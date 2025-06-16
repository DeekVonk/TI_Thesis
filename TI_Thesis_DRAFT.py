# -*- coding: utf-8 -*-
"""
TI_Thesis.py

Purpose:
    Code for the TI thesis

Version:
    1       First start 
    2       Some work on descriptive stats
    3       Basic descriptive stats + preliminary plots done
    4       Commented out dysfunctional code
    5       Summary stats on NY and NJ
    6       I don't know what I did in this version to be honest
    7       Starting aggregation function
    8       Aggregation function done
    9       Small changes
    10      Some statistics plus plots. Creates df.diff() which seems stationary
    11      Small changes, mainly started on AggregateACS.py
    12      Very close to getting AggregateACS to work, just need final year
            for general setup. Will likely run into more issues for the other
            variables...
    13      Function works after weird changes (see local .py file). Housing
            done in old way.
    14      Added removal of Census Tracts and interpolation function. Now
            starting on panel structure
    15      Panel structure imposed (can easily add regressors) and first quick
            FE OLS done. 
    16      Small changes
    17      Minute changes
    18      Added OOH and average kwh per dollar
    19      Race (White/Black) added
    20      Basically no changes
    21      First AB/AH estimator there, need to work on instrument selection
    22      Minute changes
    24      Small changes to export variables
    25      Fixing interpolate_trim() to account for lag in ACS survey variables
            Now also contains housing density, though it doesn't do much
    26      Created N2_L1 variable, extremely significant in FE_OLS and makes
            the sign of N1_L1 negative. Interesting stuff. Also added grad.py
    27      Starting MLE of Poisson. Removed commented code to comments.py
    28      Continuing Poisson MLE. Main, Objective and Poisson_RE written
            but not yet checked if it works.
    29      Poisson MLE testing. Function runs without error but no convergence.
            Might have to do with the high amount of zeroes. Check ZIP and
            built that tomorrow.
    30      Built ZIP(tau), intercept and tau parameter seem to influence each other
            but the others are relatively constant. Needs further investigation.
            May also want to build ZIP(.).
    31      Small changes (I think they are all reverted already)
    32      Creating N2 variable before exporting for GAUSS
    33      Changed variable names for GAUSS-related stuff
    34      Starting on PoissonFE
    35      Extracting cumulative number of installed panels per tract
    36      Small changes
    37      Small changes
    38      Median income and technology adjusted for inflation
    39      Currently get a singular matrix with the intercept removed. Might
            plug it back because then there are results, but I do think the
            second order neighbours should be removed. Could also do N2_L2 instead,
            works nicely as a second lag.
    40      Started on simulation function
    41      Simulation works, but means are much much lower than in the empirical
            setting. First results show that AR(1) component seems to be estimated
            quite well. Next up is adding AR(2) as control and see what happens
    42      Going to run it for first tests with N0_L0, N1_L1, N2_L1
    43      Creating more variables in panel and sim
    44      Empirical FE poisson done but standard errors don't work yet
    45      Simulations ready to run with 6 parameter configurations
    46      Simulations done. 12356 in one environment, 4 in the other. Written to 
            .csv, now hope that changing environment to extract 4 works.
    47      Starting simulation table
    48      Table code done
    49      Trying to fix CMLE standard errors
    50      Standard errors fixed and CMLE FE empirical estimation done
    51      All code done, second-to-last draft version as docstring is not done
    DRAFT   Draft version
    
Date:
    2025/06/15

Author:
    Dennis Cordes       681369
"""
###########################################################
### Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import os
import datetime
import statistics
import math
import sympy
import scipy as sp
import scipy.special as ss
import scipy.stats as st
import scipy.optimize as opt
import scipy
import statsmodels.tsa.stattools as smt
# import statsmodels.api as sm
# from statsmodels.tsa.arima.model import ARIMA
from timeit import default_timer as timer
from collections import Counter
import re
# import seaborn as sns
from linearmodels.panel import PanelOLS
from  pydynpd import regression

## Set working directory
os.chdir('C:/Users/corde/Tinbergen/Thesis')

from grad import *                  # Hessian and related functions
from AdjacencyMatrix_DRAFT import *    # Creates adjacency matrix
from AggregateACS_DRAFT import *      # Aggregates observations from ACS

### Make class
class cMaximizeLikelihood:
    def __init__(self):
            self.x0 = []
            self.x = []
            self.tx0 = []
            self.tx = []
            self.likelihoodvalue = []
            self.tcovariancematrix = []
            self.covariancematrix = []
            self.standarderrors = []
            self.tstandarderrors = []
            # self.filter = []
            self.success = False
            self.criteria = []

class cSimulations:
    def __init__(self):
            self.vTheta0 = []
            self.mTheta_mod1 = []
            self.mTheta_mod2 = []
            self.mTheta_mod3 = []
            self.vOtherParams = []
            self.NSims = []
            self.success = False


###########################################################
### readData(sFileName)
def readData(sFileName, sSource):
    """
    Purpose:  
        Loading in data for .csv file
        
    Inputs:
        sFileName       String, contains name of .csv file
        sSource         String, source of file:
                            TtS         Tracking the Sun (not worked out)
                            NYS         NYSERDA
        
    Return values
        df              Dataframe containing many, many things
    
    """
    df = pd.read_csv(sFileName)
    # df['datadate'] = pd.to_datetime(df['datadate'], format='%Y%m%d')
    
    if sSource=='TtS':
        df = df[df['customer_segment']=='RES']          # Residential
        df = df.reset_index()
        
        # # vDates = df['installation_date']
        # vYear=[]
        # for i in range(df.shape[0]):
        #     if df['installation_date'][i] != 'NaT':
        #         vYear = np.append(vYear,int(df['installation_date'][i][-4:]))
        #     else:
        #         vYear = np.append(vYear,9999)
        # vYear = int(vYear)
        
        df['year'] = df['installation_date'].str.extract('(\d{4})').fillna('9999').astype('int64')
        df = df.drop(df[(df.year < 1950)].index)
        df = df.drop(df[(df.year > 2023)].index)
        
        df.set_index(pd.to_datetime(df['installation_date']), inplace=True)
    elif sSource=='NYS':
        print(df['Sector'].unique())
        print(df.shape)
        # print(df['Community Distributed Generation'].value_counts())
        # print(df[['Community Distributed Generation','Expected KWh Annual Production']].dropna().groupby(['Community Distributed Generation'])).mean()['Expected KWh Annual Production']
        
        # print(df.groupby(['Community Distributed Generation'])).mean()['Expected KWh Annual Production']
        print(df.groupby('Community Distributed Generation')['Expected KWh Annual Production'].mean())
        dfCS = df[df['Community Distributed Generation']=='Yes']
        df = df[df['Sector']=='Residential']          # Residential
        # print(df['Community Distributed Generation'].value_counts())
        df = df.reset_index()
        df['Date'] = pd.to_datetime(df['Date Completed'])
        df.set_index(pd.to_datetime(df['Date Completed']), inplace=True)
        df['KWHpD'] = df['Expected KWh Annual Production'] / df['Project Cost']
        df['Total NYSERDA Incentive'].fillna(0)

    else:
        print('Incorrect data source specified.')

    ##
    # dfPart = df.iloc[:100,:]

    
    # vYear = np.vectorize(lambda s: s[-4:])(vDates)
    # vYearUnique, vYearCounts = np.unique(vYear, return_counts=True)
    # mYearCount = np.hstack((vYearUnique.reshape(-1,1),vYearCounts.reshape(-1,1)))

    # vStateUnique, vStateCounts = np.unique(df['state'].values, return_counts=True)
    # mStateCount = np.hstack((vStateUnique.reshape(-1,1),vStateCounts.reshape(-1,1)))
    
    
    # dfFull = df.dropna()

    return df, dfCS

###########################################################
### readData(sFileName)
def DescriptiveStats(df):
    """
    Purpose:  
        Plots the installations per year and cumulative installations. Furthermore,
        creates some descriptive statistics of key variables
        
    Inputs:
        df              Dataframe containing many, many things
        
    Return values
        dfStats         Dataframe, contains descriptive stats
    
    """
    vStates = np.sort(df['state'].unique())
    
    dfStats = pd.DataFrame(index=vStates, columns=['Observations', 'Price Mean', 'Price Std', 'Size Mean', 'Size Std', 
                                                   'Rebate', 'Rebate (% missing)', 'Grount-mounted', 'Grount-mounted (% missing)'])
        
    for s,state in enumerate(vStates):
        print('State: '+state)
        dfTemp = df[df.state==state]
        vYears = np.sort(dfTemp['year'].unique())
        vObsYear = np.zeros_like(vYears)
        for y, year in enumerate(vYears):
            # vObsYear[y] = ((df['state'] == state) & (df['year'] == year)).sum()
            vObsYear[y] = (dfTemp['year'] == year).sum()
        vObsCum = np.cumsum(vObsYear)
        
        dfStats['Observations'][s] = dfTemp.shape[0]
        dfStats['Price Mean'][s] = dfTemp[dfTemp['total_installed_price']!=-1]['total_installed_price'].mean()
        dfStats['Price Std'][s] = dfTemp[dfTemp['total_installed_price']!=-1]['total_installed_price'].std()
        dfStats['Size Mean'][s] = dfTemp[dfTemp['PV_system_size_DC']!=-1]['PV_system_size_DC'].mean()
        dfStats['Size Std'][s] = dfTemp[dfTemp['PV_system_size_DC']!=-1]['PV_system_size_DC'].std()
        # vPrice_SizeMean[s] = np.zeros_like(vStates)
        # vPrice_SizeStd[s] = np.zeros_like(vStates)
        dfStats['Rebate'][s] = dfTemp[dfTemp['rebate_or_grant']!=-1]['rebate_or_grant'].mean()
        dfStats['Rebate (% missing)'][s] = (dfTemp['rebate_or_grant']==-1).sum()/dfTemp.shape[0]
        dfStats['Grount-mounted'][s] = dfTemp[dfTemp['ground_mounted']!=-1]['ground_mounted'].mean()
        dfStats['Grount-mounted (% missing)'][s] = (dfTemp['ground_mounted']==-1).sum()/dfTemp.shape[0]
        
        fig, ax = plt.subplots(2, 1, figsize=(12,12))
        fig.suptitle('Installations for: '+state, fontsize=50)
        
        ax[0].plot(vYears,vObsYear, linewidth=5)
        ax[0].grid()
        ax[0].title.set_text('Installed systems per year')
        
        ax[1].plot(vYears,vObsCum, linewidth=5)
        ax[1].grid()
        ax[1].title.set_text('Cumulative installed systems')
        
        fig.show()

    return dfStats

###########################################################
### readData(sFileName)
def AggregateObservations(df, dfOrder, sInterval='Q'):
    """
    Purpose:  
        Counts the number of observations per census tract in a pre-specified
        time interval. Start and end date are fixed so time interval should take
        this into account. Easy solution is to make the time interval fit in a
        single calendar year.
        
    Inputs:
        df              Dataframe containing many, many things
        dfOrder         Dataframe with all census tracts, taken from adjacency matrix
        
    Return values
        dfCounts        (Quarters x Tracts) dataframe with observation counts
    
    """
    # I think the best course of action is to create a big matrix (iT x iN)
    # and just fill that with a loop. Can use aggregation function such as 
    # resample and then over a loop apply that for each census tract. Not the 
    # most beautiful solution but it should work. 
    
    tsFirst = df.index.min()
    tsLate = df.index.max()
    rng = rng = pd.date_range(tsFirst, tsLate, freq=sInterval)
    
    dfUnique = df['Census Tract'].unique()
    
    dfCounts = pd.DataFrame(index=rng, columns=dfOrder)
    dfKWHpD = pd.DataFrame(index=rng, columns=dfOrder)
    dfGrant = pd.DataFrame(index=rng, columns=dfOrder)
    dfCost = pd.DataFrame(index=rng, columns=dfOrder)
    
    for c in dfOrder:
        # print(c)
        if c in dfUnique:
            dfCounts[c] = df[df['Census Tract']==c]['Census Tract'].resample(sInterval).count()
            dfKWHpD[c] = df[(df['Census Tract']==c) & (df['Project Cost']!=0)]['KWHpD'].resample(sInterval).mean()
            # dfKWHpD[c] = df[df['Census Tract']==c]['Total NYSERDA Incentive'].resample('Q').mean()
            dfGrant[c] = df[(df['Census Tract']==c)]['Total NYSERDA Incentive'].resample(sInterval).mean()
            dfCost[c] = df[(df['Census Tract']==c)]['Project Cost'].resample(sInterval).median()
    
    dfCounts = dfCounts.fillna(0)
    
    
    return dfCounts, dfKWHpD, dfGrant, dfCost

###########################################################
### readData(sFileName)
def RemoveCTs(dfCounts, dfMedInc, dfAge, dfOOH, dfWhite, dfBlack, dfOrder, mAdj, dfALAND, dfRemove):
    """
    Purpose:  
        Removes rows from dataframe matrices, adjacency matrix, and accompanying
        order series of census tracts, based on series object dfRemove which has
        
    Inputs:
        df              Dataframe containing many, many things
        dfOrder         Dataframe with all census tracts, taken from adjacency matrix
        
    Return values
        dfCounts        (Quarters x Tracts) dataframe with observation counts
    
    """
    print(str(dfRemove.sum())+' values removed')
    dfKeep = ~dfRemove
    dfMedInc2 = dfMedInc[dfKeep.index[dfKeep]]
    dfAge2 = dfAge[dfKeep.index[dfKeep]]
    dfOOH2 = dfOOH[dfKeep.index[dfKeep]]
    dfWhite2 =  dfWhite[dfKeep.index[dfKeep]]
    dfBlack2 =  dfBlack[dfKeep.index[dfKeep]]
    dfCounts2 = dfCounts[dfKeep.index[dfKeep]]
    dfALAND2 = dfALAND[dfKeep.index[dfKeep]]
    # dfALAND2 = dfALAND    
    mAdj2 = mAdj[dfKeep,:]
    mAdj2 = mAdj2[:,dfKeep]
    
    dfOrder2 = dfOOH2.columns
    
    return dfCounts2, dfMedInc2, dfAge2, dfOOH2, dfWhite2, dfBlack2, dfOrder2, mAdj2, dfALAND2

###########################################################
### readData(sFileName)
def Interpolate_trim(dfCounts, dfKWHpD, dfMedInc, dfAge, dfOOH, dfWhite, dfBlack, sMethod = 'explorative'):
    """
    Purpose:  
        Linearly interpolates values for the given matrices. Also trims them.
        
    Inputs:
        dfCounts        Dataframe with counts
        dfs             Dataframes with variables
        sMethod         Takes the following values:
                            - 'explorative'     uses future information
                            - 'predictive'      does not use future info
        
    Return values
        dfCounts        (Quarters x Tracts) dataframe with observation counts
    
    """
    if sMethod == 'explorative':
        dfCounts2 = dfCounts.loc[(dfCounts.index >= '2010-9-30') & (dfCounts.index <= '2023-9-30')]
        
        dfKWHpD_trim = dfKWHpD.loc[(dfCounts.index >= '2010-9-30') & (dfCounts.index <= '2023-9-30')]
        dfKWHpD2 = pd.DataFrame().reindex_like(dfCounts2)
        dfKWHpD2.iloc[:,:] = np.tile(dfKWHpD_trim.values.reshape(-1,1),(1,dfKWHpD2.shape[1]))
        
        dfMedInc2 = dfMedInc.loc[(dfMedInc.index >= '2010-9-30') & (dfMedInc.index <= '2023-9-30')]
        dfMedInc2 = dfMedInc2.interpolate(method='linear', axis=0, limit=1, limit_direction='forward')
        
        dfAge2 = dfAge.loc[(dfAge.index >= '2010-9-30') & (dfAge.index <= '2023-9-30')]
        dfAge2 = dfAge2.interpolate(method='linear', axis=0, limit=1, limit_direction='forward')

        dfOOH2 = dfOOH.loc[(dfOOH.index >= '2010-9-30') & (dfOOH.index <= '2023-9-30')]
        dfOOH2 = dfOOH2.interpolate(method='linear', axis=0, limit=1, limit_direction='forward')
        
        dfWhite2 = dfWhite.loc[(dfWhite.index >= '2010-9-30') & (dfWhite.index <= '2023-9-30')]
        dfWhite2 = dfWhite2.interpolate(method='linear', axis=0, limit=1, limit_direction='forward')
        
        dfBlack2 = dfBlack.loc[(dfBlack.index >= '2010-9-30') & (dfBlack.index <= '2023-9-30')]
        dfBlack2 = dfBlack2.interpolate(method='linear', axis=0, limit=1, limit_direction='forward')
        
    if sMethod == 'explanatory':
        dfCounts_m1 = dfCounts.loc[(dfCounts.index == '2010-9-30')]
        dfCounts_0 = dfCounts.loc[(dfCounts.index == '2011-3-31')]
        dfCounts2 = dfCounts.loc[(dfCounts.index >= '2011-9-30') & (dfCounts.index <= '2024-9-30')]
        
        dfKWHpD_trim = dfKWHpD.loc[(dfCounts.index >= '2011-9-30') & (dfCounts.index <= '2024-9-30')]
        dfKWHpD2 = pd.DataFrame().reindex_like(dfCounts2)
        dfKWHpD2.iloc[:,:] = np.tile(dfKWHpD_trim.values.reshape(-1,1),(1,dfKWHpD2.shape[1]))
        
        dfMedInc2 = dfMedInc.loc[(dfMedInc.index >= '2010-9-30') & (dfMedInc.index <= '2023-9-30')]
        dfMedInc2 = dfMedInc2.interpolate(method='linear', axis=0, limit=3, limit_direction='forward')
        dfMedInc2.index = dfCounts2.index
        
        dfAge2 = dfAge.loc[(dfAge.index >= '2010-9-30') & (dfAge.index <= '2023-9-30')]
        dfAge2 = dfAge2.interpolate(method='linear', axis=0, limit=3, limit_direction='forward')
        dfAge2.index = dfCounts2.index
    
        dfOOH2 = dfOOH.loc[(dfOOH.index >= '2010-9-30') & (dfOOH.index <= '2023-9-30')]
        dfOOH2 = dfOOH2.interpolate(method='linear', axis=0, limit=3, limit_direction='forward')
        dfOOH2.index = dfCounts2.index
        
        dfWhite2 = dfWhite.loc[(dfWhite.index >= '2010-9-30') & (dfWhite.index <= '2023-9-30')]
        dfWhite2 = dfWhite2.interpolate(method='linear', axis=0, limit=3, limit_direction='forward')
        dfWhite2.index = dfCounts2.index
        
        dfBlack2 = dfBlack.loc[(dfBlack.index >= '2010-9-30') & (dfBlack.index <= '2023-9-30')]
        dfBlack2 = dfBlack2.interpolate(method='linear', axis=0, limit=3, limit_direction='forward')
        dfBlack2.index = dfCounts2.index
        
        for date in dfCounts2.index:
            dfMedInc_nan = dfMedInc2.loc[date].isnull()
            dfMedInc2.loc[date,dfMedInc_nan[dfMedInc_nan == True].index] = dfMedInc2.loc[date].median()
            
            dfAge_nan = dfAge2.loc[date].isnull()
            dfAge2.loc[date,dfAge_nan[dfAge_nan == True].index] = dfMedInc2.loc[date].median()

            dfOOH_nan = dfOOH2.loc[date].isnull()
            dfOOH2.loc[date,dfOOH_nan[dfOOH_nan == True].index] = dfOOH2.loc[date].median()

            dfWhite_nan = dfWhite2.loc[date].isnull()
            dfWhite2.loc[date,dfWhite_nan[dfWhite_nan == True].index] = dfWhite2.loc[date].median()

            dfBlack_nan = dfBlack2.loc[date].isnull()
            dfBlack2.loc[date,dfBlack_nan[dfBlack_nan == True].index] = dfBlack2.loc[date].median()
            

    return dfCounts2, dfCounts_0, dfCounts_m1, dfKWHpD2, dfMedInc2, dfAge2, dfOOH2, dfWhite2, dfBlack2

###########################################################
### readData(sFileName)
def Create_panel(lDFs):
    """
    Purpose:  
        Takes outcome variable and list of independent variables and creates a
        panel structure
        
    Inputs:
        dfCounts        Dataframe with counts
        lDFs            List with dataframes containing independent variables.
                        Dataframes need to have a .name attribute which then
                        becomes a column name
        
    Return values:
        dfIndex         Dataframe, has CT and Date as index
        df              Dataframe, has index reset
    """
    iT, iN = lDFs[0].shape
    iK = len(lDFs)
    vCols = ['Census Tract','Date', 'Time', 'Trump dummy']
    for k in lDFs:
        vCols.append(k.name)
    
    df = pd.DataFrame(index=range(iT*iN), columns=vCols)
    
    for t,c in enumerate(lDFs[0].columns):
        # print(t,c)
        df['Census Tract'].iloc[t*iT:(t+1)*iT] = c
    df['Date'] = np.tile(lDFs[0].index,iN)
    df['Timeperiod'] = np.tile(np.arange(1,iT+1),iN)
    df['Trump dummy'].loc[(df['Date'] >= '2017-01-01') & (df['Date'] <= '2021-01-01')] = 1
    df['Trump dummy'] = df['Trump dummy'].fillna(0)
    df['Intercept'] = 1
    # df.loc[(df['Date'] < '2017-01-01')] = 0
    
    for k in lDFs:
        print('Starting '+k.name)
        for t,c in enumerate(lDFs[0].columns):
            df[k.name].iloc[t*iT:(t+1)*iT] = k.iloc[:,t]
         
    df['Time x Trump'] = df['Timeperiod'] * df['Trump dummy']
    df['MedianInc'] = df['MedianInc']/1000
    
    dfIndex = df.set_index(['Census Tract', 'Date'], inplace=False)
    df.columns = df.columns.str.replace(' ', '_')
    
    return dfIndex, df


###########################################################
### readData(sFileName)
def CreateHousingDensity(dfOOH, dfALAND):
    """
    Purpose:  
        Divides housing density by land area
        
    Inputs:
        dfOOH           Dataframe, owner-occupied housing units
        dfALAND         Dataframe, land area
        
    Return values:
        dfHD            Dataframe, housing density
    """
    dfHD = pd.DataFrame().reindex_like(dfOOH)
    
    for t in range(dfOOH.shape[0]):
        # print(t)
        dfHD.iloc[t,:] = dfOOH.iloc[t,:] / dfALAND
    
    dfHD.name = 'HD'
    
    return dfHD

###########################################################
### readData(sFileName)
def Adjust_Prices(dfMedInc, dfKWHpD):
    """
    Purpose:  
        Adjusts prices back to 2010 dollars
        
    Inputs:
        dfMedInc        Dataframe, median income
        dfKWHpD         Dataframe, expected yearly output per dollar
        
    Return values:
        dfMedInc2        Dataframe, updated median income
        dfKWHpD2         Dataframe, updated expected yearly output per dollar
    """
    dfInflation = pd.read_csv('Data/inflation_data.csv')
    dfInflation.index = pd.to_datetime(dfInflation['year'], format="%d/%m/%Y")
    dfInfFull = pd.DataFrame(index=dfMedInc.index, columns=["amount"])
    # dfInfFull.merge(dfInflation, how="inner", on="")
    
    for date in dfInflation.index:
        dfInfFull.loc[date] = dfInflation.loc[date]['amount']
    dfInfFull = dfInfFull.iloc[:-2,:].astype('float64')
    dfInfFull = dfInfFull.interpolate()    
    
    dfMedInc2 = dfMedInc / dfInfFull.values
    dfKWHpD2 = (dfKWHpD * dfInfFull.values)*100
    
    return dfMedInc2, dfKWHpD2

###########################################################
### readData(sFileName)
def CumPanelsVector(df):
    """
    Purpose:  
        Creates vector with cumulative installed panels, with GEOID as index.
        Saves it as a .csv in the GIS map
        
    Inputs:
        
    Return values
        
    
    """
    dfGetCTs = pd.read_csv('GIS/NY_Neighbours.csv')['GEOID']
    dfCum = df.cumsum(axis=0)
    dfTotal = pd.DataFrame(index=np.sort(dfGetCTs))
    dfTotal['cum'] = 0
    dfTotal['GEOID'] = dfTotal.index
    
    for ct in dfCum.columns:
        dfTotal['cum'].loc[ct] = dfCum[ct].iloc[-1]
    
    dfTotal.to_csv('GIS/cumpanels.csv')
    
    return 


###########################################################
### readData(sFileName)
def first_last_missing(lDFs):
    """
    Purpose:  
        Fills in 2011 and 2024 value if 2012 or 2023 value is present
        
    Inputs:
        
    Return values
        
    
    """
    for dfr in lDFs:
        for col in dfr.columns:
            if ((pd.isna(dfr[col].loc['2011-09-30']) == True) & (pd.isna(dfr[col].loc['2012-09-30']) == False)):
                dfr[col].loc['2011-09-30'] = dfr[col].loc['2012-09-30']
            if ((pd.isna(dfr[col].loc['2024-09-30']) == True) & (pd.isna(dfr[col].loc['2023-09-30']) == False)):
                dfr[col].loc['2024-09-30'] = dfr[col].loc['2023-09-30']
    
    return lDFs


###########################################################
### readData(sFileName)
def Create_AR_Neighbours(df, mAdj, df_0, df_m1):
    """
    Purpose:  
        Takes the count dataframe and adjacency matrix and creates 0,1, and 2-
        lagged versions of nbr0, nbr1, and nbr2.
        
    Inputs:
        dfCounts        Dataframe with counts
        
    Return values:
        Loads, output speaks for itself. 
    
    """
    dfCount_L1 = pd.DataFrame().reindex_like(df)
    dfN1 = pd.DataFrame().reindex_like(df)
    dfN2 = pd.DataFrame().reindex_like(df)
    dfN1_L1 = pd.DataFrame().reindex_like(df)
    dfN2_L1 = pd.DataFrame().reindex_like(df)
    
    # new
    dfN0_L2 = pd.DataFrame().reindex_like(df)
    dfN1_L2 = pd.DataFrame().reindex_like(df)
    dfN2_L2 = pd.DataFrame().reindex_like(df)
    
    
    dfCount_L1.iloc[0,:] = df_0.values
    dfCount_L1.iloc[1:,:] = df.iloc[:-1,:]
    
    dfN0_L2.iloc[0,:] = df_m1.values
    dfN0_L2.iloc[1,:] = df_0.values
    dfN0_L2.iloc[2:,:] = df.iloc[:-2,:]
    
    dfN1_L1.iloc[0,:] = (mAdj @ df_0.values.reshape(-1,1)).T
    dfN2_L1.iloc[0,:] = (mAdj @ mAdj @ df_0.values.reshape(-1,1)).T
    
    dfN1_L2.iloc[0,:] = (mAdj @ df_m1.values.reshape(-1,1)).T
    dfN1_L2.iloc[1,:] = (mAdj @ df_0.values.reshape(-1,1)).T
    
    dfN2_L2.iloc[0,:] = (mAdj @ mAdj @ df_m1.values.reshape(-1,1)).T
    dfN2_L2.iloc[1,:] = (mAdj @ mAdj @ df_0.values.reshape(-1,1)).T

    
    for t in range(0,df.shape[0]-1):
        vX = df.iloc[t,:].values.reshape(-1,1)
        dfN1.iloc[t,:] = (mAdj @ vX).T
        dfN2.iloc[t,:] = (mAdj @ mAdj @ vX).T
        dfN1_L1.iloc[t+1,:] = (mAdj @ vX).T
        dfN2_L1.iloc[t+1,:] = (mAdj @ mAdj @ vX).T
        
        
    dfN1.iloc[-1,:] = (mAdj @ df.iloc[-1,:].values.reshape(-1,1)).T
    dfN2.iloc[-1,:] = (mAdj @ mAdj @ df.iloc[-1,:].values.reshape(-1,1)).T
    
    # dfN0_L2.iloc[2:,:] = df.iloc[:-2,:]
    dfN1_L2.iloc[2:,:] = dfN1.iloc[:-2,:]
    dfN2_L2.iloc[2:,:] = dfN2.iloc[:-2,:]
    
    dfCount_L1.name = 'N0_L1'
    dfN1.name = 'N1_L0'
    dfN2.name = 'N2_L0'
    dfN1_L1.name = 'N1_L1'
    dfN2_L1.name = 'N2_L1'
    dfN0_L2.name = 'N0_L2'
    dfN1_L2.name = 'N1_L2'
    dfN2_L2.name = 'N2_L2'
    
    return dfCount_L1, dfN1, dfN2, dfN1_L1, dfN2_L1, dfN0_L2, dfN1_L2, dfN2_L2

###########################################################
### mStats = DescriptiveStats(mLogRet)
def Poisson_RE(iT, vY, mX, vBeta, dGamma, dDelta):
    """
    Purpose:
        Calculates likelihood for a random effects conditional Poisson ML
    Inputs:
        iT                  integer, T
        vY                  iNT x 1 vector of observations, sorted per individual
        mX                  iNT x k vector of regressors
        vBeta               iK x 1 vector of parameters
        dGamma              double, first Gamma-dist parameter
        dDelta              double, second Gamma-dist parameter
        
    Returns:
        dLikelihoodValues       double, sum of loglikelihood values
    """
    
    # initialize sizes and formats
    iNT = vY.shape[0]
    iN = int(vY.shape[0] / iT)
    # vVariances = np.zeros(iT)
    # vLikelihoodValues = np.ones(iNT)
    
    dLikelihoodValue = -iN*ss.loggamma(dGamma) + iN*dGamma*np.log(dDelta)
        
    for n in range(iN):
        vY_i = vY[iT*n:iT*(n+1)].reshape(-1,1)
        mX_i = mX[iT*n:iT*(n+1),:]
        iY_sum = np.sum(vY_i)
        
        vTerm1 = np.sum(vY_i * (mX_i @ vBeta))      # n_it (Xit Beta)
        vTerm2 = -(dGamma + iY_sum)*np.log(np.sum(np.exp(mX_i @ vBeta))+dDelta)
        vTerm3 = ss.loggamma(iY_sum + dGamma)
        
        
        dLikelihoodValue += (vTerm1 + vTerm2 + vTerm3)
                
        # if math.isnan(vLikelihoodValues[i1]) == True:
        #     vLikelihoodValues = -1000000
        #     print('nan found, setting likelihood to:', vLikelihoodValues)
        #     return 0, vLikelihoodValues
    
    dLikelihoodValue = dLikelihoodValue/ iNT
    # return filtered levels and likelihood values
    return dLikelihoodValue

###########################################################
### mStats = DescriptiveStats(mLogRet)
def ZIP(iT, vY, mX, vBeta, dTau):
    """
    Purpose:
        Calculates likelihood for a zero-inflated Poisson distribution
        
    Inputs:
        iT                  integer, T
        vY                  iNT x 1 vector of observations, sorted per individual
        mX                  iNT x k vector of regressors
        vBeta               iK x 1 vector of parameters
        dTau                double, ZIP parameter
        
    Returns:
        dLikelihoodValues       double, sum of loglikelihood values
    """
    
    # initialize sizes and formats
    vZero = (vY == 0).reshape(-1,1)
    vNonZero = ~vZero
    dLikelihoodValue = 0
    
    vBiBeta = (mX @ vBeta).reshape(-1,1)
    # print(vBiBeta)
    vTerm1 = np.sum(vZero * np.log(np.exp(-dTau * vBiBeta)+np.exp(-np.exp(vBiBeta))))
    vTerm2 = np.sum(vNonZero * vY * vBiBeta - np.exp(vBiBeta))
    vTerm3 = -np.sum(np.log(1+np.exp(-dTau*vBiBeta)))
    
    dLikelihoodValue += (vTerm1 + vTerm2 + vTerm3)
        
    dLikelihoodValue = dLikelihoodValue/ vY.shape[0]
    
    # return filtered levels and likelihood values
    return dLikelihoodValue

###########################################################
### mStats = DescriptiveStats(mLogRet)
def Poisson_FE(iT, vY, mX, vBeta):
    """
    Purpose:
        Calculates likelihood for a fixed effects conditional Poisson ML
    Inputs:
        iT                  integer, T
        vY                  iNT x 1 vector of observations, sorted per individual
        mX                  iNT x k vector of regressors
        vBeta               iK x 1 vector of parameters
        
    Returns:
        vLikelihoodValues   iNT x 1 vector, loglikelihood values
    """
    
    # initialize sizes and formats
    iNT = vY.shape[0]
    iN = int(vY.shape[0] / iT)
    vNonZero = np.nonzero(vY)[0].reshape(-1,1)
    # vVariances = np.zeros(iT)
    vLikelihoodValues = np.zeros((iNT,1))
    
    vSumExpXBeta = np.zeros((iNT,1))
    vExpXBeta = np.exp(mX @ vBeta)
    # vExpMinXBeta = np.exp(-mX @ vBeta)
    
    # dLikelihoodValue = 0
        
    for n in range(iN):
        vSumExpXBeta[iT*n:iT*(n+1)] = (np.ones((iT,iT)) @ vExpXBeta[iT*n:iT*(n+1)]).reshape(-1,1)
        
    for i in vNonZero:
        vLikelihoodValues[i] -= (vY[i] * np.log((vExpXBeta[i]**(-1))*vSumExpXBeta[i]))
        if math.isnan(vLikelihoodValues[i]) == True:
            vLikelihoodValues = -1000000
            print('nan found, setting likelihood to:', vLikelihoodValues)
            return 0, vLikelihoodValues
    # vLikelihoodValues += 
    
    # vLikelihoodValue = dLikelihoodValue/ iNT
    # return filtered levels and likelihood values
    return vLikelihoodValues

###########################################################
### cMaxLikelihood = MaxLikelihood(vData)
def estim_Poisson_MLE(df, lY, lPre, lExog, iT, dTau_0, sMethod='RE'):
    """
    Purpose:
        Estimates Random effects and ZIP poisson ML
    Inputs:
        df          Dataframe containing all variables
        lY          list, y-variable
        lPre        list, pre-determined variables
        lExog       list, exogenous variables
        iT          integer, T
        dTau_0      double, starting value for tau (ZIP)
        sMethod     Method, either 'RE' or 'ZIP'
        
    Returns:
        cReturnValue    Class object containing all information about the estimation
    """
    ###########################################################
    ### cMaxLikelihood = MaxLikelihood(vData)
    def Objective(vTheta, bForAllT=False):
        """
        Purpose:
            beta | gamma | delta
        Inputs:
            vTheta      Vector, containing the parameters of interest.
            bForAllT    Boolean, True if the vector of likelihoods must be return. False for average negative log likelihood.
        Returns:
            dObjValue   Vector or Double, containin the vector of loglikelihood values or the average negative log likelihood.
        """
        # initialize the parameter values
        print('.', end='')
        vTranspar = ParameterTransform(vTheta)
        # dOmega = vTranspar[0]
        # dAlpha = vTranspar[1]
        # dBeta = vTranspar[2]
        # dNu = vTranspar[3]
        # dGamma = vTranspar[4]
        # dSig = vTranspar[5]
        
        
        if sMethod=='RE':
            vBeta = vTranspar[:-2]
            # dGamma = vTranspar[-2]
            dGamma = vTranspar[-1]
            # dDelta = vTranspar[-1]
            # run the filter
            dLikelihoodValue = Poisson_RE(iT, vY, mX, vBeta, dGamma, 1)
            # dLikelihoodValue = Poisson_RE(iT, vY, mX, vBeta, dGamma, dDelta)
        elif sMethod=='ZIP':
            vBeta = vTranspar[:-1]
            dTau = vTranspar[-1]
            dLikelihoodValue = ZIP(iT, vY, mX, vBeta, dTau)
        
        if (bForAllT == True):
            dObjValue = dLikelihoodValue.reshape(-1,)
        else:
            dObjValue = -(dLikelihoodValue/iNT)
        
        return dObjValue
    
    ###########################################################
    ### cMaxLikelihood = MaxLikelihood(vData)
    def ComputeCovarianceMatrix(vTheta):
        """
        Purpose:
            Local covariance matrix function to compute the covariance matrix 
        Inputs:
            vTheta      Vector, containing the parameters of interest.
        Returns:
            mCov        Matrix, iK x iK containing the covariance matrix of the parameters
        """
        # compute the inverse hessian of the average log likelihood
        mH= hessian_2sided(Objective, vTheta)
        # print(mH)
        mCov = np.linalg.inv(mH)
        mCov = (mCov +  mCov.T)/2       #  Force to be symmetric
        # compute the outer product of gradients of the average log likelihood
        mG = jacobian_2sided(Objective, vTheta, True)
        # mG = jacobian_2sided(Objective, vTheta)
        # print(mG)
        # mG = np.dot(mG.T, mG) / vData.shape[0]
        # mG = np.dot(mG, mCov)
        # mCov = np.dot(mCov, mG) / vData.shape[0]
        
        mG = np.dot(mG.T, mG) / iNT
        mG = np.dot(mG, mCov)
        mCov = np.dot(mCov, mG) / iNT

        return mCov
    
    ###############################################
    ### main()
    def ParameterTransform(vTheta):
        """
        Purpose:
            To transform the parameters for the robust TV-mean model such that they are restricted:
                -1 < Phi_1 < 1

            Outputs the true vector of parameters
        Inputs:
            vTheta                  Vector, containing the transformed parameters of interest
        Returns:
            vTranspar               Vector, containing the true parameters.
        """
        vTranspar = np.copy(vTheta)
        vTranspar[1] = 1/(1+np.exp(-vTheta[1]))
        if sMethod == 'RE':
            # vTranspar[-2] = np.exp(vTheta[-2])
            vTranspar[-1] = np.exp(vTheta[-1])
        # elif sMethod == 'ZIP':
        #     vTranspar[-1] = 1/(1+np.exp(-vTheta[-1]))

        return vTranspar
    
    ###############################################
    ### main()
    def CalculateCriteria(dLL):
        """
        Purpose:
            Calculates information criteria (AIC, AICc, BIC)
            
        Inputs:
            dLL         Likelihood value
    
        Returns:
            dAIC        double, AIC
            dAICc       double, AICc
            dBIC        double, BIC
        """
        # dLL = -dLL
        dAIC= 2*iK - 2*dLL
        dAICc= dAIC + ((2*iK**2 + 2*iK)/(iT-iK-1))
        dBIC= iK*np.log(iT) - 2*dLL

        return (dAIC, dAICc, dBIC)
    
    # initialize starting values and return value
    cReturnValue = cMaximizeLikelihood()
    if sMethod == 'RE':
        # iK = 1+ len(lPre) + len(lExog) + 1  # intercept + predetermined + exog + 2 (gamma)
        iK = 1+ len(lPre) + len(lExog) + 2  # intercept + predetermined + exog + 2 (gamma,delta)
    elif sMethod == 'FE':
        iK = 1+ len(lPre) + len(lExog) + 2  # intercept + predetermined + exog
    elif sMethod == 'ZIP':
        iK = 1+ len(lPre) + len(lExog) + 1  # intercept + predetermined + exog + 1 (tau)
        # iK = len(lPre) + len(lExog) + 1     # predetermined + exog + 1 (tau)
    
    vTheta = np.zeros(iK).reshape(-1,1)
    
    iNT = df.shape[0]
    # dMu0 = np.mean(vData)
    
    # Order: intercept | predetermined | exogenous | [gamma, delta]    
    # exact order is given in list arguments of function
    # first predetermined coefficient is always AR(1), of which the parameter
    # is bounded between -1 and 1. When RE is used, gamma distribution parameters
    # are strictly positive
    
    vTheta[1] = np.log(0.4/(1-0.4))                   # AR(1)
    if sMethod == 'RE':
        # vTheta[-2] = np.log(0.5)
        vTheta[-1] = np.log(0.5)
    elif sMethod == 'ZIP':
        vTheta[-1] = dTau_0

    lFull = []
    # lFull.extend(lY)
    lFull.extend(lPre)
    lFull.extend(lExog) 
    dfOrdered = df[lFull]

    vY = df[lY].values.astype('int')
    mX = np.hstack((np.ones(iNT).reshape(-1,1), dfOrdered.values.astype('float')))
    # mX = dfOrdered.values.astype('float')

    # vTheta = np.log(vTheta)
    cReturnValue.x0 = vTheta
    cReturnValue.tx0 = ParameterTransform(vTheta)
    
    # do the optimization
    tSol = opt.minimize(Objective, vTheta, method='BFGS', options={'disp': True, 'maxiter':250})
    cReturnValue.success = tSol['success']

    # check for success and store results
    if (tSol['success'] != True):
        print("*** no true convergence: ",tSol['message'])
    else:
        cReturnValue.x = tSol['x']
        cReturnValue.tx = ParameterTransform(cReturnValue.x)
        print(cReturnValue.tx0)
        print(cReturnValue.tx)
        cReturnValue.likelihoodvalue = -iNT * tSol['fun']
        cReturnValue.criteria = CalculateCriteria(cReturnValue.likelihoodvalue)
        cReturnValue.covariancematrix = ComputeCovarianceMatrix(cReturnValue.x)
        # mJ = jacobian_2sided(ParameterTransform, cReturnValue.x, True)
        mJ = jacobian_2sided(ParameterTransform, cReturnValue.x)
        cReturnValue.tcovariancematrix = np.dot(mJ, np.dot(cReturnValue.covariancematrix, mJ.T))
        cReturnValue.standarderrors = np.sqrt(np.diag(cReturnValue.covariancematrix))
        cReturnValue.tstandarderrors = np.sqrt(np.diag(cReturnValue.tcovariancematrix))
            
    return cReturnValue


###########################################################
### cMaxLikelihood = MaxLikelihood(vData)
def estim_PoissonFE_MLE(df, lY, lPre, lExog, iT):
    """
    Purpose:
        Estimates Random effects and ZIP poisson ML
    Inputs:
        df          Dataframe containing all variables
        lY          list, y-variable
        lPre        list, pre-determined variables
        lExog       list, exogenous variables
        iT          integer, T
        
    Returns:
        cReturnValue    Class object containing all information about the estimation
    """
    ###########################################################
    ### cMaxLikelihood = MaxLikelihood(vData)
    def Objective(vTheta, bForAllT=False):
        
        # initialize the parameter values
        print('.', end='')
        vBeta = ParameterTransform(vTheta)
        
        vLikelihoodValues = Poisson_FE(iT, vY, mX, vBeta)

        if (bForAllT == True):
            dObjValue = vLikelihoodValues.reshape(-1,)
        else:
            dObjValue = -np.mean(vLikelihoodValues)
            # dObjValue = -np.sum(vLikelihoodValues)
        
        return dObjValue
    
    ###########################################################
    ### cMaxLikelihood = MaxLikelihood(vData)
    def ComputeCovarianceMatrix(vTheta):
        """
        Purpose:
            Local covariance matrix function to compute the covariance matrix 
        Inputs:
            vTheta      Vector, containing the parameters of interest.
        Returns:
            mCov        Matrix, iK x iK containing the covariance matrix of the parameters of the univariate garch
        """
        # compute the inverse hessian of the average log likelihood
        mH= hessian_2sided(Objective, vTheta)
        # print(mH)
        mCov = np.linalg.inv(mH)
        mCov = (mCov +  mCov.T)/2       #  Force to be symmetric
        # compute the outer product of gradients of the average log likelihood
        mG = jacobian_2sided(Objective, vTheta, True)
        # mG = jacobian_2sided(Objective, vTheta)
        # print(mG)
        # mG = np.dot(mG.T, mG) / vData.shape[0]
        # mG = np.dot(mG, mCov)
        # mCov = np.dot(mCov, mG) / vData.shape[0]
        
        mG = np.dot(mG.T, mG) / iNT
        mG = np.dot(mG, mCov)
        mCov = np.dot(mCov, mG) / iNT

        return mCov
    
    ###############################################
    ### main()
    def ParameterTransform(vTheta):
        """
        Purpose:
            To transform the parameters for the robust TV-mean model such that they are restricted:
                -1 < Phi_1 < 1

            Outputs the true vector of parameters
        Inputs:
            vTheta                  Vector, containing the transformed parameters of interest
        Returns:
            vTranspar               Vector, containing the true parameters.
        """
        vTranspar = np.copy(vTheta)
        # vTranspar[1] = 1/(1+np.exp(-vTheta[1]))
        vTranspar[0] = 1/(1+np.exp(-vTheta[0]))     # changed for no intercept
        # elif sMethod == 'ZIP':
        #     vTranspar[-1] = 1/(1+np.exp(-vTheta[-1]))

        return vTranspar
    
    ###############################################
    ### main()
    def CalculateCriteria(dLL):
        """
        Purpose:
            Calculates information criteria (AIC, AICc, BIC)
            
        Inputs:
            dLL         Likelihood value
    
        Returns:
            dAIC        double, AIC
            dAICc       double, AICc
            dBIC        double, BIC
        """
        dLL = -dLL
        dAIC= 2*iK - 2*dLL
        dAICc= dAIC + ((2*iK**2 + 2*iK)/(iT-iK-1))
        dBIC= iK*np.log(iT) - 2*dLL

        return (dAIC, dAICc, dBIC)
    
    # initialize starting values and return value
    cReturnValue = cMaximizeLikelihood()
    # iK = 1+ len(lPre) + len(lExog)  # intercept + predetermined + exog
    iK = len(lPre) + len(lExog)  # predetermined + exog
    
    # vTheta = np.zeros(iK).reshape(-1,1)
    vTheta = np.ones(iK).reshape(-1,1)
    
    iNT = df.shape[0]
    # dMu0 = np.mean(vData)
    
    # Order: intercept | predetermined | exogenous | [gamma, delta]    
    # exact order is given in list arguments of function
    # first predetermined coefficient is always AR(1), of which the parameter
    # is bounded between -1 and 1. When RE is used, gamma distribution parameters
    # are strictly positive
    
    # vTheta[1] = np.log(0.4/(1-0.4))                   # AR(1)
    vTheta[0] = np.log(0.2/(1-0.2))                   # AR(1) - changed for no intercept
    # if sMethod == 'RE':
    #     vTheta[-2] = np.log(0.5)
    #     vTheta[-1] = np.log(0.5)
    # elif sMethod == 'ZIP':
    #     vTheta[-1] = dTau_0

    lFull = []
    # lFull.extend(lY)
    lFull.extend(lPre)
    lFull.extend(lExog) 
    dfOrdered = df[lFull]

    vY = df[lY].values.astype('int')
    # mX = np.hstack((np.ones(iNT).reshape(-1,1), dfOrdered.values.astype('float')))
    mX = dfOrdered.values.astype('float') # changed for no interept
    # mX = dfOrdered.values.astype('float')

    # vTheta = np.log(vTheta)
    cReturnValue.x0 = vTheta
    cReturnValue.tx0 = ParameterTransform(vTheta)
    
    # do the optimization
    tSol = opt.minimize(Objective, vTheta, method='BFGS', options={'disp': True, 'maxiter':250})
    cReturnValue.success = tSol['success']

    # check for success and store results
    if (tSol['success'] != True):
        print("*** no true convergence: ",tSol['message'])
    else:
        cReturnValue.x = tSol['x']
        cReturnValue.tx = ParameterTransform(cReturnValue.x)
        print(cReturnValue.tx0)
        print(cReturnValue.tx)
        cReturnValue.likelihoodvalue = -iNT * tSol['fun']
        cReturnValue.criteria = CalculateCriteria(cReturnValue.likelihoodvalue)
        cReturnValue.covariancematrix = ComputeCovarianceMatrix(cReturnValue.x)
        # mJ = jacobian_2sided(ParameterTransform, cReturnValue.x, True)
        mJ = jacobian_2sided(ParameterTransform, cReturnValue.x)
        cReturnValue.tcovariancematrix = np.dot(mJ, np.dot(cReturnValue.covariancematrix, mJ.T))
        cReturnValue.standarderrors = np.sqrt(np.diag(cReturnValue.covariancematrix))
        cReturnValue.tstandarderrors = np.sqrt(np.diag(cReturnValue.tcovariancematrix))
            
    return cReturnValue

###########################################################
### mStats = DescriptiveStats(mLogRet)
def Plot_CumMargPanels(df):
    """
    Purpose:
        Plot the cumulative and marginal number of solar panels per year
        
    Inputs:
        df          Dataframe with counts per year
        
    Returns:
        
    """
    dfPlot = pd.DataFrame(index=df.index)
    dfPlot['Cumulative'] = (df.cumsum()).sum(axis=1)-dfMarginal
    dfPlot['Marginal'] = df.sum(axis=1)
    dfPlot.index = dfPlot.index.year
    
    ax = dfPlot.plot.bar(stacked=True, figsize=(10,5))
    ax.set_ylabel('Number of residential PV projects')
    ax.set_xlabel('Year')
    ax.figure.savefig('Plots/cumplot.png',bbox_inches='tight') 
    
    dfCum = (df.cumsum()).sum(axis=1)-dfMarginal
    
    fig, ax = plt.subplots(figsize=(20, 5))
    ax.bar(dfCum.index, dfCum.values)
    # ax.bar(dfMarginal.index, dfCum.values)
    plt.show()    
    
    return 

###########################################################
### mStats = DescriptiveStats(mLogRet)
def Plot_Technology(dfCounts, dfTechnology):
    """
    Purpose:
        Plot the time series of mean KWH per dollar and avg counts per census tract
        
    Inputs:
        df          Dataframe with counts per year
        
    Returns:
        
    """
    fig, ax1 = plt.subplots(figsize=(8,4))
    
    col1 = 'firebrick'
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Expected yearly kWh per dollar', color=col1)
    ax1.plot(dfTechnology.index, dfTechnology, color=col1)
    ax1.tick_params(axis='y', labelcolor=col1)
    ax1.grid()
    
    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = 'mediumblue'
    ax2.set_ylabel('Residential PV system installations', color=color)  # we already handled the x-label with ax1
    ax2.plot(dfCounts.index, dfCounts, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.savefig('Plots/technology.png', bbox_inches='tight', dpi=300)
    plt.show()
        
    return 

###########################################################
### mStats = DescriptiveStats(mLogRet)
def Simulate_data(iN, iT, dPhi1, dPhi2, dTheta1, dTheta2, dPsi1, dEta, dGamma, dBeta, vKWHpD, mAdj, dCap):
    """
    Purpose:
        Simulate dataset in order to test FE Poisson CMLE
        
    Inputs:
        iN          integer, number of cross-sectional units
        iT          integer, time periods
        dPhi_1      double, AR(1)
        dPhi_2      double, AR(2)
        dTheta1     double, N1_L1
        dTheta2     double, N1_L2
        dPsi1       double, N2_L1
        dEta        double, gamma() parameter
        dGamma      double, gamma() parameter
        dBeta       double, KWHpD parameter
        vKWHpD      iT x 1 vector of trend variable
        mAdj        iN x iN adjacency matrix
        dCap        double, upper limit of gamma distribution random draw
        
    Returns:
        df          iT x iN dataframe
        
    """
    mY = np.zeros((iT+2,iN))
    mY_N1 = np.zeros((iT+2,iN))
    mY_N2 = np.zeros((iT+2,iN))
    vA_i = np.random.gamma(dGamma,scale=dEta,size=(1,iN))
    vA_i[vA_i > dCap] = dCap
    
    
    t=2
    while t<iT+2:
        mY_N1[t-1,:] = (mAdj @ mY[t-1,:].reshape(-1,1)).reshape(1,-1)
        mY_N2[t-1,:] = (mAdj @ mAdj @mY[t-1,:].reshape(-1,1)).reshape(1,-1)
        vXB = dPhi1*mY[t-1,:] + dPhi2*mY[t-2,:] + dTheta1*mY_N1[t-1,:] + dTheta2*mY_N1[t-2,:] + dPsi1*mY_N2[t-1,:] + dBeta * vKWHpD[t-2] 
        # print(vXB.shape)
        # mY[t,:] = np.random.poisson(vA_i*np.exp(vXB))
        
        if np.max((vA_i*np.exp(vXB))) > 999999999:
            vXB = np.ones((1,iN))*-np.inf
            vA_i = np.random.gamma(dGamma,scale=dEta,size=(1,iN))
            vA_i[vA_i > dCap] = dCap
            t = 2
            print('Large value in simulation: resetting loop')
        # print(t,np.max((vA_i*np.exp(vXB))))
        mY[t,:] = st.poisson.rvs(vA_i*np.exp(vXB))
        
        t+= 1
    
    # for t in range(2,iT+2):
    #     mY_N1[t-1,:] = (mAdj @ mY[t-1,:].reshape(-1,1)).reshape(1,-1)
    #     mY_N2[t-2,:] = (mAdj @ mAdj @mY[t-2,:].reshape(-1,1)).reshape(1,-1)
    #     vXB = dPhi*mY[t-1,:] + dTheta1*mY_N1[t-1,:] + dTheta2*mY_N2[t-2,:] + dBeta * vKWHpD[t-2] # N2_L1 third term
    #     # print(vXB.shape)
    #     # mY[t,:] = np.random.poisson(vA_i*np.exp(vXB))
        
    #     if np.max((vA_i*np.exp(vXB))) > 999999999:
    #         vXB = np.ones((1,iN))*-np.inf
    #         t = 2
    #     print(t,np.max((vA_i*np.exp(vXB))))
    #     mY[t,:] = st.poisson.rvs(vA_i*np.exp(vXB))
    #     # mY[t,:] = np.random.poisson(np.exp(vXB))
    
    # plt.plot(mY.mean(axis=1))
    # plt.show()
    
    df = pd.DataFrame()
    df['N0_L0'] = mY[2:,:].flatten('F')
    df['N0_L1'] = mY[1:-1,:].flatten('F')
    df['N0_L2'] = mY[:-2,:].flatten('F')
    df['N1_L1'] = mY_N1[1:-1,:].flatten('F')
    df['N1_L2'] = mY_N1[:-2,:].flatten('F')
    df['N2_L1'] = mY_N2[1:-1,:].flatten('F')
    df['KWH'] = np.tile(vKWHpD.reshape(-1),iN)
    
    # df['Check'] = np.array([1,2,3])    
    
    return df

###########################################################
### cMaxLikelihood = MaxLikelihood(vData)
def estim_Simulated_PoissonFE_MLE(df, iT, sModel):
    """
    Purpose:
        Roughly same as earlier PoissonFE_MLE but with fewer options for faster
        estimation. See documentation of other function.
    
    """
    ###########################################################
    ### cMaxLikelihood = MaxLikelihood(vData)
    def Objective(vTheta, bForAllT=False):

        # initialize the parameter values
        # print('.', end='')
        vBeta = ParameterTransform(vTheta)
        
        vLikelihoodValues = Poisson_FE(iT, vY, mX, vBeta)

        if (bForAllT == True):
            dObjValue = vLikelihoodValues.reshape(-1,)
        else:
            dObjValue = -np.mean(vLikelihoodValues)
        
        return dObjValue
    
    ###########################################################
    ### cMaxLikelihood = MaxLikelihood(vData)
    def ComputeCovarianceMatrix(vTheta):

        # compute the inverse hessian of the average log likelihood
        mH= hessian_2sided(Objective, vTheta)
        # print(mH)
        mCov = np.linalg.inv(mH)
        mCov = (mCov +  mCov.T)/2       #  Force to be symmetric
        # compute the outer product of gradients of the average log likelihood
        mG = jacobian_2sided(Objective, vTheta, True)
        # mG = jacobian_2sided(Objective, vTheta)
        # print(mG)
        # mG = np.dot(mG.T, mG) / vData.shape[0]
        # mG = np.dot(mG, mCov)
        # mCov = np.dot(mCov, mG) / vData.shape[0]
        
        mG = np.dot(mG.T, mG) / iNT
        mG = np.dot(mG, mCov)
        mCov = np.dot(mCov, mG) / iNT

        return mCov
    
    ###############################################
    ### main()
    def ParameterTransform(vTheta):

        vTranspar = np.copy(vTheta)
        # vTranspar[1] = 1/(1+np.exp(-vTheta[1]))
        vTranspar[0] = 1/(1+np.exp(-vTheta[0]))     # changed for no intercept
        # elif sMethod == 'ZIP':
        #     vTranspar[-1] = 1/(1+np.exp(-vTheta[-1]))

        return vTranspar
    
    ###############################################
    ### main()
    def CalculateCriteria(dLL):

        dLL = -dLL
        dAIC= 2*iK - 2*dLL
        dAICc= dAIC + ((2*iK**2 + 2*iK)/(iT-iK-1))
        dBIC= iK*np.log(iT) - 2*dLL

        return (dAIC, dAICc, dBIC)
    
    # initialize starting values and return value
    cReturnValue = cMaximizeLikelihood()
    # iK = 1+ len(lPre) + len(lExog)  # intercept + predetermined + exog
    # iK = len(lPre) + len(lExog)  # predetermined + exog
    iK = df.shape[1]-1
    # iK = 4
    
    # vTheta = np.zeros(iK).reshape(-1,1)
    vTheta = np.zeros(iK).reshape(-1,1)
    
    iNT = df.shape[0]
    # dMu0 = np.mean(vData)
    
    # Order: intercept | predetermined | exogenous | [gamma, delta]    
    # exact order is given in list arguments of function
    # first predetermined coefficient is always AR(1), of which the parameter
    # is bounded between -1 and 1. When RE is used, gamma distribution parameters
    # are strictly positive
    
    # vTheta[1] = np.log(0.4/(1-0.4))                   # AR(1)
    vTheta[0] = np.log(0.2/(1-0.2))                   # AR(1) - changed for no intercept
    # if sMethod == 'RE':
    #     vTheta[-2] = np.log(0.5)
    #     vTheta[-1] = np.log(0.5)
    # elif sMethod == 'ZIP':
    #     vTheta[-1] = dTau_0

    # lFull = []
    # # lFull.extend(lY)
    # lFull.extend(lPre)
    # lFull.extend(lExog) 
    # dfOrdered = df[lFull]

    vY = df.iloc[:,0].values.astype('int')
    # mX = np.hstack((np.ones(iNT).reshape(-1,1), dfOrdered.values.astype('float')))
    mX = df.iloc[:,1:].values.astype('float') # changed for no interept
    # mX = dfOrdered.values.astype('float')

    # vTheta = np.log(vTheta)
    cReturnValue.x0 = vTheta
    cReturnValue.tx0 = ParameterTransform(vTheta)
    
    # do the optimization
    tSol = opt.minimize(Objective, vTheta, method='BFGS', options={'disp': False, 'maxiter':250})
    cReturnValue.success = tSol['success']

    # check for success and store results
    if (tSol['success'] != True):
        print("*** no true convergence: ",tSol['message'])
        vNaN = np.zeros((iK))
        cReturnValue.tx = vNaN.fill(np.NaN)
    
    else:
        cReturnValue.x = tSol['x']
        cReturnValue.tx = ParameterTransform(cReturnValue.x)
        # print(cReturnValue.tx.shape)
        # print(cReturnValue.tx0)
        # print(cReturnValue.tx)
        # cReturnValue.likelihoodvalue = -iNT * tSol['fun']
        # cReturnValue.criteria = CalculateCriteria(cReturnValue.likelihoodvalue)
        # cReturnValue.covariancematrix = ComputeCovarianceMatrix(cReturnValue.x)
        # # mJ = jacobian_2sided(ParameterTransform, cReturnValue.x, True)
        # mJ = jacobian_2sided(ParameterTransform, cReturnValue.x)
        # cReturnValue.tcovariancematrix = np.dot(mJ, np.dot(cReturnValue.covariancematrix, mJ.T))
        # cReturnValue.standarderrors = np.sqrt(np.diag(cReturnValue.covariancematrix))
        # cReturnValue.tstandarderrors = np.sqrt(np.diag(cReturnValue.tcovariancematrix))
            
    return cReturnValue.tx

###########################################################
### mStats = DescriptiveStats(mLogRet)
def RunSims(iM, iN, iT, vParamSim, dEta, dGamma, dBeta, vKWHpD, mAdj, dCap):
    """
    Purpose:
        Runs simulations for a given dgp for three models (prespecified)
        
    Inputs:
        iM          integer, number of simulations
        iN          integer, number of cross-sectional units
        iT          integer, time periods
        vParamSim   (5,) vector of dPhi1, dPhi2, dTheta1, dTheta2, dPsi1
                    Data is simulated with these parameter values
        dEta        double, gamma() parameter
        dGamma      double, gamma() parameter
        dBeta       double, KWHpD parameter
        vKWHpD      iT x 1 vector of trend variable
        mAdj        iN x iN adjacency matrix
        dCap        double, upper limit of gamma distribution random draw
        
    Returns:
        cSimRes     class, see top of file
        
    """
    dPhi1, dPhi2, dTheta1, dTheta2, dPsi1 = vParamSim
    cSimRes = cSimulations()
    cSimRes.vTheta0 = [dPhi1, dPhi2, dTheta1, dTheta2, dPsi1, dBeta]
    cSimRes.vOtherParams = [iN, iT, dEta, dGamma]
    cSimRes.NSims = iM
    mEstim_mod1 = np.zeros((iM,3))
    mEstim_mod2 = np.zeros((iM,4))
    mEstim_mod3 = np.zeros((iM,6))
    
    for m in range(iM):
        print(vParamSim)
        print('Simulation ', str(m+1),'/',str(iM))
        dfSim = Simulate_data(iN, iT, dPhi1, dPhi2, dTheta1, dTheta2, dPsi1, dEta, dGamma, dBeta, vKWHpD, mAdj, dCap)
        mEstim_mod1[m,:] = estim_Simulated_PoissonFE_MLE(dfSim[['N0_L0','N0_L1', 'N1_L1', 'KWH']], iT, 1)
        mEstim_mod2[m,:] = estim_Simulated_PoissonFE_MLE(dfSim[['N0_L0','N0_L1', 'N1_L1', 'N2_L1', 'KWH']], iT, 2)
        mEstim_mod3[m,:] = estim_Simulated_PoissonFE_MLE(dfSim[['N0_L0','N0_L1', 'N0_L2', 'N1_L1', 'N1_L2', 'N2_L1','KWH']], iT, 3)
    
    cSimRes.mTheta_mod1 = mEstim_mod1
    cSimRes.mTheta_mod2 = mEstim_mod2
    cSimRes.mTheta_mod3 = mEstim_mod3
    
    return cSimRes

###########################################################
### mStats = DescriptiveStats(mLogRet)
def SimRes_tocsv(cSimRes, iNum):
    """
    Purpose:
        saves simulation results in .csv format. 
        
    Inputs:
        cSimRes     class, see top of file
        iNum        integer, number to name file with
        
        
    """
    for i in range(3):
        # sName = 'SimRes_p'+str(iNum)+'_mod'+str(i+1)
        np.savetxt('Simulations/SimRes_p'+str(iNum)+'_mod1v3.csv',cSimRes.mTheta_mod1, delimiter=",")
        np.savetxt('Simulations/SimRes_p'+str(iNum)+'_mod2v3.csv',cSimRes.mTheta_mod2, delimiter=",")
        np.savetxt('Simulations/SimRes_p'+str(iNum)+'_mod3v3.csv',cSimRes.mTheta_mod3, delimiter=",")
    
    return 

###########################################################
### mStats = DescriptiveStats(mLogRet)
def SimTable(lParamVecs, mSimP, dBeta):
    """
    Purpose:
        Creates table with simulation results. See Overleaf for results.
        
    Inputs:
        lParamVecs      list, has integer that denotes parameter set (see main)
        mSimP           Matrix with true parameters of simulation
        dBeta           Parameter of trend
        
    Returns:
        vY          (iN*iT) x 1 vector
        mX          (iN*iT) x 4
        
    """
    vCols = ['DGP', 'Model', 'Metric', 'phi_1', 'phi_2', 'theta_1', 'theta_2', 'psi_1', 'beta']
    df = pd.DataFrame(index=range(len(lParamVecs)*6), columns=vCols)
    for i, par in enumerate(lParamVecs):
        df.iloc[i*6,1] = '(1)'
        df.iloc[2+(i*6),1] = '(2)'
        df.iloc[4+(i*6),1] = '(3)'
        sName1 = 'Simulations/SimRes_p'+str(par)+'_mod1v3.csv'
        sName2 = 'Simulations/SimRes_p'+str(par)+'_mod2v3.csv'
        sName3 = 'Simulations/SimRes_p'+str(par)+'_mod3v3.csv'
        mSimRes1 = pd.read_csv(sName1, header=None).dropna(axis=0).to_numpy()
        mSimRes2 = pd.read_csv(sName2, header=None).dropna(axis=0).to_numpy()
        mSimRes3 = pd.read_csv(sName3, header=None).dropna(axis=0).to_numpy()
        
        vTrue1 = np.append(mSimP[par-1,[0,2]],8)
        vTrue2 = np.append(mSimP[par-1,[0,2,4]],8)
        vTrue3 = np.append(mSimP[par-1,:],8)
        
        # vBias1 = (mSimRes1-vTrue1).mean(axis=0)
        df.iloc[6*i,[3,5,8]] = (mSimRes1-vTrue1).mean(axis=0)
        df.iloc[1+6*i,[3,5,8]] = np.sqrt(np.mean((mSimRes1-vTrue1)**2,axis=0))
        
        df.iloc[2+6*i,[3,5,7,8]] = (mSimRes2-vTrue2).mean(axis=0)
        df.iloc[3+6*i,[3,5,7,8]] = np.sqrt(np.mean((mSimRes2-vTrue2)**2,axis=0))
        
        df.iloc[4+6*i,[3,4,5,6,7,8]] = (mSimRes3-vTrue3).mean(axis=0)
        df.iloc[5+6*i,[3,4,5,6,7,8]] = np.sqrt(np.mean((mSimRes3-vTrue3)**2,axis=0))

    df = df.fillna(' ')    
    df = df.style.format(decimal='.', thousands=',', precision=3)
    
    
    return df


###########################################################
### main
def main():
    ############## Below is estimation for NY only ########
    ############## 'NYS' used in syntax ###################
    sNameNYS = 'Data/Solar_Electric_Programs_Reported_by_NYSERDA__Beginning_2000_20241231.csv'
    
    dfNYS, dfCS = readData(sNameNYS, 'NYS')
    # dfCheckNYS = dfNYS.iloc[:100,:]
    
    # (dfNYS['Project Cost'] < dfNYS['Total NYSERDA Incentive']).sum() # number of cost < grant
    
    mAdj, dfOrder, dfALAND = Run_adjacency()
    
    dfALAND = (pd.DataFrame(dfALAND.values.reshape(1,-1), columns=dfOrder))/1000000
    
    # Create dataset aggregated per 6 months
    dfCounts6m, dfKWHpD6m, dfGrant, dfCost6m = AggregateObservations(dfNYS, dfOrder, '2Q')
    dfCounts, dfKWHpD, dfGrant2, dfCost = AggregateObservations(dfNYS, dfOrder, 'Y')
    
    # CumPanelsVector(dfCounts6m)
    
    dfKWHpD_mean = dfKWHpD6m.mean(axis=1)
    dfCost6m_mean = dfCost6m.mean(axis=1)
    # dfKWHpD_mean = dfKWHpD_mean.replace([np.inf, -np.inf], np.nan)
    plt.plot(dfKWHpD_mean)        # between 2010-2023, output per dollar doubled
    plt.plot(dfCost6m_mean)
    Plot_Technology(dfCounts6m.mean(axis=1), dfKWHpD_mean)
    
    # dfCostStd = dfKWHpD6m.std(axis=1)
    # dfCostStd = dfCostStd.replace([np.inf, -np.inf], np.nan)
    # plt.plot(dfCostStd)        # between 2010-2023, output per dollar doubled


    #### Matrices of explanatory data
    dfMedInc = Run_ACS(dfCounts6m, 'median income', sFillIn = '09-30')
    dfAge = Run_ACS(dfCounts6m, 'age', sFillIn = '09-30')
    dfOOH = Run_ACS(dfCounts6m, 'housing', sFillIn = '09-30')   # use this with ALAND for housing density
    dfWhite, dfBlack = Run_ACS(dfCounts6m, 'race', sFillIn = '09-30')


    # Remove values of OOH that are NaN or zero
    dfOOHmeans = dfOOH.mean(axis=0)
    # dfOOHremove = dfOOHmeans[(dfOOHmeans.isnull() == 1) | (dfOOHmeans==0)==1]
    dfOOHremove = (dfOOHmeans.isnull() == 1) | ((dfOOHmeans< 10)==1)   # has ones at NaN and <5 spots

    ### Remove entries 
    dfCounts6m, dfMedInc, dfAge, dfOOH, dfWhite, dfBlack, dfOrder, mAdj, dfALAND = RemoveCTs(dfCounts6m, dfMedInc, dfAge, dfOOH, 
                                                                                             dfWhite, dfBlack, dfOrder, mAdj, dfALAND, dfOOHremove)
    
    # Trim and interpolate
    dfCounts6m_trim, dfCounts6m_0, dfCounts6m_m1, dfKWHpD, dfMedInc, dfAge, dfOOH, dfWhite, dfBlack = Interpolate_trim(dfCounts6m, dfKWHpD_mean, dfMedInc, dfAge, 
                                                                                          dfOOH, dfWhite, dfBlack, 'explanatory')
    
    dfMedInc, dfKWHpD = Adjust_Prices(dfMedInc, dfKWHpD)    # Corrects median income and technology for inflation
    # dfCounts6m_trim.name = 'Counts'
    # dfCount_OOH = (dfCounts6m / dfOOH)*1000     # per 1000 OOHs
    
    # Create AR and neighbourhood counts
    dfCount_L1, dfN1, dfN2, dfN1_L1, dfN2_L1, dfN0_L2, dfN1_L2, dfN2_L2 = Create_AR_Neighbours(dfCounts6m_trim, mAdj, dfCounts6m_0, dfCounts6m_m1)
    
    # Per 100 OOHs to get more stable estimation
    dfOOH_100 = dfOOH/100
    dfOOH_100.name = 'OOHx100'
    dfHD = CreateHousingDensity(dfOOH_100, dfALAND)
    
    # Give them a name
    dfMedInc.name = 'MedianInc'
    dfAge.name = 'Age'
    dfOOH.name = 'OOH'
    dfCounts6m_trim.name = 'N0_L0'
    dfKWHpD.name = 'KWHpDollar'
    dfWhite.name = 'White'
    dfBlack.name = 'Black'
    
    # Put in panel format
    dfPanel, dfPanel_noindex = Create_panel([dfCounts6m_trim, dfKWHpD, dfMedInc, dfAge, dfOOH_100, 
                                             dfWhite, dfBlack, dfCount_L1, dfN1, dfN2, dfN1_L1, dfN2_L1, 
                                             dfN0_L2, dfN1_L2, dfN2_L2, dfHD,])
    vXnames_full = ['MedianInc','Age', 'OOHx100', 'KWHpDollar', 'White', 'Black',
                    'N0_L1', 'N1_L1', 'N2_L1', 'Timeperiod', 'Trump dummy', 'Time x Trump']
    
    dfPanel_noindex.to_csv('Panel_noindex.csv')
    
    # Quick test
    OLS = PanelOLS(dependent=dfPanel['N0_L0'],
             exog=dfPanel[['Median income','Age','Counts_L1', 'N1_L1', 'Timeperiod', 'Trump dummy', 'Time x Trump']],
             entity_effects=False,
             time_effects=False)
    OLS.fit(cov_type='clustered', cluster_entity=True)
    
    OLS_FE = PanelOLS(dependent=dfPanel['N0_L0'],
             exog=dfPanel[['Median income','Age', '100 OOHs', 'KWH per dollar', 'White', 'Black', 'N0_L1', 'N1_L1', 'Timeperiod', 'Trump dummy']],
             entity_effects=True,
             time_effects=False)
    OLS_FE.fit(cov_type='clustered', cluster_entity=True)
    
    OLS_FE_interaction = PanelOLS(dependent=dfPanel['N0_L0'],
             exog=dfPanel[vXnames_full],
             entity_effects=True,
             time_effects=False)
    OLS_FE_interaction.fit(cov_type='clustered', cluster_entity=True)
    
    ## Arellano Bond estimator - not used, but good to compare results with
    command_str = 'N0_L0 L1.N0_L0 MedianInc OOHx100 Age N1_L1 | gmm(N0_L0, 2:3) gmm(N1_L1,1:2) | collapse'
    mydpd = regression.abond(command_str, dfPanel_noindex, ['Census_Tract', 'Date'])

    command_str = 'N0_L0 L1.N0_L0 Median_income OOHx100 Age Neighbour_count | gmm(N0_L0, 2:3) gmm(Neighbour_count,2:3) | collapse'
    mydpd2 = regression.abond(command_str, dfPanel_noindex, ['Census_Tract', 'Date'])
    
    # Another guess
    command_str = 'N0_L0 L1.N0_L0 MedianInc OOHx100 Age N1_L1 | gmm(N0_L0, 2:3) gmm(N1_L1,1:2) iv(MedianInc OOHx100 Age) | collapse'
    mydpd = regression.abond(command_str, dfPanel_noindex, ['Census_Tract', 'Date'])
    
    
    ## POISSON REGRESSIONS
    lY = ['N0_L0']
    lPre = ['N0_L1', 'N1_L1']
    # lExog = ['KWHpDollar', 'MedianInc', 'OOHx100', 'White', 'Timeperiod', 'Trump dummy', 'Time x Trump']
    lExog = ['KWHpDollar', 'MedianInc', 'OOHx100', 'White', 'Timeperiod', 'Trump dummy']
    # lExog = ['KWHpDollar']
    
    cPoissonRE = estim_Poisson_MLE(dfPanel, lY, lPre, lExog, dfCounts6m_trim.shape[1], 0, sMethod='RE')
    # cPoissonZIP_0 = estim_Poisson_MLE(dfPanel, ['N1_L1', 'Counts_L1'], lPre, lExog, dfCounts6m_trim.shape[1], 0, sMethod='ZIP')
    cPoissonZIP_0 = estim_Poisson_MLE(dfPanel, lY, lPre, lExog, dfCounts6m_trim.shape[1], 0, sMethod='ZIP')
    cPoissonZIP_1 = estim_Poisson_MLE(dfPanel, lY, lPre, lExog, dfCounts6m_trim.shape[1], 1, sMethod='ZIP')
    cPoissonZIP_neg1 = estim_Poisson_MLE(dfPanel, lY, lPre, lExog, dfCounts6m_trim.shape[1], -1, sMethod='ZIP')    

    # # no convergence as of now, will try with counts+1 (without fixing AR variables yet)
    # dfPanel['Counts+1'] = dfPanel['N0_L0']+1
    # lY2 = ['Counts+1']
    # cPoissonRE = estim_Poisson_MLE(dfPanel, lY2, lPre, lExog, dfCounts6m_trim.shape[1], sMethod='RE')
    
    # Fixed effects (has own MLE function as well)
    cPoissonFE = estim_PoissonFE_MLE(dfPanel, lY, lPre, lExog, dfCounts6m_trim.shape[1])
    lPre2 = ['N0_L1', 'N1_L1', 'N2_L1']
    # cPoissonFE2 = estim_PoissonFE_MLE(dfPanel, lY, lPre2, lExog, dfCounts6m_trim.shape[1])
    # cPoissonFE2_ones = estim_PoissonFE_MLE(dfPanel, lY, lPre2, lExog, dfCounts6m_trim.shape[1])
    
    # Manually changing FE function to get rid of intercept, can't estimate it anyway
    lExog2 = ['KWHpDollar', 'MedianInc', 'OOHx100', 'White', 'Black', 'Trump dummy']
    lPre3 = ['N0_L1', 'N0_L2', 'N1_L1', 'N1_L2', 'N2_L1']

    cPoissonFE_mod1 = estim_PoissonFE_MLE(dfPanel, lY, lPre, lExog2, dfCounts6m_trim.shape[1])
    cPoissonFE_mod2 = estim_PoissonFE_MLE(dfPanel, lY, lPre2, lExog2, dfCounts6m_trim.shape[1])
    cPoissonFE_mod3 = estim_PoissonFE_MLE(dfPanel, lY, lPre3, lExog2, dfCounts6m_trim.shape[1])
    
    # Turning covariance functions on again 
    # cPoissonFE_test = estim_PoissonFE_MLE(dfPanel, ['N0_L0'], ['N0_L1'], ['KWHpDollar'], dfCounts6m_trim.shape[1])


    # Simulation function
    # dfSim = Simulate_data(4766, 27, 0.2, 0, 0, 0.2, 0.1, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    # cPFECMLE_sim = estim_Simulated_PoissonFE_MLE(dfSim, 27)
    
    # t0 = timer()
    # cSimRes1 = RunSims(100, 4766, 27, 0.2, 0.2, 0.2, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    # t = timer() - t0
    # print(t)

    # cSimRes2 = RunSims(100, 4766, 27, 0.2, 0.2, 0.0, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    # cSimRes3 = RunSims(100, 4766, 27, 0.2, 0.0, 0.0, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)

    # cSimRes4 = RunSims(100, 4766, 27, 0.4, 0.2, 0.2, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    # cSimRes5 = RunSims(100, 4766, 27, 0.4, 0.2, 0.0, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    # cSimRes6 = RunSims(100, 4766, 27, 0.4, 0.0, 0.0, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    
    # Below are the values for the DGP
    vSimP1 = [0.2, 0, 0, 0, 0]              # phi_1
    vSimP2 = [0.2, 0, 0.2, 0, 0]            # phi_1 theta_1
    vSimP3 = [0.2, 0, 0.2, 0, 0.2]          # phi_1 theta_1 psi_1
    vSimP4 = [0.2, 0.1, 0.2, 0.1, 0.1]      # phi_1 phi_2 theta_1 theta_2 psi_1
    vSimP5 = [0.2, 0.1, 0.2, 0.1, 0]        # phi_1 phi_2 theta_1 theta_2
    vSimP6 = [0.2, 0, 0, 0.2, 0]            # phi_1 theta_2
    
    # stack them for easy use when creating table
    mSimP = np.vstack((vSimP1, vSimP2, vSimP3, vSimP4, vSimP5, vSimP6))
    
    
    # t0 = timer()
    cSimResP1 = RunSims(80, 4766, 27, vSimP1, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    cSimResP2 = RunSims(80, 4766, 27, vSimP2, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    cSimResP3 = RunSims(80, 4766, 27, vSimP3, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    cSimResP4 = RunSims(80, 4766, 27, vSimP4, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    cSimResP5 = RunSims(80, 4766, 27, vSimP5, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    cSimResP6 = RunSims(80, 4766, 27, vSimP6, 0.01, 0.2, 8, dfKWHpD.iloc[:,0].values.reshape(-1,1)/100, mAdj, 0.5)
    # t = timer() - t0
    # print(t)
    SimRes_tocsv(cSimResP1, 1)
    SimRes_tocsv(cSimResP2, 2)
    SimRes_tocsv(cSimResP3, 3)
    SimRes_tocsv(cSimResP4, 4)
    SimRes_tocsv(cSimResP5, 5)
    SimRes_tocsv(cSimResP6, 6)

    dfSimTable = SimTable([2,3,4], mSimP, 8)    # in main text
    print(dfSimTable.to_latex())

    dfSimTable2 = SimTable([1,5,6], mSimP, 8)   # in appendix
    print(dfSimTable2.to_latex())


###########################################################
### start main
if __name__ == "__main__":
    main()
