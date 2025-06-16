# -*- coding: utf-8 -*-
"""
AggregateACS.py

Purpose:
    Aggregates ACS files into a panel setting. Idea is to call one function from
    this file that takes the solar count panel as an input to copy sizes. Need
    to think about imputation and where to place the data point (beginning or
    end of year). Need to discuss this with supervisor

Version:
    1       First start 
    2       Done for median income, but 504 missing obs from 2020 on
    3       Done for age
    4       Now stores floats instead of strings
    5       Working on filling in post 2020 values
    6       Works with cleaned difference file now, need to look into this.
            Only done for median income.
    7       (incomplete) dataframe for housing done. Changed structure to be
            able to make calculations in census data before putting it in
            the table.
    8       Changes to 'age' [fixed] and 'housing' [WIP]
    9       DO NOT USE, FUCKED UP THE HOUSING/RACE SCRIPT
    10      Housing script should be fixed, continuing race

Date:
    2025/05/15

Author:
    Dennis Cordes       681369
"""
###########################################################
### Imports
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import random
# import os
# import datetime
# import statistics
# import math
# import sympy
# import scipy.special as ss
# import scipy.stats as st
# import scipy.optimize as opt
# import scipy
# import statsmodels.api as sm
# from statsmodels.tsa.arima.model import ARIMA
# from timeit import default_timer as timer
# from collections import Counter
# import re
# import seaborn as sns
# from sklearn.metrics import roc_auc_score, roc_curve, auc

# from johansen import coint_johansen

###########################################################
### readData(sFileName)
def Run_ACS(dfOG, sType, sFillIn = '09-30'):
    """
    Purpose:  
        
        
    Inputs:
        df              DataFrame, used to copy shape and take CTs and dates
        sType           string, variable to create. Options are:
                            'median income'
                            'race'
                            'housing'
                            'age'
        sFill           string, what date to fill in value.
        
    Return values
        df              Dataframe containing many, many things
    
    """
    # df = pd.DataFrame(np.NaN, index=dfOG.index, columns=dfOG.columns).reindex_like(dfOG)
    df = pd.DataFrame(np.NaN, index=dfOG.index, columns=dfOG.columns)
    
    dfDiff = pd.read_csv('Data/Census/2020_2010_diff_cleaned.txt', sep="|")     # census differences
    vDates2020 = ['2020-09-30','2021-09-30','2022-09-30','2023-09-30']  # dates to fill in

    
    # dates = df.index()
    if sType == 'median income':        # change
        # Loop that fills in last obs of a calendar year
        sColName = 'Households!!Estimate!!Median income (dollars)'      # change
        for r in range(2010,2024):
            print(r)
            if r == 2017:
                sColName = 'Estimate!!Households!!Median income (dollars)'      # change
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/S1901_income/ACSST5Y'+str(r)+'.S1901-Data.csv'     # change
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            
            vCTs = np.intersect1d(df.columns, dfTemp.index.values)
            
            # for c in df.columns:
            for c in vCTs:
                # df[df.index==date][c]= dfTemp.loc[c,sColName]
                df.loc[date,c]= dfTemp.loc[c,sColName] 
        print((df.loc['2019-09-30']=='-').sum())
        print(df.loc['2020-09-30'].isnull().sum())
        df = df.replace('-',np.NaN)
        df = df.replace('(X)',np.NaN)
        df = df.replace('2,500-',2500)
        df = df.replace('250,000+',250000)
        df = df.astype('float64')
        
        df2020 = df.loc[vDates2020]                                     # From here is new in loop (v5)
        vToChange = df2020.loc[:, df2020.isnull().all()].columns
        
        for r in range(2020,2024):          # TEMP CHANGE TO 2023
            print(r)
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/S1901_income/ACSST5Y'+str(r)+'.S1901-Data.csv'     # change
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            dfTemp = dfTemp.replace('-',np.NaN)
            dfTemp = dfTemp.replace('(X)',np.NaN)
            dfTemp = dfTemp.replace('2,500-',2500)
            dfTemp = dfTemp.replace('250,000+',250000)

            
            # vCTs = np.intersect1d(vToChange, dfTemp.index.values)
            # print(vCTs.shape)
            
            for c in vToChange:
                # dfDiff[dfDiff['GEOID_TRACT_10']==36123150500]['GEOID_TRACT_20'].values
                df.loc[date,c] = (dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,'Estimate!!Households!!Median income (dollars)'].astype('float64').mean())

    elif sType == 'age':
        # Loop that fills in last obs of a calendar year
        sColName = 'Estimate!!SEX AND AGE!!Median age (years)'
        for r in range(2010,2024):
            print(r)
            if r == 2017:
                sColName = 'Estimate!!SEX AND AGE!!Total population!!Median age (years)'
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/DP05_NY_age/ACSDP5Y'+str(r)+'.DP05-Data.csv'
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            
            vCTs = np.intersect1d(df.columns, dfTemp.index.values)
            
            # for c in df.columns:
            for c in vCTs:
                # df[df.index==date][c]= dfTemp.loc[c,sColName]
                df.loc[date,c]= dfTemp.loc[c,sColName]
        print((df.loc['2019-09-30']=='-').sum())
        print(df.loc['2020-09-30'].isnull().sum())
        df = df.replace('-',np.NaN)
        df = df.astype('float64')
        
        df2020 = df.loc[vDates2020]                                     # From here is new in loop (v5)
        vToChange = df2020.loc[:, df2020.isnull().all()].columns
        
        for r in range(2020,2024):          # TEMP CHANGE TO 2023
            print(r)
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/DP05_NY_age/ACSDP5Y'+str(r)+'.DP05-Data.csv'     # change
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            dfTemp = dfTemp.replace('-',np.NaN)
            
            # vCTs = np.intersect1d(vToChange, dfTemp.index.values)
            # print(vCTs.shape)
            
            for c in vToChange:
                # dfDiff[dfDiff['GEOID_TRACT_10']==36123150500]['GEOID_TRACT_20'].values
                df.loc[date,c] = dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,'Estimate!!SEX AND AGE!!Total population!!Median age (years)'].astype('float64').mean()

        
    elif sType == 'housing':
        # Loop that fills in last obs of a calendar year
        sColOOH = 'Total!!Estimate!!HOUSING TENURE!!Owner-occupied housing units'
        sColTH = 'Total!!Estimate!!Total households'
        for r in range(2010,2024):
            print(r)
            if r == 2017:
                sColOOH = 'Estimate!!Total!!Total households!!HOUSING TENURE!!Owner-occupied housing units'
                sColTH = 'Estimate!!Total!!HOUSEHOLDS!!Total households'
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/S1101_NY_housing/ACSST5Y'+str(r)+'.S1101-Data.csv'
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            dfTemp = dfTemp.replace('-',np.NaN)

            vCTs = np.intersect1d(df.columns, dfTemp.index.values)
            
            # for c in df.columns:
            for c in vCTs:
                # df[df.index==date][c]= dfTemp.loc[c,sColName]
                df.loc[date,c]= (float(dfTemp.loc[c,sColOOH]) * float(dfTemp.loc[c,sColTH]))/100
        print((df.loc['2019-09-30']=='-').sum())
        print(df.loc['2020-09-30'].isnull().sum())
        df = df.replace('-',np.NaN)
        # df = df.astype('float64')
    
        df2020 = df.loc[vDates2020]                                     # From here is new in loop (v5)
        vToChange = df2020.loc[:, df2020.isnull().all()].columns
        
        for r in range(2020,2024):          # TEMP CHANGE TO 2023
            print(r)
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/S1101_NY_housing/ACSST5Y'+str(r)+'.S1101-Data.csv'     # change
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            dfTemp = dfTemp.replace('-',np.NaN)
            
            # vCTs = np.intersect1d(vToChange, dfTemp.index.values)
            # print(vCTs.shape)
            
            for c in vToChange:
                # dfDiff[dfDiff['GEOID_TRACT_10']==36123150500]['GEOID_TRACT_20'].values
                # df.loc[date,c] = dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,'Estimate!!SEX AND AGE!!Total population!!Median age (years)'].astype('float64').mean()
                df.loc[date,c] = ((dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColOOH].astype('float64')*dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColTH].astype('float64'))/100).mean()

    elif sType == 'race':
        dfB = pd.DataFrame(np.NaN, index=dfOG.index, columns=dfOG.columns) # Need a second df for black pct
        # Loop that fills in last obs of a calendar year
        sColTotal = 'Estimate!!Total'
        sColWhite = 'Estimate!!Total!!White alone'
        sColBlack = 'Estimate!!Total!!Black or African American alone'
        for r in range(2010,2024):
            print(r)
            if r == 2019:
                sColTotal = 'Estimate!!Total:'
                sColWhite = 'Estimate!!Total:!!White alone'
                sColBlack = 'Estimate!!Total:!!Black or African American alone'
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/B02001_NY_race/ACSDT5Y'+str(r)+'.B02001-Data.csv'
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            dfTemp = dfTemp.replace('-',np.NaN)

            vCTs = np.intersect1d(df.columns, dfTemp.index.values)
            
            # for c in df.columns:
            for c in vCTs:
                # df[df.index==date][c]= dfTemp.loc[c,sColName]
                if float(dfTemp.loc[c,sColTotal]) != 0:
                    df.loc[date,c]= (float(dfTemp.loc[c,sColWhite]) / float(dfTemp.loc[c,sColTotal]))*100
                    dfB.loc[date,c]= (float(dfTemp.loc[c,sColBlack]) / float(dfTemp.loc[c,sColTotal]))*100
        print((df.loc['2019-09-30']=='-').sum())
        print(df.loc['2020-09-30'].isnull().sum())
        df = df.replace('-',np.NaN)
        dfB = dfB.replace('-',np.NaN)
        # df = df.astype('float64')
    
        df2020 = df.loc[vDates2020]                                     # From here is new in loop (v5)
        vToChange = df2020.loc[:, df2020.isnull().all()].columns
        
        for r in range(2020,2024):          # TEMP CHANGE TO 2023
            print(r)
            date = str(r)+'-'+sFillIn
            sFileName = 'Data/Census/B02001_NY_race/ACSDT5Y'+str(r)+'.B02001-Data.csv'     # change
            dfTemp = pd.read_csv(sFileName, header=1)
            # dfTemp['Census Tract'] = dfTemp['Geography'].str.strip().str[-11:]
            dfTemp['Census Tract'] = (dfTemp['Geography'].str[-11:]).astype('int64')
            dfTemp = dfTemp.set_index('Census Tract')
            dfTemp = dfTemp.replace('-',np.NaN)
            
            # vCTs = np.intersect1d(vToChange, dfTemp.index.values)
            # print(vCTs.shape)
            
            for c in vToChange:
                if (dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColTotal].astype('float64')).sum()!=0:
                # if float(dfTemp.loc[c,sColTotal]) != 0:
                # dfDiff[dfDiff['GEOID_TRACT_10']==36123150500]['GEOID_TRACT_20'].values
                # df.loc[date,c] = dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,'Estimate!!SEX AND AGE!!Total population!!Median age (years)'].astype('float64').mean()
                    df.loc[date,c] = ((dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColWhite].astype('float64')).sum() / (dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColTotal].astype('float64')).sum())*100
                    dfB.loc[date,c] = ((dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColBlack].astype('float64')).sum() / (dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==c]['GEOID_TRACT_20'].values,sColTotal].astype('float64')).sum())*100

        return df, dfB

    return df


###########################################################
### readData(sFileName)
def Census_differences(df, sType):
    """
    Purpose:  
        First four lines are the necessary syntax, is now in the big function
        as well. The next two are the necessary code for the loop, applied to a 
        specific tract. The last two are again the same code, but now for a census
        tract that gave an error. In total, three errors came up and the accompanying
        lines were removed. Instead of 2010_2020_diff, it now uses 2010_2020_diff_cleaned.
        Need to look into this as this is not proper work. I did it to first make sure
        that the other parts of the code worked.
        
    Inputs:
        df              DataFrame, used to copy shape and take CTs and dates
        sType           string, variable to create. Options are:
                            'median income'
                            'race'
                            'housing'
                            'age'
        
    Return values
        df              Dataframe containing many, many things
    
    """
    dfDiff = pd.read_csv('Data/Census/2020_2010_diff.txt', sep="|")
    vDates2020 = ['2020-09-30','2021-09-30','2022-09-30','2023-09-30']
    df2020 = df.loc[vDates2020]
    vToChange = df2020.loc[:, df2020.isnull().all()].columns
    
    dfDiff[dfDiff['GEOID_TRACT_10']==36123150500]['GEOID_TRACT_20'].values
    dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==36123150500]['GEOID_TRACT_20'].values,'Estimate!!Households!!Median income (dollars)'].astype('float64').mean()
    
    dfDiff[dfDiff['GEOID_TRACT_10']==36103122403]['GEOID_TRACT_20'].values
    dfTemp.loc[dfDiff[dfDiff['GEOID_TRACT_10']==36103122403]['GEOID_TRACT_20'].values,'Estimate!!Households!!Median income (dollars)'].astype('float64').mean()

    
    
    dfMedInc.loc['2019-09-30',36123150500]
    
    return 

