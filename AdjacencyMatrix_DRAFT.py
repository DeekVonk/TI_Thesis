# -*- coding: utf-8 -*-
"""
AdjacencyMatrix.py

Purpose:
    Takes a .csv output file with neigbhour relationships from QGIS and creates
    an Adjacency Matrix.

Version:
    1       First start 
    2       Standard Adjacency matrix done
    3       Docstring updated
    4       Created new function to run in main TI_Thesis code
    5       Removes tracts with zero land area
    6       Now has weights instead of counts
    7       Experts ALAND

Date:
    2025/05/25

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


## Set working directory
# for Windows, run next line. For Mac, do it manually (how does it even work??):
# os.chdir('C:/Users/corde/Tinbergen/Thesis')


### Get hessian and related functions
# from python_scripts.grad import *

###########################################################
### readData(sFileName)
def readData(sFileName):
    """
    Purpose:  
        Loading in data for .csv file. Sorts by GEOID (ascending). Removes areas
        that do not have any land area.
        
    Inputs:
        sFileName       String, contains name of .csv file
        
    Return values
        df              Dataframe containing GEOID, NEIGHBOURS and original index
    
    """
    df = pd.read_csv(sFileName)
    df = df[df['ALAND']!=0]
    df = df[['GEOID', 'NEIGHBOURS', 'ALAND']].sort_values("GEOID")
    df = df.reset_index()
    # df = df[['GEOID', 'NEIGHBOURS']]

    return df

###########################################################
### readData(sFileName)
def Make_dict(df):
    """
    Purpose:  
        Creates dictionary for each GEOID with its neigbhours.
        
    Inputs:
        df              Dataframe containing GEOID, NEIGHBOURS and original index
        
    Return values
        dictNb          Dictionary of neigbhour relations
    
    """
    lNb = []
    dictNb = {} 
    
    for i in range(df.shape[0]):
        lNb.append(int(df['NEIGHBOURS'][i][:11]))
        for r in range(int((len(df['NEIGHBOURS'][i])-11)/13)):
            lNb.append(int(df['NEIGHBOURS'][i][(13+13*r):(11+13*(r+1))]))
        dictNb[df['GEOID'][i]] = lNb
        lNb = []
        
    return dictNb

###########################################################
### readData(sFileName)
def Fill_adjacency(df, dictNb, bDiagOne = True):
    """
    Purpose:  
        Creates an adjacency matrix: ones in the cells where the corresponding
        elements share a border. Boolean controls whether the diagonal element should
        be zero or one. This is important for later estimation as otherwise the
        autoregressive component shares its effect with the spatial lag term. 
        However, there are terms without neighbours which would yield a 
        reduced rank matrix if the diagonal term is zero. If estimation yields
        problems, should carefully go over this again.
        
    Inputs:
        df              Dataframe containing GEOID, NEIGHBOURS and original index
        dictNb          Dictionary of neigbhour relations
        bDiagOne        Boolean, true for ones on the diagonal, false for zeroes
        
    Return values
        df              Dataframe containing many, many things
    
    """
    dfAdj = pd.DataFrame(0, index= df['GEOID'], columns=df['GEOID'])
    
    for tract in df['GEOID']:
        for nb in dictNb[tract]:
            dfAdj.loc[tract,nb] = 1
    
    if bDiagOne == False:
        np.fill_diagonal(dfAdj.values, 0)   # Sets diagonal values to zero
    
    # print(dfAdj.shape)
    dfAdj = dfAdj.iloc[:dfAdj.shape[0], :dfAdj.shape[0]]
    dfAdj = dfAdj / np.tile(dfAdj.sum(axis=1).to_numpy().reshape(-1,1),dfAdj.shape[1])
    
    
    mAdj = dfAdj.to_numpy()
    dfOrder = df['GEOID']   
    dfALAND = df['ALAND']

    return mAdj, dfOrder, dfALAND

###########################################################
### readData(sFileName)
def Run_adjacency():
    """
    Purpose:  
        Creates an adjacency matrix: ones in the cells where the corresponding
        elements share a border. Boolean controls whether the diagonal element should
        be zero or one. This is important for later estimation as otherwise the
        autoregressive component shares its effect with the spatial lag term. 
        However, there are terms without neighbours which would yield a 
        reduced rank matrix if the diagonal term is zero. If estimation yields
        problems, should carefully go over this again.
        
    Inputs:
        df              Dataframe containing GEOID, NEIGHBOURS and original index
        dictNb          Dictionary of neigbhour relations
        bDiagOne        Boolean, true for ones on the diagonal, false for zeroes
        
    Return values
        df              Dataframe containing many, many things
    
    """
    sName = 'GIS/NY_Neighbours.csv' # Windows
    
    # Initialisation
    dfAdj = readData(sName)

    # Estimation
    dictNb = Make_dict(dfAdj)
    mAdj, dfOrder, dfALAND = Fill_adjacency(dfAdj, dictNb, True)
    dfALAND.index = dfOrder

    return mAdj, dfOrder, dfALAND

###########################################################
### main
def main():
    # Magic numbers
    # sName = 'Tracking the Sun 2024 Data File/TTS_LBNL_public_file_21-Aug-2024_all.csv' # Mac
    sName = 'GIS/NY_Neighbours.csv' # Windows
    
    # Initialisation
    dfAdj = readData(sName)

    # Estimation
    dictNb = Make_dict(dfAdj)
    mAdj, dfOrder = Fill_adjacency(dfAdj, dictNb, True)
    

    
###########################################################
### start main
if __name__ == "__main__":
    main()
