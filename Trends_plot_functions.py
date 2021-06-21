# %% -----------------------------------------------------------------------------------------------
# Import modules

# general functions
import xarray as xr
# import matplotlib as mat
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import pandas as pd
# For mann-kendall test and innovative Sen slope analysis
# https://pypi.org/project/pymannkendall/
import pymannkendall as mk
# For pettitt test - need to check results as function copied from stackoverflow
# import pettitt as pett
import pyhomogeneity as hg
# multiple regression
from sklearn import linear_model
import scipy as sp
from statsmodels.tsa.stattools import kpss
# EEMD 
from PyEMD import EMD, EEMD, Visualisation
# significance stuff
from sklearn.metrics import mean_squared_error
import signalz
import statsmodels.api as sm
import time
import random
import seaborn as sns

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Plot functions

# STD vs time vd depth plot

def get_variability_site(TIME,DEPTH,TEMP,DEPTH_TARGETS):
    
    DEPTH = np.array(DEPTH); TIME = np.array(TIME);
    TEMP = np.array(TEMP)
    
    # define timerange
    trange = xr.cftime_range(start='1950-01-01', end='2020-01-01', freq='Y')
    trange = trange.to_datetimeindex()    
    trangemid = trange[1:-1:2]
    trange1 = trange[0:-1:2]
    trange2 = trange[2:-1:2]

    # define lists
    var_std = np.ones((len(trangemid),len(DEPTH_TARGETS)),dtype=float)
    t_std = []    
    
    for n_depth in range(len(DEPTH_TARGETS)):
        for n_t in range(len(trangemid)):
            # select data
            c = [(DEPTH >= DEPTH_TARGETS[n_depth] - 3) & 
                 (DEPTH <= DEPTH_TARGETS[n_depth] + 3) &
                 (TIME >= trange1[n_t]) & 
                 (TIME <= trange2[n_t])]   
            print(np.sum(c))
            # get standard deviations
            if np.sum(c) > 20:
                var_std[n_t,n_depth] = np.nanstd(TEMP[c])
            else:
                var_std[n_t,n_depth] = np.nan
                
    t_std = trangemid

    # return standard deviations and time
    return t_std, var_std
               
        
        
        
        