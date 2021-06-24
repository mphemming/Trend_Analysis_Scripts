
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Import modules

import os
os.chdir('C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Trend_Analysis_Scripts')
# general functions
import xarray as xr
import matplotlib as mat
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
# Import functions from script
import Trends_functions as TF
# KPSS test
from statsmodels.tsa.stattools import kpss
# autocorrelation
import statsmodels.api as sm
# 1D interpolation
from scipy.interpolate import interp1d
# wavelet analysis
import pycwt as wavelet
from pycwt.helpers import find
# For montecarlo simulation
from scipy.stats import norm
from random import seed
from random import random
import signalz
from scipy.io import savemat

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Load data

#########################################################################
# Sydney Airport
filename = ('C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\' +
            'HC06D_Data_066037_999999999979200.txt')
SYDAIR = pd.read_csv(filename,usecols=[1,2,3,4,5,6,12,14,16,18,20,22,24,26])
SYDAIR = SYDAIR.apply(pd.to_numeric, args=('coerce',))

# sort out time
yr = np.int32(SYDAIR['Year'])
mn = np.int32(SYDAIR['Month'])
dy = np.int32(SYDAIR['Day'])
hr = np.int32(SYDAIR['Hour'])
SYDAIR_TIME = []
for n in range(len(SYDAIR['Year'])):
    # month
    if mn[n] < 10:
        mn_str = '0' + str(mn[n])
    else:
        mn_str = str(mn[n])
    # day
    if dy[n] < 10:
        dy_str = '0' + str(dy[n])
    else:
        dy_str = str(dy[n])    
    # hour
    if hr[n] < 10:
        hr_str = '0' + str(hr[n])
    else:
        hr_str = str(hr[n])           
    
    date = (str(yr[n]) + '-' + mn_str + '-' + dy_str
            + ' ' + hr_str + ':00:00')
    
    SYDAIR_TIME.append(np.datetime64(date))
    
SYDAIR_TIME = np.array(SYDAIR_TIME)
_, _, _, _, SYDAIR_DOY = TF.datevec(SYDAIR_TIME)


# grab wind data
SYDAIR_WIND = np.array(SYDAIR['Wind speed measured in km/h'])
SYDAIR_WIND[SYDAIR_WIND == 0] = np.nan
# convert to m/s
SYDAIR_WIND = (SYDAIR_WIND*1000)/3600


#########################################################################
# MAI 
filename = ('C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\' +
            'HC06D_Data_092124_999999999979200.txt')
MAI = pd.read_csv(filename,usecols=[1,2,3,4,5,6,12,14,16,18,20,22,24,26])
MAI = MAI.apply(pd.to_numeric, args=('coerce',))

# sort out time
yr = np.int32(MAI['Year'])
mn = np.int32(MAI['Month'])
dy = np.int32(MAI['Day'])
hr = np.int32(MAI['Hour'])
MAI_TIME = []
for n in range(len(MAI['Year'])):
    # month
    if mn[n] < 10:
        mn_str = '0' + str(mn[n])
    else:
        mn_str = str(mn[n])
    # day
    if dy[n] < 10:
        dy_str = '0' + str(dy[n])
    else:
        dy_str = str(dy[n])    
    # hour
    if hr[n] < 10:
        hr_str = '0' + str(hr[n])
    else:
        hr_str = str(hr[n])           
    
    date = (str(yr[n]) + '-' + mn_str + '-' + dy_str
            + ' ' + hr_str + ':00:00')
    
    MAI_TIME.append(np.datetime64(date))
    
MAI_TIME = np.array(MAI_TIME)
_, _, _, _, MAI_DOY = TF.datevec(MAI_TIME)


# grab wind data
MAI_WIND = np.array(MAI['Wind speed measured in km/h'])
MAI_WIND[MAI_WIND == 0] = np.nan
# convert to m/s
MAI_WIND = (MAI_WIND*1000)/3600

# %% -------------------------------------------------------
# get de-seasoned wind anomalies for trend analysis

#########################################################################
# Sydney Airport
# calculate simple climatology
clim = TF.calc_clim_monthly(SYDAIR_TIME,SYDAIR_WIND)  
# get de-seasoned anomalies

SYDAIR_TIME_m_deseason = []
SYDAIR_WIND_m_deseason = []

SYDAIR_TIME_m,WW = TF.bin_monthly(1953,2021,SYDAIR_TIME,np.float64(SYDAIR_WIND))
SYDAIR_WIND_m_deseason,SYDAIR_WIND_m,_,SYDAIR_WIND_m_no_fill = TF.fill_gaps(
                                    SYDAIR_TIME_m,WW,np.squeeze(clim),30*12)

#########################################################################
# MAI
# calculate simple climatology
clim = TF.calc_clim_monthly(MAI_TIME,MAI_WIND)  
# get de-seasoned anomalies

MAI_TIME_m_deseason = []
MAI_WIND_m_deseason = []

MAI_TIME_m,WW = TF.bin_monthly(1953,2021,MAI_TIME,np.float64(MAI_WIND))
MAI_WIND_m_deseason,MAI_WIND_m,_,MAI_WIND_m_no_fill = TF.fill_gaps(
                                    MAI_TIME_m,WW,np.squeeze(clim),30*12)


# %% -------------------------------------------------------
# get EEMD trend

#########################################################################
# Sydney Airport
SYDAIR_t, SYDAIR_W, SYDAIR_trend, SYDAIR_trend_EAC, SYDAIR_imfs, _, _, _ = TF.Ensemble_EMD(
    SYDAIR_TIME_m,SYDAIR_WIND_m-np.nanmean(SYDAIR_WIND_m),0,0)
#########################################################################
# Sydney Airport
MAI_t, MAI_W, MAI_trend, MAI_trend_EAC, MAI_imfs, _, _, _ = TF.Ensemble_EMD(
    MAI_TIME_m,MAI_WIND_m-np.nanmean(MAI_WIND_m),0,0)


