
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
from numpy.fft import fft, ifft, fft2, ifft2, fftshift
import scipy as sp

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Load data

# SYDAP
# main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
# Data = sp.io.loadmat(main_path + 'Data\\SYDAP_TEMP.mat');
# TIME = []
# t = Data['TIME']
# for n in range(len(t)):
#     TIME.append(np.datetime64(t[n]))
# TIME = np.array(TIME)
# TEMP = np.array(Data['T'])
filename = ('C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\' +
            'HC06D_Data_066037_999999999979200.txt')
SYDAIR = pd.read_csv(filename,usecols=[1,2,3,4,5,6,12,14,16,18,20,22,24,26])
SYDAIR = SYDAIR.apply(pd.to_numeric, args=('coerce',))
TEMP = np.array(SYDAIR['Air temperature in Degrees C'])
P = np.array(SYDAIR['Mean sea level pressure in hPa'])
# sort out time
yr = np.int32(SYDAIR['Year'])
mn = np.int32(SYDAIR['Month'])
dy = np.int32(SYDAIR['Day'])
hr = np.int32(SYDAIR['Hour'])
TIME = []
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
    
    TIME.append(np.datetime64(date))
    
TIME = np.array(TIME)

# Tasman Island
filename = ('C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\' +
            'HC06D_Data_092124_999999999979200.txt')
MAI = pd.read_csv(filename,usecols=[1,2,3,4,5,6,12,14,16,18,20,22,24,26])
MAI = MAI.apply(pd.to_numeric, args=('coerce',))
MAI_TEMP = np.array(MAI['Air temperature in Degrees C'])
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


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get climatology at same depths

# calculate simple climatology for now 
clim = TF.calc_clim_monthly(TIME,TEMP) 
clim_P = TF.calc_clim_monthly(TIME,P)
clim_MAI = TF.calc_clim_monthly(MAI_TIME,MAI_TEMP)  

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Bin the data
# This is done to get a regular time grid with daily resolution


tbin, Tbin = TF.bin_monthly(1953,2020,TIME,TEMP)
_, Pbin = TF.bin_monthly(1953,2020,TIME,P)
MAI_tbin, MAI_Tbin = TF.bin_monthly(1940,2020,MAI_TIME,MAI_TEMP)
 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

Tbin_deseason = np.array(TF.deseason(tbin,Tbin,clim))
Pbin_deseason = np.array(TF.deseason(tbin,Pbin,clim_P))
MAI_Tbin_deseason = np.array(TF.deseason(MAI_tbin,MAI_Tbin,clim_MAI))
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')

t, T, trend, trend_EAC, imfs, imfs_std, imfs_to_ave, res = \
                            TF.Ensemble_EMD(tbin,Tbin_deseason,0,0)
t, P, trend_P, trend_EAC_P, imfs_P, imfs_std_P, imfs_to_ave_P, res_P = \
                            TF.Ensemble_EMD(tbin,Pbin_deseason,0,0)

MAI_t, MAI_T, MAI_trend, MAI_trend_EAC, MAI_imfs, MAI_imfs_std, MAI_imfs_to_ave, MAI_res = \
                            TF.Ensemble_EMD(MAI_tbin,MAI_Tbin_deseason,0,0)

EEMD_t_str = []
a = []  
for n in range(len(t)):
    tt = t[n]
    a.append(tt.strftime("%Y-%m-%d %H:%M:%S"))
    EEMD_t_str.append(a)


TEMP_EEMD = {'EEMD_t':EEMD_t_str,
            'EEMD_T':T,
            'EEMD_trend':trend,
            'EEMD_trend_EAC':trend_EAC,
            'EEMD_imfs':imfs}


savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "SYDAP_TEMP_analysis.mat", TEMP_EEMD)

# %% -----------------------------------------------------------------------------------------------
# Look at summer diurnal cycle, any change over time?

yr, mn, dy, hr, yday = TF.datevec(TIME)
yr_n = [1950, 1960, 1970, 1980, 1990, 2000, 2010]
hr_n = [0, 3, 6, 9, 12, 15, 18, 21, 0]

PHB_summer_diurnal = np.ones((len(yr_n),len(hr_n)),dtype=float)*np.nan
for n in range(len(yr_n)):
    for nn in range(len(hr_n)):
        c = [(yr > yr_n[n]) & (yr <= yr_n[n]+10) & (hr == hr_n[nn])]
        cmn = [(mn == 11) | (mn == 12) | (mn == 1)]
        c = np.squeeze(np.logical_and(c, cmn))  
        PHB_summer_diurnal[n,nn] = np.nanmean(TEMP[c])
        if nn == 8:
            PHB_summer_diurnal[n,nn] = PHB_summer_diurnal[n,0]
        
plt.scatter(hr_n,PHB_summer_diurnal[0,:]-np.nanmean(PHB_summer_diurnal[0,:]),c='b'); 
plt.scatter(hr_n,PHB_summer_diurnal[6,:]-np.nanmean(PHB_summer_diurnal[6,:]),c='r');
plt.ylabel('Temperature [deg C]')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



