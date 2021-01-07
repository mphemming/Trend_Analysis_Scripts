

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
from random import random, sample, seed
import random
import signalz
from sklearn.metrics import mean_squared_error

# %% -----------------------------------------------------------------------------------------------
# Load data

# mooring data
main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
NRSPHB_clim = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_Climatology_1953-2019_BottleCTDMooringSatellite.nc')
NRSPHB_agg = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_1953-2019_aggregated.nc')

# surface
csurf = [(NRSPHB_agg.DEPTH >= 0) & (NRSPHB_agg.DEPTH <= 2)]
Dsurf = np.array(NRSPHB_agg.DEPTH)
Dsurf = Dsurf[csurf]
tsurf = np.array(NRSPHB_agg.TIME)
tsurf = tsurf[csurf]
Tsurf = np.array(NRSPHB_agg.TEMP)
Tsurf = Tsurf[csurf]

tbin,Tbin = TF.bin_daily(1953,2020,tsurf,Tsurf)

# De-season data
# select climatology at similar depth
climsurf = NRSPHB_clim.TEMP_AVE[:,0] # climatology at 20m
# get de-seasoned temperatures
Tbin_deseason = np.array(TF.deseason(tbin,Tbin,climsurf))

# Get monthly averages
# Using de-seasoned timeseries
tbin_m,Tbin_m = TF.bin_monthly(1953,2021,tbin,Tbin_deseason)


# determine leakage to get closest matching brownian noise signal to TEMP

# Autocorrelation analysis
check_nans = np.isfinite(Tbin_m)
T = Tbin_m[check_nans]
ACF_result = pd.Series(sm.tsa.acf(T, nlags=10))
tests = np.arange(0,1,0.02)
ACF_tests = []
RMSE_tests = []
for n in range(len(tests)):
    x = signalz.brownian_noise(len(Tbin_m), leak=tests[n], start=0, std=0.8, source="gaussian")
    ACF_tests.append(pd.Series(sm.tsa.acf(x, nlags=10)))
    A = ACF_tests[n]
    RMSE_tests.append(np.sqrt(mean_squared_error(ACF_result[0:3], A[0:3])))
leakage = np.float(tests[RMSE_tests == np.min(RMSE_tests)])

# %% -----------------------------------------------------------------------------------------------
# Create Synthetic Data
# With similar properties to real data at surface

# create a synthetic trend
# exponential trend
# t = np.exp([0, 1, 2, 3, 4, 5 ,6 ,7 ,8 ,9 ,10])/10000
# t = np.interp(np.arange(0,10000,1)+1, 
#               [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,10000], t)
# piecewise trend
t = []
t[0:400] = np.arange(0,0.2,0.2/400)
t[401:600] = np.arange(0.2,0.5,0.3/200)
t[601:800] = np.arange(0.5,1,0.5/200)

# create synthetic data
data = []
for n in range(0,50):
        sig = signalz.brownian_noise(len(t), leak=leakage, 
                    start=0, std=0.8, source="gaussian")
        data.append(sig + t)


# plt.plot(data[500])
# plt.plot(data[200])
# plt.plot(data[50])
# plt.plot(data[25])
# plt.plot(data[1])
# plt.show()


# %% -----------------------------------------------------------------------------------------------
# Test to see if EEMD works on 'monthly data'
# assuming this data is monthly

time = np.arange(0,800,1)+1
tr = []
tr_diff = []
for n in range(len(data)):
    print(n)
    TR,imfs = TF.Ensemble_EMD_quick(time,data[n])
    tr.append(np.array(TR))
    tr_diff.append(t-tr[n])


# get standard deviations for each month
d_mean = []
d_std = []
for n in range(len(time)):
    d = []
    for dn in range(len(data)):
        da = tr[dn]
        d.append(da[n])
    d_mean.append(np.nanmean(d))
    d_std.append(np.nanstd(d))
d_mean = np.array(d_mean)
d_std = np.array(d_std)

# Create plot 
for n in range(len(tr)):
    a = tr[n]
    plt.plot(time,tr[n]-a[0],color='grey')
plt.plot(time,t,color='k',linewidth=3)
plt.plot(time,d_mean,color='r',linewidth=3)
plt.plot(time,d_mean+d_std,color='r',linewidth=3,linestyle='--')
plt.plot(time,d_mean-d_std,color='r',linewidth=3,linestyle='--')
plt.show()

# # Create plot 
# for n in range(len(tr)):
#     a = tr[n]
#     plt.plot(time,np.abs(tr_diff[n]))
# plt.show()


# plt.plot(time,data[0])
# plt.plot(time,tr[0],2)
# plt.plot(time,t,color='k')
# plt.show()

# %% -----------------------------------------------------------------------------------------------
# Similar as above but with missing data points


time = np.arange(0,800,1)+1
tr_10 = []; # 10% missing data (random)
tr_10_chunk = []; # 10% missing data in chunk
tr_25 = []; # 25% missing data (random)
tr_25_chunk = []; # 25% missing data in chunk
tr_50 = []; # 50% missing data (random)
tr_50_chunk = []; # 50% missing data in chunk

for n in range(len(data)):
    print(n)
    d = data[n]
    d_10 = np.array(d); d_25 = np.array(d); d_50 = np.array(d);
    
    r_ind = random.sample(range(0, 799), 80)   
    d_10[r_ind] = np.nan;
    r_ind = random.sample(range(0, 799), 200)   
    d_25[r_ind] = np.nan;
    r_ind = random.sample(range(0, 799), 400)   
    d_50[r_ind] = np.nan;        
        
    d_10_chunk = np.array(d); d_10_chunk[100:140] = np.nan; d_10_chunk[600:640] = np.nan; 
    d_25_chunk = np.array(d); d_25_chunk[100:200] = np.nan; d_25_chunk[600:700] = np.nan; 
    d_50_chunk = np.array(d); d_50_chunk[100:300] = np.nan; d_50_chunk[600:800] = np.nan; 
    
    #-----------------------------------------------------
    print('10%')
    check = np.isfinite(d_10)
    TR,imfs = TF.Ensemble_EMD_quick(time[check],d_10[check])
    tr_10.append(np.array(TR))
    #-----------------------------------------------------
    print('25%')
    check = np.isfinite(d_25)
    TR,imfs = TF.Ensemble_EMD_quick(time[check],d_25[check])
    tr_25.append(np.array(TR))
    #-----------------------------------------------------
    print('50%')
    check = np.isfinite(d_50)
    TR,imfs = TF.Ensemble_EMD_quick(time[check],d_50[check])
    tr_50.append(np.array(TR))
    #-----------------------------------------------------
    print('10% chunk')
    check = np.isfinite(d_10_chunk)
    TR,imfs = TF.Ensemble_EMD_quick(time[check],d_10_chunk[check])
    tr_10_chunk.append(np.array(TR))
    #-----------------------------------------------------
    print('25% chunk')
    check = np.isfinite(d_25_chunk)
    TR,imfs = TF.Ensemble_EMD_quick(time[check],d_25_chunk[check])
    tr_25_chunk.append(np.array(TR))
    #-----------------------------------------------------
    print('50% chunk')
    check = np.isfinite(d_50_chunk)
    TR,imfs = TF.Ensemble_EMD_quick(time[check],d_50_chunk[check])
    tr_50_chunk.append(np.array(TR))   
    
    
    plt.plot(tr_10[1]-np.nanmin(tr_10[1]))
    plt.plot(tr_25[1]-np.nanmin(tr_25[1]))
    plt.plot(tr_50[1]-np.nanmin(tr_50[1]))
    plt.show()
    
    plt.plot(tr_10_chunk[1])
    plt.plot(tr_25_chunk[1])
    plt.plot(tr_50_chunk[1])
    plt.show()    
    

# %% -----------------------------------------------------------------------------------------------
# Test EEMD method

# time = np.arange(0,5000,1)+1
# d = data[200]

# # Daily
# tr_d,imfs = TF.Ensemble_EMD_quick(time,d)
# # every 30 days
# d_m = []
# # bin data
# time_m = np.arange(0,5000,30)+1
# for n in range(len(time_m)):
#     check = np.squeeze(np.logical_and([time >= time_m[n]-15], [time < time_m[n]+15]))
#     d_m.append(np.nanmean(d[check]))
# d_m = np.array(d_m)
# tr_m, imfs = TF.Ensemble_EMD_quick(time_m,d_m)
# # every year
# time_y = np.arange(0,5000,365)+1
# d_y = []
# for n in range(len(time_y)):
#     check = np.squeeze(np.logical_and([time >= time_y[n]-180], [time < time_y[n]+180]))
#     d_y.append(np.nanmean(d[check]))
# d_y = np.array(d_y)  
# tr_y, imfs = TF.Ensemble_EMD_quick(time_y,d_y)





# plt.plot(time,d)
# plt.plot(time_m,d_m)
# plt.show()


# plt.plot(tr_d-np.nanmin(tr_d))
# plt.plot((imfs[10]+imfs[9])-np.nanmin(imfs[10]+imfs[9]))
# plt.plot(t)
# plt.show()

# plt.plot(time,tr_d-np.nanmin(tr_d))
# plt.plot(time_y,(imfs[1]+imfs[2])-np.nanmin(imfs[1]+imfs[2]))
# # plt.plot(time_m,((imfs[6])+imfs[5]+imfs[4])-np.nanmin(imfs[6]+imfs[5]+imfs[4]))
# plt.plot(time,t)
# plt.show()


# plt.plot(time,t-(tr_d-np.nanmin(tr_d)))






