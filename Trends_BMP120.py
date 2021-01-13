
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

# mooring data
main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
BMP120_agg = xr.open_dataset(main_path + 
    'Data\\IMOS_ANMN-NSW_TZ_20110329_BMP120_FV01_TEMP-aggregated-timeseries_END-20200826_C-20201207.nc')

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

# code to check data distribution
# check = np.isfinite(BMP120_agg.TEMP) 
# %matplotlib qt
# plt.hist(BMP120_agg.DEPTH[check], bins = np.arange(0,120,1))
# plt.xlim(left=0, right=120)

print('Selecting data at different depths:')

depths = [18.5, 27.5, 35, 43, 50.5, 58.5, 67, 75, 83.5, 91.5, 99.5, 107.5]

D = []
t = []
T = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # index check
    c = [(BMP120_agg.DEPTH >= depths[n] - 2) & (BMP120_agg.DEPTH <= depths[n] + 2) & \
        (BMP120_agg.TEMP_quality_control > 0) & (BMP120_agg.TEMP_quality_control <3)]
    # Depth
    d = np.array(BMP120_agg.DEPTH);
    D.append(d[c])
    # time
    tt = np.array(BMP120_agg.TIME);
    t.append(tt[c])    
    # Temp
    TT = np.array(BMP120_agg.TEMP);
    T.append(TT[c])       

# plt.plot(t[0],T[0])
# plt.plot(t[10],T[10])
# plt.show()

del c, d, tt, TT

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# QC data
# n=0

# TO DO LATER - DOESN"T SEEM TO BAD TO BE HONEST

# TQC = T[0]; tQC = t[0];
# c1 = [(TQC < 20) & (tQC > np.datetime64('2010-10-29')) & (tQC < np.datetime64('2011-04-09'))]
# c2 = [(TQC < 17)]
# c3 = [(TQC < 19.4) & (tQC > np.datetime64('2016-08-15')) & (tQC < np.datetime64('2016-08-19'))]
# c4 = [(TQC < 18.06) & (tQC > np.datetime64('2016-11-05')) & (tQC < np.datetime64('2016-11-09'))]
# TQC[c1] = np.nan; TQC[c2] = np.nan; TQC[c3] = np.nan; TQC[c4] = np.nan;  
# T[0] = TQC; 
# # n=1
# TQC = T[1]; tQC = t[1];
# c1 = [(TQC < 19.6) & (tQC > np.datetime64('2010-10-29')) & (tQC < np.datetime64('2010-11-21'))]
# c2 = [(TQC < 17) & (tQC > np.datetime64('2010-12-01')) & (tQC < np.datetime64('2011-01-01'))]
# c3 = [(TQC < 20) & (tQC > np.datetime64('2011-02-01')) & (tQC < np.datetime64('2011-05-01'))]
# c4 = [(TQC < 22) & (tQC > np.datetime64('2011-04-09')) & (tQC < np.datetime64('2011-04-13'))]
# TQC[c1] = np.nan; TQC[c2] = np.nan; TQC[c3] = np.nan; TQC[c4] = np.nan;  
# T[1] = TQC; 
# # n=2
# TQC = T[2]; tQC = t[2];
# c1 = [(TQC < 16.5) & (tQC > np.datetime64('2010-12-01')) & (tQC < np.datetime64('2010-12-25'))]
# TQC[c1] = np.nan; 
# T[2] = TQC; 
# # n = 5
# TQC = T[5]; tQC = t[5];
# c1 = [(TQC < 14) & (tQC > np.datetime64('2010-12-01')) & (tQC < np.datetime64('2011-01-01'))]
# TQC[c1] = np.nan; 
# T[5] = TQC; 


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

print('Removing the season')

# select climatology at similar depth
clim = []
for n in range(len(depths)):
    clim.append(TF.calc_clim_monthly(t[n],T[n]))
# get de-seasoned temperatures
Tbin_deseason = []
for n in range(len(depths)):
    Tbin_deseason.append(np.array(TF.deseason(t[n],T[n],clim[n])))
    
del n
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get monthly averages

# print('Getting Monthly Averages')

# # Using de-seasoned timeseries
# tbin_m = []
# Tbin_m = []
# for n in range(len(depths)):
#     print(str(depths[n]) + ' m')
#     tt,TT = TF.bin_monthly(2011,2021,t[n],Tbin_deseason[n])
#     tbin_m.append(tt)
#     Tbin_m.append(TT)
    
# del tt, TT, n
    
# plt.plot(t10m,Tbin_deseason)
# plt.plot(tbin_m,Tbin_m)
# plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get daily averages

print('Getting Daily Averages')

# Using de-seasoned timeseries
tbin = []
Tbin = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # This is done to get a regular time grid with daily resolution
    tt,TT = TF.bin_daily(2011,2021,t[n],np.float64(Tbin_deseason[n]))
    tbin.append(tt)
    Tbin.append(TT)
    
del tt, TT, n    

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# KPSS test to check for stationarity
# If the result = 'Not stationary', a deterministc trend / linear regression is not suitable

print('Checking for stationarity')
KPSS_result = []
stationarity_array = []
pval_array = []
for n in range(len(depths)):
    KPSS_result.append(TF.kpss_test((Tbin[n]))) 
    a = KPSS_result[n]
    stationarity_array.append(str(depths[n]) + ' m :  ' + a.KPSS_result)       
    pval_array.append(str(depths[n]) + ' m :  ' + str(a.KPSS_p_value))      
del a, n
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Mann kendall tests
print('Estimating Sen slopes and performing Mann Kendall tests')
mk_result = []
mk_trend = []
mk_trend_per_decade = []
mk_pval = []
for n in range(len(depths)):
    mk_result.append(mk.trend_free_pre_whitening_modification_test(Tbin[n]))
    mk_pval.append(mk_result[n].p)
    mk_trend.append(range(len(tbin[n]))*mk_result[n].slope + mk_result[n].intercept)
    tr = range(0,3652)*mk_result[n].slope + mk_result[n].intercept
    mk_trend_per_decade.append(tr[-1]-tr[0])
    
del n, tr
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Innovative trend analysis

ITA_stats = []
ITA_significance = []
ITA_slope_per_decade = []
for n in range(len(depths)):
    ITA_stats.append(TF.ITA(tbin[n],Tbin[n],0,0))
    a = ITA_stats[n]
    ITA_significance.append(a.ITA_significance)
    ITA_slope_per_decade.append(a.ITA_trend_sen_per_decade)
    
plt.plot(ITA_slope_per_decade,depths)
plt.plot(mk_trend_per_decade,depths)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD (Here for testing only)
# print('Running Ensemble EMD')
# t, T, trend, imfs, res = TF.Ensemble_EMD(tbin[10],Tbin[10])

# # Autocorrelation analysis
# check = np.where(np.isnan(Tbin_m))
# ACF1 = pd.Series(sm.tsa.acf(Tbin_m[20:401], nlags=10)); # where longest streak without nans
# ACF2 = pd.Series(sm.tsa.acf(Tbin_m[671:777], nlags=10)); # where longest streak without 
# ACF_result = []
# for n in range(0,10):
#     ACF_result.append(np.nanmean([ACF1[n],ACF2[n]]))
# ACF_result = np.array(ACF_result)

# # significance
# conf_std_limit, std_array, trend_sims, x_sims =  TF.EEMD_significance(tbin_m,Tbin_m,ACF_result,200)

# # Create figure
# plt.figure(figsize=(15,8))
# for n in range(1,200):
#     tt = trend_sims[n]
#     plt.plot(tbin_m,tt-tt[0],color='grey')
# plt.plot(tbin_m,conf_std_limit,color='r')
# plt.plot(tbin_m,conf_std_limit*-1,color='r')
# plt.plot(t,trend-trend[0],color='k',linewidth=2)
# plt.xlabel('Year')
# plt.ylabel(r'$\rmTemperature Trend [^\circ C]$')
# plt.show()

    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# trend comparison plot

# plt.plot(tbin_m,Tbin_m)
# plt.plot(tbin_m,mk_trend,'r')
# # plt.plot(t,trend,'g')
# plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Save data as mat file

# convert time to string
tbin_str = []
tbin_deseason_str = []
a = []
b = []
for nn in range(len(tbin)):
    ttt = tbin[nn]
    for n in range(len(ttt)):
        tt = ttt[n]
        if 'datetime64' in str(type(tt)):
            a.append(str(tt))
        else:
            a.append(tt.strftime("%Y-%m-%d %H:%M:%S"))
        if nn == 0:    
            tbin_str.append(a)
        
    yr, mn, dy, hr, yday = TF.datevec(ttt)
    for n in range(len(yr)):
        d = dt.datetime(yr[n],mn[n],dy[n],hr[n])
        b.append(d.strftime("%Y-%m-%d %H:%M:%S"))
    tbin_deseason_str.append(b)


Trend_dict = {'MK_result': mk_result,
'MK_trend': mk_trend,
'MK_trend_per_decade': mk_trend_per_decade,
'MK_pval': mk_pval,
'KPSS_results': KPSS_result,
'ITA_stats': ITA_stats,
'ITA_significance': ITA_significance,
'ITA_trend_per_decade': ITA_slope_per_decade}

Data_dict = {'tbin': tbin_str,
'Tbin': Tbin,
't': tbin_deseason_str,
'T': T,
'D': D,
'Tbin_deseason': Tbin_deseason,
'clims': clim,
'BMP120_agg': BMP120_agg}

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "BMP120_trends.mat", Trend_dict)

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "BMP120_data.mat", Data_dict)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


