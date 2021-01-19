

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
NRSMAI_clim = xr.open_dataset(main_path + 'Data\\Maria_Island_TEMP_Climatology_1945-2020_BottleCTDMooring.nc')
NRSMAI_agg = xr.open_dataset(main_path + 'Data\\Maria_Island_TEMP_1945-2020_aggregated.nc')
# PDO data
# PDO = pd.read_csv(main_path + 'Data\\PDO.csv',index_col=0)
# PDO_IND = np.array(PDO.ind)
# PDO_IND_t = np.array(PDO.index)
# c = PDO_IND_t > 1952
# PDO_IND = PDO_IND[c]
# PDO_IND_t = PDO_IND_t[c]
# SAM data
# SAM = pd.read_csv(main_path + 'Data\\SAM.csv',index_col=0) 
# need to add code here to sort out SAM

del main_path

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

print('Selecting data at different depths:')

depths = [2, 21]

D = []
t = []
T = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # index check
    c = [(NRSMAI_agg.DEPTH >= depths[n] - 2) & (NRSMAI_agg.DEPTH <= depths[n] + 2)]
    # Depth
    d = np.array(NRSMAI_agg.DEPTH);
    D.append(d[c])
    # time
    tt = np.array(NRSMAI_agg.TIME);
    t.append(tt[c])    
    # Temp
    TT = np.array(NRSMAI_agg.TEMP);
    T.append(TT[c])       

# plt.plot(t[0],T[0])
# plt.plot(t[10],T[10])
# plt.show()

del c, d, tt, TT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get climatology at same depths


clim = np.ones((365,2),dtype=float)
for day in range(0,365):
    day_temps = NRSMAI_clim.TEMP_AVE[day,:] 
    clim[day,:] = np.interp(depths,NRSMAI_clim.DEPTH,day_temps)

# plt.plot(NRSMAI_clim.TEMP_AVE[:,1],'k')
# plt.plot(np.arange(0,364,1),clim[:,1],'r')
# plt.plot(NRSMAI_clim.TEMP_AVE[:,2],'k')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Bin the data
# This is done to get a regular time grid with daily resolution

tbin = []
Tbin = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    tt, TT = TF.bin_daily(1945,2020,t[n],T[n])
    tbin.append(tt)
    Tbin.append(TT)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

print('Removing the season')


# get de-seasoned temperatures
Tbin_deseason = []
for n in range(len(depths)):
    cl = clim[:,n]
    Tbin_deseason.append(np.array(TF.deseason(tbin[n],Tbin[n],cl)))
    
del n, cl
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get monthly averages

# Using de-seasoned timeseries
tbin_m = []
Tbin_m = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    tt,TT = TF.bin_monthly(1945,2020,tbin[n],Tbin_deseason[n])
    tbin_m.append(tt)
    Tbin_m.append(TT)
    
   
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')

EEMD_t = []
EEMD_T = []
EEMD_trend = []
EEMD_imfs = []
EEMD_res = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    t, T, trend, imfs, res = TF.Ensemble_EMD(tbin_m[n],Tbin_m[n],0)
    EEMD_t.append(t)
    EEMD_T.append(T)
    EEMD_trend.append(trend)
    EEMD_imfs.append(imfs)
    EEMD_res.append(res)
    
EEMD_IMFS = {'IMF_1':EEMD_imfs[0],
             'IMF_2':EEMD_imfs[1]}    
    
    
    
plt.plot(EEMD_t[0],EEMD_trend[0])
plt.plot(EEMD_t[1],EEMD_trend[1])


plt.plot(EEMD_t[0],EEMD_T[0],'.')
plt.plot(EEMD_t[0],EEMD_trend[0])

    

# Autocorrelation analysis and significance
print('Running autocorrelation analysis')
# Using last 10 years only

ACF_result = []
conf_std_limit = []
std_array = []
trend_sims = []
x_sims = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m') 
    check = np.where(np.logical_and([tbin_m[n] > dt.datetime(2010,1,1)], 
                   [tbin_m[n] < dt.datetime(2020,1,1)]))      
    TT = Tbin_m[n]
    TT = TT[check[1]]
    TT = TT[np.isfinite(TT)]
    # gaps = np.where(np.diff(check) > np.array(1))
    # f_times = np.where(np.diff(gaps[1]) > 10)
    # if np.max(gaps) < len(tbin_m[n])-10:
    #     if '[]' not in str(f_times):
    #         f_times = np.array(f_times)
    #         f_times[-1] = np.int64(len(tbin_m[n]))
    #     else:
    #         f_times = np.int64(len(tbin_m[n]))          
    # a = gaps[1]
    # f_times = a[f_times]
    # ACF = []
    # for ind in range(len(f_times)-1):
    #     TT = Tbin_m[n]
    #     TT = TT[f_times[ind]+1:f_times[ind+1]]
    #     TT = TT[np.where(np.isfinite(TT))]
    #     ACF.append(pd.Series(sm.tsa.acf(TT, nlags=10)));

    ACF_result.append(np.array(pd.Series(sm.tsa.acf(TT, nlags=10))))

    # significance
    csl, sa, ts, xs = \
           TF.EEMD_significance(tbin_m[n],Tbin_m[n],ACF_result[n],1000)
    conf_std_limit.append(csl)
    std_array.append(sa)
    trend_sims.append(ts)
    x_sims.append(xs)

del TT, n, check, csl, sa, ts, xs

# Create figure
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

# %% -----------------------------------------------------------------------------------------------
# Testing edge effects


# n = 0
# tt = tbin_m[n]
# TT = Tbin_m[n]
# t_0, T_0, trend_0, imfs_0, res_0 = TF.Ensemble_EMD(tt,TT,0)
# t_1, T_1, trend_1, imfs_1, res_1 = TF.Ensemble_EMD(
#     tt[np.int64(len(tt)*0.1):len(tt)-np.int64(len(tt)*0.1)],
#     TT[np.int64(len(tt)*0.1):len(tt)-np.int64(len(tt)*0.1)],0)
# t_2, T_2, trend_2, imfs_2, res_2 = TF.Ensemble_EMD(
#     tt[np.int64(len(tt)*0.25):len(tt)-np.int64(len(tt)*0.25)],
#     TT[np.int64(len(tt)*0.25):len(tt)-np.int64(len(tt)*0.25)],0)
# t_3, T_3, trend_3, imfs_3, res_3 = TF.Ensemble_EMD(
#     tt[np.int64(len(tt)*0.4):len(tt)-np.int64(len(tt)*0.4)],
#     TT[np.int64(len(tt)*0.4):len(tt)-np.int64(len(tt)*0.4)],0)
    
# plt.plot(t_0[1::],np.diff(trend_0))
# plt.plot(t_1[1::],np.diff(trend_1))
# plt.plot(t_2[1::],np.diff(trend_2))
# plt.plot(t_3[1::],np.diff(imfs_3[6]+imfs_3[5]))


# %% -----------------------------------------------------------------------------------------------
# KPSS test to check for stationarity
# If the result = 'Not stationary', a deterministc trend / linear regression is not suitable

print('Checking for stationarity')
KPSS_result = []
stationarity_array = []
for n in range(len(depths)):
    KPSS_result.append(TF.kpss_test((Tbin_m[n]))) 
    a = KPSS_result[n]
    stationarity_array.append(str(depths[n]) + ' m :  ' + a.KPSS_result)       
      
del a, n
    
print(stationarity_array)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Save data as mat file for plotting etc.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html

# convert time to string
tbin_m_str = []
tbin_deseason_str = []
for nn in range(len(tbin_m)):
    ttt = tbin_m[nn]
    a = []
    for n in range(len(ttt)):
        tt = ttt[n]
        a.append(str(tt))
    tbin_m_str.append(a)
    b = []    
    yr, mn, dy, hr, yday = TF.datevec(ttt)
    for n in range(len(yr)):
        d = dt.datetime(yr[n],mn[n],dy[n],hr[n])
        b.append(d.strftime("%Y-%m-%d %H:%M:%S"))
    tbin_deseason_str.append(b)
   
EEMD_t_str = []
for nn in range(len(EEMD_t)):
    ttt = EEMD_t[nn]
    a = []  
    for n in range(len(ttt)):
        tt = ttt[n]
        a.append(tt.strftime("%Y-%m-%d %H:%M:%S"))
    EEMD_t_str.append(a)


Trend_dict = {'ACF': ACF_result,
'KPSS_results': KPSS_result,
'EEMD_t': EEMD_t_str,
'EEMD_T': EEMD_T,
'EEMD_trend': EEMD_trend,
'EEMD_imfs': EEMD_IMFS,
'EEMD_res': EEMD_res,
'EEMD_conf_std_limit': conf_std_limit,
'EEMD_std_array': std_array,
'EEMD_trend_sims': trend_sims,
'EEMD_sims': x_sims}

Data_dict = {'tbin': tbin_m_str,
'Tbin': Tbin,
't': tbin_deseason_str,
'T': T,
'D': D,
'Tbin_deseason': Tbin_deseason,
'clims': clim,
'NRSMAI_agg': NRSMAI_agg}

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "NRSMAI_trends.mat", Trend_dict)

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "NRSMAI_data.mat", Data_dict)



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>









