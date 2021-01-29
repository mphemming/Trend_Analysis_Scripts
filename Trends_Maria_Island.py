

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

del main_path

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

print('Selecting data at different depths:')

depths = [2, 10, 20, 30, 40, 50, 85]

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


# clim = np.ones((365,7),dtype=float)
# for day in range(0,365):
#     day_temps = NRSMAI_clim.TEMP_AVE[day,:] 
#     clim[day,:] = np.interp(depths,NRSMAI_clim.DEPTH,day_temps)
  
# calculate simple climatology for now
clim = np.ones((12,7),dtype=float) 
for n in range(len(depths)):
    c = TF.calc_clim_monthly(t[n],T[n])  
    clim[:,n] = c

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

# %% -----------------------------------------------------------------------------------------------
# Get monthly averages and gap-fill

# Using de-seasoned timeseries
tbin_m = []
Tbin_m = []
Tbin_m_NG = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    tt,TT = TF.bin_monthly(1945,2021,tbin[n],Tbin[n])
    TT,TTnoDS,_ = TF.fill_gaps(tt,TT,np.squeeze(clim[:,n]),30*12)
    tbin_m.append(tt)
    Tbin_m.append(TT)    
    tt,TT = TF.bin_monthly(1945,2021,tbin[n],Tbin_deseason[n])
    Tbin_m_NG.append(TT)  
    if n == 6:
        check = tbin_m[n] < dt.datetime(2008,8,1)
        TT = Tbin_m[n]; TT[check] = np.nan
        Tbin_m[n] = TT
        TT = Tbin_m_NG[n]; TT[check] = np.nan
        Tbin_m_NG[n] = TT.astype('float64')  
    if n == 3:
        check = tbin_m[n] < dt.datetime(1956,10,1)
        TT = Tbin_m[n]; TT[check] = np.nan
        Tbin_m[n] = TT.astype('float64') 
        

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')

EEMD_t = []
EEMD_T = []
EEMD_trend = []
EEMD_trend_EAC = []
EEMD_imfs = []
EEMD_res = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    t, T, trend, trend_EAC, imfs, res = TF.Ensemble_EMD(tbin_m[n],Tbin_m[n],0)
    EEMD_t.append(t)
    EEMD_T.append(T)
    EEMD_trend.append(trend)
    EEMD_trend_EAC.append(trend_EAC)
    EEMD_imfs.append(imfs)
    EEMD_res.append(res)
    

EEMD_IMFS = {'IMF_1':EEMD_imfs[0],
             'IMF_2':EEMD_imfs[1],
             'IMF_3':EEMD_imfs[2], 
             'IMF_4':EEMD_imfs[3], 
             'IMF_5':EEMD_imfs[4], 
             'IMF_6':EEMD_imfs[5],
             'IMF_7':EEMD_imfs[6]}
    
# plt.plot(EEMD_t[0],EEMD_trend[0])
# plt.plot(EEMD_t[1],EEMD_trend[1])
# plt.plot(EEMD_t[2],EEMD_trend[2])
# plt.plot(EEMD_t[3],EEMD_trend[3])
# plt.plot(EEMD_t[4],EEMD_trend[4])
# plt.plot(EEMD_t[5],EEMD_trend[5])
# plt.plot(EEMD_t[6],EEMD_trend[6])

plt.plot(EEMD_t[0],EEMD_T[0],'.')
plt.plot(EEMD_t[0],EEMD_trend[0])
plt.plot(EEMD_t[0],EEMD_trend_EAC[0])

# Autocorrelation analysis and significance
print('Running autocorrelation analysis')
# Using last 10 years only

ACF_result = []
conf_std_limit = []
conf_std_limit_EAC = []
std_array = []
std_array_EAC = []
trend_sims = []
trend_sims_EAC = []
x_sims = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m') 
    check = np.where(np.logical_and([tbin_m[n] > dt.datetime(2010,1,1)], 
                   [tbin_m[n] < dt.datetime(2020,1,1)]))      
    TT = Tbin_m[n]
    tt = tbin_m[n]
    TT = TT[check[1]]
    TT = TT.astype('float64')
    TT = TT[np.isfinite(TT)]
    ACF_result.append(np.array(pd.Series(sm.tsa.acf(TT, nlags=10))))
    # significance
    TT = Tbin_m[n]
    TT = TT.astype('float64')
    csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
           TF.EEMD_significance(tbin_m[n],TT,ACF_result[n],1)
    conf_std_limit.append(csl)
    std_array.append(sa)
    trend_sims.append(ts)
    conf_std_limit_EAC.append(csl_EAC)
    std_array_EAC.append(sa_EAC)
    trend_sims_EAC.append(ts_EAC)    
    x_sims.append(xs)

del TT, n, check, csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs

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
# Mann kendall tests
print('Estimating Sen slopes and performing Mann Kendall tests')
mk_result = []
mk_trend = []
mk_trend_per_decade = []
mk_pval = []
for n in range(len(depths)):
    TT = Tbin_m[n]; TT = TT.astype('float64')
    mk_result.append(mk.trend_free_pre_whitening_modification_test(TT))
    mk_pval.append(mk_result[n].p)
    mk_trend.append(range(len(TT))*mk_result[n].slope + mk_result[n].intercept)
    tr = range(0,120)*mk_result[n].slope + mk_result[n].intercept
    mk_trend_per_decade.append(tr[-1]-tr[0])
    
del n, tr
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



# %% -----------------------------------------------------------------------------------------------
# Innovative trend analysis

ITA_stats = []
ITA_significance = []
ITA_slope_per_decade = []
# ITA_slope_per_decade_high = []
# ITA_slope_per_decade_low = []
for n in range(len(depths)):
    tt = TF.to_date64(tbin_m[n])
    TT = Tbin_m[n]; TT = TT.astype('float64')
    ITA_stats.append(TF.ITA(tt,TT,0,0))
    a = ITA_stats[n]
    ITA_significance.append(a.ITA_significance)
    ITA_slope_per_decade.append(a.ITA_trend_sen_per_decade)
    # ITA_slope_per_decade_high.append(a.ITA_trend_sen_high_per_decade)
    # ITA_slope_per_decade_low.append(a.ITA_trend_sen_low_per_decade)
    
# plt.plot(ITA_slope_per_decade_low,depths,color='b')
# plt.plot(ITA_slope_per_decade_high,depths,color='r')
plt.plot(ITA_slope_per_decade,depths,color='k')
plt.plot(mk_trend_per_decade,depths,color='g')

del n, a


r = np.arange(0,6,1)

for n in r:
    line = np.arange(start=-20, stop=20, step=1) 
    plt.plot(line,line,color='k')
    plt.scatter(ITA_stats[n].TEMP_half_1,ITA_stats[n].TEMP_half_2,2)
    plt.xlim(left=-4, right=4)
    plt.ylim(bottom=-4, top=4)

# %% -----------------------------------------------------------------------------------------------
# KPSS test to check for stationarity
# If the result = 'Not stationary', a deterministc trend / linear regression is not suitable

print('Checking for stationarity')
KPSS_result = []
stationarity_array = []
for n in range(len(depths)):
    TT = Tbin_m[n]; TT = TT.astype('float64')
    KPSS_result.append(TF.kpss_test((TT))) 
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


Trend_dict = {'MK_result': mk_result,
'MK_trend': mk_trend,
'MK_trend_per_decade': mk_trend_per_decade,
'MK_pval': mk_pval,
'KPSS_results': KPSS_result,
'ITA_stats': ITA_stats,
'ITA_significance': ITA_significance,
'ITA_trend_per_decade': ITA_slope_per_decade,
'ACF': ACF_result,
'KPSS_results': KPSS_result,
'EEMD_t': EEMD_t_str,
'EEMD_T': EEMD_T,
'EEMD_trend': EEMD_trend,
'EEMD_trend_EAC': EEMD_trend_EAC,
'EEMD_imfs': EEMD_IMFS,
'EEMD_res': EEMD_res,
'EEMD_conf_std_limit': conf_std_limit,
'EEMD_conf_std_limit_EAC': conf_std_limit_EAC,
'EEMD_std_array': std_array,
'EEMD_std_array_EAC': std_array_EAC,
'EEMD_trend_sims': trend_sims,
'EEMD_trend_sims_EAC': trend_sims_EAC,
'EEMD_sims': x_sims}

Data_dict = {'tbin': tbin_m_str,
'Tbin': Tbin_m,
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









