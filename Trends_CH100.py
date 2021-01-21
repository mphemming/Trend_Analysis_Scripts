
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
import scipy.stats as sts
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
CH100_agg = xr.open_dataset(main_path + 
    'Data\\IMOS_ANMN-NSW_TZ_20090815_CH100_FV01_TEMP-aggregated-timeseries_END-20201028_C-20201207.nc')

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

# code to check data distribution
# check = np.isfinite(CH100_agg.TEMP) 
# %matplotlib qt
# plt.hist(CH100_agg.DEPTH[check], bins = np.arange(0,120,1))
# plt.xlim(left=0, right=110)


print('Selecting data at different depths:')

depths = [10.5, 20, 27.5, 35.5, 43.5, 51.5, 59.5, 67.5, 75.5, 84.5, 91.5]

D = []
t = []
T = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # index check
    c = [(CH100_agg.DEPTH >= depths[n] - 2) & (CH100_agg.DEPTH <= depths[n] + 2) & \
        (CH100_agg.TEMP_quality_control > 0) & (CH100_agg.TEMP_quality_control <3)]
    # Depth
    d = np.array(CH100_agg.DEPTH);
    D.append(d[c])
    # time
    tt = np.array(CH100_agg.TIME);
    t.append(tt[c])    
    # Temp
    TT = np.array(CH100_agg.TEMP);
    T.append(TT[c])       

# plt.plot(t[0],T[0])
# plt.plot(t[10],T[10])
# plt.show()

del c, d, tt, TT, n

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# QC data
# n=0
TQC = T[0]; tQC = t[0];
c1 = [(TQC < 20) & (tQC > np.datetime64('2010-10-29')) & (tQC < np.datetime64('2011-04-09'))]
c2 = [(TQC < 17)]
c3 = [(TQC < 19.4) & (tQC > np.datetime64('2016-08-15')) & (tQC < np.datetime64('2016-08-19'))]
c4 = [(TQC < 18.06) & (tQC > np.datetime64('2016-11-05')) & (tQC < np.datetime64('2016-11-09'))]
TQC[c1] = np.nan; TQC[c2] = np.nan; TQC[c3] = np.nan; TQC[c4] = np.nan;  
T[0] = TQC; 
# n=1
TQC = T[1]; tQC = t[1];
c1 = [(TQC < 19.6) & (tQC > np.datetime64('2010-10-29')) & (tQC < np.datetime64('2010-11-21'))]
c2 = [(TQC < 17) & (tQC > np.datetime64('2010-12-01')) & (tQC < np.datetime64('2011-01-01'))]
c3 = [(TQC < 20) & (tQC > np.datetime64('2011-02-01')) & (tQC < np.datetime64('2011-05-01'))]
c4 = [(TQC < 22) & (tQC > np.datetime64('2011-04-09')) & (tQC < np.datetime64('2011-04-13'))]
TQC[c1] = np.nan; TQC[c2] = np.nan; TQC[c3] = np.nan; TQC[c4] = np.nan;  
T[1] = TQC; 
# n=2
TQC = T[2]; tQC = t[2];
c1 = [(TQC < 16.5) & (tQC > np.datetime64('2010-12-01')) & (tQC < np.datetime64('2010-12-25'))]
TQC[c1] = np.nan; 
T[2] = TQC; 
# n = 5
TQC = T[5]; tQC = t[5];
c1 = [(TQC < 14) & (tQC > np.datetime64('2010-12-01')) & (tQC < np.datetime64('2011-01-01'))]
TQC[c1] = np.nan; 
T[5] = TQC; 

del c1, c2, c3, c4

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

# print('Removing the season')

# # select climatology at similar depth
clim = []
for n in range(len(depths)):
    clim.append(TF.calc_clim_monthly(t[n],T[n]))
# # get de-seasoned temperatures
Tbin_deseason = []
for n in range(len(depths)):
    Tbin_deseason.append(np.array(TF.deseason(t[n],T[n],clim[n])))
    
del n

# interpolate climatologies to 365 days
t_months = [dt.datetime(1,1,1),
            dt.datetime(1,2,1),
            dt.datetime(1,3,1),
            dt.datetime(1,4,1),
            dt.datetime(1,5,1),
            dt.datetime(1,6,1),
            dt.datetime(1,7,1),
            dt.datetime(1,8,1),
            dt.datetime(1,9,1),
            dt.datetime(1,10,1),
            dt.datetime(1,11,1),
            dt.datetime(1,12,1),
            dt.datetime(1,12,31)]
_, _, _, _, yday = datevec(t_months)
clim_daily = []
for n in range(len(depths)):
    c = np.concatenate([clim[n],clim[n]])
    c = np.stack(c).astype(None)
    a = c[0:13]
    clim_daily.append(np.interp(np.arange(0,367,1),yday,a))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get daily averages
  
print('Getting Daily Averages')

# Using de-seasoned timeseries
tbin = []
Tbin = []
Tbin_no_deseason = []
choice = 1
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # This is done to get a regular time grid with daily resolution
    if choice == 1:
        tt,TT = TF.bin_daily(2009,2021,t[n],np.float64(T[n]))
        TT,TTnoDS,_ = TF.fill_gaps(tt,TT,np.squeeze(clim_daily[n]),2*365)
    else:
        tt,TT = TF.bin_monthly(2009,2021,tbin[n],Tbin_deseason[n])           
    tbin.append(tt)
    Tbin.append(TT)
    _,TT = TF.bin_daily(2009,2021,t[n],np.float64(T[n]))
    Tbin_no_deseason.append(TT)
    
del tt, TT, n  

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# KPSS test to check for stationarity
# If the result = 'Not stationary', a deterministc trend / linear regression is not suitable

print('Checking for stationarity')
KPSS_result = []
stationarity_array = []
for n in range(len(depths)):
    KPSS_result.append(TF.kpss_test((Tbin[n]))) 
    a = KPSS_result[n]
    stationarity_array.append(str(depths[n]) + ' m :  ' + a.KPSS_result)       
      
del a, n
    
print(stationarity_array)

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

del n, a

# # Autocorrelation analysis
# n = 0
# TT = Tbin[n]
# check = np.where(np.isnan(TT))
# ACF1 = pd.Series(sm.tsa.acf(TT[1550:1690], nlags=100)); # where longest streak without nans
# ACF2 = pd.Series(sm.tsa.acf(TT[2680:2770], nlags=100)); # where longest streak without 
# ACF3 = pd.Series(sm.tsa.acf(TT[3785:3920], nlags=100)); # where longest streak without 
# ACF_result = []
# for n in range(0,50):
#     ACF_result.append(np.nanmean([ACF1[n],ACF2[n],ACF3[n]]))
# ACF_result = np.array(ACF_result)

# # 95% confidence level for no trend case
# # confidence_ITA = TF.ITA_significance(tbin[0],Tbin[0],ACF_result,500)

# # https://www.mathsisfun.com/data/confidence-interval.html


# line = np.arange(start=-20, stop=20, step=1) 
# conf = 1.96*(np.nanstd(line)/np.sqrt(len(line)))
# # line = np.arange(start=2500, stop=5700, step=1) 
# line_upper = line+conf
# line_lower = line-conf

r = [0, 5, 10]

for n in r:
    line = np.arange(start=-20, stop=20, step=1) 
    plt.plot(line,line,color='k')
    plt.scatter(ITA_stats[n].TEMP_half_1,ITA_stats[n].TEMP_half_2,2)
    plt.xlim(left=-4, right=4)
    plt.ylim(bottom=-4, top=4)


# plt.plot(ITA_stats.TEMP_half_1,ITA_stats.TEMP_half_1*ITA_stats.ITA_trend_sen_period)



# plt.plot(ITA_stats.TEMP_half_1,
#    ITA_stats.TEMP_half_2*(ITA_stats.ITA_trend_sen_period+ITA_stats.ITA_slope_upper_conf),color='r')
# plt.plot(ITA_stats.TEMP_half_1,
#    ITA_stats.TEMP_half_2*(ITA_stats.ITA_trend_sen_period-ITA_stats.ITA_slope_upper_conf),color='r')
# plt.scatter(ITA_stats.TEMP_half_1,ITA_stats.TEMP_half_2)
# plt.xlim(left=-2, right=2)
# plt.ylim(bottom=-2, top=2)

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
    t, T, trend, imfs, res = TF.Ensemble_EMD(tbin[n],Tbin[n],0)
    EEMD_t.append(t)
    EEMD_T.append(T)
    EEMD_trend.append(trend)
    EEMD_imfs.append(imfs)
    EEMD_res.append(res)
    

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
    check = np.where(np.logical_and([tbin[n] > dt.datetime(2010,1,1)], 
                   [tbin[n] < dt.datetime(2020,1,1)]))      
    TT = Tbin[n]
    tt = tbin[n]
    TT = TT[check[1]]
    TT = TT[np.isfinite(TT)]

    ACF_result.append(np.array(pd.Series(sm.tsa.acf(TT, nlags=10))))

    # significance
    csl, sa, ts, xs = \
           TF.EEMD_significance(tbin[n],Tbin[n],ACF_result[n],1000)
    conf_std_limit.append(csl)
    std_array.append(sa)
    trend_sims.append(ts)
    x_sims.append(xs)

del TT, n, check, csl, sa, ts, xs


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
# convert time to string
tbin_str = []
tbin_deseason_str = []
for nn in range(len(tbin)):
    ttt = tbin[nn]
    a = []
    for n in range(len(ttt)):
        tt = ttt[n]
        a.append(str(tt))
    tbin_str.append(a)
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
    yr, mn, dy, hr, yday = TF.datevec(ttt)
    for n in range(len(ttt)):
        tt = ttt[n]
        d = dt.datetime(yr[n],mn[n],dy[n],hr[n])
        a.append(d.strftime("%Y-%m-%d %H:%M:%S"))        
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
'EEMD_imfs': EEMD_imfs,
'EEMD_res': EEMD_res,
'EEMD_conf_std_limit': conf_std_limit,
'EEMD_std_array': std_array,
'EEMD_trend_sims': trend_sims,
'EEMD_sims': x_sims}


Data_dict = {'tbin': tbin_str,
'Tbin': Tbin,
't': tbin_deseason_str,
'T': T,
'D': D,
'Tbin_deseason': Tbin_deseason,
'clims': clim,
'CH100_agg': CH100_agg}

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "CH100_trends.mat", Trend_dict)

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "CH100_data.mat", Data_dict)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


