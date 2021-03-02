
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
NRSNSI_agg = xr.open_dataset(main_path + 
    'Data\\IMOS_ANMN-NRS_TZ_20101213_NRSNSI_FV01_TEMP-aggregated-timeseries_END-20201028_C-20201207.nc')
# Satellite data
NRSNSI_sat = xr.open_dataset(main_path + 
    'Data\\Satellite_NRSNSI.nc')                       

# %% -----------------------------------------------------------------------------------------------
# Select satellite data and correct

SST = np.array(NRSNSI_sat.sea_surface_temperature)
sses_bias = np.array(NRSNSI_sat.sses_bias)
lon = np.array(NRSNSI_sat.lon)
lat = np.array(NRSNSI_sat.lat)
time = np.array(NRSNSI_sat.time)
dtime = np.array(NRSNSI_sat.sst_dtime)
QC = np.array(NRSNSI_sat.quality_level)

A_SST = []
A_SSES = []
A_dtime = []
for day in range(np.size(SST,0)):
    # select good QC data
    QC_day = QC[day,24:30,5:11]
    c = QC_day >= 4
    # apply QC and average data per day
    sst = SST[day,24:30,5:11]
    A_SST.append(np.nanmean(sst[c]))
    sses = sses_bias[day,24:30,5:11]
    A_SSES.append(np.nanmean(sses[c]))
    lon_sst = lon[5:11]
    lat_sst = lat[24:30]
    dtt = dtime[day,24:30,5:11]
    A_dtime.append(np.nanmedian(dtt[c]))
# Correct SST
SST_corr = np.array(A_SST) - np.array(A_SSES)
SST_corr = SST_corr + 0.17
SST_corr = SST_corr - 273.15
# Sort out time
# dt is between -0.3 and 0.2 days so won;t affect trend analysis!
# check dtime, if it won't change the time of day too much no point 
# in correcting time for dt. Daily data used for trends anyway

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

# code to check data distribution
# check = np.isfinite(NRSNSI_agg.TEMP) 
# %matplotlib qt
# plt.hist(NRSNSI_agg.DEPTH[check], bins = np.arange(0,120,1))
# plt.xlim(left=50, right=70)

print('Selecting data at different depths:')

depths = [0, 22, 61]

D = []
t = []
T = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    
    if depths[n] == 0:
        # bin satellite data 
        t.append(time)  
        T.append(SST_corr)
    else:
        # index check
        c = [(NRSNSI_agg.DEPTH >= depths[n] - 2) & (NRSNSI_agg.DEPTH <= depths[n] + 2) & \
            (NRSNSI_agg.TEMP_quality_control > 0) & (NRSNSI_agg.TEMP_quality_control <3)]
        # Depth
        d = np.array(NRSNSI_agg.DEPTH);
        D.append(d[c])
        # time
        tt = np.array(NRSNSI_agg.TIME);
        t.append(tt[c])    
        # Temp
        TT = np.array(NRSNSI_agg.TEMP);
        T.append(TT[c]) 
             

# plt.plot(t[0],T[0])
# plt.plot(t[10],T[10])
# plt.show()

del c, d, tt, TT

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# QC data

# 0m
TQC = T[0]; tQC = t[0];
c1 = [(TQC < 18.9) & (tQC > np.datetime64('2019-09-01')) &
      (tQC < np.datetime64('2019-09-15'))]
c2 = [(TQC < 19.9) & (tQC > np.datetime64('2018-08-01')) &
      (tQC < np.datetime64('2018-09-01'))]
c3 = [(TQC < 20.1) & (tQC > np.datetime64('2016-07-01')) &
      (tQC < np.datetime64('2016-08-01'))]
c4 = [(TQC < 19.6) & (tQC > np.datetime64('2013-07-01')) &
      (tQC < np.datetime64('2013-08-01'))]
TQC[c1] = np.nan; TQC[c2] = np.nan;  TQC[c3] = np.nan;  TQC[c4] = np.nan;  
T[0] = TQC; 
# 22m
TQC = T[1]; tQC = t[1];
c1 = [(TQC > 28.6) & (tQC > np.datetime64('2018-03-15')) & 
      (tQC < np.datetime64('2018-04-01'))]
TQC[c1] = np.nan; 
T[1] = TQC; 

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
_, _, _, _, yday = TF.datevec(t_months)
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
Tbin_deseason = []
Tbin_deseason_nofill = []
choice = 1
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # This is done to get a regular time grid with daily resolution
    if choice == 1:
        tt,TT = TF.bin_daily(2011,2021,t[n],np.float64(T[n]))
        TT,TTnoDS,_,non_filled = TF.fill_gaps(tt,TT,np.squeeze(clim_daily[n]),2*365)
        Tbin_deseason_nofill.append(non_filled)
    else:
        tt,TT = TF.bin_monthly(2011,2021,tbin[n],Tbin_deseason[n])           
    tbin.append(tt)
    Tbin.append(TTnoDS)
    Tbin_deseason.append(TT)
    # _,TT = TF.bin_daily(2011,2021,t[n],np.float64(T[n]))
    # Tbin_no_deseason.append(TT)
    
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
mk_trend_per_decade_Su = []
mk_trend_per_decade_Au = []
mk_trend_per_decade_Wi = []
mk_trend_per_decade_Sp = []
mk_pval = []
mk_pval_Su = []
mk_pval_Au = []
mk_pval_Wi = []
mk_pval_Sp = []

for n in range(len(depths)):
    tt = tbin[n]; TT = Tbin_deseason[n];
    mk_result.append(
        mk.trend_free_pre_whitening_modification_test(TT))
    mk_pval.append(mk_result[n].p)
    mk_trend.append(range(len(tt[np.isfinite(TT)]))
                    *mk_result[n].slope + mk_result[n].intercept)
    a = mk_trend[n]
    tr = (a[-1]-a[0]) * (3652/len(tt[np.isfinite(TT)]));
    mk_trend_per_decade.append(tr)
    # Seasons
    yr, mn, dy, hr, yday = TF.datevec(tt)
    c_summer = np.squeeze(np.logical_or([mn == 12],[mn <= 2]))
    c_autumn = np.squeeze(np.logical_and([mn > 2],[mn <= 5]))
    c_winter = np.squeeze(np.logical_and([mn > 5],[mn <= 8]))
    c_spring = np.squeeze(np.logical_and([mn > 8],[mn <= 11]))
    # Summer
    mk_res = mk.trend_free_pre_whitening_modification_test(TT[c_summer])
    a = range(len(tt[np.squeeze(
        np.logical_and([np.isfinite(TT)],[c_summer]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (3652/len(tt[np.isfinite(TT)]));
    mk_trend_per_decade_Su.append(tr)    
    mk_pval_Su.append(mk_res.p)
    # Autumn
    mk_res = mk.trend_free_pre_whitening_modification_test(TT[c_autumn])
    a = range(len(tt[np.squeeze(
        np.logical_and([np.isfinite(TT)],[c_autumn]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (3652/len(tt[np.isfinite(TT)]));
    mk_trend_per_decade_Au.append(tr)  
    mk_pval_Au.append(mk_res.p)
    # Winter
    mk_res = mk.trend_free_pre_whitening_modification_test(TT[c_winter])
    a = range(len(tt[np.squeeze(
    np.logical_and([np.isfinite(TT)],[c_winter]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (3652/len(tt[np.isfinite(TT)]));
    mk_trend_per_decade_Wi.append(tr)   
    mk_pval_Wi.append(mk_res.p)
    # Spring
    mk_res = mk.trend_free_pre_whitening_modification_test(TT[c_spring])
    a = range(len(tt[np.squeeze(
        np.logical_and([np.isfinite(TT)],[c_spring]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (3652/len(tt[np.isfinite(TT)]));
    mk_trend_per_decade_Sp.append(tr)   
    mk_pval_Sp.append(mk_res.p)     
    
del n, tr
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Innovative trend analysis

ITA_stats = []
ITA_significance = []
ITA_slope_per_decade = []

for n in range(len(depths)):
    tt = TF.to_date64(tbin[n])
    ITA_stats.append(TF.ITA(tt,Tbin_deseason_nofill[n],-1,0))
    a = ITA_stats[n]
    ITA_significance.append(a.ITA_significance)
    ITA_slope_per_decade.append(a.ITA_trend_sen_per_decade)

# plt.plot(ITA_slope_per_decade,depths,color='k')
# plt.plot(mk_trend_per_decade,depths,color='g')

del n, a


r = np.arange(0,len(ITA_slope_per_decade),1)

for n in r:
    line = np.arange(start=-20, stop=20, step=1) 
    plt.plot(line,line,color='k')
    plt.scatter(ITA_stats[n].TEMP_half_1,ITA_stats[n].TEMP_half_2,2)
    plt.xlim(left=-4, right=4)
    plt.ylim(bottom=-4, top=4)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')

EEMD_t = []
EEMD_T = []
EEMD_trend = []
EEMD_trend_Su = []
EEMD_trend_Au = []
EEMD_trend_Wi = []
EEMD_trend_Sp = []
EEMD_trend_EAC = []
EEMD_trend_EAC_Su = []
EEMD_trend_EAC_Au = []
EEMD_trend_EAC_Wi = []
EEMD_trend_EAC_Sp = []
EEMD_imfs = []
EEMD_res = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    tt = tbin[n]; TT = Tbin[n];
    t, T, trend, trend_EAC, imfs, res = TF.Ensemble_EMD(tt,TT,0)
    EEMD_t.append(t)
    EEMD_T.append(T)
    EEMD_trend.append(trend)
    EEMD_trend_EAC.append(trend_EAC)
    EEMD_imfs.append(imfs)
    EEMD_res.append(res)
    # Seasons
    yr, mn, dy, hr, yday = TF.datevec(tt)
    # Summer
    c_summer = np.squeeze(np.logical_or([mn == 12],[mn <= 2]))
    _, _, trend, trend_EAC, _, _ = TF.Ensemble_EMD(
        tt[c_summer],TT[c_summer],0)
    EEMD_trend_Su.append(trend); EEMD_trend_EAC_Su.append(trend_EAC); 
    # Autumn
    c_autumn = np.squeeze(np.logical_and([mn > 2],[mn <= 5]))
    _, _, trend, trend_EAC, _, _ = TF.Ensemble_EMD(
        tt[c_autumn],TT[c_autumn],0)
    EEMD_trend_Au.append(trend); EEMD_trend_EAC_Au.append(trend_EAC);     
    # Winter
    c_winter = np.squeeze(np.logical_and([mn > 5],[mn <= 8]))
    _, _, trend, trend_EAC, _, _ = TF.Ensemble_EMD(
        tt[c_winter],TT[c_winter],0)
    EEMD_trend_Wi.append(trend); EEMD_trend_EAC_Wi.append(trend_EAC);     
    # Spring
    c_spring = np.squeeze(np.logical_and([mn > 8],[mn <= 11]))
    _, _, trend, trend_EAC, _, _ = TF.Ensemble_EMD(
        tt[c_spring],TT[c_spring],0)
    EEMD_trend_Sp.append(trend); EEMD_trend_EAC_Sp.append(trend_EAC);     
    
EEMD_IMFS = {'IMF_1':EEMD_imfs[0],
             'IMF_2':EEMD_imfs[1],
             'IMF_3':EEMD_imfs[2]}
    
# Autocorrelation analysis and significance
print('Running autocorrelation analysis')
# Using last 10 years only

ACF_result = []
conf_std_limit = []
conf_std_limit_Su = []
conf_std_limit_Au = []
conf_std_limit_Wi = []
conf_std_limit_Sp = []
conf_std_limit_EAC = []
conf_std_limit_EAC_Su = []
conf_std_limit_EAC_Au = []
conf_std_limit_EAC_Wi = []
conf_std_limit_EAC_Sp = []
std_array = []
std_array_EAC = []
trend_sims = []
trend_sims_EAC = []
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
    # significance (using monthly values)
    del tt, TT
    tt,TT = TF.bin_monthly(2011,2021,tbin[n],Tbin[n])
    csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
           TF.EEMD_significance(tt,TT,ACF_result[n],1)
    conf_std_limit.append(csl)
    std_array.append(sa)
    trend_sims.append(ts)
    conf_std_limit_EAC.append(csl_EAC)
    std_array_EAC.append(sa_EAC)
    trend_sims_EAC.append(ts_EAC)    
    x_sims.append(xs)
    # For the seasons, changing std only for each season
    # Using same ACF result as whole time series
    # Summer
    yr, mn, dy, hr, yday = TF.datevec(tt)
    c_summer = np.squeeze(np.logical_or([mn == 12],[mn <= 2]))
    csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
           TF.EEMD_significance(tt[c_summer],TT[c_summer],ACF_result[n],1)    
    conf_std_limit_Su.append(csl)
    conf_std_limit_EAC_Su.append(csl_EAC)
    # autumn
    c_autumn = np.squeeze(np.logical_and([mn > 2],[mn <= 5]))
    csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
           TF.EEMD_significance(tt[c_autumn],TT[c_autumn],ACF_result[n],1)    
    conf_std_limit_Au.append(csl)
    conf_std_limit_EAC_Au.append(csl_EAC)
    # winter
    c_winter = np.squeeze(np.logical_and([mn > 5],[mn <= 8]))
    csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
           TF.EEMD_significance(tt[c_winter],TT[c_winter],ACF_result[n],1)    
    conf_std_limit_Wi.append(csl)
    conf_std_limit_EAC_Wi.append(csl_EAC)
    # spring
    c_spring = np.squeeze(np.logical_and([mn > 8],[mn <= 11]))
    csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
           TF.EEMD_significance(tt[c_spring],TT[c_spring],ACF_result[n],1)    
    conf_std_limit_Sp.append(csl)
    conf_std_limit_EAC_Sp.append(csl_EAC)


del TT, n, check, csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs



# %% -----------------------------------------------------------------------------------------------
# Save data as mat file

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
'MK_trend_per_decade_Su': mk_trend_per_decade_Su,
'MK_trend_per_decade_Au': mk_trend_per_decade_Au,
'MK_trend_per_decade_Wi': mk_trend_per_decade_Wi,
'MK_trend_per_decade_Sp': mk_trend_per_decade_Sp,
'MK_pval': mk_pval,
'MK_pval_Su': mk_pval_Su,
'MK_pval_Au': mk_pval_Au,
'MK_pval_Wi': mk_pval_Wi,
'MK_pval_Sp': mk_pval_Sp,
'KPSS_results': KPSS_result,
'ITA_stats': ITA_stats,
'ITA_significance': ITA_significance,
'ITA_trend_per_decade': ITA_slope_per_decade,
'ACF': ACF_result,
'KPSS_results': KPSS_result,
'EEMD_t': EEMD_t_str,
'EEMD_T': EEMD_T,
'EEMD_trend': EEMD_trend,
'EEMD_trend_Su': EEMD_trend_Su,
'EEMD_trend_Au': EEMD_trend_Au,
'EEMD_trend_Wi': EEMD_trend_Wi,
'EEMD_trend_Sp': EEMD_trend_Sp,
'EEMD_trend_EAC': EEMD_trend_EAC,
'EEMD_trend_EAC_Su': EEMD_trend_EAC_Su,
'EEMD_trend_EAC_Au': EEMD_trend_EAC_Au,
'EEMD_trend_EAC_Wi': EEMD_trend_EAC_Wi,
'EEMD_trend_EAC_Sp': EEMD_trend_EAC_Sp,
'EEMD_imfs': EEMD_IMFS,
'EEMD_res': EEMD_res,
'EEMD_conf_std_limit': conf_std_limit,
'EEMD_conf_std_limit_Su': conf_std_limit_Su,
'EEMD_conf_std_limit_Au': conf_std_limit_Au,
'EEMD_conf_std_limit_Wi': conf_std_limit_Wi,
'EEMD_conf_std_limit_Sp': conf_std_limit_Sp,
'EEMD_conf_std_limit_EAC': conf_std_limit_EAC,
'EEMD_conf_std_limit_EAC_Su': conf_std_limit_EAC_Su,
'EEMD_conf_std_limit_EAC_Au': conf_std_limit_EAC_Au,
'EEMD_conf_std_limit_EAC_Wi': conf_std_limit_EAC_Wi,
'EEMD_conf_std_limit_EAC_Sp': conf_std_limit_EAC_Sp,
'EEMD_std_array': std_array,
'EEMD_std_array_EAC': std_array_EAC,
'EEMD_trend_sims': trend_sims,
'EEMD_trend_sims_EAC': trend_sims_EAC,
'EEMD_sims': x_sims}

Data_dict = {'tbin': tbin_str,
'Tbin': Tbin,
't': tbin_deseason_str,
'T': T,
'D': D,
'Tbin_deseason': Tbin_deseason,
'clims': clim,
'NRSNSI_agg': NRSNSI_agg}

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "NRSNSI_trends.mat", Trend_dict)

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "NRSNSI_data.mat", Data_dict)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


