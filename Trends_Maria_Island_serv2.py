

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

system = 0; # for windows (1), linux (0)

if system == 1:
    # mooring data
    main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
    NRSMAI_clim = xr.open_dataset(main_path + 'Data\\MariaIsland_TEMP_Climatology_1944-2020_BottleCTDMooringSatellite.nc')
    NRSMAI_agg = xr.open_dataset(main_path + 'Data\\MariaIsland_TEMP_1944-2020_aggregated.nc')
    del main_path
else:
  # mooring data
    NRSMAI_clim = xr.open_dataset('MariaIsland_TEMP_Climatology_1944-2020_BottleCTDMooringSatellite.nc')
    NRSMAI_agg = xr.open_dataset('MariaIsland_TEMP_1944-2020_aggregated.nc')


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

print('Selecting data at different depths:')

depths = [2, 20, 50]

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

# print('Removing the season')


# # get de-seasoned temperatures
# Tbin_deseason = []
# for n in range(len(depths)):
#     cl = clim[:,n]
#     Tbin_deseason.append(np.array(TF.deseason(tbin[n],Tbin[n],cl)))
# del n, cl
    
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
Tbin_m_deseason = []
Tbin_m_deseason_nofill = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # This is done to get a regular time grid with daily resolution
    tt,TT = TF.bin_monthly(1940,2021,t[n],np.float64(T[n]))
    TT,TTnoDS,_,non_filled = TF.fill_gaps(tt,TT,np.squeeze(clim[:,n]),30*12)
    Tbin_m_deseason_nofill.append(non_filled)      
    tbin_m.append(tt)
    Tbin_m.append(TTnoDS)
    Tbin_m_deseason.append(TT)
    # _,TT = TF.bin_daily(2011,2021,t[n],np.float64(T[n]))
    # Tbin_no_deseason.append(TT)
    
del tt, TT, n  
        

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
    tt = tbin_m[n]; TT = Tbin_m_deseason[n];
    TTT = []
    for nn in range(len(TT)):
        TTT.append(TT[nn])
    TTT = np.array(TTT)
    mk_result.append(
        mk.trend_free_pre_whitening_modification_test(TTT))
    mk_pval.append(mk_result[n].p)
    mk_trend.append(range(len(tt[np.isfinite(TTT)]))
                    *mk_result[n].slope + mk_result[n].intercept)
    a = mk_trend[n]
    tr = (a[-1]-a[0]) * (120/len(tt[np.isfinite(TTT)]));
    mk_trend_per_decade.append(tr)
    # Seasons
    yr, mn, dy, hr, yday = TF.datevec(tt)
    c_summer = np.squeeze(np.logical_or([mn == 12],[mn <= 2]))
    c_autumn = np.squeeze(np.logical_and([mn > 2],[mn <= 5]))
    c_winter = np.squeeze(np.logical_and([mn > 5],[mn <= 8]))
    c_spring = np.squeeze(np.logical_and([mn > 8],[mn <= 11]))
    # Summer
    mk_res = mk.trend_free_pre_whitening_modification_test(TTT[c_summer])
    a = range(len(tt[np.squeeze(
        np.logical_and([np.isfinite(TTT)],[c_summer]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (120/len(tt[np.isfinite(TTT)]));
    mk_trend_per_decade_Su.append(tr)    
    mk_pval_Su.append(mk_res.p)
    # Autumn
    mk_res = mk.trend_free_pre_whitening_modification_test(TTT[c_autumn])
    a = range(len(tt[np.squeeze(
        np.logical_and([np.isfinite(TTT)],[c_autumn]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (120/len(tt[np.isfinite(TTT)]));
    mk_trend_per_decade_Au.append(tr)  
    mk_pval_Au.append(mk_res.p)
    # Winter
    mk_res = mk.trend_free_pre_whitening_modification_test(TTT[c_winter])
    a = range(len(tt[np.squeeze(
    np.logical_and([np.isfinite(TTT)],[c_winter]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (120/len(tt[np.isfinite(TTT)]));
    mk_trend_per_decade_Wi.append(tr)   
    mk_pval_Wi.append(mk_res.p)
    # Spring
    mk_res = mk.trend_free_pre_whitening_modification_test(TTT[c_spring])
    a = range(len(tt[np.squeeze(
        np.logical_and([np.isfinite(TTT)],[c_spring]))]))*mk_res.slope + mk_res.intercept
    tr = (a[-1]-a[0]) * (120/len(tt[np.isfinite(TTT)]));
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
# ITA_slope_per_decade_high = []
# ITA_slope_per_decade_low = []
for n in range(len(depths)):
    TT = Tbin_m_deseason_nofill[n];
    TTT = []
    for nn in range(len(TT)):
        TTT.append(TT[nn])
    TTT = np.array(TTT)    
    tt = TF.to_date64(tbin_m[n])
    ITA_stats.append(TF.ITA(tt,TTT,-1,0))
    a = ITA_stats[n]
    ITA_significance.append(a.ITA_significance)
    ITA_slope_per_decade.append(a.ITA_trend_sen_per_decade)
    # ITA_slope_per_decade_high.append(a.ITA_trend_sen_high_per_decade)
    # ITA_slope_per_decade_low.append(a.ITA_trend_sen_low_per_decade)
    
del n, a


r = np.arange(0,len(ITA_slope_per_decade),1)

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
    TT = Tbin_m_deseason_nofill[n];
    TTT = []
    for nn in range(len(TT)):
        TTT.append(TT[nn])
    TTT = np.array(TTT)        
    KPSS_result.append(TF.kpss_test((TTT))) 
    a = KPSS_result[n]
    stationarity_array.append(str(depths[n]) + ' m :  ' + a.KPSS_result)       
      
del a, n
    
print(stationarity_array)

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
# for n in range(len(depths)):
#     print(str(depths[n]) + ' m')
#     tt = tbin_m[n]; TT = Tbin_m[n];
#     TTT = []
#     for nn in range(len(TT)):
#         TTT.append(TT[nn])
#     TTT = np.array(TTT)      
#     t, T, trend, trend_EAC, imfs, imfs_std, imfs_to_ave, res = TF.Ensemble_EMD(tt,TTT,0,1)
#     EEMD_t.append(t)
#     EEMD_T.append(T)
#     EEMD_trend.append(trend)
#     EEMD_trend_EAC.append(trend_EAC)
#     EEMD_imfs.append(imfs)
#     EEMD_res.append(res)
#     # Seasons
#     yr, mn, dy, hr, yday = TF.datevec(tt)
#     # Summer
#     c_summer = np.squeeze(np.logical_or([mn == 12],[mn <= 2]))
#     _, _, trend, trend_EAC, _, _, _, _ = TF.Ensemble_EMD(
#         tt[c_summer],TT[c_summer],0,1)
#     EEMD_trend_Su.append(trend); EEMD_trend_EAC_Su.append(trend_EAC); 
#     # Autumn
#     c_autumn = np.squeeze(np.logical_and([mn > 2],[mn <= 5]))
#     _, _, trend, trend_EAC, _, _, _, _ = TF.Ensemble_EMD(
#         tt[c_autumn],TT[c_autumn],0,1)
#     EEMD_trend_Au.append(trend); EEMD_trend_EAC_Au.append(trend_EAC);     
#     # Winter
#     c_winter = np.squeeze(np.logical_and([mn > 5],[mn <= 8]))
#     _, _, trend, trend_EAC, _, _, _, _ = TF.Ensemble_EMD(
#         tt[c_winter],TT[c_winter],0,1)
#     EEMD_trend_Wi.append(trend); EEMD_trend_EAC_Wi.append(trend_EAC);     
#     # Spring
#     c_spring = np.squeeze(np.logical_and([mn > 8],[mn <= 11]))
#     _, _, trend, trend_EAC, _, _, _, _ = TF.Ensemble_EMD(
#         tt[c_spring],TT[c_spring],0,1)
#     EEMD_trend_Sp.append(trend); EEMD_trend_EAC_Sp.append(trend_EAC); 
    

# EEMD_IMFS = {'IMF_1':EEMD_imfs[0],
#              'IMF_2':EEMD_imfs[1],
#              'IMF_3':EEMD_imfs[2]}
    
# plt.plot(EEMD_t[0],EEMD_trend[0])
# plt.plot(EEMD_t[1],EEMD_trend[1])
# plt.plot(EEMD_t[2],EEMD_trend[2])
# plt.plot(EEMD_t[3],EEMD_trend[3])
# plt.plot(EEMD_t[4],EEMD_trend[4])
# plt.plot(EEMD_t[5],EEMD_trend[5])
# plt.plot(EEMD_t[6],EEMD_trend[6])

# plt.plot(EEMD_t[0],EEMD_T[0],'.')
# plt.plot(EEMD_t[0],EEMD_trend[0])
# plt.plot(EEMD_t[0],EEMD_trend_EAC[0])

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
    TT = Tbin_m[n]
    TTT = []
    for nn in range(len(TT)):
        TTT.append(TT[nn])
    TT = np.array(TTT)  
    tt = tbin_m[n]
    TT = TT[np.isfinite(TT)]
    ACF_result.append(np.array(pd.Series(sm.tsa.acf(TT, nlags=10))))
    # significance (using monthly values)
    tt,TT = TF.bin_monthly(1940,2021,tbin_m[n],Tbin_m[n])
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
            TF.EEMD_significance(tt[c_summer],TT[c_summer],ACF_result[n],1000)    
    conf_std_limit_Su.append(csl)
    conf_std_limit_EAC_Su.append(csl_EAC)
#     # autumn
#     c_autumn = np.squeeze(np.logical_and([mn > 2],[mn <= 5]))
#     csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
#            TF.EEMD_significance(tt[c_autumn],TT[c_autumn],ACF_result[n],1000)    
#     conf_std_limit_Au.append(csl)
#     conf_std_limit_EAC_Au.append(csl_EAC)
#     # winter
#     c_winter = np.squeeze(np.logical_and([mn > 5],[mn <= 8]))
#     csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
#            TF.EEMD_significance(tt[c_winter],TT[c_winter],ACF_result[n],1000)    
#     conf_std_limit_Wi.append(csl)
#     conf_std_limit_EAC_Wi.append(csl_EAC)
#     # spring
#     c_spring = np.squeeze(np.logical_and([mn > 8],[mn <= 11]))
#     csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs = \
#            TF.EEMD_significance(tt[c_spring],TT[c_spring],ACF_result[n],1000)    
#     conf_std_limit_Sp.append(csl)
#     conf_std_limit_EAC_Sp.append(csl_EAC)


# del TT, n, csl, csl_EAC, sa, sa_EAC, ts, ts_EAC, xs

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

Data_dict = {'tbin': tbin_m_str,
'Tbin': Tbin_m,
't': tbin_deseason_str,
'T': T,
'D': D,
'Tbin_m_deseason': Tbin_m_deseason,
'clims': clim,
'NRSMAI_agg': NRSMAI_agg}


system = 0; # for windows (1), linux (0)

if system == 1:
    savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
            "NRSMAI_trends.mat", Trend_dict)
    savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
            "NRSMAI_data.mat", Data_dict)
else:
    savemat("NRSMAI_trends_server2.mat", Trend_dict)
    savemat("NRSMAI_data_server2.mat", Data_dict)    



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>









