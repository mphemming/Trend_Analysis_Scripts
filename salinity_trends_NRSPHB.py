
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

# mooring data
main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
NRSPHB_clim = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_Climatology_1953-2019_BottleCTDMooringSatellite.nc')
NRSPHB_agg = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_1953-2019_aggregated.nc')
Sal = pd.read_csv(main_path + 'Data\\NRSPHB_1953_2010.csv',index_col=3)

del main_path

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

print('Selecting data at different depths:')

depths = [2, 50, 75]

D = []
Ds = []
t = []
ts = []
T = []
S = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # index check
    c = [(NRSPHB_agg.DEPTH >= depths[n] - 2) & (NRSPHB_agg.DEPTH <= depths[n] + 2)]
    cs = [(Sal.PRESSURE >= depths[n] - 2) & (Sal.PRESSURE <= depths[n] + 2)]    
    # Depth
    d = np.array(NRSPHB_agg.DEPTH);
    ds = np.array(Sal.PRESSURE);
    D.append(d[c])
    Ds.append(ds[cs])
    # time
    tt = np.array(NRSPHB_agg.TIME);
    tts = np.array(Sal.index);
    b = []
    for n in range(len(tts)):
        a = tts[n]
        cc = a[6:10] + '-' + a[3:5] + '-' + a[0:2]
        if '/' in cc or ' ' in cc:
            cc = a[5:9] + '-' + a[2:4] + '-' + a[0:1]
            if len(a[0:1]) < 2:
                cc = a[5:9] + '-' + a[2:4] + '-0' + a[0:1] 
        b.append(np.datetime64(cc)) 
    b = np.array(b)
    t.append(tt[c]) 
    ts.append(b[cs]) 
    # Temp
    TT = np.array(NRSPHB_agg.TEMP);
    SS = np.array(Sal.SALINITY_VALUE);
    T.append(TT[c])       
    S.append(SS[cs])     
# plt.plot(t[0],T[0])
# plt.plot(t[10],T[10])
# plt.show()

del c, cc, cs, d, ds, tt, tts, TT, SS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Some PSAL manual QC

# 2m
c = S[0] < 35.1
a = S[0]; a[c] = np.nan;
S[0] = a;
# 50m
c = S[1] < 35.1
a = S[1]; a[c] = np.nan;
S[1] = a;
# 75m
c = S[2] < 35
a = S[2]; a[c] = np.nan;
S[2] = a;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get climatology at same depths

# calculate simple climatology for now
clim = np.ones((12,3),dtype=float) 
clim_S = np.ones((12,3),dtype=float) 
for n in range(len(depths)):
    c = TF.calc_clim_monthly(t[n],T[n])  
    cs = TF.calc_clim_monthly(ts[n],S[n]) 
    clim[:,n] = c
    clim_S[:,n] = cs
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Bin the data
# This is done to get a regular time grid with daily resolution

tbin = []
Tbin = []
Sbin = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    tt, TT = TF.bin_daily(1953,2020,t[n],T[n])
    tbin.append(tt)
    Tbin.append(TT)
    _, SS = TF.bin_daily(1953,2020,ts[n],S[n])
    Sbin.append(SS)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

print('Removing the season')


# get de-seasoned temperatures
Tbin_deseason = []
Sbin_deseason = []
for n in range(len(depths)):
    cl = clim[:,n]
    cl_s = clim_S[:,n]
    Tbin_deseason.append(np.array(TF.deseason(tbin[n],Tbin[n],cl)))
    Sbin_deseason.append(np.array(TF.deseason(tbin[n],Sbin[n],cl_s)))
del n, cl, cl_s
    
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
Sbin_m = []
Sbin_m_NG = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # temperature
    tt,TT = TF.bin_monthly(1953,2021,tbin[n],Tbin[n])
    TT,TTnoDS,_ = TF.fill_gaps(tt,TT,np.squeeze(clim[:,n]),30*12)
    tbin_m.append(tt)
    Tbin_m.append(TT)    
    tt,TT = TF.bin_monthly(1953,2021,tbin[n],Tbin_deseason[n])
    Tbin_m_NG.append(TT) 
    # salinity
    tt,SS = TF.bin_monthly(1953,2021,tbin[n],Sbin[n])
    SS,SSnoDS,_ = TF.fill_gaps(tt,SS,np.squeeze(clim_S[:,n]),30*12)
    Sbin_m.append(SS)    
    tt,SS = TF.bin_monthly(1953,2021,tbin[n],Sbin_deseason[n])
    Sbin_m_NG.append(SS)     
        
del n, tt, TT, SS, SSnoDS, TTnoDS

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')

EEMD_t = []
EEMD_tS = []
EEMD_T = []
EEMD_S = []
EEMD_trend = []
EEMD_trend_EAC = []
EEMD_trend_S = []
EEMD_trend_EAC_S = []
EEMD_imfs = []
EEMD_res = []
EEMD_imfs_S = []
EEMD_res_S = []

for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    t, T, trend, trend_EAC, imfs, res = TF.Ensemble_EMD(tbin_m[n],Tbin_m[n],0)
    EEMD_t.append(t)
    EEMD_T.append(T)
    EEMD_trend.append(trend)
    EEMD_trend_EAC.append(trend_EAC)
    EEMD_imfs.append(imfs)
    EEMD_res.append(res)
    t, S, trend_s, trend_EAC_s, imfs_s, res_s = TF.Ensemble_EMD(
        tbin_m[n],Sbin_m[n],0)
    EEMD_tS.append(t)
    EEMD_S.append(S)
    EEMD_trend_S.append(trend_s)
    EEMD_trend_EAC_S.append(trend_EAC_s)
    EEMD_imfs_S.append(imfs_s)
    EEMD_res_S.append(res_s)

a = EEMD_imfs_S[0]


EEMD_IMFS = {'IMF_1':EEMD_imfs[0],
             'IMF_2':EEMD_imfs[1],
             'IMF_3':EEMD_imfs[2]}

EEMD_IMFS_S = {'IMF_1':EEMD_imfs_S[0],
             'IMF_2':EEMD_imfs_S[1],
             'IMF_3':EEMD_imfs_S[2]}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# regrid trends on same time grid so equal size

EEMD_t_grid = []
EEMD_trend_EAC_grid = []
EEMD_trend_EAC_grid_S = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # temperature
    tt,TT = TF.bin_monthly(1953,2021,EEMD_t[n],EEMD_trend_EAC[n])
    EEMD_t_grid.append(tt)
    EEMD_trend_EAC_grid.append(TT)
    _,SS = TF.bin_monthly(1953,2021,EEMD_tS[n],EEMD_trend_EAC_S[n])
    EEMD_trend_EAC_grid_S.append(SS)





#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Comparison plots

n = 1
plt.plot(EEMD_t[n],EEMD_trend_EAC[n]/np.nanstd(EEMD_T[n]))
plt.plot(EEMD_tS[n],EEMD_trend_EAC_S[n]/np.nanstd(EEMD_S[n])*3-0.6)
plt.ylim(bottom=-3, top=4)

# def cross_correlation_using_fft(x, y):
#     f1 = fft(x)
#     f2 = fft(np.flipud(y))
#     cc = np.real(ifft(f1 * f2))
#     return fftshift(cc)

# # shift < 0 means that y starts 'shift' time steps before x # shift > 0 means that y starts 'shift' time steps after x
# def compute_shift(x, y):
#     assert len(x) == len(y)
#     c = cross_correlation_using_fft(x, y)
#     assert len(c) == len(x)
#     zero_index = int(len(x) / 2) - 1
#     shift = zero_index - np.argmax(c)
#     return shift


# result = compute_shift(EEMD_trend_EAC_grid_S[n], EEMD_trend_EAC_grid[n])

result = 6*12
x_1 = np.arange(0,len(EEMD_trend_EAC_grid[0]),1)
x_2 = np.arange(0-result,len(EEMD_trend_EAC_grid[0])-result,1)

plt.plot(x_1,EEMD_trend_EAC_grid_S[0]/np.nanstd(EEMD_S[n]))
plt.plot(x_2,EEMD_trend_EAC_grid[0]/np.nanstd(EEMD_T[n])+0.4)

# correlation without shift
a = EEMD_trend_EAC_grid_S[0]
b = EEMD_trend_EAC_grid[0]
check_nans = np.squeeze(np.logical_and([np.isfinite(a)],[np.isfinite(b)]))
stats_no_shift = sp.stats.linregress(a[check_nans],b[check_nans])
# correlation with shift
a = a[check_nans]
b = b[check_nans]
reg_stats = []
test_stats = []
tests = np.arange(0,200,1)
for n in range(len(depths)):
    a = EEMD_trend_EAC_grid_S[n]
    b = EEMD_trend_EAC_grid[n]
    check_nans = np.squeeze(np.logical_and([np.isfinite(a)],[np.isfinite(b)]))
    a = a[check_nans]
    b = b[check_nans]
    shift_stats = []
    regress_stats = []
    for n in range(len(tests)):
        stats_shift = sp.stats.linregress(a[0:np.int64(len(a)-tests[n])],
                                         b[np.int64(tests[n])::])
        shift_stats.append(stats_shift.rvalue)
        regress_stats.append(stats_shift)
    reg_stats.append(regress_stats)
    test_stats.append(shift_stats)

best_shift = tests[np.where(shift_stats == np.nanmax(shift_stats))]


# %% -----------------------------------------------------------------------------------------------
# Save data as mat file for plotting etc.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html


EEMD_t_str = []
for nn in range(len(EEMD_t)):
    ttt = EEMD_t[nn]
    a = []  
    for n in range(len(ttt)):
        tt = ttt[n]
        a.append(tt.strftime("%Y-%m-%d %H:%M:%S"))
    EEMD_t_str.append(a)
EEMD_t_str_S = []
for nn in range(len(EEMD_tS)):
    ttt = EEMD_t[nn]
    a = []  
    for n in range(len(ttt)):
        tt = ttt[n]
        a.append(tt.strftime("%Y-%m-%d %H:%M:%S"))
    EEMD_t_str_S.append(a)    


Trend_dict = {'EEMD_t': EEMD_t_str,
'EEMD_T': EEMD_T,
'EEMD_trend': EEMD_trend,
'EEMD_trend_EAC': EEMD_trend_EAC,
'EEMD_imfs': EEMD_IMFS,
'EEMD_res': EEMD_res,
'EEMD_t_S': EEMD_t_str_S,
'EEMD_T_S': EEMD_S,
'EEMD_trend_S': EEMD_trend_S,
'EEMD_trend_EAC_S': EEMD_trend_EAC_S,
'EEMD_imfs_S': EEMD_IMFS_S,
'EEMD_res_S': EEMD_res_S,
'Regression_stats': reg_stats,
'Test_stats': test_stats
}

savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "NRSPHB_Salinity_analysis.mat", Trend_dict)




