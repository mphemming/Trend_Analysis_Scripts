

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Import modules

from rich.traceback import install
install()
import copy
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
import time
import random
import signalz
from scipy.io import savemat

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Load data

# mooring data
main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
NRSPHB_clim = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_Climatology_1953-2019_BottleCTDMooringSatellite.nc')
NRSPHB_agg = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_1953-2019_aggregated.nc')

del main_path

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

print('Selecting data at different depths:')

depths = [2]

D = []
t = []
T = []
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # index check
    c = [(NRSPHB_agg.DEPTH >= depths[n] - 2) & (NRSPHB_agg.DEPTH <= depths[n] + 2)]
    # Depth
    d = np.array(NRSPHB_agg.DEPTH);
    D.append(d[c])
    # time
    tt = np.array(NRSPHB_agg.TIME);
    t.append(tt[c])    
    # Temp
    TT = np.array(NRSPHB_agg.TEMP);
    T.append(TT[c])       

# plt.plot(t[0],T[0])
# plt.plot(t[10],T[10])
# plt.show()

del c, d, tt, TT

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Bin the data
# This is done to get a regular time grid with daily resolution

T = np.array(T)
t = np.array(t)

tbin, Tbin = TF.bin_daily(1953,2020,t,T)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get climatology at same depths


clim = np.ones((365,1),dtype=float)
for day in range(0,365):
    day_temps = NRSPHB_clim.TEMP_AVE[day,:] 
    day_std = NRSPHB_clim.TEMP_STD[day,:] 
    clim[day,:] = np.interp(depths,NRSPHB_clim.DEPTH,day_temps)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# De-season data

print('Removing the season')

#get de-seasoned temperatures

cl = clim
Tbin_deseason = np.array(TF.deseason(tbin,Tbin,cl))
    
del cl

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# prepare timeseries comparisons

# gappy timeseres, also with chunks at beginning, middle, and end
c = np.array(tbin) < np.datetime64('1965-01-01')
g_start = copy.copy(np.array(Tbin_deseason)); g_start[c] = np.nan;
c = np.squeeze(np.logical_and([np.array(tbin) > np.datetime64('1980-01-01')],
                   [np.array(tbin) < np.datetime64('1990-01-01')]))
g_mid = copy.copy(np.array(Tbin_deseason)); g_mid[c] = np.nan; 
c = np.array(tbin) > np.datetime64('2010-01-01')
g_end = copy.copy(np.array(Tbin_deseason)); g_end[c] = np.nan;
# Different resolutions
tt = tbin; TT = Tbin_deseason
tt_m,TT_m = TF.bin_monthly(1953,2021,tt,TT)
tt_a,TT_a = TF.bin_annually(1953,2021,tt,TT)
# Filled timeseries
_,TT_filled_m,_ = TF.fill_gaps(tt_m,TT_m,np.squeeze(clim[:,0]),30*12);
_,TT_filled_a,_ = TF.fill_gaps(tt_a,TT_a,np.squeeze(clim[:,0]),30);
# with random missing data different percentages
r_10 = random.sample(range(1, len(TT_filled_m)), np.int64(len(TT_filled_m)*0.9))
r_25 = random.sample(range(1, len(TT_filled_m)), np.int64(len(TT_filled_m)*0.75))
r_50 = random.sample(range(1, len(TT_filled_m)), np.int64(len(TT_filled_m)*0.5))
tt_m_10, TT_filled_m_10 = TF.bin_monthly(1953,2021,tt_m[r_10],TT_filled_m[r_10])
tt_m_25, TT_filled_m_25 = TF.bin_monthly(1953,2021,tt_m[r_25],TT_filled_m[r_25])
tt_m_50, TT_filled_m_50 = TF.bin_monthly(1953,2021,tt_m[r_50],TT_filled_m[r_50])

# change number of ensembles
# ens_1 = 500
# ens_2 = 1000
# ens_3 = 2000
# ens_4 = 3000

# save as a class
class data:
    gaps = np.array(Tbin_deseason)
    gaps_t = np.array(tbin)    
    gaps_start = np.array(g_start)
    gaps_middle = np.array(g_mid)
    gaps_end = np.array(g_end)
    daily = TT
    daily_t = tt
    monthly = TT_m
    monthly_filled = TT_filled_m
    monthly_filled_10 = TT_filled_m_10
    monthly_filled_10_t = tt_m_10
    monthly_filled_25 = TT_filled_m_25
    monthly_filled_25_t = tt_m_25
    monthly_filled_50 = TT_filled_m_50
    monthly_filled_50_t = tt_m_50
    monthly_t = tt_m
    annually = TT_a
    annually_t = tt_a
    annually_filled = TT_filled_a
    
    
del g_start, g_mid, g_end, tt, TT, tt_m, TT_m, tt_a, TT_a

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get comparisons

EEMD_t = []
EEMD_T = []
EEMD_trend = []
EEMD_imfs = []
# gaps
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.gaps_t,data.gaps,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# gaps_start
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.gaps_t,data.gaps_start,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# gaps_middle
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.gaps_t,data.gaps_middle,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# gaps_end
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.gaps_t,data.gaps_end,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# monthly
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.monthly_t,data.monthly,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# monthly filled
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.monthly_t,data.monthly_filled,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# monthly filled 10% missing
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.monthly_filled_10_t,data.monthly_filled_10,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# monthly filled 25% missing
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.monthly_filled_25_t,data.monthly_filled_25,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# monthly filled 50% missing
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.monthly_filled_50_t,data.monthly_filled_50,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# annually
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.annually_t,data.annually,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);
# annually filled
t, T, tr, _, imfs, _ = TF.Ensemble_EMD(data.annually_t,data.annually_filled,0)
EEMD_t.append(t); EEMD_T.append(T); EEMD_trend.append(tr), EEMD_imfs.append(imfs);

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Create plot comparisons

# Resolution effect
a = EEMD_imfs[0]
plt.plot(EEMD_t[0],a[10]+a[11])
a = EEMD_imfs[4]
plt.plot(EEMD_t[4],a[6]+a[7])
a = EEMD_imfs[9]
plt.plot(EEMD_t[9],a[3]+a[4])
# gaps chunks
a = EEMD_imfs[0]
plt.plot(EEMD_t[0],a[10]+a[11])
plt.plot(EEMD_t[1],EEMD_trend[1])
a = EEMD_imfs[2]
plt.plot(EEMD_t[2],a[10]+a[11])
a = EEMD_imfs[3]
plt.plot(EEMD_t[3],a[10]+a[9])
# gaps random
a = EEMD_imfs[5]
plt.plot(EEMD_t[5],a[8]+a[7]+a[6])
a = EEMD_imfs[6]
plt.plot(EEMD_t[6],a[6]+a[7])
a = EEMD_imfs[7]
plt.plot(EEMD_t[7],a[6]+a[7])
a = EEMD_imfs[8]
plt.plot(EEMD_t[8],a[6]+a[7])
# gap filling (Monthly)
# plt.plot(tbin,Tbin_deseason,'.')
a = EEMD_imfs[4]
plt.plot(EEMD_t[4],a[6]+a[7])
a = EEMD_imfs[5]
plt.plot(EEMD_t[5],a[8]+a[7]+a[6])
# gap filling (Annually)
# plt.plot(tbin,Tbin_deseason,'.')
a = EEMD_imfs[6]
plt.plot(EEMD_t[6],a[7]+a[6])
a = EEMD_imfs[7]
plt.plot(EEMD_t[7],a[7]+a[6])

# take homes from these experiments:
    # Annual measurements should not be used for this anlysis
    # Resolution has an impact on the trends
    # annual (because of gap-filling) or daily data 
    # (because of gaps) should not be used!
    # Gaps in the middle have the least effect than start and end, 
    # timeseries missing data at the beginning or end should not be used!
    # Chunk gaps have more of an effect than random gaps, anything more
    # than 25% not great!
    
# %%-----------------------------------------------------------------------------
# Explore changes in trend when changing number of ensembles

a = EEMD_imfs[7]; trend = a[7];
noise_1 = signalz.brownian_noise(len(trend), leak=0.2, start=0, std=0.8, source="gaussian")
noise_2 = signalz.brownian_noise(len(trend), leak=0.9, 
                                 start=0, std=0.8, source="gaussian")*0.5
example = trend + noise_1 + noise_2

tests = [50,500,1000]
# EEMD_imfs = []
EEMD_mean = []
EEMD_std = []
for n in range(len(tests)):
    print(tests[n])
    ensemble = np.ones((np.int32(len(example)),30))
    for ens in range(0,30):
        eemd = TF.EEMD(noise_width = 0.2, trials=tests[n], parallel=True, 
                        processes=8, max_imfs=4,include_residue=False)
        eemd.eemd(example)
        imfs, _ = eemd.get_imfs_and_residue()
        tr = imfs[-1]
        if tr[-1]-tr[0] < 0.1:
            tr = imfs[-2]+imfs[-1]
        # if tr[-1]-tr[0] < 0.5:
        #     tr = a[-2]+a[-1]+a[-3]       
        ensemble[:,ens] = tr 
    EEMD_mean.append(np.nanmean(ensemble,1))
    EEMD_std.append(np.nanstd(ensemble,1))
    
# trend_tests = []    
# for n in range(len(tests)):
#     a = EEMD_tests[n]
#     tr = a[-1]
#     if tr[-1]-tr[0] < 0.1:
#         tr = a[-2]+a[-1]
#     # if tr[-1]-tr[0] < 0.5:
#     #     tr = a[-2]+a[-1]+a[-3]       
#     trend_tests.append(tr)

# create plot 
handles = []
for n in range(len(tests)):
    a = EEMD_mean[n]
    h = plt.plot(a-a[0])
    plt.plot((a-a[0])+EEMD_std[n],'k',linestyle='dotted')
    plt.plot((a-a[0])-EEMD_std[n],'k',linestyle='dotted')
    handles.append(h)
plt.plot(trend-trend[0],'k')
plt.legend([str(tests[0]),str(tests[1]),str(tests[2]),'trend'])

# create plot 
handles = []
for n in range(len(tests)):
    a = EEMD_std[n]
    h = plt.plot(a)
    handles.append(h)
plt.legend([str(tests[0]),str(tests[1]),str(tests[2]),'trend'])


# take home messages so far:
    # Trend is no reconstructed properly
    # next time, run multiple times to get uncertainty and means
    
    

