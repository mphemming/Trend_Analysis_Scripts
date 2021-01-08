

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

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Load data

# mooring data
main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
NRSPHB_clim = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_Climatology_1953-2019_BottleCTDMooringSatellite.nc')
NRSPHB_agg = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_1953-2019_aggregated.nc')
# PDO data
PDO = pd.read_csv(main_path + 'Data\\PDO.csv',index_col=0)
PDO_IND = np.array(PDO.ind)
PDO_IND_t = np.array(PDO.index)
c = PDO_IND_t > 1952
PDO_IND = PDO_IND[c]
PDO_IND_t = PDO_IND_t[c]
# SAM data
SAM = pd.read_csv(main_path + 'Data\\SAM.csv',index_col=0) 
# need to add code here to sort out SAM

del c, main_path

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

# surface
csurf = [(NRSPHB_agg.DEPTH >= 0) & (NRSPHB_agg.DEPTH <= 2)]
Dsurf = np.array(NRSPHB_agg.DEPTH)
Dsurf = Dsurf[csurf]
tsurf = np.array(NRSPHB_agg.TIME)
tsurf = tsurf[csurf]
Tsurf = np.array(NRSPHB_agg.TEMP)
Tsurf = Tsurf[csurf]

del csurf

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Bin the data

# This is done to get a regular time grid with daily resolution
tbin,Tbin = TF.bin_daily(1953,2020,tsurf,Tsurf)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

# select climatology at similar depth
climsurf = NRSPHB_clim.TEMP_AVE[:,0] # climatology at 20m
# get de-seasoned temperatures
Tbin_deseason = np.array(TF.deseason(tbin,Tbin,climsurf))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get monthly averages

# Using de-seasoned timeseries
tbin_m,Tbin_m = TF.bin_monthly(1953,2021,tbin,Tbin_deseason)
   
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')
t, T, trend, imfs, res = TF.Ensemble_EMD(tbin_m,Tbin_m)

# Autocorrelation analysis
check = np.where(np.isnan(Tbin_m))
ACF1 = pd.Series(sm.tsa.acf(Tbin_m[20:401], nlags=10)); # where longest streak without nans
ACF2 = pd.Series(sm.tsa.acf(Tbin_m[671:777], nlags=10)); # where longest streak without 
ACF_result = []
for n in range(0,10):
    ACF_result.append(np.nanmean([ACF1[n],ACF2[n]]))
ACF_result = np.array(ACF_result)

# significance
conf_std_limit, std_array, trend_sims, x_sims =  TF.EEMD_significance(tbin_m,Tbin_m,ACF_result,200)

# Create figure
plt.figure(figsize=(15,8))
for n in range(1,200):
    tt = trend_sims[n]
    plt.plot(tbin_m,tt-tt[0],color='grey')
plt.plot(tbin_m,conf_std_limit,color='r')
plt.plot(tbin_m,conf_std_limit*-1,color='r')
plt.plot(t,trend-trend[0],color='k',linewidth=2)
plt.xlabel('Year')
plt.ylabel(r'$\rmTemperature Trend [^\circ C]$')
plt.show()



# %% -----------------------------------------------------------------------------------------------
# Pettitt tests

pett_result = hg.pettitt_test(Tbin_m)
trend_change_date = tbin_m[pett_result[1]]

# %% -----------------------------------------------------------------------------------------------
# KPSS test to check for stationarity
# If the result = 'Not stationary', a deterministc trend / linear regression is not suitable

TF.kpss_test(Tbin_m)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# estimate trend with sensitivity testing
# changing bin sizes, changing depth, changing start and end date


# Mann kendall tests


mk_result = mk.trend_free_pre_whitening_modification_test(Tbin_m)
mk_trend = range(len(tbin_m))*mk_result.slope + mk_result.intercept;

plt.plot(tbin_m,Tbin_m)
plt.plot(tbin_m,mk_trend,'r')
plt.plot(t,trend,'g')
plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Save data as mat file for plotting etc.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html





#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Innovative trend analysis

# trend_change_points = 0
# ITA_stats = TF.ITA(tbin_m,Tbin_m,trend_change_points)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Statistical significance

# run Monte Carlo simulations

# test = np.random.normal(0,0.25/np.sqrt(len(T)),len(T))

# sim = []
# for x in test:
#     sim.append(T[0]*x)
    
    
# _,_, trend_test,_,_ = TF.Ensemble_EMD(t,test)


# random walk simulation

# seed(1)
# random_walk = list()
# random_walk.append(-1 if random() < 0.5 else 1)
# for i in range(1, len(T)):
# 	movement = -1 if random() < 0.5 else 1
# 	value = random_walk[i-1] + movement
# 	random_walk.append(value)
# plt.plot(random_walk)
# plt.show()

# ACF_test = pd.Series(sm.tsa.acf(random_walk, nlags=10))
    

# https://matousc89.github.io/signalz/sources/generators/brownian_noise.html

# x = signalz.brownian_noise(len(T), leak=0.55, start=0, std=1, source="gaussian")
# ACF_test = pd.Series(sm.tsa.acf(x, nlags=10))
# _,_, trend_test,_,_ = TF.Ensemble_EMD(t,x)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>










