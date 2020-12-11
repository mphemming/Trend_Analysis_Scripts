

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
# 1D interpolation
from scipy.interpolate import interp1d


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

# 17-21 m depth
csurf = [(NRSPHB_agg.DEPTH >= 0) & (NRSPHB_agg.DEPTH <= 2)]
Dsurf = np.array(NRSPHB_agg.DEPTH)
Dsurf = Dsurf[csurf]
tsurf = np.array(NRSPHB_agg.TIME)
tsurf = tsurf[csurf]
Tsurf = np.array(NRSPHB_agg.TEMP)
Tsurf = Tsurf[csurf]

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
Tbin_deseason = TF.deseason(tbin,Tbin,climsurf)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get annual averages

tbin_ann,Tbin_ann = TF.bin_annually(1953,2020,tbin,Tbin_deseason)
   
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Determine relationship and effect of climate indices 

# reg = linear_model.LinearRegression()
# reg.fit(tbin_ann[0:-1], PDO_IND)

Tbin_reg = np.array(Tbin_ann[0:-1])
c = np.logical_not(np.isnan(Tbin_reg))
Tbin_reg = Tbin_reg[c]
PDO_IND_reg =  PDO_IND[c]

PDO_result = sp.stats.linregress(PDO_IND_reg,Tbin_reg); # no real relationship

            
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# KPSS test to check for stationarity
# If the result = 'Not stationary', a deterministc trend / linear regression is not suitable

TF.kpss_test(Tbin_deseason)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Innovative trend analysis

# 1:1 line (no trend line)
line = np.arange(start=-1.5, stop=2.5, step=0.1)

# separate the period into two
T1 = np.array(np.sort(Tbin_ann[0:33]))
T2 = np.array(np.sort(Tbin_ann[33:-2]))

T1_indices = np.argsort(Tbin_ann[0:33]) 
tb = np.array(tbin_ann[0:33])
T1_t = tb[T1_indices]
T2_indices = np.argsort(Tbin_ann[33:-2]) 
tb = np.array(tbin_ann[33:-2])
T2_t = tb[T2_indices]

# get 'low' mean slope
cl = [(T1 < -0.88) & (T2 < -0.5)]
low_stats = sp.stats.linregress(T1[cl],T2[cl])
low_line_points = np.interp(T1[cl],line,line)

# get 'medium' mean slope
cm = [(T1 > -0.88) & (T1 < -0.25) & (T2 > -0.5) & (T2 < 0.45)]
med_stats = sp.stats.linregress(T1[cm],T2[cm])
med_line_points = np.interp(T1[cm],line,line)

# get 'high' mean slope
ch = [(T1 >= -0.25) & (T2 >= 0.45)]
high_stats = sp.stats.linregress(T1[ch],T2[ch])
high_line_points = np.interp(T1[ch],line,line)

# Determine trends
#------------------------
ll = T2[cl]
low_res = []
for n in range(len(ll)):
    low_res.append(np.abs(ll[n] - low_line_points[n]))
#------------------------
lm = T2[cm]
med_res = []
for n in range(len(lm)):
    med_res.append(np.abs(lm[n] - med_line_points[n]))
#------------------------
lh = T2[ch]
high_res = []
for n in range(len(lh)):
    high_res.append(np.abs(lh[n] - high_line_points[n]))
    
    
# code from R implementation of ITA
# https://rdrr.io/cran/trendchange/src/R/innovtrend.R    
# See Sen 2017 Global Warming paper for steps and equations
# trend
trend =  ((2*(np.nanmean(T2)-np.nanmean(T1)))/66)*100 
# Calculating slope standard deviation
ssd = (2*np.sqrt(2))*np.nanstd(Tbin_ann)*np.sqrt(1-np.cov(T2,T1))/66/np.sqrt(66)
# Trend indicator calculation
D = np.nanmean((T2-T1)*10/np.nanmean(T1))
# Confidence limits (CL) of the trend slope at 95 percent
CLlower95 = 0 - 1.96*ssd
CLupper95 = 0 + 1.96*ssd

trend_mean = np.mean([np.mean(high_res),np.mean(med_res),np.mean(low_res)])
# trend_century = ((1-(len(T1)+len(T2))/100) * trend_mean_slope) + trend_mean_slope 
trend_century = (trend_mean/66)*100

a= plt.figure()
axes = a.add_axes([0.1,0.1,0.8,0.8])
axes.scatter(T1,T2)
axes.scatter(T1[cl],T2[cl])
axes.scatter(T1[cm],T2[cm])
axes.plot(T1[cl],low_stats.slope*T1[cl]+low_stats.intercept,'k')
axes.plot(T1[cm],med_stats.slope*T1[cm]+med_stats.intercept,'k')
axes.plot(T1[ch],high_stats.slope*T1[ch]+high_stats.intercept,'k')
axes.set_xlim([-1.5, 2.5])
axes.set_ylim([-1.5, 2.5])
axes.set_xlabel('1953 - 1985')
axes.set_ylabel('1986 - 2018')
axes.plot(line,line,'k')
plt.show()


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# estimate trend with sensitivity testing
# changing bin sizes, changing depth, changing start and end date


# Mann kendall tests

mk_result = mk.original_test(Tbin_ann)
mk_result_wang = mk.pre_whitening_modification_test(Tbin_ann)
mk_trend = range(0,68)*mk_result.slope + mk_result.intercept; # Theil-Sen slope
mk_trend_wang = range(0,68)*mk_result_wang.slope + mk_result_wang.intercept; # Theil-Sen slope

plt.scatter(tbin_ann,Tbin_ann)
plt.plot(tbin_ann,mk_trend,'r')
plt.plot(tbin_ann,mk_trend_wang,'g')
plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Pettitt tests

# remove NaNs from timeseries for pettitt test
# Tbin_no_nan = []
# tbin_no_nan = []
# for n in range(len(Tbin)):
#     if np.logical_not(np.isnan(Tbin[n])):
#         Tbin_no_nan.append(Tbin[n])
#         tbin_no_nan.append(tbin[n])

# pett_result = pett.pettitt_test(Tbin_no_nan)
pett_result = hg.pettitt_test(Tbin_ann)


trend_change_date = tbin_ann[pett_result[1]]




    
    












