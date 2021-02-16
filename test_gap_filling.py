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
import random 
import signalz
from scipy.io import savemat
from sklearn.metrics import mean_squared_error
import scipy.stats as sp

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
for n in range(len(depths)):
    print(str(depths[n]) + ' m')
    # This is done to get a regular time grid with daily resolution
    tt,TT = TF.bin_daily(2011,2021,t[n],np.float64(T[n]))
    tbin.append(tt)
    Tbin.append(TT)
    
del tt, TT, n  

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Create artificial data gaps

TT = Tbin[0]; tt = tbin[0];

# n% random missing data
# 10%
r_10 = random.sample(range(1, len(TT)), np.int64(len(TT)*0.9))
tt_10perc,T_ = TF.bin_daily(2011,2021,tt[r_10],TT[r_10])
_,TT_10perc,_ = TF.fill_gaps(tt_10perc,T_,np.squeeze(clim_daily[0]),2*365)
tt_10perc,TT_10perc = TF.bin_daily(2011,2021,tt_10perc,TT_10perc)
# 20%
r_20 = random.sample(range(1, len(TT)), np.int64(len(TT)*0.8))
tt_20perc,T_ = TF.bin_daily(2011,2021,tt[r_20],TT[r_20])
_,TT_20perc,_ = TF.fill_gaps(tt_20perc,T_,np.squeeze(clim_daily[0]),2*365)
tt_20perc,TT_20perc = TF.bin_daily(2011,2021,tt_20perc,TT_20perc)
# 30%
r_30 = random.sample(range(1, len(TT)), np.int64(len(TT)*0.7))
tt_30perc,T_ = TF.bin_daily(2011,2021,tt[r_30],TT[r_30])
_,TT_30perc,_ = TF.fill_gaps(tt_30perc,T_,np.squeeze(clim_daily[0]),2*365)
tt_30perc,TT_30perc = TF.bin_daily(2011,2021,tt_30perc,TT_30perc)
# 40%
r_40 = random.sample(range(1, len(TT)), np.int64(len(TT)*0.6))
tt_40perc,T_ = TF.bin_daily(2011,2021,tt[r_40],TT[r_40])
_,TT_40perc,_ = TF.fill_gaps(tt_40perc,T_,np.squeeze(clim_daily[0]),2*365)
tt_40perc,TT_40perc = TF.bin_daily(2011,2021,tt_40perc,TT_40perc)
# 50%
r_50 = random.sample(range(1, len(TT)), np.int64(len(TT)*0.5))
tt_50perc,T_ = TF.bin_daily(2011,2021,tt[r_50],TT[r_50])
_,TT_50perc,_ = TF.fill_gaps(tt_50perc,T_,np.squeeze(clim_daily[0]),2*365)
tt_50perc,TT_50perc = TF.bin_daily(2011,2021,tt_50perc,TT_50perc)

# # with data gaps
# c1 = np.logical_and([tt > np.datetime64('2016-05-01')],[tt < np.datetime64('2017-05-01')])
# # mixture of data gaps
# # one gap
# a = np.squeeze(c1 == False)
# _,TT_1gap,_ = TF.fill_gaps(tt[a],TT[a],np.squeeze(clim_daily[0]),2*365)
# tt_1gap,TT_1gap = TF.bin_daily(2011,2021,tt[a],TT_1gap)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Calculate statistics 

# 10%
a = np.squeeze(np.logical_and([np.isfinite(TT)],[np.isfinite(TT_10perc)]))
diff = TT-TT_10perc
b = diff != 0
a = np.squeeze(np.logical_and([a],[b]))
diff = diff[a]
RMSE_10 = mean_squared_error(TT[a],TT_10perc[a],squared=False)
diff = TT[a]-TT_10perc[a]; diff = diff[diff > 0]
AVE_10 = np.nanmean(diff); std_10 = np.nanstd(diff)
result_10 = sp.linregress(TT[a],TT_10perc[a])
a_10 = a
# 20%
a = np.squeeze(np.logical_and([np.isfinite(TT)],[np.isfinite(TT_20perc)]))
diff = TT-TT_20perc
b = diff != 0
a = np.squeeze(np.logical_and([a],[b]))
diff = diff[a]
RMSE_20 = mean_squared_error(TT[a],TT_20perc[a],squared=False)
diff = TT[a]-TT_20perc[a]; diff = diff[diff > 0]
AVE_20 = np.nanmean(diff); std_20 = np.nanstd(diff)
result_20 = sp.linregress(TT[a],TT_20perc[a])
a_20 = a
# 30%
a = np.squeeze(np.logical_and([np.isfinite(TT)],[np.isfinite(TT_30perc)]))
diff = TT-TT_30perc
b = diff != 0
a = np.squeeze(np.logical_and([a],[b]))
diff = diff[a]
RMSE_30 = mean_squared_error(TT[a],TT_30perc[a],squared=False)
diff = TT[a]-TT_30perc[a]; diff = diff[diff > 0]
AVE_30 = np.nanmean(diff); std_30 = np.nanstd(diff)
result_30 = sp.linregress(TT[a],TT_30perc[a])
a_30 = a
# 40%
a = np.squeeze(np.logical_and([np.isfinite(TT)],[np.isfinite(TT_40perc)]))
diff = TT-TT_40perc
b = diff != 0
a = np.squeeze(np.logical_and([a],[b]))
diff = diff[a]
RMSE_40 = mean_squared_error(TT[a],TT_40perc[a],squared=False)
diff = TT[a]-TT_40perc[a]; diff = diff[diff > 0]
AVE_40 = np.nanmean(diff); std_40 = np.nanstd(diff)
result_40 = sp.linregress(TT[a],TT_40perc[a])
a_40 = a
# 50%
a = np.squeeze(np.logical_and([np.isfinite(TT)],[np.isfinite(TT_50perc)]))
diff = TT-TT_50perc
b = diff != 0
a = np.squeeze(np.logical_and([a],[b]))
diff = diff[a]
RMSE_50 = mean_squared_error(TT[a],TT_50perc[a],squared=False)
diff = TT[a]-TT_50perc[a]; diff = diff[diff > 0]
AVE_50 = np.nanmean(diff); std_50 = np.nanstd(diff)
result_50 = sp.linregress(TT[a],TT_50perc[a])
a_50 = a
# average RMSE
AVE_RMSE = np.nanmean([RMSE_10,RMSE_20,RMSE_30,RMSE_40,RMSE_50])
AVE_diff = np.nanmean([AVE_10,AVE_20,AVE_30,AVE_40,AVE_50])
AVE_std = np.nanmean([std_10,std_20,std_30,std_40,std_50])
AVE_rval = np.nanmean(
    [result_10.rvalue,result_20.rvalue,result_30.rvalue,
     result_40.rvalue,result_50.rvalue,])


# %% -----------------------------------------------------------------------------------------------
# Create Figure

plt.plot(TT[a_50],TT_50perc[a_50],'.')
plt.plot(TT[a_40],TT_40perc[a_40],'.')
plt.plot(TT[a_30],TT_30perc[a_30],'.')
plt.plot(TT[a_20],TT_20perc[a_20],'.')
plt.plot(TT[a_10],TT_10perc[a_10],'.')
plt.plot(np.linspace(13,25,10),np.linspace(13,25,10),color='k')
plt.plot(np.linspace(13,25,10),np.linspace(13,25,10)+AVE_RMSE,color='k',linestyle=':')
plt.plot(np.linspace(13,25,10),np.linspace(13,25,10)-AVE_RMSE,color='k',linestyle=':')
plt.legend(['10 %','20 %', '30 %','40 %','50 %'])
plt.xlabel('Real Temperature [°C]')
plt.ylabel('Filled Temperature [°C]')
plt.savefig('C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\' + 
            'Plots\\gap_filling.png', dpi=300)



