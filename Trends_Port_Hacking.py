

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
# wavelet analysis
import pycwt as wavelet
from pycwt.helpers import find
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
Tbin_deseason = np.array(TF.deseason(tbin,Tbin,climsurf))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get annual averages

# Dont' need to use annual means
# tbin_ann,Tbin_ann = TF.bin_annually(1953,2020,tbin,Tbin_deseason)
   
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Determine relationship and effect of climate indices 

# reg = linear_model.LinearRegression()
# reg.fit(tbin_ann[0:-1], PDO_IND)

# Tbin_reg = np.array(Tbin_ann[0:-1])
# c = np.logical_not(np.isnan(Tbin_reg))
# Tbin_reg = Tbin_reg[c]
# PDO_IND_reg =  PDO_IND[c]

# PDO_result = sp.stats.linregress(PDO_IND_reg,Tbin_reg); # no real relationship

        
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

trend_change_points = 0
ITA_stats = TF.ITA(tbin,Tbin_deseason,trend_change_points)


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




    
    












