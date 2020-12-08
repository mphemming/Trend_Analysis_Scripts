

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Import modules

# general functions
import xarray as xr
# import matplotlib as mat
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


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

# 17-21 m depth
c19 = [(NRSPHB_agg.DEPTH >= 17) & (NRSPHB_agg.DEPTH <= 21)]
D19 = np.array(NRSPHB_agg.DEPTH)
D19 = D19[c19]
t19 = np.array(NRSPHB_agg.TIME)
t19 = t19[c19]
T19 = np.array(NRSPHB_agg.TEMP)
T19 = T19[c19]

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Bin the data

tbin,Tbin = TF.bin_daily(1953,2020,t19,T19)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# De-season data

# get climatology grid
clim_grid = range(0,365)
# convert times into year days
tbin_doy = []
for n in range(len(tbin)):
    t = time_grid[n]
    tbin_doy.append(int(t.strftime('%j')))

clim19 = NRSPHB_clim.TEMP_AVE[:,2] # climatology at 20m


Tbin_deseason = [None] * len(tbin_doy)
for n in range(len(tbin_doy)):
    if tbin_doy[n]-1 < 365:
        Tbin_deseason[n] = Tbin[n] - clim19[tbin_doy[n]-1]
    else:
        Tbin_deseason[n] = np.nan

# plot climatology
# plt.scatter(tbin_doy,Tbin)
# plt.plot(clim_grid,clim19,color='r')
# plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Get annual averages

# Create time grid
base = dt.datetime(1953, 6, 1)
time_grid_ann = np.array([base + dt.timedelta(days=i*365.25) for i in range(0,68)])
base = dt.datetime(1953, 1, 1)
time_grid_ann_lower = np.array([base + dt.timedelta(days=i*365.25) for i in range(0,68)])
base = dt.datetime(1954, 1, 1)
time_grid_ann_upper = np.array([base + dt.timedelta(days=i*365.25) for i in range(0,68)])


t_grid_ann = []; # convert to datetime
t_grid_ann_lower = []; # convert to datetime
t_grid_ann_upper = []; # convert to datetime

for n in range(len(time_grid_ann)):
    t_grid_ann.append(np.datetime64(time_grid_ann[n]))
    t_grid_ann_lower.append(np.datetime64(time_grid_ann_lower[n]))
    t_grid_ann_upper.append(np.datetime64(time_grid_ann_upper[n]))

# convert times into year and months
# https://strftime.org/
tbin_y = []
tbin_m = []
for n in range(len(tbin)):
    t = time_grid[n]
    tbin_y.append(int(t.strftime('%Y')))
    tbin_m.append(int(t.strftime('%m')))

un_yrs = np.unique(tbin_y)


# binning
Tbin_ann = []
tbin_ann = []
for n_bin in range(len(un_yrs)):
        c = tbin_y == un_yrs[n_bin]
        T_in_bin = []
        t_in_bin = []
        m_in_bin = []
        yr_in_bin = []
        for n in range(len(c)):
            if c[n] and np.isnan(Tbin_deseason[n]) == False:
                T_in_bin.append(Tbin_deseason[n])
                t_in_bin.append(tbin[n])
                m_in_bin.append(tbin_m[n])
                yr_in_bin.append(tbin_y[n])
                
        if len(np.unique(m_in_bin)) > 6:
            Tbin_ann.append(np.nanmean(T_in_bin))
            tbin_ann.append(dt.datetime(int(un_yrs[n_bin]),6,1))    
        else:
            Tbin_ann.append(np.nan)
            tbin_ann.append(dt.datetime(int(un_yrs[n_bin]),6,1))     
   

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
# estimate trend with sensitivity testing
# changing bin sizes, changing depth, changing start and end date




# Mann kendall tests

mk_result = mk.original_test(Tbin_ann)
mk_trend = range(0,68)*mk_result.slope + mk_result.intercept; # Theil-Sen slope

plt.scatter(tbin_ann,Tbin_ann)
plt.plot(tbin_ann,mk_trend,'r')
plt.show()

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




    
    












