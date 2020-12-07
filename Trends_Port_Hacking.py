

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Import modules

# general functions
import xarray as xr
import matplotlib as mat
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
# For mann-kendall test and innovative Sen slope analysis
# https://pypi.org/project/pymannkendall/
import pymannkendall as mk
# For pettitt test - need to check results as function copied from stackoverflow
import pettitt as pett
import pyhomogeneity as hg

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Load data

main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
NRSPHB_clim = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_Climatology_1953-2019_BottleCTDMooringSatellite.nc')
NRSPHB_agg = xr.open_dataset(main_path + 'Data\\Port_Hacking_TEMP_1953-2019_aggregated.nc')

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


# Create time grid
base = dt.datetime(1953, 1, 1,12,0,0)
time_grid = np.array([base + dt.timedelta(days=i) for i in range(1,24836)])
base = dt.datetime(1953, 1, 1,0,0,0)
time_grid_lower = np.array([base + dt.timedelta(days=i) for i in range(1,24836)])
base = dt.datetime(1953, 1, 2,0,0,0)
time_grid_upper = np.array([base + dt.timedelta(days=i) for i in range(1,24836)])

t_grid = []; # convert to datetime
t_grid_lower = []; # convert to datetime
t_grid_upper = []; # convert to datetime
for n in range(len(time_grid)-1):
    t_grid.append(np.datetime64(time_grid[n]))
    t_grid_lower.append(np.datetime64(time_grid_lower[n]))
    t_grid_upper.append(np.datetime64(time_grid_upper[n]))


# binning
Tbin = []
tbin = []

for n_bin in range(0,24833):
    
    c = [(t19 >= t_grid_lower[n_bin]) & (t19 <= t_grid_upper[n_bin])]
    T_med = np.median(T19[c])  
    tbin.append(t_grid[n_bin])
    Tbin.append(T_med)
    
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


for day in range(clim_grid):
    c = []
    for n in range(len(tbin_doy)):
            c.append(tbin_doy[n] == day+1)
        


# plot climatology
# plt.scatter(tbin_doy,Tbin)
# plt.plot(clim_grid,clim19,color='r')
# plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# %% -----------------------------------------------------------------------------------------------
# Mann kendall tests

mk_result = mk.original_test(Tbin)

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
pett_result = hg.pettitt_test(Tbin)


trend_change_date = tbin[pett_result[1]]



    
    












