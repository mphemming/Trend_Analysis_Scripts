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
from statsmodels.tsa.stattools import kpss

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Functions

# -----------------------------------------------------------------------------------------------
# datevec function
 
def datevec(TIME):
    """
    Convert array of datetime64 to a calendar array of year, month, day, hour,
    minute, seconds, microsecond with these quantites indexed on the last axis.

    Parameters
    ----------
    TIME : datetime64 array (...)
        numpy.ndarray of datetimes of arbitrary shape

    Returns
    -------
    cal : uint32 array (..., 7)
        calendar array with last axis representing year, month, day, hour,
        minute, second, microsecond
    """

    # allocate output 
    out = np.empty(TIME.shape + (7,), dtype="u4")
    # decompose calendar floors
    Y, M, D, h, m, s = [TIME.astype(f"M8[{x}]") for x in "YMDhms"]
    out[..., 0] = Y + 1970 # Gregorian Year
    out[..., 1] = (M - Y) + 1 # month
    out[..., 2] = (D - M) + 1 # dat
    out[..., 3] = (TIME - D).astype("m8[h]") # hour
    out[..., 4] = (TIME - h).astype("m8[m]") # minute
    out[..., 5] = (TIME - m).astype("m8[s]") # second
    out[..., 6] = (TIME - s).astype("m8[us]") # microsecond
    
    yr = out[:,0]; mn = out[:,1]; dy = out[:,2]; hr = out[:,3]; 
    yday = []
    for n in range(len(yr)):
        yday.append(dt.date(yr[n], mn[n], dy[n]).timetuple().tm_yday)

    return yr, mn, dy, hr, yday


# -----------------------------------------------------------------------------------------------
# Bin data daily


def bin_daily(start_yr,end_yr,TIME,TEMP):

    # Create time grids
    base = dt.datetime(start_yr, 1, 1,12,0,0)
    time_grid = np.array([base + dt.timedelta(days=i) for i in range(1,30000)])
    base = dt.datetime(start_yr, 1, 1,0,0,0)
    time_grid_lower = np.array([base + dt.timedelta(days=i) for i in range(1,30000)])
    base = dt.datetime(start_yr, 1, 2,0,0,0)
    time_grid_upper = np.array([base + dt.timedelta(days=i) for i in range(1,30000)])
    
    # convert to datetime
    t_grid = []; 
    t_grid_lower = []; 
    t_grid_upper = []; 
    for n in range(len(time_grid)-1):
        t_grid.append(np.datetime64(time_grid[n]))
        t_grid_lower.append(np.datetime64(time_grid_lower[n]))
        t_grid_upper.append(np.datetime64(time_grid_upper[n]))
    
    # flag time grid past end year date
    t_grid = np.array(t_grid)
    t_grid_lower = np.array(t_grid_lower)
    t_grid_upper = np.array(t_grid_upper)
    
    c = np.where(t_grid <= np.datetime64(dt.datetime(end_yr,12,31)))
    t_grid = t_grid[c]
    c = np.where(t_grid_lower <= np.datetime64(dt.datetime(end_yr,12,30)))
    t_grid_lower = t_grid_lower[c]    
    c = np.where(t_grid_upper <= np.datetime64(dt.datetime(end_yr+1,1,1)))
    t_grid_upper = t_grid_upper[c]    
    
    
    # binning
    Tbin = []
    tbin = []
    
    for n_bin in range(len(t_grid)):
        
        c = [(TIME >= t_grid_lower[n_bin]) & (TIME <= t_grid_upper[n_bin])]
        T_med = np.median(TEMP[c])  
        tbin.append(t_grid[n_bin])
        Tbin.append(T_med)
    
    return tbin, Tbin
    

# -----------------------------------------------------------------------------------------------
# De-season data


def deseason(TIME,TEMP,clim):

    # get climatology grid
    clim_grid = range(0,365)
    # get year days
    _,_,_,_,yday_bin = datevec(np.array(TIME))
    # de-season temperatures
    
    TEMP_deseason = [None] * len(yday_bin)
    for n in range(len(yday_bin)):
        if yday_bin[n]-1 < 365:
            TEMP_deseason[n] = TEMP[n] - clim[yday_bin[n]-1]
        else:
            TEMP_deseason[n] = np.nan
    
    # plot climatology
    plt.scatter(yday_bin,TEMP)
    plt.plot(clim_grid,clim,color='r')
    plt.show()
    # plot de-seasoned TEMP 
    # plt.scatter(yday_bin,TEMP_deseason)
    
    return TEMP_deseason


# -----------------------------------------------------------------------------------------------
# Bin data annually


def bin_annually(start_yr,end_yr,TIME,TEMP):

    # Create time grids
    base = dt.datetime(start_yr, 6, 1)
    time_grid_ann = np.array([base + dt.timedelta(days=i*365.25) for i in range(0,68)])
    base = dt.datetime(start_yr, 1, 1)
    time_grid_ann_lower = np.array([base + dt.timedelta(days=i*365.25) for i in range(0,68)])
    base = dt.datetime(start_yr, 1, 1)
    time_grid_ann_upper = np.array([base + dt.timedelta(days=i*365.25) for i in range(0,68)])
    
    t_grid_ann = []; # convert to datetime
    t_grid_ann_lower = []; # convert to datetime
    t_grid_ann_upper = []; # convert to datetime
    
    for n in range(len(time_grid_ann)):
        t_grid_ann.append(np.datetime64(time_grid_ann[n]))
        t_grid_ann_lower.append(np.datetime64(time_grid_ann_lower[n]))
        t_grid_ann_upper.append(np.datetime64(time_grid_ann_upper[n]))
    
    # flag time grid past end year date
    t_grid_ann = np.array(t_grid_ann)
    t_grid_ann_lower = np.array(t_grid_ann_lower)
    t_grid_ann_upper = np.array(t_grid_ann_upper)
    
    c = np.where(t_grid_ann <= np.datetime64(dt.datetime(end_yr,12,31)))
    t_grid_ann = t_grid_ann[c]
    c = np.where(t_grid_ann_lower <= np.datetime64(dt.datetime(end_yr,12,30)))
    t_grid_ann_lower = t_grid_ann_lower[c]    
    c = np.where(t_grid_ann_upper <= np.datetime64(dt.datetime(end_yr+1,1,1)))
    t_grid_ann_upper = t_grid_ann_upper[c]    
    
    # get months and years of TEMP
    yrs,mns,_,_,_ = datevec(np.array(TIME))
    # get unique years
    un_yrs = np.unique(yrs)
    
    # binning
    Tbin_ann = []
    tbin_ann = []
    for n_bin in range(len(un_yrs)):
            c = yrs == un_yrs[n_bin]
            T_in_bin = []
            t_in_bin = []
            m_in_bin = []
            yr_in_bin = []
            for n in range(len(c)):
                if c[n] and np.isnan(TEMP[n]) == False:
                    T_in_bin.append(TEMP[n])
                    t_in_bin.append(TIME[n])
                    m_in_bin.append(mns[n])
                    yr_in_bin.append(yrs[n])
            # bin only if there are data from 6 or more months       
            if len(np.unique(m_in_bin)) > 6:
                Tbin_ann.append(np.nanmean(T_in_bin))
                tbin_ann.append(dt.datetime(int(un_yrs[n_bin]),6,1))    
            else:
                Tbin_ann.append(np.nan)
                tbin_ann.append(dt.datetime(int(un_yrs[n_bin]),6,1))   
    
    return tbin_ann, Tbin_ann

# -----------------------------------------------------------------------------------------------
# KPSS test

def kpss_test(series, **kw):   
    
    # reference
    # https://www.machinelearningplus.com/time-series/kpss-test-for-stationarity/
    # remove NaNs from series
    series = np.array(series)
    series = series[np.where(np.logical_not(np.isnan(series)))]
    
    statistic, p_value, n_lags, critical_values = kpss(series, **kw)
    
    # Format Output
    print(f'KPSS Statistic: {statistic}')
    print(f'p-value: {p_value}')
    print(f'num lags: {n_lags}')
    print('Critial Values:')
    for key, value in critical_values.items():
        print(f'   {key} : {value}')
    print(f'Result: The series is {"not " if p_value < 0.05 else ""}stationary')

