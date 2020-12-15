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

# %% -----------------------------------------------------------------------------------------------
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
    
# %% -----------------------------------------------------------------------------------------------
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

# %% -----------------------------------------------------------------------------------------------
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

# %% -----------------------------------------------------------------------------------------------
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



# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# ITA function

def ITA(TIME,TEMP,trend_change_points):
    
    # 1:1 line (no trend line)
    line = np.arange(start=-10, stop=10, step=0.1)
    # separate the period into two equal parts
    split_number = np.int(np.round(len(TIME)/2)-1)
    T1 = np.array(np.sort(TEMP[0:split_number]))
    T2 = np.array(np.sort(TEMP[split_number+1::]))
    t1 = np.array(TIME[0:split_number])
    t2 = np.array(TIME[split_number+1::])    
    # sort two parts in ascending order
    T1_indices = np.argsort(T1) 
    T1_t = t1[T1_indices]
    T2_indices = np.argsort(T2) 
    T2_t = t2[T2_indices]
    check_nans = np.squeeze(np.logical_and(\
            [np.isfinite(T1)], [np.isfinite(T2)]))
    T1_nonan = T1[check_nans]
    T2_nonan = T2[check_nans]    
    # create first plot for selecting trend points
    a= plt.figure()
    axes = a.add_axes([0.1,0.1,0.8,0.8])
    axes.scatter(T1,T2)    
    axes.set_xlim([np.nanmin([np.nanmin(T1),np.nanmin(T2)])-0.5, \
                   np.nanmax([np.nanmax(T1),np.nanmax(T2)])+0.5])
    axes.set_ylim([np.nanmin([np.nanmin(T1),np.nanmin(T2)])-0.5, \
                   np.nanmax([np.nanmax(T1),np.nanmax(T2)])+0.5])
    yrs,_,_,_,_ = datevec(T1_t)
    T1_min_year = np.nanmin(yrs)
    T1_max_year = np.nanmax(yrs)
    yrs,_,_,_,_ = datevec(T2_t)
    T2_min_year = np.nanmin(yrs)
    T2_max_year = np.nanmax(yrs)    
    axes.set_xlabel(str(T1_min_year) + ' - ' + str(T1_max_year))
    axes.set_ylabel(str(T2_min_year) + ' - ' + str(T2_max_year))
    axes.plot(line,line,'k') 
    if np.size(trend_change_points) == 1:

        plt.close(fig=a)
        # get trends
        # trend 1
        check_low = T1_nonan <= np.float(trend_change_points)
        low_stats = sp.stats.linregress(T1_nonan[check_low],T2_nonan[check_low])
        low_line_points = np.interp(T1_nonan[check_low],line,line)        
        # trend 2
        check_high = T1_nonan > np.float(trend_change_points)
        high_stats = sp.stats.linregress(T1_nonan[check_high],T2_nonan[check_high])
        high_line_points = np.interp(T1_nonan[check_high],line,line)   
        # create new plot
        a= plt.figure()
        axes = a.add_axes([0.1,0.1,0.8,0.8])
        axes.scatter(T1,T2)        
        axes.scatter(T1_nonan[check_low],T2_nonan[check_low])
        axes.scatter((T1_nonan[check_high]),T2_nonan[check_high])
        axes.plot(T1_nonan[check_low],low_stats.slope*T1_nonan[check_low]+low_stats.intercept,'k')
        axes.plot(T1_nonan[check_high],high_stats.slope*T1_nonan[check_high]+high_stats.intercept,'k')
        axes.set_xlim([np.nanmin([np.nanmin(T1),np.nanmin(T2)])-0.5, \
                       np.nanmax([np.nanmax(T1),np.nanmax(T2)])+0.5])
        axes.set_ylim([np.nanmin([np.nanmin(T1),np.nanmin(T2)])-0.5, \
                       np.nanmax([np.nanmax(T1),np.nanmax(T2)])+0.5])
        yrs,_,_,_,_ = datevec(T1_t)
        T1_min_year = np.nanmin(yrs)
        T1_max_year = np.nanmax(yrs)
        yrs,_,_,_,_ = datevec(T2_t)
        T2_min_year = np.nanmin(yrs)
        T2_max_year = np.nanmax(yrs)    
        axes.set_xlabel(str(T1_min_year) + ' - ' + str(T1_max_year))
        axes.set_ylabel(str(T2_min_year) + ' - ' + str(T2_max_year))
        axes.plot(line,line,'k')
        plt.show()
        
        # Determine trends
        #------------------------
        ll = T2_nonan[check_low]
        low_res = []
        for n in range(len(ll)):
            low_res.append(np.abs(ll[n] - low_line_points[n]))
        #------------------------
        lh = T2_nonan[check_high]
        high_res = []
        for n in range(len(lh)):
            high_res.append(np.abs(lh[n] - high_line_points[n]))        
        trend_mean = np.mean([np.mean(high_res),np.mean(low_res)])
        trend_ave_per_century = (trend_mean/66)*100
        

    else:
        change_point_1 = trend_change_points[0]
        change_point_2 = trend_change_points[1]
        plt.close(fig=a)
        # get trends
        # trend 1
        check_low = T1_nonan <= np.float(change_point_1)
        low_stats = sp.stats.linregress(T1_nonan[check_low],T2_nonan[check_low])
        low_line_points = np.interp(T1_nonan[check_low],line,line)        
        # trend 2
        check_med = np.logical_and(\
                [T1_nonan > np.float(change_point_1)], [T1_nonan <= np.float(change_point_2)])
        check_med = np.squeeze(check_med)
        med_stats = sp.stats.linregress(T1_nonan[check_med],T2_nonan[check_med])
        med_line_points = np.interp(T1_nonan[check_med],line,line) 
        # trend 3
        check_high = T1_nonan > np.float(change_point_2)
        high_stats = sp.stats.linregress(T1_nonan[check_high],T2_nonan[check_high])
        high_line_points = np.interp(T1_nonan[check_high],line,line)   
        # create new plot
        a= plt.figure()
        axes = a.add_axes([0.1,0.1,0.8,0.8])
        axes.scatter(T1,T2)        
        axes.scatter(T1_nonan[check_low],T2_nonan[check_low])
        axes.scatter(T1_nonan[check_med],T2_nonan[check_med])
        axes.scatter((T1_nonan[check_high]),T2_nonan[check_high])
        axes.plot(T1_nonan[check_low],low_stats.slope*T1_nonan[check_low]+low_stats.intercept,'k')
        axes.plot(T1_nonan[check_med],med_stats.slope*T1_nonan[check_med]+med_stats.intercept,'k')
        axes.plot(T1_nonan[check_high],high_stats.slope*T1_nonan[check_high]+high_stats.intercept,'k')
        axes.set_xlim([np.nanmin([np.nanmin(T1),np.nanmin(T2)])-0.5, \
                       np.nanmax([np.nanmax(T1),np.nanmax(T2)])+0.5])
        axes.set_ylim([np.nanmin([np.nanmin(T1),np.nanmin(T2)])-0.5, \
                       np.nanmax([np.nanmax(T1),np.nanmax(T2)])+0.5])
        yrs,_,_,_,_ = datevec(T1_t)
        T1_min_year = np.nanmin(yrs)
        T1_max_year = np.nanmax(yrs)
        yrs,_,_,_,_ = datevec(T2_t)
        T2_min_year = np.nanmin(yrs)
        T2_max_year = np.nanmax(yrs)    
        axes.set_xlabel(str(T1_min_year) + ' - ' + str(T1_max_year))
        axes.set_ylabel(str(T2_min_year) + ' - ' + str(T2_max_year))
        axes.plot(line,line,'k')
        plt.show()
        
        # Determine trends
        #------------------------
        ll = T2_nonan[check_low]
        low_res = []
        for n in range(len(ll)):
            low_res.append(np.abs(ll[n] - low_line_points[n]))
        #------------------------
        lm = T2_nonan[check_med]
        med_res = []
        for n in range(len(lm)):
            med_res.append(np.abs(lm[n] - med_line_points[n]))
        #------------------------
        lh = T2_nonan[check_high]
        high_res = []
        for n in range(len(lh)):
            high_res.append(np.abs(lh[n] - high_line_points[n]))
        trend_mean = np.mean([np.mean(high_res),np.mean(med_res),np.mean(low_res)])
        trend_ave_per_century = (trend_mean/66)*100
        
    # code from R implementation of ITA
    # https://rdrr.io/cran/trendchange/src/R/innovtrend.R    
    # See Sen 2017 Global Warming paper for steps and equations
    # trend
    yr_range = (T2_max_year-T1_min_year)
    slope_sen = ((2*(np.nanmean(T2)-np.nanmean(T1)))/yr_range)
    trend_sen = slope_sen*100 
    # Calculating slope standard deviation
    # covariance with NaNs
    df = pd.DataFrame({'T1':T1,'T2':T2})
    cov_vals = df.corr()
    slope_std = (2*np.sqrt(2))*np.nanstd(TEMP)*np.sqrt(1-np.array(cov_vals))/yr_range/np.sqrt(yr_range)
    slope_std = slope_std[0,1]
    # Trend indicator calculation
    D = np.nanmean((T2-T1)*10/np.nanmean(T1))
    # Confidence limits (CL) of the trend slope at 95 percent
    CLlower95 = 0 - 1.96*slope_std
    CLupper95 = 0 + 1.96*slope_std
    
    # save information
    class ITA_stats: 
        
        ITA_trend_ave_per_century = trend_ave_per_century
        ITA_slope_sen = slope_sen
        ITA_slope_std = slope_std
        ITA_trend_sen = trend_sen
        ITA_trend_indicator = D
        ITA_slope_lower_conf = CLlower95
        ITA_slope_upper_conf = CLupper95
        
    return ITA_stats



