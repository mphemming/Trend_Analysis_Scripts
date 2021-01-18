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
# EEMD 
from PyEMD import EMD, EEMD, Visualisation
# significance stuff
from sklearn.metrics import mean_squared_error
import signalz
import statsmodels.api as sm
import time

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
    
    # If not a datetime64, convert from dt.datetime
    if '64' not in str(type(TIME)):
        t = []
        for nt in range(len(TIME)):
            o = TIME[nt]
            if '64' not in str(type(o)):
                t.append(np.datetime64(o.strftime("%Y-%m-%dT%H:%M:%S")))
            else:
                t.append(o)
        TIME = np.array(t)
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
        
    tbin = np.array(tbin)
    Tbin = np.array(Tbin)
    
    return tbin, Tbin
    
# %% -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# De-season data


def deseason(TIME,TEMP,clim):

    if np.size(clim) > 12:
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
    else:
        # get climatology grid
        clim_grid = np.arange(1,13,1)
        # get months
        _,m,_,_,_ = datevec(np.array(TIME))
        # de-season temperatures
        TEMP_deseason = [None] * len(m)
        TEMP_deseason = np.array(TEMP_deseason)
        for n in clim_grid:
            check = m == n
            TEMP_deseason[check] = TEMP[check] - clim[n-1]
    
    # plot climatology
    # plt.scatter(yday_bin,TEMP)
    # plt.plot(clim_grid,clim,color='r')
    # plt.show()
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
# Bin data monthly


def bin_monthly(start_yr,end_yr,TIME,TEMP):
    # create ranges for loop
    yrs_range = np.arange(start_yr,end_yr+1,1)
    mns_range = np.arange(1,13,1)
    # get years and months from time
    yrs,mns,_,_,_ = datevec(np.array(TIME))
    

    t_mns = []
    T_mns = []
    for yr in yrs_range:
        for mn in mns_range:
            t_mns.append(dt.datetime(yr,mn,15))
            check_bin = np.logical_and([yrs == yr], [mns == mn])
            T = TEMP[np.squeeze(check_bin)]
            if np.size(T) > 0:
                T_mns.append(np.nanmean(np.float32(TEMP[np.squeeze(check_bin)])))
            else:
                T_mns.append(np.nan)
            
    t_mns = np.array(t_mns); T_mns = np.array(T_mns);     
    return t_mns, T_mns


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
    if p_value < 0.05:
        result = 'Not Stationary'
    else:
        result = 'Stationary'    
    
    class results:
        KPSS_statistic = statistic
        KPSS_p_value = p_value
        KPSS_n_lags = n_lags
        KPSS_critical_values = critical_values
        KPSS_result = result
        
    return results


# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# ITA function

def ITA(TIME,TEMP,trend_change_points,figure):
    
    # 1:1 line (no trend line)
    line = np.arange(start=-10, stop=10, step=0.1)
    # separate the period into two equal parts
    split_number = np.int(np.round(len(TIME)/2)-1)
    T1 = np.array(np.sort(TEMP[0:split_number]))
    T2 = np.array(np.sort(TEMP[split_number+1::]))
    t1 = np.array(TIME[0:split_number])
    t2 = np.array(TIME[split_number+1::])   
    # check if segments are equal
    if len(t1) != len(t2):
        t2 = t2[0:-1]
        T2 = T2[0:-1]
    # sort two parts in ascending order
    T1_indices = np.argsort(T1) 
    T1_t = t1[T1_indices]
    T2_indices = np.argsort(T2) 
    T2_t = t2[T2_indices]
    check_nans = np.squeeze(np.logical_and(\
            [np.isfinite(T1)], [np.isfinite(T2)]))
    T1_nonan = T1[check_nans]
    T2_nonan = T2[check_nans]    
    if figure == 1:
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

        if figure == 1:
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
        if figure == 1:
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
        # trend_mean = np.mean([np.mean(high_res),np.mean(low_res)])
        date_range = np.nanmax(t2)-np.nanmin(t1)
        days = date_range.astype('timedelta64[D]')
        # trend_ave_per_decade = (trend_mean/np.int64(days))*3653
        

    else:
        change_point_1 = trend_change_points[0]
        change_point_2 = trend_change_points[1]
        if figure == 1:
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
        if figure == 1:
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
        date_range = np.nanmax(t2)-np.nanmin(t1)
        days = date_range.astype('timedelta64[D]')
        # trend_ave_per_decade = (trend_mean/np.int64(days))*3653
        
    # code from R implementation of ITA
    # https://rdrr.io/cran/trendchange/src/R/innovtrend.R    
    # See Sen 2017 Global Warming paper for steps and equations
    # trend
    yrs,_,_,_,_ = datevec(T1_t)
    T1_min_year = np.nanmin(yrs)
    T1_max_year = np.nanmax(yrs)
    yrs,_,_,_,_ = datevec(T2_t)
    T2_min_year = np.nanmin(yrs)
    T2_max_year = np.nanmax(yrs)      
    yr_range = (T2_max_year-T1_min_year)
    slope_sen = ((2*(np.nanmean(T2)-np.nanmean(T1)))/np.int64(days))
    trend_sen_period = slope_sen*np.int64(days)
    trend_sen_per_decade = slope_sen*3653
    trend_sen_period_numb_days = np.int64(days)
    # get intercept
    tt = np.interp(np.nanmean(T1),T1,T2)
    
    intercept = np.nanmean(T2) - \
        ((2*(np.nanmean(T2)-np.nanmean(T1)))/np.int64(days)) * tt 
    
    # Calculating slope standard deviation
    # covariance with NaNs
    df = pd.DataFrame({'T1':T2_nonan,'T2':T1_nonan})
    cov_vals = df.corr()
    slope_std = (2 * np.sqrt(2) ) * np.nanstd(TEMP) \
        *np.sqrt( 1 - np.array(cov_vals) ) \
            / np.int64(days) / np.sqrt(np.int64(days))
    slope_std = slope_std[0,1]
    # Trend indicator calculation
    D = np.nanmean((T2-T1)*10/np.nanmean(T1))
    # Confidence limits (CL) of the trend slope at 95 percent
    # See Sen 2017 bok for more details page 205 and around
    CLlower95 = 0 - 1.96*slope_std
    CLupper95 = 0 + 1.96*slope_std
    # Significant trend decision
    if np.abs(CLupper95) < np.abs(slope_sen):
        significance = 'significant trend (95%)'
    else:
        significance = 'Not significant (95%)'

    # save information
    class ITA_stats: 
        
        ITA_trend_sen_period = trend_sen_period
        ITA_trend_sen_period_numb_days = trend_sen_period_numb_days
        ITA_trend_sen_per_decade = trend_sen_per_decade
        ITA_slope_sen = slope_sen
        ITA_slope_std = slope_std
        ITA_trend_indicator = D
        ITA_slope_lower_conf = CLlower95
        ITA_slope_upper_conf = CLupper95
        ITA_significance = significance
        TEMP_half_1 = T1_nonan
        TEMP_half_2 = T2_nonan
        t_start_half_1 = T1_min_year
        t_end_half_1 = T1_max_year
        t_start_half_2 = T2_min_year
        t_end_half_2 = T2_max_year
        
    return ITA_stats

# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# ITA significance

def ITA_significance(TIME,TEMP,ACF_result,numb_sims):
    # determine leakage to get closest matching brownian noise signal to TEMP
    tests = np.arange(0,1,0.02)
    ACF_tests = []
    RMSE_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=tests[n], start=0, std=0.8, source="gaussian")
        ACF_tests.append(pd.Series(sm.tsa.acf(x, nlags=50)))
        A = ACF_tests[n]
        RMSE_tests.append(np.sqrt(mean_squared_error(ACF_result[0:50], A[0:50])))
    leakage = np.float(tests[RMSE_tests == np.min(RMSE_tests)])
    
    # determine standard deviation closest to reality
    real_std = np.nanstd(TEMP)
    tests = np.arange(0.1,2,0.1)
    std_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=leakage, start=0, std=tests[n], source="gaussian")
        x_std = np.nanstd(x)
        std_tests.append(real_std-x_std)
    std_chosen = np.float(tests[np.abs(std_tests) == np.nanmin(np.abs(std_tests))])  
        
        
    # Code for checking if red noise similar
    x = signalz.brownian_noise(len(TEMP), leak=leakage, start=0, std=std_chosen, source="gaussian")
    plt.plot(ACF_result)
    plt.plot(pd.Series(sm.tsa.acf(x, nlags=50)))
    plt.show()
     
    x_sims = []
    stats = []
    for n in range(0,numb_sims):
        tic = time.perf_counter()
        print('Simulation: ' + str(n))
        x_sims.append(signalz.brownian_noise(len(TEMP), leak=leakage, start=0, \
                                     std=std_chosen, source="gaussian"))
        stats.append(ITA(TIME,x_sims[n],0,0))
        toc = time.perf_counter()
        print(f"{toc - tic:0.4f} seconds")
        
    # calculate mean of slopes and mean std of slopes
    slopes = []
    std_slopes = []
    conf_slopes = []
    for n in range(0,numb_sims):
        tt = stats[n]
        slopes.append(tt.ITA_slope_sen)
        std_slopes.append(tt.ITA_slope_std)
        conf_slopes.append(tt.ITA_slope_upper_conf)

    conf_95_limit = np.nanmean(conf_slopes)

    return conf_95_limit


# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# EEMD function

def Ensemble_EMD(TIME,TEMP,figure):
    
    # https://github.com/laszukdawid/PyEMD
    # https://pyemd.readthedocs.io/en/latest/intro.html#
    
    # Ensemble empirical mode decomposition (EEMD) [Wu2009] is noise-assisted technique, 
    # which is meant to be more robust than simple Empirical Mode Decomposition (EMD). The 
    # robustness is checked by performing many decompositions on signals slightly perturbed 
    # from their initial position. In the grand average over all IMF results the noise will 
    # cancel each other out and the result is pure decomposition.
    
    # ensure is numpy array
    TIME = np.array(TIME)
    TEMP = np.array(TEMP)
    # remove nans from timeseries
    check_nans = np.isfinite(TEMP)
    T = TEMP[check_nans]
    t = TIME[check_nans]
    # perform EEMD
    eemd = EEMD(noise_width = 0.2, trials=500, parallel=True, 
                processes=5, max_imfs=4,include_residue=False) # same parameters as GMSL trends Chen et al. paper and almost same as Wu et al nature trends paper
    eemd.eemd(T)
    imfs, res = eemd.get_imfs_and_residue()
    if figure == 1:
        # visualise imfs
        vis = Visualisation()
        vis.plot_imfs(imfs=imfs, residue=res, t=t, include_residue=True)
        # vis.plot_instant_freq(np.arange(1,len(t)+1,1), imfs=imfs)
        vis.show()
    # reconstruct timeseries from imfs
    n_imfs = np.size(imfs,0)
    for n in range(n_imfs):
        if n == 0:
            recon =  imfs[n,:]
        else:
            recon = recon + imfs[n,:]
            
    if n_imfs >= 9:        
        # construct trend using last 3 imfs
        trend = imfs[n_imfs-3,:] + imfs[n_imfs-2,:] + imfs[n_imfs-1,:]
    if n_imfs == 8:        
        # construct trend using last 3 imfs
        trend = imfs[n_imfs-2,:] + imfs[n_imfs-1,:]    
        
    if figure == 1:
        # create trend figure
        plt.figure()
        plt.plot(t,T)
        plt.plot(t,trend,'k')
        plt.show()
    
    if 'trend' not in (locals()):
        trend = 0
    
    return t, T, trend, imfs, res


# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# EEMD quicker function

def Ensemble_EMD_quick(TIME,TEMP):
    
    # https://github.com/laszukdawid/PyEMD
    # https://pyemd.readthedocs.io/en/latest/intro.html#
    
    # Ensemble empirical mode decomposition (EEMD) [Wu2009] is noise-assisted technique, 
    # which is meant to be more robust than simple Empirical Mode Decomposition (EMD). The 
    # robustness is checked by performing many decompositions on signals slightly perturbed 
    # from their initial position. In the grand average over all IMF results the noise will 
    # cancel each other out and the result is pure decomposition.
    
    # ASSUMES THAT DATA HAS NO NANS AND IS NUMPY ARRAY

    # perform EEMD
    eemd = EEMD(noise_width = 0.2, trials=500, parallel=True, processes=5) # same parameters as GMSL trends Chen et al. paper and almost same as Wu et al nature trends paper
    eemd.eemd(TEMP)
    imfs, res = eemd.get_imfs_and_residue()
 
    # reconstruct timeseries from imfs
    n_imfs = np.size(imfs,0)
    for n in range(n_imfs):
        if n == 0:
            recon =  imfs[n,:]
        else:
            recon = recon + imfs[n,:]
            
    if n_imfs >= 9:        
        # construct trend using last 3 imfs
        trend = imfs[n_imfs-3,:] + imfs[n_imfs-2,:] + imfs[n_imfs-1,:]
    if n_imfs == 8:        
        # construct trend using last 3 imfs
        trend = imfs[n_imfs-2,:] + imfs[n_imfs-1,:]    

    if 'trend' not in (locals()):
        trend = 0

    return trend, imfs



# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# EEMD significance function using red noise simulations

def EEMD_significance(TIME,TEMP,ACF_result,numb_sims):

    # determine leakage to get closest matching brownian noise signal to TEMP
    tests = np.arange(0,1,0.02)
    ACF_tests = []
    RMSE_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=tests[n], start=0, std=0.8, source="gaussian")
        ACF_tests.append(pd.Series(sm.tsa.acf(x, nlags=10)))
        A = ACF_tests[n]
        RMSE_tests.append(np.sqrt(mean_squared_error(ACF_result[0:10], A[0:10])))
    leakage = np.float(tests[RMSE_tests == np.min(RMSE_tests)])
    
    # determine standard deviation closest to reality
    real_std = np.nanstd(TEMP)
    tests = np.arange(0.1,2,0.1)
    std_tests = []
    for n in range(len(tests)):
        x = signalz.brownian_noise(len(TEMP), leak=leakage, start=0, std=tests[n], source="gaussian")
        x_std = np.nanstd(x)
        std_tests.append(real_std-x_std)
    std_chosen = np.float(tests[np.abs(std_tests) == np.nanmin(np.abs(std_tests))])  
        
    
    # Code for checking if red noise similar
    # x = signalz.brownian_noise(len(TEMP), leak=0.44, start=0, std=1, source="gaussian")
    # plt.plot(ACF_result)
    # plt.plot(pd.Series(sm.tsa.acf(x, nlags=10)))
    # plt.show()

    x_sims = []
    trend_sims = []
    for n in range(0,numb_sims):
        tic = time.perf_counter()
        print('Simulation: ' + str(n))
        x_sims.append(signalz.brownian_noise(len(TEMP), leak=leakage, start=0, \
                                             std=std_chosen, source="gaussian"))
        tr,_ = Ensemble_EMD_quick(TIME,x_sims[n])
        trend_sims.append(tr)
        toc = time.perf_counter()
        print(f"{toc - tic:0.4f} seconds")
        
    # combine trends and calculate standard deviation of trends
    # get standard deviation
    std_array = []
    for n in range(len(TEMP)):
        array_for_stats = []
        for nn in range(len(trend_sims)):
            tt = trend_sims[nn]
            array_for_stats.append(tt[n] -tt[0])
        std_array.append(np.std(array_for_stats))
    conf_std_limit = (std_array * (np.ones(len(std_array))*1.96))
    
    return conf_std_limit, std_array, trend_sims, x_sims




# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# test effect of trials n on trend


def test_trials(TIME,TEMP):
    
    trial_tests = np.arange(50,1000,50)
    
    trend_results = []
    for n in range(len(trial_tests)):
        print(trial_tests[n])
        # perform EEMD
        eemd = EEMD(noise_width = 0.2, trials=trial_tests[n], parallel=True, processes=10) # same parameters as GMSL trends Chen et al. paper and almost same as Wu et al nature trends paper
        eemd.eemd(TEMP)
        imfs, res = eemd.get_imfs_and_residue()
        # get trends
        n_imfs = np.size(imfs,0)
        if n_imfs == 9:        
            # construct trend using last 3 imfs
            trend = imfs[n_imfs-3,:] + imfs[n_imfs-2,:] + imfs[n_imfs-1,:]
        if n_imfs == 8:        
            # construct trend using last 3 imfs
            trend = imfs[n_imfs-2,:] + imfs[n_imfs-1,:]   
        if n_imfs == 7:        
            # construct trend using last 3 imfs
            trend = imfs[n_imfs-2,:] + imfs[n_imfs-1,:]      
            
        trend_results.append(trend)
    
    for n in range(len(trend_results)):
        tt = np.array(trend_results[n])
        plt.scatter(x=trial_tests[n],y=tt[0],color='r')
        plt.scatter(x=trial_tests[n],y=tt[300],color='b')
        plt.scatter(x=trial_tests[n],y=tt[500],color='g')
        plt.scatter(x=trial_tests[n],y=tt[700],color='k')
    plt.show()   
    
    for n in range(len(trend_results)):
        if n < 10:
            plt.plot(TIME,trend_results[n],'r')
        else:
            plt.plot(TIME,trend_results[n],'k')
    plt.show()      
        
        

# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Simple climatology function for IMOS mooring timeseries

def calc_clim_monthly(TIME,TEMP):
    
    yr, mn, dy, hr, yday = datevec(TIME)
    clim = []
    for n_mon in np.arange(1,13,1):
        check = mn == n_mon
        clim.append(np.nanmean(TEMP[check]))
    
    return clim

# %% ----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Convert to MATLAB datenum

def datetime2matlabdn(d):
   mdn = d + dt.timedelta(days = 366)
   frac_seconds = (d-dt.datetime(d.year,d.month,d.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   frac_microseconds = d.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
   return mdn.toordinal() + frac_seconds + frac_microseconds



