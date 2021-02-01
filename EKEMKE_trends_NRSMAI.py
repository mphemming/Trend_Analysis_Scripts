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
from scipy.io import savemat
from numpy.fft import fft, ifft, fft2, ifft2, fftshift
import scipy as sp

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Load data

# mooring data
main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
CURR = pd.read_csv(main_path + 
                   'Data\\CURR_NRSMAI.csv',index_col=0)
del main_path

u = CURR.UCUR
v = CURR.VCUR
tt = CURR.index
t = []
for n in range(len(tt)):
    t.append(np.datetime64(tt[n]))
# get annual means
tbin_ann, ubin_ann = TF.bin_annually(1992,2021,t,u)
tbin_ann, vbin_ann = TF.bin_annually(1992,2021,t,v)

# calculate EKE and MKE
yr, _, _, _, _ = TF.datevec(t)
yr_tbin, _, _, _, _ = TF.datevec(tbin_ann)
u_anom = []
v_anom = []
EKE = []
MKE = []
for n in range(len(u)):
    check_yr = yr_tbin == yr[n]
    u_anom.append(u[n] - ubin_ann[check_yr])
    v_anom.append(u[n] - vbin_ann[check_yr])
    EKE.append((np.square(u_anom[n]) + np.square(v_anom[n]))/2)
    MKE.append((np.square(ubin_ann[check_yr]) + np.square(vbin_ann[check_yr]))/2)
    
EKE = np.array(EKE); MKE = np.array(MKE)
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get climatology at same depths

# calculate simple climatology for now
EKE_clim = TF.calc_clim_monthly(t,EKE)  

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% -----------------------------------------------------------------------------------------------
# Get monthly EKE

t_m,EKE_m = TF.bin_monthly(1953,2021,t,EKE)

# %% -----------------------------------------------------------------------------------------------
# Ensemble EMD
print('Running Ensemble EMD')

# NRSMAI

EEMD_t, EEMD_T, EEMD_trend, EEMD_trend_EAC, EEMD_imfs, EEMD_res = \
                                                TF.Ensemble_EMD(t_m,EKE_m,0)

EEMD_t_str = []
a = []  
for n in range(len(EEMD_t)):
    tt = EEMD_t[n]
    a.append(tt.strftime("%Y-%m-%d %H:%M:%S"))
    EEMD_t_str.append(a)


EKE_EEMD = {'EEMD_t':EEMD_t_str,
            'EEMD_T':EEMD_T,
            'EEMD_trend':EEMD_trend,
            'EEMD_trend_EAC':EEMD_trend_EAC,
            'EEMD_imfs':EEMD_imfs}


savemat("C:\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\Data\\" + 
        "NRSMAI_EKE_analysis.mat", EKE_EEMD)














