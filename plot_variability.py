
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
# Import plot functions from script
import Trends_plot_functions as TP
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
import seaborn as sns

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% load data

main_path = "\\Users\\mphem\\Documents\\Work\\UNSW\\Trends\\"
NRSPHB_agg = xr.open_dataset(main_path + 'Data\\PH100_TEMP_1953-2020_aggregated_v1.nc')
NRSMAI_agg = xr.open_dataset(main_path + 'Data\\MAI090_TEMP_1944-2020_aggregated_v1.nc')

# %% -----------------------------------------------------------------------------------------------
# Select data at specific depths

print('Selecting data at different depths:')

depths_PHB = [2, 19, 31, 40, 50, 77, 99]
depths_MAI = [2, 20, 50]

D_PHB = []; D_MAI = []
t_PHB = []; t_MAI = []
T_PHB = []; T_MAI = []
# ---------------------------------------------
# NRSPHB
for n in range(len(depths_PHB)):
    print(str(depths_PHB[n]) + ' m')
    # index check
    c = [(NRSPHB_agg.DEPTH_AGG >= depths_PHB[n] - 3) & 
         (NRSPHB_agg.DEPTH_AGG <= depths_PHB[n] + 3)]
    # Depth
    d = np.array(NRSPHB_agg.DEPTH_AGG);
    D_PHB.append(d[c])
    # time
    tt = np.array(NRSPHB_agg.TIME);
    t_PHB.append(tt[c])    
    # Temp
    TT = np.array(NRSPHB_agg.TEMP_AGG);
    T_PHB.append(TT[c])       
# ---------------------------------------------
# NRSMAI
for n in range(len(depths_MAI)):
    print(str(depths_MAI[n]) + ' m')
    # index check
    c = [(NRSMAI_agg.DEPTH_AGG >= depths_MAI[n] - 3) & 
         (NRSMAI_agg.DEPTH_AGG <= depths_MAI[n] + 3)]
    # Depth
    d = np.array(NRSMAI_agg.DEPTH_AGG);
    D_MAI.append(d[c])
    # time
    tt = np.array(NRSMAI_agg.TIME);
    t_MAI.append(tt[c])    
    # Temp
    TT = np.array(NRSMAI_agg.TEMP_AGG);
    T_MAI.append(TT[c])         


del c, d, tt, TT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# %% get standard deviations

t_std_PHB,std_PHB = TP.get_variability_site(NRSPHB_agg.TIME, 
                                            NRSPHB_agg.DEPTH_AGG, 
                                            NRSPHB_agg.TEMP_AGG, 
                                            depths_PHB)
t_std_MAI,std_MAI = TP.get_variability_site(NRSMAI_agg.TIME, 
                                            NRSMAI_agg.DEPTH_AGG, 
                                            NRSMAI_agg.TEMP_AGG, 
                                            depths_MAI)

for n in range(np.size(std_PHB,1)):
    plt.plot(t_std_PHB,std_PHB[:,n])

for n in range(np.size(std_MAI,1)):
    plt.plot(t_std_MAI,std_MAI[:,n])

# Create figure including binned profile and 
# standard deviations and std in time at each site

depth_PHB_for_MAI = [2, 20, 50]

bin_mean_PHB, _, bin_std_PHB, bin_D_PHB, bin_n_PHB = TF.bin_profile(
    NRSPHB_agg.TEMP_AGG,NRSPHB_agg.DEPTH_AGG,depth_PHB_for_MAI,2)
bin_mean_MAI, _, bin_std_MAI, bin_D_MAI, bin_n_MAI = TF.bin_profile(
                NRSMAI_agg.TEMP_AGG,NRSMAI_agg.DEPTH_AGG,depths_MAI,2)




# Seaborn violin plot

sns.set_theme()

# Create a random dataset across several variables
rs = np.random.default_rng(0)
n, p = 40, 8
d = rs.normal(0, 2, (n, p))
d += np.log(np.arange(1, p + 1)) * -5 + 10

# Show each distribution with both violins and points
sns.violinplot(data=d, palette="light:g", inner="points", orient="h")










