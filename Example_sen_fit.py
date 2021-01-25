import numpy as np

t = np.arange(1,127,1)
y = t + 2.5


a = 2.5
s = 0.25
n = 126
data = s*t + a
y1 = data[0:63]
y2 = data[63:126]


# below using equations from Sen 2017 - Innovative trend significance test nd applications
# slope
s_calc = 2*(np.nanmean(y2)-np.nanmean(y1))/n
# intercept
t_bar = np.nanmean(t)
y_bar = np.interp(t_bar,t,data)

y2_bar = np.nanmean(y2) 
y1_bar = np.nanmean(y1)

a_calc = y_bar - (2 *(y2_bar - y1_bar) / n) * t_bar

# reconstruct fit
fit_recon = s_calc*t + a_calc
