from PyEMD import EMD, EEMD, Visualisation
import numpy as np

# https://github.com/laszukdawid/PyEMD
# https://pyemd.readthedocs.io/en/latest/intro.html#

# Ensemble empirical mode decomposition (EEMD) [Wu2009] is noise-assisted technique, 
# which is meant to be more robust than simple Empirical Mode Decomposition (EMD). The 
# robustness is checked by performing many decompositions on signals slightly perturbed 
# from their initial position. In the grand average over all IMF results the noise will 
# cancel each other out and the result is pure decomposition.

check_nans = np.isfinite(Tbin_m)
Tbin_deseason_nonans = Tbin_m[check_nans]
# N = len(Tbin_m)
# t = np.arange(0, N) * dt + t0
t_nonans =tbin_m[check_nans]

# Extract imfs and residue
# In case of EMD
eemd = EEMD(noise_width = 0.2, trials=1000) # same parameters as GMSL trends Chen et al. paper
eemd.eemd(Tbin_deseason_nonans)
imfs, res = eemd.get_imfs_and_residue()

# In general:
#components = EEMD()(S)
#imfs, res = components[:-1], components[-1]

vis = Visualisation()
vis.plot_imfs(imfs=imfs, residue=res, t=t_nonans, include_residue=True)
vis.plot_instant_freq(t_nonans, imfs=imfs)
vis.show()

# re-create the signal from IMFs

recon =  imfs[0,:] + imfs[1,:] + imfs[2,:] + imfs[3,:] + imfs[4,:] + imfs[5,:] \
    + imfs[6,:] + imfs[7,:] + imfs[8,:]
trend = imfs[6,:] + imfs[7,:] + imfs[8,:]
multi_decadal = imfs[8,:] + imfs[7,:] + imfs[6,:]

a= plt.figure()
axes = a.add_axes([0.1,0.1,0.8,0.8])
plt.plot(t_nonans,Tbin_deseason_nonans)
plt.plot(t_nonans,trend,'k')
# plt.plot(t_nonans,multi_decadal,'g')
plt.plot(tbin_m,mk_trend,'r')
plt.show()

# 


