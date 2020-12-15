from PyEMD import EMD, Visualisation
import numpy as np

# https://github.com/laszukdawid/PyEMD
# https://pyemd.readthedocs.io/en/latest/intro.html#


check_nans = np.isfinite(Tbin_deseason)
Tbin_deseason_nonans = Tbin_deseason[check_nans]
N = len(Tbin_deseason_nonans)
t = np.arange(0, N) * dt + t0
t_nonans =t[check_nans]

# Extract imfs and residue
# In case of EMD
emd = EMD()
emd.emd(Tbin_deseason_nonans)
imfs, res = emd.get_imfs_and_residue()

# In general:
#components = EEMD()(S)
#imfs, res = components[:-1], components[-1]

vis = Visualisation()
vis.plot_imfs(imfs=imfs, residue=res, t=t_nonans, include_residue=True)
vis.plot_instant_freq(t_nonans, imfs=imfs)
vis.show()


# re-create the signal from IMFs

test =  imfs[0,:] + imfs[1,:] + imfs[2,:] + imfs[3,:] + imfs[4,:] + imfs[5,:] + imfs[6,:] + imfs[7,:] + \
     imfs[8,:] + imfs[9,:] 

test_multi_decadal = imfs[9,:]+imfs[8,:]

a= plt.figure()
axes = a.add_axes([0.1,0.1,0.8,0.8])
plt.scatter(t_nonans,Tbin_deseason_nonans)
plt.plot(t_nonans,imfs[9,:],'k')
plt.show()



a= plt.figure()
axes = a.add_axes([0.1,0.1,0.8,0.8])
plt.scatter(t_nonans,Tbin_deseason_nonans)
plt.plot(t_nonans,imfs[9,:],'k')
plt.show()



# same but for annual means

tbin_ann,Tbin_ann = TF.bin_annually(1953,2020,tbin,Tbin_deseason)
tbin_ann = np.array(tbin_ann)
Tbin_ann = np.array(Tbin_ann)

check_nans = np.isfinite(Tbin_ann)
Tbin_ann_nonans = Tbin_ann[check_nans]
tbin_ann_nonans = tbin_ann[check_nans]

N = len(Tbin_ann_nonans)
t = np.arange(0, N) * dt + t0
t_nonans =t[check_nans]

# Extract imfs and residue
# In case of EMD
emd = EMD()
emd.emd(Tbin_ann_nonans)
imfs, res = emd.get_imfs_and_residue()

# In general:
#components = EEMD()(S)
#imfs, res = components[:-1], components[-1]

vis = Visualisation()
vis.plot_imfs(imfs=imfs, residue=res, t=tbin_ann_nonans, include_residue=True)
vis.plot_instant_freq(tbin_ann_nonans, imfs=imfs)
vis.show()


# re-create the signal from IMFs


test_multi_decadal = imfs[4,:]+imfs[3,:]

a= plt.figure()
axes = a.add_axes([0.1,0.1,0.8,0.8])
plt.scatter(tbin_ann_nonans,Tbin_ann_nonans)
plt.plot(tbin_ann_nonans,imfs[4,:],'k')
plt.show()




