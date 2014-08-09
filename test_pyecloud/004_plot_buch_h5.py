import h5py
import pylab as pl 
import myloadmat_to_obj as mlm
import scipy as sp
import numpy as np

filename = 'bunch.h5'
#filename = 'bunch_highdens.h5'
filename = 'bunch_highdens_57_kick.h5'
#~ filename = 'bunch_noecloud_kick.h5'
#~ filename = 'bunch_Q26.h5'
filename = 'bunch_Q26_aftermerge.h5'
bunch = mlm.obj_from_dict(h5py.File(filename, 'r')['Bunch'])


pl.close('all')
pl.figure(1)
pl.plot(bunch.mean_x)
pl.plot(bunch.mean_y)





pl.figure(2)
pl.plot(bunch.epsn_x)
pl.plot(bunch.epsn_y)
pl.ylim(0, None)



spec_x = sp.fft(bunch.mean_x)
spec_y = sp.fft(bunch.mean_y)

pl.figure(3)
pl.plot(np.linspace(0, 1, len(spec_x)), np.abs(spec_x))
pl.plot(np.linspace(0, 1, len(spec_y)), np.abs(spec_y))

pl.figure(4)
pl.plot(bunch.sigma_z)
pl.plot(bunch.mean_z)
pl.show()


