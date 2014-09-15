from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


npzfile = np.load('perc_particles_lost.npz')
x = npzfile['arr_0']
print x