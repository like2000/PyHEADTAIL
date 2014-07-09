import cython_functions.stats as cp

import numpy as np

import time

x = np.random.randn(1000000)


t0 = time.clock()
cp.mean(x)

print time.clock() - t0

t0 = time.clock()
np.mean(x)

print time.clock() - t0

print 
