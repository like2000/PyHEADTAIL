import cython_functions.stats as cp

import numpy as np

import matplotlib.pyplot as plt

import time
 
l = np.linspace(1, 5000000, 100)

res = np.empty(len(l))


for i in range(len(l)):

    x = np.random.randn(l[i])
    
    t0 = time.clock()
    aa=cp.std(x)
    
    a = time.clock() - t0
    
    t0 = time.clock()
    bb=np.std(x)
    
    b = time.clock() - t0

    res[i] = a - b

print aa, bb
plt.plot(l, res)
plt.show()




