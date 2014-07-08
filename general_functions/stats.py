import numpy as np

def emittance(u, v):
    
    u2 = np.sum((u - np.mean(u)) ** 2)
    v2 = np.sum((v - np.mean(v)) ** 2)
    uv = np.sum((u - np.mean(u)) * (v - np.mean(v)))
    
    return np.sqrt(u2 * v2 - uv ** 2)