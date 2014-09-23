import numpy as np
from scipy.constants import c

beta_program_table = np.loadtxt('ps_booster_impedances/beta_program_0.05_1.4.txt')
time_table = beta_program_table[:, 0] * 1e-3
beta_table = beta_program_table[:, 1] 

time_program = []
beta_program = []

t = 0.276 
beta = 0.31367

time_program = np.append(time_program, t)
beta_program = np.append(beta_program, beta)

while t <= 0.810:
    print t
    t = t + (2 * np.pi * 25)/(beta * c)
    beta = np.interp(t, time_table, beta_table)
    time_program = np.append(time_program, t)
    beta_program = np.append(beta_program, beta)

time_program = np.array([time_program]).T
beta_program = np.array([beta_program]).T
print time_program
print beta_program
np.savetxt('real_beta_program.out', np.hstack([time_program, beta_program]), fmt=['%0.8f','%0.8f'])    





        