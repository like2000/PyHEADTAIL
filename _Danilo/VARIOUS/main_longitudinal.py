
from __future__ import division
import numpy as np
from numpy import loadtxt
from scipy.constants import c, e, m_p
import time

from input_parameters.general_parameters import *
from input_parameters.rf_parameters import *
from trackers.longitudinal_tracker import *
from beams.beams import *
from beams.longitudinal_distributions import *
from longitudinal_plots.longitudinal_plots import *
from monitors.monitors import *
from beams.slices import *
from impedances.longitudinal_impedance import *

# Simulation parameters --------------------------------------------------------
# Bunch parameters
N_b = 1.e11           # Intensity
N_p = 500000          # Macro-particles
tau_0 = 2e-9       # Initial bunch length, 4 sigma [s]
sd = .5e-4           # Initial r.m.s. momentum spread
particle_type = 'proton'

# Machine and RF parameters

gamma_transition = 18.0 #55.759505  # Transition gamma
C = 2 * np.pi * 1100.009 #26658.883        # Machine circumference [m]
p_s = 450.e9         # Synchronous momentum [eV]


# Tracking details
N_t = 1          # Number of turns to track
dt_out = 20          # Time steps between output
dt_plt = 20          # Time steps between plots

# Derived parameters
#m_p *= c**2/e
E_s = np.sqrt(p_s**2 + m_p**2 * c**4 / e**2)  # Sychronous energy [eV]
gamma = E_s/(m_p*c**2/e)          # Relativistic gamma
beta = np.sqrt(1. - 1./gamma**2)  # Relativistic beta
T0 = C / beta / c                 # Turn period
sigma_theta = 2 * np.pi * tau_0 / T0      # R.m.s. theta
sigma_dE = sd * beta**2 * E_s           # R.m.s. dE


# Monitors
bunchmonitor = BunchMonitor('bunch', N_t+1, long_gaussian_fit = "Off")


n_rf_systems = 1                                      # Number of rf systems second section
harmonic_numbers = 4620#35640                               # Harmonic number second section
voltage_program = 6.e6
       

phi_offset = 0



momentum_compaction = 1./gamma_transition**2


general_params = GeneralParameters(N_t, C, momentum_compaction, p_s, 
                                   particle_type, number_of_sections = 1)

RF_sct_par = RFSectionParameters(general_params, 1, n_rf_systems, harmonic_numbers, 
                          voltage_program, phi_offset)

ring_RF_section = RingAndRFSection(RF_sct_par)


# Simulation setup -------------------------------------------------------------
print "Setting up the simulation..."
print ""



# Define Beam
my_beam = Beam(general_params, N_p, N_b)

longitudinal_gaussian_matched(general_params, RF_sct_par, my_beam, sigma_theta*2)

number_slices = 250
slice_beam = Slices(number_slices, cut_left = 0, cut_right = 10e-9, coord = "tau", mode = 'const_space_hist')

# slice_beam.track(my_beam)
# 
# bunchmonitor.dump(my_beam, slice_beam)

temp = loadtxt('new_HQ_table.dat', comments = '!')
R_shunt = temp[0, 2] * 10**6 
f_res = temp[0, 0] * 10**9
Q_factor = temp[0, 1]

resonator = Longitudinal_resonators(R_shunt, f_res, Q_factor)
ind_volt_from_wake = Induced_voltage_from_wake(slice_beam, 'off', [resonator], my_beam)
ind_volt_from_wake2 = Induced_voltage_from_impedance(slice_beam, 'off', [resonator], 0.5e6, my_beam)


# Accelerator map
map_ = [slice_beam] + [ind_volt_from_wake] + [ring_RF_section]
print "Map set"
print ""


# Tracking ---------------------------------------------------------------------
for i in range(N_t):
    
    # Track
    t0 = time.clock()
    for m in map_:
        m.track(my_beam)
    
    bunchmonitor.dump(my_beam, slice_beam)
    t1 = time.clock()
    print t1 - t0
    plot_induced_voltage_vs_bins_centers(i+1, general_params, ind_volt_from_wake, style = '-')
    # Plot
#     if (i % dt_plt) == 0:
#         plot_long_phase_space(my_beam, general_params, ring, -0.75, 0, -1.e-3, 1.e-3, xunit='m', yunit='1')

print "Done!"
print ""

bunchmonitor.h5file.close()