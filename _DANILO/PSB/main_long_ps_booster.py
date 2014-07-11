
from __future__ import division
import numpy as np
from numpy import loadtxt
import math
from scipy.constants import c, e, m_p
import time, sys
import matplotlib.pyplot as plt

from input_parameters.general_parameters import *
from input_parameters.rf_parameters import *
from trackers.longitudinal_tracker import *
from beams.beams import *
from beams.longitudinal_distributions import *
from longitudinal_plots.longitudinal_plots import *
from monitors.monitors import *
from beams.slices import *
from impedances.longitudinal_impedance import *


# SIMULATION PARAMETERS -------------------------------------------------------

# Beam parameters
particle_type = 'proton'
N_b = 1.e11          
N_p = 5.e5           
sigma_tau = 180e-9 / 4 # [s]     
sigma_delta = .5e-4 # [1]          
kin_beam_energy = 1.4e9 # [eV]

# Machine and RF parameters
radius = 25
gamma_transition = 4.4  # [1]
C = 2 * np.pi * radius  # [m]       
      
# Tracking details
n_turns = 100          
n_turns_between_two_plots = 1          

# Derived parameters
E_0 = m_p * c**2 / e
tot_beam_energy =  E_0 + kin_beam_energy # [eV]
sync_momentum = np.sqrt(tot_beam_energy**2 - E_0**2) # [eV / c]
gamma = tot_beam_energy / E_0  # [1]        
beta = np.sqrt(1 - 1 / gamma**2)  # [1]
sigma_theta = beta * c / radius * sigma_tau # [rad]     
sigma_dE = beta**2 * tot_beam_energy * sigma_delta # [eV]
momentum_compaction = 1 / gamma_transition**2 # [1]       

# Cavities parameters
n_rf_systems = 1                                     
harmonic_numbers = 1                         
voltage_program = 8.e3
phi_offset = 0

# MONITOR----------------------------------------------------------------------

bunchmonitor = BunchMonitor('beam', n_turns+1, statistics = "Longitudinal")

# DEFINE RING------------------------------------------------------------------

section_params = RFSectionParameters(n_turns, n_rf_systems, C, harmonic_numbers, voltage_program, phi_offset, sync_momentum)
general_params = General_parameters(particle_type, n_turns, C, momentum_compaction, section_params.momentum_program)
ring = RingAndRFSection(general_params, section_params)

# DEFINE BEAM------------------------------------------------------------------

my_beam = Beam(general_params, N_p, N_b)

longitudinal_bigaussian(general_params, ring, my_beam, sigma_theta, sigma_dE)

# DEFINE SLICES----------------------------------------------------------------

number_slices = 2000
slice_beam = Slices(number_slices, cut_left = - 5.72984173562e-07 / 2, cut_right = 5.72984173562e-07 / 2, unit = 
                    "tau", mode = 'const_space_hist')


# LOAD IMPEDANCE TABLES--------------------------------------------------------

var = str(kin_beam_energy / 1e9)

# ejection kicker
Ekicker = loadtxt('ps_booster_impedances/ejection kicker/Ekicker_' + var + 'GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

Ekicker_table = Longitudinal_table(Ekicker[:,0].real, Ekicker[:,1].real, Ekicker[:,1].imag)

# ejection kicker cables
Ekicker_cables = loadtxt('ps_booster_impedances/ejection kicker cables/Ekicker_cables_' + var + 'GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

Ekicker_cables_table = Longitudinal_table(Ekicker_cables[:,0].real, Ekicker_cables[:,1].real, Ekicker_cables[:,1].imag)

# KSW magnets
KSW = loadtxt('ps_booster_impedances/KSW/KSW_' + var + 'GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

KSW_table = Longitudinal_table(KSW[:,0].real, KSW[:,1].real, KSW[:,1].imag)

# resistive wall
RW = loadtxt('ps_booster_impedances/resistive wall/RW_' + var + 'GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

RW_table = Longitudinal_table(RW[:,0].real, RW[:,1].real, RW[:,1].imag)

# indirect space charge
ISC = loadtxt('ps_booster_impedances/Indirect space charge/ISC_' + var + 'GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

ISC_table = Longitudinal_table(ISC[:,0].real, ISC[:,1].real, ISC[:,1].imag)

# steps
steps = Longitudinal_inductive_impedance(0.061, general_params)

# Finemet cavity

F_C = loadtxt('ps_booster_impedances/Finemet_cavity/Finemet.txt', dtype = float, skiprows = 1)

F_C[:, 3], F_C[:, 5], F_C[:, 7] = np.pi * F_C[:, 3] / 180, np.pi * F_C[:, 5] / 180, np.pi * F_C[:, 7] / 180
F_C[:, 2], F_C[:, 4], F_C[:, 6] = 10 ** (F_C[:, 2]/20), 10 ** (F_C[:, 4]/20), 10 ** (F_C[:, 6]/20)

option = "closed loop"

if option == "open loop":
    Re_Z = F_C[:, 4] * np.cos(F_C[:, 3])
    Im_Z = F_C[:, 4] * np.sin(F_C[:, 3])
    F_C_table = Longitudinal_table(F_C[:, 0], 13 * Re_Z, 13 * Im_Z)
elif option == "closed loop":
    Re_Z = F_C[:, 2] * np.cos(F_C[:, 5])
    Im_Z = F_C[:, 2] * np.sin(F_C[:, 5])
    F_C_table = Longitudinal_table(F_C[:, 0], 13 * Re_Z, 13 * Im_Z)
elif option == "shorted":
    Re_Z = F_C[:, 6] * np.cos(F_C[:, 7])
    Im_Z = F_C[:, 6] * np.sin(F_C[:, 7])
    F_C_table = Longitudinal_table(F_C[:, 0], 13 * Re_Z, 13 * Im_Z)
else:
    pass


# INDUCED VOLTAGE FROM IMPEDANCE------------------------------------------------

sum_impedance = [Ekicker_table] + [Ekicker_cables_table] + [KSW_table] + [RW_table] + [steps] + [F_C_table]

ind_volt_from_imp = Induced_voltage_from_impedance(slice_beam, "off", sum_impedance, 2e6, my_beam)


# ACCELERATION MAP-------------------------------------------------------------

map_ = [slice_beam] + [ind_volt_from_imp] + [ring]


# TRACKING + PLOTS-------------------------------------------------------------

plot_impedance_vs_frequency(general_params, ind_volt_from_imp, option = "single")

for i in range(n_turns):
    
    print i
    t0 = time.clock()
    for m in map_:
        m.track(my_beam)
    general_params.counter[0] += 1
    bunchmonitor.dump(my_beam, slice_beam)
    t1 = time.clock()
    print t1 - t0
    # Plots that change from turn to turn
    if (i % n_turns_between_two_plots) == 0:
        plot_long_phase_space(my_beam, general_params, ring, 
          - 5.72984173562e-07 / 2 * 1e9, 5.72984173562e-07 / 2 * 1e9, 
          - my_beam.sigma_dE * 4 * 1e-6, my_beam.sigma_dE * 4 * 1e-6, xunit = 'ns')
        
plot_bunch_length_evol(my_beam, 'beam', general_params)

print "Done!"

bunchmonitor.h5file.close()