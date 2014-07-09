
from __future__ import division
import numpy as np
from numpy import loadtxt
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


# SIMULATION PARAMETERS

# Bunch parameters
N_b = 1.e11          
N_p = 1.e5           
tau_0 = 180e-9       
sigma_delta = .5e-4           
particle_type = 'proton'

# Machine and RF parameters

gamma_transition = 4.4  
C = 2 * np.pi * 25.0         
sync_momentum = 2.12e9      

# Tracking details
N_t = 3          
dt_plt = 1          

# Derived parameters

E_s = np.sqrt(sync_momentum**2 + m_p**2 * c**4 / e**2)  #[eV]
gamma = E_s/(m_p*c**2/e)          

beta = np.sqrt(1 - 1/gamma**2)  
T0 = C / beta / c                 

sigma_theta = (np.pi * tau_0) / (2 * T0)      
sigma_dE = sigma_delta * beta**2 * E_s           

n_rf_systems = 1                                     
harmonic_numbers = 1                         
voltage_program = 8.e3
phi_offset = 0

section_params = RFSectionParameters(N_t, n_rf_systems, C, harmonic_numbers, voltage_program, phi_offset, sync_momentum)

momentum_compaction = 1./gamma_transition**2

general_params = General_parameters(particle_type, N_t, C, momentum_compaction, section_params.momentum_program)

# MONITOR

bunchmonitor = BunchMonitor('beam', N_t+1, statistics = "Longitudinal")

# DEFINE RING

ring = RingAndRFSection(general_params, section_params)

# DEFINE BEAM

my_beam = Beam(general_params, N_p, N_b)

longitudinal_bigaussian(general_params, ring, my_beam, sigma_theta, sigma_dE)

# DEFINE SLICES

number_slices = 2000
slice_beam = Slices(number_slices, cut_left = - 5.72984173562e-07 / 2, cut_right = 5.72984173562e-07 / 2, unit = 
                    "tau", mode = 'const_space_hist')

# LOAD IMPEDANCE TABLES

Ekicker = loadtxt('ps_booster_impedances/ejection kicker/Ekicker_1.4GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

Ekicker_table = Longitudinal_table(Ekicker[:,0].real, Ekicker[:,1].real, Ekicker[:,1].imag)


Ekicker_cables = loadtxt('ps_booster_impedances/ejection kicker cables/Ekicker_cables_1.4GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

Ekicker_cables_table = Longitudinal_table(Ekicker_cables[:,0].real, Ekicker_cables[:,1].real, Ekicker_cables[:,1].imag)


ISC = loadtxt('ps_booster_impedances/Indirect space charge/ISC_1.4GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

ISC_table = Longitudinal_table(ISC[:,0].real, ISC[:,1].real, ISC[:,1].imag)


KSW = loadtxt('ps_booster_impedances/KSW/KSW_1.4GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

KSW_table = Longitudinal_table(KSW[:,0].real, KSW[:,1].real, KSW[:,1].imag)


RW = loadtxt('ps_booster_impedances/resistive wall/RW_1.4GeV.txt'
        , dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

RW_table = Longitudinal_table(RW[:,0].real, RW[:,1].real, RW[:,1].imag)

steps = Longitudinal_inductive_impedance(0.061, general_params)

sum_impedance = [Ekicker_table] + [Ekicker_cables_table] + [KSW_table] + [RW_table] + [steps]

ind_volt_from_imp = Induced_voltage_from_impedance(slice_beam, "off", sum_impedance, 2e6, my_beam)


# ACCELERATION MAP

map_ = [slice_beam] + [ind_volt_from_imp] + [ring]

# TRACKING

for i in range(N_t):
    
    t0 = time.clock()
    for m in map_:
        m.track(my_beam)
    general_params.counter[0] += 1
    bunchmonitor.dump(my_beam, slice_beam)
    t1 = time.clock()
    print t1 - t0
    
    # Plot
    if (i % dt_plt) == 0:
        plot_long_phase_space(my_beam, general_params, ring, - 5.72984173562e-07 / 2 * 1e9, 5.72984173562e-07 / 2 * 1e9, -0.5, 0.5, xunit = 'ns')
        plot_impedance_vs_frequency(general_params, ind_volt_from_imp)

print "Done!"

bunchmonitor.h5file.close()