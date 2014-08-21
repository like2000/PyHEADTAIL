
from __future__ import division
import numpy as np
import math
from scipy.constants import c, e, m_p
import time, sys
import matplotlib.pyplot as plt

from input_parameters.general_parameters import *
from input_parameters.rf_parameters import *
from trackers.longitudinal_tracker import *
from beams.beams import *
from beams.longitudinal_distributions import *
from longitudinal_plots.plot_beams import *
from longitudinal_plots.plot_impedance import *
from longitudinal_plots.plot_slices import *
from monitors.monitors import *
from beams.slices import *
from impedances.longitudinal_impedance import *


# SIMULATION PARAMETERS -------------------------------------------------------

# Beam parameters
particle_type = 'proton'
n_particles = int(1.6e12)          
n_macroparticles = int(1e6)
kin_beam_energy = 1.4e9 # [eV]
emittance = 1.3 # [eVs]

# Machine and RF parameters
radius = 25
gamma_transition = 4.05  # [1]
C = 2 * np.pi * radius  # [m]       
      
# Tracking details
n_turns = 10000       
n_turns_between_two_plots = 100          

# Derived parameters
E_0 = m_p * c**2 / e    # [eV]
tot_beam_energy =  E_0 + kin_beam_energy # [eV]
sync_momentum = np.sqrt(tot_beam_energy**2 - E_0**2) # [eV / c]
momentum_compaction = 1 / gamma_transition**2 # [1]       

# Cavities parameters
n_rf_systems = 2                                     
harmonic_numbers_1 = 1    
harmonic_numbers_2 = 2                        
voltage_program_1 = 8e3
voltage_program_2 = 1.3e3
phi_offset_1 = 0
phi_offset_2 = np.pi

# DEFINE RING------------------------------------------------------------------

general_params = GeneralParameters(n_turns, C, momentum_compaction, sync_momentum, 
                                   particle_type, number_of_sections = 1)
print len([[harmonic_numbers_1], [harmonic_numbers_2]])
RF_sct_par = RFSectionParameters(general_params, n_rf_systems, [[harmonic_numbers_1], [harmonic_numbers_2]],
                          [[voltage_program_1], [voltage_program_2]], [[phi_offset_1], [phi_offset_2]])

longitudinal_tracker = RingAndRFSection(RF_sct_par)
full_tracker = FullRingAndRF([longitudinal_tracker])

# DEFINE BEAM------------------------------------------------------------------

my_beam = Beam(general_params, n_macroparticles, n_particles)


# DEFINE SLICES----------------------------------------------------------------

number_slices = 200
slice_beam = Slices(my_beam, number_slices, cut_left = - 5.72984173562e-07 / 2, 
                    cut_right = 5.72984173562e-07 / 2, mode = 'const_space_hist')

# MONITOR----------------------------------------------------------------------

bunchmonitor = BunchMonitor('6', n_turns+1, statistics = "Longitudinal")


# LOAD IMPEDANCE TABLES--------------------------------------------------------

var = str(kin_beam_energy / 1e9)

# ejection kicker
Ekicker = np.loadtxt('ps_booster_impedances/ejection kicker/Ekicker_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

Ekicker_table = InputTable(Ekicker[:,0].real, Ekicker[:,1].real, Ekicker[:,1].imag)

# ejection kicker cables
Ekicker_cables = np.loadtxt('ps_booster_impedances/ejection kicker cables/Ekicker_cables_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

Ekicker_cables_table = InputTable(Ekicker_cables[:,0].real, Ekicker_cables[:,1].real, Ekicker_cables[:,1].imag)

# KSW kickers
KSW = np.loadtxt('ps_booster_impedances/KSW/KSW_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

KSW_table = InputTable(KSW[:,0].real, KSW[:,1].real, KSW[:,1].imag)

# resistive wall
RW = np.loadtxt('ps_booster_impedances/resistive wall/RW_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

RW_table = InputTable(RW[:,0].real, RW[:,1].real, RW[:,1].imag)

# indirect space charge
ISC = np.loadtxt('ps_booster_impedances/Indirect space charge/ISC_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))

ISC_table = InputTable(ISC[:,0].real, ISC[:,1].real, ISC[:,1].imag)


# Finemet cavity

F_C = np.loadtxt('ps_booster_impedances/Finemet_cavity/Finemet.txt', dtype = float, skiprows = 1)

F_C[:, 3], F_C[:, 5], F_C[:, 7] = np.pi * F_C[:, 3] / 180, np.pi * F_C[:, 5] / 180, np.pi * F_C[:, 7] / 180

option = "closed loop"

if option == "open loop":
    Re_Z = F_C[:, 4] * np.cos(F_C[:, 3])
    Im_Z = F_C[:, 4] * np.sin(F_C[:, 3])
    F_C_table = InputTable(F_C[:, 0], 13 * Re_Z, 13 * Im_Z)
elif option == "closed loop":
    Re_Z = F_C[:, 2] * np.cos(F_C[:, 5])
    Im_Z = F_C[:, 2] * np.sin(F_C[:, 5])
    F_C_table = InputTable(F_C[:, 0], 13 * Re_Z, 13 * Im_Z)
elif option == "shorted":
    Re_Z = F_C[:, 6] * np.cos(F_C[:, 7])
    Im_Z = F_C[:, 6] * np.sin(F_C[:, 7])
    F_C_table = InputTable(F_C[:, 0], 13 * Re_Z, 13 * Im_Z)
else:
    pass


# inductive impedance from steps plus space charge

steps_plus_space_charge = InductiveImpedance(slice_beam, 
     34.6669349520904 / 10e9 * general_params.f_rev[0] 
     -376.730313462 / (general_params.beta_r[0,0] * general_params.gamma_r[0,0]**2), general_params.f_rev[0])

# cavity C02
cavity_02 = Resonators(350, 1.8*1e6, 2.8)

# cavity C04
cavity_04 = Resonators(430, 3.9*1e6, 6.8)

# cavity C16
cavity_16 = Resonators(0, 16*1e6, 30)

# INDUCED VOLTAGE FROM IMPEDANCE------------------------------------------------

imp_list = [steps_plus_space_charge, RW_table, KSW_table, Ekicker_cables_table, Ekicker_table]
ind_volt_freq = InducedVoltageFreq(slice_beam, imp_list, 1e5)

total_induced_voltage = TotalInducedVoltage(slice_beam, [ind_volt_freq])

# BEAM GENERATION

matched_from_distribution_density(my_beam, full_tracker, {'type':'parabolic', 
            'parameters':[emittance, 1.5], 'density_variable':'density_from_H'}, 
            main_harmonic_option = 'lowest_freq', TotalInducedVoltage = 
            total_induced_voltage, n_iterations_input = 100)

slice_beam.track(my_beam)
bunchmonitor.track(my_beam)


# ACCELERATION MAP-------------------------------------------------------------

map_ = [total_induced_voltage] + [full_tracker] + [slice_beam] + [bunchmonitor] 


# TRACKING + PLOTS-------------------------------------------------------------

for i in range(n_turns):
    
    print '%d turn'%i
    t0 = time.clock()
    for m in map_:
        m.track(my_beam)
    t1 = time.clock()
    print t1 - t0
    # Plots
    
    if ((i+1) % n_turns_between_two_plots) == 0:
#         
        plot_long_phase_space(my_beam, general_params, RF_sct_par, 
        - 5.72984173562e-07 / 2 * 1e9, 5.72984173562e-07 / 2 * 1e9, 
        - 1e1, 1e1, xunit = 'ns', sampling = 10, separatrix_plot = True, dirname = '6_figNoFC')
          
        plot_induced_voltage_vs_bins_centers(i+1, general_params, total_induced_voltage, 
                                             style = '-', dirname = '6_figNoFC')
           
        plot_beam_profile(i+1, general_params, slice_beam, dirname = '6_figNoFC')
    
    
# plot_impedance_vs_frequency(n_turns, general_params, ind_volt_freq, 
#            option1 = "single", style = '-', option3 = "freq_table", option2 = "spectrum", dirname = 'figCL',
#            cut_right = 0.1e8)
plot_bunch_length_evol(my_beam, '6', general_params, n_turns, dirname = '6_figNoFC')
plot_position_evol(n_turns, my_beam, '6', general_params, unit = None, dirname = '6_figNoFC')


print "Done!"

bunchmonitor.h5file.close()
