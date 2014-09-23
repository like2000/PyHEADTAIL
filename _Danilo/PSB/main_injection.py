
from __future__ import division
import numpy as np
from numpy.fft import irfft
import math
from scipy.constants import c, e, m_p
import time, sys
import matplotlib as plt
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
import os

if os.path.exists('output2'):    
    os.system('del /s/q '+ os.getcwd() +'\\output2>null')
else:
    os.makedirs('output2')

# SIMULATION PARAMETERS -------------------------------------------------------

# Beam parameters
particle_type = 'proton'
n_particles = 1e13
n_macroparticles = 5e5
kin_beam_energy = 0.05e9 # [eV]
emittance = 1.0  # [eVs]


# Machine and RF parameters
radius = 25
gamma_transition = 4.05  # [1]
C = 2 * np.pi * radius  # [m]       
      
# Tracking details
n_turns = 10
n_turns_between_two_plots = 2

# Derived parameters
E_0 = m_p * c**2 / e    # [eV]
tot_beam_energy =  E_0 + kin_beam_energy # [eV]
sync_momentum = np.sqrt(tot_beam_energy**2 - E_0**2) # [eV / c]
momentum_compaction = 1 / gamma_transition**2 # [1]       

# Cavities parameters
n_rf_systems = 2                                     
harmonic_numbers_1 = 1    
harmonic_numbers_2 = 2                        
voltage_program_1 = 4.03e3   
voltage_program_2 = 3.81e3 
phi_offset_1 = 0
phi_offset_2 = np.pi 

# DEFINE RING------------------------------------------------------------------

general_params = GeneralParameters(n_turns, C, momentum_compaction, sync_momentum, 
                                   particle_type, number_of_sections = 1)

RF_sct_par = RFSectionParameters(general_params, n_rf_systems, [[harmonic_numbers_1], [harmonic_numbers_2]],
                          [[voltage_program_1], [voltage_program_2]], [[phi_offset_1], [phi_offset_2]])

longitudinal_tracker = RingAndRFSection(RF_sct_par)
full_tracker = FullRingAndRF([longitudinal_tracker])

# DEFINE BEAM------------------------------------------------------------------

my_beam = Beam(general_params, n_macroparticles, n_particles)


# DEFINE SLICES----------------------------------------------------------------

number_slices = 200
slice_beam = Slices(my_beam, number_slices, cut_left = - np.pi, 
                    cut_right = np.pi, mode = 'const_space_hist', cuts_coord = 'theta', slicing_coord = 'tau')


# # LOAD IMPEDANCE TABLES--------------------------------------------------------
 
var = str(kin_beam_energy / 1e9)
 
# ejection kicker
Ekicker = np.loadtxt('ps_booster_impedances/ejection_kicker/Ekicker_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))
 
Ekicker_table = InputTable(Ekicker[:,0].real, Ekicker[:,1].real, Ekicker[:,1].imag)
 
# ejection kicker cables
Ekicker_cables = np.loadtxt('ps_booster_impedances/ejection_kicker_cables/Ekicker_cables_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))
 
Ekicker_cables_table = InputTable(Ekicker_cables[:,0].real, Ekicker_cables[:,1].real, Ekicker_cables[:,1].imag)
 
# KSW kickers
KSW = np.loadtxt('ps_booster_impedances/KSW/KSW_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))
 
KSW_table = InputTable(KSW[:,0].real, KSW[:,1].real, KSW[:,1].imag)
 
# resistive wall
RW = np.loadtxt('ps_booster_impedances/resistive_wall/RW_' + var + 'GeV.txt'
        , skiprows = 1, dtype=complex, converters = dict(zip((0, 1), (lambda s: 
        complex(s.replace('i', 'j')), lambda s: complex(s.replace('i', 'j'))))))
 
RW_table = InputTable(RW[:,0].real, RW[:,1].real, RW[:,1].imag)
 
# indirect space charge
ISC = np.loadtxt('ps_booster_impedances/Indirect_space_charge/ISC_' + var + 'GeV.txt'
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
     -778.85, general_params.f_rev[0])

# cavity C02
cavity_02 = Resonators(350, 1.8*1e6, 2.8)

# cavity C04
cavity_04 = Resonators(430, 3.9*1e6, 6.8)

# cavity C16
cavity_16 = Resonators(0, 16*1e6, 30)

# INDUCED VOLTAGE FROM IMPEDANCE------------------------------------------------

imp_list = [F_C_table, steps_plus_space_charge, RW_table, KSW_table, Ekicker_cables_table, Ekicker_table]

ind_volt_freq = InducedVoltageFreq(slice_beam, imp_list, 1e5)

total_induced_voltage = TotalInducedVoltage(slice_beam, [ind_volt_freq])

# SAVING DATA----------------------------------------------------------------------

save_bunch_length = np.zeros(n_turns + 1) 
save_bunch_position = np.zeros(n_turns + 1)
save_perc_number_particles_lost = np.zeros(n_turns + 1)


# BEAM GENERATION

matched_from_distribution_density(my_beam, full_tracker, 
            {'type':'parabolic_amplitude', 'emittance':emittance, 'density_variable':'density_from_H'}, 
            main_harmonic_option = 'lowest_freq', TotalInducedVoltage = 
            total_induced_voltage, n_iterations_input = 3)

my_beam.theta -= np.pi / 180
my_beam.losses_longitudinal_cut(-np.pi, np.pi)
itemindex = np.where(my_beam.id != 0)[0]
save_perc_number_particles_lost[0] = 100 - len(itemindex) * 100 / n_macroparticles
save_bunch_length[0] = 4 * np.std(my_beam.theta[itemindex]) 
save_bunch_position[0] = np.mean(my_beam.theta[itemindex])
           

# ACCELERATION MAP-------------------------------------------------------------

map_ = [slice_beam] + [total_induced_voltage] + [full_tracker]


# TRACKING + PLOTS-------------------------------------------------------------

for i in range(n_turns):
    t0 = time.clock()
    my_beam.losses_longitudinal_cut(-np.pi, np.pi)
    itemindex = np.where(my_beam.id != 0)[0]
    save_perc_number_particles_lost[i+1] = 100 - len(itemindex) * 100 / n_macroparticles
    print '%d turn'%i
    for m in map_:
        m.track(my_beam)
    save_bunch_length[i+1] = 4 * np.std(my_beam.theta[itemindex]) 
    save_bunch_position[i+1] = np.mean(my_beam.theta[itemindex])
    if (i % n_turns_between_two_plots) == 0:
        plot_long_phase_space(my_beam, general_params, RF_sct_par, 
        - np.pi, np.pi, 
        - 0.1e1, 0.1e1, sampling = 10, separatrix_plot = False, dirname = 'fig2')
        plt.figure() 
        ax = plt.axes() 
        plt.plot((ind_volt_freq.slices.bins_centers - np.min(ind_volt_freq.slices.bins_centers))*1e5, irfft(ind_volt_freq.total_impedance)[0:ind_volt_freq.slices.n_slices] * ind_volt_freq.slices.beam_spectrum_freq[1] * 2*(len(ind_volt_freq.slices.beam_spectrum)-1), '-')
        ax.set_xlabel("Time [1e5]") 
        ax.set_ylabel ("Wake")
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(10)
            tick.label1.set_fontweight('bold')
        plt.savefig('output2/%d_turn'%i+'wake.png')
        plt.clf()
    print time.clock() - t0

plot_long_phase_space(my_beam, general_params, RF_sct_par, 
        - np.pi, np.pi, 
        - 0.1e1, 0.1e1, sampling = 10, separatrix_plot = False, dirname = 'fig2')

plt.figure() 
ax = plt.axes() 
plt.plot((ind_volt_freq.slices.bins_centers - np.min(ind_volt_freq.slices.bins_centers))*1e5, irfft(ind_volt_freq.total_impedance)[0:ind_volt_freq.slices.n_slices] * ind_volt_freq.slices.beam_spectrum_freq[1] * 2*(len(ind_volt_freq.slices.beam_spectrum)-1), '-')
ax.set_xlabel("Time [1e5]") 
ax.set_ylabel ("Wake")
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(10)
    tick.label1.set_fontweight('bold')
plt.savefig('output2/%d_turn'%i+'wake.png')
plt.clf()

np.savez_compressed('output2/phase_space', my_beam.theta, my_beam.dE, n_turns)
np.savez_compressed('output2/bunch_length', save_bunch_length)
np.savez_compressed('output2/bunch_position', save_bunch_position)
np.savez_compressed('output2/perc_particles_lost', save_perc_number_particles_lost)

t = np.arange(0, n_turns+1)
plt.figure(1, figsize=(8,6)) 
ax = plt.axes([0.12, 0.1, 0.82, 0.8]) 
ax.plot(t, save_bunch_length, '-') 
ax.set_xlabel(r"No. turns [T$_0$]") 
ax.set_ylabel (r"Bunch length [rad]")
plt.savefig('output2/length_evolution.png')
plt.clf()

plt.figure(1, figsize=(8,6)) 
ax = plt.axes([0.12, 0.1, 0.82, 0.8]) 
ax.plot(t, save_bunch_position, '-') 
ax.set_xlabel(r"No. turns [T$_0$]") 
ax.set_ylabel (r"Bunch position [rad]")
plt.savefig('output2/position_evolution.png')
plt.clf()

plt.figure(1, figsize=(8,6)) 
ax = plt.axes([0.12, 0.1, 0.82, 0.8]) 
ax.plot(t, save_perc_number_particles_lost, '-') 
ax.set_xlabel(r"No. turns [T$_0$]") 
ax.set_ylabel (r"Percentage_particles_lost")
plt.savefig('output2/Percentage_particles_lost.png')
plt.clf()
print "Done!"


