'''
SPS test case
'''

# General imports
from __future__ import division
import numpy as np
import time
import matplotlib.pyplot as plt
import os
os.chdir('/users/alasheen/documents/work/eclipse/git/PyHEADTAIL_alasheen/_Alex/')

# PyHEADTAIL imports
from input_parameters.general_parameters import GeneralParameters 
from input_parameters.rf_parameters import RFSectionParameters
from beams.beams import Beam
from beams.longitudinal_distributions import matched_from_distribution_density, \
                                             matched_from_line_density   
from trackers.longitudinal_tracker import RingAndRFSection, FullRingAndRF
from beams.slices import Slices
from impedances.longitudinal_impedance import TravelingWaveCavity, Resonators, \
                                              InductiveImpedance, InputTable, \
                                              InducedVoltageTime, InducedVoltageFreq, \
                                              TotalInducedVoltage
# from beams.plot_beams import plot_long_phase_space
# 
# # Other imports
# from scipy.optimize import curve_fit, brentq
# from scipy.interpolate import interp1d
from scipy.constants import c


# Simulation parameters -------------------------------------------------------
# Simulation parameters
n_turns = 1500          # Number of turns to track
plot_step = 100           # Time steps between plots
output_step = 100

# General parameters
particle_type = 'proton'
circumference = 6911.56                         # Machine circumference [m]
gamma_transition = 1/np.sqrt(0.00192)           # Transition gamma
# gamma_transition = 18.                      # Transition gamma
momentum_compaction = 1./gamma_transition**2    # Momentum compaction array
sync_momentum = 25.92e9                         # Synchronous momentum program [eV/c]
# sync_momentum = 450e9                         # Synchronous momentum program [eV/c]

# RF parameters
n_rf_systems = 2                # Number of rf systems second section
harmonic_numbers_1 = 4620       # Harmonic number second section
# voltage_program_1 = 7.e6       # RF voltage [V] second section
voltage_program_1 = 0.9e6       # RF voltage [V] second section
phi_offset_1 = 0                # Phase offset
harmonic_numbers_2 = 4620*4     # Harmonic number second section
# voltage_program_2 = 0.65e6        # RF voltage [V] second section
voltage_program_2 = 0.e6        # RF voltage [V] second section
phi_offset_2 = np.pi            # Phase offset

# Beam parameters
# intensity = 2.6e11           # Intensity
intensity = 8.e10           # Intensity
n_macroparticles = 5e5   # Macro-particles
# tau_0 = 2.0                 # Initial bunch length, 4 sigma [ns]
# beam_loaded = np.loadtxt('test_input_files/distributions/distribution_parabolic_ESME_0.12.crd')
# n_macroparticles = len(beam_loaded[:,0])   # Macro-particles

# Slicing parameters
n_slices = 2**8
cut_left = 0.
cut_right = 2*np.pi / harmonic_numbers_1
frequency_step = 94.e3 # [Hz]

# Impedance parameters
TWC_import = np.loadtxt('test_input_files/SPS_impedance/TWC.dat', comments='!')
HOM_cavities_import = np.loadtxt('test_input_files/SPS_impedance/HOM_cavities.dat', comments='!')
Beam_scrapers_import = np.loadtxt('test_input_files/SPS_impedance/Beam_scrapers.dat', comments='!')
BPHs_import = np.loadtxt('test_input_files/SPS_impedance/BPHs.dat', comments='!')
BPVs_import = np.loadtxt('test_input_files/SPS_impedance/BPVs.dat', comments='!')
Flanges_import = np.loadtxt('test_input_files/SPS_impedance/Flanges.dat', comments='!')
Pumping_ports_import = np.loadtxt('test_input_files/SPS_impedance/Pumping_ports.dat', comments='!')
QDQDClosedFlange_import = np.loadtxt('test_input_files/SPS_impedance/QD-QD-Closed-Flange.dat', comments='!')
QDQDEnamFlange_import = np.loadtxt('test_input_files/SPS_impedance/QD-QD-Enam-Flange.dat', comments='!')
QFQFClosedFlange_import = np.loadtxt('test_input_files/SPS_impedance/QF-QF-Closed-Flange.dat', comments='!')
Y_chambers_import = np.loadtxt('test_input_files/SPS_impedance/Y_chambers.dat', comments='!')

Kickers_import = np.loadtxt('test_input_files/SPS_impedance/Kickers_impedance.dat', comments='!')
RW_import = np.loadtxt('test_input_files/SPS_impedance/Resistive_wall.dat', comments='!')

Space_charge_Z_over_n = -1.5 # Ohms
Steps_Z_over_n = 0.5 # Ohms

    
# Simulation setup ------------------------------------------------------------
# General parameters                  
general_params = GeneralParameters(n_turns, circumference, momentum_compaction, 
                                   sync_momentum, particle_type)
# general_params = GeneralParameters(n_turns, [circumference/2, circumference/2], 
#                                    [[momentum_compaction], [momentum_compaction]], 
#                                    [[sync_momentum], [sync_momentum]], particle_type,
#                                    number_of_sections = 2)

# RF parameters
# rf_params = RFSectionParameters(general_params, 1, harmonic_numbers_1, 
#                                       voltage_program_1, phi_offset_1)
rf_params = RFSectionParameters(general_params, n_rf_systems, [harmonic_numbers_1, harmonic_numbers_2], 
                                [voltage_program_1, voltage_program_2], [phi_offset_1, phi_offset_2])
# rf_params_1 = RFSectionParameters(general_params, n_rf_systems, [harmonic_numbers_1, harmonic_numbers_2], 
#                                 [voltage_program_1/2, voltage_program_2/2], [phi_offset_1, phi_offset_2], section_index = 1)
# rf_params_2 = RFSectionParameters(general_params, n_rf_systems, [harmonic_numbers_1, harmonic_numbers_2], 
#                                 [voltage_program_1/2, voltage_program_2/2], [phi_offset_1, phi_offset_2], section_index = 2)

# RF tracker
longitudinal_tracker = RingAndRFSection(rf_params)
full_tracker = FullRingAndRF([longitudinal_tracker])
# longitudinal_tracker_1 = RingAndRFSection(rf_params_1)
# longitudinal_tracker_2 = RingAndRFSection(rf_params_2)
# full_tracker = FullRingAndRF([longitudinal_tracker_1, longitudinal_tracker_2])

# Beam
SPS_beam = Beam(general_params, n_macroparticles, intensity)
# general_params_dummy = GeneralParameters(n_turns, circumference, momentum_compaction, 
#                                    sync_momentum, particle_type)
# rf_params_dummy = RFSectionParameters(general_params, 1, harmonic_numbers_1, 
#                                       voltage_program_1, phi_offset_1)
# longitudinal_gaussian_matched(general_params, rf_params_dummy, SPS_beam, tau_0, 
#                               unit='ns')
# SPS_beam.z = beam_loaded[:,0]
# SPS_beam.delta = beam_loaded[:,1]
# SPS_beam.theta = SPS_beam.theta + np.pi/harmonic_numbers_1

 
# Slicing
slicing = Slices(SPS_beam, n_slices, cut_left = cut_left, cut_right = cut_right, 
                 cuts_coord = 'theta', slicing_coord = 'tau', mode = 'const_space_hist', fit_option = 'gaussian') #
  
   
# Impedance sources
SPS_TWC = TravelingWaveCavity(TWC_import[:,1]*1e3, TWC_import[:,0]*1e9, TWC_import[:,2]*1e-6)
SPS_HOM = Resonators(HOM_cavities_import[1]*1e3, HOM_cavities_import[0]*1e9, HOM_cavities_import[2])
SPS_BS = Resonators(Beam_scrapers_import[:,1]*1e3, Beam_scrapers_import[:,0]*1e9, Beam_scrapers_import[:,2])
SPS_BPH = Resonators(BPHs_import[:,1]*1e3, BPHs_import[:,0]*1e9, BPHs_import[:,2])
SPS_BPV = Resonators(BPVs_import[:,1]*1e3, BPVs_import[:,0]*1e9, BPVs_import[:,2])
SPS_Flanges = Resonators(Flanges_import[:,1]*1e3, Flanges_import[:,0]*1e9, Flanges_import[:,2])
SPS_PP = Resonators(Pumping_ports_import[:,1]*1e3, Pumping_ports_import[:,0]*1e9, Pumping_ports_import[:,2])
SPS_Flanges2 = Resonators(QDQDClosedFlange_import[:,1]*1e3, QDQDClosedFlange_import[:,0]*1e9, QDQDClosedFlange_import[:,2])
SPS_Flanges3 = Resonators(QDQDEnamFlange_import[:,1]*1e3, QDQDEnamFlange_import[:,0]*1e9, QDQDEnamFlange_import[:,2])
SPS_Flanges4 = Resonators(QFQFClosedFlange_import[:,1]*1e3, QFQFClosedFlange_import[:,0]*1e9, QFQFClosedFlange_import[:,2])
SPS_YC = Resonators(Y_chambers_import[:,1]*1e3, Y_chambers_import[:,0]*1e9, Y_chambers_import[:,2])
    
SPS_Kickers = InputTable(Kickers_import[:,0]*1e9, Kickers_import[:,1], Kickers_import[:,2])
SPS_RW = InputTable(RW_import[:,0], RW_import[:,1]*(RW_import[:,0]/general_params.f_rev[0]), RW_import[:,2]*(RW_import[:,0]/general_params.f_rev[0]))
    
Wake_list = [SPS_TWC, SPS_HOM, SPS_BS, SPS_BPH, SPS_BPV, SPS_Flanges, SPS_PP, 
             SPS_Flanges2, SPS_Flanges3, SPS_Flanges4, SPS_YC]
SPS_intensity_time = InducedVoltageTime(slicing, Wake_list)
    
Impedance_list = [SPS_Kickers, SPS_RW]
SPS_intensity_freq = InducedVoltageFreq(slicing, Impedance_list, frequency_step)
    
SPS_inductive = InductiveImpedance(slicing, Space_charge_Z_over_n + Steps_Z_over_n, general_params.f_rev[0], deriv_mode = 'filter1d')
    
SPS_longitudinal_intensity = TotalInducedVoltage(slicing, [SPS_intensity_time, SPS_intensity_freq]) # SPS_inductive


# Beam generation
# emittance = 0.4
# coefficient = 2.
# matched_from_distribution_density(SPS_beam, full_tracker, {'type':'parabolic', 'parameters':[emittance, coefficient], 'density_variable':'density_from_H'}, main_harmonic_option = 'lowest_freq', TotalInducedVoltage = SPS_longitudinal_intensity, n_iterations_input = 50)
bunch_length = 2.e-9 / (SPS_beam.ring_radius / (SPS_beam.beta_r * c))
coefficient = 5/2
matched_from_line_density(SPS_beam, full_tracker, {'type':'parabolic', 'parameters':[bunch_length,coefficient]}, main_harmonic_option = 'lowest_freq', TotalInducedVoltage = SPS_longitudinal_intensity)


# # Initial kick
# SPS_longitudinal_intensity.induced_voltage_sum(SPS_beam)
# effective_voltage = voltage_program_1 * np.sin(harmonic_numbers_1*slicing.convert_coordinates(slicing.bins_centers, 'tau', 'theta')) + SPS_longitudinal_intensity.induced_voltage
# effective_voltage_interp = interp1d(slicing.bins_centers, effective_voltage)
# plt.plot(slicing.bins_centers, effective_voltage)
# plt.plot(slicing.bins_centers, voltage_program_1 * np.sin(harmonic_numbers_1*slicing.convert_coordinates(slicing.bins_centers, 'tau', 'theta')))
# plt.plot(slicing.bins_centers, SPS_longitudinal_intensity.induced_voltage)
# plt.show()
# sync_phase_intensity = brentq(effective_voltage_interp, 2e-9, 3e-9)
# SPS_beam.tau = SPS_beam.tau + (sync_phase_intensity - SPS_beam.bp_gauss_tau)
#        
# print (sync_phase_intensity - SPS_beam.bp_gauss)
         
# Total simulation map
sim_map = [full_tracker] + [slicing] + [SPS_longitudinal_intensity]# # + [SPS_intensity_time] + [SPS_intensity_freq] + [SPS_inductive]
         
# Tracking ---------------------------------------------------------------------

plt.ion()
plt.figure(3)
save_bl = np.zeros(n_turns)
save_bp = np.zeros(n_turns)
for i in range(n_turns):
                        
    if i % output_step == 0:
        print i
        t0 = time.clock()
                            
    # Track
    for m in sim_map:
        m.track(SPS_beam)
#         SPS_beam.longit_statistics()
                           
    if i % output_step == 0:
        t1 = time.clock()
#         print 4*SPS_beam.sigma_tau, SPS_beam.mean_tau
        print t1-t0
                                            
    if i % plot_step == 0:
#         plt.plot(SPS_beam.theta[0:10000], SPS_beam.dE[0:10000],'.')
#         plt.xlim((0,2*np.pi/harmonic_numbers_1))
#         plt.ylim((-1.50e8,1.50e8))
#         plt.pause(0.0001)
#         plt.clf()
#         plt.show()
        plt.plot(SPS_longitudinal_intensity.time_array, SPS_longitudinal_intensity.induced_voltage)
        plt.plot(SPS_longitudinal_intensity.time_array, slicing.n_macroparticles / np.max(slicing.n_macroparticles) *1e6)#* np.max(SPS_longitudinal_intensity.induced_voltage))
        plt.ylim((-2e6, 2e6))
#         plt.plot(slicing.bins_centers, slicing.n_macroparticles)
#         plt.ylim((0,6e3))
#         plt.plot(slicing.bins_centers, gauss(slicing.bins_centers, *slicing.pfit_gauss))
#         plt.plot(slicing.edges[:-1]-slicing.edges[1:])
        plt.pause(0.0001)
        plt.clf()
#         plt.show()
#         plot_long_phase_space(SPS_beam, general_params_dummy, rf_params_dummy, 
#                               0., 5., -150, 150, xunit='ns', separatrix_plot = 'True', sampling = 1e1)
       
    save_bl[i] = SPS_beam.bl_gauss_tau
    save_bp[i] = SPS_beam.bp_gauss_tau
#     save_bl[i] = 4*np.std(SPS_beam.tau)
#     save_bp[i] = np.mean(SPS_beam.tau)
       
    
plt.ioff()
     
plt.figure(1)
plt.plot(save_bl)
plt.figure(2)
plt.plot(save_bp)
plt.show()
