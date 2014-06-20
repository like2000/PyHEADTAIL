'''
Test main file
'''

from __future__ import division
import numpy as np

from input_parameters.simulation_parameters import General_parameters
from input_parameters.rf_parameters import RF_section_parameters, sum_RF_section_parameters
from trackers.ring_and_RFstation import Full_Ring_and_RF
from trackers.longitudinal_tracker import Longitudinal_tracker

from os import chdir
chdir('C:/work/git/PyHEADTAIL/ALEX_input&main_files')


# Simulation parameters --------------------------------------------------------
# Simulation parameters
n_turns = 1000000                                    # Number of turns to track

# Global parameters
particle_type = 'proton'
circumference = 6911.56                         # Machine circumference [m]
gamma_transition = 1/np.sqrt(0.00192)           # Transition gamma
momentum_compaction = 1./gamma_transition**2    # Momentum compaction array
 

# RF parameters
n_rf_systems_1 = 2                                      # Number of rf systems first section
harmonic_numbers_1_1 = np.array([4620])           # Harmonic number first section, first RF system
harmonic_numbers_1_2 = np.array([4620*4])           # Harmonic number first section, second RF system
harmonic_numbers_1_list = [harmonic_numbers_1_1, harmonic_numbers_1_2]
voltage_program_1_1 = 0.9e6 * np.ones([n_turns])        # RF voltage [V] first section, first RF system
voltage_program_1_2 = np.array([0.e6])                  # RF voltage [V] first section, second RF system
voltage_program_1_list = [voltage_program_1_1, voltage_program_1_2]
phi_offset_1_1 = 0 * np.ones([n_turns])                 # Phase offset first section, first RF system
phi_offset_1_2 = 0 * np.ones([n_turns])                 # Phase offset first section, second RF system
phi_offset_1_list = [phi_offset_1_1, phi_offset_1_2]
sync_momentum_1 = 25.92e9                               # Synchronous momentum [eV/c] first section
section_1_params = RF_section_parameters(n_turns, n_rf_systems_1, circumference/2, harmonic_numbers_1_list, voltage_program_1_list, phi_offset_1_list, sync_momentum_1)
    
n_rf_systems_2 = 1                                      # Number of rf systems second section
harmonic_numbers_2 = 4620                               # Harmonic number second section
voltage_program_2 = 0.9e6                               # RF voltage [V] second section
sync_momentum_2 = 25.92e9 * np.ones(n_turns+1)          # Synchronous momentum program [eV/c] second section
phi_offset_2 = 0
section_2_params = RF_section_parameters(n_turns, n_rf_systems_2, circumference/2, harmonic_numbers_2, voltage_program_2, phi_offset_2, sync_momentum_2)

rf_params = sum_RF_section_parameters([section_1_params, section_2_params])
global_params = General_parameters(particle_type, n_turns, circumference, momentum_compaction, rf_params.momentum_program_matrix)







