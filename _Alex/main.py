'''
SPS test case
'''

# General imports
from __future__ import division
import numpy as np
import time

# PyHEADTAIL imports
from input_parameters.general_parameters import GeneralParameters 
from input_parameters.rf_parameters import RFSectionParameters
from beams.beams import Beam
from beams.longitudinal_distributions import longitudinal_gaussian_matched
from trackers.longitudinal_tracker import RingAndRFSection
from beams.plot_beams import plot_long_phase_space

# Simulation parameters -------------------------------------------------------
# Simulation parameters
n_turns = 2000          # Number of turns to track
plot_step = 100         # Time steps between plots

# General parameters
particle_type = 'proton'
circumference = 6911.56                         # Machine circumference [m]
gamma_transition = 1/np.sqrt(0.00192)           # Transition gamma
momentum_compaction = 1./gamma_transition**2   # Momentum compaction array
sync_momentum = 25.92e9                         # Synchronous momentum program [eV/c]

# RF parameters
n_rf_systems = 2                # Number of rf systems second section
harmonic_numbers_1 = 4620       # Harmonic number second section
voltage_program_1 = 0.9e6       # RF voltage [V] second section
phi_offset_1 = 0                # Phase offset
harmonic_numbers_2 = 4620*4     # Harmonic number second section
voltage_program_2 = 0.e6        # RF voltage [V] second section
phi_offset_2 = np.pi            # Phase offset

# Beam parameters
intensity = 1.e10           # Intensity
n_macroparticles = 100000   # Macro-particles
tau_0 = 2.0                 # Initial bunch length, 4 sigma [ns]          


# Simulation setup ------------------------------------------------------------
# General parameters                  
general_params = GeneralParameters(n_turns, circumference, momentum_compaction, 
                                   sync_momentum, particle_type)

# RF parameters
rf_params = RFSectionParameters(general_params, 1, harmonic_numbers_1, 
                                      voltage_program_1, phi_offset_1)
# rf_params = RFSectionParameters(general_params, n_rf_systems, [harmonic_numbers_1, harmonic_numbers_2], 
#                                 [voltage_program_1, voltage_program_2], [phi_offset_1, phi_offset_2])

# RF tracker
longitudinal_tracker = RingAndRFSection(rf_params)

# Beam
my_beam = Beam(general_params, n_macroparticles, intensity)
rf_params_dummy = RFSectionParameters(general_params, 1, harmonic_numbers_1, 
                                      voltage_program_1, phi_offset_1)
longitudinal_gaussian_matched(general_params, rf_params_dummy, my_beam, tau_0, 
                              unit='ns')

# Total simulation map
sim_map = [longitudinal_tracker]


# Tracking ---------------------------------------------------------------------

for general_params.counter[0] in range(n_turns):
    i = general_params.counter[0]
    
    if i % 100 == 0:
        print i
        t0 = time.clock()
        
    # Track
    for m in sim_map:
        m.track(my_beam)
       
    if i % 100 == 0:
        t1 = time.clock()
        print t1-t0
       
    if i % plot_step == 0:
        plot_long_phase_space(my_beam, general_params, rf_params_dummy, 
                              0., 5., -150, 150, xunit='ns', separatrix_plot = 'True')










