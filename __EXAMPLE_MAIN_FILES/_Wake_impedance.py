# SPS simulation with intensity effects in time and frequency domains using
# a table of resonators. The input beam has been cloned to show that the two methods
# are equivalent (compare the two figure folders). Note that to create an exact 
# clone of the beam, the option seed=0 in the generation has been used. This 
# script shows also an example of how to use the class SliceMonitor (check the
# corresponding h5 files).

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
from monitors.monitors import *
from beams.slices import *
from impedances.longitudinal_impedance import *
from longitudinal_plots.plot_beams import *
from longitudinal_plots.plot_impedance import *
from longitudinal_plots.plot_slices import *


# SIMULATION PARAMETERS -------------------------------------------------------

# Beam parameters
particle_type = 'proton'
n_particles = 1.e10        
n_macroparticles = 500000
tau_0 = 2.0

# Machine and RF parameters
gamma_transition = 1/np.sqrt(0.00192)   # [1]
C = 6911.56  
      
# Tracking details
n_turns = 2          
n_turns_between_two_plots = 1          

# Derived parameters
sync_momentum = 25.92e9 # [eV / c]
momentum_compaction = 1 / gamma_transition**2 # [1]       

# Cavities parameters
n_rf_systems = 1                                     
harmonic_number = 4620                         
voltage_program = 0.9e6
phi_offset = 0



# DEFINE RING------------------------------------------------------------------

general_params = GeneralParameters(n_turns, C, momentum_compaction, sync_momentum, 
                                   particle_type, number_of_sections = 1)
general_params_copy = GeneralParameters(n_turns, C, momentum_compaction, sync_momentum, 
                                   particle_type, number_of_sections = 1)

RF_sct_par = RFSectionParameters(general_params, n_rf_systems, harmonic_number, 
                          voltage_program, phi_offset)
RF_sct_par_copy = RFSectionParameters(general_params_copy, n_rf_systems, harmonic_number, 
                          voltage_program, phi_offset)

ring_RF_section = RingAndRFSection(RF_sct_par)
ring_RF_section_copy = RingAndRFSection(RF_sct_par_copy)

# DEFINE BEAM------------------------------------------------------------------

my_beam = Beam(general_params, n_macroparticles, n_particles)

my_beam_copy = Beam(general_params_copy, n_macroparticles, n_particles)

longitudinal_gaussian_matched(general_params, RF_sct_par, my_beam, tau_0, 
                              unit='ns', seed=0)

longitudinal_gaussian_matched(general_params_copy, RF_sct_par_copy, my_beam_copy, tau_0, 
                              unit='ns', seed=0)

number_slices = 100
slice_beam = Slices(my_beam, number_slices, cut_left = 0, 
                    cut_right = 2 * np.pi / harmonic_number, mode = 
                    'const_space', cuts_coord = 'theta', slicing_coord = 'tau', 
                    statistics_option = 'on', fit_option = 'gaussian')
slice_beam_copy = Slices(my_beam_copy, number_slices, cut_left = 0, 
                    cut_right = 2 * np.pi / harmonic_number, mode = 
                    'const_space', cuts_coord = 'theta', slicing_coord = 'tau', 
                    statistics_option = 'on', fit_option = 'gaussian')


# MONITOR----------------------------------------------------------------------

bunchmonitor = BunchMonitor('bunch', n_turns+1, "Longitudinal", slice_beam)
slicesmonitor = SlicesMonitor('slices', n_turns+1, slice_beam)
bunchmonitor.track(my_beam)
slicesmonitor.track(my_beam)

bunchmonitor_copy = BunchMonitor('bunch_copy', n_turns+1, "Longitudinal", slice_beam_copy)
slicesmonitor_copy = SlicesMonitor('slices_copy', n_turns+1, slice_beam_copy)
bunchmonitor_copy.track(my_beam_copy)
slicesmonitor_copy.track(my_beam_copy)

# LOAD IMPEDANCE TABLE--------------------------------------------------------

table = np.loadtxt('new_HQ_table.dat', comments = '!')
R_shunt = table[:, 2] * 10**6 
f_res = table[:, 0] * 10**9
Q_factor = table[:, 1]
resonator = Resonators(R_shunt, f_res, Q_factor)
ind_volt_time = InducedVoltageTime(slice_beam, [resonator])
ind_volt_freq = InducedVoltageFreq(slice_beam_copy, [resonator], 1e5)
tot_vol = TotalInducedVoltage(slice_beam, [ind_volt_time])
tot_vol_copy = TotalInducedVoltage(slice_beam_copy, [ind_volt_freq])

# ACCELERATION MAP-------------------------------------------------------------

map_ = [tot_vol] + [ring_RF_section] + [slice_beam] + [bunchmonitor] + [slicesmonitor]
map_copy = [tot_vol_copy] + [ring_RF_section_copy] + [slice_beam_copy] + [bunchmonitor_copy] + [slicesmonitor_copy]

# TRACKING + PLOTS-------------------------------------------------------------

for i in range(n_turns):
    
    print i+1
    for m in map_:
        m.track(my_beam)
    for m in map_copy:
        m.track(my_beam_copy)
    
    # Plots
    if ((i+1) % n_turns_between_two_plots) == 0:
        
        plot_induced_voltage_vs_bins_centers(i+1, general_params, tot_vol, style = '-', dirname = 'fig1a')
        plot_induced_voltage_vs_bins_centers(i+1, general_params_copy, tot_vol_copy, style = '-', dirname = 'fig1b')
        
        plot_long_phase_space(my_beam, general_params, RF_sct_par, 
          0, 0.0014, - 1.5e2, 1.5e2, sampling=50, dirname = 'fig1a')
        plot_long_phase_space(my_beam_copy, general_params_copy, RF_sct_par_copy, 
          0, 0.0014, - 1.5e2, 1.5e2, sampling=50, dirname = 'fig1b')
        
print "Done!"

bunchmonitor.h5file.close()
