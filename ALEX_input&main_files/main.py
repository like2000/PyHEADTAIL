'''
Test main file
'''

from __future__ import division
import numpy as np
import time
from scipy.constants import m_p, e

from trackers.ring_and_RFstation import Ring_and_RFstation
from beams.beams import Beam
from trackers.longitudinal_tracker import Longitudinal_tracker
from beams.longitudinal_distributions import longitudinal_gaussian_matched
from longitudinal_plots.longitudinal_plots import plot_long_phase_space


# Simulation parameters --------------------------------------------------------
# Tracking details
n_turns = 2000      # Number of turns to track
plot_step = 10      # Time steps between plots
output_step = 100     # Time steps between outputs

# Bunch parameters
intensity = 1.e10           # Intensity
n_macroparticles = 100001    # Macro-particles
tau_0 = 2.                  # Initial bunch length, 4 sigma [ns]

# Machine and RF parameters
harmonic_number_array = np.array([4620])    # Harmonic number
voltage_program = np.array([0.9e6])         # RF voltage [V]
gamma_transition = 1/np.sqrt(0.00192)       # Transition gamma
circumference = 6911.56                     # Machine circumference [m]
sync_momentum = 25.92e9                     # Synchronous momentum [eV]

# Derived parameters
momentum_program = sync_momentum * np.ones(n_turns+1)
momentum_compaction = np.array([1./gamma_transition**2])


# Simulation setup -------------------------------------------------------------
print "Setting up the simulation..."
print ""

# Define Ring and RF Station
SPS_ring = Ring_and_RFstation(circumference, momentum_program, momentum_compaction, circumference, harmonic_number_array, np.array([0.4e6]))
print "Ring and RF are set..."
print ""

# Define Beam
SPS_bunch = Beam(SPS_ring, m_p, n_macroparticles, e, intensity)
longitudinal_gaussian_matched(SPS_bunch, tau_0, unit='ns')
print "Beam is set..."
print ""
 
SPS_ring.voltage =  np.array([0.9e6])
# SPS_bunch.ring = SPS_ring
 
print SPS_bunch.ring.voltage
 
# Define RF
SPS_long_tracker = Longitudinal_tracker(SPS_ring)
print "Longitudinal tracker is set..."
print ""
   
# Accelerator map
map_ = [SPS_long_tracker]
print "Map is set..."
print ""
   
  
# Tracking ---------------------------------------------------------------------
for i in range(n_turns):
    t0 = time.clock()
       
    # Track
    for m in map_:
        m.track(SPS_bunch)
         
    t1 = time.clock()
     
    # Output
    if (i % output_step) == 0:
        print t1-t0
         
    # Plot
    if (i % plot_step) == 0:
        plot_long_phase_space(SPS_bunch, i, 0, 5, -75, 75, xunit = 'ns')












