from scipy.constants import e, m_e, m_p, c
import numpy as np

# ==============================================================================================
# SIMULATION SETUP.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# PHYSICS AND MACHINE PARAMETERS.
intensity = 2.5e11                           # Number of particles (protons) per bunch.
charge    = e                               # Charge of a proton.
mass      = m_p                             # Mass of a proton.

sigma_z   = .23				                # Bunch length (RMS) [m].
gamma     = 27.7		                    # Relativistic gamma.
alpha_0   = 0.00192                      # Momentum compaction factor.
eta       = alpha_0 - 1./gamma**2           # Slippage factor.
gamma_t   = 1./np.sqrt(alpha_0)             # Transition gamma.

p0 = np.sqrt(gamma**2 - 1) * mass * c       # Momentum.

RFvoltage = 2.0e6							# RF voltage
Q_s       = 0.0059			                # Synchrotron tune.
Q_x       = 26.13                           # Betatron tune (horizontal).
Q_y       = 26.18                           # Betatron tune (vertical).

C         = 6911.                           # Ring circumference [m].
R         = C/(2.*np.pi)                    # Ring radius [m].

alpha_x   = 0.
alpha_y   = 0.

Qp_x      = 0.                              # Horizontal chromaticity.
Qp_y      = 0.                              # Vertical chromaticity.    
	
beta_x    = 42.		                    	# Horizontal beta function [m].
beta_y    = 42.                            	# Vertical beta function [m].
beta_z    = eta*R/Q_s                       # Longitudinal beta function [m].

epsn_x    = 2.5                             # Horizontal emittance [um].
epsn_y    = 2.5                             # Vertical emittance [um].
epsn_z    = 4.*np.pi*sigma_z**2 * p0 / (beta_z * e)

initial_kick_x = 0.0000                       # Initial horizontal kick of beam.
initial_kick_y = 0.0000                       # Initial vertical kick of beam.

# SIMULATION PARAMETERS.
n_macroparticles = 200000                    # Number of macroparticles per bunch (go to 1e6).
n_turns          = 256                       # Number of turn (set to 2e5 later)

n_segments = 19 

Dt_ref_ecloud = 100e-12

bunch_monitor_filename = 'test_separate_input'
