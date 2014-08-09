from __future__ import division
import numpy as np
import pylab as pl
import time


import ecloud.PyECLOUD_for_PyHEADTAIL as pyecl
from particles.particles import *
from scipy.constants import e, m_e
import numpy as np
from particles.slicer import *
from trackers.transverse_tracker_2 import *
from trackers.longitudinal_tracker import *
from monitors.monitors import *

from simulation_parameters import *


simulation_parameters_dict = {'comment': 'test instab ecloud',\
                                  'charge':           charge,\
                                  'mass':             mass,\
                                  'intensity':        intensity,\
                                  'beta_x':           beta_x,\
                                  'beta_y':           beta_y,\
                                  'beta_z':           beta_z,\
                                  'sigma_z':          sigma_z,\
                                  'epsn_x':           epsn_x,\
                                  'epsn_y':           epsn_y,\
                                  'gamma_t':          gamma_t,\
                                  'C':                C,\
                                  'n_turns':          n_turns,\
                                  'Q_x':              Q_x,\
                                  'Q_y':              Q_y,\
                                  'Q_s':              Q_s,\
                                  'Qp_x':             Qp_x,\
                                  'Qp_y':             Qp_y,\
                                  'n_macroparticles': n_macroparticles,\
                                  'initial_kick_x':   initial_kick_x,\
                                  'initial_kick_y':   initial_kick_y,\
                                  }


L_ecloud = C/n_segments



# Beam
bunch = Particles.as_gaussian(n_macroparticles, e, gamma, intensity, m_p, alpha_x, beta_x, epsn_x, alpha_y, beta_y, epsn_y, beta_z, epsn_z)

bunch.x += initial_kick_x
bunch.y += initial_kick_y


#ecloud
beamslicer = Slicer(50, nsigmaz=2)
#ecloud = pyecl.Ecloud(L_ecloud, beamslicer, Dt_ref = Dt_ref_ecloud, pyecl_input_folder='drift_for_instab_SPS')

s = np.arange(0, n_segments + 1) * C/n_segments

# BETATRON
# Loop on number of segments and create the TransverseSegmentMap for each segment.
alpha_x *= np.ones(n_segments)
beta_x  *= np.ones(n_segments)
D_x      = np.zeros(n_segments)
alpha_y *= np.ones(n_segments)
beta_y  *= np.ones(n_segments)
D_y      = np.zeros(n_segments)

# Create detuning elements.
#chromaticity       = Chromaticity(Qp_x, Qp_y)
#~ amplitude_detuning = AmplitudeDetuning.from_octupole_currents_LHC(i_oct_f, i_oct_d)

# Generate transverse map.
transverse_map = TransverseMap(s, alpha_x, beta_x, D_x, alpha_y, beta_y, D_y, Q_x, Q_y)#,
							   #chromaticity)#, amplitude_detuning)


# SYNCHROTRON
#cavity = LinearMap(C, alpha_0, Q_s)
cavity = RFSystems(C, [4620], [RFvoltage], [0.], [alpha_0], bunch.gamma)

print '{0:4d} \t {1:+3e} \t {2:+3e} \t {3:+3e} \t {4:3e} \t {5:3e} \t {6:3f} \t {7:3f} \t {8:3f} \t {9:4e} \t {10:3s}'.format(-1, bunch.mean_x(), bunch.mean_y(), bunch.mean_z(), bunch.epsn_x(), bunch.epsn_y(), bunch.epsn_z(), bunch.sigma_z(), bunch.sigma_dp(), bunch.n_macroparticles / bunch.n_macroparticles * bunch.intensity, 'c')



# build ring
elements=[]
for l in transverse_map:
    elements+=[l]#, ecloud]

elements.append(cavity)
bunch_monitor = BunchMonitor(bunch_monitor_filename, n_turns, simulation_parameters_dict)

for i in range(n_turns):
        t0 = time.clock()
	
        for ind, m in enumerate(elements):
			m.track(bunch)
			print ind, m
			

        # slice_monitor.dump(bunch)
        bunch_monitor.dump(bunch)
        # particle_monitor.dump(bunch)
      
        print '{0:4d} \t {1:+3e} \t {2:+3e} \t {3:+3e} \t {4:3e} \t {5:3e} \t {6:3f} \t {7:3f} \t {8:3f} \t {9:4e} \t {10:3s}'.format(i, bunch.mean_x(), bunch.mean_y(), bunch.mean_z(), bunch.epsn_x(), bunch.epsn_y(), bunch.epsn_z(), bunch.sigma_z(), bunch.sigma_dp(), bunch.n_macroparticles / bunch.n_macroparticles * bunch.intensity, str(time.clock() - t0))

bunch_monitor.close()          
