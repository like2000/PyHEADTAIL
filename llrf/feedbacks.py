'''
**Various feedbacks**

:Authors: **Helga Timko**
'''

from __future__ import division
from scipy.constants import c
import numpy as np



class LHCNoiseFB(object): 
    
    # Feedback based on bunch length, acting on phase noise used for blow-up 
    def __init__(self, bl_target, gain=0.1, factor=0.64, dt_update = 100, 
                 self_statistics=False ):

        self.x = 1 # multiplication factor; initially 1
        self.bl_targ = bl_target # 4 sigma, in s
        self.bl_meas = bl_target # set initial value for measured bunch length
        self.g = gain # in inverse-turns
        self.a = factor
        self.dt = dt_update # update frequency, in turns
        self.self_stat = self_statistics # using external/internal statistics
       

    def FB(self, general_params, rf_params, beam, RFnoise):
        # Update bunch length, every dt no. of turns
        if rf_params.counter[0] % self.dt:
            if self.self_stat:
                itemindex = np.where(beam.id != 0)[0]
                self.bl_meas = 4.*np.std(beam.theta[itemindex]) \
                             *beam.ring_radius/(beam.beta_r*c)
            else:
                self.bl_meas = 4.*beam.sigma_tau

        
        # Update multiplication factor
        self.x = self.a*self.x \
               + self.g*(self.bl_targ - self.bl_meas)/general_params.t_rev
        
        # Limit to range [0,1]
        if self.x < 0:
            self.x = 0
        if self.x > 1:
            self.x = 1
            
        
        # Update phase noise for the following turn
        RFnoise.dphi[rf_params.counter[0] + 1] *= self.x
        
        