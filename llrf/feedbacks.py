'''
**Various feedbacks**

:Authors: **Helga Timko**
'''

from __future__ import division
from scipy.constants import c
import numpy as np



class PhaseLoop(object): 
    '''
    One-turn phase loop for different machines with different hardware. Phase 
    loop is active only in certain turns, in which the frequency or the phase
    in the longitudinal tracker is modified. 
    '''    
    def __init__(self, general_params, gain, sampling_frequency = 100,
                 machine = 'LHC'):
        
        self.gain = gain # feedback gain, can be an array depending on machine
        self.dt = sampling_frequency # either in turns or in time, depending on machine
        self.machine = machine # machine name
        
        self.correction = 0 # PL correction in frequency or phase, depending on machine
        self.dphi = 0 # phase difference between bunch/beam and RF
        self.PL_on = False # time step when PL is active
        
        # Pre-processing
        if self.machine == 'PSB':
            self.counter = 0
            self.precalculate_time(general_params)
       
            
    def track(self, beam, tracker):
        if (self.machine == 'LHC' and tracker.counter % self.dt) \
        or (self.machine == 'PSB' and tracker.counter == self.on_time[self.counter]):
            getattr(self, self.machine)()
            self.PL_on = True
        if (self.machine == 'LHC' and tracker.counter % (self.dt+1)) \
        or (self.machine == 'PSB' and tracker.counter == (self.on_time[self.counter]+1)):
            self.PL_on = True
        else:
            self.PL_on = False

    
    def precalculate_time(self, general_params):
        n = 0
        while n < general_params.t_rev.size:
            summa = 0
            while summa < self.dt:
                summa += general_params.t_rev[n]
            np.append(self.on_time, n)
            
      
    def phase_difference(self, beam):
        self.dphi = h * beam.mean_theta - phi_s
        
    def LHC(self, beam):
        phase_difference(beam)
        self.correction = multiply by k

    def PSB(self, beam):
        another transfer function
        self.counter += 1 
    
    

class LHCNoiseFB(object): 
    
    # Feedback based on bunch length, acting on phase noise used for blow-up 
    def __init__(self, bl_target, gain = 0.1, factor = 0.64, 
                 sampling_frequency = 100, self_statistics = False ):

        self.x = 1 # multiplication factor; initially 1
        self.bl_targ = bl_target # 4 sigma, in s
        self.bl_meas = bl_target # set initial value for measured bunch length
        self.g = gain # in inverse-turns
        self.a = factor
        self.dt = sampling_frequency # bunch length sampling frequency, in turns
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
        
        