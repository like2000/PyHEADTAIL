'''
**Module gathering all the general input parameters used for the simulation**

:Authors: **Alexandre Lasheen**, **Danilo Quartullo**

'''

from __future__ import division
from scipy.constants import m_p, e, c
import numpy as np


class GeneralParameters(object):
    '''
    *Object containing all the general input parameters used for the simulation*
    '''

    def __init__(self, particle_type, n_turns, ring_circumference, alpha_array, 
                 momentum_program):
        
        #: *Particle type*
        self.particle_type = particle_type
        
        #: *Particle mass in [kg]*
        self.mass = 0
        
        #: *Particle charge in [C]*
        self.charge = 0 
        
        #: | *Counter to be incremented every turn in the tracking loop.*
        #: | *It is defined as a list in order to be passed as a reference in the other modules.*
        self.counter = [0]
        
        #: *Number of turns of the simulation*
        self.n_turns = n_turns 
        
        #: *Ring circumference in [m]*
        self.ring_circumference = ring_circumference 
        
        #: *Ring radius in [m]*
        self.ring_radius = self.ring_circumference / (2*np.pi) 
        
        #: | *Momentum program in [eV/c]*
        #: | *Can be given as a single value to be assumed constant, or as a program of (n_turns + 1) terms in case of acceleration.*
        self.momentum_program = momentum_program 
        
        #: *Momentum compation factor (up to 2nd order)*
        self.alpha_array = np.array(alpha_array) 
        
        #: *Relativistic beta program*
        self.beta_rel_program = np.sqrt(1 / (1 + (self.mass * c**2)**2 / 
                                             (self.momentum_program * e)**2)) 
        
        #: *Relativistic gamma program*
        self.gamma_rel_program = np.sqrt(1 + (self.momentum_program * e)**2 / 
                                         (self.mass * c**2)**2) 
        
        #: *Energy program in [eV]*
        self.energy_program = np.sqrt(self.momentum_program**2 + 
                                      (self.mass * c**2 / e)**2) 
        
        #: *Slippage factor (order 0)*
        self.eta0 = 0
        
        #: *Slippage factor (order 1)*
        self.eta1 = 0
        
        #: *Slippage factor (order 2)*
        self.eta2 = 0
        
        ### Pre-processing
        # Attribution of mass and charge with respect to particle_type
        if self.particle_type is 'proton':
            self.mass = m_p
            self.charge = e
        else:
            raise RuntimeError('Particle type not recognized')
        
        # Warning that higher orders for alpha_array will not be used
        if ((isinstance(self.alpha_array, np.ndarray) or 
             isinstance(self.alpha_array, list)) and len(self.alpha_array) > 3):
            print 'WARNING : Momentum compaction factor is held only up to 2nd \
            order'
        
        # Processing the momentum program for the cases where it is constant
        # or input as a program
        if isinstance(self.momentum_program, float):
            self.momentum_program = self.momentum_program * np.ones(self.n_turns + 1)
        elif isinstance(self.momentum_program, np.ndarray):
            try:
                assert self.momentum_program.shape[1] == self.n_turns + 1
            except IndexError:
                assert self.momentum_program.shape[0] == self.n_turns + 1
        
        # Processing the slippage factor
        self.eta_generation()
                
                
    def eta_generation(self):
        '''
        *Pre-processing of the slippage factor parameters with respect to the
        input momentum compaction factor (up to 2nd order) and the momentum program.
        For eta coefficients, see Lee: Accelerator Physics (Wiley).*
        '''
        
        for i in xrange(self.alpha_array.size):
            getattr(self, '_eta' + str(i))()

    
    def _eta0(self):
        '''
        *Calculation of the slippage factor (order 0) with respect to the
        momentum program and momentum compaction factor.*
        ''' 
        
        if self.alpha_array.size == 1:
            self.eta0 = self.alpha_array - self.gamma_rel_program**-2
            return self.eta0
        else:
            self.eta0 = self.alpha_array[0] - self.gamma_rel_program**-2
            return self.eta0
   
    
    def _eta1(self):
        '''
        *Calculation of the slippage factor (order 1) with respect to the
        momentum program and momentum compaction factor.*
        ''' 
                
        self.eta1 = 3 * self.beta_rel_program**2 / (2 * self.gamma_rel_program**2) + \
                    self.alpha_array[1] - self.alpha_array[0] * (self.alpha_array[0] - self.gamma_rel_program**-2)
        return self.eta1
    
    
    def _eta2(self):
        '''
        *Calculation of the slippage factor (order 2) with respect to the
        momentum program and momentum compaction factor.*
        ''' 
                
        self.eta2 = - self.beta_rel_program**2 * (5 * self.beta_rel_program**2 - 1) / (2 * self.gamma_rel_program**2) + \
                    self.alpha_array[2] - 2 * self.alpha_array[0] * self.alpha_array[1] + \
                    self.alpha_array[1] / self.gamma_rel_program**2 + \
                    self.alpha_array[0]**2 * (self.alpha_array[0] - self.gamma_rel_program**-2) - \
                    3 * self.beta_rel_program**2 * self.alpha_array[0] / (2 * self.gamma_rel_program**2)
        return self.eta2
        
