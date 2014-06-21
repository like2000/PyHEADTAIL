'''
Created on 19 juin 2014

Module gathering all the global input parameters to be called and used in other modules for the simulation

@author: Alexandre Lasheen, Danilo Quartullo
'''

from __future__ import division
from scipy.constants import m_p, e, c
import numpy as np


class General_parameters(object):
    '''
    Object containing all the global parameters
    '''

    def __init__(self, particle_type, n_turns, ring_circumference, alpha_array, momentum_program):
        
        self.particle_type = particle_type #: Defining particle type (mass in [kg] and charge in [C])
        if self.particle_type is 'proton':
            self.mass = m_p
            self.charge = e
        else:
            raise RuntimeError('Particle type not recognized')
        
        self.counter = [0] #: Counter to be incremented every turn (initiated as a list to be called as a pointer)
        self.n_turns = n_turns #: Number of turns of the simulation
        
        self.ring_circumference = ring_circumference #: Ring circumference in [m]
        self.ring_radius = self.ring_circumference / (2*np.pi) #: Ring circumference in [m]
        
        self.momentum_program = momentum_program #: Momentum program in [eV/c]
        
        self.alpha_array = np.array(alpha_array) #: Momentum compation (up to 2nd order)
        
        '''
        The momentum program can be given as a single value (initial momentum, no acceleration)
        or as a full program (array of length n_turns+1)
        For several RF sections, this program has to be redefined for each ones and only the initial
        momentum is needed here
        '''
        if type(self.momentum_program) is float:
            self.momentum_program = self.momentum_program * np.ones(self.n_turns + 1) #: Momentum program in [eV/c], if length is 1 : constant value
        elif type(self.momentum_program) is np.ndarray:
            try:
                assert self.momentum_program.shape[1] == self.n_turns + 1
            except IndexError:
                assert self.momentum_program.shape[0] == self.n_turns + 1
            
            
        self.beta_rel_program = np.sqrt( 1 / (1 + (self.mass * c**2)**2 / (self.momentum_program * e)**2) ) #: Relativistic beta program
        
        self.gamma_rel_program = np.sqrt( 1 + (self.momentum_program * e)**2 / (self.mass * c**2)**2 ) #: Relativistic gamma program
        
        self.energy_program = np.sqrt( self.momentum_program**2 + (self.mass * c**2 / e)**2 ) #: Energy program in [eV]
        
        self.eta0 = 0
        self.eta1 = 0
        self.eta2 = 0
        
        self.eta_generation()
                
                
    def eta_generation(self):
        
        """
        Pre-processing of the eta parameters with respect to the input momentum
        compaction factor (array) and the momentum program
        For eta coefficients, see Lee: Accelerator Physics (Wiley).
        """

        for i in xrange( self.alpha_array.size ):   # order = len - 1
            getattr(self, '_eta' + str(i))()

    
    def _eta0(self):
        
        if self.alpha_array.size == 1:
            self.eta0 = self.alpha_array - self.gamma_rel_program**-2 #: Slippage factor (order 0)
            return self.eta0
        else:
            self.eta0 = self.alpha_array[0] - self.gamma_rel_program**-2 #: Slippage factor (order 0)
            return self.eta0
   
    
    def _eta1(self):
        
        self.eta1 = 3 * self.beta_rel_program**2 / (2 * self.gamma_rel_program**2)  \
                + self.alpha_array[1] - self.alpha_array[0] * (self.alpha_array[0] - self.gamma_rel_program**-2)
        return self.eta1
    
    
    def _eta2(self):
        
        self.eta2 = - self.beta_rel_program**2 * (5 * self.beta_rel_program**2 - 1) / (2 * self.gamma_rel_program**2) \
                    + self.alpha_array[2] - 2 * self.alpha_array[0] * self.alpha_array[1] + self.alpha_array[1] \
                    / self.gamma_rel_program**2 + self.alpha_array[0]**2 * (self.alpha_array[0] - self.gamma_rel_program**-2) \
                    - 3 * self.beta_rel_program**2 * self.alpha_array[0] / (2 * self.gamma_rel_program**2)
        return self.eta2
        
