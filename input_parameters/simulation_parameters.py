'''
Created on 19 juin 2014

Module gathering all the global parameters to be called and used for the simulation

@author: Alexandre Lasheen
'''

from scipy.constants import m_p, e, c
import numpy as np

class Global_parameters(object):
    '''
    Object containing all the global parameters
    '''

    def __init__(self, particle_type, n_turns, ring_circumference, momentum_compaction_array, momentum_program):
        
        if particle_type is 'proton':
            self.mass = m_p
            self.charge = e
        else:
            raise RuntimeError('Particle type not recognized')
        
        self.counter = 0 #: Counter to be incremented every turn
        self.n_turns = n_turns #: Number of turns of the simulation
        
        self.ring_circumference = ring_circumference #: Ring circumference in [m]
        
        self.momentum_compaction_array = momentum_compaction_array #: Momentum compation (up to 2nd order)
        
        if type(momentum_program) is float:
            self.momentum_program = momentum_program * np.ones(self.n_turns + 1) #: Momentum program in [eV/c], if length is 1 : constant value
        else:
            assert len(momentum_program) == self.n_turns + 1
            
        self.beta_rel_program = np.sqrt( 1 / (1 + (self.mass * c**2)**2 / (self.momentum_program * e)**2) ) #: Relativistic beta program
        
        self.gamma_rel_program = np.sqrt( 1 + (self.momentum_program * e)**2 / (self.mass * c**2)**2 ) #: Relativistic gamma program
        
        self.energy_program = np.sqrt( self.momentum_program**2 + (self.mass * c**2 / e)**2 ) #: Energy program in [eV]
            
            
            
            
            