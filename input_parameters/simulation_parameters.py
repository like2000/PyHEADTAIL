'''
Created on 19 juin 2014

Module gathering all the gloabal parameters to be called and used for the simulation

@author: Alexandre Lasheen
'''

from scipy.constants import m_p, e
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
            self.momentum_program = momentum_program * np.ones(self.n_turns + 1) #: Momentum program, if length is 1 : constant value
        else:
            assert len(momentum_program) == self.n_turns + 1
            
            
            
            
            
            