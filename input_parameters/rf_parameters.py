'''
Created on 19 juin 2014

Module gathering and processing all the RF parameters to be given to the other modules

@author: Alexandre Lasheen, Danilo Quartullo
'''

from __future__ import division
import numpy as np

def input_check(input_value, expected_length):
    '''
    Function to check the length of the input
    If len(input_value) = 1, transform it to a constant array
    If len(input_value) != expected_length and != 1, raise an error
    '''
    
    if type(input_value) is float:
        return input_value * np.ones(expected_length)
    elif type(input_value) is int:
        return input_value * np.ones(expected_length)
    elif type(input_value) is np.ndarray and input_value.size == 1:
        return input_value * np.ones(expected_length)
    elif len(input_value) == expected_length:
        return input_value
    else:
        raise RuntimeError(str(input_value) + ' does not match ' + str(expected_length))
    
    
class Sum_RF_section_parameters(object):
    '''
    Method to add sections together
    '''
    
    def __init__(self, RF_section_parameters_list):
        
        self.RF_section_parameters_list = RF_section_parameters_list
        self.section_length_sum = 0
        self.total_n_sections = len(RF_section_parameters_list)
        self.momentum_program_matrix = np.zeros((self.total_n_sections, RF_section_parameters_list[0].n_turns + 1))

        for i in range(len(RF_section_parameters_list)):
            self.section_length_sum += RF_section_parameters_list[i].section_length
            self.momentum_program_matrix[i,:] = RF_section_parameters_list[i].momentum_program


class RF_section_parameters(object):
    '''
    Object containing all the RF parameters for one section
    It can be simply added to another RF_section_parameters object to have
    a complete set of parameters to be used afterwards
    '''
    
    def __init__(self, n_turns, n_rf_systems, section_length, harmonic_numbers_list, voltage_program_list, phi_offset_list, momentum_program):
        
        self.n_sections = 1
        self.n_turns = n_turns
        self.section_length = section_length #: Length of the section in [m]
        self.n_rf_systems = n_rf_systems #: Number of RF systems in the section
        self.momentum_program = momentum_program #: Momentum program in [eV/c]
        
        self.momentum_program = input_check(momentum_program, self.n_turns + 1) #: Momentum program in [eV/c]
        
        self.p_increment = self.momentum_program[1:] - self.momentum_program[0:-1] #: Momentum change (acceleration/deceleration) between two turns
        
        if n_rf_systems == 1:
            self.harmonic_numbers_list = [harmonic_numbers_list] #: Harmonic number list (according to n_rf_systems)
            self.voltage_program_list = [voltage_program_list] #: Voltage program list in [V] (according to n_rf_systems)
            self.phi_offset_list = [phi_offset_list] #: Phase offset list in [rad]
        else:
            if not n_rf_systems == len(harmonic_numbers_list) == len(voltage_program_list) == len(phi_offset_list):
                raise RuntimeError('The RF parameters to define RF_section_parameters are not homogeneous')
            self.harmonic_numbers_list = harmonic_numbers_list #: Harmonic number list (according to n_rf_systems)
            self.voltage_program_list = voltage_program_list #: Voltage program list in [V] (according to n_rf_systems)
            self.phi_offset_list = phi_offset_list #: Phase offset list in [rad]
        
        for i in range(self.n_rf_systems):
            self.harmonic_numbers_list[i] = input_check(self.harmonic_numbers_list[i], n_turns)
            self.voltage_program_list[i] = input_check(self.voltage_program_list[i], n_turns)
            self.phi_offset_list[i] = input_check(self.phi_offset_list[i], n_turns)
            

            
            

            
            
    
