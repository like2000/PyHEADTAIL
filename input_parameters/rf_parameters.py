'''
Created on 19 juin 2014

Module gathering and processing all the RF parameters to be given to the other modules

@author: Alexandre Lasheen
'''

import numpy as np

def input_check(input_value, expected_length):
    '''
    Function to check the length of the input
    If len(input_value) = 1, transform it to a constant array
    If len(input_value) != expected_length and != 1, raise an error
    '''
    
    if type(input_value) is float:
        return input_value * np.ones(expected_length)
    if type(input_value) is int:
        return input_value * np.ones(expected_length)
    elif type(input_value) is np.ndarray and input_value.size == 1:
        return input_value * np.ones(expected_length)
    elif len(input_value) == expected_length:
        return input_value
    else:
        raise RuntimeError(str(input_value) + ' does not match ' + str(expected_length))   


class RF_section_parameters(object):
    '''
    Object containing all the RF parameters for one section
    It can be simply added to another RF_section_parameters object to have
    a complete set of parameters to be used afterwards
    '''
    
    def __init__(self, n_rf_systems, harmonic_numbers_list, voltage_program_list, phi_offset_list, momentum_program):
        
        if n_rf_systems == 1:
            self.n_rf_systems = n_rf_systems #: Number of RF systems in the section
            self.harmonic_numbers_list = [harmonic_numbers_list] #: Harmonic number list (according to n_rf_systems)
            self.voltage_program_list = [voltage_program_list] #: Voltage program list in [V] (according to n_rf_systems)
            self.phi_offset_list = [phi_offset_list] #: Phase offset list in [rad]
            self.momentum_program = momentum_program #: Momentum program in [eV/c]
            
        else:
            if not n_rf_systems == len(harmonic_numbers_list) == len(voltage_program_list) == len(phi_offset_list):
                raise RuntimeError('The RF parameters to define RF_section_parameters are not homogeneous')
        
            self.n_rf_systems = n_rf_systems #: Number of RF systems in the section
            self.harmonic_numbers_list = harmonic_numbers_list #: Harmonic number list (according to n_rf_systems)
            self.voltage_program_list = voltage_program_list #: Voltage program list in [V] (according to n_rf_systems)
            self.phi_offset_list = phi_offset_list #: Phase offset list in [rad]
            self.momentum_program = momentum_program #: Momentum program in [eV/c]
        

class RF_parameters(object):
    '''
    Object containing all the RF parameters and generating arrays to be used
    in the tracking
    '''

    def __init__(self, n_turns, n_sections, RF_section_parameters_list):

        self.n_sections = n_sections #: Number of sections of the ring
        self.n_turns = n_turns #: Number of turns for the simulation
        
        # Checking the validity of the input and warning
        if len(RF_section_parameters_list) != n_sections:
            raise RuntimeError(str(RF_section_parameters_list) + ' does not match with the number of sections')
        
        self.momentum_program = np.zeros((self.n_sections, n_turns + 1)) #: Reorganizing the momentum_program as a matrix
        rf_input_section_arr = [] #: Reorganizing the rf_input as a list of np.arrays, for one section
        self.rf_input_total_arr = [] #: Gathering all the np.arrays of all the sections in one list
        
        for i in range(n_sections):
            self.momentum_program[i,:] = input_check(RF_section_parameters_list[i].momentum_program, n_turns + 1)
                
            for j in range(RF_section_parameters_list[i].n_rf_systems):
                harmonic_number_array = input_check(RF_section_parameters_list[i].harmonic_numbers_list[j], n_turns)
                voltage_array = input_check(RF_section_parameters_list[i].voltage_program_list[j], n_turns)
                phase_offset_array = input_check(RF_section_parameters_list[i].phi_offset_list[j], n_turns)
                rf_input_section_arr.append([harmonic_number_array, voltage_array, phase_offset_array])
            
            self.rf_input_total_arr.append(rf_input_section_arr)
            
            rf_input_section_arr = []
            

        
        
            
            
            
            
            
            
            
            
            
            
            
            
            