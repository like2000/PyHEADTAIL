'''
**Module gathering and processing all the RF parameters to be given to the 
other modules.**

:Authors: **Alexandre Lasheen**, **Danilo Quartullo**

'''

from __future__ import division
import numpy as np


def input_check(input_value, expected_length):
    '''
    | *Function to check the length of the input*
    | *The input can be a float, int, np.ndarray and list*
    | *If len(input_value) == 1, transform it to a constant array*
    | *If len(input_value) != expected_length and != 1, raise an error*
    '''
    
    if isinstance(input_value, float):
        return input_value * np.ones(expected_length)
    elif isinstance(input_value, int):
        return input_value * np.ones(expected_length)
    elif isinstance(input_value, np.ndarray) and input_value.size == 1:
        return input_value * np.ones(expected_length)
    elif isinstance(input_value, list) and len(input_value) == 1:
        return input_value[0] * np.ones(expected_length)
    elif len(input_value) == expected_length:
        return np.array(input_value)
    else:
        raise RuntimeError(str(input_value) + ' does not match ' + str(expected_length))
    
    
class Sum_RF_section_parameters(object):
    '''
    *Method to add RF_section_parameters objects together in order to gather
    the complete information for a longitudinal_tracker.Full_Ring_and_RF
    object*
    '''
    
    def __init__(self, RF_section_parameters_list):
        
        #: *List of RF_section_parameters objects to concatenate*
        self.RF_section_parameters_list = RF_section_parameters_list
        
        #: *Total length of the sections in [m]*
        self.section_length_sum = 0
        
        #: *The total number of sections concatenated*
        self.total_n_sections = len(RF_section_parameters_list)
        
        #: | *Momentum program matrix in [eV/c]*
        #: | *The lines of this matrix corresponds to the momentum program for one section.*
        #: | *The columns correspond to one turn of the simulation.* 
        self.momentum_program_matrix = np.zeros((self.total_n_sections, 
                                                 RF_section_parameters_list[0].n_turns + 1))
        
        ### Pre-processing the inputs
        # The length of the sections are added and the momentum program is 
        # set as a matrix.
        for i in range(len(RF_section_parameters_list)):
            self.section_length_sum += RF_section_parameters_list[i].section_length
            self.momentum_program_matrix[i,:] = RF_section_parameters_list[i].momentum_program


class RF_section_parameters(object):
    '''
    *Object gathering all the RF parameters for one section (see section
    definition in longitudinal_tracker.Ring_and_RF_section), and pre-processing 
    them in order to be used in the longitudinal_tracker.py module.
    It can be added to another RF_section_parameters object by the 
    Sum_RF_section_parameters object in order to concatenate all the parameters
    for one full ring.*
    '''
    
    def __init__(self, n_turns, n_rf_systems, section_length, 
                 harmonic_number_list, voltage_program_list, phi_offset_list, 
                 momentum_program):
        
        #: *Number of turns for the simulation*
        self.n_turns = n_turns
        
        #: *Length of the section in [m]*
        self.section_length = section_length
        
        #: *Number of RF systems in the section*
        self.n_rf_systems = n_rf_systems
        
        #: | *Momentum program in [eV/c]*
        #: | *The length of the momentum program should be n_turns + 1, check longitudinal_tracker.py for more precisions.*
        #: | *Inputing a single value will assume a constant value for all the simulation.*
        self.momentum_program = input_check(momentum_program, self.n_turns + 1)
        
        #: *Momentum increment (acceleration/deceleration) between two turns,
        #: for one section in [eV/c]*
        self.p_increment = np.diff(self.momentum_program)
        
        #: | *Harmonic number list*
        #: | *The length of the list shoud be equal to n_rf_systems.* 
        self.harmonic_number_list = 0
        
        #: | *Voltage program list in [V]*
        #: | *The length of the list shoud be equal to n_rf_systems.* 
        self.voltage_program_list = 0
        
        #: | *Phase offset list in [rad]*
        #: | *The length of the list shoud be equal to n_rf_systems.* 
        self.phi_offset_list = 0
         
        ### Pre-processing the inputs
        # The input is analysed and structured in order to have lists, which
        # length are matching the number of RF systems considered in this
        # section.
        # For one rf system, single values for h, V_rf, and phi will assume
        # that these values will remain constant for all the simulation.
        # These can be inputed directly as arrays in order to have programs
        # (the length of the arrays will then be checked)
        
        if n_rf_systems == 1:
            self.harmonic_number_list = [harmonic_number_list] 
            self.voltage_program_list = [voltage_program_list] 
            self.phi_offset_list = [phi_offset_list] 
        else:
            if not n_rf_systems == len(harmonic_number_list) == \
                   len(voltage_program_list) == len(phi_offset_list):
                raise RuntimeError('The RF parameters to define \
                                    RF_section_parameters are not homogeneous \
                                    (n_rf_systems is not matching the input)')
            self.harmonic_number_list = harmonic_number_list
            self.voltage_program_list = voltage_program_list 
            self.phi_offset_list = phi_offset_list
        
        for i in range(self.n_rf_systems):
            self.harmonic_number_list[i] = input_check(self.harmonic_number_list[i], n_turns)
            self.voltage_program_list[i] = input_check(self.voltage_program_list[i], n_turns)
            self.phi_offset_list[i] = input_check(self.phi_offset_list[i], n_turns)
            

            
            

            
            
    
