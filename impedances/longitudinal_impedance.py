'''

@author: Hannes Bartosik, Danilo Quartullo, Alexandre Lasheen

'''

from __future__ import division
import numpy as np
from numpy import convolve, interp
from scipy.constants import c, e
from scipy.constants import physical_constants
import time


class Induced_voltage_from_wake(object):
    '''
    Induced voltage derived from sum of several wake fields.
    Apart from further optimizations, these are the important results obtained
    after the benchmarking which are applied to the following code:
    1) in general update_with_interpolation is faster than update_without_interpolation
    2) in general induced_voltage_with_convolv is faster than induced_voltage_with_matrix
    3) if slices.mode == const_charge, you are obliged to use 
        induced_voltage_with_matrix, then you should use update_with_interpolation
    4) if slices.mode == const_space, i.e. you want to calculate slices statistics
        (otherwise slices.mode == const_space_hist is faster), you should use 
        induced_voltage_with_convolv and then update_with_interpolation
    5) if slices.mode == const_space_hist, use induced_voltage_with_convolv and
        update_with_interpolation
    If there is not acceleration then precalc == 'on', except for the const_charge method.
    If you have acceleration and slices.unit == z or theta, then precalc == 'off';
    if slices.unit == tau then precalc == 'on'
    '''
    
    def __init__(self, slices, acceleration, wake_object_sum):       
        '''
        Constructor
        '''
        self.slices = slices
        self.acceleration = acceleration
        self.wake_object_sum = wake_object_sum
        
        if self.slices.mode != 'const_charge':
            
            if self.acceleration == 'off' or self.slices.unit == 'tau':
                self.precalc = 'on'
                translation = self.slices.bins_centers - self.slices.bins_centers[0]
                self.wake_array = self.sum_wakes(translation, self.wake_object_sum)
            else:
                self.precalc = 'off' 
    
    
    def sum_wakes(self, translation, wake_object_sum):
        
        total_wakes = np.zeros(len(translation))
        for wake_object in wake_object_sum:
            total_wakes += wake_object.wake_calc(translation)
            
        return total_wakes
    
    def track(self, bunch):
        
        if self.slices.mode == 'const_charge':
            ind_vol = self.induced_voltage_with_matrix(bunch)
        else:
            ind_vol = self.induced_voltage_with_convolv(bunch)
            
        self.update_with_interpolation(bunch, ind_vol)
        
           
    def induced_voltage_with_matrix(self, bunch):
        
        if self.slices.unit == 'tau':
            dtau_matrix = self.slices.bins_centers - \
                            np.transpose([self.slices.bins_centers])
            self.wake_matrix = self.sum_wakes(dtau_matrix, self.wake_object_sum)
        elif self.slices.unit == 'z':
            dtau_matrix = (np.transpose([self.slices.bins_centers]) - \
                           self.slices.bins_centers) / (bunch.beta_rel * c)
            self.wake_matrix = self.sum_wakes(dtau_matrix, self.wake_object_sum)
        else:
            dtau_matrix = bunch.ring_radius / (bunch.beta_rel * c) * \
            (self.slices.bins_centers - np.transpose([self.slices.bins_centers])) 
            self.wake_matrix = self.sum_wakes(dtau_matrix, self.wake_object_sum)
        
        return - bunch.charge * bunch.intensity / bunch.n_macroparticles * \
                np.dot(self.slices.n_macroparticles, self.wake_matrix)
    
    
    def induced_voltage_with_convolv(self, bunch): 
    
        if self.precalc == 'off':
            if self.slices.unit == 'tau':
                dtau = self.slices.bins_centers - self.slices.bins_centers[0]
                self.wake_array = self.sum_wakes(dtau, self.wake_object_sum)
            elif self.slices.unit == 'z':
                dtau = (self.slices.bins_centers - self.slices.bins_centers[0])/(bunch.beta_rel * c)
                self.wake_array = self.sum_wakes(dtau, self.wake_object_sum)
                reversed_array = self.wake_array[::-1]
                return - bunch.charge * bunch.intensity / bunch.n_macroparticles * \
                    convolve(reversed_array, self.slices.n_macroparticles)[(len(reversed_array) - 1):] 
            elif self.slices.unit == 'theta':
                dtau = (self.slices.bins_centers - self.slices.bins_centers[0]) \
                       * bunch.ring_radius / (bunch.beta_rel * c)
                self.wake_array = self.sum_wakes(dtau, self.wake_object_sum)
                
        return - bunch.charge * bunch.intensity / bunch.n_macroparticles * \
            convolve(self.wake_array, self.slices.n_macroparticles)[0:len(self.wake_array)] 
    
    
    def update_without_interpolation(self, bunch, induced_voltage):
        
        for i in range(0, self.slices.n_slices):
                
                bunch.dE[self.slices.first_index_in_bin[i]:
                  self.slices.first_index_in_bin[i+1]] += induced_voltage[i]
    
    
    def update_with_interpolation(self, bunch, induced_voltage):
        
        if self.slices.unit == 'tau':
            
            temp1 = self.slices.bins_centers[0]
            temp2 = self.slices.bins_centers[-1]
            self.slices.bins_centers[0] = self.slices.edges[0]
            self.slices.bins_centers[-1] = self.slices.edges[-1]
            induced_voltage_interpolated = interp(bunch.tau, self.slices.bins_centers, induced_voltage, 0, 0)
            self.slices.bins_centers[0] = temp1
            self.slices.bins_centers[-1] = temp2
        
        elif self.slices.unit == 'z':
            
            temp1 = self.slices.bins_centers[0]
            temp2 = self.slices.bins_centers[-1]
            self.slices.bins_centers[0] = self.slices.edges[0]
            self.slices.bins_centers[-1] = self.slices.edges[-1]
            induced_voltage_interpolated = interp(bunch.z, self.slices.bins_centers, induced_voltage, 0, 0)
            self.slices.bins_centers[0] = temp1
            self.slices.bins_centers[-1] = temp2
        
        elif self.slices.unit == 'theta':
            
            temp1 = self.slices.bins_centers[0]
            temp2 = self.slices.bins_centers[-1]
            self.slices.bins_centers[0] = self.slices.edges[0]
            self.slices.bins_centers[-1] = self.slices.edges[-1]
            induced_voltage_interpolated = interp(bunch.theta, self.slices.bins_centers, induced_voltage, 0, 0)
            self.slices.bins_centers[0] = temp1
            self.slices.bins_centers[-1] = temp2

        bunch.dE += induced_voltage_interpolated
    
    
class Longit_wake_table(object):
    '''
    classdocs
    '''
    
    def __init__(self, dtau_array, wake_array):       
        '''
        Constructor
        '''
        self.dtau_array = dtau_array
        self.wake_array = wake_array
    
    def wake_calc(self, dtau):
        
        wake = interp(dtau, self.dtau_array - self.dtau_array[0], self.wake_array, left = 0, right = 0)
        
        return wake
    
    
class Longit_wake_resonators(object):
    '''
    
    '''
    def __init__(self, R_shunt, frequency, Q):
        '''
        Constructor
        '''
        self.R_shunt = np.array([R_shunt]).flatten()
        self.frequency = np.array([frequency]).flatten()
        self.Q = np.array([Q]).flatten()
        
    
    def wake_calc(self, dtau):
        
        wake = np.zeros(len(dtau))
        
        for i in range(0, len(self.R_shunt)):
       
            omega = 2 * np.pi * self.frequency[i]
            alpha = omega / (2 * self.Q[i])
            omegabar = np.sqrt(np.abs(omega ** 2 - alpha ** 2))
            
            wake += (np.sign(dtau) + 1) * self.R_shunt[i] * alpha * np.exp(-alpha * dtau) * \
              (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
       
        return wake
        
    
    
  
 
