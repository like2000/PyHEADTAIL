'''

@author: Hannes Bartosik, Danilo Quartullo, Alexandre Lasheen
'''

from __future__ import division
import numpy as np
from numpy import convolve, interp
from scipy.constants import c, e
from scipy.constants import physical_constants
import abc
import time

class Total_wake(object):
    '''
    classdocs
    '''
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def track(self, bunch):
        pass

    
class Long_wake_table(Total_wake):
    '''
    classdocs
    '''
    
    def __init__(self, array_distances, coord_distances, array_wake, slices):       
        '''
        Constructor
        '''
        self.array_distances = array_distances
        self.coord_distances = coord_distances
        self.array_wake = array_wake
        self.slices = slices
        if self.coord_distances != self.slices.unit:
            pass
            
        

    
    


class Long_BB_resonators(Total_wake):
    '''
    Induced voltage derived from resonators.
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
    def __init__(self, R_shunt, frequency, Q, slices, bunch, acceleration):
        '''
        Constructor
        '''
        self.R_shunt = np.array([R_shunt]).flatten()
        self.frequency = np.array([frequency]).flatten()
        self.Q = np.array([Q]).flatten()
        assert(len(self.R_shunt) == len(self.frequency) == len(self.Q))
        self.slices = slices
        self.bunch = bunch
        self.acceleration = acceleration
        
        if self.slices.mode != 'const_charge':
            
            if self.acceleration == 'off' or self.slices.unit == 'tau':
                self.precalc = 'on'
                translation = self.slices.bins_centers - self.slices.bins_centers[0]
                self.wake_array = self.sum_resonators(translation, bunch)
            else:
                self.precalc = 'off' 
        
    
    def sum_resonators(self, dist_betw_centers, bunch):
        return reduce(lambda x,y: x+y, [self.single_resonator(self.R_shunt[i],
         self.frequency[i], self.Q[i], dist_betw_centers, bunch) for i in np.arange(len(self.Q))])

    
    def single_resonator(self, R_shunt, frequency, Q, dtau, bunch):        
        
        omega = 2 * np.pi * frequency
        alpha = omega / (2 * Q)
        omegabar = np.sqrt(np.abs(omega ** 2 - alpha ** 2))
        
        wake = (np.sign(dtau) + 1) * R_shunt * alpha * np.exp(-alpha * dtau) * \
          (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
#         elif self.slices.unit == 'z':
#             dtau = - dist_betw_centers / (bunch.beta_rel * c)
#             wake = (np.sign(dtau) + 1) * R_shunt * alpha * np.exp(-alpha * dtau) * \
#                 (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
#         else:
#             dtau = (bunch.ring_radius * dist_betw_centers) / (bunch.beta_rel * c)
#             wake = (np.sign(dtau) + 1) * R_shunt * alpha * np.exp(-alpha * dtau) * \
#                 (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
        
        return wake
        
    
    def track(self, bunch):
        
        if self.slices.mode == 'const_charge':
            ind_vol = self.induced_voltage_with_matrix(bunch)
        else:
            ind_vol = self.induced_voltage_with_convolv(bunch)
            
        self.update_with_interpolation(bunch, ind_vol)
        
           
    def induced_voltage_with_matrix(self, bunch):
        
        dist_betw_centers = self.slices.bins_centers - np.transpose([self.slices.bins_centers])
        self.wake_matrix = self.sum_resonators(dist_betw_centers, bunch)
        
        return - bunch.charge * bunch.intensity / bunch.n_macroparticles * \
                np.dot(self.slices.n_macroparticles, self.wake_matrix)
    
    
    def induced_voltage_with_convolv(self, bunch): 
    
        if self.precalc == 'off':
            translation = self.slices.bins_centers - self.slices.bins_centers[0]
            self.wake_array = self.sum_resonators(translation, bunch)
        
        return - bunch.charge * bunch.intensity / bunch.n_macroparticles * \
            convolve(self.wake_array, self.slices.n_macroparticles)[0:len(self.wake_array)] 
    
    
    def update_without_interpolation(self, bunch, induced_voltage):
        
        for i in range(0, self.slices.n_slices):
                
                bunch.dE[self.slices.first_index_in_bin[i]:
                  self.slices.first_index_in_bin[i+1]] += induced_voltage[i]
    
    
    def update_with_interpolation(self, bunch, induced_voltage):
        
        
        if self.slices.unit == 'tau':
            
            x = (self.slices.bins_centers[0] + self.slices.bins_centers[1]) / 2
            y = (self.slices.bins_centers[-2] + self.slices.bins_centers[-1]) / 2
            self.slices.bins_centers[0] -= x
            self.slices.bins_centers[-1] += y
            induced_voltage_interpolated = interp(bunch.tau, self.slices.bins_centers, induced_voltage, 0, 0)
            self.slices.bins_centers[0] += x
            self.slices.bins_centers[-1] += y
        
        elif self.slices.unit == 'z':
            
            x = (self.slices.bins_centers[0] + self.slices.bins_centers[1]) / 2
            y = (self.slices.bins_centers[-2] + self.slices.bins_centers[-1]) / 2
            self.slices.bins_centers[0] -= x
            self.slices.bins_centers[-1] += y
            induced_voltage_interpolated = interp(bunch.z, self.slices.bins_centers, induced_voltage, 0, 0)
            self.slices.bins_centers[0] += x
            self.slices.bins_centers[-1] += y
        
        else:
            
            x = (self.slices.bins_centers[0] + self.slices.bins_centers[1]) / 2
            y = (self.slices.bins_centers[-2] + self.slices.bins_centers[-1]) / 2
            self.slices.bins_centers[0] -= x
            self.slices.bins_centers[-1] += y
            induced_voltage_interpolated = interp(bunch.theta, self.slices.bins_centers, induced_voltage, 0, 0)
            self.slices.bins_centers[0] += x
            self.slices.bins_centers[-1] += y

        bunch.dE += induced_voltage_interpolated
  
 
