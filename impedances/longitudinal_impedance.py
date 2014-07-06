'''

@author: Hannes Bartosik, Danilo Quartullo, Alexandre Lasheen

'''

from __future__ import division
import numpy as np
from numpy import convolve, interp
from scipy.constants import c, e
from scipy.constants import physical_constants
import time
from numpy.fft import rfft, irfft, rfftfreq


class Induced_voltage_from_wake(object):
    '''
    Induced voltage derived from the sum of several wake fields.
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
    
    def __init__(self, slices, acceleration, wake_sum):       
        '''
        Constructor
        '''
        self.slices = slices
        self.acceleration = acceleration
        self.wake_sum = wake_sum
        
        if self.slices.mode != 'const_charge':
            
            if self.acceleration == 'off' or self.slices.unit == 'tau':
                self.precalc = 'on'
                translation = self.slices.bins_centers - self.slices.bins_centers[0]
                self.wake_array = self.sum_wakes(translation, self.wake_object_sum)
            else:
                self.precalc = 'off' 
    
    
    def sum_wakes(self, translation, wake_object_sum):
        
        total_wake = np.zeros(len(translation))
        for wake_object in self.wake_sum:
            total_wake += wake_object.wake_calc(translation)
            
        return total_wake
    
    
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
                dtau = (self.slices.bins_centers - self.slices.bins_centers[0])\
                       /(bunch.beta_rel * c)
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
    
    
class Induced_voltage_from_impedance(object):
    '''
    Induced voltage derived from the sum of several impedances.
    '''
    
    def __init__(self, slices, acceleration, impedance_sum, frequency_step, 
                 bunch):       
        '''
        Constructor
        '''
        self.slices = slices
        self.acceleration = acceleration
        self.impedance_sum = impedance_sum
        self.frequency_step = frequency_step
        
        if self.acceleration == 'off' or self.slices.unit == 'tau':
                self.precalc = 'on'
                self.frequency_fft, self.n_sampling_fft = self.frequency_array(slices, bunch)
                self.imped_array = self.sum_impedances(self.frequency_fft, self.imped_sum)    
        else:
            self.precalc = 'off' 
    
    
    def sum_impedances(self, omega, imped_object_sum):
        
        total_impedance = np.zeros(len(omega))
        for imped_object in self.imped_sum:
            total_impedance += imped_object.imped_calc(omega)
            
        return total_impedance
    
    
    def frequency_array(self, slices, bunch):
        
        if self.slices.unit == 'tau':
                    dtau = self.slices.bins_centers[1] - self.slices.bins_centers[0]
        elif self.slices.unit == 'theta':
                    dtau = (self.slices.bins_centers[1] - self.slices.bins_centers[0]) \
                       * bunch.ring_radius / (bunch.beta_rel * c)
        elif self.slices.unit == 'z':
                    dtau = (self.slices.bins_centers[1] - self.slices.bins_centers[0])\
                       /(bunch.beta_rel * c)
        power = np.floor(np.log2(1 / (self.frequency_step * dtau))) + 2
        
        return rfftfreq(2 ** power, dtau), 2 ** power
    
    
    def track(self, bunch):
        
        if self.precalc == 'off':
            self.frequency_fft, self.n_sampling_fft = self.frequency_array(self.slices, bunch)
            self.imped_array = self.sum_impedances(self.frequency_fft, self.imped_sum)  
            
        spectrum = bunch.spectrum(self.n_sampling_fft, self.slices)
         
        ind_vol = - bunch.charge * bunch.intensity / bunch.n_macroparticles \
                    * irfft(self.impedance_array * spectrum) * \
                     self.frequency_fft[1] * 2*(len(self.frequency_fft)-1)
        ind_vol = ind_vol[0:len(self.frequency_fft)]
        self.update_with_interpolation(bunch, ind_vol)
    

def update_without_interpolation(self, bunch, induced_voltage):
    
    for i in range(0, self.slices.n_slices):
            
            bunch.dE[self.slices.first_index_in_bin[i]:
              self.slices.first_index_in_bin[i+1]] += induced_voltage[i]
    
    
def update_with_interpolation(self, bunch, induced_voltage):
    
    temp1 = self.slices.bins_centers[0]
    temp2 = self.slices.bins_centers[-1]
    self.slices.bins_centers[0] = self.slices.edges[0]
    self.slices.bins_centers[-1] = self.slices.edges[-1]
    
    if self.slices.unit == 'tau':
        
        induced_voltage_interpolated = interp(bunch.tau, 
                            self.slices.bins_centers, induced_voltage, 0, 0)
        
    elif self.slices.unit == 'z':
        
        induced_voltage_interpolated = interp(bunch.z, 
                            self.slices.bins_centers, induced_voltage, 0, 0)
        
    elif self.slices.unit == 'theta':
        
        induced_voltage_interpolated = interp(bunch.theta, 
                            self.slices.bins_centers, induced_voltage, 0, 0)
        
    self.slices.bins_centers[0] = temp1
    self.slices.bins_centers[-1] = temp2
    bunch.dE += induced_voltage_interpolated
    
    
class Longitudinal_table(object):
    '''
    classdocs
    '''
    
    def __init__(self, a, b, c = None):       
        '''
        Constructor
        '''
        if c == None:
            self.dtau_array = a
            self.wake_array = b
        else:
            self.omega_array = a
            self.Re_Z_array = b
            self.Im_Z_array = c
    
    
    def wake_calc(self, dtau):
        
        wake = interp(dtau, self.dtau_array - self.dtau_array[0], 
                      self.wake_array, left = 0, right = 0)
        
        return wake
    
    
    def imped_calc(self, omega):
        
        Re_Z = interp(omega, self.omega_array, self.Re_Z_array, right = 0)
        Im_Z = interp(omega, self.omega_array, self.Im_Z_array, right = 0)
        
        return Re_Z + 1j * Im_Z
    
    
    
class Longitudinal_resonators(object):
    '''
    
    '''
    def __init__(self, R_S, omega_R, Q):
        '''
        Constructor
        '''
        self.R_S = np.array([R_S]).flatten()
        self.omega_R = np.array([omega_R]).flatten()
        self.Q = np.array([Q]).flatten()
        self.n_resonators = len(self.R_S)
        
    
    def wake_calc(self, dtau):
        
        wake = np.zeros(len(dtau))
        
        for i in range(0, self.n_resonators):
       
            alpha = self.omega_R[i] / (2 * self.Q[i])
            omega_bar = np.sqrt(self.omega_R ** 2 - alpha ** 2)
            
            wake += (np.sign(dtau) + 1) * self.R_S[i] * alpha * np.exp(-alpha * 
                    dtau) * (np.cos(omega_bar * dtau) - alpha / omega_bar * 
                    np.sin(omega_bar * dtau))
       
        return wake
    
    
    def imped_calc(self, omega):
        
        impedance = np.zeros(len(omega)) + 0j
        
        for i in range(0, len(self.R_S)):
            
            impedance[1:] +=  self.R_S[i] / (1 + 1j * self.Q[i] * \
                    (self.omega_R[i] / omega[1:] - omega[1:] / self.omega_R[i]))
            
        return impedance
 

class Longitudinal_travelling_waves(object):
    '''
    
    '''
    def __init__(self, R_S, omega_R, a_factor):
        '''
        Constructor
        '''
        self.R_S = np.array([R_S]).flatten()
        self.omega_R = np.array([omega_R]).flatten()
        self.a_factor = np.array([a_factor]).flatten()
        self.n_twc = len(self.R_S)
        
    
    def wake_calc(self, dtau):
        
        wake = np.zeros(len(dtau))
        
        for i in range(0, self.n_twc):
       
            a_tilde = self.a_factor[i] / (2 * np.pi)
            indexes = np.where(dtau <= a_tilde)
            wake[indexes] += (np.sign(dtau[indexes]) + 1) * 2 * self.R_S[i] \
                / a_tilde * (1 - dtau[indexes] / a_tilde) * np.cos(2 * np.pi * 
                self.omega_R[i] * dtau[indexes])
                 
        return wake
    
    
    def imped_calc(self, omega):
        
        impedance = np.zeros(len(omega)) + 0j
        
        for i in range(0, self.n_twc):
            
            impedance +=  self.R_S[i] * ((np.sin(self.a_factor[i] / 2 * 
                (omega - self.omega_R[i])) / (self.a_factor[i] / 2 * (omega - 
                self.omega_R[i])))**2 - 2j*(self.a_factor[i] * (omega - 
                self.omega_R[i]) - np.sin(self.a_factor[i] * (omega - 
                self.omega_R[i]))) / (self.a_factor[i] * (omega - 
                self.omega_R[i]))**2) + self.R_S[i] * ((np.sin(self.a_factor[i] 
                / 2 * (omega + self.omega_R[i])) / (self.a_factor[i] / 2 * (
                omega + self.omega_R[i])))**2 - 2j*(self.a_factor[i] * (omega 
                + self.omega_R[i]) - np.sin(self.a_factor[i] * (omega + 
                self.omega_R[i]))) / (self.a_factor[i] * (omega + 
                self.omega_R[i]))**2)
            
        return impedance       
    
    
  
 
