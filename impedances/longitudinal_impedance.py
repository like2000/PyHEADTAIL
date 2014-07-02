'''

@author: Hannes Bartosik, Danilo Quartullo, Alexandre Lasheen
'''

from __future__ import division
import numpy as np
from numpy import convolve, interp
from scipy.constants import c, e
from scipy.constants import physical_constants
import abc


class Wakefields(object):
    '''
    classdocs
    '''
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def track(self, bunch):
        pass

    
class long_wake_table(Wakefields):
    '''
    classdocs
    '''
    def __init__(self):       
        '''
        Constructor
        '''
        self.wake_table = {}

    
    @classmethod
    def from_ASCII(cls, wake_file, keys):
        self = cls()
        table = np.loadtxt(wake_file, delimiter="\t")
        self.wake_table = dict(zip(keys, np.array(zip(*table))))
        self.unit_conversion()
        return self

    def unit_conversion(self):
        longitudinal_wakefield_keys = ['longitudinal']
        self.wake_field_keys = []
        print 'Converting wake table to correct units ... '
        self.wake_table['time'] *= 1e-9 # unit convention [ns]
        print '\t converted time from [ns] to [s]'        
        
        for wake in longitudinal_wakefield_keys:
            try: 
                self.wake_table[wake] *= - 1.e12 # unit convention [V/pC] and sign convention !!
                print '\t converted "' + wake + '" wake from [V/pC/mm] to [V/C/m]'
                self.wake_field_keys += [wake]
            except:
                print '\t "' + wake + '" wake not provided'

    def wake_longitudinal(self, bunch, z):
        time = np.array(self.wake_table['time'])
        wake = np.array(self.wake_table['longitudinal'])
        wake_interpolated = np.interp(- z / c / bunch.beta, time, wake, left=0, right=0)
        if time[0] < 0:
            return wake_interpolated
        elif time[0] == 0:
            # beam loading theorem: half value of wake at z=0; 
            return (np.sign(-z) + 1) / 2 * wake_interpolated
    
    
    def track(self, bunch):
        
        if 'longitudinal' in self.wake_field_keys:
            self.longitudinal_wakefield_kicks(bunch)


class Long_BB_resonators(Wakefields):
    '''
    classdocs
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
        if self.slices.mode == 'const_charge':
            self.mode = "matrix_no_precalc"
        elif self.acceleration == 'off':
            self.mode = 'matrix_with_precalc'
            dist_betw_centers = slices.bins_centers - np.transpose([slices.bins_centers])
            self.wake_matrix = self.wake_longitudinal(dist_betw_centers, self.bunch)
        elif self.slices.unit == 'tau':
            self.mode = 'matrix_with_precalc'
            dist_betw_centers = slices.bins_centers - np.transpose([slices.bins_centers])
            self.wake_matrix = self.wake_longitudinal(dist_betw_centers, self.bunch)
        else:
            self.mode = "matrix_no_precalc"
        
    
    def wake_longitudinal(self, dist_betw_centers, bunch):
        return reduce(lambda x,y: x+y, [self.wake_BB_resonator(self.R_shunt[i],
         self.frequency[i], self.Q[i], dist_betw_centers, bunch) for i in np.arange(len(self.Q))])

    
    def wake_BB_resonator(self, R_shunt, frequency, Q, dist_betw_centers, bunch):        
        
        omega = 2 * np.pi * frequency
        alpha = omega / (2 * Q)
        omegabar = np.sqrt(np.abs(omega ** 2 - alpha ** 2))
        
        if self.slices.unit == 'tau':
            dtau = dist_betw_centers
            wake = (np.sign(dtau) + 1) * R_shunt * alpha * np.exp(-alpha * dtau) * \
                (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
        elif self.slices.unit == 'z':
            dtau = - dist_betw_centers / (bunch.beta_rel * c)
            wake = (np.sign(dtau) + 1) * R_shunt * alpha * np.exp(-alpha * dtau) * \
                (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
        else:
            dtau = (bunch.ring_radius * dist_betw_centers) / (bunch.beta_rel * c)
            wake = (np.sign(dtau) + 1) * R_shunt * alpha * np.exp(-alpha * dtau) * \
                (np.cos(omegabar * dtau) - alpha / omegabar * np.sin(omegabar * dtau))
        
        return wake
        
    
    def track(self, bunch):
        
        ind = self.induced_voltage_with_matrix(bunch)
        self.update_without_interpolation(bunch, ind)
            
    
    def induced_voltage_with_matrix(self, bunch):
        
        if self.mode == "no_precalc":
            dist_betw_centers = self.slices.bins_centers - np.transpose([self.slices.bins_centers])
            self.wake_matrix = self.wake_longitudinal(dist_betw_centers, bunch)
        
        return - bunch.charge * np.dot(
                                self.slices.n_macroparticles, self.wake_matrix)
    
    
    def induced_voltage_with_convolv(self, bunch): 
    
        if self.mode == "no_precalc":
            translation = self.slices.bins_centers - self.slices.bins_centers[0]
            wake_array = self.wake_longitudinal(translation, bunch)
        
        return - bunch.charge * convolve(wake_array, self.slices.n_macroparticles)[0:len(self.wake_array)]
    
    
    def update_without_interpolation(self, bunch, induced_voltage):
        
        for i in range(0, self.slices.n_slices):
                
                bunch.dE[self.slices.first_index_in_bin[i]:
                  self.slices.first_index_in_bin[i+1]] += induced_voltage[i]
    
    
    def update_with_interpolation(self, bunch, induced_voltage):
        
        if self.slices.unit == 'tau':
        
            induced_voltage_interpolated = interp(bunch.tau, self.slices.bins_centers, induced_voltage)
        
        elif self.slices.unit == 'z':
            
            induced_voltage_interpolated = interp(bunch.z, self.slices.bins_centers, induced_voltage)
        
        else:
            
            induced_voltage_interpolated = interp(bunch.theta, self.slices.bins_centers, induced_voltage)

        bunch.dE += induced_voltage_interpolated
  
 
