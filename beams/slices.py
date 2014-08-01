'''
**Module to compute longitudinal beam slicing**

:Authors: **Hannes Bartosik**, **Kevin Li**, **Michael Schenk**, **Danilo Quartullo**, **Alexandre Lasheen**
'''

from __future__ import division
import numpy as np
from random import sample
from scipy.constants import c
from numpy.fft import rfft
from scipy import ndimage


class Slices(object):
    '''
    *Slices class that controls longitudinal discretisation of a Beam. This
    include the Beam profiling (including computation of Beam spectrum,
    derivative, and profile fitting) and the computation of statistics per
    slice.*
    '''

    def __init__(self, Beam, n_slices, n_sigma = None, cut_left = None, 
                 cut_right = None, coord = "theta", mode = 'const_space',
                 statistics_option = 'off', fit_option = 'off'):
        
        #: *Copy (reference) of the beam to be sliced (from Beam)*
        self.Beam = Beam
        
        #: *Number of slices*
        self.n_slices = n_slices
        
        #: | *Slicing computation mode*
        #: | *The options are: 'const_space' (default), 'const_charge'.*
        self.mode = mode
        
        #: *Left edge of the slicing (is an optionnal input, in case you use
        #: the 'const_space' mode, a default value will be set if no value is 
        #: given).*
        self.cut_left = cut_left
        
        #: *Right edge of the slicing (is an optionnal input, in case you use
        #: the 'const_space' mode, a default value will be set if no value is 
        #: given).*
        self.cut_right = cut_right
        
        #: *Optionnal input parameters, corresponding to the number of*
        #: :math:`\sigma_{RMS}` *of the Beam to slice (this will overwrite
        #: any input of cut_left and cut_right).*
        self.n_sigma = n_sigma
        
        #: | *Type of coordinates in which the slicing is done.*
        #: | *The options are: 'theta' (default), 'tau', 'z' (default if constant_charge).*
        self.coord = coord
        if self.mode is 'const_charge':
            self.coord = 'z'
        
        #: *Number of macroparticles per slice (~profile).*
        self.n_macroparticles = np.empty(n_slices)
        
        #: *Edges positions of the slicing*
        self.edges = np.empty(n_slices + 1)
        
        #: *Center of the bins*
        self.bins_centers = np.empty(n_slices)
        
        # Pre-processing the slicing (for const_space only)
        if self.cut_left is None and self.cut_right is None:
            self.set_longitudinal_cuts()
        if self.mode is 'const_space' or self.mode is 'const_space_hist':
            self.edges = np.linspace(self.cut_left, self.cut_right, self.n_slices + 1)
            self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
        
        #: *Compute statistics option allows to compute mean_theta, mean_dE, 
        #: sigma_theta and sigma_dE properties each turn.*
        self.statistics_option = statistics_option
        
        if self.statistics_option is 'on':
            #: *Average theta position of the particles in each slice (needs 
            #: the compute_statistics_option to be 'on').*
            self.mean_theta = np.empty(n_slices)
            #: *Average dE position of the particles in each slice (needs 
            #: the compute_statistics_option to be 'on').*
            self.mean_dE = np.empty(n_slices)
            #: *RMS theta position of the particles in each slice (needs 
            #: the compute_statistics_option to be 'on').*
            self.sigma_theta = np.empty(n_slices)
            #: *RMS dE position of the particles in each slice (needs 
            #: the compute_statistics_option to be 'on').*
            self.sigma_dE = np.empty(n_slices)
            
        #: *Fit option allows to fit the Beam profile, with the options
        #: 'off' (default), 'gaussian'.*
        self.fit_option = fit_option
            
        if self.fit_option is 'gaussian':
            #: *Beam length with a gaussian fit (needs fit_option to be 
            #: 'gaussian' defined as* :math:`\tau_{gauss} = 4\sigma`)
            self.bl_gauss = 0
            #: *Beam position with a gaussian fit (needs fit_option to be 
            #: 'gaussian')*
            self.bp_gauss = 0
                

    def set_longitudinal_cuts(self):
        '''
        *Method to set the self.cut_left and self.cut_right properties. This is
        done as a pre-processing if the mode is set to 'const_space', for
        'const_charge' this is calculated each turn.*
        
        *The frame is defined by :math:`n\sigma_{RMS}` or manually by the user.
        If not, a default frame consisting of taking the whole bunch +5% of the 
        maximum distance between two particles in the bunch will be taken.*
        '''
        
        if self.n_sigma is None:
            
            self.sort_particles()
            
            if self.coord == "theta":
                self.cut_left = self.Beam.theta[0] - 0.05*(self.Beam.theta[-1] - self.Beam.theta[0])
                self.cut_right = self.Beam.theta[-1] + 0.05*(self.Beam.theta[-1] - self.Beam.theta[0])
            elif self.coord == "z":
                self.cut_left = self.Beam.z[0] - 0.05*(self.Beam.z[-1] - self.Beam.z[0])
                self.cut_right = self.Beam.z[-1] + 0.05*(self.Beam.z[-1] - self.Beam.z[0])
            else:
                self.cut_left = self.Beam.tau[0] - 0.05*(self.Beam.tau[-1] - self.Beam.tau[0])
                self.cut_right = self.Beam.tau[-1] + 0.05*(self.Beam.tau[-1] - self.Beam.tau[0])
                
        else:
            if self.coord == "theta":
                mean_theta = np.mean(self.Beam.theta)
                sigma_theta = np.std(self.Beam.theta)
                self.cut_left = mean_theta - self.n_sigma * sigma_theta / 2
                self.cut_right = mean_theta + self.n_sigma * sigma_theta / 2
            elif self.coord == "z":
                mean_z = np.mean(self.Beam.z)
                sigma_z = np.std(self.Beam.z)
                self.cut_left = mean_z - self.n_sigma * sigma_z / 2
                self.cut_right = mean_z + self.n_sigma * sigma_z / 2
            else:
                mean_tau = np.mean(self.Beam.tau)
                sigma_tau = np.std(self.Beam.tau)
                self.cut_left = mean_tau - self.n_sigma * sigma_tau / 2
                self.cut_right = mean_tau + self.n_sigma * sigma_tau / 2


    def slice_constant_space(self):
        '''
        *Constant space slicing. This method consist in slicing a fixed frame
        (which length is determined in the beginning of the simulation) with
        bins of constant size. Each turn, the particles are sorted with respect
        to their longitudinal position and counted in each bin. This allows
        also to calculate the statistics of the particles for each bin (if 
        statistics_option is 'on') and fit the profile (e.g. Gaussian).*
        
        *Be careful that because the frame is not changing, a bunch with
        increasing bunch length might not be sliced properly as part of it
        might be out of the frame.*
        '''

        self.sort_particles()

        if self.coord == 'z':
            first_index_in_bin = np.searchsorted(self.Beam.z, self.edges)
        elif self.coord == 'theta':
            first_index_in_bin = np.searchsorted(self.Beam.theta, self.edges)
        else:
            first_index_in_bin = np.searchsorted(self.Beam.tau, self.edges)
            
        self.n_macroparticles = np.diff(first_index_in_bin)
        

    def slice_constant_space_histogram(self):
        '''
        *Appears to be slower than the constant_space version, to be removed...*
        '''
        
        if self.coord == 'theta':
            self.n_macroparticles = np.histogram(self.Beam.theta, self.edges)[0]
        elif self.coord == 'z':
            self.n_macroparticles = np.histogram(self.Beam.z, self.edges)[0]
        else:
            self.n_macroparticles = np.histogram(self.Beam.tau, self.edges)[0]
 
        
    def slice_constant_charge(self):
        '''
        *Constant charge slicing. This method consist in slicing with varying
        bin sizes that adapts in order to have the same number of particles
        in each bin*
        
        *Must be updated in order to take into account potential losses (in order
        for the frame size not to diverge) and for different types of coordinates
        (only z coordinate for the moment)*
        '''
        
        self.coord = 'z'
        self.cut_left = None
        self.cut_right = None
        self.n_sigma = None
        
        self.set_longitudinal_cuts()
        
        # 1. n_macroparticles - distribute macroparticles uniformly along slices.
        # Must be integer. Distribute remaining particles randomly among slices with indices 'ix'.
        n_cut_left = 0 # number of particles cut left, to be adapted for losses
        n_cut_right = 0 # number of particles cut right, to be adapted for losses
          
        q0 = self.Beam.n_macroparticles - n_cut_right - n_cut_left
         
        ix = sample(range(self.n_slices), q0 % self.n_slices)
        self.n_macroparticles = (q0 // self.n_slices) * np.ones(self.n_slices)
        self.n_macroparticles[ix] += 1
        
        # 2. edges
        # Get indices of the particles defining the bin edges
        n_macroparticles_all = np.hstack((n_cut_left, self.n_macroparticles, n_cut_right))
        first_index_in_bin = np.cumsum(n_macroparticles_all)
        first_particle_index_in_slice = first_index_in_bin[:-1]
        first_particle_index_in_slice = (first_particle_index_in_slice).astype(int)
         
        self.edges[1:-1] = (self.Beam.z[(first_particle_index_in_slice - 1)[1:-1]] + self.Beam.z[first_particle_index_in_slice[1:-1]]) / 2
        self.edges[0], self.edges[-1] = self.cut_left, self.cut_right
        self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
    

    def track(self, Beam):

        if self.mode == 'const_charge':
            self.slice_constant_charge()
        elif self.mode == 'const_space':
            self.slice_constant_space()
        elif self.mode == 'const_space_hist':
            self.slice_constant_space_histogram()
        else:
            raise RuntimeError('Choose one proper slicing mode!')


    def sort_particles(self):
        '''
        | *Sort the particles with respect to their longitudinal position.*
        | *'theta' and 'tau' type of coordinates are treated differently that 'z' as head and tail are reversed.*
        '''
               
        if self.coord == 'theta' or self.coord == 'tau':  
            argsorted = np.argsort(self.Beam.theta)
            
        elif self.coord == 'z':
            argsorted = np.argsort(self.Beam.z)
            
        self.Beam.theta = self.Beam.theta.take(argsorted)
        self.Beam.dE = self.Beam.dE.take(argsorted)
        self.Beam.id = self.Beam.id.take(argsorted)
        
        
    def gaussian_fit(self):
        '''
        *Gaussian fit of the profile, in order to get the bunch length and
        position.*
        '''
        
        pass
    
#         ##### Gaussian fit to theta-profile
#         if gaussian_fit == "On":
#             
#             if slices == None:
#                 warnings.filterwarnings("once")                
#                 warnings.warn("WARNING: The Gaussian bunch length fit cannot be calculated without slices!")
#             else:
#                 try:
#                     if slices.coord == "theta":
#                         p0 = [max(slices.n_macroparticles), self.mean_theta, self.sigma_theta]                
#                         pfit = curve_fit(gauss, slices.bins_centers, 
#                                          slices.n_macroparticles, p0)[0]
#                     elif slices.coord == "tau":
#                         p0 = [max(slices.n_macroparticles), self.mean_tau, self.sigma_tau]                
#                         pfit = curve_fit(gauss, slices.bins_centers, 
#                                          slices.n_macroparticles, p0)[0]  
#                     elif slices.coord == "z":
#                         p0 = [max(slices.n_macroparticles), self.mean_z, self.sigma_z]                
#                         pfit = curve_fit(gauss, slices.bins_centers, 
#                                          slices.n_macroparticles, p0)[0]                                    
#                     self.bl_gauss = 4 * abs(pfit[2]) 
#                 except:
#                     self.bl_gauss = 0 

    
    
    def beam_spectrum(self, n_sampling_fft):
         
        spectrum = rfft(self.n_macroparticles, n_sampling_fft)
             
        return spectrum
     
     
    def beam_profile_derivative(self, mode):      
        ''' 
        *The input is one of the two available methods for differentiating
        a function. The two outputs are the coordinate step and the discrete
        derivative of the Beam profile respectively.*
        '''
         
        dist_centers = self.bins_centers[1] - self.bins_centers[0]
         
        if mode == 1:
            x = self.bins_centers
            derivative = ndimage.gaussian_filter1d(self.n_macroparticles, sigma=1, order=1, mode='wrap') / dist_centers
        if mode == 2:
            x = self.bins_centers
            derivative = np.gradient(self.n_macroparticles, dist_centers)
             
        return x, derivative
     
     
    @property    
    def mean_z(self):
        '''*Average z position of the particles in each slice (needs 
        the compute_statistics_option to be 'on', the head and tail are
        reversed compared to theta and an tau coordinates).*'''
        return - self.mean_theta * self.Beam.ring_radius 
    @mean_z.setter
    def mean_z(self, value):
        self.mean_theta = - value / self.Beam.ring_radius 
     
    @property    
    def mean_tau(self):
        '''*Average tau position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.mean_theta * self.Beam.ring_radius / (self.Beam.beta_rel * c)
    @mean_tau.setter
    def mean_tau(self, value):
        self.mean_theta = value * self.Beam.beta_rel * c / self.Beam.ring_radius
 
    @property
    def mean_delta(self):
        '''*Average delta position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.mean_dE / (self.Beam.beta_rel**2 * self.Beam.energy)
    @mean_delta.setter
    def mean_delta(self, value):
        self.mean_dE = value * self.Beam.beta_rel**2 * self.Beam.energy
     
    @property    
    def sigma_z(self):
        return - self.sigma_theta * self.Beam.ring_radius 
    @sigma_z.setter
    def sigma_z(self, value):
        self.sigma_theta = - value / self.Beam.ring_radius 
     
    @property
    def sigma_tau(self):
        return self.sigma_theta * self.Beam.ring_radius / (self.Beam.beta_rel * c)
    @sigma_tau.setter
    def sigma_tau(self, value):
        self.sigma_theta = value * self.Beam.beta_rel * c / self.Beam.ring_radius
     
    @property
    def sigma_delta(self):
        return self.sigma_dE / (self.Beam.beta_rel**2 * self.Beam.energy)
    @sigma_delta.setter
    def sigma_delta(self, value):
        self.sigma_dE = value * self.Beam.beta_rel**2 * self.Beam.energy


def gauss(x, *p):
    '''
    *Defined as:*
    
    .. math:: A \, e^{\\frac{\\left(x-x_0\\right)^2}{2\\sigma_x^2}}
    
    '''
    
    A, x0, sx = p
    return A*np.exp(-(x-x0)**2/2./sx**2) 
    
        
