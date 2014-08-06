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
from scipy.optimize import curve_fit


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
        
        #: | *Type of coordinates in which the cuts are given.*
        #: | *The options are: 'theta' (default), 'tau', 'z'.*
        self.coord = coord
#         if self.mode is 'const_charge':
#             self.coord = 'z'
        
        #: *Number of macroparticles per slice (~profile).*
        self.n_macroparticles = np.empty(n_slices)
        
        #: *Edges positions of the slicing*
        self.edges = np.empty(n_slices + 1)
        
        #: *Center of the bins*
        self.bins_centers = np.empty(n_slices)
        
        # Pre-processing the slicing (for const_space only)
        if self.cut_left is None and self.cut_right is None:
            self.set_longitudinal_cuts()
        if self.coord == "theta":
            self.cut_left = cut_left
            self.cut_right = cut_right
        elif self.coord == "z":
            self.cut_left = cut_left / self.Beam.ring_radius 
            self.cut_right = cut_right / self.Beam.ring_radius 
        elif self.coord == 'tau':
            self.cut_left = cut_left / (self.Beam.ring_radius / (self.Beam.beta_rel * c))
            self.cut_right = cut_right / (self.Beam.ring_radius / (self.Beam.beta_rel * c))
        if self.mode is not 'const_charge':
            self.edges = np.linspace(self.cut_left, self.cut_right, self.n_slices + 1)
            self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
        
        #: *Compute statistics option allows to compute mean_theta, mean_dE, 
        #: sigma_theta and sigma_dE properties each turn.*
        self.statistics_option = statistics_option
        
        if self.statistics_option is 'on' and self.mode is 'const_space_hist':
            raise RuntimeError('Slice compute statistics does not work with \
                                the const_space_hist mode !')
        elif self.statistics_option is 'on':
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
            #: *RMS dE position of the particles in each slice (needs 
            #: the compute_statistics_option to be 'on').*
            self.eps_rms_l = np.empty(n_slices)
            
        #: *Fit option allows to fit the Beam profile, with the options
        #: 'off' (default), 'gaussian'.*
        self.fit_option = fit_option
        
        # Pre-processing Beam longitudinal statistics
        self.Beam.longit_statistics()
        
        if self.mode is not 'const_charge' and self.fit_option is 'gaussian':    
            #: *Beam length with a gaussian fit (needs fit_option to be 
            #: 'gaussian' defined as* :math:`\tau_{gauss} = 4\sigma`)
            self.bl_gauss = 0
            #: *Beam position with a gaussian fit (needs fit_option to be 
            #: 'gaussian')*
            self.bp_gauss = 0
            #: *Gaussian parameters list obtained from fit*
            self.pfit_gauss = 0
                    
        # Use of track in order to pre-process the slicing at injection
        self.track(self.Beam)
        

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
            
#             if self.coord == "theta":
#                 self.cut_left = self.Beam.theta[0] - 0.05*(self.Beam.theta[-1] - self.Beam.theta[0])
#                 self.cut_right = self.Beam.theta[-1] + 0.05*(self.Beam.theta[-1] - self.Beam.theta[0])
#             elif self.coord == "z":
#                 self.cut_left = self.Beam.z[0] - 0.05*(self.Beam.z[-1] - self.Beam.z[0])
#                 self.cut_right = self.Beam.z[-1] + 0.05*(self.Beam.z[-1] - self.Beam.z[0])
#             else:
#                 self.cut_left = self.Beam.tau[0] - 0.05*(self.Beam.tau[-1] - self.Beam.tau[0])
#                 self.cut_right = self.Beam.tau[-1] + 0.05*(self.Beam.tau[-1] - self.Beam.tau[0])

            self.cut_left = self.Beam.theta[0] - 0.05*(self.Beam.theta[-1] - self.Beam.theta[0])
            self.cut_right = self.Beam.theta[-1] + 0.05*(self.Beam.theta[-1] - self.Beam.theta[0])
                
        else:
#             if self.coord == "theta":
#                 mean_theta = np.mean(self.Beam.theta)
#                 sigma_theta = np.std(self.Beam.theta)
#                 self.cut_left = mean_theta - self.n_sigma * sigma_theta / 2
#                 self.cut_right = mean_theta + self.n_sigma * sigma_theta / 2
#             elif self.coord == "z":
#                 mean_z = np.mean(self.Beam.z)
#                 sigma_z = np.std(self.Beam.z)
#                 self.cut_left = mean_z - self.n_sigma * sigma_z / 2
#                 self.cut_right = mean_z + self.n_sigma * sigma_z / 2
#             else:
#                 mean_tau = np.mean(self.Beam.tau)
#                 sigma_tau = np.std(self.Beam.tau)
#                 self.cut_left = mean_tau - self.n_sigma * sigma_tau / 2
#                 self.cut_right = mean_tau + self.n_sigma * sigma_tau / 2
                
            mean_theta = np.mean(self.Beam.theta)
            sigma_theta = np.std(self.Beam.theta)
            self.cut_left = mean_theta - self.n_sigma * sigma_theta / 2
            self.cut_right = mean_theta + self.n_sigma * sigma_theta / 2


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

#         if self.coord == 'z':
#             first_index_in_bin = np.searchsorted(self.Beam.z, self.edges)
#         elif self.coord == 'theta':
#             first_index_in_bin = np.searchsorted(self.Beam.theta, self.edges)
#         elif self.coord == 'tau':
#             first_index_in_bin = np.searchsorted(self.Beam.tau, self.edges)
            
        first_index_in_bin = np.searchsorted(self.Beam.theta, self.edges)
            
        self.n_macroparticles = np.diff(first_index_in_bin)
        

    def slice_constant_space_histogram(self):
        '''
        *Constant space slicing with the built-in numpy histogram function,
        with a constant frame. This gives the same profile as the 
        slice_constant_space method, but no compute statistics possibilities
        (the index of the particles is needed).*
        
        *This method is faster than the classic slice_constant_space method 
        for high number of particles (~1e6).*
        '''
        
#         if self.coord == 'theta':
#             self.n_macroparticles = np.histogram(self.Beam.theta, self.edges)[0]
#         elif self.coord == 'z':
#             self.n_macroparticles = np.histogram(self.Beam.z, self.edges)[0]
#         else:
#             self.n_macroparticles = np.histogram(self.Beam.tau, self.edges)[0]
            
        self.n_macroparticles = np.histogram(self.Beam.theta, self.edges)[0]
 
        
#     def slice_constant_charge(self):
#         '''
#         *Constant charge slicing. This method consist in slicing with varying
#         bin sizes that adapts in order to have the same number of particles
#         in each bin*
#         
#         *Must be updated in order to take into account potential losses (in order
#         for the frame size not to diverge). The slicing is done in the 'z'
#         coordinate, but should be optimized in order to be done in 'theta'
#         and the results then converted in 'z'*
#         '''
#         
#         self.coord = 'z'
#         self.cut_left = None
#         self.cut_right = None
#         self.n_sigma = None
#         
#         self.set_longitudinal_cuts()
#         
#         # 1. n_macroparticles - distribute macroparticles uniformly along slices.
#         # Must be integer. Distribute remaining particles randomly among slices with indices 'ix'.
#         n_cut_left = 0 # number of particles cut left, to be adapted for losses
#         n_cut_right = 0 # number of particles cut right, to be adapted for losses
#           
#         q0 = self.Beam.n_macroparticles - n_cut_right - n_cut_left
#          
#         ix = sample(range(self.n_slices), q0 % self.n_slices)
#         self.n_macroparticles = (q0 // self.n_slices) * np.ones(self.n_slices)
#         self.n_macroparticles[ix] += 1
#         
#         # 2. edges
#         # Get indices of the particles defining the bin edges
#         n_macroparticles_all = np.hstack((n_cut_left, self.n_macroparticles, n_cut_right))
#         first_index_in_bin = np.cumsum(n_macroparticles_all)
#         first_particle_index_in_slice = first_index_in_bin[:-1]
#         first_particle_index_in_slice = (first_particle_index_in_slice).astype(int)
#          
#         self.edges[1:-1] = (self.Beam.z[(first_particle_index_in_slice - 1)[1:-1]] + 
#                             self.Beam.z[first_particle_index_in_slice[1:-1]]) / 2
#         self.edges[0], self.edges[-1] = self.cut_left, self.cut_right
#         self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
    

    def track(self, Beam):
        '''
        *Track method in order to update the slicing along with the tracker.
        This will update the beam properties (bunch length obtained from the
        fit, etc.).*
        '''

        if self.mode == 'const_charge':
            raise RuntimeError('const_charge still needs some corrections, sorry...')
#             self.slice_constant_charge()
        elif self.mode == 'const_space':
            self.slice_constant_space()
        elif self.mode == 'const_space_hist':
            self.slice_constant_space_histogram()
        else:
            raise RuntimeError('Choose one proper slicing mode!')
        
        if self.fit_option is 'gaussian':
            self.gaussian_fit()
            self.Beam.bl_gauss = self.bl_gauss
            self.Beam.bp_gauss = self.bp_gauss
            
        if self.statistics_option is 'on':
            self.compute_statistics()


    def sort_particles(self):
        '''
        | *Sort the particles with respect to their longitudinal position.*
        | *'theta' and 'tau' type of coordinates are treated differently that 'z' as head and tail are reversed.*
        '''
               
#         if self.coord == 'theta' or self.coord == 'tau':  
#             argsorted = np.argsort(self.Beam.theta)
#             
#         elif self.coord == 'z':
#             argsorted = np.argsort(self.Beam.z)
            
        argsorted = np.argsort(self.Beam.theta)
            
        self.Beam.theta = self.Beam.theta.take(argsorted)
        self.Beam.dE = self.Beam.dE.take(argsorted)
        self.Beam.id = self.Beam.id.take(argsorted)
        
        
    def gaussian_fit(self):
        '''
        *Gaussian fit of the profile, in order to get the bunch length and
        position.*
        '''
            
        if self.bl_gauss is 0 and self.bp_gauss is 0:
#             if self.coord is 'theta':
#                 p0 = [max(self.n_macroparticles), self.Beam.mean_theta, self.Beam.sigma_theta]
#             elif self.coord is 'tau':
#                 p0 = [max(self.n_macroparticles), self.Beam.mean_tau, self.Beam.sigma_tau]
#             elif self.coord is 'z':
#                 p0 = [max(self.n_macroparticles), self.Beam.mean_z, self.Beam.sigma_z]
            p0 = [max(self.n_macroparticles), self.Beam.mean_z, self.Beam.sigma_z]
        else:
            p0 = [max(self.n_macroparticles), self.bp_gauss, self.bl_gauss/4]
                                                                
        self.pfit_gauss = curve_fit(gauss, self.bins_centers, self.n_macroparticles, p0)[0] 
        self.bl_gauss = 4 * abs(self.pfit_gauss[2]) 
        self.bp_gauss = abs(self.pfit_gauss[1])

    
    
    def beam_spectrum(self, n_sampling_fft):
        '''
        *Beam spectrum calculation, to be extended (normalized profile, different
        coordinates, etc.)*
        '''
         
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
            derivative = ndimage.gaussian_filter1d(self.n_macroparticles, sigma=1, 
                                                   order=1, mode='wrap') / dist_centers
        if mode == 2:
            x = self.bins_centers
            derivative = np.gradient(self.n_macroparticles, dist_centers)
             
        return x, derivative
    
    
    def compute_statistics(self):
        '''
        *Compute statistics of each slice (average position of the particles
        in a slice and sigma_rms.*
        
        *Improvement is needed in order to include losses, and transverse
        statistics calculation. Be also careful that empty slices will
        result with NaN values for the statistics.*
        '''
        
        index = np.cumsum(np.append(0, self.n_macroparticles))

        for i in xrange(self.n_slices):
           
            theta  = self.Beam.theta[index[i]:index[i + 1]]
            dE = self.Beam.dE[index[i]:index[i + 1]]

            self.mean_theta[i] = np.mean(theta)
            self.mean_dE[i] = np.mean(dE)

            self.sigma_theta[i] = np.std(theta)
            self.sigma_dE[i] = np.std(dE)

            self.eps_rms_l[i] = np.pi * self.sigma_dE[i] * self.sigma_theta[i] \
                                * self.Beam.ring_radius / (self.Beam.beta_r * c)
     
     
    @property    
    def mean_z(self):
        '''*Average z position of the particles in each slice (needs 
        the compute_statistics_option to be 'on', the head and tail are
        reversed compared to theta and an tau coordinates).*'''
        return - self.mean_theta * self.Beam.ring_radius 
     
    @property    
    def mean_tau(self):
        '''*Average tau position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.mean_theta * self.Beam.ring_radius / (self.Beam.beta_rel * c)

 
    @property
    def mean_delta(self):
        '''*Average delta position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.mean_dE / (self.Beam.beta_rel**2 * self.Beam.energy)

     
    @property    
    def sigma_z(self):
        '''*RMS z position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return - self.sigma_theta * self.Beam.ring_radius 
     
    @property
    def sigma_tau(self):
        '''*RMS tau position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.sigma_theta * self.Beam.ring_radius / (self.Beam.beta_rel * c)
     
    @property
    def sigma_delta(self):
        '''*RMS delta position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.sigma_dE / (self.Beam.beta_rel**2 * self.Beam.energy)



def gauss(x, *p):
    '''
    *Defined as:*
    
    .. math:: A \, e^{\\frac{\\left(x-x_0\\right)^2}{2\\sigma_x^2}}
    
    '''
    
    A, x0, sx = p
    return A*np.exp(-(x-x0)**2/2./sx**2) 
    
        
