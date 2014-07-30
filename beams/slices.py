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
    *Slices class that controls longitudinal discretisation of a bunch. This
    include the bunch profiling (including computation of bunch spectrum,
    derivative, and profile fitting) and the computation of statistics per
    slice.*
    '''

    def __init__(self, n_slices, n_sigma = None, cut_left = None, 
                 cut_right = None, coord = "theta", mode = 'const_space',
                 statistics_option = 'off', fit_option = 'off'):
        
        
        #: *Number of slices*
        self.n_slices = n_slices
        
        #: | *Slicing computation mode*
        #: | *The options are: 'const_space' (default), 'const_charge'.*
        self.mode = mode
        
        #: *Left edge of the slicing (is an optionnal input, in case you use
        #: the 'const_space' mode).*
        self.cut_left = cut_left
        
        #: *Right edge of the slicing (is an optionnal input, in case you use
        #: the 'const_space' mode).*
        self.cut_right = cut_right
        
        #: *Optionnal input parameters, corresponding to the number of*
        #: :math:`\sigma_{RMS}` *of the bunch to slice (this will overwrite
        #: any input of cut_left and cut_right).*
        self.n_sigma = n_sigma
        
        #: | *Type of coordinates in which the slicing is done.*
        #: | *The options are: 'theta' (default), 'tau', 'z'.*
        self.coord = coord
        
        if mode is 'const_space':
            if self.n_sigma is not None:
                pass
            elif cut_left is not None and cut_right is not None:
                pass
        
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
            
        #: *Fit option allows to fit the bunch profile, with the options
        #: 'off' (default), 'gaussian'.*
        self.fit_option = fit_option
            
        if self.fit_option is 'gaussian':
            #: *Bunch length with a gaussian fit (needs fit_option to be 
            #: 'gaussian' defined as* :math:`\tau_{gauss} = 4\sigma`)
            self.bl_gauss = 0

#         if cut_left is not None and cut_right is not None:
#             self.cut_left = cut_left
#             self.cut_right = cut_right
#             if mode != 'const_charge':
#                 self.edges = np.linspace(cut_left, cut_right, self.n_slices + 1)
#                 self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
                
        self.sorted = False
    
    
    def set_longitudinal_cuts(self, bunch):
        
        if self.n_sigma == None:
            if self.coord == "theta":
                self.sort_particles(bunch)
                self.sorted = True
                cut_left = bunch.theta[0] - 0.05*(bunch.theta[-1] - bunch.theta[0])
                cut_right = bunch.theta[-1] + 0.05*(bunch.theta[-1] - bunch.theta[0])
            elif self.coord == "z":
                self.sort_particles(bunch)
                cut_left = bunch.z[0] - 0.05*(bunch.z[-1] - bunch.z[0])
                cut_right = bunch.z[-1] + 0.05*(bunch.z[-1] - bunch.z[0])
            else:
                self.sort_particles(bunch)
                cut_left = bunch.tau[0] - 0.05*(bunch.tau[-1] - bunch.tau[0])
                cut_right = bunch.tau[-1] + 0.05*(bunch.tau[-1] - bunch.tau[0])
        else:
            if self.coord == "theta":
                mean_theta = np.mean(bunch.theta)
                sigma_theta = np.std(bunch.theta)
                cut_left = mean_theta - self.n_sigma * sigma_theta / 2
                cut_right = mean_theta + self.n_sigma * sigma_theta / 2
            elif self.coord == "z":
                mean_z = np.mean(bunch.z)
                sigma_z = np.std(bunch.z)
                cut_left = mean_z - self.n_sigma * sigma_z / 2
                cut_right = mean_z + self.n_sigma * sigma_z / 2
            else:
                mean_tau = np.mean(bunch.tau)
                sigma_tau = np.std(bunch.tau)
                cut_left = mean_tau - self.n_sigma * sigma_tau / 2
                cut_right = mean_tau + self.n_sigma * sigma_tau / 2
            
        return cut_left, cut_right


    def slice_constant_space(self, bunch):
        
        try:
            cut_left, cut_right = self.cut_left, self.cut_right
        except AttributeError:
            cut_left, cut_right = self.set_longitudinal_cuts(bunch)
            self.edges = np.linspace(cut_left, cut_right, self.n_slices + 1)
            self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
        
        if self.sorted == False:
            self.sort_particles(bunch)
            self.sorted = True

        if self.coord == 'z':
            
            self.first_index_in_bin = np.searchsorted(bunch.z, self.edges)
            if cut_right <= bunch.z[-1] and \
                  cut_right == bunch.z[self.first_index_in_bin[-1]]:
                list_z = bunch.z[self.first_index_in_bin[-1]:].tolist()
                self.first_index_in_bin[-1] += list_z.count(bunch.z[
                                                self.first_index_in_bin[-1]])
            
        elif self.coord == 'theta':
            
            self.first_index_in_bin = np.searchsorted(bunch.theta, self.edges)
            if cut_right <= bunch.theta[-1] and \
                  cut_right == bunch.theta[self.first_index_in_bin[-1]]:
                list_theta = bunch.theta[self.first_index_in_bin[-1]:].tolist()
                self.first_index_in_bin[-1] += list_theta.count(bunch.theta[
                                                self.first_index_in_bin[-1]])
        
        else:
            
            self.first_index_in_bin = np.searchsorted(bunch.tau, self.edges)
            if cut_right <= bunch.tau[-1] and \
                  cut_right == bunch.tau[self.first_index_in_bin[-1]]:
                list_tau = bunch.tau[self.first_index_in_bin[-1]:].tolist()
                self.first_index_in_bin[-1] += list_tau.count(bunch.tau[
                                                self.first_index_in_bin[-1]])
            
        self.n_macroparticles = np.diff(self.first_index_in_bin)
        
        
    def slice_constant_space_histogram(self, bunch):
        
        try:
            cut_left, cut_right = self.cut_left, self.cut_right
        except AttributeError:
            cut_left, cut_right = self.set_longitudinal_cuts(bunch)
            self.edges = np.linspace(cut_left, cut_right, self.n_slices + 1)
            self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2
        
        if self.coord == 'theta':
            self.n_macroparticles = np.histogram(bunch.theta, self.edges)[0]
        elif self.coord == 'z':
            self.n_macroparticles = np.histogram(bunch.z, self.edges)[0]
        else:
            self.n_macroparticles = np.histogram(bunch.tau, self.edges)[0]

       
    def slice_constant_charge(self, bunch):
        
        try:
            cut_left, cut_right = self.cut_left, self.cut_right
        except AttributeError:
            cut_left, cut_right = self.set_longitudinal_cuts(bunch)
        
        if self.sorted == False:
            self.sort_particles(bunch)
            self.sorted = True

        n_cut_left = np.searchsorted(bunch.z, cut_left)
        n_cut_right = np.searchsorted(bunch.z, cut_right)
        
        q0 = self.n_macroparticles - (n_cut_right - n_cut_left)
        if (cut_right == bunch.z[n_cut_right]):
            list_z = bunch.z[n_cut_right:].tolist()
            q0 += list_z.count(bunch.z[n_cut_right])
        
        ix = sample(range(self.n_slices), q0 % self.n_slices)
        self.n_macroparticles = (q0 // self.n_slices) * np.ones(self.n_slices)
        self.n_macroparticles[ix] += 1

        n_macroparticles_all = np.hstack((n_cut_left, self.n_macroparticles))
        self.first_index_in_bin = np.cumsum(n_macroparticles_all)
        
        self.edges = (bunch.z[self.first_index_in_bin[1:-1] - 1] + 
                        bunch.z[self.first_index_in_bin[1: -1]]) / 2
        self.edges = np.hstack((cut_left, self.edges, cut_right))
        self.bins_centers = (self.edges[:-1] + self.edges[1:]) / 2


    def track(self, bunch):
        
        self.sorted = False
        self.bunch = bunch
       
        if self.mode == 'const_charge':
            self.slice_constant_charge(bunch)
        elif self.mode == 'const_space':
            self.slice_constant_space(bunch)
        elif self.mode == 'const_space_hist':
            self.slice_constant_space_histogram(bunch)
        else:
            raise RuntimeError('Choose one of the three slicing methods!')


#     def compute_statistics(self, bunch):
#          
#         if self.sorted == False:
#             self.sort_particles(bunch)
#  
#         index = self.first_index_in_bin[0] + \
#                 np.cumsum(np.append(0, self.n_macroparticles))
#  
#         for i in xrange(self.n_slices):
#             
#             x  = bunch.x[index[i]:index[i + 1]]
#             xp = bunch.xp[index[i]:index[i + 1]]
#             y  = bunch.y[index[i]:index[i + 1]]
#             yp = bunch.yp[index[i]:index[i + 1]]
#             theta  = bunch.theta[index[i]:index[i + 1]]
#             dE = bunch.dE[index[i]:index[i + 1]]
#  
#             self.mean_x[i] = np.mean(x)
#             self.mean_xp[i] = np.mean(xp)
#             self.mean_y[i] = np.mean(y)
#             self.mean_yp[i] = np.mean(yp)
#             self.mean_theta[i] = np.mean(theta)
#             self.mean_dE[i] = np.mean(dE)
#  
#             self.sigma_x[i] = np.std(x)
#             self.sigma_y[i] = np.std(y)
#             self.sigma_theta[i] = np.std(theta)
#             self.sigma_dE[i] = np.std(dE)
#  
#             self.epsn_x[i] = cp.emittance(x, xp) * bunch.gamma_rel * \
#                                 bunch.beta_rel * 1e6
#             self.epsn_y[i] = cp.emittance(y, yp) * bunch.gamma_rel * \
#                                 bunch.beta_rel * 1e6
#             self.eps_rms_l[i] = np.pi * self.sigma_dE[i] * self.sigma_theta[i] \
#                                 * bunch.ring_radius / (bunch.beta_rel * c)

    
    def sort_particles(self, bunch):
       
        if self.coord == 'theta' or self.coord == 'tau':
            
            argsorted = np.argsort(bunch.theta)
            
        elif self.coord == 'z':
            
            argsorted = np.argsort(bunch.z)
        
#         bunch.x = bunch.x.take(argsorted)
#         bunch.xp = bunch.xp.take(argsorted)
#         bunch.y = bunch.y.take(argsorted)
#         bunch.yp = bunch.yp.take(argsorted)
        bunch.theta = bunch.theta.take(argsorted)
        bunch.dE = bunch.dE.take(argsorted)
        bunch.id = bunch.id.take(argsorted)
    
    
    def beam_spectrum(self, n_sampling_fft):
        
        spectrum = rfft(self.n_macroparticles, n_sampling_fft)
            
        return spectrum
    
    
    def beam_profile_derivative(self, mode):
        
        ''' The input is one of the two available methods for differentiating
            a function. The two outputs are the coordinate step and the discrete
            derivative of the bunch profile respectively'''
        
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
        return - self.mean_theta * self.bunch.ring_radius 
    @mean_z.setter
    def mean_z(self, value):
        self.mean_theta = - value / self.bunch.ring_radius 
    
    @property    
    def mean_tau(self):
        '''*Average tau position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.mean_theta * self.ring_radius / (self.beta_rel * c)
    @mean_tau.setter
    def mean_tau(self, value):
        self.mean_theta = value * self.beta_rel * c / self.ring_radius

    @property
    def mean_delta(self):
        '''*Average delta position of the particles in each slice (needs 
        the compute_statistics_option to be 'on').*'''
        return self.mean_dE / (self.bunch.beta_rel**2 * self.bunch.energy)
    @mean_delta.setter
    def mean_delta(self, value):
        self.mean_dE = value * self.bunch.beta_rel**2 * self.bunch.energy
    
    @property    
    def sigma_z(self):
        return - self.sigma_theta * self.bunch.ring_radius 
    @sigma_z.setter
    def sigma_z(self, value):
        self.sigma_theta = - value / self.bunch.ring_radius 
    
    @property
    def sigma_tau(self):
        return self.sigma_theta * self.ring_radius / (self.beta_rel * c)
    @sigma_tau.setter
    def sigma_tau(self, value):
        self.sigma_theta = value * self.beta_rel * c / self.ring_radius
    
    @property
    def sigma_delta(self):
        return self.sigma_dE / (self.bunch.beta_rel**2 * self.bunch.energy)
    @sigma_delta.setter
    def sigma_delta(self, value):
        self.sigma_dE = value * self.bunch.beta_rel**2 * self.bunch.energy
    
        
