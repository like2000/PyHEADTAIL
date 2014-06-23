'''
Created on 12.06.2014

@author: Kevin Li, Danilo Quartullo, Helga Timko
'''

import numpy as np
#import copy, h5py, sys
import sys
from scipy.constants import c, e
import cobra_functions.stats as cp
from scipy.optimize import curve_fit


class Beam(object):
    
    def __init__(self, General_parameters, n_macroparticles, intensity):
        
        # Beam and ring-dependent properties
        self.ring_radius = General_parameters.ring_radius
        self.mass = General_parameters.mass # in kg
        self.charge = General_parameters.charge # in C
        self.intensity = intensity # total no of particles
        
        self.beta_rel = General_parameters.beta_rel_program[0] #: Relativistic beta of the synchronous particle
        self.gamma_rel = General_parameters.gamma_rel_program[0] #: Relativistic gamma of the synchronous particle
        self.energy = General_parameters.energy_program[0] #: Energy of the synchronous particle [eV]
        self.momentum = General_parameters.momentum_program[0] #: Momentum of the synchronous particle [eV/c]

        # Beam coordinates
        self.x = np.empty([n_macroparticles])
        self.xp = np.empty([n_macroparticles])
        self.y = np.empty([n_macroparticles])
        self.yp = np.empty([n_macroparticles])
        self.theta = np.empty([n_macroparticles])
        self.dE = np.empty([n_macroparticles])
        
        # Initial coordinates (e.g. for ecloud)
        self.x0 = np.empty([n_macroparticles])
        self.xp0 = np.empty([n_macroparticles])
        self.y0 = np.empty([n_macroparticles])
        self.yp0 = np.empty([n_macroparticles])
        self.theta0 = np.empty([n_macroparticles])
        self.dE0 = np.empty([n_macroparticles])
        
        # Transverse and longitudinal properties, statistics       
        self.alpha_x = 0
        self.beta_x = 0
        self.epsn_x = 0
        self.alpha_y = 0
        self.beta_y = 0
        self.epsn_y = 0
        self.sigma_theta = 0
        self.sigma_dE = 0
        
        # Particle/loss counts
        self.n_macroparticles = n_macroparticles
        self.n_macroparticles_lost = 0
        self.id = np.arange(1, self.n_macroparticles + 1, dtype=int)
        
        # Boolean value: is the beam sliced?
        self.beam_is_sliced = False
        self.slicing = None

            
    # Coordinate conversions
    @property
    def z(self):
        return - self.theta * self.ring_radius
     
    @z.setter
    def z(self, value):
        self.theta = - value / self.ring_radius
    
    @property
    def delta(self):
        return self.dE / (self.beta_rel**2 * self.ring.energy_i(self))

    @delta.setter
    def delta(self, value):
        self.dE = value * self.beta_rel**2 * self.ring.energy_i(self)

    @property
    def z0(self):
        return - self.theta0 * self.ring_radius
    @z0.setter
    def z0(self, value):
        self.theta0 = - value / self.ring_radius
    
    @property
    def delta0(self):
        return self.dE0 / (self.beta_rel**2 * self.ring.energy_i(self))
    @delta0.setter
    def delta0(self, value):
        self.dE0 = value * self.beta_rel**2 * self.ring.energy_i(self)

    def reinit(self):

        np.copyto(self.x, self.x0)
        np.copyto(self.xp, self.xp0)
        np.copyto(self.y, self.y0)
        np.copyto(self.yp, self.yp0)
        np.copyto(self.theta, self.theta0)
        np.copyto(self.dE, self.dE0)
        np.copyto(self.z, self.z0)
        np.copyto(self.delta, self.delta0)
        
    # Statistics
    @property    
    def mean_z(self):
        return - self.mean_theta * self.ring_radius
    @mean_z.setter
    def mean_z(self, value):
        self.mean_theta = - value / self.ring_radius
    
    @property
    def mean_delta(self):
        return self.mean_dE / (self.beta_rel**2 * self.energy)
    @mean_delta.setter
    def mean_delta(self, value):
        self.mean_dE = value * self.beta_rel**2 * self.energy

    @property    
    def sigma_z(self):
        return - self.sigma_theta * self.ring_radius
    @sigma_z.setter
    def sigma_z(self, value):
        self.sigma_theta = - value / self.ring_radius
    
    @property
    def sigma_delta(self):
        return self.sigma_dE / (self.beta_rel**2 * self.energy)
    @sigma_delta.setter
    def sigma_delta(self, value):
        self.sigma_dE = value * self.beta_rel**2 * self.energy

    
    def longit_statistics(self, gaussian_fit = "Off"):
        
        self.mean_theta = cp.mean(self.theta)
        self.mean_dE = cp.mean(self.dE)
        self.sigma_theta = cp.std(self.theta)
        self.sigma_dE = cp.std(self.dE)
        
        ##### R.m.s. emittance in Gaussian approximation, other emittances to be defined
        self.epsn_rms_l = np.pi * self.sigma_dE * self.sigma_theta \
                        * self.ring_radius / (self.beta_rel * c) # in eVs

        ##### Gaussian fit to theta-profile
        if self.beam_is_sliced == True and gaussian_fit == "On":
            self.slicing.compute_statistics(self)
            p0 = [max(self.slicing.n_macroparticles), self.mean_theta, self.sigma_theta] 
            def gauss(x, *p):
                A, x0, sx = p
                return A*np.exp(-(x-x0)**2/2./sx**2) 
            pfit, pvar = curve_fit(gauss, self.slicing.mean_theta[0:], 
                                   self.slicing.n_macroparticles[0:], p0 = p0)
            self.bl_gauss = 4*abs(pfit[2]) 

                                
    def transv_statistics(self):
        
        self.mean_x = cp.mean(self.x)
        self.mean_xp = cp.mean(self.xp)
        self.mean_y = cp.mean(self.y)
        self.mean_yp = cp.mean(self.yp)
        self.sigma_x = cp.std(self.x)
        self.sigma_y = cp.std(self.y)
        self.epsn_x_xp = cp.emittance(self.x, self.xp) * self.gamma_rel \
                        * self.beta_rel * 1e6
        self.epsn_y_yp = cp.emittance(self.y, self.yp) * self.gamma_rel \
                        * self.beta_rel * 1e6
    
    def losses(self, ring):
         
        for i in xrange(self.n_macroparticles):
            if not ring.is_in_separatrix(ring, self, self.theta[i], self.dE[i], self.delta[i]):
                # Set ID to zero
                self.id[i] = 0



