'''
Created on 12.06.2014

@author: Danilo Quartullo, Helga Timko, Alexandre Lasheen
'''

from __future__ import division
import numpy as np
from warnings import filterwarnings
from scipy.constants import c, e
        

class Ring_and_RFstation(object):
    '''
    Definition of an RF station and part of the ring until the next station, see figure.
    
    .. image:: https://raw.githubusercontent.com/like2000/PyHEADTAIL/PYlongitudinal/doc/source/ring_and_RFstation.png
        :align: center
        :width: 600
        :height: 600
        
    The time step is fixed to be one turn, but the tracking can consist of multiple
    ring_and_RFstation objects. In this case, the user should make sure that the lengths
    of the stations sum up exactly to the circumference.
    Each RF station may contain several RF harmonic systems which are considered to be
    in the same location. First, a kick from the cavity voltage(s) is applied, then a
    kick from the accelerating magnets, and finally a drift. 
    If one wants to do minimal longitudinal tracking, defining the RF systems in the 
    station is not necessary. 
    '''
    
    def __init__(self, Global_parameters, length=None, 
                 harmonic_list=None, voltage_list=None, phi_offset_list=None):
                    
        alpha_array = Global_parameters.momentum_compaction_array
        circumference = Global_parameters.ring_circumference
        
        if alpha_array != None and len(alpha_array) > 3:
            print "WARNING: Slippage factor implemented only till second order. Higher orders in alpha ignored. "
        self.counter = 0 # To step in the momentum program  
                      
        # Optional parameters
        if length != None and circumference != length:
            print "ATTENTION: The total length of RF stations should sum up to the circumference."
        self.length = length # in m              
        self.harmonic = harmonic_list
        self.voltage = voltage_list # in V
        self.phi_offset = phi_offset_list # in rad

    # Derived energy-related properties 
    # Energy and momentum in units of eV   
    '''The relativistic beta, gamma, and energy of the synchronous particle are
    derived from the synchronous momentum provided by the user.
    Energy and momentum are in units of eV, while all other quantities are SI.'''
    def p0_i(self):
        return self.momentum_program[self.counter]
    
    def p0_f(self):
        return self.momentum_program[self.counter + 1]

    def p0(self):
        return (self.p0_i() + self.p0_f()) / 2   
        
    def beta_i(self, beam):
        return np.sqrt( 1 / (1 + (beam.mass * c**2)**2 / (self.p0_i() * e)**2) )
        
    def beta_f(self, beam):
        return np.sqrt( 1 / (1 + (beam.mass * c**2)**2 / (self.p0_f() * e)**2) )
 
    def beta(self, beam):
        return (self.beta_i(beam) + self.beta_f(beam)) / 2
        
    def gamma_i(self, beam):
        return np.sqrt( 1 + (self.p0_i() * e)**2 / (beam.mass * c**2)**2 )
    
    def gamma_f(self, beam):
        return np.sqrt( 1 + (self.p0_f() * e)**2 / (beam.mass * c**2)**2 )
    
    def gamma(self, beam):
        return (self.gamma_i(beam) + self.gamma_f(beam)) / 2
    
    def energy_i(self, beam):
        return np.sqrt( self.p0_i()**2 + (beam.mass * c**2 / e)**2 )
    
    def energy_f(self, beam):
        return np.sqrt( self.p0_f()**2 + (beam.mass * c**2 / e)**2 )

    def energy(self, beam):
        return (self.energy_i(beam) + self.energy_f(beam)) / 2    
 
    
#    def potential(self, z, beam):
        
#        """the potential well of the rf system"""
#        phi_0 = self.accelerating_kick.calc_phi_0(beam)
#        h1 = self.accelerating_kick.harmonic
#        def fetch_potential(kick):
#            phi_0_i = kick.harmonic / h1 * phi_0
#            return kick.potential(z, beam, phi_0_i)
#        potential_list = map(fetch_potential, self.kicks)
#        return sum(potential_list)
    
    def hamiltonian(self, beam, theta, dE, delta):
        """Single RF sinusoidal Hamiltonian.
        Uses beta, energy averaged over the turn.
        To be generalized."""
        h0 = self.harmonic[0]
        V0 = self.voltage[0]
        c1 = self.eta(beam, delta) * c * np.pi / (self.circumference * 
             self.beta_i(beam) * self.energy_i(beam) )
        c2 = c * self.beta_i(beam) * V0 / (h0 * self.circumference)
        phi_s = self.calc_phi_s(beam, self.voltage)

        return c1 * dE**2 + c2 * (np.cos(h0 * theta) - np.cos(phi_s) + 
                                   (h0 * theta - phi_s) * np.sin(phi_s))

    def calc_phi_s(self, beam, voltage):
        """The synchronous phase calculated from the rate of momentum change.
        Below transition, for decelerating bucket: phi_s is in (-Pi/2,0)
        Below transition, for accelerating bucket: phi_s is in (0,Pi/2)
        Above transition, for accelerating bucket: phi_s is in (Pi/2,Pi)
        Above transition, for decelerating bucket: phi_s is in (Pi,3Pi/2)
        The synchronous phase is calculated at a certain moment.
        Uses beta, energy averaged over the turn."""
        V0 = voltage[0]
        phi_s = np.arcsin(self.beta(beam) * (self.p0_f() - self.p0_i()) / V0 )
        if self.eta(beam, 0) > 0:
            phi_s = np.pi - phi_s

        return phi_s      
    
    def separatrix(self, beam, theta):
        """Single RF sinusoidal separatrix.
        Uses beta, energy averaged over the turn.
        To be generalized."""
        h0 = self.harmonic[0]
        V0 = self.voltage[0]
        phi_s = self.calc_phi_s(beam, self.voltage)
        
        filterwarnings('ignore')
        
        separatrix_array = np.sqrt(self.beta_i(beam)**2 * self.energy_i(beam) *
                        V0 / (np.pi * self.eta(beam, 0) * h0) * 
                        (-np.cos(h0 * theta) - np.cos(phi_s) + 
                        (np.pi - phi_s - h0 * theta) * np.sin(phi_s)))
        
        filterwarnings('default')
           
        return separatrix_array

    def is_in_separatrix(self, beam, theta, dE, delta):
        """Condition for being inside the separatrix.
        Single RF sinusoidal.
        Uses beta, energy averaged over the turn.
        To be generalized."""
        h0 = self.harmonic[0]
        phi_s = self.calc_phi_s(beam, self.voltage)
        Hsep = self.hamiltonian(beam, (np.pi - phi_s) / h0, 0, 0) 
        isin = np.fabs(self.hamiltonian(beam, theta, dE, delta)) < np.fabs(Hsep)

        return isin

    def eta(self, beam, delta):
        
        """Depending on the number of entries in self.alpha_array, the slippage factor
        \eta = \sum_i \eta_i * \delta^i is calculated to the corresponding order.

        As eta is used in the tracker, it is calculated with the initial momentum
        at that time step, and the corresponding relativistic beta and gamma.
        For eta coefficients, see Lee: Accelerator Physics (Wiley)."""
        eta = 0
        for i in xrange( len(self.alpha_array) ):   # order = len - 1
            eta_i = getattr(self, '_eta' + str(i))(beam, self.alpha_array)
            eta  += eta_i * (delta**i)
        return eta

    
    def _eta0(self, beam, alpha_array):
        
        return alpha_array[0] - self.gamma_i(beam)**-2
   
    
    def _eta1(self, beam, alpha_array):
        
        return 3 * self.beta_i(beam)**2 / (2 * self.gamma_i(beam)**2) + alpha_array[1] \
            - alpha_array[0] * (alpha_array[0] - self.gamma_i(beam)**-2)
    
    
    def _eta2(self, beam, alpha_array):
        
        return - self.beta_i(beam)**2 * (5 * self.beta_i(beam)**2 - 1) / (2 * self.gamma_i(beam)**2) \
            + alpha_array[2] - 2 * alpha_array[0] * alpha_array[1] + alpha_array[1] \
            / self.gamma_i(beam)**2 + alpha_array[0]**2 * (alpha_array[0] - self.gamma_i(beam)**-2) \
            - 3 * self.beta_i(beam)**2 * alpha_array[0] / (2 * self.gamma_i(beam)**2)
    
    
