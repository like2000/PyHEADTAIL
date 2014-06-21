'''
Created on 12.06.2014

@author: Danilo Quartullo, Helga Timko, Adrian Oeftiger, Alexandre Lasheen
'''


from __future__ import division
import numpy as np
import math
import sys
from scipy.constants import c


def eta_tracking(self, General_parameters, delta, counter):
    
    """
    Depending on the number of entries in self.alpha_array, the slippage factor
    \eta = \sum_i \eta_i * \delta^i is calculated to the corresponding order.

    As eta is used in the tracker, it is calculated with the initial momentum
    at that time step, and the corresponding relativistic beta and gamma.
    """
    eta = 0
    for i in xrange( General_parameters.alpha_array.size ):   # order = len - 1
        eta_i = getattr(General_parameters, 'eta' + str(i))[counter[0]]
        eta  += eta_i * (delta**i)
    return eta


class Kick(object):
    
    """The Kick represents the kick(s) by an RF station at a certain position 
    of the ring. The kicks are summed over the different harmonic RF systems 
    in the station.

    The cavity phase can be shifted by the user via phi_offset."""

    def __init__(self, n_rf_systems, harmonic_numbers_list, voltage_program_list, phi_offset_list, counter):
        
        self.counter = counter
        self.n_rf_systems = n_rf_systems
        self.harmonic = harmonic_numbers_list
        self.voltage = voltage_program_list
        self.phi_offset = phi_offset_list
        
    def track(self, beam):
        for i in range(self.n_rf_systems):
            beam.dE += self.voltage[i][self.counter[0]] * np.sin(self.harmonic[i][self.counter[0]] * beam.theta + 
                                             self.phi_offset[i][self.counter[0]]) # in eV
    
    
class Kick_acceleration(object):
    
    """Kick_acceleration gives a single accelerating kick to the bunch. 
    The accelerating kick is defined by the change in the design momentum 
    (synchronous momentum). 
    
    The acceleration is assumed to be distributed over the length of the 
    RF station, so the average beta is used in the calculation of the kick."""
    
    def __init__(self, beta_rel_program, p_increment, counter):
        
        self.counter = counter
        self.beta = beta_rel_program
        self.p_increment = p_increment
        
    def track(self, beam):
        
        beam.dE += - self.beta_rel_program[self.counter[0]] * self.p_increment[self.counter[0]] # in eV
#         beam.beta_rel
#         beam.gamma_rel
#         beam.energy
        

class Drift(object):
    
    """The drift updates the longitudinal coordinate of the particle after 
    applying the energy kick; self.length is the drift length.  

    The correction factor \\beta_{n+1} / \\beta_n is necessary when the 
    synchronous energy is low and the range is synchronous energy is large,
    to avoid a shrinking phase space."""

    def __init__(self, drift_length, General_parameters, solver):
        
        self.drift_length = drift_length
        self.General_parameters = General_parameters
        self.ring_circumference = General_parameters.ring_circumference
        self.counter = General_parameters.counter
        self.beta_rel_program = General_parameters.beta_rel_program
        self.solver = solver
                
    def track(self, beam):  
        try: 
            beam.theta = \
            {'full' : self.beta_rel_program[self.counter[0] + 1] / self.beta_rel_program[self.counter[0]] * beam.theta \
                    + 2 * np.pi * (1 / (1 - eta_tracking(self.General_parameters, beam.delta, self.counter) \
                    * beam.delta) - 1) * self.drift_length / self.ring_circumference,
             'simple' : beam.theta + 2 * np.pi * self.General_parameters.eta0 * beam.delta * \
                    self.drift_length / self.ring_circumference
            }[self.solver]
        except KeyError:
            print "ERROR: Choice of longitudinal solver not recognized! Aborting..."
            sys.exit()
        

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
    
    def __init__(self, General_parameters, RF_parameters_section):
        
        self.kick = Kick(RF_parameters_section.n_rf_systems, 
                          RF_parameters_section.harmonic_numbers_list, RF_parameters_section.voltage_program_list, 
                          RF_parameters_section.phi_offset_list, General_parameters.counter)
        
        p_increment = RF_parameters_section.momentum_program[1:] - RF_parameters_section.momentum_program[1:]
                    
        self.kick_acceleration = Kick_acceleration(General_parameters.beta_rel_program, p_increment, General_parameters.counter)
        
        solver = 'full' # TODO : put this as a paramter
        self.drift = Drift(RF_parameters_section.section_length, General_parameters, solver)
        
        self.elements = [self.kick] + [self.kick_acceleration] + [self.drift]
        
    def track(self, beam):
        for longMap in self.elements:
            longMap.track(beam)
    
    
class Full_Ring_and_RF(object):
    '''
    Full ring object, containing the total map of Ring_and_RFstation
    '''
    
    def __init__(self, General_parameters, sum_RF_section_parameters):
        
        assert General_parameters.ring_circumference == sum_RF_section_parameters.section_length_sum
        
        self.n_sections = sum_RF_section_parameters.total_n_sections #: Passing the number of sections
        
        self.Ring_and_RFstation_list = [] #: List of ring and RF stations
         
        for i in range(self.n_sections):
            self.Ring_and_RFstation_list.append(Ring_and_RFstation(General_parameters, sum_RF_section_parameters.RF_section_parameters_list[i]))
            
    def track(self, beam):
        
        for i in range(self.n_sections):
            self.Ring_and_RFstation_list[i].track(beam)
    
        
# class Longitudinal_tracker(object):
#     
#     """
#         The Longitudinal_tracker tracks the bunch through a given RF station
#         and takes care that kicks and the drift are done in correct order.
#         
#         Different solvers can be used:
#         
#         'full' -- accurate solution of the drift
#         
#         'simple' -- drift with no correction for low energy/large energy range and zeroth order in the slippage factor
#         
#         For de-bunching, simply pass zero voltage.
#         For synchrotron radiation, energy loss term yet to be implemented.
#     """
# 
#     def __init__(self, ring, solver='full'): 
#         
#         """self.p_increment is the momentum step per turn of the synchronous 
#         particle (defined via user input, see ring_and_RFstation).
#         See the Kick_acceleration class for further details."""
#         
#         self.solver = solver
#         self.ring = ring
#         self.harmonic_list = ring.harmonic
#         self.voltage_list = ring.voltage
#         self.circumference = ring.circumference
#         
#         if ring.phi_offset is not None:
#             if not len(ring.harmonic) == len(ring.voltage) == len(ring.phi_offset):
#                 print ("Warning: parameter lists for RFSystems do not have the same length!")
# 
# 
#         """Separating the kicks from the RF and the magnets.
#         kick can contain multiple contributions:
#         self.kicks -- kick due to RF station passage
#         self.kick_acceleration -- kick due to acceleration
#         self.elements contains the full map of kicks and drift in the RF station."""
#         self.kicks = []
#         for i in xrange(len(ring.harmonic)):
#             kick = Kick(ring, i)
#             self.kicks.append(kick)
#         self.kick_acceleration = Kick_acceleration(ring, 0)
#         self.elements = self.kicks + [self.kick_acceleration] + [Drift(ring, solver)]
#         
#     def track(self, beam):
#         
#         self.kick_acceleration.p_increment = self.ring.p0_f() - self.ring.p0_i()
#         for longMap in self.elements:
#             longMap.track(beam)
        
    
class LinearMap(object):
    
    '''
    Linear Map represented by a Courant-Snyder transportation matrix.
    self.alpha is the linear momentum compaction factor.
    Qs is forced to be constant.
    '''

    def __init__(self, ring, Qs):

        """alpha is the linear momentum compaction factor,
        Qs the synchroton tune."""
        
        self.circumference = ring.circumference
        self.alpha = ring.alpha_array[0]
        self.Qs = Qs

    def track(self, beam):

        eta = self.alpha - 1 / beam.ring.gamma_i(beam)**2

        omega_0 = 2 * np.pi * beam.ring.beta_i(beam) * c / self.circumference
        omega_s = self.Qs * omega_0

        dQs = 2 * np.pi * self.Qs
        cosdQs = np.cos(dQs)
        sindQs = np.sin(dQs)

        z0 = beam.z
        delta0 = beam.delta

        beam.z = z0 * cosdQs - eta * c / omega_s * delta0 * sindQs
        beam.delta = delta0 * cosdQs + omega_s / eta / c * z0 * sindQs


