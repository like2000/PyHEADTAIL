'''
Created on 12.06.2014

@author: Danilo Quartullo, Helga Timko, Adrian Oeftiger, Alexandre Lasheen
'''


from __future__ import division
import numpy as np
import math
import sys
from scipy.constants import c
from warnings import filterwarnings


def eta_tracking(General_parameters, delta, counter, index_section = 0):
    
    """
    Depending on the number of entries in self.alpha_array, the slippage factor
    \eta = \sum_i \eta_i * \delta^i is calculated to the corresponding order.

    As eta is used in the tracker, it is calculated with the initial momentum
    at that time step, and the corresponding relativistic beta and gamma.
    """
    eta = 0
    for i in xrange( General_parameters.alpha_array.size ):   # order = len - 1
        eta_i = getattr(General_parameters, 'eta' + str(i))[index_section][counter[0]]
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
    
    def __init__(self, General_parameters, p_increment, counter, index_section = 0):
        
        self.counter = counter
        self.momentum_program = General_parameters.momentum_program
        self.p_increment = p_increment
        
        if General_parameters.beta_rel_program.size == General_parameters.n_turns + 1:
            self.beta_rel_program = General_parameters.beta_rel_program
            self.gamma_rel_program = General_parameters.gamma_rel_program
            self.energy_program = General_parameters.energy_program
        else:
            self.beta_rel_program = General_parameters.beta_rel_program[index_section]
            self.gamma_rel_program = General_parameters.gamma_rel_program[index_section]
            self.energy_program = General_parameters.energy_program[index_section]
        
    def track(self, beam):
        
        beam.dE += - self.beta_rel_program[self.counter[0]] * self.p_increment[self.counter[0]] # in eV
        
        # Updating the beam synchronous momentum
        beam.beta_rel = self.beta_rel_program[self.counter[0] + 1]
        beam.gamma_rel = self.gamma_rel_program[self.counter[0] + 1]
        beam.energy = self.energy_program[self.counter[0] + 1]
        beam.momentum = self.momentum_program[self.counter[0] + 1]
        

class Drift(object):
    
    """The drift updates the longitudinal coordinate of the particle after 
    applying the energy kick; self.length is the drift length.  

    The correction factor \\beta_{n+1} / \\beta_n is necessary when the 
    synchronous energy is low and the range is synchronous energy is large,
    to avoid a shrinking phase space."""

    def __init__(self, drift_length, General_parameters, solver, index_section = 0):
        
        
        self.index_section = index_section
        self.drift_length = drift_length
        self.General_parameters = General_parameters
        self.ring_circumference = General_parameters.ring_circumference
        self.counter = General_parameters.counter
        self.solver = solver
        self.beta_rel_program = General_parameters.beta_rel_program[index_section]

                
    def track(self, beam):  

        if self.solver == 'full': 
            beam.theta = self.beta_rel_program[self.counter[0] + 1] / self.beta_rel_program[self.counter[0]] * beam.theta \
                    + 2 * np.pi * (1 / (1 - eta_tracking(self.General_parameters, beam.delta, self.counter, self.index_section) \
                    * beam.delta) - 1) * self.drift_length / self.ring_circumference
        elif self.solver == 'simple':
            beam.theta = beam.theta + 2 * np.pi * self.General_parameters.eta0[self.index_section][self.counter[0]] * beam.delta * \
                    self.drift_length / self.ring_circumference
        else:
            raise RuntimeError("ERROR: Choice of longitudinal solver not recognized! Aborting...")
        

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
        
    def __init__(self, General_parameters, RF_parameters_section, index_section = 0):
        
        self.index_section = index_section

        self.kick = Kick(RF_parameters_section.n_rf_systems, 
                          RF_parameters_section.harmonic_numbers_list, RF_parameters_section.voltage_program_list, 
                          RF_parameters_section.phi_offset_list, General_parameters.counter)
        
        solver = 'simple' # TODO : put this as a parameter
        self.drift = Drift(RF_parameters_section.section_length, General_parameters, solver, index_section)
        
        if np.sum(RF_parameters_section.p_increment) == 0:
            self.elements = [self.kick] + [self.drift]
        else:
            self.kick_acceleration = Kick_acceleration(General_parameters, RF_parameters_section.p_increment, General_parameters.counter, index_section)
            self.elements = [self.kick] + [self.kick_acceleration] + [self.drift]
            
        self.phi_s = calc_phi_s(General_parameters, RF_parameters_section)

          
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
            self.Ring_and_RFstation_list.append(Ring_and_RFstation(General_parameters, sum_RF_section_parameters.RF_section_parameters_list[i], i))
            
    def track(self, beam):
        
        for i in range(self.n_sections):
            self.Ring_and_RFstation_list[i].track(beam)


def calc_phi_s(General_parameters, RF_parameters_section, accelerating_systems = 'all', index_section = 0):
    """The synchronous phase calculated from the rate of momentum change.
    Below transition, for decelerating bucket: phi_s is in (-Pi/2,0)
    Below transition, for accelerating bucket: phi_s is in (0,Pi/2)
    Above transition, for accelerating bucket: phi_s is in (Pi/2,Pi)
    Above transition, for decelerating bucket: phi_s is in (Pi,3Pi/2)
    The synchronous phase is calculated at a certain moment.
    Uses beta, energy averaged over the turn."""
    
#     if General_parameters.beta_rel_program.size == General_parameters.n_turns + 1:
#         beta_rel_program = General_parameters.beta_rel_program
#         eta0 = General_parameters.eta0
#     else:
    beta_rel_program = General_parameters.beta_rel_program[index_section]
    eta0 = General_parameters.eta0[index_section]
        
            
    if RF_parameters_section.n_rf_systems == 1:
        V0 = RF_parameters_section.voltage_program_list[0]
        average_beta = (beta_rel_program[1:] + beta_rel_program[0:-1])/2
        
        phi_s = np.arcsin(average_beta * RF_parameters_section.p_increment / V0 )
        
        phi_s[(eta0[1:] + eta0[0:-1])/2 > 0] = np.pi - phi_s
         
        return phi_s
    
    else:
        '''
        To be implemented
        '''
        if accelerating_systems == 'all':
            '''
            In this case, all the rf_systems are accelerating, phi_s is calculated accordingly
            with respect to the fundamental frequency
            '''
            pass
        elif accelerating_systems == 'first':
            '''
            Only the first rf_system is accelerating, so we have to correct the phi_offset of the
            other rf_systems in order that the p_increment is only caused by the first RF
            '''
            pass
        else:
            raise RuntimeError('Did not recognize the option accelerating_systems in calc_phi_s function')
        

def hamiltonian(General_parameters, RF_parameters_section, theta, dE, delta):
    """Single RF sinusoidal Hamiltonian.
    Uses beta, energy averaged over the turn.
    To be generalized."""
    
    if RF_parameters_section.section_length != General_parameters.ring_circumference:
        raise RuntimeError('WARNING : The hamiltonian is not yet properly computed for several sections !!!')
    
    
    if RF_parameters_section.n_rf_systems == 1:
        counter = General_parameters.counter[0]
        h0 = RF_parameters_section.harmonic_numbers_list[0][counter]
        V0 = RF_parameters_section.voltage_program_list[0][counter]
        
        c1 = eta_tracking(General_parameters, delta, counter) * c * np.pi / (General_parameters.ring_circumference * 
             General_parameters.beta_rel_program[counter] * General_parameters.energy_program[counter] )
        c2 = c * General_parameters.beta_rel_program[counter] * V0 / (h0 * General_parameters.ring_circumference)
        
        phi_s = RF_parameters_section.phi_s[counter]
    
        return c1 * dE**2 + c2 * (np.cos(h0 * theta) - np.cos(phi_s) + 
                                   (h0 * theta - phi_s) * np.sin(phi_s))
        
    else:
        raise RuntimeError('Hamiltonian for multiple RF is not implemeted yet')


def separatrix(General_parameters, RF_parameters_section, theta):
    """Single RF sinusoidal separatrix.
    Uses beta, energy averaged over the turn.
    To be generalized."""
    
    if RF_parameters_section.section_length != General_parameters.ring_circumference:
        print 'WARNING : The separatrix is not yet properly computed for several sections !!!'
    
    
    if RF_parameters_section.n_rf_systems == 1:
        counter = General_parameters.counter[0]
        h0 = RF_parameters_section.harmonic_numbers_list[0][counter]
        V0 = RF_parameters_section.voltage_program_list[0][counter]
        
    else:
        raise RuntimeError('Separatrix for multiple RF is not implemeted yet')

    phi_s = RF_parameters_section.phi_s[counter]
     
    filterwarnings('ignore')
    
    beta_rel_program_average = (General_parameters.beta_rel_program[counter + 1] + General_parameters.beta_rel_program[counter]) / 2
    
    energy_program_average = (General_parameters.energy_program[counter + 1] + General_parameters.energy_program[counter]) / 2
     
    separatrix_array = np.sqrt(beta_rel_program_average**2 * energy_program_average *
                    V0 / (np.pi * General_parameters.eta0 * h0) * 
                    (-np.cos(h0 * theta) - np.cos(phi_s) + 
                    (np.pi - phi_s - h0 * theta) * np.sin(phi_s)))
     
    filterwarnings('default')
        
    return separatrix_array



def is_in_separatrix(General_parameters, RF_parameters_section, theta, dE, delta):
    """Condition for being inside the separatrix.
    Single RF sinusoidal.
    Uses beta, energy averaged over the turn.
    To be generalized."""
    
    if RF_parameters_section.section_length != General_parameters.ring_circumference:
        print 'WARNING : The separatrix is not yet properly computed for several sections !!!'
    
    
    if RF_parameters_section.n_rf_systems == 1:
        counter = General_parameters.counter[0]
        h0 = RF_parameters_section.harmonic_numbers_list[counter]
        
    else:
        raise RuntimeError('is_in_separatrix for multiple RF is not implemeted yet')
        
    phi_s = RF_parameters_section.phi_s[counter]
    
    Hsep = hamiltonian(General_parameters, RF_parameters_section, (np.pi - phi_s) / h0, 0, 0) 
    isin = np.fabs(hamiltonian(General_parameters, RF_parameters_section, theta, dE, delta)) < np.fabs(Hsep)

    return isin
        
        
    
class LinearMap(object):
    
    '''
    Linear Map represented by a Courant-Snyder transportation matrix.
    self.alpha is the linear momentum compaction factor.
    Qs is forced to be constant.
    '''

    def __init__(self, General_parameters, Qs):

        """alpha is the linear momentum compaction factor,
        Qs the synchroton tune."""
        
        self.beta_rel_program = General_parameters.beta_rel_program
        self.gamma_rel_program = General_parameters.gamma_rel_program
        self.energy_program = General_parameters.energy_program
        self.momentum_program = General_parameters.momentum_program
        
        self.ring_circumference = General_parameters.ring_circumference
        self.eta = General_parameters._eta0
        self.Qs = Qs
        self.omega_0 = 2 * np.pi * General_parameters.beta_rel_program * c / self.ring_circumference
        self.omega_s = self.Qs * self.omega_0
        
        self.dQs = 2 * np.pi * self.Qs
        self.cosdQs = np.cos(self.dQs)
        self.sindQs = np.sin(self.dQs)
        
        self.counter = General_parameters.counter

    def track(self, beam):

        z0 = beam.z
        delta0 = beam.delta

        beam.z = z0 * self.cosdQs - self.eta[self.counter[0]] * c / self.omega_s[self.counter[0]] * delta0 * self.sindQs
        beam.delta = delta0 * self.cosdQs + self.omega_s[self.counter[0]] / self.eta[self.counter[0]] / c * z0 * self.sindQs
        
        # Updating the beam synchronous momentum
        beam.beta_rel = self.beta_rel_program[self.counter[0] + 1]
        beam.gamma_rel = self.gamma_rel_program[self.counter[0] + 1]
        beam.energy = self.energy_program[self.counter[0] + 1]
        beam.momentum = self.momentum_program[self.counter[0] + 1]

