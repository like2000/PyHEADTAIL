'''
Created on 12.06.2014

:Authors: **Danilo Quartullo**, **Helga Timko**, **Adrian Oeftiger**, **Alexandre Lasheen**
'''

from __future__ import division
import numpy as np
from scipy.constants import c
from warnings import filterwarnings


def eta_tracking(GeneralParameters, delta, counter, index_section = 0):
    '''
    *Depending on the number of entries in GeneralParameters.alpha_array, 
    the slippage factor is calculated depending on the delta of the beam. As eta
    is used in the tracker, it is calculated with the initial momentum at that 
    time step, and the corresponding relativistic beta and gamma. This is the
    definition of the full slippage factor:*
    
    .. math:: 
        \eta = \sum_{i}(\eta_i \, \delta^i)

    '''
    eta = 0
    for i in xrange( GeneralParameters.alpha_array.size ):   # order = len - 1
        eta_i = getattr(GeneralParameters, 'eta' + str(i))[index_section][counter[0]]
        eta  += eta_i * (delta**i)
    return eta


class Kick(object):
    '''
    *The Kick represents the kick(s) by an RF station at a certain position 
    of the ring. The kicks are summed over the different harmonic RF systems 
    in the station. The cavity phase can be shifted by the user via phi_offset.
    The increment in energy is given by the discrete equation of motion:*
    
    .. math::
        \Delta E_{n+1} = \Delta E_n + \sum_{j=0}^{n_{RF}}{V_{j,n}\,\sin{\\left(h_{j,n}\,\\theta + \phi_{j,n}\\right)}}
        
    '''
    
    def __init__(self, GeneralParameters, RFSectionParameters):
        
        self.counter = GeneralParameters.counter
        self.n_rf_systems = RFSectionParameters.n_rf_systems
        self.harmonic = RFSectionParameters.harmonic_number_list
        self.voltage = RFSectionParameters.voltage_program_list
        self.phi_offset = RFSectionParameters.phi_offset_list
        
    def manual_input(self, n_rf_systems, harmonic_number_list, 
                     voltage_program_list, phi_offset_list, counter):
        
        self.counter = counter
        self.n_rf_systems = n_rf_systems
        self.harmonic = harmonic_number_list
        self.voltage = voltage_program_list
        self.phi_offset = phi_offset_list
        
    def track(self, beam):
        for i in range(self.n_rf_systems):
            beam.dE += self.voltage[i][self.counter[0]] * \
                       np.sin(self.harmonic[i][self.counter[0]] * beam.theta + 
                              self.phi_offset[i][self.counter[0]])
    
    
class Kick_acceleration(object):
    
    """Kick_acceleration gives a single accelerating kick to the bunch. 
    The accelerating kick is defined by the change in the design momentum 
    (synchronous momentum). 
    
    The acceleration is assumed to be distributed over the length of the 
    RF station, so the average beta is used in the calculation of the kick."""
    
    def __init__(self, GeneralParameters, p_increment, counter, index_section = 0):
        
        self.counter = counter
        self.momentum_program = GeneralParameters.momentum_program
        self.p_increment = p_increment
        
        if GeneralParameters.beta_rel_program.size == GeneralParameters.n_turns + 1:
            self.beta_rel_program = GeneralParameters.beta_rel_program
            self.gamma_rel_program = GeneralParameters.gamma_rel_program
            self.energy_program = GeneralParameters.energy_program
        else:
            self.beta_rel_program = GeneralParameters.beta_rel_program[index_section]
            self.gamma_rel_program = GeneralParameters.gamma_rel_program[index_section]
            self.energy_program = GeneralParameters.energy_program[index_section]
        
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

    def __init__(self, drift_length, GeneralParameters, solver, index_section = 0):
        
        
        self.index_section = index_section
        self.drift_length = drift_length
        self.GeneralParameters = GeneralParameters
        self.ring_circumference = GeneralParameters.ring_circumference
        self.counter = GeneralParameters.counter
        self.solver = solver
        self.beta_rel_program = GeneralParameters.beta_rel_program[index_section]

                
    def track(self, beam):  

        if self.solver == 'full': 
            beam.theta = self.beta_rel_program[self.counter[0] + 1] / self.beta_rel_program[self.counter[0]] * beam.theta \
                    + 2 * np.pi * (1 / (1 - eta_tracking(self.GeneralParameters, beam.delta, self.counter, self.index_section) \
                    * beam.delta) - 1) * self.drift_length / self.ring_circumference
        elif self.solver == 'simple':
            beam.theta = beam.theta + 2 * np.pi * self.GeneralParameters.eta0[self.index_section][self.counter[0]] * beam.delta * \
                    self.drift_length / self.ring_circumference
        else:
            raise RuntimeError("ERROR: Choice of longitudinal solver not recognized! Aborting...")
        

class Ring_and_RF_section(object):
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
        
    def __init__(self, GeneralParameters, RF_section_parameters, index_section = 0):
        
        self.index_section = index_section

        self.kick = Kick(RF_section_parameters.n_rf_systems, 
                          RF_section_parameters.harmonic_numbers_list, RF_section_parameters.voltage_program_list, 
                          RF_section_parameters.phi_offset_list, GeneralParameters.counter)
        
        solver = 'simple' # TODO : put this as a parameter
        self.drift = Drift(RF_section_parameters.section_length, GeneralParameters, solver, index_section)
        
        if np.sum(RF_section_parameters.p_increment) == 0:
            self.elements = [self.kick] + [self.drift]
        else:
            self.kick_acceleration = Kick_acceleration(GeneralParameters, RF_section_parameters.p_increment, GeneralParameters.counter, index_section)
            self.elements = [self.kick] + [self.kick_acceleration] + [self.drift]
            
        self.phi_s = calc_phi_s(GeneralParameters, RF_section_parameters)

          
    def track(self, beam):
        for longMap in self.elements:
            longMap.track(beam)
    
    
class Full_Ring_and_RF(object):
    '''
    Full ring object, containing the total map of Ring_and_RFstation
    '''
    
    def __init__(self, GeneralParameters, sum_RF_section_parameters):
        
        assert GeneralParameters.ring_circumference == sum_RF_section_parameters.section_length_sum
        
        self.n_sections = sum_RF_section_parameters.total_n_sections #: Passing the number of sections
        
        self.Ring_and_RFstation_list = [] #: List of ring and RF stations
         
        for i in range(self.n_sections):
            self.Ring_and_RFstation_list.append(Ring_and_RF_section(GeneralParameters, sum_RF_section_parameters.RF_section_parameters_list[i], i))
            
    def track(self, beam):
        
        for i in range(self.n_sections):
            self.Ring_and_RFstation_list[i].track(beam)


def calc_phi_s(GeneralParameters, RF_section_parameters, accelerating_systems = 'all', index_section = 0):
    """The synchronous phase calculated from the rate of momentum change.
    Below transition, for decelerating bucket: phi_s is in (-Pi/2,0)
    Below transition, for accelerating bucket: phi_s is in (0,Pi/2)
    Above transition, for accelerating bucket: phi_s is in (Pi/2,Pi)
    Above transition, for decelerating bucket: phi_s is in (Pi,3Pi/2)
    The synchronous phase is calculated at a certain moment.
    Uses beta, energy averaged over the turn."""
    
#     if GeneralParameters.beta_rel_program.size == GeneralParameters.n_turns + 1:
#         beta_rel_program = GeneralParameters.beta_rel_program
#         eta0 = GeneralParameters.eta0
#     else:
    beta_rel_program = GeneralParameters.beta_rel_program[index_section]
    eta0 = GeneralParameters.eta0[index_section]
        
            
    if RF_section_parameters.n_rf_systems == 1:
        V0 = RF_section_parameters.voltage_program_list[0]
        average_beta = (beta_rel_program[1:] + beta_rel_program[0:-1])/2
        
        phi_s = np.arcsin(average_beta * RF_section_parameters.p_increment / V0 )
        
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
        

def hamiltonian(GeneralParameters, RF_section_parameters, theta, dE, delta):
    """Single RF sinusoidal Hamiltonian.
    Uses beta, energy averaged over the turn.
    To be generalized."""
    
    if RF_section_parameters.section_length != GeneralParameters.ring_circumference:
        raise RuntimeError('WARNING : The hamiltonian is not yet properly computed for several sections !!!')
    
    
    if RF_section_parameters.n_rf_systems == 1:
        counter = GeneralParameters.counter[0]
        h0 = RF_section_parameters.harmonic_numbers_list[0][counter]
        V0 = RF_section_parameters.voltage_program_list[0][counter]
        
        c1 = eta_tracking(GeneralParameters, delta, counter) * c * np.pi / (GeneralParameters.ring_circumference * 
             GeneralParameters.beta_rel_program[0][counter] * GeneralParameters.energy_program[0][counter] )
        c2 = c * GeneralParameters.beta_rel_program[0][counter] * V0 / (h0 * GeneralParameters.ring_circumference)
        
        phi_s = calc_phi_s(GeneralParameters, RF_section_parameters)[counter] # TODO : do not recalculate phi_s every time but preprocess it
    
        return c1 * dE**2 + c2 * (np.cos(h0 * theta) - np.cos(phi_s) + 
                                   (h0 * theta - phi_s) * np.sin(phi_s))
        
    else:
        raise RuntimeError('Hamiltonian for multiple RF is not implemeted yet')


def separatrix(GeneralParameters, RF_section_parameters, theta):
    """Single RF sinusoidal separatrix.
    Uses beta, energy averaged over the turn.
    To be generalized."""
    
    if RF_section_parameters.section_length != GeneralParameters.ring_circumference:
        print 'WARNING : The separatrix is not yet properly computed for several sections !!!'
    
    
    if RF_section_parameters.n_rf_systems == 1:
        counter = GeneralParameters.counter[0]
        h0 = RF_section_parameters.harmonic_numbers_list[0][counter]
        V0 = RF_section_parameters.voltage_program_list[0][counter]
        
    else:
        raise RuntimeError('Separatrix for multiple RF is not implemeted yet')

    phi_s = calc_phi_s(GeneralParameters, RF_section_parameters)[counter] # TODO : do not recalculate phi_s every time but preprocess it
     
    filterwarnings('ignore')
    
    beta_average = (GeneralParameters.beta_rel_program[0][counter + 1] + GeneralParameters.beta_rel_program[0][counter]) / 2
    
    energy_average = (GeneralParameters.energy_program[0][counter + 1] + GeneralParameters.energy_program[0][counter]) / 2
    
    eta0_average = (GeneralParameters.eta0[0][counter + 1] + GeneralParameters.eta0[0][counter])/2
     
    separatrix_array = np.sqrt(beta_average**2 * energy_average *
                    V0 / (np.pi * eta0_average * h0) * 
                    (-np.cos(h0 * theta) - np.cos(phi_s) + 
                    (np.pi - phi_s - h0 * theta) * np.sin(phi_s)))
     
    filterwarnings('default')
        
    return separatrix_array



def is_in_separatrix(GeneralParameters, RF_section_parameters, theta, dE, delta):
    """Condition for being inside the separatrix.
    Single RF sinusoidal.
    Uses beta, energy averaged over the turn.
    To be generalized."""
    
    if RF_section_parameters.section_length != GeneralParameters.ring_circumference:
        print 'WARNING : The separatrix is not yet properly computed for several sections !!!'
    
    
    if RF_section_parameters.n_rf_systems == 1:
        counter = GeneralParameters.counter[0]
        h0 = RF_section_parameters.harmonic_numbers_list[0][counter]
        
    else:
        raise RuntimeError('is_in_separatrix for multiple RF is not implemeted yet')
        
    phi_s = calc_phi_s(GeneralParameters, RF_section_parameters)[counter] # TODO : do not recalculate phi_s every time but preprocess it
    
    Hsep = hamiltonian(GeneralParameters, RF_section_parameters, (np.pi - phi_s) / h0, 0, 0) 
    isin = np.fabs(hamiltonian(GeneralParameters, RF_section_parameters, theta, dE, delta)) < np.fabs(Hsep)

    return isin
        
        
    
class LinearMap(object):
    
    '''
    Linear Map represented by a Courant-Snyder transportation matrix.
    self.alpha is the linear momentum compaction factor.
    Qs is forced to be constant.
    '''

    def __init__(self, GeneralParameters, Qs):

        """alpha is the linear momentum compaction factor,
        Qs the synchroton tune."""
        
        self.beta_rel_program = GeneralParameters.beta_rel_program
        self.gamma_rel_program = GeneralParameters.gamma_rel_program
        self.energy_program = GeneralParameters.energy_program
        self.momentum_program = GeneralParameters.momentum_program
        
        self.ring_circumference = GeneralParameters.ring_circumference
        self.eta = GeneralParameters._eta0
        self.Qs = Qs
        self.omega_0 = 2 * np.pi * GeneralParameters.beta_rel_program * c / self.ring_circumference
        self.omega_s = self.Qs * self.omega_0
        
        self.dQs = 2 * np.pi * self.Qs
        self.cosdQs = np.cos(self.dQs)
        self.sindQs = np.sin(self.dQs)
        
        self.counter = GeneralParameters.counter

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

