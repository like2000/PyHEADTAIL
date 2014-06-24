'''
Created on 12.06.2014

@author: Danilo Quartullo, Helga Timko, Alexandre Lasheen
'''

import numpy as np
from scipy.constants import c, e
from trackers.longitudinal_tracker import calc_phi_s

# def stationary_exponential(H, Hmax, H0, bunch):
# 
#     def psi(dz, dp):
#         result = np.exp(H(dz, dp, bunch) / H0) - np.exp(Hmax / H0)
#         return result
# 
#     return psi

# def _match_simple_gaussian_longitudinal(self, beta_z, sigma_z=None, epsn_z=None):
# 
#         if sigma_z and epsn_z:
#             sigma_delta = epsn_z / (4 * np.pi * sigma_z) * e / self.p0
#             if sigma_z / sigma_delta != beta_z:
#                 print '*** WARNING: beam mismatched in bucket. Set synchrotron tune as to obtain beta_z = ', sigma_z / sigma_delta
#         elif not sigma_z and epsn_z:
#             sigma_z = np.sqrt(beta_z * epsn_z / (4 * np.pi) * e / self.p0)
#             sigma_delta = sigma_z / beta_z
#         else:
#             sigma_delta = sigma_z / beta_z
# 
#         self.z *= sigma_z
#         self.delta *= sigma_delta
        
        
def longitudinal_bigaussian(beam, sigma_x, sigma_y, xunit=None, yunit=None):
    
    '''
    To be corrected with the proper variable names...
    '''
    
    if xunit == None or xunit == 'rad':
        sigma_theta = sigma_x
    elif xunit == 'm':
        sigma_theta = sigma_x / (- beam.ring.radius * beam.ring.harmonic[0]) 
    elif xunit == 'ns':       
        sigma_theta = sigma_x * beam.ring.beta_i(beam) * c * 1.e-9 / beam.ring.radius
        
    if yunit == None or yunit == 'eV':
        sigma_dE = sigma_y
    elif yunit == '1':    
        sigma_dE = sigma_y * beam.ring.beta_i(beam)**2 * beam.ring.energy_i(beam)
        
    
    beam.sigma_theta = sigma_theta
    beam.sigma_dE = sigma_dE
    phi_s = beam.ring.calc_phi_s(beam, beam.ring.voltage)

    for i in xrange(beam.n_macroparticles):
        beam.theta[i] = sigma_theta * np.random.randn() \
                        + phi_s/beam.ring.harmonic[0]
        beam.dE[i] = sigma_dE * np.random.randn()
        
        if not beam.ring.is_in_separatrix(beam, beam.theta[i], beam.dE[i], beam.delta[i]): 
            while not beam.ring.is_in_separatrix(beam, beam.theta[i], beam.dE[i], beam.delta[i]):
                beam.theta[i] = sigma_theta * np.random.randn() \
                        + phi_s/beam.ring.harmonic[0]
                beam.dE[i] = sigma_dE * np.random.randn()

  

def longitudinal_gaussian_matched(General_parameters, RF_section_parameters, beam, four_sigma_bunch_length, unit=None):
    
    if RF_section_parameters.section_length != General_parameters.ring_circumference:
        raise RuntimeError('WARNING : The longitudinal_gaussian_matched is not yet properly computed for several sections !!!')
        
    if RF_section_parameters.n_rf_systems != 1:
        raise RuntimeError('longitudinal_gaussian_matched for multiple RF is not implemeted yet')
    
    counter = General_parameters.counter[0]
    harmonic = RF_section_parameters.harmonic_numbers_list[0][counter]
    voltage = RF_section_parameters.voltage_program_list[0][counter]
    energy = General_parameters.energy_program[0][counter]
    beta = General_parameters.beta_rel_program[0][counter]
    eta0 = General_parameters.eta0[0][counter]
            
    if unit == None or unit == 'rad':
        sigma_theta = four_sigma_bunch_length / 4
    elif unit == 'm':
        sigma_theta = four_sigma_bunch_length / (-4 * General_parameters.ring_radius) 
    elif unit == 'ns':       
        sigma_theta = four_sigma_bunch_length * beta * c * 0.25e-9 / General_parameters.ring_radius
        
    phi_s = calc_phi_s(General_parameters, RF_section_parameters)[counter]
    phi_b = harmonic*sigma_theta + phi_s
    
    sigma_dE = np.sqrt( voltage * energy * beta**2  
             * (np.cos(phi_b) - np.cos(phi_s) + (phi_b - phi_s) * np.sin(phi_s)) 
             / (np.pi * harmonic * eta0) )
        
    beam.sigma_theta = sigma_theta
    beam.sigma_dE = sigma_dE
    
    beam.theta = sigma_theta * np.random.randn(beam.n_macroparticles) \
                        + phi_s/harmonic
    beam.dE = sigma_dE * np.random.randn(beam.n_macroparticles)

#     for i in xrange(beam.n_macroparticles):
#         beam.theta[i] = sigma_theta * np.random.randn() \
#                         + phi_s/beam.ring.harmonic[0]
#         beam.dE[i] = sigma_dE * np.random.randn()
#         
#         if not beam.ring.is_in_separatrix(beam, beam.theta[i], beam.dE[i], beam.delta[i]): 
#             while not beam.ring.is_in_separatrix(beam, beam.theta[i], beam.dE[i], beam.delta[i]):
#                 beam.theta[i] = sigma_theta * np.random.randn() \
#                         + phi_s/beam.ring.harmonic[0]
#                 beam.dE[i] = sigma_dE * np.random.randn()
    
    


