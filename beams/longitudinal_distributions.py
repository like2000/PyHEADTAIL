'''
Created on 12.06.2014

@author: Danilo Quartullo, Helga Timko, Alexandre Lasheen
'''

from __future__ import division
import numpy as np
from scipy.constants import c
from trackers.longitudinal_utilities import is_in_separatrix


def longitudinal_bigaussian(GeneralParameters, RingAndRFSection, beam, sigma_x,
                             sigma_y, xunit=None, yunit=None, reinsertion = 'off'):
    
    if GeneralParameters.n_sections > 1:
        raise RuntimeError('WARNING : The longitudinal_gaussian_matched is not\
         yet properly computed for several sections !!!')
        
    if RingAndRFSection.n_rf > 1:
        raise RuntimeError('longitudinal_gaussian_matched for multiple RF is \
        not implemeted yet')
    
    counter = RingAndRFSection.counter
    harmonic = RingAndRFSection.harmonic_list[0][counter]
    energy = RingAndRFSection.energy[counter]
    beta = RingAndRFSection.beta_r[counter]
    
    if xunit == None or xunit == 'rad':
        sigma_theta = sigma_x
    elif xunit == 'm':
        sigma_theta = sigma_x / (- beam.ring_radius * harmonic) 
    elif xunit == 'ns':       
        sigma_theta = sigma_x * beta * c * 1.e-9 / beam.ring_radius
        
    if yunit == None or yunit == 'eV':
        sigma_dE = sigma_y
    elif yunit == '1':    
        sigma_dE = sigma_y * beta**2 * energy
        
    
    beam.sigma_theta = sigma_theta
    beam.sigma_dE = sigma_dE
    phi_s = RingAndRFSection.phi_s[counter]
    
    beam.theta = sigma_theta * np.random.randn(beam.n_macroparticles) \
                        + phi_s/harmonic
    beam.dE = sigma_dE * np.random.randn(beam.n_macroparticles)
    
    if reinsertion is 'on':
    
        itemindex = np.where(is_in_separatrix(GeneralParameters, RingAndRFSection,
                                     beam.theta, beam.dE, beam.delta) == False)[0]
         
        while itemindex.size != 0:
         
            beam.theta[itemindex] = sigma_theta * np.random.randn(itemindex.size) \
                    + phi_s/harmonic
            beam.dE[itemindex] = sigma_dE * np.random.randn(itemindex.size)
            itemindex = np.where(is_in_separatrix(GeneralParameters, 
                                RingAndRFSection, beam.theta, beam.dE, beam.delta) 
                                 == False)[0]

  

def longitudinal_gaussian_matched(GeneralParameters, RingAndRFSection, beam, 
                                  four_sigma_bunch_length, unit=None, reinsertion = 'off'):
    
    
    if GeneralParameters.n_sections > 1:
        raise RuntimeError('WARNING : The longitudinal_gaussian_matched is not\
         yet properly computed for several sections !!!')
        
    if RingAndRFSection.n_rf > 1:
        raise RuntimeError('longitudinal_gaussian_matched for multiple RF is \
        not implemeted yet')
    
    counter = RingAndRFSection.counter
    harmonic = RingAndRFSection.harmonic_list[0][counter]
    energy = RingAndRFSection.energy[counter]
    voltage = RingAndRFSection.voltage[0][counter]
    beta = RingAndRFSection.beta_r[counter]
    eta0 = RingAndRFSection.eta_0[counter]
    
    if unit == None or unit == 'rad':
        sigma_theta = four_sigma_bunch_length / 4
    elif unit == 'm':
        sigma_theta = four_sigma_bunch_length / (-4 * GeneralParameters.ring_radius) 
    elif unit == 'ns':       
        sigma_theta = four_sigma_bunch_length * beta * c * \
        0.25e-9 / GeneralParameters.ring_radius
    
    phi_s = RingAndRFSection.phi_s[counter]
  
    phi_b = harmonic*sigma_theta + phi_s
    
    sigma_dE = np.sqrt( voltage * energy * beta**2  
             * (np.cos(phi_b) - np.cos(phi_s) + (phi_b - phi_s) * np.sin(phi_s)) 
             / (np.pi * harmonic * eta0) )
        
    beam.sigma_theta = sigma_theta
    beam.sigma_dE = sigma_dE
    
    beam.theta = sigma_theta * np.random.randn(beam.n_macroparticles) \
                        + phi_s/harmonic
    beam.dE = sigma_dE * np.random.randn(beam.n_macroparticles)
    
    if reinsertion is 'on':
    
        itemindex = np.where(is_in_separatrix(GeneralParameters, RingAndRFSection,
                                     beam.theta, beam.dE, beam.delta) == False)[0]
         
        while itemindex.size != 0:
         
            beam.theta[itemindex] = sigma_theta * np.random.randn(itemindex.size) \
                    + phi_s/harmonic
            beam.dE[itemindex] = sigma_dE * np.random.randn(itemindex.size)
            itemindex = np.where(is_in_separatrix(GeneralParameters, 
                                RingAndRFSection, beam.theta, beam.dE, beam.delta) 
                                 == False)[0]
    
    

