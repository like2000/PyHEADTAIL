'''
**Module to generate longitudinal distributions**

:Authors: **Danilo Quartullo**, **Helga Timko**, **Alexandre Lasheen**
'''

from __future__ import division
import numpy as np
import warnings
from scipy.constants import c
from trackers.longitudinal_utilities import is_in_separatrix



def matched_from_line_density(Beam, line_density):
    '''
    *Function to generate a beam by inputing the line density. The distribution
    density is then reconstructed with the Abel transform and the particles
    randomly generated.*
    '''
    
    pass



def matched_from_distribution_density(Beam, FullRingAndRF, TotalInducedVoltage = None, bunch_length = None, emittance = None):
    '''
    *Function to generate a beam by inputing the distribution density.*
    '''
    
    pass


def line_density_bunch_population(Beam):
    '''
    *Function to generate particles with respect to an input line density (1D grid).*
    '''
    
    pass


def distribution_density_bunch_population(Beam):
    '''
    *Function to generate particles with respect to an input distribution density
    (2D grid).*
    '''
    
    pass



def longitudinal_bigaussian(GeneralParameters, RFSectionParameters, beam, sigma_x,
                             sigma_y, xunit=None, yunit=None, reinsertion = 'off'):
    '''
    *Method to generate a bigaussian distribution by manually input the sigma values.*
    '''
    
    warnings.filterwarnings("once")
    if GeneralParameters.n_sections > 1:
        warnings.warn("WARNING: longitudinal_bigaussian is not yet properly computed for several sections!")
        
    if RFSectionParameters.n_rf > 1:
        warnings.warn("longitudinal_bigaussian for multiple RF is not yet implemented")
    
    counter = RFSectionParameters.counter[0]
    
    harmonic = RFSectionParameters.harmonic[0,counter]
    energy = RFSectionParameters.energy[counter]
    beta = RFSectionParameters.beta_r[counter]
    
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
    phi_s = RFSectionParameters.phi_s[counter]
    
    
    beam.theta = sigma_theta * np.random.randn(beam.n_macroparticles) \
                        + phi_s/harmonic
    beam.dE = sigma_dE * np.random.randn(beam.n_macroparticles)
    
    if reinsertion is 'on':
    
        itemindex = np.where(is_in_separatrix(GeneralParameters, RFSectionParameters,
                                     beam.theta, beam.dE, beam.delta) == False)[0]
         
        while itemindex.size != 0:
         
            beam.theta[itemindex] = sigma_theta * np.random.randn(itemindex.size) \
                    + phi_s/harmonic
            beam.dE[itemindex] = sigma_dE * np.random.randn(itemindex.size)
            itemindex = np.where(is_in_separatrix(GeneralParameters, 
                                RFSectionParameters, beam.theta, beam.dE, beam.delta) 
                                 == False)[0]

  

def longitudinal_gaussian_matched(GeneralParameters, RFSectionParameters, beam, 
                                  four_sigma_bunch_length, unit=None, reinsertion = 'off'):
    '''
    *Method to generate a bigaussian distribution by inputing the sigma value
    on the position coordinate and compute the sigma in the momentum coordinate
    with respect to the RF input.*
    '''
    
    warnings.filterwarnings("once")
        
    if GeneralParameters.n_sections > 1:
        warnings.warn("WARNING: longitudinal_gaussian_matched is not yet properly computed for several sections!")
        
    if RFSectionParameters.n_rf > 1:
        warnings.warn("longitudinal_gaussian_matched for multiple RF is not yet implemented!")
    
    counter = RFSectionParameters.counter[0]
    harmonic = RFSectionParameters.harmonic[0,counter]
    energy = RFSectionParameters.energy[counter]
    voltage = RFSectionParameters.voltage[0,counter]
    beta = RFSectionParameters.beta_r[counter]
    eta0 = RFSectionParameters.eta_0[counter]
    
    if unit == None or unit == 'rad':
        sigma_theta = four_sigma_bunch_length / 4
    elif unit == 'm':
        sigma_theta = four_sigma_bunch_length / (-4 * GeneralParameters.ring_radius) 
    elif unit == 'ns':       
        sigma_theta = four_sigma_bunch_length * beta * c * \
        0.25e-9 / GeneralParameters.ring_radius
    
    phi_s = RFSectionParameters.phi_s[counter]
  
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
    
        itemindex = np.where(is_in_separatrix(GeneralParameters, RFSectionParameters,
                                     beam.theta, beam.dE, beam.delta) == False)[0]
         
        while itemindex.size != 0:
         
            beam.theta[itemindex] = sigma_theta * np.random.randn(itemindex.size) \
                    + phi_s/harmonic
            beam.dE[itemindex] = sigma_dE * np.random.randn(itemindex.size)
            itemindex = np.where(is_in_separatrix(GeneralParameters, 
                                RFSectionParameters, beam.theta, beam.dE, beam.delta) 
                                 == False)[0]
    
    


