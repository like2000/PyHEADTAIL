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



def matched_from_distribution_density(Beam, FullRingAndRF, distribution_function, 
                                      TotalInducedVoltage = None, bunch_length = None, 
                                      emittance = None, several_potential_wells = True):
    '''
    *Function to generate a beam by inputing the distribution density. To 
    be implemented: iteratively converge towards the chosen bunch length/emittance 
    and add the potential due to intensity effects.
    An error will be raised if there is not a full potential well (2 max 
    and 1 min at least), or if there are several wells (more than 2 max and 
    1 min). The last error can be ignored by the user with the several_potential_wells 
    option (this might lead to generate a doublet bunch in the longitudinal 
    distribution generation). A margin of 5% is applied in order to be able
    to catch the min/max of the potential well that might be on the edge
    of the frame.*
    '''
    
    # Generate potential well
    n_points_grid = 1e5
    FullRingAndRF.potential_well_generation(n_points = n_points_grid, theta_margin_percent = 0.05)
    potential_well_array = FullRingAndRF.potential_well
    theta_coord_array = FullRingAndRF.potential_well_coordinates
    slippage_factor = abs(FullRingAndRF.RingAndRFSection_list[0].eta_0[0])
    
    # Check for the min/max of the potentiel well
    [min_theta_positions, max_theta_positions],[min_potential_values, max_potential_values] = minmax_location(theta_coord_array,potential_well_array)
    
    print [min_theta_positions, max_theta_positions],[min_potential_values, max_potential_values]
    
    import matplotlib.pyplot as plt
    plt.plot(theta_coord_array,potential_well_array)
    plt.show()
        
    # Compute deltaE coordinates
    max_potential = np.max(potential_well_array)
    max_deltaE = np.sqrt(1)
        


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


def minmax_location(x,f):
    '''
    *Function to locate the minima and maxima of the f(x) numerical function.*
    '''
    
    f_derivative = np.diff(f)
    f_derivative_second = np.diff(f_derivative)
    
    warnings.filterwarnings("ignore")
    f_derivative_zeros = np.append(np.where(f_derivative == 0), np.where(f_derivative[1:]/f_derivative[0:-1] < 0))
    min_x_position = (x[f_derivative_zeros[f_derivative_second[f_derivative_zeros]>0] + 1] + x[f_derivative_zeros[f_derivative_second[f_derivative_zeros]>0]])/2
    max_x_position = (x[f_derivative_zeros[f_derivative_second[f_derivative_zeros]<0] + 1] + x[f_derivative_zeros[f_derivative_second[f_derivative_zeros]<0]])/2
    
    min_values = np.interp(min_x_position, x, f)
    max_values = np.interp(max_x_position, x, f)

    warnings.filterwarnings("default")
                                          
    return [min_x_position, max_x_position], [min_values, max_values]


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
    
    


