'''
**Module to generate longitudinal distributions**

:Authors: **Danilo Quartullo**, **Helga Timko**, **Alexandre Lasheen**
'''

from __future__ import division
import numpy as np
import warnings
from scipy.constants import c
from trackers.longitudinal_utilities import is_in_separatrix
import matplotlib.pyplot as plt



def matched_from_line_density(Beam, line_density):
    '''
    *Function to generate a beam by inputing the line density. The distribution
    density is then reconstructed with the Abel transform and the particles
    randomly generated.*
    '''
    
    pass



def matched_from_distribution_density(Beam, FullRingAndRF, distribution_options, 
                                      emittance = None, TotalInducedVoltage = None, 
                                      bunch_length = None, 
                                      main_harmonic_option = 'lowest_freq'):
    '''
    *Function to generate a beam by inputing the distribution density (by
    choosing the type of distribution and the emittance). 
    To be implemented: iteratively converge towards the chosen bunch length 
    and add the potential due to intensity effects.
    The potential well is preprocessed to check for the min/max and center
    the frame around the separatrix.
    An error will be raised if there is not a full potential well (2 max 
    and 1 min at least), or if there are several wells (more than 2 max and 
    1 min, this case will be treated in the future).
    A margin of 5% is applied in order to be able to catch the min/max of the 
    potential well that might be on the edge of the frame. 
    The slippage factor should be updated to take the higher orders.
    Outputs should be added in order for the user to check step by step if
    his bunch is going to be well generated. More detailed 'step by step' 
    documentation should be implemented*
    '''
    
    # Initialize variables depending on the accelerator parameters
    slippage_factor = abs(FullRingAndRF.RingAndRFSection_list[0].eta_0[0])
    eom_factor_dE = (np.pi * slippage_factor * c) / \
                    (FullRingAndRF.ring_circumference * Beam.beta_r * Beam.energy)
    
    # Generate potential well
    n_points_potential = int(1e5)
    FullRingAndRF.potential_well_generation(n_points = n_points_potential, 
                                            theta_margin_percent = 0.05, 
                                            main_harmonic_option = main_harmonic_option)
    potential_well_array = FullRingAndRF.potential_well
    theta_coord_array = FullRingAndRF.potential_well_coordinates
    theta_resolution = theta_coord_array[1] - theta_coord_array[0]
        
    # Check for the min/max of the potentiel well
    minmax_positions, minmax_values = minmax_location(theta_coord_array, 
                                                      potential_well_array)
    min_theta_positions = minmax_positions[0]
    max_theta_positions = minmax_positions[1]
    max_potential_values = minmax_values[1]
    n_minima = len(min_theta_positions)
    n_maxima = len(max_theta_positions)
            
    # Process the potential well in order to take a frame around the separatrix
    if n_minima == 0:
        raise RuntimeError('The potential well has no minima...')
    if n_minima > n_maxima and n_maxima == 1:
        raise RuntimeError('The potential well has more minima than maxima, and only one maximum')
    if n_maxima == 0:
        print ('Warning: The maximum of the potential well could not be found... \
                You may reconsider the options to calculate the potential well \
                as the main harmonic is probably not the expected one. \
                You may also increase the percentage of margin to compute \
                the potentiel well. The full potential well will be taken')
    elif n_maxima == 1:
        if min_theta_positions[0] > max_theta_positions[0]:
            saved_indexes = (potential_well_array < max_potential_values[0]) * \
                            (theta_coord_array > max_theta_positions[0])
            theta_coord_array = theta_coord_array[saved_indexes]
            potential_well_array = potential_well_array[saved_indexes]
            if potential_well_array[-1] < potential_well_array[0]:
                raise RuntimeError('The potential well is not well defined. \
                                    You may reconsider the options to calculate \
                                    the potential well as the main harmonic is \
                                    probably not the expected one.')
        else:
            saved_indexes = (potential_well_array < max_potential_values[0]) * \
                            (theta_coord_array < max_theta_positions[0])
            theta_coord_array = theta_coord_array[saved_indexes]
            potential_well_array = potential_well_array[saved_indexes]
            if potential_well_array[-1] > potential_well_array[0]:
                raise RuntimeError('The potential well is not well defined. \
                                    You may reconsider the options to calculate \
                                    the potential well as the main harmonic is \
                                    probably not the expected one.')
    elif n_maxima == 2:
        lower_maximum_value = np.min(max_potential_values)
        higher_maximum_value = np.max(max_potential_values)
        lower_maximum_theta = max_theta_positions[max_potential_values == lower_maximum_value]
        higher_maximum_theta = max_theta_positions[max_potential_values == higher_maximum_value]
        if min_theta_positions[0] > lower_maximum_theta:
            saved_indexes = (potential_well_array < lower_maximum_value) * \
                            (theta_coord_array > lower_maximum_theta) * \
                            (theta_coord_array < higher_maximum_theta)
            theta_coord_array = theta_coord_array[saved_indexes]
            potential_well_array = potential_well_array[saved_indexes]
        else:
            saved_indexes = (potential_well_array < lower_maximum_value) * \
                            (theta_coord_array < lower_maximum_theta) * \
                            (theta_coord_array > higher_maximum_theta)
            theta_coord_array = theta_coord_array[saved_indexes]
            potential_well_array = potential_well_array[saved_indexes]
    elif n_maxima > 2:
        raise RuntimeError('Work in progress, case to be included in the future...')
#         left_max_theta = np.min(max_theta_positions)
#         right_max_theta = np.max(max_theta_positions)
#         saved_indexes = (theta_coord_array > left_max_theta) * (theta_coord_array < right_max_theta)
#         theta_coord_array = theta_coord_array[saved_indexes]
#         potential_well_array = potential_well_array[saved_indexes]
    
    # Potential is shifted to put the minimum on 0
    potential_well_array = potential_well_array - np.min(potential_well_array)
    n_points_potential = len(potential_well_array)
    
    # Compute deltaE frame corresponding to the separatrix
    max_potential = np.max(potential_well_array)
    max_deltaE = np.sqrt(max_potential / eom_factor_dE)

    # Saving the Hamilotian values corresponding to dE=0 (with high resolution
    # to be used in the integral to compute J further)
    H_array_dE0 = potential_well_array
    
    # Initializing the grids by reducing the resolution to a 
    # n_points_grid*n_points_grid frame.
    n_points_grid = int(1e3)
    potential_well_indexes = np.arange(0,n_points_potential)
    grid_indexes = np.arange(0,n_points_grid) * n_points_potential / n_points_grid
    theta_coord_array = np.interp(grid_indexes, potential_well_indexes, theta_coord_array)
    deltaE_coord_array = np.linspace(-max_deltaE, max_deltaE, n_points_grid)
    potential_well_array_low_res = np.interp(grid_indexes, potential_well_indexes, potential_well_array)
    theta_grid, deltaE_grid = np.meshgrid(theta_coord_array, deltaE_coord_array)
    potential_well_grid = np.meshgrid(potential_well_array_low_res, potential_well_array_low_res)[0]
    
    # Computing the action J by integrating the dE trajectories
    J_array_dE0 = np.zeros(n_points_grid)
    
    warnings.filterwarnings("ignore")
    
    for i in range(0, n_points_grid):
        dE_trajectory = np.sqrt((potential_well_array_low_res[i] - H_array_dE0)/eom_factor_dE)
        dE_trajectory[np.isnan(dE_trajectory)] = 0
        J_array_dE0[i] = 2 / (2*np.pi) * np.trapz(dE_trajectory, dx=theta_resolution * 
                                                  FullRingAndRF.ring_radius / 
                                                  (Beam.beta_r * c)) 
        
    warnings.filterwarnings("default")
    
    # Sorting the H and J functions in order to be able to interpolate the function J(H)
    H_array_dE0 = potential_well_array_low_res
    sorted_H_dE0 = H_array_dE0[H_array_dE0.argsort()]
    sorted_J_dE0 = J_array_dE0[H_array_dE0.argsort()]
    
    # Calculating the H and J grid
    H_grid = eom_factor_dE * deltaE_grid**2 + potential_well_grid
    J_grid = np.interp(H_grid, sorted_H_dE0, sorted_J_dE0, left = 0, right = np.inf)

    # Computing the density grid
    option = 'density_from_H'
    if option is 'density_from_F':
        distribution_options['parameters'][0] = distribution_options['parameters'][0] / (2*np.pi)
        density_grid = distribution_density_function(J_grid, distribution_options['type'], distribution_options['parameters'])
    elif option is 'density_from_H':
        distribution_options['parameters'][0] = np.interp(distribution_options['parameters'][0] / (2*np.pi), sorted_J_dE0, sorted_H_dE0)
        density_grid = distribution_density_function(H_grid, distribution_options['type'], distribution_options['parameters'])
    
    # Normalizing the grid
    density_grid = density_grid / np.sum(density_grid)
    
    # Populating the bunch
    indexes = np.random.choice(range(0,np.size(density_grid)), Beam.n_macroparticles, p=density_grid.flatten())
    bunch = np.zeros((2,Beam.n_macroparticles))
    bunch[0,:] = theta_grid.flatten()[indexes] + (np.random.rand(Beam.n_macroparticles) - 0.5) * (theta_coord_array[1]-theta_coord_array[0])
    bunch[1,:] = deltaE_grid.flatten()[indexes] + (np.random.rand(Beam.n_macroparticles) - 0.5) * (deltaE_coord_array[1]-deltaE_coord_array[0])
    
    Beam.theta = bunch[0,:]
    Beam.dE = bunch[1,:]



def distribution_density_function(action_array, dist_type, parameters):
    '''
    *Distribution density (formulas from Laclare).*
    '''
    
    if dist_type is 'parabolic':
        length = parameters[0]
        exponent = parameters[1]
        warnings.filterwarnings("ignore")
        density_function = (1 - action_array / length)**exponent
        warnings.filterwarnings("default")
        density_function[action_array > length] = 0
        return density_function
    
    elif dist_type is 'gaussian':
        length = parameters[0]
        exponent = parameters[1]
        density_function = np.exp(- 2 * action_array / length)
        return density_function
    
    elif dist_type is 'waterbag':
        length = parameters[0]
        density_function = np.ones(action_array.shape)
        density_function[action_array > length] = 0
        return density_function
    


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
    
    


