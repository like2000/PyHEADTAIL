'''
**Module to generate longitudinal distributions**

:Authors: **Danilo Quartullo**, **Helga Timko**, **Alexandre Lasheen**, **Juan Esteban Muller**, **Theodoros Argyropoulos**
'''

from __future__ import division
import numpy as np
import warnings
import copy
from scipy.constants import c
from trackers.longitudinal_utilities import is_in_separatrix
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt


def matched_from_line_density(Beam, FullRingAndRF, line_density_theta_coord_input, line_density_input, main_harmonic_option = 'lowest_freq', TotalInducedVoltage = None):
    '''
    *Function to generate a beam by inputing the line density. The distribution
    density is then reconstructed with the Abel transform and the particles
    randomly generated.*
    '''

    # Initialize variables depending on the accelerator parameters
    slippage_factor = abs(FullRingAndRF.RingAndRFSection_list[0].eta_0[0])
    eom_factor_dE = (np.pi * slippage_factor * c) / \
                    (FullRingAndRF.ring_circumference * Beam.beta_r * Beam.energy)
    eom_factor_potential = (Beam.beta_r * c) / (FullRingAndRF.ring_circumference)
     
    # Generate potential well
    n_points_potential = int(1e5)
    FullRingAndRF.potential_well_generation(n_points = n_points_potential, 
                                            theta_margin_percent = 0.05, 
                                            main_harmonic_option = main_harmonic_option)
    potential_well_array = FullRingAndRF.potential_well
    theta_coord_array = FullRingAndRF.potential_well_coordinates
    theta_resolution = theta_coord_array[1] - theta_coord_array[0]
    
    # Normalizing the line density
    line_density = line_density_input / np.sum(line_density_input) * Beam.n_macroparticles
    line_density_theta_coord = line_density_theta_coord_input
    
    # Induced voltage contribution
    induced_voltage_potential = 0
    total_potential = potential_well_array + induced_voltage_potential
    
    # Process the potential well in order to take a frame around the separatrix
    theta_coord_sep, potential_well_sep = potential_well_cut(theta_coord_array, total_potential)
    
    # Centering the line density into the synchronous phase
    n_iterations = 100
    if not TotalInducedVoltage:
        n_iterations = 1
        
    for i in range(0, n_iterations):
        if TotalInducedVoltage:
            pass
        minmax_positions_potential, minmax_values_potential = minmax_location(theta_coord_sep, potential_well_sep)
        minmax_positions_profile, minmax_values_profile = minmax_location(line_density_theta_coord, line_density)
        min_theta_positions_potential = minmax_positions_potential[0]
        max_theta_positions_potential = minmax_positions_potential[1]
        min_potential_values_potential = minmax_values_potential[0]
        max_potential_values_potential = minmax_values_potential[1]
        n_minima_potential = len(min_theta_positions_potential)
        n_maxima_potential = len(max_theta_positions_potential)
        min_theta_positions_profile = minmax_positions_profile[0]
        max_theta_positions_profile = minmax_positions_profile[1]
        min_potential_values_profile = minmax_values_profile[0]
        max_potential_values_profile = minmax_values_profile[1]
        n_minima_profile = len(min_potential_values_profile)
        n_maxima_profile = len(max_theta_positions_profile)
        
        if n_minima_potential > 1 or n_maxima_profile > 1:
            raise RuntimeError('n_minmax not good')
        
        line_density_theta_coord = line_density_theta_coord - (max_theta_positions_profile[0] - min_theta_positions_potential[0])
 



def matched_from_distribution_density(Beam, FullRingAndRF, distribution_options,
                                      main_harmonic_option = 'lowest_freq', 
                                      TotalInducedVoltage = None,
                                      n_iterations_input = 1):
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
    
    if not distribution_options.has_key('exponent'):  
        distribution_options['exponent'] = None
    
    # Initialize variables depending on the accelerator parameters
    slippage_factor = abs(FullRingAndRF.RingAndRFSection_list[0].eta_0[0])
    eom_factor_dE = (np.pi * slippage_factor * c) / \
                    (FullRingAndRF.ring_circumference * Beam.beta_r * Beam.energy)
    eom_factor_potential = np.sign(FullRingAndRF.RingAndRFSection_list[0].eta_0[0]) * (Beam.beta_r * c) / (FullRingAndRF.ring_circumference)
    
    # Generate potential well
    n_points_potential = int(1e5)
    FullRingAndRF.potential_well_generation(n_points = n_points_potential, 
                                            theta_margin_percent = 0.05, 
                                            main_harmonic_option = main_harmonic_option)
    potential_well_array = FullRingAndRF.potential_well
    theta_coord_array = FullRingAndRF.potential_well_coordinates
    theta_resolution = theta_coord_array[1] - theta_coord_array[0]
    
    induced_potential = 0
    total_potential = potential_well_array + induced_potential
    
    n_iterations = n_iterations_input
    if not TotalInducedVoltage:
        n_iterations = 1

        
    for i in range(0, n_iterations):
        
        old_potential = copy.deepcopy(total_potential)
        # Adding the induced potential to the RF potential
        total_potential = potential_well_array + induced_potential
        
        sse = np.sqrt(np.sum((old_potential-total_potential)**2))

        print 'Matching the bunch... (iteration: ' + str(i) + ' and sse: ' + str(sse) +')'
                
        # Process the potential well in order to take a frame around the separatrix
        theta_coord_sep, potential_well_sep = potential_well_cut(theta_coord_array, total_potential)
        
        # Potential is shifted to put the minimum on 0
        potential_well_sep = potential_well_sep - np.min(potential_well_sep)
        n_points_potential = len(potential_well_sep)
        
        # Compute deltaE frame corresponding to the separatrix
        max_potential = np.max(potential_well_sep)
        max_deltaE = np.sqrt(max_potential / eom_factor_dE)
    
        # Saving the Hamilotian values corresponding to dE=0 (with high resolution
        # to be used in the integral to compute J further)
        H_array_dE0 = potential_well_sep
        
        # Initializing the grids by reducing the resolution to a 
        # n_points_grid*n_points_grid frame.
        n_points_grid = int(1e3)
        potential_well_indexes = np.arange(0,n_points_potential)
        grid_indexes = np.arange(0,n_points_grid) * n_points_potential / n_points_grid
        theta_coord_low_res = np.interp(grid_indexes, potential_well_indexes, theta_coord_sep)
        deltaE_coord_array = np.linspace(-max_deltaE, max_deltaE, n_points_grid)
        potential_well_low_res = np.interp(grid_indexes, potential_well_indexes, potential_well_sep)
        theta_grid, deltaE_grid = np.meshgrid(theta_coord_low_res, deltaE_coord_array)
        potential_well_grid = np.meshgrid(potential_well_low_res, potential_well_low_res)[0]
        
        # Computing the action J by integrating the dE trajectories
        J_array_dE0 = np.zeros(n_points_grid)
        
        warnings.filterwarnings("ignore")
        
        for i in range(0, n_points_grid):
            dE_trajectory = np.sqrt((potential_well_low_res[i] - H_array_dE0)/eom_factor_dE)
            dE_trajectory[np.isnan(dE_trajectory)] = 0
            J_array_dE0[i] = 2 / (2*np.pi) * np.trapz(dE_trajectory, dx=theta_resolution * 
                                                      FullRingAndRF.ring_radius / 
                                                      (Beam.beta_r * c)) 

                
        warnings.filterwarnings("default")
        
        # Sorting the H and J functions in order to be able to interpolate the function J(H)
        H_array_dE0 = potential_well_low_res
        sorted_H_dE0 = H_array_dE0[H_array_dE0.argsort()]
        sorted_J_dE0 = J_array_dE0[H_array_dE0.argsort()]
        
        # Calculating the H and J grid
        H_grid = eom_factor_dE * deltaE_grid**2 + potential_well_grid
        J_grid = np.interp(H_grid, sorted_H_dE0, sorted_J_dE0, left = 0, right = np.inf)
        
        # Computing bunch length (4-rms) as a function of H/J
        density_variable_option = distribution_options['density_variable']
        if distribution_options.has_key('bunch_length'):        
            time_low_res = theta_coord_low_res * FullRingAndRF.ring_radius / (Beam.beta_r * c)
            tau = 0.0
            if density_variable_option is 'density_from_J':
                X_low = sorted_J_dE0[0]
                X_hi = sorted_J_dE0[n_points_grid - 1]
            elif density_variable_option is 'density_from_H':
                X_low = sorted_H_dE0[0]
                X_hi = sorted_H_dE0[n_points_grid - 1]

            while np.abs(distribution_options['bunch_length']-tau) > (time_low_res.max() - time_low_res.min()) / n_points_grid / 10:
                X0 = 0.5 * (X_low+X_hi)
                if density_variable_option is 'density_from_J':
                    density_grid = distribution_density_function(J_grid, distribution_options['type'], X0, distribution_options['exponent'])
                elif density_variable_option is 'density_from_H':
                    density_grid = distribution_density_function(H_grid, distribution_options['type'], X0, distribution_options['exponent'])                
                density_grid = density_grid / np.sum(density_grid)
                
                line_density = np.sum(density_grid, axis = 0)
                
                if (line_density>0).any():
                    tau = 4.0*np.sqrt(np.sum((time_low_res-np.sum(line_density*time_low_res)/np.sum(line_density))**2*line_density)/np.sum(line_density))            
                
                if tau >= distribution_options['bunch_length']:
                    X_hi = X0
                else:
                    X_low = X0
                    
            if density_variable_option is 'density_from_J':
                J0 = X0
            elif density_variable_option is 'density_from_H':
                H0 = X0
        
        # Computing the density grid
        if density_variable_option is 'density_from_J':
            if distribution_options.has_key('emittance'):
                J0 = distribution_options['emittance']/ (2*np.pi)
            density_grid = distribution_density_function(J_grid, distribution_options['type'], J0, distribution_options['exponent'])
        elif density_variable_option is 'density_from_H':
            if distribution_options.has_key('emittance'):
                H0 = np.interp(distribution_options['emittance'] / (2*np.pi), sorted_J_dE0, sorted_H_dE0)
            density_grid = distribution_density_function(H_grid, distribution_options['type'], H0, distribution_options['exponent'])
        
        # Normalizing the grid
        density_grid = density_grid / np.sum(density_grid)
        
        # Induced voltage contribution
        if TotalInducedVoltage is not None:
            # Calculating the line density
            line_density = np.sum(density_grid, axis = 0)
            line_density = line_density / np.sum(line_density) * Beam.n_macroparticles

            # Calculating the induced voltage
            induced_voltage_object = copy.deepcopy(TotalInducedVoltage)
                        
            # Inputing new line density
            induced_voltage_object.slices.n_macroparticles = line_density
            induced_voltage_object.slices.bins_centers = theta_coord_low_res * FullRingAndRF.ring_radius / (Beam.beta_r * c)
            induced_voltage_object.slices.edges = np.linspace(induced_voltage_object.slices.bins_centers[0]-(induced_voltage_object.slices.bins_centers[1]-induced_voltage_object.slices.bins_centers[0])/2,induced_voltage_object.slices.bins_centers[-1]+(induced_voltage_object.slices.bins_centers[1]-induced_voltage_object.slices.bins_centers[0])/2,len(induced_voltage_object.slices.bins_centers)+1)
            induced_voltage_object.slices.n_slices = n_points_grid
            induced_voltage_object.slices.fit_option = 'off'
            
            # Re-calculating the sources of wakes/impedances according to this slicing
            induced_voltage_object.reprocess(induced_voltage_object.slices)
            
            # Calculating the induced voltage
            induced_voltage_length_sep = int(np.ceil((theta_coord_array[-1] -  theta_coord_low_res[0]) / (theta_coord_low_res[1] - theta_coord_low_res[0])))
            induced_voltage = induced_voltage_object.induced_voltage_sum(Beam, length = induced_voltage_length_sep)
            theta_induced_voltage = np.linspace(theta_coord_low_res[0], theta_coord_low_res[0] + (induced_voltage_length_sep - 1) * (theta_coord_low_res[1] - theta_coord_low_res[0]), induced_voltage_length_sep)
            
            # Calculating the induced potential
            induced_potential_low_res = - eom_factor_potential * np.insert(cumtrapz(induced_voltage, dx=theta_induced_voltage[1]-theta_induced_voltage[0]),0,0)
            induced_potential = np.interp(theta_coord_array, theta_induced_voltage, induced_potential_low_res)
            
     
    # Populating the bunch
    indexes = np.random.choice(range(0,np.size(density_grid)), Beam.n_macroparticles, p=density_grid.flatten())
     
    Beam.theta = theta_grid.flatten()[indexes] + (np.random.rand(Beam.n_macroparticles) - 0.5) * (theta_coord_low_res[1]-theta_coord_low_res[0])
    Beam.dE = deltaE_grid.flatten()[indexes] + (np.random.rand(Beam.n_macroparticles) - 0.5) * (deltaE_coord_array[1]-deltaE_coord_array[0])
    


def distribution_density_function(action_array, dist_type, length, exponent = None):
    '''
    *Distribution density (formulas from Laclare).*
    '''
    
    if dist_type is 'binomial':
        warnings.filterwarnings("ignore")
        density_function = (1 - action_array / length)**exponent
        warnings.filterwarnings("default")
        density_function[action_array > length] = 0
        return density_function
    
    elif dist_type is 'gaussian':
        density_function = np.exp(- 2 * action_array / length)
        return density_function
    
    elif dist_type is 'waterbag':
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



def potential_well_cut(theta_coord_array, potential_array):
    '''
    *Function to cut the potential well in order to take only the separatrix
    (several cases according to the number of min/max).*
    '''
    
    # Check for the min/max of the potential well
    minmax_positions, minmax_values = minmax_location(theta_coord_array, 
                                                      potential_array)
    min_theta_positions = minmax_positions[0]
    max_theta_positions = minmax_positions[1]
    max_potential_values = minmax_values[1]
    n_minima = len(min_theta_positions)
    n_maxima = len(max_theta_positions)
    
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
            saved_indexes = (potential_array < max_potential_values[0]) * \
                            (theta_coord_array > max_theta_positions[0])
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
            if potential_array[-1] < potential_array[0]:
                raise RuntimeError('The potential well is not well defined. \
                                    You may reconsider the options to calculate \
                                    the potential well as the main harmonic is \
                                    probably not the expected one.')
        else:
            saved_indexes = (potential_array < max_potential_values[0]) * \
                            (theta_coord_array < max_theta_positions[0])
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
            if potential_array[-1] > potential_array[0]:
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
            saved_indexes = (potential_array < lower_maximum_value) * \
                            (theta_coord_array > lower_maximum_theta) * \
                            (theta_coord_array < higher_maximum_theta)
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
        else:
            saved_indexes = (potential_array < lower_maximum_value) * \
                            (theta_coord_array < lower_maximum_theta) * \
                            (theta_coord_array > higher_maximum_theta)
            theta_coord_sep = theta_coord_array[saved_indexes]
            potential_well_sep = potential_array[saved_indexes]
    elif n_maxima > 2:
#             raise RuntimeError('Work in progress, case to be included in the future...')
        left_max_theta = np.min(max_theta_positions)
        right_max_theta = np.max(max_theta_positions)
        saved_indexes = (theta_coord_array > left_max_theta) * (theta_coord_array < right_max_theta)
        theta_coord_sep = theta_coord_array[saved_indexes]
        potential_well_sep = potential_array[saved_indexes]
        
        
    return theta_coord_sep, potential_well_sep



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
    
    


