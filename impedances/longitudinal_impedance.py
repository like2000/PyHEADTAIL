'''
**Module to compute longitudinal intensity effects**

:Authors: **Danilo Quartullo**, **Alexandre Lasheen**, **Hannes Bartosik**
'''

from __future__ import division
import numpy as np
from scipy.constants import c
from numpy.fft import irfft, fftfreq
import math


class InducedVoltageTime(object):
    '''
    *Induced voltage derived from the sum of several wake fields (time domain).
    Note that if there is no acceleration then obviously precalc == 'on', 
    except for the const_charge method where the distances between the 
    slide centers change from turn to turn.
    If there is acceleration and slices.coord == z or theta, then precalc == 'off';
    If slices.coord == tau then precalc == 'on since the wake, at least for
    the analytic formulas presented in the code, doesn't depend on the energy
    of the beam.*
    '''
    
    def __init__(self, GeneralParameters, Slices, wake_source_list):       
        
        #: *Copy of the Slices object in order to access to the profile but 
        #: also some needed information about the Beam.*
        self.slices = Slices
        
        #: *Wake sources inputed as a list (eg: list of BBResonators objects)*
        self.wake_source_list = wake_source_list
        
        #: *Time array of the wake in [s]*
        self.time_array = 0
        
        #: *Total wake array of all sources in* [:math:`\Omega / s`]
        self.total_wake = 0
        
        #: *Induced voltage from the sum of the wake sources in [V]*
        self.induced_voltage = 0
        
        # Pre-processing the wakes
        if self.slices.mode is 'const_charge':
            self.precalc = 'off'
        elif np.sum(np.diff(GeneralParameters.beta_r)) == 0:
            self.precalc = 'on'
        elif self.slices.coord is 'tau':
            self.precalc = 'on'
        else:   
            self.precalc = 'off'
                
        if self.precalc is 'on':
            self.time_array_generation()
            self.sum_wakes(self.time_array)
            
            
    def time_array_generation(self):
        '''
        *Generation of the time array in [s], with respect to the slicing 
        and coordinate type.*
        '''
        
        if self.slices.coord == 'tau':
            self.time_array = self.slices.bins_centers - self.slices.bins_centers[0]
        elif self.slices.coord == 'theta':
            self.time_array = (self.slices.bins_centers - self.slices.bins_centers[0]) * self.slices.Beam.ring_radius / (self.slices.Beam.beta_r * c)
        elif self.slices.coord == 'z':
            self.time_array = (self.slices.bins_centers - self.slices.bins_centers[0]) / (self.slices.Beam.beta_r * c)
    
    
    def sum_wakes(self, time_array):
        '''
        *Summing all the wake contributions in one total wake.*
        '''
        
        self.total_wake = np.zeros(len(time_array))
        for wake_object in self.wake_source_list:
            wake_object.wake_calc(time_array)
            self.total_wake += wake_object.wake
        
           
    def induced_voltage_with_matrix(self, Beam):
        '''
        Method to calculate the induced voltage from wakes with the matrix method;
        note that if slices.coord = z one has to "fix" the usual method, since
        the head and the tail of the bunch are inverted.
        '''
        
        if self.slices.coord == 'tau':
            dtau_matrix = self.slices.bins_centers - \
                            np.transpose([self.slices.bins_centers])
            self.wake_matrix = self.sum_wakes(dtau_matrix, self.wake_object_sum)
        elif self.slices.coord == 'z':
            dtau_matrix = (np.transpose([self.slices.bins_centers]) - \
                           self.slices.bins_centers) / (Beam.beta_r * c)
            self.wake_matrix = self.sum_wakes(dtau_matrix, self.wake_object_sum)
        elif self.slices.coord == 'theta':
            dtau_matrix = Beam.ring_radius / (Beam.beta_r * c) * \
            (self.slices.bins_centers - np.transpose([self.slices.bins_centers])) 
            self.wake_matrix = self.sum_wakes(dtau_matrix, self.wake_object_sum)
        
        return - Beam.charge * Beam.intensity / Beam.n_macroparticles * \
                np.dot(self.slices.n_macroparticles, self.wake_matrix)
    
    
    def induced_voltage_with_convolv(self, Beam): 
        '''
        Method to calculate the induced voltage from wakes with convolution;
        note that if slices.coord = z one has to "fix" the usual method, since
        the head and the tail of the bunch are inverted.
        '''
        
        if self.precalc is 'off':
            self.time_array_generation()
            self.sum_wakes(self.time_array)
            
        if self.slices.coord is 'tau':
            self.induced_voltage = - Beam.charge * Beam.intensity / Beam.n_macroparticles * np.convolve(self.total_wake, self.slices.n_macroparticles)[0:len(self.total_wake)]
        
#         if self.precalc == 'off':
#             
#             if self.slices.coord == 'tau':
#                 dtau = self.slices.bins_centers - self.slices.bins_centers[0]
#                 self.wake_array = self.sum_wakes(dtau, self.wake_sum)
#             
#             elif self.slices.coord == 'z':
#                 dtau = (self.slices.bins_centers - self.slices.bins_centers[0])\
#                        /(Beam.beta_r * c)
#                 self.wake_array = self.sum_wakes(dtau, self.wake_sum)
#                 reversed_array = self.wake_array[::-1]
#                 return - Beam.charge * Beam.intensity / Beam.n_macroparticles * \
#                     np.convolve(reversed_array, self.slices.n_macroparticles)[(len(reversed_array) - 1):] 
#             
#             elif self.slices.coord == 'theta':
#                 dtau = (self.slices.bins_centers - self.slices.bins_centers[0]) \
#                        * Beam.ring_radius / (Beam.beta_r * c)
#                 self.wake_array = self.sum_wakes(dtau, self.wake_sum)
#         
#         if self.precalc == 'on' and self.slices.coord == 'z':
#                 reversed_array = self.wake_array[::-1]
#                 return - Beam.charge * Beam.intensity / Beam.n_macroparticles * \
#                     np.convolve(reversed_array, self.slices.n_macroparticles)[(len(reversed_array) - 1):]  
#         
#         return - Beam.charge * Beam.intensity / Beam.n_macroparticles * \
#             np.convolve(self.wake_array, self.slices.n_macroparticles)[0:len(self.wake_array)]
            
            
        def track(self, Beam):
            '''
            Note that if slices.mode = 'const_charge' one is obliged to use the
            matrix method for the calculation of the induced voltage; otherwise
            update_with_interpolation is faster.
            '''
            
            if self.slices.mode == 'const_charge':
                self.ind_vol = self.induced_voltage_with_matrix(Beam)
            else:
                self.ind_vol = self.induced_voltage_with_convolv(Beam)
            
            update_with_interpolation(Beam, self.ind_vol, self.slices)
    
    
class InducedVoltageFreq(object):
    '''
    *Induced voltage derived from the sum of several impedances.*
    '''
    
    def __init__(self, slices, impedance_sum, frequency_step, 
                 sum_slopes_from_induc_imp = None, deriv_mode = 2, mode = 'only_spectrum'):       
        '''
        Constructor
        
        Frequency_step is equal to 1/(dist_centers * n) where dist_centers is
        the distance between the centers of two consecutive slides and (n/2 + 1)
        is the number of sampling points for the frequency array; see the 
        frequency_array method.
        Sum_slopes_from_induc_imp is the sum of all the inductances derived from
        all the inductive impedance, included space charge; see in addition the
        ind_vol_derivative method.
        '''
        self.slices = slices
        self.impedance_sum = impedance_sum
        self.frequency_step = frequency_step
        self.sum_slopes_from_induc_imp = sum_slopes_from_induc_imp
        self.deriv_mode = deriv_mode
        self.mode = mode
        
        if self.mode != 'only_derivative':
            
            self.frequency_fft, self.n_sampling_fft = self.frequency_array(slices)
            self.impedance_array = self.sum_impedances(self.frequency_fft, \
                                           self.impedance_sum)
            
    
    def sum_impedances(self, frequency, imped_object_sum):
        
        total_impedance = np.zeros(len(frequency)) + 0j
        for imped_object in self.impedance_sum:
            imp = imped_object.imped_calc(frequency)
            total_impedance += imp
       
        return total_impedance
    
    
    def frequency_array(self, slices):
        '''
        Method to calculate the sampling frequency array through the rfftfreq
        Python command. Since this command is not included in older version
        of Numpy, the similar fftfreq has been used with some fixes to obtain
        the same result as rfftfreq.
        '''
        
        dcenters = self.slices.bins_centers[1] - self.slices.bins_centers[0]
        
        n = int(math.ceil(1 / (self.frequency_step * dcenters) ))
        
        if n/2 + 1 >= slices.n_slices:
            pass
        else:
            print 'Warning! Resolution in frequency domain is too small, \
                you can get an error or a truncated bunch profile'
            
        if n%2 == 1:
            n += 1
        
        rfftfreq = fftfreq(n, dcenters)[0:int(n/2+1)]
        rfftfreq[-1] = - rfftfreq[-1]
        
        return rfftfreq, n
        
    
    def track(self, bunch):
        '''
        Method to calculate the induced voltage through the bunch spectrum, or
        the derivative of profile, or both; these three choices are represented
        by the mode 'only_spectrum', 'only_derivative', 'spectrum + derivative'
        respectively.
        '''
        if self.mode != 'only_derivative':
        
            self.spectrum = self.slices.beam_spectrum(self.n_sampling_fft)
            
            self.ind_vol = - bunch.charge * bunch.intensity / bunch.n_macroparticles \
                * irfft(self.impedance_array * self.spectrum) * self.frequency_fft[1] \
                * 2*(len(self.frequency_fft)-1) 
            
            self.ind_vol = self.ind_vol[0:self.slices.n_slices]
            
            if self.slices.coord == 'tau':
                pass
            elif self.slices.coord == 'theta':
                self.ind_vol *= (bunch.beta_r * c / bunch.ring_radius) ** 2
            elif self.slices.coord == 'z':
                self.ind_vol *= (bunch.beta_r * c) ** 2
                self.ind_vol = self.ind_vol[::-1]
            if self.mode == 'spectrum + derivative':
                self.ind_vol += self.ind_vol_derivative(bunch)
               
        elif self.mode == 'only_derivative':
            
            self.ind_vol = self.ind_vol_derivative(bunch)
            
        update_with_interpolation(bunch, self.ind_vol, self.slices)
        
        
    def ind_vol_derivative(self, bunch):
        '''
        Method to calculate the induced voltage through the derivative of the
        profile; the impedance must be of inductive type.
        '''
        
        ind_vol_deriv = bunch.charge / (2 * np.pi) * bunch.intensity / bunch.n_macroparticles * \
                            self.sum_slopes_from_induc_imp * \
                            self.slices.beam_profile_derivative(self.deriv_mode)[1] / \
                            (self.slices.bins_centers[1] - self.slices.bins_centers[0]) 
        
        if self.slices.coord == 'tau':
            pass
        elif self.slices.coord == 'theta':
            ind_vol_deriv *= (bunch.beta_r * c /  bunch.ring_radius) ** 2
        elif self.slices.coord == 'z':
            ind_vol_deriv *= (bunch.beta_r * c) ** 2
        
        return ind_vol_deriv


def update_without_interpolation(bunch, induced_voltage, slices):
    '''
    *Other method to update the energy of the particles; this method can be used
    only if one has not used slices.mode == const_space_hist for the slicing.
    Maybe this method could be optimised through Cython or trying to avoid
    the for loop.*
    '''
    
    for i in range(0, slices.n_slices):
            
            bunch.dE[slices.first_index_in_bin[i]:
              slices.first_index_in_bin[i+1]] += induced_voltage[i]
    
    
def update_with_interpolation(bunch, induced_voltage, slices):
    '''
    *Method to update the energy of the particles through interpolation of
    the induced voltage. Note that there is a fix to prevent that one neglects
    all the particles situated between the first edge and the first slice center
    and between the last edge and the last slice center*
    '''
    
    temp1 = slices.bins_centers[0]
    temp2 = slices.bins_centers[-1]
    slices.bins_centers[0] = slices.edges[0]
    slices.bins_centers[-1] = slices.edges[-1]
    
    if slices.coord == 'tau':
        
        induced_voltage_interpolated = np.interp(bunch.tau, 
                            slices.bins_centers, induced_voltage, 0, 0)
        
    elif slices.coord == 'z':
        
        induced_voltage_interpolated = np.interp(bunch.z, 
                            slices.bins_centers, induced_voltage, 0, 0)
        
    elif slices.coord == 'theta':
        
        induced_voltage_interpolated = np.interp(bunch.theta, 
                            slices.bins_centers, induced_voltage, 0, 0)
        
    slices.bins_centers[0] = temp1
    slices.bins_centers[-1] = temp2
    bunch.dE += induced_voltage_interpolated
    
    
class InputTable(object):
    '''
    *Intensity effects from impedance and wake tables.
    If this constructor takes just two arguments, then a wake table is passed;
    if it takes three arguments, then an impedance table is passed. Be careful
    that if you input a wake, the input wake for W(t=0) should be already 
    divided by two (beam loading theorem) ; and that if you input impedance, 
    only the positive  frequencies of the impedance is needed (the impedance
    will be assumed to be Hermitian (Real part symmetric and Imaginary part
    antisymmetric).Note that we add the point (f, Z(f)) = (0, 0) to the 
    frequency and impedance arrays derived from the table.*
    '''
    
    def __init__(self, input_1, input_2, input_3 = None):       
        
        if input_3 == None:
            #: *Time array of the wake in [s]*
            self.time_array = input_1
            #: *Wake array in* [:math:`\Omega / s`]
            self.wake_array = input_2
        else:
            #: *Frequency array of the impedance in [Hz]*
            self.frequency_array = input_1
            #: *Real part of impedance in* [:math:`\Omega`]
            self.Re_Z_array = input_2
            #: *Imaginary part of impedance in* [:math:`\Omega`]
            self.Im_Z_array = input_3
            #: *Impedance array in* [:math:`\Omega`]
            self.impedance = self.Re_Z_array + 1j * self.Im_Z_array
            
            if self.frequency_array[0] != 0:
                self.frequency_array = np.hstack((0, self.frequency_array))
                self.Re_Z_array = np.hstack((0, self.Re_Z_array))
                self.Im_Z_array = np.hstack((0, self.Im_Z_array))
    
    
    def wake_calc(self, new_time_array):
        '''
        *The wake is interpolated in order to scale with the new time array.*
        '''
        
        self.wake = np.interp(new_time_array, self.time_array, self.wake_array, 
                           right = 0)
        self.time_array = new_time_array
        
    
    def imped_calc(self, new_frequency_array):
        '''
        *The impedance is interpolated in order to scale with the new frequency
        array.*
        '''
        
        Re_Z = np.interp(new_frequency_array, self.frequency_array, self.Re_Z_array, 
                      right = 0)
        Im_Z = np.interp(new_frequency_array, self.frequency_array, self.Im_Z_array, 
                      right = 0)
        self.frequency_array = new_frequency_array
        self.impedance = Re_Z + 1j * Im_Z
        
    
    
class Resonators(object):
    '''
    *Impedance contribution from resonators, analytic formulas for both wake 
    and impedance. The resonant modes (and the corresponding R and Q) 
    can be inputed as a list in case of several modes.*
    
    *The model is the following:*
    
    .. math::
    
        Z(f) = \\frac{R}{1 + j Q \\left(\\frac{f}{f_r}-\\frac{f_r}{f}\\right)}
        
    .. math::
        
        W(t>0) = 2\\alpha R e^{-\\alpha t}\\left(\\cos{\\bar{\\omega}t} - \\frac{\\alpha}{\\bar{\\omega}}\\sin{\\bar{\\omega}t}\\right)

        W(0) = \\alpha R
        
    .. math::
        
        \\omega_r = 2 \\pi f_r
        
        \\alpha = \\frac{\\omega_r}{2Q}
        
        \\bar{\\omega} = \\sqrt{\\omega_r^2 - \\alpha^2}
        
    '''
    
    def __init__(self, R_S, frequency_R, Q):
        
        #: *Shunt impepdance in* [:math:`\Omega`]
        self.R_S = np.array([R_S]).flatten()
        
        #: *Resonant frequency in [Hz]*
        self.frequency_R = np.array([frequency_R]).flatten()
        
        #: *Resonant angular frequency in [rad/s]*
        self.omega_R = 2 *np.pi * self.frequency_R
        
        #: *Quality factor*
        self.Q = np.array([Q]).flatten()
        
        #: *Number of resonant modes*
        self.n_resonators = len(self.R_S)
        
        #: *Time array of the wake in [s]*
        self.time_array = 0
        
        #: *Wake array in* [:math:`\Omega / s`]
        self.wake = 0
        
        #: *Frequency array of the impedance in [Hz]*
        self.freq_array = 0
        
        #: *Impedance array in* [:math:`\Omega`]
        self.impedance = 0


    def wake_calc(self, time_array):
        '''
        *Wake calculation method as a function of time.*
        '''
        
        self.time_array = time_array
        self.wake = np.zeros(len(self.time_array))
        
        for i in range(0, self.n_resonators):
       
            alpha = self.omega_R[i] / (2 * self.Q[i])
            omega_bar = np.sqrt(self.omega_R[i] ** 2 - alpha ** 2)
            
            self.wake += (np.sign(self.time_array) + 1) * self.R_S[i] * alpha * \
                         np.exp(-alpha * self.time_array) * \
                         (np.cos(omega_bar * self.time_array) - 
                          alpha / omega_bar * np.sin(omega_bar * self.time_array))
    
    
    def imped_calc(self, freq_array):
        '''
        *Impedance calculation method as a function of frequency.*
        '''
        
        self.freq_array = freq_array
        self.impedance = np.zeros(len(self.freq_array)) + 0j
        
        for i in range(0, self.n_resonators):
            
            self.impedance[1:] += self.R_S[i] / (1 + 1j * self.Q[i] * 
                                                 ((self.freq_array[1:] / self.frequency_R[i]) - 
                                                  (self.frequency_R[i] / self.freq_array[1:])))
 
 

class TravelingWaveCavity(object):
    '''
    *Impedance contribution from traveling wave cavities, analytic formulas for 
    both wake and impedance. The resonance modes (and the corresponding R and a) 
    can be inputed as a list in case of several modes.*
    
    *The model is the following:*
    
    .. math::
    
        Z_+(f) = R \\left[\\left(\\frac{\\sin{\\frac{a\\left(f-f_r\\right)}{2}}}{\\frac{a\\left(f-f_r\\right)}{2}}\\right)^2 - 2i \\frac{a\\left(f-f_r\\right) - \\sin{a\\left(f-f_r\\right)}}{\\left(a\\left(f-f_r\\right)\\right)^2}\\right]
        
        Z_-(f) = R \\left[\\left(\\frac{\\sin{\\frac{a\\left(f+f_r\\right)}{2}}}{\\frac{a\\left(f+f_r\\right)}{2}}\\right)^2 - 2i \\frac{a\\left(f+f_r\\right) - \\sin{a\\left(f+f_r\\right)}}{\\left(a\\left(f+f_r\\right)\\right)^2}\\right]
        
        Z = Z_+ + Z_-
        
    .. math::
        
        W(0<t<\\tilde{a}) = \\frac{4R}{\\tilde{a}}\\left(1-\\frac{t}{\\tilde{a}}\\right)\\cos{\\omega_r t} 

        W(0) = \\frac{2R}{\\tilde{a}}
        
    .. math::
        
        a = 2 \\pi \\tilde{a}
        
    '''
    
    def __init__(self, R_S, frequency_R, a_factor):
        
        #: *Shunt impepdance in* [:math:`\Omega`]
        self.R_S = np.array([R_S]).flatten()
        
        #: *Resonant frequency in [Hz]*
        self.frequency_R = np.array([frequency_R]).flatten()
        
        #: *Damping time a in [s]*
        self.a_factor = np.array([a_factor]).flatten()
        
        #: *Number of resonant modes*
        self.n_twc = len(self.R_S)
        
        #: *Time array of the wake in [s]*
        self.time_array = 0
        
        #: *Wake array in* [:math:`\Omega / s`]
        self.wake = 0
        
        #: *Frequency array of the impedance in [Hz]*
        self.freq_array = 0
        
        #: *Impedance array in* [:math:`\Omega`]
        self.impedance = 0
        
    
    def wake_calc(self, time_array):
        '''
        *Wake calculation method as a function of time.*
        '''
        
        self.time_array = time_array
        self.wake = np.zeros(len(self.time_array))
        
        for i in range(0, self.n_twc):
            a_tilde = self.a_factor[i] / (2 * np.pi)
            indexes = np.where(self.time_array <= a_tilde)
            self.wake[indexes] += (np.sign(self.time_array[indexes]) + 1) * 2 * self.R_S[i] / a_tilde * \
                                  (1 - self.time_array[indexes] / a_tilde) * \
                                  np.cos(2 * np.pi * self.frequency_R[i] * self.time_array[indexes])
    
    
    def imped_calc(self, freq_array):
        '''
        *Impedance calculation method as a function of frequency.*
        '''
        
        self.freq_array = freq_array
        self.impedance = np.zeros(len(self.freq_array)) + 0j
        
        for i in range(0, self.n_twc):
            
            Zplus = self.R_S[i] * ((np.sin(self.a_factor[i] / 2 * (self.freq_array - self.frequency_R[i])) / 
                                    (self.a_factor[i] / 2 * (self.freq_array - self.frequency_R[i])))**2 - 
                                   2j*(self.a_factor[i] * (self.freq_array - self.frequency_R[i]) - 
                                       np.sin(self.a_factor[i] * (self.freq_array - self.frequency_R[i]))) / \
                                    (self.a_factor[i] * (self.freq_array - self.frequency_R[i]))**2)
            
            Zminus = self.R_S[i] * ((np.sin(self.a_factor[i] / 2 * (self.freq_array + self.frequency_R[i])) / 
                                    (self.a_factor[i] / 2 * (self.freq_array + self.frequency_R[i])))**2 - 
                                   2j*(self.a_factor[i] * (self.freq_array + self.frequency_R[i]) - 
                                       np.sin(self.a_factor[i] * (self.freq_array + self.frequency_R[i]))) / \
                                    (self.a_factor[i] * (self.freq_array + self.frequency_R[i]))**2)
            
            self.impedance += Zplus + Zminus   
    
    

class InductiveImpedance(object):
    '''
    *Constant imaginary Z/n impedance. This needs to be extended to the
    cases where there is acceleration as the revolution frequency f0 used
    in the calculation of n=f/f0 is changing.*
    '''
    
    def __init__(self, Z_over_n, revolution_frequency):
        
        #: *Imaginary Z/n in* [:math:`\Omega / Hz`]
        self.Z_over_n = Z_over_n
        
        #: *Revolution frequency in [Hz]*
        self.revolution_frequency = revolution_frequency
        
        #: *Frequency array of the impedance in [Hz]*
        self.freq_array = 0
        
        #: *Impedance array in* [:math:`\Omega`]
        self.impedance = 0
        
        
    def imped_calc(self, freq_array):
        '''
        *Impedance calculation method as a function of frequency.*
        '''    
        
        self.freq_array = freq_array
        self.impedance = (self.freq_array / self.revolution_frequency) * \
                         self.Z_over_n * 1j


 
