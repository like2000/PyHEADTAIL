'''
**Module to plot different bunch features **

:Authors: **Helga Timko**, **Danilo Quartullo**

'''

from __future__ import division
import matplotlib.pyplot as plt
import h5py
import numpy as np
from longitudinal_plots.plot_settings import fig_folder


def plot_noise_spectrum(frequency, spectrum, sampling = 1, dirname = 'fig'):
    
    """
    Plot of the phase noise spectrum.
    For large amount of data, use "sampling" to plot a fraction of the data.
    """

    # Directory where longitudinal_plots will be stored
    fig_folder(dirname)
    
    # Plot
    plt.figure(1, figsize=(8,6))
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.set_xlim([0, 300])    
    ax.plot(frequency[::sampling], spectrum[::sampling])
    ax.set_xlabel("Frequency [Hz]")
    params = {'text.usetex': False, 'mathtext.default' : 'sf'}
    plt.rcParams.update(params)
    ax.set_ylabel (r"Noise spectrum [$\frac{rad^2}{Hz}$]")

    # Save figure
    fign = dirname +'/noise_spectrum.png'
    plt.savefig(fign)
    plt.clf()
    
    
def plot_phase_noise(time, dphi, sampling = 1, dirname = 'fig'):
    
    """
    Plot of the phase noise as a function of time.
    For large amount of data, use "sampling" to plot a fraction of the data.
    """

    # Directory where longitudinal_plots will be stored
    fig_folder(dirname)
    
    # Plot
    plt.figure(1, figsize=(8,6))
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(time[::sampling], dphi[::sampling])
    ax.set_xlabel("Time [s]")    
    ax.set_ylabel (r"Phase noise [rad]")

    # Save figure
    fign = dirname +'/phase_noise.png'
    plt.savefig(fign)
    plt.clf()     
    
    
def plot_PL_phase_corr(PhaseLoop, h5file, time_step, output_freq = 1, 
                       dirname = 'fig'):
    
    """
    Plot of the phase noise as a function of time.
    For large amount of data, use "sampling" to plot a fraction of the data.
    """

    # Directory where longitudinal_plots will be stored
    fig_folder(dirname)

    # Load/create data
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq) + 1
    t = output_freq*range(0, ndata + 1)    
    storeddata = h5py.File(h5file + '.h5', 'r')
    dphi = np.array(storeddata["/Bunch/PL_phase_corr"], dtype = np.double)
    
    # Plot
    plt.figure(1, figsize=(8,6))
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, dphi[0:ndata+1],'.')
    ax.set_xlabel(r"No. turns [T$_0$]")    
    ax.set_ylabel (r"PL $\phi$ correction [rad]")

    # Save figure
    fign = dirname +'/PL_phase_corr.png'
    plt.savefig(fign)
    plt.clf()     
           

def plot_PL_freq_corr(PhaseLoop, h5file, time_step, output_freq = 1, 
                      dirname = 'fig'):
    
    """
    Plot of the phase noise as a function of time.
    For large amount of data, use "sampling" to plot a fraction of the data.
    """

    # Directory where longitudinal_plots will be stored
    fig_folder(dirname)

    # Load/create data
    if output_freq < 1:
        output_freq = 1
    ndata = int(time_step/output_freq) + 1
    t = output_freq*range(0, ndata + 1)    
    storeddata = h5py.File(h5file + '.h5', 'r')
    dphi = np.array(storeddata["/Bunch/PL_omegaRF_corr"], dtype = np.double)
    
    # Plot
    plt.figure(1, figsize=(8,6))
    ax = plt.axes([0.15, 0.1, 0.8, 0.8])
    ax.plot(t, dphi[0:ndata+1],'.')
    ax.set_xlabel(r"No. turns [T$_0$]")    
    ax.set_ylabel (r"PL $\omega_{RF}$ correction [1/s]")
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # Save figure
    fign = dirname +'/PL_freq_corr.png'
    plt.savefig(fign)
    plt.clf()     
    
    
    