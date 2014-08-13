# Example input for longitudinal simulation with RF noise
# No intensity effects


import time 
import numpy as np

from mpi4py import MPI
from LLRF.RF_noise import *
from input_parameters.general_parameters import *
from input_parameters.rf_parameters import *
from trackers.longitudinal_tracker_mpi import *
from setup_parallel.mpi_config import *
from beams.beams import *
from beams.longitudinal_distributions import *
from beams.slices import *
from monitors.monitors import *
from longitudinal_plots.longitudinal_plots import *


# Simulation parameters --------------------------------------------------------
# Bunch parameters
N_b = 1.e9           # Intensity
N_p = 10001          # Macro-particles
tau_0 = 0.4          # Initial bunch length, 4 sigma [ns]

# Machine and RF parameters
C = 26658.883        # Machine circumference [m]
p_s = 450.e9         # Synchronous momentum [eV]
h = 35640            # Harmonic number
V = 6.e6             # RF voltage [eV]
dphi = 0             # Phase modulation/offset
gamma_t = 55.759505  # Transition gamma
alpha = 1./gamma_t/gamma_t        # First order mom. comp. factor

# Tracking details
N_t = 201           # Number of turns to track
dt_plt = 200         # Time steps between plots


# Optional parallel processing -------------------------------------------------
mpi_conf = MPI_Config()


# Pre-processing: RF phase noise -----------------------------------------------
if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
    t0 = time.clock()
    f = np.arange(0, 5.6227612455e+03, 1.12455000e-02)
    spectrum = np.concatenate((1.11100000e-07 * np.ones(4980), np.zeros(495021)))
    noise_t, noise_dphi = Phase_noise(f, spectrum, seed1=1234, seed2=7564).spectrum_to_phase_noise()
    #noise_dphi = mpi_conf.mpi_comm.Bcast(noise_dphi)

    plot_noise_spectrum(f, spectrum, sampling=100)
    plot_phase_noise(noise_t, noise_dphi, sampling=100)
    #plot_phase_noise(noise_t[0:10000], noise_dphi[0:10000], sampling=1)
    print "   Sigma of RF noise is %.4e" %np.std(noise_dphi)
    print "   Time step of RF noise is %.4e" %noise_t[1]
    print ""


# Simulation setup -------------------------------------------------------------
    print "Setting up the simulation..."
    print ""

# Define general parameters
general_params = GeneralParameters(N_t, C, alpha, p_s, 'proton')

# Define RF station parameters and corresponding tracker
rf_params = RFSectionParameters(general_params, 1, 1, h, V, dphi)
#rf_params = RFSectionParameters(general_params, 1, 1, h, V, noise_dphi)

long_tracker = RingAndRFSection(rf_params, mpi_conf=mpi_conf)

if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
    print "General and RF parameters set..."


# Define beam and distribution
beam = Beam(general_params, N_p, N_b)
# Generate new distribution
if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
    longitudinal_gaussian_matched(general_params, rf_params, beam, tau_0, 
                                  unit='ns', reinsertion='on')
    #np.savetxt('initial_long_distr.dat', np.c_[beam.theta, beam.dE], fmt='%.8e')
    # Read in old distribution
    #beam.theta, beam.dE = np.loadtxt('initial_long_distr.dat', unpack=True)

#     if mpi_conf.mpi_size > 1:
#         mpi_conf.mpi_comm.Bcast(beam.theta)
#         mpi_conf.mpi_comm.Bcast(beam.dE)

# Need slices for the Gaussian fit; slice for the first plot
    slice_beam = Slices(100)
    slice_beam.track(beam)


# Define what to save in file
    bunchmonitor = BunchMonitor('output_data', N_t+1, statistics = "Longitudinal", long_gaussian_fit = "On")

if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
    print "Statistics set..."


# Accelerator map
if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
    map_ = [long_tracker] + [slice_beam] # No intensity effects, no aperture limitations
    print "Map set"
    print ""
else:
    map_ = [long_tracker] # No intensity effects, no aperture limitations




# Tracking ---------------------------------------------------------------------
for i in range(N_t):
    t0 = time.clock()
    
    # Save data
    if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
        bunchmonitor.dump(beam, slice_beam)
        if any(np.isnan(beam.theta)):
            print "Detected NaN in beam theta at time step %d" %i    
    
    # Plot has to be done before tracking (at least for cases with separatrix)
    if (i % dt_plt) == 0:
        if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
            print "Outputting at time step %d..." %i
            print "   Beam momentum %.6e eV" %beam.momentum
            print "   Beam gamma %3.3f" %beam.gamma_r
            print "   Beam beta %3.3f" %beam.beta_r
            print "   Beam energy %.6e eV" %beam.energy
            print "   Four-times r.m.s. bunch length %.4e rad" %(4.*beam.sigma_theta)
            print "   Gaussian bunch length %.4e rad" %beam.bl_gauss
            print ""
            # In plots, you can choose following units: rad, ns, m  
            plot_long_phase_space(beam, general_params, rf_params, 0, 0.0001763, -450, 450, separatrix_plot = True)
            plot_bunch_length_evol(beam, 'output_data', general_params, i, unit='ns')
            plot_bunch_length_evol_gaussian(beam, 'output_data', general_params, slice_beam, i, unit='ns')
            plot_beam_profile(i, general_params, slice_beam)
            if i > 100:
                np.savetxt('final_distribution.dat', np.c_[beam.theta, beam.dE], fmt='%.6e')


    # Track
    for m in map_:
        m.track(beam)
    # Define losses according to separatrix and/or longitudinal position
    if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
        beam.losses_separatrix(general_params, rf_params)
        beam.losses_longitudinal_cut(0.28e-4, 0.75e-4)


if mpi_conf.mpi_comm == None or mpi_conf.mpi_rank == 0:
    print "Done!"
    print "Total computation time %.2f s" %(time.clock() - t0)
    print ""





