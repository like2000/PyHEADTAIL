Implementation in this branch (intensity effects)
=============================

- Major (feature implementation that may need further review and agreement):
  - Deleted transverse calculations from slices (need agreement with transverse development)
  - Reorganisation of the Slices module
  	- Constant_space_histogram is not faster than constant_space (and as it gives less info, can be removed)
  	- Constant_space is now the reference (and is constant frame also)
  	- Constant_charge is working, but the losses are not taken into account (the frame can diverge...)
  	- Gaussian fit inside the slice module (and the track method updates the bunch length in Beams class)
  
- Normal (small changes that are not transparent, the example main files should be adapted accordingly):
  - None
  
- Minor (can be implemented in a transparent way):
  - Corrected cut_left and cut_right calculation for n_sigma (divided by 2)
  - Added in the beam class if longitudinal/transverse is defined (only definition, not used for the moment)
  - Documentation:
  	- Slices module done
  
  
Implementation in this branch (longitudinal tracking)
=============================

- Major (feature implementation that may need further review and agreement):
  - None
  
- Normal (small changes that are not transparent, the example main files should be adapted accordingly):
  - PEP8 changes:
    - Changed LLRF module name to llrf
    - Changed GeneralParameters.T0 to t_rev + correcting the calculation method
    - Changed RFSectionParameters.sno to section_index
  - Put section_number as an option in RFSectionParameters input + changing its name to section_index
  - Optimization of the longitudinal_tracker.RingAndRFSection object
  	- Moved eta_tracking from rf_parameters to longitudinal_tracker (transparent)
  	- Changed documentation in order to refer to the RFSectionParameters documentation (transparent)
  	- Added a method that chooses the solver to be 'simple' or 'full' wrt the input order
  
- Minor (can be implemented in a transparent way):
  - Updated .gitignore
  - Documentation
  	- Small corrections in input_paramters.general_parameters
  	- Small corrections in input_paramters.rf_parameters
  	- Small corrections in trackers.longitudinal_tracker
  - Changed general_parameters to GeneralParameters as an input for RFSectionParameters
  - Changed rf_params to RFSectionParameters as an input for RingAndRFSection
  - Secured the cases where the user input momentum compaction with higher orders than 2



PYHEADTAIL LONGITUDINAL v1.0
==========

Longitudinal version of the CERN PyHeadTail code for the simulation of multi-particle 
beam dynamics with collective effects.

The structure is as follows:

1) for example main files, see __EXAMPLE_MAIN_FILES; contains examples for using
   the longitudinal package with acceleration, several RF stations, etc.;
2) 5 folders reserved for the current members of the "longitudinal team" for
   their main files, input and output data;	
3) the doc folder contains the documentation, type make html into the console 
   from the folder itself, then go to build, html and open the index file; 
   note that you need Latex and dvipng (if not present in the Latex distribution) 
   to be able to see displayed all the math formulas;
4) the various packages which are the basis of the simulation code;
5) this README.md file;
6) a setup file to compile the various cython files present in the 
   cython_functions package; this file should be run before launching any 
   simulation; from the console window type "python setup.py cleanall 
   build_ext --inplace".


VERSION CONTENTS
==========

2014-07-23
v1.1.1 Plotting routines now separated into different files:
       beams.plot_beams.py -> phase space, beam statistics plots
       beams.plot_slices.py -> profile, slice statistics
       impedances.plot_impedance.py -> induced voltage, wakes, impedances
       LLRF.plot_llrf.py -> noise in time and frequency domain

2014-07-22
v1.1   Added method in 'longitudinal_impedance' to calculate through the 
       derivative of the bunch profile the induced voltage derived from 
       inductive impedance.
       Added in 'longitudinal_plots' the feature to delete the 'fig' folder 
       for Linux users together with a method for plotting the evolution of
       the position of the bunch.
       Two new example main files showing the use of intensity effects methods
       have been included in the corresponding folder.
       The doc folder has been updated to take into account all the packages
       of the tree.
       Several bugs, mostly in the 'longitudinal_impedance' script, have been
       fixed.

2014-07-17
v1.0   Longitudinal tracker tested. Works for acceleration and multiple
       RF sections.
       Beams and slices ready for transverse features to be added.
       Basic statistics, monitors, and plotting in longitudinal plane.
       Longitudinal bunch generation options. 
       Simple versions of separatrix/Hamiltonian.
       Longitudinal impedance calculations in frequency and time domain.
       RF noise.

