Implementation in this branch (longitudinal utilities with intensity effects)
=============================

- Major
    - Implemented the synchrotron_frequency_spread method
- Normal
	- Moved potential_well_cut and minmax_location functions from longitudinal distribution generation to the longitudinal utilities
	- Added some returns in the longitudinal distribution generation to be used in the synchrotron_frequency_spread method
- Minor



Other corrections
=============================

- Input check in rf_parameters (the voltage program should be n_turns)


Implementation in this branch (bunch generation with intensity effects)
=============================

- Major
  - Bunch generation from distribution function implemented
    - Input is distribution type, emittance and bunch length
  - Bunch generation from line density implemented
    - From Abel Transform
    - Input is line density type and bunch length

- Normal
  - Option in beam spectrum calculation in order to compute the frequency_array only wrt the sampling and time step
  - Functions in the induced voltage objects to reprocess the wake/impedance sources according to a new slicing
  - Corrected InputArray object in longitudinal_impedance
  - Added returns in the induced_voltage_generation in order to return induced_voltage after the frame of the slicing
  - Initialization to 0 instead of empty for arrays
  - Initial guess for gaussian fit corrected

- Minor
  - n_macroparticles converted to int in beam class (solve some warnings)
    

Implementation in this branch (intensity effects)
=============================

- Major (feature implementation that may need further review and agreement):
  - Deleted transverse calculations from slices (need discussions with transverse development)
  - Deleted transverse coordinates from beam (need discussions with transverse development)
  - Reorganisation of the slices module
  	- The different coordinates type is redone, a function to convert values from one coordinates type to another is included to reduce the code length
  	- Constant_space_histogram is not faster than constant_space (weird correlation, does not give the same result whether you use intensity effects or not... needs proper profiling...)
  	- Constant_space is now the reference (and is constant frame also)
  	- Constant_charge is working, but the losses are not taken into account (the frame can diverge...)
  	- Gaussian fit inside the slice module (and the track method updates the bunch length in Beams class)
  - Reorganisation of the longitudinal_impedance module
  	- All calculations are now done in tau
    - The impedance coming from and impedance table is assumed to be 0 for higher frequencies
    - The wake_calc in InputTable assumes that the wake begins at t=0
    - The precalculation is always done in InducedVoltageTime unless you use constant_charge slicing
   	  - The input parameters have been changed and the constructor reorganized accordingly
   	- The number of sampling points for the fft is now a power of 2
  
- Normal (small changes that are not transparent, the example main files should be adapted accordingly):
  - PEP8 corrections:
  	- Slices module
  	- Impedance module
  	- Renamed classes/variables
  
- Minor (can be implemented in a transparent way):
  - Corrected cut_left and cut_right calculation for n_sigma (divided by 2)
  - Documentation:
  	- Slices module done (figures and more details might be added)
  	- Impedance (Work in progress)
  	
- Thoughts
  - Smoothing in the slicing can be done (by weighing the particles in the bins)
  -	Better method to convert directly any kind of value from one coordinate to another, just for this values to be called as properties (like in 
  - To be discussed : standard input/output for the functions/objects etc...
  - Should we include the normalized density in the slicing ?
  - To be implemented : varying frame for slicing with constant_space (for cases where you don't mind recalculating all the time the impedance)
  - Wake calculations can be improved by doing pure matrix calculations (this might be useful in case you don't have pre-processing)
  - Should we be able to use only one source of impedance (and have track methods for all sources) ?
  - Inductive impedance to be updated in order to take into account acceleration
  - Formulas to be included in the documentation
  - Pre-processing the induced_voltage (in order to be used directly afterwards for beam generation for example)
  - Include better filtering options
  
  
Implementation in this branch (longitudinal tracking)
=============================

- Major (feature implementation that may need further review and agreement):
  - FullRingAndRF object implementation in order to compute potential well, hamiltonian, separatrix, etc.
  
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
  - Documentation
  	- Small corrections in input_paramters.general_parameters
  	- Small corrections in input_paramters.rf_parameters
  	- Small corrections in trackers.longitudinal_tracker
  - Changed general_parameters to GeneralParameters as an input for RFSectionParameters
  - Changed rf_params to RFSectionParameters as an input for RingAndRFSection
  - Secured the cases where the user input momentum compaction with higher orders than 2
  
- Thoughts:
  - Better input check for GeneralParameters
  

Implementation in this branch (file management)
=============================
- Removed cython functions (obsolete for the moment, will be re-implemented when dedicated functions will be in use)
- Updated .gitignore


