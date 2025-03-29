## CM1 Numerical Model, release 21.1  (`cm1r21.1`)
24 March 2024

Summary of changes.

* * *

### 1.  New options/features.

#### `cm1r21.0`:

 - Officially added the "LES with Mesoscale Model" configuration, in which 
   an inner part of the domain (ie, an "inner fine mesh") uses large-eddy 
   simulation (LES) and all other parts of the domain use a mesoscale 
   model configuration in which a planetary boundary layer (PBL) 
   parameterization is used.  User must set cm1setup=4, and should turn 
   on eddy recycling.  See param17 and param18 sections of namelist.input 
   (and README.namelist) for more information.  

 - On a restart, users can now interpolate to a different model grid, which 
   can have higher (or lower) resolution.  No flags/options need to be set 
   by users; rather, if the code detects on a restart that the model grid 
   is different than it was in the previous simulation (which wrote the 
   restart file) then all data are interpolated onto the new grid 
   automatically.  See the "Changing resolution on a restart" section of 
   README.restart for more information.  

 - Added the option to use the "two-part" subgrid model with time-averages 
   (rather than only with spatial averages).  Useful for hurricane LES, among 
   other applications with horizontal heterogeneity.  See t2p_avg option in 
   namelist.input. 

 - Added new pre-configured cases in the "config_files" subroutine: 
      - Hurricane simulation using the "LES within mesoscale model" setup
        (hurricane_les_within_mm)
      - High-resolution, small-domain idealized LES with hurricane winds 
        at a coast (les_HurrCoast)  
      - Wind tunnel case with an immersed cube
        (les_ib_windtunnel) 

 - Added many more diagnostics into the azimuthal-averaging code, including: 
   max/min windspeeds; total upward/downward mass flux; fraction of grid 
   points that exceed certain thresholds (eg, 95% relative humidity, 40 dBZ
   reflectivity); hydrometeor fall velocities; PBL tendencies for radial and 
   tangential velocity components; turbulent fluctuations and variances 
   (e.g., <w'v'>, <v'v'>); and more. 
   (azimavg.F)

 - Added "domain location" as an output variable for simulations that use 
   the "adaptive moving domain" option.  (So, users easily see how far a 
   hurricane moved during a simulation, and plot its track, for example.) 

 - Added a namelist option ("outunits") so user can specify whether they 
   want the spatial variables in output files (eg, netcdf, grads) to be in 
   km or meters. 

 - Moved options for the large-scale (domain-average) nudging code to 
   namelist.input; see the param19 section. 

 - Moved options for immersed boundaries to namelist.input; see the param20 
   section.  

 - Moved options for the simple hurricane boundary layer setup to 
   namelist.input; see the param21 section.  

 - Added namelist options to use a gradual increase (i.e., "ramp up") of 
   the subgrid turbulence model at a beginning of LES runs when turbulence 
   has not yet developed.  See "ramp_sgs" and "ramp_time" in README.namelist.  

 - Added some text to the standard output that notifies users when large 
   arrays are being allocated.  Hopefully, this will make it clearer when 
   CM1 crashes because memory limits have been exceeded. 

 - Added the ability to write out fall velocities of hydrometeors from the 
   Thompson microphysics scheme. 

#### `cm1r21.1`:

 - Added the option to use a convection parameterization.  Currently only the 
   "new Tiedtke" scheme is available (cuparam=1).  

 - Added the WSM6 microphysics scheme from "MMM shared physics" (ptype=7).

### 2.  Changes, fixes, modifications, etc.

#### `cm1r21.0`:

 - Updated the eddy recycling code, based on experiences using it for several
   applications.  Among other changes, the eddy recycler now captures/injects 
   perturbations from a time-average, rather than the full fields.  Also, 
   variables related to the eddy recycling code are now included in 
   namelist.input, for convenience (see param18 section of namelist). 

 - Added some diagnostics variables to output files called recy_cap and 
   recy_inj so users can see exactly where the "capture" and "injection" 
   regions are located when using the eddy recycling code. 

 - Fixed a problem with the code hanging when there aren't enough grid points 
   to do domain decomposition.  
   (param.F)

 - Fixed values of wspa and wspan in cm1out_diag files. 
   (domaindiag.F)

 - Fixed a divide-by-zero problem when lsnudge_ramp_time is set to zero. 
   (solve1.F) 

 - Tried to fix a sneaky problem where the adaptive time stepping algorithm 
   will gradually ramp-down the time step to a very, very small value.  
   Hopefully, this bug is now fixed. 
   (misclibs.F)

 - Fixed a problem that would sometimes cause the model to crash when using 
   the isnd=17 option (ie, when ignoring the wind profile in input_sounding 
   files).  
   (base.F)

#### `cm1r21.1`:

 - Updated the Makefile; noted which options are preferred for NCAR's casper 
   and derecho computers; updated default settings for the gfortran and Cray
   fortran compilers.

 - Updated the NSSL microphysics code.  See the README.NSSLmp file for more 
   information.  Thanks to Ted Mansell of NOAA/NSSL. 

 - Updated the P3 microphysics code to the WRFv4.5.2 version.

 - Updated the MYNN PBL and surface-layer codes to the WRFv4.5.2 versions. 

 - Updated the YSU PBL and "revised WRF" surface-layer codes to the "MMM shared
   physics" versions. 

 - Modified the RRTMG radiation codes to use the same values for physical 
   constants (eg, c_p, R, g) as the CM1 solver. 

 - Fixed a problem interpolating ozone to model levels in the RRTMG radiation 
   codes.  (Thanks to Osamu Miyawaki of NSF-NCAR/ASP for pointing out this 
   problem.) 
   (module_ra_rrtmg_lw.F)

 - Fixed a major problem with the apmasscon=1 option (that adjusts the mean 
   pressure to ensure conservation of total mass) when using a large number 
   of grid points (roughly, > 10^9).  Double precision is now used in the 
   calculation to prevent the solver from blowing up with very large domains.
   (solve2.F)

 - Fixed a minor distributed-memory parallelization problem when using TKE 
   advection with the MYNN PBL code. 
   (solve2.F, cm1.F)

 - Fixed a very minor bug with the initialization of viscosity/diffusivity when 
   using DNS and some LES sgs models. 
   (cm1.F)

 - Fixed a very minor bug with the initialization of effective radii for cloud,
   ice, and snow when using RRTMG. 
   (cm1.F) 

 - Made several improvements to the TC-tracking algoritm (do_adapt_move=true).
   It is now based on the centroid of a smoothed surface-pressure field, and 
   the adaptive movement code now does a better job of keeping the TC center 
   in the middle of the domain. 
   (azimavg.F, cm1.F)

 - Fixed several problems with netcdf-format restart files (restart_format=2).
   (writeout_nc.F)

 - Fixed a bug when writing out variables from the two-part subgrid turbulent 
   model (sgsmodel=3,4); now only the actual levels are written (as opposed to the 
   entire domain, where most of the data is undefined/zero).  
   (writeout.F)

 - Fixed several bugs for LES with terrain. 
   (turb.F, turbtend.F, solve1.F)

 - Fixed a minor bug with the vertical diffusion of subgrid TKE for sgsmodel=1,3,4
   when using doimpl=1 (vertically implicit diffusion).
   (turbtend.F)

 - Fixed some bugs associated with shared-memory (OpenMP) parallelization.
   Also, sped-up OpenMP parallelization in a few places. 

 - Fixed a bug with NBA subgrid turbulence (sgsmodel=5,6) when moisture is 
   present. 
   (turb.F)

 - Fixed a few minor bugs with some "domain" LES diagnostics.
   (domaindiag.F, turb.F)

 - Fixed a few bugs with the "LES within mesoscale model" option (cm1setup=4).

 - Modified the code that calculates and writes azimuthally averaged fields 
   so that much less memory is used.  A similar change was made to the iinit=7
   option (modified Rankine initial vortex) to also use less memory. 
   (azimavg.F, cm1.F, init3d.F)

 - Fixed a few bugs associated with interpolating to a new grid on a restart. 
   Unfortunately, this exercise has exposed several flaws with the format of 
   restart files in CM1; an overhaul of restart files is planned for the next 
   major upgrade of CM1. 
