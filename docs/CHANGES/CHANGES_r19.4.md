<PRE>

 CM1 Numerical Model, release 19.4  (cm1r19.4)
 25 November 2017 

 Summary of changes.

-------------------------------------------------------------
1.  New options/features.

   -------------------
    cm1r19.2:

 - For simulations with Coriolis acceleration (icor=1), added the ability 
   to use a beta plane.  The f plane approximation is still used by 
   default.  See the variable "betaplane" in the param.F file. 

 - Added the ability to turn off the Rayleigh damper near lateral boundaries
   when hrdamp=1.  See variables hrdamp_west, hrdamp_east, hrdamp_south, and 
   hrdamp_north in the param.F file. 

 - Added subgrid-scale tke budget terms and a preliminary version of 
   resolved-scale tke budget terms to the diagnostic turbulence code
   (doturbdiag = .true.).


   -------------------
    cm1r19.3:

 - Added the updraft nudging initialization scheme of Naylor and Gilmore 
   (2012, MWR, pg 3699).  See iinit=12 in init3d.F for more information.

 - The anelastic and incompressible solvers (psolver=4,5) can now be used 
   with open lateral boundary conditions.  These codes have also been
   optimized for efficiency, and the pressure variable "phi" is now included 
   in output files. 


   -------------------
    cm1r19.4:

 - Added the option to either apply the Rayleigh damper to the perturbation 
   from the base state (irdamp,hrdamp = 1, which was the only option in 
   CM1 previously) or to apply Rayleigh damper to the perturbation from the 
   domain-average profiles (irdamp,hrdamp = 2).  The new option is useful 
   for very long simulations (roughly, >10 days). 

 - Added the RICO precipitating cumulus test case as testcase=7.  See 
   README file in the run/config_files/les_ShallowCuPrecip subdirectory
   and/or VanZanten et al (2011, JAMES) for more information. 

 - Added a radiative-convective equilibrium (RCE) test case as testcase=8.
   See README file in the run/config_files/cpm_RadConvEquil subdirectory
   and/or Bretherton et al (2005, JAS) for more information. 

 - Added a single-column version of the hurricane boundary layer test case
   (testcase=6).  See README file in the run/config_files/scm_HurrBoundLayer 
   subdirectory for more information. 

 - When writing budgets variables (output_thbudget, output_ubudget, etc)
   and an odd-ordered advection scheme is used (3rd, 5th, 7th, 9th-order 
   advection) or a WENO advection scheme is used, then the advection 
   tendencies are separated into a non-diffusive component and the component 
   attributable to implicit diffusion.  
   (Note: in some previous versions of CM1, implicit diffusion tendencies 
    were previously written when output_impdiften=1.  These output fields 
    have now been restored as part of the "budget" output variables introduced
    in cm1r19.) 

 - Explicit diffusion tendencies are now added to the budget variables, 
   including tendencies from the sixth-order horizontal diffusion scheme 
   (idiff >=1 with difforder=6). 

 - Added the option to write pressure decomposition variables to output files.  
   Three diagnostic variables are determined: buoyancy pressure perturbation, 
   non-linear dynamic pressure perturbation, and linear dynamic pressure 
   perturbation (based on the base-state wind profiles). 
   Note: this code does not work with distributed-memory parallelization (MPI)
   yet.  
   See output_pdcomp in the namelist and README.namelist for more info.

 - Added the ability to more easily write any output field to a CM1 output 
   file.  Specifically, the arbitrary arrays "out3d" and "out2d" have been 
   added to CM1.  These arrays are allocated at the beginning of a CM1 
   simulation, and have arbitrary sizes based on the "nout3d" and "nout2d" 
   variables which can be set near the top of the cm1.F file.  As long as 
   nout3d and/or nout2d are >= 1, then the contents of these arrays will be 
   added to CM1 output files, and named "out1" and "out2" and "out3" etc. 
   The only thing the users needs to do is:  1) set values for nout3d and 
   nout3d in cm1.F;  2) make sure the out(i,j,k,1) and out(i,j,k,2) etc 
   arrays are filled with the information/variable/diagnostic that the 
   use wants to output.  See also the brief description near the top of 
   the cm1.F file.  

 - Added additional variables to CM1 output files:

      - arbitrary 3d and 2d output arrays (see item above)
      - vgrad (gradient wind) when the axisymmetric model is used (axisymm=1)
      - wspd, the wind speed in the surface layer (including gust) used to 
        calculate surface fluxes for some surface-layer schemes
      - pressure decomposition variables (see "pressure decomposition" item
        above;  user must set output_pdcomp in the namelist;  (Note: does not
        work with distributed memory MPI simulations, yet)

 - Added additional variables to the azimuthal averaging code:

      - satfrac (fraction of saturated gridpoints)
      - vgrad (gradient wind)
      - buoyancy (wrt base state)
      - pressure decomposition fields (eg, buoyancy and dynamic prs perts)
      - vertical accel due to pressure decomposition fields

 - Added many, many new variables to the turbdiag (turbulence diagnostics)
   output files.  See cm1out_turbdiag* files for more information. 

 - Added max/min PBL depth as a "stats" variable.



-------------------------------------------------------------
2.  Changes, fixes, modifications, etc.

   -------------------
    cm1r19.2:

 - Updated the NSSL microphysics code.  See documentation near the top of 
   module_mp_nssl_2mom.F for more information.  Thanks to Ted Mansell of 
   NOAA/NSSL. 
   (module_mp_nssl_2mom.F)

 - Fixed a bug with the goddard microphysics module that was causing certain 
   compilers to crash.  (Does not change results.)  Thanks to Dan Stern of NRL.
   (goddard.F)

 - Fixed a bug with the netcdf-format stats output routine, which was 
   crashing for long simulations.  (Does not change results.)  Thanks to Dan 
   Kirshbaum of McGill University. 
   (writeout_nc.F)

 - Fixed a bug with the RRTMG shortwave radiation code that was causing  
   some long simulations to crash.  Thanks to Stipo Sentic of New Mexico Tech.
   (module_ra_rrtmg_sw.F)

 - Tweaked the Deardroff/Lilly TKE subgrid model for LES (sgsmodel=1) to 
   better handle the subgrid length scale in stable conditions.  Thanks to 
   Xiaoming Shi of the University of California Berkeley. 
   (turb.F) 

 - Fixed a bug with the diagnostic calculation of horizontal pressure gradient
   (output_ubudget = 1 and/or output_vbudget = 1) when using terrain.  Does 
   not affect simulations (only affects diagnostic budget output). 
   (solve.F)

 - Fixed some minor bugs with I/O when the output frequency is not exactly 
   divisible by the time step.  No change on results;  merely affects the 
   times when output is written. 
   (cm1.F)


   -------------------
    cm1r19.3:

 - Updated the Thompson microphysics scheme (ptype=3).  The version is now
   the same as that in WRF3.9.1.
   (thompson.F)

 - Updated the revised surface layer scheme (sfcmodel=3).  The version is now
   the same as that in WRF3.9.1.
   (sfclayrev.F)

 - Fixed a few minor bugs with the compressible-Boussinesq solver (psolver=6)
   and added the pressure variable "phi" to output files.
   (soundcb.F, writeout.F, restart.F)


   -------------------
    cm1r19.4:

 - Added more documentation for the input_sounding option (isnd=7).  See the 
   relevant section of code in base.F or the webpage 

       http://www2.mmm.ucar.edu/people/bryan/cm1/soundings/ 

   for more information.  
   (base.F)

 - Fixed a bug with the compressible-Boussinesq solver (psolver=6) when using
   terrain.
   (soundcb.F)

 - Modified the calculation of saturation vapor pressure for very cold 
   temperatures (< -100 C) to avoid Bad Things.  Only relevant for very
   deep domains (> 30 km).
   (many files)

 - Replaced the old crude calculation of PBL depth with a slightly less crude
   calculation based on critical Richardson number.  The new code prevents 
   ridiculously large values (> 20 km) of PBL depth.  Note1: as before, use 
   this diagnostic code with caution.  Note2: diagnosed PBL depth is 
   unchanged when using the YSU PBL scheme (ipbl=1). 
   (turb.F, sfcphys.F)

 - For calculation of squared moist Brunt-Vaisala frequency (nm array) the 
   saturation condition is changed from simply the presence of qc and/or qi 
   to a condition based on saturation state.  This changed was made primarily 
   because there are regions where the air can very dry (wrt saturation) but 
   can have cloud ice due to sedimentation of qi from above.  This change 
   affects all the LES subgrid turbulence parameterizations, as well as 
   ipbl = 1,2 (YSU and simple CM1 PBL schemes). 
   (turb.F, ysu.F)

 - Physical constants in RRTMG radiation scheme are changed to be consistent 
   with the rest of CM1. 
   (module_ra_rrtmg_lw.F, module_ra_rrtmg_sw.F, radiation_driver.F)

 - OpenMP parallelization (shared memory parallelization) can now be used for 
   the RRTMG radiation scheme. 
   (module_ra_rrtmg_lw.F, module_ra_rrtmg_sw.F, radiation_driver.F)

 - Fixed a bug that didn't allow the simple CM1 PBL scheme (ipbl=2) to be 
   used without moisture (imoist=0).  Thanks to Thomas Kloetzke of the 
   University of Queensland.  
   (param.F)

 - Fixed a bug with netcdf-format restart files and sfcmodel=1.  Thanks to 
   Dan Stern (NRL). 
   (writeout_nc.F)

 - Other minor bug fixes have been applied.  
   (solve.F, turb.F, cm1.F, etc)




</PRE>
