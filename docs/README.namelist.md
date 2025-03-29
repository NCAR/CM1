
 User's Guide to namelist.input file (which is located in the "run" directory)

 CM1 Numerical Model, Release 21.1  (cm1r21.1) 
 Last updated:  24 March 2024

 Entries with "--- NEW ---" are new (or modified) for cm1r21

 NOTE:  some namelist.input files that have been configured for certain
        idealized simulations ... e.g., squall lines, supercells, 
        hurricanes, LES, etc ... are available in the run/config_files 
        subdirectory.

        For mode information about CM1, see the website:  
          http://www2.mmm.ucar.edu/people/bryan/cm1/

-------------------------------------------------------------------------
  param0 section -- most variables are INTEGER

 nx - Total number of grid points in x direction

 ny - Total number of grid points in y direction

 nz - Total number of grid points in z direction

 ppnode - MPI processes per node  (for MPI runs only)

         NOTEs:  - ppnode is used for input/output purposes only.
                 - This is hardware dependent, so check the documentation
                   for your supercomputer.
                 - For NCAR's yellowstone, use ppnode=16.
                 - For NCAR's cheyenne, use ppnode=36.
                 - If you're not certain what to use, just take a guess;
                   this doesn't affect model performance, but it can make
                   parallel I/O a little faster and cleaner (i.e., fewer
                   restart files, for example).

 timeformat - Format for text printout of model integration time:
             1 = seconds
             2 = minutes
             3 = hours
             4 = days

 timestats - 0 = Do not provide provide timing statistics
             1 = Provide timing statistics at end of simulation
             2 = The same as 1, but include time required to complete 
                 each time step

 terrain_flag -  .true.  = With terrain
                 .false. = No terrain (flat lower boundary)

 procfiles -  .true.  = Text printout/config files for every MPI process
              .false. = Only one text printout/config file  (Default)

 --- NEW ---
 outunits  -  units of x,y,z (ie, space dimensions) in output files
                1 = km  (default for most cases)
                2 = meters


-------------------------------------------------------------------------
  param1 section -- enter REAL values

 dx - Horizontal grid spacing in x direction (m).

 dy - Horizontal grid spacing in y direction (m).

 dz - Vertical grid spacing (m).

    NOTE:
       The variables dx,dy,dz are only used when stretch_* = 0.  When 
       stretch_* >= 1, these variables should be set to an approximately 
       average value (to minimize roundoff errors).  See README.stretch 
       for more information.

 dtl - Large time step (s).

       For psolver = 2,3,4,5,6 this time step is limited by the fastest 
       nonacoustic speed.  For thunderstorm simulations, this is usually 
       the maximum vertical velocity.  Otherwise, this would be the 
       propagation speed of gravity waves.  The following is a rough 
       estimate that usually works well for convective storm simulations:
           dtl = min(dx,dy,dz)/67  (rounded to an appropriate value, of course)

       For psolver=1, this time step is limited by the propagation speed of
       sound waves.  dtl of about min(dx,dy,dz)/700 is recommended.

      When using adaptive time-stepping, set dtl to a reasonable
      "target" value (i.e., a value you think would probably be best for your
      simulation).  This will be used as the initial timestep.

 timax - Maximum integration time (s).

 run_time - Integration time (s) to run model from current time.
            NOTE: Ignored if value is less than zero.
            NOTE: Overrides timax.
            (Useful for restarts.  For example, just integrate model 
             "run_time" seconds forward from current time.)

 tapfrq - Frequency of three-dimensional model output (s).
          Output is in cm1out files.

 rstfrq - Frequency to save model restart files (s).  Set to a negative
          number if restart files are not desired.

 statfrq - Frequency for calculating some interesting output.  (seconds)
           Set to negative number to output stats every timestep.
           Output is in cm1out_stats.dat file.
           See param10 section below for the information that can 
           be requested.

 prclfrq - Frequency to output parcel data (s).  Note ... this does not 
           affect the parcel calculations themselves, which are always 
           updated every timestep;  it merely tells the model how 
           frequently to output the information.  Set to a negative number 
           to output parcel data every time step.

-------------------------------------------------------------------------
  param2 section -- enter INTEGER values
     (value in CAPS is recommended, where applicable)

 cm1setup - Overall CM1 setup, base on how turbulence is handled:

      0 = no subgrid turbulence model & no explicit diffusion
            - essentially integrates the Euler equations 
              (adiabatic and inviscid flow)
              (although, diffusion can still occur via numerical methods)
            - NOTE:  ignores sgsmodel, param7 section, etc

      1 = large-eddy simulation (LES)
            - integrates filtered Navier-Stokes equations (ie, LES equations)
            - NOTE:  user must set "sgsmodel" (and related parameters) below

      2 = mesoscale modeling with planetary boundary layer (PBL) parameterization
            - essentially uses Reynolds-averaged Navier-Stokes (RANS) equations 
            - NOTE:  user must set "ipbl" below 

      3 = direct numerical simulation (DNS)
            - integrates Navier-Stokes equations with explicit diffusion and 
              diffusivity terms
            - NOTE:  user must set the parameters in "param7" section below

               --- NEW ---
      4 = LES within mesoscale model
            - runs LES model within an inner fine mesh; runs mesoscale model 
              with PBL parameterization beyond.  See param17 section for 
              settings.



 testcase - Turns on certain simplied physics schemes and/or specified 
            tendency terms and/or special settings for well-documented 
            test cases.

    0  =  Default.  (no special physics, forcings, or configutation)
          - Most users will use testcase=0.

    1  =  Convective Boundary Layer (CBL) using Large Eddy Simulation (LES) 
          - Based on Sullivan and Patton (2011, JAS, pg 2395)
          - See namelist.input in directory run/config_files/les_ConvBoundLayer

    2  =  Sheared Boundary Layer (SBL) using Large Eddy Simulation (LES)
          - Based on Moeng and Sullivan (1994, JAS, pg 999)
          - See namelist.input in directory run/config_files/les_ShearBoundLayer

    3  =  Shallow cumulus clouds using Large Eddy Simulation (LES)
          - Based on Siebesma et al (2003, JAS, pg 1201)
          - See namelist.input in directory run/config_files/les_ShallowCu

    4  =  Nonprecipitating stratocumulus clouds using Large Eddy Simulation (LES)
          - Based on Stevens et al (2005, MWR, pg 1443)
          - See namelist.input in directory run/config_files/les_StratoCuNoPrecip

    5  =  Drizzling stratocumulus clouds using Large Eddy Simulation (LES)
          - Based on Ackerman et al (2009, MWR, pg 1083)
          - See namelist.input in directory run/config_files/les_StratoCuDrizzle

    6  =  Hurricane Boundary Layer (HBL) using Large Eddy Simulation (LES)
          - Based on Bryan et al (2017, BLM, pg 475)
          - See namelist.input in directory run/config_files/les_HurrBoundLayer

          Also: a single-column model (SCM) version using a PBL scheme is 
          available.
          - See namelist.input in directory run/config_files/scm_HurrBoundLayer

    7  =  Precipitating shallow cumulus clouds using Large Eddy Simulation (LES)
          - Based on RICO shallow Cu case (VanZanten et al 2011, JAMES)
          - See namelist.input in run/config_files/les_ShallowCuPrecip

    8  =  Convection Permitting Model (CPM) simulation of Radiative Convective 
          Equilibrium (RCE)
          - Based on Bretherton et al (2005, JAS)
          - NOTE: No diurnal cycle.  (Solar constant is fixed at 650.83 W/m2)
          - See namelist.input in run/config_files/cpm_RadConvEquil 

    9  =  Stable boundary layer using Large Eddy Simulation (LES)
          - Based on Beare et al. (2006, BLM)
          - See namelist.input in run/config_files/les_StableBoundLayer

    10 =   Hurricane boundary layer with LES or single-column modeling 
           with heat and moisture stratification.  Experimental:  uses 
           nudging to help maintain original temperature and moisture 
           profiles.  See "tqnudge" in solve.F for more information and 
           settings. 

    11 =  Convective boundary layer with moisture (but without clouds). 
          - Based on NCAR LES intercomparison case. 
          - See namelist.input in run/config_files/les_ConvPBL_moisture

    12 =  LES, wind tunnel with immersed cube.
          - Based on Martinuzzi and Tropea (1993, JFE). 
          - See run/config_files/les_ib_windtunnel

    14 =  Shallow cumulus convection over land with diurnal cycle. 
          - Based on Brown et al. (2002, QJRMS, 128, p 1075).
          - See namelist.input in run/config_files/les_ShallowCuLand

    15 =  LES, hurricane winds at a coast. 
          - See run/config_files/les_HurrCoast


 adapt_dt - Use adaptive timestep?  (0=no, 1=yes)
            Model automatically adjusts timestep to maintain stability.
        NOTE: - a reasonable value must still be assigned to dtl (above),
                which will be used as the initial timestep

 irst - Is this a restart?  (0=no, 1=yes)

 rstnum - If this is a restart, this variable specifies the number of the
          restart file.  For example, for file cm1out_0000_0002_rst.dat,
          rstnum is 2.

 iconly - Setup initial conditions only?
              1 = creates initial conditions, but does not run model.
              0 = creates initial conditions, and proceeds with integration.

 hadvordrs - Order of horizontal advection scheme for scalars.
 vadvordrs - Order of vertical advection scheme for scalars.
 hadvordrv - Order of horizontal advection scheme for velocities.
 vadvordrv - Order of vertical advection scheme for velocities.
                [Valid options are 2,3,4,5,6,7,8,9,10]

           Odd-ordered schemes have implicit diffusion.  If an odd-ordered
           scheme is used (3,5,7,9) then idiff can be set to 0 (i.e., no
           additional artifical diffusion is typically necessary).

           Even-ordered schemes (2,4,6,8,10) usually require additional artifical 
           diffusion for stability (i.e., idiff=1 is recommended).  Users can use 
           idiff=1 with difforder=6 and a value of kdiff6 between about 0.02-0.24

 advwenos - Advect scalars (except pressure) with WENO scheme?
 advwenov - Advect velocities with WENO scheme?
            0 = no
            1 = yes, apply on every Runge-Kutta step
            2 = yes, apply on final Runge-Kutta step only  (default)

 weno_order - Formulation for the WENO scheme. 
              Valid options are 3, 5, 7, 9

        References:  
        original 3rd and 5th order weno:  Jiang and Shu, 1996, 
        J. Comput. Phys., 126, pg 202

        original 7th and 9th order weno:  Balsara and Shu, 2000,
        J. Comput. Phys., 160, pg 405

              *** CM1 default formulation ***
        5th order with improved smoothness indicators:
        Borges et al, 2008, J. Comput. Phys., 227, pg 3191

 apmasscon - Adjust average pressure perturbation to ensure conservation
             of dry-air mass?    0 = no
                                 1 = yes

             Note:  This option checks the total dry-air mass
                    in the domain and adjusts the domain-average
                    pressure perturbation to ensure conservation. 
                    In general, this option is only needed for long
                    (several days or more) simulations.

 idiff - Include additional artificial diffusion?  (0=no, 1=yes)
         (in addition to any diffusion associated with cm1setup setting)

           For idiff=1, diffusion of all variables.
           For idiff=2, diffusion only applied to winds (u,v,w).

         User must also set difforder and kdiff2 or kdiff6.

 mdiff - When idiff=1 and difforder=6, apply monotonic version of 
         artificial diffusion?

           (0=no, 1=yes)

           Reference:  Xue, 2000, MWR, p 2853.

 difforder - Order of diffusion scheme.  2=second order, 6=sixth order.

           Second order diffusion is not generally recommended.  It is
           only used for certain idealized cases.  kdiff2 must be set
           appropriately when difforder=2.

           Sixth order diffusion is recommended for general use when
           diffusion of small scales (2-6 delta) is needed.  User must
           also set kdiff6 when difforder=6.

           (not available for axisymmetric simulations)

 imoist - Include moisture?  (0=no, 1=yes)

 ipbl - Use a planetary boundary layer (PBL) parameterization?

             0 = no

             1 = Yonsei University (YSU) PBL parameterization
                   Reference:  Hong et al, 2006, MWR, p 2318
                   ("bl_pbl_physics = 1" in WRF)

             2 = simple PBL parameterization (Louis-type scheme)
                 Reference:  Bryan and Rotunno (2009, MWR, pg 1773)

             3 = GFS-EDMF  (as configured in HWRF_v4.0a)
                  Reference:  Hong and Pan, 1996, MWR, pg 2322

             4 = MYNN (Mellor-Yamada-Nakanishi-Niino) level 2.5
                  Reference:  Nakanishi and Niino (2006, BLM)

             5 = MYNN (Mellor-Yamada-Nakanishi-Niino) level 3
                  Reference:  Nakanishi and Niino (2006, BLM)

             6 = MYJ (Mellor-Yamada-Janjic)
 
        Note:  ipbl >= 1 requires cm1setup = 2

 sgsmodel - Subgrid-scale turbulence model for Large Eddy Simulation (LES)
                 Note:  only used for cm1setup = 1

              1 = TKE scheme (eg, Deardorff, 1980, BLM)
                  (previously, this was iturb=1)
              2 = Smagorinsky scheme (eg, Stevens et al, 1999, JAS, see
                  Appendix B, section b, "Equilibrium models")
                  (previously, this was iturb=2)
              3 = sgsmodel=1 (TKE scheme) + Sullivan et al (1994, BLM)
                  version of two-part model
                  (see tunable parameters in param.F ... 
                   search for "2-part turbulence model")
              4 = sgsmodel=1 (TKE scheme) + Bryan (2020, in prep)
                  version of two-part model
                  (see tunable parameters in param.F ... 
                   search for "2-part turbulence model")
              5 = Nonlinear Backscatter and Anisotropy (NBA) model,
                  Deardorff-type TKE version (Mirocha et al. 2010, MWR)
              6 = Nonlinear Backscatter and Anisotropy (NBA) model,
                  Smagorinsky-type version (Mirocha et al. 2010, MWR)


 tconfig - Calculation of turbulence coefficients for sgsmodel = 1 or 2
            1 = horizontal and vertical turbulence coefficients are the
                same;  use this if dx,dy are about equal to dz (default)
            2 = horizontal turbulence coefficient is different from vertical 
                turbulence coefficient;  use this if dx,dy are much greater
                than dz.

 bcturbs - Lower/upper boundary condition for vertical diffusion of all 
           scalars.  (Applies only to sgsmodel=1,2 and ipbl=2)
                 1 = zero flux (default)
                 2 = zero gradient

 horizturb - Horizontal turbulence parameterization 
             (i.e., horizontal Smagorinsky scheme)
             Reference:  Bryan and Rotunno (2009, MWR, pg 1773)
             (note:  previously, this was part of iturb=3)

                 0 = no
                 1 = yes

        Note:  horizturb = 1 requires cm1setup = 2

 doimpl - Vertically implicit calculation for vertical turbulence tendencies 

                 0 = no   (use vertically explicit scheme)
                 1 = YES  (use vertically implicit scheme)   
                    (note:  doimpl=1 was required in cm1r18)

          Default formulation (doimpl=1) is a Crank-Nicholson scheme, which is 
          absolutely stable (i.e., numerically stable regardless of time step).

          For doimpl=0, an explicit scheme is used for vertical turbulence
          tendencies, which can severely limit the time step for simulations
          with small vertical grid spacing. 

 irdamp - Use upper-level Rayleigh damping zone?
             (acts on u,v,w, and theta only)
          (User must set rdalpha and zd below)

                  0 = no
                  1 = yes, damp towards base state
                  2 = yes, damp towards horizontal average
                        (useful for very long simulations, eg, > 10 days)

 hrdamp - Use Rayleigh damping near lateral boundaries (0=NO,  1=yes)
             (acts on u,v,w only)
          (User must set rdalpha and xhd below)

 psolver - Option for pressure solver.

    CM1 DEFAULT:  depends on model grid
        - when dx,dy,dz are approximately equal, use psolver=2
        - when dz is much smaller than dx,dy, use psolver=3

       1 = Compressible equations, integrated explicitly: 
           No time-splitting, no small time steps, no implicit numerics.
           (note: very expensive for weak wind speeds, < 10 m/s)
           (recommended only use when max wind speed is order 100 m/s)

       2 = Compressible equations, Klemp-Wilhelmson time-splitting, explicit:
           Uses K-W split time steps for acoustic modes; uses explicit 
           calculations of acoustic terms in both vertical and horizontal 
           directions.
           (use if dx,dy,dz are approximately equal)

       3 = Compressible, Klemp-Wilhelmson time-splitting, vertically implicit: 
           Uses K-W split time steps for acoustic modes, with a vetically 
           implicit solver, and horizontally explicit calculations.  
           (as in MM5, ARPS, WRF, MPAS)
           (use if dz is much smaller than dx,dy)

       4 = Anelastic solver:  
           Uses the anelastic mass continuity equation.  Pressure is retrieved 
           diagnostically.
           (Note: OpenMP parallelization only; no MPI parallelization, for now.)

       5 = Incompressible solver:  
           Uses the incompressible mass continuity equation.  Pressure is 
           retrieved diagnostically.
           (Note: OpenMP parallelization only; no MPI parallelization, for now.)

       6 = Compressible-Boussinesq:  
           Compressible equation set, using KW time splitting, fully explicit 
           (like psolver=2), but Boussinesq approx made for pressure-gradient
           terms in velocity equations.  Useful for certain types of idealized 
           modeling.  (see, eg, Bryan and Rotunno, 2014, JAS, pg 1126)
           Note: user must set value for "csound" in param.F (otherwise, default 
           value of 300 m/s will be used)

       7 = Modified compressible equations:  
           Equation set from Klemp and Wilhelmson (1978) but with modified value 
           of sound propagation speed.  Useful for certain simulations with low 
           wind speeds (<10 m/s).
           Note: user must set value for "csound" in param.F (otherwise, default 
           value of 300 m/s will be used)


       NOTE:  since cm1r19, users no longer need to set the "nsound" parameter
              for psolver=2,3,6,7   
              (It is now determined adaptively during simulations)


 ptype - Explicit moisture scheme:

             0 = no microphysics (vapor only)

             1 = Kessler scheme (water only)

             2 = NASA-Goddard version of LFO scheme

             3 = Thompson scheme

             4 = Gilmore/Straka/Rasmussen version of LFO scheme

   (default) 5 = Morrison double-moment scheme

             6 = Rotunno-Emanuel (1987) simple water-only scheme

             7 = WSM6        --- NEW ---

            (Note: options 26,27,28 use namelist nssl2mom_params, see below and 
                   README.NSSLmp
                   3-moment option can be activated with nssl_3moment = .true.)
             26 = NSSL 2-moment scheme (graupel-only, no hail); 
                  graupel density predicted
             27 = NSSL 2-moment scheme (graupel and hail); 
                  graupel and hail densities predicted
             28 = NSSL single-moment scheme (graupel-only, similar to ptype=4);
                  fixed graupel density (rho_qh)
             Ice density prediction can be turned off with nssl_density_on  (logical flag, default .true.) 
                 (ptype 26 or 27)


            (Note: P3 = Predicted Particle Property bulk microphysics scheme)
             50 = P3 1-ice category, 1-moment cloud water
             51 = P3 1-ice category plus double-moment cloud water
             52 = P3 2-ice categories plus double-moment cloud water
             53 = P3 1-ice category, 3-moment ice, plus double-moment cloud water

             55 = Jensen's ISHMAEL (Ice-Spheroids Habit Model with Aspect-ratio Evolution)

 nssl_3moment - logical (default = .false.) Works with ptype 26 and 27 to turn on 
              the reflectivity moments for rain and graupel (ptype=26/27) and for 
              hail (ptype=27) Includes option for bin-emulating melting 
              (Mansell et al. 2020; see README.NSSLmp)

 nssl_density_on  - logical flag (default .true.) to toggle graupel/hail density 
                 prediction (ptype 26 or 27)

 ihail - Use hail or graupel for large ice category when ptype=2,5.
          (Goddard-LFO and Morrison schemes only)
             1 = hail
             0 = graupel

 iautoc - Include autoconversion of qc to qr when ptype = 2?  (0=no, 1=yes)
            (Goddard-LFO scheme only)

               --- NEW ---
 cuparam - Convection parameterization:
             0 = no convection parameterization (default)
             1 = new Tiedtke convection parameterization

 icor - Include Coriolis acceleration?  (0=no,  1=yes)
        (If user chooses 1, then fcor must be set below)
        f-plane is assumed.

    NOTE: if icor=1, consider including a large-scale pressure gradient 
          acceleration term (see lspgrad below)

               --- NEW ---
 betaplane - Use beta plane (i.e., Coriolis term is a function of y)?  
             (0=no, 1=yes)
             Caution: not well tested.  May not work with some lateral
             boundary condition options.  

 lspgrad - Apply large-scale pressure gradient acceleration to u and v
           components of velocity.
        0 = no
        1 = yes, based on geostropic balance using base-state wind profiles
          (note:  lspgrad = 1 was called "pertcor" in earlier versions of cm1)
        2 = yes, based on geostropic balance using ug,vg arrays
        3 = yes, based on gradient-wind balance (Bryan et al 2017, BLM)
        4 = yes, specified values (set ulspg, vlspg in base.F)  

 eqtset - equation set for moist microphysics:
        1 = a traditional (approximate) equation set for cloud models
        2 = an energy-conserving equation set that accounts for the
            heat capacity of hydrometeors (Bryan and Fritsch 2002)
            that also conserves mass

           (note:  value is ignored if imoist = 0, because the
            equations are equivalent in a dry environment)

           (not available with ptype=4)

 idiss - Include dissipative heating?  (0=no,  1=yes)

 efall - Include energy fallout term?  (0=no,  1=yes)

 rterm - Include simple relaxation term that mimics atmospheric radiation?
            (0=no,  1=yes)
         (Note:  this is a very simple approach, and is only recommended
          for highly idealized model simulations.  See Rotunno and Emanuel
          1987, JAS, p. 546 for a description of this term.)

 wbc - West lateral boundary condition.

 ebc - East lateral boundary condition.

 sbc - South lateral boundary condition.

 nbc - North lateral boundary condition.

                        where:  1 = periodic
                                2 = open-radiative
                                3 = rigid walls, free slip
                                4 = rigid walls, no slip

 bbc - bottom boundary condition for winds

                        where:  1 = free slip
                                2 = no slip
                                3 = semi-slip (i.e., partial slip)
   (NOTE: for bbc=3, user must also set some options in param12 section below)
          [see variables that mention "bbc=3" below]
          

 tbc - top boundary condition for winds

                        where:  1 = free slip
                                2 = no slip
                                3 = semi-slip (i.e., partial slip)
              (note: tbc=3 requires sfcmodel=1, and uses cnst_znt or cnst_ust)

 irbc - For bc=2, this is the type of radiative scheme to use:
                    1 = Klemp-Wilhelmson (1978) on large steps
                    2 = Klemp-Wilhelmson (1978) on small steps
                    4 = Durran-Klemp (1983) formulation

 roflux - Restrict outward flux?  (0=no, 1=yes)

       When this option is activated, the total outward mass flux at open
       boundary conditions is not allowed to exceed total inward mass flux.
       This is a requirement for the anelastic solver.  For the compressible
       solvers, this scheme helps prevent runaway outward mass flux that can
       cause domain-total mass loss and pressure falls.

 nudgeobc - Nudge winds at inflow boundaries when using open boundary 
            conditions?    (0=no, 1=yes)

       When using open-radiative lateral boundary conditions, this option 
       nudges the horizontal winds toward the base-state fields where there 
       is inflow.  This option is useful for maintaining an inflowing wind 
       profile in long simulations.  (User must set the variable alphobc; 
       see below).

 isnd - Base-state sounding:  1 = Dry adiabatic
                              2 = Dry isothermal
                              3 = Dry, constant dT/dz
   [some variables can be     4 = Saturated neutrally stable (BF02 sounding)
    set in base.F file]       5 = Weisman-Klemp analytic sounding
                              7 = External file (named 'input_sounding')
                                  (see isnd=7 section of base.F for info)
                                  (some soundings are available at 
                                   http://www2.mmm.ucar.edu/people/bryan/cm1)
                                  (Note: wind profile is also obtained from
                                   input_sounding file; iwnd is ignored)
                              8 = Dry, constant d(theta)/dz
                              9 = Dry, constant Brunt-Vaisala frequency
                             10 = Saturated, constant Brunt-Vaisala frequency
                             11 = Saturated, constant equiv. pot. temp.
                             12 = Dry, adiabatic near surface, constant
                                  lapse rate above 
                             13 = Dry, three different layers having constant
                                  N^2 (squared Brunt-Vaisala frequency)
                             14 = Dry, profile for Convective Boundary Layer 
                                  test case.
                             15 = Moist, analytic, based on DYCOMS-II, 
                                  used for stratocumulus test cases.
                             17 = Same as isnd=7, but wind profiles are neglected.
                                  User must set wind profile using the 'iwnd' option.
                                  External file named 'input_sounding' is used,
                                  although columns 4-5 are ignored.
                                  (see isnd=7 section of base.F for more info)
                                  (some soundings are available at 
                                   http://www2.mmm.ucar.edu/people/bryan/cm1)
                             18 = Dry, sharp inversion in middle of profile, 
                                  used for sheared boundary layer test case.
                             19 = Moist analytic profiles based on BOMEX, 
                                  used for shallow cumulus test case
                             20 = Moist analytic profiles based on RICO,
                                  used for precipitating shallow cumulus test case
                             22 = Initial sounding for stable boundary-layer 
                                  test case (testcase=9).
                             23 = Initial sounding for shallow cumulus test case 
                                  over land (testcase=14). 

 iwnd - Base-state wind profile:  (ignored if isnd=7)
                                  0 = zero winds
    [additional variables         1 = RKW-type profile
     need to be set in            2 = Weisman-Klemp supercell
     base.F file]                 3 = multicell
                                  4 = Weisman-Klemp multicell
                                  5 = Dornbrack etal analytic profile
                                  6 = constant wind
                                  8 = constant or linearly decreasing wind profile, 
                                      used for simple hurricane boundary layer case
                                      (see base.F for more details)
                                  9 = linear wind profiles used for shallow 
                                      cumulus test case
                                  10 = linear wind profiles used for drizzling
                                       stratocumulus test case
                                  11 = wind profiles used for RICO precipitating 
                                       shallow cumulus case

 itern - Initial topography specifications.
         User must also set zs array in init_terrain.F:
                    0 = no terrain (zs=0)
                    1 = bell-shaped hill
                    2 = Schaer test case
                    3 = (case from T. Lane and J. Doyle)
                    4 = specified in external GrADS file

 iinit - 3D initialization option:  0 = no perturbation
                                    1 = warm bubble
                                    2 = cold pool
    [additional variables           3 = line of warm bubbles
     need to be set in              4 = initialization for moist benchmark
     init3d.F file]                 5 = cold blob
                                    7 = tropical cyclone (modified Rankine by default)
    See relevant code in            8 = line thermal with random perturbations
    init3d.F file for               9 = forced convergence (Loftus et al 2008)
    more details.                  10 = momentum forcing (Morrison et al 2015)
                                   11 = Skamarock-Klemp IG wave perturbation
                                   12 = updraft nudging (Naylor and Gilmore 2012)

 irandp - Include random potential temperature perturbations in the
          initial conditions?
             (0=no,  1=yes)
          (set magnitude of perturbations in init3d.F)

 ibalance - Specified balance assumption for initial 3D pressure field
            (ignored if iinit=7)

          0 = no balance (initial pressure perturbation is zero everywhere,
                          except for iinit=7) 
          1 = hydrostatic balance (appropriate for small aspect ratios)
          2 = anelastic balance (initial pressure perturbation is the 
              buoyancy pressure perturbation field for an anelastic 
              atmosphere).  (Does not currently work with MPI setup.)

 iorigin - Specifies location of the origin in horizontal space

          1 = At the bottom-left corner of the domain
              (x goes from 0 km to nx*dx km)
              (y goes from 0 km to ny*dy km)
          2 = At the center of the domain
              (x goes from -nx*dx/2 km to +nx*dx/2 km)
              (y goes from -ny*dy/2 km to +ny*dy/2 km)

 axisymm - Run axisymmetric version of model  (0=no,  1=yes)

        (for axisymm=1, ny must be 1, wbc must be 3, and sbc,nbc must be 1)

        (see README.axisymm for more information)

 imove - Move domain at constant speed (0=no, 1=yes)

           For imove=1, user must set umove and vmove.

 iptra - Integrate passive fluid tracer? (0=no, 1=yes)

           User must initialize "pta" array in init3d.F.

 npt - Total number of passive fluid tracers.

 pdtra - Ensure positive-definiteness for tracers?  (0=no, 1=yes)

 iprcl - Integrate passive parcels? (0=no, 1=yes)

           User must initialize "pdata" array in init3d.F.

 nparcels - Total number of parcels.

-------------------------------------------------------------------------
  param3 section -- enter REAL values

 kdiff2 - Diffusion coefficient for difforder=2.  Specified in m^2/s.

 kdiff6 - Diffusion coefficient for difforder=6.  Specified as a 
          fraction of one-dimensional stability.  A value between
          0.02-0.24 is recommended.

 fcor - Coriolis parameter (1/s).

 kdiv - Coefficient for divergence damper.  Value of ~0.1 is recommended.

        This is only used when psolver=2,3.  (The divergence damper is
        an artificial term designed to damp acoustic waves.)

 alph - Off-centering coefficient for vertically implicit acoustic
        solver.  A value of 0.5 is centered-in-time.  Slight forward-in-time 
        bias is recommended.  Default value is 0.60.
        (only used for psolver=3)

 rdalpha - Inverse e-folding time for upper-level Rayleigh damping layer
           (1/s).  Value of about 1/300 is recommended.

 zd - Height above which Rayleigh damping is applied (m).
      (when irdamp = 1)

 xhd - Distance from lateral boundaries where Rayleigh damping is applied (m).
       (when hrdamp = 1) 

 alphobc - Time scale (s) of nudging tendency when using the nudgeobc option. 

 umove - Constant speed for domain translation in x-direction (m/s)
         (for imove = 1)
         (NOTE: for imove=1 and umove not equal to 0.0, ground-relative winds 
          are umove+ua, umove+u3d, etc)

 vmove - Constant speed for domain translation in y-direction (m/s)
         (for imove = 1)
         (NOTE: for imove=1 and vmove not equal to 0.0, ground-relative winds 
          are vmove+va, vmove+v3d, etc)

 v_t - Constant terminal fall velocity of liquid water (m/s) when ptype=6

         When v_t is negative, all liquid water above a small threshold is 
         removed from the domain, ie, pseudoadiabatic thermodynamics are
         used, following Bryan and Rotunno (2009, JAS, pg 3042).

 l_h - Horizontal turbulence length scale (m) used when horizturb=1
        (ie, 2D Smagorinsky)
            Since cm1r18, this is used OVER LAND ONLY
            (see lhref1,lhref2 for settings OVER OCEAN)

 lhref1 - a reference value of l_h (m):  value for surface pressure of 1015 mb
 lhref2 - a reference value of l_h (m):  value for surface pressure of  900 mb

    notes:  - Since cm1r18, the horizontal turbulence length scale for 
              horizturb=1 is a function of surface pressure 
              (over the OCEAN ONLY).
            - This is based on studies of hurricanes (eg, Bryan 2012, MWR).
            - lhref1 and lhref2 define a linear formulation for horizontal
              turbulence length over the ocean for horizturb=1.
            - For gridpoints over land, l_h (above) is used for horizturb=1.
            - For water points above sea level (e.g., lakes) l_h is used.

 l_inf - Asymptotic vertical turbulence length scale (m) (i.e., vertical
         length scale at z = infinity) used for ipbl=2 (simple parameterized 
         turbulence / boundary layer scheme)

 ndcnst - specified cloud droplet concentration for default version of 
          Morrison microphysics scheme (units of cm-3)

        Note: typical value of ndcnst (nt_c) for maritime environments:  100 cm-3
              typical value of ndcnst (nt_c) for continental environments:  300 cm-3

 --- NEW  ---
 nt_c - same as ndcnst, but for Thompson microphysics scheme

     (NOTE: for other microphysics schemes, you will have 
            to change this value manually in the code)

 --- NEW  ---
 csound - speed of sound (m/s) for psolver=6,7
     (Note: should be roughly 5-10 times larger than maximum flow velocity)

 --- NEW  ---
 cstar - propagation speed (m/s) of outward-propagating waves at open 
         boundaries (for irbc=1,2 only)

-------------------------------------------------------------------------
  param11 section - atmospheric radiation

  NOTE:  The parameters in this section ONLY apply to atmospheric radation.
         You do not need to set anything here unless radopt >= 1.

 radopt  -  Use atmospheric radiation code?

               0 = no
               1 = yes, use the NASA-Goddard scheme
               2 = yes, use the RRTMG scheme

       Note:  the NASA-Goddard longwave and shortwave radiation codes
       were adapted from the ARPS model, courtesy of the ARPS/CAPS group 
       at the University of Oklahoma. 

       Note:  the RRTMG code was adapted from the WRF model. 

       (Note:  TIPA option is not implemented in this version of CM1)

        -----
        Note:  for the NASA-Goddard code, the interaction of radiation with 
        clouds is configured consistently for only two microphysics 
        schemes:  the NASA-Goddard LFO scheme (ptype=2) and the Morrison
        microphysics scheme (ptype=5).  A future version of CM1 might pass the 
        proper variables from all microphysics schemes into the radiation 
        code so that consistent calculations are performed.  

        That said, the radiative tendencies should still be reasonable for all 
        ice microphysics schemes, and there is no issue for clear-sky 
        conditions.  (The inconsistency arises only when radiation interacts 
        with water and ice particles, and the radiation scheme needs to be 
        sent information about hydrometeor size and distribution for 
        accurate calculations.) 
        -----

        -----
        Note:  for the RRTMG code, only the Thompson (ptype=3), Morrison 
        (ptype=5), NSSL (ptype=26/27), P3 (50-53), and Jensen ISHMAEL (55) schemes are 
        accurately coupled with the radiation calculations.
        -----

 If radopt >= 1, set the following parameters:
     
 dtrad   -  Time increment (seconds) between calculation of radiation
            tendency.  (Radiative tendencies are held fixed in-between
            calls to the atmospheric radiation subroutine.)

 ctrlat  -  Latitude (applies to entire domain, for now)

 ctrlon  -  Longitude (applies to entire domain, for now)

     NOTE:  because ctrlat and ctrlon are fixed (for now) the radiation 
            scheme is only appropriate for domains having horizontal 
            extent of order 100--1000 km or less

     (FAQ:  Why are lat and lon fixed across the entire domain?  
      It's because George doesn't have time, at the moment, to deal with 
      map projections in CM1.)

 year    -  Year (integer) at start of simulation

 month   -  Month (integer) at start of simulation

 day     -  Day (integer) at start of simulation

 hour    -  Hour (integer) at start of simulation

 minute  -  Minute (integer) at start of simulation

 second  -  Second (integer) at start of simulation

        Yet Another Note:  the radiation schemes uses three important 
        pieces of information from the surface section (param12) below:
        surface temperature, land/water flag, and land-use type.  
        Make sure you have the desired settings for your simulation below 
        (even if you are not using surface fluxes!).  

-------------------------------------------------------------------------
  param12 section:  surface model, ocean model, boundary layer:

      NOTE:  By default, surface conditions are the same everywhere at the 
      initial time.  But, users can define spatially varying initial surface 
      conditions in init_surface.F

 isfcflx - Include surface fluxes of heat and moisture (0=no, 1=yes)

 sfcmodel - Surface model:
     (Specifically, method to calculate surface fluxes and surface stress 
      over land and water)
     (NOTE:  bbc=3 requires sfcmodel >= 1 )
     (NOTE:  set_znt=1, or set_ust=1, or set_flx=1 requires sfcmodel = 1)

      List (see more info below):

            1 = original CM1 formulation

            2 = surface-layer scheme from WRF model (details below)

            3 = 'revised' surface-layer scheme from WRF model (details below)

            4 = GFDL surface layer  (as configured in HWRF_v4.0a)

            5 = Monin-Obukhov Similarity Theory (MOST) for LES 

            6 = MYNN surface layer 

            7 = MYJ surface layer 

      Further information:

        - sfcmodel=1 : Uses simple formulations wherein surface exchange 
        coefficients are specified:  see "Options for sfcmodel = 1" section 
        below.  For diagnostic surface layer calculations (such as 10-m winds),
        a neutrally stratified surface layer is assumed. 
        Notes:  - surface temperature remains fixed over time
                - surface moisture availability remains fixed over time
                - sfcmodel=1 requires oceanmodel=1

        - sfcmodel=2 : Uses the MM5/WRF similarity theory code for the surface 
        layer:  based on Monin-Obukhov with Carslon-Boland viscous sub-layer 
        and standard similarity functions from look-up tables.  
        ("sf_sfclay_physics = 1" in WRF)
        See also "Options for sfcmodel = 2" section below.
        The soil model is the "Thermal diffusion" model from MM5/WRF:  it 
        updates soil temperature only ... soil moisture availability is held 
        fixed over time.  (Same as "sf_surface_physics = 1" in WRF)
        Notes:  - sfcmodel=2 can be used with either oceanmodel=1,2

        - sfcmodel=3 : A revised version of sfcmodel=2.  See Jimenez et al 
        (2012, MWR, pg 898) for more details.
        Notes:  - sfcmodel=3 can be used with either oceanmodel=1,2

        - sfcmodel=4 : GFDL surface layer, as configured in HWRF_v4.0a

        - sfcmodel=5 : Monin-Obukhov Similarity Theory (MOST). 

        - sfcmodel=6 : MYNN surface layer (from WRFV4.2)

        - sfcmodel=7 : MYJ surface layer (from WRFV4.2)


 oceanmodel - Model for ocean/water surface:

            1 = fixed sea-surface temperature
            2 = ocean mixed layer model
                (Same as "omlcall = 1" in WRF)
                (Note:  oceanmodel=2 requires sfcmodel=2,3)
                Ref:  Pollard et al, 1973, Geophys. Fluid Dyn., 3, 381-404.

 --------
 Options for initialization of surface conditions:

 initsfc - initial surface conditions: 
              1 = constant values  (set tsk0,tmn0,xland0,lu0 below)
              2 = sea breeze test case from WRF
              3 = rough surface to west; smoother surface to east
              4 = coastline (land to west, ocean to east)

              for any other value:  you must initialize the surface conditions
              yourself in the "init_surface.F" file.

 tsk0 - default initial value for "skin temperature" (K) of soil/water  
        (~1 cm deep)
        NOTE:  this replaces sea surface temperature (tsurf) in cm1r15

 tmn0 - default initial value for deep-layer temperature (K) of soil
        (Note:  remains fixed throughout simulation)
        (only used if sfcmodel=2)
        (only used over land ... ignored over water)

 xland0 - default initial value for land/water flag: 
             1 for land,   2 for water 

 lu0 - default initial value for land-use index   (see LANDUSE.TBL file)
        (NOTE:  for water/ocean, use lu0 = 16)

 season - which set of land-use conditions to use from LANDUSE.TBL file:
            1 = summer values
            2 = winter values

    c-------------------------------------------------------------c
     To reiterate:  if you want to use spatially varying values of 
     tsk,tmn,xland,lu then you must code it up yourself in the 
     "init_surface.F" file.
    c-------------------------------------------------------------c

 --------
 Options for sfcmodel = 1:

 cecd - When bbc=3 and/or isfcflx=1, this allows the user to choose the 
        formulation for the surface exchange coefficients for enthalphy (Ce)
        and momentum (Cd).  Options are:
                1 = constant value:  user must set cnstce and/or cnstcd below
                    (applies to land and water)
                2 = Deacon's formula  [eg, Rotunno and Emanuel (1987, JAS)]
                    (over water only) 
                    (WRF LANDUSE.TBL used over land)
       default: 3 = Cd based roughly on Fairall et al (2003) at low wind speeds
                                    and Donelan (2004, GRL) at high wind speeds
                    Ce constant, based on Drennan et al. (2007, JAS)
                    (over water only) 
                    (WRF LANDUSE.TBL used over land)

 pertflx - Use only perturbation winds for calculation of surface fluxes?
             (0=NO,  1=yes)
           (Only available with sfcmodel=1)

 cnstce - Constant value of Ce (surface exchange coefficient for enthalpy)
          if isfcflx=1 and cecd=1

 cnstcd - Constant value of Cd (surface exchange coefficient for momentum)
          if bbc=3 and cecd=1

 --------
 Options for sfcmodel = 2,3,6:

 isftcflx - Use alternative Ck and Cd for tropical storm applications:
              (0=off) (For Cd: 1,2 = Donelan)
                      (For Ce: 1=constant Z0q, 2=Garratt)
            MYNN has other options: see module_sf_mynn.F for details
            (Cannot be used with sfcmodel=1,4,5)

 iz0tlnd  - When using sfcmodel=2, option for thermal roughness length:
                 0 = Carlson-Boland  (original mm5/wrf version)
                 1 = Czil_new        (depends on vegetation height)

 --------
 Options for oceanmodel = 2:

 oml_hml0  - default ocean mixed layer depth (m) at initial time

 oml_gamma - default deep water lapse rate (K m-1)

 --------
 Further options for sfcmodel = 1 only:

 set_flx - impose constant surface heat fluxes  (0=no, 1=yes)

 cnst_shflx  - value for surface sensible heat flux (K m/s) if set_flx=1 

 cnst_lhflx - value for surface latent heat flux (g/g m/s) if set_flx=1

 set_znt - impose constant surface roughness length  (0=no, 1=yes)

 cnst_znt - value of surface roughness length (z0, meters) if set_znt=1

 set_ust - impose constant surface friction velocity (0=no, 1=yes)

 cnst_ust - value of surface friction velocity (u-star, m/s) if set_ust=1

 --- NEW  ---
 ramp_sgs - gradually turn on (ie, "ramp up") the subgrid-scale model for 
            LES  (0=no, 1=yes)   (for sgsmodel >=1 only)

 --- NEW  ---
 ramp_time - for ramp_sgs=1 only: this is the time over which the LES 
             subgrid-scale model is linearly ramped up. 

 --- NEW  ---
 t2p_avg - for two-part models (sgsmodel=3,4) use either:
              1 - spatial average (over entire domain at each model level)
              2 - time average (each grid point has a different time-avg)


-------------------------------------------------------------------------
  param4 section -- Horizontally stretched (x) grid options.

          See README.stretch for more information.

 stretch_x - Use horizontally stretched grid in x?
                0 = no
                1 = yes, stretching on both west and east sides of domain
                2 = yes, stretching on west side of domain only
                3 = Arbitrary, specified by user in input_grid_x file
                   (note: options below are ignored for stretch_x = 3)

 dx_inner - Smallest grid spacing (m).

 dx_outer - Largest grid spacing, at edge of domain (m).

 nos_x_len - Length of the no-stretching part of domain (m).

 tot_x_len - Total length of the domain (m).

-------------------------------------------------------------------------
  param5 section -- Horizontally stretched (y) grid options.

          See README.stretch for more information.

 stretch_y - Use horizontally stretched grid in y?
                0 = no
                1 = yes, stretching on both south and north sides of domain
                2 = yes, stretching on south side of domain only
                3 = Arbitrary, specified by user in input_grid_y file
                   (note: options below are ignored for stretch_y = 3)

 dy_inner - Smallest grid spacing (m).

 dy_outer - Largest grid spacing, at edge of domain (m).

 nos_y_len - Length of the no-stretching part of domain (m).

 tot_y_len - Total length of the domain (m).

-------------------------------------------------------------------------
  param6 section -- Vertically stretched grid options.

          See README.stretch for more information.

 stretch_z - Use vertically stretched grid spacing?  
               0 = no
               1 = Wilhelmson and Chen
               2 = Smooth geometric (L. Wicker) Stretch factor, 
                          dz(k+1) = stretch*dz(k)
                   is determined iteratively. Max stretch factor 
                   of 1.1 to limit finite difference errors. If an error
                   is issued, either increase nz or adjust other parameters
                   and try again.
               3 = Arbitrary, specified by user in input_grid_z file.
                   Note: this option specifies the heights of scalar
                   levels (also known as "half levels").
                   For example:  0.5*dz , 1.5*dz , 2.5*dz , etc.
                   Total number of levels specified in file is nz. 
                   See run/config_files/les_StratoCuDrizzle/input_grid_z 
                   for an example.
                   (note: options below are ignored for stretch_z = 3)
               4 = Arbitrary, specified by user in input_grid_z file.
                   Note: this option specifies the heights of w levels 
                   (also known as "full levels").
                   For example:  0.0*dz , 1.0*dz , 2.0*dz , 3.0*dz , etc.
                   Total number of levels specified in file is nz+1. 
                   (note: options below are ignored for stretch_z = 4)

 ztop - Total depth of the domain (i.e., the height of the top of the
        domain) (m).

 str_bot - Level where stretching begins (m). For stretch_z=2, this sets the
           value of nbndlyr = Max(0,Int( str_bot/dz_bot + 0.01) - 1)

 str_top - Level where stretching ends (m).  (diagnosed iteratively for stretch_z = 2)

 dz_bot - Grid spacing at (and below) str_bot (m).

 dz_top - Grid spacing at (and above) str_top (m). For stretch_z=2, this sets
          the maximum dz (which might not be reached, depending on parameters)

-------------------------------------------------------------------------
  param7 section -- Options relating to Direct Numerical Simulation (DNS)
        (used ONLY when turbconf = 3)

 bc_temp - top/bottom boundary condition for potential temperature

             1 = constant theta is specified at boundaries
             2 = constant flux is specified at boundaries

 ptc_top - potential temperature closure for top of model

        if bc_temp = 1, this is theta at the top boundary (K)
        if bc_temp = 2, this is the flux at the top boundary (K m/s)

 ptc_bot - potential temperature closure for bottom of model

        if bc_temp = 1, this is theta at the bottom boundary (K)
        if bc_temp = 2, this is the flux at the bottom boundary (K m/s)

 viscosity - value for kinematic viscosity (m^2/s)

 pr_num - value for Prandtl number (unitless)

-------------------------------------------------------------------------
  param8 section - flex variables

 var1,var2,var3 ... var20 - Use these variables to easily change parameters
 in the model without re-compiling the code.

 Example:  A user wants to run a series of simulations in which the low-level
 shear is changed.  In file "base.F", the user can set "uconst2 = var1", then
 compile the code.  Then, var1 can be changed in the namelist.input file
 and the initial wind profile will change without needing to re-compile.

 Example: A user wants to change the location, size, and amplitude of the 
 initial thermal bubble.  Using the flex vars, the code could be modified in 
 this manner:
        ric     =    var1
        rjc     =    var2
        zc      =    var3
        bhrad   =    var4
        bvrad   =    var5
        bptpert =    var6

 and then the values in namelist.input might be set as follows:
        var1    =  50000.0,
        var2    =  25000.0,
        var3    =   2000.0,
        var4    =  10000.0,
        var5    =   1000.0,
        var6    =      1.5,

 Notes: - the flex variables are all real (float) variables.
        - the flex variables are automatically passed into most subroutines
          in the model


-------------------------------------------------------------------------
  param9 section - output options


 output_format - Specifies the format of the output files:

                      1 = GrADS format
                      2 = netcdf

          !---  IMPORTANT NOTE about large runs  ---!
          ! NOTE:  output_format=2 (netcdf-format output files) tends
                   to be inefficient in CM1 for large processor/core 
                   counts, ie, when using >2,000 (roughly) cores.
                   For many thousands of processors/cores, we recommend 
                   using output_format=1 (also known as unformatted,
                   direct-access binary format) and then converting the 
                   subsequent output files to netcdf format (if desired) 
                   using a conversion program; CDO (climate data operators)
                   is useful to this end.  On NCAR's derecho, use the two
                   following commands:
                            module load cdo
                            cdo -f nc4 import_binary cm1out_s.ctl cm1out_s.nc
                  (see https://code.mpimet.mpg.de/projects/cdo for more info on CDO)


 output_filetype - Type of output file:
                      1 = all output goes into one file
                          (note:  the file size can become very big!)
                      2 = one output file per output time
                          (note:  good to use if you want to keep
                           individual file sizes relatively small ...
                           but this produces many files)
                      3 = for MPI runs only:  one output file per output time
                          ...AND... one output file per MPI process
                          (note:  creates many output files that need to 
                           be combined together using special code ... see, eg,
                           http://www2.mmm.ucar.edu/people/bryan/cm1/programs/)
                          but this is the most efficient way to write output
                          for MPI runs with very large number of processes
                          (say, > 200).

 output_interp - For simulations with terrain, this option (0=no, 1=yes)
                 will generate a second set of output files, with "_i"
                 as part of the filename, wherein the output has been
                 interpolated to the nominal model height levels.  This
                 differs from all other model output files, for which the
                 output is on the native terrain-following coordinate
                 surfaces.

     ------ For the remaining variables, 0=no and 1=yes ------

 output_rain     - surface rainfall.  If imove = 1, two output fields are 
                   generated.  The first (rn) is the accumulated rainfall
                   at model grid points.  The second (rn2) is the translated
                   rainfall pattern, assuming a lower surface is moving at 
                   umove and vmove.

 output_sws      - maximum surface wind speed (aka, surface wind swath, sws)
                   If imove = 1, two output fields are generated.  The first 
                   (sws) is the max sws at model grid points.  The second 
                   (sws2) is the translated sws pattern, assuming a lower 
                   surface is moving at umove and vmove.

   Similar to output_sws:

 output_svs - maximum vertical vorticity at lowest model level

 output_sps - minumum pressure perturbation at lowest model level

 output_srs - maximum rainwater mixing ratio at lowest model level

 output_sgs - maximum graupel/hail mixing ratio at lowest model level

 output_sus - maximum w at 5 km AGL (i.e., maximum updraft swath)

 output_shs - maximum integrated updraft helicity swath

 output_coldpool - properties of surface-based cold pools:

      cpc = cold pool intensity, C
      cph = cold pool depth, h

            NOTE:  the reference profile (for calculation of buoyancy) is 
                   simply the initial sounding.

 output_sfcflx   - surface fluxes of potential temperature and water vapor,
                   and Ce/Cd (if isfcflx=1)

 output_sfcparams - parameters used in the surface/soil/ocean models

                  - includes two-dimensional (surface AND top-of-atmosphere)
                    parameters from radiation scheme (if radopt >= 1)

 output_sfcdiags - diagnostics from surface-layer parameterization
                   (eg, 10-m winds, 2-m temp/moisture, roughness length, etc)

 output_psfc     - surface pressure (z=0 ... not lowest model level)?

 output_zs       - terrain height?

 output_zh       - height on model levels?

 output_basestate - output the base-state arrays?

 output_th       - potential temperature?

 output_thpert   - potential temperature perturbation?

 output_prs      - pressure?

 output_prspert  - pressure perturbation?

 output_pi       - nondimensional pressure (ie, Exner function)?

 output_pipert   - nondimensional pressure perturbation?

 output_rho      - dry air density?

 output_rhopert  - dry air density perturbation?

 output_tke      - subgrid turbulence kinetic energy?

 output_km       - subgrid eddy viscosity?

 output_kh       - subgrid eddy diffusivity?

 output_qv       - water vapor mixing ratio?

 output_qvpert   - perturbation water vapor mixing ratio?

 output_q        - liquid and solid water mixing ratios?
                   (and number concentrations for double-moment schemes)

 output_dbz      - reflectivity (dBZ)?  Only available for ptype=2,3,5,26-28,50-53,55
                   (Goddard-LFO, Thompson, Morrison, NSSL, P3, Jensen ISHMAEL schemes,
                       respectively) also outputs composite reflectivity (cref)
                                (i.e., max reflectivity in the column)

 output_buoyancy - buoyancy (relative to model's base state)?

      ---- all variables above here are in the scalar file (_s) ---

 output_u        - u-velocity?  (_u file)

 output_upert    - u-velocity perturbation?  (_u file)

 output_uinterp  - u-velocity interpolated to scalar points?  (_s file)

 output_v        - v-velocity?  (_v file)

 output_vpert    - v-velocity perturbation?  (_v file)

 output_vinterp  - v-velocity interpolated to scalar points?  (_s file)

 output_w        - w-velocity?  (_w file)

 output_winterp  - w-velocity interpolated to scalar points?  (_s file)

 output_vort     - vorticity?  (all three components, interpolated
                                to scalar points)  (_s file)
                   (not available when using terrain ... at least for now)

 output_pv       - (dry) potential vorticity?

 output_uh       - vertically integrated (2-5 km AGL) updraft helicity?
                   (Kain et al, 2008, WAF, p 931)

 output_pblten   - tendencies from PBL scheme?  (only if ipbl >= 1)
                   (Note:  since cm1r19, pbl tendencies for theta, qv, 
                    u, and v are including with the "budget" output 
                    variables;  see below)

 output_dissten  - dissipation rate?

 output_fallvel  - fall velocities of hydrometeors?
                   (currently only available for ptype=3,5)
                   (ie, Thompson and Morrison microphysics)

 output_nm       - squared Brunt-Vasala frequency?

 output_def      - deformation?

 output_radten   - radiative tendencies?  (only if radopt >= 1)
                 - temperature tendencies from radiation subroutines
                   (potential temperature tendencies due to radiation
                    are included with the "budget" output variables below)
                 - includes surface and top-of-atmosphere fluxes (eg, OLR)

 output_cape -  convective available potential energy (CAPE)?

 output_cin - convective inhibition (CIN)?

 output_lcl - lifted condensation level (LCL)?

 output_lfc - level of free convection (LFC)?

 output_pwat - precipitable water?

 output_lwp - liquid water path and cloud water path?


      Budget variables: includes tendencies due to advection, turbulence, 
                        microphysics, radiation, pbl scheme, buoyancy, 
                        pressure-gradient acceleration, etc, as calculated
                        within CM1.  See output files for more information.

 output_thbudget - Budget terms for potential temperature (theta)?

 output_qvbudget - Budget terms for water vapor mixing ratio (qv)?

 output_ubudget - Budget terms for u velocity?

 output_vbudget - Budget terms for v velocity?

 output_wbudget - Budget terms for w velocity?

 output_pdcomp  -  Pressure decomposition variables?
   Three diagnostic variables are determined: buoyancy pressure perturbation, 
   non-linear dynamic pressure perturbation, and linear dynamic pressure 
   perturbation (based on the base-state wind profiles). 
   Note: this code does not work with distributed-memory parallelization (MPI)
   yet. 

----------------------------------------------------------------------------
            Variables/settings for restart files

 restart_format - Specifies the format of the restart files:

                      1 = binary  (default; highly recommended)
                      2 = netcdf  (note: inefficient for large number of processors)

 restart_filetype - Type of restart file:
                    (since cm1r18.3:  very similar to output_filetype)

                      1 = all data goes into one file
                          (note:  the file size can become very big!)
                          (only available for netcdf format, at the moment)

                      2 = one restart file per restart time
                          (note:  good to use if you want to keep
                           individual file sizes relatively small ...
                           but this produces many files)

                      3 = for MPI runs only:  one restart file per restart time
                          ...AND... one restart file per node (base on ppnode)
                          (faster when using large numbers of procs, ie > 1000)
                          (requires restart_format=1, at the moment)

              Summary:  For restart_format=1 (binary format), 
                        only restart_filetype=2,3 are available

                        For restart_format=2 (netcdf format), 
                        only restart_filetype=1,2 are available


 restart_reset_frqtim  -  Reset next output/stat/restart time based on 
                          parameters in namelist.input file?

                          .true. means the next time CM1 writes 
                          output/stats/restart files will be current time 
                          (ie, the time at restart) plus tapfrq/statfrq/rstfrq.
                          In other words, the user decides when the next file
                          write will be based on namelist parameters.
                          (This option should be used in most instances.)
                          (note:  this is the default)

                          .false. means the next time CM1 writes
                          output/stats/restart files will be whatever time
                          is specified in the restart files, meaning whatever
                          time would have occurred in the original simulation 
                          (i.e., the simulation that wrote the restart file).
                          (This option is useful for so-called "true" 
                           restarts, meaning the user has not changed anything
                           in the namelist, and merely wants to continue a
                           previous simulation.)
                          (This option may be needed to obtain bit-identical 
                           restarts in some instances.)


        --- Additional variables in restart files: ---

         (NOTE: None of the following are required for bit-identical restarts.
                These options have been made available for specialized applications 
                only, eg data assimilation systems.)

 restart_file_theta  -  potential temperature (base state + perturbation)?

 restart_file_dbz  -  reflectivity?

 restart_file_th0  -  base-state potential temperature?

 restart_file_prs0  -  base-state pressure?

 restart_file_pi0  -  base-state nondimensional pressure?

 restart_file_rho0  -  base-state density?

 restart_file_qv0  -  base-state water vapor mixing ratio?

 restart_file_u0  -  base-state x-component velocity (u)?

 restart_file_v0  -  base-state y-component velocity (v)?

 restart_file_zs  -  terrain height?

 restart_file_zh  -  height of half (scalar) grid points (3d array)?

 restart_file_zf  -  height of full (w) grid points (3d array)?

 restart_file_diags  -  theta and qv diags?
                        (dissipative heating if output_thbudget=1 and idiss=1)
                        (theta tndcy from microphys if output_thbudget=1)
                        (hydrometeor fall velocities if output_fallvel=1)


        --- For reading restart files (i.e., when irst=1): ---

 restart_use_theta  -  Use total potential temperature on restarts? 

                       .false. means that perturbation potential temperature
                       is read from the restart file and used (as is) for the
                       simulation.  (Note:  this is needed to obtain 
                       bit-identical restarts with CM1.) 
                       (note:  this is the default)

                       .true. means that the total (i.e. base-state + 
                       perturbation) potential temperature is read from 
                       the restart file, then perturbation potential temp
                       is calculated and then used for the simulation.
                       (Can only be used if restart_file_theta was set to .true.
                        when the restart file was written.)
                       (Note:  bit-identical restarts are usually not 
                        possible for this option, because of round-off errors
                        associated with adding, then subtracting, the base-state
                        field.)
           

                       Default is .false.  

                       Note that the .true. option shouldn't be used in most
                       cases, but it can be useful for some applications, such 
                       as ensemble data assimilation.


-------------------------------------------------------------------------
  param10 section - statistical output options
                      (0=no,  1=yes)

 stat_w        - max/min vertical velocity

 stat_wlevs    - max/min vertical velocity at five different levels
                 (0.5, 1, 2.5, 5, 10 km AGL)

 stat_u        - max/min horizontal (x-direction) velocity
                 (radial velocity if axisymm=1)

 stat_v        - max/min horizontal (y-direction) velocity
                 (azimuthal velocity if axisymm=1)

 stat_rmw      - radius of maximum horizontal velocity
                 (for axisymmetric simulations only)

 stat_pipert   - max/min perturbation nondimensional pressure

 stat_prspert  - max/min perturbation pressure

 stat_thpert   - max/min perturbation potential temperature

 stat_q        - max/min moisture variables

 stat_tke      - max/min subgrid turbulence kinetic energy

 stat_km       - max/min turbulence coefficient for momentum

 stat_kh       - max/min turbulence coefficient for scalars

 stat_div      - max/min divergence

 stat_rh       - max/min relative humidity (wrt liquid)

 stat_rhi      - max/min relative humidity (wrt ice)

 stat_the      - max/min equivalent potential temperature

 stat_cloud    - max/min cloud top/bottom

 stat_sfcprs   - max/min pressure at lowest model level

 stat_wsp      - max/min wind speed at lowest model level
                   (includes 10-m windspeed if bbc=3)

 stat_cfl      - max Courant number

               - also prints KSHMAX and KSVMAX
          (analysis of numerical stability for diffusion terms)
          (NOTE:  KSHMAX should ideally be less than 0.125 for 3d runs)

 stat_vort     - max vertical vorticity at several levels

 stat_tmass    - total dry-air mass

 stat_tmois    - total moisture

 stat_qmass    - total mass of moisture variables

 stat_tenerg   - total energy

 stat_mo       - total momentum

 stat_tmf      - total downward/upward mass flux

 stat_pcn      - precipitation/moisture statistics

 stat_qsrc     - sources of moisture mass


-------------------------------------------------------------------------

  param13 section  -  parcel output options (for iprcl = 1)
                      (0=no,  1=yes)

 prcl_th    -  potential temperature

 prcl_t     -  temperature

 prcl_prs   -  pressure

 prcl_ptra  -  passive tracer concentration (if iptra = 1)

 prcl_q     -  mixing ratios of moisture variables (if imoist = 1)

 prcl_nc    -  number concentrations (for double-moment microphysic schemes)

 prcl_km    -  viscosity

 prcl_kh    -  diffusivity

 prcl_tke   -  subgrid turbulence kinetic energy (if sgsmodel = 1)

 prcl_dbz   -  reflectivity

 prcl_b     -  buoyancy

 prcl_vpg   -  vertical perturbation pressure gradient acceleration

 prcl_vort  -  vertical vorticity

 prcl_rho   -  dry-air density

 prcl_qsat  -  saturation vapor pressure

 prcl_sfc   -  surface variables  (roughness length and friction velocity)


-------------------------------------------------------------------------
  nssl2mom_params section - parameters/options for NSSL microphysics scheme
    (used only if ptype=26,27,28)
    Default values are noted. There is an additional internal namelist that
    can be accessed for more options. See README.NSSLmp
  
  ihlcnh       - Graupel -> hail conversion option (old default was 1, new default is 3)
  rho_qr       - Rain density (1000 kg m**-3)
  rho_qs       - Snow density (100 kg m**-3)
  rho_qh       - Graupel density (500 kg m**-3) 
                 [Only for ptype=28 (single moment) or with nssl_density_on=.false.]
  rho_qhl      - Hail density (800 kg m**-3) 
                 [Only for ptype=28 (single moment) or with nssl_density_on=.false.]
  cnor         - Rain intercept (8.e5 m**-4) 
                 [Only for ptype=28 (single moment)]
  cnos         - Snow intercept (3.e6 m**-4) 
                 [Only for ptype=28 (single moment)]
  cnoh         - Graupel/hail intercept (4.e5 m**-4) 
                 [Only for ptype=28 (single moment)]

  ccn          - Initial concentration of cloud condensation nuclei
                   0.25e+9 maritime
                   0.6e+9 "low" continental (default)
                   1.0e+9 "med-high" continental
                   1.5e+09 - high-extreme continental CCN)
                 Value sets the concentration at MSL, and an initially
                 homogeneous number mixing ratio (ccn/1.225) is assumed throughout the depth of
                 the domain.
  
  irenuc (new option) : 2 = ccn field is UNactivated aerosol (previous default; old droplet activation)
                        7 = ccn field is ACTVIATED aerosol (new default 2023) (new droplet activation)

  infall       - Two-moment sedimation options (default infall=4 recommended)
                          ! 0 -> uses number-wgt for N; NO correction applied 
                                 (results in excessive size sorting)
                          ! 1 -> uses mass-weighted fallspeed for N ALWAYS
                                 (prevents size sorting)
                          ! 2 -> uses number-wgt for N and mass-weighted
                                 correction for N 
                                 (Method II in Mansell, 2010 JAS)
                          ! 3 -> uses number-wgt for N and Z-weighted
                                 correction for N 
                                 (Method I in Mansell, 2010 JAS)
                          ! 4 -> Hybrid of 2 and 3: Uses minimum N from each
                                 method (z-wgt and m-wgt corrections) 
                                 (Method I+II in Mansell, 2010 JAS)
  alphah       - Shape parameter for graupel (0.0)
  alphahl      - Shape parameter for hail (1.0)
  
  Less-used parameters (see code)
  imurain      - DSD function option for rain (1 - gamma of diameter)
  icdx         - fall speed option for graupel (default was 3, now is 6)
  icdxhl       - fall speed option for hail (default was 3, now is 6)

  dfrz         - Minimum diameter of new frozen drops (graupel);
                 default 0.15e-3, but can be set larger (up to 0.5e-5) to produce
                 larger graupel
  hldnmn       - Minimum hail particle density (500 kg m**-3) (changed from 750 in r17)
  iehw,iehlw   - Graupel (iehw) and Hail (iehlw) droplet collection efficiency
                 options 
  ehw0,ehlw0   - Maximim droplet collection efficiencies for graupel (default was ehw0=0.75, now 0.9)
                 and hail (default was ehlw0=0.75, now 0.9)
  dmrauto      - Options for limiting rain autoconversion

-------------------------------------------------------------------------

  param14 section:  Options related to domain-wide diagnostics.


  dodomaindiag - Calculate and write domain-wide / domain-averaged diagnostics 
               (eg, domain-average vertical profiles, variances, vertical
                fluxes, etc)
              (logical:  .true. or .false.)

                     Can now be written in CF-compliant netcdf format! 
                     Use output_format = 2

  diagfrq - Frequency (seconds) for calculating/writing diagnostics.


-------------------------------------------------------------------------

  param15 section:   Options related to azimuthally averaged output.


  doazimavg  -  Calculate and write azimuthally averaged vertical cross 
                sections?
                (logical:  .true. or .false.)

                      Can now be written in CF-compliant netcdf format! 
                      Use output_format = 2

      Note:  To define the center-point for these calculatins, there is 
             only one option, at the moment:  the code searches for the 
             point that yields maximum possible azimuthally averaged  
             tangential velocity, and defines this point as "storm center."
             Other mehtods could be coded in the future

  azimavgfrq  -  Time between writes to azimavg files (s).

  rlen  -  Length (m) of analysis grid (in radial direction).

  do_adapt_move  -  Automatically identify center of a storm (using min 
           surface pressure), and adaptively adjust umove,vmove to keep 
           the storm near the center of the domain.

           Caution: this is an experimental capability in CM1 
           (although it seems to work well in our tests).

  adapt_move_frq  -  Frequency (seconds) to identify center of storm, 
           and update umove,vmove (when do_adapt_move = .true.)



-------------------------------------------------------------------------
                        --- NEW ---

  param17 section:  options for cm1setup=4 only 
                    (ie, for LES-within-Mesoscale-Model framework)


 les_subdomain_shape  -  1 = square-shaped LES subdomain 
                         2 = circular-shaped LES subdomain (not implemented yet)

 les_subdomain_xlen  -  size of LES subdomain in x direction 
                        (for les_subdomain_shape=1 only)
                    (Note: for stretch_x=1, this should be similar to nos_x_len)

 les_subdomain_ylen  -  size of LES subdomain in y direction
                        (for les_subdomain_shape=1 only)
                    (Note: for stretch_y=1, this should be similar to nos_y_len)

 les_subdomain_dlen  -  diameter of LES subdomain for les_subdomain_shape=2 only
                           (not implemented yet)

 les_subdomain_trnslen  -  a transition scale (meters) for expected development of 
                           resolved-scale turbulence.  Specifically, this is how 
                           far into the LES subdomain at which the LES subgrid 
                           model is completely turned on. 


-------------------------------------------------------------------------
                        --- NEW ---

  param18 section:  options for eddy recycling
                    (ie, turbulent fluctuations are nudged into the flow 
                     at either open boundaries when cm1setup=0,1,3, or at 
                     the edges of an LES subdomain when cm1setup=4)

 do_recycle_w  -  use eddy recycling near western boundary of an LES domain?
                  (T/F)

 do_recycle_s  -  use eddy recycling near southern boundary of an LES domain?
                  (T/F)

 do_recycle_e  -  use eddy recycling near eastern boundary of an LES domain?
                  (T/F)

 do_recycle_n  -  use eddy recycling near northern boundary of an LES domain?
                  (T/F)

 recycle_width_dx  -  width of the recycled data, in number of grid points 
                      (eg, a value of 6.0 means 6*dx)

 recycle_depth_m  -  depth of recycled data, from the surface (meters)

 recycle_cap_loc_m  -  distance (meters) from the edge of the LES domain 
                       where turbulent fluctuations are "captured"

 recycle_inj_loc_m  -  distance (meters) from the edge of the LES domain 
                       where turbulent fluctuations are "injected"


-------------------------------------------------------------------------
                        --- NEW ---

  param19 section:   Options related to large-scale (i.e., domain-average) nudging.

 (for details, see Appendix of Alland et al. 2021, doi:10.1175/JAS-D-20-0054.1)


  ! use large-scale wind nudging?

  do_lsnudge  -  use large-scale nudging?   (T/F)

   ! Note: see the documentation at the top of lsnudge.F for some important info


        ! for do_lsnudge=.true., select which variables are nudged:

  do_lsnudge_u    -   nudge domain-avg u wind profile?
  do_lsnudge_v    -   nudge domain-avg v wind profile?
  do_lsnudge_th   -   nudge domain-avg potential temperature (theta) profile?
  do_lsnudge_qv   -   nudge domain-avg water vapor mixing ratio (qv) profile?



  lsnudge_tau   -   time scale (seconds) for damping term


  ! NOTE:  CM1 only applies large-scale nudging when t > lsnudge_start and 
  !        t < lsnudge_end;  the user must set these two variables below.  
  !        For any other times, large-scale nudging is not applied, regardless 
  !        of the times provided for lsnudge_time1 and lsnudge_time2 in the 
  !        lsnudge_xxxx.dat files.

  lsnudge_start  -   time (seconds) to begin large-scale nudging

  lsnudge_end    -   time (seconds) to end large-scale nudging

  lsnudge_ramp_time  -   time (seconds) over which to gradually introduce 
     nudging; that is, for  
         lsnudge_start  <  t (seconds)  <  lsnudge_start+nudge_ramp_time 

     Nudging term is linearly increased from 0 to full nudging term
     over this time period.

     (set lsnudge_ramp_time to zero for instantaneous introduction of nudging)


-------------------------------------------------------------------------
                        --- NEW ---

  param20 section:  options related to immersed boundaries 


 do_ib  -  use immersed boundaries?  (T/F)

 ib_init  -  distribution of immersed boundaries 
             (see ib_module.F for more details)

 top_cd  -  drag coefficient (unitless) on top of immersed boundaries

 side_cd  -  drag coefficient (unitless) on sides of immersed boundaries 


-------------------------------------------------------------------------
                        --- NEW ---

  param21 section:  options related to idealized hurricane test cases

          see Bryan et al 2017 (doi:10.1007/s10546-016-0207-0)
              Chen et al 2021 (doi:10.1175/JAS-D-20-0227.1)


 hurr_vg  -  Gradient wind (m/s)

             This is also the default initial wind speed in idealized 
             initial conditions.
             (assumes iwnd=8)

 hurr_rad  -  Radius (meters) from center of tropical cyclone

 hurr_vgpl  -  The "power law" (or "radial decay" parameter) for the 
               gradient wind  (see Eqn 20 in Bryan et al. 2017). 
               (Note: should be a negative number)

 hurr_rotate  -  How much to rotate winds and pressure gradient from the 
                 convention used by Bryan et al (2017). 

                 Here, positive is clockwise, negative is counterclockwise.  

                 Example 1: for hurr_rotate = 0.0, winds are southerly 
                 (ie, from south-to-north) and the TC center is to the west 
                 (see Fig 1 of B17). 

                 Example 2: for hurr_rotate = -90.0, winds are easterly 
                 (ie, from east-to-west) and the TC center is to the south. 

                 Example 3: for hurr_rotate = 90.0, winds are westerly 
                 (ie, from west-to-east) and the TC center is to the north.

                 (Note: Nothern hemisphere is assumed.) 


-------------------------------------------------------------------------


Questions, comments, suggestions:  send email to gbryan@ucar.edu

