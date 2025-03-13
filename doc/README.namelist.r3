
User's Guide to namelist.input file (located in "run" directory)

Bryan Cloud Model, Version 1, Release 3
16 March 2003

-------------------------------------------------------------------------
  param1 section -- enter REAL values

 dx - Horizontal grid spacing (m).

 dz - Vertical grid spacing (m).

 dtl - Large time step (s).

       For kwss=1, this time step is limited by the fastest nonacoustic 
       speed.  For thunderstorm simulations, this is usually the maximum 
       vertical velocity.  Otherwise, this would be the propagation speed
       of gravity waves.
       The following usually works well:
           dtl = min(dx,dz)/67  (rounded to an appropriate value, of course)

       For kwss=0, this time step is limited by the propagation speed of
       sound waves.  dtl of about min(dx,dz)/700 is recommended.

 timax - Maximum forecast time (s).

 tapfrq - Frequency of three-dimensional model output (s).
          Output is in cm1out files.

 rstfrq - Frequency to save model restart files (s).  Set to a negative
          number if restart files are not desired.

-------------------------------------------------------------------------
  param2 section -- enter INTEGER values
     (value in CAPS is recommended, where applicable)

 statfrq - Frequency for calculating some interesting output.
           Output is in stats.dat file.
           Set as number of time steps.  (e.g., 1 = every model
           time step, 5 = every fifth time step)

 irst - Is this a restart?  (0=no, 1=yes)

 iconly - Setup initial conditions only?
              1 = creates initial conditions, but does not run model.
              0 = creates initial conditions, and proceeds with integration.

 advorder - Order of advection scheme. 5=fifth order, 6=sixth order

           Fifth order scheme has implicit diffusion.  If fifth order
           scheme is used, idiff should be set to 0 (i.e., no additional
           artifical diffusion is typically necessary).

           Sixth order scheme requires additional artifical diffusion
           for stability (i.e., idiff=1 is recommended).  User should
           use difforder=6 and value of kdiff6 between about 0.02-0.06

 idiff - Include additional artificial diffusion?  (0=no, 1=yes)

           For idiff=1, user must set difforder and kdiff2 or kdiff6.

 difforder - Order of diffusion scheme.  2=second order, 6=sixth order.

           Second order diffusion is not generally recommended.  It is
           only used for certain idealized cases.  kdiff2 must be set
           appropriately when difforder=2.

           Sixth order diffusion is recommended for general use when
           diffusion of small scales (2-6 delta) is needed.  User must
           also set kdiff6 when difforder=6.

 imoist - Include moisture?  (0=no, 1=yes)

 iturb - Include subgrid turbulence scheme?  (0=no, 1=TKE, 2=Smagorinsky)

 smix - Use new conserved-variable saturated mixing scheme?  (not available)

 irdamp - Use upper-level Rayleigh damping zone?  (0=no,  1=YES)
          (User must set rdalpha and zd below)

 implsound - Use vertical implicit scheme for acoustic modes?
                 (Only used if kwss = 1)
              0 = explicit treament of acoustic terms in both vertical 
                  and horizontal directions
              1 = vertical implicit sound (only horizontal terms are explicit,
                  as in MM5, ARPS, WRF)

 kwss - Use Klemp-Wilhelmson split time steps?  (0=no,  1=YES)
         (User must set nsound if kwss=1)

 nsound - Number of small (acoustic) time steps per large time step.

       Only used for kwss=1.  (Ignored if kwss=0.)

       NOTE!  nsound must be an EVEN integer (exactly) for the
       split time steps to work properly.  There is no check for
       this ... so please remember.  Typically, nsound of
       4,6,8,10,12 is appropriate.  For nsound less than 4,
       you should probably use kwss=0 instead.  A value of
       nsound greater than 12 is not recommended.

       Also keep in mind the small time step (dts=dtl/nsound) that
       is required for stability.

         For implsound=0, dts should be about min(dx,dz)/700.
         For implsound=1, dts should be about dx/700.

 ptype - Explicit moisture scheme.
             1 = Kessler scheme (water only)
             2 = Goddard version of Lin et al. scheme (includes ice)

 ihail - Use hail or graupel for large ice category when ptype=2.

             1 = hail
             0 = graupel

 idrag - Simple surface drag parameterization (0=no, 1=yes)

 isfcflx - Simple surface fluxes of heat and moisture (0=no, 1=yes)

 icor - Include Coriolis force?  (0=no,  1=yes)
        (If user chooses 1, then fcor must be set below)
        f-plane is assumed.

 neweqts - Use new equations (from Bryan and Fritsch 2002) ?
               (1 = YES,  0 = no)

 wbc - West lateral boundary condition.

 ebc - East lateral boundary condition.

 sbc - South lateral boundary condition.

 nbc - North lateral boundary condition.

                        where:  1 = periodic
                                2 = open-radiative
                                3 = rigid walls

 irbc - For bc=2, this is the type of radiative scheme to use:
                    1 = Klemp-Wilhelmson on large steps
                    2 = Klemp-Wilhelmson on small steps
                    (no other options available at this time)

 isnd - Base-state sounding:  1 = Dry adiabatic
                              2 = Dry isothermal
                              3 = Dry constant lapse rate
                              4 = Saturated neutrally stable
                              5 = Weisman-Klemp analytic sounding
                              6 = Sounding for tropical cyclone
                              7 = External file

 iwnd - Base-state wind profile:  0 = zero winds
                                  1 = RKW-type profile
                                  2 = Weisman-Klemp supercell
                                  3 = multicell
                                  4 = Weisman-Klemp multicell

 iinit - 3D initialization option:  1 = warm bubble
                                    2 = cold pool
                                    3 = line of warm bubbles
                                    4 = initialization for moist benchmark
                                    5 = cold blob
                                    6 = idealized tropical cyclone

 ibalance - Specified balance assumption for initial 3D pressure field
            (not used for idealized tropical cyclone)

          0 = no balance (initial pressure perturbation is zero everywhere)
          1 = hydrostatic balance (appropriate for small aspect ratios)
          2 = anelastic balance (initial pressure perturbation is the 
              buoyancy pressure perturbation field for an anelastic 
              atmosphere).  CAUTION:  relatively untested.  Also slow.

 imove - Move domain at constant speed (0=no, 1=yes)

           For imove=1, user must set umove and vmove.

 iptra - integrate passive fluid tracer (0=no, 1=yes)

           User must initialize "pta" array in init.F.
           Currently, only one tracer is supported. 

-------------------------------------------------------------------------
  param3 section -- enter REAL values

 kdiff2 - Diffusion coefficient for difforder=2.  Specified in m^2/s.

 kdiff6 - Diffusion coefficient for difforder=6.  Specified as a 
          fraction of one-dimensional stability.  A value between
          0.02-0.06 is recommended.

 fcor - Coriolis parameter (1/s).

 kdiv - Coefficient for divergence damper.  Value between 0.02-0.06 is
        recommended.

        This is only used when kwss=1.  The divergence damper is
        an artificial term designed to attenuate acoustic waves.

 rdalpha - Inverse folding time for upper-level Rayleigh damping layer
           (1/s).  Value of about 1/300 is recommended.

 zd - Height above which Rayleigh damping is applied (m).

 umove - Constant speed for domain translation in x-direction (m/s)

 vmove - Constant speed for domain translation in y-direction (m/s)


