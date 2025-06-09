!WRF:MODEL_LAYER:PHYSICS
!
! Current version is tagged as V3.6
! Reference is Morrison et al. (2005), JAS, Morrison et al. (2009), MWR

! changes with respect to V1.4

! V1.5
! 1) more pathways to allow hail to form (only affects IHAIL=1 option), from collisions of snow/cloud water
! 2) bug fix to PGAM calculation (multiplication instead of division by air density)

! V1.6
! 1) added parameter TMELT for all calculations involving melting point
! 2) replaced hard-wired gas constant for air with parameter value 'R'

! V1.7
! 1) modification to minimum mixing ratio in dry conditions, change from 10^-6 to 10^-8 kg/kg
!   to improve reflectivity at low mixing ratio amounts
! 2) bug fix to prevent possible division by zero error involving LAMI
! 3) change for liquid saturation vapor pressure, replace old formula with Flatau et al. 1992

! V2
! 1) bug fix to maximum-allowed particle fallspeeds (air density correction factor considered)
! 2) change to comments

! V2.1
! 1) addition of rain drop breakup following Verlinde and Cotton (1993)
! 2) change to minimum allowed lambda (slope parameter) for rain
! 3) addition of accelerated melting of graupel/hail/snow due to collision with rain

! V3
! minor revisions by Andy Ackerman
! 1) replaced kinematic with dynamic viscosity 
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! Additional fixes
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NO SUBTRACTED SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SIN
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS 
! 7) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 8) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 9) BUG FIX TO IGRAUP SWITCH FOR NO GRAUPEL/HAIL

! hm bug fix 3/16/12

! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted

! hm, changes 3/4/13 for V3.3

! 1) removed second initialization of evpms (non-answer-changing)
! 2) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass 
!    collected by rain
! 3) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!    cm-3. This was done to address the problem of excessive and persistent
!    anvil cirrus produced by the scheme, and was found to greatly improve forecasts over
!    at convection-permitting scales over the central U.S. in summertime.

! hm, changes 7/25/13 for V3.4

! 1) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!    Hailong Wang (PNNL)
! 2) very minor fix to immersion freezing rate formulation (negligible impact)
! 3) clarifications to code comments
! 4) minor change to shedding of rain, remove limit so that the number of
!    collected drops can smaller than number of shed drops
! 5) change of specific heat of liquid water from 4218 to 4187 J/kg/K

! hm, changes 1/20/15 for version 3.5

! 1) minor bug fix to diagnostic steady-state supersaturation equation for droplet activation in
!    cloud interior for option IBASE=1, minus sign is missing from the pressure term in this equation,
!    note this only has ~10% effect relative to the DQSDT term so the impact of the fix is small.
!    Bug reported by Xiaowen Li (NASA Goddard).
! 2) minor bug fix to melting of snow and graupel, an extra factor of air density (RHO) was removed
!    from the calculation of PSMLT and PGMLT
! 3) redundant initialization of PSMLT (non answer-changing)

! hm, changes 7/12/16 for version 3.6

! 1) changes and cleanup of code comments and variable names (non-answer changing)
! 2) correction to universal gas constant (very small answer change)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_mp_GRAUPEL
!!!   USE     module_wrf_error
!!!      USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
!!!      USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT

!  USE module_state_description

   implicit none

   private
   public :: graupel_init,mp_graupel

   real, parameter :: PI = 3.1415926535897932384626434
   real, parameter :: xxx = 0.9189385332046727417803297

!!!   PUBLIC  ::  MP_GRAUPEL
!!!   PUBLIC  ::  POLYSVP

   private :: GAMMA, DERF1
   private :: PI, xxx
   private :: M2005MICRO_GRAUPEL

! SWITCHES FOR MICROPHYSICS SCHEME
! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

     integer, private ::  IACT

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

     integer, private ::  INUM

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
     real, private ::      NDCNST

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

     integer, private ::  ILIQ

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS

     integer, private ::  INUC

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

     integer, private ::  IBASE

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

     integer, private ::  ISUB      

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

     integer, private ::  IGRAUP

! HM ADDED NEW OPTION FOR HAIL V1.3
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

     integer, private ::  IHAIL

! HM ADDED 8/1/08, v1.4
! SWITCH FOR WARM RAIN SCHEME
! IRAIN = 0, WARM RAIN (AUTO, ACC, SELF-COLL) FROM KHAIROUTIDNOV AND KOGAN (2000)
! IRAIN = 1, WARM RAIN (AUTO, ACC, SELF-COLL) FROM SEIFERT AND BEHENG (2001)

     integer, private ::  IRAIN      

! CLOUD MICROPHYSICS CONSTANTS

     real, private ::      AI,AC,AS,AR,AG ! 'A' parameter in FALLSPEED-DIAM RELATIONSHIP
     real, private ::      BI,BC,BS,BR,BG ! 'B' parameter in FALLSPEED-DIAM RELATIONSHIP
     real, private ::      R           ! GAS CONSTANT FOR AIR
     real, private ::      RV          ! GAS CONSTANT FOR WATER VAPOR
     real, private ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
     real, private ::      EP_2
     ! GHB: Needed for eqtset=2 in CM1
     real, private ::      cv
     real, private ::      cvv
     real, private ::      cpv
     real, private ::      cpl
     real, private ::      cpi
     real, private ::      xlv0
     real, private ::      xls0
     real, private ::      lv1
     real, private ::      lv2
     real, private ::      ls1
     real, private ::      ls2
     ! GHB
     real, private ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
     real, private ::      RHOW        ! DENSITY OF LIQUID WATER
     real, private ::      RHOI        ! BULK DENSITY OF CLOUD ICE
     real, private ::      RHOSN       ! BULK DENSITY OF SNOW
     real, private ::      RHOG        ! BULK DENSITY OF GRAUPEL
     real, private ::      AIMM        ! parameter in BIGG IMMERSION FREEZING
     real, private ::      BIMM        ! parameter in BIGG IMMERSION FREEZING
     real, private ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
     real, private ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
     real, private ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
     real, private ::      MG0         ! MASS OF EMBRYO GRAUPEL
     real, private ::      F1S         ! VENTILATION parameter FOR SNOW
     real, private ::      F2S         ! VENTILATION parameter FOR SNOW
     real, private ::      F1R         ! VENTILATION parameter FOR RAIN
     real, private ::      F2R         ! VENTILATION parameter FOR RAIN
     real, private ::      G           ! GRAVITATIONAL ACCELERATION
     real, private ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
     real, private ::      QTHR        ! ghb: SMALL HYDROMETEOR MIXING RATIO THRESHOLD in SUBSATURATED CONDITIONS
     real, private ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
     real, private ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
     real, private ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
     real, private ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
! V1.6
     real, private ::      TMELT     ! melting temp (K)
! hm, add for V2.1
     real, private ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

! CCN SPECTRA FOR IACT = 1

     real, private ::      C1     ! 'C' in NCCN = CS^K (CM-3)
     real, private ::      K1     ! 'K' in NCCN = CS^K

! AEROSOL PARAMETERS FOR IACT = 2

     real, private ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
     real, private ::      OSM     ! OSMOTIC COEFFICIENT
     real, private ::      VI      ! NUMBER OF ION DISSOCIATED in SOLUTION
     real, private ::      EPSM    ! AEROSOL SOLUBLE FRACTION
     real, private ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
     real, private ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
     real, private ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
     real, private ::      RR      ! UNIVERSAL GAS CONSTANT
     real, private ::      BACT    ! ACTIVATION parameter
     real, private ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
     real, private ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
     real, private ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
     real, private ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
     real, private ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
     real, private ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
     real, private ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     real, private ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     real, private ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     real, private ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2     
     real, private ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
     real, private ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING
     real, private ::      D_b     ! cutoff diameter for breakup

! CONSTANTS TO IMPROVE EFFICIENCY

     real, private :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
     real, private :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
     real, private :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
     real, private :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
     real, private :: CONS41

! v1.4
     real, private :: dnu(16)

!..Various radar related variables, from GT

!..Lookup table dimensions
      integer, parameter, private:: nbins = 100
      integer, parameter, private:: nbr = nbins
      integer, parameter, private:: nbs = nbins
      integer, parameter, private:: nbg = nbins
      double precision, dimension(nbins+1):: ddx
      double precision, dimension(nbr):: Dr, dtr
      double precision, dimension(nbs):: Dds, dts
      double precision, dimension(nbg):: Ddg, dtg

      double precision, parameter, private:: lamda_radar = 0.10         ! in meters
      double precision, private:: K_w, PI5, lamda4
      complex*16, private:: m_w_0, m_i_0
      double precision, dimension(nbins+1), private:: simpson
      double precision, dimension(3), parameter, private:: basis =      &
                           (/1.d0/3.d0, 4.d0/3.d0, 1.d0/3.d0/)

      integer, parameter, private:: slen = 20
      character(len=slen), private::                                    &
              mixingrulestring_s, matrixstring_s, inclusionstring_s,    &
              hoststring_s, hostmatrixstring_s, hostinclusionstring_s,  &
              mixingrulestring_g, matrixstring_g, inclusionstring_g,    &
              hoststring_g, hostmatrixstring_g, hostinclusionstring_g

      real, parameter, private:: D0r = 50.E-6
      real, parameter, private:: D0s = 100.E-6
      real, parameter, private:: D0g = 100.E-6
      character*256:: mp_debug

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GRAUPEL_INIT(cm1hail,cm1inum,cm1ndcnst)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: cm1hail
      integer, intent(inout) :: cm1inum
      real,    intent(in) :: cm1ndcnst
      integer n,i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
! SET PRIOR TO CODE COMPILATION

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

      INUM = 1

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)

      NDCNST = cm1ndcnst

! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      IACT = 2

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      IBASE = 1

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      ISUB = 1      

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      ILIQ = 0

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

      INUC = 0

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

      IGRAUP = 0

! HM ADDED 11/7/07, V1.3
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL

      IHAIL = cm1hail
 
! HM ADDED 8/1/08, v1.4
! SWITCH FOR WARM RAIN SCHEME
! IRAIN = 0, WARM RAIN (AUTO, ACC, SELF-COLL) FROM KHAIROUTIDNOV AND KOGAN (2000)
! IRAIN = 1, WARM RAIN (AUTO, ACC, SELF-COLL) FROM SEIFERT AND BEHENG (2001)

      IRAIN = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET PHYSICAL CONSTANTS

! FALLSPEED PARAMETERS (V=AD^B)
         AI = 700.
         AC = 3.E7
         AS = 11.72
         AR = 841.99667
         BI = 1.
         BC = 2.
         BS = 0.41
         BR = 0.8
! V1.3
         if (IHAIL.EQ.0) then
	 AG = 19.3
	 BG = 0.37
         else ! (MATSUN AND HUGGINS 1980)
         AG = 114.5 
         BG = 0.5
         end if

! CONSTANTS AND PARAMETERS
         R = 287.04
         RV = 461.5
         CP = 1005.7
         EP_2 = R/RV
         ! GHB: Needed for eqtset=2 in CM1
         cv = CP - R
         cvv = 1408.5
         cpv = 1870.0
         cpl = 4190.0
         cpi = 2106.0
         xlv0 = 2501000.0
         xls0 = 2834000.0
         lv1 = xlv0+(cpl-cpv)*273.15
         lv2 = cpl-cpv
         ls1 = xls0+(cpi-cpv)*273.15
         ls2 = cpi-cpv
         ! GHB
! V1.6
         TMELT = 273.15
! V1.6
         RHOSU = 85000./(R*TMELT)
         RHOW = 997.
         RHOI = 500.
         RHOSN = 100.
! V1.3
         if (IHAIL.EQ.0) then
	 RHOG = 400.
         else
         RHOG = 900.
         end if
         AIMM = 0.66
         BIMM = 100.
         ECR = 1.
         DCS = 125.E-6
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3
	 MG0 = 1.6E-10
         F1S = 0.86
         F2S = 0.28
         F1R = 0.78
! V3 5/27/11
!         F2R = 0.32
! AA revision 4/1/11
         F2R = 0.308
         G = 9.81
         QSMALL = 1.E-14
         QTHR = 1.E-10
         EII = 0.1
         ECI = 0.7
! HM, ADD FOR V3.2
! hm, 7/23/13
!         CPW = 4218.
!         CPW = 4187.
         CPW = cpl

! SIZE DISTRIBUTION PARAMETERS

         CI = RHOI*PI/6.
         DI = 3.
         CS = RHOSN*PI/6.
         DS = 3.
         CG = RHOG*PI/6.
         DG = 3.

! RADIUS OF CONTACT NUCLEI
         RIN = 0.1E-6

         MMULT = 4./3.*PI*RHOI*(5.E-6)**3

! SIZE LIMITS FOR LAMBDA

         LAMMAXI = 1./1.E-6
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./20.E-6
!         LAMMINR = 1./500.E-6
! V2.1
         LAMMINR = 1./2800.E-6

         D_b = 300.e-6  ! cutoff diameter for breakup

         LAMMAXS = 1./10.E-6
         LAMMINS = 1./2000.E-6
         LAMMAXG = 1./20.E-6
         LAMMING = 1./2000.E-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! note: these parameters only used by the non-wrf-chem version of the
!       scheme with predicted droplet number
! CCN SPECTRA FOR IACT = 1

! MARITIME
! MODIFIED FROM RASMUSSEN ET AL. 2002
! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

              K1 = 0.4
              C1 = 120. 

! CONTINENTAL

!              K1 = 0.5
!              C1 = 1000. 

! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

         MW = 0.018
         OSM = 1.
         VI = 3.
         EPSM = 0.7
         RHOA = 1777.
         MAP = 0.132
         MA = 0.0284
         RR = 8.3187
! hm fix 6/23/16
!         RR = 8.3187
         RR = 8.3145
         BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE
! MODE 1

         RM1 = 0.052E-6
         SIG1 = 2.04
         NANEW1 = 72.2E6
         F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
         F21 = 1.+0.25*LOG(SIG1)

! MODE 2

         RM2 = 1.3E-6
         SIG2 = 2.5
         NANEW2 = 1.8E6
         F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
         F22 = 1.+0.25*LOG(SIG2)

! CONSTANTS FOR EFFICIENCY

         CONS1=GAMMA(1.+DS)*CS
         CONS2=GAMMA(1.+DG)*CG
         CONS3=GAMMA(4.+BS)/6.
         CONS4=GAMMA(4.+BR)/6.
         CONS5=GAMMA(1.+BS)
         CONS6=GAMMA(1.+BR)
         CONS7=GAMMA(4.+BG)/6.
         CONS8=GAMMA(1.+BG)
         CONS9=GAMMA(5./2.+BR/2.)
         CONS10=GAMMA(5./2.+BS/2.)
         CONS11=GAMMA(5./2.+BG/2.)
         CONS12=GAMMA(1.+DI)*CI
         CONS13=GAMMA(BS+3.)*PI/4.*ECI
         CONS14=GAMMA(BG+3.)*PI/4.*ECI
         CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
         CONS16=GAMMA(BI+3.)*PI/4.*ECI
         CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
         CONS18=RHOSN*RHOSN
         CONS19=RHOW*RHOW
         CONS20=20.*PI*PI*RHOW*BIMM
         CONS21=4./(DCS*RHOI)
         CONS22=PI*RHOI*DCS**3/6.
         CONS23=PI/4.*EII*GAMMA(BS+3.)
         CONS24=PI/4.*ECR*GAMMA(BR+3.)
         CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
         CONS26=PI/6.*RHOW
         CONS27=GAMMA(1.+BI)
         CONS28=GAMMA(4.+BI)/6.
         CONS29=4./3.*PI*RHOW*(25.E-6)**3
         CONS30=4./3.*PI*RHOW
         CONS31=PI*PI*ECR*RHOSN
         CONS32=PI/2.*ECR
         CONS33=PI*PI*ECR*RHOG
         CONS34=5./2.+BR/2.
         CONS35=5./2.+BS/2.
         CONS36=5./2.+BG/2.
         CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
         CONS38=PI*PI/3.*RHOW
         CONS39=PI*PI/36.*RHOW*BIMM
         CONS40=PI/6.*BIMM
         CONS41=PI*PI*ECR*RHOW

! v1.4
         dnu(1) = -0.557
         dnu(2) = -0.557
         dnu(3) = -0.430
         dnu(4) = -0.307
         dnu(5) = -0.186
         dnu(6) = -0.067
         dnu(7) = 0.050
         dnu(8) = 0.167
         dnu(9) = 0.282
         dnu(10) = 0.397
         dnu(11) = 0.512
         dnu(12) = 0.626
         dnu(13) = 0.739
         dnu(14) = 0.853
         dnu(15) = 0.966
         dnu(16) = 0.966

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables for radar reflecitivity calculations
!..Create bins of rain (from min diameter up to 5 mm).
      ddx(1) = D0r*1.0d0
      ddx(nbr+1) = 0.005d0
      do n = 2, nbr
         ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbr) &
                  *DLOG(ddx(nbr+1)/ddx(1)) +DLOG(ddx(1)))
      enddo
      do n = 1, nbr
         Dr(n) = DSQRT(ddx(n)*ddx(n+1))
         dtr(n) = ddx(n+1) - ddx(n)
      enddo

!..Create bins of snow (from min diameter up to 2 cm).
      Ddx(1) = D0s*1.0d0
      Ddx(nbs+1) = 0.02d0
      do n = 2, nbs
         Ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbs) &
                  *DLOG(Ddx(nbs+1)/Ddx(1)) +DLOG(Ddx(1)))
      enddo
      do n = 1, nbs
         Dds(n) = DSQRT(Ddx(n)*Ddx(n+1))
         dts(n) = Ddx(n+1) - Ddx(n)
      enddo

!..Create bins of graupel (from min diameter up to 5 cm).
      Ddx(1) = D0g*1.0d0
      Ddx(nbg+1) = 0.05d0
      do n = 2, nbg
         Ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbg) &
                  *DLOG(Ddx(nbg+1)/Ddx(1)) +DLOG(Ddx(1)))
      enddo   
      do n = 1, nbg
         Ddg(n) = DSQRT(Ddx(n)*Ddx(n+1))
         dtg(n) = Ddx(n+1) - Ddx(n)
      enddo

      do i = 1, 256
         mp_debug(i:i) = char(0)
      enddo

      call radar_init

      cm1inum   = inum

end subroutine GRAUPEL_INIT

!interface copied from new thompson interface
!and added NC, NS, NR, and NG variables.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE M2005MICRO_GRAUPEL) 
! WHICH OPERATES ON 1D VERTICAL COLUMNS.
! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE

! ******IMPORTANT******
! THIS CODE ASSUMES THE DRIVER MODEL USES PROCESS-SPLITTING FOR SOLVING THE TIME-DEPENDENT EQS.
! THUS, MODEL VARIABLES ARE UPDATED WITH MICROPHYSICS TENDENCIES INSIDE OF THE MICROPHYSICS
! SCHEME. THESE UPDATED VARIABLES ARE PASSED BACK TO DRIVER MODEL. THIS IS WHY THERE
! ARE NO TENDENCIES PASSED BACK AND FORTH BETWEEN DRIVER AND THE INTERFACE SUBROUTINE

! AN EXCEPTION IS THE TURBULENT MIXING TENDENCIES FOR DROPLET AND CLOUD ICE NUMBER CONCENTRATIONS
! (NCTEND, NITEND BELOW). FOR APPLICATION IN MODELS OTHER THAN WRF, TURBULENT MIXING TENDENCIES
! CAN BE ADDED TO THE VARIABLES ELSEWHERE (IN DRIVER OR PBL ROUTINE), AND THEN DON'T
! NEED TO BE PASSED INTO THE SUBROUTINE HERE.....

! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MP_GRAUPEL(ITIMESTEP,                       &
                T3D, pii ,QV, QC, QR, QI, QS, QG, NNI, NS, NR, NG, NC,   &
                     P, RHO, DT_IN, DDZ, W, RAINNC, prate ,         &
                effc,effi,effs,effr,effg,effis,                    &
                tcond,tevac,tevar,train,                           &
                ruh,rvh,rmh,rr,getdbz,                             &
                vtc,vtr,vti,vts,vtg,getvt,dorad,                   &
                dotbud,doqbud,tdiag,qdiag,out3d)
 
! QV - water vapor mixing ratio (kg/kg)
! QC - cloud water mixing ratio (kg/kg)
! QR - rain water mixing ratio (kg/kg)
! QI - cloud ice mixing ratio (kg/kg)
! QS - snow mixing ratio (kg/kg)
! QG - graupel mixing ratio (KG/KG)
! NI - cloud ice number concentration (1/kg)
! NC - Droplet Number concentration (1/kg)
! NS - Snow Number concentration (1/kg)
! NR - Rain Number concentration (1/kg)
! NG - Graupel number concentration (1/kg)
! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
! P - AIR PRESSURE (PA)
! W - VERTICAL AIR VELOCITY (M/S)
! TH - POTENTIAL TEMPERATURE (K)
! PII - exner function - used to convert potential temp to temp
! DZ - difference in height over interface (m)
! DT_IN - model time step (sec)
! ITIMESTEP - time step counter
! RAINNC - accumulated grid-scale precipitation (mm)
! RAINNCV - one time step grid scale precipitation (mm/time step)
! SR - one time step mass ratio of snow to total precip
! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
!................................
! GRID_CLOCK, GRID_ALARMS - parameters to limit radar reflectivity calculation only when needed
! otherwise radar reflectivity calculation every time step is too slow
! only needed for coupling with WRF, see code below for details

! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

! ADDITIONAL INPUT NEEDED BY MICRO
! ********NOTE: WVAR IS SHOULD BE USED IN DROPLET ACTIVATION
! FOR CASES WHEN UPDRAFT IS NOT RESOLVED, EITHER BECAUSE OF
! LOW MODEL RESOLUTION OR CLOUD TYPE
! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

   use input
   use constants , only : rovcp,rd,rp00
   implicit none

!!!   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
!!!                                       ims, ime, jms, jme, kms, kme , &
!!!                                       its, ite, jts, jte, kts, kte
! Temporary changed from INOUT to IN

!!!   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
!!!                          qv, qc, qr, qi, qs, qg, ni, ns, nr, TH, NG
!, effcs, effis

!!!   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
!!!                          pii, p, dz, rho, w !, tke, nctend, nitend,kzh
   real, dimension(ib:ie,jb:je,kb:ke), intent(inout):: &
                          qv, qc, qr, qi, qs, qg, nni, nc, ns, nr, T3D, NG,   &
                          vtc,vtr,vti,vts,vtg,pii
   real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: effc,effi,effs,effr,effg,effis

   real, dimension(ib:ie,jb:je,kb:ke), intent(in):: &
                          rmh, rr,      p, rho, ddz
   real, dimension(ib:ie,jb:je,kb:ke+1), intent(in):: &
                          w
   real, intent(in):: dt_in
   integer, intent(in):: ITIMESTEP

   real, dimension(ib:ie,jb:je,nrain), intent(inout) :: rainnc
   real, dimension(ib:ie,jb:je),       intent(inout) :: prate
   double precision :: tcond,tevac,tevar,train
   real, dimension(ib:ie) :: ruh
   real, dimension(jb:je) :: rvh
   logical, intent(in) :: getdbz,getvt,dorad

      logical, intent(in) :: dotbud,doqbud
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

!!!   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
!!!                          RAINNC, RAINNCV, SR

!   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &  ! GT
!                          refl_10cm

!!!   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

!      TYPE (WRFU_Clock):: grid_clock                  ! GT
!      TYPE (WRFU_Alarm), POINTER:: grid_alarms(:)     ! GT

   ! LOCAL VARIABLES

!!!   REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
!!!                      effi, effs, effr, EFFG

!!!   REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
!!!                      T, WVAR, EFFC

   real, dimension(nk) ::                                                                & 
                            QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, NC_TEND1D,           &
                            NI_TEND1D, NS_TEND1D, NR_TEND1D,                                  &
                            QC1D, QI1D, QR1D, NC1D,NI1D, NS1D, NR1D, QS1D,                    &
                            T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, RHO1D, W1D, WVAR1D,         &
                            EFFC1D, EFFI1D, EFFS1D, EFFR1D,DZ1D,   &
   ! HM ADD GRAUPEL
                            QG_TEND1D, NG_TEND1D, QG1D, NG1D, EFFG1D, &

! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
                            QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN, &

! HM add reflectivity      
                            dbz,effis1d  ! hm 8/20/12
                          
   real PRECPRT1D, SNOWRT1D

   integer I,K,J,n

   real :: DT,rdt,rpii

   real, dimension(1:nk) :: PRE,EVPMS,EVPMG,PCC,FR,FI,FC,FS,FG,  &
                          PRD,PRDS,PRDG,EPRD,EPRDS,EPRDG,MNUCCD, &
                          tbcond,tbevar,tbmelt,tbfrz,tbdep,tbsub,torig,qvorig
   real :: tem1,tem2,tfoo,thfoo,qvfoo,thorig
   double precision, dimension(nj) :: bud1,bud2,bud3,bud4

   logical:: dBZ_tstep ! GT


   ! Initialize tendencies (all set to 0) and transfer
   ! array to local variables
   DT = DT_IN   
   rdt = 1.0/dt

      dBZ_tstep = getdbz

    tem1 = dx*dy*dz

!   DO I=ITS,ITE
!   DO J=JTS,JTE
!   DO K=KTS,KTE
!       T(I,K,J)        = TH(i,k,j)*PII(i,k,j)
!
!! wvar is the ST. DEV. OF sub-grid vertical velocity, used for calculating droplet 
!! activation rates.
!! WVAR CAN BE DERIVED EITHER FROM PREDICTED TKE (AS IN MYJ PBL SCHEME),
!! OR FROM EDDY DIFFUSION COEFFICIENT KZH (AS IN YSU PBL SCHEME),
!! DEPENDING ON THE PARTICULAR pbl SCHEME DRIVER MODEL IS COUPLED WITH
!! NOTE: IF MODEL HAS HIGH ENOUGH RESOLUTION TO RESOLVE UPDRAFTS, WVAR MAY 
!! NOT BE NEEDED 
!
!! currently assign wvar to 0.5 m/s (not coupled with PBL scheme)
!
!       WVAR(I,K,J)     = 0.5
!
!! currently mixing of number concentrations also is neglected (not coupled with PBL schemes)
!
!   END DO
!   END DO
!   END DO

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,QC_TEND1D,QI_TEND1D,QNI_TEND1D,QR_TEND1D,NC_TEND1D,     &
!$omp NI_TEND1D,NS_TEND1D,NR_TEND1D,T_TEND1D,QV_TEND1D,QC1D,QI1D,QS1D,QR1D,   &
!$omp NC1D,NI1D,NS1D,NR1D,QG1D,NG1D,QG_TEND1D,NG_TEND1D,T1D,QV1D,P1D,RHO1D,DZ1D, &
!$omp W1D,WVAR1D,PRE,EVPMS,EVPMG,PCC,dBZ,PRECPRT1D,SNOWRT1D,EFFC1D,EFFI1D,EFFS1D,EFFR1D,  &
!$omp EFFG1D,QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN,FR,FI,FC,FS,FG,tem2,effis1d,  &
!$omp PRD,PRDS,PRDG,EPRD,EPRDS,EPRDG,MNUCCD,tbcond,tbevar,tbmelt,tbfrz,tbdep,tbsub,torig,qvorig,  &
!$omp tfoo,thfoo,qvfoo,thorig,rpii)
   do j=1,nj         ! j loop (north-south)
   bud1(j)=0.0d0
   bud2(j)=0.0d0
   bud3(j)=0.0d0
   bud4(j)=0.0d0
   do i=1,ni         ! i loop (east-west)
   !
   ! Transfer 3D arrays into 1D for microphysical calculations
   !

! hm , initialize 1d tendency arrays to zero

          PRECPRT1D = 0.0

          QC_TEND1D  = 0.
          QI_TEND1D  = 0.
          QNI_TEND1D = 0.
          QR_TEND1D  = 0.
          NC_TEND1D  = 0.
          NI_TEND1D  = 0.
          NS_TEND1D  = 0.
          NR_TEND1D  = 0.
          T_TEND1D   = 0.
          QV_TEND1D  = 0.

          QG_TEND1D  = 0.
          NG_TEND1D  = 0.

          QGSTEN     = 0.
          QRSTEN     = 0.
          QISTEN     = 0.
          QNISTEN    = 0.
          QCSTEN     = 0.

          ! wvar is neglected in cm1
          WVAR1D     = 0.0

          NC1D = 0.0

      if(inum.eq.0)then
        do k=1,nk
          NC1D(k)       = NC(i,j,k)
        enddo
      endif

      do k=1,nk      ! k loop (vertical)

          QC1D(k)       = QC(i,j,k)
          QI1D(k)       = QI(i,j,k)
          QS1D(k)       = QS(i,j,k)
          QR1D(k)       = QR(i,j,k)

          NI1D(k)       = NNI(i,j,k)

          NS1D(k)       = NS(i,j,k)
          NR1D(k)       = NR(i,j,k)
! HM ADD GRAUPEL
          QG1D(K)       = QG(i,j,k)
          NG1D(K)       = NG(i,j,k)

!!!          T1D(k)        = T(i,k,j)
!!!          T1D(k)        = TH(i,j,k)*PII(i,j,k)
          T1D(k)        = T3D(i,j,k)
          QV1D(k)       = QV(i,j,k)
          P1D(k)        = P(i,j,k)
          RHO1D(k)      = RHO(i,j,k)
          DZ1D(k)       = DDZ(i,j,k)
!!!          W1D(k)        = W(i,k,j)
          W1D(k)        = 0.5*( w(i,j,k) + w(i,j,k+1) )

        PRE(k) = 0.0
        EVPMS(k) = 0.0
        EVPMG(k) = 0.0
        PCC(k) = 0.0
        FR(k) = 0.0
        FI(k) = 0.0
        FC(k) = 0.0
        FS(k) = 0.0
        FG(k) = 0.0

        PRD(k) = 0.0
        PRDS(k) = 0.0
        PRDG(k) = 0.0
        EPRD(k) = 0.0
        EPRDS(k) = 0.0
        EPRDG(k) = 0.0
        MNUCCD(k) = 0.0

      end do

      if( dorad )then
      EFFC1D = 0.0
      EFFI1D = 0.0
      EFFS1D = 0.0
      EFFR1D = 0.0
      EFFG1D = 0.0
      effis1d = 0.0
      endif

      tbcond = 0.0
      tbevar = 0.0
      tbmelt = 0.0
      tbfrz = 0.0
      tbdep = 0.0
      tbsub = 0.0

         if( dotbud )then
           do k=1,nk
             torig(k) = T1D(K)
             qvorig(k) = QV1D(K)
           enddo
         endif

      call m2005micro_graupel(QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, NC_TEND1D,            &
       NI_TEND1D, NS_TEND1D, NR_TEND1D,                                                  &
       QC1D, QI1D, QS1D, QR1D, NC1D,NI1D, NS1D, NR1D,                                    &
       T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, RHO1D, DZ1D, W1D, WVAR1D,                   &
       PRECPRT1D,SNOWRT1D,                                                               &
       EFFC1D,EFFI1D,EFFS1D,EFFR1D,DT,RDT,                                               &
                                              1,ni ,   1,nj ,   1,nk ,                   &
                                              1,ni ,   1,nj ,   1,nk ,                   &
                                    QG_TEND1D,NG_TEND1D,QG1D,NG1D,EFFG1D, &
! ADD SEDIMENTATION TENDENCIES
                                  QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN,     &
                                  PRE,EVPMS,EVPMG,PCC,FR,FI,FC,FS,FG,      &
                                  PRD,PRDS,PRDG,EPRD,EPRDS,EPRDG,MNUCCD,   &
                                  tbcond,tbevar,tbmelt,tbfrz,tbdep,tbsub,  &
                                  eqtset,dorad,tsmall,effis1d,testcase)


   !
   ! Transfer 1D arrays back into 3D arrays
   !
      do k=1,nk

! hm, add tendencies to update global variables 
! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

          QC(i,j,k)        = QC1D(k)
          QI(i,j,k)        = QI1D(k)
          QS(i,j,k)        = QS1D(k)
          QR(i,j,k)        = QR1D(k)
          NNI(i,j,k)        = NI1D(k)
          NS(i,j,k)        = NS1D(k)          
          NR(i,j,k)        = NR1D(k)
	  QG(i,j,k)        = QG1D(K)
          NG(i,j,k)        = NG1D(K)

!!!          T(i,k,j)         = T1D(k)
!!!          TH(I,K,J)        = T(i,k,j)/PII(i,k,j) ! CONVERT TEMP BACK TO POTENTIAL TEMP
!!!          TH(i,j,k)        = T1D(k)/PII(i,j,k) ! CONVERT TEMP BACK TO POTENTIAL TEMP
          T3D(i,j,k)        = T1D(k)
          QV(i,j,k)        = QV1D(k)

      end do

      if( dorad )then
        if( radopt.eq.1 )then
          ! for goddard radiation scheme:
          do k=1,nk
! hm 8/20/12
            EFFC(i,j,k)      = EFFC1D(k)
            EFFI(i,j,k)      = EFFI1D(k)
            EFFS(i,j,k)      = EFFS1D(k)
            EFFR(i,j,k)      = EFFR1D(k)
            EFFG(i,j,k)      = EFFG1D(K)
            effis(i,j,k)      = effis1D(K)  ! combined cloud ice plus snow, with correct weighting
            ! GHB, 130904: still seeing blowups in goddard code, so add an upper limit:
            EFFI(i,j,k) = min( EFFI(i,j,k) , 400.0 )
            EFFS(i,j,k) = min( EFFS(i,j,k) , 400.0 )
            effis(i,j,k) = min( effis(i,j,k) , 400.0 )
          end do
        elseif( radopt.eq.2 )then
          ! for rrtmg radiation scheme:
          do k=1,nk
            EFFC(i,j,k) = EFFC1D(k)*1.0e-6
            EFFI(i,j,k) = min( 140.0 , EFFI1D(k) )*1.0e-6
            EFFS(i,j,k) = EFFS1D(k)*1.0e-6
          enddo
        endif
      endif

      if(inum.eq.0)then
        do k=1,nk
          NC(i,j,k)        = NC1D(k)
        enddo
      endif

      if( getvt )then
        do k=1,nk
          vtc(i,j,k) = fc(k)
          vtr(i,j,k) = fr(k)
          vti(i,j,k) = fi(k)
          vts(i,j,k) = fs(k)
          vtg(i,j,k) = fg(k)
        enddo
      endif

! hm modified so that m2005 precip variables correctly match wrf precip variables
!!!      RAINNC(i,j) = RAINNC(I,J)+PRECPRT1D
!!!      RAINNCV(i,j) = PRECPRT1D
!!!      SR(i,j) = SNOWRT1D/(PRECPRT1D+1.E-12)


         do n = 1,nrain
           RAINNC(i,j,n) = RAINNC(i,j,n) + 0.1*PRECPRT1D
                                        ! mult by 100 cm/m  (convert from m to cm)
                                        ! div by 1000 kg/m^3  (density of water)
         enddo
         prate(i,j) = PRECPRT1D*rdt
         if( axisymm.eq.1 )then
           bud1(j) = bud1(j) + PRECPRT1D*ruh(i)*rvh(j)*dx*dy*rr(i,j,1)/rho(i,j,1)
         else
           bud1(j) = bud1(j) + PRECPRT1D*ruh(i)*rvh(j)*dx*dy
         endif

         do k=1,nk
           tem2 = dt*rr(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*tem1
           bud2(j) = bud2(j) - tem2*( PRE(K)+EVPMS(K)+EVPMG(K) )
           bud3(j) = bud3(j) + tem2*max(PCC(K),0.0)
           bud4(j) = bud4(j) - tem2*min(PCC(K),0.0)
         enddo

         tbuds:  &
         if( dotbud )then

           eqs2:  &
           if(eqtset.eq.2)then

           if( td_cond.ge.1 )then
             do k=1,nk
               rpii = 1.0/pii(i,j,k)
               thorig = torig(k)*rpii
               tfoo = torig(k)+dt*max(tbcond(k),0.0)
               qvfoo = qvorig(k)-dt*max(PCC(k),0.0)
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_cond) = (thfoo-thorig)*rdt
             enddo
           endif
           if( td_evac.ge.1 )then
             do k=1,nk
               tfoo = torig(k)+dt*min(tbcond(k),0.0)
               qvfoo = qvorig(k)-dt*min(PCC(k),0.0)
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_evac) = (thfoo-thorig)*rdt
             enddo
           endif
           if( td_evar.ge.1 )then
             do k=1,nk
               tfoo = torig(k)+dt*tbevar(k)
               qvfoo = qvorig(k)+dt*( -PRE(K)-EVPMS(K)-EVPMG(K) )
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_evar) = (thfoo-thorig)*rdt
             enddo
           endif
           if( td_melt.ge.1 )then
             do k=1,nk
               tfoo = torig(k)+dt*tbmelt(k)
               qvfoo = qvorig(k)
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_melt) = (thfoo-thorig)*rdt
             enddo
           endif
           if( td_frz .ge.1 )then
             do k=1,nk
               tfoo = torig(k)+dt*tbfrz(k)
               qvfoo = qvorig(k)
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_frz) = (thfoo-thorig)*rdt
             enddo
           endif
           if( td_dep .ge.1 )then
             do k=1,nk
               tfoo = torig(k)+dt*tbdep(k)
               qvfoo = qvorig(k)+dt*( -PRD(k)-PRDS(k)-PRDG(k)-MNUCCD(K) )
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_dep) = (thfoo-thorig)*rdt
             enddo
           endif
           if( td_subl.ge.1 )then
             do k=1,nk
               tfoo = torig(k)+dt*tbsub(k)
               qvfoo = qvorig(k)+dt*( -EPRD(k)-EPRDS(k)-EPRDG(k) )
               thfoo = tfoo/( (rho(i,j,k)*(rd+rv*qvfoo)*tfoo*rp00)**rovcp )
               tdiag(i,j,k,td_subl) = (thfoo-thorig)*rdt
             enddo
           endif

           else  eqs2

           if( td_cond.ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_cond) = max(tbcond(k),0.0)/pii(i,j,k)
             enddo
           endif
           if( td_evac.ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_evac) = min(tbcond(k),0.0)/pii(i,j,k)
             enddo
           endif
           if( td_evar.ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_evar) = tbevar(k)/pii(i,j,k)
             enddo
           endif
           if( td_melt.ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_melt) = tbmelt(k)/pii(i,j,k)
             enddo
           endif
           if( td_frz .ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_frz)  = tbfrz(k)/pii(i,j,k)
             enddo
           endif
           if( td_subl.ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_subl)  = tbsub(k)/pii(i,j,k)
             enddo
           endif
           if( td_dep .ge.1 )then
             do k=1,nk
               tdiag(i,j,k,td_dep)  = tbdep(k)/pii(i,j,k)
             enddo
           endif

           endif  eqs2
         endif  tbuds

         qbuds:  &
         if( doqbud )then
         if( qd_cond.ge.1 )then
           do k=1,nk
             qdiag(i,j,k,qd_cond) = -max(PCC(k),0.0)
           enddo
         endif
         if( qd_evac.ge.1 )then
           do k=1,nk
             qdiag(i,j,k,qd_evac) = -min(PCC(k),0.0)
           enddo
         endif
         if( qd_evar.ge.1 )then
           do k=1,nk
             qdiag(i,j,k,qd_evar) = -PRE(K)-EVPMS(K)-EVPMG(K)
           enddo
         endif
         if( qd_dep.ge.1 )then
           do k=1,nk
             qdiag(i,j,k,qd_dep) = -PRD(k)-PRDS(k)-PRDG(k)-MNUCCD(K)
           enddo
         endif
         if( qd_subl.ge.1 )then
           do k=1,nk
             qdiag(i,j,k,qd_subl) = -EPRD(k)-EPRDS(k)-EPRDG(k)
           enddo
         endif
         endif  qbuds

!!!         if( nout3d.ge.1 )then
!!!           do k=1,nk
!!!             out3d(i,j,k,1) = PRE(K)
!!!             out3d(i,j,k,2) = EVPMS(K)
!!!             out3d(i,j,k,3) = EVPMG(K)
!!!           enddo
!!!         endif

! add reflectivity calculations
! only calculate if logical parameter dbz_tstep = .true.

         if (dBZ_tstep .and. qd_dbz.ge.1) then
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          call calc_refl10cm (qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ,    &
                      1, nk, i, j, nr1d, ns1d, ng1d, rho1d)
          do k = 1, nk
             qdiag(i,j,k,qd_dbz) = dBZ(k)
          enddo
          if(timestats.ge.1) time_dbz=time_dbz+mytime()
         endif

   end do
   end do   

   do j=1,nj
     train=train+bud1(j)
     tevar=tevar+bud2(j)
     tcond=tcond+bud3(j)
     tevac=tevac+bud4(j)
   enddo

end subroutine MP_GRAUPEL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine M2005MICRO_GRAUPEL(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,NC3DTEN,    &
       NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NC3D,NI3D,NS3D,NR3D,         &
       T3DTEN,QV3DTEN,T3D,QV3D,PRES,RHO3D,DZQ,W3D,WVAR,PRECRT,SNOWRT,            &
       EFFC,EFFI,EFFS,EFFR,DT,RDT,                                               &
                                            IMS,IME, JMS,JME, KMS,KME,           &
                                            ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
                        QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN,   &
                        PRE,EVPMS,EVPMG,PCC,FR,FI,FC,FS,FG,      &
                        PRD,PRDS,PRDG,EPRD,EPRDS,EPRDG,MNUCCD,   &
                        tbcond,tbevar,tbmelt,tbfrz,tbdep,tbsub,  &
                        eqtset,dorad,tsmall,effis,testcase)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
! MORRISON ET AL. 2005 JAS AND MORRISON ET AL. 2009 MWR

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

! CODE STRUCTURE: MAIN SUBROUTINE IS 'M2005MICRO_GRAUPEL'. ALSO INCLUDED IN THIS FILE IS
! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
! 'FUNCTION GAMMA'.

! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! DECLARATIONS

      implicit none

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
! DEFINE ARRAY SIZES

! INPUT NUMBER OF GRID CELLS

! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
      integer, intent( in)  :: IMS,IME, JMS,JME, KMS,KME,          &
                               ITS,ITE, JTS,JTE, KTS,KTE

      real, dimension(KMS:KME) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
      real, dimension(KMS:KME) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
      real, dimension(KMS:KME) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
      real, dimension(KMS:KME) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
      real, dimension(KMS:KME) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
      real, dimension(KMS:KME) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
      real, dimension(KMS:KME) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
      real, dimension(KMS:KME) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)
      real, dimension(KMS:KME) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
      real, dimension(KMS:KME) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
      real, dimension(KMS:KME) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
      real, dimension(KMS:KME) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)
      real, dimension(KMS:KME) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
      real, dimension(KMS:KME) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
      real, dimension(KMS:KME) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
      real, dimension(KMS:KME) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)
      real, dimension(KMS:KME) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
      real, dimension(KMS:KME) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
      real, dimension(KMS:KME) ::  T3D                ! TEMPERATURE (K)
      real, dimension(KMS:KME) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
      real, dimension(KMS:KME) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
      real, dimension(KMS:KME) ::  RHO3D              ! DENSITY (KG/M^3)
      real, dimension(KMS:KME) ::  DZQ                ! DIFFERENCE in HEIGHT ACROSS LEVEL (m)
      real, dimension(KMS:KME) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
      real, dimension(KMS:KME) ::  WVAR               ! SUB-GRID VERTICAL VELOCITY (M/S)

! HM ADDED GRAUPEL VARIABLES
      real, dimension(KMS:KME) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
      real, dimension(KMS:KME) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
      real, dimension(KMS:KME) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
      real, dimension(KMS:KME) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

      real, dimension(KMS:KME) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
      real, dimension(KMS:KME) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
      real, dimension(KMS:KME) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
      real, dimension(KMS:KME) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
      real, dimension(KMS:KME) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)      

      integer, intent(in) :: eqtset
      logical, intent(in) :: dorad
      real, intent(in) :: tsmall

! OUTPUT VARIABLES

        real PRECRT                ! TOTAL PRECIP PER TIME STEP (mm)
        real SNOWRT                ! SNOW PER TIME STEP (mm)

        real, dimension(KMS:KME) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
        real, dimension(KMS:KME) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
        real, dimension(KMS:KME) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
        real, dimension(KMS:KME) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
        real, dimension(KMS:KME) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)
! hm 8/20/12
        real, dimension(KMS:KME) ::   effis           ! GRAUPEL EFFECTIVE RADIUS (MICRON)

! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

        real, intent(in) :: DT,RDT         ! MODEL TIME STEP (SEC)
        integer, intent(in) :: testcase

!.....................................................................................................
! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
! REST OF THE MODEL.

! SIZE PARAMETER VARIABLES

     real, dimension(KMS:KME) :: LAMC          ! SLOPE parameter FOR DROPLETS (M-1)
     real, dimension(KMS:KME) :: LAMI          ! SLOPE parameter FOR CLOUD ICE (M-1)
     real, dimension(KMS:KME) :: LAMS          ! SLOPE parameter FOR SNOW (M-1)
     real, dimension(KMS:KME) :: LAMR          ! SLOPE parameter FOR RAIN (M-1)
     real, dimension(KMS:KME) :: LAMG          ! SLOPE parameter FOR GRAUPEL (M-1)
     real, dimension(KMS:KME) :: CDIST1        ! PSD parameter FOR DROPLETS
     real, dimension(KMS:KME) :: N0I           ! INTERCEPT parameter FOR CLOUD ICE (KG-1 M-1)
     real, dimension(KMS:KME) :: N0S           ! INTERCEPT parameter FOR SNOW (KG-1 M-1)
     real, dimension(KMS:KME) :: N0RR          ! INTERCEPT parameter FOR RAIN (KG-1 M-1)
     real, dimension(KMS:KME) :: N0G           ! INTERCEPT parameter FOR GRAUPEL (KG-1 M-1)
     real, dimension(KMS:KME) :: PGAM          ! SPECTRAL SHAPE parameter FOR DROPLETS

! MICROPHYSICAL PROCESSES

     real, dimension(KMS:KME) ::  NSUBC     ! LOSS OF NC DURING EVAP
     real, dimension(KMS:KME) ::  NSUBI     ! LOSS OF NI DURING SUB.
     real, dimension(KMS:KME) ::  NSUBS     ! LOSS OF NS DURING SUB.
     real, dimension(KMS:KME) ::  NSUBR     ! LOSS OF NR DURING EVAP
     real, dimension(KMS:KME), intent(inout) ::  PRD       ! DEP CLOUD ICE
     real, dimension(KMS:KME), intent(inout) ::  PRE       ! EVAP OF RAIN
     real, dimension(KMS:KME), intent(inout) ::  PRDS      ! DEP SNOW
     real, dimension(KMS:KME) ::  NNUCCC    ! CHANGE N DUE to CONTACT FREEZ DROPLETS
     real, dimension(KMS:KME) ::  MNUCCC    ! CHANGE Q DUE to CONTACT FREEZ DROPLETS
     real, dimension(KMS:KME) ::  PRA       ! ACCRETION DROPLETS BY RAIN
     real, dimension(KMS:KME) ::  PRC       ! AUTOCONVERSION DROPLETS
     real, intent(inout), dimension(KMS:KME) ::  PCC       ! COND/EVAP DROPLETS
     real, dimension(KMS:KME) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     real, dimension(KMS:KME), intent(inout) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     real, dimension(KMS:KME) ::  MNUCCR    ! CHANGE Q DUE to CONTACT FREEZ RAIN
     real, dimension(KMS:KME) ::  NNUCCR    ! CHANGE N DUE to CONTACT FREEZ RAIN
     real, dimension(KMS:KME) ::  NPRA      ! CHANGE in N DUE to DROPLET ACC BY RAIN
     real, dimension(KMS:KME) ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
     real, dimension(KMS:KME) ::  NSAGG     ! SELF-COLLECTION OF SNOW
     real, dimension(KMS:KME) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
     real, dimension(KMS:KME) ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
     real, dimension(KMS:KME) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE
     real, dimension(KMS:KME) ::  PRCI      ! CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW
     real, dimension(KMS:KME) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
     real, dimension(KMS:KME) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
     real, dimension(KMS:KME) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     real, dimension(KMS:KME) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     real, dimension(KMS:KME) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     real, dimension(KMS:KME) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
     real, dimension(KMS:KME) ::  NMULTS    ! ICE MULT DUE to RIMING DROPLETS BY SNOW
     real, dimension(KMS:KME) ::  NMULTR    ! ICE MULT DUE to RIMING RAIN BY SNOW
     real, dimension(KMS:KME) ::  QMULTS    ! CHANGE Q DUE to ICE MULT DROPLETS/SNOW
     real, dimension(KMS:KME) ::  QMULTR    ! CHANGE Q DUE to ICE RAIN/SNOW
     real, dimension(KMS:KME) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
     real, dimension(KMS:KME) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
     real, dimension(KMS:KME) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
     real, dimension(KMS:KME) ::  PSMLT     ! CHANGE Q MELTING SNOW to RAIN
     real, dimension(KMS:KME), intent(inout) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
     real, dimension(KMS:KME) ::  NSMLTS    ! CHANGE N MELTING SNOW
     real, dimension(KMS:KME) ::  NSMLTR    ! CHANGE N MELTING SNOW to RAIN
! HM ADDED 12/13/06
     real, dimension(KMS:KME) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
     real, dimension(KMS:KME) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
     real, dimension(KMS:KME) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
     real, dimension(KMS:KME) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED to SNOW
     real, dimension(KMS:KME) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED to SNOW
     real, dimension(KMS:KME) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED to SNOW
     real, dimension(KMS:KME), intent(inout) ::  EPRD      ! SUBLIMATION CLOUD ICE
     real, dimension(KMS:KME), intent(inout) ::  EPRDS     ! SUBLIMATION SNOW
! HM ADDED GRAUPEL PROCESSES
     real, dimension(KMS:KME) ::  PRACG    ! CHANGE in Q COLLECTION RAIN BY GRAUPEL
     real, dimension(KMS:KME) ::  PSACWG    ! CHANGE in Q COLLECTION DROPLETS BY GRAUPEL
     real, dimension(KMS:KME) ::  PGSACW    ! CONVERSION Q to GRAUPEL DUE to COLLECTION DROPLETS BY SNOW
     real, dimension(KMS:KME) ::  PGRACS    ! CONVERSION Q to GRAUPEL DUE to COLLECTION RAIN BY SNOW
     real, dimension(KMS:KME), intent(inout) ::  PRDG    ! DEP OF GRAUPEL
     real, dimension(KMS:KME), intent(inout) ::  EPRDG    ! SUB OF GRAUPEL
     real, dimension(KMS:KME), intent(inout) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     real, dimension(KMS:KME) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
     real, dimension(KMS:KME) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
     real, dimension(KMS:KME) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     real, dimension(KMS:KME) ::  NSCNG    ! CHANGE N CONVERSION to GRAUPEL DUE to COLLECTION DROPLETS BY SNOW
     real, dimension(KMS:KME) ::  NGRACS    ! CHANGE N CONVERSION to GRAUPEL DUE to COLLECTION RAIN BY SNOW
     real, dimension(KMS:KME) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
     real, dimension(KMS:KME) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL to RAIN
     real, dimension(KMS:KME) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
     real, dimension(KMS:KME) ::  PSACR    ! CONVERSION DUE to COLL OF SNOW BY RAIN
     real, dimension(KMS:KME) ::  NMULTG    ! ICE MULT DUE to ACC DROPLETS BY GRAUPEL
     real, dimension(KMS:KME) ::  NMULTRG    ! ICE MULT DUE to ACC RAIN BY GRAUPEL
     real, dimension(KMS:KME) ::  QMULTG    ! CHANGE Q DUE to ICE MULT DROPLETS/GRAUPEL
     real, dimension(KMS:KME) ::  QMULTRG    ! CHANGE Q DUE to ICE MULT RAIN/GRAUPEL

! TIME-VARYING ATMOSPHERIC PARAMETERS

     real, dimension(KMS:KME) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
     real, dimension(KMS:KME) ::   EVS   ! SATURATION VAPOR PRESSURE
     real, dimension(KMS:KME) ::   EIS   ! ICE SATURATION VAPOR PRESSURE
     real, dimension(KMS:KME) ::   QVS   ! SATURATION MIXING RATIO
     real, dimension(KMS:KME) ::   QVI   ! ICE SATURATION MIXING RATIO
     real, dimension(KMS:KME) ::   QVQVS ! SAUTRATION RATIO
     real, dimension(KMS:KME) ::   QVQVSI! ICE SATURAION RATIO
     real, dimension(KMS:KME) ::   DV    ! DIFFUSIVITY OF WATER VAPOR in AIR
     real, dimension(KMS:KME) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
     real, dimension(KMS:KME) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
     real, dimension(KMS:KME) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
     ! GHB
     real :: cvm
     real, dimension(KMS:KME) ::   ecnd,edep,efrz
     ! GHB
     real, dimension(KMS:KME) ::   MU    ! VISCOCITY OF AIR
     real, dimension(KMS:KME) ::   SC    ! SCHMIDT NUMBER
     real, dimension(KMS:KME) ::   XLF   ! LATENT HEAT OF FREEZING
     real, dimension(KMS:KME) ::   RHO   ! AIR DENSITY
     real, dimension(KMS:KME) ::   AB    ! CORRECTION to CONDENSATION RATE DUE to LATENT HEATING
     real, dimension(KMS:KME) ::   ABI    ! CORRECTION to DEPOSITION RATE DUE to LATENT HEATING

! TIME-VARYING MICROPHYSICS PARAMETERS

     real, dimension(KMS:KME) ::   DAP    ! DIFFUSIVITY OF AEROSOL
     real    NACNT                    ! NUMBER OF CONTACT in
     real    FMULT                    ! TEMP.-DEP. parameter FOR RIME-SPLINTERING
     real    COFFI                    ! ICE AUTOCONVERSION parameter

! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

      real, dimension(KMS:KME), intent(inout) :: FR,FI,FC,FS,FG
      real, dimension(KMS:KME) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
      real UNI, UMI,UMR
      real, dimension(KMS:KME) ::    FNI,FNG
      real RGVM
      real, dimension(KMS:KME) ::   FALOUTR,FALOUTI,FALOUTNI
      real FALTNDR,FALTNDI,FALTNDNI,RHO2
      real, dimension(KMS:KME) ::   DUMQS,DUMFNS
      real UMS,UNS
      real, dimension(KMS:KME) ::   FNS, FALOUTS,FALOUTNS,FALOUTG,FALOUTNG
      real FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
      real, dimension(KMS:KME) ::    DUMC,DUMFNC
      real UNC,UMC,UNG,UMG
      real, dimension(KMS:KME) ::   FALOUTC,FALOUTNC
      real FALTNDC,FALTNDNC
      real, dimension(KMS:KME) ::   FNC,DUMFNR,FALOUTNR
      real FALTNDNR
      real, dimension(KMS:KME) ::   FNR(KMS:KME)

! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

      real, dimension(KMS:KME) ::    AIN,ARN,ASN,ACN,AGN

      real, dimension(KMS:KME), intent(inout) :: tbcond,tbevar,tbmelt,tbfrz,tbdep,tbsub

! EXTERNAL FUNCTION CALL RETURN VARIABLES

!      REAL GAMMA,      ! EULER GAMMA FUNCTION
!      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
!      REAL DERF1        ! ERROR FUNCTION

! DUMMY VARIABLES

     real DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS
     real DUMMQS,DUMQG,SAVEQC,SAVET

! PROGNOSTIC SUPERSATURATION

     real DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
     real DQSIDT   ! CHANGE in ICE SAT. MIXING RAT. WITH T
     real EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
     real EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
     real EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
     real EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

! NEW DROPLET ACTIVATION VARIABLES
     real TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
     real TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
     real TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
     real TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
     real TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
     real DUMACT,DUM3

! COUNTING/INDEX VARIABLES

! V1.3 DIFFERENT NSTEP FOR EACH SPECIES
     integer K,N,NSTEPR,NSTEPI,NSTEPS,NSTEPC,NSTEPG ! ,I

! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN, 
!               = 1, HYDROMETEORS IN COLUMN

      integer LTRUE

! DROPLET ACTIVATION/FREEZING AEROSOL


     real    CT      ! DROPLET ACTIVATION parameter
     real    TEMP1   ! DUMMY TEMPERATURE
     real    SAT1    ! DUMMY SATURATION
     real    SIGVL   ! SURFACE TENSION LIQ/VAPOR
     real    KEL     ! KELVIN parameter
     real    KC2     ! TOTAL ICE NUCLEATION RATE

       real CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

! MORE WORKING/DUMMY VARIABLES

     real DUMQI,DUMNI,DC0,DS0,DG0
     real DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

! EFFECTIVE VERTICAL VELOCITY  (M/S)
     real WEF

! WORKING PARAMETERS FOR ICE NUCLEATION

      real ANUC,BNUC

! WORKING PARAMETERS FOR AEROSOL ACTIVATION

        real AACT,GAMM,GG,PSI,ETA1,ETA2,SM1,SM2,SMAX,UU1,UU2,ALPHA

! DUMMY SIZE DISTRIBUTION PARAMETERS

        real DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN

        integer IDROP

! v1.4
! new variables for seifert and beheng warm rain scheme
      real, dimension(KMS:KME) :: nu
      integer dumii

      real :: dttmp
      real, dimension(KMS:KME) :: temr,temz

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! SET LTRUE INITIALLY TO 0

         LTRUE = 0

      if( dorad )then
! V.13 initialize effective radii to default values (from P. Blossey)
! hm 8/20/12
         effc(kts:kte) = 25.         
         effi(kts:kte) = 25.         
         effs(kts:kte) = 25.         
         effr(kts:kte) = 25.         
         effg(kts:kte) = 25.         
	 effis(kts:kte) = 25.
      endif
!---------------------------------!

! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
         kloop1:  &
         do K = KTS,KTE

! LATENT HEAT OF VAPORATION

!!!            XXLV(K) = 3.1484E6-2370.*T3D(K)
            ! GHB, 130722:
            XXLV(K) = lv1-lv2*T3D(K)

! LATENT HEAT OF SUBLIMATION

!!!            XXLS(K) = 3.15E6-2370.*T3D(K)+0.3337E6
            ! GHB, 130722:
            XXLS(K) = ls1-ls2*T3D(K)

            CPM(K) = CP*(1.+0.887*QV3D(K))

            ! GHB
            if(eqtset.eq.2)then
              cvm = 1.0/( cv+cvv*qv3d(k)+cpl*(qc3d(k)+qr3d(k))   &
                                        +cpi*(qi3d(k)+qni3d(k)+qg3d(k)) )
              ecnd(k) = (XXLV(K)-RV*T3D(K))*cvm
              edep(k) = (XXLS(K)-RV*T3D(K))*cvm
              efrz(k) = (XXLS(K)-XXLV(K))*cvm
            else
              ecnd(k) = XXLV(K)/CPM(K)
              edep(k) = XXLS(K)/CPM(K)
              efrz(k) = (XXLS(K)-XXLV(K))/CPM(K)
            endif
            ! GHB

! SATURATION VAPOR PRESSURE AND MIXING RATIO

!!!            EVS(K) = POLYSVP(T3D(K),0)   ! PA
!!!            EIS(K) = POLYSVP(T3D(K),1)   ! PA
            EVS(K) = min( PRES(K)*0.5 , 611.2 * EXP( 17.67 * ( T3D(K) - 273.15 ) / ( T3D(K) - 29.65 ) ) )
            EIS(K) = min( PRES(K)*0.5 , 611.2 * EXP( 21.8745584 * ( T3D(K) - 273.15 ) / ( T3D(K) - 7.66 ) ) )

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            if (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AIR DENSITY

!!!            RHO(K) = PRES(K)/(R*T3D(K))
            RHO(K) = RHO3D(K)

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! V1.3, change limit from 10^-7 to 10^-6
! V1.7 7/9/09 change limit from 10^-6 to 10^-8
! this improves reflectivity at low mixing ratios
! ghb 160428:  use QTHR

             if (QVQVS(K).LT.0.9) then
               if (QR3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QR3D(K)
!!!                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  T3D(K)=T3D(K)-QR3D(K)*ecnd(k)
                  QR3D(K)=0.
               end if
               if (QC3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QC3D(K)
!!!                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  T3D(K)=T3D(K)-QC3D(K)*ecnd(k)
                  QC3D(K)=0.
               end if
             end if

             if (QVQVSI(K).LT.0.9) then
               if (QI3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QI3D(K)
!!!                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  T3D(K)=T3D(K)-QI3D(K)*edep(k)
                  QI3D(K)=0.
               end if
               if (QNI3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QNI3D(K)
!!!                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  T3D(K)=T3D(K)-QNI3D(K)*edep(k)
                  QNI3D(K)=0.
               end if
               if (QG3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QG3D(K)
!!!                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  T3D(K)=T3D(K)-QG3D(K)*edep(k)
                  QG3D(K)=0.
               end if
             end if

! HEAT OF FUSION

            XLF(K) = XXLS(K)-XXLV(K)

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       if (QC3D(K).LT.QSMALL) then
         QC3D(K) = 0.
         NC3D(K) = 0.
       end if
       if (QR3D(K).LT.QSMALL) then
         QR3D(K) = 0.
         NR3D(K) = 0.
       end if
       if (QI3D(K).LT.QSMALL) then
         QI3D(K) = 0.
         NI3D(K) = 0.
       end if
       if (QNI3D(K).LT.QSMALL) then
         QNI3D(K) = 0.
         NS3D(K) = 0.
       end if
       if (QG3D(K).LT.QSMALL) then
         QG3D(K) = 0.
         NG3D(K) = 0.
       end if

! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

      QRSTEN(K) = 0.
      QISTEN(K) = 0.
      QNISTEN(K) = 0.
      QCSTEN(K) = 0.
      QGSTEN(K) = 0.

!..................................................................
! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

! DYNAMIC VISCOSITY OF AIR

            MU(K) = 1.496E-6*T3D(K)**1.5/(T3D(K)+120.)

! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

            DUM = (RHOSU/RHO(K))**0.54

! v3 5/27/11
!            AIN(K) = DUM*AI
! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction 
! hm fix 11/18/11
            AIN(K) = (RHOSU/RHO(K))**0.35*AI
            ARN(K) = DUM*AR
            ASN(K) = DUM*AS
!            ACN(K) = DUM*AC
! AA revision 4/1/11: temperature-dependent Stokes fall speed
            ACN(K) = G*RHOW/(18.*MU(K))
! HM ADD GRAUPEL 8/28/06
            AGN(K) = DUM*AG

! V1.7
! bug fix 7/10/09 
!hm 4/15/09 bug fix, initialize lami to prevent later division by zero
            LAMI(K)=0.

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
! FOR THIS LEVEL

            if (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) then
                 if (T3D(K).LT.TMELT.AND.QVQVSI(K).LT.0.999) GOTO 200
                 if (T3D(K).GE.TMELT.AND.QVQVS(K).LT.0.999) GOTO 200
            end if

! THERMAL CONDUCTIVITY FOR AIR

! v3 5/27/11
            KAP(K) = 1.414E3*MU(K)

! DIFFUSIVITY OF WATER VAPOR

            DV(K) = 8.794E-5*T3D(K)**1.81/PRES(K)

! SCHMIT NUMBER

! v3 5/27/11
            SC(K) = MU(K)/(RHO(K)*DV(K))

! PSYCHOMETIC CORRECTIONS

! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

            DUM = (RV*T3D(K)**2)

            DQSDT = XXLV(K)*QVS(K)/DUM
            DQSIDT =  XXLS(K)*QVI(K)/DUM

            ABI(K) = 1.+DQSIDT*XXLS(K)/CPM(K)
            AB(K) = 1.+DQSDT*XXLV(K)/CPM(K)

! 
!.....................................................................
!.....................................................................
! CASE FOR TEMPERATURE ABOVE FREEZING

            if (T3D(K).GE.TMELT) then

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         if (INUM.EQ.1) then
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         end if

! GET SIZE DISTRIBUTION PARAMETERS

! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
       if (QNI3D(K).LT.1.E-6) then
          QR3D(K)=QR3D(K)+QNI3D(K)
          NR3D(K)=NR3D(K)+NS3D(K)
!!!          T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
          T3D(K)=T3D(K)-QNI3D(K)*efrz(k)
          QNI3D(K) = 0.
          NS3D(K) = 0.
       end if
       if (QG3D(K).LT.1.E-6) then
          QR3D(K)=QR3D(K)+QG3D(K)
          NR3D(K)=NR3D(K)+NG3D(K)
!!!          T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
          T3D(K)=T3D(K)-QG3D(K)*efrz(k)
          QG3D(K) = 0.
          NG3D(K) = 0.
       end if

       if (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8) GOTO 300

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NS3D(K) = MAX(0.,NS3D(K))
      if(INUM.EQ.0)  &
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! RAIN

      if (QR3D(K).GE.QSMALL) then
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMR(K).LT.LAMMINR) then

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      else if (LAMR(K).GT.LAMMAXR) then
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      end if
      end if

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      if (QC3D(K).GE.QSMALL) then

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)
! v1.4
! interpolate
         dumii=int(pgam(k))
         nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam(k)-real(dumii))

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      if (LAMC(K).LT.LAMMIN) then
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      else if (LAMC(K).GT.LAMMAX) then
      LAMC(K) = LAMMAX

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      end if

      end if

!......................................................................
! SNOW

      if (QNI3D(K).GE.QSMALL) then
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMS(K).LT.LAMMINS) then
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      else if (LAMS(K).GT.LAMMAXS) then

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      end if
      end if

!......................................................................
! GRAUPEL

      if (QG3D(K).GE.QSMALL) then
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! ADJUST VARS

      if (LAMG(K).LT.LAMMING) then
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      else if (LAMG(K).GT.LAMMAXG) then

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      end if
      end if

!.....................................................................
! ZERO OUT PROCESS RATES

            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            NSMLTS(K) = 0.
            NSMLTR(K) = 0.
            EVPMS(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            NSUBC(K) = 0.
            NSUBR(K) = 0.
            PRACG(K) = 0.
            NPRACG(K) = 0.
            PSMLT(K) = 0.
            PGMLT(K) = 0.
            EVPMG(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         if (QC3D(K).GE.1.E-6) then

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

            if (IRAIN.EQ.0) then

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

            else if (IRAIN.EQ.1) then

! v1.4
! replace with seifert and beheng

        dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
        dum1 = 600.*dum**0.68*(1.-dum**0.68)**3

        prc(k) = 9.44e9/(20.*2.6e-7)* &
        (nu(k)+2.)*(nu(k)+4.)/(nu(k)+1.)**2* &
        (rho(k)*qc3d(k)/1000.)**4/(rho(k)*nc3d(k)/1.e6)**2* &
        (1.+dum1/(1.-dum)**2)*1000./rho(k)

        nprc(k) = prc(k)*2./2.6e-7*1000.
        nprc1(k) = 0.5*nprc(k)

        end if
         end if

!.......................................................................
! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING

         if (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) then

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! hm fix, 3/4/13
! for above freezing conditions to get accelerated melting of snow,
! we need collection of rain by snow (following Lin et al. 1983)
!            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
!                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
!                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
!                  (5./(LAMS(K)**3*LAMR(K))+                    &
!                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
!                  0.5/(LAMS(K)*LAMR(K)**3)))

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

! v3 5/27/11 npracs no longer used
!            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
!                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
!                (1./(LAMR(K)**3*LAMS(K))+                      &
!                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
!                 1./(LAMR(K)*LAMS(K)**3))

         end if

! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
! ASSUME SHED DROPS ARE 1 MM IN SIZE

         if (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) then

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

! ASSUME 1 MM DROPS ARE SHED, GET NUMBER CONC (KG-1) SHED PER SEC

            DUM = PRACG(K)/5.2E-7

! GET NUMBER CONC OF RAIN DROPS COLLECTED

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! hm 7/15/13, remove limit so that the number of collected drops can smaller than
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
            NPRACG(K)=NPRACG(K)-DUM

	    end if

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         if (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) then

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

            if (IRAIN.EQ.0) then

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

           else if (IRAIN.EQ.1) then

! v1.4
! seifert and beheng (2001) formulation

           dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
           dum1 = (dum/(dum+5.e-4))**4
           pra(k) = 5.78e3*rho(k)/1000.*qc3d(k)*qr3d(k)*dum1
           npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
           (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)

         end if
         end if
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

! v1.4, replace with seifert and beheng (2001)

         if (QR3D(K).GE.1.E-8) then
! v1.4
! seifert and beheng
! include breakup, V2.1
            dum1=D_b
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
            nragg(k) = -5.78*dum*qr3d(k)*nr3d(k)*rho(k)

         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF RAIN

      if (QR3D(K).GE.QSMALL) then
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      else
      EPSR = 0.
      end if

! NO CONDENSATION ONTO RAIN, ONLY EVAP

           if (QV3D(K).LT.QVS(K)) then
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
              ! dont allow evap to exceed saturation:  GHB, 130712:
              dum = QV3D(K)-PRE(K)*DT
              if( dum.gt.QVS(K) )then
                PRE(K)=(QV3D(K)-QVS(K))/DT
              endif
           else
              PRE(K) = 0.
           end if

!.......................................................................
! MELTING OF SNOW

! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

          if (QNI3D(K).GE.1.E-8) then

! v3 5/27/11 bug fix
!             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACS(K)

! hm fix 1/20/15
!             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(TMELT-T3D(K))/       &
!                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
!                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS10/                   &
!                   (LAMS(K)**CONS35))+DUM
             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(TMELT-T3D(K))/       &
                    XLF(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
                   (LAMS(K)**CONS35))+DUM

! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

      if (QVQVS(K).LT.1.) then
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
! bug fix V1.4
        EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)    
        EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
        PSMLT(K) = PSMLT(K)-EVPMS(K)
      end if
      end if

!.......................................................................
! MELTING OF GRAUPEL

! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

          if (QG3D(K).GE.1.E-8) then

! v3 5/27/11 bug fix
!             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACG(K)

! hm fix 1/20/15
!             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(TMELT-T3D(K))/ 		 &
!                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
!                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS11/                   &
!                   (LAMG(K)**CONS36))+DUM
             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(TMELT-T3D(K))/ 		 &
                    XLF(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
                   (LAMG(K)**CONS36))+DUM

! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

      if (QVQVS(K).LT.1.) then
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))
! bug fix V1.4
        EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
        EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
        PGMLT(K) = PGMLT(K)-EVPMG(K)
      end if
      end if

! HM, V2.1
! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

      PRACG(K) = 0.
      PRACS(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
! CALCULATION

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K))*DT

      if (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) then

        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO

        end if

! CONSERVATION OF SNOW

        DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

        if (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) then

! NO SOURCE TERMS FOR SNOW AT T > FREEZING
        RATIO = QNI3D(K)/DUM

        PSMLT(K) = PSMLT(K)*RATIO
        EVPMS(K) = EVPMS(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO

        end if

! CONSERVATION OF GRAUPEL

        DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT

        if (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) then

! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
        RATIO = QG3D(K)/DUM

        PGMLT(K) = PGMLT(K)*RATIO
        EVPMG(K) = EVPMG(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        end if

! CONSERVATION OF QR
! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

        DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT

        if (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) then

        RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                        (-PRE(K))
        PRE(K) = PRE(K)*RATIO
        
        end if

!....................................

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

!!!      T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
!!!                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)
      T3DTEN(K) = T3DTEN(K)+(PRE(K)*ecnd(k)+(EVPMS(K)+EVPMG(K))*edep(k)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*efrz(k))

      ! cm1 thbud !
      tbevar(k) = tbevar(k)+PRE(K)*ecnd(k)+(EVPMS(K)+EVPMG(K))*edep(k)
      tbmelt(k) = tbmelt(k)+(PSMLT(K)+PGMLT(K))*efrz(k)
      tbfrz(k) = tbfrz(k)+(-PRACS(K)-PRACG(K))*efrz(k)

      QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
      QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
      QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
! v3 5/27/11
!      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
!      NG3DTEN(K) = NG3DTEN(K)
      NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
      NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))

      if (PRE(K).LT.0.) then
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      end if

! V1.3 move code below to before saturation adjustment
        if (EVPMS(K)+PSMLT(K).LT.0.) then
         DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSMLTS(K) = DUM*NS3D(K)/DT
        end if
        if (PSMLT(K).LT.0.) then
          DUM = PSMLT(K)*DT/QNI3D(K)
          DUM = MAX(-1.0,DUM)
          NSMLTR(K) = DUM*NS3D(K)/DT
        end if
        if (EVPMG(K)+PGMLT(K).LT.0.) then
         DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NGMLTG(K) = DUM*NG3D(K)/DT
        end if
        if (PGMLT(K).LT.0.) then
          DUM = PGMLT(K)*DT/QG3D(K)
          DUM = MAX(-1.0,DUM)
          NGMLTR(K) = DUM*NG3D(K)/DT
        end if

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

         NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
         NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
         NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))

 300  continue

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQR = QR3D(K)+DT*QR3DTEN(K)
      DUMMQS = QNI3D(K)+DT*QNI3DTEN(K)
      DUMQG = QG3D(K)+DT*QG3DTEN(K)
      DUMQI = 0.0
      if(eqtset.eq.2) PRES(K) = RHO(K)*DUMT*(R+RV*DUMQV)
      SAVEQC = DUMQC
      SAVET = DUMT
      call satadj2(DUMT,PRES(K),DUMQV,DUMQC,DUMQR,DUMQI,DUMMQS,DUMQG,RHO(K),eqtset,tsmall)

!!!      DUMQSS = EP_2*POLYSVP(DUMT,0)/ (PRES(K)-POLYSVP(DUMT,0))
!!!      DUMQC = MAX(DUMQC,0.)
!!!
!!!! SATURATION ADJUSTMENT FOR LIQUID
!!!
!!!      DUMS = DUMQV-DUMQSS
!!!      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
!!!      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
!!!           PCC(K) = -DUMQC/DT
!!!      END IF
!!!
!!!      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
!!!      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
!!!      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
      ! GHB2

      if( abs(DUMT-SAVET).ge.0.0001 )then
        PCC(K) = (DUMQC-SAVEQC)*RDT
        QV3DTEN(K) = QV3DTEN(K)-PCC(K)
        QC3DTEN(K) = QC3DTEN(K)+PCC(K)
        T3DTEN(K) = T3DTEN(K)+(DUMT-SAVET)*RDT
        tbcond(k) = tbcond(k)+(DUMT-SAVET)*RDT
      endif

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS

    ! GHB, 120124:
    if(INUM.EQ.0)then
      if (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) then

! EFFECTIVE VERTICAL VELOCITY (M/S)

      if (ISUB.EQ.0) then
! ADD SUB-GRID VERTICAL VELOCITY
         DUM = W3D(K)+WVAR(K)

! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
         DUM = MAX(DUM,0.10)

      else if (ISUB.EQ.1) then
         DUM=W3D(K)
      end if

! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
      if (DUM.GE.0.001) then

      if (IBASE.EQ.1) then

! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

         IDROP=0

! V1.3 USE CURRENT VALUE OF QC FOR IDROP
         if (QC3D(K).LE.0.05E-3/RHO(K)) then
            IDROP=1
         end if
         if (K.EQ.1) then
            IDROP=1
         else if (K.GE.2) then
            if (QC3D(K).GT.0.05E-3/RHO(K).AND. &
             QC3D(K-1).LE.0.05E-3/RHO(K-1)) then
            IDROP=1
            end if
         end if

         if (IDROP.EQ.1) then
! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

           if (IACT.EQ.1) then
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S to CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 to M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 to KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           else if (IACT.EQ.2) then
! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2
           end if  ! IACT

!.............................................................................
        else if (IDROP.EQ.0) then
! ACTIVATE IN CLOUD INTERIOR
! FIND EQUILIBRIUM SUPERSATURATION

           TAUC=1./(2.*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.)/LAMC(K))
           if (EPSR.GT.1.E-8) then
             TAUR=1./EPSR
           else
             TAUR=1.E8
           end if

! hm fix 1/20/15
!           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=(-QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=DUM3*TAUC*TAUR/(TAUC+TAUR)

           if (DUM3/QVS(K).GE.1.E-6) then
           if (IACT.EQ.1) then

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

            DUM=DUM*100.  ! CONVERT FROM M/S to CM/S
            DUMACT = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))

! USE POWER LAW CCN SPECTRA

! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
            DUM3=DUM3/QVS(K)*100.

            DUM2=C1*DUM3**K1
! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
            DUM2=MIN(DUM2,DUMACT)
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 to M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 to KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           else if (IACT.EQ.2) then

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

           DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! USE LOGNORMAL AEROSOL
           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
           SMAX = DUM3/QVS(K)

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
            DUM2=MIN(DUM2,DUMACT)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           end if ! IACT
           end if ! DUM3/QVS > 1.E-6
        end if  ! IDROP = 1

!.......................................................................
      else if (IBASE.EQ.2) then

           if (IACT.EQ.1) then
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S to CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 to M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 to KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           else if (IACT.EQ.2) then

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2
           end if  ! IACT
        end if  ! IBASE
        end if  ! W > 0.001
        end if  ! QC3D > QSMALL
    ! GHB, 120124:
    endif     ! INUM=0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

!.....................................................................
!.....................................................................
         else  ! TEMPERATURE < 273.15

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         if (INUM.EQ.1) then
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         end if

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      if(INUM.EQ.0)  &
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

      if (QI3D(K).GE.QSMALL) then
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)
         N0I(K) = NI3D(K)*LAMI(K)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMI(K).LT.LAMMINI) then

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      else if (LAMI(K).GT.LAMMAXI) then
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      end if
      end if

!......................................................................
! RAIN

      if (QR3D(K).GE.QSMALL) then
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMR(K).LT.LAMMINR) then

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      else if (LAMR(K).GT.LAMMAXR) then
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      end if
      end if

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      if (QC3D(K).GE.QSMALL) then

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! v1.4
! interpolate
         dumii=int(pgam(k))
         nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam(k)-real(dumii))

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      if (LAMC(K).LT.LAMMIN) then
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      else if (LAMC(K).GT.LAMMAX) then
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      end if

! TO CALCULATE DROPLET FREEZING

        CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)

      end if

!......................................................................
! SNOW

      if (QNI3D(K).GE.QSMALL) then
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMS(K).LT.LAMMINS) then
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      else if (LAMS(K).GT.LAMMAXS) then

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      end if
      end if

!......................................................................
! GRAUPEL

      if (QG3D(K).GE.QSMALL) then
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMG(K).LT.LAMMING) then
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      else if (LAMG(K).GT.LAMMAXG) then

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      end if
      end if

!.....................................................................
! ZERO OUT PROCESS RATES

            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PGSACW(K) = 0.
            PGRACS(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES
! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
!.......................................................................
! FREEZING OF CLOUD DROPLETS
! ONLY ALLOWED BELOW -4 C
        if (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15) then

! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

! MEYERS CURVE

           NACNT = EXP(-2.80+0.262*(TMELT-T3D(K)))*1000.

! COOPER CURVE
!        NACNT =  5.*EXP(0.304*(TMELT-T3D(K)))

! FLECTHER
!     NACNT = 0.01*EXP(0.6*(TMELT-T3D(K)))

! CONTACT FREEZING

! MEAN FREE PATH

            DUM = 7.37*T3D(K)/(288.*10.*PRES(K))/100.

! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
! BASED ON BROWNIAN DIFFUSION

            DAP(K) = CONS37*T3D(K)*(1.+DUM/RIN)/MU(K)
 
           MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                   LOG(GAMMA(PGAM(K)+5.))-4.*LOG(LAMC(K)))
           NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
                    GAMMA(PGAM(K)+2.)/                         &
                    LAMC(K)

! IMMERSION FREEZING (BIGG 1953)

!           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
!                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
!                   EXP(AIMM*(273.15-T3D(K)))

!           NNUCCC(K) = NNUCCC(K)+                                  &
!            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
!                *EXP(AIMM*(273.15-T3D(K)))

! hm 7/15/13 fix for consistency w/ original formula
           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                   (EXP(AIMM*(273.15-T3D(K)))-1.)

           NNUCCC(K) = NNUCCC(K)+                                  &
            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                *(EXP(AIMM*(273.15-T3D(K)))-1.)

! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
! MIXING RATIO, SINCE CONSERVATION NOT APPLIED TO NUMBER CONC

           NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

        end if

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         if (QC3D(K).GE.1.E-6) then

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

            if (IRAIN.EQ.0) then

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

             else if (IRAIN.EQ.1) then

! v1.4
! replace with seifert and beheng

        dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
        dum1 = 600.*dum**0.68*(1.-dum**0.68)**3

        prc(k) = 9.44e9/(20.*2.6e-7)* &
        (nu(k)+2.)*(nu(k)+4.)/(nu(k)+1.)**2* &
        (rho(k)*qc3d(k)/1000.)**4/(rho(k)*nc3d(k)/1.e6)**2* &
        (1.+dum1/(1.-dum)**2)*1000./rho(k)

        nprc(k) = prc(k)*2./2.6e-7*1000.
        nprc1(k) = 0.5*nprc(k)

         end if
         end if

!.......................................................................
! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

         if (QNI3D(K).GE.1.E-8) then
             NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
            ((2.+BS)/3.)*QNI3D(K)**((2.+BS)/3.)*                  &
            (NS3D(K)*RHO(K))**((4.-BS)/3.)/                       &
            (RHO(K))
         end if

!.......................................................................
! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
! HERE USE CONTINUOUS COLLECTION EQUATION WITH
! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

! SNOW

         if (QNI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) then

           PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)
           NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)

         end if

!............................................................................
! COLLECTION OF CLOUD WATER BY GRAUPEL

         if (QG3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) then

           PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
           NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
	    end if

!.......................................................................
! HM, ADD 12/13/06
! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
! BEFORE RIMING CAN OCCUR
! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
! TO HALLET-MOSSOP SPLINTERING

         if (QI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) then

! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
! FROM THOMPSON ET AL. 2004, MWR

            if (1./LAMI(K).GE.100.E-6) then

           PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           end if
         end if

!.......................................................................
! ACCRETION OF RAIN WATER BY SNOW
! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

         if (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) then

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                (1./(LAMR(K)**3*LAMS(K))+                      &
                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
                 1./(LAMR(K)*LAMS(K)**3))

! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

! V1.3
! V1.5
!            IF (IHAIL.EQ.0) THEN
            if (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) then
            PSACR(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                 N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                  (5./(LAMS(K)**3*LAMR(K))+                    &
                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                  0.5/(LAMS(K)*LAMR(K)**3)))            
            end if
!            END IF

         end if

!.......................................................................

! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990, 
! USED BY REISNER ET AL 1998
         if (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) then

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

	    end if

!.......................................................................
! RIME-SPLINTERING - SNOW
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

!v1.4
         if (QNI3D(K).GE.0.1E-3) then
         if (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) then
         if (PSACWS(K).GT.0..OR.PRACS(K).GT.0.) then
            if (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) then

               if (T3D(K).GT.270.16) then
                  FMULT = 0.
               else if (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  then
                  FMULT = (270.16-T3D(K))/2.
               else if (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   then
                  FMULT = (T3D(K)-265.16)/3.
               else if (T3D(K).LT.265.16) then
                  FMULT = 0.
               end if

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

               if (PSACWS(K).GT.0.) then
                  NMULTS(K) = 35.E4*PSACWS(K)*FMULT*1000.
                  QMULTS(K) = NMULTS(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                  QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                  PSACWS(K) = PSACWS(K)-QMULTS(K)

               end if

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               if (PRACS(K).GT.0.) then
                   NMULTR(K) = 35.E4*PRACS(K)*FMULT*1000.
                   QMULTR(K) = NMULTR(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                   QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                   PRACS(K) = PRACS(K)-QMULTR(K)

               end if

            end if
            end if
         end if
         end if

!.......................................................................
! RIME-SPLINTERING - GRAUPEL 
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

! V1.3
! V1.5
!         IF (IHAIL.EQ.0) THEN
! v1.4
         if (QG3D(K).GE.0.1E-3) then
         if (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) then
         if (PSACWG(K).GT.0..OR.PRACG(K).GT.0.) then
            if (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) then

               if (T3D(K).GT.270.16) then
                  FMULT = 0.
               else if (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  then
                  FMULT = (270.16-T3D(K))/2.
               else if (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   then
                  FMULT = (T3D(K)-265.16)/3.
               else if (T3D(K).LT.265.16) then
                  FMULT = 0.
               end if

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

               if (PSACWG(K).GT.0.) then
                  NMULTG(K) = 35.E4*PSACWG(K)*FMULT*1000.
                  QMULTG(K) = NMULTG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                  QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                  PSACWG(K) = PSACWG(K)-QMULTG(K)

               end if

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               if (PRACG(K).GT.0.) then
                   NMULTRG(K) = 35.E4*PRACG(K)*FMULT*1000.
                   QMULTRG(K) = NMULTRG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                   QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                   PRACG(K) = PRACG(K)-QMULTRG(K)

               end if

            end if
            end if
            end if
            end if
!         END IF

!........................................................................
! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL

! V1.3
! V1.5
!           IF (IHAIL.EQ.0) THEN
	   if (PSACWS(K).GT.0.) then
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              if (QNI3D(K).GE.0.1E-3.AND.QC3D(K).GE.0.5E-3) then

! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
	     PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
                          ASN(K)*ASN(K)/ &
                           (RHO(K)*LAMS(K)**(2.*BS+2.))) 

! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
	     DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.) 

! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
	     NSCNG(K) = DUM/MG0*RHO(K)
! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
             NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

! PORTION OF RIMING LEFT FOR SNOW
             PSACWS(K) = PSACWS(K) - PGSACW(K)
             end if
	   end if

! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

	   if (PRACS(K).GT.0.) then
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              if (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) then
! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
	      DUM = CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3 &    
                   /(CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3+ &  
                   CONS19*(4./LAMR(K))**3*(4./LAMR(K))**3)
              DUM=MIN(DUM,1.)
              DUM=MAX(DUM,0.)
	      PGRACS(K) = (1.-DUM)*PRACS(K)
            NGRACS(K) = (1.-DUM)*NPRACS(K)
! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
            NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
            NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

! AMOUNT LEFT FOR SNOW PRODUCTION
            PRACS(K) = PRACS(K) - PGRACS(K)
            NPRACS(K) = NPRACS(K) - NGRACS(K)
! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
            PSACR(K)=PSACR(K)*(1.-DUM)
            end if
	   end if
!           END IF

!.......................................................................
! FREEZING OF RAIN DROPS
! FREEZING ALLOWED BELOW -4 C

         if (T3D(K).LT.269.15.AND.QR3D(K).GE.QSMALL) then

! IMMERSION FREEZING (BIGG 1953)
!            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3 &
!                 /LAMR(K)**3

!            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3

! hm fix 7/15/13 for consistency w/ original formula
            MNUCCR(K) = CONS20*NR3D(K)*((EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3 &
                 /LAMR(K)**3)

            NNUCCR(K) = PI*NR3D(K)*BIMM*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3

! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
            NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

         end if

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         if (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) then

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

            if (IRAIN.EQ.0) then

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

           else if (IRAIN.EQ.1) then

! v1.4
! seifert and beheng (2001) formulation

           dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
           dum1 = (dum/(dum+5.e-4))**4
           pra(k) = 5.78e3*rho(k)/1000.*qc3d(k)*qr3d(k)*dum1
           npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
           (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)

         end if
         end if
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

! v1.4 replace with seifert and beheng (2001)

         if (QR3D(K).GE.1.E-8) then
! v1.4
! seifert and beheng
! include breakup, V2.1
            dum1=D_b
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
            nragg(k) = -5.78*dum*qr3d(k)*nr3d(k)*rho(k)

         end if

!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

         if (QI3D(K).GE.1.E-8 .AND.QVQVSI(K).GE.1.) then

!           COFFI = 2./LAMI(K)
!           IF (COFFI.GE.DCS) THEN
              NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
              PRCI(K) = CONS22*NPRCI(K)
              NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

!           END IF
         end if

!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

         if (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) then
            PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                     LAMS(K)**(BS+3.)
            NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                  RHO(K)*N0S(K)/                                 &
                  LAMS(K)**(BS+3.)
            NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
         end if

!.......................................................................
! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
! FOLLOWS REISNER ET AL. 1998
! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

         if (QR3D(K).GE.1.E-8.AND.QI3D(K).GE.1.E-8.AND.T3D(K).LE.TMELT) then

! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
! OTHERWISE ADD TO SNOW

            if (QR3D(K).GE.0.1E-3) then
            NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
            NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
            else 
            NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
            NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
            end if
         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

         if (INUC.EQ.0) then

! FREEZING OF AEROSOL ONLY ALLOWED BELOW -5 C
! AND ABOVE DELIQUESCENCE THRESHOLD OF 80%
! AND ABOVE ICE SATURATION

! add threshold according to Greg Thomspon

         if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
              QVQVSI(K).ge.1.25) then

! hm, modify dec. 5, 2006, replace with cooper curve
      kc2 = 0.005*exp(0.304*(TMELT-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
      kc2 = min(kc2,500.e3)
      kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1

          if (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) then
             NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD(K) = NNUCCD(K)*MI0
          end if

          end if

          else if (INUC.EQ.1) then

          if (T3D(K).LT.TMELT.AND.QVQVSI(K).GT.1.) then

             KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 to KG-1
          if (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) then
             NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD(K) = NNUCCD(K)*MI0
          end if
          end if

         end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 101      continue

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

! NO VENTILATION FOR CLOUD ICE

        if (QI3D(K).GE.QSMALL) then

         EPSI = 2.*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

      else
         EPSI = 0.
      end if

      if (QNI3D(K).GE.QSMALL) then
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
      else
      EPSS = 0.
      end if

      if (QG3D(K).GE.QSMALL) then
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))


      else
      EPSG = 0.
      end if

      if (QR3D(K).GE.QSMALL) then
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      else
      EPSR = 0.
      end if

! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

              if (QI3D(K).GE.QSMALL) then              
              DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
              PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
              else
              DUM=0.
              end if
! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
              if (QNI3D(K).GE.QSMALL) then
              PRDS(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
                EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
! OTHERWISE ADD TO CLOUD ICE
              else
              PRD(K) = PRD(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
              end if
! VAPOR DPEOSITION ON GRAUPEL
              PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)

! NO CONDENSATION ONTO RAIN, ONLY EVAP

           if (QV3D(K).LT.QVS(K)) then
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
              ! dont allow evap to exceed saturation:  GHB, 130712:
              dum = QV3D(K)-PRE(K)*DT
              if( dum.gt.QVS(K) )then
                PRE(K)=(QV3D(K)-QVS(K))/DT
              endif
           else
              PRE(K) = 0.
           end if

! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
! FORMULA FROM REISNER 2 SCHEME

           DUM = (QV3D(K)-QVI(K))/DT

           FUDGEF = 0.9999
           SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)

           if( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
               (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) then
               MNUCCD(K) = FUDGEF*MNUCCD(K)*DUM/SUM_DEP
               PRD(K) = FUDGEF*PRD(K)*DUM/SUM_DEP
               PRDS(K) = FUDGEF*PRDS(K)*DUM/SUM_DEP
	       PRDG(K) = FUDGEF*PRDG(K)*DUM/SUM_DEP
           endif

! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

           if (PRD(K).LT.0.) then
              EPRD(K)=PRD(K)
              PRD(K)=0.
           end if
           if (PRDS(K).LT.0.) then
              EPRDS(K)=PRDS(K)
              PRDS(K)=0.
           end if
           if (PRDG(K).LT.0.) then
              EPRDG(K)=PRDG(K)
              PRDG(K)=0.
           end if

!.......................................................................
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CONSERVATION OF WATER
! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.
! THIS SECTION IS SEPARATED INTO TWO PARTS, IF T < 0 C, T > 0 C
! DUE TO DIFFERENT PROCESSES THAT ACT DEPENDING ON FREEZING/ABOVE FREEZING

! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

! NOTE: CONSERVATION NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
! STEP

!****SENSITIVITY - NO ICE

      if (ILIQ.EQ.1) then
      MNUCCC(K)=0.
      NNUCCC(K)=0.
      MNUCCR(K)=0.
      NNUCCR(K)=0.
      MNUCCD(K)=0.
      NNUCCD(K)=0.
      end if

! ****SENSITIVITY - NO GRAUPEL
      if (IGRAUP.EQ.1) then
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PGSACW(K) = 0.
            PGRACS(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
            EVPMG(K) = 0.
            PGMLT(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.
	    NGMLTG(K) = 0.
            NGMLTR(K) = 0.
! v3 5/27/11
            PIACRS(K)=PIACRS(K)+PIACR(K)
            PIACR(K) = 0.
! fix 070713
	    PRACIS(K)=PRACIS(K)+PRACI(K)
	    PRACI(K) = 0.
	    PSACWS(K)=PSACWS(K)+PGSACW(K)
	    PGSACW(K) = 0.
	    PRACS(K)=PRACS(K)+PGRACS(K)
	    PGRACS(K) = 0.
       end if

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT

      if (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) then
        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO
        MNUCCC(K) = MNUCCC(K)*RATIO
        PSACWS(K) = PSACWS(K)*RATIO
        PSACWI(K) = PSACWI(K)*RATIO
        QMULTS(K) = QMULTS(K)*RATIO
        QMULTG(K) = QMULTG(K)*RATIO
        PSACWG(K) = PSACWG(K)*RATIO
	PGSACW(K) = PGSACW(K)*RATIO
        end if
 
! CONSERVATION OF QI

      DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
                -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT

      if (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) then

        RATIO = (QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
                     MNUCCD(K)+PSACWI(K))/ &
                      (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))

        PRCI(K) = PRCI(K)*RATIO
        PRAI(K) = PRAI(K)*RATIO
        PRACI(K) = PRACI(K)*RATIO
        PRACIS(K) = PRACIS(K)*RATIO
        EPRD(K) = EPRD(K)*RATIO

        end if

! CONSERVATION OF QR

      DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
             PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT

      if (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) then

        RATIO = (QR3D(K)/DT+PRC(K)+PRA(K))/ &
             (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))

        PRE(K) = PRE(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO
        QMULTR(K) = QMULTR(K)*RATIO
        QMULTRG(K) = QMULTRG(K)*RATIO
        MNUCCR(K) = MNUCCR(K)*RATIO
        PIACR(K) = PIACR(K)*RATIO
        PIACRS(K) = PIACRS(K)*RATIO
        PGRACS(K) = PGRACS(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        end if

! CONSERVATION OF QNI
! CONSERVATION FOR GRAUPEL SCHEME

        if (IGRAUP.EQ.0) then

      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT

      if (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) then

        RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K))

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       end if

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       else if (IGRAUP.EQ.1) then

      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT

      if (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) then

       RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K))

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       end if

       end if

! CONSERVATION OF QG

      DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT

      if (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) then

        RATIO = (QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
                  PIACR(K)+PRACI(K))/(-EPRDG(K))

       EPRDG(K) = EPRDG(K)*RATIO

      end if

! TENDENCIES

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))

! v3 5/27/11 bug fix
!!!      T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
!!!               *XXLV(K)+(PRD(K)+PRDS(K)+                            &
!!!                MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
!!!               (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
!!!                QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
!!!                +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K))*XLF(K))/CPM(K)
      T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
               *ecnd(k)+(PRD(K)+PRDS(K)+                            &
                MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*edep(k)+         &
               (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
                QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
                +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*efrz(k))

      ! cm1 thbud !
      tbevar(k) = tbevar(k)+PRE(K)*ecnd(k)
      tbfrz(k) = tbfrz(k)+( MNUCCC(K)+MNUCCR(K)  &
                           +PSACWS(K)+PSACWI(K)+PSACWG(K)  &
                           +QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)  &
                           +PRACS(K)+PRACG(K)+PGSACW(K)+PGRACS(K)  &
                           +PIACR(K)+PIACRS(K)  &
                          )*efrz(k)
      tbdep(k) = tbdep(k)+( PRD(K)+PRDS(K)+PRDG(K)+MNUCCD(K) )*edep(k)
      tbsub(k) = tbsub(k)+( EPRD(K)+EPRDS(K)+EPRDG(K) )*edep(k)


      QC3DTEN(K) = QC3DTEN(K)+                                      &
                 (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
                  PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
      QI3DTEN(K) = QI3DTEN(K)+                                      &
         (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
                  PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
      QR3DTEN(K) = QR3DTEN(K)+                                      &
                 (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
             -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))

      if (IGRAUP.EQ.0) then

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
                    PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
      NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
      else if (IGRAUP.EQ.1) then

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))

      end if

      NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
            -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))

      NI3DTEN(K) = NI3DTEN(K)+                                      &
       (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
               NNUCCD(K)-NIACR(K)-NIACRS(K))

      NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
                   +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))

! V1.3 move code below to before saturation adjustment
      if (EPRD(K).LT.0.) then
         DUM = EPRD(K)*DT/QI3D(K)
            DUM = MAX(-1.,DUM)
         NSUBI(K) = DUM*NI3D(K)/DT
      end if
      if (EPRDS(K).LT.0.) then
         DUM = EPRDS(K)*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSUBS(K) = DUM*NS3D(K)/DT
      end if
      if (PRE(K).LT.0.) then
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      end if
      if (EPRDG(K).LT.0.) then
         DUM = EPRDG(K)*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NSUBG(K) = DUM*NG3D(K)/DT
      end if

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

         NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
         NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
         NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
         NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQR = QR3D(K)+DT*QR3DTEN(K)
      DUMMQS = QNI3D(K)+DT*QNI3DTEN(K)
      DUMQG = QG3D(K)+DT*QG3DTEN(K)
      DUMQI = QI3D(K)+DT*QI3DTEN(K)
      if(eqtset.eq.2) PRES(K) = RHO(K)*DUMT*(R+RV*DUMQV)
      SAVEQC = DUMQC
      SAVET = DUMT
      call satadj2(DUMT,PRES(K),DUMQV,DUMQC,DUMQR,DUMQI,DUMMQS,DUMQG,RHO(K),eqtset,tsmall)

!!!      DUMQSS = EP_2*POLYSVP(DUMT,0)/ (PRES(K)-POLYSVP(DUMT,0))
!!!      DUMQC = MAX(DUMQC,0.)
!!!
!!!! SATURATION ADJUSTMENT FOR LIQUID
!!!
!!!      DUMS = DUMQV-DUMQSS
!!!      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
!!!      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
!!!           PCC(K) = -DUMQC/DT
!!!      END IF
!!!
!!!      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
!!!      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
!!!      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
      ! GHB2

      if( abs(DUMT-SAVET).ge.0.0001 )then
        PCC(K) = (DUMQC-SAVEQC)*RDT
        QV3DTEN(K) = QV3DTEN(K)-PCC(K)
        QC3DTEN(K) = QC3DTEN(K)+PCC(K)
        T3DTEN(K) = T3DTEN(K)+(DUMT-SAVET)*RDT
        tbcond(k) = tbcond(k)+(DUMT-SAVET)*RDT
      endif

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS

    ! GHB, 120124:
    if(INUM.EQ.0)then
      if (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) then

! EFFECTIVE VERTICAL VELOCITY (M/S)

      if (ISUB.EQ.0) then
! ADD SUB-GRID VERTICAL VELOCITY
         DUM = W3D(K)+WVAR(K)

! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
         DUM = MAX(DUM,0.10)

      else if (ISUB.EQ.1) then
         DUM=W3D(K)
      end if

! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
      if (DUM.GE.0.001) then

      if (IBASE.EQ.1) then

! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

         IDROP=0

! V1.3 USE CURRENT VALUE OF QC FOR IDROP
         if (QC3D(K).LE.0.05E-3/RHO(K)) then
            IDROP=1
         end if
         if (K.EQ.1) then
            IDROP=1
         else if (K.GE.2) then
            if (QC3D(K).GT.0.05E-3/RHO(K).AND. &
             QC3D(K-1).LE.0.05E-3/RHO(K-1)) then
            IDROP=1
            end if
         end if

         if (IDROP.EQ.1) then
! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

           if (IACT.EQ.1) then
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S to CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 to M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 to KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           else if (IACT.EQ.2) then
! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2
           end if  ! IACT

!.............................................................................
        else if (IDROP.EQ.0) then
! ACTIVATE IN CLOUD INTERIOR
! FIND EQUILIBRIUM SUPERSATURATION

           TAUC=1./(2.*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.)/LAMC(K))
           if (EPSR.GT.1.E-8) then
             TAUR=1./EPSR
           else
             TAUR=1.E8
           end if
           if (EPSI.GT.1.E-8) then
             TAUI=1./EPSI
           else
             TAUI=1.E8
           end if
           if (EPSS.GT.1.E-8) then
             TAUS=1./EPSS
           else
             TAUS=1.E8
           end if
           if (EPSG.GT.1.E-8) then
             TAUG=1./EPSG
           else
             TAUG=1.E8
           end if

! EQUILIBRIUM SS INCLUDING BERGERON EFFECT

! hm fix 1/20/15
!           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=(-QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=(DUM3*TAUC*TAUR*TAUI*TAUS*TAUG- &
           (QVS(K)-QVI(K))*(TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS))/ &
           (TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS+ &
            TAUR*TAUI*TAUS*TAUG+TAUC*TAUI*TAUS*TAUG)

           if (DUM3/QVS(K).GE.1.E-6) then
           if (IACT.EQ.1) then

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

            DUM=DUM*100.  ! CONVERT FROM M/S to CM/S
            DUMACT = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))

! USE POWER LAW CCN SPECTRA

! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
            DUM3=DUM3/QVS(K)*100.

            DUM2=C1*DUM3**K1
! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
            DUM2=MIN(DUM2,DUMACT)
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 to M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 to KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           else if (IACT.EQ.2) then

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

           DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! USE LOGNORMAL AEROSOL
           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
           SMAX = DUM3/QVS(K)

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
            DUM2=MIN(DUM2,DUMACT)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           end if ! IACT
           end if ! DUM3/QVS > 1.E-6
        end if  ! IDROP = 1

!.......................................................................
      else if (IBASE.EQ.2) then

           if (IACT.EQ.1) then
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S to CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 to M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 to KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2

           else if (IACT.EQ.2) then

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT to KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
            NC3DTEN(K) = NC3DTEN(K)+DUM2
           end if  ! IACT
        end if  ! IBASE
        end if  ! W > 0.001
        end if  ! QC3D > QSMALL
    ! GHB, 120124:
    endif       ! INUM=0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

         end if !!!!!! TEMPERATURE

! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
         LTRUE = 1

 200     continue

         enddo  kloop1

! V1.3 move precip initialization to here
! INITIALIZE PRECIP AND SNOW RATES

      PRECRT = 0.
      SNOWRT = 0.

! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

        if (LTRUE.EQ.0) GOTO 400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!.......................................................................
! CALCULATE SEDIMENATION
! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
! STABILITY, I.E. COURANT# < 1

!.......................................................................

      NSTEPR=1
      NSTEPI=1
      NSTEPS=1
      NSTEPC=1
      NSTEPG=1

      do K = KTS,KTE
        temz(k) = 1.0/DZQ(K)
      enddo

! v3 5/27/11
      kloop2:  &
      do K = KTE,KTS,-1

        DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
        DUMQS(K) = QNI3D(K)+QNI3DTEN(K)*DT
        DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
        DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
        DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
        DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
        DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
        DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
	DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
	DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT

! SWITCH FOR CONSTANT DROPLET NUMBER
        if (INUM.EQ.1) then
        DUMFNC(K) = NC3D(K)
        end if

! GET DUMMY LAMDA

! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
      DUMFNI(K) = MAX(0.,DUMFNI(K))
      DUMFNS(K) = MAX(0.,DUMFNS(K))
      DUMFNC(K) = MAX(0.,DUMFNC(K))
      DUMFNR(K) = MAX(0.,DUMFNR(K))
      DUMFNG(K) = MAX(0.,DUMFNG(K))

!......................................................................
! CLOUD ICE

      if (DUMI(K).GE.QSMALL) then
        DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
        DLAMI=MAX(DLAMI,LAMMINI)
        DLAMI=MIN(DLAMI,LAMMAXI)
      end if
!......................................................................
! RAIN

      if (DUMR(K).GE.QSMALL) then
        DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1./3.)
        DLAMR=MAX(DLAMR,LAMMINR)
        DLAMR=MIN(DLAMR,LAMMAXR)
      end if
!......................................................................
! CLOUD DROPLETS

      if (DUMC(K).GE.QSMALL) then
!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

        DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.)/(DUMC(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
        LAMMIN = (PGAM(K)+1.)/60.E-6
        LAMMAX = (PGAM(K)+1.)/1.E-6
        DLAMC=MAX(DLAMC,LAMMIN)
        DLAMC=MIN(DLAMC,LAMMAX)
      end if
!......................................................................
! SNOW

      if (DUMQS(K).GE.QSMALL) then
        DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1./DS)
        DLAMS=MAX(DLAMS,LAMMINS)
        DLAMS=MIN(DLAMS,LAMMAXS)
      end if
!......................................................................
! GRAUPEL

      if (DUMG(K).GE.QSMALL) then
        DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1./DG)
        DLAMG=MAX(DLAMG,LAMMING)
        DLAMG=MIN(DLAMG,LAMMAXG)
      end if

!......................................................................
! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

! CLOUD WATER

      if (DUMC(K).GE.QSMALL) then
      UNC =  ACN(K)*GAMMA(1.+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.))
      UMC = ACN(K)*GAMMA(4.+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.))
      else
      UMC = 0.
      UNC = 0.
      end if

      if( testcase.eq.5 )then
        ! drizzling stratoCu (Ackerman et al, 2009, MWR)
        if( dumc(k).ge.qsmall )then
          umc = 1.19e8                                       &
               *(3.0/(4.0*pi*rhow*ndcnst*1.0e6))**(2.0/3.0)  &
               *(rho(k)*dumc(k))**(2.0/3.0)                  &
               *exp(5.0*log(1.5)**2)
        else
          umc = 0.0
        endif
      endif

      if (DUMI(K).GE.QSMALL) then
      UNI =  AIN(K)*CONS27/DLAMI**BI
      UMI = AIN(K)*CONS28/(DLAMI**BI)
      else
      UMI = 0.
      UNI = 0.
      end if

      if (DUMR(K).GE.QSMALL) then
      UNR = ARN(K)*CONS6/DLAMR**BR
      UMR = ARN(K)*CONS4/(DLAMR**BR)
      else
      UMR = 0.
      UNR = 0.
      end if

      if (DUMQS(K).GE.QSMALL) then
      UMS = ASN(K)*CONS3/(DLAMS**BS)
      UNS = ASN(K)*CONS5/DLAMS**BS
      else
      UMS = 0.
      UNS = 0.
      end if

      if (DUMG(K).GE.QSMALL) then
      UMG = AGN(K)*CONS7/(DLAMG**BG)
      UNG = AGN(K)*CONS8/DLAMG**BG
      else
      UMG = 0.
      UNG = 0.
      end if

! SET REALISTIC LIMITS ON FALLSPEED

! bug fix, 10/08/09
        dum=(rhosu/rho(k))**0.54
        UMS=MIN(UMS,1.2*dum)
        UNS=MIN(UNS,1.2*dum)
! v3 5/27/11
! fix for correction by AA 4/6/11
        UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)
        UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)
        UMR=MIN(UMR,9.1*dum)
        UNR=MIN(UNR,9.1*dum)
        UMG=MIN(UMG,20.*dum)
        UNG=MIN(UNG,20.*dum)

      FR(K) = UMR
      FI(K) = UMI
      FNI(K) = UNI
      FS(K) = UMS
      FNS(K) = UNS
      FNR(K) = UNR
      FC(K) = UMC
      FNC(K) = UNC
      FG(K) = UMG
      FNG(K) = UNG

! v3 5/27/11 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

	if (K.LE.KTE-1) then
        if (FR(K).LT.1.E-10) then
	FR(K)=FR(K+1)
	end if
        if (FI(K).LT.1.E-10) then
	FI(K)=FI(K+1)
	end if
        if (FNI(K).LT.1.E-10) then
	FNI(K)=FNI(K+1)
	end if
        if (FS(K).LT.1.E-10) then
	FS(K)=FS(K+1)
	end if
        if (FNS(K).LT.1.E-10) then
	FNS(K)=FNS(K+1)
	end if
        if (FNR(K).LT.1.E-10) then
	FNR(K)=FNR(K+1)
	end if
        if (FC(K).LT.1.E-10) then
	FC(K)=FC(K+1)
	end if
        if (FNC(K).LT.1.E-10) then
	FNC(K)=FNC(K+1)
	end if
        if (FG(K).LT.1.E-10) then
	FG(K)=FG(K+1)
	end if
        if (FNG(K).LT.1.E-10) then
	FNG(K)=FNG(K+1)
	end if
	end if ! K LE KTE-1

      RGVM = MAX(FR(K),FNR(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPR = MAX(INT(RGVM*DT*temz(K)+1.),NSTEPR)

      RGVM = MAX(FI(K),FNI(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPI = MAX(INT(RGVM*DT*temz(K)+1.),NSTEPI)

      RGVM = MAX(FS(K),FNS(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPS = MAX(INT(RGVM*DT*temz(K)+1.),NSTEPS)

      RGVM = MAX(FC(K),FNC(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPC = MAX(INT(RGVM*DT*temz(K)+1.),NSTEPC)

      RGVM = MAX(FG(K),FNG(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEPG = MAX(INT(RGVM*DT*temz(K)+1.),NSTEPG)

! MULTIPLY VARIABLES BY RHO
      DUMR(k) = DUMR(k)*RHO(K)
      DUMI(k) = DUMI(k)*RHO(K)
      DUMFNI(k) = DUMFNI(K)*RHO(K)
      DUMQS(k) = DUMQS(K)*RHO(K)
      DUMFNS(k) = DUMFNS(K)*RHO(K)
      DUMFNR(k) = DUMFNR(K)*RHO(K)
      DUMC(k) = DUMC(K)*RHO(K)
      DUMFNC(k) = DUMFNC(K)*RHO(K)
      DUMG(k) = DUMG(K)*RHO(K)
      DUMFNG(k) = DUMFNG(K)*RHO(K)

      enddo  kloop2

! V1.3, change so that sub-stepping is done
! individually for each species

! RAIN

      dttmp = DT/NSTEPR
      do K = KTS,KTE
        temr(k) = 1.0/(NSTEPR*RHO(k))
      enddo

      do N = 1,NSTEPR

      do K = KTS,KTE
      FALOUTR(K) = FR(K)*DUMR(K)
      FALOUTNR(K) = FNR(K)*DUMFNR(K)
      end do
! TOP OF MODEL
      K = KTE
      FALTNDR = FALOUTR(K)*temz(k)
      FALTNDNR = FALOUTNR(K)*temz(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
      QRSTEN(K) = QRSTEN(K)-FALTNDR*temr(k)
      NR3DTEN(K) = NR3DTEN(K)-FALTNDNR*temr(k)
      DUMR(K) = DUMR(K)-FALTNDR*dttmp
      DUMFNR(K) = DUMFNR(K)-FALTNDNR*dttmp
      do K = KTE-1,KTS,-1
      FALTNDR = (FALOUTR(K+1)-FALOUTR(K))*temz(K)
      FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))*temz(K)
      QRSTEN(K) = QRSTEN(K)+FALTNDR*temr(k)
      NR3DTEN(K) = NR3DTEN(K)+FALTNDNR*temr(k)
      DUMR(K) = DUMR(K)+FALTNDR*dttmp
      DUMFNR(K) = DUMFNR(K)+FALTNDNR*dttmp
      end do
      PRECRT = PRECRT+(FALOUTR(KTS))  &
                     *dttmp
      end do

! CLOUD ICE

      dttmp = DT/NSTEPI
      do K = KTS,KTE
        temr(k) = 1.0/(NSTEPI*RHO(k))
      enddo

      do N = 1,NSTEPI

      do K = KTS,KTE
      FALOUTI(K) = FI(K)*DUMI(K)
      FALOUTNI(K) = FNI(K)*DUMFNI(K)
      end do
! TOP OF MODEL
      K = KTE
      FALTNDI = FALOUTI(K)*temz(k)
      FALTNDNI = FALOUTNI(K)*temz(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
      QISTEN(K) = QISTEN(K)-FALTNDI*temr(k)
      NI3DTEN(K) = NI3DTEN(K)-FALTNDNI*temr(k)
      DUMI(K) = DUMI(K)-FALTNDI*dttmp
      DUMFNI(K) = DUMFNI(K)-FALTNDNI*dttmp
      do K = KTE-1,KTS,-1
      FALTNDI = (FALOUTI(K+1)-FALOUTI(K))*temz(K)
      FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))*temz(K)
      QISTEN(K) = QISTEN(K)+FALTNDI*temr(k)
      NI3DTEN(K) = NI3DTEN(K)+FALTNDNI*temr(k)
      DUMI(K) = DUMI(K)+FALTNDI*dttmp
      DUMFNI(K) = DUMFNI(K)+FALTNDNI*dttmp
      end do
      PRECRT = PRECRT+(FALOUTI(KTS))  &
                     *dttmp
      SNOWRT = SNOWRT+(FALOUTI(KTS))*dttmp
      end do

! CLOUD DROPLETS

      dttmp = DT/NSTEPC
      do K = KTS,KTE
        temr(k) = 1.0/(NSTEPC*RHO(k))
      enddo

      do N = 1,NSTEPC

      do K = KTS,KTE
      FALOUTC(K) = FC(K)*DUMC(K)
      FALOUTNC(K) = FNC(K)*DUMFNC(K)
      end do
! TOP OF MODEL
      K = KTE
      FALTNDC = FALOUTC(K)*temz(k)
      FALTNDNC = FALOUTNC(K)*temz(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
      QCSTEN(K) = QCSTEN(K)-FALTNDC*temr(k)
      NC3DTEN(K) = NC3DTEN(K)-FALTNDNC*temr(k)
      DUMC(K) = DUMC(K)-FALTNDC*dttmp
      DUMFNC(K) = DUMFNC(K)-FALTNDNC*dttmp
      do K = KTE-1,KTS,-1
      FALTNDC = (FALOUTC(K+1)-FALOUTC(K))*temz(K)
      FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))*temz(K)
      QCSTEN(K) = QCSTEN(K)+FALTNDC*temr(k)
      NC3DTEN(K) = NC3DTEN(K)+FALTNDNC*temr(k)
      DUMC(K) = DUMC(K)+FALTNDC*dttmp
      DUMFNC(K) = DUMFNC(K)+FALTNDNC*dttmp
      end do
      PRECRT = PRECRT+(FALOUTC(KTS))  &
                     *dttmp
      end do

! SNOW

      dttmp = DT/NSTEPS
      do K = KTS,KTE
        temr(k) = 1.0/(NSTEPS*RHO(k))
      enddo

      do N = 1,NSTEPS

      do K = KTS,KTE
      FALOUTS(K) = FS(K)*DUMQS(K)
      FALOUTNS(K) = FNS(K)*DUMFNS(K)
      end do
! TOP OF MODEL
      K = KTE
      FALTNDS = FALOUTS(K)*temz(k)
      FALTNDNS = FALOUTNS(K)*temz(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
      QNISTEN(K) = QNISTEN(K)-FALTNDS*temr(k)
      NS3DTEN(K) = NS3DTEN(K)-FALTNDNS*temr(k)
      DUMQS(K) = DUMQS(K)-FALTNDS*dttmp
      DUMFNS(K) = DUMFNS(K)-FALTNDNS*dttmp
      do K = KTE-1,KTS,-1
      FALTNDS = (FALOUTS(K+1)-FALOUTS(K))*temz(K)
      FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))*temz(K)
      QNISTEN(K) = QNISTEN(K)+FALTNDS*temr(k)
      NS3DTEN(K) = NS3DTEN(K)+FALTNDNS*temr(k)
      DUMQS(K) = DUMQS(K)+FALTNDS*dttmp
      DUMFNS(K) = DUMFNS(K)+FALTNDNS*dttmp
      end do
      PRECRT = PRECRT+(FALOUTS(KTS))  &
                     *dttmp
      SNOWRT = SNOWRT+(FALOUTS(KTS))*dttmp
      end do

! GRAUPEL

      dttmp = DT/NSTEPG
      do K = KTS,KTE
        temr(k) = 1.0/(NSTEPG*RHO(k))
      enddo

      do N = 1,NSTEPG

      do K = KTS,KTE
      FALOUTG(K) = FG(K)*DUMG(K)
      FALOUTNG(K) = FNG(K)*DUMFNG(K)
      end do
! TOP OF MODEL
      K = KTE
      FALTNDG = FALOUTG(K)*temz(k)
      FALTNDNG = FALOUTNG(K)*temz(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
      QGSTEN(K) = QGSTEN(K)-FALTNDG*temr(k)
      NG3DTEN(K) = NG3DTEN(K)-FALTNDNG*temr(k)
      DUMG(K) = DUMG(K)-FALTNDG*dttmp
      DUMFNG(K) = DUMFNG(K)-FALTNDNG*dttmp
      do K = KTE-1,KTS,-1
      FALTNDG = (FALOUTG(K+1)-FALOUTG(K))*temz(K)
      FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))*temz(K)
      QGSTEN(K) = QGSTEN(K)+FALTNDG*temr(k)
      NG3DTEN(K) = NG3DTEN(K)+FALTNDNG*temr(k)
      DUMG(K) = DUMG(K)+FALTNDG*dttmp
      DUMFNG(K) = DUMFNG(K)+FALTNDNG*dttmp
      end do
      PRECRT = PRECRT+(FALOUTG(KTS))  &
                     *dttmp
      SNOWRT = SNOWRT+(FALOUTG(KTS))*dttmp
      end do

      kloop:  &
      do K=KTS,KTE

! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

        QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
        QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
        QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
        QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
        QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)

! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

! V1.7
!hm 7/9/09 bug fix
!        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
        if (QI3D(K).GE.QSMALL.AND.T3D(K).LT.TMELT.AND.LAMI(K).GE.1.E-10) then

        if (1./LAMI(K).GE.2.*DCS) then
           QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
           NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
           QI3DTEN(K) = -QI3D(K)/DT
           NI3DTEN(K) = -NI3D(K)/DT
        end if
        end if

! hm add tendencies here, then call sizeparameter
! to ensure consisitency between mixing ratio and number concentration

          QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
          QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
          QNI3D(k)        = QNI3D(k)+QNI3DTEN(k)*DT
          QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
          NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
          NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
          NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
          NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT

          if (IGRAUP.EQ.0) then
          QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
          NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
          end if

! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
          T3D(K)         = T3D(K)+T3DTEN(k)*DT
          QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT
          if(eqtset.eq.2) PRES(K) = RHO(K)*T3D(K)*(R+RV*QV3D(K))

! SATURATION VAPOR PRESSURE AND MIXING RATIO

!!!            EVS(K) = POLYSVP(T3D(K),0)   ! PA
!!!            EIS(K) = POLYSVP(T3D(K),1)   ! PA
            EVS(K) = min( PRES(K)*0.5 , 611.2 * EXP( 17.67 * ( T3D(K) - 273.15 ) / ( T3D(K) - 29.65 ) ) )
            EIS(K) = min( PRES(K)*0.5 , 611.2 * EXP( 21.8745584 * ( T3D(K) - 273.15 ) / ( T3D(K) - 7.66 ) ) )

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            if (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! V1.3, change limit from 10^-7 to 10^-6
! V1.7 7/9/09 change limit from 10^-6 to 10^-8
! this improves reflectivity at low mixing ratios
! ghb 160428:  use QTHR

             if (QVQVS(K).LT.0.9) then
               if (QR3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QR3D(K)
!!!                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  T3D(K)=T3D(K)-QR3D(K)*ecnd(k)
                  QR3D(K)=0.
               end if
               if (QC3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QC3D(K)
!!!                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  T3D(K)=T3D(K)-QC3D(K)*ecnd(k)
                  QC3D(K)=0.
               end if
             end if

             if (QVQVSI(K).LT.0.9) then
               if (QI3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QI3D(K)
!!!                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  T3D(K)=T3D(K)-QI3D(K)*edep(k)
                  QI3D(K)=0.
               end if
               if (QNI3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QNI3D(K)
!!!                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  T3D(K)=T3D(K)-QNI3D(K)*edep(k)
                  QNI3D(K)=0.
               end if
               if (QG3D(K).LT.QTHR) then
                  QV3D(K)=QV3D(K)+QG3D(K)
!!!                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  T3D(K)=T3D(K)-QG3D(K)*edep(k)
                  QG3D(K)=0.
               end if
             end if

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       if (QC3D(K).LT.QSMALL) then
         QC3D(K) = 0.
         NC3D(K) = 0.
       end if
       if (QR3D(K).LT.QSMALL) then
         QR3D(K) = 0.
         NR3D(K) = 0.
       end if
       if (QI3D(K).LT.QSMALL) then
         QI3D(K) = 0.
         NI3D(K) = 0.
       end if
       if (QNI3D(K).LT.QSMALL) then
         QNI3D(K) = 0.
         NS3D(K) = 0.
       end if
       if (QG3D(K).LT.QSMALL) then
         QG3D(K) = 0.
         NG3D(K) = 0.
       end if

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

            if (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) GOTO 500

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE INSTANTANEOUS PROCESSES

! ADD MELTING OF CLOUD ICE TO FORM RAIN

        if (QI3D(K).GE.QSMALL.AND.T3D(K).GE.TMELT) then
           QR3D(K) = QR3D(K)+QI3D(K)
!!!           T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
           T3D(K) = T3D(K)-QI3D(K)*efrz(k)
           tbmelt(k) = tbmelt(k)-QI3D(K)*efrz(k)*rdt
           QI3D(K) = 0.
           NR3D(K) = NR3D(K)+NI3D(K)
           NI3D(K) = 0.
        end if

! ****SENSITIVITY - NO ICE
        if (ILIQ.EQ.1) GOTO 778

! HOMOGENEOUS FREEZING OF CLOUD WATER

        if (T3D(K).LE.233.15.AND.QC3D(K).GE.QSMALL) then
           QI3D(K)=QI3D(K)+QC3D(K)
!!!           T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
           T3D(K)=T3D(K)+QC3D(K)*efrz(k)
           tbfrz(k) = tbfrz(k)+QC3D(K)*efrz(k)*rdt
           QC3D(K)=0.
           NI3D(K)=NI3D(K)+NC3D(K)
           NC3D(K)=0.
        end if

! HOMOGENEOUS FREEZING OF RAIN

        if (IGRAUP.EQ.0) then

        if (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) then
           QG3D(K) = QG3D(K)+QR3D(K)
!!!           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
           T3D(K) = T3D(K)+QR3D(K)*efrz(k)
           tbfrz(k) = tbfrz(k)+QR3D(K)*efrz(k)*rdt
           QR3D(K) = 0.
           NG3D(K) = NG3D(K)+ NR3D(K)
           NR3D(K) = 0.
        end if

        else if (IGRAUP.EQ.1) then

        if (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) then
           QNI3D(K) = QNI3D(K)+QR3D(K)
!!!           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
           T3D(K) = T3D(K)+QR3D(K)*efrz(k)
           tbfrz(k) = tbfrz(k)+QR3D(K)*efrz(k)*rdt
           QR3D(K) = 0.
           NS3D(K) = NS3D(K)+NR3D(K)
           NR3D(K) = 0.
        end if

        end if

 778    continue

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

      if (QI3D(K).GE.QSMALL) then
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMI(K).LT.LAMMINI) then

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      else if (LAMI(K).GT.LAMMAXI) then
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      end if
      end if

!......................................................................
! RAIN

      if (QR3D(K).GE.QSMALL) then
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMR(K).LT.LAMMINR) then

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      else if (LAMR(K).GT.LAMMAXR) then
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      end if

      end if

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      if (QC3D(K).GE.QSMALL) then

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      if (LAMC(K).LT.LAMMIN) then
      LAMC(K) = LAMMIN
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      else if (LAMC(K).GT.LAMMAX) then
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      end if

      end if

!......................................................................
! SNOW

      if (QNI3D(K).GE.QSMALL) then
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMS(K).LT.LAMMINS) then
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      else if (LAMS(K).GT.LAMMAXS) then

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1
      NS3D(K) = N0S(K)/LAMS(K)
      end if

      end if

!......................................................................
! GRAUPEL

      if (QG3D(K).GE.QSMALL) then
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)

! CHECK FOR SLOPE

! ADJUST VARS

      if (LAMG(K).LT.LAMMING) then
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      else if (LAMG(K).GT.LAMMAXG) then

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      end if

      end if

 500  continue

      ifrad:  if( dorad )then
!! CALCULATE EFFECTIVE RADIUS

! hm 8/20/12
      if (QI3D(K).GE.QSMALL) then
         EFFI(K) = 3./LAMI(K)/2.*1.E6
      else
         EFFI(K) = 25.
      end if

      if (QNI3D(K).GE.QSMALL) then
         EFFS(K) = 3./LAMS(K)/2.*1.E6
      else
         EFFS(K) = 25.
      end if

      if (QR3D(K).GE.QSMALL) then
         EFFR(K) = 3./LAMR(K)/2.*1.E6
      else
         EFFR(K) = 25.
      end if

      if (QC3D(K).GE.QSMALL) then
      EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
             GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
      else
      EFFC(K) = 25.
      end if

      if (QG3D(K).GE.QSMALL) then
         EFFG(K) = 3./LAMG(K)/2.*1.E6
      else
         EFFG(K) = 25.
      end if

! add combined cloud ice + snow effective radius, used by Goddard radiation

      if (qi3d(k).ge.qsmall.and.qni3d(k).ge.qsmall) then
         effis(k) = 3./2.*1.e6*(1./lami(k)**4+1./lams(k)**4)/        &
                    (1./lami(k)**3+1./lams(k)**3)
      end if
! if there is no cloud ice but snow is present, use snow eff radius, and vice versa
      if (qi3d(k).ge.qsmall.and.qni3d(k).lt.qsmall) then
         effis(k) = effi(k)
      end if
      if (qi3d(k).lt.qsmall.and.qni3d(k).ge.qsmall) then
         effis(k) = effs(k)
      end if
      if (qi3d(k).lt.qsmall.and.qni3d(k).lt.qsmall) then
         effis(k) = 25.
      end if
      endif  ifrad

!---------------------------------!

! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
!          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
! HM, 3/4/13, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
! OF EXCESSIVE AND PERSISTENT ANVIL
! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
          NI3D(K) = MIN(NI3D(K),0.3E6/RHO(K))

! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
          if (INUM.EQ.0.AND.IACT.EQ.2) then
          NC3D(K) = MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K))
          end if
! SWITCH FOR CONSTANT DROPLET NUMBER
          if (INUM.EQ.1) then
! CHANGE NDCNST FROM CM-3 TO KG-1
             NC3D(K) = NDCNST*1.E6/RHO(K)
          end if

      end do  kloop !!! K LOOP

 400         continue

! ALL DONE !!!!!!!!!!!
      return
      end  subroutine M2005MICRO_GRAUPEL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!      REAL FUNCTION POLYSVP (T,TYPE)
!
!!-------------------------------------------
!
!!  COMPUTE SATURATION VAPOR PRESSURE
!
!!  POLYSVP RETURNED IN UNITS OF PA.
!!  T IS INPUT IN UNITS OF K.
!!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!
!      IMPLICIT NONE
!
!      REAL DUM
!      REAL T
!      INTEGER TYPE
!
!! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
!
!! ice
!      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
!      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
!	6.11147274, 0.503160820, 0.188439774e-1, &
!        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
!        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	
!
!! liquid
!      real a0,a1,a2,a3,a4,a5,a6,a7,a8 
!
!! V1.7
!      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
!	6.11239921, 0.443987641, 0.142986287e-1, &
!        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
!        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
!      real dt
!
!! ICE
!
!      IF (TYPE.EQ.1) THEN
!
!!         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
!!          LOG10(6.1071))*100.
!
!
!!!!      dt = max(-80.,t-273.16)
!!!!      polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
!!!!      polysvp = polysvp*100.
!
!        ! GHB ... needed for consistency with CM1 satadj scheme
!        polysvp = 611.2 * EXP( 21.8745584 * ( T  - 273.15 ) / ( T  - 7.66 ) )
!
!
!      END IF
!
!! LIQUID
!
!      IF (TYPE.EQ.0) THEN
!
!!!!       dt = max(-80.,t-273.16)
!!!!       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
!!!!       polysvp = polysvp*100.
!
!!         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
!!             5.02808*LOG10(373.16/T)-									&
!!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
!!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
!!             LOG10(1013.246))*100.
!
!           ! GHB ... needed for consistency with CM1 satadj scheme
!           polysvp = 611.2 * EXP( 17.67 * ( T  - 273.15 ) / ( T  - 29.65 ) )
!
!         END IF
!
!
!      END FUNCTION POLYSVP

!------------------------------------------------------------------------------

      real function GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      integer I,N
      logical PARITY
      real                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      real, dimension(7) :: C
      real, dimension(8) :: P
      real, dimension(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      data ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      data XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      data P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      data Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      data C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,				   &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = real(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      if(Y.LE.ZERO)then
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        if(RES.NE.ZERO)then
          if(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        else
          RES=XINF
          GOTO 900
        endif
      endif
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      if(Y.LT.EPS)then
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        if(Y.GE.XMININ)then
          RES=ONE/Y
        else
          RES=XINF
          GOTO 900
        endif
      elseif(Y.LT.TWELVE)then
        Y1=Y
        if(Y.LT.ONE)then
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        else
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        endif
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        do I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        end do
        RES=XNUM/XDEN+ONE
        if(Y1.LT.Y)then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        elseif(Y1.GT.Y)then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          do I=1,N
            RES=RES*Y
            Y=Y+ONE
          end do
        endif
      else
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        if(Y.LE.XBIG)then
          YSQ=Y*Y
          SUM=C(7)
          do I=1,6
            SUM=SUM/YSQ+C(I)
          end do
          SUM=SUM/Y-Y+xxx
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        else
          RES=XINF
          GOTO 900
        endif
      endif
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      if(PARITY)RES=-RES
      if(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      return
! ---------- LAST LINE OF GAMMA ----------
      end function GAMMA


      real function DERF1(X)
      implicit none
      real X
      real, dimension(0 : 64) :: A, B
      real W,T,Y
      integer K,I
      data A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,	                          &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,							  &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0, 							  &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,							  &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      data (B(I), I = 0, 12) /                                  &
         -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
         -0.00000000640788827665E0, -0.00000001667820132100E0,  &
         -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
          0.00001612722157047886E0, -0.00025616361025506629E0, 	&
          0.00015380842432375365E0,  0.00815533022524927908E0, 	&
         -0.01402283663896319337E0, -0.19746892495383021487E0,  &
          0.71511720328842845913E0 /
      data (B(I), I = 13, 25) /                                 &
         -0.00000000001951073787E0, -0.00000000032302692214E0,  &
          0.00000000522461866919E0,  0.00000000342940918551E0, 	&
         -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
          0.00002687044575042908E0, -0.00011843240273775776E0, 	&
         -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
          0.00909530922354827295E0, -0.20160072778491013140E0, 	&
          0.51169696718727644908E0 /
      data (B(I), I = 26, 38) /                                 &
         0.00000000003147682272E0, -0.00000000048465972408E0,   &
         0.00000000063675740242E0,  0.00000003377623323271E0, 	&
        -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
         0.00001947204525295057E0,  0.00002854147231653228E0, 	&
        -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
         0.02328095035422810727E0, -0.16725021123116877197E0, 	&
         0.32490054966649436974E0 /
      data (B(I), I = 39, 51) /                                 &
         0.00000000002319363370E0, -0.00000000006303206648E0,   &
        -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
         0.00000011371857327578E0, -0.00000211211337219663E0, 	&
         0.00000368797328322935E0,  0.00009823686253424796E0, 	&
        -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
         0.02585434424202960464E0, -0.11637092784486193258E0, 	&
         0.18267336775296612024E0 /
      data (B(I), I = 52, 64) /                                 &
        -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
        -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
         0.00000018006992266137E0, -0.00000078441223763969E0, 	&
        -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
        -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
         0.02064129023876022970E0, -0.06905562880005864105E0,   &
         0.09084526782065478489E0 /
      W = ABS(X)
      if (W .LT. 2.2D0) then
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
              A(K + 11)) * T + A(K + 12)) * W
      else if (W .LT. 6.9D0) then
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      else
          Y = 1
      end if
      if (X .LT. 0) Y = -Y
      DERF1 = Y
      end function DERF1

!+---+-----------------------------------------------------------------+
!
      subroutine radar_init

      implicit none
      integer:: n
      PI5 = PI*PI*PI*PI*PI
      lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
      m_w_0 = m_complex_water_ray (lamda_radar, 0.0d0)
      m_i_0 = m_complex_ice_maetzler (lamda_radar, 0.0d0)
      K_w = (ABS( (m_w_0*m_w_0 - 1.0) /(m_w_0*m_w_0 + 2.0) ))**2

      do n = 1, nbins+1
         simpson(n) = 0.0d0
      enddo
      do n = 1, nbins-1, 2
         simpson(n) = simpson(n) + basis(1)
         simpson(n+1) = simpson(n+1) + basis(2)
         simpson(n+2) = simpson(n+2) + basis(3)
      enddo

      do n = 1, slen
         mixingrulestring_s(n:n) = char(0)
         matrixstring_s(n:n) = char(0)
         inclusionstring_s(n:n) = char(0)
         hoststring_s(n:n) = char(0)
         hostmatrixstring_s(n:n) = char(0)
         hostinclusionstring_s(n:n) = char(0)
         mixingrulestring_g(n:n) = char(0)
         matrixstring_g(n:n) = char(0)
         inclusionstring_g(n:n) = char(0)
         hoststring_g(n:n) = char(0)
         hostmatrixstring_g(n:n) = char(0)
         hostinclusionstring_g(n:n) = char(0)
      enddo

      mixingrulestring_s = 'maxwellgarnett'
      hoststring_s = 'air'
      matrixstring_s = 'water'
      inclusionstring_s = 'spheroidal'
      hostmatrixstring_s = 'icewater'
      hostinclusionstring_s = 'spheroidal'

      mixingrulestring_g = 'maxwellgarnett'
      hoststring_g = 'air'
      matrixstring_g = 'water'
      inclusionstring_g = 'spheroidal'
      hostmatrixstring_g = 'icewater'
      hostinclusionstring_g = 'spheroidal'

      end subroutine radar_init
!+---+-----------------------------------------------------------------+

      complex*16 function m_complex_water_ray(lambda,T)

!      Complex refractive Index of Water as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
!      after Ray (1972)

      implicit none
      double precision, intent(in):: T,lambda
      double precision:: epsinf,epss,epsr,epsi
      double precision:: alpha,lambdas,sigma,nenner
      complex*16, parameter:: i = (0d0,1d0)

      epsinf  = 5.27137d0 + 0.02164740d0 * T - 0.00131198d0 * T*T
      epss    = 78.54d+0 * (1.0 - 4.579d-3 * (T - 25.0)                 &
              + 1.190d-5 * (T - 25.0)*(T - 25.0)                        &
              - 2.800d-8 * (T - 25.0)*(T - 25.0)*(T - 25.0))
      alpha   = -16.8129d0/(T+273.16) + 0.0609265d0
      lambdas = 0.00033836d0 * exp(2513.98d0/(T+273.16)) * 1e-2

      nenner = 1.d0+2.d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PI*0.5) &
             + (lambdas/lambda)**(2d0-2d0*alpha)
      epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)   &
           * sin(alpha*PI*0.5)+1d0)) / nenner
      epsi = ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)            &
           * cos(alpha*PI*0.5)+0d0)) / nenner                           &
           + lambda*1.25664/1.88496
      
      m_complex_water_ray = SQRT(CMPLX(epsr,-epsi))
      
      end function m_complex_water_ray

!+---+-----------------------------------------------------------------+
      
      complex*16 function m_complex_ice_maetzler(lambda,T)
      
!      complex refractive index of ice as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.0001,30] m; T in [-250.0,0.0] C
!      Original comment from the Matlab-routine of Prof. Maetzler:
!      Function for calculating the relative permittivity of pure ice in
!      the microwave region, according to C. Maetzler, "Microwave
!      properties of ice and snow", in B. Schmitt et al. (eds.) Solar
!      System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
!      Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
!      TK = temperature (K), range 20 to 273.15
!      f = frequency in GHz, range 0.01 to 3000
         
      implicit none
      double precision, intent(in):: T,lambda
      double precision:: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa

      c = 2.99d8
      TK = T + 273.16
      f = c / lambda * 1d-9

      B1 = 0.0207
      B2 = 1.16d-11
      b = 335.0d0
      deltabeta = EXP(-10.02 + 0.0364*(TK-273.16))
      betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1)**2) ) + B2*f*f
      beta = betam + deltabeta
      theta = 300. / TK - 1.
      alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
      m_complex_ice_maetzler = 3.1884 + 9.1e-4*(TK-273.16)
      m_complex_ice_maetzler = m_complex_ice_maetzler                   &
                             + CMPLX(0.0d0, (alfa/f + beta*f)) 
      m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))
      
      end function m_complex_ice_maetzler
!+---+-----------------------------------------------------------------+

      subroutine rayleigh_soak_wetgraupel (x_g, a_geo, b_geo, fmelt,    &
                     meltratio_outside, m_w, m_i, lambda, C_back,       &
                     mixingrule,matrix,inclusion,                       &
                     host,hostmatrix,hostinclusion)

      implicit none

      double precision, intent(in):: x_g, a_geo, b_geo, fmelt, lambda,  &
                                     meltratio_outside
      double precision, intent(out):: C_back
      complex*16, intent(in):: m_w, m_i
      character(len=*), intent(in):: mixingrule, matrix, inclusion,     &
                                     host, hostmatrix, hostinclusion

      complex*16:: m_core, m_air
      double precision:: D_large, D_g, rhog, x_w, xw_a, fm, fmgrenz,    &
                         volg, vg, volair, volice, volwater,            &
                         meltratio_outside_grenz, mra
      integer:: error
      real :: rho_i, rho_w

      rho_i = 900.
      rho_w = 1000.


!     refractive index of air:
      m_air = (1.0d0,0.0d0)

!     Limiting the degree of melting --- for safety: 
      fm = DMAX1(DMIN1(fmelt, 1.0d0), 0.0d0)
!     Limiting the ratio of (melting on outside)/(melting on inside):
      mra = DMAX1(DMIN1(meltratio_outside, 1.0d0), 0.0d0)

!    ! The relative portion of meltwater melting at outside should increase
!    ! from the given input value (between 0 and 1)
!    ! to 1 as the degree of melting approaches 1,
!    ! so that the melting particle "converges" to a water drop.
!    ! Simplest assumption is linear:
      mra = mra + (1.0d0-mra)*fm

      x_w = x_g * fm

      D_g = a_geo * x_g**b_geo

      if (D_g .ge. 1d-12) then

       vg = PI/6. * D_g**3
       rhog = DMAX1(DMIN1(x_g / vg, DBLE(rho_i)), 10.0d0)
       vg = x_g / rhog
      
       meltratio_outside_grenz = 1.0d0 - rhog / rho_w

       if (mra .le. meltratio_outside_grenz) then
        !..In this case, it cannot happen that, during melting, all the
        !.. air inclusions within the ice particle get filled with
        !.. meltwater. This only happens at the end of all melting.
        volg = vg * (1.0d0 - mra * fm)
 
       else
        !..In this case, at some melting degree fm, all the air
        !.. inclusions get filled with meltwater.
        fmgrenz=(rho_i-rhog)/(mra*rho_i-rhog+rho_i*rhog/rho_w)

        if (fm .le. fmgrenz) then
         !.. not all air pockets are filled:
         volg = (1.0 - mra * fm) * vg
        else
         !..all air pockets are filled with meltwater, now the
         !.. entire ice sceleton melts homogeneously:
         volg = (x_g - x_w) / rho_i + x_w / rho_w
        endif

       endif

       D_large  = (6.0 / PI * volg) ** (1./3.)
       volice = (x_g - x_w) / (volg * rho_i)
       volwater = x_w / (rho_w * volg)
       volair = 1.0 - volice - volwater
      
       !..complex index of refraction for the ice-air-water mixture
       !.. of the particle:
       m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,      &
                         volwater, mixingrule, host, matrix, inclusion, &
                         hostmatrix, hostinclusion, error)
       if (error .ne. 0) then
        C_back = 0.0d0
        return
       endif

       !..Rayleigh-backscattering coefficient of melting particle: 
       C_back = (ABS((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2           &
                * PI5 * D_large**6 / lamda4

      else
       C_back = 0.0d0
      endif

      end subroutine rayleigh_soak_wetgraupel
!+---+-----------------------------------------------------------------+

      complex*16 function get_m_mix_nested (m_a, m_i, m_w, volair,      &
                     volice, volwater, mixingrule, host, matrix,        &
                     inclusion, hostmatrix, hostinclusion, cumulerror)

      implicit none

      double precision, intent(in):: volice, volair, volwater
      complex*16, intent(in):: m_a, m_i, m_w
      character(len=*), intent(in):: mixingrule, host, matrix,          &
                     inclusion, hostmatrix, hostinclusion
      integer, intent(out):: cumulerror

      double precision:: vol1, vol2
      complex*16:: mtmp
      integer:: error

      !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
      !.. air, ice, or water

      cumulerror = 0
      get_m_mix_nested = CMPLX(1.0d0,0.0d0)

      if (host .eq. 'air') then

       if (matrix .eq. 'air') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
!!!        CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, 0.0d0, vol1, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error
          
        if (hostmatrix .eq. 'air') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'icewater') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
!!!         CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'ice') then

       if (matrix .eq. 'ice') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
!!!        CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0d0, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'ice') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airwater') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         'air', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
!!!         CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'water') then

       if (matrix .eq. 'water') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
!!!        CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0d0,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'water') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,            &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airice') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,            &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
                           hostmatrix
!!!         CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'none') then

       get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                     &
                       volair, volice, volwater, mixingrule,            &
                       matrix, inclusion, error)
       cumulerror = cumulerror + error
        
      else
       write(mp_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
!!!       CALL wrf_debug(150, mp_debug)
       cumulerror = cumulerror + 1
      endif

      if (cumulerror .ne. 0) then
       write(mp_debug,*) 'GET_M_MIX_NESTED: error encountered'
!!!       CALL wrf_debug(150, mp_debug)
       get_m_mix_nested = CMPLX(1.0d0,0.0d0)    
      endif

      end function get_m_mix_nested

!+---+-----------------------------------------------------------------+

      complex*16 function get_m_mix (m_a, m_i, m_w, volair, volice,     &
                     volwater, mixingrule, matrix, inclusion, error)

      implicit none

      double precision, intent(in):: volice, volair, volwater
      complex*16, intent(in):: m_a, m_i, m_w
      character(len=*), intent(in):: mixingrule, matrix, inclusion
      integer, intent(out):: error

      error = 0
      get_m_mix = CMPLX(1.0d0,0.0d0)

      if (mixingrule .eq. 'maxwellgarnett') then
       if (matrix .eq. 'ice') then
        get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,  &
                           m_i, m_a, m_w, inclusion, error)
       elseif (matrix .eq. 'water') then
        get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,  &
                           m_w, m_a, m_i, inclusion, error)
       elseif (matrix .eq. 'air') then
        get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,  &
                           m_a, m_w, m_i, inclusion, error)
       else
        write(mp_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
!!!        CALL wrf_debug(150, mp_debug)
        error = 1
       endif

      else
       write(mp_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
!!!       CALL wrf_debug(150, mp_debug)
       error = 2
      endif

      if (error .ne. 0) then
       write(mp_debug,*) 'GET_M_MIX: error encountered'
!!!       CALL wrf_debug(150, mp_debug)
      endif

      end function get_m_mix

!+---+-----------------------------------------------------------------+

      complex*16 function m_complex_maxwellgarnett(vol1, vol2, vol3,    &
                     m1, m2, m3, inclusion, error)

      implicit none

      complex*16 :: m1, m2, m3
      double precision :: vol1, vol2, vol3
      character(len=*) :: inclusion

      complex*16 :: beta2, beta3, m1t, m2t, m3t
      integer, intent(out) :: error

      error = 0

      if (DABS(vol1+vol2+vol3-1.0d0) .gt. 1d-6) then
       write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',       &
              'partial volume fractions is not 1...ERROR'
!!!       CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=CMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m1t = m1**2
      m2t = m2**2
      m3t = m3**2

      if (inclusion .eq. 'spherical') then
       beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
       beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
      elseif (inclusion .eq. 'spheroidal') then
       beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
       beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
      else
       write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ',                  &
                         'unknown inclusion: ', inclusion
!!!       CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=DCMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m_complex_maxwellgarnett = &
       SQRT(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
       (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

      end function m_complex_maxwellgarnett

!+---+-----------------------------------------------------------------+
!..Compute radar reflectivity assuming 10 cm wavelength radar and using
!.. Rayleigh approximation.  Only complication is melted snow/graupel
!.. which we treat as water-coated ice spheres and use Uli Blahak's
!.. library of routines.  The meltwater fraction is simply the amount
!.. of frozen species remaining from what initially existed at the
!.. melting level interface.
!+---+-----------------------------------------------------------------+
      subroutine calc_refl10cm (qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ,  &
                          kts, kte, ii, jj, nr1d, ns1d, ng1d, rho1d)

      implicit none

!..Sub arguments
      integer, intent(in):: kts, kte, ii, jj
      real, dimension(kts:kte), intent(in)::                            &
                qv1d, qr1d, qs1d, qg1d, t1d, p1d, nr1d, ns1d, ng1d, rho1d
      real, dimension(kts:kte), intent(inout):: dBZ

!..Local variables
      real, dimension(kts:kte):: temp, pres, qv, rho
      real, dimension(kts:kte):: rr, rs, rg,rnr,rns,rng

      double precision, dimension(kts:kte):: ilamr, ilamg, N0_r, N0_g,ilams,n0_s

      real, dimension(kts:kte):: ze_rain, ze_snow, ze_graupel

      double precision:: lamg
      double precision:: fmelt_s, fmelt_g

      integer:: i, k, k_0
      logical:: melti
      logical, dimension(kts:kte):: L_qr, L_qs, L_qg

!..Single melting snow/graupel particle 70% meltwater on external sfc
      double precision, parameter:: melt_outside_s = 0.7d0
      double precision, parameter:: melt_outside_g = 0.7d0

      double precision:: cback, x, eta, f_d

! hm added parameter
      real R1,t_0,dumlams,dumlamr,dumlamg,dumn0s,dumn0r,dumn0g,ocms,obms,ocmg,obmg

      integer n

      R1 = 1.E-12
      t_0 = 273.15

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
!!!         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rho(k) = rho1d(K)
         if (qr1d(k) .gt. R1) then
            rr(k) = qr1d(k)*rho(k)
            L_qr(k) = .true.
         else
            rr(k) = R1
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R1) then
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R1) then
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif

! hm add number concentration
         if (nr1d(k) .gt. R1) then
            rnr(k) = nr1d(k)*rho(k)
         else
            rnr(k) = R1
         endif
         if (ns1d(k) .gt. R1) then
            rns(k) = ns1d(k)*rho(k)
         else
            rns(k) = R1
         endif
         if (ng1d(k) .gt. R1) then
            rng(k) = ng1d(k)*rho(k)
         else
            rng(k) = R1
         endif

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte

! compute moments for snow

! calculate slope and intercept parameter

      dumLAMS = (CONS1*rns(K)/rs(K))**(1./DS)
      dumN0S = rns(K)*dumLAMS/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      if (dumLAMS.LT.LAMMINS) then
      dumLAMS = LAMMINS
      dumN0S = dumLAMS**4*rs(K)/CONS1
      else if (dumLAMS.GT.LAMMAXS) then
      dumLAMS = LAMMAXS
      dumN0S = dumLAMS**4*rs(k)/CONS1
      end if

      ilams(k)=1./dumlams
      n0_s(k)=dumn0s

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

      do k = kte, kts, -1


! calculate slope and intercept parameter

      dumLAMg = (CONS2*rng(K)/rg(K))**(1./Dg)
      dumN0g = rng(K)*dumLAMg/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      if (dumLAMg.LT.LAMMINg) then
      dumLAMg = LAMMINg
      dumN0g = dumLAMg**4*rg(K)/CONS2
      else if (dumLAMg.GT.LAMMAXg) then
      dumLAMg = LAMMAXg
      dumN0g = dumLAMg**4*rg(k)/CONS2
      end if

      ilamg(k)=1./dumlamg
      n0_g(k)=dumn0g

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept & slope values for rain.
!+---+-----------------------------------------------------------------+

      do k = kte, kts, -1

! calculate slope and intercept parameter

      dumLAMr = (PI*RHOW*rnr(K)/rr(K))**(1./3.)
      dumN0r = rnr(K)*dumLAMr/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      if (dumLAMr.LT.LAMMINr) then
      dumLAMr = LAMMINr
      dumN0r = dumLAMr**4*rr(K)/(PI*RHOW)
      else if (dumLAMr.GT.LAMMAXr) then
      dumLAMr = LAMMAXr
      dumN0r = dumLAMr**4*rr(k)/(PI*RHOW)
      end if

      ilamr(k)=1./dumlamr
      n0_r(k)=dumn0r

      enddo

      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt. T_0) .and. (rr(k).gt. 0.001e-3) &
                   .and. ((rs(k+1)+rg(k+1)).gt. 0.01e-3) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*720.*ilamr(k)**7

         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (pi*rhosn/6./900.)*(pi*rhosn/6./900.) &
                                    * N0_s(k)*720.*ilams(k)**7
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (pi*rhog/6./900.)* (pi*rhog/6./900.)        &
                                    * N0_g(k)*720.*ilamg(k)**7
      enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.2) then
       do k = k_0-1, 1, -1

!..Reflectivity contributed by melting snow
          fmelt_s = DMIN1(1.0d0-rs(k)/rs(k_0), 1.0d0)
          if (fmelt_s.gt.0.01d0 .and. fmelt_s.lt.0.99d0 .and.           &
                         rs(k).gt.R1) then
           eta = 0.d0
           obms = 1./ds
           ocms = (1./(pi*rhosn/6.))**obms
           do n = 1, nbs
              x = pi*rhosn/6. * Dds(n)**3
              call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)* DEXP(-Dds(n)/ilams(k))
              eta = eta + f_d * CBACK * simpson(n) * dts(n)

           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif


!..Reflectivity contributed by melting graupel

          fmelt_g = DMIN1(1.0d0-rg(k)/rg(k_0), 1.0d0)
          if (fmelt_g.gt.0.01d0 .and. fmelt_g.lt.0.99d0 .and.           &
                         rg(k).gt.R1) then
           eta = 0.d0
           lamg = 1./ilamg(k)
           obmg = 1./dg
           ocmg = (1./(pi*rhog/6.))**obmg
           do n = 1, nbg
              x = pi*rhog/6. * Ddg(n)**3
              call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)* DEXP(-lamg*Ddg(n))
              eta = eta + f_d * CBACK * simpson(n) * dtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo


      end subroutine calc_refl10cm

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine satadj2(tnew,prs,qvnew,qcnew,qrnew,qinew,qsnew,qgnew,rho,eqtset,tsmall)
      implicit none

      ! CM1 iterative saturation adjustment scheme for Morrison microphysics

      real, intent(inout) :: tnew,prs,qvnew,qcnew
      real, intent(in) :: qrnew,qinew,qsnew,qgnew,rho
      integer, intent(in) :: eqtset
      real, intent(in) :: tsmall

      integer :: n,iflag
      real :: esl,qvs,cvml,rm,lhv,tlast,dqv
      real :: converge,t1,d1,tem,ql,qi
      logical :: doit

!--------------------------------------------------------------------
!  iterative sat adj.

    iflag=0

    if(eqtset.eq.2)then

        converge=2.0*tsmall

        esl=min( prs*0.5 , 611.2*exp( 17.67 * ( tnew - 273.15 ) / ( tnew - 29.65 ) ) )
        qvs=ep_2*esl/(prs-esl)

        if(qcnew.gt.qsmall .or. qvnew.gt.qvs)then

          ql=max(0.0,qcnew+qrnew)
          qi=max(0.0,qinew+qsnew+qgnew)
          cvml=cv+cvv*qvnew+cpl*ql+cpi*qi
          lhv=lv1-lv2*tnew

          t1=(lhv-rv*tnew)/cvml
          d1=t1*17.67*243.5

          n=0
          tlast=tnew
          doit=.true.

          do while( doit )
            n=n+1
            dqv=(qvs-qvnew)/(1.0+d1*qvs/((tnew-29.65)**2) )
            dqv=min(dqv,qcnew)
            if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew

            qvnew=qvnew+dqv
            qcnew=qcnew-dqv
            tnew=tnew-dqv*t1
            prs=rho*(R+rv*qvnew)*tnew

            doit = .false.
            if( abs(tnew-tlast).gt.converge )then
              tlast=tnew
              esl=min( prs*0.5 , 611.2*exp( 17.67 * ( tnew - 273.15 ) / ( tnew - 29.65 ) ) )
              qvs=ep_2*esl/(prs-esl)
              doit = .true.
            endif

            if(n.gt.50) print *,'  satadj2:',n,tnew,prs
            if(n.eq.100)then
              print *,'  infinite loop!'
              iflag=1
              doit=.false.
            endif

          enddo

        endif

    else

        esl=min( prs*0.5 , 611.2*exp( 17.67 * ( tnew - 273.15 ) / ( tnew - 29.65 ) ) )
        qvs=ep_2*esl/(prs-esl)
        lhv=lv1-lv2*tnew

        dqv=(qvs-qvnew)/(1.0+lhv*qvs*17.67*243.5/(cp*((tnew-29.65)**2)))
        if(  (qvnew+dqv).lt.1.0e-20 ) dqv=1.0e-20-qvnew
        dqv=min(dqv,max(0.0,qcnew))

        qvnew=qvnew+dqv
        qcnew=qcnew-dqv
        tnew=tnew-dqv*lhv/(cp*(1.+0.887*qvnew))

    endif

      if(iflag.ne.0)then
        print *
        print *,' Convergence cannot be reached in satadj2 subroutine.'
        print *
        print *,' This may be a problem with the algorithm in satadj2.'
        print *,' However, the model may have became unstable somewhere'
        print *,' else and the symptoms first appeared here.'
        print *
        print *,' Try decreasing the timestep.'
        print *
        print *,'  ... stopping cm1 ... '
        print *
        call stopcm1
      endif

      end subroutine satadj2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module module_mp_GRAUPEL
