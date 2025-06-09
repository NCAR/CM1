!COMMAS/1:MODEL_LAYER:PHYSICS

! prepocessed on "Aug 14 2023" at "11:12:05" - update Oct 29 2023







!---------------------------------------------------------------------
! IMPORTANT: Best results are attained using the 5th-order WENO (Weighted Essentially Non-Oscillatory) advection option (4) for scalars:
! moist_adv_opt                       = 4,
! scalar_adv_opt                      = 4, (can also use option 3, which is WENO without the positive definite filter)
! The WENO-5 scheme provides a 5th-order (horizontal and vertical) adaptive weighting of components that 
! better preserve monotinicity in strong gradients. The standard 5th-order formulation is prone to undershoots 
! (negative values) of mass and number concentrations at cloud edges. The WENO scheme helps 
! to prevent undershoots and results in less noise at cloud and reflectivity boundaries. This is particularly
! useful for multi-moment schemes to preserve relationships between mass and number concentration. An option is also available
! for WENO-5 advection of momentum, but this can result in excessive damping of poorly-resolved features. For both scalar and momentum
! the steps 1 and 2 of the Runge-Kutta time integration use standare 5th-order advection, and the WENO-5 is applied on the 3rd (final)
! RK step. Option 3 applies the WENO-5, and option 4 adds the positive definite filter (as also used in option 1). 
!
! WENO references: Jiang and Shu, 1996, J. Comp. Phys. v. 126, 202-223; Shu 2003, Int. J. Comp. Fluid Dyn. v. 17 107-118;
!
!! This module provides a 1/2/3-moment bulk microphysics scheme based on a combination of
!! Straka and Mansell (2005, JAM) and Zeigler (1985, JAS) and modified/upgraded in
!! in Mansell, Zeigler, and Bruning (2010, JAS).  Two-moment adaptive sedimentation
!! follows Mansell (2010, JAS), using parameter infall = 4.
!!
!! Added info on graupel density and soaking is in Mansell and Ziegler (2013, JAS)
!!
!! Average graupel and hail particle densities are predicted, which affects fall speed as well.
!!
!! Maintainer: Ted Mansell, National Severe Storms Laboratory <ted.mansell@noaa.gov>
!!
!! Microphysics References:
!!
!! Mansell, E. R., C. L. Ziegler, and E. C. Bruning, 2010: Simulated electrification of a small
!!   thunderstorm with two-moment bulk microphysics. J. Atmos. Sci., 67, 171-194, doi:10. 1175/2009JAS2965.1.
!!
!!  Mansell, E. R. and C. L. Ziegler, 2013: Aerosol effects on simulated storm electrification and
!!     precipitation in a two-moment bulk microphysics model. J. Atmos. Sci., 70 (7), 2032-2050,
!!     doi:10.1175/JAS-D-12-0264.1.
!!
!! Ziegler, C. L., 1985: Retrieval of thermal and microphysical variables in observed convective storms.
!!    Part I: Model development and preliminary testing. J. Atmos. Sci., 42, 1487-1509.
!!
!! Sedimentation reference:
!!
!! Mansell, E. R., 2010: On sedimentation and advection in multimoment bulk microphysics.
!!    J. Atmos. Sci., 67, 3084-3094, doi:10.1175/2010JAS3341.1.
!
! Possible parameters to adjust:
!
!  ccn : base cloud condensation nuclei concentration (use namelist.input value "nssl_cccn")
!  alphah, alphahl : Size distribution shape parameters for graupel (h) and hail (hl)
!  infall : changes sedimentation options to see effects (see below)
!
! lightning model references:
!
!    Fierro, A. O., E.R. Mansell, C. Ziegler and D. R. MacGorman 2013: The
!    implementation of an explicit charging and discharge lightning scheme
!    within the WRF-ARW model: Benchmark simulations of a continental squall line, a
!    tropical cyclone and a winter storm. Monthly Weather Review, Volume 141, 2390-2415
!
!    Mansell et al. 2005: Charge structure and lightning sensitivity in a simulated
!     multicell thunderstorm. J. Geophys. Res., 110, D12101, doi:10.1029/2004JD005287
!
! Note: Some parameters below apply to unreleased features.
!
!
!---------------------------------------------------------------------
! Apr. 2023
!  - Update to 3-moment for rain, graupel, and hail
!  - Change default graupel/hail fall speeds to icdx/icdxhl=6 (Milbrandt & Morrison 2013)
!     and also set default ehw0=0.9 and ehlw0=0.9 to compensate for lower fall speeds.
!  - Change default hail conversion to ihlcnh=-1, and then =1 for 2-mom or =3 for 3-mom,
!     using wet growth diameter to convert large graupel
!---------------------------------------------------------------------
! Sept. 2021:
! Fixes:
!   Restored previous formulation of snow reflectivity, as it was realized that the last change incorrectly assumed a fixed
!     density independent of size. Generally lower snow reflectivity values as a result (no effect on microphysics)
! Other:
!   Generic fall speed coeffecients (axx,bxx) to accomodate future frozen drops category (no effect)
!   Reordered collection coefficients (dab1lh) to be consistent (no effect)
!   Switched to full calculation of rain number loss via collection by graupel (chacr; to be consisted with collection by hail) (minor effects)
!---------------------------------------------------------------------
! April 2021:
! Fixes:
!  Fall speed air density factor limited to air density of 0.05 (for very high model top) to mitigate excessive fall speeds
!  Fixed issue of spurious creation of large concentrations of very small droplets and transient large condensation (also increased minimum droplet size)
!  Fixed issue of negligible "seed" values of graupel from Bigg freezing at relatively high temperatures (thanks to S. Lasher-Trapp)
!  Minor bug fix in effective radius calculation of snow. (thanks to T. Iguchi)
! Updates:
!  Enabled regeneration of CCN by droplet evaporation and background restore (default time constant of 3600s)
!  Updated the routine that handles single-moment variables on the first time step. This sets a higher threshold for meaningful mixing ratios and sets a more realistic droplet concentration (also activating CCN as needed).
!  Enabled radar reflectivity from cloud ice (new formulation) ( idbzci = 1 )
!  Added internal option for ice crystal nucleation by DeMott et al. (2010, PNAS) (inucopt=4)
!  Allow greater fraction of hail to melt in one time step
!  Reduced minimum number concentration from 1e-4 to 1e-8 (based on CAPS input)
!  Added internal namelist for easier access to internal variables for development/testing and easier setup for ensemble microphysics diversity
!    (namelist read is disabled by default)
!  Increased resolution of lookup table for incomplete gamma functions
!
!---------------------------------------------------------------------
! Sept. 2019:
! Bug fixes:
!  - Effective radius calculation was only done at history times. Now every time step (though should be just before radiation is called)
!  - Snow reflectivity: Previous "fix" was incorrect and yields snow dBZ that is too low. Reverted to old version which was correct
!  - Incorrectly updated a state value in the reflectivity code. (Could cause small differences if reflectivity is not calculated)
! Updates:
!  - Added code hints to use the "axtra2d" array to communicate rates from the microphysics routine into any 3d arrays that are passed in to the driver.
!  - Graupel and hail drag coefficients are returned from fall speed subroutine to use in ventilation coeffs. for consistency (minor change)
!  - Added (compile) option flag to turn on diagnosis of cloud droplet shape parameter based on number concentration
!  - Added (compile) option flag icracr to turn off rain self-collection
!  - Added compile options 'depfac' and 'meltfac' to adjust deposition/sublimation and melting (not freezing) rates of graupel/hail by a constant factor (for experimentation). Default value is 1.0
!  - Put limit on snow volume (2 cm) in aggregation rate
!---------------------------------------------------------------------
! WRF 4.0 update:
!  Major:
!   Fixed excessive sublimation that could occur in very strong downdrafts (3.9.1.1 update)
!
!  Minor:
!    icefallopt=3 : New ice crystal fall speed that has faster speeds for small ice particles. Main effect
!                   is on anvil clouds to help them decay a bit faster. Old behavior can be recovered with icefallopt=1
!    Cosmetic: removed stray single quotes because some preprocessors complain about unclosed quotes even in comments
!
!---------------------------------------------------------------------
! WRF 3.9.1.1 update:
!
!  Added a check on overdepletion of ice by sublimation, which could sometimes result in water supersaturation
!  Bug fix: setting of t7 used 'dn' instead of 'dn1' (Thanks to Chunxi Zhang)
!
!---------------------------------------------------------------------
! WRF 3.9 updates:
!
!   2-moment scheme now creates number concentration tendencies from cumulus scheme mass mixing ratio rates
!   Renamed internal gamma function routine from 'gamma' to 'gamma_sp' to avoid name conflicts
!   Restored older settings that allow snow aggregation starting at T > -25C
!   Adjusted Meyers number of activated nuclei by the local air density to compensate for using data at surface
!   Minor updates to rain-ice crystal and hail-rain collection efficiencies
!
!   
!   Reduced minimum mean snow diameter from 100 microns to 10 microns
!
!---------------------------------------------------------------------
! WRF 3.8 updates:
!   Fixed issue with reflectivity conservation for graupel melting into rain. Rain number concentrations were too low,
!       resulting in excessive reflectivity of a couple dBZ
!   Changed default value of iusewetgraupel to 1 (turns off diagnostic meltwater on graupel for reflectivity)
!   Apply a 70 m/s fall speed limit for sedimentation
!   Changed vapor ice nucleation to Meyers-Ferrier method (original scheme)
!   New method for Bigg freezing (ibiggopt=2)
!   Reduced snow aggregration efficiency and restricted aggregation to higher temperatures (assuming dendrites and mechanical aggregation)
!   Increased maximum graupel-droplet collection efficiency when hail is turned off (nssl_2momg)
!   Updates for compatibility with WRF-NMM
!   Added calculation of hail number concentration in calcnfromq (creates number concentration from mixing ratio
!       when starting from an analysis). And fixed error in graupel intercept
!   Bug fix in snow fall speeds
!   Further fix in snow reflectivity
!   Use diameter of maximum mass rather than mean diamter when checking maximum size
!   Helped performance in sedimentation with flag "do_accurate_sedimentation" to control recalculation of fall speeds when
!       more than one sub-time step is needed (often happens with large time steps and small dz near the ground):
!        = .true. : recalculates fall speed after each substep (more accurate)
!        = .false. : (default) reuses fall speeds calculated on the first substep (typical for most schemes), theoretically could cause an occasional glitch, but none seen in practice
!   Increased maximum mean droplet radius from 40 to 60 microns, which alleviates spurious number concentration increases at low CCN concentration.
!   Removed a duplicate factor from hail reflectivity that was causing a loss of about 6 dBZ (since WRF 3.5).
!
!---------------------------------------------------------------------



module module_mp_nssl_2mom
  implicit none
  
  public nssl_2mom_driver
  public nssl_2mom_init
  private gamma_sp,gamxinf,GAML02, GAML02d300, GAML02d500, fqvs, fqis
  private gamma_dp, gamxinfdp, gamma_dpr
  private delbk, delabk
  private gammadp
  
  logical, private :: cleardiag = .false.
  private

#if ( WRF_CHEM == 1 )
  integer, parameter :: wrfchem_flag = 1
#else
  integer, parameter :: wrfchem_flag = 0
#endif

   logical, private:: is_aerosol_aware = .false.
      
      logical, private :: turn_on_cin = .false.
  
  integer, private :: myrankcm1 = 0
  integer, private :: eqtset = 1 ! Flag for use with cm1 to use alternate equation set (changes latent heating rates)
                                 ! value of > 2 invokes the equivalent version of eqtset=2 that applies updates to both theta and Pi.
   double precision, parameter, public :: zscale = 1.0d0 ! 1.000e-10
   double precision, parameter, public :: zscaleinv = 1.0d0/zscale ! 1.000e-10

  
  real, parameter :: warmonly = 0.0 ! testing parameter, set to 1.0 to reduce to warm-rain physics (ice variables stay zero)
  
  logical, parameter :: lwsm6 = .false. ! act like wsm6 for some single moment interactions

! some constants from WSM6
  real, parameter  :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
  real, parameter  :: roqimax = 2.08e22*dimax**8
  
! Params for dbz:
  integer  :: iuseferrier = 1  ! =1: use dry graupel only from Ferrier 1994; = 0: use Smith (wet graupel)
  integer  :: idbzci      = 1
  integer  :: iusewetgraupel = 1 ! =1 to turn on use of QHW for graupel reflectivity (only for ZVDM -- mixedphase)
                                 ! =2 turn on for graupel density less than 300. only 
  integer  :: iusewethail = 0 ! =1 to turn on use of QHW for graupel reflectivity (only for ZVDM -- mixedphase)
  integer  :: iusewetsnow = 0 ! =1 to turn on diagnosed bright band; =2 'old' snow reflectivity (dry), =3 'old' snow dbz + brightband
! microphysics

  real, public :: rho_qr = 1000., cnor = 8.0e5  ! cnor is set in namelist!!  rain params
  real, public :: rho_qs =  100., cnos = 3.0e6  ! set in namelist!!  snow params
  real, public :: rho_qh =  500., cnoh = 4.0e5  ! set in namelist!!  graupel params
  real, public :: rho_qhl=  800., cnohl = 4.0e4 ! set in namelist!!  hail params

  real, public :: hdnmn  = 170.0  ! minimum graupel density (for variable density graupel)
  real, public :: hldnmn = 500.0  ! minimum hail density (for variable density hail)

  real :: cnohmn  = 1.e-2 ! minimum intercept for 2-moment graupel (alphah < 0.5)
  real :: cnohlmn = 1.e-2 ! minimum intercept for 2-moment hail (alphahl < 0.5)
  
! Autoconversion parameters

  real   , private :: qcmincwrn      = 2.0e-3    ! qc threshold for autonconversion (LFO; for 10ICE use qminrncw for ircnw != 5)
  real   , private :: cwdiap         = 20.0e-6   ! threshold diameter of cloud drops (Ferrier 1994 autoconversion)
  real   , private :: cwdisp         = 0.15      ! assume droplet dispersion parameter (can be 0.3 for maritime)
  real   , public :: ccn            = 0.6e+09   ! set in namelist!! Central plains CCN value
  real   , public :: ccnuf          = 0        ! set in namelist!! Central plains CCN value
  real   , public  :: qccn, qccnuf               ! ccn "mixing ratio"
  real   , private :: old_qccn = -1.0
  integer, private :: iauttim        = 1         ! 10-ice rain delay flag
  real   , private :: auttim         = 300.      ! 10-ice rain delay time
  real   , private :: qcwmntim       = 1.0e-5    ! 10-ice rain delay min qc for time accrual

#if (NMM_CORE == 1)
! NMM WRF core does not have special boundary conditions for CCN, therefore set invertccn to true
      logical, parameter :: invertccn = .true. ! =true for base state of ccn=0, =false for ccn initialized in the base state
#else
      logical, private :: invertccn = .false. ! =true for base state of ccn=0, =false for ccn initialized in the base state
#endif
  logical :: switchccn = .false.
  real    :: old_cccn = -1.0
  logical :: restoreccn = .true. ! whether or not to nudge CCN back to base state (qccn) (only applies if CCNA is NOT predicted)
  real    :: ccntimeconst = 3600.  ! time constant for CCN restore (either for CCNA or when restoreccn = true)
  real, public  :: restoreccnfrac = 1.0  ! fraction of evaporated droplets that restore CCN
  real    :: ufccntimeconst = 6.*3600.  ! time constant for UFCCN decay (Blossey et al. 2018)
  real    :: ufbackground = 0.1e9       ! background ccnuf value (Blossey et al.)
  logical :: decayufccn = .false.
  integer :: i_uf_or_ccn = 0      ! 0 = ship adds UF; 1 = treat UF as regular ccn (add to qccn)

! sedimentation flags
! itfall -> 0 = 1st order fallout (other options removed)
! iscfall, infall -> fallout options for charge and number concentration, respectively
!                    1 = mass-weighted fall speed; 2 = number-weighted fallspeed.
  integer, private :: itfall = 0
  integer, private :: iscfall = 1
  integer, public :: irfall = -1
  integer, public :: isfall =  2 ! default limit with method II (more restrictive)
  logical, public :: do_accurate_sedimentation = .false. ! if true, recalculate fall speeds on sub time steps; (more expensive)
                                                         ! if false, reuse fall speeds on multiple steps (can have a noticeable speedup)
                                                         ! Mainly is an issue for small dz near the surface. 
  integer, public :: interval_sedi_vt = 1 ! interval for recalculating Vt in sedimentation subloop (only when do_accurate_sedimentation = .true.)
  integer, public :: infall = 4   ! 0 -> uses number-wgt for N; NO correction applied (results in excessive size sorting)
                          ! 1 -> uses mass-weighted fallspeed for N ALWAYS
                          ! 2 -> uses number-wgt for N and mass-weighted correction for N (Method II in Mansell, 2010 JAS)
                          ! 3 -> uses number-wgt for N and Z-weighted correction for N (Method I in Mansell, 2010 JAS)
                          ! 4 -> Hybrid of 2 and 3: Uses minimum N from each method (z-wgt and m-wgt corrections) (Method I+II in Mansell, 2010 JAS)
                          ! 5 -> uses number-wgt for N and uses average of N-wgt and q-wgt instead of Max.
  integer :: imydiagalpha = 0 ! apply MY diagnostic shape parameter for fall speeds (1=for fall speed only; 2=also for microphysics rates)
  real, public    :: rainfallfac = 1.0 ! factor to adjust rain fall speed (single moment only)
  real, public    :: icefallfac = 1.0 ! factor to adjust ice fall speed
  real, public    :: snowfallfac = 1.0 ! factor to adjust snow fall speed
  real, public    :: graupelfallfac = 1.0 ! factor to adjust graupel fall speed
  real, public    :: hailfallfac = 1.0 ! factor to adjust hail fall speed
  integer, public :: icefallopt = 3 ! 1= default, 2 = Ferrier ice fall speed; 3 = adjusted Ferrier (slightly high Vt)
  integer, public :: icdx = 6 ! (graupel) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.
                               ! 6= Milbrandt and Morrison (2013) density-based fall speed
  integer, public :: icdxhl = 6 ! (hail) 0=Ferrier; 1=leave drag coef. cd fixed; 2=vary by density, 4=set by user with cdxmin,cdxmax,etc.
                               ! 6= Milbrandt and Morrison (2013) density-based fall speed
  real    :: axh = 75.7149, bxh = 0.5
  real    :: axf = 75.7149, bxf = 0.5
  real    :: axhl = 206.984, bxhl = 0.6384
  real   , private :: cdhmin = 0.45, cdhmax = 0.8        ! defaults for graupel (icdx=4)
  real   , private :: cdhdnmin = 500., cdhdnmax = 800.0  ! defaults for graupel (icdx=4)
  real   , private :: cdhlmin = 0.45, cdhlmax = 0.6      ! defaults for hail (icdx=4)
  real   , private :: cdhldnmin = 500., cdhldnmax = 800.0  ! defaults for hail (icdx=4)
  real   , public :: vtmaxsed = 70. ! Limit on fall speed (m/s, all moments) for sedimentation calculations. Not applied to fall speeds for microphysical rates
  
  integer :: rssflg = 1   ! Rain size-sorting allowed (1, default), or disallowed (0).  if 0, sets N and Z-weighted fall speeds to q-weighted value
  integer :: sssflg = 1   ! As above but for snow
  integer :: hssflg = 1   ! As above but for graupel
  integer :: hlssflg = 1  ! As above but for hail

! input flags

  integer, public :: ndebug = -1, ncdebug = 0
  integer, private :: ipconc = 5
  integer, public :: inucopt = 0
  integer, private :: ichaff = 0
  integer, parameter :: ilimit = 0
  
  real, private :: constccw = -1.

  real, private :: cimn = 1.0e3, cimx = 1.0e6

  real   , public :: rhofrz = 900 ! density of freezing drops
  real   , public :: ifrzg = 1.0 ! fraction of frozen drops (Bigg freezing) going to graupel. 1=freeze all rain to graupel, 0=freeze all to hail
  real   , public :: ifiacrg = 1.0 ! fraction of frozen drops (3-component freezing qiacr) going to graupel. 1=freeze all rain to graupel, 0=freeze all to hail
  real   , public :: ifrzs = 1.0 ! fraction of small frozen drops going to snow. 1=freeze rain to snow, 0=freeze to cloud ice
  real   , public :: ffrzs = 0.0 ! fraction of other initiated cloud ice going to snow. 1=freeze rain to snow, 0=freeze to cloud ice
  real   , public :: f2h = 1.0 ! fraction of cloud ice conversion going to graupel (vs. frozen drops). For testing
  integer, private :: irwfrz = 1 ! compute total rain that can freeze (checks heat budget)
  integer, private :: irimtim = 0 ! future use
!  integer, private :: infdo = 1   ! 1 = calculate number-weighted fall speeds

  integer, public :: irimdenopt = 1 ! = 1 for default Heymsfield and Pflaum (1985); = 2 for experimental Cober and List (1993); = 3 Macklin
  real   , private :: rimc1 = 300.0, rimc2 = 0.44  ! rime density coeff. and power (Default Heymsfield and Pflaum, 1985)
  real   , private :: rimc3 = 170.0                ! minimum rime density
  real    :: rimc4 = 900.0                ! maximum rime density
  real   , private :: rimtim = 120.0               ! cut-off rime time (10ICE)
  real   , private :: eqtot = 1.0e-9               ! threshold for mass budget reporting
  real, public :: rimdenvwgt = 0.0 ! weight (0-1) given to number-weighted fall speed when calculating rime density

  integer, private :: ireadmic = 0

  integer, public :: idiagnosecnu = 0 ! =1 to diagnose cnu based on Chandrakar et al. 2016 data; =2 for Geoffroy et al. (2010, ACP)
  integer, private :: iccwflg = 1     ! sets max size of first droplets in parcel to 4 micron radius (in two-moment liquid)
                             ! (first nucleation is done with a KW sat. adj. step)
  integer, private :: issfilt = 0     ! flag to turn on filtering of supersaturation field
  integer, public :: icnuclimit = 0  ! limit droplet nucleation based on Konwar et al. (2012) and Chandrakar et al. (2016)
  integer, public :: irenuc = 7      ! =1 to always allow renucleation of droplets within the cloud (do no use, obsolete)
                                      ! =2 renucleation following Twomey/Cohard&Pinty
                                      ! =7 New renucleation that requires prediction of the number of activated nuclei
                             ! i.e., not only at cloud base
  integer, private :: irenuc3d = 0      ! =1 to include horizontal gradient in renucleation of droplets within the cloud
  real    :: renucfrac = 0.0 ! = 0 : cnuc = cwccn
                             ! = 1 : cnuc = actual available CCN
                             ! otherwise cnuc = cwccn*(1. - renufrac) + ccnc(1:ngscnt)*renucfrac
  real    :: ssf2kmax = 10. ! max value for ssf**cck in irenuc=4 or 5
  real   , private :: cck = 0.6       ! exponent in Twomey expression
  real   , private :: ciintmx = 1.0e6 ! limit on ice concentration from primary nucleation

  real   , private :: cwccn ! , cwmasn,cwmasx
  real   , private :: ccwmx

  integer, private :: idocw = 1, idorw = 1, idoci = 1, idoir = 1, idoip = 1, idosw = 1
  integer, private :: idogl = 1, idogm = 1, idogh = 1, idofw = 1, idohw = 1, idohl = 1
!  integer, private :: ido(3:14) = / 12*1 /


! 0,2, 5.00e-10, 1, 0, 0, 0      : itype1,itype2,cimas0,icfn,ihrn,ibfc,iacr
  integer, public :: itype1 = 0, itype2 = 2  ! controls Hallett-Mossop process
  integer, private :: in_freeze_rain_first = 0 ! =1 use in to freezed rain drops (if none, then freeze droplets)
  integer, public :: icenucopt = 1       ! =1 Meyers/Ferrier primary ice nucleation; =2 Thompson/Cooper, =3 Phillips (Meyers/Demott), =4 DeMott (2010)
  real, public :: naer = 1.0e6  ! background large aerosol conc. for DeMott
  integer, private :: icfn = 2                ! contact freezing: 0 = off; 1 = hack (ok for single moment); 2 = full Cotton/Meyers version
  integer, private :: ihrn = 0            ! Hobbs-Rangno ice multiplication (Ferrier, 1994; use in 10-ice only)
  integer, public :: ibfc = 1            ! Flag to use Bigg freezing on droplets (0 = off (uses alternate freezing), 1 = on)
  real, private  :: cwfrz2snowfrac = 0.0 ! fraction of freezing droplet mass to send to snow
  real, private  :: cwfrz2snowratio = 5. ! Assumed number of frozen droplets in a cluster
  integer, public :: iremoveqwfrz = 1    ! Whether to remove (=1) or not (=0) the newly-frozen cloud droplets (ibfc=1) from the CWC used for charge separation
  integer, public :: iacr = 2            ! Flag for drop contact freezing with crytals
                                 ! (0=off; 1=drops > 500micron diameter; 2 = > 300micron)
  integer, public :: icrcev = 1          ! 1 = old crcev; 2 = crcev scaled by vtrain ratio (num/mass); 3 = set to zero
  integer, public :: icracr = 1          ! Flag to turn rain self-collection on/off (=0 to turn off)
  integer, public :: icracrthresh = 1    ! For rain self-coll. thresh. use: 1 = mean diam of 2mm; 2 = rain median volume diam of 1.9mm
  integer, private :: ibfr = 2            ! Flag for Bigg freezing conversion of freezing drops to graupel
                                 ! (1=min graupel size is vr1mm; 2=use min size of dfrz, 5= as for 2 and apply dbz conservation)
  integer, public :: ibiggopt = 2        ! 1 = old Bigg; 2 = experimental Bigg (only for imurain = 1, however)
  integer :: ibiggsmallrain = 0  ! 1 = When rain is too small, freeze none to graupel and send all to snow (experimental)
  integer, public :: iacrsize = 5        ! assumed min size of drops freezing by capture
                                 !  1: > 500 micron diam
                                 !  2: > 300 micron
                                 !  3: > 40 micron
                                 !  4: all sizes
                                 !  5: > 150 micron (only for imurain = 1)
  real   , private :: cimas0 = 6.62e-11   ! default mass of Hallett-Mossop crystals
                                 ! 6.62e-11kg results in half the diam. (60 microns) of old default value of 5.0e-10
  real   , private :: cimas1 = 6.88e-13   ! default mass of new ice crystals
  real   , private :: splintermass = 6.88e-13
  real   , private :: cfnfac = 0.1        ! Hack factor that goes with icfn=1
  integer, private :: iscni = 4           ! default option for ice crystal aggregation/conversion to snow
  real   , public :: fscni = 1.0         ! factor for calculating cscni
  logical, private :: imeyers5 = .false.  ! .false.=off, true=on for Meyers ice nucleation for temp > -5 C
  real   , private :: dmincw = 15.0e-6    ! minimum droplet diameter for collection for iehw=3
  integer, public :: iehw = 1            ! 0 -> ehw=ehw0; 1 -> old ehw; 2 -> test ehw with Mason table data
  integer, public :: iefw = 1            ! 0 -> ehw=ehw0; 1 -> old ehw; 2 -> test ehw with Mason table data
  integer, public :: iehlw = 1           ! 0 -> ehlw=ehlw0; 1 -> old ehlw; 2 -> test ehlw with Mason table data
                                 ! For ehw/ehlw = 1, ehw0/ehlw0 act as maximum limit on collection efficiency (defaults are 1.0)
  integer, private :: ierw = 1            ! for single-moment rain (LFO/Z)
  integer, public :: iehr0c = 0          ! 0 -> no collection for T > 0C;  1 -> turn on collection/shedding for T > 0C
  integer, public :: iehlr0c = 0         ! 0 -> no collection for T > 0C;  1 -> turn on collection/shedding for T > 0C
  real   , public :: ehw0 = 0.9 ! 0.5          ! constant or max assumed graupel-droplet collection efficiency
  real   , private :: erw0 = 1.0          ! constant assumed rain-droplet collection efficiency
  real   , public :: ehlw0 = 0.9 ! 0.75        ! constant or max assumed hail-droplet collection efficiency
  real   , public :: efw0 = 0.5          ! constant or max assumed graupel-droplet collection efficiency
  real    :: ehr0 = 1.0          ! constant or max assumed graupel-rain collection efficiency
  real    :: efr0 = 1.0          ! constant or max assumed graupel-rain collection efficiency
  real    :: ehlr0 = 1.0         ! constant or max assumed hail-rain collection efficiency
  real   , public :: exwmindiam = 0.0    ! minimum diameter of droplets for riming. if set > 0, will exclude that fraction of mass/number from accretion (idea from Furtado and Field 2017 JAS but also Fierro and Mansell 2017)
  

  real   , private :: esilfo0 = 1.0       ! factor for LFO collection efficiency of snow for cloud ice.
  real   , private :: ehslfo0 = 1.0       ! factor for LFO collection efficiency of hail/graupel for snow.

  integer, private :: ircnw    = 5        ! single-moment warm-rain autoconversion option.  5= Ferrier 1994.
  real   , private :: qminrncw = 2.0e-3   ! qc threshold for rain autoconversion (NA for ircnw=5)

  integer, private :: iqcinit = 2         ! For ZVDxx schemes, flag to choose which way to initialize droplets
                                 ! 1 = Soong-Ogura adjustment
                                 ! 2 = Saturation adjustment to value of ssmxinit
                                 ! 3 = KW adjustment

  real   , private :: ssmxinit = 0.4      ! saturation percentage to adjust down to for initial cloud
                                 ! formation (ZVDxx scheme only)

  real   , private :: ewfac = 1.0         ! hack factor applied to graupel and hail collection eff. for droplets
  real   , private :: eii0 = 0.1 ,eii1 = 0.1  ! graupel-crystal coll. eff. parameters: eii0*exp(eii1*min(temcg(mgs),0.0))
                                     ! set eii1 = 0 to get a constant value of eii0
  real   , public :: eii0hl = 0.2 ,eii1hl = 0.0  ! hail-crystal coll. eff. parameters: eii0hl*exp(eii1hl*min(temcg(mgs),0.0))
                                     ! set eii1hl = 0 to get a constant value of eii0hl
  real, public :: ewi_dcmin = 15.0e-06 ! minimum droplet diameter for nonzero ewi
  real, public :: ewi_dimin = 30.0e-06 ! minimum ice crystal diameter for nonzero ewi
  real   , public :: eri0 = 0.1   ! rain efficiency to collect ice crystals
  real   , public :: eri_cimin = 10.e-6      ! minimum ice crystal diameter for collection by rain
  real   , public :: esi0 = 0.1              ! linear factor in snow-ice collection efficiency
  real   , private :: ehs0 = 0.1, ehs1 = 0.1  ! graupel-snow coll. eff. parameters: ehs0*exp(ehs1*min(temcg(mgs),0.0))
                                     ! set ehs1 = 0 to get a constant value of ehs0
  integer :: iessopt = 1  ! 1 = Original (no factor); 2 = factor based on wvel; 3 = factor based on SSI
                          ! 4 = as 3 but sets min factor of 0.1 and goes to full value at 0.5% SSI
  real   , private :: ess0 = 0.5, ess1 = 0.05 ! snow aggregation coefficients: ess0*exp(ess1*min(temcg(mgs),0.0))
                                     ! set ess1 = 0 to get a constant value of ess0
  real   , public :: esstem1 = -25.  ! lower temperature where snow aggregation turns on
  real   , public :: esstem2 = -20.  ! higher temperature for linear ramp of ess from zero at esstem1 to formula value at esstem2
  real   , public :: essrmax = 0.02  ! maximum snow radius (meters) for csacs
  real   , public :: essfrac1 = 0.5  ! snow mass fraction 1 for aggregation roll-off
  real   , public :: essfrac2 = 0.75 ! snow mass fraction 2 for aggregation roll-off
  integer, public :: iessec0flag = 0 ! flag to activate aggregation roll-off
  real   , private :: ehsfrac = 1.0           ! multiplier for graupel collection efficiency in wet growth
  real   , private :: ehimin = 0.0 ! Minimum collection efficiency (graupel - ice crystal)
  real   , private :: ehimax = 1.0 ! Maximum collection efficiency (graupel - ice crystal)
  real   , private :: ehsmax = 0.5 ! Maximum collection efficiency (graupel - snow)
  real   , private :: ecollmx = 0.5 ! Maximum collision efficiency for graup/hail with ice; used only for charging rates
  integer, public :: iglcnvi = 1  ! flag for riming conversion from cloud ice to rimed ice/graupel
  integer, public :: iglcnvs = 2  ! flag for conversion from snow to rimed ice/graupel

  real   , private :: rz          ! reflectivity conservation factor for graupel/rain
                         ! now calculated in icezvd_dr.F from alphah and rnu
                         ! currently only used for graupel melting to rain
  real   , private :: rzhl        ! reflectivity conservation factor for hail/rain
                         ! now calculated in icezvd_dr.F from alphahl and rnu

  real   , private :: rzs     ! reflectivity conservation factor for snow(imusnow=3) with rain (imurain=1)

  real   , private :: alphahacx = 0.0 ! assumed minimum shape parameter for zhacw and zhacr

  real   , private :: fconv = 1.0  ! factor to boost max graupel depletion by riming conversions in 10ICE

  real   , private :: rg0 = 400.0  ! reference graupel density for graupel fall speed

  integer, private :: rcond = 2    ! (Z only) rcond = 2 includes rain condensation in loop with droplet condensation
                                   ! 0 = no condensation on rain; 1 = bulk condensation on rain
  integer, parameter, private :: icond = 1    ! (Z only) icond = 1 calculates ice deposition (crystals and snow) BEFORE droplet condensation
                          ! icond = 2 does not work (intended to calc. dep in loop with droplet cond.)
  integer, public :: iqis0 = 2    ! = 1 for normal qis; = 2 to set qis to use T = 0C when T > 0C  
  
  real   , public :: dfrz = 0.15e-3 ! 0.25e-3  ! minimum diameter of frozen drops from Bigg freezing (used for vfrz) for iacr > 1
                            ! and for ciacrf for iacr=4
  real   , public :: dmlt = 3.0e-3  ! maximum diameter for rain melting from graupel and hail
  real   , public :: dshd = 1.0e-3  ! nominal diameter for rain drops shed from graupel/hail
  integer, public :: ivshdgs   = 1  ! 0 = use 1mm for all shedding (non-mixedphase); 1 = use vshdgs with sheddiam
  integer, private :: ished2cld = 0  ! 1: Send shed liquid (from wet growth) to cloud droplets

  integer, public :: ihmlt = 2      ! 1=old melting with vmlt; 2=new melting using mean volume diam of graupel/hail
  integer, public :: imltshddmr = 2 ! 0 (default)=mean diameter of drops produced during melting+shedding as before (using mean diameter of graupel/hail
                            ! and max mean diameter of rain)
                            ! 1=new method where mean diameter of rain during melting is adjusted linearly downward 
                            ! toward 3 mm for large (> sheddiam) graupel and hail, to take into account shedding of 
                            ! smaller drops.  sheddiam0 controls the size of graupel/hail above which the assumed 
                            ! mean diameter of rain is set to 3 mm
                            ! Only valid for ihmlt = 2 for ZVD(H) but also applies to ZVD(H)M
                            ! 2 = method that sets the resulting rain size ( vshdgs ) according to the mass-weighted diameter of the ice

   real  :: mltdiam1 = 9.0e-3, mltdiam2 = 16.0e-3, mltdiam3 = 19.0e-3, mltdiam4 = 200.0e-3, mltdiam05 = 4.5e-3

  integer, public :: nsplinter = 0  ! number of ice splinters per freezing drop, if negative, then per resulting graupel particle
  real,    private :: lawson_splinter_fac = 2.5e-11  ! constant in Lawson et al. (2015, JAS) for ice particle production from freezing drops
  integer, private :: isnwfrac = 0   ! 0= no snow fragmentation; 1 = turn on snow fragmentation (Schuur, 2000)

!  integer, private :: denscale = 1  ! 1=scale num. conc. and charge by air density for advection, 0=turn off for comparison

  real, public  :: qhdpvdn = -1.
  real, public  :: qhacidn = -1.

  integer, private :: iraintypes = 0
  logical, private :: mixedphase = .false.   ! .false.=off, true=on to include mixed phase graupel
  integer, private :: imixedphase = 0
  logical, private :: qsdenmod = .false.     ! true = modify snow density by linear interpolation of snow and rain density
  logical, private :: qhdenmod = .false.     ! true = modify graupel density by linear interpolation of graupel and rain density
  logical, private :: qsvtmod = .false.      ! true = modify snow fall speed by linear interpolation of snow and rain vt
  real   , private :: sheddiam   = 8.0e-03  ! minimum diameter of graupel before shedding occurs
  real    :: sheddiamlg = 10.0e-03  ! diameter of hail to use fwmlarge
  real    :: sheddiam0  = 20.0e-03  ! diameter of hail at which all water is shed
  
  integer :: ifwmhopt = 2 ! option for calculating maximum liquid fraction when fwmh and/or fwmhl is set to -1
                          ! 1 = maximum based on size of maximum mass diameter
                          ! 2 = integrate over spectrum for maximum liquid (experimental)

  integer :: ihxw2rain = 0 ! = 0 no transfer
                           ! = 1 transfer completely melted (99.5%) graupel/hail to rain when fwmh/fwmhl is set to -1.

  real   , private :: fwms = 0.5 ! maximum liquid water fraction on snow
  real   , private :: fwmh = 0.5 ! maximum liquid water fraction on graupel
  real   , private :: fwmhl = 0.5 ! maximum liquid water fraction on hail
  real    :: fwmlarge = 0.2 ! maximum liquid water fraction on hail larger than sheddiam
  integer :: ifwmfall = 0   ! whether to interpolate toward rain fall speed for graupel and hail
                            ! when diam < sheddiam and liquid fraction is predicted (0=no, 1=yes)
  
  logical :: rescale_high_alpha = .false.  ! whether to rescale number. conc. when alpha = alphamax (3-moment only)
  logical :: rescale_low_alpha = .true.    ! whether to rescale Z (graupel/hail) when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphar = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphah = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)
  logical :: rescale_low_alphahl = .true.    ! whether to rescale Z for rain when alpha = alphamin (3-moment only)

  real, parameter :: alpharmax = 8. ! limited for rwvent calculation
  
  integer, public ::  ihlcnh = -1  ! which graupel -> hail conversion to use
                          ! 1 = Milbrandt and Yau (2005) using Ziegler 1985 wet growth diameter
                          ! 2 = Straka and Mansell (2005) conversion using size threshold
                          ! 3 = Conversion using wet growth diameter
  real, public :: hlcnhdia = 1.e-3 ! threshold diameter for graupel -> hail conversion for ihlcnh = 1 option.
  real, public :: hlcnhqmin = 0.1e-3 ! minimum graupel mass content for graupel -> hail conversion (ihlcnh = 1)
  real   , private :: hldia1 = 10.0e-3  ! threshold diameter for graupel -> hail conversion for ihlcnh = 2 option.
  integer, public  :: incwet = 0    ! flag to do wet growth only on D > D_wet
  integer, public  :: iusedw = 0    ! flag to use experimental wet growth ice diameter for gr -> hl conversion (=1 turns on)
  real   , private  :: dwmin   = 5.0e-3  ! Minimum diameter with iusedw (can stay at 0 or be set to something larger)
  real   , private  :: dwetmin = 5.0e-3  ! Minimum diameter with iusedw (can stay at 0 or be set to something larger)
  real   , public  :: dwmax  = 15.e-3 ! for ihlcnh, always convert this size and larger whether or not there is wet growth
  real   , public  :: dwtempmin = 242. ! lowest temperature to allow wet growth conversion to hail
  real   , public  :: dwehwmin = 0.   ! Minimum ehw to use to find wet growth diameter (if > ehw0, then wet growth diam becomes smaller)
  real   , public  :: dg0thresh = 0.15 ! graupel wet growth diameter above which we say do not bother
  integer :: ifddenfac = 0  ! = 1 to use density threshold to count FD as GR when converting to HL
  real    :: fddenthresh = 500. ! if ifddenfac > 0, then hail from FD with lower density are considered to come from graupel
  integer :: icvhl2h = 0   ! allow conversion of hail back to graupel when hail density gets close to minimum allowed

  integer, public :: imurain = 1 ! 3 for gamma-volume, 1 for gamma-diameter DSD for rain.
  integer, public :: imusnow = 3 ! 3 for gamma-volume, 1 for gamma-diameter DSD for snow (=1 NOT IMPLEMENTED!!).
  integer, public :: iturbenhance = 0 ! warm-rain collision enhancement
                              ! 1 = enhance autoconversion only
                              ! 2 = add rain collection of cloud
                              ! 3 = add rain self-collection
  integer, public :: isedonly = 0 ! 1= only do sedimentation and skip other microphysics
  integer, public :: iferwisventr = 2 ! =1 for Ferrier rwvent, =2 for Wisner rwvent (imurain=1)
  integer, public :: izwisventr   = 2 ! =1 for old Ziegler rwvent, =2 for Wisner-style rwvent (imurain=3)
  integer :: iresetmoments = 0 ! if >0, then set all moments to zero when one of them is zero (3-moment only)
  integer, public :: imaxdiaopt = 3 
                               ! = 1 use mean diameter for breakup
                               ! = 2 use maximum mass diameter for breakup
                               ! = 3 use mass-weighted diameter for breakup
  integer :: iraintailbreak = 0 ! 1 = on
  real    :: draintail      = 8.e-3 ! starting size for rain breakup
  integer, public :: dmrauto       = 0 
                              ! = -1 no limiter on crcnw
                              ! =  0 limit crcnw when qr > 1.2*L (Cohard-Pinty 2002)
                              ! =  1 DTD version based on MY code
                              ! =  2 DTD mass-weighted version based on MY code
                              ! =  3 Milbrandt version (from Cohard and Pinty code
  integer :: dmropt = 0 ! extra option for crcnw
  integer :: dmhlopt = 0 ! options for graupel -> hail conversion
  integer :: irescalerainopt = 3 ! 0 = default option
                                 ! 1 = qx(mgs,lc) > qxmin(lc) 
                                 ! 2 = qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < 3.0
                                 ! 3 = temcg(mgs) > 0.0.and. qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < 3.0 
  real    :: rescale_wthresh = 3.0
  real    :: rescale_tempthresh = 0.0
  real, parameter :: alpharaut = 0.0 ! MY2005 for autoconversion
  real    :: cxmin = 1.e-8  ! threshold cutoff for number concentration
  real    :: zxmin = 1.e-28 ! threshold cutoff for reflectivity moment
  
  integer :: ithompsoncnoh = 0 ! For single moment graupel only
                           ! 0 = fixed intercept
                           ! 1 = intercept based on graupel mass

  integer :: ivhmltsoak = 1   ! 0=off, 1=on : flag to simulate soaking (graupel/hail) during melting 
                         ! when liquid fraction is not predicted
  logical, public :: iwetsoak = .true. ! soak and freeze during wet growth or not
  integer, public :: ioldlimiter = 0 ! test switch for new(=0) or old(=1) size limiter at the end of GS for 3-moment categories
  integer, public :: isnowfall = 2   ! Option for choosing between snow fall speed parameters
                         ! 1 = original Zrnic et al. (Mansell et al. 2010)
                         ! 2 = Ferrier 1994 (results in slower fall speeds)

  integer, public :: isnowdens = 1   ! Option for choosing between snow density options
                             ! 1 = constant of 100 kg m^-3
                             ! 2 = Option based on Cox 
  
  integer, public  :: ibiggsnow   = 3 ! 1 = switch conversion over to snow for small frozen drops from Bigg freezing
                                       ! 2 = switch conversion over to snow for small frozen drops from rain-ice interaction
                                       ! 3 = switch conversion over to snow for small frozen drops from both
  real    :: biggsnowdiam = -1.0 ! if >0, use for ibiggsnow threshold
  
  integer, public  :: ixtaltype = 1 ! =1 column, =2 disk (similar to Takahashi)

  real, public  :: takshedsize1 = 0.15 ! diameter (cm) of drop shed from ice with D > 1.9 cm
  real, public  :: takshedsize2 = 0.3 ! diameter (cm) of drop shed from ice with D < 1.9 cm and D > 0.8 cm
  real, public  :: takshedsize3 = 0.45 ! diameter (cm) of drop shed from ice with D < 1.6 cm and D > 0.8 cm
  integer, public  :: numshedregimes = 3
  
  real, public     :: evapfac     = 1.0 ! Multiplier on rain evaporation rate
  real, public     :: depfac      = 1.0 ! Multiplier on graupel/hail deposition/sublimation rate
  real,private,parameter :: meltfac     = 1.0 ! Multiplier on graupel/hail melting rate

  integer, public :: ibinhmlr = 0  ! =1 use incomplete gammas to determine melting from larger and smaller sizes of graupel, and appropriate shed drop sizes 
                           ! =2 to test melting by temporary bins
  integer, public :: ibinhlmlr = 0  ! =1 use incomplete gammas to determine melting from larger and smaller sizes of hail, and appropriate shed drop sizes 
                            ! =2 to test melting by temporary bins
  integer, public :: ibinnum   = 2  ! number of bins for melting of smaller ice (for ibinhmlr = 1)
  integer, public :: iqhacrmlr = 1  ! turn on/off qhacrmlr
  integer, public :: iqhlacrmlr = 1  ! turn on/off qhlacrmlr
  integer, public :: iqhacwshr = 1  ! turn on/off qhacw for T > 0
  integer, public :: iqhlacwshr = 1  ! turn on/off qhlacw for T > 0
  real, public :: binmlrmxdia = 40.e-3 ! threshold diameter (graupel/hail) to switch bin-bulk melting to use standard chmlr
  real, public :: binmlrzrrfac = 1.0 ! factor for reflectivity change ice that sheds while melting
  real, public :: snowmeltdia = 0 ! if nonzero, sets the size of rain drops from melting snow.
  real, public :: alphasmlr0 = 14.0 ! shape parameter for drops formed from melting/shedding snow
  real, public :: delta_alphamlr = 0.5 ! offset from alphamax at which melting does not further collapse the shape parameter
  
  integer :: iqvsopt = 0 ! =0 use old default for tabqvs; =1 use Bolton formulation (Rogers and Yau)

  integer :: imaxsupopt = 4 ! how to treat saturation adjustment in two-moment droplets
                            ! 1 = add droplets with same mean mass as current droplets
                            ! 2 = add droplets with minimum radius of 30 microns
                            ! 3 = only add 1.5*cxmin to number concentration (allow max size to apply)
                            ! 4 = add droplets with minimum radius of 20 microns
  real    :: maxsupersat = 1.9 ! maximum supersaturation ratio, above which a saturation adustment is done
  real    :: maxlowtempss = 1.08 ! Sat. ratio threshold for allowing droplet nucleation at T < tfrh
  real    :: ssmxuf = 4.0 ! supersaturation at which to start using "ultrafine" CCN (if ccnuf > 0.)
  

  integer, parameter :: icespheres = 0 ! turn ice spheres (frozen droplets) on (1) or off (0). NOT COMPLETE in WRF/ARPS/1 CODE!
  integer, parameter :: lqmx = 30
  integer, parameter :: lt = 1
  integer, parameter :: lv = 2
  integer, parameter :: lc = 3
  integer, parameter :: lr = 4
  integer, parameter :: li = 5
  integer, private :: lis = 0
  integer, private :: ls = 6
  integer, private :: lh = 7
  integer, private :: lf = 0
  integer, private :: lhl = 0

  integer, public  :: lccn = 9 ! 0 or 9, other indices adjusted accordingly
  integer, private :: lccnuf = 0
  integer, private :: lccna = 0
  integer, private :: lcina = 0
  integer, private :: lcin = 0
  integer, private :: lnc = 9
  integer, private :: lnr = 10
  integer, private :: lni = 11
  integer, private :: lnis = 0
  integer, private :: lns = 12
  integer, private :: lnh = 13
  integer, private :: lnf = 0
  integer, private :: lnhl = 0
  integer, private :: lnhf = 0
  integer, private :: lnhlf = 0
  integer, private :: lss = 0
  integer :: lvh = 15

  integer, private :: lhab = 8
  integer, private :: lg = 7

! Particle volume

  integer :: lvi = 0
  integer :: lvs = 0
  integer :: lvgl = 0
  integer :: lvgm = 0
  integer :: lvgh = 0
  integer :: lvf = 0
!  integer :: lvh = 16
  integer :: lvhl = 0

! liquid water fraction (not predicted here but tested for)
  integer :: lhw = 0
  integer :: lfw = 0
  integer :: lsw = 0
  integer :: lhlw = 0
  integer :: lhwlg = 0
  integer :: lhlwlg = 0

! reflectivity (6th moment) ! not predicted here but may be tested against

  integer :: lzr = 0
  integer :: lzi = 0
  integer :: lzs = 0
  integer :: lzgl = 0
  integer :: lzgm = 0
  integer :: lzgh = 0
  integer :: lzf = 0
  integer :: lzh = 0
  integer :: lzhl = 0

! Space charge

  integer :: lscw = 0
  integer :: lscr = 0
  integer :: lsci = 0
  integer :: lscis = 0
  integer :: lscs = 0
  integer :: lsch = 0
  integer :: lscf = 0
  integer :: lschl = 0
  integer :: lscwi = 0
  integer :: lscpi = 0
  integer :: lscni = 0
  integer :: lscpli = 0
  integer :: lscnli = 0
  integer :: lschab = 0

  integer :: lscb = 0
  integer :: lsce = 0
  integer :: lsceq = 0

!  integer, parameter :: lscmx = 100

  integer :: lne = 0 ! last varible for transforming

  real :: cnoh0 = 4.0e+5
  real :: hwdn1 = 700.0

  real    :: alphai  = 0.0 ! shape parameter for ZIEG ice crystals ! not currently used
  real    :: alphas  = 0.0 ! shape parameter for ZIEG snow         ! used only for single moment
  real    :: alphar  = 0.0 ! shape parameter for rain (imurain=1 only)
  real, public    :: alphah  = 0.0 ! set in namelist!! shape parameter for ZIEG graupel
  real, public    :: alphahl = 1.0 ! set in namelist!! shape parameter for ZIEG hail

  real    :: dmuh    = 1.0  ! power in exponential part (graupel)
  real    :: dmuhl   = 1.0  ! power in exponential part (hail)

  real, private   :: alphamax = 15.
  real, private   :: alphamin = 0.
  real, parameter :: rnumin = -0.8
  real, parameter :: rnumax = 15.0

  
  real            :: cnu = 0.0 ! default value of droplet shape parameter. Can be diagnosed by setting idiagnosecnu=1
  real, parameter :: rnu = -0.8, snu = -0.8, cinu = 0.0
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
  
  real xnu(lc:lqmx) ! 1st shape parameter (mass)
  real xmu(lc:lqmx) ! 2nd shape parameter (mass)
  real dnu(lc:lqmx) ! 1st shape parameter (diameter)
  real dmu(lc:lqmx) ! 2nd shape parameter (diameter)
  
  real ax(lc:lqmx)
  real bx(lc:lqmx)
  real fx(lc:lqmx)

      real da0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real dab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real dab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real da1 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real bb  (lc:lqmx)


! put ipelec here for now....
  integer :: ipelec = 0
  integer :: isaund = 0
  logical :: idoniconly = .false.
  integer, public :: elec_on_time = -1     ! time (seconds) to turn on charge separation.
  integer, public :: elec_ramp_time = 0   ! time (interval) for linear ramp after elec_on_time 
                                   ! (i.e., linear factor on chg sep to smoothly turn on elec)
                                   ! full charging rate is achieved at time = elec_on_time + elec_ramp_time
  integer :: jchgs = 3  ! number of points near boundary where charging is turned off (to keep lightning from getting wonky)
  integer :: jchgn = 2
  integer :: ichge = 3
  integer :: ichgw = 2
  real    :: charging_border = 4000. ! width of no-charging zone from boundary
      real, private    :: delqnw = -1.0e-10!-1.0e-12 !
      real, private    :: delqxw =  1.0e-10! 1.0e-12 !
      real :: tindmn = 233, tindmx = 298.0  ! min and max temperatures where inductive charging is allowed

!
!  gamma function lookup table
!
      integer ngm0,ngm1,ngm2
      parameter (ngm0=3001,ngm1=500,ngm2=500)
      double precision, parameter :: dgam = 0.01, dgami = 100.
      double precision gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)

      integer, parameter :: nqiacralpha =  300 !480 ! 240 ! 120 ! 15
      integer, parameter :: nqiacrratio =  400 ! 500 !50  ! 25
!      real,    parameter :: maxratiolu = 25.
      real,    parameter :: maxratiolu = 100. ! 25.
      real,    parameter :: maxalphalu = 15.
      real,    parameter :: minalphalu = -0.95
      real,    parameter :: dqiacralpha = maxalphalu/Float(nqiacralpha), dqiacrratio = maxratiolu/Float(nqiacrratio) 
      real,    parameter :: dqiacrratioinv = 1./dqiacrratio, dqiacralphainv = 1./dqiacralpha
      integer, parameter :: ialpstart = minalphalu*dqiacralphainv
      real :: ciacrratio(0:nqiacrratio,ialpstart:nqiacralpha)
      real :: qiacrratio(0:nqiacrratio,ialpstart:nqiacralpha)
      real :: ziacrratio(0:nqiacrratio,ialpstart:nqiacralpha)
      double precision :: gamxinflu(0:nqiacrratio,ialpstart:nqiacralpha,12,2) ! last index for graupel (1) or hail (2)
!      real :: ciacrratio(0:nqiacrratio,0:nqiacralpha)
!      real :: qiacrratio(0:nqiacrratio,0:nqiacralpha)
!      real :: ziacrratio(0:nqiacrratio,0:nqiacralpha)
!      double precision :: gamxinflu(0:nqiacrratio,0:nqiacralpha,12,2) ! last index for graupel (1) or hail (2)

! for 3-moment collection coefficients
      real, save :: dab0lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real, save :: dab1lu(ialpstart:nqiacralpha,ialpstart:nqiacralpha,lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005

    integer, parameter :: ngdnmm = 9
    real :: mmgraupvt(ngdnmm,3)  ! Milbrandt and Morrison (2013) fall speed coefficients for graupel/hail

    data mmgraupvt(:,1) / 50., 150., 250., 350., 450., 550., 650., 750., 850./
    data mmgraupvt(:,2) / 62.923, 94.122, 114.74, 131.21, 145.26, 157.71, 168.98, 179.36, 189.02 /
    data mmgraupvt(:,3) / 0.67819, 0.63789, 0.62197, 0.61240, 0.60572, 0.60066, 0.59663, 0.59330, 0.59048 /

      integer lsc(lc:lqmx)
      integer ln(lc:lqmx)
      integer ipc(lc:lqmx)
      integer lvol(lc:lqmx)
      integer lz(lc:lqmx)
      integer lliq(li:lqmx)
      integer denscale(lc:lqmx) ! flag for density scaling (mixing ratio conversion)

      integer ido(lc:lqmx)
      logical ldovol

      real xdn0(lc:lqmx)
      real xdnmx(lc:lqmx), xdnmn(lc:lqmx)
      real cdx(lc:lqmx)
      real cno(lc:lqmx)
      real xvmn(lc:lqmx), xvmx(lc:lqmx)
      real qxmin(lc:lqmx)
      real qxmin_init(lc:lqmx)

      integer nqsat
      parameter (nqsat=1000001) ! (nqsat=20001)
      real fqsat,fqsati
      parameter (fqsat=0.002,fqsati=1./fqsat)
      real tabqvs(nqsat),tabqis(nqsat),dtabqvs(nqsat),dtabqis(nqsat)

!
!  constants
!
      real, parameter :: ar = 841.99666         ! rain terminal velocity power law coefficient (LFO)
      real, parameter :: br = 0.8               ! rain terminal velocity power law coefficient (LFO)
      real, parameter :: aradcw = -0.27544      !
      real, parameter :: bradcw = 0.26249e+06   !
      real, parameter :: cradcw = -1.8896e+10   !
      real, parameter :: dradcw = 4.4626e+14    !
      real, parameter :: bta1 = 0.6             ! beta-1 constant used for ice nucleation by deposition (Ferrier 94, among others)
      real, parameter :: cnit = 1.0e-02         ! No for ice nucleation by deposition (Cotton et al. 86)
      real, parameter :: dragh = 0.60           ! coefficient used to adjust fall speed for hail versus graupel (Pruppacher and Klett 78)
      real, parameter :: dnz00 = 1.225          ! reference/MSL air density
      real, parameter :: rho00 = 1.225          ! reference/MSL air density
!      cs = 4.83607122       ! snow terminal velocity power law coefficient (LFO)
!      ds = 0.25             ! snow terminal velocity power law coefficient (LFO)
!  new values for  cs and ds
      real, parameter :: cs = 12.42             ! snow terminal velocity power law coefficient 
      real, parameter :: ds = 0.42              ! snow terminal velocity power law coefficient 
      real, parameter :: cp608 = 0.608          ! constant used in conversion of T to Tv

      real, parameter :: gr = 9.8

      real, parameter :: pi = 3.141592653589793
      real, parameter :: piinv = 1./pi
      real, parameter :: pid4 = pi/4.0

!
! max and min mean volumes
!
      real xvrmn, xvrmx0  ! min, max rain volumes
      real xvsmn, xvsmx  ! min, max snow volumes
      real xvfmn, xvfmx  ! min, max frozen drop volumes
      real xvgmn, xvgmx  ! min, max graupel volumes
      real xvhmn, xvhmn0, xvhmx, xvhmx0  ! min, max hail volumes
      real xvhlmn, xvhlmx  ! min, max lg hail volumes

      real, parameter :: dhlmn = 0.3e-3, dhlmx = 40.e-3
      real, parameter :: dhmn0 = 0.3e-3
      real, public :: dhmn = dhmn0, dhmx = -1.

      real, parameter :: cwradn = 2.0e-6, xcradmn = cwradn    ! minimum radius
      real, parameter :: cwradx = 60.e-6, xcradmx = cwradx    ! maximum radius
      real, parameter :: cwc1 = 6.0/(pi*1000.)

!      parameter( xvcmn=4.188e-18 )   ! mks  min volume = 3 micron radius
      real, parameter :: xvcmn=0.523599*(2.*cwradn)**3    ! mks  min volume = 2.5 micron radius
      real, parameter :: xvcmx=0.523599*(2.*xcradmx)**3    ! mks  max volume = 60 micron radius
      real, parameter :: cwmasn = 1000.*xvcmn   ! minimum mass, defined by radius of 5.0e-6
      real, parameter :: cwmasx = 1000.*xvcmx   ! maximum mass, defined by radius of 50.0e-6
      real, parameter :: cwmasn5 = 1000.*0.523599*(2.*5.0e-6)**3 !  5.23e-13

      real, parameter :: xvimn=0.523599*(2.*5.e-6)**3    ! mks  min volume = 5 micron radius
      real, parameter :: xvimx=0.523599*(2.*1.e-3)**3    ! mks  max volume = 1 mm radius (solid sphere approx)
      
      real, public   :: xvdmx = -1.0 ! 3.0e-3
      real     :: xvrmx
      parameter( xvrmn=0.523599*(80.e-6)**3, xvrmx0=0.523599*(6.e-3)**3 ) !( was 4.1887e-9 )  ! mks
      parameter( xvsmn=0.523599*(0.01e-3)**3, xvsmx=0.523599*(10.e-3)**3 ) !( was 4.1887e-9 )  ! mks
      parameter( xvfmn=0.523599*(0.1e-3)**3, xvfmx=0.523599*(10.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvgmn=0.523599*(0.1e-3)**3, xvgmx=0.523599*(10.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvhmn0=0.523599*(0.3e-3)**3, xvhmx0=0.523599*(20.e-3)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3
      parameter( xvhlmn=0.523599*(dhlmn)**3, xvhlmx=0.523599*(dhlmx)**3 )  ! mks xvfmx = (pi/6)*(10mm)**3

!
!  electrical permitivity of air C / (N m**2) -  check the units
!
      real eperao
      parameter (eperao  = 8.8592e-12 )

      real ec,eci  ! fundamental unit of charge
      parameter (ec = 1.602e-19)
      parameter (eci = 1.0/ec)

      real    :: scwppmx = 20.0e-12
      real    :: scippmx = 20.0e-12
!
!  constants
!
      real, parameter :: c1f3 = 1.0/3.0

      real, parameter :: cai = 21.87455
      real, parameter :: caw = 17.2693882
      real, parameter :: cbi = 7.66
      real, parameter :: cbw = 35.86

      real, parameter :: cbwbolton = 29.65 ! constants for Bolton formulation
      real, parameter :: cawbolton = 17.67

      real, parameter :: tfrh = 233.15
      real, parameter :: tfr = 273.15

      real, parameter :: cp = 1004.0, rd = 287.04
      real, parameter :: rw = 461.5              ! gas const. for water vapor
      real, private, parameter ::      cpl = 4190.0
      real, private, parameter ::      cpigb = 2106.0
      real, parameter :: cpi = 1./cp
      real, parameter :: cap = rd/cp
      real, parameter :: tfrcbw = tfr - cbw
      real, parameter :: tfrcbi = tfr - cbi
      real, parameter :: rovcp = rd/cp
      real, public :: rdorv = 0.622

      real, parameter :: poo = 1.0e+05
      real, parameter :: advisc0 = 1.832e-05     ! reference dynamic viscosity (SMT; see Beard & Pruppacher 71)
      real, parameter :: advisc1 = 1.718e-05     ! dynamic viscosity constant used in thermal conductivity calc
      real, parameter :: tka0 = 2.43e-02         ! reference thermal conductivity

     ! GHB: Needed for eqtset=2 in cm1
!     REAL, PRIVATE ::      cv = cp - rd
      real, private, parameter ::      cv = 717.0             ! specific heat at constant volume - air
      real, private, parameter ::      cvv = 1408.5
     ! GHB

      real, parameter ::  bfnu0 = (rnu + 2.0)/(rnu + 1.0) 
      real :: ventr, ventrn, ventc, c1sw


      real :: cckm,ccne,ccnefac,cnexp,CCNE0

      integer :: na = 9
      integer :: nxtra = 1
      real gf4p5, gf4ds, gf4br
      real gsnow1, gsnow53, gsnow73
      real gfcinu1, gfcinu1p47, gfcinu2p47
      real gfcinu1p22,gfcinu2p22
      real gfcinu1p18,gfcinu2p18

      real :: cwchtmp0 = 1.0
      real :: cwchltmp0 = 1.0

      real    :: esctot = 1.0e-13

      integer iexy(lc:lqmx,lc:lqmx)
      integer :: ieswi = 1,  ieswc = 1, ieswr = 0
      integer :: iehlsw = 1, iehli = 1,  iehlc = 1, iehlr = 0
      integer :: iehwsw = 1, iehwi = 1,  iehwc = 1, iehwr = 0

      logical, parameter :: do_satadj_for_wrfchem = .true.

      integer, parameter :: ac_opt = 0 ! option flag for alternate aerosol (for NUWRF only)
      logical, private :: nuaccoinp = .false.

! Note to users: Many of these options are for development and not guaranteed to perform well.
! Some may not be functional depending on the version of the code.
! Some may be useful for ensemble physics diversity. Feel free to contact Ted Mansell if you have questions
! in that regard.
  NAMELIST /nssl_mp_params/               &
                        ndebug, ncdebug,&
                        iusewetgraupel, &
                        iusewethail,    &
                        iusewetsnow,    &
                        idbzci,         &
                        vtmaxsed,       &
                        itfall,iscfall, &
                        infall,irfall,isfall,  &
                        rssflg,         &
                        sssflg,         &
                        hssflg,         &
                        hlssflg,        &
                        irimdenopt,rimdenvwgt,     &
                        rimc1, rimc2, rimc3, rimc4,   &
                        idiagnosecnu,   &
                        icnuclimit,     &
                        irenuc,         &
                        restoreccn, ccntimeconst, cck, &
                        decayufccn, ufccntimeconst, &
                        switchccn, old_cccn,  &
                        ciintmx,        &
                        itype1, itype2, &
                        icenucopt, in_freeze_rain_first,     &
                        naer,           &
                        icfn,           &
                        ibfc, iacr, icracr, &
                        icracrthresh,   &
                        cwfrz2snowfrac, cwfrz2snowratio, &
                        ibfr,           &
                        ibiggopt,       &
                        ibiggsmallrain, &
                        ifrzg,ifiacrg,  &
                        ifrzs,ffrzs,    &
                        iacrsize,       &
                        cimas0, cimas1, cfnfac, &
                        splintermass,   &
                        ewfac,          &
                        eii0, eii1,     &
                        eri0, esi0,     &
                        eri_cimin,      &
                        eii0hl, eii1hl, &
                        ehs0, ehs1,     &
                        ess0, ess1, iessopt,    &
                        esstem1,esstem2, &
                        ircnw, qminrncw,& ! single-moment only
                        iglcnvi,        &
                        iglcnvs,        &
                        alphahacx,      &
                        fconv,          &
                        eqtot,          &
                        imeyers5,       &
                        iehw,           &
                        ierw,           &
                        iehr0c,iehlr0c, &
                        alphai,         &
                        alphar,         &
                        alphas,         & ! note that alphah and alphahl come through physics namelist
                        cnu,            &
                        iscni,fscni,    &
                        dfrz,           &
                        dmlt,           &
                        rainfallfac,    &
                        icefallfac,     &
                        snowfallfac,    &
                        graupelfallfac,    &
                        hailfallfac,    &
                        icefallopt,     &
                        icdx,icdxhl,    &
                        axh,bxh,axf,bxf,axhl,bxhl, &
                        cdhmin, cdhmax,       &
                        cdhdnmin, cdhdnmax,   &
                        cdhlmin, cdhlmax,     &
                        cdhldnmin, cdhldnmax, &
                        ihmlt,          &
                        ehimin,         &
                        ehimax,         &
                        ehsmax,         &
                        ecollmx,        &
                        ehw0, ehlw0,    &
                        ehr0, ehlr0,    &
                        erw0,           &
                        exwmindiam,     &
                        nsplinter,      &
                        lawson_splinter_fac, &
                        iqcinit,        &
                        ssmxinit,      &
                        xvdmx,          &
                        dhmn, dhmx,     &
                        fwms,fwmh,fwmhl,  &
                        ifwmhopt,         &
                        ihxw2rain,        &
                        fwmlarge,         &
                        ifwmfall,         &
                        iturbenhance,     &
                        qsdenmod,qhdenmod, &
                        qsvtmod,          &
                        alphamin,alphamax, &
                        isnwfrac,          &
                        rescale_low_alpha, &
                        rescale_low_alphar, &
                        rescale_low_alphah, &
                        rescale_low_alphahl, &
                        rescale_high_alpha, &
                        ihlcnh, hldia1,iusedw, dwehwmin, dwmin, dwmax, dwtempmin, dg0thresh, &
                        icvhl2h, hldnmn,hdnmn,    &
                        hlcnhdia, hlcnhqmin, &
                        isedonly,           &
                        iresetmoments,      &
                        cxmin, zxmin,       &
                        imurain,            &
                        iferwisventr,       &
                        izwisventr,         &
                        qhdpvdn,            &
                        qhacidn,            &
                        sheddiam,sheddiamlg, &
                        sheddiam0,           &
                        mltdiam1,mltdiam2,mltdiam3,mltdiam4,mltdiam05, &
                        imaxdiaopt,          &
                        ithompsoncnoh,       &
                        cnohmn,             &
                        ivhmltsoak,         &
                        ioldlimiter,        &
                        isnowfall,          &
                        isnowdens,          &
                        ibiggsnow,          &
                        ixtaltype,          &
                        evapfac,    &
                        depfac,             &
                        dmrauto,irescalerainopt, dmropt,dmhlopt,     &
                        rescale_tempthresh, rescale_wthresh, &
                        ibinhmlr,ibinhlmlr,imltshddmr, binmlrmxdia, binmlrzrrfac,ibinnum,   &
                        iqhacrmlr, iqhlacrmlr, &
                        snowmeltdia,    &
                        delta_alphamlr, &
                        iqvsopt,     &
                        maxsupersat, &
                        do_accurate_sedimentation, interval_sedi_vt
! #####################################################################
! #####################################################################

 contains

! #####################################################################
! #####################################################################


 real function fqvs(t)
  implicit none
  real :: t
  fqvs = exp(caw*(t-273.15)/(t-cbw))
 end function fqvs

 real function fqis(t)
  implicit none
  real :: t
  fqis = exp(cai*(t-273.15)/(t-cbi))
 end function fqis
 





! #####################################################################
! #####################################################################
       subroutine nssl_2mom_init(  &
     & ipctmp, mixphase,ihvol,idoniconlytmp,eqtset_tmp, index_ccna, index_cina, outfile, &
     & nssl_graupelfallfac, &
     & nssl_hailfallfac, &
     & nssl_ehw0, &
     & nssl_ehlw0, &
     & nssl_icdx, &
     & nssl_icdxhl, &
     & nssl_icefallfac, &
     & nssl_snowfallfac, &
     & nssl_cccn,   &
     & nssl_ufccn,  &
     & nssl_alphah, &
     & nssl_alphahl, &
     & nssl_alphar, &
     & nssl_density_on, nssl_hail_on, nssl_ccn_on, nssl_icecrystals_on, ccn_is_ccna, &
     & infileunit, &
     & myrank, mpiroot &
     )

  implicit none
  
   real, intent(in), optional ::  &
     & nssl_graupelfallfac, &
     & nssl_hailfallfac, &
     & nssl_ehw0, &
     & nssl_ehlw0, &
     & nssl_icefallfac, &
     & nssl_snowfallfac, &
     & nssl_cccn,   &
     & nssl_alphah, &
     & nssl_alphahl, &
     & nssl_alphar
   integer, intent(in), optional ::  &
     & nssl_icdx, &
     & nssl_icdxhl, myrank, mpiroot, &
     & nssl_ufccn
   logical, intent(in), optional :: nssl_density_on, nssl_hail_on, nssl_ccn_on, nssl_icecrystals_on
   integer, intent(inout), optional :: ccn_is_ccna

  integer, intent(in),optional      :: infileunit

   integer, optional, intent(in) :: eqtset_tmp
   integer, optional, intent(in) :: index_ccna, index_cina, outfile


   integer, intent(in) :: ipctmp,mixphase
   integer, optional, intent(in) :: ihvol
   logical, optional, intent(in) :: idoniconlytmp

    integer :: igvol_local = 1
    logical :: wrote_namelist = .false.
    logical :: wrf_dm_on_monitor
    integer :: hail_on = -1, density_on = -1, icecrystals_on = 1
    integer :: ccn_on = -1

     double precision :: arg
     real    :: temq
     integer :: igam
     integer :: i,il,j,l
     integer :: ltmp
     integer :: isub
     real    :: bxh1,bxhl1

      real    :: alp,ratio
      double precision  :: x,y,y2,y7
      logical :: turn_on_ccna, turn_on_cina
      integer :: iufccn = 0
      integer :: istat

      real :: alpjj, alpii, xnuii, xnujj
      integer :: ii, jj
     

     turn_on_ccna = .false.
     turn_on_cina = .false.

!      IF ( present( igvol ) ) THEN
!        igvol_local = igvol
!      ENDIF
      
      if ( present( nssl_hail_on ) ) then
        if ( nssl_hail_on ) then
          hail_on = 1
        else
          hail_on = 0
        endif
      endif

      if ( present( nssl_density_on ) ) then
        if ( nssl_density_on ) then
          density_on = 1
        else
          density_on = 0
        endif
      endif
      
      if ( present( nssl_icecrystals_on ) ) then
        if ( nssl_icecrystals_on ) then
          icecrystals_on = 1
        else
          icecrystals_on = 0
          ! renucfrac = 1.0 ! why was this set to 1?
          ffrzs = 1.0
        endif
      endif


!
! set some global values from namelist input
!
!      open(unit=20,file='namelist.input',form='formatted',status='old',    &
!           access='sequential')
!         read(20,nml=micro_params)
!      close(unit=20)

      if ( present( eqtset_tmp ) ) eqtset = eqtset_tmp
      if ( present( index_ccna ) ) then
        if ( index_ccna > 1 ) turn_on_ccna = .true.
      endif
      if ( present( index_cina ) ) then
        if ( index_cina > 1 ) turn_on_cina = .true.
      endif



      if ( present(nssl_graupelfallfac) ) graupelfallfac = nssl_graupelfallfac
      if ( present(nssl_hailfallfac) ) hailfallfac = nssl_hailfallfac
      if ( present(nssl_ehw0) ) then
        if ( nssl_ehw0 > 0.0 ) ehw0 = nssl_ehw0
      endif
      if ( present(nssl_ehlw0) ) then
        if ( nssl_ehlw0 > 0.0 ) ehlw0 = nssl_ehlw0
      endif
      if ( present(nssl_icdx) ) icdx = nssl_icdx
      if ( present(nssl_icdxhl) ) icdxhl = nssl_icdxhl
      if ( present(nssl_icefallfac) ) icefallfac = nssl_icefallfac
      if ( present(nssl_snowfallfac) ) snowfallfac = nssl_snowfallfac
      if ( present(nssl_cccn) ) then
        if (nssl_cccn > 1 ) ccn = nssl_cccn
      endif
      if ( present(nssl_alphah) ) then
        if ( nssl_alphah > -1. ) alphah = nssl_alphah
      endif
      if ( present(nssl_alphahl) ) then
        if ( nssl_alphahl > -1. ) alphahl = nssl_alphahl
      endif
      if ( present(nssl_alphar) ) then
        if ( nssl_alphar > -1.0 ) alphar = nssl_alphar
      endif


    ipconc = ipctmp
    
    if ( ipconc < 5 ) then
       ihlcnh = 0
    endif

    if ( ihlcnh <= 0 ) then
      if ( ipconc == 5 ) then
       ihlcnh = 3
      elseif ( ipconc >= 6 ) then
       ihlcnh = 3
      endif
    endif

      



      open(15,file='namelist.input',status='old',form='formatted',action='read')
      rewind(15)
      read(15,NML=nssl_mp_params,iostat=istat)
      close(15)
      if ( istat /= 0 ) then
#ifdef WRF_ELEC
        if ( wrf_dm_on_monitor() ) then
        write(0,*) 'NSSL_2MOM_INIT: PROBLEM WITH NSSL_MP_PARAMS namelist: not found or bad token'
        endif
#else
       ! write(0,*) 'NSSL_2MOM_INIT: PROBLEM WITH NSSL_MP_PARAMS namelist: not found or bad token'
#endif
      endif
      if ( istat /= 0 ) then
        if ( present(myrank) ) then
          if ( myrank == 0 ) write(0,*) 'READ_NAMELIST: PROBLEM WITH NSSL_MP_PARAMS namelist: not found or bad token'
        endif
      else
        if ( present(myrank) ) then
          if ( myrank == 0 ) write(0,*) 'read nssl_mp_params'
        endif
      endif
!        IF ( present( outfile ) ) THEN
!          IF ( outfile >= 0 ) THEN
!          write(outfile,NML=nssl_mp_params)
!          ENDIF
!        ENDIF
        if ( present(myrank) ) then
          myrankcm1 = myrank
          if ( myrank == 0 ) then
           open(15,file='nssl_mp_params.out',status='unknown',form='formatted')
           write(15,NML=nssl_mp_params)
           close(15)
          endif
        endif



      if ( iufccn > 0 ) then ! make sure to use option that uses UF ccn
        irenuc = 7
        if ( ccnuf <= 0.0 ) decayufccn = .true. ! assume surface emission and need decay
        if ( i_uf_or_ccn > 0 ) then
          ufbackground = 0.0
          ccntimeconst = ufccntimeconst
        endif
      endif

        if ( present( nssl_ccn_on ) ) then
          if ( nssl_ccn_on ) then
            ccn_on = 1
          else
            ccn_on = 0
            irenuc = 2
          endif
        endif

      if ( irenuc >= 5 ) then
        turn_on_ccna = .true.
        if ( present( nssl_ccn_on ) ) then
          if ( .not. nssl_ccn_on ) then
      write(0,*) 'NSSL_MP Error: Must have nssl_ccn_on=1 for irenuc >= 5!'
        stop
          endif
        endif
      endif

      if ( present( ccn_is_ccna ) .and. ccn_on == 1 ) then
        if ( ccn_is_ccna > 0 ) then
          turn_on_ccna = .true.
        else
          if ( irenuc >= 5 ) then
            ccn_is_ccna = 1
          endif
        endif
      endif

      cwccn = ccn

      lhab = 8
      lhl = 8
      if ( icespheres >= 1 ) then
        lhab = lhab + 1
        lis = li + 1
        ls = ls + 1
        lh = lh + 1
        lhl = lhl + 1
      endif
      if ( hail_on == -1 ) then ! hail_on is not set
        hail_on = 1
        if ( ihvol <= -1 .or. ihvol == 2 ) then
          if ( ihvol == -1 .or. ihvol == -2 ) then
          lhab = lhab - 1  ! turns off hail 
          lhl = 0
          hail_on = 0
          ! past me thought it would be a good idea to change graupel factors when hail is off....
          ! ehw0 = 0.75
          ! iehw = 2
          ! dfrz = Max( dfrz, 0.5e-3 )
          endif
          if ( ihvol == -2 .or. ihvol == 2 .or. icecrystals_on == 0 ) then ! ice crystals are turned off
           ! a value of 2? means to turn off ice crystals but turn on hail
           ! renucfrac = 1.0 ! why?
            ffrzs = 1.0
           ! idoci = 0 ! try this later
          endif
        endif
      
      else ! hail_on is set
        if ( hail_on == 0 ) then
          lhab = lhab - 1  ! turns off hail 
          lhl = 0
        else
          ! assume default that hail is on
        endif
      endif
      
      if ( density_on == -1 ) then ! density flag not set, so default is to predict it
        density_on = 1
      endif


!      write(0,*) 'wrf_init: lhab,lhl,hail_on,density_on = ',lhab,lhl,hail_on,density_on

!      IF ( ipelec > 0 ) idonic = .true.

!
! Build lookup table for saturation mixing ratio (Soong and Ogura 73)
!

      do l = 1,nqsat
      temq = 163.15 + (l-1)*fqsat
      if ( iqvsopt == 0 ) then
      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
      dtabqvs(l) = ((-caw*(-273.15 + temq))/(temq - cbw)**2 + &
     &                 caw/(temq - cbw))*tabqvs(l)
      else
      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
      dtabqvs(l) = ((-cawbolton*(-273.15 + temq))/(temq - cbwbolton)**2 + &
     &                 cawbolton/(temq - cbwbolton))*tabqvs(l)
      endif
      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
      dtabqis(l) = ((-cai*(-273.15 + temq))/(temq - cbi)**2 + &
     &                 cai/(temq - cbi))*tabqis(l)
      end do

      bx(lr) = 0.85
      ax(lr) = 1647.81
      fx(lr) = 135.477

      
      if ( icdx == 6 ) then
        bx(lh) = 0.6 ! Milbrandt and Morrison (2013) for density of 550.
        ax(lh) = 157.71
!      ELSEIF ( icdx == 1 ) THEN
!        bx(lh) = bxh
!        ax(lh) = axh
      elseif ( icdx > 1 ) then
        bx(lh) = 0.5
        ax(lh) = 75.7149
      elseif ( icdx == 0 ) then
        bx(lh) = 0.37 ! 0.6  ! Ferrier 1994 graupel
        ax(lh) = 19.3
      else ! icdx < 0
!        ax(lh) = 206.984 ! Ferrier 1994 hail/frozen drops
!        bx(lh) = 0.6384
        bx(lh) = bxh
        ax(lh) = axh
      endif

!      bx(lh) = 0.6

      if ( lhl .gt. 1 ) then
        if ( icdxhl == 6 ) then
          bx(lhl) = 0.593 ! Milbrandt and Morrison (2013) for density of 750.
          ax(lhl) = 179.36
        elseif (icdxhl == 0 ) then
          ax(lhl) = 206.984 ! Ferrier 1994
          bx(lhl) = 0.6384
        elseif (icdxhl > 0 ) then
         bx(lhl) = 0.5
         ax(lhl) = 75.7149
        else
         bx(lhl) = bxhl
         ax(lhl) = axhl
        endif
      endif

! fill in the complete gamma function lookup table
     gmoi(0) = 1.d32
     do igam = 1,ngm0
      arg = dgam*igam
      gmoi(igam) = gamma_dp(arg)
     end do

     ! build lookup table to compute the number and mass fractions of rain drops 
     ! (imurain=1) greater than a given diameter. Used for qiacr and ciacr
     ! Uses incomplete gamma functions
     ! The terms with bxh or bxhl will be off if the actual bxh or bxhl is different from the base value (icdx=6 option)
      
      bxh1 = bx(lh)
      bxhl1 = bx(Max(lh,lhl))
      
!      DO j = 0,nqiacralpha
      do j = ialpstart,nqiacralpha
      alp = float(j)*dqiacralpha
      y = gamma_dpr(1.+alp)
      y2 = gamma_dpr(2.+alp)
      do i = 0,nqiacrratio
        ratio = float(i)*dqiacrratio
        x = gamxinfdp( 1.+alp, ratio )
!        write(0,*) 'i, x/y = ',i, x/y
        ciacrratio(i,j) = x/y

        ! graupel (.,.,.,1)
        gamxinflu(i,j,1,1) = x/y
        gamxinflu(i,j,2,1) = gamxinfdp( 2.0+alp, ratio )/y
        gamxinflu(i,j,3,1) = gamxinfdp( 2.5+alp+0.5*bxh1, ratio )/y
        gamxinflu(i,j,5,1) = (gamma_dpr(5.0+alp) - gamxinfdp( 5.0+alp, ratio ))/y
        gamxinflu(i,j,6,1) = (gamma_dpr(5.5+alp+0.5*bxh1) - gamxinfdp( 5.5+alp+0.5*bxh1, ratio ))/y
        gamxinflu(i,j,9,1) = gamxinfdp( 1.0+alp, ratio )/y
        gamxinflu(i,j,10,1)= gamxinfdp( 4.0+alp, ratio )/y

        gamxinflu(i,j,12,1) = gamxinfdp( 2.0+alp, ratio )/y2
       
        ! hail (.,.,.,2)
        gamxinflu(i,j,1,2) = gamxinflu(i,j,1,1)
        gamxinflu(i,j,2,2) = gamxinflu(i,j,2,1)
        gamxinflu(i,j,3,2) = gamxinfdp( 2.5+alp+0.5*bxhl1, ratio )/y
        gamxinflu(i,j,5,2) = gamxinflu(i,j,5,1)
        gamxinflu(i,j,6,2) = (gamma_dpr(5.5+alp+0.5*bxhl1) - gamxinfdp( 5.5+alp+0.5*bxhl1, ratio ))/y
        gamxinflu(i,j,9,2) = gamxinflu(i,j,9,1)
        gamxinflu(i,j,10,2)= gamxinflu(i,j,10,1)

      if ( alp > 1.1 ) then
!       gamxinflu(i,j,7,1) = gamxinfdp( alp - 1., ratio )/y
       gamxinflu(i,j,7,1) = (gamma_dpr(alp - 1.) - gamxinfdp( alp - 1., ratio ))/y
!       gamxinflu(i,j,8,1) = gamxinfdp( alp - 0.5 + 0.5*bxh, ratio )/y
       gamxinflu(i,j,8,1) = (gamma_dpr(alp - 0.5 + 0.5*bxh1) - gamxinfdp( alp - 0.5 + 0.5*bxh1, ratio ))/y
!       gamxinflu(i,j,8,2) = gamxinfdp( alp - 0.5 + 0.5*bxhl1, ratio )/y
       gamxinflu(i,j,8,2) = (gamma_dpr(alp - 0.5 + 0.5*bxhl1) - gamxinfdp( alp - 0.5 + 0.5*bxhl1, ratio ))/y
      else
!       gamxinflu(i,j,7,1) = gamxinfdp( .1, ratio )/y
       gamxinflu(i,j,7,1) = (gamma_dpr(0.1) - gamxinfdp( 0.1, ratio ) )/y
!       gamxinflu(i,j,8,1) = gamxinfdp( 1.1 - 0.5 + 0.5*bxh1, ratio )/y
!       gamxinflu(i,j,8,2) = gamxinfdp( 1.1 - 0.5 + 0.5*bxhl1, ratio )/y
       gamxinflu(i,j,8,1) = (gamma_dpr(1.1 - 0.5 + 0.5*bxh1) - gamxinfdp( 1.1 - 0.5 + 0.5*bxh1, ratio ) )/y
       gamxinflu(i,j,8,2) = (gamma_dpr(1.1 - 0.5 + 0.5*bxhl1) - gamxinfdp( 1.1 - 0.5 + 0.5*bxhl1, ratio ) )/y
      endif
        
        gamxinflu(i,j,7,2) = gamxinflu(i,j,7,1)

      enddo
      enddo
      ciacrratio(0,:) = 1.0

      do j = ialpstart,nqiacralpha
      alp = float(j)*dqiacralpha
      y = gamma_sp(4.+alp)
      y7 = gamma_sp(7.+alp)
      do i = 0,nqiacrratio
        ratio = float(i)*dqiacrratio
        
        ! mass fraction
        x = gamxinfdp( 4.+alp, ratio )
!        write(0,*) 'i, x/y = ',i, x/y
        qiacrratio(i,j) = x/y
        gamxinflu(i,j,4,1) = x/y
        gamxinflu(i,j,4,2) = x/y

        ! reflectivity fraction
        x = gamxinfdp( 7.+alp, ratio )
        ziacrratio(i,j) = x/y7
        gamxinflu(i,j,11,1) = x/y7
        gamxinflu(i,j,11,2) = x/y7

      enddo
      enddo
      qiacrratio(0,:) = 1.0


      lccn = 0
      lccnuf = 0
      lccna = 0
      lnc = 0
      lnr = 0
      lni = 0
      lnis = 0
      lns = 0
      lnh = 0
      lnhl = 0
      lvh = 0
      lvhl = 0
      lzr = 0
      lzh = 0
      lzhl = 0
      lsw = 0
      lhw = 0
      lhlw = 0

      denscale(:) = 0
      
!      lccn = 9


    if ( ipconc == 0 ) then
       if ( hail_on == 1 ) then ! turn on graupel density for 1-moment scheme
       lvh = 9
       ltmp = 9
       denscale(lvh) = 1
       else ! no hail, 'LFO' scheme
       ltmp = lhab
       lhl = 0
       endif
    elseif ( ipconc == 5 ) then
      ltmp = lhab
      if ( iufccn > 0 ) then
        ltmp = ltmp+1
        lccnuf = ltmp
        denscale(lccnuf) = 1
      endif
      lccn= ltmp+1 ! 9
      lnc = ltmp+2 ! 10
      lnr = ltmp+3 ! 11
      lni = ltmp+4 !12
      lns = ltmp+5 !13
      lnh = ltmp+6 !14
      ltmp = lnh
      if ( hail_on == 1 ) then
      ltmp = ltmp + 1
      lnhl = ltmp ! lhab+7 ! 15
      endif
      if ( density_on >= 1 ) then
      ltmp = ltmp + 1
      lvh = ltmp ! lhab+8 + isub ! 16 + isub ! isub adjusts to 15 if hail is off
!      ltmp = lvh
      endif
      denscale(lccn:ltmp) = 1
      if ( density_on == 1 .and. hail_on == 1 ) then
       ltmp = ltmp + 1
       lvhl = ltmp
!       ltmp = lvhl
       denscale(lvhl) = 1
      endif
      if ( mixedphase ) then
      ltmp = ltmp + 1
      lsw  = ltmp
      ltmp = ltmp + 1
      lhw  = ltmp
        if ( lhl > 1 ) then
          ltmp = ltmp + 1
          lhlw = ltmp
        endif
!      ltmp = lhlw
      endif
    elseif ( ipconc >= 6 ) then
      ltmp = lhab
      if ( iufccn > 0 ) then
        ltmp = ltmp+1
        lccnuf = ltmp
        denscale(lccnuf) = 1
      endif

      lccn= ltmp+1 ! 9
      lnc = ltmp+2 ! 10
      lnr = ltmp+3 ! 11
      lni = ltmp+4 !12
      lns = ltmp+5 !13
      lnh = ltmp+6 !14
      ltmp = lnh
      if ( lhl > 0 ) then
      ltmp = ltmp + 1
      lnhl = ltmp ! lhab+7 ! 15
      endif
      if ( density_on == 1 ) then
      ltmp = ltmp + 1
      lvh = ltmp ! lhab+8 + isub ! 16 + isub ! isub adjusts to 15 if hail is off
      endif
!      ltmp = lvh
      denscale(lccn:ltmp) = 1
      if ( density_on == 1 .and. hail_on == 1 ) then
       ltmp = ltmp + 1
       lvhl = ltmp
!       ltmp = lvhl
       denscale(lvhl) = 1
      endif

      if ( ipconc == 6 ) then
       ltmp = ltmp + 1
       lzh = ltmp
      elseif ( ipconc == 7 ) then
       ltmp = ltmp + 1
       lzh = ltmp
       ltmp = ltmp + 1
       lzr = ltmp
      elseif ( ipconc == 8 ) then
       ltmp = ltmp + 1
       lzh = ltmp
       ltmp = ltmp + 1
       lzr = ltmp
       if ( lhl > 1 ) then
         ltmp = ltmp + 1
         lzhl = ltmp
       endif
      ! write(0,*) 'ipcon,lzr = ',ipconc,lzr,lzh,lzhl
      endif
!      ltmp = lvh
 !     denscale(lccn:lvh) = 1
      if ( mixedphase ) then
      ltmp = ltmp + 1
      lsw  = ltmp
      ltmp = ltmp + 1
      lhw  = ltmp
        if ( lhl > 1 ) then
          ltmp = ltmp + 1
          lhlw = ltmp
        endif
!      ltmp = lhlw
      endif
    else
      write(0,*) 'nssl_2mom_init: Invalid value of ipctmp'
      stop
    endif



      ! write(0,*) 'wrf_init: lh,lhl,lzh,lzhl = ',lh,lhl,lzh,lzhl 
      ! write(0,*) 'wrf_init: ipconc = ',ipconc
      ! write(0,*) 'wrf_init: irenuc, turn_on_ccna = ',irenuc, turn_on_ccna 
      if ( turn_on_ccna ) then
        ltmp = ltmp + 1
        lccna = ltmp
        denscale(ltmp) = 1
      endif

      if ( turn_on_cina ) then
        ltmp = ltmp + 1
        lcina = ltmp
        denscale(ltmp) = 1
      endif

      if ( turn_on_cin .or. is_aerosol_aware ) then
        ltmp = ltmp + 1
        lcin = ltmp
        denscale(ltmp) = 1
!debug        write(0,*) 'Setting lcin to ',lcin
      endif
      na = ltmp
      
      ln(lc) = lnc
      ln(lr) = lnr
      ln(li) = lni
      ln(ls) = lns
      ln(lh) = lnh
      if ( lhl .gt. 1 ) ln(lhl) = lnhl

      ipc(lc) = 2
      ipc(lr) = 3
      ipc(li) = 1
      ipc(ls) = 4
      ipc(lh) = 5
      if ( lhl .gt. 1 ) ipc(lhl) = 5
      
      ldovol = .false.
      lvol(:) = 0
      lvol(li) = lvi
      lvol(ls) = lvs
      lvol(lh) = lvh
      if ( lhl .gt. 1 .and. lvhl .gt. 1 ) lvol(lhl) = lvhl
      
      lne = Max(lnh,lnhl)
      lne = Max(lne,lvh)
      lne = Max(lne,lvhl)
      lne = Max(lne,na)

      lsc(:) = 0
      lsc(lc) = lscw
      lsc(lr) = lscr
      lsc(li) = lsci
      lsc(ls) = lscs
      lsc(lh) = lsch
      if ( lhl .gt. 1 ) lsc(lhl) = lschl


      do il = lc,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      enddo

!      write(0,*) 'nssl_2mom_init: ldovol = ',ldovol

      lz(:) = 0
      lz(lr) = lzr
      lz(li) = lzi
      lz(ls) = lzs
      lz(lh) = lzh
      if ( lhl .gt. 1 .and. lzhl > 1 ) lz(lhl) = lzhl

      lliq(:) = 0
      lliq(ls) = lsw
      lliq(lh) = lhw
      if ( lhl .gt. 1 ) lliq(lhl) = lhlw
      if ( mixedphase ) then
!       write(0,*) 'lsw,lhw,lhlw = ',lsw,lhw,lhlw
      endif



      xnu(lc) = cnu
      xmu(lc) = 1.
      
      if ( imurain == 3 ) then
        xnu(lr) = rnu
        xmu(lr) = 1.
      elseif ( imurain == 1 ) then
        xnu(lr) = (alphar - 2.0)/3.0
        xmu(lr) = 1./3.
      endif

      xnu(li) = cinu
      xmu(li) = 1.

      if ( lis >= 1 ) then
      xnu(lis) = 0.0
      xmu(lis) = 1.
      endif

      dnu(lc) = 3.*xnu(lc) + 2. ! alphac
      dmu(lc) = 3.*xmu(lc)

      dnu(lr) = 3.*xnu(lr) + 2. ! alphar
      dmu(lr) = 3.*xmu(lr)

      xnu(ls) = snu
      xmu(ls) = 1.

      dnu(ls) = 3.*xnu(ls) + 2.  ! -0.4 ! alphas
      dmu(ls) = 3.*xmu(ls)


      dnu(lh) = alphah
      dmu(lh) = dmuh

      xnu(lh) = (dnu(lh) - 2.)/3.
      xmu(lh) = dmuh/3.


      if ( imurain == 3 ) then ! rain is gamma of volume
      rz =  ((4. + alphah)*(5. + alphah)*(6. + alphah)*(1. + xnu(lr)))/ & 
     &  ((1 + alphah)*(2 + alphah)*(3 + alphah)*(2. + xnu(lr)))

!      IF ( ipconc .lt. 5 ) alphahl = alphah
      
      rzhl =  ((4. + alphahl)*(5. + alphahl)*(6. + alphahl)*(1. + xnu(lr)))/ & 
     &  ((1. + alphahl)*(2. + alphahl)*(3. + alphahl)*(2. + xnu(lr)))

      rzs =  1. ! assume rain and snow are both gamma volume

      else ! rain is gamma of diameter
      
      rz =  ((4. + alphah)*(5. + alphah)*(6. + alphah)*(1. + alphar)*(2. + alphar)*(3. + alphar))/ & 
     &  ((1 + alphah)*(2 + alphah)*(3 + alphah)*(4. + alphar)*(5. + alphar)*(6. + alphar))
      
      rzhl =  ((4. + alphahl)*(5. + alphahl)*(6. + alphahl)*(1. + alphar)*(2. + alphar)*(3. + alphar))/ & 
     &  ((1 + alphahl)*(2 + alphahl)*(3 + alphahl)*(4. + alphar)*(5. + alphar)*(6. + alphar))

      
      rzs =   & 
     &  ((1. + alphar)*(2. + alphar)*(3. + alphar)*(2. + xnu(ls)))/  &
     &  ((4. + alphar)*(5. + alphar)*(6. + alphar)*(1. + xnu(ls)))
       

      endif

      if ( ipconc <= 5 ) then 
        imltshddmr = Min(1, imltshddmr)
        ibinhmlr = 0
        ibinhlmlr = 0
      endif

      if ( ipconc > 5 .and. (ibinhmlr == 0 .and. ibinhlmlr == 0 ) ) then 
        imltshddmr = Min(1, imltshddmr)
      endif

!      write(0,*) 'rz,rzhl = ', rz,rzhl

      if ( ipconc .lt. 4 ) then

      dnu(ls) = alphas
      dmu(ls) = 1.

      xnu(ls) = (dnu(ls) - 2.)/3.
      xmu(ls) = 1./3.


      endif

      if ( lhl .gt. 1 ) then

      dnu(lhl) = alphahl
      dmu(lhl) = dmuhl

      xnu(lhl) = (dnu(lhl) - 2.)/3.
      xmu(lhl) = dmuhl/3.

      endif

      cno(lc)  = 1.0e+08
      if ( li .gt. 1 ) cno(li)  = 1.0e+08
      cno(lr)  = cnor
      if ( ls .gt. 1 ) cno(ls)  = cnos ! 8.0e+06
      if ( lh .gt. 1 ) cno(lh)  = cnoh ! 4.0e+05
      if ( lhl .gt. 1 ) cno(lhl)  = cnohl ! 4.0e+05
!
!  density maximums and minimums
!
      xdnmx(:) = 900.0

      xdnmx(lr) = 1000.0
      xdnmx(lc) = 1000.0
      xdnmx(li) =  917.0
      xdnmx(ls) =  300.0
      xdnmx(lh) =  900.0
      if ( lhl .gt. 1 ) xdnmx(lhl) = 900.0
!
      xdnmn(:) = 900.0

      xdnmn(lr) = 1000.0
      xdnmn(lc) = 1000.0
      xdnmn(li) =  100.0
      xdnmn(ls) =  100.0
      xdnmn(lh) =  hdnmn
      if ( lhl .gt. 1 ) xdnmn(lhl) = hldnmn

      xdn0(:) = 900.0

      xdn0(lc) = 1000.0
      xdn0(li) = 900.0
      xdn0(lr) = 1000.0
      xdn0(ls) = rho_qs ! 100.0
      xdn0(lh) = rho_qh ! (0.5)*(xdnmn(lh)+xdnmx(lh))
      if ( lhl .gt. 1 ) xdn0(lhl) = rho_qhl ! 800.0

!
!  Set terminal velocities...
!    also set drag coefficients
!
      cdx(lr) = 0.60
      cdx(lh) = 0.8 ! 1.0 ! 0.45
      cdx(ls) = 2.00
      if ( lhl .gt. 1 ) cdx(lhl) = 0.45

      ido(lc) = idocw
      ido(lr) = idorw
      ido(li) = idoci
      ido(ls) = idosw
      ido(lh)  = idohw
      if ( lhl .gt. 1 ) ido(lhl) = idohl

      if ( irfall .lt. 0 ) irfall = infall
      if ( isfall .lt. 0 ) isfall = infall
      if ( lzr > 0 ) irfall = 0

      qccn = ccn/rho00
      qccnuf = ccnuf/rho00
      if ( old_cccn > 0.0 ) then
         old_qccn = old_cccn/rho00
      else
         old_qccn = qccn
      endif
!      xvcmx = (4./3.)*pi*xcradmx**3

! set max rain diameter
      if ( xvdmx .gt. 0.0 ) then
        xvrmx = 0.523599*(xvdmx)**3
      else
        xvrmx = xvrmx0
      endif

         if ( dhmn <= 0.0 ) then
           xvhmn = xvhmn0
!           xvhmn = Min(xvhmn0, 0.523599*(dfrz)**3 )
         else
           xvhmn = 0.523599*(dhmn)**3
!           xvhmn = 0.523599*(Min(dhmn,dfrz))**3
         endif

         if ( dhmx <= 0.0 ) then
           xvhmx = xvhmx0
         else
           xvhmx = 0.523599*(dhmx)**3
         endif
         
         if ( qhdpvdn < 0. ) qhdpvdn = xdnmn(lh)
         if ( qhacidn < 0. ) qhacidn = xdnmn(lh)

! load max/min diameters
      xvmn(lc) = xvcmn
      xvmn(li) = xvimn
      xvmn(lr) = xvrmn
      xvmn(ls) = xvsmn
      xvmn(lh) = xvhmn

      xvmx(lc) = xvcmx
      xvmx(li) = xvimx
      xvmx(lr) = xvrmx
      xvmx(ls) = xvsmx
      xvmx(lh) = xvhmx

      if ( lhl .gt. 1 ) then
      xvmn(lhl) = xvhlmn
      xvmx(lhl) = xvhlmx
      endif

!
!  cloud water constants in mks units
!
!      cwmasn = 4.25e-15  ! radius of 1.0e-6
!      cwmasn = 5.23e-13   ! minimum mass, defined by radius of 5.0e-6
!      cwmasn5 =  5.23e-13
!      cwradn = 5.0e-6     ! minimum radius
!      cwmasx = 5.25e-10   ! maximum mass, defined by radius of 50.0e-6
!      mwfac = 6.0**(1./3.)
      if ( ipconc .ge. 2 ) then
!        cwmasn = xvmn(lc)*1000.  ! minimum mass, defined by minimum droplet volume
!        cwradn = 1.0e-6          ! minimum radius
!        cwmasx = xvmx(lc)*1000.  ! maximum mass, defined by maximum droplet volume
        
      endif
!        rwmasn = xvmn(lr)*1000.  ! minimum mass, defined by minimum rain volume
!        rwmasx = xvmx(lr)*1000.  ! maximum mass, defined by maximum rain volume

      if ( lhl < 1 ) ifrzg = 1

      ventr = 1.
      if ( imurain == 3 ) then
!       IF ( izwisventr == 1 ) THEN
        ventr = Gamma_sp(rnu + 4./3.)/((rnu + 1.)**(1./3.)*Gamma_sp(rnu + 1.)) ! Ziegler 1985
!       ELSE
        ventrn =  Gamma_sp(rnu + 1.5 + br/6.)/(Gamma_sp(rnu + 1.)*(rnu + 1.)**((1.+br)/6. + 1./3.) ) ! adapted from Wisner et al. 1972; for second term in rwvent
!        ventr = Gamma_sp(rnu + 4./3.)/((rnu + 1.)**(1./3.)*Gamma_sp(rnu + 1.)) ! Ziegler 1985, still use for first term in rwvent
!        ventr  = Gamma_sp(rnu + 4./3.)/Gamma_sp(rnu + 1.) 
!       ENDIF
      else ! imurain == 1
!       IF ( iferwisventr == 1 ) THEN
        ventr = Gamma_sp(2. + alphar)  ! Ferrier 1994
!       ELSEIF ( iferwisventr == 2 ) THEN
        ventrn =  Gamma_sp(alphar + 2.5 + br/2.)/Gamma_sp(alphar + 1.) ! adapted from Wisner et al. 1972
!       ENDIF
      endif
      ventc   = Gamma_sp(cnu + 4./3.)/(cnu + 1.)**(1./3.)/Gamma_sp(cnu + 1.)
      c1sw = Gamma_sp(snu + 4./3.)*(snu + 1.0)**(-1./3.)/gamma_sp(snu + 1.0)

  ! set threshold mixing ratios

      qxmin(:) = 1.0e-12

      qxmin(lc) = 1.e-9
      qxmin(lr) = 1.e-7
      if ( li > 1 ) qxmin(li) = 1.e-12
      if ( ls > 1 ) qxmin(ls) = 1.e-7
      if ( lh > 1 ) qxmin(lh) = 1.e-7
      if ( lhl .gt. 1 ) qxmin(lhl) = 1.e-7

      if ( lc .gt. 1 .and. lnc .gt. 1 ) qxmin(lc) = 1.0e-13
      if ( lr .gt. 1 .and. lnr .gt. 1 ) qxmin(lr) = 1.0e-12

      if ( li .gt. 1 .and. lni .gt. 1 ) qxmin(li ) = 1.0e-13
      if ( ls .gt. 1 .and. lns .gt. 1 ) qxmin(ls ) = 1.0e-13
      if ( lh .gt. 1 .and. lnh .gt. 1 ) qxmin(lh ) = 1.0e-12
      if ( lhl.gt. 1 .and. lnhl.gt. 1 ) qxmin(lhl) = 1.0e-12

      qxmin_init(:) = 1.0e-8 ! threshold for considering single-moment initial condition mixing ratios
  ! constants for droplet nucleation

      cckm = cck-1.
      ccnefac =  (1.63/(cck * beta(3./2., cck/2.)))**(cck/(cck + 2.0))
      cnexp   = (3./2.)*cck/(cck+2.0)
! ccne is all the factors with w in eq. A7 in Mansell et al. 2010 (JAS).  The constant changes
! if k (cck) is changed!
      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
      ccne0 = ccnefac*1.e6*(1.e-6)**(2./(2.+cck))
!      write(0,*) 'cwccn, cck, ccne = ',cwccn,cck,ccne,ccnefac,cnexp
      if ( cwccn .lt. 0.0 ) then
      cwccn = Abs(cwccn)
      ccwmx = 50.e9 ! cwccn
      else
      ccwmx = 50.e9 ! cwccn ! *1.4
      endif

!
!
!  Set collection coefficients (Seifert and Beheng 05)
!
      bb(:) = 1.0/3.0
      bb(li) = 0.3429
      do il = lc,lhab
        da0(il) = delbk(bb(il), xnu(il), xmu(il), 0)
        da1(il) = delbk(bb(il), xnu(il), xmu(il), 1)

!        write(0,*) 'il, da0, da1, xnu, xmu = ', il, da0(il), da1(il), xnu(il), xmu(il)
      enddo

      dab0(:,:) = 0.0
      dab1(:,:) = 0.0

      do il = lc,lhab
        do j = lc,lhab
          if ( il .ne. j ) then

            dab0(il,j) = delabk(bb(il), bb(j), xnu(il), xnu(j), xmu(il), xmu(j), 0)
            dab1(il,j) = delabk(bb(il), bb(j), xnu(il), xnu(j), xmu(il), xmu(j), 1)

!           write(0,*) 'il, j, dab0, dab1 = ',il, j, dab0(il,j), dab1(il,j)
          endif
        enddo
      enddo

      dab0lu(:,:,:,:) = 0.0
      dab1lu(:,:,:,:) = 0.0
      
      if ( ipconc >= 6 ) then
      do il = lc,lhab ! collector
        do j = lc,lhab ! collected
          if ( il .ne. j ) then

            do jj = ialpstart,nqiacralpha
                alpjj = float(jj)*dqiacralpha
                xnujj = (alpjj - 2.)/3.
            do ii = ialpstart,nqiacralpha
                alpii = float(ii)*dqiacralpha
                xnuii = (alpii - 2.)/3.
          
            dab0lu(ii,jj,il,j) = delabk(bb(il), bb(j), xnuii, xnujj, xmu(il), xmu(j), 0)
            dab1lu(ii,jj,il,j) = delabk(bb(il), bb(j), xnuii, xnujj, xmu(il), xmu(j), 1)
          
            enddo
            enddo
!           write(0,*) 'il, j, dab0, dab1 = ',il, j, dab0(il,j), dab1(il,j)
          endif
        enddo
      enddo
      
      endif

        gf4br = gamma_sp(4.0+br)
        gf4ds = gamma_sp(4.0+ds)
        gf4p5 = gamma_sp(4.0+0.5)
        gfcinu1 = gamma_sp(cinu + 1.0)
        gfcinu1p47 = gamma_sp(cinu + 1.47167)
        gfcinu2p47 = gamma_sp(cinu + 2.47167)
        gfcinu1p22 = gamma_sp(cinu + 1.22117)
        gfcinu2p22 = gamma_sp(cinu + 2.22117)
        gfcinu1p18 = gamma_sp(cinu + 1.18333)
        gfcinu2p18 = gamma_sp(cinu + 2.18333)
        
        gsnow1 = gamma_sp(snu + 1.0)
        gsnow53 = gamma_sp(snu + 5./3.)
        gsnow73 = gamma_sp(snu + 7./3.)

        if ( lh  .gt. 1 ) cwchtmp0 = 6.0/pi*gamma_sp( (xnu(lh) + 1.)/xmu(lh) )/gamma_sp( (xnu(lh) + 2.)/xmu(lh) )
        if ( lhl .gt. 1 ) cwchltmp0 = 6.0/pi*gamma_sp( (xnu(lhl) + 1)/xmu(lhl) )/gamma_sp( (xnu(lhl) + 2)/xmu(lhl) )


      iexy(:,:)=0; ! sets to zero the ones Imight have forgotten

!     snow
      iexy(ls,li) = ieswi
      iexy(ls,lc) = ieswc ; iexy(ls,lr) = ieswr ;

!     graupel
      iexy(lh,ls)  = iehwsw ; iexy(lh,li) = iehwi ;
      iexy(lh,lc) = iehwc ; iexy(lh,lr)  = iehwr ;

!     hail
      if (lhl .gt. 1 ) then
      iexy(lhl,ls)  = iehlsw ; iexy(lhl,li) = iehli ;
      iexy(lhl,lc) = iehlc ; iexy(lhl,lr)  = iehlr ;
      endif
      
!      IF ( icefallfac /= 1.0 ) write(0,*) 'icefallfac = ',icefallfac
!      IF ( snowfallfac /= 1.0 ) write(0,*) 'snowfallfac = ',snowfallfac


  return
end subroutine nssl_2mom_init

! #####################################################################
! #####################################################################

subroutine nssl_2mom_driver(qv, qc, qr, qi, qs, qh, qhl, ccw, crw, cci, csw, chw, chl,  &
                              cn, vhw, vhl, cna, cni, f_cn, f_cna, f_cina,              &
                              f_qc, f_qr, f_qi, f_qs, f_qh, f_qhl,                      &
                              cnuf, f_cnuf,                                             &
                              zrw, zhw, zhl, f_zrw, f_zhw, f_zhl, f_vhw, f_vhl,         &
                              qsw, qhw, qhlw,                                           &
                              tt, th, pii, p, w, dn, dz, dtp, itimestep,                &
                              is_theta_or_temp,                                         &
                              ntmul, ntcnt, lastloop,                                   &
                              RAIN,nrain,prate,                                         &
                              tcond,tevac,tevar,train,                                  &
                              ruh,rvh,rmh,rr,                                           &
                              dx, dy,                                                   &
                              ib3d,ie3d,jb3d,je3d,kb3d,ke3d,                            &
                              nout3d,out3d,                                             &
                              axtra, dorad,                                             &
                              SNOWNC, SNOWNCV, GRPLNC, GRPLNCV,                         &
                              SR,HAILNC, HAILNCV,                                       &
                              hail_maxk1, hail_max2d, nwp_diagnostics,                  &
                              tkediss,                                                  &
                              re_cloud, re_ice, re_snow, re_rain,                       &
                              re_graup, re_hail,                                        &
                              has_reqc, has_reqi, has_reqs, has_reqr,                   &
                              has_reqg, has_reqh,                                       &
                              rainncw2, rainnci2,                                       &
                              dbz, vzf,compdbz,                                         &
                              rscghis_2d,rscghis_2dp,rscghis_2dn,                       &
                              scr,scw,sci,scs,sch,schl,sctot,                           &
                              elec_physics,                                             &
                              induc,elecz,scion,sciona,                                 &
                              noninduc,noninducp,noninducn,                             &
                              pcc2, pre2, depsubr,      &
                              mnucf2, melr2, ctr2,     &
                              rim1_2, rim2_2,rim3_2, &
                              nctr2, nnuccd2, nnucf2, &
                              effc2,effr2,effi2,       &
                              effs2, effg2,                       &
                              fc2, fr2,fi2,fs2,fg2, &
                              fnc2, fnr2,fni2,fns2,fng2, &
!                              qcond,qdep,qfrz,qrauto,qhcnvi,qhcollw,qscollw,            &
!                              ncauto, niinit,nifrz,                                     &
!                              re_liquid, re_graupel, re_hail, re_icesnow,               &
!                              vtcloud, vtrain, vtsnow, vtgraupel, vthail,               &
                              ipelectmp,                                                &
                              diagflag,ke_diag,                                         &
                              nssl_progn,                                              & ! wrf-chem 
! 20130903 acd_mb_washout start
                              wetscav_on, rainprod, evapprod,                           & ! wrf-chem 
! 20130903 acd_mb_washout end
                              cu_used, qrcuten, qscuten, qicuten, qccuten,              & ! hm added
                              ids,ide, jds,jde, kds,kde,                                &  ! domain dims
                              ims,ime, jms,jme, kms,kme,                                &  ! memory dims
                              its,ite, jts,jte, kts,kte)                                   ! tile dims




! this line intentionally commented
#ifdef MPI
      use mpi
! this line intentionally commented
#endif

      implicit none


 !Subroutine arguments:
      integer :: mytask = 0
      integer, intent(in), optional ::                                                   &
                            ids,ide, jds,jde, kds,kde

      integer, intent(in)::                                                             &
                            ims,ime, jms,jme, kms,kme,                                   &
                            its,ite, jts,jte, kts,kte
      real, dimension(ims:ime, jms:jme, kms:kme), intent(inout)::                        &
                            qv,qc,qr,qs,qh
      ! tt is air temperature -- used by CCPP instead of th (theta)
      real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(inout)::                       &
                              th, tt,                                                   &
                              zrw, zhw, zhl,                                            &
                              qsw, qhw, qhlw,                                           &
                            qi,qhl,ccw,crw,cci,csw,chw,chl,vhw,vhl
      integer, optional, intent(in) :: is_theta_or_temp
      logical, optional, intent(in) ::  f_zrw, f_zhw, f_zhl, f_vhw, f_vhl ! not used yet
      real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(inout):: dbz, vzf, cn, cna, cni, cnuf
      real, dimension(ims:ime, jms:jme), optional, intent(inout):: compdbz
      real, dimension(ims:ime, jms:jme), optional, intent(inout):: rscghis_2d,  & ! 2D accumulation arrays for vertically-integrated charging rate
                                                                   rscghis_2dp, & ! 2D accumulation arrays for vertically-integrated charging rate (positive only)
                                                                   rscghis_2dn    ! 2D accumulation arrays for vertically-integrated charging rate (negative only)
!      real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(inout)::rscghis_3d
      integer, optional, intent(in) :: elec_physics
      real, dimension(ims:ime, jms:jme, kms:kme),  optional, intent(inout)::                   &
                            scr,scw,sci,scs,sch,schl,sciona,sctot  ! space charge
      real, dimension(ims:ime, jms:jme, kms:kme),  optional, intent(inout)::                   &
                            induc,noninduc,noninducp,noninducn  ! charging rates: inductive, noninductive (all, positive, negative to graupel)
      real, dimension(ims:ime, jms:jme, kms:kme),  optional, intent(in) :: elecz ! elecsave = Ez
      real, dimension(ims:ime, jms:jme, kms:kme,2),optional, intent(inout) :: scion
      real, dimension(ims:ime, jms:jme, kms:kme), intent(in)::  p,w,dz,dn

      real, dimension(ims:ime, jms:jme, kms:kme), intent(inout)::  pii
      real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(inout)::   &
                              pcc2, pre2, depsubr,      &
                              mnucf2, melr2, ctr2,     &
                              rim1_2, rim2_2,rim3_2, &
                              nctr2, nnuccd2, nnucf2, &
                              effc2,effr2,effi2,       &
                              effs2, effg2,                       &
                              fc2, fr2,fi2,fs2,fg2, &
                              fnc2, fnr2,fni2,fns2,fng2
!                              qcond,qdep,qfrz,qrauto,qhcnvi,qhcollw,qscollw,            &
!                              ncauto, niinit,nifrz,                                     &
!                              re_liquid, re_graupel, re_hail, re_icesnow,               &
!                              vtcloud, vtrain, vtsnow, vtgraupel, vthail               

       real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(inout) :: axtra

! COMMAS vars
      integer, intent(in) :: nrain
      real, dimension(ims:ime, jms:jme,nrain), intent(inout)::                                 &
                            rain                          ! accumulated precip (NC) and rate (NCV)
      real, dimension(ims:ime, jms:jme), optional, intent(inout) :: prate
      real, dimension(ims:ime, jms:jme) ::                                 & ! used as local temp array in 1/COMMAS
                            RAINNC,RAINNCV          ! accumulated precip (NC) and rate (NCV)
      real, dimension(ims:ime, jms:jme), optional, intent(inout) ::                                 &
                            SNOWNC,SNOWNCV,GRPLNC,GRPLNCV,SR,HAILNC,HAILNCV         ! accumulated precip (NC) and rate (NCV)
!      logical, parameter :: f_THPROC = .false.
      
! extra vars for cm1
      real*8, intent(inout) :: tcond,tevac,tevar,train
      real, dimension(ims:ime, jms:jme, kms:kme), intent(in)::  rmh, rr
      real, intent(in), dimension(ims:ime) :: ruh
      real, intent(in), dimension(jms:jme) :: rvh
      logical, intent(in) :: dorad
      real*8, dimension(jte) :: bud1,bud2
      real, dimension(ims:ime, jms:jme), optional, intent(inout) :: hail_maxk1, hail_max2d
      integer, optional, intent(in) :: nwp_diagnostics
!     for cm1, set nproctot=44 (or as needed) to get domain total rates
      integer, parameter :: nproc = 1
      double precision :: proctot(nproc),proctotmpi(nproc)
      real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(inout)::  re_cloud, re_ice, re_snow, &
                                                                   re_rain, re_graup, re_hail
      real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(in):: tkediss
      integer, intent(in), optional :: has_reqc, has_reqi, has_reqs, has_reqr, has_reqg, has_reqh
      real, dimension(ims:ime, jms:jme), intent(out), optional ::                                 &
                            rainncw2, rainnci2       ! liquid rain, ice, accumulation rates
      real, optional, intent(in) :: dx,dy
      real, intent(in)::    dtp
      integer, intent(in):: itimestep !, ccntype
      integer, intent(in), optional :: ntmul, ntcnt
      logical, optional, intent(in) :: lastloop
      logical, optional, intent(in) :: diagflag, f_cna, f_cn, f_cina, f_cnuf
      logical, optional, intent(in) :: f_qc, f_qr, f_qi, f_qs, f_qh, f_qhl
      integer, optional, intent(in) :: ipelectmp, ke_diag


  logical, intent(in), optional ::    nssl_progn   ! flags for wrf-chem 
  
!   REAL, DIMENSION(ims:ime, jms:jme, kms:kme), optional,INTENT(INOUT):: qndrop
  logical :: flag_qndrop  ! wrf-chem
  logical :: flag_qnifa , flag_qnwfa
  logical :: flag_cnuf = .false.
  logical :: flag_ccn = .false.
  logical :: flag_qi  = .true.
  logical :: has_reqg_local = .false., has_reqh_local = .false.
  logical :: flag
  logical :: nwp_diagflag = .false.
  real :: cinchange, t7max,testmax,wmax

! 20130903 acd_ck_washout start
! rainprod - total tendency of conversion of cloud water/ice and graupel to rain (kg kg-1 s-1)
! evapprod - tendency of evaporation of rain (kg kg-1 s-1)
! 20130903 acd_ck_washout end
   real, dimension(ims:ime, jms:jme, kms:kme), optional,intent(inout)::  rainprod, evapprod

! qrcuten, rain tendency from parameterized cumulus convection
! qscuten, snow tendency from parameterized cumulus convection
! qicuten, cloud ice tendency from parameterized cumulus convection
! mu : air mass in column
   real, dimension(ims:ime, jms:jme, kms:kme), optional, intent(in):: qrcuten, qscuten, qicuten, qccuten
   integer, optional, intent(in) :: cu_used
   logical, optional, intent(in) :: wetscav_on

      integer, intent(in) :: ib3d,ie3d,jb3d,je3d,kb3d,ke3d,nout3d
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
!
! local variables
!
     real, dimension(its:ite, 1, kts:kte) :: elec2 ! ez = elecsave slab
!     real, dimension(its:ite, 1, kts:kte,2) :: scion2 ! 1=- , 2=+
     real, dimension(its:ite, kts:kte) :: rainprod2d, evapprod2d,tke2d
     real, dimension(its:ite, 1, kts:kte, na) :: an, ancuten
     real, dimension(its:ite, 1, kts:kte, nxtra) :: axtra2d
     real, dimension(its:ite, 1, kts:kte, 3) :: alpha2d
     real, dimension(its:ite, 1, kts:kte) :: t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
     real, dimension(its:ite, 1, kts:kte) :: dn1,t00,t77,ssat,pn,wn,dz2d,dz2dinv,dbz2d,vzf2d
     real, dimension(its:ite, 1, na) :: xfall
     real, dimension(its:ite, 1) :: hailmax1d,hailmaxk1
     real, dimension(kts:kte, nproc) :: thproclocal
     integer, parameter :: nor = 0, ng = 0
     integer :: nx,ny,nz,ngs
     integer ix,jy,kz,i,j,k,il,n
     integer :: infdo
     real :: ssival, ssifac, t8s, t9s, qvapor
     integer :: ltemq
     double precision :: dp1
     integer :: jye, lnb
     integer :: imx,kmx
     real    :: dbzmx,refl
     integer :: vzflag0 = 0
     logical :: makediag
     real    :: dx1,dy1
      real, parameter :: cnin20 = 1.0e3
      real, parameter :: cnin10 = 5.0e1
      real, parameter :: cnin1a = 4.5
      real, parameter :: cnin2a = 12.96
      real, parameter :: cnin2b = 0.639

      double precision :: cwmass1,cwmass2
      double precision :: rwmass1,rwmass2
      double precision :: icemass1,icemass2
      double precision :: swmass1,swmass2
      double precision :: grmass1,grmass2
      double precision :: hlmass1,hlmass2
      double precision :: wvol5,wvol10
      real :: tmp,dv,dv1,tmpchg
      real :: rdt
      
      double precision :: dt1,dt2
      double precision :: timesed,timesed1,timesed2,timesed3, timegs, timenucond, timedbz,zmaxsed
      double precision :: timevtcalc,timesetvt
      
      logical :: f_cnatmp, f_cinatmp
      logical :: has_wetscav

      integer :: kediagloc
      integer :: iunit

      real :: ycent, y, emissrate, emissrate0, emissrate1, z, fac, factot
      real :: fach(kts:kte)
      
      logical, parameter :: debugdriver = .false.
      
      integer :: loopcnt, loopmax, outerloopcnt
      logical :: lastlooptmp

!
#if defined(MPI) 
      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
      integer :: mpi_error_code = 1
#endif

! -------------------------------------------------------------------


      rdt = 1.0/dtp
      
     if ( debugdriver ) write(0,*) 'N2M: entering routine'

     flag_qndrop = .false.
     flag_qnifa = .false.
     flag_qnwfa = .false.
     flag_cnuf = .false.
     flag_ccn = .false.
     nwp_diagflag = .false.
     
     if ( PRESENT ( nssl_progn ) ) flag_qndrop = nssl_progn
     if ( present ( f_cnuf ) ) flag_cnuf = f_cnuf
     if ( present ( nwp_diagnostics ) ) nwp_diagflag = ( nwp_diagnostics > 0 )

     if ( present ( f_cn ) .and. present( cn ) ) then 
       flag_ccn = f_cn
     elseif ( present( cn ) ) then
       flag_ccn = .true.
     endif
     
     if ( present( f_qi ) ) then
       flag_qi = f_qi
     else
       if ( ffrzs < 1.0 ) then
         flag_qi = .true.
       else
         flag_qi = .false.
       endif
     endif
     
     if ( .not. flag_qi .and. ffrzs < 1.0 ) ffrzs = 1.0

     
     if ( PRESENT ( has_reqg ) ) has_reqg_local = has_reqg > 0
     if ( PRESENT ( has_reqh ) ) has_reqh_local = has_reqh > 0
     
     loopmax = 1
     outerloopcnt = 1
     lastlooptmp = .true.
     if ( present( ntmul ) .and. present( ntcnt ) .and. present( lastloop ) ) then
       loopmax = ntmul
       outerloopcnt = ntcnt
       lastlooptmp = lastloop
     endif
          

         has_wetscav = .false.
         if ( wrfchem_flag > 0 ) then
           if ( PRESENT( wetscav_on ) ) then
             has_wetscav = wetscav_on
           endif
         endif

     if ( present( f_cna ) ) then
       f_cnatmp = f_cna
     else 
       f_cnatmp = .false.
     endif

     if ( present( f_cina ) ) then
       f_cinatmp = f_cina
     else 
       f_cinatmp = .false.
     endif
       
     if ( present( vzf ) ) vzflag0 = 1
     
     if ( present( ipelectmp ) ) then
       ipelec = ipelectmp
     else
       ipelec = 0
     endif
!       IF ( present( dbz ) ) THEN
!       DO jy = jts,jte
!         DO kz = kts,kte
!           DO ix = its,ite
!             dbz(ix,jy,kz) = 0.0
!           ENDDO
!         ENDDO
!       ENDDO
!       ENDIF

     if ( present( dx ) .and. present( dy ) ) then
        dx1 = dx
        dy1 = dy
     else
        dx1 = 1.0
        dy1 = 1.0
     endif

     
     makediag = .true.
     if ( present( diagflag ) ) then
      makediag = diagflag .or. itimestep == 1
      if ( diagflag ) cleardiag = .true.
     endif

     if ( debugdriver ) write(0,*) 'N2M: makediag = ',makediag
     
     
     nx = ite-its+1
     ny = 1         ! set up as 2D slabs
     nz = kte-kts+1
     ngs = 64
     
     if ( .not. flag_ccn ) then
       renucfrac = 1.0
     endif

     


!     ENDIF ! itimestep == 1


! sedimentation settings

      infdo = 2
      
      if ( infall .ne. 1 .or. iscfall .ge. 2 ) then
         infdo = 1
      else
         infdo = 0
      endif

      if ( infall .ge. 3 .or. ipconc .ge. 6 ) then
         infdo = 2
      endif
     

      if (present( HAILNCV ) ) HAILNCV(its:ite,jts:jte) = 0.

      tke2d(:,:) = 0.0 ! initialize if not used

     lnb = Max(lh,lhl)+1 ! lnc
!     IF ( lccn > 1 ) lnb = lccn

       jye = jte

     if ( present( compdbz ) .and. makediag ) then
     do jy = jts,jye
       do ix = its,ite
        compdbz(ix,jy) = -3.0
       enddo
     enddo
     endif

      zmaxsed = 0.0d0
      timevtcalc = 0.0d0
      timesetvt = 0.0d0
      timesed = 0.0d0
      timesed1 = 0.0d0
      timesed2 = 0.0d0
      timesed3 = 0.0d0
      timegs = 0.0d0
      timenucond = 0.0d0



     if ( debugdriver ) write(0,*) 'N2M: jy loop 1, lhl,na = ',lhl,na,present(qhl)

          ancuten(its:ite,1,kts:kte,:) = 0.0
          thproclocal(:,:) = 0.0
          proctot(:) = 0

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ix,jy,kz,xfall,axtra2d,an,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9, &
!$OMP t00,t77,dbz2d,vzf2d,dn1,pn,wn,dz2d,dz2dinv,ltemq,ssival,dp1,    &
!$OMP elec2,thproclocal,t8s,t9s,il,n,ssat,loopcnt,   &
!$OMP qvapor,ssifac,ancuten,rainprod2d,evapprod2d,kediagloc,alpha2d, &
!$OMP hailmax1d,hailmaxk1)

! need to figure out how to treat thproc: with a critical section?

     do jy = jts,jye
     
!     write(0,*) 'N2M: load an, jy,lccn = ',jy,lccn,qccn

     if ( present( pcc2 ) .and. makediag ) then
         axtra2d(its:ite,1,kts:kte,:) = 0.0
     endif

     if ( nwp_diagflag ) then
        alpha2d(its:ite,1,kts:kte,1) = alphar
        alpha2d(its:ite,1,kts:kte,2) = alphah
        alpha2d(its:ite,1,kts:kte,3) = alphahl
     endif


! for now, only using thproclocal for domain totals
       thproclocal(:,:) = 0.0
   ! copy from 3D array to 2D slab
   
       do kz = kts,kte
        do ix = its,ite
          an(ix,1,kz,lt) = th(ix,jy,kz)
          an(ix,1,kz,lv)   = qv(ix,jy,kz)
          an(ix,1,kz,lc)   = qc(ix,jy,kz)
          an(ix,1,kz,lr)   = qr(ix,jy,kz)
          if ( flag_qi ) then
            an(ix,1,kz,li)   = qi(ix,jy,kz)
          else
            an(ix,1,kz,li) = 0.0
          endif
          an(ix,1,kz,ls)   = qs(ix,jy,kz)
          an(ix,1,kz,lh)   = qh(ix,jy,kz)
          if ( lhl > 1 ) an(ix,1,kz,lhl)  = qhl(ix,jy,kz)
          if ( lccn > 1 ) then
           if ( is_aerosol_aware .and. flag_qnwfa ) then
            ! 
           elseif ( flag_ccn ) then
             if ( lccna > 1 .and. .not. ( present( cna ) .and. f_cnatmp ) ) then
               an(ix,1,kz,lccna) = cn(ix,jy,kz)
               an(ix,1,kz,lccn) = qccn ! cn(ix,jy,kz)
             else
               an(ix,1,kz,lccn) = cn(ix,jy,kz)
             endif
             if ( i_uf_or_ccn > 0 .and. lccnuf > 1 ) then ! UF ccn are extra regular ccn
               an(ix,1,kz,lccn) = an(ix,1,kz,lccn) + cnuf(ix,jy,kz)
             endif
           else
            if ( lccna == 0 .and. ( .not. f_cnatmp ) ) then
              an(ix,1,kz,lccn) = qccn - ccw(ix,jy,kz)
            else
              an(ix,1,kz,lccn) = qccn 
            endif
           
           endif
          endif

          if ( lccnuf > 0 .and. flag_cnuf ) then
            if ( i_uf_or_ccn == 0 ) then ! UF are UF
              an(ix,1,kz,lccnuf) = Max(0.0, cnuf(ix,jy,kz) )
            else ! UF were added to lccn
              an(ix,1,kz,lccnuf) = 0.0
            endif
          endif

          if ( lccna > 1 ) then
            if ( present( cna ) .and. f_cnatmp ) then
              an(ix,1,kz,lccna) = cna(ix,jy,kz)
            endif
          endif

          if ( lcina > 1 ) then
            if ( present( cni ) .and. f_cinatmp ) then
              an(ix,1,kz,lcina) = cni(ix,jy,kz)
            endif
          endif
          
          if ( ipconc >= 5 ) then
             an(ix,1,kz,lnc)  = ccw(ix,jy,kz)
          if ( constccw > 0.0 ) then
            an(ix,1,kz,lnc)  = constccw
          endif
          an(ix,1,kz,lnr)  = crw(ix,jy,kz)
          if ( present( cci ) ) then
            an(ix,1,kz,lni)  = cci(ix,jy,kz)
          else
            an(ix,1,kz,lni) = 0.0
          endif
          an(ix,1,kz,lns)  = csw(ix,jy,kz)
          an(ix,1,kz,lnh)  = chw(ix,jy,kz)
          if ( lhl > 1 ) an(ix,1,kz,lnhl) = chl(ix,jy,kz)
          endif
          if ( lvh > 0 ) an(ix,1,kz,lvh)  = vhw(ix,jy,kz)
          if ( lvhl > 0 .and. present( vhl ) ) an(ix,1,kz,lvhl)  = vhl(ix,jy,kz)

          if ( ipconc >= 6 ) then
            if ( lzr > 0 )  an(ix,1,kz,lzr)  = zrw(ix,jy,kz)*zscale
            if ( lzh > 0 )  an(ix,1,kz,lzh)  = zhw(ix,jy,kz)*zscale
            if ( lzhl > 0 ) an(ix,1,kz,lzhl) = zhl(ix,jy,kz)*zscale
          endif
          


        enddo
       enddo
       
       do kz = kts,kte
        do ix = its,ite

          
          t0(ix,1,kz) = th(ix,jy,kz)*pii(ix,jy,kz) ! temperature (Kelvin)
          t00(ix,1,kz) = 380.0/p(ix,jy,kz)
          t77(ix,1,kz) = pii(ix,jy,kz)
          dbz2d(ix,1,kz) = 0.0
          vzf2d(ix,1,kz) = 0.0
        enddo
       enddo
       
       do ix = its,ite
         RAINNCV(ix,jy) = 0.0
         if ( present( GRPLNCV ) ) GRPLNCV(ix,jy) = 0.0
         if ( present( HAILNCV ) ) HAILNCV(ix,jy) = 0.0
         if ( present( SNOWNCV ) ) SNOWNCV(ix,jy) = 0.0
       enddo

      do loopcnt = 1,loopmax
       
       do kz = kts,kte
        do ix = its,ite

          
          t1(ix,1,kz) = 0.0
          t2(ix,1,kz) = 0.0
          t3(ix,1,kz) = 0.0
          t4(ix,1,kz) = 0.0
          t5(ix,1,kz) = 0.0
          t6(ix,1,kz) = 0.0
          t7(ix,1,kz) = 0.0
          t8(ix,1,kz) = 0.0
          t9(ix,1,kz) = 0.0

          pn(ix,1,kz) = p(ix,jy,kz)
          wn(ix,1,kz) = w(ix,jy,kz)
          dn1(ix,1,kz) = dn(ix,jy,kz)
!          wmax = Max(wmax,wn(ix,1,kz))
          dz2d(ix,1,kz) = dz(ix,jy,kz)
          dz2dinv(ix,1,kz) = 1./dz(ix,jy,kz)
          
         ltemq = Int( (t0(ix,1,kz)-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )
!
! saturation mixing ratio
!
      t8s = t00(ix,1,kz)*tabqvs(ltemq)  !saturation mixing ratio wrt water
      t9s = t00(ix,1,kz)*tabqis(ltemq)  !saturation mixing ratio wrt ice

!
!  calculate rate of nucleation
!
      ssival = Min(t8s,max(an(ix,1,kz,lv),0.0))/t9s  ! qv/qvi


      if ( ssival .gt. 1.0 ) then
!
      if ( icenucopt == 1 ) then

      if ( t0(ix,1,kz).le.268.15 ) then

       dp1 = dn1(ix,1,kz)/rho00*cnin20*exp( Min( 57.0 ,(cnin2a*(ssival-1.0)-cnin2b) ) )
       t7(ix,1,kz) = Min(dp1, 1.0d30)
      end if

!
!   Default value of imeyers5 turns off nucleation by Meyer at higher temperatures
!  This is really from Ferrier (1994), eq. 4.31 - 4.34
      if ( imeyers5 ) then
      if ( t0(ix,1,kz).lt.tfr .and. t0(ix,1,kz).gt.268.15 ) then
      qvapor = max(an(ix,1,kz,lv),0.0)
      ssifac = 0.0
      if ( (qvapor-t9s) .gt. 1.0e-5 ) then
      if ( (t8s-t9s) .gt. 1.0e-5 ) then
      ssifac = (qvapor-t9s) /(t8s-t9s)
      ssifac = ssifac**cnin1a
      end if
      end if
      t7(ix,1,kz) = dn1(ix,1,kz)/rho00*cnin10*ssifac*exp(-(t0(ix,1,kz)-tfr)*bta1)
      end if
      endif
      
!       t7max = Max(t7max,  t7(ix,1,kz) )

      elseif ( icenucopt == 2 ) then ! Thompson/Cooper; Note Thompson 2004 has constants of
                                     ! 0.005 and 0.304 because the line function was estimated from Cooper plot
                                     ! Here, the fit line values from Cooper 1986 are converted. Very little difference 
                                     ! in practice
      
        t7(ix,1,kz) = 1000.*0.00446684*exp(0.3108*(273.16 - Max(233.0, t0(ix,1,kz) ) ) ) ! factor of 1000 to convert L**-1 to m**-3

!        write(0,*) 'Cooper t7,ssival = ',ix,kz,t7(ix,1,kz),ssival
      
      elseif ( icenucopt == 3 ) then ! Phillips (Meyers/DeMott)

      if ( t0(ix,1,kz).le.268.15 .and.  t0(ix,1,kz) > 243.15 ) then ! Meyers with factor of Psi=0.06
        
       dp1 = 0.06*cnin20*exp( Min( 57.0 ,(cnin2a*(ssival-1.0)-cnin2b) ) )
       t7(ix,1,kz) = Min(dp1, 1.0d30)
      elseif ( t0(ix,1,kz) <= 243.15 ) then ! Phillips estimate of DeMott et al (2003) data
       dp1 = 1000.*( exp( Min( 57.0 ,cnin2a*(ssival-1.1) ) ) )**0.3
       t7(ix,1,kz) = Min(dp1, 1.0d30)
      
      end if

      elseif ( icenucopt == 4 ) then ! DeMott 2010

        if ( t0(ix,1,kz) < 268.16 .and.  t0(ix,1,kz) > 223.15 .and. ssival > 1.001 ) then ! 
      
        ! a = 0.0000594, b = 3.33, c = 0.0264, d = 0.0033,
        ! nint = a*(-Tc)**b * naer**(c*(-Tc) + d)
        ! nint has units of per (standard) liter, so mult by 1.e3 and scale by dn/rho00
        ! naer needs units of cm**-3, so mult by 1.e-6
        
        !  dp1 = 1.e3*0.0000594*(273.16 - t0(ix,1,kz))**3.33 * (1.e-6*cin*dn(ix,1,kz))**(0.0264*(273.16 - t0(ix,1,kz)) + 0.0033)
            tmp = 1.e-6*naer
          dp1 = 1.e3*dn1(ix,1,kz)/rho00*0.0000594*(273.16 - t0(ix,1,kz))**3.33 * tmp**(0.0264*(273.16 - t0(ix,1,kz)) + 0.0033)
          t7(ix,1,kz) = Min(dp1, 1.0d30)
      
        else
       !   t7(ix,1,kz) = 0.0
        endif
      
      endif ! icenucopt


!
      end if ! ( ssival .gt. 1.0 )
!

        enddo ! ix
       enddo ! kz

      if ( wrfchem_flag > 0 ) then
          if ( has_wetscav ) then
            if ( PRESENT( rainprod ) ) rainprod2d(its:ite,kts:kte) = 0
            if ( PRESENT( evapprod ) ) evapprod2d(its:ite,kts:kte) = 0
          endif
      endif
         

   ! transform from number mixing ratios to number conc.
     
    if ( loopcnt == 1 ) then
     do il = lnb,na
       if ( denscale(il) == 1 ) then
         do kz = kts,kte
          do ix = its,ite
           an(ix,1,kz,il) = an(ix,1,kz,il)*dn1(ix,1,kz) ! dn(ix,jy,kz)
          enddo
         enddo
       endif
     enddo ! il
    endif

        
! sedimentation
      xfall(:,:,:) = 0.0
       

!      IF ( .true. ) THEN


! #ifndef 1
! for real cases when hydrometeor mixing ratios have been initialized without concentrations
       if ( itimestep == 1 .and. ipconc > 0 .and. loopcnt == 1 ) then
         call calcnfromq(nx,ny,nz,an,na,nor,nor,dn1)
       endif
! #endif

      if ( present(cu_used) .and.         &
           ( present( qrcuten ) .or. present( qscuten ) .or.  &
             present( qicuten ) .or. present( qccuten ) ) ) then !{

       if ( cu_used == 1 ) then !{
       do kz = kts,kte
        do ix = its,ite

         if ( present( qrcuten ) ) ancuten(ix,1,kz,lr) = dtp*qrcuten(ix,jy,kz)
         if ( present( qscuten ) ) ancuten(ix,1,kz,ls) = dtp*qscuten(ix,jy,kz)
         if ( present( qicuten ) ) ancuten(ix,1,kz,li) = dtp*qicuten(ix,jy,kz)
         if ( present( qccuten ) ) ancuten(ix,1,kz,lc) = dtp*qccuten(ix,jy,kz)
         
        enddo
       enddo
       
         call calcnfromcuten(nx,ny,nz,ancuten,an,na,nor,nor,dn1)

       do kz = kts,kte
        do ix = its,ite


          if ( ipconc >= 6 ) then
!            IF ( lzr > 0 )  an(ix,1,kz,lzr)  = an(ix,1,kz,lzr) + ancuten(ix,1,kz,lzr) 
          endif
         
        enddo
       enddo
       
       endif !}
       
      endif !}
      
      


      call sediment1d(dtp,nx,ny,nz,an,na,nor,nor,xfall,dn1,dz2d,dz2dinv, &
     &                    t0,t7,infdo,jy,its,jts &
     &   ,timesed1,timesed2,timesed3,zmaxsed,timesetvt)


! copy xfall to appropriate places...

     if ( debugdriver ) write(0,*) 'N2M: end sediment, jy = ',jy

       do ix = its,ite
         if ( lhl > 1 ) then
         RAINNCV(ix,jy) = RAINNCV(ix,jy) + &
                          dtp*dn1(ix,1,1)*(xfall(ix,1,lr) + xfall(ix,1,ls)*1000./xdn0(lr) + &
              &            xfall(ix,1,lh)*1000./xdn0(lr) + xfall(ix,1,lhl)*1000./xdn0(lr) )
         else
         RAINNCV(ix,jy) = RAINNCV(ix,jy) + &
                           dtp*dn1(ix,1,1)*(xfall(ix,1,lr) + xfall(ix,1,ls)*1000./xdn0(lr) + &
              &            xfall(ix,1,lh)*1000./xdn0(lr) )
         endif
         if ( present ( rainncw2 ) ) then ! rain only
           rainncw2(ix,jy) =  rainncw2(ix,jy) +  dtp*dn1(ix,1,1)*xfall(ix,1,lr)
         endif
         if ( present ( rainnci2 ) ) then ! ice only
           if ( lhl > 1 ) then
             rainnci2(ix,jy) =rainnci2(ix,jy) + dtp*dn1(ix,1,1)*(xfall(ix,1,ls)*1000./xdn0(lr) + &
     &            xfall(ix,1,lh)*1000./xdn0(lr) + xfall(ix,1,lhl)*1000./xdn0(lr) )
           else
             rainnci2(ix,jy) = rainnci2(ix,jy) + dtp*dn1(ix,1,1)*(xfall(ix,1,ls)*1000./xdn0(lr) + &
     &            xfall(ix,1,lh)*1000./xdn0(lr) )
           endif
         endif
         if ( present( SNOWNCV ) ) SNOWNCV(ix,jy) = SNOWNCV(ix,jy) + dtp*dn1(ix,1,1)*xfall(ix,1,ls)*1000./xdn0(lr)
         if ( present( GRPLNCV ) ) then 
           if ( lhl > 1 .and. .not. present( HAILNC) ) then ! if no separate hail accum, then add to graupel
             GRPLNCV(ix,jy) = GRPLNCV(ix,jy) + dtp*dn1(ix,1,1)*(xfall(ix,1,lh) + xfall(ix,1,lhl)) *1000./xdn0(lr)
           else
             GRPLNCV(ix,jy) = GRPLNCV(ix,jy) + dtp*dn1(ix,1,1)*xfall(ix,1,lh)*1000./xdn0(lr)
           endif
         endif
         if ( loopcnt == loopmax ) RAINNC(ix,jy)  = RAINNC(ix,jy) + RAINNCV(ix,jy)

         if ( loopcnt == loopmax ) then
         do n = 1,nrain
           RAIN(ix,jy,n) = RAIN(ix,jy,n) + 0.1*RAINNCV(ix,jy)
         enddo
         endif

         if ( present( prate ) ) prate(ix,jy) = rdt*RAINNCV(ix,jy)

!         bud1(jy) = bud1(jy) + RAINNCV(ix,jy)*ruh(ix)*rvh(jy)*dx*dy
!  Moved train sum to end of jy loop
!         do k=1,kte
!           bud2(jy) = bud2(jy) - rr(ix,jy,kz)*PRE(kz)*ruh(ix)*rvh(jy)*rmh(ix,jy,kz)*tem1
!         enddo

         if ( present (SNOWNC) .and. present (SNOWNCV) .and. loopcnt == loopmax ) then
           SNOWNC(ix,jy)  = SNOWNC(ix,jy) + SNOWNCV(ix,jy)
         endif
         if ( lhl > 1 ) then
!#ifdef 1
!           IF ( .true. ) THEN
!#else
           if ( present( HAILNC ) ) then
!#endif
             HAILNCV(ix,jy) = dtp*dn1(ix,1,1)*xfall(ix,1,lhl)*1000./xdn0(lr)
             if ( loopcnt == loopmax ) HAILNC(ix,jy)  = HAILNC(ix,jy) + HAILNCV(ix,jy)
!           ELSEIF ( present( GRPLNCV ) ) THEN ! if no separate hail accum, then add to graupel
!             GRPLNCV(ix,jy) = GRPLNCV(ix,jy) + dtp*dn1(ix,1,1)*xfall(ix,1,lhl)*1000./xdn0(lr)
           endif
         endif
         if ( present( GRPLNCV ) .and. loopcnt == loopmax ) then
           GRPLNC(ix,jy)  = GRPLNC(ix,jy) + GRPLNCV(ix,jy)
         endif
       enddo
       
!      ENDIF ! .false.
 
      if ( isedonly /= 1 ) then
   ! call nssl_2mom_gs: main gather-scatter routine to calculate microphysics

     if ( debugdriver ) write(0,*) 'N2M: gs, jy = ',jy
!      IF ( isedonly /= 2 ) THEN


      call nssl_2mom_gs   &
     &  (nx,ny,nz,na,jy   &
     &  ,nor,nor          &
     &  ,dtp,dz2d       &
     &  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9     &
     &  ,an,dn1,t77                        &
     &  ,pn,wn,0                           &
     &  ,t00,t77,                          &
     &   ventr,ventc,c1sw,1,ido,           &
     &   xdnmx,xdnmn,                      &
!     &   ln,ipc,lvol,lz,lliq,              &
     &   cdx,                              &
     &   xdn0,dbz2d,tke2d,                 &
     &   thproclocal,nproc,dx1,dy1,ngs,    &
     &   timevtcalc,axtra2d, makediag        &
     &   ,has_wetscav, rainprod2d, evapprod2d, alpha2d  &
     &   ,elec2,its,ids,ide,jds,jde          &
     & )





   endif ! isedonly /= 1
   
 ! droplet nucleation/condensation/evaporation
   if ( .true. ) then
   call NUCOND    &
     &  (nx,ny,nz,na,jy & 
     &  ,nor,nor,dtp,nx  &
     &  ,dz2d & 
     &  ,t0,t9 & 
     &  ,an,dn1,t77 & 
     &  ,pn,wn &
     &  ,ngs   &
     &  ,axtra2d, makediag  &
     &  ,ssat,t00,t77,flag_qndrop)



   endif




     enddo ! loopcnt=1,loopmax
     if ( present( pcc2 ) .and. makediag ) then
         do kz = kts,kte
          do ix = its,ite
! example of using the 'axtra2d' array to get rates out of the microphysics routine for output.
! Search for 'axtra' to find example code below
!            pcc2(ix,jy,kz)    = axtra2d(ix,1,kz,1)
          enddo
         enddo
     endif


! compute diagnostic S-band reflectivity if needed
     if ( present( dbz ) .and. makediag .and. lastlooptmp ) then
   ! calc dbz
      
      if ( .true. ) then
      if ( present(ke_diag) ) then
        kediagloc = ke_diag
      else
        kediagloc = nz
      endif
      call radardd02(nx,ny,nz,nor,na,an,t0,         &
     &    dbz2d,dn1,nz,cnoh,rho_qh,ipconc,kediagloc, 0)
      endif ! .false.

     
       do kz = kts,kediagloc ! kte
        do ix = its,ite
         dbz(ix,jy,kz) = dbz2d(ix,1,kz)
         if ( present( vzf ) ) then
           vzf(ix,jy,kz) = vzf2d(ix,1,kz)
           if ( dbz2d(ix,1,kz) <= 0.0 ) then
             vzf(ix,jy,kz) = 0.0
           elseif ( dbz2d(ix,1,kz) <= 15.0 ) then
             refl = 10**(0.1*dbz2d(ix,1,kz))
             vzf(ix,jy,kz) = Min( vzf2d(ix,1,kz), 2.6 * Max(0.0,refl)**0.107 * (1.2/dn1(ix,1,kz))**0.4 )
           endif
         endif
          if ( present( compdbz ) ) then
            compdbz(ix,jy) = Max( compdbz(ix,jy), dbz2d(ix,1,kz) )
          endif
        enddo
       enddo

       endif



! Following Greg Thompson, calculation for effective radii. Used by RRTMG LW/SW schemes if enabled in module_physics_init.F
      if ( dorad ) then
      if ( present( has_reqc ).and. present( has_reqi ) .and. present( has_reqs ) .and.  &
           present( re_cloud ).and. present( re_ice ) .and. present( re_snow )    .and.  &
           lastlooptmp) then
       if ( has_reqc.ne.0 .or. has_reqi.ne.0 .or. has_reqs.ne.0) then
         do kz = kts,kte
          do ix = its,ite
             re_cloud(ix,jy,kz)  = 2.51E-6
             re_ice(ix,jy,kz)    = 10.01E-6
             re_snow(ix,jy,kz)   = 25.E-6
             t1(ix,1,kz) = 2.51E-6
             t2(ix,1,kz) = 10.01E-6
             t3(ix,1,kz) = 25.E-6
             t4(ix,1,kz) = 50.e-6
          enddo
         enddo


          call calc_eff_radius   &
     &         (nx,ny,nz,na,jy & 
     &          ,nor,nor & 
     &          ,t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,f_t5=has_reqg_local, f_t6=has_reqh_local  & 
     &          ,an=an,dn=dn1 )

        do kz = kts,kte
          do ix = its,ite
             re_cloud(ix,jy,kz) = MAX(2.51E-6, MIN(t1(ix,1,kz), 50.E-6))
             re_ice(ix,jy,kz)   = MAX(10.01E-6, MIN(t2(ix,1,kz), 125.E-6))
             re_snow(ix,jy,kz)  = MAX(25.E-6, MIN(t3(ix,1,kz), 999.E-6))
             ! check for case where snow needs to be treated as cloud ice (for rrtmg radiation)
             if ( .not. present(qi) ) re_ice(ix,jy,kz)  = MAX(10.E-6, MIN(t3(ix,1,kz), 125.E-6))
          enddo
         enddo

       if ( present(has_reqr) .and. present( re_rain ) ) then
         if ( has_reqr /= 0 ) then
         do kz = kts,kte
           do ix = its,ite
            re_rain(ix,jy,kz)  = MAX(50.E-6, MIN(t4(ix,1,kz), 2999.E-6))
           enddo
         enddo
         endif
       endif

       if ( present(has_reqg) .and. present( re_graup ) ) then
         if ( has_reqg /= 0 ) then
         do kz = kts,kte
           do ix = its,ite
            re_graup(ix,jy,kz)  = MAX(50.E-6, MIN(t5(ix,1,kz), 10.E-3))
           enddo
         enddo
         endif
       endif

       if ( present(has_reqh) .and. present( re_hail ) ) then
         if ( has_reqh /= 0 ) then
         do kz = kts,kte
           do ix = its,ite
            re_hail(ix,jy,kz)  = MAX(50.E-6, MIN(t5(ix,1,kz), 40.E-3))
           enddo
         enddo
         endif
       endif
       
         endif
        endif
      endif ! dorad 


     if ( present( hail_maxk1 ) .and. present( hail_max2d ) .and. nwp_diagflag ) then
         do ix = its,ite
            hailmax1d(ix,1) = hail_max2d(ix,jy)
            hailmaxk1(ix,1) = hail_maxk1(ix,jy)
         enddo

         call hailmaxd(dtp,nx,ny,nz,an,na,nor,nor,alpha2d,dn1,   &
                          hailmax1d,hailmaxk1,1 )

         do ix = its,ite
           hail_max2d(ix,jy) = hailmax1d(ix,1)
           hail_maxk1(ix,jy) = hailmaxk1(ix,1)
         enddo
!       ENDIF
     endif
   
! transform concentrations back to mixing ratios
     do il = lnb,na
      if ( denscale(il) == 1 ) then
       do kz = kts,kte
        do ix = its,ite
         an(ix,1,kz,il) = an(ix,1,kz,il)/dn1(ix,1,kz) ! dn(ix,jy,kz)
        enddo
       enddo
      endif
     enddo ! il
   
   ! copy 2D slabs back to 3D

   
       do kz = kts,kte
        do ix = its,ite
        
         th(ix,jy,kz)  = an(ix,1,kz,lt)

         if ( eqtset > 2 ) then
          pii(ix,jy,kz) = t77(ix,1,kz)
         endif
         qv(ix,jy,kz)  = an(ix,1,kz,lv)
         qc(ix,jy,kz)  = an(ix,1,kz,lc)
         qr(ix,jy,kz)  = an(ix,1,kz,lr)
         if ( flag_qi ) qi(ix,jy,kz)  = an(ix,1,kz,li)
         qs(ix,jy,kz)  = an(ix,1,kz,ls)
         qh(ix,jy,kz)  = an(ix,1,kz,lh)
         if ( lhl > 1 ) qhl(ix,jy,kz) = an(ix,1,kz,lhl)
         
         if ( lccn > 1 .and. is_aerosol_aware .and. flag_qnwfa ) then
           ! not used here
         elseif ( flag_ccn .and. lccn > 1 .and. .not. flag_qndrop) then
            if ( lccna > 1 .and. .not. ( present( cna ) .and. f_cnatmp ) ) then
              cn(ix,jy,kz) = Max(0.0, an(ix,1,kz,lccna) )
            else
              cn(ix,jy,kz) = an(ix,1,kz,lccn)
            endif
         endif
         if ( lccna > 1 ) then
           if ( present( cna ) .and. f_cnatmp ) then
              cna(ix,jy,kz) = Max(0.0, an(ix,1,kz,lccna) )
           endif
         endif

         if ( lcina > 1 ) then
           if ( present( cni ) .and. f_cinatmp ) then
              cni(ix,jy,kz) = Max(0.0, an(ix,1,kz,lcina) )
           endif
         endif

         if ( lccnuf > 0 .and. flag_cnuf ) then
           if ( i_uf_or_ccn > 0 ) then ! UF are ccn and lccnuf is zero, so put cnuf into lccnuf to do decay
             an(ix,1,kz,lccnuf) = Max(0.0, cnuf(ix,jy,kz) )
           endif
           if ( decayufccn ) then
             if ( an(ix,1,kz,lccnuf) > ufbackground ) then
               an(ix,1,kz,lccnuf) = an(ix,1,kz,lccnuf) - (an(ix,1,kz,lccnuf) - &
                            ufbackground)*(1.0 - exp(-dtp/ufccntimeconst))
             endif
           endif
           cnuf(ix,jy,kz) = an(ix,1,kz,lccnuf)
         endif



         if ( ipconc >= 5 ) then

          ccw(ix,jy,kz) = an(ix,1,kz,lnc)
          crw(ix,jy,kz) = an(ix,1,kz,lnr)
          if ( present( cci ) ) cci(ix,jy,kz) = an(ix,1,kz,lni)
          csw(ix,jy,kz) = an(ix,1,kz,lns)
          chw(ix,jy,kz) = an(ix,1,kz,lnh)
          if ( lhl > 1 ) chl(ix,jy,kz) = an(ix,1,kz,lnhl)
         endif

         if ( ipconc >= 6 ) then
            if ( lzr > 0 )  zrw(ix,jy,kz) = an(ix,1,kz,lzr) *zscaleinv
            if ( lzh > 0 )  zhw(ix,jy,kz) = an(ix,1,kz,lzh) *zscaleinv
            if ( lzhl > 0 ) zhl(ix,jy,kz) = an(ix,1,kz,lzhl)*zscaleinv
         endif



         if ( lvh > 0 )  vhw(ix,jy,kz) = an(ix,1,kz,lvh)
         if ( lvhl > 0 .and. present( vhl ) ) vhl(ix,jy,kz) = an(ix,1,kz,lvhl)

#if ( WRF_CHEM == 1 )
         if ( has_wetscav ) then
           if ( loopmax > 1 ) then
             ! wrferror not supported
           endif
           if ( PRESENT( rainprod ) ) rainprod(ix,jy,kz) = rainprod2d(ix,kz)
           if ( PRESENT( evapprod ) ) evapprod(ix,jy,kz) = evapprod2d(ix,kz)
         endif
#endif

        enddo
       enddo

      if ( nproc > 1 ) then
!$OMP CRITICAL      

       do k = kms,kme
        do n = 1,nproc
          proctot(n) = proctot(n) + thproclocal(k,n)
        enddo
       enddo

!$OMP END CRITICAL
      endif

     enddo ! jy
     
       do jy = jts,jte
       do ix = its,ite
         train = train + RAINNCV(ix,jy)*ruh(ix)*rvh(jy)*dx*dy
       enddo
       enddo
     

     




 11   format('processes: ',100(1x,1pe12.5))
      if ( nproc > 1 ) then
#ifdef MPI
      call MPI_Reduce(proctot, proctotmpi, nproc, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_error_code)
       if ( myrankcm1 == 0 ) then
        proctot(1:nproc) = proctotmpi(1:nproc)
       endif
#endif
       if ( myrankcm1 == 0 ) then
         write(6,11) (proctot(i), i=1,nproc)
       endif

      endif

  return
end subroutine nssl_2mom_driver

! #####################################################################
! #####################################################################

      real function GAMMA_SP(xx)

      implicit none
      real xx
      integer j

! Double precision ser,stp,tmp,x,y,cof(6)

      real*8 ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d+0,  &
     &            -86.50532032941677d0,   &
     &             24.01409824083091d0,   &
     &             -1.231739572450155d0,  &
     &              0.1208650973866179d-2,&
     &             -0.5395239384953d-5,   &
     &              2.5066282746310005d0/

      if ( xx <= 0.0 ) then
        write(0,*) 'Argument to gamma must be > 0!! xx = ',xx
        stop
      endif
      
      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x + 0.5d0)*Log(tmp) - tmp
      ser = 1.000000000190015d0
      do j=1,6
        y = y + 1.0d0
        ser = ser + cof(j)/y
      end do
      gamma_sp = Exp(tmp + log(stp*ser/x))

      return
      end function GAMMA_SP

! #####################################################################

      double precision function GAMMA_DPR(x)
      ! dp gamma with real input
        implicit none
        real :: x
        double precision :: xx
        
        xx = x
        
        gamma_dpr = gamma_dp(xx)
        
        return
        end function GAMMA_DPR
        



! #####################################################################

        real function GAMXINF(A1,X1)

!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                from x to infinity
!       Input :  a   --- Parameter ( a  170 )
!                x   --- Argument 
!       Output:  GIM --- gamma(a,x) t=x,Infinity
!       Routine called: GAMMA for computing gamma(x)
!       ===================================================

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        real :: a1,x1
        double precision :: xam,dlog,s,r,ga,t0,a,x
        integer :: k
        double precision :: gin, gim
        
        a = a1
        x = x1
        if ( x1 <= 0.0 ) then
           gamxinf = GAMMA_SP(A1)
           return
        endif
        XAM=-X+A*DLOG(X)
        if (XAM.GT.700.0.OR.A.GT.170.0) then
           write(*,*)'a and/or x too large'
           stop
        endif
        if (X.EQ.0.0) then
           GIN=0.0
           GIM = GAMMA_SP(A1)
        else if (X.LE.1.0+A) then
           S=1.0D0/A
           R=S
           do 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              if (DABS(R/S).LT.1.0D-15) go to 15
10         continue
15         GIN=DEXP(XAM)*S
           ga = GAMMA_SP(A1)
           GIM=GA-GIN
        else if (X.GT.1.0+A) then
           T0=0.0D0
           do 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         continue
           GIM=DEXP(XAM)/(X+T0)
!           GA = GAMMA_SP(A1)
!           GIN=GA-GIM
        endif
        
        gamxinf = GIM
        return
        end function GAMXINF

! #####################################################################

        double precision function GAMXINFDP(A1,X1)

!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                from x to infinity
!       Input :  a   --- Parameter ( a < 170 )
!                x   --- Argument 
!       Output:  GIM --- Gamma(a,x) t=x,Infinity
!       Routine called: GAMMA for computing gamma_dp(x)
!       ===================================================

!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        real :: a1,x1
! dont declare gamma_dp because it is within the module
!        double precision :: gamma_dp
        double precision :: xam,dlog,s,r,ga,t0,a,x
        integer :: k
        double precision :: gin, gim
        
        a = a1
        x = x1
        if ( x1 <= 0.0 ) then
           gamxinfdp = GAMMA_DP(A)
           return
        endif
        XAM=-X+A*DLOG(X)
        if (XAM.GT.700.0.OR.A.GT.170.0) then
           write(*,*)'a and/or x too large'
           stop
        endif
        if (X.EQ.0.0) then
           GIN=0.0
           GIM = GAMMA_dp(A)
        else if (X.LE.1.0+A) then
           S=1.0D0/A
           R=S
           do 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              if (DABS(R/S).LT.1.0D-15) go to 15
10         continue
15         GIN=DEXP(XAM)*S
           ga = GAMMA_DP(A)
           GIM=GA-GIN
        else if (X.GT.1.0+A) then
           T0=0.0D0
           do 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         continue
           GIM=DEXP(XAM)/(X+T0)
!           GA = GAMMA_dp(A)
!           GIN=GA-GIM
        endif
        
        gamxinfdp = GIM
        return
        end function GAMXINFDP


! #####################################################################

      real function gaminterp(ratio, alp, luindex, ilh)
      
      implicit none

      real, intent(in) :: ratio, alp
      integer, intent(in) :: ilh  ! 1 = graupel, 2 = hail
      integer, intent(in) :: luindex ! which argument: 
                         ! gamxinflu(i,j,1,1) = x/y
                          ! gamxinflu(i,j,2,1) = gamxinf( 2.0+alp, ratio )/y
                          ! gamxinflu(i,j,3,1) = gamxinf( 2.5+alp+0.5*bxh, ratio )/y
                          ! gamxinflu(i,j,5,1) = gamxinf( 5.0+alp, ratio )/y
                          ! gamxinflu(i,j,6,1) = gamxinf( 5.5+alp+0.5*bxh, ratio )/y

      
      real :: delx, dely, tmp1, tmp2, temp3
      integer :: i,j,ip1,jp1 !,ilh
      
!      ilh = Abs(ilh0)


           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           delx = Min(maxratiolu,ratio) - float(i)*dqiacrratio
           dely = alp - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,luindex,ilh) + delx*dqiacrratioinv*         &
     &                 (gamxinflu(ip1,j,luindex,ilh) - gamxinflu(i,j,luindex,ilh))
           tmp2 = gamxinflu(i,jp1,luindex,ilh) + delx*dqiacrratioinv*       &
     &                 (gamxinflu(ip1,jp1,luindex,ilh) - gamxinflu(i,jp1,luindex,ilh))
           
           ! interpolate along alpha; 
           
           gaminterp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
           ! debug
!           IF ( ilh0 < 0 ) THEN
!             write(0,*) 'gaminterp: ',i,j,ilh,ratio,delx,dely,gamxinflu(i,j,luindex,ilh),tmp1,tmp2
!           ENDIF
           
        end function gaminterp
! #####################################################################

!**************************** GAML02 *********************** 
!  This calculates Gamma(0.2,x)/Gamma[0.2], where is a ratio
!   It is used for qiacr with the gamma of volume to calculate what 
!   fraction of drops exceed a certain size (this version is for 40 micron drops)
! **********************************************************
      real function GAML02(x) 
      implicit none
      integer ig, i, ii, n, np
      real x
      integer ng
      parameter(ng=12)
      real gamxg(ng), xg(ng)
      data xg/0.01,0.02,0.025,0.04,0.075,0.1,0.25,0.5,0.75,1.,2.,10./ 
      data gamxg/  &
     &  7.391019203578037e-8,0.02212726874591478,0.06959352407989682, &
     &  0.2355654024970809,0.46135930387500346,0.545435791452399,     &
     &  0.7371571313308203,                                           &
     &  0.8265676632204345,0.8640182781845841,0.8855756211304151,     &
     &  0.9245079225301251,                                           &
     &  0.9712578342732681/
      if ( x .ge. xg(ng) ) then
        gaml02 = xg(ng)
        return
      endif
      if ( x .lt. xg(1) ) then
        gaml02 = 0.0
        return
      endif
      do ii = 1,ng-1
        i = ng - ii
        n = i
        np = n + 1
        if ( x .ge. xg(i) ) then
!         GOTO 2 
          gaml02 = gamxg(N)+((X-XG(N))/(XG(NP)-XG(N)))* &
     &            ( gamxg(NP) - gamxg(N) ) 
          return
        endif
      enddo
      return
      end function GAML02

!**************************** GAML02d300 *********************** 
!  This calculates Gamma(0.2,x)/Gamma[0.2], where is a ratio
!   It is used for qiacr with the gamma of volume to calculate what 
!   fraction of drops exceed a certain size (this version is for 300 micron drops) (see zieglerstuff.nb)
! **********************************************************
     real function GAML02d300(x)
      implicit none
      integer ig, i, ii, n, np
      real x
      integer ng
      parameter(ng=9)
      real gamxg(ng), xg(ng)
      data xg/0.04,0.075,0.1,0.25,0.5,0.75,1.,2.,10./ 
      data gamxg/                           &
     &  0.0,                                  &
     &  7.391019203578011e-8,0.0002260640810600053,  &
     &  0.16567071824457152,                         &
     &  0.4231369044918005,0.5454357914523988,       &
     &  0.6170290936864555,                           &
     &  0.7471346054110058,0.9037156157718299 /
      if ( x .ge. xg(ng) ) then
        GAML02d300 = xg(ng)
        return
      endif
      if ( x .lt. xg(1) ) then
        GAML02d300 = 0.0
        return
      endif
      do ii = 1,ng-1
        i = ng - ii
        n = i
        np = n + 1
        if ( x .ge. xg(i) ) then
!         GOTO 2 
          GAML02d300 = gamxg(N)+((X-XG(N))/(XG(NP)-XG(N)))*  &
     &            ( gamxg(NP) - gamxg(N) ) 
          return
        endif
      enddo
      return
      end function GAML02d300
!c

! #####################################################################
! #####################################################################

!**************************** GAML02 *********************** 
!  This calculates Gamma(0.2,x)/Gamma[0.2], where is a ratio
!   It is used for qiacr with the gamma of volume to calculate what 
!   fraction of drops exceed a certain size (this version is for 500 micron drops) (see zieglerstuff.nb)
! **********************************************************
      real function GAML02d500(x) 
      implicit none
      integer ig, i, ii, n, np
      real x
      integer ng
      parameter(ng=9)
      real gamxg(ng), xg(ng)
      data xg/0.04,0.075,0.1,0.25,0.5,0.75,1.,2.,10./ 
      data gamxg/  &
     &  0.0,0.0,   &
     &  2.2346039e-13, 0.0221272687459,  &
     &  0.23556540,  0.38710348,         &
     &  0.48136183,0.6565833,            &
     &  0.86918315 /
      if ( x .ge. xg(ng) ) then
        GAML02d500 = xg(ng)
        return
      endif
      if ( x .lt. xg(1) ) then
        GAML02d500 = 0.0
        return
      endif
      do ii = 1,ng-1
        i = ng - ii
        n = i
        np = n + 1
        if ( x .ge. xg(i) ) then
!         GOTO 2 
          GAML02d500 = gamxg(N)+((X-XG(N))/(XG(NP)-XG(N)))*  &
     &            ( gamxg(NP) - gamxg(N) ) 
          return
        endif
      enddo
      return
      end function GAML02d500
!c

! #####################################################################

! #####################################################################


        real function BETA(P,Q)
!
!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA for computing gamma(x)
!       ==========================================
!
!        IMPLICIT real (A-H,O-Z)
        implicit none
        double precision p1,gp,q1,gq, ppq,gpq
        real p,q
        
        p1 = p
        q1 = q
        call GAMMADP(P1,GP)
        call GAMMADP(Q1,GQ)
        PPQ=P1+Q1
        call GAMMADP(PPQ,GPQ)
        beta=GP*GQ/GPQ
        return
        end function BETA

! #####################################################################
! #####################################################################

      double precision function GAMMA_DP(xx)

      implicit none
      double precision xx
      integer j

! Double precision ser,stp,tmp,x,y,cof(6)

      real*8 ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d+0,  &
     &            -86.50532032941677d0,   &
     &             24.01409824083091d0,   &
     &             -1.231739572450155d0,  &
     &              0.1208650973866179d-2,&
     &             -0.5395239384953d-5,   &
     &              2.5066282746310005d0/

      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x + 0.5d0)*Log(tmp) - tmp
      ser = 1.000000000190015d0
      do j=1,6
        y = y + 1.0d0
        ser = ser + cof(j)/y
      end do
      gamma_dp = Exp(tmp + log(stp*ser/x))

      return
      end function gamma_dp
! #####################################################################

        subroutine GAMMADP(X,GA)
!
!       ==================================================
!       Purpose: Compute gamma function Gamma(x)
!       Input :  x  --- Argument of Gamma(x)
!                       ( x is not equal to 0,-1,-2,...)
!       Output:  GA --- gamma(x)
!       ==================================================
!
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit none
        
        double precision, parameter :: PI=3.141592653589793D0
        double precision :: x,ga,z,r,gr
        integer :: k,m1,m
        
        double precision :: G(26)
        
        if (X.EQ.INT(X)) then
           if (X.GT.0.0D0) then
              GA=1.0D0
              M1=X-1
              do K=2,M1
                GA=GA*K
              enddo
           else
              GA=1.0D+300
           endif
        else
           if (DABS(X).GT.1.0D0) then
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              do K=1,M
                 R=R*(Z-K)
              enddo
              Z=Z-M
           else
              Z=X
           endif
           data G/1.0D0,0.5772156649015329D0,                  &
     &          -0.6558780715202538D0, -0.420026350340952D-1,  &
     &          0.1665386113822915D0,-.421977345555443D-1,     &
     &          -.96219715278770D-2, .72189432466630D-2,       &
     &          -.11651675918591D-2, -.2152416741149D-3,       &
     &          .1280502823882D-3, -.201348547807D-4,          &
     &          -.12504934821D-5, .11330272320D-5,             &
     &          -.2056338417D-6, .61160950D-8,                 &
     &          .50020075D-8, -.11812746D-8,                   &
     &          .1043427D-9, .77823D-11,                       &
     &          -.36968D-11, .51D-12,                          &
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           do K=25,1,-1
             GR=GR*Z+G(K)
           enddo
           GA=1.0D0/(GR*Z)
           if (DABS(X).GT.1.0D0) then
              GA=GA*R
              if (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           endif
        endif
        return
        end subroutine GAMMADP


! #####################################################################
! #####################################################################
!
!
! #####################################################################
      function delbk(bb,nu,mu,k)
!   
!  Purpose: Caluculates collection coefficients following Siefert (2006)
!
!  delbk is equation (90) (b collecting b -- self-collection)
!  mass-diameter relationship: D = a*x**(b), where x = particle mass
!  general distribution: n(x) = A*x**(nu)*Exp(-lam*x**(mu))
!  where
!      A = mu*N/(Gamma((nu+1)/mu)) *lam**((nu+1)/mu)
!
!      lam = ( Gamma((nu+1)/mu)/Gamma((nu+2)/mu) * xbar )**(-mu)
!
!     where  xbar = L/N  (mass content)/(number concentration) = q*rhoa/N
!

      implicit none
      real delbk
      real nu, mu, bb
      integer k
      
      real tmp, del
      real x1, x2, x3, x4
      integer i

        tmp = ((1.0 + nu)/mu)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((2.0 + nu)/mu)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1.0 + 2.0*bb + k + nu)/mu)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x3 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
      
!      delbk =  &
!     &  ((Gamma_sp((1.0 + nu)/mu)/Gamma_sp((2.0 + nu)/mu))**(2.0*bb + k)* &
!     &    Gamma_sp((1.0 + 2.0*bb + k + nu)/mu))/Gamma_sp((1.0 + nu)/mu)

      delbk =  &
     &  ((x1/x2)**(2.0*bb + k)* &
     &    x3)/x1
      
      return
      end  function delbk
      
! #####################################################################
!
!
! #####################################################################
! Equation (91) in Seifert and Beheng (2006) ("a" collecting "b")
      function delabk(ba,bb,nua,nub,mua,mub,k)
      
      implicit none
      real delabk
      real nua, mua, ba
      integer k
      real nub, mub, bb
      
      integer i
      real tmp,del
      
      real g1pnua, g2pnua, g1pbapnua, g1pbbpk, g1pnub, g2pnub
      
        tmp = (1. + nua)/mua
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        if ( i+1 > ngm0 ) then
          write(0,*) 'delabk: i+1 > ngm0!!!!',i,ngm0,nua,mua,tmp
           stop
        endif
        g1pnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!        write(91,*) 'delabk: g1pnua,gamma = ',g1pnua,Gamma_sp((1. + nua)/mua)

        tmp = ((2. + nua)/mua)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g2pnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1. + ba + nua)/mua)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1pbapnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1. + nub)/mub)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1pnub = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((2 + nub)/mub)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g2pnub = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = ((1. + bb + k + nub)/mub)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1pbbpk = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

      delabk =  &
     &  (2.*(g1pnua/g2pnua)**ba*     &
     &    g1pbapnua*                                               &
     &    (g1pnub/g2pnub)**(bb + k)*                                &
     &    g1pbbpk)/                                                &
     &  (g1pnua*g1pnub)              
      
      return
      end function delabk


! #####################################################################
!
!
! #####################################################################

     real function cnudiag(ccw)
   ! function to diagnose cloud drop DSD shape parameter (cnu) based on Chandrakar et al. 2016 (PNAS) data
     
     implicit none
     
     real :: ccw
     
     ! local
     real :: vard, tmp
     
     real :: vard_nu(9,2)  ! table of variance vs. shape parameter
     real :: vard_ccw(4,3) ! table of Chandrakar et al. 2016

      data vard_nu(:,1) / 0.406873, 0.381485, 0.360077, 0.34174, 0.325827, 0.311862, 0.299489, 0.288433, 0.278483 /
      data vard_nu(:,2) / -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5 /
     
      data vard_ccw(:,1) / 0.40, 0.38, 0.35, 0.28 /
      data vard_ccw(:,2) / 21.3e6, 76.9e6, 201.2e6, 564.6e6 /
      data vard_ccw(:,3) / -0.274661, -0.193569, -0.0468808, 0.484082 /
      
! ----------------

      if ( ccw <= vard_ccw(1,2) ) then
        vard = 0.4
        cnudiag = -0.274661
      elseif ( ccw >= vard_ccw(4,2) ) then
        vard = 0.28
        cnudiag = 0.484082
      else
        ! interpolate
        if ( ccw <= vard_ccw(2,2) ) then
         ! range of 21.3e6 to 76.9e6
          
          cnudiag = vard_ccw(1,3) + (vard_ccw(2,3)-vard_ccw(1,3))*(ccw - vard_ccw(1,2))/( vard_ccw(2,2) - vard_ccw(1,2) )
          
        elseif ( ccw <= vard_ccw(3,2) ) then
         ! range of 76.9e6 to 201.2e6

          cnudiag = vard_ccw(2,3) + (vard_ccw(3,3)-vard_ccw(2,3))*(ccw - vard_ccw(2,2))/( vard_ccw(3,2) - vard_ccw(2,2) )
        
        else
         ! range of 201.2e6 to 564.6e6

          cnudiag = vard_ccw(3,3) + (vard_ccw(4,3)-vard_ccw(3,3))*(ccw - vard_ccw(3,2))/( vard_ccw(4,2) - vard_ccw(3,2) )
        
        endif
      
      
      endif
     return
     end function cnudiag

! #####################################################################
!
! #####################################################################
!--------------------------------------------------------------------------
      subroutine cld_cpu(string)

      implicit none
      character( LEN = * ) string
      
      return
      
      end subroutine cld_cpu

!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
! #######################################################################
!  HAILMAXD - calculated maximum expected hail size
! #######################################################################
     subroutine hailmaxd(dtp,nx,ny,nz,an,na,nor,norz,alpha2d,dn,  &
     &                    hailmax1d,hailmaxk1,jslab ) 
!
! Calculate maximum hail size from the tail of of the distribution. The value
! of thresh_conc sets the minimum concentration in the integral over (Dmax, Inf).
! This uses the lookup tables for incomplete gamma functions and simply search for
! the expected value (and linearly interpolate) on D.
!
!  Written by ERM 7/2023
!
!
!
      implicit none

      integer nx,ny,nz,nor,norz,ngt,jgs,na,ia
      integer id ! =1 use density, =0 no density
!      integer :: its,ite ! x-range to calculate
      
      integer ng1
      parameter(ng1 = 1)

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

!      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
      real dtp
      real alpha2d(-nor+1:nx+nor,1,-norz+1:nz+norz,3)  ! array for PSD shape parameters
      real  :: hailmax1d(nx,ny),hailmaxk1(nx,ny)
      integer infdo
      integer jslab ! which line of xfall to use
            
      integer ix,jy,kz,ndfall,n,k,il,in
      double precision :: tmp, ratio, del, g1palp
      real, parameter :: dz = 200.

      real :: db1(nx,nz+1),dtz1(nz+1,nx,0:1),dz2dinv(nz+1,nx),db1inv(nx,nz+1)
      
      real :: rhovtzx(nz,nx)

      real :: alp, diam, diam1, hwdn
      
!      real, parameter :: cmin = 0.001 ! threshold number per m^3 for maximum diamter (threshold from diag_nwp)
      double precision, parameter:: thresh_conc = 0.0005d0                 ! number conc. of graupel/hail per cubic meter
      real :: cwchtmp,cwchltmp, maxdia

!-----------------------------------------------------------------------------

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: plo, phi
      integer :: ialp, i, j

      logical :: debug_mpi = .TRUE.

! ###################################################################


      if ( lh > 1 ) then
        cwchtmp  = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
      endif
      if ( lhl > 1 ) then
        cwchltmp = ((3. + dnu(lhl))*(2. + dnu(lhl))*(1.0 + dnu(lhl)))**(-1./3.)
      endif


      kzb = 1
      kze = nz

      ixb = 1  ! aliased its
      ixe = nx ! aliased ite


      jy = jslab
      jgs = jy


!      hailmax1d(:,jy) = 0.0
!      hailmaxk1(:,jy) = 0.0

      if ( ndebug .gt. 0 ) write(0,*) 'dbg = 3a'


! first graupel, even if hail is also predicted, since graupel can sometime be large on its own
      if ( lh > 1 .and. lnh > 1 ) then
      do kz = kzb,kze
      do ix = ixb,ixe
        if ( an(ix,jy,kz,lh) .gt. qxmin(lh) .and. an(ix,jy,kz,lnh) .gt. thresh_conc ) then
          if ( lvh .gt. 1 ) then
            hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
          else
            hwdn = rho_qh
          endif

          tmp = 1. + alpha2d(ix,1,kz,2)
          i = Int(dgami*(tmp))
          del = tmp - dgam*i
          g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

          tmp = dn(ix,jy,kz)*an(ix,jy,kz,lh)/(hwdn*an(ix,jy,kz,lnh))
          diam = (6.0*tmp/pi)**(1./3.)
          if ( lzh > 1 ) then ! 3moment
            cwchtmp = ((3. + alpha2d(ix,1,kz,2))*(2. + alpha2d(ix,1,kz,2))*(1.0 + alpha2d(ix,1,kz,2)))**(-1./3.)
          endif
          diam1 = diam*cwchtmp ! characteristic diameter, i.e., 1/lambda
         ! want cxd1 = thresh_conc
         !  tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
         ! cxd1 = cx(mgs,lh)*(tmp)/g1palp
         ! tmp = thresh_conc*g1palp/cx
         ! 
         tmp = thresh_conc*g1palp/an(ix,jy,kz,lnh)
         alp = alpha2d(ix,1,kz,2)
         ! gamxinflu(i,j,luindex,ilh)
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           ratio = 0.0
           maxdia = 0.0
           ! eventually could replace with bisection search, but final value of i is usually small
           ! compared to nqiacrratio
           do i = 0,nqiacrratio-1
              if ( gamxinflu(i,j,1,1) >= tmp .and. tmp >= gamxinflu(i+1,j,1,1) ) then
               !  interpolate here for FWIW
                ratio = i*dqiacrratio
                del = tmp - gamxinflu(i,j,1,1)
                ratio = (float(i) + del/(gamxinflu(i+1,j,1,1) - gamxinflu(i,j,1,1)))*dqiacrratio
                exit
              endif
           enddo
           
           if ( ratio > 0.0 ) then
              maxdia = ratio*diam1 ! units of m
           endif

           if ( kz == kzb ) then
             hailmaxk1(ix,jy) = Max( maxdia, hailmaxk1(ix,jy) )
!             IF ( maxdia > 0.1 ) THEN
!            IF ( an(ix,jy,kz,lh) > 1.e-4 ) THEN
!              write(0,*) 'maxdia,tmp,alp,ratio,diam,diam1= ',maxdia,tmp,alp,ratio,diam*100.,diam1*100.
!              write(0,*) 'hwdn, cxhl, qx, g1palp = ',hwdn, an(ix,jy,kz,lnhl), an(ix,jy,kz,lhl), g1palp
!              write(0,*) 'j,gamxinflu(0,2,4) = ',j,gamxinflu(0,j,1,1),gamxinflu(2,j,1,1), &
!                gamxinflu(4,j,1,1)
!            ENDIF
           endif
           
           hailmax1d(ix,jy) = Max(maxdia, hailmax1d(ix,jy) )

        ! 

        endif

      enddo
      enddo

      endif ! lh

! And diam for hail if present
      if ( lhl > 1 .and. lnhl > 1 ) then
      do kz = kzb,kze
      do ix = ixb,ixe
        if ( an(ix,jy,kz,lhl) .gt. qxmin(lhl) .and. an(ix,jy,kz,lnhl) .gt. thresh_conc ) then
          if ( lvhl .gt. 1 ) then
            hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
          else
            hwdn = rho_qhl
          endif

          tmp = 1. + alpha2d(ix,1,kz,3)
          i = Int(dgami*(tmp))
          del = tmp - dgam*i
          g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

          tmp = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/(hwdn*an(ix,jy,kz,lnhl))
          diam = (6.0*tmp/pi)**(1./3.)
          if ( lzhl > 1 ) then ! 3moment
            cwchltmp = ((3. + alpha2d(ix,1,kz,3))*(2. + alpha2d(ix,1,kz,3))*(1.0 + alpha2d(ix,1,kz,3)))**(-1./3.)
          endif
          diam1 = diam*cwchltmp ! characteristic diameter, i.e., 1/lambda
         ! want cxd1 = thresh_conc
         !  tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
         ! cxd1 = cx(mgs,lh)*(tmp)/g1palp
         ! tmp = thresh_conc*g1palp/cx
         ! 
         tmp = thresh_conc*g1palp/an(ix,jy,kz,lnhl)
         alp = alpha2d(ix,1,kz,3)
         ! gamxinflu(i,j,luindex,ilh)
           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           ratio = 0.0
           maxdia = 0.0
           ! eventually could replace with bisection search, but final value of i is usually small
           ! compared to nqiacrratio
           do i = 0,nqiacrratio-1
              if ( gamxinflu(i,j,1,1) >= tmp .and. tmp >= gamxinflu(i+1,j,1,1) ) then
               !  interpolate here for FWIW
                ratio = i*dqiacrratio
                del = tmp - gamxinflu(i,j,1,1)
                ratio = (float(i) + del/(gamxinflu(i+1,j,1,1) - gamxinflu(i,j,1,1)))*dqiacrratio
                exit
              endif
           enddo
           
           if ( ratio > 0.0 ) then
              maxdia = ratio*diam1 ! units of m
           endif

           if ( kz == kzb ) then
             hailmaxk1(ix,jy) = Max( maxdia, hailmaxk1(ix,jy) )
!             IF ( maxdia > 0.1 ) THEN
!            IF ( an(ix,jy,kz,lhl) > 1.e-4 ) THEN
!              write(0,*) 'maxdia,tmp,alp,ratio,diam,diam1= ',maxdia,tmp,alp,ratio,diam*100.,diam1*100.
!              write(0,*) 'hwdn, cxhl, qx, g1palp = ',hwdn, an(ix,jy,kz,lnhl), an(ix,jy,kz,lhl), g1palp
!              write(0,*) 'j,gamxinflu(0,2,4) = ',j,gamxinflu(0,j,1,1),gamxinflu(2,j,1,1), &
!                gamxinflu(4,j,1,1)
!            ENDIF
           endif
           
           hailmax1d(ix,jy) = Max(maxdia, hailmax1d(ix,jy) )

        ! 

        endif

      enddo
      enddo

      endif


     end subroutine HAILMAXD
! #######################################################################
! #######################################################################
     subroutine sediment1d(dtp,nx,ny,nz,an,na,nor,norz,xfall,dn,dz3d,dz3dinv, &
     &                    t0,t7,infdo,jslab,its,jts,  &
     &   timesed1,timesed2,timesed3,zmaxsed,timesetvt) ! used for timing
!
! Sedimentation driver -- column by column
!
!  Written by ERM 10/2011
!
!
!
      implicit none

      integer nx,ny,nz,nor,norz,ngt,jgs,na,ia
      integer id ! =1 use density, =0 no density
      integer :: its,jts ! SW point of local tile
      
      integer ng1
      parameter(ng1 = 1)

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real dz3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real dz3dinv(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

!      real gz(-nor+ng1:nz+nor),z1d(-nor+ng1:nz+nor,4)
      real dtp
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      integer infdo
      integer jslab ! which line of xfall to use
            
      integer ix,jy,kz,ndfall,n,k,il,in
      real tmp, vtmax, dtptmp, dtfrac
      real, parameter :: dz = 200.


      real, allocatable :: db1(:,:), dtz1(:,:,:),dz2dinv(:,:),db1inv(:,:) ! db1(nx,nz+1),dtz1(nz+1,nx,0:1),dz2dinv(nz+1,nx),db1inv(nx,nz+1)
      real, allocatable :: rhovtzx(:,:)
      real, allocatable :: xfall0(:,:), xvt(:,:,:,:),tmpn(:,:,:),tmpn2(:,:,:),z(:,:,:)
      
      double precision :: timesed1,timesed2,timesed3, zmaxsed,timesetvt,dummy
      double precision :: dt1,dt2,dt3,dt4

      integer :: ngs ! = 512
      integer :: ngscnt,mgs,ipconc0

      real, allocatable ::  qx(:,:)
      real, allocatable ::  qxw(:,:)
      real, allocatable ::  cx(:,:)
      real, allocatable ::  xv(:,:)
      real, allocatable ::  vtxbar(:,:,:)
      real, allocatable ::  xmas(:,:)
      real, allocatable ::  xdn(:,:)
      real, allocatable ::  xdia(:,:,:)
      real, allocatable ::  vx(:,:)
      real, allocatable ::  alpha(:,:)
      real, allocatable ::  zx(:,:)
      logical, allocatable :: hasmass(:,:)

      integer, allocatable :: igs(:),kgs(:)
      
      real, allocatable :: rho0(:),temcg(:)

      real, allocatable :: temg(:)
      
      real, allocatable :: rhovt(:)
      
      real, allocatable :: cwnc(:),cinc(:)
      real, allocatable :: fadvisc(:),cwdia(:),cipmas(:)
      
      real, allocatable :: cnina(:),cimas(:)
      
      real, allocatable :: cnostmp(:)
      
      real :: cimasn,cimasx
      

!-----------------------------------------------------------------------------

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: plo, phi

      logical :: debug_mpi = .TRUE.

! ###################################################################


      allocate( db1(nx,nz+1),dtz1(nz+1,nx,0:1),dz2dinv(nz+1,nx),db1inv(nx,nz+1),rhovtzx(nz,nx) )
      allocate( xfall0(nx,ny), xvt(nz+1,nx,3,lc:lhab), tmpn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) )
      allocate( tmpn2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz), z(-nor+ng1:nx+nor,-norz+ng1:nz+norz,lr:lhab))

      ngs = nz+3
      
      allocate( qx(ngs,lv:lhab),  &
                qxw(ngs,ls:lhab),  &
                cx(ngs,lc:lhab),  &
                xv(ngs,lc:lhab),  &
                vtxbar(ngs,lc:lhab,3),  &
                xmas(ngs,lc:lhab),  &
                xdn(ngs,lc:lhab),  &
                xdia(ngs,lc:lhab,3),  &
                vx(ngs,li:lhab),  &
                alpha(ngs,lc:lhab),  &
                zx(ngs,lr:lhab),     &
                hasmass(nx,lc+1:lhab), &
                igs(ngs),kgs(ngs), &
                rho0(ngs),temcg(ngs),temg(ngs), rhovt(ngs), &
                cwnc(ngs),cinc(ngs), &
                fadvisc(ngs),cwdia(ngs),cipmas(ngs), &
                cnina(ngs),cimas(ngs), &
                cnostmp(ngs) )

      kzb = 1
      kze = nz

      ixb = 1
      ixe = nx


      jy = 1
      jgs = jy


!
!  zero the precip flux arrays (2d)
!

      xvt(:,:,:,:) = 0.0

      if ( ndebug .gt. 0 ) write(0,*) 'dbg = 3a'


      do kz = kzb,kze
      do ix = ixb,ixe
       db1(ix,kz) = dn(ix,jy,kz)
       db1inv(ix,kz) = 1./dn(ix,jy,kz)
       rhovtzx(kz,ix) = Sqrt(rho00*Min(1.0/0.05, db1inv(ix,kz))) ! prevent excessive rhovt
      enddo
      enddo

      do kz = kzb,kze
      do ix = ixb,ixe
       dtz1(kz,ix,0) = dz3dinv(ix,jy,kz)
       dtz1(kz,ix,1) = dz3dinv(ix,jy,kz)*db1inv(ix,kz) 
       dz2dinv(kz,ix) = dz3dinv(ix,jy,kz)
      enddo
      enddo

      if ( lzh .gt. 1 ) then
      do kz = kzb,kze
      do ix = ixb,ixe
        an(ix,jy,kz,lzh) = Max( 0., an(ix,jy,kz,lzh) )
      enddo
      enddo
      endif

      
      do il = lc+1,lhab
       do ix = ixb,ixe
!        hasmass(ix,il) = Any( an(ix,jy,:,il) > qxmin(il) )
       enddo
      enddo




      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3a2'

! loop over columns
      do ix = ixb,ixe
      
      dummy = 0.d0

      
      call ziegfall1d(nx,ny,nz,nor,norz,na,dtp,jgs,ix, & 
     &  xvt, rhovtzx, & 
     &  an,dn,ipconc,t0,t7,cwmasn,cwmasx, & 
     &  cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,xvmn,xvmx, & 
     &  ngs,qx,qxw,cx,xv,vtxbar,xmas,xdn,xdia,vx,alpha,zx,igs,kgs, &
     &  rho0,temcg,temg,rhovt,cwnc,cinc,fadvisc,cwdia,cipmas,cnina,cimas, &
     &  cnostmp,              &
     &  infdo,0               &
     & )


! loop over each species and do sedimentation for all moments
     do il = lc,lhab
       if ( ido(il) == 0 ) cycle

!       IF ( .not. hasmass(ix,il) ) CYCLE

!      plo = nz
!      phi = 0


      vtmax = 0.0
      
      do kz = kzb,kze

      ! apply limit vtmaxsed (08/20/2015)
      xvt(kz,ix,1,il) = Min( vtmaxsed,  xvt(kz,ix,1,il) )
      xvt(kz,ix,2,il) = Min( vtmaxsed,  xvt(kz,ix,2,il) )
      xvt(kz,ix,3,il) = Min( vtmaxsed,  xvt(kz,ix,3,il) )
      
      vtmax = Max(vtmax,xvt(kz,ix,1,il)*dz2dinv(kz,ix))
      vtmax = Max(vtmax,xvt(kz,ix,2,il)*dz2dinv(kz,ix))
      vtmax = Max(vtmax,xvt(kz,ix,3,il)*dz2dinv(kz,ix))

!      IF ( dtp*xvt(kz,ix,1,il)*dz2dinv(kz,ix) >= 0.7 .or. &
!     &     dtp*xvt(kz,ix,2,il)*dz2dinv(kz,ix) >= 0.7 .or. &
!     &     dtp*xvt(kz,ix,3,il)*dz2dinv(kz,ix) >= 0.7 ) THEN
!          
!          zmaxsed = Max(zmaxsed, float(kz) )
!!          plo = Min(plo,kz)
!!          phi = Max(phi,kz)
!           
!      ENDIF
      
      enddo
      
      if ( vtmax == 0.0 ) cycle


      
      if ( dtp*vtmax .lt. 0.7 ) then ! check whether multiple steps are needed.
        ndfall = 1
      else
       if ( dtp > 20.0 ) then ! more stringent subdivision for large time steps
         ndfall = Max(2, Int(dtp*vtmax/0.7) + 1)
       else ! more relaxed for small time steps, but might still be a problem for very thin vertical layers near the ground
         ndfall = 1+Int(dtp*vtmax + 0.301)
       endif
      endif
      
      if ( ndfall .gt. 1 ) then
        dtptmp = dtp/real(ndfall)
!        write(0,*) 'subdivide fallout on its,jts,ix,plo,phi = ',its,jts,ix,plo,phi
!        write(0,*) 'for il,jsblab,c,ndfall = ',il,jslab,dtp*vtmax,ndfall
      else
        dtptmp = dtp
      endif
      
      dtfrac = dtptmp/dtp


      do n = 1,ndfall

      if ( do_accurate_sedimentation .and. n .ge. 2 .and. ( n == interval_sedi_vt*(n/interval_sedi_vt) ) ) then
!
!  zero the precip flux arrays (2d)
!
      
      dummy = 0.d0

      xvt(kzb:kze,ix,1:3,il) = 0.0 ! reset to zero because routine will only compute points with q > qmin

      call ziegfall1d(nx,ny,nz,nor,norz,na,dtp,jgs,ix, & 
     &  xvt, rhovtzx, & 
     &  an,dn,ipconc,t0,t7,cwmasn,cwmasx, & 
     &  cwradn, & 
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,xvmn,xvmx, & 
     &  ngs,qx,qxw,cx,xv,vtxbar,xmas,xdn,xdia,vx,alpha,zx,igs,kgs, &
     &  rho0,temcg,temg,rhovt,cwnc,cinc,fadvisc,cwdia,cipmas,cnina,cimas, &
     &  cnostmp,             &
     &  infdo,il)


      do kz = kzb,kze
      ! apply limit vtmaxsed (08/20/2015)
        xvt(kz,ix,1,il) = Min( vtmaxsed,  xvt(kz,ix,1,il) )
        xvt(kz,ix,2,il) = Min( vtmaxsed,  xvt(kz,ix,2,il) )
        xvt(kz,ix,3,il) = Min( vtmaxsed,  xvt(kz,ix,3,il) )
      enddo




      endif ! (n .ge. 2)


        if ( il >= lr .and. ( infall .eq. 3 .or. infall .eq. 4 ) .and. ln(il) > 0 ) then
           if ( (il .eq. lr .and. irfall .eq. infall .and. lzr < 1) .or. & 
                 (il .ge. lh .and. lz(il) .lt. 1 ) .or. (il == ls .and. isfall == infall ) ) then
            call calczgr1d(nx,ny,nz,nor,na,an,ixe,kze, & 
     &         z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il), lvol(il), xdn0(il), ix )
           endif
        endif

      if (ndebug .gt. 0 ) write(0,*) 'dbg = 1b'

! mixing ratio

      call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &             an,db1,il,1,xfall,dtz1,ix)


      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3c'

! volume

      if ( ldovol .and. il >= li ) then
        if ( lvol(il) .gt. 1 ) then
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &              an,db1,lvol(il),0,xfall,dtz1,ix)
        endif
      endif

! reflectivity

      if ( ipconc .ge. 6 ) then
        if ( lz(il) .gt. 1 ) then
         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,3,il), & 
     &              an,db1,lz(il),0,xfall,dtz1,ix)
        endif
      endif

      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3d'

      
      if ( ipconc .gt. 0 ) then !{
        if ( ipconc .ge. ipc(il) ) then

      if ( ( infall .ge. 2 .or. (infall .eq. 0 .and. il .lt. lh) ) .and. lz(il) .lt. 1) then !{
!
! load number conc. into tmpn to do fallout by mass-weighted mean fall speed
!  to put a lower bound on number conc.
!

        if ( ( infall .eq. 3 .or. infall .eq. 4 ) .and. (  (il == ls .and. isfall .eq. infall ) &
     &        .or. il .eq. lh .or. il .eq. lhl .or.  il == lf .or. & 
     &      ( il .eq. lr .and. irfall .eq. infall) ) ) then

          ! set up for method I+II
          do kz = kzb,kze
!            DO ix = ixb,ixe
              tmpn2(ix,jy,kz) = z(ix,kz,il)
!            ENDDO
          enddo
          do kz = kzb,kze
!            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
!            ENDDO
          enddo
        
        else
          ! set up for method II only
          do kz = kzb,kze
!            DO ix = ixb,ixe
              tmpn(ix,jy,kz) = an(ix,jy,kz,ln(il))
!            ENDDO
          enddo

        endif

      endif !}


      if (ndebug .gt. 0 ) write(0,*) 'dbg = 3f'

       in = 2
       if ( infall .eq. 1 ) in = 1

         call fallout1d(nx,ny,nz,nor,na,dtptmp,dtfrac,jgs,xvt(1,1,in,il), & 
     &        an,db1,ln(il),0,xfall,dtz1,ix)


         if ( lz(il) .lt. 1 ) then ! if not 3-moment, run one of the correction schemes
         if ( (infall .ge. 2 .or. infall .eq. 3) .and. .not. (infall .eq. 0 .and. il .ge. lh) & 
     &       .and. ( il .eq. lr .or. (il .ge. li .and. il .le. lhab) )) then
!     :        .or. il .eq. lhl )) THEN
           
           xfall0(:,jgs) = 0.0

           if ( ( infall .eq. 3 .or. infall .eq. 4 ) .and.  & 
     &        ( il .ge. lh .or. (il .eq. lr .and. irfall .eq. infall)    &
                   .or. (il .eq. ls .and. isfall .eq. infall) ) ) then
             call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,3,il), & 
     &         tmpn2,db1,1,0,xfall0,dtz1,ix)
             call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &         tmpn,db1,1,0,xfall0,dtz1,ix)
           else
             call fallout1d(nx,ny,nz,nor,1,dtptmp,dtfrac,jgs,xvt(1,1,1,il), & 
     &         tmpn,db1,1,0,xfall0,dtz1,ix)
           endif

           if ( ( infall .eq. 3 .or. infall .eq. 4 ) .and. ( (il .eq. lr .and. irfall .eq. infall) & 
     &            .or. il .ge. lh .or. (il == ls .and. isfall .eq. infall ) ) ) then
! "Method I" - dbz correction

             call calcnfromz1d(nx,ny,nz,nor,na,an,tmpn2,ixe,kze, & 
     &       z,db1,jgs,ipconc, dnu(il), il, ln(il), qxmin(il), xvmn(il), xvmx(il),tmpn,  & 
     &       lvol(il), xdn0(il), infall, ix)

           elseif ( infall .eq. 5 .and. il .ge. lh .or. ( il == lr .and. irfall == 5 ) ) then

             do kz = kzb,kze
!              DO ix = ixb,ixe
               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), 0.5* ( an(ix,jgs,kz,ln(il)) + tmpn(ix,jy,kz) ))
              
!              ENDDO
             enddo           

           elseif ( .not. (il .eq. lr .and. irfall .eq. 0) .and. .not. (il .eq. ls .and. isfall .eq. 0) ) then
! "Method II" M-wgt N-fallout correction

             do kz = kzb,kze
!              DO ix = ixb,ixe

               an(ix,jgs,kz,ln(il)) = Max( an(ix,jgs,kz,ln(il)), tmpn(ix,jy,kz) )
              
!              ENDDO
             enddo
           endif 
           endif ! lz(il) .lt. 1
           

         endif
        endif


      endif !}


      enddo ! n=1,ndfall
      enddo ! il
      
      enddo ! ix


      deallocate( db1,dtz1,dz2dinv,db1inv,rhovtzx )
      deallocate( xfall0, xvt, tmpn )
      deallocate( tmpn2, z)

      deallocate( qx,  &
                qxw,  &
                cx,  &
                xv,  &
                vtxbar,  &
                xmas,  &
                xdn,  &
                xdia,  &
                vx,  &
                alpha,  &
                zx,     &
                hasmass, &
                igs,kgs, &
                rho0,temcg,temg, rhovt, &
                cwnc,cinc, &
                fadvisc,cwdia,cipmas, &
                cnina,cimas, &
                cnostmp )

      return
      end subroutine SEDIMENT1D


! #####################################################################

!
! #####################################################################


!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
      subroutine fallout1d(nx,ny,nz,nor,na,dtp,dtfrac,jgs,vt,   &
     &  a,db1,ia,id,xfall,dtz1,ixcol)
!
! First-order, upwind fallout scheme
!
!  Written by ERM 6/10/2011
!
!
!
      implicit none

      integer nx,ny,nz,nor,ngt,jgs,na,ia
      integer id ! =1 use density, =0 no density
      integer ng1
      parameter(ng1 = 1)
      integer :: ixcol

!      real dz3dinv(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
!      real a(nx,ny,nz,na)
      real a(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na) ! quantity to be 'advected'
      real vt(nz+1,nx)  ! terminal speed for a
      real dtp,dtfrac
      real cmax
      real xfall(nx,ny,na)  ! array for stuff landing on the ground
      real db1(nx,nz+1),dtz1(nz+1,nx,0:1)

! Local
           
      integer ix,jy,kz,n,k
      integer iv1,iv2
      real tmp
      integer imn,imx,kmn,kmx
      real qtmp1(nz+1)

!-----------------------------------------------------------------------------

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .TRUE.

! ###################################################################

      jy = 1

      iv1 = 0
      iv2 = 0

      imn = nx
      imx = 1
      kmn = nz
      kmx = 1

      cmax = 0.0

      kzb = 1
      kze = nz

      ixb = ixcol
      ixe = ixcol
      ix  = ixcol

      qtmp1(nz+1) = 0.0
      
      do kz = kzb,kze
!        DO ix = ixb,ixe
!         cmax = Max(cmax, vt(ix,kz)*dz3dinv(ix,jy,kz)) 
         
         if ( id == 1 ) then
         qtmp1(kz) = a(ix,jgs,kz,ia)*vt(kz,ix)*db1(ix,kz)
         else
         qtmp1(kz) = a(ix,jgs,kz,ia)*vt(kz,ix)
         endif
         
         if ( a(ix,jgs,kz,ia) .ne. 0.0 ) then
!           imn = Min(ix,imn)
!           imx = Max(ix,imx)
           kmn = Min(kz,kmn)
           kmx = Max(kz,kmx)
         endif
!        ENDDO
      enddo
      
      kmn = Max(1,kmn-1)
      
! first check if fallout is worth doing
!      IF ( cmax .eq. 0.0 .or. imn .gt. imx ) THEN
!        RETURN
!      ENDIF
      
      if ( kmn == 1 ) then
      
      kz = 1
!      do ix = imn,imx ! 1,nx-1
         xfall(ix,jy,ia) = xfall(ix,jy,ia) + a(ix,jgs,kz,ia)*vt(kz,ix)*dtfrac
!      enddo
      
      endif

      do kz = 1,nz
!      do ix = 1,nx
        a(ix,jgs,kz,ia) =  a(ix,jgs,kz,ia) + dtp*dtz1(kz,ix,id)*(qtmp1(kz+1) - qtmp1(kz) )
!      enddo
      enddo

      
      return
      end subroutine FALLOUT1D

! ##############################################################################
! ##############################################################################

      subroutine calczgr1d(nx,ny,nz,nor,na,a,ixe,kze,              &
     &    z,db,jgs,ipconc, alpha, l,ln, qmin, xvmn,xvmx, lvol, rho_qx, ixcol)


      implicit none

      integer nx,ny,nz,nor,na,ngt,jgs
      integer :: ixcol
      integer, parameter :: norz = 3
      real a(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)
      real z(-nor+1:nx+nor,-nor+1:nz+nor,lr:lhab)   ! reflectivity
      real db(nx,nz+1)  ! air density
!      real gt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,ngt)

      integer ixe,kze
      real    alpha
      real    qmin
      real    xvmn,xvmx
      integer ipconc
      integer l   ! index for q
      integer ln  ! index for N
      integer lvol ! index for volume
      real    rho_qx


      integer ix,jy,kz
      real vr,qr,nrx,rd,xv,g1,zx,chw,xdn,ynu
      
      
      jy = jgs
      ix = ixcol
      
      if ( l .eq. lh .or. l .eq. lhl .or. ( l .eq. lr .and. imurain == 1 )  &
           .or. ( l .eq. ls .and. imusnow == 1 ) ) then
      
      
      do kz = 1,kze
          
          
          
          if (  a(ix,jy,kz,l) .gt. qmin .and. a(ix,jy,kz,ln) .gt. 1.e-15 ) then
            
            if ( lvol .gt. 1 ) then
                if ( a(ix,jy,kz,lvol) .gt. 0.0 ) then
                  xdn = db(ix,kz)*a(ix,jy,kz,l)/a(ix,jy,kz,lvol)
                  xdn = Min( 900., Max( hdnmn, xdn ) )
                else
                  xdn = rho_qx
                endif
            else
                xdn = rho_qx
            endif

            if ( l == lr ) xdn = 1000.

            qr = a(ix,jy,kz,l)
            xv = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
            chw = a(ix,jy,kz,ln)

             if ( xv .lt. xvmn .or. xv .gt. xvmx ) then
              xv = Min( xvmx, Max( xvmn,xv ) )
              chw = db(ix,kz)*a(ix,jy,kz,l)/(xv*xdn)
             endif

             g1 = (6.0 + alpha)*(5.0 + alpha)*(4.0 + alpha)/  &
     &            ((3.0 + alpha)*(2.0 + alpha)*(1.0 + alpha))
             zx = g1*db(ix,kz)**2*(a(ix,jy,kz,l))*a(ix,jy,kz,l)/chw
!             z(ix,kz,l)  = 1.e18*zx*(6./(pi*1000.))**2
             z(ix,kz,l)  = zx*(6./(pi*1000.))**2


!          IF ( ny.eq.2 .and. kz .ge. 25 .and. kz .le. 29 .and. z(ix,kz,l) .gt. 0. ) THEN
!             write(*,*) 'calczgr: z,dbz,xdn = ',ix,kz,z(ix,kz,l),10*log10(z(ix,kz,l)),xdn
!          ENDIF
          
          else
           
            z(ix,kz,l) = 0.0
           
          endif
          
      enddo
      
      elseif ( (l == ls .and. imusnow == 3) .or. ( l .eq. lr .and. imurain == 3 ) ) then

      xdn = rho_qx ! 1000.
      if ( l == ls ) ynu = snu
      if ( l == lr ) ynu = rnu
      
      do kz = 1,kze

          if (  a(ix,jy,kz,l) .gt. qmin .and. a(ix,jy,kz,ln) .gt. 1.e-15 ) then

            vr = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
!            z(ix,kz,l) = 3.6e18*(ynu+2.0)*a(ix,jy,kz,ln)*vr**2/(ynu+1.0)
            z(ix,kz,l) = 3.6*(ynu+2.0)*a(ix,jy,kz,ln)*vr**2/(ynu+1.0)
!            qr = a(ix,jy,kz,lr)
!            nrx = a(ix,jy,kz,lnr)
          
          else
           
            z(ix,kz,l) = 0.0
           
          endif
      
          
      enddo
      
      endif
      
      return
      
      end subroutine calczgr1d

! ##############################################################################
! ##############################################################################
!
!  Subroutine to correct number concentration to prevent reflectivity growth by 
!  sedimentation in 2-moment ZXX scheme.
!  Calculation is in a slab (constant jgs)
!

      subroutine calcnfromz1d(nx,ny,nz,nor,na,a,t0,ixe,kze,    &
     &    z0,db,jgs,ipconc, alpha, l,ln, qmin, xvmn,xvmx,t1, &
     &    lvol, rho_qx, infall, ixcol)

      
      implicit none

      integer nx,ny,nz,nor,na,ngt,jgs,ixcol

      real a(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,na)  ! sedimented N and q
      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)    ! sedimented reflectivity
      real t1(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor)    ! sedimented N (by Vm)
!      real gt(-nor+1:nx+nor,-nor+1:ny+nor,-nor+1:nz+nor,ngt)
      real z0(-nor+1:nx+nor,-nor+1:nz+nor,lr:lhab)   ! initial reflectivity

      real db(nx,nz+1)  ! air density
      
      integer ixe,kze
      real    alpha
      real    qmin
      real    xvmn,xvmx
      integer ipconc
      integer l   ! index for q
      integer ln  ! index for N
      integer lvol ! index for volume
      real    rho_qx
      integer infall
      
      
      integer ix,jy,kz
      double precision vr,qr,nrx,rd,g1,zx,chw,z,znew,zt,zxt
      real xv,xdn
      integer :: ndbz, nmwgt, nnwgt, nwlessthanz
      
      ndbz = 0
      nmwgt = 0
      nnwgt = 0
      nwlessthanz = 0
      

      
      jy = jgs
      ix = ixcol
      
      if ( l .eq. lh .or. l .eq. lhl .or. ( l == lr .and. imurain == 1 ) ) then
      
             g1 = (6.0 + alpha)*(5.0 + alpha)*(4.0 + alpha)/  &
     &            ((3.0 + alpha)*(2.0 + alpha)*(1.0 + alpha))
      
      do kz = 1,kze

         
          if (   t0(ix,jy,kz) .gt. 0. ) then ! {
            
            if ( lvol .gt. 1 ) then
               if ( a(ix,jy,kz,lvol) .gt. 0.0 ) then
                 xdn = db(ix,kz)*a(ix,jy,kz,l)/a(ix,jy,kz,lvol)
                 xdn = Min( 900., Max( hdnmn, xdn ) )
               else 
                 xdn = rho_qx
               endif
            else
               xdn = rho_qx
            endif
            
            if ( l == lr ) xdn = 1000.
          
            qr = a(ix,jy,kz,l)
            xv = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
            chw = a(ix,jy,kz,ln)

             if ( xv .lt. xvmn .or. xv .gt. xvmx ) then
              xv = Min( xvmx, Max( xvmn,xv ) )
              chw = db(ix,kz)*a(ix,jy,kz,l)/(xv*xdn)
             endif

             zx = g1*db(ix,kz)**2*( a(ix,jy,kz,l))*a(ix,jy,kz,l)/chw
             z  = zx*(6./(pi*1000.))**2

            
           if ( (z .gt. t0(ix,jy,kz) .and. z .gt. 0.0 .and.  &
     &           t0(ix,jy,kz) .gt. z0(ix,kz,l) )) then !{
           
            zx = t0(ix,jy,kz)/((6./(pi*1000.))**2)
            
            nrx =  g1*db(ix,kz)**2*( a(ix,jy,kz,l))*a(ix,jy,kz,l)/zx
            if ( infall .eq. 3 ) then
              if ( nrx .gt. a(ix,jy,kz,ln) ) then
                ndbz = ndbz + 1
                if ( t1(ix,jy,kz) .lt. ndbz ) nwlessthanz = nwlessthanz + 1
              else
                nnwgt = nnwgt + 1
              endif
              a(ix,jy,kz,ln) = Max( real(nrx), a(ix,jy,kz,ln) )
            else
             if (  nrx .gt. a(ix,jy,kz,ln) .and. t1(ix,jy,kz) .gt. a(ix,jy,kz,ln) ) then
              if ( nrx .lt. t1(ix,jy,kz)  ) then
                ndbz = ndbz + 1
              else
                nmwgt = nmwgt + 1
                if ( t1(ix,jy,kz) .lt. ndbz ) nwlessthanz = nwlessthanz + 1
              endif
             else
              nnwgt = nnwgt + 1
             endif
              
              a(ix,jy,kz,ln) = Max(Min( real(nrx), t1(ix,jy,kz) ), a(ix,jy,kz,ln) )
            endif

           else ! } {
             if ( t1(ix,jy,kz) .gt. 0 .and. a(ix,jy,kz,ln) .gt. 0 ) then
              if ( t1(ix,jy,kz) .gt. a(ix,jy,kz,ln) ) then
                nmwgt = nmwgt + 1
              else
                nnwgt = nnwgt + 1
              endif
            endif
            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )
            nrx = a(ix,jy,kz,ln)



           endif ! }

           ! }
          else ! {
            if ( t1(ix,jy,kz) .gt. 0 .and. a(ix,jy,kz,ln) .gt. 0 ) then
              if ( t1(ix,jy,kz) .gt. a(ix,jy,kz,ln) ) then
                nmwgt = nmwgt + 1
              else
                nnwgt = nnwgt + 1
              endif
            endif
          endif! }
          
      enddo
      
      
      elseif ( l .eq. lr .and. imurain == 3) then

      xdn = 1000.
      
      do kz = 1,kze
          if (  t0(ix,jy,kz) .gt. 0. ) then

            vr = db(ix,kz)*a(ix,jy,kz,l)/(xdn*a(ix,jy,kz,ln))
            z = 3.6*(rnu+2.0)*a(ix,jy,kz,ln)*vr**2/(rnu+1.0)
          
             if ( z .gt. t0(ix,jy,kz) .and. z .gt. 0.0 .and.  &
     &          t0(ix,jy,kz) .gt. 0.0                         &
     &          .and. t0(ix,jy,kz) .gt. z0(ix,kz,l) ) then

            vr = db(ix,kz)*a(ix,jy,kz,l)/(xdn)
             chw =  a(ix,jy,kz,ln)
            nrx =   3.6*(rnu+2.0)*vr**2/((rnu+1.0)*t0(ix,jy,kz))
             if ( infall .eq. 3 ) then
              a(ix,jy,kz,ln) = Max( real(nrx), a(ix,jy,kz,ln) )
            elseif ( infall .eq. 4 ) then
              a(ix,jy,kz,ln) = Max( Min( real(nrx), t1(ix,jy,kz)), a(ix,jy,kz,ln) )
            endif

           else

            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )

           endif

          else

            a(ix,jy,kz,ln) = Max(t1(ix,jy,kz), a(ix,jy,kz,ln) )

          endif


      enddo

      endif

      return

      end subroutine calcnfromz1d


! ##############################################################################
! ##############################################################################
!
!  Subroutine to calculate number concentrations from initial state that has only mixing ratio.
!  Output N will be in #/m^3 in 'an' array, since sedimentation is done next.
!  Output ccw,cci etc. will be in #/kg

!
! 10.27.2015: Added hail calculation
!
      subroutine calcnfromq(nx,ny,nz,an,na,nor,norz,dn, &
     &  qcw,qci,qsw,qrw,qhw,qhl, &
     &  ccw,cci,csw,crw,chw,chl, &
     &  cccn,cccna, vhw,vhl,qv,spechum, invertccn_flag, cwmasin )
      

      
      implicit none

      integer nx,ny,nz,nor,norz,na,ngt,jgs,ixcol

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)  ! scalars (q, N, Z)

      real dn(nx,nz+1)  ! air density

      real, optional, dimension(nx,nz), intent(inout) :: qcw,qci,qsw,qrw,qhw,qhl, &
                                          ccw,cci,csw,crw,chw,chl, &
                                          cccn,cccna,vhw,vhl,qv, spechum
      logical, optional, intent(in) :: invertccn_flag
      real, optional :: cwmasin
      
      integer ixe,kze
      real    alpha
      real    qmin
      real    xvmn,xvmx
      integer ipconc
      integer lvol ! index for volume
      integer infall
      
      
      integer ix,jy,kz
      double precision vr,q,nrx,nrx2,rd,g1h,g1hl,g1r,g1s,zx,z,znew,zt,zxt,n1,laminv1
      double precision :: zr, zs, zh, dninv
      real, parameter :: xn0s = 3.0e6, xn0r = 8.0e6, xn0h = 2.0e5, xn0hl = 4.0e4
      real, parameter :: xdnr = 1000., xdns = 100. ,xdnh = 700.0, xdnhl = 900.0
      real, parameter :: zhlfac = 1./(pi*xdnhl*xn0hl)
      real, parameter :: zhfac = 1./(pi*xdnh*xn0h)
      real, parameter :: zrfac = 1./(pi*xdnr*xn0r)
      real, parameter :: zsfac = 1./(pi*xdns*xn0s)
      real, parameter :: g0 = (6.0)*(5.0)*(4.0)/((3.0)*(2.0)*(1.0))
      real, parameter :: xims=900.*0.523599*(2.*50.e-6)**3    ! mks   (100 micron diam solid sphere approx)
      real, parameter :: xgms=xdnh*0.523599*(300.e-6)**3    ! mks   (300 micron diam  sphere approx)
      real, parameter :: cwmas09 = 1000.*0.523599*(2.*9.e-6)**3 ! mass of 9-micron radius droplet

      real xv,xdn,cwmasinv
      integer :: ndbz, nmwgt, nnwgt, nwlessthanz
      double precision :: mixconv, mixconvqv, qsmax,qsmax2,qsmax3,qsmax4
      logical :: invertccn_local

! ------------------------------------------------------------------
      
      if ( present( invertccn_flag ) ) then
        invertccn_local = invertccn_flag
      else
        invertccn_local = .false.
      endif
      
      if ( present( cwmasin ) ) then
        cwmasinv = 1.0/cwmasin
      else
        cwmasinv = 1.0/cwmas09
      endif
      
      jy = 1
      
      
         g1h = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/  &
     &        ((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))

         g1hl = (6.0 + alphahl)*(5.0 + alphahl)*(4.0 + alphahl)/  &
     &        ((3.0 + alphahl)*(2.0 + alphahl)*(1.0 + alphahl))
     
         if ( imurain == 3 ) then
         g1r = (rnu+2.0)/(rnu+1.0)
         else ! imurain == 1
         g1r = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/  &
     &        ((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))
         endif

         g1s = (snu+2.0)/(snu+1.0)
      qsmax = 0
      qsmax2 = 0
      qsmax3 = 0
      qsmax4 = 0
!      IF ( .not. present( qcw ) ) THEN
      do kz = 1,nz
       do ix = 1,nx ! ixcol

!         qv_mp = spechum/(1.0_kind_phys-spechum)
!         IF ( convertdry ) THEN
!         qc_mp = qc/(1.0_kind_phys-spechum)
        mixconv = 1
        if ( present( spechum ) ) then ! convert to "dry" mixing ratios
          an(ix,jy,kz,lv)  = spechum(ix,kz)/(1.0d0 - spechum(ix,kz))
          mixconv = 1.0d0/(1.0d0 - spechum(ix,kz))
        else
          mixconv = 1.0d0
        endif
        if ( present( qv ) ) an(ix,jy,kz,lv) = qv(ix,kz) ! assume qv is "dry" mixing ratio if passed in
        if ( present( qcw ) ) an(ix,jy,kz,lc) = qcw(ix,kz)*mixconv
        if ( present( qrw ) ) an(ix,jy,kz,lr) = qrw(ix,kz)*mixconv
        if ( present( qci ) ) an(ix,jy,kz,li) = qci(ix,kz)*mixconv
        if ( present( qsw ) ) then
          an(ix,jy,kz,ls) = qsw(ix,kz)*mixconv
!          qsmax = Max( qsmax, qsw(ix,kz) )
!          qsmax2 = Max( qsmax2, an(ix,jy,kz,ls) )
        endif
        if ( present( qhw ) ) an(ix,jy,kz,lh) = qhw(ix,kz)*mixconv
        if ( lhl > 1 .and. present( qhl ) ) an(ix,jy,kz,lhl) = qhl(ix,kz)*mixconv
        if ( present( ccw ) ) an(ix,jy,kz,lnc) = ccw(ix,kz)*mixconv*dn(ix,kz)
        if ( present( crw ) ) an(ix,jy,kz,lnr) = crw(ix,kz)*mixconv*dn(ix,kz)
        if ( present( cci ) ) an(ix,jy,kz,lni) = cci(ix,kz)*mixconv*dn(ix,kz)
        if ( present( csw ) ) an(ix,jy,kz,lns) = csw(ix,kz)*mixconv*dn(ix,kz)
        if ( present( chw ) ) an(ix,jy,kz,lnh) = chw(ix,kz)*mixconv*dn(ix,kz)
        if ( lhl > 1 .and. present( chl ) ) an(ix,jy,kz,lnhl) = chl(ix,kz)*mixconv*dn(ix,kz)
        if ( lvh > 1 .and. present( vhw ) ) an(ix,jy,kz,lvh) = vhw(ix,kz)*mixconv
        if ( lvhl > 1 .and. present( vhl ) ) an(ix,jy,kz,lvhl) = vhl(ix,kz)*mixconv
        if ( lccn > 1 .and. present( cccn ) ) an(ix,jy,kz,lccn) = cccn(ix,kz)*mixconv*dn(ix,kz)
        if ( lccna > 1 .and. present( cccna ) ) an(ix,jy,kz,lccna) = cccna(ix,kz)*mixconv

         dninv = 1./dn(ix,kz)
         
!        IF ( .not. present( qcw ) ) THEN
   !  Cloud droplets
         
         if ( lnc > 1 ) then
           if ( an(ix,jy,kz,lnc) <= cxmin .and. an(ix,jy,kz,lc) > qxmin_init(lc) ) then
             
             an(ix,jy,kz,lnc) = Min(qccn, an(ix,jy,kz,lc)*cwmasinv )*dn(ix,kz)
             
             if ( invertccn_local ) then
               an(ix,jy,kz,lccn) = an(ix,jy,kz,lccn) + an(ix,jy,kz,lnc)
             else
             
             if ( lccn > 1 .and. lccna < 1 ) then
                an(ix,jy,kz,lccn) = an(ix,jy,kz,lccn) - an(ix,jy,kz,lnc)
             endif
             if ( lccna > 1 ) then
                an(ix,jy,kz,lccna) = an(ix,jy,kz,lccna) + an(ix,jy,kz,lnc)
             endif
             endif

           elseif ( an(ix,jy,kz,lc) <= qxmin(lc) .or.  &
                    ( an(ix,jy,kz,lnc) <= cxmin .and. an(ix,jy,kz,lc) <= qxmin_init(lc)) ) then
             
             an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
             an(ix,jy,kz,lnc) = 0.0
             an(ix,jy,kz,lc) = 0.0
           
           endif
         endif

   !  Cloud ice
         
         if ( lni > 1 ) then
           if ( an(ix,jy,kz,lni) <= cxmin .and. an(ix,jy,kz,li) > qxmin_init(li) ) then
             an(ix,jy,kz,lni) = dn(ix,kz)*an(ix,jy,kz,li)/xims
           
           elseif ( an(ix,jy,kz,li) <= qxmin(li) .or. &
                    ( an(ix,jy,kz,lni) <= cxmin .and. an(ix,jy,kz,li) <= qxmin_init(li)) ) then
             an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
             an(ix,jy,kz,lni) = 0.0
             an(ix,jy,kz,li) = 0.0
           endif
         endif

   !  rain
         
         if ( lnr > 1 ) then
           if ( an(ix,jy,kz,lnr) <= 0.1*cxmin .and. an(ix,jy,kz,lr) > qxmin_init(lr) ) then

             q = an(ix,jy,kz,lr)
             
             laminv1 = (dn(ix,kz) * q * zrfac)**(0.25)  ! inverse of slope
             
             n1 = laminv1*xn0r  ! number concentration for inv. exponential single moment input
             
             nrx =  n1*g1r/g0   ! number concentration for different shape parameter

             an(ix,jy,kz,lnr) = nrx ! *dninv ! convert to number mixing ratio
             
           elseif ( an(ix,jy,kz,lr) <= qxmin(lr) .or. &
                    ( an(ix,jy,kz,lnr) <= cxmin .and. an(ix,jy,kz,lr) <= qxmin_init(lr)) ) then
             an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lr)
             an(ix,jy,kz,lnr) = 0.0
             an(ix,jy,kz,lr) = 0.0
           endif
         endif

             if ( lzr > 1 ) then ! set reflectivity moment
               if ( an(ix,jy,kz,lr) > qxmin_init(lr) .and. an(ix,jy,kz,lzr) < zxmin .and. &
                    an(ix,jy,kz,lnr) > cxmin ) then
                  q = an(ix,jy,kz,lr)
                  nrx = an(ix,jy,kz,lnr)
                  an(ix,jy,kz,lzr) = 36.*g1r*dn(ix,kz)**2*q**2/(pi**2*xdnr**2*nrx) ! *dninv
               endif
             endif

  ! snow
         if ( lns > 1 ) then
           if ( an(ix,jy,kz,lns) <= 0.1*cxmin .and. an(ix,jy,kz,ls) > qxmin_init(ls) ) then

             q = an(ix,jy,kz,ls)
             
             laminv1 = (dn(ix,kz) * q * zsfac)**(0.25)  ! inverse of slope
             
             n1 = laminv1*xn0s  ! number concentration for inv. exponential single moment input
             
             nrx =  n1*g1s/g0   ! number concentration for different shape parameter

             an(ix,jy,kz,lns) = nrx ! *dninv ! convert to number mixing ratio

           elseif ( an(ix,jy,kz,ls) <= qxmin(ls) .or. &
                    ( an(ix,jy,kz,lns) <= cxmin .and. an(ix,jy,kz,ls) <= qxmin_init(ls)) ) then
             an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
             an(ix,jy,kz,lns) = 0.0
             an(ix,jy,kz,ls) = 0.0
             
           endif
         endif
         
    ! graupel

         if ( lnh > 1 ) then
           if ( an(ix,jy,kz,lnh) <= 0.1*cxmin .and. an(ix,jy,kz,lh) > qxmin_init(lh) ) then
             if ( lvh > 1 ) then
               if ( an(ix,jy,kz,lvh) <= 0.0 ) then
                 an(ix,jy,kz,lvh) = an(ix,jy,kz,lh)/xdnh
               endif
             endif

             q = an(ix,jy,kz,lh)
             
             laminv1 = (dn(ix,kz) * q * zhfac)**(0.25)  ! inverse of slope
             
             n1 = laminv1*xn0h  ! number concentration for inv. exponential single moment input
             
             nrx =  n1*g1h/g0   ! number concentration for different shape parameter

             nrx2 = dn(ix,kz) * q / xgms
             
             nrx = Min( nrx, nrx2 )

             if ( nrx > cxmin ) then
               an(ix,jy,kz,lnh) = nrx ! *dninv ! convert to number mixing ratio
             else
               an(ix,jy,kz,lh) = 0.0
               an(ix,jy,kz,lnh) = 0.0
               an(ix,jy,kz,lvh) = 0.0
             endif

           elseif ( an(ix,jy,kz,lh) <= qxmin(lh) .or. &
                    ( an(ix,jy,kz,lnh) <= cxmin .and. an(ix,jy,kz,lh) <= qxmin_init(lh)) ) then
           
              an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lh)
              an(ix,jy,kz,lh) = 0.0
           
           endif
         endif

             if ( lzh > 1 ) then ! set reflectivity moment
               if ( an(ix,jy,kz,lh) > qxmin_init(lh) .and. an(ix,jy,kz,lzh) < zxmin .and. &
                    an(ix,jy,kz,lnh) > cxmin ) then
                  q = an(ix,jy,kz,lh)
                  nrx = an(ix,jy,kz,lnh)
                  an(ix,jy,kz,lzh) = 36.*g1h*dn(ix,kz)**2*q**2/(pi**2*xdnh**2*nrx) ! *dninv
               endif
             endif

    ! hail

         if ( lnhl > 1 .and. lhl > 1 ) then
           if ( an(ix,jy,kz,lnhl) <= 0.1*cxmin .and. an(ix,jy,kz,lhl) > qxmin_init(lhl) ) then
             if ( lvhl > 1 ) then
               if ( an(ix,jy,kz,lvhl) <= 0.0 ) then
                 an(ix,jy,kz,lvhl) = an(ix,jy,kz,lhl)/xdnhl
               endif
             endif

             q = an(ix,jy,kz,lhl)
             
             laminv1 = (dn(ix,kz) * q * zhlfac)**(0.25)  ! inverse of slope
             
             n1 = laminv1*xn0hl  ! number concentration for inv. exponential single moment input
             
             nrx =  n1*g1hl/g0   ! number concentration for different shape parameter

             an(ix,jy,kz,lnhl) = nrx ! *dninv ! convert to number mixing ratio

           elseif ( an(ix,jy,kz,lhl) <= qxmin(lhl) .or.  &
                   ( an(ix,jy,kz,lnhl) <= cxmin .and. an(ix,jy,kz,lhl) <= qxmin_init(lhl)) ) then
           
              an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lhl)
              an(ix,jy,kz,lhl) = 0.0
           
           endif
         endif

             if ( lzhl > 1 ) then ! set reflectivity moment
               if ( an(ix,jy,kz,lhl) > qxmin_init(lhl) .and. an(ix,jy,kz,lzhl) < zxmin .and. &
                    an(ix,jy,kz,lnhl) > cxmin ) then
                  q = an(ix,jy,kz,lhl)
                  nrx = an(ix,jy,kz,lnhl)
                  an(ix,jy,kz,lzhl) = 36.*g1hl*dn(ix,kz)**2*q**2/(pi**2*xdnhl**2*nrx) ! *dninv
               endif
             endif
         
         
!         ENDIF

!         spechum = qv_mp/(1.0_kind_phys+qv_mp)
!         IF ( convertdry ) THEN
!         qc      = qc_mp/(1.0_kind_phys+qv_mp)
        mixconvqv = 1
        if ( present( spechum ) ) then ! convert back to "dry+vapor" mixing ratios
          !an(ix,jy,kz,lv)  = spechum(ix,kz)/(1.0d0 - spechum(ix,kz))
          mixconvqv = 1.0d0/(1.0d0 + an(ix,jy,kz,lv))
          spechum(ix,kz) = an(ix,jy,kz,lv)*mixconvqv
        else
          mixconvqv = 1.0d0
        endif

        if ( present( qv ) )  qv(ix,kz)  = an(ix,jy,kz,lv)
        if ( present( qcw ) ) qcw(ix,kz) = an(ix,jy,kz,lc)*mixconvqv
        if ( present( qrw ) ) qrw(ix,kz) = an(ix,jy,kz,lr)*mixconvqv
        if ( present( qci ) ) qci(ix,kz) = an(ix,jy,kz,li)*mixconvqv
        if ( present( qsw ) ) then
          qsw(ix,kz) = an(ix,jy,kz,ls)*mixconvqv
!          qsmax3 = Max( qsmax3, qsw(ix,kz) )
!          qsmax4 = Max( qsmax4, an(ix,jy,kz,ls) )
        endif
        if ( present( qhw ) ) qhw(ix,kz) = an(ix,jy,kz,lh)*mixconvqv
        if ( lhl > 1 .and. present( qhl ) ) qhl(ix,kz) = an(ix,jy,kz,lhl)*mixconvqv
        if ( present( ccw ) ) ccw(ix,kz) = an(ix,jy,kz,lnc)*mixconvqv*dninv
        if ( present( crw ) ) crw(ix,kz) = an(ix,jy,kz,lnr)*mixconvqv*dninv
        if ( present( cci ) ) cci(ix,kz) = an(ix,jy,kz,lni)*mixconvqv*dninv
        if ( present( csw ) ) csw(ix,kz) = an(ix,jy,kz,lns)*mixconvqv*dninv
        if ( present( chw ) ) chw(ix,kz) = an(ix,jy,kz,lnh)*mixconvqv*dninv
        if ( lhl > 1 .and. present( chl ) ) chl(ix,kz) = an(ix,jy,kz,lnhl)*mixconvqv*dninv
        if ( lvh > 1 .and. present( vhw ) ) vhw(ix,kz) = an(ix,jy,kz,lvh)*mixconvqv
        if ( lvhl > 1 .and. present( vhl ) ) vhl(ix,kz) = an(ix,jy,kz,lvhl)*mixconvqv
        if ( lccn > 1 .and. present( cccn ) ) cccn(ix,kz) = an(ix,jy,kz,lccn)*mixconvqv*dninv
        if ( lccna > 1 .and. present( cccna ) ) cccna(ix,kz) = an(ix,jy,kz,lccna)*mixconvqv


      enddo ! ix
      enddo ! kz
!      ELSE
!        write(0,*) 'calcnfromq: lv = ',lv,lc,lr,li,ls,lh,lvh,lhl,lccn,lccna
!        write(0,*) 'calcnfromq: nx,ny,nz,na = ',nx,ny,nz,na
!      
!      ENDIF
      
!      IF ( present( qsw ) ) THEN
!      write(0,*) 'calcnfromq: qsmax = ',qsmax,qsmax2,qsmax3,qsmax4
!      ENDIF
      
      return
      
      end subroutine calcnfromq

! ##############################################################################
! ##############################################################################
!
!  Subroutine to calculate number concentrations from convection parameterization rates that have only mixing ratio.
!  N will be in #/kg, NOT #/m^3, since sedimentation is done next.
!

!
! 10.27.2015: Added hail calculation
!
      subroutine calcnfromcuten(nx,ny,nz,an,anold,na,nor,norz,dn)

      
      implicit none

      integer nx,ny,nz,nor,norz,na,ngt,jgs,ixcol

      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)  ! scalars (q, N, Z) from CUTEN arrays
      real anold(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)  ! scalars (q, N, Z)

      real dn(nx,nz+1)  ! air density
      
      integer ixe,kze
      real    alpha
      real    qmin
      real    xvmn,xvmx
      integer ipconc
      integer lvol ! index for volume
      integer infall
      
      
      integer ix,jy,kz
      double precision vr,q,nrx,rd,g1h,g1hl,g1r,g1s,zx,chw,z,znew,zt,zxt,n1,laminv1
      double precision :: zr, zs, zh, dninv
      real, parameter :: xn0s = 3.0e6, xn0r = 8.0e6, xn0h = 4.0e4, xn0hl = 4.0e4
      real, parameter :: xdnr = 1000., xdns = 100. ,xdnh = 700.0, xdnhl = 900.0
      real, parameter :: zhlfac = 1./(pi*xdnhl*xn0hl)
      real, parameter :: zhfac = 1./(pi*xdnh*xn0h)
      real, parameter :: zrfac = 1./(pi*xdnr*xn0r)
      real, parameter :: zsfac = 1./(pi*xdns*xn0s)
      real, parameter :: g0 = (6.0)*(5.0)*(4.0)/((3.0)*(2.0)*(1.0))
      real, parameter :: xims=900.*0.523599*(2.*50.e-6)**3    ! mks   (100 micron diam solid sphere approx)
      real, parameter :: xcms=1000.*0.523599*(2.*7.5e-6)**3    ! mks   (100 micron diam solid sphere approx)

      real :: xmass,xv,xdn
      integer :: ndbz, nmwgt, nnwgt, nwlessthanz

! ------------------------------------------------------------------
      
      
      jy = 1
      
      
         g1h = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/  &
     &        ((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))

         g1hl = (6.0 + alphahl)*(5.0 + alphahl)*(4.0 + alphahl)/  &
     &        ((3.0 + alphahl)*(2.0 + alphahl)*(1.0 + alphahl))
     
         if ( imurain == 3 ) then
         g1r = (rnu+2.0)/(rnu+1.0)
         else ! imurain == 1
         g1r = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/  &
     &        ((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))
         endif

         g1s = (snu+2.0)/(snu+1.0)
      
      do kz = 1,nz
       do ix = 1,nx ! ixcol

         dninv = 1./dn(ix,kz)
         
   !  Cloud droplets
         
         if ( lnc > 1 ) then
!           IF ( an(ix,jy,kz,lnc) <= 0.1*cxmin .and. an(ix,jy,kz,lc) > qxmin(lc) ) THEN
           if ( an(ix,jy,kz,lnc) > qxmin(lc) ) then
             anold(ix,jy,kz,lnc) = anold(ix,jy,kz,lnc) + an(ix,jy,kz,lc)/xcms
           endif
         endif

   !  Cloud ice
         
         if ( lni > 1 ) then
           if ( an(ix,jy,kz,lni) > qxmin(li) ) then
             anold(ix,jy,kz,lni) = anold(ix,jy,kz,lni) + an(ix,jy,kz,li)/xims
           endif
         endif

   !  rain
         
         if ( lnr > 1 ) then
           if ( an(ix,jy,kz,lr) > qxmin(lr) ) then ! adding rain mass from CU scheme

            if ( .true. .or. (anold(ix,jy,kz,lr) - an(ix,jy,kz,lr)) < qxmin(lr) .or. anold(ix,jy,kz,lnr) < cxmin ) then 

             q = an(ix,jy,kz,lr)
             
             laminv1 = (dn(ix,kz) * q * zrfac)**(0.25)  ! inverse of slope
             
             n1 = laminv1*xn0r  ! number concentration for inv. exponential single moment input
             
             nrx =  n1*g1r/g0   ! number concentration for different shape parameter

             anold(ix,jy,kz,lnr) = anold(ix,jy,kz,lnr) + nrx ! *dninv ! convert to number mixing ratio

            else
             ! assume mean particle mass of pre-existing snow
                xmass = anold(ix,jy,kz,lr)/anold(ix,jy,kz,lnr)
                anold(ix,jy,kz,lnr) = anold(ix,jy,kz,lnr) + an(ix,jy,kz,lr)/xmass
            endif
             
             if ( lzr > 1 ) then ! set reflectivity moment
               an(ix,jy,kz,lzr) = 36.*g1r*dn(ix,kz)**2*q**2/(pi**2*xdnr**2*nrx) ! *dninv
             endif
           endif
         endif

  ! snow
         if ( lns > 1 ) then
           if ( an(ix,jy,kz,ls) > qxmin(ls) ) then ! adding snow mass from CU scheme

             if ( .true. .or. (anold(ix,jy,kz,ls) - an(ix,jy,kz,ls)) < qxmin(ls) .or. anold(ix,jy,kz,lns) < cxmin ) then 
             
             ! assume that there was no snow before this
             
             q = an(ix,jy,kz,ls)
             
             laminv1 = (dn(ix,kz) * q * zsfac)**(0.25)  ! inverse of slope
             
             n1 = laminv1*xn0s  ! number concentration for inv. exponential single moment input
             
             nrx =  n1*g1s/g0   ! number concentration for different shape parameter

             anold(ix,jy,kz,lns) = anold(ix,jy,kz,lns) + nrx ! *dninv ! convert to number mixing ratio
             
             else
             ! assume mean particle mass of pre-existing snow
                xmass = anold(ix,jy,kz,ls)/anold(ix,jy,kz,lns)
                anold(ix,jy,kz,lns) = anold(ix,jy,kz,lns) + an(ix,jy,kz,ls)/xmass
             endif
             
           endif
         endif
         
    ! graupel

!         IF ( lnh > 1 ) THEN
!           IF ( an(ix,jy,kz,lnh) <= 0.1*cxmin .and. an(ix,jy,kz,lh) > qxmin(lh) ) THEN
!             IF ( lvh > 1 ) THEN
!               IF ( an(ix,jy,kz,lvh) <= 0.0 ) THEN
!                 an(ix,jy,kz,lvh) = an(ix,jy,kz,lh)/xdnh
!               ENDIF
!             ENDIF
!
!             q = an(ix,jy,kz,lh)
!             
!             laminv1 = (dn(ix,kz) * q * zhfac)**(0.25)  ! inverse of slope
!             
!             n1 = laminv1*xn0h  ! number concentration for inv. exponential single moment input
!             
!             nrx =  n1*g1h/g0   ! number concentration for different shape parameter
!
!             an(ix,jy,kz,lnh) = nrx ! *dninv ! convert to number mixing ratio
!
!             IF ( lzh > 1 ) THEN ! set reflectivity moment
!               an(ix,jy,kz,lzh) = 36.*g1h*dn(ix,kz)**2*q**2/(pi**2*xdnh**2*nrx) ! *dninv
!             ENDIF
!           ENDIF
!         ENDIF
!
!    ! hail
!
!         IF ( lnhl > 1 .and. lhl > 1 ) THEN
!           IF ( an(ix,jy,kz,lnhl) <= 0.1*cxmin .and. an(ix,jy,kz,lhl) > qxmin(lhl) ) THEN
!             IF ( lvhl > 1 ) THEN
!               IF ( an(ix,jy,kz,lvhl) <= 0.0 ) THEN
!                 an(ix,jy,kz,lvhl) = an(ix,jy,kz,lhl)/xdnhl
!               ENDIF
!             ENDIF
!
!             q = an(ix,jy,kz,lhl)
!             
!             laminv1 = (dn(ix,kz) * q * zhlfac)**(0.25)  ! inverse of slope
!             
!             n1 = laminv1*xn0hl  ! number concentration for inv. exponential single moment input
!             
!             nrx =  n1*g1hl/g0   ! number concentration for different shape parameter
!
!             an(ix,jy,kz,lnhl) = nrx ! *dninv ! convert to number mixing ratio
!
!             IF ( lzhl > 1 ) THEN ! set reflectivity moment
!               an(ix,jy,kz,lzhl) = 36.*g1hl*dn(ix,kz)**2*q**2/(pi**2*xdnhl**2*nrx) ! *dninv
!             ENDIF
!           ENDIF
!         ENDIF
 
      enddo ! ix
      enddo ! kz
      
      return
      
      end subroutine calcnfromcuten

! #####################################################################
! #####################################################################

   subroutine calc_eff_radius    &
     &  (nx,ny,nz,na,jyslab & 
     &  ,nor,norz & 
     &  ,t1,t2,t3,t4,t5,t6, f_t5,f_t6  & 
     &  ,qcw,qci,qsw,qrw &
     &  ,ccw,cci,csw,crw &
     &  ,an,dn )

   implicit none

      integer, parameter :: ng1 = 1
      integer :: nx,ny,nz,na
      integer :: ng
      integer :: nor,norz, jyslab ! ,nht,ngt,igsr
      real    :: dtp  ! time step


!
! external temporary arrays
!

      real,optional :: t1(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real,optional :: t2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real,optional :: t3(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real,optional :: t4(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real,optional :: t5(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real,optional :: t6(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      logical, optional :: f_t5, f_t6 ! flags to fill t5/t6 for graupel/hail

      real, optional :: an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real, optional, dimension(nx,nz) :: qcw,qci,qsw,qrw,ccw,cci,csw,crw 

      


      
    ! local
    
      real pb(-norz+ng1:nz+norz)
      real pinit(-norz+ng1:nz+norz)

! 
!  declarations microphysics and for gather/scatter
!
      integer nxmpb,nzmpb,nxz
      integer mgs,ngs,numgs,inumgs
      parameter (ngs=1)
      integer ngscnt,igs(ngs),kgs(ngs)
      real rho0(ngs)

      integer ix,kz,i,n, kp1
      integer :: jy, jgs
      integer ixb,ixe,jyb,jye,kzb,kze
    
      integer itile,jtile,ktile
      integer ixend,jyend,kzend,kzbeg
      integer nxend,nyend,nzend,nzbeg

      real :: qx(ngs,lv:lhab)
      real :: cx(ngs,lc:lhab)
      real :: xv(ngs,lc:lhab)
      real :: xmas(ngs,lc:lhab)
      real :: xdn(ngs,lc:lhab)
      real :: xdia(ngs,lc:lhab,3)
      real :: alpha(ngs,lc:lhab)
      
      real :: gamc1,gamc2,gami1,gami2,gams1,gams2,gamr1,gamr2
      real :: factor_c, factor_i, factor_s, factor_r
      real :: lam_c, lam_i, lam_s, lam_r
      integer :: il


! -------------------------------------------------------------------------------
      itile = nx
      jtile = ny
      ktile = nz
      ixend = nx
      jyend = ny
      kzend = nz
      nxend = nx + 1
      nyend = ny + 1
      nzend = nz
      kzbeg = 1
      nzbeg = 1

       jy = 1
       pb(:) = 0.0
       pinit(:) = 0.0

     gamc1 = Gamma_sp(2. + cnu)
     gamc2 = 1. ! Gamma[1 + alphac]
     gami1 = Gamma_sp(2. + cinu)
     gami2 = 1. ! Gamma[1 + alphac]
     gams1 = Gamma_sp(2. + snu)
     gams2 = Gamma_sp(1. + snu)
     gamr1 = Gamma_sp(2. + rnu)
     gamr2 = Gamma_sp(1. + rnu)

     factor_c = (1. + cnu)*Gamma_sp(1. + cnu)/Gamma_sp(5./3. + cnu)
     factor_i = (1. + cinu)*Gamma_sp(1. + cinu)/Gamma_sp(5./3. + cinu)
     factor_s = (1. + snu)*Gamma_sp(1. + snu)/Gamma_sp(5./3. + snu)

     if ( present(t4) ) then
     if ( imurain == 3 ) then
       factor_r = (1. + rnu)*Gamma_sp(1. + rnu)/Gamma_sp(5./3. + rnu)
     else
       factor_r = ((Pi*(alphar+3.)*(alphar+1.)*(alphar+1.))/6.)**(1./3.)
     endif
     endif

!
!     jy = 1 ! working on a 2d slab
!!  VERY IMPORTANT:  SET jgs = jy

      jgs = jy

      mgs = 1
      do kz = 1,nz
       do ix = 1,nx ! ixcol

         rho0(mgs) = dn(ix,jy,kz)
         if ( present( an ) ) then
         do il = lc,ls
          qx(mgs,il) = max(an(ix,jy,kz,il), 0.0) 
          cx(mgs,il) = max(an(ix,jy,kz,ln(il)), 0.0) 
         enddo
         else
          qx(mgs,:) = 0.0
          cx(mgs,:) = 0.0
          if ( present(qcw) ) qx(mgs,lc) = qcw(ix,kz)
          if ( present(qci) ) qx(mgs,li) = qci(ix,kz)
          if ( present(qsw) ) qx(mgs,ls) = qsw(ix,kz)
          if ( present(qrw) ) qx(mgs,lr) = qrw(ix,kz)
          if ( present(ccw) ) cx(mgs,lc) = ccw(ix,kz)*rho0(mgs)
          if ( present(cci) ) cx(mgs,li) = cci(ix,kz)*rho0(mgs)
          if ( present(csw) ) cx(mgs,ls) = csw(ix,kz)*rho0(mgs)
          if ( present(crw) ) cx(mgs,lr) = crw(ix,kz)*rho0(mgs)
         
         endif
         
         if ( present( t1 ) .and. qx(mgs,lc) > qxmin(lc) .and. cx(mgs,lc) > cxmin ) then
! Lambda for cloud droplets 
         lam_c = ((cx(mgs,lc)*(Pi/6.)*xdn0(lc)*Gamc1)/(qx(mgs,lc)*rho0(mgs)*Gamc2))**(1./3.)
          t1(ix,jy,kz) = 0.5*factor_c/lam_c
         endif

         if ( present( t2 ) .and. qx(mgs,li) > qxmin(li) .and. cx(mgs,li) > cxmin ) then
! Lambda for cloud ice 
         lam_i = ((cx(mgs,li)*(Pi/6.)*xdn0(li)*Gami1)/(qx(mgs,li)*rho0(mgs)*Gami2))**(1./3.)
          t2(ix,jy,kz) = 0.5*factor_i/lam_i
         endif

         if ( present( t3 ) .and. qx(mgs,ls) > qxmin(ls) .and. cx(mgs,ls) > cxmin  ) then
! Lambda for snow
         lam_s = ((cx(mgs,ls)*(Pi/6.)*xdn0(ls)*Gams1)/(qx(mgs,ls)*rho0(mgs)*Gams2))**(1./3.)
          t3(ix,jy,kz) = 0.5*factor_s/lam_s
         endif

         if ( present( t4 ) .and. present(qrw) .and. present(crw) ) then
         if ( qx(mgs,lr) > Max(1.e-8,qxmin(lr))  .and. cx(mgs,lr) > cxmin ) then
           if ( imurain == 1 ) then ! gamma-diameter
! Lambda for rain
             lam_r = factor_r *((xdn0(lr)*cx(mgs,lr))/(qx(mgs,lr)*rho0(mgs)))**(1./3.)
             t4(ix,jy,kz) = 0.5*(alphar+3.)/lam_r
           else ! gamma-volume
! Lambda for rain
             lam_r = ((cx(mgs,lr)*(Pi/6.)*xdn0(lr)*Gamr1)/(qx(mgs,lr)*rho0(mgs)*Gamr2))**(1./3.)
             t4(ix,jy,kz) = 0.5*factor_r/lam_r
           endif
         endif
         endif

      
       enddo ! ix
      enddo ! kz

   return
   end subroutine calc_eff_radius


! #####################################################################
! #####################################################################

      subroutine QVEXCESS(ngs,mgs,qwvp0,qv0,qcw1,pres,thetap0,theta0, &
     &    qvex,pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felvcp,ss1,pk,ngscnt)
      
!#####################################################################
!  Purpose: find the amount of vapor that can be condensed to liquid
!#####################################################################

      implicit none

      integer ngs,mgs,ngscnt
      
      real theta2temp
      
      real qvex
      
      integer nqsat
      real fqsat, cbw
      
      real ss1  ! 'target' supersaturation
!
!  input arrays
!
      real qv0(ngs), qcw1(ngscnt), pres(ngs), qwvp0(mgs)
      real thetap0(ngs), theta0(ngs)
      real fcqv1(ngs), felvcp(ngs), pi0(ngs)
      real pk(ngs)
      
      real tabqvs(nqsat)
!
! Local stuff
!
      
      integer itertd
      integer ltemq
      real gamss
      real theta(ngs), qvap(ngs), pqs(ngs), qcw(ngs), qwv(ngs)
      real qcwtmp(ngs), qss(ngs), qvs(ngs), qwvp(ngs)
      real dqcw(ngs), dqwv(ngs), dqvcnd(ngs)
      real temg(ngs), temcg(ngs), thetap(ngs)
      
      real tfr
      parameter ( tfr = 273.15 )
            
!      real poo,cap
!      parameter ( cap = rd/cp, poo = 1.0e+05 )
!
!
!  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
!
!
!
!  set up temperature and vapor arrays
!
      pqs(mgs) = (380.0)/(pres(mgs))
      thetap(mgs) = thetap0(mgs)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      qwvp(mgs) = qwvp0(mgs)
      qvap(mgs) = max( (qwvp0(mgs) + qv0(mgs)), 0.0 )
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
!      temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
!
!
!
!  reset temporaries for cloud particles and vapor
!
      
      qwv(mgs) = max( 0.0, qvap(mgs) )
      qcw(mgs) = max( 0.0, qcw1(mgs) )
!
!
      qcwtmp(mgs) = qcw(mgs)
      temcg(mgs) = temg(mgs) - tfr
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qss(mgs) = (0.01*ss1 + 1.0)*qvs(mgs)
!
!  iterate  adjustment
!
      do itertd = 1,2
!
!
!  calculate super-saturation
!
      dqcw(mgs) = 0.0
      dqwv(mgs) = ( qwv(mgs) - qss(mgs) )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then           !  subsaturated
        if( qcw(mgs) .gt. -dqwv(mgs) ) then  ! check if qc can make up all of the deficit
          dqcw(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                 !  otherwise make all qc available for evap
          dqcw(mgs) = -qcw(mgs)
          dqwv(mgs) = dqwv(mgs) + qcw(mgs)
        end if
!
        qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs)  )  ! add to perturbation vapor
!
        qcw(mgs) = qcw(mgs) + dqcw(mgs)

        thetap(mgs) = thetap(mgs) +  &
     &                1./pi0(mgs)*  &
     &                (felvcp(mgs)*dqcw(mgs) )

      end if  ! dqwv(mgs) .lt. 0. (end of evap/sublim)
!
! condensation/deposition
!
      if ( dqwv(mgs) .ge. 0. ) then
!
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/  &
     &  ((temg(mgs)-cbw)**2))
!
!
      dqcw(mgs) = dqvcnd(mgs)
!
      thetap(mgs) = thetap(mgs) +  &
     &   (felvcp(mgs)*dqcw(mgs) )    &
     & / (pi0(mgs))
      qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )
      qcw(mgs) = qcw(mgs) + dqcw(mgs)
!
      end if !  dqwv(mgs) .ge. 0.

      theta(mgs) = thetap(mgs) + theta0(mgs)
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
!      temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
      qvap(mgs) = Max((qwvp(mgs) + qv0(mgs)), 0.0)
      temcg(mgs) = temg(mgs) - tfr
!      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qcw(mgs) = max( 0.0, qcw(mgs) )
      qwv(mgs) = max( 0.0, qvap(mgs))
      qss(mgs) = (0.01*ss1 + 1.0)*qvs(mgs)
      end do
!
!  end the saturation adjustment iteration loop
!
!
      qvex = Max(0.0, qcw(mgs) - qcw1(mgs) )

      return
      end subroutine QVEXCESS

! #####################################################################
! #####################################################################





!
! ##############################################################################
!
      subroutine setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,cnostmp, &
     &                 xmas,vtxbar,xdn,xvmn0,xvmx0,xv,cdx,cdxgs,            &
     &                 ipconc1,ndebug1,ngs,nz,kgs,fadvisc,   &
     &                 cwmasn,cwmasx,cwradn,cnina,cimna,cimxa,      &
     &                 itype1a,itype2a,temcg,infdo,alpha,ildo,axx,bxx)
!     &                 itype1a,itype2a,temcg,infdo,alpha,ildo,axh,bxh,axhl,bxhl)


      implicit none
      
      integer ngscnt,ngs0,ngs,nz
!      integer infall    ! whether to calculate number-weighted fall speeds
      
      real xv(ngs,lc:lhab)
      real qx(ngs,lv:lhab)
      real qxw(ngs,ls:lhab)
      real cx(ngs,lc:lhab)
      real vtxbar(ngs,lc:lhab,3)
      real xmas(ngs,lc:lhab)
      real xdn(ngs,lc:lhab)
      real cdxgs(ngs,lc:lhab)
      real xdia(ngs,lc:lhab,3)
      real xvmn0(lc:lhab), xvmx0(lc:lhab)
      real qxmin(lc:lhab)
      real cdx(lc:lhab)
      real alpha(ngs,lc:lhab)
      
      real rho0(ngs),rhovt(ngs),temcg(ngs)
      real cno(lc:lhab)
      real cnostmp(ngs)
      
      real cwc1, cimna, cimxa
      real cnina(ngs)
      integer kgs(ngs)
      real fadvisc(ngs)
      real fsw
      
      integer ipconc1
      integer ndebug1
      
      integer, intent (in) :: itype1a,itype2a,infdo
      integer, intent (in) :: ildo ! which species to do, or all if ildo=0

      real :: axx(ngs,lh:lhab),bxx(ngs,lh:lhab)
!!      real :: axh(ngs),bxh(ngs)
!      real :: axhl(ngs),bxhl(ngs)
      
! Local vars

      
      
      real swmasmx, dtmp
      real cd
      real cwc0 ! ,cwc1
      real :: cwch(ngscnt), cwchl(ngscnt)
      real :: cwchtmp,cwchltmp,xnutmp
      real pii
      real cimasx,cimasn
      real cwmasn,cwmasx,cwradn
      real cwrad
      real vr,rnux
      real alp
      
      real ccimx

      integer mgs
      
      real arx,frx,vtrain,fw
      real fwlo,fwhi,rfwdiff
      real ar,br,cs,ds
!      real gf4p5, gf4ds, gf4br, ifirst, gf1ds
!      real gfcinu1, gfcinu1p47, gfcinu2p47
      real gr
      real rwrad,rwdia
      real mwfac
      integer il

!      save gf4p5, gf4ds, gf4br, ifirst, gf1ds
!      save gfcinu1, gfcinu1p47, gfcinu2p47
!      data ifirst /0/
      
      real bta1,cnit
      parameter ( bta1 = 0.6, cnit = 1.0e-02 )
      real x,y,tmp,del
      real aax,bbx,delrho
      integer :: indxr
      real mwt, nwt, zwt
      real, parameter :: rho00 = 1.225
      integer i
      real xvbarmax

      integer l1, l2


!
! set values
!
!      cwmasn = 5.23e-13  ! radius of 5.0e-6
!      cwradn = 5.0e-6
!      cwmasx = 5.25e-10  ! radius of 50.0e-6

      fwlo = 0.2                ! water fraction to start weighting toward rain fall speed
      fwhi = 0.4                ! water fraction at which rain fall speed only is used
      rfwdiff = 1./(fwhi - fwlo)
      
!      pi = 4.0*atan(1.0)
      pii = piinv ! 1.0/pi

      arx = 10.
      frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.

      ar = 841.99666  
      br = 0.8
      gr = 9.8
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42

      if ( ildo == 0 ) then
        l1 = lc
        l2 = lhab
      else
        l1 = ildo
        l2 = ildo
      endif

!      IF ( ifirst .eq. 0 ) THEN
!        ifirst = 1
!        gf4br = gamma(4.0+br)
!        gf4ds = gamma(4.0+ds)
!!        gf1ds = gamma(1.0+ds)
!        gf4p5 = gamma(4.0+0.5)
!        gfcinu1 = gamma(cinu + 1.0)
!        gfcinu1p47 = gamma(cinu + 1.47167)
!        gfcinu2p47 = gamma(cinu + 2.47167)
        
        if ( lh  .gt. 1 ) then
          if ( dmuh == 1.0 ) then
            cwchtmp = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
          else
            cwchtmp = 6.0*pii*gamma_sp( (xnu(lh) + 1.)/xmu(lh) )/gamma_sp( (xnu(lh) + 2.)/xmu(lh) )
          endif
        endif
        if ( lhl .gt. 1 ) then
          if ( dmuhl == 1.0 ) then
            cwchltmp = ((3. + dnu(lhl))*(2. + dnu(lhl))*(1.0 + dnu(lhl)))**(-1./3.)
          else
            cwchltmp = 6.0*pii*gamma_sp( (xnu(lhl) + 1)/xmu(lhl) )/gamma_sp( (xnu(lhl) + 2)/xmu(lhl) )
          endif
        endif

        if ( ipconc .le. 5 ) then
          if ( lh  .gt. 1 ) cwch(:) =  cwchtmp 
          if ( lhl .gt. 1 ) cwchl(:) = cwchltmp
        else
          do mgs = 1,ngscnt
          
          if ( lh  .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) then
           if ( qx(mgs,lh) .gt. qxmin(lh) ) then
            if ( dmuh == 1.0 ) then
              cwch(mgs) = ((3. + alpha(mgs,lh))*(2. + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))**(-1./3.)
             else
             xnutmp = (alpha(mgs,lh) - 2.0)/3.0
             cwch(mgs) =  6.0*pii*gamma_sp( (xnutmp + 1.)/xmu(lh) )/gamma_sp( (xnutmp + 2.)/xmu(lh) )
            endif
           else
             cwch(mgs) = cwchtmp
           endif
          endif
          if ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) then
           if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
            if ( dmuhl == 1.0 ) then
              cwchl(mgs) = ((3. + alpha(mgs,lhl))*(2. + alpha(mgs,lhl))*(1.0 + alpha(mgs,lhl)))**(-1./3.)
             else
             xnutmp = (alpha(mgs,lhl) - 2.0)/3.0
             cwchl(mgs) = 6.0*pii*gamma_sp( (xnutmp + 1)/xmu(lhl) )/gamma_sp( (xnutmp + 2)/xmu(lhl) )
            endif
           else
             cwchl(mgs) = cwchltmp
           endif
          endif
          
          enddo
        
        endif
       

      cimasn = Min( cimas0, 6.88e-13)
      cimasx = 1.0e-8
      ccimx = 5000.0e3   ! max of 5000 per liter

      cwc1 = 6.0/(pi*1000.)
      cwc0 = pii ! 6.0*pii
      mwfac = 6.0**(1./3.)

      
      if (ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set scale diameter'
!


!
!  cloud water variables
! ################################################################
!
!  DROPLETS
!
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set cloud water variables'
      
      if ( ildo == 0 .or. ildo == lc ) then
      
      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then !{
      
      if ( ipconc .ge. 2 ) then
        if ( cx(mgs,lc) .gt. cxmin) then !{
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        else
         cx(mgs,lc) = Max( cxmin, rho0(mgs)*qx(mgs,lc)/cwmasx )
         xmas(mgs,lc) = Min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
         xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        
        endif
      else
       if ( ipconc .lt. 2 ) then
         cx(mgs,lc) = rho0(mgs)*ccn/rho00 ! scales to local density, relative to standard air density
       endif
       if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 0.01 ) then !{
        xmas(mgs,lc) =  &
     &     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)), &
     &      xdn(mgs,lc)*xvmx(lc) )
        
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)
        
       elseif ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 1.0e-9 ) then
        cx(mgs,lc) = Max( cxmin, rho0(mgs)*qx(mgs,lc)/cwmasx )
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)

       elseif ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. 0.01 ) then
        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
        
       else
        xmas(mgs,lc) = cwmasn
        xv(mgs,lc) = xmas(mgs,lc)/1000.
! do not define ccw here! it can feed back to ccn!!!    cx(mgs,lc) = 0.0 ! cwnc(mgs)
       endif !}
      endif !}
!      IF ( ipconc .lt. 2 ) THEN
!        xmas(mgs,lc) = &
!     &    min( max(qx(mgs,lc)*rho0(mgs)/cwnc(mgs),cwmasn),cwmasx )
!        cx(mgs,lc) = Max(1.0,qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc))
!      ELSE
!        cwnc(mgs) = an(igs(mgs),jgs,kgs(mgs),lnc)
!        cx(mgs,lc) = cwnc(mgs)
!      ENDIF
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**(1./3.)
      xdia(mgs,lc,2) = xdia(mgs,lc,1)**2
      xdia(mgs,lc,3) = xdia(mgs,lc,1)
      cwrad = 0.5*xdia(mgs,lc,1)
      if ( fadvisc(mgs) > 0.0 ) then
      vtxbar(mgs,lc,1) =  &
     &   (2.0*gr*xdn(mgs,lc) *(cwrad**2)) &
     &  /(9.0*fadvisc(mgs))
      else
       vtxbar(mgs,lc,1) = 0.0
      endif

      
      else
       xmas(mgs,lc) = cwmasn
       xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
       if ( qx(mgs,lc) <= 0.0 ) cx(mgs,lc) = 0.0
       if ( ipconc .le. 1 ) cx(mgs,lc) = 0.01
       xdia(mgs,lc,1) = 2.*cwradn
       xdia(mgs,lc,2) = 4.*cwradn**2
       xdia(mgs,lc,3) = xdia(mgs,lc,1)
       vtxbar(mgs,lc,1) = 0.0
       
      endif !} qcw .gt. qxmin(lc)
      
      end do
      
      endif



!
! cloud ice variables
! columns
!
! ################################################################
!
!  CLOUD ICE
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set cip'
      
      if ( li .gt. 1 .and. ( ildo == 0 .or. ildo == li ) ) then
      do mgs = 1,ngscnt
       xdn(mgs,li)  = 900.0
      if ( ipconc .eq. 0 ) then
!       cx(mgs,li) = min(cnit*exp(-temcg(mgs)*bta1),1.e+09)
        cx(mgs,li) = cnina(mgs)
       if ( cimna .gt. 1.0 ) then
         cx(mgs,li) = Max(cimna,cx(mgs,li))
       endif
       if ( cimxa .gt. 1.0 ) then
         cx(mgs,li) = Min(cimxa,cx(mgs,li))
       endif
! erm 3/28/2002
       if ( itype1a .ge. 1 .or. itype2a .ge. 1 ) then
        cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
        cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)
       endif
!
       cx(mgs,li) = max(1.0e-20,cx(mgs,li))
!       cx(mgs,li) = Min(ccimx, cx(mgs,li))

      
      elseif ( ipconc .ge. 1 ) then
        if ( qx(mgs,li) .gt. qxmin(li) ) then
         cx(mgs,li) = Max(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasx)
         cx(mgs,li) = Min(cx(mgs,li),qx(mgs,li)*rho0(mgs)/cimasn)
!         cx(mgs,li) = Max(1.0,cx(mgs,li))
        endif
      endif
      
      if ( qx(mgs,li) .gt. qxmin(li) ) then
      xmas(mgs,li) = &
     &     max( qx(mgs,li)*rho0(mgs)/cx(mgs,li), cimasn )
!     &  min( max(qx(mgs,li)*rho0(mgs)/cx(mgs,li),cimasn),cimasx )
      
!      if ( temcg(mgs) .gt. 0.0 ) then
!      xdia(mgs,li,1) = 0.0
!      else
      if ( xmas(mgs,li) .gt. 0.0 ) then ! cimasn ) then
!c      xdia(mgs,li,1) = 0.4892*(xmas(mgs,li)**(0.4554))
!       xdia(mgs,li,1) = 0.1871*(xmas(mgs,li)**(0.3429))

!       xdia(mgs,li,1) = (132.694*5.40662/xmas(mgs,li))**(-1./2.9163)  ! for inverse exponential distribution
       if ( ixtaltype == 1 ) then ! column
       xdia(mgs,li,1) = 0.1871*(xmas(mgs,li)**(0.3429))
       xdia(mgs,li,3) = 0.1871*(xmas(mgs,li)**(0.3429))
       elseif  ( ixtaltype == 2 ) then ! disk
        xdia(mgs,li,1) = 0.277823*xmas(mgs,li)**0.359971
        xdia(mgs,li,3) = 0.277823*xmas(mgs,li)**0.359971
       endif
      end if
!      end if
!      xdia(mgs,li,1) = max(xdia(mgs,li,1), 5.e-6)
!      xdia(mgs,li,1) = min(xdia(mgs,li,1), 1000.e-6)

       if ( ipconc .ge. 0 ) then
!      vtxbar(mgs,li,1) = rhovt(mgs)*49420.*40.0005/5.40662*xdia(mgs,li,1)**(1.415) ! mass-weighted
!      vtxbar(mgs,li,1) = (4.942e4)*(xdia(mgs,li,1)**(1.4150))
        xv(mgs,li) = xmas(mgs,li)/xdn(mgs,li)
        if ( icefallopt == 1 ) then ! default ice fall
          if ( ixtaltype == 1 ) then ! column
          tmp = (67056.6300748612*rhovt(mgs))/  &
     &     (((1.0 + cinu)/xv(mgs,li))**0.4716666666666667*gfcinu1)
          vtxbar(mgs,li,2) = tmp*gfcinu1p47
          vtxbar(mgs,li,1) = tmp*gfcinu2p47/(1. + cinu)
          vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 
        elseif  ( ixtaltype == 2 ) then ! disk -- but just use Ferrier (1994) snow fall speeds for now
            vtxbar(mgs,li,1) = 11.9495*rhovt(mgs)*(xv(mgs,li))**(0.14)
            vtxbar(mgs,li,2) = 7.02909*rhovt(mgs)*(xv(mgs,li))**(0.14)
           vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 
        
          endif
          
       elseif ( icefallopt == 2 ) then !   ! Ferrier ice fall speed
          tmp = (82.3166*rhovt(mgs))/  &
     &     (((1.0 + cinu)/xv(mgs,li))**0.22117*gfcinu1)
          vtxbar(mgs,li,2) = tmp*gfcinu1p22
          vtxbar(mgs,li,1) = tmp*gfcinu2p22/(1. + cinu)
          vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 

       elseif ( icefallopt == 3 ) then !   ! Adjusted Ferrier (smaller exponent of 0.55 instead of 0.6635)
       
          tmp = (47.6273*rhovt(mgs))/  &
     &     (((1.0 + cinu)/xv(mgs,li))**0.18333*gfcinu1)
          vtxbar(mgs,li,2) = tmp*gfcinu1p18
          vtxbar(mgs,li,1) = tmp*gfcinu2p18/(1. + cinu)
          vtxbar(mgs,li,3) = vtxbar(mgs,li,1) 
       
       endif
!      vtxbar(mgs,li,1) = vtxbar(mgs,li,2)*(1.+cinu)/(1. + cinu)
!      xdn(mgs,li)   = min(max(769.8*xdia(mgs,li,1)**(-0.0140),300.0),900.0)
!      xdn(mgs,li) = 900.0
        xdia(mgs,li,2) = xdia(mgs,li,1)**2
!      vtxbar(mgs,li,1) = vtxbar(mgs,li,1)*rhovt(mgs)
       else
         xdia(mgs,li,1) = max(xdia(mgs,li,1), 10.e-6)
         xdia(mgs,li,1) = min(xdia(mgs,li,1), 1000.e-6)
         vtxbar(mgs,li,1) = (4.942e4)*(xdia(mgs,li,1)**(1.4150))
!      xdn(mgs,li)   = min(max(769.8*xdia(mgs,li,1)**(-0.0140),300.0),900.0)
         xdn(mgs,li) = 900.0
         xdia(mgs,li,2) = xdia(mgs,li,1)**2
         vtxbar(mgs,li,1) = vtxbar(mgs,li,1)*rhovt(mgs)
         xv(mgs,li) = xmas(mgs,li)/xdn(mgs,li)
       endif ! ipconc gt 3
      else
       xmas(mgs,li) = 1.e-13
       if ( qx(mgs,li) <= 0.0 ) cx(mgs,li) = 0.0
       xdn(mgs,li)  = 900.0
       xdia(mgs,li,1) = 1.e-7
       xdia(mgs,li,2) = (1.e-14)
       xdia(mgs,li,3) = 1.e-7
       vtxbar(mgs,li,1) = 0.0
!       cicap(mgs) = 0.0
!       ciat(mgs) = 0.0
      endif
      
      if ( icefallfac /= 1.0 ) then
        vtxbar(mgs,li,1) = icefallfac*vtxbar(mgs,li,1)
        vtxbar(mgs,li,2) = icefallfac*vtxbar(mgs,li,2)
        vtxbar(mgs,li,3) = icefallfac*vtxbar(mgs,li,3)
      endif

      
      
      end do
      
      endif ! li .gt. 1


! ################################################################
!
!  RAIN
!
      
!
      if ( ildo == 0 .or. ildo == lr ) then
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      
!      IF ( qx(mgs,lr) .gt. 10.0e-3 ) &
!     &  write(0,*)  'RAIN1: ',igs(mgs),kgs(mgs),qx(mgs,lr)
      
      if ( ipconc .ge. 3 ) then
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        xvbarmax = xvmx(lr)
        if ( imaxdiaopt == 1 ) then
          xvbarmax = xvmx(lr)
        elseif ( imaxdiaopt == 2 ) then ! test against maximum mass diameter
         if ( imurain == 1 ) then
           xvbarmax = xvmx(lr)/((3. + alpha(mgs,lr))**3/((3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr))))
         elseif ( imurain == 3 ) then
           
         endif
        elseif ( imaxdiaopt == 3 ) then ! test against mass-weighted diameter
         if ( imurain == 1 ) then
           xvbarmax = xvmx(lr)/((4. + alpha(mgs,lr))**3/((3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr))))
         elseif ( imurain == 3 ) then
           
         endif
        endif
       
        if ( xv(mgs,lr) .gt. xvbarmax ) then
          xv(mgs,lr) = xvbarmax
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvbarmax*xdn(mgs,lr))
        elseif ( xv(mgs,lr) .lt. xvmn(lr) ) then
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        endif


        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.) ! xdia(mgs,lr,1)
        if ( imurain == 3 ) then
!          xdia(mgs,lr,1) = (6.*pii*xv(mgs,lr)/(alpha(mgs,lr)+1.))**(1./3.)
          xdia(mgs,lr,1) = xdia(mgs,lr,3) ! formulae for Ziegler (1985) use mean volume diameter, not lambda**(-1)
        else ! imurain == 1, Characteristic diameter (1/lambda)
          xdia(mgs,lr,1) = (6.*pii*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)
        endif
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)

! Inverse exponential version:
!        xdia(mgs,lr,1) =
!     &  (qx(mgs,lr)*rho0(mgs)
!     & /(pi*xdn(mgs,lr)*cx(mgs,lr)))**(0.333333)
      else
        xdia(mgs,lr,1) = &
     &  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25) 
        xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.)
        cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
      end if
      else
        xdia(mgs,lr,1) = 1.e-9
        xdia(mgs,lr,3) = 1.e-9
        xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)
      end if
      xdia(mgs,lr,2) = xdia(mgs,lr,1)**2
!      xmas(mgs,lr) = xdn(mgs,lr)*(pi/6.)*xdia(mgs,lr,1)**3
      end do
      
      endif
! ################################################################
!
!  SNOW
!

      if ( ls .gt. 1 .and. ( ildo == 0 .or. ildo == ls ) ) then
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      if ( ipconc .ge. 4 ) then ! 

        xmas(mgs,ls) =  rho0(mgs)*qx(mgs,ls)/(Max(1.0e-9,cx(mgs,ls)))
        swmasmx = 13.7e-6
!       IF ( xmas(mgs,ls) > swmasmx ) THEN
!          xmas(mgs,ls) = swmasmx
!          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
!        ENDIF

        if ( isnowdens == 2 ) then ! Set values according to Cox relationship
        
          xdn(mgs,ls) = 0.0346159*Sqrt(cx(mgs,ls)/(qx(mgs,ls)*rho0(mgs)) )
          xdn(mgs,ls) = Max( 100.0, xdn(mgs,ls) )  ! limit snow to 100. to keep other equations in line
          
          if ( xdn(mgs,ls) <= 900. ) then
             dtmp = Sqrt( xmas(mgs,ls)/0.069 ) ! diameter (meters) of mean mass particle using Cox 1998 relation (m = p d^2)
             xv(mgs,ls) = 28.8887*xmas(mgs,ls)**(3./2.)
          else ! at small sizes, assume ice spheres
             xdn(mgs,ls) = 900.
             xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
             dtmp = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
          endif
          
        else ! leave xdn(ls) at default value
             xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*Max(1.0e-9,cx(mgs,ls)))
             dtmp = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
        endif

        xdia(mgs,ls,1) = dtmp ! (xv(mgs,ls)*cwc0*6.0)**(1./3.)

        if ( xv(mgs,ls) .lt. xvmn(ls) .and. isnowdens == 1) then
          xv(mgs,ls) = Max( xvmn(ls),xv(mgs,ls) )
          xmas(mgs,ls) = xv(mgs,ls)*xdn(mgs,ls)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
          xdia(mgs,ls,1) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
        endif

        if ( xv(mgs,ls) .gt. xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls))) ) then
          xv(mgs,ls) = Min( xvmx(ls), Max( xvmn(ls),xv(mgs,ls) ) )
          xmas(mgs,ls) = 0.106214*xv(mgs,ls)**(2./3.)
          cx(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xmas(mgs,ls))
          xdn(mgs,ls) = 0.0346159*Sqrt(cx(mgs,ls)/(qx(mgs,ls)*rho0(mgs)) )
          xdia(mgs,ls,1) = Sqrt( xmas(mgs,ls)/0.069 ) 
        endif

        xdia(mgs,ls,3) = xdia(mgs,ls,1)

      else
        xdia(mgs,ls,1) =  &
     &    (qx(mgs,ls)*rho0(mgs)/(pi*xdn(mgs,ls)*cnostmp(mgs)))**(0.25) 
        cx(mgs,ls) = cnostmp(mgs)*xdia(mgs,ls,1)
        xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*cx(mgs,ls))
        xdia(mgs,ls,3) = (xv(mgs,ls)*cwc0*6.0)**(1./3.)
      end if
      else
      xdia(mgs,ls,1) = 1.e-9
      xdia(mgs,ls,3) = 1.e-9
      cx(mgs,ls) = 0.0
      
       if ( isnowdens == 2 ) then ! Set values according to Cox relationship
         xdn(mgs,ls) = 90.
       endif

      end if
      xdia(mgs,ls,2) = xdia(mgs,ls,1)**2
!      swdia3(mgs) = xdia(mgs,ls,2)*xdia(mgs,ls,1)
!      xmas(mgs,ls) = xdn(mgs,ls)*(pi/6.)*swdia3(mgs)
      end do
      
      endif ! ls .gt 1
!
!
! ################################################################
!
!  GRAUPEL
!

      if ( lh .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) then
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*Max(1.0e-9,cx(mgs,lh)))
        xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)

        if ( xv(mgs,lh) .lt. xvmn(lh) .or. xv(mgs,lh) .gt. xvmx(lh) ) then
          xv(mgs,lh) = Min( xvmx(lh), Max( xvmn(lh),xv(mgs,lh) ) )
          xmas(mgs,lh) = xv(mgs,lh)*xdn(mgs,lh)
          cx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xmas(mgs,lh))
        endif

         xdia(mgs,lh,3) = (xv(mgs,lh)*6.*pii)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
         if ( dmuh == 1.0 ) then
           xdia(mgs,lh,1) = cwch(mgs)*xdia(mgs,lh,3)
         else
           xdia(mgs,lh,1) = (xv(mgs,lh)*cwch(mgs))**(1./3.)
         endif

      else
      xdia(mgs,lh,1) =  &
     &  (qx(mgs,lh)*rho0(mgs)/(pi*xdn(mgs,lh)*cno(lh)))**(0.25) 
      cx(mgs,lh) = cno(lh)*xdia(mgs,lh,1)
      xv(mgs,lh) = Max(xvmn(lh), rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh)) )
      xdia(mgs,lh,3) = (xv(mgs,lh)*6./pi)**(1./3.) 
      end if
      else
      xdia(mgs,lh,1) = 1.e-9
      xdia(mgs,lh,3) = 1.e-9
      end if
      xdia(mgs,lh,2) = xdia(mgs,lh,1)**2
!      hwdia3(mgs) = xdia(mgs,lh,2)*xdia(mgs,lh,1)
!      xmas(mgs,lh) = xdn(mgs,lh)*(pi/6.)*hwdia3(mgs)
      end do
      
      endif

!
! ################################################################
!
!  HAIL
!

      if ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) then
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      if ( ipconc .ge. 5 ) then

        xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*Max(1.0e-9,cx(mgs,lhl)))
        xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
!        write(0,*) 'setvt: xv = ',xv(mgs,lhl),xdn(mgs,lhl),cx(mgs,lhl),xmas(mgs,lhl),qx(mgs,lhl)

        if ( xv(mgs,lhl) .lt. xvmn(lhl) .or. xv(mgs,lhl) .gt. xvmx(lhl) ) then
          xv(mgs,lhl) = Min( xvmx(lhl), Max( xvmn(lhl),xv(mgs,lhl) ) )
          xmas(mgs,lhl) = xv(mgs,lhl)*xdn(mgs,lhl)
          cx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xmas(mgs,lhl))
        endif

        xdia(mgs,lhl,3) = (xv(mgs,lhl)*6./pi)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
         if ( dmuhl == 1.0 ) then
           xdia(mgs,lhl,1) = cwchl(mgs)*xdia(mgs,lhl,3)
         else
           xdia(mgs,lhl,1) = (xv(mgs,lhl)*cwchl(mgs))**(1./3.)
         endif
        
!        write(0,*) 'setvt: xv = ',xv(mgs,lhl),xdn(mgs,lhl),cx(mgs,lhl),xdia(mgs,lhl,3)
      else
      xdia(mgs,lhl,1) = &
     &  (qx(mgs,lhl)*rho0(mgs)/(pi*xdn(mgs,lhl)*cno(lhl)))**(0.25) 
      cx(mgs,lhl) = cno(lhl)*xdia(mgs,lhl,1)
      xv(mgs,lhl) = Max(xvmn(lhl), rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl)) )
      xdia(mgs,lhl,3) = (xv(mgs,lhl)*6./pi)**(1./3.) 
      end if
      else
      xdia(mgs,lhl,1) = 1.e-9
      xdia(mgs,lhl,3) = 1.e-9
      end if
      xdia(mgs,lhl,2) = xdia(mgs,lhl,1)**2
!      hwdia3(mgs) = xdia(mgs,lh,2)*xdia(mgs,lh,1)
!      xmas(mgs,lh) = xdn(mgs,lh)*(pi/6.)*hwdia3(mgs)
      end do
      
      endif
!      
!
!
!  Set terminal velocities...
!    also set drag coefficients (moved to start of subroutine)
!
!      cdx(lr) = 0.60
!      cdx(lh) = 0.45
!      cdx(lhl) = 0.45
!      cdx(lf) = 0.45
!      cdx(lgh) = 0.60
!      cdx(lgm) = 0.80
!      cdx(lgl) = 0.80
!      cdx(lir) = 2.00
!
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set terminal velocities'
!
!
! ################################################################
!
!  RAIN
!
      if ( ildo == 0 .or. ildo == lr ) then
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
      if ( ipconc .lt. 3 ) then
        vtxbar(mgs,lr,1) = rainfallfac*(ar*gf4br/6.0)*(xdia(mgs,lr,1)**br)*rhovt(mgs)
!        write(91,*) 'vtxbar: ',vtxbar(mgs,lr,1),mgs,gf4br,xdia(mgs,lr,1),rhovt(mgs)
      else
        
        if ( imurain == 1 ) then ! DSD of Diameter
        
        ! using functional form of  arx*(1 - Exp(-frx*diameter) ), with arx =       arx = 10.
        !  and frx = 516.575 ! raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.
        ! Similar form as in Atlas et al. (1973), who had 9.65 - 10.3*Exp[-600 * d]

        
          alp = alpha(mgs,lr)
          
          vtxbar(mgs,lr,1) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 4.0) ) ! mass weighted
          
          if ( infdo .ge. 1 .and. rssflg == 1 ) then
            vtxbar(mgs,lr,2) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 1.0) ) ! number weighted
          else
            vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
          endif
          
          if ( infdo .ge. 2 .and. rssflg == 1 ) then
            vtxbar(mgs,lr,3) = rhovt(mgs)*arx*(1.0 - (1.0 + frx*xdia(mgs,lr,1))**(-alp - 7.0) ) ! z-weighted
          else
            vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
          endif
          
!          write(91,*) 'setvt: alp,vn,vm,vz = ',alp,vtxbar(mgs,lr,2), vtxbar(mgs,lr,1), vtxbar(mgs,lr,3),alpha(mgs,lr)

        elseif ( imurain == 3 ) then ! DSD of Volume
        
        if ( lzr < 1 ) then ! not 3-moment rain
        rwdia = Min( xdia(mgs,lr,1), 8.0e-3 )
        
         vtxbar(mgs,lr,1) = rhovt(mgs)*6.0*pii*( 0.04771 + 3788.0*rwdia -  &
     &        1.105e6*rwdia**2 + 1.412e8*rwdia**3 - 6.527e9*rwdia**4)
        
        if ( infdo .ge. 1 ) then
          if (  rssflg >= 1 ) then
         vtxbar(mgs,lr,2) = (0.09112 + 2714.0*rwdia - 4.872e5*rwdia**2 +  &
     &            4.495e7*rwdia**3 - 1.626e9*rwdia**4)*rhovt(mgs)
          else
            vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
          endif
        endif
        
        if ( infdo .ge. 2 ) then ! Z-weighted fall speed
        vtxbar(mgs,lr,3)  = rhovt(mgs)*(  &
     &       0.0911229 +                  &
     &  9246.494*(rwdia) -               &
     &  3.2839926e6*(rwdia**2) +          &
     &  4.944093e8*(rwdia**3) -          &
     &  2.631718e10*(rwdia**4) )
        endif
        
        else ! 3-moment rain, gamma-volume

        vr = xv(mgs,lr)
        rnux = alpha(mgs,lr)
        
        if ( infdo .ge. 1 .and. rssflg == 1) then ! number-weighted; DTD: added size-sorting flag
        vtxbar(mgs,lr,2) = rhovt(mgs)*                             &
     &     (((1. + rnux)/vr)**(-1.333333)*                         &
     &    (0.0911229*((1. + rnux)/vr)**1.333333*Gamma_sp(1. + rnux) + &
     &      (5430.3131*(1. + rnux)*Gamma_sp(4./3. + rnux))/           &
     &       vr - 1.0732802e6*((1. + rnux)/vr)**0.6666667*         &
     &       Gamma_sp(1.666667 + rnux) +                              &
     &      8.584110982429507e7*((1. + rnux)/vr)**(1./3.)*         &
     &       Gamma_sp(2. + rnux) -                                    &
     &      2.3303765697228556e9*Gamma_sp(7./3. + rnux)))/            &
     &  Gamma_sp(1. + rnux)
        endif

!  mass-weighted
       vtxbar(mgs,lr,1)  = rhovt(mgs)*                                                 &
     &   (0.0911229*(1 + rnux)**1.3333333333333333*Gamma_sp(2. + rnux) +                  &
     &    5430.313059683277*(1 + rnux)*vr**0.3333333333333333*                         &
     &     Gamma_sp(2.333333333333333 + rnux) -                                           &
     &    1.0732802065650471e6*(1 + rnux)**0.6666666666666666*vr**0.6666666666666666*  &
     &     Gamma_sp(2.6666666666666667 + rnux) +                                          &
     &    8.584110982429507e7*(1 + rnux)**0.3333333333333333*vr*Gamma_sp(3 + rnux) -      &
     &    2.3303765697228556e9*vr**1.3333333333333333*                                 &
     &     Gamma_sp(3.333333333333333 + rnux))/                                           &
     &  ((1 + rnux)**2.333333333333333*Gamma_sp(1 + rnux)) 
     
        if(infdo .ge. 1 .and. rssflg == 0) then ! No size-sorting, set N-weighted fall speed to mass-weighted
          vtxbar(mgs,lr,2) = vtxbar(mgs,lr,1)
        endif     
      
        if ( infdo .ge. 2 .and. rssflg == 1) then ! Z-weighted fall speed
        vtxbar(mgs,lr,3)  =   rhovt(mgs)*                                          &
     &  ((1. + rnux)*(0.0911229*(1 + rnux)**1.3333333333333333*Gamma_sp(3. + rnux) +  &
     &      5430.313059683277*(1 + rnux)*vr**0.3333333333333333*                   &
     &       Gamma_sp(3.3333333333333335 + rnux) -                                    &
     &      1.0732802065650471e6*(1 + rnux)**0.6666666666666666*                   &
     &       vr**0.6666666666666666*Gamma_sp(3.6666666666666665 + rnux) +             &
     &      8.5841109824295e7*(1 + rnux)**0.3333333333333333*vr*Gamma_sp(4. + rnux) - &
     &      2.3303765697228556e9*vr**1.3333333333333333*                           &
     &       Gamma_sp(4.333333333333333 + rnux)))/                                    &
     &  ((1 + rnux)**3.3333333333333335*(2 + rnux)*Gamma_sp(1 + rnux))
        
!         write(0,*) 'setvt: mgs,lzr,infdo = ',mgs,lzr,infdo
!         write(0,*) 'vt1,2,3 = ',vtxbar(mgs,lr,1),vtxbar(mgs,lr,2),vtxbar(mgs,lr,3)
        
        elseif (infdo .ge. 2) then ! No size-sorting, set Z-weighted fall speed to mass-weighted
          vtxbar(mgs,lr,3) = vtxbar(mgs,lr,1)
        endif
        
        
        endif
       endif ! imurain

!        IF ( rwrad*mwfac .gt. 6.0e-4  ) THEN
!          vtxbar(mgs,lr,1) = 20.1*Sqrt(100.*rwrad*mwfac)*rhovt(mgs)
!        ELSE
!          vtxbar(mgs,lr,1) = 80.0e2*rwrad*rhovt(mgs)*mwfac
!        ENDIF
!        IF ( rwrad .gt. 6.0e-4  ) THEN
!          vtxbar(mgs,lr,2) = 20.1*Sqrt(100.*rwrad)*rhovt(mgs)
!        ELSE
!          vtxbar(mgs,lr,2) = 80.0e2*rwrad*rhovt(mgs)
!        ENDIF
      endif ! ipconc
      else  ! qr < qrmin
      vtxbar(mgs,lr,1) = 0.0
      vtxbar(mgs,lr,2) = 0.0
      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set rain vt'
      
      endif
!
! ################################################################
!
!  SNOW !Zrnic et al. (1993)
!
      if ( ls .gt. 1 .and. ( ildo == 0 .or. ildo == ls ) ) then
      do mgs = 1,ngscnt
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
        if ( ipconc .ge. 4 ) then
         if ( mixedphase .and. qsvtmod ) then
         else
          if ( isnowfall == 1 ) then
           ! original (Zrnic et al. 1993)
           vtxbar(mgs,ls,1) = 5.72462*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
          elseif ( isnowfall == 2 ) then
          ! Ferrier:
            if ( isnowdens == 1 ) then
              vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls))**(0.14)
            else
              vtxbar(mgs,ls,1) = 11.9495*rhovt(mgs)*(xv(mgs,ls)*xdn(mgs,ls)/100.)**(0.14) 
            endif
          elseif ( isnowfall == 3 ) then
          ! Cox, mass distrib:
            vtxbar(mgs,ls,1) = 50.092*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
          endif
          
          if(Abs(sssflg) >= 1) then
            if ( isnowfall == 1 ) then
              vtxbar(mgs,ls,2) = 4.04091*rhovt(mgs)*(xv(mgs,ls))**(1./12.)
            elseif ( isnowfall == 2 ) then
            ! Ferrier:
              if ( isnowdens == 1 ) then
                vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls))**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
              else
                vtxbar(mgs,ls,2) = 7.02909*rhovt(mgs)*(xv(mgs,ls)*xdn(mgs,ls)/100.)**(0.14)  ! bug fix 11/15/2015: was rewriting to mass fall speed vtxbar(mgs,ls,1)
              endif
            elseif ( isnowfall == 3 ) then
            ! Cox, mass distrib:
              vtxbar(mgs,ls,2) = 21.6147*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
            endif
          else
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
          endif
           if ( infdo  >= 2 ) then
            if ( isnowfall == 1 ) then
             vtxbar(mgs,ls,3) = 6.12217*rhovt(mgs)*(xv(mgs,ls))**(1./12.) ! Zrnic et al 93
            elseif ( isnowfall == 2 ) then
             vtxbar(mgs,ls,3) = 13.3436*rhovt(mgs)*(xv(mgs,ls))**(0.14)   ! Ferrier 94
            elseif ( isnowfall == 3 ) then
            ! Cox, mass distrib:
              vtxbar(mgs,ls,3) = 61.0914*rhovt(mgs)*(xmas(mgs,ls))**(0.2635)
            endif
           endif
         
         if ( sssflg < 0 .and. temcg(mgs) > Abs(sssflg) ) then ! above a given temperature, effectively turn off size sorting
            vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
            vtxbar(mgs,ls,3) = vtxbar(mgs,ls,1)
         endif
         
         endif
        else ! single-moment:
         vtxbar(mgs,ls,1) = (cs*gf4ds/6.0)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
         vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        endif
      else
      vtxbar(mgs,ls,1) = 0.0
      end if

      if ( snowfallfac /= 1.0 ) then
        vtxbar(mgs,ls,1) = snowfallfac*vtxbar(mgs,ls,1)
        vtxbar(mgs,ls,2) = snowfallfac*vtxbar(mgs,ls,2)
        vtxbar(mgs,ls,3) = snowfallfac*vtxbar(mgs,ls,3)
      endif


      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set snow vt'
      
      endif ! ls .gt. 1
!
!
! ################################################################
!
!  GRAUPEL !Wisner et al. (1972)
!
      if ( lh .gt. 1 .and. ( ildo == 0 .or. ildo == lh ) ) then
      
      do mgs = 1,ngscnt
      vtxbar(mgs,lh,1) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
         cd = cdx(lh)
        if ( icdx .eq. 1 ) then
         cd = cdx(lh)
       elseif ( icdx .eq. 2 ) then
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lh) - xdn(mgs,lh))/(xdnmx(lh)-xdnmn(lh)) ) )
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lh))/(900. - 300.) ) )
         cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
!         cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
       elseif ( icdx .eq. 3 ) then
!         cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 300.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
       elseif ( icdx .eq. 4 ) then
         cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
     &        (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lh) ) ) )/(cdhdnmax - cdhdnmin) ) )
       elseif ( icdx .eq. 5 ) then
         cd = cdx(lh)*(xdn(mgs,lh)/rho_qh)**(2./3.)
       elseif ( icdx .eq. 6 ) then ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lh)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )
         
         
         delrho = Max( 0.0, 0.01*(xdn(mgs,lh) - mmgraupvt(indxr,1)) )
         if ( indxr < ngdnmm ) then
          
          axx(mgs,lh) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxx(mgs,lh) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )

          
         else
          axx(mgs,lh) = mmgraupvt(indxr,2)
          bxx(mgs,lh) = mmgraupvt(indxr,3)
         endif
         
         aax = axx(mgs,lh)
         bbx = bxx(mgs,lh)

         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
         
       elseif ( icdx <= 0 ) then ! 
         aax = ax(lh)
         bbx = bx(lh)
          cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
       else
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
       endif
       
       cdxgs(mgs,lh) = cd
      if ( alpha(mgs,lh) .eq. 0.0 .and. icdx > 0 .and. icdx /= 6 ) then
!      axx(mgs,lh) =  (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lh)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      axx(mgs,lh) = Sqrt(4.0*xdn(mgs,lh)*gr/(3.0*cd*rho00))
      bxx(mgs,lh) = 0.5
      vtxbar(mgs,lh,1) = (gf4p5/6.0)* rhovt(mgs)*axx(mgs,lh) * Sqrt(xdia(mgs,lh,1)) 
!      vtxbar(mgs,lh,1) = (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lh)*xdia(mgs,lh,1)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      else
        if ( icdx /= 6 ) bbx = bx(lh)
        tmp = 4. + alpha(mgs,lh) + bbx
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
!        aax = Max( 1.0, Min(2.0, (xdn(mgs,lh)/400.) ) )
!        vtxbar(mgs,lh,1) =  rhovt(mgs)*aax*ax(lh)*(xdia(mgs,lh,1)**bx(lh)*x)/y
        
        if ( icdx > 0 .and. icdx /= 6) then
          aax = Sqrt(4.0*xdn(mgs,lh)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lh,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lh,1)) * x/y
          axx(mgs,lh) = aax
          bxx(mgs,lh) = bbx
        elseif (icdx == 6 ) then
          vtxbar(mgs,lh,1) =  rhovt(mgs)*aax* xdia(mgs,lh,1)**bbx * x/y
        else ! icdx < 0
          axx(mgs,lh) = ax(lh)
          bxx(mgs,lh) = bx(lh)
          vtxbar(mgs,lh,1) =  rhovt(mgs)*ax(lh)*(xdia(mgs,lh,1)**bx(lh)*x)/y          
        endif

!     &    Gamma_sp(4.0 + dnu(lh) + 0.6))/Gamma_sp(4. + dnu(lh))
      endif

      if ( lwsm6 .and. ipconc == 0 ) then
!         vtxbar(mgs,lh,1) = (330.*gf4ds/6.0)*(xdia(mgs,ls,1)**ds)*rhovt(mgs)
         vtxbar(mgs,lh,1) = (330.*gf4br/6.0)*(xdia(mgs,lh,1)**br)*rhovt(mgs)
      endif
      
      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'
      
      endif ! lh .gt. 1
!
!
! ################################################################
!
!  HAIL
!
      if ( lhl .gt. 1 .and. ( ildo == 0 .or. ildo == lhl ) ) then
      
      do mgs = 1,ngscnt
      vtxbar(mgs,lhl,1) = 0.0
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then

       if ( icdxhl .eq. 1 ) then
         cd = cdx(lhl)
       elseif ( icdxhl .eq. 3 ) then
!         cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 300.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
       elseif ( icdxhl .eq. 4 ) then
         cd = Max(cdhlmin, Min(cdhlmax, cdhlmin + (cdhlmax-cdhlmin)*  &
     &       (cdhldnmax - Max( cdhldnmin, Min( cdhldnmax, xdn(mgs,lhl) ) ) )/(cdhldnmax - cdhldnmin) ) )
       elseif ( icdxhl .eq. 5 ) then
         cd = cdx(lh)*(xdn(mgs,lhl)/rho_qh)**(2./3.)
       elseif ( icdxhl .eq. 6 ) then ! Milbrandt and Morrison (2013)
         indxr = Int( (xdn(mgs,lhl)-50.)/100. ) + 1
         indxr = Min( ngdnmm, Max(1,indxr) )
         
         
         delrho = Max( 0.0, 0.01*(xdn(mgs,lhl) - mmgraupvt(indxr,1)) )
         if ( indxr < ngdnmm ) then
          
          axx(mgs,lhl) = mmgraupvt(indxr,2) + delrho*(mmgraupvt(indxr+1,2) - mmgraupvt(indxr,2) )
          bxx(mgs,lhl) = mmgraupvt(indxr,3) + delrho*(mmgraupvt(indxr+1,3) - mmgraupvt(indxr,3) )

          
         else
          axx(mgs,lhl) = mmgraupvt(indxr,2)
          bxx(mgs,lhl) = mmgraupvt(indxr,3)
         endif
         
         aax = axx(mgs,lhl)
         bbx = bxx(mgs,lhl)

         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
         
       else
!         cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lhl))/(900. - 300.) ) )
!        cd = Max(0.5, Min(0.8, 0.5 + 0.3*(xdnmx(lhl) - xdn(mgs,lhl))/(xdnmx(lhl)-xdnmn(lhl)) ) )
!         cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 500.) ) )
         cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
       endif

       cdxgs(mgs,lhl) = cd

      if ( alpha(mgs,lhl) .eq. 0.0 .and. icdxhl > 0 .and. icdxhl /= 6) then
!      axx(mgs,lhl) =  (gf4p5/6.0)*  &
!     &  Sqrt( (xdn(mgs,lhl)*4.0*gr) /  &
!     &    (3.0*cd*rho0(mgs)) )
      axx(mgs,lhl) = Sqrt(4.0*xdn(mgs,lhl)*gr/(3.0*cd*rho00))
      bxx(mgs,lhl) = 0.5
      vtxbar(mgs,lhl,1) = (gf4p5/6.0)* rhovt(mgs)*axx(mgs,lhl) * Sqrt(xdia(mgs,lhl,1)) 
      else
        if ( icdxhl /= 6 ) bbx = bx(lhl)
        tmp = 4. + alpha(mgs,lhl) + bbx
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        if ( icdxhl > 0 .and. icdxhl /= 6) then
          aax = Sqrt(4.0*xdn(mgs,lhl)*gr/(3.0*cd*rho00))
          vtxbar(mgs,lhl,1) =  rhovt(mgs)*aax* Sqrt(xdia(mgs,lhl,1)) * x/y
          axx(mgs,lhl) = aax
          bxx(mgs,lhl) = bbx
        elseif ( icdxhl == 6 ) then
          vtxbar(mgs,lhl,1) =  rhovt(mgs)*aax* (xdia(mgs,lhl,1))**bbx * x/y
        else
          axx(mgs,lhl) = ax(lhl)
          bxx(mgs,lhl) = bx(lhl)
         vtxbar(mgs,lhl,1) =  rhovt(mgs)*(ax(lhl)*xdia(mgs,lhl,1)**bx(lhl)*x)/y
        endif
        
!     &    Gamma_sp(4.0 + dnu(lh) + 0.6))/Gamma_sp(4. + dnu(lh))
      endif


      end if
      end do
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt'
      
      endif ! lhl .gt. 1


      if ( infdo .ge. 1 ) then

!      DO il = lc,lhab
!      IF ( il .ne. lr ) THEN
        do mgs = 1,ngscnt
          if ( ildo == 0 .or. ildo == lc ) then
            vtxbar(mgs,lc,2) = vtxbar(mgs,lc,1)
          endif
        if ( li .gt. 1 ) then
!          vtxbar(mgs,li,2) = rhovt(mgs)*49420.*1.25447*xdia(mgs,li,1)**(1.415) ! n-wgt (Ferrier 94)
!          vtxbar(mgs,li,2) = vtxbar(mgs,li,1)

! test print stuff...
!          IF ( xdia(mgs,li,1) .gt. 200.e-6 ) THEN
!            tmp = (xv(mgs,li)*cwc0)**(1./3.)
!            x = rhovt(mgs)*49420.*40.0005/5.40662*tmp**(1.415)
!            y = rhovt(mgs)*49420.*1.25447*tmp**(1.415)
!            write(6,*) 'Ice fall: ',vtxbar(mgs,li,1),x,y,tmp,xdia(mgs,li,1)
!          ENDIF
        endif
!          vtxbar(mgs,ls,2) = vtxbar(mgs,ls,1)
        enddo

        if ( lg .gt. lr ) then

        do il = lg,lhab
         if ( ildo == 0 .or. ildo == il ) then

            do mgs = 1,ngscnt
             if ( qx(mgs,il) .gt. qxmin(il) ) then
              if ( (il .eq. lh .and. hssflg == 1) .or. ( lhl .gt. 1 .and. il .eq. lhl .and. hlssflg == 1) ) then ! DTD: added flag for size-sorting
              
              ! DTD: allow for setting of number-weighted and z-weighted fall speeds to the mass-weighted value,
              ! effectively turning off size-sorting

              if ( il .eq. lh ) then ! {
             
               if ( icdx .eq. 1 ) then
                 cd = cdx(lh)
               elseif ( icdx .eq. 2 ) then
!                 cd = Max(0.6, Min(1.0, 0.6 + 0.4*(xdnmx(lh) - xdn(mgs,lh))/(xdnmx(lh)-xdnmn(lh)) ) )
!                 cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lh))/(900. - 300.) ) )
                 cd = Max(0.45, Min(1.0, 0.45 + 0.35*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
!                 cd = Max(0.55, Min(1.0, 0.55 + 0.25*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 500.) ) )
               elseif ( icdx .eq. 3 ) then
!                 cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 170.0, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
                 cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hdnmn, Min( 800.0, xdn(mgs,lh) ) ) )/(800. - 170.0) ) )
               elseif ( icdx .eq. 4 ) then
                 cd = Max(cdhmin, Min(cdhmax, cdhmin + (cdhmax-cdhmin)* &
     &            (cdhdnmax - Max( cdhdnmin, Min( cdhdnmax, xdn(mgs,lh) ) ) )/(cdhdnmax - cdhdnmin) ) )
               elseif ( icdx .eq. 5 ) then
                 cd = cdx(lh)*(xdn(mgs,lh)/rho_qh)**(2./3.)
               elseif ( icdx .eq. 6 ) then ! Milbrandt and Morrison (2013)
                  aax = axx(mgs,lh)
                  bbx = bxx(mgs,lh)
               elseif ( icdx <= 0 ) then ! 
                  aax = ax(lh)
                  bbx = bx(lh)
               endif
               
              elseif ( lhl .gt. 1 .and. il .eq. lhl ) then
             
               if ( icdxhl .eq. 1 ) then
                 cd = cdx(lhl)
               elseif ( icdxhl .eq. 3 ) then
!               cd = Max(0.45, Min(1.0, 0.45 + 0.55*(800.0 - Max( 300., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 300.) ) )
                cd = Max(0.45, Min(1.2, 0.45 + 0.55*(800.0 - Max( hldnmn, Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 170.0) ) )
               elseif ( icdxhl .eq. 4 ) then
                cd = Max(cdhlmin, Min(cdhlmax, cdhlmin + (cdhlmax-cdhlmin)*  &
     &               (cdhldnmax - Max( cdhldnmin, Min( cdhldnmax, xdn(mgs,lhl) ) ) )/(cdhldnmax - cdhldnmin) ) )
               elseif ( icdxhl == 5 ) then
!                cd = Max(0.6, Min(1.0, 0.6 + 0.4*(900.0 - xdn(mgs,lhl))/(900. - 300.) ) )
!                cd = Max(0.5, Min(0.8, 0.5 + 0.3*(xdnmx(lhl) - xdn(mgs,lhl))/(xdnmx(lhl)-xdnmn(lhl)) ) )
                 cd = Max(0.45, Min(0.6, 0.45 + 0.15*(800.0 - Max( 500., Min( 800.0, xdn(mgs,lhl) ) ) )/(800. - 500.) ) )
               elseif ( icdxhl .eq. 6 ) then ! Milbrandt and Morrison (2013)
                  aax = axx(mgs,lhl)
                  bbx = bxx(mgs,lhl)
               elseif ( icdxhl <= 0 ) then ! 
                  aax = ax(lhl)
                  bbx = bx(lhl)
               endif
               
              endif ! }

               if ( alpha(mgs,il) .eq. 0. .and. infdo .lt. 2 .and.   &
               ( ( il==lh .and. icdx > 0 .and. icdx /= 6) .or. ( il==lhl .and. icdxhl > 0 .and. icdxhl /= 6 ) ) ) then ! {
                 vtxbar(mgs,il,2) =   &
     &              Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) / &
     &                (3.0*cd*Max(0.05,rho0(mgs))) )

               else
               if ( il == lh  .and. icdx   /= 6 ) bbx = bx(il)
               if ( il == lhl .and. icdxhl /= 6 ) bbx = bx(il)
               tmp = 1. + alpha(mgs,il) + bbx
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
  
               tmp = 1. + alpha(mgs,il)
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

                 if ( il .eq. lh  .or. il .eq. lhl) then ! {
                   if ( ( il==lh .and. icdx > 0 ) ) then
                     if ( icdx /= 6 ) then
                      aax = Sqrt(4.0*xdn(mgs,il)*gr/(3.0*cd*rho00))
                      vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bx(il) * x/y
                     else !  (icdx == 6 ) then
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bbx * x/y
                     endif

                   elseif ( ( il==lhl .and. icdxhl > 0 ) ) then
                     if ( icdxhl /= 6 ) then
                       aax = Sqrt(4.0*xdn(mgs,il)*gr/(3.0*cd*rho00))
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bx(il) * x/y
                     else ! ( icdxhl == 6 )
                       vtxbar(mgs,il,2) =  rhovt(mgs)*aax* xdia(mgs,il,1)**bbx * x/y
                     endif
                   else ! get here if il==lh and icdx < 0 -- or -- il==lhl and icdxhl < 0
                     aax = ax(il)
                     vtxbar(mgs,il,2) =  rhovt(mgs)*ax(il)*(xdia(mgs,il,1)**bx(il)*x)/y
                   endif
!                  vtxbar(mgs,il,2) =  &
!     &               rhovt(mgs)*(xdn(mgs,il)/400.)*(75.715*xdia(mgs,il,1)**0.6* &
!     &               x)/y
!                  vtxbar(mgs,il,2) =  &
!     &               rhovt(mgs)*(xdn(mgs,il)/400.)*(ax(il)*xdia(mgs,il,1)**bx(il)* &
!     &               x)/y
                  if ( infdo .ge. 2 ) then ! Z-weighted

               tmp = 7. + alpha(mgs,il) + bbx
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
  
               tmp = 7. + alpha(mgs,il)
               i = Int(dgami*(tmp))
               del = tmp - dgam*i
               y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

                   vtxbar(mgs,il,3) = rhovt(mgs)*                 &
     &                (aax*(xdia(mgs,il,1) )**bbx *  &
     &                 x)/y
!     &                 Gamma(7.0 + alpha(mgs,il) + bbx)/Gamma(7. + alpha(mgs,il))
          if ( .not. (vtxbar(mgs,il,1) > -1. .and. vtxbar(mgs,il,1) < 200. ) .or. &
               .not. (vtxbar(mgs,il,3) > -1. .and. vtxbar(mgs,il,3) < 200. ) ) then
           write(0,*) 'Setvtz: problem with vtxbar1/3: ',il,vtxbar(mgs,il,1),vtxbar(mgs,il,3),aax,bbx,x,y
           write(0,*) 'q, number, diam1,3(mm) = ', qx(mgs,il),cx(mgs,il),1000.*xdia(mgs,il,1),1000.*xdia(mgs,il,3)
           ! call commasmpi_abort()
          endif
!     &                (aax*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &                 Gamma_sp(7.0 + alpha(mgs,il) + bx(il)))/Gamma_sp(7. + alpha(mgs,il))
                  endif

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt3'

                 else ! hail
                  vtxbar(mgs,il,2) =  &
     &               rhovt(mgs)*(ax(il)*xdia(mgs,il,1)**bx(il)* &
     &               x)/y

                 if ( infdo .ge. 2 ) then ! Z-weighted
                  vtxbar(mgs,il,3) = rhovt(mgs)*                 &
     &              (aax*(1.0/xdia(mgs,il,1) )**(- bbx)*  &
     &               Gamma_sp(7.0 + alpha(mgs,il) + bbx))/Gamma_sp(7. + alpha(mgs,il))
!     &              (ax(il)*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &               Gamma_sp(7.0 + alpha(mgs,il) + bx(il)))/Gamma_sp(7. + alpha(mgs,il))
                  endif

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set hail vt4'

                 endif ! }
!     &             Gamma_sp(1.0 + dnu(il) + 0.6)/Gamma_sp(1. + dnu(il))
               endif ! }

!              IF ( infdo .ge. 2 ) THEN ! Z-weighted
!               vtxbar(mgs,il,3) = rhovt(mgs)*                 &
!     &            (ax(il)*(1.0/xdia(mgs,il,1) )**(- bx(il))*  &
!     &             Gamma_sp(7.0 + alpha(mgs,il) + bx(il)))/Gamma_sp(7. + alpha(mgs,il))
!              ENDIF

!               IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!                write(0,*) 'setvt: ',qx(mgs,il),xdia(mgs,il,1),xdia(mgs,il,3),dnu(il),ax(il),bx(il)
!               ENDIF
             elseif ( (il .eq. lh .and. hssflg == 0) .or. ( lhl .gt. 1 .and. il .eq. lhl .and. hlssflg == 0) ) then ! no size-sorting for graupel or hail
              vtxbar(mgs,il,2) = vtxbar(mgs,il,1)
              vtxbar(mgs,il,3) = vtxbar(mgs,il,1)
             else ! not lh or lhl
              vtxbar(mgs,il,2) = &
     &            Sqrt( (xdn(mgs,il)*xdia(mgs,il,1)*pi*gr) /  &
     &              (3.0*cdx(il)*Max(0.05,rho0(mgs))) )
              vtxbar(mgs,il,3) = vtxbar(mgs,il,1)

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt5'


              endif
             else ! qx < qxmin
              vtxbar(mgs,il,2) = 0.0

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt6'

             endif
           enddo ! mgs

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt7'

        endif
        enddo ! il

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt8'

        endif ! lg .gt. 1 
        
!      ENDIF
!      ENDDO

      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: Set graupel vt9'

!       DO mgs = 1,ngscnt
!        IF ( qx(mgs,lr) > qxmin(lr) ) THEN
!         write(0,*) 'setvt2: mgs,lzr,infdo = ',mgs,lzr,infdo
!         write(0,*) 'vt1,2,3 = ',vtxbar(mgs,lr,1),vtxbar(mgs,lr,2),vtxbar(mgs,lr,3)
!        ENDIF
!       ENDDO

      endif ! infdo .ge. 1 

        if (  lh > 0 .and. graupelfallfac /= 1.0 ) then
          do mgs = 1,ngscnt
            vtxbar(mgs,lh,1) = graupelfallfac*vtxbar(mgs,lh,1)
            vtxbar(mgs,lh,2) = graupelfallfac*vtxbar(mgs,lh,2)
            vtxbar(mgs,lh,3) = graupelfallfac*vtxbar(mgs,lh,3)
            axx(mgs,lh) = graupelfallfac*axx(mgs,lh)
          enddo
        endif

        if ( lhl > 0 .and. hailfallfac /= 1.0 ) then
          do mgs = 1,ngscnt
            vtxbar(mgs,lhl,1) = hailfallfac*vtxbar(mgs,lhl,1)
            vtxbar(mgs,lhl,2) = hailfallfac*vtxbar(mgs,lhl,2)
            vtxbar(mgs,lhl,3) = hailfallfac*vtxbar(mgs,lhl,3)
            axx(mgs,lhl) = hailfallfac*axx(mgs,lhl)
          enddo
        endif
      
      if ( ndebug1 .gt. 0 ) write(0,*) 'SETVTZ: end OF ROUTINE'

!############ SETVTZ ############################

      return
      end subroutine setvtz
!--------------------------------------------------------------------------

!
! ##############################################################################

!
!  subroutine to calculate fall speeds of hydrometeors
!

      subroutine ziegfall1d(nx,ny,nz,nor,norz,na,dtp,jgs,ixcol, &
     &  xvt, rhovtzx,                                           &
     &  an,dn,ipconc0,t0,t7,cwmasn,cwmasx,       &
     &  cwradn,                                   &
     &  qxmin,xdnmx,xdnmn,cdx,cno,xdn0,xvmn,xvmx,  &
     &  ngs,qx,qxw,cx,xv,vtxbar,xmas,xdn,xdia,vx,alpha,zx,igs,kgs, &
     &  rho0,temcg,temg,rhovt,cwnc,cinc,fadvisc,cwdia,cipmas,cnina,cimas, &
     &  cnostmp,                     &
     &  infdo,ildo,timesetvt)

! 12.16.2005: .F version use in transitional SWM model
!
! 10.10.2003: Added cimn and cimx to setting for cci and cip.
!
! TO DO LIST:
!
! need to set up values for:
!     :  cipdia,cidia,cwdia,cwmas,vtwbar,
!     :  rho0,temcg,cip,cci
!
! and need to put fallspeed values in cwvt etc.
!
      
      implicit none
      integer ng1
      parameter(ng1 = 1)
      
      integer, intent(in) :: ixcol ! which column to return
      integer, intent(in) :: ildo
      
      integer nx,ny,nz,nor,norz,ngt,jgs,na
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real dn(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real dtp,dtz1
      
      real :: rhovtzx(nz,nx)
      
      integer ndebugzf
      parameter (ndebugzf = 0)

      integer ix,jy,kz,i,j,k,il
      integer infdo
!
!
      real xvt(nz+1,nx,3,lc:lhab) ! 1=mass-weighted, 2=number-weighted

      real qxmin(lc:lhab)
      real xdn0(lc:lhab)
      real xvmn(lc:lhab), xvmx(lc:lhab)
      double precision,optional :: timesetvt

      integer :: ngs
      integer :: ngscnt,mgs,ipconc0
!      parameter ( ngs=200 )
      
      real ::  qx(ngs,lv:lhab) 
      real ::  qxw(ngs,ls:lhab) 
      real ::  cx(ngs,lc:lhab) 
      real ::  xv(ngs,lc:lhab) 
      real ::  vtxbar(ngs,lc:lhab,3) 
      real ::  xmas(ngs,lc:lhab) 
      real ::  xdn(ngs,lc:lhab) 
      real ::  cdxgs(ngs,lc:lhab) 
      real ::  xdia(ngs,lc:lhab,3) 
      real ::  vx(ngs,li:lhab) 
      real ::  alpha(ngs,lc:lhab) 
      real ::  zx(ngs,lr:lhab) 

      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real :: axx(ngs,lh:lhab), bxx(ngs,lh:lhab)
!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)

!
!   drag coefficients
!
      real cdx(lc:lhab)
!
! Fixed intercept values for single moment scheme
!
      real cno(lc:lhab)
      
      real cwccn0,cwmasn,cwmasx,cwradn
!      real cwc0

      integer nxmpb,nzmpb,nxz,numgs,inumgs
      integer kstag
      parameter (kstag=1)

      integer igs(ngs),kgs(ngs)
      
      real rho0(ngs),temcg(ngs)

      real temg(ngs)
      
      real rhovt(ngs)
      
      real cwnc(ngs),cinc(ngs)
      real fadvisc(ngs),cwdia(ngs),cipmas(ngs)
      
!      real cimasn,cimasx,
      real :: cnina(ngs),cimas(ngs)
      
      real :: cnostmp(ngs)

!      real pii
!
!
!  general constants for microphysics
!

! 
! Miscellaneous
!
      
      logical flag
      logical ldoliq
      
    
      real chw, qr, z, rd, alp, z1, g1, vr, nrx, tmp
      
      real vtmax
      real xvbarmax

      real, parameter ::  c1r=19.0, c2r=0.6, c3r=1.8, c4r=17.0   ! rain
      real, parameter ::  c1h=5.5, c2h=0.7, c3h=4.5, c4h=8.5   ! Graupel
      real, parameter ::  c1hl=3.7, c2hl=0.3, c3hl=9.0, c4hl=6.5, c5hl=1.0, c6hl=6.5 ! Hail

      integer l1, l2
      
      double precision :: dpt1, dpt2


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES 

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze

      logical :: debug_mpi = .false.


      if (ndebugzf .gt. 0 ) write(0,*) "ZIEGFALL: ENTERED subroutine"

! #####################################################################
! BEGIN EXECUTABLE
! #####################################################################
!

!  constants
!

      ldoliq = .false.
      if ( ls .gt. 1 ) then
      do il = ls,lhab
        ldoliq = ldoliq .or. ( lliq(il) .gt. 1 )
      enddo
      endif
      
!      poo = 1.0e+05
!      cp608 = 0.608
!      cp = 1004.0
!      cv = 717.0
!      dnz00 = 1.225
!      rho00 = 1.225
!      cs = 4.83607122
!      ds = 0.25
!  new values for  cs and ds
!      cs = 12.42
!      ds = 0.42
!      pi = 4.0*atan(1.0)
!      pii = piinv ! 1./pi
!      pid4 = pi/4.0 
!      qccrit = 2.0e-03
!      qscrit = 6.0e-04
!      cwc0 = pii
      
!
!
!  general constants for microphysics
!
      
!
!  ci constants in mks units
!
!      cimasn = 6.88e-13 
!      cimasx = 1.0e-8
!
!  Set terminal velocities...
!    also set drag coefficients
!
      jy = jgs
      nxmpb = ixcol
      nzmpb = 1
      nxz = 1*nz
!      ngs = nz
      numgs = 1

      if ( ildo == 0 ) then
        l1 = lc
        l2 = lhab
      else
        l1 = ildo
        l2 = ildo
      endif


      do inumgs = 1,numgs
       ngscnt = 0


       do kz = nzmpb,nz
        do ix = ixcol,ixcol
        flag = .false.

        
        do il = l1,l2
          flag =  flag .or. ( an(ix,jy,kz,il)  .gt. qxmin(il) ) 
        enddo

        if ( flag ) then
! load temp quantities

        ngscnt = ngscnt + 1
        igs(ngscnt) = ix
        kgs(ngscnt) = kz
        if ( ngscnt .eq. ngs ) goto 1100
        end if
        end do !!ix
        nxmpb = 1
       end do !! kz

!      if ( jy .eq. (ny-jstag) ) iend = 1

 1100 continue

      if ( ngscnt .eq. 0 ) go to 9998
!
!  set temporaries for microphysics variables
!


!
!  Reconstruct various quantities 
!
      do mgs = 1,ngscnt

       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhovt(mgs) = rhovtzx(kgs(mgs),ixcol) !  Sqrt(rho00/rho0(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
       temcg(mgs) = temg(mgs) - tfr

        
!
      end do
!
! only need fadvisc for 
      if ( lc .gt. 1 .and. (ildo == 0 .or. ildo == lc ) ) then
        do mgs = 1,ngscnt
         fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* &
     &   (temg(mgs)/296.0)**(1.5)
        end do
      endif

      if ( ipconc .eq. 0 ) then
      do mgs = 1,ngscnt
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      end do
      endif


      if ( ildo > 0 ) then
        vtxbar(:,ildo,:) = 0.0
      else
        vtxbar(:,:,:) = 0.0
      endif
      
!      do mgs = 1,ngscnt
!        qx(mgs,lv) = max(an(igs(mgs),jy,kgs(mgs),lv), 0.0) 
!      ENDDO
      do il = l1,l2
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      enddo
      end do
      
      cnostmp(:) = cno(ls)
      if ( ipconc < 1 .and. lwsm6 .and. (ildo == 0 .or. ildo == ls )) then
        do mgs = 1,ngscnt
          tmp = Min( 0.0, temcg(mgs) )
          cnostmp(mgs) = Min( 2.e8, 2.e6*exp(0.12*tmp) )
        enddo
      endif


!
!  set concentrations
!
      cx(:,:) = 0.0
      
      if ( ipconc .ge. 1 .and. li .gt. 1 .and. (ildo == 0 .or. ildo == li ) ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 .and. lc .gt. 1 .and. (ildo == 0 .or. ildo == lc ) ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
!        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
       end do
      end if
      if ( ipconc .ge. 3 .and. lr .gt. 1 .and. (ildo == 0 .or. ildo == lr ) ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
!        IF ( qx(mgs,lr) .le. qxmin(lr) ) THEN
!        ELSE
!          cx(mgs,lr) = Max( 0.0, cx(mgs,lr) )
!        ENDIF
       end do
      end if
      if ( ipconc .ge. 4  .and. ls .gt. 1 .and. (ildo == 0 .or. ildo == ls ) ) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
!        IF ( qx(mgs,ls) .le. qxmin(ls) ) THEN
!        ELSE
!          cx(mgs,ls) = Max( 0.0, cx(mgs,ls) )
!        ENDIF
       end do
      end if

      if ( ipconc .ge. 5  .and. lh .gt. 1 .and. (ildo == 0 .or. ildo == lh ) ) then
       do mgs = 1,ngscnt

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
!        IF ( qx(mgs,lh) .le. qxmin(lh) ) THEN
!        ELSE
!          cx(mgs,lh) = Max( 0.0, cx(mgs,lh) )
!        ENDIF

       end do
      endif

      if ( ipconc .ge. 5  .and. lhl .gt. 1 .and. (ildo == 0 .or. ildo == lhl ) ) then
       do mgs = 1,ngscnt

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
!        IF ( qx(mgs,lhl) .le. qxmin(lhl) ) THEN
!          cx(mgs,lhl) = 0.0
!        ELSEIF ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lhl) ) THEN
!          qx(mgs,lhl) = 0.0
!        ELSE
!          cx(mgs,lhl) = Max( 0.0, cx(mgs,lhl) )
!        ENDIF

       end do
      end if
       
      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
!        IF ( ls .gt. 1 .and. lvs .eq. 0 ) xdn(mgs,ls) = xdn0(ls)
!        IF ( lh .gt. 1 .and. lvh .eq. 0 ) xdn(mgs,lh) = xdn0(lh)
        if ( li .gt. 1 )  xdn(mgs,li) = xdn0(li)
        if ( ls .gt. 1 )  xdn(mgs,ls) = xdn0(ls)
        if ( lh .gt. 1 )  xdn(mgs,lh) = xdn0(lh)
        if ( lhl .gt. 1 ) xdn(mgs,lhl) = xdn0(lhl)
      end do

!
! Set mean particle volume
!
      if ( ldovol .and. (ildo == 0 .or. ildo >= li ) ) then
      
      vx(:,:) = 0.0
      
       do il = l1,l2
        
        if ( lvol(il) .ge. 1 ) then
        
          do mgs = 1,ngscnt
            vx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lvol(il)), 0.0)
            if ( vx(mgs,il) .gt. rho0(mgs)*qxmin(il)*1.e-3 .and. qx(mgs,il) .gt. qxmin(il) ) then
              xdn(mgs,il) = Min( xdnmx(il), Max( xdnmn(il), rho0(mgs)*qx(mgs,il)/vx(mgs,il) ) )
            endif
          enddo
          
        endif
      
       enddo
      
      endif

      do il = lg,lhab
      do mgs = 1,ngscnt
        alpha(mgs,il) = dnu(il)
      enddo
      enddo
      
      if ( imurain == 1 ) then
        alpha(:,lr) = alphar
      elseif ( imurain == 3 ) then
        alpha(:,lr) = xnu(lr)
      endif


      if ( ipconc == 5 .and. imydiagalpha > 0 ) then
        do mgs = 1,ngscnt
          if ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > cxmin ) then
             xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))            ! 
             xdia(mgs,lr,3) = (xv(mgs,lr)*6.0*cwc1)**(1./3.) 
             alpha(mgs,lr) = Min(alphamax, c1r*tanh(c2r*(xdia(mgs,lr,3)*1000. - c3r)) + c4r)
          endif
          if ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) > cxmin ) then
             xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh))            ! 
             xdia(mgs,lh,3) = (xv(mgs,lh)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
             alpha(mgs,lh) = Min(alphamax, c1h*tanh(c2h*(xdia(mgs,lh,3)*1000. - c3h)) + c4h)
          endif
!        alpha(:,lr) = 0. ! 10.
!        alpha(:,lh) = 0. ! 10.
          if ( lhl > 0 ) then
          if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) > cxmin ) then
             xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl))            ! 
             xdia(mgs,lhl,3) = (xv(mgs,lhl)*6.*piinv)**(1./3.)
             if ( xdia(mgs,lhl,3) < 0.008 ) then
               alpha(mgs,lhl) = Min(alphamax, c1hl*tanh(c2hl*(xdia(mgs,lhl,3)*1000. - c3hl)) + c4hl)
             else
               alpha(mgs,lhl) = Min(alphamax, c5hl*xdia(mgs,lhl,3)*1000. + c6hl)
             endif
          endif
          endif
        enddo
      endif


!
! Set 6th moments
!
      if ( ipconc .ge. 6 .or. lzr > 1) then
      
      zx(:,:) = 0.0
      
!      DO il = lr,lhab
       do il = l1,l2
        
        if ( lz(il) .ge. 1 ) then
        
          do mgs = 1,ngscnt
            zx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0)
          enddo
          
        
        endif
      
       enddo
      
      endif
       



       
!  Find shape parameter rain


     if ( lz(lr) > 1 .and. (ildo == 0 .or. ildo == lr ) .and. imurain == 3  ) then ! { RAIN SHAPE PARAM
          il = lr
          do mgs = 1,ngscnt
         
         if ( iresetmoments == 1 .or. iresetmoments == il  ) then
!         IF (  .false. .and. zx(mgs,lr) <= zxmin ) THEN
         if ( zx(mgs,lr) <= zxmin ) then
           qx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
!         ELSEIF ( zx(mgs,lr) <= 0.0 .and. cx(mgs,lr) > 0.0 .and. qx(mgs,il) .gt. qxmin(il)) THEN
!           write(91,*) 'ZF: overdepletion of Zr: z,c,q = ',zx(mgs,il),cx(mgs,il),qx(mgs,il)
         elseif ( cx(mgs,lr) <= cxmin ) then
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         endif
         endif
         
          
         
         if ( qx(mgs,lr) .gt. qxmin(lr) ) then

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        if ( xv(mgs,lr) .gt. xvmx(lr) ) then
!          tmp = cx(mgs,lr)
!          xv(mgs,lr) = xvmx(lr)
!          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
!          an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
!          IF ( tmp < cx(mgs,il) ) THEN ! breakup
!             g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
!!             zx(mgs,lr) = zx(mgs,lr) + g1*(rho0(mgs)/(1000.))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
!!             an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
!          ENDIF
        elseif ( xv(mgs,lr) .lt. xvmn(lr) ) then
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
          an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
        endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)

            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)

!            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-9,cx(mgs,lr)))
!            vr = xv(mgs,lr)

!             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z

            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              write(91,*) 'alpha = ',alpha(mgs,il)
             if ( qx(mgs,il) < 1.e-8 ) then
             qx(mgs,il) = 0.0
             an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
             an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
             else
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
             endif
          endif
          
          if ( zx(mgs,lr) > 0.0 ) then
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*Max(1.0e-9,cx(mgs,lr)))
            vr = xv(mgs,lr)
!            z = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           if ( z .gt. 0.0 ) then
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           do i = 1,20
!            IF ( 100.*Abs(alp - alpha(mgs,lr))/Abs(alpha(mgs,lr)) .lt. 1. ) EXIT
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
!           write(0,*) 'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo
!           write(0,*) 'kz, alp, alpha(kz) = ',kz,alp,alpha(mgs,lr),qr*1000,z*1.e18,vr,nrx


! check for artificial breakup (rain larger than allowed max size)
        if (  xv(mgs,il) .gt. xvmx(il) ) then
          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          if ( tmp < cx(mgs,il) ) then ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           do i = 1,20
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo

            
          endif
        endif

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
!           IF ( alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax ) THEN
           if ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) then

            if ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( rescale_low_alphar .and. alp <= rnumin ) then

             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             
             endif
           endif

           endif
          endif
           
          else
          
           zx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
          
          endif
          
          enddo
        endif ! }
        

      if ( ipconc .ge. 6 ) then

!  Find shape parameters for graupel,hail

        do il = lr,lhab
        
        if ( lz(il) .gt. 1 .and. (ildo == 0 .or. ildo == il ) .and. ( .not. ( il == lr .and. imurain == 3 )) ) then
        
        do mgs = 1,ngscnt

         if ( iresetmoments == 1 .or. iresetmoments == il  .or. iresetmoments == -1 ) then
         if ( zx(mgs,il) <= zxmin ) then !  .and. qx(mgs,il) > 0.05e-3 ) then
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         elseif ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         elseif ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) then !  .and. qx(mgs,il) > 0.05e-3  ) then
!!            write(91,*) 'cx=0; qx,zx = ',1000.*qx(mgs,il),1.e18*zx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif
         endif

         if (  zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif

         if ( qx(mgs,il) .gt. qxmin(il) ) then

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        if ( xv(mgs,il) .lt. xvmn(il)  ) then
!          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
!          IF ( tmp < cx(mgs,il) ) THEN ! breakup
!            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
!             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
!             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
!            
!          ENDIF
        endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
            an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
!              write(91,*) 'ziegfall: something screwy with moments: il = ',il
!              write(91,*) 'q,n,z = ', 1.e3*qx(mgs,il),cx(mgs,il),zx(mgs,il)
!              write(91,*) 'alpha = ',alpha(mgs,il)

             if ( qx(mgs,il) < 1.e-8 ) then
             qx(mgs,il) = 0.0
             an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
             an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
             else
!              write(0,*) 'alpha = ',alpha(mgs,il)
         ! set values according to dBZ of -10
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
            endif
          endif
         endif

        if ( qx(mgs,il) .gt. qxmin(il) .and. cx(mgs,il) .gt. cxmin ) then
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          if ( zx(mgs,il) .gt. 0. ) then
           
!            rd = z*(pi/6.*1000.)**2*chw/(0.224*(dn(igs(mgs),jy,kgs(mgs))*qr)**2)
            rd = z*(pi/6.*xdn(mgs,il))**2*chw/((dn(igs(mgs),jy,kgs(mgs))*qr)**2)

           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           do i = 1,10
            if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
!           write(0,*) 'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo



! check for artificial breakup (graupel/hail larger than allowed max size)
        
        if ( imaxdiaopt == 1 ) then
          xvbarmax = xvmx(il) 
        elseif ( imaxdiaopt == 2 ) then ! test against maximum mass diameter
          xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        elseif ( imaxdiaopt == 3 ) then ! test against mass-weighted diameter
          xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        endif
        
        if (  xv(mgs,il) .gt. xvbarmax ) then
          tmp = cx(mgs,il)
          xv(mgs,il) = Min( xvbarmax, Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          if ( tmp < cx(mgs,il) ) then ! breakup
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

            rd = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
           do i = 1,10
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo

            
          endif
        endif
           
!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
           if ( (rescale_low_alpha .or. rescale_high_alpha ) .and.  &
     &        ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) then

             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))

            if ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( rescale_low_alpha .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) ) then

!!             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*( 0.224*qr)*qr/chw
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
            endif
           endif
          else
          endif
        endif
        enddo ! mgs
        
        endif ! lz(il) .gt. 1
        
        enddo ! il

!      CALL cld_cpu('Z-MOMENT-ZFAll')  
          
      endif

      if ( lzhl > 1 ) then
        if ( lhl .gt. 1 ) then
        
        endif
      endif



!
!  Set density
!
      if (ndebugzf .gt. 0 ) write(0,*)  'ZIEGFALL: call setvtz'
!
      
      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,cnostmp,   &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,        &
     &                 ipconc,ndebugzf,ngs,nz,kgs,fadvisc, &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,    &
     &                 itype1,itype2,temcg,infdo,alpha,ildo,axx,bxx)
!     &                 itype1,itype2,temcg,infdo,alpha,ildo,axh,bxh,axhl,bxhl)



!
! put fall speeds into the x-z arrays
!
      do il = l1,l2
      do mgs = 1,ngscnt
       
       vtmax = 150.0

       
       if ( vtxbar(mgs,il,2) .gt. vtxbar(mgs,il,1)  .or. &
     &      ( vtxbar(mgs,il,1) .gt. vtxbar(mgs,il,3) .and. vtxbar(mgs,il,3) > 0.0) ) then
          
          
!          IF ( qx(mgs,il) > 1.e-4 .and.  &
!     &        .not. ( il == lr .and. 1.e3*xdia(mgs,il,3) > 5.0 ) ) THEN
!          write(0,*) 'infdo,mgs = ',infdo,lzr,mgs
!          write(0,*) 'Moment problem with vtxbar for il at i,j,k = ',il,igs(mgs),jy,kgs(mgs)
!          write(0,*) 'nx,ny,nz,ng = ',nx,ny,nz,nor
!          write(0,*) 'cwmasn,cwmasx = ',cwmasn,cwmasx
!          write(0,*) 'vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
!          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
!          IF ( il .ge. lr  .and. lz(il) > 1 ) write(0,*) 'z = ', zx(mgs,il)
!          IF ( il .ge. lg .or. il == lr ) THEN
!            write(0,*) 'alpha = ',alpha(mgs,il)
!          ENDIF
!          ENDIF
          
          vtxbar(mgs,il,1) = Max( vtxbar(mgs,il,1), vtxbar(mgs,il,2) )
          vtxbar(mgs,il,3) = Max( vtxbar(mgs,il,3), vtxbar(mgs,il,1) )
          
       endif

       
       if ( vtxbar(mgs,il,1) .gt. vtmax .or. vtxbar(mgs,il,2) .gt. vtmax .or. &
     &      vtxbar(mgs,il,3) .gt. vtmax ) then
       
!        IF ( ndebugzf >= 0 .and.  1.e3*qx(mgs,il) > 0.1 ) THEN
!          write(0,*) 'infdo = ',infdo
!          write(0,*) 'Problem with vtxbar for il at i,j,k = ',il,igs(mgs),jy,kgs(mgs)
!          write(0,*) 'nx,ny,nz,ng = ',nx,ny,nz,nor
!          write(0,*) 'cwmasn,cwmasx = ',cwmasn,cwmasx
!          write(0,*) 'vt1,2,3 = ',vtxbar(mgs,il,1),vtxbar(mgs,il,2),vtxbar(mgs,il,3)
!          write(0,*) 'q,n,d = ', 1.e3*qx(mgs,il),cx(mgs,il),1.e3*xdia(mgs,il,3)
!          IF ( il .ge. lr  .and. lz(il) > 1 ) write(0,*) 'z = ', zx(mgs,il)
!          IF ( il .ge. lg ) THEN
!            write(0,*) 'alpha = ',alpha(mgs,il)
!          ENDIF
!        ENDIF
        vtxbar(mgs,il,1) = Min(vtmax,vtxbar(mgs,il,1) )
        vtxbar(mgs,il,2) = Min(vtmax,vtxbar(mgs,il,2) )
        vtxbar(mgs,il,3) = Min(vtmax,vtxbar(mgs,il,3) )
        
!        call commasmpi_abort()
       endif


       xvt(kgs(mgs),igs(mgs),1,il) = vtxbar(mgs,il,1)
       xvt(kgs(mgs),igs(mgs),2,il) = vtxbar(mgs,il,2)
       if ( infdo .ge. 2 ) then
       xvt(kgs(mgs),igs(mgs),3,il) = vtxbar(mgs,il,3)
       else
       xvt(kgs(mgs),igs(mgs),3,il) = 0.0
       endif

!       xvt(kgs(mgs),igs(mgs),2,il) = xvt(kgs(mgs),igs(mgs),1,il)

      enddo
      enddo


      if (ndebugzf .gt. 0 ) write(0,*)  'ZIEGFALL: COPIED FALL SPEEDS'



 9998 continue

      if (ndebugzf .gt. 0 ) write(0,*)  'ZIEGFALL: DONE WITH LOOP'

      if ( kz .gt. nz-1 ) then
        go to 1200
      else
        nzmpb = kz 
      end if

      if (ndebugzf .gt. 0 ) write(0,*) 'ZIEGFALL: SET NZMPB'

      end do !! inumgs

      if (ndebugzf .gt. 0 ) write(0,*) 'ZIEGFALL: SET NXMPB'

 1200 continue


!       ENDDO ! ix
!      ENDDO ! kz


      if (ndebugzf .gt. 0 ) write(0,*) "ZIEGFALL: EXITING subroutine"


      return
      end subroutine ziegfall1d

! #####################################################################
! #####################################################################


! #####################################################################
! #####################################################################

! ##############################################################################
      subroutine radardd02(nx,ny,nz,nor,na,an,temk,         &
     &    dbz,db,nzdbz,cnoh0t,hwdn1t,ipconc,ke_diag, iunit)
!
! 11.13.2005: Changed values of indices for reordering of lip
!
! 07.13.2005: Fixed an error where cnoh was being used for graupel and frozen drops
!
! 01.24.2005: add ice crystal reflectivity using parameterization of
!             Heymsfield (JAS, 1977).  Could also try Ferrier for this, too.
!
!  09.28.2002 Test alterations for dry ice following Ferrier (1994)
!      for equivalent melted diameter reflectivity.
!      Converted to Fortran by ERM.
!      
!Date: Tue, 21 Nov 2000 10:13:36 -0600 (CST)
!From: Matthew Gilmore <gilmore@hesston.met.tamu.edu>
!
!PRO RF_SPEC ; Computes Radar Reflectivity
!COMMON MAINB, data, x1d, y1d, z1d, iconst, rconst, labels, nx, ny, nz, dshft
!
!;MODIFICATION HISTORY
!; 5/99  -Svelta Veleva introduces variable dielf (const_ki_x) as a (weak)
!;   function of density.  This leads to slight modification of dielf such
!;   that the snow reflectivity is slightly increased - not a big effect.
!;   This is believed to be more accurate than assuming the dielectric
!;   constant for snow is the same as for hail in previous versions.
!
!;On 6/13/99 I added the VIL computation (k=0 in vil array)
!;On 6/15/99 I removed the number concentration dependencies as a function
!;           of temperature (only use for ferrier!)
!;On 6/15/99 I added the Composite reflectivity (k=1 in VIL array)
!;On 6/15/99 I added the Severe Hail Index computation (k=2 in vil array)
!;
!; 6/99 - Veleva and Seo argue that since graupel is more similar to
!;   snow (in number conc and size density) than it is to hail, we
!;   should not weight wetted graupel with the .95 exponent correction
!;   factor as in the case of hail.  An if-statement checks the size
!;   density for wet hail/graupel and treats them appropriately.
!;
!; 6/22/99 - Added function to compute height of max rf and 40 dbz echo top
!;           Also added vilqr which is the model vertical integrated liquid only
!;           using qr.  Will need to check...does not seem consistent with vilZ
!;


      implicit none
      
      character(LEN=15), parameter :: microp = 'ZVD'
      integer nx,ny,nz,nor,na,ngt
      integer nzdbz    !  how many levels actually to process
      
      integer ng1,n10
      integer iunit
      integer, parameter :: printyn = 0

      parameter( ng1 = 1 )
      
      real cnoh0t,hwdn1t
      integer ke_diag
      integer ipconc
      real vr


      integer imapz,mzdist
      
      integer vzflag
      integer, parameter :: norz = 3
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,na)
      real db(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)  ! air density
!      real gt(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor,ngt)
      real temk(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)  ! air temperature (kelvin)
      real dbz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)   ! reflectivity
      real gz(-nor+1:nz+nor) ! ,z1d(-nor+1:nz+nor,4)
      
!      real g,rgas,eta,inveta
      real cr1, cr2 ,  hwdnsq,swdnsq
      real rwdnsq, dhmin, qrmin, qsmin, qhmin, qhlmin, tfr, tfrh, zrc
      real reflectmin,  kw_sq
      real const_ki_sn, const_ki_h, ki_sq_sn
      real ki_sq_h, dielf_sn, dielf_h
      real pi
      logical ltest

!  Other data arrays
       real gtmp     (nx,nz)
       real dtmp     (nx,nz)
       real tmp

       real*8 dtmps, dtmpr, dtmph, dtmphl, g1, zx, ze, x

       integer i,j,k,ix,jy,kz,ihcnt

        real*8 xcnoh, xcnos, dadh, dads, zhdryc, zsdryc, zhwetc,zswetc
        real*8 dadr
        real dbzmax,dbzmin
        parameter ( dbzmin = 0 )

      real cnow,cnoi,cnoip,cnoir,cnor,cnos
      real cnogl,cnogm,cnogh,cnof,cnoh,cnohl

      real swdn, rwdn ,hwdn,gldn,gmdn,ghdn,fwdn,hldn
      real swdn0

      real rwdnmx,cwdnmx,cidnmx,xidnmx,swdnmx,gldnmx,gmdnmx
      real ghdnmx,fwdnmx,hwdnmx,hldnmx
      real rwdnmn,cwdnmn,cidnmn,xidnmn,swdnmn,gldnmn,gmdnmn
      real ghdnmn,fwdnmn,hwdnmn,hldnmn
 
      real gldnsq,gmdnsq,ghdnsq,fwdnsq,hldnsq

      real dadgl,dadgm,dadgh,dadhl,dadf
      real zgldryc,zglwetc,zgmdryc, zgmwetc,zghdryc,zghwetc
      real zhldryc,zhlwetc,zfdryc,zfwetc

      real dielf_gl,dielf_gm,dielf_gh,dielf_hl,dielf_fw
      
      integer imx,jmx,kmx
      
      real swdia,gldia,gmdia,ghdia,fwdia,hwdia,hldia
      
      real csw,cgl,cgm,cgh,cfw,chw,chl
      real xvs,xvgl,xvgm,xvgh,xvf,xvh,xvhl
      
      real cwc0
      integer izieg
      integer ice10
      real rhos
      parameter ( rhos = 0.1 )
      
      real qxw,qxw1    ! temp value for liquid water on ice mixing ratio
      real :: dnsnow
      real qh

      real, parameter :: cwmasn = 5.23e-13   ! minimum mass, defined by radius of 5.0e-6
      real, parameter :: cwmasx = 5.25e-10   ! maximum mass, defined by radius of 50.0e-6
      real, parameter :: cwradn = 5.0e-6     ! minimum radius

      real cwnccn(nz)
      
      real :: vzsnow, vzrain, vzgraupel, vzhail
      real :: ksq
      real :: dtp


! #########################################################################      

      vzflag = 0
      
      izieg = 0
      ice10 = 0
!      g=9.806                 ! g: gravity constant
!      rgas=287.04             ! rgas: gas constant for dry air
!      rcp=rgas/cp             ! rcp: gamma constant
!      eta=0.622
!      inveta = 1./eta
!      rcpinv = 1./rcp
!      cpr=cp/rgas
!      cvr=cv/rgas
      pi = 4.0*ATan(1.)
      cwc0 = piinv ! 1./pi ! 6.0/pi
      
      cnoh = cnoh0t
      hwdn = hwdn1t

      rwdn = 1000.0
      swdn = 100.0

      qrmin = 1.0e-05
      qsmin = 1.0e-06
      qhmin = 1.0e-05

!
!  default slope intercepts
!
      cnow  = 1.0e+08
      cnoi  = 1.0e+08
      cnoip = 1.0e+08 
      cnoir = 1.0e+08 
      cnor  = 8.0e+06 
      cnos  = 8.0e+06 
      cnogl = 4.0e+05 
      cnogm = 4.0e+05 
      cnogh = 4.0e+05 
      cnof  = 4.0e+05
      cnohl = 1.0e+03


      imx = 1
      jmx = 1
      kmx = 1
      i = 1


       if ( microp(1:4) .eq. 'ZIEG' ) then !  na .ge. 14 .and. ipconc .ge. 3 ) then 

!        write(0,*)  'Set reflectivity for ZIEG'
         izieg = 1

         hwdn = hwdn1t ! 500.


         cnor  = cno(lr)
         cnos  = cno(ls)
         cnoh  = cno(lh)
         qrmin = qxmin(lr)
         qsmin = qxmin(ls)
         qhmin = qxmin(lh)
         if ( lhl .gt. 1 ) then
            cnohl  = cno(lhl)
            qhlmin = qxmin(lhl)
         endif

       elseif ( microp(1:3) .eq. 'ZVD' ) then !  na .ge. 14 .and. ipconc .ge. 3 ) then 

         izieg = 1
         
         swdn0 = swdn

         cnor  = cno(lr)
         cnos  = cno(ls)
         cnoh  = cno(lh)
         
         qrmin = qxmin(lr)
         qsmin = qxmin(ls)
         qhmin = qxmin(lh)
         if ( lhl .gt. 1 ) then
            cnohl  = cno(lhl)
            qhlmin = qxmin(lhl)
         endif
!         write(*,*) 'radardbz: ',db(1,1,1),temk(1,1,1),an(1,1,1,lr),an(1,1,1,ls),an(1,1,1,lh)


        endif


!      cdx(lr) = 0.60
!      
!      IF ( lh > 1 ) THEN
!      cdx(lh) = 0.8 ! 1.0 ! 0.45
!      cdx(ls) = 2.00
!      ENDIF
!
!      IF ( lhl .gt. 1 ) cdx(lhl) = 0.45
!
!      xvmn(lc) = xvcmn
!      xvmn(lr) = xvrmn
!
!      xvmx(lc) = xvcmx
!      xvmx(lr) = xvrmx
!
!      IF ( lh > 1 ) THEN
!      xvmn(ls) = xvsmn
!      xvmn(lh) = xvhmn
!      xvmx(ls) = xvsmx
!      xvmx(lh) = xvhmx
!      ENDIF
!
!      IF ( lhl .gt. 1 ) THEN
!      xvmn(lhl) = xvhlmn
!      xvmx(lhl) = xvhlmx
!      ENDIF
!
!      xdnmx(lr) = 1000.0
!      xdnmx(lc) = 1000.0
!      IF ( lh > 1 ) THEN
!      xdnmx(li) =  917.0
!      xdnmx(ls) =  300.0
!      xdnmx(lh) =  900.0
!      ENDIF
!      IF ( lhl .gt. 1 ) xdnmx(lhl) = 900.0
!!
!      xdnmn(:) = 900.0
!      
!      xdnmn(lr) = 1000.0
!      xdnmn(lc) = 1000.0
!      IF ( lh > 1 ) THEN
!      xdnmn(li) =  100.0
!      xdnmn(ls) =  100.0
!      xdnmn(lh) =  hdnmn
!      ENDIF
!      IF ( lhl .gt. 1 ) xdnmn(lhl) = 500.0
!
!      xdn0(:) = 900.0
!      
!      xdn0(lc) = 1000.0
!      xdn0(lr) = 1000.0
!      IF ( lh > 1 ) THEN
!      xdn0(li) = 900.0
!      xdn0(ls) = 100.0 ! 100.0
!      xdn0(lh) = hwdn1t ! (0.5)*(xdnmn(lh)+xdnmx(lh))
!      ENDIF
!      IF ( lhl .gt. 1 ) xdn0(lhl) = 800.0

!
!  slope intercepts
!
!      cnow  = 1.0e+08
!      cnoi  = 1.0e+08
!      cnoip = 1.0e+08 
!      cnoir = 1.0e+08 
!      cnor  = 8.0e+06 
!      cnos  = 8.0e+06 
!      cnogl = 4.0e+05 
!      cnogm = 4.0e+05 
!      cnogh = 4.0e+05 
!      cnof  = 4.0e+05
!c      cnoh  = 4.0e+04
!      cnohl = 1.0e+03
!
!
!  density maximums and minimums
!
      rwdnmx = 1000.0
      cwdnmx = 1000.0
      cidnmx =  917.0
      xidnmx =  917.0
      swdnmx =  200.0
      gldnmx =  400.0
      gmdnmx =  600.0
      ghdnmx =  800.0
      fwdnmx =  900.0
      hwdnmx =  900.0
      hldnmx =  900.0
!
      rwdnmn = 1000.0
      cwdnmn = 1000.0
      xidnmn =  001.0
      cidnmn =  001.0
      swdnmn =  001.0
      gldnmn =  200.0
      gmdnmn =  400.0
      ghdnmn =  600.0
      fwdnmn =  700.0
      hwdnmn =  700.0
      hldnmn =  900.0

      
      gldn = (0.5)*(gldnmn+gldnmx)  ! 300.
      gmdn = (0.5)*(gmdnmn+gmdnmx)  ! 500.
      ghdn = (0.5)*(ghdnmn+ghdnmx)  ! 700.
      fwdn = (0.5)*(fwdnmn+fwdnmx)  ! 800.
      hldn = (0.5)*(hldnmn+hldnmx)  ! 900.


      cr1  = 7.2e+20
      cr2  = 7.295e+19
      hwdnsq = hwdn**2
      swdnsq = swdn**2
      rwdnsq = rwdn**2

      gldnsq = gldn**2
      gmdnsq = gmdn**2
      ghdnsq = ghdn**2
      fwdnsq = fwdn**2
      hldnsq = hldn**2
      
      dhmin = 0.005
      tfr   = 273.16
      tfrh  = tfr - 8.0
      zrc   = cr1*cnor
      reflectmin = 0.0
      kw_sq = 0.93
      dbzmax = dbzmin
      
      ihcnt=0

            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Dielectric Factor  - Formulas implemented by Svetla Veleva
!                       following Battan, "Radar Meteorology" - p. 40
!  The result of these calculations is that the dielf numerator (ki_sq) without
!  the density ratio is  .2116 for hail if using 917 density and .25 for
!  snow if using 220 density.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      const_ki_sn = 0.5 - (0.5-0.46)/(917.-220.)*(swdn-220.)
      const_ki_h  = 0.5 - (0.5-0.46)/(917.-220.)*(hwdn-220.)
      ki_sq_sn = (swdnsq/rwdnsq) * const_ki_sn**2
      ki_sq_h  = (hwdnsq/rwdnsq) * const_ki_h**2
      dielf_sn = ki_sq_sn / kw_sq
      dielf_h  = ki_sq_h  / kw_sq
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Use the next line if you want to hardwire dielf for dry hail for both dry
!  snow and dry hail.
!  This would be equivalent to what Straka had originally. (i.e, .21/.93)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dielf_sn = (swdnsq/rwdnsq)*.21/ kw_sq
      dielf_h  = (hwdnsq/rwdnsq)*.21/ kw_sq

      dielf_gl  = (gldnsq/rwdnsq)*.21/ kw_sq
      dielf_gm  = (gmdnsq/rwdnsq)*.21/ kw_sq
      dielf_gh  = (ghdnsq/rwdnsq)*.21/ kw_sq
      dielf_hl  = (hldnsq/rwdnsq)*.21/ kw_sq
      dielf_fw  = (fwdnsq/rwdnsq)*.21/ kw_sq

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Notes on dielectric factors  - from Eun-Kyoung Seo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! constants for both snow and hail would be (x=s,h).....
!       xwdnsq/rwdnsq *0.21/kw_sq   ! Straka/Smith - the original
!       xwdnsq/rwdnsq *0.224        ! Ferrier - for particle sizes in equiv. drop diam
!       xwdnsq/rwdnsq *0.176/kw_sq  ! =0.189 in Smith - for particle sizes in equiv 
!                       ice spheres
!       xwdnsq/rwdnsq *0.208/kw_sq  ! Smith 1984 - for particle sizes in equiv melted drop diameter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! VIL algorithm constants
!      Ztop = 10.**(56./10)           !56 dbz is the max rf used by WATADS in cell vil


! Hail detection algorithm constants
!      ZL = 40.
!      ZU = 50.
!      Ho = 3400.  !WATADS Defaults
!      Hm20 = 6200.      !WATADS Defaults

!      DO kz = 1,Min(nzdbz,nz-1)

      do jy=1,1

        do kz = 1,ke_diag ! nz
         
          do ix=1,nx
            dbz(ix,jy,kz) = 0.0
                      
          vzsnow = 0.0
          vzrain = 0.0
          vzgraupel = 0.0
          vzhail = 0.0
          
          dtmph = 0.0
          dtmps = 0.0
          dtmphl = 0.0
          dtmpr = 0.0
           dadr = (db(ix,jy,kz)/(pi*rwdn*cnor))**(0.25)
!-----------------------------------------------------------------------
! Compute Rain Radar Reflectivity
!-----------------------------------------------------------------------
           
           dtmp(ix,kz) = 0.0
           gtmp(ix,kz) = 0.0
           if ( an(ix,jy,kz,lr) .ge. qrmin ) then
             if ( ipconc .le. 2 ) then
               gtmp(ix,kz) = dadr*an(ix,jy,kz,lr)**(0.25)
               dtmp(ix,kz) = zrc*gtmp(ix,kz)**7
             elseif ( lzr .gt. 1 ) then
               dtmp(ix,kz) = 1e18*an(ix,jy,kz,lzr)
             elseif ( an(ix,jy,kz,lnr) .gt. 1.e-3 ) then
               if ( imurain == 3 ) then
                 vr = db(ix,jy,kz)*an(ix,jy,kz,lr)/(1000.*an(ix,jy,kz,lnr))
                 dtmp(ix,kz) = 3.6e18*(rnu+2.)*an(ix,jy,kz,lnr)*vr**2/(rnu+1.)
               else ! imurain == 1
                g1 = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))
                zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lr))**2/an(ix,jy,kz,lnr)
                ze =1.e18*zx*(6./(pi*1000.))**2 ! note: using 1000. here for water density
                dtmp(ix,kz) = ze
               endif
             endif
             dtmpr = dtmp(ix,kz)
           endif
           
!-----------------------------------------------------------------------
! Compute snow and graupel reflectivity
!
! Lou modified to look at parcel temperature rather than base state
!-----------------------------------------------------------------------

          if( lhab .gt. lr ) then

!    qs2d   = reform(data[*,*,k,10],[nx*ny])
!    qh2d   = reform(data[*,*,k,11],[nx*ny])

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Only use the following lines if running Straka GEMS microphysics
!  (Sam 1-d version modified by L Wicker does not use this)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ;xcnoh    = cnoh*exp(-0.025*(temp-tfr))
!    ;xcnos    = cnos*exp(-0.038*(temp-tfr))
!    ;good = where(temp GT tfr, n_elements)
!    ;IF n_elements NE 0 THEN xcnoh(good) = cnoh*exp(-0.075*(temp(good)-tfr))
!    ;IF n_elements NE 0 THEN xcnos(good) = cnos*exp(-0.088*(temp(good)-tfr))

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Only use the following lines if running Ferrier micro with No=No(T)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    ;  NOSE = -.15
!    ;  NOGE =  .0
!    ;  xcnoh = cnoh*(1.>exp(NOGE*(temp-tfr)) )
!    ;  xcnos = cnos*(1.>exp(NOSE*(temp-tfr)) )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Use the following lines if Nos and Noh are constant
!  (As in Svetla version of Ferrier, GCE Tao, and SAM 1-d)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xcnoh    = cnoh
        xcnos    = cnos

!
! Temporary fix for predicted number concentration -- need a 
! more appropriate reflectivity equation!
!
!        IF ( an(ix,jy,kz,lns) .lt. 0.1 ) THEN
!         swdia = (xvrmn*cwc0)**(1./3.)
!         xcnos = an(ix,jy,kz,ls)*db(ix,jy,kz)/(xvrmn*swdn*swdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         swdia =
!     >  (an(ix,jy,kz,ls)*db(ix,jy,kz)
!     > /(pi*swdn*an(ix,jy,kz,lns)))**(1./3.)
!
!        xcnos = an(ix,jy,kz,lns)/swdia
!        ENDIF

        if ( ls .gt. 1 ) then ! {
        
        if ( lvs .gt. 1 ) then
          if ( an(ix,jy,kz,lvs) .gt. 0.0 ) then
            swdn = db(ix,jy,kz)*an(ix,jy,kz,ls)/an(ix,jy,kz,lvs)
            swdn = Min( 300., Max( 100., swdn ) )
          else 
            swdn = swdn0
          endif
        
        endif 
        
        if ( ipconc .ge. 5 ) then ! {

        xvs = db(ix,jy,kz)*an(ix,jy,kz,ls)/  &
     &      (swdn*Max(1.0e-3,an(ix,jy,kz,lns)))
        if ( xvs .lt. xvsmn .or. xvs .gt. xvsmx ) then
          xvs = Min( xvsmx, Max( xvsmn,xvs ) )
          csw = db(ix,jy,kz)*an(ix,jy,kz,ls)/(xvs*swdn)
        endif

         swdia = (xvs*cwc0)**(1./3.)
         xcnos = an(ix,jy,kz,ls)*db(ix,jy,kz)/(xvs*swdn*swdia)
         
         endif ! }
         endif  ! }

!        IF ( an(ix,jy,kz,lnh) .lt. 0.1 ) THEN
!         hwdia = (xvrmn*cwc0)**(1./3.)
!         xcnoh = an(ix,jy,kz,lh)*db(ix,jy,kz)/(xvrmn*hwdn*hwdia)
!        ELSE
!      ! changed back to diameter of mean volume!!!
!         hwdia =
!     >  (an(ix,jy,kz,lh)*db(ix,jy,kz)
!     > /(pi*hwdn*an(ix,jy,kz,lnh)))**(1./3.)
!        
!         xcnoh = an(ix,jy,kz,lnh)/hwdia
!        ENDIF

        if ( lh .gt. 1 ) then ! {

        if ( lvh .gt. 1 ) then
          if ( an(ix,jy,kz,lvh) .gt. 0.0 ) then
            hwdn = db(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
            hwdn = Min( 900., Max( hdnmn, hwdn ) )
          else 
            hwdn = 500. ! hwdn1t
          endif
        else
          hwdn = hwdn1t
        endif 
        
        if ( ipconc .ge. 5 ) then ! {

        xvh = db(ix,jy,kz)*an(ix,jy,kz,lh)/       &
     &      (hwdn*Max(1.0e-3,an(ix,jy,kz,lnh)))
        if ( xvh .lt. xvhmn .or. xvh .gt. xvhmx ) then
          xvh = Min( xvhmx, Max( xvhmn,xvh ) )
          chw = db(ix,jy,kz)*an(ix,jy,kz,lh)/(xvh*hwdn)
        endif

         hwdia = (xvh*cwc0)**(1./3.)
         xcnoh = an(ix,jy,kz,lh)*db(ix,jy,kz)/(xvh*hwdn*hwdia)
         
        endif ! } ipconc .ge. 5
 
        endif ! }

        dadh = 0.0
        dadhl = 0.0
        dads = 0.0
        if ( xcnoh .gt. 0.0 ) then 
          dadh = ( db(ix,jy,kz) /(pi*hwdn*xcnoh) )**(.25)
          zhdryc = 0.224*cr2*(db(ix,jy,kz)/rwdn)**2/xcnoh ! dielf_h*cr1*xcnoh          ! SV - equiv formula as before but
                                        ! ratio of densities included in
                                        ! dielf_h rather than here following
                                        ! Battan.
        else
          dadh = 0.0
          zhdryc = 0.0
        endif
        
        if ( xcnos .gt. 0.0 ) then
          dads = ( db(ix,jy,kz) /(pi*swdn*xcnos) )**(.25)
          zsdryc = 0.224*cr2*(db(ix,jy,kz)/rwdn)**2/xcnos ! dielf_sn*cr1*xcnos         ! SV - similar change as above
        else
          dads = 0.0
          zsdryc = 0.0
        endif
        zhwetc = zhdryc ! cr1*xcnoh      !Hail/graupel version with .95 power bug removed
        zswetc = zsdryc ! cr1*xcnos
!           
! snow contribution
!
          if ( ls .gt. 1 ) then
          
          gtmp(ix,kz) = 0.0 
          qxw = 0.0 
          qxw1 = 0.0
          dtmps = 0.0
           if ( an(ix,jy,kz,ls) .ge. qsmin ) then !{
            if ( ipconc .ge. 4 ) then  ! (Ferrier 94) !{

             if (lsw .gt. 1) then 
               qxw = an(ix,jy,kz,lsw)
               qxw1 = 0.0
             elseif ( ( iusewetsnow == 1 .or. iusewetsnow == 3) .and. temk(ix,jy,kz) .gt. tfr+1. & 
     &                  .and. an(ix,jy,kz,ls) > an(ix,jy,kz,lr) .and. an(ix,jy,kz,lr) > qsmin) then
               qxw = Min(0.5*an(ix,jy,kz,ls), an(ix,jy,kz,lr))
               qxw1 = qxw
             endif

             vr = xvs ! db(ix,jy,kz)*an(ix,jy,kz,lr)/(1000.*an(ix,jy,kz,lnr))
!             gtmp(ix,kz) = 3.6e18*(0.243*rhos**2/0.93)*(snu+2.)*an(ix,jy,kz,lns)*vr**2/(snu+1.)
             
             ksq = 0.189 ! Smith (1984, JAMC) for equiv. ice sphere
             if ( an(ix,jy,kz,lns) .gt. 1.e-7 ) then
     !          IF ( .true. ) THEN
               if ( qxw > qsmin .or. iusewetsnow >= 2 ) then ! old version
!                gtmp(ix,kz) = 3.6e18*(snu+2.)*( 0.224*an(ix,jy,kz,ls) + 0.776*qxw)*an(ix,jy,kz,ls)/ &
!     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2
                gtmp(ix,kz) = 3.6e18*(snu+2.)*( 0.224*(an(ix,jy,kz,ls)+qxw1) + 0.776*qxw)*(an(ix,jy,kz,ls)+qxw1)/ &
     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2

               else ! new form using a mass relationship m = p d^2 (instead of d^3 -- Cox 1988 QJRMS) so that density depends on size
                    ! p = 0.106214 for m = p v^(2/3)
                 dnsnow = 0.0346159*sqrt(an(ix,jy,kz,lns)/(an(ix,jy,kz,ls)*db(ix,jy,kz)) )
                 if ( .true. .or. dnsnow < 900. ) then
                  gtmp(ix,kz) = 1.e18*323.3226* 0.106214**2*(ksq*an(ix,jy,kz,ls) + &
     &             (1.-ksq)*qxw)*an(ix,jy,kz,ls)*db(ix,jy,kz)**2*gsnow73/         &
     &                   (an(ix,jy,kz,lns)*(917.)**2* gsnow1*(1.0+snu)**(4./3.))
                 else ! otherwise small enough to assume ice spheres?
                  gtmp(ix,kz) = (36./pi**2) * 1.e18*(snu+2.)*( 0.224*(an(ix,jy,kz,ls)+qxw1) + 0.776*qxw)*(an(ix,jy,kz,ls)+qxw1)/ &
     &              (an(ix,jy,kz,lns)*(snu+1.)*rwdn**2)*db(ix,jy,kz)**2
                 endif

               endif
             
             endif
             
!             tmp = Min(1.0,1.e3*(an(ix,jy,kz,ls))*db(ix,jy,kz))
!             gtmp(ix,kz) = Max( 1.0*gtmp(ix,kz), 750.0*(tmp)**1.98)
             dtmps = gtmp(ix,kz)
             dtmp(ix,kz) = dtmp(ix,kz) + gtmp(ix,kz)
            else ! }{ single-moment snow:
             gtmp(ix,kz) = dads*an(ix,jy,kz,ls)**(0.25)
             
             if ( gtmp(ix,kz) .gt. 0.0 ) then !{
             dtmps = zsdryc*an(ix,jy,kz,ls)**2/gtmp(ix,kz)
             if ( temk(ix,jy,kz) .lt. tfr ) then
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                   zsdryc*an(ix,jy,kz,ls)**2/gtmp(ix,kz)
             else
               dtmp(ix,kz) = dtmp(ix,kz) +          &
     &                  zswetc*an(ix,jy,kz,ls)**2/gtmp(ix,kz)
             endif
             endif !}
            endif !}
           
           endif !}
           
           endif


!
! ice crystal contribution (Heymsfield, 1977, JAS)
!
         if ( li .gt. 1 .and. idbzci .ne. 0 ) then
          
           if ( idbzci == 1 .and. lni > 0 ) then
          ! assume spherical ice with density of 900 for dbz calc
            if ( an(ix,jy,kz,li) > qxmin(li) .and. an(ix,jy,kz,lni) > 1.0 ) then
                 vr = db(ix,jy,kz)*an(ix,jy,kz,li)/(900.*an(ix,jy,kz,lni))
                 dtmp(ix,kz) = dtmp(ix,kz) +  &
     &                 0.224*3.6e18*(cinu+2.)*an(ix,jy,kz,lni)*vr**2/(cinu+1.)*(900./1000.)**2
            endif

          elseif ( idbzci == 2 ) then
!
! ice crystal contribution (Heymsfield, 1977, JAS)
!
         gtmp(ix,kz) = 0.0 
           if ( an(ix,jy,kz,li) .ge. 0.1e-3 ) then
             gtmp(ix,kz) = Min(1.0,1.e3*(an(ix,jy,kz,li))*db(ix,jy,kz))
             dtmp(ix,kz) = dtmp(ix,kz) + 750.0*(gtmp(ix,kz))**1.98
           endif
           
          endif
        
        endif
          
!           
! graupel/hail contribution
!
         if ( lh .gt. 1 ) then ! {
           gtmp(ix,kz) = 0.0 
           dtmph = 0.0
           qxw = 0.0

          if ( izieg .ge. 1 .and. ipconc .ge. 5 ) then

           ltest = .false.
           if ( lzh > 1 ) then
             if ( an(ix,jy,kz,lzh) > 0.0 .and. an(ix,jy,kz,lh) > qhmin .and. &
                  an(ix,jy,kz,lnh) >= cxmin ) ltest = .true.
           endif
           
           if ( ltest .or. (an(ix,jy,kz,lh) .ge. qhmin .and. an(ix,jy,kz,lnh) .ge. cxmin )) then
            
            if ( lvh .gt. 1 ) then
             
             if ( an(ix,jy,kz,lvh) .gt. 0.0 ) then
               hwdn = db(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
               hwdn = Min( 900., Max( 100., hwdn ) )
              else 
               hwdn = 500. ! hwdn1t
              endif

             endif

             chw = an(ix,jy,kz,lnh)
            if ( chw .gt. 0.0 ) then                                         ! (Ferrier 94)
             xvh = db(ix,jy,kz)*an(ix,jy,kz,lh)/(hwdn*Max(1.0e-3,chw))
             if ( xvh .lt. xvhmn .or. xvh .gt. xvhmx ) then
              xvh = Min( xvhmx, Max( xvhmn,xvh ) )
              chw = db(ix,jy,kz)*an(ix,jy,kz,lh)/(xvh*hwdn)
             endif
             
             qh = an(ix,jy,kz,lh)
             
             if ( lhw .gt. 1 ) then
               if ( iusewetgraupel .eq. 1 ) then
                  qxw = an(ix,jy,kz,lhw)
               elseif ( iusewetgraupel .eq. 2 ) then
                  if ( hwdn .lt. 300. ) then
                    qxw = an(ix,jy,kz,lhw)
                  endif
               endif
             elseif ( iusewetgraupel .eq. 3 ) then
                  if ( hwdn .lt. 300. .and. temk(ix,jy,kz) > tfr .and. an(ix,jy,kz,lr) > qhmin ) then
                    qxw = Min( an(ix,jy,kz,lh), an(ix,jy,kz,lr))
                    qh = qh + qxw
                  endif
             elseif ( iusewetgraupel == 4 .and. temk(ix,jy,kz) .gt. tfr+0.25 .and. an(ix,jy,kz,lh) > an(ix,jy,kz,lr) &
     &              .and. an(ix,jy,kz,lr) > qhmin) then
               qxw = Min(0.5*an(ix,jy,kz,lh), an(ix,jy,kz,lr))
               qh = qh + qxw

             endif
             
             if ( lzh .gt. 1 ) then
              x = (0.224*qh +  0.776*qxw)/an(ix,jy,kz,lh)  ! weighted average of dielectric const
              dtmph = 1.e18*x*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
              dtmp(ix,kz) = dtmp(ix,kz) + dtmph
             else
             g1 = (6.0 + alphah)*(5.0 + alphah)*(4.0 + alphah)/((3.0 + alphah)*(2.0 + alphah)*(1.0 + alphah))
!             zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lh))**2/chw
!             ze = 0.224*1.e18*zx*(6./(pi*1000.))**2
             zx = g1*db(ix,jy,kz)**2*( 0.224*qh + 0.776*qxw)*qh/chw
             ze =1.e18*zx*(6./(pi*1000.))**2
             dtmp(ix,kz) = dtmp(ix,kz) + ze
             dtmph = ze
             endif
             
            endif
             
        !     IF ( an(ix,jy,kz,lh) .gt. 1.0e-3 ) write(0,*)  'Graupel Z : ',dtmph,ze
           endif
          
          else
          
          dtmph = 0.0
          
           if ( an(ix,jy,kz,lh) .ge. qhmin ) then
             gtmp(ix,kz) = dadh*an(ix,jy,kz,lh)**(0.25)
             if ( gtmp(ix,kz) .gt. 0.0 ) then
             dtmph =  zhdryc*an(ix,jy,kz,lh)**2/gtmp(ix,kz)
             if ( temk(ix,jy,kz) .lt. tfr ) then
               dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhdryc*an(ix,jy,kz,lh)**2/gtmp(ix,kz)
             else
!               IF ( hwdn .gt. 700.0 ) THEN
                 dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhdryc*an(ix,jy,kz,lh)**2/gtmp(ix,kz)
! 
!     &                               (zhwetc*gtmp(ix,kz)**7)**0.95
!               ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + zhwetc*gtmp(ix,kz)**7
!               ENDIF
             endif
             endif
           endif
          
         
          
          endif
 

          endif ! }
          
          endif ! na .gt. 5

        
        if ( izieg .ge. 1 .and. lhl .gt. 1 ) then

        hldn = 900.0
        gtmp(ix,kz) = 0.0
        dtmphl = 0.0
        qxw = 0.0
        

        if ( lvhl .gt. 1 ) then
          if ( an(ix,jy,kz,lvhl) .gt. 0.0 ) then
            hldn = db(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
            hldn = Min( 900., Max( 300., hldn ) )
          else 
            hldn = 900. 
          endif
        else
          hldn = rho_qhl
        endif 


        if ( ipconc .ge. 5 ) then

           ltest = .false.
           if ( lzhl > 1 ) then
             if ( an(ix,jy,kz,lzhl) > 0.0 .and. an(ix,jy,kz,lhl) > qhlmin .and. &
                  an(ix,jy,kz,lnhl) > 0.0 ) ltest = .true.
           endif

          if ( ltest .or. ( an(ix,jy,kz,lhl) .ge. qhlmin .and. an(ix,jy,kz,lnhl) .gt. 0.) ) then !{
            chl = an(ix,jy,kz,lnhl)
            if ( chl .gt. 0.0 ) then !{
             xvhl = db(ix,jy,kz)*an(ix,jy,kz,lhl)/         &
     &        (hldn*Max(1.0e-9,an(ix,jy,kz,lnhl)))
            if ( xvhl .lt. xvhlmn .or. xvhl .gt. xvhlmx ) then ! {
              xvhl = Min( xvhlmx, Max( xvhlmn,xvhl ) )
              chl = db(ix,jy,kz)*an(ix,jy,kz,lhl)/(xvhl*hldn)
              ! do not update state in dbz calc. ! an(ix,jy,kz,lnhl) = chl
            endif ! }

             if ( lhlw .gt. 1 ) then
               if ( iusewethail .eq. 1 ) then
                  qxw = an(ix,jy,kz,lhlw)
               elseif ( iusewethail .eq. 2 ) then
                  if ( hldn .lt. 300. ) then
                    qxw = an(ix,jy,kz,lhlw)
                  endif
               endif
             endif
            
             if ( lzhl .gt. 1 ) then !{
              x = (0.224*an(ix,jy,kz,lhl) +  0.776*qxw)/an(ix,jy,kz,lhl)  ! weighted average of dielectric const
              dtmphl = 1.e18*x*an(ix,jy,kz,lzhl)*(hldn/rwdn)**2
              dtmp(ix,kz) = dtmp(ix,kz) + dtmphl
             else !}

             g1 = (6.0 + alphahl)*(5.0 + alphahl)*(4.0 + alphahl)/((3.0 + alphahl)*(2.0 + alphahl)*(1.0 + alphahl))
             zx = g1*db(ix,jy,kz)**2*( 0.224*an(ix,jy,kz,lhl) + 0.776*qxw)*an(ix,jy,kz,lhl)/chl
!             zx = g1*(db(ix,jy,kz)*an(ix,jy,kz,lhl))**2/chl
             ze = 1.e18*zx*(6./(pi*1000.))**2 ! 3/28/2016 removed extra factor of 0.224
             dtmp(ix,kz) = dtmp(ix,kz) + ze
             dtmphl = ze
             
             endif !}
            endif!}
        !     IF ( an(ix,jy,kz,lh) .gt. 1.0e-3 ) write(0,*)  'Graupel Z : ',dtmph,ze
           endif

          
          else
          
          
           if ( an(ix,jy,kz,lhl) .ge. qhlmin ) then ! {
            dadhl = ( db(ix,jy,kz) /(pi*hldn*cnohl) )**(.25)
             gtmp(ix,kz) = dadhl*an(ix,jy,kz,lhl)**(0.25)
             if ( gtmp(ix,kz) .gt. 0.0 ) then ! {

              zhldryc = 0.224*cr2*( db(ix,jy,kz)/rwdn)**2/cnohl 

             dtmphl =  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)

             if ( temk(ix,jy,kz) .lt. tfr ) then
               dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)
             else
!               IF ( hwdn .gt. 700.0 ) THEN
                 dtmp(ix,kz) = dtmp(ix,kz) +                   &
     &                  zhldryc*an(ix,jy,kz,lhl)**2/gtmp(ix,kz)
! 
!     :                               (zhwetc*gtmp(ix,kz)**7)**0.95
!               ELSE
!                 dtmp(ix,kz) = dtmp(ix,kz) + zhwetc*gtmp(ix,kz)**7
!               ENDIF
             endif
             endif ! }
           
           endif ! }
          
         endif ! ipconc .ge. 5


        endif ! izieg .ge. 1 .and. lhl .gt. 1 

          
           
          if ( dtmp(ix,kz) .gt. 0.0 ) then
            dbz(ix,jy,kz) = Max(dbzmin, 10.0*Log10(dtmp(ix,kz)) )
            
            if ( dbz(ix,jy,kz) .gt. dbzmax ) then
              dbzmax = Max(dbzmax,dbz(ix,jy,kz))
              imx = ix
              jmx = jy
              kmx = kz
            endif
          else 
             dbz(ix,jy,kz) = dbzmin
             if ( lh > 1 .and. lhl > 1) then
               if ( an(ix,jy,kz,lh) > 1.0e-3 ) then
                 write(0,*) 'radardbz: qr,qh,qhl = ',an(ix,jy,kz,lr), an(ix,jy,kz,lh),an(ix,jy,kz,lhl)
                 write(0,*) 'radardbz: dtmps,dtmph,dadh,dadhl,dtmphl = ',dtmps,dtmph,dadh,dadhl,dtmphl
                 
                 if ( lzh>1 .and. lzhl>1 ) write(0,*) 'radardbz: zh, zhl = ',an(ix,jy,kz,lzh),an(ix,jy,kz,lzhl)
               endif
             endif
          endif

!         IF ( an(ix,jy,kz,lh) .gt. 1.e-4 .and. 
!     &        dbz(ix,jy,kz) .le. 0.0 ) THEN
!          write(0,*) 'dbz = ',dbz(ix,jy,kz)
!          write(0,*) 'Hail intercept: ',xcnoh,ix,kz
!          write(0,*) 'Hail,snow q: ',an(ix,jy,kz,lh),an(ix,jy,kz,ls)
!          write(0,*) 'Hail,snow c: ',an(ix,jy,kz,lnh),an(ix,jy,kz,lns)
!          write(0,*) 'dtmps,dtmph = ',dtmps,dtmph
!         ENDIF
        if ( .not. dtmp(ix,kz) .lt. 1.e30 .or. dbz(ix,jy,kz) > 190.0 ) then
!        IF ( ix == 31 .and. kz == 20 .and. jy == 23 ) THEN
!          write(0,*) 'my_rank = ',my_rank
          write(0,*) 'ix,jy,kz = ',ix,jy,kz
          write(0,*) 'dbz = ',dbz(ix,jy,kz)
          write(0,*) 'db, zhdryc = ',db(ix,jy,kz),zhdryc
          write(0,*) 'Hail intercept: ',xcnoh,ix,kz
          write(0,*) 'Hail,snow q: ',an(ix,jy,kz,lh),an(ix,jy,kz,ls)
          write(0,*) 'graupel density hwdn = ',hwdn
          write(0,*) 'rain q: ',an(ix,jy,kz,lr)
          write(0,*) 'ice q: ',an(ix,jy,kz,li)
          if ( lhl .gt. 1 ) write(0,*) 'Hail (lhl): ',an(ix,jy,kz,lhl)
          if (ipconc .ge. 3 ) write(0,*) 'rain c: ',an(ix,jy,kz,lnr)
          if ( lzr > 1 ) write(0,*) 'rain Z: ',an(ix,jy,kz,lzr)
          if ( ipconc .ge. 5 ) then
          write(0,*) 'Hail,snow c: ',an(ix,jy,kz,lnh),an(ix,jy,kz,lns)
          if ( lhl .gt. 1 ) write(0,*) 'Hail (lnhl): ',an(ix,jy,kz,lnhl)
          if ( lzhl .gt. 1 ) then 
            write(0,*) 'Hail (lzhl): ',an(ix,jy,kz,lzhl)
            write(0,*) 'chl,xvhl,dhl = ',chl,xvhl,(xvhl*6./3.14159)**(1./3.)
            write(0,*) 'xvhlmn,xvhlmx = ',xvhlmn,xvhlmx
          endif
          endif
          write(0,*) 'chw,xvh = ', chw,xvh
          write(0,*) 'dtmps,dtmph,dadh,dadhl,dtmphl = ',dtmps,dtmph,dadh,dadhl,dtmphl
          write(0,*) 'dtmpr = ',dtmpr
          write(0,*) 'gtmp = ',gtmp(ix,kz),dtmp(ix,kz)
          if ( .not. (dbz(ix,jy,kz) .gt. -100 .and. dbz(ix,jy,kz) .lt. 200 ) ) then
            write(0,*) 'dbz out of bounds!'
!           STOP
          endif
         endif

           
          enddo ! ix
         enddo ! kz
      enddo ! jy
            
      
      
      
!      write(0,*)  'na,lr = ',na,lr
      if ( printyn .eq. 1 ) then
!      IF ( dbzmax .gt. dbzmin ) THEN
        write(iunit,*) 'maxdbz,ijk = ',dbzmax,imx,jmx,kmx
        write(iunit,*) 'qrw = ',an(imx,jmx,kmx,lr)
        
        if ( lh .gt. 1 ) then
          write(iunit,*) 'qi  = ',an(imx,jmx,kmx,li)
          write(iunit,*) 'qsw = ',an(imx,jmx,kmx,ls)
          write(iunit,*) 'qhw = ',an(imx,jmx,kmx,lh)
          if ( lhl .gt. 1 ) write(iunit,*) 'qhl = ',an(imx,jmx,kmx,lhl)
        endif

      
      endif
      
      
      return
      end subroutine radardd02
      

! ##############################################################################
! ##############################################################################


! #####################################################################
! #####################################################################
!
! Subroutine for explicit cloud condensation and droplet nucleation
!
! 11/30/2022: Fixed droplet evaporation heating term for 1 eqtset=2 (was only doing eqtset=1)
!
   subroutine NUCOND    &
     &  (nx,ny,nz,na,jyslab & 
     &  ,nor,norz,dtp,nxi & 
     &  ,dz3d & 
     &  ,t0,t9 & 
     &  ,an,dn,p2 & 
     &  ,pn,w & 
     &  ,ngs   &
     &  ,axtra,io_flag &
     &  ,ssfilt,t00,t77,flag_qndrop  &
     & )


   implicit none

!      real :: cwmasn = 1000.*0.523599*(2.*2.e-6)**3 
      integer :: nx,ny,nz,na,nxi
      integer :: nor,norz, jyslab ! ,nht,ngt,igsr
      real    :: dtp  ! time step
      logical :: flag_qndrop

      integer, parameter :: ng1 = 1


!
! external temporary arrays
!
      real t00(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t77(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real t0(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t1(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t3(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t4(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t5(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t6(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t7(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real t8(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real t9(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      

      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)  ! perturbation Pi
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
!      real qv(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      real ssfilt(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)
      

      real pb(-norz+ng1:nz+norz)
      real pinit(-norz+ng1:nz+norz)

      real dz3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

      
    ! local


      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)
      logical :: io_flag
      
      real :: dv
      real :: ccnefactwo, sstmp, cn1, cnuctmp

! 
!  declarations microphysics and for gather/scatter
!
      real, parameter :: cwmas30 = 1000.*0.523599*(2.*30.e-6)**3 ! mass of 30-micron radius droplet, for sat. adj.
      real, parameter :: cwmas20 = 1000.*0.523599*(2.*20.e-6)**3 ! mass of 20-micron radius droplet, for sat. adj.
      integer nxmpb,nzmpb,nxz
      integer mgs,ngs,numgs,inumgs
      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs)
      integer nsvcnt
      
      integer ix,kz,i,n, kp1, km1
      integer :: jy, jgs
      integer ixb,ixe,jyb,jye,kzb,kze
    
      integer itile,jtile,ktile
      integer ixend,jyend,kzend,kzbeg
      integer nxend,nyend,nzend,nzbeg

!
! Variables for Ziegler warm rain microphysics
!      


      real ccnc(ngs), ccna(ngs), cnuc(ngs), cwnccn(ngs)
      real :: ccnc_nu(ngs), ccnc_ac(ngs), ccnc_co(ngs)
      real ccncuf(ngs)
      real sscb  ! 'cloud base' SS threshold
      parameter ( sscb = 2.0 )
      integer idecss  ! flag to turn on (=1) decay of ssmax when no cloud or ice crystals
      parameter ( idecss = 1 )
      integer iba ! flag to do condensation/nucleation in 1st or 2nd loop
                  ! =0 to use ad to calculate SS
                  ! =1 to use an at end of main jy loop to calculate SS
      parameter (iba = 1)
      integer ifilt   ! =1 to filter ssat, =0 to set ssfilt=ssat
      parameter ( ifilt = 0 ) 
      real temp1,temp2 ! ,ssold
      real :: ssmax(ngs)      ! maximum SS experienced by a parcel
      real ssmx
      real dnnet,dqnet
!      real cnu,rnu,snu,cinu
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
      real ventrx(ngs)
      real ventrxn(ngs)
      real volb, t2s
      real, parameter :: aa1 = 9.44e15, aa2 = 5.78e3  ! a1 in Ziegler

      real ec0, ex1, ft, rhoinv(ngs)
      
      real chw, g1, rd1

      real ac1,bc, taus, c1,d1,e1,f1,p380,tmp,tmp2 ! , sstdy, super
      real tmpmx, fw, qctmp
      real x,y,del,r,alpr
      double precision :: vent1,vent2
      real g1palp
      real bs
      real v1, v2
      real d1r, d1i, d1s, e1i
      integer nc ! condensation step
      real dtcon,dtcon1,dtcon2 ! condensation time step (dtcon*nc = dtp)
      real delta
      integer ltemq1,ltemq1m ! ,ltemq1m2
      real dqv,qv1,ss1,ss2,qvs1,dqvs,dtemp,dt1   ! temporaries for condensation

      real ssi1, ssi2, dqvi, dqvis, dqvii,qis1
      real dqvr, dqc, dqr, dqi, dqs
      real qv1m,qvs1m,ss1m,ssi1m,qis1m
      real cwmastmp 
      real  dcloud,dcloud2 ! ,as, bs
      real dcrit
      real cn(ngs), cnuf(ngs)
      real :: ccwmax

      integer ltemq
      
      integer il

      real  es(ngs) ! ss(ngs),
!      real  eis(ngs)
      real ssf(ngs),ssfkp1(ngs),ssfkm1(ngs),ssat0(ngs)
      real, parameter :: ssfcut = 4.0
      real ssfjp1(ngs),ssfjm1(ngs)
      real ssfip1(ngs),ssfim1(ngs)

      real supcb, supmx
      parameter (supcb=0.5,supmx=238.0)
      real r2dxm, r2dym, r2dzm
      real dssdz, dssdy, dssdx
!      real tqvcon
      real epsi,d
      parameter (epsi = 0.622, d = 0.266)
      real r1,qevap ! ,slv
      
      real vr,nrx,qr,z1,z2,rdi,alp,xnutmp,xnuc
      real ctmp, ccwtmp
      real f5, qvs0  ! Kessler condensation factor
      real    :: t0p1, t0p3
      real qvex
      
!      real, dimension(ngs) :: temp, tempc, elv, elf, els, pqs, theta, temg, temcg
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs)
      real temp(ngs),tempc(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) ! ,tembzg(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real felv(ngs),felf(ngs),fels(ngs)
      real felvcp(ngs),felvpi(ngs)
      real gamw(ngs),gams(ngs)   !   qciavl(ngs),
      real tsqr(ngs),ssi(ngs),ssw(ngs)
      real cc3(ngs),cqv1(ngs),cqv2(ngs)
      real qcwtmp(ngs),qtmp

      real fvent(ngs) !,fraci(ngs),fracl(ngs)
      real fwvdf(ngs),ftka(ngs),fthdf(ngs)
      real fadvisc(ngs),fakvisc(ngs)
      real fci(ngs),fcw(ngs)
      real fschm(ngs),fpndl(ngs)

      real pres(ngs),pipert(ngs)
      real pk(ngs)
      real rho0(ngs),pi0(ngs)
      real rhovt(ngs)
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
      real thsave(ngs)
      real qss0(ngs)
      real fcqv1(ngs)
      real wvel(ngs),wvelkm1(ngs)

      real wvdf(ngs),tka(ngs)
      real advisc(ngs)

      real rwvent(ngs)
      

      real :: qx(ngs,lv:lhab)
      real :: cx(ngs,lc:lhab)
      real :: xv(ngs,lc:lhab)
      real :: xmas(ngs,lc:lhab)
      real :: xdn(ngs,lc:lhab)
      real :: xdia(ngs,lc:lhab,3)
      real :: alpha(ngs,lc:lhab)
      real :: zx(ngs,lr:lhab)


      logical zerocx(lc:lqmx)
      
      logical :: lprint

      integer, parameter :: iunit = 0
      
      real :: frac, hwdn, tmpg
      
      real :: cvm,cpm,rmm

      real, parameter ::      cpv = 1885.0       ! specific heat of water vapor at constant pressure
      
      integer :: kstag
      
      integer :: count
      
! -------------------------------------------------------------------------------
      itile = nxi
      jtile = ny
      ktile = nz
      ixend = nxi
      jyend = ny
      kzend = nz
      nxend = nxi + 1
      nyend = ny + 1
      nzend = nz
      kzbeg = 1
      nzbeg = 1

      if ( ac_opt > 0 )  ccnefactwo =  (1.63e-3/(cck * beta(3./2., cck/2.)))**(1.0/(cck + 2.0))
      f5 = 237.3 * 17.27 * 2.5e6 / cp ! combined constants for rain condensation (Soong and Ogura 73)

       jy = 1
       kstag = 0
       pb(:) = 0.0
       pinit(:) = 0.0
      
      if ( ipconc <= 1 .or. isedonly == 2 ) GOTO 2200

!
!  Ziegler nucleation 
!

!      ssfilt(:,:,:) = 0.0
      ssmx = 0
      count = 0

      do kz = 1,nz-kstag
        do ix = 1,nxi

         temp1 = an(ix,jy,kz,lt)*t77(ix,jy,kz)
          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          c1 = t00(ix,jy,kz)*tabqvs(ltemq)

          if ( c1 > 0. ) then
            ssfilt(ix,jy,kz) = 100.*(an(ix,jy,kz,lv)/c1 - 1.0)  ! from "new" values
          endif

        enddo
      enddo


!
!     jy = 1 ! working on a 2d slab
!!  VERY IMPORTANT:  SET jgs = jy

      jgs = jy

!
!..Gather microphysics
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: Gather stage'

      nxmpb = 1
      nzmpb = 1
      nxz = nxi*nz
      numgs = nxz/ngs + 1


      do 2000 inumgs = 1,numgs

      ngscnt = 0


      kzb = nzmpb
      kze = nz-kstag
 !     if (kzbeg .le. nzmpb .and. kzend .gt. nzmpb) kzb = nzmpb

      ixb = nxmpb
      ixe = itile

      do kz = kzb,kze
      do ix = nxmpb,nxi

      pqs(1) = 380.0/(pn(ix,jy,kz) + pb(kz))
      theta(1) = an(ix,jy,kz,lt)
      temg(1) = t0(ix,jy,kz)

      temcg(1) = temg(1) - tfr
      ltemq = (temg(1)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvs(1) = pqs(1)*tabqvs(ltemq)
      qis(1) = pqs(1)*tabqis(ltemq)

      qss(1) = qvs(1)


      if ( temg(1) .lt. tfr ) then
      end if
!
      if ( (temg(1) .gt. tfrh .or. an(ix,jy,kz,lv)/qvs(1) > maxlowtempss ) .and.  &
     &   ( an(ix,jy,kz,lv)  .gt. qss(1) .or. &
     &     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or.  &
     &     ( an(ix,jy,kz,lr)  .gt. qxmin(lr) .and. rcond == 2 )  &
     &     )) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 2100
      end if

      end do  !ix

      nxmpb = 1
      end do  !kz
!      if ( jy .eq. (ny-jstag) ) iend = 1
 2100 continue

      if ( ngscnt .eq. 0 ) go to 29998

      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: dbg = 8'
      
!      write(0,*) 'NUCOND: dbg = 8, ngscnt,ssmx = ',ngscnt,ssmx

      
      qx(:,:) = 0.0
      cx(:,:) = 0.0
      zx(:,:) = 0.0

      xv(:,:) = 0.0
      xmas(:,:) = 0.0

      if ( imurain == 1 ) then
        alpha(:,lr) = alphar
      elseif ( imurain == 3 ) then
        alpha(:,lr) = xnu(lr)
      endif

!
!  define temporaries for state variables to be used in calculations
!
      do mgs = 1,ngscnt
      qx(mgs,lv) = an(igs(mgs),jy,kgs(mgs),lv)
       do il = lc,lhab
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0)
       enddo

       qcwtmp(mgs) = qx(mgs,lc)


      theta0(mgs) = an(igs(mgs),jy,kgs(mgs),lt) !
      thetap(mgs) = 0.0
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      qv0(mgs) =  qx(mgs,lv)
      qwvp(mgs) = qx(mgs,lv) - qv0(mgs)

       pres(mgs) = pn(igs(mgs),jy,kgs(mgs)) + pb(kgs(mgs))
       pipert(mgs) = p2(igs(mgs),jy,kgs(mgs))
       rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
       rhoinv(mgs) = 1.0/rho0(mgs)
       rhovt(mgs) = Sqrt(rho00/rho0(mgs))
       pi0(mgs) = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs))
       temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
!       pk(mgs) = t77(igs(mgs),jy,kgs(mgs)) ! ( pres(mgs) / poo ) ** cap
       pk(mgs)   = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs)) ! t77(igs(mgs),jy,kgs(mgs))
       temcg(mgs) = temg(mgs) - tfr
       qss0(mgs) = (380.0)/(pres(mgs))
       pqs(mgs) = (380.0)/(pres(mgs))
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
       qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       qis(mgs) = pqs(mgs)*tabqis(ltemq)
!
        qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
        es(mgs) = 6.1078e2*tabqvs(ltemq)
        qss(mgs) = qvs(mgs)


        temgx(mgs) = min(temg(mgs),313.15)
        temgx(mgs) = max(temgx(mgs),233.15)
        felv(mgs) = 2500837.367 * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
!
        if ( eqtset <= 1 ) then
          felvcp(mgs) = felv(mgs)*cpi
        else ! equation set 2 in cm1
          tmp = qx(mgs,li)+qx(mgs,ls)+qx(mgs,lh)
          if ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          if ( lf > 1 ) tmp = tmp + qx(mgs,lf)
          cvm = cv+cvv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)
          cpm = cp+cpv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)
          rmm=rd+rw*qx(mgs,lv)
          
          if ( eqtset == 2 ) then

           felvcp(mgs) = (felv(mgs)-rw*temg(mgs))/cvm

          else
            felvcp(mgs) = (felv(mgs)*cv/(cp) - rw*temg(mgs)*(1.0-rovcp*cpm/rmm))/cvm
            felvpi(mgs) = pi0(mgs)*rovcp*(felv(mgs)/(temg(mgs)) - rw*cpm/rmm)/cvm
          endif

        endif

        temcgx(mgs) = min(temg(mgs),273.15)
        temcgx(mgs) = max(temcgx(mgs),223.15)
        temcgx(mgs) = temcgx(mgs)-273.15
        felf(mgs) = 333690.6098 + (2030.61425)*temcgx(mgs) - (10.46708312)*temcgx(mgs)**2
!
        fels(mgs) = felv(mgs) + felf(mgs)
        fcqv1(mgs) = 4098.0258*felv(mgs)*cpi

      wvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)* &
     &  (101325.0/(pb(kgs(mgs)) + pn(igs(mgs),jgs,kgs(mgs))))                            ! diffusivity of water vapor, Hall and Pruppacher (76)
      advisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* &
     &  (temg(mgs)/296.0)**(1.5)                         ! dynamic viscosity (SMT; see Beard & Pruppacher 71)
      tka(mgs) = tka0*advisc(mgs)/advisc1                 ! thermal conductivity


      enddo



!
! load concentrations
!
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
        cwnccn(mgs) = cwccn*rho0(mgs)/rho00 ! background ccn count
        cn(mgs) = 0.0
        if ( lss > 1 ) then 
          ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        else
          ssmax(mgs) = 0.0
        endif
        if ( lccn .gt. 1 .and. ac_opt == 0 ) then
          if ( lccnuf .gt. 1 .and. i_uf_or_ccn > 0 ) then
             ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn) + an(igs(mgs),jy,kgs(mgs),lccnuf)
          else
             ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
          endif
        else
          ccnc(mgs) = cwnccn(mgs)
        endif
        if ( lccnuf .gt. 1 .and. i_uf_or_ccn == 0 ) then
          ccncuf(mgs) = an(igs(mgs),jy,kgs(mgs),lccnuf)
        else
          ccncuf(mgs) = 0.0
        endif
        cnuf(mgs) = 0.0
        if ( lccna > 1 ) then
          ccna(mgs) = an(igs(mgs),jy,kgs(mgs),lccna) ! predicted count of activated ccn
        else
          if ( lccn > 1 ) then
            ccna(mgs) = cwnccn(mgs) - ccnc(mgs) ! diagnose activated ccn as background value - remaining unactivated ccn
          else
            ccna(mgs) = cx(mgs,lc) ! approximation of number of activated ccn
          endif
        endif
       end do
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
       end do
      end if

!        cnuc(1:ngscnt) = cwccn*rho0(mgs)/rho00*(1. - renucfrac) + ccnc(1:ngscnt)*renucfrac
       do mgs = 1,ngscnt
       ! default value of renucfrac is 0.0
        if ( irenuc /= 6 ) then
        cnuc(mgs) = Max(ccnc(mgs),cwnccn(mgs))*(1. - renucfrac) + ccnc(mgs)*renucfrac
        else
        cnuc(mgs) = Max(ccnc(mgs),cwnccn(mgs))*(1. - renucfrac) + Max(0.0,ccnc(mgs) - ccna(mgs))*renucfrac
        endif
        if ( renucfrac >= 0.999 ) then
          if ( temg(mgs) < 265. ) then
            if ( qx(mgs,lc) > 10.*qxmin(lc) .and. w(igs(mgs),jgs,kgs(mgs)) > 2.0 ) then
             cnuc(mgs) = 0.0 !  Min(cnuc(mgs), 0.5*cx(mgs,lc) ) ! Hack to reduce nucleation at low temp in updraft when ccn are not predicted
            else
             cnuc(mgs) = 0.1*cnuc(mgs)
            endif
          endif
        endif
       enddo

!  Set density
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: Set density'

      do mgs = 1,ngscnt
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
      end do

      ventrx(:) = ventr
      ventrxn(:) = ventrn
      

!  Find shape parameter rain

      if ( lzr > 1 .and. rcond == 2 ) then ! { RAIN SHAPE PARAM
      do mgs = 1,ngscnt
         zx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lzr), 0.0)
      enddo

!      CALL cld_cpu('Z-MOMENT-1r2')
          il = lr
          do mgs = 1,ngscnt

         if ( zx(mgs,il) <= zxmin ) then
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         elseif ( cx(mgs,il) <= 0.0 ) then
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif

         if ( qx(mgs,lr) .gt. qxmin(lr) ) then

          xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
          if ( xv(mgs,lr) .gt. xvmx(lr) ) then
            xv(mgs,lr) = xvmx(lr)
            cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
          elseif ( xv(mgs,lr) .lt. xvmn(lr) ) then
            xv(mgs,lr) = xvmn(lr)
            cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
          endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            if ( imurain == 3 ) then
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z1   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
            else
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
            z1   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
            
            endif
!            an(igs(mgs),jgs,kgs(mgs),ln(il)) = zx(mgs,il)
           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            if ( imurain == 3 ) then
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(chw*1000.*1000)
            else
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(chw*1000.*1000)
            
            endif

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
              if ( imurain == 3 ) then
               g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z1   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z1*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              elseif ( imurain == 1 ) then
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z1   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(z1*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
              
              endif
          endif
        
          if ( zx(mgs,lr) > 0.0 ) then
            vr = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
!            z1 = 36.*(alpha(kz)+2.0)*a(ix,jy,kz,lnr)*vr**2/((alpha(kz)+1.0)*pi**2)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z1 = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z1*(pi/6.*1000.)**2/xv


! determine shape parameter alpha by iteration
        if ( z1 .gt. 0.0 ) then

          if ( imurain == 3 ) then
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
!           write(0,*) 'kz, alp, alpha(kz) = ',kz,alp,alpha(kz),rd,z1,xv
           do i = 1,20
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z1*pi**2) - 1.
!           write(0,*) 'i,alp = ',i,alp
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo

         else ! imurain == 1
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)

            rd1 = z1*(pi/6.*xdn(mgs,il))**2*nrx/(rho0(mgs)*qr)**2

           alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd1) - 1.0

           do i = 1,10
            if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )

             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/ &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rd1) - 1.0

             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo

         
         endif
!         ENDIF

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
          if ( imurain == 3 ) then
           if ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) then

             if ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
               g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z1*(1./(xdn(mgs,il)))**2
               an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
             elseif ( rescale_low_alphar .and. alp <= rnumin ) then

              z1  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
              zx(mgs,il) = z1
             endif
           endif
           
          elseif ( imurain == 1 ) then
          
             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))

           if ( (rescale_low_alpha .or. rescale_high_alpha ) .and.  &
     &          ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) then



            if ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*rho0(mgs)**2*(qr)*qr/zx(mgs,lr)*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( rescale_low_alpha .and. alp <= alphamin ) then ! alpha = alphamin, so reset Z to prevent growth in C
             z1 = g1*rho0(mgs)**2*(qr)*qr/nrx
             z2  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z2
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z2
            endif
          endif ! imurain

          endif ! z > 0

           tmp = alpha(mgs,lr) + 4./3.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,lr) + 1.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma(alpha(mgs,lr) + 1.)
           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))

           if ( imurain == 3 .and. izwisventr == 2 ) then

           tmp = alpha(mgs,lr) + 1.5 + br/6.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 1.5 + br/6.)/Gamma(alpha(mgs,lr) + 1.)
           ventrxn(mgs) = x/(y*(alpha(mgs,lr) + 1.)**((1.+br)/6. + 1./3.))
           
           elseif ( imurain == 1 .and.  iferwisventr == 2 ) then

           tmp = alpha(mgs,lr) + 2.5 + br/2.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma(alpha(mgs,lr) + 1.5 + br/6.)/Gamma(alpha(mgs,lr) + 1.)
           ventrxn(mgs) = x/y
           
           
           endif

           
           endif
          endif
          
          endif
          
          enddo
!        CALL cld_cpu('Z-MOMENT-1r2')  
        endif ! }


!       write(0,*) 'NUCOND: Set ssf variables, ssmxinit =',ssmxinit
      ssmx = 0.0
      do mgs = 1,ngscnt
      
      kp1 = Min(nz, kgs(mgs)+1 )
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kp1) & 
     &                  +w(igs(mgs),jgs,kgs(mgs)))
      wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs)) & 
     &                  +w(igs(mgs),jgs,Max(1,kgs(mgs)-1)))

      ssat0(mgs)  = ssfilt(igs(mgs),jgs,kgs(mgs))
      ssf(mgs)    = ssfilt(igs(mgs),jgs,kgs(mgs))
!      ssmx = Max( ssmx, ssf(mgs) )

      
      ssfkp1(mgs) = ssfilt(igs(mgs),jgs,Min(nz-1,kgs(mgs)+1))
      ssfkm1(mgs) = ssfilt(igs(mgs),jgs,Max(1,kgs(mgs)-1))

!        IF ( wvel(mgs) /= 0.0 ) write(0,*) 'mgs,wvel1,ssf = ',mgs,wvel(mgs),ssf(mgs)


      enddo



!
!  cloud water variables
!

      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_DR: Set cloud water variables'

      do mgs = 1,ngscnt
      xv(mgs,lc) = 0.0
      if ( ipconc .ge. 2 .and. cx(mgs,lc) .gt. 1.0e6 ) then
        xmas(mgs,lc) = &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)
      else
       if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. cxmin ) then
        xmas(mgs,lc) = &
     &     min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),xdn(mgs,lc)*xvmn(lc)), &
     &      xdn(mgs,lc)*xvmx(lc) )

        cx(mgs,lc) = qx(mgs,lc)*rho0(mgs)/xmas(mgs,lc)

       elseif ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .le. cxmin ) then
!        xmas(mgs,lc) = xdn(mgs,lc)*4.*pi/3.*(5.0e-6)**3
!        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
        cx(mgs,lc) = Max( cxmin, rho0(mgs)*qx(mgs,lc)/cwmasx )
        xmas(mgs,lc) =  &
     &    min( max(qx(mgs,lc)*rho0(mgs)/cx(mgs,lc),cwmasn),cwmasx )
        xv(mgs,lc) = xmas(mgs,lc)/xdn(mgs,lc)

       else
        xmas(mgs,lc) = cwmasn
       endif
      endif
      xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3


      end do
!
! rain
!
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then

      if ( ipconc .ge. 3 ) then
        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-9,cx(mgs,lr)))
!      parameter( xvmn(lr)=2.8866e-13, xvmx(lr)=4.1887e-9 )  ! mks
        if ( xv(mgs,lr) .gt. xvmx(lr) ) then
          xv(mgs,lr) = xvmx(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        elseif ( xv(mgs,lr) .lt. xvmn(lr) ) then
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        endif

        xmas(mgs,lr) = xv(mgs,lr)*xdn(mgs,lr)
        xdia(mgs,lr,3) = (xmas(mgs,lr)*cwc1)**(1./3.) ! xdia(mgs,lr,1)
        if ( imurain == 3 ) then
!          xdia(mgs,lr,1) = (6.*pii*xv(mgs,lr)/(alpha(mgs,lr)+1.))**(1./3.)
          xdia(mgs,lr,1) = xdia(mgs,lr,3) ! formulae for Ziegler (1985) use mean volume diameter, not lambda**(-1)
        else ! imurain == 1, Characteristic diameter (1/lambda)
          xdia(mgs,lr,1) = (6.*piinv*xv(mgs,lr)/((alpha(mgs,lr)+3.)*(alpha(mgs,lr)+2.)*(alpha(mgs,lr)+1.)))**(1./3.)
        endif
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)

! Inverse exponential version:
!        xdia(mgs,lr,1) =
!     >  (qx(mgs,lr)*rho0(mgs)
!     > /(pi*xdn(mgs,lr)*cx(mgs,lr)))**(0.333333)
      else
        xdia(mgs,lr,1) = &
     &  (qx(mgs,lr)*rho0(mgs)/(pi*xdn(mgs,lr)*cno(lr)))**(0.25)
      end if
      else
        xdia(mgs,lr,1) = 1.e-9
!        rwrad(mgs) = 0.5*xdia(mgs,lr,1)
      end if

      end do


!
!  Ventilation coefficients

      do mgs = 1,ngscnt


      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))* & 
     &  (temg(mgs)/296.0)**(1.5)

      fakvisc(mgs) = fadvisc(mgs)*rhoinv(mgs)

      fwvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)* & 
     &  (101325.0/(pres(mgs)))
      
      fschm(mgs) = (fakvisc(mgs)/fwvdf(mgs))

      fvent(mgs) = (fschm(mgs)**(1./3.)) * (fakvisc(mgs)**(-0.5))

      end do
!
!
!  Ziegler nucleation 
!
!
! cloud evaporation, condensation, and nucleation
!  sqsat -> qss(mgs)

      do mgs=1,ngscnt
        dcloud = 0.0
        ! Skip points at low temperature if SS stays less than 1.08, 
        ! otherwise allow nucleation at low temp (will freeze at next time step)
        if ( temg(mgs) .le. tfrh .and. qx(mgs,lv)/qvs(mgs) < maxlowtempss ) then 
         cycle
        endif

      if( ssat0(mgs) .GT. 0. .OR. ssf(mgs) .GT. 0. ) go to 620
!6/4      IF( qvap(mgs) .EQ. qss(mgs) ) GO TO 631
!
!.... EVAPORATION. QV IS LESS THAN qss(mgs).
!.... EVAPORATE CLOUD FIRST
!
      if ( qx(mgs,lc) .LE. 0. ) go to 631
!.... CLOUD EVAPORATION.
! convert input 'cp' to cgs
      R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/ &
     &            (cp*(temg(mgs) - cbw)**2))
      QEVAP= Min( qx(mgs,lc), R1*(qss(mgs)-qvap(mgs)) )


      if ( qx(mgs,lc) <= QEVAP ) then ! go to 63
        qwvp(mgs) = qwvp(mgs) + qx(mgs,lc)
        thetap(mgs) = thetap(mgs) - felvcp(mgs)*qx(mgs,lc)/(pi0(mgs))
        if ( io_flag .and. nxtra > 1 ) then
           axtra(igs(mgs),jy,kgs(mgs),1) = -qx(mgs,lc)/dtp
        endif
        qx(mgs,lc) = 0.
        if ( restoreccn ) then
           if ( lccna > 1 ) then
              ccna(mgs) = ccna(mgs) - restoreccnfrac*cx(mgs,lc)
           elseif ( irenuc <= 2 ) then
              if ( .not. invertccn ) then
               ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + restoreccnfrac*cx(mgs,lc) ) )
              else
               ccnc(mgs) = ccnc(mgs) + restoreccnfrac*cx(mgs,lc)
              endif
          endif
        endif
        cx(mgs,lc) = 0.
      else
        qctmp = qx(mgs,lc)
        qwvp(mgs) = qwvp(mgs) + QEVAP
        qx(mgs,lc) = qx(mgs,lc) - QEVAP
        if ( qx(mgs,lc) .le. 0. ) then
          if ( restoreccn ) then
            if ( lccna > 1 ) then
              ccna(mgs) = ccna(mgs) - restoreccnfrac*cx(mgs,lc)
            elseif ( irenuc <= 2 ) then
!              ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + cx(mgs,lc) ) )
!              ccnc(mgs) = ccnc(mgs) + cx(mgs,lc)
              if ( .not. invertccn ) then
               ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + restoreccnfrac*cx(mgs,lc) ) )
              else
               ccnc(mgs) = ccnc(mgs) + restoreccnfrac*cx(mgs,lc)
              endif
            endif
          endif
          cx(mgs,lc) = 0.
        else
          tmp = 0.9*QEVAP*cx(mgs,lc)/qctmp ! let droplets get smaller but also remove some. A factor of 1.0 would maintain same size
          if ( restoreccn ) then
            if ( lccna > 1 ) then
              ccna(mgs) = ccna(mgs) - restoreccnfrac*tmp
            elseif ( irenuc <= 2 ) then
 !             ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + tmp ) )
!              ccnc(mgs) = ccnc(mgs) + tmp
              if ( .not. invertccn ) then
               ccnc(mgs) = Max( ccnc(mgs), Min( qccn*rho0(mgs), ccnc(mgs) + restoreccnfrac*tmp ) )
              else
               ccnc(mgs) = ccnc(mgs) + restoreccnfrac*tmp
              endif
            endif
          endif
          cx(mgs,lc) = cx(mgs,lc) - tmp
        endif
        thetap(mgs) = thetap(mgs) - felvcp(mgs)*QEVAP/(pi0(mgs))
        if ( io_flag .and. nxtra > 1 ) then
           axtra(igs(mgs),jy,kgs(mgs),1) = -QEVAP/dtp
        endif

      endif

      go to 631


  620 continue

!.... CLOUD CONDENSATION

        if ( qx(mgs,lc) .GT. qxmin(lc) .and. cx(mgs,lc) .ge. 1. ) then



!       ac1 =  xdn(mgs,lc)*elv(kgs(mgs))**2*epsi/
!     :        (tka(kgs(mgs))*rw*temg(mgs)**2)
! took out xdn factor because it cancels later...
       ac1 =  felv(mgs)**2/(tka(mgs)*rw*temg(mgs)**2)


!       bc = xdn(mgs,lc)*rw*temg(mgs)/
!     :       (epsi*wvdf(kgs(mgs))*es(mgs))
! took out xdn factor because it cancels later...
       bc =   rw*temg(mgs)/(wvdf(mgs)*es(mgs))

!       bs = rho0(mgs)*((rd*temg(mgs)/(epsi*es(mgs)))+
!     :             (epsi*elv(kgs(mgs))**2/(pres(mgs)*temg(mgs)*cp)))

!       taus = Min(dtp, xdn(mgs,lc)*rho0(mgs)*(ac1+bc)/
!     :        (4*pi*0.89298*BS*0.5*xdia(mgs,lc,1)*cx(mgs,lc)*xdn(mgs,lc)))

!
      if ( ssf(mgs) .gt. 0.0 .or. ssat0(mgs) .gt. 0.0 ) then
       if ( ny .le. 2 ) then
!        write(0,*)  'undershoot: ',ssf(mgs),
!     :   ( (qx(mgs,lv) - dcloud)/c1 - 1.0)*100.
       endif



       if ( qx(mgs,lc) .gt. qxmin(lc) ) then

         if ( xdia(mgs,lc,1) .le. 0.0 ) then
          xmas(mgs,lc) = cwmasn
          xdia(mgs,lc,1) = (xmas(mgs,lc)*cwc1)**c1f3
         endif
        d1 = (1./(ac1 + bc))*4.0*pi*ventc &
     &        *0.5*xdia(mgs,lc,1)*cx(mgs,lc)*rhoinv(mgs)

       else
         d1 = 0.0
       endif

       if ( rcond .eq. 2 .and. qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > 1.e-9 ) then
          if ( imurain == 3 ) then
           if ( izwisventr == 1 ) then
            rwvent(mgs) = ventrx(mgs)*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)
           else ! izwisventr = 2
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
          rwvent(mgs) =   &
     &  (0.78*ventrx(mgs) + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )
           endif

          else ! imurain == 1

           if ( iferwisventr == 1 ) then
             alpr = Min(alpharmax,alpha(mgs,lr) )
!             alpr = alpha(mgs,lr)
             x =  1. + alpr

              tmp = 1 + alpr
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = 2.5 + alpr + 0.5*bx(lr)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions

!         vent1 = dble(xdia(mgs,lr,1))**(-2. - alpr) ! Actually OK
!         vent2 = dble(1./xdia(mgs,lr,1) + 0.5*fx(lr))**dble(2.5+alpr+0.5*bx(lr))  ! Actually OK
         vent1 = dble(xdia(mgs,lr,1))**(0.5 + 0.5*bx(lr)) ! 2016.2.26 Changed for consistency with derivation (recast formula)
         vent2 = dble(1. + 0.5*fx(lr)*xdia(mgs,lr,1))**dble(2.5+alpr+0.5*bx(lr))
        
        
        rwvent(mgs) =    &
     &    0.78*x +    &
     &    0.308*fvent(mgs)*y*   &
     &            Sqrt(ax(lr)*rhovt(mgs))*(vent1/vent2)

           elseif ( iferwisventr == 2 ) then
          
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
            x =  1. + alpha(mgs,lr)

            rwvent(mgs) =   &
     &        (0.78*x + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &         *Sqrt((ar*rhovt(mgs)))   &
     &         *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )

          
          endif ! iferwisventr
          
       endif ! imurain

       d1r = (1./(ac1 + bc))*4.0*pi*rwvent(mgs) & 
     &        *0.5*xdia(mgs,lr,1)*cx(mgs,lr)*rhoinv(mgs)
       else
       d1r = 0.0
       endif
       
       
       e1  = felvcp(mgs)/(pi0(mgs))
       f1 = pk(mgs) ! (pres(mgs)/poo)**cap

!
!  fifth trial to see what happens:
!
       ltemq = (temg(mgs)-163.15)/fqsat+1.5
       ltemq = Min( nqsat, Max(1,ltemq) )
       ltemq1 = ltemq
       temp1 = temg(mgs)
       p380 = 380.0/pres(mgs)

!       taus = Max( 0.05*dtp, Min(taus, 0.25*dtp ) )
!       nc = NInt(dtp/Min(1.0,0.5*taus))
!       dtcon = dtp/float(nc)
       ss1 = qx(mgs,lv)/qvs(mgs)
       ss2 = ss1
       temp2 = temp1
       qv1 = qx(mgs,lv)
       qvs1 = qvs(mgs)
       qis1 = qis(mgs)
       dt1 = 0.0


!          dtcon = Max(dtcon,0.2)
!          nc = Nint(dtp/dtcon)

       ltemq1 = ltemq
! want to start out with a small time step to handle the steep slope
! and fast changes, then can switch to a larger step (dtcon2) for the
! rest of the big time step.
! base the initial time step (dtcon1) on the slope (delta)
       if ( Abs(ss1 - 1.0) .gt. 1.e-5 ) then
         delta = 0.5*(qv1-qvs1)/(d1*(ss1 - 1.0))
       else
         delta = 0.1*dtp
       endif
! delta is the extrapolated time to get halfway from qv1 to qvs1
! want at least 5 time steps to the halfway point, so multiply by 0.2
! for the initial time step
       dtcon1 = Min(0.05,0.2*delta)
       nc = Max(5,2*NInt( (dtp-4.0*dtcon1)/delta))
       dtcon2 = (dtp-4.0*dtcon1)/nc

       n = 1
       dt1 = 0.0
       nc = 0
       dqc = 0.0
       dqr = 0.0
       dqi = 0.0
       dqs = 0.0
       dqvii = 0.0
       dqvis = 0.0

       RK2c: do WHILE ( dt1 .lt. dtp )
          nc = 0
          if ( n .le. 4 ) then
            dtcon = dtcon1
          else
            dtcon = dtcon2
          endif
 609       dqv  = -(ss1 - 1.)*d1*dtcon
           dqvr = -(ss1 - 1.)*d1r*dtcon
            dtemp = -0.5*e1*f1*(dqv + dqvr)
!          write(0,*) 'RK2c dqv1 = ',dqv
! calculate midpoint values:
     !      ltemq1m = ltemq1 + Nint(dtemp*fqsat + 0.5)

         ! 7.6.2016: Test full calc of ltemq
           ltemq1m = (temp1+dtemp-163.15)*fqsati+1.5
           ltemq1m = Min( nqsat, Max(1,ltemq1m) )

           if ( ltemq1m .lt. 1 .or. ltemq1m .gt. nqsat ) then
             write(0,*) 'stop in nucond line 1192 '
             write(0,*) ' ltemq1m,icond = ',ltemq1m,icond
             write(0,*) ' dtemp,e1,f1,dqv,dqvr = ', dtemp,e1,f1,dqv,dqvr
             write(0,*) ' d1,d1r,dtcon,ss1 = ',d1,d1r,dtcon,ss1
             write(0,*) ' dqc, dqr = ',dqc,dqr
             write(0,*) ' qv,qc,qr = ',qx(mgs,lv)*1000.,qx(mgs,lc)*1000.,qx(mgs,lr)*1000.
             write(0,*) ' i, j, k = ',igs(mgs),jy,kgs(mgs)
             write(0,*) ' dtcon1,dtcon2,delta = ',dtcon1,dtcon2,delta
             write(0,*) ' nc,dtp = ',nc,dtp
             write(0,*) ' rwvent,xdia,crw,ccw = ', rwvent(mgs),xdia(mgs,lr,1),cx(mgs,lr),cx(mgs,lc)
             write(0,*) ' fvent,alphar = ',fvent(mgs),alpha(mgs,lr)
             write(0,*) ' xvr,xmasr,xdnr,cwc1 = ',xv(mgs,lr),xmas(mgs,lr),xdn(mgs,lr),cwc1
           endif
            dqvs = dtemp*p380*dtabqvs(ltemq1m)
            qv1m = qv1 + dqv + dqvr
!          qv1mr = qv1r + dqvr

            qvs1m = qvs1 + dqvs
            ss1m = qv1m/qvs1m

    ! check for undersaturation when no ice is present, if so, then reduce time step
          if ( ss1m .lt. 1.  .and. (dqvii + dqvis) .eq. 0.0 ) then
            dtcon = (0.5*dtcon)
            if ( dtcon .ge. dtcon1 ) then
             GOTO 609
            else
             exit
            endif
          endif
! calculate full step:
          dqv  = -(ss1m - 1.)*d1*dtcon
          dqvr = -(ss1m - 1.)*d1r*dtcon


!          write(0,*) 'RK2a dqv1m = ',dqv
          dtemp = -e1*f1*(dqv + dqvr)
          
         ! ltemq1 = ltemq1 + Nint(dtemp*fqsat + 0.5)

         ! 7.6.2016: Test full calc of ltemq
           ltemq1 = (temp1+dtemp-163.15)*fqsati+1.5
           ltemq1 = Min( nqsat, Max(1,ltemq1) )

           if ( ltemq1 .lt. 1 .or. ltemq1 .gt. nqsat ) then
             write(0,*) 'stop in nucond line 1230 '
             write(0,*) ' ltemq1m,icond = ',ltemq1m,icond
             write(0,*) ' dtemp,e1,dqv,dqvr = ', dtemp,e1,dqv,dqvr
           endif
          dqvs = dtemp*p380*dtabqvs(ltemq1)

          qv1 = qv1 + dqv + dqvr

          dqc = dqc - dqv
          dqr = dqr - dqvr

          qvs1 = qvs1 + dqvs
          ss1 = qv1/qvs1
          temp1 = temp1 + dtemp
          if ( temp2 .eq. temp1 .or. ss2 .eq. ss1 .or.  &
     &           ss1 .eq. 1.00 .or.  &
     &      ( n .gt. 10 .and. ss1 .lt. 1.0005 ) ) then
!           write(0,*) 'RK2c break'
           exit
          else
           ss2 = ss1
           temp2 = temp1
           dt1 = dt1 + dtcon
           n = n + 1
          endif
       enddo RK2c


        dcloud = dqc ! qx(mgs,lv) - qv1
        thetap(mgs) = thetap(mgs) + e1*(DCLOUD + dqr)


        if ( eqtset > 2 ) then
           pipert(mgs) = pipert(mgs) + felvpi(mgs)*(DCLOUD + dqr)
        endif
        if ( io_flag .and. nxtra > 1 ) then
           axtra(igs(mgs),jy,kgs(mgs),1) = DCLOUD/dtp
           axtra(igs(mgs),jy,kgs(mgs),2) = axtra(igs(mgs),jy,kgs(mgs),2) + dqr/dtp
        endif
        qwvp(mgs) = qwvp(mgs) - (DCLOUD + dqr)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        qx(mgs,lr) = qx(mgs,lr) + dqr
!        t9(igs(mgs),jy,kgs(mgs)) = t9(igs(mgs),jy,kgs(mgs)) + (DCLOUD + dqr)/dtp*felv(mgs)/(cp*pi0(mgs)) !* &
!!     &                 dx*dy*dz3d(igs(mgs),jy,kgs(mgs))


        if ( lzr > 1 .and. rcond == 2 .and. qx(mgs,lr) .gt. qxmin(lr)   &
     &       .and. cx(mgs,lr) .gt. 1.e-9 ) then
          tmp = qx(mgs,lr)/cx(mgs,lr)
          if ( imurain == 3 ) then
          g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
          else
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
          
          endif
          zx(mgs,lr) = zx(mgs,lr) + g1*(rho0(mgs)/(xdn(mgs,lr)))**2*( 2.*( tmp ) * dqr )
        endif

        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*f1
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
!        es(mgs) = 6.1078e2*tabqvs(ltemq)

!

      endif  ! dcloud .gt. 0.


      else  ! qc .le. qxmin(lc)

!        IF ( ssf(mgs) .gt. 0.0 .and. .not. flag_qndrop ) THEN ! flag_qndrop turns off primary nucleation when using wrf-chem with progn=1
        if ( ssf(mgs) .gt. 0.0 ) then ! .and.  ssmax(mgs) .lt. sscb ) then ! except that wrf-chem does not seem to initialize qc for activated aerosols, so keep this, after all

          if ( iqcinit == 1 ) then

         qvs0   = 380.*exp(17.27*(temg(mgs)-273.)/(temg(mgs)- 36.))/pk(mgs)

         dcloud = Max(0.0, (qx(mgs,lv)-qvs0) / (1.+qvs0*f5/(temg(mgs)-36.)**2) )

          elseif ( iqcinit == 3 ) then
              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felvcp(mgs)/ & 
     &             ((temg(mgs) - cbw)**2))
            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
          
          elseif ( iqcinit == 2 ) then
!              R1=1./(1. + caw*(273.15 - cbw)*qss(mgs)*felv(mgs)/
!     :             (cp*(temg(mgs) - cbw)**2))
!            DCLOUD=R1*(qvap(mgs) - qvs(mgs))  ! KW model adjustment; 
                              ! this will put mass into qc if qv > sqsat exists
         ssmx = ssmxinit

!          IF ( ssf(mgs) > ssmx  .and. ssmax(mgs) < 3.0 ) THEN
!          IF ( ssf(mgs) > ssmx  .and. ccnc(mgs) > 1.0 ) THEN
!          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 5.0 .and. ccnc(mgs) > 0.1*cwnccn(mgs) ) THEN ! this one works
!          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 20.0 ) THEN ! test -- fails
!          IF ( ssf(mgs) > ssmx  .and. ssf(mgs) < 20.0 .and. ccnc(mgs) > 0.1*cwnccn(mgs)) THEN ! test -- is OK
          if ( ssf(mgs) > ssmx  .and. ssf(mgs) < 20.0 .and.  &
             ( ccnc(mgs) > 0.05*cwnccn(mgs) .or. ( ac_opt > 0 .and. ccnc_ac(mgs) - cx(mgs,lc) > 0.0 ) ) ) then ! test
!          IF ( ssf(mgs) > ssmx ) THEN ! original condition
           call QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,dcloud, & 
     &      pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felvcp,ssmx,pk,ngscnt)
          else
            dcloud = 0.0
          endif
         endif
        else
            dcloud = 0.0
        endif
        
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        if ( io_flag .and. nxtra > 1 ) then
           axtra(igs(mgs),jy,kgs(mgs),1) = DCLOUD/dtp
        endif
        theta(mgs) = thetap(mgs) + theta0(mgs)
        temg(mgs) = theta(mgs)*pk(mgs) !( pres(mgs) / poo ) ** cap
!        temg(mgs) = theta2temp( theta(mgs), pres(mgs) )
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )
        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
!        es(mgs) = 6.1078e2*tabqvs(ltemq)

!.... S. TWOMEY (1959)
! Note: get here if there is no previous cloud water and w > 0.
      cn(mgs) = 0.0
      
      if ( ncdebug .ge. 1 ) then
        write(iunit,*) 'at 613: ',qx(mgs,lc),cx(mgs,lc),wvel(mgs),ssmax(mgs),kgs(mgs)
      endif
      
      if (  .not. flag_qndrop ) then ! { do not calculate number of droplets if using wrf-chem

      if ( ac_opt == 0 ) then
        cnuctmp = cnuc(mgs)
      else
        cnuctmp = ccnc_ac(mgs)
      endif
      
!      IF ( ssmax(mgs) .lt. sscb .and. qx(mgs,lc) .gt. qxmin(lc)) THEN
      if ( dcloud .gt. qxmin(lc) .and. wvel(mgs) > 0.0) then
!       CN(mgs) =   CCNE*wvel(mgs)**cnexp ! *Min(1.0,1./dtp) ! 0.3465
       CN(mgs) =   CCNE0*cnuctmp**(2./(2.+cck))*wvel(mgs)**cnexp ! *Min(1.0,1./dtp) ! 0.3465
        if ( ny .le. 2 .and. cn(mgs) .gt. 0.0    &
     &                    .and. ncdebug .ge. 1 ) then 
          write(iunit,*) 'CN: ',cn(mgs)*1.e-6, cx(mgs,lc)*1.e-6, qx(mgs,lc)*1.e3,   &
     &       wvel(mgs), dcloud*1.e3
          if ( cn(mgs) .gt. 1.0 ) write(iunit,*) 'cwrad = ',   &
     &       1.e6*(rho0(mgs)*qx(mgs,lc)/cn(mgs)*cwc1)**c1f3,   &
     &   igs(mgs),kgs(mgs),temcg(mgs),    &
     &   1.e3*an(igs(mgs),jgs,kgs(mgs)-1,lc)
        endif
        if ( iccwflg .eq. 1 ) then
          cn(mgs) = Min(cwccn*rho0(mgs)/rho00, Max(cn(mgs),   &
     &       rho0(mgs)*qx(mgs,lc)/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3)))
        endif
      else
       cn(mgs) = 0.0
       dcloud = 0.0
!          cn(mgs) = Min(cwccn,    &
!     &       rho0(mgs)*dcloud/(xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3) )
      endif

      if ( cn(mgs) .gt. 0.0 ) then
       if ( ac_opt == 0 ) then
         if ( cn(mgs) .gt. ccnc(mgs) ) then
           cn(mgs) = ccnc(mgs)
!          ccnc(mgs) = 0.0
         endif
       else 
         cn(mgs) = Min( cn(mgs), ccnc_ac(mgs) )
       endif
!      cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
      if ( irenuc <= 2 .and. lccna < 1  ) ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
      ccna(mgs) = ccna(mgs) + cn(mgs)
      endif

!       write(91,*) 'nuc1: cn, ix, kz = ',cn(mgs),igs(mgs),kgs(mgs),wvel(mgs),cnexp,ccnc(mgs)

      if( CN(mgs) .GT. cx(mgs,lc) ) cx(mgs,lc) = CN(mgs)
      if( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .le. qxmin(lc) ) then
        cx(mgs,lc) = 0.
      else
        cx(mgs,lc) = Min(cx(mgs,lc),rho0(mgs)*Max(0.0,qx(mgs,lc))/cwmasn)
      endif
      
      endif ! }.not. flag_qndrop

        GOTO 613
        
        end if ! qc .gt. 0.

!        ES=EES(PIB(K)*PT)
!        SQSAT=EPSI*ES/(PB(K)*1000.-ES)

!.... CLOUD NUCLEATION
!      T=PIB(K)*PT
!      ES=1.E3*PB(K)*QV/EPSI

      if ( wvel(mgs) .le. 0. ) go to 616
      if ( cx(mgs,lc) .le. 0. )  go to 613                             !TWOMEY (1959) Nucleation
      if ( kzbeg-1+kgs(mgs) .GT. 1 .and. qx(mgs,lc) .le. qxmin(lc)) go to 613  !TWOMEY (1959) Nucleation
      if ( kzbeg-1+kgs(mgs) .eq. 1 .and. wvel(mgs) .gt. 0. ) go to 613         !TWOMEY (1959) Nucleation
!.... ATTEMPT ZIEGLER CLOUD NUCLEATION IN CLOUD INTERIOR UNLESS...
  616 if ( ssf(mgs) .LE. SUPCB .AND. wvel(mgs) .GT. 0. ) go to 631 !... weakly saturated updraft
      if ( kzbeg-1+kgs(mgs) .GT. 1 .AND. kzbeg-1+kgs(mgs) .LT. nzend-1 .AND.  &
     &    (ssfkp1(mgs) .GE. SUPMX .OR. &
     &     ssf(mgs)    .GE. SUPMX .OR. &
     &     ssfkm1(mgs) .GE. SUPMX)) go to 631                      !... too much vapour
      if (ssf(mgs) .LT. 1.E-10 .OR. ssf(mgs) .GE. SUPMX) go to 631 !... at the extremes for ss

!
! get here if ( qc > 0 and ss > supcb) or (w < 0)
!

      if (ndebug .gt. 0) write(0,*) "ICEZVD_DR: Entered Ziegler Cloud Nucleation" !mpidebug

      DSSDZ=0.
      r2dzm=0.50/dz3d(igs(mgs),jy,kgs(mgs))

      if ( irenuc >= 0 .and. ac_opt == 0 .and. .not. flag_qndrop ) then ! turn off nucleation when flag_qndrop (using WRF-CHEM for activation)

      if ( irenuc < 2 ) then !{

        if ( kzend == nzend ) then
          t0p3 = t0(igs(mgs),jgs,Min(kze,kgs(mgs)+3))
          t0p1 = t0(igs(mgs),jgs,Min(kze,kgs(mgs)+1))
        else
          t0p3 = t0(igs(mgs),jgs,kgs(mgs)+3)
          t0p1 = t0(igs(mgs),jgs,kgs(mgs)+1)
        endif

      if ( ( ssf(mgs) .gt. ssmax(mgs) .or.  irenuc .eq. 1 ) &
     &   .and.  ( ( lccn .lt. 1 .and.  &
     &            cx(mgs,lc) .lt. cwccn*(Min(1.0,rho0(mgs)))) .or. &
     &    ( lccn .gt. 1 .and. ccnc(mgs) .gt. 0. )   ) &
     &    ) then
      if( kzbeg-1+kgs(mgs) .GT. 1 .AND. kzbeg-1+kgs(mgs) .LT. nzend-1 &
     &  .and. ssf(mgs) .gt. 0.0 &
     &  .and. ssfkp1(mgs) .LT. SUPMX .and. ssfkp1(mgs) .ge. 0.0  &
     &  .AND. ssfkm1(mgs) .LT. SUPMX .AND. ssfkm1(mgs) .ge. 0.0  &
     &  .AND. ssfkp1(mgs) .gt. ssfkm1(mgs)  &
     &  .and. t0p3 .gt. 233.2) then
          DSSDZ = (ssfkp1(mgs) - ssfkm1(mgs))*R2DZM
!
! otherwise check for cloud base condition with updraft:
!
        elseif( kzbeg-1+kgs(mgs) .GT. 1 .AND. kzbeg-1+kgs(mgs) .LT. nzend-1 &
!        IF( kgs(mgs) .GT. 1 .AND. kgs(mgs) .LT. NZ-1 & !)
     &  .and. ssf(mgs) .gt. 0.0  .and. wvel(mgs) .gt. 0.0 &
     &  .and. ssfkp1(mgs) .gt. 0.0   &
     &  .AND. ssfkm1(mgs) .le. 0.0 .and. wvelkm1(mgs) .gt. 0.0 &
     &  .AND. ssf(mgs) .gt. ssfkm1(mgs)  &
     &  .and. t0p1 .gt. 233.2) then
         DSSDZ = 2.*(ssf(mgs) - ssfkm1(mgs))*R2DZM  ! 1-sided difference
        endif

       endif
!
!CLZ  IF(wijk.LE.0.) CN=CCN*ssfilt(ix,jy,kz)**CCK
! note: CCN -> cwccn, DELT -> dtp
      c1 = Max(0.0, rho0(mgs)*(qx(mgs,lv) - qss(mgs))/ &
     &        (xdn(mgs,lc)*(4.*pi/3.)*(4.e-6)**3))
      if ( lccn .lt. 1 ) then
       CN(mgs) = cwccn*rho0(mgs)/rho00*CCK*ssf(mgs)**CCKM*dtp*   &
     & Max(0.0,    &
     &         (wvel(mgs)*DSSDZ) )      ! probably the vertical gradient dominates
      else
      CN(mgs) =  &
     &    Min(ccnc(mgs), cnuc(mgs)*CCK*ssf(mgs)**CCKM*dtp*   &
     & Max(0.0,    &
     &         ( wvel(mgs)*DSSDZ) )  )
!      IF ( cn(mgs) .gt. 0 ) ccnc(mgs) = ccnc(mgs) - cn(mgs)
      endif

      if ( cn(mgs) .gt. 0.0 ) then
       if ( ccnc(mgs) .lt. 5.e7 .and. cn(mgs) .ge. 5.e7 ) then
          cn(mgs) = 5.e7
          ccnc(mgs) = 0.0
       elseif ( cn(mgs) .gt. ccnc(mgs) ) then
         cn(mgs) = ccnc(mgs)
         ccnc(mgs) = 0.0
       endif
      cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
      endif

      elseif ( irenuc == 2 ) then !} { 
      ! simple Twomey scheme
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       CN(mgs) =   CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp ! *Min(1.0,1./dtp) ! 0.3465
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       CN(mgs) = Min( CN(mgs), Max(0.0, (cnuc(mgs) - ccna(mgs) )) )
       
        if ( .false. .and. ny <= 2 ) then
          write(0,*) 'i,k, cwmasn = ',igs(mgs),kgs(mgs),cwmasn
          write(0,*) 'wvel, cnuc, cn = ',wvel(mgs),cnuc(mgs),cn(mgs)
          write(0,*) 'ccne0,cnexp,cck = ',ccne0,cnexp,cck
          write(0,*) 'part1, part2 = ',CCNE0*cnuc(mgs)**(2./(2.+cck)), Max(0.0,wvel(mgs))**cnexp
          write(0,*) 'ccnc, dqc, dqc/cwmasn = ',ccnc(mgs), dqc, 0.5*dqc/cwmasn
        endif
       
       if ( icnuclimit > 0 ) then 
         tmp = ccnc(mgs) + cx(mgs,lc)
         if ( tmp < 330.34e6 ) then
           ccwmax = 1.1173e6 * (1.e-6*tmp)**0.9504
         else
           ccwmax = 21.57e6 * (1.e-6*tmp)**0.44
         endif
         
!         IF ( cn(mgs) > 0. ) THEN
!          write(0,*) 'cn,tmp,ccwmax,cx,c-cx = ',cn(mgs),tmp,ccwmax,cx(mgs,lc),ccwmax - cx(mgs,lc) 
!         ENDIF
         
        cn(mgs) = Max( 0.0, Min( cn(mgs), ccwmax - cx(mgs,lc) ) )
       
       endif
       
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       if ( lccna < 1 ) ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))

      elseif ( irenuc == 3 ) then !} { 
      ! Phillips Donner Garner 2007
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
!       CN(mgs) =   cwccn*Min(ssf(mgs),ssfcut)**cck 

! Need to calculate new ssf since condensation has happened:
         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          c1= pqs(mgs)*tabqvs(ltemq)

          ssf(mgs) = 0.0
          if ( c1 > 0. ) then
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          endif
       CN(mgs) =   cnuc(mgs)*Min(1.0, (ssf(mgs))**cck ) ! 

       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       ! 6/13/2016: Phillips et al. appears not to decrement CCN, but only increments CCNa.
       ! This would allow an initially non-homogeneous (vertically, e.g.) initial value of CCN/rho_air
        ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs)) 
       
      elseif ( irenuc == 4 ) then !} { 
      ! modification of Phillips Donner Garner 2007
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
!       CN(mgs) =   CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp
!       cn(mgs) = Min( cn(mgs), cnuc(mgs) )
! Need to calculate new ssf since condensation has happened:
         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          c1= pqs(mgs)*tabqvs(ltemq)
          if ( c1 > 0. ) then
            ssf(mgs) = Max(0.0, 100.*((qv0(mgs) + qwvp(mgs))/c1 - 1.0) )  ! from "new" values
          else
            ssf(mgs) = 0.0
          endif
       CN(mgs) =   cnuc(mgs)*Min(ssf2kmax, ssf(mgs)**cck) ! this allows cn(mgs) > cnuc(mgs)

       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       
       if ( cn(mgs) > 0.0 ) then
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       ! ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs)) 
       
       dcrit = 2.0*2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
        endif
       ! 6/13/2016: Phillips et al. appears not to decrement CCN, but only increments CCNa.
       ! This would allow an initially non-homogeneous (vertically, e.g.) initial value of CCN/rho_air
!        ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
       


      elseif ( irenuc == 6 ) then !} { 

      ! simple Twomey scheme but limit activation to try to do most activation near cloud base, but keep some CCN available for renuclation
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       cn(mgs) = 0.0
!       IF ( ccna(mgs) < 0.7*cnuc(mgs) .and. ccnc(mgs) > 0.69*cnuc(mgs) - ccna(mgs)) THEN ! here, assume we are near cloud base and use Twomey formulation
       if ( ccna(mgs) < 0.7*cnuc(mgs) ) then ! here, assume we are near cloud base and use Twomey formulation
         CN(mgs) =  Min( 0.9*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
!         IF ( cn(mgs) + ccna(mgs) > 0.71*cnuc ) THEN
         ! prevent this branch from activating more than 70% of CCN
           CN(mgs) = Min( CN(mgs), Max(0.0, (0.7*cnuc(mgs) - ccna(mgs) )) )
!           CN(mgs) = Min( CN(mgs), Max(0.0, 0.71*ccnc(mgs) - ccna(mgs) ) )
           
       else
        ! if a large fraction of CCN have been activated, then assume we are in the cloud interior and use local SSw as in Phillips et al. 2007.

         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
!          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

!          c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
          c1= pqs(mgs)*tabqvs(ltemq)
          if ( c1 > 0. ) then
            ssf(mgs) = Max(0.0, 100.*((qv0(mgs) + qwvp(mgs))/c1 - 1.0) )  ! from "new" values
          else
            ssf(mgs) = 0.0
          endif

!        CN(mgs) = cnuc(mgs)*Min(0.99, Min(ssf(mgs),ssfcut)**cck ) ! 
         CN(mgs) =   cnuc(mgs)*Min(2.0, Max(0.0,ssf(mgs))**cck ) ! 
!         CN(mgs) =   cnuc(mgs)*Min(ssf(mgs),ssfcut)**cck ! 


        CN(mgs) = Min(0.01*cnuc(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from
!        cn(mgs) = 0.0
       endif
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
!       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       
       if ( cn(mgs) > 0.0 ) then
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       dcrit = 2.0*2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
        endif
      elseif ( irenuc == 5 ) then !} { 

      ! modification of Phillips Donner Garner 2007
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
!      CN(mgs) =  Min( 0.91*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
       CN(mgs) =  Min( cnuc(mgs),  CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )

         
        if ( ccna(mgs) >= cnuc(mgs) ) then ! apply limit after all "base" CCN have been depleted
        temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

          c1= pqs(mgs)*tabqvs(ltemq)
          if ( c1 > 0. ) then
            ssf(mgs) = Max(0.0, 100.*((qv0(mgs) + qwvp(mgs))/c1 - 1.0) )  ! from "new" values
          else
            ssf(mgs) = 0.0
          endif
          

       CN(mgs) =   Max( cn(mgs), cnuc(mgs)*Min(ssf2kmax, ssf(mgs)**cck) ) ! this allows cn(mgs) > cnuc(mgs)

   !    cn(mgs) = Min( cn(mgs), cnuc(mgs) )

!       IF ( ccna(mgs) >= cnuc(mgs) ) THEN ! apply limit after all "base" CCN have been depleted
       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
       
       else
         CN(mgs) =  Min( cn(mgs), cnuc(mgs) - ccna(mgs) ) ! no more than remaining "base" CCN
       endif
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
!       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       dcrit = 2.0*2.0e-6
       dcloud = 1000.*dcrit**3*Pi/6.
 !      cn(mgs) = Min(cn(mgs), 0.5*dqc/dcloud) ! limit the nucleation mass to half of the condensation mass
       ! check new droplet size:
         ! tmp is number of droplets at diameter dcrit
         tmp = Max(0.0,  rho0(mgs)*qx(mgs,lc)/dcloud - cx(mgs,lc)) ! (cx(mgs,lc) + cn(mgs))
         cn(mgs) = Min(tmp, cn(mgs) )

      
       if ( cn(mgs) > 0.0 ) then
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       dcrit = 2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
        endif
       ! 6/13/2016: Phillips et al. appears not to decrement CCN, but only increments CCNa.
       ! This would allow an initially non-homogeneous (vertically, e.g.) initial value of CCN/rho_air
       ! ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
      elseif ( irenuc == 7 .or. irenuc == 17 ) then !} { 

      ! simple Twomey scheme but limit activation to try to do most activation near cloud base, but keep some CCN available for renuclation
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       cn(mgs) = 0.0
       if ( irenuc == 7 ) then
         frac = 0.9
       else
         frac = 0.98
       endif
!       IF ( ccna(mgs) < 0.7*cnuc(mgs) .and. ccnc(mgs) > 0.69*cnuc(mgs) - ccna(mgs)) THEN ! here, assume we are near cloud base and use Twomey formulation
       if ( ccna(mgs) < frac*cnuc(mgs) ) then ! { here, assume we are near cloud base and use Twomey formulation
         CN(mgs) =  Min( (frac+0.01)*cnuc(mgs), CCNE0*cnuc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
!         IF ( cn(mgs) + ccna(mgs) > 0.71*cnuc ) THEN
         ! prevent this branch from activating more than 70% of CCN
           CN(mgs) = Min( CN(mgs), Max(0.0, (frac*cnuc(mgs) - ccna(mgs) )) )
!           CN(mgs) = Min( CN(mgs), Max(0.0, 0.71*ccnc(mgs) - ccna(mgs) ) )
         !  write(0,*) '1: k,cn = ',kgs(mgs),cn(mgs),ssf(mgs)
!!           IF ( ccncuf(mgs) > 0.0 .and. cn(mgs) < 1.e-3 .and. ssmax(mgs) > 1.0 ) THEN
!           IF ( ccncuf(mgs) > 0.0 .and. ssf(mgs) > ssmxuf .and. ssmax(mgs) > ssmxuf ) THEN
!            CNuf(mgs) =  Min( ccncuf(mgs), CCNE0*ccncuf(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
          !  IF ( cnuf(mgs) >= 0.0 ) write(0,*) '1: cnuf, k = ',cnuf(mgs),ccncuf(mgs),kgs(mgs)
!           ENDIF

           
       else ! }{
        ! if a large fraction of CCN have been activated, then assume we are in the cloud interior and use local SSw as in Phillips et al. 2007.

         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
!          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

        !  c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
          c1= pqs(mgs)*tabqvs(ltemq)

          ssf(mgs) = 0.0
          if ( c1 > 0. ) then
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          endif

!          IF ( ssf(mgs) <= 1.0 .or. cnuc(mgs) > ccna(mgs) ) THEN
          if ( ssf(mgs) <= 1.0 ) then
          CN(mgs) =   cnuc(mgs)*Min(1.0, Max(0.0,ssf(mgs))**cck ) ! 
          else
          CN(mgs) =   cnuc(mgs)*Min(2.0, Max(0.0,0.03*(ssf(mgs)-1.0)+1.)**cck ) !           
!          write(0,*) 'iren7: ssf,ssmx = ',ssf(mgs),ssmax(mgs),cn(mgs),ccna(mgs),cnuc(mgs)
!          write(0,*) 'c1,qv = ',c1,qx(mgs,lv),temp1,ltemq
          endif

         !  write(0,*) 'k,cn = ',kgs(mgs),cn(mgs),ssf(mgs)
         !  write(0,*) 'ccn-ccna = ',cnuc(mgs) - ccna(mgs),ccnc(mgs) - ccna(mgs)
!           IF ( ccncuf(mgs) > 0.0 .and. cn(mgs) < 1.e-3 .and. ssmax(mgs) > 1.0 ) THEN
           if ( ccncuf(mgs) > 0.0 .and. ssf(mgs) > ssmxuf .and. ( ssmax(mgs) > ssmxuf .or. lss < 1 ) ) then
            CNuf(mgs) =  Min( ccncuf(mgs), CCNE0*ccncuf(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp )! *Min(1.0,1./dtp) ! 0.3465
          !  IF ( cnuf(mgs) >= 0.0 ) write(0,*) 'cnuf, k = ',cnuf(mgs),ccncuf(mgs),kgs(mgs)
           endif
          

!        CN(mgs) = Min( Min(0.1,ssf(mgs)-1.)*cnuc(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from
!        CN(mgs) = Min( Min(0.5*cx(mgs,lc), Min(0.1,ssf(mgs)/100.)*cnuc(mgs)), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from
        
        CN(mgs) = Min(0.01*cnuc(mgs), Max( 0.0, CN(mgs) - ccna(mgs) ) ) ! this was from

       endif ! }
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
!       CN(mgs) = Min(cn(mgs), ccnc(mgs))
!       cn(mgs) = Min(cn(mgs), 0.5*dqc/cwmasn) ! limit the nucleation mass to half of the condensation mass
       

        if ( icnuclimit > 0 ) then
! max droplet conc. based on Chandrakar et al. (2016) and Konwar et al. (2012)
           tmp = ccnc(mgs) - ccna(mgs) + cx(mgs,lc)
           if ( tmp < 330.34e6 ) then
             ccwmax = 1.1173e6 * (1.e-6*tmp)**0.9504
           else
             ccwmax = 21.57e6 * (1.e-6*tmp)**0.44
           endif
          
           cn(mgs) = Max( 0.0, Min( cn(mgs), ccwmax - cx(mgs,lc) ) )
           
        endif

       if ( cn(mgs) + cnuf(mgs) > 0.0 ) then

       dcrit = 2.0*2.0e-6
       dcloud = 1000.*dcrit**3*Pi/6.
 !      cn(mgs) = Min(cn(mgs), 0.5*dqc/dcloud) ! limit the nucleation mass to half of the condensation mass
       ! check new droplet size:
         ! tmp is number of droplets at diameter dcrit
         tmp = Max(0.0,  rho0(mgs)*qx(mgs,lc)/dcloud - cx(mgs,lc)) ! (cx(mgs,lc) + cn(mgs))
         cn(mgs) = Min(tmp, cn(mgs) )

        cx(mgs,lc) = cx(mgs,lc) + cn(mgs) + cnuf(mgs)
 
 
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       
       dcrit = 2.0*2.5e-7
       dcloud = 1000.*dcrit**3*Pi/6.*(cn(mgs) + cnuf(mgs) )
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
         ccncuf(mgs) = Max(0.0, ccncuf(mgs) - cnuf(mgs))
        endif

      elseif ( irenuc == 8 ) then !} { 
      ! simple Twomey scheme
!      if (ndebug .gt. 0) write(0,*) 'ICEZVD_DR:  Cloud reNucleation, wvel = ',wvel(mgs)
       
       cn(mgs) = 0.0
       
       if ( ccnc(mgs) > 0. ) then
       CN(mgs) =   CCNE0*ccnc(mgs)**(2./(2.+cck))*Max(0.0,wvel(mgs))**cnexp ! *Min(1.0,1./dtp) ! 0.3465
!      ccne = ccnefac*1.e6*(1.e-6*Abs(cwccn))**(2./(2.+cck))
!!!       CN(mgs) = Max( 0.0, CN(mgs) - ccna(mgs) ) ! this was from
               ! Philips, Donner et al. 2007, but results in too much limitation of
               ! nucleation
       CN(mgs) = Min(cn(mgs), ccnc(mgs))
       
       elseif ( cx(mgs,lc) < 0.01e9 ) then

        ! if a large fraction of CCN have been activated, then assume we are in the cloud interior and use local SSw as in Phillips et al. 2007.

         temp1 = (theta0(mgs)+thetap(mgs))*pk(mgs) ! t77(ix,jy,kz)
!          t0(ix,jy,kz) = temp1
          ltemq = Int( (temp1-163.15)/fqsat+1.5 )
         ltemq = Min( nqsat, Max(1,ltemq) )

        !  c1 = t00(igs(mgs),jy,kgs(mgs))*tabqvs(ltemq)
          c1= pqs(mgs)*tabqvs(ltemq)

          ssf(mgs) = 0.0
          if ( c1 > 0. ) then
            ssf(mgs) = 100.*(qx(mgs,lv)/c1 - 1.0)  ! from "new" values
          endif

!          IF ( ssf(mgs) <= 1.0 .or. cnuc(mgs) > ccna(mgs) ) THEN
          if ( ssf(mgs) <= 1.0 ) then
          CN(mgs) = 0.0
          else
!           CN(mgs) = 0.01e9*rho0(mgs)/rho00*Min(2.0, Max(0.0,0.03*(ssf(mgs)-1.0)+1.)**cck ) - cx(mgs,lc) !           
           CN(mgs) = 0.01e9*Min(2.0, Max(0.0,0.03*(ssf(mgs)-1.0)+1.)**cck ) - cx(mgs,lc) !           
          endif
       
       endif

       if ( cn(mgs) > 0.0 ) then
       cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
       
       ! ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
       
       ! create some small droplets at minimum size (CP 2000), although it adds very little liquid
       
       dcrit = 2.0*2.5e-7
       
       dcloud = 1000.*dcrit**3*Pi/6.*cn(mgs)
        qx(mgs,lc) = qx(mgs,lc) + DCLOUD
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*DCLOUD/(pi0(mgs))
        qwvp(mgs) = qwvp(mgs) - DCLOUD
  !      ccnc(mgs) = Max(0.0, ccnc(mgs) - cn(mgs))
        endif
       


      endif ! }

      ccna(mgs) = ccna(mgs) + cn(mgs)

      endif ! irenuc >= 0 .and. .not. flag_qndrop

      if( cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .LE. qxmin(lc)) cx(mgs,lc)=0.
      go to 631
!.... NUCLEATION ON CLOUD INFLOW BOUNDARY POINT

  613 continue

  631  continue

!
! Check for supersaturation greater than ssmx and adjust down
!
       ssmx = maxsupersat
       qv1 = qv0(mgs) + qwvp(mgs)
       qvs1 = qvs(mgs)
       
!       IF ( flag_qndrop .and. do_satadj_for_wrfchem ) ssmx = 1.04 ! set lower threshold for progn=1 when using WRF-CHEM

       if ( qv1 .gt. (ssmx*qvs1) ) then
! use line below to disable saturation adjustment when flag_qndrop is true
!       IF ( qv1 .gt. (ssmx*qvs1) .and. .not. flag_qndrop ) THEN
        
         ss1 = qv1/qvs1

        ssmx = 100.*(ssmx - 1.0)
        
        qvex = 0.0

        call QVEXCESS(ngs,mgs,qwvp,qv0,qx(1,lc),pres,thetap,theta0,qvex,   &
     &    pi0,tabqvs,nqsat,fqsat,cbw,fcqv1,felvcp,ssmx,pk,ngscnt)



        if ( qvex .gt. 0.0 ) then
        thetap(mgs) = thetap(mgs) + felvcp(mgs)*qvex/(pi0(mgs))
        if ( io_flag .and. nxtra > 1 ) then
           axtra(igs(mgs),jy,kgs(mgs),1) = axtra(igs(mgs),jy,kgs(mgs),1) + qvex/dtp
        endif
        qwvp(mgs) = qwvp(mgs) - qvex
        qx(mgs,lc) = qx(mgs,lc) + qvex
        if ( .not. flag_qndrop) then
          if ( imaxsupopt == 1 ) then
            cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, xmas(mgs,lc) )  )
          elseif ( imaxsupopt == 2 ) then
            cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, Max(cwmas30,xmas(mgs,lc)) )  )
          elseif ( imaxsupopt == 3 ) then
            cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, Max(cwmasx,xmas(mgs,lc)) )  )
!            cn(mgs) = 1.5*cxmin
          elseif ( imaxsupopt == 4 ) then
            cn(mgs) = Min( Max(ccnc(mgs),cwnccn(mgs)), rho0(mgs)*qvex/Max( cwmasn5, Max(cwmas20,xmas(mgs,lc)) )  )
          endif
        if ( lccna > 1 ) then
          ccna(mgs) = ccna(mgs) + cn(mgs)
        else
          ccnc(mgs) = Max( 0.0, ccnc(mgs) - cn(mgs) )
        endif
        cx(mgs,lc) = cx(mgs,lc) + cn(mgs)
        endif
        
!        write(iunit,*) 'theta = ',theta0(mgs) + thetap(mgs)

!        temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** cap

        endif

       
       endif

!
! Calculate droplet volume and check if it is within bounds.
!  Adjust if necessary
!  
!      if (ndebug .gt. 0) write(0,*) "ICEZVD_DR: check droplet volume" 


!      cx(mgs,lc) = Min( cwnccn(mgs), cx(mgs,lc) )
      if( cx(mgs,lc) > cxmin .AND. qx(mgs,lc) .GT. qxmin(lc)) then
!        SVC(mgs) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc)*xdn(mgs,lc))
        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/(cx(mgs,lc))
        
       if (  xmas(mgs,lc) < cwmasn .or.  xmas(mgs,lc) > cwmasx ) then
        tmp = cx(mgs,lc)
        xmas(mgs,lc) = Min( xmas(mgs,lc), cwmasx )
        xmas(mgs,lc) = Max( xmas(mgs,lc), cwmasn )
        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc)
!        IF ( cx(mgs,lc) > tmp*1.1 ) THEN
!          write(0,*) 'nucond: kgs, ccw1,2 = ',kgs(mgs),tmp,cx(mgs,lc)
!        ENDIF
       endif
      endif


!      IF( cx(mgs,lc) .GT. 10.e6 .AND. qx(mgs,lc) .GT. qxmin(lc) ) GO TO 681
!        ccwtmp = cx(mgs,lc)
!        cwmastmp = xmas(mgs,lc)
!       xmas(mgs,lc) = Max(xmas(mgs,lc), cwmasn)
!       IF (qx(mgs,lc) .GT. qxmin(lc) .AND. cx(mgs,lc) .le. 0.) THEN
!          cx(mgs,lc) = Min(0.5*cwccn,rho0(mgs)*qx(mgs,lc)/xmas(mgs,lc))
!          xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
!       ENDIF
!      IF (cx(mgs,lc) .GT. 0. .AND. qx(mgs,lc) .GT. qxmin(lc))    &
!     &        xmas(mgs,lc) = rho0(mgs)*qx(mgs,lc)/cx(mgs,lc)
!      IF (qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .LT. cwmasn)    &
!     &          xmas(mgs,lc) = cwmasn
!      IF (qx(mgs,lc) .GT. qxmin(lc) .AND. xmas(mgs,lc) .GT. cwmasx)    &
!     &    xmas(mgs,lc) = cwmasx
!      IF ( qx(mgs,lc) .gt. qxmin(lc) ) THEN
!        cx(mgs,lc) = rho0(mgs)*qx(mgs,lc)/Max(cwmasn,xmas(mgs,lc))
!      ENDIF
!        
!
! 681  CONTINUE

        
      if ( ipconc .ge. 3 .and. rcond == 2 ) then

        
        if (cx(mgs,lr) .GT. 0. .AND. qx(mgs,lr) .GT. qxmin(lr))    &
     &       xv(mgs,lr)=rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
        if (xv(mgs,lr) .GT. xvmx(lr)) xv(mgs,lr) = xvmx(lr)
        if (xv(mgs,lr) .LT. xvmn(lr)) xv(mgs,lr) = xvmn(lr)

      endif



      enddo ! mgs


! ################################################################
      do mgs=1,ngscnt
      if ( lss > 1 .and. ssf(mgs) .gt. ssmax(mgs)    &
     &  .and. ( idecss .eq. 0 .or. qx(mgs,lc) .gt. qxmin(lc)) ) then
        ssmax(mgs) = ssf(mgs)
      endif
      enddo
!

      do mgs = 1,ngscnt
      an(igs(mgs),jy,kgs(mgs),lt) = theta0(mgs) + thetap(mgs)
      an(igs(mgs),jy,kgs(mgs),lv) =  qv0(mgs) + qwvp(mgs)
!      tmp3d(igs(mgs),jy,kgs(mgs)) = tmp3d(igs(mgs),jy,kgs(mgs)) + t9(igs(mgs),jy,kgs(mgs)) !  pi0(mgs) ! wvdf(mgs) ! ssf(mgs) ! cn(mgs)
!
      if ( eqtset > 2 ) then
        p2(igs(mgs),jy,kgs(mgs)) = pipert(mgs)
      endif

       if ( ido(lc) .eq. 1 )  then
        an(igs(mgs),jy,kgs(mgs),lc) = qx(mgs,lc) +    &
     &    min( an(igs(mgs),jy,kgs(mgs),lc), 0.0 )
!        qx(mgs,lc) = an(igs(mgs),jy,kgs(mgs),lc)
       end if
!

       if ( ido(lr) .eq. 1 .and. rcond == 2 )  then
        an(igs(mgs),jy,kgs(mgs),lr) = qx(mgs,lr) +    &
     &    min( an(igs(mgs),jy,kgs(mgs),lr), 0.0 )
!        qx(mgs,lr) = an(igs(mgs),jy,kgs(mgs),lr)
       end if

        if ( lzr > 1 .and. rcond == 2 ) then
        an(igs(mgs),jy,kgs(mgs),lzr) = zx(mgs,lr) +  &
     &    min( an(igs(mgs),jy,kgs(mgs),lzr), 0.0 )
        endif


       if (  ipconc .ge. 2 ) then
        an(igs(mgs),jy,kgs(mgs),lnc) = Max(cx(mgs,lc) , 0.0)
        if ( lss > 1 ) an(igs(mgs),jy,kgs(mgs),lss) = Max( 0.0, ssmax(mgs) )
        if ( ac_opt == 0 ) then
          if ( lccn .gt. 1 .and. lccna .lt. 1  ) then
            an(igs(mgs),jy,kgs(mgs),lccn) = Max(0.0,  ccnc(mgs) )
          endif
        endif
        if ( lccnuf .gt. 1 .and. .not. ( lccna .gt. 1 .and. i_uf_or_ccn > 0 ) ) then
          an(igs(mgs),jy,kgs(mgs),lccnuf) = Max(0.0,  ccncuf(mgs) )
        endif
        if ( lccna .gt. 1 ) then
          an(igs(mgs),jy,kgs(mgs),lccna) = Max(0.0, ccna(mgs) )
        endif
       endif
       if (  ipconc .ge. 3 .and. rcond == 2 ) then
        an(igs(mgs),jy,kgs(mgs),lnr) = Max(cx(mgs,lr) , 0.0)
       endif
      end do


29998 continue


      if ( kz .gt. nz-1 .and. ix .ge. nxi) then
        if ( ix .ge. nxi ) then
         go to 2200 ! exit gather scatter
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz
      end if

      if ( ix .ge. nxi ) then
        nxmpb = 1
        nzmpb = kz+1
      else
       nxmpb = ix+1
      end if

 2000 continue ! inumgs
 2200 continue
!
!  end of gather scatter (for this jy slice)


!#ifdef COMMAS
!    GOTO 9999
!#endif

! Redistribute inappreciable cloud particles and charge
!
! Redistribution everywhere in the domain...
!
    if ( .true. ) then
      
      frac = 1.0 ! 0.25 ! 1.0 ! 0.2
!
!  alternate test version for ipconc .ge. 3
!  just vaporize stuff to prevent noise in the number concentrations


      do kz = 1,nz
!      do jy = 1,1
      do ix = 1,nxi
      
      t0(ix,jy,kz) = an(ix,jy,kz,lt)*t77(ix,jy,kz)
      
      zerocx(:) = .false.
      do il = lc,lhab
       if ( iresetmoments == 1 .or. iresetmoments == il ) then
        if ( ln(il) > 1 ) zerocx(il) = ( an(ix,jy,kz,ln(il)) < cxmin )
        if ( lz(il) > 1 ) zerocx(il) = ( zerocx(il) .or. an(ix,jy,kz,lz(il)) < zxmin )
       else
        if ( il == lc ) then
          if ( ln(il) > 1 ) zerocx(il) = ( an(ix,jy,kz,ln(il)) <= 0 ) .and. .not. flag_qndrop ! do not reset if progn=1 (WRF-CHEM)
        else
         if ( ln(il) > 1 ) zerocx(il) = ( an(ix,jy,kz,ln(il)) <= 0 )
        endif
       endif
      enddo

      if ( lhl .gt. 1 ) then
      
      if ( lzhl .gt. 1 ) then

        an(ix,jy,kz,lzhl) = Max(0.0, an(ix,jy,kz,lzhl) )
        
        if ( an(ix,jy,kz,lhl) .ge. frac*qxmin(lhl) .and. rescale_low_alpha ) then ! check 6th moment
          
          if ( an(ix,jy,kz,lnhl) .gt. 0.0 ) then

           if ( lvhl .gt. 1 ) then
             if ( an(ix,jy,kz,lvhl) .gt. 0.0 ) then
               hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
             else
               hwdn = xdn0(lhl)
             endif
             hwdn = Max( xdnmn(lhl), hwdn )
           else
             hwdn = xdn0(lhl)
           endif

             chw = an(ix,jy,kz,lnhl)
             g1 = (6.0+alphamin)*(5.0+alphamin)*(4.0+alphamin)/   &
     &            ((3.0+alphamin)*(2.0+alphamin)*(1.0+alphamin))
             z1 = g1*dn(ix,jy,kz)**2*( an(ix,jy,kz,lhl) )*an(ix,jy,kz,lhl)/chw
             z1 = z1*(6./(pi*hwdn))**2
          else
             z1 = 0.0
          endif
          
          an(ix,jy,kz,lzhl) = Min( z1, an(ix,jy,kz,lzhl) )
          
          if (  an(ix,jy,kz,lnhl) .lt. 1.e-5 ) then
!            an(ix,jy,kz,lzhl) = 0.9*an(ix,jy,kz,lzhl)
          endif
        endif
        
      endif !lzhl
      
      if ( an(ix,jy,kz,lhl) .lt. frac*qxmin(lhl) .or. zerocx(lhl) ) then

!        IF ( an(ix,jy,kz,lhl) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lhl)
          an(ix,jy,kz,lhl) = 0.0
!        ENDIF

        if ( ipconc .ge. 5 ) then ! .and. an(ix,jy,kz,lnh) .gt. 0.0 ) then
          an(ix,jy,kz,lnhl) = 0.0
        endif

        if ( lvhl .gt. 1 ) then
           an(ix,jy,kz,lvhl) = 0.0
        endif

        if ( lhlw .gt. 1 ) then
           an(ix,jy,kz,lhlw) = 0.0
        endif

        if ( lnhlf .gt. 1 ) then
           an(ix,jy,kz,lnhlf) = 0.0
        endif
      
        if ( lzhl .gt. 1 ) then
           an(ix,jy,kz,lzhl) = 0.0
        endif

      else
       if ( lvol(lhl) .gt. 1 ) then  ! check density
        if ( an(ix,jy,kz,lvhl) .gt. 0.0 ) then
         tmp = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
        else ! in case volume is zero but mass is above threshold (should not happen, of course)
          tmp = rho_qhl
          an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
        endif

        if (  tmp .lt. xdnmn(lhl) ) then
          tmp = Max( xdnmn(lhl), tmp )
          an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
        endif

        if ( tmp .gt. xdnmx(lhl) .and. lhlw .le. 0 ) then ! no liquid allowed on hail
          tmp = Min( xdnmx(lhl), tmp )
          an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
        elseif ( tmp .gt. xdnmx(lhl) .and. lhlw .gt. 1 ) then  ! allow for liquid on hail
          fw = an(ix,jy,kz,lhlw)/an(ix,jy,kz,lhl)
!          tmpmx = xdnmx(lhl) + fw*(xdnmx(lr) - xdnmx(lhl)) ! maximum possible average density
                                                           ! it is not exactly linear, but approx. is close enough for this
!          tmpmx = 1./( (1. - fw)/900. + fw/1000. ) is exact max, where 900 is xdnmx

          tmpmx = xdnmx(lhl)/( 1. - fw*(1. - xdnmx(lhl)/xdnmx(lr) )) 

          if ( tmp .gt. tmpmx  ) then
            an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmpmx
          endif

!          IF ( tmp .gt. xdnmx(lhl) .and. an(ix,jy,kz,lhlw) .lt. qxmin(lhl) ) THEN
!            tmp = Min( xdnmx(lhl), tmp )
!            an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
!          ELSEIF ( tmp .gt. xdnmx(lr) ) THEN
!            tmp =  xdnmx(lr)
!            an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
!          ENDIF
        endif

        if ( lhlw .gt. 1 ) then ! check if basically pure water
          if ( an(ix,jy,kz,lhlw) .gt. 0.98*an(ix,jy,kz,lhl) ) then
           tmp = xdnmx(lr)
           an(ix,jy,kz,lvhl) = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/tmp
          endif
        endif
        
       endif
       
       
!  CHECK INTERCEPT
       if ( ipconc == 5 .and.  an(ix,jy,kz,lhl) .gt. qxmin(lhl) .and.  alphahl .le. 0.1 .and. lnhl .gt. 1 .and. lzhl == 0 ) then
       
         if ( lvhl .gt. 1 ) then
           hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lhl)/an(ix,jy,kz,lvhl)
         else
           hwdn = xdn0(lhl)
         endif
           tmp = (hwdn*an(ix,jy,kz,lnhl))/(dn(ix,jy,kz)*an(ix,jy,kz,lhl))
           tmpg = an(ix,jy,kz,lnhl)*(tmp*(3.14159))**(1./3.)
           if ( tmpg .lt. cnohlmn ) then
             tmp = ( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lhl))*(3.14159))**(1./3.)
              an(ix,jy,kz,lnhl) = (cnohlmn/tmp)**(3./4.)
           endif
       
       endif
!      ELSE  ! check mean size here?

      end if

      endif !lhl



      if ( lzh .gt. 1 ) then

        an(ix,jy,kz,lzh) = Max(0.0, an(ix,jy,kz,lzh) )
        
        if ( .false. .and. an(ix,jy,kz,lh) .ge. frac*qxmin(lh) .and. rescale_low_alpha ) then
          
          if ( an(ix,jy,kz,lnh) .gt. 0.0 ) then

           if ( lvh .gt. 1 ) then
             if ( an(ix,jy,kz,lvh) .gt. 0.0 ) then
               hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
             else
               hwdn = xdn0(lh)
             endif
             hwdn = Max( xdnmn(lh), hwdn )
           else
             hwdn = xdn0(lh)
           endif

             chw = an(ix,jy,kz,lnh)
             g1 = (6.0+alphamin)*(5.0+alphamin)*(4.0+alphamin)/   &
     &            ((3.0+alphamin)*(2.0+alphamin)*(1.0+alphamin))
             z1 = g1*dn(ix,jy,kz)**2*( an(ix,jy,kz,lh) )*an(ix,jy,kz,lh)/chw
             z1  = z1*(6./(pi*hwdn))**2
          else
             z1 = 0.0
          endif
          
          an(ix,jy,kz,lzh) = Min( z1, an(ix,jy,kz,lzh) )
          
          if (  an(ix,jy,kz,lnh) .lt. 1.e-5 ) then
!            an(ix,jy,kz,lzh) = 0.9*an(ix,jy,kz,lzh)
          endif
        endif
        
      endif

      if ( an(ix,jy,kz,lh) .lt. frac*qxmin(lh) .or. zerocx(lh) ) then

!        IF ( an(ix,jy,kz,lh) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lh)
          an(ix,jy,kz,lh) = 0.0
!        ENDIF

        if ( ipconc .ge. 5 ) then ! .and. an(ix,jy,kz,lnh) .gt. 0.0 ) then
          an(ix,jy,kz,lnh) = 0.0
        endif

        if ( lvh .gt. 1 ) then
           an(ix,jy,kz,lvh) = 0.0
        endif
      
        if ( lhw .gt. 1 ) then
           an(ix,jy,kz,lhw) = 0.0
        endif

        if ( lnhf .gt. 1 ) then
           an(ix,jy,kz,lnhf) = 0.0
        endif
      
        if ( lzh .gt. 1 ) then
           an(ix,jy,kz,lzh) = 0.0
        endif

      else
       if ( lvol(lh) .gt. 1 ) then  ! check density
        if ( an(ix,jy,kz,lvh) .gt. 0.0 ) then
         tmp = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
        else
         tmp = rho_qh
          an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
        endif

        if (  tmp .lt. xdnmn(lh) ) then
          tmp = Max( xdnmn(lh), tmp )
          an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
        endif

        if ( tmp .gt. xdnmx(lh) .and. lhw .le. 0 ) then ! no liquid allowed on graupel
          tmp = Min( xdnmx(lh), tmp )
          an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
        elseif ( tmp .gt. xdnmx(lh) .and. lhw .gt. 1 ) then  ! allow for liquid on graupel
          fw = an(ix,jy,kz,lhw)/an(ix,jy,kz,lh)
!          tmpmx = xdnmx(lh) + fw*(xdnmx(lr) - xdnmx(lh)) ! maximum possible average density
                                                           ! it is not exactly linear, but approx. is close enough for this
!          tmpmx = 1./( (1. - fw)/900. + fw/1000. ) is exact max, where 900 is xdnmx
          tmpmx = xdnmx(lh)/( 1. - fw*(1. - xdnmx(lh)/xdnmx(lr) )) 

          if ( tmp .gt. tmpmx  ) then
            an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmpmx
          endif

!          IF ( tmp .gt. xdnmx(lh) .and. an(ix,jy,kz,lhw) .lt. qxmin(lh) ) THEN
!            tmp = Min( xdnmx(lh), tmp )
!            an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
!          ELSEIF ( tmp .gt. xdnmx(lr) ) THEN
!            tmp =  xdnmx(lr)
!            an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
!          ENDIF

        endif

        if ( lhw .gt. 1 ) then ! check if basically pure water
          if ( an(ix,jy,kz,lhw) .gt. 0.98*an(ix,jy,kz,lh) ) then
           tmp = xdnmx(lr)
           an(ix,jy,kz,lvh) = dn(ix,jy,kz)*an(ix,jy,kz,lh)/tmp
          endif
        endif
        
       endif

!  CHECK INTERCEPT
       if ( ipconc == 5 .and.  an(ix,jy,kz,lh) .gt. qxmin(lh) .and.  alphah .le. 0.1 .and. lnh .gt. 1 .and. lzh == 0 ) then
       
         if ( lvh .gt. 1 ) then
           if ( an(ix,jy,kz,lvh) .gt. 0.0 ) then
             hwdn = dn(ix,jy,kz)*an(ix,jy,kz,lh)/an(ix,jy,kz,lvh)
           else
             hwdn = xdn0(lh)
           endif
           hwdn = Max( xdnmn(lh), hwdn )
         else
           hwdn = xdn0(lh)
         endif
           tmp = (hwdn*an(ix,jy,kz,lnh))/(dn(ix,jy,kz)*an(ix,jy,kz,lh))
           tmpg = an(ix,jy,kz,lnh)*(tmp*(3.14159))**(1./3.)
           if ( tmpg .lt. cnohmn ) then
!           tmpg = an(ix,jy,kz,lnh)*( (hwdn*an(ix,jy,kz,lnh))/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
!           tmpg = an(ix,jy,kz,lnh)**(4./3.)*( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
             tmp = ( (hwdn)/(dn(ix,jy,kz)*an(ix,jy,kz,lh))*(3.14159))**(1./3.)
              an(ix,jy,kz,lnh) = (cnohmn/tmp)**(3./4.)
           endif
       
       endif
        
      end if


      if ( an(ix,jy,kz,ls) .lt.  frac*qxmin(ls)  .or. zerocx(ls)  & ! .or.  an(ix,jy,kz,lns) .lt. 0.1 ! .and.
     &         ) then
      if ( t0(ix,jy,kz) .lt. 273.15 ) then
!        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
!        ENDIF
      
        if ( ipconc .ge. 4 ) then ! .and. an(ix,jy,kz,lns) .gt. 0.0  ) then ! 
!          an(ix,jy,kz,lni) = an(ix,jy,kz,lni) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        endif
        
        if ( lvs .gt. 1 ) then
           an(ix,jy,kz,lvs) = 0.0
        endif

        if ( lsw .gt. 1 ) then
           an(ix,jy,kz,lsw) = 0.0
        endif

      else
!        IF ( an(ix,jy,kz,ls) .gt. 0 ) THEN
          an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,ls)
          an(ix,jy,kz,ls) = 0.0
!        ENDIF

        if ( lvs .gt. 1 ) then
           an(ix,jy,kz,lvs) = 0.0
        endif

        if ( lsw .gt. 1 ) then
           an(ix,jy,kz,lsw) = 0.0
        endif

        if ( ipconc .ge. 4 ) then ! .and. an(ix,jy,kz,lns) .gt. 0.0  ) then ! 
!          an(ix,jy,kz,lnr) = an(ix,jy,kz,lnr) + an(ix,jy,kz,lns)
          an(ix,jy,kz,lns) = 0.0
        endif

      endif
      

      elseif ( lvol(ls) .gt. 1 ) then  ! check density
        if ( an(ix,jy,kz,lvs) .gt. 0.0 ) then
          tmp = dn(ix,jy,kz)*an(ix,jy,kz,ls)/an(ix,jy,kz,lvs)
          if ( tmp .gt. xdnmx(ls) .or. tmp .lt. xdnmn(ls) ) then
            tmp = Min( xdnmx(ls), Max( xdnmn(ls), tmp ) )
            an(ix,jy,kz,lvs) = dn(ix,jy,kz)*an(ix,jy,kz,ls)/tmp
          endif
        else
          tmp = rho_qs
          an(ix,jy,kz,lvs) = dn(ix,jy,kz)*an(ix,jy,kz,ls)/tmp
        endif


      end if

        if ( lzr > 1 ) then
          an(ix,jy,kz,lzr) = Max(0.0, an(ix,jy,kz,lzr) )
        endif

      if ( an(ix,jy,kz,lr) .lt. frac*qxmin(lr)  .or. zerocx(lr)  &
     &  ) then
        an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lr)
        an(ix,jy,kz,lr) = 0.0
        if ( ipconc .ge. 3 ) then
!          an(ix,jy,kz,lnc) = an(ix,jy,kz,lnc) + an(ix,jy,kz,lnr)
          an(ix,jy,kz,lnr) = 0.0
        endif
        
        if ( lzr > 1 ) then
          an(ix,jy,kz,lzr) = 0.0
        endif

      end if

!
!  for qci
!
      if ( an(ix,jy,kz,li) .le. frac*qxmin(li) .or. zerocx(li)   & ! .or.  an(ix,jy,kz,lni) .lt. 0.1
     &    ) then
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,li)
      an(ix,jy,kz,li)= 0.0
       if ( ipconc .ge. 1 ) then
         an(ix,jy,kz,lni) = 0.0
       endif
      endif

!
!  for qis
!
      if ( lis > 1 ) then ! {
      if ( an(ix,jy,kz,lis) .le. frac*qxmin(lis) .or. zerocx(lis)   & ! .or.  an(ix,jy,kz,lni) .lt. 0.1
     &    ) then ! { {
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lis)
      an(ix,jy,kz,lis)= 0.0
       if ( ipconc .ge. 1 ) then
         an(ix,jy,kz,lnis) = 0.0
       endif
      
      elseif ( icespheres >= 2 ) then ! } {
       km1 = Max(1, kz-1)
       if ( 0.5*( w(ix,jy,kz) + w(ix,jy,kz+1)) < -1.0 .or.    &
     &      (icespheres == 3 .and. ( t0(ix,jy,kz) < 232.15 .or. an(ix,jy,kz,lc) < qxmin(lc) ) ) .or. &
     &      (icespheres == 5 .and. ( t0(ix,jy,kz) < 232.15 .or. &
     &         ( an(ix,jy,kz,lc) < qxmin(lc) .and. an(ix,jy,km1,lc) < qxmin(lc)  )) ) .or. &
     &      (icespheres == 4 .and. ( t0(ix,jy,kz) < 235.15 )) ) then ! transfer to regular ice crystals in downdraft or at low temp
         an(ix,jy,kz,li) = an(ix,jy,kz,li) + an(ix,jy,kz,lis)
         an(ix,jy,kz,lni) = an(ix,jy,kz,lni) + an(ix,jy,kz,lnis)
         an(ix,jy,kz,lis)= 0.0
         an(ix,jy,kz,lnis)= 0.0
         
       endif
       
      endif ! } }
      endif ! }

!
!  for qcw
!

      if ( an(ix,jy,kz,lc) .le. frac*qxmin(lc) .or. zerocx(lc)   &
     &       ) then
      an(ix,jy,kz,lv) = an(ix,jy,kz,lv) + an(ix,jy,kz,lc)
      an(ix,jy,kz,lc)= 0.0
       if ( ipconc .ge. 2 ) then
        if ( lccn .gt. 1 .or. ac_opt == 1 ) then
          if ( irenuc < 5 .and. lccna <= 1 ) then
            if ( ac_opt == 0 ) then
               an(ix,jy,kz,lccn) = an(ix,jy,kz,lccn) + Max(0.0,an(ix,jy,kz,lnc))
            endif
          elseif ( lccna > 1 ) then
            an(ix,jy,kz,lccna) = Max( 0.0, an(ix,jy,kz,lccna) - Max(0.0,an(ix,jy,kz,lnc)) )
          endif
        endif
         an(ix,jy,kz,lnc) = 0.0
         if ( lccn > 1 ) an(ix,jy,kz,lccn) = Max( 0.0, an(ix,jy,kz,lccn) )
         
         if ( lccna > 0 .and. ac_opt == 0  ) then ! apply exponential decay to activated CCN to restore to environmental value
           if ( restoreccn ) then
           tmp = an(ix,jy,kz,li) + an(ix,jy,kz,ls)  
           
           if ( an(ix,jy,kz,lccna) > 1. .and. tmp < qxmin(li) ) an(ix,jy,kz,lccna) = an(ix,jy,kz,lccna)*Exp(-dtp/ccntimeconst)
           endif
         elseif ( lccn > 1 .and. restoreccn .and. ac_opt == 0  ) then
           ! in this case, we are treating the ccn field as ccna
           tmp = an(ix,jy,kz,li) + an(ix,jy,kz,ls)  
!           IF ( ny == 2 .and. ix == nx/2 ) THEN
!             write(0,*) 'restore: k, qccn,exp = ',kz,qccn,dn(ix,jy,kz)*qccn,Exp(-dtp/ccntimeconst)
!             write(0,*) 'ccn1,ccn2 = ',an(ix,jy,kz,lccn),dn(ix,jy,kz)*qccn - Max(0.0 , dn(ix,jy,kz)*qccn - an(ix,jy,kz,lccn))*Exp(-dtp/ccntimeconst)
!           ENDIF
           if ( an(ix,jy,kz,lccn) > 1. .and. tmp < qxmin(li) .and. ( an(ix,jy,kz,lccn) < dn(ix,jy,kz)*qccn .or. .not. invertccn ) ) then 
        !      an(ix,jy,kz,lccn) =  &
        !            an(ix,jy,kz,lccn) +  Max(0.0 , dn(ix,jy,kz)*qccn - an(ix,jy,kz,lccn))*(1.0 - Exp(-dtp/ccntimeconst))
        ! Equivalent form after expanding last term:
               an(ix,jy,kz,lccn) =  &
                    dn(ix,jy,kz)*qccn - Max(0.0 , dn(ix,jy,kz)*qccn - an(ix,jy,kz,lccn))*Exp(-dtp/ccntimeconst)
           endif
         
         endif

       endif

      endif

      end do
!      end do
      end do
      
      endif ! true/false
      
      if ( ndebug .ge. 1 ) write(6,*) 'end OF ICEZVD_DR'
!
!
   
   
   9999 return
   
   end subroutine NUCOND


! #####################################################################
! #####################################################################




!c--------------------------------------------------------------------------
!
!
!--------------------------------------------------------------------------
!

      subroutine nssl_2mom_gs   &
     &  (nx,ny,nz,na,jyslab  &
     &  ,nor,norz          &
     &  ,dtp,gz       &
     &  ,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9      &
     &  ,an,dn,p2                  &
     &  ,pn,w,iunit                   &
     &  ,t00,t77,                             &
     &   ventr,ventc,c1sw,jgs,ido,    &
     &   xdnmx,xdnmn,               &
!     &   ln,ipc,lvol,lz,lliq,   &
     &   cdx,                              &
     &   xdn0,tmp3d,tkediss  &
     &  ,thproc,numproc,dx1,dy1,ngs     &
     & ,timevtcalc,axtra,io_flag  &
     & , has_wetscav,rainprod2d, evapprod2d, alpha2d &
     & ,elec,its,ids,ide,jds,jde &
     & )


!
!--------------------------------------------------------------------------
!                                
!     Ziegler 1985 parameterized microphysics (also Zrnic et al. 1993)
!     1)  cloud water
!     2)  rain
!     3)  column ice 
!     6)  snow
!     11) graupel/hail
!
!--------------------------------------------------------------------------
!
! Notes:
!
!  4/27/2009: allows for liquid water to be advected on snow and graupel particles using flag "mixedphase"
!
!  3/14/2007: (APS) added qproc temp to make microphysic process timeseries
!
!  10/17/2006: added flag (iehw) to select how to calculate ehw
!
!  10/5/2006: switched chacr to integrated version rather than assuming that average rain
!             drop mass does not change.  This acts to reduce rain size somewhat via graupel
!             collection.
!             Use Mason data for ehw, with scaling toward ehw=1 as air density decreases.
!
!  10/3/2006: Turned off Meyers nucleation for T > -5 (can turn on with imeyers5 flag)
!             Turned off contact nucleation in updrafts
!
!  7/24/2006:  Turned on Meyers nucleation for -5 < T < 0
!
!  5/12/2006:  Converted qsacw/csacw and qsaci/csaci to Z93
!
!  5/12/2006:  Put a threshold on Bigg rain freezing.  If the frozen drops
!              have an average volume less than xvhmn, then the drops are put
!              into snow instead of graupel/hail.
!
!              Fixed bug when vapor deposition was limited.
!
!  5/13/2006:  Note that qhacr has a large effect, but Z85 did not include it.
!              Turned off qsacr (set to zero).
!
!  9/14/2007: erm: recalculate vx(lh) after setting xdn(lh) in case xdn was out of allowed range.
!             added parameter rimc3 for minimum rime density.  Default value set at 170. kg/m**3
!             instead of previous use of 100.  (Farley, 1987)
!
!--------------------------------------------------------------------------
!
!  general declarations
!
!--------------------------------------------------------------------------
!
!
!


      implicit none
!
!      integer icond 
!      parameter ( icond = 2 )

      integer, parameter :: ng1 = 1

      integer nx,ny,nz,na,nba,nv
      integer nor,norz,istag,jstag,kstag ! ,nht,ngt,igsr
      integer iwrite
      real dtp,dx,dy,dz

      logical, intent(in) :: io_flag

      integer itile,jtile,ktile
      integer ixbeg,jybeg
      integer ixend,jyend,kzend,kzbeg
      integer nxend,nyend,nzend,nzbeg
      integer :: my_rank = 0
      integer, parameter :: myprock = 1, nprock = 1
      logical, intent(in) :: has_wetscav
      integer, intent(in) :: numproc
      real, intent(inout) :: thproc(nz,numproc)
      real, intent(in) :: dx1,dy1
      real rainprod2d(-nor+1:nx+nor,-norz+ng1:nz+norz)
      real evapprod2d(-nor+1:nx+nor,-norz+ng1:nz+norz)

      real alpha2d(-nor+1:nx+nor,-norz+ng1:nz+norz,3)

      real, parameter :: tfrdry = 243.15

      logical lrescalelow(lc:lhab)
      real tkediss(-nor+1:nx+nor,-norz+ng1:nz+norz)
      real axtra(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,nxtra)

      real :: galpharaut
      real :: xvbarmax
      
      integer jyslab,its,ids,ide,jds,jde ! domain boundaries
      integer, intent(in) :: iunit !,iunit0
      real qvex
      integer iraincv, icgxconv
      parameter ( iraincv = 1, icgxconv = 1)
      real ffrz
      real :: ffrzh = 1.0

      real qcitmp,cirdiatmp ! ,qiptmp,qirtmp
      real ccwtmp,ccitmp ! ,ciptmp,cirtmp
      real cpqc,cpci ! ,cpip,cpir
      real cpqc0,cpci0 ! ,cpip0,cpir0
      real scfac ! ,cpip1
      
      double precision dp1
      
      double precision frac, frach, xvfrz, xvbiggsnow
      
      double precision :: timevtcalc
      double precision :: dpt1,dpt2
            
      logical, parameter :: gammacheck = .false.
      integer :: luindex
      double precision :: tmpgam
      logical, parameter :: usegamxinfcnu = .false.
      logical, parameter :: usegamxinf = .false.
      logical, parameter :: usegamxinf2 = .false.
      logical, parameter :: usegamxinf3 = .false.
!      real rar  ! rime accretion rate as calculated from qxacw

! a few vars for time-split fallout      
      real vtmax
      integer n,ndfall
      
      double precision chgneg,chgpos,sctot
      
      real temgtmp

      real pb(-norz+ng1:nz+norz)
      real pinit(-norz+ng1:nz+norz)

      real gz(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) ! dz
      
      real qimax,xni0,roqi0


      real dv

      real dtptmp
      integer itest,nidx,id1,jd1,kd1
      parameter (itest=1)
      parameter (nidx=10)
      parameter (id1=1,jd1=1,kd1=1)
      integer ierr
      integer iend

      integer ix,kz, il, ic, ir, icp1, irp1, ip1,jp1,kp1
      integer :: jy
      integer i,j,k,i1
      integer kzb,kze
      real slope1, slope2
      real x1, x2, x3
      real eps,eps2
      parameter (eps=1.e-20,eps2=1.e-5)
!
!  Other elec. vars
!
      real  temele
      real  trev
      
      logical ldovol, ishail, ltest, wtest
      logical , parameter :: alp0flag = .false.
!
!
!  wind indicies
!
      integer mu,mv,mw
      parameter (mu=1,mv=2,mw=3)
!
!  conversion parameters
!
      integer mqcw,mqxw,mtem,mrho,mtim
      parameter (mqcw=21,mqxw=21,mtem=21,mrho=5,mtim=6)

      real xftim,xftimi,yftim, xftem,yftem, xfqcw,yfqcw, xfqxw,yfqxw
      parameter (xftim=0.05,xftimi = 1./xftim,yftim=1.)
      parameter (xftem=0.5,yftem=1.)
      parameter (xfqcw=2000.,yfqcw=1.)
      parameter (xfqxw=2000.,yfqxw=1.)
      real dtfac
      parameter ( dtfac = 1.0 )
      integer ido(lc:lqmx)

!      integer iexy(lc:lqmx,lc:lqmx)
!      integer ieswi, ieswir, ieswip, ieswc, ieswr
!      integer ieglsw, iegli, ieglir, ieglip, ieglc, ieglr
!      integer iegmsw, iegmi, iegmir, iegmip, iegmc, iegmr
!      integer ieghsw, ieghi, ieghir, ieghip, ieghc, ieghr
!      integer iefwsw, iefwi, iefwir, iefwip, iefwc, iefwr
!      integer iehwsw, iehwi, iehwir, iehwip, iehwc, iehwr
!      integer iehlsw, iehli, iehlir, iehlip, iehlc, iehlr
!      real delqnsa, delqxsa, delqnsb, delqxsb, delqnia, delqxia
!      real delqnra, delqxra

       real delqnxa(lc:lqmx)
       real delqxxa(lc:lqmx)
!
! external temporary arrays
!
      real t00(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t77(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t1(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t2(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t3(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t4(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t5(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t6(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t7(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t8(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)
      real t9(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

      real p2(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)  ! perturbation Pi
      real pn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real an(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz,na)
      real dn(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)
      real w(-nor+1:nx+nor,-nor+1:ny+nor,-norz+ng1:nz+norz)

      real tmp3d(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz)

! 
!  declarations microphyscs and for gather/scatter
!
      integer nxmpb,nzmpb,nxz
      integer jgs,mgs,ngs,numgs
      integer, parameter :: ngsz = 500
      integer ntt
      parameter (ntt=300)

      real dvmgs(ngs)
      
      integer ngscnt,igs(ngs),kgs(ngs)
      integer kgsp(ngs),kgsm(ngs),kgsm2(ngs)
      integer ncuse
      parameter (ncuse=0)
      integer il0(ngs),il5(ngs),il2(ngs),il3(ngs)
!      integer il1m(ngs),il2m(ngs),il3m(ngs),il4m(ngs),il5m(ngs)
!
      real tdtol,temsav,tfrcbw,tfrcbi
      real, parameter :: thnuc = 235.15
!
!  Ice Multiplication Arrays.
!
      real  fimt1(ngs),fimta(ngs),fimt2(ngs) !,qmul1(ngs),qmul2(ngs)
      real xcwmas
!
!
! Variables for Ziegler warm rain microphysics
!      


      real ccnc(ngs),ccin(ngs),cina(ngs),ccna(ngs)
      real cwnccn(ngs)
      real sscb  ! 'cloud base' SS threshold
      parameter ( sscb = 2.0 )
      integer idecss  ! flag to turn on (=1) decay of ssmax when no cloud or ice crystals
      parameter ( idecss = 1 )
      integer iba ! flag to do condensation/nucleation in 1st or 2nd loop
                  ! =0 to use ad to calculate SS
                  ! =1 to use an at end of main jy loop to calculate SS
      parameter (iba = 1)
      integer ifilt   ! =1 to filter ssat, =0 to set ssfilt=ssat
      parameter ( ifilt = 0 ) 
      real temp1,temp2 ! ,ssold
      real :: mwat, mice, dice, mwshed, fwmax, fw, mwcrit, massfactor, tmpdiam
      real, parameter :: shedalp = 3.  ! set 3 for maximum mass diameter (same as area-weighted diameter), 4 for mass-weighted diameter
      real ssmax(ngs)       ! maximum SS experienced by a parcel
      real ssmx
      real dnnet,dqnet
!      real cnu,rnu,snu,cinu
!      parameter ( cnu = 0.0, rnu = -0.8, snu = -0.8, cinu = 0.0 )
      real bfnu, bfnu0, bfnu1
      parameter ( bfnu0 = (rnu + 2.0)/(rnu + 1.0)  )
      real ventr, ventc
      real volb
      double precision t2s, xdp
      double precision xl2p(ngs),rb(ngs)
      real, parameter :: aa1 = 9.44e15, aa2 = 5.78e3  ! a1 in Ziegler
! snow parameters:
      real, parameter :: cexs = 0.1, cecs = 0.5 
      real, parameter :: rvt = 0.104  ! ratio of collection kernels (Zrnic et al, 1993)
      real, parameter :: kfrag = 1.0e-6 ! rate coefficent for collisional splintering (Schuur & Rutledge 00b)
      real, parameter :: mfrag = 1.0e-10 ! assumed ice fragment mass for collisional splintering (Schuur & Rutledge 00b)
      double precision cautn(ngs), rh(ngs), nh(ngs)
      real ex1, ft, rhoinv(ngs)
      double precision ec0(ngs)
      
      real ac1,bc, taus, c1,d1,e1,f1,p380,tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,temp3 ! , sstdy, super
      real :: flim
      real dw,dwr
      double precision :: tmpz, tmpzmlt
      real ratio, delx, dely
      real dbigg,volt
      real chgtmp,fac,mixedphasefac
      real x,y,y2,del,r,rtmp,alpr
      double precision :: vent1,vent2
      double precision :: g1palp,g4palp
      double precision :: g1palpinf,g4palpinf
      real fqt !charge separation as fn of temperature from Dong and Hallett 1992
      real bs
      real v1, v2
      real d1r, d1i, d1s, e1i
      real c1sw   ! integration factor for snow melting with snu = -0.8
      real, parameter :: vr1mm = 5.23599e-10 ! volume of 1mm diameter sphere (m**3)
      real, parameter :: vr3mm = 5.23599e-10*(3.0/1.)**3   ! volume of a 3 mm diameter sphere (m**3) (Rasmussen et al. 1984b, JAS)
      real, parameter :: vr4p5mm = 5.23599e-10*(4.5/1.)**3 ! volume of 4.5mm diameter sphere (m**3) (Rasmussen et al. 1984b, JAS)
      real vmlt,vshd, vshdgs(ngs,lh:lhab), maxmassfac(lc:lhab)
      real rhosm
      parameter ( rhosm = 500. )
      integer nc ! condensation step
      real dtcon,dtcon1,dtcon2 ! condensation time step (dtcon*nc = dtp)
      real delta
      integer ltemq1,ltemq1m ! ,ltemq1m2
      real dqv,qv1,ss1,ss2,qvs1,dqvs,dtemp,dt1   ! temporaries for condensation
      real ssi1, ssi2, dqvi, dqvis, dqvii,qis1
      real dqvr, dqc, dqr, dqi, dqs
      real qv1m,qvs1m,ss1m,ssi1m,qis1m
      real cwmastmp
      real  dcloud,dcloud2 ! ,as, bs
      real cn(ngs)
      double precision xvc, xvr
      real mwfac
!      real  es(ngs) ! ss(ngs),
!      real  eis(ngs)

      real rwmasn,rwmasx

      real vgra,vfrz
      parameter ( vgra = 0.523599*(1.0e-3)**3 )
     
!      real, parameter :: epsi = 0.622
!      real, parameter :: d = 0.266
      real :: d, dold, denom,denominv,vth
      double precision :: h1, h2, h3, h4,denomdp, denominvdp
      real r1,qevap ! ,slv
      
      real vr,nrx,chw,g1,qr,z,z1,rdi,alp,xnutmp,xnuc,g1r,rd1,rdia,rmas
      real :: snowmeltmass = 0
      
!      real, parameter :: rhofrz = 900.   ! density of graupel from newly-frozen rain
      real, parameter :: rimedens = 500. ! default rime density

!      real svc(ngs)  !  droplet volume
!
!  contact freezing nucleation
!
      real raero,kaero !assumd aerosol radius, thermal conductivity
      parameter ( raero = 3.e-7, kaero = 5.39e-3 )
      real kb   ! Boltzman constant  J K-1
      parameter (kb = 1.3807e-23)
      
      real knud(ngs),knuda(ngs) !knudsen number and correction factor
      real gtp(ngs)  !G(T,p) = 1/(a' + b')  Cotton 72b
      real dfar(ngs) !aerosol diffusivity
      real fn1(ngs),fn2(ngs),fnft(ngs)
      
      real ccia(ngs)
      real ctfzbd(ngs),ctfzth(ngs),ctfzdi(ngs)
!
!  misc
!
      real ni,nis,nr,d0
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs),dqcitmp(ngs),dqwvtmp(ngs)
      real tempc(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),qvap(ngs) 
      real temgkm1(ngs), temgkm2(ngs)
      real temgx(ngs),temcgx(ngs)
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real elv(ngs),elf(ngs),els(ngs)
      real tsqr(ngs),ssi(ngs),ssw(ngs),ssi0(ngs)
      real qcwtmp(ngs),qtmp,qtot(ngs) 
      real qcond(ngs)
      real ctmp, sctmp
      real cimasn,cimasx,ccimx
      real pid4
      real cs,ds,gf7,gf6,gf5,gf4,gf3,gf2,gf1
      real gcnup1,gcnup2
      real gf73rds, gf83rds
      real gamice73fac, gamsnow73fac
      real gf43rds, gf53rds
      real aradcw,bradcw,cradcw,dradcw,cwrad,rwrad,rwradmn
      parameter ( rwradmn = 50.e-6 )
      real dh0
      real dg0(ngs),df0(ngs)
      real dhwet(ngs),dhlwet(ngs),dfwet(ngs)
      
      real clionpmx,clionnmx
      parameter (clionpmx=1.e9,clionnmx=1.e9) ! Takahashi 84
!
!  other arrays

      real fwet1(ngs),fwet2(ngs)   
      real fmlt1(ngs),fmlt2(ngs),fmlt1e(ngs)
      real fvds(ngs),fvce(ngs),fiinit(ngs) 
      real fvent(ngs),fraci(ngs),fracl(ngs)
!
      real fai(ngs),fav(ngs),fbi(ngs),fbv(ngs)
      real felv(ngs),fels(ngs),felf(ngs)
      real felvcp(ngs),felscp(ngs),felfcp(ngs)
      real felvpi(ngs),felspi(ngs),felfpi(ngs)
      real felvs(ngs),felss(ngs)      !   ,felfs(ngs)
      real fwvdf(ngs),ftka(ngs),fthdf(ngs)
      real fadvisc(ngs),fakvisc(ngs)
      real fci(ngs),fcw(ngs) ! heat capacities of ice and liquid
      real fschm(ngs),fpndl(ngs)
      real fgamw(ngs),fgams(ngs)
      real fcqv1(ngs),fcqv2(ngs),fcc3(ngs) 
      
      real cvm,cpm,rmm

      real, parameter ::      cpv = 1885.0       ! specific heat of water vapor at constant pressure
!
      real fcci(ngs), fcip(ngs)
!
      real :: sfm1(ngs),sfm2(ngs)
      real :: gfm1(ngs),gfm2(ngs)
      real :: ffm1(ngs),ffm2(ngs)
      real :: hfm1(ngs),hfm2(ngs)

      logical :: wetsfc(ngs),wetsfchl(ngs),wetsfcf(ngs)
      logical :: wetgrowth(ngs), wetgrowthhl(ngs), wetgrowthf(ngs)

      real qitmp(ngs),qistmp(ngs)
       
      real rzxh(ngs), rzxhl(ngs), rzxhlh(ngs), rzxhlf(ngs)
      real rzxs(ngs), rzxf(ngs)
!      real axh(ngs),bxh(ngs),axhl(ngs),bxhl(ngs)
      real cdh(ngs),cdhl(ngs)
      real :: axx(ngs,lh:lhab),bxx(ngs,lh:lhab)
      real vt2ave(ngs)

      real :: qcwresv(ngs), ccwresv(ngs) ! "reserved" droplet mass and number that are too small for accretion
      
      real ::  lfsave(ngs,6)
      real ::  qx(ngs,lv:lhab)
      real ::  qxw(ngs,ls:lhab)
      real ::  qxwlg(ngs,lh:lhab)
      real ::  chxf(ngs,lh:lhab)
      real ::  cx(ngs,lc:lhab)
      real ::  cxmxd(ngs,lc:lhab)
      real ::  qxmxd(ngs,lv:lhab)
      real ::  scx(ngs,lc:lhab)
      real ::  xv(ngs,lc:lhab)
      real ::  vtxbar(ngs,lc:lhab,3)
      real ::  xmas(ngs,lc:lhab)
      real ::  xdn(ngs,lc:lhab)
      real ::  xdntmp(ngs,lc:lhab)
      real ::  cdxgs(ngs,lc:lhab)
      real ::  xdia(ngs,lc:lhab,3)
      real ::  vtwtdia(ngs,lr:lhab) ! sweep-out volume weighted diameter
      real ::  rarx(ngs,ls:lhab)
      real ::  vx(ngs,li:lhab)
      real ::  rimdn(ngs,li:lhab)
      real ::  raindn(ngs,li:lhab)
      real ::  alpha(ngs,lc:lhab)
      real ::  dab0lh(ngs,lc:lhab,lc:lhab)
      real ::  dab1lh(ngs,lc:lhab,lc:lhab)
      real ::  zx(ngs,lr:lhab)
      real ::  zxmxd(ngs,lr:lhab)
      real ::  g1x(ngs,lr:lhab)
      

      real :: qsimxdep(ngs) ! max sublimation of qi+qs+qis
      real :: qsimxsub(ngs) ! max depositionof qi+qs+qis
      logical,parameter :: DoSublimationFix = .true.
      real :: qrtmp(ngs),qvtmp(ngs),qctmp(ngs)
      real :: felvcptmp,felscptmp,qsstmp
      real :: thetatmp, thetaptmp, temcgtmp,qvaptmp
      real :: qvstmp, qisstmp, qvptmp, qitmp1, qctmp1
      
      real :: galphrout
      
      real ventrx(ngs)
      real ventrxn(ngs)
      real g1shr, alphashr
      real g1mlr, alphamlr
      real g1smlr, alphasmlr
      real massfacshr, massfacmlr
      
      real :: qhgt8mm ! ice mass greater than 8mm
      real :: qhwgt8mm ! ice + max water mass greater than 8mm
      real :: qhgt10mm ! mass greater than 10mm
      real :: qhgt20mm ! mass greater than 20mm
      real :: fwmhtmp
      real, parameter :: fwmhtmptem = -15. ! temperature at which fwmhtmp fully switches to liquid water only being on large particles
      real, parameter :: d1t = (6.0 * 0.268e-3/(917.* pi))**(1./3.) ! d1t is the diameter of the ice sphere with the mass (0.268e-3 kg) of an 8mm spherical drop
      real, parameter :: srasheym = 0.1389 ! slope fraction from Rasmussen and Heymsfield 
!
      real swvent(ngs),hwvent(ngs),rwvent(ngs),hlvent(ngs),hwventy(ngs),hlventy(ngs),rwventz(ngs)
      real hxventtmp
      real hlventinc(ngs),hwventinc(ngs)
      integer, parameter :: ndiam = 10
      integer :: numdiam
      real hwvent0(ndiam+4),hlvent0 ! 0 to d1
      real hwvent1,hlvent1 ! d1 to infinity
      real hwvent2,hlvent2 ! d2 to infinity
      real gama0,gamb0
      real gama1,gamb1
      real gama2,gamb2
!      real, parameter :: mltdiam1 = 9.0e-3, mltdiam1p5 = 16.0e-3, mltdiam2 = 19.0e-3, mltdiam3 = 200.0e-3, mltdiam05 = 4.5e-3
      real :: mltdiam(ndiam+4)
      real mltmass0inv,mltmass1inv,mltmass2inv, mltmass1cgs, mltmass2cgs,mltmass3inv, mltmass3cgs
      real qhmlr0, qhmlr05, qhmlr1, qhmlr2,qhmlr3, qhmlr12, qhmlr23
      real qhlmlr0, qhlmlr05, qhlmlr1, qhlmlr2,qhlmlr3, qhlmlr12, qhlmlr23
      real qxd1, cxd1, zxd1 ! mass and number up to mltdiam1
      real qxd05, cxd05 ! mass and number up to mltdiam1/2
      
      real :: qxd(ndiam+4), cxd(ndiam+4), qhml(ndiam+4), qhml0(ndiam+4)
      real :: dqxd(ndiam+4), dcxd(ndiam+4), dqhml(ndiam+4)
      
      
      real civent(ngs)
      real isvent(ngs)
!
      real xmascw(ngs)
      real xdnmx(lc:lhab), xdnmn(lc:lhab)
      real dnmx
      real :: xdiamxmas(ngs,lc:lhab)
!
      real cilen(ngs) ! ,ciplen(ngs)
!
!
      real rwcap(ngs),swcap(ngs)
      real hwcap(ngs)
      real hlcap(ngs)
      real cicap(ngs)
      real iscap(ngs)

      real qvimxd(ngs)
      real qimxd(ngs),qismxd(ngs),qcmxd(ngs),qrmxd(ngs),qsmxd(ngs),qhmxd(ngs),qhlmxd(ngs)
      real cimxd(ngs),ccmxd(ngs),crmxd(ngs),csmxd(ngs),chmxd(ngs)
      real cionpmxd(ngs),cionnmxd(ngs)
      real clionpmxd(ngs),clionnmxd(ngs)


      real elec(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-norz+ng1:nz+norz) ! Ez (elecsave)

!
!
      ! Hallett-Mossop arrays
      real chmul1(ngs),chlmul1(ngs),csmul1(ngs),csmul(ngs)
      real qhmul1(ngs),qhlmul1(ngs),qsmul1(ngs),qsmul(ngs)
      
      ! splinters from drop freezing
      real csplinter(ngs),qsplinter(ngs)
      real csplinter2(ngs),qsplinter2(ngs)
!
!
!  concentration arrays...
!
      real :: chlcnh(ngs), vhlcnh(ngs), vhlcnhl(ngs)
      real :: chlcnhhl(ngs) ! number of new hail particles (may be different from number of lost graupel)
      real cracif(ngs), ciacrf(ngs)
      real cracr(ngs)

!
      real ciint(ngs), crfrz(ngs), crfrzf(ngs), crfrzs(ngs)
      real cicint(ngs)
      real cipint(ngs)
      real ciacw(ngs), cwacii(ngs) 
      real ciacr(ngs), craci(ngs)
      real csacw(ngs)
      real csacr(ngs)
      real csaci(ngs),   csacs(ngs)
      real cracw(ngs) 
      real chacw(ngs), chacr(ngs)
      real :: chlacw(ngs) 
      real chaci(ngs), chacs(ngs)
!
      real :: chlacr(ngs)
      real :: chlaci(ngs), chlacs(ngs)
      real crcnw(ngs) 
      real cidpv(ngs),cisbv(ngs)
      real cisdpv(ngs),cissbv(ngs)
      real cimlr(ngs),cismlr(ngs)

      real chlsbv(ngs), chldpv(ngs)
      real chlmlr(ngs), chlmlrr(ngs)
      real chlfmlr(ngs)
!      real chlmlrsave(ngs),chlsave(ngs),qhlsave(ngs)
      real chlshr(ngs), chlshrr(ngs)


      real chdpv(ngs),chsbv(ngs)
      real chmlr(ngs),chcev(ngs)
      real chmlrr(ngs)
      real chshr(ngs), chshrr(ngs)

      real csdpv(ngs),cssbv(ngs)
      real csmlr(ngs),csmlrr(ngs),cscev(ngs)
      real csshr(ngs), csshrr(ngs)

      real crcev(ngs)
      real crshr(ngs)
      real cwshw(ngs), qwshw(ngs)
!
!
! arrays for w-ac-x ;  x-ac-w
!
!
!
      real qrcnw(ngs), qwcnr(ngs)
      real zrcnw(ngs),zracr(ngs),zracw(ngs),zrcev(ngs)

      real qracw(ngs) ! qwacr(ngs),
      real qiacw(ngs) !, qwaci(ngs)

      real qsacw(ngs) ! ,qwacs(ngs),
      real qhacw(ngs) ! qwach(ngs),
      real :: qhlacw(ngs), qxacwtmp, qxacrtmp, qxacitmp, qxacstmp !
      real vhacw(ngs), vsacw(ngs), vhlacw(ngs), vhlacr(ngs)

      real qfcev(ngs)
      real qfmul1(ngs),cfmul1(ngs)
!
      real qsacws(ngs)

!
!  arrays for x-ac-r and r-ac-x; 
!
      real qsacr(ngs),qracs(ngs)
      real qhacr(ngs),qhacrmlr(ngs),qhacwmlr(ngs) ! ,qrach(ngs)
      real vhacr(ngs), zhacr(ngs), zhacrf(ngs), zrach(ngs), zrachl(ngs)
      real qiacr(ngs),qraci(ngs)
      
      real ziacr(ngs)

      real qracif(ngs),qiacrf(ngs),qiacrs(ngs),ciacrs(ngs)

      real :: qhlacr(ngs),qhlacrmlr(ngs), qhlacwmlr(ngs)
      real qsacrs(ngs) !,qracss(ngs)
!
!  ice - ice interactions
!
      real qsaci(ngs)
      real qsacis(ngs)
      real qhaci(ngs)
      real qhacs(ngs)

      real :: qhacis(ngs) 
      real :: chacis(ngs) 
      real :: chacis0(ngs)

      real :: csaci0(ngs) ! collision rate only
      real :: chaci0(ngs) ! collision rate only
      real :: chacs0(ngs) ! collision rate only
      real :: chlaci0(ngs)
      real :: chlacis(ngs)
      real :: chlacis0(ngs)
      real :: chlacs0(ngs) 

      real :: qsaci0(ngs) ! collision rate only
      real :: qsacis0(ngs) ! collision rate only
      real :: qhaci0(ngs) ! collision rate only
      real :: qhacis0(ngs) ! collision rate only
      real :: qhacs0(ngs) ! collision rate only
      real :: qhlaci0(ngs)  
      real :: qhlacis0(ngs)
      real :: qhlacs0(ngs) 

      real :: qhlaci(ngs)  
      real :: qhlacis(ngs)
      real :: qhlacs(ngs)
!
!  conversions
!
      real qrfrz(ngs) ! , qirirhr(ngs)
      real zrfrz(ngs), zrfrzf(ngs), zrfrzs(ngs)
      real ziacrf(ngs), zhcnsh(ngs), zhcnih(ngs)
      real zhacw(ngs), zhacs(ngs), zhaci(ngs)
      real zhmlr(ngs), zhdsv(ngs), zhsbv(ngs), zhlcnh(ngs), zhshr(ngs)
      real zfacw(ngs), zfacs(ngs), zfaci(ngs)
      real zfmlr(ngs), zfdsv(ngs), zfsbv(ngs), zhlcnf(ngs), zfshr(ngs), zfshrr(ngs)
      real zhmlrtmp,zhmlr0inf,zhlmlr0inf
      real zhmlrr(ngs),zhlmlrr(ngs),zhshrr(ngs),zhlshrr(ngs),zfmlrr(ngs)
!      real zsmlr(ngs)
      real zsmlrr(ngs), zsshr(ngs), zsshrr(ngs)
      real zhcns(ngs), zhcni(ngs)
      real zhwdn(ngs), zfwdn(ngs) ! change in Z due to density changes
      real zhldn(ngs) ! change in Z due to density changes

      real zhlacw(ngs), zhlacs(ngs), zhlacr(ngs)
      real zhlmlr(ngs), zhldsv(ngs), zhlsbv(ngs), zhlshr(ngs)

      
      real vrfrzf(ngs), viacrf(ngs)
      real qrfrzs(ngs), qrfrzf(ngs)
      real qwfrz(ngs), qwctfz(ngs)
      real cwfrz(ngs), cwctfz(ngs)
      real qwfrzis(ngs), qwctfzis(ngs) ! droplet freezing to ice spheres
      real cwfrzis(ngs), cwctfzis(ngs)
      real qwfrzc(ngs), qwctfzc(ngs) ! droplet freezing to columns
      real cwfrzc(ngs), cwctfzc(ngs)
      real qwfrzp(ngs), qwctfzp(ngs) ! droplet freezing to plates
      real cwfrzp(ngs), cwctfzp(ngs)
      real xcolmn(ngs), xplate(ngs)
      real ciihr(ngs), qiihr(ngs)
      real cicichr(ngs), qicichr(ngs)
      real cipiphr(ngs), qipiphr(ngs)
      real qscni(ngs), cscni(ngs), cscnis(ngs)
      real qscnvi(ngs), cscnvi(ngs), cscnvis(ngs)
      real qhcns(ngs), chcns(ngs), chcnsh(ngs), vhcns(ngs)
      real qscnh(ngs), cscnh(ngs), vscnh(ngs)
      real qhcni(ngs), chcni(ngs), chcnih(ngs), vhcni(ngs)
      real qiint(ngs),qipipnt(ngs),qicicnt(ngs)
      real cninm(ngs),cnina(ngs),cninp(ngs),wvel(ngs),wvelkm1(ngs)
      real tke(ngs)
      real uvel(ngs),vvel(ngs)
!
      real qidpv(ngs),qisbv(ngs) ! qicnv(ngs),qievv(ngs),
      real qimlr(ngs),qidsv(ngs),qisdsv(ngs),qidsvp(ngs) ! ,qicev(ngs)
      real qismlr(ngs)

!
!
      real :: qhldpv(ngs), qhlsbv(ngs) ! qhlcnv(ngs),qhlevv(ngs),
      real :: qhlmlr(ngs), qhldsv(ngs), qhlmlrsave(ngs)
      real :: qhlwet(ngs), qhldry(ngs), qhlshr(ngs), qxwettmp
!
      real :: qrfz(ngs),qsfz(ngs),qhfz(ngs),qhlfz(ngs)
      real :: qffz(ngs)
!
      real qhdpv(ngs),qhsbv(ngs) ! qhcnv(ngs),qhevv(ngs),
      real qhmlr(ngs),qhdsv(ngs),qhcev(ngs),qhcndv(ngs),qhevv(ngs)
      real qhlcev(ngs), chlcev(ngs)
      real qhwet(ngs),qhdry(ngs),qhshr(ngs)
      real qhshrp(ngs)
      real qhshh(ngs) !accreted water that remains on graupel
      real qhmlh(ngs) !melt water that remains on graupel
      real qhfzh(ngs) !water that freezes on mixed-phase graupel
      real qffzf(ngs) !water that freezes on mixed-phase FD
      real qhlfzhl(ngs) !water that freezes on mixed-phase hail
      
      real qhmlrlg(ngs),qhlmlrlg(ngs) ! melting from the larger diameters
      real qhfzhlg(ngs) !water that freezes on mixed-phase graupel (large sizes)
      real qhlfzhllg(ngs) !water that freezes on mixed-phase hail (large sizes)
      real qhlcevlg(ngs), chlcevlg(ngs)
      real qhcevlg(ngs), chcevlg(ngs)

      real vhfzh(ngs), vffzf(ngs) ! change in volume from water that freezes on mixed-phase graupel, frozen drops
      real vhlfzhl(ngs) !  change in volume from water that freezes on mixed-phase hail

      real vhshdr(ngs) !accreted water that leaves on graupel (mixedphase)
      real vhlshdr(ngs) !accreted water that leaves on hail (mixedphase)
      real vhmlr(ngs) !melt water that leaves graupel (single phase)
      real vhlmlr(ngs) !melt water that leaves hail (single phase)
      real vhsoak(ngs) !  aquired water that seeps into graupel.
      real vhlsoak(ngs) !  aquired water that seeps into hail.
      
!
      real qsdpv(ngs),qssbv(ngs) ! qscnv(ngs),qsevv(ngs),
      real qsmlr(ngs),qsdsv(ngs),qscev(ngs),qscndv(ngs),qsevv(ngs)
      real qswet(ngs),qsdry(ngs),qsshr(ngs)
      real qsshrp(ngs)
      real qsfzs(ngs)
!
!
      real qipdpv(ngs),qipsbv(ngs)
      real qipmlr(ngs),qipdsv(ngs)
!
      real qirdpv(ngs),qirsbv(ngs)
      real qirmlr(ngs),qirdsv(ngs),qirmlw(ngs)
!
      real qgldpv(ngs),qglsbv(ngs)
      real qglmlr(ngs),qgldsv(ngs)
      real qglwet(ngs),qgldry(ngs),qglshr(ngs)
      real qglshrp(ngs)
!
      real qgmdpv(ngs),qgmsbv(ngs)
      real qgmmlr(ngs),qgmdsv(ngs)
      real qgmwet(ngs),qgmdry(ngs),qgmshr(ngs)
      real qgmshrp(ngs)
      real qghdpv(ngs),qghsbv(ngs)
      real qghmlr(ngs),qghdsv(ngs) 
      real qghwet(ngs),qghdry(ngs),qghshr(ngs)
      real qghshrp(ngs)
!
      real qrztot(ngs),qrzmax(ngs),qrzfac(ngs)
      real qrcev(ngs)
      real qrshr(ngs)
      real fsw(ngs),fhw(ngs),fhlw(ngs),ffw(ngs) !liquid water fractions
      real fswmax(ngs),fhwmax(ngs),fhlwmax(ngs) !liquid water fractions
      real ffwmax(ngs)
      real qhcnf(ngs) 
      real :: qhlcnh(ngs)
      real qhcngh(ngs),qhcngm(ngs),qhcngl(ngs)
      
      real :: qhcnhl(ngs), chcnhl(ngs), zhcnhl(ngs), vhcnhl(ngs) ! conversion of low-density hail back to graupel

      real eiw(ngs),eii(ngs),eiri(ngs),eipir(ngs),eisw(ngs)
      real erw(ngs),esw(ngs),eglw(ngs),eghw(ngs),efw(ngs)
      real ehxw(ngs),ehlw(ngs),egmw(ngs),ehw(ngs)
      real err(ngs),esr(ngs),eglr(ngs),eghr(ngs),efr(ngs)
      real ehxr(ngs),ehlr(ngs),egmr(ngs) 
      real eri(ngs),esi(ngs),egli(ngs),eghi(ngs),efi(ngs),efis(ngs)
      real ehxi(ngs),ehli(ngs),egmi(ngs),ehi(ngs),ehis(ngs),ehlis(ngs)
      real ers(ngs),ess(ngs),egls(ngs),eghs(ngs),efs(ngs),ehs(ngs),ehsfac(ngs)
      real ehscnv(ngs)
      real ehxs(ngs),ehls(ngs),egms(ngs),egmip(ngs) 

      real ehsclsn(ngs),ehiclsn(ngs),ehisclsn(ngs)
      real efsclsn(ngs),eficlsn(ngs),efisclsn(ngs)
      real ehlsclsn(ngs),ehliclsn(ngs),ehlisclsn(ngs)
      real esiclsn(ngs)

      real :: ehs_collsn = 0.5, ehi_collsn = 1.0
      real :: efs_collsn = 0.5, efi_collsn = 1.0
      real :: ehls_collsn = 1.0, ehli_collsn = 1.0
      real :: esi_collsn = 1.0
      
      real ew(8,6)
      real cwr(8,2)  ! radius and inverse of interval
      data cwr / 2.0, 3.0, 4.0, 6.0,  8.0,  10.0, 15.0,  20.0 , & ! radius
     &           1.0, 1.0, 0.5, 0.5,  0.5,   0.2,  0.2,  1.  /   ! inverse of interval
      integer icwr(ngs), igwr(ngs), irwr(ngs), ihlr(ngs), ifwr(ngs)
      real grad(6,2) ! graupel radius and inverse of interval
      data grad / 100., 200., 300., 400., 600., 1000.,   &
     &            1.e-2,1.e-2,1.e-2,5.e-3,2.5e-3, 1.    /
!droplet radius: 2     3     4     6     8    10    15    20
      data ew /0.03, 0.07, 0.17, 0.41, 0.58, 0.69, 0.82, 0.88,  & ! 100
!     :         0.07, 0.13, 0.27, 0.48, 0.65, 0.73, 0.84, 0.91,  ! 150
     &         0.10, 0.20, 0.34, 0.58, 0.70, 0.78, 0.88, 0.92,  & ! 200
     &         0.15, 0.31, 0.44, 0.65, 0.75, 0.83, 0.96, 0.91,  & ! 300
     &         0.17, 0.37, 0.50, 0.70, 0.81, 0.87, 0.93, 0.96,  & ! 400
     &         0.17, 0.40, 0.54, 0.71, 0.83, 0.88, 0.94, 0.98,  & ! 600
     &         0.15, 0.37, 0.52, 0.74, 0.82, 0.88, 0.94, 0.98 / ! 1000
!     :         0.11, 0.34, 0.49, 0.71, 0.83, 0.88, 0.94, 0.95 / ! 1400


      real da0lr(ngs),da1lr(ngs)
      real da0lc(ngs),da1lc(ngs)
      real da0lh(ngs)
      real da0lhl(ngs)
      real da0lf(ngs)
      real :: da0lx(ngs,lr:lhab)
      
      real va0 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real vab0(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real vab1(lc:lqmx,lc:lqmx)  ! collection coefficients from Seifert 2005
      real va1 (lc:lqmx)          ! collection coefficients from Seifert 2005
      real ehip(ngs),ehlip(ngs),ehlir(ngs)
      real erir(ngs),esir(ngs),eglir(ngs),egmir(ngs),eghir(ngs)
      real efir(ngs),ehir(ngs),eirw(ngs),eirir(ngs),ehr(ngs)
      real erip(ngs),esip(ngs),eglip(ngs),eghip(ngs)
      real efip(ngs),eipi(ngs),eipw(ngs),eipip(ngs)
!
!  arrays for production terms
!
      real ptotal(ngs) ! , pqtot(ngs)
!
      real pqcwi(ngs),pqcii(ngs),pqrwi(ngs),pqisi(ngs)
      real pqswi(ngs),pqhwi(ngs),pqwvi(ngs)
      real pqgli(ngs),pqghi(ngs),pqfwi(ngs)
      real pqgmi(ngs),pqhli(ngs) ! ,pqhxi(ngs)
      real pqiri(ngs),pqipi(ngs) ! pqwai(ngs),
      real pqlwsi(ngs),pqlwhi(ngs),pqlwhli(ngs),pqlwfi(ngs)
      
      real pqlwlghi(ngs),pqlwlghli(ngs)
      real pqlwlghd(ngs),pqlwlghld(ngs)
      
      
      

      real pvhwi(ngs), pvhwd(ngs)
      real pvfwi(ngs), pvfwd(ngs)
      real pvhli(ngs), pvhld(ngs)
      real pvswi(ngs), pvswd(ngs)
!
      real pqcwd(ngs),pqcid(ngs),pqrwd(ngs),pqisd(ngs), pqcwdacc(ngs)
      real pqswd(ngs),pqhwd(ngs),pqwvd(ngs)
      real pqgld(ngs),pqghd(ngs),pqfwd(ngs)
      real pqgmd(ngs),pqhld(ngs) ! ,pqhxd(ngs)
      real pqird(ngs),pqipd(ngs) ! pqwad(ngs),
      real pqlwsd(ngs),pqlwhd(ngs),pqlwhld(ngs),pqlwfd(ngs)
!
!      real pqxii(ngs,nhab),pqxid(ngs,nhab)
!
      real  pctot(ngs)
      real  pcipi(ngs), pcipd(ngs)
      real  pciri(ngs), pcird(ngs)
      real  pccwi(ngs), pccwd(ngs), pccwdacc(ngs)
      real  pccii(ngs), pccid(ngs)
      real  pcisi(ngs), pcisd(ngs)
      real  pccin(ngs)
      real  pcrwi(ngs), pcrwd(ngs)
      real  pcswi(ngs), pcswd(ngs)
      real  pchwi(ngs), pchwd(ngs)
      real  pchli(ngs), pchld(ngs)
      real  pcfwi(ngs), pcfwd(ngs)
      real  pcgli(ngs), pcgld(ngs)
      real  pcgmi(ngs), pcgmd(ngs)
      real  pcghi(ngs), pcghd(ngs)

      real  pzrwi(ngs), pzrwd(ngs)
      real  pzhwi(ngs), pzhwd(ngs)
      real  pzfwi(ngs), pzfwd(ngs)
      real  pzhli(ngs), pzhld(ngs)
      real  pzswi(ngs), pzswd(ngs)

!
!  other arrays
!
      real dqisdt(ngs) !,advisc(ngs) !dqwsdt(ngs), ,schm(ngs),pndl(ngs)

      real qss0(ngs)

      real qsacip(ngs)
      real pres(ngs),pipert(ngs)
      real pk(ngs)
      real rho0(ngs),pi0(ngs)
      real rhovt(ngs),sqrtrhovt
      real thetap(ngs),theta0(ngs),qwvp(ngs),qv0(ngs)
      real thsave(ngs)
      real ptwfzi(ngs),ptimlw(ngs)
      real psub(ngs),pvap(ngs),pfrz(ngs),ptem(ngs),pmlt(ngs),pevap(ngs),pdep(ngs),ptem2(ngs)
      
      real cnostmp(ngs)   ! for diagnosed snow intercept
!
!  iholef = 1 to do hole filling technique version 1
!  which uses all hydrometerors to do hole filling of all hydrometeors
!  iholef = 2 to do hole filling technique version 2
!  which uses an individual hydrometeror species to do hole
!  filling of a species of a hydrometeor
!
!  iholen = interval that hole filling is done
!
      integer  iholef
      integer  iholen
      parameter (iholef = 1)
      parameter (iholen = 1)
      real  cqtotn,cqtotn1
      real  cctotn
      real  citotn
      real  crtotn
      real  cstotn
      real  cvtotn
      real  cftotn
      real  cgltotn
      real  cghtotn
      real  chtotn
      real  cqtotp,cqtotp1
      real  cctotp
      real  citotp
      real  ciptotp
      real  crtotp
      real  cstotp
      real  cvtotp
      real  cftotp
      real  chltotp
      real  cgltotp
      real  cgmtotp
      real  cghtotp
      real  chtotp
      real  cqfac
      real  ccfac
      real  cifac
      real  cipfac
      real  crfac
      real  csfac
      real  cvfac
      real  cffac
      real  cglfac
      real  cghfac
      real  chfac
      
      real ssifac, qvapor
!
!   Miscellaneous variables
!
      real, parameter :: cwmas30 = 1000.*0.523599*(2.*30.e-6)**3 ! mass of 30-micron radius droplet, for sat. adj.
      real, parameter :: cwmas20 = 1000.*0.523599*(2.*20.e-6)**3 ! mass of 20-micron radius droplet, for sat. adj.
      integer ireadqf,lrho,lqsw,lqgl,lqgm ,lqgh 
      integer lqrw
      real vt
      real arg  ! gamma is a function
      real erbnd1, fdgt1, costhe1
      real qeps
      real dyi2,dzi2,bta1,cnit,dragh,dnz00,pii ! ,cp608
      real qccrit,gf4br,gf4ds,gf4p5, gf3ds, gf1ds
      real gf1palp(ngs) ! for storing Gamma[1.0 + alphar]

      
      real xdn0(lc:lhab)
      real xdn_new,drhodt
      
      integer l ,ltemq,inumgs, idelq

      real brz,arz,temq

      real ssival,tqvcon
      real cdx(lc:lhab)
      real cnox
      real cval,aval,eval,fval,gval ,qsign,ftelwc,qconkq,elecfac,altelecfac
      real qconm,qconn,cfce15,gf8,gf4i,gf3p5,gf1a,gf1p5,qdiff,argrcnw
      real c4,bradp,bl2,bt2,dthr,hrifac, hdia0,hdia1,civenta,civentb
      real civentc,civentd,civente,civentf,civentg,cireyn,xcivent
      real cipventa,cipventb,cipventc,cipventd,cipreyn,cirventa
      real cirventb
      integer igmrwa,igmrwb,igmswa, igmswb,igmfwa,igmfwb,igmhwa,igmhwb
      real rwventa ,rwventb,swventa,swventb,fwventa,fwventb,fwventc
      real hwventa,hwventb
      real    hwventc, hlventa, hlventb,  hlventc
      real  glventa, glventb, glventc
      real   gmventa, gmventb,  gmventc, ghventa, ghventb, ghventc
      real  dzfacp,  dzfacm,  cmassin,  cwdiar 
      real  rimmas, rhobar
      real   argtim, argqcw, argqxw, argtem
      real   frcswsw, frcswgl, frcswgm, frcswgh, frcswfw, frcswsw1
      real   frcglgl, frcglgm, frcglgh,  frcglfw, frcglgl1
      real   frcgmgl, frcgmgm, frcgmgh,  frcgmfw, frcgmgm1
      real   frcghgl, frcghgm, frcghgh,  frcghfw,  frcghgh1
      real   frcfwgl, frcfwgm, frcfwgh, frcfwfw,  frcfwfw1
      real   frcswrsw, frcswrgl,  frcswrgm,  frcswrgh, frcswrfw
      real   frcswrsw1
      real   frcrswsw, frcrswgl, frcrswgm, frcrswgh, frcrswfw
      real  frcrswsw1
      real  frcglrgl, frcglrgm, frcglrgh,  frcglrfw, frcglrgl1
      real  frcrglgl
      real  frcrglgm,  frcrglgh, frcrglfw, frcrglgl1
      real  frcgmrgl, frcgmrgm, frcgmrgh, frcgmrfw,  frcgmrgm1
      real  frcrgmgl, frcrgmgm,  frcrgmgh, frcrgmfw, frcrgmgm1
      real  total,  qweps,  gf2a, gf4a, dqldt, dqidt, dqdt
      real frcghrgl, frcghrgm, frcghrgh, frcghrfw, frcghrgh1, frcrghgl
      real frcrghgm, frcrghgh,  frcrghfw, frcrghgh1
      real    a1,a2,a3,a4,a5,a6
      real   gamss
      real cdw, cdi, denom1, denom2, delqci1, delqip1
      real cirtotn,  ciptotn, cgmtotn, chltotn,  cirtotp
      real  cgmfac, chlfac,  cirfac
      integer igmhla, igmhlb, igmgla, igmglb, igmgma,  igmgmb
      integer igmgha, igmghb
      integer idqis, item, itim0 
      integer  iqgl, iqgm, iqgh, iqrw, iqsw 
      integer  itertd, ia
      
      integer :: infdo
      
      real tau, ewtmp
      
      integer cntnic_noliq
      real     q_noliqmn, q_noliqmx
      real     scsacimn, scsacimx
      
      real :: dtpinv
      
!   arrays for temporary bin space

      real :: xden,xmlt,cmlt,cmlttot,fventm,fventh,am,ah,felfinv,dmwdt

      real :: qhmlrtmp,qhmlrtmp2, chmlrtmp, chmlrtmpd1inf, chlmlrtmp, zhlmlrtmp, zhlmlrrtmp, qvs0,tmpcmlt

      real :: term1,term2,term3,term4
      real :: qaacw ! combined qsacw-qhacw for WSM6 variation
      real :: cwchtmp

      real, parameter ::  c1r=19.0, c2r=0.6, c3r=1.8, c4r=17.0   ! rain
      real, parameter ::  c1h=5.5, c2h=0.7, c3h=4.5, c4h=8.5   ! Graupel
      real, parameter ::  c1hl=3.7, c2hl=0.3, c3hl=9.0, c4hl=6.5, c5hl=1.0, c6hl=6.5 ! Hail


! inline functions for Newton method
       real :: galpha, dgalpha
       real :: a_in
       logical, parameter :: newton = .false.


      galpha(a_in) = ((4. + a_in)*(5. + a_in)*(6. + a_in))/((1. + a_in)*(2. + a_in)*(3. + a_in))
      dgalpha(a_in) = (876. + 1260.*a_in + 621.*a_in**2 + 126.*a_in**3 + 9.*a_in**4)/            &
     &  (36. + 132.*a_in + 193.*a_in**2 + 144.*a_in**3 + 58.*a_in**4 + 12.*a_in**5 + a_in**6)
!
! ####################################################################
!
!  Start routine
!
! ####################################################################



!

       pb(:) = 0.0
       pinit(:) = 0.0
      itile = nx
      jtile = ny
      ktile = nz
      ixend = nx
      jyend = ny
      kzend = nz
      nxend = nx + 1
      nyend = ny + 1
      nzend = nz
      kzbeg = 1
      nzbeg = 1

      istag = 0
      jstag = 0
      kstag = 1

      lrescalelow(:) = rescale_low_alpha
      lrescalelow(lr) = rescale_low_alphar .and. rescale_low_alpha
      lrescalelow(lh) = rescale_low_alphah .and. rescale_low_alpha
      if ( lf > 1 ) lrescalelow(lf) = rescale_low_alphah .and. rescale_low_alpha
      if ( lhl > 1 ) lrescalelow(lhl) = rescale_low_alphahl .and. rescale_low_alpha


!
!  slope intercepts
!

      if ( ngs .lt. nz ) then
!       write(0,*) 'Error in ICEZVD: Must have ngs .ge. nz!'
!       STOP
      endif

      cntnic_noliq = 0
      q_noliqmn = 0.0
      q_noliqmx = 0.0
      scsacimn = 0.0
      scsacimx = 0.0

      ldovol = .false.

      do il = lc,lhab
        ldovol = ldovol .or. ( lvol(il) .gt. 1 )
      enddo


      ffrzh = 1
!      DO il = lc,lhab
!        write(iunit,*) 'delqnxa(',il,') = ',delqnxa(il)
!      ENDDO
      
!
!  density maximums and minimums
!

!
!  Set terminal velocities...
!    also set drag coefficients
!

      dtpinv = 1.d0/dtp

!

!
!  electricity constants
!
!  mixing ratio epsilon
!
      qeps  = 1.0e-20

!  rebound efficiency (erbnd)
!
!
!
!  constants
!

!      cp608 = 0.608
      aradcw = -0.27544
      bradcw = 0.26249e+06
      cradcw = -1.8896e+10
      dradcw = 4.4626e+14
      bta1 = 0.6
      cnit = 1.0e-02
      dragh = 0.60
      dnz00 = 1.225
!      cs = 4.83607122
!      ds = 0.25
!  new values for  cs and ds
      cs = 12.42
      ds = 0.42
      pii = piinv ! 1./pi
      pid4 = pi/4.0 
!      qscrit = 6.0e-04
      gf1 = 1.0 ! gamma(1.0)
      gf1p5 = 0.8862269255  ! gamma(1.5)
      gf2 = 1.0 ! gamma(2.0)
      gf3 = 2.0 ! gamma(3.0)
      gf3p5 = 3.32335097 ! gamma(3.5)
      gf4 = 6.00 ! gamma(4.0)
      gf5 = 24.0 ! gamma(5.0)
      gf6 = 120.0 ! gamma(6.0)
      gf7 = 720.0 ! gamma(7.0)
      gf4br = 17.837861981813607 ! gamma(4.0+br)
      gf4ds = 10.41688578110938 ! gamma(4.0+ds)
      gf4p5 = 11.63172839656745 ! gamma(4.0+0.5)
      gf3ds = 3.0458730354120997 ! gamma(3.0+ds)
      gf1ds = 0.8863557896089221 ! gamma(1.0+ds)

      gf43rds = 0.8929795116 ! gamma(4./3.)
      gf53rds = 0.9027452930 ! gamma(5./3.)
      gf73rds = 1.190639349 ! gamma(7./3.)
      gf83rds = 1.504575488 ! gamma(8./3.)
      
      gamice73fac =  (Gamma_sp(7./3. + cinu))**3/ (Gamma_sp(1. + cinu)**3 * (1. + cinu)**4)
      gamsnow73fac =  (Gamma_sp(7./3. + snu))**3/ (Gamma_sp(1. + snu)**3 * (1. + snu)**4)
      
!      gcnup1 = Gamma_sp(cnu + 1.)
!      gcnup2 = Gamma_sp(cnu + 2.)
!
!  constants
!
!
!  general constants for microphysics
!
      brz = 100.0
      arz = 0.66
      
      bfnu1 = (4. + alphar)*(5. + alphar)*(6. + alphar)/ &
     &       ((1. + alphar)*(2. + alphar)*(3. + alphar))

       galpharaut = (6.+alpharaut)*(5.+alpharaut)*(4.+alpharaut)/ &
     &             ((3.+alpharaut)*(2.+alpharaut)*(1.+alpharaut))
      
      vfrz = 0.523599*(dfrz)**3 
      vmlt = Min(xvmx(lr), 0.523599*(dmlt)**3 )
      vshd = Min(xvmx(lr), 0.523599*(dshd)**3 )

      if ( snowmeltdia > 0.0 ) then
        snowmeltmass = pi/6.0 * 1000. * snowmeltdia**3  ! maximum rain particle mass from melting snow (if snowmeltdia > 0)
      endif

      tdtol = 1.0e-05
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi
      
      if ( mixedphase ) then
       ibinhmlr = 0
       ibinhlmlr = 0
      endif
!
!
! #ifdef COMMAS
!      print*,'ventr,ventc = ',ventr,ventc

!
!  Set up look up tables for supersaturation w.r.t. liq and ice
!
!VD$L SKIP
!      do l = 1,nqsat
!      temq = 163.15 + (l-1)*fqsat
!      tabqvs(l) = exp(caw*(temq-273.15)/(temq-cbw))
!      tabqis(l) = exp(cai*(temq-273.15)/(temq-cbi))
!      end do

      mltmass0inv = 1.0/( 1000.0* xvmx(lr) ) ! for drops melting from ice with diameter > 1.9cm
      mltmass1inv = 1.0/( 1000.0*(4.0*pi/3.0)*((0.01*0.5*takshedsize1)**3) ) ! for drops melting from ice with diameter > 1.9cm; 0.01 converts cm to m, 0.5 conv. diam to radius
      mltmass2inv = 1.0/( 1000.0*(4.0*pi/3.0)*((0.01*0.5*takshedsize2)**3) ) ! for drops melting from ice with 0.9cm < d < 1.9cm (or 1.6cm to 1.9cm)
      mltmass3inv = 1.0/( 1000.0*(4.0*pi/3.0)*((0.01*0.5*takshedsize3)**3) ) ! for drops melting from ice with 0.9cm < d < 1.6cm
      mltmass1cgs =  1.0*(4.0*pi/3.0)*((0.5*takshedsize1)**3) 
      mltmass2cgs =  1.0*(4.0*pi/3.0)*((0.5*takshedsize2)**3) 
      mltmass3cgs =  1.0*(4.0*pi/3.0)*((0.5*takshedsize3)**3) 
      
!      real, parameter :: mltdiam1 = 9.0e-3, mltdiam2 = 19.0e-3, mltdiam05 = 4.5e-3

      if ( ibinnum == 1 ) then
        numdiam = 1 ! must have numdiam < ndiam because numdiam+1 holds values for the interval of mltdiam(numdiam) to mltdiam(ndiam+1)
        mltdiam(1) = 4.5e-3
      elseif ( ibinnum == 2 ) then
        numdiam = 2 ! must have numdiam < ndiam because numdiam+1 holds values for the interval of mltdiam(numdiam) to mltdiam(ndiam+1)
        mltdiam(1) = mltdiam1/6. ! 1.5e-3
        mltdiam(2) = mltdiam1/2. ! 4.5e-3
      elseif ( ibinnum > 2 ) then
        numdiam = Min(ibinnum, ndiam)
        do k = 1,numdiam
          mltdiam(k) = (k - 0.5)*mltdiam1/float(numdiam)
        enddo
      
      else
        numdiam = 5 ! must have numdiam < ndiam because numdiam+1 holds values for the interval of mltdiam(numdiam) to mltdiam(ndiam+1)
        mltdiam(1) = 0.5e-3
        mltdiam(2) = 1.0e-3
        mltdiam(3) = 2.0e-3
        mltdiam(4) = 4.0e-3
        mltdiam(5) = 6.0e-3
      endif


      if ( numshedregimes == 2 ) then
        mltdiam(ndiam+1) = mltdiam1 !  9.0e-3
        mltdiam(ndiam+2) = mltdiam3 ! 19.0e-3
        mltdiam(ndiam+3) = mltdiam4 !100.0e-3
      elseif ( numshedregimes == 3 ) then
        mltdiam(ndiam+1) = mltdiam1 !  9.0e-3
        mltdiam(ndiam+2) = mltdiam2 ! 16.0e-3
        mltdiam(ndiam+3) = mltdiam3 ! 19.0e-3
        mltdiam(ndiam+4) = mltdiam4 !200.0e-3
      endif

      kzb = 1
      kze = ktile
!      if (kzend .eq. nzend) kze = kzend-kzbeg+1-kstag

!
!  cw constants in mks units
!
!      cwmasn = 4.25e-15  ! radius of 1.0e-6
      mwfac = 6.0**(1./3.)
      if ( ipconc .ge. 2 ) then
!        cwmasn = xvmn(lc)*1000.
!        cwradn = 1.0e-6
!        cwmasx = xvmx(lc)*1000.
      endif
        rwmasn = xvmn(lr)*1000.
        rwmasx = xvmx(lr)*1000.

      if ( biggsnowdiam > 0.0 ) then
        xvbiggsnow = (pi/6.0)*biggsnowdiam**3
      else
        xvbiggsnow = xvmn(lh)
      endif

!
!  ci constants in mks units
!
      cimasn = Min(cimas0, cimas1) ! 12 microns for  0.1871*(xmas(mgs,li)**(0.3429))
      cimasx = 1.0e-8   ! 338 microns
      ccimx = 5000.0e3   ! max of 5000 per liter

!
!  constants for paramerization
!
!
!  set save counter (number of saves):  nsvcnt
!
!      nsvcnt = 0
      iend = 0


!      timetd1 = etime(tarray)
!      timetd1 = tarray(1)

!
!***********************************************************
!  start jy loop
!***********************************************************
!

!      do 9999 jy = 1,ny-jstag
!
!  VERY IMPORTANT:  SET jy = jgs
!
      jy = jgs
     
     
!      t1(:,:,:) = 0
!      t2(:,:,:) = 0
!      t3(:,:,:) = 0
!      t4(:,:,:) = 0
!      t5(:,:,:) = 0
!      t6(:,:,:) = 0
!      t8(:,:,:) = 0
      
      if ( ipconc < 2 ) then ! Make a copy of cloud droplet mixing ratio to use for homogeneous freezing
        do kz = 1,kze
         do ix = 1,itile
           t9(ix,jy,kz) = an(ix,jy,kz,lc)
         enddo
        enddo
      endif
      
!
!..Gather microphysics  
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: ENTER GATHER STAGE'


      
      nxmpb = 1
      nzmpb = 1
      nxz = itile*nz
      numgs = nxz/ngs + 1
!      write(0,*) 'ICEZVD_GS: ENTER GATHER STAGE: nx,nz,nxz,numgs,ngs = ',nx,nz,nxz,numgs,ngs

      do 1000 inumgs = 1,numgs
      ngscnt = 0
      
      do kz = nzmpb,kze
      do ix = nxmpb,itile

      pqs(1) = t00(ix,jy,kz)

      theta(1) = an(ix,jy,kz,lt)
      temg(1) = t0(ix,jy,kz)
      temcg(1) = temg(1) - tfr
      tqvcon = temg(1)-cbw
      ltemq = (temg(1)-163.15)/fqsat + 1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvs(1) = pqs(1)*tabqvs(ltemq)
      if ( iqis0 == 1 .or. temg(1) <= tfr+0.5 ) then
        qis(1) = pqs(1)*tabqis(ltemq)
      else
        ltemq = (tfr - 163.15)/fqsat + 1.5
        qis(1) = pqs(1)*tabqis(ltemq)
      endif

      qss(1) = qvs(1)

      if ( temg(1) .lt. tfr ) then
        qss(1) = qis(1)
      end if
!
      ishail = .false.
      if ( lhl > 1 ) then
        if ( an(ix,jy,kz,lhl)  .gt. qxmin(lhl) ) ishail = .true.
      endif


      
      if ( an(ix,jy,kz,lv)  .gt. qss(1) .or.   &
     &     an(ix,jy,kz,lc)  .gt. qxmin(lc)   .or.    &
     &     an(ix,jy,kz,li)  .gt. qxmin(li)   .or.   &
     &     an(ix,jy,kz,lr)  .gt. qxmin(lr)   .or.   &
     &     an(ix,jy,kz,ls)  .gt. qxmin(ls)   .or.   &
     &     an(ix,jy,kz,lh)  .gt. qxmin(lh)   .or.  ishail ) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 1100
      end if
      enddo !ix
      nxmpb = 1
      enddo !kz
 1100 continue

      if ( ngscnt .eq. 0 ) go to 9998

      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: dbg = 5, ngscnt = ',ngscnt
      
!      write(0,*) 'allocating qc'

      
      xv(:,:) = 0.0
      xmas(:,:) = 0.0
      vtxbar(:,:,:) = 0.0
      xdia(:,:,:) = 0.0
      raindn(:,:) = 900.
      cx(:,:) = 0.0
      if ( lnhf > 1 .or. lnhlf > 1 ) chxf(:,:) = 0.0
      alpha(:,:) = 0.0
      do il = li,lhab
        do mgs = 1,ngscnt
          rimdn(mgs,il)  = rimedens ! xdn0(il)
        enddo
      enddo
!
!  define temporaries for state variables to be used in calculations
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: dbg = def temps'
      do mgs = 1,ngscnt
      kgsm(mgs) = max(kgs(mgs)-1,1)
      kgsp(mgs) = min(kgs(mgs)+1,nz-1)
      kgsm2(mgs) = Max(kgs(mgs)-2,1)
      theta0(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      thetap(mgs) = an(igs(mgs),jy,kgs(mgs),lt) - theta0(mgs)
      theta(mgs) = an(igs(mgs),jy,kgs(mgs),lt)
      qv0(mgs) = an(igs(mgs),jy,kgs(mgs),lv)
      qwvp(mgs) = an(igs(mgs),jy,kgs(mgs),lv)  - qv0(mgs) ! qv0(mgs) is full qv, so qwvp starts as zero!

      pres(mgs) = pn(igs(mgs),jy,kgs(mgs)) + pb(kgs(mgs))
      pipert(mgs) = p2(igs(mgs),jy,kgs(mgs))
      rho0(mgs) = dn(igs(mgs),jy,kgs(mgs))
      rhoinv(mgs) = 1.0/rho0(mgs)
      rhovt(mgs) = Sqrt(rho00/Max(0.05,rho0(mgs))) ! prevent excessive rhovt
      pi0(mgs) = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs))
      temg(mgs) = t0(igs(mgs),jy,kgs(mgs))
      temgkm1(mgs) = t0(igs(mgs),jy,kgsm(mgs))
      temgkm2(mgs) = t0(igs(mgs),jy,kgsm2(mgs))
      pk(mgs)   = p2(igs(mgs),jy,kgs(mgs)) + pinit(kgs(mgs)) ! t77(igs(mgs),jy,kgs(mgs))
      temcg(mgs) = temg(mgs) - tfr
      qss0(mgs) = (380.0)/(pres(mgs))
      pqs(mgs) = (380.0)/(pres(mgs))
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      if ( iqis0 == 1 .or. temg(mgs) <= tfr+0.5 ) then
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
      else
        ltemq = (tfr - 163.15)/fqsat + 1.5
        qis(mgs) = pqs(mgs)*tabqis(ltemq)
      endif
      qss(mgs) = qvs(mgs)
!      es(mgs)  = 6.1078e2*tabqvs(ltemq)
!      eis(mgs) = 6.1078e2*tabqis(ltemq)
      cnostmp(mgs) = cno(ls)
!

      il5(mgs) = 0
      if ( temg(mgs) .lt. tfr ) then
      il5(mgs) = 1
      end if
      enddo !mgs
      
      if ( ipconc < 1 .and. lwsm6 ) then
        do mgs = 1,ngscnt
          tmp = Min( 0.0, temcg(mgs) )
          cnostmp(mgs) = Min( 2.e8, 2.e6*exp(0.12*tmp) )
        enddo
      endif


!
! zero arrays that are used but not otherwise set (tm)
!
      do mgs = 1,ngscnt
         qhshr(mgs) = 0.0 
       end do
!
!  set temporaries for microphysics variables
!
      do il = lv,lhab
      do mgs = 1,ngscnt
        qx(mgs,il) = max(an(igs(mgs),jy,kgs(mgs),il), 0.0) 
      enddo
      end do

      qxw(:,:) = 0.0
      qxwlg(:,:) = 0.0




!
!  set concentrations
!
!      ssmax = 0.0
      
      
      if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*)  'ICEZVD_GS: dbg = 5b'
      
      if ( ipconc .ge. 1 ) then
       do mgs = 1,ngscnt
        cx(mgs,li) = Max(an(igs(mgs),jy,kgs(mgs),lni), 0.0)
          if ( qx(mgs,li) .le. qxmin(li) ) then
            cx(mgs,li) = 0.0
          endif

        if ( lcina .gt. 1 ) then
         cina(mgs) = an(igs(mgs),jy,kgs(mgs),lcina)
        else
         cina(mgs) = cx(mgs,li)
        endif
        if ( lcin > 1 ) then
         ccin(mgs) = an(igs(mgs),jy,kgs(mgs),lcin)
        endif
       end do
      end if
      if ( ipconc .ge. 2 ) then
       do mgs = 1,ngscnt
        cx(mgs,lc) = Max(an(igs(mgs),jy,kgs(mgs),lnc), 0.0)
!        cx(mgs,lc) = Min( ccwmx, cx(mgs,lc) )
        if ( qx(mgs,lc) .le. qxmin(lc) ) then
          cx(mgs,lc) = 0.0
        endif
        if ( lss > 1 ) then
        ssmax(mgs) = an(igs(mgs),jy,kgs(mgs),lss)
        endif
        if ( lccn .gt. 1 ) then
         ccnc(mgs) = an(igs(mgs),jy,kgs(mgs),lccn)
        else
         ccnc(mgs) = 0.0
        endif
        if ( lccna .gt. 1 ) then
         ccna(mgs) = an(igs(mgs),jy,kgs(mgs),lccna)
        else
         ccna(mgs) = cx(mgs,lc)
        endif
       end do
!       ELSE
!       cx(mgs,lc) = Abs(ccn)
      end if
      if ( ipconc .ge. 3 ) then
       do mgs = 1,ngscnt
        cx(mgs,lr) = Max(an(igs(mgs),jy,kgs(mgs),lnr), 0.0)
        if ( qx(mgs,lr) .le. qxmin(lr) ) then
!          cx(mgs,lr) = 0.0
        elseif ( cx(mgs,lr) .eq. 0.0 .and. qx(mgs,lr) .lt. 3.0*qxmin(lr) ) then
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lr)
          qx(mgs,lr) = 0.0
        else
          cx(mgs,lr) = Max( 1.e-9, cx(mgs,lr) )
        endif
       end do
      end if
      if ( ipconc .ge. 4 ) then
       do mgs = 1,ngscnt
        cx(mgs,ls) = Max(an(igs(mgs),jy,kgs(mgs),lns), 0.0)
        if ( qx(mgs,ls) .le. qxmin(ls) ) then
!          cx(mgs,ls) = 0.0
        elseif ( cx(mgs,ls) .eq. 0.0 .and. qx(mgs,ls) .lt. 3.0*qxmin(ls) ) then
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,ls)
          qx(mgs,ls) = 0.0
        else
          cx(mgs,ls) = Max( 1.e-9, cx(mgs,ls) )

         if ( ilimit .ge. ipc(ls) ) then
            tmp = (xdn0(ls)*cx(mgs,ls))/(rho0(mgs)*qx(mgs,ls))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,ls)*(tmp2)
         if ( cnox .gt. 3.0*cno(ls) ) then
           cx(mgs,ls) = 3.0*cno(ls)/tmp2
         endif
         endif
        endif
       end do
      end if
      if ( ipconc .ge. 5 ) then
       do mgs = 1,ngscnt

        cx(mgs,lh) = Max(an(igs(mgs),jy,kgs(mgs),lnh), 0.0)
        if ( qx(mgs,lh) .le. qxmin(lh) ) then
!          cx(mgs,lh) = 0.0
        elseif ( cx(mgs,lh) .eq. 0.0 .and. qx(mgs,lh) .lt. 3.0*qxmin(lh) ) then
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lh) 
          qx(mgs,lh) = 0.0
        else
          cx(mgs,lh) = Max( 1.e-9, cx(mgs,lh) )
         if ( ilimit .ge. ipc(lh) ) then
            tmp = (xdn0(lh)*cx(mgs,lh))/(rho0(mgs)*qx(mgs,lh))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lh)*(tmp2)
         if ( cnox .gt. 3.0*cno(lh) ) then
           cx(mgs,lh) = 3.0*cno(lh)/tmp2
         endif
         endif
        endif


       end do


      end if

      if ( lhl .gt. 1 .and. ipconc .ge. 5 ) then
       do mgs = 1,ngscnt

        cx(mgs,lhl) = Max(an(igs(mgs),jy,kgs(mgs),lnhl), 0.0)
        if ( qx(mgs,lhl) .le. qxmin(lhl) ) then
          cx(mgs,lhl) = 0.0
        elseif ( cx(mgs,lhl) .eq. 0.0 .and. qx(mgs,lhl) .lt. 3.0*qxmin(lhl) ) then
          qx(mgs,lv) = qx(mgs,lv) + qx(mgs,lhl) 
          qx(mgs,lhl) = 0.0
        else
          cx(mgs,lhl) = Max( 1.e-9, cx(mgs,lhl) )
         if ( ilimit .ge. ipc(lhl) ) then
            tmp = (xdn0(lhl)*cx(mgs,lhl))/(rho0(mgs)*qx(mgs,lhl))
            tmp2 = (tmp*(3.14159))**(1./3.)
            cnox = cx(mgs,lhl)*(tmp2)
         if ( cnox .gt. 3.0*cno(lhl) ) then
           cx(mgs,lhl) = 3.0*cno(lhl)/tmp2
         endif
         endif
        endif


       end do
      end if

!
! Set mean particle volume
!
      if ( ldovol ) then
      
      vx(:,:) = 0.0
      
       do il = li,lhab
        
        if ( lvol(il) .ge. 1 ) then
        
          do mgs = 1,ngscnt
            vx(mgs,il) = Max(an(igs(mgs),jy,kgs(mgs),lvol(il)), 0.0)
          enddo

        endif

       enddo

      endif


!
! Set liquid water fraction
!
      fhw(:) = 0.0
      fsw(:) = 0.0
      fhlw(:) = 0.0



!
!  6th moments
!

      if ( ipconc .ge. 6 ) then
       zx(:,:) = 0.0
       do il = lr,lhab
        if ( lz(il) .gt. 1 ) then
         do mgs = 1,ngscnt
          zx(mgs,il) = Max( an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0 )
         enddo
        endif
       enddo

      endif

      if ( ipconc .ge. 6 ) then

       if ( lz(lr) .lt. 1 ) then
         g1x(:,lr) = (6.0 + alphar)*(5.0 + alphar)*(4.0 + alphar)/ &
     &            ((3.0 + alphar)*(2.0 + alphar)*(1.0 + alphar))

         
         do mgs = 1,ngscnt
           if ( cx(mgs,lr) .gt. 0.0 .and. qx(mgs,lr) .gt. qxmin(lr)  ) then
            
            vr = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            if ( lzr < 1 ) then
             if ( imurain == 3 ) then
               zx(mgs,lr) = 3.6476*(rnu+2.0)*cx(mgs,lr)*vr**2/(rnu+1.0)
             else ! imurain == 1
               zx(mgs,lr) = 3.6476*g1x(mgs,lr)*cx(mgs,lr)*vr**2
             endif
            endif
             
           endif
         enddo
       endif
      
      endif


        scx(:,:) = 0.0
!
!  set shape parameters
!
       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = set alpha'
      if ( imurain == 1 ) then
        alpha(:,lr) = alphar
      elseif ( imurain == 3 ) then
        alpha(:,lr) = xnu(lr)
      endif
      
      alpha(:,li) = xnu(li)
      alpha(:,lc) = xnu(lc)
        if ( idiagnosecnu > 0 ) then
          do mgs = 1,ngscnt
            if ( cx(mgs,lc) > cxmin ) then
              if ( idiagnosecnu == 1 ) then ! diagnose cloud drop DSD shape parameter (cnu) based on Chandrakar et al. 2016 (PNAS) data
                alpha(mgs,lc) = cnudiag( cx(mgs,lc) )
              elseif ( idiagnosecnu > 1 ) then ! Geoffrey et al. (2010, ACP, eq. 12 for alpha = nu - 1)
                alpha(mgs,lc) = 1.58*rho0(mgs)*qx(mgs,lc)*1.e6 + 0.72 - 1.0
              endif
            endif
          enddo
        endif

      if ( imusnow == 1 ) then
        alpha(:,ls) = alphas
      elseif ( imusnow == 3 ) then
        alpha(:,ls) = xnu(ls)
      endif

       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = set dab'
      
      do il = lr,lhab
      do mgs = 1,ngscnt
        if ( il .ge. lg ) alpha(mgs,il) = dnu(il)


        do ic = lc,lhab
        dab0lh(mgs,il,ic) =  dab0(il,ic) ! dab0(ic,il)
        dab1lh(mgs,il,ic) =  dab1(il,ic) ! dab1(ic,il)
        enddo
      end do
      enddo

      
!      DO mgs = 1,ngscnt
        do il = lr,lhab
          da0lx(:,il) = da0(il)
        enddo
        da0lh(:) = da0(lh)
        da0lr(:) = da0(lr)
        da1lr(:) = da1(lr)
        da0lc(:) = da0(lc)
        da1lc(:) = da1(lc)
        if ( idiagnosecnu > 0 ) then
          do mgs = 1,ngscnt
            if ( cx(mgs,lc) > cxmin ) then
              da0lc(mgs) =  delbk(bb(lc), alpha(mgs,lc), xmu(lc), 0)
              da1lc(mgs) =  delbk(bb(lc), alpha(mgs,lc), xmu(lc), 1)
            endif
          enddo
        endif


       if ( ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank,  'ICEZVD_GS: dbg = set rz'

        if ( lzh < 1 .or. lzhl < 1 ) then
          rzxhlh(:) = rzhl/rz
        elseif ( lzh > 1 .and. lzhl > 1 ) then
          rzxhlh(:) = 1.
        endif
        if ( lzr > 1 ) then
          rzxh(:) = 1.
          rzxhl(:) = 1.
        else
          rzxh(:) = rz
          rzxhl(:) = rzhl
        endif
        
        if ( imurain == 1 .and. imusnow == 3 .and. lzr < 1 ) then
          rzxs(:) = rzs
        elseif ( imurain == imusnow .or. lzr > 1 ) then
          rzxs(:) = 1.
        endif
 !     ENDDO
      
      if ( lhl .gt. 1 ) then
      do mgs = 1,ngscnt
        da0lhl(mgs) = da0(lhl)
      enddo
      endif
      
      ventrx(:) = ventr
      ventrxn(:) = ventrn
      gf1palp(:) = gamma_sp(1.0 + alphar)

!
!  set factors
!
      do mgs = 1,ngscnt
!
      ssi(mgs) = qx(mgs,lv)/qis(mgs)
      ssw(mgs) = qx(mgs,lv)/qvs(mgs)
!
      tsqr(mgs) = temg(mgs)**2
!
      temgx(mgs) = min(temg(mgs),313.15)
      temgx(mgs) = max(temgx(mgs),233.15)
      felv(mgs) = 2500837.367 * (273.15/temgx(mgs))**((0.167)+(3.67e-4)*temgx(mgs))
!
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),223.15)
      temcgx(mgs) = temcgx(mgs)-273.15

! felf = latent heat of fusion, fels = LH of sublimation, felv = LH of vaporization
      felf(mgs) = 333690.6098 + (2030.61425)*temcgx(mgs) - (10.46708312)*temcgx(mgs)**2
!
      fels(mgs) = felv(mgs) + felf(mgs)
!
      felvs(mgs) = felv(mgs)*felv(mgs)
      felss(mgs) = fels(mgs)*fels(mgs)
      
        if ( eqtset <= 1 ) then
          felvcp(mgs) = felv(mgs)*cpi
          felscp(mgs) = fels(mgs)*cpi
          felfcp(mgs) = felf(mgs)*cpi
        else
          
          ! equations from appendix in Bryan and Morrison (2012, MWR)
          ! note that rw is Rv in the paper, and rd is R.
          
          tmp = qx(mgs,li)+qx(mgs,ls)+qx(mgs,lh)
          if ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          if ( lf > 1 ) tmp = tmp + qx(mgs,lf)
          cvm = cv+cvv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)

          if ( eqtset == 2 ) then ! compact form from treating dT/dt = theta*d(pi)/dt + pi*d(theta)dt and then applied to theta assuming constant pi
          felvcp(mgs) = (felv(mgs)-rw*temg(mgs))/cvm
          felscp(mgs) = (fels(mgs)-rw*temg(mgs))/cvm
          felfcp(mgs) = felf(mgs)/cvm
          
          else
           ! equivalent version that applies separate updates of latent heating to theta and pi, when both are returned.

          cpm = cp+cpv*qx(mgs,lv)+cpl*(qx(mgs,lc)+qx(mgs,lr))   &
                                  +cpigb*(tmp)
          rmm=rd+rw*qx(mgs,lv)
          
          felvcp(mgs) = (felv(mgs)*cv/(cp) - rw*temg(mgs)*(1.0-rovcp*cpm/rmm))/cvm
          felscp(mgs) = (fels(mgs)*cv/(cp) - rw*temg(mgs)*(1.0-rovcp*cpm/rmm))/cvm
          felfcp(mgs) = felf(mgs)*cv/(cp*cvm)

          felvpi(mgs) = pi0(mgs)*rovcp*(felv(mgs)/(temg(mgs)) - rw*cpm/rmm)/cvm
          felspi(mgs) = pi0(mgs)*rovcp*(fels(mgs)/(temg(mgs)) - rw*cpm/rmm)/cvm 
          felfpi(mgs) = pi0(mgs)*rovcp*(felf(mgs)/(cvm*temg(mgs)))
          
          endif

        endif
!
      fgamw(mgs) = felvcp(mgs)/pi0(mgs)
      fgams(mgs) = felscp(mgs)/pi0(mgs)
!
      fcqv1(mgs) = 4098.0258*pi0(mgs)*fgamw(mgs)
      fcqv2(mgs) = 5807.6953*pi0(mgs)*fgams(mgs)
      fcc3(mgs) = felfcp(mgs)/pi0(mgs)
!
!  fwvdf = water vapor diffusivity
      fwvdf(mgs) = (2.11e-05)*((temg(mgs)/tfr)**1.94)*(101325.0/(pres(mgs)))
!
! fadvisc = 'd' for dynamic viscosity
! fakvisc = 'k' for kinematic viscosity
      fadvisc(mgs) = advisc0*(416.16/(temg(mgs)+120.0))*(temg(mgs)/296.0)**(1.5) ! dynamic visc.
!
      fakvisc(mgs) = fadvisc(mgs)*rhoinv(mgs) ! divide by rho_air to get kinematic visc. (note the 'k' vs. 'd')
!
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),233.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fci(mgs) = (2.118636 + 0.007371*(temcgx(mgs)))*(1.0e+03)
!
      if ( temg(mgs) .lt. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),273.15)
      temcgx(mgs) = max(temcgx(mgs),233.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4203.1548  + (1.30572e-2)*((temcgx(mgs)-35.)**2)   &
     &                 + (1.60056e-5)*((temcgx(mgs)-35.)**4)
      end if
      if ( temg(mgs) .ge. 273.15 ) then
      temcgx(mgs) = min(temg(mgs),308.15)
      temcgx(mgs) = max(temcgx(mgs),273.15)
      temcgx(mgs) = temcgx(mgs)-273.15
      fcw(mgs) = 4243.1688  + (3.47104e-1)*(temcgx(mgs)**2)
      end if
!
      ftka(mgs) = tka0*fadvisc(mgs)/advisc1  ! thermal conductivity: proportional to dynamic viscosity
      fthdf(mgs) = ftka(mgs)*cpi*rhoinv(mgs)
!
      fschm(mgs) = (fakvisc(mgs)/fwvdf(mgs))  ! Schmidt number
      fpndl(mgs) = (fakvisc(mgs)/fthdf(mgs))  ! Prandl number (only used for bin melting)
!
      fai(mgs) = (fels(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbi(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qis(mgs)))
      fav(mgs) = (felv(mgs)**2)/(ftka(mgs)*rw*temg(mgs)**2)
      fbv(mgs) = (1.0/(rho0(mgs)*fwvdf(mgs)*qvs(mgs)))

      kp1 = Min(nz, kgs(mgs)+1 )
      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kp1)   &
     &                  +w(igs(mgs),jgs,kgs(mgs)))

!
      end do
!
!
!   ice habit fractions
!
!
!
!  Set density
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: Set density'
!

      do mgs = 1,ngscnt
        xdn(mgs,li) = xdn0(li)
        xdn(mgs,lc) = xdn0(lc)
        xdn(mgs,lr) = xdn0(lr)
        xdn(mgs,ls) = xdn0(ls)
        xdn(mgs,lh) = xdn0(lh)
        if ( lvol(ls) .gt. 1 ) then
         if ( vx(mgs,ls) .gt. 0.0 .and. qx(mgs,ls) .gt. qxmin(ls) ) then
           xdn(mgs,ls) = Min( xdnmx(ls), Max( xdnmn(ls), rho0(mgs)*qx(mgs,ls)/vx(mgs,ls) ) )
         endif
        endif

        if ( lvol(lh) .gt. 1 ) then
         if ( vx(mgs,lh) .gt. 0.0 .and. qx(mgs,lh) .gt. qxmin(lh) ) then
           if ( mixedphase ) then
           else
             dnmx = xdnmx(lh)
           endif
           xdn(mgs,lh) = Min( dnmx, Max( xdnmn(lh), rho0(mgs)*qx(mgs,lh)/vx(mgs,lh) ) )
           vx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/xdn(mgs,lh)
         
         elseif ( vx(mgs,lh) == 0.0 .and. qx(mgs,lh) .gt. qxmin(lh) ) then ! if volume is zero, need to initialize the default value

           vx(mgs,lh) = rho0(mgs)*qx(mgs,lh)/xdn(mgs,lh)
         
         endif
        endif


        if ( lhl .gt. 1 ) then

          xdn(mgs,lhl) = xdn0(lhl)
          xdntmp(mgs,lhl) = xdn0(lhl)

          if ( lvol(lhl) .gt. 1 ) then
           if ( vx(mgs,lhl) .gt. 0.0 .and. qx(mgs,lhl) .gt. qxmin(lhl) ) then

           if ( mixedphase .and. lhlw > 1 ) then
           else
             dnmx = xdnmx(lhl)
           endif

             xdn(mgs,lhl) = Min( dnmx, Max( xdnmn(lhl), rho0(mgs)*qx(mgs,lhl)/vx(mgs,lhl) ) )
             vx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/xdn(mgs,lhl)
             xdntmp(mgs,lhl) = xdn(mgs,lhl)
         
           elseif ( vx(mgs,lhl) == 0.0 .and. qx(mgs,lhl) .gt. qxmin(lhl) ) then ! if volume is zero, need to initialize the default value

             vx(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/xdn(mgs,lhl)
         
           endif
          endif

        endif


      end do

      if ( ipconc == 5 .and. imydiagalpha == 2 ) then

        cwchtmp = ((3. + dnu(lh))*(2. + dnu(lh))*(1.0 + dnu(lh)))**(-1./3.)
        
        do mgs = 1,ngscnt
          !IF ( igs(mgs) == 19 ) write(0,*) 'k,qr,qh,cr,ch = ',kgs(mgs),qx(mgs,lr),cx(mgs,lr),qx(mgs,lh),cx(mgs,lh)
          if ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) > cxmin ) then
             xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))            ! 
             xdia(mgs,lr,3) = (xv(mgs,lr)*6.0*cwc1)**(1./3.) 
           !  alpha(mgs,lr) = Min(alphamax, c1r*tanh(c2r*(xdia(mgs,lr,3)*1000. - c3r)) + c4r)
           ! IF ( igs(mgs) == 19 ) write(0,*) 'imy: i,k,alpr,xdia = ',igs(mgs),kgs(mgs),alpha(mgs,lr),xdia(mgs,lr,3)*1000.

            ! M&M-C 2010:
             tmp = 4. + alphar
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = 1. + alphar
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = (x/y)**(1./3.)*xdia(mgs,lr,3)*cwchtmp

             alpha(mgs,lr) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
          endif
          if ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) > cxmin ) then
!      MY 2005:
             xv(mgs,lh) = rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh))            ! 
             xdia(mgs,lh,3) = (xv(mgs,lh)*6.*piinv)**(1./3.) ! mwfac*xdia(mgs,lh,1) ! (xv(mgs,lh)*cwc0*6.0)**(1./3.)
!             alpha(mgs,lh) = Min(alphamax, c1h*tanh(c2h*(xdia(mgs,lh,3)*1000. - c3h)) + c4h)

            ! M&M-C 2010:
             tmp = 4. + dnu(lh)
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = 1. + dnu(lh)
             i = Int(dgami*(tmp))
             del = tmp - dgam*i
             y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

             tmp = (x/y)**(1./3.)*xdia(mgs,lh,3)*cwchtmp

             alpha(mgs,lh) = Min(15., 11.8*(1000.*tmp - 0.7)**2 + 2.)
            ! alphan(mgs,lh) = alpha(mgs,lh)
            
           ! IF ( igs(mgs) == 19 ) write(0,*) 'imy: i,k,alph,xdia = ',igs(mgs),kgs(mgs),alpha(mgs,lh),xdia(mgs,lh,3)*1000.
            il = lh
            do ic = lc,lh-1 ! lhab
               i = Nint( alpha(mgs,il)*dqiacralphainv )
               if ( ic == lc .or. ic == li .or. ic == ls .or. (ic == lr .and. imurain == 3) ) then
                 alp = (3.*alpha(mgs,ic) + 2.)
                 j = Nint( (3.*alpha(mgs,ic) + 2.)*dqiacralphainv )
               else ! if ( ic == lr .and. imurain == 1 ) ! rain
                 alp = alpha(mgs,ic)
                 j = Nint( alpha(mgs,ic)*dqiacralphainv )
               endif
             
               dab0lh(mgs,ic,il) = dab0lu(j,i,ic,il)
               dab1lh(mgs,ic,il) = dab1lu(j,i,ic,il)
               dab0lh(mgs,il,ic) = dab0lu(i,j,il,ic)
               dab1lh(mgs,il,ic) = dab1lu(i,j,il,ic)
             enddo
          endif
!        alpha(:,lr) = 0. ! 10.
!        alpha(:,lh) = 0. ! 10.
          if ( lhl > 0 ) then
          if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) > cxmin ) then
             xv(mgs,lhl) = rho0(mgs)*qx(mgs,lhl)/(xdn(mgs,lhl)*cx(mgs,lhl))            ! 
             xdia(mgs,lhl,3) = (xv(mgs,lhl)*6.*piinv)**(1./3.)
             if ( xdia(mgs,lhl,3) < 0.008 ) then
               alpha(mgs,lhl) = Min(alphamax, c1hl*tanh(c2hl*(xdia(mgs,lhl,3)*1000. - c3hl)) + c4hl)
             else
               alpha(mgs,lhl) = Min(alphamax, c5hl*xdia(mgs,lhl,3)*1000. + c6hl)
             endif

            il = lhl
            do ic = lc,lh-1 ! lhab
               i = Nint( alpha(mgs,il)*dqiacralphainv )
               if ( ic == lc .or. ic == li .or. ic == ls .or. (ic == lr .and. imurain == 3) ) then
                 alp = (3.*alpha(mgs,ic) + 2.)
                 j = Nint( (3.*alpha(mgs,ic) + 2.)*dqiacralphainv )
               else ! if ( ic == lr .and. imurain == 1 ) ! rain
                 alp = alpha(mgs,ic)
                 j = Nint( alpha(mgs,ic)*dqiacralphainv )
               endif
             
               dab0lh(mgs,ic,il) = dab0lu(j,i,ic,il)
               dab1lh(mgs,ic,il) = dab1lu(j,i,ic,il)
               dab0lh(mgs,il,ic) = dab0lu(i,j,il,ic)
               dab1lh(mgs,il,ic) = dab1lu(i,j,il,ic)
             enddo

          endif
          endif



        enddo
      endif
      

       if ( imurain == 3 ) then
         if ( lzr > 1 ) then
           alphashr = 0.0
           alphamlr = -2.0/3.0
           alphasmlr = -2.0/3.0
         else
           alphashr = xnu(lr)
           alphamlr = xnu(lr)
           alphasmlr = xnu(lr)
         endif
!         massfacshr = ( (2. + 3.*(1. +alphashr) )/( 3.*(1. + alphashr) ) )**(1./3.) ! this is the diameter factor
!         massfacmlr = ( (2. + 3.*(1. +alphamlr) )/( 3.*(1. + alphamlr) ) )**(1./3.)
         massfacshr = ( (2. + 3.*(1. +alphashr) )**3/( 3.*(1. + alphashr) ) )  ! this is the mass or volume factor
         massfacmlr = ( (2. + 3.*(1. +alphamlr) )**3/( 3.*(1. + alphamlr) ) )
       elseif ( imurain == 1 ) then
         if ( lzr > 1 ) then
           alphashr = 4.0
           alphamlr = 4.0
           alphasmlr = alphasmlr0
         else
           alphashr = alphar
           alphamlr = alphar
           alphasmlr = alphar
         endif
!         massfacshr = (3.0 + alphashr)*((3.+alphashr)*(2.+alphashr)*(1. + alphashr) )**(-1./3.) ! this is the diameter factor
!         massfacmlr = (3.0 + alphamlr)*((3.+alphamlr)*(2.+alphamlr)*(1. + alphamlr) )**(-1./3.)
         massfacshr = (3.0 + alphashr)**3/((3.+alphashr)*(2.+alphashr)*(1. + alphashr) ) ! this is the mass or volume factor
         massfacmlr = (3.0 + alphamlr)**3/((3.+alphamlr)*(2.+alphamlr)*(1. + alphamlr) )
       endif
       
!  Find shape parameter rain

      g1shr = 1.0
      g1mlr = 1.0
      g1smlr = 1.0
 
!      CALL cld_cpu('Z-MOMENT-1')  
      
      if ( ipconc >= 6 ) then
      
      ! set base g1x in case rain is not 3-moment
       if ( ipconc >= 6 .and. imurain == 3 ) then
         il = lr
         do mgs = 1,ngscnt
!           g1x(mgs,il) = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
           g1x(mgs,il) = (alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0))
         enddo
       endif

      if (lzr > 1 ) then
       if ( imurain == 3 ) then
         g1shr = (alphashr+2.0)/((alphashr+1.0))
         g1mlr = (alphamlr+2.0)/((alphamlr+1.0))
         g1smlr = (alphasmlr+2.0)/((alphasmlr+1.0))
       elseif ( imurain == 1 ) then
!         g1shr = 36.*(6.0 + alphashr)*(5.0 + alphashr)*(4.0 + alphashr)/ &
!     &            (pi**2*(3.0 + alphashr)*(2.0 + alphashr)*(1.0 + alphashr))
         g1shr = (6.0 + alphashr)*(5.0 + alphashr)*(4.0 + alphashr)/ &
     &            ((3.0 + alphashr)*(2.0 + alphashr)*(1.0 + alphashr))
!         g1mlr = 36.*(6.0 + alphamlr)*(5.0 + alphamlr)*(4.0 + alphamlr)/ &
!     &            (pi**2*(3.0 + alphamlr)*(2.0 + alphamlr)*(1.0 + alphamlr))
         g1mlr = (6.0 + alphamlr)*(5.0 + alphamlr)*(4.0 + alphamlr)/ &
     &            ((3.0 + alphamlr)*(2.0 + alphamlr)*(1.0 + alphamlr))
         g1smlr = (6.0 + alphasmlr)*(5.0 + alphasmlr)*(4.0 + alphasmlr)/ &
     &            ((3.0 + alphasmlr)*(2.0 + alphasmlr)*(1.0 + alphasmlr))
       endif
      endif

      if ( lzr > 1 .and. imurain == 3 ) then ! { RAIN SHAPE PARAM
      
      
!      CALL cld_cpu('Z-MOMENT-1r')  
          il = lr
          do mgs = 1,ngscnt
          

         if ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) then ! .or. qx(mgs,il) <= qxmin(il)  then
         if ( zx(mgs,il) <= zxmin ) then !  .and. qx(mgs,il) > 0.05e-3  then
!!            write(91,*) 'zx=0; qx,cx = ',1000.*qx(mgs,il),cx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         elseif ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         elseif ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) then !  .and. qx(mgs,il) > 0.05e-3   then
         
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         endif
         endif

         if ( .false. .and. zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif
         
         if ( qx(mgs,lr) .gt. qxmin(lr) ) then

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        if ( xv(mgs,lr) .gt. xvmx(lr) ) then
!          xv(mgs,lr) = xvmx(lr)
!          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        elseif ( xv(mgs,lr) .lt. xvmn(lr) ) then
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*xdn(mgs,lr)**2)
!            an(igs(mgs),jgs,kgs(mgs),ln(il)) = zx(mgs,il)
           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
          endif
        
          if ( zx(mgs,lr) > 0.0 ) then
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           if ( z .gt. 0.0 ) then
!           alpha(mgs,lr) = 3.
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           do i = 1,20
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo

! check for artificial breakup (rain larger than allowed max size)
        if (  (xv(mgs,il) .gt. xvmx(il) .or. (ioldlimiter >= 2 .and. xv(mgs,il) .gt. xvmx(il)/8.) )) then
          tmp = cx(mgs,il)
          if ( ioldlimiter >= 2 ) then ! MY-style active breakup
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.222e3*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.222e3*x1**2), x3)
          else ! simple cutoff 
            xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          endif
            !xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            !cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))

          if ( tmp < cx(mgs,il) ) then ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           do i = 1,20
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo

            
          endif
        endif

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
           if ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) then

            if ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( rescale_low_alphar .and. alp <= rnumin ) then
             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = zx(mgs,il)
            endif
           endif
           
         ! set g1x to use as G factor later. If alpha is in the range ( rnumin < alpha < rnumax ), then 
         ! this will be the same as computing G from alpha.  If alpha = rnumax, however, it probably means that
         ! the moments are not matched correctly, so we compute G from the moments instead so that the dZ/dt rates
         ! stay consistent with dN/dt and dq/dt.
           if ( alp >= rnumax - 0.01 ) then
!             g1x(mgs,il) = 6**2*zx(mgs,il)/(cx(mgs,il)*(pi*xv(mgs,lr))**2)
!             g1x(mgs,il) = xdn(mgs,il)*zx(mgs,il)*cx(mgs,il)/((rho0(mgs)*qx(mgs,lr))**2)
             g1x(mgs,il) = (pi*xdn(mgs,il))**2*zx(mgs,il)*cx(mgs,il)/((6.*rho0(mgs)*qx(mgs,il))**2)
           else
             g1x(mgs,il) = g1
           endif
           
           tmp = alpha(mgs,lr) + 4./3.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           tmp = alpha(mgs,lr) + 1.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
           
           gf1palp(mgs) = y

!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma_sp(alpha(mgs,lr) + 1.)
           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))

           if ( imurain == 3 .and. izwisventr == 2 ) then

           tmp = alpha(mgs,lr) + 1.5 + br/6.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 1.5 + br/6.)/Gamma_sp(alpha(mgs,lr) + 1.)
           ventrxn(mgs) = x/(y*(alpha(mgs,lr) + 1.)**((1.+br)/6. + 1./3.))
           
! This whole section is imurain == 3, so this branch never runs
!           ELSEIF ( imurain == 1 .and.  iferwisventr == 2 ) THEN
!
!           tmp = alpha(mgs,lr) + 2.5 + br/2.
!           i = Int(dgami*(tmp))
!           del = tmp - dgam*i
!           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!
!!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 1.5 + br/6.)/Gamma_sp(alpha(mgs,lr) + 1.)
!           ventrxn(mgs) = x/y
           
           
           endif
           
           endif
          endif
          
          endif
          
          enddo
!        CALL cld_cpu('Z-MOMENT-1r')  
        endif ! }
        
      endif ! ipconc >= 6

!  Find shape parameters for graupel and hail
      if ( ipconc .ge. 6 ) then
            
        do il = lr,lhab
          
        ! set base values of g1x
          if ( (.not. ( il == lr .and. imurain == 3 )) .and. ( il == lr .or. il == lh .or. il == lhl .or. il == lf ) ) then
          do mgs = 1,ngscnt
            g1x(mgs,il) = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
          enddo
          endif
        
        if ( lz(il) .gt. 1   .and. ( .not. ( il == lr .and. imurain == 3 )) ) then
        
        do mgs = 1,ngscnt


         if ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) then ! .or. qx(mgs,il) <= qxmin(il) ) then
         if ( zx(mgs,il) <= zxmin ) then !  .and. qx(mgs,il) > 0.05e-3 ) then
!!            write(91,*) 'zx=0; qx,cx = ',1000.*qx(mgs,il),cx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           zx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         elseif ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         elseif ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) then !  .and. qx(mgs,il) > 0.05e-3  ) then
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         endif
         endif

         if ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif
        
        if ( qx(mgs,il) .gt. qxmin(il) ) then

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        if ( xv(mgs,il) .lt. xvmn(il) ) then
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
        endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
!            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > cxmin ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
!            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            g1 = (6.0 + alphamax)*(5.0 + alphamax)*(4.0 + alphamax)/ &
     &            ((3.0 + alphamax)*(2.0 + alphamax)*(1.0 + alphamax))
            zx(mgs,il) = Max(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
!               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           else
          
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          if ( zx(mgs,il) .gt. 0. ) then
           
!            rdi = z*(pi/6.*1000.)**2*chw/((rho0(mgs)*qr)**2)
            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

!           alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           alp = Max( alphamin, Min( alphamax, alp ) )
           
         if ( newton ) then
           do i = 1,10
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = alp + ( galpha(alp) - rdi )/dgalpha(alp)
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo
           
         else
           do i = 1,10
!            IF ( 100.*Abs(alp - alpha(mgs,il))/(Abs(alpha(mgs,il))+1.e-5) .lt. 1. ) EXIT
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo
          endif


! check for artificial breakup (graupel/hail larger than allowed max size)
        if ( imaxdiaopt == 1 ) then
          xvbarmax = xvmx(il) 
        elseif ( imaxdiaopt == 2 ) then ! test against maximum mass diameter
          xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        elseif ( imaxdiaopt == 3 ) then ! test against mass-weighted diameter
          xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
        else
          xvbarmax = xvmx(il) 
        endif

        if (  xv(mgs,il) .gt. xvbarmax .or. (il == lr .and. ioldlimiter >= 2 .and. xv(mgs,il) .gt. xvmx(il)/8.)) then
          tmp = cx(mgs,il)
          if( ioldlimiter >= 2 .and. il == lr) then ! MY-style drop limiter for rain
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.222e3*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.222e3*x1**2), x3)
          else
            xv(mgs,il) = Min( xvbarmax, Max( xvmn(il),xv(mgs,il) ) )
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          endif
          if ( tmp < cx(mgs,il) ) then ! artificial breakup has happened, so need to adjust reflectivity and find new shape parameter
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           do i = 1,10
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo

            
          endif
        endif

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
 
           if ( ( lrescalelow(il) .or. rescale_high_alpha ) .and.  &
     &          ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) then



            if ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( lrescalelow(il) .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) .and. &
                     .not. ( il == lr .and. .not. rescale_low_alphar ) ) then ! alpha = alphamin, so reset Z to prevent growth in C
             wtest = .false.
             if ( irescalerainopt == 0 ) then
               wtest = .false.
             elseif ( irescalerainopt == 1 ) then
               wtest = qx(mgs,lc) > qxmin(lc) 
             elseif ( irescalerainopt == 2 ) then
               wtest = qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             elseif ( irescalerainopt == 3 ) then
               wtest = temcg(mgs) > rescale_tempthresh .and. qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             endif
             
             if ( il == lr .and. ( wtest ) ) then
!             IF ( temcg(mgs) > 0.0 .and. il == lr .and. qx(mgs,lc) > qxmin(lc) ) THEN
             ! certain situations where rain number is adjusted instead of Z. Helps avoid rain being 'zapped' by autoconverted 
             ! drops (i.e., favor preserving Z when alpha tries to go negative)
             chw = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2 ! g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z1
             cx(mgs,il) = chw
             an(igs(mgs),jy,kgs(mgs),ln(il)) = chw
             else
             
             ! Usual resetting of reflectivity moment to force consisntency between Q, N, Z, and alpha when alpha = alphamin
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             endif
            endif
           endif
          
          
         ! set g1x to use as G factor later. If alpha is in the range ( rnumin < alpha < rnumax ), then 
         ! this will be the same as computing G from alpha.  If alpha = rnumax, however, it probably means that
         ! the moments are not matched correctly, so we compute G from the moments instead so that the dZ/dt rates
         ! stay consistent with dN/dt and dq/dt.
!          g1x(mgs,il) = zx(mgs,il)*chw*(pi*xdn(mgs,il))**2/(6.*qr*dn(igs(mgs),jy,kgs(mgs)))**2
!          g1x(mgs,il) = g1 ! zx(mgs,il)*cx(mgs,il)/(qr)**2
           if ( alp >= alphamax - 0.5 ) then
!             g1x(mgs,il) = 6**2*zx(mgs,il)/(cx(mgs,il)*(pi*xv(mgs,lr))**2)
!             g1x(mgs,il) = (xdn(mgs,il))**2*zx(mgs,il)*cx(mgs,il)/((rho0(mgs)*qx(mgs,il))**2)
             g1x(mgs,il) = (pi*xdn(mgs,il))**2*zx(mgs,il)*cx(mgs,il)/((6.*rho0(mgs)*qx(mgs,il))**2)
           else
             g1x(mgs,il) = g1
           endif
          
           endif
          
!          IF ( ny .eq. 2 ) THEN
!          IF ( qr .gt. 1.e-3 ) THEN
!           write(0,*) 'alphah at nstep,i,k = ',dtp*(nstep-1),igs(mgs),kgs(mgs),alpha(mgs,il),qr*1000.
!          ENDIF
!          ENDIF
          
           
           endif ! .true.

          if ( il == lr ) then
           
!           tmp = alpha(mgs,lr) + 4./3.
!           i = Int(dgami*(tmp))
!           del = tmp - dgam*i
!           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!
!           tmp = alpha(mgs,lr) + 1.
!           i = Int(dgami*(tmp))
!           del = tmp - dgam*i
!           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
!
!!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 4./3.)/(alpha(mgs,lr) + 1.)**(1./3.)/Gamma_sp(alpha(mgs,lr) + 1.)
!           ventrx(mgs) = x/(y*(alpha(mgs,lr) + 1.)**(1./3.))


           tmp = alpha(mgs,lr) + 1.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

           gf1palp(mgs) = y

           if (   iferwisventr == 2 ) then
           tmp = alpha(mgs,lr) + 2.5 + br/2.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!           ventrx(mgs) = Gamma_sp(alpha(mgs,lr) + 1.5 + br/6.)/Gamma_sp(alpha(mgs,lr) + 1.)

           ventrxn(mgs) = x/y
           
           endif
           
          endif ! il==lr
 
          
          else ! below mass threshold
!             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
!             z1 = g1*rho0(mgs)**2*(qr)*qr/chw
!             z  = 1.e18*z1*(6./(pi*1000.))**2
!             z  = z1*(6./(pi*1000.))**2
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
          endif ! ( qx(mgs,il) .gt. qxmin(il) )
        
        
        
!        ENDIF
        enddo ! mgs

!         CALL cld_cpu('Z-DELABK')  
        
!        IF ( il == lr ) THEN
!          xnutmp = (alpha(mgs,il) - 2.)/3.
!           da0lr(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
!        ENDIF
        
        if ( .not. ( il == lr .and. imurain == 3 ) ) then
!          CALL cld_cpu('Z-DELABK')  
        do mgs = 1,ngscnt
          if ( qx(mgs,il) > qxmin(il) ) then
          xnutmp = (alpha(mgs,il) - 2.)/3.
          
!          IF ( .true. ) THEN
          do ic = lc,lh-1 ! lhab
           if ( il .ne. ic .and.  qx(mgs,ic) .gt. qxmin(ic)) then
             xnuc = xnu(ic)
             if ( ic == lc .and. idiagnosecnu > 0 ) xnuc = alpha(mgs,lc) ! alpha for droplets is actually nu
             if ( il /= lr .and. ic == lr .and. lzr > 1 ) then
               if ( imurain == 3 ) then
                 xnuc = alpha(mgs,lr) ! alpha is nu already
               else
                 xnuc = ( alpha(mgs,lr) - 2. )/3. ! convert alpha to nu
               endif
             endif
                                 ! delabk(ba,bb,nua,nub,mua,mub,k), where a (il)  is collector and b (ic) is collected
             if ( .false. ) then
             dab0lh(mgs,ic,il) =  delabk(bb(ic), bb(il), xnuc, xnutmp, xmu(ic), xmu(il), 0) !dab0(il,ic)
             dab1lh(mgs,ic,il) =  delabk(bb(ic), bb(il), xnuc, xnutmp, xmu(ic), xmu(il), 1) !dab1(il,ic)
             dab0lh(mgs,il,ic) =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 0) !dab0(il,ic)
             dab1lh(mgs,il,ic) =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 1) !dab1(il,ic)
             else ! use lookup table -- not interpolating yet because table resolution of 0.05 is good enough
               i = Nint( alpha(mgs,il)*dqiacralphainv )
               if ( ic == lc .or. ic == li .or. ic == ls .or. (ic == lr .and. imurain == 3) ) then
                 alp = (3.*alpha(mgs,ic) + 2.)
                 j = Nint( (3.*alpha(mgs,ic) + 2.)*dqiacralphainv )
               else ! if ( ic == lr .and. imurain == 1 ) ! rain
                 alp = alpha(mgs,ic)
                 j = Nint( alpha(mgs,ic)*dqiacralphainv )
               endif
             
               dab0lh(mgs,ic,il) = dab0lu(j,i,ic,il)
               dab1lh(mgs,ic,il) = dab1lu(j,i,ic,il)
               dab0lh(mgs,il,ic) = dab0lu(i,j,il,ic)
               dab1lh(mgs,il,ic) = dab1lu(i,j,il,ic)

!               tmp1 = dab0lu(j,i,ic,il)
!               tmp2 = dab1lu(j,i,ic,il)
!               tmp3 = dab0lu(i,j,il,ic)
!               tmp4 = dab1lu(i,j,il,ic)
!               tmp5 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(ic), xmu(il), 0) !dab0(il,ic)
!               tmp6 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(ic), xmu(il), 1) !dab1(il,ic)
!               tmp5 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 0) !dab0(il,ic)
!               tmp6 =  delabk(bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic), 1) !dab1(il,ic)
               
               if ( .false. .and. ny <= 2 ) then
                 write(0,*)
                 write(0,*) 'bb: ', bb(il), bb(ic), xnutmp, xnuc, xmu(il), xmu(ic)
                 write(0,*) 'il,ic = ',il,ic,alpha(mgs,il),i,xnuc,alp,j
                 write(0,*) 'dab0lh,tmp1 = ',dab0lh(mgs,ic,il),tmp1
                 write(0,*) 'dab1lh,tmp2 = ',dab1lh(mgs,ic,il),tmp2
                 write(0,*) 'dab0lh,tmp3 = ',dab0lh(mgs,il,ic),tmp3,tmp5
                 write(0,*) 'dab1lh,tmp4 = ',dab1lh(mgs,il,ic),tmp4,tmp6
               
               endif
             
             endif
             
           endif
          enddo

!          ENDIF
           
             da0lx(mgs,il) = delbk(bb(il), xnutmp, xmu(il), 0)
           if ( il .eq. lh ) then
             da0lh(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
            if ( lzr > 1 ) then
             rzxh(mgs) = 1.
            else
             rzxh(mgs) = ((4. + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr)))
            endif
            
            if ( lzhl < 1 ) then
              rzxhlh(mgs) = rzxhl(mgs)/(((4. + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr))))
            endif
           elseif ( il .eq. lhl ) then
             da0lhl(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
            if ( lzr > 1 ) then
             rzxhl(mgs) = 1.
            else
             rzxhl(mgs) = ((4.0 + alpha(mgs,il))*(5. + alpha(mgs,il))*(6. + alpha(mgs,il))*(1. + xnu(lr)))/   &
     &  ((1. + alpha(mgs,il))*(2. + alpha(mgs,il))*(3. + alpha(mgs,il))*(2. + xnu(lr)))
            endif
           elseif ( il == lr ) then
             xnutmp = (alpha(mgs,il) - 2.)/3.
             da0lr(mgs) = delbk(bb(il), xnutmp, xmu(il), 0)
             da1lr(mgs) = delbk(bb(il), xnutmp, xmu(il), 1)
           endif
          
          endif ! ( qx(mgs,il) > qxmin(il) )
        enddo ! mgs
!          CALL cld_cpu('Z-DELABK')  
        endif ! il /= lr

!         CALL cld_cpu('Z-DELABK')  
        
        endif ! lz(il) .gt. 1
        
        enddo ! il
          
      endif ! ipconc .ge. 6

!      CALL cld_cpu('Z-MOMENT-1')  

!
!  set some values for ice nucleation
!
      do mgs = 1,ngscnt
      kp1 = Min(nz, kgs(mgs)+1 )
!      wvel(mgs) = (0.5)*(w(igs(mgs),jgs,kp1)   &
!     &                  +w(igs(mgs),jgs,kgs(mgs)))

      
        wvelkm1(mgs) = (0.5)*(w(igs(mgs),jgs,kgs(mgs))   &
     &                    +w(igs(mgs),jgs,kgsm(mgs)))
      cninm(mgs) = t7(igs(mgs),jgs,kgsm(mgs))
      cnina(mgs) = t7(igs(mgs),jgs,kgs(mgs))
      cninp(mgs) = t7(igs(mgs),jgs,kgsp(mgs))
      end do

!
!  Set a couple of cloud variables...
!

!      SUBROUTINE setvt(ngscnt,qx,qxmin,cx,rho0,rhovt,xdia,cno,
!     :                 xmas,xdn,xvmn,xvmx,xv,cdx,
!     :                 ipconc,ndebug)
!      SUBROUTINE setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno, &
!     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,            &
!     &                 ipconc1,ndebug1,ngs,nz,kgs,cwnccn,fadvisc,   &
!     &                 cwmasn,cwmasx,cwradn,cnina,cimna,cimxa,      &
!     &                 itype1a,itype2a,temcg,infdo,alpha)


      infdo = 1
      if ( rimdenvwgt > 0 ) infdo = 1

      call setvtz(ngscnt,qx,qxmin,qxw,cx,rho0,rhovt,xdia,cno,cnostmp,   &
     &                 xmas,vtxbar,xdn,xvmn,xvmx,xv,cdx,cdxgs,   &
     &                 ipconc,ndebug,ngs,nz,kgs,fadvisc,   &
     &                 cwmasn,cwmasx,cwradn,cnina,cimn,cimx,   &
     &                 itype1,itype2,temcg,infdo,alpha,0,axx,bxx) ! ,cdh,cdhl)
!     &                 itype1,itype2,temcg,infdo,alpha,0,axh,bxh,axhl,bxhl) ! ,cdh,cdhl)


       if ( lwsm6 .and. ipconc == 0 ) then
         tmp = Max(qxmin(lh), qxmin(ls))
         do mgs = 1,ngscnt
           total = qx(mgs,lh) + qx(mgs,ls)
           if ( total > tmp ) then
             vt2ave(mgs) = (qx(mgs,lh)*vtxbar(mgs,lh,1) + qx(mgs,ls)*vtxbar(mgs,ls,1))/total
           else
             vt2ave(mgs) = 0.0
           endif
         enddo
       endif


!
!  Set number concentrations (need xdia from setvt)
!
      if ( ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: Set concentration'
      if ( ipconc .lt. 1 ) then
         cina(1:ngscnt) = cx(1:ngscnt,li)
      endif
      if ( ipconc .lt. 5 ) then
      do mgs = 1,ngscnt


      if ( ipconc .lt. 3 ) then
!      cx(mgs,lr) = 0.0
      if ( qx(mgs,lr) .gt. qxmin(lh) )  then
!      cx(mgs,lr) = cno(lr)*xdia(mgs,lr,1)
!      xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*cx(mgs,lr))
      end if
      endif

      if ( ipconc .lt. 4 ) then
!      tmp = cx(mgs,ls)
!      cx(mgs,ls) = 0.0
      if ( qx(mgs,ls) .gt. qxmin(ls) )  then
!      cx(mgs,ls) = cno(ls)*xdia(mgs,ls,1)
!      xv(mgs,ls) = rho0(mgs)*qx(mgs,ls)/(xdn(mgs,ls)*cx(mgs,ls))
      end if
      endif ! ( ipconc .lt. 4 )

      if ( ipconc .lt. 5 ) then


!      cx(mgs,lh) = 0.0
      if ( qx(mgs,lh) .gt. qxmin(lh) )  then
!      cx(mgs,lh) = cno(lh)*xdia(mgs,lh,1)
!      xv(mgs,lh) = Max(xvmn(lh), rho0(mgs)*qx(mgs,lh)/(xdn(mgs,lh)*cx(mgs,lh)) )
!      xdia(mgs,lh,3) = (xv(mgs,lh)*6./pi)**(1./3.) 
      end if

      endif ! ( ipconc .lt. 5 )

      end do
      end if
      
      if ( ipconc .ge. 2 ) then
      do mgs = 1,ngscnt
        
        rb(mgs) = 0.5*xdia(mgs,lc,1)*(1./(1.+alpha(mgs,lc)))**(1./6.)
        xl2p(mgs) = Max(0.0d0, 2.7e-2*xdn(mgs,lc)*cx(mgs,lc)*xv(mgs,lc)*   &
     &           ((0.5e20*rb(mgs)**3*xdia(mgs,lc,1))-0.4) )
        if ( rb(mgs) .gt. 3.51e-6 ) then
!          rh(mgs) = Max( 0.5d0*xdia(mgs,lc,1), 6.3d-4/(1.d6*(rb(mgs) - 3.5d-6)) )
          rh(mgs) = Max( 41.d-6, 6.3d-4/(1.d6*(rb(mgs) - 3.5d-6)) )
        else
          rh(mgs) = 41.d-6
        endif
        if ( xl2p(mgs) .gt. 0.0 ) then
          nh(mgs) = 4.2d9*xl2p(mgs)
        else
          nh(mgs) = 1.e30
        endif
      enddo
      endif
      
!
!
!              
!
!  maximum depletion tendency by any one source
!
!
      if( ndebug .ge. 0 ) then
!mpi!        write(0,*) 'Set depletion max/min1'
      endif
      do mgs = 1,ngscnt
      qvimxd(mgs) = 0.70*(qx(mgs,lv)-qis(mgs))*dtpinv ! depletion by all vap. dep to ice.
      
      if ( qx(mgs,lc) < qxmin(lc) ) qvimxd(mgs) = 0.99*(qx(mgs,lv)-qis(mgs))*dtpinv ! this makes virtually no difference whatsoever, but what the heck
      
      qvimxd(mgs) = max(qvimxd(mgs), 0.0)

      frac = 0.1d0
      qimxd(mgs)  = frac*qx(mgs,li)*dtpinv
      qcmxd(mgs)  = frac*qx(mgs,lc)*dtpinv
      qrmxd(mgs)  = frac*qx(mgs,lr)*dtpinv
      qsmxd(mgs)  = frac*qx(mgs,ls)*dtpinv
      qhmxd(mgs)  = frac*qx(mgs,lh)*dtpinv
      if ( lhl > 1 ) qhlmxd(mgs)  = frac*qx(mgs,lhl)*dtpinv
      end do
!
      if( ndebug .ge. 0 ) then
!mpi!        write(0,*) 'Set depletion max/min2'
      endif

      do mgs = 1,ngscnt
!
      if ( qx(mgs,lc) .le. qxmin(lc) ) then
      ccmxd(mgs)  = 0.20*cx(mgs,lc)*dtpinv
      else
      if ( ipconc .ge. 2 ) then
        ccmxd(mgs)  = frac*cx(mgs,lc)*dtpinv
      else
        ccmxd(mgs)  = frac*qx(mgs,lc)/(xmas(mgs,lc)*rho0(mgs)*dtp)
      endif
      end if
!
      if ( qx(mgs,li) .le. qxmin(li) ) then
      cimxd(mgs)  = frac*cx(mgs,li)*dtpinv
      else
      if ( ipconc .ge. 1 ) then
        cimxd(mgs)  = frac*cx(mgs,li)*dtpinv
      else
        cimxd(mgs)  = frac*qx(mgs,li)/(xmas(mgs,li)*rho0(mgs)*dtp)
      endif
      end if
!
!
      crmxd(mgs)  = 0.10*cx(mgs,lr)*dtpinv
      csmxd(mgs)  = frac*cx(mgs,ls)*dtpinv
      chmxd(mgs)  = frac*cx(mgs,lh)*dtpinv

      ccmxd(mgs)  = frac*cx(mgs,lc)*dtpinv
      cimxd(mgs)  = frac*cx(mgs,li)*dtpinv
      crmxd(mgs)  = frac*cx(mgs,lr)*dtpinv
      csmxd(mgs)  = frac*cx(mgs,ls)*dtpinv
      chmxd(mgs)  = frac*cx(mgs,lh)*dtpinv

      qxmxd(mgs,lv) = Max(0.0, 0.1*(qx(mgs,lv) - qvs(mgs))*dtpinv)

      do il = lc,lhab
       qxmxd(mgs,il) = frac*qx(mgs,il)*dtpinv
       cxmxd(mgs,il) = frac*cx(mgs,il)*dtpinv
      enddo

      end do




      if ( ipconc >= 6 ) then
      frac = 0.4d0
      zxmxd(:,:) = 0.0
      do il = lr,lhab
       if ( lz(il) > 0 .or. ( il == lr ) ) then
         do mgs = 1,ngscnt
           zxmxd(mgs,il) = frac*zx(mgs,il)*dtpinv
         enddo
       endif
      enddo
      endif




    ! default factors between mean volume and maximum mass volume
      maxmassfac(lc) = ( (2. + 3.*(1. + xnu(lc)) )**3/( 3.*(1. + xnu(lc)) ) )
      maxmassfac(li) = ( (2. + 3.*(1. + xnu(li)) )**3/( 3.*(1. + xnu(li)) ) )

      if ( imurain == 3 ) then
        maxmassfac(lr) = ( (2. + 3.*(1. + xnu(lr)) )**3/( 3.*(1. + xnu(lr)) ) )
      else
        maxmassfac(lr) =  (3.0 + alphar)**3/    &
     &                  ((3.+alphar)*(2.+alphar)*(1. + alphar) )
      endif

      if ( imusnow == 3 ) then
        maxmassfac(ls) = ( (2. + 3.*(1. + alphas) )**3/( 3.*(1. + alphas) ) )
      else
        maxmassfac(ls) =  (3.0 + alphas)**3/    &
     &                  ((3.+alphas)*(2.+alphas)*(1. + alphas) )
      endif
      
        maxmassfac(lh) =  (3.0 + alphah)**3/    &
     &                  ((3.+alphah)*(2.+alphah)*(1. + alphah) )

       if ( lhl > 1 ) then
        maxmassfac(lhl) =  (3.0 + alphahl)**3/    &
     &                  ((3.+alphahl)*(2.+alphahl)*(1. + alphahl) )
       endif
      


       do mgs = 1,ngscnt
          do il = lh,lhab ! graupel and hail only (and frozen drops)
            
            vshdgs(mgs,il) = vshd ! base value
            
            if ( qx(mgs,il) > qxmin(il) .and. ivshdgs > 0 ) then
              
              ! tmpdiam is weighted diameter of d^(shedalp-1), so for shedalp=3, this is the area-weighted diameter or maximum mass diameter.
              tmpdiam = (shedalp+alpha(mgs,il))*xdia(mgs,il,1) ! *( xdn(mgs,il)/917. )**(1./3.) ! erm added density factor for equiv. solid ice sphere 10.12.2015
              
              if ( tmpdiam > sheddiam0 ) then
                vshdgs(mgs,il) = 0.523599*(1.5e-3)**3/massfacshr ! 1.5mm drops from very large ice
              elseif ( tmpdiam > sheddiam ) then ! intermediate size
                vshdgs(mgs,il) = 0.523599*(3.0e-3)**3/massfacshr ! 3.0mm drops from medium-large ice
              else
!                vshdgs(mgs,il) = Min( xvmx(lr), xv(mgs,il)*xdn(mgs,il)*0.001 ) ! size of drop from melted mean ice particle
                vshdgs(mgs,il) = Min( xvmx(lr), 6./pi*xdn(mgs,il)*0.001*tmpdiam**3 )/massfacshr ! size of drop from melted mean ice particle; 0.001 is 1/rhow
              endif
            endif
          enddo
       enddo

!
!
!  microphysics source terms (1/s) for mixing ratios 
!
!
!
!  Collection efficiencies:
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: Set collection efficiencies'
!
      do mgs = 1,ngscnt
!
!
!
      qcwresv(mgs) = 0.0
      ccwresv(mgs) = 0.0
      
      erw(mgs) = 0.0
      esw(mgs) = 0.0
      ehw(mgs) = 0.0
      efw(mgs) = 0.0
      ehlw(mgs) = 0.0
!      ehxw(mgs) = 0.0
!
      err(mgs) = 0.0
      esr(mgs) = 0.0
      il2(mgs) = 0
      il3(mgs) = 0
      ehr(mgs) = 0.0
      ehlr(mgs) = 0.0
!      ehxr(mgs) = 0.0
!
      eri(mgs) = 0.0
      esi(mgs) = 0.0
      ehi(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehi*ehiclsn
      ehis(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehi*ehiclsn
      ehli(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehli*ehliclsn
      ehlis(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehli*ehliclsn
!      ehxi(mgs) = 0.0
!
      ers(mgs) = 0.0
      ess(mgs) = 0.0
      ehs(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehs*ehsclsn
      ehsfac(mgs) = 1.0 ! factor based on ice saturation
      ehls(mgs) = 0.0 ! used as sticking efficiency, so collection efficiency is ehls*ehlsclsn
      ehscnv(mgs) = 0.0
!      ehxs(mgs) = 0.0
!
      eiw(mgs) = 0.0
      eii(mgs) = 0.0
      ehsclsn(mgs) = 0.0
      ehiclsn(mgs) = 0.0
      ehlsclsn(mgs) = 0.0
      ehliclsn(mgs) = 0.0
      esiclsn(mgs) = 0.0


! reserve droplets
         if ( exwmindiam > 0 .and. qx(mgs,lc) > qxmin(lc) ) then
           tmp = cx(mgs,lc)*Exp(- (exwmindiam/xdia(mgs,lc,1))**3 )
           ccwresv(mgs) =  Min( cx(mgs,lc), Max( 2.e6,  cx(mgs,lc) -  tmp ) )
           
           tmp = cx(mgs,lc) - ccwresv(mgs)

           volt = pi/6.*(exwmindiam)**3
           qcwresv(mgs) = qx(mgs,lc) - tmp*xdn0(lc)*rhoinv(mgs)*(volt + xv(mgs,lc))
           
           
           if ( .false. .and. qx(mgs,lc) > 0.1e-3 ) then
           
             write(0,*) 'cx,qx,crsv,qrsv = ',cx(mgs,lc),qx(mgs,lc),ccwresv(mgs),qcwresv(mgs)
           
           endif

         endif


      icwr(mgs) = 1
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then
       cwrad = 0.5*xdia(mgs,lc,1)
      do il = 1,8
         if ( cwrad .ge. 1.e-6*cwr(il,1) ) icwr(mgs) = il
      enddo
      endif


      irwr(mgs) = 1
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
         rwrad = 0.5*xdia(mgs,lr,3)  ! changed to mean volume diameter (10/6/06)
      do il = 1,6
         if ( rwrad .ge. 1.e-6*grad(il,1) ) irwr(mgs) = il
      enddo
      endif


      igwr(mgs) = 1
!      IF ( qx(mgs,lr) .gt. qxmin(lr) ) THEN
!         rwrad = 0.5*xdia(mgs,lr,1)
! setting erw = 1 always, so now use igwr for graupel
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
         rwrad = 0.5*xdia(mgs,lh,3)  ! changed to mean volume diameter (10/6/06)
      do il = 1,6
         if ( rwrad .ge. 1.e-6*grad(il,1) ) igwr(mgs) = il
      enddo
      endif


      if ( lhl .gt. 1 ) then ! hail is turned on
      ihlr(mgs) = 1
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
         rwrad = 0.5*xdia(mgs,lhl,3)  ! changed to mean volume diameter (10/6/06)
      do il = 1,6
         if ( rwrad .ge. 1.e-6*grad(il,1) ) ihlr(mgs) = il
      enddo
      endif
      endif

!
!
!  Ice-Ice: Collection (cxc) efficiencies
!
!
      if ( qx(mgs,li) .gt. qxmin(li) ) then
!      IF ( ipconc .ge. 14 ) THEN
!       eii(mgs)=0.1*exp(0.1*temcg(mgs))
!       if ( temg(mgs) .lt. 243.15 .and. qx(mgs,lc) .gt. 1.e-6 ) then
!        eii(mgs)=0.1
!       end if
!      
!      ELSE
        eii(mgs) = exp(0.025*Min(temcg(mgs),0.0))  ! alpha1 from LFO83 (21)
!      ENDIF
      if ( temg(mgs) .gt. 273.15 ) eii(mgs) = 1.0
      end if
!
!
!
!  Ice-cloud water: Collection (cxc) efficiencies
!
!
      eiw(mgs) = 0.0
      if ( qx(mgs,li).gt.qxmin(li) .and. qx(mgs,lc).gt.qxmin(lc) ) then
      
      
      if (xdia(mgs,lc,1).gt.ewi_dcmin .and. xdia(mgs,li,1).gt.ewi_dimin) then
! erm 5/10/2007 test following change:
!      if (xdia(mgs,lc,1).gt.12.0e-06 .and. xdia(mgs,li,1).gt.50.0e-06) then
      eiw(mgs) = 0.5
      end if
      if ( temg(mgs) .ge. 273.15 ) eiw(mgs) = 0.0
      end if

!
!
!
!  Rain: Collection (cxc) efficiencies
!
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lc).gt.qxmin(lc) ) then

       if ( lnr .gt. 1 ) then
       erw(mgs) = 1.0

       else

!      cwrad = 0.5*xdia(mgs,lc,1)
!      erw(mgs) =
!     >  min((aradcw + cwrad*(bradcw + cwrad*
!     <  (cradcw + cwrad*(dradcw)))), 1.0)
!       IF ( xdia(mgs,lc,1) .lt. 2.4e-06 .or. xdia(mgs,lr,1) .le. 50.0e-6 ) THEN
!          erw(mgs)=0.0
!       ENDIF
!       erw(mgs) = ew(icwr(mgs),igwr(mgs))
! interpolate along droplet radius
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = irwr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,3)
       rwrad = 0.5*xdia(mgs,lr,3)
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)

!       write(iunit,*) 'slop1: ',slope1,slope2,ew(ic,ir),cwr(ic,2)

       x1 = ew(ic,  ir) + slope1*Max(0.0, (cwrad - cwr(ic,1)) )
       x2 = ew(icp1,ir) + slope2*Max(0.0, (cwrad - cwr(ic,1)) )

       slope1 = (x2 - x1)*grad(ir,2)

       erw(mgs) = Max(0.0, x1 + slope1*Max(0.0, (rwrad - grad(ir,1)) ))

!       write(iunit,*) 'erw: ',erw(mgs),1.e6*cwrad,1.e6*rwrad,ic,ir,x1,x2
!       write(iunit,*)

       erw(mgs) = Max(0.0, erw(mgs) )
       if ( rwrad .lt. 50.e-6 ) then
         erw(mgs) = 0.0
       elseif (  rwrad .lt. 100.e-6 ) then  ! linear change from zero at 50 to erw at 100 microns
         erw(mgs) = erw(mgs)*(rwrad - 50.e-6)/50.e-6
       endif

       endif
      end if
      if ( cx(mgs,lc) .le. 0.0 ) erw(mgs) = 0.0
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lr).gt.qxmin(lr) ) then
      err(mgs)=1.0
      end if
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,ls).gt.qxmin(ls) ) then
      ers(mgs)=1.0
      end if
!
      if ( qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,li).gt.qxmin(li) ) then
!        IF ( vtxbar(mgs,lr,1) .gt. vtxbar(mgs,li,1) .and.
!     :       xdia(mgs,lr,3) .gt. 200.e-6 .and. xdia(mgs,li,3) .gt. 100.e-6 ) THEN
         eri(mgs) = eri0
!      cwrad = 0.5*xdia(mgs,li,3)
!      eri(mgs) =
!     >  1.0*min((aradcw + cwrad*(bradcw + cwrad*
!     <  (cradcw + cwrad*(dradcw)))), 1.0)
!         ENDIF
!       if ( xdia(mgs,li,1) .lt. 10.e-6 ) eri(mgs)=0.0
       if ( xdia(mgs,li,3) .lt. eri_cimin ) eri(mgs)=0.0
      end if
!
!
!  Snow aggregates: Collection (cxc) efficiencies
!
! Modified by ERM with a linear function for small droplets and large
! snow agg. based numerical data from Wang and Ji (1992) in P&K 1997 (Fig. 14-13), which
! allows collection of very small droplets, albeit at low efficiency.  But slow
! fall speeds of snow make up for the efficiency.
!
      esw(mgs) = 0.0
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lc).gt.qxmin(lc) ) then
        esw(mgs) = 0.5
        if ( xdia(mgs,lc,1) .gt. 15.e-6 .and. xdia(mgs,ls,1) .gt. 100.e-6) then
          esw(mgs) = 0.5
        elseif ( xdia(mgs,ls,1) .ge. 500.e-6 ) then
          esw(mgs) = Min(0.5, 0.05 + (0.8-0.05)/(40.e-6)*xdia(mgs,lc,1) )
        endif
      end if
!
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,lr).gt.qxmin(lr)  &
     &     .and. temg(mgs) .lt. tfr - 1.   &
     &                               ) then
      esr(mgs)=Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,ls,1))
      if ( qx(mgs,ls) < 1.e-4 .and. qx(mgs,lr) < 1.e-4 ) il2(mgs) = 1
      end if
      
      if ( ipconc < 3 .and. temg(mgs) < tfr .and. qx(mgs,lr).gt.qxmin(lr) .and. qx(mgs,lr) < 1.e-4 ) then
        il3(mgs) = 1
      endif
!
!      if ( qx(mgs,ls).gt.qxmin(ls) ) then
      if ( temcg(mgs) < 0.0 ) then
            
      if ( ipconc .lt. 4 .or. temcg(mgs) < esstem1 ) then
        ess(mgs) = 0.0
!        ess(mgs)=0.1*exp(0.1*min(temcg(mgs),0.0))
!        ess(mgs)=min(0.1,ess(mgs))
      
      else
      
        fac = Abs(ess0)
        if ( iessopt == 2 ) then ! experimental code
!         IF ( wvel(mgs) > 2.0 .or. wvel(mgs) < -0.5 .or. ssi(mgs) < 1.0 ) THEN
         if ( wvel(mgs) > 2.0 ) then
          ! assume convective cell or downdraft
           fac = 0.0
         elseif ( wvel(mgs) > 1.0 ) then ! transition to stratiform range of values
           fac = Max(0.0, 2.0 - wvel(mgs))*fac
         endif
        elseif ( iessopt == 3 ) then ! factor based on ice supersat
           if ( ssi(mgs) <= 1.0 ) then
             fac = 0.0
             ehsfac(mgs) = 0.0
           elseif ( ssi(mgs) <= 1.02 ) then
             fac = fac*(ssi(mgs) - 1.0)/0.02
             ehsfac(mgs) = (ssi(mgs) - 1.0)/0.02
           endif
        elseif ( iessopt == 4 ) then ! factor based on ice supersat
           if ( ssi(mgs) <= 1.0 ) then
             fac = 0.1
             ehsfac(mgs) = 0.1
           elseif ( ssi(mgs) <= 1.005 ) then
             fac = Max(0.1, fac*(ssi(mgs) - 1.0)/0.005)
             ehsfac(mgs) = Max(0.1, (ssi(mgs) - 1.0)/0.005)
           endif
        endif
        
        if ( temcg(mgs) > esstem1 .and. temcg(mgs) < esstem2 ) then  ! only nonzero for T > esstem1
          ess(mgs) = fac*Exp(ess1*(esstem2) )*(temcg(mgs) - esstem1)/(esstem2 - esstem1) ! linear ramp up from zero at esstem1 to value at esstem2
        elseif ( temcg(mgs) >= esstem2 ) then
          ess(mgs) = fac*Exp(ess1*Min( temcg(mgs), 0.0 ) )
        endif
        
      endif
      end if
!
      if ( qx(mgs,ls).gt.qxmin(ls) .and. qx(mgs,li).gt.qxmin(li) ) then
       esiclsn(mgs) = esi_collsn
!      IF ( ipconc .lt. 4 ) THEN
      if ( ipconc < 1 .and. lwsm6 ) then
        esi(mgs) = exp(0.7*min(temcg(mgs),0.0))
      else
        esi(mgs) = esi0*exp(0.1*min(temcg(mgs),0.0))
        esi(mgs) = Min(0.1,esi(mgs))
      endif
      if ( ipconc .le. 3 ) then
       esi(mgs) =  exp(0.025*min(temcg(mgs),0.0)) ! LFO
!       esi(mgs) =  Min(0.5, exp(0.025*min(temcg(mgs),0.0)) ) ! LFO
!       esi(mgs)=0.5*exp(0.1*min(temcg(mgs),0.0))  ! 10ice
      endif
!      ELSE ! zrnic/ziegler 1993
!      esi(mgs)= 0.1 ! 0.5*exp(0.1*min(temcg(mgs),0.0))
!      ENDIF
      if ( temg(mgs) .gt. 273.15 ) esi(mgs) = 0.0
      end if
!
!
!
!
!  Graupel: Collection (cxc) efficiencies
!
!
       xmascw(mgs) = xmas(mgs,lc)
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc).gt.qxmin(lc) ) then !{
       ehw(mgs) = 1.0
       if ( iehw .eq. 0 ) then
       ehw(mgs) = ehw0  ! default value is 1.0
       elseif ( iehw .eq. 1 .or. iehw .eq. 10 ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      ehw(mgs) = Min( ehw0,    &
     &  ewfac*min((aradcw + cwrad*(bradcw + cwrad*   &
     &  (cradcw + cwrad*(dradcw)))), 1.0) )
      
       elseif ( iehw .eq. 2 .or. iehw .eq. 10 ) then
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = igwr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,1)
       rwrad = 0.5*xdia(mgs,lh,3)  ! changed to mean volume diameter
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)
 
!        write(iunit,*) 'slop1: ',slope1,slope2,ew(ic,ir),cwr(ic,2)

       x1 = ew(ic,  ir) + slope1*Max(0.0, (cwrad - cwr(ic,1)) )
       x2 = ew(icp1,ir) + slope2*Max(0.0, (cwrad - cwr(ic,1)) )
       
       slope1 = (x2 - x1)*grad(ir,2)
       
       tmp = Max( 0.0, Min( 1.0, x1 + slope1*Max(0.0, (rwrad - grad(ir,1)) ) ) )
       ehw(mgs) = Min( ehw(mgs), tmp )

!       write(iunit,*) 'ehw: ',ehw(mgs),1.e6*cwrad,1.e6*rwrad,ic,ir,x1,x2
!       write(iunit,*)

!       ehw(mgs) = Max( 0.2, ehw(mgs) )
!  assume that ehw = 1 for zero air resistance (rho0 = 0.0) and extrapolate toward that
!      ehw(mgs) = ehw(mgs) + (ehw(mgs) - 1.0)*(rho0(mgs) - rho00)/rho00
!      ehw(mgs) = ehw(mgs) + (1.0 - ehw(mgs))*((Max(0.0,rho00 - rho0(mgs)))/rho00)**2

       elseif ( iehw .eq. 3 .or. iehw .eq. 10 ) then ! use fraction of droplets greater than dmincw diameter
         tmp = Exp(- (dmincw/xdia(mgs,lc,1))**3)
         xmascw(mgs) = xmas(mgs,lc) + xdn0(lc)*(pi*dmincw**3/6.0) ! this is the average mass of the droplets with d > dmincw
         ehw(mgs) = Min( ehw(mgs), tmp )
       elseif ( iehw .eq. 4 .or. iehw .eq. 10 ) then ! Cober and List 1993, eq. 19-20
         tmp =  &
     &   2.0*xdn(mgs,lc)*vtxbar(mgs,lh,1)*(0.5*xdia(mgs,lc,1))**2 &
     &  /(9.0*fadvisc(mgs)*0.5*xdia(mgs,lh,3))
         tmp = Max( 1.5, Min(10.0, tmp) )
         ehw(mgs) = Min( ehw(mgs), 0.55*Log10(2.51*tmp) )
       endif
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) ehw(mgs)=0.0

       ehw(mgs) = Min( ehw0, ehw(mgs) )
       
       if ( ibfc == -1 .and. temcg(mgs) < -41.0 ) then
        ehw(mgs) = 0.0
       endif 

      end if !}
!
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lr).gt.qxmin(lr)    &
!     &     .and. temg(mgs) .lt. tfr    &
     &                               ) then
!      ehr(mgs) = Exp(-(40.e-6)**3/xv(mgs,lr))*Exp(-40.e-6/xdia(mgs,lh,1))
!      ehr(mgs) = 1.0
       ehr(mgs) = Exp(-(40.e-6)/xdia(mgs,lr,3))*Exp(-40.e-6/xdia(mgs,lh,3))
       ehr(mgs) = Min( ehr0, ehr(mgs) )
      end if
!
      if ( qx(mgs,ls).gt.qxmin(ls) ) then
        if ( ipconc .ge. 4 ) then
        ehscnv(mgs) = ehs0*exp(ehs1*min(temcg(mgs),0.0)) ! for 2-moment, used as default for ehs and ehls. Otherwise not used for snow->graupel conversion
        else
        ehscnv(mgs) = exp(0.09*min(temcg(mgs),0.0))
        endif
        
        if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc) >= qxmin(lc)  ) then
!          ehsclsn(mgs) = ehs_collsn
!          ehs(mgs) = ehscnv(mgs)*ehsfac(mgs)*Min(1.0, Max(0.0,xdn(mgs,lh) - 300.)/300.  )
!        ELSEIF ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lc) >= qxmin(lc)  ) then
          ehsclsn(mgs) = ehs_collsn
          if ( xdia(mgs,ls,3) < 40.e-6 ) then
            ehsclsn(mgs) = 0.0
          elseif ( xdia(mgs,ls,3) < 150.e-6 ) then
            ehsclsn(mgs) =  ehs_collsn*(xdia(mgs,ls,3) - 40.e-6)/(150.e-6 - 40.e-6)
          else
            ehsclsn(mgs) = ehs_collsn
          endif
!          ehs(mgs) = ehscnv(mgs)*Min(1.0, Max(0., xdn(mgs,lh) - xdnmn(lh)*1.2)/xdnmn(lh)  ) ! shut off qhacs as graupel goes to lowest density
          ehs(mgs) = ehscnv(mgs)*Min(1.0, Max(0.0,xdn(mgs,lh) - 300.)/300.  ) ! shut off qhacs as graupel goes to low density; limits scavenging of snow in bright band
!          ehs(mgs) = ehscnv(mgs) ! *Min(1.0, Max(0.0,xdn(mgs,lh) - 300.)/300.  ) ! shut off qhacs as graupel goes to low density
          ehs(mgs) = Min(ehs(mgs),ehsmax)
        end if
      endif
!
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,li).gt.qxmin(li) ) then
      ehiclsn(mgs) = ehi_collsn
      ehi(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      ehi(mgs) = Min( ehimax, Max( ehi(mgs), ehimin ) )
!      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehi(mgs) = 0.0
      end if

      if ( lis > 1 ) then
      if ( qx(mgs,lh).gt.qxmin(lh) .and. qx(mgs,lis).gt.qxmin(lis) ) then
      ehisclsn(mgs) = ehi_collsn
      ehis(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      ehis(mgs) = Min( ehimax, Max( ehis(mgs), ehimin ) )
!      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehis(mgs) = 0.0
      end if
      endif

 
!
!
!  Hail: Collection (cxc) efficiencies
!
!
      if ( lhl .gt. 1 ) then

      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,lc).gt.qxmin(lc) ) then
       if ( iehw == 3 ) iehlw = 3
       if ( iehw == 4 ) iehlw = 4
       ehlw(mgs) = ehlw0
       if ( iehlw .eq. 0 ) then
       ehlw(mgs) = ehlw0  ! default value is 1.0
       elseif ( iehlw .eq. 1 .or. iehlw .eq. 10 ) then
      cwrad = 0.5*xdia(mgs,lc,1)
      ehlw(mgs) = Min( ehlw0,    &
     &  ewfac*min((aradcw + cwrad*(bradcw + cwrad*   &
     &  (cradcw + cwrad*(dradcw)))), 1.0) )
      
       elseif ( iehlw .eq. 2 .or. iehlw .eq. 10 ) then
       ic = icwr(mgs)
       icp1 = Min( 8, ic+1 )
       ir = ihlr(mgs)
       irp1 = Min( 6, ir+1 )
       cwrad = 0.5*xdia(mgs,lc,1)
       rwrad = 0.5*xdia(mgs,lhl,3)  ! changed to mean volume diameter
       
       slope1 = (ew(icp1, ir  ) - ew(ic,ir  ))*cwr(ic,2)
       slope2 = (ew(icp1, irp1) - ew(ic,irp1))*cwr(ic,2)
       
       x1 = ew(ic,  ir) + slope1*(cwrad - cwr(ic,1))
       x2 = ew(icp1,ir) + slope2*(cwrad - cwr(ic,1))
       
       slope1 = (x2 - x1)*grad(ir,2)
       
       tmp = Max( 0.0, Min( 1.0, x1 + slope1*(rwrad - grad(ir,1)) ) )
         ehlw(mgs) = Min( ehlw(mgs), tmp )
       ehlw(mgs) = Min( ehlw0, ehlw(mgs) )
!       ehw(mgs) = Max( 0.2, ehw(mgs) )
!  assume that ehw = 1 for zero air resistance (rho0 = 0.0) and extrapolate toward that
!      ehw(mgs) = ehw(mgs) + (ehw(mgs) - 1.0)*(rho0(mgs) - rho00)/rho00
!      ehlw(mgs) = ehlw(mgs) + (1.0 - ehlw(mgs))*((Max(0.0,rho00 - rho0(mgs)))/rho00)**2

       elseif ( iehlw .eq. 3 .or. iehlw .eq. 10 ) then ! use fraction of droplets greater than 15 micron diameter
         tmp = Exp(- (dmincw/xdia(mgs,lc,1))**3)
         ehlw(mgs) = Min( ehlw(mgs), tmp )
       elseif ( iehlw .eq. 4 .or. iehlw .eq. 10 ) then ! Cober and List 1993
         tmp =  &
     &   2.0*xdn(mgs,lc)*vtxbar(mgs,lhl,1)*(0.5*xdia(mgs,lc,1))**2 &
     &  /(9.0*fadvisc(mgs)*0.5*xdia(mgs,lhl,3))
         tmp = Max( 1.5, Min(10.0, tmp) )
         ehlw(mgs) = Min( ehlw(mgs), 0.55*Log10(2.51*tmp) )
       endif
      if ( xdia(mgs,lc,1) .lt. 2.4e-06 ) ehlw(mgs)=0.0
       ehlw(mgs) = Min( ehlw0, ehlw(mgs) )

       if ( ibfc == -1 .and. temcg(mgs) < -41.0 ) then 
        ehlw(mgs) = 0.0
       endif 

      end if
!
      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,lr).gt.qxmin(lr)    &
!     &     .and. temg(mgs) .lt. tfr    &
     &                               ) then
        ehlr(mgs) = 1.0
       ehlr(mgs) = Min( ehlr0, ehlr(mgs) )
      end if
!
      if ( qx(mgs,ls).gt.qxmin(ls) ) then
        if ( qx(mgs,lhl).gt.qxmin(lhl)  ) then
          ehlsclsn(mgs) = ehls_collsn
          ehls(mgs) = ehscnv(mgs)
          ehls(mgs) = Min(ehls(mgs),ehsmax)
        end if
      endif
!
      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,li).gt.qxmin(li) ) then
      ehliclsn(mgs) = ehli_collsn
      ehli(mgs)=eii0hl*exp(eii1hl*min(temcg(mgs),0.0))
      ehli(mgs) = Min( ehimax, Max( ehli(mgs), ehimin ) )
      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehli(mgs) = 0.0
      end if

      if ( lis > 1 ) then
      if ( qx(mgs,lhl).gt.qxmin(lhl) .and. qx(mgs,lis).gt.qxmin(lis) ) then
      ehlisclsn(mgs) = ehli_collsn
      ehlis(mgs)=eii0*exp(eii1*min(temcg(mgs),0.0))
      ehlis(mgs) = Min( ehimax, Max( ehlis(mgs), ehimin ) )
      if ( temg(mgs) .gt. 273.15 .or. ( qx(mgs,lc) < qxmin(lc)) ) ehlis(mgs) = 0.0
      end if
      endif


      endif ! lhl .gt. 1

      enddo  ! mgs loop for collection efficiencies

!
!
!
!  Set flags for plates vs. columns
!
!
      do mgs = 1,ngscnt
!
      xplate(mgs) = 0.0
      xcolmn(mgs) = 1.0
!
!      if ( temcg(mgs) .lt. 0. .and. temcg(mgs) .ge. -4. ) then
!      xplate(mgs) = 1.0
!      xcolmn(mgs) = 0.0
!      end if
!c
!      if ( temcg(mgs) .lt. -4. .and. temcg(mgs) .ge. -9. ) then
!      xplate(mgs) = 0.0
!      xcolmn(mgs) = 1.0
!      end if
!c
!      if ( temcg(mgs) .lt. -9. .and. temcg(mgs) .ge. -22.5 ) then
!      xplate(mgs) = 1.0
!      xcolmn(mgs) = 0.0
!      end if
!c
!      if ( temcg(mgs) .lt. -22.5 .and. temcg(mgs) .ge. -90. ) then
!      xplate(mgs) = 0.0
!      xcolmn(mgs) = 1.0
!      end if
!
      end do
      
      

!
!
!
!  Collection growth equations....
!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: rain collects xxxxx'
!
      do mgs = 1,ngscnt
      qracw(mgs) =  0.0
      if ( qx(mgs,lr) .gt. qxmin(lr) .and. erw(mgs) .gt. 0.0 ) then
      if ( ipconc .lt. 3 ) then
       if ( erw(mgs) .gt. 0.0 .and. qx(mgs,lr) .gt. 1.e-7 ) then
       vt = (ar*(xdia(mgs,lc,1)**br))*rhovt(mgs)
       qracw(mgs) =    &
     &   (0.25)*pi*erw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lr) &
!     >  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))   &
     &  *Max(0.0, vtxbar(mgs,lr,1)-vt)   &
     &  *(  gf3*xdia(mgs,lr,2)    &
     &    + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,lc,1)    &
     &    + gf1*xdia(mgs,lc,2) )
!       qracw(mgs) = 0.0
!      write(iunit,*) 'qracw,cx =',qracw(mgs),1.e6*xdia(mgs,lr,1),erw(mgs)
!      write(iunit,*) 'qracw,cx =',qracw(mgs),cx(mgs,lc),kgs(mgs),cx(mgs,lr),1.e6*xdia(mgs,lr,1),vtxbar(mgs,lr,1),vt
!      write(iunit,*) 'vtr: ',vtxbar(mgs,lr,1), ar*gf4br/6.0*xdia(mgs,lr,1)**br, rhovt(mgs),
!     :         ar*gf4br/6.0*xdia(mgs,lr,1)**br * rhovt(mgs)
       endif
      else

      if ( dmrauto <= 0 .or.  rho0(mgs)*qx(mgs,lr) > 1.2*xl2p(mgs) ) then 
       rwrad = 0.5*xdia(mgs,lr,3)
        if ( rwrad .gt. rh(mgs) ) then ! .or. cx(mgs,lr) .gt. nh(mgs) ) then
         if ( rwrad .gt. rwradmn ) then
!      DM1CCC=A2*XNC*XNR*XVC*(((CNU+2.)/(CNU+1.))*XVC+XVR)       ! (A12)
!     NOTE: Result is independent of imurain, assumes mucloud = 3
           qracw(mgs) = erw(mgs)*aa2*cx(mgs,lr)*cx(mgs,lc)*xmas(mgs,lc)*   &
     &        ((alpha(mgs,lc) + 2.)*xv(mgs,lc)/(alpha(mgs,lc) + 1.) + xv(mgs,lr))/rho0(mgs) !*rhoinv(mgs)
         else

          if ( imurain == 3 ) then

!      DM1CCC=A1*XNC*XNR*(((CNU+3.)*(CNU+2.)/(CNU+1.)**2)*XVC**3+ ! (A14)
!     1 ((RNU+2.)/(RNU+1.))*XVC*XVR**2)

!           qracw(mgs) = aa1*cx(mgs,lr)*cx(mgs,lc)*xdn(mgs,lc)*   &
!     &        ((cnu + 3.)*(cnu + 2.)*xv(mgs,lc)**3/(cnu + 1.)**2 +    &
!     &         (alpha(mgs,lr) + 2.)*xv(mgs,lc)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.))/rho0(mgs) !*rhoinv(mgs)
! save multiplies by converting cx*xdn*xv/rho0 to qx
           qracw(mgs) = aa1*cx(mgs,lr)*(qx(mgs,lc)-qcwresv(mgs))*   &
     &        ((alpha(mgs,lc) + 3.)*(alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.)**2 +    &
     &         (alpha(mgs,lr) + 2.)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.)) 
           
           else ! imurain == 1

           qracw(mgs) = aa1*cx(mgs,lr)*(qx(mgs,lc)-qcwresv(mgs))*   &
     &        ((alpha(mgs,lc) + 3.)*(alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.)**2 +    &
     &         (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)*xv(mgs,lr)**2/ &
     &          ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.))) 
           
           endif
           
         endif
        endif
        endif
       endif
!       qracw(mgs) = Min(qracw(mgs), qx(mgs,lc))
       qracw(mgs) = Min(qracw(mgs), qcmxd(mgs))
       endif
      end do
!
      do mgs = 1,ngscnt
      qraci(mgs) = 0.0
      craci(mgs) = 0.0
      qracs(mgs) = 0.0
      if ( eri(mgs) .gt. 0.0 .and. iacr .ge. 1 .and. xdia(mgs,lr,3) .gt. 2.*rwradmn ) then
        if ( ipconc .ge. 3 ) then

           tmp = eri(mgs)*aa2*cx(mgs,lr)*cx(mgs,li)*   &
     &        ((cinu + 2.)*xv(mgs,li)/(cinu + 1.) + xv(mgs,lr))

        qraci(mgs) = Min( qxmxd(mgs,li), tmp*xmas(mgs,li)*rhoinv(mgs) )
        craci(mgs) = Min( cxmxd(mgs,li), tmp )

!       vt = Sqrt((vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))**2 +
!     :            0.04*vtxbar(mgs,lr,1)*vtxbar(mgs,li,1) )
!
!          qraci(mgs) = 0.25*pi*eri(mgs)*cx(mgs,lr)*qx(mgs,li)*vt*
!     :         (  da0(lr)*xdia(mgs,lr,3)**2 +
!     :            dab1(lr,li)*xdia(mgs,lr,3)*xdia(mgs,li,3) +
!     :            da1(li)*xdia(mgs,li,3)**2 )
!
!
!       vt = Sqrt((vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))**2 +
!     :            0.04*vtxbar(mgs,lr,1)*vtxbar(mgs,li,1) )
!
!          craci(mgs) = 0.25*pi*eri(mgs)*cx(mgs,lr)*cx(mgs,li)*vt*
!     :         (  da0(lr)*xdia(mgs,lr,3)**2 +
!     :            dab0(lr,li)*xdia(mgs,lr,3)*xdia(mgs,li,3) +
!     :            da0(li)*xdia(mgs,li,3)**2 )
!
!          qraci(mgs) = Min( qraci(mgs), qxmxd(mgs,li) )
!          craci(mgs) = Min( craci(mgs), cxmxd(mgs,li) )

        else
          qraci(mgs) =    &
     &     min(   &
     &     (0.25)*pi*eri(mgs)*qx(mgs,li)*cx(mgs,lr)   &
     &    *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))   &
     &    *(  gf3*xdia(mgs,lr,2)    &
     &      + 2.0*gf2*xdia(mgs,lr,1)*xdia(mgs,li,1)    &
     &      + gf1*xdia(mgs,li,2) )     &
     &    , qimxd(mgs))
        endif
      if ( temg(mgs) .gt. 268.15 ) then
      qraci(mgs) = 0.0
      end if
      endif
      end do
!
      if ( ipconc < 3 ) then
      do mgs = 1,ngscnt
      qracs(mgs) = 0.0
      if ( ers(mgs) .gt. 0.0 .and. ipconc < 3 ) then
       if ( lwsm6 .and. ipconc == 0 ) then
         vt = vt2ave(mgs)
       else
         vt = vtxbar(mgs,ls,1)
       endif
      qracs(mgs) =      &
     &   min(     &
     &   ((0.25)*pi/gf4)*ers(mgs)*qx(mgs,ls)*cx(mgs,lr)     &
     &  *abs(vtxbar(mgs,lr,1)-vt)     &
     &  *(  gf6*gf1*xdia(mgs,ls,2)     &
     &    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lr,1)      &
     &    + gf4*gf3*xdia(mgs,lr,2) )      &
     &  , qsmxd(mgs))
      endif
      end do
      endif

!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: snow collects xxxxx'
!
      do mgs = 1,ngscnt
      qsacw(mgs) =  0.0
      csacw(mgs) =  0.0
      vsacw(mgs) =  0.0
      if ( esw(mgs) .gt. 0.0 ) then

       if ( ipconc .ge. 4 ) then
!      QSACC=CECS*RVT*A2*XNC*XNS*XVC*ROS*
!     *    (((CNU+2.)/(CNU+1.))*XVC+XVS)/RO

!        tmp = esw(mgs)*rvt*aa2*cx(mgs,ls)*cx(mgs,lc)*
!     :        ((cnu + 2.)*xv(mgs,lc)/(cnu + 1.) + xv(mgs,ls))
        tmp = 1.0*rvt*aa2*cx(mgs,ls)*cx(mgs,lc)*   &
     &        ((alpha(mgs,lc) + 2.)*xv(mgs,lc)/(alpha(mgs,lc) + 1.) + xv(mgs,ls))

        qsacw(mgs) = Min( qxmxd(mgs,lc), tmp*xmas(mgs,lc)*rhoinv(mgs) )
        csacw(mgs) = Min( cxmxd(mgs,lc), tmp )

          if ( lvol(ls) .gt. 1 ) then
             if ( temg(mgs) .lt. 273.15) then
             rimdn(mgs,ls) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,ls,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
             rimdn(mgs,ls) = Min( Max( rimc3, rimdn(mgs,ls) ), rimc4 )
             else
             rimdn(mgs,ls) = 1000.
             endif

           vsacw(mgs) = rho0(mgs)*qsacw(mgs)/rimdn(mgs,ls)

          endif


!        qsacw(mgs) = cecs*aa2*cx(mgs,ls)*cx(mgs,lc)*xmas(mgs,lc)*
!     :        ((alpha(mgs,lc) + 2.)*xv(mgs,lc)/(alpha(mgs,lc) + 1.) + xv(mgs,ls))*rhoinv(mgs)
       else
!      qsacw(mgs) =
!     >   min(
!     >   ((0.25)*pi)*esw(mgs)*qx(mgs,lc)*cx(mgs,ls)
!     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))
!     >  *(  gf3*xdia(mgs,ls,2)
!     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,lc,1)
!     >    + gf1*xdia(mgs,lc,2) )
!     <  , qcmxd(mgs))

            vt = abs(vtxbar(mgs,ls,1)-vtxbar(mgs,lc,1))

          qsacw(mgs) = 0.25*pi*esw(mgs)*cx(mgs,ls)*qx(mgs,lc)*vt*   &
     &         (  da0(ls)*xdia(mgs,ls,3)**2 +     &
     &            dab1(ls,lc)*xdia(mgs,ls,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 )
        qsacw(mgs) = Min( qsacw(mgs), qxmxd(mgs,ls) )
        csacw(mgs) = rho0(mgs)*qsacw(mgs)/xmas(mgs,lc)
       endif
      endif
      end do
!
!
      do mgs = 1,ngscnt
      qsaci(mgs) = 0.0
      csaci(mgs) = 0.0
      csaci0(mgs) = 0.0
      if ( ipconc .ge. 4 ) then
      if ( esi(mgs) .gt. 0.0 .or. ( ipelec > 0 .and. esiclsn(mgs) > 0.0 )) then
!      QSCOI=CEXS*RVT*A2*XNCI*XNS*XVCI*ROS*
!     *  (((CINU+2.)/(CINU+1.))*VCIP+XVS)/RO

        tmp = esiclsn(mgs)*rvt*aa2*cx(mgs,ls)*cx(mgs,li)*   &
     &        ((cinu + 2.)*xv(mgs,li)/(cinu + 1.) + xv(mgs,ls))

        qsaci(mgs) = Min( qxmxd(mgs,li), esi(mgs)*tmp*xmas(mgs,li)*rhoinv(mgs) )
        csaci0(mgs) = tmp
        csaci(mgs) = Min(cxmxd(mgs,li), esi(mgs)*tmp )

!      qsaci(mgs) =
!     >   min(
!     >   ((0.25)*pi)*esi(mgs)*qx(mgs,li)*cx(mgs,ls)
!     >  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))
!     >  *(  gf3*xdia(mgs,ls,2)
!     >    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,li,1)
!     >    + gf1*xdia(mgs,li,2) )
!     <  , qimxd(mgs))
      endif
      else ! 
      if ( esi(mgs) .gt. 0.0 ) then
         qsaci(mgs) =    &
     &   min(   &
     &   ((0.25)*pi)*esi(mgs)*qx(mgs,li)*cx(mgs,ls)   &
     &  *abs(vtxbar(mgs,ls,1)-vtxbar(mgs,li,1))   &
     &  *(  gf3*xdia(mgs,ls,2)    &
     &    + 2.0*gf2*xdia(mgs,ls,1)*xdia(mgs,li,1)    &
     &    + gf1*xdia(mgs,li,2) )     &
     &  , qimxd(mgs))
      endif
      endif
      end do
!
!
!
      do mgs = 1,ngscnt
      qsacr(mgs) = 0.0
      qsacrs(mgs) = 0.0
      csacr(mgs) = 0.0
      if ( esr(mgs) .gt. 0.0 ) then
      if ( ipconc .ge. 3 ) then
!       vt = Sqrt((vtxbar(mgs,ls,1)-vtxbar(mgs,lr,1))**2 + 
!     :            0.04*vtxbar(mgs,ls,1)*vtxbar(mgs,lr,1) )
!       qsacr(mgs) = esr(mgs)*cx(mgs,ls)*vt*
!     :     qx(mgs,lr)*0.25*pi*
!     :      (3.02787*xdia(mgs,lr,2) + 
!     :       3.30669*xdia(mgs,ls,1)*xdia(mgs,lr,1) + 
!     :       2.*xdia(mgs,ls,2))
!        qsacr(mgs) = Min( qsacr(mgs), qrmxd(mgs) )
!        csacr(mgs) = qsacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
!        csacr(mgs) = min(csacr(mgs),crmxd(mgs))
      else
       if ( lwsm6 .and. ipconc == 0 ) then
         vt = vt2ave(mgs)
       else
         vt = vtxbar(mgs,ls,1)
       endif
       
       qsacr(mgs) =   &
     &   min(   &
     &   ((0.25)*pi/gf4)*esr(mgs)*qx(mgs,lr)*cx(mgs,ls)   &
     &  *abs(vtxbar(mgs,lr,1)-vt)   &
     &  *(  gf6*gf1*xdia(mgs,lr,2)   &
     &    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,ls,1)    &
     &    + gf4*gf3*xdia(mgs,ls,2) )    &
     &  , qrmxd(mgs))
      endif
      endif
      end do
!
!
!

      if (ndebug .gt. 0 ) write(0,*) 'Collection: graupel collects xxxxx'
!
      do mgs = 1,ngscnt
      qhacw(mgs) = 0.0
      qhacwmlr(mgs) = 0.0
      rarx(mgs,lh) = 0.0
      vhacw(mgs) = 0.0
      vhsoak(mgs) = 0.0
      zhacw(mgs) = 0.0
      
      if ( .false. ) then
        vtmax = (gz(igs(mgs),jgs,kgs(mgs))*dtpinv)
        vtxbar(mgs,lh,1) = Min( vtmax, vtxbar(mgs,lh,1))
        vtxbar(mgs,lh,2) = Min( vtmax, vtxbar(mgs,lh,2))
        vtxbar(mgs,lh,3) = Min( vtmax, vtxbar(mgs,lh,3))
      endif
      if ( ehw(mgs) .gt. 0.0 ) then

        if ( ipconc .ge. 2 ) then

        if ( .false. ) then  
        qhacw(mgs) = (ehw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lh)*pi*   &
     &    abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))*   &
     &    (2.0*xdia(mgs,lh,1)*(xdia(mgs,lh,1) +    &
     &         xdia(mgs,lc,1)*gf73rds) +    &
     &      xdia(mgs,lc,2)*gf83rds))/4.     
     
         else  ! using Seifert coefficients
            vt = abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1)) 

          qhacw(mgs) = 0.25*pi*ehw(mgs)*cx(mgs,lh)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,lc)*xdia(mgs,lh,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 ) 
         
         endif
          qhacw(mgs) = Min( qhacw(mgs), 0.5*qx(mgs,lc)*dtpinv )
        
         if ( lzh .gt. 1 ) then
          tmp = qx(mgs,lh)/cx(mgs,lh)
          
!!          g1 = (6.0 + alpha(mgs,lh))*(5.0 + alpha(mgs,lh))*(4.0 + alpha(mgs,lh))/
!!     :         ((3.0 + alpha(mgs,lh))*(2.0 + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))
!          alp = Max( 1.0, alpha(mgs,lh)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
         endif
        
        else
         qhacw(mgs) =    &
     &   min(   &
     &   ((0.25)*pi)*ehw(mgs)*(qx(mgs,lc)-qcwresv(mgs))*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf3*xdia(mgs,lh,2)    &
     &    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,lc,1)    &
     &    + gf1*xdia(mgs,lc,2) )     &
     &    , 0.5*(qx(mgs,lc)-qcwresv(mgs))*dtpinv)
!     <  , qxmxd(mgs,lc))
!     <  , qcmxd(mgs))
       
       
         if ( lwsm6 .and. qsacw(mgs) > 0.0 .and.  qhacw(mgs) > 0.0) then
           qaacw = ( qx(mgs,ls)*qsacw(mgs) + qx(mgs,lh)*qhacw(mgs) )/(qx(mgs,ls) + qx(mgs,lh))
!           qaacw = Min( qaacw, 0.5*(qsacw(mgs) + qhacw(mgs) ) )
           qsacw(mgs) = qaacw
           qhacw(mgs) = qaacw
         endif
         
       endif

          qhacwmlr(mgs) = qhacw(mgs)
          if ( temg(mgs) > tfr .and. iqhacwshr == 0 ) then
            qhacw(mgs) = 0.0
          endif
          
          if ( lvol(lh) .gt. 1 .or. lhl .gt. 1 ) then ! calculate rime density for graupel volume and/or for graupel conversion to hail
             
             if ( temg(mgs) .lt. 273.15) then
               if ( irimdenopt == 1 ) then ! Heymsfield and Pflaum (1985)
               vt = ( (1.0-rimdenvwgt)*vtxbar(mgs,lh,1) + rimdenvwgt*vtxbar(mgs,lh,2) )
               
             rimdn(mgs,lh) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vt )   &
     &                /(temg(mgs)-273.15))**(rimc2)
!             rimdn(mgs,lh) = Min( Max( hdnmn, rimc3, rimdn(mgs,lh) ), rimc4 )
             rimdn(mgs,lh) = Min( Max( rimc3, rimdn(mgs,lh) ), rimc4 )

!               IF ( igs(mgs) == 30 ) THEN
!                 write(0,*) 'k,vt: ',kgs(mgs),vt, vtxbar(mgs,lh,1),vtxbar(mgs,lh,2), rhovt(mgs)*axx(mgs,lh)*( (alpha(mgs,lh)+3.)*xdia(mgs,lh,1) )**bxx(mgs,lh)
!                 write(0,*) 'diam: char, mean, maxmass = ',xdia(mgs,lh,1),xdia(mgs,lh,3),(alpha(mgs,lh)+3.)*xdia(mgs,lh,1)
!                 write(0,*) 'ax,bx,cd,xdn = ',axx(mgs,lh),bxx(mgs,lh),cdxgs(mgs,lh),xdn(mgs,lh)
!                 write(0,*) 'vt_char,vt_mean = ',rhovt(mgs)*axx(mgs,lh)*( xdia(mgs,lh,1) )**bxx(mgs,lh),rhovt(mgs)*axx(mgs,lh)*( xdia(mgs,lh,3) )**bxx(mgs,lh)
!                 write(0,*) 'rimdn,alpha = ',rimdn(mgs,lh),alpha(mgs,lh)
!               ENDIF

               elseif ( irimdenopt == 2 ) then ! Cober and List (1993)

                tmp = (-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lh,1) + rimdenvwgt*vtxbar(mgs,lh,2) )   &
     &                /(temg(mgs)-273.15))
                tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) ) ! have to limit range of "R" because quadratic function starts to decrease (unphysically) at higher values
                
                rimdn(mgs,lh) = 1000.*(0.051 + 0.114*tmp - 0.0055*tmp**2)

               elseif ( irimdenopt == 3 .or. irimdenopt == 4) then ! Macklin (3) or Saunders and Hosseini 2001

                tmp = (-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lh,1) + rimdenvwgt*vtxbar(mgs,lh,2) )   &
     &                /(temg(mgs)-273.15))
              !  tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                if ( irimdenopt == 3 ) then
                  rimdn(mgs,lh) =  Min(900., Max( 170., 110.*tmp**0.76 ) )
                elseif ( irimdenopt == 4 ) then ! Saunders and Hosseini
                  rimdn(mgs,lh) =  Min(917., Max( 10.,  900.0*(1.0 - 0.905**tmp ) ) )
                endif
               
               endif
             else
             rimdn(mgs,lh) = 1000.
             endif
             
             if ( lvol(lh) > 1 ) vhacw(mgs) = rho0(mgs)*qhacw(mgs)/rimdn(mgs,lh)

          endif
      
        if ( qx(mgs,lh) .gt. qxmin(lh) .and. ipelec .ge. 1 ) then
         rarx(mgs,lh) =     &
     &    qhacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lh,2)*cx(mgs,lh))
        endif
      
      endif  
      end do   
!
!
      do mgs = 1,ngscnt
      qhaci(mgs) = 0.0
      qhaci0(mgs) = 0.0
      if ( ehi(mgs) .gt. 0.0 ) then
       if (  ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,li,1) )

          qhaci0(mgs) = 0.25*pi*ehiclsn(mgs)*cx(mgs,lh)*qx(mgs,li)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,li)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            da1(li)*xdia(mgs,li,3)**2 ) 
          qhaci(mgs) = Min( ehi(mgs)*qhaci0(mgs), qimxd(mgs) )
       else
        qhaci(mgs) =    &
     &  min(   &
     &  ((0.25)*pi)*ehi(mgs)*ehiclsn(mgs)*qx(mgs,li)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))   &
     &  *(  gf3*xdia(mgs,lh,2)    &
     &    + 2.0*gf2*xdia(mgs,lh,1)*xdia(mgs,li,1)    &
     &    + gf1*xdia(mgs,li,2) )     &
     &  , qimxd(mgs))
       endif
      endif
      end do   


      if ( lis > 1 .and. ipconc >= 5 ) then
      do mgs = 1,ngscnt
      qhacis(mgs) = 0.0
      qhacis0(mgs) = 0.0
      if ( ehis(mgs) .gt. 0.0 ) then

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,lis,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,lis,1) )

          qhacis0(mgs) = 0.25*pi*ehisclsn(mgs)*cx(mgs,lh)*qx(mgs,lis)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,lis)*xdia(mgs,lh,3)*xdia(mgs,lis,3) +    &
     &            da1(li)*xdia(mgs,lis,3)**2 ) 
          qhacis(mgs) = Min( ehis(mgs)*qhacis0(mgs), qxmxd(mgs,lis) )
      endif
      end do
      endif

!
!
      do mgs = 1,ngscnt
      qhacs(mgs) = 0.0
      qhacs0(mgs) = 0.0
      if ( ehs(mgs) .gt. 0.0 ) then
       if ( ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,ls,1) )

          qhacs0(mgs) = 0.25*pi*ehsclsn(mgs)*cx(mgs,lh)*qx(mgs,ls)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,ls)*xdia(mgs,lh,3)*xdia(mgs,ls,3) +    &
     &            da1(ls)*xdia(mgs,ls,3)**2 ) 
      
          qhacs(mgs) = Min( ehs(mgs)*qhacs0(mgs), qsmxd(mgs) )

       else
         qhacs(mgs) =   &
     &   min(   &
     &   ((0.25)*pi/gf4)*ehs(mgs)*ehsclsn(mgs)*qx(mgs,ls)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))   &
     &  *(  gf6*gf1*xdia(mgs,ls,2)   &
     &    + 2.0*gf5*gf2*xdia(mgs,ls,1)*xdia(mgs,lh,1)   &
     &    + gf4*gf3*xdia(mgs,lh,2) )   &
     &  , qsmxd(mgs))
        endif
      endif
      end do   
!
      do mgs = 1,ngscnt
      qhacr(mgs) = 0.0
      qhacrmlr(mgs) = 0.0
      vhacr(mgs) = 0.0
      chacr(mgs) = 0.0
      zhacr(mgs) = 0.0
      if ( temg(mgs) .gt. tfr ) raindn(mgs,lh) = 1000.0

      if ( ehr(mgs) .gt. 0.0 ) then
      if ( ipconc .ge. 3 ) then
       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,lr,1) )
!       qhacr(mgs) = ehr(mgs)*cx(mgs,lh)*vt*
!     :     qx(mgs,lr)*0.25*pi*
!     :      (3.02787*xdia(mgs,lr,2) + 
!     :       3.30669*xdia(mgs,lh,1)*xdia(mgs,lr,1) + 
!     :       2.*xdia(mgs,lh,2))
     
       qhacr(mgs) = 0.25*pi*ehr(mgs)*cx(mgs,lh)*qx(mgs,lr)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab1lh(mgs,lh,lr)*xdia(mgs,lh,3)*xdia(mgs,lr,3) +    &
     &            da1lr(mgs)*xdia(mgs,lr,3)**2 )
!     &            da1(lr)*xdia(mgs,lr,3)**2 )
!       IF ( qhacr(mgs) .gt. 0. .or. tmp .gt. 0.0 ) write(0,*) 'qhacr= ',qhacr(mgs),tmp
!!        qhacr(mgs) = Min( qhacr(mgs), qrmxd(mgs) )
!!        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
!!        chacr(mgs) = min(chacr(mgs),crmxd(mgs))

        qhacr(mgs) = Min( qhacr(mgs), qxmxd(mgs,lr) )

            qhacrmlr(mgs) = qhacr(mgs)
        
        if ( temg(mgs) > tfr .and. iehr0c == 0 ) then
          qhacr(mgs) = 0.0

          if ( iqhacrmlr == 0 ) then
              qhacrmlr(mgs) = -qhacw(mgs)
          endif

        else
!        chacr(mgs) = Min( qhacr(mgs)*rho0(mgs)/xmas(mgs,lr), cxmxd(mgs,lr) )

!       chacr(mgs) = ehr(mgs)*cx(mgs,lh)*vt*
!     :     cx(mgs,lr)*0.25*pi*
!     :      (0.69874*xdia(mgs,lr,2) +
!     :       1.24001*xdia(mgs,lh,1)*xdia(mgs,lr,1) +
!     :       2.*xdia(mgs,lh,2))

        chacr(mgs) = 0.25*pi*ehr(mgs)*cx(mgs,lh)*cx(mgs,lr)*vt*      &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +                     &
     &            dab0lh(mgs,lh,lr)*xdia(mgs,lh,3)*xdia(mgs,lr,3) +  &
     &            da0lr(mgs)*xdia(mgs,lr,3)**2 )

!       IF ( qhacr(mgs) .gt. 0. .or. tmp .gt. 0.0 ) write(0,*) 'chacr= ',chacr(mgs),tmp

!        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
        chacr(mgs) = min(chacr(mgs),crmxd(mgs))

      if ( lzh .gt. 1 ) then
          tmp = qx(mgs,lh)/cx(mgs,lh)

!          g1 = (6.0 + alpha(mgs,lh))*(5.0 + alpha(mgs,lh))*(4.0 + alpha(mgs,lh))/
!     :         ((3.0 + alpha(mgs,lh))*(2.0 + alpha(mgs,lh))*(1.0 + alpha(mgs,lh)))
!          alp = Max( 1.0, alpha(mgs,lh)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) - tmp**2 * chacr(mgs) )
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhacr(mgs) )
      endif
      endif ! temg > tfr
      
      else
       if ( lwsm6 .and. ipconc == 0 ) then
         vt = vt2ave(mgs)
       else
         vt = vtxbar(mgs,lh,1)
       endif

      qhacr(mgs) =   &
     &   min(   &
     &   ((0.25)*pi/gf4)*ehr(mgs)*qx(mgs,lr)*cx(mgs,lh)   &
     &  *abs(vt-vtxbar(mgs,lr,1))   &
     &  *(  gf6*gf1*xdia(mgs,lr,2)   &
     &    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,lh,1)   &
     &    + gf4*gf3*xdia(mgs,lh,2) )   &
     &  , qrmxd(mgs))
      
        if ( temg(mgs) > tfr ) then
          if ( iqhacrmlr >= 1 ) qhacrmlr(mgs) = qhacr(mgs)
          qhacr(mgs) = 0.0
        endif
      
      endif
          if ( lvol(lh) .gt. 1 .or. lhl .gt. 1 ) then ! calculate rime density for graupel volume and/or for graupel conversion to hail
             
             if ( temg(mgs) .lt. 273.15) then
             raindn(mgs,lh) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lr,3))   &
     &                *((0.60)*vt)   &
     &                /(temg(mgs)-273.15))**(rimc2)

             raindn(mgs,lh) = Min( Max( rimc3, rimdn(mgs,lh) ), rimc4 )
             else
             raindn(mgs,lh) = 1000.
             endif
             
             if ( lvol(lh) > 1 )  vhacr(mgs) = rho0(mgs)*qhacr(mgs)/raindn(mgs,lh)
        endif
      endif
      end do

!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: hail collects xxxxx'
!

      do mgs = 1,ngscnt
      qhlacw(mgs) = 0.0
      qhlacwmlr(mgs) = 0.0
      vhlacw(mgs) = 0.0
      vhlsoak(mgs) = 0.0
      if ( lhl > 1 .and. .true.) then
        vtmax = (gz(igs(mgs),jgs,kgs(mgs))*dtpinv)
        vtxbar(mgs,lhl,1) = Min( vtmax, vtxbar(mgs,lhl,1))
        vtxbar(mgs,lhl,2) = Min( vtmax, vtxbar(mgs,lhl,2))
        vtxbar(mgs,lhl,3) = Min( vtmax, vtxbar(mgs,lhl,3))
      endif

      if ( lhl > 0 ) then
      rarx(mgs,lhl) = 0.0
      endif

      if ( lhl .gt. 1 .and. ehlw(mgs) .gt. 0.0 ) then


!        IF ( ipconc .ge. 2 ) THEN

            vt = abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))

          qhlacw(mgs) = 0.25*pi*ehlw(mgs)*cx(mgs,lhl)*(qx(mgs,lc)-qcwresv(mgs))*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,lc)*xdia(mgs,lhl,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 )


          qhlacw(mgs) = Min( qhlacw(mgs), 0.5*qx(mgs,lc)*dtpinv )

          qhlacwmlr(mgs) = qhlacw(mgs)
          if ( temg(mgs) > tfr .and. iqhlacwshr == 0 ) then
            qhlacw(mgs) = 0.0
          endif

          if ( lvol(lhl) .gt. 1 ) then

             if ( temg(mgs) .lt. 273.15) then
               if ( irimdenopt == 1 ) then ! Heymsfeld and Pflaum (1985)
             rimdn(mgs,lhl) = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*( (1.0-rimdenvwgt)*vtxbar(mgs,lhl,1) + rimdenvwgt*vtxbar(mgs,lhl,2) ))   &
     &                /(temg(mgs)-273.15))**(rimc2)
             rimdn(mgs,lhl) = Min( Max( hldnmn, rimc3, rimdn(mgs,lhl) ), rimc4 )
               
               elseif ( irimdenopt == 2 ) then ! Cober and List (1993)
                tmp = -0.5*(1.e+06)*xdia(mgs,lc,1)   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lhl,1) + rimdenvwgt*vtxbar(mgs,lhl,2) )   &
     &                /(temg(mgs)-273.15)
                tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                rimdn(mgs,lhl) = 1000.*(0.051 + 0.114*tmp - 0.005*tmp**2)
               
               elseif ( irimdenopt == 3 .or. irimdenopt == 4) then ! Macklin (3) or Saunders and Hosseini 2001
                tmp = -0.5*(1.e+06)*xdia(mgs,lc,1)   &
     &                *( (1.0-rimdenvwgt)*vtxbar(mgs,lhl,1) + rimdenvwgt*vtxbar(mgs,lhl,2) )  &
     &                /(temg(mgs)-273.15)
              !  tmp = Min( 5.5/0.6, Max( 0.3/0.6, tmp ) )
                
                if ( irimdenopt == 3 ) then
                  rimdn(mgs,lhl) =  Min(900., Max( 170., 110.*tmp**0.76 ) )
                elseif ( irimdenopt == 4 ) then ! Saunders and Hosseini
                  rimdn(mgs,lhl) =  Min(917., Max( 10.,  900.0*(1.0 - 0.905**tmp ) ) )
                endif
               
               endif
             else
             rimdn(mgs,lhl) = 1000.
             endif

             vhlacw(mgs) = rho0(mgs)*qhlacw(mgs)/rimdn(mgs,lhl)

          endif


        if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. ipelec .ge. 1 ) then
         rarx(mgs,lhl) =     &
     &    qhlacw(mgs)*1.0e3*rho0(mgs)/((pi/2.0)*xdia(mgs,lhl,2)*cx(mgs,lhl))
        endif

      endif
      end do

      qhlaci(:) = 0.0
      qhlaci0(:) = 0.0
      if ( lhl .gt. 1  ) then
      do mgs = 1,ngscnt
      if ( ehli(mgs) .gt. 0.0 ) then
       if (  ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,li,1) )

          qhlaci0(mgs) = 0.25*pi*ehliclsn(mgs)*cx(mgs,lhl)*qx(mgs,li)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,li)*xdia(mgs,lhl,3)*xdia(mgs,li,3) +    &
     &            da1(li)*xdia(mgs,li,3)**2 )
        ! qhlaci(mgs) = Min( qhlaci(mgs), qimxd(mgs) )
          qhlaci(mgs) = Min( ehli(mgs)*qhlaci0(mgs), qimxd(mgs) )
       endif
      endif
      end do
      endif
!
      qhlacs(:) = 0.0
      qhlacs0(:) = 0.0
      if ( lhl .gt. 1 ) then
      do mgs = 1,ngscnt
      if ( ehls(mgs) .gt. 0.0) then
       if ( ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,ls,1) )

          qhlacs0(mgs) = 0.25*pi*ehlsclsn(mgs)*cx(mgs,lhl)*qx(mgs,ls)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,ls)*xdia(mgs,lhl,3)*xdia(mgs,ls,3) +    &
     &            da1(ls)*xdia(mgs,ls,3)**2 )

          qhlacs(mgs) = Min( ehls(mgs)*qhlacs0(mgs), qsmxd(mgs) )
        endif
      endif
      end do
      endif


      do mgs = 1,ngscnt
      qhlacr(mgs) = 0.0
      qhlacrmlr(mgs) = 0.0
      chlacr(mgs) = 0.0
      vhlacr(mgs) = 0.0
      if ( lhl .gt. 1 .and. temg(mgs) .gt. tfr ) raindn(mgs,lhl) = 1000.0

      if ( lhl .gt. 1 .and. ehlr(mgs) .gt. 0.0 ) then
      if ( ipconc .ge. 3 ) then
       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,lr,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,lr,1) )

       qhlacr(mgs) = 0.25*pi*ehlr(mgs)*cx(mgs,lhl)*qx(mgs,lr)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab1lh(mgs,lhl,lr)*xdia(mgs,lhl,3)*xdia(mgs,lr,3) +    &
     &            da1lr(mgs)*xdia(mgs,lr,3)**2 )
!     &            da1(lr)*xdia(mgs,lr,3)**2 )
!       IF ( qhacr(mgs) .gt. 0. .or. tmp .gt. 0.0 ) write(0,*) 'qhacr= ',qhacr(mgs),tmp
!!        qhacr(mgs) = Min( qhacr(mgs), qrmxd(mgs) )
!!        chacr(mgs) = qhacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
!!        chacr(mgs) = min(chacr(mgs),crmxd(mgs))

        qhlacr(mgs) = Min( qhlacr(mgs), qxmxd(mgs,lr) )

     
        if ( iqhlacrmlr >= 1 ) qhlacrmlr(mgs) = qhlacr(mgs)
        
        if ( temg(mgs) > tfr .and. iehlr0c == 0) then
          qhlacr(mgs) = 0.0
          if ( iqhlacrmlr == 0 ) then
              qhlacrmlr(mgs) = -qhlacw(mgs)
          endif
        else
        chlacr(mgs) = 0.25*pi*ehlr(mgs)*cx(mgs,lhl)*cx(mgs,lr)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0lh(mgs,lhl,lr)*xdia(mgs,lhl,3)*xdia(mgs,lr,3) +    &
     &            da0lr(mgs)*xdia(mgs,lr,3)**2 )

        chlacr(mgs) = min(chlacr(mgs),crmxd(mgs))

        if ( lvol(lhl) .gt. 1 ) then
         vhlacr(mgs) = rho0(mgs)*qhlacr(mgs)/raindn(mgs,lhl)
        endif
        endif
      endif
      endif
      end do



!
!
!
!
!      if (ndebug .gt. 0 ) write(0,*) 'Collection: Cloud collects xxxxx'

      if (ndebug .gt. 0 ) write(0,*) 'Collection: cloud ice collects xxxx2'
!
      do mgs = 1,ngscnt
      qiacw(mgs) = 0.0
      if ( eiw(mgs) .gt. 0.0 ) then

       vt = Sqrt((vtxbar(mgs,li,1)-vtxbar(mgs,lc,1))**2 +    &
     &            0.04*vtxbar(mgs,li,1)*vtxbar(mgs,lc,1) )

          qiacw(mgs) = 0.25*pi*eiw(mgs)*cx(mgs,li)*qx(mgs,lc)*vt*   &
     &         (  da0(li)*xdia(mgs,li,3)**2 +     &
     &            dab1(li,lc)*xdia(mgs,li,3)*xdia(mgs,lc,3) +    &
     &            da1lc(mgs)*xdia(mgs,lc,3)**2 )

       qiacw(mgs) = Min( qiacw(mgs), qxmxd(mgs,lc) )
      endif
      end do


!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: cloud ice collects xxxx8'
!
      do mgs = 1,ngscnt
      qiacr(mgs) = 0.0
      qiacrf(mgs) = 0.0
      qiacrs(mgs) = 0.0
      ciacrs(mgs) = 0.0
      ciacr(mgs) = 0.0
      ciacrf(mgs) = 0.0
      viacrf(mgs) = 0.0
      csplinter(mgs) = 0.0
      qsplinter(mgs) = 0.0
      csplinter2(mgs) = 0.0
      qsplinter2(mgs) = 0.0
      if ( iacr .ge. 1 .and. eri(mgs) .gt. 0.0    &
     &     .and. temg(mgs) .le. 270.15 ) then
      if ( ipconc .ge. 3 ) then
       ni = 0.0
         if ( xdia(mgs,li,1) .ge. 10.e-6 ) then
          ni = ni + cx(mgs,li)*Exp(- (40.e-6/xdia(mgs,li,1))**3 )
         endif
       if ( imurain == 1 ) then ! gamma of diameter
           if ( iacrsize /= 4 ) then
           if ( iacrsize .eq. 1 ) then
             ratio = 500.e-6/xdia(mgs,lr,1)
           elseif ( iacrsize .eq. 2 ) then
             ratio = 300.e-6/xdia(mgs,lr,1)
           elseif ( iacrsize .eq. 3 ) then
             ratio = 40.e-6/xdia(mgs,lr,1)
           elseif ( iacrsize .eq. 5 ) then
             ratio = 150.e-6/xdia(mgs,lr,1)
           endif
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
!           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lr)))*dqiacralphainv)
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lr) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio
           tmp1 = ciacrratio(i,j) + delx*dqiacrratioinv*(ciacrratio(ip1,j) - ciacrratio(i,j))
           tmp2 = ciacrratio(i,jp1) + delx*dqiacrratioinv*(ciacrratio(ip1,jp1) - ciacrratio(i,jp1))
           
           ! interpolate along alpha
           
           nr = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*cx(mgs,lr)
           
           ! interpolate along x, i.e., ratio; 
           tmp1 = qiacrratio(i,j) + delx*dqiacrratioinv*(qiacrratio(ip1,j) - qiacrratio(i,j))
           tmp2 = qiacrratio(i,jp1) + delx*dqiacrratioinv*(qiacrratio(ip1,jp1) - qiacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           qr = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*qx(mgs,lr)
           
           else ! iacrsize == 4 : use all
             nr = cx(mgs,lr)
             qr = qx(mgs,lr)
           endif

          vt = Sqrt((vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))**2 +     &
     &            0.04*vtxbar(mgs,lr,1)*vtxbar(mgs,li,1) )

          qiacr(mgs) = 0.25*pi*eri(mgs)*ni*qr*vt*   &
     &         (  da0(li)*xdia(mgs,li,3)**2 +     &
     &            dab1lh(mgs,li,lr)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            da1(lr)*xdia(mgs,lr,3)**2 ) 
          
          qiacr(mgs) = Min( qrmxd(mgs), qiacr(mgs) )
          

          ciacr(mgs) = 0.25*pi*eri(mgs)*ni*nr*vt*   &
     &         (  da0(li)*xdia(mgs,li,3)**2 +     &
     &            dab0lh(mgs,li,lr)*xdia(mgs,lr,3)*xdia(mgs,li,3) +    &
     &            da0(lr)*xdia(mgs,lr,3)**2 ) 

          ciacr(mgs) = Min( crmxd(mgs), ciacr(mgs) )
          
!          write(iunit,*) 'qiacr: ',cx(mgs,lr),nr,qx(mgs,lr),qr,qiacr(mgs),ciacr(mgs)
!          write(iunit,*) 'xdia r li = ',xdia(mgs,lr,3),xdia(mgs,li,3),xdia(mgs,lr,1),xdia(mgs,li,1)
!          write(iunit,*) 'i,j,ratio = ',i,j,ciacrratio(i,j),qiacrratio(i,j)
!          write(iunit,*) 'ni,ci = ',ni,cx(mgs,li),qx(mgs,li)

       elseif ( imurain == 3 ) then ! gamma of volume
!   Set nr to the number of drops greater than 40 microns.
         arg = 1000.*xdia(mgs,lr,3)
!         nr = cx(mgs,lr)*gaml02( arg )
!        IF ( iacr .eq. 1 ) THEN
         if ( ipconc .ge. 3 ) then
           if ( iacrsize .eq. 1 ) then
            nr = cx(mgs,lr)*gaml02d500( arg )  ! number greater than 500 microns in diameter
           elseif ( iacrsize .eq. 2 .or. iacrsize .eq. 5 ) then
            nr = cx(mgs,lr)*gaml02d300( arg )  ! number greater than 300 microns in diameter
           elseif ( iacrsize .eq. 3 ) then
            nr = cx(mgs,lr)*gaml02( arg ) ! number greater than 40 microns in diameter
           elseif ( iacrsize .eq. 4 ) then
            nr = cx(mgs,lr) ! all raindrops
           endif
         else
         nr = cx(mgs,lr)*gaml02( arg )
         endif
!        ELSEIF ( iacr .eq. 2 ) THEN
!         nr = cx(mgs,lr)*gaml02d300( arg )  ! number greater than 300 microns in diameter
!        ENDIF
       if ( ni .gt. 0.0 .and. nr .gt. 0.0 ) then
       d0 = xdia(mgs,lr,3)
       qiacr(mgs) = xdn(mgs,lr)*rhoinv(mgs)*   &
     &     (0.217239*(0.522295*(d0**5) +    &
     &      49711.81*(d0**6) -    &
     &      1.673016e7*(d0**7)+    &
     &      2.404471e9*(d0**8) -    &
     &      1.22872e11*(d0**9))*ni*nr)
      qiacr(mgs) = Min( qrmxd(mgs), qiacr(mgs) )
      ciacr(mgs) =   &
     &   (0.217239*(0.2301947*(d0**2) +    &
     &      15823.76*(d0**3) -    &
     &      4.167685e6*(d0**4) +    &
     &      4.920215e8*(d0**5) -    &
     &      2.133344e10*(d0**6))*ni*nr)
      ciacr(mgs) = Min( crmxd(mgs), ciacr(mgs) )
!      ciacr(mgs) = qiacr(mgs)*cx(mgs,lr)/qx(mgs,lr)
      endif
      endif
       if ( iacr .eq. 1 .or. iacr .eq. 3 ) then
         ciacrf(mgs) = Min(ciacr(mgs), qiacr(mgs)/(1.0*vr1mm*1000.0)*rho0(mgs) ) ! *rzxh(mgs)
       elseif ( iacr .eq. 2 ) then
         ciacrf(mgs) = ciacr(mgs) ! *rzxh(mgs)
       elseif ( iacr .eq. 4 ) then
         ciacrf(mgs) = Min(ciacr(mgs), qiacr(mgs)/(1.0*vfrz*1000.0)*rho0(mgs) ) ! *rzxh(mgs)
       elseif ( iacr .eq. 5 ) then
         ciacrf(mgs) = ciacr(mgs)*rzxh(mgs)
       endif 
!      crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*27.0*vr1mm*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
       endif
      
      
      else ! single-moment rain
      qiacr(mgs) =    &
     &  min(        &
     &   ((0.25/gf4)*pi)*eri(mgs)*cx(mgs,li)*qx(mgs,lr)   &
     &  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,li,1))   &
     &  *(  gf6*gf1*xdia(mgs,lr,2)    &
     &    + 2.0*gf5*gf2*xdia(mgs,lr,1)*xdia(mgs,li,1)    &
     &    + gf4*gf3*xdia(mgs,li,2) )     &
     &  , qrmxd(mgs))
      endif
!      if ( temg(mgs) .gt. 268.15 ) then
!      qiacr(mgs) = 0.0
!      ciacr(mgs) = 0.0
!      end if

      if ( ipconc .ge. 1 ) then
        if ( nsplinter .ge. 1000 ) then
        ! Lawson et al. 2015 JAS
         ! ave. diam of freezing drops in microns
           if ( qiacr(mgs)*dtp > qxmin(lh) .and. ciacr(mgs) > 1.e-3 ) then
             tmpdiam = 1.e6*( 6.*qiacr(mgs)/(1000.*pi*ciacr(mgs) ) )**(1./3.) ! avg. diameter of newly frozen drops in microns
             csplinter(mgs) = lawson_splinter_fac*tmpdiam**4*ciacr(mgs)
           endif
        elseif ( nsplinter .ge. 0 ) then
          csplinter(mgs) = nsplinter*ciacr(mgs)
        else
          csplinter(mgs) = -nsplinter*ciacrf(mgs)
        endif
        qsplinter(mgs) = Min(0.1*qiacr(mgs), csplinter(mgs)*splintermass/rho0(mgs) ) ! makes splinters smaller if too much mass is taken from graupel
      endif
      
      frach = 1.0
           if ( ibiggsnow == 2 .or. ibiggsnow == 3 ) then
           if ( ciacr(mgs) > qxmin(lh) ) then
           xvfrz = rho0(mgs)*qiacr(mgs)/(ciacr(mgs)*900.) ! mean volume of frozen drops; 900. for frozen drop density
           frach = 0.5 *(1. +  Tanh(0.2e12 *( xvfrz - 1.15*xvbiggsnow)))

             qiacrs(mgs) = (1.-frach)*qiacr(mgs)
             ciacrs(mgs) = (1.-frach)*ciacrf(mgs) ! *rzxh(mgs)
           
           endif
           endif

      qiacrf(mgs) = frach*qiacr(mgs)
      ciacrf(mgs) = frach*ciacrf(mgs)

      if ( lvol(lh) > 1 ) then
         viacrf(mgs) = rho0(mgs)*qiacrf(mgs)/rhofrz
      endif
      
      end do
!
!
!
!

! snow aggregation here
      if ( ipconc .ge. 4 ) then !
      do mgs = 1,ngscnt
      csacs(mgs) = 0.0
      if ( qx(mgs,ls) > qxmin(ls) .and. ess(mgs) .gt. 0.0 ) then ! .and. xv(mgs,ls) < 0.25*xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls)))  ) then

        if ( iessec0flag == 0 ) then
          ec0(mgs) = 1.0
        else
          tmp = xv(mgs,ls)/(xvmx(ls)*Max(1.,100./Min(100.,xdn(mgs,ls)))) ! fraction of max snow mass
          if ( tmp .lt. essfrac1 ) then
            ec0(mgs) = 1.0
          elseif ( tmp .ge. essfrac2 ) then
            ec0(mgs) = 0.0
          else
            ec0(mgs) = (essfrac2 - tmp)/(essfrac2 - essfrac1)
          endif
        endif

      csacs(mgs) = ec0(mgs)*rvt*aa2*ess(mgs)*cx(mgs,ls)**2*Min( xv(mgs,ls), 4.*pii/3.*essrmax**3 ) ! *Min(1.,xdn(mgs,ls)/100. ) ! Min func tries to recalibrate for low diagnosed density 
!      csacs(mgs) = rvt*aa2*ess(mgs)*cx(mgs,ls)**2*Min( xv(mgs,ls), 4.*pii/3.*0.02**3 ) ! *Min(1.,xdn(mgs,ls)/100. ) ! Min func tries to recalibrate for low diagnosed density 
      csacs(mgs) = Min(csacs(mgs),csmxd(mgs))
      endif
      end do
      end if
!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 11'
      if ( ipconc .ge. 2 .or. ipelec .ge. 9 ) then
      do mgs = 1,ngscnt
      ciacw(mgs) = 0.0
      if ( eiw(mgs) .gt. 0.0 ) then
        ciacw(mgs) = qiacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        ciacw(mgs) = min(ciacw(mgs),ccmxd(mgs))
      endif
      end do

      end if

      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 18'
      if ( ipconc .ge. 2 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
       cracw(mgs) = 0.0
       cracr(mgs) = 0.0
       ec0(mgs) = 1.e9
      if ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr)    &
     &      .and. qracw(mgs) .gt. 0.0 ) then

       if ( ipconc .lt. 3 ) then
        if ( erw(mgs) .gt. 0.0 ) then
        cracw(mgs) =   &
     &   ((0.25)*pi)*erw(mgs)*(cx(mgs,lc) - ccwresv(mgs))*cx(mgs,lr)   &
     &  *abs(vtxbar(mgs,lr,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf1*xdia(mgs,lc,2)   &
     &    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lr,1)   &
     &    + gf3*xdia(mgs,lr,2) )
        endif
       else ! if ( ipconc .ge. 3 .and. 
        if ( dmrauto <= 0 .or.  rho0(mgs)*qx(mgs,lr) > 1.2*xl2p(mgs) ) then  !{
        if ( 0.5*xdia(mgs,lr,3) .gt. rh(mgs) ) then ! { .or. cx(mgs,lr) .gt. nh(mgs) 
!        IF ( qx(mgs,lc) .gt. qxmin(lc) .and. qx(mgs,lr) .gt. qxmin(lr) ) THEN
          if ( 0.5*xdia(mgs,lr,3) .gt. rwradmn ) then ! r > 50.e-6 
!          DM0CCC=A2*XNC*XNR*(XVC+XVR)                               ! (A11)
!         NOTE: murain drops out, so same result for imurain = 1 and 3
            cracw(mgs) = aa2*cx(mgs,lr)*(cx(mgs,lc) - ccwresv(mgs))*(xv(mgs,lc) + xv(mgs,lr))
          else
            if ( imurain == 3 ) then
!          DM0CCC=A1*XNC*XNR*(((CNU+2.)/(CNU+1.))*XVC**2+((RNU+2.)/(RNU+1.))*XVR**2) ! (A13)
            cracw(mgs) = aa1*cx(mgs,lr)*(cx(mgs,lc) - ccwresv(mgs))*   &
     &          ((alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.) +    &
     &          (alpha(mgs,lr) + 2.)*xv(mgs,lr)**2/(alpha(mgs,lr) + 1.))
            else ! imurain == 1 use CP00 for rain DSD in diameter
            cracw(mgs) = aa1*cx(mgs,lr)*(cx(mgs,lc) - ccwresv(mgs))*   &
     &          ((alpha(mgs,lc) + 2.)*xv(mgs,lc)**2/(alpha(mgs,lc) + 1.) +    &
     &          (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)*xv(mgs,lr)**2/  &
     &             ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.)) )
            endif ! imurain
          endif
        endif ! } rh
        endif ! } dmrauto
       endif ! ipconc
      endif ! qc > qcmin & qr > qrmin
        
! Rain self collection (cracr) and break-up (factor of ec0)
!
!       
        ec0(mgs) = 2.e9
        if ( qx(mgs,lr) .gt. qxmin(lr) ) then
        rwrad = 0.5*xdia(mgs,lr,3)
        
        
        ! check median volume diameter
        if ( icracrthresh > 1 ) then
         if ( imurain == 1 ) then
           tmp =  (3.67+alpha(mgs,lr))*xdia(mgs,lr,1) ! median volume diameter; units of mm (Ulbrich 1983, JCAM)
         else ! imurain == 3, 
           tmp =  (1.678+alpha(mgs,lr))**(1./3.)*xdia(mgs,lr,1) ! units of mm (using method of Ulbrich 1983. See ventillation_stuff.nb)
         endif
        else
          tmp = xdia(mgs,lr,3) - 0.1e-3
        endif
         
!        IF ( xdia(mgs,lr,3) .gt. 2.0e-3 .or. icracr <= 0 ) THEN
        if ( tmp .gt. 1.9e-3 .or. icracr <= 0 ) then
          ec0(mgs) = 0.0
          cracr(mgs) = 0.0
        else
         if ( dmrauto <= 0 .or.  rho0(mgs)*qx(mgs,lr) > 1.2*xl2p(mgs) ) then 
          if ( xdia(mgs,lr,3) .lt. 6.1e-4 ) then
            ec0(mgs) = 1.0
          else
            ec0(mgs) = Exp(-50.0*(50.0*(xdia(mgs,lr,3) - 6.0e-4)))
          endif
          

          if ( rwrad .ge. 50.e-6 ) then
              cracr(mgs) = ec0(mgs)*aa2*cx(mgs,lr)**2*xv(mgs,lr)
          else
            if ( imurain == 3 ) then
             cracr(mgs) = ec0(mgs)*aa1*(cx(mgs,lr)*xv(mgs,lr))**2*   &
     &                   (alpha(mgs,lr) + 2.)/(alpha(mgs,lr) + 1.)
            else ! imurain == 1
             cracr(mgs) = ec0(mgs)*aa1*(cx(mgs,lr)*xv(mgs,lr))**2*   &
     &                   (alpha(mgs,lr) + 6.)*(alpha(mgs,lr) + 5.)*(alpha(mgs,lr) + 4.)/ &
     &                  ((alpha(mgs,lr) + 3.)*(alpha(mgs,lr) + 2.)*(alpha(mgs,lr) + 1.))
              
            endif
          endif
!          cracr(mgs) = Min(cracr(mgs),crmxd(mgs))
         endif
        endif
        endif

!      cracw(mgs) = min(cracw(mgs),cxmxd(mgs,lc)) 
      end do
      end if

!
!
!
!  Graupel
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22ii'
      chacw(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt

      if ( ipconc .ge. 5 ) then
       if ( qhacw(mgs) .gt. 0.0 .and. xmas(mgs,lc) .gt. 0.0 ) then

!  This is the explict version of chacw, which turns out to be very close to the
!  approximation that the droplet size does not change, to within a few percent.
!  This may _not_ be the case for cnu other than zero!
!          chacw(mgs) = (ehw(mgs)*cx(mgs,lc)*cx(mgs,lh)*(pi/4.)*
!     :    abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))*
!     :    (2.0*xdia(mgs,lh,1)*(xdia(mgs,lh,1) +
!     :         xdia(mgs,lc,1)*gf43rds) +
!     :      xdia(mgs,lc,2)*gf53rds))

!          chacw(mgs) = Min( chacw(mgs), 0.6*cx(mgs,lc)*dtpinv )

!        chacw(mgs) = qhacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        chacw(mgs) = qhacw(mgs)*rho0(mgs)/xmascw(mgs)
!        chacw(mgs) = min(chacw(mgs),cxmxd(mgs,lc))
        chacw(mgs) = Min( chacw(mgs), 0.5*(cx(mgs,lc) - ccwresv(mgs))*dtpinv )
       else
        qhacw(mgs) = 0.0
       endif
      else
      ! single-moment
      chacw(mgs) =   &
     &   ((0.25)*pi)*ehw(mgs)*cx(mgs,lc)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,lc,1))   &
     &  *(  gf1*xdia(mgs,lc,2)   &
     &    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lh,1)   &
     &    + gf3*xdia(mgs,lh,2) )
      chacw(mgs) = min(chacw(mgs),0.5*cx(mgs,lc)*dtpinv)
!      chacw(mgs) = min(chacw(mgs),cxmxd(mgs,lc))
!      chacw(mgs) = min(chacw(mgs),ccmxd(mgs))
      endif
      end do
      end if
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22kk'
      chaci(:) = 0.0
      chaci0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      if ( ehi(mgs) .gt. 0.0 .or. ( ehiclsn(mgs) > 0.0 .and. ipelec > 0 )) then
       if ( ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,li,1) )

          chaci0(mgs) = 0.25*pi*ehiclsn(mgs)*cx(mgs,lh)*cx(mgs,li)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab0lh(mgs,lh,li)*xdia(mgs,lh,3)*xdia(mgs,li,3) +    &
     &            da0(li)*xdia(mgs,li,3)**2 )

       else
        chaci0(mgs) =   &
     &   ((0.25)*pi)*ehiclsn(mgs)*cx(mgs,li)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,li,1))   &
     &  *(  gf1*xdia(mgs,li,2)   &
     &    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lh,1)   &
     &    + gf3*xdia(mgs,lh,2) )
        endif

        chaci(mgs) = min(ehi(mgs)*chaci0(mgs),cimxd(mgs))
       endif
      end do
      end if


      chacis(:) = 0.0
      if ( lis > 1 .and. ipconc .ge. 5 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      if ( ehis(mgs) .gt. 0.0 .or. ( ehisclsn(mgs) > 0.0 .and. ipelec > 0 )) then

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,lis,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,lis,1) )

          chacis0(mgs) = 0.25*pi*ehisclsn(mgs)*cx(mgs,lh)*cx(mgs,lis)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab0lh(mgs,lh,lis)*xdia(mgs,lh,3)*xdia(mgs,lis,3) +    &
     &            da0(lis)*xdia(mgs,lis,3)**2 )


        chacis(mgs) = min(ehis(mgs)*chacis0(mgs),cxmxd(mgs,lis))
       endif
      end do
      end if
!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22nn'
      chacs(:) = 0.0
      chacs0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      if ( ehs(mgs) .gt. 0 ) then
       if ( ipconc .ge. 5 .or. ( ehsclsn(mgs) > 0.0 .and. ipelec > 0 ) ) then

       vt = Sqrt((vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lh,1)*vtxbar(mgs,ls,1) )

          chacs0(mgs) = 0.25*pi*ehsclsn(mgs)*cx(mgs,lh)*cx(mgs,ls)*vt*   &
     &         (  da0lh(mgs)*xdia(mgs,lh,3)**2 +     &
     &            dab0lh(mgs,lh,ls)*xdia(mgs,lh,3)*xdia(mgs,ls,3) +    &
     &            da0(ls)*xdia(mgs,ls,3)**2 )

       else
      chacs0(mgs) =   &
     &   ((0.25)*pi)*ehsclsn(mgs)*cx(mgs,ls)*cx(mgs,lh)   &
     &  *abs(vtxbar(mgs,lh,1)-vtxbar(mgs,ls,1))   &
     &  *(  gf3*gf1*xdia(mgs,ls,2)   &
     &    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lh,1)   &
     &    + gf1*gf3*xdia(mgs,lh,2) )
      endif
      chacs(mgs) = min(ehs(mgs)*chacs0(mgs),csmxd(mgs))
      endif
      end do
      end if


!
!
!  Hail
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22ii'
      chlacw(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt

      if ( lhl .gt. 1 .and. ipconc .ge. 5 ) then
       if ( qhlacw(mgs) .gt. 0.0 .and. xmas(mgs,lc) .gt. 0.0 ) then

!  This is the explict version of chacw, which turns out to be very close to the
!  approximation that the droplet size does not change, to within a few percent.
!  This may _not_ be the case for cnu other than zero!
!          chlacw(mgs) = (ehlw(mgs)*cx(mgs,lc)*cx(mgs,lhl)*(pi/4.)*
!     :    abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))*
!     :    (2.0*xdia(mgs,lhl,1)*(xdia(mgs,lhl,1) +
!     :         xdia(mgs,lc,1)*gf43rds) +
!     :      xdia(mgs,lc,2)*gf53rds))

!          chlacw(mgs) = Min( chlacw(mgs), 0.6*cx(mgs,lc)*dtpinv )

!        chlacw(mgs) = qhlacw(mgs)*rho0(mgs)/xmas(mgs,lc)
        chlacw(mgs) = qhlacw(mgs)*rho0(mgs)/xmascw(mgs)
!        chlacw(mgs) = min(chlacw(mgs),cxmxd(mgs,lc))
        chlacw(mgs) = Min( chlacw(mgs), 0.5*cx(mgs,lc)*dtpinv )
       else
        qhlacw(mgs) = 0.0
       endif
!      ELSE
!      chlacw(mgs) =
!     >   ((0.25)*pi)*ehlw(mgs)*cx(mgs,lc)*cx(mgs,lhl)
!     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,lc,1))
!     >  *(  gf1*xdia(mgs,lc,2)
!     >    + 2.0*gf2*xdia(mgs,lc,1)*xdia(mgs,lhl,1)
!     >    + gf3*xdia(mgs,lhl,2) )
!      chlacw(mgs) = min(chlacw(mgs),0.5*cx(mgs,lc)*dtpinv)
!      chlacw(mgs) = min(chlacw(mgs),cxmxd(mgs,lc))
!      chlacw(mgs) = min(chlacw(mgs),ccmxd(mgs))
      endif
      end do
      end if
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22kk'
      chlaci(:) = 0.0
      chlaci0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      if ( lhl .gt. 1 .and. ( ehli(mgs) .gt. 0.0 .or. (ipelec > 0 .and. ehliclsn(mgs) > 0.0) )  ) then
       if ( ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,li,1) )

          chlaci0(mgs) = 0.25*pi*ehliclsn(mgs)*cx(mgs,lhl)*cx(mgs,li)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0(lhl,li)*xdia(mgs,lhl,3)*xdia(mgs,li,3) +    &
     &            da0(li)*xdia(mgs,li,3)**2 )

!       ELSE
!        chlaci(mgs) =
!     >   ((0.25)*pi)*ehli(mgs)*cx(mgs,li)*cx(mgs,lhl)
!     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,li,1))
!     >  *(  gf1*xdia(mgs,li,2)
!     >    + 2.0*gf2*xdia(mgs,li,1)*xdia(mgs,lhl,1)
!     >    + gf3*xdia(mgs,lhl,2) )
        endif

        chlaci(mgs) = min(ehli(mgs)*chlaci0(mgs),cimxd(mgs))
       endif
      end do
      end if


      if ( lis > 1 .and. ipconc .ge. 5) then
      
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22kk'
      chlacis(:) = 0.0
      chlacis0(:) = 0.0
       do mgs = 1,ngscnt
      if ( lhl .gt. 1 .and. ( ehlis(mgs) .gt. 0.0 .or. (ipelec > 0 .and. ehlisclsn(mgs) > 0.0) )  ) then

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,lis,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,lis,1) )

          chlacis0(mgs) = 0.25*pi*ehlisclsn(mgs)*cx(mgs,lhl)*cx(mgs,lis)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0(lhl,lis)*xdia(mgs,lhl,3)*xdia(mgs,lis,3) +    &
     &            da0(lis)*xdia(mgs,lis,3)**2 )


        chlacis(mgs) = min(ehlis(mgs)*chlacis0(mgs),cxmxd(mgs,lis))
       endif
      end do
      endif

!
!
      if (ndebug .gt. 0 ) write(0,*) 'ICEZVD_GS: conc 22jj'
      chlacs(:) = 0.0
      chlacs0(:) = 0.0
      if ( ipconc .ge. 1 .or. ipelec .ge. 1 ) then
      do mgs = 1,ngscnt
      if ( lhl .gt. 1 .and. ( ehls(mgs) .gt. 0.0 .or. (ipelec > 0 .and. ehlsclsn(mgs) > 0.0) ) ) then
       if ( ipconc .ge. 5 ) then

       vt = Sqrt((vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))**2 +    &
     &            0.04*vtxbar(mgs,lhl,1)*vtxbar(mgs,ls,1) )

          chlacs0(mgs) = 0.25*pi*ehlsclsn(mgs)*cx(mgs,lhl)*cx(mgs,ls)*vt*   &
     &         (  da0lhl(mgs)*xdia(mgs,lhl,3)**2 +     &
     &            dab0(lhl,ls)*xdia(mgs,lhl,3)*xdia(mgs,ls,3) +    &
     &            da0(ls)*xdia(mgs,ls,3)**2 )

!       ELSE
!      chlacs(mgs) =
!     >   ((0.25)*pi)*ehls(mgs)*cx(mgs,ls)*cx(mgs,lhl)
!     >  *abs(vtxbar(mgs,lhl,1)-vtxbar(mgs,ls,1))
!     >  *(  gf3*gf1*xdia(mgs,ls,2)
!     >    + 2.0*gf2*gf2*xdia(mgs,ls,1)*xdia(mgs,lhl,1)
!     >    + gf1*gf3*xdia(mgs,lhl,2) )
      endif
      chlacs(mgs) = min(ehls(mgs)*chlacs0(mgs),csmxd(mgs))
      endif
      end do
      end if

!
! Ziegler (1985) autoconversion
!
!
      if ( ipconc .ge. 2 ) then
      if (ndebug .gt. 0 ) write(0,*) 'conc 26a'
      
      do mgs = 1,ngscnt
        zrcnw(mgs) = 0.0
        qrcnw(mgs) = 0.0
        crcnw(mgs) = 0.0
        cautn(mgs) = 0.0
      enddo
      
      if ( dmrauto >= -1 ) then !{
      do mgs = 1,ngscnt
!      qracw(mgs) = 0.0
!      cracw(mgs) = 0.0
       if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1000. .and. temg(mgs) .gt. tfrh+4.) then
       !( .and. w(igs(mgs),jgs,kgs(mgs)) > 5.0) THEN ! DTD: added w threshold for testing                                                                                                            
         volb = xv(mgs,lc)*(1./(1.+alpha(mgs,lc)))**(1./2.)
         cautn(mgs) = Min(ccmxd(mgs),   &
     &      ((alpha(mgs,lc)+2.)/(alpha(mgs,lc)+1.))*aa1*cx(mgs,lc)**2*xv(mgs,lc)**2)
         cautn(mgs) = Max( 0.0d0, cautn(mgs) )
         if ( rb(mgs) .le. 7.51d-6 .or. dmrauto == -1) then
           t2s = 1.d30
!           cautn(mgs) = 0.0
         else
!         XL2P=2.7E-2*XNC*XVC*((1.E12*RB**3*RC)-0.4)
         
!        T2S=3.72E-3/(((1.E4*RB)-7.5)*XNC*XVC) 
!           t2s = 3.72E-3/(((1.e6*rb)-7.5)*cx(mgs,lc)*xv(mgs,lc))
!           t2s = 3.72/(((1.e6*rb(mgs))-7.5)*rho0(mgs)*qx(mgs,lc))
           t2s = 3.72/(1.e6*(rb(mgs)-7.500d-6)*rho0(mgs)*qx(mgs,lc))

           qrcnw(mgs) = Max( 0.0d0, xl2p(mgs)/(t2s*rho0(mgs)) )
           crcnw(mgs) = Max( 0.0d0, Min(3.5e9*xl2p(mgs)/t2s,0.5*cautn(mgs)) )
           
           if ( dmrauto == 0 ) then
             if ( qx(mgs,lr)*rho0(mgs) > 1.2*xl2p(mgs) .and. cx(mgs,lr) > cxmin ) then ! Cohard and Pinty (2000a) switch over from (18) to (19)
               crcnw(mgs) = cx(mgs,lr)/qx(mgs,lr)*qrcnw(mgs)
             elseif ( ( dmropt == 1 .or. dmropt == 3 ) .and. qx(mgs,lr) > qxmin(lr) ) then
               tmp = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               crcnw(mgs) = Min(tmp,crcnw(mgs) )
             elseif ( ( dmropt == 4 ) .and. qx(mgs,lr) > qxmin(lr) ) then
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass-weighted average of old and new Dmr using converted qc mass
               crcnw(mgs) = (tmp*qrcnw(mgs)+tmp2*qx(mgs,lr))/(qrcnw(mgs)+qx(mgs,lr))
             elseif ( ( dmropt == 5 ) .and. qx(mgs,lr) > qxmin(lr) ) then
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass-weighted average of old and new Dmr using full qc mass
               crcnw(mgs) = (tmp*qx(mgs,lc)+tmp2*qx(mgs,lr))/(qx(mgs,lc)+qx(mgs,lr))
             elseif ( ( dmropt == 6 ) .and. qx(mgs,lr) > qxmin(lr) ) then
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass*diameter-weighted average of old and new Dmr (using full qc mass)
               crcnw(mgs) = (tmp*xdia(mgs,lc,3)*qx(mgs,lc)+tmp2*xdia(mgs,lr,3)*qx(mgs,lr))/(xdia(mgs,lc,3)*qx(mgs,lc)+xdia(mgs,lr,3)*qx(mgs,lr))
             elseif ( ( dmropt == 7 ) .and. qx(mgs,lr) > qxmin(lr) ) then
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try diameter-weighted average of old and new Dmr
               crcnw(mgs) = (tmp*xdia(mgs,lc,3)+tmp2*xdia(mgs,lr,3))/(xdia(mgs,lc,3)+xdia(mgs,lr,3))
             elseif ( ( dmropt == 8 ) .and. qx(mgs,lr) > qxmin(lr) ) then
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try sqrt(diameter)-weighted average of old and new Dmr
               crcnw(mgs) = (tmp*sqrt(xdia(mgs,lc,3))+tmp2*sqrt(xdia(mgs,lr,3)))/(sqrt(xdia(mgs,lc,3))+sqrt(xdia(mgs,lr,3)))
             endif
           elseif ( dmrauto == 1  .and. cx(mgs,lr) > cxmin) then
             if ( qx(mgs,lr) > qxmin(lr) ) then
               tmp = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               crcnw(mgs) = Min(tmp,crcnw(mgs) )
             endif
           elseif ( dmrauto == 2  .and. cx(mgs,lr) > cxmin) then
               tmp = crcnw(mgs)
               tmp2 = qrcnw(mgs)*cx(mgs,lr)/qx(mgs,lr)
               ! try mass-weighted average of old and new Dmr
               crcnw(mgs) = (tmp*qrcnw(mgs)+tmp2*qx(mgs,lr))/(qrcnw(mgs)+qx(mgs,lr))
           elseif ( dmrauto == 3  .and. cx(mgs,lr) > cxmin) then ! adapted from MY/CP code
              tmp = Max( 2.d0*rh(mgs), dble( xdia(mgs,lr,3) ) )
              crcnw(mgs) = rho0(mgs)*qrcnw(mgs)/(pi/6.*1000.*tmp**3)
           endif
           
           if ( crcnw(mgs) < 1.e-30 ) qrcnw(mgs) = 0.0

           if ( ipconc >= 6 ) then
           if ( lzr > 1 .and. qrcnw(mgs) > 0.0 ) then
!            vr = rho0(mgs)*qrcnw(mgs)/(1000.*crcnw(mgs))
!            zrcnw(mgs) = 36.*(xnu(lr)+2.0)*crcnw(mgs)*vr**2/((xnu(lr)+1.0)*pi**2)
             ! DTD: If rain exists at a grid point already either use the alpha-preserving Z-rate eqn. (dmrauto == 1)
             ! or a mass-weighted average of the alpha-preserving Z-rate eqn. and the init. rate eqn. (dmrauto == 2)
             ! or the original initiation rate equation (dmrauto == 0).  Not sure if this is the correct way to go but seems to work ok.
             if (qx(mgs,lr) .gt. qxmin(lr) .and. ( dmrauto == 1 .or. dmrauto ==2 ) ) then
              tmp3 = qx(mgs,lr)/cx(mgs,lr)
              tmp4 =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2* &
     &                 ( 2.*tmp3 * qrcnw(mgs) - tmp3**2 * crcnw(mgs)  )
              if (imurain == 3) then
                vr = rho0(mgs)*qrcnw(mgs)/(1000.)
                tmp3 = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
              else
                tmp3 = galpharaut*(6.*rho0(mgs)*qrcnw(mgs)/(pi*xdn0(lr)))**2/crcnw(mgs)
              endif
              if ( dmrauto == 1 ) then ! Preserve alpha
                zrcnw(mgs) = tmp4
              elseif ( dmrauto == 2 ) then ! Mass-weighted average
                zrcnw(mgs) = (tmp3*qrcnw(mgs)+tmp4*qx(mgs,lr))/(qrcnw(mgs)+qx(mgs,lr))
              endif
             else ! original formulation
              if ( imurain == 3 ) then
                vr = rho0(mgs)*qrcnw(mgs)/(1000.) ! crcnw(mgs) not divided here but is in next line, cancels one factor in the numerator
                zrcnw(mgs) = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
              else ! rain in gamma of diameter
                if ( dmropt <= 1 .or. dmropt >= 4 .or. ( qx(mgs,lr) < qxmin(lr) .and. cx(mgs,lr) < cxmin ) ) then
                  zrcnw(mgs) = galpharaut*(6.*rho0(mgs)*qrcnw(mgs)/(pi*xdn0(lr)))**2/crcnw(mgs)
                else
                  tmp3 = qx(mgs,lr)/cx(mgs,lr)
                  zrcnw(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2* &
     &                 ( 2.*tmp3 * qrcnw(mgs) - tmp3**2 * crcnw(mgs)  )
                endif
!             vr = rho0(mgs)*qrcnw(mgs)/(1000.) ! crcnw(mgs) not divided here but is in next line, cancels one factor in the numerator
!             zrcnw(mgs) = 36.*(xnu(lc)+2.0)*vr**2/(crcnw(mgs)*(xnu(lc)+1.0)*pi**2)
              endif
             endif
!             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
           endif 
           endif ! ipconc >= 6
!           IF (  crcnw(mgs) .gt. cautn(mgs) .and. crcnw(mgs) .gt. 1.0 )
!     :          THEN
!             write(0,*)  'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
!     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s,qx(mgs,lr)
!             write(0,*)  '            ',qx(mgs,lc),cx(mgs,lc),0.5e6*xdia(mgs,lc,1)
!             write(0,*)  '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
!     :         1.e6*(( 3/(4.*pi))*rho0(mgs)*qrcnw(mgs)/
!     :       (crcnw(mgs)*xdn(mgs,lr)))**(1./3.),rh(mgs)*1.e6,rwrad(mgs)
!           ELSEIF ( crcnw(mgs) .gt. 1.0 .and. cautn(mgs) .gt. 0.) THEN
!             write(0,*)  'crcnw,cautn ',crcnw(mgs)/cautn(mgs),
!     :          crcnw(mgs),cautn(mgs),igs(mgs),kgs(mgs),t2s
!             write(0,*)  '            ',rho0(mgs)*qrcnw(mgs)/crcnw(mgs),
!     :  1.e6*(( 3*pi/4.)*rho0(mgs)*qrcnw(mgs)/
!     :   (crcnw(mgs)*xdn(mgs,lr)))**(1./3.)
!           ENDIF
!           crcnw(mgs) = Min(cautn(mgs),3.5e9*xl2p(mgs)/t2s)

!           IF ( qrcnw(mgs) .gt. 0.3e-2 ) THEN
!            write(0,*)  'QRCNW'
!            write(0,*)  qrcnw(mgs),crcnw(mgs),cautn(mgs)
!            write(0,*)  xl2p,t2s,rho0(mgs),xv(mgs,lc),cx(mgs,lc),qx(mgs,lc)
!            write(0,*)  rb,0.5*xdia(mgs,lc,1),mgs,igs(mgs),kgs(mgs)
!           ENDIF
!           qrcnw(mgs) = Min(qrcnw(mgs),qcmxd(mgs))
         endif


       endif
      enddo
      
      endif !} dmrauto >= 0



      else

!
!  Berry 1968 auto conversion for rain (Orville & Kopp 1977)
!
!
      if ( ircnw .eq. 4 ) then
      do mgs = 1,ngscnt
!      sconvmix(lcw,mgs) = 0.0
      qrcnw(mgs) =  0.0
      qdiff = max((qx(mgs,lc)-qminrncw),0.0)
      if ( qdiff .gt. 0.0 .and. xdia(mgs,lc,1) .gt. 20.0e-6 ) then
      argrcnw =   &
     &  ((1.2e-4)+(1.596e-12)*(cx(mgs,lc)*1.0e-6)   &
     &  /(cwdisp*qdiff*1.0e-3*rho0(mgs)))
      qrcnw(mgs) = (rho0(mgs)*1e-3)*(qdiff**2)/argrcnw
!      sconvmix(lcw,mgs) = max(sconvmix(lcw,mgs),0.0)
      qrcnw(mgs) = (max(qrcnw(mgs),0.0))
      end if
      end do

      endif
!
!
!
!  Berry 1968 auto conversion for rain (Ferrier 1994)
!
!
      if ( ircnw .eq. 5 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      qrcnw(mgs) =  0.0
      qccrit = (pi/6.)*(cx(mgs,lc)*cwdiap**3)*xdn(mgs,lc)/rho0(mgs)
      qdiff = max((qx(mgs,lc)-qccrit),0.)
      if ( qdiff .gt. 0.0 .and. cx(mgs,lc) .gt. 1.0 ) then
      argrcnw = &
!     >  ((1.2e-4)+(1.596e-12)*cx(mgs,lc)/(cwdisp*rho0(mgs)*qdiff))   &
     &  ((1.2e-4)+(1.596e-12)*cx(mgs,lc)*1.0e-3/(cwdisp*rho0(mgs)*qdiff))
      qrcnw(mgs) = &
!     >  timflg(mgs)*rho0(mgs)*(qdiff**2)/argrcnw   &
     &  1.0e-3*rho0(mgs)*(qdiff**2)/argrcnw
      qrcnw(mgs) = Min(qxmxd(mgs,lc), (max(qrcnw(mgs),0.0)) )

!      write(iunit,*) 'qrcnw,cx =',qrcnw(mgs),cx(mgs,lc),mgs,1.e3*qx(mgs,lc),cno(lr)
      end if
      end do
      end if

!
!
!  kessler auto conversion for rain.
!
      if ( ircnw .eq. 2 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      qrcnw(mgs) = (0.001)*max((qx(mgs,lc)-qminrncw),0.0)
      end do
      end if
!
!  c4 = pi/6
!  c1 = 0.12-0.32 for colorado storms...typically 0.3-0.4
!  berry reinhart type conversion (proctor 1988)
!
      if ( ircnw .eq. 1 ) then
      do mgs = 1,ngscnt
      qrcnw(mgs) = 0.0
      c1 = 0.2
      c4 = pi/(6.0)
      bradp =    &
     & (1.e+06) * ((c1/(0.38))**(1./3.)) * (xdia(mgs,lc,1)*(0.5))
      bl2 =   &
     & (0.027) * ((100.0)*(bradp**3)*(xdia(mgs,lc,1)*(0.5)) - (0.4))
      bt2 = (bradp -7.5) / (3.72)
      qrcnw(mgs) = 0.0
      if ( bl2 .gt. 0.0 .and. bt2 .gt. 0.0 ) then
      qrcnw(mgs) = bl2 * bt2 * rho0(mgs)   &
     &  * qx(mgs,lc) * qx(mgs,lc)
      end if
      end do
      end if



      endif  !  ( ipconc .ge. 2 )

!
!
!
!  Bigg Freezing of Rain
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 27a'
      qrfrz(:) = 0.0
      qrfrzs(:) = 0.0
      qrfrzf(:) = 0.0
      vrfrzf(:) = 0.0
      crfrz(:) = 0.0
      crfrzs(:) = 0.0
      crfrzf(:) = 0.0
      zrfrz(:)  = 0.0
      zrfrzs(:)  = 0.0
      zrfrzf(:)  = 0.0
      qwcnr(:) = 0.0
      
      if ( .not. ( ipconc == 0 .and. lwsm6 ) ) then
      
      do mgs = 1,ngscnt 
      if ( qx(mgs,lr) .gt. qxmin(lr) .and. temcg(mgs) .lt. -5. .and. ibiggopt > 0 ) then
!      brz = 100.0
!      arz = 0.66
       if ( ipconc .lt. 3 ) then
       qrfrz(mgs) =    &
     &  min(   &
     &  (20.0)*(pi**2)*brz*(xdn(mgs,lr)/rho0(mgs))   &
     &   *cx(mgs,lr)*(xdia(mgs,lr,1)**6)   &
     &   *(exp(max(-arz*temcg(mgs), 0.0))-1.0)   &
     &  , qrmxd(mgs))
        qrfrzf(mgs) = qrfrz(mgs)

!       ELSEIF ( ipconc .ge. 3 .and. xv(mgs,lr) .gt. 1.1*xvmn(lr) ) THEN
       elseif ( ipconc .ge. 3 ) then
!         tmp = brz*cx(mgs,lr)*(Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
!         crfrz(mgs) = xv(mgs,lr)*tmp

         frach = 1.0d0
         
!         IF ( ibiggopt == 2 .and. imurain == 1 .and. lzr < 1 ) THEN ! lzr check because results are weird for 3-moment
         if ( ibiggopt == 2 .and. imurain == 1 ) then !
         ! integrate from Bigg diameter (for given supercooling Ts) to infinity
           
           volt = exp( 16.2 + 1.0*temcg(mgs) )* 1.0e-6 !  Ts == -temcg ; volt comes from the fit in Fig. 1 in Bigg 1953 (Proc. Phys. Soc. London) 
                                               ! for mean temperature for freezing: -ln (V) = a*Ts - b, where a = 6.9/6.8, or approx a = 1.0, and b = 16.2
                                               ! volt is given in cm**3, so convert to m**3
           dbigg = (6./pi* volt )**(1./3.) 
           
           ! perhaps should also test that W > V_t_dbigg, i.e., that drops the size of dbigg are being lifted and cooled. 
           if ( dbigg < 8.e-3 ) then !{ only bother if freezing diameter is reasonable
           
             ratio = Min(maxratiolu, dbigg/xdia(mgs,lr,1) )
           
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lr)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lr) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = ciacrratio(i,j) + delx*dqiacrratioinv*(ciacrratio(ip1,j) - ciacrratio(i,j))
           tmp2 = ciacrratio(i,jp1) + delx*dqiacrratioinv*(ciacrratio(ip1,jp1) - ciacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           crfrz(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*cx(mgs,lr)*dtpinv
           crfrzf(mgs) = crfrz(mgs)
           ! interpolate along x, i.e., ratio; 
           tmp1 = qiacrratio(i,j) + delx*dqiacrratioinv*(qiacrratio(ip1,j) - qiacrratio(i,j))
           tmp2 = qiacrratio(i,jp1) + delx*dqiacrratioinv*(qiacrratio(ip1,jp1) - qiacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           qrfrz(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*qx(mgs,lr)*dtpinv
           qrfrzf(mgs) = qrfrz(mgs)

           if ( qrfrz(mgs)*dtp < qxmin(lh) .or. crfrz(mgs)*dtp < cxmin ) then
           
             crfrz(mgs) = 0.0
             qrfrz(mgs) = 0.0
             qrfrzf(mgs) = 0.0
            
           else !{

            
           if ( ipconc >= 6 .and. lzr > 1 ) then
           ! interpolate along x, i.e., ratio; 
            tmp1 = ziacrratio(i,j) + delx*dqiacrratioinv*(ziacrratio(ip1,j) - ziacrratio(i,j))
            tmp2 = ziacrratio(i,jp1) + delx*dqiacrratioinv*(ziacrratio(ip1,jp1) - ziacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
            zrfrz(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*zx(mgs,lr)*dtpinv
           endif
           
            if ( ibiggsmallrain > 0 .and. xv(mgs,lr) < 2.*xvmn(lr) .and. ( ibiggsnow == 1 .or. ibiggsnow == 3 ) ) then
!            IF ( ibiggsmallrain > 0 .and. xv(mgs,lr) < xvbiggsnow .and. ( ibiggsnow == 1 .or. ibiggsnow == 3 ) ) THEN
             ! rain drops are so small that they cannot be pushed smaller, so put into snow (or cloud ice, depending on ifrzs)
              crfrzf(mgs) = 0.0
              qrfrzf(mgs) = 0.0
              crfrzs(mgs) = crfrz(mgs)
              qrfrzs(mgs) = qrfrz(mgs)

              if ( ipconc >= 6 .and. lzr > 1 ) then
                zrfrzs(mgs) = zrfrz(mgs)
                zrfrzf(mgs) = 0.
              endif
           elseif ( dbigg < Max( biggsnowdiam, Max(dfrz,dhmn)) .and. ( ibiggsnow == 1 .or. ibiggsnow == 3 ) ) then ! { convert some to snow or ice crystals
            ! temporarily store qrfrz and crfrz in snow terms and caclulate new crfrzf, qrfrzf, and zrfrzf. Leave crfrz etc. alone!
            
            crfrzs(mgs) = crfrz(mgs)
            qrfrzs(mgs) = qrfrz(mgs)
            
            if ( ibiggsmallrain > 0 .and. xv(mgs,lr) < 1.2*xvmn(lr) ) then
             ! rain drops are so small that they cannot be pushed smaller, so put into snow (or cloud ice, depending on ifrzs)
            crfrzf(mgs) = 0.0
            qrfrzf(mgs) = 0.0

             if (ipconc >= 6 .and. lzr > 1 ) then
               zrfrzs(mgs) = zrfrz(mgs)
               zrfrzf(mgs) = 0.
             endif
            else !{
            
           ! recalculate using dhmn for ratio
           ratio = Min( maxratiolu, Max(dfrz,dhmn)/xdia(mgs,lr,1) )
           
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
!           j = Int(Max(alphamin,Min(alphamax,alpha(mgs,lr)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lr)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lr)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lr) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = ciacrratio(i,j) + delx*dqiacrratioinv*(ciacrratio(ip1,j) - ciacrratio(i,j))
           tmp2 = ciacrratio(i,jp1) + delx*dqiacrratioinv*(ciacrratio(ip1,jp1) - ciacrratio(i,jp1))


           ! interpolate along alpha; 
           
           crfrzf(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*cx(mgs,lr)*dtpinv
           
           ! interpolate along x, i.e., ratio; 
           tmp1 = qiacrratio(i,j) + delx*dqiacrratioinv*(qiacrratio(ip1,j) - qiacrratio(i,j))
           tmp2 = qiacrratio(i,jp1) + delx*dqiacrratioinv*(qiacrratio(ip1,jp1) - qiacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
           qrfrzf(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*qx(mgs,lr)*dtpinv

           ! now subtract off the difference
            crfrzs(mgs) = crfrzs(mgs) - crfrzf(mgs)
            qrfrzs(mgs) = qrfrzs(mgs) - qrfrzf(mgs)

           if ( ipconc >= 6 .and. lzr > 1 ) then
            zrfrzs(mgs) = zrfrz(mgs)
           ! interpolate along x, i.e., ratio; 
            tmp1 = ziacrratio(i,j) + delx*dqiacrratioinv*(ziacrratio(ip1,j) - ziacrratio(i,j))
            tmp2 = ziacrratio(i,jp1) + delx*dqiacrratioinv*(ziacrratio(ip1,jp1) - ziacrratio(i,jp1))
           
           ! interpolate along alpha; 
           
            zrfrzf(mgs) = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))*zx(mgs,lr)*dtpinv
            zrfrzs(mgs) = zrfrzs(mgs) - zrfrzf(mgs)
            zrfrzf(mgs) = (1000./900.)**2*zrfrzf(mgs)
           endif
            endif ! }
           else
            crfrzs(mgs) = 0.0
            qrfrzs(mgs) = 0.0
            zrfrzs(mgs) = 0.0
           endif ! }
           
           endif !}
           
           if ( (qrfrz(mgs))*dtp > qx(mgs,lr) ) then
             fac = ( qrfrz(mgs) )*dtp/qx(mgs,lr)
             qrfrz(mgs) = fac*qrfrz(mgs)
             qrfrzs(mgs) = fac*qrfrzs(mgs)
             qrfrzf(mgs) = fac*qrfrzf(mgs)
             crfrz(mgs) = fac*crfrz(mgs)
             crfrzs(mgs) = fac*crfrzs(mgs)
             crfrzf(mgs) = fac*crfrzf(mgs)
             if ( ipconc >= 6 .and. lzr > 1 ) then
               zrfrz(mgs) = fac*zrfrz(mgs)
               zrfrzf(mgs) = fac*zrfrzf(mgs)
             endif
           endif
           
            endif !}

!           IF ( (crfrzs(mgs) + crfrz(mgs))*dtp > cx(mgs,lr) ) THEN
!             fac = ( crfrzs(mgs) + crfrz(mgs) )*dtp/cx(mgs,lr)
!             crfrz(mgs) = fac*crfrz(mgs)
!             crfrzs(mgs) = fac*crfrzs(mgs)
!           ENDIF
           
!           qrfrzf(mgs) = qrfrz(mgs)
!           crfrzf(mgs) = crfrz(mgs)
           
   !        qrfrz(mgs) = qrfrzf(mgs) + qrfrzs(mgs)
   !        crfrz(mgs) = crfrzf(mgs) + crfrzs(mgs)

           
         elseif ( ibiggopt == 1 ) then
         ! Z85, eq. A34
         tmp = xv(mgs,lr)*brz*cx(mgs,lr)*(Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
         if ( .false. .and. tmp .gt. cxmxd(mgs,lr) ) then ! {
!           write(iunit,*) 'Bigg Freezing problem!',mgs,igs(mgs),kgs(mgs)
!           write(iunit,*)  'tmp, cx(lr), xv = ',tmp, cx(mgs,lr), xv(mgs,lr), (Exp(Max( -arz*temcg(mgs), 0.0 )) - 1.0)
!           write(iunit,*)  'qr,temcg = ',qx(mgs,lr)*1000.,temcg(mgs)
           crfrz(mgs) = cxmxd(mgs,lr) ! cx(mgs,lr)*dtpinv
           qrfrz(mgs) = qxmxd(mgs,lr) ! qx(mgs,lr)*dtpinv
!           STOP
         else ! } {
         crfrz(mgs) = tmp
 !        crfrzfmx = cx(mgs,lr)*Exp(-4./3.*pi*(40.e-6)**3/xv(mgs,lr))
 !        IF ( crfrz(mgs) .gt. crfrzmx ) THEN
 !          crfrz(mgs) = crfrzmx
 !          qrfrz(mgs) = bfnu*xmas(mgs,lr)*rhoinv(mgs)*crfrzmx
 !          qwcnr(mgs) = cx(mgs,lr) - crfrzmx
 !        ELSE
         if ( lzr < 1 ) then
           if ( imurain == 3 ) then
             bfnu = bfnu0
           else !imurain == 1
             bfnu = bfnu1
           endif
         else
 !         bfnu = 1.0 ! (alpha(mgs,lr)+2.0)/(alpha(mgs,lr)+1.)
           if ( imurain == 3 ) then
             bfnu = (alpha(mgs,lr)+2.0)/(alpha(mgs,lr)+1.)
           else !imurain == 1
!             bfnu = bfnu1
            bfnu = (4. + alpha(mgs,lr))*(5. + alpha(mgs,lr))*(6. + alpha(mgs,lr))/  &
     &            ((1. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(3. + alpha(mgs,lr)))
!            bfnu = 1.
           endif
         endif 
         qrfrz(mgs) = bfnu*xmas(mgs,lr)*rhoinv(mgs)*crfrz(mgs)

         qrfrz(mgs) = Min( qrfrz(mgs), 1.*qx(mgs,lr)*dtpinv ) ! qxmxd(mgs,lr) 
         crfrz(mgs) = Min( crfrz(mgs), 1.*cx(mgs,lr)*dtpinv ) !cxmxd(mgs,lr) 
         qrfrz(mgs) = Min( qrfrz(mgs), qx(mgs,lr) )
         qrfrzf(mgs) = qrfrz(mgs)
         endif !}

         
         
         
         if ( crfrz(mgs) .gt. qxmin(lh) ) then !{ Yes, it compares cx and qxmin, but this is just to be sure that 
                                                  ! crfrz is greater than zero in the division
!          IF ( xdia(mgs,lr,1) .lt. 200.e-6 ) THEN
!           IF ( xv(mgs,lr) .lt. xvmn(lh) ) THEN
           
           if ( (ibiggsnow == 1 .or. ibiggsnow == 3 ) .and. ibiggopt /= 2 ) then
           xvfrz = rho0(mgs)*qrfrz(mgs)/(crfrz(mgs)*900.) ! mean volume of frozen drops; 900. for frozen drop density
           frach = 0.5 *(1. +  Tanh(0.2e12 *( xvfrz - 1.15*xvmn(lh))))

             qrfrzs(mgs) = (1.-frach)*qrfrz(mgs)
             crfrzs(mgs) = (1.-frach)*crfrz(mgs) ! *rzxh(mgs)
!             qrfrzf(mgs) = frach*qrfrz(mgs)
           
           endif
           
           if ( ipconc .ge. 14 .and. 1.e-3*rho0(mgs)*qrfrz(mgs)/crfrz(mgs) .lt. xvmn(lh) ) then
             qrfrzs(mgs) = qrfrz(mgs)
             crfrzs(mgs) = crfrz(mgs) ! *rzxh(mgs)
           else
!           crfrz(mgs) = Min( crfrz(mgs), 0.1*cx(mgs,lr)*dtpinv ) ! cxmxd(mgs,lr) 
!           qrfrz(mgs) = Min( qrfrz(mgs), 0.1*qx(mgs,lr)*dtpinv ) ! qxmxd(mgs,lr) 
             qrfrzf(mgs) = frach*qrfrz(mgs)
!             crfrzf(mgs) = Min( qrfrz(mgs)*rho0(mgs)/(xdn(mgs,lh)*vgra), crfrz(mgs) )
            if ( ibfr .le. 1 ) then
             crfrzf(mgs) = frach*Min(crfrz(mgs), qrfrz(mgs)/(bfnu*1.0*vr1mm*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
            elseif ( ibfr .eq. 5 ) then
             crfrzf(mgs) = frach*Min(crfrz(mgs), qrfrz(mgs)/(bfnu*vfrz*1000.0)*rho0(mgs) )*rzxh(mgs)  !*crfrz(mgs)
            elseif ( ibfr .eq. 2 ) then
             crfrzf(mgs) = frach*Min(crfrz(mgs), qrfrz(mgs)/(bfnu*vfrz*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
            elseif ( ibfr .eq. 6 ) then
             crfrzf(mgs) = frach*Max(crfrz(mgs), qrfrz(mgs)/(bfnu*9.*xv(mgs,lr)*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
            else
             crfrzf(mgs) = frach*crfrz(mgs)
            endif 
!             crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*xvmn(lh)*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
!            IF ( lz(lr) > 1 .and. lz(lh) > 1 ) THEN
!              crfrzf(mgs) = crfrz(mgs)
!            ENDIF
            
           endif
!         crfrz(mgs) = Min( cxmxd(mgs,lr), rho0(mgs)*qrfrz(mgs)/xmas(mgs,lr) )
         else
          crfrz(mgs) = 0.0
          qrfrz(mgs) = 0.0
         endif !}

         endif ! ibiggopt

          if ( lvol(lh) .gt. 1 ) then
           vrfrzf(mgs) = rho0(mgs)*qrfrzf(mgs)/rhofrz
          endif

        
        if ( nsplinter .ne. 0 ) then
          if ( nsplinter .ge. 1000 ) then
           ! Lawson et al. 2015 JAS
           ! ave. diam of freezing drops in microns
            tmp = 0
            if ( qrfrz(mgs)*dtp > qxmin(lh) .and. crfrz(mgs) > 1.e-3 ) then
              tmpdiam = 1.e6*( 6.*qrfrz(mgs)/(1000.*pi*crfrz(mgs) ))**(1./3.)  ! avg. diameter of newly frozen drops in microns
              tmp = lawson_splinter_fac*tmpdiam**4*crfrz(mgs)
            endif
          elseif ( nsplinter .gt. 0 ) then
            tmp = nsplinter*crfrz(mgs)
          else
            tmp = -nsplinter*crfrzf(mgs)
          endif
          csplinter2(mgs) = tmp
          qsplinter2(mgs) = Min(0.1*qrfrz(mgs), tmp*splintermass/rho0(mgs) ) ! makes splinters smaller if too much mass is taken from graupel

!          csplinter(mgs) = csplinter(mgs) + tmp
!          qsplinter(mgs) = qsplinter(mgs) + Min(0.1*qrfrz(mgs), tmp*splintermass/rho0(mgs) ) ! makes splinters smaller if too much mass is taken from graupel
        endif
!         IF ( temcg(mgs) .lt. -31.0 ) THEN
!           qrfrz(mgs) = qx(mgs,lr)*dtpinv + qrcnw(mgs)
!           qrfrzf(mgs) = qrfrz(mgs)
!           crfrz(mgs) = cx(mgs,lr)*dtpinv + crcnw(mgs)
!           crfrzf(mgs) = Min(crfrz(mgs), qrfrz(mgs)/(bfnu*1.0*vr1mm*1000.0)*rho0(mgs) ) ! rzxh(mgs)*crfrz(mgs)
!         ENDIF
!         qrfrz(mgs) = 6.0*xdn(mgs,lr)*xv(mgs,lr)**2*tmp*rhoinv(mgs)
!         qrfrz(mgs) = Min( qrfrz(mgs), ffrz*qrmxd(mgs) )
!         crfrz(mgs) = Min( crmxd(mgs), ffrz*crfrz(mgs))
!         crfrz(mgs) = Min(crmxd(mgs),qrfrz(mgs)*rho0(mgs)/xmas(mgs,lr))
       endif
!      if ( temg(mgs) .gt. 268.15 ) then
      else
!      end if
      end if
      end do
      
      endif
!
!  Homogeneous freezing of cloud drops to ice crystals
!  following Bigg (1953) and Ferrier (1994).
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 25b'
      do mgs = 1,ngscnt
      qwfrz(mgs) = 0.0
      cwfrz(mgs) = 0.0
      qwfrzc(mgs) = 0.0
      cwfrzc(mgs) = 0.0
      qwfrzp(mgs) = 0.0
      cwfrzp(mgs) = 0.0
      if ( ibfc .ge. 1 .and. ibfc /= 3 .and. temg(mgs) < 268.15 ) then
!      if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. 1.  .and.   &
!     &     .not. (ipconc .ge. 2 .and. xdia(mgs,lc,1) .lt. 10.e-6) ) then
      if ( qx(mgs,lc) .gt. qxmin(lc) .and. cx(mgs,lc) .gt. cxmin ) then
      if ( ipconc < 2 ) then
      qwfrz(mgs) = ((2.0)*(brz)/(xdn(mgs,lc)*cx(mgs,lc)))   &
     &  *(exp(max(-arz*temcg(mgs), 0.0))-1.0)   &
     &  *rho0(mgs)*(qx(mgs,lc)**2)
      qwfrz(mgs) = max(qwfrz(mgs), 0.0)
      qwfrz(mgs) = min(qwfrz(mgs),qcmxd(mgs))
         cwfrz(mgs) = qwfrz(mgs)*rho0(mgs)/xmas(mgs,li)
       elseif ( ipconc .ge. 2 ) then
         if ( xdia(mgs,lc,3) > 0.e-6 ) then
          volt = exp( 16.2 + 1.0*temcg(mgs) )* 1.0e-6 !  Ts == -temcg ; volt comes from the fit in Fig. 1 in Bigg 1953 
                                               ! for mean temperature for freezing: -ln (V) = a*Ts - b
                                               ! volt is given in cm**3, so factor of 1.e-6 to convert to m**3
!           dbigg = (6./pi* volt )**(1./3.) 

         if (  alpha(mgs,lc) == 0.0 ) then
         cwfrz(mgs) = cx(mgs,lc)*Exp(-volt/xv(mgs,lc))*dtpinv ! number of droplets with volume greater than volt
!turn off limit so that all can freeze at low temp
!!!       cwfrz(mgs) = Min(cwfrz(mgs),ccmxd(mgs))

         qwfrz(mgs) = cwfrz(mgs)*xdn0(lc)*rhoinv(mgs)*(volt + xv(mgs,lc))
          else
            ratio = (1. + alpha(mgs,lc))*volt/xv(mgs,lc)
            
            if ( .false. .and. usegamxinfcnu ) then
              i = Nint(dgami*(1. + alpha(mgs,lc)))
              gcnup1 = gmoi(i)
              i = Nint(dgami*(2. + alpha(mgs,lc)))
              gcnup2 = gmoi(i)

              cwfrz(mgs) = cx(mgs,lc)*Gamxinf(1.+alpha(mgs,lc), ratio)/(dtp*gcnup1) ! gamxinflu(i,j,1,1)

              qwfrz(mgs) = cx(mgs,lc)*xdn0(lc)*xv(mgs,lc)*rhoinv(mgs)*Gamxinf(2.+alpha(mgs,lc), ratio)/(dtp*gcnup2) !  gamxinflu(i,j,12,1)
            
            else
            
              ratio = Min( maxratiolu, ratio )
!              write(0,*) 'cwfrz: temp,ratio = ',temcg(mgs),ratio
!              write(0,*) 'cwfrz: xv,volt,qx = ',xv(mgs,lc),volt,qx(mgs,lc)
!              write(0,*) 'cwfrz: i,j,k = ',igs(mgs),jgs,kgs(mgs)
              tmp = gaminterp(ratio,alpha(mgs,lc),1,1)
!              write(0,*) 'cwfrz: tmp1 = ',tmp
              cwfrz(mgs) = cx(mgs,lc)*tmp*dtpinv ! Gamxinf(1.+alpha(mgs,lc), ratio)/(dtp*gcnup1) ! gamxinflu(i,j,1,1)

              tmp = gaminterp(ratio,alpha(mgs,lc),12,1)
!              write(0,*) 'cwfrz: tmp2 = ',tmp
              qwfrz(mgs) = cx(mgs,lc)*xdn0(lc)*xv(mgs,lc)*rhoinv(mgs)*dtpinv*tmp ! Gamxinf(2.+alpha(mgs,lc), ratio)/(dtp*gcnup2) !  gamxinflu(i,j,12,1)
            
            endif
          
          endif

         endif
       endif
      if ( temg(mgs) .gt. 268.15 ) then
      qwfrz(mgs) = 0.0
      cwfrz(mgs) = 0.0
      end if
      end if
      endif
!
        if ( xplate(mgs) .eq. 1 ) then
          qwfrzp(mgs) = qwfrz(mgs)
          cwfrzp(mgs) = cwfrz(mgs)
        end if
!  
        if ( xcolmn(mgs) .eq. 1 ) then
          qwfrzc(mgs) = qwfrz(mgs)
          cwfrzc(mgs) = cwfrz(mgs)
        end if
      
!
!     qwfrzp(mgs) = 0.0
!     qwfrzc(mgs) = qwfrz(mgs)
!
      end do
!
!
!  Contact freezing nucleation:  factor is to convert from L-1
!  T < -2C:  via Meyers et al. JAM July, 1992 (31, 708-721)
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 25a'
      do mgs = 1,ngscnt

       ccia(mgs) = 0.0

       cwctfz(mgs) = 0.0
       qwctfz(mgs) = 0.0
       ctfzbd(mgs) = 0.0
       ctfzth(mgs) = 0.0
       ctfzdi(mgs) = 0.0

       cwctfzc(mgs) = 0.0
       qwctfzc(mgs) = 0.0
       cwctfzp(mgs) = 0.0
       qwctfzp(mgs) = 0.0
       if ( icfn .ge. 1 ) then

       if ( temg(mgs) .lt. 271.15  .and. qx(mgs,lc) .gt. qxmin(lc)) then

!       find available # of ice nuclei & limit value to max depletion of cloud water

        if ( icfn .ge. 2 ) then
         ccia(mgs) = exp( 4.11 - (0.262)*temcg(mgs) )  ! in m-3, see Walko et al. 1995; 1000*exp(-2.8 -b*t) = exp(6.91)*exp(-2.8 - b*t) = exp(4.11 -b*t)
         !ccia(mgs) = Min(cwctfz(mgs), ccmxd(mgs) )

!       now find how many of these collect cloud water to form IN
!       Cotton et al 1986

         knud(mgs) = 2.28e-5 * temg(mgs) / ( pres(mgs)*raero ) !Walko et al. 1995
         knuda(mgs) = 1.257 + 0.4*exp(-1.1/knud(mgs))          !Pruppacher & Klett 1997 eqn 11-16
         gtp(mgs) = 1. / ( fai(mgs) + fbi(mgs) )               !Byers 65 / Cotton 72b
         dfar(mgs) = kb*temg(mgs)*(1.+knuda(mgs)*knud(mgs))/(6.*pi*fadvisc(mgs)*raero) !P&K 1997 eqn 11-15
         fn1(mgs) = 2.*pi*xdia(mgs,lc,1)*cx(mgs,lc)*ccia(mgs)
         fn2(mgs) = -gtp(mgs)*(ssw(mgs)-1.)*felv(mgs)/pres(mgs)
         fnft(mgs) = 0.4*(1.+1.45*knud(mgs)+0.4*knud(mgs)*exp(-1./knud(mgs)))*(ftka(mgs)+2.5*knud(mgs)*kaero)      &
     &              / ( (1.+3.*knud(mgs))*(2*ftka(mgs)+5.*knud(mgs)*kaero+kaero) )


!      Brownian diffusion
         ctfzbd(mgs) = fn1(mgs)*dfar(mgs)

!      Thermophoretic contact nucleation
         ctfzth(mgs) = fn1(mgs)*fn2(mgs)*fnft(mgs)/rho0(mgs)

!      Diffusiophoretic contact nucleation
         ctfzdi(mgs) = fn1(mgs)*fn2(mgs)*rw*temg(mgs)/(felv(mgs)*rho0(mgs))

         cwctfz(mgs) = max( ctfzbd(mgs) + ctfzth(mgs) + ctfzdi(mgs) , 0.)

!      Sum of the contact nucleation processes
!         IF ( cx(mgs,lc) .gt. 1.e6) write(0,*) 'ctfzbd,etc = ',cwctfz(mgs),ctfzbd(mgs),ctfzth(mgs),ctfzdi(mgs)
!         IF ( wvel(mgs) .lt. -0.05 ) write(6,*) 'ctfzbd,etc = ',ctfzbd(mgs),ctfzth(mgs),ctfzdi(mgs),cx(mgs,lc)*1e-6,wvel(mgs)
!         IF ( ssw(mgs) .lt. 1.0 .and. cx(mgs,lc) .gt. 1.e6 .and. cwctfz(mgs) .gt. 1. ) THEN
!          write(6,*) 'ctfzbd,etc = ',ctfzbd(mgs),ctfzth(mgs),ctfzdi(mgs),cx(mgs,lc)*1e-6,wvel(mgs),fn1(mgs),fn2(mgs)
!          write(6,*) 'more = ',nstep,ssw(mgs),dfar(mgs),gtp(mgs),felv(mgs),pres(mgs)
!         ENDIF

        elseif ( icfn .eq. 1 ) then
         if ( wvel(mgs) .lt. -0.05 ) then ! older kludgy version
           cwctfz(mgs) = cfnfac*exp( (-2.80) - (0.262)*temcg(mgs) )
           cwctfz(mgs) = Min((1.0e3)*cwctfz(mgs), ccmxd(mgs) )  !convert to m-3
         endif
        endif   ! icfn

        if ( ipconc .ge. 2 ) then
         cwctfz(mgs) = Min( cwctfz(mgs)*dtpinv, ccmxd(mgs) )
         qwctfz(mgs) = xmas(mgs,lc)*cwctfz(mgs)/rho0(mgs)
        else
         qwctfz(mgs) = (cimasn)*cwctfz(mgs)/(dtp*rho0(mgs))
         qwctfz(mgs) = max(qwctfz(mgs), 0.0)
         qwctfz(mgs) = min(qwctfz(mgs),qcmxd(mgs))
        endif

!
        if ( xplate(mgs) .eq. 1 ) then
         qwctfzp(mgs) = qwctfz(mgs)
         cwctfzp(mgs) = cwctfz(mgs)
        end if
!
        if ( xcolmn(mgs) .eq. 1 ) then
         qwctfzc(mgs) = qwctfz(mgs)
         cwctfzc(mgs) = cwctfz(mgs)
        end if
        
!        IF ( cwctfz(mgs)*dtp > 0.5 .and. dtp*qwctfz(mgs) > qxmin(li) ) THEN
!          write(91,*) 'cwctfz: ',cwctfz(mgs),qwctfz(mgs) ! ,cwctfzc(mgs),qwctfzc(mgs)
!        ENDIF
       
!
!     qwctfzc(mgs) = qwctfz(mgs)
!     qwctfzp(mgs) = 0.0
!
       end if

       endif ! icfn

      end do
!
!
!
! Hobbs-Rangno ice enhancement (Ferrier, 1994)
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 23a'
      dthr = 300.0
      hrifac = (1.e-3)*((0.044)*(0.01**3))
      do mgs = 1,ngscnt
      ciihr(mgs) = 0.0
      qiihr(mgs) = 0.0
      cicichr(mgs) = 0.0
      qicichr(mgs) = 0.0
      cipiphr(mgs) = 0.0
      qipiphr(mgs) = 0.0
      if ( ihrn .ge. 1 ) then
      if ( qx(mgs,lc) .gt. qxmin(lc) ) then
      if ( temg(mgs) .lt. 273.15 ) then
!      write(iunit,'(3(1x,i3),3(1x,1pe12.5))')
!     : igs(mgs),jgs,kgs(mgs),cx(mgs,lc),rho0(mgs),qx(mgs,lc)
!      write(iunit,'(1pe15.6)')
!     :  log(cx(mgs,lc)*(1.e-6)/(3.0)),
!     :  ((1.e-3)*rho0(mgs)*qx(mgs,lc)),
!     :  (cx(mgs,lc)*(1.e-6)),
!     : ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)),
!     : (alog(cx(mgs,lc)*(1.e-6)/(3.0)) *
!     >  ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)))

      if ( Log(cx(mgs,lc)*(1.e-6)/(3.0)) .gt. 0.0 ) then
      ciihr(mgs) = ((1.69e17)/dthr)   &
     & *(log(cx(mgs,lc)*(1.e-6)/(3.0)) *   &
     &  ((1.e-3)*rho0(mgs)*qx(mgs,lc))/(cx(mgs,lc)*(1.e-6)))**(7./3.)
      ciihr(mgs) = ciihr(mgs)*(1.0e6)
      qiihr(mgs) = hrifac*ciihr(mgs)/rho0(mgs)
      qiihr(mgs) = max(qiihr(mgs), 0.0)
      qiihr(mgs) = min(qiihr(mgs),qcmxd(mgs))
      endif
!
      if ( xplate(mgs) .eq. 1 ) then
      qipiphr(mgs) = qiihr(mgs)
      cipiphr(mgs) = ciihr(mgs)
      end if
!
      if ( xcolmn(mgs) .eq. 1 ) then
      qicichr(mgs) = qiihr(mgs)
      cicichr(mgs) = ciihr(mgs)
      end if
!
!     qipiphr(mgs) = 0.0
!     qicichr(mgs) = qiihr(mgs)
!
      end if
      end if
      endif ! ihrn
      end do
!
!
!
!  simple frozen rain to hail conversion.  All of the
!  frozen rain larger than 5.0e-3 m in diameter are converted
!  to hail.  This is done by considering the equation for
!  frozen rain mixing ratio:
!
!
!  qfw = [ cno(lf) * pi * fwdn / (6 rhoair) ]
!
!         /inf
!      *  |     fwdia*3 exp(-dia/fwdia) d(dia)
!         /Do
!
!  The amount to be reclassified as hail is the integral above from
!  Do to inf where Do is 5.0e-3 m.
!
!
!  qfauh = [ cno(lf) * pi * fwdn / (6 rhoair) ]
!
!


      hdia0 = 300.0e-6
      do mgs = 1,ngscnt
      qscnvi(mgs) = 0.0
      cscnvi(mgs) = 0.0
      cscnvis(mgs) = 0.0
!      IF ( .false. ) THEN
!      IF ( temg(mgs) .lt. tfr .and. ssi(mgs) .gt. 1.01 .and. qx(mgs,li) .gt. qxmin(li) ) THEN
      if ( temg(mgs) .lt. tfr .and. qx(mgs,li) .gt. qxmin(li) ) then
        if ( ipconc .ge. 4 .and. .false. ) then
         if ( cx(mgs,li) .gt. 10. .and. xdia(mgs,li,1) .gt. 50.e-6 ) then !{
         cirdiatmp =   &
     &  (qx(mgs,li)*rho0(mgs)   &
     & /(pi*xdn(mgs,li)*cx(mgs,li)))**(1./3.)
          if ( cirdiatmp .gt. 100.e-6 ) then !{
          qscnvi(mgs) =   &
     &  ((pi*xdn(mgs,li)*cx(mgs,li)) / (6.0*rho0(mgs)*dtp))   &
     & *exp(-hdia0/cirdiatmp)   &
     & *( (hdia0**3) + 3.0*(hdia0**2)*cirdiatmp   &
     &  + 6.0*(hdia0)*(cirdiatmp**2) + 6.0*(cirdiatmp**3) )
      qscnvi(mgs) =   &
     &  min(qscnvi(mgs),qimxd(mgs))
          if ( ipconc .ge. 4 ) then
            cscnvi(mgs) = Min( cimxd(mgs), cx(mgs,li)*Exp(-hdia0/cirdiatmp))
          endif
         endif  ! }
        end if ! }

       elseif ( ipconc .lt. 4 ) then

        qscnvi(mgs) = 0.001*eii(mgs)*max((qx(mgs,li)-1.e-3),0.0)
        qscnvi(mgs) = min(qscnvi(mgs),qxmxd(mgs,li))
        cscnvi(mgs) = qscnvi(mgs)*rho0(mgs)/xmas(mgs,li)
        cscnvis(mgs) = 0.5*cscnvi(mgs)

       endif
      endif
!      ENDIF
      end do



!
!  Ventilation coeficients
!
      do mgs = 1,ngscnt
      fvent(mgs) = (fschm(mgs)**(1./3.)) * (fakvisc(mgs)**(-0.5))
      end do
!
!
      if ( ndebug .gt. 0 ) write(0,*) 'civent'
!
      civenta = 1.258e4
      civentb = 2.331
      civentc = 5.662e4
      civentd = 2.373
      civente = 0.8241
      civentf = -0.042
      civentg = 1.70

      do mgs = 1,ngscnt
      if ( icond .eq. 1 .or. temg(mgs) .le. tfrh    &
     &      .or. (qx(mgs,lr) .le. qxmin(lr) .and. qx(mgs,lc) .le. qxmin(lc)) ) then
      if ( qx(mgs,li) .gt. qxmin(li) ) then
      cireyn =   &
     &  (civenta*xdia(mgs,li,1)**civentb   &
     &  +civentc*xdia(mgs,li,1)**civentd)   &
     &  /   &
     &  (civente*xdia(mgs,li,1)**civentf+civentg)
      xcivent = (fschm(mgs)**(1./3.))*((cireyn/fakvisc(mgs))**0.5)
      if ( xcivent .lt. 1.0 ) then
      civent(mgs) = 1.0 + 0.14*xcivent**2
      end if
      if ( xcivent .ge. 1.0 ) then
      civent(mgs) = 0.86 + 0.28*xcivent
      end if
      else
       civent(mgs) = 0.0
      endif


      endif ! icond .eq. 1
      end do

!
!
      igmrwa = 100.0*2.0
      igmrwb = 100.*((5.0+br)/2.0)
      rwventa = (0.78)*gmoi(igmrwa)  ! 0.78
      rwventb = (0.308)*gmoi(igmrwb) ! 0.562825
      do mgs = 1,ngscnt
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then
        if ( ipconc .ge. 3 ) then
          if ( imurain == 3 ) then
           if ( izwisventr == 1 ) then
            rwvent(mgs) = ventrx(mgs)*(1.6 + 124.9*(1.e-3*rho0(mgs)*qx(mgs,lr))**.2046)
           else ! izwisventr = 2
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
          rwvent(mgs) =   &
     &  (0.78*ventrx(mgs) + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )
           endif

          else ! imurain == 1
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lr)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        if ( iferwisventr == 1 ) then

  ! Ferrier fall speed in the ventillation term [uses fx(lr) ]
  
        alpr = Min(alpharmax,alpha(mgs,lr) )

        x =  1. + alpha(mgs,lr)

        if ( ipconc >= 6 .and. lzr > 1 ) then ! 3 moment
        tmp = 1. + alpr ! alpha(mgs,lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lr) + 0.5*bx(lr)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions
        else
         y = ventrxn(mgs)
        endif

!         vent1 = dble(xdia(mgs,lr,1))**(-2. - alpr) ! Actually OK
!         vent2 = dble(1./xdia(mgs,lr,1) + 0.5*fx(lr))**dble(2.5+alpr+0.5*bx(lr))  ! Actually OK
         vent1 = dble(xdia(mgs,lr,1))**(0.5 + 0.5*bx(lr)) ! 2016.2.26 Changed for consistency with derivation (recast formula -- should be equivalent)
         vent2 = dble(1. + 0.5*fx(lr)*xdia(mgs,lr,1))**dble(2.5+alpr+0.5*bx(lr))
        
        
        rwvent(mgs) =    &
     &    0.78*x +    &
     &    0.308*fvent(mgs)*y*   &
     &            Sqrt(ax(lr)*rhovt(mgs))*(vent1/vent2)
       
        rwventz(mgs) = 0.0

!        rwventz(mgs) =    &
!     &    0.78*x +    &
!     &    0.308*fvent(mgs)*y*   &
!     &            Sqrt(ax(lr)*rhovt(mgs))*(vent1/vent2)


        elseif ( iferwisventr == 2 ) then
          
!  Following Wisner et al. (1972) but using gamma of volume. Note that Ferrier rain fall speed does not integrate with gamma of volume, so using Vr = ar*d^br
         x =  1. + alpha(mgs,lr)

           rwvent(mgs) =   &
     &  (0.78*x + 0.308*ventrxn(mgs)*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )


        if ( ipconc >= 7 ) then
        alpr = Min(alpharmax,alpha(mgs,lr) )

           tmp = alpr + 5.5 + br/2.
           i = Int(dgami*(tmp))
           del = tmp - dgam*i
           y = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

!        rwventz(mgs) =    &
!     &    0.78*(4. + alpha(mgs,lr))*(3. + alpha(mgs,lr))*(2. + alpha(mgs,lr))*(1. + alpha(mgs,lr)) +    &
        rwventz(mgs) =    &
     &    0.78*(4. + alpr)*(3. + alpr)*(2. + alpr)*(1. + alpr) +    &
     &    0.308*fvent(mgs)*   &
     &            Sqrt(ax(lr)*rhovt(mgs))*(y/gf1palp(mgs))*(xdia(mgs,lr,1)**((1.0+br)/2.0))

        endif

          
          endif ! iferwisventr
          
          endif ! imurain
        else
         rwvent(mgs) =   &
     &  (rwventa + rwventb*fvent(mgs)   &
     &   *Sqrt((ar*rhovt(mgs)))   &
     &    *(xdia(mgs,lr,1)**((1.0+br)/2.0)) )
        endif
      else
       rwvent(mgs) = 0.0
      endif
      end do
!
      igmswa = 100.0*2.0
      igmswb = 100.*((5.0+ds)/2.0)
      swventa = (0.78)*gmoi(igmswa)
      swventb = (0.308)*gmoi(igmswb)
      do mgs = 1,ngscnt
      if ( qx(mgs,ls) .gt. qxmin(ls) ) then
      if ( ipconc .ge. 4 ) then
      swvent(mgs) = 0.65 + 0.44*fvent(mgs)*Sqrt(vtxbar(mgs,ls,1)*xdia(mgs,ls,1))
      else
! 10-ice version:
       swvent(mgs) =   &
     &  (swventa + swventb*fvent(mgs)   &
     &   *Sqrt((cs*rhovt(mgs)))   &
     &   *(xdia(mgs,ls,1)**((1.0+ds)/2.0)) )
      endif
      else
      swvent(mgs) = 0.0
      endif
      end do
!
!

      igmhwa = 100.0*2.0
      igmhwb = 100.0*2.75
      hwventa = (0.78)*gmoi(igmhwa)
      hwventb = (0.308)*gmoi(igmhwb)
!      hwventc = (4.0*gr/(3.0*cdx(lh)))**(0.25)
      hwvent(:) = 0.0
      hwventy(:) = 0.0

      do mgs = 1,ngscnt
      if ( qx(mgs,lh) .gt. qxmin(lh) ) then
       hwventc = (4.0*gr/(3.0*cdxgs(mgs,lh)))**(0.25)
       if ( .false. .or. alpha(mgs,lh) .eq. 0.0 ) then
        hwvent(mgs) =   &
     &  ( hwventa + hwventb*hwventc*fvent(mgs)   &
     &    *((xdn(mgs,lh)/rho0(mgs))**(0.25))   &
     &    *(xdia(mgs,lh,1)**(0.75)))
       else ! Ferrier 1994, eq. B.36
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lh)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
        
! note that hwvent includes a division by Gamma(1+alpha), so Gamma(2+alpha)/Gamma(1+alpha) = 1 + alpha
! and g1palp = Gamma(1+alpha) divides into y
        x =  1. + alpha(mgs,lh)

        tmp = 1 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        
        hwventy(mgs) = 0.308*fvent(mgs)*(xdia(mgs,lh,1)**(0.5 + 0.5*bxx(mgs,lh)))*Sqrt(axx(mgs,lh)*rhovt(mgs)) 
        hwvent(mgs) =    &
     &  ( 0.78*x +  y*hwventy(mgs) ) !   &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lh,1)**(0.5 + 0.5*bxx(mgs,lh)))*   &
!     &            Sqrt(axx(mgs,lh)*rhovt(mgs)) )
       
       endif
      else
      hwvent(mgs) = 0.0
      hwventy(mgs) = 0.0
      endif
      end do
      

      hlvent(:) = 0.0
      hlventy(:) = 0.0

      if ( lhl .gt. 1 ) then
      igmhwa = 100.0*2.0
      igmhwb = 100.0*2.75
      hwventa = (0.78)*gmoi(igmhwa)
      hwventb = (0.308)*gmoi(igmhwb)
!      hwventc = (4.0*gr/(3.0*cdx(lhl)))**(0.25)
      do mgs = 1,ngscnt
      if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
      hwventc = (4.0*gr/(3.0*cdxgs(mgs,lhl)))**(0.25)

       if ( .false. .or. alpha(mgs,lhl) .eq. 0.0 ) then
        hlvent(mgs) =   &
     &  ( hwventa + hwventb*hwventc*fvent(mgs)   &
     &    *((xdn(mgs,lhl)/rho0(mgs))**(0.25))   &
     &    *(xdia(mgs,lhl,1)**(0.75)))
       else ! Ferrier 1994, eq. B.36
       ! linear interpolation of complete gamma function
!        tmp = 2. + alpha(mgs,lhl)
!        i = Int(dgami*(tmp))
!        del = tmp - dgam*i
!        x = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

! note that hlvent includes a division by Gamma(1+alpha), so x = Gamma(2+alpha)/Gamma(1+alpha) = 1 + alpha
! and g1palp = Gamma(1+alpha) divides into y

        x =  1. + alpha(mgs,lhl)

        tmp = 1 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp ! ratio of gamma functions

        hlventy(mgs) = 0.308*fvent(mgs)*(xdia(mgs,lhl,1)**(0.5 + 0.5*bxx(mgs,lhl)))*Sqrt(axx(mgs,lhl)*rhovt(mgs)) 
        
        hlvent(mgs) =  0.78*x + y*hlventy(mgs)  !   &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lhl,1)**(0.5 + 0.5*bxx(mgs,lhl)))*   &
!     &            Sqrt(axx(mgs,lhl)*rhovt(mgs)))
!     :            Sqrt(xdn(mgs,lhl)*ax(lhl)*rhovt(mgs)/rg0))/tmp

        endif
       endif
      end do
      endif

!
!
!
!  Wet growth constants
!
      do mgs = 1,ngscnt
      fwet1(mgs) =   &
     & (2.0*pi)*   &
     & ( felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qss0(mgs)-qx(mgs,lv))   &
     &  -ftka(mgs)*temcg(mgs) )   &
     & / ( rho0(mgs)*(felf(mgs)+fcw(mgs)*temcg(mgs)) )
      fwet2(mgs) =   &
     &  (1.0)-fci(mgs)*temcg(mgs)   &
     & / ( felf(mgs)+fcw(mgs)*temcg(mgs) )
      end do
!
!  Melting constants
!
      do mgs = 1,ngscnt
      fmlt1(mgs) = (2.0*pi)*   &
     &  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))   &
     &   -ftka(mgs)*temcg(mgs)/rho0(mgs) )    &
     &  / (felf(mgs))
      fmlt2(mgs) = -fcw(mgs)*temcg(mgs)/felf(mgs)
      fmlt1e(mgs) = (2.0*pi)*   &
     &  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))  ) / (felf(mgs))
      end do
!
!  Vapor Deposition constants
!
      do mgs = 1,ngscnt
      fvds(mgs) =    &
     &  (4.0*pi/rho0(mgs))*(ssi(mgs)-1.0)*   &
     &  (1.0/(fai(mgs)+fbi(mgs)))
      end do
      do mgs = 1,ngscnt
      fvce(mgs) =    &
     &  (4.0*pi/rho0(mgs))*(ssw(mgs)-1.0)*   &
     &  (1.0/(fav(mgs)+fbv(mgs)))
      end do

!
!  deposition, sublimation, and melting of snow, graupel and hail
!
      qsmlr(:) = 0.0
      qimlr(:) = 0.0 ! this is not used. qi melts to qc way down in the code.
      qhmlr(:) = 0.0
      qhlmlr(:) = 0.0
      if ( lhwlg > 1 ) then
        qhmlrlg(:) = 0.0
        qhlmlrlg(:) = 0.0
      endif
      qhfzh(:) = 0.0
      qffzf(:) = 0.0
      qhlfzhl(:) = 0.0
      qhfzhlg(:) = 0.0
      qhlfzhllg(:) = 0.0
      vhfzh(:) = 0.0
      vffzf(:) = 0.0
      vhlfzhl(:) = 0.0
      qsfzs(:) = 0.0
!      zsmlr(:) = 0.0
      zhmlr(:) = 0.0
      zhmlrr(:) = 0.0
      zsmlrr(:) = 0.0
      zhshr(:) = 0.0
      zhlmlr(:) = 0.0
      zhlshr(:) = 0.0

      zhshrr(:) = 0.0
      zhlmlrr(:) = 0.0
      zhlshrr(:) = 0.0

      csmlr(:) = 0.0
      csmlrr(:) = 0.0
      chmlr(:) = 0.0
      chmlrr(:) = 0.0
      chlmlr(:) = 0.0
      chlfmlr(:) = 0.0
!      chlmlrsave(:) = 0.0
!      qhlmlrsave(:) = 0.0
!      chlsave(:) = 0.0
!      qhlsave(:) = 0.0
      chlmlrr(:) = 0.0


      if ( .not. mixedphase ) then !{
      do mgs = 1,ngscnt
!
      if ( temg(mgs) .gt. tfr ) then
      
      if (  qx(mgs,ls) .gt. qxmin(ls) ) then
      qsmlr(mgs) =   &
     &   min(   &
     &  (c1sw*fmlt1(mgs)*cx(mgs,ls)*swvent(mgs)*xdia(mgs,ls,1) ) & ! /rhosm    &
     &   , 0.0 )
      endif

      
!       IF ( qx(mgs,ls) .gt. 0.1e-4 ) write(0,*) 'qsmlr: ',qsmlr(mgs),qx(mgs,ls),cx(mgs,ls),fmlt1(mgs),
!     :        temcg(mgs),swvent(mgs),xdia(mgs,ls,1),qss0(mgs)-qx(mgs,lv)
!      ELSE
!       qsmlr(mgs) = 0.0
!      ENDIF
! 10ice version:
!     >   min(
!     >  (fmlt1(mgs)*cx(mgs,ls)*swvent(mgs)*xdia(mgs,ls,1) +
!     >   fmlt2(mgs)*(qsacr(mgs)+qsacw(mgs)) )
!     <   , 0.0 )

      if (  qx(mgs,lh) .gt. qxmin(lh) ) then

      if ( ibinhmlr == 0 .or. lzh < 1 ) then
       qhmlr(mgs) =   &
     &   meltfac*min(   &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1)   &
     &  + fmlt2(mgs)*(qhacrmlr(mgs)+qhacwmlr(mgs))    &
     &   , 0.0 )
       elseif ( ibinhmlr == 1 ) then ! use incomplete gamma functions to approximate the bin results


       qhmlr(mgs) =   &
     &   min(   &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1)   &
!     &  + fmlt2(mgs)*(qhacrmlr(mgs)+qhacw(mgs))    & ! turn off the collection part for now
     &   , 0.0 )

! first part to integrate from mltdiam1 and mltdiam2 to infinity (k loop)

        tmp = 1. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g4palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
 
       do k = 1,numshedregimes ! do not need mltdiam4 -- treat as infinity; goes to 3 to cover the top 3 diameters (independent of ndiam, the number of smaller diameters)
        
!        ratio = mltdiam1/xdia(mgs,lh,1)
        ratio = Min( maxratiolu, mltdiam(ndiam+k)/xdia(mgs,lh,1) )

        
        if ( usegamxinf2 ) then
           x = gamxinfdp(2. + alpha(mgs,lh), ratio)/g1palp
        else
           x = gaminterp(ratio,alpha(mgs,lh),2,1)
        endif ! usegamxinf

! qxd(ndiam+2), cxd(ndiam+2), qhml(ndiam+2)        
        if ( usegamxinf2 ) then
           cxd1 = cx(mgs,lh)*(1. - gamxinfdp(1. + alpha(mgs,lh), ratio)/g1palp)
        else

          tmp = gaminterp(ratio,alpha(mgs,lh),1,1)

!          IF ( tmp > 0.95 ) THEN
!           tmp = gamxinfdp(1. + alpha(mgs,lh), ratio)/g1palp
!          ENDIF
          cxd1 = cx(mgs,lh)*(1. - tmp)
        endif

         cxd(ndiam+k) = cxd1

        if ( usegamxinf2 ) then
         qxd1 = qx(mgs,lh)*(1. - gamxinfdp(4. + alpha(mgs,lh), ratio)/g4palp)
        else

          tmp = gaminterp(ratio,alpha(mgs,lh),4,1)

!           IF ( tmp > 0.95 ) THEN
!            tmp = gamxinfdp(4. + alpha(mgs,lh), ratio)/g4palp
!           ENDIF
           
           qxd1 = qx(mgs,lh)*(1. - tmp)
        endif

         qxd(ndiam+k) = qxd1

        if ( usegamxinf2 ) then
          y =  gamxinfdp(2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio)/g1palp
        else

          y = gaminterp(ratio,alpha(mgs,lh),3,1)

        endif
        
          hwvent1 =  0.78*x + y*hwventy(mgs) 

          qhmlr1 = min( fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1), 0.0 )

! zero to mltdiam1 is then the difference between qhmlr and qhmlr1
       qhmlr0 = qhmlr(mgs) - qhmlr1
       qhml(ndiam+k) = qhmlr1 ! from diam to infinity
       qhml0(ndiam+k) = qhmlr0 ! from zero to diam

!        IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
!          write(iunit,*) 'diam loop k,qhmlr1 = ',k,qhmlr1,qxd1,cxd1
!          write(iunit,*) 'fmlt1 = ',fmlt1(mgs),cx(mgs,lh),hwvent1,x,y,hwventy(mgs)
!        
!        ENDIF

       
       enddo ! k

         if ( numshedregimes == 2 ) then

           dqhml(ndiam+2) = qhml(ndiam+2) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam4
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhmlr1 = qhml(ndiam+1)
           qhmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2

         
         elseif ( numshedregimes == 3 ) then
         
           dqhml(ndiam+3) = qhml(ndiam+3) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam4
           dqhml(ndiam+2) = qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
           qhmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhmlr1 = qhml(ndiam+1)
           qhmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           qhmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
         
         endif


         if ( .not. mixedphase ) qhmlrlg(mgs) = qhmlr1
         

        if ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) then
 !         write(iunit,*) 'k,qhml0(n+1,+2,+3) = ',kgs(mgs),qhml0(ndiam+1),qhml0(ndiam+2),qhml0(ndiam+3)
 !         write(iunit,*) 'k,qhml(n+1,+2,+3) = ',kgs(mgs),qhml(ndiam+1),qhml(ndiam+2),qhml(ndiam+3)
 !         write(iunit,*) 'qhmlr1,qhmlr2,qhmlr12 = ',qhmlr1,qhmlr2,qhmlr12,qhml(ndiam+1)-qhml(ndiam+2)
 !         write(iunit,*) 'qhmlr,qhmlrlg = ',qhmlr(mgs),qhmlrlg(mgs)
        endif


! First calc of number rate from 0 to D1
! numdiam is the number of diameter "bins" for finding the number loss at diameters less than mltdiam1 (9mm)
! If alpha > 1.5, this will be replaced with integrated value

        do k = 1,numdiam !{
        
!        ratio = mltdiam05/xdia(mgs,lh,1)
        ratio = Min(maxratiolu, mltdiam(k)/xdia(mgs,lh,1) )

        
        if ( usegamxinf2 ) then
          x =  gamxinf(2. + alpha(mgs,lh), ratio)/g1palp
        else

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh))))
!           j = Int(Max(alphamin,Min(alphamax,alpha(mgs,lh)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,2,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,2,1) - gamxinflu(i,j,2,1))
           tmp2 = gamxinflu(i,jp1,2,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,2,1) - gamxinflu(i,jp1,2,1))
           
           ! interpolate along alpha; 

          x = gaminterp(ratio,alpha(mgs,lh),2,1)

         
           x = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
        endif
        ! number and mass from 0 to mltdiam05

!        IF ( usegamxinf2 .or. k == 1 ) THEN
        if ( usegamxinf2 ) then
          cxd05 = cx(mgs,lh)*(1. - gamxinf(1. + alpha(mgs,lh), ratio)/g1palp)

          if ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) then
           tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
           tmp2 = cx(mgs,lh)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
!           IF ( tmp3 > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lh)*(1. - tmp)
           cxd05 = tmp2
          endif
        else

           ! interpolate along x, i.e., ratio; 
!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
!           tmp1 = gamxinflu(i,j,1,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,1,1) - gamxinflu(i,j,1,1))
!           tmp2 = gamxinflu(i,jp1,1,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,1,1) - gamxinflu(i,jp1,1,1))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
          if ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) then
          cxd05 = cx(mgs,lh)*(1. - gamxinf(1. + alpha(mgs,lh), ratio)/g1palp)
           tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
           tmp2 = cx(mgs,lh)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
!           IF ( tmp3 > 5. ) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
          
          endif

           tmp = gaminterp(ratio,alpha(mgs,lh),1,1)

           if ( tmp > 0.98 ) then ! Interpolated values can be flaky close to 1.0, so revert to expensive gamxinf subroutine
             cxd05 = cx(mgs,lh)*(1. - gamxinf(1. + alpha(mgs,lh), ratio)/g1palp)
           else
             cxd05 = cx(mgs,lh)*(1. - tmp)
           endif

        endif
        
         cxd(k) = cxd05


!        IF ( usegamxinf2 .or. k == 1 ) THEN
        if ( usegamxinf2 ) then
          qxd05 = qx(mgs,lh)*(1. - gamxinf(4. + alpha(mgs,lh), ratio)/g4palp)

!           tmp = gaminterp(ratio,alpha(mgs,lh),4,1)
!           tmp2 = qx(mgs,lh)*(1. - tmp)
!           tmp3 = 2.*100.*(qxd05 - tmp2)/(qxd05 + tmp2)
        
           if ( Abs(tmp3) > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) then
!            write(0,*) 'k,qxd05,lu = ',k,igs(mgs),kgs(mgs),qxd05,tmp2, tmp3
!            write(0,*) 'ratio,alph,gam,interp = ', ratio,alpha(mgs,lh),1.0 - gamxinf(4. + alpha(mgs,lh), ratio)/g4palp, 1.0 - tmp
           endif
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lh)*(1. - tmp)
!          qxd05 = tmp2


        else
           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,4,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,4,1) - gamxinflu(i,j,4,1))
!           tmp2 = gamxinflu(i,jp1,4,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,4,1) - gamxinflu(i,jp1,4,1))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lh),4,1)

           if ( tmp > 0.98 ) then ! Interpolated values can be flaky close to 1.0, so revert to expensive gamxinf subroutine
             qxd05 = qx(mgs,lh)*(1. - gamxinf(4. + alpha(mgs,lh), ratio)/g4palp)
           else
             qxd05 = qx(mgs,lh)*(1. - tmp)
           endif

        endif
        
        qxd(k) = qxd05

        if ( usegamxinf ) then
          y = gamxinf(2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio)/g1palp
        else
           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,3,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,3,1) - gamxinflu(i,j,3,1))
           tmp2 = gamxinflu(i,jp1,3,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,3,1) - gamxinflu(i,jp1,3,1))
           
           ! interpolate along alpha; 
           
           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lh),3,1)

            y =  tmp ! 
        endif


        
        ! reuse hwvent1
        hwvent1 =  0.78*x + y*hwventy(mgs)


       qhml(k) = Min( fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1), 0.0 ) ! from diam to inf.
       qhml0(k) = qhmlr(mgs) - qhml(k)                                    ! from zero to diam
       
!       qhml0(k) = qhmlr05
! zero to 0.5*mltdiam1 is then the difference between qhmlr and qhmlr05
!       qhmlr05 = qhmlr(mgs) - qhmlr05
        if ( k == 1 ) then
         dqhml(k) = qhml0(k) ! qhmlr(mgs) - qhml0(k)
        else 
         dqhml(k) = qhml0(k) - qhml0(k-1)
        endif
       
       enddo ! } k loop over diameters less than 9mm

! compute between last diameter and 9mm
       dqhml(numdiam+1) = qhml0(ndiam+1) - qhml0(numdiam)
!       cxd(numdiam+1) = 
       cxd(numdiam+1) = cxd(ndiam+1) ! - cxd(numdiam)
       qxd(numdiam+1) = qxd(ndiam+1) ! - qxd(numdiam)

! compute qx and cx for each 'bin'
! sum the number of melted particle in tmp

     dcxd(1) = cxd(1)
     dqxd(1) = qxd(1)
     tmpcmlt =  dqhml(1)*dcxd(1)/dqxd(1)
     
     if ( .false. .and. igs(mgs) == 52 ) then
       write(0,*) 'HML: kgs, tmpcmlt = ',kgs(mgs),tmpcmlt,dqhml(1),dcxd(1),dqxd(1)
       write(0,*) 'dqhml+1,+2 = ',dqhml(ndiam+1),dqhml(ndiam+2)
       write(0,*) 'qhml0+1,+2 = ',qhml0(ndiam+1),qhml0(ndiam+2)
       write(0,*) 'qhml+1,+2 = ',qhml(ndiam+1),qhml(ndiam+2)
       write(0,*) 'qhmlr(mgs) = ',qhmlr(mgs)
     endif
     
     do k = 2,numdiam+1
      dcxd(k) = cxd(k) - cxd(k-1)
      dqxd(k) = qxd(k) - qxd(k-1)
!     IF ( .false. .and. igs(mgs) == 52 ) THEN
!       write(0,*) 'HML: kgs, k,delta- tmpcmlt = ',kgs(mgs),k,dqhml(k)*dcxd(k)/dqxd(k), dqhml(k),dcxd(k),dqxd(k),cxd(k),cxd(k-1)
!       write(0,*) cxd(ndiam+1), cxd(ndiam+1)-cxd(k-1)
!       write(0,*) qxd(ndiam+1), qxd(ndiam+1)-qxd(k-1)
!     ENDIF
      if ( dqxd(k) > qxmin(lh) ) then
        tmpcmlt = tmpcmlt + dqhml(k)*dcxd(k)/dqxd(k)
      endif
     enddo

     ! now, qhmlr1 and qhmlr2 lose mass but not number
     ! qhmlr0 loses both mass and number

         ! first term is melting for d < d1 divided by mean mass of particles with d < d1
         ! second and third terms are the shedding from larger particles
!         chmlrr(mgs) = qhmlr0/(qxd1/cxd1) + rho0(mgs)*(qhmlr2*mltmass1inv + qhmlr12*mltmass2inv)



! test integrated N rate for graupel

        if ( alpha(mgs,lh) > 1.5 ) then
        
        ratio = Min( maxratiolu, mltdiam1/xdia(mgs,lh,1) )

        tmp = 1.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami



        tmp = -1. + alpha(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf ) then
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(-1. + alpha(mgs,lh), ratio))/g1palp
          x = (Gamma_sp(-1. + alpha(mgs,lh)) - gamxinf(-1. + alpha(mgs,lh), ratio))/g1palp
        else

          
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,7,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,7,1) - gamxinflu(i,j,7,1))
           tmp2 = gamxinflu(i,jp1,7,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,7,1) - gamxinflu(i,jp1,7,1))
           
           ! interpolate along alpha; 
           
!           tmp3 = (tmp1 + dely*(tmp2 - tmp1))*Gamma_sp(5. + alpha(mgs,lh))
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
!           IF ( tmp3 > 0.98 ) THEN
!            tmp3 = gamxinf(5. + alpha(mgs,lh), ratio)/g1palp
!           ENDIF


          !  x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
            
            x = tmp3  ! changed lookup table to be the difference (gam - gamxinf)/g1palp
        
        endif

        tmp = -0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf ) then
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(-0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
           y = (Gamma_sp(-0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)) - gamxinf(-0.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
        else
          ! FIX THIS FOR NUMBER!

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,8,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,8,1) - gamxinflu(i,j,8,1))
           tmp2 = gamxinflu(i,jp1,8,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,8,1) - gamxinflu(i,jp1,8,1))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))
           
        !   y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3
           
        endif

        hwvent1 =  0.78*x + y*hwventy(mgs)
        
        chmlrtmp =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1/xdia(mgs,lh,1)**(2)   &
     &   , 0.0 )/xdn(mgs,lh)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lh))/pi comes from removing the mass/volume factors; 
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
         
!          IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
!           write(iunit,*) 'k,chmlrtmp',kgs(mgs),chmlrtmp,hwvent1,cx(mgs,lh),fmlt1(mgs),xdia(mgs,lh,1)
!           write(iunit,*) 'x,y,hwventy ',x,y,hwventy(mgs),bxx(mgs,lh),bx(lh)
!         ENDIF
        
        tmp = (chmlrtmp )
        
        tmpcmlt = chmlrtmp
        
!        IF ( igs(mgs) == 40 ) THEN
!          write(iunit,*) 'chmlrtmp,chmlrtmpd1inf = ',chmlrtmp
!          write(iunit,*) 'tmpcmlt old/new, alph = ',tmpcmlt,tmp,alpha(mgs,lh)
!        ENDIF
       
       endif

! test using mltdiam05
!         chmlrr(mgs) = (qhml0(ndiam+1)-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhmlr05/(qxd05/cxd05) + rho0(mgs)*(qhmlr2*mltmass1inv + qhmlr12*mltmass2inv)
!         tmpcmlt = (qhml0(ndiam+1)-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhmlr05/(qxd05/cxd05)
         if ( numshedregimes == 2 ) then
           chmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhmlr2*mltmass1inv + qhmlr12*mltmass2inv)

         elseif ( numshedregimes == 3 ) then
!           qhmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
!           qhmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
!           qhmlr1 = qhml(ndiam+1)
!           qhmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
!           qhmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           chmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhmlr3*mltmass1inv + qhmlr23*mltmass2inv + qhmlr12*mltmass3inv)
         endif

!         IF ( igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) THEN
!           write(iunit,*) 'k,qh,alp,qhmlr,chmlrr(mgs),xdia ',kgs(mgs),1.e3*qx(mgs,lh),alpha(mgs,lh),qhmlr(mgs),chmlrr(mgs),xdia(mgs,lh,3)
!           write(iunit,*) 'zx,tmpcmlt,qhmlr2,qhmlr12,xdn ',zx(mgs,lh),tmpcmlt,qhmlr2,qhmlr12,xdn(mgs,lh),qhmlr12*mltmass2inv,qhmlr2*mltmass1inv
!         ENDIF

         if ( .true. .or. alpha(mgs,lh) < alphamax - delta_alphamlr .or. xdia(mgs,lh,3) > binmlrmxdia ) then ! { allow spectrum to get narrower
     !     chmlr(mgs)  = qhmlr0/(rho0(mgs)*qxd1/cxd1) 
          chmlr(mgs)  =  tmpcmlt ! (qhml0(ndiam+1)-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhmlr05/(qxd05/cxd05)
          

! test integrated Z rate for graupel

        tmp = 1.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp

        tmp = 5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        hwvent1 =  0.78*x + y*hwventy(mgs)
!        hwvent1 =    &
!     &  ( 0.78*x +    &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lh,1)**(0.5 + 0.5*bxx(mgs,lh)))*   &
!     &            Sqrt(axx(mgs,lh)*rhovt(mgs)) )
       
       if ( .false. ) then
       ! full shedding version
       zhmlr(mgs) =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
     &   , 0.0 )/xdn(mgs,lh)*(6.*2./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lh))*2/pi comes from removing the mass/volume factors; 
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 24
         ! factor of rho0 to nullify the 1/rho0 in fmlt1

! alternate zhmlr assumes Z = N D**6, from which dZ/dt = (dN/dt)*D**6 = (1/M) dM/dt D**6
! turns out to be exactly a factor of 2 smaller than the other zhmlr (and matches with bin results when set to remove full particles
! from each bin, i.e., particles do not shrink)
       else
       ! whole particle loss version: factor of 1/2 different from shedding version
       zhmlr(mgs) =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
     &   , 0.0 )/xdn(mgs,lh)*(6./pi)  
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
       endif
       
       zhmlr0inf = zhmlr(mgs)

! note that fmlt1 has the factor of 2*pi
!        fmlt1(mgs) = (2.0*pi)*   &
!     &  ( felv(mgs)*fwvdf(mgs)*(qss0(mgs)-qx(mgs,lv))   &
!     &   -ftka(mgs)*temcg(mgs)/rho0(mgs) )    &
!     &  / (felf(mgs))


          if ( .false. .and. igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) then
          tmp1 = zhmlr0inf
          
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
!          zhmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
          zhmlrtmp =  g1*rho0(mgs)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs)   - tmp**2 * qhmlr(mgs)/tmp  )
          
          write(iunit,*) 'zhmlrtmp,zhmlr,1.5*zhmlr: ',zhmlrtmp,zhmlr(mgs),1.5*zhmlr(mgs),zhmlrtmp/(zhmlr(mgs)),alpha(mgs,lh)
          
         ! try average:
          
    !      zhmlr(mgs) = 0.5*(zhmlr(mgs) + zhmlrtmp)
          
          endif

!        IF ( .false. .and. igs(mgs) == 40 ) THEN !.and. zhmlr(mgs) < -zxmin  THEN
! #include 'bintest.F90'  
!        ENDIF


          !  try using the new N and q with the old Z to compute a new alpha, then use new alpha (assuming it is larger)
          ! to compute delta-Z?

          il = lh
          chw = cx(mgs,il) + chmlr(mgs)
          qr  = qx(mgs,il) + qhmlr(mgs)
          z   = zx(mgs,il)  ! + zhmlr(mgs)

          if ( .false. .and. z .gt. zxmin .and. qr > qxmin(lh) .and. chw > cxmin ) then
           
           il = lh
           
           tmp2 = alpha(mgs,lh) ! save alpha
           
            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = Max( alphamin, Min( alphamax, alp ) )
           
         if ( newton ) then
           do i = 1,10
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = alp + ( galpha(alp) - rdi )/dgalpha(alp)
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo
           
         else
           do i = 1,10
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo
          endif
          
            
            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            
            tmp = qx(mgs,lh)/cx(mgs,lh)
            
            zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
            
            alpha(mgs,lh) = tmp2 ! restore alpha
            
          endif ! everything greater than min values

      endif !}
!
! not sure this is right -- should it use g1x(mgs,lh) or some assumed value for melt rain?
!
! should separately add the reflecitivities from the 2 sizes of shed drops, then add the completely melted particles

!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs)  ) )


! test integrated Z rate for rain for size range up to mltdiam1
! integrate from d1 to infinity, then subtract that from the full value (zhmlr)
        ratio = Min( maxratiolu, mltdiam1/xdia(mgs,lh,1) )


        tmp = 1. + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5. + alpha(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf ) then
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lh), ratio))/g1palp
          x = (Gamma_sp(5. + alpha(mgs,lh)) - gamxinf(5. + alpha(mgs,lh), ratio))/g1palp
        else

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,5,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,5,1) - gamxinflu(i,j,5,1))
           tmp2 = gamxinflu(i,jp1,5,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,5,1) - gamxinflu(i,jp1,5,1))
           
           ! interpolate along alpha; 
           
!           tmp3 = (tmp1 + dely*(tmp2 - tmp1))*Gamma_sp(5. + alpha(mgs,lh))
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
           tmp3 = gaminterp(ratio,alpha(mgs,lh),5,1)
           
!           IF ( tmp3 > 0.98 ) THEN
!            tmp3 = gamxinf(5. + alpha(mgs,lh), ratio)/g1palp
!           ENDIF


!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - tmp3)/g1palp
!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)*(1. - tmp3)
       ! x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
         x = tmp3 ! changed lookup table to be the difference (gam - gamxinf)/g1palp
       !  x = Max(0.0,x)
         
        endif

        tmp = 5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf ) then
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
           y = (Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))  - gamxinf(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio))/g1palp
        else
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,1) - gamxinflu(i,j,6,1))
           tmp2 = gamxinflu(i,jp1,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,1) - gamxinflu(i,jp1,6,1))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))
           
           
           tmp3 = gaminterp(ratio,alpha(mgs,lh),6,1)
           
         !  y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3 ! changed lookup table to be the difference (gam - gamxinf)/g1palp
           y = Max(0.0, y)
        endif
        if ( .false. .and. igs(mgs) == 38 .and. qhmlr(mgs) < 0.0 ) then
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lh)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lh)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lh) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,1) - gamxinflu(i,j,6,1))
           tmp2 = gamxinflu(i,jp1,6,1) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,1) - gamxinflu(i,jp1,6,1))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1)) ! *Gamma_sp(5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh))
           
          write(iunit,*) 'ibinhmlr = 1: i,k = ',igs(mgs),kgs(mgs)
          write(iunit,*) 'y,ratio,alpha,mltdiam1,xdia(mgs,lh,1) = ',y,ratio,alpha(mgs,lh),mltdiam1,xdia(mgs,lh,1)
          write(iunit,*) 'full gamx = ',gamxinf(5.5 +0.5*bxx(mgs,lh) + alpha(mgs,lh), ratio),bxx(mgs,lh),bx(lh)
           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
          write(iunit,*) 'interpolated y = ', y,tmp1,tmp2,tmp3, &
               gamxinflu(i,j,6,1),gamxinflu(ip1,j,6,1),gamxinflu(i,jp1,6,1),gamxinflu(ip1,jp1,6,1)
           
        endif

        hwvent1 =  0.78*x + y*hwventy(mgs)
       
       zhmlrr(mgs) =   &  !  -- THIS IS THE LOWER INCOMPLETE GAMMA... 0 to ratio
     &   rho0(mgs)*min(   &
     &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
     &   , 0.0 )/xdn(mgs,lh)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lh))/pi comes from removing the mass/volume factors; 

! zhmlrr here is whole particle loss for d < mltdiam1
!
! so now have zhmlrr as the Z loss (whole particle) for graupel for 0 <= d <= mltdiam1, and zhmlr0inf as the integrated whole-particle loss
! We want whole particle loss for 0 < d < mltdiam1 which is (zhmlrr(mgs)), plus shedding loss for d > mltdiam1, which is 2*(zhmlr0inf - zhmlrr(mgs)) 
! Shedding loss is simply double the whole particle loss

       zhmlr(mgs) = zhmlrr(mgs) + binmlrzrrfac*2.*(zhmlr0inf - zhmlrr(mgs)) ! factor of 2 to convert whole particle loss to shedding loss
       ! should that be zhmlr(mgs) = 2.*zhmlrr(mgs) + (zhmlr0inf - zhmlrr(mgs)) = zhmlrr(mgs) + zhmlr0inf
       ! yes, I think it should -- NOPE!
       ! zhmlr(mgs) = 2.*zhmlrr(mgs) + (zhmlr0inf - zhmlrr(mgs))

          if ( .false. .and. igs(mgs) == 40 .and. qhmlr(mgs) < -qxmin(lh) ) then
          
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
!          zhmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
          zhmlrtmp =  g1*rho0(mgs)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs)   - tmp**2 * qhmlr(mgs)/tmp  )

        if ( alpha(mgs,lh) >= alphamax - delta_alphamlr ) then
        tmp = 1.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5.0 + alpha(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp

        tmp = 5.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        hwvent1 =  0.78*x + y*hwventy(mgs)

       ! whole particle loss version: factor of 1/2 different from shedding version
 !      zhmlr0inf =   &
 !    &   min(   rho0(mgs)* &
 !    &  fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1)**4   &
 !    &   , 0.0 )/xdn(mgs,lh)*(6./pi)  
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
        
        tmp = zhmlrr(mgs) + 2.*(zhmlr0inf - zhmlrr(mgs)) 
        
       endif
          
!          write(iunit,*) 'zhmlrtmp,zhmlr,zhmlr0inf: ',zhmlrtmp,zhmlr(mgs),zhmlr0inf,zhmlrtmp/(zhmlr(mgs)),alpha(mgs,lh),tmp,zhmlrtmp/(zhmlr0inf)
                    
          endif

!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr0 - tmp**2 * chmlr(mgs)  ) )

         if ( numshedregimes == 2 ) then
           zhmlrr(mgs) = (xdn(mgs,lh)/1000.)**2*zhmlrr(mgs) + rho0(mgs)*   &
     &                   (qhmlr2*mltmass1inv*(0.01*takshedsize1)**6 + qhmlr12*mltmass2inv*(0.01*takshedsize2)**6 )
         elseif ( numshedregimes == 3 ) then
!           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
          zhmlrr(mgs) = (xdn(mgs,lh)/1000.)**2 *zhmlrr(mgs) + rho0(mgs)*  &
     &              (qhmlr3*mltmass1inv*(0.01*takshedsize1)**6 + qhmlr23*mltmass2inv*(0.01*takshedsize2)**6  + &
     &               qhmlr12*mltmass3inv*(0.01*takshedsize3)**6 )
         endif


         if ( alpha(mgs,lh) >= alphamax - delta_alphamlr ) then ! { cannot allow spectrum to get any narrower, 
         
          chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
         
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
          zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )

!          chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh)))  ! into rain 

!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs)  ) )

!          IF(imltshddmr == 1 ) THEN
!            ! DTD: If Dmg < sheddiam, then assume complete melting into
!            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
!            tmp = -rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh))) ! Min of Maximum raindrop size/mean hail size
!            tmp2 = -rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
!            chmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lh,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lh,3)-sheddiam)/(sheddiam0-sheddiam)
!            chmlrr(mgs) = -Max(tmp,Min(tmp2,chmlrr(mgs)))
!          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
!            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
!            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
!            chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lh))  ! into rain 
!          ENDIF

         elseif (  xdia(mgs,lh,3) > binmlrmxdia ) then ! 
         
          chmlrtmp  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lh)/cx(mgs,lh)
          zhmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrtmp )
          
          tmp1 = 1.e-3
          if ( xdia(mgs,lh,3) > binmlrmxdia + tmp1 ) then 
             chmlr(mgs) = chmlrtmp
             zhmlr(mgs) = zhmlrtmp
          else ! weighted average
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'k,xdia,chmlr,zhmlr = ',kgs(mgs),xdia(mgs,lh,3)*1.e3,chmlr(mgs),zhmlr(mgs)
!              write(0,*) 'weights 1,2 = ',(xdia(mgs,lh,3) - binmlrmxdia )/tmp1 , (binmlrmxdia + tmp1 - xdia(mgs,lh,3))/tmp1
!            ENDIF
            chmlr(mgs) = chmlrtmp*(xdia(mgs,lh,3) - binmlrmxdia )/tmp1 + chmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lh,3))/tmp1
            zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
            !zhlmlr(mgs) = zhlmlrtmp*(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 + zhlmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'tmp,chlmlr,= ',chlmlrtmp,chlmlr(mgs)
!              write(0,*) 'tmp,zhlmlr,= ',zhlmlrtmp,zhlmlr(mgs)
!            ENDIF
          endif
          

         endif ! }

           
!         IF ( igs(mgs) == 38 ) THEN
!           x = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
!           write(iunit,*) 'k,chmlrr: ',kgs(mgs),chmlrr(mgs),qx(mgs,lh),qhmlr(mgs),zhmlrr(mgs),qxd05,qxd1,cxd05,cxd1
!!           write(iunit,*) 'k,zhmlr: ',kgs(mgs),zhmlr(mgs),x,qx(mgs,lh),qhmlr(mgs),zhmlrr(mgs),xdia(mgs,lh,1),zx(mgs,lh)
!!           write(iunit,*) 'zhmlrr parts: ',  2.*tmp * qhmlr(mgs), - tmp**2 * chmlrr(mgs) , 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs) 
!           write(iunit,*) 'chmlrr parts: ',(qhmlr0-qhmlr05)/((qxd1-qxd05)/(cxd1-cxd05)), qhmlr05/(qxd05/cxd05), rho0(mgs)*(qhmlr2*mltmass1inv), rho0(mgs)*( qhmlr12*mltmass2inv)
!           write(iunit,*) 'qhmlr0: ',qhmlr0,qhmlr05,qhmlr2,qhmlr12
!         ENDIF
         
       elseif ( ibinhmlr == 2 .or. ibinhmlr == 3 ) then

       endif
       
       
       if ( ivhmltsoak > 0 .and. qhmlr(mgs) < 0.0 .and. lvol(lh) > 1 .and. xdn(mgs,lh) .lt. xdnmx(lh) ) then
         ! act as if 100% of the meltwater were soaked into the graupel
           v1 = (1. - xdn(mgs,lh)/xdnmx(lh))*(vx(mgs,lh) + rho0(mgs)*qhmlr(mgs)/xdn(mgs,lh) )/(dtp) ! volume available for filling
           v2 = -1.0*rho0(mgs)*qhmlr(mgs)/xdnmx(lh)  ! volume of melted ice if it were refrozen in the matrix
           
           vhsoak(mgs) = Min(v1,v2)
           
       endif

      endif !  qx(mgs,lh) .gt. qxmin(lh)

      
      if ( lhl .gt. 1  .and. lhlw < 1 ) then

       if ( qx(mgs,lhl) .gt. qxmin(lhl) ) then
         if ( ibinhlmlr == 0  .or. lzhl < 1) then
       qhlmlr(mgs) =   &
     &   meltfac*min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*hlvent(mgs)*xdia(mgs,lhl,1)   &
     &  + fmlt2(mgs)*(qhlacrmlr(mgs)+qhlacwmlr(mgs))    &
     &   , 0.0 )

       elseif ( ibinhlmlr == 1 ) then ! use incomplete gamma functions to approximate the bin results

! #ifdef 1
! #if (defined 1) && defined( COMMAS ) || defined( COMMASTMP )

       qhlmlr(mgs) =   &
     &   min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*hlvent(mgs)*xdia(mgs,lhl,1)   &
!     &  + fmlt2(mgs)*(qhacrmlr(mgs)+qhacw(mgs))    & ! turn off the collection part for now
     &   , 0.0 )

! first part to integrate from mltdiam1 and mltdiam2 to infinity (k loop)

        tmp = 1. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 4. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g4palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
 
       do k = 1,numshedregimes  ! do not need k=3, treat as infinity; goes to 3 to cover the top 3 diameters (independent of ndiam, the number of smaller diameters)
        
!        ratio = mltdiam1/xdia(mgs,lhl,1)
        ratio = Min(maxratiolu, mltdiam(ndiam+k)/xdia(mgs,lhl,1) )

        
        if ( usegamxinf3 ) then
        x =  gamxinfdp(2. + alpha(mgs,lhl), ratio)/g1palp
        else

          x = gaminterp(ratio,alpha(mgs,lhl),2,2)
          
          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = gamxinfdp(2. + alpha(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck01: i,k,xint,x = ',igs(mgs),kgs(mgs),x,tmpgam
            endif
          endif

       endif ! usegamxinf

! qxd(ndiam+2), cxd(ndiam+2), qhml(ndiam+2)        
        if ( usegamxinf3 ) then
         cxd1 = cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
        else

          tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)

          
!          IF ( tmp > 0.95 ) THEN
!           tmp = gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp
!          ENDIF

          cxd1 = cx(mgs,lhl)*(1. - tmp)

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
              write(6,*) 'gamcheck02: i,k,cxd1,x = ',igs(mgs),kgs(mgs),cxd1,tmpgam
            endif
          endif


        endif

         cxd(ndiam+k) = cxd1

        if ( usegamxinf3 ) then
         qxd1 = qx(mgs,lhl)*(1. - gamxinfdp(4. + alpha(mgs,lhl), ratio)/g4palp)
        else

          tmp = gaminterp(ratio,alpha(mgs,lhl),4,2)

!           IF ( tmp > 0.95 ) THEN
!            tmp = gamxinfdp(4. + alpha(mgs,lhl), ratio)/g4palp
!           ENDIF

           qxd1 = qx(mgs,lhl)*(1. - tmp)

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = qx(mgs,lhl)*(1. - gamxinfdp(4. + alpha(mgs,lhl), ratio)/g4palp)
             ! tmp1 =  qx(mgs,lhl)*(1. - gaminterp(ratio,alpha(mgs,lhl),10,2))
              write(6,*) 'gamcheck03: i,k,qxd1,x = ',igs(mgs),kgs(mgs),qxd1,tmpgam !,tmp1
             if ( .false. ) then
              write(6,*) 'alpha,gaminterp = ',alpha(mgs,lhl),tmp
              write(6,*) 'ratio,ratio2,mltdiam,xdia = ',ratio,mltdiam(ndiam+k)/xdia(mgs,lhl,1),mltdiam(ndiam+k),xdia(mgs,lhl,1) 
              tmp1=gamxinfdp(4. + alpha(mgs,lhl), ratio)
              tmp2=gamma_dpr(4. + alpha(mgs,lhl))
              write(6,*) 'gam,gam4,gam4alt,gam/gam4 = ',tmp1,g4palp,tmp2,tmp1/g4palp
           alp = alpha(mgs,lhl)
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(maxalphalu,alp))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           endif
           delx = Min(maxratiolu,ratio) - float(i)*dqiacrratio
           dely = Min(maxalphalu, alp ) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )
           
           luindex = 4
           il = 2

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,luindex,il) + delx*dqiacrratioinv*         &
     &                 (gamxinflu(ip1,j,luindex,il) - gamxinflu(i,j,luindex,il))
           tmp2 = gamxinflu(i,jp1,luindex,il) + delx*dqiacrratioinv*       &
     &                 (gamxinflu(ip1,jp1,luindex,il) - gamxinflu(i,jp1,luindex,il))
           
           ! interpolate along alpha; 
           
           temp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
              write(6,*) 'i,j,temp3,gam1,',i,j,temp3, gamxinflu(i,j,luindex,il)
            endif ! false/true
            
            endif
          endif

        endif

         qxd(ndiam+k) = qxd1

        if ( usegamxinf3 ) then
          y =  gamxinfdp(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
        else

          y = gaminterp(ratio,alpha(mgs,lhl),3,2)

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = gamxinfdp(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck04: i,k,yint,y = ',igs(mgs),kgs(mgs),y,tmpgam
            endif
          endif

        endif
        
          hwvent1 =  0.78*x + y*hlventy(mgs) 

          qhlmlr1 = min( fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1), 0.0 )

        if ( igs(mgs) == 40 ) then
!          write(0,*) 'kgs,k,usegamxinf/2/3 = ',kgs(mgs),k,usegamxinf,usegamxinf2,usegamxinf3
!          write(0,*) 'diam loop k,qhmlr1 = ',k,qhmlr1,qxd1,cxd1
!          write(0,*) 'fmlt1 = ',fmlt1(mgs),cx(mgs,lhl),hwvent1,x,y,hwventy(mgs)
!          write(0,*) 'ratio0,ratio,xdia(1) = ',mltdiam(ndiam+k)/xdia(mgs,lhl,1),ratio,xdia(mgs,lhl,1)
        
        endif

! zero to mltdiam1 is then the difference between qhlmlr and qhlmlr1
       qhlmlr0 = qhlmlr(mgs) - qhlmlr1
       qhml(ndiam+k) = qhlmlr1 ! from diam to infinity
       qhml0(ndiam+k) = qhlmlr0 ! from zero to diam
        
!        IF ( igs(mgs) == 40 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!          write(iunit,*) 'diam loop k,qhlmlr1 = ',k,qhlmlr1,qxd1,cxd1
!          write(iunit,*) 'fmlt1 = ',fmlt1(mgs),cx(mgs,lhl),hwvent1,x,tmp3,y,hlventy(mgs)
!          write(iunit,*) 'lookup: ',i,j,gamxinflu(i,j,2,1),gamxinflu(i,j,2,2) !, gaminterp(ratio,alpha(mgs,lhl),2,-2)
!        
!        ENDIF
       
       enddo ! k = 1,2

         if ( numshedregimes == 2 ) then
           dqhml(ndiam+2) = qhml(ndiam+2) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhlmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhlmlr1 = qhml(ndiam+1)
           qhlmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
         
         elseif ( numshedregimes == 3 ) then
         
           dqhml(ndiam+3) = qhml(ndiam+3) ! qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+2) = qhml0(ndiam+3) - qhml0(ndiam+2) ! amount of melting for  19mm < D < mltdiam3
           dqhml(ndiam+1) = qhml0(ndiam+2) - qhml0(ndiam+1) ! amount of melting for 9mm < D < 19mm
         
           qhlmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
           qhlmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
           qhlmlr1 = qhml(ndiam+1)
           qhlmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           qhlmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
         
         endif

         if ( .not. mixedphase ) qhlmlrlg(mgs) = qhlmlr1
         

!  integrate from mltdiam1/2 to infinity

        do k = 1,numdiam
        
!        ratio = mltdiam05/xdia(mgs,lhl,1)
        ratio = Min( maxratiolu,  mltdiam(k)/xdia(mgs,lhl,1) )

        
        if ( usegamxinf3 ) then
          x =  gamxinf(2. + alpha(mgs,lhl), ratio)/g1palp
        else

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl))))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,2,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,2,2) - gamxinflu(i,j,2,2))
           tmp2 = gamxinflu(i,jp1,2,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,2,2) - gamxinflu(i,jp1,2,2))
           
           ! interpolate along alpha; 

         !  x = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

          x = gaminterp(ratio,alpha(mgs,lhl),2,2)

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = gamxinfdp(2.0 + alpha(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck05: i,k,xint,x = ',igs(mgs),kgs(mgs),x,tmpgam
            endif
          endif

         
        endif
        ! number and mass from 0 to mltdiam05

!        IF ( usegamxinf2 .or. k == 1 ) THEN
        if ( usegamxinf3 ) then
          cxd05 = cx(mgs,lhl)*(1. - gamxinf(1. + alpha(mgs,lhl), ratio)/g1palp)

          if ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) then
           tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)
           tmp2 = cx(mgs,lhl)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
       !    IF ( tmp3 > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lhl)*(1. - tmp)
          cxd05 = tmp2
          endif
        else

           ! interpolate along x, i.e., ratio; 
!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
!           tmp1 = gamxinflu(i,j,1,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,1,2) - gamxinflu(i,j,1,2))
!           tmp2 = gamxinflu(i,jp1,1,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,1,2) - gamxinflu(i,jp1,1,2))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
          if ( .false. .and. igs(mgs) == 52 .and. kgs(mgs) == 6 ) then
          cxd05 = cx(mgs,lhl)*(1. - gamxinf(1. + alpha(mgs,lhl), ratio)/g1palp)
           tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)
           tmp2 = cx(mgs,lhl)*(1. - tmp)
           tmp3 = 2.*100.*Abs(cxd05 - tmp2)/(cxd05 + tmp2)
        
!           IF ( tmp3 > 5. ) write(0,*) 'k,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05,tmp2, tmp3
          
          endif

           tmp = gaminterp(ratio,alpha(mgs,lhl),1,2)

           if ( tmp > 0.98 ) then
             cxd05 = cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
           else
             cxd05 = cx(mgs,lhl)*(1. - tmp)
           endif

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = gamxinfdp(1.0 + alpha(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck06: i,k,tmp,x = ',igs(mgs),kgs(mgs),tmp,tmpgam
              tmp1 =  cx(mgs,lhl)*(1. - gamxinfdp(1. + alpha(mgs,lhl), ratio)/g1palp)
              write(6,*) 'cxd05,gam: ',cxd05,tmp1
            endif
          endif
          
          cxd05 = Max(0.0, cxd05)

        endif
        
         cxd(k) = cxd05


!        IF ( usegamxinf2 .or. k == 1 ) THEN
        if ( usegamxinf3 ) then
          qxd05 = qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)

!           tmp = gaminterp(ratio,alpha(mgs,lhl),4,2)
!           tmp2 = qx(mgs,lhl)*(1. - tmp)
!           tmp3 = 2.*100.*(qxd05 - tmp2)/(qxd05 + tmp2)
        
           if ( Abs(tmp3) > 50. .and. (cxd05 > 0.1 .or. tmp2 > 0.1 )) then
!            write(0,*) 'k,qxd05,lu = ',k,igs(mgs),kgs(mgs),qxd05,tmp2, tmp3
!            write(0,*) 'ratio,alph,gam,interp = ', ratio,alpha(mgs,lhl),1.0 - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp, 1.0 - tmp
           endif
!            write(0,*) 'k,igs,kgs,cxd05,lu = ',k,igs(mgs),kgs(mgs),cxd05, cx(mgs,lhl)*(1. - tmp)
!          qxd05 = tmp2


        else
           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,4,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,4,2) - gamxinflu(i,j,4,2))
!           tmp2 = gamxinflu(i,jp1,4,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,4,2) - gamxinflu(i,jp1,4,2))
           
           ! interpolate along alpha; 
           
!           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lhl),4,2)

           if ( tmp > 0.98 ) then
             qxd05 = qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
           else
             qxd05 = qx(mgs,lhl)*(1. - tmp)
           endif

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp
              write(6,*) 'gamcheck07: i,k,tmp,x = ',igs(mgs),kgs(mgs),tmp,tmpgam
              tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
              write(6,*) 'qxd05,gam: ',qxd05,tmp1
            endif
          endif


        endif
        
        qxd(k) = qxd05

        if ( usegamxinf3 ) then
          y = gamxinf(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
        else
           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,3,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,3,2) - gamxinflu(i,j,3,2))
           tmp2 = gamxinflu(i,jp1,3,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,3,2) - gamxinflu(i,jp1,3,2))
           
           ! interpolate along alpha; 
           
           tmp = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))

           tmp = gaminterp(ratio,alpha(mgs,lhl),3,2)

            y =  tmp ! 

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = gamxinf(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
              write(6,*) 'gamcheck08: i,k,tmp,y = ',igs(mgs),kgs(mgs),tmp,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            endif
          endif

        endif


        
        ! reuse hwvent1
        hwvent1 =  0.78*x + y*hlventy(mgs)


       qhml(k) = Min( fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1), 0.0 ) ! from diam to inf.
       qhml0(k) = qhlmlr(mgs) - qhml(k)                                    ! from zero to diam
 
 
        if ( igs(mgs) == 40 ) then
!          write(0,*) '2 kgs,k,usegamxinf/2/3 = ',kgs(mgs),k,usegamxinf,usegamxinf2,usegamxinf3
!          write(0,*) '2 diam loop k,qhml,qhml0 = ',k,qhml(k),qhml0(k) ,qxd05,cxd05
!          write(0,*) '2 fmlt1 = ',fmlt1(mgs),cx(mgs,lhl),hwvent1,x,y,hwventy(mgs)
!          write(0,*) '2 ratio0,ratio,xdia(1) = ',mltdiam(k)/xdia(mgs,lhl,1),ratio,xdia(mgs,lhl,1),alpha(mgs,lhl)
        
        endif
      
!       qhml0(k) = qhlmlr05
! zero to 0.5*mltdiam1 is then the difference between qhlmlr and qhlmlr05
!       qhlmlr05 = qhlmlr(mgs) - qhlmlr05
        if ( k == 1 ) then
         dqhml(k) = qhml0(k) ! qhlmlr(mgs) - qhml0(k)
        else 
         dqhml(k) = qhml0(k) - qhml0(k-1)
        endif
       
       enddo ! k loop over diameters less than 9mm

! compute between last diameter and 9mm
       dqhml(numdiam+1) = qhml0(ndiam+1) - qhml0(numdiam)
!       cxd(numdiam+1) = 
       cxd(numdiam+1) = cxd(ndiam+1) ! - cxd(numdiam)
       qxd(numdiam+1) = qxd(ndiam+1) ! - qxd(numdiam)

! compute qx and cx for each 'bin'
! sum the number of melted particle in tmp

     dcxd(1) = cxd(1)
     dqxd(1) = qxd(1)
     tmpcmlt =  dqhml(1)*dcxd(1)/dqxd(1)
     
     if ( .false. .and. igs(mgs) == 52 ) then
       write(0,*) 'HML: kgs, tmpcmlt = ',kgs(mgs),tmpcmlt,dqhml(1),dcxd(1),dqxd(1)
       write(0,*) 'dqhml+1,+2 = ',dqhml(ndiam+1),dqhml(ndiam+2)
       write(0,*) 'qhml0+1,+2 = ',qhml0(ndiam+1),qhml0(ndiam+2)
       write(0,*) 'qhml+1,+2 = ',qhml(ndiam+1),qhml(ndiam+2)
       write(0,*) 'qhlmlr(mgs) = ',qhlmlr(mgs)
     endif
     
     do k = 2,numdiam+1
      dcxd(k) = cxd(k) - cxd(k-1)
      dqxd(k) = qxd(k) - qxd(k-1)
     if ( .false. .and. igs(mgs) == 52 ) then
       write(0,*) 'HML: kgs, k,delta- tmpcmlt = ',kgs(mgs),k,dqhml(k)*dcxd(k)/dqxd(k), dqhml(k),dcxd(k),dqxd(k),cxd(k),cxd(k-1)
       write(0,*) cxd(ndiam+1), cxd(ndiam+1)-cxd(k-1)
       write(0,*) qxd(ndiam+1), qxd(ndiam+1)-qxd(k-1)
     endif
      if ( dqxd(k) > qxmin(lhl) ) then
        tmpcmlt = tmpcmlt + dqhml(k)*dcxd(k)/dqxd(k) ! this value of tmpcmlt is used if alpha(lhl) <= 1.5
      endif
     enddo

     ! now, qhlmlr1 and qhlmlr2 lose mass but not number
     ! qhlmlr0 loses both mass and number

         ! first term is melting for d < d1 divided by mean mass of particles with d < d1
         ! second and third terms are the shedding from larger particles
!         chlmlrr(mgs) = qhlmlr0/(qxd1/cxd1) + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)

! test using mltdiam05
!         chlmlrr(mgs) = (qhml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05) + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)
!         tmp = (qhml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05)
!         chlmlrr(mgs) = tmp + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)



! test integrated N rate for Hail

        if ( alpha(mgs,lhl) > 1.5 ) then ! For lower alpha, the estimate of tmpcmlt from above is used
        
        ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lhl,1) )

        tmp = 1.0 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = -1. + alpha(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf3 ) then
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(-1. + alpha(mgs,lhl), ratio))/g1palp
          x = (Gamma_sp(-1. + alpha(mgs,lhl)) - gamxinf(-1. + alpha(mgs,lhl), ratio))/g1palp
        else

          
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,7,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,7,2) - gamxinflu(i,j,7,2))
           tmp2 = gamxinflu(i,jp1,7,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,7,2) - gamxinflu(i,jp1,7,2))
           
           ! interpolate along alpha; 
           
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
           
!           IF ( tmp3 > 0.98 ) THEN ! use direct calculation in this nonlinear region
!            tmp3 = gamxinf(5. + alpha(mgs,lhl), ratio)/g1palp
!           ENDIF

          !  x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
            x = tmp3 ! changed lookup table to be the difference already:  Gamma_sp(-1+alp) - Gamma_sp(-1+alp,ratio)

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = (Gamma_sp(-1. + alpha(mgs,lhl)) - gamxinf(-1. + alpha(mgs,lhl), ratio))/g1palp
              tmp1 = gaminterp(ratio,alpha(mgs,lhl),7,2)
              write(6,*) 'gamcheck09: i,k,tmp3,interp,x = ',igs(mgs),kgs(mgs),tmp3,tmp1,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            endif
          endif
        
        endif

        tmp = -0.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf3 ) then
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(-0.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
           y = (Gamma_sp( alpha(mgs,lhl) - 0.5 + 0.5*bxx(mgs,lhl) )  - gamxinf( alpha(mgs,lhl) - 0.5 + 0.5*bxx(mgs,lhl), ratio))/g1palp
        else

           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,8,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,8,2) - gamxinflu(i,j,8,2))
           tmp2 = gamxinflu(i,jp1,8,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,8,2) - gamxinflu(i,jp1,8,2))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lh))
           
        !   y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3
           
          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = (Gamma_sp(-0.5 + alpha(mgs,lhl)) - gamxinf(-0.5 + alpha(mgs,lhl), ratio))/g1palp
              tmp1 = gaminterp(ratio,alpha(mgs,lhl),8,2)
              write(6,*) 'gamcheck10: i,k,tmp3,interp,x = ',igs(mgs),kgs(mgs),tmp3,tmp1,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            endif
          endif
        endif

        hwvent1 =  0.78*x + y*hlventy(mgs)
        
        chlmlrtmp =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lhl)*hwvent1/xdia(mgs,lhl,1)**(2)   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lhl))/pi comes from removing the mass/volume factors; 
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1
         
!          IF ( igs(mgs) == 40 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!           write(iunit,*) 'k,chlmlrtmp',kgs(mgs),chlmlrtmp,hwvent1,cx(mgs,lhl),fmlt1(mgs),xdia(mgs,lhl,1)
!           write(iunit,*) 'x,y,hlventy ',x,y,hlventy(mgs),bxx(mgs,lhl),bx(lhl)
!         ENDIF
       
        tmp = (chlmlrtmp )
        
        tmpcmlt = chlmlrtmp
        
!        IF ( igs(mgs) == 40 ) THEN
!          write(iunit,*) 'chlmlrtmp,chlmlrtmpd1inf = ',chlmlrtmp
!          write(iunit,*) 'tmpcmlt old/new, alph = ',tmpcmlt,tmp,alpha(mgs,lhl)
!        ENDIF
       
       endif

! test using mltdiam05
!         chlmlrr(mgs) = (qhlml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05) + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)
!         tmpcmlt = (qhlml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05)
         if ( numshedregimes == 2 ) then
           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr2*mltmass1inv + qhlmlr12*mltmass2inv)
         elseif ( numshedregimes == 3 ) then
!           qhlmlr3 = qhml(ndiam+3) ! dqhml(ndiam+2)
!           qhlmlr2 = qhml(ndiam+2) ! dqhml(ndiam+2)
!           qhlmlr1 = qhml(ndiam+1)
!           qhlmlr12 = qhml(ndiam+1)-qhml(ndiam+2) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
!           qhlmlr23 = qhml(ndiam+2)-qhml(ndiam+3) ! dqhml(ndiam+1) ! qhmlr1 - qhmlr2
           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
         endif
         
         ! hack!
         !  chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)*mltmass0inv
         !   IF ( igs(mgs) == nx/2 ) THEN
         !     write(6,*) 'i,k,tmpcmlt,othercmlt = ',igs(mgs),kgs(mgs),tmpcmlt,rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
         !     write(6,*) 'mlr0,1,2,3: ', qhlmlr(mgs)-qhlmlr1,qhlmlr12,qhlmlr23,qhlmlr3  !
         !     tmp = (6.*rho0(mgs)*(qhlmlr(mgs)-qhlmlr1)/(xdn0(lr)*tmpcmlt*pi) )**(1./3.)
         !     tmp1 = (6.*rho0(mgs)*(qhlmlr(mgs))/(xdn0(lr)*chlmlrr(mgs)*pi) )**(1./3.)
         !     write(6,*) 'diam0, diamave = ',tmp*1000.,tmp1*1000.
         !   ENDIF
         
!         IF ( igs(mgs) == 40 .and. qhlmlr(mgs) < -qxmin(lhl) ) THEN
!           write(iunit,*) 'k,qhl,alp,qhlmlr,chlmlrr(mgs),xdia ',kgs(mgs),1.e3*qx(mgs,lhl),alpha(mgs,lhl),qhlmlr(mgs),chlmlrr(mgs),xdia(mgs,lhl,3)
!           write(iunit,*) 'zx,tmpcmlt,qhlmlr2,qhlmlr12,xdn',zx(mgs,lhl),tmpcmlt,qhlmlr2,qhlmlr12,xdn(mgs,lhl),qhlmlr12*mltmass2inv,qhlmlr2*mltmass1inv
!         ENDIF

         if ( .true. .or. alpha(mgs,lhl) < alphamax - delta_alphamlr ) then ! { allow spectrum to get narrower, BUT still want to calculate chlmlrr
     !     chlmlr(mgs)  = qhlmlr0/(rho0(mgs)*qxd1/cxd1) 
          chlmlr(mgs)  =  tmpcmlt ! (qhlml0(ndiam+1)-qhlmlr05)/((qxd1-qxd05)/(cxd1-cxd05)) + qhlmlr05/(qxd05/cxd05)
          
          ! check for negative value of cx, and, if so revert to typical number loss rate.
          if ( cx(mgs,lhl) + dtp*chlmlr(mgs) < 0.0 ) then
             chlmlr(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
          endif
          
! test integrated Z rate for hail

        tmp = 1.0 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5.0 + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        x = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp

        tmp = 5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        y = (gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami)/g1palp
        
        hwvent1 =  0.78*x + y*hlventy(mgs)
       
       ! whole particle loss version: factor of 1/2 different from shedding version -- THIS IS THE COMPLETE GAMMA... 0 to infinity
       zhlmlr(mgs) =   &
     &   min(   rho0(mgs)* &
     &  fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1)**4   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  
         ! also note that fmlt1 includes a factor of 2*pi, resulting in a correct numerical factor of 12
         ! factor of rho0 to nullify the 1/rho0 in fmlt1

       zhlmlr0inf = zhlmlr(mgs)

      endif !}


!
! not sure this is right -- should it use g1x(mgs,lhl) or some assumed value for melt rain?
!
! should separately add the reflecitivities from the 2 sizes of shed drops, then add the completely melted particles

!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrr(mgs)  ) )


! test integrated Z rate for rain for size range up to mltdiam1
! integrate from d1 to infinity, then subtract that from the full value (zhlmlr)
        ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lhl,1) )


        tmp = 1. + alpha(mgs,lhl)
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

        tmp = 5. + alpha(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        if ( usegamxinf3 ) then
!          x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
          x = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
        else

!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
!           delx = ratio - float(i)*dqiacrratio
!           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
!           ip1 = Min( i+1, nqiacrratio )
!           jp1 = Min( j+1, nqiacralpha )
!
!           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,5,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,5,2) - gamxinflu(i,j,5,2))
!           tmp2 = gamxinflu(i,jp1,5,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,5,2) - gamxinflu(i,jp1,5,2))
!           
!           ! interpolate along alpha; 
!           
!!           tmp3 = (tmp1 + dely*(tmp2 - tmp1))*Gamma_sp(5. + alpha(mgs,lhl))
!           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1))
!           tmp5 = tmp3

           tmp3 = gaminterp(ratio,alpha(mgs,lhl),5,2)

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
            !  tmp1 = gaminterp(ratio,alpha(mgs,lhl),8,2)
              write(6,*) 'gamcheck11: i,k,tmp3,x = ',igs(mgs),kgs(mgs),tmp3,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            endif
          endif

!           IF ( tmp3 > 0.98 ) THEN
!            tmp3 = gamxinf(5. + alpha(mgs,lhl), ratio)/g1palp
!           ENDIF


!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - tmp3)/g1palp
!        x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)*(1. - tmp3)
   !     x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp - tmp3
           x = tmp3
           x2 = x

          if ( x <= 0.0 ) then ! check for unphysical negative result from interpolation
       !    x = 0.0
    !       x = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
           x = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
          else
    !       x2 = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
    !       x2 = (Gamma_sp(5. + alpha(mgs,lhl)) - gamxinf(5. + alpha(mgs,lhl), ratio))/g1palp
          endif

          if ( .false. .and. igs(mgs) == 40 ) then
            write(0,*) '3 x: tmp3,tmp5 = ',tmp3,tmp5
          endif

        endif

        tmp = 5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl)
        i1 = Int(dgami*(tmp))
        del = tmp - dgam*i1

        ! problem with interpolation of this one....
        if ( usegamxinf3 ) then
!           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
           y = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
        else
!           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
!           delx = ratio - float(i)*dqiacrratio
!           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
!           ip1 = Min( i+1, nqiacrratio )
!           jp1 = Min( j+1, nqiacralpha )
!
!           ! interpolate along x, i.e., ratio; 
!           tmp1 = gamxinflu(i,j,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,2) - gamxinflu(i,j,6,2))
!           tmp2 = gamxinflu(i,jp1,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,2) - gamxinflu(i,jp1,6,2))
!           
!           ! interpolate along alpha; 
!           tmp3 = (tmp1 + dely*(tmp2 - tmp1)*dqiacralphainv) ! *Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))
!
!           tmp4 = tmp3
           
!           IF ( igs(mgs) == 40 ) THEN
!             tmp3 = gaminterp(ratio,alpha(mgs,lhl),6,-2)
!           ELSE
             tmp3 = gaminterp(ratio,alpha(mgs,lhl),6,2)
!           ENDIF

          if ( gammacheck ) then
            if ( igs(mgs) == nx/2 ) then
              tmpgam = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - &
                        gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
            !  tmp1 = gaminterp(ratio,alpha(mgs,lhl),8,2)
              write(6,*) 'gamcheck12: i,k,tmp3,y = ',igs(mgs),kgs(mgs),tmp3,tmpgam
          !    tmp1 =  qx(mgs,lhl)*(1. - gamxinf(4. + alpha(mgs,lhl), ratio)/g4palp)
          !    write(6,*) 'qxd05,gam: ',qxd05,tmp1
            endif
          endif
           
      !     y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
           y = tmp3
           y2 = y

          if ( y <= 0.0 ) then ! check for unphysical negative result from interpolation
          ! y = 0.0
          ! y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
           y = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
          else
!           y2 = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
!           y2 = (Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))  - gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
          endif
           
        if ( .false. .and. igs(mgs) == nx/2 ) then
          write(0,*) '3 kgs,k,usegamxinf/2/3 = ',kgs(mgs),k,usegamxinf,usegamxinf2,usegamxinf3
          write(0,*) '3 y-interp, tmp = ',y,tmp
          tmp = (gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio))/g1palp
          write(0,*) '3 tmp3,tmp4 gamx = ',tmp3,tmp4,tmp
!          call gaminterpsub(tmp3,ratio,alpha(mgs,lhl),6,2)
          write(0,*) 'tmp3 from subroutine: ',tmp3
          write(0,*) 'i,j,jp1,delx,dely,gam1,gamjp1 = ',i,j,jp1,delx,dely,gamxinflu(i,j,6,2),gamxinflu(i,jp1,6,2) 
          write(0,*) '3 term1,2 = ', (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami),  &
                       gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)
          write(0,*) '3 term1,2/g1palp = ', (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp, &
                       gamxinf(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp
          write(0,*) '3 ratio0,ratio,xdia(1) = ',mltdiam1/xdia(mgs,lhl,1),ratio,xdia(mgs,lhl,1),alpha(mgs,lhl)
        
        endif
           
        endif
        if ( .false. .and. igs(mgs) == 38 .and. qhlmlr(mgs) < 0.0 ) then
           i = Min(nqiacrratio,Int(ratio*dqiacrratioinv))
!           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           if ( alp0flag ) then
           j = Int(Max(0.0,Min(15.,alpha(mgs,lhl)))*dqiacralphainv)
           else
           j = Int(Max(minalphalu,Min(maxalphalu,alpha(mgs,lhl)))*dqiacralphainv)
           endif
           delx = ratio - float(i)*dqiacrratio
           dely = alpha(mgs,lhl) - float(j)*dqiacralpha
           ip1 = Min( i+1, nqiacrratio )
           jp1 = Min( j+1, nqiacralpha )

           ! interpolate along x, i.e., ratio; 
           tmp1 = gamxinflu(i,j,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,j,6,2) - gamxinflu(i,j,6,2))
           tmp2 = gamxinflu(i,jp1,6,2) + delx*dqiacrratioinv*(gamxinflu(ip1,jp1,6,2) - gamxinflu(i,jp1,6,2))
           
           ! interpolate along alpha; 
           tmp3 = (tmp1 + dely*dqiacralphainv*(tmp2 - tmp1)) ! *Gamma_sp(5.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl))
           
          write(iunit,*) 'ibinhmlr = 1: i,k = ',igs(mgs),kgs(mgs)
          write(iunit,*) 'y,ratio,alpha,mltdiam1,xdia(mgs,lhl,1) = ',y,ratio,alpha(mgs,lhl),mltdiam1,xdia(mgs,lhl,1)
          write(iunit,*) 'full gamx = ',gamxinf(5.5 +0.5*bxx(mgs,lhl) + alpha(mgs,lhl), ratio),bxx(mgs,lhl),bx(lhl)
           y = (gmoi(i1) + (gmoi(i1+1) - gmoi(i1))*del*dgami)/g1palp  - tmp3 
          write(iunit,*) 'interpolated y = ', y,tmp1,tmp2,tmp3, &
               gamxinflu(i,j,6,2),gamxinflu(ip1,j,6,2),gamxinflu(i,jp1,6,2),gamxinflu(ip1,jp1,6,2)
           
        endif

        hwvent1 =  0.78*x + y*hlventy(mgs)
!        hwvent1 =    &
!     &  ( 0.78*x +    &
!     &    0.308*fvent(mgs)*y*(xdia(mgs,lhl,1)**(0.5 + 0.5*bxx(mgs,lhl)))*   &
!     &            Sqrt(axx(mgs,lh)*rhovt(mgs)) )
       
       zhlmlrr(mgs) =   &
     &   rho0(mgs)*min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1)**4   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lhl))/pi comes from removing the mass/volume factors; 

        if ( .false. .and. igs(mgs) == nx/2 ) then
         tmp3 = &
     &   rho0(mgs)*min(   &
     &  fmlt1(mgs)*cx(mgs,lhl)*(0.78*x2 + y2*hlventy(mgs))*xdia(mgs,lhl,1)**4   &
     &   , 0.0 )/xdn(mgs,lhl)*(6./pi)  ! factor of 6 comes from d/dt(D^6) = 6*D^5; 1/(xdn(mgs,lhl))/pi comes from removing the mass/volume factors; 

          write(0,*) '3 zhlmlrr(mgs),zhlmlrr2 = ',zhlmlrr(mgs),tmp3
          write(0,*) '3 x,x2,y,y2 = ',x,x2,y,y2,hwvent1,0.78*x + y2*hlventy(mgs)
        
        endif

! so now have zhmlrr as the Z whole-particle loss for hail for 0 to mltdiam1, and zhmlr0inf as the integrated whole-particle loss
! We want whole particle loss for 0 < d < mltdiam1, plus shedding loss for d > mltdiam1. Shedding loss is simple double the whole particle loss

       zhlmlrr(mgs) = Min( 0.0, zhlmlrr(mgs) )
       zhlmlr0inf   = Min( 0.0, zhlmlr0inf)
       zhlmlr(mgs) = zhlmlrr(mgs) + binmlrzrrfac*2.*(zhlmlr0inf - zhlmlrr(mgs)) ! factor of 2 to convert whole particle loss to shedding loss
       ! should that be zhlmlr(mgs) = 2.*zhlmlrr(mgs) + (zhlmlr0inf - zhlmlrr(mgs)) = zhlmlrr(mgs) + zhlmlr0inf
       ! yes, I think it should -- NOPE!
     !  zhlmlr(mgs) = 2.*zhlmlrr(mgs) + Min(0.0, (zhlmlr0inf - zhlmlrr(mgs)) )

        if ( .false. .and. igs(mgs) == nx/2 ) then

          write(0,'(a,i3,10(1pe18.10))') 'kz, zhlmlr0inf,z0-zrr,zhlmlrr,zhlmlr  = ',kgs(mgs),zhlmlr0inf, &
                     zhlmlr0inf - zhlmlrr(mgs),zhlmlrr(mgs),zhlmlr(mgs)
          write(0,'(a,10(1pe18.10))') 'old zhlmlr = ',zhlmlrr(mgs) + 2.*(zhlmlr0inf - zhlmlrr(mgs))
          write(0,'(a,10(1pe18.10))' ) 'test zhlmlr = ',2.*zhlmlrr(mgs) + Min(0.0, (zhlmlr0inf - zhlmlrr(mgs)) )
          write(0,'(a,10(1pe18.10))') 'new - old = ',zhlmlr(mgs) - ( zhlmlrr(mgs) + 2.*(zhlmlr0inf - zhlmlrr(mgs)) )
        !  write(0,*) '3 x,x2,y,y2 = ',x,x2,y,y2,hwvent1,0.78*x + y2*hlventy(mgs)
        
        endif

!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr0 - tmp**2 * chlmlr(mgs)  ) )

         tmp1 = (xdn(mgs,lhl)/1000.)**2 * zhlmlrr(mgs) 
         if ( numshedregimes == 2 ) then
          zhlmlrr(mgs) = (xdn(mgs,lhl)/1000.)**2 *zhlmlrr(mgs) + rho0(mgs)*  &
     &              (qhlmlr2*mltmass1inv*(0.01*takshedsize1)**6 + qhlmlr12*mltmass2inv*(0.01*takshedsize2)**6 )
         elseif ( numshedregimes == 3 ) then
!           chlmlrr(mgs) = tmpcmlt + rho0(mgs)*(qhlmlr3*mltmass1inv + qhlmlr23*mltmass2inv + qhlmlr12*mltmass3inv)
          zhlmlrr(mgs) = (xdn(mgs,lhl)/1000.)**2 *zhlmlrr(mgs) + rho0(mgs)*  &
     &              (qhlmlr3*mltmass1inv*(0.01*takshedsize1)**6 + qhlmlr23*mltmass2inv*(0.01*takshedsize2)**6  + &
     &               qhlmlr12*mltmass3inv*(0.01*takshedsize3)**6 )
!     &              (qhlmlr3*mltmass1inv*(0.01*0.5*takshedsize1)**6 + qhlmlr23*mltmass2inv*(0.01*0.5*takshedsize2)**6  + &
!     &               qhlmlr12*mltmass3inv*(0.01*0.5*takshedsize3)**6 )
         endif
         
         if ( .false. .and. igs(mgs) == nx/2 ) then
          write(0,'(a,i3,10(1pe18.10))') 'kz, zhlmlrr1,zhlmlrr,zhlmlr  = ',kgs(mgs),tmp1,zhlmlrr(mgs),zhlmlr(mgs)
          write(0,'(a,10(1pe18.10))') 'z1,z2,z3 = ', rho0(mgs)*qhlmlr3*mltmass1inv*(0.01*takshedsize1)**6 , &
     &               rho0(mgs)*qhlmlr23*mltmass2inv*(0.01*takshedsize2)**6 , &
     &               rho0(mgs)*qhlmlr12*mltmass3inv*(0.01*takshedsize3)**6 
!          write(0,'(a,10(1pe18.10))') 'z1,z2,z3 = ', rho0(mgs)*qhlmlr3*mltmass1inv*(0.01*0.5*takshedsize1)**6 , rho0(mgs)*qhlmlr23*mltmass2inv*(0.01*0.5*takshedsize2)**6 , &
!     &               rho0(mgs)*qhlmlr12*mltmass3inv*(0.01*0.5*takshedsize3)**6 
          write(0,*) 'chlmlrr = ',chlmlrr(mgs)
         endif

         ! chlmlr(mgs)  = Min( chlmlr(mgs), (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)) ! Make sure hail size does not go *down*

!          tmp = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
         if ( alpha(mgs,lhl) >= alphamax - delta_alphamlr ) then ! { cannot allow spectrum to get any narrower, revert to 2-moment method
         
          chlmlr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
!          chlmlrsave(mgs) = chlmlr(mgs)
!          qhlmlrsave(mgs) = qhlmlr(mgs)
!          chlsave(mgs) = cx(mgs,lhl)
!          qhlsave(mgs) = qx(mgs,lhl)
 
         

          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )

!        chmlrr(mgs) = Min( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lh)*xv(mgs,lh)) ) ! into rain 
! guess what, this is the same as chmlr: rho0*qhmlr/xmas(lh) --> cx/qx = rho0/xmas
!          IF(imltshddmr == 1 ) THEN
!            ! DTD: If Dmg < sheddiam, then assume complete melting into
!            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
!            tmp = -rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl))) ! Min of Maximum raindrop size/mean hail size
!            tmp2 = -rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
!            chlmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lh,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lh,3)-sheddiam)/(sheddiam0-sheddiam)
!            chlmlrr(mgs) = -Max(tmp,Min(tmp2,chlmlrr(mgs)))
!          ELSEIF ( imltshddmr == 2 .or. imltshddmr == 3 ) THEN
!            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
!            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
!            chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lhl))  ! into rain 
!          ENDIF


!          chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl)))  ! into rain 

!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrr(mgs)  ) )

         
         elseif (  xdia(mgs,lhl,3) > binmlrmxdia ) then ! 
         
          chlmlrtmp  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          zhlmlrtmp =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrtmp )
          
          tmp1 = 1.e-3
          if ( xdia(mgs,lhl,3) > binmlrmxdia + tmp1 ) then 
             chlmlr(mgs) = chlmlrtmp
             zhlmlr(mgs) = zhlmlrtmp
          else ! weighted average
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'k,xdia,chlmlr,zhlmlr = ',kgs(mgs),xdia(mgs,lhl,3)*1.e3,chlmlr(mgs),zhlmlr(mgs)
!              write(0,*) 'weights 1,2 = ',(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 , (binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
!            ENDIF
            chlmlr(mgs) = chlmlrtmp*(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 + chlmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
            zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )
            !zhlmlr(mgs) = zhlmlrtmp*(xdia(mgs,lhl,3) - binmlrmxdia )/tmp1 + zhlmlr(mgs)*(binmlrmxdia + tmp1 - xdia(mgs,lhl,3))/tmp1
!            IF ( igs(mgs) == 3 ) THEN
!              write(0,*) 'tmp,chlmlr,= ',chlmlrtmp,chlmlr(mgs)
!              write(0,*) 'tmp,zhlmlr,= ',zhlmlrtmp,zhlmlr(mgs)
!            ENDIF
          endif

          ! check for negative value of cx, and, if so revert to typical number loss rate.
          if ( cx(mgs,lhl) + dtp*chlmlr(mgs) < 0.0 ) then
             chlmlr(mgs) = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
             zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )
          endif
          
         
         endif !} 


       elseif ( ibinhlmlr == -1 ) then ! OLD VERSION use incomplete gamma functions to approximate the bin results

        endif ! ibinhlmlr


       if ( ivhmltsoak > 0 .and.  qhlmlr(mgs) < 0.0 .and. lvol(lhl) > 1 .and. xdn(mgs,lhl) .lt. xdnmx(lhl) ) then
         ! act as if 50% of the meltwater were soaked into the graupel
           v1 = (1. - xdn(mgs,lhl)/xdnmx(lhl))*(vx(mgs,lhl) + rho0(mgs)*qhlmlr(mgs)/xdn(mgs,lhl) )/(dtp) ! volume available for filling
           v2 = -1.0*rho0(mgs)*qhlmlr(mgs)/xdnmx(lhl)  ! volume of melted ice if it were refrozen in the matrix
           
           vhlsoak(mgs) = Min(v1,v2)
           
       endif
        
        endif
       endif

      endif
      
!
!      qimlr(mgs)  = max( qimlr(mgs), -qimxd(mgs) ) 
!      qsmlr(mgs)  = max( qsmlr(mgs),  -qsmxd(mgs) ) 
! erm 5/10/2007 changed to next line:
      if ( .not. mixedphase ) qsmlr(mgs)  = max( qsmlr(mgs),  Min( -qsmxd(mgs), -0.7*qx(mgs,ls)*dtpinv ) ) 
      if ( .not. mixedphase ) then
        qhmlr(mgs)  = max( qhmlr(mgs),  Min( -qhmxd(mgs), -0.95*qx(mgs,lh)*dtpinv ) ) 
        chmlr(mgs)  = max( chmlr(mgs),  Min( -chmxd(mgs), -0.95*cx(mgs,lh)*dtpinv ) ) 
      endif
!      qhmlr(mgs)  = max( max( qhmlr(mgs),  -qhmxd(mgs) ) , -0.5*qx(mgs,lh)*dtpinv ) !limits to 1/2 qh or max depletion
      qhmlh(mgs)  = 0. ! not used


      ! Rasmussen and Heymsfield say melt water remains on graupel up to 9 mm before shedding


      if ( lhl .gt. 1 .and. lhlw < 1 ) then
        qhlmlr(mgs)  = max( qhlmlr(mgs),  Min( -qxmxd(mgs,lhl), -0.95*qx(mgs,lhl)*dtpinv ) )
        chlmlr(mgs)  = max( chlmlr(mgs),  Min( -cxmxd(mgs,lhl), -0.95*cx(mgs,lhl)*dtpinv ) )
      endif

!
      end do

      endif  ! } not mixedphase
!
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      cimlr(mgs)  = (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qimlr(mgs)
      if ( .not. mixedphase ) then !{
        if ( xdia(mgs,ls,1) .gt. 1.e-6 .and. -qsmlr(mgs) .ge. 0.5*qxmin(ls) .and. ipconc .ge. 4 ) then 
!         csmlr(mgs)  = rho0(mgs)*qsmlr(mgs)/(xv(mgs,ls)*rhosm)
         csmlr(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)))*qsmlr(mgs)
        elseif ( qx(mgs,ls) > qxmin(ls) ) then
         csmlr(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)))*qsmlr(mgs)
        endif
        
        csmlrr(mgs) = csmlr(mgs)/rzxs(mgs)
         if ( -csmlrr(mgs)*dtp > cxmin .and. -qsmlr(mgs)*dtp > qxmin(lr) .and. snowmeltdia > 0.0 ) then
           rmas = rho0(mgs)*qsmlr(mgs)/csmlrr(mgs)
           if ( rmas > snowmeltmass ) then
             csmlrr(mgs) = rho0(mgs)*qsmlr(mgs)/snowmeltmass
           endif
         endif
           


!        IF ( xdia(mgs,lh,1) .gt. 1.e-6 .and. Abs(qhmlr(mgs)) .ge. qxmin(lh) ) THEN
!          chmlr(mgs) = rho0(mgs)*qhmlr(mgs)/(pi*xdn(mgs,lh)*xdia(mgs,lh,1)**3)  ! out of hail
!          chmlr(mgs) = Max( chmlr(mgs), -chmxd(mgs) )
!        ELSE
         if ( ibinhmlr == 0 .or. lzh < 1 ) then
           chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
           if ( imltshddmr == 3 .and. qhmlr(mgs) < -qxmin(lh) ) then
            !  tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            !  
            !  IF ( tmpdiam > sheddiam ) THEN ! let size get smaller until it reaches sheddiam
            !   chmlr(mgs) = 0.0
            !  ENDIF
            
            ! test to remove the part of the melting associated with large ice particles so they get smaller

            tmp = 1. + alpha(mgs,lh)
            i = Int(dgami*(tmp))
            del = tmp - dgam*i
            g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

            ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lh,1) )

            x =  gamxinfdp(2. + alpha(mgs,lh), ratio)/g1palp
            y =  gamxinfdp(2.5 + alpha(mgs,lh) + 0.5*bxx(mgs,lh), ratio)/g1palp

            hwvent1 =  0.78*x + y*hwventy(mgs) 

            qhlmlr1 = min( fmlt1(mgs)*cx(mgs,lh)*hwvent1*xdia(mgs,lh,1), 0.0 )

            chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*(qhmlr(mgs) - qhlmlr1)
           
           
           endif
!           IF ( igs(mgs) == 40 ) THEN
!             write(0,*) 'is this running? chmlr = ',kgs(mgs), chmlr(mgs)
!           ENDIF
         endif
!        ENDIF


     if ( chmlr(mgs) < 0.0 .and. (ibinhmlr < 1 .or. lzh < 1) ) then ! { already done if ibinhmlr > 0
      if (  ipconc >= 6 .and. lzr .gt. 1 .and. lzh < 1  .and. qx(mgs,lh) > qxmin(lh) ) then ! Only compute if rain is 3-moment but graupel is not, otherwise is computed later
          tmp = qx(mgs,lh)/cx(mgs,lh)
          alp = alpha(mgs,lh)
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
        
        zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )

      endif
      
      if ( ibinhmlr == 0 .or. lzh < 1 ) then
      if ( ihmlt .eq. 1 ) then
        chmlrr(mgs)  = Min( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
      elseif ( ihmlt .eq. 2 ) then
        if ( xv(mgs,lh) .gt. 0.0 .and. chmlr(mgs) .lt. 0.0 ) then
!        chmlrr(mgs) = Min( chmlr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lh)*xv(mgs,lh)) ) ! into rain 
! guess what, this is the same as chmlr: rho0*qhmlr/xmas(lh) --> cx/qx = rho0/xmas
          if(imltshddmr == 1) then
            ! DTD: If Dmg < sheddiam, then assume complete melting into
            ! maximal raindrop.  Between sheddiam and sheddiam0 mm, linearly ramp down to a 3 mm shed drop
            tmp = -rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh))) ! Min of Maximum raindrop size/mean hail size
            tmp2 = -rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
            
            chmlrr(mgs) = tmp*(sheddiam0-xdia(mgs,lh,3))/(sheddiam0-sheddiam)+tmp2*(xdia(mgs,lh,3)-sheddiam)/(sheddiam0-sheddiam)  ! old version
            chmlrr(mgs) = -Max(tmp,Min(tmp2,chmlrr(mgs)))
          elseif ( imltshddmr == 2 .or. imltshddmr == 3 ) then
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lh))  ! into rain 
          else ! Old method
            chmlrr(mgs) =  rho0(mgs)*qhmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lh)*xv(mgs,lh)))  ! into rain 
         endif
        else
        chmlrr(mgs) = chmlr(mgs)
        endif
      elseif ( ihmlt .eq. 0 ) then
        chmlrr(mgs) = chmlr(mgs)
      endif

      else ! ibinhmlr < 0? Already have an outer if test for ibinhmlr < 1
        chmlrr(mgs)  = Min( chmlrr(mgs), rho0(mgs)*qhmlr(mgs)/(xdn(mgs,lr)*xvmx(lr)) ) ! into rain 
      endif
      
      endif ! } ( chmlr(mgs) < 0.0 .and. ibinhmlr < 1)

      if ( lhl .gt. 1 .and. lhlw < 1 .and. .not. mixedphase .and. qhlmlr(mgs) < 0.0 ) then ! {
      
      if ( ibinhlmlr == 0 .or. lzhl < 1 ) then
!      IF ( xdia(mgs,lhl,1) .gt. 1.e-6 .and. Abs(qhlmlr(mgs)) .ge. qxmin(lhl) ) THEN
!      chlmlr(mgs) = rho0(mgs)*qhlmlr(mgs)/(pi*xdn(mgs,lhl)*xdia(mgs,lhl,1)**3)  ! out of hail
!      chlmlr(mgs) = Max( chlmlr(mgs), -cxmxd(mgs,lhl) )
!      ELSE
      chlmlr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlmlr(mgs)
           if ( imltshddmr == 3 .and. qhlmlr(mgs) < -qxmin(lhl) ) then
!           IF ( .false. .and. imltshddmr == 3  ) THEN
!              tmpdiam = (shedalp+alpha(mgs,lhl))*xdia(mgs,lhl,1)
!              
!              IF ( tmpdiam > sheddiam ) THEN ! let size get smaller until it reaches sheddiam
!                chlmlr(mgs) = 0.0
!              ENDIF
 
            ! test to remove the part of the melting associated with large ice particles so they get smaller
!
            tmp = 1. + alpha(mgs,lhl)
            i = Int(dgami*(tmp))
            del = tmp - dgam*i
            g1palp = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

            ratio = Min( maxratiolu,  mltdiam1/xdia(mgs,lhl,1) )

            x =  gamxinfdp(2. + alpha(mgs,lhl), ratio)/g1palp
            y =  gamxinfdp(2.5 + alpha(mgs,lhl) + 0.5*bxx(mgs,lhl), ratio)/g1palp

            hwvent1 =  0.78*x + y*hlventy(mgs) 

            qhlmlr1 = min( fmlt1(mgs)*cx(mgs,lhl)*hwvent1*xdia(mgs,lhl,1), 0.0 )

            chlmlr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*Min(0.0, qhlmlr(mgs) - qhlmlr1)

            endif
!      ENDIF
      endif
      
      if ( ibinhlmlr == 0 .or. lzhl < 1 ) then !{
      if ( ihmlt .eq. 1 ) then
        chlmlrr(mgs)  = Min( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vmlt) ) ! into rain 
      elseif ( ihmlt .eq. 2 ) then
        if ( xv(mgs,lhl) .gt. 0.0 .and. chlmlr(mgs) .lt. 0.0 ) then
!        chlmlrr(mgs) = rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl)))  ! into rain 
!        chlmlrr(mgs) = Min( chlmlr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lhl)*xv(mgs,lhl)) ) ! into rain 
          if(imltshddmr == 1 ) then
            tmp = -rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl))) ! Min of Maximum raindrop size/mean hail size
            tmp2 = -rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vr3mm) ! conc. change for a 3 mm mean drop diameter
            chlmlrr(mgs) = tmp*(20.e-3-xdia(mgs,lhl,3))/(20.e-3-sheddiam)+tmp2*(xdia(mgs,lhl,3)-sheddiam)/(20.e-3-sheddiam)
            chlmlrr(mgs) = -Max(tmp,Min(tmp2,chlmlrr(mgs)))
          elseif ( imltshddmr == 2 .or. imltshddmr == 3 ) then
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chlmlrr(mgs) =  rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lhl))  ! into rain 
          else ! old method
            chlmlrr(mgs) = rho0(mgs)*qhlmlr(mgs)/(Min(xdn(mgs,lr)*xvmx(lr), xdn(mgs,lhl)*xv(mgs,lhl)))  ! into rain 
          endif
        else
        chlmlrr(mgs) = chlmlr(mgs)
        endif
      elseif ( ihmlt .eq. 0 ) then
        chlmlrr(mgs) = chlmlr(mgs)
      endif

      else ! } { ibinhlmlr > 0
        chlmlrr(mgs)  = Min( chlmlrr(mgs), rho0(mgs)*qhlmlr(mgs)/(xdn(mgs,lr)*xvmx(lr)) ) ! into rain 
      endif !}
      
        
       if ( ipconc >= 8 .and. lzhl .gt. 1 .and. ibinhlmlr <= 0 ) then
        if ( cx(mgs,lhl) > 0.0 ) then

          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          alp = alpha(mgs,lhl)
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
        
        zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( tmp * qhlmlr(mgs) )
       endif
      endif
      endif ! }

      endif ! }.not. mixedphase 

! 10ice versions:
!      chmlr(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhmlr(mgs)
!      chmlrr(mgs) = chmlr(mgs)
      end do
      end if

!
!  deposition/sublimation of ice
!
      do mgs = 1,ngscnt

      rwcap(mgs) = (0.5)*xdia(mgs,lr,1)
      swcap(mgs) = (0.5)*xdia(mgs,ls,1)
      hwcap(mgs) = (0.5)*xdia(mgs,lh,1)
      if ( lhl .gt. 1 ) hlcap(mgs) = (0.5)*xdia(mgs,lhl,1)

      if ( qx(mgs,li).gt.qxmin(li) .and. xdia(mgs,li,1) .gt. 0.0 ) then
!
! from Cotton, 1972 (Part II)
!
        cilen(mgs)   = 0.4764*(xdia(mgs,li,1))**(0.958)
        cval = xdia(mgs,li,1)
        aval = cilen(mgs)
        eval = Sqrt(1.0-(aval**2)/(cval**2))
        fval = min(0.99,eval)
        gval = alog( abs( (1.+fval)/(1.-fval) ) )
        cicap(mgs) = cval*fval / gval
      else
       cicap(mgs) = 0.0
      end if
      enddo
!
!
      qhdsv(:) = 0.0
      qhldsv(:) = 0.0

      do mgs = 1,ngscnt
      if ( icond .eq. 1 .or. temg(mgs) .le. tfrh    &
     &      .or. (qx(mgs,lr) .le. qxmin(lr) .and. qx(mgs,lc) .le. qxmin(lc)) ) then
        qidsv(mgs) =   &
     &    fvds(mgs)*cx(mgs,li)*civent(mgs)*cicap(mgs)*depfac
        qsdsv(mgs) =   &
     &    fvds(mgs)*cx(mgs,ls)*swvent(mgs)*swcap(mgs)*depfac

!        IF ( ny .eq. 2 .and. igs(mgs) .eq. 302 .and. temg(mgs) .le. tfrh+10 .and. qx(mgs,lv) .gt. qis(mgs)
!     :       .and. qx(mgs,li) .gt. qxmin(li) ) THEN
!         write(0,*) 'qidsv = ',nstep,kgs(mgs),qidsv(mgs),temg(mgs)-tfrh,100.*(qx(mgs,lv)/qis(mgs) - 1.),1.e6*xdia(mgs,li,1),
!     :            fvds(mgs),civent(mgs),cicap(mgs)
!        ENDIF
      else
        qidsv(mgs) = 0.0
        qsdsv(mgs) = 0.0
      endif
        qhdsv(mgs) =   &
     &    fvds(mgs)*cx(mgs,lh)*hwvent(mgs)*hwcap(mgs)*depfac

        if ( lhl .gt. 1 ) qhldsv(mgs) = fvds(mgs)*cx(mgs,lhl)*hlvent(mgs)*hlcap(mgs)*depfac
!
!
      end do
!


! #include "nssl.qlimit.F"

!
!  Use a test saturation adjustment to set limits on ice deposition/sublimation
!  and rain evaporation
!
!
      if ( DoSublimationFix ) then
      
      do mgs = 1,ngscnt

        qitmp(mgs) = qx(mgs,li) + qx(mgs,ls) + qx(mgs,lh)
        if ( lis > 1 ) qitmp(mgs) = qitmp(mgs) + qx(mgs,lis)
        if ( lhl > 1 ) qitmp(mgs) = qitmp(mgs) + qx(mgs,lhl)
        qrtmp(mgs) = qx(mgs,lr)
        qctmp(mgs) = qx(mgs,lc)
        qsimxdep(mgs) = 0.0
        qsimxsub(mgs) = 0.0
        dqcitmp(mgs) = 0.0
        

!      IF ( ( qitmp(mgs) > qxmin(li) .or. qrtmp(mgs) > qxmin(lr) ) ) THEN
      if ( qitmp(mgs) > qxmin(li)  ) then
      
        qitmp1    = qitmp(mgs)
        qctmp1    = qctmp(mgs)
        felvcptmp = felvcp(mgs)
        felscptmp = felscp(mgs)
        qvtmp(mgs) = qx(mgs,lv)
        qss(mgs) = qvs(mgs)
        qsstmp = qvs(mgs)
        qvstmp = qvs(mgs)
        qisstmp = qis(mgs)
        thetatmp  = theta(mgs)
        thetaptmp = thetap(mgs)
        temgtmp   = temg(mgs)
        temcgtmp  = temcg(mgs)
        qvaptmp   = qx(mgs,lv) ! qwvp(mgs) + qv0(mgs)
        qvptmp    = 0.0 ! qwvp(mgs)  ! qv pertubation

        qsstmp = qisstmp

      
       dqwvtmp(mgs) = ( qvtmp(mgs) - qsstmp )

      do itertd = 1,2
      
!
!  calculate super-saturation
!
      if ( itertd == 1 ) then
      
      else
        dqcitmp(mgs) = dqci(mgs)
   !     dqwvtmp(mgs) = dqwv(mgs)
      endif

      dqcw(mgs) = 0.0
      dqci(mgs) = 0.0
      dqwv(mgs) = ( qvtmp(mgs) - qsstmp )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then           ! { subsaturated
        if( qitmp(mgs) .gt. -dqwv(mgs) ) then  ! check if qi can make up all the deficit
          dqci(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                  ! otherwise make all ice available for sublimation
          dqci(mgs) = -qitmp(mgs)
          dqwv(mgs) = dqwv(mgs) + qitmp(mgs)
        end if
!
       qvptmp = qvptmp - ( dqcw(mgs) + dqci(mgs) )  ! add to perturbation vapor

       if ( itertd == 2 .and. eqtset > 1 ) then
       ! if eqtset == 2, then need to update the latent heats for change in hydrometeor content
          tmp = qitmp(mgs) !+ qx(mgs,lh)
!          IF ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          cvm = cv+cvv*qvtmp(mgs)+cpl*(qx(mgs,lc)+qrtmp(mgs))   &
                                  +cpigb*(tmp)

          felvcptmp = (felv(mgs)-rw*temg(mgs))/cvm
          felscptmp = (fels(mgs)-rw*temg(mgs))/cvm
       endif


!      qitmp(mgs) = qx(mgs,li)
      qctmp(mgs) = qctmp(mgs) + dqcw(mgs) ! dqcw is zero
      qitmp(mgs) = qitmp(mgs) + dqci(mgs)
      thetaptmp = thetaptmp +   &
     &  1./pi0(mgs)*   &
     &  (felvcp(mgs)*dqcw(mgs) +felscp(mgs)*dqci(mgs))


      end if  ! } dqwv(mgs) .lt. 0. (end of evap/sublim)
!
! condensation/deposition
!
      if ( dqwv(mgs) .ge. 0. ) then ! {
      
!      write(iunit,*) 'satadj: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qx(mgs,lv),qx(mgs,lc)
!
!        qitmp(mgs) = qx(mgs,li)
        fracl(mgs) = 0.0
        fraci(mgs) = 1.0
        if ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) then
!          fracl(mgs) = max(min(1.,(temg(mgs)-233.15)/(20.)),0.0)
!          fraci(mgs) = 1.0-fracl(mgs)
        end if
        if ( temg(mgs) .le. thnuc ) then
           fraci(mgs) = 1.0
           fracl(mgs) = 0.0
         end if
!        fraci(mgs) = 1.0-fracl(mgs)

       gamss = (felvcp(mgs)*fracl(mgs) + felscp(mgs)*fraci(mgs))   &
     &      / (pi0(mgs))

          dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv2(mgs)*qsstmp/   &
     &  ((temg(mgs)-cbi)**2))

      if ( temg(mgs) .ge. tfr ) then
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qsstmp/   &
     &  ((temg(mgs)-cbw)**2))
      end if

      delqci1=qx(mgs,li)


      dqcw(mgs) = dqvcnd(mgs)*fracl(mgs) ! is zero
      dqci(mgs) = dqvcnd(mgs)*fraci(mgs)

      thetaptmp = thetaptmp +   &
     &   (felvcp(mgs)*dqcw(mgs) + felscp(mgs)*dqci(mgs))   &
     & / (pi0(mgs))

      qvptmp = qvptmp - ( dqvcnd(mgs) )
      qctmp(mgs) = qctmp(mgs) + dqcw(mgs)
      qitmp(mgs) = qitmp(mgs) + dqci(mgs)

       if ( itertd == 2 .and. eqtset > 1 ) then
       ! if eqtset == 2, then need to update the latent heats for change in hydrometeor content
          tmp = qitmp(mgs) ! + qx(mgs,lh)
!          IF ( lhl > 1 ) tmp = tmp + qx(mgs,lhl)
          cvm = cv+cvv*qvtmp(mgs)+cpl*(qctmp(mgs) +qrtmp(mgs))   &
                                  +cpigb*(tmp)

          felvcptmp = (felv(mgs)-rw*temg(mgs))/cvm
          felscptmp = (fels(mgs)-rw*temg(mgs))/cvm
       endif

      if ( eqtset > 2 ) then
        pipert(mgs) = pipert(mgs) + (0   &
     &  +felspi(mgs)*dqci(mgs)    &
     &  +felvpi(mgs)*dqcw(mgs))*dtp
      endif

!
!
      end if ! } dqwv(mgs) .ge. 0.


!
      if ( itertd == 1 ) then
      ! update temporary saturation values

      thetatmp = thetaptmp + theta0(mgs)
      temgtmp = thetatmp*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      qvaptmp = Max((qvptmp + qv0(mgs)), 0.0)
      temcgtmp = temgtmp - tfr
      tqvcon = temgtmp-cbw
      ltemq = (temgtmp-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvstmp = pqs(mgs)*tabqvs(ltemq)
      qisstmp = pqs(mgs)*tabqis(ltemq)
      qctmp(mgs) = max( 0.0, qctmp(mgs) )
      qitmp(mgs) = max( 0.0, qitmp(mgs) )
      qvtmp(mgs) = max( 0.0, qvaptmp )
      
!      qsstmp = qvstmp
      qsstmp = qisstmp
      
      else
       ! set max depletion
        qctmp(mgs) = max( 0.0, qctmp(mgs) )
        qitmp(mgs) = max( 0.0, qitmp(mgs) )
       
        if ( qitmp(mgs) < qitmp1 ) then
          qsimxsub(mgs) = (qitmp1 - qitmp(mgs))*dtpinv
        elseif ( qitmp(mgs) > qitmp1 ) then
          qsimxdep(mgs) = (qitmp(mgs) - qitmp1)*dtpinv
        endif
       
      
      endif
!      pceds(mgs) = (thetap(mgs) - thsave(mgs))*dtpinv
!      write(iunit,*) 'satadj2: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qxtmp,qctmp(mgs)
!
!  end the saturation adjustment iteration loop
!
      end do ! itertd
      
      endif
      
      end do ! mgs
      
      else
      
       do mgs = 1,ngscnt
         qsimxdep(mgs) = qvimxd(mgs)
         qsimxsub(mgs) = 1.e20
       enddo
      
      endif

! end of qlimit

      qhcev(:) = 0.0
      chcev(:) = 0.0
      qhlcev(:) = 0.0
      chlcev(:) = 0.0
      qfcev(:) = 0.0

      do mgs = 1,ngscnt
      qisbv(mgs) = 0.0
      qssbv(mgs) = 0.0
      qidpv(mgs) = 0.0
      qsdpv(mgs) = 0.0
      qhsbv(mgs) = 0.0
      qscev(mgs) = 0.0
      cscev(mgs) = 0.0
      if ( icond .eq. 1 .or. temg(mgs) .le. tfrh    &
     &      .or. (qx(mgs,lr) .le. qxmin(lr) .and. qx(mgs,lc) .le. qxmin(lc)) ) then ! last condition (qr<qmin & qc<qmin) for case icond=0
!        qisbv(mgs) = max( min(qidsv(mgs), 0.0), -qimxd(mgs) )
!        qssbv(mgs) = max( min(qsdsv(mgs), 0.0), -qsmxd(mgs) )
! erm 5/10/2007:
        qisbv(mgs) = max( min(qidsv(mgs), 0.0),  Min( -qimxd(mgs), -0.5*qx(mgs,li)*dtpinv ) )
        if ( temg(mgs) < tfr .or. .not. qsmlr(mgs) < 0.0 ) then
        qssbv(mgs) = max( min(qsdsv(mgs), 0.0),  Min( -qsmxd(mgs), -0.5*qx(mgs,ls)*dtpinv ) )
        endif
        qidpv(mgs) = Max(qidsv(mgs), 0.0)
        qsdpv(mgs) = Max(qsdsv(mgs), 0.0)
        
        if ( qsmlr(mgs) < 0.0 .and. .not. mixedphase ) then ! switch snow sublimation to evaporation if there is melting

          qscev(mgs) = evapfac*   &
     &  4.0*pi*(qx(mgs,lv)-qss0(mgs))*cx(mgs,ls)*swcap(mgs)*swvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))
          qscev(mgs) = Max( Min(0.0,qscev(mgs)),  Min( -qsmxd(mgs), -0.5*qx(mgs,ls)*dtpinv ) )
        else

        endif



      else
        qisbv(mgs) = 0.0
        qssbv(mgs) = 0.0
        qidpv(mgs) = 0.0
        qsdpv(mgs) = 0.0
      endif

      qhsbv(mgs) = 0.0
      qhdpv(mgs) = 0.0
      if ( qx(mgs,lh) > qxmin(lh) ) then
      if ( temg(mgs) < tfr .or. .not. qhmlr(mgs) < 0.0 ) then
      ! no liquid from melting, so evaporation is greater. Thus can calculate sublimation rate
      qhsbv(mgs) = max( min(qhdsv(mgs), 0.0), -qhmxd(mgs) )
      qhdpv(mgs) = Max(qhdsv(mgs), 0.0)
      endif
      
      if ( .true. .and. qhmlr(mgs) < 0.0 .and. .not. mixedphase ) then
        ! Liquid is forming, so find the evaporation that was subtracted from melting (if it is not condensing)
!       qhcev(mgs) =   &
!     &   evapfac*min(   &
!     &  fmlt1e(mgs)*cx(mgs,lh)*hwvent(mgs)*xdia(mgs,lh,1), 0.0 )
        
      qhcev(mgs) =  evapfac*  &
     &  2.0*pi*(qx(mgs,lv)-qss0(mgs))*cx(mgs,lh)*xdia(mgs,lh,1)*hwvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))

      qhcev(mgs)  = max(qhcev(mgs), -qhmxd(mgs))
        
      endif
      endif


      qhlsbv(mgs) = 0.0
      qhldpv(mgs) = 0.0
      if ( lhl .gt. 1 ) then
      if ( qx(mgs,lhl) > qxmin(lhl) ) then
        if ( temg(mgs) < tfr .or. .not. qhlmlr(mgs) < 0.0 ) then
        qhlsbv(mgs) = max( min(qhldsv(mgs), 0.0), -qxmxd(mgs,lhl) )
        qhldpv(mgs) = Max(qhldsv(mgs), 0.0)
        endif
        if ( qhlmlr(mgs) < 0.0 .and. .not. mixedphase ) then
        ! Liquid is forming, so find the evaporation that was subtracted from melting (if it is not condensing)
         qhlcev(mgs) =  evapfac*  &
     &      2.0*pi*(qx(mgs,lv)-qss0(mgs))*cx(mgs,lhl)*xdia(mgs,lhl,1)*hlvent(mgs)/(qss0(mgs)*(fav(mgs)+fbv(mgs)))

         qhlcev(mgs)  = max(qhlcev(mgs), -qhlmxd(mgs))
        
      endif
      endif
      endif
      
      temp1 = qidpv(mgs) + qsdpv(mgs) + qhdpv(mgs) + qhldpv(mgs)

!      IF ( temp1 .gt. qvimxd(mgs) ) THEN

!      frac = qvimxd(mgs)/temp1

      if ( temp1 .gt. qsimxdep(mgs) ) then
      frac = qsimxdep(mgs)/temp1

      qidpv(mgs) = frac*qidpv(mgs)
      qsdpv(mgs) = frac*qsdpv(mgs)
      qhdpv(mgs) = frac*qhdpv(mgs)
      qhldpv(mgs) = frac*qhldpv(mgs)

!        IF ( ny .eq. 2 .and. igs(mgs) .eq. 302 .and. temg(mgs) .le. tfrh+10 .and. qx(mgs,lv) .gt. qis(mgs)
!     :       .and. qx(mgs,li) .gt. qxmin(li) ) THEN
!         write(0,*) 'qidpv,frac = ',kgs(mgs),qidpv(mgs),frac
!        ENDIF

      endif

      temp1 = qisbv(mgs) + qssbv(mgs) + qhsbv(mgs) + qhlsbv(mgs)


      if ( temp1 <  -qsimxsub(mgs) ) then
      frac = -qsimxsub(mgs)/temp1

      qisbv(mgs) = frac*qisbv(mgs)
      qssbv(mgs) = frac*qssbv(mgs)
      qhsbv(mgs) = frac*qhsbv(mgs)
      qhlsbv(mgs) = frac*qhlsbv(mgs)

!        IF ( ny .eq. 2 .and. igs(mgs) .eq. 302 .and. temg(mgs) .le. tfrh+10 .and. qx(mgs,lv) .gt. qis(mgs)
!     :       .and. qx(mgs,li) .gt. qxmin(li) ) THEN
!         write(0,*) 'qidpv,frac = ',kgs(mgs),qidpv(mgs),frac
!        ENDIF

      endif


      end do
!
!
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      cssbv(mgs)  = (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qssbv(mgs)
      cisbv(mgs)  = (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qisbv(mgs)
      chsbv(mgs)  = (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhsbv(mgs)
      if ( lhl .gt. 1 ) chlsbv(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlsbv(mgs)
      csdpv(mgs)  = 0.0 ! (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*qsdpv(mgs)
      cidpv(mgs) =  0.0 ! (cx(mgs,li)/(qx(mgs,li)+1.e-20))*qidpv(mgs)
      cisdpv(mgs) = 0.0
      chdpv(mgs)  = 0.0 ! (cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhdpv(mgs)
      chldpv(mgs) = 0.0
      end do
      end if

!
!  Aggregation or size conversion of small crystals to snow
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 29a'
      do mgs = 1,ngscnt
      qscni(mgs) =  0.0
      cscni(mgs) = 0.0
      cscnis(mgs) = 0.0
      if ( ipconc .ge. 4 .and. iscni .ge. 1 .and. qx(mgs,li) .gt. qxmin(li) ) then
        if ( iscni .eq. 1 ) then
         qscni(mgs) =    &
     &      pi*rho0(mgs)*((0.25)/(6.0))   &
     &      *eii(mgs)*(qx(mgs,li)**2)*(xdia(mgs,li,2))   &
     &      *vtxbar(mgs,li,1)/xmas(mgs,li)
         cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li))
         cscnis(mgs) = 0.5*cscni(mgs)
        elseif ( iscni .eq. 2 .or. iscni .eq. 4 .or. iscni .eq. 5 ) then  ! Zeigler 1985/Zrnic 1993, sort of
          if ( iscni .ne. 5 .and. qidpv(mgs) .gt. 0.0 .and.  xdia(mgs,li,3) .ge. 100.e-6 ) then
          ! convert larger crystals to snow
!            IF ( xdia(mgs,ls,3) .gt. xdia(mgs,li,3) ) THEN
!              qscni(mgs) = Max(0.1,xdia(mgs,li,3)/xdia(mgs,ls,3))*qidpv(mgs)
! erm 9/5/08 changed max to min
              qscni(mgs) = Min(0.5, xdia(mgs,li,3)/200.e-6)*qidpv(mgs)
!            ELSE
!              qscni(mgs) = 0.1*qidpv(mgs)
!            ENDIF
            cscni(mgs) = fscni*qscni(mgs)*rho0(mgs)/Max(rho_qs*xvmn(ls),xmas(mgs,li))
!            cscni(mgs) = fscni*Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/Max(xdn(mgs,ls)*xvmn(ls),xmas(mgs,li)))
!            cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li) )
!            IF ( xdia(mgs,ls,3) .le. 200.e-6 ) THEN
              cscnis(mgs) = cscni(mgs)
!            ELSE
!              cscnis(mgs) = 0.0
!            ENDIF
          endif

           if ( iscni .ne. 4 ) then
           ! crystal aggregation to become snow
! erm 9/5/08 commented second line and added xv to 1st line (zrnic et al 1993)
             tmp = ess(mgs)*rvt*aa2*cx(mgs,li)*cx(mgs,li)*xv(mgs,li)
!     :         ((cinu + 2.)*xv(mgs,li)/(cinu + 1.) + xv(mgs,li))

!           csacs(mgs) = rvt*aa2*ess(mgs)*cx(mgs,ls)**2*xv(mgs,ls)

             qscni(mgs) = qscni(mgs) + Min( qxmxd(mgs,li), 2.0*tmp*xmas(mgs,li)*rhoinv(mgs) )
             cscni(mgs) = cscni(mgs) + Min( cxmxd(mgs,li), 2.0*tmp )
             cscnis(mgs) = cscnis(mgs) + Min( cxmxd(mgs,li), tmp )
           endif
        elseif ( iscni .eq. 3 ) then ! LFO
           qscni(mgs) = 0.001*eii(mgs)*max((qx(mgs,li)-1.e-3),0.0)
           qscni(mgs) = min(qscni(mgs),qxmxd(mgs,li))
           cscni(mgs) = qscni(mgs)*rho0(mgs)/xmas(mgs,li)
           cscnis(mgs) = 0.5*cscni(mgs)
!           write(iunit,*) 'qscni, qi = ',qscni(mgs),qx(mgs,li),igs(mgs),kgs(mgs)
        endif

      elseif ( ipconc < 4 ) then ! LFO
           if ( lwsm6 ) then
             qimax = rhoinv(mgs)*roqimax
             qscni(mgs) = Min(0.90*qx(mgs,li), Max( 0.0, (qx(mgs,li) - qimax)*dtpinv ) )
           else
             qscni(mgs) = 0.001*eii(mgs)*max((qx(mgs,li)-1.e-3),0.0)
             qscni(mgs) = min(qscni(mgs),qxmxd(mgs,li))
           endif
      else ! 10-ice version
      if ( iscni > 0 .and. qx(mgs,li) .gt. qxmin(li) ) then
          qscni(mgs) =    &
     &    pi*rho0(mgs)*((0.25)/(6.0))   &
     &    *eii(mgs)*(qx(mgs,li)**2)*(xdia(mgs,li,2))   &
     &    *vtxbar(mgs,li,1)/xmas(mgs,li)
         cscni(mgs) = Min(cimxd(mgs),qscni(mgs)*rho0(mgs)/xmas(mgs,li))
        end if

      end if
      end do





!
!
!  compute dry growth rate of snow, graupel, and hail
!
      do mgs = 1,ngscnt
!
      qsdry(mgs)  = qsacr(mgs)    + qsacw(mgs)   &
     &            + qsaci(mgs)
!
      qhdry(mgs)  = qhaci(mgs)    + qhacs(mgs)   &
     &            + qhacr(mgs)   &
     &            + qhacw(mgs)
!

      qhldry(mgs) = 0.0
      if ( lhl .gt. 1 ) then
      qhldry(mgs)  = qhlaci(mgs)    + qhlacs(mgs)   &
     &               + qhlacr(mgs)   &
     &               + qhlacw(mgs)
      endif
      end do
!
!  set wet growth and shedding
!
      do mgs = 1,ngscnt
      
      if ( tfrdry < temg(mgs) .and. temg(mgs) < tfr ) then
!
!      qswet(mgs) =
!     >  ( xdia(mgs,ls,1)*swvent(mgs)*cx(mgs,ls)*fwet1(mgs)
!     >  + fwet2(mgs)*(qsaci(mgs)+qsacir(mgs)
!     >               +qsacip(mgs)) )
!      qswet(mgs) = max( 0.0, qswet(mgs))
!
!      IF ( dnu(lh) .ne. 0. ) THEN
!        qhwet(mgs) = qhdry(mgs)
!      ELSE
       if ( incwet == 0 ) then
        qhwet(mgs) =   &
     &    ( xdia(mgs,lh,1)*hwvent(mgs)*cx(mgs,lh)*fwet1(mgs)   &
     &   + fwet2(mgs)*(qhaci(mgs) + qhacs(mgs)) )
       qhwet(mgs) = max( 0.0, qhwet(mgs))
         else
         endif

!      ENDIF


       qhlwet(mgs) = 0.0
       if ( lhl .gt. 1 ) then
         if ( incwet == 0 ) then
         qhlwet(mgs) =   &
     &     ( xdia(mgs,lhl,1)*hlvent(mgs)*cx(mgs,lhl)*fwet1(mgs)   &
     &     + fwet2(mgs)*(qhlaci(mgs) + qhlacs(mgs)) )
         qhlwet(mgs) = max( 0.0, qhlwet(mgs))
         
         else
         endif ! incwet
       endif
       
       else
       
        qhwet(mgs) = qhdry(mgs)
        qhlwet(mgs) = qhldry(mgs)
       endif
!
!      qhlwet(mgs) = qhldry(mgs)

      end do

!
! shedding rate
!
      qsshr(:)  =  0.0
      qhshr(:)  =  0.0
      qhlshr(:) =  0.0
      qhshh(:)  =  0.0
      csshr(:)  =  0.0
      csshrr(:) = 0.0
      chshr(:)  =  0.0
      chlshr(:)  =  0.0
      chshrr(:)  =  0.0
      chlshrr(:)  =  0.0
      vhshdr(:)  = 0.0
      vhlshdr(:)  = 0.0
      wetsfc(:)  = .false.
      wetgrowth(:)  = .false.
      wetsfchl(:)  = .false.
      wetgrowthhl(:)  = .false.

      do mgs = 1,ngscnt
!
!
!
      qhshr(mgs)  = Min( 0.0, qhwet(mgs) - qhdry(mgs) )  ! water that freezes should never be more than what sheds
      


      qhlshr(mgs)  =  Min( 0.0, qhlwet(mgs) - qhldry(mgs) )

!
! limit wet growth to only higher density particles
!
      qsshr(mgs)  =  0.0
!
!
!  no shedding for temperatures < 243.15 
!
      if ( temg(mgs) .lt. 243.15 ) then
       qsshr(mgs)  =  0.0
       qhshr(mgs)  =  0.0
       qhlshr(mgs) =  0.0
       vhshdr(mgs)  = 0.0
       vhlshdr(mgs)  = 0.0
       wetsfc(mgs) = .false.
       wetgrowth(mgs) = .false.
       wetsfchl(mgs) = .false.
       wetgrowthhl(mgs) = .false.
      end if
!
!  shed all at temperatures > 273.15
!
      if ( temg(mgs) .gt. tfr ) then

       if ( .false. ) then ! old and incorrect -- Thanks to Shaofeng Hua for noticing this error (9/17/2017)
       qsshr(mgs)  = -qsdry(mgs)
       qhshr(mgs)  = -qhdry(mgs)
       qhlshr(mgs) = -qhldry(mgs)
       else ! new and correct
       ! note that the qxacr terms should be zero here, so shedding at T > 0 is all from the droplets
       qsshr(mgs)   = - qsacr(mgs) - qsacw(mgs) ! -qsdry(mgs)
       qhlshr(mgs)  = - qhlacw(mgs) - qhlacr(mgs) ! -qhldry(mgs)
       qhshr(mgs)  = - qhacw(mgs) - qhacr(mgs) ! -qhdry(mgs)

       endif

       vhshdr(mgs)  = -vhacw(mgs) - vhacr(mgs)
       vhlshdr(mgs)  = -vhlacw(mgs) - vhlacr(mgs)
       qhwet(mgs)  = 0.0
       qhlwet(mgs) = 0.0
      end if
!
!      if (qhshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr  ) THEN
        wetsfc(mgs) =  (qhshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) .or. ( qhmlr(mgs) < -qxmin(lh) .and.  temg(mgs) > tfr )
        wetgrowth(mgs) = (qhshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr )
!      ENDIF
      if (qhlshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) then
        wetsfchl(mgs) = (qhlshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr ) .or. ( qhlmlr(mgs) < -qxmin(lhl) .and.  temg(mgs) > tfr )
        wetgrowthhl(mgs) = (qhlshr(mgs) .lt. 0.0 .and. temg(mgs) < tfr )
      endif

      end do
!
      if ( ipconc .ge. 1 ) then
      do mgs = 1,ngscnt
      csshr(mgs)  = 0.0 ! (cx(mgs,ls)/(qx(mgs,ls)+1.e-20))*Min(0.0,qsshr(mgs))
       
       chshr(mgs) = 0.0 ! no change to graupel number concentration for wet-growth shedding
       
      !   tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
        ! Base the drop size on the shedding regime
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chshrr(mgs) =  rho0(mgs)*qhshr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lh))  ! into rain 
      
      
      
      chlshr(mgs) = 0.0
      chlshrr(mgs) = 0.0
      if ( lhl .gt. 1 ) then 
!         chlshr(mgs)  = (cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlshr(mgs)


       chlshr(mgs) = 0.0 ! no change to hail number concentration for wet-growth shedding
       
      !   tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
        ! Base the drop size on the shedding regime
            ! 8/26/2015 ERM updated to use shedalp and tmpdiam
            ! tmpdiam = (shedalp+alpha(mgs,lh))*xdia(mgs,lh,1)
            chlshrr(mgs) =  rho0(mgs)*qhlshr(mgs)/(xdn(mgs,lr)*vshdgs(mgs,lhl))  ! into rain 

      endif ! ( lhl > 1 )

      
      end do
      end if



!
!  final decisions
!
      do mgs = 1,ngscnt
!
!  Snow
!
      if ( qsshr(mgs) .lt. 0.0 ) then
      qsdpv(mgs) = 0.0
      qssbv(mgs) = 0.0
      else
      qsshr(mgs) = 0.0
      end if
!
!     if ( qsdry(mgs) .lt. qswet(mgs) ) then
!     qswet(mgs) = 0.0
!     else
!     qsdry(mgs) = 0.0
!     end if
!

!  graupel
!
!
      if ( wetgrowth(mgs) .or. (mixedphase .and. fhw(mgs) .gt. 0.05 .and. temg(mgs) .gt. 243.15) ) then
      

! soaking (when not advected liquid water film with graupel)

        if ( lvol(lh) .gt. 1 .and. .not. mixedphase) then
        ! rescale volumes to maximum density
         if ( iwetsoak ) then 

         rimdn(mgs,lh) = xdnmx(lh)
         raindn(mgs,lh) = xdnmx(lh)
         vhacw(mgs) = qhacw(mgs)*rho0(mgs)/rimdn(mgs,lh)
         vhacr(mgs) = qhacr(mgs)*rho0(mgs)/raindn(mgs,lh)
!        IF ( lvol(lh) .gt. 1 .and. wetgrowth(mgs) ) THEN
         if ( xdn(mgs,lh) .lt. xdnmx(lh) ) then
         ! soak some liquid into the graupel
!           v1 = xdnmx(lh)*vx(mgs,lh)/(xdn(mgs,lh)*dtp) ! volume available for filling
           v1 = (1. - xdn(mgs,lh)/xdnmx(lh))*vx(mgs,lh)/(dtp) ! volume available for filling
!            tmp = (vx(mgs,lh)/rho0(mgs))*(xdnmx(lh) - xdn(mgs,lh)) ! max mixing ratio of liquid water that can be added
           v2 = rho0(mgs)*qhwet(mgs)/xdnmx(lh)  ! volume of frozen accretion
           
           vhsoak(mgs) = Min(v1,v2)

           
         endif
         
         endif

         vhshdr(mgs) = Min(0.0, rho0(mgs)*qhwet(mgs)/xdnmx(lh) - vhacw(mgs) - vhacr(mgs) )
         
        elseif ( lvol(lh) .gt. 1  .and. mixedphase ) then
!         vhacw(mgs) = rho0(mgs)*qhacw(mgs)/xdn0(lr)
!         vhacr(mgs) = rho0(mgs)*qhacr(mgs)/xdn0(lr)
        endif
        

      qhdpv(mgs) = 0.0
!      qhsbv(mgs) = 0.0
      chdpv(mgs) = 0.0
!      chsbv(mgs) = 0.0

! collection efficiency modification

      if ( ehi(mgs) .gt. 0.0 ) then
        qhaci(mgs) = Min(qimxd(mgs),qhaci0(mgs))  ! effectively sets collection eff to 1
        chaci(mgs) = Min(cimxd(mgs),chaci0(mgs))  ! effectively sets collection eff to 1
      endif
      if ( ehs(mgs) .gt. 0.0 ) then
!        qhacs(mgs) = Min(qsmxd(mgs),qhacs(mgs)/ehs(mgs))  ! effectively sets collection eff to 1
        qhacs(mgs) = Min(qsmxd(mgs),qhacs0(mgs)) !/ehs(mgs)                   ! divide out the collection efficiency
        chacs(mgs) = Min(csmxd(mgs),chacs0(mgs)) !/ehs(mgs)                   ! divide out the collection efficiency
        ehs(mgs) = ehsmax ! 1.0 ! min(ehsfrac*ehs(mgs),ehsmax)            ! modify it
        qhacs(mgs) = Min(qsmxd(mgs),qhacs(mgs))   ! plug it back in
      endif

! be sure to catch particles with wet surfaces but not in wet growth to turn off Hallett-Mossop
      wetsfc(mgs) = .true.

      else
!        qhshr(mgs) = 0.0
      end if
!
!
!  hail
!
!      if ( lhl .gt. 1 .and. qhlshr(mgs) .lt. 0.0 ) then
      if ( lhl > 1 .and. ( wetgrowthhl(mgs) .or. (mixedphase .and. fhlw(mgs) .gt. 0.05 .and. temg(mgs) .gt. 243.15) ) ) then
!      if ( wetgrowthhl(mgs) ) then
       

      qhldpv(mgs) = 0.0
!      qhlsbv(mgs) = 0.0
      chldpv(mgs) = 0.0
!      chlsbv(mgs) = 0.0




        if ( lvol(lhl) .gt. 1  .and. .not. mixedphase ) then
!        IF ( lvol(lhl) .gt. 1 .and. wetgrowthhl(mgs) ) THEN

         if ( iwetsoak ) then 

         rimdn(mgs,lhl) = xdnmx(lhl) 
         raindn(mgs,lhl) = xdnmx(lhl) 
         vhlacw(mgs) = qhlacw(mgs)*rho0(mgs)/rimdn(mgs,lhl)
         vhlacr(mgs) = qhlacr(mgs)*rho0(mgs)/raindn(mgs,lhl)

         if ( xdn(mgs,lhl) .lt. xdnmx(lhl) ) then
         ! soak some liquid into the hail
!           v1 = xdnmx(lhl)*vx(mgs,lhl)/(xdn(mgs,lhl)*dtp) ! volume available for filling
           v1 = (1. - xdn(mgs,lhl)/xdnmx(lhl))*vx(mgs,lhl)/(dtp) ! volume available for filling
!            tmp = (vx(mgs,lhl)/rho0(mgs))*(xdnmx(lhl) - xdn(mgs,lhl)) ! max mixing ratio of liquid water that can be added
           v2 = rho0(mgs)*qhlwet(mgs)/xdnmx(lhl)  ! volume of frozen accretion
           if ( v1 > v2 ) then ! all the frozen stuff fits in
             vhlsoak(mgs) = v2
           else  ! fill up the available space
             vhlsoak(mgs) = v1
           endif
!           vhlacw(mgs) = 0.0
!           vhlacr(mgs) = Max( 0.0, v2 - v1 )
         else
           vhlsoak(mgs) = 0.0
!           vhlacw(mgs) = 0.0
!           vhlacr(mgs) = rho0(mgs)*qhlwet(mgs)/raindn(mgs,lhl)
         
         endif
         
         endif

         vhlshdr(mgs) = Min(0.0, rho0(mgs)*qhlwet(mgs)/xdnmx(lhl) - vhlacw(mgs) - vhlacr(mgs) )


        elseif ( lvol(lhl) .gt. 1  .and. mixedphase ) then
!         vhlacw(mgs) = rho0(mgs)*qhlacw(mgs)/xdn0(lr)
!         vhlacr(mgs) = rho0(mgs)*qhlacr(mgs)/xdn0(lr)
        endif

      if ( ehli(mgs) .gt. 0.0 ) then
        qhlaci(mgs) = Min(qimxd(mgs),qhlaci0(mgs))  ! effectively sets collection eff to 1
        chlaci(mgs) = Min(cimxd(mgs),chlaci0(mgs))  ! effectively sets collection eff to 1
      endif

!      IF ( ehls(mgs) .gt. 0.0 ) THEN
!        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs(mgs)/ehls(mgs))
!      ENDIF
      if ( ehls(mgs) .gt. 0.0 ) then
        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs0(mgs)) !/ehls(mgs)                   ! divide out the collection efficiency
        chlacs(mgs) = Min(csmxd(mgs),chlacs0(mgs)) !/ehls(mgs)                   ! divide out the collection efficiency
        ehls(mgs) = ehsmax ! 1.0 ! min(ehsfrac*ehs(mgs),ehsmax)            ! modify it
!        qhlacs(mgs) = Min(qsmxd(mgs),qhlacs(mgs))   ! plug it back in
      endif

      
!      qhlwet(mgs) = 1.0

! be sure to catch particles with wet surfaces but not in wet growth to turn off Hallett-Mossop
      wetsfchl(mgs) = .true.


      else
!      qhlshr(mgs) = 0.0
!      qhlwet(mgs) = 0.0
      end if

      end do
!
! Ice -> graupel conversion
!
      do mgs = 1,ngscnt
      
      qhcni(mgs) = 0.0
      chcni(mgs) = 0.0
      chcnih(mgs) = 0.0
      vhcni(mgs) = 0.0
      
      if ( iglcnvi .ge. 1 ) then
      if ( temg(mgs) .lt. 273.0 .and. qiacw(mgs) - qidpv(mgs) .gt. 0.0 ) then
      
        
        tmp = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,li,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
        tmp = Min( Max( rimc3, tmp ), 900.0 )
        
        !  Assume that half the volume of the embryo is rime with density 'tmp'
        !  m = rhoi*(V/2) + rhorime*(V/2) = (rhoi + rhorime)*V/2
        !  V = 2*m/(rhoi + rhorime)
        
!        write(0,*)  'rime dens = ',tmp
        
        if ( tmp .ge. 200.0 .or. iglcnvi >= 2 ) then
          r = Max( 0.5*(xdn(mgs,li) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcni(mgs) = (qiacw(mgs) - qidpv(mgs)) ! *float(iglcnvi)
          chcni(mgs) = cx(mgs,li)*qhcni(mgs)/qx(mgs,li)
!          chcnih(mgs) = rho0(mgs)*qhcni(mgs)/(1.6e-10)
          chcnih(mgs) = Min(chcni(mgs), rho0(mgs)*qhcni(mgs)/(r*xvmn(lh)) )
!          vhcni(mgs) = rho0(mgs)*2.0*qhcni(mgs)/(xdn(mgs,li) + tmp)
          vhcni(mgs) = rho0(mgs)*qhcni(mgs)/r
        endif
      
      elseif ( iglcnvi == 3 ) then

       if ( temg(mgs) .lt. 273.0 .and. qiacw(mgs)*dtp > 2.*qxmin(lh) .and. gamice73fac*xmas(mgs,li) > xdnmn(lh)*xvmn(lh) ) then
      
        
        tmp = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,li,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
        tmp = Min( Max( rimc3, tmp ), 900.0 )
        
        !  Assume that half the volume of the embryo is rime with density 'tmp'
        !  m = rhoi*(V/2) + rhorime*(V/2) = (rhoi + rhorime)*V/2
        !  V = 2*m/(rhoi + rhorime)
        
!        write(0,*)  'rime dens = ',tmp
        ! convert to particles with the mass of the mass-weighted diameter
      !  massofmwr = gamice73fac*xmas(mgs,li)
        
        if ( tmp .ge. xdnmn(lh)  ) then
          r = Max( 0.5*(xdn(mgs,li) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcni(mgs) = 0.5*qiacw(mgs)
          chcni(mgs) = qhcni(mgs)/(gamice73fac*xmas(mgs,li))
          chcnih(mgs) = Min(chcni(mgs), rho0(mgs)*qhcni(mgs)/(r*xvmn(lh)) )
!          vhcni(mgs) = rho0(mgs)*2.0*qhcni(mgs)/(xdn(mgs,li) + tmp)
          vhcni(mgs) = rho0(mgs)*qhcni(mgs)/r
        endif
      
      endif

      
      endif
      endif
      
      
      enddo
      
      
      qhlcnh(:) = 0.0
      chlcnh(:) = 0.0
      chlcnhhl(:) = 0.0
      vhlcnh(:) = 0.0
      vhlcnhl(:) = 0.0
      zhlcnh(:) = 0.0

      qhcnhl(:) = 0.0
      chcnhl(:) = 0.0
      vhcnhl(:) = 0.0
      zhcnhl(:) = 0.0


      if ( lhl .gt. 1  ) then
      
      if ( ihlcnh == 1 .or. ihlcnh == 3 ) then

!
!  Graupel (h) conversion to hail (hl) based on Milbrandt and Yau 2005b
!
      do mgs = 1,ngscnt

!        IF ( lhl .gt. 1 .and. ipconc .ge. 5 .and. qx(mgs,lh) .gt. 1.0e-3 .and.
!     :        xdn(mgs,lh) .gt. 750. .and. qhshr(mgs) .lt. 0.0 .and.
!     :        xdia(mgs,lh,3) .gt. 1.e-3 ) THEN
        if ( hlcnhdia > 0 ) then
          ltest = xdia(mgs,lh,3) .gt. hlcnhdia  ! test on mean volume diameter
        else 
!          ltest =  xdia(mgs,lh,1)*(3. + alpha(mgs,lh)) > Abs( hlcnhdia ) ! test on maximum mass diameter
          ltest =  xdia(mgs,lh,1)*(4. + alpha(mgs,lh)) > Abs( hlcnhdia ) ! test on mass-weighted diameter
        endif

         if ( iusedw == 0 .and. ihlcnh == 1 ) then
           dg0(mgs) = -1.
         else
         if (((qhacw(mgs) + qhacr(mgs))*dtp > qxmin(lh) .and. qx(mgs,lh) > hlcnhqmin .and. temg(mgs) .le. tfr-2.0  &
               .and.  temg(mgs) .gt. dwtempmin ) .or. ( wetgrowth(mgs) .and. qx(mgs,lh) > hlcnhqmin ) ) then
!         dw = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*ehw(mgs)*qx(mgs,lc) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
!         dwr = 0.01*( Exp( -temcg(mgs)/( 1.1e4 * rho0(mgs)*(ehw(mgs)*qx(mgs,lc)+ehr(mgs)*qx(mgs,lr)) - &
!                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 ) ) - 1.0 )
            x =   1.1e4 * rho0(mgs)*(ehw(mgs)*qx(mgs,lc)+ehr(mgs)*qx(mgs,lr)) - &
                1.3e3*rho0(mgs)*qx(mgs,li) + 1.0 
            if ( x > 1.e-20 ) then
              arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
              dwr = 0.01*(exp(arg) - 1.0)
            else
              dwr = 1.e30
            endif
          d = dwr
           if ( dwr < 0.2 .and. dwr > 0.0 .and. rho0(mgs)*(qx(mgs,lc)+qx(mgs,lr)) > 1.e-4 ) then
                      sqrtrhovt = Sqrt( rhovt(mgs) )
                      fventh = sqrtrhovt*(fpndl(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5) 
                      fventm = sqrtrhovt*(fschm(mgs)**(1./3.)) * (fakvisc(mgs))**(-0.5)
                      ltemq = (tfr-163.15)/fqsat+1.5
                      qvs0 = pqs(mgs)*tabqvs(ltemq)
                      denomdp = felf(mgs) + fcw(mgs)*temcg(mgs)
                      denominvdp = 1.d0/(felf(mgs) + fcw(mgs)*temcg(mgs))

!                      write(91,*) 'dw,dwr,temcg = ',100.*dw,100.*dwr,temcg(mgs)
                      h1 = ( -ftka(mgs)*temcg(mgs) - felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qx(mgs,lv) - qvs0) )
                      h2 = ehi(mgs)*qx(mgs,li)*rho0(mgs)*fci(mgs)*temcg(mgs)
                      h3 = Max(dwehwmin, ehw(mgs))*qx(mgs,lc) 
                      h4 = ehr(mgs)* qx(mgs,lr)
                      ! iterate to find minimum diameter for wet growth. Start with value of dwr
                      do n = 1,10
                        d = Max(d, 1.e-4)
                        dold = d
                        vth = axx(mgs,lh)*d**bxx(mgs,lh) 
                        x2 = fventh*sqrtrhovt*Sqrt(d*vth)
                       if ( x2 > 1.4 ) then
                         ah = 0.78 + 0.308*x2  ! heat ventillation
                       else
                         ah = 1.0 + 0.108*x2**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                       endif

                       if ( .false. ) then ! this option includes 'am' separate from ah, which makes only small differences. Otherwise equivalent to second option
                        x1 = fventm*sqrtrhovt*Sqrt(d*vth)
                        if ( x1 > 1.4 ) then
                          am = 0.78 + 0.308*x1 ! mass ventillation (Beard and Pruppacher 1971, eq. 8)
                        else
                          am = 1.0 + 0.108*x1**2 ! mass ventillation (Beard and Pruppacher 1971, eq. 9)
                        endif
                        
                        d = 8.*denominvdp*( am*felv(mgs)*fwvdf(mgs)*rho0(mgs)*(qvs0 - qx(mgs,lv)) - ah*ftka(mgs)*temcg(mgs)  )/ &
                           (dtp* ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs) +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2*denominvdp))
                       
                        else

                        ! Based on Farley and Orville (1986), eq. 5-9 but neglecting the Ci*(T0-Ts) term in (8) since we want Ts=T0
                        ! Simplified mass rates as dm_w/dt = pi/4*d**2*(Vh - Vc)*rhoair*qc*ehw, etc.
                        d = 8.*ah*h1/ &
                            ( ( Max(0.001,vth - vtxbar(mgs,lc,1))*h3 +                              &
                            Max(0.001,vth - vtxbar(mgs,lr,1))*h4) *rho0(mgs)*denomdp +               &
                            Max(0.001,vth - vtxbar(mgs,li,1))*h2)
                            
                        endif
                        if ( Abs(dold - d)/dold < 0.05 .or. ( n > 3 .and. d > dg0thresh ) ) exit
                        
                      enddo
              endif
              
              dg0(mgs) = Max( d, dwmin )
          else
            dg0(mgs) = dg0thresh + 0.0001
          endif
          
            if ( ihlcnh == 3 .and. (qhacw(mgs) + qhacr(mgs))*dtp > qxmin(lh) .and. qx(mgs,lh) > hlcnhqmin &
                   .and. temg(mgs) .le. tfr-2.0 ) then
           ! set a secondary condition on to capture large graupel that is riming but not in wet growth
                dg0(mgs) = Min( dg0(mgs), dg0thresh - 0.0001 )
            endif
            
          endif

        wtest = (dg0(mgs) > 0.0 .and. dg0(mgs) < dg0thresh )
        
        if ( ihlcnh == 1 ) then ! .or. iusedw == 0  then
        
        if ( ( wetgrowth(mgs) .and. (xdn(mgs,lh) .gt. hldnmn .or. lvh < 1 ) .and.  & ! correct this when hail gets turned on
     &        rimdn(mgs,lh) .gt. 800. .and.   &
     &        ltest .and. qx(mgs,lh) .gt. hlcnhqmin ) .or. wtest ) then ! {
!     :        xdia(mgs,lh,3) .gt. 2.e-3 .and. qx(mgs,lh) .gt. 1.0e-3  THEN ! 0823.2008 erm test
!        IF ( xdia(mgs,lh,3) .gt. 1.e-3 ) THEN
        if ( qhacw(mgs) .gt. 0.0 .and. qhacw(mgs) .gt. qhaci(mgs) .and. temg(mgs) .le. tfr-2.0 ) then ! {
        ! dh0 is the diameter dividing wet growth from dry growth (Ziegler 1985), modified by MY05
!          dh0 = 0.01*(exp(temcg(mgs)/(1.1e4*(qx(mgs,lc)+qx(mgs,lr)) - 
!     :           1.3e3*qx(mgs,li) + 1.0e-3 ) ) - 1.0)
          if ( wtest ) then
            dh0 = dg0(mgs)
          else
            x = (1.1e4*(rho0(mgs)*qx(mgs,lc)) - 1.3e3*rho0(mgs)*qx(mgs,li) + 1.0e-3 )
            if ( x > 1.e-20 ) then
            arg = Min(70.0, (-temcg(mgs)/x )) ! prevent overflow of the exp function in 32 bit
            dh0 = 0.01*(exp(arg) - 1.0)
            else
             dh0 = 1.e30
            endif
          endif ! wtest
!          dh0 = Max( dh0, 5.e-3 )
          
!         IF ( dh0 .gt. 0.0 ) write(0,*) 'dh0 = ',dh0
!         IF ( dh0 .gt. 1.0e-4 ) THEN
         if ( xdia(mgs,lh,3)/dh0 .gt. 0.1 ) then !{
!         IF ( xdia(mgs,lh,3) .lt. 20.*dh0 .and. dh0 .lt. 2.0*xdia(mgs,lh,3) ) THEN 
           tmp = qhacw(mgs) + qhacr(mgs) + qhaci(mgs) + qhacs(mgs)
!           qtmp = Min( 1.0, xdia(mgs,lh,3)/(2.0*dh0) )*(tmp)
           qtmp = Min( 100.0, xdia(mgs,lh,3)/(2.0*dh0) )*(tmp)
           qhlcnh(mgs) = Min(  qxmxd(mgs,lh), qtmp )
           
           if ( ipconc .ge. 5 ) then !{
!           dh0 = Max( xdia(mgs,lh,3), Min( dh0, 5.e-3 ) ) ! do not create hail greater than 5mm diam. unless the graupel is larger
           if ( .not. wtest ) dh0 = Min( dh0, 10.e-3 ) ! do not create hail greater than 10mm diam., which is the max graupel size
           if ( qx(mgs,lhl) > 0.1e-3 ) dh0 = Max( dh0, xdia(mgs,lhl,3) ) ! when enough hail is established, do not dilute the size
           chlcnhhl(mgs) = Min( cxmxd(mgs,lh), rho0(mgs)*qhlcnh(mgs)/(pi*xdn(mgs,lh)*dh0**3/6.0) )

           r = rho0(mgs)*qhlcnh(mgs)/(xdn(mgs,lh)*xv(mgs,lh))  ! number of graupel particles at mean volume diameter
           chlcnh(mgs) = Max( chlcnhhl(mgs), r )
           endif !}
           
           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))

          endif !}

        endif ! }
        endif ! }
        
        elseif ( ihlcnh == 3 ) then !{
         

          if ( wtest  .and. &
               ( qhacw(mgs)*dtp > qxmin(lh) .and. temg(mgs) .lt. tfr-2. .and. qx(mgs,lh) > hlcnhqmin ) ) then
        ! convert number, mass, and reflectivity for d > dw
           if ( ipconc == 5 ) then
            ! dg0(mgs) = Min( dg0(mgs), hldia1 )
             !dg0(mgs) = hldia1
           endif
           
           ratio = Min( maxratiolu, dg0(mgs)/xdia(mgs,lh,1) )


           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lh),4,1)
           if ( ipconc == 5 ) then
       !      tmp2 = Min( 0.25, tmp2 )
           endif
            qxd1 = qx(mgs,lh)*(tmp2)
            qhlcnh(mgs) = dtpinv*qxd1
            flim = 1.0
!            tmp3 = Min( dtp*(qfacw(mgs) + qfacr(mgs) ), qxmxd(mgs,lf) )
            tmp3 = qxmxd(mgs,lh)
            if (qxd1 > tmp3 ) then
              flim = tmp3/(qxd1)
              qhlcnh(mgs) = flim*qhlcnh(mgs)
            endif

            
            
            if ( ( qxd1 > qxmin(lhl) .and. ipconc > 5 ) .or. ( qxd1 > 10.*qxmin(lhl) .and. ipconc == 5) ) then
            
           ! number
            tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
             if ( ipconc == 5 ) then
          !     tmp = Min( 0.2, tmp )
             endif
            cxd1 = flim*cx(mgs,lh)*( tmp)
            chlcnh(mgs) = dtpinv*cxd1
            chlcnhhl(mgs) = chlcnh(mgs)

           if ( qx(mgs,lhl) > qxmin(lhl) .and. dmhlopt > 0 ) then
             dh0 = rho0(mgs)*qhlcnh(mgs)/chlcnhhl(mgs)
             if ( dh0 < xmas(mgs,lhl) ) then
               ! dh0 = ( qxd1*dh0 + qx(mgs,lhl)*xmas(mgs,lhl))/( qxd1 + qx(mgs,lhl))  ! weighted average
               dh0 = (( qxd1*dh0**(1./3.) + qx(mgs,lhl)*xmas(mgs,lhl)**(1./3.))/( qxd1 + qx(mgs,lhl)))**3  ! weighted average
               chlcnhhl(mgs) = Min( chlcnhhl(mgs), rho0(mgs)*qhlcnh(mgs)/dh0 )
             else
!               dh0 = Max( dh0, xmas(mgs,lhl) ) ! when enough hail is established, do not dilute the size
             endif
           endif


           ! reflectivity
           if ( ipconc >= 6 .and. lzh > 1 .and. lzhl > 1 ) then
            tmp3 = gaminterp(ratio,alpha(mgs,lh),11,1)
            zxd1 = flim*zx(mgs,lh)*(tmp3)
            zhlcnh(mgs) = dtpinv*zxd1
           else
            zxd1 = 0
           endif

            else
               qhlcnh(mgs) = 0.0
            endif

           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))
           
           endif


        endif !}
      
      enddo
      
      elseif ( ihlcnh == 2 ) then ! 10-ice type conversion 

!
! Staka and Mansell (2005) type conversion
!
!      hldia1 is set in micro_module and namelist
!      IF ( .true. ) THEN
      
        ! convert number, mass, and reflectivity for d > hldia1,
        ! regardless of wet growth status, but as long as riming > 0
        do mgs = 1,ngscnt
        if ( qhacw(mgs)*dtp > qxmin(lh) .and. temg(mgs) .lt. tfr-2. .and. qx(mgs,lh) > qxmin(lh) ) then
           ratio = Min( maxratiolu, hldia1/xdia(mgs,lh,1) )

           ! number
            tmp = gaminterp(ratio,alpha(mgs,lh),1,1)
            cxd1 = cx(mgs,lh)*( tmp)
            chlcnh(mgs) = dtpinv*cxd1
            chlcnhhl(mgs) = chlcnh(mgs)

           ! mass
            tmp2 = gaminterp(ratio,alpha(mgs,lh),4,1)
            qxd1 = qx(mgs,lh)*(tmp2)
            qhlcnh(mgs) = dtpinv*qxd1

           ! reflectivity
           if ( lzh > 1 .and. lzhl > 1 ) then
            tmp3 = gaminterp(ratio,alpha(mgs,lh),11,1)
            zxd1 = zx(mgs,lh)*(tmp3)
            zhlcnh(mgs) = dtpinv*zxd1
           else
            zxd1 = 0
           endif
           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))
           
         endif
         
         enddo
!        ENDIF
      elseif ( ihlcnh == 0 ) then

      do mgs = 1,ngscnt
!      qhlcnh(mgs) = 0.0
!      chlcnh(mgs) = 0.0
      if ( wetgrowth(mgs) .and. temg(mgs) .lt. tfr-5. .and. qx(mgs,lh) > qxmin(lh) ) then
      if ( qhacw(mgs).gt.1.e-6 .and. xdn(mgs,lh) > 700. ) then
      qhlcnh(mgs) =                                                   &
        ((pi*xdn(mgs,lh)*cx(mgs,lh)) / (6.0*rho0(mgs)*dtp))           &
       *exp(-hldia1/xdia(mgs,lh,1))                                    &
       *( (hldia1**3) + 3.0*(hldia1**2)*xdia(mgs,lh,1)                  &
        + 6.0*(hldia1)*(xdia(mgs,lh,1)**2) + 6.0*(xdia(mgs,lh,1)**3) )
      qhlcnh(mgs) =   min(qhlcnh(mgs),qhmxd(mgs))
      if ( ipconc .ge. 5 ) then
        chlcnh(mgs) = Min( cxmxd(mgs,lh), cx(mgs,lh)*Exp(-hldia1/xdia(mgs,lh,1)))
        chlcnhhl(mgs) = chlcnh(mgs)
!        chlcnh(mgs) = Min( cxmxd(mgs,lh), rho0(mgs)*qhlcnh(mgs)/(2.0*xmas(mgs,lh) ))
      endif
           vhlcnh(mgs) = rho0(mgs)*qhlcnh(mgs)/xdn(mgs,lh)
           vhlcnhl(mgs) = rho0(mgs)*qhlcnh(mgs)/Max(xdnmn(lhl), xdn(mgs,lh))
      end if
      end if
      end do
      
!      ENDIF ! true
      
      endif ! ihlcnh options
      
     ! convert low-density hail to graupel
      if ( icvhl2h >= 1 ) then
      do mgs = 1,ngscnt
        if (  qx(mgs,lhl) > qxmin(lhl) .and. xdn(mgs,lhl) < 0.5*(xdnmn(lhl) + xdnmx(lhl)) ) then
          tmp = Min(0.95, 1. - 0.5*(1. + tanh(0.125*(xdn(mgs,lhl) - 1.01*xdnmn(lhl) )) ))
          qhcnhl(mgs) = tmp*qx(mgs,lhl)*dtpinv
          chcnhl(mgs) = cx(mgs,lhl)*qhcnhl(mgs)/qx(mgs,lhl)
          vhcnhl(mgs) = vx(mgs,lhl)*qhcnhl(mgs)/qx(mgs,lhl)
          
        endif
      enddo
      
      endif
      
      endif ! lhl > 1

  
  

!
! Ziegler snow conversion to graupel
!
      do mgs = 1,ngscnt

      qhcns(mgs) = 0.0
      chcns(mgs) = 0.0
      chcnsh(mgs) = 0.0
      vhcns(mgs) = 0.0

      qscnh(mgs) = 0.0
      cscnh(mgs) = 0.0
      vscnh(mgs) = 0.0

      if ( ipconc .ge. 5 ) then

        ! test attempt at converting graupel to snow when not riming but growing by deposition
        if ( temg(mgs) < tfr .and. qx(mgs,lh) .gt. qxmin(lh) .and. qhdpv(mgs) > qxmin(lh)*dtpinv  &
     &       .and. qhacw(mgs) < qxmin(lh)*dtpinv ) then
          if ( xdn(mgs,lh) < 290. ) then
!          qscnh(mgs) = 2.*qhdpv(mgs)
!          cscnh(mgs) = cx(mgs,lh)*qscnh(mgs)/qx(mgs,lh)
!          vscnh(mgs) = rho0(mgs)*qscnh(mgs)/xdn(mgs,lh)
          endif
        endif


        if ( qx(mgs,ls) .gt. qxmin(ls) .and. qsacw(mgs) .gt. 0.0 ) then

!      DATA VGRA/1.413E-2/  ! this is the volume (cm**3) of a 3mm diam. sphere
!    vgra = 1.4137e-8 m**3

!      DNNET=DNCNV-DNAGG
!      DQNET=QXCON+QSACC+SDEP
!
!      DNSCNV=EXP(-(ROS*XNS*VGRA/(RO*QI)))*((1.-(XNS*VGRA*ROS/
!     / (RO*QI)))*DNNET + (XNS**2*VGRA*ROS/(RO*QI**2))*DQNET)
!      IF(DNSCNV.LT.0.) DNSCNV=0.
!
!      QIHC=(ROS*VGRA/RO)*DNSCNV
!
!      QH=QH+DT*QIHC
!      QI=QI-DT*QIHC
!      XNH=XNH+DT*DNSCNV
!      XNS=XNS-DT*DNSCNV

        if ( iglcnvs .eq. 1 ) then  ! Zrnic, Ziegler et al (1993)

        dnnet = cscnvis(mgs) + cscnis(mgs) - csacs(mgs)
        dqnet = qscnvi(mgs) + qscni(mgs) + qsacw(mgs) + qsdpv(mgs) + qssbv(mgs)

        a3 = 1./(rho0(mgs)*qx(mgs,ls))
        a1 = Exp( - xdn(mgs,ls)*cx(mgs,ls)*vgra*a3 ) !! EXP(-(ROS*XNS*VGRA/(RO*QI)))
! (1.-(XNS*VGRA*ROS/(RO*QI)))*DNNET
        a2 =  (1.-(cx(mgs,ls)*vgra*xdn(mgs,ls)*a3))*dnnet
! (XNS**2*VGRA*ROS/(RO*QI**2))*DQNET
        a4 = cx(mgs,ls)**2*vgra*xdn(mgs,ls)*a3/qx(mgs,ls)*dqnet

        chcns(mgs) = Max( 0.0, a1*(a2 + a4) )
        chcns(mgs) = Min( chcns(mgs), cxmxd(mgs,ls) )
        chcnsh(mgs) = chcns(mgs)

        qhcns(mgs) = Min( xdn(mgs,ls)*vgra*rhoinv(mgs)*chcns(mgs), qxmxd(mgs,ls) )
        vhcns(mgs) = rho0(mgs)*qhcns(mgs)/Max(xdn(mgs,ls),xdnmn(lh))
!        vhcns(mgs) = rho0(mgs)*qhcns(mgs)/Max(xdn(mgs,ls),400.)

        elseif ( iglcnvs .ge. 2  ) then  ! treat like ice crystals, i.e., check for rime density (ERM)

          if ( temg(mgs) .lt. 273.0 .and. ( qsacw(mgs) - qsdpv(mgs) .gt. 0.0 .or. &
              ( iglcnvs >= 3 .and. qsacw(mgs)*dtp > 2.*qxmin(lh) .and. gamsnow73fac*xmas(mgs,ls) > xdnmn(lh)*xvmn(lh)  ) ) ) then !{


        tmp = rimc1*(-((0.5)*(1.e+06)*xdia(mgs,lc,1))   &
     &                *((0.60)*vtxbar(mgs,ls,1))   &
     &                /(temg(mgs)-273.15))**(rimc2)
!        tmp = Min( Max( rimc3, tmp ), 900.0 )
        tmp = Min( tmp , 900.0 )

        !  Assume that half the volume of the embryo is rime with density 'tmp'
        !  m = rhoi*(V/2) + rhorime*(V/2) = (rhoi + rhorime)*V/2
        !  V = 2*m/(rhoi + rhorime)

!        write(0,*)  'rime dens = ',tmp

        if ( iglcnvs == 2 ) then !{
        if ( tmp .ge. 200.0  ) then
          r = Max( 0.5*(xdn(mgs,ls) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcns(mgs) = (qsacw(mgs) - qsdpv(mgs))
          chcns(mgs) = cx(mgs,ls)*qhcns(mgs)/qx(mgs,ls)
!          chcnih(mgs) = rho0(mgs)*qhcni(mgs)/(1.6e-10)
          chcnsh(mgs) = Min(chcns(mgs), rho0(mgs)*qhcns(mgs)/(r*xvmn(lh)) )
!          vhcni(mgs) = rho0(mgs)*2.0*qhcni(mgs)/(xdn(mgs,li) + tmp)
          vhcns(mgs) = rho0(mgs)*qhcns(mgs)/r
        endif
        
        elseif ( iglcnvs == 3 ) then
 
         ! convert to particles with the mass of the mass-weighted diameter
      !  massofmwr = gamice73fac*xmas(mgs,li)
        
        if ( tmp > xdnmn(lh) ) then
          r = Max( 0.5*(xdn(mgs,ls) + tmp), xdnmn(lh) )
!          r = Max( r, 400. )
          qhcns(mgs) = 0.5*qsacw(mgs)
          chcns(mgs) = qhcns(mgs)/(gamsnow73fac*xmas(mgs,ls))
          chcns(mgs) = Min( chcns(mgs), cx(mgs,ls)*qhcns(mgs)/qx(mgs,ls))
          chcnsh(mgs) = Min(chcns(mgs), rho0(mgs)*qhcns(mgs)/(r*xvmn(lh)) )
          vhcns(mgs) = rho0(mgs)*qhcns(mgs)/r
        endif

        endif !}

      endif !}

        endif


        endif

       else ! single moment lfo

        qhcns(mgs) = 0.001*ehscnv(mgs)*max((qx(mgs,ls)-6.e-4),0.0)
        qhcns(mgs) = min(qhcns(mgs),qxmxd(mgs,ls))
        if ( lvol(lh) .ge. 1 ) vhcns(mgs) = rho0(mgs)*qhcns(mgs)/Max(xdn(mgs,ls),400.)

       endif
      enddo
!
!
!  heat budget for rain---not all rain that collects ice can freeze
!
!
!
      if ( irwfrz .gt. 0 .and. .not. mixedphase) then
!
      do mgs = 1,ngscnt
!
!  compute total rain that freeze when it interacts with cloud ice
!
      qrztot(mgs) = qrfrz(mgs) + qiacr(mgs) + qsacr(mgs)
!
!  compute the maximum amount of rain that can freeze
!  Used to limit freezing to 4*qrmxd, but now allow all rain to freeze if possible
!
      qrzmax(mgs) =   &
     &  ( xdia(mgs,lr,1)*rwvent(mgs)*cx(mgs,lr)*fwet1(mgs) )
      qrzmax(mgs) = max(qrzmax(mgs), 0.0)
      qrzmax(mgs) = min(qrztot(mgs), qrzmax(mgs))
      qrzmax(mgs) = min(qx(mgs,lr)*dtpinv, qrzmax(mgs))

      if ( temcg(mgs) < -30. ) then ! allow all to freeze if T < -30 because fwet becomes invalid (negative)
        qrzmax(mgs) = qx(mgs,lr)*dtpinv
      endif
!      qrzmax(mgs) = min(4.*qrmxd(mgs), qrzmax(mgs))
!
!  compute the correction factor
!
!      IF ( qrztot(mgs) .gt. qxmin(lr) ) THEN
      if ( qrztot(mgs) .gt. qrzmax(mgs) .and. qrztot(mgs) .gt. qxmin(lr) ) then
        qrzfac(mgs) = qrzmax(mgs)/(qrztot(mgs))
      else
        qrzfac(mgs) = 1.0
      endif
      qrzfac(mgs) = min(1.0, qrzfac(mgs))
!
      end do
!
!
! now correct the above sources
!
!
      do mgs = 1,ngscnt
      if ( temg(mgs) .le. 273.15 .and. qrzfac(mgs) .lt. 1.0 ) then
      qrfrz(mgs)   = qrzfac(mgs)*qrfrz(mgs)
      qrfrzs(mgs)  = qrzfac(mgs)*qrfrzs(mgs)
      qrfrzf(mgs)  = qrzfac(mgs)*qrfrzf(mgs)
      qiacr(mgs)   = qrzfac(mgs)*qiacr(mgs)
      qsacr(mgs)   = qrzfac(mgs)*qsacr(mgs)
      qiacrf(mgs)  = qrzfac(mgs)*qiacrf(mgs)
      qiacrs(mgs)  = qrzfac(mgs)*qiacrs(mgs)
      crfrz(mgs)   = qrzfac(mgs)*crfrz(mgs)
      crfrzf(mgs)  = qrzfac(mgs)*crfrzf(mgs)
      crfrzs(mgs)  = qrzfac(mgs)*crfrzs(mgs)
      ciacr(mgs)   = qrzfac(mgs)*ciacr(mgs)
      ciacrf(mgs)  = qrzfac(mgs)*ciacrf(mgs)
      ciacrs(mgs)  = qrzfac(mgs)*ciacrs(mgs)

!      IF ( lzh .gt. 1 ) THEN
!        zrfrzf(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.)) * &
!        ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
!      ENDIF
      
       vrfrzf(mgs)  = qrzfac(mgs)*vrfrzf(mgs)
       viacrf(mgs)  = qrzfac(mgs)*viacrf(mgs)
      end if
      end do
!
!
!
      end if
!
!
!
!  evaporation of rain
!
!
!
      qrcev(:) = 0.0
      crcev(:) = 0.0


      do mgs = 1,ngscnt
!
      if ( qx(mgs,lr) .gt. qxmin(lr) ) then

      qrcev(mgs) =   &
     &  fvce(mgs)*cx(mgs,lr)*rwvent(mgs)*rwcap(mgs)*evapfac
! this line to allow condensation on rain:
      if ( rcond .eq. 1 ) then
        qrcev(mgs) = min(qrcev(mgs), qxmxd(mgs,lv))
! this line to have evaporation only:
      else
        qrcev(mgs) = min(qrcev(mgs), 0.0)
      endif

      qrcev(mgs) = max(qrcev(mgs), -qrmxd(mgs))
!      if ( temg(mgs) .lt. 273.15 ) qrcev(mgs) = 0.0
      if ( qrcev(mgs) .lt. 0. .and. lnr > 1 ) then
!        qrcev(mgs) =   -qrmxd(mgs)
!        crcev(mgs) = (rho0(mgs)/(xmas(mgs,lr)+1.e-20))*qrcev(mgs)
        if ( icrcev == 1 ) then
          crcev(mgs) = (cx(mgs,lr)/(qx(mgs,lr)))*qrcev(mgs)
        elseif ( icrcev == 2 ) then
          crcev(mgs) = (cx(mgs,lr)/(qx(mgs,lr)))*qrcev(mgs)*vtxbar(mgs,lr,2)/vtxbar(mgs,lr,1)
        else
          crcev(mgs) = 0.0
        endif
      else
         crcev(mgs) = 0.0
      endif
!      if ( temg(mgs) .lt. 273.15 ) crcev(mgs) = 0.0
!
      endif

      end do
!
! evaporation/condensation of wet graupel and snow
!
      if ( lhwlg > 1 ) then
      qhcevlg(:) = 0.0
      chcevlg(:) = 0.0
      endif
      if ( lhlwlg > 1 ) then
      qhlcevlg(:) = 0.0
      chlcevlg(:) = 0.0
      endif


!
!
!
!  ICE MULTIPLICATION: Two modes (rimpa, and rimpb)
!  (following Cotton et al. 1986)
!
 
      chmul1(:) =  0.0
      chlmul1(:) =  0.0
      csmul1(:) = 0.0
!
      qhmul1(:) =  0.0
      qhlmul1(:) =  0.0
      qsmul1(:) =  0.0
      do mgs = 1,ngscnt
 
       ltest =  qx(mgs,lh) .gt. qxmin(lh)
       if ( lhl > 1 )  ltest =  ltest .or. qx(mgs,lhl) .gt. qxmin(lhl)
       
      if ( (itype1 .ge. 1 .or. itype2 .ge. 1 )   &
     &              .and. qx(mgs,lc) .gt. qxmin(lc)) then
      if ( temg(mgs) .ge. 265.15 .and. temg(mgs) .le. 271.15 ) then
       if ( ipconc .ge. 2 ) then
        if ( xv(mgs,lc) .gt. 0.0     &
     &     .and.  ltest &
!     .and. itype2 .ge. 2    &
     &       ) then
!
!  Ziegler et al. 1986 Hallett-Mossop process.  VSTAR = 7.23e-15 (vol of 12micron radius)
!
         if ( alpha(mgs,lc) == 0.0 ) then
           ex1 = (1./250.)*Exp(-7.23e-15/xv(mgs,lc))
         else
           
           ratio = (1. + alpha(mgs,lc))*(7.23e-15)/xv(mgs,lc)

           if ( usegamxinfcnu ) then
            i = Nint(dgami*(1. + alpha(mgs,lc)))
            gcnup1 = gmoi(i)
            ex1 = (1./250.)*Gamxinf(1.+alpha(mgs,lc), ratio)/(gcnup1)
           else
             ratio = Min( maxratiolu, ratio )
             tmp = gaminterp(ratio,alpha(mgs,lc),1,1)
             ex1 = (1./250.)*tmp
           endif
         endif
       if ( itype2 .le. 2 ) then
         ft = Max(0.0,Min(1.0,-0.11*temcg(mgs)**2 - 1.1*temcg(mgs)-1.7))
       else
        if ( temg(mgs) .ge. 265.15 .and. temg(mgs) .le. 267.15 ) then
          ft = 0.5
        elseif (temg(mgs) .ge. 267.15 .and. temg(mgs) .le. 269.15 ) then
          ft = 1.0
        elseif (temg(mgs) .ge. 269.15 .and. temg(mgs) .le. 271.15 ) then
          ft = 0.5
        else 
          ft = 0.0
        endif
       endif
!        rhoinv = 1./rho0(mgs)
!        DNSTAR = ex1*cglacw(mgs)
        
       if ( ft > 0.0 ) then
        
        if ( itype2 > 0 ) then
         if ( qx(mgs,lh) .gt. qxmin(lh) .and. (.not. wetsfc(mgs))  ) then
          chmul1(mgs) = ft*ex1*chacw(mgs)
!          chmul1(mgs) = Min( ft*ex1*chacw(mgs), ft*(30.*1.e+06)*rho0(mgs)*qhacw(mgs) ) ! 1.e+6 converts kg to mg; Saunders & Hosseini (2001) average of about 30 crystals per mg
          qhmul1(mgs) = cimas0*chmul1(mgs)*rhoinv(mgs)
         endif
         if ( lhl .gt. 1 ) then
           if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. (.not. wetsfchl(mgs))  ) then
            chlmul1(mgs) = (ft*ex1*chlacw(mgs))
            qhlmul1(mgs) = cimas0*chlmul1(mgs)*rhoinv(mgs)
           endif
         endif
        endif ! itype2

        if ( itype1 > 0 ) then
         if ( qx(mgs,lh) .gt. qxmin(lh) .and. (.not. wetsfc(mgs))  ) then
          tmp = ft*(3.5e+08)*rho0(mgs)*qhacw(mgs)
          chmul1(mgs) = chmul1(mgs) + tmp
          qhmul1(mgs) = qhmul1(mgs) + cimas0*tmp*rhoinv(mgs)
         endif
         if ( lhl .gt. 1 ) then
           if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. (.not. wetsfchl(mgs)) ) then
            tmp = ft*(3.5e+08)*rho0(mgs)*qhlacw(mgs)
            chlmul1(mgs) = chlmul1(mgs) + tmp
            qhlmul1(mgs) = qhlmul1(mgs) + cimas0*tmp*rhoinv(mgs)
           endif
         endif
        endif ! itype1

        
        endif ! ft

        endif ! xv(mgs,lc) .gt. 0.0 .and.

       else ! ipconc .lt. 2
!
!  define the temperature function
!
      fimt1(mgs) = 0.0
!
! Cotton et al. (1986) version
!
      if ( temg(mgs) .ge. 268.15 .and. temg(mgs) .le. 270.15 ) then
        fimt1(mgs) = 1.0 -(temg(mgs)-268.15)/2.0
      elseif (temg(mgs) .le. 268.15 .and. temg(mgs) .ge. 265.15 ) then
        fimt1(mgs) = 1.0 +(temg(mgs)-268.15)/3.0
      else 
        fimt1(mgs) = 0.0
      end if
!
! Ferrier (1994) version
!
      if ( temg(mgs) .ge. 265.15 .and. temg(mgs) .le. 267.15 ) then
        fimt1(mgs) = 0.5
      elseif (temg(mgs) .ge. 267.15 .and. temg(mgs) .le. 269.15 ) then
        fimt1(mgs) = 1.0
      elseif (temg(mgs) .ge. 269.15 .and. temg(mgs) .le. 271.15 ) then
        fimt1(mgs) = 0.5
      else 
        fimt1(mgs) = 0.0
      end if
!
!
!   type I:  350 splinters are formed for every 1e-3 grams of cloud
!            water accreted by graupel/hail (note converted to MKS units)
!            3.5e+8 has units of 1/kg
!
      if ( itype1 .ge. 1 ) then
       fimta(mgs) = (3.5e+08)*rho0(mgs)
      else
       fimta(mgs) = 0.0
      endif

!
!
!   type II:  1 splinter formed for every 250 cloud droplets larger than
!             24 micons in diameter (12 microns in radius) accreted by
!             graupel/hail
!
!
      fimt2(mgs) = 0.0
      xcwmas = xmas(mgs,lc) * 1000.
!
      if ( itype2 .ge. 1 ) then
      if ( xcwmas.lt.1.26e-9 ) then
        fimt2(mgs) = 0.0
      end if
      if ( xcwmas .le. 3.55e-9 .and. xcwmas .ge. 1.26e-9 ) then
        fimt2(mgs) = (2.27)*alog(xcwmas) + 13.39
      end if
      if ( xcwmas .gt. 3.55e-9 ) then
        fimt2(mgs) = 1.0
      end if

      fimt2(mgs) = min(fimt2(mgs),1.0)
      fimt2(mgs) = max(fimt2(mgs),0.0)
      
      endif
!
!     qhmul2 = 0.0
!     qsmul2 = 0.0
!
!     qhmul2 =
!    >  (4.0e-03)*fimt1(mgs)*fimt2(mgs)*qhacw(mgs)
!     qsmul2 =
!    >  (4.0e-03)*fimt1(mgs)*fimt2(mgs)*qsacw(mgs)
!
!      cimas0 = (1.0e-12)
!      cimas0 = 2.5e-10
      if ( .not. wetsfc(mgs) ) then
      chmul1(mgs) =  fimt1(mgs)*(fimta(mgs) +   &
     &                           (4.0e-03)*fimt2(mgs))*qhacw(mgs)
      endif
!
      qhmul1(mgs) =  chmul1(mgs)*(cimas0/rho0(mgs))

         if ( lhl .gt. 1 ) then
           if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. (.not. wetsfchl(mgs)) ) then
            tmp = fimt1(mgs)*(fimta(mgs) +   &
     &                           (4.0e-03)*fimt2(mgs))*qhlacw(mgs)
            chlmul1(mgs) =  tmp
            qhlmul1(mgs) = cimas0*tmp*rhoinv(mgs)
           endif
         endif

!      qsmul1(mgs) =  csmul1(mgs)*(cimas0/rho0(mgs))
!
      endif ! ( ipconc .ge. 2 )
      
      end if ! (in temperature range)
      
      endif ! ( itype1 .eq. 1 .or. itype2 .eq. 1)
!
      end do
!
!
!
!     end if
!
!     end do
!
!
! ICE MULTIPLICATION FROM SNOW
!   Lo and Passarelli 82 / Willis and Heymsfield 89 / Schuur and Rutledge 00b
!   using kfrag as fragmentation rate (s-1) / 500 microns as char mean diam for max snow mix ratio
!
      csmul(:) = 0.0
      qsmul(:) = 0.0
      
      if ( isnwfrac /= 0 ) then
      do mgs = 1,ngscnt
       if (temg(mgs) .gt. 265.0) then !{
        if (xdia(mgs,ls,1) .gt. 100.e-6 .and. xdia(mgs,ls,1) .lt. 2.0e-3) then  ! equiv diameter 100microns to 2mm

        tmp = rhoinv(mgs)*pi*xdn(mgs,ls)*cx(mgs,ls)*(500.e-6)**3
        qsmul(mgs) = Max( kfrag*( qx(mgs,ls) - tmp ) , 0.0 )

        qsmul(mgs) = Min( qxmxd(mgs,li), qsmul(mgs) )
        csmul(mgs) = Min( cxmxd(mgs,li), rho0(mgs)*qsmul(mgs)/mfrag )

        endif
       endif !}
      enddo
      endif

!
!  frozen rain-rain interaction....
!
!
!
!
!  rain-ice interaction
!
!
      do mgs = 1,ngscnt
      qracif(mgs) = qraci(mgs)
      cracif(mgs) = craci(mgs)
!      ciacrf(mgs) = ciacr(mgs)
      end do
!
! 
!  vapor to pristine ice crystals   UP
!
!
!
!  compute the nucleation rate
!
!     do mgs = 1,ngscnt
!     idqis = 0
!     if ( ssi(mgs) .gt. 1.0 ) idqis = 1
!     fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
!     dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/
!    >  (1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
!     qidsvp(mgs) = dqisdt(mgs)
!     cnnt = min(cnit*exp(-temcg(mgs)*bta1),1.0e+09)
!     qiint(mgs) = 
!    >  il5(mgs)*idqis*(1.0*dtpinv)
!    <  *min((6.88e-13)*cnnt/rho0(mgs), 0.25*dqisdt(mgs)) 
!     end do
!
!  Meyers et al. (1992; JAS) and Ferrier (1994) primary ice nucleation
!
      cmassin = cimasn  ! 6.88e-13
      do mgs = 1,ngscnt
      qiint(mgs) = 0.0
      ciint(mgs) = 0.0
      qicicnt(mgs) = 0.0
      cicint(mgs) = 0.0
      qipipnt(mgs) = 0.0
      cipint(mgs) = 0.0
      ccitmp = 0.0
      if ( icenucopt == 1 .or. icenucopt == -10 .or. icenucopt == -11 ) then
      if ( ( temg(mgs) .lt. 268.15 .or.  &
!     : ( imeyers5 .and. temg(mgs) .lt.  273.0) ) .and.    &
     & ( imeyers5 .and. temg(mgs) .lt.  272.0 .and. temgkm2(mgs) .lt. tfr) ) .and.    &
     &    ciintmx .gt. (cx(mgs,li)+ccitmp)  &
!     :    .and. cninm(mgs) .gt. 0.   &
     &     ) then
      fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
      dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/   &
     &  (1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
!      qidsvp(mgs) = dqisdt(mgs)
      idqis = 0
      if ( ssi(mgs) .gt. 1.0 ) then
      idqis = 1 
      dzfacp = max( float(kgsp(mgs)-kgs(mgs)), 0.0 )
      dzfacm = max( float(kgs(mgs)-kgsm(mgs)), 0.0 )
      qiint(mgs) =   &
     &  idqis*il5(mgs)   &
     &  *(cmassin/rho0(mgs))   &
     &  *max(0.0,wvel(mgs))   &
     &  *max((cninp(mgs)-cninm(mgs)),0.0)/gz(igs(mgs),jgs,kgs(mgs))   &
     &  /((dzfacp+dzfacm))

      qiint(mgs) = min(qiint(mgs), max(0.25*dqisdt(mgs),0.0)) 
      ciint(mgs) = qiint(mgs)*rho0(mgs)/cmassin
      
!
! limit new crystals so it does not increase the current concentration
!  above ciintmx 20,000 per liter (2.e7 per m**3)
!
!      ciintmx = 1.e9
!      ciintmx = 1.e9
      if ( icenucopt /= -10 ) then
      
        if ( lcin > 1 ) then
          ciint(mgs) = Min(ciint(mgs), ccin(mgs)*dtpinv) ! because ciint is a *rate*
          ccin(mgs) = ccin(mgs) - ciint(mgs)*dtp
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
        elseif ( lcina > 1 ) then
          ciint(mgs) = Max(0.0, Min( ciint(mgs), Min( cnina(mgs), ciintmx ) - cina(mgs) ))
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
      
        elseif ( icenucopt == 1 .and. ciint(mgs) .gt. Max(0.0, ciintmx - cx(mgs,li) - ccitmp )*dtpinv  ) then
          ciint(mgs) = Max(0.0, ciintmx - (cx(mgs,li)) )*dtpinv 
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)

        elseif ( icenucopt == -11 .and. dtp*ciint(mgs) .gt. ( cnina(mgs) - (cx(mgs,li) - ccitmp))) then
          ciint(mgs) = Max(0.0,  cnina(mgs) - (cx(mgs,li)+ccitmp)*dtpinv )
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)

        endif
      endif
      
      end if
      endif

      elseif ( icenucopt == 2 .or. icenucopt == -1 .or. icenucopt == -2 ) then
      
        if ( ( temg(mgs) .lt. 268.15 .and. ssw(mgs) > 1.0 ) .or. ssi(mgs) > 1.25 ) then
          if ( lcin > 1 ) then
           ciint(mgs) = Min(cnina(mgs), ccin(mgs))
           ciint(mgs) = Min( ciint(mgs), Max(0.0, ciintmx - (cx(mgs,li) - ccitmp) ) ) ! do not initiate ice beyond concentration of ciintmx
           ccin(mgs) = ccin(mgs) - ciint(mgs)
           ciint(mgs) = ciint(mgs)*dtpinv ! convert total initiation to a rate
          else
           ciint(mgs) = Max( 0.0, cnina(mgs) - cina(mgs) )*dtpinv
          endif
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)

          fiinit(mgs) = (felv(mgs)**2)/(cp*rw)
          dqisdt(mgs) = (qx(mgs,lv)-qis(mgs))/(1.0 + fiinit(mgs)*qis(mgs)/tsqr(mgs))
          qiint(mgs) = min(qiint(mgs), max(0.25*dqisdt(mgs),0.0))
          ciint(mgs) = qiint(mgs)*rho0(mgs)/cmassin
        endif
      
      
      
      elseif ( icenucopt == 3 .or. icenucopt == 4 .or. icenucopt == 10 ) then
        if (  temg(mgs) .lt. 268.15 ) then
          if ( lcin > 1 ) then
           ciint(mgs) = Min(cnina(mgs), ccin(mgs))
           ciint(mgs) = Min( ciint(mgs), Max(0.0, ciintmx - (cx(mgs,li) + ccitmp) ) ) ! do not initiate ice beyond concentration of ciintmx
           ccin(mgs) = ccin(mgs) - ciint(mgs)
           ciint(mgs) = ciint(mgs)*dtpinv ! convert total initiation to a rate
          else
           ciint(mgs) = Max( 0.0, cnina(mgs) - cina(mgs) )*dtpinv
          endif
          qiint(mgs) = ciint(mgs)*cmassin/rho0(mgs)
        endif

      endif
!
      if ( xplate(mgs) .eq. 1 ) then
      qipipnt(mgs) = qiint(mgs)
      cipint(mgs) = ciint(mgs)
      end if
!
      if ( xcolmn(mgs) .eq. 1 ) then
      qicicnt(mgs) = qiint(mgs)
      cicint(mgs) = ciint(mgs)
      end if
!
!     qipipnt(mgs) = 0.0
!     qicicnt(mgs) = qiint(mgs)
!
      end do
!
! 

!
!  vapor to cloud droplets   UP
!
      if (ndebug .gt. 0 ) write(0,*) 'dbg = 8'
!
!
      if (ndebug .gt. 0 ) write(0,*) 'Collection: set 3-component'
!
!  time for riming....
!
!     rimtim = 240.0
!     dtrim = rimtim
!     xacrtim  = 120.0
!     tranfr = 0.50
!     tranfw = 0.50
!
!  coefficients for riming
!
!     rimc1 = 300.00
!     rimc2 = 0.44
!
! 
!  zero some arrays
!
!
      do mgs = 1,ngscnt
      qrshr(mgs) = 0.0
      qwshw(mgs) = 0.0
      cwshw(mgs) = 0.0
      qsshrp(mgs) = 0.0
      qhshrp(mgs) = 0.0
      end do
!
!
!  first sum all of the shed rain
!
!
      do mgs = 1,ngscnt
      qrshr(mgs) = qsshr(mgs) + qhshr(mgs) + qhlshr(mgs)
      crshr(mgs) = chshrr(mgs)/rzxh(mgs) + chlshrr(mgs)/rzxhl(mgs)
      
      
      if ( ipconc .ge. 3 ) then
!       crshr(mgs) = Max(crshr(mgs), rho0(mgs)*qrshr(mgs)/(xdn(mgs,lr)*vr1mm) )
      endif
      end do 
!
!
!

!
!
!
!
      if ( ipconc .ge. 1 ) then
!
!
!  concentration production terms
!
!  YYY
!
!
!       DO mgs = 1,ngscnt
       pccwi(:) = 0.0
       pccwd(:) = 0.0
       pccwdacc(:) = 0.0
       pccii(:) = 0.0
       pccin(:) = 0.0
       pccid(:) = 0.0
       pcisi(:) = 0.0
       pcisd(:) = 0.0
       pcrwi(:) = 0.0
       pcrwd(:) = 0.0
       pcswi(:) = 0.0
       pcswd(:) = 0.0
       pchwi(:) = 0.0
       pchwd(:) = 0.0
       pchli(:) = 0.0
       pchld(:) = 0.0
!       ENDDO
!
!  Cloud ice
!
!      IF ( ipconc .ge. 1 ) THEN

      if ( warmonly < 0.5 ) then
      if ( ffrzs < 1.0 ) then
      do mgs = 1,ngscnt
      pccii(mgs) =   &
     &   il5(mgs)*cicint(mgs)  &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*cwfrzc(mgs)+cwctfzc(mgs)   &
     &  +cicichr(mgs))   &
     &  +chmul1(mgs)   &
     &  +chlmul1(mgs)    &
     &  + csplinter(mgs) + csplinter2(mgs)   &
     &  +csmul(mgs)
     
       pccii(mgs) = pccii(mgs)*(1.0 - ffrzs)
       
!     >  + nsplinter*(crfrzf(mgs) + crfrz(mgs))
      pccid(mgs) =   &
     &   il5(mgs)*(-cscni(mgs) - cscnvi(mgs) & ! - cwaci(mgs)   &
     &  -craci(mgs)    &
     &  -csaci(mgs)   &
     &  -chaci(mgs) - chlaci(mgs)   &
     &  -chcni(mgs))   &
     &  +il5(mgs)*cisbv(mgs)   &
     &  -(1.-il5(mgs))*cimlr(mgs)

      pccin(mgs) = ciint(mgs)
      

      end do
      endif ! ffrzs
      elseif ( warmonly < 0.8 ) then
      do mgs = 1,ngscnt
      
!      qiint(mgs) = 0.0
!      cicint(mgs) = 0.0
!      qicicnt(mgs) = 0.0
      
      pccii(mgs) =   &
     &   il5(mgs)*cicint(mgs)   &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*cwfrzc(mgs)+cwctfzc(mgs)   &
     &  +cicichr(mgs))   &
     &  +chmul1(mgs)   &
     &  +chlmul1(mgs)    &
     &  + csplinter(mgs) + csplinter2(mgs)   &
     &  +csmul(mgs)
     
       pccii(mgs) = pccii(mgs)*(1. - ffrzs)
      pccid(mgs) =   &
!     &   il5(mgs)*(-cscni(mgs) - cscnvi(mgs) & ! - cwaci(mgs)   &
!     &  -craci(mgs)    &
!     &  -csaci(mgs)   &
!     &  -chaci(mgs) - chlaci(mgs)   &
!     &  -chcni(mgs))   &
     &  +il5(mgs)*cisbv(mgs)   &
     &  -(1.-il5(mgs))*cimlr(mgs)

      pccin(mgs) = ciint(mgs)

      end do
      endif ! warmonly

      
!      ENDIF ! ( ipconc .ge. 1 )
!
!  Cloud water
!
      if ( ipconc .ge. 2 ) then
      
      do mgs = 1,ngscnt
      pccwi(mgs) =  (0.0) - cwshw(mgs) ! + (1-il5(mgs))*(-cirmlw(mgs))
      
      if ( warmonly < 0.5 ) then
      pccwd(mgs) =    &
     &  - cautn(mgs) +   &
     &  il5(mgs)*(-ciacw(mgs)-cwfrz(mgs)-cwctfzp(mgs)   &
     &  -cwctfzc(mgs)   &
     &   )   &
     &  -cracw(mgs) -csacw(mgs)  -chacw(mgs) - chlacw(mgs)


      elseif ( warmonly < 0.8 ) then
      pccwd(mgs) =    &
     &  - cautn(mgs) +   &
     &  il5(mgs)*(  &
     & -ciacw(mgs)-cwfrz(mgs)-cwctfzp(mgs)   &
     &  -cwctfzc(mgs)   &
     &   )   &
     &  -cracw(mgs) -chacw(mgs) -chlacw(mgs) 
      else
      
!       tmp3d(igs(mgs),jy,kgs(mgs)) = crcnw(mgs)

!       cracw(mgs) = 0.0 ! turn off accretion
!       qracw(mgs) = 0.0
!       crcev(mgs) = 0.0 ! turn off evap
!       qrcev(mgs) = 0.0 ! turn off evap
!       cracr(mgs) = 0.0 ! turn off self collection
       
       
!       cautn(mgs) = 0.0 
!       crcnw(mgs) = 0.0
!       qrcnw(mgs) = 0.0

      pccwd(mgs) =    &
     &  - cautn(mgs) -cracw(mgs)
      endif


      if ( .false. .and. exwmindiam > 0.0 .and. ccwresv(mgs) > 0.0 ) then
      pccwdacc(mgs) =    &
     &  il5(mgs)*(-ciacw(mgs)  &
     &   )   &
     &  -cracw(mgs) -csacw(mgs)  -chacw(mgs) - chlacw(mgs)

      if ( -pccwdacc(mgs)*dtp .gt. cx(mgs,lc) - ccwresv(mgs) ) then

       frac = -(cx(mgs,lc) - ccwresv(mgs) )/(pccwdacc(mgs)*dtp)
       pccwdacc(mgs) = -(cx(mgs,lc) - ccwresv(mgs) )*dtpinv

        ciacw(mgs)   = frac*ciacw(mgs)
        cracw(mgs)   = frac*cracw(mgs)
        csacw(mgs)   = frac*csacw(mgs)
        chacw(mgs)   = frac*chacw(mgs)
        cautn(mgs)   = frac*cautn(mgs)
       
        if ( lhl .gt. 1 ) chlacw(mgs)   = frac*chlacw(mgs)

! resum
      pccwd(mgs) =    &
     &  - cautn(mgs) +   &
     &  il5(mgs)*(-ciacw(mgs)-cwfrzp(mgs)-cwctfzp(mgs)   &
     &  -cwfrzc(mgs)-cwctfzc(mgs)   &
     &  -il5(mgs)*(ciihr(mgs))   &
     &   )   &
     &  -cracw(mgs) -csacw(mgs)  -chacw(mgs) - chlacw(mgs)

      endif

      endif


      if ( -pccwd(mgs)*dtp .gt. cx(mgs,lc) ) then
!       write(0,*) 'OUCH! pccwd(mgs)*dtp .gt. ccw(mgs) ',pccwd(mgs),cx(mgs,lc)
!       write(0,*) 'qc = ',qx(mgs,lc)
!       write(0,*) -ciacw(mgs)-cwfrzp(mgs)-cwctfzp(mgs)-cwfrzc(mgs)-cwctfzc(mgs)
!       write(0,*)  -cracw(mgs) -csacw(mgs)  -chacw(mgs)
!       write(0,*) - cautn(mgs)

       frac = -cx(mgs,lc)/(pccwd(mgs)*dtp)
       pccwd(mgs) = -cx(mgs,lc)*dtpinv

        ciacw(mgs)   = frac*ciacw(mgs)
        cwfrz(mgs)  = frac*cwfrz(mgs)
        cwfrzp(mgs)  = frac*cwfrzp(mgs)
        cwctfzp(mgs) = frac*cwctfzp(mgs)
        cwfrzc(mgs)  = frac*cwfrzc(mgs)
        cwctfzc(mgs) = frac*cwctfzc(mgs)
        cwctfz(mgs) = frac*cwctfz(mgs)
        cracw(mgs)   = frac*cracw(mgs)
        csacw(mgs)   = frac*csacw(mgs)
        chacw(mgs)   = frac*chacw(mgs)
        cautn(mgs)   = frac*cautn(mgs)
       
        pccii(mgs) = pccii(mgs) - (1.-frac)*il5(mgs)*(cwfrzc(mgs)+cwctfzc(mgs))*(1. - ffrzs)
        if ( lhl .gt. 1 ) chlacw(mgs)   = frac*chlacw(mgs)

!       STOP
      endif

      end do

      endif ! ipconc

!
!  Rain
!
      if ( ipconc .ge. 3 ) then

      do mgs = 1,ngscnt

      if ( warmonly < 0.5 ) then
      pcrwi(mgs) = &
!     >   cracw(mgs) +    &
     &   crcnw(mgs)   &
     &  +(1-il5(mgs))*(   &
     &    -chmlrr(mgs)/rzxh(mgs)   &
     &    -chlmlrr(mgs)/rzxhl(mgs)   &
!     &    -csmlr(mgs)/rzxs(mgs)     &
     &   -csmlrr(mgs)     &
     &   - cimlr(mgs) )   &
     &  -crshr(mgs)             !null at this point when wet snow/graupel included
      pcrwd(mgs) =   &
     &   il5(mgs)*(-ciacr(mgs) - crfrz(mgs) ) & ! - cipacr(mgs))
!     >  -csacr(mgs)   &
     &  - chacr(mgs) - chlacr(mgs)   &
     &  +crcev(mgs)   &
     &  - cracr(mgs)
!     >  -il5(mgs)*ciracr(mgs)


      elseif ( warmonly < 0.8 ) then
       pcrwi(mgs) = &
     &   crcnw(mgs)   &
     &  +(1-il5(mgs))*(   &
     &    -chmlrr(mgs)/rzxh(mgs)    &
     &    -chlmlrr(mgs)/rzxhl(mgs)   &
!     &    -csmlr(mgs)     &
     &   -csmlrr(mgs)     &
     &   - cimlr(mgs) )   &
     &  -crshr(mgs)             !null at this point when wet snow/graupel included
      pcrwd(mgs) =   &
     &   il5(mgs)*( - crfrz(mgs) ) & ! - cipacr(mgs))
     &  - chacr(mgs)    &
     &  - chlacr(mgs)    &
     &  +crcev(mgs)   &
     &  - cracr(mgs)
      else
      pcrwi(mgs) =   &
     &   crcnw(mgs)
      pcrwd(mgs) =   &
     &  +crcev(mgs)   &
     &  - cracr(mgs)

!        tmp3d(igs(mgs),jy,kgs(mgs)) = vtxbar(mgs,lr,1) ! crcnw(mgs) ! (pcrwi(mgs) + pcrwd(mgs))
!        pcrwi(mgs) = 0.0
!        pcrwd(mgs) = 0.0
!        qrcnw(mgs) = 0.0

      endif


      frac = 0.0
      if ( -pcrwd(mgs)*dtp .gt. cx(mgs,lr) ) then
!       write(0,*) 'OUCH! pcrwd(mgs)*dtp .gt. crw(mgs) ',pcrwd(mgs)*dtp,cx(mgs,lr),mgs,igs(mgs),kgs(mgs)
!       write(0,*) -ciacr(mgs)
!       write(0,*) -crfrz(mgs)
!       write(0,*) -chacr(mgs)
!       write(0,*)  crcev(mgs)
!       write(0,*)  -cracr(mgs)

       frac =  -cx(mgs,lr)/(pcrwd(mgs)*dtp)
       pcrwd(mgs) = -cx(mgs,lr)*dtpinv

        ciacr(mgs) = frac*ciacr(mgs)
        ciacrf(mgs) = frac*ciacrf(mgs)
        ciacrs(mgs) = frac*ciacrs(mgs)
        crfrz(mgs) = frac*crfrz(mgs)
        crfrzf(mgs) = frac*crfrzf(mgs)
        crfrzs(mgs) = frac*crfrzs(mgs)
        chacr(mgs) = frac*chacr(mgs)
        chlacr(mgs) = frac*chlacr(mgs)
        crcev(mgs) = frac*crcev(mgs)
        cracr(mgs) = frac*cracr(mgs)

!       STOP
      endif

      end do

      endif


      if ( warmonly < 0.5 ) then

!
!  Snow
!
      if ( ipconc .ge. 4 ) then !

      do mgs = 1,ngscnt
      pcswi(mgs) =   &
     &   il5(mgs)*(cscnis(mgs) + cscnvis(mgs) )    &
     &  + cwfrz2snowfrac*cwfrz(mgs)/cwfrz2snowratio  &
     &  + cscnh(mgs)
      
      if (  ffrzs > 0.0 ) then
       pcswi(mgs) =  pcswi(mgs) + ffrzs* (  &
     &   il5(mgs)*cicint(mgs)   &
     &  +il5(mgs)*(cwfrzc(mgs)+cwctfzc(mgs)   &
     &  +cicichr(mgs))  &
     &  +chmul1(mgs)   &
     &  +chlmul1(mgs)    &
     &  + csplinter(mgs) + csplinter2(mgs)   &
     &  +csmul(mgs) )
      endif

      
      if ( ess0 < 0.0 ) then
         csacs(mgs) = Max(0.0, csacs(mgs) - (ifrzs)*(crfrzs(mgs) + ciacrs(mgs)))
      endif
      
      pcswd(mgs) = &
!     :  cracs(mgs)     &
     &  -chacs(mgs) - chlacs(mgs)   &
     &  -chcns(mgs)   &
     &  +(1-il5(mgs))*csmlr(mgs) + csshr(mgs) & ! + csshrp(mgs)
!     >  +il5(mgs)*(cssbv(mgs))   &
     &   + cssbv(mgs)   &
     &  - csacs(mgs)

      frac = 0.0
      if ( imixedphase == 0 ) then
        if ( cx(mgs,ls) + dtp*(pcswi(mgs) + pcswd(mgs)) < 0.0 ) then
         frac = (-cx(mgs,ls) + pcswi(mgs)*dtp)/(pcswd(mgs)*dtp)
         
           pcswd(mgs) = frac*pcswd(mgs)
           
           chacs(mgs)  = frac*chacs(mgs) 
           chlacs(mgs) = frac*chlacs(mgs)
           chcns(mgs)  = frac*chcns(mgs) 
           csmlr(mgs)  = frac*csmlr(mgs) 
           csshr(mgs)  = frac*csshr(mgs) 
           cssbv(mgs)  = frac*cssbv(mgs) 
           csacs(mgs)  = frac*csacs(mgs)
      
        endif
      endif


      
      pccii(mgs) =  pccii(mgs) &
     &  + (1. - ifrzs)*crfrzs(mgs) &
     &  + (1. - ifrzs)*ciacrs(mgs)

      pcswi(mgs) =  pcswi(mgs) &
     &  + (ifrzs)*crfrzs(mgs) &
     &  + (ifrzs)*ciacrs(mgs)

      end do

      endif

!
!  Graupel
!
      if ( ipconc .ge. 5 ) then !
      do mgs = 1,ngscnt
      pchwi(mgs) =   &
     &  +(ffrzh*ifrzg*crfrzf(mgs)   &
     & +il5(mgs)*ffrzh*ifiacrg*(ciacrf(mgs) ))    &
     & + f2h*chcnsh(mgs) + f2h*chcnih(mgs) + chcnhl(mgs)

      pchwd(mgs) =   &
     &  (1-il5(mgs))*chmlr(mgs) &
!     >  + il5(mgs)*chsbv(mgs)   &
     &  + chsbv(mgs)   &
     &  - il5(mgs)*chlcnh(mgs) &
     &  - cscnh(mgs)

      end do



!

!
!  Hail
!
      if ( lhl .gt. 1 .and. lnhl > 1 ) then !
      do mgs = 1,ngscnt
      pchli(mgs) = (ffrzh*(1.0-ifrzg)*crfrzf(mgs) +il5(mgs)*ffrzh*(1.0-ifiacrg)*(ciacrf(mgs) ))  &
     & + chlcnhhl(mgs) *rzxhlh(mgs)

      pchld(mgs) =   &
     &  (1-il5(mgs))*chlmlr(mgs)   &
!     >  + il5(mgs)*chlsbv(mgs)   &
     &  + chlsbv(mgs) - chcnhl(mgs)
      
      if ( imixedphase == 0 ) then
      frac = 0.0
      if ( cx(mgs,lhl) + dtp*(pchli(mgs) + pchld(mgs)) < 0.0 ) then
        ! rescale depletion

         frac = (-cx(mgs,lhl) + pchli(mgs)*dtp)/(pchld(mgs)*dtp)
         
         chlmlr(mgs) = frac*chlmlr(mgs)
         chlsbv(mgs) = frac*chlsbv(mgs)
         chcnhl(mgs) = frac*chcnhl(mgs)
           
         pchld(mgs) = frac*pchld(mgs)
           
      endif
      endif

      end do
      
      endif
!

      endif ! (ipconc .ge. 5 )

      elseif ( warmonly < 0.8 ) then

!
!  Graupel
!
      if ( ipconc .ge. 5 ) then !
      do mgs = 1,ngscnt
      pchwi(mgs) =   &
     &  +ifrzg*(crfrzf(mgs) ) ! +il5(mgs)*(ciacrf(mgs) ))

      pchwd(mgs) =   &
     &  (1-il5(mgs))*chmlr(mgs) &
     &  - il5(mgs)*chlcnh(mgs)
      end do
!
!  Hail
!
      if ( lhl .gt. 1 ) then !
      do mgs = 1,ngscnt
      pchli(mgs) = (1.0-ifrzg)*(crfrzf(mgs)) & ! +il5(mgs)*(ciacrf(mgs) ))  &
     & + chlcnhhl(mgs) *rzxhl(mgs)/rzxh(mgs)

      pchld(mgs) =   &
     &  (1-il5(mgs))*chlmlr(mgs) !  &
!     >  + il5(mgs)*chlsbv(mgs)   &
!     &  + chlsbv(mgs)

!      IF ( pchli(mgs) .ne. 0. .or. pchld(mgs) .ne. 0 ) THEN
!       write(0,*) 'dr: pchli,pchld = ', pchli(mgs),pchld(mgs), igs(mgs),kgs(mgs)
!      ENDIF
      end do

      endif

      endif ! ipconc >= 5

      endif ! warmonly

!

!
!  Balance and checks for continuity.....within machine precision...
!
      do mgs = 1,ngscnt
      pctot(mgs)   = pccwi(mgs) +pccwd(mgs) +   &
     &               pccii(mgs) +pccid(mgs) +   &
     &               pcrwi(mgs) +pcrwd(mgs) +   &
     &               pcswi(mgs) +pcswd(mgs) +   &
     &               pchwi(mgs) +pchwd(mgs) +   &
     &               pchli(mgs) +pchld(mgs)
      end do
!
!
      endif ! ( ipconc .ge. 1 )
!
!
!
!
!
!  GOGO
!  production terms for mass
!
!
       pqwvi(:) = 0.0
       pqwvd(:) = 0.0
       pqcwi(:) = 0.0
       pqcwd(:) = 0.0
       pqcwdacc(:) = 0.0
       pqcii(:) = 0.0
       pqcid(:) = 0.0
       pqrwi(:) = 0.0
       pqrwd(:) = 0.0
       pqswi(:) = 0.0
       pqswd(:) = 0.0
       pqhwi(:) = 0.0
       pqhwd(:) = 0.0
       pqhli(:) = 0.0
       pqhld(:) = 0.0
       pqlwsi(:) = 0.0
       pqlwsd(:) = 0.0
       pqlwhi(:) = 0.0
       pqlwhd(:) = 0.0
       pqlwlghi(:) = 0.0
       pqlwlghd(:) = 0.0
       pqlwlghli(:) = 0.0
       pqlwlghld(:) = 0.0
       pqlwhli(:) = 0.0
       pqlwhld(:) = 0.0
       if ( ipconc > 5 ) then
       pzhwi(:) = 0.0
       pzhwd(:) = 0.0
       pzrwi(:) = 0.0
       pzrwd(:) = 0.0
       pzhli(:) = 0.0
       pzhld(:) = 0.0
       endif


!
!  Vapor
!
      if ( warmonly < 0.5 ) then
      do mgs = 1,ngscnt
      
! NOTE: ANY CHANGES HERE ALSO NEED TO GO INTO THE RESUM FARTHER DOWN!
      pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs))   &
     &  -Min(0.0, qhcev(mgs))   &
     &  -Min(0.0, qfcev(mgs))   &
     &  -Min(0.0, qhlcev(mgs))   &
     &  -Min(0.0, qscev(mgs))   &
!     >  +il5(mgs)*(-qhsbv(mgs) - qhlsbv(mgs) )   &
     &  -qhsbv(mgs) - qhlsbv(mgs)   &
     &  -qssbv(mgs)    &
     &  -il5(mgs)*qisbv(mgs)
      
      pqwvd(mgs) =     &
     &  -Max(0.0, qrcev(mgs))   &
     &  -Max(0.0, qhcev(mgs))   &
     &  -Max(0.0, qfcev(mgs))   &
     &  -Max(0.0, qhlcev(mgs))   &
     &  -Max(0.0, qscev(mgs))   &
     &  +il5(mgs)*(-qiint(mgs)   &
     &  -qhdpv(mgs) -qsdpv(mgs) - qhldpv(mgs))   &
     &  -il5(mgs)*qidpv(mgs)  
      
      end do

      elseif ( warmonly < 0.8 ) then
      do mgs = 1,ngscnt
      pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs)) &
     &  -il5(mgs)*qisbv(mgs)
      pqwvd(mgs) =     &
     &  +il5(mgs)*(-qiint(mgs)   &
!     &  -qhdpv(mgs) ) & !- qhldpv(mgs))   &
     &  -qhdpv(mgs) - qhldpv(mgs))   &
!     &  -qhdpv(mgs) -qsdpv(mgs) - qhldpv(mgs))   &
     &  -Max(0.0, qrcev(mgs))     &
     &  -il5(mgs)*qidpv(mgs)  
      end do

      else
      do mgs = 1,ngscnt
      pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs))
      pqwvd(mgs) =     &
     &  -Max(0.0, qrcev(mgs))
      end do

      endif ! warmonly
!
!  Cloud water
!
      do mgs = 1,ngscnt

      pqcwi(mgs) =  (0.0) + qwcnr(mgs) - qwshw(mgs)

      if ( warmonly < 0.5 ) then
      pqcwd(mgs) =    &
     &  il5(mgs)*(-qiacw(mgs)-qwfrz(mgs)-qwctfz(mgs))   &
     &  -il5(mgs)*(qiihr(mgs))   &
     &  -qracw(mgs) -qsacw(mgs) -qrcnw(mgs) -qhacw(mgs) - qhlacw(mgs)  !&
!     &  -il5(mgs)*(qwfrzp(mgs))
      elseif ( warmonly < 0.8 ) then
      pqcwd(mgs) =    &
     &  il5(mgs)*(-qiacw(mgs)-qwfrz(mgs)-qwctfz(mgs))   &
     &  -il5(mgs)*(qiihr(mgs))   &
     &  -qracw(mgs) -qrcnw(mgs) -qhacw(mgs) -qhlacw(mgs)
      else
      pqcwd(mgs) =    &
     &  -qracw(mgs) - qrcnw(mgs)
      endif


      if ( pqcwd(mgs) .lt. 0.0 .and. -pqcwd(mgs)*dtp .gt. qx(mgs,lc) ) then

       frac = -Max(0.0,qx(mgs,lc))/(pqcwd(mgs)*dtp)
       pqcwd(mgs) = -qx(mgs,lc)*dtpinv

        qiacw(mgs)   = frac*qiacw(mgs)
!        qwfrzp(mgs)  = frac*qwfrzp(mgs)
!        qwctfzp(mgs) = frac*qwctfzp(mgs)
        qwfrzc(mgs)  = frac*qwfrzc(mgs)
        qwfrz(mgs)  = frac*qwfrz(mgs)
        qwctfzc(mgs) = frac*qwctfzc(mgs)
        qwctfz(mgs) = frac*qwctfz(mgs)
        qracw(mgs)   = frac*qracw(mgs)
        qsacw(mgs)   = frac*qsacw(mgs)
        qhacw(mgs)   = frac*qhacw(mgs)
        vhacw(mgs)   = frac*vhacw(mgs)
        qrcnw(mgs)   = frac*qrcnw(mgs)
        qwfrzp(mgs)  = frac*qwfrzp(mgs)
        if ( lhl .gt. 1 ) then
          qhlacw(mgs)   = frac*qhlacw(mgs)
          vhlacw(mgs)   = frac*vhlacw(mgs)
        endif
!        IF ( lzh .gt. 1 ) zhacw(mgs) = frac*zhacw(mgs)

!       STOP
      endif
      

      end do
!
!  Cloud ice
!
      if ( warmonly < 0.5 ) then

      do mgs = 1,ngscnt
      if ( ffrzs < 1.0 ) then
      pqcii(mgs) =     &
     &   il5(mgs)*qicicnt(mgs)    &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*qwfrzc(mgs)+qwctfzc(mgs))   &
     &  +il5(mgs)*(qicichr(mgs))  &
     &  +qsmul(mgs)               &
     &  +qhmul1(mgs) + qhlmul1(mgs)   &
     & + qsplinter(mgs) + qsplinter2(mgs)
!     > + cimas0*nsplinter*(crfrzf(mgs) + crfrz(mgs))/rho0(mgs)
      endif
       
       pqcii(mgs) = pqcii(mgs)*(1.0 - ffrzs) &
     &  +il5(mgs)*qidpv(mgs)    &
     &  +il5(mgs)*qiacw(mgs)
       
      pqcid(mgs) =     &
     &   il5(mgs)*(-qscni(mgs) - qscnvi(mgs)    & ! -qwaci(mgs)    &
     &  -qraci(mgs)    &
     &  -qsaci(mgs) )   &
     &  -qhaci(mgs)   &
     &  -qhlaci(mgs)    &
     &  +il5(mgs)*qisbv(mgs)    &
     &  +(1.-il5(mgs))*qimlr(mgs)   &
     &  - qhcni(mgs)
      end do

      
      elseif ( warmonly < 0.8 ) then

      do mgs = 1,ngscnt
      pqcii(mgs) =     &
     &   il5(mgs)*qicicnt(mgs)*(1. - ffrzs)    &
     &  +il5(mgs)*((1.0-cwfrz2snowfrac)*qwfrzc(mgs)+qwctfzc(mgs))*(1. - ffrzs)   &
     &  +il5(mgs)*(qicichr(mgs))*(1. - ffrzs)   &
!     &  +il5(mgs)*(qicichr(mgs))   &
!     &  +qsmul(mgs)               &
     &  +qhmul1(mgs) + qhlmul1(mgs)   &
     & + qsplinter(mgs) + qsplinter2(mgs) &
     &  +il5(mgs)*qidpv(mgs)    &
     &  +il5(mgs)*qiacw(mgs)  ! & ! (qiacwi(mgs)+qwacii(mgs))   &
!     &  +il5(mgs)*(qwfrzc(mgs)+qwctfzc(mgs))   &
!     &  +il5(mgs)*(qicichr(mgs))   &
!     &  +qsmul(mgs)               &
!     &  +qhmul1(mgs) + qhlmul1(mgs)   &
!     & + qsplinter(mgs) + qsplinter2(mgs)

      pqcid(mgs) =     &
!     &   il5(mgs)*(-qscni(mgs) - qscnvi(mgs)    & ! -qwaci(mgs)    &
!     &  -qraci(mgs)    &
!     &  -qsaci(mgs) )   &
!     &  -qhaci(mgs)   &
!     &  -qhlaci(mgs)    &
     &  +il5(mgs)*qisbv(mgs)    &
     &  +(1.-il5(mgs))*qimlr(mgs)  ! &
!     &  - qhcni(mgs)
      end do

      endif
!
!  Rain
!

      do mgs = 1,ngscnt
      if ( warmonly < 0.5 ) then
      pqrwi(mgs) =     &
     &   qracw(mgs) + qrcnw(mgs) + Max(0.0, qrcev(mgs))   &
     &  +(1-il5(mgs))*(   &
     &    -qhmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qsmlr(mgs)  - qhlmlr(mgs)     &
     &    -qimlr(mgs))   &
!     &    -qsshr(mgs)       &                      !null at this point when wet snow/graupel included
!     &    -qhshr(mgs)       &                      !null at this point when wet snow/graupel included
!     &    -qhlshr(mgs)      &
     & - qrshr(mgs)

      pqrwd(mgs) =     &
     &  il5(mgs)*(-qiacr(mgs)-qrfrz(mgs))    &
     &  - qsacr(mgs) - qhacr(mgs) - qhlacr(mgs) - qwcnr(mgs)   &
     &  + Min(0.0,qrcev(mgs))
      elseif ( warmonly < 0.8 ) then
      pqrwi(mgs) =     &
     &   qracw(mgs) + qrcnw(mgs) + Max(0.0, qrcev(mgs))   &
     &  +(1-il5(mgs))*(   &
     &    -qhlmlr(mgs)                 &            !null at this point when wet snow/graupel included
     &    -qhmlr(mgs)  )               &            !null at this point when wet snow/graupel included
     &    -qhshr(mgs)                 &           !null at this point when wet snow/graupel included
     &    -qhlshr(mgs)                            !null at this point when wet snow/graupel included
      pqrwd(mgs) =     &
     &  il5(mgs)*(-qrfrz(mgs))    &
     &   - qhacr(mgs)    &
     &   - qhlacr(mgs)    &
     &  + Min(0.0,qrcev(mgs))
      else
      pqrwi(mgs) =     &
     &   qracw(mgs) + qrcnw(mgs) + Max(0.0, qrcev(mgs))
      pqrwd(mgs) =  Min(0.0,qrcev(mgs))
      endif ! warmonly


 !      IF ( pqrwd(mgs) .lt. 0.0 .and. -(pqrwd(mgs) + pqrwi(mgs))*dtp .gt. qx(mgs,lr)  ) THEN
      if ( pqrwd(mgs) .lt. 0.0 .and. -(pqrwd(mgs) + pqrwi(mgs))*dtp .gt. qx(mgs,lr)  ) then

       frac = (-qx(mgs,lr) + pqrwi(mgs)*dtp)/(pqrwd(mgs)*dtp)
!       pqrwd(mgs) = -qx(mgs,lr)*dtpinv  + pqrwi(mgs)

       pqwvi(mgs) = pqwvi(mgs)    &
     &  + Min(0.0, qrcev(mgs))   &
     &  - frac*Min(0.0, qrcev(mgs))
       pqwvd(mgs) =  pqwvd(mgs)   &
     &  + Max(0.0, qrcev(mgs))   &
     &  - frac*Max(0.0, qrcev(mgs))

       qiacr(mgs)  = frac*qiacr(mgs)
       qiacrf(mgs) = frac*qiacrf(mgs)
       qiacrs(mgs) = frac*qiacrs(mgs)
       viacrf(mgs) = frac*viacrf(mgs)
       qrfrz(mgs)  = frac*qrfrz(mgs) 
       qrfrzs(mgs) = frac*qrfrzs(mgs) 
       qrfrzf(mgs) = frac*qrfrzf(mgs)
       vrfrzf(mgs) = frac*vrfrzf(mgs)
       qsacr(mgs)  = frac*qsacr(mgs)
       qhacr(mgs)  = frac*qhacr(mgs)
       vhacr(mgs)  = frac*vhacr(mgs)
       qrcev(mgs)  = frac*qrcev(mgs)
       qhlacr(mgs) = frac*qhlacr(mgs)
       vhlacr(mgs) = frac*vhlacr(mgs)
       qhcev(mgs)  = frac*qhcev(mgs)
       qhlcev(mgs)  = frac*qhlcev(mgs)


      if ( warmonly < 0.5 ) then
       pqrwd(mgs) =     &
     &  il5(mgs)*(-qiacr(mgs)-qrfrz(mgs) - qsacr(mgs))    &
     &  - qhacr(mgs) - qhlacr(mgs) - qwcnr(mgs)   &
     &  + Min(0.0,qrcev(mgs))
      elseif ( warmonly < 0.8 ) then
      pqrwd(mgs) =     &
     &  il5(mgs)*(-qrfrz(mgs))    &
     &   - qhacr(mgs)    &
     &   - qhlacr(mgs)    &
     &  + Min(0.0,qrcev(mgs))
      else
       pqrwd(mgs) =  Min(0.0,qrcev(mgs))
      endif ! warmonly

!
! Resum for vapor since qrcev has changed
!
      if ( qrcev(mgs) .ne. 0.0 ) then
       pqwvi(mgs) =    &
     &  -Min(0.0, qrcev(mgs))   &
     &  -Min(0.0, qhcev(mgs))   &
     &  -Min(0.0, qhlcev(mgs))   &
     &  -Min(0.0, qscev(mgs))   &
!     >  +il5(mgs)*(-qhsbv(mgs)  - qhlsbv(mgs) )   &
     &  -qhsbv(mgs)  - qhlsbv(mgs)   &
     &  -qssbv(mgs)    &
     &  -il5(mgs)*qisbv(mgs) 
     
       pqwvd(mgs) =     &
     &  -Max(0.0, qrcev(mgs))   &
     &  -Max(0.0, qhcev(mgs))   &
     &  -Max(0.0, qhlcev(mgs))   &
     &  -Max(0.0, qscev(mgs))   &
     &  +il5(mgs)*(-qiint(mgs)   &
     &  -qhdpv(mgs) -qsdpv(mgs) - qhldpv(mgs))   &
     &  -il5(mgs)*qidpv(mgs)  

       endif


!       STOP
      endif


      end do

      if ( warmonly < 0.5 ) then

!
!  Snow
!
      do mgs = 1,ngscnt
      pqswi(mgs) =     &
     &   il5(mgs)*(qscni(mgs)+qsaci(mgs)+qsdpv(mgs)   &
     &   + qscnvi(mgs)                        &
     &   + ifrzs*(qiacrs(mgs) + qrfrzs(mgs))  &
     &   + il5(mgs)*(( qwfrzc(mgs) + qwctfzc(mgs) + qicichr(mgs) )*ffrzs   &
     &   +  (1.0 - ffrzs)*cwfrz2snowfrac*qwfrz(mgs) ) &
     &   + il2(mgs)*qsacr(mgs))   &
     &   + il5(mgs)*qicicnt(mgs)*ffrzs        &
     &   + il3(mgs)*(qiacrf(mgs)+qracif(mgs)) & ! only applies for ipconc <= 3
     &   + Max(0.0, qscev(mgs))   &
     &   + qsacw(mgs) + qscnh(mgs) &
     &  + ffrzs*(qsmul(mgs)               &
     &  +qhmul1(mgs) + qhlmul1(mgs)   &
     & + qsplinter(mgs) + qsplinter2(mgs))
      pqswd(mgs) =    &
!     >  -qfacs(mgs) ! -qwacs(mgs)   &
     &  -qracs(mgs)*(1-il2(mgs)) -qhacs(mgs) - qhlacs(mgs)   &
     &  -qhcns(mgs)   &
     &  +(1-il5(mgs))*qsmlr(mgs) + qsshr(mgs)    &    !null at this point when wet snow included
!     >  +il5(mgs)*(qssbv(mgs))   &
     &  + qssbv(mgs)   &
     &  + Min(0.0, qscev(mgs))  &
     &  -qsmul(mgs)
      
      
      if ( imixedphase == 0 .and. pqswd(mgs) .lt. 0.0  ) then
        if ( qx(mgs,ls) + dtp*(pqswi(mgs) + pqswd(mgs)) < 0.0 ) then
         frac = (-qx(mgs,ls) + pqswi(mgs)*dtp)/(pqswd(mgs)*dtp)
         
           pqswd(mgs) = frac*pqswd(mgs)
           
           qracs(mgs)  = frac*qracs(mgs) ! only used for single moment at this time
           qhacs(mgs)  = frac*qhacs(mgs) 
           qhlacs(mgs) = frac*qhlacs(mgs)
           qhcns(mgs)  = frac*qhcns(mgs) 
           qsmlr(mgs)  = frac*qsmlr(mgs) 
           qsshr(mgs)  = frac*qsshr(mgs) 
           qssbv(mgs)  = frac*qssbv(mgs) 
           qsmul(mgs)  = frac*qsmul(mgs) 
           if ( qscev(mgs) < 0.0 ) qscev(mgs) = frac*qscev(mgs)

        endif
      endif
      
      pqcii(mgs) =  pqcii(mgs) &
     &  + (1. - ifrzs)*qrfrzs(mgs) &
     &  + (1. - ifrzs)*qiacrs(mgs)
      
      end do 
      
!
!  Graupel
!
      do mgs = 1,ngscnt
      pqhwi(mgs) =    &
     &  +il5(mgs)*(ffrzh*ifrzg*qrfrzf(mgs)  + (1-il3(mgs))*ffrzh*ifiacrg*(qiacrf(mgs)+qracif(mgs)))   &
     &  + (1-il2(mgs))*(qracs(mgs) + qsacr(mgs))  & ! only used for ipconc < 3
     &  +il5(mgs)*(qhdpv(mgs))   &
     &  +Max(0.0, qhcev(mgs))   &
     &  +qhacr(mgs)+qhacw(mgs)   &
     &  +qhacs(mgs)+qhaci(mgs)   &
     &  + f2h*qhcns(mgs) + f2h*qhcni(mgs) + qhcnhl(mgs)
      pqhwd(mgs) =     &
     &   qhshr(mgs)                &    !null at this point when wet graupel included
     &  +(1-il5(mgs))*qhmlr(mgs)   &    !null at this point when wet graupel included
!     >  +il5(mgs)*qhsbv(mgs)   &
     &  + qhsbv(mgs)   &
     &  + Min(0.0, qhcev(mgs))   &
     &  -qhmul1(mgs) - qhlcnh(mgs) - qscnh(mgs)  &
     &  - ffrzh*(qsplinter(mgs) + qsplinter2(mgs))
!     > - cimas0*nsplinter*(crfrzf(mgs) + crfrz(mgs))/rho0(mgs)

      end do


!
!  Hail
!
      if ( lhl .gt. 1 ) then

      do mgs = 1,ngscnt
      pqhli(mgs) =    &
     &  +il5(mgs)*(qhldpv(mgs) + ((1.0-ifrzg)*qrfrzf(mgs) + (1.0-ifiacrg)*(qiacrf(mgs)+ qracif(mgs))))   &
     &  +Max(0.0, qhlcev(mgs))   &
     &  +qhlacr(mgs)+qhlacw(mgs)   &
     &  +qhlacs(mgs)+qhlaci(mgs)   &
     &  + qhlcnh(mgs)
      pqhld(mgs) =     &
     &   qhlshr(mgs)    &
     &  +(1-il5(mgs))*qhlmlr(mgs)    &
!     >  +il5(mgs)*qhlsbv(mgs)   &
     &  + qhlsbv(mgs)   &
     &  + Min(0.0, qhlcev(mgs))   &
     &  -qhlmul1(mgs) - qhcnhl(mgs)

      if ( imixedphase == 0 ) then
      frac = 0.0
      if ( qx(mgs,lhl) + dtp*(pqhli(mgs) + pqhld(mgs)) < 0.0 ) then
        ! rescale depletion

         frac = (-qx(mgs,lhl) + pqhli(mgs)*dtp)/(pqhld(mgs)*dtp)
         
         qhlmlr(mgs) = frac*qhlmlr(mgs)
         qhlsbv(mgs) = frac*qhlsbv(mgs)
         qhcnhl(mgs) = frac*qhcnhl(mgs)
         qhlmul1(mgs) = frac*qhlmul1(mgs)
         if ( qhlcev(mgs) < 0.0 ) qhlcev(mgs) = frac*qhlcev(mgs)
           
         pqhld(mgs) = frac*pqhld(mgs)
           
      endif
      endif


      end do
      
      endif ! lhl

      elseif ( warmonly < 0.8 ) then
!
!  Graupel
!
      do mgs = 1,ngscnt
      pqhwi(mgs) =    &
     &  +il5(mgs)*ifrzg*(qrfrzf(mgs) )   &
     &  +il5(mgs)*(qhdpv(mgs))   &
     &  +qhacr(mgs)+qhacw(mgs)   
      pqhwd(mgs) =     &
     &   qhshr(mgs)                &    !null at this point when wet graupel included
     &  - qhlcnh(mgs)   &
     &  - qhmul1(mgs)   &
     &  - qsplinter(mgs) - qsplinter2(mgs) &
     &  +(1-il5(mgs))*qhmlr(mgs)        !null at this point when wet graupel included
       end do

!
!  Hail
!
      if ( lhl .gt. 1 ) then

      do mgs = 1,ngscnt
      pqhli(mgs) =    &
     &  +il5(mgs)*(qhldpv(mgs) ) & ! + (1.0-ifrzg)*(qiacrf(mgs)+qrfrzf(mgs)  + qracif(mgs)))   &
     &  +il5(mgs)*(1.0-ifrzg)*(qrfrzf(mgs) )  &
     &  +qhlacr(mgs)+qhlacw(mgs)   &
!     &  +qhlacs(mgs)+qhlaci(mgs)   &
     &  + qhlcnh(mgs)
      pqhld(mgs) =     &
     &   qhlshr(mgs)    &
     &  +(1-il5(mgs))*qhlmlr(mgs)    &
!     >  +il5(mgs)*qhlsbv(mgs)   &
     &  + qhlsbv(mgs)   &
     &  -qhlmul1(mgs) - qhcnhl(mgs)

      end do

      endif ! lhl

      endif ! warmonly

!
!  Liquid water on snow and graupel 
!

      vhmlr(:) = 0.0
      vhlmlr(:) = 0.0
      vhfzh(:) = 0.0
      vhlfzhl(:) = 0.0

      if ( mixedphase ) then
      else ! set arrays for non-mixedphase graupel
      
!        vhshdr(:) = 0.0
        vhmlr(:) = qhmlr(:) ! not actually volume, but treated as q in rate equation
!        vhsoak(:) = 0.0

!        vhlshdr(:) = 0.0
        vhlmlr(:) = qhlmlr(:) ! not actually volume, but treated as q in rate equation
!        vhlmlr(:) = rho0(:)*qhlmlr(:)/xdn(:,lhl) 
!        vhlsoak(:) = 0.0

      endif  ! mixedphase



!
!  Graupel reflectivity
!
      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'graupel reflectivity'

      do mgs = 1,ngscnt
      
!      zhmlr(mgs) = 0.0
!      zhshr(mgs) = 0.0
!      zhmlrr(mgs) = 0.0
!      zhshrr(mgs) = 0.0
      zhdsv(mgs) = 0.0
!      IF ( lf < 1 ) THEN
      if ( ffrzh > 0.0 ) then
      ziacr(mgs) = 0.0
      ziacrf(mgs) = 0.0
      endif
!      ENDIF
      zhcns(mgs) = 0.0
      zhcni(mgs) = 0.0
      zhacs(mgs) = 0.0
      zhaci(mgs) = 0.0
      
      enddo

      if ( lzh .gt. 1 ) then ! 
      do mgs = 1,ngscnt
      
      
      if ( qx(mgs,lh) .gt. qxmin(lh) .and. cx(mgs,lh) .gt. 0.0 ) then
          tmp = qx(mgs,lh)/cx(mgs,lh)
          alp = Max( alphamin, alpha(mgs,lh) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
!          g1r = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)

           zhaci(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhaci(mgs) )
           zhacs(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhacs(mgs) )
        
        if ( .not. mixedphase  .and. ibinhmlr < 1 ) then
        zhmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlr(mgs)  )
        endif
        
        zhshr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs)  )

!        IF ( lzr > 0 .and. qhshr(mgs) /= 0.0 .and. chshrr(mgs) /= 0.0 .and. ibinhmlr < 1 ) THEN
        if ( lzr > 0 .and. qhshr(mgs) /= 0.0 .and. chshrr(mgs) /= 0.0 ) then
!         IF ( temg(mgs) > tfr + 2.0 ) THEN
!           zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshrr(mgs)  )
!           IF ( zhshrr(mgs) > 0. ) THEN
!             zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs) )
!           ENDIF
!           z1 = g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) ! should this be g1shr?
!           zhshrr(mgs) = Max( z1, zhshrr(mgs))
!         ELSE
!          zhshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  )


         if ( temg(mgs) >= tfr ) then
 !           zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshrr(mgs)  )
 !           IF ( zhshrr(mgs) > 0.0 ) THEN
 !             zhshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs)  )
 !           ENDIF
           if ( (shedalp + alpha(mgs,lh))*xdia(mgs,lh,1) < sheddiam ) then ! if not shedding small drops, then use alpha of hail
             z1 = g1*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) 
           else
             z1 = g1shr*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) ! should this be g1shr?
           endif
           zhshrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  ) ! should this be g1shr?
!           zhshrr(mgs) = Max( z1, zhshrr(mgs))
         else
          zhshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhshr(mgs)**2/ chshrr(mgs)  )
         endif
         
         zhshrr(mgs) = Min( 0.0, zhshrr(mgs) )
        endif

        if ( zhshr(mgs) > 0.0 ) then
          write(0,*) 'Problem with zhshr! zhshr,qhshr,chshr = ',zhshr(mgs),qhshr(mgs),chshr(mgs)
          write(0,*) 'g1,tmp, qx,cx,zx = ',g1,tmp,qx(mgs,lh),cx(mgs,lh),zx(mgs,lh)
          write(0,*) ( 2.*tmp * qhshr(mgs) - tmp**2 * chshr(mgs)  ),  2.*tmp * qhshr(mgs), - tmp**2 * chshr(mgs)
          write(0,*) 'temcg = ',temcg(mgs),'chshr recalc = ',(cx(mgs,lh)/(qx(mgs,lh)+1.e-20))*qhshr(mgs)
          
          stop
        endif


!        zhshr(mgs) =  (xdn0(lr)/(xdn(mgs,lh)))**2*( zx(mgs,lh) * qhshr(mgs) )
        
        qtmp = qhdpv(mgs) + qhcev(mgs) + qhsbv(mgs)
        ctmp = chdpv(mgs) + chcev(mgs) + chsbv(mgs)

        zhdsv(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

          alp = Max( alphahacx, alpha(mgs,lh) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lh) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))

          if ( .true. ) then  ! {
          if ( qhacr(mgs) .gt. 0.0 ) then
!          zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacr(mgs) )

!          g1r = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
!          zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacr(mgs) )
          zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacr(mgs) )
!          zhacrf(mgs) = g1*zhacr


!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lh)+dtp*qhacr(mgs))**2)/(cx(mgs,lh))

          if ( z > zx(mgs,lh) ) then
!            zhacr(mgs) = (z - zx(mgs,lh))*dtpinv
          else
!            zhacr(mgs) = 0.0
          endif
          endif

!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) )
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) - tmp**2 * chacr(mgs) )

!          alp = Max( 1.0, alpha(mgs,lh)+1. )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/
!     :         ((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          if ( qhacw(mgs) .gt. 0.0 ) then
!          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )

!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lh)+dtp*(qhacw(mgs)-qhmul1(mgs)))**2)/(cx(mgs,lh))
          if ( z > zx(mgs,lh) ) then
!            zhacw(mgs) = (z - zx(mgs,lh))*dtpinv
          endif
          endif

          else ! } { ! this is not used because of the 'true' above

          if ( qhacw(mgs) .gt. 0.0 .or. qhacr(mgs) .gt. 0.0 ) then
          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lh)+dtp*(qhacr(mgs) + qhacw(mgs)-qhmul1(mgs)))**2)/(cx(mgs,lh))
!          zhacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lh)/cx(mgs,lh)) * qhacw(mgs) )
          if ( z > zx(mgs,lh) ) then
            zhacw(mgs) = (z - zx(mgs,lh))*dtpinv
          endif
          endif

          endif ! }

          if ( qhlcnh(mgs) .gt. 0.0 .and. ihlcnh < 2  ) then
           zhlcnh(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*( tmp ) * qhlcnh(mgs) - tmp**2 * chlcnh(mgs) )
          endif
      endif
! qsplinter(mgs)
      if ( ffrzh*qiacrf(mgs) .gt. 0.0 .and. cx(mgs,lr) .gt. 0.0 .and. qx(mgs,lr) .gt. qxmin(lr) ) then
            tmp = qx(mgs,lr)/cx(mgs,lr)
!            alp = 3.0
!            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            if ( imurain == 3 ) then
            ! note that 3.6476 = (6/pi)**2
            ziacr(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.))*  &
     &           ( 2.*tmp * qiacrf(mgs) - tmp**2 * ciacrf(mgs)  )
            else ! imurain == 1 
            ziacr(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2)*  &
     &           ( 2.*tmp * qiacrf(mgs) - tmp**2 * ciacrf(mgs)  )
            endif
            ziacr(mgs) = Min( ziacr(mgs), zxmxd(mgs,lr) )
!            ziacrf(mgs) = (xdn(mgs,lr)/xdn(mgs,lh))**2 * ziacr(mgs)
            ziacrf(mgs) = (xdn(mgs,lr)/xdnmx(lh))**2 * ziacr(mgs)
!            z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * (qiacrf(mgs) - qsplinter(mgs)) - tmp**2 * ciacrf(mgs)  )
!            ziacrf(mgs) = Min(  ziacrf(mgs), z )
      endif
      
      
      
      if ( ffrzh*qrfrzf(mgs) .gt. 0.0 .and. cx(mgs,lr) .gt. 0.0 ) then
            tmp = qx(mgs,lr)/cx(mgs,lr)
!            alp = 3.0
!            g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
            if ( imurain == 3 ) then
            zrfrz(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,lr)+2.)/(xdn0(lr)**2*(alpha(mgs,lr)+1.)) * &
     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
            zrfrzf(mgs) = (xdn(mgs,lr)/xdn(mgs,lh))**2 * zrfrz(mgs)
            elseif ( imurain == 1 .and. ibiggopt /= 2 ) then
!            zrfrz(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2) * &
!     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrz(mgs)  )
            zrfrz(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(xdn0(lr)**2) * &
     &         ( 2.*tmp * qrfrz(mgs) - tmp**2 * crfrz(mgs)  )
            zrfrzf(mgs) = 3.6476*rho0(mgs)**2*g1x(mgs,lr)/(rhofrz**2) * &
     &         ( 2.*tmp * qrfrzf(mgs) - tmp**2 * crfrzf(mgs)  )
            endif
            zrfrz(mgs) = Min( zrfrz(mgs), Max(0.4,qrfrz(mgs)/qx(mgs,lr))*zx(mgs,lr)*dtpinv )
!            zrfrzf(mgs) = (xdn(mgs,lr)/xdn(mgs,lh))**2 * zrfrz(mgs)
!            zrfrzf(mgs) = (xdn(mgs,lr)/xdnmx(lh))**2 * zrfrz(mgs)
!            z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * (qrfrzf(mgs)-qsplinter2(mgs)) - tmp**2 * crfrzf(mgs)  )
!             zrfrzf(mgs) = Min(  zrfrzf(mgs), z )
      ! change this to be alpha=0?
      endif
      
      if ( lhl > 1 .and. qhcnhl(mgs) .gt. 0.0 ) then
        tmp = qx(mgs,lhl)/cx(mgs,lhl)
        zhcnhl(mgs) = g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*( tmp ) * qhcnhl(mgs) - tmp**2 * chcnhl(mgs) )
        
      endif
      
      if ( qhcns(mgs) > 0.0 .and. chcns(mgs) > 0.0 .and. cx(mgs,ls) > cxmin .and. vhcns(mgs) > 0 ) then
        tmp = qx(mgs,ls)/cx(mgs,ls)
        r = rho0(mgs)*qhcns(mgs)/vhcns(mgs) ! density of new graupel particles
        if ( imusnow == 3 ) then
        zhcns(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,ls)+2.)/(r**2*(alpha(mgs,ls)+1.)) * &
     &         ( 2.*tmp * qhcns(mgs) - tmp**2 * chcns(mgs)  )
        else
         write(0,*) 'value of imusnow not valid. Must be 3 (fix me for =1). imusnow = ',imusnow
        stop
        endif
      endif

      if ( qhcni(mgs) > 0.0 .and. chcni(mgs) > 0.0 .and. cx(mgs,li) > cxmin .and. vhcni(mgs) > 0 ) then
        tmp = qx(mgs,li)/cx(mgs,li)
        r = rho0(mgs)*qhcni(mgs)/vhcni(mgs) ! density of new graupel particles
        zhcni(mgs) = 3.6476*rho0(mgs)**2*(alpha(mgs,li)+2.)/(r**2*(alpha(mgs,li)+1.)) * &
     &         ( 2.*tmp * qhcni(mgs) - tmp**2 * chcni(mgs)  )
      endif
 

      pzhwi(mgs) =   &
     &  +ifrzg*ffrzh*(zrfrzf(mgs)   &
     & +il5(mgs)*ifiacrg*(ziacrf(mgs) ) )   &
!     : + zhcnsh(mgs) + zhcnih(mgs)   &
     & + zhacw(mgs)   &
     & + zhacr(mgs)   &
     & + zhcnhl(mgs)  &
     & + zhacs(mgs)   &
     & + zhaci(mgs)   &
     &  + f2h*zhcni(mgs) + f2h*zhcns(mgs) &
     & + Max( 0.0, zhdsv(mgs) )

      pzhwd(mgs) = 0.0   &
     & + (1-il5(mgs))*zhmlr(mgs)   &
     & + zhshr(mgs)   &
!     >  + il5(mgs)*chsbv(mgs)   &
     &  + Min( 0.0, zhdsv(mgs) )   &
     &  - il5(mgs)*zhlcnh(mgs)


           if ( igs(mgs) == 44 .and. kgs(mgs) == 23 .or. dtp*( pqhwi(mgs) + pqhwd(mgs) ) > qxmin(lh) ) then
!             write(0,*)  'i,k,time = ',igs(mgs),kgs(mgs),time_real
!             write(0,*) 'pzhwi,d = ',pzhwi(mgs),pzhwd(mgs),dtp*( pzhwi(mgs) + pzhwd(mgs) ),zx(mgs,lh)
!             write(0,*) 'pqhwi,d = ',pqhwi(mgs),pqhwd(mgs),dtp*( pqhwi(mgs) + pqhwd(mgs) ),qx(mgs,lh)
!             write(0,*) 'pchwi,d = ',pchwi(mgs),pchwd(mgs),dtp*( pchwi(mgs) + pchwd(mgs) ),cx(mgs,lh)
           endif


!        IF ( zhcnhl(mgs) < 0.0 ) THEN
!          write(0,*) 'Problem with zhcnhl! zhcnhl,qhcnhl,chcnhl = ',zhcnhl(mgs),qhcnhl(mgs),chcnhl(mgs)
!          write(0,*) 'g1,tmp = ',g1x(mgs,lhl),tmp
!          write(0,*) ( 2.*( tmp ) * qhcnhl(mgs) - tmp**2 * chcnhl(mgs) )
!          
!!          STOP
!        ENDIF
      end do

      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'end graupel reflectivity'
      
      endif

!
!  Hail reflectivity
!

      do mgs = 1,ngscnt
      
      zhldsv(mgs) = 0.0
      zhlacr(mgs) = 0.0
      zhlacw(mgs) = 0.0
      
      enddo

      if ( lzhl .gt. 1 .or. ( lzr > 1 .and. lnhl > 1 ) ) then ! also run for 2-moment hail for 3-moment rain sources

      if (ndebug .gt. 0 .and. my_rank>=0 ) write(0,*) my_rank, 'hail reflectivity'

      do mgs = 1,ngscnt
      
      if ( qx(mgs,lhl) .gt. qxmin(lhl) .and. cx(mgs,lhl) .gt. 0.0 ) then
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
          alp = Max( alphamin, alpha(mgs,lhl) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
        
        if ( .not. mixedphase .and. qhlmlr(mgs) /= 0.0 .and. chlmlr(mgs) /= 0.0 .and. ibinhlmlr < 1 ) then
         zhlmlr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlr(mgs)  )
        endif
        
        zhlshr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshr(mgs)  )
        if ( lzr > 1 .and. qhlshr(mgs) /= 0.0 .and. chlshrr(mgs) /= 0.0 ) then
         if ( temg(mgs) >= tfr ) then
 !           zhlshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshrr(mgs)  )
 !           IF ( zhlshrr(mgs) > 0.0 ) THEN
 !             zhlshrr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn0(lr)))**2*( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshr(mgs)  )
 !           ENDIF
           if ( (shedalp + alpha(mgs,lhl))*xdia(mgs,lhl,1) < sheddiam ) then ! if not shedding small drops, then use alpha of hail
             z1 = g1*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  ) 
           else
             z1 = g1shr*(6.0*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  ) ! should this be g1shr?
           endif
           zhlshrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  ) ! should this be g1shr?
!           zhlshrr(mgs) = Max( z1, zhlshrr(mgs))
         else
          zhlshrr(mgs) =  g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlshr(mgs)**2/ chlshrr(mgs)  )
         endif

          zhlshrr(mgs) = Min( 0.0, zhlshrr(mgs) )
        endif

        if ( zhlshr(mgs) > 0.0 ) then
          write(0,*) 'Problem with zhlshr! zhlshr,qhlshr,chlshr = ',zhlshr(mgs),qhlshr(mgs),chlshr(mgs)
          write(0,*) 'g1,tmp, qx,cx,zx = ',g1,tmp,qx(mgs,lhl),cx(mgs,lhl),zx(mgs,lhl)
          write(0,*) ( 2.*tmp * qhlshr(mgs) - tmp**2 * chlshr(mgs)  ),  2.*tmp * qhlshr(mgs), - tmp**2 * chlshr(mgs)
          write(0,*) 'temcg = ',temcg(mgs),'chlshr recalc = ',(cx(mgs,lhl)/(qx(mgs,lhl)+1.e-20))*qhlshr(mgs)
          
          stop
        endif
!        zhlshr(mgs) = Min( 0.0, zhlshr(mgs) )

!        zhlshr(mgs) =  (xdn0(lr)/(xdn(mgs,lhl)))**2*( zx(mgs,lhl) * qhlshr(mgs) )
        
        qtmp = qhldpv(mgs) + qhlcev(mgs)
        ctmp = chldpv(mgs) + chlcev(mgs)
        
        zhldsv(mgs) = g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

          alp = Max( alphahacx, alpha(mgs,lhl) )
!          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))

          if ( .true. ) then ! {
          if ( qhlacr(mgs) .gt. 0.0 ) then
!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lhl)+dtp*qhlacr(mgs))**2)/(cx(mgs,lhl))
          zhlacr(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*( tmp ) * qhlacr(mgs) )
!          zhlacr(mgs) = Min( zxmxd(mgs,lr), zhlacr(mgs) )
          
!          IF ( z > zx(mgs,lhl) ) THEN
!            zhlacr(mgs) = (z - zx(mgs,lhl))*dtpinv
!          ELSE
!            zhlacr(mgs) = 0.0
!          ENDIF
          endif

!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) )
!        zhacr(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( tmp ) * qhacr(mgs) - tmp**2 * chacr(mgs) )

          if ( qhlacw(mgs) .gt. 0.0 ) then
          alp = Max( 3.0, alpha(mgs,lhl)+1. )
          g1 = (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          
!          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lhl)+dtp*(qhlacw(mgs)-qhlmul1(mgs)))**2)/(cx(mgs,lhl))
!          zhlacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lhl)/cx(mgs,lhl)) * qhlacw(mgs) )
          zhlacw(mgs) =  g1*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlacw(mgs) )

!          IF ( z > zx(mgs,lhl) ) THEN
!            zhlacw(mgs) = (z - zx(mgs,lhl))*dtpinv
!          ENDIF
          g1 = g1x(mgs,lhl) ! (6.0 + alp)*(5.0 + alp)*(4.0 + alp)/((3.0 + alp)*(2.0 + alp)*(1.0 + alp))
          endif
          
          else ! }  .false. {

          if ( qhlacw(mgs) .gt. 0.0 .or. qhlacr(mgs) .gt. 0.0 ) then
          z = g1*(6.*rho0(mgs)/(pi*1000.))**2*( (qx(mgs,lhl)+dtp*(qhlacr(mgs) + qhlacw(mgs)-qhlmul1(mgs)))**2)/(cx(mgs,lhl))
!          zhlacw(mgs) =  g1*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*( qx(mgs,lhl)/cx(mgs,lhl)) * qhlacw(mgs) )
          if ( z > zx(mgs,lhl) ) then
            zhlacw(mgs) = (z - zx(mgs,lhl))*dtpinv
          endif
          endif
          
          endif ! }
        
      endif
! qsplinter(mgs)
      
      if ( lzhl > 1 ) then
      pzhli(mgs) = ffrzh*(((1.0-ifrzg)*zrfrzf(mgs)   &
     & +il5(mgs)*(1.0-ifiacrg)*ziacrf(mgs) )) &
     &  + il5(mgs)*zhlcnh(mgs)   &
     & + zhlacw(mgs)   &
     & + zhlacr(mgs)   &
!     : + zhlacs(mgs)   &
     & + Max( 0.0, zhldsv(mgs) )

      pzhld(mgs) = 0.0   &
     & + (1-il5(mgs))*zhlmlr(mgs)   &
     & + zhlshr(mgs)   &
     & - zhcnhl(mgs)   &
     &  + Min( 0.0, zhldsv(mgs) )
      

       if ( .not. ( -1.0 < pzhli(mgs) .and. pzhli(mgs) < 1.e20 ) ) then
         write(iunit,*) 'Problem with pzhli!'
         write(iunit,*) 'zhlcnh,zhlacw,zhlacr,zhldsv = ',zhlcnh(mgs),zhlacw(mgs),zhlacr(mgs),zhldsv(mgs)
       endif

       if ( .not. ( -1.0e20 < pzhld(mgs) .and. pzhld(mgs) < 1. ) ) then
         write(iunit,*) 'Problem with pzhld!'
         write(iunit,*) 'zhlmlr,zhlshr,zhldsv = ',zhlmlr(mgs),zhlshr(mgs),zhldsv(mgs)
       endif
       
      endif ! lzhl > 1
      
      end do
      
      endif

!
!  rain reflectivity
!
      if (ndebug .gt. 0 ) write(0,*) 'WARMZIEG: dbg = 11'

      if ( lzr .gt. 1 ) then ! 
       
        do mgs = 1,ngscnt
        
        zracw(mgs) = 0.0
        zracr(mgs) = 0.0
        zrcev(mgs) = 0.0
        zrach(mgs) = 0.0
        zrachl(mgs) = 0.0
        zsshr(mgs) = 0.0
        zsshrr(mgs) = 0.0
!        zsmlr(mgs) = 0.0
        zsmlrr(mgs) = 0.0

        if ( qx(mgs,ls) .gt. qxmin(ls) .and. ( csmlr(mgs) /= 0.0 .or. csshr(mgs) /= 0.0 .or. &
              csmlrr(mgs) /= 0.0 .or. csshrr(mgs) /= 0.0) ) then !{
         tmp = qx(mgs,ls)/cx(mgs,ls)
         g1 = 36.*(xnu(ls)+2.0)/((xnu(ls)+1.0)*pi**2)
        if ( .not. mixedphase ) then
!          zsmlr(mgs) =  (xdn(mgs,ls)/xdn(mgs,lr))**2*g1*(rho0(mgs)/(xdn(mgs,ls)))**2* &
!     &                 ( 2.*tmp * qsmlr(mgs) - tmp**2 * csmlr(mgs)  )

          if ( csmlrr(mgs) /= 0.0 ) then
            z1 = g1smlr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qsmlr(mgs)**2/ csmlrr(mgs)  )
            zsmlrr(mgs) = z1
          endif
        endif
        
!        zsshr(mgs) =  (xdn(mgs,ls)/xdn(mgs,lr))**2*g1*(rho0(mgs)/(xdn(mgs,ls)))**2*  &
!     &                 ( 2.*tmp * qsshr(mgs) - tmp**2 * csshr(mgs)  )

         if ( csshrr(mgs) /= 0.0 ) then
          z1 = g1smlr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qsshr(mgs)**2/ csshrr(mgs)  )
          zsshrr(mgs) = z1
         endif
        
        endif !}
        
        if ( .not. mixedphase ) then !{
          if ( zhmlr(mgs) < 0.0 .and. chmlrr(mgs) /= 0.0 .and. ibinhmlr == 0 ) then !{
          tmp = qx(mgs,lh)/cx(mgs,lh)
!          zhmlrr(mgs) =  Min(0.0, (xdn(mgs,lh)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lh)))**2*( 2.*tmp * qhmlr(mgs) - tmp**2 * chmlrr(mgs)  ) )
            
!            IF ( zhmlrr(mgs) >= 0. ) THEN
!              zhmlrr(mgs) =  (xdn(mgs,lh)/xdn(mgs,lr))**2 * zhmlr(mgs)
!            ENDIF
           if ( (shedalp + alpha(mgs,lh))*xdia(mgs,lh,1) < sheddiam ) then ! if not shedding small drops, then use alpha of graupel
             z1 = g1x(mgs,lh)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhmlr(mgs)**2/ chmlrr(mgs)  ) 
           else ! assume drops are shed off, so use either alpha for shedding or graupel alpha, whichever gives the lower g-factor (i.e., larger alpha)
             z1 = Min(g1x(mgs,lh),g1shr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhmlr(mgs)**2/ chmlrr(mgs)  )
           endif
           zhmlrr(mgs) = z1
!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhmlr(mgs)**2/ chmlrr(mgs)  ) 
!           zhmlrr(mgs) = Max( z1, zhmlrr(mgs))
          endif !}


!          zhshrr(mgs) =  (xdn(mgs,lh)/xdn(mgs,lr))**2 * zhshr(mgs)

         if ( lhl > 1 .and. qhlmlr(mgs) /= 0 .and. ibinhlmlr == 0) then
          tmp = qx(mgs,lhl)/cx(mgs,lhl)
!          zhlmlrr(mgs) =  Min(0.0, (xdn(mgs,lhl)/xdn(mgs,lr))**2 * &
!     &       g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lhl)))**2*( 2.*tmp * qhlmlr(mgs) - tmp**2 * chlmlrr(mgs)  ) )

!          IF ( zhlmlrr(mgs) >= 0. ) THEN ! should be negative, if not, then use alternate calculation
!           zhlmlrr(mgs) =  (xdn(mgs,lhl)/xdn(mgs,lr))**2 * zhlmlr(mgs)
!          ENDIF

           if ( (shedalp + alpha(mgs,lhl))*xdia(mgs,lhl,1) < sheddiam ) then ! if not shedding small drops, then use alpha of hail
             z1 = g1x(mgs,lhl)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  ) 
           else ! assume drops are shed off, so use either alpha for shedding or graupel alpha, whichever gives the lower g-factor (i.e., larger alpha)
             z1 = Min(g1x(mgs,lhl),g1shr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  )
!             z1 = g1shr*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  )
           endif
           zhlmlrr(mgs) = z1

!           z1 = g1mlr*(rho0(mgs)/(xdn(mgs,lr)))**2*( qhlmlr(mgs)**2/ chlmlrr(mgs)  )
!           zhlmlrr(mgs) = Max( z1, zhlmlrr(mgs))
!         zhlmlr(mgs) =
!          zhlshrr(mgs) =  (xdn(mgs,lhl)/xdn(mgs,lr))**2 * zhlshr(mgs)
         endif
         
         endif ! }

        if ( qx(mgs,lr) .gt. qxmin(lr) .and. cx(mgs,lr) .gt. 0.0 ) then

          tmp = qx(mgs,lr)/cx(mgs,lr)
          g1 = g1x(mgs,lr) ! 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)


        if ( qracw(mgs) > 0.0 .and. cx(mgs,lr) > 0.0 ) then
         zracw(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*1000.))**2*( 2.*tmp * qracw(mgs) )
        endif
        
        if ( cracr(mgs) > 0.0 .and. cx(mgs,lr) > 0.0  ) then
         zracr(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*1000.))**2*( tmp**2 * cracr(mgs) )
        endif

        qtmp = qrcev(mgs)
        ctmp = crcev(mgs)
        
!        IF ( .false. .or. iferwisventr == 2 ) THEN
!        zrcev(mgs) = Min(0.0, (12./(pii*xdn(mgs,lr)))*xdia(mgs,lr,1)**3*fvce(mgs)*rwcap(mgs)*rwventz(mgs) )
!        ELSE
        zrcev(mgs) = g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*( 2.*( tmp ) * qtmp - tmp**2 * ctmp )

        
        if (  iferwisventr == 2 ) then
          vent1 = Min(0.0, (12./(pii*xdn(mgs,lr)))*xdia(mgs,lr,1)**3*fvce(mgs)*rwcap(mgs)*rwventz(mgs))
          zrcev(mgs) = Max( zrcev(mgs), vent1 )
        endif
!        IF ( ny == 2 .and. igs(mgs) == 20 ) THEN
!          write(0,*) 'k,zrcevold,new,maxdep : ',kgs(mgs),zrcev(mgs),vent1,-zxmxd(mgs,lr),alpha(mgs,lr),cx(mgs,lr)
!        ENDIF


!        ENDIF
        zrcev(mgs) = Max( zrcev(mgs), -zxmxd(mgs,lr) )

        if ( qhacr(mgs) > 0.0 ) then 
          zrach(mgs) =  g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2* &
     &     ( 2.*( qx(mgs,lr)/cx(mgs,lr)) * qhacr(mgs) - tmp**2 * chacr(mgs) )
          zrach(mgs) = Min( zrach(mgs), zxmxd(mgs,lr) )
         
         endif

        if ( lhl > 1 .and. qhlacr(mgs) > 0.0 ) then 
          zrachl(mgs) = g1x(mgs,lr)*(6.*rho0(mgs)/(pi*xdn(mgs,lr)))**2*   &
     &     ( 2.*( qx(mgs,lr)/cx(mgs,lr)) * qhlacr(mgs) - tmp**2 * chlacr(mgs) )
          zrachl(mgs) = Min( zrachl(mgs), zxmxd(mgs,lr) )
         endif


        
        endif

         pzrwi(mgs) = zrcnw(mgs) + zracw(mgs) + zracr(mgs) &
     &    + Max( 0.,zrcev(mgs) )  &
     &  - (1-il5(mgs))*zsmlrr(mgs)   &
     &  - zsshrr(mgs)   &
     &  - (1-il5(mgs))*zhmlrr(mgs)   &
     &  - zhshrr(mgs)   &
     &  - (1-il5(mgs))*zhlmlrr(mgs)   &
     &  - zhlshrr(mgs)   


         pzrwd(mgs) = 0.0   &
     &   +  Min(0.,zrcev(mgs) )  &
     &    - zrach(mgs)  &
     &    - zrachl(mgs)  &
     &    - zrfrz(mgs)  &
     &    - il5(mgs)*(ziacr(mgs) ) 


         if ( zx(mgs,lr) + dtp*(pzrwi(mgs)+pzrwd(mgs))  <= 0.0  &
              .and. qx(mgs,lr) > qxmin(lr) ) then
           pzrwd(mgs) =  -zx(mgs,lr)*dtpinv - pzrwi(mgs)
         endif

        enddo

      endif



!
!  Snow volume
!
      if ( lvol(ls) .gt. 1 ) then
      do mgs = 1,ngscnt
!      pvswi(mgs) = rho0(mgs)*( pqswi(mgs) )/xdn0(ls)

      pvswi(mgs) = rho0(mgs)*(    &
!aps     >   il5*qsfzs(mgs)/xdn(mgs,ls)   &
!aps     >  -il5*qsfzs(mgs)/xdn(mgs,lr)   &
     &  +il5(mgs)*(qscni(mgs)+qsaci(mgs)+qsdpv(mgs)   &
     &   + qscnvi(mgs) + (1. - ifrzs)*qiacrs(mgs) &
     &   + (1. - ifrzs)*qrfrzs(mgs)  &
     &  )/xdn0(ls)   &
     &    + (qsacr(mgs))/rimdn(mgs,ls) ) + vsacw(mgs)
!     >   + (qsacw(mgs) + qsacr(mgs))/rimdn(mgs,ls) )
      pvswd(mgs) = rho0(mgs)*( pqswd(mgs) )/xdn0(ls)  &
!     >  -qhacs(mgs)
!     >  -qhcns(mgs)
!     >  +(1-il5(mgs))*qsmlr(mgs) + qsshr(mgs)
!     >  +il5(mgs)*(qssbv(mgs))
     &   -rho0(mgs)*qsmul(mgs)/xdn0(ls)
!aps     >   +rho0(mgs)*(1-il5(mgs))*(
!aps     >             qsmlr(mgs)/xdn(mgs,ls)
!aps     >    +(qscev-qsmlr(mgs))/xdn(mgs,lr) )
      end do

!aps      IF (mixedphase) THEN
!aps        pvswd(mgs) = pvswd(mgs)
!aps     >   + rho0(mgs)*qsshr(mgs)/xdn(mgs,lr)
!aps      ENDIF

      endif
!
!  Graupel volume
!
      if ( lvol(lh) .gt. 1 ) then
      do mgs = 1,ngscnt
!      pvhwi(mgs) = rho0(mgs)*( (pqhwi(mgs) )/xdn0(lh) )

!      pvhwi(mgs) = rho0(mgs)*( (pqhwi(mgs) - il5(mgs)*qrfrzf(mgs) )/xdn0(lh) !
!     :  +  il5(mgs)*qrfrzf(mgs)/rhofrz )

      pvhwi(mgs) = rho0(mgs)*(   &
     &  +il5(mgs)*( ifiacrg*ffrzh*qracif(mgs))/rhofrz   &
!erm     >  + il5(mgs)*qhfzh(mgs)/rhofrz !aps: or use xdnmx(lh)?   &
     &  + (  il5(mgs)*qhdpv(mgs)/qhdpvdn   &
     &     + (qhacs(mgs) + qhaci(mgs))/qhacidn ) )   &
     &  +   rho0(mgs)*Max(0.0, qhcev(mgs))/1000.   & ! only used in mixed phase: evaporation/condensation of liquid water coating
!     >     + qhacs(mgs) + qhaci(mgs) )/xdn0(ls) )   &
     &  + f2h*vhcns(mgs)   &
     &  + vhacr(mgs) + vhacw(mgs)  + vhfzh(mgs)   & ! qhacw(mgs)/rimdn(mgs,lh)
!     >  + vhfrh(mgs)   &
     &  + f2h*vhcni(mgs) + (ifiacrg*viacrf(mgs) + ifrzg*vrfrzf(mgs))*ffrzh
!     >  +qhacr(mgs)/raindn(mgs,lh) + qhacw(mgs)/rimdn(mgs,lh)
      
!      pvhwd(mgs) = rho0(mgs)*(pqhwd(mgs) )/xdn0(lh)

      pvhwd(mgs) = rho0(mgs)*(   &
!     >   qhshr(mgs)/xdn0(lr)   &
!     >  - il5(mgs)*qhfzh(mgs)/xdn(mgs,lr)   &
     &  +( (1-il5(mgs))*vhmlr(mgs)    &
!     >     +il5(mgs)*qhsbv(mgs)   &
     &     + qhsbv(mgs)   &
     &     + Min(0.0, qhcev(mgs))   &
     &     -qhmul1(mgs) )/xdn(mgs,lh) )   &
     &  - vhlcnh(mgs) + vhshdr(mgs) - vhsoak(mgs) - vscnh(mgs)

!      IF (mixedphase) THEN
!       pvhwd(mgs) = pvhwd(mgs) 
!     >  + rho0(mgs)*qhshr(mgs)/xdn(mgs,lh) !xdn(mgs,lr)
!      ENDIF

       if ( lzh > 1 .and. qx(mgs,lh) > qxmin(lh) ) then
!       Calculate change in reflectivity due to density changes

        xdn_new = rho0(mgs)*(qx(mgs,lh) + dtp*(pqhwi(mgs) + pqhwd(mgs) ))/   &
     &   (vx(mgs,lh) + dtp*(pvhwi(mgs) + pvhwd(mgs))  )

           if ( mixedphase ) then 
             if ( qxw(mgs,lh) .gt. 0.0 ) then
               dnmx = xdnmx(lr)
             else
               dnmx = xdnmx(lh)
             endif
           else
             dnmx = xdnmx(lh)
           endif

        xdn_new = Max( Min( xdn_new, dnmx ), xdnmn(lh) )
        
        drhodt = (xdn_new - xdn(mgs,lh))*dtpinv
        
        zhwdn(mgs) = -2.*g1x(mgs,lh)*(rho0(mgs)*qx(mgs,lh)*6.*pii )**2/(cx(mgs,lh)*xdn(mgs,lh)**3)*drhodt
        
        pzhwi(mgs) = pzhwi(mgs) + Max(0.0, zhwdn(mgs))
        pzhwd(mgs) = pzhwd(mgs) + Min(0.0, zhwdn(mgs))
        
       
       endif
      if ( .false. .and. ny .eq. 2 .and. kgs(mgs) .eq. 9 .and. igs(mgs) .eq. 19 ) then

      write(iunit,*)
      write(iunit,*)   'Graupel at ',igs(mgs),kgs(mgs)
!
      write(iunit,*)   il5(mgs)*qrfrzf(mgs), qrfrzf(mgs) - qrfrz(mgs)
      write(iunit,*)   il5(mgs)*qiacrf(mgs)
      write(iunit,*)   il5(mgs)*qracif(mgs)
      write(iunit,*)   'qhcns',qhcns(mgs)
      write(iunit,*)   'qhcni',qhcni(mgs)
      write(iunit,*)   il5(mgs)*(qhdpv(mgs))
      write(iunit,*)   'qhacr ',qhacr(mgs)
      write(iunit,*)   'qhacw', qhacw(mgs)
      write(iunit,*)   'qhacs', qhacs(mgs)
      write(iunit,*)   'qhaci', qhaci(mgs)
      write(iunit,*)   'pqhwi = ',pqhwi(mgs)
      write(iunit,*)
      write(iunit,*) 'qhcev',qhcev(mgs)
      write(iunit,*)
      write(iunit,*)   'qhshr',qhshr(mgs)
      write(iunit,*)  'qhmlr', (1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)   'qhsbv', qhsbv(mgs)
      write(iunit,*)   'qhlcnh',-qhlcnh(mgs)
      write(iunit,*)   'qhmul1',-qhmul1(mgs)
      write(iunit,*)   'pqhwd = ', pqhwd(mgs)
      write(iunit,*)
      write(iunit,*)  'Volume'
      write(iunit,*)
      write(iunit,*)  'pvhwi',pvhwi(mgs)
      write(iunit,*)   'vhcns', vhcns(mgs)
      write(iunit,*)  'vhacr,vhacw',vhacr(mgs), vhacw(mgs) ! qhacw(mgs)/rimdn(mgs,lh)
      write(iunit,*)  'vhcni',vhcni(mgs)
      write(iunit,*)
      write(iunit,*)  'pvhwd',pvhwd(mgs)
      write(iunit,*)  'vhlcnh,vhshdr,vhsoak ', vhlcnh(mgs),  vhshdr(mgs), vhsoak(mgs)
      write(iunit,*)  'vhmlr', vhmlr(mgs)
      write(iunit,*)
!      write(iunit,*)
!      write(iunit,*)
!      write(iunit,*)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pchwi(mgs),pchwd(mgs)
      write(iunit,*)  crfrzf(mgs)
      write(iunit,*)  chcns(mgs)
      write(iunit,*)  ciacrf(mgs)


      endif


      enddo

      endif
!
!
!

!
!  Hail volume
!
      if ( lhl .gt. 1 ) then
      if ( lvol(lhl) .gt. 1 ) then
      do mgs = 1,ngscnt

      pvhli(mgs) = rho0(mgs)*(   &
     &  + (  il5(mgs)*(((1.0-ifiacrg)*ffrzh*qracif(mgs))/rhofrz  + qhldpv(mgs) )   &
!     &  +    Max(0.0, qhlcev(mgs))   &
!     &     + qhlacs(mgs) + qhlaci(mgs) )/xdnmn(lhl) )   & ! xdn0(ls) )   &
!     &     + qhlacs(mgs) + qhlaci(mgs) )/xdnmn(lh) )   &  ! yes, this is 'lh' on purpose
     &     + qhlacs(mgs) + qhlaci(mgs) )/500. )   &  ! changed to 500 instead of min graupel density to keep hail density from dropping too much
     &  +   rho0(mgs)*Max(0.0, qhlcev(mgs))/1000.   &
     &  + vhlcnhl(mgs) + ((1.0-ifiacrg)*ffrzh*viacrf(mgs) + (1.0-ifrzg)*ffrzh*vrfrzf(mgs))  & 
     &  + vhlacr(mgs) + vhlacw(mgs) + vhlfzhl(mgs) ! qhlacw(mgs)/rimdn(mgs,lhl)
      
      pvhld(mgs) = rho0(mgs)*(   &
     &  +(  qhlsbv(mgs)   &
     &     + Min(0.0, qhlcev(mgs))   &
     &     -qhlmul1(mgs) )/xdn(mgs,lhl) ) &
!     &   + vhlmlr(mgs)                    &
     &   + rho0(mgs)*(1-il5(mgs))*vhlmlr(mgs)/xdn(mgs,lhl)  &
     &   + vhlshdr(mgs) - vhlsoak(mgs)

       if ( lzhl > 1 .and. qx(mgs,lhl) > qxmin(lhl) ) then
!       Calculate change in reflectivity due to density changes

        xdn_new = rho0(mgs)*(qx(mgs,lhl) + dtp*(pqhli(mgs) + pqhld(mgs) ))/   &
     &   (vx(mgs,lhl) + dtp*(pvhli(mgs) + pvhld(mgs))  )
        
           if ( mixedphase ) then 
             if ( qxw(mgs,lhl) .gt. 0.0 ) then
               dnmx = xdnmx(lr)
             else
               dnmx = xdnmx(lhl)
             endif
           else
             dnmx = xdnmx(lhl)
           endif
        xdn_new = Max( Min( xdn_new, dnmx ), xdnmn(lhl) )
        
        drhodt = (xdn_new - xdn(mgs,lhl))*dtpinv
        
        zhldn(mgs) = -2.*g1x(mgs,lhl)*(rho0(mgs)*qx(mgs,lhl)*6.*pii )**2/(cx(mgs,lhl)*xdn(mgs,lhl)**3)*drhodt
        
        pzhli(mgs) = pzhli(mgs) + Max(0.0, zhldn(mgs))
        pzhld(mgs) = pzhld(mgs) + Min(0.0, zhldn(mgs))
        
       
       endif

      enddo
      
      endif
      endif


      if ( ndebug .ge. 1 ) then
      do mgs = 1,ngscnt
!
      ptotal(mgs) = 0.
      ptotal(mgs) = ptotal(mgs)     &
     &  + pqwvi(mgs) + pqwvd(mgs)   &
     &  + pqcwi(mgs) + pqcwd(mgs)   &
     &  + pqcii(mgs) + pqcid(mgs)   &
     &  + pqrwi(mgs) + pqrwd(mgs)   &
     &  + pqswi(mgs) + pqswd(mgs)   &
     &  + pqhwi(mgs) + pqhwd(mgs)   &
     &  + pqhli(mgs) + pqhld(mgs)
!

      
      
      enddo
      
      do mgs = 1,ngscnt

      if ( ( (ndebug .ge. 0  ) .and. abs(ptotal(mgs)) .gt. eqtot )   &
!      if ( (  abs(ptotal(mgs)) .gt. eqtot )
!     :    .or. pqswi(mgs)*dtp .gt. 1.e-3
!     :    .or. pqhwi(mgs)*dtp .gt. 1.e-3
!     :     .or. dtp*(pqrwi(mgs)+pqrwd(mgs)) .gt. 10.0e-3
!     :     .or. dtp*(pccii(mgs)+pccid(mgs)) .gt. 1.e7
!     :     .or. dtp*(pcipi(mgs)+pcipd(mgs)) .gt. 1.e7    &
     &  .or.  .not. (ptotal(mgs) .lt. 1.0 .and.  ptotal(mgs) .gt. -1.0)   & ! this line is basically checking for NaNs
     &              ) then
      write(iunit,*) 'YIKES! ','ptotal1',mgs,igs(mgs),jgs,   &
     &       kgs(mgs),ptotal(mgs)

      write(iunit,*) 't7: ', t7(igs(mgs),jgs,kgs(mgs))
      write(iunit,*)  'cci,ccw,crw,rdia: ',cx(mgs,li),cx(mgs,lc),cx(mgs,lr),0.5*xdia(mgs,lr,1)
      write(iunit,*)  'qc,qi,qr : ',qx(mgs,lc),qx(mgs,li),qx(mgs,lr)
      write(iunit,*)  'rmas, qrcalc : ',xmas(mgs,lr),xmas(mgs,lr)*cx(mgs,lr)/rho0(mgs)
      write(iunit,*)  'vti,vtc,eiw,vtr: ',vtxbar(mgs,li,1),vtxbar(mgs,lc,1),eiw(mgs),vtxbar(mgs,lr,1)
      write(iunit,*)  'cidia,cwdia,qcmxd: ', xdia(mgs,li,1),xdia(mgs,lc,1),qcmxd(mgs)
      write(iunit,*)  'snow: ',qx(mgs,ls),cx(mgs,ls),swvent(mgs),vtxbar(mgs,ls,1),xdia(mgs,ls,1)
      write(iunit,*)  'graupel: ',qx(mgs,lh),cx(mgs,lh),hwvent(mgs),vtxbar(mgs,lh,1),xdia(mgs,lh,1)
      if ( lhl .gt. 1 ) write(iunit,*)  'hail: ',qx(mgs,lhl),cx(mgs,lhl),hlvent(mgs),vtxbar(mgs,lhl,1),xdia(mgs,lhl,1)


      write(iunit,*)  'li: ',xdia(mgs,li,1),xdia(mgs,li,2),xmas(mgs,li),qx(mgs,li),   &
     &         vtxbar(mgs,li,1)


      write(iunit,*)  'rain cx,xv : ',cx(mgs,lr),xv(mgs,lr)
      write(iunit,*)  'temcg = ', temcg(mgs)

      write(iunit,*) 'v ', pqwvi(mgs) ,pqwvd(mgs)
      write(iunit,*) 'c ', pqcwi(mgs) ,pqcwd(mgs)
      write(iunit,*) 'ci', pqcii(mgs) ,pqcid(mgs)
      write(iunit,*) 'r ', pqrwi(mgs) ,pqrwd(mgs)
      write(iunit,*) 's ', pqswi(mgs) ,pqswd(mgs)
      write(iunit,*) 'h ', pqhwi(mgs) ,pqhwd(mgs)
      write(iunit,*) 'hl', pqhli(mgs) ,pqhld(mgs)
       tmp =  pqwvi(mgs) + pqwvd(mgs)   &
     &  + pqcwi(mgs) + pqcwd(mgs)   &
     &  + pqcii(mgs) + pqcid(mgs)   &
     &  + pqrwi(mgs) + pqrwd(mgs)   &
     &  + pqswi(mgs) + pqswd(mgs)   &
     &  + pqhwi(mgs) + pqhwd(mgs)   &
     &  + pqhli(mgs) + pqhld(mgs)

      write(iunit,*) 'total = ',tmp
      write(iunit,*) 'end OF OUTPUT OF SOURCE AND SINK'

!
!  print production terms
!
      write(iunit,*)
      write(iunit,*)   'Vapor'
!
      write(iunit,*)   -Min(0.0,qrcev(mgs))
      write(iunit,*)   -il5(mgs)*qhsbv(mgs)
      write(iunit,*)   -il5(mgs)*qhlsbv(mgs)
      write(iunit,*)   -il5(mgs)*qssbv(mgs)
      write(iunit,*)   -il5(mgs)*qisbv(mgs)
      write(iunit,*)    'pqwvi= ', pqwvi(mgs)
      write(iunit,*)   -Max(0.0,qrcev(mgs))
      write(iunit,*)   -Max(0.0,qhcev(mgs))
      write(iunit,*)   -Max(0.0,qhlcev(mgs))
      write(iunit,*)   -Max(0.0,qscev(mgs))
      write(iunit,*)   -il5(mgs)*qiint(mgs)
      write(iunit,*)   -il5(mgs)*qhdpv(mgs)
      write(iunit,*)   -il5(mgs)*qhldpv(mgs)
      write(iunit,*)   -il5(mgs)*qsdpv(mgs)
      write(iunit,*)   -il5(mgs)*qidpv(mgs)
      write(iunit,*)    'pqwvd = ', pqwvd(mgs)
!
      write(iunit,*)
      write(iunit,*)   'Cloud ice'
!
      write(iunit,*)   il5(mgs)*qicicnt(mgs)
      write(iunit,*)   il5(mgs)*qidpv(mgs)
      write(iunit,*)   il5(mgs)*qiacw(mgs)
      write(iunit,*)   il5(mgs)*qwfrzc(mgs)
      write(iunit,*)   il5(mgs)*qwctfzc(mgs)
      write(iunit,*)   il5(mgs)*qicichr(mgs)
      write(iunit,*)   qhmul1(mgs)
      write(iunit,*)   qhlmul1(mgs)
      write(iunit,*)   'pqcii = ', pqcii(mgs)
      write(iunit,*)   -il5(mgs)*qscni(mgs)
      write(iunit,*)   -il5(mgs)*qscnvi(mgs)
      write(iunit,*)   -il5(mgs)*qraci(mgs)
      write(iunit,*)   -il5(mgs)*qsaci(mgs)
      write(iunit,*)   -il5(mgs)*qhaci(mgs)
      write(iunit,*)   -il5(mgs)*qhlaci(mgs)
      write(iunit,*)   il5(mgs)*qisbv(mgs)
      write(iunit,*)   (1.-il5(mgs))*qimlr(mgs)
      write(iunit,*)   -il5(mgs)*qhcni(mgs)
      write(iunit,*)   'pqcid = ', pqcid(mgs)
      write(iunit,*)   ' Conc:'
      write(iunit,*)   pccii(mgs),pccid(mgs)
      write(iunit,*)   il5(mgs),cicint(mgs)
      write(iunit,*)   cwfrzc(mgs),cwctfzc(mgs)
      write(iunit,*)   cicichr(mgs)
      write(iunit,*)   chmul1(mgs)
      write(iunit,*)   chlmul1(mgs)
      write(iunit,*)   csmul(mgs)
!
!
!
!
      write(iunit,*)
      write(iunit,*)   'Cloud water'
!
      write(iunit,*)   'pqcwi =', pqcwi(mgs)
      write(iunit,*)   -il5(mgs)*qiacw(mgs)
      write(iunit,*)   -il5(mgs)*qwfrzc(mgs)
      write(iunit,*)   -il5(mgs)*qwctfzc(mgs)
!      write(iunit,*)   -il5(mgs)*qwfrzp(mgs)
!      write(iunit,*)   -il5(mgs)*qwctfzp(mgs)
      write(iunit,*)   -il5(mgs)*qiihr(mgs)
      write(iunit,*)   -il5(mgs)*qicichr(mgs)
      write(iunit,*)   -il5(mgs)*qipiphr(mgs)
      write(iunit,*)   -qracw(mgs)
      write(iunit,*)   -qsacw(mgs)
      write(iunit,*)   -qrcnw(mgs)
      write(iunit,*)   -qhacw(mgs)
      write(iunit,*)   -qhlacw(mgs)
      write(iunit,*)   'pqcwd = ', pqcwd(mgs)


      write(iunit,*)
      write(iunit,*)  'Concentration:'
      write(iunit,*)   -cautn(mgs)
      write(iunit,*)   -cracw(mgs)
      write(iunit,*)   -csacw(mgs)
      write(iunit,*)   -chacw(mgs)
      write(iunit,*)  -ciacw(mgs)
      write(iunit,*)  -cwfrzp(mgs)
      write(iunit,*)  -cwctfzp(mgs)
      write(iunit,*)  -cwfrzc(mgs)
      write(iunit,*)  -cwctfzc(mgs)
      write(iunit,*)   pccwd(mgs)
!
      write(iunit,*)
      write(iunit,*)      'Rain '
!
      write(iunit,*)      qracw(mgs)
      write(iunit,*)      qrcnw(mgs)
      write(iunit,*)      Max(0.0, qrcev(mgs))
      write(iunit,*)       -(1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qhlmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qsmlr(mgs)
      write(iunit,*)       -(1-il5(mgs))*qimlr(mgs)
      write(iunit,*)       -qrshr(mgs)
      write(iunit,*)       'pqrwi = ', pqrwi(mgs)    
      write(iunit,*)        -qsshr(mgs)     
      write(iunit,*)        -qhshr(mgs)     
      write(iunit,*)        -qhlshr(mgs)
      write(iunit,*)        -il5(mgs)*qiacr(mgs),qiacr(mgs), qiacrf(mgs)
      write(iunit,*)        -il5(mgs)*qrfrz(mgs)
      write(iunit,*)        -qsacr(mgs)
      write(iunit,*)        -qhacr(mgs)
      write(iunit,*)        -qhlacr(mgs)
      write(iunit,*)        qrcev(mgs)
      write(iunit,*)       'pqrwd = ', pqrwd(mgs) 
      write(iunit,*)        'qrzfac = ', qrzfac(mgs)
!
      
      write(iunit,*)
      write(iunit,*)  'Rain concentration'
      write(iunit,*)  pcrwi(mgs) 
      write(iunit,*)    crcnw(mgs)
      write(iunit,*)    1-il5(mgs)
      write(iunit,*)   -chmlr(mgs),-csmlr(mgs)
      write(iunit,*)     -crshr(mgs)
      write(iunit,*)  pcrwd(mgs) 
      write(iunit,*)    il5(mgs)
      write(iunit,*)   -ciacr(mgs),-crfrz(mgs) 
      write(iunit,*)   -csacr(mgs),-chacr(mgs)
      write(iunit,*)   +crcev(mgs)
      write(iunit,*)   cracr(mgs)
!      write(iunit,*)   -il5(mgs)*ciracr(mgs)


      write(iunit,*)
      write(iunit,*)   'Snow'
!
      write(iunit,*)        il5(mgs)*qscni(mgs), qscnvi(mgs)
      write(iunit,*)        il5(mgs)*qsaci(mgs)
      write(iunit,*)        il5(mgs)*qrfrzs(mgs)
      write(iunit,*)        il5(mgs)*qiacrs(mgs),il3(mgs)*(qiacrf(mgs)+qracif(mgs)),il3(mgs),qiacrf(mgs),qracif(mgs)
      write(iunit,*)        il5(mgs)*qsdpv(mgs), qscev(mgs)
      write(iunit,*)        qsacw(mgs)
      write(iunit,*)        qsacr(mgs), qscnh(mgs)
       write(iunit,*)        'pqswi = ',pqswi(mgs)
      write(iunit,*)        -qhcns(mgs)
      write(iunit,*)        -qracs(mgs)
      write(iunit,*)        -qhacs(mgs)
      write(iunit,*)        -qhlacs(mgs)
      write(iunit,*)       (1-il5(mgs))*qsmlr(mgs)
      write(iunit,*)       qsshr(mgs)
!      write(iunit,*)       qsshrp(mgs)
      write(iunit,*)       il5(mgs)*(qssbv(mgs))
      write(iunit,*)       'pqswd = ', pqswd(mgs)
      write(iunit,*)   -qracs(mgs)*(1-il2(mgs)) , qhacs(mgs) , qhlacs(mgs)   
      write(iunit,*)   -qhcns(mgs)   
      write(iunit,*)   +(1-il5(mgs))*qsmlr(mgs) , qsshr(mgs)     
      write(iunit,*)    qssbv(mgs)
      write(iunit,*)   Min(0.0, qscev(mgs))  
      write(iunit,*)   -qsmul(mgs)
!
!
      write(iunit,*)
      write(iunit,*)   'Graupel'
!
      write(iunit,*)   il5(mgs)*qrfrzf(mgs), qrfrzf(mgs) - qrfrz(mgs)
      write(iunit,*)   il5(mgs)*qiacrf(mgs)
      write(iunit,*)   il5(mgs)*qracif(mgs)
      write(iunit,*)   qhcns(mgs)
      write(iunit,*)   qhcni(mgs)
      write(iunit,*)   il5(mgs)*(qhdpv(mgs))
      write(iunit,*)   qhacr(mgs)
      write(iunit,*)   qhacw(mgs)
      write(iunit,*)   qhacs(mgs)
      write(iunit,*)   qhaci(mgs)
      write(iunit,*)   'pqhwi = ',pqhwi(mgs)
      write(iunit,*)
      write(iunit,*)   qhshr(mgs)
      write(iunit,*)   (1-il5(mgs))*qhmlr(mgs)
      write(iunit,*)   il5(mgs),qhsbv(mgs)
      write(iunit,*)   -qhlcnh(mgs)
      write(iunit,*)   -qhmul1(mgs)
      write(iunit,*)   'pqhwd = ', pqhwd(mgs)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pchwi(mgs),pchwd(mgs)
      write(iunit,*)  crfrzf(mgs)
      write(iunit,*)  chcns(mgs)
      write(iunit,*)  ciacrf(mgs)

!
      write(iunit,*)
      write(iunit,*)   'Hail'
!
      write(iunit,*)   qhlcnh(mgs)
      write(iunit,*)   il5(mgs)*(qhldpv(mgs))
      write(iunit,*)   qhlacr(mgs)
      write(iunit,*)   qhlacw(mgs)
      write(iunit,*)   qhlacs(mgs)
      write(iunit,*)   qhlaci(mgs)
      write(iunit,*)   pqhli(mgs)
      write(iunit,*)
      write(iunit,*)   qhlshr(mgs)
      write(iunit,*)   (1-il5(mgs))*qhlmlr(mgs)
      write(iunit,*)   il5(mgs)*qhlsbv(mgs)
      write(iunit,*)   pqhld(mgs)
      write(iunit,*)  'Concentration'
      write(iunit,*)   pchli(mgs),pchld(mgs)
      write(iunit,*)  chlcnh(mgs)
!
!  Balance and checks for continuity.....within machine precision...
!
!
      write(iunit,*) 'end OF OUTPUT OF SOURCE AND SINK'
      write(iunit,*) 'PTOTAL',ptotal(mgs)
!
      end if ! ptotal out of bounds or NaN
!
      end do
!

      end if ! ( nstep/12*12 .eq. nstep )

!
!  latent heating from phase changes (except qcw, qci cond, and evap)
!
      do mgs = 1,ngscnt
      if ( warmonly < 0.5 ) then
      pfrz(mgs) =    &
     &  (1-il5(mgs))*   &
     &  (qhmlr(mgs)+    &
     &   qsmlr(mgs)+qhlmlr(mgs))   & !+qhmlh(mgs))   &
     &  +il5(mgs)*(1-imixedphase)*(   &
     &   qsacw(mgs)+qhacw(mgs) + qhlacw(mgs)   &
     &  +qsacr(mgs)+qhacr(mgs) + qhlacr(mgs)   &
     &  +qsshr(mgs)   &
     &  +qhshr(mgs)   &
     &  +qhlshr(mgs)  &
     &  +qrfrz(mgs)+qiacr(mgs)  &
     &  )  &
     &  +il5(mgs)*(qwfrz(mgs)    &
     &  +qwctfz(mgs)+qiihr(mgs)   &
     &  +qiacw(mgs))
      pmlt(mgs) =    &
     &  (1-il5(mgs))*   &
     &  (qhmlr(mgs)+qsmlr(mgs)+  &
     &   qhlmlr(mgs))    !+qhmlh(mgs))   
      ! NOTE: psub is sum of sublimation and deposition
      psub(mgs) =    &
     &   il5(mgs)*(   &
     &  + qsdpv(mgs) + qhdpv(mgs)   &
     &  + qhldpv(mgs)    &
     &  + qidpv(mgs) + qisbv(mgs) )   &
     &   + qssbv(mgs)  + qhsbv(mgs) &
     &  + qhlsbv(mgs)   &
     &  +il5(mgs)*(qiint(mgs))
      pvap(mgs) =    &
     &   qrcev(mgs) + qhcev(mgs) + qscev(mgs) + qhlcev(mgs) + qfcev(mgs)
      pevap(mgs) =    &
     &   Min(0.0,qrcev(mgs)) + Min(0.0,qhcev(mgs)) + Min(0.0,qscev(mgs)) + Min(0.0,qhlcev(mgs)) &
         +  Min(0.0,qfcev(mgs))
      ! NOTE: pdep is the deposition part only
      pdep(mgs) =    &
     &   il5(mgs)*(   &
     &  + qsdpv(mgs) + qhdpv(mgs)   &
     &  + qhldpv(mgs)    &
     &  + qidpv(mgs)  )   & 
     &  +il5(mgs)*(qiint(mgs))
      elseif ( warmonly < 0.8 ) then
      pfrz(mgs) =    &
     &  (1-il5(mgs))*   &
     &  (qhmlr(mgs)+qhlmlr(mgs))   & !+qhmlh(mgs))   &
     &  +il5(mgs)*(qhfzh(mgs)+qhlfzhl(mgs))   &
     &  +il5(mgs)*(   &
     &  +qhshr(mgs)   &
     &  +qhlshr(mgs)   &
     &  +qrfrz(mgs)+qwfrz(mgs)   &
     &  +qwctfz(mgs)+qiihr(mgs)   &
     &  +qiacw(mgs)                &
     & +qhacw(mgs) + qhlacw(mgs)   &
     & +qhacr(mgs) + qhlacr(mgs)  ) 
      psub(mgs) =  0.0 +  &
     &   il5(mgs)*(   &
     &  + qhdpv(mgs)   &
     &  + qhldpv(mgs)    &
     &  + qidpv(mgs) + qisbv(mgs) )   &
     &  +il5(mgs)*(qiint(mgs))
      pvap(mgs) =    &
     &   qrcev(mgs) + qhcev(mgs) + qhlcev(mgs) + qfcev(mgs) 
      else
      pfrz(mgs) = 0.0
      psub(mgs) = 0.0
      pvap(mgs) = qrcev(mgs)
      endif ! warmonly
      ptem(mgs) =    &
     &  (1./pi0(mgs))*   &
     &  (felfcp(mgs)*pfrz(mgs)   &
     &  +felscp(mgs)*psub(mgs)    &
     &  +felvcp(mgs)*pvap(mgs))
      thetap(mgs) = thetap(mgs) + dtp*ptem(mgs)
      ptem2(mgs) = ptem(mgs)
      if ( eqtset > 2 ) then
        pipert(mgs) = pipert(mgs) + (felfpi(mgs)*pfrz(mgs)   &
     &  +felspi(mgs)*psub(mgs)    &
     &  +felvpi(mgs)*pvap(mgs))*dtp
      endif
      end do




!
!  sum the sources and sinks for qwvp, qcw, qci, qrw, qsw
!
!
      do mgs = 1,ngscnt


      qwvp(mgs) = qwvp(mgs) +        &
     &   dtp*(pqwvi(mgs)+pqwvd(mgs))
      qx(mgs,lc) = qx(mgs,lc) +   &
     &   dtp*(pqcwi(mgs)+pqcwd(mgs))
      qx(mgs,lr) = qx(mgs,lr) +   &
     &   dtp*(pqrwi(mgs)+pqrwd(mgs))
      qx(mgs,li) = qx(mgs,li) +   &
     &   dtp*(pqcii(mgs)+pqcid(mgs))
      qx(mgs,ls) = qx(mgs,ls) +   &
     &   dtp*(pqswi(mgs)+pqswd(mgs))
      qx(mgs,lh) = qx(mgs,lh) +    &
     &   dtp*(pqhwi(mgs)+pqhwd(mgs))

      if ( lhl .gt. 1 ) then
      qx(mgs,lhl) = qx(mgs,lhl) +    &
     &   dtp*(pqhli(mgs)+pqhld(mgs))
      endif


      end do

! sum sources for particle volume

      if ( ldovol ) then

      do mgs = 1,ngscnt

      if ( lvol(ls) .gt. 1 ) then
      vx(mgs,ls) = vx(mgs,ls) +    &
     &   dtp*(pvswi(mgs)+pvswd(mgs))
      endif

      if ( lvol(lh) .gt. 1 ) then
      vx(mgs,lh) = vx(mgs,lh) +    &
     &   dtp*(pvhwi(mgs)+pvhwd(mgs))
!     >   rho0(mgs)*dtp*(pqhwi(mgs)+pqhwd(mgs))/xdn0(lh)
      endif

      if ( lhl .gt. 1 ) then
      if ( lvol(lhl) .gt. 1 ) then
      vx(mgs,lhl) = vx(mgs,lhl) +    &
     &   dtp*(pvhli(mgs)+pvhld(mgs))
!     >   rho0(mgs)*dtp*(pqhwi(mgs)+pqhwd(mgs))/xdn0(lh)
      endif
      endif

      enddo

      endif  ! ldovol

!
!
!
! concentrations
!
      if ( ipconc .ge. 1  ) then
      do mgs = 1,ngscnt
      cx(mgs,li) = cx(mgs,li) +   &
     &   dtp*(pccii(mgs)+pccid(mgs)) 
      cina(mgs) = cina(mgs) + pccin(mgs)*dtp
      if ( ipconc .ge. 2 ) then
      cx(mgs,lc) = cx(mgs,lc) +   &
     &   dtp*(pccwi(mgs)+pccwd(mgs))
      endif
      if ( ipconc .ge. 3 ) then
      cx(mgs,lr) = cx(mgs,lr) +   &
     &   dtp*(pcrwi(mgs)+pcrwd(mgs))
      endif
      if ( ipconc .ge. 4 ) then
      cx(mgs,ls) = cx(mgs,ls) +   &
     &   dtp*(pcswi(mgs)+pcswd(mgs))
      endif
      if ( ipconc .ge. 5 ) then
      cx(mgs,lh) = cx(mgs,lh) +    &
     &   dtp*(pchwi(mgs)+pchwd(mgs))
       if ( lhl .gt. 1 ) then
        cx(mgs,lhl) = cx(mgs,lhl) +    &
     &     dtp*(pchli(mgs)+pchld(mgs))


        
        
       endif
      endif
      if ( ipconc .ge. 6 ) then
       if ( lzr .gt. 1 ) then
       zx(mgs,lr) = zx(mgs,lr) +    &
     &   dtp*(pzrwi(mgs)+pzrwd(mgs))
       endif
       if ( lzs .gt. 1 ) then
       zx(mgs,ls) = zx(mgs,ls) +    &
     &   dtp*(pzswi(mgs)+pzswd(mgs))
       endif
       if ( lzh .gt. 1 ) then
       zx(mgs,lh) = zx(mgs,lh) +    &
     &   dtp*(pzhwi(mgs)+pzhwd(mgs))
       endif
       if ( lzhl .gt. 1 ) then
        zx(mgs,lhl) = zx(mgs,lhl) +    &
     &     dtp*(pzhli(mgs)+pzhld(mgs))
!      IF ( pchli(mgs) .ne. 0. .or. pchld(mgs) .ne. 0 ) THEN
!       write(0,*) 'dr: cx,pchli,pchld = ', cx(mgs,lhl),pchli(mgs),pchld(mgs), igs(mgs),kgs(mgs)
!      ENDIF
       endif
      endif
      end do
      end if

!
!
!
! start saturation adjustment
!
      if (ndebug .gt. 0 ) write(0,*) 'conc 30a'
!      include 'sam.jms.satadj.sgi'
!
!
!
!  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
!
!
!
!  set up temperature and vapor arrays
!
      do mgs = 1,ngscnt
      pqs(mgs) = (380.0)/(pres(mgs))
      theta(mgs) = thetap(mgs) + theta0(mgs)
      qvap(mgs) = max( (qwvp(mgs) + qv0(mgs)), 0.0 )
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      end do
!
!  melting of cloud ice
!
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      ptimlw(mgs) = 0.0
      end do
!
      do mgs = 1,ngscnt
      qitmp(mgs) = qx(mgs,li)
      if( temg(mgs) .gt. tfr .and.   &
     &    qitmp(mgs) .gt. 0.0 ) then
      qx(mgs,lc) = qx(mgs,lc) + qitmp(mgs)
!      pfrz(mgs) = pfrz(mgs) - qitmp(mgs)*dtpinv
      ptem(mgs) =  ptem(mgs) +   &
     &  (1./pi0(mgs))*   &
     &  felfcp(mgs)*(- qitmp(mgs)*dtpinv)  
      if ( eqtset > 2 ) then
        pipert(mgs) = pipert(mgs) - (felfpi(mgs)*qitmp(mgs))
      endif
      pmlt(mgs) = pmlt(mgs) - qitmp(mgs)*dtpinv
      scx(mgs,lc) = scx(mgs,lc) + scx(mgs,li)
      thetap(mgs) = thetap(mgs) -   &
     &  fcc3(mgs)*qitmp(mgs)
      ptimlw(mgs) = -fcc3(mgs)*qitmp(mgs)*dtpinv
      cx(mgs,lc) = cx(mgs,lc) + cx(mgs,li)
      qx(mgs,li) = 0.0
      cx(mgs,li) = 0.0
      scx(mgs,li) = 0.0
      vx(mgs,li) = 0.0
      qitmp(mgs) = 0.0
      end if
      end do

!
!


!      do mgs = 1,ngscnt
!      qimlw(mgs) = (qcwtmp(mgs)-qx(mgs,lc))*dtpinv
!      end do
!
!  homogeneous freezing of cloud water
!
      if ( warmonly < 0.8 ) then

      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      ptwfzi(mgs) = 0.0
      end do
!
      do mgs = 1,ngscnt

!      if( temg(mgs) .lt. tfrh ) THEN
!       write(0,*) 'GS: mgs,temp,qc,qi = ',mgs,temg(mgs),temcg(mgs),qx(mgs,lc),qx(mgs,li)
!      ENDIF

      ctmp = 0.0
      frac = 0.0
      qtmp = 0.0
      
!      if( ( temg(mgs) .lt. thnuc + 2. .or. (ibfc == 2 .and. temg(mgs) < thnuc + 10. ) ) .and.    &
!     &  qx(mgs,lc) .gt. qxmin(lc) .and. (ipconc < 2 .or. ibfc == 0 .or. ibfc == 2 )) then
! commented for test (12/01/2015):
!      if( temg(mgs) .lt. thnuc + 0. .and.    &
!     &  qx(mgs,lc) .gt. 0.0 .and. (ipconc < 2 .or. ibfc == 0 )) then
      if( ( ( temg(mgs) .lt. thnuc + 0.) .or. (temg(mgs) .lt. thnuc + 2. .and. ibfc >= 3) ) .and.    &
     &  qx(mgs,lc) .gt. 0.0 .and. (ipconc < 2 .or. ibfc == 0 .or. ibfc == 2)) then

      if ( ibfc >= 3 ) then
        frac = Max( 0.25, Min( 1., ((thnuc + 2.) - temg(mgs) )/4.0 ) )
      elseif ( ibfc /= 2 .or. ipconc < 2 ) then
        frac = Max( 0.25, Min( 1., ((thnuc + 1.) - temg(mgs) )/4.0 ) )
      else
          volt = exp( 16.2 + 1.0*temcg(mgs) )* 1.0e-6 !  Ts == -temcg ; volt comes from the fit in Fig. 1 in Bigg 1953 
                                               ! for mean temperature for freezing: -ln (V) = a*Ts - b
                                               ! volt is given in cm**3, so factor of 1.e-6 to convert to m**3
         
         cwfrz(mgs) = cx(mgs,lc)*Exp(-volt/xv(mgs,lc)) ! number of droplets with volume greater than volt

         qtmp = cwfrz(mgs)*xdn0(lc)*rhoinv(mgs)*(volt + xv(mgs,lc))
         frac = qtmp/qx(mgs,lc) ! reset number frozen to same fraction as mass. This makes 
                                                       ! sure that cwfrz and qwfrz are consistent and prevents 
                                                       ! spurious creation of ice crystals.
      
      endif
      qtmp = frac*qx(mgs,lc)

      if ( ibfc == 4 .and. lis >= 1 ) then
        qx(mgs,lis) = qx(mgs,lis) + qtmp
      else
        qx(mgs,li) = qx(mgs,li) + qtmp ! qx(mgs,lc)
      endif
      pfrz(mgs) = pfrz(mgs) + qtmp*dtpinv
      ptem(mgs) =  ptem(mgs) +   &
     &  (1./pi0(mgs))*   &
     &  felfcp(mgs)*(qtmp*dtpinv)  

      if ( eqtset > 2 ) then
        pipert(mgs) = pipert(mgs) + felfpi(mgs)*qtmp
      endif

!      IF ( lvol(li) .gt. 1 ) vx(mgs,li) = vx(mgs,li) + rho0(mgs)*qx(mgs,lc)/xdn0(li)
      if ( lvol(li) .gt. 1 ) vx(mgs,li) = vx(mgs,li) + rho0(mgs)*qtmp/xdn0(li)

      if ( ipconc .ge. 2 ) then
        ctmp = frac*cx(mgs,lc)
!        cx(mgs,li) = cx(mgs,li) + cx(mgs,lc)
        if ( ibfc == 4 .and. lis >= 1 ) then
          cx(mgs,lis) = cx(mgs,lis) + ctmp
        else
          cx(mgs,li) = cx(mgs,li) + ctmp
        endif
      else ! (ipconc .lt. 2 )
        ctmp = 0.0
        if ( t9(igs(mgs),jgs,kgs(mgs)-1) .gt. qx(mgs,lc) ) then
           qtmp = frac*t9(igs(mgs),jgs,kgs(mgs)-1)  

!           cx(mgs,lc) = cx(mgs,lc)*qx(mgs,lc)*rho0(mgs)/qtmp
           ctmp = cx(mgs,lc)*qx(mgs,lc)*rho0(mgs)/qtmp
        else
           cx(mgs,lc) = Max(0.0,wvel(mgs))*dtp*cwccn   &
     &      /gz(igs(mgs),jgs,kgs(mgs))
          cx(mgs,lc) = cwccn
        endif

       if ( ipconc .ge. 1 ) cx(mgs,li) = Min(ccimx, cx(mgs,li) + cx(mgs,lc))
      endif

      sctmp = frac*scx(mgs,lc)
!      scx(mgs,li) = scx(mgs,li) + scx(mgs,lc)
      scx(mgs,li) = scx(mgs,li) + sctmp
!      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qx(mgs,lc)
!      ptwfzi(mgs) = fcc3(mgs)*qx(mgs,lc)*dtpinv
!      qx(mgs,lc) = 0.0
!      cx(mgs,lc) = 0.0
!      scx(mgs,lc) = 0.0
      thetap(mgs) = thetap(mgs) + fcc3(mgs)*qtmp
      ptwfzi(mgs) = fcc3(mgs)*qtmp*dtpinv
      qx(mgs,lc) = qx(mgs,lc) - qtmp
      cx(mgs,lc) = cx(mgs,lc) - ctmp
      scx(mgs,lc) = scx(mgs,lc) - sctmp
      end if
      end do

      endif ! warmonly
!
!      do mgs = 1,ngscnt
!      qwfzi(mgs) = (qcwtmp(mgs)-qx(mgs,lc))*dtpinv   ! Not used?? (ERM)
!      end do
!
!  reset temporaries for cloud particles and vapor
!
      qcond(:) = 0.0
      
      if ( ipconc .le. 1 .and.  lwsm6 ) then ! Explicit cloud condensation/evaporation (Rutledge and Hobbs 1983)
       do mgs = 1,ngscnt

        qcwtmp(mgs) = qx(mgs,lc)
        theta(mgs) = thetap(mgs) + theta0(mgs)
        temgtmp = temg(mgs)
!        temg(mgs) = theta(mgs)*(p2(igs(mgs),jgs,kgs(mgs)) ) ! *pk(mgs) ! ( pres(mgs) / poo ) ** cap
!        temsav = temg(mgs)
!        thsave(mgs) = thetap(mgs)
        temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
        temcg(mgs) = temg(mgs) - tfr
        ltemq = (temg(mgs)-163.15)/fqsat+1.5
        ltemq = Min( nqsat, Max(1,ltemq) )

        qvs(mgs) = pqs(mgs)*tabqvs(ltemq)

        if ( ( qvap(mgs) > qvs(mgs) .or. qx(mgs,lc) > qxmin(lc) ) .and. temg(mgs) > tfrh ) then
          tmp = (qvap(mgs) - qvs(mgs))/(1. + qvs(mgs)*felv(mgs)**2/(cp*rw*temg(mgs)**2) )
          qcond(mgs) = Min( Max( 0.0, tmp ), (qvap(mgs)-qvs(mgs)) )
          if ( qx(mgs,lc) > qxmin(lc) .and. tmp < 0.0 ) then ! evaporation
            qcond(mgs) = Max( tmp, -qx(mgs,lc) )
          endif
          qwvp(mgs) = qwvp(mgs) - qcond(mgs)
          qvap(mgs) = qvap(mgs) - qcond(mgs)
          qx(mgs,lc) = Max( 0.0, qx(mgs,lc) + qcond(mgs) )
          thetap(mgs) = thetap(mgs) + felvcp(mgs)*qcond(mgs)/(pi0(mgs))
          
        endif
        
        enddo
      
      endif
      
      
      if ( ipconc .le. 1 .and. .not. lwsm6 ) then
!      IF ( ipconc .le. 1  ) THEN
      
      do mgs = 1,ngscnt
      qx(mgs,lv) = max( 0.0, qvap(mgs) )
      qx(mgs,lc) = max( 0.0, qx(mgs,lc) )
      qx(mgs,li) = max( 0.0, qx(mgs,li) )
      qitmp(mgs) = qx(mgs,li)
      end do
!
!
      do mgs = 1,ngscnt
      qcwtmp(mgs) = qx(mgs,lc)
      qitmp(mgs) = qx(mgs,li)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temgtmp = temg(mgs)
      temg(mgs) = theta(mgs)*(pinit(kgs(mgs)) + p2(igs(mgs),jgs,kgs(mgs)) ) ! *pk(mgs) ! ( pres(mgs) / poo ) ** cap
      temsav = temg(mgs)
      thsave(mgs) = thetap(mgs)
      temcg(mgs) = temg(mgs) - tfr
      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )

      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
      qss(mgs) = qvs(mgs)
      if ( temg(mgs) .lt. tfr ) then
      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )   &
     &  qss(mgs) = qvs(mgs)
      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &  qss(mgs) = qis(mgs)
      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /   &
     &   (qx(mgs,lc) + qitmp(mgs))
      end if
      end do
!
!  iterate  adjustment
!
      do itertd = 1,2
!
      do mgs = 1,ngscnt
!
!  calculate super-saturation
!
      qitmp(mgs) = qx(mgs,li)
      fcci(mgs) = 0.0
      fcip(mgs) = 0.0
      dqcw(mgs) = 0.0
      dqci(mgs) = 0.0
      dqwv(mgs) = ( qx(mgs,lv) - qss(mgs) )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then           !  subsaturated
        if( qx(mgs,lc) .gt. -dqwv(mgs) ) then  ! check if qc can make up all of the deficit
          dqcw(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                 !  otherwise make all qc available for evap
          dqcw(mgs) = -qx(mgs,lc)
          dqwv(mgs) = dqwv(mgs) + qx(mgs,lc)
        end if
!
        if( qitmp(mgs) .gt. -dqwv(mgs) ) then  ! check if qi can make up all the deficit
          dqci(mgs) = dqwv(mgs)
          dqwv(mgs) = 0.
        else                                  ! otherwise make all ice available for sublimation
          dqci(mgs) = -qitmp(mgs)
          dqwv(mgs) = dqwv(mgs) + qitmp(mgs)
        end if
!
       qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs) + dqci(mgs) )  ! add to perturbation vapor
!
! This next line removed 3/19/2003 thanks to Adam Houston,
!  who found the bug in the 3-ICE code
!      qwvp(mgs) = max(qwvp(mgs), 0.0)
      qitmp(mgs) = qx(mgs,li)
      if ( qitmp(mgs) .ge. qxmin(li) ) then
        fcci(mgs) = qx(mgs,li)/(qitmp(mgs))
      else
        fcci(mgs) = 1.0
      endif
      qx(mgs,lc) = qx(mgs,lc) + dqcw(mgs)
      qx(mgs,li) = qx(mgs,li) + dqci(mgs) * fcci(mgs)
      thetap(mgs) = thetap(mgs) +   &
     &  1./pi0(mgs)*   &
     &  (felvcp(mgs)*dqcw(mgs) +felscp(mgs)*dqci(mgs))

      if ( eqtset > 2 ) then
        pipert(mgs) = pipert(mgs)   &
     &  +(felspi(mgs)*dqci(mgs)    &
     &  +felvpi(mgs)*dqcw(mgs))*dtp
      endif

      end if  ! dqwv(mgs) .lt. 0. (end of evap/sublim)
!
! condensation/deposition
!
      if ( dqwv(mgs) .ge. 0. ) then
      
!      write(iunit,*) 'satadj: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qx(mgs,lv),qx(mgs,lc)
!
        qitmp(mgs) = qx(mgs,li)
        fracl(mgs) = 1.0
        fraci(mgs) = 0.0
        if ( temg(mgs) .lt. tfr .and. temg(mgs) .gt. thnuc ) then
          fracl(mgs) = max(min(1.,(temg(mgs)-233.15)/(20.)),0.0)
          fraci(mgs) = 1.0-fracl(mgs)
        end if
        if ( temg(mgs) .le. thnuc ) then
           fraci(mgs) = 1.0
           fracl(mgs) = 0.0
         end if
        fraci(mgs) = 1.0-fracl(mgs)
!
       gamss = (felvcp(mgs)*fracl(mgs) + felscp(mgs)*fraci(mgs))   &
     &      / (pi0(mgs))
!
      if ( temg(mgs) .lt. tfr ) then
        if (qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) ) then
         dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/   &
     &  ((temg(mgs)-cbw)**2))
        end if
        if ( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li) ) then
          dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv2(mgs)*qss(mgs)/   &
     &  ((temg(mgs)-cbi)**2))
        end if
        if ( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li) ) then
         cdw = caw*pi0(mgs)*tfrcbw/((temg(mgs)-cbw)**2)
         cdi = cai*pi0(mgs)*tfrcbi/((temg(mgs)-cbi)**2)
         denom1 = qx(mgs,lc) + qitmp(mgs)
         denom2 = 1.0 + gamss*   &
     &    (qx(mgs,lc)*qvs(mgs)*cdw + qitmp(mgs)*qis(mgs)*cdi) / denom1
         dqvcnd(mgs) =  dqwv(mgs) / denom2
        end if 

      endif  !  temg(mgs) .lt. tfr
!
      if ( temg(mgs) .ge. tfr ) then
      dqvcnd(mgs) = dqwv(mgs)/(1. + fcqv1(mgs)*qss(mgs)/   &
     &  ((temg(mgs)-cbw)**2))
      end if
!
      delqci1=qx(mgs,li)
!
      if ( qitmp(mgs) .gt. qxmin(li) ) then
        fcci(mgs) = qx(mgs,li)/(qitmp(mgs))
      else
        fcci(mgs) = 1.0
      endif
!
      dqcw(mgs) = dqvcnd(mgs)*fracl(mgs)
      dqci(mgs) = dqvcnd(mgs)*fraci(mgs)
!
      thetap(mgs) = thetap(mgs) +   &
     &   (felvcp(mgs)*dqcw(mgs) + felscp(mgs)*dqci(mgs))   &
     & / (pi0(mgs))

      if ( eqtset > 2 ) then
        pipert(mgs) = pipert(mgs) + (0   &
     &  +felspi(mgs)*dqci(mgs)    &
     &  +felvpi(mgs)*dqcw(mgs))*dtp
      endif

      qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )
      qx(mgs,lc) = qx(mgs,lc) + dqcw(mgs)
!      IF ( qitmp(mgs) .gt. qxmin(li) ) THEN
        qx(mgs,li) = qx(mgs,li) + dqci(mgs)*fcci(mgs)
        qitmp(mgs) = qx(mgs,li)
!      ENDIF
!
!      delqci(mgs) =  dqci(mgs)*fcci(mgs)
!
      end if !  dqwv(mgs) .ge. 0.
      end do
!
      do mgs = 1,ngscnt
      qitmp(mgs) = qx(mgs,li)
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temg(mgs) = theta(mgs)*pk(mgs) ! ( pres(mgs) / poo ) ** cap
      qvap(mgs) = Max((qwvp(mgs) + qv0(mgs)), 0.0)
      temcg(mgs) = temg(mgs) - tfr
      tqvcon = temg(mgs)-cbw
      ltemq = (temg(mgs)-163.15)/fqsat+1.5
      ltemq = Min( nqsat, Max(1,ltemq) )
      qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
      qis(mgs) = pqs(mgs)*tabqis(ltemq)
      qx(mgs,lc) = max( 0.0, qx(mgs,lc) )
      qitmp(mgs) = max( 0.0, qitmp(mgs) )
      qx(mgs,lv) = max( 0.0, qvap(mgs))
!      if ( temg(mgs) .lt. tfr ) then
!      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )
!     >  qss(mgs) = qvs(mgs)
!c      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
!      if( qx(mgs,lc) .eq. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
!     >  qss(mgs) = qis(mgs)
!c      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))
!      if( qx(mgs,lc) .gt. 0.0 .and. qitmp(mgs) .gt. qxmin(li))
!     >  qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /
!     > (qx(mgs,lc) + qitmp(mgs))
!      else
!      qss(mgs) = qvs(mgs)
!      end if
      qss(mgs) = qvs(mgs)
      if ( temg(mgs) .lt. tfr ) then
      if( qx(mgs,lc) .ge. 0.0 .and. qitmp(mgs) .le. qxmin(li) )   &
     &  qss(mgs) = qvs(mgs)
      if( qx(mgs,lc) .le. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &  qss(mgs) = qis(mgs)
      if( qx(mgs,lc) .gt. qxmin(lc) .and. qitmp(mgs) .gt. qxmin(li))   &
     &   qss(mgs) = (qx(mgs,lc)*qvs(mgs) + qitmp(mgs)*qis(mgs)) /   &
     &   (qx(mgs,lc) + qitmp(mgs))
      end if
!      pceds(mgs) = (thetap(mgs) - thsave(mgs))*dtpinv
!      write(iunit,*) 'satadj2: mgs,iter = ',mgs,itertd,dqwv(mgs),qss(mgs),qx(mgs,lv),qx(mgs,lc)
      end do
!
!  end the saturation adjustment iteration loop
!
      end do

     endif ! ( ipconc .le. 1 )

!
!  spread the growth owing to vapor diffusion onto the
!  ice crystal categories using the
!
!  END OF SATURATION ADJUSTMENT
!

      if (ndebug .gt. 0 ) write(0,*) 'conc 30b'
!
!
!  end of saturation adjustment

!
!
! !DIR$ IVDEP
      do mgs = 1,ngscnt
      t0(igs(mgs),jy,kgs(mgs)) =  temg(mgs)
      end do
!
! Load the save arrays
!
      if ( numproc > 1 ) then
      do mgs = 1,ngscnt
       dv = dx1*dy1*gz(igs(mgs),1,kgs(mgs))
       if ( ipconc > 2 ) then
       thproc(kzbeg-1+kgs(mgs),1) = thproc(kzbeg-1+kgs(mgs),1) + crfrzf(mgs)*dtp*dv
       else
       thproc(kzbeg-1+kgs(mgs),1) = thproc(kzbeg-1+kgs(mgs),1) + qrfrzf(mgs)*rho0(mgs)*dtp*dv
       endif
       thproc(kzbeg-1+kgs(mgs),2) = thproc(kzbeg-1+kgs(mgs),2) + il5(mgs)*ciacrf(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),3) = thproc(kzbeg-1+kgs(mgs),3) + chcnsh(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),4) = thproc(kzbeg-1+kgs(mgs),4) + chcnih(mgs)*dtp*dv
       if (  qhacw(mgs)+qhacr(mgs) > 0.0 .and. temg(mgs) < tfr ) then
       thproc(kzbeg-1+kgs(mgs),5) = thproc(kzbeg-1+kgs(mgs),5) + (qhacw(mgs)+qhacr(mgs)+qhshr(mgs))*rho0(mgs)*dtp*dv
       endif
       thproc(kzbeg-1+kgs(mgs),6) = thproc(kzbeg-1+kgs(mgs),6) + qracw(mgs)*rho0(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),7) = thproc(kzbeg-1+kgs(mgs),7) + qrcnw(mgs)*rho0(mgs)*dtp*dv
       if ( qhacw(mgs) > 0.0 .and. temg(mgs) < tfr ) then
       thproc(kzbeg-1+kgs(mgs),8) = thproc(kzbeg-1+kgs(mgs),8) + (vhacw(mgs)+vhacr(mgs)+vhshdr(mgs))*dtp*dv
!       thproc(kzbeg-1+kgs(mgs),8) = thproc(kzbeg-1+kgs(mgs),8) + qhacw(mgs)*rho0(mgs)/rimdn(mgs,lh)*dtp*dv
       endif
       thproc(kzbeg-1+kgs(mgs),9) = thproc(kzbeg-1+kgs(mgs),9) + pi0(mgs)*ptem(mgs)*dtp*dv  ! latent heating -- mult by pi0 to get T rather than theta
       thproc(kzbeg-1+kgs(mgs),10) = thproc(kzbeg-1+kgs(mgs),10) +  &
     &                             ( chmul1(mgs) + chlmul1(mgs)  )*dtp*dv
       if ( lf > 1 ) then
       else
       thproc(kzbeg-1+kgs(mgs),11) = thproc(kzbeg-1+kgs(mgs),11) +  &
     &                             ( csplinter(mgs) + csplinter2(mgs) )*dtp*dv
       endif
       thproc(kzbeg-1+kgs(mgs),12) = thproc(kzbeg-1+kgs(mgs),12) + qrfrzf(mgs)*rho0(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),13) = thproc(kzbeg-1+kgs(mgs),13) + il5(mgs)*qiacrf(mgs)*rho0(mgs)*dtp*dv ! mass of rain freezing by ice crystal capture
       thproc(kzbeg-1+kgs(mgs),14) = thproc(kzbeg-1+kgs(mgs),14) + crcnw(mgs)*dtp*dv    ! rain drop prod. by autoconv.
       thproc(kzbeg-1+kgs(mgs),15) = thproc(kzbeg-1+kgs(mgs),15) + (pcrwi(mgs)-crcnw(mgs))*dtp*dv ! rain drop prod by melting/shedding (i.e., everything but autoconv.)
!       thproc(kzbeg-1+kgs(mgs),18) = thproc(kzbeg-1+kgs(mgs),18) + pevap(mgs)*rho0(mgs)*dv ! rain evaporation rate
       thproc(kzbeg-1+kgs(mgs),19) = thproc(kzbeg-1+kgs(mgs),19) + pmlt(mgs)*rho0(mgs)*dv  ! melting rate
       thproc(kzbeg-1+kgs(mgs),20) = thproc(kzbeg-1+kgs(mgs),20) + pdep(mgs)*rho0(mgs)*dv  ! deposition rate
       thproc(kzbeg-1+kgs(mgs),21) = thproc(kzbeg-1+kgs(mgs),21) + (psub(mgs)-pdep(mgs))*rho0(mgs)*dv  ! sublimation rate
       thproc(kzbeg-1+kgs(mgs),22) = thproc(kzbeg-1+kgs(mgs),22) + (pfrz(mgs)-pmlt(mgs))*rho0(mgs)*dv  ! freezing rate

!       thproc(kzbeg-1+kgs(mgs),20) = thproc(kzbeg-1+kgs(mgs),20) + (1./pi0(mgs))*felfcp(mgs)*pvap(mgs)*rho0(mgs)*dv  ! deposition rate
!       thproc(kzbeg-1+kgs(mgs),21) = thproc(kzbeg-1+kgs(mgs),21) + (1./pi0(mgs))*felscp(mgs)*psub(mgs)*rho0(mgs)*dv  ! sublimation rate
!       thproc(kzbeg-1+kgs(mgs),22) = thproc(kzbeg-1+kgs(mgs),22) + (1./pi0(mgs))*felfcp(mgs)*pfrz(mgs)*rho0(mgs)*dv ! (pfrz(mgs)-pmlt(mgs))*rho0(mgs)*dv  ! freezing rate

       thproc(kzbeg-1+kgs(mgs),23) = thproc(kzbeg-1+kgs(mgs),23) + crfrzs(mgs)*dtp*dv
       thproc(kzbeg-1+kgs(mgs),24) = thproc(kzbeg-1+kgs(mgs),24) + il5(mgs)*ciacrs(mgs)*dtp*dv

       thproc(kzbeg-1+kgs(mgs),25) = thproc(kzbeg-1+kgs(mgs),25) + qhmlr(mgs)*rho0(mgs)*dv  ! melting rate
       thproc(kzbeg-1+kgs(mgs),26) = thproc(kzbeg-1+kgs(mgs),26) + qhlmlr(mgs)*rho0(mgs)*dv  ! melting rate

       if (  qhlacw(mgs)+qhlacr(mgs) > 0.0 .and. temg(mgs) < tfr ) then
        thproc(kzbeg-1+kgs(mgs),27) = thproc(kzbeg-1+kgs(mgs),27) + (qhlacw(mgs)+qhlacr(mgs)+qhlshr(mgs))*rho0(mgs)*dtp*dv
        thproc(kzbeg-1+kgs(mgs),28) = thproc(kzbeg-1+kgs(mgs),28) + (qhlacw(mgs))*rho0(mgs)*dtp*dv
        thproc(kzbeg-1+kgs(mgs),29) = thproc(kzbeg-1+kgs(mgs),29) + (qhlacr(mgs))*rho0(mgs)*dtp*dv
       endif

       if ( temg(mgs) < tfr ) then
        thproc(kzbeg-1+kgs(mgs),30) = thproc(kzbeg-1+kgs(mgs),30) + (qhacw(mgs))*rho0(mgs)*dtp*dv
        thproc(kzbeg-1+kgs(mgs),31) = thproc(kzbeg-1+kgs(mgs),31) + (qhacr(mgs))*rho0(mgs)*dtp*dv
       endif

       thproc(kzbeg-1+kgs(mgs),32) = thproc(kzbeg-1+kgs(mgs),32) + qhlcnh(mgs)*rho0(mgs)*dtp*dv ! graupel mass conversion to hail

       if ( ihrn > 0 ) then
       thproc(kzbeg-1+kgs(mgs),33) = thproc(kzbeg-1+kgs(mgs),33) + ciihr(mgs)*dtp*dv ! contact freezing of droplets
       else
       if ( qwctfz(mgs)*dtp >= qxmin(li) ) then
       thproc(kzbeg-1+kgs(mgs),33) = thproc(kzbeg-1+kgs(mgs),33) + cwctfz(mgs)*dtp*dv ! contact freezing of droplets
       endif
       endif
       thproc(kzbeg-1+kgs(mgs),34) = thproc(kzbeg-1+kgs(mgs),34) + pevap(mgs)*rho0(mgs)*dv ! rain evaporation rate
       if ( qiint(mgs)*dtp >= qxmin(li) ) then
       thproc(kzbeg-1+kgs(mgs),35) = thproc(kzbeg-1+kgs(mgs),35) + ciint(mgs)*dtp*dv ! primary ice initiation
       endif
       if ( lf > 1 ) then
       else
          thproc(kzbeg-1+kgs(mgs),38) = thproc(kzbeg-1+kgs(mgs),38) + (vhacw(mgs)+vhacr(mgs)+vhshdr(mgs))*dtp*dv
        endif
      if ( lhl > 1 ) then
        thproc(kzbeg-1+kgs(mgs),36) = thproc(kzbeg-1+kgs(mgs),36) + chlcnhhl(mgs)*dtp*dv
          thproc(kzbeg-1+kgs(mgs),37) = thproc(kzbeg-1+kgs(mgs),37) + (vhlacw(mgs)+vhlacr(mgs)+vhlshdr(mgs))*dtp*dv
       else
         if ( lf > 1 ) then 
         else
           thproc(kzbeg-1+kgs(mgs),36) = thproc(kzbeg-1+kgs(mgs),36) + (pchwi(mgs))*dtp*dv
           thproc(kzbeg-1+kgs(mgs),39) = thproc(kzbeg-1+kgs(mgs),39) + (pchwd(mgs))*dtp*dv
           thproc(kzbeg-1+kgs(mgs),37) = thproc(kzbeg-1+kgs(mgs),37) + (chmlr(mgs))*dtp*dv
         endif
       endif

!        IF ( numproc >= 44 ) THEN
          thproc(kzbeg-1+kgs(mgs),42) = thproc(kzbeg-1+kgs(mgs),42) - qhshr(mgs)*dtp*dv
          thproc(kzbeg-1+kgs(mgs),43) = 0.0
          thproc(kzbeg-1+kgs(mgs),44) = thproc(kzbeg-1+kgs(mgs),44) - qhlshr(mgs)*dtp*dv
!        ENDIF
       
       
!       thproc(kzbeg-1+kgs(mgs),35) = thproc(kzbeg-1+kgs(mgs),35) + pevap(mgs)*rho0(mgs)*dv ! rain evaporation rate


!      ptem(mgs) =    &
!     &  (1./pi0(mgs))*   &
!     &  (felfcp(mgs)*pfrz(mgs)   &
!     &  +felscp(mgs)*psub(mgs)    &
!     &  +felvcp(mgs)*pvap(mgs))

      enddo
      endif


! Sample code for using the axtra array to load microphysical rates or quantities for output
!
! Note that indices 1 and 2 are used in the nucond subroutine for condensation/evap of droplets (1) and
!    condensation of rain (2)
!
!      IF ( io_flag .and. nxtra > 1 ) THEN
!        DO mgs = 1,ngscnt
!          axtra(igs(mgs),jy,kgs(mgs),3)  = pfrz(mgs) !
!          axtra(igs(mgs),jy,kgs(mgs),4)  = qrcev(mgs) ! pre2
!          axtra(igs(mgs),jy,kgs(mgs),5)  = psub(mgs) ! depsubr
!          axtra(igs(mgs),jy,kgs(mgs),6)  = qrfrz(mgs) ! rain freezing (Bigg)
!          axtra(igs(mgs),jy,kgs(mgs),7)  = pmlt(mgs) ! melr2
!        ENDDO
!      ENDIF



      if (ndebug .gt. 0 ) write(0,*) 'gs 11'

      do mgs = 1,ngscnt
!
      an(igs(mgs),jy,kgs(mgs),lt) =    &
     &  theta0(mgs) + thetap(mgs) 
      an(igs(mgs),jy,kgs(mgs),lv) = qwvp(mgs) + qv0(mgs) !

      if ( eqtset > 2 ) then
        p2(igs(mgs),jy,kgs(mgs)) = pipert(mgs)
      endif
!
      
      do il = lc,lhab
        if ( ido(il) .eq. 1 ) then
        if ( lf > 1 .and. il == lf ) then 
           lfsave(mgs,1) = an(igs(mgs),jy,kgs(mgs),il)
           lfsave(mgs,2) = qx(mgs,il)
        endif
         an(igs(mgs),jy,kgs(mgs),il) = qx(mgs,il) +   &
     &     min( an(igs(mgs),jy,kgs(mgs),il), 0.0 )
         qx(mgs,il) = an(igs(mgs),jy,kgs(mgs),il)
        endif
      enddo

      if ( lcina > 1 ) then
        an(igs(mgs),jy,kgs(mgs),lcina) = cina(mgs)
      endif





!
!  6th moments
!

      if ( ipconc .ge. 6 ) then
       do il = lr,lhab
        if ( lz(il) .gt. 1 ) then
        if ( lf > 1 .and. il == lf ) then 
           lfsave(mgs,3) = an(igs(mgs),jy,kgs(mgs),lz(il))
           lfsave(mgs,4) = zx(mgs,il)
        endif

         an(igs(mgs),jy,kgs(mgs),lz(il)) = zx(mgs,il) +   &
     &     min( an(igs(mgs),jy,kgs(mgs),lz(il)), 0.0 )
         zx(mgs,il) = an(igs(mgs),jy,kgs(mgs),lz(il))
         
        endif
       enddo
       
      endif
!
      end do
!

      if ( ipconc .ge. 1 ) then
      do il = lc,lhab !{

!        write(0,*) 'limiter loop: il,ipc,lz: ',il,ipc(il),lz(il),ipconc

       if ( ipconc .ge. ipc(il) .and. ido(il) > 0 ) then ! {

         if (  ipconc .ge. 4 .and. ipc(il) .ge. 1 ) then ! {

!            write(0,*) 'MY limiter: il,ipc,lz: ',il,ipc(il),lz(il),lr,lzr
!            STOP

          if ( lz(il) <= 1 .or. ioldlimiter == 1 ) then ! { { is a two-moment category so dont worry about reflectivity
          

           do mgs = 1,ngscnt
            if ( qx(mgs,il) .le. 0.0 ) then
              cx(mgs,il) = 0.0
            else !{
              if ( cx(mgs,il) .gt. cxmin ) then !{
!              xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
!              xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(cxmin,cx(mgs,il)))
                xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*cx(mgs,il))
              
!              IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!               write(0,*) 'dr: xv,cx,qx,xdn,ln = ',xv(mgs,il),cx(mgs,il),qx(mgs,il),xdn(mgs,il),ln(il)
!              ENDIF

               ! 8/26/2015 erm: apply imaxdiaopt for 2-moment also
               if ( imaxdiaopt == 1 .or. il == lc .or. il == li .or. (il == lr .and. imurain == 3) .or. &
     &              (il == ls .and. imusnow == 3 ) ) then
                 xvbarmax = xvmx(il)
               elseif ( imaxdiaopt == 2 ) then ! test against maximum mass diameter
                 xvbarmax = xvmx(il) /((3. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
               elseif ( imaxdiaopt == 3 ) then ! test against mass-weighted diameter
                 xvbarmax = xvmx(il) /((4. + alpha(mgs,il))**3/((3. + alpha(mgs,il))*(2. + alpha(mgs,il))*(1. + alpha(mgs,il))))
               else
                 xvbarmax = xvmx(il)
               endif

               tmp = 1.0
               if ( il == ls ) then
                 xvbarmax = xvbarmax*Max(1.,100./Min(100.,xdn(mgs,ls)))
               endif
               
               if ( xv(mgs,il) .lt. xvmn(il) .or. xv(mgs,il) .gt. xvbarmax ) then
                xv(mgs,il) = Min( xvbarmax, xv(mgs,il) )
                xv(mgs,il) = Max( xvmn(il), xv(mgs,il) )
                cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xv(mgs,il)*xdn(mgs,il))
               endif
              
             endif !}

!              IF ( lhl .gt. 1 .and. il .eq. lhl ) THEN
!               write(0,*) 'dr: xv,cx,= ',xv(mgs,il),cx(mgs,il)
!              ENDIF

            endif !}
           enddo ! mgs
          
          else ! } { is three-moment, so have to adjust Z if size is too large
           if ( il == lr .and. imurain == 3 ) then ! { { RAIN

!          rdmx = 
!          rdmn = 

          do mgs = 1,ngscnt
          

         if ( iresetmoments == 1 .or. iresetmoments == il  ) then
         if ( zx(mgs,lr) <= zxmin ) then
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           qx(mgs,lr) = 0.0
           cx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),ln(lr)) = cx(mgs,lr)
         elseif ( cx(mgs,lr) <= cxmin ) then
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,lr) = 0.0
           qx(mgs,lr) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),lr)
           an(igs(mgs),jgs,kgs(mgs),lr) = qx(mgs,lr)
           an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)
         endif
         endif
         
         if ( qx(mgs,lr) .gt. qxmin(lr) ) then

        xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xdn(mgs,lr)*Max(1.0e-11,cx(mgs,lr)))
        if ( xv(mgs,lr) .gt. xvmx(lr) ) then
!          xv(mgs,lr) = xvmx(lr)
!          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmx(lr)*xdn(mgs,lr))
        elseif ( xv(mgs,lr) .lt. xvmn(lr) ) then
          xv(mgs,lr) = xvmn(lr)
          cx(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(xvmn(lr)*xdn(mgs,lr))
        endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*xdn(mgs,lr)**2)
!            an(igs(mgs),jgs,kgs(mgs),ln(il)) = zx(mgs,il)
           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(xdn(mgs,lr)**2*chw)
            an(igs(mgs),jgs,kgs(mgs),lz(lr)) = zx(mgs,lr)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/(z*1000.*1000)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
          endif
        
          if ( zx(mgs,lr) > 0.0 ) then
            xv(mgs,lr) = rho0(mgs)*qx(mgs,lr)/(1000.*cx(mgs,lr))
            vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)

!           xv = (db(1,kz)*a(1,1,kz,lr))**2/(a(1,1,kz,lnr))
!           rd = z*(pi/6.*1000.)**2/xv

! determine shape parameter alpha by iteration
           if ( z .gt. 0.0 ) then
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           do i = 1,20
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo

! check for artificial breakup (rain larger than allowed max size)
        if (  xv(mgs,il) .gt. xvmx(il) .or. (ioldlimiter == 2 .and. xv(mgs,il) .gt. xvmx(il)/8.) ) then
          tmp = cx(mgs,il)
!            write(0,*) 'MY limiter: xv: ',xv(mgs,il), xv(mgs,il)/(xvmx(il)/8.)
!            STOP
          if ( ioldlimiter == 2 ) then ! MY-style active breakup
            x = (6.*rho0(mgs)*qx(mgs,il)/(pi*xdn(mgs,il)*cx(mgs,il)))**(1./3.)
            x1 = Max(0.0e-3, x - 3.0e-3)
            x2 = Max(0.5, x/6.0e-3)
            x3 = x2**3
            cx(mgs,il) = cx(mgs,il)*Max((1.+2.222e3*x1**2), x3)
            xv(mgs,il) = xv(mgs,il)/Max((1.+2.222e3*x1**2), x3)
          else ! simple cutoff 
            xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
            xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          endif
            !xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
            !cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          
          
          if ( tmp < cx(mgs,il) ) then ! breakup

            g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
            zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           vr = xv(mgs,lr)
           qr = qx(mgs,lr)
           nrx = cx(mgs,lr)
           z = zx(mgs,lr)


! determine shape parameter alpha by iteration
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
           do i = 1,20
            if ( Abs(alp - alpha(mgs,lr)) .lt. 0.01 ) exit
             alpha(mgs,lr) = Max( rnumin, Min( rnumax, alp ) )
           alp = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/(z*pi**2) - 1.
             alp = Max( rnumin, Min( rnumax, alp ) )
           enddo

            
          endif
        endif

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
              g1 = 36.*(alpha(mgs,lr)+2.0)/((alpha(mgs,lr)+1.0)*pi**2)
           if ( .true. .and. (alpha(mgs,il) <= rnumin .or. alp == rnumin .or. alp == rnumax) ) then

            if ( rescale_high_alpha .and. alp >= rnumax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(1./(xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( rescale_low_alphar .and. alp <= rnumin ) then
             z  = 36.*(alpha(mgs,lr)+2.0)*nrx*vr**2/((alpha(mgs,lr)+1.0)*pi**2)
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = zx(mgs,il)
            endif
           endif
           

           
           endif
          endif
          
          endif
          
          enddo
!        CALL cld_cpu('Z-MOMENT-1r')  
           
           
           elseif ( il == lh .or. il == lhl .or. il == lf .or. (il == lr .and. imurain == 1 )) then ! } { Rain, GRAUPEL OR HAIL

        
        
        do mgs = 1,ngscnt

        if ( lf > 1 .and. il == lf ) then 
           lfsave(mgs,5) = an(igs(mgs),jy,kgs(mgs),ln(il))
           lfsave(mgs,6) = cx(mgs,il)
        endif
        
        if ( il == lhl .and. lnhlf > 1 ) then
          if ( cx(mgs,lhl) > cxmin ) then
            frac = chxf(mgs,lhl)/cx(mgs,lhl)
          else
            frac = 0.0
          endif
        endif

        if ( il == lh .and. lnhf > 1 ) then
          if ( cx(mgs,lh) > cxmin ) then
            frach = chxf(mgs,lh)/cx(mgs,lh)
          else
            frach = 0.0
          endif
        endif



         if ( iresetmoments == 1 .or. iresetmoments == il .or. iresetmoments == -1  ) then ! { .or. qx(mgs,il) <= qxmin(il) 
         if ( zx(mgs,il) <= zxmin ) then !  .and. qx(mgs,il) > 0.05e-3 
!!            write(91,*) 'zx=0; qx,cx = ',1000.*qx(mgs,il),cx(mgs,il)
           qx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
         elseif ( iresetmoments == -1 .and. qx(mgs,il) < qxmin(il) ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)

           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         
         elseif ( cx(mgs,il) <= cxmin .and. iresetmoments /= -1 ) then !  .and. qx(mgs,il) > 0.05e-3  
           qx(mgs,lv) = qx(mgs,lv) + qx(mgs,il)
           zx(mgs,il) = 0.0
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif
         else
            if ( zx(mgs,il) < 0.0 ) then !  .and. qx(mgs,il) > 0.05e-3 
               zx(mgs,il) = 0.0
             endif
         endif !}


         if (  zx(mgs,il) <= zxmin .and. cx(mgs,il) <= cxmin ) then
           zx(mgs,il) = 0.0
           cx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),lv) = an(igs(mgs),jgs,kgs(mgs),lv) + an(igs(mgs),jgs,kgs(mgs),il)
           qx(mgs,il) = 0.0
           an(igs(mgs),jgs,kgs(mgs),il) = qx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
           an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
         endif
        
        if ( qx(mgs,il) .gt. qxmin(il) ) then !{

        xv(mgs,il) = rho0(mgs)*qx(mgs,il)/(xdn(mgs,il)*Max(1.0e-9,cx(mgs,il)))
        xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)

        if ( xv(mgs,il) .lt. xvmn(il) ) then
          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
        endif

          if ( zx(mgs,il) > 0.0 .and. cx(mgs,il) <= 0.0 ) then !{
!  have mass and reflectivity but no concentration, so set concentration, using default alpha
            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            z   = zx(mgs,il)
            qr  = qx(mgs,il)
!            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
            cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)


           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) > 0.0 ) then
!  have mass and concentration but no reflectivity, so set reflectivity, using default alpha
!            g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
!     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
            chw = cx(mgs,il)
            qr  = qx(mgs,il)
!            zx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!            zx(mgs,il) = Min(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            g1 = (6.0 + alphamax)*(5.0 + alphamax)*(4.0 + alphamax)/ &
     &            ((3.0 + alphamax)*(2.0 + alphamax)*(1.0 + alphamax))
            zx(mgs,il) = Max(zxmin*1.1, g1*dn(igs(mgs),jy,kgs(mgs))**2*(6*qr)**2/(chw*(pi*xdn(mgs,il))**2) )
            an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

           elseif ( zx(mgs,il) <= zxmin .and. cx(mgs,il) <= 0.0 ) then
!   How did this happen?
         ! set values according to dBZ of -10, or Z = 0.1
!              0.1 = 1.e18*0.224*an(ix,jy,kz,lzh)*(hwdn/rwdn)**2

!               write(0,*) 'GS: moment problem! il,c,z,q = ',il,cx(mgs,il),zx(mgs,il),qx(mgs,il)
               
               zx(mgs,il) = 1.e-19/0.224*(xdn0(lr)/xdn0(il))**2
               an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)
               
               g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
               z   = zx(mgs,il)
               qr  = qx(mgs,il)
!               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z
               cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(6.*qr)**2/(z*(pi*xdn(mgs,il))**2)
               an(igs(mgs),jgs,kgs(mgs),ln(il)) = cx(mgs,il)
               
!               write(0,*) 'GS: moment problem! reset il,c,z,q = ',il,cx(mgs,il),zx(mgs,il),qx(mgs,il)
               
           else
          ! have all valid moments, so find shape parameter
          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

          if ( zx(mgs,il) .gt. 0. ) then !{
           
!            rdi = z*(pi/6.*1000.)**2*chw/((rho0(mgs)*qr)**2)
            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)

!           alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'kz, alp, alpha(mgs,il) = ',kz,alp,alpha(mgs,il),rdi,z,xv
           do i = 1,10
!            IF ( 100.*Abs(alp - alpha(mgs,il))/(Abs(alpha(mgs,il))+1.e-5) .lt. 1. ) EXIT
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
!             alp = 1.e18*(6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/
!     :            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
!           print*,'i,alp = ',i,alp
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo


! check for artificial breakup (graupel/hail larger than allowed max size)
        if (  xv(mgs,il) .gt. xvmx(il) ) then !{
          tmp = cx(mgs,il)


          xv(mgs,il) = Min( xvmx(il), Max( xvmn(il),xv(mgs,il) ) )
          xmas(mgs,il) = xv(mgs,il)*xdn(mgs,il)
          cx(mgs,il) = rho0(mgs)*qx(mgs,il)/(xmas(mgs,il))
          if ( tmp < cx(mgs,il) ) then ! breakup
            g1 = 36.*(6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il))*pi**2)
             zx(mgs,il) = zx(mgs,il) + g1*(rho0(mgs)/xdn(mgs,il))**2*( (qx(mgs,il)/tmp)**2 * (tmp-cx(mgs,il)) )
             an(igs(mgs),jgs,kgs(mgs),lz(il)) = zx(mgs,il)

          chw = cx(mgs,il)
          qr  = qx(mgs,il)
          z   = zx(mgs,il)

            rdi = z*(pi/6.*xdn(mgs,il))**2*chw/((rho0(mgs)*qr)**2)
            alp = (6.0+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
           do i = 1,10
             if ( Abs(alp - alpha(mgs,il)) .lt. 0.01 ) exit
             alpha(mgs,il) = Max( alphamin, Min( alphamax, alp ) )
             alp = (6.+alpha(mgs,il))*(5.0+alpha(mgs,il))*(4.0+alpha(mgs,il))/   &
     &            ((3.0+alpha(mgs,il))*(2.0+alpha(mgs,il))*rdi) - 1.0
             alp = Max( alphamin, Min( alphamax, alp ) )
           enddo

            
          endif
        endif !}

!
! Check whether the shape parameter is at or less than the minimum, and if it is, reset the 
! concentration or reflectivity to match (prevents reflectivity from being out of balance with Q and N)
!
             g1 = (6.0 + alpha(mgs,il))*(5.0 + alpha(mgs,il))*(4.0 + alpha(mgs,il))/ &
     &            ((3.0 + alpha(mgs,il))*(2.0 + alpha(mgs,il))*(1.0 + alpha(mgs,il)))
 
           if ( ( lrescalelow(il) .or. rescale_high_alpha ) .and.  &
     &          ( alpha(mgs,il) <= alphamin .or. alp == alphamin .or. alp == alphamax ) ) then !{

            if ( rescale_high_alpha .and. alp >= alphamax - 0.01  ) then  ! reset c at high alpha to prevent growth in Z
              cx(mgs,il) = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2
              an(igs(mgs),jy,kgs(mgs),ln(il)) = cx(mgs,il)
            
            elseif ( lrescalelow(il) .and. alp <= alphamin .and. .not. (il == lh .and. icvhl2h > 0 ) .and. &
                     .not. ( il == lr .and. .not. rescale_low_alphar ) ) then ! alpha = alphamin, so reset Z to prevent growth in C
             
             wtest = .false.
             if ( irescalerainopt == 0 ) then
               wtest = .false.
             elseif ( irescalerainopt == 1 ) then
               wtest = qx(mgs,lc) > qxmin(lc) 
             elseif ( irescalerainopt == 2 ) then
               wtest = qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             elseif ( irescalerainopt == 3 ) then
               wtest = temcg(mgs) > rescale_tempthresh .and. qx(mgs,lc) > qxmin(lc) .and. wvel(mgs) < rescale_wthresh
             endif
             
             if ( il == lr .and. ( wtest .or. .not. rescale_low_alphar ) ) then
             ! certain situations where rain number is adjusted instead of Z. Helps avoid rain being 'zapped' by autoconverted 
             ! drops (i.e., favor preserving Z when alpha tries to go negative)
             chw = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z*(6./(pi*xdn(mgs,il)))**2 ! g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/z1
             cx(mgs,il) = chw
             an(igs(mgs),jy,kgs(mgs),ln(il)) = chw
             else
             ! Usual resetting of reflectivity moment to force consisntency between Q, N, Z, and alpha when alpha = alphamin
             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
             z  = z1*(6./(pi*xdn(mgs,il)))**2
             zx(mgs,il) = z
             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
             endif

!             z1 = g1*dn(igs(mgs),jy,kgs(mgs))**2*(qr)*qr/chw
!             z  = z1*(6./(pi*xdn(mgs,il)))**2
!             zx(mgs,il) = z
!             an(igs(mgs),jy,kgs(mgs),lz(il)) = z
            endif

           endif !}
          
          
           endif !}
          
           
           endif ! !}
 
          
          
          endif !}

! #if !defined( COMMAS ) && !defined( 1 )
        if ( lzr > 1 ) then
          alpha2d(igs(mgs),kgs(mgs),1) = alpha(mgs,lr)
        endif
        if ( lzh > 1 ) then
          alpha2d(igs(mgs),kgs(mgs),2) = alpha(mgs,lh)
        endif
        if ( lzhl > 1 ) then
          alpha2d(igs(mgs),kgs(mgs),3) = alpha(mgs,lhl)
        endif
! #endif

        if ( il == lhl .and. lnhlf > 1 ) then
        ! update chxf in case cx has changed
          chxf(mgs,lhl) = frac*cx(mgs,lhl)
        endif
        if ( il == lh .and. lnhf > 1 ) then
        ! update chxf in case cx has changed
          chxf(mgs,lh) = frach*cx(mgs,lh)
        endif


!      IF ( lf > 0 .and. il == lf .and. kgs(mgs) <= 20 .and. ( cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ) > 200. .or. cx(mgs,lf) > 400. )) THEN
!        write(0,*) 'ix,jy, kz, cf = ',igs(mgs)+ixbeg,jy+jybeg,kgs(mgs), an(igs(mgs),jy,kgs(mgs),ln(lf)),lfsave(mgs,5),lfsave(mgs,6)
!        write(0,*) 'qold,qxold,zold,zxold = ',lfsave(mgs,1),lfsave(mgs,2),lfsave(mgs,3),lfsave(mgs,4)
!        write(0,*) 'cf_new,pcfwi,pcfwd = ',cx(mgs,lf),cx(mgs,lf) + dtp*( pcfwi(mgs) + pcfwd(mgs) ),pcfwi(mgs) + pcfwd(mgs)
!      
!      ENDIF
        
        enddo ! mgs

!         CALL cld_cpu('Z-DELABK')  
        

!         CALL cld_cpu('Z-DELABK')  
        
        
 
           
           endif ! } }

          endif ! }}
          endif ! }

          do mgs = 1,ngscnt

            if ( il == lh ) then
            if ( lnhf > 1 ) then ! number of graupel from frozen drops
              an(igs(mgs),jy,kgs(mgs),lnhf) = Max( chxf(mgs,lh), 0.0)
            endif
            endif

            if ( il == lhl ) then
            
            if ( lnhlf > 1 ) then ! number of hail from frozen drops
!              an(igs(mgs),jy,kgs(mgs),lnhlf) = Min( cx(mgs,lhl), Max( chxf(mgs,lhl), 0.0) )
              an(igs(mgs),jy,kgs(mgs),lnhlf) = Max( chxf(mgs,lhl), 0.0)
            endif
            endif
            an(igs(mgs),jy,kgs(mgs),ln(il)) = Max(cx(mgs,il), 0.0)
          enddo
        endif ! }
      enddo ! il }

      if ( lcin > 1 ) then
      do mgs = 1,ngscnt
        an(igs(mgs),jy,kgs(mgs),lcin) = Max(0.0, ccin(mgs))
      end do
      endif

      if ( ipconc .ge. 2 ) then
      do mgs = 1,ngscnt
        if ( lss > 1 ) then
          an(igs(mgs),jy,kgs(mgs),lss) = Max(0.0, ssmax(mgs) )
        endif

        if ( lccn > 1 ) then
          an(igs(mgs),jy,kgs(mgs),lccn) = Max(0.0, ccnc(mgs) )
        endif
      end do
      endif
      
      elseif ( ipconc .eq. 0 .and. lni .gt. 1 ) then
      
          do mgs = 1,ngscnt
            an(igs(mgs),jy,kgs(mgs),lni) = Max(cx(mgs,li), 0.0)
          enddo


      end if

      if ( ldovol ) then

       do il = li,lhab

        if ( lvol(il) .ge. 1 ) then

          do mgs = 1,ngscnt

           an(igs(mgs),jy,kgs(mgs),lvol(il)) = Max( 0.0, vx(mgs,il) )
          enddo
          
        endif
      
       enddo
      
      endif
!
!
!
!
!
      if (ndebug .gt. 0 ) write(0,*) 'gs 12'



      if (ndebug .gt. 0 ) write(0,*) 'gs 13'

 9998 continue

      if ( kz .gt. nz-1 .and. ix .ge. itile) then
        if ( ix .ge. itile ) then
         go to 1200 ! exit gather scatter
        else
         nzmpb = kz
        endif
      else
        nzmpb = kz
      end if

      if ( ix .ge. itile ) then
        nxmpb = 1
        nzmpb = kz+1
      else
       nxmpb = ix+1
      end if

 1000 continue
 1200 continue
!
!  end of gather scatter (for this jy slice)
!
!

      return
      end subroutine nssl_2mom_gs
!
!--------------------------------------------------------------------------
!



!
!--------------------------------------------------------------------------
!


end module module_mp_nssl_2mom
