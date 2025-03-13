  MODULE constants

  implicit none

  public

    real :: g,rd,rv,cp,cpinv,cv,cpv,cvv,rcp,                                     &
            cpdcv,cvdcp,rovcp,rddcp,rddcv,rddrv,cvdrd,cpdrd,eps,reps,repsm1,     &
            cpt,cvt,pnum,xlv,lathv,xls,xlf,lvdcp,condc,cpl,cpi,lv1,lv2,ls1,ls2,  &
            rhow,svp1,svp2,svp3,svpt0,c_e1,c_e2,c_s,rcs,earth_radius,govtwo

      !----------------

      real, parameter ::  pi     = 3.1415926535897932384626433
      real, parameter ::  th0r   = 300.0
      real, parameter ::  to     = 273.15
      real, parameter ::  p00    = 1.0e5
      real, parameter ::  rp00   = 1.0/p00
      real, parameter ::  piddeg = pi/180.0
      real, parameter ::  degdpi = 180.0/pi
      real, parameter ::  clwsat = 1.e-6
      real, parameter ::  cmemin = 1.e-4
      real, parameter ::  karman = 0.40
      real, parameter ::  rkarman = 1.0/karman
      real, parameter ::  omega  = 2.0*pi/86400.0

      !----------------
      ! Implicit vertical diffusion.  For vialpha:
      !      0.0 = explicit forward-in-time (unstable if K dt / (dz^2) > 0.5)
      !      0.5 = implicit centered-in-time (Crank-Nicholson)
      !      1.0 = implicit backward-in-time (implicit Euler)
      real, parameter ::  vialpha = 1.0

      ! DO NOT CHANGE VIBETA:
      real, parameter ::  vibeta = 1.0 - vialpha

      !----------------

      ! open bc phase velocity:
      ! 211210: set this in param.F
      ! 220227: now set this via namelist, with check in param.F
!!!      real, parameter ::  cstar = 30.0

      ! open bc max velocity:
      real, parameter ::  csmax = 350.0

      !----------------

      ! parameters for LES subgrid turbulence:
      real, parameter ::  c_m  = 0.10
      real, parameter ::  c_l  = 0.8165
      real, parameter ::  ri_c = 0.25


      !----------------
      ! stuff for weno / advection:
      real, parameter :: epsilon = 1.0e-18
      !----------------


      integer, parameter :: nrkmax = 3

      real, parameter :: grads_undef  =  -999999.9

  CONTAINS

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine set_constants( testcase )
      implicit none

      integer, intent(in) :: testcase


      !--------------------
      !  For test cases:

      IF( testcase.eq.4 )THEN

        !  nonprecipitating stratocumulus (Stevens et al, 2005, MWR)
        g      = 9.81
        rd     = 287.0
        cp     = 1015.0
        cpv    = 1870.0
        xlv    = 2.47e6
        xls    = 2834000.0
        earth_radius  =  6370000.0

      ELSEIF( testcase.eq.5 )THEN

        !  drizzling stratocumulus (Ackerman et al, 2009, MWR)
        g      = 9.81
        rd     = 287.0
        cp     = 1004.0
        cpv    = 1870.0
        xlv    = 2.50e6
        xls    = 2834000.0
        earth_radius  =  6370000.0

      ELSE

        !  RCEMIP
!!!        g      = 9.79764
!!!        rd     = 287.04
!!!        cp     = 1004.64
!!!        cpv    = 1846.0
!!!        xlv    = 2501000.0
!!!        xls    = 2834700.0
!!!        earth_radius  =  6371000.0

        !  cm1 defaults:
        g      = 9.81
        rd     = 287.04
        cp     = 1005.7
        cpv    = 1870.0
        xlv    = 2501000.0
        xls    = 2834000.0
        earth_radius  =  6370000.0

      ENDIF

      !--------------------

      govtwo = 0.5*g
      rv     = 461.5
      cpinv  = 1.0/cp
      cv     = cp-rd
      cvv    = cpv-rv
      rcp    = 1.0/cp
      cpdcv  = cp/cv
      cvdcp  = cv/cp
      rovcp  = rd/cp
      rddcp  = rd/cp
      rddcv  = rd/cv
      rddrv  = rd/rv
      cvdrd  = cv/rd
      cpdrd  = cp/rd
      eps    = rd/rv
      reps   = rv/rd
      repsm1 = rv/rd-1.0
      cpt    = (cpv/cp)-1.0
      cvt    = (cvv/cv)-1.0
      pnum   = (cp*rv)-(cpv*rd)
      lathv  = xlv
      lvdcp  = xlv/cp
      condc  = xlv*xlv/(rv*cp)
      cpl    = 4190.0
      cpi    = 2106.0
      lv1    = xlv+(cpl-cpv)*to
      lv2    = cpl-cpv
      ls1    = xls+(cpi-cpv)*to
      ls2    = cpi-cpv
      xlf    = xls-xlv
      rhow   = 1.0e3
      svp1   = 0.6112
      svp2   = 17.67
      svp3   = 29.65
      svpt0  = 273.15


      IF( testcase.eq.4 .or. testcase.eq.5 )THEN
        lv1 = xlv
        lv2 = 0.0
      ENDIF


      ! DO NOT CHANGE THESE:
      c_e1 = c_m * c_l * c_l * ( 1.0 / ri_c - 1.0 )
      c_e2 = max( 0.0 , c_m * pi * pi - c_e1 )
      c_s = ( c_m * c_m * c_m / ( c_e1 + c_e2 ) )**0.25   ! Smagorinsky constant
      rcs = 1.0/c_s


      end subroutine set_constants

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE constants
  
! ---------------------------------------------------------------------------  

    ! GHB: index_module is for NSSL microphysics
  
  MODULE index_module
  
  implicit none

! Mixing Ratio
  integer :: lv = 2
  integer :: lc = 3
  integer :: lr = 4
  integer :: li = 5
  integer :: ls = 6
  integer :: lgh = 7
  integer :: lhl = 0
  integer :: lhw = 0
  integer :: lsw = 0
  integer :: lhlw = 0
  integer :: lhab = 7

! Number concentration

  integer :: lccna = 0 ! activated ccn
  integer :: lcina = 0 ! activated ice nuclei
  integer :: lccn = 8
  integer :: lccnuf = 0

  integer :: lnc = 1
  integer :: lnr = 1
  integer :: lni = 1
  integer :: lns = 1
  integer :: lnh = 1
  integer :: lnhl = 1

! Particle volume

  integer :: lvh = 0
  integer :: lvhl = 0

! reflectivity (6th moment)

  integer :: lzr = 0
  integer :: lzh = 0
  integer :: lzhl = 0
  
! flags
  
  logical :: f_ccna = .false.
  logical :: f_cina = .false.
  
  
  END MODULE index_module
 
