  MODULE radtrns3d_module

  use radlib3d_module

  implicit none

  private
  public :: radtrns,setradwrk,zenangl,fito3,julday

!-----------------------------------------------------------------------

  integer, public, parameter :: radstgr  = 0  ! this must be 0 in CM1
  integer, public, parameter :: radshade = 0  ! this must be 0 in CM1
  integer, public, parameter :: rlwopt   = 1  ! 1 = "high"  (recommended)
                                              ! 0 = "low"

!-----------------------------------------------------------------------
!  for zenangl:

  real, parameter :: pi  = 3.14159265358979323
  real, public ::     hour0, pi2, deg2rad, yrday, rad2deg
  integer, public ::  nxmid, nymid, source

!-----------------------------------------------------------------------
!
!##################################################################
!##################################################################
!######                                                      ######
!######                      RADCST.INC                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'radcst.inc' for ARPS
!
!  This file initializes values of some constants used for radiation
!  calculation that remain unchanged throughout a model run
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/18/1996
!
!  MODIFICATION HISTORY:
!
!  08/14/1996 (Yuhe Liu)
!    Combined radims.inc and radcst.inc into radcst.inc. Defined the
!    radiation working arrays in the ARPS main program
!
!  04/09/1997 (Yuhe Liu)
!    Changed the indeces of working arrays from parameters to
!    variables and listed them in common blocks. This is to cooperate
!    with nested runs which have different dimension sizes for
!    different grids.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Define buffer array for option radopt=2. For option radopt=1, the
!  buffer size should be set to 1 to reduce the memory allocation.
!  Otherwise, the dimension sizes should be the same as nx, ny, and
!  nz, and the buffer size should be larger than the total size of 27
!  2-d arrays and 44 3-d arrays.
!
!    n2d_radiat  number of 2-d arrays in the buffer
!    n3d_radiat  number of 3-d arrays in the buffer
!
!  The 2-d arrays should be always at the beginning of radbuf and
!  the 3-d then follow.
!
!-----------------------------------------------------------------------
!
  INTEGER, public :: n2d_radiat  ! number of 2-d arrays in the buffer
  INTEGER, public :: n3d_radiat  ! number of 3-d arrays in the buffer

  PARAMETER ( n2d_radiat = 27 )
  PARAMETER ( n3d_radiat = 44 )
!
!-----------------------------------------------------------------------
!
!  Define the indeces that specify the positions of 2-D arrays in a
!  buffer array which carrys some of the radiation variables
!  calculated and/or used for both radiation options, radopt=1, or 2.
!
!-----------------------------------------------------------------------
!
  INTEGER, public :: nrad2d   ! Number of 2-D arrays in radiation buffer array
  PARAMETER ( nrad2d = 10 )
           ! Buffur array to carry the variables calculated or used in
! radiation calculation. The last index defines the variables:
! 1  = nrsirbm,  Solar IR surface albedo for beam
! 2  = nrsirdf,  Solar IR surface albedo for diffuse
! 3  = nrsuvbm,  Solar UV surface albedo for beam
! 4  = nrsuvdf,  Solar UV surface albedo for diffuse
! 5  = ncosz,    Cosine of zenith
! 6  = ncosss,   Cosine of angle between sun light and
!                terrain slope
! 7  = nfdirir,  all-sky direct downward IR flux
!                (0.7-10 micron) at the surface
! 8  = nfdifir,  all-sky diffuse downward IR flux
           !                at the surface
! 9  = nfdirpar, all-sky direct downward par flux
!                (0.4-0.7 micron) at the surface
! 10 = nfdifpar, all-sky diffuse downward par flux
!                at the surface

  INTEGER, public :: nrsirbm  ! Index for solar IR surface albedo for beam
  PARAMETER ( nrsirbm  = 1 )

  INTEGER, public :: nrsirdf  ! Index for solar IR surface albedo for diffuse
  PARAMETER ( nrsirdf  = 2 )

  INTEGER, public :: nrsuvbm  ! Index for solar UV surface albedo for beam
  PARAMETER ( nrsuvbm  = 3 )

  INTEGER, public :: nrsuvdf  ! Index for solar UV surface albedo for diffuse
  PARAMETER ( nrsuvdf  = 4 )

  INTEGER, public :: ncosz    ! Index for cosine of solar zenith angle
  PARAMETER ( ncosz = 5 )

  INTEGER, public :: ncosss   ! Index for cosine of angle between sun light and
                      ! terrain slope
  PARAMETER ( ncosss = 6 )

  INTEGER, public :: nfdirir  ! Index for all-sky direct downward IR flux
                      !           (0.7-10 micron) at the surface
  PARAMETER ( nfdirir = 7 )

  INTEGER, public :: nfdifir  ! Index for all-sky diffuse downward IR flux
                      !           at the surface
  PARAMETER ( nfdifir = 8 )

  INTEGER, public :: nfdirpar ! Index for all-sky direct downward and par flux
                      !           (0.4-0.7 micron) at the surface
  PARAMETER ( nfdirpar = 9 )

  INTEGER, public :: nfdifpar ! Index for all-sky diffuse downward and par flux
                      !           at the surface
  PARAMETER ( nfdifpar = 10 )
!
!-----------------------------------------------------------------------
!
!  Define co2 concentration in mppv
!
!-----------------------------------------------------------------------
!
  REAL :: co2

  PARAMETER ( co2 = 370.0E-6 )  ! in mppv
!
!-----------------------------------------------------------------------
!
!  Define parameters for solid water components
!
!-----------------------------------------------------------------------
!
  REAL, public :: roqr,tnw,roqs,tns,roqg,tng
!!!  COMMON /radmicro/ roqr,tnw,roqs,tns,roqg,tng

!  These are specific to Goddard LFO scheme
!  Moved to init_physics in CM1:  GHB, 101202
!  PARAMETER ( roqr = 1.00   )
!  PARAMETER ( tnw  = 0.08   )
!  PARAMETER ( roqs = 0.10   )
!  PARAMETER ( tns  = 0.03   )
!  PARAMETER ( roqg = 0.913  )
!  PARAMETER ( tng  = 0.0004 )
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine IRRAD and the subroutine IRRAD calls.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ir2d1,ir2d2,ir2d3,ir2d4,ir2d5,ir2d6     ! for 2-d arrays

!!!  COMMON /ir2dcmn/ ir2d1,ir2d2,ir2d3,ir2d4,ir2d5,ir2d6

  INTEGER :: ir3d1, ir3d2, ir3d3, ir3d4, ir3d5     ! for 3-d arrays
  INTEGER :: ir3d6, ir3d7, ir3d8, ir3d9, ir3d10    ! for 3-d arrays
  INTEGER :: ir3d11,ir3d12,ir3d13,ir3d14,ir3d15    ! for 3-d arrays
  INTEGER :: ir3d16

!!!  COMMON /ir3dcmn/ ir3d1, ir3d2, ir3d3, ir3d4, ir3d5,                   &
!!!                   ir3d6, ir3d7, ir3d8, ir3d9, ir3d10,                  &
!!!                   ir3d11,ir3d12,ir3d13,ir3d14,ir3d15,                  &
!!!                   ir3d16

  INTEGER :: ir4d1, ir4d2                          ! for 4-d arrays

!!!  COMMON /ir4dcmn/ ir4d1, ir4d2

  INTEGER :: ir5d1                                 ! for 5-d arrays

!!!  COMMON /ir5dcmn/ ir5d1
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine SOLIR, SOLUV, and CLDFLX.
!
!-----------------------------------------------------------------------
!
  INTEGER :: so2d1, so2d2, so2d3, so2d4, so2d5     ! for 2-d arrays
  INTEGER :: so2d6, so2d7, so2d8, so2d9, so2d10
  INTEGER :: so2d11,so2d12,so2d13,so2d14,so2d15
  INTEGER :: so2d16,so2d17,so2d18,so2d19,so2d20
  INTEGER :: so2d21,so2d22,so2d23

!!!  COMMON /so2dcmn/ so2d1, so2d2, so2d3, so2d4, so2d5,                   &
!!!                   so2d6, so2d7, so2d8, so2d9, so2d10,                  &
!!!                   so2d11,so2d12,so2d13,so2d14,so2d15,                  &
!!!                   so2d16,so2d17,so2d18,so2d19,so2d20,                  &
!!!                   so2d21,so2d22,so2d23

  INTEGER :: so3d1, so3d2, so3d3, so3d4, so3d5     ! for 3-d arrays
  INTEGER :: so3d6, so3d7, so3d8, so3d9, so3d10    ! for 3-d arrays
  INTEGER :: so3d11,so3d12,so3d13,so3d14           ! for 3-d arrays

!!!  COMMON /so3dcmn/ so3d1, so3d2, so3d3, so3d4, so3d5,                   &
!!!                   so3d6, so3d7, so3d8, so3d9, so3d10,                  &
!!!                   so3d11,so3d12,so3d13,so3d14

  INTEGER :: so4d1, so4d2, so4d3, so4d4, so4d5     ! for 4-d arrays

!!!  COMMON /so4dcmn/ so4d1, so4d2, so4d3, so4d4, so4d5

  INTEGER :: so5d1, so5d2, so5d3, so5d4, so5d5     ! for 5-d arrays

!!!  COMMON /so5dcmn/ so5d1, so5d2, so5d3, so5d4, so5d5

!-----------------------------------------------------------------------
!
!  Parameters for radiation
!
!-----------------------------------------------------------------------
!
  INTEGER :: nradstp         ! Time steps to update radiation
                             ! calculation

  REAL :: cldh2m             ! Height in meter that separates high cloud
                             ! and middle cloud
  REAL :: cldm2l             ! Height in meter that separates middle cloud
                             ! and low cloud

!!!  INTEGER :: ict             ! Vertical indices of height cldh2m
!!!  INTEGER :: icb             ! Vertical indices of height cldm2l

  PARAMETER ( cldh2m = 5500.0, cldm2l = 3000.0 )   ! in meter

  REAL :: solarc    ! Solar constant
  PARAMETER ( solarc = 1353.0 )

!augustin
!!!  COMMON  /arps500/ nradstp

!-----------------------------------------------------------------------

  CONTAINS

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RADTRNS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radtrns(nx,ny,nz,rbufsz,tipa,myid,dx,dy,                     &
           ib,ie,jb,je,kb,ke,xh,yh,prs0s,olr,dsr,                       &  ! MS add olr, dsr
           ptprt,pprt,qv,qc,qr,qi,qs,qh,cvm,                            &
           ptbar,pbar,ppi, o31, rhostr, tsfc, zp,                       &
           radsw,rnflx,radswnet,radlwin, cosss,                         &
           rsirbm,rsirdf,rsuvbm,rsuvdf, cosz, az,                       &
           fdirir,fdifir,fdirpar,fdifpar,                               &
           plinv,tinv,qvinv,o3a,ccld,                                   &
           flxir,flcir,flxuv,flcuv,dfdts,                               &
           tauir,taual, tauswi,tauswl,reffi,reffl,                      &
           radbuf, tem1, swfrc,lwfrc,doirrad,dosorad,                   &
           effc,effi,effs,effr,effg,effis,                              &
           cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,ptype,g,cp,eqtset)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to compute the atmospheric radiation forcing
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/11/1996
!
!  MODIFICATION HISTORY:
!
!  04/09/1997 (Yuhe Liu)
!  Removed dimension size check statement. For normal ARPS run, the
!  dimension size check is now done in the main driver, arps##. For
!  nested runs, no need to check the dimension size because the
!  working array is allocated automatically from existing space.
!
!  Added call of subroutine SETRADWRK to set the indeces of working
!  arrays in the radiation buffer. Those indeces are passed by
!  common blocks in radcst.inc.
!
!  10/11/1998 (Keith Brewster)
!  Added option for using RH in cloud optical depth calculation.
!
!  10/30/1998 (Keith Brewster)
!  Added calculation of aerosol density.
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  07/16/03  (J. Brotzge)
!  Fixed new radswnet bug (detected by T. Katopodes)
!
!-----------------------------------------------------------------------
!
  use irrad3d_module, only : irrad
  use sorad3d_module, only : sorad
  IMPLICIT NONE

  INTEGER, intent(in) :: nx,ny,nz
  INTEGER, intent(in) :: rbufsz
  integer, intent(in) :: tipa,myid
  real, intent(in) :: dx,dy
  integer, intent(in) :: ib,ie,jb,je,kb,ke
  real, dimension(ib:ie), intent(in) :: xh
  real, dimension(jb:je), intent(in) :: yh
  real, intent(in) :: prs0s
  logical, intent(in) :: doirrad,dosorad
!
!-----------------------------------------------------------------------
!
!  Define ARPS variables
!
!-----------------------------------------------------------------------
!
  REAL :: x(nx+1)
  REAL :: y(ny+1)
!!!  REAL :: z(nz)
!!!  REAL :: za(nz)
  REAL :: zp(nx,ny,nz)

  REAL :: ptprt (nx,ny,nz)
  REAL :: pprt  (nx,ny,nz)
  REAL :: qv    (nx,ny,nz)
  REAL :: qc    (nx,ny,nz)
  REAL :: qr    (nx,ny,nz)
  REAL :: qi    (nx,ny,nz)
  REAL :: qs    (nx,ny,nz)
  REAL :: qh    (nx,ny,nz)
  REAL :: cvm   (nx,ny,nz)

  REAL :: ptbar (nx,ny,nz)
  REAL :: pbar  (nx,ny,nz)
  REAL :: ppi   (nx,ny,nz)
  REAL :: o31   (nx,ny,nz)

  REAL :: rhostr(nx,ny,nz)

  REAL :: tsfc  (nx,ny)     ! Surface temperature (K)

  REAL, intent(inout) ::  swfrc(nx,ny,nz)  ! Radiation forcing (K/s)
  REAL, intent(inout) ::  lwfrc(nx,ny,nz)  ! Radiation forcing (K/s)

  REAL :: radsw (nx,ny)     ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)     ! Net upward radiation flux
  REAL :: radswnet(nx,ny)   ! Net solar radiation at surface
  REAL :: radlwin(nx,ny)    ! Incoming longwave radiation at surface
  REAL :: flxd(nx,ny)         
  ! MS addition
  REAL :: olr(nx,ny)        ! TOA outgoing longwave
  REAL :: dsr(nx,ny)        ! TOA outgoing shortwave

  REAL :: cosz  (nx,ny)     ! Cosine of zenith
  REAL :: cosss (nx,ny)     ! Cosine of angle between sun light and
                            ! surface terrain slope

  REAL :: az, zen           ! Solar angles

  REAL :: sh(nx,ny)         ! augustin add sh
!
!-----------------------------------------------------------------------
!
!  Define 2-D variables for radiation calculation.
!
!-----------------------------------------------------------------------
!
  REAL :: rsirbm(nx,ny)     ! Solar IR surface albedo for beam radiation
  REAL :: rsirdf(nx,ny)     ! Solar IR surface albedo for diffuse radiation
  REAL :: rsuvbm(nx,ny)     ! Solar UV surface albedo for beam radiation
  REAL :: rsuvdf(nx,ny)     ! Solar UV surface albedo for diffuse radiation

  REAL :: fdirir (nx,ny)    ! all-sky direct downward IR flux
                            ! (0.7-10 micron) at the surface
  REAL :: fdifir (nx,ny)    ! all-sky diffuse downward IR flux
                            ! at the surface
  REAL :: fdirpar(nx,ny)    ! all-sky direct downward par flux
                            ! (0.4-0.7 micron) at the surface
  REAL :: fdifpar(nx,ny)    ! all-sky diffuse downward par flux
                            ! at the surface
! MS add
  REAL :: fdiruv(nx,ny)    ! all-sky direct downward uv flux
                            ! (<0.4 micron) at the surface
  REAL :: fdifuv(nx,ny)    ! all-sky diffuse downward uv flux
                            ! at the surface

!!!  REAL :: st4(nx,ny)        ! Emission by the surface
!
!-----------------------------------------------------------------------
!
!  Arrays which have the vertical coordinate inversed, that
!  is, k=1 is for top while k=nz is at the surface.
!
!-----------------------------------------------------------------------
!
  REAL :: plinv (nx,ny,nz)  ! Pressure in mb at scalar points
                 ! GHB:  I think plinv is really p at w points
  REAL :: tinv  (nx,ny,nz)  ! Temperature
  REAL :: qvinv (nx,ny,nz)  ! Water vapor mixing ratio (g/g)

  REAL :: o3a   (nx,ny,nz)  ! Ozone (o3) mixing ratio (g/g)
  REAL :: ccld  (nx,ny,nz)  ! Cloud coverage (fraction)

  REAL :: flxir (nx,ny,nz)  ! all-sky net downward flux
  REAL :: flcir (nx,ny,nz)  ! clear-sky net downward flux

  REAL :: flxuv (nx,ny,nz)  ! all-sky solar flux (downward minus upward)
  REAL :: flcuv (nx,ny,nz)  ! clear-sky solar flux (downward minus upward)

  REAL :: dfdts (nx,ny,nz)  ! Sensitivity of net downward flux to surface
                            ! temperature

  REAL :: tauir (nx,ny,nz)  ! Cloud optical depth for LW IR
  REAL :: taual (nx,ny,nz)  ! Aerosol optical thickness

  REAL :: tauswi(nx,ny,nz)  ! Cloud optical depth for solar IR for
                            ! ice particles
  REAL :: tauswl(nx,ny,nz)  ! Cloud optical depth for solar IR for
                            ! liquid particles

  REAL :: reffi (nx,ny,nz)  ! Effective cloud-particle size for
                            ! ice particles
  REAL :: reffl (nx,ny,nz)  ! Effective cloud-particle size for
                            ! liquid particles

  REAL :: tem1 (nx,ny,nz)   ! Work array for message passing

! hm 8/20/12
  real, intent(in), dimension(nz) :: effc,effi,effs,effr,effg,effis
  real, intent(in) :: cgs1,cgs2,cgs3
  real, intent(in) :: cgt1,cgt2,cgt3
  integer, intent(in) :: ptype,eqtset
  real, intent(in) :: g,cp

!    don't bother saving arrays in current CM1 implementation:  GHB, 100729
!  REAL :: temp (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: tempp (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: tempt (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temptp (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temppi (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temqv (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temqc (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temqr (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temqi (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temqs (nx,ny,nz)   ! Work array for new coordinate
!  REAL :: temqh (nx,ny,nz)   ! Work array for new coordinate
!!!  REAL :: dz (nz)            ! Vertical grid spacing
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!!!  INCLUDE 'globcst.inc'
!!!  INCLUDE 'phycst.inc'
!!!  INCLUDE 'radcst.incl'
!!!  INCLUDE 'radmore.incl'

!!!  INCLUDE 'bndry.inc'                !  add by Yunheng
!!!  INCLUDE 'mp.inc'                   
!
!-----------------------------------------------------------------------
!
!  Include file radcst.inc which contains the definition of dimension
!  sizes for 2-d and 3-d temporary arrays represented by a buffer,
!  radbuf. When radopt is NOT set to 1, the dimensions and buffer
!  sizes can be 1. Otherwise, the dimension sizes should be the same
!  as nx, ny, and nz, and the buffer size should be larger than the
!  total size of 27 2-d arrays and 44 3-d arrays.
!
!  integer n2d_rad  ! number of 2-d arrays in the buffer
!  integer n3d_rad  ! number of 3-d arrays in the buffer
!
!  integer rbufsz   ! nx*ny*(n2d_rad+n3d_rad*nz)
!  real radbuf( rbufsz )
!
!  The 2-d arrays should be always at the beginning of radbuf and
!  the 3-d arrays then following.
!
!-----------------------------------------------------------------------
!
  REAL :: radbuf( rbufsz )
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxy, nxyz
  INTEGER :: i,j,k
  INTEGER :: im,jm,km, ij, ijm

  INTEGER :: istgr, jstgr
  INTEGER :: nxodd, nyodd, jeven, jodd

  INTEGER :: m, n

  INTEGER :: night           ! Flag for night time

  REAL :: tqe, dp, coolrate,heatrate
  REAL :: pk,psfc,tk,qvsat,rh,hgtagl,ccld1,ccld2,aden,psqc,aersum
  REAL :: gcp,grav

!  REAL, ALLOCATABLE :: difir(:,:)
!  REAL, ALLOCATABLE :: difpa(:,:)
!  REAL, ALLOCATABLE :: ptmp(:,:,:)
!  REAL, ALLOCATABLE :: temc(:,:,:)
!  REAL, ALLOCATABLE :: temal(:,:,:)
!  REAL, ALLOCATABLE :: temsi(:,:,:)
!  REAL, ALLOCATABLE :: temsl(:,:,:)
!  REAL, ALLOCATABLE :: temri(:,:,:)
!  REAL, ALLOCATABLE :: temrl(:,:,:)
!  REAL, ALLOCATABLE :: h_rate(:,:,:)
!  REAL, ALLOCATABLE :: c_rate(:,:,:)
  real, dimension(nx,ny,nz) :: h_rate,c_rate

  LOGICAL :: high

  INTEGER :: rh2cldopt
  REAL :: rhcldwgt,qcwgt
  PARAMETER (rh2cldopt=0,rhcldwgt=0.667)

  INTEGER :: ict             ! Vertical indices of height cldh2m
  INTEGER :: icb             ! Vertical indices of height cldm2l

!-----------------------------------------------------------------------

!  REAL :: fclr(1,1)
!  REAL :: dbs(1,1)
!  REAL :: trant(1,1)
!  REAL :: th2o(1,1,6)
!  REAL :: tcon(1,1,3)
!  REAL :: tco2(1,1,6,2)

!  REAL :: pa(1,1,nz-1)
!  REAL :: dt(1,1,nz-1)
!  REAL :: sh2o(1,1,nz)
!  REAL :: swpre(1,1,nz)
!  REAL :: swtem(1,1,nz)
!  REAL :: sco3(1,1,nz)
!  REAL :: scopre(1,1,nz)
!  REAL :: scotem(1,1,nz)
!  REAL :: dh2o(1,1,nz-1)
!  REAL :: dcont(1,1,nz-1)
!  REAL :: dco2(1,1,nz-1)
!  REAL :: do3(1,1,nz-1)
!  REAL :: flxu(1,1,nz)
!  REAL :: flxdd(1,1,nz)
!  REAL :: clr(1,1,0:nz)
!  REAL :: blayer(1,1,0:nz)

!  REAL :: h2oexp(1,1,nz-1,6)
!  REAL :: conexp(1,1,nz-1,3)

!  REAL :: co2exp(1,1,nz-1,6,2)

!-----------------------------------------------------------------------

 
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
!!!  REAL :: aeroden

!-----------------------------------------------------------------------
!
!  Define additional layer to the top of model atmosphere
!
!-----------------------------------------------------------------------
!
  REAL :: padd                ! Additional pressure level (mb) to TOA
  REAL :: tadd                ! Temperature for additional layer to TOA

  PARAMETER (padd =  1.0)
  PARAMETER (tadd =223.0)

 

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!  ALLOCATE(difir(nx,ny))
!  ALLOCATE(difpa(nx,ny))
!  ALLOCATE(ptmp(nx,ny,nz))
!  ALLOCATE(temc(nx,ny,nz))
!  ALLOCATE(temal(nx,ny,nz))
!  ALLOCATE(temsi(nx,ny,nz))
!  ALLOCATE(temsl(nx,ny,nz))
!  ALLOCATE(temri(nx,ny,nz))
!  ALLOCATE(temrl(nx,ny,nz))
!  ALLOCATE(h_rate(nx,ny,nz))
!  ALLOCATE(c_rate(nx,ny,nz))

  IF ( rlwopt == 0 ) THEN
    high = .false.
  ELSE
    high = .true.
  END IF

  sh = 1.0

  gcp = 0.01*g/cp
  grav = 0.01*g

!!!  CALL setradwrk(nx,ny,nz)

  nxy  = nx*ny
  nxyz = nxy*nz
  qcwgt=1.-rhcldwgt

  IF ( (pbar(1,1,nz-1)*0.01) <= padd ) THEN
!    WRITE (6,'(a/a)')                                                   &
!        'The pressure at the top of atmosphere was too low to add ',    &
!        'additional levels. Check the sounding profile.',               &
!        'Program stopped in RADTRNS.'
!    CALL arpsstop('arpsstop stopped RADTRNS pressure to low at the top ',1)
    print *,'Program stopped in RADTRNS.'
    call stopcm1
  END IF

!!! 700   CONTINUE   
  ! replace go-to with do loop:   GHB, 100729

!!!      zen = acos(cosz(max(1,nx/2),max(1,ny/2)))
      zen = acos(cosz(1,1))
!---------------


  DO j=1,ny
    DO i=1,nx
      DO k=1,nz-1
        o3a(i,j,k)=qc(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  In the staggering case, only those points of (i,j) = (even,even)
!  and (odd,odd) will be calculated. Do defragment to move the points
!  to the left half domain, i=1,nx/2 and j=1,ny-1
!
!-----------------------------------------------------------------------
!
  nyodd = MOD(ny,2)
  nxodd = MOD(nx,2)

  ict = 0
  icb = 0

  m = nx
  n = ny

  IF ( radstgr == 0 ) THEN
!    WARNING: this loop is repeated for radstgr=1, changes made here must
!    also be made there as well.
      DO j=1,n
        DO i=1,m
          psfc = prs0s + pprt(i,j,1)
    DO km=3,nz-1
      k = nz+1-km
          pk = pbar(i,j,k)+pprt(i,j,k)
          ! plinv is pressure (in mb) at w points:
          plinv(i,j,km) = 0.005*((pbar(i,j,k  ) + pprt(i,j,k  ))        &
                                +(pbar(i,j,k+1) + pprt(i,j,k+1)))
          tk = (ptbar(i,j,k)+ptprt(i,j,k))*ppi(i,j,k)
          tinv (i,j,km) = MAX(tk, 190.)
          qvinv(i,j,km) = MAX(qv(i,j,k), 1.0E-6)

          tqe = qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qh(i,j,k)
          tqe = MAX(1.0E-20, tqe)

          if( ict.eq.0 )then
            if( (zp(i,j,k)-zp(i,j,2)).le.cldh2m )  ict = km
          endif
          if( icb.eq.0 )then
            if( (zp(i,j,k)-zp(i,j,2)).le.cldm2l )  icb = km
          endif

          IF( rh2cldopt > 0 ) THEN
            stop 22222
          ELSE
            ccld(i,j,km) = MIN(1.0, MAX(0.0, 0.25*ALOG10(tqe)+1.5))
!  borrow o3a to store psuedo-cloud
            o3a(i,j,k) = qc(i,j,k)
          END IF

!!!          psfc = 0.5* (pbar(i,j,1)+pprt(i,j,1) +                        &
!!!                       pbar(i,j,2)+pprt(i,j,2) )
!!!          psfc = prs0s + pprt(i,j,1)
          aden=2.0*aeroden(pk,psfc)
          taual(i,j,km)=aden*0.5*((pbar(i,j,k-1)+pprt(i,j,k-1))         &
                                 -(pbar(i,j,k+1)+pprt(i,j,k+1)))
        END DO
      END DO
    END DO

    if( ict.eq.0 .or. icb.eq.0 )then
      print *
      print *,'  Could not find ict,icb '
      print *
      print *,'  ict,icb = ',ict,icb
      print *
      call stopcm1
    endif

!    WARNING: this loop is repeated for radstgr=1, changes made here must
!    also be made there as well.
    DO j=1,n
      DO i=1,m
        ij = nx*(j-1) + i
        radbuf(ij) = tsfc(i,j)

        plinv(i,j,1) = padd
        tinv (i,j,1) = tadd
        qvinv(i,j,1) = 1.0E-6
        ccld (i,j,1) = 0.0
        taual(i,j,1) = 0.0

        plinv(i,j,3) = 0.01*( cgt1*(pbar(i,j,nz-2) + pprt(i,j,nz-2)) &
                             +cgt2*(pbar(i,j,nz-3) + pprt(i,j,nz-3)) &
                             +cgt3*(pbar(i,j,nz-4) + pprt(i,j,nz-4)) )
        plinv(i,j,2) = 0.5*(plinv(i,j,1)+plinv(i,j,3))
        tinv (i,j,2) = MAX(0.5*(tinv(i,j,1)+tinv(i,j,3)), 190.0)
        qvinv(i,j,2) = MAX(qv(i,j,nz-1), 1.0E-6)
        ccld (i,j,2) = 0.0
        taual(i,j,2) = 0.0

!!!        plinv(i,j,nz) = 0.005*((pbar(i,j,1) + pprt(i,j,1))              &
!!!                            +(pbar(i,j,2) + pprt(i,j,2)))
!!!        plinv(i,j,nz) = 0.01*( prs0s + pprt(i,j,1) )
        plinv(i,j,nz) = 0.01*( cgs1*(pbar(i,j,2) + pprt(i,j,2)) &
                              +cgs2*(pbar(i,j,3) + pprt(i,j,3)) &
                              +cgs3*(pbar(i,j,4) + pprt(i,j,4)) )
      END DO
    END DO

  ELSE IF ( radstgr == 1 ) THEN
    stop 77771

  END IF

  DO km=1,nz-1
    DO j=1,n
      DO i=1,m
        dp = plinv(i,j,km+1) - plinv(i,j,km)
        IF ( dp <= 0.0 ) THEN
        print*, 'dp is=',dp
        print*, 'location is=', i  ,j,  km
        print*, 'plinv =',plinv(i,j,:)
!          WRITE (6,'(a,i3,a,i3,a,i3,a,i3,a/a,a/a)')                     &
!              'ERROR: The pressure gradient between level k = ',nz+1-km, &
!              ' and ',nz-km, ' at i = ',i,' and j = ',j,' was <=0.',    &
!              'Please check the sounding file, ',                       &
!              runname(1:lfnkey)//'.sound, or the data sets.',           &
!              'Program stopped in RADTRNS.'
!          CALL arpsstop('arpsstop stopped RADTRNS problem with sounding',1)
          print *,'Program stopped in RADTRNS.'
          call stopcm1
        END IF
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Notes: The arguments in subroutine CLDOPTD have different vertical
!         coordinates orders:
!
!         plinv  -- 1 for top
!         tinv   -- 1 for top
!
!         qc     -- 1 for bottom
!         qr     -- 1 for bottom
!         qi     -- 1 for bottom
!         qs     -- 1 for bottom
!         qh     -- 1 for bottom
!
!         tauir  -- 1 for top
!         tauswi -- 1 for top
!         tauswl -- 1 for top
!         reffi  -- 1 for top
!         reffl  -- 1 for top
!
!-----------------------------------------------------------------------
!
!  note borrowed o3a to store psuedo-cloud in place of qc.
  CALL cldoptd(nx,ny,m,n,nz, radstgr,                                   &
               plinv,tinv,o3a,qr,qi,qs,qh, rhostr, zp,                  &
               tauir,tauswi,tauswl,reffi,reffl,zen,                     &
               effc,effi,effs,effr,effg,effis,ptype)   ! hm 8/20/12

  IF ( radstgr == 1 ) THEN
    stop 77771
  END IF
!
!-----------------------------------------------------------------------
!
!  Fit the ozone concentration by interpolating a standard o3 profile
!  to ARPS pressure levels.
!
!-----------------------------------------------------------------------
!
!!!  CALL fito3(nx,ny,m,n, nz-1, plinv,o3a)
    do k=1,nz-1
    do j=1,ny
    do i=1,nx
      o3a(i,j,k) = o31(i,j,k)
    enddo
    enddo
    enddo
!
!-----------------------------------------------------------------------
!
!  Calculate the downward longwave IR radiation.
!
!  Positions of 2-d arrays in the buffer:
!
!    fclr (m,n)     -- radbuf(1+ 1*nxy)
!    dbs  (m,n)     -- radbuf(1+ 2*nxy)
!    trant(m,n)     -- radbuf(1+ 3*nxy)
!
!    th2o (m,n,6)   -- radbuf(1+ 4*nxy)
!    tcon (m,n,3)   -- radbuf(1+10*nxy)
!    tco2 (m,n,6,2) -- radbuf(1+13*nxy)
!
!  Positions of 3-d arrays in the buffer:
!
!    pa    (m,n,np)     -- radbuf(1+25*nxy)
!    dt    (m,n,np)     -- radbuf(1+25*nxy+ 1*nxyz)
!    sh2o  (m,n,np+1)   -- radbuf(1+25*nxy+ 2*nxyz)
!    swpre (m,n,np+1)   -- radbuf(1+25*nxy+ 3*nxyz)
!    swtem (m,n,np+1)   -- radbuf(1+25*nxy+ 4*nxyz)
!    sco3  (m,n,np+1)   -- radbuf(1+25*nxy+ 5*nxyz)
!    scopre(m,n,np+1)   -- radbuf(1+25*nxy+ 6*nxyz)
!    scotem(m,n,np+1)   -- radbuf(1+25*nxy+ 7*nxyz)
!    dh2o  (m,n,np)     -- radbuf(1+25*nxy+ 8*nxyz)
!    dcont (m,n,np)     -- radbuf(1+25*nxy+ 9*nxyz)
!    dco2  (m,n,np)     -- radbuf(1+25*nxy+10*nxyz)
!    do3   (m,n,np)     -- radbuf(1+25*nxy+11*nxyz)
!    flxu  (m,n,np+1)   -- radbuf(1+25*nxy+12*nxyz)
!    flxd  (m,n,np+1)   -- radbuf(1+25*nxy+13*nxyz)
!    clr   (m,n,0:np+1) -- radbuf(1+25*nxy+14*nxyz) ! nz+1
!    blayer(m,n,0:np+1) -- radbuf(1+26*nxy+15*nxyz) ! nz+1
!
!    h2oexp(m,n,np,6)   -- radbuf(1+27*nxy+16*nxyz)
!    conexp(m,n,np,3)   -- radbuf(1+27*nxy+22*nxyz)
!
!    co2exp(m,n,np,6,2) -- radbuf(1+27*nxy+25*nxyz)
!
!-----------------------------------------------------------------------
!
 IF (doirrad) THEN
  ! (arps code is absurd ... interpolate, then intepolate back?  ugh.)
  ! give irrad the pressure at scalar pts directly:
  do j=1,ny
  do i=1,nx
    km = 1
    radbuf(ir3d1-1+km) = 0.5*(plinv(i,j,km)+plinv(i,j,km+1))
    km = 2
    radbuf(ir3d1-1+km) = 0.5*(plinv(i,j,km)+plinv(i,j,km+1))
    do km=3,nz-1
      k = nz+1-km
      radbuf(ir3d1-1+km) = 0.01*(pbar(i,j,k)+pprt(i,j,k))
    enddo
  enddo
  enddo
  m = nx
  n = ny
  CALL irrad(nx,ny,m,n,nz-1,                                            &
             tauir(1,1,1),ccld(1,1,1), plinv(1,1,1),tinv(1,1,1),qvinv(1,1,1),o3a(1,1,1), co2,radbuf(1),           &
             high,flxd(1,1),flxir(1,1,1),flcir(1,1,1),dfdts(1,1,1),rnflx(1,1),                        &
!!!           fclr,dbs,trant,th2o,tcon,tco2,                               &
!!!           pa,dt,sh2o,swpre,swtem,sco3,scopre,scotem,                   &
!!!           dh2o,dcont,dco2,do3,flxu,flxdd,clr,blayer,                    &
!!!           h2oexp,conexp,co2exp)
             radbuf(ir2d1),radbuf(ir2d2),radbuf(ir2d3),                 &
             radbuf(ir2d4),radbuf(ir2d5),radbuf(ir2d6),                 &
             radbuf(ir3d1),radbuf(ir3d2),radbuf(ir3d3),                 &
             radbuf(ir3d4),radbuf(ir3d5),radbuf(ir3d6),                 &
             radbuf(ir3d7),radbuf(ir3d8),radbuf(ir3d9),                 &
             radbuf(ir3d10),radbuf(ir3d11),radbuf(ir3d12),              &
             radbuf(ir3d13),radbuf(ir3d14),radbuf(ir3d15),              &
             radbuf(ir3d16),radbuf(ir4d1),radbuf(ir4d2),                &
             radbuf(ir5d1))
  IF( eqtset.eq.1 )THEN
    DO k=2,nz-2
      km=nz+1-k          ! inverse vertical coordinates to ARPS grid
      DO j=1,ny
        DO i=1,nx
          ! NOTE:  plinv is pressure (in mb) at w points:
!!!          c_rate(i,j,k) = 9.770687E-05               & ! = g/cp in cgs unit
!!!               * ( flxir(i,j,km+1) - flxir(i,j,km) )                    &
!!!                   / (ppi(i,j,k) * ( plinv(i,j,km) - plinv(i,j,km+1) ))
          ! cm1r17:  c_rate is now a TEMPERATURE tendency: 
          c_rate(i,j,k) = gcp                        & ! = g/cp in cgs unit
                   * ( flxir(i,j,km+1) - flxir(i,j,km) )                    &
                   / ( plinv(i,j,km) - plinv(i,j,km+1) )
        END DO
      END DO
    END DO
  ELSEIF( eqtset.eq.2 )THEN
    DO k=2,nz-2
      km=nz+1-k          ! inverse vertical coordinates to ARPS grid
      DO j=1,ny
        DO i=1,nx
          ! NOTE:  plinv is pressure (in mb) at w points:
!!!          c_rate(i,j,k) = 9.770687E-05               & ! = g/cp in cgs unit
!!!               * ( flxir(i,j,km+1) - flxir(i,j,km) )                    &
!!!                   / (ppi(i,j,k) * ( plinv(i,j,km) - plinv(i,j,km+1) ))
          ! cm1r17:  c_rate is now a TEMPERATURE tendency: 
          c_rate(i,j,k) = grav                       & ! = g/cp in cgs unit
                   * ( flxir(i,j,km+1) - flxir(i,j,km) )                    &
                   / ( cvm(i,j,k)*( plinv(i,j,km) - plinv(i,j,km+1) ) )
        END DO
      END DO
    END DO
  ENDIF

 END IF
!
!-----------------------------------------------------------------------
!
!  Calculate solar radiation fluxes.
!
!  Output flxuv and flcuv are the fractions to incoming solar flux
!  at the top of atmosphere
!
!  Positions of 2-d arrays used in subroutine SORAD
!
!    sdf  (m,n)     -- radbuf(1+ 1*nxy    )
!    sclr (m,n)     -- radbuf(1+ 2*nxy)
!    csm  (m,n)     -- radbuf(1+ 3*nxy)
!    cc   (m,n,3)   -- radbuf(1+ 4*nxy)
!
!  Positions of 3-d arrays used in subroutine SORAD
!
!    tauclb(m,n,np)   -- radbuf(1+ 7*nxy            )
!    tauclf(m,n,np)   -- radbuf(1+ 7*nxy+ 1*nxyz)
!    dp    (m,n,np)   -- radbuf(1+ 7*nxy+ 2*nxyz)
!    wh    (m,n,np)   -- radbuf(1+ 7*nxy+ 3*nxyz)
!    oh    (m,n,np)   -- radbuf(1+ 7*nxy+ 4*nxyz)
!    scal  (m,n,np)   -- radbuf(1+ 7*nxy+ 5*nxyz)
!    swh   (m,n,np+1) -- radbuf(1+ 7*nxy+ 6*nxyz)
!    so2   (m,n,np+1) -- radbuf(1+ 7*nxy+ 7*nxyz)
!    df    (m,n,np+1) -- radbuf(1+ 7*nxy+ 8*nxyz)
!
!  Positions of temporary arrays used in subroutine SOLIR, SOLUV, and
!  CLDFLX:
!
!    tem2d1 (m,n)         -- radbuf(1+ 7*nxy+ 9*nxyz)
!    tem2d2 (m,n)         -- radbuf(1+ 8*nxy+ 9*nxyz)
!    tem2d3 (m,n)         -- radbuf(1+ 9*nxy+ 9*nxyz)
!    tem2d4 (m,n)         -- radbuf(1+10*nxy+ 9*nxyz)
!    tem2d5 (m,n)         -- radbuf(1+11*nxy+ 9*nxyz)
!    tem2d6 (m,n)         -- radbuf(1+12*nxy+ 9*nxyz)
!    tem2d7 (m,n)         -- radbuf(1+13*nxy+ 9*nxyz)
!    tem2d8 (m,n)         -- radbuf(1+14*nxy+ 9*nxyz)
!    tem2d9 (m,n)         -- radbuf(1+15*nxy+ 9*nxyz)
!    tem2d10(m,n)         -- radbuf(1+16*nxy+ 9*nxyz)
!    tem2d11(m,n)         -- radbuf(1+17*nxy+ 9*nxyz)
!    tem2d12(m,n)         -- radbuf(1+18*nxy+ 9*nxyz)
!    tem2d13(m,n)         -- radbuf(1+19*nxy+ 9*nxyz)
!    tem2d14(m,n)         -- radbuf(1+20*nxy+ 9*nxyz)
!    tem2d15(m,n)         -- radbuf(1+21*nxy+ 9*nxyz)
!    tem2d16(m,n)         -- radbuf(1+22*nxy+ 9*nxyz)
!    tem2d17(m,n)         -- radbuf(1+23*nxy+ 9*nxyz)
!    tem2d18(m,n)         -- radbuf(1+24*nxy+ 9*nxyz)
!    tem2d19(m,n)         -- radbuf(1+25*nxy+ 9*nxyz)
!
!    tem3d1 (m,n,np+1)    -- radbuf(1+26*nxy+ 9*nxyz)
!    tem3d2 (m,n,np+1)    -- radbuf(1+26*nxy+10*nxyz)
!    tem3d3 (m,n,np+1)    -- radbuf(1+26*nxy+11*nxyz)
!    tem3d4 (m,n,np+1)    -- radbuf(1+26*nxy+12*nxyz)
!    tem3d5 (m,n,np+1)    -- radbuf(1+26*nxy+13*nxyz)
!
!    tem4d1 (m,n,np+1,2)  -- radbuf(1+26*nxy+14*nxyz)
!    tem4d2 (m,n,np+1,2)  -- radbuf(1+26*nxy+16*nxyz)
!    tem4d3 (m,n,np+1,2)  -- radbuf(1+26*nxy+18*nxyz)
!    tem4d4 (m,n,np+1,2)  -- radbuf(1+26*nxy+20*nxyz)
!    tem4d5 (m,n,np+1,2)  -- radbuf(1+26*nxy+22*nxyz)
!
!    tem5d1(m,n,np+1,2,2) -- radbuf(1+26*nxy+24*nxyz)
!    tem5d2(m,n,np+1,2,2) -- radbuf(1+26*nxy+28*nxyz)
!    tem5d3(m,n,np+1,2,2) -- radbuf(1+26*nxy+32*nxyz)
!    tem5d4(m,n,np+1,2,2) -- radbuf(1+26*nxy+36*nxyz)
!    tem5d5(m,n,np+1,2,2) -- radbuf(1+26*nxy+40*nxyz)
!
!-----------------------------------------------------------------------
!
 ! GHB, 100729:
 IF( dosorad )THEN

  night = 1
  DO j=1,n
    DO i=1,m
      IF ( cosz(i,j) > 0.0 ) THEN
        night = 0
        GO TO 500
      END IF
    END DO
  END DO

  500   CONTINUE

!  aersum=0.
!  DO 505 k=1,nz-1
!    aersum=aersum+taual(2,2,k)
! 505 CONTINUE
!  write(6,'(a,f10.4))') ' Total aerosol optical depth: ',aersum

  IF ( night == 0 ) THEN

    m = nx
    n = ny

    CALL sorad(nx,ny,m,n,nz-1, plinv,tinv,qvinv,o3a,co2,                &
         tauswi,tauswl,reffi,reffl,ccld,ict,icb,                        &
         taual,rsirbm,rsirdf,rsuvbm,rsuvdf,cosz,                        &
         flxuv,flcuv,fdirir,fdifir,fdirpar,fdifpar,                     &
         fdiruv,fdifuv,                                                 & ! MS add
         radbuf(so2d1),radbuf(so2d2),radbuf(so2d3),radbuf(so2d4),       &
         radbuf(so3d1),radbuf(so3d2),radbuf(so3d3),radbuf(so3d4),       &
         radbuf(so3d5),radbuf(so3d6),radbuf(so3d7),radbuf(so3d8),       &
         radbuf(so3d9),radbuf(so2d5),radbuf(so2d6),radbuf(so2d7),       &
         radbuf(so2d8),radbuf(so2d9),radbuf(so2d10),radbuf(so2d11),     &
         radbuf(so2d12),radbuf(so2d13),radbuf(so2d14),                  &
         radbuf(so2d15),radbuf(so2d16),radbuf(so2d17),                  &
         radbuf(so2d18),radbuf(so2d19),radbuf(so2d20),                  &
         radbuf(so2d21),radbuf(so2d22),radbuf(so2d23),                  &
         radbuf(so3d10),radbuf(so3d11),radbuf(so3d12),                  &
         radbuf(so3d13),radbuf(so3d14),radbuf(so4d1),                   &
         radbuf(so4d2),radbuf(so4d3),radbuf(so4d4),radbuf(so4d5),       &
         radbuf(so5d1), radbuf(so5d2),                                  &
         radbuf(so5d3), radbuf(so5d4), radbuf(so5d5))

  ELSE

    DO k=1,nz
      DO j=1,n
        DO i=1,m
          flxuv(i,j,k) = 0.0
        END DO
      END DO
    END DO

    DO j=1,n
      DO i=1,m
        fdirir (i,j) = 0.0
        fdifir (i,j) = 0.0
        fdirpar(i,j) = 0.0
        fdifpar(i,j) = 0.0
        fdiruv(i,j) = 0.0 ! MS add
        fdifuv(i,j) = 0.0 ! MS add
      END DO
    END DO

  ENDIF

 ENDIF

!------------------------------------------------------------------------------


!   IF ((tipa == 1) .AND. (repeat == 0))THEN !.AND. (abs(cosz(1,1)) > 0.5))THEN
!    don't bother saving arrays in current CM1 implementation:  GHB, 100729
!      ptmp = plinv
!      temc = ccld
!      temal = taual
!      temsi = tauswi
!      temsl = tauswl
!      temri = reffi
!      temrl = reffl
!      temp = pbar
!      tempp = pprt
!      tempt = ptbar
!      temptp = ptprt
!      temppi = ppi
!      temqv = qv
!      temqc = qc
!      temqr = qr
!      temqi = qi
!      temqs = qs
!      temqh = qh
!      tem1 = rhostr
!
!      repeat = 1
!
!      zen = acos(cosz(nx/2,ny/2))
!
!      GOTO 700
!
!   END IF


!   IF (tipa == 1) THEN

!    don't bother saving arrays in current CM1 implementation:  GHB, 100729
!      pbar = temp
!      pprt = tempp
!      ptbar = tempt
!      ptprt = temptp
!      qv = temqv
!      qc = temqc
!      qr = temqr
!      qi = temqi
!      qs = temqs
!      qh = temqh
!      rhostr = tem1
!      ccld = temc
!      taual = temal
!      tauswi = temsi
!      tauswl = temsl
!      reffi = temri
!      reffl = temrl

      !Fill-in edges of SW fluxes

!    DO k=1,nz
!        DO j=1,ny
!          flxuv(1,j,k) = flxuv(2,j,k)
!          flxuv(nx,j,k) = flxuv(nx-1,j,k)
!          flcuv(1,j,k) = flcuv(2,j,k)
!          flcuv(nx,j,k) = flcuv(nx-1,j,k)
!        END DO
!
!        DO i=1,nx
!          flxuv(i,1,k) = flxuv(i,2,k)
!          flxuv(i,ny,k) = flxuv(i,ny-1,k)
!          flcuv(i,1,k) = flcuv(i,2,k)
!          flcuv(i,ny,k) = flcuv(i,ny-1,k)
!        END DO
!
!        flxuv(1,1,k) = flxuv(2,2,k)
!        flxuv(nx,1,k) = flxuv(nx-1,2,k)
!        flxuv(nx,ny,k) = flxuv(nx-1,ny-1,k)
!        flxuv(1,ny,k) = flxuv(2,ny-1,k)
!        flcuv(1,1,k) = flcuv(2,2,k)
!        flcuv(nx,1,k) = flcuv(nx-1,2,k)
!        flcuv(nx,ny,k) = flcuv(nx-1,ny-1,k)
!        flcuv(1,ny,k) = flcuv(2,ny-1,k)
!
!    END DO
!   END IF


!  ELSEIF ((repeat == 0) .AND. (tipa == 1)) THEN 
!
!    ! not needed:  GHB, 100729
!    DO k=1,nz-1
!      DO j=1,n
!        DO i=1,m
!          flxuv(i,j,k) = 0.0
!        END DO
!      END DO
!    END DO
!
!    ! moved back to original location (just after sorad):  GHB, 100729
!    DO j=1,n
!      DO i=1,m
!        fdirir (i,j) = 0.0
!        fdifir (i,j) = 0.0
!        fdirpar(i,j) = 0.0
!        fdifpar(i,j) = 0.0
!      END DO
!    END DO

!  END IF
!
!-----------------------------------------------------------------------
!
!  Added the heating rate of solar radiation to total radiation
!  forcing (K/s)
!
!  Constant 9.770687e-05 is equal to g/cp in cgs unit, where g = 980
!  and cp = 1.003e7.
!
!  Outputs from SORAD such as flxuv, flcuv, etc., are fractions of
!  solar flux at the top of atmosphere. Therefore we need to multipy
!  solar constant and cosine of zenith angle to obtain the solar
!  radiation flux.
!
!-----------------------------------------------------------------------
!

!    don't bother saving arrays in current CM1 implementation:  GHB, 100729
!  IF (repeat == 0) THEN
!    temppi = ppi
!    ptmp = plinv
!  END IF

  IF ( radstgr == 0 ) THEN

  IF( night.eq.0 .and. dosorad )THEN

! moved calculation of coolrate up near irrad to avoid unnecessary
! calculations when using tipa:  GHB, 100726

  IF( eqtset.eq.1 )THEN
    DO k=2,nz-2
      km=nz+1-k          ! inverse vertical coordinates to ARPS grid
      DO j=1,ny
        DO i=1,nx
          ! NOTE:  plinv is pressure (in mb) at w points:
!!!          h_rate(i,j,k) = solarc * cosz(i,j) * 9.770687E-05                  &
!!!                   * ( flxuv(i,j,km+1) - flxuv(i,j,km) )                &
!!!                   / (ppi(i,j,k) * ( plinv(i,j,km) - plinv(i,j,km+1) ))
          ! cm1r17:  h_rate is now a TEMPERATURE tendency: 
          h_rate(i,j,k) = solarc * cosz(i,j) * gcp                           &
                   * ( flxuv(i,j,km+1) - flxuv(i,j,km) )                &
                   / ( plinv(i,j,km) - plinv(i,j,km+1) )
        END DO
      END DO
    END DO
  ELSEIF( eqtset.eq.2 )THEN
    DO k=2,nz-2
      km=nz+1-k          ! inverse vertical coordinates to ARPS grid
      DO j=1,ny
        DO i=1,nx
          ! NOTE:  plinv is pressure (in mb) at w points:
!!!          h_rate(i,j,k) = solarc * cosz(i,j) * 9.770687E-05                  &
!!!                   * ( flxuv(i,j,km+1) - flxuv(i,j,km) )                &
!!!                   / (ppi(i,j,k) * ( plinv(i,j,km) - plinv(i,j,km+1) ))
          ! cm1r17:  h_rate is now a TEMPERATURE tendency: 
          h_rate(i,j,k) = solarc * cosz(i,j) * grav                          &
                   * ( flxuv(i,j,km+1) - flxuv(i,j,km) )                &
                   / ( cvm(i,j,k)*( plinv(i,j,km) - plinv(i,j,km+1) ) )
        END DO
      END DO
    END DO
  ENDIF

  ELSE
    h_rate = 0.0
  ENDIF


!!!    radfrc = 0

  IF( doirrad )THEN
    DO k=2,nz-2
      DO j=1,ny
        DO i=1,nx
          lwfrc(i,j,k) = c_rate(i,j,k)
        END DO
      END DO
    END DO
    DO j=1,ny
      DO i=1,nx
        radlwin(i,j) = flxd(i,j) 
      END DO
    END DO
  ENDIF

  IF( dosorad )THEN
    DO k=2,nz-2
      DO j=1,ny
        DO i=1,nx
          swfrc(i,j,k) = h_rate(i,j,k)  
        END DO
      END DO
    END DO

    ! MS altered code begin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! MS altered this calculation to include UV (<0.4 micron) radiation
    ! don't know what 'sh' is, but in current code is set to unity.
    ! also changed to use cosz instead of cosss. In my view this is more correct in 
    ! the presence of terrain than cosss. 

    DO j=1,ny
      DO i=1,nx
        ! net radiation at the surface
        rnflx(i,j) = solarc * radsw(i,j) * cosz(i,j)                    & ! radsw = a2dr2
                 * ( sh(i,j)*(1.0-rsirbm(i,j)) * fdirir(i,j)            &
                       + sh(i,j)*(1.0-rsuvbm(i,j)) * fdirpar(i,j)       &
                       + sh(i,j)*(1.0-rsuvbm(i,j)) * fdiruv(i,j)        &
                       + (1.0-rsirdf(i,j)) * fdifir(i,j)                &
                       + (1.0-rsuvdf(i,j)) * fdifpar(i,j)               &
                       + (1.0-rsuvdf(i,j)) * fdifuv(i,j) )              &
                       + flxir(i,j,nz)       ! net downward LW flux at sfc
	! net shortwave radiation at the surface
        radswnet(i,j) = solarc * radsw(i,j) * cosz(i,j)       & ! radsw = a2dr2
                 * (sh(i,j)*(1.0-rsirbm(i,j)) * fdirir(i,j)             &
                 + sh(i,j)*(1.0-rsuvbm(i,j)) * fdirpar(i,j)             &
                 + sh(i,j)*(1.0-rsuvbm(i,j)) * fdiruv(i,j)              &
                          + (1.0-rsirdf(i,j)) * fdifir(i,j)             &
                          + (1.0-rsuvdf(i,j)) * fdifpar(i,j)		&
                          + (1.0-rsuvdf(i,j)) * fdifuv(i,j) )
        ! MS add
        ! TOA fluxes - currently calculated at third level from top
        !              This is to be consistent with the heating rates
        !              which ignore the top two flux levels 
        olr(i,j) = -flxir(i,j,3)
        dsr(i,j) = solarc*radsw(i,j)*cosz(i,j)*flxuv(i,j,3)
        radsw(i,j) = solarc * radsw(i,j) * cosz(i,j)                      &
                     * ( sh(i,j)*(fdirir(i,j) + fdirpar(i,j)+ fdiruv(i,j))  &
                       + fdifir(i,j) + fdifpar(i,j) + fdifuv(i,j) )

      END DO
    END DO
    ! MS end altered code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENDIF

  ELSE IF ( radstgr == 1 ) THEN
    stop 77773

  END IF

!---------------------------------------------------------------------
!
! Added by Yunheng to update the fake zone for radfrc, radsw, rnflx.
! two more variables (radswnet and radlwin) since IHOP_3
!
!--------------------------------------------------------------------

!  MPI stuff?  Commented out, for now:  GHB, 100720
!  IF (mp_opt > 0) THEN
!    CALL acct_interrupt(mp_acct)
!    CALL mpsendrecv2dew(radfrc, nx, ny, nz, ebc, wbc, 0, tem1)
!    CALL mpsendrecv2dns(radfrc, nx, ny, nz, nbc, sbc, 0, tem1)
!
!    CALL mpsendrecv1dew(radsw,  nx, ny, ebc, wbc, 0, tem1)
!    CALL mpsendrecv1dns(radsw,  nx, ny, nbc, sbc, 0, tem1)
!
!    CALL mpsendrecv1dew(rnflx,  nx, ny, ebc, wbc, 0, tem1)
!    CALL mpsendrecv1dns(rnflx,  nx, ny, nbc, sbc, 0, tem1)
!
!    CALL mpsendrecv1dew(radswnet, nx, ny, ebc, wbc, 0, tem1)
!    CALL mpsendrecv1dns(radswnet, nx, ny, nbc, sbc, 0, tem1)
!
!    CALL mpsendrecv1dew(radlwin,  nx, ny, ebc, wbc, 0, tem1)
!    CALL mpsendrecv1dns(radlwin,  nx, ny, nbc, sbc, 0, tem1)
!
!    CALL acct_stop_inter
!  END IF

!  IF ( raddiag == 1 ) THEN
!    WRITE(6,'(a,i8,a,f10.2,a)')                                         &
!        ' Dump radiation variables at time step,', nstep,               &
!        ', model time=',curtim,' (s)'
!!
!!-----------------------------------------------------------------------
!!
!!  Write out results to GrADS file for display
!!
!!-----------------------------------------------------------------------
!!
!    CALL wrtrad(nx,ny,nz,m,n,x,y,z,                                     &
!                plinv,tinv,qvinv,qc,qr,qi,qs,qh, o3a,radbuf(1),         &
!                ccld, tauir,taual,tauswi,tauswl,reffi,reffl,            &
!                rsirbm,rsirdf,rsuvbm,rsuvdf,                            &
!                fdirir,fdifir,fdirpar,fdifpar,                          &
!                dfdts, radsw,rnflx, cosz,                               &
!                flxir,flcir, flxuv,flcuv,                               &
!                radfrc)
!!    :              radfrc, coolr,heatr)
!
!  END IF

!  DEALLOCATE(difir)
!  DEALLOCATE(difpa)
!  DEALLOCATE(ptmp)
!  DEALLOCATE(temc)
!  DEALLOCATE(temal)
!  DEALLOCATE(temsi)
!  DEALLOCATE(temsl)
!  DEALLOCATE(temri)
!  DEALLOCATE(temrl)
!  DEALLOCATE(h_rate)
!  DEALLOCATE(c_rate)

  RETURN

END SUBROUTINE radtrns
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SETRADWRK                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setradwrk( nx,ny,nz )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the indeces for radiation working arrays
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/09/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! The number grid points in 3 directions
  INTEGER :: nxy, nxyz
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
!!!  INCLUDE 'radcst.incl'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxy  = nx*ny
  nxyz = nxy*nz
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine IRRAD and the subroutine IRRAD calls.
!
!-----------------------------------------------------------------------
!
  ir2d1  = 1+ 1*nxy
  ir2d2  = 1+ 2*nxy
  ir2d3  = 1+ 3*nxy
  ir2d4  = 1+ 4*nxy
  ir2d5  = 1+10*nxy
  ir2d6  = 1+13*nxy

  ir3d1  = 1+25*nxy+ 0*nxyz
  ir3d2  = 1+25*nxy+ 1*nxyz
  ir3d3  = 1+25*nxy+ 2*nxyz
  ir3d4  = 1+25*nxy+ 3*nxyz
  ir3d5  = 1+25*nxy+ 4*nxyz
  ir3d6  = 1+25*nxy+ 5*nxyz
  ir3d7  = 1+25*nxy+ 6*nxyz
  ir3d8  = 1+25*nxy+ 7*nxyz
  ir3d9  = 1+25*nxy+ 8*nxyz
  ir3d10 = 1+25*nxy+ 9*nxyz
  ir3d11 = 1+25*nxy+10*nxyz
  ir3d12 = 1+25*nxy+11*nxyz
  ir3d13 = 1+25*nxy+12*nxyz
  ir3d14 = 1+25*nxy+13*nxyz
  ir3d15 = 1+25*nxy+14*nxyz
  ir3d16 = 1+26*nxy+15*nxyz

  ir4d1  = 1+27*nxy+16*nxyz
  ir4d2  = 1+27*nxy+22*nxyz

  ir5d1  = 1+27*nxy+25*nxyz
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine SOLIR, SOLUV, and CLDFLX.
!
!-----------------------------------------------------------------------
!
  so2d1  = 1+ 1*nxy
  so2d2  = 1+ 2*nxy
  so2d3  = 1+ 3*nxy
  so2d4  = 1+ 4*nxy
  so2d5  = 1+ 7*nxy
  so2d6  = 1+ 8*nxy
  so2d7  = 1+ 9*nxy
  so2d8  = 1+10*nxy
  so2d9  = 1+11*nxy
  so2d10 = 1+12*nxy
  so2d11 = 1+13*nxy
  so2d12 = 1+14*nxy
  so2d13 = 1+15*nxy
  so2d14 = 1+16*nxy
  so2d15 = 1+17*nxy
  so2d16 = 1+18*nxy
  so2d17 = 1+19*nxy
  so2d18 = 1+20*nxy
  so2d19 = 1+21*nxy
  so2d20 = 1+22*nxy
  so2d21 = 1+23*nxy
  so2d22 = 1+24*nxy
  so2d23 = 1+25*nxy

  so3d1  = 1+26*nxy+ 0*nxyz
  so3d2  = 1+26*nxy+ 1*nxyz
  so3d3  = 1+26*nxy+ 2*nxyz
  so3d4  = 1+26*nxy+ 3*nxyz
  so3d5  = 1+26*nxy+ 4*nxyz
  so3d6  = 1+26*nxy+ 5*nxyz
  so3d7  = 1+26*nxy+ 6*nxyz
  so3d8  = 1+26*nxy+ 7*nxyz
  so3d9  = 1+26*nxy+ 8*nxyz
  so3d10 = 1+26*nxy+ 9*nxyz
  so3d11 = 1+26*nxy+10*nxyz
  so3d12 = 1+26*nxy+11*nxyz
  so3d13 = 1+26*nxy+12*nxyz
  so3d14 = 1+26*nxy+13*nxyz

  so4d1  = 1+26*nxy+14*nxyz
  so4d2  = 1+26*nxy+16*nxyz
  so4d3  = 1+26*nxy+18*nxyz
  so4d4  = 1+26*nxy+20*nxyz
  so4d5  = 1+26*nxy+22*nxyz

  so5d1  = 1+26*nxy+24*nxyz
  so5d2  = 1+26*nxy+28*nxyz
  so5d3  = 1+26*nxy+32*nxyz
  so5d4  = 1+26*nxy+36*nxyz
  so5d5  = 1+26*nxy+40*nxyz

  RETURN
END SUBROUTINE setradwrk
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTRAD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  REMOVED ... GHB 100720
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CLDOPTD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cldoptd(nx,ny,m,n,nz, radstgr,                               &
           pres,temp,qc,qr,qi,qs,qh, rhostr, zp,                        &
           tauir,tauswi,tauswl,reffi,reffl,zen,                         &
           effc,effi,effs,effr,effg,effis,ptype)  ! hm 8/20/12
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the optical depth
!
!-----------------------------------------------------------------------
!
!  AUTHOR: NASA Goddard Center
!
!  MODIFICATION:
!
!  03/11/1996 (Yuhe Liu)
!  Modified the original code from 1-D to 3-D
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!                                               Vertical index order
!    pres     Pressure (mb)                       1 for top
!    temp     Temperature (K),                    1 for top
!    qc       Cloud water mixing ratio (g/g),     1 for bottom
!    qr       Rain water mixing ratio (g/g),      1 for bottom
!    qi       cloud ice mixing ratio (g/g),       1 for bottom
!    qs       Snow mixing ratio (g/g),            1 for bottom
!    qh       Hail mixing ratio (g/g),            1 for bottom
!    rhostr   Density multiply by j3,             1 for bottom
!    j3inv    1/j3,                               1 for bottom
!    zen      Solar zenith angle
!
!  OUTPUT:
!
!    tauir    Cloud optical depth for longwave,   1 for bottom
!    tauswi   Cloud optical depth for ice cloud
!             for shortwave,                      1 for bottom
!    tauswl   Cloud optical depth for liquid cloud
!             for shortwave,                      1 for bottom
!    reffi    Effective cloud-particle size (mm)
!             for ice cloud for shortwave,        1 for bottom
!    reffl    Effective cloud-particle size (mm)
!             for liquid cloud for shortwave,     1 for bottom
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: m,n
  INTEGER :: radstgr

  REAL :: pres  (nx,ny,nz)      ! in mb
  REAL :: temp  (nx,ny,nz)      ! in K
  REAL :: qc    (nx,ny,nz)      ! in g/g
  REAL :: qr    (nx,ny,nz)      ! in g/g
  REAL :: qi    (nx,ny,nz)      ! in g/g
  REAL :: qs    (nx,ny,nz)      ! in g/g
  REAL :: qh    (nx,ny,nz)      ! in g/g
  REAL :: rhostr(nx,ny,nz)      ! in kg/m**3
  REAL :: zp    (nx,ny,nz)      ! in kg/m**3

  REAL :: tauir (nx,ny,nz)
  REAL :: tauswi(nx,ny,nz)
  REAL :: tauswl(nx,ny,nz)
  REAL :: reffi (nx,ny,nz)
  REAL :: reffl (nx,ny,nz)

  REAL :: tauqc
  REAL :: tauqr
  REAL :: tauqi
  REAL :: tauqs
  REAL :: tauqh
  REAL :: reff1
  REAL :: reff2
  REAL :: zen

  REAL :: w1, effrad

  REAL :: cpi,twco, dpg, rho

! hm 8/20/12
  real, intent(in), dimension(nz) :: effc,effi,effs,effr,effg,effis
  integer, intent(in) :: ptype

  INTEGER :: i,j,k, im,jm,km
!!!  INTEGER :: istgr, jstgr, jodd

      real, parameter :: g      = 9.81
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!!!  INCLUDE 'phycst.inc'
!!!  INCLUDE 'radcst.incl'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  cpi  = 4.*ATAN(1.)
  twco = 1.e-6

!!!  IF ( radstgr == 0 ) THEN
!!!    istgr = 0
!!!    jstgr = 1
!!!  ELSE
!!!    stop 77774
!!!  END IF

  if( nx.ne.1 .or. ny.ne.1 ) stop 23323

  DO km=3,nz-1
    k = nz+1-km

    DO j=1,ny
      DO i=1,nx

        im = i
        jm = j

        rho = rhostr(i,j,k)*1.0E-3        ! g/cm**3
!!!        dpg = 10.0*(pres(im,jm,km+1)-pres(im,jm,km))/g ! g/cm**2
        ! cm1r17:  instead of using (dry?) hydrostatic equation, 
        !          use actual rho and actual delta-z
        dpg = 0.1*rhostr(i,j,k)*(zp(i,j,k+1)-zp(i,j,k))

        IF ( qc(i,j,k) >= twco ) THEN
          w1     = dpg*qc(i,j,k)
          effrad = 0.0015
! hm 8/20/12
          if(ptype.eq.5) effrad = effc(k)*1.e-4 !convert from micron to cm
          tauqc  = w1/effrad
          reff2  = effrad*1.0E4
        ELSE
          tauqc = 0.0
          reff2 = 0.0
        END IF

        IF ( qr(i,j,k) >= twco ) THEN
          w1     = dpg*qr(i,j,k)
          effrad = 3./((cpi*tnw*roqr/(rho*qr(i,j,k)))**.25)
! hm 8/20/12
          if(ptype.eq.5) effrad = effr(k)*1.e-4 !convert from micron to cm 
          tauqr  = w1/effrad
        ELSE
          tauqr=0.0
        END IF

        IF ( qi(i,j,k)+qs(i,j,k) >= twco ) THEN
          w1 = 1.e4*dpg*(qi(i,j,k)+qs(i,j,k))

          IF ( temp(im,jm,km) > 243.16 ) THEN
            effrad = 0.0125
          ELSE IF ( temp(im,jm,km) < 223.16 ) THEN
            effrad = 0.0025
          ELSE
            effrad = 0.0125+(temp(im,jm,km)-243.16)*0.00050
          END IF

! hm 8/20/12 overwrite with predicted effective radii
          if(ptype.eq.5) effrad = effs(k)*1.e-4
          tauqs = w1*(-0.006656 +  3.686E-4/effrad)
          if(ptype.eq.5) effrad = effi(k)*1.e-4
          tauqi = w1*(-0.011500 +  4.110E-4/effrad                      &
                                + 17.300E-8/(effrad*effrad))
! hm 8/20/12 combined cloud ice + snow effective radius
          if(ptype.eq.5) effrad = effis(k)*1.e-4 ! micron to cm
          reff1 = effrad*1.0E4
        ELSE
          tauqi = 0.0
          tauqs = 0.0
          reff1 = 0.0
        END IF

        IF ( qh(i,j,k) >= twco ) THEN
          w1     = dpg*qh(i,j,k)
          effrad = 3./((cpi*tng*roqg/(rho*qh(i,j,k)))**.25)
! hm 8/20/12
          if(ptype.eq.5) effrad = effg(k)*1.e-4 ! micron to cm
          tauqh  = w1/effrad
        ELSE
          tauqh  = 0.0
        END IF

!!!        tauswi(im,jm,km) = tauqs + tauqh
        ! GHB, bug fix (I think):  include qi  (130903)
        tauswi(im,jm,km) = tauqs + tauqh + tauqi
        tauswl(im,jm,km) = 1.5 * ( tauqc + tauqr )
        reffi (im,jm,km) = reff1
        reffl (im,jm,km) = reff2
        tauir (im,jm,km)  = 0.5 * tauswl(im,jm,km) + tauqi + tauqh
      END DO
    END DO
  END DO

  IF ( radstgr /= 0 .AND. MOD(ny,2) == 0 ) THEN
    stop 77775
  END IF

  DO jm=1,n
    DO im=1,m
      tauswi(im,jm,1) = 0.0
      tauswl(im,jm,1) = 0.0
      reffi (im,jm,1) = 0.0
      reffl (im,jm,1) = 0.0
      tauir (im,jm,1) = 0.0

      tauswi(im,jm,2) = 0.0
      tauswl(im,jm,2) = 0.0
      reffi (im,jm,2) = 0.0
      reffl (im,jm,2) = 0.0
      tauir (im,jm,2) = 0.0
    END DO
  END DO

  RETURN
END SUBROUTINE cldoptd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FITO3                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fito3(nx,ny,m,n,np,pl,ao,prs0,o30,ib,ie,jb,je,kb,ke,nk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to fit o3 to the model grid
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Modified the subroutine from 1-D to 3-D
!
!-----------------------------------------------------------------------
!
!fpp$ expand (terp1)
!!dir$ inline always terp1
!*$*  inline routine (terp1)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,np
  INTEGER :: m,n

  REAL :: pl(nx,ny,np+1)      ! Model pressure (mb)
  REAL :: ao(nx,ny,np+1)        ! Model o3 mixing ratio (g/g)

  integer, intent(in) :: ib,ie,jb,je,kb,ke,nk
  real, dimension(ib:ie,jb:je,kb:ke) :: prs0,o30
!
!-----------------------------------------------------------------------
!
!  Local definitions
!
!-----------------------------------------------------------------------
!
  INTEGER :: lay
  PARAMETER (lay=75)

!  integer iop(2),itab(3),iflag(5)
  INTEGER :: iop(2),itab(3)

  REAL :: pa(lay)       ! Local layer pressure (mb)
  REAL :: ta(lay)       ! Local layer temperature (K)
  REAL :: wa(lay)       ! Local layer water vapor mixing ratio (g/g)
  REAL :: oa(lay)       ! Local layer o3 mixing ratio (g/g)

!  real tsi(5)
!  real w(lay)
  REAL :: tab(3)
  REAL :: wk(lay,4)

  INTEGER :: i,j,k
  INTEGER :: INT, ix
!  integer lun

!  real y
  REAL :: p
!
!----- ix=        1:trp; 2:mls; 3:mlw; 4:sas; 5:saw
!  data iflag/  1   ,  0   ,  0   ,  0   ,  0   /
!  data   tsi/ 300.0, 294.0, 272.2, 287.0, 257.1/
!
!-----------------------------------------------------------------------
!
!  Local data bank for pressure, temperature, moisture, and o3
!
!-----------------------------------------------------------------------
!
  DATA (pa(k),k=1,lay)/                                                 &
           .0003,    .0008,    .0011,    .0015,     .0021,              &
           .0029,    .0041,    .0058,    .0081,     .0113,              &
           .0158,    .0221,    .0310,    .0435,     .0609,              &
           .0855,    .1200,    .1700,    .2400,     .3350,              &
           .4650,    .6500,    .9150,   1.2850,    1.8000,              &
          2.5250,   3.5450,   4.9700,   6.9700,    9.7800,              &
         13.7150,  19.2350,  26.9850,  37.8550,   53.1000,              &
         73.8900,  97.6650, 121.4350, 145.2100,  168.9900,              &
        192.7650, 216.5400, 240.3150, 264.0900,  287.8650,              &
        311.6350, 335.4100, 359.1900, 382.9650,  406.7400,              &
        430.5150, 454.2850, 478.0600, 501.8350,  525.6100,              &
        549.3900, 573.1650, 596.9400, 620.7150,  644.4900,              &
        668.2650, 692.0350, 715.8100, 739.5850,  763.3600,              &
        787.1400, 810.9150, 834.6900, 858.4650,  882.2400,              &
        906.0150, 929.7850, 953.5600, 977.3350, 1001.1100/

  DATA (ta(k),k=1,lay)/                                                 &
        209.86, 210.20, 210.73, 211.27, 211.81,                         &
        212.35, 212.89, 213.44, 213.98, 214.53,                         &
        215.08, 215.62, 216.17, 216.74, 218.11,                         &
        223.20, 230.04, 237.14, 244.46, 252.00,                         &
        259.76, 267.70, 274.93, 274.60, 269.38,                         &
        262.94, 256.45, 250.12, 244.31, 238.96,                         &
        233.74, 228.69, 224.59, 221.75, 219.10,                         &
        216.64, 215.76, 215.75, 215.78, 216.22,                         &
        219.15, 223.79, 228.29, 232.45, 236.33,                         &
        239.92, 243.32, 246.53, 249.56, 252.43,                         &
        255.14, 257.69, 260.11, 262.39, 264.57,                         &
        266.66, 268.67, 270.60, 272.48, 274.29,                         &
        276.05, 277.75, 279.41, 281.02, 282.59,                         &
        284.09, 285.53, 286.86, 288.06, 289.13,                         &
        290.11, 291.03, 291.91, 292.76, 293.59/


  DATA (wa(k),k=1,lay)/                                                 &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.406E-05, 0.520E-05,          &
        0.115E-04, 0.275E-04, 0.572E-04, 0.107E-03, 0.166E-03,          &
        0.223E-03, 0.285E-03, 0.360E-03, 0.446E-03, 0.547E-03,          &
        0.655E-03, 0.767E-03, 0.890E-03, 0.103E-02, 0.118E-02,          &
        0.136E-02, 0.159E-02, 0.190E-02, 0.225E-02, 0.264E-02,          &
        0.306E-02, 0.351E-02, 0.399E-02, 0.450E-02, 0.504E-02,          &
        0.560E-02, 0.619E-02, 0.680E-02, 0.742E-02, 0.805E-02,          &
        0.869E-02, 0.935E-02, 0.100E-01, 0.107E-01, 0.113E-01/

  DATA (oa(k),k=1,lay)/                                                 &
       .643E-07, .202E-06, .246E-06, .290E-06, .334E-06,                &
       .378E-06, .422E-06, .467E-06, .512E-06, .557E-06,                &
       .603E-06, .648E-06, .694E-06, .740E-06, .793E-06,                &
       .101E-05, .131E-05, .164E-05, .198E-05, .234E-05,                &
       .272E-05, .312E-05, .359E-05, .465E-05, .590E-05,                &
       .765E-05, .910E-05, .960E-05, .994E-05, .101E-04,                &
       .990E-05, .853E-05, .710E-05, .576E-05, .423E-05,                &
       .260E-05, .152E-05, .102E-05, .786E-06, .598E-06,                &
       .448E-06, .352E-06, .302E-06, .252E-06, .212E-06,                &
       .193E-06, .176E-06, .160E-06, .147E-06, .137E-06,                &
       .127E-06, .118E-06, .109E-06, .103E-06, .975E-07,                &
       .924E-07, .883E-07, .846E-07, .810E-07, .778E-07,                &
       .749E-07, .721E-07, .694E-07, .671E-07, .648E-07,                &
       .626E-07, .607E-07, .593E-07, .579E-07, .565E-07,                &
       .552E-07, .540E-07, .528E-07, .517E-07, .505E-07/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ix=1

  iop(1)=4
  iop(2)=4
  INT=1

  CALL coeff(lay,pa,oa,wa,iop,INT,wk)

!$omp parallel do default(shared)  &
!$omp private(i,j,k,p,tab,itab)
  DO k=kb,ke
    DO j=jb,je
      DO i=ib,ie
        itab(1)=1
        itab(2)=0
        itab(3)=0
!!!        p = 0.5 * ( pl(i,j,k+1) + pl(i,j,k) )
        p = 0.01*prs0(i,j,k)
        CALL terp1(lay,pa,oa,wa,p,INT,tab,itab)
!!!        ao(i,j,k)=tab(1)
        o30(i,j,nk+1-k)=tab(1)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE fito3
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION AERODEN                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  REAL FUNCTION aeroden(p,psfc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign aerosol optical density (per Pa)
!  For now we use an estimate based on the aerosol density used
!  in simple climate models, where 0.05 of optical depth is
!  assumed uniformly distributed over mass sfc to 800 mb and
!  0.025 is uniformly distributed over mass 800 mb to 225 mb.
!
!  Here we normalize that distribution according to the
!  surface pressure so that the mountains don't lose aerosols.
!
!  Aerosol optical depth is computed by calling routine as
!  layer depth (Pa) times this aerosol density.
!
!  AUTHOR:   Keith Brewster
!
!  MODIFICATION HISTORY
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: p,psfc
!!!  REAL :: aeroden
!
  IF (p < (.225*psfc)) THEN
    aeroden=0.
  ELSE IF (p < (.80*psfc)) THEN
    aeroden=0.025/(0.575*psfc)
  ELSE
    aeroden=0.050/(0.2*psfc)
  END IF
  RETURN
  END FUNCTION aeroden

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE JULDAY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE julday( year, month, day, jday )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute Julian day from year, month, and day
!
!  Start from 1 (Jan. 1) to 365, or 366 for leap year (Dec. 31)
!
!  The rule is that a year will be a leap year if
!
!    the year can be divided by 400, or
!    the year can by divided by 4, but not by 100
!
!  Form this rule year 1972, 1984, 1996, and 2000 are leap years,
!  but 1700, 1800 and 1900 are not.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  07/29/93
!
!  MODIFICATIONS:
!
!  05/06/1998 (Yuhe Liu)
!  Corrected the leap year calculation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    year       Reference calendar year
!    month      Reference monthe of the year
!    day        Reference day of the month
!
!    OUTPUT:
!
!    jday       Julian day, start from 1 -- Jan. 1 to 365 -- Dec. 31
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: year, month, day, jday
  INTEGER :: lpyear, lp

  INTEGER :: mndys(12)     ! Day numbers for each month
  DATA mndys/0,31,59,90,120,151,181,212,243,273,304,334/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( MOD(year,400) == 0 .OR.                                          &
         (MOD(year,4) == 0 .AND. MOD(year,100) /= 0 ) ) THEN
    lpyear = 1
  ELSE
    lpyear = 0
  END IF

  lp = 0
  IF ( month > 2 ) lp = lpyear

  jday = mndys(month) + day + lp

  RETURN
END SUBROUTINE julday

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ZENANGL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE zenangl(nx,ny,      hterain, cosz, cosss, a2dr2,             &
           rjday,tloc, latscl,lonscl, slpmag,slpdir,                    &
           tem1,tem2,saltitude,sazimuth,dx,dy,dtbig,curtim,             &
           ctrlat,ctrlon,year,month,day,hour,minute,second,jday)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate cosine of solar zenith angle.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  11/16/93
!
!  MODIFICATION HISTORY:
!  2/2/99  Vince Wong and Jik Leong
!  This modification calculates the solar declination angle and
!  equation of time using a method found on page C24 of the
!  1996 Astronomical Almanac.
!  The mothod is good to 0.01 degrees in the sky over the
!  period 1950 to 2050.
!
! augustin
!  8/23/01 Augustin Colette EFML/ Stanford University
!
! Computation of the solar altitude and azimuth
! saltitude and sazimuth are outputs of zenangl to be used in shade
! sources:
! http://www.usc.edu/dept/architecture/mbs/tools/vrsolar/Help/ &
! solar_concepts.html
! http://www.uwinnipeg.ca/~blair/physclim/lab2.htm
! http://ra.stsci.edu/cgi-bin/gethelp.cgi?altaz.src
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    x        X coordinates at scalar points
!    y        Y coordinates at scalar points
!    hterain  Surface terrain
!
!  OUTPUT:
!
!    cosz     Cosine of zenith
!    cosss    Cosine of angle between sun light and terrain slope
!    a2dr2    Square ratio of average distance to the time
!             dependent distance from the earth to the sun
!
!augustin
!    sazimuth  solar azimuth
!    saltitude solar altitude
!
!  WORK ARRAY:
!
!    rjday    Julian day at each grid point
!    tloc     Local time at each grid point
!    latscl   Latitudes  at scalar points
!    lonscl   Longitudes at scalar points
!    slpmag   Surface terrain slope magnitude
!    slpdir   Surface terrain slope direction
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny
  real :: dx,dy
  real :: dtbig,curtim
  real :: ctrlat,ctrlon
  integer :: year,month,day,hour,minute,second,jday

  REAL :: hterain(nx,ny)

  REAL :: cosz(nx,ny),   cosss(nx,ny), a2dr2(nx,ny)
  REAL :: rjday(nx,ny),  tloc(nx,ny)
  REAL :: latscl(nx,ny), lonscl(nx,ny)
  REAL :: slpmag(nx,ny), slpdir(nx,ny)
!augustin add saltitude and sazimuth in the outputs of zenanlg
  REAL :: saltitude
  REAL :: sazimuth

  REAL :: tem1(nx,ny), tem2(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Include file:
!
!-----------------------------------------------------------------------
!
!!!  INCLUDE 'radcst.incl'
!!!  INCLUDE 'radmore.incl'
!!!  INCLUDE 'globcst.inc'
!!!  INCLUDE 'grid.inc'          ! Grid & map parameters.
!!!  INCLUDE 'phycst.inc'
!!!  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j

  REAL :: xs, ys

!!!  REAL :: hour0, yrday
!!!  REAL :: deg2rad, pi, pi2

  REAL :: etau, shrangl, sdeclin
  REAL :: azimuth, sinz
  REAL :: dpsi, sinpsi, cospsi

  REAL :: anncyc

  REAL :: hr, days2k, lsun, gsun, obliq, lambda, xsun, ysun
!!!  REAL :: asun, alpha, rad2deg
  REAL :: asun, alpha

  LOGICAL :: firstcall        ! First call flag of this subroutine

!!!  include 'radzen.incl'

!!!  SAVE firstcall, hour0, pi, pi2, deg2rad, yrday, rad2deg
!!!  DATA firstcall/.true./

!!!  ! added to calculate saltitude and sazimuth at the middle domain
!!!  INTEGER :: nxmid, nymid, source 
!!!  SAVE nxmid, nymid, source

  REAL :: shrangl_mid
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!----- moved this section of code to initrad -----!
!  IF (firstcall) THEN
!    pi  = 3.14159265358979
!    pi2 = 2.0 * pi
!    deg2rad = pi/180.0
!    rad2deg = 1./deg2rad
!
!    hour0 = FLOAT(hour)                                                 &
!          + FLOAT(minute)/60.0                                          &
!          + FLOAT(second)/3600.0
!
!    IF ( MOD(year, 4) == 0 ) THEN
!      yrday = 366.
!    ELSE
!      yrday = 365.
!    END IF
!
!    nxmid = CEILING( 0.5*((nx-3)*nproc_x + 3) )   ! Middle point index at
!    nymid = CEILING( 0.5*((ny-3)*nproc_y + 3) )   ! global domain
!!!! Not using arps MPI code:  GHB, 100720
!    source = 0
!!!!    source = proc( (nxmid-2)/(nx-3)+1 + ( (nymid-2)/(ny-3) )*nproc_x )
!!!!                                ! source processor contain the middle point.
!    nxmid = MOD( (nxmid-2), (nx-3) ) + 2    ! local index of central domain
!    nymid = MOD( (nymid-2), (ny-3) ) + 2
!
!    firstcall = .false.
!  END IF

        slpmag = 0.
        slpdir = 0.
!!!  CALL sfcslp( nx,ny, hterain, slpmag,slpdir, tem1,tem2,dx,dy )

!!! not considering map projections, for now:  GHB, 100720
!  IF ( mapproj == 0 ) THEN
    DO j=1,ny
      DO i=1,nx
        latscl(i,j) = ctrlat
        lonscl(i,j) = ctrlon
      END DO
    END DO
!  ELSE
!    DO j=1,ny
!      ys = 0.5*(y(j)+y(j+1))
!      DO i=1,nx
!        xs = 0.5*(x(i)+x(i+1))
!        CALL xytoll(1,1, xs,ys, latscl(i,j), lonscl(i,j))
!      END DO
!    END DO
!  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the local time at each grid point. The
!  following formula is based on that the input time is the GMT
!  time at the reference grid point of the center.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx

      latscl(i,j) = deg2rad * latscl(i,j)        ! lat: -90 to 90

      tloc(i,j) = hour0 + (curtim-dtbig)/3600.0                         &
                + lonscl(i,j)/15.0

      rjday(i,j) = jday + INT( tloc(i,j)/24.0 )

      tloc(i,j) = MOD( tloc(i,j), 24.0 )

      IF ( tloc(i,j) < 0. ) THEN
        tloc(i,j) = tloc(i,j) + 24.0            ! Local time
        rjday(i,j) = MOD( rjday(i,j)-1, yrday ) ! Julian day at each pts
      END IF

    END DO
  END DO

!!!  print *,'  curtim,dtbig = ',curtim,dtbig
!!!  print *,'  tloc,rjday,latscl = ',tloc(2,2),rjday(2,2),latscl(2,2)

  DO j=1,ny
    DO i=1,nx
      anncyc = pi2 * ( rjday(i,j) - 1.0 ) / yrday

      a2dr2(i,j) = 1.000110                                             &
                 + 0.034221 * COS(anncyc)                               &
                 + 0.00128  * SIN(anncyc)                               &
                 + 0.000719 * COS(2.*anncyc)                            &
                 + 0.000077 * SIN(2.*anncyc)  ! PX, Eq. 17

      hr = hour + minute / 60.

!  days before (-ve) or after (+ve) 1/1/2000
      days2k = 367 * year - 7 * ( year + ( month + 9 ) / 12 ) / 4       &
              + 275 * month / 9 + day - 730531.5 + hr / 24.

      lsun = 280.461 + 0.9856474 * days2k     ! Mean Longitude of the Sun
      950     IF ( lsun < 0 ) THEN
        lsun = lsun + 360.
        GO TO 950
      ELSE IF ( lsun > 360 ) THEN
        lsun = lsun - 360.
        GO TO 950
      END IF

      gsun = 357.528 + 0.9856003 * days2k     ! Mean anomaly of the Sun
      960     IF ( gsun < 0 ) THEN
        gsun = gsun + 360.
        GO TO 960
      ELSE IF ( gsun > 360 ) THEN
        gsun = gsun - 360.
        GO TO 960
      END IF

      lambda = lsun + 1.915 * SIN(gsun*deg2rad)  & ! Ecliptic longitude
          + 0.02 * SIN(2*gsun*deg2rad)
      970     IF ( lambda < 0 ) THEN
        lambda = lambda + 360.
        GO TO 970
      ELSE IF ( lambda > 360 ) THEN
        lambda = lambda - 360.
        GO TO 970
      END IF

      obliq = 23.439 - 0.0000004 * days2k     ! Obliquity of the ecliptic

      xsun = COS(lambda*deg2rad)
      ysun = COS(obliq*deg2rad) * SIN(lambda*deg2rad)
      asun = ATAN(ysun/xsun)*rad2deg
      IF ( xsun < 0. ) THEN
        alpha = asun + 180   ! Right Ascension (RA)
      ELSE IF ( ( ysun < 0. ) .AND. ( xsun > 0. ) ) THEN
        alpha = asun + 360
      ELSE
        alpha = asun
      END IF

      etau = ( lsun - alpha ) * 4. / 60.      ! Equation of time in hour

!    etau = 0.158 * sin( pi*(rjday(i,j)+10.)/91.25 ) ! Equation of time
!    :       + 0.125 * sin( pi*rjday(i,j)/182.5 )       ! Wong, Eq. 8

      shrangl = 15.0 * deg2rad                        & ! Hour angle
                     * ( tloc(i,j) + etau - 12.0)       ! Wong, Eq. 7

!    sdeclin = 23.5 * deg2rad
!    :          * cos( 2.0*pi*(rjday(i,j)-173.)/yrday ) ! Wong, Eq. 6
      sdeclin = ASIN(SIN(obliq*deg2rad)*SIN(lambda*deg2rad))
                                                ! Declination (in radian)

      cosz(i,j) = COS(latscl(i,j)) * COS(sdeclin) * COS(shrangl)        &
                + SIN(latscl(i,j)) * SIN(sdeclin)

!    print *, cos(latscl(i,j)),cos(sdeclin),cos(shrangl)
!    print *, sin(latscl(i,j)),sin(sdeclin)
!    print *,sdeclin,shrangl

      sinz = SIN ( ACOS(cosz(i,j)) )
!
!-----------------------------------------------------------------------
!
!  Consider the effects of the terrain slope on the solar radiation.
!  The slope magnitude and direction has been computed by subroutine
!  SFCSLP and passed in by slpmag and slpdir.
!
!-----------------------------------------------------------------------
!

      sinpsi = COS(sdeclin) * SIN(shrangl) *COS(latscl(i,j))
      cospsi = COSZ(i,j) * SIN( latscl(i,j) ) - SIN( sdeclin )
      azimuth = ATAN2( sinpsi, cospsi)

      dpsi = azimuth - slpdir(i,j)

      cosss(i,j) = COS( slpmag(i,j) ) * cosz(i,j)                       &
                 + SIN( slpmag(i,j) ) * sinz * COS( dpsi )

      cosz (i,j) = MAX( cosz (i,j), 0.0 )
      cosss(i,j) = MAX( cosss(i,j), 0.0 )

! added to calculate saltitude and sazimuth later for shading
!     IF( radshade /= 0) THEN  want azimuth regardless of shading

        IF (i == nxmid .AND. j == nymid) THEN
  
! augustin
! computes the solar altitudes and azimuth,
! sazimuth: 0=North; pi/2=East; Pi=South; 3*pi/2=West
! see the reference listed in the header of subroutine ZENANGL for more 
! information

        saltitude = pi/2-ACOS(cosz(i,j))

!  the following code is not computation safe (division by zero) DBW
!       sazimuth  = ACOS ( (SIN(sdeclin) * COS (latscl(i,j))             &
!                   - COS(sdeclin) * SIN (latscl(i,j))                   &
!                      * COS (shrangl))/sinz)

!  the above is replaced with:   DBW

        sinpsi = COS(sdeclin) * SIN(shrangl) * COS(latscl(i,j))
        cospsi = - COSZ(i,j) * SIN( latscl(i,j) ) + SIN( sdeclin )
        sazimuth = ABS( ATAN2(sinpsi, cospsi) )

!  end or replacement code  DBW

        shrangl_mid = shrangl

       END IF

 !     END IF  !  end of radshade if block..!want azimuth regardless of rshade

    END DO
  END DO

 ! IF( radshade /= 0) THEN

! broadcast these value from source to all other processors

! GHB, 100720
!!!    CALL mpbcastr(saltitude,   source)
!!!    CALL mpbcastr(sazimuth,    source)
!!!    CALL mpbcastr(shrangl_mid, source)

! modification of shrangl

! In the definition of the hour angle, shrangl  should be
! between -pi and 0 before solar noon and between 0 and pi after.
! Before this error, zenangl gave positive values
! between pi and 2*pi of shrangl
! before solar noon for some latitudes and day.
! for instance, +45 deg north March 21 2001
! it causes not problem in the original ARPS since it is just 
! a problem of module but subroutine shade needs to have shrangl 
! between -pi and pi


    IF ((shrangl_mid > pi) .AND. (shrangl_mid < pi2)) THEN
      shrangl_mid = shrangl_mid - pi2
    END IF
    IF ((shrangl_mid > pi2) .AND. (shrangl_mid < pi2+pi)) THEN
      shrangl_mid = shrangl_mid - pi2
    END IF

    IF (shrangl_mid > 0) THEN
      sazimuth = pi2 - sazimuth
    END IF

 ! END IF  !  end of radshade if block...

  RETURN
END SUBROUTINE zenangl

  END MODULE radtrns3d_module
