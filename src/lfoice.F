  MODULE lfoice_module

  implicit none

  private
  public :: lfo_ice_drive,lfoice_init

!-------------------------------------------------------------------------------
! This is from the 091406 version of lfoice.F which defines all constants within the wrapper
! It seems to work without problems.

!
! User-Changeable parameters
!
      real cwdiap, qcmincwrn, cwdisp ; parameter( cwdiap=20.e-6, qcmincwrn = 2.0e-3, cwdisp = 0.15 )
      integer autoconversion ; parameter( autoconversion = 0 )
      real rwdn, cnor ; parameter(rwdn = 1000., cnor = 8.0e6) ! rain params
      real swdn, cnos ; parameter(swdn =  100., cnos = 3.0e6) ! snow params
      real hwdn, cnoh ; parameter(hwdn =  900., cnoh = 4.0e4) ! hail params
!

 
! Air pressure, density, & Exner function
      real rho00

! Biggs
      real arz, brz

! Bergeron Process

      real cs10(32),cs11(32)
      real bfa1(32),bfa2(32),cbtim(32)
      real cmn, cmi40, cmi50, ri50, vti50, atemp, cm50a, cm40b, cm50b, bsfw
      integer ibb
!   Koenig (1971) was corrected at T=-14 C (fifteenth index of bfa1) from 0.1725e-6 to 0.1725e-4 after a reviewer noted a typo.
      data bfa1/0., 0.7939E-7, 0.7841E-6, 0.3369E-5, 0.4336E-5, 0.5285E-5, 0.3728E-5, 0.1852E-5, 0.2991E-6, 0.4248E-6,   &
         0.7434E-6, 0.1812E-5, 0.4394E-5, 0.9145E-5, 0.1725E-4, 0.3348E-4, 0.1725E-4, 0.9175E-5, 0.4412E-5, 0.2252E-5,   &
         0.9115E-6, 0.4876E-6, 0.3473E-6, 0.4758E-6, 0.6306E-6, 0.8573E-6, 0.7868E-6, 0.7192E-6, 0.6513E-6, 0.5956E-6,   &
         0.5333E-6, 0.4834E-6/
      data bfa2/0., 0.4006,    0.4831,    0.5320,    0.5307,    0.5319,    0.5249,    0.4888,    0.3894,    0.4047,      &
         0.4318,    0.4771,    0.5183,    0.5463,    0.5651,    0.5813,    0.5655,    0.5478,    0.5203,    0.4906,      &
         0.4447,    0.4126,    0.3960,    0.4149,    0.4320,    0.4506,    0.4483,    0.4460,    0.4433,    0.4413,      &
         0.4382,    0.4361/

! Accretions

      real xacwi
      real xcnos,xcnoh,xcnor
      real xxacx

! Evaporation/Deposition/Sublimation

      real xrcev1

! Initiation of cloud ice

      real xiint

! Initiation

      real cwmasn,cwmasx
      real cimasn,cimasx

! Melting/Wet growth of hail, Melting snow

      real xvth3
      real xmlt3
      real xxmlt1

!  Distribution parameters, Fallout, Mean diameter's, mass, and mixing ratio's

      real xslop
 

! Fallout, fall velocity parameters/vars

      real xvtr,xvts,xvth1


! Misc constants

      real advisc0,advisc1
      real elv,elf,els
      real ar,br,bta1,cs,ds
      real cbi,cbw,cnit,c1f3,cwc1
      real gf6,gf5,gf4,gf3,gf2,gf1,gf4ds,gf4br,gf4p5,gf5ds,gf5br,gf5p5,gf2p75
      real pid4
      real dragh
      real tfr,thnuc,tfrcbw,tfrcbi,tka0
      real cw, cpinv


! Saturation lookup table, vapor pressures, ratio's

      integer nqsat, lll
      parameter(nqsat=20001)
      real tabqvs(nqsat),tabqis(nqsat)
      real fqsat, temq
      parameter (fqsat=0.01)
      real cai, caw


! Saturation Adjustment

      real cqv1,cqv2


! Species threshold,  mass,    mass constraints

      real qcmin,qhmin,qimin,qrmin,qsmin
      real qccrit, qscrit, qicrit

! Functions

      real gamma
      external gamma

!-----------------------------------------
!  stuff added by George Bryan, 060929:
      real :: poo,pi,dnz00,cpi,cwdn
      integer :: outfile
      real g, cp, cv, rd, rw, rcp ; parameter (g=9.81, cp=1004.0, rd=287.0, cv=cp-rd, rw=461.5, rcp=rd/cp)
!-----------------------------------------


  CONTAINS

!-------------------------------------------------------------------------------
!   lfo_ice_drive
!
!   PURPOSE: To transform the CM1 variables to what 3-ICE expects
!   DATE:   Sept 12, 2006 
!   AUTHOR: Matt Gilmore
!
!-------------------------------------------------------------------------------
!           call lfo_ice_drive(dt, mf, pi0, pp3d, prs, th0, th3d, qv0, rho0, q3d) 

      subroutine lfo_ice_drive(dt, mf, pi0, prs0, pp3d, prs, th0, th3d,    &
                               qv0, rho0, q3d, vq, pn)

      use input
      implicit none


! Variable declarations for variables passed from CM1

!(CM1 input.incl common block)
!       ni       # grid points (for tile) in x-direction (east/west)
!       nj       # grid points (for tile) in y-direction (north/south)
!       nk       # grid points (for tile) in vertical
!       dt       model large timestep     (s)
!       rdz      inverse of constant model dz (m^-1)
 
! CM1 INPUT (passed)
!       mf     = 1/[rdz*(zh(k)-zh(k-1)]   (No Dim)    (i,j,k)
!  Note rdz*mf is 1/dzactual for the stretched grid ! zh is scalar ht (m)
!
!       pi0      base state total Exner   (No Dim)    (i,j,k)
!       pp3d     perturbation Exner pressure          (i,j,k)
!       prs      full pressure            (Pascals)   (i,j,k)
!       th0      base state potential temp(K)         (i,j,k)
!       th3d     potential temperature perturb. (K)   (i,j,k)
!       rho0     base state dry air density     (kg/m^3)    (i,j,k)   ! replaced db(kz) in NCOMMAS version
!       qv0      base state vapor mixing ratio        (i,j,k)
!
!       q3d      qv, qc, qr, qi, qs, qh  (kg/kg)      (i,j,k,numq)

!      include 'input.incl'     !  dt, ib, ie, jb, je, kb, ke, ibm, ibi, jbm, jbi, kbm, kbi, iem ,iei, jem, jei, kem, kei, ni, nj, nk, numq
!       ni        number of grid pts in the x-direction (east/west)
!       nj        number of grid pts in the y-direction (north/south)
!       nk        number of grid pts in the vertical
!       ib    -2  starting point including 3-point buffer
!       ie  ni+3
!       jb    -2
!       je  nj+3
!       kb     0
!       ke  nk+1


      real :: dt
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,pp3d,prs,th3d,th0,   &
                                            rho0,qv0,pn
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d,vq
!
! The following are 3-ICE variables/terminology
! We set the dimensions to that of the tile without buffer (not the full domain)
!
      integer nstep          ! Current Model timestep
      parameter (nstep = 1)  !    for now, set to a constant
      real dzc(nk)           ! 1/dz for the stretched grid 
      real rinit(nk)         ! base state density (kg/m^3)
      real pinit(nk)         ! base state pressure (Pa)
      real sb(nk, numq)        ! base state array of th, qv, qc, etc

!-------------------------------------------------------------------------------
!     misc. local variables
!-------------------------------------------------------------------------------

      integer i,j,k,n

!-------------------------------------------------------------------------------

!
! Set up Conversions between CM1 and SAM
!
      sb(:,:)  = 0.0
      rinit(:) = 0.0
      pinit(:) = 0.0
      dzc(:)   = 0.0

!     Set up 1D vertical arrays.  Use corner of domain for no-terrain case.

      i = 1
      j = 1
      DO k = 1,nk                   !MSG might need to change to nk instead of nk-1
       dzc(k)  =   mf(i,j,k)*rdz
       sb(k,1) =  th0(i,j,k)
       sb(k,2) =  qv0(i,j,k)
       rinit(k)= rho0(i,j,k)
       pinit(k)= prs0(i,j,k)
      ENDDO

!     Copy 3-D arrays from CM1 to 3-ICE.

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      DO k = 1, nk
      DO j = 1, nj
      DO i = 1, ni
       pn(i,j,k)   = prs(i,j,k) - pinit(k)       !perturbation pressure (Pascals)
      ENDDO
      ENDDO
      ENDDO

!     Make sure vq array is zero, by default

      DO n = 1, numq
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        DO k = 1, nk
        DO j = 1, nj
        DO i = 1, ni
          vq(i,j,k,n) = 0.0
        ENDDO
        ENDDO
        ENDDO
      ENDDO

! DO MICROPHYSICS

           CALL LFO_ICE (ni,nj,nk,numq,dt,dzc,sb,th3d,q3d,vq,rinit,pn,pinit,nstep, &
                         ib,ie,jb,je,kb,ke,ibm,iem,jbm,jem,kbm,kem)

!-------------------------------------------------------------------------------

! Uncomment the next lines if micro is allowed to feedback on pressure
!  This could also be put directly into LFO_ICE
!     DO k = 1, nk
!     DO j = 1, nj
!     DO i = 1, ni
!       prs(i,j,k) = pn(i,j,k) + prs0(i,j,k)    !total pressure (Pascals)
!       pp3d(i,j,k) = (prs(i,j,k)*rp00)**rovcp  !total pressure (Ekner)
!     ENDDO
!     ENDDO
!     ENDDO

      if(timestats.ge.1) time_microphy=time_microphy+mytime()


      END subroutine lfo_ice_drive
!--------------------------------------------------------------------------
!
!
!  3-ICE MODEL
!
!  VERSION: 1.3 with Rates (2/10/05)
!
!  LIN-FARELY-ORVILLE-like "Simple Ice and Liquid Microphysics Scheme"
!   based on Gilmore et al. (2004) Mon. Wea. Rev.
!
!   Copyright (C) <2004>  <Jerry Straka and Matthew Gilmore>
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Lesser General Public
!   License as published by the Free Software Foundation; either
!   version 2.1 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Lesser General Public License for more details.
!
!   If you find this code useful and publish results using it, please reference:
!
!     Gilmore M. S., J. M. Straka, and E. N. Rasmussen,
!           Monthly Weather Review: Vol. 132, No. 8, pp. 1897-1916.
!
!--------------------------------------------------------------------------
!          CALL LFO_ICE (ni,nj,nk,numq,dt,dzc,sb,th3d,q3d,rinit,pn,pinit,nstep, &
!                        ib,ie,jb,je,kb,ke,ibm,iem,jbm,jem,kbm,kem)
      SUBROUTINE LFO_ICE(nx,ny,nz,numq,dtp,dzc,ab,th3d,q3d,vq,   db,pn,   pb,nstep, &
                         ib,ie,jb,je,kb,ke,ibm,iem,jbm,jem,kbm,kem)
!--------------------------------------------------------------------------
! RELEASE NOTES
!
! MSG - 2/10/05 Ported the SAM gather-scatter version to NCOMMAS
!       This is most similar to that actually used in the published manuscripts.
!       This version also provides rate output (labels consistent with Gilmore et al) for post-analysis.
!       Known inconsistencies/omissions: 
!          1) min q criteria applied to some (but not all) processes
!          2) internal rate names and signs differ compared to Gilmore et al. (2004),
!               however, rate output is consistent
!          3) time splitting for fallout omitted
!
!
! MSG - 5/5/05  Added time-splitting fallout.
!
! MSG - 9/12/06c Ported code to CM1 with wrapper routine.  Strategy is to *not* use any common blocks
!                or include files.  Thus, slight mismatch between CM1 constants and herein.  
!                Also, here we assume 1-D base state (faster code) but CM1 is really 3-D base state.
!                See the alternate 9/13/06 version if you wish to use the full 3D base state.
!                Note that this code has internal saturation adjustment and fallout.
!
! GHB - 061021: - Now using common block, using include file 'lfoice.incl'
!                 Several constants are now defined in 'lfoice.incl'
!                 Call to "lfoice_init" has been added to "param"
!
!--------------------------------------------------------------------------
!                                
!
!  general declarations
!
!
      implicit none

!!!      include 'lfoice.incl'

!     include 'param.h'
!--------------------------------------------------------------------------

      integer ix,jy,kz,nx,ny,nz

      integer ib, ie, jb, je, kb, ke
      integer ibm,iem,jbm,jem,kbm,kem
      integer numq

      integer nstep
      real dtp
      integer istag,jstag,kstag
      parameter (istag=0, jstag=0, kstag=0)            !MSG Set these to 1 for NCOMMAS; set to 0 for CM1
      integer      lv,  lc,   lr,   li,   ls,   lh
      parameter (  lv=1,lc=2, lr=3, li=4, ls=5, lh=6)  !MSG indicee numbers consistent with CM1

      real dzc(nz)                                     !MSG weren't defined before

!     real  an(nx,ny,nz,numq)                                !MSG scalars (NCOMMAS)
      real, dimension(ib:ie,jb:je,kb:ke) :: th3d             !MSG TH scalar (CM1)
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d   !MSG Q scalars (CM1)
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: vq    !GHB Q fallspeeds (CM1)
      real  ab(nz,numq)                                      !MSG base state for scalars
      real, dimension(ib:ie,jb:je,kb:ke) :: pn               !MSG perturb pressure
      real pb(nz), db(nz)                                    !MSG base state pressure, density
!
!---
! 
!  declarations microphysics and for gather/scatter
!
!   (MSG: In this version, if ngs, the number of gather/scatter points, is smaller than the number of
!         microphysically active gridpoints in an XZ slice, then the code does multiple
!         gather/scatter iterations within a single XZ slice using inumgs and the pointers nxmpb,nzmpb.)
!
      integer jgs,mgs,ngs,numgs,inumgs,nxmpb,nzmpb,nxz      !MSG added innumgs
      parameter ( ngs=500 )
      integer ngscnt,igs(ngs),kgs(ngs)

! Air and particle temperatures (equivalent to each other in LFO version)

      real temp(ngs)
      real temg(ngs),temcg(ngs),theta(ngs),thetap(ngs),theta0(ngs)

! Air pressure, density, & Exner function

      real pbz(ngs),pres(ngs),presp(ngs),pres0(ngs)
      real rho0(ngs),dnz(ngs),pi0n(ngs),piz(ngs)

! Accretions

      real qiacr(ngs)
      real qracw(ngs), qraci(ngs), qracs(ngs)
      real qsacw(ngs), qsaci(ngs),            qsacr(ngs)
      real qhacw(ngs), qhaci(ngs), qhacs(ngs),qhacr(ngs)
      real eri(ngs),erw(ngs),ers(ngs)
      real esw(ngs),esi(ngs)
      real ehw(ngs),ehr(ngs),ehi(ngs),ehs(ngs)

! Biggs Freezing

      real qrfrz(ngs), xrfrz(ngs)

! Bergeron Process

      real qsfw(ngs),qsfi(ngs), eic(ngs)
      real cs9

! Conversions

      real qhcns(ngs), qrcnw(ngs), qscni(ngs), qdiff,argrcnw
      real ehscnv(ngs), esicnv(ngs)

! Evaporation/Deposition/Sublimation

      real qhdsv(ngs),qhdpv(ngs),qhsbv(ngs)
      real qsdsv(ngs),qsdpv(ngs),qssbv(ngs)
      real qrcev(ngs)
      real xav(ngs),xbv,xas(ngs),xbs
      real xrcev2(ngs)
      real xce,xrv,xds
      real xxdsv,xxcev

! Initiation of cloud ice

      real dqisdt(ngs),qiint(ngs),cnnt

! Melting/Wet growth of hail, Melting snow

      real qhmlr(ngs), qsmlr(ngs)
      real xmlt1,xmlt2
      real xhmlt2(ngs), xsmlt2(ngs)
      real xsv, xhv, xhsw

! Wet/dry growth, shedding

      real qhdry(ngs),qhwet(ngs)
      real xhwet1,xhwet2, xcwt,xwt1
      real qhacip(ngs),qhacsp(ngs),qhshr(ngs)

! Water Budgets

      real ptotal(ngs),ptotsat(ngs) 
      real pqcwi(ngs),pqcii(ngs),pqrwi(ngs)
      real pqswi(ngs),pqhwi(ngs),pqwvi(ngs)
      real pqcwd(ngs),pqcid(ngs),pqrwd(ngs)
      real pqswd(ngs),pqhwd(ngs),pqwvd(ngs)
      integer il2(ngs),il3(ngs),il5(ngs)

! Latent Heating Computation

      real psub(ngs),pvap(ngs),pfrz(ngs),ptem(ngs)

! Maximum Depletion Tendencies

      real qc5dt,qi5dt,qr5dt,qs5dt,qh5dt

! Flags
!!!      integer imake                    !prints intercept/density to output file
!!!      save imake
!!!      data imake / 0 /
      integer  ndebug, nrates          !prints debug stuff to output file, prints rates to out file
      integer itfall
      parameter (ndebug = 0, nrates=0) !prints debugging, rate flags
        ! GHB, 061013:  moved fallout to CM1 solve
      parameter (itfall = 0)           !0 no flux fallout, 2: timesplit flux fallout, 1 flux on regular dt

! Fallout, fall velocity parameters/vars

      real dtz1
      ! GHB, 061013:  fallout calculations have been moved to CM1 solve
!!!      real cwflx(ngs), cflux(nx,nz)
!!!      real piflx(ngs), pflux(nx,nz)
!!!      real rwflx(ngs), rflux(nx,nz)
!!!      real swflx(ngs), sflux(nx,nz)
!!!      real hwflx(ngs), hflux(nx,nz)

      real maxfall            ! max fallspeed of qg, qh, qf, qr, qd, qm
      real vtrdzmax           ! input courant                         ! MSG
      real dtsplit            ! small timestep for rain sedimentation ! MSG
      integer nrnstp,inrnstp  ! number of small timesteps, counter    ! MSG

 
!  Distribution parameters, Fallout, Mean diameter's, mass, and mixing ratio's
 
!
!
      real qwv(ngs),qcw(ngs),qci(ngs),qrw(ngs),qsw(ngs),qhw(ngs)
      real ccw(ngs),cci(ngs),crw(ngs),csw(ngs),chw(ngs)
!
      real vtwbar(ngs),cwmas(ngs)
      real vtibar(ngs),cimas(ngs)
      real vtrbar(ngs)
      real vtsbar(ngs)
      real vthbar(ngs)
!
      real cwdia(ngs),cwdia2(ngs)
      real cidia(ngs),cidia2(ngs)
      real rwdia(ngs),rwdia2(ngs)
      real swdia(ngs),swdia2(ngs)
      real hwdia(ngs),hwdia2(ngs)

! Saturation Adjustment

      real qvap(ngs)
      real dqvcnd(ngs),dqwv(ngs),dqcw(ngs),dqci(ngs)
      real gamss,denom1,denom2
      real gamw(ngs),gams(ngs)
      real qwvp(ngs), qv0n(ngs)
      real fraci(ngs),fracl(ngs)
      real cdw, cdi
      integer itertd
      real qwfzi(ngs),qimlw(ngs)
      real qidep(ngs),qisub(ngs)
      real qcevp(ngs), qccnd(ngs)
      real qcevpcnd(ngs), qisubdep(ngs)

! Saturation lookup table, vapor pressures, ratio's

      integer ltemq
      real qvs(ngs),qis(ngs),qss(ngs),pqs(ngs)
      real tsqr(ngs),ssi(ngs),ssw(ngs)

! Misc Variables

      real advisc,schm,tkarwinv,rhoratio,mlttemp
      real wvdf(ngs),akvisc,ci(ngs),tka(ngs)
      real cc3(ngs),cc4(ngs),cc5(ngs)

!  Rate Output (Domain-total g/m^3)

      real hfrz, hdep, hcnd, cevap, cmelt, csub		!MSG added total cooling/heating vars
      real tqva,  tqia,  tqca
      real tqvb,  tqib,  tqcb
      real tqvap, tqsap, tqiap, tqrap, tqcap, tqhap
      real tqvbp, tqsbp, tqibp, tqrbp, tqcbp, tqhbp
      real tvsum, tssum, tisum, trsum, tcsum, thsum, tqc,tqi,tqv,tsumall 
      real suma,  sumb,  psum 
      real trqsacw, trqhacr, trqhshr, trqhmlr, trqsmlr, tiqiint, tiqidep, trqhacs
      real trqhacw, trqhaci, tvqhsbv, tvqcevp, tvqssbv, tvqrcev, trqrcnw, trqracw
      real tcqcmli, tvqisub, tcqccnd, tiqifzc, thqhacs, thqsacr, thqhaci, thqhacr
      real thqhacw, thqhdpv, thqrfrz, thqiacr, thqracs, thqraci, tsqsacw, tsqsacr
      real tsqscni,  tsqsfi,  tsqsfw, tsqsdpv, thqhcns, tsqsaci, tsqraci, tsqiacr, thqhwet 
!     
!
!  read in constants from 'inmicro.jmslfo'
!  [deleted this part - MSG]

      if ( ndebug .eq. 1 ) print*,'Just entered micro...'

!
!  ZERO
!
!
!  totals for source / sink terms
!
!
!  vapor
!
      tvqrcev = 0.0
      tvqssbv = 0.0
      tvqhsbv = 0.0
      tvqcevp = 0.0
      tvqisub = 0.0
!
!  cloud water
!
      tcqccnd = 0.0 
      tcqcmli = 0.0
!
!  rain
!
      trqrcnw = 0.0
      trqracw = 0.0
      trqhmlr = 0.0
      trqsmlr = 0.0
      trqhshr = 0.0
      trqsacw = 0.0
      trqhacr = 0.0
      trqhacw = 0.0
      trqhaci = 0.0
      trqhacs = 0.0

!
!  cloud ice
!
      tiqiint = 0.0
      tiqidep = 0.0
      tiqifzc = 0.0
!
!  snow
!
      tsqsfi  = 0.0
      tsqsfw  = 0.0
      tsqscni = 0.0
      tsqsacw = 0.0
      tsqsacr = 0.0
      tsqraci = 0.0
      tsqiacr = 0.0
      tsqsaci = 0.0
      tsqsdpv = 0.0
!
!  hail
!
      thqhcns = 0.0
      thqhacr = 0.0
      thqhacw = 0.0
      thqhaci = 0.0
      thqhacs = 0.0
      thqsacr = 0.0
      thqracs = 0.0
      thqraci = 0.0
      thqiacr = 0.0
      thqhdpv = 0.0
      thqrfrz = 0.0
      thqhwet = 0.0
!
!  total heating and cooling rates
!
      hfrz    = 0.0
      hdep    = 0.0
      hcnd    = 0.0
      cevap   = 0.0
      cmelt   = 0.0
      csub    = 0.0
!
!
!  various rate budgets
!
      tqvap     = 0.0
      tqcap     = 0.0
      tqiap     = 0.0
      tqrap     = 0.0
      tqsap     = 0.0
      tqhap     = 0.0
!
      tqvbp     = 0.0
      tqcbp     = 0.0
      tqibp     = 0.0
      tqrbp     = 0.0
      tqsbp     = 0.0
      tqhbp     = 0.0
!
      tqva     = 0.0
      tqca     = 0.0
      tqia     = 0.0
!
      tqvb     = 0.0
      tqcb     = 0.0
      tqib     = 0.0
!
      suma      = 0.0
      psum      = 0.0
!
      sumb      = 0.0
!
      tvsum    = 0.0
      tcsum    = 0.0
      tisum    = 0.0
      trsum    = 0.0
      tssum    = 0.0
      thsum    = 0.0
!
      tsumall   = 0.0
      tqv       = 0.0
      tqc       = 0.0
      tqi       = 0.0
!
!
!  end of totals
!
!
!  other constants
!


!      write(6,*) '-----------------------------------------------------------------------'
!
      if ( ndebug .eq. 1 ) print*,'dbg = 0b'

!
!  MSG: These 3 things a function of base state (eventually will be moved to 3-D base state)
!
      do 10 kz = 1,nz-kstag
       pbz(kz) = pb(kz)
       piz(kz) = (pbz(kz)/poo)**rcp
       temp(kz) = piz(kz)*ab(kz,1)
  10  continue
      
      if (ndebug .eq. 1 ) print*,'dbg = 1'
!
      if (ndebug .eq. 1 ) print*,'dbg = 2'
!
!  start jy loop  (for doing XZ slabs)
!
!
      do 9999 jy = 1,ny-jstag
!
!  VERY IMPORTANT:  SET jgs
!
      jgs = jy
!
!  zero precip flux arrays
!
      ! GHB, 061013:  fallout calculations have been moved to CM1 solve
!!!      vtrdzmax = 0.0                     !MSG Courant for the slice
!!!      maxfall = 0.0                      !MSG max fallspeed of any particle for the slice
!!!      if (ndebug .eq. 1 ) print*,'dbg = 3'
!!!      do 95 kz = 1,nz-kstag
!!!      do 96 ix = 1,nx-istag
!!!      hflux(ix,kz) = 0.0
!!!      cflux(ix,kz) = 0.0
!!!      pflux(ix,kz) = 0.0
!!!      rflux(ix,kz) = 0.0
!!!      sflux(ix,kz) = 0.0
!!!  96  continue
!!!  95  continue
!
!..gather microphysics  
!
      if (ndebug .eq. 1 ) print*,'dbg = 4'
      nxmpb = 1
      nzmpb = 1
      nxz = nx*nz
      numgs = nxz/ngs + 1
      do 1000 inumgs = 1,numgs 
      ngscnt = 0
! 061023, GHB:  dunno why the "kstag-1" is here ... removing the -1
!!!      do kz = nzmpb,nz-kstag-1 
      do kz = nzmpb,nz-kstag
      do ix = nxmpb,nx-istag

      theta(kz) = th3d(ix,jy,kz) + ab(kz,1)
      temg(kz) = theta(kz)*( (pn(ix,jy,kz)+pbz(kz)) / poo ) ** rcp
      ltemq = nint((temg(kz)-163.15)/fqsat+1.5)
      ltemq = min(max(ltemq,1),nqsat)         
      pqs(kz) = 380.0/(pn(ix,jy,kz)+pbz(kz))
      qvs(kz) = pqs(kz)*tabqvs(ltemq)        
      qis(kz) = pqs(kz)*tabqis(ltemq)       
      
      if ( temg(kz) .lt. tfr ) then

       qcw(kz) = max(q3d(ix,jy,kz,lc) ,0.0)
       qci(kz) = max(q3d(ix,jy,kz,li) ,0.0)

       if( qcw(kz) .ge. 0.0 .and. qci(kz) .eq. 0.0 ) qss(kz) = qvs(kz)
       if( qcw(kz) .eq. 0.0 .and. qci(kz) .gt. 0.0)  qss(kz) = qis(kz)
       if( qcw(kz) .gt. 0.0 .and. qci(kz) .gt. 0.0)  qss(kz) = (qcw(kz)*qvs(kz) + qci(kz)*qis(kz)) /(qcw(kz) + qci(kz))
      else
       qss(kz) = qvs(kz)
      end if
!
      if ( q3d(ix,jy,kz,lv) .gt. qss(kz) .or.   &
           q3d(ix,jy,kz,lc) .gt. qcmin  .or.    &
           q3d(ix,jy,kz,li) .gt. qimin .or.     &
           q3d(ix,jy,kz,lr) .gt. qrmin .or.     &
           q3d(ix,jy,kz,ls) .gt. qsmin .or.     &
           q3d(ix,jy,kz,lh) .gt. qhmin ) then
      ngscnt = ngscnt + 1
      igs(ngscnt) = ix
      kgs(ngscnt) = kz
      if ( ngscnt .eq. ngs ) goto 1100
      end if
      ENDDO     !MSG - i loop
      nxmpb = 1
      ENDDO     !MSG - k loop
 1100 continue
      if ( ngscnt .eq. 0 ) go to 9998
      if ( ndebug .eq. 1 ) print*,'dbg = 5'
!
!  define temporaries to be used in calculations
!
      do 1010 mgs = 1,ngscnt
      dnz(mgs)   = db(kgs(mgs))
      pres0(mgs) = pbz(kgs(mgs))                    !MSG need to make 3-D base state eventually
      presp(mgs) = pn(igs(mgs),jy,kgs(mgs))
      pres(mgs)  = presp(mgs) + pres0(mgs)
      pi0n(mgs)  = piz(kgs(mgs))                    !MSG need to make 3-D base state eventually
       cc3(mgs)  = cpi*elf/pi0n(mgs)
       cc4(mgs)  = cpi*elv/pi0n(mgs)
       cc5(mgs)  = cpi*els/pi0n(mgs)
      theta0(mgs)= ab(kgs(mgs),1)                   !MSG need to make 3-D base state eventually
      thetap(mgs)= th3d(igs(mgs),jy,kgs(mgs))
      theta(mgs) = thetap(mgs) + theta0(mgs)
      temg(mgs)  = theta(mgs)*( pres(mgs) / poo ) ** rcp      

!      if (IEEE_IS_NAN(temg(mgs))) call stopcm1

      temcg(mgs) = temg(mgs) - tfr
      pqs(mgs)   = 380.0/pres(mgs)
      ltemq      = nint((temg(mgs)-163.15)/fqsat+1.5)
      ltemq      = min(max(ltemq,1),nqsat)
      qvs(mgs)   = pqs(mgs)*tabqvs(ltemq)
      qis(mgs)   = pqs(mgs)*tabqis(ltemq)
      qwv(mgs)   = q3d(igs(mgs),jy,kgs(mgs),lv)
      qv0n(mgs)  = ab(kgs(mgs),2)                   !MSG need to make 3-D base state eventually
      qwvp(mgs)  = qwv(mgs) - qv0n(mgs) 
      qcw(mgs) = max(q3d(igs(mgs),jy,kgs(mgs),lc), 0.0) 
      qci(mgs) = max(q3d(igs(mgs),jy,kgs(mgs),li), 0.0) 
      qrw(mgs) = max(q3d(igs(mgs),jy,kgs(mgs),lr), 0.0) 
      qsw(mgs) = max(q3d(igs(mgs),jy,kgs(mgs),ls), 0.0) 
      qhw(mgs) = max(q3d(igs(mgs),jy,kgs(mgs),lh), 0.0) 
      il2(mgs) = 0
      il3(mgs) = 0
      il5(mgs) = 0
      if ( temg(mgs) .lt. tfr ) then 
       il5(mgs) = 1
       if ( qrw(mgs) .lt. 1.0e-04 .and. qsw(mgs) .lt. 1.0e-04 ) il2(mgs) = 1
       if ( qrw(mgs) .lt. 1.0e-04 ) il3(mgs) = 1
      end if
 1010 continue
!
!
! 
!  other constants for paramerization
!
      do mgs = 1,ngscnt
      advisc     = advisc0*(416.16/(temp(kgs(mgs))+120.0))*(temp(kgs(mgs))/296.0)**(1.5)
      akvisc     = advisc/dnz(mgs)
      ci(mgs)    = (2.118636 + 0.007371*(temp(kgs(mgs))-tfr))*(1.0e+03)     !MSG used in qhwet
      tka(mgs)   = tka0*advisc/advisc1
      wvdf(mgs)  = (2.11e-05)*((temp(kgs(mgs))/tfr)**1.94)*(101325.0/(pbz(kgs(mgs))))
      schm       = akvisc/wvdf(mgs)
      tkarwinv   = 1./(tka(mgs)*rw)
      xav(mgs)   = (elv**2)*tkarwinv
      xas(mgs)   = (els**2)*tkarwinv
      rhoratio   = (dnz00/dnz(mgs))**0.25
      mlttemp    = 0.308*(schm**(1./3.))*(akvisc**(-0.5))
      xhmlt2(mgs)= mlttemp*gf2p75
      xsmlt2(mgs)= mlttemp*gf5ds*(cs**(0.5))*rhoratio
      xrcev2(mgs)= mlttemp*gf5br*(ar**(0.5))*rhoratio
      enddo 

!
      if (ndebug .eq. 1 ) print*,'dbg = 6'
!
! cloud water variables
!
      do 4101 mgs = 1,ngscnt
!     ccw(mgs) = 1.e9    !LFO default (Western Plains)
      ccw(mgs) = .6e9    !Central plains CCN value
!     ccw(mgs) = .3e9    !Maritime CCN value
      cwmas(mgs) = min( max(qcw(mgs)*dnz(mgs)/ccw(mgs),cwmasn),cwmasx )
      cwdia(mgs) = (cwmas(mgs)*cwc1)**c1f3 
      cwdia2(mgs) = cwdia(mgs)**2
      vtwbar(mgs) = (ar*(cwdia(mgs)**br))*(rho00/dnz(mgs))**0.5
 4101 continue
!
! cloud ice variables
!
      do 4201 mgs = 1,ngscnt
!
      cimasx     = 3.23e-8
      cci(mgs)   = max(min(cnit*exp(-temcg(mgs)*bta1),1.e+09),1.0)       !Fletcher's formula
      cimas(mgs) = min( max(qci(mgs)*dnz(mgs)/cci(mgs),cimasn),cimasx )
      if ( temcg(mgs) .gt. 0 ) then
        cidia(mgs) = 0.0
      else
        cidia(mgs) = 16.7*(cimas(mgs)**(0.5))
        cidia(mgs) = max(cidia(mgs), 1.e-5)
        cidia(mgs) = min(cidia(mgs), 3.e-3)
      endif
      cidia2(mgs) = cidia(mgs)**2
      vtibar(mgs) = (cs*(cidia(mgs)**ds))*(rho00/dnz(mgs))**0.5
!
 4201 continue
!
!  mp-distribution information for rain, snow agg's, and graupel/hail
!
!
!  definitions for marshall palmer distribution variables
!  (rain, snow, hail) when mixing ratio only is predicted
!
      if (ndebug .eq. 1 ) print*,'dbg = 7a'
!
      do 4301 mgs = 1,ngscnt
!
      rwdia(mgs) = 1.e-20
      swdia(mgs) = 1.e-20
      hwdia(mgs) = 1.e-20
!
      if ( qrw(mgs) .gt. 1.0e-10 ) rwdia(mgs) = xslop*(dnz(mgs)*qrw(mgs)/(rwdn*xcnor))**(0.25) 
      if ( qsw(mgs) .gt. 1.0e-10 ) swdia(mgs) = xslop*(dnz(mgs)*qsw(mgs)/(swdn*xcnos))**(0.25)
      if ( qhw(mgs) .gt. 1.0e-10 ) hwdia(mgs) = xslop*(dnz(mgs)*qhw(mgs)/(hwdn*xcnoh))**(0.25)
!
      rwdia2(mgs) = rwdia(mgs)**2
      swdia2(mgs) = swdia(mgs)**2
      hwdia2(mgs) = hwdia(mgs)**2
!
      vtrbar(mgs) = xvtr*(dnz(mgs)**(-0.5))*(rwdia(mgs)**br)
      vtsbar(mgs) = xvts*(dnz(mgs)**(-0.5))*(swdia(mgs)**ds)
      vthbar(mgs) = xvth1*(dnz(mgs)**(-0.5))*(xvth3**2)*((hwdn*hwdia(mgs))**(0.5))
!     
      crw(mgs) = xcnor*rwdia(mgs)
      csw(mgs) = xcnos*swdia(mgs)
      chw(mgs) = xcnoh*hwdia(mgs)
!
 4301 continue

      ! GHB, 061013:  fallout calculations have been moved to CM1 solve
!!!      if (itfall .eq. 2) then
!!!        do mgs = 1,ngscnt                                        !MSG Find max fallspeed for later
!!!         maxfall = max(maxfall,vtrbar(mgs))
!!!         maxfall = max(maxfall,vthbar(mgs))
!!!        enddo
!!!      endif

!
      do 5001 mgs = 1,ngscnt 
!
!  maximum depletion tendency by any one source
!
      qc5dt = 0.20*qcw(mgs)/dtp
      qi5dt = 0.20*qci(mgs)/dtp
      qr5dt = 0.20*qrw(mgs)/dtp
      qs5dt = 0.20*qsw(mgs)/dtp
      qh5dt = 0.20*qhw(mgs)/dtp

!
!  collection efficiencies
!
      eic(mgs) = 1.0
      eri(mgs) = 1.0
      erw(mgs) = 1.0
      esi(mgs) = exp(0.025*min(temcg(mgs),0.0))
      esicnv(mgs) = esi(mgs)
      esw(mgs) = 1.0
      ers(mgs) = 1.0
      ehw(mgs) = 1.0
      ehr(mgs) = 1.0
      ehs(mgs) = exp(0.09*min(temcg(mgs),0.0))
      ehscnv(mgs) = ehs(mgs)
      if ( temcg(mgs) .gt. 0.0 ) ehs(mgs) = 1.0
      ehi(mgs) = 0.1

      if ( qcw(mgs) .lt. qcmin .or. qci(mgs) .lt. qimin ) eic(mgs) = 0.0
      if ( qrw(mgs) .lt. qrmin .or. qci(mgs) .lt. qimin ) eri(mgs) = 0.0
      if ( qrw(mgs) .lt. qrmin .or. qcw(mgs) .lt. qcmin ) erw(mgs) = 0.0
      if ( qsw(mgs) .lt. qsmin .or. qci(mgs) .lt. qimin ) esi(mgs) = 0.0
      if ( qsw(mgs) .lt. qsmin .or. qcw(mgs) .lt. qcmin ) esw(mgs) = 0.0
      if ( qsw(mgs) .lt. qsmin .or. qrw(mgs) .lt. qrmin ) ers(mgs) = 0.0
      if ( qhw(mgs) .lt. qhmin .or. qcw(mgs) .lt. qcmin ) ehw(mgs) = 0.0
      if ( qhw(mgs) .lt. qhmin .or. qrw(mgs) .lt. qrmin ) ehr(mgs) = 0.0
      if ( qhw(mgs) .lt. qhmin .or. qsw(mgs) .lt. qsmin ) ehs(mgs) = 0.0
      if ( qhw(mgs) .lt. qhmin .or. qci(mgs) .lt. qimin ) ehi(mgs) = 0.0
!
!  accretions:
!    marshall-palmer size distribution collection 
!    of constant size distribution
!      1)  sink for constant size distribution
!      2)  source for marshall-palmer size distribution
!
      qracw(mgs) =  &
         min(erw(mgs)*qcw(mgs)*xacwi*crw(mgs)*abs(vtrbar(mgs)-vtwbar(mgs))  &
        *(  gf3*rwdia2(mgs) + 2.0*gf2*rwdia(mgs)*cwdia(mgs) + gf1*cwdia2(mgs) )      , qc5dt)
      qsacw(mgs) =  &
         min(esw(mgs)*qcw(mgs)*xacwi*csw(mgs)*abs(vtsbar(mgs)-vtwbar(mgs))  &
        *(  gf3*swdia2(mgs) + 2.0*gf2*swdia(mgs)*cwdia(mgs) + gf1*cwdia2(mgs) )      , qc5dt)
      qhacw(mgs) =  &
         min(ehw(mgs)*qcw(mgs)*xacwi*chw(mgs)*abs(vthbar(mgs)-vtwbar(mgs))  &
        *(  gf3*hwdia2(mgs) + 2.0*gf2*hwdia(mgs)*cwdia(mgs) + gf1*cwdia2(mgs) )      , qc5dt)
      qraci(mgs) =  &
         min(eri(mgs)*qci(mgs)*xacwi*crw(mgs)*abs(vtrbar(mgs)-vtibar(mgs))  &
        *(  gf3*rwdia2(mgs) + 2.0*gf2*rwdia(mgs)*cidia(mgs) + gf1*cidia2(mgs) )      , qi5dt)
      qsaci(mgs) =  &
         min(esi(mgs)*qci(mgs)*xacwi*csw(mgs)*abs(vtsbar(mgs)-vtibar(mgs))  &
        *(  gf3*swdia2(mgs) + 2.0*gf2*swdia(mgs)*cidia(mgs) + gf1*cidia2(mgs) )      , qi5dt)
      qhaci(mgs) =  &
        min(ehi(mgs)*qci(mgs)*xacwi*chw(mgs)*abs(vthbar(mgs)-vtibar(mgs))   &
        *(  gf3*hwdia2(mgs) + 2.0*gf2*hwdia(mgs)*cidia(mgs) + gf1*cidia2(mgs) )      , qi5dt)
      qiacr(mgs) =  &
        min(eri(mgs)*qrw(mgs)* xxacx*cci(mgs)*abs(vtrbar(mgs)-vtibar(mgs))   &
        *(  gf6*rwdia2(mgs) + 2.0*gf5*rwdia(mgs)*cidia(mgs) + gf4*cidia2(mgs) )      , qr5dt) 
!
!  accretions:
!    marshall-palmer size distribution collecting marshall-palmer size
!    distribution
!
      qhacs(mgs) =  &
         min( xxacx*abs(vthbar(mgs)-vtsbar(mgs))*ehs(mgs)*qsw(mgs)*chw(mgs)  &
        *(  gf6*gf1*swdia2(mgs) + 2.0*gf5*gf2*swdia(mgs)*hwdia(mgs) + gf4*gf3*hwdia2(mgs) ) , qs5dt)
      qhacr(mgs) =  &
         min( xxacx*abs(vthbar(mgs)-vtrbar(mgs))*ehr(mgs)*qrw(mgs)*chw(mgs)  &
        *(  gf6*gf1*rwdia2(mgs) + 2.0*gf5*gf2*rwdia(mgs)*hwdia(mgs) + gf4*gf3*hwdia2(mgs) ) , qr5dt)
      qracs(mgs) =  &
         min( xxacx*abs(vtrbar(mgs)-vtsbar(mgs))*ers(mgs)*qsw(mgs)*crw(mgs)  &
        *(  gf6*gf1*swdia2(mgs) + 2.0*gf5*gf2*swdia(mgs)*rwdia(mgs) + gf4*gf3*rwdia2(mgs) ) , qs5dt)
      qsacr(mgs) =  &
         min( xxacx*abs(vtrbar(mgs)-vtsbar(mgs))*ers(mgs)*qrw(mgs)*csw(mgs)  &
        *(  gf6*gf1*rwdia2(mgs) + 2.0*gf5*gf2*rwdia(mgs)*swdia(mgs) + gf4*gf3*swdia2(mgs) ) , qr5dt)
! 
!  bergeron process for snow
!
      ibb = min(max(1,int(-temcg(mgs))),32)
        cs9     = bsfw*dnz(mgs)*(0.001) 
      qsfw(mgs) = qci(mgs)*cs11(ibb)*(cs10(ibb) + eic(mgs)*cs9*qcw(mgs))
      qsfw(mgs) = min(qsfw(mgs),qc5dt)
      qsfi(mgs) = qci(mgs)/cbtim(ibb) 
      qsfi(mgs) = min(qsfi(mgs),qi5dt)
!
!  conversions
!
      qscni(mgs) = 0.001*esicnv(mgs)*max((qci(mgs)-qicrit),0.0)
      qscni(mgs) = min(qscni(mgs),qi5dt)
!
      qhcns(mgs) = 0.001*ehscnv(mgs)*max((qsw(mgs)-qscrit),0.0)
      qhcns(mgs) = min(qhcns(mgs),qs5dt)
!
      xrfrz(mgs) = 20.0*(pi**2)*brz*cwdn/dnz(mgs)
      qrfrz(mgs) = min(xrfrz(mgs)*crw(mgs)*(rwdia(mgs)**6)*(exp(max(-arz*temcg(mgs), 0.0))-1.0), qr5dt)

! Berry (1968) Autoconversion == 0 (critical qc) or ==1 (critical diameter from Ferrier 1994)

      IF( autoconversion .eq. 0 ) THEN
       qdiff  = max((qcw(mgs)-qcmincwrn),0.)
      ELSE
       qccrit = (pi/6.)*((ccw(mgs)*cwdiap**3)*cwdn)/dnz(mgs)
       qdiff  = max((qcw(mgs)-qccrit),0.)
      ENDIF

      qrcnw(mgs) =  0.0
      if ( qdiff .gt. 0.0 ) then
       argrcnw = ((1.2e-4)+(1.596e-12)*ccw(mgs)*(1e-6)/(dnz(mgs)*1.e-3*cwdisp*qdiff))
       qrcnw(mgs) = dnz(mgs)*1e-3*(qdiff**2)/argrcnw
       qrcnw(mgs) = (max(qrcnw(mgs),0.0))
      end if
      qrcnw(mgs) = min(qrcnw(mgs),qc5dt)
!
!  constants for hydrometeor-vapor interactions
!
      ssi(mgs) = qwv(mgs)/qis(mgs)
      ssw(mgs) = qwv(mgs)/qvs(mgs)
      tsqr(mgs) = temg(mgs)**2
!
!  melting of snow and hail
!
      xsv   = (xxmlt1*(swdia(mgs)) + xsmlt2(mgs)*(swdia(mgs)**((3.0+ds)/2.0)))
      xhv   = (xxmlt1*(hwdia(mgs)) + (hwdn**(0.25))*xvth3*xhmlt2(mgs)*(hwdia(mgs)**(1.75)))
      xmlt2 = wvdf(mgs)*elv*dnz(mgs)
      xhsw =(tka(mgs)*temcg(mgs) + xmlt2*(qwv(mgs)-pqs(mgs)))
      xmlt1 = -2.0*pi/(elf*dnz(mgs))
      qsmlr(mgs) = min( (xmlt1*csw(mgs)*xsv*xhsw + temcg(mgs)*xmlt3*(qsacr(mgs)+qsacw(mgs)) ) , 0.0 )
      qhmlr(mgs) = min( (xmlt1*chw(mgs)*xhv*xhsw + temcg(mgs)*xmlt3*(qhacr(mgs)+qhacw(mgs)) ) , 0.0 )
      qsmlr(mgs) = max( qsmlr(mgs), -qs5dt ) 
      qhmlr(mgs) = max( qhmlr(mgs), -qh5dt ) 
!
!  deposition/sublimation of snow and hail
!
      xbs   = (1.0/(dnz(mgs)*wvdf(mgs)))
      xxdsv = 2.0*pi/dnz(mgs)
      xds = xxdsv*(ssi(mgs)-1.0)*(1.0/(xas(mgs)/tsqr(mgs)+xbs/qis(mgs)))
      qsdsv(mgs) =   xds*csw(mgs)*xsv
      qhdsv(mgs) =   xds*chw(mgs)*xhv
      qhsbv(mgs) = max( min(qhdsv(mgs), 0.0), -qh5dt )
      qhdpv(mgs) = max( qhdsv(mgs), 0.0 )
      qssbv(mgs) = max( min(qsdsv(mgs), 0.0), -qs5dt )
      qsdpv(mgs) = max( qsdsv(mgs), 0.0 )

!
! SHEDDING CALCULATION
! New version by MSG closer to JMS original  - Last modified 4/6/03
!       
!
!  compute dry growth rate of hail regardless of location
!
      qhdry(mgs) = qhacr(mgs) + qhacw(mgs) + qhaci(mgs) + qhacs(mgs)
!
!  compute wet growth rate of hail regardless of location
!
      qhacip(mgs)= qhaci(mgs)		!ehi=0 case
      qhacsp(mgs)= qhacs(mgs)		!ehs=0 case
      IF ( ehi(mgs) .gt. 0.0 ) qhacip(mgs) = min(qhaci(mgs)/ehi(mgs),qi5dt)
      IF ( ehs(mgs) .gt. 0.0 ) qhacsp(mgs) = min(qhacs(mgs)/ehs(mgs),qs5dt)

      xcwt = 1.0/( elf +cw*temcg(mgs) )

      xhwet1 = 2.0*pi/dnz(mgs)
      xhwet2 = dnz(mgs)*elv*wvdf(mgs)
      xwt1   = xhwet2*(pqs(mgs)-qwv(mgs)) -tka(mgs)*temcg(mgs)

      qhwet(mgs) =  max( 0.0, ( xhv*chw(mgs)*xwt1*xhwet1*xcwt   &
                                + ( 1.0 -ci(mgs)*temcg(mgs)*xcwt )* ( qhacip(mgs)+qhacsp(mgs) )  )  )

!
!  evaluate shedding rate (effective range is 243 < T < 273 due to other "if" checks below)
!
      qhshr(mgs) = 0.0
      if ( qhwet(mgs) .lt. qhdry(mgs) .and. qhwet(mgs) .gt. 0.0 ) then
        qhdry(mgs) = 0.0                                        ! Wet growth
        qhshr(mgs) = qhwet(mgs) -(qhacw(mgs) +qhacr(mgs))
      else					                ! Dry growth (defaults here if qhwet<0)
        qhwet(mgs) = 0.0
        qhshr(mgs) = 0.0
      endif
!
!  Special shedding case when warmer than freezing
!
      if ( temg(mgs) .gt. tfr ) then 
        qhwet(mgs) = 0.0
        qhdry(mgs) = 0.0
        qhshr(mgs) =  -qhacr(mgs) -qhacw(mgs)-qhacip(mgs)-qhacsp(mgs)
      end if
!
!  Special no-shedding (dry) case when T<243....
!
      if ( temg(mgs) .lt. 243.15 ) then
       qhwet(mgs) = 0.0
       qhshr(mgs) = 0.0
      end if
!
!  Reset some vars if wet particle surface due to shedding....
!
      if ( qhshr(mgs) .lt. 0.0 ) then
       qhaci(mgs) = qhacip(mgs)
       qhacs(mgs) = qhacsp(mgs)
       qhdpv(mgs) = 0.0
       qhsbv(mgs) = 0.0
      end if

!
!  evaporation/condensation on wet snow and hail (NOT USED)
! 

!
!  evaporation of rain
!
      xxcev = 2.0*pi/dnz(mgs)
      xbv   = (1.0/(dnz(mgs)*wvdf(mgs)))
      xce = xxcev*(ssw(mgs)-1.0)*(1.0/(xav(mgs)/tsqr(mgs)+xbv/qvs(mgs)))
      xrv = (xrcev1*(rwdia(mgs)) +  xrcev2(mgs)*(rwdia(mgs)**((3.0+br)/2.0)))
      qrcev(mgs) = max(min(xce*crw(mgs)*xrv, 0.0), -qr5dt)
! 
!  vapor to pristine ice crystals 
!
      qiint(mgs) = 0.0
      IF ( ssi(mgs) .gt. 1.0 ) THEN
        dqisdt(mgs)= (qwv(mgs)-qis(mgs))/ (1.0 + xiint*qis(mgs)/tsqr(mgs))
        cnnt       = cci(mgs)
        qiint(mgs) = (1.0/dtp) *min((1.0e-12)*cnnt/dnz(mgs), 0.50*dqisdt(mgs)) 
      ENDIF
!
!  Domain totals for source terms 
!
!
!  vapor
!
      tvqssbv = tvqssbv - il5(mgs)*qssbv(mgs)*dnz(mgs)   
      tvqhsbv = tvqhsbv - il5(mgs)*qhsbv(mgs)*dnz(mgs)     
      tvqrcev = tvqrcev - qrcev(mgs)*dnz(mgs)             
      tvqcevp = tvqcevp + 0.0	
!     tvqisub = tvqisub + 0.0
!
!  cloud water
!
      tcqccnd = tcqccnd	+ 0.0	
      tcqcmli = tcqcmli	+ 0.0		
!
!  rain
!
      trqracw = trqracw + qracw(mgs)*dnz(mgs)           
      trqrcnw = trqrcnw + qrcnw(mgs)*dnz(mgs)            
      trqsacw = trqsacw + (1-il5(mgs))*qsacw(mgs)*dnz(mgs)
      trqsmlr = trqsmlr - (1-il5(mgs))*qsmlr(mgs)*dnz(mgs) 
      trqhmlr = trqhmlr - (1-il5(mgs))*qhmlr(mgs)*dnz(mgs)
      trqhshr = trqhshr - qhshr(mgs)*dnz(mgs)             
!
!  cloud ice
!
      tiqiint = tiqiint + il5(mgs)*qiint(mgs)*dnz(mgs) 	
      tiqidep = tiqidep + 0.0	
      tiqifzc = tiqifzc + 0.0		
!
!  snow
!
      tsqsacw = tsqsacw + il5(mgs)*qsacw(mgs)*dnz(mgs)    
      tsqscni = tsqscni + il5(mgs)*qscni(mgs)*dnz(mgs)     
      tsqsaci = tsqsaci + il5(mgs)*qsaci(mgs)*dnz(mgs)       
      tsqsfi  = tsqsfi  + il5(mgs)*qsfi(mgs) *dnz(mgs)      
      tsqsfw  = tsqsfw  + il5(mgs)*qsfw(mgs) *dnz(mgs)        
      tsqraci = tsqraci + il5(mgs)*il3(mgs)*qraci(mgs)*dnz(mgs)
      tsqiacr = tsqiacr + il5(mgs)*il3(mgs)*qiacr(mgs)*dnz(mgs)
      tsqsacr = tsqsacr + il5(mgs)*il2(mgs)*qsacr(mgs)*dnz(mgs)
      tsqsdpv = tsqsdpv + il5(mgs)*qsdpv(mgs)*dnz(mgs)         
!
!  hail/graupel
!
      thqhcns = thqhcns + qhcns(mgs)*dnz(mgs)			
      thqiacr = thqiacr + il5(mgs)*(1-il3(mgs))*qiacr(mgs)*dnz(mgs) 
      thqraci = thqraci + il5(mgs)*(1-il3(mgs))*qraci(mgs)*dnz(mgs) 
      thqracs = thqracs + il5(mgs)*(1-il2(mgs))*qracs(mgs)*dnz(mgs) 
      thqsacr = thqsacr + il5(mgs)*(1-il2(mgs))*qsacr(mgs)*dnz(mgs) 
      thqhdpv = thqhdpv + il5(mgs)*qhdpv(mgs)*dnz(mgs)             
      thqrfrz = thqrfrz + il5(mgs)*qrfrz(mgs)*dnz(mgs)   
!     thqhacr = thqhacr + qhacr(mgs)*dnz(mgs)           ! MSG see below instead 2/6/04    
!     thqhacw = thqhacw + qhacw(mgs)*dnz(mgs)           ! MSG see below instead 2/6/04   
!     thqhacs = thqhacs + qhacs(mgs)*dnz(mgs)           ! MSG see below instead 2/6/04
!     thqhaci = thqhaci + qhaci(mgs)*dnz(mgs)           ! MSG see below instead 2/6/04

!
!  hail/graupel and rain (based upon wet growth budget)  !MSG added on 2/6/04
!
!--     
      if ( temg(mgs) .ge. 273.15 ) then
       trqhaci = trqhaci + qhaci(mgs)*dnz(mgs)
       trqhacs = trqhacs + qhacs(mgs)*dnz(mgs)
      else
       thqhaci = thqhaci + qhaci(mgs)*dnz(mgs)
       thqhacs = thqhacs + qhacs(mgs)*dnz(mgs)
      endif

      if ((qhwet(mgs) .gt. 0.0).or.( temg(mgs) .ge. 273.15)) then 
       trqhacw = trqhacw + qhacw(mgs)*dnz(mgs)
       trqhacr = trqhacr + qhacr(mgs)*dnz(mgs)
      else
       thqhacw = thqhacw + qhacw(mgs)*dnz(mgs)
       thqhacr = thqhacr + qhacr(mgs)*dnz(mgs)
      endif

      thqhwet  = thqhwet + qhwet(mgs)*dnz(mgs)    !MSG qhwet is positive or zero here

!--     
!  end of totals
!
 5001 continue
!
 5002 continue
!
!
      if (ndebug .eq. 1 ) print*,'dbg = 8'

 
!  rain, snow, hail fluxes due to gravity    !MSG q3d() was ad() in SAM version
!
!  061013: Moved fallout code to CM1 solve.  Store vtbar info in vq array.
!
      do 5020 mgs = 1,ngscnt 
!!!      hwflx(mgs) = dnz(mgs)*q3d(igs(mgs),jgs,kgs(mgs),lh)*vthbar(mgs)
!!!      piflx(mgs) = dnz(mgs)*q3d(igs(mgs),jgs,kgs(mgs),li)*vtibar(mgs)
!!!      cwflx(mgs) = dnz(mgs)*q3d(igs(mgs),jgs,kgs(mgs),lc)*vtwbar(mgs)
!!!      rwflx(mgs) = dnz(mgs)*q3d(igs(mgs),jgs,kgs(mgs),lr)*vtrbar(mgs)
!!!      swflx(mgs) = dnz(mgs)*q3d(igs(mgs),jgs,kgs(mgs),ls)*vtsbar(mgs)
      vq(igs(mgs),jgs,kgs(mgs),lh) = vthbar(mgs)
      vq(igs(mgs),jgs,kgs(mgs),li) = vtibar(mgs)
      vq(igs(mgs),jgs,kgs(mgs),lc) = vtwbar(mgs)
      vq(igs(mgs),jgs,kgs(mgs),lr) = vtrbar(mgs)
      vq(igs(mgs),jgs,kgs(mgs),ls) = vtsbar(mgs)
 5020 continue

!
!  Compute total-domain content (g/m^3) before production rates
!
      do mgs = 1,ngscnt                        !MSG domain total of each species before microphysics
      tqvbp = tqvbp + qwvp(mgs)*dnz(mgs)
      tqcbp = tqcbp + qcw(mgs)*dnz(mgs)
      tqibp = tqibp + qci(mgs)*dnz(mgs)
      tqrbp = tqrbp + qrw(mgs)*dnz(mgs)
      tqsbp = tqsbp + qsw(mgs)*dnz(mgs)
      tqhbp = tqhbp + qhw(mgs)*dnz(mgs)
      sumb = sumb+tqvbp+tqcbp+tqibp+tqrbp+tqsbp+tqhbp
      end do

!
! CALCULATE RATE TOTALS
!
      do 9000 mgs = 1,ngscnt

!
      pqwvi(mgs) =  il5(mgs)*( -qhsbv(mgs) -qssbv(mgs)              )           - qrcev(mgs)
      pqwvd(mgs) =  il5(mgs)*( -qhdpv(mgs) -qsdpv(mgs) - qiint(mgs) ) 
!
      pqcii(mgs) =  il5(mgs)*qiint(mgs) 
      pqcid(mgs) =  il5(mgs)*( -qscni(mgs) -qsaci(mgs) -qraci(mgs) -qsfi(mgs))  - qhaci(mgs)  
!
      pqcwi(mgs) =  0.0
      pqcwd(mgs) =  (-il5(mgs)*qsfw(mgs))  -qracw(mgs) -qsacw(mgs) -qrcnw(mgs) -qhacw(mgs)
!
      pqrwi(mgs) =  qracw(mgs) +qrcnw(mgs) +(1-il5(mgs))*(qsacw(mgs)-qhmlr(mgs) -qsmlr(mgs)) -qhshr(mgs)
      pqrwd(mgs) =  il5(mgs)*(-qiacr(mgs) -qrfrz(mgs)  -qsacr(mgs)) +qrcev(mgs) -qhacr(mgs)
!
      pqswi(mgs) =  il5(mgs)*( qsacw(mgs) +qscni(mgs) +qsaci(mgs) + qsfi(mgs) +qsfw(mgs)    &
                               +il3(mgs)*(qraci(mgs) +qiacr(mgs)) +il2(mgs)*qsacr(mgs) +qsdpv(mgs)     )
      pqswd(mgs) = -qhcns(mgs) -qhacs(mgs) +(1-il5(mgs))*qsmlr(mgs)  +il5(mgs)*(qssbv(mgs) -(1-il2(mgs))*qracs(mgs)) 
!
      pqhwi(mgs) =  qhcns(mgs) +qhacr(mgs) +qhacw(mgs) +qhacs(mgs) +qhaci(mgs)   &
         + il5(mgs)*( (1-il3(mgs))*(qraci(mgs)+qiacr(mgs)) +(1-il2(mgs))*(qsacr(mgs)+qracs(mgs))+qhdpv(mgs)+qrfrz(mgs))
      pqhwd(mgs) =  qhshr(mgs) +(1-il5(mgs))*qhmlr(mgs) + il5(mgs)*qhsbv(mgs)
!
      ptotal(mgs) = pqwvi(mgs) +pqwvd(mgs) + pqcwi(mgs) +pqcwd(mgs) + pqcii(mgs) +pqcid(mgs) +   &
                    pqrwi(mgs) +pqrwd(mgs) + pqswi(mgs) +pqswd(mgs) + pqhwi(mgs) +pqhwd(mgs) 
!
      if (ndebug .eq. 1) then
       if(abs(ptotal(mgs)).gt.1.e-7)then
         print*,'NOTICE:PTOTAL>1e-7 ', mgs, kgs(mgs), ptotal(mgs)
       end if
      endif
!
      psum = psum + ptotal(mgs)
!
 9000 continue
!
!
!  latent heating from phase changes (except qcw, qci cond, and evap)
!   (22 processes involve phase changes, 10 do not)
!
      do 9010 mgs = 1,ngscnt
      pfrz(mgs) = (1.-il5(mgs))*(qhmlr(mgs) + qsmlr(mgs))   &
                + ( il5(mgs)  )*(qiacr(mgs)+qsacr(mgs)+ qsfw(mgs)+qrfrz(mgs)+qsacw(mgs)+qhacw(mgs)+qhacr(mgs)+qhshr(mgs))
      psub(mgs) = ( il5(mgs)  )*(qhdpv(mgs)+qhsbv(mgs)+qiint(mgs)+qsdpv(mgs)+qssbv(mgs))
      pvap(mgs) = qrcev(mgs)
      ptem(mgs) = cc3(mgs)*pfrz(mgs) + cc5(mgs)*psub(mgs) + cc4(mgs)*pvap(mgs)
      thetap(mgs) = thetap(mgs) + dtp*ptem(mgs)

!
!  partitioned domain-total heating and cooling rates 
!  (all are adjusted again later within saturation adjustment)
!

      hfrz  = hfrz  + ( qiacr(mgs)+qsacr(mgs)+qsfw(mgs)+qsacw(mgs)   &
                       +qhacw(mgs)+qhacr(mgs)+qhshr(mgs)+qrfrz(mgs))*il5(mgs)*cc3(mgs)
      hdep  = hdep  + (qiint(mgs)+qhdpv(mgs)+qsdpv(mgs))*il5(mgs)*cc5(mgs)
      hcnd  = hcnd  + 0.0
      cevap = cevap +  qrcev(mgs)                          *cc4(mgs)
      cmelt = cmelt + (qhmlr(mgs)+qsmlr(mgs))*(1.-il5(mgs))*cc3(mgs)
      csub  = csub  + (qhsbv(mgs)+qssbv(mgs))*il5(mgs)     *cc5(mgs)

 9010 continue
 9004 continue
!
!  sum the sources and sinks for qwvp, qcw, qci, qrw, qsw
!
      do 9100 mgs = 1,ngscnt
      qwvp(mgs)= qwvp(mgs) + dtp*(pqwvi(mgs)+pqwvd(mgs))      !initial qwvp is being adjusted by all source/sink
      qcw(mgs) = qcw(mgs) +  dtp*(pqcwi(mgs)+pqcwd(mgs)) 
      qci(mgs) = qci(mgs) +  dtp*(pqcii(mgs)+pqcid(mgs)) 
      qrw(mgs) = qrw(mgs) +  dtp*(pqrwi(mgs)+pqrwd(mgs)) 
      qsw(mgs) = qsw(mgs) +  dtp*(pqswi(mgs)+pqswd(mgs)) 
      qhw(mgs) = qhw(mgs) +  dtp*(pqhwi(mgs)+pqhwd(mgs)) 
      
 9100 continue
!
!  domain-total content (g/m^3) before saturation adjustment
!
      do mgs = 1,ngscnt
      tqvap = tqvap + qwvp(mgs)*dnz(mgs)
      tqcap = tqcap + qcw(mgs)*dnz(mgs)
      tqiap = tqiap + qci(mgs)*dnz(mgs)
      tqrap = tqrap + qrw(mgs)*dnz(mgs)
      tqsap = tqsap + qsw(mgs)*dnz(mgs)
      tqhap = tqhap + qhw(mgs)*dnz(mgs)
      suma = suma+tqvap+tqcap+tqiap+tqrap+tqsap+tqhap
      end do


! 
      if (ndebug .eq. 1 ) print*,'dbg = 10a'
      
!
!  set up temperature and vapor arrays
!
      if ( ndebug .eq. 1 ) print*,'dbg = 10.1'
!
      do mgs = 1,ngscnt
       pqs(mgs) = (380.0)/(pres(mgs))
       theta(mgs) = thetap(mgs) + theta0(mgs)
       qvap(mgs) = max( (qwvp(mgs) + qv0n(mgs)), 0.0 )         !MSG Current total qwv
       temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp
      end do
!
!  melting of cloud ice
!
      if ( ndebug .eq. 1 ) print*,'dbg = 10.2'
!
      do mgs = 1,ngscnt
      if( temg(mgs) .gt. tfr .and. qci(mgs) .gt. 0.0 ) then
        qimlw(mgs) = - qci(mgs)/dtp                         !MSG Rate of cloudice melting
        tcqcmli = tcqcmli - qimlw(mgs)*dnz(mgs)        !MSG updated 11/22/03 (domain-total rate of qc increase)
        thetap(mgs) = thetap(mgs) - cc3(mgs)*qci(mgs)  !MSG heat decrease
        cmelt   = cmelt   + cc3(mgs)*qimlw(mgs)        !MGS cooling rate
        qcw(mgs) = qcw(mgs) + qci(mgs)
        qci(mgs) = 0.0
      end if
      end do
!
!
!  homogeneous freezing of cloud water
!
      if ( ndebug .eq. 1 ) print*,'dbg = 10.3'
!
      do mgs = 1,ngscnt
      if( temg(mgs) .lt. thnuc .and. qcw(mgs) .gt. 0.0 ) then
        qwfzi(mgs)  = -qcw(mgs)/dtp                         ! MSG Rate of clouwater freezing
        tiqifzc = tiqifzc - qwfzi(mgs)*dnz(mgs)        ! MSG updated 11/22/03 (domain-total rate of qi increase)
        thetap(mgs) = thetap(mgs) + cc3(mgs)*qcw(mgs)  ! MSG heat increase
        hfrz    = hfrz    + cc3(mgs)*(-qwfzi(mgs))     ! MSG heating rate
        qci(mgs) = qci(mgs) + qcw(mgs)
        qcw(mgs) = 0.0
      end if
      end do

!
!  Saturation adjustment iteration procedure
!
!  Modified Straka adjustment (nearly identical to Tao et al. 1989 MWR)
!
!

!
!  reset temporaries for cloud particles and vapor
!

      if ( ndebug .eq. 1 ) print*,'dbg = 10.4'
      do mgs = 1,ngscnt
       ptotsat(mgs) = 0.0
       qwv(mgs) = max( 0.0, qvap(mgs) )
       qcw(mgs) = max( 0.0, qcw(mgs) )
       qci(mgs) = max( 0.0, qci(mgs) )
       ptotsat(mgs) = qwv(mgs)+qci(mgs)+qcw(mgs)       !MSG updated just before sat adj. (qwv+qci+qcw)
       qcevpcnd(mgs) = 0.0
       qisubdep(mgs) = 0.0
      end do
!
      tqvb = tqvap              !MSG domain-total vapor perturb. prior to sat adj.
      do mgs = 1,ngscnt         !MSG domain-total qcw and qci prior to sat adj.
       tqcb = tqcb + qcw(mgs)*dnz(mgs)
       tqib = tqib + qci(mgs)*dnz(mgs)
      enddo

      do mgs = 1,ngscnt
       theta(mgs) = thetap(mgs) + theta0(mgs)
       temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp
       temcg(mgs) = temg(mgs) - tfr
       ltemq = nint((temg(mgs)-163.15)/fqsat+1.5)
       ltemq = min(max(ltemq,1),nqsat) 
       qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       qis(mgs) = pqs(mgs)*tabqis(ltemq)
       if ( temg(mgs) .lt. tfr ) then
        if( qcw(mgs) .ge. 0.0 .and. qci(mgs) .eq. 0.0 ) qss(mgs) = qvs(mgs)
        if( qcw(mgs) .eq. 0.0 .and. qci(mgs) .gt. 0.0)  qss(mgs) = qis(mgs)
        if( qcw(mgs) .gt. 0.0 .and. qci(mgs) .gt. 0.0)  qss(mgs) = (qcw(mgs)*qvs(mgs) + qci(mgs)*qis(mgs))   &	
                                                                 / (qcw(mgs) + qci(mgs))
       else
        qss(mgs) = qvs(mgs)
       end if
      end do
!
!  iterate  adjustment
!
      if ( ndebug .eq. 1 ) print*,'dbg = 10.5'
      do itertd = 1,2
!
      do mgs = 1,ngscnt
!
!  calculate super-saturation
!
      dqcw(mgs) = 0.0
      dqci(mgs) = 0.0
      dqwv(mgs) = ( qwv(mgs) - qss(mgs) )
!
!  evaporation and sublimation adjustment
!
      if( dqwv(mgs) .lt. 0. ) then
       if( qcw(mgs) .gt. -dqwv(mgs) ) then		!Evap some of qc
         dqcw(mgs) = dqwv(mgs)
         dqwv(mgs) = 0.
       else						!Evap all of qc
         dqcw(mgs) = -qcw(mgs)
         dqwv(mgs) = dqwv(mgs) + qcw(mgs)
       end if
!
       if( qci(mgs) .gt. -dqwv(mgs) ) then		!Sublimate some of qi
         dqci(mgs) = dqwv(mgs)
         dqwv(mgs) = 0.
       else						!Sublimate all of qi
         dqci(mgs) = -qci(mgs)
         dqwv(mgs) = dqwv(mgs) + qci(mgs)
       end if
!
       qwvp(mgs) = qwvp(mgs) - ( dqcw(mgs) + dqci(mgs) )	!Increase vapor
      
       qcw(mgs) = qcw(mgs) + dqcw(mgs)			!Decrease cloudwater (dqcw<0)
       qci(mgs) = qci(mgs) + dqci(mgs) 			!Decrease cloudice   (dqci<0)
       thetap(mgs) = thetap(mgs) + cpi/pi0n(mgs)*(elv*dqcw(mgs) +els*dqci(mgs))
       qcevpcnd(mgs) = qcevpcnd(mgs) + (dqcw(mgs)/dtp)
       qisubdep(mgs) = qisubdep(mgs) + (dqci(mgs)/dtp)
      end if
!
! condensation/deposition
!
      if( dqwv(mgs) .ge. 0. ) then
!
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
      gamss = (elv*fracl(mgs) + els*fraci(mgs))/ (pi0n(mgs)*cp)
!
      if ( temg(mgs) .lt. tfr ) then
       if (qcw(mgs) .ge. 0.0 .and. qci(mgs) .le. 0.0 ) then
         dqvcnd(mgs) = dqwv(mgs)/(1. + cqv1*qss(mgs)/((temg(mgs)-cbw)**2))
       end if
       if( qcw(mgs) .eq. 0.0 .and. qci(mgs) .gt. 0.0 ) then
         dqvcnd(mgs) = dqwv(mgs)/(1. + cqv2*qss(mgs)/((temg(mgs)-cbi)**2))
       end if
       if ( qcw(mgs) .gt. 0.0 .and. qci(mgs) .gt. 0.0 ) then
        cdw = caw*pi0n(mgs)*tfrcbw/((temg(mgs)-cbw)**2)
        cdi = cai*pi0n(mgs)*tfrcbi/((temg(mgs)-cbi)**2)
        denom1 = qcw(mgs) + qci(mgs)
        denom2 = 1.0 + gamss*(qcw(mgs)*qvs(mgs)*cdw + qci(mgs)*qis(mgs)*cdi) / denom1
        dqvcnd(mgs) =  dqwv(mgs) / denom2
       end if
      end if

      if ( temg(mgs) .ge. tfr ) then
        dqvcnd(mgs) = dqwv(mgs)/(1. + cqv1*qss(mgs)/ ((temg(mgs)-cbw)**2))
      end if
!
      dqcw(mgs) = dqvcnd(mgs)*fracl(mgs)
      dqci(mgs) = dqvcnd(mgs)*fraci(mgs)
!
      thetap(mgs) = thetap(mgs) + (elv*dqcw(mgs) + els*dqci(mgs))/ (pi0n(mgs)*cp)
      qwvp(mgs) = qwvp(mgs) - ( dqvcnd(mgs) )		!Decrease vapor
      qcw(mgs) = qcw(mgs) + dqcw(mgs)			!Increase cloudwater (dqcw>0)
      qci(mgs) = qci(mgs) + dqci(mgs)			!Increase cloudice   (dqci>0)

      qcevpcnd(mgs) = qcevpcnd(mgs) + (dqcw(mgs)/dtp)
      qisubdep(mgs) = qisubdep(mgs) + (dqci(mgs)/dtp)
      
!
      end if
      end do
!
      if ( ndebug .eq. 1 ) print*,'dbg = 10.51'
      do mgs = 1,ngscnt
       theta(mgs) = thetap(mgs) + theta0(mgs)
       temg(mgs) = theta(mgs)*( pres(mgs) / poo ) ** rcp
       qvap(mgs) =max((qwvp(mgs) + qv0n(mgs)), 0.0)
       temcg(mgs) = temg(mgs) - tfr

       ltemq = nint((temg(mgs)-163.15)/fqsat+1.5)
       ltemq = min(max(ltemq,1),nqsat)
       qvs(mgs) = pqs(mgs)*tabqvs(ltemq)
       qis(mgs) = pqs(mgs)*tabqis(ltemq)
       qcw(mgs) = max( 0.0, qcw(mgs) )
       qwv(mgs) = max( 0.0, qvap(mgs))
       qci(mgs) = max( 0.0, qci(mgs) )      !MSG 
      
       if ( temg(mgs) .lt. tfr ) then
        if( qcw(mgs) .ge. 0.0 .and. qci(mgs) .eq. 0.0 ) qss(mgs) = qvs(mgs)
        if( qcw(mgs) .eq. 0.0 .and. qci(mgs) .gt. 0.0)  qss(mgs) = qis(mgs)
        if( qcw(mgs) .gt. 0.0 .and. qci(mgs) .gt. 0.0)  qss(mgs) = (qcw(mgs)*qvs(mgs) + qci(mgs)*qis(mgs))   &
                                                               / (qcw(mgs) + qci(mgs))
       else
        qss(mgs) = qvs(mgs)
       end if
      end do
      if ( ndebug .eq. 1 ) print*,'dbg = 10.52'
!
!  end the saturation adjustment iteration loop
!
      end do
      if ( ndebug .eq. 1 ) print*,'dbg = 10.6'

      do mgs = 1,ngscnt                     !MSG net at each gpt after all sat adj. iterations are finished
       qcevp(mgs) = min(qcevpcnd(mgs),0.)  ! qcevp <=0
       qisub(mgs) = min(qisubdep(mgs),0.)  ! qisub <=0
       qccnd(mgs) = max(qcevpcnd(mgs),0.)  ! qccnd >=0
       qidep(mgs) = max(qisubdep(mgs),0.)  ! qidep >=0
      end do

      do mgs = 1,ngscnt
       tcqccnd = tcqccnd + qccnd(mgs)*dnz(mgs)    ! MSG updated 2/12/05 domain-total condensation
       tvqcevp = tvqcevp - qcevp(mgs)*dnz(mgs)    ! MSG updated 2/12/05 domain-total qc evaporation
       tvqisub = tvqisub - qisub(mgs)*dnz(mgs)    ! MSG updated 2/12/05 domain-total sublimation
       tiqidep = tiqidep + qidep(mgs)*dnz(mgs)    ! MSG updated 2/12/05 domain-total deposition

          hcnd = hcnd    + qccnd(mgs)*cc4(mgs)    ! MSG Update domain-total heating rate via qc condensation
         cevap = cevap   + qcevp(mgs)*cc4(mgs)    ! MSG Update domain-total cooling rate via qc evap
          csub = csub    + qisub(mgs)*cc5(mgs)    ! MSG Update domain-total cooling rate via qi sublim
          hdep = hdep    + qidep(mgs)*cc5(mgs)    ! MSG Update domain-total heating rate via qi deposition
      end do

!
!  Compute vapor, ice, and cloud totals after saturation adjustment.  
!            
      if (ndebug .eq. 1 ) then
       do mgs = 1,ngscnt
        if(abs(ptotsat(mgs)-qwv(mgs)-qci(mgs)-qcw(mgs)).gt.1.e-7)then
         print*,'NOTICE:PTOTSAT>1e-7 ', mgs, kgs(mgs), ptotsat(mgs),qwv(mgs),qci(mgs),qcw(mgs)
        end if
       end do
      end if

      if (ndebug .eq. 1 ) print*,'dbg = 10b'
! 
      do mgs = 1,ngscnt
        tqva = tqva + qwvp(mgs)*dnz(mgs)
        tqca = tqca +  qcw(mgs)*dnz(mgs)
        tqia = tqia +  qci(mgs)*dnz(mgs)
      end do

      tqv=tqva-tqvb          ! Change in vapor due to sat adj. (Should equal +tvqcevp+tvqisub-tiqidep-tcqccnd  )
      tqc=tqca-tqcb          ! Change in cloud due to sat adj. (Should equal +tcqccnd+tcqcmli-tvqcevp-tiqifzc  )
      tqi=tqia-tqib          ! Change in ice due to sat adj.   (Should equal +tiqidep+tiqifzc-tvqisub-tcqcmli  )
!
!
!
!  end of saturation adjustment
!
!  scatter precipitation fluxes, and thetap, and hydrometeors
!
!DIR$ IVDEP
      do 4001 mgs = 1,ngscnt
      ! GHB, 061013:  fallout calculations have been moved to CM1 solve
!!!      pflux(igs(mgs),kgs(mgs)) = piflx(mgs)
!!!      cflux(igs(mgs),kgs(mgs)) = cwflx(mgs)
!!!      rflux(igs(mgs),kgs(mgs)) = rwflx(mgs)
!!!      sflux(igs(mgs),kgs(mgs)) = swflx(mgs)
!!!      hflux(igs(mgs),kgs(mgs)) = hwflx(mgs)
      th3d(igs(mgs),jy,kgs(mgs))  =  thetap(mgs) 
      q3d(igs(mgs),jy,kgs(mgs),lv) = ab(kgs(mgs),2) +   qwvp(mgs) 
      q3d(igs(mgs),jy,kgs(mgs),lc) = qcw(mgs)  + min( q3d(igs(mgs),jy,kgs(mgs),lc), 0.0 ) !MSG putting any neg Gibbs values back
      q3d(igs(mgs),jy,kgs(mgs),li) = qci(mgs)  + min( q3d(igs(mgs),jy,kgs(mgs),li), 0.0 )  
      q3d(igs(mgs),jy,kgs(mgs),lr) = qrw(mgs)  + min( q3d(igs(mgs),jy,kgs(mgs),lr), 0.0 )  
      q3d(igs(mgs),jy,kgs(mgs),ls) = qsw(mgs)  + min( q3d(igs(mgs),jy,kgs(mgs),ls), 0.0 )  
      q3d(igs(mgs),jy,kgs(mgs),lh) = qhw(mgs)  + min( q3d(igs(mgs),jy,kgs(mgs),lh), 0.0 )  
 4001 continue
!
!
 9998 continue
!
!---------------------------------
! new, 071008, GHB (from MSG):
      if ( (ix+1 .gt. nx-istag) .and. &   ! x&z at end. Exit.       !MSG Corrected logic on 5 Oct 2007
           (kz+1 .gt. nz-kstag-1) ) then
        go to 1200
      else if (ix+1 .gt. nx-istag ) then  ! x at end. Reset x & inc z
        nzmpb = kz+1 
        nxmpb = 1
      else                                ! none at the end.  Inc x.
        nzmpb = kz 
        nxmpb = ix+1
      end if
!---------------------------------
! 061023, GHB:  dunno why the "kstag-1" is here ... removing the -1
!!!      if ( kz .gt. nz-kstag-1 .and. ix .gt. nx-istag ) then
!  old code:
!      if ( kz .gt. nz-kstag .and. ix .gt. nx-istag ) then
!        go to 1200
!      else
!        nzmpb = kz 
!      end if
!      
!      if ( ix+1 .gt. nx ) then
!        nxmpb = 1
!      else
!        nxmpb = ix+1
!      end if
!---------------------------------
      
 1000 continue      !MSG  end of numgs loop
 1200 continue
!
!  end of gather scatter
!
!  precipitation fallout contributions 
!  MSG Technically we probably should really be re-computing the fluxes based upon new mixing ratio's but we don't.
!
      ! GHB, 061013:  fallout is now computed in CM1 solve.

!      if (itfall .gt. 0) then
!       dtsplit = dtp
!       nrnstp  = 1
!
!       if (itfall .eq. 2) then  ! COMPUTE SPLIT-EXPLICIT TIME STEP USING MAX FALLSPEED OF ANY SPECIES IN THE SLICE
!        vtrdzmax = dtp*dzc(kz)*maxfall   ! Courant Number
!        nrnstp = 1 + aint(2.*vtrdzmax)   ! MSG Define number of small steps based upon
!                                         ! fastest falling particles (rain or hail/graupel)
!                                         ! nrnstp for 1 + aint(2*0.45) = 1 (use original timestep)
!                                         ! whereas 1+aint(2*0.55) = 2 (requiring dtp be split into
!                                         ! two smaller timesteps.) Thus if courant number, vtrdzmax>=0.5,
!                                         ! then automatically uses 2 small steps (0.25 courant each).
!                                         ! If vtrdzmax=1.0, then uses 3 small steps (0.33 courant each).
!                                         ! If vtrdzmax=1.5, then uses 4 small steps (0.25 courant, Etc.)
!                                         ! Thus the actual fallout courant ranges 0.25 to 0.5
!        dtsplit = dtp / float( nrnstp )  ! MSG Define small timestep (sec)
!       endif
!
!       DO inrnstp = 1,nrnstp            ! MSG BEGIN OF TIMESPLITTING
!        if (ndebug .eq. 1 ) print*,'dbg = 10g'
! 061023, GHB:  dunno why the "kstag-1" is here ... removing the -1
!!!!        do kz = 1,nz-kstag-1 
!        do kz = 1,nz-kstag
!        do ix = 1,nx-istag
!!         dtz1 = dzc(kz)*dtp/1.0        !MSG old version before timesplitting
!          dtz1 = dzc(kz)*dtsplit
!          q3d(ix,jy,kz,li) = q3d(ix,jy,kz,li) + dtz1*(pflux(ix,kz+1)-pflux(ix,kz))/db(kz)
!          q3d(ix,jy,kz,lc) = q3d(ix,jy,kz,lc) + dtz1*(cflux(ix,kz+1)-cflux(ix,kz))/db(kz)
!          q3d(ix,jy,kz,lr) = q3d(ix,jy,kz,lr) + dtz1*(rflux(ix,kz+1)-rflux(ix,kz))/db(kz)
!          q3d(ix,jy,kz,ls) = q3d(ix,jy,kz,ls) + dtz1*(sflux(ix,kz+1)-sflux(ix,kz))/db(kz)
!          q3d(ix,jy,kz,lh) = q3d(ix,jy,kz,lh) + dtz1*(hflux(ix,kz+1)-hflux(ix,kz))/db(kz)
!        enddo
!        enddo
!       ENDDO      !MSG End of timesplitting on fallout
!
!      endif     !end of itfall check


!
!
!  end of jy loop
!
 9999 continue
!
!
      if (ndebug .eq. 1 ) print*,'dbg = 10h'
!
!
!  WRITE totals for source / sink terms
      if (nrates .eq. 1 ) then
!
        tvsum  = tvqrcev -tiqiint -thqhdpv -tsqsdpv +tvqhsbv +tvqssbv
        tcsum  = -tsqsfw -thqhacw -trqrcnw -trqracw -tsqsacw -trqsacw
        trsum  = -tvqrcev -thqhacr +trqrcnw -thqrfrz +trqhshr +trqhmlr +trqsmlr -tsqiacr -thqiacr -tsqsacr -thqsacr   &
                 +trqracw +trqsacw
        tisum  = -tsqsfi -thqhaci -tsqscni +tiqiint -tsqraci -thqraci  -tsqsaci
        tssum  = tsqsfw +tsqsfi -thqhacs +tsqscni -thqhcns -trqsmlr +tsqiacr +tsqraci +tsqsacr +tsqsdpv -tvqssbv      &
                +tsqsaci +tsqsacw -thqracs
        thsum  = thqhacw +thqhacr +thqhacs +thqhaci +thqhcns +thqrfrz -trqhshr -trqhmlr +thqiacr +thqraci +thqsacr +thqhdpv   &
                 +thqracs -tvqhsbv
        tsumall =  tvsum + tcsum + tisum + trsum + tssum + thsum        !MSG gives ~1e-9 (machine precision)
!
        if (ndebug .eq. 1) then    ! MSG only print details if debugging on
!
        write(6,*) 'Sum species source/sink domain totals for all but Sat. Adj. (kg s^-1 m^-3)'
        write(6,*) 'qv', tvsum            ! Total rate for qv (not including sat adj.)
        write(6,*) 'qc', tcsum
        write(6,*) 'qi', tisum
        write(6,*) 'qr', trsum
        write(6,*) 'qs', tssum
        write(6,*) 'qh', thsum
!
        write(6,*) 'Sum only for saturation adjustment (kg s^-1 m^-3)'
        write(6,*) 'qv', tqv              !MSG total change in vapor only due to sat adj.
        write(6,*) 'qc', tqc              !MSG change in cloud only due to sat adj.
        write(6,*) 'qi', tqi              !MSG change in ice only due to sat adj.
!
        write(6,*) 'Sum all species but not sat. adj. (kg s^-1 m^-3)'
        write(6,*) 'sum', tsumall         !MSG - all sources/sink rates except for sat adj.(comes from tsum rates)
!       write(6,*) 'psum', psum           !MSG - all sources/sink rates except for sat adj.(comes from model rates)
!       write(6,*) 'sumchange', suma-sumb !MSG - change in content (g/m^3) over timestep (not incl. sat adj.)
!
        write(6,*) 'Sum of all rates (kg s^-1 m^-3)'
        write(6,*) 'tsa', tsumall+tqv+tqc+tqi  !MSG all source/sink rates including sat adj.
!       write(6,*) 'psa', psum+tqv+tqc+tqi     !MSG all source/sink rates including sat adj.
!
!       write(6,*) 'tqvp',tqvap-tqvbp          !Changes due to all (except sat adj or fallout)
!       write(6,*) 'tqcp',tqcap-tqcbp
!       write(6,*) 'tqip',tqiap-tqibp
!       write(6,*) 'tqrp',tqrap-tqrbp
!       write(6,*) 'tqsp',tqsap-tqsbp
!       write(6,*) 'tqhp',tqhap-tqhbp
!
        endif
!
!
        write(6,*) 'Individual domain total rates (kg s^-1 m^-3)'
        write(6,*) 'Using Gilmore et al. (2004b) terminology'
        write(6,*) 'vapor sources'
!
        write(6,*) 'qvevr ' , tvqrcev 
        write(6,*) 'qvevw ' , tvqcevp 
        write(6,*) 'qvsbi ' , tvqisub 
        write(6,*) 'qvsbs ' , tvqssbv 
        write(6,*) 'qvsbh ' , tvqhsbv 
!
        write(6,*) 'cloud water sources'
!
        write(6,*) 'qwcdv ' , tcqccnd 
        write(6,*) 'qwmli ' , tcqcmli 
!
        write(6,*) 'rain sources'
!
        write(6,*) 'qrhacr' , trqhacr     !Added 2/6/04  (rain accreted and shed during same timestep)
        write(6,*) 'qrmlh ' , trqhmlr 
        write(6,*) 'qracw ' , trqracw 
        write(6,*) 'qrcnw ' , trqrcnw 
!       write(6,*) 'qrshh ' , trqhshr      !MSG ambiguous since sign can switch. Instead, use 4 individual terms (2/6/04)
        write(6,*) 'qrhacw' , trqhacw     !Added 2/6/04  (from shedding)
        write(6,*) 'qrhacs' , trqhacs     !Added 2/6/04  (from shedding)
        write(6,*) 'qrmls ' , trqsmlr 
        write(6,*) 'qrsacw' , trqsacw
        write(6,*) 'qrhaci' , trqhaci     !Added 2/6/04  (from shedding)
!
        write(6,*) 'cloud ice sources'
!
        write(6,*) 'qidpv ' , tiqidep
        write(6,*) 'qiint ' , tiqiint 
        write(6,*) 'qifzw ' , tiqifzc 
!
        write(6,*) 'snow sources'
!
        write(6,*) 'qsfi  ' , tsqsfi 
        write(6,*) 'qsacw ' , tsqsacw
        write(6,*) 'qsaci ' , tsqsaci
        write(6,*) 'qsdpv ' , tsqsdpv
        write(6,*) 'qsiacr' , tsqiacr
        write(6,*) 'qsacr ' , tsqsacr
        write(6,*) 'qsfw  ' , tsqsfw 
        write(6,*) 'qsraci' , tsqraci
        write(6,*) 'qscni ' , tsqscni
!
        write(6,*) 'hail sources'
!
        write(6,*) 'qhacw ' , thqhacw   ! 2/6/04 (only that which is not shed as rain)
        write(6,*) 'qhacr ' , thqhacr   ! 2/6/04 (only that which is not shed as rain)
        write(6,*) 'qhdpv ' , thqhdpv
        write(6,*) 'qhsacr' , thqsacr
        write(6,*) 'qhwtr ' , thqhwet   ! 2/6/04 (only that which is not shed as rain)
        write(6,*) 'qhacs ' , thqhacs   ! 2/6/04 (only that which is not shed as rain)
        write(6,*) 'qhfzr ' , thqrfrz
        write(6,*) 'qhaci ' , thqhaci   ! 2/6/04 (only that which is not shed as rain)
        write(6,*) 'qhracs' , thqracs
        write(6,*) 'qhcns ' , thqhcns
        write(6,*) 'qhiacr' , thqiacr
        write(6,*) 'qhraci' , thqraci
!
        write(6,*) 'total heating and cooling rates'
!
        write(6,*) 'hfrz'    , hfrz
        write(6,*) 'hdep'    , hdep
        write(6,*) 'hcnd'    , hcnd
        write(6,*) 'cmelt'   , cmelt
        write(6,*) 'cevap'   , cevap
        write(6,*) 'csub'    , csub
      endif
! 
!
!
      if (ndebug .eq. 1 ) print*,'dbg = 11a'

      end subroutine lfo_ice
!
!  end of subroutine
!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
! Routine from Numerical Recipes to replace other gamma function
!  using 32-bit reals, this is accurate to 6th decimal place.

!     REAL FUNCTION GAMMA(xx)

!     implicit none
!     real xx
!     integer j

! Double precision ser,stp,tmp,x,y,cof(6)

!     real*8 ser,stp,tmp,x,y,cof(6)
!     SAVE cof,stp
!     DATA cof,stp/76.18009172947146d+0,
!    $            -86.50532032941677d0,
!    $             24.01409824083091d0,
!    $             -1.231739572450155d0,
!    $              0.1208650973866179d-2,
!    $             -0.5395239384953d-5,
!    $              2.5066282746310005d0/

!     x = xx
!     y = x
!     tmp = x + 5.5d0
!     tmp = (x + 0.5d0)*Log(tmp) - tmp
!     ser = 1.000000000190015d0
!     DO j=1,6
!       y = y + 1.0d0
!       ser = ser + cof(j)/y
!     END DO
!     gamma = Exp(tmp + log(stp*ser/x))

!     RETURN
!     END

!--------------------------------------------------------------------------

      subroutine lfoice_init(dtp)
      implicit none

!!!      include 'lfoice.incl'

      real, intent(in) :: dtp

      outfile=6
#ifdef MPI
      outfile=10
#endif

!
!  constants
!
      poo = 1.0e+05
      ar = 841.99666  
      br = 0.8
      bta1 = 0.6
      cnit = 1.0e-02
      dnz00 = 1.225
      rho00 = 1.225
      cs = 4.83607122
      ds = 0.25
      pi = 4.0*atan(1.0)
      pid4 = pi/4.0 
      qccrit = 2.0e-03
      qscrit = 6.0e-04
      qicrit = 1.0e-03
!      
! Define gamma functions
!
      gf1 = 1.0      !gamma(1.0)
      gf2 = 1.0      !gamma(2.0)
      gf3 = 2.0      !gamma(3.0)
      gf4 = 6.0      !gamma(4.0)
      gf5 = 24.0     !gamma(5.0)
      gf6 = 120.0    !gamma(6.0)
      gf4br = 17.837862  !gamma(4.0+br)
      gf4ds = 8.2850851  !gamma(4.0+ds)
      gf5br = 1.8273551  !gamma((5.0+br)/2.)
      gf5ds = 1.4569332  !gamma((5.0+ds)/2.)
      gf2p75= 1.6083594  !gamma(2.75)
      gf4p5 = 11.631728  !gamma(4.0+0.5)
      
      if ( hwdn .lt. 600.0 ) then
      dragh = 1.00
      else
      dragh = 0.60
      end if

      xcnor = cnor
      xcnos = cnos
      xcnoh = cnoh

       write(outfile,*) '-----------------------------------------------------------------------'
       write(outfile,*)
       write(outfile,*) '2004 STRAKA(GILMORE) SAM microphysics'
       write(outfile,*)
       write(outfile,*) '3-ICE MICROPHYSICAL CONSTANTS'
101    Format(1x,a,6(g12.5,2x))
       write(outfile,101) 'CNOR/DENR:        ',cnor, rwdn
       write(outfile,101) 'CNOS/DENS:        ',cnos, swdn
       write(outfile,101) 'CNOH/DENH/DRAGH:  ',cnoh, hwdn,dragh
       write(outfile,*)   'TIME STEP:        ', dtp
         IF( autoconversion .eq. 0 ) THEN
        write(outfile,*) 'Berry (1968) Critical qc g/g for autoconversion= ',qcmincwrn
         ELSE
        write(outfile,*) 'Berry (1968) Critical qc diam for autoconversion= ',cwdiap*1e6,' microns'
         ENDIF
       write(outfile,*) '-----------------------------------------------------------------------'
       write(outfile,*)

!  constants
!
      c1f3 = 1.0/3.0
!
!  general constants for microphysics
!
      brz = 100.0
      arz = 0.66
      cai = 21.87455
      caw = 17.2693882
      cbi = 7.66
      cbw = 35.86
      qcmin = 1.0e-09
      qimin = 1.0e-12
      qrmin = 1.0e-07
      qsmin = 1.0e-07
      qhmin = 1.0e-07

      tfr = 273.15
      thnuc = 233.15
      advisc0 = 1.832e-05
      advisc1 = 1.718e-05
      tka0 = 2.43e-02
      cpi = 1.0/cp
      tfrcbw = tfr - cbw
      tfrcbi = tfr - cbi

      elv = 2500300. 
      elf = 335717. 
      els = elv + elf
      cw = 4218.0   
      xmlt3 = -cw/elf
      xiint = (elv**2)/(cp*rw)

      xacwi = pid4          
      xxacx = pid4/gf4
      xvth1 = gf4p5/(6.0)
      xxmlt1 = 0.78*gf2
      xrcev1 = 0.78*gf2
      xvth3 = (4.0*g/(3.0*dragh))**(0.25)
      xslop = (1./pi)**(0.25)
      xvtr  = ar*gf4br*(dnz00**0.5)/6.0
      xvts  = cs*gf4ds*(dnz00**0.5)/6.0

      cqv1 = 4097.8531*elv*cpi
      cqv2 = 5807.4743*els*cpi
!
!
!  Saturation Vapor Pressure Lookup Table
!
      do lll = 1,nqsat
      temq = 163.15 + (lll-1)*fqsat
      tabqvs(lll) = exp(caw*(temq-273.15)/(temq-cbw))
      tabqis(lll) = exp(cai*(temq-273.15)/(temq-cbi))
      end do
!
!  cw constants in mks units
!
      cwmasn = 4.25e-15
      cwmasx = 5.25e-10
      cwdn = 1000.0
      cwc1 = 6.0/(pi*cwdn)
!
!  ci constants in mks units
!
      cimasx = 3.23e-8
      cimasn = 4.25e-15
!
!  constants for bergeron process, note in cgs units
!
      cmn = 1.05e-15
      cmi40 = 2.4546e-07
      cmi50 = 4.8e-07
      ri50 = 5.0e-03
      vti50 = 100.0
      bsfw = (ri50**2)*vti50*pi
      cbtim(1) = 1.0e+30
      cs10(1)  = 0.0
      cs11(1)  = 0.0
      do ibb = 2,32
        cm50a = cmi50**bfa2(ibb)
        atemp = 1.0-bfa2(ibb)
        cm40b = cmi40**atemp
        cm50b = cmi50**atemp
        cbtim(ibb) = (cm50b-cm40b)/(bfa1(ibb)*atemp)
        cs10(ibb) = bfa1(ibb)*cm50a
        cs11(ibb) = dtp/(cbtim(ibb)*cmi50)
!       write(6,*) 'ibb, cbtim(ibb)=',ibb, cbtim(ibb)
      ENDDO

      end subroutine lfoice_init

!--------------------------------------------------------------------------

  END MODULE lfoice_module
