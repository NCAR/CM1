  MODULE base_module

  implicit none

  private
  public :: base

  CONTAINS

!-----------------------------------------------------------------------------
!
!  subroutine BASE:  Specify base-state arrays for CM1. 
!
!    Note:   The base state in CM1 is assumed to be in hydrostatic balance.
!            See "The governing equations for CM1" for more details:
!               http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!
!    Note:  The variables psurf, tsurf, and qsurf refer to very-near-surface 
!           ATMOSPHERIC conditions (not land/ocean surface conditions).
!           They are used primarily for diagnostics purposes, such as for 
!           calculation of surface-based CAPE and CIN.
!
!-----------------------------------------------------------------------------

      subroutine base(zh,mh,rmh,c1,c2,zf,mf,rho0s,pi0s,prs0s,rth0s,         &
                      wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,               &
                      uavg,vavg,thavg,pavg,qavg,ulspg,vlspg,                &
                      pi0,prs0,rho0,thv0,th0,rth0,qv0,u0,v0,thrd,           &
                      qc0,qi0,rr0,rf0,rrf0,t0  ,rh0 ,                       &
                      reqs_u,reqs_v,reqs_s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2, &
                      n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,              &
                      uw31,uw32,ue31,ue32,us31,us32,un31,un32,              &
                      vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,              &
                      sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
      use input
      use constants
      use bc_module
      use comm_module
      use goddard_module, only : T0K,T00K,RT0
      use cm1libs , only : rslf,rsif
      use getcape_module
#ifdef MPI
      use mpi
#endif
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s

      real, intent(inout), dimension(kb:ke) :: wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,  &
                                               uavg,vavg,thavg,pavg
      real, intent(inout), dimension(kb:ke,numq) :: qavg
      real, intent(inout), dimension(kb:ke) :: ulspg,vlspg
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,rho0,thv0,th0,rth0,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: qc0,qi0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: t0,rh0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ibb2:ibe2,jbb2:jbe2,kbb2:kbe2) :: thrd
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, intent(inout), dimension(cmp,jmp,kmp) :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp) :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp) :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp) :: ss31,ss32,sn31,sn32

!-----------------------------------------------------------------------

      integer i,j,k,m,n,nn,irec,niter,nsnd,kbot,ktop,tflag,nmax
      real zu,zv
      real z_trop,th_trop,th_sfc,t_trop,prs_sfc,qv_pbl,pi_sfc,t_sfc,rh_sfc,rh_pbl
      real qv_sfc,thv_sfc,psurf,tsurf,qsurf,thsurf,thvsurf
      real tlast,ql2,thbar,qvbar,tlcl
      real ns,ns1,ns2,ns3,zl1,zl2,zsfc
      real qcval,qtval,p_sfc,ql_sfc,thn,thlast,qtm,qtp,tavg,lhv,     &
           qvavg,qlavg,qtavg,desdt,gamma,tp,delz,pim,thvm,qvm,qlm,   &
           tm,pm,thm,drdt,prsavg,qvl,qvi,fliq,fice,cpml,qim
      real thex,es,qvs,aaa
      real udep,uconst1,uconst2
      real udep1,udep2,umax1,umax2,vmax1,angle,a1,a2,wsp
      real tmp,ql,t1,t2,pitmp
      real hs,lapse
      real alpha,umax,nm,dudz,dvdz,rinum
      real, dimension(:), allocatable :: zsnd,thsnd,qvsnd,usnd,vsnd,   &
                                         thvsnd,pisnd,psnd,tsnd,rhsnd
      real, dimension(:), allocatable :: thinterp,qvinterp,uinterp,vinterp, &
                                         pinterp,tinterp,rhinterp
      integer :: kk,kup,kdn
      real :: interp_frac

      integer :: flag,ttype
      real :: pisfc,the_sfc,thv1,thv2,pi1,pi2,p2,theq,qt_sfc

      real :: ztrop,zmix,qv_mix,zmin,dtheta,thv_trop,qv_trop,pi_trop,p_trop,rhexp
      real :: depth_layer_01,depth_layer_02,shear_layer_01,shear_layer_02
      real :: du,dv

      real, dimension(:), allocatable :: pfoo,tfoo,qvfoo
      real :: cape,cin,zlcl,zlfc,zel,psource,tsource,thsource,qvsource

      integer :: k1,k2
      real :: lr1,lr2,ths1
      real :: tem
      logical :: dothis
      real :: wsp10,cd,znt,ust,wspr,gam,wsp1,wsp2,wtop,wref,radsfc,qvsfc,thsfc
      real :: uval,vval,vmax,wmin
      real :: z1,z2,z3,z4,z5,z6,z7,z8,z9
      real :: th1,th2,th3,th4,th5,th6,th7,th8,th9
      real :: qv1,qv2,qv3,qv4,qv5,qv6,qv7,qv8,qv9
      real :: pilast,thvlast,zlast,qt,thl,z_inv,th0last,lwp
      real :: qsfc,zq1,zq2,tv0,tsfc,tvsfc,psfc,lr,tvtrop,ptrop
      real :: rhval

!------------------------------------------------------------------

      if(dowr) write(outfile,*) 'Inside BASE'

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  Start definition of base state sounding (isnd opton)  cccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      psurf = 0.0
      tsurf = 0.0
      qsurf = 0.0

!-----------------------------------------------------------------------
!  isnd = 1
!  Dry adiabatic base state

      IF(isnd.eq.1)THEN

        ! Set these two variables for dry adiabatic sounding

        th_sfc   =   300.0   ! Potential temperature of atmosphere (K)
        pi_sfc   =   1.0     ! Exner function at surface

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          th0(i,j,k)=th_sfc
          pi0(i,j,k)=pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo

        thsurf = th_sfc
        psurf = p00*(pi_sfc**cpdrd)
        tsurf = thsurf*((psurf*rp00)**rovcp)
        qsurf = 0.0

!-----------------------------------------------------------------------
!  isnd = 2
!  Dry isothermal base state
            
      ELSEIF(isnd.eq.2)THEN
          
        ! Set these two variables for dry isothermal sounding

        t_sfc    =   250.0    ! Temperature of atmosphere (K)
        prs_sfc  =   p00      ! Pressure at surface (Pa)


        hs=rd*t_sfc/g        ! scale height of atmosphere

        do k=kb,ke
        do j=jb,je
        do i=ib,ie

          ! calculate pressure field
          prs0(i,j,k)=prs_sfc*EXP(-zh(i,j,k)/hs)

          ! using the pressure field, calculate the exner pressure
          pi0(i,j,k)=(prs0(i,j,k)/p00)**(rd/cp)

          ! using exner pressure, determine the potential temperature
          th0(i,j,k)=t_sfc/pi0(i,j,k)

        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo

        psurf = prs_sfc
        tsurf = t_sfc
        qsurf = 0.0


!-----------------------------------------------------------------------
!  isnd = 3
!  Dry, constant dT/dz sounding.
!  Lapse rate of 0.0065  =  standard dry atmosphere.

      ELSEIF(isnd.eq.3)THEN

        ! Set these three variables for dry constant lapse rate sounding

        th_sfc   =  300.0     ! theta at surface (K)
        prs_sfc  =  p00       ! pressure at surface (Pa)
        lapse    =  0.0065    ! dT/dz (K m^-1)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie

          ! Calculate the temperature using the specified lapse rate.
          t0(i,j,k)=th_sfc-lapse*zh(i,j,k)

          ! Calculate the pressure from the temperature field.
          prs0(i,j,k)=p00*(t0(i,j,k)/th_sfc)**(g/(lapse*rd))

          ! Calculate the exner pressure from the pressure field.
          pi0(i,j,k)=(prs0(i,j,k)/p00)**(rd/cp)

          ! Calculate the theta field from temperature and the 
          ! specified lapse rate.
          th0(i,j,k)=th_sfc*(t0(i,j,k)/th_sfc)**(1-g/(lapse*cp))

        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo

        psurf = prs_sfc
        tsurf = th_sfc*((psurf*rp00)**rovcp)
        qsurf = 0.0


!------------------------------------------------------------------
!  isnd = 4
!  Saturated, neutrally-stable sounding for moist benchmark simulation.
!  reference:  Bryan and Fritsch, 2002, MWR, 130, 2917-2928.

      ELSEIF(isnd.eq.4)THEN

        ! these two parameters define the sounding

        thec_mb   =   320.0     ! wet equivalent potential temp (K)
        qt_mb     =   0.020     ! total water mixing ratio (unitless)


        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Saturated, neutrally-stable sounding'
        if(dowr) write(outfile,*) '    thec,qt=',thec_mb,qt_mb
        if(dowr) write(outfile,*)

        prs_sfc = 100000.0

      do j=jb,je
      do i=ib,ie

        !! First guesses at lowest model level
        prs0(i,j,1)=prs_sfc
        qv1=qt_mb
        tmp=thec_mb
        tlast=tmp

        do n=1,20
          thbar=(tmp*(p00/prs0(i,j,1))**(rd/cp))*(1.0+qv1*reps)/(1.0+qt_mb)
          pi0(i,j,1)=1.0-g*zh(i,j,1)/(cp*thbar)
          prs0(i,j,1)=p00*(pi0(i,j,1)**(cpdrd))
          qv1=rslf(prs0(i,j,1),tmp)
          tmp=thec_mb*((prs0(i,j,1)/(1.0+qv1/eps)/p00)**(rd/(cp+cpl*qt_mb)))  &
             /exp((lv1-lv2*tmp)*qv1/(tmp*(cp+cpl*qt_mb)))
          tmp=tlast+0.3*(tmp-tlast)
          tlast=tmp
        enddo
 
        t0(i,j,1)=tmp
        th0(i,j,1)=tmp/pi0(i,j,1)
        rh0(i,j,1)=1.0
        qv0(i,j,1)=rslf(prs0(i,j,1),tmp)
        t1=tmp
        qv1=qv0(i,j,1)
        th1=th0(i,j,1)

        do k=2,nk
          tlast=th1
          t2=t1
          qv2=qv1
          th2=th1
          n=0
          pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*th2)
100       continue
            n=n+1
            th2=tlast
            t2=th2*pi0(i,j,k)
            thbar=0.5*( th1*(1.0+qv1*reps)/(1.0+qt_mb)    &
                       +th2*(1.0+qv2*reps)/(1.0+qt_mb) )
            pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*thbar)
            prs0(i,j,k)=p00*(pi0(i,j,k)**(cpdrd))
            qv2=rslf(prs0(i,j,k),t2)
            t2=thec_mb*((prs0(i,j,k)/(1.0+qv2/eps)/p00)**(rd/(cp+cpl*qt_mb))) &
             /exp((lv1-lv2*t2)*qv2/(t2*(cp+cpl*qt_mb)))
            th2=t2/pi0(i,j,k)

            if(n.gt.50.and.dowr) write(outfile,*) n,th2
            if(abs(th2-tlast).gt.0.0001 .and. n.lt.100)then
              tlast=tlast+0.3*(th2-tlast)
              go to 100
            elseif(n.ge.100)then
              if(dowr) write(outfile,*) '  stuck in loop!'
              call stopcm1
            endif
 
          t0(i,j,k)=t2
          th0(i,j,k)=th2
          qv0(i,j,k)=rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k))
          rh0(i,j,k)=1.0
 
          th1=th2
          qv1=qv2
          t1=t2
 
        enddo

        do k=1,nk
          qc0(i,j,k)=qt_mb-qv0(i,j,k)
        enddo

      enddo
      enddo

        psurf = prs_sfc
        tsurf = cgs1*t0(1,1,1)+cgs2*t0(1,1,2)+cgs3*t0(1,1,3)
        qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)


!-----------------------------------------------------------------------
!  isnd = 5
!  Weisman-Klemp analytic sounding.
!  reference:  Weisman and Klemp, 1982, MWR, 110, 504-520.

      ELSEIF(isnd.eq.5)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  WK sounding'

!  variables related to Weisman-Klemp analytic sounding
        z_trop   = 12000.0      ! height of tropopause (m)
        th_trop  = 343.0        ! theta at tropopause (K)
        t_trop   = 213.0        ! temp at tropopause (K)
        th_sfc   = 300.0        ! theta at surface (K)
        prs_sfc  = 100000.0     ! pressure at surface (Pa)
        qv_pbl   = 0.014        ! constant value of mixing ratio in PBL

!--------------

        pi_sfc  = (prs_sfc/p00)**(rd/cp)
        qv_sfc  = rslf(prs_sfc,th_sfc*pi_sfc)
        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)

      do j=jb,je
      do i=ib,ie

        do k=kb,ke
          rh0(i,j,k)=0.0
        enddo

        do k=1,nk
          if(zh(i,j,k).lt.z_trop)then
            th0(i,j,k)=th_sfc+(th_trop-th_sfc)*((zh(i,j,k)/z_trop)**1.25)
            if(imoist.eq.1) rh0(i,j,k)=1.0-0.75*((zh(i,j,k)/z_trop)**1.25)
          else
            th0(i,j,k)=th_trop*exp((g/(t_trop*cp))*(zh(i,j,k)-z_trop))
            if(imoist.eq.1) rh0(i,j,k)=0.25
          endif
        enddo

        th0(i,j,0)=th0(i,j,1)
        th0(i,j,nk+1)=th0(i,j,nk)


!  Get pressure, temperature, and mixing ratio using hydrostatic eqt.

        do k=kb,ke
          qv0(i,j,k)=0.0
        enddo

        do n=1,20
! virtual potential temperature
          do k=kb,ke
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
          enddo

          pi0(i,j,1)=pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
          do k=2,nk
            pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
          enddo

! pressure
          do k=1,nk
            prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          enddo

! mixing ratio
          do k=1,nk
            qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k))
            if(qv0(i,j,k).gt.qv_pbl) qv0(i,j,k)=qv_pbl
          enddo

        enddo

        do k=1,nk
          rh0(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
        enddo

      enddo
      enddo

        psurf = prs_sfc
        tsurf = th_sfc*((psurf*rp00)**rovcp)
      IF( imoist.eq.1 )THEN
        if( qv0(1,1,1).lt.qv_pbl )then
          qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)
        else
          qsurf = qv_pbl
        endif
      ELSE
        qsurf = 0.0
      ENDIF

!------------------------------------------------------------------
!  isnd = 6

      ELSEIF(isnd.eq.6)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' isnd = 6 code has been removed (for now)'
        if(dowr) write(outfile,*)
        call stopcm1

!------------------------------------------------------------------
!  isnd = 7
!  Read base-state sounding from an external text file.
!  Assumes file name = "input_sounding"
!  NOTE:  for isnd=7, iwnd is ignored.
!
!  The format is the same as that for the WRF Model.
!
!  One-line header containing:   sfc pres (mb)    sfc theta (K)    sfc qv (g/kg)
!   
!   (Note1: here, "sfc" refers to near-surface atmospheric conditions. 
!    Technically, this should be z = 0, but in practice is obtained from the   
!    standard reporting height of 2 m AGL/ASL from observations)
!   (Note2: land-surface temperature and/or sea-surface temperature (SST) are 
!    specified elsewhere: see tsk0 in namelist.input and/or tsk array in 
!    init_surface.F)
!
! Then, the following lines are:   z (m)    theta (K)   qv (g/kg)    u (m/s)    v (m/s)
!
!   (Note3: # of levels is arbitrary)
!
!     Index:   sfc    =  surface (technically z=0, but typically from 2 m AGL/ASL obs)
!              z      =  height AGL/ASL
!              pres   =  pressure
!              theta  =  potential temperature
!              qv     =  mixing ratio
!              u      =  west-east component of velocity
!              v      =  south-north component of velocity
!
! Note4:  For final line of input_sounding file, z (m) must be greater than the model top
!         (which is nz * dz when stretch_z=0, or ztop when stretch_z=1,  etc)
!

      ELSEIF( isnd.eq.7 .or. isnd.eq.17 )then

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' Reading sounding from external file, input_sounding'
        if(dowr) write(outfile,*)

        nmax = 1000000  !should be enough!

        allocate(   zsnd(nmax) )
        zsnd = 0.0
        allocate(  thsnd(nmax) )
        thsnd = 0.0
        allocate(  qvsnd(nmax) )
        qvsnd = 0.0
        allocate(   usnd(nmax) )
        usnd = 0.0
        allocate(   vsnd(nmax) )
        vsnd = 0.0
        allocate( thvsnd(nmax) )
        thvsnd = 0.0
        allocate(  pisnd(nmax) )
        pisnd = 0.0
        allocate(   psnd(nmax) )
        psnd = 0.0
        allocate(   tsnd(nmax) )
        tsnd = 0.0
        allocate(  rhsnd(nmax) )
        rhsnd = 0.0

        allocate( thinterp(nk) )
        thinterp = 0.0
        allocate( qvinterp(nk) )
        qvinterp = 0.0
        allocate(  uinterp(nk) )
        uinterp = 0.0
        allocate(  vinterp(nk) )
        vinterp = 0.0
        allocate(  pinterp(nk) )
        pinterp = 0.0
        allocate(  tinterp(nk) )
        tinterp = 0.0
        allocate( rhinterp(nk) )
        rhinterp = 0.0

      prc0:  &
      IF( myid.eq.0 )THEN
        ! only processor 0 does this:

        open(unit=40,file='input_sounding',status='old',err=9009)

        ! read surface parameters:
        read(40,*,err=9008) p_sfc, th_sfc, qv_sfc
        if(dowr) write(outfile,*) ' p_sfc, th_sfc, qv_sfc = ',p_sfc, th_sfc, qv_sfc

        ! now, read entire sounding until end of file is discovered
        nsnd=1
        do k=1,nmax
          if( isnd.eq.7 )then
            read(40,*,end=445,err=9008) zsnd(k+1),thsnd(k+1),qvsnd(k+1),usnd(k+1),vsnd(k+1)
          elseif( isnd.eq.17 )then
            read(40,*,end=445,err=9008) zsnd(k+1),thsnd(k+1),qvsnd(k+1)
          else
            print *,'  23985 '
            call stopcm1
          endif
          ! put qv in g/g
          qvsnd(k+1) = qvsnd(k+1)/1000.0
          nsnd=nsnd+1
        enddo
445     continue
        if(dowr) write(outfile,*) '  Found ',nsnd,'  levels (including surface)'
        if(dowr) write(outfile,*)
        close(unit=40)

        if( nsnd.le.2 )then
          print *
          print *,'  There is only one u,v level.  At least two u,v levels are needed. '
          print *
          print *,'  ..... stopping cm1 ..... '
          print *
          call stopcm1
        endif

      ENDIF  prc0

#ifdef MPI
        call MPI_BCAST(nsnd    ,  1 ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(p_sfc   ,  1 ,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(th_sfc  ,  1 ,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qv_sfc  ,  1 ,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zsnd(1) ,nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(thsnd(1),nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qvsnd(1),nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(usnd(1) ,nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vsnd(1) ,nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
#endif

        if(imoist.eq.0) qv_sfc = 0.0

        p_sfc = p_sfc * 100.0
        ! put qv in g/g
        qv_sfc = qv_sfc / 1000.0
        pi_sfc  = (p_sfc/p00)**(rd/cp)
        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)

        psurf  = p_sfc
        thsurf = th_sfc
        tsurf  = th_sfc * pi_sfc
        qsurf  = qv_sfc

        zsnd(1) = 0.0
        thsnd(1) = th_sfc
        qvsnd(1) = qv_sfc

        if(imoist.eq.0) qvsnd = 0.0

      if( isnd.eq.7 )then
        ! dont remember where this formulation comes from:
!!!        usnd(1) = 1.75*usnd(2)-usnd(3)+0.25*usnd(4)
!!!        vsnd(1) = 1.75*vsnd(2)-vsnd(3)+0.25*vsnd(4)

        ! 190324:  use simple extrapolation
        usnd(1) = usnd(2) - (zsnd(2)-zsnd(1))*(usnd(3)-usnd(2))/(zsnd(3)-zsnd(2))
        vsnd(1) = vsnd(2) - (zsnd(2)-zsnd(1))*(vsnd(3)-vsnd(2))/(zsnd(3)-zsnd(2))
      endif

!--------------------------------------------------------------------
!  Added by GHB, 061021:
!  Get thv and prs. Check if qv is too small; if so, set rh to 0.1%
!  (This code has no effect on the sounding if qv > 1e-12 everywere)
        do k=1,nsnd
          thvsnd(k)=thsnd(k)*(1.0+reps*qvsnd(k))/(1.0+qvsnd(k))
        enddo
        pisnd(1)=pi_sfc
        do k=2,nsnd
          pisnd(k)=pisnd(k-1)-g*(zsnd(k)-zsnd(k-1))   &
                               /(cp*0.5*(thvsnd(k)+thvsnd(k-1)))
        enddo
        do k=1,nsnd
          psnd(k)=p00*(pisnd(k)**(cp/rd))
          tsnd(k)=thsnd(k)*pisnd(k)
        enddo
      if(imoist.eq.1)then
        do k=1,nsnd
          if(qvsnd(k).lt.1.0e-12)then
            if(dowr) write(outfile,*) '  Qv is too small.  Setting rh to 0.01%.  k,zsnd=',k,zsnd(k)
            qvsnd(k)=0.0001*rslf(psnd(k),thsnd(k)*pisnd(k))
          endif
          rhsnd(k)=qvsnd(k)/rslf(psnd(k),tsnd(k))
        enddo
      endif
!--------------------------------------------------------------------

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' Data from input_sounding file: '
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '    k,z,th,qv:'
        do k=1,nsnd
          if(dowr) write(outfile,*) k,zsnd(k),thsnd(k),1000.0*qvsnd(k)
        enddo
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '    k,z,u,v:'
        do k=1,nsnd
          if(dowr) write(outfile,*) k,zsnd(k),usnd(k),vsnd(k)
        enddo
        if(dowr) write(outfile,*)

! check to make sure sounding levels span computational grid (WICKER)
!!!        061011, GHB:  Commented this out ... we'll use the surface info.
!!!        if (zsnd(1) .gt. zh(1,1,1)) then
!!!          write(*,*) 'zmin of sounding > zmin of grid!'
!!!          write(*,*) 'zmin of sounding = ',zsnd(1)
!!!          write(*,*) 'zmin of grid = ',zh(1,1,1)
!!!          call stopcm1
!!!        endif

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) 'interpolating sounding to grid: '
        if(dowr) write(outfile,*)

      DO j=jb,je
      DO i=ib,ie

        if (zsnd(nsnd) .lt. zh(i,j,nk)) then
          if(dowr) write(outfile,*) 'zmax of sounding < zmax of grid!'
          if(dowr) write(outfile,*) 'zmax of sounding = ',zsnd(nsnd)
          if(dowr) write(outfile,*) 'zmax of grid  = ',zh(i,j,nk)
          call stopcm1
        endif

        DO k=1,nk

            kk = 1
            do while( zsnd(kk) .lt. zh(i,j,k) )
              kk = kk+1
            enddo
            kdn = kk-1
            kup = kk

            zu=0.5*(zh(max(ib,i-1),j,k)+zh(i,j,k))
            interp_frac = (   zu        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            uinterp(k) =  usnd(kdn) + ( usnd(kup)- usnd(kdn))*interp_frac

            zv=0.5*(zh(i,max(jb,j-1),k)+zh(i,j,k))
            interp_frac = (   zv        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            vinterp(k) =  vsnd(kdn) + ( vsnd(kup)- vsnd(kdn))*interp_frac

!!!            ! if this is first grid point, utilize surface values of th,qv:
!!!            if( k.eq.1 ) kdn = 1

            interp_frac = (   zh(i,j,k) - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            thinterp(k) = thsnd(kdn) + (thsnd(kup)-thsnd(kdn))*interp_frac
            qvinterp(k) = qvsnd(kdn) + (qvsnd(kup)-qvsnd(kdn))*interp_frac
             pinterp(k) =  psnd(kdn) + ( psnd(kup)- psnd(kdn))*interp_frac
             tinterp(k) =  tsnd(kdn) + ( tsnd(kup)- tsnd(kdn))*interp_frac
            rhinterp(k) = rhsnd(kdn) + (rhsnd(kup)-rhsnd(kdn))*interp_frac

            if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '       ',zsnd(kdn),zh(i,j,k),zsnd(kup),interp_frac

        ENDDO
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*)

        do k=1,nk
           u0(i,j,k) =  uinterp(k)
           v0(i,j,k) =  vinterp(k)
          qv0(i,j,k) = qvinterp(k)
          th0(i,j,k) = thinterp(k)
         prs0(i,j,k) =  pinterp(k)
           t0(i,j,k) =  tinterp(k)
          rh0(i,j,k) = rhinterp(k)
        enddo

        ! get pi0 and prs0 from thv0, using hydrostatic equation

        do k=1,nk
          ! get qv from linear interpolation of rh:
          qv0(i,j,k) = rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
        enddo

        do k=1,nk
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
        enddo

        pi0(i,j,1)=pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
        do k=2,nk
          pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))   &
                                   /(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
        enddo

        do k=1,nk
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
        enddo

        ! rh, just in case we want/need it later

      if(imoist.eq.1)then
        do k=1,nk
          rh0(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
        enddo
      endif

      ENDDO    ! enddo for i loop
      ENDDO    ! enddo for j loop


        ! deallocate temporary 1D arrays
        deallocate(   zsnd )
        deallocate(  thsnd )
        deallocate(  qvsnd )
        deallocate(   usnd )
        deallocate(   vsnd )
        deallocate( thvsnd )
        deallocate(  pisnd )
        deallocate(   psnd )
        deallocate(   tsnd )
        deallocate(  rhsnd )
        deallocate( thinterp )
        deallocate( qvinterp )
        deallocate(  uinterp )
        deallocate(  vinterp )
        deallocate(  pinterp )
        deallocate(  tinterp )
        deallocate( rhinterp )


!-----------------------------------------------------------------------
!  isnd = 8
!  Dry, constant d(theta)/dz sounding

      ELSEIF(isnd.eq.8)THEN

        ! Set these three variables for dry, constant d(theta)/dz sounding

        th_sfc   =  300.0     ! theta at surface (K)
        pi_sfc   =    1.0     ! Exner function at surface
        lapse    = -0.0030581 ! potential temperature lapse rate (K/m)

        rhval = 0.0

        IF( testcase.eq.6 .or. testcase.eq.10 )THEN
          !  hurricane boundary layer (Bryan et al, 2017, BLM)
          lapse  =  0.005
        ENDIF
        IF( testcase.eq.10 )THEN
          rhval = 0.90
        ENDIF
        IF( testcase.eq.15 )THEN
          lapse = 0.005
          rhval = 0.80
        ENDIF


        do k=kb,ke
        do j=jb,je
        do i=ib,ie

          ! Calculate theta using the specified lapse rate.
          th0(i,j,k)=th_sfc+lapse*zh(i,j,k)

          ! Calculate pi from theta
          pi0(i,j,k)=pi_sfc-(g/(cp*lapse))*alog(th0(i,j,k)/th_sfc)

          ! Calculate pressure from pi
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))

          ! Calculate temperature from theta and pi
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)

          rh0(i,j,k) = rhval
          qv0(i,j,k) = rhval*rslf(prs0(i,j,k),t0(i,j,k))

        enddo
        enddo
        enddo

        thsurf = th_sfc
        psurf = p00*(pi_sfc**(cp/rd))
        tsurf = th_sfc * pi_sfc
        qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)

!------------------------------------------------------------------
!  isnd = 9
!  Constant Brunt-Vaisala frequency

      ELSEIF(isnd.eq.9)then

        ns1 =  0.0001
        ns2 =  0.0000

        zl1 = 40000.0

      do j=jb,je
      do i=ib,ie
        th_sfc   =  288.00
        pi_sfc   =    1.0
        !-----
        thsurf = th_sfc
        psurf = p00*(pi_sfc**(cp/rd))
        tsurf = th_sfc * pi_sfc
        !-----
        if(zh(i,j,1).lt.zl1)then
          zsfc = 0.0
          ns   = ns1
        else
          if(abs(ns1).lt.1.0e-6)then
            pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
          else
            pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
            th_sfc=th_sfc*exp(ns1*zl1/g)
          endif
          zsfc = zl1
          ns   = ns2
        endif
        do k=1,nk
          qv0(i,j,k)=0.0
          qc0(i,j,k)=0.0
          rh0(i,j,k)=0.0
          if(abs(ns).lt.1.0e-6)then
            thv0(i,j,k)=th_sfc
            pi0(i,j,k)=pi_sfc-g*(zh(i,j,k)-zsfc)/(cp*th_sfc)
          else
            thv0(i,j,k)=th_sfc*exp(ns*(zh(i,j,k)-zsfc)/g)
            pi0(i,j,k)=pi_sfc+g*g/(cp*ns*th_sfc)   &
                           *(exp(-ns*(zh(i,j,k)-zsfc)/g)-1.0)
          endif
          prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
          th0(i,j,k)=thv0(i,j,k)
          do n=1,20
            t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
            qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            th0(i,j,k)=thv0(i,j,k)*(1.0+qv0(i,j,k)+qc0(i,j,k))/(1.0+reps*qv0(i,j,k))
!!!            if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,th0(i,j,k),qv0(i,j,k)
          enddo
          if(zh(i,j,k+1).gt.zl1.and.zsfc.lt.1.0)then
            if(abs(ns1).lt.1.0e-6)then
              pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
            else
              pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
              th_sfc=th_sfc*exp(ns1*zl1/g)
            endif
            zsfc=zl1
            ns=ns2
          endif
        enddo

      enddo
      enddo

        qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)

!------------------------------------------------------------------
!  isnd = 10
!  Moist, constant Brunt-Vaisala frequency
!    (assuming reversible moist microphysics)

      ELSEIF(isnd.eq.10)then

        if(dowr) write(outfile,*) '  isnd = 10'

        qcval    =   0.005
        qtval    =   0.00000

        ns1 =  0.0
        ns2 =  4.0e-4

        zl1 = 111700.0

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        ! DEFINE RH here:  saturated (RH = 100%) is assumed by default
        rh0(i,j,k)=1.0
        qc0(i,j,k)=qcval
      enddo
      enddo
      enddo

      if(dowr) write(outfile,*)
      do j=jb,je
      do i=ib,ie
        th_sfc   =  288.0
        pi_sfc   =    1.0
        !-----
        thsurf = th_sfc
        psurf = p00*(pi_sfc**(cp/rd))
        tsurf = th_sfc * pi_sfc
        !-----
        ns=ns1
        p_sfc=p00*(pi_sfc**(cp/rd))
        t_sfc=th_sfc*pi_sfc
        qv_sfc=rh0(i,j,1)*rslf(p_sfc,t_sfc)
        ql_sfc=qcval
!!!        ql_sfc=qtval-qv_sfc
        thv_sfc=th_sfc*(1.0+reps*qv_sfc)/(1.0+qv_sfc+ql_sfc)
        thn=th_sfc
        thlast=thn
      !---------- do k=1 first ----------!
        k=1
        n=1
        qv0(i,j,k)=qv_sfc
        qc0(i,j,k)=ql_sfc
 571      continue
          th0(i,j,k)=thn
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))   &
                                /(1.0+qv0(i,j,k)+qc0(i,j,k))
          pi0(i,j,k)=pi_sfc-g*zh(i,j,k)/(cp*0.5*(thv0(i,j,k)+thv_sfc))
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
          qv0(i,j,k)=qvs
          qc0(i,j,k)=qcval
!!!          qc0(i,j,k)=qtval-qv0(i,j,k)
          qtm=qv_sfc+ql_sfc
          qtp=qv0(i,j,k)+qc0(i,j,k)
          tavg=0.5*( t0(i,j,k)+t_sfc )
          lhv=lv1-lv2*tavg
          qvavg=0.5*( qv_sfc + qv0(i,j,k) )
          qlavg=0.5*( ql_sfc + qc0(i,j,k) )
          qtavg=qvavg+qlavg
          desdt=17.67*(273.15-29.65)/((tavg-29.65)**2)
          drdt=17.67*(273.15-29.65)*qvavg/((tavg-29.65)**2)
        if(eqtset.eq.2)then
          gamma=g*(1.0+qtavg)*(1.0+lhv*qvavg/rd/tavg)   &
               /( cp+cpv*qvavg+cpl*qlavg          &
                 +lhv*drdt )
        else
          gamma=g*(1.0+qtavg)*(rd/(rd+rv*qvavg)+lhv*qvavg/rd/tavg)   &
               /( cp+lhv*(1.0+qvavg*reps)*qvavg*desdt )
        endif
          tp=t_sfc*exp( zh(i,j,k)*(                     &
                ((ns/g)+alog((1.0+qtp)/(1.0+qtm))/zh(i,j,k))  &
                   /(1.0+tavg*drdt/(eps+qvavg))-gamma/tavg                  &
                                                      ) )
          thn=tp/pi0(i,j,k)

          n=n+1

          if(n.gt.180.and.dowr) write(outfile,*) n,thn
          if(abs(thn-thlast).gt.0.0001 .and. n.lt.200)then
            thn=thlast+0.3*(thn-thlast)
            thlast=thn
            go to 571
          elseif(n.ge.200)then
            if(dowr) write(outfile,*) '  stuck in loop (111)!'
            call stopcm1
          endif

        th0(i,j,k)=thn
        prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
        t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
        qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
        qv0(i,j,k)=qvs
        qc0(i,j,k)=qcval
!!!        qc0(i,j,k)=qtval-qv0(i,j,k)
        if(qc0(i,j,1).lt.0.0)then
          print *,'  calling stopcm1 (222) '
          call stopcm1
        endif
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,th0(i,j,k),prs0(i,j,k)
      !---------- done with k=1, do all other levels ----------!
        qvs=qv0(i,j,1)
        ns=ns1
        do k=2,nk
        tflag=0
        do niter=1,2
        if(tflag.lt.2)then
          if(zh(i,j,k).gt.zl1.and.ns.eq.ns1.and.tflag.eq.0.and.zh(i,j,k-1).lt.zl1)then
            tflag=1
            delz=zl1-zh(i,j,k-1)
            pim=pi0(i,j,k-1)
            thm=th0(i,j,k-1)
            thvm=thv0(i,j,k-1)
            qtm=qv0(i,j,k-1)+qc0(i,j,k-1)+qi0(i,j,k-1)
            qvm=qv0(i,j,k-1)
            qlm=qc0(i,j,k-1)
            qim=qi0(i,j,k-1)
            tm=t0(i,j,k-1)
            pm=prs0(i,j,k-1)
          elseif(tflag.eq.1)then
            tflag=0
            ns=ns2
            delz=zh(i,j,k)-zl1
            pim=pi0(i,j,k)
            thm=th0(i,j,k)
            thvm=thv0(i,j,k)
            qtm=qv0(i,j,k)+qc0(i,j,k)+qi0(i,j,k)
            qvm=qv0(i,j,k)
            qlm=qc0(i,j,k)
            qim=qi0(i,j,k)
            tm=t0(i,j,k)
            pm=prs0(i,j,k)
          else
            tflag=2
            delz=zh(i,j,k)-zh(i,j,k-1)
            pim=pi0(i,j,k-1)
            thm=th0(i,j,k-1)
            thvm=thv0(i,j,k-1)
            qtm=qv0(i,j,k-1)+qc0(i,j,k-1)+qi0(i,j,k-1)
            qvm=qv0(i,j,k-1)
            qlm=qc0(i,j,k-1)
            qim=qi0(i,j,k-1)
            tm=t0(i,j,k-1)
            pm=prs0(i,j,k-1)
          endif
!!!          t0(i,j,k)=tm
          th0(i,j,k)=thm
          qv0(i,j,k)=qvm
          qc0(i,j,k)=qlm
          qi0(i,j,k)=qim
          thv0(i,j,k)=thvm
!!!          thlast=tm*pim
          thlast=thm
          n=0
          thn=thlast
          thlast=thn
 572      continue
            th0(i,j,k)=thn
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))   &
                                  /(1.0+qv0(i,j,k)+qc0(i,j,k)+qi0(i,j,k))
            pi0(i,j,k)=pim-g*delz/(cp*0.5*(thv0(i,j,k)+thvm))
            prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
            t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          if(iice.eq.0)then
            qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qv0(i,j,k)=qvs
            qc0(i,j,k)=qcval
!!!            qc0(i,j,k)=qtval-qv0(i,j,k)
            qtp=qv0(i,j,k)+qc0(i,j,k)
            prsavg=0.5*(prs0(i,j,k-1)+prs0(i,j,k))
            tavg=0.5*( t0(i,j,k)+tm )
            qvavg=rh0(i,j,k)*rslf(prsavg,tavg)
            qlavg=0.5*( qlm + qc0(i,j,k) )
            qtavg=0.5*( qtm + qtp )
            drdt=17.67*(273.15-29.65)*qvavg/((tavg-29.65)**2)
            lhv=lv1-lv2*tavg
            cpml=cp+cpv*qvavg+cpl*qlavg
          else
            qvl=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qvi=rh0(i,j,k)*rsif(prs0(i,j,k),t0(i,j,k))
            fliq=max(min((t0(i,j,k)-t00k)*rt0,1.0),0.0)
            fice=1.0-fliq
            qvs=fliq*qvl+fice*qvi
            qv0(i,j,k)=qvs
            qc0(i,j,k)=fliq*qcval
            qi0(i,j,k)=fice*qcval
            qtp=qv0(i,j,k)+qc0(i,j,k)+qi0(i,j,k)
            prsavg=0.5*(prs0(i,j,k-1)+prs0(i,j,k))
            tavg=0.5*( t0(i,j,k)+tm )
            qvl=rh0(i,j,k)*rslf(prsavg,tavg)
            qvi=rh0(i,j,k)*rsif(prsavg,tavg)
            fliq=max(min((tavg-t00k)*rt0,1.0),0.0)
            fice=1.0-fliq
            qvavg=fliq*qvl+fice*qvi
            qlavg=0.5*( qlm + qc0(i,j,k) )
            qtavg=0.5*( qtm + qtp )
            drdt=fliq*17.67*(273.15-29.65)*qvl/((tavg-29.65)**2)    &
                +fice*21.8745584*(273.15-7.66)*qvi/((tavg-7.66)**2)
            if(tavg.gt.t00k.and.tavg.lt.t0k)then
              drdt=drdt+(qvl-qvi)*rt0
            endif
            lhv=fliq*(lv1-lv2*tavg)+fice*(ls1-ls2*tavg)
            cpml=cp+cpv*qvavg+cpl*qlavg+cpi*(qtavg-qlavg-qvavg)
          endif
          if(eqtset.eq.2)then
            gamma=g*(1.0+qtavg)*(1.0+lhv*qvavg/rd/tavg)   &
                 /( cpml+lhv*drdt )
          else
            gamma=g*(1.0+qtavg)*(rd/(rd+rv*qvavg)+lhv*qvavg/rd/tavg)   &
                 /( cp+lhv*(1.0+qvavg*reps)*qvavg*desdt )
          endif
            tp=tm*exp( delz*(                                 &
                  ((ns/g)+alog((1.0+qtp)/(1.0+qtm))/delz)     &
                     /(1.0+tavg*drdt/(eps+qvavg))-gamma/tavg  &
                                                        ) )
            thn=tp/pi0(i,j,k)

            n=n+1

            if(n.gt.180.and.dowr) write(outfile,*) n,tp,tm,delz
            if(abs(thn-thlast).gt.0.0001 .and. n.lt.200)then
              thn=thlast+0.3*(thn-thlast)
              thlast=thn
              go to 572
            elseif(n.ge.200)then
              if(dowr) write(outfile,*) '  stuck in loop (333)!'
              call stopcm1
            endif

          t0(i,j,k)=tp
          th0(i,j,k)=thn
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          if(iice.eq.0)then
            qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qv0(i,j,k)=qvs
            fliq = 1.0
            fice = 0.0
          else
            qvl=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qvi=rh0(i,j,k)*rsif(prs0(i,j,k),t0(i,j,k))
            fliq=max(min((t0(i,j,k)-t00k)*rt0,1.0),0.0)
            fice=1.0-fliq
            qvs=fliq*qvl+fice*qvi
            qv0(i,j,k)=qvs
            qc0(i,j,k)=fliq*qcval
            qi0(i,j,k)=fice*qcval
          endif
          if(qc0(i,j,1).lt.0.0)then
            print *,'  calling stopcm1 (444) '
            call stopcm1
          endif
!!!          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,t0(i,j,k),fliq,fice
        endif
        enddo  ! enddo for iteration loop
        enddo  ! enddo for k loop

      enddo  ! enddo for i loop
      enddo  ! enddo for j loop
      if(dowr) write(outfile,*)

        qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)

!------------------------------------------------------------------
!  constant theta-e, saturated
!  Reference:  Bryan and Rotunno, 2009, JAS, v10, pp. 3042-3060

      ELSEIF(isnd.eq.11)THEN

        IF(imoist.eq.0)THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  isnd=11 requires imoist=1'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  stopping model ...'
          if(dowr) write(outfile,*)
          call stopcm1
        ENDIF

        zl1 = 15000.0    ! tropopause height

        ttype =  2    ! 1 = reversible
                      ! 2 = pseudoadiabatic

        psurf = 101510.0    ! surface pressure (Pa)

      DO j=jb,je
      DO i=ib,ie

        ! use psurf:
        pi_sfc   =  (psurf*rp00)**rovcp
        ! use tsk from namelist.input, make 3 K cooler:
        th_sfc   =  ( tsk0 - 3.0 )/pi_sfc

        !-----
        thsurf = th_sfc
        tsurf = th_sfc * pi_sfc
        !-----

        t_sfc=th_sfc*pi_sfc
        p_sfc=p00*(pi_sfc**(cp/rd))
        qv_sfc=rslf(p_sfc,t_sfc)
        qt_sfc=qv_sfc
        thv_sfc=th_sfc*(1.0+reps*qv_sfc)/(1.0+qv_sfc)

        the_sfc=getthe(ttype,p_sfc,t_sfc,qv_sfc,qt_sfc)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '  p_sfc,t_sfc,qv_sfc = ',p_sfc,t_sfc,qv_sfc
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '  the_sfc,qt_sfc = ',the_sfc,qt_sfc

        th_sfc=0.0
        qv_sfc=0.0
        call revthe(ttype,the_sfc,p_sfc,qt_sfc,t_sfc,qv_sfc)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '  p_sfc,t_sfc,qv_sfc = ',p_sfc,t_sfc,qv_sfc

        pi1=pi_sfc
        z1=0.0
        t1=t_sfc
        thv1=thv_sfc
        thv2=thv1
        flag=0

        do k=1,nk

          z2=zh(i,j,k)
          t2=t1
          tlast=0.0
          n=0

        if(z2.le.zl1)then
          ! troposphere
          theq=the_sfc
          ! iterate:
          do while( abs(t2-tlast).gt.0.0001 )
            tlast=t2
            n=n+1
            pi2=pi1-g*(z2-z1)/(cp*0.5*(thv1+thv2))
            p2=p00*(pi2**(cp/rd))
            call revthe(ttype,theq,p2,qt_sfc,t2,qv2)
            th2=t2/pi2
            if(ttype.eq.1)then
              ql2=qt_sfc-qv2
            elseif(ttype.eq.2)then
              ql2=0.0
            endif
            thv2=th2*(1.0+reps*qv2)/(1.0+qv2+ql2)
          enddo
        else
          ! stratosphere
          ns=4.0e-4
          if(flag.eq.0)then
            ! first time in this section ... get sfc params
            !-------------------
            z2=zl1
            t2=t1
            tlast=0.0
            n=0
            theq=the_sfc
            do while( abs(t2-tlast).gt.0.0001 )
              tlast=t2
              n=n+1
              pi2=pi1-g*(z2-z1)/(cp*0.5*(thv1+thv2))
              p2=p00*(pi2**(cp/rd))
              call revthe(ttype,theq,p2,qt_sfc,t2,qv2)
              th2=t2/pi2
              if(ttype.eq.1)then
                ql2=qt_sfc-qv2
              elseif(ttype.eq.2)then
                ql2=0.0
              endif
              thv2=th2*(1.0+reps*qv2)/(1.0+qv2+ql2)
            enddo
            flag=1
            th_sfc=thv2
            pi_sfc=pi2
            zsfc=z2
            !-------------------
            z2=zh(i,j,k)
            t2=t1
            tlast=0.0
            n=0
          endif
          if(abs(ns).lt.1.0e-6)then
            thv2=th_sfc
            pi2=pi_sfc-g*(z2-zsfc)/(cp*th_sfc)
          else
            thv2=th_sfc*exp(ns*(z2-zsfc)/g)
            pi2=pi_sfc+g*g/(cp*ns*th_sfc)   &
                           *(exp(-ns*(z2-zsfc)/g)-1.0)
          endif
          p2=p00*(pi2**cpdrd)
          th2=thv2
          do n=1,20
            t2=th2*pi2
            qv2=rslf(p2,t2)
            if(ttype.eq.1)then
              ql2=qt_sfc-qv2
            elseif(ttype.eq.2)then
              ql2=0.0
            endif
            th2=thv2*(1.0+qv2+ql2)/(1.0+reps*qv2)
          enddo
        endif
!!!          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) n,p2,th2

          t1=t2
          thv1=thv2
          z1=z2
          pi1=pi2

          pi0(i,j,k)=pi2
          prs0(i,j,k)=p2
          t0(i,j,k)=t2
          th0(i,j,k)=th2
          thv0(i,j,k)=thv2
          rh0(i,j,k)=1.0
          qv0(i,j,k)=qv2
          qc0(i,j,k)=ql2

        enddo

      ENDDO
      ENDDO

        qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)

!------------------------------------------------------------------
!  PBL simulation:  assumed dry

      ELSEIF(isnd.eq.12)THEN

        pi_sfc = 1.0
        th_sfc = 300.0
        zl1    = 960.0
        lapse  = 0.010

        pisfc = pi_sfc-g*zl1/(cp*th_sfc)

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          IF(zh(i,j,k).le.zl1)THEN
            th0(i,j,k) = th_sfc
            pi0(i,j,k)=pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          ELSE
            th0(i,j,k) = th_sfc+lapse*(zh(i,j,k)-zl1)
            pi0(i,j,k)=pisfc-(g/(cp*lapse))*alog(th0(i,j,k)/th_sfc)
          ENDIF
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          rh0(i,j,k)=0.0
          qv0(i,j,k)=0.0
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  from Klemp (2011)

      ELSEIF(isnd.eq.13)then

        ns1 =  0.0001
        ns2 =  0.0004
        ns3 =  0.0001

        zl1 =  2000.0
        zl2 =  3000.0

      do j=jb,je
      do i=ib,ie
        th_sfc   =  288.00
        pi_sfc   =    1.0
        !-----
        thsurf = th_sfc
        psurf = p00*(pi_sfc**(cp/rd))
        tsurf = th_sfc * pi_sfc
        !-----
        if(zh(i,j,1).lt.zl1)then
          zsfc = 0.0
          ns   = ns1
        else
          if(abs(ns1).lt.1.0e-6)then
            pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
          else
            pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
            th_sfc=th_sfc*exp(ns1*zl1/g)
          endif
          zsfc = zl1
          ns   = ns2
        endif
        do k=1,nk
          if(abs(ns).lt.1.0e-6)then
            th0(i,j,k)=th_sfc
            pi0(i,j,k)=pi_sfc-g*(zh(i,j,k)-zsfc)/(cp*th_sfc)
          else
            th0(i,j,k)=th_sfc*exp(ns*(zh(i,j,k)-zsfc)/g)
            pi0(i,j,k)=pi_sfc+g*g/(cp*ns*th_sfc)   &
                           *(exp(-ns*(zh(i,j,k)-zsfc)/g)-1.0)
          endif
          prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
          thv0(i,j,k)=th0(i,j,k)
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          if(zh(i,j,k+1).gt.zl1.and.zsfc.lt.1.0)then
            if(abs(ns1).lt.1.0e-6)then
              pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
            else
              pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
              th_sfc=th_sfc*exp(ns1*zl1/g)
            endif
            zsfc=zl1
            ns=ns2
          elseif(zh(i,j,k+1).gt.zl2.and.zsfc.lt.zl2)then
            if(abs(ns1).lt.1.0e-6)then
              pi_sfc=pi_sfc-g*(zl2-zl1)/(cp*th_sfc)
            else
              pi_sfc=pi_sfc+g*g/(cp*ns2*th_sfc)*(exp(-ns2*(zl2-zl1)/g)-1.0)
              th_sfc=th_sfc*exp(ns2*(zl2-zl1)/g)
            endif
            zsfc=zl2
            ns=ns3
          endif
        enddo

      enddo
      enddo

!------------------------------------------------------------------
!  Convective boundary layer  (Sullivan and Patton, 2011, JAS)

      ELSEIF(isnd.eq.14)THEN

        pi_sfc = 1.0
        th_sfc = 300.0

        lr1 =  0.08       ! lapse rate within inversion
        lr2 =  0.003      ! lapse rate above inversion

        zl1 =   974.0     ! bottom of inversion
        zl2 =  1074.0     ! top of inversion

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          IF(zh(i,j,k).le.zl1)THEN
            th0(i,j,k) = th_sfc
            pi0(i,j,k) = pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          ELSEIF(zh(i,j,k).le.zl2)THEN
            th0(i,j,k) = th_sfc + (zh(i,j,k)-zl1)*lr1
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(th0(i,j,k)+th0(i,j,k-1)))
          ELSE
            th0(i,j,k) = th_sfc + (zl2-zl1)*lr1 + (zh(i,j,k)-zl2)*lr2
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(th0(i,j,k)+th0(i,j,k-1)))
          ENDIF
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          thv0(i,j,k)=th0(i,j,k)
          rh0(i,j,k)=0.0
          qv0(i,j,k)=0.0
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------
! isnd = 15
      
      ! base state for the simulation of DYCOMS-II RF01,
      ! nonprecipitating stratocumulus (Stevens et al, 2005, MWR)            

      ! and DYCOMS-II RF02
      ! drizzling stratocumulus (Ackerman et al, 2009, MWR)

      ELSEIF( isnd.eq.15 )THEN

        if( testcase.lt.4 .or. testcase.gt.5 )then
          if( myid.eq.0 )then
            print *
            print *,'  isnd       = ',isnd
            print *,'  testcase = ',testcase
            print *
            print *,'  testcase must be 4 or 5 for this case '
            print *
            print *,'   stopping model .... '
            print *
            print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
         
        do j = 1,nj
        do i = 1,ni

          if( testcase .eq. 4 )then
            th_sfc  =  289.0                   ! surface (air) potential temperature (K)
            qv_sfc  =  0.009                   ! surface mixing ratio (kg/kg)
            prs_sfc = 101780.0                 ! surface pressure (Pa)
            z_inv   =  840.0                   ! base of the inversion (m)
          elseif( testcase .eq. 5 )then
            th_sfc  =  288.3                   ! surface (air) potential temperature (K)
            qv_sfc  =  0.00945                 ! surface mixing ratio (kg/kg)
            prs_sfc = 101780.0                 ! surface pressure (Pa)
            z_inv   =  795.0                   ! base of the inversion (m)
          endif

          pi_sfc  = (prs_sfc/p00)**(rd/cp)                 ! Exner function at surface
          thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)  ! surface density potential temperature (K)

          thvlast = thv_sfc
          th0last = th_sfc
          zlast = 0.0
          pilast = pi_sfc

          if( myid.eq.0 .and. i.eq.1 .and. j.eq.1 ) print *,'  pi_sfc,z_inv = ',pi_sfc,z_inv

          kloop1:  do k = 1,nk

            if( zh(i,j,k) .le. z_inv )then
              thl = th_sfc
              qt = qv_sfc
            else
              if( testcase .eq. 4 )then
                thl = 297.5 + (zh(i,j,k)-z_inv)**(1.0/3.0)
                qt = 0.0015
                ! convert specific humidity to mixing ratio:
                qt = qt/(1.0-qt)
              elseif( testcase .eq. 5 )then
                thl = 295.0 + (zh(i,j,k)-z_inv)**(1.0/3.0)
                qt = 0.001*( 5.0-3.0*(1.0-exp((z_inv-zh(i,j,k))/500.0)) )
              endif
            endif

            if( myid.eq.0 .and. i.eq.1 .and. j.eq.1 ) print *
            if( myid.eq.0 .and. i.eq.1 .and. j.eq.1 ) print *,'  k,zh = ',k,zh(i,j,k)
            if( myid.eq.0 .and. i.eq.1 .and. j.eq.1 ) print *,'  pi0,th0,qc0:'

            ! first guess:
            thlast = 0.0
            th0(i,j,k) = thl
            qv0(i,j,k) = qt
            qc0(i,j,k) = 0.0
            n = 0

            ! iterate:
            do while( abs(th0(i,j,k)-thlast).ge.0.0001 )

              thlast = th0(i,j,k)
              thv0(i,j,k) = th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k)+qc0(i,j,k))
              pi0(i,j,k) = pilast - g*(zh(i,j,k)-zlast)/(cp*0.5*(thvlast+thv0(i,j,k)))
              prs0(i,j,k) = p00*(pi0(i,j,k)**(cp/rd))
              t0(i,j,k) = th0(i,j,k)*pi0(i,j,k)
              qvs = rslf( prs0(i,j,k) , t0(i,j,k) )
              if( qt .lt. qvs )then
                qv0(i,j,k) = qt
                qc0(i,j,k) = 0.0
                th0(i,j,k) = thl
              else
                qv0(i,j,k) = qvs
                qc0(i,j,k) = qt-qvs
                th0(i,j,k) = thl/(1.0-xlv*qc0(i,j,k)/(cp*t0(i,j,k)))
              endif

              th0(i,j,k) = thlast + 0.3*(th0(i,j,k)-thlast)

              n = n+1
              if( myid.eq.0 .and. i.eq.1 .and. j.eq.1 ) print *,n,pi0(i,j,k),th0(i,j,k),qc0(i,j,k)
              if( n.gt.100.0 ) stop 1111

            enddo

            t0(i,j,k) = th0(i,j,k)*pi0(i,j,k)
            rho0(i,j,k) = prs0(i,j,k)/(rd*t0(i,j,k)*(1.0+qv0(i,j,k)*reps))

            thvlast = thv0(i,j,k)
            th0last = th0(i,j,k)
            pilast = pi0(i,j,k)
            zlast = zh(i,j,k)

          enddo  kloop1

        enddo
        enddo
 
        psurf = prs_sfc
        qsurf = qv_sfc
        tsurf = th_sfc * pi_sfc

        lwp = 0.0
        do k=1,nk
          lwp = lwp + rho0(1,1,k)*qc0(1,1,k)*dz*rmh(1,1,k)
        enddo

        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  lwp = ',lwp
        if( myid.eq.0 ) print *
 

!------------------------------------------------------------------
!   Moeng and Sullivan (1994, JAS) PBL

      ELSEIF(isnd.eq.18)THEN

        pi_sfc = 1.0
        th_sfc = 300.0

        delz = 60.0

        lr1 =  8.0/delz   ! lapse rate within inversion
        lr2 =  0.003      ! lapse rate above inversion

        zl1 =   465.0 - delz/2.0   ! bottom of inversion
        zl2 =   465.0 + delz/2.0   ! top of inversion

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          IF(zh(i,j,k).le.zl1)THEN
            th0(i,j,k) = th_sfc
            pi0(i,j,k) = pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          ELSEIF(zh(i,j,k).le.zl2)THEN
            th0(i,j,k) = th_sfc + (zh(i,j,k)-zl1)*lr1
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(th0(i,j,k)+th0(i,j,k-1)))
          ELSE
            th0(i,j,k) = th_sfc + (zl2-zl1)*lr1 + (zh(i,j,k)-zl2)*lr2
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(th0(i,j,k)+th0(i,j,k-1)))
          ENDIF
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          thv0(i,j,k)=th0(i,j,k)
          rh0(i,j,k)=0.0
          qv0(i,j,k)=0.0
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!   Siebesma et al (2003, JAS) shallow Cu case:

      ELSEIF(isnd.eq.19)THEN

        ! (m)
        z1   =     0.0
        z2   =   520.0
        z3   =  1480.0
        z4   =  2000.0
        z5   =  3000.0

        ! (K)
        th1  =   298.7
        th2  =   298.7
        th3  =   302.4
        th4  =   308.2
        th5  =   311.85

        ! specific humidity (convert to mixing ratio below)
        ! (g/g)
        qv1  =  0.0170
        qv2  =  0.0163
        qv3  =  0.0107
        qv4  =  0.0042
        qv5  =  0.0030

        prs_sfc  =  101500.0   ! (Pa)

      IF( testcase.eq.11 )THEN
        ! SAS convective PBL:
        ! (m)
        z1   =     0.0
        z2   =   352.5
        z3   =   442.5
        z4   =  3000.0
        z5   = 20000.0

        ! (K)
        th1  =   296.6
        th2  =   296.6
        th3  =   298.1
        th4  =   th3 + 0.003*(z4-z3)
        th5  =   th3 + 0.003*(z5-z3)

        ! specific humidity (convert to mixing ratio below)
        ! (g/g)
        qv1  =  11.8e-3
        qv2  =  11.8e-3
        qv3  =   7.8e-3
        qv4  =  qv3 - 0.004e-3*(z4-z3)
        qv5  =  qv3 - 0.004e-3*(z5-z3)

        prs_sfc  =  100000.0   ! (Pa)
      ENDIF

        pi_sfc   =  (prs_sfc*rp00)**(rovcp)
        thv_sfc  = th1*(1.0+reps*qv1)/(1.0+qv1)

        do j=1,nj
        do i=1,ni

          do k=kb,ke
            if(     zh(i,j,k).le.z2 )then
              th0(i,j,k) = th1+(th2-th1)*(zh(i,j,k)-z1)/(z2-z1)
              qv0(i,j,k) = qv1+(qv2-qv1)*(zh(i,j,k)-z1)/(z2-z1)
            elseif( zh(i,j,k).le.z3 )then
              th0(i,j,k) = th2+(th3-th2)*(zh(i,j,k)-z2)/(z3-z2)
              qv0(i,j,k) = qv2+(qv3-qv2)*(zh(i,j,k)-z2)/(z3-z2)
            elseif( zh(i,j,k).le.z4 )then
              th0(i,j,k) = th3+(th4-th3)*(zh(i,j,k)-z3)/(z4-z3)
              qv0(i,j,k) = qv3+(qv4-qv3)*(zh(i,j,k)-z3)/(z4-z3)
            elseif( zh(i,j,k).le.z5 )then
              th0(i,j,k) = th4+(th5-th4)*(zh(i,j,k)-z4)/(z5-z4)
              qv0(i,j,k) = qv4+(qv5-qv4)*(zh(i,j,k)-z4)/(z5-z4)
            endif
            ! convert specific humidity to mixing ratio:
            qv0(i,j,k) = max( 0.0 , qv0(i,j,k)/(1.0-qv0(i,j,k)) )
            thv0(i,j,k) = th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
          enddo

          pi0(i,j,1) = pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
          do k=2,nk
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
          enddo

          do k=1,nk
            prs0(i,j,k) = p00*(pi0(i,j,k)**(cp/rd))
            rh0(i,j,k) = qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
          enddo

      enddo
      enddo


        psurf = prs_sfc
        tsurf = th1*((psurf*rp00)**rovcp)
      IF( imoist.eq.1 )THEN
        if( qv0(1,1,1).lt.qv1 )then
          qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)
        else
          qsurf = qv1
        endif
      ELSE
        qsurf = 0.0
      ENDIF

!------------------------------------------------------------------
!   RICO, precipitating shallow Cu  (VanZanten et al 2011, JAMES)

      ELSEIF(isnd.eq.20)THEN

        ! (m)
        z1   =     0.0
        z2   =   740.0
        z3   =  3260.0
        z4   =  4000.0

        ! (K)
        th1  =   297.9
        th2  =   297.9
        th3  =   312.6644
        th4  =   317.0

        ! specific humidity (convert to mixing ratio below)
        ! (g/g)
        qv1  =  0.0160
        qv2  =  0.0138
        qv3  =  0.0024
        qv4  =  0.0018

        ! (Pa)
        prs_sfc  =  101540.0

        pi_sfc   =  (prs_sfc*rp00)**(rovcp)
        thv_sfc  = th1*(1.0+reps*qv1)/(1.0+qv1)

        do j=1,nj
        do i=1,ni

          do k=kb,ke
            if(     zh(i,j,k).le.z2 )then
              th0(i,j,k) = th1+(th2-th1)*(zh(i,j,k)-z1)/(z2-z1)
              qv0(i,j,k) = qv1+(qv2-qv1)*(zh(i,j,k)-z1)/(z2-z1)
            elseif( zh(i,j,k).le.z3 )then
              th0(i,j,k) = th2+(th3-th2)*(zh(i,j,k)-z2)/(z3-z2)
              qv0(i,j,k) = qv2+(qv3-qv2)*(zh(i,j,k)-z2)/(z3-z2)
            elseif( zh(i,j,k).le.z4 )then
              th0(i,j,k) = th3+(th4-th3)*(zh(i,j,k)-z3)/(z4-z3)
              qv0(i,j,k) = qv3+(qv4-qv3)*(zh(i,j,k)-z3)/(z4-z3)
            endif
            ! convert specific humidity to mixing ratio:
            qv0(i,j,k) = qv0(i,j,k)/(1.0-qv0(i,j,k))
            thv0(i,j,k) = th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
          enddo

          pi0(i,j,1) = pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
          do k=2,nk
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
          enddo

          do k=1,nk
            prs0(i,j,k) = p00*(pi0(i,j,k)**(cp/rd))
            rh0(i,j,k) = qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
          enddo

      enddo
      enddo


        psurf = prs_sfc
        tsurf = th1*((psurf*rp00)**rovcp)
      IF( imoist.eq.1 )THEN
        if( qv0(1,1,1).lt.qv1 )then
          qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)
        else
          qsurf = qv1
        endif
      ELSE
        qsurf = 0.0
      ENDIF

!------------------------------------------------------------------
!   RCEMIP  (Wing et al. 2018, GMD)

      ELSEIF(isnd.eq.21)THEN

        ztrop = 15000.0

        qsfc = 0.01865
        qt = 10.0**(-11)
        zq1 = 4000.0
        zq2 = 7500.0

        tsfc = tsk0           ! set air temperature equal to SST 
        tvsfc = tsfc*(1.0+0.608*qsfc)

        psfc = 101480.0
        pisfc = (psfc/p00)**(rd/cp)

        lr = 0.0067

        tvtrop = tvsfc - lr*ztrop
        ptrop = psfc*(tvtrop/tvsfc)**(g/(rd*lr))

        th_sfc = tsfc/pisfc
        qv_sfc = qsfc/(1.0-qsfc)
        thv_sfc = th_sfc*(1.0+reps*qv_sfc)/(1.0+qv_sfc)

        if(dowr) write(outfile,*) '  tvsfc,ptrop = ',tvsfc,ptrop

        do j=jb,je
        do i=ib,ie

          do k=1,nk
            if( zh(i,j,k).lt.ztrop )then
              qv0(i,j,k) = qsfc*exp(-zh(i,j,k)/zq1)*exp(-(zh(i,j,k)/zq2)**2)
              tv0 = tvsfc - lr*zh(i,j,k)
              prs0(i,j,k) = psfc*( ((tvsfc-lr*zh(i,j,k))/tvsfc)**(g/(rd*lr)) )
            else
              qv0(i,j,k) = qt
              tv0 = tvtrop
              prs0(i,j,k) = ptrop*exp( -(g*(zh(i,j,k)-ztrop)/(rd*tvtrop)) )
            endif
            t0(i,j,k) = tv0/(1.0+0.608*qv0(i,j,k))
            if(dowr.and.i.eq.1.and.j.eq.1) write(outfile,*) k,zh(i,j,k),qv0(i,j,k),t0(i,j,k),prs0(i,j,k)
            pi0(i,j,k) = (prs0(i,j,k)*rp00)**rovcp
            th0(i,j,k) = t0(i,j,k)/pi0(i,j,k)
            ! convert specific humidity to mixing ratio:
            qv0(i,j,k) = qv0(i,j,k)/(1.0-qv0(i,j,k))
            thv0(i,j,k) = th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
          enddo

        enddo
        enddo

!------------------------------------------------------------------
!   stable boundary layer test case (Beare et al. 2006, BLM, pg 247)

      ELSEIF( isnd.eq.22 )THEN

        pi_sfc = 1.0
        th_sfc = 265.0
        zl1 = 100.0
        lr1 = 0.01

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          IF(zh(i,j,k).le.zl1)THEN
            th0(i,j,k) = th_sfc
            pi0(i,j,k) = pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          ELSE
            th0(i,j,k) = th_sfc + (zh(i,j,k)-zl1)*lr1
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(th0(i,j,k)+th0(i,j,k-1)))
          ENDIF
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          thv0(i,j,k)=th0(i,j,k)
          rh0(i,j,k)=0.0
          qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!   Diurnal cycle of shallow cumulus (Brown et al, 2002, QJRMS)

      ELSEIF( isnd.eq.23 )THEN

        pi_sfc = 1.0
        th_sfc = 299.0
        qv_sfc = 0.0152

        ! (m)
        z1   =     0.0
        z2   =    50.0
        z3   =   350.0
        z4   =   650.0
        z5   =   700.0
        z6   =  1300.0
        z7   =  2500.0
        z8   =  5500.0

        ! (K)
        th1  =   299.0
        th2  =   301.5
        th3  =   302.5
        th4  =   303.53
        th5  =   303.7
        th6  =   307.13
        th7  =   314.0
        th8  =   343.2

        ! specific humidity (convert to mixing ratio below)
        ! (g/g)
        qv1  =  0.0152
        qv2  =  0.01517
        qv3  =  0.01498
        qv4  =  0.0148
        qv5  =  0.0147
        qv6  =  0.0135
        qv7  =  0.003
        qv8  =  0.003

        prs_sfc  = 97000.0

        pi_sfc   =  (prs_sfc*rp00)**(rovcp)
        thv_sfc  = th1*(1.0+reps*qv1)/(1.0+qv1)

        z9 = 10000.0
        th9 = th8+(z9-z8)*(th8-th7)/(z8-z7)
        qv9 = qv8+(z9-z8)*(qv8-qv7)/(z8-z7)

        do j=1,nj
        do i=1,ni

          do k=kb,ke
            if(     zh(i,j,k).le.z2 )then
              th0(i,j,k) = th1+(th2-th1)*(zh(i,j,k)-z1)/(z2-z1)
              qv0(i,j,k) = qv1+(qv2-qv1)*(zh(i,j,k)-z1)/(z2-z1)
            elseif( zh(i,j,k).le.z3 )then
              th0(i,j,k) = th2+(th3-th2)*(zh(i,j,k)-z2)/(z3-z2)
              qv0(i,j,k) = qv2+(qv3-qv2)*(zh(i,j,k)-z2)/(z3-z2)
            elseif( zh(i,j,k).le.z4 )then
              th0(i,j,k) = th3+(th4-th3)*(zh(i,j,k)-z3)/(z4-z3)
              qv0(i,j,k) = qv3+(qv4-qv3)*(zh(i,j,k)-z3)/(z4-z3)
            elseif( zh(i,j,k).le.z5 )then
              th0(i,j,k) = th4+(th5-th4)*(zh(i,j,k)-z4)/(z5-z4)
              qv0(i,j,k) = qv4+(qv5-qv4)*(zh(i,j,k)-z4)/(z5-z4)
            elseif( zh(i,j,k).le.z6 )then
              th0(i,j,k) = th5+(th6-th5)*(zh(i,j,k)-z5)/(z6-z5)
              qv0(i,j,k) = qv5+(qv6-qv5)*(zh(i,j,k)-z5)/(z6-z5)
            elseif( zh(i,j,k).le.z7 )then
              th0(i,j,k) = th6+(th7-th6)*(zh(i,j,k)-z6)/(z7-z6)
              qv0(i,j,k) = qv6+(qv7-qv6)*(zh(i,j,k)-z6)/(z7-z6)
            elseif( zh(i,j,k).le.z8 )then
              th0(i,j,k) = th7+(th8-th7)*(zh(i,j,k)-z7)/(z8-z7)
              qv0(i,j,k) = qv7+(qv8-qv7)*(zh(i,j,k)-z7)/(z8-z7)
            elseif( zh(i,j,k).le.z9 )then
              th0(i,j,k) = th8+(th9-th8)*(zh(i,j,k)-z8)/(z9-z8)
              qv0(i,j,k) = qv8+(qv9-qv8)*(zh(i,j,k)-z8)/(z9-z8)
            endif
            ! convert specific humidity to mixing ratio:
            qv0(i,j,k) = qv0(i,j,k)/(1.0-qv0(i,j,k))
            thv0(i,j,k) = th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
          enddo

          pi0(i,j,1) = pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
          do k=2,nk
            pi0(i,j,k) = pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
          enddo

          do k=1,nk
            prs0(i,j,k) = p00*(pi0(i,j,k)**(cp/rd))
            rh0(i,j,k) = qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
          enddo

      enddo
      enddo


        psurf = prs_sfc
        tsurf = th1*((psurf*rp00)**rovcp)
      IF( imoist.eq.1 )THEN
        if( qv0(1,1,1).lt.qv1 )then
          qsurf = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)
        else
          qsurf = qv1
        endif
      ELSE
        qsurf = 0.0
      ENDIF

!------------------------------------------------------------------

      ENDIF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  End definition of base state sounding (isnd opton)  cccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!------------------------------------------------------------------
!  fill in ghost cells

      call bcs(pi0)
      call bcs(prs0)
      call bcs(th0)
      call bcs(qv0)
      call bcs(qc0)
      call bcs(qi0)
      call bcs(rh0)

#ifdef MPI
      nf=0
      nu=0
      nv=0
      nw=0
      call comm_all_s( pi0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
      call comm_all_s(prs0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
      call comm_all_s( th0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
      call comm_all_s( qv0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
      call comm_all_s( qc0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
      call comm_all_s( qi0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
      call comm_all_s( rh0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
#endif

      call extrapbcs(pi0)
      call extrapbcs(prs0)
      call extrapbcs(th0)
      call extrapbcs(qv0)
      call extrapbcs(qc0)
      call extrapbcs(qi0)
      call extrapbcs(rh0)

    do j=jb,je
    do i=ib,ie

!------------------------------------------------------------------
!  check thv0
!  Assumes th0, qv0 are accurate

      if(imoist.eq.1)then
        do k=kb,ke
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k)+qc0(i,j,k)+qi0(i,j,k))
        enddo
      else
        do k=kb,ke
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
          thv0(i,j,k)=th0(i,j,k)
        enddo
      endif

!----------------------------
!  calculate pressure, density, and temperature

      do k=kb,ke
        prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
        rho0(i,j,k)=prs0(i,j,k)/(rd*th0(i,j,k)*pi0(i,j,k)*(1.0+qv0(i,j,k)*reps))
        t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
      enddo

!----------------------------
!  This reduces errors associated with buoyancy term
!  (seems kind of redundant ... but it works)

      !  This qv0 must match specification of qva array in 
      !  the INIT3D subroutine
      !  (i.e., identical bit-for-bit calculation)

!    IF(imoist.eq.1)THEN
!      do k=kb,ke
!        qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k))
!      enddo
!    ENDIF


      !  This thv0 must exactly match the manner in which thv
      !  is calculated in the SOLVE subroutine
      !  (i.e., identical bit-for-bit calculation)

      do k=kb,ke
        if(imoist.eq.1)then
          thv0(i,j,k)=(th0(i,j,k)+0.0)*(1.0+reps*max(0.0,qv0(i,j,k)))  &
                                      /(1.0+max(0.0,qv0(i,j,k))+max(0.0,qc0(i,j,k)+qi0(i,j,k)))
        else
          thv0(i,j,k)=th0(i,j,k)
        endif
      enddo

    enddo
    enddo

      call bcs(rho0)
#ifdef MPI
      call comm_all_s(rho0,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                           n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
#endif
      call extrapbcs(rho0)

      IF( psolver.eq.5 .or. psolver.eq.6 )THEN   
        ! incompressible or compressible-Boussinesq:  set rho0 to a constant

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          rho0(i,j,k) = 1.0
        enddo
        enddo
        enddo

      ENDIF

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rr0(i,j,k)=1.0/rho0(i,j,k)
      enddo
      enddo
      enddo

      do k=2,nk
      do j=jb,je
      do i=ib,ie
        rf0(i,j,k)=c1(i,j,k)*rho0(i,j,k-1)+c2(i,j,k)*rho0(i,j,k)
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,1) = cgs1*rho0(i,j,1)+cgs2*rho0(i,j,2)+cgs3*rho0(i,j,3)
        rf0(i,j,0)=rf0(i,j,1)
        rho0s(i,j) = rf0(i,j,1)
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,nk+1) = cgt1*rho0(i,j,nk)+cgt2*rho0(i,j,nk-1)+cgt3*rho0(i,j,nk-2)
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rrf0(i,j,k)=1.0/rf0(i,j,k)
      enddo
      enddo
      enddo

!        i = 1
!        j = 1
!        print *
!        print *,'  nk = ',nk
!        k = 1
!        print *,k,rf0(i,j,k)
!      do k=2,nk+1
!        print *,k,rf0(i,j,k),rf0(i,j,k)-rf0(i,j,k-1)
!      enddo
!        print *
!      stop 11111

!-----------------------------------------------------------------------
!  values at surface:

      ! Get surface p/T/q (for surface models and for CAPE calculation):

      IF( psurf.ge.tsmall .and. tsurf.ge.tsmall .and. (.not.terrain_flag) )THEN
        ! this section of code only if no terrain
#ifdef MPI
        call MPI_BCAST(psurf,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(tsurf,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qsurf,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
        thsurf = tsurf/( (psurf*rp00)**rovcp )
        do j=jb,je
        do i=ib,ie
          prs0s(i,j) = psurf
          pi0s(i,j) = (psurf*rp00)**rovcp
          rth0s(i,j) = thsurf**(-1)
        enddo
        enddo
      ELSE
        do j=jb,je
        do i=ib,ie
          thsurf = cgs1*th0(i,j,1)+cgs2*th0(i,j,2)+cgs3*th0(i,j,3)
           qsurf = cgs1*qv0(i,j,1)+cgs2*qv0(i,j,2)+cgs3*qv0(i,j,3)
          thvsurf = thsurf*(1.0+qsurf*reps)/(1.0+qsurf)
          ! use hydrostatic equation:
          pi0s(i,j) = pi0(i,j,1)+zh(i,j,1)*g/(cp*0.5*(thvsurf+thv0(i,j,1)))
          prs0s(i,j) = p00*( pi0s(i,j)**(cp/rd) )
          rth0s(i,j) = thsurf**(-1)
        enddo
        enddo
      ENDIF

      i = 1
      j = 1
      if( psurf.lt.tsmall ) psurf = prs0s(i,j)
      if( tsurf.le.tsmall ) tsurf = cgs1*t0(i,j,1)+cgs2*t0(i,j,2)+cgs3*t0(i,j,3)
      if( qsurf.le.tsmall ) qsurf = cgs1*qv0(i,j,1)+cgs2*qv0(i,j,2)+cgs3*qv0(i,j,3)
      thsurf = tsurf/( (psurf*rp00)**rovcp )


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  Start definition of base state wind (iwnd option)  ccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!-----------------------------------------------------------------------
!  Get wind profiles  ...  assume zero wind to start

!--------------------------
! 061012:
! Ignore this section if isnd = 7;  in that case, wind profile has
! already been retrieved from the input_sounding file
!--------------------------

    IF(isnd.ne.7)THEN

      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        u0(i,j,k)= 0.0
      enddo
      enddo
      enddo

      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        v0(i,j,k)= 0.0
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!  iwnd = 1
!  RKW-type wind profile
!  reference: Rotunno, Klemp, and Weisman, 1988, JAS, 463-485.

      if(iwnd.eq.1)then

        udep1   =     0.0    ! height of bottom of shear layer (m)
        udep2   =  2500.0    ! height of top of shear layer (m)
        uconst1 =     0.0    ! u at bottom of shear layer
        uconst2 =    10.0    ! u at top of shear layer

        do k=1,nk
        do j=0,nj+1
        do i=0,ni+2
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.lt.udep1)then
            u0(i,j,k)=uconst1
          elseif(zu.gt.udep1 .and. zu.lt.udep2)then
            u0(i,j,k)=(uconst2-uconst1)*(zu-udep1)/(udep2-udep1)+uconst1
          else
            u0(i,j,k)=uconst2
          endif
        enddo
        enddo
        enddo

        do k=1,nk
        do j=0,nj+2
        do i=0,ni+1
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          v0(i,j,k)=0.0
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 2
!  Weisman-Klemp type supercell profile

      elseif(iwnd.eq.2)then

        udep1=2000.0
        udep2=6000.0
        umax1=7.0
        umax2=31.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.le.udep1)THEN
            ANGLE=90.0*(zu/udep1)*(pi/180.0)
            u0(i,j,k)=umax1-umax1*cos(ANGLE)
          elseif(zu.gt.udep1 .and. zu.le.udep2)THEN
            u0(i,j,k)=umax1+(zu-udep1)*(umax2-umax1)/(udep2-udep1)
          ELSE
            u0(i,j,k)=umax2
          ENDIF
        enddo
        enddo
        enddo

        vmax1=umax1

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          if(zv.le.udep1)THEN
            ANGLE=90.0*(zv/udep1)*(pi/180.0)
            v0(i,j,k)=vmax1*SIN(ANGLE)
          elseif(zv.gt.udep1 .and. zv.le.udep2)THEN
            v0(i,j,k)=vmax1
          ELSE
            v0(i,j,k)=vmax1
          ENDIF
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 3
!  Mulit-cell type profile (?)

      elseif(iwnd.eq.3)then

        udep1=0.0
        udep2=7500.0
        umax1=-40.0/pi
        umax2=40.0/pi+40.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.le.udep2)then
            u0(i,j,k)=umax1+(zu-udep1)*(umax2-umax1)/(udep2-udep1)
          else
            u0(i,j,k)=umax2
          endif
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          v0(i,j,k)=40.0/pi
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 4
!  Multi-cell
!  reference:  Weisman and Klemp, 1982, MWR, 110, 504-520.

      elseif(iwnd.eq.4)then

        umax1=35.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          u0(i,j,k)=umax1*tanh(0.5*(zh(i-1,j,k)+zh(i,j,k))/3000.0)
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          v0(i,j,k)=0.
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 5
!  reference:  Dornbrack et al., 2005, Atmos. Sci. Let., 6, 118-122

      elseif(iwnd.eq.5)then

        umax  =   15.0
        udep1 = 4000.0
        udep2 = 6000.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.lt.udep1)then
            u0(i,j,k) = umax
          elseif(zu.lt.udep2)then
            alpha=0.25*pi*(1.0+cos(pi*(zu-udep1)/(udep2-udep1)))
            u0(i,j,k) = umax*sin(alpha)
          else
            u0(i,j,k) = 0.0
          endif
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          if(zv.lt.udep1)then
            v0(i,j,k) = 0.0
          elseif(zv.lt.udep2)then
            alpha=0.25*pi*(1.0+cos(pi*(zv-udep1)/(udep2-udep1)))
            v0(i,j,k) = umax*cos(alpha)
          else
            v0(i,j,k) = umax
          endif
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 6
!  constant wind

      ELSEIF(iwnd.eq.6)THEN

        if( testcase.eq.1 )then
          ! Sullivan and Patton (2011, JAS)
          uval = 1.0
          vval = 0.0
        elseif( testcase.eq.2 )then
          ! Moeng and Sullivan (1994, JAS)
          uval = 15.0
          vval =  0.0
        elseif( testcase.eq.4 )then
          ! nonprecipitating stratocumulus case (Stevens et al. MWR 2005)            
          uval =  7.0
          vval = -5.5
        elseif( testcase.eq.9 )then
          ! stable boundary layer (Beare et al. 2006, BLM)
          uval =  8.0
          vval =  0.0
        elseif( testcase.eq.11 )then
          ! SAS
          uval =  2.0
          vval =  0.0
        elseif( testcase.eq.12 )then
          ! Martinuzzi and Tropea (1993) wind tunnel
          uval = 20.0
          vval =  0.0
        elseif( testcase.eq.14 )then
          ! shallow Cu over land (Brown et al. 2002)
          uval = 10.0
          vval =  0.0
        else
          ! default:
          uval = 10.0
          vval =  0.0
        endif

        do k=kb,ke
        do j=jb,je
        do i=ib,ie+1
          u0(i,j,k) = uval
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je+1
        do i=ib,ie
          v0(i,j,k) = vval
        enddo
        enddo
        enddo


!-----------------------------------------------------------------------
!  iwnd = 8
!  Wind speed decreasing with height
!  Simple hurricane boundary layer (Bryan et al 2017, BLM)

      ELSEIF(iwnd.eq.8)THEN

!!!        vmax = 40.0        ! this is now set in namelist.input
                              ! (see hurr_vg in param21 section)

        ! velocity in cylindrical coords:
        do k=1,nk

!!!          ! linearly decreasing with height:
!!!          ug(k) = 0.0
!!!          vg(k) = hurr_vg * max( 0.0 , (1.0-zh(1,1,k)/18000.0) )

          ! constant with height:
          ug(k) = 0.0
          vg(k) = hurr_vg

        enddo

        ! convert to Cartesian components:
        do k=1,nk
        do j=0,nj+2
        do i=0,ni+2
          u0(i,j,k) = - vg(k) * sin(-hurr_angle)
          v0(i,j,k) =   vg(k) * cos(-hurr_angle)
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 9
!  Siebesma et al (2003, JAS) shallow Cu case:

      ELSEIF(iwnd.eq.9)THEN

        udep1   =   700.0    ! height of bottom of shear layer (m)
        udep2   =  3000.0    ! height of top of shear layer (m)
        uconst1 =    -8.75   ! u at bottom of shear layer
        uconst2 =    -4.61   ! u at top of shear layer

        do k=1,nk

          do j=0,nj+1
          do i=0,ni+2
            zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
            if(zu.le.udep1)then
              u0(i,j,k)=uconst1
            elseif(zu.gt.udep1 .and. zu.lt.udep2)then
              u0(i,j,k)=(uconst2-uconst1)*(zu-udep1)/(udep2-udep1)+uconst1
            else
              u0(i,j,k)=uconst2
            endif
          enddo
          enddo

          do j=0,nj+2
          do i=0,ni+1
            zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
            v0(i,j,k)=0.0
          enddo
          enddo

        enddo


!-----------------------------------------------------------------------
!  iwnd = 10
!  Ackerman et al (2009, MWR) drizzling Sc case

      ELSEIF(iwnd.eq.10)THEN

        do k=1,nk

          do j=0,nj+1
          do i=0,ni+2
            zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
            u0(i,j,k) =  3.0 + 4.3*zu/1000.0
          enddo
          enddo

          do j=0,nj+2
          do i=0,ni+1
            zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
            v0(i,j,k) = -9.0 + 5.6*zv/1000.0
          enddo
          enddo

        enddo

!-----------------------------------------------------------------------
!  iwnd = 11
!   RICO, precipitating shallow Cu  (VanZanten et al 2011, JAMES)

      ELSEIF(iwnd.eq.11)THEN

        udep1   =      0.0   ! height of bottom of shear layer (m)
        udep2   =   4000.0   ! height of top of shear layer (m)
        uconst1 =     -9.9   ! u at bottom of shear layer
        uconst2 =     -1.9   ! u at top of shear layer

        do k=1,nk

          do j=0,nj+1
          do i=0,ni+2
            zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
            if(zu.le.udep1)then
              u0(i,j,k)=uconst1
            elseif(zu.gt.udep1 .and. zu.lt.udep2)then
              u0(i,j,k)=(uconst2-uconst1)*(zu-udep1)/(udep2-udep1)+uconst1
            else
              u0(i,j,k)=uconst2
            endif
          enddo
          enddo

          do j=0,nj+2
          do i=0,ni+1
            zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
            v0(i,j,k) = -3.8
          enddo
          enddo

        enddo

!-----------------------------------------------------------------------

      ENDIF    ! endif for iwnd options

    ENDIF   ! endif for isnd=7 check


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!-----------------------------------------------------------------------
!  subtract off umove and vmove (if applicable)

      if(imove.eq.1)then
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          u0(i,j,k)=u0(i,j,k)-umove
          v0(i,j,k)=v0(i,j,k)-vmove
        enddo
        enddo
        enddo
      else
        umove=0.0
        vmove=0.0
      endif

!-----------------------------------------------------------------------
!  Fill in ghost cells

      call bcu(u0)
      call bcv(v0)

      !--------
#ifdef MPI
      call comm_3u_start(u0,uw31,uw32,ue31,ue32,   &
                            us31,us32,un31,un32,reqs_u)
      call comm_3u_end(u0,uw31,uw32,ue31,ue32,   &
                          us31,us32,un31,un32,reqs_u)
      call getcorneru3(u0,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                          s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      call bcu2(u0)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=jb,je
      do i=ib,ie+1
        u0(i,j,0)    = cgs1*u0(i,j,1)+cgs2*u0(i,j,2)+cgs3*u0(i,j,3)
        u0(i,j,nk+1) = cgt1*u0(i,j,nk)+cgt2*u0(i,j,nk-1)+cgt3*u0(i,j,nk-2)
      enddo
      enddo
      !--------
#ifdef MPI
      call comm_3v_start(v0,vw31,vw32,ve31,ve32,   &
                            vs31,vs32,vn31,vn32,reqs_v)
      call comm_3v_end(v0,vw31,vw32,ve31,ve32,   &
                          vs31,vs32,vn31,vn32,reqs_v)
      call getcornerv3(v0,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                          s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      call bcv2(v0)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=jb,je+1
      do i=ib,ie
        v0(i,j,0)    = cgs1*v0(i,j,1)+cgs2*v0(i,j,2)+cgs3*v0(i,j,3)
        v0(i,j,nk+1) = cgt1*v0(i,j,nk)+cgt2*v0(i,j,nk-1)+cgt3*v0(i,j,nk-2)
      enddo
      enddo
      !--------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  End definition of base state wind (iwnd option)  ccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  forcing tendencies and settings for certain configurations:
!
    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase .eq. 3 )THEN

      ! shallow Cu case  (Siebesma et al, 2003, JAS)

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  k,zf,w:'

      wmin  =  -0.0065
      z1    =   1500.0
      z2    =   2100.0

      do k=1,nk
        ! note:  wprof is on w levels (ie, full levels).   Use zf.
        if( zf(1,1,k).le.z1 )then
          wprof(k) = 0.0 + ( wmin - 0.0 )*(zf(1,1,k)-0.0)/(z1-0.0)
        elseif( zf(1,1,k).le.z2 )then
          wprof(k) = wmin + ( 0.0 - wmin )*(zf(1,1,k)-z1)/(z2-z1)
        else
          wprof(k) = 0.0
        endif
        if( myid.eq.0 ) print *,k,zf(1,1,k),wprof(k)
      enddo
 
      radsfc = -2.0 /( 3600.0 * 24.0 )
      z1     =  1500.0
      z2     =  2100.0

      do k=1,nk
        if( zh(1,1,k).le.z1 )then
          thfrc(k) = radsfc
        elseif( zh(1,1,k).le.z2 )then
          thfrc(k) = radsfc + ( 0.0 - radsfc )*(zh(1,1,k)-z1)/(z2-z1)
        else
          thfrc(k) = 0.0
        endif
      enddo

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  k,zh,thfrc,qvfrc '

      qvsfc  =  -1.2e-8
      z1     =   300.0
      z2     =   500.0

      do k=1,nk
        if( zh(1,1,k).le.z1 )then
          qvfrc(k) = -1.2e-8
        elseif( zh(1,1,k).le.z2 )then
          qvfrc(k) = qvsfc + ( 0.0 - qvsfc )*(zh(1,1,k)-z1)/(z2-z1)
        else
          qvfrc(k) = 0.0
        endif
        if( myid.eq.0 ) print *,k,zh(1,1,k),thfrc(k)*3600.0*24.0,qvfrc(k)
      enddo


      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  k,ug,vg '
      do k=1,nk
        ug(k) = -10.0 + 0.0018 * zh(1,1,k)
        vg(k) =   0.0
      enddo
      if( myid.eq.0 ) print *

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase .eq. 4 )THEN

      ! nonprecipitating stratoCu (Stevens et al, 2005, MWR)

      do k=2,nk
        ! note:  wprof is on w levels (ie, full levels):   Use zf.
        wprof(k) = -3.75e-6 * zf(1,1,k)
      enddo

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase .eq. 5 )THEN

      ! drizzling stratoCu (Ackerman et al, 2009, MWR)

      do k=2,nk
        ! note:  wprof is on w levels (ie, full levels):   Use zf.
        wprof(k) = -3.75e-6 * zf(1,1,k)
      enddo

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase.eq.6 .or. testcase.eq.10 .or. testcase.eq.15 )THEN

      ! hurricane boundary layer (Bryan et al, 2017, BLM)

!!!      hurr_rad = 40000.0      ! radius from center of hurricane (m) 
                                 ! (cm1r21.0: now set this in namelist.input)

      wmin  =   -0.01    ! value of w at z = z1
!!!      wmin  =   var1     ! value of w at z = z1
      z1    =  1000.0    ! level of min w
      z2    =  4000.0    ! level where w = 0

      do k=1,nk

        ! radial gradient of tangential velocity (defined on scalar levels)
!!!        dvdr(k) = -0.00085*max(0.0,1.0-zh(1,1,k)/18000.0)
        ! set this in namelist.input:
        ! hurr_vgpl is the power law of the gradient wind
        dvdr(k) = hurr_vgpl * ( vg(k) / hurr_rad )

        ! subsidence (defined on w levels)
        if( zf(1,1,k).le.z1 )then
          wprof(k) = 0.0 + ( wmin - 0.0 )*(zf(1,1,k)-0.0)/(z1-0.0)
        elseif( zf(1,1,k).le.z2 )then
          wprof(k) = wmin + ( 0.0 - wmin )*(zf(1,1,k)-z1)/(z2-z1)
        else
          wprof(k) = 0.0
        endif

      enddo

      if( dowr ) write(outfile,*)
      if( dowr ) write(outfile,*) '  hurr_rad (m) = ',hurr_rad
      if( dowr ) write(outfile,*)
      if( dowr ) write(outfile,*) '      zh(m),vg(m/s),dvdr(1/s),n = '
      do k=1,nk
        if( dowr ) write(outfile,*) zh(1,1,k),vg(k),dvdr(k),-dvdr(k)*hurr_rad/max(1.0e-12,vg(k))
      enddo
      if( dowr ) write(outfile,*)

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase .eq. 7 )THEN

      !   RICO, precipitating shallow Cu  
      !   (VanZanten et al 2011, JAMES)

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  k,zf,w:'

      wmin  =  -0.005
      z1    =   2260.0

      do k=1,nk
        ! note:  wprof is on w levels (ie, full levels).   Use zf.
        if( zf(1,1,k).le.z1 )then
          wprof(k) = 0.0 + ( wmin - 0.0 )*(zf(1,1,k)-0.0)/(z1-0.0)
        else
          wprof(k) = wmin
        endif
        if( myid.eq.0 ) print *,k,zf(1,1,k),wprof(k)
      enddo
 
      thsfc  = -2.5 /( 3600.0 * 24.0 )

      do k=1,nk
        thfrc(k) = thsfc
      enddo

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  k,zh,thfrc,qvfrc '

      qv1    =    -1.0 / 1000.0 / ( 3600.0 * 24.0 )
      qv2    =  0.3456 / 1000.0 / ( 3600.0 * 24.0 )
      z1     =  2980.0

      do k=1,nk
        ! humidity tendency:
        if( zh(1,1,k).le.z1 )then
          qvfrc(k) = qv1 + ( qv2 - qv1 )*(zh(1,1,k)-0.0)/(z1-0.0)
        else
          qvfrc(k) = qv2
        endif
        if( myid.eq.0 ) print *,k,zh(1,1,k),thfrc(k)*3600.0*24.0,qvfrc(k)
      enddo


      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  k,ug,vg '
      do k=1,nk
        ug(k) = u0(1,1,k)
        vg(k) = v0(1,1,k)
      enddo
      if( myid.eq.0 ) print *

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase .eq. 9 .or. testcase.eq.14 )THEN


    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase.eq.11 )THEN
      ! SAS NCAR PBL intercomparison

!!!      ug = 2.0
!!!      vg = 0.0
      ug = u0(1,1,1)+umove
      vg = v0(1,1,1)+vmove

      do k=2,nk
        ! note:  wprof is on w levels (ie, full levels):   Use zf.
        wprof(k) = -9.0e-6 * zf(1,1,k)
      enddo

      ufrc = 0.0
      vfrc = 0.0
!!!      thfrc = 6.4e-4
!!!      qvfrc = 1.5e-4 / 2.5e6

      thfrc = 0.0
      qvfrc = 0.0
      wprof = 0.0

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase.eq.12 )THEN

      ! Martinuzzi and Tropea (1993, JFE) wind tunnel case:

        ulspg = 36.0
        vlspg = 0.0

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    IF( testcase.eq.13 )THEN

      ! UF Boundary Layer Wind Tunnel case:

      ! BLWT:
      if( isnd.ne.7 )then
        u0 = 18.0
        v0 =  0.0
        ug = 18.0
        vg =  0.0
      endif
      ulspg = 0.000545
      vlspg = 0.0

    ENDIF

    !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

    do k=kb,ke
      ug(k) = ug(k) - umove
      vg(k) = vg(k) - vmove
    enddo


    ! boundary conditions:
    wprof(1) = 0.0
    wprof(nk+1) = 0.0


    if( lspgrad.eq.4 )then
    if( myid.eq.0 )then
      print *
      print *,'  zh,ulspg,vlspg = '
      do k=1,nk
        print *,k,zh(1,1,k),ulspg(k),vlspg(k)
      enddo
      print *
    endif
    endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do k=1,nk
        thavg(k) = th0(1,1,k)
        pavg(k) = prs0(1,1,k)
      enddo

      if( imoist.eq.1 )then
        do k=1,nk
          qavg(k,nqv) = qv0(1,1,k)
        enddo
      endif


!---------------------------------------------------------------------

!  Print out base state

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,101)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(dowr) write(outfile,103) k,zh(1,1,k),prs0(1,1,k),pi0(1,1,k),rho0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

101   format(7x,'k      zh (m)    prs0 (Pa)        pi0        rho0 (kg/m^3)')
102   format(4x,'---------------------------------------------------------------------')
103   format(4x,i4,4x,f8.2,4x,f9.2,4x,f10.7,4x,f10.7)

!-----

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,104)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(dowr) write(outfile,105) k,zh(1,1,k),th0(1,1,k),thv0(1,1,k),t0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

104   format(7x,'k      zh (m)    th0 (K)     thv0 (K)     t0 (K) ')
105   format(4x,i4,4x,f8.2,4x,f8.4,4x,f8.4,4x,f8.4)

!-----

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,114)
      if(dowr) write(outfile,102)
      do k=2,nk
        if(dowr) write(outfile,115) k,zf(1,1,k),  &
                   -1000*(t0(1,1,k)-t0(1,1,k-1))*rdz*mf(1,1,k),  &
                    1.0e4*g*alog(thv0(1,1,k)/thv0(1,1,k-1))*rdz*mf(1,1,k)
      enddo
      if(dowr) write(outfile,*)

114   format(7x,'k      zf (m)     l.r. (K/km)     N^2')
115   format(4x,i4,4x,f8.2,4x,f10.4,4x,f8.4)

!-----
! Ri check ... warn user if Ri < 0.25, but do not stop model
! 061021: added N^2 check ... stop model if N^2 < 0

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Richardson number:  NOTE!!!  Ri should be > 0.25 for most applications!'
      if(dowr) write(outfile,*)
      do k=2,nk
        nm = g*alog(thv0(1,1,k)/thv0(1,1,k-1))*rdz*mf(1,1,k)
        if(nm .lt. -1.0e-12)then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  k,zf,N^2:',k,zf(1,1,k),nm
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' Warning.  N^2 (Brunt-Vaisala frequency squared) is less than zero!'
          if(dowr) write(outfile,*) ' This is really, really not recommended for the base state.'
!!!          if(dowr) write(outfile,*) ' Stopping model ....'
          if(dowr) write(outfile,*)
!!!          call stopcm1
        endif
        dudz = (u0(1,1,k)-u0(1,1,k-1))*rdz*mf(1,1,k)
        dvdz = (v0(1,1,k)-v0(1,1,k-1))*rdz*mf(1,1,k)
        rinum = nm/(1.0e-12+dudz*dudz+dvdz*dvdz)
        if(rinum.gt.0.25)then
          if(dowr) write(outfile,*) '  k,zf,Ri:',k,zf(1,1,k),rinum
        else
          if(dowr) write(outfile,*) '  k,zf,Ri:',k,zf(1,1,k),rinum,'<---- NOTE!  Ri < 0.25'
        endif
      enddo
      if(dowr) write(outfile,*)

!-----

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,106)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(rh0(1,1,k).gt.0.999 .or. rh0(1,1,k).lt.0.0001 .or. imoist.eq.0)then
          tlcl=t0(1,1,k)
        else
          tlcl=55.0+(2840./(3.5*alog(t0(1,1,k))-    &
                   alog(0.01*prs0(1,1,k)*qv0(1,1,k)/(0.622+qv0(1,1,k)))-4.805))
        endif
        if(dowr) write(outfile,107) k,zh(1,1,k),                                 &
            th0(1,1,k)*exp((3376./tlcl-2.54)*qv0(1,1,k)*(1.0+0.81*qv0(1,1,k))),  &
            rh0(1,1,k),qv0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

106   format(7x,'k      zh (m)      theta-e (K)       rh0          qv0 ')
107   format(4x,i4,4x,f8.2,4x,4x,f8.4,4x,f9.6,4x,f9.6)

!-----

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,126)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(rh0(1,1,k).gt.0.999 .or. rh0(1,1,k).lt.0.0001 .or. imoist.eq.0)then
          tlcl=t0(1,1,k)
        else
          tlcl=55.0+(2840./(3.5*alog(t0(1,1,k))-    &
                   alog(0.01*prs0(1,1,k)*qv0(1,1,k)/(0.622+qv0(1,1,k)))-4.805))
        endif
        if(dowr) write(outfile,127) k,zh(1,1,k),qc0(1,1,k),qi0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

126   format(7x,'k      zh (m)        qc0          qi0')
127   format(4x,i4,4x,f8.2,4x,f9.6,4x,f9.6)

!------------------------------------------------------------------
!  Get CAPE,CIN,etc:

  IF(imoist.eq.1)THEN

    allocate(  pfoo(nk+1) )
    allocate(  tfoo(nk+1) )
    allocate( qvfoo(nk+1) )

    pfoo(1) = psurf * 0.01
    tfoo(1) = tsurf - 273.15
   qvfoo(1) = max( qsurf , 1.0e-20 )

    do k=1,nk
      pfoo(k+1) = 0.01*prs0(1,1,k)
      tfoo(k+1) = t0(1,1,k) - 273.15
      qvfoo(k+1) = max( qv0(1,1,k) , 1.0e-20 )
    enddo

    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) '  Thermodynamic properties of base-state sounding:'
  IF(terrain_flag)THEN
    if(dowr) write(outfile,*) '    (for lower-left corner of domain:  i=1,j=1)  '
  ENDIF
    if(dowr) write(outfile,*)
    do n=1,3
      call getcape( n , nk+1 , pfoo , tfoo , qvfoo , cape , cin ,   &
                    zlcl, zlfc, zel , psource , tsource , qvsource )
      if(n.eq.1)then
        if(dowr) write(outfile,*) '    for surface parcel:'
      elseif(n.eq.2)then
        if(dowr) write(outfile,*) '    for most-unstable parcel:'
      elseif(n.eq.3)then
        if(dowr) write(outfile,*) '    for mixed-layer parcel:'
      endif
      if(dowr) write(outfile,116) 0.01*psource,tsource,1000.0*qvsource
      if(dowr) write(outfile,118) zlcl,zlfc,zel
      if(dowr) write(outfile,117) cape,cin
116   format('        source p(mb),T(K),qv(g/kg) = ',3(4x,f6.1))
118   format('        LCL,LFC,EL (m AGL)         = ',3(3x,f7.1))
117   format('        CAPE,CIN (J/kg)            = ',2(4x,f6.1))
      if(dowr) write(outfile,*)
    enddo
    if(dowr) write(outfile,*)

    deallocate(  pfoo )
    deallocate(  tfoo )
    deallocate( qvfoo )

  ENDIF

!------------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,108)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(dowr) write(outfile,109) k,zh(1,1,k),u0(1,1,k),v0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

108   format(7x,'k      zh (m)     u0 (m/s)     v0 (m/s)')
109   format(4x,i4,4x,f8.2,4x,f9.4,4x,f9.4)

!------------------------------------------------------------------

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          rth0(i,j,k)=1.0/th0(i,j,k)
        enddo
        enddo
        enddo

        do k=kbb2,kbe2
        do j=jbb2,jbe2
        do i=ibb2,ibe2
          thrd(i,j,k)=th0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------

      pmin = 1.0e30

      do k=1,nk
      do j=1,nj
      do i=1,ni
        pmin = min( pmin , prs0(i,j,k) )
      enddo
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(mpi_in_place,pmin,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
#endif

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  pmin = ',pmin
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  psurf,tsurf  = ',psurf,tsurf
      if(dowr) write(outfile,*) '  qsurf,thsurf = ',qsurf,thsurf
      if(dowr) write(outfile,*)

      dothis = .false.
      IF( dothis )THEN
      IF( myid.eq.0 )THEN

        open(unit=21,file='input_sounding_grid')
        i = 1
        j = 1
        write(21,123) 0.01*psurf,thsurf,1000.0*qsurf
123     format(2x,f11.4,2x,f11.4,2x,f10.6)
124     format(2x,f11.4,2x,f11.4,2x,f10.6,2x,f11.4,2x,f11.4)
        do k=1,nk
          write(21,124) zh(i,j,k),th0(i,j,k),1000.0*qv0(i,j,k),u0(i,j,k),v0(i,j,k)
        enddo
          write(21,124) cgt1*zh(i,j,nk)+cgt2*zh(i,j,nk-1)+cgt3*zh(i,j,nk-2),     &
                        cgt1*th0(i,j,nk)+cgt2*th0(i,j,nk-1)+cgt3*th0(i,j,nk-2),  &
                1000.0*(cgt1*qv0(i,j,nk)+cgt2*qv0(i,j,nk-1)+cgt3*qv0(i,j,nk-2)), &
                        cgt1*u0(i,j,nk)+cgt2*u0(i,j,nk-1)+cgt3*u0(i,j,nk-2),     &
                        cgt1*v0(i,j,nk)+cgt2*v0(i,j,nk-1)+cgt3*v0(i,j,nk-2)
        close(unit=21)

      ENDIF
      ENDIF

!--------------------------------------------------------------

      dothis = .false.
      IF( dothis )THEN
      IF( myid.eq.0 )THEN

        open(unit=21,file='lsnudge_0001.dat')
        write(21,*) '*  Header:   lsnudge_time1 (s)   lsnudge_time2 (s)'
        write(21,*) '                    0.0                1.0e30 '
        write(21,*) '*  Profile.   Note: values will be ignored if nudging for that variable is off.'
        write(21,*) '*    z (m)   theta (K)    qv (g/kg)  u (m/s, grnd-reltv) v (m/s, grnd-reltv)'
        i = 1
        j = 1
!!!        write(21,123) 0.01*psurf,thsurf,1000.0*qsurf
        do k=1,nk
          write(21,124) zh(i,j,k),th0(i,j,k),1000.0*qv0(i,j,k),u0(i,j,k),v0(i,j,k)
        enddo
          write(21,124) cgt1*zh(i,j,nk)+cgt2*zh(i,j,nk-1)+cgt3*zh(i,j,nk-2),     &
                        cgt1*th0(i,j,nk)+cgt2*th0(i,j,nk-1)+cgt3*th0(i,j,nk-2),  &
                1000.0*(cgt1*qv0(i,j,nk)+cgt2*qv0(i,j,nk-1)+cgt3*qv0(i,j,nk-2)), &
                        cgt1*u0(i,j,nk)+cgt2*u0(i,j,nk-1)+cgt3*u0(i,j,nk-2),     &
                        cgt1*v0(i,j,nk)+cgt2*v0(i,j,nk-1)+cgt3*v0(i,j,nk-2)
        close(unit=21)

      ENDIF
      ENDIF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      base_pbot = cgs1*prs0(1,1,1)+cgs2*prs0(1,1,2)+cgs3*prs0(1,1,3)
      base_ptop = cgt1*prs0(1,1,nk)+cgt2*prs0(1,1,nk-1)+cgt3*prs0(1,1,nk-2)

      base_pibot = cgs1*pi0(1,1,1)+cgs2*pi0(1,1,2)+cgs3*pi0(1,1,3)
      base_pitop = cgt1*pi0(1,1,nk)+cgt2*pi0(1,1,nk-1)+cgt3*pi0(1,1,nk-2)

      base_tbot = cgs1*t0(1,1,1)+cgs2*t0(1,1,2)+cgs3*t0(1,1,3)
      base_ttop = cgt1*t0(1,1,nk)+cgt2*t0(1,1,nk-1)+cgt3*t0(1,1,nk-2)

      base_thbot = cgs1*th0(1,1,1)+cgs2*th0(1,1,2)+cgs3*th0(1,1,3)
      base_thtop = cgt1*th0(1,1,nk)+cgt2*th0(1,1,nk-1)+cgt3*th0(1,1,nk-2)

      base_qvbot = cgs1*qv0(1,1,1)+cgs2*qv0(1,1,2)+cgs3*qv0(1,1,3)
      base_qvtop = cgt1*qv0(1,1,nk)+cgt2*qv0(1,1,nk-1)+cgt3*qv0(1,1,nk-2)

      if( myid.eq.0 )then
        print *
        print *,'  base_pbot, base_ptop  = ',base_pbot,base_ptop
        print *,'  base_pibot,base_pitop = ',base_pibot,base_pitop
        print *,'  base_tbot, base_ttop  = ',base_tbot,base_ttop
        print *,'  base_thbot,base_thtop = ',base_thbot,base_thtop
        print *,'  base_qvbot,base_qvtop = ',base_qvbot,base_qvtop
        print *
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(dowr) write(outfile,*) 'Leaving BASE'

      return

    !--------------------------------------------------------------

9008  print *
      print *,'  9008: error reading input_sounding file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9009  print *
      print *,'  9009: error opening input_sounding file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

    !--------------------------------------------------------------

      return
      end subroutine base

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    real function getthe(ttype,p,t,qv,qt)
    use constants
    implicit none
    integer ttype
    real p,t,qv,qt

    ! Assumes air is saturated

    real :: pd,cpm,lhv

    getthe = 0.0

  IF(ttype.eq.1)THEN
    pd=p/(1.0+reps*qv)
    cpm=cp+cpl*qt
    lhv=lv1-lv2*t
    getthe=t*((p00/pd)**(rd/cpm))*exp(lhv*qv/(cpm*t))
  ELSEIF(ttype.eq.2)THEN
    getthe=t*( (p00/p)**(0.2854*(1.0-0.28*qv)) )   &
            *exp( ((3376.0/t)-2.54)*qv*(1.0+0.81*qv) )
  ENDIF
    end function getthe

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine revthe(ttype,the,p,qt,t,qv)
    use constants
    use cm1libs , only : rslf
    implicit none
    integer ttype
    real the,p,qt,rh,t,qv

    ! Assumes air is saturated
    ! Input:   the,p,qt  (t is first guess value upon input)
    ! Output:  t,qv

    integer n
    real tlast,tinc,pd,cpm,lhv,thx,diff

    n=1
    tlast=t
    tinc=0.0
250 continue
      t=tlast+tinc
      qv=rslf(p,t)
    IF(ttype.eq.1)THEN
      pd=p/(1.0+reps*qv)
      cpm=cp+cpl*qt
      lhv=lv1-lv2*t
      thx=t*((p00/pd)**(rd/cpm))*exp(lhv*qv/(cpm*t))
    ELSEIF(ttype.eq.2)THEN
      thx=t*( (p00/p)**(0.2854*(1.0-0.28*qv)) )   &
           *exp( ((3376.0/t)-2.54)*qv*(1.0+0.81*qv) )
    ENDIF

      diff=the-thx
      if(n.ge.40) print *,n,p,the,thx
      if(abs(diff).gt.0.0001 .and. n.lt.50)then
        n=n+1
        tinc=0.30*diff
        tlast=t
        if(abs(tinc).ge.0.0001) go to 250
      elseif(n.ge.50)then
        print *,'n exceeded 50!'
        print *,'n=',n
        print *,'  1222 '
        call stopcm1
      endif

    return
    end subroutine revthe

  END MODULE base_module
