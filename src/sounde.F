  MODULE sounde_module

  implicit none

  private
  public :: sounde

  CONTAINS

      subroutine sounde(dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,zf,      &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,thv0,zs,           &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,dtu,dtu0,dtv,dtv0,    &
                        dum1,dum2,dum3,dum4,dum5,dum6,                    &
                        ppd ,fpk ,qk ,pk1,pk2,ftk,                        &
                        u0,rru,ua,u3d,uten,                               &
                        v0,rrv,va,v3d,vten,                               &
                        rrw,wa,w3d,wten,                                  &
                        ppi,pp3d,piadv,ppten,ppx,                         &
                        thv,ppterm,nrk,dttmp,rtime,mtime,get_time_avg,    &
                        bndy,kbdy,                                        &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      use input
      use constants
      use misclibs , only : convinitu,convinitv,get_wnudge
      use bc_module
      use comm_module
      use ib_module
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,rr0,rf0,rrf0,th0,rth0,thv0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je) :: dtu,dtv
      real, intent(in),    dimension(ib:ie,jb:je) :: dtu0,dtv0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppd,fpk,qk,pk1,pk2,ftk
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,ua,u3d,uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,va,v3d,vten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,piadv,ppten,ppx
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: thv,ppterm
      integer, intent(in) :: nrk
      real, intent(in)  :: dttmp,rtime
      double precision, intent(in) :: mtime
      logical, intent(in) :: get_time_avg
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      integer, intent(inout), dimension(rmp) :: reqs_p

!-----

      integer :: i,j,k,n,nloop
      real :: tem,tem1,tem2,tem3,tem4,r1,r2,rr1,rr2,dts,delp

      real :: temx,temy,u1,u2,v1,v2,w1,w2,ww,div,tavg,c1a,c1b,c2a,c2b,p1,p2,p3,p4

!---------------------------------------------------------------------

    IF( nrkmax.eq.3 )THEN
      if(nrk.eq.1)then
!!!        nloop=1
!!!        dts=dt/3.
        nloop=nint(float(nsound)/3.0)
        dts=dt/(nloop*3.0)
        if( dts.gt.(dt/nsound) )then
          nloop=nloop+1
          dts=dt/(nloop*3.0)
        endif
      elseif(nrk.eq.2)then
        nloop=0.5*nsound
        dts=dt/nsound
      elseif(nrk.eq.3)then
        nloop=nsound
        dts=dt/nsound
      endif
    ELSE
      stop 97393
    ENDIF

!!!      print *,'  nloop,dts,dttmp = ',nloop,dts,nloop*dts

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        dtu(i,j) = dts*dtu0(i,j)
        dtv(i,j) = dts*dtv0(i,j)
      enddo
      enddo

!---------------------------------------------------------------------
!  Arrays for vadv:

    IF( psolver.eq.2 )THEN

      IF(.not.terrain_flag)THEN

        ! without terrain:
        ! "s" velocities ARE NOT coupled with reference density
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem,tem1,r1,r2,c1a,c1b,c2a,c2b,p1,p2,p3,p4)
        do k=1,nk
          p1 =  dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k+1)*(1.0-c1(1,1,k+1))
          p2 = -dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k+1)*c2(1,1,k+1)
          p3 =  dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k)*c1(1,1,k)
          p4 =  dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k)*(c2(1,1,k)-1.0)
          do j=1,nj
          do i=1,ni
            pk2(i,j,k) = p2*piadv(i,j,k+1)+p1*piadv(i,j,k)
            pk1(i,j,k) = p3*piadv(i,j,k-1)+p4*piadv(i,j,k)
          enddo
          enddo
        enddo

      ELSE

        ! with terrain:
        ! "s" velocities ARE coupled with reference density
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem,tem1,tem2,r1,r2)
        do k=1,nk
          do j=1,nj
          do i=1,ni
            tem2 = dts*gz(i,j)*rdsf(k)*rr0(i,j,k)
            pk2(i,j,k) = tem2*rf0(i,j,k+1)*( -c2(i,j,k+1)*piadv(i,j,k+1)+(1.0-c1(i,j,k+1))*piadv(i,j,k) )
            pk1(i,j,k) = tem2*rf0(i,j,k+1)*( +c1(i,j,k  )*piadv(i,j,k-1)+(c2(i,j,k  )-1.0)*piadv(i,j,k) )
          enddo
          enddo
          IF( k.eq.1 )THEN
            do j=1,nj
            do i=1,ni
              dum3(i,j,1)=0.0
              dum3(i,j,nk+1)=0.0
            enddo
            enddo
          ENDIF
        enddo

      ENDIF

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        pk1(i,j,1) = 0.0
        pk2(i,j,nk) = 0.0
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------
!  Prepare for acoustic steps

      if( nrk.eq.1 )then

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          ppd(i,j,k)=ppx(i,j,k)
          rru(i,j,k)=0.0
          rrv(i,j,k)=0.0
          rrw(i,j,k)=0.0
        enddo
        enddo
        enddo

      else

        IF( axisymm.eq.0 )THEN

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            u3d(i,j,k)=ua(i,j,k)
            rru(i,j,k)=0.0
            v3d(i,j,k)=va(i,j,k)
            rrv(i,j,k)=0.0
            w3d(i,j,k)=wa(i,j,k)
            rrw(i,j,k)=0.0
            pp3d(i,j,k)=ppi(i,j,k)
            ppd(i,j,k)=ppx(i,j,k)
          enddo
          enddo
          enddo

        ELSEIF( axisymm.eq.1 )THEN

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            u3d(i,j,k)=ua(i,j,k)
            rru(i,j,k)=0.0
            w3d(i,j,k)=wa(i,j,k)
            rrw(i,j,k)=0.0
            pp3d(i,j,k)=ppi(i,j,k)
            ppd(i,j,k)=ppx(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

      endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Begin small steps:

      small_step_loop:  DO N=1,NLOOP

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-----

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,u3d)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,v3d)
 
        endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Open boundary conditions:

        IF(wbc.eq.2.and.ibw.eq.1)THEN
          ! west open bc tendency:
          call   ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dts)
        ENDIF
        IF(ebc.eq.2.and.ibe.eq.1)THEN
          ! east open bc tendency:
          call   ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dts)
        ENDIF

        IF(roflux.eq.1)THEN
          call restrict_openbc_we(rvh,rmh,rho0,u3d)
        ENDIF

!-----

      IF(axisymm.eq.0)THEN
        IF(sbc.eq.2.and.ibs.eq.1)THEN
          ! south open bc tendency:
          call   ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dts)
        ENDIF
        IF(nbc.eq.2.and.ibn.eq.1)THEN
          ! north open bc tendency:
          call   ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dts)
        ENDIF

        IF(roflux.eq.1)THEN
          call restrict_openbc_sn(ruh,rmh,rho0,v3d)
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  integrate u,v forward in time:

#ifdef MPI
        if( n.ne.1 )then
          if(timestats.ge.1) time_sound=time_sound+mytime()
          call comm_1p_end1(ppd,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        endif
#endif

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
!  psolver = 2

  ps27a:  &
  IF(psolver.eq.2)THEN

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid without terrain:

        tem1 = rdx*cp*0.5
        tem2 = rdy*cp*0.5

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            u3d(i,j,k)=u3d(i,j,k)+dtu(i,j)*( uten(i,j,k)     &
                       -tem1*uf(i)*(ppd(i,j,k)-ppd(i-1,j,k)) &
                                  *(thv(i,j,k)+thv(i-1,j,k)) )
            v3d(i,j,k)=v3d(i,j,k)+dtv(i,j)*( vten(i,j,k)     &
                       -tem2*vf(j)*(ppd(i,j,k)-ppd(i,j-1,k)) &
                                  *(thv(i,j,k)+thv(i,j-1,k)) )
          enddo
          enddo
        enddo

        IF( do_ib )THEN
          call zero_out_uv(bndy,kbdy,u3d,v3d)
        ENDIF

      ELSE
        ! axisymmetric grid:

        tem1 = rdx*cp*0.5

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dtu(i,j)*( uten(i,j,k)    &
                   -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)    &
                        *(thv(i,j,k)+thv(i-1,j,k)) )
          enddo
          enddo
        enddo

      ENDIF

    ELSE

        ! Cartesian grid with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=0,nj+1
          do k=2,nk
          do i=0,ni+1
            dum1(i,j,k) = (ppd(i,j,k)-ppd(i,j,k-1))*rds(k)
          enddo
          enddo
          do i=0,ni+1
            dum1(i,j,1) = 0.0
            dum1(i,j,nk+1) = 0.0
          enddo
        enddo

        tem = cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          ! x-dir
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                   -tem*(thv(i,j,k)+thv(i-1,j,k))*(         &
                     (ppd(i,j,k)-ppd(i-1,j,k))*rdx*uf(i)    &
              +0.125*( (dum1(i,j,k+1)+dum1(i-1,j,k+1))      &
                      +(dum1(i,j,k  )+dum1(i-1,j,k  )) )    &
                    *(gxu(i,j,k)+gxu(i,j,k+1))    ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)         &
                   -tem*(thv(i,j,k)+thv(i,j-1,k))*(         &
                     (ppd(i,j,k)-ppd(i,j-1,k))*rdy*vf(j)    &
              +0.125*( (dum1(i,j,k+1)+dum1(i,j-1,k+1))      &
                      +(dum1(i,j,k  )+dum1(i,j-1,k  )) )    &
                    *(gyv(i,j,k)+gyv(i,j,k+1))    ) )
          enddo
          enddo
        enddo


    ENDIF

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
!  psolver = 7

  ELSEIF(psolver.eq.7)THEN  ps27a

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem1,tem2)
        do k=1,nk
          tem1 = rdx*cp*thv0(1,1,k)
          tem2 = rdy*cp*thv0(1,1,k)
          do j=1,nj+1
          do i=1,ni+1
            u3d(i,j,k)=u3d(i,j,k)+dtu(i,j)*( uten(i,j,k)      &
                    -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i) )
            v3d(i,j,k)=v3d(i,j,k)+dtv(i,j)*( vten(i,j,k)      &
                    -tem2*(ppd(i,j,k)-ppd(i,j-1,k))*vf(j) )
          enddo
          enddo
        enddo

        IF( do_ib )THEN
          call zero_out_uv(bndy,kbdy,u3d,v3d)
        ENDIF

      ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem1)
        do k=1,nk
          tem1 = rdx*cp*thv0(1,1,k)
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)      &
                    -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i) )
          enddo
          enddo
        enddo

      ENDIF

    ELSE

        ! Cartesian grid with terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,r1,r2)
        do j=0,nj+1
          do k=2,nk
          do i=0,ni+1
            dum1(i,j,k) = (ppd(i,j,k)-ppd(i,j,k-1))*rds(k)
          enddo
          enddo
          do i=0,ni+1
            dum1(i,j,1) = 0.0
            dum1(i,j,nk+1) = 0.0
          enddo
        enddo

        tem = cp*0.5

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
          ! x-dir
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                   -tem*(thv0(i,j,k)+thv0(i-1,j,k))*(         &
                     (ppd(i,j,k)-ppd(i-1,j,k))*rdx*uf(i)    &
              +0.125*( (dum1(i,j,k+1)+dum1(i-1,j,k+1))      &
                      +(dum1(i,j,k  )+dum1(i-1,j,k  )) )    &
                    *(gxu(i,j,k)+gxu(i,j,k+1))    ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)         &
                   -tem*(thv0(i,j,k)+thv0(i,j-1,k))*(         &
                     (ppd(i,j,k)-ppd(i,j-1,k))*rdy*vf(j)    &
              +0.125*( (dum1(i,j,k+1)+dum1(i,j-1,k+1))      &
                      +(dum1(i,j,k  )+dum1(i,j-1,k  )) )    &
                    *(gyv(i,j,k)+gyv(i,j,k+1))    ) )
          enddo
          enddo
        enddo

    ENDIF

  ENDIF  ps27a

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
          ENDIF
        ENDIF

!----------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

!----------------------------------------------

      ENDIF

      if(timestats.ge.1) time_sound=time_sound+mytime()

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  integrate w forward in time:


      IF( wnudge.eq.1 )THEN
        !  updraft nudging tendency:
        IF( (mtime+dt).le.t2_wnudge )THEN
          call get_wnudge(mtime,dts,xh,yh,zf,w3d,dum1)
        ENDIF
      ENDIF

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    ps27b:  &
    IF( psolver.eq.2 )THEN

      IF(.not.terrain_flag)THEN
        ! without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem1,c1a,c1b,c2a,c2b)
        do k=2,nk
          tem1 = rdz*cp*mf(1,1,k)
          c1a = c1(1,1,k)
          c2a = c2(1,1,k)
          do j=1,nj
          do i=1,ni
            w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)                     &
                    -tem1*(ppd(i,j,k)-ppd(i,j,k-1))                     &
                         *(c2a*thv(i,j,k)+c1a*thv(i,j,k-1)) )
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

      ELSE
        ! with terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rds(k)*cp
          do j=1,nj
          do i=1,ni
            w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)                     &
                    -tem1*(ppd(i,j,k)-ppd(i,j,k-1))*gz(i,j)             &
                         *(c2(i,j,k)*thv(i,j,k)+c1(i,j,k)*thv(i,j,k-1)) )
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)

      ENDIF

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    ELSEIF( psolver.eq.7 )THEN  ps27b

      IF(.not.terrain_flag)THEN
        ! without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem1)
        do k=2,nk
          tem1 = rdz*mf(1,1,k)*cp*(c2(1,1,k)*thv0(1,1,k)+c1(1,1,k)*thv0(1,1,k-1))
          do j=1,nj
          do i=1,ni
            w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)            &
                    -tem1*(ppd(i,j,k)-ppd(i,j,k-1)) )
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

      ELSE
        ! with terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem1)
        do k=2,nk
          do j=1,nj
          do i=1,ni
            tem1 = rds(k)*cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))
            w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)                     &
                  -tem1*(ppd(i,j,k)-ppd(i,j,k-1))*gz(i,j) )
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)

      ENDIF

    ENDIF  ps27b

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( wnudge.eq.1 )THEN
        !  apply updraft nudging:
        IF( (mtime+dt).le.t2_wnudge )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            w3d(i,j,k)=w3d(i,j,k)+dts*dum1(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
      ENDIF

      IF( do_ib )THEN
          call zero_out_w(bndy,kbdy,w3d)
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get terms for div,vadv (terrain only):

    doter:  &
    IF(terrain_flag)THEN
      ! Cartesian grid with terrain:

      ps27d:  &
      IF(psolver.eq.2)THEN 

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)*rgzu(i,j)
            dum4(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*dum1(i,j,k)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            dum2(i,j,k)=v3d(i,j,k)*rgzv(i,j)
            dum5(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*dum2(i,j,k)
          enddo
          enddo
        ENDDO
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,r1,r2)
        DO k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do j=1,nj
          do i=1,ni
            dum3(i,j,k)=w3d(i,j,k)                                               &
                       +0.5*( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                               +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                               +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                            )*(sigmaf(k)-zt)*gz(i,j)*rzt
            ! do not couple dum6 with rho0:
            ! (note formulation of pk1,pk2)
            dum6(i,j,k)=w3d(i,j,k)                                               &
                       +0.5*( ( r2*(dum4(i,j,k  )+dum4(i+1,j,k  ))               &
                               +r1*(dum4(i,j,k-1)+dum4(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum5(i,j,k  )+dum5(i,j+1,k  ))               &
                               +r1*(dum5(i,j,k-1)+dum5(i,j+1,k-1)) )*dzdy(i,j)   &
                            )*(sigmaf(k)-zt)*gz(i,j)*rzt * rrf0(i,j,k)
          enddo
          enddo
        ENDDO

      ELSEIF( psolver.eq.7 )THEN  ps27d

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)*rgzu(i,j)*0.25*(thv0(i-1,j,k)+thv0(i,j,k))*(rho0(i-1,j,k)+rho0(i,j,k))
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            dum2(i,j,k)=v3d(i,j,k)*rgzv(i,j)*0.25*(thv0(i,j-1,k)+thv0(i,j,k))*(rho0(i,j-1,k)+rho0(i,j,k))
          enddo
          enddo
        ENDDO

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,r1,r2,rr1,rr2)
        DO k=1,nk
          IF(k.eq.1)THEN
            do j=1,nj
            do i=1,ni
              dum3(i,j,1)=0.0
              dum3(i,j,nk+1)=0.0
            enddo
            enddo
          ELSE
            rr2 = (sigmaf(k)-sigma(k-1))*rds(k)
            rr1 = 1.0-rr2
            r1 = 0.5*rr1
            r2 = 0.5*rr2
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=w3d(i,j,k)*(rr1*thv0(i,j,k-1)+rr2*thv0(i,j,k))*rf0(i,j,k)  &
                         +( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                             +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                           +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                             +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                             )*(sigmaf(k)-zt)*gz(i,j)*rzt
            enddo
            enddo
          ENDIF
        ENDDO

      ENDIF  ps27d

    ENDIF  doter

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get new pp,th

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

  ps27c:  &
  IF( psolver.eq.2 )THEN

      temx = dts*0.5*rdx
      temy = dts*0.5*rdy

    IF( axisymm.eq.0 )THEN

      IF(.not.terrain_flag)THEN
        ! Cartesian grid, without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div,delp)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          div=( (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)    &
               +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j) )  &
               +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
          if(abs(div).lt.smeps) div=0.0
          delp = dts*( ppten(i,j,k)-ppterm(i,j,k)*div )                         &
             +( -( u3d(i  ,j,k)*(piadv(i  ,j,k)-piadv(i-1,j,k))                 &
                  +u3d(i+1,j,k)*(piadv(i+1,j,k)-piadv(i  ,j,k)) )*temx*uh(i)    &
                -( v3d(i,j  ,k)*(piadv(i,j  ,k)-piadv(i,j-1,k))                 &
                  +v3d(i,j+1,k)*(piadv(i,j+1,k)-piadv(i,j  ,k)) )*temy*vh(j) )  &
                +( w3d(i,j,k  )*pk1(i,j,k)  &
                  +w3d(i,j,k+1)*pk2(i,j,k) )
          pp3d(i,j,k)=pp3d(i,j,k)+delp
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*delp
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
        enddo
        enddo
        enddo

      ELSE
        ! Cartesian grid, with terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div,u1,u2,v1,v2,w1,w2,tem)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          div = gz(i,j)*( ( (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i)    &
                           +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j) )  &
                           +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k) )
          if(abs(div).lt.smeps) div=0.0
          u2 = temx*dum4(i+1,j,k)*uh(i)
          u1 = temx*dum4(i  ,j,k)*uh(i)
          v2 = temy*dum5(i,j+1,k)*vh(j)
          v1 = temy*dum5(i,j  ,k)*vh(j)
          w2 = dum6(i,j,k+1)
          w1 = dum6(i,j,k  )
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+dts*( ppten(i,j,k)-ppterm(i,j,k)*div )  &
                 +( -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))                 &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )               &
                    -( v2*(piadv(i,j+1,k)-piadv(i,j  ,k))                 &
                      +v1*(piadv(i,j  ,k)-piadv(i,j-1,k)) ) )*rr0(i,j,k)*gz(i,j) &
                    +( w1*pk1(i,j,k)+w2*pk2(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
        enddo
        enddo
        enddo

      ENDIF

    ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div,u1,u2,v1,v2,w1,w2,tem)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          div=(arh2(i)*u3d(i+1,j,k)-arh1(i)*u3d(i,j,k))*rdx*uh(i)   &
             +(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
          if(abs(div).lt.smeps) div=0.0
          u2 = temx*u3d(i+1,j,k)*uh(i)*arh2(i)
          u1 = temx*u3d(i  ,j,k)*uh(i)*arh1(i)
          w2 = w3d(i,j,k+1)
          w1 = w3d(i,j,k  )
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+dts*( ppten(i,j,k)-ppterm(i,j,k)*div )  &
                    -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))                 &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )               &
                    +( w1*pk1(i,j,k)+w2*pk2(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
        enddo
        enddo
        enddo

    ENDIF

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

  ELSEIF( psolver.eq.7 )THEN  ps27c

    IF( axisymm.eq.0 )THEN

      IF(.not.terrain_flag)THEN
        ! Cartesian grid, without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,r1,r2,div,tem)
        do k=1,nk
          r2 = rf0(1,1,k+1)*(c2(1,1,k+1)*thv0(1,1,k+1)+c1(1,1,k+1)*thv0(1,1,k  ))  &
              *rdz*mh(1,1,k)/(rho0(1,1,k)*thv0(1,1,k))
          r1 = rf0(1,1,k  )*(c2(1,1,k  )*thv0(1,1,k  )+c1(1,1,k  )*thv0(1,1,k-1))  &
              *rdz*mh(1,1,k)/(rho0(1,1,k)*thv0(1,1,k))
          tem = dts*csound*csound/(cp*thv0(1,1,k))
          do j=1,nj
          do i=1,ni
            !-----
            div=( (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)        &
                 +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j) )      &
                 +(r2*w3d(i,j,k+1)-r1*w3d(i,j,k))
            if(abs(div).lt.smeps) div=0.0
            !-----
            ppd(i,j,k)=pp3d(i,j,k)
            pp3d(i,j,k)=pp3d(i,j,k)-tem*div
            if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
            dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
            ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
            !-----
          enddo
        enddo
        enddo

      ELSE
        ! Cartesian grid, with terrain:

        tem = dts*csound*csound/cp

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          !-----
          div=( ( (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i)        &
                 +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j) )      &
                 +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k) )/(rho0(i,j,k)*thv0(i,j,k))
          if(abs(div).lt.smeps) div=0.0
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)-tem*div/thv0(i,j,k)
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
        enddo
        enddo
        enddo

      ENDIF

    ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,r1,r2,div,tem)
        do k=1,nk
          r2 = rf0(1,1,k+1)*(c2(1,1,k+1)*thv0(1,1,k+1)+c1(1,1,k+1)*thv0(1,1,k  ))  &
              *rdz*mh(1,1,k)/(rho0(1,1,k)*thv0(1,1,k))
          r1 = rf0(1,1,k  )*(c2(1,1,k  )*thv0(1,1,k  )+c1(1,1,k  )*thv0(1,1,k-1))  &
              *rdz*mh(1,1,k)/(rho0(1,1,k)*thv0(1,1,k))
          tem = dts*csound*csound/(cp*thv0(1,1,k))
          do j=1,nj
          do i=1,ni
            !-----
            div=(arh2(i)*u3d(i+1,j,k)-arh1(i)*u3d(i,j,k))*rdx*uh(i)        &
               +(r2*w3d(i,j,k+1)-r1*w3d(i,j,k))
            if(abs(div).lt.smeps) div=0.0
            !-----
            ppd(i,j,k)=pp3d(i,j,k)
            pp3d(i,j,k)=pp3d(i,j,k)-tem*div
            if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
            dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
            ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
            !-----
          enddo
        enddo
        enddo

    ENDIF

  ENDIF  ps27c

!c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

        IF( n.lt.nloop )THEN
          if(timestats.ge.1) time_sound=time_sound+mytime()
          call bcp(ppd)
        ENDIF
#ifdef MPI
        IF( n.ne.1 )THEN
          if(timestats.ge.1) time_sound=time_sound+mytime()
          call comm_1p_end2(reqs_p)
        ENDIF
        IF( n.lt.nloop )THEN
          call comm_1p_start(ppd,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        ENDIF
#endif

!--------------------------------------------------------------------
!  time-averaged velocities:

      if( n.lt.nloop )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          rru(i,j,k)=rru(i,j,k)+u3d(i,j,k)
          rrv(i,j,k)=rrv(i,j,k)+v3d(i,j,k)
          rrw(i,j,k)=rrw(i,j,k)+w3d(i,j,k)
        enddo
        enddo
        ENDDO
      else
        tavg = 1.0/float(nloop)
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          rru(i,j,k)=(rru(i,j,k)+u3d(i,j,k))*tavg
          rrv(i,j,k)=(rrv(i,j,k)+v3d(i,j,k))*tavg
          rrw(i,j,k)=(rrw(i,j,k)+w3d(i,j,k))*tavg
        enddo
        enddo
        ENDDO
      endif

!--------------------------------------------------------------------

      ENDDO  small_step_loop

!  end of small steps
!--------------------------------------------------------------------

      IF( nrk.eq.nrkmax )THEN
        ! pressure tendency term: save for next timestep:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppx(i,j,k)=dum1(i,j,k)
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_sound=time_sound+mytime()


      end subroutine sounde

  END MODULE sounde_module
