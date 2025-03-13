  MODULE soundcb_module

  implicit none

  private
  public :: soundcb

  CONTAINS

      subroutine soundcb(dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,    &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,zf,      &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,                      &
                        dum1,dum2,dum3,dum4,dum5,dum6,                    &
                        ppd ,fpk ,qk ,pk1,pk2,ftk,                        &
                        u0,rru,ua,u3d,uten,                               &
                        v0,rrv,va,v3d,vten,                               &
                        rrw,wa,w3d,wten,                                  &
                        ppi,pp3d,piadv,ppten,ppx,phi1,phi2,               &
                        thv,ppterm,nrk,dttmp,rtime,mtime,get_time_avg,    &
                        bndy,kbdy,                                        &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      use input
      use constants
      use misclibs , only : convinitu,convinitv,get_wnudge,getdiv
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
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,rr0,rf0,rrf0,th0,rth0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppd,fpk,qk,pk1,pk2,ftk
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,ua,u3d,uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,va,v3d,vten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,piadv,ppten,ppx
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
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
      real :: tem,tem1,tem2,tem3,tem4,r1,r2,dts

      real :: temx,temy,tavg

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
      stop 97394
    ENDIF

!!!      print *,'  nloop,dts,dttmp = ',nloop,dts,nloop*dts

!---------------------------------------------------------------------
!  Arrays for vadv:

      IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem,tem1,r1,r2)
            do j=1,nj
            do i=1,ni
              dum3(i,j,1)=0.0
              dum3(i,j,nk+1)=0.0
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
        enddo
        enddo
        enddo

      else

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1,ni+1
            u3d(i,j,k)=ua(i,j,k)
          enddo
          enddo
          IF(axisymm.eq.0)THEN
            ! Cartesian grid:
            do j=1,nj+1
            do i=1,ni
              v3d(i,j,k)=va(i,j,k)
            enddo
            enddo
          ENDIF
          IF(k.ge.2)THEN
            do j=1,nj
            do i=1,ni
              w3d(i,j,k)=wa(i,j,k)
            enddo
            enddo
          ENDIF
          do j=1,nj
          do i=1,ni
            phi2(i,j,k)=phi1(i,j,k)
          enddo
          enddo
          do j=0,nj+1
          do i=0,ni+1
            ppd(i,j,k)=ppx(i,j,k)
          enddo
          enddo
        enddo

      endif

!---------------------------------------------------------------------
!  time-averaged velocities:

      IF( get_time_avg )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            rru(i,j,k)=0.0
          enddo
          enddo
          IF(axisymm.eq.0)THEN
            do j=1,nj+1
            do i=1,ni
              rrv(i,j,k)=0.0
            enddo
            enddo
          ENDIF
          IF(k.ge.2)THEN
            do j=1,nj
            do i=1,ni
              rrw(i,j,k)=0.0
            enddo
            enddo
          ENDIF
        ENDDO
      ENDIF
      if(timestats.ge.1) time_sound=time_sound+mytime()

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Begin small steps:

      small_step_loop:  DO N=1,NLOOP

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

!-----

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)      &
                    -rdx*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)      &
                    -rdy*(ppd(i,j,k)-ppd(i,j-1,k))*vf(j) )
          enddo
          enddo
        enddo

        IF( do_ib )THEN
          call zero_out_uv(bndy,kbdy,u3d,v3d)
        ENDIF

      ELSE
        ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)      &
                    -rdx*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i) )
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

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          ! x-dir
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                  -( (ppd(i,j,k)-ppd(i-1,j,k))*rdx*uf(i)    &
              +0.125*( (dum1(i,j,k+1)+dum1(i-1,j,k+1))      &
                      +(dum1(i,j,k  )+dum1(i-1,j,k  )) )    &
                    *(gxu(i,j,k)+gxu(i,j,k+1))    ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)         &
                  -( (ppd(i,j,k)-ppd(i,j-1,k))*rdy*vf(j)    &
              +0.125*( (dum1(i,j,k+1)+dum1(i,j-1,k+1))      &
                      +(dum1(i,j,k  )+dum1(i,j-1,k  )) )    &
                    *(gyv(i,j,k)+gyv(i,j,k+1))    ) )
          enddo
          enddo
        enddo

    ENDIF

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
!  convergence forcing:

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:
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

      IF(.not.terrain_flag)THEN
        ! without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rdz*mf(1,1,k)
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)            &
                  -tem1*(ppd(i,j,k)-ppd(i,j,k-1)) )
        enddo
        enddo
        enddo

      ELSE
        ! with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=w3d(i,j,k)+dts*( wten(i,j,k)                     &
                -rds(k)*(ppd(i,j,k)-ppd(i,j,k-1))*gz(i,j) )
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)

      ENDIF

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
!  get divergence:


        call     getdiv(arh1,arh2,uh,vh,mh,u3d,v3d,w3d,dum1,dum2,dum3,dum4,  &
                        rds,rdsf,sigma,sigmaf,gz,rgzu,rgzv,dzdx,dzdy)


!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get new pp,th

      tem = dts*csound*csound

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
        do j=1,nj
        do i=1,ni
          ppd(i,j,k)=phi2(i,j,k)
          phi2(i,j,k)=phi2(i,j,k)-tem*dum4(i,j,k)
          if(abs(phi2(i,j,k)).lt.smeps) phi2(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( phi2(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=phi2(i,j,k)+dum1(i,j,k)
        enddo
        enddo
      ENDDO
      if(timestats.ge.1) time_sound=time_sound+mytime()

        IF( n.lt.nloop )THEN
          call bcp(ppd)
        ENDIF
#ifdef MPI
        IF( n.ne.1 )THEN
          call comm_1p_end2(reqs_p)
        ENDIF
        IF( n.lt.nloop )THEN
          call comm_1p_start(ppd,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        ENDIF
#endif

!--------------------------------------------------------------------
!  time-averaged velocities:

      IF( get_time_avg )THEN

      if( n.lt.nloop )then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            rru(i,j,k)=rru(i,j,k)+u3d(i,j,k)
          enddo
          enddo
          IF( axisymm.eq.0 )THEN
            do j=1,nj+1
            do i=1,ni
              rrv(i,j,k)=rrv(i,j,k)+v3d(i,j,k)
            enddo
            enddo
          ENDIF
          IF( k.ge.2 )THEN
            do j=1,nj
            do i=1,ni
              rrw(i,j,k)=rrw(i,j,k)+w3d(i,j,k)
            enddo
            enddo
          ENDIF
        ENDDO
      else
        tavg = 1.0/float(nloop)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            rru(i,j,k)=(rru(i,j,k)+u3d(i,j,k))*tavg
          enddo
          enddo
          IF( axisymm.eq.0 )THEN
            do j=1,nj+1
            do i=1,ni
              rrv(i,j,k)=(rrv(i,j,k)+v3d(i,j,k))*tavg
            enddo
            enddo
          ENDIF
          IF( k.ge.2 )THEN
            do j=1,nj
            do i=1,ni
              rrw(i,j,k)=(rrw(i,j,k)+w3d(i,j,k))*tavg
            enddo
            enddo
          ENDIF
        ENDDO
      endif

      ENDIF

!--------------------------------------------------------------------

      ENDDO  small_step_loop

!  end of small steps
!--------------------------------------------------------------------


      if(timestats.ge.1) time_sound=time_sound+mytime()
      call bcs(phi2)
#ifdef MPI
      call comm_1s_start(phi2,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      call comm_1s_end(  phi2,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
#endif


      IF( nrk.eq.nrkmax )THEN

        ! pressure tendency term: save for next timestep:
        ! 170823:  also, step forward phi

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1,ni
            ppx(i,j,k)=dum1(i,j,k)
          enddo
          enddo
          do j=0,nj+1
          do i=0,ni+1
            phi1(i,j,k)=phi2(i,j,k)
          enddo
          enddo
        enddo

      ENDIF
      if(timestats.ge.1) time_sound=time_sound+mytime()


      end subroutine soundcb

  END MODULE soundcb_module
