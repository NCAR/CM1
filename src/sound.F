  MODULE sound_module

  implicit none

  private
  public :: sound

  CONTAINS

      subroutine sound( dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,     &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,zf,      &
                        pi0,rho0,rr0,rf0,rrf0,th0,rth0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,dtu,dtu0,dtv,dtv0,    &
                        dum1,dum2,dum3,dum4,mm  ,bk  ,                    &
                        ppd ,fpk ,qk ,pk1,pk2,wtens,                      &
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
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,rr0,rf0,rrf0,th0,rth0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je) :: dtu,dtv
      real, intent(in),    dimension(ib:ie,jb:je) :: dtu0,dtv0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,mm,bk
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppd,fpk,qk,pk1,pk2,wtens
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
      real :: tem,tem1,tem2,tem3,tem4,r1,r2,dts,rdt
      real :: beta,aa,bb,cc,dd,asq,fwk

      real :: temx,temy,u1,u2,v1,v2,w1,w2,ww,div,tavg,c1a,c1b,c2a,c2b

!---------------------------------------------------------------------
      ! beta = backward weight
      ! alph = forward weight

      beta=1.0-alph
      asq=alph*alph

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
      stop 97392
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
!  Arrays for vertically implicit solver

      IF(.not.terrain_flag)THEN

        ! without terrain:
        ! "s" velocities ARE NOT coupled with reference density

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem,tem1,r1,r2,c1a,c1b,c2a,c2b)
        do k=1,nk
          tem = dts*rdz*mh(1,1,k)
          r2 = dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k+1)
          r1 = dts*rdz*mh(1,1,k)*rr0(1,1,k)*rf0(1,1,k)
          tem1 = dts*rdz*cp*mf(1,1,k)
          c1a = c1(1,1,k  )
          c1b = c1(1,1,k+1)
          c2a = c2(1,1,k  )
          c2b = c2(1,1,k+1)
          do j=1,nj
          do i=1,ni
            qk(i,j,k) = tem*ppterm(i,j,k)
            pk2(i,j,k) = r2*( -c2b*piadv(i,j,k+1)+(1.0-c1b)*piadv(i,j,k) )
            pk1(i,j,k) = r1*( +c1a*piadv(i,j,k-1)+(c2a-1.0)*piadv(i,j,k) )
            mm(i,j,k) = tem1*(c2a*thv(i,j,k)+c1a*thv(i,j,k-1))
          enddo
          enddo
        enddo

      ELSE

        ! with terrain:
        ! "s" velocities ARE coupled with reference density

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tem,tem1,tem2,r1,r2)
        do k=1,nk
          tem = dts*rdsf(k)
          do j=1,nj
          do i=1,ni
            qk(i,j,k) = tem*ppterm(i,j,k)*gz(i,j)
            tem2 = dts*gz(i,j)*rdsf(k)*rr0(i,j,k)
            pk2(i,j,k) = tem2*rf0(i,j,k+1)*( -c2(i,j,k+1)*piadv(i,j,k+1)+(1.0-c1(i,j,k+1))*piadv(i,j,k) )
            pk1(i,j,k) = tem2*rf0(i,j,k  )*( +c1(i,j,k  )*piadv(i,j,k-1)+(c2(i,j,k  )-1.0)*piadv(i,j,k) )
          enddo
          enddo
          IF( k.ge.2 )THEN
            tem1 = dts*rds(k)*cp
            do j=1,nj
            do i=1,ni
              mm(i,j,k) = tem1*(c2(i,j,k)*thv(i,j,k)+c1(i,j,k)*thv(i,j,k-1))*gz(i,j)
            enddo
            enddo
          ELSE
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
          rrw(i,j,k)=wa(i,j,k)*beta
        enddo
        enddo
        enddo

      else

        if( axisymm.eq.0 )then

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
            rrw(i,j,k)=wa(i,j,k)*beta
            pp3d(i,j,k)=ppi(i,j,k)
            ppd(i,j,k)=ppx(i,j,k)
          enddo
          enddo
          enddo

        else

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            u3d(i,j,k)=ua(i,j,k)
            rru(i,j,k)=0.0
            w3d(i,j,k)=wa(i,j,k)
            rrw(i,j,k)=wa(i,j,k)*beta
            pp3d(i,j,k)=ppi(i,j,k)
            ppd(i,j,k)=ppx(i,j,k)
          enddo
          enddo
          enddo

        endif

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

!-----

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
            u3d(i,j,k)=u3d(i,j,k)+dtu(i,j)*( uten(i,j,k)    &
                   -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)    &
                        *(thv(i,j,k)+thv(i-1,j,k)) )
            v3d(i,j,k)=v3d(i,j,k)+dtv(i,j)*( vten(i,j,k)    &
                   -tem2*(ppd(i,j,k)-ppd(i,j-1,k))*vf(j)    &
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
          do i=2,ni+1
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

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get terms for horiz div (terrain only):

      IF(terrain_flag)THEN
        ! Cartesian grid with terrain:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=u3d(i,j,k)*rgzu(i,j)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            dum2(i,j,k)=v3d(i,j,k)*rgzv(i,j)
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
            dum3(i,j,k)=0.5*( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                               +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                               +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                            )*(sigmaf(k)-zt)*gz(i,j)*rzt
          enddo
          enddo
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  update arrays for vertically implicit solver

      temx = dts*0.5*rdx
      temy = dts*0.5*rdy

    IF(axisymm.eq.0)THEN

      IF(.not.terrain_flag)THEN
        ! Cartesian grid without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div,u1,u2,v1,v2,w1,w2)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          div = (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)  &
               +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)
          if(abs(div).lt.smeps) div=0.0
          fpk(i,j,k)=dts*( ppten(i,j,k)-ppterm(i,j,k)*div )                        &
                 +( -( u3d(i  ,j,k)*(piadv(i  ,j,k)-piadv(i-1,j,k))                &
                      +u3d(i+1,j,k)*(piadv(i+1,j,k)-piadv(i  ,j,k)) )*temx*uh(i)   &
                    -( v3d(i,j  ,k)*(piadv(i,j  ,k)-piadv(i,j-1,k))                &
                      +v3d(i,j+1,k)*(piadv(i,j+1,k)-piadv(i,j  ,k)) )*temy*vh(j) ) &
              +beta*( w3d(i,j,k  )*(pk1(i,j,k)+qk(i,j,k))                          &
                     +w3d(i,j,k+1)*(pk2(i,j,k)-qk(i,j,k)) )
        enddo
        enddo
        ENDDO

      ELSE
        ! Cartesian grid with terrain:
        ! (note: see below for advection)

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div,u1,u2,v1,v2,w1,w2)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          div = gz(i,j)*( ( (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i)    &
                           +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j) )  &
                           +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k) )
          if(abs(div).lt.smeps) div=0.0
          w1 = w3d(i,j,k)
          w2 = w3d(i,j,k+1)
          if(k.eq. 1) w1 = 0.0
          if(k.eq.nk) w2 = 0.0
          fpk(i,j,k)=dts*( ppten(i,j,k)-ppterm(i,j,k)*div )  &
              +beta*( -qk(i,j,k)*(w2-w1) )
        enddo
        enddo
        ENDDO

      ENDIF

    ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,div,u1,u2,v1,v2,w1,w2)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          div = (arh2(i)*u3d(i+1,j,k)-arh1(i)*u3d(i,j,k))*rdx*uh(i)
          if(abs(div).lt.smeps) div=0.0
          u2 = temx*u3d(i+1,j,k)*uh(i)*arh2(i)
          u1 = temx*u3d(i  ,j,k)*uh(i)*arh1(i)
          fpk(i,j,k)=dts*( ppten(i,j,k)-ppterm(i,j,k)*div )     &
                    -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))       &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )     &
              +beta*( w3d(i,j,k  )*(pk1(i,j,k)+qk(i,j,k))       &
                     +w3d(i,j,k+1)*(pk2(i,j,k)-qk(i,j,k)) )
        enddo
        enddo
        ENDDO

    ENDIF

      if(timestats.ge.1) time_sound=time_sound+mytime()


!----------------------------------------------
!  updraft nudging:


      IF( wnudge.eq.1 )THEN
        !  updraft nudging tendency:
        !  here in sound.F, add to fwk array:
        IF( (mtime+dt).le.t2_wnudge )THEN
          if(n.eq.1)then
            do k=2,nk
            do j=1,nj
            do i=1,ni
              wtens(i,j,k) = wten(i,j,k)
            enddo
            enddo
            enddo
          endif
          call get_wnudge(mtime,dts,xh,yh,zf,w3d,dum4)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = wtens(i,j,k) + dum4(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
      ENDIF


!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  advection for terrain:

      IF(terrain_flag)THEN
        ! Cartesian grid with terrain:
        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
        ! Cartesian grid with terrain:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni+1
            dum1(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*dum1(i,j,k)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            dum2(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*dum2(i,j,k)
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
            ! NOTE:  dum3 is NOT coupled with density
            dum3(i,j,k)=0.5*( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                               +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                               +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                            )*(sigmaf(k)-zt)*gz(i,j)*rzt * rrf0(i,j,k)
          enddo
          enddo
        ENDDO
        temx = dts*0.5*rdx
        temy = dts*0.5*rdy
!$omp parallel do default(shared)   &
!$omp private(i,j,k,u1,u2,v1,v2,w1,w2)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          u2 = temx*dum1(i+1,j,k)*uh(i)
          u1 = temx*dum1(i  ,j,k)*uh(i)
          v2 = temy*dum2(i,j+1,k)*vh(j)
          v1 = temy*dum2(i,j  ,k)*vh(j)
          w2 = dum3(i,j,k+1)+beta*w3d(i,j,k+1)
          w1 = dum3(i,j,k  )+beta*w3d(i,j,k  )
          fpk(i,j,k)=fpk(i,j,k)                                           &
                 +( -( u2*(piadv(i+1,j,k)-piadv(i  ,j,k))                 &
                      +u1*(piadv(i  ,j,k)-piadv(i-1,j,k)) )               &
                    -( v2*(piadv(i,j+1,k)-piadv(i,j  ,k))                 &
                      +v1*(piadv(i,j  ,k)-piadv(i,j-1,k)) ) )*rr0(i,j,k)*gz(i,j) &
                    +( w1*pk1(i,j,k)+w2*pk2(i,j,k) )
        enddo
        enddo
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  the vertcally implicit solver .... get new w

        k=2

        !$omp parallel do default(shared)   &
        !$omp private(i,j,aa,bb,cc,dd,fwk,r1)
        do j=1,nj
        do i=1,ni
          cc =      -asq*( -mm(i,j,k)*( pk2(i,j,k  )-qk(i,j,k  ) ) )
          bb = 1.0-asq*( +mm(i,j,k)*(-(pk1(i,j,k  )+qk(i,j,k  ))   &
                                                +(pk2(i,j,k-1)-qk(i,j,k-1)) ) )
          fwk=dts*wten(i,j,k)+beta*( -mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1)) )
          dd = w3d(i,j,k)+fwk       +alph*(                             &
                 -mm(i,j,k)*( (pp3d(i,j,k  )+fpk(i,j,k  ))              &
                             -(pp3d(i,j,k-1)+fpk(i,j,k-1)) ) )
          r1 = 1.0/bb
          dum1(i,j,k) = -cc*r1
          dum2(i,j,k) =  dd*r1
        enddo
        enddo

        do k=3,(nk-1)

        !$omp parallel do default(shared)   &
        !$omp private(i,j,aa,bb,cc,dd,fwk,r1)
        do j=1,nj
        do i=1,ni
          aa =      -asq*( +mm(i,j,k)*( pk1(i,j,k-1)+qk(i,j,k-1) ) )
          cc =      -asq*( -mm(i,j,k)*( pk2(i,j,k  )-qk(i,j,k  ) ) )
          bb = 1.0-asq*( +mm(i,j,k)*(-(pk1(i,j,k  )+qk(i,j,k  ))   &
                                                +(pk2(i,j,k-1)-qk(i,j,k-1)) ) )
          fwk=dts*wten(i,j,k)+beta*( -mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1)) )
          dd = w3d(i,j,k)+fwk       +alph*(                             &
                 -mm(i,j,k)*( (pp3d(i,j,k  )+fpk(i,j,k  ))              &
                             -(pp3d(i,j,k-1)+fpk(i,j,k-1)) ) )
          r1 = 1.0/(aa*dum1(i,j,k-1)+bb)
          dum1(i,j,k) = -cc*r1
          dum2(i,j,k) = (dd-aa*dum2(i,j,k-1))*r1
        enddo
        enddo

        enddo

        k = nk

        !$omp parallel do default(shared)   &
        !$omp private(i,j,aa,bb,cc,dd,fwk,r1)
        do j=1,nj
        do i=1,ni
          aa =      -asq*( +mm(i,j,k)*( pk1(i,j,k-1)+qk(i,j,k-1) ) )
          bb = 1.0-asq*( +mm(i,j,k)*(-(pk1(i,j,k  )+qk(i,j,k  ))   &
                                                +(pk2(i,j,k-1)-qk(i,j,k-1)) ) )
          fwk=dts*wten(i,j,k)+beta*( -mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1)) )
          dd = w3d(i,j,k)+fwk       +alph*(                             &
                 -mm(i,j,k)*( (pp3d(i,j,k  )+fpk(i,j,k  ))              &
                             -(pp3d(i,j,k-1)+fpk(i,j,k-1)) ) )
          r1 = 1.0/(aa*dum1(i,j,k-1)+bb)
          w3d(i,j,k) = (dd-aa*dum2(i,j,k-1))*r1
        enddo
        enddo

        do k=(nk-1),2,-1

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=dum1(i,j,k)*w3d(i,j,k+1)+dum2(i,j,k)
        enddo
        enddo

        enddo

      IF( do_ib )THEN
          call zero_out_w(bndy,kbdy,w3d)
      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get new pp,th


      IF(.not.terrain_flag)THEN
        ! without terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,w1,w2)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+fpk(i,j,k)              &
              +alph*( w3d(i,j,k  )*(pk1(i,j,k)+qk(i,j,k)) &
                     +w3d(i,j,k+1)*(pk2(i,j,k)-qk(i,j,k)) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
        enddo
        enddo
        ENDDO

      ELSE
        ! with terrain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,w1,w2)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          !-----
          w1 = w3d(i,j,k)
          w2 = w3d(i,j,k+1)
          if(k.eq. 1) w1 = 0.0
          if(k.eq.nk) w2 = 0.0
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+fpk(i,j,k)    &
              +alph*( w1*(pk1(i,j,k)+qk(i,j,k)) &
                     +w2*(pk2(i,j,k)-qk(i,j,k)) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          dum1(i,j,k)=kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          ppd(i,j,k)=pp3d(i,j,k)+dum1(i,j,k)
          !-----
        enddo
        enddo
        ENDDO

      ENDIF

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
          rrw(i,j,k)=(rrw(i,j,k)+w3d(i,j,k)*alph)*tavg
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


      end subroutine sound

  END MODULE sound_module
