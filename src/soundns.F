  MODULE soundns_module

  implicit none

  private
  public :: soundns

  CONTAINS

      subroutine soundns(xh,rxh,arh1,arh2,uh,xf,uf,yh,vh,yf,vf,           &
                         zh,mh,c1,c2,mf,zf,pi0,thv0,rr0,rf0,              &
                         rds,sigma,rdsf,sigmaf,                           &
                         zs,gz,rgz,gzu,rgzu,gzv,rgzv,                     &
                         dzdx,dzdy,gx,gxu,gy,gyv,                         &
                         radbcw,radbce,radbcs,radbcn,dtu,dtu0,dtv,dtv0,   &
                         dum1,dum2,dum3,div ,                             &
                         u0,ua,u3d,uten,                                  &
                         v0,va,v3d,vten,                                  &
                         wa,w3d,wten,                                     &
                         ppi,pp3d,ppten,thv,ppterm,dttmp,nrk,rtime,mtime, &
                         bndy,kbdy)
      use input
      use constants
      use misclibs , only : convinitu,convinitv,get_wnudge,getdiv
      use bc_module
      use ib_module
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rr0,rf0
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(in), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(in), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je) :: dtu,dtv
      real, intent(in),    dimension(ib:ie,jb:je) :: dtu0,dtv0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,div
      real, intent(in   ), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d,uten
      real, intent(in   ), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d,vten
      real, intent(in   ), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d,wten
      real, intent(in   ), dimension(ib:ie,jb:je,kb:ke) :: ppi,ppten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: pp3d
      real, intent(in   ), dimension(ib:ie,jb:je,kb:ke) :: thv,ppterm
      real, intent(in) :: dttmp
      integer, intent(in) :: nrk
      real, intent(in) :: rtime
      double precision, intent(in) :: mtime
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy

!-----

      integer :: i,j,k
      real :: tem,tem1,tem2,r1,r2

      integer :: pgform

      do j=1,nj+1
      do i=1,ni+1
        dtu(i,j) = dttmp*dtu0(i,j)
        dtv(i,j) = dttmp*dtv0(i,j)
      enddo
      enddo

!---------------------------------------------------------------------

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)
 
        endif

!-----

        if(wbc.eq.2.and.ibw.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(1,j,k)=ua(1,j,k)
          enddo
          enddo
          call   ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dttmp)
        endif

        if(ebc.eq.2.and.ibe.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(ni+1,j,k)=ua(ni+1,j,k)
          enddo
          enddo
          call   ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dttmp)
        endif

!-----

      IF(axisymm.eq.0)THEN

        if(sbc.eq.2.and.ibs.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,1,k)=va(i,1,k)
          enddo
          enddo
          call   ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dttmp)
        endif
 
        if(nbc.eq.2.and.ibn.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,nj+1,k)=va(i,nj+1,k)
          enddo
          enddo
          call   ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dttmp)
        endif

      ENDIF

!-----

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN

        tem1=dttmp*rdx*cp*0.5
        tem2=dttmp*rdy*cp*0.5

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            u3d(i,j,k)=ua(i,j,k)+dtu(i,j)*uten(i,j,k)          &
                   -(tem1*(pp3d(i,j,k)-pp3d(i-1,j,k))*uf(i)    &
                         *(thv(i,j,k)+thv(i-1,j,k)))
            v3d(i,j,k)=va(i,j,k)+dtv(i,j)*vten(i,j,k)          &
                   -(tem2*(pp3d(i,j,k)-pp3d(i,j-1,k))*vf(j)    &
                         *(thv(i,j,k)+thv(i,j-1,k)))
          enddo
          enddo
        enddo

        IF( do_ib )THEN
          call zero_out_uv(bndy,kbdy,u3d,v3d)
        ENDIF

      ELSE

        tem1=dttmp*rdx*cp*0.5

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1+ibw,ni+1-ibe
          u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)             &
                 -(tem1*(pp3d(i,j,k)-pp3d(i-1,j,k))*uf(i)    &
                       *(thv(i,j,k)+thv(i-1,j,k)))
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
            dum1(i,j,k) = (pp3d(i,j,k)-pp3d(i,j,k-1))*rds(k)
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
            u3d(i,j,k)=ua(i,j,k)+dttmp*( uten(i,j,k)        &
                   -tem*(thv(i,j,k)+thv(i-1,j,k))*(         &
                     (pp3d(i,j,k)-pp3d(i-1,j,k))*rdx*uf(i)  &
              +0.125*( (dum1(i,j,k+1)+dum1(i-1,j,k+1))      &
                      +(dum1(i,j,k  )+dum1(i-1,j,k  )) )    &
                    *(gxu(i,j,k)+gxu(i,j,k+1))    ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=va(i,j,k)+dttmp*( vten(i,j,k)        &
                   -tem*(thv(i,j,k)+thv(i,j-1,k))*(         &
                     (pp3d(i,j,k)-pp3d(i,j-1,k))*rdy*vf(j)  &
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

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

      ENDIF

!-----

      IF( wnudge.eq.1 )THEN
        !  updraft nudging tendency:
        IF( (mtime+dttmp).le.t2_wnudge )THEN
          call get_wnudge(mtime,dttmp,xh,yh,zf,w3d,dum1)
        ENDIF
      ENDIF

    IF(.not.terrain_flag)THEN
        ! without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rdz*cp*mf(1,1,k)
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=wa(i,j,k)+dttmp*( wten(i,j,k)                    &
                  -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))                   &
                       *(c2(1,1,k)*thv(i,j,k)+c1(1,1,k)*thv(i,j,k-1)) )
        enddo
        enddo
        enddo

      ELSE
        ! with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
        do k=2,nk
        tem1 = rds(k)*cp
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=wa(i,j,k)+dttmp*( wten(i,j,k)                    &
                  -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))*gz(i,j)           &
                       *(c2(i,j,k)*thv(i,j,k)+c1(i,j,k)*thv(i,j,k-1)) )
        enddo
        enddo
        enddo

      ENDIF

      IF( wnudge.eq.1 )THEN
        !  apply updraft nudging:
        IF( (mtime+dttmp).le.t2_wnudge )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            w3d(i,j,k)=w3d(i,j,k)+dttmp*dum1(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
      ENDIF

      IF( do_ib )THEN
          call zero_out_w(bndy,kbdy,w3d)
      ENDIF

        if(timestats.ge.1) time_sound=time_sound+mytime()

!----------------------------------------------

        call     getdiv(arh1,arh2,uh,vh,mh,u3d,v3d,w3d,dum1,dum2,dum3,div,   &
                        rds,rdsf,sigma,sigmaf,gz,rgzu,rgzv,dzdx,dzdy)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pp3d(i,j,k)=ppi(i,j,k)+dttmp*( ppten(i,j,k)-ppterm(i,j,k)*div(i,j,k) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
        enddo
        enddo
        enddo

        if(timestats.ge.1) time_sound=time_sound+mytime()
 
!-------------------------------------------------------------------- 

      end subroutine soundns

  END MODULE soundns_module
