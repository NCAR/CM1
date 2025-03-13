  MODULE anelp_module

  implicit none

  private
  public :: anelp

  CONTAINS

      subroutine anelp(xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,             &
                       zh,mh,rmh,mf,rmf,zf,pi0,thv0,rho0,prs0,rf0,  &
                       rds,sigma,rdsf,sigmaf,                       &
                       gz,rgz,gzu,rgzu,gzv,rgzv,                    &
                       dzdx,dzdy,gx,gxu,gy,gyv,                     &
                       radbcw,radbce,radbcs,radbcn,                 &
                       dum1,dum2,dum3,dum4,divx,                    &
                       u0,ua,u3d,uten,                              &
                       v0,va,v3d,vten,                              &
                       wa,w3d,wten,                                 &
                       ppi,pp3d,phi1,phi2,cfb,cfa,cfc,              &
                        d1, d2,pdt,lgbth,lgbph,rhs,trans,dttmp,nrk,rtime,mtime)
      use input
      use constants
      use misclibs , only : convinitu,convinitv,get_wnudge
      use bc_module
      use poiss_module
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf,zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rho0,prs0,rf0
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(in), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(in), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in   ), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d,uten
      real, intent(in   ), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d,vten
      real, intent(in   ), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(in), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, intent(in), dimension(kpb:kpe) :: cfa,cfc,d1,d2
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,lgbth,lgbph
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      real, intent(in) :: dttmp
      integer, intent(in) :: nrk
      real, intent(in) :: rtime
      double precision, intent(in) :: mtime

!-----

      integer :: i,j,k
      real :: tem,temx,temy,temz,rdt
      double precision :: fluxout,fluxin,u1,v1

!---------------------------------------------------------------------

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)
 
        endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Open boundary conditions:

        IF( (wbc.eq.2.and.ibw.eq.1) .or. (ebc.eq.2.and.ibe.eq.1) )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(1,j,k)=ua(1,j,k)
            u3d(ni+1,j,k)=ua(ni+1,j,k)
          enddo
          enddo
        ENDIF

        IF(wbc.eq.2.and.ibw.eq.1)THEN
          ! west open bc tendency:
          call   ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dttmp)
        ENDIF
        IF(ebc.eq.2.and.ibe.eq.1)THEN
          ! east open bc tendency:
          call   ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dttmp)
        ENDIF


        IF( (wbc.eq.2.and.ibw.eq.1) .or. (ebc.eq.2.and.ibe.eq.1) )THEN

          call restrict_openbc_we(rvh,rmh,rho0,u3d)

          rdt = 1.0/dttmp
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj
            uten(1,j,k)=(u3d(1,j,k)-ua(1,j,k))*rdt
            uten(ni+1,j,k)=(u3d(ni+1,j,k)-ua(ni+1,j,k))*rdt
          enddo
          enddo

        ENDIF

!-----

      IF(axisymm.eq.0)THEN

        IF( (sbc.eq.2.and.ibs.eq.1) .or. (nbc.eq.2.and.ibn.eq.1) )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,1,k)=va(i,1,k)
            v3d(i,nj+1,k)=va(i,nj+1,k)
          enddo
          enddo
        ENDIF

        IF(sbc.eq.2.and.ibs.eq.1)THEN
          ! south open bc tendency:
          call   ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dttmp)
        ENDIF
        IF(nbc.eq.2.and.ibn.eq.1)THEN
          ! north open bc tendency:
          call   ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dttmp)
        ENDIF


        IF( (sbc.eq.2.and.ibs.eq.1) .or. (nbc.eq.2.and.ibn.eq.1) )THEN

          call restrict_openbc_sn(ruh,rmh,rho0,v3d)

          rdt = 1.0/dttmp
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni
            vten(i,1,k)=(v3d(i,1,k)-va(i,1,k))*rdt
            vten(i,nj+1,k)=(v3d(i,nj+1,k)-va(i,nj+1,k))*rdt
          enddo
          enddo

        ENDIF

      ENDIF

!---------------------------------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni+1
              u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)
            enddo
            enddo
            enddo
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
            tem=1.0/dttmp
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni+1
              uten(i,j,k) = (u3d(i,j,k)-ua(i,j,k))*tem
            enddo
            enddo
            enddo
          ENDIF
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj+1
            do i=1,ni
              v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)
            enddo
            enddo
            enddo
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
            tem=1.0/dttmp
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj+1
            do i=1,ni
              vten(i,j,k) = (v3d(i,j,k)-va(i,j,k))*tem
            enddo
            enddo
            enddo
          ENDIF
        ENDIF

!----------------------------------------------
!  updraft nudging:

      IF( wnudge.eq.1 )THEN
        IF( (mtime+dttmp).le.t2_wnudge )THEN
          call get_wnudge(mtime,dttmp,xh,yh,zf,w3d,dum1)
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = wten(i,j,k)+dum1(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
      ENDIF

        if(timestats.ge.1) time_sound=time_sound+mytime()

!---------------------------------------------------------------------
!  Get pressure

        call     poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                       dum1,dum2,dum3,dum4,divx,phi2,uten,vten,wten,  &
                       cfb,cfa,cfc,d1,d2,pdt,lgbth,lgbph,rhs,trans,dttmp)

        call bcs(phi2)

!---------------------------------------------------------------------

        temx=dttmp*rdx
        temy=dttmp*rdy
        temz=dttmp*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)             &
                  -temx*(phi2(i,j,k)-phi2(i-1,j,k))*uf(i)
          v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)             &
                  -temy*(phi2(i,j,k)-phi2(i,j-1,k))*vf(j)
          w3d(i,j,k)=wa(i,j,k)+dttmp*wten(i,j,k)             &
                  -temz*(phi2(i,j,k)-phi2(i,j,k-1))*mf(i,j,k)
        enddo
        enddo
        enddo

!-------------------------------------------------------------------- 

      IF( nrk.eq.nrkmax )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          phi1(i,j,k) = phi2(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-------------------------------------------------------------------- 

      end subroutine anelp

  END MODULE anelp_module
