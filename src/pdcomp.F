  MODULE pdcomp_module

  implicit none

  private
  public :: pidcomp

  CONTAINS

      subroutine pidcomp(dt,xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,          &
                         gz,rgz,gzu,gzv,mh,rmh,mf,rmf,rds,rdsf,c1,c2,f2d,wprof,     &
                         pi0,th0,rth0,thv0,qv0,qc0,qi0,rho0,rr0,rf0,rrf0,u0,v0,     &
                         dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,divx,              &
                         u3d,rru,uten,uten1,v3d,rrv,vten,vten1,w3d,rrw,wten,wten1,  &
                         rho,pp3d,th3d,q3d,udiag,vdiag,wdiag,pdiag,                 &
                         cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans,             &
                         bndy,kbdy,hflxw,hflxe,hflxs,hflxn)

      use input
      use constants
      use adv_module
      use poiss_module
      use ib_module

      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(kb:ke) :: rds
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je) :: f2d
      real, intent(in), dimension(kb:ke) :: wprof
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,rth0,thv0,qv0,qc0,qi0,rho0,rr0,rf0,rrf0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,divx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,uten,uten1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,vten,vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wten,wten1
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,pp3d,th3d
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
      real, intent(inout) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag
      real, intent(in), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, intent(in), dimension(kpb:kpe) :: cfa,cfc,ad1,ad2
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,lgbth,lgbph
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn

      integer :: i,j,k,n
      real :: tem,dudz,dvdz
      double precision :: pavg1,pavg2

      !--------------------------------------------------------------!
      ! buoyancy pressure:

        IF( imoist.eq.1 )THEN
          ! buoyancy (with moisture terms):

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,n)
          do j=1,nj
            !---
            do k=1,nk
            do i=1,ni
              dum1(i,j,k) = 0.0
            enddo
            enddo
            !---
            do n=nql1,nql2
            do k=1,nk
            do i=1,ni
              dum1(i,j,k) = dum1(i,j,k)+q3d(i,j,k,n)
            enddo
            enddo
            enddo
            !---
            if( iice.eq.1 )then
            do n=nqs1,nqs2
            do k=1,nk
            do i=1,ni
              dum1(i,j,k) = dum1(i,j,k)+q3d(i,j,k,n)
            enddo
            enddo
            enddo
            endif
            !---
            do k=1,nk
            do i=1,ni
              dum6(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k) + repsm1*(q3d(i,j,k,nqv)-qv0(i,j,k)) - (dum1(i,j,k)-qc0(i,j,k)-qi0(i,j,k)) )
            enddo
            enddo
            !---
            do k=2,nk
            do i=1,ni
              wten(i,j,k)=(c1(i,j,k)*dum6(i,j,k-1)+c2(i,j,k)*dum6(i,j,k))
            enddo
            enddo
            !---
            do i=1,ni
              wten(i,j,1) = 0.0
              wten(i,j,nk+1) = 0.0
            enddo
            !---
          enddo

        ELSE
          ! buoyancy (dry):

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
            !---
            do k=1,nk
            do i=1,ni
              dum6(i,j,k) = g*th3d(i,j,k)*rth0(i,j,k)
            enddo
            enddo
            !---
            do k=2,nk
            do i=1,ni
              wten(i,j,k)=(c1(i,j,k)*dum6(i,j,k-1)+c2(i,j,k)*dum6(i,j,k))
            enddo
            enddo
            !---
            do i=1,ni
              wten(i,j,1) = 0.0
              wten(i,j,nk+1) = 0.0
            enddo
            !---
          enddo

        ENDIF

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          uten(i,j,k) = 0.0
          vten(i,j,k) = 0.0
          divx(i,j,k) = 0.0
        enddo
        enddo
        enddo

        call   poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                     dum2,dum3,dum4,dum5,divx,pdiag(ibdp,jbdp,kbdp,1),uten,vten,wten,  &
                     cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans,dt)

        ! Adjust mean pressure ... set average value at upper-most level to 
        ! mean of actual pp3d field:

        pavg1 = 0.0
        pavg2 = 0.0

        do j=1,nj
        do i=1,ni
          pavg1 = pavg1+pp3d(i,j,nk)
          pavg2 = pavg2+pdiag(i,j,nk,1)
        enddo
        enddo

        pavg1 = pavg1/dble(nx*ny)
        pavg2 = pavg2/dble(nx*ny)

        tem = sngl(pavg1-pavg2)

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pdiag(i,j,k,1) = pdiag(i,j,k,1)+tem
        enddo
        enddo
        enddo

      !--------------------------------------------------------------!
      ! linear dynamic pressure:

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,dudz,dvdz)
        do k=2,nk
          dudz = -0.5*rf0(1,1,k)*(u0(1,1,k)-u0(1,1,k-1))*rdz*mf(1,1,k)
          dvdz = -0.5*rf0(1,1,k)*(v0(1,1,k)-v0(1,1,k-1))*rdz*mf(1,1,k)
          do j=1,nj+1
          do i=1,ni+1
            dum1(i,j,k) = dudz*(w3d(i-1,j,k)+w3d(i,j,k))
            dum2(i,j,k) = dvdz*(w3d(i,j-1,k)+w3d(i,j,k))
          enddo
          enddo
        enddo

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni+1
          dum1(i,j,1) = 0.0
          dum1(i,j,nk+1) = 0.0
          dum2(i,j,1) = 0.0
          dum2(i,j,nk+1) = 0.0
        enddo
        enddo

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,tem)
        do k=1,nk
          tem = 0.5*rr0(1,1,k)
          do j=1,nj+1
          do i=1,ni+1
            uten(i,j,k) = tem*(dum1(i,j,k)+dum1(i,j,k+1))
            vten(i,j,k) = tem*(dum2(i,j,k)+dum2(i,j,k+1))
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = 0.0
            divx(i,j,k) = 0.0
          enddo
          enddo
        enddo

        call   poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                     dum2,dum3,dum4,dum5,divx,pdiag(ibdp,jbdp,kbdp,2),uten,vten,wten,  &
                     cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans,dt)

        ! Adjust mean pressure ... set average value at upper-most level to zero:

        pavg2 = 0.0

        do j=1,nj
        do i=1,ni
          pavg2 = pavg2+pdiag(i,j,nk,2)
        enddo
        enddo

        pavg2 = pavg2/dble(nx*ny)

        tem = sngl(0.0d0-pavg2)

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pdiag(i,j,k,2) = pdiag(i,j,k,2)+tem
        enddo
        enddo
        enddo

      !--------------------------------------------------------------!
      ! nonlinear dynamic pressure:

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        DO k=1,nk
          ! store perturbation u,v in uten1,vten1 arrays:
          do j=jb,je
          do i=ib,ie+1
            rru(i,j,k)=rho0(1,1,k)*u3d(i,j,k)
            uten1(i,j,k)=u3d(i,j,k)-u0(i,j,k)
          enddo
          enddo
          do j=jb,je+1
          do i=ib,ie
            rrv(i,j,k)=rho0(1,1,k)*v3d(i,j,k)
            vten1(i,j,k)=v3d(i,j,k)-v0(i,j,k)
          enddo
          enddo
          IF(k.eq.1)THEN
            do j=0,nj+1
            do i=0,ni+1
              rrw(i,j,   1) = 0.0
              rrw(i,j,nk+1) = 0.0
            enddo
            enddo
          ELSE
            do j=0,nj+1
            do i=0,ni+1
              rrw(i,j,k)=rf0(1,1,k)*w3d(i,j,k)
            enddo
            enddo
          ENDIF
        ENDDO

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          divx(i,j,k)=( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )      &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
          if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          uten(i,j,k) = 0.0
          vten(i,j,k) = 0.0
          wten(i,j,k) = 0.0
        enddo
        enddo
        enddo

        call advu(nrkmax,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,gz,rgz,gzu,mh,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,divx, &
                   rru,uten1,uten,rrv,rrw,rdsf,c1,c2,rho, dt  ,.false.,udiag,wprof,bndy,kbdy,hflxw,hflxe,hflxs,hflxn,vf,mf,vten1,w3d)

        call advv(nrkmax,xh,rxh,arh1,arh2,uh,xf,vh,vf,gz,rgz,gzv,mh,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,divx, &
                   rru,rrv,vten1,vten,rrw,rdsf,c1,c2,rho, dt  ,.false.,vdiag,wprof,bndy,kbdy,hflxw,hflxe,hflxs,hflxn,uf,mf,uten1,w3d)

        call   advw(nrkmax,xh,rxh,arh1,arh2,uh,xf,vh,gz,rgz,mh,mf,rho0,rr0,rf0,rrf0,  &
                    dum1,dum2,dum3,dum4,dum5,dum6,dum7,divx,                       &
                    rru,rrv,rrw,w3d  ,wten,rds,rdsf,c1,c2,rho, dt  ,               &
                    .false.,wdiag,hadvordrv,vadvordrv,advwenov,bndy,kbdy,uf,vf,uten1,vten1,hflxw,hflxe,hflxs,hflxn,.true.)

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          divx(i,j,k) = 0.0
        enddo
        enddo
        enddo

        call   poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                     dum2,dum3,dum4,dum5,divx,pdiag(ibdp,jbdp,kbdp,3),uten,vten,wten,  &
                     cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans,dt)

        ! Adjust mean pressure ... set average value at upper-most level to zero:

        pavg2 = 0.0

        do j=1,nj
        do i=1,ni
          pavg2 = pavg2+pdiag(i,j,nk,3)
        enddo
        enddo

        pavg2 = pavg2/dble(nx*ny)

        tem = sngl(0.0d0-pavg2)

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pdiag(i,j,k,3) = pdiag(i,j,k,3)+tem
        enddo
        enddo
        enddo

      !--------------------------------------------------------------!
      ! Coriolis pressure:

      IF( icor.eq.1 )THEN

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          uten(i,j,k) =  0.125*(f2d(i,j)+f2d(i-1,j))           &
                              *( (v3d(i  ,j,k)+v3d(i  ,j+1,k)) &
                                +(v3d(i-1,j,k)+v3d(i-1,j+1,k)) )
          vten(i,j,k) = -0.125*(f2d(i,j)+f2d(i,j-1))           &
                              *( (u3d(i,j  ,k)+u3d(i+1,j  ,k)) &
                                +(u3d(i,j-1,k)+u3d(i+1,j-1,k)) )
          wten(i,j,k) = 0.0
          divx(i,j,k) = 0.0
        enddo
        enddo
        enddo

        call   poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                     dum2,dum3,dum4,dum5,divx,pdiag(ibdp,jbdp,kbdp,4),uten,vten,wten,  &
                     cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans,dt)

        ! Adjust mean pressure ... set average value at upper-most level to zero:

        pavg2 = 0.0

        do j=1,nj
        do i=1,ni
          pavg2 = pavg2+pdiag(i,j,nk,4)
        enddo
        enddo

        pavg2 = pavg2/dble(nx*ny)

        tem = sngl(0.0d0-pavg2)

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          pdiag(i,j,k,4) = pdiag(i,j,k,4)+tem
        enddo
        enddo
        enddo

      ENDIF

      !--------------------------------------------------------------!

      end subroutine pidcomp

  END MODULE pdcomp_module
