  MODULE diff2_module

  implicit none

  private
  public :: diff2u,diff2v,diff2w,diff2s,diff2def

  CONTAINS


      subroutine diff2u(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                        dum,diffx,diffy,diffz,uten,ust,rho,rr,rf,divx,t11,t12,t13, &
                        doubud,udiag)
      use input
      use constants
      implicit none
 
      real, intent(in), dimension(ib:ie) :: rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,diffx,diffy,diffz
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in), dimension(ib:ie,jb:je) :: ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf,divx,t11,t12,t13
      logical, intent(in) :: doubud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag

      integer :: i,j,k,i2

      real, parameter :: twdth = 2.0/3.0

!--------------------------
!
!  2nd order artificial diffusion
!
!-----------------------------------------------------------------------

      IF(terrain_flag)THEN
        if(myid.eq.0) print *,'  cannot use diff2u with terrain (yet) '
        call stopcm1
      ENDIF

!-----------------------------------------------------------------------

    IF(axisymm.eq.0)THEN
      ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        !  x-direction
        do j=1,nj
        do i=0,ni+1
          dum(i,j,k)=2.0*t11(i,j,k)-rho(i,j,k)*twdth*divx(i,j,k)
        enddo
        enddo
        do j=1,nj
        do i=1,ni+1
          diffx(i,j,k)=(dum(i,j,k)-dum(i-1,j,k))*rdx*uf(i)
        enddo
        enddo
        !  y-direction
        do j=1,nj
        do i=1,ni+1
          diffy(i,j,k)=(t12(i,j+1,k)-t12(i,j,k))*rdy*vh(j)
        enddo
        enddo
      enddo

    ELSE
      ! axisymmetric grid:

      i2 = ni+1
      IF( ebc.eq.3 .or. ebc.eq.4 ) i2 = ni

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        do i=2,i2
          diffx(i,j,k)=2.0*(arf2(i)*arf2(i)*t11(i,j,k)-arf1(i)*arf1(i)*t11(i-1,j,k))*rdx*uf(i) &
                      -twdth*(rho(i,j,k)*divx(i,j,k)-rho(i-1,j,k)*divx(i-1,j,k))*rdx*uf(i)
        enddo
        IF( ebc.eq.3 .or. ebc.eq.4 ) diffx(ni+1,j,k) = 0.0
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO j=1,nj

        do k=1,nk
        do i=1,ni+1
          diffz(i,j,k)=(t13(i,j,k+1)-t13(i,j,k))*rdz*mh(1,1,k)
        enddo
        enddo

        IF( axisymm.eq.1 .and. (ebc.eq.3 .or. ebc.eq.4) )THEN
          do k=1,nk
            diffz(ni+1,j,k)=0.0
          enddo
        ENDIF

      ENDDO

!-----------------------------------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten(i,j,k)=uten(i,j,k)+((diffx(i,j,k)+diffy(i,j,k))+diffz(i,j,k))  &
                                 /(0.5*(rho(i-1,j,k)+rho(i,j,k)))
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=2,ni+1
          uten(i,j,k)=uten(i,j,k)+(diffx(i,j,k)+diffz(i,j,k))  &
                                 /(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
        enddo
        enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------
!  Diagnostics:

    IF( doubud .and. ud_hediff.ge.1 .and. ud_vediff.ge.1 )THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          udiag(i,j,k,ud_hediff)=(diffx(i,j,k)+diffy(i,j,k))/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          udiag(i,j,k,ud_vediff)=diffz(i,j,k)/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=2,ni+1
          udiag(i,j,k,ud_hediff)=diffx(i,j,k)/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          udiag(i,j,k,ud_vediff)=diffz(i,j,k)/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()

      end subroutine diff2u


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine diff2v(xh,arh1,arh2,uh,rxf,arf1,arf2,uf,vh,vf,mh,mf,  &
                        dum,diffx,diffy,diffz,vten,ust,rho,rr,rf,divx,t22,t12,t23, &
                        dovbud,vdiag)
      use input
      use constants
      implicit none
 
      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: rxf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,diffx,diffy,diffz
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in), dimension(ib:ie,jb:je) :: ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf,divx,t22,t12,t23
      logical, intent(in) :: dovbud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
 
      integer :: i,j,k,j1,j2

      real, parameter :: twdth = 2.0/3.0
 
!--------------------------
!
!  2nd order artificial diffusion
!
!-----------------------------------------------------------------------

      IF(terrain_flag)THEN
        if(myid.eq.0) print *,'  cannot use diff2v with terrain (yet) '
        call stopcm1
      ENDIF

!-----------------------------------------------------------------------

    IF(axisymm.eq.0)THEN
      ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj+1
        do i=1,ni
          diffx(i,j,k)=(t12(i+1,j,k)-t12(i,j,k))*rdx*uh(i)
        enddo
        enddo
        !  y-direction
        do j=0,nj+1
        do i=1,ni
          dum(i,j,k)=2.0*t22(i,j,k)-rho(i,j,k)*twdth*divx(i,j,k)
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          diffy(i,j,k)=(dum(i,j,k)-dum(i,j-1,k))*rdy*vf(j)
        enddo
        enddo
      enddo

    ELSE
      ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(j,k)
      do k=1,nk
      do j=1,nj
        do i=1,ni
          diffx(i,j,k)=(arh2(i)*arh2(i)*t12(i+1,j,k)-arh1(i)*arh1(i)*t12(i,j,k))*rdx*uh(i)
        enddo
      enddo
      enddo

    ENDIF
 
!-----------------------------------------------------------------------
!  z-direction

      j1 = 1
      j2 = nj+1
      IF( axisymm.eq.1 ) j2 = nj

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO j=j1,j2

        do k=1,nk
        do i=1,ni
          diffz(i,j,k)=(t23(i,j,k+1)-t23(i,j,k))*rdz*mh(1,1,k)
        enddo
        enddo

      ENDDO

!-----------------------------------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten(i,j,k)=vten(i,j,k)+((diffx(i,j,k)+diffy(i,j,k))+diffz(i,j,k))  &
                                 /(0.5*(rho(i,j-1,k)+rho(i,j,k)))
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          vten(i,j,k)=vten(i,j,k)+(diffx(i,j,k)+diffz(i,j,k))*rr(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------
!  Diagnostics:

    IF( dovbud .and. vd_hediff.ge.1 .and. vd_vediff.ge.1 )THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vdiag(i,j,k,vd_hediff)=(diffx(i,j,k)+diffy(i,j,k))/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          vdiag(i,j,k,vd_vediff)=diffz(i,j,k)/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          vdiag(i,j,k,vd_hediff)=diffx(i,j,k)*rr(i,j,k)
          vdiag(i,j,k,vd_vediff)=diffz(i,j,k)*rr(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      end subroutine diff2v


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine diff2w(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                        dum,diffx,diffy,diffz,wten,rho,rr,rf,divx,t33,t13,t23, &
                        dowbud,wdiag)
      use input
      use constants
      implicit none
 
      real, intent(in), dimension(ib:ie) :: rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,diffx,diffy,diffz
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf,divx,t33,t13,t23
      logical, intent(in) :: dowbud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
 
      integer :: i,j,k

      real, parameter :: twdth = 2.0/3.0

!--------------------------
!
!  2nd order artificial diffusion
!
!-----------------------------------------------------------------------

      IF(terrain_flag)THEN
        if(myid.eq.0) print *,'  cannot use diff2w with terrain (yet) '
        call stopcm1
      ENDIF

!-----------------------------------------------------------------------

    IF(axisymm.eq.0)THEN
      ! Cartesian grid:

      !  x-direction
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
        do j=1,nj
        do i=1,ni
          diffx(i,j,k)=(t13(i+1,j,k)-t13(i,j,k))*rdx*uh(i)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          diffy(i,j,k)=(t23(i,j+1,k)-t23(i,j,k))*rdy*vh(j)
        enddo
        enddo
      enddo

    ELSE
      ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
        do i=1,ni
          diffx(i,j,k)=(arh2(i)*t13(i+1,j,k)-arh1(i)*t13(i,j,k))*rdx*uh(i)
        enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO j=1,nj

        do k=1,nk
        do i=1,ni
          dum(i,j,k)=2.0*t33(i,j,k)-rho(i,j,k)*twdth*divx(i,j,k)
        enddo
        enddo

        do k=2,nk
        do i=1,ni
          diffz(i,j,k)=(dum(i,j,k)-dum(i,j,k-1))*rdz*mf(1,1,k)
        enddo
        enddo

      ENDDO

!-----------------------------------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=wten(i,j,k)+((diffx(i,j,k)+diffy(i,j,k))+diffz(i,j,k))/rf(i,j,k)
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=wten(i,j,k)+(diffx(i,j,k)+diffz(i,j,k))/rf(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------
!  Diagnostics:

    IF( dowbud .and. wd_hediff.ge.1 .and. wd_vediff.ge.1 )THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_hediff)=(diffx(i,j,k)+diffy(i,j,k))/rf(i,j,k)
          wdiag(i,j,k,wd_vediff)=diffz(i,j,k)/rf(i,j,k)
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_hediff)=diffx(i,j,k)/rf(i,j,k)
          wdiag(i,j,k,wd_vediff)=diffz(i,j,k)/rf(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      end subroutine diff2w


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
 
      subroutine diff2s(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                        dum,diffx,diffy,diffz,s,sten,rho,rr,rf,  &
                        dobud,ibd,ied,jbd,jed,kbd,ked,ndiag,diag,sd_hediff,sd_vediff)
      use input
      use constants
      implicit none
 
      real, dimension(ib:ie) :: rxh,arh1,arh2,uh
      real, dimension(ib:ie+1) :: xf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,diffx,diffy,diffz
      real, dimension(ib:ie,jb:je,kb:ke) :: s,sten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      logical, intent(in) :: dobud
      integer, intent(in) :: ibd,ied,jbd,jed,kbd,ked,ndiag,sd_hediff,sd_vediff
      real, intent(inout) , dimension(ibd:ied,jbd:jed,kbd:ked,ndiag) :: diag
 
      integer :: i,j,k
      real :: coef

!--------------------------
!
!  2nd order artificial diffusion
!
!-----------------------------------------------------------------------

      IF(terrain_flag)THEN
        if(myid.eq.0) print *,'  cannot use diff2s with terrain (yet) '
        call stopcm1
      ENDIF

      coef = kdiff2

!-----------------------------------------------------------------------
!  x-direction
!  y-direction

    IF(axisymm.eq.0)THEN
      ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni+1
          dum(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*(s(i,j,k)-s(i-1,j,k))*rdx*uf(i)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          diffx(i,j,k)=coef*(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          dum(i,j,k)=0.5*(rho(i,j,k)+rho(i,j-1,k))*(s(i,j,k)-s(i,j-1,k))*rdy*vf(j)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          diffy(i,j,k)=coef*(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
        enddo
        enddo
      enddo

    ELSE
      ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
        dum(1,j,k)=0.0
        do i=2,ni+1
          dum(i,j,k)=(s(i,j,k)-s(i-1,j,k))*rdx*uf(i)   &
                *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
        enddo
        !-----
        ! assume zero flux:
        if(wbc.eq.3.or.wbc.eq.4)then
          dum(1,j,k) = 0.0
        endif
        if(ebc.eq.3.or.ebc.eq.4)then
          dum(ni+1,j,k) = 0.0
        endif
        !-----
        do i=1,ni
          diffx(i,j,k)=coef*(arh2(i)*dum(i+1,j,k)-arh1(i)*dum(i,j,k))*rdx*uh(i)
        enddo
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO j=1,nj

        do k=2,nk
        do i=1,ni
          dum(i,j,k)=(s(i,j,k)-s(i,j,k-1))*rdz*mf(1,1,k)
        enddo
        enddo

        ! assume zero flux boundary conditions:
        do i=1,ni
          dum(i,j,1)=0.0
          dum(i,j,nk+1)=0.0
        enddo

        do k=1,nk
        do i=1,ni
          diffz(i,j,k)=coef*(rf(i,j,k+1)*dum(i,j,k+1)-rf(i,j,k)*dum(i,j,k))*rdz*mh(1,1,k)
        enddo
        enddo

      ENDDO

!-----------------------------------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          sten(i,j,k)=sten(i,j,k)+((diffx(i,j,k)+diffy(i,j,k))+diffz(i,j,k))*rr(i,j,k)
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          sten(i,j,k)=sten(i,j,k)+(diffx(i,j,k)+diffz(i,j,k))*rr(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------
!  Diagnostics:

    IF( dobud .and. sd_hediff.ge.1 .and. sd_vediff.ge.1 )THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          diag(i,j,k,sd_hediff)=(diffx(i,j,k)+diffy(i,j,k))*rr(i,j,k)
          diag(i,j,k,sd_vediff)=diffz(i,j,k)*rr(i,j,k)
        enddo
        enddo
        enddo

      ELSE
        ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          diag(i,j,k,sd_hediff)=diffx(i,j,k)*rr(i,j,k)
          diag(i,j,k,sd_vediff)=diffz(i,j,k)*rr(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      end subroutine diff2s


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine diff2def(uh,arh1,arh2,uf,arf1,arf2,vh,vf,mh,c1,c2,mf,ust,znt,u1,v1,s1,  &
                          divx,rho,rr,rf,t11,t12,t13,t22,t23,t33,ua,va,wa,dissten)
      use input
      use constants
      implicit none
 
      real, intent(in), dimension(ib:ie) :: uh,arh1,arh2
      real, intent(in), dimension(ib:ie+1) :: uf,arf1,arf2
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je) :: ust,znt,u1,v1,s1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rr,rf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in),    dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in),    dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: dissten

      integer :: i,j,k
      real :: coef,rcoef,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,tem

      real, parameter :: twodthree  = 2.0/3.0

!--------------------------------------------------------------------
!  NOTE:  diff2 arrays all assume no terrain
!--------------------------------------------------------------------

      coef = kdiff2
      rcoef = 1.0/coef

!-----------------------------------------------------------------------
!  Cartesian grid:

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        DO k=1,nk

          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)        &
                       +(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)        &
                       +(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo

          do j=0,nj+1
          do i=0,ni+1
            t11(i,j,k) = coef*rho(i,j,k)*(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
            t22(i,j,k) = coef*rho(i,j,k)*(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
            t33(i,j,k) = coef*rho(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)
            t12(i,j,k) = coef*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)    &
                               +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) )  &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
          enddo
          enddo

          !-----
          ! lateral boundary conditions:
          if(wbc.eq.3.and.ibw.eq.1)then
            ! free slip b.c.
            do j=1,nj+1
              t12(1,j,k) = t12(2,j,k)
            enddo
          elseif(wbc.eq.4.and.ibw.eq.1)then
            ! no slip b.c.
            i = 1
            do j=1,nj+1
              t12(1,j,k) = coef*2.0*va(1,j,k)*rdx*uf(1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          if(ebc.eq.3.and.ibe.eq.1)then
            ! free slip b.c.
            do j=1,nj+1
              t12(ni+1,j,k) = t12(ni,j,k)
            enddo
          elseif(ebc.eq.4.and.ibe.eq.1)then
            ! no slip b.c.
            i = ni+1
            do j=1,nj+1
              t12(ni+1,j,k) = -coef*2.0*va(ni,j,k)*rdx*uf(ni+1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          !-----
          !-----
          if(sbc.eq.3.and.ibs.eq.1)then
            ! free slip b.c.
            do i=1,ni+1
              t12(i,1,k) = t12(i,2,k)
            enddo
          elseif(sbc.eq.4.and.ibs.eq.1)then
            ! no slip b.c.
            j = 1
            do i=1,ni+1
              t12(i,1,k) = coef*2.0*ua(i,1,k)*rdy*vf(1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          if(nbc.eq.3.and.ibn.eq.1)then
            ! free slip b.c.
            do i=1,ni+1
              t12(i,nj+1,k) = t12(i,nj,k)
            enddo
          elseif(nbc.eq.4.and.ibn.eq.1)then
            ! no slip b.c.
            j = nj+1
            do i=1,ni+1
              t12(i,nj+1,k) = -coef*2.0*ua(i,nj,k)*rdy*vf(nj+1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
          endif
          !-----

          IF(k.ne.1)THEN
            do j=1,nj+1
            do i=1,ni+1
              t13(i,j,k)=coef*( (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)  &
                               +(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i) )
              t23(i,j,k)=coef*( (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)  &
                               +(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j) )
            enddo
            enddo
            !-----
            ! lateral boundary conditions:
            if(wbc.eq.3.and.ibw.eq.1)then
              ! free slip b.c.
              do j=1,nj
                t13(1,j,k) = t13(2,j,k)
              enddo
            elseif(wbc.eq.4.and.ibw.eq.1)then
              ! no slip b.c.
              do j=1,nj
                t13(1,j,k) = coef*2.0*wa(1,j,k)*rdx*uf(1)
              enddo
            endif
            if(ebc.eq.3.and.ibe.eq.1)then
              ! free slip b.c.
              do j=1,nj
                t13(ni+1,j,k) = t13(ni,j,k)
              enddo
            elseif(ebc.eq.4.and.ibe.eq.1)then
              ! no slip b.c.
              do j=1,nj
                t13(ni+1,j,k) = -coef*2.0*wa(ni,j,k)*rdx*uf(ni+1)
              enddo
            endif
            !-----
            !-----
            if(sbc.eq.3.and.ibs.eq.1)then
              ! free slip b.c.
              do i=1,ni
                t23(i,1,k) = t23(i,2,k)
              enddo
            elseif(sbc.eq.4.and.ibs.eq.1)then
              ! no slip b.c.
              do i=1,ni
                t23(i,1,k) = coef*2.0*wa(i,1,k)*rdy*vf(1)
              enddo
            endif
            if(nbc.eq.3.and.ibn.eq.1)then
              ! free slip b.c.
              do i=1,ni
                t23(i,nj+1,k) = t23(i,nj,k)
              enddo
            elseif(nbc.eq.4.and.ibn.eq.1)then
              ! no slip b.c.
              do i=1,ni
                t23(i,nj+1,k) = -coef*2.0*wa(i,nj,k)*rdy*vf(nj+1)
              enddo
            endif
            !-----
          ENDIF  ! endif for k.ne.1

        ENDDO  ! enddo for k loop

        !---------------------------------
        ! bottom/top boundary conditions:
        !---------------------------------

        if(bbc.eq.1)then
          ! free slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,1)=t13(i,j,2)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            t23(i,j,1)=t23(i,j,2)
          enddo
          enddo
        elseif(bbc.eq.2)then
          ! no slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,1)=coef*2.0*(ua(i,j,1)+umove)*rdz*mf(1,1,1)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            t23(i,j,1)=coef*2.0*(va(i,j,1)+vmove)*rdz*mf(1,1,1)
          enddo
          enddo
        elseif(bbc.eq.3)then
          ! semi-slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,1) = 0.5*( (ust(i-1,j)**2)*(u1(i-1,j)/max(s1(i-1,j),0.01)) &
                              +(ust(i  ,j)**2)*(u1(i  ,j)/max(s1(i  ,j),0.01)) )
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            t23(i,j,1) = 0.5*( (ust(i,j-1)**2)*(v1(i,j-1)/max(s1(i,j-1),0.01)) &
                              +(ust(i,j  )**2)*(v1(i,j  )/max(s1(i,j  ),0.01)) )
          enddo
          enddo
        endif

        if(tbc.eq.1)then
          ! free slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,nk+1)=t13(i,j,nk)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            t23(i,j,nk+1)=t23(i,j,nk)
          enddo
          enddo
        elseif(tbc.eq.2)then
          ! no slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,nk+1)=-coef*2.0*ua(i,j,nk)*rdz*mf(1,1,nk+1)
          enddo
          enddo
          do j=1,nj+1
          do i=1,ni
            t23(i,j,nk+1)=-coef*2.0*va(i,j,nk)*rdz*mf(1,1,nk+1)
          enddo
          enddo
        endif

        ! mult t13,t23 by density:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,k)=t13(i,j,k)*0.5*(rf(i-1,j,k)+rf(i,j,k))
          t23(i,j,k)=t23(i,j,k)*0.5*(rf(i,j-1,k)+rf(i,j,k))
        enddo
        enddo
        enddo

        ! get dissipation rate:
        IF(idiss.eq.1.or.output_dissten.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,tem)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            tmp11=( c1(i,j,k)*t11(i,j,k-1)**2 + c2(i,j,k)*t11(i,j,k)**2 )
            tmp22=( c1(i,j,k)*t22(i,j,k-1)**2 + c2(i,j,k)*t22(i,j,k)**2 )
            tmp33=( c1(i,j,k)*t33(i,j,k-1)**2 + c2(i,j,k)*t33(i,j,k)**2 )
            tmp12=0.25*( c1(i,j,k)*( ( t12(i,j  ,k-1)**2 + t12(i+1,j+1,k-1)**2 )     &
                                   + ( t12(i,j+1,k-1)**2 + t12(i+1,j  ,k-1)**2 ) )   &
                        +c2(i,j,k)*( ( t12(i,j  ,k  )**2 + t12(i+1,j+1,k  )**2 )     &
                                   + ( t12(i,j+1,k  )**2 + t12(i+1,j  ,k  )**2 ) ) ) 
            tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )
            tmp23=0.5*( t23(i,j,k)**2 + t23(i,j+1,k)**2 )
            dissten(i,j,k)= rcoef*( ( 2.0*( tmp33 ) + ( tmp13 + tmp23 )               &
                                     +2.0*( tmp11 + tmp22 ) + tmp12 )/(rf(i,j,k)**2)  &
                   -twodthree*(divx(i,j,k)**2) )
          enddo
          enddo
          enddo
        ENDIF

!-----------------------------------------------------------------------
!  axisymmetric grid:

      ELSEIF(axisymm.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        DO k=1,nk

          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(arh2(i)*ua(i+1,j,k)-arh1(i)*ua(i,j,k))*rdx*uh(i)   &
                       +(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)
            if(abs(divx(i,j,k)).lt.smeps) divx(i,j,k)=0.0
          enddo
          enddo

          do j=1,nj
          do i=0,ni+1
            t11(i,j,k) = coef*rho(i,j,k)*(arf1(i+1)*ua(i+1,j,k)-arf2(i)*ua(i,j,k))*rdx*uh(i)
            t33(i,j,k) = coef*rho(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)
            t12(i,j,k) = coef*(arh1(i)*va(i,j,k)-arh2(i-1)*va(i-1,j,k))*rdx*uf(i)  &
                        *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
          enddo
          enddo

          !-----
          ! lateral boundary conditions:
          j = 1
          if(wbc.eq.3)then
            ! free slip b.c.
            t12(1,j,k) = 0.0
          elseif(wbc.eq.4)then
            ! no slip b.c.
            i = 1
            t12(1,j,k) = coef*2.0*va(1,j,k)*rdx*uf(1)   &
                      *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
          endif
          if(ebc.eq.3)then
            ! free slip b.c.
            t12(ni+1,j,k) = t12(ni,j,k)
          elseif(ebc.eq.4)then
            ! no slip b.c.
            i = ni+1
            t12(ni+1,j,k) = -coef*2.0*va(ni,j,k)*rdx*uf(ni+1)   &
                      *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
          endif
          !-----

          IF(k.ne.1)THEN
            do j=1,nj
            do i=1,ni+1
              t13(i,j,k)=coef*( (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)  &
                               +(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i) )
              t23(i,j,k)=coef*(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
            enddo
            enddo
            !-----
            ! lateral boundary conditions:
            j = 1
            if(wbc.eq.3)then
              ! free slip b.c.
              t13(1,j,k) = 0.0
            elseif(wbc.eq.4)then
              ! no slip b.c.
              t13(1,j,k) = coef*2.0*wa(1,j,k)*rdx*uf(1)
            endif
            if(ebc.eq.3)then
              ! free slip b.c.
              t13(ni+1,j,k) = t13(ni,j,k)
            elseif(ebc.eq.4)then
              ! no slip b.c.
              t13(ni+1,j,k) = -coef*2.0*wa(ni,j,k)*rdx*uf(ni+1)
            endif
            !-----
          ENDIF  ! endif for k.ne.1

        ENDDO  ! enddo for k-loop

        !---------------------------------
        ! bottom/top boundary conditions:
        !---------------------------------

        if(bbc.eq.1)then
          ! free slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,1)=t13(i,j,2)
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            t23(i,j,1)=t23(i,j,2)
          enddo
          enddo
        elseif(bbc.eq.2)then
          ! no slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,1)=coef*2.0*(ua(i,j,1)+umove)*rdz*mf(1,1,1)
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            t23(i,j,1)=coef*2.0*(va(i,j,1)+vmove)*rdz*mf(1,1,1)
          enddo
          enddo
        elseif(bbc.eq.3)then
          ! semi-slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,1) = 0.5*( arf1(i)*(ust(i-1,j)**2)*(u1(i-1,j)/max(s1(i-1,j),0.01)) &
                              +arf2(i)*(ust(i  ,j)**2)*(u1(i  ,j)/max(s1(i  ,j),0.01)) )
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            t23(i,j,1) = (ust(i,j)**2)*(v1(i,j)/max(s1(i,j),0.01))
          enddo
          enddo
        endif

        if(tbc.eq.1)then
          ! free slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,nk+1)=t13(i,j,nk)
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            t23(i,j,nk+1)=t23(i,j,nk)
          enddo
          enddo
        elseif(tbc.eq.2)then
          ! no slip b.c.
          do j=1,nj
          do i=1,ni+1
            t13(i,j,nk+1)=-coef*2.0*ua(i,j,nk)*rdz*mf(1,1,nk+1)
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            t23(i,j,nk+1)=-coef*2.0*va(i,j,nk)*rdz*mf(1,1,nk+1)
          enddo
          enddo
        endif

        ! mult t13,t23 by density:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni+1
          t13(i,j,k)=t13(i,j,k)*0.5*(arf1(i)*rf(i-1,j,k)+arf2(i)*rf(i,j,k))
          t23(i,j,k)=t23(i,j,k)*rf(i,j,k)
        enddo
        enddo
        enddo

        ! get dissipation rate:
        IF(idiss.eq.1.or.output_dissten.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k,tmp11,tmp33,tmp12,tmp13,tmp23,tem)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            tmp11=( c1(1,1,k)*t11(i,j,k-1)**2 + c2(1,1,k)*t11(i,j,k)**2 )
            tmp33=( c1(1,1,k)*t33(i,j,k-1)**2 + c2(1,1,k)*t33(i,j,k)**2 )
            tmp12=0.5*( c1(1,1,k)*( t12(i,j  ,k-1)**2 + t12(i+1,j  ,k-1)**2 ) &
                       +c2(1,1,k)*( t12(i,j  ,k  )**2 + t12(i+1,j  ,k  )**2 ) ) 
            tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )
            tmp23=      t23(i,j,k)**2
            dissten(i,j,k)= rcoef*( ( 2.0*( tmp33 ) + ( tmp13 + tmp23 )               &
                                     +2.0*( tmp11 ) + tmp12 )/(rf(i,j,k)**2)  &
                   -twodthree*(divx(i,j,k)**2) )
          enddo
          enddo
          enddo
        ENDIF

!-----------------------------------------------------------------------

      ELSE
        stop 333
      ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_diffu=time_diffu+mytime()
 
      end subroutine diff2def


  END MODULE diff2_module
