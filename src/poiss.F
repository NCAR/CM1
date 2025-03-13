  MODULE poiss_module

  implicit none

  private
  public :: poiss

  CONTAINS

      subroutine poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,   &
                       dum1,dum2,dum3,def,divx,ppi,uten,vten,wten,  &
                       cfb,cfa,cfc,d1,d2,pdt,lgbth,lgbph,rhs,trans,dttmp)

      use input
      use singleton , only : fft
      implicit none

      real, intent(in), dimension(ib:ie) :: uh
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rho0,rf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,def,ppi
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wten
      real, intent(in), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, intent(in), dimension(kpb:kpe) :: cfa,cfc,d1,d2
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,lgbth,lgbph
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      real, intent(in) :: dttmp

      integer :: i,j,k
      real :: tem
      real, dimension(nk) :: r1,deft2
      complex :: temc

      IF(axisymm.eq.1)THEN
        if(myid.eq.0)then
          print *
          print *,'  The anelastic/incompressible solver cannot be '
          print *,'  used with the axisymmetric model(yet)'
          print *
          print *,'  Stopping model ...'
          print *
        endif
        call stopcm1
      ENDIF

!-----------------------------------------------------------------------
!!!
!!!  Get the forcing from the exact numerical representation
!!!  of the right-hand side of the momentum equation
!!!
!!!  The divx term reduces accumulation of numerical errors over time
!!!
!!!
!!! fourier transform the total forcing, then forward part of tridiagonal solver:
!!!

      DO k=1,nk

        tem = 1.0/dttmp

!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          def(i,j,k)=rho0(1,1,k)*(                                  &
                     (uten(i+1,j,k)-uten(i,j,k))*rdx*uh(i)          &
                    +(vten(i,j+1,k)-vten(i,j,k))*rdy*vh(j) )        &
                    +( rf0(i,j,k+1)*wten(i,j,k+1)                   &
                      -rf0(i,j,k  )*wten(i,j,k  ) )*rdz*mh(i,j,k)   &
                    +divx(i,j,k)*tem
          rhs(i,j)=cmplx(def(i,j,k)*d1(k),0.0)
        enddo
        enddo

        if(imirror.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            rhs(ipe+1-i,j)=rhs(i,j)
          enddo
          enddo

        endif

        if(jmirror.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            rhs(i,jpe+1-j)=rhs(i,j)
          enddo
          enddo

        endif

        if(imirror.eq.1.and.jmirror.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            rhs(ipe+1-i,jpe+1-j)=rhs(i,j)
          enddo
          enddo

        endif

        trans=fft(rhs)

        IF( k.eq.1 )THEN

          !$omp parallel do default(shared)  &
          !$omp private(i,j,tem)
          do j=jpb,jpe
          do i=ipb,ipe
            tem = 1.0/cfb(i,j,1)
            lgbth(i,j,1) = -cfc(1)*tem
            lgbph(i,j,1) = trans(i,j)*tem
          enddo
          enddo

        ELSEIF( k.lt.nk )THEN

          !$omp parallel do default(shared)  &
          !$omp private(i,j,temc)
          do j=jpb,jpe
          do i=ipb,ipe
            temc = 1.0/(cfa(k)*lgbth(i,j,k-1)+cfb(i,j,k))
            lgbth(i,j,k) = -cfc(k)*temc
            lgbph(i,j,k) = (trans(i,j)-cfa(k)*lgbph(i,j,k-1))*temc
          enddo
          enddo

        ELSE

          !$omp parallel do default(shared)  &
          !$omp private(i,j,temc)
          do j=jpb,jpe
          do i=ipb,ipe
            temc = 1.0/(cfa(k)*lgbth(i,j,k-1)+cfb(i,j,k))
            lgbth(i,j,k) = -cfc(k)*temc
            lgbph(i,j,k) = (trans(i,j)-cfa(k)*lgbph(i,j,k-1))*temc
            pdt(i,j,k) = lgbph(i,j,k)
          enddo
          enddo

        ENDIF

        deft2(k) = real(trans(1,1))

      ENDDO

!-----------------------------------------------------------------------
!!! backward part of tridiagonal solver:

      DO k=(nk-1),1,-1

        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=jpb,jpe
        do i=ipb,ipe
          pdt(i,j,k) = lgbth(i,j,k)*pdt(i,j,k+1)+lgbph(i,j,k)
        enddo
        enddo

      ENDDO

      r1(1) = 0.0
      pdt(1,1,1) = 0.0

      r1(2) = (deft2(1)-cfb(1,1,1)*r1(1))/cfc(1)
      pdt(1,1,2) = cmplx( r1(2) , 0.0 )

      do k=2,(nk-1)
        r1(k+1) = (deft2(k)-cfa(k)*r1(k-1)-cfb(1,1,k)*r1(k))/cfc(k)
        pdt(1,1,k+1) = cmplx( r1(k+1) , 0.0 )
      enddo

!---------------------------------------------

!!!
!!! reverse fourier transform:
!!!

      DO k=1,nk

        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=jpb,jpe
        do i=ipb,ipe
          rhs(i,j)=pdt(i,j,k)
        enddo
        enddo

        trans=fft(rhs,inv=.true.)

        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          ppi(i,j,k)=real(trans(i,j))*d2(k)
        enddo
        enddo

      ENDDO

      if(timestats.ge.1) time_poiss=time_poiss+mytime()

!-----------------------------------------------------------------------
!  All done.

      end subroutine poiss


  END MODULE poiss_module
