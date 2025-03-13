  MODULE bc_module

  implicit none

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bc2d(s)
      use input
      implicit none

      real, dimension(ib:ie,jb:je) :: s

      integer i,j

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jb,je
        do i=1,ngxy
          s(1-i,j)=s(ni+1-i,j)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jb,je
        do i=1-ngxy,0
          s(i,j)=s(1,j)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jb,je
        do i=1,ngxy
          s(1-i,j)=s(i,j)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jb,je
        do i=1,ngxy
          s(ni+i,j)=s(i,j)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jb,je
        do i=ni+1,ni+ngxy
          s(i,j)=s(ni,j)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=jb,je
        do i=1,ngxy
          s(ni+i,j)=s(ni+1-i,j)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=ib,ie
          s(i,1-j)=s(i,nj+1-j)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1-ngxy,0
        do i=ib,ie
          s(i,j)=s(i,1)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=ib,ie
          s(i,1-j)=s(i,j)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=ib,ie
          s(i,nj+j)=s(i,j)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=nj+1,nj+ngxy
        do i=ib,ie
          s(i,j)=s(i,nj)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=ib,ie
          s(i,nj+j)=s(i,nj+1-j)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bc2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcs(s)
      use input
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: s

      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=1,ngxy
          s(1-i,j,k)=s(ni+1-i,j,k)
        enddo
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1-ngxy,0
          s(i,j,k)=s(1,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ngxy
          s(1-i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=1,ngxy
          s(ni+i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=ni+1,ni+ngxy
          s(i,j,k)=s(ni,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ngxy
          s(ni+i,j,k)=s(ni+1-i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=ib,ie
          s(i,1-j,k)=s(i,nj+1-j,k)
        enddo
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1-ngxy,0
        do i=0,ni+1
          s(i,j,k)=s(i,1,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=0,ni+1
          s(i,1-j,k)=s(i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=ib,ie
          s(i,nj+j,k)=s(i,j,k)
        enddo
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=nj+1,nj+ngxy
        do i=0,ni+1
          s(i,j,k)=s(i,nj,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=0,ni+1
          s(i,nj+j,k)=s(i,nj+1-j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcs


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcp(p)
      use input
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: p

      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared) private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(0,j,k)=p(ni,j,k)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(0,j,k)=p(1,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
        !$omp parallel do default(shared) private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(0,j,k)=p(1,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
        !$omp parallel do default(shared) private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(ni+1,j,k)=p(1,j,k)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
        !$omp parallel do default(shared) private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(ni+1,j,k)=p(ni,j,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
        !$omp parallel do default(shared) private(j,k)
        do k=1,nk
        do j=0,nj+1
          p(ni+1,j,k)=p(ni,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
        !$omp parallel do default(shared) private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,0,k)=p(i,nj,k)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,0,k)=p(i,1,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
        !$omp parallel do default(shared) private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,0,k)=p(i,1,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
        !$omp parallel do default(shared) private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,nj+1,k)=p(i,1,k)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,nj+1,k)=p(i,nj,k)
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
        !$omp parallel do default(shared) private(i,k)
        do k=1,nk
        do i=0,ni+1
          p(i,nj+1,k)=p(i,nj,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcu(u)
      use input
      implicit none

      real, dimension(ib:ie+1,jb:je,kb:ke) :: u

      integer i,j,k  

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=1,ngxy
          u(1-i,j,k)=u(ni+1-i,j,k)
        enddo
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1-ngxy,0
          u(i,j,k)=u(1,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
          u(1,j,k)=0.
        do i=1,ngxy
          u(1-i,j,k)=-u(1+i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=1,ngxy
          u(ni+1+i,j,k)=u(1+i,j,k)
        enddo
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=ni+2,ni+1+ngxy
          u(i,j,k)=u(ni+1,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
          u(ni+1,j,k)=0.0
        do i=1,ngxy
          u(ni+1+i,j,k)=-u(ni+1-i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=ib,ie+1
          u(i,1-j,k)=u(i,nj+1-j,k)
        enddo
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1-ngxy,0
        do i=0,ni+2
          u(i,j,k)=u(i,1,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=0,ni+2
          u(i,1-j,k)=u(i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=ib,ie+1
          u(i,nj+j,k)=u(i,j,k)
        enddo
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=nj+1,nj+ngxy
        do i=0,ni+2
          u(i,j,k)=u(i,nj,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=0,ni+2
          u(i,nj+j,k)=u(i,nj+1-j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine bcu


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcv(v)
      use input
      implicit none

      real, dimension(ib:ie,jb:je+1,kb:ke) :: v
 
      integer i,j,k

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=ib,ie
          v(i,1-j,k)=v(i,nj+1-j,k)
        enddo
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1-ngxy,0
        do i=0,ni+1
          v(i,j,k)=v(i,1,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=0,ni+1
          v(i,1-j,k)=-v(i,j+1,k)
        enddo
        enddo
        do i=0,ni+1
          v(i,   1,k)=0.
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

      IF(axisymm.eq.1)THEN
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do i=ib,ie
          v(i,2,k)=v(i,1,k)
        enddo
        enddo
      ENDIF

#ifndef MPI
      if(nbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=ib,ie
          v(i,nj+1+j,k)=v(i,1+j,k)
        enddo
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=nj+2,nj+1+ngxy
        do i=0,ni+1
          v(i,j,k)=v(i,nj+1,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,ngxy
        do i=0,ni+1
          v(i,nj+1+j,k)=-v(i,nj+1-j,k)
        enddo
        enddo
        do i=0,ni+1
          v(i,nj+1,k)=0.
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je+1
        do i=1,ngxy
          v(1-i,j,k)=v(ni+1-i,j,k)
        enddo
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+2
        do i=1-ngxy,0
          v(i,j,k)=v(1,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+2
        do i=1,ngxy
          v(1-i,j,k)=v(i,j,k)
        enddo
        enddo
        enddo
      endif

      IF(axisymm.eq.1)THEN
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+2
        do i=1,ngxy
          v(1-i,j,k) = -v(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je+1
        do i=1,ngxy
          v(ni+i,j,k)=v(i,j,k)
        enddo
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+2
        do i=ni+1,ni+ngxy
          v(i,j,k)=v(ni,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+2
        do i=1,ngxy
          v(ni+i,j,k)=v(ni+1-i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine bcv


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcw(w,flag)
      use input
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke+1) :: w
      integer flag
 
      integer i,j,k

!-----------------------------------------------------------------------
!  west boundary condition
 
#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=jb,je
        do i=1,ngxy
          w(1-i,j,k)=w(ni+1-i,j,k)
        enddo
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=1-ngxy,0
          w(i,j,k)=w(1,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=1,ngxy
          w(1-i,j,k)=w(i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=jb,je
        do i=1,ngxy
          w(ni+i,j,k)=w(i,j,k)
        enddo
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=ni+1,ni+ngxy
          w(i,j,k)=w(ni,j,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=1,ngxy
          w(ni+i,j,k)=w(ni+1-i,j,k)
        enddo
        enddo
        enddo 
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,ngxy
        do i=ib,ie
          w(i,1-j,k)=w(i,nj+1-j,k)
        enddo
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1-ngxy,0
        do i=0,ni+1
          w(i,j,k)=w(i,1,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,ngxy
        do i=0,ni+1
          w(i,1-j,k)=w(i,j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,ngxy
        do i=ib,ie
          w(i,nj+j,k)=w(i,j,k)
        enddo
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=nj+1,nj+ngxy
        do i=0,ni+1
          w(i,j,k)=w(i,nj,k)
        enddo
        enddo
        enddo
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,ngxy
        do i=0,ni+1
          w(i,nj+j,k)=w(i,nj+1-j,k)
        enddo
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  top/bottom boundary condition

    IF(flag.eq.1)THEN

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        w(i,j, 1)=0.0
        w(i,j,nk+1)=0.0
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine bcw


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcwsfc(gz,dzdx,dzdy,u,v,w)
      use input
      implicit none

      real, intent(in), dimension(itb:ite,jtb:jte) :: gz
      real, intent(in), dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w

      integer :: i,j

!-----------------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        w(i,j,1) = 0.5*( ( u(i,j,1)+u(i+1,j,1) )*dzdx(i,j) &
                        +( v(i,j,1)+v(i,j+1,1) )*dzdy(i,j) )*gz(i,j)
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcwsfc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine extrapbcs(s)
      use input
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s

      integer :: i,j

      ! cm1r18 extrapolation formulation:
      ! assumes zh(0) is same as zf(1), and zh(nk+1) is same as zf(nk+1)

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=jb,je
      do i=ib,ie
        s(i,j,0)    = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        s(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
      enddo
      enddo

      end subroutine extrapbcs


#ifdef MPI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcs2(s)
      use input
      implicit none

      real s(ib:ie,jb:je,kb:ke)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1-ngxy,0
            s(i,j,k)=s(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1-ngxy,0
            s(i,j,k)=s(1,j,k)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=ni+1,ni+ngxy
            s(i,j,k)=s(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=ni+1,ni+ngxy
            s(i,j,k)=s(ni,j,k)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1-ngxy,0
            s(i,j,k)=s(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1-ngxy,0
            s(i,j,k)=s(i,1,k)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=nj+1,nj+ngxy
            s(i,j,k)=s(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=nj+1,nj+ngxy
            s(i,j,k)=s(i,nj,k)
          enddo
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcs2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bct2(t)
      use input
      implicit none

      real t(ib:ie,jb:je,kb:ke+1)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=1-ngxy,0
            t(i,j,k)=t(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=1-ngxy,0
            t(i,j,k)=t(1,j,k)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=ni+1,ni+ngxy
            t(i,j,k)=t(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=ni+1,ni+ngxy
            t(i,j,k)=t(ni,j,k)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=1-ngxy,0
            t(i,j,k)=t(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=1-ngxy,0
            t(i,j,k)=t(i,1,k)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=nj+1,nj+ngxy
            t(i,j,k)=t(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=nj+1,nj+ngxy
            t(i,j,k)=t(i,nj,k)
          enddo
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bct2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcu2(u)
      use input
      implicit none

      real u(ib:ie+1,jb:je,kb:ke)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1-ngxy,0
            u(i,j,k)=u(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1-ngxy,0
            u(i,j,k)=u(1,j,k)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=ni+2,ni+1+ngxy
            u(i,j,k)=u(ni+1,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=ni+2,ni+1+ngxy
            u(i,j,k)=u(ni+1,j,k)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1-ngxy,0
            u(i,j,k)=u(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+2
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1-ngxy,0
            u(i,j,k)=u(i,1,k)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=nj+1,nj+ngxy
            u(i,j,k)=u(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+2
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=nj+1,nj+ngxy
            u(i,j,k)=u(i,nj,k)
          enddo
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcu2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcv2(v)
      use input
      implicit none

      real v(ib:ie,jb:je+1,kb:ke)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1-ngxy,0
            v(i,j,k)=v(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+2
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1-ngxy,0
            v(i,j,k)=v(1,j,k)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=ni+1,ni+ngxy
            v(i,j,k)=v(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+2
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=ni+1,ni+ngxy
            v(i,j,k)=v(ni,j,k)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1-ngxy,0
            v(i,j,k)=v(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1-ngxy,0
            v(i,j,k)=v(i,1,k)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=nj+2,nj+1+ngxy
            v(i,j,k)=v(i,nj+1,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=nj+2,nj+1+ngxy
            v(i,j,k)=v(i,nj+1,k)
          enddo
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcv2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcw2(w)
      use input
      implicit none

      real w(ib:ie,jb:je,kb:ke+1)

      integer i,j,k

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=1-ngxy,0
            w(i,j,k)=w(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=1-ngxy,0
            w(i,j,k)=w(1,j,k)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=ni+1,ni+ngxy
            w(i,j,k)=w(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk+1
          do i=ni+1,ni+ngxy
            w(i,j,k)=w(ni,j,k)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=1-ngxy,0
            w(i,j,k)=w(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=1-ngxy,0
            w(i,j,k)=w(i,1,k)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=nj+1,nj+ngxy
            w(i,j,k)=w(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk+1
          do j=nj+1,nj+ngxy
            w(i,j,k)=w(i,nj,k)
          enddo
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcw2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcs2_2d(s)
      use input
      implicit none

      real s(ib:ie,jb:je)

      integer i,j

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          do j=-2,0
          do i=1-ngxy,0
          s(i,j)=s(1,j)
          enddo
          enddo
        endif

        if(patchnww)then
          do j=nj+1,nj+3
          do i=1-ngxy,0
          s(i,j)=s(1,j)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          do j=-2,0
          do i=ni+1,ni+ngxy
          s(i,j)=s(ni,j)
          enddo
          enddo
        endif

        if(patchnee)then
          do j=nj+1,nj+3
          do i=ni+1,ni+ngxy
          s(i,j)=s(ni,j)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          do j=1-ngxy,0
          do i=-2,0
          s(i,j)=s(i,1)
          enddo
          enddo
        endif

        if(patchses)then
          do j=1-ngxy,0
          do i=ni+1,ni+3
          s(i,j)=s(i,1)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          do j=nj+1,nj+ngxy
          do i=-2,0
          s(i,j)=s(i,nj)
          enddo
          enddo
        endif

        if(patchnen)then
          do j=nj+1,nj+ngxy
          do i=ni+1,ni+3
          s(i,j)=s(i,nj)
          enddo
          enddo
        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bcs2_2d


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcew(radbcw,radbce,ua)
      use input
      use constants
      implicit none

      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
 
      integer j,k
      real cbcw,cbce
 
      if(ibw.eq.1.and.wbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k,cbcw)
        do k=1,nk
        do j=1,nj
          cbcw=ua(1,j,k)-cstar
          if(cbcw.lt.0.0)then
            radbcw(j,k)=cbcw
          else
            radbcw(j,k)=0.
          endif
        enddo
        enddo
      endif
 
      if(ibe.eq.1.and.ebc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(j,k,cbce)
        do k=1,nk
        do j=1,nj
          cbce=ua(ni+1,j,k)+cstar
          if(cbce.gt.0.0)then
            radbce(j,k)=cbce
          else
            radbce(j,k)=0.
          endif
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine radbcew


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcns(radbcs,radbcn,va)
      use input
      use constants
      implicit none

      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
 
      integer i,k
      real cbcs,cbcn

      if(ibs.eq.1.and.sbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k,cbcs)
        do k=1,nk
        do i=1,ni
          cbcs=va(i,1,k)-cstar
          if(cbcs.lt.0.0)then
            radbcs(i,k)=cbcs
          else
            radbcs(i,k)=0.
          endif
        enddo
        enddo
      endif
 
      if(ibn.eq.1.and.nbc.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,k,cbcn)
        do k=1,nk
        do i=1,ni
          cbcn=va(i,nj+1,k)+cstar
          if(cbcn.gt.0.0)then
            radbcn(i,k)=cbcn
          else
            radbcn(i,k)=0.
          endif
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine radbcns


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcew4(ruf,radbcw,radbce,u1,u2,dt)
      use input
      use constants
      implicit none

      real, dimension(ib:ie+1) :: ruf
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u1,u2
      real dt

      integer j,k
      real umax,avgw,avge,foo,cbcw,cbce

      umax=csmax

      if(ibw.eq.1.and.wbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(j,k,foo,avgw,cbcw)
        do j=1,nj
          avgw=0.
          do k=1,nk
            foo=(u1(3,j,k)-u1(2,j,k))
            cbcw=dx*ruf(2)*(u1(2,j,k)-u2(2,j,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbcw=max(min(cbcw,0.0),-umax)
            avgw=avgw+cbcw
          enddo
          avgw=avgw/float(nk)
          do k=1,nk
            radbcw(j,k)=avgw
          enddo
        enddo
      endif

      if(ibe.eq.1.and.ebc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(j,k,foo,avge,cbce)
        do j=1,nj
          avge=0.
          do k=1,nk
            foo=(u1(ni+1-1,j,k)-u1(ni+1-2,j,k))
            cbce=dx*ruf(ni+1-1)*(u1(ni+1-1,j,k)-u2(ni+1-1,j,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbce=min(max(cbce,0.0),umax)
            avge=avge+cbce
          enddo
          avge=avge/float(nk)
          do k=1,nk
            radbce(j,k)=avge
          enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine radbcew4


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine radbcns4(rvf,radbcs,radbcn,v1,v2,dt)
      use input
      use constants
      implicit none

      real, dimension(jb:je+1) :: rvf
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v1,v2
      real dt

      integer i,k
      real umax,avgs,avgn,foo,cbcs,cbcn

      umax=csmax

      if(ibs.eq.1.and.sbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,k,avgs,foo,cbcs)
        do i=1,ni
          avgs=0.
          do k=1,nk
            foo=(v1(i,3,k)-v1(i,2,k))
            cbcs=dy*rvf(2)*(v1(i,2,k)-v2(i,2,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbcs=max(min(cbcs,0.0),-umax)
            avgs=avgs+cbcs
          enddo
          avgs=avgs/float(nk)
          do k=1,nk
            radbcs(i,k)=avgs
          enddo
        enddo
      endif

      if(ibn.eq.1.and.nbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,k,avgn,foo,cbcn)
        do i=1,ni
          avgn=0.
          do k=1,nk
            foo=(v1(i,nj+1-1,k)-v1(i,nj+1-2,k))
            cbcn=dy*rvf(nj+1-1)*(v1(i,nj+1-1,k)-v2(i,nj+1-1,k))   &
                   /(dt*(sign(1.e-10,foo)+foo))
            cbcn=min(max(cbcn,0.0),umax)
            avgn=avgn+cbcn
          enddo
          avgn=avgn/float(nk)
          do k=1,nk
            radbcn(i,k)=avgn
          enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine radbcns4


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restrict_openbc_we(rvh,rmh,rho0,u3d)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none

      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u3d

      integer i,j,k
      double precision :: fluxout,fluxin,tem,u1,t3
      double precision, dimension(nk) :: temout,temin

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        temout(k) = 0.0d0
        temin(k)  = 0.0d0
      enddo

      if(wbc.eq.2.and.ibw.eq.1)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          temout(k)=temout(k)-min(0.0,rho0(1,j,k)*u3d(i,j,k)*rvh(j)*rmh(1,j,k))
          temin(k) =temin(k) +max(0.0,rho0(1,j,k)*u3d(i,j,k)*rvh(j)*rmh(1,j,k))
        enddo
        enddo
      endif

      if(ebc.eq.2.and.ibe.eq.1)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          temout(k)=temout(k)+max(0.0,rho0(ni,j,k)*u3d(i,j,k)*rvh(j)*rmh(ni,j,k))
          temin(k) =temin(k) -min(0.0,rho0(ni,j,k)*u3d(i,j,k)*rvh(j)*rmh(ni,j,k))
        enddo
        enddo
      endif

      fluxout = 0.0d0
      fluxin  = 0.0d0

      do k=1,nk
        fluxout = fluxout + temout(k)
        fluxin  = fluxin  + temin(k)
      enddo

#ifdef MPI
      tem=0.0d0
      call MPI_ALLREDUCE(fluxout,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxout=tem
      tem=0.0d0
      call MPI_ALLREDUCE(fluxin,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxin=tem
#endif

      t3=(fluxin+1.0d-20)/(fluxout+1.0d-20)

      if(wbc.eq.2.and.ibw.eq.1)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k,u1)
        do k=1,nk
        do j=1,nj
          u1=rho0(1,j,k)*u3d(i,j,k)
          if(u1.lt.0.0)then
            u3d(i,j,k)=u1*t3/rho0(1,j,k)
          endif
        enddo
        enddo
      endif

      if(ebc.eq.2.and.ibe.eq.1)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k,u1)
        do k=1,nk
        do j=1,nj
          u1=rho0(ni,j,k)*u3d(i,j,k)
          if(u1.gt.0.0)then
            u3d(i,j,k)=u1*t3/rho0(ni,j,k)
          endif
        enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine restrict_openbc_we


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restrict_openbc_sn(ruh,rmh,rho0,v3d)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none

      real, dimension(ib:ie) :: ruh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d

      integer i,j,k
      double precision :: fluxout,fluxin,tem,u1,t3
      double precision, dimension(nk) :: temout,temin

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        temout(k) = 0.0d0
        temin(k)  = 0.0d0
      enddo

      if(sbc.eq.2.and.ibs.eq.1)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          temout(k)=temout(k)-min(0.0,rho0(i,1,k)*v3d(i,j,k)*ruh(i)*rmh(i,1,k))
          temin(k) =temin(k) +max(0.0,rho0(i,1,k)*v3d(i,j,k)*ruh(i)*rmh(i,1,k))
        enddo
        enddo
      endif

      if(nbc.eq.2.and.ibn.eq.1)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          temout(k)=temout(k)+max(0.0,rho0(i,nj,k)*v3d(i,j,k)*ruh(i)*rmh(i,nj,k))
          temin(k) =temin(k) -min(0.0,rho0(i,nj,k)*v3d(i,j,k)*ruh(i)*rmh(i,nj,k))
        enddo
        enddo
      endif

      fluxout = 0.0d0
      fluxin  = 0.0d0

      do k=1,nk
        fluxout = fluxout + temout(k)
        fluxin  = fluxin  + temin(k)
      enddo

#ifdef MPI
      tem=0.0d0
      call MPI_ALLREDUCE(fluxout,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxout=tem
      tem=0.0d0
      call MPI_ALLREDUCE(fluxin,tem,1,MPI_DOUBLE_PRECISION,MPI_SUM,   &
                         MPI_COMM_WORLD,ierr)
      fluxin=tem
#endif

      t3=(fluxin+1.0d-20)/(fluxout+1.0d-20)

      if(sbc.eq.2.and.ibs.eq.1)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k,u1)
        do k=1,nk
        do i=1,ni
          u1=rho0(i,1,k)*v3d(i,j,k)
          if(u1.lt.0.0)then
            v3d(i,j,k)=u1*t3/rho0(i,1,k)
          endif
        enddo
        enddo
      endif

      if(nbc.eq.2.and.ibn.eq.1)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k,u1)
        do k=1,nk
        do i=1,ni
          u1=rho0(i,nj,k)*v3d(i,j,k)
          if(u1.gt.0.0)then
            v3d(i,j,k)=u1*t3/rho0(i,nj,k)
          endif
        enddo
        enddo
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine restrict_openbc_sn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dts)
      use input
      implicit none

      real, intent(in),    dimension(ib:ie) :: uh
      real, intent(in),    dimension(kb:ke) :: rds,sigma
      real, intent(in),    dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in),    dimension(itb:ite,jtb:jte) :: gz,rgzu
      real, intent(in),    dimension(itb:ite,jtb:jte,ktb:kte) :: gx
      real, intent(in),    dimension(jb:je,kb:ke) :: radbcw
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in) :: dts

      integer :: i,j,k
      real :: r1,r2

          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(j,k)
            do k=1,nk
            do j=1,nj
              u3d(1,j,k)=u3d(1,j,k)+dts*( -radbcw(j,k)       &
                        *(u3d(2,j,k)-u3d(1,j,k))*rdx*uh(1)   &
                           +uten(1,j,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores u at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=1,nj
              ! lowest model level:
              do i=1,2
                dum1(i,j,1) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=1,2
                dum1(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=1,2
                dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
              u3d(1,j,k)=u3d(1,j,k)+dts*( -radbcw(j,k)*(                            &
                     gz(1,j)*(u3d(2,j,k)*rgzu(2,j)-u3d(1,j,k)*rgzu(1,j))*rdx*uh(1)  &
                    +0.5*( gx(1,j,k+1)*(dum1(2,j,k+1)+dum1(1,j,k+1))                &
                          -gx(1,j,k  )*(dum1(2,j,k  )+dum1(1,j,k  )) )*rdsf(k)      &
                                                       )+uten(1,j,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine ssopenbcw


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dts)
      use input
      implicit none

      real, intent(in),    dimension(ib:ie) :: uh
      real, intent(in),    dimension(kb:ke) :: rds,sigma
      real, intent(in),    dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in),    dimension(itb:ite,jtb:jte) :: gz,rgzu
      real, intent(in),    dimension(itb:ite,jtb:jte,ktb:kte) :: gx
      real, intent(in),    dimension(jb:je,kb:ke) :: radbce
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in) :: dts

      integer :: i,j,k
      real :: r1,r2

          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(j,k)
            do k=1,nk
            do j=1,nj
              u3d(ni+1,j,k)=u3d(ni+1,j,k)+dts*( -radbce(j,k)           &
                           *(u3d(ni+1,j,k)-u3d(ni  ,j,k))*rdx*uh(ni)   &
                           +uten(ni+1,j,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores u at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=1,nj
              ! lowest model level:
              do i=ni,ni+1
                dum1(i,j,1) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=ni,ni+1
                dum1(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=ni,ni+1
                dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
              u3d(ni+1,j,k)=u3d(ni+1,j,k)+dts*( -radbce(j,k)*(                                &
                     gz(ni,j)*(u3d(ni+1,j,k)*rgzu(ni+1,j)-u3d(ni,j,k)*rgzu(ni,j))*rdx*uh(ni)  &
                    +0.5*( gx(ni,j,k+1)*(dum1(ni+1,j,k+1)+dum1(ni,j,k+1))                     &
                          -gx(ni,j,k  )*(dum1(ni+1,j,k  )+dum1(ni,j,k  )) )*rdsf(k)           &
                                                             )+uten(ni+1,j,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine ssopenbce


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dts)
      use input
      implicit none

      real, intent(in),    dimension(jb:je) :: vh
      real, intent(in),    dimension(kb:ke) :: rds,sigma
      real, intent(in),    dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in),    dimension(itb:ite,jtb:jte) :: gz,rgzv
      real, intent(in),    dimension(itb:ite,jtb:jte,ktb:kte) :: gy
      real, intent(in),    dimension(ib:ie,kb:ke) :: radbcs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in) :: dts

      integer :: i,j,k
      real :: r1,r2

          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(i,k)
            do k=1,nk
            do i=1,ni
              v3d(i,1,k)=v3d(i,1,k)+dts*( -radbcs(i,k)      &
                        *(v3d(i,2,k)-v3d(i,1,k))*rdy*vh(1)  &
                        +vten(i,1,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores v at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=1,2
              ! lowest model level:
              do i=1,ni
                dum1(i,j,1) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=1,ni
                dum1(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=1,ni
                dum1(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do i=1,ni
              v3d(i,1,k)=v3d(i,1,k)+dts*( -radbcs(i,k)*(                            &
                     gz(i,1)*(v3d(i,2,k)*rgzv(i,2)-v3d(i,1,k)*rgzv(i,1))*rdy*vh(1)  &
                    +0.5*( gy(i,1,k+1)*(dum1(i,2,k+1)+dum1(i,1,k+1))                &
                          -gy(i,1,k  )*(dum1(i,2,k  )+dum1(i,1,k  )) )*rdsf(k)      &
                                                       )+vten(i,1,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine ssopenbcs


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dts)
      use input
      implicit none

      real, intent(in),    dimension(jb:je) :: vh
      real, intent(in),    dimension(kb:ke) :: rds,sigma
      real, intent(in),    dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in),    dimension(itb:ite,jtb:jte) :: gz,rgzv
      real, intent(in),    dimension(itb:ite,jtb:jte,ktb:kte) :: gy
      real, intent(in),    dimension(ib:ie,kb:ke) :: radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in) :: dts

      integer :: i,j,k
      real :: r1,r2

          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(i,k)
            do k=1,nk
            do i=1,ni
              v3d(i,nj+1,k)=v3d(i,nj+1,k)+dts*( -radbcn(i,k)        &
                           *(v3d(i,nj+1,k)-v3d(i,nj,k))*rdy*vh(nj)  &
                           +vten(i,nj+1,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores v at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=nj,nj+1
              ! lowest model level:
              do i=1,ni
                dum1(i,j,1) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=1,ni
                dum1(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=1,ni
                dum1(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do i=1,ni
              v3d(i,nj+1,k)=v3d(i,nj+1,k)+dts*( -radbcn(i,k)*(                                &
                     gz(i,nj)*(v3d(i,nj+1,k)*rgzv(i,nj+1)-v3d(i,nj,k)*rgzv(i,nj))*rdy*vh(nj)  &
                    +0.5*( gy(i,nj,k+1)*(dum1(i,nj+1,k+1)+dum1(i,nj,k+1))                     &
                          -gy(i,nj,k  )*(dum1(i,nj+1,k  )+dum1(i,nj,k  )) )*rdsf(k)           &
                                                             )+vten(i,nj+1,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine ssopenbcn


  END MODULE bc_module
