  MODULE bc_module

  implicit none
  private

  !GPU enabled
  public :: bc2d, bcu, bcv, bcw, bcs, bcp, bcwsfc
  public :: radbcew, radbcns, radbcew4, radbcns4, extrapbcs
  public :: restrict_openbc_sn, restrict_openbc_we
  public :: ssopenbcw, ssopenbce, ssopenbcs, ssopenbcn, bcs_tend_halo

  !MPI + GPU enabled
  public :: bcs2, bct2, bcu2, bcv2, bcw2, bcs2_2d

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine bcp(p,device)
      use input, only : ib,ie,jb,je,kb,ke,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,ni,nj,nk,timestats,time_bc,mytime
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: p
      logical, optional, intent(in) :: device

      integer i,j,k
      logical :: openacc 

      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared) private(j,k)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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

      subroutine bcwsfc(gz,dzdx,dzdy,u,v,w,device)
      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ni,nj, &
          timestats,time_bc,mytime

      implicit none

      real, intent(in), dimension(itb:ite,jtb:jte) :: gz
      real, intent(in), dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w
      logical, optional, intent(in) :: device


      integer :: i,j
      logical :: openacc

      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------

      !$omp parallel do default(shared) private(i,j)
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine bcs_tend_halo(s,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,wbc,ebc,sbc,nbc, &
          time_bcst,timestats,mytime
      implicit none

      double precision, dimension(ib:ie,jb:je,kb:ke) :: s
      logical, optional, intent(in) :: device

      integer i,j,k
      logical :: openacc

      ! For periodic lateral boundary conditions only:  
      ! account for tendencies from halo zones. 

#ifndef MPI

      call SetGPUParams(openacc,device)

!-----------------------------------------------------------------------
!  west boundary condition

      if(wbc.eq.1)then
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do k=1,nk
        do j=1,nj
          s(1,j,k)=s(1,j,k)+s(ni+1,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  east boundary condition

      if(ebc.eq.1)then
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do k=1,nk
        do j=1,nj
          s(ni,j,k)=s(ni,j,k)+s(0,j,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  south boundary condition

      if(sbc.eq.1)then
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do k=1,nk
        do i=1,ni
          s(i,1,k)=s(i,1,k)+s(i,nj+1,k)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  north boundary condition

      if(nbc.eq.1)then
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do k=1,nk
        do i=1,ni
          s(i,nj,k)=s(i,nj,k)+s(i,0,k)
        enddo
        enddo
      endif

#endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bcst=time_bcst+mytime()

      end subroutine bcs_tend_halo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine extrapbcs(s,device)
      use input, only : ib,ie,jb,je,kb,ke,nk,cgs1,cgs2,cgs3,cgt1,cgt2,cgt3
      implicit none

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s
      logical, optional, intent(in) :: device

      integer :: i,j
      logical :: openacc

      ! cm1r18 extrapolation formulation:
      ! assumes zh(0) is same as zf(1), and zh(nk+1) is same as zf(nk+1)

      call SetGPUParams(openacc,device)

      !$omp parallel do default(shared) private(i,j)
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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

      subroutine bct2(t,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,ngxy, &
          ibw,ibe,ibs,ibn,patchsww,patchnww,patchsee,patchnee, &
          patchsws,patchses,patchnwn,patchnen,timestats,time_bc,mytime
      implicit none

      real, intent(inout) ::  t(ib:ie,jb:je,kb:ke+1)
      logical, optional, intent(in) :: device

      integer i,j,k
      logical :: openacc

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------
      call SetGPUParams(openacc,device)

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk+1
          do i=1-ngxy,0
            t(i,j,k)=t(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
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
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk+1
          do i=ni+1,ni+ngxy
            t(i,j,k)=t(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk+1
          do j=1-ngxy,0
            t(i,j,k)=t(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk+1
          do j=nj+1,nj+ngxy
            t(i,j,k)=t(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
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


      subroutine bcu2(u,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,ngxy, &
          ibw,ibe,ibs,ibn,patchsww,patchnww,patchsee,patchnee, &
          patchsws,patchses,patchnwn,patchnen,timestats,time_bc,mytime
      implicit none

      real u(ib:ie+1,jb:je,kb:ke)
      logical, optional, intent(in) :: device

      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)
!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do i=1-ngxy,0
            u(i,j,k)=u(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do i=ni+2,ni+1+ngxy
            u(i,j,k)=u(ni+1,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do j=1-ngxy,0
            u(i,j,k)=u(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+2
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do j=nj+1,nj+ngxy
            u(i,j,k)=u(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+2
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine bcv2(v,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,ngxy, &
          ibw,ibe,ibs,ibn,patchsww,patchnww,patchsee,patchnee, &
          patchsws,patchses,patchnwn,patchnen,timestats,time_bc,mytime
      implicit none

      real v(ib:ie,jb:je+1,kb:ke)
      logical, optional, intent(in) :: device

      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)
!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do i=1-ngxy,0
            v(i,j,k)=v(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+2
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do i=ni+1,ni+ngxy
            v(i,j,k)=v(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+2
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do j=1-ngxy,0
            v(i,j,k)=v(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk
          do j=nj+2,nj+1+ngxy
            v(i,j,k)=v(i,nj+1,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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

      subroutine bcw2(w,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,ngxy, &
          ibw,ibe,ibs,ibn,patchsww,patchnww,patchsee,patchnee, &
          patchsws,patchses,patchnwn,patchnen,timestats,time_bc,mytime
      implicit none

      real w(ib:ie,jb:je,kb:ke+1)
      logical, optional, intent(in) :: device

      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)
!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk+1
          do i=1-ngxy,0
            w(i,j,k)=w(1,j,k)
          enddo
          enddo
        endif

        if(patchnww)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk+1
          do i=ni+1,ni+ngxy
            w(i,j,k)=w(ni,j,k)
          enddo
          enddo
        endif

        if(patchnee)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk+1
          do j=1-ngxy,0
            w(i,j,k)=w(i,1,k)
          enddo
          enddo
        endif

        if(patchses)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do k=1,nk+1
          do j=nj+1,nj+ngxy
            w(i,j,k)=w(i,nj,k)
          enddo
          enddo
        endif

        if(patchnen)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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

      subroutine bcs2_2d(s,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk,ngxy, &
          ibw,ibe,ibs,ibn,patchsww,patchnww,patchsee,patchnee, &
          patchsws,patchses,patchnwn,patchnen,timestats,time_bc,mytime
      implicit none

      real s(ib:ie,jb:je)
      logical, optional, intent(in) :: device

      integer i,j
      logical :: openacc

      !Setup GPU specific parameters
      call SetGPUParams(openacc,device)

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=-2,0
          do i=1-ngxy,0
          s(i,j)=s(1,j)
          enddo
          enddo
        endif

        if(patchnww)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=nj+1,nj+3
          do i=1-ngxy,0
          s(i,j)=s(1,j)
          enddo
          enddo
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=-2,0
          do i=ni+1,ni+ngxy
          s(i,j)=s(ni,j)
          enddo
          enddo
        endif

        if(patchnee)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=nj+1,nj+3
          do i=ni+1,ni+ngxy
          s(i,j)=s(ni,j)
          enddo
          enddo
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=1-ngxy,0
          do i=-2,0
          s(i,j)=s(i,1)
          enddo
          enddo
        endif

        if(patchses)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=1-ngxy,0
          do i=ni+1,ni+3
          s(i,j)=s(i,1)
          enddo
          enddo
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
          do j=nj+1,nj+ngxy
          do i=-2,0
          s(i,j)=s(i,nj)
          enddo
          enddo
        endif

        if(patchnen)then
          !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine radbcew(radbcw,radbce,ua,device)
      use input, only : ib,ie,jb,je,kb,ke,wbc,ebc,sbc,nbc,ni,nj,nk, &
          ibw,ibe,ibs,ibn,timestats,time_bc,mytime
      use constants
      implicit none

      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      logical, optional, intent(in) :: device
 
      integer j,k
      real cbcw,cbce
      logical :: openacc

 
      call SetGPUParams(openacc,device)
      if(ibw.eq.1.and.wbc.eq.2)then
        !$omp parallel do default(shared) private(j,k,cbcw)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(j,k,cbce)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine radbcns(radbcs,radbcn,va,device)
      use input, only : ib,ie,jb,je,kb,ke,wbc,ebc,sbc,nbc,ni,nj,nk, &
          ibw,ibe,ibs,ibn,timestats,time_bc,mytime

      use constants
      implicit none

      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      logical, optional, intent(in) :: device
 
      integer i,k
      real cbcs,cbcn
      logical :: openacc


      call SetGPUParams(openacc,device)
      if(ibs.eq.1.and.sbc.eq.2)then
        !$omp parallel do default(shared) private(i,k,cbcs)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc) 
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
        !$omp parallel do default(shared) private(i,k,cbcn)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine radbcew4(ruf,radbcw,radbce,u1,u2,dt,device)
      use input, only : ib,ie,jb,je,kb,ke,wbc,ebc,sbc,nbc,ni,nj,nk, &
          ibw,ibe,ibs,ibn,dx,timestats,time_bc,mytime

      use constants
      implicit none

      real, dimension(ib:ie+1) :: ruf
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u1,u2
      real dt
      logical, optional, intent(in) :: device

      integer j,k
      real umax,avgw,avge,foo,cbcw,cbce
      logical :: openacc

      call SetGPUParams(openacc,device)
      umax=csmax

      if(ibw.eq.1.and.wbc.eq.2)then
       !$omp parallel do default(shared) private(j,k,foo,avgw,cbcw)
       !$acc parallel loop gang vector default(present) private(j,k,foo,avgw,cbcw) if(openacc)
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
        !$omp parallel do default(shared) private(j,k,foo,avge,cbce)
        !$acc parallel loop gang vector default(present) private(j,k,foo,avge,cbce) if(openacc)
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


      subroutine radbcns4(rvf,radbcs,radbcn,v1,v2,dt,device)
      use input, only : ib,ie,jb,je,kb,ke,wbc,ebc,sbc,nbc,ni,nj,nk, &
          ibw,ibe,ibs,ibn,dy,timestats,time_bc,mytime
      use constants
      implicit none

      real, dimension(jb:je+1) :: rvf
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v1,v2
      real dt
      logical, optional, intent(in) :: device

      integer i,k
      real umax,avgs,avgn,foo,cbcs,cbcn
      logical :: openacc


      call SetGPUParams(openacc,device)
      umax=csmax

      if(ibs.eq.1.and.sbc.eq.2)then
        !$omp parallel do default(shared) private(i,k,avgs,foo,cbcs)
        !$acc parallel loop gang vector default(present) private(i,k,avgs,foo,cbcs) if(openacc)
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
        !$omp parallel do default(shared) private(i,k,avgn,foo,cbcn)
        !$acc parallel loop gang vector default(present) private(i,k,avgn,foo,cbcn) if(openacc)
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


      subroutine restrict_openbc_we(rvh,rmh,rho0,u3d,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk, &
          wbc,ebc,sbc,nbc,ibw,ibe,ibs,ibn,time_bc, &
          ierr,timestats,mytime
#ifdef MPI
      use mpi
#endif
      implicit none

      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      logical, optional, intent(in) :: device

      integer i,j,k
      double precision :: fluxout,fluxin,tem,u1,t3
      double precision :: temp1,temp2
      double precision, dimension(nk) :: temout,temin
      logical :: openacc

      call SetGPUParams(openacc,device)
      !$acc data create(temout,temin)

      !$omp parallel do default(shared) private(k)
      !$acc parallel loop gang vector default(present) if(openacc)
      do k=1,nk
        temout(k) = 0.0d0
        temin(k)  = 0.0d0
      enddo

      if(wbc.eq.2.and.ibw.eq.1)then
        i=1
        !$omp parallel do default(shared) private(j,k,temp1,temp2)
        !$acc parallel default(present) reduction(+:temp1,temp2) if(openacc)
        !$acc loop gang
        do k=1,nk
          temp1 = temout(k)
          temp2 = temin(k)
          !$acc loop vector reduction(+:temp1,temp2) 
          do j=1,nj
            temp1 = temp1 - min(0.0,rho0(1,j,k)*u3d(i,j,k)*rvh(j)*rmh(1,j,k))
            temp2 = temp2 + max(0.0,rho0(1,j,k)*u3d(i,j,k)*rvh(j)*rmh(1,j,k))
          enddo
          temout(k) = temp1
          temin(k)  = temp2
        enddo
        !$acc end parallel
      endif

      if(ebc.eq.2.and.ibe.eq.1)then
        i=ni+1
        !$omp parallel do default(shared) private(j,k)
        !$acc parallel default(present) reduction(+:temp1,temp2) if(openacc)
        !$acc loop gang
        do k=1,nk
          temp1 = temout(k)
          temp2 = temin(k)
          !$acc loop vector reduction(+:temp1,temp2)
          do j=1,nj
            temp1 = temp1 + max(0.0,rho0(ni,j,k)*u3d(i,j,k)*rvh(j)*rmh(ni,j,k))
            temp2 = temp2 - min(0.0,rho0(ni,j,k)*u3d(i,j,k)*rvh(j)*rmh(ni,j,k))
          enddo
          temout(k) = temp1
          temin(k)  = temp2
        enddo
        !$acc end parallel
      endif

      fluxout = 0.0d0
      fluxin  = 0.0d0

      !$acc parallel loop gang vector default(present) reduction(+:fluxout,fluxin) if(openacc)
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
        !$omp parallel do default(shared) private(j,k,u1)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(j,k,u1)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do k=1,nk
        do j=1,nj
          u1=rho0(ni,j,k)*u3d(i,j,k)
          if(u1.gt.0.0)then
            u3d(i,j,k)=u1*t3/rho0(ni,j,k)
          endif
        enddo
        enddo
      endif

      !$acc end data

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine restrict_openbc_we


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restrict_openbc_sn(ruh,rmh,rho0,v3d,device)
      use input, only : ib,ie,jb,je,kb,ke,ni,nj,nk, &
          wbc,ebc,sbc,nbc,ibw,ibe,ibs,ibn,time_bc, &
          ierr,timestats,mytime
#ifdef MPI
      use mpi
#endif
      implicit none

      real, dimension(ib:ie) :: ruh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      logical, optional, intent(in) :: device

      integer i,j,k
      double precision :: fluxout,fluxin,tem,u1,t3,tmp1,tmp2
      double precision, dimension(nk) :: temout,temin
      logical :: openacc

      call SetGPUParams(openacc,device)
      if(sbc.eq.2.and.ibs.eq.1)then
        j=1
!$omp parallel do default(shared) private(i,k,tmp1,tmp2)
        !$acc parallel default(present) reduction(+:tmp1,tmp2) if(openacc)
        !$acc loop gang
        do k=1,nk
        tmp1 = 0.D0
        tmp2 = 0.D0
        !$acc loop vector reduction (+:tmp1,tmp2)
        do i=1,ni
          tmp1=tmp1-min(0.0,rho0(i,1,k)*v3d(i,j,k)*ruh(i)*rmh(i,1,k))
          tmp2=tmp2+max(0.0,rho0(i,1,k)*v3d(i,j,k)*ruh(i)*rmh(i,1,k))
        enddo
        temout(k)=tmp1
        temin(k)=tmp2
        enddo
        !$acc end parallel
      endif

      if(nbc.eq.2.and.ibn.eq.1)then
        j=nj+1
!$omp parallel do default(shared) private(i,k,tmp1,tmp2)
        !$acc parallel default(present) reduction(+:tmp1,tmp2) if(openacc)
        !$acc loop gang
        do k=1,nk
        tmp1 = 0.D0
        tmp2 = 0.D0
        !$acc loop vector reduction (+:tmp1,tmp2)
        do i=1,ni
          tmp1=tmp1+max(0.0,rho0(i,nj,k)*v3d(i,j,k)*ruh(i)*rmh(i,nj,k))
          tmp2=tmp2-min(0.0,rho0(i,nj,k)*v3d(i,j,k)*ruh(i)*rmh(i,nj,k))
        enddo
        temout(k)=tmp1
        temin(k)=tmp2
        enddo
        !$acc end parallel
      endif

      fluxout = 0.0d0
      fluxin  = 0.0d0

      !$acc parallel loop gang vector default(present) reduction(+:fluxout,fluxin) if(openacc)
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
        !$omp parallel do default(shared) private(i,k,u1)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
!$omp parallel do default(shared) private(i,k,u1)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine ssopenbcw(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbcw,dum1,u3d,uten,dts,device)
      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ktb,kte,ni,nj,nk, &
          cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,rdx,terrain_flag,timestats,time_bc,mytime
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
      logical, optional, intent(in) :: device

      integer :: i,j,k
      real :: r1,r2
      logical :: openacc

      call SetGPUParams(openacc,device)

          IF(.not.terrain_flag)THEN
            ! no terrain:
            !$omp parallel do default(shared) private(j,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
            !$omp parallel do default(shared) private(i,j,k,r1,r2)
            !JMD what to do about nested directives
            !$acc parallel default(present)
            !$acc loop gang
            do j=1,nj
              ! lowest model level:
              !$acc loop vector
              do i=1,2
                dum1(i,j,1) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
              enddo
              ! upper-most model level:
              !$acc loop vector
              do i=1,2
                dum1(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
              enddo
              ! interior:
              !$acc loop vector
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              !$acc loop seq
              do i=1,2
                dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
              enddo
              enddo
            enddo
            !$acc end parallel

            !$omp parallel do default(shared) private(i,j,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine ssopenbce(uh,rds,sigma,rdsf,sigmaf,gz,rgzu,gx,radbce,dum1,u3d,uten,dts,device)
      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ktb,kte,ni,nj,nk, &
          cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,rdx,terrain_flag,timestats,time_bc,mytime
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
      logical, optional, intent(in) :: device

      integer :: i,j,k
      real :: r1,r2
      logical :: openacc

      call SetGPUParams(openacc,device)
          IF(.not.terrain_flag)THEN
            ! no terrain:
            !$omp parallel do default(shared) private(j,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
            !$omp parallel do default(shared) private(i,j,k,r1,r2)
            !JMD nested-openacc directives
            !$acc parallel default(present)
            !$acc loop gang
            do j=1,nj
              ! lowest model level:
              !$acc loop vector
              do i=ni,ni+1
                dum1(i,j,1) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
              enddo
              ! upper-most model level:
              !$acc loop vector
              do i=ni,ni+1
                dum1(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
              enddo
              ! interior:
              !$acc loop vector
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              !$acc loop seq
              do i=ni,ni+1
                dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
              enddo
              enddo
            enddo
            !$acc end parallel

            !$omp parallel do default(shared) private(i,j,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine ssopenbcs(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcs,dum1,v3d,vten,dts,device)
      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ktb,kte,ni,nj,nk, &
          cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,rdy,terrain_flag,timestats,time_bc,mytime
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
      logical, optional, intent(in) :: device

      integer :: i,j,k
      real :: r1,r2
      logical :: openacc


      call SetGPUParams(openacc,device)

          IF(.not.terrain_flag)THEN
            ! no terrain:
            !$omp parallel do default(shared) private(i,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
            !$omp parallel do default(shared) private(i,j,k,r1,r2)
            !$acc parallel default(present) if(openacc)
            !$acc loop gang 
            do j=1,2
              ! lowest model level:
              !$acc loop vector
              do i=1,ni
                dum1(i,j,1) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
              enddo
              ! upper-most model level:
              !$acc loop vector
              do i=1,ni
                dum1(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
              enddo
              ! interior:
              !$acc loop vector
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              !$acc loop seq 
              do i=1,ni
                dum1(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
              enddo
              enddo
            enddo
            !$acc end parallel

            !$omp parallel do default(shared) private(i,j,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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


      subroutine ssopenbcn(vh,rds,sigma,rdsf,sigmaf,gz,rgzv,gy,radbcn,dum1,v3d,vten,dts,device)
      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ktb,kte,ni,nj,nk, &
          cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,rdy,terrain_flag,timestats,time_bc,mytime
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
      logical, optional, intent(in) :: device

      integer :: i,j,k
      real :: r1,r2
      logical :: openacc

      call SetGPUParams(openacc,device)
          IF(.not.terrain_flag)THEN
            ! no terrain:
            !$omp parallel do default(shared) private(i,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
            !$omp parallel do default(shared) private(i,j,k,r1,r2)
            !JMD nested openacc 
            !$acc parallel default(present) if(openacc)
            !$acc loop gang
            do j=nj,nj+1
              ! lowest model level:
              !$acc loop vector
              do i=1,ni
                dum1(i,j,1) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
              enddo
              ! upper-most model level:
              !$acc loop vector
              do i=1,ni
                dum1(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
              enddo
              ! interior:
              !$acc loop vector
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              !$acc loop seq
              do i=1,ni
                dum1(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
              enddo
              enddo
            enddo
            !$acc end parallel

            !$omp parallel do default(shared) private(i,j,k)
            !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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

!GPU================================================================================================
!GPU ported subroutine below here 
!GPU================================================================================================

      subroutine bc2d(s,device)
      use input, only : ib,ie,jb,je,ngxy,ni,nj,nk,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,timestats,time_bc,mytime 
      implicit none

      real, dimension(ib:ie,jb:je) :: s
      logical, intent(in), optional :: device

      integer i,j
      logical :: openacc


      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do j=jb,je
        do i=1,ngxy
          s(1-i,j)=s(ni+1-i,j)
        enddo
        enddo
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do j=jb,je
        do i=1,ngxy
          s(ni+i,j)=s(i,j)
        enddo
        enddo
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do j=1,ngxy
        do i=ib,ie
          s(i,1-j)=s(i,nj+1-j)
        enddo
        enddo
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do j=1,ngxy
        do i=ib,ie
          s(i,nj+j)=s(i,j)
        enddo
        enddo
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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
        !$omp parallel do default(shared) private(i,j)
        !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
        do j=1,ngxy
        do i=ib,ie
          s(i,nj+j)=s(i,nj+1-j)
        enddo
        enddo
      endif

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()

      end subroutine bc2d

      subroutine bcs(s,device)
      use input, only : ib,ie,jb,je,kb,ke,ngxy,ni,nj,nk,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,timestats,time_bcs,mytime
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke) :: s
      logical, optional :: device

      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
      if(wbc.eq.1)then
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
        do k=1,nk
           do j=jb,je
              do i=1,ngxy
                 s(1-i,j,k)=s(ni+1-i,j,k)
              end do
           end do
        end do
      elseif(wbc.eq.2)then
#else
      if(ibw.eq.1.and.wbc.eq.2)then
#endif
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
        do k=1,nk
           do j=0,nj+1
              do i=1-ngxy,0
                 s(i,j,k)=s(1,j,k)
              end do
           end do
        end do
#ifdef MPI
      elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
      elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=0,nj+1
             do i=1,ngxy
                s(1-i,j,k)=s(i,j,k)
             enddo
          enddo
       end do
      endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
      if(ebc.eq.1)then
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
        do k=1,nk
           do j=jb,je
              do i=1,ngxy
                 s(ni+i,j,k)=s(i,j,k)
              end do
           end do
        end do
      elseif(ebc.eq.2)then
#else
      if(ibe.eq.1.and.ebc.eq.2)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=0,nj+1
             do i=ni+1,ni+ngxy
                s(i,j,k)=s(ni,j,k)
             end do
          end do
       end do
#ifdef MPI
      elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
      elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=0,nj+1
             do i=1,ngxy
                s(ni+i,j,k)=s(ni+1-i,j,k)
             end do
          end do
       end do 
      endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
      if(sbc.eq.1)then
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=1,ngxy
             do i=ib,ie
                s(i,1-j,k)=s(i,nj+1-j,k)
             end do
          end do
       end do 
      elseif(sbc.eq.2)then
#else
      if(ibs.eq.1.and.sbc.eq.2)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=1-ngxy,0
             do i=0,ni+1
                s(i,j,k)=s(i,1,k)
             end do
          end do
       end do 
#ifdef MPI
      elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
      elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=1,ngxy
             do i=0,ni+1
                s(i,1-j,k)=s(i,j,k)
             end do
          end do
       end do 
      endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
      if(nbc.eq.1)then
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=1,ngxy
             do i=ib,ie
                s(i,nj+j,k)=s(i,j,k)
             end do
          end do
       end do 
      elseif(nbc.eq.2)then
#else
      if(ibn.eq.1.and.nbc.eq.2)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=nj+1,nj+ngxy
             do i=0,ni+1
                s(i,j,k)=s(i,nj,k)
             end do
          end do
       end do
#ifdef MPI
      elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
      elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector default(present) collapse(3) if(openacc)
       do k=1,nk
          do j=1,ngxy
             do i=0,ni+1
                s(i,nj+j,k)=s(i,nj+1-j,k)
             end do
          end do
       end do 
      endif

!-----------------------------------------------------------------------

      !!$acc end data

      if(timestats.ge.1) time_bcs=time_bcs+mytime()

      end subroutine bcs

      subroutine bcu(u,device)
      use input, only : ib,ie,jb,je,kb,ke,ngxy,ni,nj,nk,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,timestats,time_bc,mytime
      implicit none

      real, dimension(ib:ie+1,jb:je,kb:ke) :: u
      logical, optional, intent(in) :: device

      integer i,j,k  
      logical :: openacc

      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
    if(wbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=jb,je
            do i=1,ngxy
               u(1-i,j,k)=u(ni+1-i,j,k)
            end do
         end do
      end do
    elseif(wbc.eq.2)then
#else
    if(ibw.eq.1.and.wbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=0,nj+1
            do i=1-ngxy,0
               u(i,j,k)=u(1,j,k)
            end do
         end do
      end do
#ifdef MPI
    elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
    elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
      do k=1,nk
         do j=0,nj+1
            u(1,j,k)=0.
            !$acc loop seq
            do i=1,ngxy
               u(1-i,j,k)=-u(1+i,j,k)
            end do
         end do
      end do
    endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
    if(ebc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=jb,je
            do i=1,ngxy
               u(ni+1+i,j,k)=u(1+i,j,k)
            end do
         end do
      end do
    elseif(ebc.eq.2)then
#else
    if(ibe.eq.1.and.ebc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=0,nj+1
            do i=ni+2,ni+1+ngxy
               u(i,j,k)=u(ni+1,j,k)
            end do
         end do
      end do
#ifdef MPI
    elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
    elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
      do k=1,nk
         do j=0,nj+1
            u(ni+1,j,k)=0.0
            !$acc loop seq
            do i=1,ngxy
               u(ni+1+i,j,k)=-u(ni+1-i,j,k)
            end do
        end do
      end do
    endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
    if(sbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=1,ngxy
            do i=ib,ie+1
               u(i,1-j,k)=u(i,nj+1-j,k)
            end do
         end do
      end do
    elseif(sbc.eq.2)then
#else
    if(ibs.eq.1.and.sbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=1-ngxy,0
            do i=0,ni+2
               u(i,j,k)=u(i,1,k)
            end do
         end do
      end do
#ifdef MPI
    elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
    elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=1,ngxy
            do i=0,ni+2
               u(i,1-j,k)=u(i,j,k)
            end do
         end do
      end do
    endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
    if(nbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=1,ngxy
            do i=ib,ie+1
               u(i,nj+j,k)=u(i,j,k)
            end do
        end do
      end do
    elseif(nbc.eq.2)then
#else
    if(ibn.eq.1.and.nbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=nj+1,nj+ngxy
            do i=0,ni+2
               u(i,j,k)=u(i,nj,k)
            end do
         end do
      end do
#ifdef MPI
    elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
    elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      do k=1,nk
         do j=1,ngxy
            do i=0,ni+2
               u(i,nj+j,k)=u(i,nj+1-j,k)
            end do
         end do
      end do
    endif


!-----------------------------------------------------------------------

      if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine bcu

      subroutine bcv(v,device)
      use input, only : ib,ie,jb,je,kb,ke,ngxy,ni,nj,nk,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,timestats,time_bc,mytime,axisymm
      implicit none

      real, dimension(ib:ie,jb:je+1,kb:ke) :: v
      logical, optional, intent(in) :: device
 
      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------
!  south boundary condition


#ifndef MPI
    if(sbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=1,ngxy
        do i=ib,ie
          v(i,1-j,k)=v(i,nj+1-j,k)
        enddo
        enddo
      ENDDO
    elseif(sbc.eq.2)then
#else
    if(ibs.eq.1.and.sbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=1-ngxy,0
        do i=0,ni+1
          v(i,j,k)=v(i,1,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
    elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=1,ngxy
        do i=0,ni+1
          v(i,1-j,k)=-v(i,j+1,k)
        enddo
        enddo
      ENDDO
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
      DO k=1,nk
        do i=0,ni+1
          v(i,   1,k)=0.
        enddo
      ENDDO
    endif

!-----------------------------------------------------------------------
!  north boundary condition

    IF(axisymm.eq.1)THEN
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
      DO k=1,nk
        do i=ib,ie
          v(i,2,k)=v(i,1,k)
        enddo
      ENDDO
    ENDIF

#ifndef MPI
    if(nbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=1,ngxy
        do i=ib,ie
          v(i,nj+1+j,k)=v(i,1+j,k)
        enddo
        enddo
      ENDDO
    elseif(nbc.eq.2)then
#else
    if(ibn.eq.1.and.nbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=nj+2,nj+1+ngxy
        do i=0,ni+1
          v(i,j,k)=v(i,nj+1,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
    elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=1,ngxy
        do i=0,ni+1
          v(i,nj+1+j,k)=-v(i,nj+1-j,k)
        enddo
        enddo
      ENDDO
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
      DO k=1,nk
        do i=0,ni+1
          v(i,nj+1,k)=0.
        enddo
      ENDDO
    endif

!-----------------------------------------------------------------------
!  west boundary condition

#ifndef MPI
    if(wbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=jb,je+1
        do i=1,ngxy
          v(1-i,j,k)=v(ni+1-i,j,k)
        enddo
        enddo
      ENDDO
    elseif(wbc.eq.2)then
#else
    if(ibw.eq.1.and.wbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=0,nj+2
        do i=1-ngxy,0
          v(i,j,k)=v(1,j,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
    elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=0,nj+2
        do i=1,ngxy
          v(1-i,j,k)=v(i,j,k)
        enddo
        enddo
      ENDDO
    endif

    IF(axisymm.eq.1)THEN
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=0,nj+2
        do i=1,ngxy
          v(1-i,j,k) = -v(i,j,k)
        enddo
        enddo
      ENDDO
    ENDIF

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
    if(ebc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=jb,je+1
        do i=1,ngxy
          v(ni+i,j,k)=v(i,j,k)
        enddo
        enddo
      ENDDO
    elseif(ebc.eq.2)then
#else
    if(ibe.eq.1.and.ebc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=0,nj+2
        do i=ni+1,ni+ngxy
          v(i,j,k)=v(ni,j,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
    elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=1,nk
        do j=0,nj+2
        do i=1,ngxy
          v(ni+i,j,k)=v(ni+1-i,j,k)
        enddo
        enddo
      ENDDO
    endif


!-----------------------------------------------------------------------

    if(timestats.ge.1) time_bc=time_bc+mytime()
 
      end subroutine bcv

      subroutine bcw(w,flag,device)
      use input, only : ib,ie,jb,je,kb,ke,ngxy,ni,nj,nk,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,timestats,time_bc,mytime
      implicit none

      real, dimension(ib:ie,jb:je,kb:ke+1) :: w
      integer flag
      logical, optional, intent(in) :: device
 
      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)
!-----------------------------------------------------------------------
!  west boundary condition

 
#ifndef MPI
    if(wbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=jb,je
        do i=1,ngxy
          w(1-i,j,k)=w(ni+1-i,j,k)
        enddo
        enddo
      ENDDO
    elseif(wbc.eq.2)then
#else
    if(ibw.eq.1.and.wbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=0,nj+1
        do i=1-ngxy,0
          w(i,j,k)=w(1,j,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibw.eq.1.and.(wbc.eq.3.or.wbc.eq.4))then
#else
    elseif(wbc.eq.3.or.wbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=0,nj+1
        do i=1,ngxy
          w(1-i,j,k)=w(i,j,k)
        enddo
        enddo
      ENDDO
    endif

!-----------------------------------------------------------------------
!  east boundary condition

#ifndef MPI
    if(ebc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=jb,je
        do i=1,ngxy
          w(ni+i,j,k)=w(i,j,k)
        enddo
        enddo
      ENDDO
    elseif(ebc.eq.2)then
#else
    if(ibe.eq.1.and.ebc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=0,nj+1
        do i=ni+1,ni+ngxy
          w(i,j,k)=w(ni,j,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibe.eq.1.and.(ebc.eq.3.or.ebc.eq.4))then
#else
    elseif(ebc.eq.3.or.ebc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=0,nj+1
        do i=1,ngxy
          w(ni+i,j,k)=w(ni+1-i,j,k)
        enddo
        enddo 
      ENDDO
    endif

!-----------------------------------------------------------------------
!  south boundary condition

#ifndef MPI
    if(sbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=1,ngxy
        do i=ib,ie
          w(i,1-j,k)=w(i,nj+1-j,k)
        enddo
        enddo
      ENDDO
    elseif(sbc.eq.2)then
#else
    if(ibs.eq.1.and.sbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=1-ngxy,0
        do i=0,ni+1
          w(i,j,k)=w(i,1,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibs.eq.1.and.(sbc.eq.3.or.sbc.eq.4))then
#else
    elseif(sbc.eq.3.or.sbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=1,ngxy
        do i=0,ni+1
          w(i,1-j,k)=w(i,j,k)
        enddo
        enddo
      ENDDO
    endif

!-----------------------------------------------------------------------
!  north boundary condition

#ifndef MPI
    if(nbc.eq.1)then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=1,ngxy
        do i=ib,ie
          w(i,nj+j,k)=w(i,j,k)
        enddo
        enddo
      ENDDO
    elseif(nbc.eq.2)then
#else
    if(ibn.eq.1.and.nbc.eq.2)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=nj+1,nj+ngxy
        do i=0,ni+1
          w(i,j,k)=w(i,nj,k)
        enddo
        enddo
      ENDDO
#ifdef MPI
    elseif(ibn.eq.1.and.(nbc.eq.3.or.nbc.eq.4))then
#else
    elseif(nbc.eq.3.or.nbc.eq.4)then
#endif
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) if(openacc)
      DO k=2,nk
        do j=1,ngxy
        do i=0,ni+1
          w(i,nj+j,k)=w(i,nj+1-j,k)
        enddo
        enddo
      ENDDO
    endif


!-----------------------------------------------------------------------
!  top/bottom boundary condition

    IF(flag.eq.1)THEN

      !$omp parallel do default(shared) private(i,j)
      !$acc parallel loop gang vector collapse(2) default(present) if(openacc)
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

      subroutine bcs2(s,device)
      use input, only : ib,ie,jb,je,kb,ke,ngxy,ni,nj,nk,wbc,ebc,sbc,nbc, &
          ibw,ibe,ibs,ibn,timestats,time_bcs2,mytime,patchsww,patchnww, &
          patchsee,patchnee,patchsws,patchses,patchnwn,patchnen
      implicit none

      real, intent(inout) ::  s(ib:ie,jb:je,kb:ke)
      logical, optional :: device

      integer i,j,k
      logical :: openacc

      call SetGPUParams(openacc,device)

!---------------------------------------------------------
!  This subroutine sets the corner points
!---------------------------------------------------------

      if(ibw.eq.1)then

        if(patchsww)then
          j=0
          !$omp parallel do default(shared)  private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do i=1-ngxy,0
                s(i,j,k)=s(1,j,k)
             end do
          end do
        endif

        if(patchnww)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do i=1-ngxy,0
                s(i,j,k)=s(1,j,k)
             end do
          end do
        endif

      endif

      if(ibe.eq.1)then

        if(patchsee)then
          j=0
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do i=ni+1,ni+ngxy
                s(i,j,k)=s(ni,j,k)
             end do
          end do
        endif

        if(patchnee)then
          j=nj+1
          !$omp parallel do default(shared) private(i,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do i=ni+1,ni+ngxy
                s(i,j,k)=s(ni,j,k)
             end do
          end do
        endif

      endif

      if(ibs.eq.1)then

        if(patchsws)then
          i=0
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do j=1-ngxy,0
                s(i,j,k)=s(i,1,k)
             end do
          end do
        endif

        if(patchses)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do j=1-ngxy,0
                s(i,j,k)=s(i,1,k)
             end do
          end do
        endif

      endif

      if(ibn.eq.1)then

        if(patchnwn)then
          i=0
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do j=nj+1,nj+ngxy
                s(i,j,k)=s(i,nj,k)
             end do
          end do
        endif

        if(patchnen)then
          i=ni+1
          !$omp parallel do default(shared) private(j,k)
          !$acc parallel loop gang vector default(present) collapse(2) if(openacc)
          do k=1,nk
             do j=nj+1,nj+ngxy
                s(i,j,k)=s(i,nj,k)
             end do
          end do
        endif

      endif
      
      if(timestats.ge.1) time_bcs2=time_bcs2+mytime()

      end subroutine bcs2

      subroutine SetGPUParams(openacc,device)

         logical, intent(inout) :: openacc
         logical, optional :: device

      ! Set logical flag if subroutine is called on GPU or CPU resident
      ! data.  If the device flag is not set, it is assumed that it is
      ! GPU resident.
      if(present(device)) then
         if(device) then
            openacc = .true.
         else
            openacc = .false.
         endif
      else
         openacc = .true.
      endif

      end subroutine SetGPUParams


  END MODULE bc_module
