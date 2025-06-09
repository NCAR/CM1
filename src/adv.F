  module adv_module

  implicit none

  private
  public :: advs,advu,advv,advw

  integer, parameter ::  gk00fix  =  3   !  1 = uvw terms
                                         !  2 = uv terms
                                         !  3 = uvw & uv terms

  contains

      subroutine advs(nrk,wflag,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,  &
                       rho0,rr0,rf0,rrf0,advx,advy,advz,dumx,dumy,dumz,mass,subs,divx,      &
                       rru,rrv,rrw,s0,s,sten,pdef,pdefweno,dt,weps,                         &
                       flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,ri,diffit, &
                       dobud,ibd,ied,jbd,jed,kbd,ked,ndiag,diag,sd_hadv,sd_vadv,sd_lsw,     &
                       sd_hidiff,sd_vidiff,sd_hediff,wprof,dumk1,dumk2,hadvorder,vadvorder,kbdy,bndy,hflxw,hflxe,hflxs,hflxn,out3d)
      use input
      use constants
      use pdef_module
      use adv_routines
      use ib_module
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: nrk
      integer, intent(in) :: wflag,bflag
      double precision, intent(inout) :: bsq
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: vh,rvh
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh,rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advx,advy,advz,dumx,dumy,dumz,mass,subs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s0,s
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: sten
      integer, intent(in) :: pdef,pdefweno
      real, intent(in) :: dt
      double precision, intent(in) :: weps
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,ri
      integer, intent(in) :: diffit
      logical, intent(in) :: dobud
      integer, intent(in) :: ibd,ied,jbd,jed,kbd,ked,ndiag,sd_hadv,sd_vadv,sd_lsw,sd_hidiff,sd_vidiff,sd_hediff
      real, intent(inout) , dimension(ibd:ied,jbd:jed,kbd:ked,ndiag) :: diag
      real, intent(in), dimension(kb:ke) :: wprof
      double precision, intent(inout), dimension(kb:ke) :: dumk1,dumk2
      integer, intent(in) :: hadvorder,vadvorder
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
 
      integer :: i,j,k,hadv
      logical :: doitw,doite,doits,doitn
      logical :: doweno
      real :: tem0,coef,tot,ndiff

      integer, dimension(4) :: reqsx,reqsy

!----------------------------------------------------------------

      doweno = .false.
      if( wflag.eq.1 )then
        if( (advwenos.eq.1) .or. (advwenos.eq.2.and.nrk.eq.nrkmax) ) doweno = .true.
      endif

      if(diffit.eq.1)then
        coef = kdiff6/64.0/dt
      endif

!-----------------

      doitw = wbc.eq.2 .and. ibw.eq.1
      doite = ebc.eq.2 .and. ibe.eq.1

!-----------------

      doits = sbc.eq.2 .and. ibs.eq.1
      doitn = nbc.eq.2 .and. ibn.eq.1

!-----------------

      if( bflag.eq.1 )then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          dumk1(k) = 0.0d0
          dumk2(k) = 0.0d0
        enddo
      endif

!-----------------

    hadvsection:  if(axisymm.eq.1)then
      call advsaxi(doweno,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,rmh,gz,rgz, &
                   rho0,rr0,rf0,rrf0,advx,dumx,mass,rru,s0,s,pdef,dt,weps,   &
                   hadvorder,flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    else

!----------------------------------------------------------------
! Advection in horizontal directions

    if( doweno )then

      if(     weno_order.eq.3 )then

        call     hadv_weno3(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s,pdef,weps)

      elseif( weno_order.eq.5 )then

        call     hadv_weno5(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s,pdef,weps)

      elseif( weno_order.eq.7 )then

        call     hadv_weno7(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s,pdef,weps)

      elseif( weno_order.eq.9 )then

        call     hadv_weno9(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s,pdef,weps)

      else

        print *,' 12941 '
        call stopcm1

      endif

    else

      if(     hadvorder.eq.2 )then

        call     hadv_flx2(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.3 )then

        call     hadv_flx3(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.4 )then

        call     hadv_flx4(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.5 )then

        call     hadv_flx5(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.6 )then

        call     hadv_flx6(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.7 )then

        call     hadv_flx7(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.8 )then

        call     hadv_flx8(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.9 )then

        call     hadv_flx9(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      elseif( hadvorder.eq.10 )then

        call     hadv_flx10(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)

      else

        print *,' 98611 '
        call stopcm1

      endif

    endif

    if( do_ib .and. hadvorder.ge.3 )then
      if( doweno )then
        call   ib_side_weno(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s,hadvorder,bndy,hflxw,hflxe,hflxs,hflxn,weps,out3d)
      else
        call   ib_side_flx(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s,hadvorder,bndy,hflxw,hflxe,hflxs,hflxn,out3d)
      endif
    endif

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow

      if(doitw)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k)=dumx(i+1,j,k)
          endif
          dumk1(k) = dumk1(k)+dumx(1,j,k)*rvh(j)*rmh(1,j,k)
        enddo
        enddo
      endif

      if(doite)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          if(rru(i,j,k).le.0.0)then
            dumx(i,j,k)=dumx(i-1,j,k)
          endif
          dumk1(k) = dumk1(k)-dumx(ni+1,j,k)*rvh(j)*rmh(ni+1,j,k)
        enddo
        enddo
      endif

      if(doits)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k)=dumy(i,j+1,k)
          endif
          dumk2(k) = dumk2(k)+dumy(i,1,k)*ruh(i)*rmh(i,1,k)
        enddo
        enddo
      endif

      if(doitn)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          if(rrv(i,j,k).le.0.0)then
            dumy(i,j,k)=dumy(i,j-1,k)
          endif
          dumk2(k) = dumk2(k)-dumy(i,nj+1,k)*ruh(i)*rmh(i,nj+1,k)
        enddo
        enddo
      endif

    !-------------------------------------------------------
    !  hadv tendencies:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
!!!        if( j.eq.41 .and. k.eq.20 ) print *,bndy(i,j,k),dumx(i,j,k),advx(i,j,k)
      enddo
      enddo
      enddo

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow
    !            (here, we subtract-off the divx piece)

      if(doitw)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          if(rru(1,j,k).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        enddo
      endif

      if(doite)then
        i=ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          if(rru(ni+1,j,k).le.0.0)then
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        enddo
      endif

      if(doits)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          if(rrv(i,1,k).ge.0.0)then
            advy(i,j,k)=advy(i,j,k)-s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
          endif
        enddo
        enddo
      endif

      if(doitn)then
        j=nj
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          if(rrv(i,nj+1,k).le.0.0)then
            advy(i,j,k)=advy(i,j,k)-s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    ! 6th-order diffusion-s:
    if(diffit.eq.1)then
      if( dobud .and. nrk.eq.nrkmax .and. sd_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          diag(i,j,k,sd_hediff) = (advx(i,j,k)+advy(i,j,k))
        enddo
        enddo
        enddo
      endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=1,nk
      do j=1,nj
      do i=1,ni+1
        dumx(i,j,k)=( 10.0*(s(i  ,j,k)-s(i-1,j,k))     &
                      -5.0*(s(i+1,j,k)-s(i-2,j,k))     &
                          +(s(i+2,j,k)-s(i-3,j,k)) )   &
                   *0.5*(rho(i-1,j,k)+rho(i,j,k))
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj
        do i=1,ni+1
          if( dumx(i,j,k)*(s(i,j,k)-s(i-1,j,k)).le.0.0 )then
            dumx(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=advx(i,j,k)+coef*(dumx(i+1,j,k)-dumx(i,j,k))*ri(i,j,k)*rho0(i,j,k)
      enddo
      enddo
      do j=1,nj+1
      do i=1,ni
        dumy(i,j,k)=( 10.0*(s(i,j  ,k)-s(i,j-1,k))     &
                      -5.0*(s(i,j+1,k)-s(i,j-2,k))     &
                          +(s(i,j+2,k)-s(i,j-3,k)) )   &
                   *0.5*(rho(i,j-1,k)+rho(i,j,k))
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni
          if( dumy(i,j,k)*(s(i,j,k)-s(i,j-1,k)).le.0.0 )then
            dumy(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=advy(i,j,k)+coef*(dumy(i,j+1,k)-dumy(i,j,k))*ri(i,j,k)*rho0(i,j,k)
      enddo
      enddo
    enddo
      if( dobud .and. nrk.eq.nrkmax .and. sd_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          diag(i,j,k,sd_hediff) = (advx(i,j,k)+advy(i,j,k))-diag(i,j,k,sd_hediff)
        enddo
        enddo
        enddo
      endif
    endif
    !-------------------------------------------------------


!----------------------------------------------------------------
!  Misc for x-direction

      if(stat_qsrc.eq.1.and.(wbc.eq.2.or.ebc.eq.2).and.bflag.eq.1)then
        tem0=dt*dy*dz
        do k=1,nk
          bsq=bsq+dumk1(k)*tem0
        enddo
      endif

      if(pdscheme.eq.1 .and. pdef.eq.1)then
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefx1(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dumx,mass,s0,s,dt,flag,sw31,sw32,se31,se32,reqsx)
      endif

!----------------------------------------------------------------
!  Misc for y-direction

      if(stat_qsrc.eq.1.and.(sbc.eq.2.or.nbc.eq.2).and.bflag.eq.1)then
        tem0=dt*dx*dz
        do k=1,nk
          bsq=bsq+dumk2(k)*tem0
        enddo
      endif

      if(pdscheme.eq.1 .and. pdef.eq.1)then
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefy1(vh,rho0,gz,rgz,rrv,advy,dumy,mass,s0,s,dt,flag,ss31,ss32,sn31,sn32,reqsy)
      endif

    endif  hadvsection

!----------------------------------------------------------------
! Advection in z-direction

    if( doweno )then

      if(     weno_order.eq.3 )then

        call     vadv_weno3(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s,pdef,weps)

      elseif( weno_order.eq.5 )then

        call     vadv_weno5(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s,pdef,pdefweno,weps)

      elseif( weno_order.eq.7 )then

        call     vadv_weno7(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s,pdef,weps)

      elseif( weno_order.eq.9 )then

        call     vadv_weno9(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s,pdef,weps)

      else

        print *,' 12942 '
        call stopcm1

      endif

    else

      if(     vadvorder.eq.2 )then

        call     vadv_flx2(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.3 )then

        call     vadv_flx3(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.4 )then

        call     vadv_flx4(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.5 )then

        call     vadv_flx5(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.6 )then

        call     vadv_flx6(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.7 )then

        call     vadv_flx7(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.8 )then

        call     vadv_flx8(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.9 )then

        call     vadv_flx9(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      elseif( vadvorder.eq.10 )then

        call     vadv_flx10(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)

      else

        print *,' 98612 '
        call stopcm1

      endif

    endif

    if( do_ib .and. vadvorder.ge.3 )then
      if( doweno )then
        call   ib_lwr_weno( 1 ,ni,nj,nk,c1,c2,rrw,dumz,s,kbdy,vadvorder,weps)
      else
        call   ib_lwr_flx(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s,kbdy,vadvorder)
      endif
    endif

!------

    if(terrain_flag)then

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        advz(i,j,1) = -dumz(i,j,2)*rdsf(1)
        advz(i,j,nk) = +dumz(i,j,nk)*rdsf(nk)
      enddo
      enddo

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk-1
      do j=1,nj
      do i=1,ni
        advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdsf(k)
      enddo
      enddo
      enddo

    else

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        advz(i,j,1) = -dumz(i,j,2)*rdz*mh(1,1,1)
        advz(i,j,nk) = +dumz(i,j,nk)*rdz*mh(1,1,nk)
      enddo
      enddo

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk-1
      do j=1,nj
      do i=1,ni
        advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
      enddo
      enddo
      enddo

    endif

!----------------------------------------------------------------
!  Large-scale subsidence:

    if( dobud .and. nrk.eq.nrkmax )then
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        subs(i,j,k) = 0.0
      enddo
      enddo
      enddo
    endif

    if( dolsw )then
      ! vertical advection from specified large-scale w profile:

      call     wsub(ni  ,nj  ,nk  ,s  ,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        ! Add subsidence tendency to advz array so it can be included 
        ! in pdefz calculation
        advz(i,j,k) = advz(i,j,k)+subs(i,j,k)*rho0(1,1,k)
      enddo
      enddo
      enddo

    endif

!----------------------------------------------------------------
!  Misc for z-direction

      if(pdscheme.eq.1 .and. pdef.eq.1)then
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefz(mh,rho0,gz,rgz,rdsf,rrw,advz,dumz,mass,s0,s,dt,flag)
      endif

!----------------------------------------------------------------
!  Finish pdefxy:

      if(pdscheme.eq.1 .and. pdef.eq.1 .and. axisymm.eq.0)then
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefx2(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dumx,mass,s0,s,dt,flag,sw31,sw32,se31,se32,reqsx)
        call pdefy2(vh,rho0,gz,rgz,rrv,advy,dumy,mass,s0,s,dt,flag,ss31,ss32,sn31,sn32,reqsy)
      endif

!----------------------------------------------------------------
!  Total advection tendency:

    if(terrain_flag)then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                                 +s(i,j,k)*divx(i,j,k) )*rr0(i,j,k)*gz(i,j)
      enddo
      enddo
      enddo

    else

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                                 +s(i,j,k)*divx(i,j,k) )*rr0(1,1,k)
      enddo
      enddo
      enddo

    endif

!----------------------------------------------------------------
!  Diagnostics:

  doingsbud:  &
  if( dobud .and. nrk.eq.nrkmax )then
  if( sd_hadv.ge.1 .and. sd_vadv.ge.1 )then
    !--------
    ! advective tendencies:
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          diag(i,j,k,sd_hadv) = ( advx(i,j,k)+advy(i,j,k)             &
                + s(i,j,k)*( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) &
                            +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) ) )*rr0(1,1,k)
          diag(i,j,k,sd_vadv) = ( (advz(i,j,k)-subs(i,j,k)*rho0(1,1,k))         &
                + s(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) ) )*rr0(1,1,k)
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k)
        do k=1,nk
        do i=1,ni
          diag(i,j,k,sd_hadv) = ( advx(i,j,k)                         &
                + s(i,j,k)*( (arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i) )  )*rr0(1,1,k)
          diag(i,j,k,sd_vadv) = ( advz(i,j,k)                         &
                + s(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) )  )*rr0(1,1,k)
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          diag(i,j,k,sd_hadv) = ( advx(i,j,k)+advy(i,j,k)             &
                + s(i,j,k)*( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) &
                            +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) ) )*rr0(i,j,k)*gz(i,j)
          diag(i,j,k,sd_vadv) = ( advz(i,j,k)                         &
                + s(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k) ) )*rr0(i,j,k)*gz(i,j)
        enddo
        enddo
        enddo
    endif
  endif
    !--------
    if( sd_lsw.ge.1 )then
      ! large-scale subsidence:
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        diag(i,j,k,sd_lsw) = subs(i,j,k)
      enddo
      enddo
      enddo
    endif
    !--------
    if( diffit.eq.1 .and. sd_hediff.ge.1 .and. sd_hadv.ge.1 )then
      ! subtract-off diffusion from advection:
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        diag(i,j,k,sd_hadv) = diag(i,j,k,sd_hadv)-diag(i,j,k,sd_hediff)
      enddo
      enddo
      enddo
    endif
    !--------
    gethidiffs:  &
    if( sd_hidiff.ge.1 .and. sd_hadv.ge.1 )then
      ! horiz implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh31a !
    if( axisymm.eq.0 )then
      if(     hadvorder.eq.3 .or. ( advwenos.ge.1 .and. weno_order.eq.3 ) )then
        call     hadv_flx4(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)
      elseif( hadvorder.eq.5 .or. ( advwenos.ge.1 .and. weno_order.eq.5 ) )then
        call     hadv_flx6(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)
      elseif( hadvorder.eq.7 .or. ( advwenos.ge.1 .and. weno_order.eq.7 ) )then
        call     hadv_flx8(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)
      elseif( hadvorder.eq.9 .or. ( advwenos.ge.1 .and. weno_order.eq.9 ) )then
        call     hadv_flx10(  1 ,ni,nj,nk,c1,c2,rru,rrv,dumx,dumy,s)
      else
        print *,' 13941 '
        call stopcm1
      endif
    else
      if(     hadvorder.eq.3 .or. ( advwenos.ge.1 .and. weno_order.eq.3 ) )then
        hadv = 4
      elseif( hadvorder.eq.5 .or. ( advwenos.ge.1 .and. weno_order.eq.5 ) )then
        hadv = 6
      elseif( hadvorder.eq.7 .or. ( advwenos.ge.1 .and. weno_order.eq.7 ) )then
        hadv = 8
      elseif( hadvorder.eq.9 .or. ( advwenos.ge.1 .and. weno_order.eq.9 ) )then
        hadv = 10
      else
        print *,' 13981 '
        call stopcm1
      endif
      call advsaxi(.false.,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,rmh,gz,rgz, &
                   rho0,rr0,rf0,rrf0,advx,dumx,mass,rru,s0,s, 0  ,dt,weps,   &
                   hadv     ,flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
    endif
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh41 !
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tot = diag(i,j,k,sd_hadv)
          advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
          ndiff               = ( advx(i,j,k)+advy(i,j,k)             &
                + s(i,j,k)*( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) &
                            +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) ) )*rr0(1,1,k)
          diag(i,j,k,sd_hadv) = ndiff
          diag(i,j,k,sd_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k,tot,ndiff)
        do k=1,nk
        do i=1,ni
          tot = diag(i,j,k,sd_hadv)
          ndiff               = ( advx(i,j,k)                         &
                + s(i,j,k)*( (arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i) )  )*rr0(1,1,k)
          diag(i,j,k,sd_hadv) = ndiff
          diag(i,j,k,sd_hidiff) = tot-ndiff
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tot = diag(i,j,k,sd_hadv)
          advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
          ndiff               = ( advx(i,j,k)+advy(i,j,k)             &
                + s(i,j,k)*( (rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) &
                            +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) ) )*rr0(i,j,k)*gz(i,j)
          diag(i,j,k,sd_hadv) = ndiff
          diag(i,j,k,sd_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  gethidiffs
    !--------
    getvidiffs:  &
    if( sd_vidiff.ge.1 .and. sd_vadv.ge.1 )then
      ! vert implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh32a !
      if(     vadvorder.eq.3 .or. ( advwenos.ge.1 .and. weno_order.eq.3 ) )then
        call     vadv_flx4(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)
      elseif( vadvorder.eq.5 .or. ( advwenos.ge.1 .and. weno_order.eq.5 ) )then
        call     vadv_flx6(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)
      elseif( vadvorder.eq.7 .or. ( advwenos.ge.1 .and. weno_order.eq.7 ) )then
        call     vadv_flx8(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)
      elseif( vadvorder.eq.9 .or. ( advwenos.ge.1 .and. weno_order.eq.9 ) )then
        call     vadv_flx10(  1 ,ni,nj,nk,c1,c2,rrw,dumz,s)
      else
        print *,' 13942 '
        call stopcm1
      endif
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        dumz(i,j,1) = 0.0
        dumz(i,j,nk+1) = 0.0
      enddo
      enddo
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh42a !
    if(.not.terrain_flag)then
        ! no terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tot = diag(i,j,k,sd_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          ndiff               = ( (advz(i,j,k)-subs(i,j,k)*rho0(1,1,k))         &
                + s(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) ) )*rr0(1,1,k)
          diag(i,j,k,sd_vadv) = ndiff
          diag(i,j,k,sd_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tot = diag(i,j,k,sd_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdsf(k)
          ndiff               = ( advz(i,j,k)                         &
                + s(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k) ) )*rr0(i,j,k)*gz(i,j)
          diag(i,j,k,sd_vadv) = ndiff
          diag(i,j,k,sd_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  getvidiffs
    !--------
  endif  doingsbud

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advs=time_advs+mytime()
 
      end subroutine advs


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advu(nrk,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,gz,rgz,gzu,mh,rho0,rr0,rf0,rrf0,dumx,dumy,dumz,advx,advy,advz,subs,divx, &
                       rru,u3d,uten,rrv,rrw,rdsf,c1,c2,rho,dt,doubud,udiag,wprof,bndy,kbdy,hflxw,hflxe,hflxs,hflxn,vf,mf,v3d,w3d)
      use input
      use constants
      use adv_routines
      use ib_module
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: nrk
      real, intent(in), dimension(ib:ie) :: arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,dumz,advx,advy,advz,subs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(in) :: dt
      logical, intent(in) :: doubud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
      real, intent(in), dimension(kb:ke) :: wprof
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
 
      integer :: i,j,k,i1,i2,j1,j2,id1,id2,hadv
      real :: ubar,vbar,cc1,cc2
      logical :: doitw,doite,doits,doitn
      logical :: doweno
      double precision :: weps
      real :: coef,tot,ndiff

!------------------------------------------------------------

      doweno = .false.
      if( (advwenov.eq.1) .or. (advwenov.eq.2.and.nrk.eq.nrkmax) ) doweno = .true.
      weps = 100.0*epsilon

      if( idiff.ge.1 .and. difforder.eq.6 )then
        coef = kdiff6/64.0/dt
      endif

!-----------------

      if(ibw.eq.1)then
        i1=2
      else
        i1=1
      endif
 
      if(ibe.eq.1)then
        i2=ni+1-1
      else
        i2=ni+1
      endif

      id1 = i1-1
      id2 = i2

      doitw = .false.
      doite = .false.

      if(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      if(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

      j1 = 1
      j2 = nj+1

      doits = .false.
      doitn = .false.

      if(sbc.eq.2 .and. ibs.eq.1) doits = .true.
      if(nbc.eq.2 .and. ibn.eq.1) doitn = .true.

!----------------------------------------------------------------

    hadvsection:  if(axisymm.eq.1)then
      call advuaxi(doweno,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,rho0,rr0,rf0,rrf0,dumx,advx,rru,u3d,hadvordrv)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    else

!----------------------------------------------------------------
! Advection in horizontal directions

    if( doweno )then

      if(     weno_order.eq.3 )then

        call     hadv_weno3(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d,0,weps)

      elseif( weno_order.eq.5 )then

        call     hadv_weno5(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d,0,weps)

      elseif( weno_order.eq.7 )then

        call     hadv_weno7(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d,0,weps)

      elseif( weno_order.eq.9 )then

        call     hadv_weno9(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d,0,weps)

      else

        print *,' 12943 '
        call stopcm1

      endif

    else

      if(     hadvordrv.eq.2 )then

        call     hadv_flx2(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.3 )then

        call     hadv_flx3(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.4 )then

        call     hadv_flx4(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.5 )then

        call     hadv_flx5(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.6 )then

        call     hadv_flx6(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.7 )then

        call     hadv_flx7(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.8 )then

        call     hadv_flx8(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.9 )then

        call     hadv_flx9(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      elseif( hadvordrv.eq.10 )then

        call     hadv_flx10(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)

      else

        print *,' 98613 '
        call stopcm1

      endif

    endif

    if( do_ib .and. hadvordrv.ge.3 )then
      if( doweno )then
        call   ib_side_weno(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d,hadvordrv,bndy,hflxw,hflxe,hflxs,hflxn,weps)
      else
        call   ib_side_flx(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d,hadvordrv,bndy,hflxw,hflxe,hflxs,hflxn)
      endif
    endif

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow

      if(doits)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=i1,i2
          if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
            dumy(i,j,k)=dumy(i,j+1,k)
          endif
        enddo
        enddo
      endif

      if(doitn)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=i1,i2
          if((rrv(i,j,k)+rrv(i-1,j,k)).le.0.0)then
            dumy(i,j,k)=dumy(i,j-1,k)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    !  hadv tendencies:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advx(i,j,k) = -(dumx(i,j,k)-dumx(i-1,j,k))*rdx*uf(i)
        advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow
    !            (here, we subtract-off the divx piece)

      if(doits)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=i1,i2
          if((rrv(i,1,k)+rrv(i-1,1,k)).ge.0.0)then
            advy(i,j,k)=advy(i,j,k)-u3d(i,j,k)*0.5*(                    &
                            (rrv(i-1,j+1,k)-rrv(i-1,j,k))               &
                           +(rrv(i  ,j+1,k)-rrv(i  ,j,k)) )*rdy*vh(j)
          endif
        enddo
        enddo
      endif

      if(doitn)then
        j=nj
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=i1,i2
          if((rrv(i,nj+1,k)+rrv(i-1,nj+1,k)).le.0.0)then
            advy(i,j,k)=advy(i,j,k)-u3d(i,j,k)*0.5*(                    &
                            (rrv(i-1,j+1,k)-rrv(i-1,j,k))               &
                           +(rrv(i  ,j+1,k)-rrv(i  ,j,k)) )*rdy*vh(j)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    ! 6th-order diffusion-u:
    if( idiff.ge.1 .and. difforder.eq.6 )then
      if( doubud .and. nrk.eq.nrkmax .and. ud_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          udiag(i,j,k,ud_hediff) = (advx(i,j,k)+advy(i,j,k))
        enddo
        enddo
        enddo
      endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=1,nk
      do j=1,nj
      do i=1,ni+2
        dumx(i,j,k)=( 10.0*(u3d(i  ,j,k)-u3d(i-1,j,k))     &
                      -5.0*(u3d(i+1,j,k)-u3d(i-2,j,k))     &
                          +(u3d(i+2,j,k)-u3d(i-3,j,k)) )*rho(i-1,j,k)
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj
        do i=1,ni+2
          if( dumx(i,j,k)*(u3d(i,j,k)-u3d(i-1,j,k)).le.0.0 )then
            dumx(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni+1
        advx(i,j,k)=advx(i,j,k)+coef*(dumx(i+1,j,k)-dumx(i,j,k))*(rho0(i-1,j,k)+rho0(i,j,k))/(rho(i-1,j,k)+rho(i,j,k))
      enddo
      enddo
      do j=1,nj+1
      do i=1,ni+1
        dumy(i,j,k)=( 10.0*(u3d(i,j  ,k)-u3d(i,j-1,k))     &
                      -5.0*(u3d(i,j+1,k)-u3d(i,j-2,k))     &
                          +(u3d(i,j+2,k)-u3d(i,j-3,k)) )   &
                  *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))     &
                         +(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni+1
          if( dumy(i,j,k)*(u3d(i,j,k)-u3d(i,j-1,k)).le.0.0 )then
            dumy(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni+1
        advy(i,j,k)=advy(i,j,k)+coef*(dumy(i,j+1,k)-dumy(i,j,k))*(rho0(i-1,j,k)+rho0(i,j,k))/(rho(i-1,j,k)+rho(i,j,k))
      enddo
      enddo
    enddo
      if( doubud .and. nrk.eq.nrkmax .and. ud_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          udiag(i,j,k,ud_hediff) = (advx(i,j,k)+advy(i,j,k))-udiag(i,j,k,ud_hediff)
        enddo
        enddo
        enddo
      endif
    endif
    !-------------------------------------------------------


    endif  hadvsection

!----------------------------------------------------------------
! Advection in z-direction  (Cartesian grid)

  vadvu:  if(axisymm.eq.0)then

    if( doweno )then

      if(     weno_order.eq.3 )then

        call     vadv_weno3(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d,0,weps)

      elseif( weno_order.eq.5 )then

        call     vadv_weno5(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d,0,0,weps)

      elseif( weno_order.eq.7 )then

        call     vadv_weno7(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d,0,weps)

      elseif( weno_order.eq.9 )then

        call     vadv_weno9(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d,0,weps)

      else

        print *,' 12944 '
        call stopcm1

      endif

    else

      if(     vadvordrv.eq.2 )then

        call     vadv_flx2(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.3 )then

        call     vadv_flx3(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.4 )then

        call     vadv_flx4(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.5 )then

        call     vadv_flx5(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.6 )then

        call     vadv_flx6(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.7 )then

        call     vadv_flx7(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.8 )then

        call     vadv_flx8(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.9 )then

        call     vadv_flx9(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      elseif( vadvordrv.eq.10 )then

        call     vadv_flx10(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)

      else

        print *,' 98614 '
        call stopcm1

      endif

    endif

    if( do_ib .and. vadvordrv.ge.3 )then
      if( doweno )then
        call   ib_lwr_weno( 2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d,kbdy,vadvordrv,weps)
      else
        call   ib_lwr_flx(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d,kbdy,vadvordrv)
      endif
    endif

!------

    if(terrain_flag)then

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do j=1,nj
        do i=i1,i2
          k=1
          advz(i,j,k) = -dumz(i,j,k+1)*rdsf(k)
          uten(i,j,k) = uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )         &
                  *gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
          k=nk
          advz(i,j,k) = +dumz(i,j,k)*rdsf(k)
          uten(i,j,k) = uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )         &
                  *gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj
        do i=i1,i2
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdsf(k)
          uten(i,j,k) = uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )         &
                  *gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
        enddo
        enddo
        enddo

    else

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do j=1,nj
        do i=i1,i2
          k=1
          advz(i,j,k) = -dumz(i,j,k+1)*rdz*mh(1,1,k)
          k=nk
          advz(i,j,k) = +dumz(i,j,k)*rdz*mh(1,1,k)
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj
        do i=i1,i2
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
        enddo
        enddo
        enddo

          !-----------------------------------------------------------------------------------------------------------
        if( do_ib .and. ( gk00fix.eq.1 .or. gk00fix.eq.3 ) )then
          ! apply fix from Gallus and Klemp (2000, MWR, pg 1160)
          ! replace -w du/dz with -w dw/dx
          do k=1,kmaxib
          do j=1,nj
          do i=i1,i2
            if( bndy(i,j,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i-1,j,k)) .and. (.not.bndy(i-1,j,k-1)) )then
              ! upper-left corner:
              ! Note: advective form is used here, so remove the piece of advz from divx
              advz(i,j,k) = 0.25*( -(rrw(i-1,j,k+1)+rrw(i,j,k+1))*(u3d(i,j,k+1)-u3d(i,j,k  ))*rdz*mf(1,1,k+1)    &
                                   -(rrw(i-1,j,k  )+rrw(i,j,k  ))*(w3d(i,j,k)-w3d(i-1,j,k))*rdx*uf(i) )  &
                            -u3d(i,j,k)*0.5*( (rrw(i  ,j,k+1)-rrw(i  ,j,k))*rdz*mh(1,1,k) &
                                             +(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdz*mh(1,1,k) )
            endif
            if( bndy(i-1,j,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i-1,j,k)) .and. (.not.bndy(i,j,k-1)) )then
              ! upper-right corner:
              ! Note: advective form is used here, so remove the piece of advz from divx
              advz(i,j,k) = 0.25*( -(rrw(i-1,j,k+1)+rrw(i,j,k+1))*(u3d(i,j,k+1)-u3d(i,j,k  ))*rdz*mf(1,1,k+1)    &
                                   -(rrw(i-1,j,k  )+rrw(i,j,k  ))*(w3d(i,j,k)-w3d(i-1,j,k))*rdx*uf(i) )  &
                            -u3d(i,j,k)*0.5*( (rrw(i  ,j,k+1)-rrw(i  ,j,k))*rdz*mh(1,1,k) &
                                             +(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdz*mh(1,1,k) )
            endif
          enddo
          enddo
          enddo
        endif
          !-----------------------------------------------------------------------------------------------------------
        if( do_ib .and. ( gk00fix.eq.2 .or. gk00fix.eq.3 ) )then
          ! apply fix from Gallus and Klemp (2000, MWR, pg 1160)
          ! replace -w du/dz with -w dw/dx
          do k=1,kmaxib
          do j=1,nj
          do i=i1,i2
            if( bndy(i,j+1,k) .and. (.not.bndy(i-1,j+1,k)) .and. (.not.bndy(i-1,j,k)) .and. (.not.bndy(i,j,k)) )then
              ! SW corner:  replace -v du/dy with -v dv/dx:
!!!              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i-1,j+1,k))*(v3d(i,j+1,k)-v3d(i-1,j+1,k))*rdx*uf(i)     &
              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i-1,j+1,k))*(u3d(i,j+1,k)-u3d(i,j  ,k))*rdy*vf(j+1)     &
                                   -(rrv(i,j  ,k)+rrv(i-1,j  ,k))*(u3d(i,j  ,k)-u3d(i,j-1,k))*rdy*vf(j  )  )  &
                            -u3d(i,j,k)*0.5*( (rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j)  &
                                             +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j) )
            endif
            if( (.not.bndy(i,j+1,k)) .and. bndy(i-1,j+1,k) .and. (.not.bndy(i-1,j,k)) .and. (.not.bndy(i,j,k)) )then
              ! SE corner:  replace -v du/dy with -v dv/dx:
!!!              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i-1,j+1,k))*(v3d(i,j+1,k)-v3d(i-1,j+1,k))*rdx*uf(i)     &
              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i-1,j+1,k))*(u3d(i,j+1,k)-u3d(i,j  ,k))*rdy*vf(j+1)     &
                                   -(rrv(i,j  ,k)+rrv(i-1,j  ,k))*(u3d(i,j  ,k)-u3d(i,j-1,k))*rdy*vf(j  )  )  &
                            -u3d(i,j,k)*0.5*( (rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j)  &
                                             +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j) )
            endif
            if( (.not.bndy(i-1,j,k)) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i-1,j-1,k)) .and. bndy(i,j-1,k) )then
              ! NW corner:  replace -v du/dy with -v dv/dx:
              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i-1,j+1,k))*(u3d(i,j+1,k)-u3d(i,j  ,k))*rdy*vf(j+1)     &
                                   -(rrv(i,j  ,k)+rrv(i-1,j  ,k))*(u3d(i,j  ,k)-u3d(i,j-1,k))*rdy*vf(j  )  )  &
!!!                                   -(rrv(i,j  ,k)+rrv(i-1,j  ,k))*(v3d(i,j  ,k)-v3d(i-1,j,k))*rdx*uf(i)    )  &
                            -u3d(i,j,k)*0.5*( (rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j)  &
                                             +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j) )
            endif
            if( (.not.bndy(i-1,j,k)) .and. (.not.bndy(i,j,k)) .and. bndy(i-1,j-1,k) .and. (.not.bndy(i,j-1,k)) )then
              ! NE corner:  replace -v du/dy with -v dv/dx:
              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i-1,j+1,k))*(u3d(i,j+1,k)-u3d(i,j  ,k))*rdy*vf(j+1)     &
                                   -(rrv(i,j  ,k)+rrv(i-1,j  ,k))*(u3d(i,j  ,k)-u3d(i,j-1,k))*rdy*vf(j  )  )  &
!!!                                   -(rrv(i,j  ,k)+rrv(i-1,j  ,k))*(v3d(i,j  ,k)-v3d(i-1,j,k))*rdx*uf(i)    )  &
                            -u3d(i,j,k)*0.5*( (rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j)  &
                                             +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j) )
            endif
          enddo
          enddo
          enddo
        endif
          !-----------------------------------------------------------------------------------------------------------

      ! cm1ib: moved uten

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          uten(i,j,k) = uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        enddo
        enddo

    endif

!----------------------------------------------------------------

    if( dolsw )then
      ! vertical advection from specified large-scale w profile:

      call     wsub(ni+1,nj  ,nk  ,u3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)

    endif

!  end vadvu for Cartesian grid
!----------------------------------------------------------------
! Advection in z-direction  (axisymmetric grid)

  elseif(axisymm.eq.1)then  vadvu

      if(ebc.eq.3.or.ebc.eq.4) i2 = ni

      call       vadv_axiu(doweno,arf1,arf2,c1,c2,rrw,u3d,dumz,weps)

!------

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do j=1,nj
        do i=2,i2
          k=1
          advz(i,j,k) = -dumz(i,j,k+1)*rdz*mh(1,1,k)
          uten(i,j,k) = uten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(arf2(i)*divx(i,j,k)+arf1(i)*divx(i-1,j,k)) )*rr0(1,1,k)
          k=nk
          advz(i,j,k) = +dumz(i,j,k)*rdz*mh(1,1,k)
          uten(i,j,k) = uten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(arf2(i)*divx(i,j,k)+arf1(i)*divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk-1
        do j=1,nj
        do i=2,i2
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          uten(i,j,k) = uten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(arf2(i)*divx(i,j,k)+arf1(i)*divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        enddo
        enddo

  else
    print *,' 54525 '
    call stopcm1
  endif  vadvu

!  end vadvu for axisymmetric grid
!----------------------------------------------------------------

!----------------------------------------------------------------
!  Diagnostics:

  doingubud:  &
  if( doubud .and. nrk.eq.nrkmax )then
  if( ud_hadv.ge.1 .and. ud_vadv.ge.1 )then
    !--------
    ! advective tendencies:
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          udiag(i,j,k,ud_hadv) = ( advx(i,j,k)+advy(i,j,k)                        &
                + u3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i  ,j,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j  ) )  &
                                  +( (rru(i  ,j  ,k)-rru(i-1,j,k))*rdx*uh(i-1)    &
                                    +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j  ) )  &
                                 ) )*rr0(1,1,k)
          udiag(i,j,k,ud_vadv) = ( advz(i,j,k)                                  &
                + u3d(i,j,k)*0.5*( (rrw(i  ,j,k+1)-rrw(i  ,j,k))*rdz*mh(1,1,k)  &
                                  +(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdz*mh(1,1,k)  &
                                 ) )*rr0(1,1,k)
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k)
        do k=1,nk
        do i=2,i2
          udiag(i,j,k,ud_hadv) = ( advx(i,j,k)                         &
                + u3d(i,j,k)*0.5*( arf2(i)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)     &
                                  +arf1(i)*(arh2(i-1)*rru(i,j,k)-arh1(i-1)*rru(i-1,j,k))*rdx*uh(i-1) )    &
                                 )*rr0(1,1,k)
          udiag(i,j,k,ud_vadv) = ( advz(i,j,k)                         &
                + u3d(i,j,k)*0.5*( arf2(i)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)     &
                                  +arf1(i)*(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdz*mh(1,1,k) )    &
                                 )*rr0(1,1,k)
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          udiag(i,j,k,ud_hadv) = ( advx(i,j,k)+advy(i,j,k)                        &
                + u3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i  ,j,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j  ) )  &
                                  +( (rru(i  ,j  ,k)-rru(i-1,j,k))*rdx*uh(i-1)    &
                                    +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j  ) )  &
                                 ) )*gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
          udiag(i,j,k,ud_vadv) = ( advz(i,j,k)                                  &
                + u3d(i,j,k)*0.5*( (rrw(i  ,j,k+1)-rrw(i  ,j,k))*rdsf(k)        &
                                  +(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdsf(k)        &
                                 ) )*gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
        enddo
        enddo
        enddo
    endif
  endif
    !--------
    if( ud_lsw.ge.1 )then
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        udiag(i,j,k,ud_lsw) = subs(i,j,k)
      enddo
      enddo
      enddo
    endif
    !--------
    if( ud_hediff.ge.1 .and. ud_hadv.ge.1 )then
    if( idiff.ge.1 .and. difforder.eq.6 )then
      ! subtract-off diffusion from advection:
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        udiag(i,j,k,ud_hadv) = udiag(i,j,k,ud_hadv)-udiag(i,j,k,ud_hediff)
      enddo
      enddo
      enddo
    endif
    endif
    !--------
    gethidiffu:  &
    if( ud_hidiff.ge.1 .and. ud_hadv.ge.1 )then
      ! horiz implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh31a !
    if( axisymm.eq.0 )then
      if(     hadvordrv.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        call     hadv_flx4(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)
      elseif( hadvordrv.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        call     hadv_flx6(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)
      elseif( hadvordrv.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        call     hadv_flx8(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)
      elseif( hadvordrv.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        call     hadv_flx10(  2 ,ni+1,nj,nk,c1,c2,rru,rrv,dumx,dumy,u3d)
      else
        print *,' 13951 '
        call stopcm1
      endif
    else
      if(     hadvordrv.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        hadv = 4
      elseif( hadvordrv.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        hadv = 6
      elseif( hadvordrv.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        hadv = 8
      elseif( hadvordrv.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        hadv = 10
      else
        print *,' 13981 '
        call stopcm1
      endif
      call advuaxi(.false.,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,rho0,rr0,rf0,rrf0,dumx,advx,rru,u3d,hadv)
    endif
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh41 !
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          tot = udiag(i,j,k,ud_hadv)
          advx(i,j,k) = -(dumx(i,j,k)-dumx(i-1,j,k))*rdx*uf(i)
          advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
          ndiff                = ( advx(i,j,k)+advy(i,j,k)                        &
                + u3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i  ,j,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j  ) )  &
                                  +( (rru(i  ,j  ,k)-rru(i-1,j,k))*rdx*uh(i-1)    &
                                    +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j  ) )  &
                                 ) )*rr0(1,1,k)
          udiag(i,j,k,ud_hadv) = ndiff
          udiag(i,j,k,ud_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k,tot,ndiff)
        do k=1,nk
        do i=2,i2
          tot = udiag(i,j,k,ud_hadv)
          ndiff                = ( advx(i,j,k)                         &
                + u3d(i,j,k)*0.5*( arf2(i)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)     &
                                  +arf1(i)*(arh2(i-1)*rru(i,j,k)-arh1(i-1)*rru(i-1,j,k))*rdx*uh(i-1) )    &
                                 )*rr0(1,1,k)
          udiag(i,j,k,ud_hadv) = ndiff
          udiag(i,j,k,ud_hidiff) = tot-ndiff
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          tot = udiag(i,j,k,ud_hadv)
          advx(i,j,k) = -(dumx(i,j,k)-dumx(i-1,j,k))*rdx*uf(i)
          advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
          ndiff                = ( advx(i,j,k)+advy(i,j,k)                        &
                + u3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i  ,j,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i  ,j,k))*rdy*vh(j  ) )  &
                                  +( (rru(i  ,j  ,k)-rru(i-1,j,k))*rdx*uh(i-1)    &
                                    +(rrv(i-1,j+1,k)-rrv(i-1,j,k))*rdy*vh(j  ) )  &
                                 ) )*gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
          udiag(i,j,k,ud_hadv) = ndiff
          udiag(i,j,k,ud_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  gethidiffu
    !--------
    getvidiffu:  &
    if( ud_vidiff.ge.1 .and. ud_vadv.ge.1 )then
      ! vert implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh32a !
      if(     vadvordrv.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        call     vadv_flx4(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)
      elseif( vadvordrv.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        call     vadv_flx6(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)
      elseif( vadvordrv.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        call     vadv_flx8(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)
      elseif( vadvordrv.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        call     vadv_flx10(  2 ,ni+1,nj,nk,c1,c2,rrw,dumz,u3d)
      else
        print *,' 13952 '
        call stopcm1
      endif
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      do i=i1,i2
        dumz(i,j,1) = 0.0
        dumz(i,j,nk+1) = 0.0
      enddo
      enddo
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh42a !
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          tot = udiag(i,j,k,ud_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          ndiff                = ( advz(i,j,k)                                  &
                + u3d(i,j,k)*0.5*( (rrw(i  ,j,k+1)-rrw(i  ,j,k))*rdz*mh(1,1,k)  &
                                  +(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdz*mh(1,1,k)  &
                                 ) )*rr0(1,1,k)
          udiag(i,j,k,ud_vadv) = ndiff
          udiag(i,j,k,ud_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k,tot,ndiff)
        do k=1,nk
        do i=2,i2
          tot = udiag(i,j,k,ud_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          ndiff                = ( advz(i,j,k)                         &
                + u3d(i,j,k)*0.5*( arf2(i)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)     &
                                  +arf1(i)*(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdz*mh(1,1,k) )    &
                                 )*rr0(1,1,k)
          udiag(i,j,k,ud_vadv) = ndiff
          udiag(i,j,k,ud_vidiff) = tot-ndiff
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=1,nj
        do i=i1,i2
          tot = udiag(i,j,k,ud_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdsf(k)
          ndiff                = ( advz(i,j,k)                                  &
                + u3d(i,j,k)*0.5*( (rrw(i  ,j,k+1)-rrw(i  ,j,k))*rdsf(k)        &
                                  +(rrw(i-1,j,k+1)-rrw(i-1,j,k))*rdsf(k)        &
                                 ) )*gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
          udiag(i,j,k,ud_vadv) = ndiff
          udiag(i,j,k,ud_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  getvidiffu
    !--------
  endif  doingubud

!----------------------------------------------------------------

    if( dolsw )then
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        uten(i,j,k) = uten(i,j,k)+subs(i,j,k)
      enddo
      enddo
      enddo
    endif

!----------------------------------------------------------------

      if(timestats.ge.1) time_advu=time_advu+mytime()
 
      end subroutine advu


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advv(nrk,xh,rxh,arh1,arh2,uh,xf,vh,vf,gz,rgz,gzv,mh,rho0,rr0,rf0,rrf0,dumx,dumy,dumz,advx,advy,advz,subs,divx, &
                       rru,rrv,v3d,vten,rrw,rdsf,c1,c2,rho,dt,dovbud,vdiag,wprof,bndy,kbdy,hflxw,hflxe,hflxs,hflxn,uf,mf,u3d,w3d)
      use input
      use constants
      use adv_routines
      use ib_module
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: nrk
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,dumz,advx,advy,advz,subs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(in) :: dt
      logical, intent(in) :: dovbud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
      real, intent(in), dimension(kb:ke) :: wprof
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
 
      integer :: i,j,k,j1,j2,jd1,jd2,hadv
      real :: ubar,vbar,cc1,cc2
      logical :: doitw,doite,doits,doitn
      logical :: doweno
      double precision :: weps
      real :: coef,tot,ndiff

!------------------------------------------------------------

      doweno = .false.
      if( (advwenov.eq.1) .or. (advwenov.eq.2.and.nrk.eq.nrkmax) ) doweno = .true.
      weps = 100.0*epsilon

      if( idiff.ge.1 .and. difforder.eq.6 )then
        coef = kdiff6/64.0/dt
      endif

!-----------------

      doitw = .false.
      doite = .false.

      if(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      if(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

      if(ibs.eq.1)then
        j1=2
      else
        j1=1
      endif
 
      if(ibn.eq.1)then
        j2=nj+1-1
      else
        j2=nj+1
      endif

      jd1 = j1-1
      jd2 = j2

      doits = .false.
      doitn = .false.

      if(sbc.eq.2 .and. ibs.eq.1) doits = .true.
      if(nbc.eq.2 .and. ibn.eq.1) doitn = .true.

!----------------------------------------------------------------

    hadvsection:  if(axisymm.eq.1)then

      j1 = 1
      j2 = 1
      ! advz stores M
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,1
      do i=0,ni+1
        advz(i,j,k) = xh(i)*( v3d(i,j,k) + 0.5*fcor*xh(i) )
      enddo
      enddo
      enddo
      call advvaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vf,rho0,rr0,rf0,rrf0,dumx,advx,advz,rru,hadvordrv)

    else

!-----------------

    if( doweno )then

      if(     weno_order.eq.3 )then

        call     hadv_weno3(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d,0,weps)

      elseif( weno_order.eq.5 )then

        call     hadv_weno5(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d,0,weps)

      elseif( weno_order.eq.7 )then

        call     hadv_weno7(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d,0,weps)

      elseif( weno_order.eq.9 )then

        call     hadv_weno9(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d,0,weps)

      else

        print *,' 12945 '
        call stopcm1

      endif

    else

      if(     hadvordrv.eq.2 )then

        call     hadv_flx2(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.3 )then

        call     hadv_flx3(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.4 )then

        call     hadv_flx4(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.5 )then

        call     hadv_flx5(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.6 )then

        call     hadv_flx6(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.7 )then

        call     hadv_flx7(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.8 )then

        call     hadv_flx8(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.9 )then

        call     hadv_flx9(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      elseif( hadvordrv.eq.10 )then

        call     hadv_flx10(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)

      else

        print *,' 98615 '
        call stopcm1

      endif

    endif

    if( do_ib .and. hadvordrv.ge.3 )then
      if( doweno )then
        call   ib_side_weno(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d,hadvordrv,bndy,hflxw,hflxe,hflxs,hflxn,weps)
      else
        call   ib_side_flx(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d,hadvordrv,bndy,hflxw,hflxe,hflxs,hflxn)
      endif
    endif

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow

      if(doitw)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=j1,j2
          if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
            dumx(i,j,k)=dumx(i+1,j,k)
          endif
        enddo
        enddo
      endif

      if(doite)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=j1,j2
          if((rru(i,j,k)+rru(i,j-1,k)).le.0.0)then
            dumx(i,j,k)=dumx(i-1,j,k)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    !  hadv tendencies:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        advy(i,j,k) = -(dumy(i,j,k)-dumy(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow
    !            (here, we subtract-off the divx piece)

      if(doitw)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=j1,j2
          if((rru(1,j,k)+rru(1,j-1,k)).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
        enddo
      endif

      if(doite)then
        i=ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=j1,j2
          if((rru(ni+1,j,k)+rru(ni+1,j-1,k)).le.0.0)then
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    ! 6th-order diffusion-v:
    if( idiff.ge.1 .and. difforder.eq.6 )then
      if( dovbud .and. nrk.eq.nrkmax .and. vd_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vdiag(i,j,k,vd_hediff) = (advx(i,j,k)+advy(i,j,k))
        enddo
        enddo
        enddo
      endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        dumx(i,j,k)=( 10.0*(v3d(i  ,j,k)-v3d(i-1,j,k))     &
                      -5.0*(v3d(i+1,j,k)-v3d(i-2,j,k))     &
                          +(v3d(i+2,j,k)-v3d(i-3,j,k)) )   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))     &
                          +(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni+1
          if( dumx(i,j,k)*(v3d(i,j,k)-v3d(i-1,j,k)).le.0.0 )then
            dumx(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj+1
      do i=1,ni
        advx(i,j,k)=advx(i,j,k)+coef*(dumx(i+1,j,k)-dumx(i,j,k))*(rho0(i,j-1,k)+rho0(i,j,k))/(rho(i,j-1,k)+rho(i,j,k))
      enddo
      enddo
      do j=1,nj+2
      do i=1,ni
        dumy(i,j,k)=( 10.0*(v3d(i,j  ,k)-v3d(i,j-1,k))     &
                      -5.0*(v3d(i,j+1,k)-v3d(i,j-2,k))     &
                          +(v3d(i,j+2,k)-v3d(i,j-3,k)) )*rho(i,j-1,k)
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+2
        do i=1,ni
          if( dumy(i,j,k)*(v3d(i,j,k)-v3d(i,j-1,k)).le.0.0 )then
            dumy(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj+1
      do i=1,ni
        advy(i,j,k)=advy(i,j,k)+coef*(dumy(i,j+1,k)-dumy(i,j,k))*(rho0(i,j-1,k)+rho0(i,j,k))/(rho(i,j-1,k)+rho(i,j,k))
      enddo
      enddo
    enddo
      if( dovbud .and. nrk.eq.nrkmax .and. vd_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vdiag(i,j,k,vd_hediff) = (advx(i,j,k)+advy(i,j,k))-vdiag(i,j,k,vd_hediff)
        enddo
        enddo
        enddo
      endif
    endif
    !-------------------------------------------------------


    endif  hadvsection

!----------------------------------------------------------------
! Advection in z-direction


    if( doweno )then

      if(     weno_order.eq.3 )then

        call     vadv_weno3(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d,0,weps)

      elseif( weno_order.eq.5 )then

        call     vadv_weno5(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d,0,0,weps)

      elseif( weno_order.eq.7 )then

        call     vadv_weno7(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d,0,weps)

      elseif( weno_order.eq.9 )then

        call     vadv_weno9(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d,0,weps)

      else

        print *,' 12946 '
        call stopcm1

      endif

    else

      if(     vadvordrv.eq.2 )then

        call     vadv_flx2(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.3 )then

        call     vadv_flx3(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.4 )then

        call     vadv_flx4(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.5 )then

        call     vadv_flx5(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.6 )then

        call     vadv_flx6(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.7 )then

        call     vadv_flx7(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.8 )then

        call     vadv_flx8(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.9 )then

        call     vadv_flx9(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      elseif( vadvordrv.eq.10 )then

        call     vadv_flx10(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)

      else

        print *,' 98616 '
        call stopcm1

      endif

    endif

    if( do_ib .and. vadvordrv.ge.3 )then
      if( doweno )then
        call   ib_lwr_weno( 3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d,kbdy,vadvordrv,weps)
      else
        call   ib_lwr_flx(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d,kbdy,vadvordrv)
      endif
    endif

!------

    if(terrain_flag)then

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do j=j1,j2
        do i=1,ni
          k=1
          advz(i,j,k) = -dumz(i,j,k+1)*rdsf(k)
          vten(i,j,k) = vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )         &
                  *gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
          k=nk
          advz(i,j,k) = +dumz(i,j,k)*rdsf(k)
          vten(i,j,k) = vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )         &
                  *gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk-1
        do j=j1,j2
        do i=1,ni
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdsf(k)
          vten(i,j,k) = vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )         &
                  *gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
        enddo
        enddo
        enddo

    else

    !--------
    if( axisymm.eq.0 )then
      ! Cartesian grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do j=j1,j2
        do i=1,ni
          k=1
          advz(i,j,k) = -dumz(i,j,k+1)*rdz*mh(1,1,k)
          k=nk
          advz(i,j,k) = +dumz(i,j,k)*rdz*mh(1,1,k)
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk-1
        do j=j1,j2
        do i=1,ni
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
        enddo
        enddo
        enddo

          !-----------------------------------------------------------------------------------------------------------
        if( do_ib .and. ( gk00fix.eq.1 .or. gk00fix.eq.3 ) )then
          ! apply fix from Gallus and Klemp (2000, MWR, pg 1160)
          do k=1,kmaxib
          do j=j1,j2
          do i=1,ni
            if( bndy(i,j,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i,j-1,k)) .and. (.not.bndy(i,j-1,k-1)) )then
              ! upper-left corner:
              ! Note: advective form is used here, so remove the piece of advz from divx
              advz(i,j,k) = 0.25*( -(rrw(i,j-1,k+1)+rrw(i,j,k+1))*(v3d(i,j,k+1)-v3d(i,j,k  ))*rdz*mf(1,1,k+1)    &
                                   -(rrw(i,j-1,k  )+rrw(i,j,k  ))*(w3d(i,j,k)-w3d(i,j-1,k))*rdy*vf(j) )  &
                            -v3d(i,j,k)*0.5*( (rrw(i,j  ,k+1)-rrw(i,j  ,k))*rdz*mh(1,1,k) &
                                             +(rrw(i,j-1,k+1)-rrw(i,j-1,k))*rdz*mh(1,1,k) )
            endif
            if( bndy(i,j-1,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i,j-1,k)) .and. (.not.bndy(i,j,k-1)) )then
              ! upper-right corner:
              ! Note: advective form is used here, so remove the piece of advz from divx
              advz(i,j,k) = 0.25*( -(rrw(i,j-1,k+1)+rrw(i,j,k+1))*(v3d(i,j,k+1)-v3d(i,j,k  ))*rdz*mf(1,1,k+1)    &
                                   -(rrw(i,j-1,k  )+rrw(i,j,k  ))*(w3d(i,j,k)-w3d(i,j-1,k))*rdy*vf(j) )  &
                            -v3d(i,j,k)*0.5*( (rrw(i,j  ,k+1)-rrw(i,j  ,k))*rdz*mh(1,1,k) &
                                             +(rrw(i,j-1,k+1)-rrw(i,j-1,k))*rdz*mh(1,1,k) )
            endif
          enddo
          enddo
          enddo
        endif
          !-----------------------------------------------------------------------------------------------------------
        if( do_ib .and. ( gk00fix.eq.2 .or. gk00fix.eq.3 ) )then
          ! apply fix from Gallus and Klemp (2000, MWR, pg 1160)
          do k=1,kmaxib
          do j=j1,j2
          do i=1,ni
            if( bndy(i+1,j,k) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i+1,j-1,k)) .and. (.not.bndy(i,j-1,k)) )then
              ! SW corner:  replace -u dv/dx with -u du/dy:
              advx(i,j,k) = 0.25*( -(rru(i+1,j,k)+rru(i+1,j-1,k))*(v3d(i+1,j,k)-v3d(i  ,j,k))*rdx*uf(i+1)     &
                                   -(rru(i  ,j,k)+rru(i  ,j-1,k))*(v3d(i  ,j,k)-v3d(i-1,j,k))*rdx*uf(i  )  )  &
                            -v3d(i,j,k)*0.5*( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i) &
                                             +(rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i) )
            endif
            if( (.not.bndy(i,j,k)) .and. bndy(i-1,j,k) .and. (.not.bndy(i,j-1,k)) .and. (.not.bndy(i-1,j-1,k)) )then
              ! SE corner:  replace -u dv/dx with -u du/dy:
              advx(i,j,k) = 0.25*( -(rru(i+1,j,k)+rru(i+1,j-1,k))*(v3d(i+1,j,k)-v3d(i  ,j,k))*rdx*uf(i+1)     &
                                   -(rru(i  ,j,k)+rru(i  ,j-1,k))*(v3d(i  ,j,k)-v3d(i-1,j,k))*rdx*uf(i  )  )  &
                            -v3d(i,j,k)*0.5*( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i) &
                                             +(rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i) )
            endif
            if( (.not.bndy(i+1,j,k)) .and. (.not.bndy(i,j,k)) .and. bndy(i+1,j-1,k) .and. (.not.bndy(i,j-1,k)) )then
              ! NW corner:  replace -u dv/dx with -u du/dy:
              advx(i,j,k) = 0.25*( -(rru(i+1,j,k)+rru(i+1,j-1,k))*(v3d(i+1,j,k)-v3d(i  ,j,k))*rdx*uf(i+1)     &
                                   -(rru(i  ,j,k)+rru(i  ,j-1,k))*(v3d(i  ,j,k)-v3d(i-1,j,k))*rdx*uf(i  )  )  &
                            -v3d(i,j,k)*0.5*( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i) &
                                             +(rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i) )
            endif
            if( (.not.bndy(i,j,k)) .and. (.not.bndy(i-1,j,k)) .and. (.not.bndy(i,j-1,k)) .and. bndy(i-1,j-1,k) )then
              ! NE corner:  replace -u dv/dx with -u du/dy:
              advx(i,j,k) = 0.25*( -(rru(i+1,j,k)+rru(i+1,j-1,k))*(v3d(i+1,j,k)-v3d(i  ,j,k))*rdx*uf(i+1)     &
                                   -(rru(i  ,j,k)+rru(i  ,j-1,k))*(v3d(i  ,j,k)-v3d(i-1,j,k))*rdx*uf(i  )  )  &
                            -v3d(i,j,k)*0.5*( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i) &
                                             +(rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i) )
            endif
          enddo
          enddo
          enddo
        endif
          !-----------------------------------------------------------------------------------------------------------

      ! cm1ib: moved vten

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          vten(i,j,k) = vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )*rr0(1,1,k)
        enddo
        enddo
        enddo

    !--------
    elseif( axisymm.eq.1 )then
      ! axisymmetric grid:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do j=j1,j2
        do i=1,ni
          k=1
          advz(i,j,k) = -dumz(i,j,k+1)*rdz*mh(1,1,k)
          vten(i,j,k) = vten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +v3d(i,j,k)*divx(i,j,k) )*rr0(1,1,k)
          k=nk
          advz(i,j,k) = +dumz(i,j,k)*rdz*mh(1,1,k)
          vten(i,j,k) = vten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +v3d(i,j,k)*divx(i,j,k) )*rr0(1,1,k)
        enddo
        enddo

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk-1
        do j=j1,j2
        do i=1,ni
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          vten(i,j,k) = vten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +v3d(i,j,k)*divx(i,j,k) )*rr0(1,1,k)
        enddo
        enddo
        enddo

    endif
        !--------

    endif

!----------------------------------------------------------------

    if( dolsw )then
      ! vertical advection from specified large-scale w profile:

      call     wsub(ni  ,nj+1,nk  ,v3d,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)

    endif

!----------------------------------------------------------------
!  Diagnostics:

  doingvbud:  &
  if( dovbud .and. nrk.eq.nrkmax )then
  if( vd_hadv.ge.1 .and. vd_vadv.ge.1 )then
    !--------
    ! advective tendencies:
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          vdiag(i,j,k,vd_hadv) = ( advx(i,j,k)+advy(i,j,k)                        &
                + v3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i,j  ,k))*rdy*vh(j  ) )  &
                                  +( (rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j  ,k)-rrv(i,j-1,k))*rdy*vh(j-1) )  & 
                                 ) )*rr0(1,1,k)
          vdiag(i,j,k,vd_vadv) = ( advz(i,j,k)                                  &
                + v3d(i,j,k)*0.5*( (rrw(i,j  ,k+1)-rrw(i,j  ,k))*rdz*mh(1,1,k)  &
                                  +(rrw(i,j-1,k+1)-rrw(i,j-1,k))*rdz*mh(1,1,k)  &
                                 ) )*rr0(1,1,k)
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k)
        do k=1,nk
        do i=1,ni
          vdiag(i,j,k,vd_hadv) = ( advx(i,j,k)                         &
                + v3d(i,j,k)*( (arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i) )  )*rr0(1,1,k)
          vdiag(i,j,k,vd_vadv) = ( advz(i,j,k)                         &
                + v3d(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) )  )*rr0(1,1,k)
          vdiag(i,2,k,vd_hadv) = vdiag(i,1,k,vd_hadv)
          vdiag(i,2,k,vd_vadv) = vdiag(i,1,k,vd_vadv)
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          vdiag(i,j,k,vd_hadv) = ( advx(i,j,k)+advy(i,j,k)                        &
                + v3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i,j  ,k))*rdy*vh(j  ) )  &
                                  +( (rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j  ,k)-rrv(i,j-1,k))*rdy*vh(j-1) )  & 
                                 ) )*gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
          vdiag(i,j,k,vd_vadv) = ( advz(i,j,k)                                  &
                + v3d(i,j,k)*0.5*( (rrw(i,j  ,k+1)-rrw(i,j  ,k))*rdsf(k)        &
                                  +(rrw(i,j-1,k+1)-rrw(i,j-1,k))*rdsf(k)        &
                                 ) )*gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
        enddo
        enddo
        enddo
    endif
  endif
    !--------
    if( vd_lsw.ge.1 )then
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vdiag(i,j,k,vd_lsw) = subs(i,j,k)
      enddo
      enddo
      enddo
    endif
    !--------
    if( idiff.ge.1 .and. difforder.eq.6 )then
    if( vd_hediff.ge.1 .and. vd_hadv.ge.1 )then
      ! subtract-off diffusion from advection:
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vdiag(i,j,k,vd_hadv) = vdiag(i,j,k,vd_hadv)-vdiag(i,j,k,vd_hediff)
      enddo
      enddo
      enddo
    endif
    endif
    !--------
    gethidiffv:  &
    if( vd_hidiff.ge.1 .and. vd_hadv.ge.1 )then
      ! horiz implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh31a !
    if( axisymm.eq.0 )then
      if(     hadvordrv.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        call     hadv_flx4(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)
      elseif( hadvordrv.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        call     hadv_flx6(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)
      elseif( hadvordrv.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        call     hadv_flx8(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)
      elseif( hadvordrv.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        call     hadv_flx10(  3 ,ni,nj+1,nk,c1,c2,rru,rrv,dumx,dumy,v3d)
      else
        print *,' 13961 '
        call stopcm1
      endif
    else
      if(     hadvordrv.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        hadv = 4
      elseif( hadvordrv.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        hadv = 6
      elseif( hadvordrv.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        hadv = 8
      elseif( hadvordrv.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        hadv = 10
      else
        print *,' 13981 '
        call stopcm1
      endif
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,1
      do i=0,ni+1
        ! advz stores M
        advz(i,j,k) = xh(i)*( v3d(i,j,k) + 0.5*fcor*xh(i) )
      enddo
      enddo
      enddo
      call advvaxi(.false.,xh,rxh,arh1,arh2,uh,xf,vf,rho0,rr0,rf0,rrf0,dumx,advx,advz,rru,hadv)
    endif
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh41 !
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          tot = vdiag(i,j,k,vd_hadv)
          advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          advy(i,j,k) = -(dumy(i,j,k)-dumy(i,j-1,k))*rdy*vf(j)
          ndiff                = ( advx(i,j,k)+advy(i,j,k)                        &
                + v3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i,j  ,k))*rdy*vh(j  ) )  &
                                  +( (rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j  ,k)-rrv(i,j-1,k))*rdy*vh(j-1) )  & 
                                 ) )*rr0(1,1,k)
          vdiag(i,j,k,vd_hadv) = ndiff
          vdiag(i,j,k,vd_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k,tot,ndiff)
        do k=1,nk
        do i=1,ni
          tot = vdiag(i,j,k,vd_hadv)
          ndiff                = ( advx(i,j,k)                         &
                + v3d(i,j,k)*( (arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i) )  )*rr0(1,1,k)
          vdiag(i,j,k,vd_hadv) = ndiff
          vdiag(i,j,k,vd_hidiff) = tot-ndiff
          !----
          vdiag(i,2,k,vd_hadv) = vdiag(i,1,k,vd_hadv)
          vdiag(i,2,k,vd_hidiff) = vdiag(i,1,k,vd_hidiff)
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          tot = vdiag(i,j,k,vd_hadv)
          advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          advy(i,j,k) = -(dumy(i,j,k)-dumy(i,j-1,k))*rdy*vf(j)
          ndiff                = ( advx(i,j,k)+advy(i,j,k)                        &
                + v3d(i,j,k)*0.5*( ( (rru(i+1,j  ,k)-rru(i,j  ,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j+1,k)-rrv(i,j  ,k))*rdy*vh(j  ) )  &
                                  +( (rru(i+1,j-1,k)-rru(i,j-1,k))*rdx*uh(i  )    &
                                    +(rrv(i  ,j  ,k)-rrv(i,j-1,k))*rdy*vh(j-1) )  & 
                                 ) )*gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
          vdiag(i,j,k,vd_hadv) = ndiff
          vdiag(i,j,k,vd_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  gethidiffv
    !--------
    getvidiffv:  &
    if( vd_vidiff.ge.1 .and. vd_vadv.ge.1 )then
      ! vert implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh32a !
      if(     vadvordrv.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        call     vadv_flx4(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)
      elseif( vadvordrv.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        call     vadv_flx6(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)
      elseif( vadvordrv.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        call     vadv_flx8(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)
      elseif( vadvordrv.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        call     vadv_flx10(  3 ,ni,nj+1,nk,c1,c2,rrw,dumz,v3d)
      else
        print *,' 13962 '
        call stopcm1
      endif
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=j1,j2
      do i=1,ni
        dumz(i,j,1) = 0.0
        dumz(i,j,nk+1) = 0.0
      enddo
      enddo
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh42a !
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          tot = vdiag(i,j,k,vd_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          ndiff                = ( advz(i,j,k)                                  &
                + v3d(i,j,k)*0.5*( (rrw(i,j  ,k+1)-rrw(i,j  ,k))*rdz*mh(1,1,k)  &
                                  +(rrw(i,j-1,k+1)-rrw(i,j-1,k))*rdz*mh(1,1,k)  &
                                 ) )*rr0(1,1,k)
          vdiag(i,j,k,vd_vadv) = ndiff
          vdiag(i,j,k,vd_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k,tot,ndiff)
        do k=1,nk
        do i=1,ni
          tot = vdiag(i,j,k,vd_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)
          ndiff                = ( advz(i,j,k)                         &
                + v3d(i,j,k)*( (rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) )  )*rr0(1,1,k)
          vdiag(i,j,k,vd_vadv) = ndiff
          vdiag(i,j,k,vd_vidiff) = tot-ndiff
          vdiag(i,2,k,vd_vadv) = vdiag(i,1,k,vd_vadv)
          vdiag(i,2,k,vd_vidiff) = vdiag(i,1,k,vd_vidiff)
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=1,nk
        do j=j1,j2
        do i=1,ni
          tot = vdiag(i,j,k,vd_vadv)
          advz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdsf(k)
          ndiff                = ( advz(i,j,k)                                  &
                + v3d(i,j,k)*0.5*( (rrw(i,j  ,k+1)-rrw(i,j  ,k))*rdsf(k)        &
                                  +(rrw(i,j-1,k+1)-rrw(i,j-1,k))*rdsf(k)        &
                                 ) )*gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
          vdiag(i,j,k,vd_vadv) = ndiff
          vdiag(i,j,k,vd_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  getvidiffv
    !--------
  endif  doingvbud

!----------------------------------------------------------------

    if( dolsw )then
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        vten(i,j,k) = vten(i,j,k)+subs(i,j,k)
      enddo
      enddo
      enddo
    endif

!----------------------------------------------------------------

      if(timestats.ge.1) time_advv=time_advv+mytime()
 
      end subroutine advv


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advw(nrk,xh,rxh,arh1,arh2,uh,xf,vh,gz,rgz,mh,mf,rho0,rr0,rf0,rrf0,  &
                      dumx,dumy,dumz,advx,advy,advz,subs,divx,                       &
                      rru,rrv,rrw,w3d  ,wten,rds,rdsf,c1,c2,rho,dt,                  &
                      dowbud ,wdiag,hadvorder,vadvorder,advweno ,bndy,kbdy,uf,vf,u3d,v3d,hflxw,hflxe,hflxs,hflxn,iamw)
      use input
      use constants
      use adv_routines
      use ib_module
      implicit none

      integer, intent(in) :: nrk
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy,dumz,advx,advy,advz,subs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: divx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten
      real, intent(in), dimension(kb:ke) :: rds
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(in) :: dt
      logical, intent(in) :: dowbud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
      integer, intent(in) :: hadvorder,vadvorder,advweno
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
      logical, intent(in) :: iamw
 
      integer :: i,j,k,hadv
      real :: ubar,vbar,cc1,cc2
      logical :: doitw,doite,doits,doitn
      logical :: doweno
      double precision :: weps
      real :: coef,tot,ndiff

!----------------------------------------------------------------

      doweno = .false.
      if( (advweno.eq.1) .or. (advweno.eq.2.and.nrk.eq.nrkmax) ) doweno = .true.
      weps = 100.0*epsilon

      if( idiff.ge.1 .and. difforder.eq.6 )then
        coef = kdiff6/64.0/dt
      endif

!-----------------

      doitw = wbc.eq.2 .and. ibw.eq.1
      doite = ebc.eq.2 .and. ibe.eq.1

!-----------------

      doits = sbc.eq.2 .and. ibs.eq.1
      doitn = nbc.eq.2 .and. ibn.eq.1

!----------------------------------------------------------------

    hadvsection:  if(axisymm.eq.1)then
      call advwaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vh,rho0,rr0,rf0,rrf0,dumx,advx,rru,w3d,c1,c2,hadvorder)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    else

!----------------------------------------------------------------
! Advection in horizontal directions

    if( doweno )then

      if(     weno_order.eq.3 )then

        call     hadv_weno3(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d,0,weps)

      elseif( weno_order.eq.5 )then

        call     hadv_weno5(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d,0,weps)

      elseif( weno_order.eq.7 )then

        call     hadv_weno7(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d,0,weps)

      elseif( weno_order.eq.9 )then

        call     hadv_weno9(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d,0,weps)

      else

        print *,' 12947 '
        call stopcm1

      endif

    else

      if(     hadvorder.eq.2 )then

        call     hadv_flx2(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.3 )then

        call     hadv_flx3(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.4 )then

        call     hadv_flx4(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.5 )then

        call     hadv_flx5(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.6 )then

        call     hadv_flx6(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.7 )then

        call     hadv_flx7(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.8 )then

        call     hadv_flx8(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.9 )then

        call     hadv_flx9(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      elseif( hadvorder.eq.10 )then

        call     hadv_flx10(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)

      else

        print *,' 98617 '
        call stopcm1

      endif

    endif

    if( do_ib .and. hadvorder.ge.3 )then
      if( doweno )then
        call   ib_side_weno(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d,hadvorder,bndy,hflxw,hflxe,hflxs,hflxn,weps)
      else
        call   ib_side_flx(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d,hadvorder,bndy,hflxw,hflxe,hflxs,hflxn)
      endif
    endif

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow

      if(doitw)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do j=1,nj
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k)=dumx(i+1,j,k)
          endif
        enddo
        enddo
      endif

      if(doite)then
        i=ni+1
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do j=1,nj
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.le.0.0)then
           dumx(i,j,k)=dumx(i-1,j,k)
          endif
        enddo
        enddo
      endif

      if(doits)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do i=1,ni
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k)=dumy(i,j+1,k)
          endif
        enddo
        enddo
      endif

      if(doitn)then
        j=nj+1
!$omp parallel do default(shared)   &
!$omp private(i,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do i=1,ni
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.le.0.0)then
            dumy(i,j,k)=dumy(i,j-1,k)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    !  hadv tendencies:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
      enddo
      enddo
      enddo

    !-------------------------------------------------------
    !  open bc:  set hadv to zero at inflow
    !            (here, we subtract-off the divx piece)

      if(doitw)then
        i=1
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do j=1,nj
          cc2 = 0.5*(c2(0,j,k)+c2(1,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(1,j,k)+cc1*rru(1,j,k-1)
          if(ubar.ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*(                    &
                    c1(i,j,k)*(rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +c2(i,j,k)*(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
        enddo
      endif

      if(doite)then
        i=ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do j=1,nj
          cc2 = 0.5*(c2(ni,j,k)+c2(ni+1,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(ni+1,j,k)+cc1*rru(ni+1,j,k-1)
          if(ubar.le.0.0)then
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*(                    &
                    c1(i,j,k)*(rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +c2(i,j,k)*(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
        enddo
      endif

      if(doits)then
        j=1
!$omp parallel do default(shared)   &
!$omp private(i,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do i=1,ni
          cc2 = 0.5*(c2(i,0,k)+c2(i,1,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,1,k)+cc1*rrv(i,1,k-1)
          if(vbar.ge.0.0)then
            advy(i,j,k)=advy(i,j,k)-w3d(i,j,k)*(                     &
                           c1(i,j,k)*(rrv(i,j+1,k-1)-rrv(i,j,k-1))   &
                          +c2(i,j,k)*(rrv(i,j+1,k  )-rrv(i,j,k  )) )*rdy*vh(j)
          endif
        enddo
        enddo
      endif

      if(doitn)then
        j=nj
!$omp parallel do default(shared)   &
!$omp private(i,k,ubar,vbar,cc1,cc2)
        do k=2,nk
        do i=1,ni
          cc2 = 0.5*(c2(i,nj,k)+c2(i,nj+1,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,nj+1,k)+cc1*rrv(i,nj+1,k-1)
          if(vbar.le.0.0)then
            advy(i,j,k)=advy(i,j,k)-w3d(i,j,k)*(                     &
                           c1(i,j,k)*(rrv(i,j+1,k-1)-rrv(i,j,k-1))   &
                          +c2(i,j,k)*(rrv(i,j+1,k  )-rrv(i,j,k  )) )*rdy*vh(j)
          endif
        enddo
        enddo
      endif

    !-------------------------------------------------------
    ! 6th-order diffusion-w:
    if( idiff.ge.1 .and. difforder.eq.6 )then
      if( dowbud .and. nrk.eq.nrkmax .and. wd_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_hediff) = (advx(i,j,k)+advy(i,j,k))
        enddo
        enddo
        enddo
      endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k,cc1,cc2)
    do k=2,nk
      do j=1,nj
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumx(i,j,k)=( 10.0*(w3d(i  ,j,k)-w3d(i-1,j,k))     &
                      -5.0*(w3d(i+1,j,k)-w3d(i-2,j,k))     &
                          +(w3d(i+2,j,k)-w3d(i-3,j,k)) )   &
              *0.5*( cc2*(rho(i-1,j,k  )+rho(i,j,k  ))    &
                    +cc1*(rho(i-1,j,k-1)+rho(i,j,k-1)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj
        do i=1,ni+1
          if( dumx(i,j,k)*(w3d(i,j,k)-w3d(i-1,j,k)).le.0.0 )then
            dumx(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=advx(i,j,k)+coef*(dumx(i+1,j,k)-dumx(i,j,k))*rf0(i,j,k)/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
      enddo
      enddo
      do j=1,nj+1
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumy(i,j,k)=( 10.0*(w3d(i,j  ,k)-w3d(i,j-1,k))     &
                     -5.0*(w3d(i,j+1,k)-w3d(i,j-2,k))     &
                         +(w3d(i,j+2,k)-w3d(i,j-3,k)) )   &
              *0.5*( cc2*(rho(i,j-1,k  )+rho(i,j,k  ))    &
                    +cc1*(rho(i,j-1,k-1)+rho(i,j,k-1)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni
          if( dumy(i,j,k)*(w3d(i,j,k)-w3d(i,j-1,k)).le.0.0 )then
            dumy(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=advy(i,j,k)+coef*(dumy(i,j+1,k)-dumy(i,j,k))*rf0(i,j,k)/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
      enddo
      enddo
    enddo
      if( dowbud .and. nrk.eq.nrkmax .and. wd_hediff.ge.1 )then
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_hediff) = (advx(i,j,k)+advy(i,j,k))-wdiag(i,j,k,wd_hediff)
        enddo
        enddo
        enddo
      endif
    endif
    !-------------------------------------------------------


    endif  hadvsection

!----------------------------------------------------------------
! Advection in z-direction


    if( doweno )then

      if(     weno_order.eq.3 )then

        call     vadv_weno3(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d,0,weps)

      elseif( weno_order.eq.5 )then

        call     vadv_weno5(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d,0,0,weps)

      elseif( weno_order.eq.7 )then

        call     vadv_weno7(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d,0,weps)

      elseif( weno_order.eq.9 )then

        call     vadv_weno9(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d,0,weps)

      else

        print *,' 12948 '
        call stopcm1

      endif

    else

      if(     vadvorder.eq.2 )then

        call     vadv_flx2(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.3 )then

        call     vadv_flx3(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.4 )then

        call     vadv_flx4(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.5 )then

        call     vadv_flx5(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.6 )then

        call     vadv_flx6(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.7 )then

        call     vadv_flx7(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.8 )then

        call     vadv_flx8(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.9 )then

        call     vadv_flx9(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      elseif( vadvorder.eq.10 )then

        call     vadv_flx10(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)

      else

        print *,' 98618 '
        call stopcm1

      endif

    endif

    if( do_ib .and. vadvorder.ge.3 )then
      if( doweno )then
        call   ib_lwr_weno( 4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d,kbdy,vadvorder,weps)
      else
        call   ib_lwr_flx(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d,kbdy,vadvorder)
      endif
    endif

!------

      if(terrain_flag)then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advz(i,j,k) = -(dumz(i,j,k)-dumz(i,j,k-1))*rds(k)
        wten(i,j,k) = wten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                   +w3d(i,j,k)*(c2(i,j,k)*divx(i,j,k)+c1(i,j,k)*divx(i,j,k-1)) )*rrf0(i,j,k)*gz(i,j)
      enddo
      enddo
      enddo

      else

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advz(i,j,k) = -(dumz(i,j,k)-dumz(i,j,k-1))*rdz*mf(1,1,k)
      enddo
      enddo
      enddo

      if( iamw )then
        if( do_ib .and. ( gk00fix.eq.1 .or. gk00fix.eq.3 ) )then
          ! apply fix from Gallus and Klemp (2000, MWR, pg 1160)
          do k=2,kmaxib
          do j=1,nj
          do i=1,ni
            if( bndy(i+1,j,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i,j,k-1)) .and. (.not.bndy(i+1,j,k)) )then
              ! upper-left corner, x dir:
              ! Note: advective form is used here, so remove the piece of advx from divx
              advx(i,j,k) = 0.25*( -(rru(i+1,j,k)+rru(i+1,j,k-1))*(u3d(i+1,j,k)-u3d(i+1,j,k-1))*rdz*mf(1,1,k)   &
                                   -(rru(i  ,j,k)+rru(i  ,j,k-1))*(w3d(i  ,j,k)-w3d(i-1,j,k))*rdx*uf(i)   )  &
                            -w3d(i,j,k)*( c2(1,1,k)*(rru(i+1,j,k  )-rru(i,j,k  ))*rdx*uh(i)    &
                                         +c1(1,1,k)*(rru(i+1,j,k-1)-rru(i,j,k-1))*rdx*uh(i)  )
            endif
            if( bndy(i-1,j,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i,j,k-1)) .and. (.not.bndy(i-1,j,k)) )then
              ! upper-right corner, x dir:
              ! Note: advective form is used here, so remove the piece of advx from divx
              advx(i,j,k) = 0.25*( -(rru(i+1,j,k)+rru(i+1,j,k-1))*(w3d(i+1,j,k)-w3d(i  ,j,k))*rdx*uf(i+1)    &
                                   -(rru(i  ,j,k)+rru(i  ,j,k-1))*(u3d(i,j,k)-u3d(i,j,k-1))*rdz*mf(1,1,k)  )  &
                            -w3d(i,j,k)*( c2(1,1,k)*(rru(i+1,j,k  )-rru(i,j,k  ))*rdx*uh(i)    &
                                         +c1(1,1,k)*(rru(i+1,j,k-1)-rru(i,j,k-1))*rdx*uh(i)  )
            endif
            if( bndy(i,j+1,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i,j,k-1)) .and. (.not.bndy(i,j+1,k)) )then
              ! upper-left corner, y dir:
              ! Note: advective form is used here, so remove the piece of advx from divx
              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i,j+1,k-1))*(v3d(i,j+1,k)-v3d(i,j+1,k-1))*rdz*mf(1,1,k)   &
                                   -(rrv(i,j  ,k)+rrv(i,j  ,k-1))*(w3d(i  ,j,k)-w3d(i,j-1,k))*rdy*vf(j)   )  &
                            -w3d(i,j,k)*( c2(1,1,k)*(rrv(i,j+1,k  )-rrv(i,j,k  ))*rdy*vh(j)    &
                                         +c1(1,1,k)*(rrv(i,j+1,k-1)-rrv(i,j,k-1))*rdy*vh(j)  )
            endif
            if( bndy(i,j-1,k-1) .and. (.not.bndy(i,j,k)) .and. (.not.bndy(i,j,k-1)) .and. (.not.bndy(i,j-1,k)) )then
              ! upper-right corner, y dir:
              ! Note: advective form is used here, so remove the piece of advx from divx
              advy(i,j,k) = 0.25*( -(rrv(i,j+1,k)+rrv(i,j+1,k-1))*(w3d(i,j+1,k)-w3d(i,j  ,k))*rdy*vf(j+1)    &
                                   -(rrv(i,j  ,k)+rrv(i,j  ,k-1))*(v3d(i,j,k)-v3d(i,j,k-1))*rdz*mf(1,1,k)  )  &
                            -w3d(i,j,k)*( c2(1,1,k)*(rrv(i,j+1,k  )-rrv(i,j,k  ))*rdy*vh(j)    &
                                         +c1(1,1,k)*(rrv(i,j+1,k-1)-rrv(i,j,k-1))*rdy*vh(j)  )
            endif
          enddo
          enddo
          enddo
        endif
      endif

      ! cm1ib: moved wten

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k) = wten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                   +w3d(i,j,k)*(c2(1,1,k)*divx(i,j,k)+c1(1,1,k)*divx(i,j,k-1)) )*rrf0(1,1,k)
      enddo
      enddo
      enddo

      endif

!----------------------------------------------------------------
!  Diagnostics:

  doingwbud:  &
  if( dowbud .and. nrk.eq.nrkmax )then
  if( wd_hadv.ge.1 .and. wd_vadv.ge.1 )then
    !--------
    ! advective tendencies:
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_hadv) = ( advx(i,j,k)+advy(i,j,k)                            &
                + w3d(i,j,k)*( c2(i,j,k)*( (rru(i+1,j,k  )-rru(i,j,k  ))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k  )-rrv(i,j,k  ))*rdy*vh(j) )  &
                              +c1(i,j,k)*( (rru(i+1,j,k-1)-rru(i,j,k-1))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k-1)-rrv(i,j,k-1))*rdy*vh(j) )  &
                             ) )*rrf0(1,1,k)
          wdiag(i,j,k,wd_vadv) = ( advz(i,j,k)                                        &
                + w3d(i,j,k)*( c2(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k  ))*rdz*mh(1,1,k  )  &
                              +c1(i,j,k)*(rrw(i,j,k  )-rrw(i,j,k-1))*rdz*mh(1,1,k-1)  &
                             ) )*rrf0(1,1,k)
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k)
        do k=2,nk
        do i=1,ni
          wdiag(i,j,k,wd_hadv) = ( advx(i,j,k)                         &
                + w3d(i,j,k)*( c2(1,1,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)  &
                              +c1(1,1,k)*(arh2(i)*rru(i+1,j,k-1)-arh1(i)*rru(i,j,k-1))*rdx*uh(i)  &
                             )   )*rrf0(1,1,k)
          wdiag(i,j,k,wd_vadv) = ( advz(i,j,k)                         &
                + w3d(i,j,k)*( c2(1,1,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)  &
                              +c1(1,1,k)*(rrw(i,j,k)-rrw(i,j,k-1))*rdz*mh(1,1,k-1)  &
                             )   )*rrf0(1,1,k)
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_hadv) = ( advx(i,j,k)+advy(i,j,k)                            &
                + w3d(i,j,k)*( c2(i,j,k)*( (rru(i+1,j,k  )-rru(i,j,k  ))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k  )-rrv(i,j,k  ))*rdy*vh(j) )  &
                              +c1(i,j,k)*( (rru(i+1,j,k-1)-rru(i,j,k-1))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k-1)-rrv(i,j,k-1))*rdy*vh(j) )  &
                             ) )*rrf0(i,j,k)*gz(i,j)
          wdiag(i,j,k,wd_vadv) = ( advz(i,j,k)                                        &
                + w3d(i,j,k)*( c2(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k  ))*rdsf(k)  &
                              +c1(i,j,k)*(rrw(i,j,k  )-rrw(i,j,k-1))*rdsf(k)  &
                             ) )*rrf0(i,j,k)*gz(i,j)
        enddo
        enddo
        enddo
    endif
  endif
    !--------
    if( idiff.ge.1 .and. difforder.eq.6 )then
    if( wd_hediff.ge.1 .and. wd_hadv.ge.1 )then
      ! subtract-off diffusion from advection:
      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wdiag(i,j,k,wd_hadv) = wdiag(i,j,k,wd_hadv)-wdiag(i,j,k,wd_hediff)
      enddo
      enddo
      enddo
    endif
    endif
    !--------
    gethidiffw:  &
    if( wd_hidiff.ge.1 .and. wd_hadv.ge.1 )then
      ! horiz implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh31 !
    if( axisymm.eq.0 )then
      if(     hadvorder.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        call     hadv_flx4(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)
      elseif( hadvorder.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        call     hadv_flx6(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)
      elseif( hadvorder.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        call     hadv_flx8(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)
      elseif( hadvorder.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        call     hadv_flx10(  4 ,ni,nj,nk+1,c1,c2,rru,rrv,dumx,dumy,w3d)
      else
        print *,' 13971 '
        call stopcm1
      endif
    else
      if(     hadvorder.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        hadv = 4
      elseif( hadvorder.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        hadv = 6
      elseif( hadvorder.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        hadv = 8
      elseif( hadvorder.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        hadv = 10
      else
        print *,' 13981 '
        call stopcm1
      endif
      call advwaxi(.false.,xh,rxh,arh1,arh2,uh,xf,vh,rho0,rr0,rf0,rrf0,dumx,advx,rru,w3d,c1,c2,hadv)
    endif
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh41 !
    if(.not.terrain_flag)then
      if( axisymm.eq.0 )then
        ! Cartesian grid:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          tot = wdiag(i,j,k,wd_hadv)
          advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
          ndiff                = ( advx(i,j,k)+advy(i,j,k)                            &
                + w3d(i,j,k)*( c2(i,j,k)*( (rru(i+1,j,k  )-rru(i,j,k  ))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k  )-rrv(i,j,k  ))*rdy*vh(j) )  &
                              +c1(i,j,k)*( (rru(i+1,j,k-1)-rru(i,j,k-1))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k-1)-rrv(i,j,k-1))*rdy*vh(j) )  &
                             ) )*rrf0(1,1,k)
          wdiag(i,j,k,wd_hadv) = ndiff
          wdiag(i,j,k,wd_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
      else
        ! axisymmetric grid:
        j=1
        !$omp parallel do default(shared)   &
        !$omp private(i,k,tot,ndiff)
        do k=2,nk
        do i=1,ni
          tot = wdiag(i,j,k,wd_hadv)
          ndiff                = ( advx(i,j,k)                         &
                + w3d(i,j,k)*( c2(1,1,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)  &
                              +c1(1,1,k)*(arh2(i)*rru(i+1,j,k-1)-arh1(i)*rru(i,j,k-1))*rdx*uh(i)  &
                             )   )*rrf0(1,1,k)
          wdiag(i,j,k,wd_hadv) = ndiff
          wdiag(i,j,k,wd_hidiff) = tot-ndiff
        enddo
        enddo
      endif
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          tot = wdiag(i,j,k,wd_hadv)
          advx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          advy(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
          ndiff                = ( advx(i,j,k)+advy(i,j,k)                            &
                + w3d(i,j,k)*( c2(i,j,k)*( (rru(i+1,j,k  )-rru(i,j,k  ))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k  )-rrv(i,j,k  ))*rdy*vh(j) )  &
                              +c1(i,j,k)*( (rru(i+1,j,k-1)-rru(i,j,k-1))*rdx*uh(i)    &
                                          +(rrv(i,j+1,k-1)-rrv(i,j,k-1))*rdy*vh(j) )  &
                             ) )*rrf0(i,j,k)*gz(i,j)
          wdiag(i,j,k,wd_hadv) = ndiff
          wdiag(i,j,k,wd_hidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  gethidiffw
    !--------
    getvidiffw:  &
    if( wd_vidiff.ge.1 .and. wd_vadv.ge.1 )then
      ! vert implicit diffusion tendency:
      ! step1: get non-diffusive advective fluxes:
      ! buh32 !
      if(     vadvorder.eq.3 .or. ( advwenov.ge.1 .and. weno_order.eq.3 ) )then
        call     vadv_flx4(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)
      elseif( vadvorder.eq.5 .or. ( advwenov.ge.1 .and. weno_order.eq.5 ) )then
        call     vadv_flx6(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)
      elseif( vadvorder.eq.7 .or. ( advwenov.ge.1 .and. weno_order.eq.7 ) )then
        call     vadv_flx8(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)
      elseif( vadvorder.eq.9 .or. ( advwenov.ge.1 .and. weno_order.eq.9 ) )then
        call     vadv_flx10(  4 ,ni,nj,nk+1,c1,c2,rrw,dumz,w3d)
      else
        print *,' 13972 '
        call stopcm1
      endif
      ! step2: get non-diffusive and diffusive components of advection:
      ! buh42a !
    if(.not.terrain_flag)then
        ! without terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          tot = wdiag(i,j,k,wd_vadv)
          advz(i,j,k) = -(dumz(i,j,k)-dumz(i,j,k-1))*rdz*mf(1,1,k)
          ndiff                = ( advz(i,j,k)                                        &
                + w3d(i,j,k)*( c2(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k  ))*rdz*mh(1,1,k  )  &
                              +c1(i,j,k)*(rrw(i,j,k  )-rrw(i,j,k-1))*rdz*mh(1,1,k-1)  &
                             ) )*rrf0(1,1,k)
          wdiag(i,j,k,wd_vadv) = ndiff
          wdiag(i,j,k,wd_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
    else
        ! Cartesian with terrain:
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tot,ndiff)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          tot = wdiag(i,j,k,wd_vadv)
          advz(i,j,k) = -(dumz(i,j,k)-dumz(i,j,k-1))*rds(k)
          ndiff                = ( advz(i,j,k)                                        &
                + w3d(i,j,k)*( c2(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k  ))*rdsf(k)  &
                              +c1(i,j,k)*(rrw(i,j,k  )-rrw(i,j,k-1))*rdsf(k)  &
                             ) )*rrf0(i,j,k)*gz(i,j)
          wdiag(i,j,k,wd_vadv) = ndiff
          wdiag(i,j,k,wd_vidiff) = tot-ndiff
        enddo
        enddo
        enddo
    endif
    endif  getvidiffw
    !--------
  endif  doingwbud

!----------------------------------------------------------------

      if(timestats.ge.1) time_advw=time_advw+mytime()

      end subroutine advw


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  end module adv_module
