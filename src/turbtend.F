  MODULE turbtend_module

  implicit none

  private
  public :: turbs,turbt,turbu,turbv,turbw,turbsz,turbuz,turbvz

  real :: fac,fac2

  CONTAINS


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbs(iflux,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,sflux,  &
                       rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                       turbx,turby,turbz,dumx,dumy,dumz,rho,rr,rf,s,sten,khh,khv,cm0,dum7,dum8, &
                       dobud,ibd,ied,jbd,jed,kbd,ked,ndiag,diag,sd_hturb,sd_vturb,ivar)
      use input
      use constants
      use cm1libs , only : rslf,rsif
      implicit none

      integer, intent(in) :: iflux
      real, intent(in) :: dt
      logical, intent(in) :: dosfcflx
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je) :: sflux
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dumx,dumy,dumz,sten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf,s
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: khh,khv
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: dum7,dum8
      logical, intent(in) :: dobud
      integer, intent(in) :: ibd,ied,jbd,jed,kbd,ked,ndiag,sd_hturb,sd_vturb
      real, intent(inout) , dimension(ibd:ied,jbd:jed,kbd:ked,ndiag) :: diag
      integer, intent(in) :: ivar    ! 1 = potential temperature
                                     ! 2 = water vapor mixing ratio

      integer :: i,j,k
      real :: r1,r2

!---------------------------------------------------------------

  dohoriz:  &
  IF( dohturb )THEN

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        do j=1,nj+1
        do i=1,ni+1
          !  x-direction
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                             *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
          !  y-direction
          dumy(i,j,k)= -0.125*( rho(i,j,k)+rho(i,j-1,k) )           &
                             *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                                +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )   &
                            *(    s(i,j,k)-   s(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            dumx(1,j,k) = dumx(2,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            dumx(ni+1,j,k) = dumx(ni,j,k)
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            dumy(i,1,k) = dumy(i,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            dumy(i,nj+1,k) = dumy(i,nj,k)
          enddo
        ENDIF

        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          turby(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE
      ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj

        do i=1,ni+1
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                             *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
        enddo

        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          dumx(ni+1,j,k) = arh1(ni)*dumx(ni,j,k)/arh2(ni)
        ENDIF

        !-----
        if(wbc.eq.3.or.wbc.eq.4)then
          ! assume zero flux:
          dumx(1,j,k) = 0.0
        endif
        if(ebc.eq.3.or.ebc.eq.4)then
          ! assume zero flux:
          dumx(ni+1,j,k) = 0.0
        endif
        !-----

        do i=1,ni
          turbx(i,j,k)=-(arh2(i)*dumx(i+1,j,k)-arh1(i)*dumx(i,j,k))*rdx*uh(i)
          turby(i,j,k)=0.0
        enddo

      enddo
      enddo

    ENDIF   ! endif for axisymm check

!---------------------------------------------------------------

  ELSE
      ! Cartesian with terrain:

      ! use turbz as a temporary array for s at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=0,nj+1

        ! lowest model level:
        do i=0,ni+1
          turbz(i,j,1) = cgs1*s(i,j,1)+cgs2*s(i,j,2)+cgs3*s(i,j,3)
        enddo

        ! upper-most model level:
        do i=0,ni+1
          turbz(i,j,nk+1) = cgt1*s(i,j,nk)+cgt2*s(i,j,nk-1)+cgt3*s(i,j,nk-2)
        enddo

        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+1
          turbz(i,j,k) = r1*s(i,j,k-1)+r2*s(i,j,k)
        enddo
        enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        ! x-flux
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.125*( gz(i,j)*rho(i,j,k)+gz(i-1,j)*rho(i-1,j,k) )  &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )*( &
                  (s(i,j,k)*rgz(i,j)-s(i-1,j,k)*rgz(i-1,j))         &
                   *rdx*uf(i)                                       &
              +0.5*( gxu(i,j,k+1)*(turbz(i,j,k+1)+turbz(i-1,j,k+1)) &
                    -gxu(i,j,k  )*(turbz(i,j,k  )+turbz(i-1,j,k  )) &
                   )*rdsf(k)*rgzu(i,j) )
        enddo
        enddo

        ! y-flux
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.125*( gz(i,j)*rho(i,j,k)+gz(i,j-1)*rho(i,j-1,k) )  &
                             *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                                +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )*( &
                  (s(i,j,k)*rgz(i,j)-s(i,j-1,k)*rgz(i,j-1))         &
                   *rdy*vf(j)                                       &
              +0.5*( gyv(i,j,k+1)*(turbz(i,j,k+1)+turbz(i,j-1,k+1)) &
                    -gyv(i,j,k  )*(turbz(i,j,k  )+turbz(i,j-1,k  )) &
                   )*rdsf(k)*rgzv(i,j) )
        enddo
        enddo

      enddo

      ! use turbz,dumz as temporary arrays for fluxes at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=1,nj+1
        ! lowest model level:
        do i=1,ni+1
          turbz(i,j,1) = cgs1*dumx(i,j,1)+cgs2*dumx(i,j,2)+cgs3*dumx(i,j,3)
           dumz(i,j,1) = cgs1*dumy(i,j,1)+cgs2*dumy(i,j,2)+cgs3*dumy(i,j,3)
        enddo

        ! upper-most model level:
        do i=1,ni+1
          turbz(i,j,nk+1) = cgt1*dumx(i,j,nk)+cgt2*dumx(i,j,nk-1)+cgt3*dumx(i,j,nk-2)
           dumz(i,j,nk+1) = cgt1*dumy(i,j,nk)+cgt2*dumy(i,j,nk-1)+cgt3*dumy(i,j,nk-2)
        enddo

        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=1,ni+1
          turbz(i,j,k) = r1*dumx(i,j,k-1)+r2*dumx(i,j,k)
           dumz(i,j,k) = r1*dumy(i,j,k-1)+r2*dumy(i,j,k)
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        ! x-tendency
        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -gz(i,j)*( dumx(i+1,j,k)*rgzu(i+1,j)             &
                                   -dumx(i  ,j,k)*rgzu(i  ,j) )*rdx*uh(i) &
                -( ( gx(i,j,k+1)*(turbz(i,j,k+1)+turbz(i+1,j,k+1))        &
                    -gx(i,j,k  )*(turbz(i,j,k  )+turbz(i+1,j,k  )) )      &
                 )*0.5*rdsf(k)
        enddo
        enddo

        ! y-tendency
        do j=1,nj
        do i=1,ni
          turby(i,j,k) = -gz(i,j)*( dumy(i,j+1,k)*rgzv(i,j+1)             &
                                   -dumy(i,j  ,k)*rgzv(i,j  ) )*rdy*vh(j) &
                -( ( gy(i,j,k+1)*( dumz(i,j,k+1)+ dumz(i,j+1,k+1))        &
                    -gy(i,j,k  )*( dumz(i,j,k  )+ dumz(i,j+1,k  )) )      &
                 )*0.5*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

  ELSE  dohoriz

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=0.0
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dohoriz

!---------------------------------------------------------------------
!  z-direction

  dovert:  &
  IF( dovturb )THEN

      call       turbsz(iflux,dt,dosfcflx,sflux,mh,mf,  &
                       turbz,dumx,dumy,dumz,rho,rr,rf,s,khv,dum7,dum8,ivar)

  ELSE  dovert

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dovert


!---------------------------------------------------------------------
!  Tendencies:

    IF(axisymm.eq.0)THEN

      IF( cm1setup.eq.4 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          if( cm0(i,j) .le. cmemin )then
            ! zero-out turb tendencies outside of LES domain:
            turbx(i,j,k) = 0.0
            turby(i,j,k) = 0.0
            turbz(i,j,k) = 0.0
          endif

          sten(i,j,k)=sten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rr(i,j,k)
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni

          sten(i,j,k)=sten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rr(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

    ELSE

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))*rr(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------
!  Diagnostics:

      IF( dobud )THEN
      if( sd_hturb.ge.1 .and. sd_vturb.ge.1 )then
        if( axisymm.eq.0 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            diag(i,j,k,sd_hturb) = (turbx(i,j,k)+turby(i,j,k))*rr(i,j,k)
            diag(i,j,k,sd_vturb) = turbz(i,j,k)*rr(i,j,k)
          enddo
          enddo
          enddo
        else
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            diag(i,j,k,sd_hturb) = turbx(i,j,k)*rr(i,j,k)
            diag(i,j,k,sd_vturb) = turbz(i,j,k)*rr(i,j,k)
          enddo
          enddo
          enddo
        endif
      endif
      ENDIF

!---------------------------------------------------------------------

      if(timestats.ge.1) time_ttend=time_ttend+mytime()

      return
      end subroutine turbs


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbsz(iflux,dt,dosfcflx,sflux,mh,mf,  &
                       turbz,dum1,dum2,dumz,rho,rr,rf,s,khv,dum7,dum8,ivar)
      use input
      use constants
      use cm1libs , only : rslf,rsif
      implicit none

      integer, intent(in) :: iflux
      real, intent(in) :: dt
      logical, intent(in) :: dosfcflx
      real, intent(in), dimension(ib:ie,jb:je) :: sflux
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbz,dum1,dum2,dumz
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf,s
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: khv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: dum7,dum8
      integer, intent(in) :: ivar

      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,cfa,cfb,cfc,cfd,kappa,ptcb,ptct


    if( cm1setup.eq.3 )then
      kappa = viscosity/pr_num
      ! bottom/top boundary conditions:
      if( ivar.eq.1 )then
        ! potential temperature:
        ptcb = ptc_bot - th0r
        ptct = ptc_top - th0r
      elseif( ivar.eq.2 )then
        ! water vapor mixing ratio:
        ! (assume saturation)
        ptcb = rslf( base_pbot , ptc_bot*base_pibot )
        ptct = rslf( base_ptop , ptc_top*base_pitop )
      else
        ptcb = 0.0
        ptct = 0.0
      endif
    endif

  ifimpls:  &
  IF( doimpl.eq.0 )THEN
      ! explicit vertical turbulence:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO j=1,nj

        do k=2,nk
        do i=1,ni
          dumz(i,j,k)= -khv(i,j,k)*(s(i,j,k)-s(i,j,k-1))*rdz*mf(i,j,k)*rf(i,j,k)
        enddo
        enddo

        IF( cm1setup.eq.1 .or. cm1setup.eq.2 .or. cm1setup.eq.4 )THEN
          ! LES or mesoscale modeling:

          IF(bcturbs.eq.1)THEN
 
            do i=1,ni
              dumz(i,j,1)=0.0
              dumz(i,j,nk+1)=0.0
            enddo

          ELSEIF(bcturbs.eq.2)THEN
 
            do i=1,ni
              dumz(i,j,1)=dumz(i,j,2)
              dumz(i,j,nk+1)=dumz(i,j,nk)
            enddo

          ENDIF

          if(iflux.eq.1 .and. dosfcflx)then
            do i=1,ni
              dumz(i,j,1) = sflux(i,j)*rf(i,j,1)
            enddo
          endif

        ELSEIF( cm1setup.eq.3 )THEN

          if(bc_temp.eq.1)then
            ! specified theta at boundary

            do i=1,ni
              dumz(i,j,1) = -kappa*2.0*(s(i,j,1)-ptcb)*rdz*mf(i,j,1)*rf(i,j,1)
              dumz(i,j,nk+1) = -kappa*2.0*(ptct-s(i,j,nk))*rdz*mf(i,j,nk+1)*rf(i,j,nk+1)
            enddo

          elseif(bc_temp.eq.2)then
            ! specified flux at boundary

            do i=1,ni
              dumz(i,j,1) = kappa*ptc_bot*rf(i,j,1)
              dumz(i,j,nk+1) = kappa*ptc_top*rf(i,j,nk+1)
            enddo

          endif

        ELSE

          print *,'  21086 '
          call stopcm1

        ENDIF

        do k=1,nk
        do i=1,ni
          turbz(i,j,k) = -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(i,j,k)
        enddo
        enddo

      ENDDO

  ELSE  ifimpls

      ! implicit vertical turbulence:

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

      ! boundary conditions:
      IF( cm1setup.eq.1 .or. cm1setup.eq.2 .or. cm1setup.eq.4 )THEN
        ! LES or mesoscale modeling:
        IF(bcturbs.eq.1)THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dumz(i,j,1)=0.0
            dumz(i,j,nk+1)=0.0
          enddo
          enddo
        ELSEIF(bcturbs.eq.2)THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dumz(i,j,1) = -khv(i,j,2)*(s(i,j,2)-s(i,j,1))*rdz*mf(i,j,2)*rf(i,j,2)
            dumz(i,j,nk+1) = -khv(i,j,nk)*(s(i,j,nk)-s(i,j,nk-1))*rdz*mf(i,j,nk)*rf(i,j,nk)
          enddo
          enddo
        ENDIF
        if(iflux.eq.1 .and. dosfcflx)then
          ! surface heat/moisture flux:
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dumz(i,j,1)=sflux(i,j)*rf(i,j,1)
          enddo
          enddo
        endif
      ELSEIF( cm1setup.eq.3 )THEN
        ! DNS bc:
        if(bc_temp.eq.1)then
          ! specified theta at boundary
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dumz(i,j,1) = -kappa*2.0*(s(i,j,1)-ptcb)*rdz*mf(i,j,1)*rf(i,j,1)
            dumz(i,j,nk+1) = -kappa*2.0*(ptct-s(i,j,nk))*rdz*mf(i,j,nk+1)*rf(i,j,nk+1)
          enddo
          enddo
        elseif(bc_temp.eq.2)then
          ! specified flux at boundary
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dumz(i,j,1) = kappa*ptc_bot*rf(i,j,1)
            dumz(i,j,nk+1) = kappa*ptc_top*rf(i,j,nk+1)
          enddo
          enddo
        endif
      ENDIF

      k = 1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem)
      do j=1,nj
      DO i=1,ni
          r2 = dum8(i,j,k)
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = s(i,j,k) + temb*( -r2*s(i,j,k)+r2*s(i,j,k+1) )  &
                         + temc*dumz(i,j,1)*mh(i,j,1)*rr(i,j,1)
        tem = 1.0/cfb
        dum1(i,j,1)=-cfc*tem
        dum2(i,j,1)= cfd*tem
      ENDDO
      enddo

        do k=2,nk-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem)
        do j=1,nj
        do i=1,ni
          r1 = dum7(i,j,k)
          r2 = dum8(i,j,k)
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = s(i,j,k) + temb*(r1*s(i,j,k-1)-(r1+r2)*s(i,j,k)+r2*s(i,j,k+1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
        enddo
        enddo

        enddo

        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem)
        do j=1,nj
        do i=1,ni
          r1 = dum7(i,j,k)
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = s(i,j,k) + temb*( r1*s(i,j,k-1)-r1*s(i,j,k) )  &
                         - temc*dumz(i,j,nk+1)*mh(i,j,nk)*rr(i,j,nk)
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dumz(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = rho(i,j,k)*(dumz(i,j,k)-s(i,j,k))*rdt
        enddo
        enddo

        do k=nk-1,1,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        DO i=1,ni
          dumz(i,j,k)=dum1(i,j,k)*dumz(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = rho(i,j,k)*(dumz(i,j,k)-s(i,j,k))*rdt
        ENDDO
        enddo

        enddo

  ENDIF  ifimpls

      end subroutine turbsz


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,          &
                       rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                       turbx,turby,turbz,dumx,dumy,dumz,t,tten,kmh,kmv)
      use input
      use constants
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dumx,dumy,dumz
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: t
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tten
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv

      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrf
      real :: cfa,cfb,cfc,cfd

!---------------------------------------------------------------

    IF(.not.terrain_flag)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj+1
        do i=1,ni+1
          !  x-direction
          ! note:  K is multiplied by 2:
          dumx(i,j,k)= -0.25*( rf(i,j,k)+rf(i-1,j,k) )   &
                       *2.0*( kmh(i,j,k)+kmh(i-1,j,k) )   &
                           *(   t(i,j,k)-  t(i-1,j,k) )*rdx*uf(i)
          !  y-direction
          ! note:  K is multiplied by 2:
          dumy(i,j,k)= -0.25*( rf(i,j,k)+rf(i,j-1,k) )   &
                       *2.0*( kmh(i,j,k)+kmh(i,j-1,k) )   &
                           *(   t(i,j,k)-  t(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            dumx(1,j,k) = dumx(2,j,k)
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            dumx(ni+1,j,k) = dumx(ni,j,k)
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            dumy(i,1,k) = dumy(i,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            dumy(i,nj+1,k) = dumy(i,nj,k)
          enddo
        ENDIF

        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
          turby(i,j,k) = -(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

!---------------------------------------------------------------
!  Cartesian with terrain:

    ELSE

      ! turbz stores t at s-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          turbz(i,j,k) = 0.5*(t(i,j,k)+t(i,j,k+1))
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        ! x-flux:
        do j=1,nj
        do i=1,ni+1
          ! note:  K is multiplied by 2:
          dumx(i,j,k)= -0.25*( gz(i,j)*rf(i,j,k)+gz(i-1,j)*rf(i-1,j,k) )                 &
                       *2.0*( kmh(i,j,k)+kmh(i-1,j,k) )*(                                &
                            (t(i,j,k)*rgz(i,j)-t(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)          &
                     +0.5*( (zt-sigma(k  ))*(turbz(i-1,j,k  )+turbz(i,j,k  ))            &
                           -(zt-sigma(k-1))*(turbz(i-1,j,k-1)+turbz(i,j,k-1))            &
                          )*rds(k)*(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)                       &
                                                        )
        enddo
        enddo

        ! y-flux:
        do j=1,nj+1
        do i=1,ni
          ! note:  K is multiplied by 2:
          dumy(i,j,k)= -0.25*( gz(i,j)*rf(i,j,k)+gz(i,j-1)*rf(i,j-1,k) )                 &
                       *2.0*( kmh(i,j,k)+kmh(i,j-1,k) )*(                                &
                            (t(i,j,k)*rgz(i,j)-t(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)          &
                     +0.5*( (zt-sigma(k  ))*(turbz(i,j-1,k  )+turbz(i,j,k  ))            &
                           -(zt-sigma(k-1))*(turbz(i,j-1,k-1)+turbz(i,j,k-1))            &
                          )*rds(k)*(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)                       &
                                                        )
        enddo
        enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          dumx(i,j,   1)=0.0
          dumx(i,j,nk+1)=0.0
          dumy(i,j,   1)=0.0
          dumy(i,j,nk+1)=0.0
        enddo
        enddo

      ! turbz stores dumx at s-pts:
      !  dumz stores dumy at s-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          turbz(i,j,k)=0.5*(dumx(i,j,k)+dumx(i,j,k+1))
           dumz(i,j,k)=0.5*(dumy(i,j,k)+dumy(i,j,k+1))
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        ! x-tendency:
        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -(dumx(i+1,j,k)*rgzu(i+1,j)-dumx(i,j,k)*rgzu(i,j))*gz(i,j)*rdx*uh(i) &
                         -0.5*( (zt-sigma(k  ))*(turbz(i,j,k  )+turbz(i+1,j,k  ))             &
                               -(zt-sigma(k-1))*(turbz(i,j,k-1)+turbz(i+1,j,k-1))             &
                              )*rds(k)*(rgzu(i+1,j)-rgzu(i,j))*gz(i,j)*rdx*uh(i)
        enddo
        enddo

        ! y-tendency:
        do j=1,nj
        do i=1,ni
          turby(i,j,k) = -(dumy(i,j+1,k)*rgzv(i,j+1)-dumy(i,j,k)*rgzv(i,j))*gz(i,j)*rdy*vh(j) &
                         -0.5*( (zt-sigma(k  ))*( dumz(i,j,k  )+ dumz(i,j+1,k  ))             &
                               -(zt-sigma(k-1))*( dumz(i,j,k-1)+ dumz(i,j+1,k-1))             &
                              )*rds(k)*(rgzv(i,j+1)-rgzv(i,j))*gz(i,j)*rdy*vh(j)
        enddo
        enddo

      enddo

    ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=2,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

!---------------------------------------------------------------------
!  z-direction

  ifimplt:  &
  IF( doimpl.eq.0 )THEN
      ! explicit vertical turbulence:

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        ! note:  K is multiplied by 2:
        dumz(i,j,k) = -(kmv(i,j,k)+kmv(i,j,k+1))*(t(i,j,k+1)-t(i,j,k))*rdz*mh(i,j,k)*rho(i,j,k)
      enddo
      enddo
      enddo

      if( bbc.eq.3 )then
      ! 210504: tke is approximately constant in the surface layer (thus dt/dz=0)
      do j=1,nj
      do i=1,ni
        dumz(i,j,1) = 0.0
      enddo
      enddo
      endif

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k) = -(dumz(i,j,k)-dumz(i,j,k-1))*rdz*mf(i,j,k)
        tten(i,j,k)=tten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
      enddo
      enddo
      enddo

  ELSE

      ! implicit vertical turbulence:

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb =      dt*vibeta*rdz*rdz
      temc = dt*rdz

        k=2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r2 = (kmv(i,j,k  )+kmv(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = t(i,j,k) + temb*( r2*t(i,j,k+1)-r2*t(i,j,k) )
!!!          tem = -(kmv(i,j,k-1)+kmv(i,j,k))*(t(i,j,k)-t(i,j,k-1))*rdz*mh(i,j,k-1)*rho(i,j,k-1)
!!!          cfd = cfd + temc*tem*rrf
          tem = 1.0/cfb
          dumx(i,j,k) = -cfc*tem
          dumy(i,j,k) =  cfd*tem
        enddo
        enddo

        do k=3,(nk-1)

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmv(i,j,k-1)+kmv(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          r2 = (kmv(i,j,k  )+kmv(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = t(i,j,k) + temb*(r2*t(i,j,k+1)-(r1+r2)*t(i,j,k)+r1*t(i,j,k-1))
          tem = 1.0/(cfa*dumx(i,j,k-1)+cfb)
          dumx(i,j,k) = -cfc*tem
          dumy(i,j,k) = (cfd-cfa*dumy(i,j,k-1))*tem
        enddo
        enddo

        enddo

        k = nk

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmv(i,j,k-1)+kmv(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = t(i,j,k) + temb*( -r1*t(i,j,k)+r1*t(i,j,k-1) )
!!!          tem = -(kmv(i,j,k)+kmv(i,j,k+1))*(t(i,j,k+1)-t(i,j,k))*rdz*mh(i,j,k)*rho(i,j,k)
!!!          cfd = cfd - temc*tem*rrf
          tem = 1.0/(cfa*dumx(i,j,k-1)+cfb)
          dumz(i,j,k) = (cfd-cfa*dumy(i,j,k-1))*tem
          turbz(i,j,k) = rf(i,j,k)*(dumz(i,j,k)-t(i,j,k))*rdt
          tten(i,j,k)=tten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
        enddo
        enddo

        do k=(nk-1),2,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dumz(i,j,k) = dumx(i,j,k)*dumz(i,j,k+1)+dumy(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dumz(i,j,k)-t(i,j,k))*rdt
          tten(i,j,k)=tten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
        enddo
        enddo

        enddo

  ENDIF  ifimplt

!---------------------------------------------------------------------

      if(timestats.ge.1) time_ttend=time_ttend+mytime()

      return
      end subroutine turbt


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                       turbx,turby,turbz,dum1,dum2,dum3,dum7,dum8,dum9,u,uten,w,t11,t12,t13,t22,kmv,cm0, &
                       kmw,ufw,u1b,u2pt,ufwk,doubud,udiag)
      use input
      use constants
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2,dum3,dum7,dum8,dum9
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(kb:ke) :: kmw,ufw,u1b
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt
      real, intent(in),    dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: ufwk
      logical, intent(in) :: doubud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag

      integer :: i,j,k,i1,i2
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rru0
      real :: cfa,cfb,cfc,cfd

!---------------------------------------------------------------

  dohoriz:  &
  IF( dohturb )THEN

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        !  y-direction
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=(t11(i,j,k)-t11(i-1,j,k))*rdx*uf(i)
          turby(i,j,k)=(t12(i,j+1,k)-t12(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(j,k)
      do k=1,nk

        do j=1,nj
        turbx(1,j,k)=0.0
        do i=2,ni+1
          turbx(i,j,k) = ( arf2(i)*arf2(i)*t11(i,j,k) - arf1(i)*arf1(i)*t11(i-1,j,k) )*rdx*uf(i)
        enddo
        IF(ebc.eq.3.or.ebc.eq.4)THEN
          turbx(ni+1,j,k)=0.0
        ENDIF
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------
!  Terrain:

  ELSE

      ! dum1 stores t11 at w-pts:
      ! dum2 stores t12 at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=1,nj+1

          ! lowest model level:
          do i=0,ni+1
            dum1(i,j,1) = cgs1*t11(i,j,1)+cgs2*t11(i,j,2)+cgs3*t11(i,j,3)
            dum2(i,j,1) = cgs1*t12(i,j,1)+cgs2*t12(i,j,2)+cgs3*t12(i,j,3)
          enddo

          ! upper-most model level:
          do i=0,ni+1
            dum1(i,j,nk+1) = cgt1*t11(i,j,nk)+cgt2*t11(i,j,nk-1)+cgt3*t11(i,j,nk-2)
            dum2(i,j,nk+1) = cgt1*t12(i,j,nk)+cgt2*t12(i,j,nk-1)+cgt3*t12(i,j,nk-2)
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+1
            dum1(i,j,k) = r1*t11(i,j,k-1)+r2*t11(i,j,k)
            dum2(i,j,k) = r1*t12(i,j,k-1)+r2*t12(i,j,k)
          enddo
          enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=gzu(i,j)*(t11(i,j,k)*rgz(i,j)-t11(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)  &
                      +0.5*( gxu(i,j,k+1)*(dum1(i-1,j,k+1)+dum1(i,j,k+1))                &
                            -gxu(i,j,k  )*(dum1(i-1,j,k  )+dum1(i,j,k  )) )*rdsf(k)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni+1
          r1 = 0.25*((rgz(i-1,j-1)+rgz(i,j))+(rgz(i-1,j)+rgz(i,j-1)))
          r2 = 0.25*((rgz(i-1,j+1)+rgz(i,j))+(rgz(i-1,j)+rgz(i,j+1)))
          turby(i,j,k)=gzu(i,j)*(t12(i,j+1,k)*r2-t12(i,j,k)*r1)*rdy*vh(j)      &
                      +0.5*( (zt-sigmaf(k+1))*(dum2(i,j,k+1)+dum2(i,j+1,k+1))  &
                            -(zt-sigmaf(k  ))*(dum2(i,j,k  )+dum2(i,j+1,k  ))  &
                           )*gzu(i,j)*(r2-r1)*rdy*vh(j)*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni+1,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni+1
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni+1
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

  ELSE  dohoriz

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbx(i,j,k)=0.0
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dohoriz

!-----------------------------------------------------------------
!  z-direction

  dovert:  &
  IF( dovturb )THEN

      call       turbuz(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                       turbz,dum1,dum2,dum3,dum7,dum9,u ,w ,t11,t12,t13,t22,kmv, &
                       kmw,ufw,u1b,u2pt,ufwk)

  ELSE  dovert

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dovert


!-----------------------------------------------------------------
!  Tendencies:

    IF(axisymm.eq.0)THEN

      IF( cm1setup.eq.4 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,rru0)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          if( 0.5*(cm0(i-1,j)+cm0(i,j)) .le. cmemin )then
            ! zero-out turb tendencies outside of LES domain:
            turbx(i,j,k) = 0.0
            turby(i,j,k) = 0.0
            turbz(i,j,k) = 0.0
          endif
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          uten(i,j,k)=uten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rru0
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,rru0)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          uten(i,j,k)=uten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rru0
        enddo
        enddo
        enddo
      ENDIF

    ELSE

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,rru0)
      do k=1,nk
      do j=1,nj
      do i=2,ni+1
        rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
        uten(i,j,k)=uten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))*rru0
      enddo
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------
!  Diagnostics:

      IF( doubud )THEN
      if( ud_hturb.ge.1 .and. ud_vturb.ge.1 )then
        if( axisymm.eq.0 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k,rru0)
          do k=1,nk
          do j=1,nj
          do i=1,ni+1
            rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
            udiag(i,j,k,ud_hturb) = (turbx(i,j,k)+turby(i,j,k))*rru0
            udiag(i,j,k,ud_vturb) = turbz(i,j,k)*rru0
          enddo
          enddo
          enddo
        else
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k,rru0)
          do k=1,nk
          do j=1,nj
          do i=2,ni+1
            rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
            udiag(i,j,k,ud_hturb) = turbx(i,j,k)*rru0
            udiag(i,j,k,ud_vturb) = turbz(i,j,k)*rru0
          enddo
          enddo
          enddo
        endif
      endif
      ENDIF

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_ttend=time_ttend+mytime()

      return
      end subroutine turbu


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbuz(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                       turbz,dum1,dum2,dum3,dum7,dum9,u ,w ,t11,t12,t13,t22,kmv, &
                       kmw,ufw,u1b,u2pt,ufwk)
      use input
      use constants
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbz,dum1,dum2,dum3,dum7,dum9
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, intent(in), dimension(kb:ke) :: kmw,ufw,u1b
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt
      real, intent(in),    dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: ufwk

      integer :: i,j,k,i1,i2
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rru0
      real :: cfa,cfb,cfc,cfd

  ifimplu:  &
  IF( doimpl.eq.0 )THEN
      ! explicit vertical turbulence:

      tem = rdz*0.5

      i1 = 1
      i2 = ni+1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbz(i,j,k)=(t13(i,j,k+1)-t13(i,j,k))*tem*(mh(i-1,j,k)+mh(i,j,k))
      enddo
      enddo
      enddo

  ELSE

      ! implicit vertical turbulence:

    check_grid:  &
    IF(axisymm.eq.0)THEN
      ! Cartesian grid:

      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temc =  0.5*dt*rdz

      i1 = 1
      i2 = ni+1

      IF( .not. terrain_flag )THEN
        ! without terrain:

        !--------
        k = 1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-r2*u(i,j,k) )  &
                         - temc*t13(i,j,1)*(mh(i-1,j,1)+mh(i,j,1))*rru0
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
          dum7(i,j,1) = 0.0
        enddo
        enddo
        !--------
        do k=2,nk-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          ! explicit piece ... dwdx term
          dum7(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)   &
                     *0.25*( kmv(i-1,j,k)+kmv(i,j,k) )  &
                          *( rf(i-1,j,k)+rf(i,j,k) )
        enddo
        enddo

        enddo
        !--------
        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = u(i,j,k) + temb*( -r1*u(i,j,k)+r1*u(i,j,k-1) )  &
                         + temc*t13(i,j,nk+1)*(mh(i-1,j,nk)+mh(i,j,nk))*rru0
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
          ! explicit piece ... dwdx term
          dum7(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)   &
                     *0.25*( kmv(i-1,j,k)+kmv(i,j,k) )  &
                          *( rf(i-1,j,k)+rf(i,j,k) )
          turbz(i,j,k)=turbz(i,j,k)+(0.0-dum7(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo
        !--------

        do k=nk-1,1,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        do i=i1,i2
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
          ! explicit piece ... dwdx term
          turbz(i,j,k)=turbz(i,j,k)+(dum7(i,j,k+1)-dum7(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo

        enddo

      ELSE
        ! with terrain:

        ! dum9 stores w at scalar-pts:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          dum9(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
        enddo
        enddo
        ENDDO

        !--------
        k = 1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-r2*u(i,j,k) )  &
                         - temc*t13(i,j,1)*(mh(i-1,j,1)+mh(i,j,1))*rru0
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
          dum7(i,j,1) = 0.0
        enddo
        enddo
        !--------
        do k=2,nk-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          ! explicit piece ... dwdx term
          dum7(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)          &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum9(i,j,k  )+dum9(i-1,j,k  ))    &
                               -(zt-sigma(k-1))*(dum9(i,j,k-1)+dum9(i-1,j,k-1)) )  &
                             *(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)
          dum7(i,j,k)=dum7(i,j,k)*0.25*( kmv(i-1,j,k)+kmv(i,j,k) )                 &
                                      *( gz(i-1,j)*rf(i-1,j,k)+gz(i,j)*rf(i,j,k) )
        enddo
        enddo

        enddo
        !--------
        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = u(i,j,k) + temb*( -r1*u(i,j,k)+r1*u(i,j,k-1) )  &
                         + temc*t13(i,j,nk+1)*(mh(i-1,j,nk)+mh(i,j,nk))*rru0
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
          ! explicit piece ... dwdx term
          dum7(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)          &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum9(i,j,k  )+dum9(i-1,j,k  ))    &
                               -(zt-sigma(k-1))*(dum9(i,j,k-1)+dum9(i-1,j,k-1)) )  &
                             *(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)
          dum7(i,j,k)=dum7(i,j,k)*0.25*( kmv(i-1,j,k)+kmv(i,j,k) )                 &
                                      *( gz(i-1,j)*rf(i-1,j,k)+gz(i,j)*rf(i,j,k) )
          turbz(i,j,k)=turbz(i,j,k)+(0.0-dum7(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo
        !--------

        do k=nk-1,1,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        do i=i1,i2
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
          ! explicit piece ... dwdx term
          turbz(i,j,k)=turbz(i,j,k)+(dum7(i,j,k+1)-dum7(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo

        enddo

      ENDIF

      !------------------------------------------------------------
      !------------------------------------------------------------
      !------------------------------------------------------------

    ELSEIF(axisymm.eq.1)THEN
      ! axisymmetric grid:

      rdt = 0.5/dt
      tema = -0.25*dt*vialpha*rdz*rdz
      temb =  0.25*dt*vibeta*rdz*rdz
      temc =  dt*rdz

      i1 = 2
      i2 = ni+1

        !--------
        k = 1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          tem = mh(1,1,k)*rru0
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*mf(1,1,k+1)   &
              *(arf1(i)*rf(i-1,j,k+1)+arf2(i)*rf(i,j,k+1))*tem
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-r2*u(i,j,k) )  &
                         - temc*t13(i,j,1)*mh(1,1,1)*rru0
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
          dum7(i,j,1) = 0.0
        enddo
        enddo
        !--------
        do k=2,nk-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          tem = mh(1,1,k)*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*mf(1,1,k  )   &
              *(arf1(i)*rf(i-1,j,k  )+arf2(i)*rf(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*mf(1,1,k+1)   &
              *(arf1(i)*rf(i-1,j,k+1)+arf2(i)*rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = u(i,j,k) + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          ! explicit piece ... dwdx term
          dum7(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)
          dum7(i,j,k)=dum7(i,j,k)*0.25*( arf1(i)*kmv(i-1,j,k)+arf2(i)*kmv(i,j,k) )  &
                                      *(arf1(i)*rf(i-1,j,k)+arf2(i)*rf(i,j,k))
        enddo
        enddo

        enddo
        !--------
        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rru0)
        do j=1,nj
        do i=i1,i2
          rru0 = 1.0/(0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k)))
          tem = mh(1,1,k)*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*mf(1,1,k  )   &
              *(arf1(i)*rf(i-1,j,k  )+arf2(i)*rf(i,j,k  ))*tem
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = u(i,j,k) + temb*( -r1*u(i,j,k)+r1*u(i,j,k-1) )  &
                         + temc*t13(i,j,nk+1)*mh(1,1,nk)*rru0
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = (arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
          ! explicit piece ... dwdx term
          dum7(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)
          dum7(i,j,k)=dum7(i,j,k)*0.25*( arf1(i)*kmv(i-1,j,k)+arf2(i)*kmv(i,j,k) )  &
                                      *(arf1(i)*rf(i-1,j,k)+arf2(i)*rf(i,j,k))
          turbz(i,j,k)=turbz(i,j,k)+(0.0-dum7(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo
        !--------

        do k=nk-1,1,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        do i=i1,i2
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))*(dum3(i,j,k)-u(i,j,k))*rdt
          ! explicit piece ... dwdx term
          turbz(i,j,k)=turbz(i,j,k)+(dum7(i,j,k+1)-dum7(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
        enddo
        enddo

        enddo

    ENDIF  check_grid

    !------------------------------------------------------------


  ENDIF  ifimplu

      IF(axisymm.eq.1)THEN
        !$omp parallel do default(shared)   &
        !$omp private(k)
        DO k=1,nk
          turbz(1,1,k) = 0.0
        ENDDO
        IF( ebc.eq.3 .or. ebc.eq.4 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            turbz(ni+1,1,k)=0.0
          enddo
        ENDIF
      ENDIF

!-----------------------------------------------------------------
!  2nd part of 2-part model:

      IF( dot2p )THEN

        do k=2,ntwk
!!!          tem = 0.5*rdz*mf(1,1,k)
          do j=1,nj
          do i=i1,i2
!!!            dum2(i,j,k) = tem*kmw(k)*(u1b(k)-u1b(k-1))*(rf(i-1,j,k)+rf(i,j,k))
!!!            dum2(i,j,k) = ufw(k)*0.5*(rf(i-1,j,k)+rf(i,j,k))
            dum2(i,j,k) = 0.5*(rf(i-1,j,k)*ufwk(i-1,j,k)+rf(i,j,k)*ufwk(i,j,k))
          enddo
          enddo
        enddo

        do j=1,nj
        do i=i1,i2
          dum2(i,j,1) = 0.0
          dum2(i,j,ntwk+1) = 0.0
        enddo
        enddo

        do k=1,ntwk
          tem = rdz*mh(1,1,k)
          do j=1,nj
          do i=i1,i2
            turbz(i,j,k) = turbz(i,j,k)+tem*(dum2(i,j,k+1)-dum2(i,j,k))
            ! 211207: save 2-pt tendency
            rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
            u2pt(i,j,k)=tem*(dum2(i,j,k+1)-dum2(i,j,k))*rru0
          enddo
          enddo
        enddo

      ENDIF

!-----------------------------------------------------------------

      return
      end subroutine turbuz


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                       turbx,turby,turbz,dum1,dum2,dum3,dum7,dum8,dum9,v,vten,w,t12,t22,t23,kmv,cm0, &
                       kmw,vfw,v1b,v2pt,vfwk,dovbud,vdiag)
      use input
      use constants
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gyv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2,dum3,dum7,dum8,dum9
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t12,t22,t23
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(in), dimension(kb:ke) :: kmw,vfw,v1b
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: v2pt
      real, intent(in),    dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: vfwk
      logical, intent(in) :: dovbud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
 
      integer :: i,j,k,j1,j2
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrv0
      real :: cfa,cfb,cfc,cfd

!---------------------------------------------------------------

  dohoriz:  &
  IF( dohturb )THEN

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        !  y-direction
        do j=1,nj+1
        do i=1,ni
          turbx(i,j,k)=(t12(i+1,j,k)-t12(i,j,k))*rdx*uh(i)
          turby(i,j,k)=(t22(i,j,k)-t22(i,j-1,k))*rdy*vf(j)
        enddo
        enddo

      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(arh2(i)*arh2(i)*t12(i+1,j,k)-arh1(i)*arh1(i)*t12(i,j,k))*rdx*uh(i)
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------
!  Terrain:

  ELSE

      ! dum1 stores t12 at w-pts:
      ! dum2 stores t22 at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=0,nj+1

          ! lowest model level:
          do i=1,ni+1
            dum1(i,j,1) = cgs1*t12(i,j,1)+cgs2*t12(i,j,2)+cgs3*t12(i,j,3)
            dum2(i,j,1) = cgs1*t22(i,j,1)+cgs2*t22(i,j,2)+cgs3*t22(i,j,3)
          enddo

          ! upper-most model level:
          do i=1,ni+1
            dum1(i,j,nk+1) = cgt1*t12(i,j,nk)+cgt2*t12(i,j,nk-1)+cgt3*t12(i,j,nk-2)
            dum2(i,j,nk+1) = cgt1*t22(i,j,nk)+cgt2*t22(i,j,nk-1)+cgt3*t22(i,j,nk-2)
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=1,ni+1
            dum1(i,j,k) = r1*t12(i,j,k-1)+r2*t12(i,j,k)
            dum2(i,j,k) = r1*t22(i,j,k-1)+r2*t22(i,j,k)
          enddo
          enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do k=1,nk

        !  x-direction
        do j=1,nj+1
        do i=1,ni
          r1 = 0.25*((rgz(i-1,j-1)+rgz(i,j))+(rgz(i-1,j)+rgz(i,j-1)))
          r2 = 0.25*((rgz(i+1,j-1)+rgz(i,j))+(rgz(i+1,j)+rgz(i,j-1)))
          turbx(i,j,k)=gzv(i,j)*(t12(i+1,j,k)*r2-t12(i,j,k)*r1)*rdx*uh(i)      &
                      +0.5*( (zt-sigmaf(k+1))*(dum1(i,j,k+1)+dum1(i+1,j,k+1))  &
                            -(zt-sigmaf(k  ))*(dum1(i,j,k  )+dum1(i+1,j,k  ))  &
                           )*gzv(i,j)*(r2-r1)*rdx*uh(i)*rdsf(k)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          turby(i,j,k)=gzv(i,j)*(t22(i,j,k)*rgz(i,j)-t22(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)  &
                      +0.5*( gyv(i,j,k+1)*(dum2(i,j-1,k+1)+dum2(i,j,k+1))                &
                            -gyv(i,j,k  )*(dum2(i,j-1,k  )+dum2(i,j,k  )) )*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj+1
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj+1
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj+1,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

  ELSE  dohoriz

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbx(i,j,k)=0.0
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dohoriz

!-----------------------------------------------------------------
!  z-direction

  dovert:  &
  IF( dovturb )THEN

      call       turbvz(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                       turbz,dum1,dum2,dum3,dum7,dum9,v ,w ,t12,t22,t23,kmv, &
                       kmw,vfw,v1b,v2pt,vfwk)

  ELSE  dovert

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dovert


!-----------------------------------------------------------------
!  Tendencies:

    IF(axisymm.eq.0)THEN

      IF( cm1setup.eq.4 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,rrv0)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          if( 0.5*(cm0(i,j-1)+cm0(i,j)) .le. cmemin )then
            ! zero-out turb tendencies outside of LES domain:
            turbx(i,j,k) = 0.0
            turby(i,j,k) = 0.0
            turbz(i,j,k) = 0.0
          endif
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          vten(i,j,k)=vten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rrv0
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,rrv0)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          vten(i,j,k)=vten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rrv0
        enddo
        enddo
        enddo
      ENDIF

    ELSE

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        vten(i,j,k)=vten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))*rr(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------
!  Diagnostics:

      IF( dovbud )THEN
      if( vd_hturb.ge.1 .and. vd_vturb.ge.1 )then
        if( axisymm.eq.0 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k,rrv0)
          do k=1,nk
          do j=1,nj+1
          do i=1,ni
            rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
            vdiag(i,j,k,vd_hturb) = (turbx(i,j,k)+turby(i,j,k))*rrv0
            vdiag(i,j,k,vd_vturb) = turbz(i,j,k)*rrv0
          enddo
          enddo
          enddo
        else
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            vdiag(i,j,k,vd_hturb) = turbx(i,j,k)*rr(i,j,k)
            vdiag(i,j,k,vd_vturb) = turbz(i,j,k)*rr(i,j,k)
            vdiag(i,2,k,vd_hturb) = vdiag(i,1,k,vd_hturb)
            vdiag(i,2,k,vd_vturb) = vdiag(i,1,k,vd_vturb)
          enddo
          enddo
          enddo
        endif
      endif
      ENDIF

!-------------------------------------------------------------------
!  All done
 
      if(timestats.ge.1) time_ttend=time_ttend+mytime()
 
      return
      end subroutine turbv


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbvz(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                       turbz,dum1,dum2,dum3,dum7,dum9,v ,w ,t12,t22,t23,kmv, &
                       kmw,vfw,v1b,v2pt,vfwk)
      use input
      use constants
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gyv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbz,dum1,dum2,dum3,dum7,dum9
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t12,t22,t23
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      real, intent(in), dimension(kb:ke) :: kmw,vfw,v1b
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: v2pt
      real, intent(in),    dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: vfwk
 
      integer :: i,j,k,j1,j2
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrv0
      real :: cfa,cfb,cfc,cfd

  ifimplv:  &
  IF( doimpl.eq.0 )THEN
      ! explicit vertical turbulence:

      tem = rdz*0.5

      j1 = 1
      j2 = nj+1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbz(i,j,k)=(t23(i,j,k+1)-t23(i,j,k))*tem*(mh(i,j-1,k)+mh(i,j,k))
      enddo
      enddo
      enddo

  ELSE

      ! implicit vertical turbulence:

      IF( .not. terrain_flag )THEN

        ! without terrain:

        rdt = 0.5/dt
        tema = -0.0625*dt*vialpha*rdz*rdz
        temb =  0.0625*dt*vibeta*rdz*rdz
        temc =  0.5*dt*rdz

        if( axisymm.eq.1 )then
          j1 = 1
          j2 = 1
        else
          j1 = 1
          j2 = nj+1
        endif

        !--------
        k = 1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
        do j=j1,j2
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = v(i,j,k) + temb*( r2*v(i,j,k+1)-r2*v(i,j,k) )  &
                         - temc*t23(i,j,1)*(mh(i,j-1,1)+mh(i,j,1))*rrv0
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
          dum7(i,j,1) = 0.0
        enddo
        enddo
        !--------
        do k=2,nk-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
        do j=j1,j2
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = v(i,j,k) + temb*( r2*v(i,j,k+1)-(r1+r2)*v(i,j,k)+r1*v(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          ! explicit piece ... dwdy term
          dum7(i,j,k)=(w(i,j,k)-w(i,j-1,k))*rdy*vf(j)   &
                     *0.25*( kmv(i,j-1,k)+kmv(i,j,k) )  &
                          *( rf(i,j-1,k)+rf(i,j,k) )
        enddo
        enddo

        enddo
        !--------
        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
        do j=j1,j2
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = v(i,j,k) + temb*( -r1*v(i,j,k)+r1*v(i,j,k-1) )  &
                         + temc*t23(i,j,nk+1)*(mh(i,j-1,nk)+mh(i,j,nk))*rrv0
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(dum3(i,j,k)-v(i,j,k))*rdt
          ! explicit piece ... dwdy term
          dum7(i,j,k)=(w(i,j,k)-w(i,j-1,k))*rdy*vf(j)   &
                     *0.25*( kmv(i,j-1,k)+kmv(i,j,k) )  &
                          *( rf(i,j-1,k)+rf(i,j,k) )
          turbz(i,j,k)=turbz(i,j,k)+(0.0-dum7(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
        enddo
        enddo
        !--------

        do k=nk-1,1,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=j1,j2
        do i=1,ni
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(dum3(i,j,k)-v(i,j,k))*rdt
          ! explicit piece ... dwdy term
          turbz(i,j,k)=turbz(i,j,k)+(dum7(i,j,k+1)-dum7(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
        enddo
        enddo

        enddo

    !------------------------------------------------------------

      ELSE

        ! with terrain:

        ! dum9 stores w at scalar-pts:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          dum9(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
        enddo
        enddo
        ENDDO

        rdt = 0.5/dt
        tema = -0.0625*dt*vialpha*rdz*rdz
        temb =  0.0625*dt*vibeta*rdz*rdz
        temc =  0.5*dt*rdz

        j1 = 1
        j2 = nj+1

        !--------
        k = 1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
        do j=j1,j2
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = v(i,j,k) + temb*( r2*v(i,j,k+1)-r2*v(i,j,k) )  &
                         - temc*t23(i,j,1)*(mh(i,j-1,1)+mh(i,j,1))*rrv0
          tem = 1.0/cfb
          dum1(i,j,1)=-cfc*tem
          dum2(i,j,1)= cfd*tem
          dum7(i,j,1) = 0.0
        enddo
        enddo
        !--------
        do k=2,nk-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
        do j=j1,j2
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = v(i,j,k) + temb*( r2*v(i,j,k+1)-(r1+r2)*v(i,j,k)+r1*v(i,j,k-1) )
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k)=-cfc*tem
          dum2(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          dum7(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)          &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum9(i,j,k  )+dum9(i,j-1,k  ))    &
                               -(zt-sigma(k-1))*(dum9(i,j,k-1)+dum9(i,j-1,k-1)) )  &
                             *(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)
          ! explicit piece ... dwdy term
          dum7(i,j,k)=dum7(i,j,k)*0.25*( kmv(i,j-1,k)+kmv(i,j,k) )                 &
                                      *( gz(i,j-1)*rf(i,j-1,k)+gz(i,j)*rf(i,j,k) )
        enddo
        enddo

        enddo
        !--------
        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrv0)
        do j=j1,j2
        do i=1,ni
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = v(i,j,k) + temb*( -r1*v(i,j,k)+r1*v(i,j,k-1) )  &
                         + temc*t23(i,j,nk+1)*(mh(i,j-1,nk)+mh(i,j,nk))*rrv0
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k)=(cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(dum3(i,j,k)-v(i,j,k))*rdt
          ! explicit piece ... dwdy term
          dum7(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)          &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum9(i,j,k  )+dum9(i,j-1,k  ))    &
                               -(zt-sigma(k-1))*(dum9(i,j,k-1)+dum9(i,j-1,k-1)) )  &
                             *(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)
          dum7(i,j,k)=dum7(i,j,k)*0.25*( kmv(i,j-1,k)+kmv(i,j,k) )                 &
                                      *( gz(i,j-1)*rf(i,j-1,k)+gz(i,j)*rf(i,j,k) )
          turbz(i,j,k)=turbz(i,j,k)+(0.0-dum7(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
        enddo
        enddo
        !--------

        do k=nk-1,1,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=j1,j2
        do i=1,ni
          dum3(i,j,k)=dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(dum3(i,j,k)-v(i,j,k))*rdt
          ! explicit piece ... dwdy term
          turbz(i,j,k)=turbz(i,j,k)+(dum7(i,j,k+1)-dum7(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
        enddo
        enddo

        enddo

      ENDIF

  ENDIF  ifimplv

!-----------------------------------------------------------------
!  2nd part of 2-part model:

      IF( dot2p )THEN

        do k=2,ntwk
!!!          tem = 0.5*rdz*mf(1,1,k)
          do j=j1,j2
          do i=1,ni
!!!            dum2(i,j,k) = tem*kmw(k)*(v1b(k)-v1b(k-1))*(rf(i,j-1,k)+rf(i,j,k))
!!!            dum2(i,j,k) = vfw(k)*0.5*(rf(i,j-1,k)+rf(i,j,k))
            dum2(i,j,k) = 0.5*(rf(i,j-1,k)*vfwk(i,j-1,k)+rf(i,j,k)*vfwk(i,j,k))
          enddo
          enddo
        enddo

        do j=j1,j2
        do i=1,ni
          dum2(i,j,1) = 0.0
          dum2(i,j,ntwk+1) = 0.0
        enddo
        enddo

        do k=1,ntwk
          tem = rdz*mh(1,1,k)
          do j=j1,j2
          do i=1,ni
            turbz(i,j,k) = turbz(i,j,k)+tem*(dum2(i,j,k+1)-dum2(i,j,k))
            ! 211207: save 2-pt tendency
            rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
            v2pt(i,j,k)=tem*(dum2(i,j,k+1)-dum2(i,j,k))*rrv0
          enddo
          enddo
        enddo

      ENDIF

!-----------------------------------------------------------------

      return
      end subroutine turbvz


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
      subroutine turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                       turbx,turby,turbz,dum1,dum2,dum3,w,wten,t13,t23,t33,t22,kmh,cm0,  &
                       dowbud,wdiag)
      use input
      use constants
      implicit none

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgzu,rgzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2,dum3
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t13,t23,t33,t22
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      logical, intent(in) :: dowbud
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
 
      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrf
      real :: cfa,cfb,cfc,cfd

!----------------------------------------------------------------

  dohoriz:  &
  IF( dohturb )THEN

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        !  x-direction
        !  y-direction
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(t13(i+1,j,k)-t13(i,j,k))*rdx*uh(i)
          turby(i,j,k)=(t23(i,j+1,k)-t23(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE
      ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(arh2(i)*t13(i+1,j,k)-arh1(i)*t13(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!----------------------------------------------------------------

  ELSE
      ! Cartesian with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.25*( (t13(i,j,k+1)+t13(i+1,j,k+1)) &
                              +(t13(i,j,k  )+t13(i+1,j,k  )) )
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = 0.25*( (t23(i,j,k+1)+t23(i,j+1,k+1)) &
                              +(t23(i,j,k  )+t23(i,j+1,k  )) )
        enddo
        enddo
      enddo


!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=gz(i,j)*( t13(i+1,j,k)*rgzu(i+1,j)             &
                                -t13(i  ,j,k)*rgzu(i  ,j) )*rdx*uh(i) &
              +( (zt-sigma(k  ))*dum1(i,j,k  )                        &
                -(zt-sigma(k-1))*dum1(i,j,k-1)                        &
               )*gz(i,j)*(rgzu(i+1,j)-rgzu(i,j))*rdx*uh(i)*rds(k)
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          turby(i,j,k)=gz(i,j)*( t23(i,j+1,k)*rgzv(i,j+1)             &
                                -t23(i,j  ,k)*rgzv(i,j  ) )*rdy*vh(j) &
              +( (zt-sigma(k  ))*dum2(i,j,k  )                        &
                -(zt-sigma(k-1))*dum2(i,j,k-1)                        &
               )*gz(i,j)*(rgzv(i,j+1)-rgzv(i,j))*rdy*vh(j)*rds(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!-----------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=2,nk

        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          do j=1,nj
            turbx(1,j,k) = 0.0
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          do j=1,nj
            turbx(ni,j,k) = 0.0
          enddo
        ENDIF

        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          do i=1,ni
            turby(i,1,k) = 0.0
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          do i=1,ni
            turby(i,nj,k) = 0.0
          enddo
        ENDIF

      ENDDO
    ENDIF

  ELSE  dohoriz

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbx(i,j,k)=0.0
        turby(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dohoriz

!-----------------------------------------------------------------
!  z-direction

  dovert:  &
  IF( dovturb .or. ( cm1setup.eq.2 .and. ipbl.eq.2 ) )THEN

  ifimplw:  &
  IF( doimpl.eq.0 )THEN
      ! explicit vertical turbulence:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*rdz*mf(i,j,k)
      enddo
      enddo
      enddo

  ELSE

      ! implicit vertical turbulence:

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb =      dt*vibeta*rdz*rdz
      temc = dt*rdz

        k=2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r2 = (kmh(i,j,k  )+kmh(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfc = tema*r2
          cfb = 1.0 - cfc
          cfd = w(i,j,k) + temb*( r2*w(i,j,k+1)-r2*w(i,j,k) )  &
                         - temc*t33(i,j,k-1)*rrf
          tem = 1.0/cfb
          dum1(i,j,k) = -cfc*tem
          dum2(i,j,k) =  cfd*tem
        enddo
        enddo

        do k=3,(nk-1)

      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          r2 = (kmh(i,j,k  )+kmh(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*rrf
          cfa = tema*r1
          cfc = tema*r2
          cfb = 1.0 - cfa - cfc
          cfd = w(i,j,k) + temb*(r2*w(i,j,k+1)-(r1+r2)*w(i,j,k)+r1*w(i,j,k-1))
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum1(i,j,k) = -cfc*tem
          dum2(i,j,k) = (cfd-cfa*dum2(i,j,k-1))*tem
        enddo
        enddo

        enddo

      IF( axisymm.eq.0 )THEN
        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = w(i,j,k) + temb*( -r1*w(i,j,k)+r1*w(i,j,k-1) )  &
                         + temc*t33(i,j,k)*rrf
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k) = (cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = rf(i,j,k)*(dum3(i,j,k)-w(i,j,k))*rdt
        enddo
        enddo

        do k=(nk-1),2,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dum3(i,j,k)-w(i,j,k))*rdt
        enddo
        enddo

        enddo

      ELSE

        k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,r1,r2,cfa,cfb,cfc,cfd,tem,rrf)
        do j=1,nj
        do i=1,ni
          rrf = mf(i,j,k)/rf(i,j,k)
          r1 = (kmh(i,j,k-1)+kmh(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*rrf
          cfa = tema*r1
          cfb = 1.0 - cfa
          cfd = w(i,j,k) + temb*( -r1*w(i,j,k)+r1*w(i,j,k-1) )  &
                         + temc*t33(i,j,k)*rrf
          tem = 1.0/(cfa*dum1(i,j,k-1)+cfb)
          dum3(i,j,k) = (cfd-cfa*dum2(i,j,k-1))*tem
          turbz(i,j,k) = rf(i,j,k)*(dum3(i,j,k)-w(i,j,k))*rdt
        enddo
        enddo

        do k=(nk-1),2,-1

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = dum1(i,j,k)*dum3(i,j,k+1)+dum2(i,j,k)
          turbz(i,j,k) = rf(i,j,k)*(dum3(i,j,k)-w(i,j,k))*rdt
        enddo
        enddo

        enddo

      ENDIF

  ENDIF  ifimplw

  ELSE  dovert

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

  ENDIF  dovert

!---------------------------------------------------------------------
!  Tendencies:

    IF(axisymm.eq.0)THEN

      IF( cm1setup.eq.4 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          if( cm0(i,j) .le. cmemin )then
            ! zero-out turb tendencies outside of LES domain:
            turbx(i,j,k) = 0.0
            turby(i,j,k) = 0.0
            turbz(i,j,k) = 0.0
          endif

          wten(i,j,k)=wten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni

          wten(i,j,k)=wten(i,j,k)+((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))/rf(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

    ELSE

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        wten(i,j,k)=wten(i,j,k)+(turbx(i,j,k)+turbz(i,j,k))/rf(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!---------------------------------------------------------------------
!  Diagnostics:

      IF( dowbud )THEN
      if( wd_hturb.ge.1 .and. wd_vturb.ge.1 )then
        if( axisymm.eq.0 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wdiag(i,j,k,wd_hturb) = (turbx(i,j,k)+turby(i,j,k))/rf(i,j,k)
            wdiag(i,j,k,wd_vturb) = turbz(i,j,k)/rf(i,j,k)
          enddo
          enddo
          enddo
        else
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wdiag(i,j,k,wd_hturb) = turbx(i,j,k)/rf(i,j,k)
            wdiag(i,j,k,wd_vturb) = turbz(i,j,k)/rf(i,j,k)
          enddo
          enddo
          enddo
        endif
      endif
      ENDIF

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_ttend=time_ttend+mytime()
 
      return
      end subroutine turbw


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  END MODULE turbtend_module
