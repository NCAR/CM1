  MODULE interp_routines

  implicit none

  private

  public :: vinterp_weno5
  public :: vinterp_flx5,vinterp_flx6


    !-----------------------------------------------------!
    !       WENO only:                                    !
    ! formulation for weno "smoothness indicators"        !
          !  1 = original (eg, Jiang and Shu, 1996, JCP)
          !  2 = Borges et al. (2008, JCP)
    integer, parameter :: siform = 2
    !-----------------------------------------------------!


  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vinterp_weno5(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,pdefweno,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef,pdefweno
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5
    real :: bmin,bmax
    logical :: doit

    bmin = weps

!ccccccccccccccccccccccccccccccccccc
!  s-staggered variable:

  IF( stag.eq.1 )THEN

  pdefs:  &
  if( pdefweno.eq.1 )then

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,s1,s2,s3,s4,s5,bmax,doit)
    DO j=1,nj
      do k=3,nk-1
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          s1=a(i,j,k-3)
          s2=a(i,j,k-2)
          s3=a(i,j,k-1)
          s4=a(i,j,k  )
          s5=a(i,j,k+1)
        else
          s1=a(i,j,k+2)
          s2=a(i,j,k+1)
          s3=a(i,j,k  )
          s4=a(i,j,k-1)
          s5=a(i,j,k-2)
        endif
        doit = .false.
        bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5))
        if( bmax.gt.bmin ) doit = .true.
        IF( k.eq.3 .and. rrw(i,j,k).gt.0.0 ) doit = .false.
        IF( k.eq.(nk-1) .and. rrw(i,j,k).lt.0.0 ) doit = .false.
        IF(doit)THEN
          dumz(i,j,k) = weno5(s1,s2,s3,s4,s5,weps)
        ELSE
          dumz(i,j,k) = 0.0
        ENDIF
      enddo
      enddo

      !----

      ! flux at k=3 if w > 0
      k = 3
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).lt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        endif
      enddo

    ENDDO

  !--------------------------------------!
  else  pdefs

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    DO j=1,nj
      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      !----

      ! flux at k=3 if w > 0
      k = 3
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dumz(i,j,k) = weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        endif
      enddo

    ENDDO

  endif  pdefs

!ccccccccccccccccccccccccccccccccccc
!  u-staggered variable:

  ELSEIF( stag.eq.2 )THEN
    ! u-staggered:

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

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=4,nk-2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      !----

      k = 3
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo

    ENDDO

!ccccccccccccccccccccccccccccccccccc
!  v-staggered variable:

  ELSEIF( stag.eq.3 )THEN
    ! v-staggered:

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

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      !----

      k = 3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo

    ENDDO

!ccccccccccccccccccccccccccccccccccc
!  w-staggered variable:

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=3,nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = weno5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = weno5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = weno5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),weps)
        endif
      enddo

      !----

      k = 1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = 0.5*(a(i,j,k)+a(i,j,k+1))
        else
          dumz(i,j,k) = upstrpd(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = upstrpd(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = 0.5*(a(i,j,k)+a(i,j,k+1))
        endif
      enddo

    ENDDO

!ccccccccccccccccccccccccccccccccccc

  ENDIF

    end subroutine vinterp_weno5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vinterp_flx5(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=1,nj

      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

      !----

      k = 3
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        endif
      enddo

    ENDDO

  ELSEIF( stag.eq.2 )THEN
    ! u-staggered:

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

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=4,nk-2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

      !----

      k = 3
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo

    ENDDO


  ELSEIF( stag.eq.3 )THEN
    ! v-staggered:

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

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

      !----

      k = 3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=3,nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = flx5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

      !----

      k = 2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = flx5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

      k = nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = flx3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ))
        endif
      enddo

      !----

      k = 1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = 0.5*(a(i,j,k)+a(i,j,k+1))
        else
          dumz(i,j,k) = flx3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ))
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = flx3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = 0.5*(a(i,j,k)+a(i,j,k+1))
        endif
      enddo

    ENDDO

  ENDIF

    end subroutine vinterp_flx5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vinterp_flx6(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=1,nj

      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo

      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = (c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
      enddo
      enddo

    ENDDO

  ELSEIF( stag.eq.2 )THEN
    ! u-staggered:

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

!$omp parallel do default(shared)   &
!$omp private(i,j,k,cc1,cc2)
    DO j=1,nj

      do k=4,nk-2
      !dir$ vector always
      do i=i1,i2
        dumz(i,j,k) = flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=i1,i2
        dumz(i,j,k) = flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo

      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=i1,i2
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
      enddo
      enddo

    ENDDO

  ELSEIF( stag.eq.3 )THEN
    ! v-staggered:

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

!$omp parallel do default(shared)   &
!$omp private(i,j,k,cc1,cc2)
    DO j=j1,j2

      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo

      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumz(i,j,k) = (cc1*a(i,j,k-1)+cc2*a(i,j,k))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=1,nj

      do k=3,nk-2
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = flx6(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo

      do k=2,(nk-1),(nk-3)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = flx4(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

      do k=1,nk,(nk-1)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = 0.5*(a(i,j,k)+a(i,j,k+1))
      enddo
      enddo

    ENDDO

  ENDIF

    end subroutine vinterp_flx6

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx3
      real function flx3(s1,s2,s3)
      implicit none

      real, intent(in) :: s1,s2,s3

      ! 3rd-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx3 = (  (-1.0/6.0)*s1  &
               +( 5.0/6.0)*s2  &
               +( 2.0/6.0)*s3  )

      end function flx3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx4
      real function flx4(s1,s2,s3,s4)
      implicit none

      real, intent(in) :: s1,s2,s3,s4

      ! 4th-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx4 = (  (7.0/12.0)*(s3+s2)  &
               -(1.0/12.0)*(s4+s1)  )

      end function flx4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx5
      real function flx5(s1,s2,s3,s4,s5)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5

      ! 5th-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx5 = (  (  2.0/60.0)*s1  &
               +(-13.0/60.0)*s2  &
               +( 47.0/60.0)*s3  &
               +( 27.0/60.0)*s4  &
               +( -3.0/60.0)*s5  )

      end function flx5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx6
      real function flx6(s1,s2,s3,s4,s5,s6)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6

      ! 6th-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx6 = (  (37.0/60.0)*(s4+s3)  &
               +(-8.0/60.0)*(s5+s2)  &
               +( 1.0/60.0)*(s6+s1)  )

      end function flx6

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: weno5
      real function weno5(s1,s2,s3,s4,s5,weps)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5
      double precision, intent(in) :: weps

      double precision :: b1,b2,b3
      double precision :: w1,w2,w3

      ! 5th-order weighted essentially non-oscillatory (weno)
      ! Jiang and Shu, 1996, JCP

      b1 = (13.0/12.0)*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
      b2 = (13.0/12.0)*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
      b3 = (13.0/12.0)*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

      if( siform.eq.1 )then
        ! original WENO (eg, Jiang and Shu, 1996, JCP)
        w1 = 0.1/(weps+b1)**2
        w2 = 0.6/(weps+b2)**2
        w3 = 0.3/(weps+b3)**2
      elseif( siform.eq.2 )then
        ! improved smoothness indicators (Borges et al, 2008, JCP)
        w1 = 0.1*(1.0+min(1.0d30,abs(b1-b3)/(b1+weps))**2)
        w2 = 0.6*(1.0+min(1.0d30,abs(b1-b3)/(b2+weps))**2)
        w3 = 0.3*(1.0+min(1.0d30,abs(b1-b3)/(b3+weps))**2)
      endif

      weno5 = ( w1*( ( 2.0/6.0)*s1 + (-7.0/6.0)*s2 + (11.0/6.0)*s3 )  &
               +w2*( (-1.0/6.0)*s2 + ( 5.0/6.0)*s3 + ( 2.0/6.0)*s4 )  &
               +w3*( ( 2.0/6.0)*s3 + ( 5.0/6.0)*s4 + (-1.0/6.0)*s5 )  &
              )/( w1+w2+w3 )

      end function weno5


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: upstrpd
      real function upstrpd(s1,s2,s3,weps)
      implicit none

      real, intent(in) :: s1,s2,s3
      double precision, intent(in) :: weps

      real :: dd,rr,phi

      ! Positive-definite upstream scheme of Beets & Koren (1996, 
      ! Department of Numerical Mathematics Rep. NM-R9601, Utrecht 
      ! University, 24 pp).

      dd = s2-s1
      rr = (s3-s2)/(sign(sngl(weps),dd)+dd)
      phi = max(0.0,min(rr,min( (1.0/6.0)+(2.0/6.0)*rr , 1.0 ) ) )
      upstrpd = s2 + phi*(s2-s1)

      end function upstrpd

!-----------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------------------------------

  END MODULE interp_routines
