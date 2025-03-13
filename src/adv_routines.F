  MODULE adv_routines

  implicit none

  private
  public :: hadv_weno3,hadv_weno5,hadv_weno7,hadv_weno9
  public :: vadv_weno3,vadv_weno5,vadv_weno7,vadv_weno9
  public :: hadv_flx2,hadv_flx3,hadv_flx4,hadv_flx5,hadv_flx6,hadv_flx7,hadv_flx8,hadv_flx9,hadv_flx10
  public :: vadv_flx2,vadv_flx3,vadv_flx4,vadv_flx5,vadv_flx6,vadv_flx7,vadv_flx8,vadv_flx9,vadv_flx10
  public :: advsaxi,advuaxi,advvaxi,advwaxi,vadv_axiu
  public :: movesfc,wsub,zsgrad


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

    subroutine hadv_weno3(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2
    real :: s1,s2,s3,s4,s5
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,bmax,doit)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          dumx(i,j,k) = rru(i,j,k)*weno3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
        else
          dumx(i,j,k) = rru(i,j,k)*weno3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          dumy(i,j,k) = rrv(i,j,k)*weno3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
        else
          dumy(i,j,k) = rrv(i,j,k)*weno3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
        endif
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno3(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno3(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
        else
          dumy(i,j,k) = vbar*weno3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
        endif
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
        endif
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno3(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),weps)
        else
          dumy(i,j,k) = vbar*weno3(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),weps)
        endif
      enddo
      enddo
    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
        else
          dumy(i,j,k) = vbar*weno3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_weno3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_weno5(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2
    real :: s1,s2,s3,s4,s5
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

  pdcheck:  &
  if( pdef.eq.1 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,bmax,doit)
    DO k=1,nk
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1

        if(rru(i,j,k).ge.0.0)then
          s1=a(i-3,j,k)
          s2=a(i-2,j,k)
          s3=a(i-1,j,k)
          s4=a(i  ,j,k)
          s5=a(i+1,j,k)
        else
          s1=a(i+2,j,k)
          s2=a(i+1,j,k)
          s3=a(i  ,j,k)
          s4=a(i-1,j,k)
          s5=a(i-2,j,k)
        endif
        doit = .true.
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif


        IF(doit)THEN
          dumx(i,j,k) = rru(i,j,k)*weno5(s1,s2,s3,s4,s5,weps)
        ELSE
          dumx(i,j,k) = 0.0
        ENDIF

        if(rrv(i,j,k).ge.0.0)then
          s1=a(i,j-3,k)
          s2=a(i,j-2,k)
          s3=a(i,j-1,k)
          s4=a(i,j  ,k)
          s5=a(i,j+1,k)
        else
          s1=a(i,j+2,k)
          s2=a(i,j+1,k)
          s3=a(i,j  ,k)
          s4=a(i,j-1,k)
          s5=a(i,j-2,k)
        endif
        doit = .true.
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif


        IF(doit)THEN
          dumy(i,j,k) = rrv(i,j,k)*weno5(s1,s2,s3,s4,s5,weps)
        ELSE
          dumy(i,j,k) = 0.0
        ENDIF

      enddo
      enddo
    ENDDO

  else  pdcheck

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1

        if(rru(i,j,k).ge.0.0)then
          dumx(i,j,k) = rru(i,j,k)*weno5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),weps)
        else
          dumx(i,j,k) = rru(i,j,k)*weno5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),weps)
        endif

        if(rrv(i,j,k).ge.0.0)then
          dumy(i,j,k) = rrv(i,j,k)*weno5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),weps)
        else
          dumy(i,j,k) = rrv(i,j,k)*weno5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),weps)
        endif

      enddo
      enddo
    ENDDO

  endif  pdcheck

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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=1,nj+1
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno5(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno5(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
        endif
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),weps)
        else
          dumy(i,j,k) = vbar*weno5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),weps)
        endif
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),weps)
        endif
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno5(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),weps)
        else
          dumy(i,j,k) = vbar*weno5(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
        endif
      enddo
      enddo
    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),weps)
        endif
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),weps)
        else
          dumy(i,j,k) = vbar*weno5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),weps)
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_weno5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_weno7(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,s6,s7,bmax,doit)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          s1=a(i-4,j,k)
          s2=a(i-3,j,k)
          s3=a(i-2,j,k)
          s4=a(i-1,j,k)
          s5=a(i  ,j,k)
          s6=a(i+1,j,k)
          s7=a(i+2,j,k)
        else
          s1=a(i+3,j,k)
          s2=a(i+2,j,k)
          s3=a(i+1,j,k)
          s4=a(i  ,j,k)
          s5=a(i-1,j,k)
          s6=a(i-2,j,k)
          s7=a(i-3,j,k)
        endif
        doit = .true.
        IF( pdef.eq.1 )THEN
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5),abs(s6),abs(s7))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif
        ENDIF
        IF(doit)THEN
          dumx(i,j,k) = rru(i,j,k)*weno7(s1,s2,s3,s4,s5,s6,s7,weps)
        ELSE
          dumx(i,j,k) = 0.0
        ENDIF
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          s1=a(i,j-4,k)
          s2=a(i,j-3,k)
          s3=a(i,j-2,k)
          s4=a(i,j-1,k)
          s5=a(i,j  ,k)
          s6=a(i,j+1,k)
          s7=a(i,j+2,k)
        else
          s1=a(i,j+3,k)
          s2=a(i,j+2,k)
          s3=a(i,j+1,k)
          s4=a(i,j  ,k)
          s5=a(i,j-1,k)
          s6=a(i,j-2,k)
          s7=a(i,j-3,k)
        endif
        doit = .true.
        IF( pdef.eq.1 )THEN
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5),abs(s6),abs(s7))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif
        ENDIF
        IF(doit)THEN
          dumy(i,j,k) = rrv(i,j,k)*weno7(s1,s2,s3,s4,s5,s6,s7,weps)
        ELSE
          dumy(i,j,k) = 0.0
        ENDIF
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno7(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno7(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno7(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),weps)
        else
          dumy(i,j,k) = vbar*weno7(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),weps)
        endif
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno7(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno7(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),weps)
        endif
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno7(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),weps)
        else
          dumy(i,j,k) = vbar*weno7(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),weps)
        endif
      enddo
      enddo
    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno7(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno7(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno7(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),weps)
        else
          dumy(i,j,k) = vbar*weno7(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),weps)
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_weno7

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_weno9(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,s6,s7,s8,s9,bmax,doit)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          s1=a(i-5,j,k)
          s2=a(i-4,j,k)
          s3=a(i-3,j,k)
          s4=a(i-2,j,k)
          s5=a(i-1,j,k)
          s6=a(i  ,j,k)
          s7=a(i+1,j,k)
          s8=a(i+2,j,k)
          s9=a(i+3,j,k)
        else
          s1=a(i+4,j,k)
          s2=a(i+3,j,k)
          s3=a(i+2,j,k)
          s4=a(i+1,j,k)
          s5=a(i  ,j,k)
          s6=a(i-1,j,k)
          s7=a(i-2,j,k)
          s8=a(i-3,j,k)
          s9=a(i-4,j,k)
        endif
        doit = .true.
        IF( pdef.eq.1 )THEN
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5),abs(s6),abs(s7),abs(s8),abs(s9))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif
        ENDIF
        IF(doit)THEN
          dumx(i,j,k) = rru(i,j,k)*weno9(s1,s2,s3,s4,s5,s6,s7,s8,s9,weps)
        ELSE
          dumx(i,j,k) = 0.0
        ENDIF
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          s1=a(i,j-5,k)
          s2=a(i,j-4,k)
          s3=a(i,j-3,k)
          s4=a(i,j-2,k)
          s5=a(i,j-1,k)
          s6=a(i,j  ,k)
          s7=a(i,j+1,k)
          s8=a(i,j+2,k)
          s9=a(i,j+3,k)
        else
          s1=a(i,j+4,k)
          s2=a(i,j+3,k)
          s3=a(i,j+2,k)
          s4=a(i,j+1,k)
          s5=a(i,j  ,k)
          s6=a(i,j-1,k)
          s7=a(i,j-2,k)
          s8=a(i,j-3,k)
          s9=a(i,j-4,k)
        endif
        doit = .true.
        IF( pdef.eq.1 )THEN
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5),abs(s6),abs(s7),abs(s8),abs(s9))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif
        ENDIF
        IF(doit)THEN
          dumy(i,j,k) = rrv(i,j,k)*weno9(s1,s2,s3,s4,s5,s6,s7,s8,s9,weps)
        ELSE
          dumy(i,j,k) = 0.0
        ENDIF
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno9(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno9(a(i+5,j,k),a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno9(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),weps)
        else
          dumy(i,j,k) = vbar*weno9(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),a(i,j-4,k),weps)
        endif
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
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno9(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno9(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),a(i-4,j,k),weps)
        endif
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno9(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k),weps)
        else
          dumy(i,j,k) = vbar*weno9(a(i,j+5,k),a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),weps)
        endif
      enddo
      enddo
    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*weno9(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),weps)
        else
          dumx(i,j,k) = ubar*weno9(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),a(i-4,j,k),weps)
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*weno9(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),weps)
        else
          dumy(i,j,k) = vbar*weno9(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),a(i,j-4,k),weps)
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_weno9

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx2(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        dumx(i,j,k) = rru(i,j,k)*0.5*(a(i-1,j,k)+a(i  ,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        dumy(i,j,k) = rrv(i,j,k)*0.5*(a(i,j-1,k)+a(i,j  ,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        dumx(i,j,k) = ubar*0.5*(a(i  ,j,k)+a(i+1,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i  ,j,k))
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        dumy(i,j,k) = vbar*0.5*(a(i,j  ,k)+a(i,j+1,k))
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i  ,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx3(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          dumx(i,j,k) = rru(i,j,k)*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
        else
          dumx(i,j,k) = rru(i,j,k)*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          dumy(i,j,k) = rrv(i,j,k)*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
        else
          dumy(i,j,k) = rrv(i,j,k)*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx3(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
        else
          dumx(i,j,k) = ubar*flx3(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
        else
          dumy(i,j,k) = vbar*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
        else
          dumx(i,j,k) = ubar*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
        endif
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx3(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
        else
          dumy(i,j,k) = vbar*flx3(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k))
        endif
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
        else
          dumx(i,j,k) = ubar*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
        else
          dumy(i,j,k) = vbar*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx4(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        dumx(i,j,k) = rru(i,j,k)*flx4(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        dumy(i,j,k) = rrv(i,j,k)*flx4(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        dumx(i,j,k) = ubar*flx4(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        dumy(i,j,k) = vbar*flx4(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        dumx(i,j,k) = ubar*flx4(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        dumy(i,j,k) = vbar*flx4(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        dumx(i,j,k) = ubar*flx4(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        dumy(i,j,k) = vbar*flx4(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx5(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          dumx(i,j,k) = rru(i,j,k)*flx5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
        else
          dumx(i,j,k) = rru(i,j,k)*flx5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k))
        endif
        if(rrv(i,j,k).ge.0.0)then
          dumy(i,j,k) = rrv(i,j,k)*flx5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
        else
          dumy(i,j,k) = rrv(i,j,k)*flx5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj+1
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx5(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
        else
          dumx(i,j,k) = ubar*flx5(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
        endif
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
        else
          dumy(i,j,k) = vbar*flx5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
        else
          dumx(i,j,k) = ubar*flx5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k))
        endif
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx5(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
        else
          dumy(i,j,k) = vbar*flx5(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
        endif
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
        else
          dumx(i,j,k) = ubar*flx5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k))
        endif
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
        else
          dumy(i,j,k) = vbar*flx5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k))
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx6(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        dumx(i,j,k) = rru(i,j,k)*flx6(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        dumy(i,j,k) = rrv(i,j,k)*flx6(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        dumx(i,j,k) = ubar*flx6(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        dumy(i,j,k) = vbar*flx6(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        dumx(i,j,k) = ubar*flx6(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        dumy(i,j,k) = vbar*flx6(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        dumx(i,j,k) = ubar*flx6(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        dumy(i,j,k) = vbar*flx6(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx6

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx7(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          dumx(i,j,k) = rru(i,j,k)*flx7(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
        else
          dumx(i,j,k) = rru(i,j,k)*flx7(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          dumy(i,j,k) = rrv(i,j,k)*flx7(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
        else
          dumy(i,j,k) = rrv(i,j,k)*flx7(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx7(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
        else
          dumx(i,j,k) = ubar*flx7(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx7(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
        else
          dumy(i,j,k) = vbar*flx7(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx7(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
        else
          dumx(i,j,k) = ubar*flx7(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k))
        endif
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx7(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
        else
          dumy(i,j,k) = vbar*flx7(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k))
        endif
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx7(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k))
        else
          dumx(i,j,k) = ubar*flx7(a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx7(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k))
        else
          dumy(i,j,k) = vbar*flx7(a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k))
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx7

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx8(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        dumx(i,j,k) = rru(i,j,k)*flx8(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        dumy(i,j,k) = rrv(i,j,k)*flx8(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        dumx(i,j,k) = ubar*flx8(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        dumy(i,j,k) = vbar*flx8(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        dumx(i,j,k) = ubar*flx8(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        dumy(i,j,k) = vbar*flx8(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k))
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        dumx(i,j,k) = ubar*flx8(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        dumy(i,j,k) = vbar*flx8(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx8

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx9(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(rru(i,j,k).ge.0.0)then
          dumx(i,j,k) = rru(i,j,k)*flx9(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
        else
          dumx(i,j,k) = rru(i,j,k)*flx9(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),a(i-4,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          dumy(i,j,k) = rrv(i,j,k)*flx9(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
        else
          dumy(i,j,k) = rrv(i,j,k)*flx9(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),a(i,j-4,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx9(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k))
        else
          dumx(i,j,k) = ubar*flx9(a(i+5,j,k),a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx9(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
        else
          dumy(i,j,k) = vbar*flx9(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),a(i,j-4,k))
        endif
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx9(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
        else
          dumx(i,j,k) = ubar*flx9(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),a(i-4,j,k))
        endif
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx9(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k))
        else
          dumy(i,j,k) = vbar*flx9(a(i,j+5,k),a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k))
        endif
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dumx(i,j,k) = ubar*flx9(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k))
        else
          dumx(i,j,k) = ubar*flx9(a(i+4,j,k),a(i+3,j,k),a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k),a(i-3,j,k),a(i-4,j,k))
        endif
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dumy(i,j,k) = vbar*flx9(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k))
        else
          dumy(i,j,k) = vbar*flx9(a(i,j+4,k),a(i,j+3,k),a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k),a(i,j-3,k),a(i,j-4,k))
        endif
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx9

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine hadv_flx10(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        dumx(i,j,k) = rru(i,j,k)*flx10(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        dumy(i,j,k) = rrv(i,j,k)*flx10(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=(i1-1),i2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        dumx(i,j,k) = ubar*flx10(a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k),a(i+5,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        dumy(i,j,k) = vbar*flx10(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k))
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
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk
      do j=j1,j2
      !dir$ vector always
      do i=1,ni+1
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        dumx(i,j,k) = ubar*flx10(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k))
      enddo
      enddo
      do j=(j1-1),j2
      !dir$ vector always
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        dumy(i,j,k) = vbar*flx10(a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k),a(i,j+5,k))
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        dumx(i,j,k) = ubar*flx10(a(i-5,j,k),a(i-4,j,k),a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),a(i+2,j,k),a(i+3,j,k),a(i+4,j,k))
      enddo
      enddo
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        dumy(i,j,k) = vbar*flx10(a(i,j-5,k),a(i,j-4,k),a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),a(i,j+2,k),a(i,j+3,k),a(i,j+4,k))
      enddo
      enddo
    ENDDO

  ENDIF

    end subroutine hadv_flx10

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_weno3(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,bmax,doit)
    DO j=1,nj
      do k=3,nk-1
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*weno3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=3,nk-1
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*weno3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=3,nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*weno3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=2,nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),weps)
        endif
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

    end subroutine vadv_weno3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_weno5(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,pdefweno,weps)
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
      do k=3,nk-1
      do j=1,nj
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
          dumz(i,j,k) = rrw(i,j,k)*weno5(s1,s2,s3,s4,s5,weps)
        ELSE
          dumz(i,j,k) = 0.0
        ENDIF
      enddo
      enddo
      enddo

      call     vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_35_b(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

  !--------------------------------------!
  else  pdefs

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k)
      do k=4,nk-2
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo
      enddo

      call     vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

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
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
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
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
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
          dumz(i,j,k) = wbar*weno5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

    ENDDO

!ccccccccccccccccccccccccccccccccccc

  ENDIF

    IF( stag.ne.1 )THEN
      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

    ENDIF

    end subroutine vadv_weno5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_weno7(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,s1,s2,s3,s4,s5,s6,s7,bmax,doit)
    DO j=1,nj
      do k=5,nk-3
      !dir$ vector always
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          s1=a(i,j,k-4)
          s2=a(i,j,k-3)
          s3=a(i,j,k-2)
          s4=a(i,j,k-1)
          s5=a(i,j,k  )
          s6=a(i,j,k+1)
          s7=a(i,j,k+2)
        else
          s1=a(i,j,k+3)
          s2=a(i,j,k+2)
          s3=a(i,j,k+1)
          s4=a(i,j,k  )
          s5=a(i,j,k-1)
          s6=a(i,j,k-2)
          s7=a(i,j,k-3)
        endif
        doit = .true.
        IF( pdef.eq.1 )THEN
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5),abs(s6),abs(s7))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif
        ENDIF
        IF(doit)THEN
          dumz(i,j,k) = rrw(i,j,k)*weno7(s1,s2,s3,s4,s5,s6,s7,weps)
        ELSE
          dumz(i,j,k) = 0.0
        ENDIF
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
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=5,nk-3
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
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
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=5,nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=4,nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_57(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

    end subroutine vadv_weno7

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_weno9(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,s1,s2,s3,s4,s5,s6,s7,s8,s9,bmax,doit)
    DO j=1,nj
      do k=6,nk-4
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          s1=a(i,j,k-5)
          s2=a(i,j,k-4)
          s3=a(i,j,k-3)
          s4=a(i,j,k-2)
          s5=a(i,j,k-1)
          s6=a(i,j,k  )
          s7=a(i,j,k+1)
          s8=a(i,j,k+2)
          s9=a(i,j,k+3)
        else
          s1=a(i,j,k+4)
          s2=a(i,j,k+3)
          s3=a(i,j,k+2)
          s4=a(i,j,k+1)
          s5=a(i,j,k  )
          s6=a(i,j,k-1)
          s7=a(i,j,k-2)
          s8=a(i,j,k-3)
          s9=a(i,j,k-4)
        endif
        doit = .true.
        IF( pdef.eq.1 )THEN
          bmax = max(abs(s1),abs(s2),abs(s3),abs(s4),abs(s5),abs(s6),abs(s7),abs(s8),abs(s9))
          if( bmax.gt.bmin )then
            doit = .true.
          else
            doit = .false.
          endif
        ENDIF
        IF(doit)THEN
          dumz(i,j,k) = rrw(i,j,k)*weno9(s1,s2,s3,s4,s5,s6,s7,s8,s9,weps)
        ELSE
          dumz(i,j,k) = 0.0
        ENDIF
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
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=6,nk-4
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
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
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=6,nk-4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=5,nk-4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+5),a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_57(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

      call     vadv_lwr_weno_79(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)

    end subroutine vadv_weno9

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx2(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=2,nk
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=2,nk
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=2,nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=1,nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
      enddo
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_flx2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx3(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=3,nk-1
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=3,nk-1
      !dir$ vector always
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=3,nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=2,nk-1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ))
        endif
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx4(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=3,(nk-1)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=3,(nk-1)
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=3,(nk-1)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=2,(nk-1)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_2(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx5(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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
      do k=4,nk-2
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo
      enddo

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
      !$omp private(i,j,k,wbar)
      do k=4,nk-2
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo
      enddo

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
      !$omp private(i,j,k,wbar)
      do k=4,nk-2
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo
      enddo

  ELSEIF( stag.eq.4 )THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,wbar)
      do k=3,nk-2
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo
      enddo

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx6(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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
        dumz(i,j,k) = rrw(i,j,k)*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=4,nk-2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=4,nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
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
        dumz(i,j,k) = wbar*flx6(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_4(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_2(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx6

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx7(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=1,nj

      do k=5,nk-3
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
      enddo

    ENDDO

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      do k=5,nk-3
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      do k=5,nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=4,nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_57(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx7

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx8(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=5,nk-3
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=5,nk-3
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=5,nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=4,nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_6(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_4(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_2(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx8

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx9(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=1,nj

      do k=6,nk-4
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo
      enddo

    ENDDO

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      do k=6,nk-4
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      do k=6,nk-4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=5,nk-4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+5),a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_57(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_79(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx9

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_flx10(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=6,nk-4
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*flx10(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=6,nk-4
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx10(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=6,nk-4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx10(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=5,nk-4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*flx10(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4),a(i,j,k+5))
      enddo
      enddo

    ENDDO

  ENDIF

      ! Get lower-order fluxes near bottom/top of domain:

      call     vadv_lwr_flx_8(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_6(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_4(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

      call     vadv_lwr_flx_2(stag,ix,jy,kz,c1,c2,rrw,dumz,a)

    end subroutine vadv_flx10

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_weno_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,s1,s2,s3,s4,s5,s6,s7,s8,s9,bmax,doit)
    DO j=1,nj

      k = 2
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
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

      k = 2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
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

      k = 2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      k = 1
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),weps)
        endif
      enddo

      k = nk
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
        endif
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_weno_23

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_weno_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

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

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

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

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo


  ELSEIF( stag.eq.4 )THEN

      k = 2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),weps)
        endif
      enddo
      enddo

  ENDIF

    end subroutine vadv_lwr_weno_35

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_weno_35_b(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).lt.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

    end subroutine vadv_lwr_weno_35_b

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_weno_57(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

      k = 4
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

      k = nk-2
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

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

      k = 4
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

      k = nk-2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

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

      k = 4
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

      k = nk-2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo


  ELSEIF( stag.eq.4 )THEN

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

      k = nk-2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

  ENDIF

    end subroutine vadv_lwr_weno_57

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_weno_79(stag,ix,jy,kz,c1,c2,rrw,dumz,a,pdef,weps)
    use input
    implicit none

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: pdef
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2
    real :: wbar,cc1,cc2
    real :: s1,s2,s3,s4,s5,s6,s7,s8,s9
    real :: bmin,bmax
    logical :: doit

    bmin = weps

  IF( stag.eq.1 )THEN

      k = 5
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
      enddo
      enddo

      k = nk-3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dumz(i,j,k) = rrw(i,j,k)*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = rrw(i,j,k)*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

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

      k = 5
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
      enddo
      enddo

      k = nk-3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

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

      k = 5
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
      enddo
      enddo

      k = nk-3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo


  ELSEIF( stag.eq.4 )THEN

      k = 4
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+5),a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

      k = nk-3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

  ENDIF

    end subroutine vadv_lwr_weno_79

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_23(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      k = 2
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

      k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        endif
      enddo
      enddo

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

      k = 2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

      k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo
      enddo


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

      k = 2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

      k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        endif
      enddo
      enddo


  ELSEIF( stag.eq.4 )THEN

      k = 1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ))
        endif
      enddo
      enddo

      k = nk
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
        endif
      enddo
      enddo

  ENDIF

    end subroutine vadv_lwr_flx_23

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_35(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

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

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

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

      k = 3
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=j1,j2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

  ELSEIF( stag.eq.4 )THEN

      k = 2
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo

      k = nk-1
      !$omp parallel do default(shared)   &
      !$omp private(i,j,cc1,cc2,wbar)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ))
        endif
      enddo
      enddo

  ENDIF

    end subroutine vadv_lwr_flx_35

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_57(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      k = 4
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

      k = nk-2
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

    ENDDO

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      k = 4
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

      k = nk-2
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      k = 4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

      k = nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      k = 3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

      k = nk-2
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_flx_57

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_79(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  IF( stag.eq.1 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=1,nj

      k = 5
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo

      k = nk-3
      !dir$ vector always
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

    ENDDO

  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      k = 5
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo

      k = nk-3
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
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

      k = 5
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo

      k = nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

    ENDDO


  !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      k = 4
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+5),a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo

      k = nk-3
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_flx_79

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_2(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=1,nk,(nk-1)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
      enddo
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_flx_2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_4(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=2,(nk-1),(nk-3)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_flx_4

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_6(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=4,nk-2,(nk-6)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=4,nk-2,(nk-6)
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=4,nk-2,(nk-6)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=3,nk-2,(nk-5)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_flx_6

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine vadv_lwr_flx_8(stag,ix,jy,kz,c1,c2,rrw,dumz,a)
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

      do k=5,nk-3,(nk-8)
      !dir$ vector always
      do i=1,ni
        dumz(i,j,k) = rrw(i,j,k)*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=1,nj

      do k=5,nk-3,(nk-8)
      !dir$ vector always
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
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
!$omp private(i,j,k,wbar,cc1,cc2)
    DO j=j1,j2

      do k=5,nk-3,(nk-8)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo

    ENDDO


  ELSEIF( stag.eq.4 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
    DO j=1,nj

      do k=4,nk-3,(nk-7)
      !dir$ vector always
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
      enddo
      enddo

    ENDDO

  ENDIF

    end subroutine vadv_lwr_flx_8

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

      !dir$ attributes forceinline :: flx7
      real function flx7(s1,s2,s3,s4,s5,s6,s7)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6,s7

      ! 7th-order flux

      flx7 = (  (  -3.0/420.0)*s1  &
               +(  25.0/420.0)*s2  &
               +(-101.0/420.0)*s3  &
               +( 319.0/420.0)*s4  &
               +( 214.0/420.0)*s5  &
               +( -38.0/420.0)*s6  &
               +(   4.0/420.0)*s7  )

      end function flx7

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx8
      real function flx8(s1,s2,s3,s4,s5,s6,s7,s8)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6,s7,s8

      ! 8th-order flux

      flx8 = (  ( 533.0/840.0)*(s5+s4)  &
               +(-139.0/840.0)*(s6+s3)  &
               +(  29.0/840.0)*(s7+s2)  &
               +(  -3.0/840.0)*(s8+s1)  )

      end function flx8

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx9
      real function flx9(s1,s2,s3,s4,s5,s6,s7,s8,s9)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6,s7,s8,s9

      ! 9th-order flux

      flx9 = (  (   12.0/7560.0)*s1  &
               +( -123.0/7560.0)*s2  &
               +(  597.0/7560.0)*s3  &
               +(-1923.0/7560.0)*s4  &
               +( 5637.0/7560.0)*s5  &
               +( 4125.0/7560.0)*s6  &
               +( -915.0/7560.0)*s7  &
               +(  165.0/7560.0)*s8  &
               +(  -15.0/7560.0)*s9  )

      end function flx9


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx10
      real function flx10(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10

      ! 10th-order flux

      flx10 = (  ( 4881.0/7560.0)*( s6+s5)  &
                +(-1419.0/7560.0)*( s7+s4)  &
                +(  381.0/7560.0)*( s8+s3)  &
                +(  -69.0/7560.0)*( s9+s2)  &
                +(   +6.0/7560.0)*(s10+s1)  )

      end function flx10

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: weno3
      real function weno3(s1,s2,s3,weps)
      implicit none

      real, intent(in) :: s1,s2,s3
      double precision, intent(in) :: weps

      double precision :: b1,b2
      double precision :: w1,w2

      ! 3rd-order weighted essentially non-oscillatory (weno)
      ! Jiang and Shu, 1996, JCP

      b1 = (s1-s2)**2
      b2 = (s2-s3)**2

      if( siform.eq.1 )then
        ! original WENO (eg, Jiang and Shu, 1996, JCP)
        w1 = (1.0/3.0)/(weps+b1)**2
        w2 = (2.0/3.0)/(weps+b2)**2
      elseif( siform.eq.2 )then
        ! improved smoothness indicators (Borges et al, 2008, JCP)
        w1 = (1.0/3.0)*(1.0+min(1.0d30,abs(b1-b2)/(b1+weps))**2)
        w2 = (2.0/3.0)*(1.0+min(1.0d30,abs(b1-b2)/(b2+weps))**2)
      endif

      weno3 = ( w1*( (-1.0/2.0)*s1 + ( 3.0/2.0)*s2 )  &
               +w2*( ( 1.0/2.0)*s2 + ( 1.0/2.0)*s3 )  &
              )/( w1+w2 )

      end function weno3

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

      !dir$ attributes forceinline :: weno7
      real function weno7(s1,s2,s3,s4,s5,s6,s7,weps)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6,s7
      double precision, intent(in) :: weps

      double precision :: b0,b1,b2,b3
      double precision :: w0,w1,w2,w3

      ! 7th-order weighted essentially non-oscillatory (weno)
      ! Balsara and Shu, 2000, JCP

      b0 = s1*(   547.0*s1 -  3882.0*s2 + 4642.0*s3 - 1854.0*s4 )  &
          +s2*(  7043.0*s2 - 17246.0*s3 + 7042.0*s4 )              &
          +s3*( 11003.0*s3 -  9402.0*s4 )                          &
          +s4*s4*2107.0
      b1 = s2*(   267.0*s2 -  1642.0*s3 + 1602.0*s4 - 494.0*s5 )   &
          +s3*(  2843.0*s3 -  5966.0*s4 + 1922.0*s5 )              &
          +s4*(  3443.0*s4 -  2522.0*s5 )                          &
          +s5*s5*547.0
      b2 = s3*(   547.0*s3 -  2522.0*s4 + 1922.0*s5 - 494.0*s6 )   &
          +s4*(  3443.0*s4 -  5966.0*s5 + 1602.0*s6 )              &
          +s5*(  2843.0*s5 -  1642.0*s6 )                          &
          +s6*s6*267.0
      b3 = s4*(  2107.0*s4 -  9402.0*s5 + 7042.0*s6 - 1854.0*s7 )  &
          +s5*( 11003.0*s5 - 17246.0*s6 + 4642.0*s7 )              &
          +s6*(  7043.0*s6 -  3882.0*s7 )                          &
          +s7*s7*547.0

      if( siform.eq.1 )then
        ! original WENO (eg, Balsara and Shu, 2000, JCP)
        w0 = ( 1.0/35.0)*min(1.0d30,(weps+b0)**(-2))
        w1 = (12.0/35.0)*min(1.0d30,(weps+b1)**(-2))
        w2 = (18.0/35.0)*min(1.0d30,(weps+b2)**(-2))
        w3 = ( 4.0/35.0)*min(1.0d30,(weps+b3)**(-2))
      elseif( siform.eq.2 )then
        ! improved smoothness indicators (Borges et al, 2008, JCP)
        w0 = ( 1.0/35.0)*(1.0+min(1.0d30,abs(b0-b3)/(b0+weps))**2)
        w1 = (12.0/35.0)*(1.0+min(1.0d30,abs(b0-b3)/(b1+weps))**2)
        w2 = (18.0/35.0)*(1.0+min(1.0d30,abs(b0-b3)/(b2+weps))**2)
        w3 = ( 4.0/35.0)*(1.0+min(1.0d30,abs(b0-b3)/(b3+weps))**2)
      endif

      weno7 = ( w0*( -(3.0/12.0)*s1 +(13.0/12.0)*s2 -(23.0/12.0)*s3 +(25.0/12.0)*s4 )  &
               +w1*( +(1.0/12.0)*s2 -( 5.0/12.0)*s3 +(13.0/12.0)*s4 +( 3.0/12.0)*s5 )  &
               +w2*( -(1.0/12.0)*s3 +( 7.0/12.0)*s4 +( 7.0/12.0)*s5 -( 1.0/12.0)*s6 )  &
               +w3*( +(3.0/12.0)*s4 +(13.0/12.0)*s5 -( 5.0/12.0)*s6 +( 1.0/12.0)*s7 )  &
              )/(w0+w1+w2+w3)

      end function weno7

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: weno9
      real function weno9(s1,s2,s3,s4,s5,s6,s7,s8,s9,weps)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6,s7,s8,s9
      double precision, intent(in) :: weps

      double precision :: b0,b1,b2,b3,b4
      double precision :: w0,w1,w2,w3,w4

      ! 9th-order weighted essentially non-oscillatory (weno)
      ! Balsara and Shu, 2000, JCP

      b0 = s1*(   22658.0*s1 -  208501.0*s2 + 364863.0*s3 - 288007.0*s4 + 86329.0*s5 )  &
          +s2*(  482963.0*s2 - 1704396.0*s3 + 1358458.0*s4 - 411487.0*s5 )              &
          +s3*( 1521393.0*s3 - 2462076.0*s4 + 758823.0*s5 )                             &
          +s4*( 1020563.0*s4 -  649501.0*s5 )                                           &
          +s5*s5*107918.0
      b1 = s2*(    6908.0*s2 -   60871.0*s3 +  99213.0*s4 -  70237.0*s5 + 18079.0*s6 )  &
          +s3*(  138563.0*s3 -  464976.0*s4 + 337018.0*s5 -  88297.0*s6 )               &
          +s4*(  406293.0*s4 -  611976.0*s5 + 165153.0*s6 )                             &
          +s5*(  242723.0*s5 -  140251.0*s6 )                                           &
          +s6*s6*22658.0
      b2 = s3*(    6908.0*s3 -   51001.0*s4 +  67923.0*s5 - 38947.0*s6 + 8209.0*s7 )    &
          +s4*(  104963.0*s4 -  299076.0*s5 + 179098.0*s6 - 38947.0*s7 )                &
          +s5*(  251153.0*s5 -  299076.0*s6 +  67923.0*s7 )                             &
          +s6*(  104963.0*s6 -   51001.0*s7 )                                           &
          +s7*s7*6908.0
      b3 = s4*(   22658.0*s4 -  140251.0*s5 + 165153.0*s6 -  88297.0*s7 + 18079.0*s8 )  &
          +s5*(  242723.0*s5 -  611976.0*s6 + 337018.0*s7 -  70237.0*s8 )               &
          +s6*(  406293.0*s6 -  464976.0*s7 +  99213.0*s8 )                             &
          +s7*(  138563.0*s7 -   60871.0*s8 )                                           &
          +s8*s8*6908.0
      b4 = s5*(  107918.0*s5 -  649501.0*s6 + 758823.0*s7 - 411487.0*s8 + 86329.0*s9 )  &
          +s6*( 1020563.0*s6 - 2462076.0*s7 + 1358458.0*s8 - 288007.0*s9 )              &
          +s7*( 1521393.0*s7 - 1704396.0*s8 + 364863.0*s9 )                             &
          +s8*(  482963.0*s8 -  208501.0*s9 )                                           &
          +s9*s9*22658.0

      if( siform.eq.1 )then
        ! original WENO (eg, Balsara and Shu, 2000, JCP)
        w0 = ( 1.0/126.0)*min(1.0d28,(weps+b0)**(-2))
        w1 = (20.0/126.0)*min(1.0d28,(weps+b1)**(-2))
        w2 = (60.0/126.0)*min(1.0d28,(weps+b2)**(-2))
        w3 = (40.0/126.0)*min(1.0d28,(weps+b3)**(-2))
        w4 = ( 5.0/126.0)*min(1.0d28,(weps+b4)**(-2))
      elseif( siform.eq.2 )then
        ! improved smoothness indicators (Borges et al, 2008, JCP)
        w0 = ( 1.0/126.0)*(1.0+min(1.0d30,abs(b0-b4)/(b0+weps))**2)
        w1 = (20.0/126.0)*(1.0+min(1.0d30,abs(b0-b4)/(b1+weps))**2)
        w2 = (60.0/126.0)*(1.0+min(1.0d30,abs(b0-b4)/(b2+weps))**2)
        w3 = (40.0/126.0)*(1.0+min(1.0d30,abs(b0-b4)/(b3+weps))**2)
        w4 = ( 5.0/126.0)*(1.0+min(1.0d30,abs(b0-b4)/(b4+weps))**2)
      endif

      weno9 = ( w0*( +(12.0/60.0)*s1  -(63.0/60.0)*s2 +(137.0/60.0)*s3 -(163.0/60.0)*s4 +(137.0/60.0)*s5 )  &
               +w1*( -( 3.0/60.0)*s2  +(17.0/60.0)*s3 -( 43.0/60.0)*s4 +( 77.0/60.0)*s5 +( 12.0/60.0)*s6 )  &
               +w2*( +( 2.0/60.0)*s3  -(13.0/60.0)*s4 +( 47.0/60.0)*s5 +( 27.0/60.0)*s6 -(  3.0/60.0)*s7 )  &
               +w3*( -( 3.0/60.0)*s4  +(27.0/60.0)*s5 +( 47.0/60.0)*s6 -( 13.0/60.0)*s7 +(  2.0/60.0)*s8 )  &
               +w4*( +(12.0/60.0)*s5  +(77.0/60.0)*s6 -( 43.0/60.0)*s7 +( 17.0/60.0)*s8 -(  3.0/60.0)*s9 )  &
              )/(w0+w1+w2+w3+w4)

      end function weno9

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
!
!                Advection for axisymmetric model:
!
!-----------------------------------------------------------------------------------------------


      subroutine vadv_axiu(doweno,arf1,arf2,c1,c2,rrw,a,dumz,weps)
      use input
      use constants
      implicit none

      ! vertical advection for u, axisymmetric solver:

      logical, intent(in) :: doweno
      real, intent(in), dimension(ib:ie+1) :: arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: a
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
      double precision, intent(in) :: weps
 
      integer :: i,k,i1,i2
      real :: wbar

      integer, parameter :: j = 1

!----------------------------------------------------------------
    ! u-staggered, axisymmetric grid:

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

      IF(ebc.eq.3.or.ebc.eq.4) i2 = ni

  wenozu:  &
  IF( doweno )THEN

    !-------------------------------------------------------------------
    !  begin weno section:

    IF( weno_order.eq.3 )THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=3,nk-1
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*weno3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      enddo

    ELSEIF( weno_order.eq.5 )THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=4,nk-2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      enddo

    ELSEIF( weno_order.eq.7 )THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=5,nk-3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      enddo

    ELSEIF( weno_order.eq.9 )THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=6,nk-4
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
      enddo
      enddo

    ELSE

      stop 23598

    ENDIF

    !----------------
    ! near boundaries:

    IF( weno_order.eq.9 )THEN
      k = 5
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4),weps)
        endif
      enddo
      k = nk-3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
    ENDIF

    IF( weno_order.eq.9 .or. weno_order.eq.7 )THEN
      k = 4
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*weno7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),weps)
        endif
      enddo
      k = nk-2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
    ENDIF

    IF( weno_order.eq.9 .or. weno_order.eq.7 .or. weno_order.eq.5 )THEN
      k = 3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
        endif
      enddo
      k = nk-1
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*weno5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
    ENDIF

    ! all weno orders do this:
      k = 2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dumz(i,j,k) = wbar*(c1(1,1,k)*a(i,j,k-1)+c2(1,1,k)*a(i,j,k))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
      enddo
      k = nk
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
        else
          dumz(i,j,k) = wbar*(c1(1,1,k)*a(i,j,k-1)+c2(1,1,k)*a(i,j,k))
        endif
      enddo

    ! end of weno section

  ELSE  wenozu

    !-------------------------------------------------------------------
    ! begin not-weno section:

    IF(     vadvordrv.eq.2 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=2,nk
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*0.5*(a(i,j,k-1)+a(i,j,k))
      enddo
      enddo
    ELSEIF( vadvordrv.eq.3 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=3,nk-1
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      enddo
    ELSEIF( vadvordrv.eq.4 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=3,nk-1
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo
    ELSEIF( vadvordrv.eq.5 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=4,nk-2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      enddo
    ELSEIF( vadvordrv.eq.6 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=4,nk-2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo
    ELSEIF( vadvordrv.eq.7 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=5,nk-3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
      enddo
    ELSEIF( vadvordrv.eq.8 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=5,nk-3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo
    ELSEIF( vadvordrv.eq.9 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=6,nk-4
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo
      enddo
    ELSEIF( vadvordrv.eq.10 )THEN
      !$omp parallel do default(shared)   &
      !$omp private(i,k,wbar)
      do k=6,nk-4
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx10(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3),a(i,j,k+4))
      enddo
      enddo
    ENDIF

    !--------
    ! near boundaries (odd-ordered schemes):

    IF( vadvordrv.eq.9 )THEN
      k = 5
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx9(a(i,j,k+4),a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3),a(i,j,k-4))
        endif
      enddo
      k = nk-3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx9(a(i,j,k-5),a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
    ENDIF
    IF( vadvordrv.eq.9 .or. vadvordrv.eq.7 )THEN
      k = 4
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx7(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),a(i,j,k-3))
        endif
      enddo
      k = nk-2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx7(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
    ENDIF
    IF( vadvordrv.eq.9 .or. vadvordrv.eq.7 .or. vadvordrv.eq.5 )THEN
      k = 3
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
        endif
      enddo
      k = nk-1
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx5(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
    ENDIF
    IF( vadvordrv.eq.9 .or. vadvordrv.eq.7 .or. vadvordrv.eq.5 .or. vadvordrv.eq.3 )THEN
      k = 2
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*(c1(1,1,k)*a(i,j,k-1)+c2(1,1,k)*a(i,j,k))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
      enddo
      k = nk
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
        else
          dumz(i,j,k) = wbar*(c1(1,1,k)*a(i,j,k-1)+c2(1,1,k)*a(i,j,k))
        endif
      enddo
    ENDIF


    !--------
    ! near boundaries (even-ordered schemes):

    IF( vadvordrv.eq.10 )THEN
      do k=5,(nk-3),(nk-8)
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx8(a(i,j,k-4),a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2),a(i,j,k+3))
      enddo
      enddo
    ENDIF
    IF( vadvordrv.eq.10 .or. vadvordrv.eq.8 )THEN
      do k=4,(nk-2),(nk-6)
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx6(a(i,j,k-3),a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),a(i,j,k+2))
      enddo
      enddo
    ENDIF
    IF( vadvordrv.eq.10 .or. vadvordrv.eq.8 .or. vadvordrv.eq.6 )THEN
      do k=3,(nk-1),(nk-4)
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*flx4(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
      enddo
      enddo
    ENDIF
    IF( vadvordrv.eq.10 .or. vadvordrv.eq.8 .or. vadvordrv.eq.6 .or. vadvordrv.eq.4 )THEN
      do k=2,nk,(nk-2)
      !dir$ vector always
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dumz(i,j,k) = wbar*(c1(1,1,k)*a(i,j,k-1)+c2(1,1,k)*a(i,j,k))
      enddo
      enddo
    ENDIF

    !--------

  ENDIF  wenozu

      end subroutine vadv_axiu


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advsaxi(doweno,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,rmh,gz,rgz, &
                  rho0,rr0,rf0,rrf0,advx,dum,mass,rru,s0,s,pdef,dt,weps, &
                  hadvorder,flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
      use input
      use constants
      use pdef_module
      implicit none

      ! horizontal advection for s, axisymmetric solver:

      logical, intent(in) :: doweno
      integer, intent(in) :: bflag
      double precision, intent(inout) :: bsq
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: vh,rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh,rho0,rr0,rf0,rrf0
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s0,s
      integer, intent(in) :: pdef
      real, intent(in) :: dt
      double precision, intent(in) :: weps
      integer, intent(in) :: hadvorder
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
 
      integer :: i,j,k
      real :: tem0
      double precision, dimension(nk) :: budx
      integer, dimension(4) :: reqsx

!!!      budx(k) = 0.0d0

!----------------------------------------------------------------
! Advection in r-direction

  wenochecksaxi:  &
  if(doweno)then

    if( weno_order.eq.3 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno3(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno3(s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.5 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno5(s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno5(s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.7 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno7(s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno7(s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.9 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno9(s(i-5,j,k),s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),s(i+3,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno9(s(i+4,j,k),s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k),s(i-4,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    else

      stop 9389

    endif


    ! Near boundaries ...

    if( weno_order.eq.9 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno7(s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno9(s(i+4,j,k),s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k),s(i-4,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno5(s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno7(s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 .or. weno_order.eq.5 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*upstrpd(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno5(s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),weps)
        endif
      enddo
      ENDDO
    endif

      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k) = rru(i,j,k)*upstrpd(s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),weps)
        endif
        dum(1,j,k) = 0.0
      enddo
      ENDDO

      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*upstrpd(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*s(i,j,k)
        endif
        dum(ni+1,j,k) = 0.0
      enddo
      ENDDO
      ENDIF


  else  wenochecksaxi


    if(hadvorder.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      dum(1,j,k) = 0.0
      !dir$ vector always
      do i=2,ni+1
        dum(i,j,k) = rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.3)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx3(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx3(s(i+1,j,k),s(i  ,j,k),s(i-1,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.4)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        dum(i,j,k) = rru(i,j,k)*flx4(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.5)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx5(s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx5(s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.6)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        dum(i,j,k) = rru(i,j,k)*flx6(s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.7)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx7(s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx7(s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.8)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        dum(i,j,k) = rru(i,j,k)*flx8(s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),s(i+3,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.9)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx9(s(i-5,j,k),s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),s(i+3,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx9(s(i+4,j,k),s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k),s(i-4,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadvorder.eq.10)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        dum(i,j,k) = rru(i,j,k)*flx10(s(i-5,j,k),s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),s(i+3,j,k),s(i+4,j,k))
      enddo
      enddo
      ENDDO
    else
      stop 34788
    endif


    ! Near boundaries ...

    if( hadvorder.eq.9 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx7(s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx9(s(i+4,j,k),s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k),s(i-4,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadvorder.eq.9 .or. hadvorder.eq.7 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx5(s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx7(s(i+3,j,k),s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k),s(i-3,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadvorder.eq.9 .or. hadvorder.eq.7 .or. hadvorder.eq.5 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx3(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx5(s(i+2,j,k),s(i+1,j,k),s(i  ,j,k),s(i-1,j,k),s(i-2,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadvorder.eq.9 .or. hadvorder.eq.7 .or. hadvorder.eq.5 .or. hadvorder.eq.3 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k) = rru(i,j,k)*flx3(s(i+1,j,k),s(i  ,j,k),s(i-1,j,k))
        endif
        dum(1,j,k) = 0.0
      enddo
      ENDDO
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
        DO k=1,nk
        do j=1,nj
          if(rru(i,j,k).ge.0.0)then
            dum(i,j,k) = rru(i,j,k)*flx3(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k))
          else
            dum(i,j,k) = rru(i,j,k)*s(i,j,k)
          endif
          dum(ni+1,j,k) = 0.0
        enddo
        ENDDO
      ENDIF
    endif


    if( hadvorder.eq.10 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*flx8(s(i-4,j,k),s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k),s(i+3,j,k))
      enddo
      ENDDO
    endif
    if( hadvorder.eq.10 .or. hadvorder.eq.8 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*flx6(s(i-3,j,k),s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k),s(i+2,j,k))
      enddo
      ENDDO
    endif
    if( hadvorder.eq.10 .or. hadvorder.eq.8 .or. hadvorder.eq.6 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*flx4(s(i-2,j,k),s(i-1,j,k),s(i  ,j,k),s(i+1,j,k))
      enddo
      ENDDO
    endif
    if( hadvorder.eq.10 .or. hadvorder.eq.8 .or. hadvorder.eq.6 .or. hadvorder.eq.4 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        dum(1,j,k) = 0.0
      enddo
      ENDDO
    endif

    if( hadvorder.eq.10 .or. hadvorder.eq.8 .or. hadvorder.eq.6 .or. hadvorder.eq.4 .or. hadvorder.eq.2 )then
    IF(ebc.eq.3.or.ebc.eq.4)THEN
      i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        dum(ni+1,j,k) = 0.0
      enddo
      ENDDO
    ENDIF
    endif

  endif  wenochecksaxi

    !-------------------------------------

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).le.0.0)then
          dum(i,j,k) = dum(i-1,j,k)*arh1(i-1)/arh2(i-1)
        endif
        budx(k) = budx(k)-dum(ni+1,j,k)*rvh(j)*rmh(ni+1,j,k)
      enddo
      ENDDO
      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        advx(i,j,k) = -(arh2(i)*dum(i+1,j,k)-arh1(i)*dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      ENDDO

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        i=ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
        DO k=1,nk
        do j=1,nj
          if(rru(ni+1,j,k).le.0.0)then
            advx(i,j,k) = advx(i,j,k)-s(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        ENDDO
      ENDIF


!----------------------------------------------------------------
!  Misc for r-direction

!      IF(stat_qsrc.eq.1.and.(wbc.eq.2.or.ebc.eq.2).and.bflag.eq.1)THEN
!        tem0=dt*dy*dz
!        do k=1,nk
!          bsq=bsq+budx(k)*tem0
!        enddo
!      ENDIF

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefx1(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dum,mass,s0,s,dt,flag,sw31,sw32,se31,se32,reqsx)
        call pdefx2(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dum,mass,s0,s,dt,flag,sw31,sw32,se31,se32,reqsx)
      ENDIF

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advs=time_advs+mytime()
 
      end subroutine advsaxi


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advuaxi(doweno,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,u3d,hadv)
      use input
      use constants
      implicit none

      ! horizontal advection for u, axisymmetric solver:

      logical, intent(in) :: doweno
      real, intent(in), dimension(ib:ie) :: arh1,arh2
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d
      integer, intent(in) :: hadv
 
      integer :: i,j,k,i1,i2,j1,j2,id1,id2
      real :: ubar
      double precision :: weps

!------------------------------------------------------------

      weps = 100.0*epsilon

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

      i1 = max(2,i1)
      id1 = max(1,id1)

!----------------------------------------------------------------
! Advection in r-direction

  wenocheckuaxi:  &
  if(doweno)then

    if( weno_order.eq.3 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=2,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno3(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),weps)
        else
          dum(i,j,k) = ubar*weno3(u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.5 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno5(u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),weps)
        else
          dum(i,j,k) = ubar*weno5(u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.7 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno7(u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),weps)
        else
          dum(i,j,k) = ubar*weno7(u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.9 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno9(u3d(i-4,j,k),u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),u3d(i+4,j,k),weps)
        else
          dum(i,j,k) = ubar*weno9(u3d(i+5,j,k),u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k),u3d(i-3,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    endif


    ! near boundaries ...

    if( weno_order.eq.9 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno7(u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),weps)
        else
          dum(i,j,k) = ubar*weno9(u3d(i+5,j,k),u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k),u3d(i-3,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno5(u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),weps)
        else
          dum(i,j,k) = ubar*weno7(u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 .or. weno_order.eq.5 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*upstrpd(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),weps)
        else
          dum(i,j,k) = ubar*weno5(u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),weps)
        endif
      enddo
      ENDDO
    endif

      i = 1
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*0.5*(u3d(i,j,k)+u3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*upstrpd(u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),weps)
        endif
      enddo
      ENDDO

      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*upstrpd(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),weps)
        else
          dum(i,j,k) = ubar*u3d(i+1,j,k)
        endif
      enddo
      ENDDO
      ENDIF


  else  wenocheckuaxi


    if(hadv.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*0.5*(u3d(i,j,k)+u3d(i+1,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.3)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=2,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx3(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*flx3(u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.4)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=2,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx4(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.5)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx5(u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k))
        else
          dum(i,j,k) = ubar*flx5(u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.6)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx6(u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.7)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx7(u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k))
        else
          dum(i,j,k) = ubar*flx7(u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.8)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx8(u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),u3d(i+4,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.9)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx9(u3d(i-4,j,k),u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),u3d(i+4,j,k))
        else
          dum(i,j,k) = ubar*flx9(u3d(i+5,j,k),u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k),u3d(i-3,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.10)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx10(u3d(i-4,j,k),u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),u3d(i+4,j,k),u3d(i+5,j,k))
      enddo
      enddo
      ENDDO
    else
      stop 34789
    endif


    ! Near boundaries ...

    if( hadv.eq.9 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx7(u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k))
        else
          dum(i,j,k) = ubar*flx9(u3d(i+5,j,k),u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k),u3d(i-3,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx5(u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k))
        else
          dum(i,j,k) = ubar*flx7(u3d(i+4,j,k),u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k),u3d(i-2,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 .or. hadv.eq.5 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx3(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*flx5(u3d(i+3,j,k),u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k),u3d(i-1,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 .or. hadv.eq.5 .or. hadv.eq.3 )then
      i = 1
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*0.5*(u3d(i,j,k)+u3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*flx3(u3d(i+2,j,k),u3d(i+1,j,k),u3d(i  ,j,k))
        endif
      enddo
      ENDDO
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx3(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*u3d(i+1,j,k)
        endif
      enddo
      ENDDO
      ENDIF
    endif


    if( hadv.eq.10 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx8(u3d(i-3,j,k),u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k),u3d(i+4,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx6(u3d(i-2,j,k),u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k),u3d(i+3,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*flx4(u3d(i-1,j,k),u3d(i  ,j,k),u3d(i+1,j,k),u3d(i+2,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 .or. hadv.eq.4 )then
      i = 1
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*0.5*(u3d(i+1,j,k)+u3d(i,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 .or. hadv.eq.4 .or. hadv.eq.2 )then
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=1,nk
      do j=1,nj
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*0.5*(u3d(i+1,j,k)+u3d(i,j,k))
      enddo
      ENDDO
      ENDIF
    endif


  endif  wenocheckuaxi


    !-------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
        advx(1,j,k)=0.0
        !dir$ vector always
        do i=2,i2
          advx(i,j,k) = -(arf2(i)*dum(i,j,k)-arf1(i)*dum(i-1,j,k))*rdx*uf(i)
        enddo
      enddo
      ENDDO


      IF(ebc.eq.3.or.ebc.eq.4)THEN
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
        DO k=1,nk
        do j=1,nj
          advx(ni+1,j,k)=0.0
        enddo
        ENDDO
      ENDIF

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advu=time_advu+mytime()
 
      end subroutine advuaxi


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advvaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vf,rho0,rr0,rf0,rrf0,dum,advx,v3d,rru,hadv)
      use input
      use constants
      implicit none

      ! horizontal advection for v, axisymmetric solver:

      logical, intent(in) :: doweno
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: v3d
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      integer, intent(in) :: hadv
 
      integer :: i,j,k
      double precision :: weps

      !----------------------------------------------------------!
      !  cm1r18:   NOTE:  v3d is actually M  (angular momentum)  !
      !----------------------------------------------------------!
 
!------------------------------------------------------------

      weps = 100.0*epsilon

!----------------------------------------------------------------
! Advection in r-direction

  wenocheckvaxi:  &
  if(doweno)then

    if( weno_order.eq.3 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno3(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno3(v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.5 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno5(v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno5(v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.7 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno7(v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno7(v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.9 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno9(v3d(i-5,j,k),v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),v3d(i+3,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno9(v3d(i+4,j,k),v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k),v3d(i-4,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    endif


    ! Near boundaries ...

    if( weno_order.eq.9 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno7(v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno9(v3d(i+4,j,k),v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k),v3d(i-4,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*weno5(v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno7(v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 .or. weno_order.eq.5 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*upstrpd(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*weno5(v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),weps)
        endif
      enddo
      ENDDO
    endif

      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*v3d(i-1,j,k)
        else
          dum(i,j,k) = rru(i,j,k)*upstrpd(v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),weps)
        endif
        dum(1,j,k) = 0.0
      enddo
      ENDDO

      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*upstrpd(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),weps)
        else
          dum(i,j,k) = rru(i,j,k)*v3d(i,j,k)
        endif
        dum(ni+1,j,k) = 0.0
      enddo
      ENDDO
      ENDIF


  else  wenocheckvaxi


    if(hadv.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      dum(1,j,k) = 0.0
      !dir$ vector always
      do i=2,ni+1
        dum(i,j,k) = rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.3)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx3(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx3(v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.4)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        dum(i,j,k) = rru(i,j,k)*flx4(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.5)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx5(v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx5(v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.6)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        dum(i,j,k) = rru(i,j,k)*flx6(v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.7)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx7(v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx7(v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.8)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        dum(i,j,k) = rru(i,j,k)*flx8(v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),v3d(i+3,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.9)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx9(v3d(i-5,j,k),v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),v3d(i+3,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx9(v3d(i+4,j,k),v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k),v3d(i-4,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.10)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        dum(i,j,k) = rru(i,j,k)*flx10(v3d(i-5,j,k),v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),v3d(i+3,j,k),v3d(i+4,j,k))
      enddo
      enddo
      ENDDO
    else
      stop 34790
    endif


    ! Near boundaries ...

    if( hadv.eq.9 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx7(v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx9(v3d(i+4,j,k),v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k),v3d(i-4,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx5(v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx7(v3d(i+3,j,k),v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k),v3d(i-3,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 .or. hadv.eq.5 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*flx3(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k))
        else
          dum(i,j,k) = rru(i,j,k)*flx5(v3d(i+2,j,k),v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k),v3d(i-2,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 .or. hadv.eq.5 .or. hadv.eq.3 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k) = rru(i,j,k)*v3d(i-1,j,k)
        else
          dum(i,j,k) = rru(i,j,k)*flx3(v3d(i+1,j,k),v3d(i  ,j,k),v3d(i-1,j,k))
        endif
        dum(1,j,k) = 0.0
      enddo
      ENDDO
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
        DO k=1,nk
        do j=1,nj
          if(rru(i,j,k).ge.0.0)then
            dum(i,j,k) = rru(i,j,k)*flx3(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k))
          else
            dum(i,j,k) = rru(i,j,k)*v3d(i,j,k)
          endif
          dum(ni+1,j,k) = 0.0
        enddo
        ENDDO
      ENDIF
    endif


    if( hadv.eq.10 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*flx8(v3d(i-4,j,k),v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k),v3d(i+3,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*flx6(v3d(i-3,j,k),v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k),v3d(i+2,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*flx4(v3d(i-2,j,k),v3d(i-1,j,k),v3d(i  ,j,k),v3d(i+1,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 .or. hadv.eq.4 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        dum(1,j,k) = 0.0
      enddo
      ENDDO
    endif

    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 .or. hadv.eq.4 .or. hadv.eq.2 )then
    IF(ebc.eq.3.or.ebc.eq.4)THEN
      i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
      DO k=1,nk
      do j=1,nj
        dum(i,j,k) = rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        dum(ni+1,j,k) = 0.0
      enddo
      ENDDO
    ENDIF
    endif

  endif  wenocheckvaxi

    !-------------------------------------

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).le.0.0)then
          i = ni
          dum(ni+1,j,k)=arh1(i)*arh1(i)*dum(i,j,k)/(arh2(i)*arh2(i))
        endif
      enddo
      ENDDO
      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        ! cm1r17: include centrifugal accel term here
!!!        advx(i,j,k)=-(arh2(i)*arh2(i)*dum(i+1,j,k)-arh1(i)*arh1(i)*dum(i,j,k))*rdx*uh(i)
        ! cm1r18: include centrifugal accel term + Coriolis term here
        advx(i,j,k) = -(arh2(i)*dum(i+1,j,k)-arh1(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)  &
               +0.5*fcor*(xf(i+1)*rru(i+1,j,k)-xf(i)*rru(i,j,k))*rdx*uh(i)
      enddo
      enddo
      ENDDO

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        i=ni
!$omp parallel do default(shared)   &
!$omp private(j,k)
        DO k=1,nk
        do j=1,nj
          if(rru(ni+1,j,k).le.0.0)then
            advx(i,j,k) = advx(i,j,k)-v3d(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        ENDDO
      ENDIF


!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advv=time_advv+mytime()
 
      end subroutine advvaxi


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advwaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,w3d,c1,c2,hadv)
      use input
      use constants
      implicit none

      logical, intent(in) :: doweno
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      integer, intent(in) :: hadv

      ! horizontal advection for w, axisymmetric solver:
 
      integer :: i,j,k
      real :: ubar
      double precision :: weps

!----------------------------------------------------------------

      weps = 100.0*epsilon

!----------------------------------------------------------------
! Advection in r-direction

  wenocheckwaxi:  &
  if(doweno)then

    if( weno_order.eq.3 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno3(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),weps)
        else
          dum(i,j,k) = ubar*weno3(w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.5 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno5(w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),weps)
        else
          dum(i,j,k) = ubar*weno5(w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.7 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno7(w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),weps)
        else
          dum(i,j,k) = ubar*weno7(w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    elseif( weno_order.eq.9 )then

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno9(w3d(i-5,j,k),w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),w3d(i+3,j,k),weps)
        else
          dum(i,j,k) = ubar*weno9(w3d(i+4,j,k),w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k),w3d(i-4,j,k),weps)
        endif
      enddo
      enddo
      ENDDO

    endif


    ! Near boundaries ...

    if( weno_order.eq.9 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno7(w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),weps)
        else
          dum(i,j,k) = ubar*weno9(w3d(i+4,j,k),w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k),w3d(i-4,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*weno5(w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),weps)
        else
          dum(i,j,k) = ubar*weno7(w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k),weps)
        endif
      enddo
      ENDDO
    endif
    if( weno_order.eq.9 .or. weno_order.eq.7 .or. weno_order.eq.5 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*upstrpd(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),weps)
        else
          dum(i,j,k) = ubar*weno5(w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),weps)
        endif
      enddo
      ENDDO
    endif

      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*w3d(i-1,j,k)
        else
          dum(i,j,k) = ubar*upstrpd(w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),weps)
        endif
        dum(1,j,k) = 0.0
      enddo
      ENDDO

      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*upstrpd(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),weps)
        else
          dum(i,j,k) = ubar*w3d(i,j,k)
        endif
        dum(ni+1,j,k) = 0.0
      enddo
      ENDDO
      ENDIF


  else  wenocheckwaxi


    if(hadv.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      dum(1,j,k) = 0.0
      !dir$ vector always
      do i=2,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.3)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx3(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k))
        else
          dum(i,j,k) = ubar*flx3(w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.4)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=3,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx4(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.5)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx5(w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*flx5(w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.6)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=4,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx6(w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.7)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx7(w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k))
        else
          dum(i,j,k) = ubar*flx7(w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.8)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=5,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx8(w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),w3d(i+3,j,k))
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.9)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx9(w3d(i-5,j,k),w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),w3d(i+3,j,k))
        else
          dum(i,j,k) = ubar*flx9(w3d(i+4,j,k),w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k),w3d(i-4,j,k))
        endif
      enddo
      enddo
      ENDDO
    elseif(hadv.eq.10)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=6,ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx10(w3d(i-5,j,k),w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),w3d(i+3,j,k),w3d(i+4,j,k))
      enddo
      enddo
      ENDDO
    else
      stop 34791
    endif


    ! Near boundaries ...

    if( hadv.eq.9 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx7(w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k))
        else
          dum(i,j,k) = ubar*flx9(w3d(i+4,j,k),w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k),w3d(i-4,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx5(w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k))
        else
          dum(i,j,k) = ubar*flx7(w3d(i+3,j,k),w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k),w3d(i-3,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 .or. hadv.eq.5 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*flx3(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k))
        else
          dum(i,j,k) = ubar*flx5(w3d(i+2,j,k),w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k),w3d(i-2,j,k))
        endif
      enddo
      ENDDO
    endif
    if( hadv.eq.9 .or. hadv.eq.7 .or. hadv.eq.5 .or. hadv.eq.3 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k) = ubar*w3d(i-1,j,k)
        else
          dum(i,j,k) = ubar*flx3(w3d(i+1,j,k),w3d(i  ,j,k),w3d(i-1,j,k))
        endif
        dum(1,j,k) = 0.0
      enddo
      ENDDO
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
        DO k=2,nk
        do j=1,nj
          ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dum(i,j,k) = ubar*flx3(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k))
          else
            dum(i,j,k) = ubar*w3d(i,j,k)
          endif
          dum(ni+1,j,k) = 0.0
        enddo
        ENDDO
      ENDIF
    endif



    if( hadv.eq.10 )then
      i = 5
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx8(w3d(i-4,j,k),w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k),w3d(i+3,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 )then
      i = 4
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx6(w3d(i-3,j,k),w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k),w3d(i+2,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 )then
      i = 3
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*flx4(w3d(i-2,j,k),w3d(i-1,j,k),w3d(i  ,j,k),w3d(i+1,j,k))
      enddo
      ENDDO
    endif
    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 .or. hadv.eq.4 )then
      i = 2
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
      DO k=2,nk
      do j=1,nj
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k) = ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        dum(1,j,k) = 0.0
      enddo
      ENDDO
    endif

    if( hadv.eq.10 .or. hadv.eq.8 .or. hadv.eq.6 .or. hadv.eq.4 .or. hadv.eq.2 )then
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        i = ni
!$omp parallel do default(shared)   &
!$omp private(j,k,ubar)
        DO k=2,nk
        do j=1,nj
          ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
          dum(i,j,k) = ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
          dum(ni+1,j,k) = 0.0
        enddo
        ENDDO
      ENDIF
    endif


  endif  wenocheckwaxi


    !-------------------------------------

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
      DO k=2,nk
      do j=1,nj
        i=ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.le.0.0)then
          i=ni
          dum(i+1,j,k)=dum(i,j,k)*arh1(i)/arh2(i)
        endif
      enddo
      ENDDO
      ENDIF


!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=2,nk
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        advx(i,j,k) = -(arh2(i)*dum(i+1,j,k)-arh1(i)*dum(i,j,k))*rdx*uh(i)
      enddo
      enddo
      ENDDO


      IF(ebc.eq.2 .and. ibe.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
        DO k=2,nk
        do j=1,nj
          i=ni+1
          ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
          if(ubar.le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*(                                  &
                    c1(1,1,k)*(arh2(i)*rru(i+1,j,k-1)-arh1(i)*rru(i,j,k-1))         &
                   +c2(1,1,k)*(arh2(i)*rru(i+1,j,k  )-arh1(i)*rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
        ENDDO
      ENDIF

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advw=time_advw+mytime()
 
      end subroutine advwaxi


!-----------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------------------------------
!
!                Advection for lower surface (imove)
!
!-----------------------------------------------------------------------------------------------


      subroutine movesfc(rmax,dt,weps,uh,vh,sfc,s,dum1,dum2,   &
                         reqs,west,newwest,east,neweast,       &
                         south,newsouth,north,newnorth)
      use input
#ifdef MPI
      use comm_module
#endif
      implicit none

      real, intent(in) :: rmax,dt
      double precision, intent(in) :: weps
      real, intent(in), dimension(ib:ie) :: uh
      real, intent(in), dimension(jb:je) :: vh
      real, intent(inout), dimension(ib:ie,jb:je) :: sfc,s,dum1,dum2
      integer, intent(inout), dimension(8) :: reqs
      real, intent(inout), dimension(cmp,jmp) :: west,newwest,east,neweast
      real, intent(inout), dimension(imp,cmp) :: south,newsouth,north,newnorth

      integer :: i,j,nrk
      real :: tem

!------------------------------------------------------------

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=1,nj
    do i=1,ni
      s(i,j)=sfc(i,j)
    enddo
    enddo

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=0,nj+1
    do i=0,ni+1
      dum1(i,j)=0.0
      dum2(i,j)=0.0
    enddo
    enddo

    rkloop:  &
    DO nrk=1,3

!-------------------------
!  set boundary conditions

#ifdef MPI
      if(timestats.ge.1) time_swath=time_swath+mytime()
      call comm_2d_start(s,west,newwest,east,neweast,   &
                           south,newsouth,north,newnorth,reqs)
#else
      if(wbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ngxy
          s(1-i,j)=s(ni+1-i,j)
        enddo
        enddo
      endif
      if(ebc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ngxy
          s(ni+i,j)=s(i,j)
        enddo
        enddo
      endif
      if(sbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=1,ni
          s(i,1-j)=s(i,nj+1-j)
        enddo
        enddo
      endif
      if(nbc.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=1,ni
          s(i,nj+j)=s(i,j)
        enddo
        enddo
      endif
#endif

      if(ibw.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ngxy
          s(1-i,j)=s(1,j)
        enddo
        enddo
      endif

      if(ibe.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ngxy
          s(ni+i,j)=s(ni,j)
        enddo
        enddo
      endif

      if(ibs.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=1,ni
          s(i,1-j)=s(i,1)
        enddo
        enddo
      endif

      if(ibn.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,ngxy
        do i=1,ni
          s(i,nj+j)=s(i,nj)
        enddo
        enddo
      endif

!-------------------------

#ifdef MPI
      if(timestats.ge.1) time_swath=time_swath+mytime()
      call comm_2dew_end(s,west,newwest,east,neweast,reqs)
#endif

!-------------------------

    if(abs(umove).gt.0.01)then

    if( hadvordrs.le.4 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(umove.ge.0.0)then
          dum1(i,j) = weno3(s(i+1,j),s(i  ,j),s(i-1,j),weps)
        else
          dum1(i,j) = weno3(s(i-2,j),s(i-1,j),s(i  ,j),weps)
        endif
      enddo
      enddo
    elseif( hadvordrs.eq.5 .or. hadvordrs.eq.6 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(umove.ge.0.0)then
          dum1(i,j) = weno5(s(i+2,j),s(i+1,j),s(i  ,j),s(i-1,j),s(i-2,j),weps)
        else
          dum1(i,j) = weno5(s(i-3,j),s(i-2,j),s(i-1,j),s(i  ,j),s(i+1,j),weps)
        endif
      enddo
      enddo
    elseif( hadvordrs.eq.7 .or. hadvordrs.eq.8 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(umove.ge.0.0)then
          dum1(i,j) = weno7(s(i+3,j),s(i+2,j),s(i+1,j),s(i  ,j),s(i-1,j),s(i-2,j),s(i-3,j),weps)
        else
          dum1(i,j) = weno7(s(i-4,j),s(i-3,j),s(i-2,j),s(i-1,j),s(i  ,j),s(i+1,j),s(i+2,j),weps)
        endif
      enddo
      enddo
    elseif( hadvordrs.eq.9 .or. hadvordrs.eq.10 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni+1
        if(umove.ge.0.0)then
          dum1(i,j) = weno9(s(i+4,j),s(i+3,j),s(i+2,j),s(i+1,j),s(i  ,j),s(i-1,j),s(i-2,j),s(i-3,j),s(i-4,j),weps)
        else
          dum1(i,j) = weno9(s(i-5,j),s(i-4,j),s(i-3,j),s(i-2,j),s(i-1,j),s(i  ,j),s(i+1,j),s(i+2,j),s(i+3,j),weps)
        endif
      enddo
      enddo
    endif

    endif

!-------------------------

#ifdef MPI
      if(timestats.ge.1) time_swath=time_swath+mytime()
      call comm_2dns_end(s,south,newsouth,north,newnorth,reqs)
#endif

!-------------------------

    if(abs(vmove).gt.0.01)then

    if( hadvordrs.le.4 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(vmove.ge.0.0)then
          dum2(i,j) = weno3(s(i,j+1),s(i,j  ),s(i,j-1),weps)
        else
          dum2(i,j) = weno3(s(i,j-2),s(i,j-1),s(i,j  ),weps)
        endif
      enddo
      enddo
    elseif( hadvordrs.eq.5 .or. hadvordrs.eq.6 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(vmove.ge.0.0)then
          dum2(i,j) = weno5(s(i,j+2),s(i,j+1),s(i,j  ),s(i,j-1),s(i,j-2),weps)
        else
          dum2(i,j) = weno5(s(i,j-3),s(i,j-2),s(i,j-1),s(i,j  ),s(i,j+1),weps)
        endif
      enddo
      enddo
    elseif( hadvordrs.eq.7 .or. hadvordrs.eq.8 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(vmove.ge.0.0)then
          dum2(i,j) = weno7(s(i,j+3),s(i,j+2),s(i,j+1),s(i,j  ),s(i,j-1),s(i,j-2),s(i,j-3),weps)
        else
          dum2(i,j) = weno7(s(i,j-4),s(i,j-3),s(i,j-2),s(i,j-1),s(i,j  ),s(i,j+1),s(i,j+2),weps)
        endif
      enddo
      enddo
    elseif( hadvordrs.eq.9 .or. hadvordrs.eq.10 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni
        if(vmove.ge.0.0)then
          dum2(i,j) = weno9(s(i,j+4),s(i,j+3),s(i,j+2),s(i,j+1),s(i,j  ),s(i,j-1),s(i,j-2),s(i,j-3),s(i,j-4),weps)
        else
          dum2(i,j) = weno9(s(i,j-5),s(i,j-4),s(i,j-3),s(i,j-2),s(i,j-1),s(i,j  ),s(i,j+1),s(i,j+2),s(i,j+3),weps)
        endif
      enddo
      enddo
    endif

    endif

!-------------------------

      tem=dt/(4-nrk)

!$omp parallel do default(shared)  &
!$omp private(i,j)
      do j=1,nj
      !dir$ vector always
      do i=1,ni
        s(i,j)=max( rmax , sfc(i,j)+tem*(                          &
                          umove*(dum1(i+1,j)-dum1(i,j))*rdx*uh(i)  &
                         +vmove*(dum2(i,j+1)-dum2(i,j))*rdy*vh(j) ) )
      enddo
      enddo

!-------------------------

    ENDDO  rkloop

!$omp parallel do default(shared)  &
!$omp private(i,j)
    do j=1,nj
    do i=1,ni
      sfc(i,j)=s(i,j)
    enddo
    enddo

!----------------------------------------------------------------

      if(timestats.ge.1) time_swath=time_swath+mytime()

      end subroutine movesfc


!-----------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------------------------------
!
!                Advection from specified w profile (eg, subsidence):
!
!-----------------------------------------------------------------------------------------------


    subroutine wsub(ix  ,jy  ,kz  ,a  ,wprof,c1,c2,mh,rr0,rf0,weps,dumz,subs)
    use input
    implicit none

    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    real, intent(in), dimension(kb:ke) :: wprof
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,mh,rr0,rf0
    double precision, intent(in) :: weps
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz,subs

    integer :: i,j,k
    real :: div

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,div)
    DO j=1,jy
      !-----
      do k=2,(nk-1)
        if( wprof(k).le.0.0 )then
          !dir$ vector always
          do i=1,ix
            dumz(i,j,k) = rf0(1,1,k)*wprof(k)*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
          enddo
        else
          if( k.eq.2 )then
            !dir$ vector always
            do i=1,ix
              dumz(i,j,k) = rf0(1,1,k)*wprof(k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
            enddo
          else
            !dir$ vector always
            do i=1,ix
              dumz(i,j,k) = rf0(1,1,k)*wprof(k)*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
            enddo
          endif
        endif
      enddo
      !-----
      !dir$ vector always
      do i=1,ix
        dumz(i,j,nk) = rf0(1,1,nk)*wprof(nk)*(c1(i,j,nk)*a(i,j,nk-1)+c2(i,j,nk)*a(i,j,nk))
        dumz(i,j,1) = 0.0
        dumz(i,j,nk+1) = 0.0
      enddo
      !-----
      do k=1,nk
        div = (rf0(1,1,k+1)*wprof(k+1)-rf0(1,1,k)*wprof(k))*rdz*mh(1,1,k)
        !dir$ vector always
        do i=1,ix
          subs(i,j,k) = ( -(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(1,1,k)   &
                          +a(i,j,k)*div )*rr0(1,1,k)
        enddo
      enddo
      !-----
    ENDDO

    end subroutine wsub


!-----------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------------------------------


    subroutine zsgrad(zs,gz,dzdx,dzdy,uh,vh)
    use input
    implicit none

    real, intent(in), dimension(ib:ie,jb:je) :: zs
    real, intent(inout), dimension(itb:ite,jtb:jte) :: gz,dzdx,dzdy
    real, intent(in), dimension(ib:ie) :: uh
    real, intent(in), dimension(jb:je) :: vh

    integer :: i,j,n

    n = max( hadvordrs , hadvordrv )

    do j=1,nj
      if( n.le.2 )then
        do i=1,ni+1
          gz(i,j) = 0.5*(zs(i-1,j)+zs(i,j))
        enddo
      elseif( n.eq.3 .or. n.eq.4 )then
        do i=1,ni+1
          gz(i,j) = flx4(zs(i-2,j),zs(i-1,j),zs(i  ,j),zs(i+1,j))
        enddo
      elseif( n.eq.5 .or. n.eq.6 )then
        do i=1,ni+1
          gz(i,j) = flx6(zs(i-3,j),zs(i-2,j),zs(i-1,j),zs(i  ,j),zs(i+1,j),zs(i+2,j))
        enddo
      elseif( n.eq.7 .or. n.eq.8 )then
        do i=1,ni+1
          gz(i,j) = flx8(zs(i-4,j),zs(i-3,j),zs(i-2,j),zs(i-1,j),zs(i  ,j),zs(i+1,j),zs(i+2,j),zs(i+3,j))
        enddo
      elseif( n.ge.9 )then
        do i=1,ni+1
          gz(i,j) = flx10(zs(i-5,j),zs(i-4,j),zs(i-3,j),zs(i-2,j),zs(i-1,j),zs(i  ,j),zs(i+1,j),zs(i+2,j),zs(i+3,j),zs(i+4,j))
        enddo
      endif
      do i=1,ni
        dzdx(i,j) = (gz(i+1,j)-gz(i,j))*rdx*uh(i)
      enddo
    enddo

    do i=1,ni
      if( n.le.2 )then
        do j=1,nj+1
          gz(i,j) = 0.5*(zs(i,j-1)+zs(i,j))
        enddo
      elseif( n.eq.3 .or. n.eq.4 )then
        do j=1,nj+1
          gz(i,j) = flx4(zs(i,j-2),zs(i,j-1),zs(i,j  ),zs(i,j+1))
        enddo
      elseif( n.eq.5 .or. n.eq.6 )then
        do j=1,nj+1
          gz(i,j) = flx6(zs(i,j-3),zs(i,j-2),zs(i,j-1),zs(i,j  ),zs(i,j+1),zs(i,j+2))
        enddo
      elseif( n.eq.7 .or. n.eq.8 )then
        do j=1,nj+1
          gz(i,j) = flx8(zs(i,j-4),zs(i,j-3),zs(i,j-2),zs(i,j-1),zs(i,j  ),zs(i,j+1),zs(i,j+2),zs(i,j+3))
        enddo
      elseif( n.ge.9 )then
        do j=1,nj+1
          gz(i,j) = flx10(zs(i,j-5),zs(i,j-4),zs(i,j-3),zs(i,j-2),zs(i,j-1),zs(i,j  ),zs(i,j+1),zs(i,j+2),zs(i,j+3),zs(i,j+4))
        enddo
      endif
      do j=1,nj
        dzdy(i,j) = (gz(i,j+1)-gz(i,j))*rdy*vh(j)
      enddo
    enddo

    end subroutine zsgrad


!-----------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------------------------------


  END MODULE adv_routines
