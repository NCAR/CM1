  module ib_module
  implicit none

  public

    integer :: ibib,ieib,jbib,jeib,kbib,keib
    integer :: kmaxib


  CONTAINS

!-----------------------------------------------------------------------

    subroutine ib_setup

      !-----------------------------------------
      ! USE IMMERSED BOUNDARY TECHNIQUE?

!   since cm1r21.0 ... this is now set in the namelist (param20 section)
!!!        do_ib  =  .false.    

      !-----------------------------------------

    end subroutine ib_setup

!-----------------------------------------------------------------------

    subroutine init_immersed_boundaries(                                &
                       xh,yh,xf,yf,sigma,sigmaf,zs,zh,bndy,kbdy,out3d,  &
                       west,newwest,east,neweast,                       &
                       south,newsouth,north,newnorth,reqs_p)

    use input
    use constants
    use bc_module
    use comm_module
#ifdef MPI
    use mpi
#endif

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      logical, intent(inout), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(inout), dimension(ibib:ieib,jbib:jeib) :: kbdy
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(cmp,jmp) :: west,newwest,east,neweast
      real, intent(inout), dimension(imp,cmp) :: south,newsouth,north,newnorth
      real, intent(inout), optional, dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

      integer :: i,j,k,n
      real :: aa,hh,xc,xloc,yloc,foox,fooy,shift
      real, dimension(:,:), allocatable :: zsfoo


      ! Define immersed grid cells when do_ib = .true.


      ! Set case here:

      ib_init  =  4



      allocate( zsfoo(ib:ie,jb:je) )
      zsfoo = 0.0

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( ib_init.eq.1 )THEN

        !  standard nh mountain wave case:

        hh =      400.0              ! max. height (m)
        aa =     1000.0              ! half width (m)
        xc =        0.0 + 0.5*dx     ! x-location (m)

        do j=jb,je
        do i=ib,ie
          zsfoo(i,j) = hh/( 1.0+( (xh(i)-xc)/aa )**2 )
        enddo
        enddo

      ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( ib_init.eq.2 )THEN

        !  2D block for tests of advection scheme:

        hh =      400.0              ! max. height (m)

        if( nx.eq.1 )then
          do j=jb,je
          do i=ib,ie
            if( abs(yh(j)-centery).le.400.0 ) zsfoo(i,j) = hh
          enddo
          enddo
        else
          do j=jb,je
          do i=ib,ie
            if( abs(xh(i)-centerx).le.400.0 ) zsfoo(i,j) = hh
          enddo
          enddo
        endif

      ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

      IF( ib_init.eq.3 )THEN

        ! multiple cubes

      do n=1,1
        if( n.eq.1 ) yloc = 0.5*maxy
        if( n.eq.2 ) yloc = 0.5*maxy + 3.0*100.0
        if( n.eq.3 ) yloc = 0.5*maxy - 3.0*100.0
        print *,'  n,yloc = ',n,yloc
        do j=jb,je
        do i=ib,ie
          ! place obstacle roughly 1/3 across the domain:
          if( abs(xh(i)-(minx+(maxx-minx)/3.0+360.0)).le.60.0 .and. abs(yh(j)-centery).le.60.0 ) zsfoo(i,j) = 120.0
        enddo
        enddo
      enddo

      ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    IF( ib_init.eq.4 )THEN

      ! Martinuzzi and Tropea (1993, JFE) wind tunnel channel case:

    IF( ny.eq.1 )THEN
      do j=jb,je
      do i=ib,ie
        if( abs(xh(i)-(minx+0.50*(maxx-minx)) ).le.(0.25*maxz) )then
          zsfoo(i,j) = 0.5*maxz
        endif
      enddo
      enddo
    ELSEIF( nx.eq.1 )THEN
      stop 12321
    ELSE
      do j=jb,je
      do i=ib,ie
        if( abs(xh(i)-(minx+0.50*(maxx-minx)) ).le.(0.25*maxz) .and.  &
            abs(yh(j)-centery).le.(0.25*maxz) )then
          zsfoo(i,j) = 0.5*maxz
        endif
      enddo
      enddo
    ENDIF

    ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    IF( ib_init.eq.5 )THEN

      ! cube in center of domain:

      do j=jb,je
      do i=ib,ie
        if( abs(xh(i)-centerx).le.250.0 .and. abs(yh(j)-centery).le.250.0 ) zsfoo(i,j) = 500.0
      enddo
      enddo

    ENDIF

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c

    !---  DO NOT CHANGE ANYTHING BELOW HERE  ---!

        call bc2d(zsfoo)
#ifdef MPI
      nf=0
      nu=0
      nv=0
      nw=0
      call comm_2d_start(zsfoo,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs_p)
      call comm_2dew_end(zsfoo,west,newwest,east,neweast,reqs_p)
      call comm_2dns_end(zsfoo,south,newsouth,north,newnorth,reqs_p)
      call bcs2_2d(zsfoo)
      call bc2d(zsfoo)
      call getcorner3_2d(zsfoo)
#endif

      do k=1,nk
      do j=jb,je
      do i=ib,ie
        if( sigma(k).le.zsfoo(i,j) ) bndy(i,j,k) = .true.
      enddo
      enddo
      enddo

      ! last step:  get kbdy

      kmaxib = 0

      do k=1,nk
      do j=jb,je
      do i=ib,ie
        if( bndy(i,j,k) )then
          kbdy(i,j) = k+1
        endif
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        kmaxib = max( kmaxib , kbdy(i,j) )
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(mpi_in_place,kmaxib,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '    kmaxib    = ',kmaxib
      if(dowr) write(outfile,*)

      IF( kmaxib .gt. (nk-3) )THEN
        if(myid.eq.0)then
        print *
        print *,'  nk     = ',nk
        print *,'  kmaxib = ',kmaxib
        print *
        print *,'  kmaxib must be < (nk-3) '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

      deallocate( zsfoo )

    end subroutine init_immersed_boundaries

!-----------------------------------------------------------------------

    subroutine zero_out_uv(bndy,kbdy,u3d,v3d)

    use input
    use constants

      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d

      integer :: i,j,k

          ! set u to zero on east/west faces of immersed gridpoints:
          do k=1,kmaxib
          do j=0,nj+1
          do i=0,ni+1
            if( bndy(i,j,k) )then
              if( .not. bndy(i-1,j,k) )  u3d(i  ,j,k) = 0.0
              if( .not. bndy(i+1,j,k) )  u3d(i+1,j,k) = 0.0
              if( .not. bndy(i,j-1,k) )  v3d(i,j  ,k) = 0.0
              if( .not. bndy(i,j+1,k) )  v3d(i,j+1,k) = 0.0
            endif
          enddo
          enddo
          enddo

          ! set u to zero on east/west faces of immersed gridpoints:
          do k=1,kmaxib
          do j=1,nj
          do i=1,ni+1
            if( bndy(i-1,j,k).and.bndy(i,j,k) ) u3d(i,j,k) = 0.0
          enddo
          enddo
          enddo

          ! set v to zero on south/north faces of immersed gridpoints:
          do k=1,kmaxib
          do j=1,nj+1
          do i=1,ni
            if( bndy(i,j-1,k).and.bndy(i,j,k) ) v3d(i,j,k) = 0.0
          enddo
          enddo
          enddo

    end subroutine zero_out_uv

!-----------------------------------------------------------------------

    subroutine zero_out_w(bndy,kbdy,w3d)

    use input
    use constants

      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d

      integer :: i,j,k

        ! set gridpoints in/on immersed gridpoints to zero:
        do j=1,nj
        do i=1,ni
          do k = 1,kbdy(i,j)
            w3d(i,j,k) = 0.0
          enddo
        enddo
        enddo

    end subroutine zero_out_w

!-----------------------------------------------------------------------

    subroutine drag_obstacles(xh,yh,zh,zf,rho,rf,dum1,dum2,dum3,dum4,dum5,dum6,t11,t12,t13,t22,t23,t33,ua,va,wa,kbdy)

    use input
    use constants

    real, intent(in), dimension(ib:ie) :: xh
    real, intent(in), dimension(jb:je) :: yh
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
    integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy

    integer :: i,j,k,kk
    real :: tem,tem2,rhox
    logical :: doit

!-----------------------------

!  drag on top & sides of obstacles


!  NOTE: if you are looking to set top_cd or side_cd, they are now specified 
!        in namelist.input (param20 section)


  doit = .true.
!!!  IF( doit )THEN
!!!    stop 4444

    doib:  &
    IF( do_ib )THEN

!!!      if( myid.eq.0 ) print *,'  top_cd,side_cd = ',top_cd,side_cd

      do k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k) = 0.5*(ua(i,j,k)+ua(i+1,j,k))
        dum2(i,j,k) = 0.5*(va(i,j,k)+va(i,j+1,k))
        dum4(i,j,k) = 0.5*(wa(i,j,k)+wa(i,j,k+1))
        dum3(i,j,k) = sqrt( dum1(i,j,k)**2 + dum2(i,j,k)**2 )
        dum5(i,j,k) = sqrt( dum4(i,j,k)**2 + dum1(i,j,k)**2 )  ! for south/north bndys
        dum6(i,j,k) = sqrt( dum4(i,j,k)**2 + dum2(i,j,k)**2 )  ! for west/east bndys
      enddo
      enddo
      enddo

!-----------------------------
!  sides of obstacles:

      tem = 0.25*side_cd

      do j=0,nj+1
      do i=0,ni+1
        kbdy_gt_1:  &
        if( kbdy(i,j).gt.1 )then
          do k=1,(kbdy(i,j)-1)
            t11(i,j,k) = 0.0
            t22(i,j,k) = 0.0
            t33(i,j,k) = 0.0
          enddo
          !-----
          if( kbdy(i-1,j).le.1 )then
            ! west side of obstacle:
            do k=1,(kbdy(i,j)-1)
              t12(i,j  ,k) = -tem * va(i-1,j  ,k) * (dum6(i-1,j,k)+dum6(i-1,j-1,k)) * (rho(i-1,j,k)+rho(i-1,j-1,k))
              t12(i,j+1,k) = -tem * va(i-1,j+1,k) * (dum6(i-1,j,k)+dum6(i-1,j+1,k)) * (rho(i-1,j,k)+rho(i-1,j+1,k))
              t13(i,j,k+1) = -tem * wa(i-1,j,k+1) * (dum6(i-1,j,k)+dum6(i-1,j,k+1)) * (rho(i-1,j,k)+rho(i-1,j,k+1))
            enddo
          endif
          if( kbdy(i+1,j).le.1 )then
            ! east side of obstacle:
            do k=1,(kbdy(i,j)-1)
              t12(i+1,j  ,k) =  tem * va(i+1,j  ,k) * (dum6(i+1,j,k)+dum6(i+1,j-1,k)) * (rho(i+1,j,k)+rho(i+1,j-1,k))
              t12(i+1,j+1,k) =  tem * va(i+1,j+1,k) * (dum6(i+1,j,k)+dum6(i+1,j+1,k)) * (rho(i+1,j,k)+rho(i+1,j+1,k))
              t13(i+1,j,k+1) =  tem * wa(i+1,j,k+1) * (dum6(i+1,j,k)+dum6(i+1,j,k+1)) * (rho(i+1,j,k)+rho(i+1,j,k+1))
            enddo
          endif
          !-----
          if( kbdy(i,j-1).le.1 )then
            ! south side of obstacle:
            do k=1,(kbdy(i,j)-1)
              t12(i  ,j,k) = -tem * ua(i  ,j-1,k) * (dum5(i,j-1,k)+dum5(i-1,j-1,k)) * (rho(i,j-1,k)+rho(i-1,j-1,k))
              t12(i+1,j,k) = -tem * ua(i+1,j-1,k) * (dum5(i,j-1,k)+dum5(i+1,j-1,k)) * (rho(i,j-1,k)+rho(i+1,j-1,k))
              t23(i,j,k+1) = -tem * wa(i,j-1,k+1) * (dum5(i,j-1,k)+dum5(i,j-1,k+1)) * (rho(i,j-1,k)+rho(i,j-1,k+1))
            enddo
          endif
          if( kbdy(i,j+1).le.1 )then
            ! north side of obstacle:
            do k=1,(kbdy(i,j)-1)
              t12(i  ,j+1,k) =  tem * ua(i  ,j+1,k) * (dum5(i,j+1,k)+dum5(i-1,j+1,k)) * (rho(i,j+1,k)+rho(i-1,j+1,k))
              t12(i+1,j+1,k) =  tem * ua(i+1,j+1,k) * (dum5(i,j+1,k)+dum5(i+1,j+1,k)) * (rho(i,j+1,k)+rho(i+1,j+1,k))
              t23(i,j+1,k+1) =  tem * wa(i,j+1,k+1) * (dum5(i,j+1,k)+dum5(i,j+1,k+1)) * (rho(i,j+1,k)+rho(i,j+1,k+1))
            enddo
          endif
          !-----
          ! check for corners:
          ! i am not sure what to do, so just set t12 to zero (for now)
          if( kbdy(i-1,j).le.1 .and. kbdy(i,j+1).le.1 .and. kbdy(i-1,j+1).le.1 )then
            ! nw corner:
            do k=1,(kbdy(i,j)-1)
              t12(i,j+1,k) = 0.0
            enddo
          endif
          if( kbdy(i+1,j).le.1 .and. kbdy(i,j+1).le.1 .and. kbdy(i+1,j+1).le.1 )then
            ! ne corner:
            do k=1,(kbdy(i,j)-1)
              t12(i+1,j+1,k) = 0.0
            enddo
          endif
          if( kbdy(i-1,j).le.1 .and. kbdy(i,j-1).le.1 .and. kbdy(i-1,j-1).le.1 )then
            ! sw corner:
            do k=1,(kbdy(i,j)-1)
              t12(i,j,k) = 0.0
            enddo
          endif
          if( kbdy(i+1,j).le.1 .and. kbdy(i,j-1).le.1 .and. kbdy(i+1,j-1).le.1 )then
            ! se corner:
            do k=1,(kbdy(i,j)-1)
              t12(i+1,j,k) = 0.0
            enddo
          endif
          !-----
          if( kbdy(i-1,j).gt.1 .and. kbdy(i,j-1).gt.1 .and. kbdy(i-1,j-1).gt.1 )then
            ! embedded on all sides:
            do k=1,(kbdy(i,j)-1)
              t12(i,j,k) = 0.0
            enddo
          endif
          !-----
        endif  kbdy_gt_1
      enddo
      enddo

!-----------------------------
!  tops of obstacles:

      ! t13:
      do j=1,nj
      do i=1,ni+1
        if( kbdy(i-1,j).gt.1 .or. kbdy(i,j).gt.1 )then
          k = max( kbdy(i-1,j) , kbdy(i,j) )
          rhox = 0.5*(   &
      (rho(i-1,j,k)+(rho(i-1,j,k+1)-rho(i-1,j,k))*(zf(i-1,j,k)-zh(i-1,j,k))/(zh(i-1,j,k+1)-zh(i-1,j,k))) &
     +(rho(i  ,j,k)+(rho(i  ,j,k+1)-rho(i  ,j,k))*(zf(i  ,j,k)-zh(i  ,j,k))/(zh(i  ,j,k+1)-zh(i  ,j,k))) )
          t13(i,j,k) = top_cd * rhox * ua(i,j,k) * 0.5*(dum3(i-1,j,k)+dum3(i,j,k))
        endif
        if( kbdy(i-1,j).gt.1 .and. kbdy(i,j).gt.1 )then
          kk = max( kbdy(i-1,j) , kbdy(i,j) ) - 1
          do k=1,kk
            t13(i,j,k) = 0.0
          enddo
        endif
      enddo
      enddo

      ! t23:
      do j=1,nj+1
      do i=1,ni
        if( kbdy(i,j-1).gt.1 .or. kbdy(i,j).gt.1 )then
          k = max( kbdy(i,j-1),kbdy(i,j) )
          rhox = 0.5*(   &
      (rho(i,j-1,k)+(rho(i,j-1,k+1)-rho(i,j-1,k))*(zf(i,j-1,k)-zh(i,j-1,k))/(zh(i,j-1,k+1)-zh(i,j-1,k))) &
     +(rho(i,j  ,k)+(rho(i,j  ,k+1)-rho(i,j  ,k))*(zf(i,j  ,k)-zh(i,j  ,k))/(zh(i,j  ,k+1)-zh(i,j  ,k))) )
          t23(i,j,k) = top_cd * rhox * va(i,j,k) * 0.5*(dum3(i,j-1,k)+dum3(i,j,k))
        endif
        if( kbdy(i,j-1).gt.1 .and. kbdy(i,j).gt.1 )then
          kk = max( kbdy(i,j-1),kbdy(i,j) ) - 1
          do k=1,kk
            t23(i,j,k) = 0.0
          enddo
        endif
      enddo
      enddo

    ENDIF  doib

!!!  ENDIF

    end subroutine drag_obstacles

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

      !dir$ attributes forceinline :: weno3
      real function weno3(s1,s2,s3,weps)
      implicit none

      real, intent(in) :: s1,s2,s3
      double precision, intent(in) :: weps

      double precision :: b1,b2
      double precision :: w1,w2

      integer, parameter :: siform = 2

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

      integer, parameter :: siform = 2

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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ib_flx_init(bndy,hflxw,hflxe,hflxs,hflxn,dum1,dum2,     &
                           uw31,uw32,ue31,ue32,us31,us32,un31,un32,    &
                           vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,reqs_u,reqs_v)

    use input
    use constants
    use bc_module
    use comm_module
#ifdef MPI
    use mpi
#endif

    logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
    integer, intent(inout), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
    real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
    real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
    real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
    real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
    integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v

    integer :: i,j,k
    logical :: doit

    DO k=1,kmaxib
      do j=jb,je
      do i=ib,ie-1
        if( bndy(i,j,k) .and. ( .not. bndy(i+1,j,k) ) )then
          !---  east face  ---!
          if( i+1 .le. ie                      ) hflxe(i+1,j,k) = 1
          if( i+2 .le. ie .and. hadvordrs.ge.3 ) hflxe(i+2,j,k) = 3
          if( i+3 .le. ie .and. hadvordrs.ge.5 ) hflxe(i+3,j,k) = 5
        endif
      enddo
      enddo
    ENDDO
    DO k=1,kmaxib
      do j=jb,je
      do i=ib+1,ie
        if( bndy(i,j,k) .and. ( .not. bndy(i-1,j,k) ) )then
          !---  west face  ---!
          if( i   .ge. ib                      ) hflxw(i  ,j,k) = 1
          if( i-1 .ge. ib .and. hadvordrs.ge.3 ) hflxw(i-1,j,k) = 3
          if( i-2 .ge. ib .and. hadvordrs.ge.5 ) hflxw(i-2,j,k) = 5
        endif
      enddo
      enddo
    ENDDO
    DO k=1,kmaxib
      do j=jb,je-1
      do i=ib,ie
        if( bndy(i,j,k) .and. ( .not. bndy(i,j+1,k) ) )then
          !---  north face  ---!
          if( j+1 .le. je                      ) hflxn(i,j+1,k) = 1
          if( j+2 .le. je .and. hadvordrs.ge.3 ) hflxn(i,j+2,k) = 3
          if( j+3 .le. je .and. hadvordrs.ge.5 ) hflxn(i,j+3,k) = 5
        endif
      enddo
      enddo
    ENDDO
    DO k=1,kmaxib
      do j=jb+1,je
      do i=ib,ie
        if( bndy(i,j,k) .and. ( .not. bndy(i,j-1,k) ) )then
          !---  south face  ---!
          if( j   .ge. jb                      ) hflxs(i,j  ,k) = 1
          if( j-1 .ge. jb .and. hadvordrs.ge.3 ) hflxs(i,j-1,k) = 3
          if( j-2 .ge. jb .and. hadvordrs.ge.5 ) hflxs(i,j-2,k) = 5
        endif
      enddo
      enddo
    ENDDO

      end subroutine ib_flx_init

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ib_side_flx(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a,hadvorder,bndy,hflxw,hflxe,hflxs,hflxn,out3d)
    use input
    implicit none

    ! reduce order of horizontal advective fluxes near sides of obstacles:

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: hadvorder
    logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
    integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
    real, intent(inout), optional, dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2
    logical :: doit

    doit = .true.
!!!    if( doit )then

  IF( stag.eq.1 )THEN

    DO k=1,(kmaxib-1)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxe(i,j,k).eq.3 )then
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k) = rru(i,j,k)*0.5*(a(i-1,j,k)+a(i,j,k))
          else
            dumx(i,j,k) = rru(i,j,k)*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
          endif
        endif
        if( hflxe(i,j,k).eq.5 )then
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k) = rru(i,j,k)*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxw(i,j,k).eq.3 )then
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k) = rru(i,j,k)*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
          else
            dumx(i,j,k) = rru(i,j,k)*0.5*(a(i-1,j,k)+a(i,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 )then
          if(rru(i,j,k).lt.0.0)then
            dumx(i,j,k) = rru(i,j,k)*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxn(i,j,k).eq.3 )then
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*0.5*(a(i,j-1,k)+a(i,j,k))
          else
            dumy(i,j,k) = rrv(i,j,k)*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
          endif
        endif
        if( hflxn(i,j,k).eq.5 )then
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxs(i,j,k).eq.3 )then
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
          else
            dumy(i,j,k) = rrv(i,j,k)*0.5*(a(i,j-1,k)+a(i,j  ,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 )then
          if(rrv(i,j,k).lt.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.2 )THEN

    if( doit )then
    DO k=1,(kmaxib-1)
      do j=1,nj+1
      !dir$ vector always
      do i=0,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.1 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*0.5*(a(i  ,j,k)+a(i+1,j,k))
          else
            dumx(i,j,k) = ubar*flx3(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k))
          endif
        endif
        if( hflxe(i,j,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
          else
            dumx(i,j,k) = ubar*0.5*(a(i  ,j,k)+a(i+1,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.lt.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k))
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.3 .or. hflxn(i-1,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
          else
            dumy(i,j,k) = vbar*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
          endif
        endif
        if( hflxn(i,j,k).eq.5 .or. hflxn(i-1,j,k).eq.5 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.3 .or. hflxs(i-1,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
          else
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 .or. hflxs(i-1,j,k).eq.5 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if( vbar.lt.0.0 )then
            dumy(i,j,k) = vbar*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO
    endif

  ELSEIF( stag.eq.3 )THEN

    if( doit )then
    DO k=1,(kmaxib-1)
      do j=0,nj+1
      !dir$ vector always
      do i=1,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.3 .or. hflxe(i,j-1,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i  ,j,k))
          else
            dumx(i,j,k) = ubar*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
          endif
        endif
        if( hflxe(i,j,k).eq.5 .or. hflxe(i,j-1,k).eq.5 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.3 .or. hflxw(i,j-1,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
          else
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i  ,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 .or. hflxw(i,j-1,k).eq.5 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.lt.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.1 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*0.5*(a(i,j  ,k)+a(i,j+1,k))
          else
            dumy(i,j,k) = vbar*flx3(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k))
          endif
        endif
        if( hflxn(i,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*flx3(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*flx3(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
          else
            dumy(i,j,k) = vbar*0.5*(a(i,j  ,k)+a(i,j+1,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if( vbar.lt.0.0 )then
            dumy(i,j,k) = vbar*flx3(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k))
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO
    endif

  ELSEIF( stag.eq.4 )THEN

  if( doit )then
    DO k=1,(kmaxib-1)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxe(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i,j,k))
          else
            dumx(i,j,k) = ubar*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
          endif
        endif
        if( hflxe(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
          else
            dumx(i,j,k) = ubar*flx5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k))
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxw(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx3(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k))
          else
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
          else
            dumx(i,j,k) = ubar*flx3(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k))
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxn(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j,k))
          else
            dumy(i,j,k) = vbar*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
          endif
        endif
        if( hflxn(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
          else
            dumy(i,j,k) = vbar*flx5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k))
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxs(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*flx3(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k))
          else
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*flx5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
          else
            dumy(i,j,k) = vbar*flx3(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k))
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO
  endif

  ENDIF

!!!    endif

    end subroutine ib_side_flx

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ib_side_weno(stag,ix,jy,kz,c1,c2,rru,rrv,dumx,dumy,a,hadvorder,bndy,hflxw,hflxe,hflxs,hflxn,weps,out3d)
    use input
    implicit none

    ! reduce order of horizontal advective fluxes near sides of obstacles:

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumx,dumy
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in) :: hadvorder
    logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
    integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
    double precision, intent(in) :: weps
    real, intent(inout), optional, dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

    integer :: i,j,k,i1,i2,j1,j2
    real :: ubar,vbar,cc1,cc2
    logical :: doit

    doit = .true.
!!!    if( doit )then

  IF( stag.eq.1 )THEN

    DO k=1,(kmaxib-1)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxe(i,j,k).eq.3 )then
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k) = rru(i,j,k)*0.5*(a(i-1,j,k)+a(i,j,k))
          else
            dumx(i,j,k) = rru(i,j,k)*upstrpd(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
          endif
        endif
        if( hflxe(i,j,k).eq.5 )then
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k) = rru(i,j,k)*upstrpd(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxw(i,j,k).eq.3 )then
          if(rru(i,j,k).ge.0.0)then
            dumx(i,j,k) = rru(i,j,k)*upstrpd(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
          else
            dumx(i,j,k) = rru(i,j,k)*0.5*(a(i-1,j,k)+a(i,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 )then
          if(rru(i,j,k).lt.0.0)then
            dumx(i,j,k) = rru(i,j,k)*upstrpd(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxn(i,j,k).eq.3 )then
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*0.5*(a(i,j-1,k)+a(i,j,k))
          else
            dumy(i,j,k) = rrv(i,j,k)*upstrpd(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
          endif
        endif
        if( hflxn(i,j,k).eq.5 )then
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*upstrpd(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxs(i,j,k).eq.3 )then
          if(rrv(i,j,k).ge.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*upstrpd(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
          else
            dumy(i,j,k) = rrv(i,j,k)*0.5*(a(i,j-1,k)+a(i,j  ,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 )then
          if(rrv(i,j,k).lt.0.0)then
            dumy(i,j,k) = rrv(i,j,k)*upstrpd(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO

  ELSEIF( stag.eq.2 )THEN

    if( doit )then
    DO k=1,(kmaxib-1)
      do j=1,nj+1
      !dir$ vector always
      do i=0,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.1 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*0.5*(a(i  ,j,k)+a(i+1,j,k))
          else
            dumx(i,j,k) = ubar*upstrpd(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),weps)
          endif
        endif
        if( hflxe(i,j,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i-1,j,k),a(i  ,j,k),a(i+1,j,k),weps)
          else
            dumx(i,j,k) = ubar*0.5*(a(i  ,j,k)+a(i+1,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 )then
          ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
          if(ubar.lt.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.3 .or. hflxn(i-1,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
          else
            dumy(i,j,k) = vbar*upstrpd(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
          endif
        endif
        if( hflxn(i,j,k).eq.5 .or. hflxn(i-1,j,k).eq.5 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*upstrpd(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.3 .or. hflxs(i-1,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*upstrpd(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
          else
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 .or. hflxs(i-1,j,k).eq.5 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
          if( vbar.lt.0.0 )then
            dumy(i,j,k) = vbar*upstrpd(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO
    endif

  ELSEIF( stag.eq.3 )THEN

    if( doit )then
    DO k=1,(kmaxib-1)
      do j=0,nj+1
      !dir$ vector always
      do i=1,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.3 .or. hflxe(i,j-1,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i  ,j,k))
          else
            dumx(i,j,k) = ubar*upstrpd(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
          endif
        endif
        if( hflxe(i,j,k).eq.5 .or. hflxe(i,j-1,k).eq.5 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.3 .or. hflxw(i,j-1,k).eq.3 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
          else
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i  ,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 .or. hflxw(i,j-1,k).eq.5 )then
          ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
          if(ubar.lt.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.1 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*0.5*(a(i,j  ,k)+a(i,j+1,k))
          else
            dumy(i,j,k) = vbar*upstrpd(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),weps)
          endif
        endif
        if( hflxn(i,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*upstrpd(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),weps)
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.3 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if( vbar.ge.0.0 )then
            dumy(i,j,k) = vbar*upstrpd(a(i,j-1,k),a(i,j  ,k),a(i,j+1,k),weps)
          else
            dumy(i,j,k) = vbar*0.5*(a(i,j  ,k)+a(i,j+1,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 )then
          vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
          if( vbar.lt.0.0 )then
            dumy(i,j,k) = vbar*upstrpd(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),weps)
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO
    endif

  ELSEIF( stag.eq.4 )THEN

  if( doit )then
    DO k=1,(kmaxib-1)
      do j=1,nj+1
      !dir$ vector always
      do i=1,ni+1
        !---------------------------------------
          !---  east face  ---!
        if( hflxe(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxe(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i,j,k))
          else
            dumx(i,j,k) = ubar*upstrpd(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
          endif
        endif
        if( hflxe(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
          else
            dumx(i,j,k) = ubar*flx5(a(i+2,j,k),a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),a(i-2,j,k))
          endif
        endif
        !---------------------------------------
          !---  west face  ---!
        if( hflxw(i,j,k).eq.1 )then
          dumx(i,j,k) = 0.0
        endif
        if( hflxw(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*upstrpd(a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),weps)
          else
            dumx(i,j,k) = ubar*0.5*(a(i-1,j,k)+a(i,j,k))
          endif
        endif
        if( hflxw(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            dumx(i,j,k) = ubar*flx5(a(i-3,j,k),a(i-2,j,k),a(i-1,j,k),a(i  ,j,k),a(i+1,j,k))
          else
            dumx(i,j,k) = ubar*upstrpd(a(i+1,j,k),a(i  ,j,k),a(i-1,j,k),weps)
          endif
        endif
        !---------------------------------------
          !---  north face  ---!
        if( hflxn(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxn(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j,k))
          else
            dumy(i,j,k) = vbar*upstrpd(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
          endif
        endif
        if( hflxn(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*upstrpd(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
          else
            dumy(i,j,k) = vbar*flx5(a(i,j+2,k),a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),a(i,j-2,k))
          endif
        endif
        !---------------------------------------
          !---  south face  ---!
        if( hflxs(i,j,k).eq.1 )then
          dumy(i,j,k) = 0.0
        endif
        if( hflxs(i,j,k).eq.3 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*upstrpd(a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),weps)
          else
            dumy(i,j,k) = vbar*0.5*(a(i,j-1,k)+a(i,j  ,k))
          endif
        endif
        if( hflxs(i,j,k).eq.5 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            dumy(i,j,k) = vbar*flx5(a(i,j-3,k),a(i,j-2,k),a(i,j-1,k),a(i,j  ,k),a(i,j+1,k))
          else
            dumy(i,j,k) = vbar*upstrpd(a(i,j+1,k),a(i,j  ,k),a(i,j-1,k),weps)
          endif
        endif
        !---------------------------------------
      enddo
      enddo
    ENDDO
  endif

  ENDIF

!!!    endif

    end subroutine ib_side_weno

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ib_lwr_flx(stag,ix,jy,kz,c1,c2,rrw,dumz,a,kbdy,vadvorder)
    use input
    implicit none

    ! reduce order of vertical advective fluxes near top of obstacles:

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
    integer, intent(in) :: vadvorder

    integer :: i,j,k,i1,i2,j1,j2,kval
    real :: wbar,cc1,cc2

  IF( stag.eq.1 )THEN

    do j=1,nj
    !dir$ vector always
    do i=1,ni
      if( kbdy(i,j).gt.1 )then
        do k=1,kbdy(i,j)
          dumz(i,j,k) = 0.0
        enddo
        k = kbdy(i,j)+1
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
        if( vadvorder.ge.5 )then
          k = kbdy(i,j)+2
          if(rrw(i,j,k).ge.0.0)then
            dumz(i,j,k) = rrw(i,j,k)*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
          else
            dumz(i,j,k) = rrw(i,j,k)*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
          endif
        endif
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

    do j=1,nj
    !dir$ vector always
    do i=i1,i2
      kval = max(kbdy(i-1,j),kbdy(i,j))
      if( kval.gt.1 )then
        do k=1,kval
          dumz(i,j,k) = 0.0
        enddo
        k = kval+1
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
        if( vadvorder.ge.5 )then
          k = kval+2
          wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
          if(wbar.ge.0.0)then
            dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
          else
            dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
          endif
        endif
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

    do j=j1,j2
    !dir$ vector always
    do i=1,ni
      kval = max(kbdy(i,j-1),kbdy(i,j))
      if( kval.gt.1 )then
        do k=1,kval
          dumz(i,j,k) = 0.0
        enddo
        k = kval+1
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
        endif
        if( vadvorder.ge.5 )then
          k = kval+2
          wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
          if(wbar.ge.0.0)then
            dumz(i,j,k) = wbar*flx3(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ))
          else
            dumz(i,j,k) = wbar*flx5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2))
          endif
        endif
      endif
    enddo
    enddo

  ELSEIF( stag.eq.4 )THEN

    do j=1,nj
    !dir$ vector always
    do i=1,ni
      if( kbdy(i,j).gt.1 )then
        do k=1,(kbdy(i,j)-1)
          dumz(i,j,k) = 0.0
        enddo
        k = kbdy(i,j)
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*flx3(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ))
        endif
        if( vadvorder.ge.5 )then
          k = kbdy(i,j)+1
          wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
          if(wbar.ge.0.0)then
            dumz(i,j,k) = wbar*flx3(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1))
          else
            dumz(i,j,k) = wbar*flx5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1))
          endif
        endif
      endif
    enddo
    enddo

  ENDIF

    end subroutine ib_lwr_flx

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ib_lwr_weno(stag,ix,jy,kz,c1,c2,rrw,dumz,a,kbdy,vadvorder,weps)
    use input
    implicit none

    ! reduce order of vertical advective fluxes near top of obstacles (for weno):

    integer, intent(in) :: stag
    integer, intent(in) :: ix,jy,kz
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dumz
    real, intent(in), dimension(1-ngxy:ix+ngxy,1-ngxy:jy+ngxy,1-ngz:kz+ngz)   :: a
    integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
    integer, intent(in) :: vadvorder
    double precision, intent(in) :: weps

    integer :: i,j,k,i1,i2,j1,j2,kval
    real :: wbar,cc1,cc2

  IF( stag.eq.1 )THEN

    do j=1,nj
    !dir$ vector always
    do i=1,ni
      if( kbdy(i,j).gt.1 )then
        do k=1,kbdy(i,j)
          dumz(i,j,k) = 0.0
        enddo
        k = kbdy(i,j)+1
        if(rrw(i,j,k).ge.0.0)then
          dumz(i,j,k) = rrw(i,j,k)*(c1(i,j,k)*a(i,j,k-1)+c2(i,j,k)*a(i,j,k))
        else
          dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
        if( vadvorder.ge.5 )then
          k = kbdy(i,j)+2
          if(rrw(i,j,k).ge.0.0)then
            dumz(i,j,k) = rrw(i,j,k)*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
          else
            dumz(i,j,k) = rrw(i,j,k)*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
          endif
        endif
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

    do j=1,nj
    !dir$ vector always
    do i=i1,i2
      kval = max(kbdy(i-1,j),kbdy(i,j))
      if( kval.gt.1 )then
        do k=1,kval
          dumz(i,j,k) = 0.0
        enddo
        k = kval+1
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
        if( vadvorder.ge.5 )then
          k = kval+2
          wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
          if( wbar.ge.0.0 )then
            dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
          else
            dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
          endif
        endif
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

    do j=j1,j2
    !dir$ vector always
    do i=1,ni
      kval = max(kbdy(i,j-1),kbdy(i,j))
      if( kval.gt.1 )then
        do k=1,kval
          dumz(i,j,k) = 0.0
        enddo
        k = kval+1
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dumz(i,j,k) = wbar*(cc1*a(i,j,k-1)+cc2*a(i,j,k))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
        endif
        if( vadvorder.ge.5 )then
          k = kval+2
          wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
          if( wbar.ge.0.0 )then
            dumz(i,j,k) = wbar*upstrpd(a(i,j,k-2),a(i,j,k-1),a(i,j,k  ),weps)
          else
            dumz(i,j,k) = wbar*weno5(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),a(i,j,k-2),weps)
          endif
        endif
      endif
    enddo
    enddo

  ELSEIF( stag.eq.4 )THEN

    do j=1,nj
    !dir$ vector always
    do i=1,ni
      if( kbdy(i,j).gt.1 )then
        do k=1,(kbdy(i,j)-1)
          dumz(i,j,k) = 0.0
        enddo
        k = kbdy(i,j)
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          dumz(i,j,k) = wbar*0.5*(a(i,j,k)+a(i,j,k+1))
        else
          dumz(i,j,k) = wbar*upstrpd(a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),weps)
        endif
        if( vadvorder.ge.5 )then
          k = kbdy(i,j)+1
          wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
          if( wbar.ge.0.0 )then
            dumz(i,j,k) = wbar*upstrpd(a(i,j,k-1),a(i,j,k  ),a(i,j,k+1),weps)
          else
            dumz(i,j,k) = wbar*weno5(a(i,j,k+3),a(i,j,k+2),a(i,j,k+1),a(i,j,k  ),a(i,j,k-1),weps)
          endif
        endif
      endif
    enddo
    enddo

  ENDIF

    end subroutine ib_lwr_weno

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  end module ib_module
