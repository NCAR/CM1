  module init_surface_module

  implicit none

  contains

      subroutine init_surface(num_soil_layers,dosfcflx,xh,ruh,xf,yh,rvh,yf,   &
           lu_index,xland,tsk,slab_zs,slab_dzs,tslb, &
           emiss,thc,albd,znt,rznt,mavail,dsxy,prs0s,prs0,   &
           tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,z0base)

      use input
      use constants
      use module_sf_oml

      implicit none

      integer, intent(in) :: num_soil_layers
      logical, intent(inout) :: dosfcflx
      real, intent(in), dimension(ib:ie) :: xh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,rznt,xland
      real, intent(inout), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: emiss,thc,albd,mavail,dsxy
      real, intent(in), dimension(ib:ie,jb:je) :: prs0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs0
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: z0base

      integer :: i,j,k,l
      real :: x1,x2,xcoast

!-----------------------------------------------------------------------
!  BEGIN user settings

    if( (sfcmodel.ge.1) .or. (bbc.eq.3) )then

      ! Initialize the surface:

      ! Key:
      ! tsk      = "skin temperature" (K) of soil/water  (~1 cm deep)
      ! tmn      = deep-layer temperature (K) of soil  (sfcmodel=2,3,4,6,7 only)
      ! xland    = land/water flag:   1 for land,   2 for water
      ! lu_index = land use index   (1-33) (see LANDUSE.TBL file)

!----------
!  Deafult ... fill in arrays with values from Namelist:

      do j=jb,je
      do i=ib,ie
        tsk(i,j)      = tsk0
      enddo
      enddo

      do j=jbl,jel
      do i=ibl,iel
        tmn(i,j)      = tmn0
        xland(i,j)    = xland0
        lu_index(i,j) = lu0
      enddo
      enddo

!----------
!  Initialize the sea breeze test case from WRF:

      if( initsfc.eq.2 )then 

        ! sea breeze test case from WRF:
        ! (modified on 150830 to match WRF3.7 setup)

        x1 = minx + (3.0/8.0)*(maxx-minx)
        x2 = minx + (5.0/8.0)*(maxx-minx)

        do j=jb,je
        do i=ib,ie
          if( xh(i).ge.x1 .and. xh(i).le.x2 )then
            ! land:
            xland(i,j) = 1.0
            lu_index(i,j) = 18
            tsk(i,j) = 280.0
          else
            ! water:
            xland(i,j) = 2.0
            lu_index(i,j) = 16
            tsk(i,j) = 287.0
          endif
          tmn(i,j) = tsk(i,j)-0.5
        enddo
        enddo

!----------

      elseif( initsfc.eq.3 ) then

        do j=jb,je
        do i=ib,ie
          if( xh(i).le.centerx )then
            znt(i,j) = 0.010
          else
            znt(i,j) = 0.001
          endif
        enddo
        enddo

!----------

      elseif( initsfc.eq.4 ) then

        xcoast = minx + 0.5*(maxx-minx)

        do j=jb,je
        do i=ib,ie
          if( xh(i) .lt. xcoast )then
            xland(i,j) = 1.0
            lu_index(i,j) = 18
          else
            xland(i,j) = 2.0
            lu_index(i,j) = 16
          endif
        enddo
        enddo

!----------
!  if initsfc is not 1,2:

      elseif( initsfc.ne.1 .and. initsfc.ne.2 )then

        ! build your own initial conditions here:

!!!        do j=jb,je
!!!        do i=ib,ie
!!!          tsk(i,j)      = ?
!!!          tmn(i,j)      = ?
!!!          xland(i,j)    = ?
!!!          lu_index(i,j) = ?
!!!        enddo
!!!        enddo

!----------

      endif     ! endif for initsfc options


!     END of user settings
!-----------------------------------------------------------------------
!     NO NOT CHANGE ANYTHING BELOW HERE

    if( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.4 .or. sfcmodel.eq.6 .or. sfcmodel.eq.7 )then

      if (num_soil_layers.NE.1)then
      slab_dzs(1)=.01
      slab_zs(1)=.5*slab_dzs(1)
      do l=2,num_soil_layers
         slab_dzs(l)=2*slab_dzs(l-1)
         slab_zs(l)=slab_zs(l-1)+.5*slab_dzs(l-1)+.5*slab_dzs(l)
      enddo
      do l=1,num_soil_layers
         do j=jb,je
            do i=ib,ie
              tslb(i,j,l)=( tsk(i,j)*(slab_zs(num_soil_layers)-slab_zs(l))   &
                          + tmn(i,j)*(slab_zs(l)-slab_zs(1)) ) / &
                          ( slab_zs(num_soil_layers)-slab_zs(1) )
            enddo
         enddo
      enddo
      endif

    endif

    if( sfcmodel.ge.1 )then

      if( axisymm.eq.1 .or. ny.eq.1 )then

        do j=jb,je
        do i=ib,ie
          dsxy(i,j) = dx*ruh(i)
        enddo
        enddo

      elseif( nx.eq.1 )then

        do j=jb,je
        do i=ib,ie
          dsxy(i,j) = dy*rvh(j)
        enddo
        enddo

      else

        do j=jb,je
        do i=ib,ie
!!!          dsxy(i,j) = sqrt( (dx*ruh(i))**2 + (dy*rvh(j))**2 )
          dsxy(i,j) = max( dx*ruh(i) , dy*rvh(j) )
        enddo
        enddo

      endif

    endif

  endif    ! endif for sfcmodel >= 1

!-----------------------------------------------------------------------
! ensure interoperability of surface schemes and turbulence schemes:
!   NOTE:  DO NOT modify anything here unless you really, really
!          know what you are doing

      ! parameters needed in turb subroutines:
      dosfcflx = .false.

      if( cm1setup.ge.1 .or. ipbl.eq.2 )then
        if( isfcflx.eq.1 ) dosfcflx = .true.
        if( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.4 .or. sfcmodel.eq.6 .or. sfcmodel.eq.7 )then
          dosfcflx = .true.
        endif
      endif
      if( ipbl.eq.1 .or. ipbl.eq.3 .or. ipbl.eq.4 .or. ipbl.eq.5 )then
        ! sfc flux/drag is handled by ysu/gfsedmf/mynn pbl subroutines:
        dosfcflx = .false.
      endif

      ! 230424: bug fix
      if( cm1setup.eq.4 .and. isfcflx.eq.1 ) dosfcflx = .true.

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Settings used in this simulation:'
      if(dowr) write(outfile,*) '  isfcflx  = ',isfcflx
      if(dowr) write(outfile,*) '  bbc      = ',bbc
      if(dowr) write(outfile,*) '  dosfcflx = ',dosfcflx
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!  Initialize land-surface parameters:

      if( sfcmodel.ge.1 .or. bbc.eq.3 )then
        ! initialize MM5/WRF landuse categories:
        call getlanduse(season,myid,ib,ie,jb,je,ibl,iel,jbl,jel,   &
                        lu_index,xland,emiss,thc,albd,znt,mavail)
      endif

!-----------------------------------------------------------------------
!  special cases:

      if( set_znt.eq.1 )then
        do j=jb,je
        do i=ib,ie
          znt(i,j) = max(1.0e-12,cnst_znt)
        enddo
        enddo
      endif

      if( testcase.eq.13 .and. (nx.gt.ny) )then
        ! UF Boundary Layer Wind Tunnel case:
        do j=jb,je
        do i=ib,ie
          if( xh(i).le.minx+(maxx-minx)/3.0 )then
            znt(i,j) = 0.011
          else
            znt(i,j) = 0.001
          endif
        enddo
        enddo
      endif

!-----------------------------------------------------------------------
!  reality checks:

      if( sfcmodel.ge.1 .or. bbc.eq.3 )then

        do j=jb,je
        do i=ib,ie

          ! 171020:
          if( xland(i,j).gt.1.5 ) mavail(i,j) = 1.0

          if( nint(xland(i,j)).le.0 .or. nint(xland(i,j)).ge.3 )then
            print *
            print *,'  Invalid value for xland '
            print *
            print *,'    xland must be 1 (for land) or 2 (for water) '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          endif
          if( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.4 .or. sfcmodel.eq.6 .or. sfcmodel.eq.7 )then
          if( lu_index(i,j).le.0 .or. lu_index(i,j).ge.34 )then
            print *
            print *,'  Invalid value for lu_index '
            print *
            print *,'    lu_index must be between 1 and 33 '
            print *,'    (see LANDUSE.TBL for more info) '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          endif
          endif
          if( tsk(i,j).lt.100.0 )then
            print *
            print *,'  WARNING:  tsk is less than 100 K '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          endif
          if( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.4 .or. sfcmodel.eq.6 .or. sfcmodel.eq.7 )then
          if( tmn(i,j).lt.100.0 )then
            print *
            print *,'  WARNING:  tmn is less than 100 K '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          endif
          endif
        enddo
        enddo

      endif

!-----------------------------------------------------------------------
!  Initialize the ocean mixed layer model:

      if(oceanmodel.eq.2)then
        call omlinit(outfile,oml_hml0, tsk,                     &
                       tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                       .false.,         .true. ,                &
                         1, ni+1,   1, nj+1,   1, nk+1,         &
                        ib, ie,  jb, je,  kb, ke,               &
                         1, ni,   1, nj,   1, nk                )
      endif

!-----------------------------------------------------------------------

      do j=jb,je
      do i=ib,ie
        rznt(i,j) = 1.0/max(1.0e-20,znt(i,j))
      enddo
      enddo

      if( sfcmodel.eq.7 )then
        do j=jb,je
        do i=ib,ie
          z0base(i,j) = znt(i,j)
        enddo
        enddo
        print *,'  znt,z0base = ',znt(1,1),z0base(1,1)
      endif

!-----------------------------------------------------------------------

      end subroutine init_surface


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getlanduse(season,myid,ib,ie,jb,je,ibl,iel,jbl,jel,   &
                            lu_index,xland,emiss,thc,albedo,znt,mavail)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: season,myid,ib,ie,jb,je,ibl,iel,jbl,jel
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: emiss,thc,albedo,mavail
      real, intent(inout), dimension(ib:ie,jb:je) :: znt,xland

      integer, parameter :: max_cats = 33    ! max categories
      integer, parameter :: max_seas =  2    ! max seasons

      integer, parameter :: iswater = 16  ! water must be category 16

      real, dimension(:,:), allocatable :: albd,slmo,sfem,sfz0,therin,scfx,sfhc
      integer :: cat,sea

      integer :: i,j,is,isn,ierr

      allocate(   albd(max_cats,max_seas) )
      allocate(   slmo(max_cats,max_seas) )
      allocate(   sfem(max_cats,max_seas) )
      allocate(   sfz0(max_cats,max_seas) )
      allocate( therin(max_cats,max_seas) )
      allocate(   scfx(max_cats,max_seas) )
      allocate(   sfhc(max_cats,max_seas) )

!-----------------------------------------------------------------------

      if( myid.eq.0 )then
        open(unit=11,file='LANDUSE.TBL',status='old',err=888)
        read(11,*,err=808)
        read(11,*,err=808)
        do isn=1,max_seas
          read(11,*,err=808)
          do is=1,max_cats
            read(11,*,err=808) i,albd(is,isn),slmo(is,isn),sfem(is,isn),sfz0(is,isn), &
                       therin(is,isn),scfx(is,isn),sfhc(is,isn)
          enddo
        enddo
        close(unit=11)
!        print *
!        print *,'  summer:'
!        print *,'  albd = ',albd(:,1)
!        print *,'  slmo = ',slmo(:,1)
!        print *,'  sfem = ',sfem(:,1)
!        print *,'  sfz0 = ',sfz0(:,1)
!        print *,'  ther = ',therin(:,1)
!        print *,'  scfx = ',scfx(:,1)
!        print *,'  sfhc = ',sfhc(:,1)
!        print *
!        print *,'  winter:'
!        print *,'  albd = ',albd(:,2)
!        print *,'  slmo = ',slmo(:,2)
!        print *,'  sfem = ',sfem(:,2)
!        print *,'  sfz0 = ',sfz0(:,2)
!        print *,'  ther = ',therin(:,2)
!        print *,'  scfx = ',scfx(:,2)
!        print *,'  sfhc = ',sfhc(:,2)
      endif

#ifdef MPI
      call MPI_BCAST(albd  ,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(slmo  ,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sfem  ,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sfz0  ,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(therin,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(scfx  ,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sfhc  ,max_cats*max_seas,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

!-----------------------------------------------------------------------

      ! ISN = season:    summer = 1    winter = 2
      ISN = season

      do j=jb,je
      do i=ib,ie
        IS = lu_index(i,j)
        ! SET NO-DATA POINTS (IS=0) TO WATER
        if(IS.EQ.0)then
          IS=ISWATER
        endif
        if( albd(is,isn).le.0.0 )then
          print *,'  category not found '
          call stopcm1
        endif
        ALBEDO(I,J) = ALBD(IS,ISN)/100.
        THC(i,j) = THERIN(IS,ISN)/100.
        EMISS(I,J) = SFEM(IS,ISN)
        ZNT(I,J) = SFZ0(IS,ISN)/100.
        MAVAIL(I,J) = SLMO(IS,ISN)
        if(IS.NE.ISWATER)then
          XLAND(I,J)=1.0
        else
          XLAND(I,J)=2.0
        endif
      enddo
      enddo

!-----------------------------------------------------------------------

      deallocate(   albd )
      deallocate(   slmo )
      deallocate(   sfem )
      deallocate(   sfz0 )
      deallocate( therin )
      deallocate(   scfx )
      deallocate(   sfhc )


      return

!-----------------------------------------------------------------------

808   print *
      print *,'  There was an error reading the LANDUSE.TBL file '
      print *
      call stopcm1

888   print *
      print *,'  There was an error opening the LANDUSE.TBL file '
      print *
      print *,'  Please make sure that LANDUSE.TBL is in the same directory '
      print *,'  as cm1.exe ... it is distributed with CM1 in the "run" directory '
      print *
      call stopcm1

      end subroutine getlanduse

  end module init_surface_module
