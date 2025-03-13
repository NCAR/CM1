  MODULE lsnudge_module

  ! code for large-scale nudging 
  ! (nudging of domain-average profiles of wind, temperature, and/or humidity)

  implicit none

  public


  ! Notes:
  !
  !   - See variables in the param19 section of namelist.input for settings related 
  !     to large-scale nudging.
  !
  !   - CM1 expects to find a series of "large-scale" data files in the same
  !     location as the cm1.exe file.  These must be named lsnudge_xxxx.dat, where 
  !     "xxxx" denotes an integer (eg, 0001, 0002, 0003, etc).
  !
  !   - For the contents of lsnudge_xxxx.dat files, see the example in the run directory.
  !     Note: the lsnudge_xxxx.dat files describe the domain-avg sounding that you want 
  !           to nudge towards. 
  !
  !   - If you only provide one lsnudge file (lsnudge_0001.dat), then lsnudge_time1 
  !     should be <= lsnudge_start, and lsnudge_time2 should be >= lsnudge_end.  
  !


  ! other variables (do not change):
    integer, parameter :: lsnudge_nmax  =  100
    integer :: lsnudge_unit,lsnudge_count,lsnudge_maxcount,lsnudge_maxlevels
    real :: lsnudgefac
    real, dimension(lsnudge_nmax) :: lsnudge_time1,lsnudge_time2

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine read_lsnudge(lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,timenudge,zh)
    use input
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(inout), dimension(lsnudge_maxlevels,lsnudge_nmax) :: lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv
    real, intent(in) :: timenudge
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh

    integer :: i,j,k,kk,kup,kdn,n,nloop,lsnudge_levels,iprint
    real :: interp_frac
    character(len=60) :: fname
    real, dimension(:), allocatable :: ztmp,utmp,vtmp,ttmp,qtmp
    logical :: added_level,domore

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Inside READ_LSNUDGE '
      if(dowr) write(outfile,*)

      if( terrain_flag )then
        if(myid.eq.0)then
        print *
        print *,'  cannot use lsnudge with terrain (for now)  '
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

      if( lsnudge_ramp_time .gt. (lsnudge_end-lsnudge_start) )then
        if( myid.eq.0 ) print *,'  12870: lsnudge_ramp_time is too large '
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

    iprint = 0
#ifdef MPI
    iprint = 2
    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
    if( myid.eq.iprint ) print *,'  iprint = ',iprint


    allocate( ztmp(lsnudge_maxlevels) )
    ztmp = 0.0
    allocate( utmp(lsnudge_maxlevels) )
    utmp = 0.0
    allocate( vtmp(lsnudge_maxlevels) )
    vtmp = 0.0
    allocate( ttmp(lsnudge_maxlevels) )
    ttmp = 0.0
    allocate( qtmp(lsnudge_maxlevels) )
    qtmp = 0.0

    nloop = 0
    lsnudge_count = 0

    ! GHB, 200214:  loop through all possible files at once. 

    doloop:  &
    DO while( nloop.lt.10000 )

      nloop = nloop+1

      domore = .false.

      ! open and read file (only processor 0 does this):
      myid0:  &
      if( myid.eq.0 )then

        fname = 'lsnudge_XXXX.dat'
        write(fname(9:12),101) nloop
101     format(i4.4)

        print *
        print *,' ------------------------------------------------------- '
        print *
        print *,'  LSNUDGE: nloop = ',nloop

        open(unit=lsnudge_unit,file=fname,status='old',err=9090)

        print *
        print *,'  Found fname = ',fname
        lsnudge_count = lsnudge_count+1

        read(lsnudge_unit,*)
        read(lsnudge_unit,*) lsnudge_time1(lsnudge_count),lsnudge_time2(lsnudge_count)
        read(lsnudge_unit,*)
        read(lsnudge_unit,*)

        print *,'  lsnudge_time1 = ',lsnudge_time1(lsnudge_count)
        print *,'  lsnudge_time2 = ',lsnudge_time2(lsnudge_count)
        print *

        if( nloop.ge.1000 )then
          if( myid.eq.0 )then
            print *
            print *,'  nloop = ',nloop
            print *
            print *,'  Can not find appropriate lsnudge file.  Stopping cm1 .... '
            print *
          endif
          call stopcm1
        endif

        lsnudge_levels = 0

        do k=1,lsnudge_maxlevels
          read(lsnudge_unit,*,err=8888,end=9999) ztmp(k),ttmp(k),qtmp(k),utmp(k),vtmp(k)
          ! convert qv from g/kg to kg/kg:
          qtmp(k) = 0.001*qtmp(k)
          lsnudge_levels = lsnudge_levels + 1
          if( lsnudge_levels .gt. lsnudge_maxlevels )then
            print *
            print *,'  lsnudge_maxlevels is too small '
            print *,'  (increase value and recompile) '
            print *
            print *,'  ... stopping cm1 ... '
            print *
            call stopcm1
          endif
        enddo

9999    print *,'  Found ',lsnudge_levels,' levels '

        close(unit=lsnudge_unit)

      if( do_lsnudge_u )then
        print *
        print *,'  ztmp(m),utmp(m/s): '
        do k=1,lsnudge_levels
          print *,k,ztmp(k),utmp(k)
        enddo
      endif

      if( do_lsnudge_v )then
        print *
        print *,'  ztmp(m),vtmp(m/s): '
        do k=1,lsnudge_levels
          print *,k,ztmp(k),vtmp(k)
        enddo
      endif

      if( do_lsnudge_th )then
        print *
        print *,'  ztmp(m),ttmp(K): '
        do k=1,lsnudge_levels
          print *,k,ztmp(k),ttmp(k)
        enddo
      endif

      if( do_lsnudge_qv )then
        print *
        print *,'  ztmp(m),qtmp(g/kg): '
        do k=1,lsnudge_levels
          print *,k,ztmp(k),1000.0*qtmp(k)
        enddo
      endif

        added_level = .false.

        if( ztmp(1).gt.0.0 )then
          added_level = .true.
          lsnudge_levels = lsnudge_levels + 1
          do k=lsnudge_levels,2,-1
            ztmp(k) = ztmp(k-1)
            utmp(k) = utmp(k-1)
            vtmp(k) = vtmp(k-1)
            ttmp(k) = ttmp(k-1)
            qtmp(k) = qtmp(k-1)
          enddo
          ztmp(1) = 0.0
          utmp(1) = utmp(2)-ztmp(2)*(utmp(3)-utmp(2))/(ztmp(3)-ztmp(2))
          vtmp(1) = vtmp(2)-ztmp(2)*(vtmp(3)-vtmp(2))/(ztmp(3)-ztmp(2))
          ttmp(1) = ttmp(2)-ztmp(2)*(ttmp(3)-ttmp(2))/(ztmp(3)-ztmp(2))
          qtmp(1) = qtmp(2)-ztmp(2)*(qtmp(3)-qtmp(2))/(ztmp(3)-ztmp(2))
        endif

        if( ztmp(lsnudge_levels).lt.maxz )then
          added_level = .true.
          lsnudge_levels = lsnudge_levels + 1
          ztmp(lsnudge_levels) = maxz
          utmp(lsnudge_levels) = utmp(lsnudge_levels-1)+(ztmp(lsnudge_levels)-ztmp(lsnudge_levels-1))  &
                                                     *(utmp(lsnudge_levels-1)-utmp(lsnudge_levels-2))  &
                                                     /(ztmp(lsnudge_levels-1)-ztmp(lsnudge_levels-2))
          vtmp(lsnudge_levels) = vtmp(lsnudge_levels-1)+(ztmp(lsnudge_levels)-ztmp(lsnudge_levels-1))  &
                                                     *(vtmp(lsnudge_levels-1)-vtmp(lsnudge_levels-2))  &
                                                     /(ztmp(lsnudge_levels-1)-ztmp(lsnudge_levels-2))
          ttmp(lsnudge_levels) = ttmp(lsnudge_levels-1)+(ztmp(lsnudge_levels)-ztmp(lsnudge_levels-1))  &
                                                     *(ttmp(lsnudge_levels-1)-ttmp(lsnudge_levels-2))  &
                                                     /(ztmp(lsnudge_levels-1)-ztmp(lsnudge_levels-2))
          qtmp(lsnudge_levels) = qtmp(lsnudge_levels-1)+(ztmp(lsnudge_levels)-ztmp(lsnudge_levels-1))  &
                                                     *(qtmp(lsnudge_levels-1)-qtmp(lsnudge_levels-2))  &
                                                     /(ztmp(lsnudge_levels-1)-ztmp(lsnudge_levels-2))
        endif

!        IF( added_level )THEN
!          print *
!          print *,'  After adding levels (at bottom, top, or both): '
!          print *,'  ztmp,utmp,vtmp: '
!          do k=1,lsnudge_levels
!            print *,k,ztmp(k),utmp(k),vtmp(k)
!          enddo
!          print *
!        ENDIF

      endif  myid0

      domore = .true.

9090  continue

      IF( myid.eq.0 .and. ( .not. domore ) )THEN
        print *
        print *,'  did not find fname = ',fname
        print *
      ENDIF

#ifdef MPI
      call MPI_BCAST(domore,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif


    do_mas:  &
    IF( domore ) THEN

#ifdef MPI
      ! communicate the new variables to all processors:
      call MPI_BCAST(lsnudge_levels,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_time1,lsnudge_nmax,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_time2,lsnudge_nmax,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ztmp(1),lsnudge_levels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(utmp(1),lsnudge_levels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vtmp(1),lsnudge_levels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ttmp(1),lsnudge_levels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qtmp(1),lsnudge_levels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif


      ! interpolate to the actual model levels:

      if( myid.eq.iprint ) print *
      if( myid.eq.iprint ) print *,'  Interpolating to model grid: '

      i = 1
      j = 1

!!!      if( myid.eq.iprint ) print *
!!!      if( myid.eq.iprint ) print *,'  myid,zk-1,z,zk,interp_frac:'

        DO k=1,nk

          kk = 1
          do while( ztmp(kk) .lt. zh(i,j,k) )
            kk = kk+1
          enddo
          kdn = kk-1
          kup = kk

          IF( kdn.lt.1 )THEN

            lsnudge_u(k,lsnudge_count) =  utmp(1)
            lsnudge_v(k,lsnudge_count) =  vtmp(1)

          ELSE

            interp_frac = (   zh(1,1,k) - ztmp(kdn) )   &
                        / ( ztmp( kup ) - ztmp(kdn) )

            lsnudge_u(k,lsnudge_count) =  utmp(kdn) + ( utmp(kup)- utmp(kdn))*interp_frac
            lsnudge_v(k,lsnudge_count) =  vtmp(kdn) + ( vtmp(kup)- vtmp(kdn))*interp_frac
            lsnudge_th(k,lsnudge_count) =  ttmp(kdn) + ( ttmp(kup)- ttmp(kdn))*interp_frac
            lsnudge_qv(k,lsnudge_count) =  qtmp(kdn) + ( qtmp(kup)- qtmp(kdn))*interp_frac

          ENDIF

!!!          if(myid.eq.iprint) write(outfile,*) '       ',ztmp(kdn),zh(i,j,k),ztmp(kup),interp_frac

        ENDDO


    if( do_lsnudge_u )then
      if( myid.eq.iprint ) print *
      if( myid.eq.iprint ) print *,'  zh(m), lsnudge_u(m/s):'
      if( myid.eq.iprint )then
        do k=1,nk
          if( myid.eq.iprint ) print *,k,zh(1,1,k),lsnudge_u(k,lsnudge_count)
        enddo
      endif
      if( myid.eq.iprint ) print *
    endif

    if( do_lsnudge_v )then
      if( myid.eq.iprint ) print *
      if( myid.eq.iprint ) print *,'  zh(m), lsnudge_v(m/s):'
      if( myid.eq.iprint )then
        do k=1,nk
          if( myid.eq.iprint ) print *,k,zh(1,1,k),lsnudge_v(k,lsnudge_count)
        enddo
      endif
      if( myid.eq.iprint ) print *
    endif

    if( do_lsnudge_th )then
      if( myid.eq.iprint ) print *
      if( myid.eq.iprint ) print *,'  zh(m), lsnudge_th(K):'
      if( myid.eq.iprint )then
        do k=1,nk
          if( myid.eq.iprint ) print *,k,zh(1,1,k),lsnudge_th(k,lsnudge_count)
        enddo
      endif
      if( myid.eq.iprint ) print *
    endif

    if( do_lsnudge_qv )then
      if( myid.eq.iprint ) print *
      if( myid.eq.iprint ) print *,'  zh(m), lsnudge_qv(kg/kg):'
      if( myid.eq.iprint )then
        do k=1,nk
          if( myid.eq.iprint ) print *,k,zh(1,1,k),lsnudge_qv(k,lsnudge_count)
        enddo
      endif
      if( myid.eq.iprint ) print *
    endif


    ELSE

      nloop = 10000000

    ENDIF  do_mas

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

    ENDDO  doloop

      if( myid.eq.iprint ) print *
      if( myid.eq.iprint ) print *,'  loop ended.  lsnudge_count = ',lsnudge_count
      if( myid.eq.iprint ) print *

      if( lsnudge_count .eq. 0 )then
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  Did not find any lsnudge files '
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  ... stopping CM1 ... '
        if(myid.eq.0) print *
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

      deallocate( ztmp )
      deallocate( utmp )
      deallocate( vtmp )
      deallocate( ttmp )
      deallocate( qtmp )

      lsnudge_maxcount = lsnudge_count
      lsnudge_count = 1

    ! all done
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Leaving READ_LSWIND_NUDGE '
      if(dowr) write(outfile,*)

    return

8888 print *
     print *,'  There was an error reading lsnudge_unit = ',lsnudge_unit
     print *
     print *,'  ... stopping cm1 ... '
     print *
     call stopcm1

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Leaving READ_LSWIND_NUDGE '
      if(dowr) write(outfile,*)

    end subroutine read_lsnudge

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE lsnudge_module
