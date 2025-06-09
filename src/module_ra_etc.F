
 !-------------------------------------------------------------------!
 !---  cm1r19:  misc stuff needed to adapt WRF RRTMG code to cm1  ---!
 !-------------------------------------------------------------------!

 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

module module_wrf_error
  integer           :: wrf_debug_level = 0
  character*256     :: wrf_err_message

  ! LOGICAL silence -- if TRUE (non-zero), this MPI rank does not send
  !   messages via wrf_message, end_timing, wrf_debug, atm_announce,
  !   cmp_announce, non-fatal glob_abort, or the like.  If FALSE, this
  !   MPI rank DOES send messages.  Regardless of this setting, fatal
  !   errors (wrf_error_fatal or fatal glob_aborts) and anything sent to
  !   write or print will be sent.
#if defined(DM_PARALLEL)
  integer, save :: silence=0
#else
  integer, parameter :: silence=0 ! Per-rank silence requires MPI
#endif

  ! LOGICAL buffered -- if TRUE, messages are buffered via clog_write.
  !   Once the buffer is full, messages are sent to stdout.  This does
  !   not apply to WRF_MESSAGE2, WRF_ERROR_FATAL, or anything sent to
  !   write or print.  The buffering implementation will not write
  !   partial lines, and buffer size is specified via namelist (see
  !   init_module_wrf_error).
  !   If FALSE, messages are send directly to WRITE.
  !
  !   This must be enabled at compile time by setting $WRF_LOG_BUFFERING

#if defined(WRF_LOG_BUFFERING)
  integer :: buffered=0
#else
  integer, parameter :: buffered=0 ! buffering disabled at compile time
#endif

  ! LOGICAL stderrlog -- if TRUE, messages are sent to stderr via
  !   write(0,...).  If FALSE, messages are not sent to stderr.
  !   This is set to FALSE automatically when buffering is enabled.

  ! Defaults: Non-MPI configurations and HWRF turn OFF stderr.
  !    MPI configurations other than HWRF turn ON stderr.

#if defined( DM_PARALLEL ) && ! defined( STUBMPI ) && !(HWRF == 1)
  integer :: stderrlog=1 ! 1/T = send to write(0,...) if buffered=0
#else
  integer :: stderrlog=0! 1/T = send to write(0,...) if buffered=0
#endif

  integer, parameter :: wrf_log_flush=0, wrf_log_set_buffer_size=1, &
                        wrf_log_write=2

  !NOTE: Make sure silence, buffered and stderrlog defaults here match
  ! the namelist defaults in init_module_wrf_error.

! min_allowed_buffer_size: requested buffer sizes smaller than this
! will simply result in disabling of log file buffering.  This number
! should be larger than any line WRF prints frequently.  If you set it 
! too small, the buffering code will still work.  However, any line 
! that is larger than the buffer may result in two writes: one for 
! the message and one for the end-of-line character at the end (if the
! message didn't already have one).
  integer, parameter :: min_allowed_buffer_size=200

!$OMP THREADPRIVATE (wrf_err_message)
contains

! ------------------------------------------------------------------------------

  logical function wrf_at_debug_level ( level )
    implicit none
    integer , intent(in) :: level
    wrf_at_debug_level = ( level .LE. wrf_debug_level )
    return
  end function wrf_at_debug_level

! ------------------------------------------------------------------------------

  subroutine init_module_wrf_error(on_io_server)
    implicit none
    logical,optional,intent(in) :: on_io_server
#if defined(DM_PARALLEL)
    logical, external :: wrf_dm_on_monitor
#endif
    logical :: compute_tasks_silent
    logical :: io_servers_silent
    integer :: buffer_size,iostat,stderr_logging
    namelist /logging/ buffer_size,compute_tasks_silent, &
                       io_servers_silent,stderr_logging

    ! MAKE SURE THE NAMELIST DEFAULTS MATCH THE DEFAULT VALUES
    ! AT THE MODULE LEVEL

    ! Default: original behavior.  No buffering, all ranks talk
    compute_tasks_silent=.false.
    io_servers_silent=.false.
    buffer_size=0

    ! MPI configurations default to stderr logging, except for HWRF.
    ! Non-MPI does not log to stderr.  (Note that fatal errors always
    ! are sent to both stdout and stderr regardless of config.)
#if defined( DM_PARALLEL ) && ! defined( STUBMPI ) && !(HWRF == 1)
    stderr_logging=1
#else
    stderr_logging=0
#endif
500 format(A)
    ! Open namelist.input using the same unit used by module_io_wrf 
    ! since we know nobody else will use that unit:
    open(unit=27, file="namelist.input", form="formatted", status="old")
    read(27,nml=logging,iostat=iostat)
    if(iostat /= 0) then
#if (DA_CORE!=1)
       call wrf_debug ( 1 , 'Namelist logging not found in namelist.input. ' )
       call wrf_debug ( 1 , ' --> Using registry defaults for variables in logging.' )
#else
       write(0,*) 'Namelist logging not found in namelist.input. Using registry defaults for variables in logging.'
       write(6,*) 'Namelist logging not found in namelist.input. Using registry defaults for variables in logging.'
#endif
#      ifdef _WIN32
          FLUSH(0)
#      endif
       close(27)
       return
    endif
    close(27)

#if defined(WRF_LOG_BUFFERING)
    ! Forbid small buffers.  See the comment above for min_allowed_buffer_size:
    if(buffer_size>=min_allowed_buffer_size) then
       call wrf_log_action(wrf_log_set_buffer_size,buffer_size,' ')
       buffered=1
    else
       buffered=0
    endif
#else
    if(buffer_size>=min_allowed_buffer_size) then
       write(0,500) 'Forcing disabling of buffering due to compile-time configuration.'
       write(6,500) 'Forcing disabling of buffering due to compile-time configuration.'
    endif
#endif

    stderrlog=stderr_logging
    if(buffered/=0 .and. stderrlog/=0) then
       write(0,500) 'Disabling stderr logging since buffering is enabled.'
       write(6,500) 'Disabling stderr logging since buffering is enabled.'
#      ifdef _WIN32
          FLUSH(0)
#      endif
       stderrlog=0
    endif

#if defined(DM_PARALLEL)
    silence=0
    if(present(on_io_server)) then
       if(on_io_server) then
          if(io_servers_silent) &
               silence=1
          return
       endif
    endif
    if(compute_tasks_silent) then
       if(wrf_dm_on_monitor()) then
          silence=0
       else
          silence=1
       endif
    endif
#endif
  end subroutine init_module_wrf_error

end module module_wrf_error

! ------------------------------------------------------------------------------
! ------------------------  GLOBAL SCOPE SUBROUTINES  --------------------------
! ------------------------------------------------------------------------------
#if defined(WRF_LOG_BUFFERING)
subroutine wrf_log_action( act,int,str )
  ! The underlying clog.c is not thread-safe, so this wrapper subroutine
  ! ensures that only one thread accesses clog.c at a time.

  ! NOTE: This routine only exists if WRF_LOG_BUFFERING is defined at
  ! compile time.
  use module_wrf_error
  implicit none
  integer, intent(in) :: int,act
  character(*), intent(in) :: str
!$OMP CRITICAL(wrf_log_action_critical)
  if(act==wrf_log_flush) then
     call clog_flush(int)
  elseif(act==wrf_log_set_buffer_size) then
     call clog_set_buffer_len(int)
  elseif(act==wrf_log_write) then
     call clog_write(int,str)
  endif
!$OMP END CRITICAL(wrf_log_action_critical)
end subroutine wrf_log_action
#endif
! ------------------------------------------------------------------------------

! wrf_message: ordinary message
!   Write to stderr if stderrlog=T to ensure immediate output
!   Write to stdout for buffered output.
subroutine wrf_message( str )
#ifdef ESMFIO
  use ESMF
#endif
  use module_wrf_error, only: silence, buffered, stderrlog, wrf_log_write
  implicit none

  character*(*) str
  if(silence/=0) return
  if(buffered/=0) then
#if defined(WRF_LOG_BUFFERING)
     call wrf_log_action(wrf_log_write,len_trim(str),str)
#endif
  else
!$OMP MASTER
     if(stderrlog/=0) then
300     format(A)
        write(0,300) trim(str)
# ifdef _WIN32
  FLUSH(0)
# endif
     endif
     print 300,trim(str)
!$OMP END MASTER
  endif

#ifdef ESMFIO
  call ESMF_LogWrite(TRIM(str),ESMF_LOGMSG_INFO)
#endif
end subroutine wrf_message

! ------------------------------------------------------------------------------

! Intentionally write to stderr only
! This is set to stderr, even in silent mode, because
! it is used for potentially fatal error or warning messages and
! we want the message to get to the log file before any crash 
! or MPI_Abort happens.
subroutine wrf_message2( str )
#ifdef ESMFIO
  use ESMF
#endif
  implicit none
  character*(*) str
!$OMP MASTER
400 format(A)
  write(0,400) str
# ifdef _WIN32
  FLUSH(0)
# endif
!$OMP END MASTER
#ifdef ESMFIO
  call ESMF_LogWrite(TRIM(str),ESMF_LOGMSG_INFO)
#endif
end subroutine wrf_message2

! ------------------------------------------------------------------------------

subroutine wrf_error_fatal3( file_str, line, str )
  use module_wrf_error
#ifdef ESMFIO
! 5.2.0r  USE ESMF_Mod
  use ESMF
#endif
  implicit none
  character*(*) file_str
  integer , intent (in) :: line  ! only print file and line if line > 0
  character*(*) str
  character*256 :: line_str

  write(line_str,'(i6)') line

  ! Fatal errors are printed to stdout and stderr regardless of
  ! any &logging namelist settings.

  call wrf_message( '-------------- FATAL CALLED ---------------' )
  ! only print file and line if line is positive
  if ( line > 0 ) then
    call wrf_message( 'FATAL CALLED FROM FILE:  '//file_str//'  LINE:  '//TRIM(line_str) )
  endif
  call wrf_message( str )
  call wrf_message( '-------------------------------------------' )

  force_stderr: if(stderrlog==0) then
  call wrf_message2( '-------------- FATAL CALLED ---------------' )
  ! only print file and line if line is positive
  if ( line > 0 ) then
        call wrf_message2( 'FATAL CALLED FROM FILE:  '//file_str//'  LINE:  '//TRIM(line_str) )
  endif
     call wrf_message2( trim(str) )
  call wrf_message2( '-------------------------------------------' )
  endif force_stderr

  ! Flush all streams.
  flush(6)
#if defined(WRF_LOG_BUFFERING)
  if(buffered/=0) call wrf_log_action(wrf_log_flush,1,' ')
# endif
  flush(0)

#ifdef ESMFIO
! 5.2.0r  CALL esmf_finalize(terminationflag=ESMF_ABORT)
  call esmf_finalize(endflag=ESMF_END_ABORT)
#endif

#ifdef TRACEBACKQQ
  call tracebackqq
#endif
  call wrf_abort
end subroutine wrf_error_fatal3

! ------------------------------------------------------------------------------

subroutine wrf_error_fatal( str )
  use module_wrf_error
  implicit none
  character*(*) str
  call wrf_error_fatal3 ( ' ', 0, str )
end subroutine wrf_error_fatal

! ------------------------------------------------------------------------------

! Check to see if expected value == actual value
! If not, print message and exit.  
subroutine wrf_check_error( expected, actual, str, file_str, line )
  use module_wrf_error
  implicit none
  integer , intent (in) :: expected
  integer , intent (in) :: actual
  character*(*) str
  character*(*) file_str
  integer , intent (in) :: line
  character (LEN=512)   :: rc_str
  character (LEN=512)   :: str_with_rc

  if ( expected .ne. actual ) then
    write (rc_str,*) '  Routine returned error code = ',actual
    str_with_rc = TRIM(str // rc_str)
    call wrf_error_fatal3 ( file_str, line, str_with_rc )
  endif
end subroutine wrf_check_error

! ------------------------------------------------------------------------------

!  Some compilers do not yet support the entirety of the Fortran 2003 standard.
!  This is a small patch to pick up the two most common events.  Most xlf 
!  compilers have an extension fflush.  That is available here.  For other older
!  compilers with no flush capability at all, we just stub it out completely.
!  These CPP ifdefs are defined in the configure file.

#ifdef USE_FFLUSH
subroutine flush ( iunit )
  implicit none
  integer :: iunit
  call fflush ( iunit ) 
end subroutine flush
#endif

#ifdef NO_FLUSH_SUPPORT
subroutine flush ( iunit )
  implicit none
  integer :: iunit
  return
end subroutine flush
#endif


 !-------------------------------------------------------------------!
 !---  cm1r19:  misc stuff needed to adapt WRF RRTMG code to cm1  ---!
 !-------------------------------------------------------------------!

 module module_state_description

   implicit none

    integer, parameter :: KFCUPSCHEME          =  10
    integer, parameter :: KFETASCHEME          =   1
    integer, parameter :: FER_MP_HIRES         =   5
    integer, parameter :: FER_MP_HIRES_ADVECT  =  15

    integer, parameter :: param_first_scalar = 1
    integer, parameter :: p_qc = 1
    integer, parameter :: p_qr = 1
    integer, parameter :: p_qi = 1
    integer, parameter :: p_qs = 1
    integer, parameter :: p_qg = 1
    integer, parameter :: p_qnc = 1
    integer, parameter :: p_qni = 1

 end module module_state_description

 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

 module module_model_constants

   implicit none

    !  240213: use the same values as CM1

   real    , parameter :: r_d          = 287.04
   real    , parameter :: cp           = 1005.7


   ! needed for mynn:
   real    , parameter ::  KARMAN=0.4
   real    , parameter :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
   real    , parameter :: p1000mb      = 100000.
   real    , parameter :: r_v          = 461.5
   real    , parameter :: rcp          = r_d/cp
   real    , parameter :: XLS          = 2.834E6
   real    , parameter :: XLV          = 2.501E6
   real    , parameter :: XLF          = XLS-XLV
   real    , parameter ::  SVP1=0.6112
   real    , parameter ::  SVP2=17.67
   real    , parameter ::  SVP3=29.65
   real    , parameter ::  SVPT0=273.15
   real    , parameter ::  EP_1=R_v/R_d-1.
   real    , parameter ::  EP_2=R_d/R_v
   real    , parameter :: rvovrd       = r_v/r_d
   real    , parameter :: cpv          = 4.*r_v
   real    , parameter :: cliq         = 4190.
   real    , parameter :: cice         = 2106.

   ! needed for myj:
   real , parameter ::  pq0=379.90516     ! 
   real , parameter ::  epsq2=0.2         ! initial TKE for camuw PBL scheme (m2 s^-2) 
   real , parameter ::  a2=17.2693882
   real , parameter ::  a3=273.16
   real , parameter ::  a4=35.86
   real , parameter ::  p608=rvovrd-1.

 end module module_model_constants

 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

subroutine set_wrf_debug_level ( level ) 
  use module_wrf_error 
  implicit none 
  integer , intent(in) :: level 
  wrf_debug_level = level 
  return 
end subroutine set_wrf_debug_level 
 
subroutine get_wrf_debug_level ( level ) 
  use module_wrf_error 
  implicit none 
  integer , intent(out) :: level 
  level = wrf_debug_level 
  return 
end subroutine get_wrf_debug_level 
 
subroutine wrf_debug( level , str ) 
  use module_wrf_error 
  implicit none 
  character*(*) str 
  integer , intent (in) :: level 
  integer :: debug_level 
  character (LEN=256) :: time_str 
  character (LEN=256) :: grid_str 
  character (LEN=512) :: out_str 
  if(silence/=0) return
  call get_wrf_debug_level( debug_level ) 
  if ( level .LE. debug_level ) then 
#ifdef _OPENMP 
  ! old behavior 
  call wrf_message( str ) 
#else 
  ! TBH: This fails on pgf90 6.1-4 when built with OpenMP and using more 
  ! TBH: than one thread. It works fine multi-threaded on AIX with xlf 
  ! TBH: 10.1.0.0 . Hence the cpp nastiness. 
  ! new behavior: include domain name and time-stamp 
!!!  CALL get_current_time_string( time_str ) 
!!!  CALL get_current_grid_name( grid_str ) 
  out_str = TRIM(grid_str)//' '//TRIM(time_str)//' '//TRIM(str) 
  call wrf_message( TRIM(out_str) ) 
#endif 
  endif 
  return 
end subroutine wrf_debug 

 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

!=========================================================================

! These are stub functions that do the right thing (usually nothing)
! in case DM_PARALLEL is not compiled for.
! This file, src/module_dm_stubs.F is copied to src/module_dm.F  when
! the code is built.
! If, on the other hand, a DM package is specified, the module_dm.F 
! provided with that package (e.g. RSL) is copied from /external/RSL/module_dm.F
! into src/module_dm.F.
! It is important to recognize this, because changes directly to src/module_dm.F
! will be lost!

logical function wrf_dm_on_monitor()
  wrf_dm_on_monitor = .true.
end function wrf_dm_on_monitor

integer function wrf_dm_monitor_rank()
  wrf_dm_monitor_rank = 0
end function wrf_dm_monitor_rank

subroutine wrf_get_myproc( myproc )
  implicit none
  integer myproc
  myproc = 0
  return
end subroutine wrf_get_myproc

subroutine wrf_get_nproc( nprocs )
  implicit none
  integer nprocs
  nprocs = 1
  return
end subroutine wrf_get_nproc

subroutine wrf_get_nprocx( nprocs )
  implicit none
  integer nprocs
  nprocs = 1
  return
end subroutine wrf_get_nprocx

subroutine wrf_get_nprocy( nprocs )
  implicit none
  integer nprocs
  nprocs = 1
  return
end subroutine wrf_get_nprocy

subroutine wrf_dm_bcast_string ( buf , size )
  implicit none
  integer size
  integer BUF(*)
  return
end subroutine wrf_dm_bcast_string

subroutine wrf_dm_bcast_bytes ( buf , size )
  implicit none
  integer size
  integer BUF(*)
  return
end subroutine wrf_dm_bcast_bytes

subroutine wrf_dm_bcast_integer( BUF, N1 )
   implicit none
   integer n1
   integer  buf(*)
   return
end subroutine wrf_dm_bcast_integer

subroutine wrf_dm_bcast_real( BUF, N1 )
   implicit none
   integer n1
   real  buf(*)
   return
end subroutine wrf_dm_bcast_real

subroutine wrf_dm_bcast_logical( BUF, N1 )
   implicit none
   integer n1
   logical  buf(*)
   return
end subroutine wrf_dm_bcast_logical

subroutine wrf_dm_halo ( domdesc , comms , stencil_id )
   implicit none
   integer domdesc , comms(*) , stencil_id
   return
end subroutine wrf_dm_halo

subroutine wrf_dm_boundary ( domdesc , comms , period_id , &
                             periodic_x , periodic_y )
   implicit none
   integer domdesc , comms(*) , period_id
   logical , intent(in)      :: periodic_x, periodic_y
   return
end subroutine wrf_dm_boundary

subroutine wrf_dm_xpose_z2x ( domdesc , comms , xpose_id  )
   implicit none
   integer domdesc , comms(*), xpose_id
   return
end subroutine wrf_dm_xpose_z2x
subroutine wrf_dm_xpose_x2y ( domdesc , comms , xpose_id  )
   implicit none
   integer domdesc , comms(*), xpose_id
   return
end subroutine wrf_dm_xpose_x2y
subroutine wrf_dm_xpose_y2z ( domdesc , comms , xpose_id  )
   implicit none
   integer domdesc , comms(*), xpose_id
   return
end subroutine wrf_dm_xpose_y2z

!-------

subroutine wrf_get_dm_communicator ( communicator )
   implicit none
   integer , intent(out) :: communicator
   communicator = 0
   return
end subroutine wrf_get_dm_communicator

subroutine wrf_get_dm_iocommunicator ( iocommunicator )
   implicit none
   integer , intent(out) :: iocommunicator
   iocommunicator = 0
   return
end subroutine wrf_get_dm_iocommunicator

subroutine wrf_dm_shutdown
      return
end subroutine wrf_dm_shutdown
subroutine wrf_abort
      stop 'wrf_abort'
end subroutine wrf_abort

subroutine wrf_patch_to_global_real (buf,globbuf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   real globbuf(*)
   real buf(*)
   return
end subroutine wrf_patch_to_global_real

subroutine wrf_global_to_patch_real (globbuf,buf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   real globbuf(*)
   real buf(*)
   return
end subroutine wrf_global_to_patch_real


subroutine wrf_patch_to_global_double (buf,globbuf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   double precision globbuf(*)
   double precision buf(*)
   return
end subroutine wrf_patch_to_global_double

subroutine wrf_global_to_patch_double (globbuf,buf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   double precision globbuf(*)
   double precision buf(*)
   return
end subroutine wrf_global_to_patch_double

subroutine wrf_patch_to_global_integer (buf,globbuf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   integer globbuf(*)
   integer buf(*)
   return
end subroutine wrf_patch_to_global_integer

subroutine wrf_global_to_patch_integer (globbuf,buf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   integer globbuf(*)
   integer buf(*)
   return
end subroutine wrf_global_to_patch_integer

subroutine wrf_patch_to_global_logical (buf,globbuf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   logical globbuf(*)
   logical buf(*)
   return
end subroutine wrf_patch_to_global_logical

subroutine wrf_global_to_patch_LOGICAL (globbuf,buf,domdesc,ndim,&
                                       ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe )
   implicit none
   integer                             ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
   integer fid,domdesc,ndim,glen(3),llen(3)
   logical globbuf(*)
   logical buf(*)
   return
end subroutine wrf_global_to_patch_LOGICAL

#if ( HWRF == 1 )
   subroutine hwrf_coupler_init
   end subroutine hwrf_coupler_init
#endif

   subroutine push_communicators_for_domain( id )
      implicit none
      integer, optional, intent(in) :: id   ! if specified also does an instate for grid id
   end subroutine push_communicators_for_domain
   subroutine pop_communicators_for_domain
   end subroutine pop_communicators_for_domain
   subroutine instate_communicators_for_domain( id )
      implicit none
      integer, intent(in) :: id
   end subroutine instate_communicators_for_domain
   subroutine store_communicators_for_domain( id )
      implicit none
      integer, intent(in) :: id
   end subroutine store_communicators_for_domain

!------------------------------------------------

integer function get_unused_unit()

    implicit none

    integer, parameter :: min_unit_number = 30
    integer, parameter :: max_unit_number = 99

    logical :: opened

    do get_unused_unit = min_unit_number, max_unit_number
       INQUIRE(UNIT=get_unused_unit, OPENED=opened)
       if ( .NOT. opened ) return
    end do

    get_unused_unit = -1

    return

end function get_unused_unit

 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!


!-----------------------------------------------------------------------
!
! Contains all subroutines to get time-varying mixing ratios of CO2, 
! N2O, CH4, CFC11 and CFC12 into radiation schemes.
! These subroutines enable the user to specify the mixing ratios
! in the file CAMtr_volume_mixing_ratio, giving the user an easy way
! to specify trace gases concentrations.
!
! The subroutines were first developed by and put together in this module
! by C. Carouge.
!-----------------------------------------------------------------------
module module_ra_clWRF_support

  use module_wrf_error

  implicit none
  private
  integer, parameter               :: r8 = 8
  integer, parameter               :: cyr = 1800     ! maximum num. of lines in 'CAMtr_volume_mixing_ratio' file
  integer, dimension(1:cyr), save  :: yrdata
  real, dimension(1:cyr), save  :: juldata
  real(r8), dimension(1:cyr), save :: co2r, n2or, ch4r, cfc11r, cfc12r


  ! Same values as in module_ra_cam_support.F
  integer, parameter               :: VARCAM_in_years = 233
  integer                          :: yr_CAM(1:VARCAM_in_years) = &
       (/ 1869, 1870, 1871, 1872, 1873, 1874, 1875, &
       1876, 1877, 1878, 1879, 1880, 1881, 1882, &
       1883, 1884, 1885, 1886, 1887, 1888, 1889, &
       1890, 1891, 1892, 1893, 1894, 1895, 1896, &
       1897, 1898, 1899, 1900, 1901, 1902, 1903, &
       1904, 1905, 1906, 1907, 1908, 1909, 1910, &
       1911, 1912, 1913, 1914, 1915, 1916, 1917, &
       1918, 1919, 1920, 1921, 1922, 1923, 1924, &
       1925, 1926, 1927, 1928, 1929, 1930, 1931, &
       1932, 1933, 1934, 1935, 1936, 1937, 1938, &
       1939, 1940, 1941, 1942, 1943, 1944, 1945, &
       1946, 1947, 1948, 1949, 1950, 1951, 1952, &
       1953, 1954, 1955, 1956, 1957, 1958, 1959, &
       1960, 1961, 1962, 1963, 1964, 1965, 1966, &
       1967, 1968, 1969, 1970, 1971, 1972, 1973, &
       1974, 1975, 1976, 1977, 1978, 1979, 1980, &
       1981, 1982, 1983, 1984, 1985, 1986, 1987, &
       1988, 1989, 1990, 1991, 1992, 1993, 1994, &
       1995, 1996, 1997, 1998, 1999, 2000, 2001, &
       2002, 2003, 2004, 2005, 2006, 2007, 2008, &
       2009, 2010, 2011, 2012, 2013, 2014, 2015, &
       2016, 2017, 2018, 2019, 2020, 2021, 2022, &
       2023, 2024, 2025, 2026, 2027, 2028, 2029, &
       2030, 2031, 2032, 2033, 2034, 2035, 2036, &
       2037, 2038, 2039, 2040, 2041, 2042, 2043, &
       2044, 2045, 2046, 2047, 2048, 2049, 2050, &
       2051, 2052, 2053, 2054, 2055, 2056, 2057, &
       2058, 2059, 2060, 2061, 2062, 2063, 2064, &
       2065, 2066, 2067, 2068, 2069, 2070, 2071, &
       2072, 2073, 2074, 2075, 2076, 2077, 2078, &
       2079, 2080, 2081, 2082, 2083, 2084, 2085, &
       2086, 2087, 2088, 2089, 2090, 2091, 2092, &
       2093, 2094, 2095, 2096, 2097, 2098, 2099, &
       2100, 2101                               /)
  real                             :: co2r_CAM(1:VARCAM_in_years) = &
       (/ 289.263, 289.263, 289.416, 289.577, 289.745, 289.919, 290.102, &
       290.293, 290.491, 290.696, 290.909, 291.129, 291.355, 291.587, 291.824, &
       292.066, 292.313, 292.563, 292.815, 293.071, 293.328, 293.586, 293.843, &
       294.098, 294.35, 294.598, 294.842, 295.082, 295.32, 295.558, 295.797,   &
       296.038, 296.284, 296.535, 296.794, 297.062, 297.338, 297.62, 297.91,   &
       298.204, 298.504, 298.806, 299.111, 299.419, 299.729, 300.04, 300.352,  &
       300.666, 300.98, 301.294, 301.608, 301.923, 302.237, 302.551, 302.863,  &
       303.172, 303.478, 303.779, 304.075, 304.366, 304.651, 304.93, 305.206,  &
       305.478, 305.746, 306.013, 306.28, 306.546, 306.815, 307.087, 307.365,  &
       307.65, 307.943, 308.246, 308.56, 308.887, 309.228, 309.584, 309.956,   &
       310.344, 310.749, 311.172, 311.614, 312.077, 312.561, 313.068, 313.599, &
       314.154, 314.737, 315.347, 315.984, 316.646, 317.328, 318.026, 318.742, &
       319.489, 320.282, 321.133, 322.045, 323.021, 324.06, 325.155, 326.299,  &
       327.484, 328.698, 329.933, 331.194, 332.499, 333.854, 335.254, 336.69,  &
       338.15, 339.628, 341.125, 342.65, 344.206, 345.797, 347.397, 348.98,    &
       350.551, 352.1, 354.3637, 355.7772, 357.1601, 358.5306, 359.9046,       &
       361.4157, 363.0445, 364.7761, 366.6064, 368.5322, 370.534, 372.5798,    &
       374.6564, 376.7656, 378.9087, 381.0864, 383.2994, 385.548, 387.8326,    &
       390.1536, 392.523, 394.9625, 397.4806, 400.075, 402.7444, 405.4875,     &
       408.3035, 411.1918, 414.1518, 417.1831, 420.2806, 423.4355, 426.6442,   &
       429.9076, 433.2261, 436.6002, 440.0303, 443.5168, 447.06, 450.6603,     &
       454.3059, 457.9756, 461.6612, 465.3649, 469.0886, 472.8335, 476.6008,   &
       480.3916, 484.2069, 488.0473, 491.9184, 495.8295, 499.7849, 503.7843,   &
       507.8278, 511.9155, 516.0476, 520.2243, 524.4459, 528.7127, 533.0213,   &
       537.3655, 541.7429, 546.1544, 550.6005, 555.0819, 559.5991, 564.1525,   &
       568.7429, 573.3701, 578.0399, 582.7611, 587.5379, 592.3701, 597.2572,   &
       602.1997, 607.1975, 612.2507, 617.3596, 622.524, 627.7528, 633.0616,    &
       638.457, 643.9384, 649.505, 655.1568, 660.8936, 666.7153, 672.6219,     &
       678.6133, 684.6945, 690.8745, 697.1569, 703.5416, 710.0284, 716.6172,   &
       723.308, 730.1008, 736.9958, 743.993, 751.0975, 758.3183, 765.6594,     &
       773.1207, 780.702, 788.4033, 796.2249, 804.1667, 812.2289, 820.4118,    &
       828.6444, 828.6444 /)

  public :: read_CAMgases

contains
  
  subroutine read_CAMgases(yr, julian, model, co2vmr, n2ovmr, ch4vmr, cfc11vmr, cfc12vmr) 

    integer, intent(in)            :: yr
    real, intent(in)               :: julian
    character(LEN=*), intent(in)   :: model           ! Radiation scheme name
    real(r8), optional, intent(out)    :: co2vmr, n2ovmr, ch4vmr, cfc11vmr, cfc12vmr

!Local
    
    integer                                          :: yearIN, found_yearIN, iyear  &
                                                       ,yr1,yr2
    integer                                         :: mondata(1:cyr)
    logical, external                                :: wrf_dm_on_monitor
    integer, external                                :: get_unused_unit
      
    integer                                          :: istatus, iunit, idata        
    integer, save                             :: VARCAM_in_years
    integer                                          :: nyrm, nyrp, njulm, njulp
    logical                                          :: exists
    logical, save                                    :: READtrFILE=.FALSE.
    character(LEN=256)                               :: message
    integer                                   :: monday(13)=(/0,31,28,31,30,31,30,31,31,30,31,30,31/)
    integer                                   :: mondayi(13)
    integer                                   :: my1,my2,my3, tot_valid
    
! CLWRF-UC June.09  (Copy from share/wrf_tsin.F)
    if ( .NOT. READtrFILE ) then
       READtrFILE= .TRUE.
       
       INQUIRE(FILE='CAMtr_volume_mixing_ratio', EXIST=exists)
       
       if (exists) then
          iunit = get_unused_unit()
          if ( iunit <= 0 ) then
             if ( wrf_dm_on_monitor() ) then
                call wrf_error_fatal('Error in module_ra_rrtm: could not find a free Fortran unit.')
             end if
          end if
          
          ! Read volume mixing ratio 
          open(UNIT=iunit, FILE='CAMtr_volume_mixing_ratio', FORM='formatted', &
               STATUS='old', IOSTAT=istatus)
          
          if (istatus == 0) then
             ! Ignore first two lines which constitute a header
             read(UNIT=iunit, FMT='(1X)')
             read(UNIT=iunit, FMT='(1X)')
             
             istatus = 0
             idata = 1
             do WHILE (istatus == 0)
                read(UNIT=iunit, FMT='(I4, 1x, F8.3,1x, 4(F10.3,1x))', IOSTAT=istatus)    &
                     yrdata(idata), co2r(idata), n2or(idata), ch4r(idata), cfc11r(idata), &
                     cfc12r(idata)
                if ( wrf_dm_on_monitor() ) then
                   write(message,*)'CLWRF reading...: istatus:',istatus,' idata:',idata,   &
                        ' year:', yrdata(idata), ' co2: ',co2r(idata), ' n2o: ',&
                        n2or(idata),' ch4:',ch4r(idata)
                   call wrf_debug( 0, message) 
                endif
                mondata(idata) = 6
                
                idata=idata+1
             end do

             if (istatus /= -1) then
                print *,'CLWRF -- clwrf -- CLWRF ALERT!'
                print *,"   Not normal ending of 'CAMtr_volume_mixing_ratio' file"
                print *,"   Lecture ends with 'IOSTAT'=",istatus
             end if
             VARCAM_in_years = idata - 1
             close(iunit)

             ! Calculate the julian day for each month.
             do idata=1,VARCAM_in_years
                mondayi = monday
                MY1=MOD(yrdata(idata),4)
                MY2=MOD(yrdata(idata),100)
                MY3=MOD(yrdata(idata),400)
                if(MY1.EQ.0.AND.MY2.NE.0.OR.MY3.EQ.0) mondayi(3)=29
                juldata(idata) = sum(mondayi(1:mondata(idata)))+real(mondayi(mondata(idata)+1))/2.-0.5   ! 1st Jan 00:00 = 0 julian day
             enddo
          endif
       else
          print *,"CLWRF -- clwrf -- CLWRF ALERT!"
          print *,"   'CAMtr_volume_mixing_ratio' does not exists"
          print *,'   Recovering original concentrations'
          
          !ccc Fake yrdata. Won't be use but a value > yr is needed so that
          ! everything work.
          ! With negative mixing ratios in co2r, n2or and ch4r, no interpolation
          ! is done and we use the constant mixing ratios that were used before
          ! CLWRF modifications.
          yrdata = yr + 1
          juldata = 1
          co2r   = -9999.999
          n2or   = -9999.999
          ch4r   = -9999.999
          cfc11r = -9999.999
          cfc12r = -9999.999

          ! For CAM model, recovers original data sets.
          if (model .eq. "CAM") then
             yrdata(1:VARCAM_in_years) = yr_CAM
             co2r(1:VARCAM_in_years)   = co2r_CAM
          endif
       endif ! CAMtr_volume_mixing_ratio exists
       
    endif ! File already opened and read 
    
    found_yearIN=0
    iyear=1
    !ccc Crash if iyear get > cyr (max. # of years in the mixing ratio file) ?
    !DO WHILE (found_yearIN == 0) 
    do WHILE (found_yearIN == 0 .and. iyear <= cyr)
       if (yrdata(iyear) .GT. yr )  then
          yearIN=iyear
          found_yearIN=1
       else if ((yrdata(iyear) .EQ. yr) .AND. (juldata(iyear) .GT. julian)) then
          yearIN=iyear
          found_yearIN = 1
       endif
       iyear=iyear+1
    enddo
    
    ! Prevent yr > last year in data
!    IF (yearIN .ge. VARCAM_in_years) yearIN=VARCAM_in_years-1
    if (iyear .ge. VARCAM_in_years) then
       yearIN=VARCAM_in_years-1
       found_yearIN = 1
    endif
    
    if (found_yearIN .NE. 0 ) then
       if (yearIN .eq. 1) yearIN = yearIN + 1    ! to take 2 first lines of the file
       nyrm = yrdata(yearIN-1)
       njulm = juldata(yearIN-1)
       nyrp = yrdata(yearIN)
       njulp = juldata(yearIN)
    endif

    if (PRESENT(co2vmr)) then
       co2vmr=-9999.999
       if (found_yearIN /= 0) then
          ! Interpolate data only if we have at least 2 valid concentrations.
          tot_valid = count(co2r(1:VARCAM_in_years) > 0)
          if (tot_valid >= 2 ) then
             call valid_years(yearIN, co2r, VARCAM_in_years,yr1, yr2)

             ! Set nyrm, njulm, nyrp, njulp
             nyrm = yrdata(yr1)
             njulm = juldata(yr1)
             nyrp = yrdata(yr2)
             njulp = juldata(yr2)

             call interpolate_CAMgases(yr, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, co2r, co2vmr)
          endif
       endif
       ! Verification of interpolated values. In case of no value 
       ! original values extracted from ghg_surfvals.F90 module

       if (co2vmr < 0. .or. found_yearIN == 0) then
          call orig_val("CO2",model,co2vmr)
       else
          ! If extrapolation, need to bound the data to pre-industrial concentrations
          if (co2vmr < 270.) co2vmr = 270.
          co2vmr=co2vmr*1.e-06
       end if
    endif

    if (PRESENT(n2ovmr)) then
       n2ovmr=-9999.999
       if (found_yearIN /= 0) then
          tot_valid = count(n2or(1:VARCAM_in_years) > 0)
          if (tot_valid >= 2 ) then
             call valid_years(yearIN, n2or, VARCAM_in_years,yr1, yr2)

             ! Set nyrm, njulm, nyrp, njulp
             nyrm = yrdata(yr1)
             njulm = juldata(yr1)
             nyrp = yrdata(yr2)
             njulp = juldata(yr2)

             call interpolate_CAMgases(yr, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, n2or, n2ovmr)
          endif
       endif
       
       if (n2ovmr < 0. .or. found_yearIN == 0) then
          call orig_val("N2O",model,n2ovmr)
       else
          ! If extrapolation, need to bound the data to pre-industrial concentrations
          if (n2ovmr < 270.) n2ovmr = 270.          
          n2ovmr=n2ovmr*1.e-09
       endif
       
    endif
    
    if (PRESENT(ch4vmr)) then
       ch4vmr=-9999.999
       if (found_yearIN /= 0) then
          tot_valid = count(ch4r(1:VARCAM_in_years) > 0)
          if (tot_valid >= 2 ) then
             call valid_years(yearIN, ch4r, VARCAM_in_years,yr1, yr2)

             ! Set nyrm, njulm, nyrp, njulp
             nyrm = yrdata(yr1)
             njulm = juldata(yr1)
             nyrp = yrdata(yr2)
             njulp = juldata(yr2)

             call interpolate_CAMgases(yr, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, ch4r, ch4vmr)
          endif
       endif
       
       if (ch4vmr < 0. .or. found_yearIN == 0) then
          call orig_val("CH4",model,ch4vmr)
       else
          ! If extrapolation, need to bound the data to pre-industrial concentrations
          if (ch4vmr < 700. ) ch4vmr = 700.
          ch4vmr=ch4vmr*1.e-09
       endif
    endif
    
    if (PRESENT(cfc11vmr)) then
       cfc11vmr = -9999.999
       if (found_yearIN /= 0) then
          tot_valid = count(cfc11r(1:VARCAM_in_years) > 0)
          if (tot_valid >= 2 ) then
             call valid_years(yearIN, cfc11r, VARCAM_in_years,yr1, yr2)

             ! Set nyrm, njulm, nyrp, njulp
             nyrm = yrdata(yr1)
             njulm = juldata(yr1)
             nyrp = yrdata(yr2)
             njulp = juldata(yr2)

             call interpolate_CAMgases(yr, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, cfc11r, cfc11vmr)
          endif
       endif
       
       if (cfc11vmr < 0. .or. found_yearIN == 0) then
          call orig_val("CFC11",model,cfc11vmr)
       else
          cfc11vmr=cfc11vmr*1.e-12
       endif
    endif
    
    if (PRESENT(cfc12vmr)) then
       cfc12vmr = -9999.999
       if (found_yearIN /= 0) then
          tot_valid = count(cfc12r(1:VARCAM_in_years) > 0)
          if (tot_valid >= 2 ) then
             call valid_years(yearIN, cfc12r, VARCAM_in_years,yr1, yr2)

             ! Set nyrm, njulm, nyrp, njulp
             nyrm = yrdata(yr1)
             njulm = juldata(yr1)
             nyrp = yrdata(yr2)
             njulp = juldata(yr2)

             call interpolate_CAMgases(yr, julian, nyrm, njulm, yr1, yr2, nyrp, njulp, cfc12r, cfc12vmr)
          endif
       endif
       
       if (cfc12vmr < 0. .or. found_yearIN == 0) then
          call orig_val("CFC12",model,cfc12vmr)
       else
          cfc12vmr=cfc12vmr*1.e-12
       endif
    endif
       
  end subroutine read_CAMgases
  
  subroutine valid_years(yearIN, gas, tot_years, yr1, yr2)

! Find 
    integer, intent(in) :: yearIN, tot_years
    integer, intent(out) :: yr2, yr1
    real(r8), intent(in) :: gas(:)

    ! Local variables
    integer :: yr_loc, idata
    

    yr_loc = yearIN
    yr2 = yr_loc
    yr1 = yr_loc-1
    
    ! If all valid dates are > yearIN then find the 2 lowest dates with
    ! valid data.
    if (count(gas(1:yr_loc-1) > 0.) == 0) then
!ccc       DO idata = yr_loc-1, tot_years-1
       do idata = yr_loc, tot_years-1
          if (gas(idata) > 0.) then
             yr1 = idata
             exit
          endif
       enddo
!ccc       DO idata = yr1, tot_years
       do idata = yr1+1, tot_years
          if (gas(idata) > 0.) then
             yr2 = idata
             exit
          endif
       enddo
    else  ! There is at least 1 valid year below yearIN
       if (gas(yr_loc) < 0.) then
          
          ! Find the closest valid year, look for higher year first
          if (any(gas(yr_loc:tot_years) > 0)) then
             do idata=yr_loc+1, tot_years
                if (gas(idata) > 0) then
                   yr2 = idata
                   exit
                endif
             enddo
          else  ! Need to take lower years and extrapolate data.
             do idata=yr_loc-1,2,-1
                if (gas(idata) > 0) then
                   yr2 = idata
                   exit
                endif
             enddo
          endif
       endif
       
       yr_loc = min(yr_loc-1, yr2-1)
       if (gas(yr_loc) < 0.) then
          
          ! Find the closest valid year, lower than yr1
          if (any(gas(1:yr_loc-1) > 0)) then
             do idata=yr_loc-1,1,-1
                if (gas(idata) > 0) then
                   yr1 = idata
                   exit
                endif
             enddo
          else  ! Need to take higher years and extrapolate data.
             print*, 'Problem: this case should never happen'
          endif
       else ! then yr1 is yr_loc (first valid date < yr2)
          yr1 = yr_loc
       endif
    endif
  end subroutine valid_years

  subroutine interpolate_CAMgases(yr, julian, yeari, juli, yr1, yr2, yearf, julf, gas, interp_gas)
    implicit none
! These subroutine interpolates a trace gas concentration from a non-homogeneously 
! distributed gas concentration evolution
    integer, intent (in)                      :: yr, yeari, yr1, yr2, yearf, juli, julf
    real, intent (in)                         :: julian
    real(r8), dimension(500), intent (in)         :: gas
    real(r8), intent (out)                    :: interp_gas
!Local
    real(r8)                                  :: yearini, yearend, gas1, gas2
    real                                      :: doymodel, doydatam, doydatap,    &
                                                 deltat, fact1, fact2, x
    integer                                   :: ny, my1,my2,my3,nday, maxyear, minyear
 
    
    ! Add support for leap-years

    ! Find smallest and largest year: yearf >= yeari since the file is ordered with increasing dates.
    minyear = MIN(yr,yeari)
    maxyear = MAX(yr,yearf)

    ! Initialise with the day in the year for each date.
    fact2 = juli
    fact1 = julf
    x = julian

    ! Calculate the julian day (day 0 = 1 Jan minyear at 00:00)
    do ny=minyear, maxyear-1
       nday = 365
       ! Leap-year?
       MY1=MOD(ny,4)
       MY2=MOD(ny,100)
       MY3=MOD(ny,400)
       if(MY1.EQ.0.AND.MY2.NE.0.OR.MY3.EQ.0) nday=366
       
       if (ny < yeari) then
          fact2 = fact2+nday
       endif
       if (ny < yr) then
          x = x+nday
       endif
       if (ny < yearf) then
          fact1 = fact1+nday
       endif
    enddo

    deltat = fact1-fact2
    fact1 = (fact1 - x)/deltat
    fact2 = (x - fact2)/deltat

    interp_gas = gas(yr1)*fact1+gas(yr2)*fact2 

    if (interp_gas .LT. 0. ) then
       interp_gas=-99999.
    endif
    
  end subroutine interpolate_CAMgases

  subroutine interpolate_lin(x1,y1,x2,y2,x0,y)
    implicit none
! Program to interpolate values y=a+b*x with:
!       a=y1
!       b=(y2-y1)/(x2-x1)
!       x=abs(x1-x0)

    real, intent (in)                       :: x1,x2,x0
    real(r8), intent (in)                   :: y1,y2
    real(r8), intent (out)                  :: y
    real(r8)                                :: a,b,x
    
    a=y1
    b=(y2-y1)/(x2-x1)
    
    if (x0 .GE. x1) then
       x=x0-x1
    else
       x=x1-x0
       b=-b
    endif
    
    y=a+b*x
    
  end subroutine interpolate_lin


  subroutine orig_val(tracer,model,out)

    character(LEN=*), intent(in) :: tracer  ! The trace gas name
    character(LEN=*), intent(in) :: model   ! The radiation scheme name
    real(r8), intent(inout)      :: out     ! The output value
!Local
    logical, external            :: wrf_dm_on_monitor
    character(LEN=256)           :: message


    if (model .eq. "CAM") then
       if (tracer .eq. "CO2") then
          out = 3.55e-4

       else if (tracer .eq. "N2O") then
          out = 0.311e-6

       else if (tracer .eq. "CH4") then
          out = 1.714e-6

       else if (tracer .eq. "CFC11") then
          out = 0.280e-9

       else if (tracer .eq. "CFC12") then
          out = 0.503e-9

       else
          if ( wrf_dm_on_monitor() ) then
             write(message,*) 'CLWRF : Trace gas ',tracer,' not valid for scheme ',model
             call wrf_error_fatal(message)
          endif
       endif

    else if (model .eq. "RRTM") then
       if (tracer .eq. "CO2") then
          out = 3.300e-4

       else if ((tracer .eq. "N2O") .or. &
                (tracer .eq. "CH4")) then
          out = 0.
       else
          if ( wrf_dm_on_monitor() ) then
             write(message,*) 'CLWRF : Trace gas ',tracer,' not valid for scheme ',model
             call wrf_error_fatal(message)
          endif
       endif

    else if (model .eq. "RRTMG") then
       if (tracer .eq. "CO2") then
          out = 379.e-6

       else if (tracer .eq. "N2O") then
          out = 319e-9

       else if (tracer .eq. "CH4") then
          out = 1774.e-9

       else if (tracer .eq. "CFC11") then
          out = 0.251e-9

       else if (tracer .eq. "CFC12") then
          out = 0.538e-9

       else
          if ( wrf_dm_on_monitor() ) then
             write(message,*) 'CLWRF : Trace gas ',tracer,' not valid for scheme ',model
             call wrf_error_fatal(message)
          endif
       endif

    else
       if ( wrf_dm_on_monitor() ) then
          write(message,*) 'CLWRF not implemented for the ',model,' radiative scheme.'
          call wrf_error_fatal(message)
       endif
    endif

  end subroutine orig_val

end module module_ra_clWRF_support


 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

module module_ra_rrtmg_wrf

    implicit none

    integer :: stderrlog=1

contains

! ------------------------------------------------------------------------------

subroutine wrf_error_fatal3( file_str, line, str )
!!!  USE module_wrf_error
#ifdef ESMFIO
! 5.2.0r  USE ESMF_Mod
  use ESMF
#endif
  implicit none
  character*(*) file_str
  integer , intent (in) :: line  ! only print file and line if line > 0
  character*(*) str
  character*256 :: line_str

  write(line_str,'(i6)') line

  ! Fatal errors are printed to stdout and stderr regardless of
  ! any &logging namelist settings.

  call wrf_message( '-------------- FATAL CALLED ---------------' )
  ! only print file and line if line is positive
  if ( line > 0 ) then
    call wrf_message( 'FATAL CALLED FROM FILE:  '//file_str//'  LINE:  '//TRIM(line_str) )
  endif
  call wrf_message( str )
  call wrf_message( '-------------------------------------------' )

  force_stderr: if(stderrlog==0) then
  call wrf_message2( '-------------- FATAL CALLED ---------------' )
  ! only print file and line if line is positive
  if ( line > 0 ) then
        call wrf_message2( 'FATAL CALLED FROM FILE:  '//file_str//'  LINE:  '//TRIM(line_str) )
  endif
     call wrf_message2( trim(str) )
  call wrf_message2( '-------------------------------------------' )
  endif force_stderr

  ! Flush all streams.
  flush(6)
#if defined(WRF_LOG_BUFFERING)
  if(buffered/=0) call wrf_log_action(wrf_log_flush,1,' ')
# endif
  flush(0)

#ifdef ESMFIO
! 5.2.0r  CALL esmf_finalize(terminationflag=ESMF_ABORT)
  call esmf_finalize(endflag=ESMF_END_ABORT)
#endif

#ifdef TRACEBACKQQ
  call tracebackqq
#endif
  call wrf_abort
end subroutine wrf_error_fatal3

! ------------------------------------------------------------------------------

subroutine wrf_error_fatal( str )
!!!  USE module_wrf_error
  implicit none
  character*(*) str
  call wrf_error_fatal3 ( ' ', 0, str )
end subroutine wrf_error_fatal

! ------------------------------------------------------------------------------

end module module_ra_rrtmg_wrf

 !-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-!

