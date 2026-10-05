  MODULE writeout_nc_module

  implicit none

  private
  public :: netcdf_prelim
#ifdef NETCDF
  public :: restart_prelim,writestat_nc,writepdata_nc,writehifrq_nc,writediag_nc,writeazim_nc,writecenter_nc,writemove_nc,disp_err
#endif

  CONTAINS

!-------------------------------------------------------------
!
!  This subroutine writes data in NetCDF files.
!
!  Code originally written by Daniel Kirshbaum
!  Code converted to netcdf4 (f90) by George Bryan, May 2013
!  Code last modified by George Bryan, 130910
!
!-------------------------------------------------------------


      subroutine netcdf_prelim(rtime,nwrite,fnum,ncid,time_index,qname,                      &
                               name_output,desc_output,unit_output,grid_output,cmpr_output,  &
                               xh,xf,yh,yf,xfref,yfref,sigma,sigmaf,zs,zh,zf,                &
                               d2d,ds,du,dv,                                                 &
                               dumz,dumx,dumy)
      use input
      use constants
#ifdef NETCDF
      use netcdf
#endif
      implicit none

      real, intent(in) :: rtime
      integer, intent(in) :: nwrite,fnum
      integer, intent(inout) :: ncid,time_index
      character(len=3), dimension(maxq), intent(in) :: qname
      character(len=60), intent(in), dimension(maxvars) :: desc_output
      character(len=40), intent(in), dimension(maxvars) :: name_output,unit_output
      character(len=1),  intent(in), dimension(maxvars) :: grid_output
      logical, intent(in), dimension(maxvars) :: cmpr_output
      real, dimension(ib:ie),   intent(in) :: xh
      real, dimension(ib:ie+1), intent(in) :: xf
      real, dimension(jb:je),   intent(in) :: yh
      real, dimension(jb:je+1), intent(in) :: yf
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, dimension(kb:ke)  , intent(in) :: sigma
      real, dimension(kb:ke+1), intent(in) :: sigmaf
      real, dimension(ib:ie,jb:je), intent(in) :: zs
      real, dimension(ib:ie,jb:je,kb:ke),   intent(in) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf
      real, dimension(ni,nj) :: d2d
      real, dimension(ni,nj,nk) :: ds
      real, dimension(ni+1,nj,nk) :: du
      real, dimension(ni,nj+1,nk) :: dv
      real, intent(inout), dimension(nz+1) :: dumz
      real, intent(inout), dimension(d2i) :: dumx
      real, intent(inout), dimension(d2j) :: dumy

#ifdef NETCDF
      integer i,j,k,n,irec

      integer :: cdfid    ! ID for the netCDF file to be created
      integer :: status,dimid,varid
      integer :: timeid,oneid,tfile,num_write
      integer :: xhid,yhid,zhid,xfid,yfid,zfid
      character(len=8) chid

!-------------------------------------------------------------
! Declare and set integer values for the netCDF dimensions 
!-------------------------------------------------------------

      integer :: ival,jval,kval,ivalp1,jvalp1,kvalp1
      real :: actual_time

      logical :: allinfo

      integer, parameter :: shuffle       = 1
      integer, parameter :: deflate       = 1
      integer, parameter :: deflate_level = 1

      integer :: chkx,chky,chkxp1,chkyp1

      logical, parameter :: stop_on_error = .true.

      if( myid.eq.0 ) print *,'  Entering netcdf_prelim '

!--------------------------------------------------------------
! Initializing some things
!--------------------------------------------------------------

      actual_time = rtime

!--------------------------------------------------------------
!  Write data to cdf file
!--------------------------------------------------------------

    do i=1,maxstring
      string(i:i) = ' '
    enddo

  IF( fnum.eq.51 )THEN
    IF(output_filetype.eq.1)THEN
      string = 'cm1out.nc'
    ELSEIF(output_filetype.eq.2)THEN
      string = 'cm1out_XXXXXX.nc'
      write(string(8:13),100) nwrite
    ELSEIF(output_filetype.eq.3)THEN
      string = 'cm1out_XXXXXX_YYYYYY.nc'
      write(string( 8:13),100) myid
      write(string(15:20),100) nwrite
    ELSE
      if(dowr) write(outfile,*) '  for netcdf output, output_filetype must be either 1,2, or 3 '
      call stopcm1
    ENDIF
  ELSEIF( fnum.eq.71 )THEN
    IF(output_filetype.eq.1)THEN
      string = 'cm1out_i.nc'
    ELSEIF(output_filetype.eq.2)THEN
      string = 'cm1out_XXXXXX_i.nc'
      write(string(8:13),100) nwrite
    ELSEIF(output_filetype.eq.3)THEN
      string = 'cm1out_XXXXXX_YYYYYY_i.nc'
      write(string( 8:13),100) myid
      write(string(15:20),100) nwrite
    ELSE
      if(dowr) write(outfile,*) '  for netcdf output, output_filetype must be either 1,2, or 3 '
      call stopcm1
    ENDIF
  ELSE
    print *,'  Unrecognized value for fnum in netcdf_prelim '
    call stopcm1
  ENDIF

      if(myid.eq.0) print *,string

100   format(i6.6)

!--------------------------------------------------------------
!  Dimensions of data:

    IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN
      ival = nx
      jval = ny
      chkx = nx
      chky = ny
      chkxp1 = nx+1
      chkyp1 = ny+1
    ELSEIF( output_filetype.eq.3 )THEN
      ival = ni
      jval = nj
      chkx = ni
      chky = nj
      chkxp1 = ni+1
      chkyp1 = nj+1
    ELSE
      print *,'  unrecognized value for output_filetype '
      call stopcm1
    ENDIF

    ivalp1 = ival+1
    jvalp1 = jval+1

    kval = min(maxk,nk)
    kvalp1 = min(maxk+1,nk+1)

!--------------------------------------------------------------
!  if this is the start of a file, then do this stuff:

    num_write = nwrite

    allinfo = .false.
    IF(num_write.eq.1) allinfo=.true.
    IF(output_filetype.ge.2)THEN
      allinfo=.true.
      num_write=1
    ENDIF

    IF( num_write.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set num_write to 1 and write all info.
      status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.num_write )then
          if(myid.eq.0) print *,'  tfile,num_write = ',tfile,num_write
          num_write = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        num_write = 1
      endif
    ENDIF

    time_index = num_write

    ifallinfo: IF(allinfo)THEN
!!!      print *,'  allinfo ... '

!-----------------------------------------------------------------------
!  BEGIN NEW:

      if( myid.eq.0 ) print *,'  calling nf90_create '
      if( myid.eq.0 ) print *,'    string = ',string
#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

!ORF set our 1d mesh arrays as our dimensions - this will preserve the mesh for
!analysis/viz programs that follow this convention (e.g., VAPOR3, VisIt, Paraview)
!Also include the axis attribute for VAPOR3, setting xh, yh, zh to the X, Y and Z axes.

      call disp_err( nf90_def_dim(ncid,'xh',ival,xhid) , .true. )
      call disp_err( nf90_def_dim(ncid,'yh',jval,yhid) , .true. )
      call disp_err( nf90_def_dim(ncid,'zh',kval,zhid) , .true. )
      call disp_err( nf90_def_dim(ncid,'xf',chkxp1,xfid) , .true. )
      call disp_err( nf90_def_dim(ncid,'yf',chkyp1,yfid) , .true. )
      call disp_err( nf90_def_dim(ncid,'zf',kvalp1,zfid) , .true. )

      call disp_err( nf90_def_dim(ncid,'time',nf90_unlimited,timeid) , .true. )
      call disp_err( nf90_def_dim(ncid,'one',1,oneid) , .true. )

    IF(icor.eq.1)THEN
      call disp_err( nf90_def_var(ncid,"f_cor",nf90_float,oneid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","Coriolis parameter") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","1/s") , .true. )
    ENDIF

      call disp_err( nf90_def_var(ncid,"ztop",nf90_float,oneid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","height at top of model") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )

      call disp_err( nf90_def_var(ncid,"time",nf90_float,timeid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","time since beginning of simulation") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

      call disp_err( nf90_def_var(ncid,"xh",nf90_float,xhid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

      call disp_err( nf90_def_var(ncid,"xf",nf90_float,xfid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location of staggered u grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

      call disp_err( nf90_def_var(ncid,"yh",nf90_float,yhid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

      call disp_err( nf90_def_var(ncid,"yf",nf90_float,yfid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location of staggered v grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

      call disp_err( nf90_def_var(ncid,"zh",nf90_float,zhid,varid) , .true. ) !ORF make zh
      call disp_err( nf90_put_att(ncid,varid,"long_name","nominal height of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

      call disp_err( nf90_def_var(ncid,"zf",nf90_float,zfid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","nominal height of staggered w grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )


    if( imove.ge.1 )then

      call disp_err( nf90_def_var(ncid,"umove",nf90_float,timeid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","umove") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"vmove",nf90_float,timeid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","vmove") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

    endif

!--------------------------------------------------------
!  Just to be sure:

        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inq_dimid(ncid,'xh',xhid) , .true. )
        call disp_err( nf90_inq_dimid(ncid,'yh',yhid) , .true. )
        call disp_err( nf90_inq_dimid(ncid,'zh',zhid) , .true. )
        call disp_err( nf90_inq_dimid(ncid,'xf',xfid) , .true. )
        call disp_err( nf90_inq_dimid(ncid,'yf',yfid) , .true. )
        call disp_err( nf90_inq_dimid(ncid,'zf',zfid) , .true. )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 2D and 3D vars:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! new for cm1r19:  use "_output" arrays, which are
    !                  defined in writeout.F

      varloop:  &
      do n = 1,n_out
        print *,'  n,name_output = ',n,name_output(n)
!!!        if( myid.eq.0 ) print *,'  n = ',n,trim(name_output(n))
        if(     grid_output(n).eq.'s' )then
!ORF keep original calls for reference
!         call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
          call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        elseif( grid_output(n).eq.'u' .and. fnum.eq.51 )then
!         call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/nip1id,njid,nkid,timeid/),varid) , .true. )
          call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/xfid,yhid,zhid,timeid/),varid) , .true. )
        elseif( grid_output(n).eq.'v' .and. fnum.eq.51 )then
!         call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/niid,njp1id,nkid,timeid/),varid) , .true. )
          call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/xhid,yfid,zhid,timeid/),varid) , .true. )
        elseif( grid_output(n).eq.'w' .and. fnum.eq.51 )then
!         call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/niid,njid,nkp1id,timeid/),varid) , .true. )
          call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/xhid,yhid,zfid,timeid/),varid) , .true. )
        elseif( grid_output(n).eq.'2' )then
          call disp_err( nf90_def_var(ncid,trim(name_output(n)),nf90_float,(/xhid,yhid,timeid/),varid) , .true. )
        else
          IF( fnum.eq.71 .and. ( grid_output(n).eq.'u' .or. grid_output(n).eq.'v' .or. grid_output(n).eq.'w' ) )THEN
            ! no error for u,v,w with fnum=71
            print *,'  fnum=71, skipping: ',trim(name_output(n))
          ELSE
            print *,'  Unknown setting for grid_output = ',grid_output(n)
            print *,'    67329 '
            call stopcm1
          ENDIF
        endif
        call disp_err( nf90_put_att(ncid,varid,"long_name",trim(desc_output(n))) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units",trim(unit_output(n))) , .true. )
#ifdef NCFPLUS
        if(     grid_output(n).eq.'s' )then
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        elseif( grid_output(n).eq.'u' .and. fnum.eq.51 )then
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkxp1,chky,1,1/)) , .true. )
        elseif( grid_output(n).eq.'v' .and. fnum.eq.51 )then
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chkyp1,1,1/)) , .true. )
        elseif( grid_output(n).eq.'w' .and. fnum.eq.51 )then
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        elseif( grid_output(n).eq.'2' )then
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        else
          IF( fnum.eq.71 .and. ( grid_output(n).eq.'u' .or. grid_output(n).eq.'v' .or. grid_output(n).eq.'w' ) )THEN
            ! no error for u,v,w with fnum=71
            print *,'  fnum=71, skipping: ',trim(name_output(n))
          ELSE
            print *,'  Unknown setting for grid_output = ',grid_output(n)
            print *,'    67330 '
            call stopcm1
          ENDIF
        endif
        if( cmpr_output(n) )then
!!!          if( myid.eq.0 ) print *,'    ... compressing ... '
          call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
        endif
#endif
      enddo  varloop

!--------------------------------------------------
      ! global attributes:

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

      call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
      call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'x_units',trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'x_label','x') , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'y_units',trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'y_label','y') , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'z_units',trim(aunit)) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'z_label','z') , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nx",nx) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ny",ny) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nz",nz) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"imoist",imoist) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"sgsmodel",sgsmodel) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"tconfig",tconfig) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"bcturbs",bcturbs) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ptype",ptype) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"wbc",wbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ebc",ebc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"sbc",sbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nbc",nbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"bbc",bbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"tbc",tbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"iorigin",iorigin) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"axisymm",axisymm) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"iptra",iptra) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"npt",npt) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"fcor",fcor) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"radopt",radopt) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dtrad" ,dtrad) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ctrlat",ctrlat) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ctrlon",ctrlon) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"year"  ,year) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"month" ,month) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"day"   ,day) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"hour"  ,hour) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"minute",minute) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"second",second) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"sfcmodel",sfcmodel) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"oceanmodel",oceanmodel) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ipbl",ipbl) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"iice",iice) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"idm",idm) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"idmplus",idmplus) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"numq",numq) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nql1",nql1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nql2",nql2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nqs1",nqs1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nqs2",nqs2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nnc1",nnc1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nnc2",nnc2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nzl1",nzl1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nzl2",nzl2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nvl1",nvl1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nvl2",nvl2) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_m",c_m) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_e1",c_e1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_e2",c_e2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_s",c_s) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgs1",cgs1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgs2",cgs2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgs3",cgs3) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgs1",dgs1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgs2",dgs2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgs3",dgs3) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgt1",cgt1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgt2",cgt2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgt3",cgt3) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgt1",dgt1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgt2",dgt2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgt3",dgt3) , .true. )

      call disp_err( nf90_enddef(ncid) , .true. )

! ... end of defs
!--------------------------------------------------
! begin data ... initial time ...

    IF(icor.eq.1)THEN
      call disp_err( nf90_inq_varid(ncid,'f_cor',varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,fcor) , .true. )
    ENDIF

      call disp_err( nf90_inq_varid(ncid,'ztop',varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,outunitconv*ztop) , .true. )

        call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
    if( output_filetype.eq.1 .or. output_filetype.eq.2 )then
      do i=1,nx
        dumx(i) = outunitconv*0.5*(xfref(i)+xfref(i+1))
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/nx/)) , .true. )
    elseif( output_filetype.eq.3 )then
      do i=1,ni
        dumx(i) = outunitconv*xh(i)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/ni/)) , .true. )
    endif

        call disp_err( nf90_inq_varid(ncid,'xf',varid) , .true. )
    if( output_filetype.eq.1 .or. output_filetype.eq.2 )then
      do i=1,nx+1
        dumx(i) = outunitconv*xfref(i)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/nx+1/)) , .true. )
    elseif( output_filetype.eq.3 )then
      do i=1,ni+1
        dumx(i) = outunitconv*xfref(myi1-1+i)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/ni+1/)) , .true. )
    endif

        call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
    if( output_filetype.eq.1 .or. output_filetype.eq.2 )then
      do j=1,ny
        dumy(j) = outunitconv*0.5*(yfref(j)+yfref(j+1))
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/ny/)) , .true. )
    elseif( output_filetype.eq.3 )then
      do j=1,nj
        dumy(j) = outunitconv*yh(j)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/nj/)) , .true. )
    endif

        call disp_err( nf90_inq_varid(ncid,'yf',varid) , .true. )
    if( output_filetype.eq.1 .or. output_filetype.eq.2 )then
      do j=1,ny+1
        dumy(j) = outunitconv*yfref(j)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/ny+1/)) , .true. )
    elseif( output_filetype.eq.3 )then
      do j=1,nj+1
        dumy(j) = outunitconv*yfref(myj1-1+j)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/nj+1/)) , .true. )
    endif

        call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
      if(terrain_flag)then
        do k=1,kval
          dumz(k) = outunitconv*sigma(k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kval/)) , .true. )
      else
        do k=1,kval
          dumz(k) = outunitconv*zh(1,1,k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kval/)) , .true. )
      endif

        call disp_err( nf90_inq_varid(ncid,'zf',varid) , .true. )
      if(terrain_flag)then
        do k=1,kvalp1
          dumz(k) = outunitconv*sigmaf(k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/)) , .true. )
      else
        do k=1,kvalp1
          dumz(k) = outunitconv*zf(1,1,k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/)) , .true. )
      endif

      ! ... end if info at initial time only

!----------------------------------------------------------

    ENDIF ifallinfo

      call disp_err( nf90_inq_varid(ncid,'time',varid) , stop_on_error )
      call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , stop_on_error )

    if( imove.ge.1 )then

      call disp_err( nf90_inq_varid(ncid,'umove',varid) , stop_on_error )
      call disp_err( nf90_put_var(ncid,varid,umove,(/time_index/)) , stop_on_error )

      call disp_err( nf90_inq_varid(ncid,'vmove',varid) , stop_on_error )
      call disp_err( nf90_put_var(ncid,varid,vmove,(/time_index/)) , stop_on_error )

    endif

#endif

      if( myid.eq.0 ) print *,'  Leaving netcdf_prelim '

    end subroutine netcdf_prelim


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef NETCDF


    ! cm1r18:  netcdf restart files
    subroutine restart_prelim(nrst,ncid,mtime,xfref,yfref,zh,zf,sigma,sigmaf,  &
                              qname,num_soil_layers,nrad2d,dumx,dumy,dumz,time_index)

    use input
    use constants
    use netcdf
    implicit none

    integer, intent(in) :: nrst
    integer, intent(inout) :: ncid
    double precision, intent(in) :: mtime
    real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
    real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
    real, intent(in), dimension(ib:ie,jb:je,kb:ke)   :: zh
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
    real, intent(in), dimension(kb:ke)   :: sigma
    real, intent(in), dimension(kb:ke+1) :: sigmaf
    character(len=3), intent(in), dimension(maxq) :: qname
    integer, intent(in) :: num_soil_layers,nrad2d
    real, intent(inout), dimension(d2i) :: dumx
    real, intent(inout), dimension(d2j) :: dumy
    real, intent(inout), dimension(nz+1) :: dumz
    integer, intent(inout) :: time_index

    integer :: i,j,k,n
    integer :: ival,jval,kval,ivalp1,jvalp1,kvalp1
    integer :: chkx,chky,chkxp1,chkyp1
    integer :: status,varid
    real :: actual_time
    integer :: niid,njid,nkid,nip1id,njp1id,nkp1id,timeid
    integer :: nbudgetid,numqid,plocid,nparcelsid
    character(len=8) :: text1
    integer :: num_write,tfile
    logical :: allinfo

    chkx = nx
    chky = ny

    chkxp1 = nx+1
    chkyp1 = ny+1

    time_index = 1

    actual_time = mtime

    do i=1,maxstring
      string(i:i) = ' '
    enddo

  IF(     restart_filetype.eq.1 )THEN
    string = 'cm1rst.nc'
  ELSEIF( restart_filetype.eq.2 )THEN
    string = 'cm1rst_XXXXXX.nc'
    write(string(8:13),100) nrst
100 format(i6.6)
  ELSE
    print *,'  invalid value for restart_filetype '
    stop 33443
  ENDIF

    if(myid.eq.0) print *,'  string = ',string

!--------------------------------------------------------------
!  Dimensions of data:

    ival = nx
    jval = ny
    kval = nk

    ivalp1 = ival+1
    jvalp1 = jval+1
    kvalp1 = kval+1

!--------------------------------------------------------------
!  if this is the start of a file, then do this stuff:

    num_write = nrst

    allinfo = .false.
    IF( num_write.eq.1 ) allinfo = .true.
    IF( restart_filetype.ge.2 )THEN
      allinfo = .true.
      num_write = 1
    ENDIF

    IF( num_write.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set num_write to 1 and write all info.
      status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.num_write )then
          if(myid.eq.0) print *,'  tfile,num_write = ',tfile,num_write
          num_write = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        num_write = 1
      endif
    ENDIF

    time_index = num_write

!-----------------------------------------------------------------------
!  BEGIN

    ifallinfo:  IF( allinfo )THEN

#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

      !------------------
      ! define dims:

      call disp_err( nf90_def_dim(ncid,'ni',ival,niid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nj',jval,njid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nk',kval,nkid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nip1',ivalp1,nip1id) , .true. )
      call disp_err( nf90_def_dim(ncid,'njp1',jvalp1,njp1id) , .true. )
      call disp_err( nf90_def_dim(ncid,'nkp1',kvalp1,nkp1id) , .true. )
      call disp_err( nf90_def_dim(ncid,'time',nf90_unlimited,timeid) , .true. )
      call disp_err( nf90_def_dim(ncid,'nbudget',nbudget,nbudgetid) , .true. )
      call disp_err( nf90_def_dim(ncid,'numq',numq,numqid) , .true. )
    if( iprcl.eq.1 )then
      n = 3
      call disp_err( nf90_def_dim(ncid,'nploc',n,plocid) , .true. )
      n = max( 1 , nparcels )
      call disp_err( nf90_def_dim(ncid,'nparcels',n,nparcelsid) , .true. )
    endif

      !------------------
      ! define vars:

      call disp_err( nf90_def_var(ncid,"time",nf90_float,timeid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","time since beginning of simulation") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )

      call disp_err( nf90_def_var(ncid,"xh",nf90_float,niid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"xf",nf90_float,nip1id,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location of staggered u grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"yh",nf90_float,njid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"yf",nf90_float,njp1id,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location of staggered v grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"zh",nf90_float,nkid,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","nominal height of scalar grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      call disp_err( nf90_def_var(ncid,"zf",nf90_float,nkp1id,varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","nominal height of staggered w grid points") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

      !--------------------------------------------------------
      !  Just to be sure:

      call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )

      !------------------
      ! vars:

      call disp_err( nf90_def_var(ncid,"cm1rversion",nf90_float,(/timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"nstep"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"srec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"sirec"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"urec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vrec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"wrec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nrec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"mrec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"prec"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"trecs"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"trecw"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"arecs"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"arecw"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"mvrec"   ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nwrite"  ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nwritet" ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nwritea" ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nwritemv",nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nwriteh" ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nrst"    ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nstatout",nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ndt"     ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"icenter" ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"jcenter" ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"old_format",nf90_int,(/timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"npt"     ,nf90_int,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"numparcels",nf90_int,(/timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"dt"      ,nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","timestep") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s") , .true. )

      call disp_err( nf90_def_var(ncid,"dtlast"  ,nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","previous timestep") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s") , .true. )

      call disp_err( nf90_def_var(ncid,"xcenter" ,nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ycenter" ,nf90_float,(/timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"umove"  ,nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"umovem1",nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"umovem2",nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vmove"  ,nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vmovem1",nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vmovem2",nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"domainlocx"   ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"domainlocy"   ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"adaptmovetim" ,nf90_double,(/timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"cflmax"  ,nf90_float,(/timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max Courant number from previous timestep") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( nf90_def_var(ncid,"mtime"   ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","model time (i.e., time since beginning of simulation)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s") , .true. )

      call disp_err( nf90_def_var(ncid,"stattim" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"taptim"  ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"rsttim"  ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"radtim"  ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"prcltim" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"adt"     ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"adtlast" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"acfl"    ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"dbldt"   ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"mass1"   ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcu" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcv" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcs" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcsu" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcsv" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfct" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcq" ,nf90_double,(/timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"avgsfcp" ,nf90_double,(/timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"qbudget" ,nf90_double,(/nbudgetid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"asq"     ,nf90_double,(/numqid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"bsq"     ,nf90_double,(/numqid,timeid/),varid) , .true. )

      call disp_err( nf90_def_var(ncid,"ug"      ,nf90_float,(/nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vg"      ,nf90_float,(/nkid,timeid/),varid) , .true. )

      !------------------
      ! 2d/3d vars:

      call disp_err( nf90_def_var(ncid,"rain"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","accumulated surface rainfall") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","cm") , .true. )

      call disp_err( nf90_def_var(ncid,"hail"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","accumulated surface hailfall") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","cm") , .true. )

      call disp_err( nf90_def_var(ncid,"sws"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max windspeed at lowest level") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"svs"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max vert vorticity at lowest level") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s^(-1)") , .true. )

      call disp_err( nf90_def_var(ncid,"sps"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","min pressure at lowest level") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( nf90_def_var(ncid,"srs"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface rainwater") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( nf90_def_var(ncid,"sgs"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface graupel/hail") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( nf90_def_var(ncid,"sus"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max w at 5 km AGL") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"shs"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max integrated updraft helicity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m2/s2") , .true. )

    if( nrain.eq.2 )then

      call disp_err( nf90_def_var(ncid,"rain2"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","accumulated surface rainfall (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","cm") , .true. )

      call disp_err( nf90_def_var(ncid,"hail2"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","accumulated surface hailfall (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","cm") , .true. )

      call disp_err( nf90_def_var(ncid,"sws2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max windspeed at lowest level (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"svs2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max vert vorticity at lowest level (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s^(-1)") , .true. )

      call disp_err( nf90_def_var(ncid,"sps2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","min pressure at lowest level (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( nf90_def_var(ncid,"srs2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface rainwater (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( nf90_def_var(ncid,"sgs2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max surface graupel/hail (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

      call disp_err( nf90_def_var(ncid,"sus2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max w at 5 km AGL (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"shs2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","max integrated updraft helicity (translated)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m2/s2") , .true. )

    endif

      call disp_err( nf90_def_var(ncid,"tsk"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","soil/ocean temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

      call disp_err( nf90_def_var(ncid,"rho"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","dry-air density") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/m^3") , .true. )

      call disp_err( nf90_def_var(ncid,"prs"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

      call disp_err( nf90_def_var(ncid,"ua"      ,nf90_float,(/nip1id,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","west-east velocity (at u points)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"va"      ,nf90_float,(/niid,njp1id,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","south-north velocity (at v points)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"wa"      ,nf90_float,(/niid,njid,nkp1id,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","vertical velocity (at w points)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

      call disp_err( nf90_def_var(ncid,"ppi"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","perturbation non-dimensional pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( nf90_def_var(ncid,"tha"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","perturbation potential temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

      call disp_err( nf90_def_var(ncid,"ppx"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","change in nondimensional pressure used for forward-time-weighting on small steps") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

    if( psolver.eq.6 )then
      call disp_err( nf90_def_var(ncid,"phi1",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","phi1") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

      call disp_err( nf90_def_var(ncid,"phi2",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","phi2") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )
    endif

    if( imoist.eq.1 )then
      do n=1,numq
        text1 = '        '
        write(text1(1:3),156) qname(n)
156     format(a3)
        call disp_err( nf90_def_var(ncid,text1     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        if( n.eq.nqv )then
          call disp_err( nf90_put_att(ncid,varid,"long_name","water vapor mixing ratio") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )
        elseif(idm.eq.1.and.n.ge.nnc1.and.n.le.nnc2)then
          call disp_err( nf90_put_att(ncid,varid,"long_name","number concentration") , .true. )
          IF ( output_n_as_m3 == 0 ) THEN
          call disp_err( nf90_put_att(ncid,varid,"units","kg^{-1}") , .true. )
          ELSE
          ! for NSSL microphysics
          call disp_err( nf90_put_att(ncid,varid,"units","m^{-3}") , .true. )
          ENDIF
        elseif(idm .eq. 1 .and. n.ge.nzl1 .and. n .le. nzl2)then
          call disp_err( nf90_put_att(ncid,varid,"long_name","reflectivity moment") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","Z m^{-3}kg^{-1}") , .true. )
        elseif(idm .eq. 1 .and. n.ge.nvl1 .and. n .le. nvl2)then
          call disp_err( nf90_put_att(ncid,varid,"long_name","particle volume") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","m^{3}kg^{-1}") , .true. )
        else
          call disp_err( nf90_put_att(ncid,varid,"long_name","mixing ratio") , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )
        endif
      enddo
    endif
    if(imoist.eq.1.and.eqtset.eq.2)then
      call disp_err( nf90_def_var(ncid,"qpten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","ppi tendency from microphysics on previous timestep (related to h_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","s^(-1)") , .true. )
      !---
      call disp_err( nf90_def_var(ncid,"qtten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","theta tendency from microphysics on previous timestep (related to h_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K/s") , .true. )
      !---
      call disp_err( nf90_def_var(ncid,"qvten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","qv tendency from microphysics on previous timestep (related to qv_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg/s") , .true. )
      !---
      call disp_err( nf90_def_var(ncid,"qcten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","qc tendency from microphysics on previous timestep (related to qc_diabatic in WRF)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg/s") , .true. )
    endif
    if(sgsmodel.eq.1.or.sgsmodel.eq.3.or.sgsmodel.eq.4.or.sgsmodel.eq.5)then
      call disp_err( nf90_def_var(ncid,"tkea"    ,nf90_float,(/niid,njid,nkp1id,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","subgrid turbulence kinetic energy") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m^2/s^2") , .true. )
    endif

    if( radopt.ge.1 )then
      call disp_err( nf90_def_var(ncid,"lwten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"swten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"thraten" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"radsw"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"rnflx"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"radswnet",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"radlwin" ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      do n=1,nrad2d
        if( n.lt.10 )then
          text1 = 'radX    '
          write(text1(4:4),181) n
181       format(i1.1)
        elseif( n.lt.100 )then
          text1 = 'radXX   '
          write(text1(4:5),182) n
182       format(i2.2)
        elseif( n.lt.1000 )then
          text1 = 'radXXX  '
          write(text1(4:6),183) n
183       format(i3.3)
        else
          stop 11611
        endif
        call disp_err( nf90_def_var(ncid,text1     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      enddo
     endif

    if( radopt.ge.1 .and. ptype.eq.5 )then
      call disp_err( nf90_def_var(ncid,"effc"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"effi"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"effs"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"effr"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"effg"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"effis"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
    endif

    if((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1))then
      if(sfcmodel.ge.1)then
        call disp_err( nf90_def_var(ncid,"ust"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","friction velocity") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"znt"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","roughness length") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"cd"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface drag coefficient") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

        call disp_err( nf90_def_var(ncid,"ch"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface exchange coefficient for sensible heat") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

        call disp_err( nf90_def_var(ncid,"cq"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface exchange coefficient for moisture") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )

        call disp_err( nf90_def_var(ncid,"u1"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","u component of velocity at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"v1"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","v component of velocity at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"s1"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","horizontal windspeed at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"t1"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature at lowest model level") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

        call disp_err( nf90_def_var(ncid,"u10"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","u component of windspeed at z=10 m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"v10"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","v component of windspeed at z=10 m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"s10"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","horizontal windspeed at z=10 m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"xland"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","land/water flag (1=land,2=water)") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","integer flag") , .true. )

        call disp_err( nf90_def_var(ncid,"thflux"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface potential temperature flux") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K m s^{-1}") , .true. )

        call disp_err( nf90_def_var(ncid,"qvflux"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface water vapor flux") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","kg kg^{-1} m s^{-1}") , .true. )

        call disp_err( nf90_def_var(ncid,"psfc"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","surface pressure") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )
      if( tbc.eq.3 )then

        call disp_err( nf90_def_var(ncid,"ustt"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","ustt") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

        call disp_err( nf90_def_var(ncid,"ut"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","ut") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

        call disp_err( nf90_def_var(ncid,"vt"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","vt") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )

        call disp_err( nf90_def_var(ncid,"st"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","st") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )
      endif
      endif
      if(sfcmodel.ge.1)then
        call disp_err( nf90_def_var(ncid,"lu_index",nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"kpbl2d"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hfx"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"qfx"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hpbl"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"wspd"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"psim"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"psih"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"gz1oz0"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"br"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"chs"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"chs2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cqs2"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cpmm"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"zol"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"mavail"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"mol"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"rmol"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"regime"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"lh"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"tmn"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"flhc"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"flqc"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"qgh"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"ck"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cka"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"cda"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"ustm"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"qsfc"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )

        call disp_err( nf90_def_var(ncid,"t2"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","temperature at z=2m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

        call disp_err( nf90_def_var(ncid,"q2"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","water vapor mixing ratio at z=2m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )

        call disp_err( nf90_def_var(ncid,"th2"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature at z=2m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )

        call disp_err( nf90_def_var(ncid,"emiss"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"thc"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"albd"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"gsw"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"glw"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"chklowq" ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"capg"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"snowc"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"fm"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"fh"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"mznt"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"swspd"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"wstar"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"delta"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        do n=1,num_soil_layers
          if( n.lt.10 )then
            text1 = 'tslbX   '
            write(text1(5:5),171) n
171         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'tslbXX  '
            write(text1(5:6),172) n
172         format(i2.2)
          else
            stop 22122
          endif
          call disp_err( nf90_def_var(ncid,text1     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        enddo
      endif
      if(oceanmodel.eq.2)then
        call disp_err( nf90_def_var(ncid,"tml"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"t0ml"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hml"     ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"h0ml"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"huml"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"hvml"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
        call disp_err( nf90_def_var(ncid,"tmoml"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      endif
    endif

  !------------------
  !  mynn:

    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. sfcmodel.eq.6 )THEN
      call disp_err( nf90_def_var(ncid,"qke"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qke3d"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"tsq"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qsq"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"cov"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"sh3d"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"sm3d"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"el_pbl"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qc_bl"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qi_bl"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"cldfra_b",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"edmf_a"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"edmf_w"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"edmf_qt" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"edmf_thl",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"edmf_ent",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"edmf_qc" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"sub_thl3",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"sub_sqv3",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"det_thl3",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"det_sqv3",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"thpten"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qvpten"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qcpten"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qipten"  ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"upten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vpten"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qnipten" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qncpten" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vdfg"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"maxmf"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"nupdraft",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ktop_plu",nf90_float,(/niid,njid,timeid/),varid) , .true. )
    ENDIF

  !------------------
  ! myj:

    IF(ipbl.eq.6)THEN
      call disp_err( nf90_def_var(ncid,"tke_myj" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"el_myj " ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
    ENDIF         
    IF(sfcmodel.eq.7)THEN
      call disp_err( nf90_def_var(ncid,"mixht" ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"akhs"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"akms"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"elflx" ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ct"    ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"snow"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"sice"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"thz0"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qz0"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"uz0"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vz0"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"th10"  ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"q10"   ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"z0base",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"zntmyj",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"lowlyr",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ivgtyp",nf90_float,(/niid,njid,timeid/),varid) , .true. )
    ENDIF

  !------------------

    IF( cuparam.ge.1 )THEN
      call disp_err( nf90_def_var(ncid,"raincv",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      if(nrain.eq.2)  &
      call disp_err( nf90_def_var(ncid,"raincv2",nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"thcuten",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qvcuten",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qccuten",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"qicuten",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ucuten" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vcuten" ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
    ENDIF

  !------------------

      call disp_err( nf90_def_var(ncid,"gamk",nf90_float,(/nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"kmw" ,nf90_float,(/nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"ufw" ,nf90_float,(/nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"vfw" ,nf90_float,(/nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"u1b" ,nf90_float,(/nkid,timeid/),varid) , .true. )
      call disp_err( nf90_def_var(ncid,"v1b" ,nf90_float,(/nkid,timeid/),varid) , .true. )

  !------------------

    if( iptra.eq.1 )then
      do n=1,npt
        if( n.le.9 )then
          text1 = 'ptX     '
          write(text1(3:3),161) n
161       format(i1.1)
        elseif( n.le.99 )then
          text1 = 'ptXX    '
          write(text1(3:4),162) n
162       format(i2.2)
        elseif( n.le.999 )then
          text1 = 'ptXXX   '
          write(text1(3:5),163) n
163       format(i3.3)
        else
          stop 11511
        endif
        call disp_err( nf90_def_var(ncid,text1,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      enddo
    endif

    if( iprcl.eq.1 )then
      call disp_err( nf90_def_var(ncid,"ploc"    ,nf90_float,(/nparcelsid,plocid,timeid/),varid) , .true. )
    endif

    if(irbc.eq.4)then
      if( wbc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbcw"  ,nf90_float,(/njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on west boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( ebc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbce"  ,nf90_float,(/njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on east boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( sbc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbcs"  ,nf90_float,(/niid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on south boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( nbc.eq.2 )then
        call disp_err( nf90_def_var(ncid,"radbcn"  ,nf90_float,(/niid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","estimated gravity wave phase speed on north boundary") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
    endif

      !------------------
      ! 150820:  optionals

    !-----
    IF( restart_file_theta )THEN
      call disp_err( nf90_def_var(ncid,"theta"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )
    ENDIF
    IF( restart_file_u0 .or. do_lsnudge .or. do_adapt_move )THEN
      call disp_err( nf90_def_var(ncid,"u0"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
    ENDIF
    IF( restart_file_v0 .or. do_lsnudge .or. do_adapt_move )THEN
      call disp_err( nf90_def_var(ncid,"v0"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
    ENDIF
    IF( restart_file_dbz )THEN
      call disp_err( nf90_def_var(ncid,"dbz"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","reflectivity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","dBZ") , .true. )
    ENDIF
    !-----
    IF( restart_file_th0 )THEN
      call disp_err( nf90_def_var(ncid,"th0"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state potential temperature") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","K") , .true. )
    ENDIF
    IF( restart_file_prs0 )THEN
      call disp_err( nf90_def_var(ncid,"prs0"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","Pa") , .true. )
    ENDIF
    IF( restart_file_pi0 )THEN
      call disp_err( nf90_def_var(ncid,"pi0"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state nondimensional pressure") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","nondimensional") , .true. )
    ENDIF
    IF( restart_file_rho0 )THEN
      call disp_err( nf90_def_var(ncid,"rho0"    ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state dry-air density") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/m^3") , .true. )
    ENDIF
    IF( restart_file_qv0 )THEN
      call disp_err( nf90_def_var(ncid,"qv0"     ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state water vapor mixing ratio") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","kg/kg") , .true. )
    ENDIF
    IF( restart_file_u0 )THEN
      call disp_err( nf90_def_var(ncid,"u0"      ,nf90_float,(/nip1id,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state x-component velocity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
    ENDIF
    IF( restart_file_v0 )THEN
      call disp_err( nf90_def_var(ncid,"v0"      ,nf90_float,(/niid,njp1id,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","base-state y-component velocity") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
    ENDIF
    !-----
    IF( restart_file_zs )THEN
      call disp_err( nf90_def_var(ncid,"zs"      ,nf90_float,(/niid,njid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","terrain height") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
    ENDIF
    IF( restart_file_zh )THEN
      call disp_err( nf90_def_var(ncid,"zhalf"   ,nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","height of half (scalar) grid points (3d array)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
    ENDIF
    IF( restart_file_zf )THEN
      call disp_err( nf90_def_var(ncid,"zfull"   ,nf90_float,(/niid,njid,nkp1id,timeid/),varid) , .true. )
      call disp_err( nf90_put_att(ncid,varid,"long_name","height of full (w) grid points (3d array)") , .true. )
      call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
    ENDIF
    !-----
    IF( restart_file_diags )THEN
      if( td_diss.gt.0 )then
        call disp_err( nf90_def_var(ncid,"dissheat",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","dissipative heating (potential temperature tendency)") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K/s") , .true. )
      endif
      if( td_mp.gt.0 )then
        call disp_err( nf90_def_var(ncid,"mptend",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","potential temperature tendency from microphysics") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","K/s") , .true. )
      endif
      if( qd_vtc.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vtc",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qc") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( qd_vtr.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vtr",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qr") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( qd_vts.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vts",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qs") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( qd_vtg.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vtg",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qg") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
      if( qd_vti.gt.0 )then
        call disp_err( nf90_def_var(ncid,"vti",nf90_float,(/niid,njid,nkid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","terminal fall velocity of qi") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )
      endif
    ENDIF
    !-----

      !------------------
      ! global attributes:

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nx",nx) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ny",ny) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nz",nz) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"imoist",imoist) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"sgsmodel",sgsmodel) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"tconfig",tconfig) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"bcturbs",bcturbs) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ptype",ptype) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"wbc",wbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ebc",ebc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"sbc",sbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nbc",nbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"bbc",bbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"tbc",tbc) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"iorigin",iorigin) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"axisymm",axisymm) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"iptra",iptra) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"npt",npt) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"fcor",fcor) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"radopt",radopt) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dtrad" ,dtrad) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ctrlat",ctrlat) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ctrlon",ctrlon) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"year"  ,year) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"month" ,month) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"day"   ,day) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"hour"  ,hour) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"minute",minute) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"second",second) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"sfcmodel",sfcmodel) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"oceanmodel",oceanmodel) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"ipbl",ipbl) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"iice",iice) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"idm",idm) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"idmplus",idmplus) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"numq",numq) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nql1",nql1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nql2",nql2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nqs1",nqs1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nqs2",nqs2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nnc1",nnc1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nnc2",nnc2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nzl1",nzl1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nzl2",nzl2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nvl1",nvl1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"nvl2",nvl2) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_m",c_m) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_e1",c_e1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_e2",c_e2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"c_s",c_s) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgs1",cgs1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgs2",cgs2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgs3",cgs3) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgs1",dgs1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgs2",dgs2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgs3",dgs3) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgt1",cgt1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgt2",cgt2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"cgt3",cgt3) , .true. )

      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgt1",dgt1) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgt2",dgt2) , .true. )
      call disp_err( nf90_put_att(ncid,NF90_GLOBAL,"dgt3",dgt3) , .true. )

      !------------------
      ! end definitions:

      call disp_err( nf90_enddef(ncid) , .true. )

      !---------------------------------------------
      !ccccccccccccccccccccccccccccccccccccccccccccc
      !---------------------------------------------


      call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
      do i=1,nx
        dumx(i) = 0.5*(xfref(i)+xfref(i+1))
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/nx/)) , .true. )


      call disp_err( nf90_inq_varid(ncid,'xf',varid) , .true. )
      do i=1,nx+1
        dumx(i) = xfref(i)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/nx+1/)) , .true. )


      call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
      do j=1,ny
        dumy(j) = 0.5*(yfref(j)+yfref(j+1))
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/ny/)) , .true. )


      call disp_err( nf90_inq_varid(ncid,'yf',varid) , .true. )
      do j=1,ny+1
        dumy(j) = yfref(j)
      enddo
      call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/ny+1/)) , .true. )


      call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
      if(terrain_flag)then
        do k=1,kval
          dumz(k) = sigma(k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kval/)) , .true. )
      else
        do k=1,kval
          dumz(k) = zh(1,1,k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kval/)) , .true. )
      endif


      call disp_err( nf90_inq_varid(ncid,'zf',varid) , .true. )
      if(terrain_flag)then
        do k=1,kvalp1
          dumz(k) = sigmaf(k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/)) , .true. )
      else
        do k=1,kvalp1
          dumz(k) = zf(1,1,k)
        enddo
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/kvalp1/)) , .true. )
      endif


    ENDIF  ifallinfo

!  END
!-----------------------------------------------------------------------

      call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

      print *,'  leaving restart_prelim '

    end subroutine restart_prelim


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writestat_nc(nrec,rtime,nstat,rstat,qname,budname,name_stat,desc_stat,unit_stat)
      use input
      use constants
      use netcdf
      implicit none

      integer, intent(inout) :: nrec
      real,    intent(in)    :: rtime
      integer, intent(in)    :: nstat
      real, dimension(stat_out), intent(in) :: rstat
      character(len=3), dimension(maxq), intent(in) :: qname
      character(len=6), dimension(maxq), intent(in) :: budname
      character(len=40), intent(in), dimension(maxvars) :: name_stat,desc_stat,unit_stat

      integer :: i,n,ncid,status,dimid,varid,time_index,timeid,tfile
      integer :: xhid,yhid,zhid
      logical :: allinfo

      if( myid.eq.0 ) print *,'  Entering writestat_nc '

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string = 'cm1out_stats.nc'

    if(myid.eq.0) print *,string

    allinfo = .false.
    IF(nrec.eq.1) allinfo=.true.

    ncid = 1

    IF( nrec.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set nrec to 1 and write all info.
      status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.nrec )then
          if(myid.eq.0) print *,'  tfile,nrec = ',tfile,nrec
          nrec = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        nrec = 1
      endif
    ENDIF

    if( myid.eq.0 ) print *,'  nrec = ',nrec


  allinfo2:  IF( allinfo )THEN
    ! Definitions/descriptions:

      if( myid.eq.0 ) print *,'  calling nf90_create '
#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

    call disp_err( nf90_def_dim(ncid,"xh",1,xhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"zh",1,zhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"time",nf90_unlimited,timeid) , .true. )

    call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

    call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

    call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","height") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

    call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

  !---------------------------------------------------------------------------!

    ! new for cm1r19:  use "_stat" arrays, which are
    !                  defined in statpack.F

      statloop:  &
      do n = 1,stat_out
        call disp_err( nf90_def_var(ncid,trim(name_stat(n)),nf90_float,timeid,varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name",trim(desc_stat(n))) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units",trim(unit_stat(n))) , .true. )
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/1/)) , .true. )
#endif
      enddo  statloop

  !---------------------------------------------------------------------------!

    call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

    call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
    call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

    call disp_err( nf90_enddef(ncid=ncid) , .true. )

    call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
    call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

    call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
    call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

    call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
    call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

  ENDIF  allinfo2

  !---------------------------------------------------------------------
    ! Write data:

    time_index = nrec

    call disp_err( nf90_inq_varid(ncid,'time',timeid) , .true. )
    call disp_err( nf90_put_var(ncid,timeid,rtime,(/time_index/)) , .true. )

    DO n=1,nstat
      call disp_err( nf90_inq_varid(ncid,trim(name_stat(n)),varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,rstat(n),(/time_index/)) , .true. )
    ENDDO

    ! close file

    call disp_err( nf90_close(ncid) , .true. )

    ! all done

      if( myid.eq.0 ) print *,'  Leaving writestat_nc '

      end subroutine writestat_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine disp_err( status , stop_on_error )
      use netcdf
      implicit none

      integer, intent(in) :: status
      logical, intent(in) :: stop_on_error

      IF( status.ne.nf90_noerr )THEN
        IF( stop_on_error )THEN
          print *,'  netcdf status returned an error: ', status,' ... stopping program'
          print *
          print *,nf90_strerror(status)
          print *
          call stopcm1
        ENDIF
      ENDIF

      end subroutine disp_err

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writepdata_nc(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,pdata2)
      use input
      use constants
      use netcdf
      implicit none

      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels) :: pdata2

      integer :: i,ncid,status,dimid,varid,time_index,n,n2,np,timeid,tfile,xid
      integer :: npid,yhid,zhid
      logical :: allinfo

!-----------------------------------------------------------------------

      if( myid.eq.0 ) print *,'  Entering writepdata_nc '

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string(1:22) = 'cm1out_pdata.nc'

    allinfo = .false.
    IF(prec.eq.1) allinfo=.true.

    ncid = 1

    IF( prec.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set prec to 1 and write all info.
      status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.prec )then
!!!          if(myid.eq.0) print *,'  tfile,prec = ',tfile,prec
          prec = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        prec = 1
      endif
    ENDIF

    if( myid.eq.0 ) print *,'  pdata prec = ',prec


  allinfo3:  IF( allinfo )THEN
    ! Definitions/descriptions:

#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

    call disp_err( nf90_def_dim(ncid,"xh",nparcels,npid) , .true. )
    call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"zh",1,zhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"time",nf90_unlimited,timeid) , .true. )

    call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/npid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location ... actually, really parcel ID number") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )

    call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )

    call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","height") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

    call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

  !---------------------------------------------------------------------------!

    ! new for cm1r19:  use "_prcl" arrays, which are
    !                  defined in parcel.F

      prclloop:  &
      do n = 1,prcl_out
        if(myid.eq.0) print *,n,trim(name_prcl(n))
        call disp_err( nf90_def_var(ncid,trim(name_prcl(n)),nf90_float,(/npid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name",trim(desc_prcl(n))) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units",trim(unit_prcl(n))) , .true. )
!!!#ifdef NCFPLUS
!!!        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/nparcels,1/)) , .true. )
!!!#endif
      enddo  prclloop

  !---------------------------------------------------------------------------!

    call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

    call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
    call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

    call disp_err( nf90_enddef(ncid) , .true. )

  do np=1,nparcels
    call disp_err( nf90_put_var(ncid,npid,float(np),(/np/)) , .true. )
  enddo
    call disp_err( nf90_put_var(ncid,yhid,0.0) , .true. )
    call disp_err( nf90_put_var(ncid,zhid,0.0) , .true. )

!------------------------

  ENDIF  allinfo3

      ! Write data:

      time_index = prec

      call disp_err( nf90_inq_varid(ncid,'time',timeid) , .true. )
      call disp_err( nf90_put_var(ncid,timeid,rtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,'x',xid) , .true. )

      DO n=1,prcl_out
        call disp_err( nf90_inq_varid(ncid,trim(name_prcl(n)),varid) , .true. )
        do np=1,nparcels
          pdata2(np) = pdata(np,n)
        enddo
        call disp_err( nf90_put_var(ncid,varid,pdata2,(/1,time_index/),(/nparcels,1/)) , .true. )
      ENDDO

      ! close file

      call disp_err( nf90_close(ncid) , .true. )

      prec = prec + 1

      ! all done

      if( myid.eq.0 ) print *,'  Leaving writepdata_nc '

      end subroutine writepdata_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writehifrq_nc(wloop,nlevels,nvar,varmax,varname,vardesc,varunit,varlvls,ncid,xhref,yhref,sigma,klev,zlev,mtime,dumx,dumy,dumz)
      use input
      use netcdf
      use constants , only : grads_undef
      implicit none

      integer, intent(in) :: wloop,nlevels,nvar,varmax
      character(len=80), dimension(varmax), intent(in) :: varname,vardesc,varunit
      integer, dimension(varmax), intent(in) :: varlvls
      integer, intent(inout) :: ncid
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xhref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yhref
      real, intent(in), dimension(kb:ke) :: sigma
      integer, intent(in), dimension(nk) :: klev
      real, intent(in), dimension(nk) :: zlev
      double precision, intent(in) :: mtime
      real, intent(inout), dimension(nx) :: dumx
      real, intent(inout), dimension(ny) :: dumy
      real, intent(inout), dimension(nlevels) :: dumz

      integer :: i,j,k,n,xhid,yhid,zhid,zfid,timeid,status,dimid,varid,time_index
      real :: actual_time

      if( wloop.eq.1 )then
        ! open file, define vars:

        ncid = 1

#ifdef NCFPLUS
        call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
        call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

        call disp_err( nf90_def_dim(ncid,"xh",nx,xhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"yh",ny,yhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zh",nlevels,zhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"time",1,timeid) , .true. )

        call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

        call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

        call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of scalar levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

        do n = 1,nvar
          if( varlvls(n).eq.0 )then
            print *,'    n,varname,0 = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,timeid/),varid) , .true. )
          else
            print *,'    n,varname,s = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
          endif
          call disp_err( nf90_put_att(ncid,varid,"long_name",trim(vardesc(n))) , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units",trim(varunit(n))) , .true. )
        enddo

        call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

        call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
        call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

        call disp_err( nf90_enddef(ncid) , .true. )

        do i=1,nx
          dumx(i) = xhref(i)
        enddo
        call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/nx/)) , .true. )

        do j=1,ny
          dumy(j) = yhref(j)
        enddo
        call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dumy,(/1/),(/ny/)) , .true. )

        do k=1,nlevels
          dumz(k) = sigma(klev(k))
        enddo
        call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dumz,(/1/),(/nlevels/)) , .true. )

        time_index = 1
        actual_time = mtime

        call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

      elseif( wloop.eq.2 )then
        ! close file:

        call disp_err( nf90_close(ncid) , .true. )

      else
        stop 32123
      endif

      end subroutine writehifrq_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writediag_nc(wloop,varmax,ncid,varname,vardesc,varunit,varlvls,vargrid,nvar,mtime,zh,zf,dum1,dum2)
      use input
      use netcdf
      use constants , only : grads_undef
      implicit none

      integer, intent(in) :: wloop,varmax
      integer, intent(inout) :: ncid
      character(len=80), dimension(varmax), intent(in) :: varname,vardesc,varunit
      integer, dimension(varmax), intent(in) :: varlvls
      character(len=1),  dimension(varmax), intent(in) :: vargrid
      integer, intent(in) :: nvar
      double precision, intent(in) :: mtime
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(nk) :: dum1
      real, intent(inout), dimension(nk+1) :: dum2

      integer :: i,j,k,n,xhid,yhid,zhid,zfid,timeid,status,dimid,varid,time_index
      real :: actual_time

      if( wloop.eq.1 )then
        ! open file:

#ifdef NCFPLUS
        call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
        call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

        call disp_err( nf90_def_dim(ncid,"xh",1,xhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zh",nk,zhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zf",nk+1,zfid) , .true. )
        call disp_err( nf90_def_dim(ncid,"time",1,timeid) , .true. )

        call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

        call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

        call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of scalar levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"zf",nf90_float,(/zfid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of w levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

        do n = 1,nvar
          if( vargrid(n).eq.'2' )then
            print *,'    n,varname,2 = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,timeid/),varid) , .true. )
          elseif( vargrid(n).eq.'s' )then
            print *,'    n,varname,s = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
          elseif( vargrid(n).eq.'w' )then
            print *,'    n,varname,w = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,zfid,timeid/),varid) , .true. )
          else
            stop 23987
          endif
          call disp_err( nf90_put_att(ncid,varid,"long_name",trim(vardesc(n))) , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units",trim(varunit(n))) , .true. )
        enddo

        call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

        call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
        call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

        call disp_err( nf90_enddef(ncid) , .true. )

        call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        do k=1,nk
          dum1(k) = zh(1,1,k)
        enddo
        call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dum1,(/1/),(/nk/)) , .true. )

        do k=1,nk+1
          dum2(k) = zf(1,1,k)
        enddo
        call disp_err( nf90_inq_varid(ncid,'zf',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dum2,(/1/),(/nk+1/)) , .true. )

        time_index = 1
        actual_time = mtime

        call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

      elseif( wloop.eq.2 )then
        ! close file:

        call disp_err( nf90_close(ncid) , .true. )

      else
        stop 32123
      endif


      end subroutine writediag_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writeazim_nc(wloop,icrs,nvar,varmax,varname,vardesc,varunit,varlvls,vargrid,ncid,mtime,zh,zf,dumx,dum1,dum2)
      use input
      use netcdf
      use constants , only : grads_undef
      implicit none

      integer, intent(in) :: wloop,icrs,nvar,varmax
      character(len=80), dimension(varmax), intent(in) :: varname,vardesc,varunit
      integer, dimension(varmax), intent(in) :: varlvls
      character(len=1),  dimension(varmax), intent(in) :: vargrid
      integer, intent(inout) :: ncid
      double precision, intent(in) :: mtime
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(icrs) :: dumx
      real, intent(inout), dimension(nk)   :: dum1
      real, intent(inout), dimension(nk+1) :: dum2

      integer :: i,j,k,n,xhid,yhid,zhid,zfid,timeid,status,dimid,varid,time_index
      real :: actual_time

      if( wloop.eq.1 )then
        ! open file:

        ncid = 1

#ifdef NCFPLUS
        call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
        call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

        call disp_err( nf90_def_dim(ncid,"xh",icrs,xhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zh",nk,zhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zf",nk+1,zfid) , .true. )
        call disp_err( nf90_def_dim(ncid,"time",1,timeid) , .true. )

        call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","radius") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units",trim(aunit)) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

        call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

        call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of scalar levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"zf",nf90_float,(/zfid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of w levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

        do n = 1,nvar
          if( varlvls(n).eq.0 )then
            print *,'    n,varname,0 = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,timeid/),varid) , .true. )
          elseif( vargrid(n).eq.'s' )then
            print *,'    n,varname,s = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
          elseif( vargrid(n).eq.'w' )then
            print *,'    n,varname,w = ',n,trim(varname(n))
            call disp_err( nf90_def_var(ncid,trim(varname(n)),nf90_float,(/xhid,yhid,zfid,timeid/),varid) , .true. )
          else
            stop 23987
          endif
          call disp_err( nf90_put_att(ncid,varid,"long_name",trim(vardesc(n))) , .true. )
          call disp_err( nf90_put_att(ncid,varid,"units",trim(varunit(n))) , .true. )
        enddo

        call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

        call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
        call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

        call disp_err( nf90_enddef(ncid) , .true. )

        do i=1,icrs
          dumx(i) = 0.5*ddr + (i-1)*ddr
        enddo
        call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dumx,(/1/),(/icrs/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        do k=1,nk
          dum1(k) = zh(1,1,k)
        enddo
        call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dum1,(/1/),(/nk/)) , .true. )

        do k=1,nk+1
          dum2(k) = zf(1,1,k)
        enddo
        call disp_err( nf90_inq_varid(ncid,'zf',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,dum2,(/1/),(/nk+1/)) , .true. )

        time_index = 1
        actual_time = mtime

        call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

      elseif( wloop.eq.2 )then
        ! close file:

        call disp_err( nf90_close(ncid) , .true. )

      else
        stop 32123
      endif



      end subroutine writeazim_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writecenter_nc(nwritea,ncid,mtime,icenter,jcenter,xcenter,ycenter,vmax,rmax,zmax)
      use input
      use netcdf
      use constants , only : grads_undef
      implicit none

      integer, intent(in) :: nwritea
      integer, intent(inout) :: ncid
      double precision, intent(in) :: mtime
      integer, intent(in) :: icenter,jcenter
      real,    intent(in) :: xcenter,ycenter,vmax,rmax,zmax

      integer :: i,j,k,n,xhid,yhid,zhid,zfid,timeid,status,dimid,varid,time_index
      integer :: tfile,num_write
      logical :: allinfo
      real :: actual_time
      integer, dimension(1,1,1,1) :: ivar
      real,    dimension(1,1,1,1) :: rvar

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_azimavg_stats.nc'
        print *,'  string = ',string

        allinfo = .false.
        IF(nwritea.eq.1) allinfo=.true.
        num_write = nwritea

        ncid = 1

        IF( nwritea.ne.1 )THEN
          !   Try to open file.
          !   If error, set num_write to 1 and write all info.
          status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
          if( status.eq.nf90_noerr )then
            ! no error, file exists.  Get number of time levels in file:
            call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
            call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
            if( (tfile+1).lt.num_write )then
              if(myid.eq.0) print *,'  tfile,num_write = ',tfile,num_write
              num_write = tfile+1
            endif
          else
            ! if error opening file, then write all info:
            if(myid.eq.0) print *,'  status = ',status
            allinfo = .true.
            num_write = 1
          endif
        ENDIF

        time_index = num_write

        print *,'  allinfo,time_index = ',allinfo,time_index

    ifallinfo:  &
    IF(allinfo)THEN

#ifdef NCFPLUS
        call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
        call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

        call disp_err( nf90_def_dim(ncid,"xh",1,xhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zh",1,zhid) , .true. )
        call disp_err( nf90_def_dim(ncid,'time',nf90_unlimited,timeid) , .true. )

        call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

        call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

        call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of scalar levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

        call disp_err( nf90_def_var(ncid,"icenter",nf90_int,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","i value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","unitless") , .true. )

        call disp_err( nf90_def_var(ncid,"jcenter",nf90_int,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","j value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","unitless") , .true. )

        call disp_err( nf90_def_var(ncid,"xcenter",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","x value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"ycenter",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","y value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"vmax",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","max azimuthally averaged tangential velocity"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"rmax",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","radius of Vmax"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"zmax",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height (ASL) of Vmax"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

        call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
        call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

        call disp_err( nf90_enddef(ncid) , .true. )

        call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

    ENDIF  ifallinfo

        actual_time = mtime
        call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
        ivar(1,1,1,1) = icenter
        call disp_err( nf90_put_var(ncid,varid,ivar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
        ivar(1,1,1,1) = jcenter
        call disp_err( nf90_put_var(ncid,varid,ivar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
        rvar(1,1,1,1) = xcenter
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
        rvar(1,1,1,1) = ycenter
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"vmax",varid) , .true. )
        rvar(1,1,1,1) = vmax
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"rmax",varid) , .true. )
        rvar(1,1,1,1) = rmax
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"zmax",varid) , .true. )
        rvar(1,1,1,1) = zmax
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_close(ncid) , .true. )

      end subroutine writecenter_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writemove_nc(nwritemv,ncid,mtime,icenter,jcenter,xcenter,ycenter,umovex,vmovex)
      use input
      use netcdf
      use constants , only : grads_undef
      implicit none

      integer, intent(in) :: nwritemv
      integer, intent(inout) :: ncid
      double precision, intent(in) :: mtime
      integer, intent(in) :: icenter,jcenter
      real,    intent(in) :: xcenter,ycenter,umovex,vmovex

      integer :: i,j,k,n,xhid,yhid,zhid,zfid,timeid,status,dimid,varid,time_index
      integer :: tfile,num_write
      logical :: allinfo
      real :: actual_time
      integer, dimension(1,1,1,1) :: ivar
      real,    dimension(1,1,1,1) :: rvar

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_move.nc'
        print *,'  string = ',string

        allinfo = .false.
        IF(nwritemv.eq.1) allinfo=.true.
        num_write = nwritemv

        ncid = 1

        IF( nwritemv.ne.1 )THEN
          !   Try to open file.
          !   If error, set num_write to 1 and write all info.
          status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
          if( status.eq.nf90_noerr )then
            ! no error, file exists.  Get number of time levels in file:
            call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
            call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
            if( (tfile+1).lt.num_write )then
              if(myid.eq.0) print *,'  tfile,num_write = ',tfile,num_write
              num_write = tfile+1
            endif
          else
            ! if error opening file, then write all info:
            if(myid.eq.0) print *,'  status = ',status
            allinfo = .true.
            num_write = 1
          endif
        ENDIF

        time_index = num_write

        print *,'  allinfo,time_index = ',allinfo,time_index

    ifallinfo:  &
    IF(allinfo)THEN

#ifdef NCFPLUS
        call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
        call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

        call disp_err( nf90_def_dim(ncid,"xh",1,xhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
        call disp_err( nf90_def_dim(ncid,"zh",1,zhid) , .true. )
        call disp_err( nf90_def_dim(ncid,'time',nf90_unlimited,timeid) , .true. )

        call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/xhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","X") , .true. )

        call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Y") , .true. )

        call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","height of scalar levels") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","Z") , .true. )

        call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
        call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

        call disp_err( nf90_def_var(ncid,"icenter",nf90_int,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","i value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","unitless") , .true. )

        call disp_err( nf90_def_var(ncid,"jcenter",nf90_int,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","j value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","unitless") , .true. )

        call disp_err( nf90_def_var(ncid,"xcenter",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","x value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"ycenter",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","y value of center grid point"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

        call disp_err( nf90_def_var(ncid,"umove",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","umove"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_def_var(ncid,"vmove",nf90_float,(/xhid,yhid,zhid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name","vmove"), .true. )
        call disp_err( nf90_put_att(ncid,varid,"units","m/s") , .true. )

        call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

        call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
        call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

        call disp_err( nf90_enddef(ncid) , .true. )

        call disp_err( nf90_inq_varid(ncid,'xh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        call disp_err( nf90_inq_varid(ncid,'yh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

        call disp_err( nf90_inq_varid(ncid,'zh',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,0.0) , .true. )

    ENDIF  ifallinfo

        actual_time = mtime
        call disp_err( nf90_inq_varid(ncid,'time',varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,actual_time,(/time_index/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
        ivar(1,1,1,1) = icenter
        call disp_err( nf90_put_var(ncid,varid,ivar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
        ivar(1,1,1,1) = jcenter
        call disp_err( nf90_put_var(ncid,varid,ivar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
        rvar(1,1,1,1) = xcenter
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
        rvar(1,1,1,1) = ycenter
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
        rvar(1,1,1,1) = umovex
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
        rvar(1,1,1,1) = vmovex
        call disp_err( nf90_put_var(ncid,varid,rvar,(/1,1,1,time_index/),(/1,1,1,1/)) , .true. )

        call disp_err( nf90_close(ncid) , .true. )

      end subroutine writemove_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif

  END MODULE writeout_nc_module
