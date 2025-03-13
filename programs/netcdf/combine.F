!------------------------------------------------------------------------
! This program combines netCDF data from different processors from
! George Bryan's cm1 numerical model.
!
! Provided by:  Daniel Kirshbaum, University of Reading
! 
! Last modified:  2 April 2008
!------------------------------------------------------------------------

      program combine
      implicit none

      include '/usr/local/netcdf/include/netcdf.inc'  ! My machine
!      include '/usr/local/netcdf-pgi/include/netcdf.inc' ! MMM at NCAR

!-----------------------------------------
! input and auxiliary variables
!-----------------------------------------
      integer :: nfile,xtype,nwrite
      character(len=120), dimension(200) :: fname
      integer, dimension(200) :: ncid
      real,dimension(200) :: x_min,x_max,y_min,y_max
      integer, dimension(200) :: nx,nxp1,ny,nyp1

      character(len=120) :: fname1
      character(len=15) :: fname2
      character(len=2) :: fid
      integer :: ncid1,ncid2

      integer :: i, j, k, n,lenf, nVars2

  ! global vars for the original file
      integer :: status, nDims, nVars, nGlobalAtts

  ! dimensions vars
      integer :: DimID
      integer, dimension(nf_max_var_dims) :: DimLen
      character(nf_max_name) :: DimName

  ! attributes vars
      character(nf_max_name) :: AttName
  ! variables vars
      character(nf_max_name) :: VarName
      integer :: VarID, VarID1, nDimsVar, nAttsVar, VarType
      integer, dimension(4) :: VarDim
      integer, dimension(nf_max_var_dims) :: VarDimIDs

      character(nf_max_name) :: TimeDimName
      integer :: TimeDimLength, TimeDimID

      integer :: mx,mxp1,my,myp1,nz,nzp1,nt
      integer :: istart,istop,jstart,jstop
      real :: x_min1,x_max1,y_min1,y_max1,z_min,z_max
      real :: dx,dy,dz,ztop

      real :: var1D
      real, dimension(:), allocatable :: time,sgz,wgz
      real, dimension(:,:), allocatable :: zs,var2D
      real, dimension(:,:,:), allocatable :: var3D
!-----------------------------------------------------------------------

  ! read in the netCDF file name

      write(*,91)
 91   format(/,'Input number of netCDF files to be processed :')
      read(*,*) nfile

      write(*,92)
 92   format(/,'Input output time index (nwrite) :')
      read(*,*) nwrite

      print *, 'nwrite=',nwrite

      do i=1,nfile
        fname2 = 'out00XX_001.cdf'
        write(fname2(6:7),93) i-1
        write(fname2(9:11),94) nwrite
 93     format(i2.2)
 94     format(i3.3)
        fname(i) = fname2
        write(*,95) fname(i)
 95     format(/,'Reading ',A)
  ! open the netCDF file to read
        status=nf_open(fname(i),nf_nowrite,ncid(i))
        if (status /= nf_noerr) call disp_err(status)

 ! get dimensions here
        status = nf_inq_dimid(ncid(i),"ni", DimID)
        status = nf_inq_dimlen(ncid(i), DimID,nx(i))
        status = nf_inq_dimid(ncid(i),"nip1", DimID)
        status = nf_inq_dimlen(ncid(i), DimID,nxp1(i))
        status = nf_inq_dimid(ncid(i),"nj", DimID)
        status = nf_inq_dimlen(ncid(i), DimID,ny(i))
        status = nf_inq_dimid(ncid(i),"njp1", DimID)
        status = nf_inq_dimlen(ncid(i), DimID,nyp1(i))
        status = nf_inq_dimid(ncid(i),"nk", DimID)
        status = nf_inq_dimlen(ncid(i), DimID,nz)
        status = nf_inq_dimid(ncid(i),"nkp1", DimID)
        status = nf_inq_dimlen(ncid(i), DimID,nzp1)
        status = nf_inq_dimid(ncid(i),"time",DimID)
        status = nf_inq_dimlen(ncid(i), DimID,nt)
        write(*,96) nx(i),nxp1(i),ny(i),nyp1(i),nz,nzp1,nt
 96     format('nx = ',I5,5x,'nxp1 = ',I5,/,'ny = ',I5,5x,'nyp1 = ',
     *       I5,/,'nz = ',I5,5x,'nzp1 = ',I5,/,'nt = ',I5)
      enddo
!-----------------------------------------------------------------------
! get dx,dy,dz
      status = nf_get_att_real(ncid(1),nf_global,"x_delta",dx)  
      if (status /= nf_noerr) call disp_err(status)

      status = nf_get_att_real(ncid(1),nf_global,"y_delta",dy)
      if (status /= nf_noerr) call disp_err(status)

      write(*,97) dx,dy,dz
 97   format('dx =',f15.5,5x,'dy = ',f15.5,5x,'dz = ',f15.5)

! get x_min, x_max, y_min, y_max, z_min, z_max
      do i=1,nfile
        status = nf_get_att_real(ncid(i),nf_global,"x_min",x_min(i))
        status = nf_get_att_real(ncid(i),nf_global,"x_max",x_max(i))
        status = nf_get_att_real(ncid(i),nf_global,"y_min",y_min(i))
        status = nf_get_att_real(ncid(i),nf_global,"y_max",y_max(i))
        status = nf_get_att_real(ncid(i),nf_global,"z_min",z_min)
        status = nf_get_att_real(ncid(i),nf_global,"z_max",z_max)
      enddo

      x_min1=minval(x_min)
      x_max1=maxval(x_max)
      y_min1=minval(y_min)
      y_max1=maxval(y_max)

      write(*,98) x_min1,x_max1
 98   format('x_min =',f15.5,5x,'x_max = ',f15.5)
      write(*,98) y_min1,y_max1
 99   format('y_min =',f15.5,5x,'y_max = ',f15.5)
      write(*,10) z_min,z_max
 10   format('z_min =',f15.5,5x,'z_max = ',f15.5)

! total dimensions
      mx=(x_max1-x_min1)/dx
      mxp1=mx+1
      my=(y_max1-y_min1)/dx
      myp1=my+1
          
! times
      allocate(time(nt))
      time(:)=0.0
      status = nf_inq_varid(ncid(1), 'time', varid)
      status = nf_get_var_real(ncid(1), varid,time(1:nt))

! sgz
      allocate(sgz(nz))
      sgz(:)=0.0
      status = nf_inq_varid(ncid(1), 'sgz', varid)
      status = nf_get_var_real(ncid(1), varid, sgz(1:nz))

! wgz
      allocate(wgz(nzp1))
      wgz(:)=0.0
      status = nf_inq_varid(ncid(1), 'wgz', varid)
      status = nf_get_var_real(ncid(1), varid, wgz(1:nzp1))
      ztop = wgz(nzp1)

! zs
      allocate(zs(mx,my))
      zs(:,:)=0.0
      do i=1,nfile
        istart=(x_min(i)-x_min1)/dx+1
        istop=istart+nx(i)-1
        jstart=(y_min(i)-y_min1)/dy+1
        jstop=jstart+ny(i)-1
        status = nf_inq_varid(ncid(i),'zbot_p',varid)
        status = nf_get_vara_real(ncid(i),varid,(/1,1/),
     *              (/nx(i),ny(i)/),zs(istart:istop,jstart:jstop))
      enddo

!-----------------------------------------------------------------------
! write output, necessary variables
!-----------------------------------------------------------------------
      write(*,81) 
 81   format(/,'Input name of the new netCDF file (name.cdf):')
      read(*,*) fname1
        
      status = nf_create(fname1,nf_noclobber,ncid1)
      if (status /= nf_noerr) call disp_err(status)

      status = nf_def_dim(ncid1,"nx",mx,DimID)
      status = nf_def_dim(ncid1,"ny",my,DimID)
      status = nf_def_dim(ncid1,"nz",nz,DimID)
      status = nf_def_dim(ncid1,"nxp1",mxp1,DimID)
      status = nf_def_dim(ncid1,"nyp1",myp1,DimID)
      status = nf_def_dim(ncid1,"nzp1",nzp1,DimID)
      status = nf_def_dim(ncid1,"time",nf_UNLIMITED,DimID)
      status = nf_def_dim(ncid1,"one",1,DimID)

      write(*,82)
 82   format(/,'writing global attributes')

      status = nf_put_att_real(ncid1,nf_global,"x_min",
&              nf_real,1,x_min1)
      status = nf_put_att_real(ncid1,nf_global,"x_max",
&              nf_real,1,x_max1)
      status = nf_put_att_real(ncid1,nf_global,"x_delta",
&              nf_real,1,dx)
      status = nf_put_att_text(ncid1,nf_global,"x_units",2,"km")
      status = nf_put_att_text(ncid1,nf_global,"x_label",1,"x")
      status = nf_put_att_text(ncid1,nf_global,"x_display_units",
&              2,"km")
      status = nf_put_att_real(ncid1,nf_global,"y_min",
&              nf_real,1,y_min1)
      status = nf_put_att_real(ncid1,nf_global,"y_max",
&              nf_real,1,y_max1)
      status = nf_put_att_real(ncid1,nf_global,"y_delta",
&              nf_real,1,dy)
      status = nf_put_att_text(ncid1,nf_global,"y_units",2,"km")
      status = nf_put_att_text(ncid1,nf_global,"y_label",1,"y")
      status = nf_put_att_text(ncid1,nf_global,"y_display_units",
&              2,"km")
      status = nf_put_att_real(ncid1,nf_global,"z_min",
&              nf_real,1,z_min)
      status = nf_put_att_real(ncid1,nf_global,"z_max",
&              nf_real,1,z_max)
      status = nf_put_att_text(ncid1,nf_global,"z_units",2,"km")
      status = nf_put_att_text(ncid1,nf_global,"z_label",1,"z")
      status = nf_put_att_text(ncid1,nf_global,"z_display_units",
&              2,"km")
      status = nf_put_att_text(ncid1,nf_global,"runname",1,fname(1))

      status = nf_redef(ncid1)
      status = nf_def_var(ncid1,"f_cor",nf_real,1,8,varid)
      status = nf_put_att_text(ncid1,varid,"units",3,"1/s")
      status = nf_put_att_text(ncid1,varid,"def",
&              18,"coriolis parameter")
      status = nf_put_att_int(ncid1,varid,"no_button",
&              nf_int,1,int(1))
      status = nf_enddef(ncid1)
      var1D=0.0
      status = nf_inq_varid(ncid(1),"f_cor",varid1)
      status = nf_get_var_real(ncid(1),varid1,var1D)
      status = nf_put_var_real(ncid1,varid,var1D)
      status = nf_redef(ncid1)

      status = nf_def_var(ncid1,"ztop",nf_real,1,8,varid)
      status = nf_put_att_text(ncid1,varid,"units",2,"km")
      status = nf_put_att_text(ncid1,varid,"def",12,"domain depth")
      status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
      status = nf_enddef(ncid1)
      status = nf_put_var_real(ncid1,varid,ztop)
      status = nf_redef(ncid1)

      status = nf_def_var(ncid1,"time",nf_float,1,7,varid)
      status = nf_put_att_text(ncid1,varid,"units",1,"s")
      status = nf_enddef(ncid1)
      print *, 'nt,time(nt)',nt,time(nt)
      status = nf_put_vara_real(ncid1,varid,1,nt,time)
      status = nf_redef(ncid1)

      status = nf_def_var(ncid1,"sgz",nf_float,1,3,varid)
      status = nf_put_att_text(ncid1,varid,"units",2,"km")
      status = nf_put_att_text(ncid1,varid,"def",11,"staggered z") 
      status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
      status = nf_enddef(ncid1)
      status = nf_put_vara_real(ncid1,varid,1,nz,sgz)
      status = nf_redef(ncid1)

      status = nf_def_var(ncid1,"wgz",nf_float,1,6,varid)
      status = nf_put_att_text(ncid1,varid,"units",2,"km")
      status = nf_put_att_text(ncid1,varid,"def",13,"unstaggered z") 
      status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
      status = nf_enddef(ncid1)
      status = nf_put_vara_real(ncid1,varid,1,nzp1,wgz)
      status = nf_redef(ncid1)

      do i = 1,nfile

        istart=(x_min(i)-x_min1)/dx+1
        istop=istart+nx(i)-1
        jstart=(y_min(i)-y_min1)/dy+1
        jstop=jstart+ny(i)-1

        status = nf_def_var(ncid1,"zbot_p",nf_float,2,(/1,2/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1,"m")
        status = nf_put_att_text(ncid1,varid,"def",12,"ter at p pts")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                         1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                         1,y_min1+0.5*dy)
        status = nf_enddef(ncid1)
        status = nf_put_var_real(ncid1,varid,zs(1:mx,1:my))
        status = nf_redef(ncid1)

      enddo

      status = nf_inq_varid(ncid(1),'thflux',varid)
      if (status.eq.nf_noerr) then
        print *,'defining thflux'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"thflux",nf_float,2,(/1,2/),varid)
        status = nf_put_att_text(ncid1,varid,"units",7,"K/m/m/s")
        status = nf_put_att_text(ncid1,varid,"def",10,"theta flux")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_enddef(ncid1)
        status = nf_inq_varid(ncid1,"thflux",varid)
        print *,'writing thflux'
        allocate(var2D(mx,my))
        var2D(:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1
          status = nf_inq_varid(ncid(i),"thflux",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var2D(istart:istop,jstart:jstop))
        enddo
        status = nf_put_var_real(ncid1,varid,var2D(1:mx,1:my))
        deallocate(var2D)
      endif

      status = nf_inq_varid(ncid(1),'qvflux',varid)
      if (status.eq.nf_noerr) then
        print *,'defining qvflux'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"qvflux",nf_float,2,(/1,2/),varid)
        status = nf_put_att_text(ncid1,varid,"units",7,"1/m/m/s")
        status = nf_put_att_text(ncid1,varid,"def",7,"qv flux") 
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_enddef(ncid1)
        status = nf_inq_varid(ncid1,"qvflux",varid)
        print *,'writing qvflux'
        allocate(var2D(mx,my))
        var2D(:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1
          status = nf_inq_varid(ncid(i),"qvflux",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var2D(istart:istop,jstart:jstop))
        enddo
        status = nf_put_var_real(ncid1,varid,var2D(1:mx,1:my))
        deallocate(var2D)
      endif

      status = nf_inq_varid(ncid(1),'rain',varid)
      if (status.eq.nf_noerr) then
        print *,'defining rain'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"rain",nf_float,2,(/1,2/),varid)
        status = nf_put_att_text(ncid1,varid,"units",2,"mm")
        status = nf_put_att_text(ncid1,varid,"def",4,"rain") 
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_enddef(ncid1)
        status = nf_inq_varid(ncid1,"rain",varid)
        print *,'writing qvflux'
        allocate(var2D(mx,my))
        var2D(:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1
          status = nf_inq_varid(ncid(i),"rain",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var2D(istart:istop,jstart:jstop))
        enddo
        status = nf_put_var_real(ncid1,varid,var2D(1:mx,1:my))
        deallocate(var2D)
      endif

      print *, 'Done with 1D and 2D'

!-----------------------------------------------------------------------
!  Define 3D variables
!-----------------------------------------------------------------------
      status = nf_inq_varid(ncid(1),'zh',varid)
      if (status.eq.nf_noerr) then
        print *,'defining zh'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"zh",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",8,"km")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing zh'
        status = nf_inq_varid(ncid1,"zh",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"zh",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,
     *              var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      endif

      status = nf_inq_varid(ncid(1),'U0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining u0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"u0",nf_float,3,(/4,2,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,1,x_min1)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
        status = nf_enddef(ncid1)
        print *,'writing u0'
        status = nf_inq_varid(ncid1,"u0",varid)
        allocate(var3D(mxp1,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"u0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop+1,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mxp1,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'V0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining v0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"v0",nf_float,3,(/1,5,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
        status = nf_enddef(ncid1)
        print *,'writing v0'
        status = nf_inq_varid(ncid1,"v0",varid)
        allocate(var3D(mx,myp1,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"v0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop+1,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:myp1,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'pi0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining pi0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"pi0",nf_float,3,(/1,2,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1," ")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                         1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                         1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
        status = nf_enddef(ncid1)
        print *,'writing pi0'
        status = nf_inq_varid(ncid1,"pi0",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"pi0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *            var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'th0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining th0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"th0",nf_float,3,(/1,2,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1,"K")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing th0'
        status = nf_inq_varid(ncid1,"th0",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"th0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,
     *             var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'thet0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining th0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"th0",nf_float,3,(/1,2,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1,"K")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing th0'
        status = nf_inq_varid(ncid1,"th0",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"thet0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,
     *             var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'rho0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining rho0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"rho0",nf_float,3,(/1,2,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",8,"kg/m/m/m")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
        status = nf_enddef(ncid1)
        print *,'writing rho0'
        status = nf_inq_varid(ncid1,"rho0",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"rho0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'qv0',varid)
      if (status.eq.nf_noerr) then
        print *,'defining qv0'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"qv0",nf_float,3,(/1,2,3/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_put_att_int(ncid1,varid,"no_button",nf_int,1,int(1))
        status = nf_enddef(ncid1)
        print *,'writing qv0'
        status = nf_inq_varid(ncid1,"qv0",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"qv0",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'u',varid)
      if (status.eq.nf_noerr) then
        print *,'defining u'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"u",nf_float,4,(/4,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing u'
        status = nf_inq_varid(ncid1,"u",varid)
        allocate(var3D(mxp1,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"u",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop+1,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,
     *              var3D(1:mxp1,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'v',varid)
      if (status.eq.nf_noerr) then
        print *,'defining v'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"v",nf_float,4,(/1,5,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                         1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing v'
        status = nf_inq_varid(ncid1,"v",varid)
        allocate(var3D(mx,myp1,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"v",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *            var3D(istart:istop,jstart:jstop+1,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:myp1,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'w',varid)
      if (status.eq.nf_noerr) then
        print *,'defining w'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"w",nf_float,4,(/1,2,6,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",3,"m/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,wgz(1))
        status = nf_enddef(ncid1)
        print *,'writing w'
        status = nf_inq_varid(ncid1,"w",varid)
        allocate(var3D(mx,my,nzp1))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"w",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nzp1))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nzp1))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'p',varid)
      if (status.eq.nf_noerr) then
        print *,'defining p'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"p",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1,"K")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing p'
        status = nf_inq_varid(ncid1,"p",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"p",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'th',varid)
      if (status.eq.nf_noerr) then
        print *,'defining th'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"th",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1,"K")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing th'
        status = nf_inq_varid(ncid1,"th",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"th",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'thet',varid)
      if (status.eq.nf_noerr) then
        print *,'defining th'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"th",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1,"K")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing th'
        status = nf_inq_varid(ncid1,"th",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"thet",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'ppi',varid)
      if (status.eq.nf_noerr) then
        print *,'defining ppi'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"ppi",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",1," ")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing ppi'
        status = nf_inq_varid(ncid1,"ppi",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"ppi",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *              var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'rho',varid)
      if (status.eq.nf_noerr) then
        print *,'defining rho'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"rho",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",8,"kg/m/m/m")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing rho'
        status = nf_inq_varid(ncid1,"rho",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"rho",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'qv',varid)
      if (status.eq.nf_noerr) then
        print *,'defining qv'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"qv",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing qv'
        status = nf_inq_varid(ncid1,"qv",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"qv",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'qc',varid)
      if (status.eq.nf_noerr) then
        print *,'defining qc'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"qc",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing qc'
        status = nf_inq_varid(ncid1,"qc",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"qc",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'qr',varid)
      if (status.eq.nf_noerr) then
        print *,'defining qr'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"qr",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"kg/kg")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing qr'
        status = nf_inq_varid(ncid1,"qr",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"qr",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'km',varid)
      if (status.eq.nf_noerr) then
        print *,'defining km'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"km",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"m^2/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing km'
        status = nf_inq_varid(ncid1,"km",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"km",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'kh',varid)
      if (status.eq.nf_noerr) then
        print *,'defining kh'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"kh",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"m^2/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing kh'
        status = nf_inq_varid(ncid1,"kh",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"kh",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      end if

      status = nf_inq_varid(ncid(1),'tke',varid)
      if (status.eq.nf_noerr) then
        print *,'defining tke'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"tke",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",5,"m/s/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing tke'
        status = nf_inq_varid(ncid1,"tke",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"tke",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *              var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      endif

      status = nf_inq_varid(ncid(1),'fp',varid)
      if (status.eq.nf_noerr) then
        print *,'defining fp'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"fp",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",8,"kg/m/m/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing fp'
        status = nf_inq_varid(ncid1,"fp",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"fp",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,
     *              var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      endif

      status = nf_inq_varid(ncid(1),'lh',varid)
      if (status.eq.nf_noerr) then
        print *,'defining lh'
        status = nf_redef(ncid1)
        status = nf_def_var(ncid1,"lh",nf_float,4,(/1,2,3,7/),varid)
        status = nf_put_att_text(ncid1,varid,"units",3,"K/s")
        status = nf_put_att_real(ncid1,varid,"x_min",nf_real,
     *                           1,x_min1+0.5*dx)
        status = nf_put_att_real(ncid1,varid,"y_min",nf_real,
     *                           1,y_min1+0.5*dy)
        status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
        status = nf_enddef(ncid1)
        print *,'writing lh'
        status = nf_inq_varid(ncid1,"lh",varid)
        allocate(var3D(mx,my,nz))
        var3D(:,:,:)=0.0
        do i=1,nfile
          istart=(x_min(i)-x_min1)/dx+1
          istop=istart+nx(i)-1
          jstart=(y_min(i)-y_min1)/dy+1
          jstop=jstart+ny(i)-1

          status = nf_inq_varid(ncid(i),"lh",varid1)
          status = nf_get_var_real(ncid(i),varid1,
     *             var3D(istart:istop,jstart:jstop,1:nz))
        enddo
        status = nf_put_var_real(ncid1,varid,var3D(1:mx,1:my,1:nz))
        deallocate(var3D)
      endif

!-----------------------------------------------------------------------
! close files
!-----------------------------------------------------------------------
      print *, 'Closing files'
      do i=1,nfile
      status = nf_close(ncid(i)) 
      if (status /= nf_noerr) call disp_err(status)
      enddo

      status = nf_close(ncid1) 
      if (status /= nf_noerr) call disp_err(status)

      end program combine

!-----------------------------------------------------------------------
! subroutines
!-----------------------------------------------------------------------
      subroutine disp_err(status)
  ! converts error message to text form and prints it

      implicit none
!      include '/usr/local/netcdf-pgi/include/netcdf.inc'
      include '/usr/local/netcdf/include/netcdf.inc'

      integer, intent (in) :: status
  
      print *, trim(nf_strerror(status))
      stop "Stopped"
      end subroutine disp_err
