!----------------------------------------------------------------------------------------
! This program combines netCDF data from different processors from
! George Bryan's cm1 numerical model.
!
! Provided by:  Daniel Kirshbaum, University of Reading
! 
! Last modified:  2 April 2008
!
! Updated:  Early Sept 2008 by Lou Wicker, National Severe Storms Lab
!           ==> Converted it to F90 freeform style (needs compiling for 132 char lines)
!           ==> Created subroutines for reading in variables
!           ==> Made it more compatible with cm1r12 (vertical staggering of tke/kh/kv)
!           OTHER CHANGES
!           * I removed the "no_button" attributes (this used for GUI???)
!           * I added in dumping out all of the coordinates / horiz grid stretching not yet supported
!           * There is no terrain coordinate dump - its easy to add in
!----------------------------------------------------------------------------------------

 PROGRAM COMBINE_V2

  implicit none

  include '/usr/local/netcdf3-32/include/netcdf.inc' 

!-----------------------------------------
! input and auxiliary variables

  integer :: nfile,xtype,nwrite
  character(len=120), dimension(200) :: fname
  integer, dimension(200) :: ncid
  real,dimension(200) :: x_min,x_max,y_min,y_max
  integer, dimension(200) :: nx,nxp1,ny,nyp1

  character(len=120) :: fname1
  character(len=19)  :: fname2
  character(len=2)   :: fid
  integer            :: ncid1,ncid2
  integer            :: i, j, k, n,lenf, nVars2

!-----------------------------------------
! global vars for the original file

  integer            :: status, nDims, nVars, nGlobalAtts

!-----------------------------------------
! dimensions vars

  integer            :: DimID
  integer, dimension(nf_max_var_dims) :: DimLen
  character(nf_max_name) :: DimName

!-----------------------------------------
! attributes vars

  character(nf_max_name) :: AttName

!-----------------------------------------
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
  real, dimension(:), allocatable :: time,sgz,wgz,sgx,ugx,sgy,vgy
  real, dimension(:,:), allocatable :: zs,var2D
  real, dimension(:,:,:), allocatable :: var3D

!-----------------------------------------------------------------------
! read in the netCDF file name

  write(*,"(/,'Input number of netCDF files to be processed :')")
  read(*,*) nfile

  write(*,"(/,'Input output time index (nwrite) :')")
  read(*,*) nwrite

  print *, 'NFILE  = ',nfile
  print *, 'NWRITE = ',nwrite

  DO i = 1,nfile

    fname2 = 'cm1out_XXXX_XXXX.nc'
    write(fname2(8:11),"(i4.4)") i-1
    write(fname2(13:16),"(i4.4)") nwrite
    fname(i) = fname2
    write(6,"(/,'Reading ',A)") fname(i)

! open the netCDF file to read, added in call to disp_err to write out error message always

    CALL DISP_ERR( nf_open(fname(i),nf_nowrite,ncid(i)), .true.)

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
  ENDDO

96  format('nx = ',I5,5x,'nxp1 = ',I5,/,'ny = ',I5,5x,'nyp1 = ', I5,/,'nz = ',I5,5x,'nzp1 = ',I5,/,'nt = ',I5)

!-----------------------------------------------------------------------
! get dx,dy

   call disp_err( nf_get_att_real(ncid(1),nf_global,"x_delta",dx), .true. )
   call disp_err( nf_get_att_real(ncid(1),nf_global,"y_delta",dy), .true. )

   write(6,"('dx =',f15.5,5x,'dy = ',f15.5,5x)") dx,dy

! get x_min, x_max, y_min, y_max, z_min, z_max

   DO i=1,nfile
    status = nf_get_att_real(ncid(i),nf_global,"x_min",x_min(i))
    status = nf_get_att_real(ncid(i),nf_global,"x_max",x_max(i))
    status = nf_get_att_real(ncid(i),nf_global,"y_min",y_min(i))
    status = nf_get_att_real(ncid(i),nf_global,"y_max",y_max(i))
    status = nf_get_att_real(ncid(i),nf_global,"z_min",z_min)
    status = nf_get_att_real(ncid(i),nf_global,"z_max",z_max)
   ENDDO

   x_min1=minval(x_min)
   x_max1=maxval(x_max)
   y_min1=minval(y_min)
   y_max1=maxval(y_max)

   write(6,"('x_min =',f15.5,5x,'x_max = ',f15.5)") x_min1,x_max1
   write(6,"('y_min =',f15.5,5x,'y_max = ',f15.5)") y_min1,y_max1

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
   write(6,"(/,'nt / time(nt)', 2x,i3,2x,f8.1)") nt,time(nt)

! sgx/ugx  - just create a straight grid here

   allocate(sgx(mx))
   allocate(ugx(mxp1))
   sgx(:) = 0.0
   ugx(:) = 0.0
   DO i = 1,mx
    sgx(i) = x_min1 + (float(i)-0.5)*dx
    ugx(i) = x_min1 + (float(i)-1.0)*dx
   ENDDO
   ugx(mxp1) = (float(mxp1)-1.0)*dx

! sgy/vgy  - just create a straight grid here

   allocate(sgy(my))
   allocate(vgy(myp1))
   sgy(:) = 0.0
   vgy(:) = 0.0
   DO i = 1,my
    sgy(i) = y_min1 + (float(i)-0.5)*dy
    vgy(i) = y_min1 + (float(i)-1.0)*dy
   ENDDO
   vgy(myp1) = (float(myp1)-1.0)*dy

! sgz
   allocate(sgz(nz))
   sgz(:) = 0.0
   status = nf_inq_varid(ncid(1), 'zh', varid)
   status = nf_get_var_real(ncid(1), varid, sgz(1:nz))

! wgz
   allocate(wgz(nzp1))
   wgz(:) = 0.0
   status = nf_inq_varid(ncid(1), 'zf', varid)
   status = nf_get_var_real(ncid(1), varid, wgz(1:nzp1))
   ztop = wgz(nzp1)

!-----------------------------------------------------------------------
! Write output, necessary variables

   write(*,"(/,'Input name of the new netCDF file (name.nc):')") 
   read(*,*) fname1
        
   call disp_err( nf_create(fname1,nf_clobber,ncid1), .true. )

   status = nf_def_dim(ncid1,"nx",mx,DimID)
   status = nf_def_dim(ncid1,"ny",my,DimID)
   status = nf_def_dim(ncid1,"nz",nz,DimID)
   status = nf_def_dim(ncid1,"nxp1",mxp1,DimID)
   status = nf_def_dim(ncid1,"nyp1",myp1,DimID)
   status = nf_def_dim(ncid1,"nzp1",nzp1,DimID)
   status = nf_def_dim(ncid1,"time",nf_UNLIMITED,DimID)
   status = nf_def_dim(ncid1,"one",1,DimID)

   write(*,"(/,'Writing global attributes')")

   status = nf_put_att_real(ncid1,nf_global,"x_min",nf_real,1,x_min1)
   status = nf_put_att_real(ncid1,nf_global,"x_max",nf_real,1,x_max1)
   status = nf_put_att_real(ncid1,nf_global,"x_delta",nf_real,1,dx)
   status = nf_put_att_text(ncid1,nf_global,"x_units",2,"km")
   status = nf_put_att_text(ncid1,nf_global,"x_label",1,"x")
   status = nf_put_att_text(ncid1,nf_global,"x_display_units",2,"km")
   status = nf_put_att_real(ncid1,nf_global,"y_min",nf_real,1,y_min1)
   status = nf_put_att_real(ncid1,nf_global,"y_max",nf_real,1,y_max1)
   status = nf_put_att_real(ncid1,nf_global,"y_delta",nf_real,1,dy)
   status = nf_put_att_text(ncid1,nf_global,"y_units",2,"km")
   status = nf_put_att_text(ncid1,nf_global,"y_label",1,"y")
   status = nf_put_att_text(ncid1,nf_global,"y_display_units",2,"km")
   status = nf_put_att_real(ncid1,nf_global,"z_min",nf_real,1,z_min)
   status = nf_put_att_real(ncid1,nf_global,"z_max",nf_real,1,z_max)
   status = nf_put_att_text(ncid1,nf_global,"z_units",2,"km")
   status = nf_put_att_text(ncid1,nf_global,"z_label",1,"z")
   status = nf_put_att_text(ncid1,nf_global,"z_display_units",2,"km")
   status = nf_put_att_text(ncid1,nf_global,"runname",1,fname(1))

   status = nf_redef(ncid1)
   status = nf_def_var(ncid1,"f_cor",nf_real,1,8,varid)
   status = nf_put_att_text(ncid1,varid,"units",3,"1/s")
   status = nf_put_att_text(ncid1,varid,"def",18,"coriolis parameter")
   status = nf_enddef(ncid1)
   var1D  = 0.0
   status = nf_inq_varid(ncid(1),"f_cor",varid1)
   status = nf_get_var_real(ncid(1),varid1,var1D)
   status = nf_put_var_real(ncid1,varid,var1D)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"ztop",nf_real,1,8,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",12,"domain depth")
   status = nf_enddef(ncid1)
   status = nf_put_var_real(ncid1,varid,ztop)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"time",nf_float,1,7,varid)
   status = nf_put_att_text(ncid1,varid,"units",1,"s")
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,nt,time)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"xh",nf_float,1,1,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",11,"staggered x") 
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,mx,sgx)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"xf",nf_float,1,4,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",13,"unstaggered x") 
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,mxp1,ugx)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"yh",nf_float,1,2,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",11,"staggered y") 
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,my,sgy)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"yf",nf_float,1,5,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",13,"unstaggered y") 
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,myp1,vgy)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"zh",nf_float,1,3,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",11,"staggered z") 
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,nz,sgz)
   status = nf_redef(ncid1)

   status = nf_def_var(ncid1,"zf",nf_float,1,6,varid)
   status = nf_put_att_text(ncid1,varid,"units",2,"km")
   status = nf_put_att_text(ncid1,varid,"def",13,"unstaggered z") 
   status = nf_enddef(ncid1)
   status = nf_put_vara_real(ncid1,varid,1,nzp1,wgz)
   status = nf_redef(ncid1)

!-----------------------------------------------------------------------
! Read/Write 2D variables

   CALL READ_WRITE_VAR2D("zbot_p", 2, (/1,2/), "m",   0, 0)
   CALL READ_WRITE_VAR2D("thflux", 2, (/1,2/), "K/m/m/s",   0, 0)
   CALL READ_WRITE_VAR2D("qvflux", 2, (/1,2/), "1/m/m/s",   0, 0)
   CALL READ_WRITE_VAR2D("rain",   2, (/1,2/), "mm",   0, 0)

   write(6,"(/,'Completed the 1D and 2D fields',/)")

!-----------------------------------------------------------------------
!  Read/Write 3D variables
! 
! LJW:  I created a subroutine to handle to do the read/write to create
!       code that was a bit more manageable

   CALL READ_WRITE_VAR3D("u0",      3, (/4,2,3/),   "m/s",      1, 0, 0)
   CALL READ_WRITE_VAR3D("v0",      3, (/1,5,3/),   "m/s",      0, 1, 0)
   CALL READ_WRITE_VAR3D("pi0",     3, (/1,2,3/),   " ",        0, 0, 0)
   CALL READ_WRITE_VAR3D("th0",     3, (/1,2,3/),   "K",        0, 0, 0)
   CALL READ_WRITE_VAR3D("rho0",    3, (/1,2,3/),   "kg/m/m/m", 0, 0, 0)
   CALL READ_WRITE_VAR3D("qv0",     3, (/1,2,3/),   "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("u",       4, (/4,2,3,7/), "m/s",      1, 0, 0)
   CALL READ_WRITE_VAR3D("v",       4, (/1,5,3,7/), "m/s",      0, 1, 0)
   CALL READ_WRITE_VAR3D("w",       4, (/1,2,6,7/), "m/s",      0, 0, 1)
   CALL READ_WRITE_VAR3D("uinterp", 4, (/1,2,3,7/), "m/s",      0, 0, 0)
   CALL READ_WRITE_VAR3D("vinterp", 4, (/1,2,3,7/), "m/s",      0, 0, 0)
   CALL READ_WRITE_VAR3D("winterp", 4, (/1,2,3,7/), "m/s",      0, 0, 0)
   CALL READ_WRITE_VAR3D("prspert", 4, (/1,2,3,7/), " ",        0, 0, 0)
   CALL READ_WRITE_VAR3D("th",      4, (/1,2,3,7/), "K",        0, 0, 0)
   CALL READ_WRITE_VAR3D("thpert",  4, (/1,2,3,7/), "K",        0, 0, 0)
   CALL READ_WRITE_VAR3D("ppi",     4, (/1,2,3,7/), " ",        0, 0, 0)
   CALL READ_WRITE_VAR3D("rho",     4, (/1,2,3,7/), "kg/m/m/m", 0, 0, 0)
   CALL READ_WRITE_VAR3D("qv",      4, (/1,2,3,7/), "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("qc",      4, (/1,2,3,7/), "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("qr",      4, (/1,2,3,7/), "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("qi",      4, (/1,2,3,7/), "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("qs",      4, (/1,2,3,7/), "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("qg",      4, (/1,2,3,7/), "kg/kg",    0, 0, 0)
   CALL READ_WRITE_VAR3D("km",      4, (/1,2,6,7/), "m^2/s",    0, 0, 1)
   CALL READ_WRITE_VAR3D("tke",     4, (/1,2,6,7/), "m/s/s",    0, 0, 1)

!-----------------------------------------------------------------------
! close files
!-----------------------------------------------------------------------
   DO i=1,nfile
     CALL DISP_ERR(nf_close(ncid(i)),.true.)
   ENDDO

   CALL DISP_ERR( nf_close(ncid1) , .true. )

   write(6,"(/,'Closed the files',/)")

  CONTAINS

   SUBROUTINE READ_WRITE_VAR2D(name, dim, shape, units, is, js)   ! Currently only works correctly on X-Y variables

     character*(*) :: name
     integer       :: dim
     integer       :: shape(:)
     character*(*) :: units
     integer       :: is, js

     IF( nf_inq_varid(ncid(1),name,varid) == nf_noerr) THEN
       write(6,"(/,'Read/Write VAR2D:  name/units = ',a, 2x, a, 2x, i3,/)") name, units
       status = nf_redef(ncid1)
       status = nf_def_var(ncid1,name,nf_float,dim,shape,varid)
       status = nf_put_att_text(ncid1,varid,"units",len(units),units)
       IF( is == 0 ) THEN
         status = nf_put_att_real(ncid1,varid,"x_min",nf_real,1,x_min1+0.5*dx)
       ELSE
         status = nf_put_att_real(ncid1,varid,"x_min",nf_real,1,x_min1)
       ENDIF
       IF( js == 0 ) THEN
         status = nf_put_att_real(ncid1,varid,"y_min",nf_real,1,y_min1+0.5*dz)
       ELSE
         status = nf_put_att_real(ncid1,varid,"y_min",nf_real,1,y_min1)
       ENDIF

       status = nf_enddef(ncid1)
       status = nf_inq_varid(ncid1,name,varid)
       allocate(var2D(mx+is,my+js))
       var2D(:,:)=0.0

       DO i = 1,nfile
        istart=(x_min(i)-x_min1)/dx+1
        istop=istart+nx(i)-1
        jstart=(y_min(i)-y_min1)/dy+1
        jstop=jstart+ny(i)-1
        status = nf_inq_varid(ncid(i),name,varid1)
        status = nf_get_var_real(ncid(i),varid1,var2D(istart:istop+is,jstart:jstop+js))
       ENDDO
        status = nf_put_var_real(ncid1,varid,var2D(1:mx+is,1:my+js))
        deallocate(var2D)
      ELSE

       write(6,"(/,'Read/Write VAR2D ',a,' DID NOT FIND VARIABLE',/)") name

      ENDIF

   RETURN
   END SUBROUTINE

   SUBROUTINE READ_WRITE_VAR3D(name, dim, shape, units, is, js, ks)

     character*(*) :: name
     integer       :: dim
     integer       :: shape(:)
     character*(*) :: units
     integer       :: is, js, ks

     IF( nf_inq_varid(ncid(1),name,varid) == nf_noerr) THEN
       write(6,"(/,'Read/Write VAR3D:  name/units = ',a, 2x, a, 2x, i3,/)") name, units
       status = nf_redef(ncid1)
       status = nf_def_var(ncid1,name,nf_float,dim,shape,varid)
       status = nf_put_att_text(ncid1,varid,"units",len(units),units)
       IF( is == 0 ) THEN
         status = nf_put_att_real(ncid1,varid,"x_min",nf_real,1,x_min1+0.5*dx)
       ELSE
         status = nf_put_att_real(ncid1,varid,"x_min",nf_real,1,x_min1)
       ENDIF
       IF( js == 0 ) THEN
         status = nf_put_att_real(ncid1,varid,"y_min",nf_real,1,y_min1+0.5*dz)
       ELSE
         status = nf_put_att_real(ncid1,varid,"y_min",nf_real,1,y_min1)
       ENDIF
       IF( ks == 0 ) THEN
         status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,sgz(1))
       ELSE
         status = nf_put_att_real(ncid1,varid,"z_min",nf_real,1,wgz(1))
       ENDIF

       status = nf_enddef(ncid1)
       status = nf_inq_varid(ncid1,name,varid)
       allocate(var3D(mx+is,my+js,nz+ks))
       var3D(:,:,:)=0.0

       DO i = 1,nfile
        istart=(x_min(i)-x_min1)/dx+1
        istop=istart+nx(i)-1
        jstart=(y_min(i)-y_min1)/dy+1
        jstop=jstart+ny(i)-1
        status = nf_inq_varid(ncid(i),name,varid1)
        status = nf_get_var_real(ncid(i),varid1,var3D(istart:istop+is,jstart:jstop+js,1:nz+ks))
       ENDDO
        status = nf_put_var_real(ncid1,varid,var3D(1:mx+is,1:my+js,1:nz+ks))
        deallocate(var3D)
      ELSE

       write(6,"(/,'Read/Write VAR3D ',a,' DID NOT FIND VARIABLE',/)") name

      ENDIF

   RETURN
   END SUBROUTINE

  END PROGRAM COMBINE_V2

!-----------------------------------------------------------------------
! SUBROUTINE DISP_ERROR
! converts error message to text form and prints it
!-----------------------------------------------------------------------
  SUBROUTINE DISP_ERR(status,die_on_error)

   implicit none
   include '/usr/local/netcdf3-32/include/netcdf.inc'

   integer, intent (in) :: status
   logical, intent (in) :: die_on_error
  
   IF (status /= nf_noerr) THEN
     print *, trim(nf_strerror(status))
    IF( die_on_error) THEN
     stop "STATUS returned an error, program stopped"
    ENDIF
   ENDIF

  END SUBROUTINE DISP_ERR