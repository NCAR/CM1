  MODULE restart_read_module

  implicit none

  private
  public :: restart_read


  logical :: dointerp


  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine restart_read(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
                              trecs,trecw,arecs,arecw,                             &
                              nwrite,nwritet,nwritea,nwriteh,nrst,nstatout,        &
                              num_soil_layers,nrad2d,                              &
                              avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp, &
                              dt,dtlast,mtime,ndt,adt,adtlast,acfl,dbldt,mass1,    &
                              stattim,taptim,rsttim,radtim,prcltim,                &
                              qbudget,asq,bsq,qname,                               &
                              xfref,xhref,yfref,yhref,xh,xf,yh,yf,zh,zf,sigma,sigmaf,zs, &
                              th0,prs0,pi0,rho0,qv0,u0,v0,                         &
                              rain,raincv,sws,svs,sps,srs,sgs,sus,shs,             &
                              tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,thflux,qvflux,      &
                              prate,ustt,ut,vt,st,                                 &
                              radbcw,radbce,radbcs,radbcn,                         &
                              rho,prs,ua,dumu,va,dumv,wa,ppi,tha,qa,tkea,          &
                              swten,lwten,thraten,radsw,rnflx,radswnet,radlwin,rad2d,      &
                              effc,effi,effs,effr,effg,effis,                      &
                              thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,       &
                              lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,xland,      &
                              hpbl,wspd,psim,psih,gz1oz0,br,                       &
                              CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                              MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                              CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                              gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta,tslb,    &
                              tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                              qpten,qtten,qvten,qcten,pta,pdata,ploc,ppx,          &
                              tdiag,qdiag,phi1,phi2,                               &
                   tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl,       &
                   qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,            &
                   edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,          &
                   sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D,                  &
                   vdfg,maxmf,nupdraft,ktop_plume,                           &
                              tke_myj,el_myj,mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj,lowlyr,ivgtyp, &
                              thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten, &
                              icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy,adaptmovetim,mvrec,nwritemv,ug,vg, &
                              gamk,kmw,ufw,vfw,u1b,v1b,                                &
                              ntavg,rtavg,tavg,timavg,sfctavg,sfctimavg,dumsfc       , &
                              dum1,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,restarted,restart_prcl)
          ! end_restart_read
      use input
      use constants
      use lsnudge_module
      use goddard_module, only : consat,consat2
      use lfoice_module, only : lfoice_init
#ifdef MPI
      use mpi
#endif
#ifdef NETCDF
      use netcdf
      use writeout_nc_module, only : disp_err
#endif

      use restart_write_module

      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the reading of restart files
      !----------------------------------------------------------

      integer, intent(inout) :: nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer, intent(inout) :: nwrite,nwritet,nwritea,nwriteh,nrst,nstatout
      integer, intent(in) :: num_soil_layers,nrad2d
      double precision, intent(inout) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(inout) :: dt,dtlast
      integer, intent(inout) :: ndt
      double precision, intent(inout) :: adt,adtlast,acfl,dbldt
      double precision, intent(inout) :: mass1
      double precision, intent(inout) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      character(len=3), intent(in), dimension(maxq) :: qname
      real, dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
      real, dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, dimension(kb:ke) :: sigma
      real, dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,prs0,pi0,rho0,qv0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,raincv,sws,svs,sps,srs,sgs,sus,shs
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,xland,psfc,thflux,qvflux,prate,ustt,ut,vt,st
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: dumu
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: dumv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,thraten,effc,effi,effs,effr,effg,effis
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(inout), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppx
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl, &
           qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,  &
           sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: vdfg,maxmf
      integer, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: nupdraft,ktop_plume
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj,el_myj
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj
      integer, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: lowlyr,ivgtyp
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, intent(in   ) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(in   ) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      integer, intent(inout) :: icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      double precision, intent(inout) :: domainlocx,domainlocy
      double precision, intent(inout) :: adaptmovetim
      integer, intent(inout) :: mvrec,nwritemv
      real, intent(inout), dimension(kb:ke) :: ug,vg
      real, intent(inout), dimension(kb:ke) :: gamk,kmw,ufw,vfw,u1b,v1b
      integer, intent(inout), dimension(ntim) :: ntavg
      double precision, intent(inout), dimension(ntim) :: rtavg
      double precision, intent(inout), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntim,ntavr) :: tavg
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      double precision, intent(inout), dimension(ibta:ieta,jbta:jeta,ntim,nsfctavr) :: sfctavg
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg
      real, intent(inout), dimension(ib:ie,jb:je) :: dumsfc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
      logical, intent(in) :: restarted
      logical, intent(inout) :: restart_prcl

      character(len=maxstring) fname
      character(len=8) :: text1
      character(len=6) :: aname
      integer :: i,j,k,kk,n,np,nt,nvar,nread,reqs,orecs,orecu,orecv,orecw,ndum,ifoo
      integer :: ncid,time_index,old_format
      double precision, dimension(nbudget,0:numprocs-1) :: sbudget
      double precision, dimension(numq,0:numprocs-1) :: csq,dsq
      real, dimension(:,:), allocatable :: pfoo
      real, dimension(:), allocatable :: dumx,dumy,dumz
#ifdef MPI
      integer :: proc,index,count,req1,req2,req3,reqp
#endif
#ifdef NETCDF
      integer :: varid,ncstatus
#endif
      integer :: nxold,nyold,nzold,ngxyold,iold,jold,kold,nxr,nyr,ketaold
      integer :: wbcold,ebcold,sbcold,nbcold,bbcold,tbcold
      real :: maxxold,maxyold,maxzold
      real, dimension(:), allocatable :: xhref0,xfref0,yhref0,yfref0,sigma0,sigmaf0
      real, dimension(:,:), allocatable :: datk1,datk2
      integer, dimension(:), allocatable :: ix,iy,iz,ixu,iyv,izw
      real, dimension(:), allocatable :: rx,ry,rz,rxu,ryv,rzw

!-----------------------------------------------------------------------

  rformat:  &
  IF( restart_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

  rftype:  &
  IF( restart_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_x.dat'

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='old',err=778)

      if(dowr) write(outfile,*)

    ENDIF

  ELSEIF( restart_filetype.eq.2 )THEN  rftype

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_x.dat'
      write(fname( 8:13),101) rstnum
101   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='old',err=778)

      if(dowr) write(outfile,*)

    ENDIF

  ELSEIF( restart_filetype.eq.3 )THEN  rftype

    !------------------
    ! one restart file per node (cm1r17 format):
    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_YYYYYY.dat'

      write(fname( 8:13),102) mynode
      write(fname(15:20),102) rstnum
102   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname
      if(dowr) write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='old')
    ENDIF
  ELSE  rftype
    stop 12389
  ENDIF  rftype

#ifdef NETCDF
  ELSEIF( restart_format.eq.2 )THEN  rformat
    ! netcdf format:

    if( myid.eq.0 )then

      do i=1,maxstring
        string(i:i) = ' '
      enddo

    IF(     restart_filetype.eq.1 )THEN
      string = 'cm1rst.nc'
      time_index = rstnum
    ELSEIF( restart_filetype.eq.2 )THEN
      string = 'cm1rst_XXXXXX.nc'
      write(string(8:13),100) rstnum
100   format(i6.6)
      time_index = 1
    ENDIF
      if(myid.eq.0) print *,'  string = ',string

      call disp_err( nf90_open( path=string , mode=nf90_nowrite , ncid=ncid ) , .true. )

    endif

#endif
  ELSE  rformat

    if( myid.eq.0 )then
      print *
      print *,'  unrecognized value for restart_format '
      print *
      print *,'      restart_format = ',restart_format
      print *
    endif
#ifdef MPI
    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
    call stopcm1

  ENDIF  rformat

!---------------------------------------------------------------
! metadata:

  interp_on_restart = .false.
  dointerp = .false.
  nxold = 1
  nyold = 1
  nzold = 1
  ngxyold = 1

  IF( restart_format.eq.1 )THEN
    myid0only:  &
    IF(myid.eq.0)THEN
      ! only processor 0 has these variables:
      read(50,err=7001) cm1rversion
      print *,'  cm1rversion = ',cm1rversion
! kludge !
!!!      if( cm1rversion.lt.20.999 ) go to 7001
      read(50) nxold
      read(50) nyold
      read(50) nzold
      read(50) ngxyold
      read(50) wbcold
      read(50) ebcold
      read(50) sbcold
      read(50) nbcold
      read(50) bbcold
      read(50) tbcold
      if( nxold.ne.nx .or. nyold.ne.ny .or. nzold.ne.nz )then
        interp_on_restart = .true.
        dointerp = .true.
        print *
        print *,'  Different grids detected: '
        print *
        print *,'  nxold,nx = ',nxold,nx
        print *,'  nyold,ny = ',nyold,ny
        print *,'  nzold,nz = ',nzold,nz
        print *
        print *,'  Interpolating to new grid ... '
        print *
      else
        interp_on_restart = .false.
        nxold = 1
        nyold = 1
        nzold = 1
      endif
      read(50) maxxold
      read(50) maxyold
      read(50) maxzold
      if( interp_on_restart )then
        if( (maxx-minx).gt.maxxold .or.  &
            (maxy-miny).gt.maxyold .or.  &
             maxz      .gt.maxzold )then
          print *
          print *,'  new domain is too large '
          print *
          print *,'  old,new x size = ',maxxold,maxx-minx
          print *,'  old,new y size = ',maxyold,maxy-miny
          print *,'  old,new z size = ',maxzold,maxz
          print *
          print *,'  ...  stopping cm1  ... '
          print *
          call stopcm1
        endif
      endif
      read(50) nstep
      read(50) srec
      read(50) sirec
      read(50) urec
      read(50) vrec
      read(50) wrec
      read(50) nrec
      read(50) mrec
      read(50) prec
      read(50) trecs
      read(50) trecw
      read(50) arecs
      read(50) arecw
      read(50) mvrec
      read(50) nwrite
      read(50) nwritet
      read(50) nwritea
      read(50) nwritemv
      read(50) nwriteh
      read(50) nrst
      read(50) nstatout
      read(50) ndt
      read(50) icenter
      read(50) jcenter
      read(50) old_format
      read(50) dt
      read(50) dtlast
      read(50) xcenter
      read(50) ycenter
      read(50) umove
      read(50) umovem1
      read(50) umovem2
      read(50) vmove
      read(50) vmovem1
      read(50) vmovem2
      read(50) domainlocx
      read(50) domainlocy
      read(50) adaptmovetim
      read(50) cflmax
      read(50) mtime
      read(50) stattim
      read(50) taptim
      read(50) rsttim
      read(50) radtim
      read(50) prcltim
      read(50) adt
      read(50) adtlast
      read(50) acfl
      read(50) dbldt
      read(50) mass1
      read(50) avgsfcu
      read(50) avgsfcv
      read(50) avgsfcs
      read(50) avgsfcsu
      read(50) avgsfcsv
      read(50) avgsfct
      read(50) avgsfcq
      read(50) avgsfcp
      allocate( xhref0(1-ngxyold:nxold+ngxyold+1) )
      xhref0 = 0.0
      allocate( xfref0(1-ngxyold:nxold+ngxyold+1) )
      xfref0 = 0.0
      allocate( yhref0(1-ngxyold:nyold+ngxyold+1) )
      yhref0 = 0.0
      allocate( yfref0(1-ngxyold:nyold+ngxyold+1) )
      yfref0 = 0.0
      allocate( sigma0(0:nzold+1) )
      sigma0 = 0.0
      allocate( sigmaf0(0:nzold+2) )
      sigmaf0 = 0.0
      read(50) xhref0
      read(50) xfref0
      read(50) yhref0
      read(50) yfref0
      read(50) sigma0
      read(50) sigmaf0

      print *,'  Depth of old domain: '
      print *,'  sigma0 = ',sigma0(1),sigma0(nzold)
      print *,'  sigmaf0 = ',sigmaf0(1),sigmaf0(nzold+1)
    ENDIF  myid0only

#ifdef MPI
      call MPI_BCAST(interp_on_restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dointerp         ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ngxyold,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nxold  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nyold  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nzold  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if( myid.ne.0 )then
        allocate( xhref0(1-ngxyold:nxold+ngxyold+1) )
        xhref0 = 0.0
        allocate( xfref0(1-ngxyold:nxold+ngxyold+1) )
        xfref0 = 0.0
        allocate( yhref0(1-ngxyold:nyold+ngxyold+1) )
        yhref0 = 0.0
        allocate( yfref0(1-ngxyold:nyold+ngxyold+1) )
        yfref0 = 0.0
        allocate( sigma0(0:nzold+1) )
        sigma0 = 0.0
        allocate( sigmaf0(0:nzold+2) )
        sigmaf0 = 0.0
      endif
      call MPI_BCAST(xhref0,nxold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xfref0,nxold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(yhref0,nyold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(yfref0,nyold+2*ngxyold+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigma0 ,nzold+2,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigmaf0,nzold+3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

      if(myid.eq.0) print *,'  allocating datk1 ... '
      allocate( datk1(0:nxold+2,0:nyold+2) )
      datk1 = 0.0
      if(myid.eq.0) print *,'  allocating datk2 ... '
      allocate( datk2(0:nxold+2,0:nyold+2) )
      datk2 = 0.0

      allocate( ix(ni) )
      ix = 0
      allocate( iy(nj) )
      iy = 0
      allocate( iz(nk) )
      iz = 0
      allocate( ixu(ni+1) )
      ixu = 0
      allocate( iyv(nj+1) )
      iyv = 0
      allocate( izw(nk+1) )
      izw = 0

      allocate( rx(ni) )
      rx = 0.0
      allocate( ry(nj) )
      ry = 0.0
      allocate( rz(nk) )
      rz = 0.0
      allocate( rxu(ni+1) )
      rxu = 0.0
      allocate( ryv(nj+1) )
      ryv = 0.0
      allocate( rzw(nk+1) )
      rzw = 0.0

      if(myid.eq.0) print *,'  ... done allocating data '

      if( interp_on_restart .and. apmasscon.ge.1 )then
        if(myid.eq.0)then
          print *
          print *,'  NOTE: setting apmasscon to 0 '
          print *
        endif
        apmasscon = 0
      endif

      if( interp_on_restart )then
        !...................................................................
        ! Get interpolation info for scalar points:
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  xold,xnew,xold,rx: '
        iold = 0
        do i=1,ni
          do while( xhref0(iold+1).lt.xh(i) )
            iold = iold + 1
          enddo
          ix(i) = iold
          rx(i) = (xh(i)-xhref0(iold))/(xhref0(iold+1)-xhref0(iold))
          if(myid.eq.0) print *,iold,xhref0(iold),xh(i),xhref0(iold+1),rx(i)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  yold,ynew,yold,ry: '
        jold = 0
        do j=1,nj
          do while( yhref0(jold+1).lt.yh(j) )
            jold = jold + 1
          enddo
          iy(j) = jold
          ry(j) = (yh(j)-yhref0(jold))/(yhref0(jold+1)-yhref0(jold))
          if(myid.eq.0) print *,jold,yhref0(jold),yh(j),yhref0(jold+1),ry(j)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  zold,znew,zold,rz: '
        kold = 0
        do k=1,nk
          do while( sigma0(kold+1).lt.sigma(k) )
            kold = kold + 1
          enddo
          iz(k) = kold
          rz(k) = (sigma(k)-sigma0(kold))/(sigma0(kold+1)-sigma0(kold))
          if(myid.eq.0) print *,kold,sigma0(kold),sigma(k),sigma0(kold+1),rz(k)
        enddo
        if(myid.eq.0) print *
        !...................................................................
        ! Get interpolation info for velocity points:
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  xold,xnew,xold,rx: '
        iold = 0
        do i=1,ni+1
          do while( xfref0(iold+1).lt.xf(i) )
            iold = iold + 1
          enddo
          ixu(i) = iold
          rxu(i) = (xf(i)-xfref0(iold))/(xfref0(iold+1)-xfref0(iold))
          if(myid.eq.0) print *,iold,xfref0(iold),xf(i),xfref0(iold+1),rxu(i)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  yold,ynew,yold,ry: '
        jold = 0
        do j=1,nj+1
          do while( yfref0(jold+1).lt.yf(j) )
            jold = jold + 1
          enddo
          iyv(j) = jold
          ryv(j) = (yf(j)-yfref0(jold))/(yfref0(jold+1)-yfref0(jold))
          if(myid.eq.0) print *,jold,yfref0(jold),yf(j),yfref0(jold+1),ryv(j)
        enddo
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  zold,znew,zold,rz: '
        kold = 0
        do k=1,nk+1
          do while( sigmaf0(kold+1).lt.sigmaf(k) )
            kold = kold + 1
          enddo
          izw(k) = kold
          rzw(k) = (sigmaf(k)-sigmaf0(kold))/(sigmaf0(kold+1)-sigmaf0(kold))
          if(myid.eq.0) print *,kold,sigmaf0(kold),sigmaf(k),sigmaf0(kold+1),rzw(k)
        enddo
        if(myid.eq.0) print *
        !...................................................................
      endif

#ifdef NETCDF
  ELSEIF( restart_format.eq.2 )THEN

    IF(myid.eq.0)THEN

      call disp_err( nf90_inq_varid(ncid,"cm1rversion",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,cm1rversion,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstep",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nstep,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"srec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,srec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"sirec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,sirec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"urec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,urec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"wrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,wrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,prec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,trecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecw",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,trecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,arecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecw",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,arecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mvrec",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mvrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwrite",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwrite,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritet",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritet,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritea",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritea,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritemv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwritemv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrst",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nrst,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstatout",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,nstatout,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ndt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,ndt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,icenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,jcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"old_format",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,old_format,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dtlast",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,xcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,ycenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,umove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umovem1",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,umovem1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umovem2",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,umovem2,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vmove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmovem1",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vmovem1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmovem2",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,vmovem2,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocx",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,domainlocx,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocy",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,domainlocy,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adaptmovetim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adaptmovetim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"cflmax",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,cflmax,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mtime",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"stattim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,stattim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"taptim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,taptim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"rsttim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,rsttim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"radtim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,radtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prcltim",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,prcltim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adtlast",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,adtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"acfl",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,acfl,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dbldt",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,dbldt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mass1",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,mass1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcu",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcs",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsu",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcsu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsv",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcsv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfct",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfct,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcq",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcq,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcp",varid) , .true. )
      call disp_err( nf90_get_var(ncid,varid,avgsfcp,(/time_index/)) , .true. )

    ENDIF
#endif
  ENDIF

#ifdef MPI
      ! communicate to all other processors:
      call MPI_BCAST(nstep  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(srec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sirec  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(urec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trecs  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trecw  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(arecs  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(arecw  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mvrec  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwrite ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwritet,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwritea,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwritemv,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwriteh ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nrst   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nstatout,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ndt    ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icenter,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jcenter,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(old_format,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt     ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtlast ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xcenter,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ycenter,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umove  ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umovem1,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umovem2,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmove  ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmovem1,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmovem2,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(domainlocx  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(domainlocy  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adaptmovetim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cflmax ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mtime  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stattim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(taptim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rsttim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(radtim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcltim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adt    ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adtlast,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(acfl   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dbldt  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mass1  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcsu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcsv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfct,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(avgsfcp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
!---------------------------------------------------------------

  if(myid.eq.0) print *,'  nxold,nyold,nzold = ',nxold,nyold,nzold

  if( .not. allocated(datk1) )    allocate( datk1(0:nxold+1,0:nyold+1) )
  if( .not. allocated(datk2) )    allocate( datk2(0:nxold+1,0:nyold+1) )

  if( .not. allocated(sigma0) )   allocate( sigma0(0:nzold) )
  if( .not. allocated(sigmaf0) )  allocate( sigmaf0(0:nzold+2) )

  if( .not. allocated(ix) )       allocate( ix(ni) )
  if( .not. allocated(iy) )       allocate( iy(nj) )
  if( .not. allocated(iz) )       allocate( iz(nk) )

  if( .not. allocated(ixu) )      allocate( ixu(ni+1) )
  if( .not. allocated(iyv) )      allocate( iyv(nj+1) )
  if( .not. allocated(izw) )      allocate( izw(nk+1) )

  if( .not. allocated(rx) )       allocate( rx(ni) )
  if( .not. allocated(ry) )       allocate( ry(nj) )
  if( .not. allocated(rz) )       allocate( rz(nk) )

  if( .not. allocated(rxu) )      allocate( rxu(ni+1) )
  if( .not. allocated(ryv) )      allocate( ryv(nj+1) )
  if( .not. allocated(rzw) )      allocate( rzw(nk+1) )

!---------------------------------------------------------------
! budget variables:

    IF( myid.eq.0 )THEN

      IF( restart_format.eq.1 )THEN
        read(50) qbudget
        read(50) asq
        read(50) bsq
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"qbudget",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,qbudget,(/1,time_index/),(/nbudget,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"asq",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,asq,(/1,time_index/),(/numq,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"bsq",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,bsq,(/1,time_index/),(/numq,1/)) , .true. )
#endif
      ENDIF

    ELSE

      qbudget = 0.0
      asq = 0.0
      bsq = 0.0

    ENDIF

!---------------------------------------------------------------

  if( myid.eq.0 )then
    if( interp_on_restart )then
      allocate( dumz(0:nzold+1) )
      dumz = 0.0
      IF( restart_format.eq.1 )THEN
        read(50) dumz
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        stop 1111
#endif
      ENDIF
      do k=1,nk
        kk = iz(k)
        ug(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
      enddo
      IF( restart_format.eq.1 )THEN
        read(50) dumz
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
#endif
      ENDIF
      do k=1,nk
        kk = iz(k)
        vg(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
      enddo
      deallocate( dumz )
    else
      IF( restart_format.eq.1 )THEN
        read(50) ug
        read(50) vg
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"ug",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,ug(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"vg",varid) , .true. )
        call disp_err( nf90_get_var(ncid,varid,vg(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
#endif
      ENDIF
    endif
  endif

#ifdef MPI
    call MPI_BCAST(ug(kb),ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vg(kb),ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

!-----------------------------------------------------------------------

  rf2:  &
  IF( restart_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

    if( interp_on_restart )then
      nxr = nxold
      nyr = nyold
    else
      nxr = nx
      nyr = ny
    endif

    if( myid.eq.0 ) print *,'  nxr,nyr = ',nxr,nyr

  rft2:  &
  IF( restart_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_s.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nxr*nyr,status='old',err=778)
      orecs = 1

      fname = 'cm1rst_u.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nxr+1)*nyr,status='old',err=778)
      orecu = 1

      fname = 'cm1rst_v.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nxr*(nyr+1),status='old',err=778)
      orecv = 1

      fname = 'cm1rst_w.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nxr*nyr,status='old',err=778)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.2 )THEN  rft2

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_s.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nxr*nyr,status='old',err=778)
      orecs = 1

      fname = 'cm1rst_XXXXXX_u.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nxr+1)*nyr,status='old',err=778)
      orecu = 1

      fname = 'cm1rst_XXXXXX_v.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nxr*(nyr+1),status='old',err=778)
      orecv = 1

      fname = 'cm1rst_XXXXXX_w.dat'
      write(fname( 8:13),101) rstnum
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nxr*nyr,status='old',err=778)
      orecw = 1

      if(dowr) write(outfile,*)

    ENDIF

  ENDIF  rft2

  ENDIF  rf2

!---------------------------------------------------------------
! standard 2D:

      n = 1
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    if( nrain.eq.2 )then
      n = 2
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
      call  readr(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------
! standard 3D:

      call  readr(ni,nj,1,nk,nx,ny,rho(ib,jb,1),'rho     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,prs(ib,jb,1),'prs     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'ua      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ixu,iy,iz,rxu,ry,rz,nxold+1,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
      call  readr(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'va      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iyv,iz,rx,ryv,rz,nxold,nyold+1,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
      call  readr(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'wa      ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz+1,ix,iy,izw,rx,ry,rzw,nxold,nyold,nzold+1,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'ppi     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'tha     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,ppx(ib,jb,1),'ppx     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    if( psolver.eq.6 )then
      call  readr(ni,nj,1,nk,nx,ny,phi1(ib,jb,1),'phi1    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,phi2(ib,jb,1),'phi2    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
    IF(imoist.eq.1)THEN
    do n=1,numq
      text1 = '        '
      write(text1(1:3),156) qname(n)
156   format(a3)
      call  readr(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),text1     ,       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call  readr(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),'qpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),'qtten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),'qvten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),'qcten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
    if( idoles .and. iusetke )then
    if( .not. interp_on_restart )then  ! cm1r21.1: do not read tke when interpolating to a new grid !
      call  readr(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tkea    ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz+1,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold+1,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
    endif

!---------------------------------------------------------------
!  radiation:

      if(radopt.ge.1)then
        call  readr(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,swten(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,thraten(ib,jb,1),'thraten ',     &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radsw(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rnflx(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radswnet(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radlwin(i,j) = dum1(i,j,1)
        enddo
        enddo
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
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rad2d(i,j,n) = dum1(i,j,1)
        enddo
        enddo
        enddo
      endif
      if( radopt.ge.1 .and. ptype.eq.5 )then
        call  readr(ni,nj,1,nk,nx,ny,effc(ib,jb,1),'effc    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effi(ib,jb,1),'effi    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effs(ib,jb,1),'effs    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effr(ib,jb,1),'effr    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effg(ib,jb,1),'effg    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,nk,nx,ny,effis(ib,jb,1),'effis   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      if((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1))then
        !---- (1) ----!
      if(sfcmodel.ge.1)then
        call  readr(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,u1(ib,jb),'u1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,v1(ib,jb),'v1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,s1(ib,jb),'s1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,t1(ib,jb),'t1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,s10(ib,jb),'s10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psfc(ib,jb),'psfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( tbc.eq.3 )then
        call  readr(ni,nj,1,1,nx,ny,ustt(ib,jb),'ustt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  ut(ib,jb),'ut      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  vt(ib,jb),'vt      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,  st(ib,jb),'st      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      endif
      endif


      if(sfcmodel.ge.1)then
        !---- (2) ----!
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu_index',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          lu_index(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'kpbl2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          kpbl2d(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,wspd(ib,jb),'wspd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psim(ib,jb),'psim    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,psih(ib,jb),'psih    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),'gz1oz0  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CHS(ib,jb),'chs     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CHS2(ib,jb),'chs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CQS2(ib,jb),'cqs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CPMM(ib,jb),'cpmm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ZOL(ib,jb),'zol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),'mavail  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,MOL(ib,jb),'mol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,RMOL(ib,jb),'rmol    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,REGIME(ib,jb),'regime  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,LH(ib,jb),'lh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,FLHC(ib,jb),'flhc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,FLQC(ib,jb),'flqc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,QGH(ib,jb),'qgh     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CK(ib,jb),'ck      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CKA(ib,jb),'cka     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,CDA(ib,jb),'cda     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,USTM(ib,jb),'ustm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,QSFC(ib,jb),'qsfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,T2(ib,jb),'t2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,Q2(ib,jb),'q2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,TH2(ib,jb),'th2     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,EMISS(ib,jb),'emiss   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,THC(ib,jb),'thc     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,ALBD(ib,jb),'albd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,chklowq(ib,jb),'chklowq ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,capg(ib,jb),'capg    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,snowc(ib,jb),'snowc   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,fm(ib,jb),'fm      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,fh(ib,jb),'fh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,mznt(ib,jb),'mznt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,swspd(ib,jb),'swspd   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,wstar(ib,jb),'wstar   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,delta(ib,jb),'delta   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
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
        call  readr(ni,nj,1,1,nx,ny,tslb(ib,jb,n),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
      endif
      endif

      if(oceanmodel.eq.2)then
        call  readr(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,t0ml(ib,jb),'t0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,h0ml(ib,jb),'h0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call  readr(ni,nj,1,1,nx,ny,tmoml(ib,jb),'tmoml   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------

    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. sfcmodel.eq.6 )THEN
      if(myid.eq.0) print *,'  reading mynn vars ... '
      call  readr(ni,nj,1,nk,nx,ny,qke(ib,jb,1),'qke     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qke3d(ib,jb,1),'qke3d   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,tsq(ib,jb,1),'tsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qsq(ib,jb,1),'qsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,cov(ib,jb,1),'cov     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sh3d(ib,jb,1),'sh3d    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sm3d(ib,jb,1),'sm3d    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,el_pbl(ib,jb,1),'el_pbl  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qc_bl(ib,jb,1),'qc_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qi_bl(ib,jb,1),'qi_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,cldfra_bl(ib,jb,1),'cldfra_b',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_a(ib,jb,1),'edmf_a  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_w(ib,jb,1),'edmf_w  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_qt(ib,jb,1),'edmf_qt ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_thl(ib,jb,1),'edmf_thl',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_ent(ib,jb,1),'edmf_ent',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,edmf_qc(ib,jb,1),'edmf_qc ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sub_thl3d(ib,jb,1),'sub_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,sub_sqv3d(ib,jb,1),'sub_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,det_thl3d(ib,jb,1),'det_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,det_sqv3d(ib,jb,1),'det_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,thpten(ib,jb,1),'thpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvpten(ib,jb,1),'qvpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qcpten(ib,jb,1),'qcpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qipten(ib,jb,1),'qipten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,upten(ib,jb,1),'upten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,vpten(ib,jb,1),'vpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qnipten(ib,jb,1),'qnipten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qncpten(ib,jb,1),'qncpten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,vdfg(ib,jb),'vdfg    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,maxmf(ib,jb),'maxmf   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)






      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'nupdraft',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ktop_plu',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        nupdraft(i,j) = nint(dum1(i,j,1))
        ktop_plume(i,j) = nint(dum1(i,j,2))
      enddo
      enddo
      if(myid.eq.0) print *,'  ... done reading mynn vars '
    ENDIF

!---------------------------------------------------------------

    IF(ipbl.eq.6)THEN
      if(myid.eq.0) print *,'  reading myj vars ... '
      call  readr(ni,nj,1,nk,nx,ny,tke_myj(ib,jb,1),'tke_myj ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,el_myj(ib,jb,1),'el_myj  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF(sfcmodel.eq.7)THEN
      call  readr(ni,nj,1,1,nx,ny,mixht(ib,jb),'mixht   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,akhs(ib,jb),'akhs    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,akms(ib,jb),'akms    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,elflx(ib,jb),'elflx   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,ct(ib,jb),'ct      ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,snow(ib,jb),'snow    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,sice(ib,jb),'sice    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,thz0(ib,jb),'thz0    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,qz0(ib,jb),'qz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,uz0(ib,jb),'uz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,vz0(ib,jb),'vz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,th10(ib,jb),'th10    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,q10(ib,jb),'q10     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,z0base(ib,jb),'z0base  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,zntmyj(ib,jb),'zntmyj  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)






      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lowlyr  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ivgtyp  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        lowlyr(i,j) = nint(dum1(i,j,1))
        ivgtyp(i,j) = nint(dum1(i,j,2))
      enddo
      enddo
      if(myid.eq.0) print *,'  ... done reading myj vars '
    endif

!---------------------------------------------------------------
!  convection param.

    IF( cuparam.ge.1 )THEN
      call  readr(ni,nj,1,1,nx,ny,raincv(ib,jb,1),'raincv  ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if(nrain.eq.2)  &
      call  readr(ni,nj,1,1,nx,ny,raincv(ib,jb,2),'raincv2 ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,thcuten(ib,jb,1),'thcuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qvcuten(ib,jb,1),'qvcuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qccuten(ib,jb,1),'qccuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,qicuten(ib,jb,1),'qicuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,ucuten(ib,jb,1),'ucuten  ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call  readr(ni,nj,1,nk,nx,ny,vcuten(ib,jb,1),'vcuten  ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF

!---------------------------------------------------------------
!  2-part turbulence param:

      if( myid.eq.0 )then
        if( interp_on_restart )then

        if( restart_format.eq.1 )then
          allocate( dumz(0:nzold+1) )
          dumz = 0.0
          !---
          read(50) dumz
          do k=1,nk
            kk = iz(k)
            gamk(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(50) dumz
          do k=1,nk
            kk = iz(k)
            kmw(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(50) dumz
          do k=1,nk
            kk = iz(k)
            ufw(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(50) dumz
          do k=1,nk
            kk = iz(k)
            vfw(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(50) dumz
          do k=1,nk
            kk = iz(k)
            u1b(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          read(50) dumz
          do k=1,nk
            kk = iz(k)
            v1b(k) = (1.0-rz(k))*dumz(kk) + rz(k)*dumz(kk+1)
          enddo
          !---
          deallocate( dumz )
        else
          stop 2222
        endif

        else
          if( restart_format.eq.1 )then
            read(50) gamk
            read(50) kmw
            read(50) ufw
            read(50) vfw
            read(50) u1b
            read(50) v1b
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
#endif
          endif
        endif

      endif

#ifdef MPI
      call MPI_BCAST(gamk,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmw ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ufw ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vfw ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(u1b ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(v1b ,ke-kb+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

!---------------------------------------------------------------
!  passive tracers:

      if( myid.eq.0 )then
        if( restart_format.eq.1 )then
          read(50) nvar
#ifdef NETCDF
        elseif( restart_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( iptra.eq.0 ) nvar = 0
#endif
        endif
        print *,'  nvar_npt = ',nvar
      endif

#ifdef MPI
      call MPI_BCAST(nvar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      if( iptra.eq.1 .or. nvar.gt.0 )then
        if( nvar.gt.0 )then
          nread = 0
          if( iptra.eq.1 )then
            do n=1,min(nvar,npt)
              if( n.lt.10 )then
                text1 = 'ptX     '
                write(text1(3:3),161) n
161             format(i1.1)
              elseif( n.lt.100 )then
                text1 = 'ptXX    '
                write(text1(3:4),162) n
162             format(i2.2)
              else
                stop 11512
              endif
              call  readr(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),text1,           &
                          ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                          ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                          dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
              nread = nread+1
            enddo
          endif
          if( nread .lt. nvar )then
            ! need to read more data ....
            do n=nread+1,nvar
              if( n.lt.10 )then
                text1 = 'ptX     '
                write(text1(3:3),161) n
              elseif( n.lt.100 )then
                text1 = 'ptXX    '
                write(text1(3:4),162) n
              else
                stop 11513
              endif
              call  readr(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),text1,            &
                          ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                          ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                          dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
            enddo
          endif
        else
          if( myid.eq.0 ) print *
          if( myid.eq.0 ) print *,'  Note:  no passive tracer data in the restart file '
          if( myid.eq.0 ) print *
        endif
      endif

!---------------------------------------------------------------
!  time-average vars:

      IF( dotimeavg )THEN
      if( .not. interp_on_restart )then  ! cm1r21.1: do not read timeavg fields when interpolating to a new grid !

        if( myid.eq.0 )then
           read(50) ifoo
           read(50) ifoo
           read(50) ifoo
           read(50) ketaold
           read(50) ntavg
           read(50) rtavg
           print *,'    ketaold = ',ketaold
           print *,'    ntavg   = ',ntavg
           print *,'    rtavg   = ',rtavg
        endif
#ifdef MPI
        call MPI_BCAST(ntavg,ntim,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rtavg,ntim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
        do n=1,ntavr
          if(     n.eq.utav .or. n.eq.uutav )then
            call  readr(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'timavgu ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ixu,iy,iz,rxu,ry,rz,nxold+1,nyold,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
            do k=kbta,keta
            do j=1,nj
            do i=1,ni+1
              timavg(i,j,k,n) = dumu(i,j,k)
            enddo
            enddo
            enddo
          elseif( n.eq.vtav .or. n.eq.vvtav )then
            call  readr(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'timavgv ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iyv,iz,rx,ryv,rz,nxold,nyold+1,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
            do k=kbta,keta
            do j=1,nj+1
            do i=1,ni
              timavg(i,j,k,n) = dumv(i,j,k)
            enddo
            enddo
            enddo
          else
            call  readr(ni,nj,kbta,keta,nx,ny,timavg(ibta,jbta,kbta,n),'timavg  ',               &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          endif
        enddo
        do n=1,nsfctavr
          call  readr(ni,nj,1,1,nx,ny,sfctimavg(ibta,jbta,n),'sfctimav',                       &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
        do n=1,ntavr
        do nt=1,ntim
          if(     n.eq.utav .or. n.eq.uutav )then
            call  readr(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'tavgu   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ixu,iy,iz,rxu,ry,rz,nxold+1,nyold,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
            do k=kbta,keta
            do j=1,nj
            do i=1,ni+1
              tavg(i,j,k,nt,n) = dumu(i,j,k)
            enddo
            enddo
            enddo
          elseif( n.eq.vtav .or. n.eq.vvtav )then
            call  readr(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'tavgv   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iyv,iz,rx,ryv,rz,nxold,nyold+1,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
            do k=kbta,keta
            do j=1,nj+1
            do i=1,ni
              tavg(i,j,k,nt,n) = dumv(i,j,k)
            enddo
            enddo
            enddo
          else
            call  readr(ni,nj,kbta,keta,nx,ny,dum1(ib,jb,kbta),'tavg    ',                       &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
            do k=kbta,keta
            do j=1,nj
            do i=1,ni
              tavg(i,j,k,nt,n) = dum1(i,j,k)
            enddo
            enddo
            enddo
          endif
        enddo
        enddo
        do n=1,nsfctavr
        do nt=1,ntim
          call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'sfctavg ',                                &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,ketaold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          do j=1,nj
          do i=1,ni
            sfctavg(i,j,nt,n) = dum1(i,j,1)
          enddo
          enddo
        enddo
        enddo

      endif
      ENDIF

!---------------------------------------------------------------
!  parcels:

      if( myid.eq.0 )then
        if( restart_format.eq.1 )then
          read(50) nvar
#ifdef NETCDF
        elseif( restart_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
          call disp_err( nf90_get_var(ncid,varid,nvar,(/time_index/)) , .true. )
          if( iprcl.eq.0 ) nvar = 0
#endif
        endif
        print *,'  nvar_parcels = ',nvar
      endif

#ifdef MPI
      call MPI_BCAST(nvar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      if( iprcl.eq.1 .or. nvar.gt.0 )then
        if( nvar.gt.0 )then
          ! only read position info:
          if( myid.eq.0 )then
            IF( nvar.eq.nparcels )THEN
              ! easy:  restart file matches current config
              if( restart_format.eq.1 )then
                read(50) ploc
#ifdef NETCDF
              elseif( restart_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
                n = 3
                call disp_err( nf90_get_var(ncid,varid,ploc,(/1,1,time_index/),(/nparcels,n,1/)) , .true. )
#endif
              endif
            ELSE
              ! annoying:  restart file has different nparcels than current config
              IF( iprcl.eq.1 )THEN
                if( .not. terrain_flag )then
                  do np=1,nparcels
                    ploc(np,1) = pdata(np,prx)
                    ploc(np,2) = pdata(np,pry)
                    ploc(np,3) = pdata(np,prz)
                  enddo
                else
                  do np=1,nparcels
                    ploc(np,1) = pdata(np,prx)
                    ploc(np,2) = pdata(np,pry)
                    ploc(np,3) = pdata(np,prsig)
                  enddo
                endif
              ENDIF
              if( myid.eq.0 ) print *,'  start pfoo ' 
              allocate( pfoo(nvar,3) )
              if( restart_format.eq.1 )then
                read(50) pfoo
#ifdef NETCDF
              elseif( restart_format.eq.2 )then
                call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
                n = 3
                call disp_err( nf90_get_var(ncid,varid,pfoo,(/1,1,time_index/),(/nvar,n,1/)) , .true. )
#endif
              endif
              IF( iprcl.eq.1 )THEN
                do n=1,3
                do np=1,min(nvar,nparcels)
                  ploc(np,n) = pfoo(np,n)
                enddo
                enddo
              ENDIF
              deallocate( pfoo )
              if( myid.eq.0 ) print *,'  end pfoo ' 
            ENDIF
          endif
          IF( iprcl.eq.1 )THEN
#ifdef MPI
            call MPI_BCAST(ploc,3*nparcels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
            if( .not. terrain_flag )then
              DO np=1,nparcels
                pdata(np,prx)=ploc(np,1)
                pdata(np,pry)=ploc(np,2)
                pdata(np,prz)=ploc(np,3)
              ENDDO
            else
              DO np=1,nparcels
                pdata(np,prx)=ploc(np,1)
                pdata(np,pry)=ploc(np,2)
                pdata(np,prsig)=ploc(np,3)
              ENDDO
            endif
            restart_prcl = .true.
          ENDIF
        else
          if( myid.eq.0 ) print *
          if( myid.eq.0 ) print *,'  Note:  no parcel data in the restart file '
          if( myid.eq.0 ) print *
          restart_prcl = .false.
        endif
      endif

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  open bc:

      if(irbc.eq.4)then
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        if(myid.eq.0)then
          ndum = ny
        else
          ndum = 1
        endif
        allocate( dumy(ndum) )
        !----------------------
      if( wbc.eq.2 )then
        aname = 'radbcw'
#ifdef NETCDF
        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumy
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
          call readrbcwe(radbcw,aname,ndum,dumy,ibw,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      if( ebc.eq.2 )then
        aname = 'radbce'
#ifdef NETCDF
        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumy
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
          call readrbcwe(radbce,aname,ndum,dumy,ibe,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
        deallocate( dumy )
        if(myid.eq.0)then
          ndum = nx
        else
          ndum = 1
        endif
        allocate( dumx(ndum) )
        !----------------------
      if( sbc.eq.2 )then
        aname = 'radbcs'
#ifdef NETCDF
        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumx
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
          call readrbcsn(radbcs,aname,ndum,dumx,ibs,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
        enddo
      endif
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      if( nbc.eq.2 )then
        aname = 'radbcn'
#ifdef NETCDF
        if( restart_format.eq.2 .and. myid.eq.0 )then
          ncstatus = nf90_inq_varid(ncid,aname,varid)
          if(ncstatus.ne.nf90_noerr)then
            print *,'  Error1 in readrbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              read(50) dumx
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_get_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in readrbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
          call readrbcsn(radbcn,aname,ndum,dumx,ibn,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
        enddo
      endif
        !----------------------
        deallocate( dumx )
        !----------------------
        !cccccccccccccccccccccc
        !----------------------
      endif

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------
!  151001:  use theta  (over-rides perturbation value that was read-in above)


    IF( restart_use_theta )THEN
      call  readr(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'theta   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        tha(i,j,k) = dum1(i,j,k)-th0(i,j,k)
      enddo
      enddo
      enddo
    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

    IF( do_lsnudge .or. do_adapt_move )THEN
      call  readr(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ixu,iy,iz,rxu,ry,rz,nxold+1,nyold,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
      call  readr(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),datk1(0,0),datk2(0,0),nx,ny,nz,ix,iyv,iz,rx,ryv,rz,nxold,nyold+1,nzold,sigma(0),sigma0(0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
    ENDIF

!---------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!---------------------------------------------------------------

    IF( restart_format.eq.1 )THEN
      IF(myid.eq.0) close(unit=50)
      IF(myid.eq.nodeleader) close(unit=51)
      IF(myid.eq.nodeleader) close(unit=52)
      IF(myid.eq.nodeleader) close(unit=53)
      IF(myid.eq.nodeleader) close(unit=54)
#ifdef NETCDF
    ELSEIF( restart_format.eq.2 )THEN
      if( myid.eq.0 )then
        call disp_err( nf90_close(ncid) , .true. )
      endif
#endif
    ENDIF

    if( restarted ) nrst = nrst+1

!---------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  From restart file: '
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mtime   = ',mtime
      if(dowr) write(outfile,*) '   stattim = ',stattim
      if(dowr) write(outfile,*) '   taptim  = ',taptim
      if(dowr) write(outfile,*) '   rsttim  = ',rsttim
      if(dowr) write(outfile,*) '   radtim  = ',radtim
      if(dowr) write(outfile,*) '   prcltim = ',prcltim
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   nstep   = ',nstep
      if(dowr) write(outfile,*) '   srec    = ',srec
      if(dowr) write(outfile,*) '   sirec   = ',sirec
      if(dowr) write(outfile,*) '   urec    = ',urec
      if(dowr) write(outfile,*) '   vrec    = ',vrec
      if(dowr) write(outfile,*) '   wrec    = ',wrec
      if(dowr) write(outfile,*) '   nrec    = ',nrec
      if(dowr) write(outfile,*) '   mrec    = ',mrec
      if(dowr) write(outfile,*) '   prec    = ',prec
      if(dowr) write(outfile,*) '   nwrite  = ',nwrite
      if(dowr) write(outfile,*) '   nrst    = ',nrst
      if(dowr) write(outfile,*) '   nstatout= ',nstatout
      if( dodomaindiag )then
      if(dowr) write(outfile,*) '   trecs   = ',trecs
      if(dowr) write(outfile,*) '   trecw   = ',trecw
      if(dowr) write(outfile,*) '   nwritet = ',nwritet
      endif
      if( doazimavg )then
      if(dowr) write(outfile,*) '   arecs   = ',arecs
      if(dowr) write(outfile,*) '   arecw   = ',arecw
      if(dowr) write(outfile,*) '   nwritea = ',nwritea
      endif
      if( dohifrq )then
      if(dowr) write(outfile,*) '   nwriteh = ',nwriteh
      endif
      if( do_adapt_move )then
      if(dowr) write(outfile,*) '   nwritemv      = ',nwritemv
      if(dowr) write(outfile,*) '   mvrec         = ',mvrec
      if(dowr) write(outfile,*) '   umove         = ',umove
      if(dowr) write(outfile,*) '   vmove         = ',vmove
      if(dowr) write(outfile,*) '   adaptmovetime = ',adaptmovetim
      endif
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mass1 = ',mass1
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   dt      = ',dt
      if(dowr) write(outfile,*) '   dtlast  = ',dtlast
      if(dowr) write(outfile,*) '   dbldt   = ',dbldt
      if(dowr) write(outfile,*) '   adt     = ',adt
      if(dowr) write(outfile,*) '   adtlast = ',adtlast

!---------

      if( adapt_dt.eq.0 )then
        dt = dtl
        dbldt = dtl
      endif

      ! this is needed for stats files:
      nrec=nrec-1
      nstatout=nstatout-1

      IF( output_format .ne. old_format )THEN
        srec = 1
        sirec = 1
        urec = 1
        vrec = 1
        wrec = 1
        nrec = 1
        mrec = 1
        nwrite = 1
        prec = 1
      ENDIF

!---------

        IF( (imoist.eq.1).and.(ptype.eq.2) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call consat2(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF
        IF( (imoist.eq.1).and.(ptype.eq.4) )then
          if(timestats.ge.1) time_misc=time_misc+mytime()
          call lfoice_init(dt)
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF

!---------

#ifdef MPI
      if(timestats.ge.1)then
        ! this is needed for proper accounting of timing:
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
#endif

      if(allocated(xhref0))  deallocate( xhref0 )
      if(allocated(xfref0))  deallocate( xfref0 )
      if(allocated(yhref0))  deallocate( yhref0 )
      if(allocated(yfref0))  deallocate( yfref0 )
      if(allocated(sigma0))  deallocate( sigma0 )
      if(allocated(sigmaf0)) deallocate( sigmaf0 )
      if(allocated(datk1))   deallocate( datk1 )
      if(allocated(datk2))   deallocate( datk2 )
      if(allocated(ix))      deallocate( ix )
      if(allocated(iy))      deallocate( iy )
      if(allocated(iz))      deallocate( iz )
      if(allocated(ixu))     deallocate( ixu )
      if(allocated(iyv))     deallocate( iyv )
      if(allocated(izw))     deallocate( izw )
      if(allocated(rx))      deallocate( rx )
      if(allocated(ry))      deallocate( ry )
      if(allocated(rz))      deallocate( rz )
      if(allocated(rxu))     deallocate( rxu )
      if(allocated(ryv))     deallocate( ryv )
      if(allocated(rzw))     deallocate( rzw )

      if( interp_on_restart )then
        if(myid.eq.0) print *,'  setting qpten to zero '
        qpten = 0.0
        qtten = 0.0
        qvten = 0.0
        qcten = 0.0
      endif

      return

778   print *,'  error opening restart file '
      print *,'    ... stopping cm1 ... '
      call stopcm1

7001  print *
      print *,'  7001: error opening restart files '
      print *
      print *,'    Error attempting to read cm1rversion '
      print *
      print *,'  Note that cm1r21.0 and later versions require restart files '
      print *,'  written by cm1r21.0 or later. '
      print *
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

      end subroutine restart_read


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine  readr(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,datk1,datk2,nx,ny,nzr,ix,iy,iz,rx,ry,rz,nxold,nyold,nzold,sigma,sigma0,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
#ifdef MPI
    use mpi
#endif
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine reads restart files and then passes data 
    ! to other processors if this is a MPI run. 
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr,nxold,nyold,nzold,nx,ny,nzr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format,restart_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(0:nxold+1,0:nyold+1) :: datk1,datk2
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodeleader,nodes,nfile
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    real, intent(in), dimension(0:numk2) :: sigma
    real, intent(in), dimension(0:nzold) :: sigma0
    integer, intent(in), dimension(numi) :: ix
    integer, intent(in), dimension(numj) :: iy
    integer, intent(in), dimension(nzr)  :: iz
    real,    intent(in), dimension(numi) :: rx
    real,    intent(in), dimension(numj) :: ry
    real,    intent(in), dimension(nzr)  :: rz

    integer :: msk

!-------------------------------------------------------------------------------

    rf1:  IF( restart_filetype.eq.1 .or. restart_filetype.eq.2 )THEN

    if(myid.eq.0) print *,aname,numk1,numk2

    msk = 0

#ifndef MPI
        call   readr_sp(numk1,numk2,nxold,nyold,nxr,nyr,d2i,d2j,numi,numj,ngxy,dat2,var,orec,restart_format,nfile,ncid,aname,time_index)
#else
    IF( .not. dointerp )THEN
      ! no interpolation:
      IF(myid.ne.nodeleader)THEN
        call   readr_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
      ELSE
        IF(myid.ne.msk)THEN
          call readr_comm2(ppnode,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat3,var)
        ELSE
          call readr_comm3(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat2,dat3,var,  &
                           nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_format,nfile,orec,  &
                           msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2,  &
                           ncid,time_index,aname)
        ENDIF
      ENDIF
    ELSE
      ! with interpolation:
      call readr_interp(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat2,dat3,var,  &
                        nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_format,nfile,orec,  &
                        msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2)
    ENDIF
#endif


    ENDIF  rf1

!-------------------------------------------------------------------------------

#ifdef MPI
    rf2:  IF( restart_filetype.eq.3 )THEN

      call     readr2(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1(1,1),dat2(1,1),dat3(1,1,0),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)

    ENDIF  rf2
#endif

!-------------------------------------------------------------------------------
!ccccc  done  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-------------------------------------------------------------------------------

!!!#ifdef MPI
!!!    ! helps with memory:
!!!    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!    !----------------- end MPI section -----------------!
!!!#endif

    return
    end subroutine  readr


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifndef MPI


    ! single processor runs:
    subroutine readr_sp(numk1,numk2,nxold,nyold,nxr,nyr,d2i,d2j,numi,numj,ngxy,dat2,var,orec,restart_format,nfile,ncid,aname,time_index)
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: numk1,numk2,nxold,nyold,nxr,nyr,d2i,d2j,numi,numj,ngxy
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: restart_format,nfile
    character(len=8), intent(in) :: aname
    integer, intent(in) :: time_index

    integer :: i,j,k

#ifdef NETCDF
    integer :: varid,status
#endif

#ifdef NETCDF
    if( restart_format.eq.2 )then
      status = nf90_inq_varid(ncid,aname,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in  readr, aname = ',aname
        print *,nf90_strerror(status)
        call stopcm1
      endif
    endif
#endif


    kloop:  &
    DO k=numk1,numk2

      IF( restart_format.eq.1 )THEN
        ! ..... binary format .....
          ! not doing interp:
          call read_binary(nfile,nxr,nyr,1,d2i,1,d2j,dat2,orec)
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        ! ..... netcdf format .....
        if(numk1.eq.numk2)then
          status = nf90_get_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
        else
          status = nf90_get_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
        endif
        if(status.ne.nf90_noerr)then
          print *,'  Error2 in  readr, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
#endif
      ENDIF


      if( dointerp )then
  stop 43987
      endif


      ! all done ... put contents of dat2 into var array and move on !

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        var(i,j,k)=dat2(i,j)
      enddo
      enddo

    ENDDO  kloop



    end subroutine readr_sp


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! ordinary processor ... send data to nodeleader:
    subroutine readr_comm1(numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader,dat1,var)
    use mpi
    implicit none

    integer, intent(in) :: numk1,numk2,d3i,d3j,ngxy,numi,numj,nodeleader
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

    integer :: i,j,k,tag,reqs,ierr

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      ! recv data from nodeleader:
      call MPI_IRECV(dat1(1,1),d3i*d3j,MPI_REAL,nodeleader,tag,MPI_COMM_WORLD,reqs,ierr)
      call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        var(i,j,k)=dat1(i,j)
      enddo
      enddo
      ! DONE

      tag = tag+2

    ENDDO  kloop

    end subroutine readr_comm1
    ! done, ordinary processor


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader but not proc 0:
    subroutine readr_comm2(ppnode,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid,reqt,dat3,var)
    use mpi
    implicit none

    integer, intent(in) :: ppnode,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,msk,myid
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

    integer :: i,j,k,tag,reqs,ierr,proc
    integer, dimension(mpi_status_size,ppnode-1) :: status1

    tag = 1

    kloop:  &
    DO k=numk1,numk2

      ! get data from msk:
      call MPI_IRECV(dat3(1,1,myid),d3i*d3j*ppnode,MPI_REAL,msk,tag+1,MPI_COMM_WORLD,reqs,ierr)
      ! wait for data to arrive:
      call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
      ! start sends to other processors on a node:
      do proc=myid+1,myid+(ppnode-1)
        call MPI_ISEND(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
      enddo
      ! my data:
      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,numj
      do i=1,numi
        var(i,j,k)=dat3(i,j,myid)
      enddo
      enddo
      ! wait for sends to finish:
      call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status1,ierr)
      ! DONE

      tag = tag+2

    ENDDO  kloop

    end subroutine readr_comm2
    ! done, nodeleader but not proc 0:


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader and proc 0:
    subroutine readr_comm3(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat2,dat3,var,  &
                           nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_format,nfile,orec,  &
                           msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2,  &
                           ncid,time_index,aname)
    use mpi
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: nxold,nyold,nzold,nxr,nyr,nzr
    integer, intent(in), dimension(nxr) :: ix
    integer, intent(in), dimension(nyr) :: iy
    integer, intent(in), dimension(nzr) :: iz
    real,    intent(in), dimension(nxr) :: rx
    real,    intent(in), dimension(nyr) :: ry
    real,    intent(in), dimension(nzr) :: rz
    integer, intent(in) :: restart_format,nfile
    integer, intent(inout) :: orec
    integer, intent(in) :: msk,myid,numprocs,mynode,nodes,ni,nj
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    real, intent(inout), dimension(0:nxold+1,0:nyold+1) :: datk1,datk2
    integer, intent(in) :: ncid,time_index
    character(len=8), intent(in) :: aname

    integer :: i,j,k,tag,reqs,ierr,kold1,kold2
    integer :: fooi,fooj,n,nn,nnn,ntot,n1,n2,nitmp,njtmp,proc,index,index2
    integer, dimension(mpi_status_size,ppnode-1) :: status1
    integer, dimension(mpi_status_size,ppnode-1 + nodes-1) :: status2
    logical :: alldone
#ifdef NETCDF
    integer :: varid,status
#endif

#ifdef NETCDF
    if( restart_format.eq.2 )then
      status = nf90_inq_varid(ncid,aname,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in  readr, aname = ',aname
        print *,nf90_strerror(status)
        call stopcm1
      endif
    endif
#endif

    tag = 1
    kold1 = -1
    kold2 =  0
    alldone = .false.

    kloop:  &
    DO k=numk1,numk2

          ! read data:
      IF( restart_format.eq.1 )THEN
        ! ..... binary format .....
        if( dointerp )then
          ! doing interp:
          call   prep_for_interp(k,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid,iz,kold1,kold2,orec,alldone,datk1,datk2)
        else
          ! not doing interp:
          call read_binary(nfile,nxr,nyr,1,d2i,1,d2j,dat2,orec)
        endif
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        ! ..... netcdf format .....
        if(numk1.eq.numk2)then
          status = nf90_get_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
        else
          status = nf90_get_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
        endif
        if(status.ne.nf90_noerr)then
          print *,'  Error2 in  readr, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
#endif
      ENDIF


      if( dointerp )then
        stop 3333
      endif



          ! send data:
          do nn=1,( nodes-1 )
              ! send data to other nodeleaders:
              index2 = nn
              if( index2.le.mynode )then
                index2 = index2-1
              endif
              n1 = index2*ppnode
              n2 = (index2+1)*ppnode-1
              do nnn=n1,n2
                proc = nnn
                fooi = myi1p(proc+1)-1
                fooj = myj1p(proc+1)-1
                nitmp = myi2p(proc+1)-myi1p(proc+1)+1
                njtmp = myj2p(proc+1)-myj1p(proc+1)+1
                if( numi.gt.ni ) nitmp = nitmp+1
                if( numj.gt.nj ) njtmp = njtmp+1
                !$omp parallel do default(shared)  &
                !$omp private(i,j)
                do j=1,njtmp
                do i=1,nitmp
                  dat3(i,j,proc) = dat2(fooi+i,fooj+j)
                enddo
                enddo
              enddo
              proc = index2*ppnode
              call MPI_ISEND(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+1,MPI_COMM_WORLD,reqt(ppnode-1+nn),ierr)
          enddo
          do nn=1,( ppnode-1 )
              ! send data to ordinary procs on this node:
              proc = myid+nn
              fooi = myi1p(proc+1)-1
              fooj = myj1p(proc+1)-1
              nitmp = myi2p(proc+1)-myi1p(proc+1)+1
              njtmp = myj2p(proc+1)-myj1p(proc+1)+1
              if( numi.gt.ni ) nitmp = nitmp+1
              if( numj.gt.nj ) njtmp = njtmp+1
              !$omp parallel do default(shared)  &
              !$omp private(i,j)
              do j=1,njtmp
              do i=1,nitmp
                dat3(i,j,proc) = dat2(fooi+i,fooj+j)
              enddo
              enddo
              call MPI_ISEND(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(nn),ierr)
          enddo
          ! my data:
          if( myid.eq.0 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j)
            do j=1,numj
            do i=1,numi
              var(i,j,k) = dat2(i,j)
            enddo
            enddo
          else
            proc = myid
            fooi = myi1p(proc+1)-1
            fooj = myj1p(proc+1)-1
            nitmp = myi2p(proc+1)-myi1p(proc+1)+1
            njtmp = myj2p(proc+1)-myj1p(proc+1)+1
            if( numi.gt.ni ) nitmp = nitmp+1
            if( numj.gt.nj ) njtmp = njtmp+1
            !$omp parallel do default(shared)  &
            !$omp private(i,j)
            do j=1,njtmp
            do i=1,nitmp
              var(i,j,k) = dat2(fooi+i,fooj+j)
            enddo
            enddo
          endif
          ntot = ppnode-1 + nodes-1
          call mpi_waitall(ntot,reqt(1:ntot),status2,ierr)


      tag = tag+2

    ENDDO  kloop

    end subroutine readr_comm3
    ! done nodeleader and proc 0


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    subroutine readr_interp(ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2,reqt,dat2,dat3,var,  &
                           nxold,nyold,nzold,nxr,nyr,nzr,ix,iy,iz,rx,ry,rz,restart_format,nfile,orec,  &
                           msk,myid,numprocs,mynode,nodes,ni,nj,myi1p,myi2p,myj1p,myj2p,datk1,datk2)
    use mpi
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    integer, intent(in) :: ppnode,d2i,d2j,d3i,d3j,d3n,d3t,ngxy,numi,numj,numk1,numk2
    integer, intent(inout), dimension(d3t) :: reqt
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: nxold,nyold,nzold,nxr,nyr,nzr
    integer, intent(in), dimension(numi) :: ix
    integer, intent(in), dimension(numj) :: iy
    integer, intent(in), dimension(nzr)  :: iz
    real,    intent(in), dimension(numi) :: rx
    real,    intent(in), dimension(numj) :: ry
    real,    intent(in), dimension(nzr)  :: rz
    integer, intent(in) :: restart_format,nfile
    integer, intent(inout) :: orec
    integer, intent(in) :: msk,myid,numprocs,mynode,nodes,ni,nj
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    real, intent(inout), dimension(0:nxold+1,0:nyold+1) :: datk1,datk2

    integer :: i,j,k,tag,reqs,ierr,kold1,kold2
    integer :: fooi,fooj,n,nn,nnn,ntot,n1,n2,nitmp,njtmp,proc,index,index2
    integer, dimension(mpi_status_size,ppnode-1) :: status1
    integer, dimension(mpi_status_size,ppnode-1 + nodes-1) :: status2
    logical :: alldone
#ifdef NETCDF
    integer :: varid,status
#endif

#ifdef NETCDF
!    if( myid.eq.0 )then
!    if( restart_format.eq.2 )then
!      status = nf90_inq_varid(ncid,aname,varid)
!      if(status.ne.nf90_noerr)then
!        print *,'  Error1 in  readr, aname = ',aname
!        print *,nf90_strerror(status)
!        call stopcm1
!      endif
!    endif
!    endif
#endif

    tag = 1
    kold1 = -1
    kold2 =  0
    alldone = .false.

    kloop:  &
    DO k=numk1,numk2

    if(myid.eq.msk)then
!!!      print *,'  k = ',k
          ! read data:
      IF( restart_format.eq.1 )THEN
        ! ..... binary format .....
          call   prep_for_interp(k,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid,iz,kold1,kold2,orec,alldone,datk1,datk2)
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        ! ..... netcdf format .....
        stop 12871
#endif
      ENDIF
    endif

        call MPI_BCAST(datk1,(nxold+2)*(nyold+2),MPI_REAL,msk,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(datk2,(nxold+2)*(nyold+2),MPI_REAL,msk,MPI_COMM_WORLD,ierr)

      if( dointerp )then
        call     do_interp(k,myid,nxr,nyr,nzr,nxold,nyold,d2i,d2j,ix,iy,iz,rx,ry,rz,datk1,datk2,numi,numj,numk1,numk2,ngxy,var)
      endif


    ENDDO  kloop

    end subroutine readr_interp


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine read_binary(nfile,nxr,nyr,i1,i2,j1,j2,datk2,orec)
      implicit none

      integer, intent(in) :: nfile,nxr,nyr,i1,i2,j1,j2
      real, intent(inout) :: datk2(i1:i2,j1:j2)
      integer, intent(inout) :: orec

      integer :: i,j

        read(nfile,rec=orec) ((datk2(i,j),i=1,nxr),j=1,nyr)
        orec = orec+1

      end subroutine read_binary


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine setbc(nxold,nyold,datk)
    implicit none

    integer, intent(in) :: nxold,nyold
    real, intent(inout), dimension(0:nxold+1,0:nyold+1) :: datk

    integer :: i,j

    do j=1,nyold
      datk(0,j) = datk(1,j)
      datk(nxold+1,j) = datk(nxold,j)
    enddo
    do i=0,nxold+1
      datk(i,0) = datk(i,1)
      datk(i,nyold+1) = datk(i,nyold)
    enddo

    end subroutine setbc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine prep_for_interp(k,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid,iz,kold1,kold2,orec,alldone,datk1,datk2)
      implicit none

      integer, intent(in) :: k,nzr,nxold,nyold,nzold,numk1,numk2,nfile,myid
      integer, intent(in), dimension(nzr) :: iz
      integer, intent(inout) :: kold1,kold2,orec
      logical, intent(inout) :: alldone
      real, intent(inout), dimension(0:nxold+1,0:nyold+1) :: datk1,datk2

      integer :: i,j,kk,kk1,kk2


    !   Layout:
    !
    !    old grid     new grid     old grid
    !      kold1                   kold2
    !        iz          k          iz+1
    !

!!!            print *,'    iz(k)+1,kold2 = ',iz(k)+1,kold2

            ! check to see if we need to increment kold:
            neednew:  &
            IF( iz(k)+1.gt.kold2 )THEN

              ! yes, we need a new level

              ! is more data available?
              inside_old_domain:  &
              IF( (iz(k)+1).le.nzold )THEN

                  ! yes, more levels are available ... so, get to work:

                  ! determine how many levels to read from the old grid:
                  if( numk1.eq.numk2 )then
                    ! 2D variable:
                    kk1 = 1
                    kk2 = 1
                  else
                    ! 3D variable:
                    kk1 = kold2+1
                    kk2 = iz(k)+1
                  endif

!!!                  print *,'    kk1,kk2 = ',kk1,kk2

                  ! read through old levels:
                  do kk=kk1,kk2

                    kold1 = kold1+1
                    kold2 = kold2+1

                    ! save old level:

                    do j=0,nyold+1
                    do i=0,nxold+1
                      datk1(i,j) = datk2(i,j)
                    enddo
                    enddo

!!!                    print *,'      ...... reading kold ...... ',kk
                    call read_binary(nfile,nxold,nyold,0,nxold+1,0,nyold+1,datk2,orec)
                  enddo

                  call setbc(nxold,nyold,datk2)

                  if( kold1.eq.0 )then
                    ! we just read the first level ... extrapolate to lower levels

!!!                    print *,'    below first level: '
                    ! (NOTE: need to replace with extrapolation)

                    do j=0,nyold+1
                    do i=0,nxold+1
                      datk1(i,j) = datk2(i,j)
                    enddo
                    enddo

                  endif


              ELSE  inside_old_domain

                  ! we have already read all levels from old grid:

                  if( .not. alldone )then
                    ! First time here ... equate the two old levels (for now)

                    ! (NOTE: need to replace with extrapolation)

!!!                    print *,'    above last level: '
                    do j=0,nyold+1
                    do i=0,nxold+1
                      datk1(i,j) = datk2(i,j)
                    enddo
                    enddo
                  endif

                  alldone = .true.

              ENDIF  inside_old_domain

            ENDIF  neednew

      end subroutine prep_for_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine do_interp(k,myid,nxr,nyr,nzr,nxold,nyold,d2i,d2j,ix,iy,iz,rx,ry,rz,datk1,datk2,numi,numj,numk1,numk2,ngxy,var)
      implicit none

      integer, intent(in) :: k,myid,nxr,nyr,nzr,nxold,nyold,d2i,d2j
      integer, intent(in), dimension(numi) :: ix
      integer, intent(in), dimension(numj) :: iy
      integer, intent(in), dimension(nzr)  :: iz
      real, intent(in), dimension(numi) :: rx
      real, intent(in), dimension(numj) :: ry
      real, intent(in), dimension(nzr)  :: rz
      real, intent(in), dimension(0:nxold+1,0:nyold+1) :: datk1,datk2
      integer, intent(in) :: numi,numj,numk1,numk2,ngxy
      real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

      double precision :: w1,w2,w3,w4,w5,w6,w7,w8,wsum

      integer :: i,j,ii,jj

        ! this is where the interpolation happens:

        do j=1,numj
        do i=1,numi

          w1=(1.0-rx(i))*(1.0-ry(j))*(1.0-rz(k))
          w2=rx(i)*(1.0-ry(j))*(1.0-rz(k))
          w3=(1.0-rx(i))*ry(j)*(1.0-rz(k))
          w4=(1.0-rx(i))*(1.0-ry(j))*rz(k)
          w5=rx(i)*(1.0-ry(j))*rz(k)
          w6=(1.0-rx(i))*ry(j)*rz(k)
          w7=rx(i)*ry(j)*(1.0-rz(k))
          w8=rx(i)*ry(j)*rz(k)
          wsum = w1+w2+w3+w4+w5+w6+w7+w8

!          ! debug !
!          if( rx(i).lt.-0.0001 .or. rx(i).gt.1.0001 .or.  &
!              ry(j).lt.-0.0001 .or. ry(j).gt.1.0001 .or.  &
!              rz(k).lt.-0.0001 .or. rz(k).gt.1.0001 .or.  &
!              wsum.le.0.99999 .or.                  &
!              wsum.ge.1.00001 )then
!            print *,'  245987 '
!            call stopcm1
!          endif
!          ! debug !

          var(i,j,k)=( datk1(ix(i)  ,iy(j)  )*w1    &
                      +datk1(ix(i)+1,iy(j)  )*w2    &
                      +datk1(ix(i)  ,iy(j)+1)*w3    &
                      +datk2(ix(i)  ,iy(j)  )*w4    &
                      +datk2(ix(i)+1,iy(j)  )*w5    &
                      +datk2(ix(i)  ,iy(j)+1)*w6    &
                      +datk1(ix(i)+1,iy(j)+1)*w7    &
                      +datk2(ix(i)+1,iy(j)+1)*w8  )

        enddo
        enddo

      end subroutine do_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! cm1r17-format restart files !
    subroutine  readr2(numi,numj,numk1,numk2,nxr,nyr,var,aname,          &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
    use mpi
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine reads restart files and then passes data 
    ! to other processors if this is a MPI run. 
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format,restart_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d3i*ppnode,d3j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodeleader,nodes,nfile
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,msk
    integer :: reqs,index,index2,n,nn,nnn,fooi,fooj,proc,ierr,ntot,n1,n2
    integer :: tag
    integer, dimension(mpi_status_size,ppnode-1) :: status1

    DO k=numk1,numk2
      IF(myid.ne.nodeleader)THEN
        call MPI_IRECV(dat1,d3i*d3j,MPI_REAL,nodeleader,k,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,mpi_status_ignore,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          var(i,j,k) = dat1(i,j)
        enddo
        enddo
      ELSE
        read(50) dat2
        do proc=myid+1,myid+(ppnode-1)
          fooi = numi*(proc-myid)
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,numj
          do i=1,numi
            dat3(i,j,proc)=dat2(fooi+i,j)
          enddo
          enddo
          call MPI_ISEND(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,k,MPI_COMM_WORLD,reqt(proc-myid),ierr)
        enddo
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          var(i,j,k)=dat2(i,j)
        enddo
        enddo
        call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status1,ierr)
      ENDIF
    ENDDO

    return
    end subroutine  readr2
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE restart_read_module
