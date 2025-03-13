  MODULE restart_write_module

  implicit none

  private
  public :: restart_write,readrbcwe,readrbcsn

  CONTAINS

      subroutine restart_write(nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,      &
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
                               dum1,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p)
          ! end_restart_write
      use input
      use constants
      use lsnudge_module
#ifdef MPI
      use mpi
#endif
#ifdef NETCDF
      use netcdf
      use writeout_nc_module, only : restart_prelim,disp_err
#endif
      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the writing of restart files
      !----------------------------------------------------------

      integer, intent(in) :: nstep,srec,sirec,urec,vrec,wrec,nrec,mrec,prec,trecs,trecw,arecs,arecw
      integer, intent(in) :: nwrite,nwritet,nwritea,nwriteh,nrst,nstatout
      integer, intent(in) :: num_soil_layers,nrad2d
      double precision, intent(in)    :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(in) :: dt,dtlast
      integer, intent(in) :: ndt
      double precision, intent(in) :: adt,adtlast,acfl,dbldt
      double precision, intent(in) :: mass1
      double precision, intent(in) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      character(len=3), intent(in), dimension(maxq) :: qname
      real, dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
      real, dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0,prs0,pi0,rho0,qv0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain,raincv,sws,svs,sps,srs,sgs,sus,shs
      real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,s1,t1,xland,psfc,thflux,qvflux,prate,ustt,ut,vt,st
      real, intent(in), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(in), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(in),    dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: dumu
      real, intent(in),    dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: dumv
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,thraten,effc,effi,effs,effr,effg,effis
      real, intent(in), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(in), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(in), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx,    &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,fm,fh,mznt,swspd,wstar,delta
      real, intent(in), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(in), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(in), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppx
      real, intent(in),    dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, intent(in)   , dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl, &
           qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,  &
           sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D
      real, intent(in)   , dimension(ibmynn:iemynn,jbmynn:jemynn) :: vdfg,maxmf
      integer, intent(in)   , dimension(ibmynn:iemynn,jbmynn:jemynn) :: nupdraft,ktop_plume
      real, intent(in),    dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj,el_myj
      real, intent(in),    dimension(ibmyj:iemyj,jbmyj:jemyj) :: mixht,akhs,akms,elflx,ct,snow,sice,thz0,qz0,uz0,vz0,th10,q10,z0base,zntmyj
      integer, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: lowlyr,ivgtyp
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, intent(in   ) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(in   ) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      integer, intent(in   ) :: icenter,jcenter
      real, intent(in   ) :: xcenter,ycenter
      double precision, intent(in) :: domainlocx,domainlocy
      double precision, intent(in   ) :: adaptmovetim
      integer, intent(in   ) :: mvrec,nwritemv
      real, intent(in),    dimension(kb:ke) :: ug,vg
      real, intent(in),    dimension(kb:ke) :: gamk,kmw,ufw,vfw,u1b,v1b
      integer, intent(in   ), dimension(ntim) :: ntavg
      double precision, intent(in   ), dimension(ntim) :: rtavg
      double precision, intent(in   ), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntim,ntavr) :: tavg
      real, intent(in   ), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      double precision, intent(in   ), dimension(ibta:ieta,jbta:jeta,ntim,nsfctavr) :: sfctavg
      real, intent(in   ), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg
      real, intent(inout), dimension(ib:ie,jb:je) :: dumsfc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

      character(len=maxstring) :: fname
      character(len=8) :: text1
      character(len=6) :: aname
      integer :: i,j,k,n,np,nt,nvar,reqs,orecs,orecu,orecv,orecw,ndum
      integer :: ncid,time_index
      real, dimension(:), allocatable :: dumx,dumy
#ifdef MPI
      integer :: proc,index,count,req1,req2,req3,reqp
      double precision, dimension(nbudget) :: cfoo
      double precision, dimension(numq) :: afoo,bfoo
#endif
#ifdef NETCDF
      integer :: varid,ncstatus
#endif

!-----------------------------------------------------------------------

  IF( restart_format.eq.1 )THEN
    ! unformatted direct-access (grads) format:

  IF( restart_filetype.eq.1 )THEN

    !------------------
    ! one restart file (per stagger type):
    IF(myid.eq.nodeleader)THEN
      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_x.dat'

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='unknown')

      fname = 'cm1rst_s.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecs = 1

      fname = 'cm1rst_u.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny)
      orecu = 1

      fname = 'cm1rst_v.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1))
      orecv = 1


      fname = 'cm1rst_w.dat'
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.2 )THEN

    !------------------
    ! one restart file (per restart time):
    IF(myid.eq.nodeleader)THEN
      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_x.dat'
      write(fname( 8:13),101) nrst
101   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname

      if( myid.eq.0 )  &
      open(unit=50,file=fname,form='unformatted',status='unknown')

      fname = 'cm1rst_XXXXXX_s.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=51,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecs = 1

      fname = 'cm1rst_XXXXXX_u.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=52,file=fname,form='unformatted',access='direct',recl=4*(nx+1)*ny)
      orecu = 1

      fname = 'cm1rst_XXXXXX_v.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=53,file=fname,form='unformatted',access='direct',recl=4*nx*(ny+1))
      orecv = 1

      fname = 'cm1rst_XXXXXX_w.dat'
      write(fname( 8:13),101) nrst
      if(dowr) write(outfile,*) '  fname=',fname
      open(unit=54,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
      orecw = 1

      if(dowr) write(outfile,*)
    ENDIF

  ELSEIF( restart_filetype.eq.3 )THEN

    !------------------
    ! one restart file per node (cm1r17 format):
    IF(myid.eq.nodeleader)THEN

      do i=1,maxstring
        fname(i:i) = ' '
      enddo

      fname = 'cm1rst_XXXXXX_YYYYYY.dat'

      write(fname( 8:13),102) mynode
      write(fname(15:20),102) nrst
102   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname
      if(dowr) write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='unknown')
    ENDIF
  ELSE
    stop 12388
  ENDIF

#ifdef NETCDF
  ELSEIF( restart_format.eq.2 )THEN
    ! netcdf format:

    if( myid.eq.0 )then
      call     restart_prelim(nrst,ncid,mtime,xfref,yfref,zh,zf,sigma,sigmaf,  &
                              qname,num_soil_layers,nrad2d,dat2(1,1),dat2(1,2),dum1(ib,jb,kb),time_index)
    endif

#endif
  ELSE

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

  ENDIF

!---------------------------------------------------------------
! metadata:

  IF( restart_format.eq.1 )THEN
    IF(myid.eq.0)THEN
      ! only processor 0 does this:
      write(50) cm1rversion
      write(50) nx
      write(50) ny
      write(50) nz
      write(50) ngxy
      write(50) wbc
      write(50) ebc
      write(50) sbc
      write(50) nbc
      write(50) bbc
      write(50) tbc










      write(50) (maxx-minx)
      write(50) (maxy-miny)
      write(50) maxz
















      write(50) nstep
      write(50) srec
      write(50) sirec
      write(50) urec
      write(50) vrec
      write(50) wrec
      write(50) nrec
      write(50) mrec
      write(50) prec
      write(50) trecs
      write(50) trecw
      write(50) arecs
      write(50) arecw
      write(50) mvrec
      write(50) nwrite
      write(50) nwritet
      write(50) nwritea
      write(50) nwritemv
      write(50) nwriteh
      write(50) nrst
      write(50) nstatout
      write(50) ndt
      write(50) icenter
      write(50) jcenter
      write(50) output_format
      write(50) dt
      write(50) dtlast
      write(50) xcenter
      write(50) ycenter
      write(50) umove
      write(50) umovem1
      write(50) umovem2
      write(50) vmove
      write(50) vmovem1
      write(50) vmovem2
      write(50) domainlocx
      write(50) domainlocy
      write(50) adaptmovetim
      write(50) cflmax
      write(50) mtime
      write(50) stattim
      write(50) taptim
      write(50) rsttim
      write(50) radtim
      write(50) prcltim
      write(50) adt
      write(50) adtlast
      write(50) acfl
      write(50) dbldt
      write(50) mass1
      write(50) avgsfcu
      write(50) avgsfcv
      write(50) avgsfcs
      write(50) avgsfcsu
      write(50) avgsfcsv
      write(50) avgsfct
      write(50) avgsfcq
      write(50) avgsfcp
      write(50) xhref
      write(50) xfref
      write(50) yhref
      write(50) yfref
      write(50) sigma
      write(50) sigmaf
    ENDIF
#ifdef NETCDF
  ELSEIF( restart_format.eq.2 )THEN

    IF(myid.eq.0)THEN

      call disp_err( nf90_inq_varid(ncid,"cm1rversion",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,cm1rversion,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstep",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nstep,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"srec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,srec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"sirec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,sirec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"urec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,urec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"wrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,wrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,prec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,trecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"trecw",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,trecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,arecs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"arecw",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,arecw,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mvrec",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mvrec,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwrite",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwrite,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritet",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritet,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritea",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritea,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwritemv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwritemv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nwriteh",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nwriteh,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nrst",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nrst,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"nstatout",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,nstatout,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ndt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,ndt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"icenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,icenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"jcenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,jcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"old_format",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,output_format,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dtlast",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"xcenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,xcenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"ycenter",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,ycenter,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umove",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,umove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umovem1",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,umovem1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"umovem2",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,umovem2,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmove",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vmove,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmovem1",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vmovem1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"vmovem2",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,vmovem2,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocx",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,domainlocx,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"domainlocy",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,domainlocy,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adaptmovetim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adaptmovetim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"cflmax",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,cflmax,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mtime",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"stattim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,stattim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"taptim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,taptim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"rsttim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,rsttim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"radtim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,radtim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"prcltim",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,prcltim,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"adtlast",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,adtlast,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"acfl",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,acfl,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"dbldt",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,dbldt,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"mass1",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,mass1,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcu",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcs",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcs,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsu",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcsu,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcsv",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcsv,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfct",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfct,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcq",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcq,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,"avgsfcp",varid) , .true. )
      call disp_err( nf90_put_var(ncid,varid,avgsfcp,(/time_index/)) , .true. )

    ENDIF
#endif
  ENDIF

!-----------------------------------------------------------------------
#ifdef MPI
      cfoo = 0.0
      call MPI_REDUCE(qbudget(1),cfoo(1),nbudget,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                      MPI_COMM_WORLD,ierr)
      if( myid.eq.0 )then
        do n=1,nbudget
          qbudget(n)=cfoo(n)
        enddo
      else
        qbudget = 0.0
      endif
      if( imoist.eq.1 )then
        afoo = 0.0
        call MPI_REDUCE(asq(1),afoo(1),numq,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                        MPI_COMM_WORLD,ierr)
        if( myid.eq.0 )then
          do n=1,numq
            asq(n)=afoo(n)
          enddo
        else
          asq = 0.0
        endif
        bfoo = 0.0
        call MPI_REDUCE(bsq(1),bfoo(1),numq,MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                        MPI_COMM_WORLD,ierr)
        if( myid.eq.0 )then
          do n=1,numq
            bsq(n)=bfoo(n)
          enddo
        else
          bsq = 0.0
        endif
      endif
#endif
!-----------------------------------------------------------------------
! budget variables:

    IF( myid.eq.0 )THEN

      IF( restart_format.eq.1 )THEN
        write(50) qbudget
        write(50) asq
        write(50) bsq
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"qbudget",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,qbudget,(/1,time_index/),(/nbudget,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"asq",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,asq,(/1,time_index/),(/numq,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"bsq",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,bsq,(/1,time_index/),(/numq,1/)) , .true. )
#endif
      ENDIF

    ENDIF

!---------------------------------------------------------------

    if( myid.eq.0 )then
      IF( restart_format.eq.1 )THEN
        write(50) ug
        write(50) vg
#ifdef NETCDF
      ELSEIF( restart_format.eq.2 )THEN
        call disp_err( nf90_inq_varid(ncid,"ug",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,ug(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
        call disp_err( nf90_inq_varid(ncid,"vg",varid) , .true. )
        call disp_err( nf90_put_var(ncid,varid,vg(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
#endif
      ENDIF
    endif

!---------------------------------------------------------------
! standard 2D:

      n = 1
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs     ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    if( nrain.eq.2 )then
      n = 2
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),'rain2   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),'sws2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),'svs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),'sps2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),'srs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),'sgs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),'sus2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),'shs2    ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
      call writer(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

!---------------------------------------------------------------
! standard 3D:

      call writer(ni,nj,1,nk,nx,ny,rho(ib,jb,1),'rho     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,prs(ib,jb,1),'prs     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'ua      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
      call writer(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'va      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
      call writer(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'wa      ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),'ppi     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,tha(ib,jb,1),'tha     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ppx(ib,jb,1),'ppx     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    if( psolver.eq.6 )then
      call writer(ni,nj,1,nk,nx,ny,phi1(ib,jb,1),'phi1    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,phi2(ib,jb,1),'phi2    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
    IF(imoist.eq.1)THEN
    do n=1,numq
      text1 = '        '
      write(text1(1:3),156) qname(n)
156   format(a3)
      call writer(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),text1     ,       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call writer(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),'qpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),'qtten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),'qvten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),'qcten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif
    if( idoles .and. iusetke )then
      call writer(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tkea    ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    endif

!---------------------------------------------------------------
!  radiation:

      if(radopt.ge.1)then
        call writer(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,swten(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,thraten(ib,jb,1),'thraten ',     &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radsw(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rnflx(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radswnet(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radlwin(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        do n=1,nrad2d
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rad2d(i,j,n)
        enddo
        enddo
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
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),text1,             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
      endif

      if( radopt.ge.1 .and. ptype.eq.5 )then
        call writer(ni,nj,1,nk,nx,ny,effc(ib,jb,1),'effc    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effi(ib,jb,1),'effi    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effs(ib,jb,1),'effs    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effr(ib,jb,1),'effr    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effg(ib,jb,1),'effg    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,nk,nx,ny,effis(ib,jb,1),'effis   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      if((oceanmodel.eq.2).or.(ipbl.ge.1).or.(sfcmodel.ge.1))then
        !---- (1) ----!
      if(sfcmodel.ge.1)then
        call writer(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,u1(ib,jb),'u1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,v1(ib,jb),'v1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,s1(ib,jb),'s1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,t1(ib,jb),'t1      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,s10(ib,jb),'s10     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psfc(ib,jb),'psfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( tbc.eq.3 )then
        call writer(ni,nj,1,1,nx,ny,ustt(ib,jb),'ustt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  ut(ib,jb),'ut      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  vt(ib,jb),'vt      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,  st(ib,jb),'st      ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      endif
      endif


      if(sfcmodel.ge.1)then
        !---- (2) ----!
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = lu_index(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu_index',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = kpbl2d(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'kpbl2d  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,wspd(ib,jb),'wspd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psim(ib,jb),'psim    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,psih(ib,jb),'psih    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),'gz1oz0  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CHS(ib,jb),'chs     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CHS2(ib,jb),'chs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CQS2(ib,jb),'cqs2    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CPMM(ib,jb),'cpmm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ZOL(ib,jb),'zol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),'mavail  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,MOL(ib,jb),'mol     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,RMOL(ib,jb),'rmol    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,REGIME(ib,jb),'regime  ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,LH(ib,jb),'lh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,FLHC(ib,jb),'flhc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,FLQC(ib,jb),'flqc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,QGH(ib,jb),'qgh     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CK(ib,jb),'ck      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CKA(ib,jb),'cka     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,CDA(ib,jb),'cda     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,USTM(ib,jb),'ustm    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,QSFC(ib,jb),'qsfc    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,T2(ib,jb),'t2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,Q2(ib,jb),'q2      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,TH2(ib,jb),'th2     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,EMISS(ib,jb),'emiss   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,THC(ib,jb),'thc     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,ALBD(ib,jb),'albd    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,chklowq(ib,jb),'chklowq ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,capg(ib,jb),'capg    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,snowc(ib,jb),'snowc   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,fm(ib,jb),'fm      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,fh(ib,jb),'fh      ',            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,mznt(ib,jb),'mznt    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,swspd(ib,jb),'swspd   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,wstar(ib,jb),'wstar   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,delta(ib,jb),'delta   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
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
          call writer(ni,nj,1,1,nx,ny,tslb(ib,jb,n),text1,             &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
      endif
      endif

      if(oceanmodel.eq.2)then
        call writer(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,t0ml(ib,jb),'t0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,h0ml(ib,jb),'h0ml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        call writer(ni,nj,1,1,nx,ny,tmoml(ib,jb),'tmoml   ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                    ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      endif

!---------------------------------------------------------------

    IF( ipbl.eq.4 .or. ipbl.eq.5 .or. sfcmodel.eq.6 )THEN
      if(myid.eq.0) print *,'  writing mynn vars ... '
      call writer(ni,nj,1,nk,nx,ny,qke(ib,jb,1),'qke     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qke3d(ib,jb,1),'qke3d   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,tsq(ib,jb,1),'tsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qsq(ib,jb,1),'qsq     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,cov(ib,jb,1),'cov     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sh3d(ib,jb,1),'sh3d    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sm3d(ib,jb,1),'sm3d    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,el_pbl(ib,jb,1),'el_pbl  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qc_bl(ib,jb,1),'qc_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qi_bl(ib,jb,1),'qi_bl   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,cldfra_bl(ib,jb,1),'cldfra_b',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_a(ib,jb,1),'edmf_a  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_w(ib,jb,1),'edmf_w  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_qt(ib,jb,1),'edmf_qt ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_thl(ib,jb,1),'edmf_thl',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_ent(ib,jb,1),'edmf_ent',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,edmf_qc(ib,jb,1),'edmf_qc ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sub_thl3d(ib,jb,1),'sub_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,sub_sqv3d(ib,jb,1),'sub_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,det_thl3d(ib,jb,1),'det_thl3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,det_sqv3d(ib,jb,1),'det_sqv3',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,thpten(ib,jb,1),'thpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvpten(ib,jb,1),'qvpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qcpten(ib,jb,1),'qcpten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qipten(ib,jb,1),'qipten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,upten(ib,jb,1),'upten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,vpten(ib,jb,1),'vpten   ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qnipten(ib,jb,1),'qnipten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qncpten(ib,jb,1),'qncpten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,vdfg(ib,jb),'vdfg    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,maxmf(ib,jb),'maxmf   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = nupdraft(i,j)
        dum1(i,j,2) = ktop_plume(i,j)
      enddo
      enddo
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'nupdraft',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ktop_plu',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)






      if(myid.eq.0) print *,'  ... done writing mynn vars '
    ENDIF

!---------------------------------------------------------------

    IF(ipbl.eq.6)THEN
      if(myid.eq.0) print *,'  writing myj vars ... '
      call writer(ni,nj,1,nk,nx,ny,tke_myj(ib,jb,1),'tke_myj ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,el_myj(ib,jb,1),'el_myj  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF(sfcmodel.eq.7)THEN
      call writer(ni,nj,1,1,nx,ny,mixht(ib,jb),'mixht   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,akhs(ib,jb),'akhs    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,akms(ib,jb),'akms    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,elflx(ib,jb),'elflx   ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,ct(ib,jb),'ct      ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,snow(ib,jb),'snow    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,sice(ib,jb),'sice    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,thz0(ib,jb),'thz0    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,qz0(ib,jb),'qz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,uz0(ib,jb),'uz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,vz0(ib,jb),'vz0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,th10(ib,jb),'th10    ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,q10(ib,jb),'q10     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,z0base(ib,jb),'z0base  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,zntmyj(ib,jb),'zntmyj  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = lowlyr(i,j)
        dum1(i,j,2) = ivgtyp(i,j)
      enddo
      enddo
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lowlyr  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'ivgtyp  ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)






      if(myid.eq.0) print *,'  ... done writing myj vars '
    endif

!---------------------------------------------------------------
!  convection param.

    IF( cuparam.ge.1 )THEN
      call writer(ni,nj,1,1,nx,ny,raincv(ib,jb,1),'raincv  ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if(nrain.eq.2)  &
      call writer(ni,nj,1,1,nx,ny,raincv(ib,jb,2),'raincv2 ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,thcuten(ib,jb,1),'thcuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qvcuten(ib,jb,1),'qvcuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qccuten(ib,jb,1),'qccuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,qicuten(ib,jb,1),'qicuten ',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,ucuten(ib,jb,1),'ucuten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      call writer(ni,nj,1,nk,nx,ny,vcuten(ib,jb,1),'vcuten  ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF

!---------------------------------------------------------------
!  2-part turbulence param:

      if( myid.eq.0 )then
        if( restart_format.eq.1 )then
          write(50) gamk
          write(50) kmw
          write(50) ufw
          write(50) vfw
          write(50) u1b
          write(50) v1b
#ifdef NETCDF
        elseif( restart_format.eq.2 )then
          call disp_err( nf90_inq_varid(ncid,"gamk",varid) , .true. )
          call disp_err( nf90_put_var(ncid,varid,gamk(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
          call disp_err( nf90_inq_varid(ncid,"kmw",varid) , .true. )
          call disp_err( nf90_put_var(ncid,varid,kmw(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
          call disp_err( nf90_inq_varid(ncid,"ufw",varid) , .true. )
          call disp_err( nf90_put_var(ncid,varid,ufw(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
          call disp_err( nf90_inq_varid(ncid,"vfw",varid) , .true. )
          call disp_err( nf90_put_var(ncid,varid,vfw(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
          call disp_err( nf90_inq_varid(ncid,"u1b",varid) , .true. )
          call disp_err( nf90_put_var(ncid,varid,u1b(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
          call disp_err( nf90_inq_varid(ncid,"v1b",varid) , .true. )
          call disp_err( nf90_put_var(ncid,varid,v1b(1:nk),(/1,time_index/),(/nk,1/)) , .true. )
#endif
        endif
      endif

!---------------------------------------------------------------
!  passive tracers:

      if(iptra.eq.1)then

        if(myid.eq.0)then
          if( restart_format.eq.1 )then
            write(50) npt
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,npt,(/time_index/)) , .true. )
#endif
          endif
        endif
        do n=1,npt
          if( n.lt.10 )then
            text1 = 'ptX     '
            write(text1(3:3),161) n
161         format(i1.1)
          elseif( n.lt.100 )then
            text1 = 'ptXX    '
            write(text1(3:4),162) n
162         format(i2.2)
          else
            stop 11512
          endif
          call writer(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),text1,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
      else
        if(myid.eq.0)then
          nvar = 0
          if( restart_format.eq.1 )then
            write(50) nvar
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"npt",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nvar,(/time_index/)) , .true. )
#endif
          endif
        endif
      endif

!---------------------------------------------------------------
!  time-average vars:

      IF( dotimeavg )THEN

        if( myid.eq.0 )then
          write(50) ntim
          write(50) ntavr
          write(50) nsfctavr
          write(50) keta
          write(50) ntavg
          write(50) rtavg
        endif
        do n=1,ntavr
          if(     n.eq.utav .or. n.eq.uutav )then
            do k=kbta,keta
            do j=1,nj
            do i=1,ni+1
              dumu(i,j,k) = timavg(i,j,k,n)
            enddo
            enddo
            enddo
            call writer(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'timavgu ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
          elseif( n.eq.vtav .or. n.eq.vvtav )then
            do k=kbta,keta
            do j=1,nj+1
            do i=1,ni
              dumv(i,j,k) = timavg(i,j,k,n)
            enddo
            enddo
            enddo
            call writer(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'timavgv ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
          else
            call writer(ni,nj,kbta,keta,nx,ny,timavg(ibta,jbta,kbta,n),'timavg  ',               &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          endif
        enddo
        do n=1,nsfctavr
          call writer(ni,nj,1,1,nx,ny,sfctimavg(ibta,jbta,n),'sfctimav',                       &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
        do n=1,ntavr
        do nt=1,ntim
          if(     n.eq.utav .or. n.eq.uutav )then
            do k=kbta,keta
            do j=1,nj
            do i=1,ni+1
              dumu(i,j,k) = tavg(i,j,k,nt,n)
            enddo
            enddo
            enddo
            call writer(ni+1,nj,kbta,keta,nx+1,ny,dumu(ib,jb,1),'tavgu   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
          elseif( n.eq.vtav .or. n.eq.vvtav )then
            do k=kbta,keta
            do j=1,nj+1
            do i=1,ni
              dumv(i,j,k) = tavg(i,j,k,nt,n)
            enddo
            enddo
            enddo
            call writer(ni,nj+1,kbta,keta,nx,ny+1,dumv(ib,jb,1),'tavgv   ',     &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
          else
            do k=kbta,keta
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = tavg(i,j,k,nt,n)
            enddo
            enddo
            enddo
            call writer(ni,nj,kbta,keta,nx,ny,dum1(ib,jb,kbta),'tavg    ',                       &
                        ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                        ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                        dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
          endif
        enddo
        enddo
        do n=1,nsfctavr
        do nt=1,ntim
          do j=1,nj
          do i=1,ni
            dum1(i,j,1) = sfctavg(i,j,nt,n)
          enddo
          enddo
          call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'sfctavg ',                                &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,                           &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
        enddo
        enddo

      ENDIF

!---------------------------------------------------------------
!  parcels:

      if(iprcl.eq.1)then
        !-----------------
        ! with parcels:
        if(myid.eq.0)then
          if( restart_format.eq.1 )then
            write(50) nparcels
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nparcels,(/time_index/)) , .true. )
#endif
          endif
        endif
        ! only write position info:
        if(myid.eq.0)then
          if( .not. terrain_flag )then
            DO np=1,nparcels
              ploc(np,1)=pdata(np,prx)
              ploc(np,2)=pdata(np,pry)
              ploc(np,3)=pdata(np,prz)
            ENDDO
          else
            DO np=1,nparcels
              ploc(np,1)=pdata(np,prx)
              ploc(np,2)=pdata(np,pry)
              ploc(np,3)=pdata(np,prsig)
            ENDDO
          endif
          if( restart_format.eq.1 )then
            write(50) ploc
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"ploc",varid) , .true. )
            n = 3
            call disp_err( nf90_put_var(ncid,varid,ploc,(/1,1,time_index/),(/nparcels,n,1/)) , .true. )
#endif
          endif
        endif
      else
        !-----------------
        ! without parcels:
        if(myid.eq.0)then
          nvar = 0
          if( restart_format.eq.1 )then
            write(50) nvar
#ifdef NETCDF
          elseif( restart_format.eq.2 )then
            call disp_err( nf90_inq_varid(ncid,"numparcels",varid) , .true. )
            call disp_err( nf90_put_var(ncid,varid,nvar,(/time_index/)) , .true. )
#endif
          endif
        endif
        !-----------------
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
            print *,'  Error1 in writerbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          call writerbcwe(radbcw,aname,ndum,dumy,ibw,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumy
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
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
            print *,'  Error1 in writerbcwe, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          call writerbcwe(radbce,aname,ndum,dumy,ibe,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumy
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumy,(/1,k,time_index/),(/ny,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcwe, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
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
            print *,'  Error1 in writerbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          call writerbcsn(radbcs,aname,ndum,dumx,ibs,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumx
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
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
            print *,'  Error1 in writerbcsn, aname = ',aname
            print *,nf90_strerror(ncstatus)
            call stopcm1
          endif
        endif
#endif
        do k=1,nk
          call writerbcsn(radbcn,aname,ndum,dumx,ibn,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
          if( myid.eq.0 )then
            if( restart_format.eq.1 )then
              write(50) dumx
#ifdef NETCDF
            elseif( restart_format.eq.2 )then
              ncstatus = nf90_put_var(ncid,varid,dumx,(/1,k,time_index/),(/nx,1,1/))
              if(ncstatus.ne.nf90_noerr)then
                print *,'  Error2 in writerbcsn, aname = ',aname
                print *,nf90_strerror(ncstatus)
                call stopcm1
              endif
#endif
            endif
          endif
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
!  150820:  optional variables


    IF( restart_file_theta )THEN
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = th0(i,j,k)+tha(i,j,k)
      enddo
      enddo
      enddo
      call writer(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'theta   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_u0 .or. do_lsnudge .or. do_adapt_move )THEN
      call writer(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecu,52,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
    ENDIF
    IF( restart_file_v0 .or. do_lsnudge .or. do_adapt_move )THEN
      call writer(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecv,53,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
    ENDIF
    IF( restart_file_dbz .and. qd_dbz.ge.1 )THEN
      ! cm1r19:
      call writer(ni,nj,1,nk,nx,ny,qdiag(ibdq,jbdq,1,qd_dbz),'dbz     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    !-----
    IF( restart_file_th0 )THEN
      call writer(ni,nj,1,nk,nx,ny,th0(ib,jb,1),'th0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_prs0 )THEN
      call writer(ni,nj,1,nk,nx,ny,prs0(ib,jb,1),'prs0    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_pi0 )THEN
      call writer(ni,nj,1,nk,nx,ny,pi0(ib,jb,1),'pi0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_rho0 )THEN
      call writer(ni,nj,1,nk,nx,ny,rho0(ib,jb,1),'rho0    ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_qv0 )THEN
      call writer(ni,nj,1,nk,nx,ny,qv0(ib,jb,1),'qv0     ',        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    !-----
    IF( restart_file_zs )THEN
      call writer(ni,nj,1,1,nx,ny,zs(ib,jb),'zs      ',            &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_zh )THEN
      call writer(ni,nj,1,nk,nx,ny,zh(ib,jb,1),'zhalf   ',         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF
    IF( restart_file_zf )THEN
      call writer(ni,nj,1,nk+1,nx,ny,zf(ib,jb,1),'zfull   ',       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecw,54,   &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p, &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
    ENDIF

    IF( restart_file_diags )THEN
      if( td_diss.gt.0 )                                                 &
      call writer(ni,nj,1,nk,nx,ny,tdiag(ib,jb,1,td_diss),'dissheat',    &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( td_mp.gt.0 )                                                   &
      call writer(ni,nj,1,nk,nx,ny,tdiag(ib,jb,1,td_mp),'mptend  ',      &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( qd_vtc.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtc),'vtc     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( qd_vtr.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtr),'vtr     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( qd_vts.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vts),'vts     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( qd_vtg.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vtg),'vtg     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
      if( qd_vti.gt.0 )                                                  &
      call writer(ni,nj,1,nk,nx,ny,qdiag(ib,jb,1,qd_vti),'vti     ',     &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,orecs,51,         &
                  ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,       &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
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

#ifdef MPI
      if(timestats.ge.1)then
        ! this is needed for proper accounting of timing:
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
#endif

      return
      end subroutine restart_write


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  !  211215:  move restart_read to its own file


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writerbcwe(radbc,aname,ndum,dumy,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: ndum,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,k
      real, intent(in), dimension(jb:je,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumy
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myj1p

      integer :: j,j1,j2
#ifdef MPI
      integer :: fooj,proc,reqs,ierr,req1,njrecv
#endif

#ifndef MPI
      do j=1,ny
        dumy(j) = radbc(j,k)
      enddo
#else
      IF(myid.ne.0)THEN
        if( ibndy.eq.1 )then
          call MPI_ISEND(nj,1,MPI_INTEGER,0,30,MPI_COMM_WORLD,req1,ierr)
          call MPI_ISEND(radbc(1,k),nj,MPI_REAL,0,31,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(req1,mpi_status_ignore,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        endif
      ELSE
        if( (aname.eq.'radbcw') .or. (aname.eq.'radbce' .and. nodex.eq.1) )then
          do j=1,nj
            dumy(j) = radbc(j,k)
          enddo
          j1 = 2
          j2 = nodey
        else
          j1 = 1
          j2 = nodey
        endif
        do j=j1,j2
          if( aname.eq.'radbcw' )then
            proc = (j-1)*nodex
          else
            proc = (j-1)*nodex + (nodex-1)
          endif
          call MPI_IRECV(njrecv,1,MPI_INTEGER,proc,30,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
          fooj = myj1p(proc+1)
          call MPI_IRECV(dumy(fooj),njrecv,MPI_REAL,proc,31,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        enddo
      ENDIF
#endif

      end subroutine writerbcwe


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writerbcsn(radbc,aname,ndum,dumx,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: ndum,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,k
      real, intent(in), dimension(ib:ie,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumx
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myi1p

      integer :: i,i1,i2
#ifdef MPI
      integer :: fooi,proc,reqs,ierr,req1,nirecv
#endif

#ifndef MPI
      do i=1,nx
        dumx(i) = radbc(i,k)
      enddo
#else
      IF(myid.ne.0)THEN
        if( ibndy.eq.1 )then
          call MPI_ISEND(ni,1,MPI_INTEGER,0,33,MPI_COMM_WORLD,req1,ierr)
          call MPI_ISEND(radbc(1,k),ni,MPI_REAL,0,32,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(req1,mpi_status_ignore,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        endif
      ELSE
        if( (aname.eq.'radbcs') .or. (aname.eq.'radbcn' .and. nodey.eq.1) )then
          do i=1,ni
            dumx(i) = radbc(i,k)
          enddo
          i1 = 2
          i2 = nodex
        else
          i1 = 1
          i2 = nodex
        endif
        do i=i1,i2
          if( aname.eq.'radbcs' )then
            proc = (i-1)
          else
            proc = (i-1) + nodex*(nodey-1)
          endif
          call MPI_IRECV(nirecv,1,MPI_INTEGER,proc,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
          fooi = myi1p(proc+1)
          call MPI_IRECV(dumx(fooi),nirecv,MPI_REAL,proc,32,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        enddo
      ENDIF
#endif

      end subroutine writerbcsn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine  readrbcwe(radbc,aname,ndum,dumy,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myj1p)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: ndum,ibndy,jb,je,kb,ke,ny,ni,nj,nk,nodex,nodey,k
      real, intent(inout), dimension(jb:je,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumy
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myj1p

      integer :: j,j1,j2
#ifdef MPI
      integer :: fooi,fooj,proc,reqs,ierr,req1,njrecv
#endif

#ifndef MPI
      do j=1,ny
        radbc(j,k) = dumy(j)
      enddo
#else
      IF(myid.ne.0)THEN
        if( ibndy.eq.1 )then
          call MPI_ISEND(nj,1,MPI_INTEGER,0,30,MPI_COMM_WORLD,req1,ierr)
          call MPI_IRECV(radbc(1,k),nj,MPI_REAL,0,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(req1,mpi_status_ignore,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        endif
      ELSE
        if( (aname.eq.'radbcw') .or. (aname.eq.'radbce' .and. nodex.eq.1) )then
          do j=1,nj
            radbc(j,k) = dumy(j)
          enddo
          j1 = 2
          j2 = nodey
        else
          j1 = 1
          j2 = nodey
        endif
        do j=j1,j2
          if( aname.eq.'radbcw' )then
            proc = (j-1)*nodex
          else
            proc = (j-1)*nodex + (nodex-1)
          endif
          call MPI_IRECV(njrecv,1,MPI_INTEGER,proc,30,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
          call MPI_ISEND(dumy(myj1p(proc+1)),njrecv,MPI_REAL,proc,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        enddo
      ENDIF
#endif

      end subroutine  readrbcwe


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine  readrbcsn(radbc,aname,ndum,dumx,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,restart_format,myid,k,numprocs,myi1p)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(in) :: ndum,ibndy,ib,ie,kb,ke,nx,ni,nj,nk,nodex,nodey,k
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbc
      character(len=6), intent(in) :: aname
      real, intent(inout), dimension(ndum) :: dumx
      integer, intent(in) :: restart_format,myid
      integer, intent(in) :: numprocs
      integer, intent(in), dimension(numprocs) :: myi1p

      integer :: i,i1,i2
#ifdef MPI
      integer :: fooi,fooj,proc,reqs,ierr,req1,nirecv
#endif

#ifndef MPI
      do i=1,nx
        radbc(i,k) = dumx(i)
      enddo
#else
      IF(myid.ne.0)THEN
        if( ibndy.eq.1 )then
          call MPI_ISEND(ni,1,MPI_INTEGER,0,33,MPI_COMM_WORLD,req1,ierr)
          call MPI_IRECV(radbc(1,k),ni,MPI_REAL,0,34,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(req1,mpi_status_ignore,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        endif
      ELSE
        if( (aname.eq.'radbcs') .or. (aname.eq.'radbcn' .and. nodey.eq.1) )then
          do i=1,ni
            radbc(i,k) = dumx(i)
          enddo
          i1 = 2
          i2 = nodex
        else
          i1 = 1
          i2 = nodex
        endif
        do i=i1,i2
          if( aname.eq.'radbcs' )then
            proc = (i-1)
          else
            proc = (i-1) + nodex*(nodey-1)
          endif
          call MPI_IRECV(nirecv,1,MPI_INTEGER,proc,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
          call MPI_ISEND(dumx(myi1p(proc+1)),nirecv,MPI_REAL,proc,34,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,mpi_status_ignore,ierr)
        enddo
      ENDIF
#endif

      end subroutine  readrbcsn


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writer(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
#ifdef MPI
    use mpi
#endif
#ifdef NETCDF
    use netcdf
#endif
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing of restart files.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=8), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    integer, intent(inout) :: orec,ncid
    integer, intent(in) :: time_index,restart_format,restart_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodeleader,nodes,nfile
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,msk,nitmp,njtmp
#ifdef MPI
    integer :: reqs,index,index2,n,nn,nnn,fooi,fooj,proc,ierr,ntot,n1,n2,tag
    logical :: recv1,recv2
    integer, dimension(mpi_status_size,ppnode-1) :: status1
#endif
#ifdef NETCDF
    integer :: varid,status
#endif

!-------------------------------------------------------------------------------

    rf1:  IF( restart_filetype.eq.1 .or. restart_filetype.eq.2 )THEN

    if(myid.eq.0) print *,aname

  ! 220315: lets clean this code up, following what was done in writeout:

    msk = 0

#ifdef MPI
    !----------------- MPI section -----------------!
    recv1 = .true.
    recv2 = .true.
    tag = 1
#endif

    kloop:  DO k=numk1,numk2

#ifndef MPI
      !-------------------- non-MPI section --------------------!
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,numj
      do i=1,numi
        dat2(i,j)=var(i,j,k)
      enddo
      enddo
#ifdef NETCDF
      if( restart_format.eq.2 )then
        status = nf90_inq_varid(ncid,aname,varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error1 in writer, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
      endif
#endif
#else
      iamnodeleader:  IF(myid.ne.nodeleader)THEN
        ! ordinary processor ... send data to nodeleader:
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j)=var(i,j,k)
        enddo
        enddo
        call MPI_ISEND(dat1(1,1),d3i*d3j,MPI_REAL,nodeleader,tag,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
        ! DONE, ordinary processors
      ELSE
        ! begin nodeleader section:
        if( recv1 )then
          ! start receives from all other processors on a node:
          do proc=myid+1,myid+(ppnode-1)
            call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
          enddo
        endif
        iammsk:  IF(myid.ne.msk)THEN
          ! nodeleader, not proc msk:
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,numj
          do i=1,numi
            dat3(i,j,myid)=var(i,j,k)
          enddo
          enddo
          ! wait for receives to finish:
          call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status1,ierr)
          ! send data to processor msk:
          call MPI_ISEND(dat3(1,1,myid),d3i*d3j*ppnode,MPI_REAL,msk,tag+1,MPI_COMM_WORLD,reqs,ierr)
          ! wait for send to finish:
          call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
          recv1 = .true.
          ! DONE, nodeleader (not proc msk)
        ELSE
          ! proc msk:
          if( recv2 )then
            ! start receives from other nodeleaders:
            do n = 1,(nodes-1)
              if( n.le.mynode )then
                proc = (n-1)*ppnode
              else
                proc = n*ppnode
              endif
              call MPI_IRECV(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+1,MPI_COMM_WORLD,reqt(ppnode-1+n),ierr)
            enddo
          endif
#ifdef NETCDF
          if( restart_format.eq.2 .and. k.eq.numk1 )then
            status = nf90_inq_varid(ncid,aname,varid)
            if(status.ne.nf90_noerr)then
              print *,'  Error1 in writer, aname = ',aname
              print *,nf90_strerror(status)
              call stopcm1
            endif
          endif
#endif
          ! my data:
          if( myid.eq.0 )then
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,numj
            do i=1,numi
              dat2(i,j)=var(i,j,k)
            enddo
            enddo
          else
            fooj = myid / nodex + 1
            fooi = myid - (fooj-1)*nodex  + 1
            fooi = (fooi-1)*ni
            fooj = (fooj-1)*nj
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,numj
            do i=1,numi
              dat2(fooi+i,fooj+j)=var(i,j,k)
            enddo
            enddo
          endif
          ! wait for data to arrive:
          ntot = ppnode-1 + nodes-1
          do nn=1,ntot
            call mpi_waitany(ntot,reqt(1:ntot),index,MPI_STATUS_IGNORE,ierr)
            if( index.le.(ppnode-1) )then
              ! data from ordinary procs on node:
              proc = myid+index
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
                dat2(fooi+i,fooj+j) = dat3(i,j,proc)
              enddo
              enddo
            else
              ! data from other nodeleaders:
              index2 = index-(ppnode-1)
              if( index2.le.mynode )then
                index2 = index2-1
              endif
              n1 = index2*ppnode
              n2 = (index2+1)*ppnode-1
              do nnn = n1,n2
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
                  dat2(fooi+i,fooj+j) = dat3(i,j,proc)
                enddo
                enddo
              enddo
            endif
          enddo
          ! DONE, proc msk
          ! processor is ready to write.
          IF( k.lt.numk2 )THEN
            ! start receives for next level:
            do proc=myid+1,myid+(ppnode-1)
              call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag+2,MPI_COMM_WORLD,reqt(proc-myid),ierr)
            enddo
            recv1 = .false.
!!!#ifdef NETCDF
!!!            IF( restart_format.eq.2 )THEN
              do n = 1,(nodes-1)
                proc = n*ppnode
                call MPI_IRECV(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+3,MPI_COMM_WORLD,reqt(ppnode-1+n),ierr)
              enddo
              recv2 = .false.
!!!            ENDIF
!!!#endif
          ENDIF
        ENDIF  iammsk
      ENDIF  iamnodeleader
#endif

        ! WRITE DATA:
        IF( myid.eq.msk )THEN
          !---   write data   ------------------!
          IF( restart_format.eq.1 )THEN
            write(nfile,rec=orec) ((dat2(i,j),i=1,nxr),j=1,nyr)
#ifdef NETCDF
          ELSEIF( restart_format.eq.2 )THEN
            ! ----- netcdf format -----
            if(numk1.eq.numk2)then
              status = nf90_put_var(ncid,varid,dat2,(/1,1,time_index/),(/d2i,d2j,1/))
            else
              status = nf90_put_var(ncid,varid,dat2,(/1,1,k,time_index/),(/d2i,d2j,1,1/))
            endif
            if(status.ne.nf90_noerr)then
              print *,'  Error2 in writer, aname = ',aname
              print *,nf90_strerror(status)
              call stopcm1
            endif
#endif
          ENDIF
        ENDIF

      !---  prepare for next level   -------!
      IF( restart_format.eq.1 )THEN
        orec = orec+1
!!!#ifdef MPI
!!!        msk = msk+ppnode
!!!        if( msk.ge.numprocs ) msk = msk-numprocs
!!!#endif
      ENDIF
#ifdef MPI
      tag = tag+2
#endif
      !---  done with this level   ---------!
    ENDDO  kloop

    ENDIF  rf1

!-------------------------------------------------------------------------------

#ifdef MPI
    rf2:  IF( restart_filetype.eq.3 )THEN

      call    writer2(numi,numj,numk1,numk2,nxr,nyr,var,aname,           &
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
    end subroutine writer


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! ordinary processor ... send data to nodeleader:
    subroutine writer_comm1
    use mpi
    implicit none

    end subroutine writer_comm1
    ! done, ordinary processor


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader but not proc 0:
    subroutine writer_comm2
    use mpi
    implicit none

    end subroutine writer_comm2
    ! done, nodeleader but not proc 0:


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


    ! nodeleader and proc 0:
    subroutine writer_comm3
    use mpi
    implicit none


    end subroutine writer_comm3
    ! done nodeleader and proc 0


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI

    ! cm1r17-format restart files !
    subroutine writer2(numi,numj,numk1,numk2,nxr,nyr,var,aname,          &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,orec,nfile,   &
                      ncid,time_index,restart_format,restart_filetype,myi1p,myi2p,myj1p,myj2p,   &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
    use mpi
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing of restart files.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(in   ), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
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
    integer :: reqs,index,index2,n,nn,nnn,fooi,fooj,proc,ierr,ntot,n1,n2,tag
    logical :: recv1,recv2

    DO k=numk1,numk2
      IF(myid.ne.nodeleader)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j) = var(i,j,k)
        enddo
        enddo
        call MPI_ISEND(dat1,d3i*d3j,MPI_REAL,nodeleader,k,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,mpi_status_ignore,ierr)
      ELSE
        do proc=myid+1,myid+(ppnode-1)
          call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,k,MPI_COMM_WORLD,reqt(proc-myid),ierr)
        enddo
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat2(i,j)=var(i,j,k)
        enddo
        enddo
        nn = 1
        do while( nn.le.(ppnode-1) )
          nn = nn + 1
          call mpi_waitany(ppnode-1,reqt(1:ppnode-1),index,MPI_STATUS_IGNORE,ierr)
          fooi = numi*index
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,numj
          do i=1,numi
            dat2(fooi+i,j)=dat3(i,j,nodeleader+index)
          enddo
          enddo
        enddo
        write(50) dat2
      ENDIF
    ENDDO

    return
    end subroutine writer2


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE restart_write_module
