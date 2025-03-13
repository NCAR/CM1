  MODULE writeout_module

  implicit none

  private
  public :: writeout

  integer :: nout
  logical :: opens,openu,openv,openw,dointerp
  integer, parameter :: varmax = 10000

  CONTAINS


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


             subroutine writeout(srec,urec,vrec,wrec,mrec,rtime,dt,fnum,nwrite,qname,qunit,    &
                        nstep,mtime,adt,ndt,domainlocx,domainlocy,                             &
                        name_output,desc_output,unit_output,grid_output,levs_output,cmpr_output, &
                        xh,xf,uf,yh,yf,vf,xfref,yfref,                                         &
                        rds,sigma,rdsf,sigmaf,zh,zf,mf,gx,gy,wprof,                            &
                        pi0,prs0,rho0,rr0,rf0,rrf0,th0,qv0,u0,v0,thv0,rth0,qc0,qi0,            &
                        zs,rgzu,rgzv,rain,raincv,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,psfc,      &
                        rxh,arh1,arh2,uh,ruh,rxf,arf1,arf2,vh,rvh,mh,rmh,rmf,rr,rf,            &
                        gz,rgz,gzu,gzv,gxu,gyv,dzdx,dzdy,c1,c2,                                &
                        cd,ch,cq,tlh,f2d,psmth,prate,ustt,cm0,                                 &
                        dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,                          &
                        t11,t12,t13,t22,t23,t33,rho,prs,divx,                                  &
                        rru,ua ,dumu,ugr  ,rrv,va ,dumv,vgr  ,rrw,wa ,dumw,ppi ,tha ,phi2,     &
                        sadv,thten,nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2,               &
                        thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,xkzh,xkzq,xkzm, &
                        lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,  &
                        qa ,p3a,p3o,kmh,kmv,khh,khv,cme,tkea ,swten,lwten,cldfra,o30,          &
                        radsw,rnflx,radswnet,radlwin,dsr,olr,pta,                              &
                        effc,effi,effs,effr,effg,effis,                                        &
                        lwupt,lwdnt,lwupb,lwdnb,                                               &
                        swupt,swdnt,swupb,swdnb,lwcf,swcf,coszen,                              &
                        num_soil_layers,u10,v10,s10,t2,q2,znt,ust,stau,tst,qst,z0t,z0q,u1,v1,s1,     &
                        hpbl,zol,mol,rmol,br,brcr,wscale,wscaleu,phim,phih,psim,psih,psiq,wspd,qsfc,wstar,delta,prkpp,fm,fh, &
                        mznt,taux,tauy,                                                        &
                        z0base,thz0,qz0,uz0,vz0,u10e,v10e,tke_myj,el_myj,  &
                        tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl,                    &
                        qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,                               &
                        edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,                       &
                        vdfg,maxmf,nupdraft,ktop_plume,                                        &
                        dat1,dat2,dat3,reqt,dum2d          ,dumk1,dumk2,                       &
                        tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out2d,out3d,cir,crr,cangle,        &
                        bndy,kbdy,hflxw,hflxe,hflxs,hflxn,recy_cap,recy_inj,timavg,sfctimavg,kmwk, &
                        thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,                         &
                        nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,myi1p,myi2p,myj1p,myj2p)
                        ! end_writeout
      use input
      use constants
      use bc_module
      use comm_module
      use misclibs
      use getcape_module
      use cm1libs , only : rslf,rsif
      use ib_module
      use sfcphys_module , only : stabil_funcs
      use eddy_recycle
#ifdef MPI
      use mpi
#endif
#ifdef NETCDF
      use netcdf
      use writeout_nc_module, only : disp_err,netcdf_prelim
#endif
      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes writeouts for GrADS-format and
      ! netcdf-format output.
      !----------------------------------------------------------

      integer, intent(inout) :: srec,urec,vrec,wrec,mrec
      real, intent(in) :: rtime,dt
      integer, intent(in) :: fnum,nwrite
      character(len=3), dimension(maxq), intent(in) :: qname
      character(len=20), intent(in), dimension(maxq) :: qunit
      integer, intent(in) :: nstep,ndt
      double precision, intent(in) :: mtime,adt,domainlocx,domainlocy
      character(len=60), intent(inout), dimension(maxvars) :: desc_output
      character(len=40), intent(inout), dimension(maxvars) :: name_output,unit_output
      character(len=1),  intent(inout), dimension(maxvars) :: grid_output
      integer, intent(inout), dimension(maxvars) :: levs_output
      logical, intent(inout), dimension(maxvars) :: cmpr_output
      real, dimension(ib:ie), intent(in) :: xh
      real, dimension(ib:ie+1), intent(in) :: xf,uf
      real, dimension(jb:je), intent(in) :: yh
      real, dimension(jb:je+1), intent(in) :: yf,vf
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, dimension(kb:ke), intent(in) :: rds,sigma
      real, dimension(kb:ke+1), intent(in) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf,mf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gy
      real, intent(in), dimension(kb:ke) :: wprof
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: pi0,prs0,rho0,rr0,rf0,rrf0,th0,qv0,thv0,rth0,qc0,qi0
      real, dimension(ib:ie,jb:je), intent(in) :: zs
      real, dimension(itb:ite,jtb:jte), intent(in) :: rgzu,rgzv
      real, dimension(ib:ie,jb:je,nrain), intent(in) :: rain,raincv,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je), intent(in) :: xland,psfc,psmth,thflux,qvflux,cd,ch,cq,tlh,f2d,prate,ustt,cm0
      real, intent(in), dimension(ib:ie) :: rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: rxf,arf1,arf2
      real, intent(in), dimension(jb:je) :: vh,rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rr,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu,gyv
      real, intent(in), dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, dimension(ib:ie,jb:je,kb:ke), intent(inout) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: rho,prs
      real, dimension(ib:ie,jb:je,kb:ke), intent(inout) :: divx
      real, dimension(ib:ie+1,jb:je,kb:ke), intent(in) :: u0,ua
      real, dimension(ib:ie+1,jb:je,kb:ke), intent(inout) :: rru,dumu,ugr
      real, dimension(ib:ie,jb:je+1,kb:ke), intent(in) :: v0,va
      real, dimension(ib:ie,jb:je+1,kb:ke), intent(inout) :: rrv,dumv,vgr
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: wa
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(inout) :: rrw,dumw
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: ppi,tha
      real, intent(in), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi2
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: sadv,thten
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2
      real, dimension(ibb:ieb,jbb:jeb,kbb:keb), intent(in) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, dimension(ibb:ieb,jbb:jeb,kbb:keb), intent(in) :: xkzh,xkzq,xkzm
      integer, dimension(ibl:iel,jbl:jel), intent(in) :: lu_index
      real, dimension(ib:ie,jb:je), intent(in) :: tsk
      real, dimension(ibl:iel,jbl:jel), intent(in) :: mavail,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw
      real, dimension(ibl:iel,jbl:jel,num_soil_layers), intent(in) :: tslb
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(in) :: qa
      real, intent(in), dimension(ibp3:iep3,kbp3:kep3,np3a) :: p3a
      real, intent(in), dimension(ibp3:iep3,jbp3:jep3,kbp3:kep3,np3o) :: p3o
      real, dimension(ibc:iec,jbc:jec,kbc:kec), intent(in) :: kmh,kmv,khh,khv
      real, dimension(ibc:iec,jbc:jec,kbc:kec), intent(in) :: cme
      real, dimension(ibt:iet,jbt:jet,kbt:ket), intent(in) :: tkea
      real, dimension(ibr:ier,jbr:jer,kbr:ker), intent(in) :: swten,lwten,cldfra
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      real, dimension(ni,nj), intent(in) :: radsw,rnflx,radswnet,radlwin,dsr,olr
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt), intent(in) :: pta
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: effc,effi,effs,effr,effg,effis
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: lwupt,lwdnt,lwupb,lwdnb
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: swupt,swdnt,swupb,swdnb
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: lwcf,swcf,coszen
      integer, intent(in) :: num_soil_layers
      real, dimension(ibl:iel,jbl:jel), intent(in) :: u10,v10,s10,t2,q2,hpbl,zol,mol,rmol,br,brcr,wscale,wscaleu,phim,phih,psim,psih,psiq,wspd,qsfc,wstar,delta,prkpp,fm,fh
      real, dimension(ibl:iel,jbl:jel), intent(in) :: mznt,taux,tauy
      real, intent(in), dimension(ibmyj:iemyj,jbmyj:jemyj) :: z0base,thz0,qz0,uz0,vz0,u10e,v10e
      real, intent(in), dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj,el_myj
      real, dimension(ib:ie,jb:je), intent(in) :: znt,ust,stau,tst,qst,z0t,z0q,u1,v1,s1
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl, &
                                                                  qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,  &
                                                                  edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: vdfg,maxmf
      integer, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: nupdraft,ktop_plume
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
      integer, intent(inout), dimension(d3t) :: reqt
      real, intent(inout), dimension(ib:ie,jb:je) :: dum2d
      double precision, intent(inout), dimension(kb:ke) :: dumk1,dumk2
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
      real, intent(inout) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag
      real, intent(inout) , dimension(ib2d:ie2d,jb2d:je2d,nout2d) :: out2d
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

      integer, intent(in), dimension(ib:ie,jb:je) :: cir
      real, intent(in), dimension(ib:ie,jb:je) :: crr,cangle

      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib,kmaxib) :: hflxw,hflxe,hflxs,hflxn
      integer, intent(in), dimension(ib:ie,jb:je) :: recy_cap,recy_inj

      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      real, intent(in), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg

      real, intent(in), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: kmwk

      real, intent(in), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten

      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

      integer :: i,j,k,n,nn,im,ip,nloop,nmax
      integer :: ncid,time_index,varid
      real :: tnew,pnew,thold,thnew,rdt,qv,ql,thv,ucrit
      real :: tem,r1,r2,epsd,pint,plast,dtdz,tmax,zibar
      character(len=8) :: text1
      character(len=30) :: text2
      logical, parameter :: dosfcflx = .true.
      logical :: dothis,doit
#ifdef MPI
      character(len=80) sname,uname,vname,wname
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, parameter :: nlim = 1000
      integer :: reqs
#endif
      real, dimension(:), allocatable :: pfoo,tfoo,qfoo
      real :: zlcl, zlfc, zel , psource , tsource , qvsource

      integer :: hloop
      real :: bri,uavg,vavg,thvflux,tmpprt,thf1,thf2,gamfac,ws
      real, dimension(:,:), allocatable :: tmpqsfc,govthv,thv1,thvsfc

      real :: tpsim,tpsih,tphim,tphih,zeta,ustbar,zntbar,molbar

      integer, parameter :: unum = 52
      integer, parameter :: vnum = 53
      integer, parameter :: wnum = 54

      ! write effective radii:
      logical, parameter :: doeff = .false.

!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!#endif
      if( myid.eq.0 ) print *,'  Entering writeout ... '
!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!#endif

!--------------------------------------------------------------
!  writeout data

      levs_output = maxk

211   format(i1.1)
212   format(i2.2)
213   format(i3.3)
214   format(i4.4)
215   format(i5.5)

      opens = .false.
      openu = .false.
      openv = .false.
      openw = .false.

      ncid = 1
      time_index = 1

      if( myid.eq.0 ) print *,'  nwrite = ',nwrite

#ifdef MPI
      !  limit to "nlim" writes at a time:
      IF( output_filetype.eq.3 )THEN
      IF( numprocs.gt.nlim )THEN
        doit = .false.
        IF( myid.ge.nlim )THEN
          call MPI_IRECV(doit,1,mpi_logical,myid-nlim,999999,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ENDIF
      ENDIF
      ENDIF
#endif


  if(output_filetype.ge.2)then
    srec=1
    urec=1
    vrec=1
    wrec=1
  endif


!-----------------------------------------------------------------------

      if( output_format.eq.1 )then

        ! write metadata:
        call     write_grads_metadata(nwrite,mrec,nstep,mtime,dt,adt,ndt,domainlocx,domainlocy)

      endif

!-----------------------------------------------------------------------

      IF( imove.eq.1 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ! get ground-relative winds:
          ugr(i,j,k) = ua(i,j,k)+umove
          vgr(i,j,k) = va(i,j,k)+vmove
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ugr(i,j,k) = ua(i,j,k)
          vgr(i,j,k) = va(i,j,k)
        enddo
        enddo
        enddo
      ENDIF


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    nmax = 1
    if( output_format.eq.2 ) nmax = 2


    bigloop:  &
    DO nloop = 1 , nmax

      if( myid.eq.0 ) print *,'  nloop = ',nloop

      if( fnum.eq.71 )then
        dointerp = .true.
      else
        dointerp = .false.
      endif

#ifdef NETCDF
      IF(output_format.eq.2)THEN
      IF( nloop.eq.2 )THEN
        ! netcdf stuff:
        if( output_filetype.eq.3 .or. myid.eq.0 )then
            if(dowr) write(outfile,*) '  calling netcdf_prelim ... '
            call netcdf_prelim(rtime,nwrite,fnum,ncid,time_index,qname,                      &
                               name_output,desc_output,unit_output,grid_output,cmpr_output,  &
                               xh,xf,yh,yf,xfref,yfref,sigma,sigmaf,zs,zh,zf,                &
                               dum1(ib,jb,kb),dum2(ib,jb,kb),dum3(ib,jb,kb),dum4(ib,jb,kb),  &
                               dum5(ib,jb,kb),dat2(1,1),dat2(1,2))
            if(dowr) write(outfile,*) '  ... done '
          opens = .true.
        endif
      ENDIF
      ENDIF
#endif

      n_out = 0

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! 2d s-staggered variables:

  !.............................................

    if(output_rain   .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'rain    '
      desc_output(n_out) = 'accumulated surface rainfall'
      unit_output(n_out) = 'cm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      cmpr_output(n_out) = .true.

        call   write2d(rain(ib,jb,1),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_rain   .eq.1 .and. cuparam.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'raincv  '
      desc_output(n_out) = 'accumulated surface rainfall (convective param.)'
      unit_output(n_out) = 'cm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      cmpr_output(n_out) = .true.

        call   write2d(raincv(ib,jb,1),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_rain   .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'prate'
      desc_output(n_out) = 'surface precipitation rate'
      unit_output(n_out) = 'kg/m2/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      cmpr_output(n_out) = .true.

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = prate(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sws    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'sws     '
      desc_output(n_out) = 'max horiz wind speed at lowest model level'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sws(i,j,1)
          enddo
          enddo

        call   write2d(sws(ib,jb,1),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_svs    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'svs     '
      desc_output(n_out) = 'max vert vorticity at lowest model level'
      unit_output(n_out) = '1/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = svs(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sps    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'sps     '
      desc_output(n_out) = 'min pressure at lowest model level'
      unit_output(n_out) = 'Pa'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sps(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_srs    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'srs     '
      desc_output(n_out) = 'max qr at lowest model level'
      unit_output(n_out) = 'kg/kg'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = srs(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sgs    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'sgs     '
      desc_output(n_out) = 'max qg at lowest model level'
      unit_output(n_out) = 'kg/kg'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sgs(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sus    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'sus     '
      desc_output(n_out) = 'max w at 5 km AGL'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sus(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_shs    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'shs     '
      desc_output(n_out) = 'max integrated updraft helicity'
      unit_output(n_out) = 'm2/s2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = shs(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_rain   .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'rain2   '
        desc_output(n_out) = 'translated surface rainfall'
        unit_output(n_out) = 'cm'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = rain(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_sws    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'sws2    '
        desc_output(n_out) = 'translated max horiz wspd at lowest model level'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sws(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_svs    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'svs2    '
        desc_output(n_out) = 'translated max vert vort at lowest model level'
        unit_output(n_out) = '1/s'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = svs(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_sps    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'sps2    '
        desc_output(n_out) = 'translated min pressure at lowest model level'
        unit_output(n_out) = 'Pa'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sps(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_srs    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'srs2    '
        desc_output(n_out) = 'translated max qr at lowest model level'
        unit_output(n_out) = 'kg/kg'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = srs(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_sgs    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'sgs2    '
        desc_output(n_out) = 'translated max qg at lowest model level'
        unit_output(n_out) = 'kg/kg'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sgs(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_sus    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'sus2    '
        desc_output(n_out) = 'translated max w at 5 km AGL'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sus(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(nrain.eq.2)then
      if(output_shs    .eq.1)then
        n_out = n_out + 1
        name_output(n_out) = 'shs2    '
        desc_output(n_out) = 'translated max integrated updraft helicity'
        unit_output(n_out) = 'm2/s2'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = shs(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    endif

  !.............................................

    if(output_uh.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'uh      '
      desc_output(n_out) = 'integrated (2-5 km) AGL) updraft helicity'
      unit_output(n_out) = 'm2/s2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          ! get height AGL:
          if( terrain_flag )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk+1
            do j=1,nj
            do i=1,ni
              dum3(i,j,k) = zh(i,j,k)-zs(i,j)
              dumw(i,j,k) = zf(i,j,k)-zs(i,j)
            enddo
            enddo
            enddo
          else
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk+1
            do j=1,nj
            do i=1,ni
              dum3(i,j,k) = zh(i,j,k)
              dumw(i,j,k) = zf(i,j,k)
            enddo
            enddo
            enddo
          endif
          if(timestats.ge.1) time_write=time_write+mytime()
          call calcuh(uf,vf,dum3,dumw,ua,va,wa,dum1(ib,jb,1),dum2,dum5,dum6, &
                      zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf)

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_coldpool.eq.1)then
          if(timestats.ge.1) time_write=time_write+mytime()
          call calccpch(zh,zf,th0,qv0,dum1(ib,jb,1),dum1(ib,jb,2),tha,qa)
    endif

  !.............................................

    if(output_coldpool.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'cpc     '
      desc_output(n_out) = 'cold pool intensity C'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_coldpool.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'cph     '
      desc_output(n_out) = 'cold pool depth h'
      unit_output(n_out) = 'm AGL'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcflx .eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'thflux  '
      desc_output(n_out) = 'surface potential temperature flux'
      unit_output(n_out) = 'K m/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = thflux(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcflx .eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'qvflux  '
      desc_output(n_out) = 'surface water vapor mixing ratio flux'
      unit_output(n_out) = 'g/g m/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = qvflux(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcflx .eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'tsk     '
      desc_output(n_out) = 'soil/ocean temperature'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tsk(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcparams.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'cd      '
      desc_output(n_out) = 'sfc exchange coeff for momentum'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = cd(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcparams.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'ch      '
      desc_output(n_out) = 'sfc exchange coeff for sensible heat'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = ch(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcparams.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'cq      '
      desc_output(n_out) = 'sfc exchange coeff for moisture'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = cq(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcparams.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'tlh     '
      desc_output(n_out) = 'horiz lengthscale for turbulence scheme'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tlh(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if( idoles .and. cm1setup.eq.4 )then

      n_out = n_out + 1
      name_output(n_out) = 'cm0'
      desc_output(n_out) = 'cm0'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   write2d(cm0,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if( betaplane.eq.1 )then

      n_out = n_out + 1
      name_output(n_out) = 'f2d'
      desc_output(n_out) = 'Coriolis parameter'
      unit_output(n_out) = '1/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = f2d(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_psfc   .eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'psfc    '
      desc_output(n_out) = 'surface pressure'
      unit_output(n_out) = 'Pa'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = psfc(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    IF( do_adapt_move )THEN

      n_out = n_out + 1
      name_output(n_out) = 'psmth   '
      desc_output(n_out) = 'smooth surface pressure for getcenter'
      unit_output(n_out) = 'Pa'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = psmth(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    if(output_zs     .eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'zs      '
      desc_output(n_out) = 'terrain height'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = zs(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_dbz    .eq.1 .and. qd_dbz.ge.1 )then

      n_out = n_out + 1
      name_output(n_out) = 'cref    '
      desc_output(n_out) = 'composite reflectivity'
      unit_output(n_out) = 'dBZ'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      cmpr_output(n_out) = .true.

          if(timestats.ge.1) time_write=time_write+mytime()
          call calccref(dum1(ib,jb,1),qdiag(ibdq,jbdq,kbdq,qd_dbz))

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif

  !.............................................

    if(output_sfcparams.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'xland   '
      desc_output(n_out) = 'land/water flag (1=land,2=water)'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = xland(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sfcparams.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'lu      '
      desc_output(n_out) = 'land use index'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = float( lu_index(i,j) )
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sfcparams.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'mavail  '
      desc_output(n_out) = 'surface moisture availability '
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = mavail(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4.or.oceanmodel.eq.2))then
      n_out = n_out + 1
      name_output(n_out) = 'tmn     '
      desc_output(n_out) = 'deep-layer soil temperature'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tmn(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.ge.1.or.oceanmodel.eq.2))then
      n_out = n_out + 1
      name_output(n_out) = 'hfx     '
      desc_output(n_out) = 'surface sensible heat flux'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = hfx(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.ge.1.or.oceanmodel.eq.2))then
      n_out = n_out + 1
      name_output(n_out) = 'qfx     '
      desc_output(n_out) = 'surface moisture flux'
      unit_output(n_out) = 'kg/m^2/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = qfx(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4.or.oceanmodel.eq.2))then
      n_out = n_out + 1
      name_output(n_out) = 'gsw     '
      desc_output(n_out) = 'downward SW flux at surface'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = gsw(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4.or.oceanmodel.eq.2))then
      n_out = n_out + 1
      name_output(n_out) = 'glw     '
      desc_output(n_out) = 'downward LW flux at surface'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = glw(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4))then
      n_out = n_out + 1
      name_output(n_out) = 'tslb1   '
      desc_output(n_out) = 'soil temperature, layer 1'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tslb(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4))then
      n_out = n_out + 1
      name_output(n_out) = 'tslb2   '
      desc_output(n_out) = 'soil temperature, layer 2'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tslb(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4))then
      n_out = n_out + 1
      name_output(n_out) = 'tslb3   '
      desc_output(n_out) = 'soil temperature, layer 3'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tslb(i,j,3)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4))then
      n_out = n_out + 1
      name_output(n_out) = 'tslb4   '
      desc_output(n_out) = 'soil temperature, layer 4'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tslb(i,j,4)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4))then
      n_out = n_out + 1
      name_output(n_out) = 'tslb5   '
      desc_output(n_out) = 'soil temperature, layer 5'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tslb(i,j,5)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
      n_out = n_out + 1
      name_output(n_out) = 'tml     '
      desc_output(n_out) = 'ocean mixed layer temperature'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tml(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
      n_out = n_out + 1
      name_output(n_out) = 'hml     '
      desc_output(n_out) = 'ocean mixed layer depth'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = hml(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
      n_out = n_out + 1
      name_output(n_out) = 'huml    '
      desc_output(n_out) = 'ocean mixed layer u velocity'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = huml(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
      n_out = n_out + 1
      name_output(n_out) = 'hvml    '
      desc_output(n_out) = 'ocean mixed layer v velocity'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = hvml(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( radopt.eq.1 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      ! nasa-goddard vars:
      n_out = n_out + 1
      name_output(n_out) = 'radsw   '
      desc_output(n_out) = 'solar radiation at surface'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = radsw(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.eq.1 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'rnflx   '
      desc_output(n_out) = 'net radiation absorbed by surface'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = rnflx(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.eq.1 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'radswnet'
      desc_output(n_out) = 'net solar radiation'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = radswnet(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.eq.1 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'radlwin '
      desc_output(n_out) = 'incoming longwave radiation'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = radlwin(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.eq.1 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'olr     '
      desc_output(n_out) = 'TOA net outgoing longwave radiation'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = olr(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.eq.1 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'dsr     '
      desc_output(n_out) = 'TOA net incoming solar radiation'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dsr(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'lwupt'
      desc_output(n_out) = 'lw flux, upward, top of atmosphere (OLR)'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = lwupt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'lwdnt'
      desc_output(n_out) = 'lw flux, downward, top of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = lwdnt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'lwupb'
      desc_output(n_out) = 'lw flux, upward, bottom of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = lwupb(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'lwdnb'
      desc_output(n_out) = 'lw flux, downward, bottom of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = lwdnb(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'swupt'
      desc_output(n_out) = 'sw flux, upward, top of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = swupt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'swdnt'
      desc_output(n_out) = 'sw flux, downward, top of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = swdnt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'swupb'
      desc_output(n_out) = 'sw flux, upward, bottom of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = swupb(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'swdnb'
      desc_output(n_out) = 'sw flux, downward, bottom of atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = swdnb(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      ! cloud forcing vars:
      n_out = n_out + 1
      name_output(n_out) = 'lwcf'
      desc_output(n_out) = 'longwave cloud forcing at top-of-atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = lwcf(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN
    if( output_sfcparams.eq.1 .and. radopt.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'swcf'
      desc_output(n_out) = 'shortwave cloud forcing at top-of-atmosphere'
      unit_output(n_out) = 'W/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = swcf(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF( radopt.ge.2 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'coszen'
      desc_output(n_out) = 'cosine of solar zenith angle'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

      call   write2d(coszen(ibr,jbr),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'u10     '
    if( imove.eq.1 )then
      desc_output(n_out) = 'diagnostic 10m u wind speed (ground-rel.)'
    else
      desc_output(n_out) = 'diagnostic 10m u wind speed'
    endif
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = u10(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'v10     '
    if( imove.eq.1 )then
      desc_output(n_out) = 'diagnostic 10m v wind speed (ground-rel.)'
    else
      desc_output(n_out) = 'diagnostic 10m v wind speed'
    endif
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = v10(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 's10'
      desc_output(n_out) = 's10'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = s10(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 't2      '
      desc_output(n_out) = 'diagnostic 2m temperature'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = t2(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'q2      '
      desc_output(n_out) = 'diagnostic 2m mixing ratio'
      unit_output(n_out) = 'g/g'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = q2(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'znt     '
    if( sfcmodel.ne.4 )then
      desc_output(n_out) = 'roughness length'
    else
      desc_output(n_out) = 'thermal roughness length (m) (NOTE: unusual naming convention)'
    endif
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = znt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. dotimeavg .and. dot2p )THEN
      n_out = n_out + 1
      name_output(n_out) = 'znttavg '
      desc_output(n_out) = 'time-average roughness length'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sfctimavg(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
    if( sfcmodel.eq.2 .or. sfcmodel.eq.3 )then
      n_out = n_out + 1
      name_output(n_out) = 'z0t     '
      desc_output(n_out) = 'roughness length for temperature'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = z0t(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
    if( sfcmodel.eq.2 .or. sfcmodel.eq.3 )then
      n_out = n_out + 1
      name_output(n_out) = 'z0q     '
      desc_output(n_out) = 'roughness length for moisture'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = z0q(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'ust     '
      desc_output(n_out) = 'friction velocity'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = ust(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. dotimeavg .and. dot2p )THEN
      n_out = n_out + 1
      name_output(n_out) = 'usttavg '
      desc_output(n_out) = 'time-average friction velocity'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sfctimavg(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. tbc.eq.3 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'ustt    '
      desc_output(n_out) = 'friction velocity at model top'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = ustt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'stau'
      desc_output(n_out) = 'surface stress'
      unit_output(n_out) = 'm^2/s^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = stau(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
    if( dosfcflx .or. use_pbl )then
    if( sfcmodel.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'tst'
      desc_output(n_out) = 'theta-star (pot temp scaling parameter in similarity theory)'
      unit_output(n_out) = 'K'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tst(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    endif
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
    if( dosfcflx .or. use_pbl )then
    if( sfcmodel.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qst'
      desc_output(n_out) = 'q-star (water vapor scaling parameter in similarity theory)'
      unit_output(n_out) = 'g/g'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = qst(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    endif
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'hpbl    '
  IF( testcase.ge.1 .and. testcase.le.7 )THEN
      desc_output(n_out) = 'PBL height (using max theta gradient)'
  ELSE
      desc_output(n_out) = 'diagnosed PBL height'
  ENDIF
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = hpbl(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
    IF( ( testcase.ge.1 .and. testcase.le.7 ) .or. testcase.eq.9 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'zi'
      desc_output(n_out) = 'estimate of boundary-layer depth (max gradient method)'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

      !--------
      !  Get zi from avg theta:

      if( imoist.eq.1 )then
        do k=1,nk
          dumk1(k) = 0.0
          do j=1,nj
          do i=1,ni
            dumk1(k) = dumk1(k) + (th0(i,j,k)+tha(i,j,k))*(1.0+repsm1*qa(i,j,k,nqv))
          enddo
          enddo
        enddo
      else
        do k=1,nk
          dumk1(k) = 0.0
          do j=1,nj
          do i=1,ni
            dumk1(k) = dumk1(k) + (th0(i,j,k)+tha(i,j,k))
          enddo
          enddo
        enddo
      endif

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,dumk1(1),nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        do k=1,nk
          dumk1(k) = dumk1(k)/dble(nx*ny)
        enddo

        tmax = -1.0e30
        do k=nk,2,-1
          dtdz = (dumk1(k)-dumk1(k-1))*rdz*mf(1,1,k)
          if( dtdz .gt. tmax )then
            tmax = dtdz
            zibar = zf(1,1,k)
          endif
        enddo

      !--------
      !  Get zi:

        do j=1,nj
        do i=1,ni
          dum2d(i,j) = 0.0
          tmax = -1.0e30
          do k=nk,2,-1
            dtdz = ((th0(i,j,k)+tha(i,j,k))-(th0(i,j,k-1)+tha(i,j,k-1)))*rdz*mf(i,j,k)
            if( dtdz .ge. tmax .and. abs(zf(i,j,k)-zibar).le.(0.5*zibar) )then
              tmax = dtdz
              dum2d(i,j) = zf(i,j,k)
            endif
          enddo
        enddo
        enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN

      n_out = n_out + 1
      name_output(n_out) = 'zol     '
      desc_output(n_out) = 'z/L (z over Monin-Obukhov length)'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = zol(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'mol'
      desc_output(n_out) = 'Monin-Obukhov length (L)'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            if( abs(rmol(i,j)).le.1.0e-10 )then
              dum2d(i,j) = sign( 1.0e10 , rmol(i,j) )
            else
              dum2d(i,j) = 1.0/rmol(i,j)
            endif
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. dotimeavg .and. dot2p )THEN
      n_out = n_out + 1
      name_output(n_out) = 'moltavg '
      desc_output(n_out) = 'time-average Monin-Obukhov length (L)'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = sfctimavg(i,j,3)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'br      '
      desc_output(n_out) = 'bulk Richardson number in surface layer'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = br(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. use_pbl )THEN
      n_out = n_out + 1
      name_output(n_out) = 'brcr'
      desc_output(n_out) = 'critical bulk Richardson number for PBL'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = brcr(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. use_pbl )THEN
      n_out = n_out + 1
      name_output(n_out) = 'wstar'
      desc_output(n_out) = 'w-star from pbl code'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = wstar(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. (ipbl.eq.1.or.ipbl.eq.4.or.ipbl.eq.5) )THEN
      n_out = n_out + 1
      name_output(n_out) = 'delta'
      desc_output(n_out) = 'thickness of entrainment zone in pbl scheme'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = delta(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1 .and. (ipbl.eq.1.or.ipbl.eq.3) )THEN
      n_out = n_out + 1
      name_output(n_out) = 'prkpp'
      desc_output(n_out) = 'Prandtl number in KPP pbl'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = prkpp(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1 .and. ipbl.eq.3 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'wscale'
      desc_output(n_out) = 'wscale from GFS-EDMF'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = wscale(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1 .and. ipbl.eq.3 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'wscaleu'
      desc_output(n_out) = 'wscaleu from GFS-EDMF'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = wscaleu(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'psim'
      desc_output(n_out) = 'similarity nondimen wind shear at lowest model level'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = psim(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'psih'
      desc_output(n_out) = 'similarity nondimen temp grad at lowest model level'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = psih(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'psiq    '
      desc_output(n_out) = 'similarity stability function (moisture) at lowest model level'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = psiq(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'fm'
      desc_output(n_out) = 'fm (from surface layer: ln(z/z0)-psim)'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = fm(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'fh'
      desc_output(n_out) = 'fh (from surface layer: ln(z/z0t)-psih)'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = fh(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'qsfc    '
      desc_output(n_out) = 'land/ocean water vapor mixing ratio'
      unit_output(n_out) = 'g/g'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = qsfc(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF(output_sfcdiags.eq.1)THEN
      n_out = n_out + 1
      name_output(n_out) = 'wspd    '
      desc_output(n_out) = 'sfc layer wind speed (with gust)   '
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = wspd(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.4 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'mznt'
      desc_output(n_out) = 'momentum roughness length (m) (NOTE: unusual naming convention)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = mznt(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.4 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'taux'
      desc_output(n_out) = 'Instantaneous stress along X direction'
      unit_output(n_out) = 'Kg/m/s^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = taux(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.4 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'tauy'
      desc_output(n_out) = 'Instantaneous stress along Y direction'
      unit_output(n_out) = 'Kg/m/s^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = tauy(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'z0base'
      desc_output(n_out) = 'z0base'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = z0base(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'thz0'
      desc_output(n_out) = 'thz0'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = thz0(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'qz0'
      desc_output(n_out) = 'qz0'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = qz0(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'uz0'
      desc_output(n_out) = 'uz0'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = uz0(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'vz0'
      desc_output(n_out) = 'vz0'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = vz0(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'u10e'
      desc_output(n_out) = 'u10e'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = u10e(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................

    IF( output_sfcdiags.eq.1 .and. sfcmodel.eq.7 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'v10e'
      desc_output(n_out) = 'v10e'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = v10e(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    ENDIF

  !.............................................


    IF( imoist.eq.1 .and. output_lwp.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'cwp     '
      desc_output(n_out) = 'cloud water path'
      unit_output(n_out) = 'kg/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
            do i=1,ni
              dum1(i,j,1) = 0.0
            enddo
            if( nqc.ge.1 )then
              do k=1,nk
              do i=1,ni
                dum1(i,j,1) = dum1(i,j,1) + rho(i,j,k)*qa(i,j,k,nqc)*dz*rmh(i,j,k)
              enddo
              enddo
            endif
          enddo

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. output_lwp.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'lwp     '
      desc_output(n_out) = 'liquid water path'
      unit_output(n_out) = 'kg/m^2'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
            do i=1,ni
              dum1(i,j,1) = 0.0
            enddo
            if( nqc.ge.1 .and. nqr.ge.1 )then
              do k=1,nk
              do i=1,ni
                dum1(i,j,1) = dum1(i,j,1) + rho(i,j,k)*(qa(i,j,k,nqc)+qa(i,j,k,nqr))*dz*rmh(i,j,k)
              enddo
              enddo
            endif
          enddo

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. output_pwat.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'pwat    '
      desc_output(n_out) = 'precipitable water'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
            do i=1,ni
              dum1(i,j,1) = 0.0
            enddo
            if( nqv.ge.1 )then
              do k=1,nk
              do i=1,ni
                                                                                ! 1000 kg/m3
                dum1(i,j,1) = dum1(i,j,1) + rho(i,j,k)*qa(i,j,k,nqv)*dz*rmh(i,j,k)/1000.0
              enddo
              enddo
            endif
          enddo

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    IF( imoist.eq.1 )THEN
    IF( output_cape.eq.1 .or. output_cin.eq.1 .or. output_lcl.eq.1 .or. output_lfc.eq.1 )THEN

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,pfoo,tfoo,qfoo,zel,psource,tsource,qvsource)
          DO j=1,nj
          DO i=1,ni

            allocate( pfoo(nk+1) )
            allocate( tfoo(nk+1) )
            allocate( qfoo(nk+1) )

            do k=1,nk
              pfoo(k+1) = 0.01*prs(i,j,k)
              tfoo(k+1) = (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k)) - 273.15
              qfoo(k+1) = qa(i,j,k,nqv)
            enddo

            pfoo(1) = cgs1*pfoo(2)+cgs2*pfoo(3)+cgs3*pfoo(4)
            tfoo(1) = cgs1*tfoo(2)+cgs2*tfoo(3)+cgs3*tfoo(4)
            qfoo(1) = cgs1*qfoo(2)+cgs2*qfoo(3)+cgs3*qfoo(4)

            ! dum1(1) = cape
            ! dum1(2) = cin
            ! dum2(1) = lcl
            ! dum2(2) = lfc

            call getcape( 3 , nk+1 , pfoo , tfoo , qfoo , dum1(i,j,1) , dum1(i,j,2) ,   &
                          dum2(i,j,1), dum2(i,j,2), zel , psource , tsource , qvsource )

            deallocate( pfoo )
            deallocate( tfoo )
            deallocate( qfoo )

          ENDDO
          ENDDO

    ENDIF
    ENDIF

  !.............................................

    IF( imoist.eq.1 )THEN
    IF( output_cape.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'cape    '
      desc_output(n_out) = 'convective available potential energy'
      unit_output(n_out) = 'J/kg'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF
    ENDIF

  !.............................................

    IF( imoist.eq.1 )THEN
    IF( output_cin.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'cin     '
      desc_output(n_out) = 'convective inhibition'
      unit_output(n_out) = 'J/kg'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum1(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF
    ENDIF

  !.............................................

    IF( imoist.eq.1 )THEN
    IF( output_lcl.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'lcl     '
      desc_output(n_out) = 'lifted condensation level'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum2(i,j,1)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF
    ENDIF

  !.............................................

    IF( imoist.eq.1 )THEN
    IF( output_lfc.eq.1 )THEN

      n_out = n_out + 1
      name_output(n_out) = 'lfc     '
      desc_output(n_out) = 'level of free convection'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = dum2(i,j,2)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF
    ENDIF

  !.............................................

!    if( pmin.lt.40000.0 )then
!
!      n_out = n_out + 1
!      name_output(n_out) = 'wa500   '
!      desc_output(n_out) = 'vertical velocity at 500 mb'
!      unit_output(n_out) = 'm/s'
!      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
!
!          !$omp parallel do default(shared)  &
!          !$omp private(i,j,k,pint,plast)
!          do j=1,nj
!          do i=1,ni
!            pint = 1.0e30
!            k = 1
!            do while( pint.gt.50000.0 .and. k.lt.nk )
!              plast = pint
!              k = k + 1
!              pint = 0.5*(prs(i,j,k-1)+prs(i,j,k))
!            enddo
!            dum2d(i,j) = wa(i,j,k-1)+(wa(i,j,k)-wa(i,j,k-1))  &
!                                    *(50000.0-plast)  &
!                                    /(pint-plast)
!          enddo
!          enddo
!
!        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
!
!    endif

  !.............................................

    IF( do_ib )THEN

      n_out = n_out + 1
      name_output(n_out) = 'kbdy'
      desc_output(n_out) = 'kbdy'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = kbdy(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    IF( do_recycle )THEN

      n_out = n_out + 1
      name_output(n_out) = 'recy_cap'
      desc_output(n_out) = 'recy_cap'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = recy_cap(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    IF( do_recycle )THEN

      n_out = n_out + 1
      name_output(n_out) = 'recy_inj'
      desc_output(n_out) = 'recy_inj'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dum2d(i,j) = recy_inj(i,j)
          enddo
          enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    ENDIF

  !.............................................

    out2dcheck:  &
    IF( nout2d.ge.1 .and. ie2d.gt.1 .and. je2d.gt.1 )THEN

      ! arbitrary output (out2d array)

      do n=1,nout2d

        n_out = n_out + 1
        text1 = 'out2d   '
        if(n.lt.10)then
          write(text1(6:6),211) n
        elseif(n.lt.100)then
          write(text1(6:7),212) n
        elseif(n.lt.1000)then
          write(text1(6:8),213) n
        else
          print *,'  nout2d is too large '
          call stopcm1
        endif
        name_output(n_out) = text1
        desc_output(n_out) = '2d output'
        unit_output(n_out) = 'unknown'
        grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

                !$omp parallel do default(shared)  &
                !$omp private(i,j)
                do j=1,nj
                do i=1,ni
                  dum2d(i,j) = out2d(i,j,n)
                enddo
                enddo

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

      enddo

    ENDIF  out2dcheck

  !.............................................

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! 3d s-staggered variables:

  !.............................................

    if(output_zh     .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'zhval   '
      desc_output(n_out) = 'height on model levels'
      unit_output(n_out) = 'm'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          if( fnum.eq.71 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = sigma(k)-zs(i,j)
            enddo
            enddo
            enddo
          else
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = zh(i,j,k)
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_th     .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'th      '
      desc_output(n_out) = 'potential temperature'
      unit_output(n_out) = 'K'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = th0(i,j,k)+tha(i,j,k)
          enddo
          enddo
          enddo

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    IF( dotimeavg .and. ttav.ge.1 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'thtavg'
      desc_output(n_out) = 'time-average theta'
      unit_output(n_out) = 'K'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dum1 = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = th0(i,j,k)+timavg(i,j,k,ttav)
          enddo
          enddo
          enddo

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    ENDIF

  !.............................................

    if(output_thpert .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'thpert  '
      desc_output(n_out) = 'potential temperature perturbation'
      unit_output(n_out) = 'K'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tha(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_prs    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'prs     '
      desc_output(n_out) = 'pressure'
      unit_output(n_out) = 'Pa'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
           dum1(i,j,k) = prs(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_prspert.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'prspert '
      desc_output(n_out) = 'pressure perturbation'
      unit_output(n_out) = 'Pa'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = prs(i,j,k)-prs0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pi     .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'pi      '
      desc_output(n_out) = 'nondimensional pressure'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          if( psolver.eq.6 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = (prs(i,j,k)*rp00)**rovcp
            enddo
            enddo
            enddo
          else
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = pi0(i,j,k)+ppi(i,j,k)
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pipert .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'pipert  '
      desc_output(n_out) = 'nondimensional pressure perturbation'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          if( psolver.eq.6 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = (prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
            enddo
            enddo
            enddo
          else
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = ppi(i,j,k)
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if( psolver.eq.4 .or. psolver.eq.5 .or. psolver.eq.6 )then
      n_out = n_out + 1
      name_output(n_out) = 'phi'

      if(psolver.eq.4 )  &
      desc_output(n_out) = 'pressure variable for anelastic equations'

      if(psolver.eq.5 )  &
      desc_output(n_out) = 'pressure variable for incompressible equations'

      if(psolver.eq.6 )  &
      desc_output(n_out) = 'pressure variable for compr.-Bouss. equations'

      unit_output(n_out) = 'm2/s2'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = phi2(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_rho    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'rho     '
      desc_output(n_out) = 'dry-air density'
      unit_output(n_out) = 'kg/m^3'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = rho(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_rhopert.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'rhopert '
      desc_output(n_out) = 'dry-air density perturbation'
      unit_output(n_out) = 'kg/m^3'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = rho(i,j,k)-rho0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(iptra         .eq.1)then

      do n=1,npt

        text1='pt      '
        if(n.le.9)then
          write(text1(3:3),155) n
155       format(i1.1)
        elseif(n.le.99)then
          write(text1(3:4),154) n
154       format(i2.2)
        else
          write(text1(3:5),153) n
153       format(i3.3)
        endif
        n_out = n_out + 1
        name_output(n_out) = text1
        desc_output(n_out) = 'passive tracer mixing ratio'
        unit_output(n_out) = 'kg/kg'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        cmpr_output(n_out) = .true.

                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk
                do j=1,nj
                do i=1,ni
                  dum1(i,j,k) = pta(i,j,k,n)
                enddo
                enddo
                enddo

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)

      enddo

    endif

  !.............................................

    if(output_qv     .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'qv      '
      desc_output(n_out) = 'water vapor mixing ratio'
      unit_output(n_out) = 'kg/kg'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qa(i,j,k,nqv)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    IF( dotimeavg .and. qtav.ge.1 .and. imoist.eq.1 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'qvtavg'
      desc_output(n_out) = 'time-average qv'
      unit_output(n_out) = 'kg/kg'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dum1 = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = timavg(i,j,k,qtav)
          enddo
          enddo
          enddo

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    ENDIF

  !.............................................

    if(output_qvpert .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'qvpert  '
      desc_output(n_out) = 'water vapor mixing ratio perturbation'
      unit_output(n_out) = 'kg/kg'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qa(i,j,k,nqv)-qv0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    ! all moisture variables (except qv):

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then

      do n=1,numq

        if(n.ne.nqv)then

          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)
          n_out = n_out + 1
          name_output(n_out) = text1
          desc_output(n_out) = text2
          unit_output(n_out) = qunit(n)
          grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
          cmpr_output(n_out) = .true.


                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk
                do j=1,nj
                do i=1,ni
                  dum1(i,j,k) = qa(i,j,k,n)
                enddo
                enddo
                enddo

                ! for NSSL microphysics:
                IF( ptype.ge.26 .and. ptype.le.28 )THEN
                IF ( output_n_as_m3 > 0 ) THEN
                ! convert output to number per m^3
                IF ( trim(qunit(n)) == '#/kg' ) THEN
                  unit_output(n_out) = '#/m3'

                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk
                do j=1,nj
                do i=1,ni
                  dum1(i,j,k) = qa(i,j,k,n)*rho(i,j,k)
                enddo
                enddo
                enddo
                ENDIF
                ENDIF
                ENDIF

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)


        endif

      enddo

    endif
    ENDIF

!--------------------------------------------------------------------------------
!  begin P3 diagnostics

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 )then
        n_out = n_out + 1
        name_output(n_out) = 'p3_vmi'
        desc_output(n_out) = 'P3: mean mass weighted ice fallspeed'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,1)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 )then
        n_out = n_out + 1
        name_output(n_out) = 'p3_di'
        desc_output(n_out) = 'P3: mean mass weighted ice size'
        unit_output(n_out) = 'm'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,2)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 )then
        n_out = n_out + 1
        name_output(n_out) = 'p3_rhopo'
        desc_output(n_out) = 'P3: mean mass weighted ice density'
        unit_output(n_out) = 'kg/m3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,3)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.52 )then
        n_out = n_out + 1
        name_output(n_out) = 'p3_vmi2'
        desc_output(n_out) = 'P3: mean mass weighted ice fallspeed, 2nd ice cat.'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,4)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.52 )then
        n_out = n_out + 1
        name_output(n_out) = 'p3_di2'
        desc_output(n_out) = 'P3: mean mass weighted ice size, 2nd ice cat.'
        unit_output(n_out) = 'm'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,5)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.52 )then
        n_out = n_out + 1
        name_output(n_out) = 'p3_rhopo2'
        desc_output(n_out) = 'P3: mean mass weighted ice density, 2nd ice cat.'
        unit_output(n_out) = 'kg/m3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,6)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

!  end P3 diagnostics
!--------------------------------------------------------------------------------
!  begin ISHMAEL diagnostics


  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'vmi3d_1'
        desc_output(n_out) = 'ISHMAEL: planar-nucleated mass-weighted fall speeds'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,1)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'di3d_1'
        desc_output(n_out) = 'ISHMAEL: planar-nucleated mass-weighted maximum diameter'
        unit_output(n_out) = 'm'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,2)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'rhopo3d_1'
        desc_output(n_out) = 'ISHMAEL: planar-nucleated mass-weighted effective density'
        unit_output(n_out) = 'kg/m3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,3)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'phii3d_1'
        desc_output(n_out) = 'ISHMAEL: planar-nucleated number-weighted aspect ratio'
        unit_output(n_out) = 'unitless'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,4)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'vmi3d_2'
        desc_output(n_out) = 'ISHMAEL: columnar-nucleated mass-weighted fall speeds'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,5)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'di3d_2'
        desc_output(n_out) = 'ISHMAEL: columnar-nucleated mass-weighted maximum diameter'
        unit_output(n_out) = 'm'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,6)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'rhopo3d_2'
        desc_output(n_out) = 'ISHMAEL: columnar-nucleated mass-weighted effective density'
        unit_output(n_out) = 'kg/m3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,7)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'phii3d_2'
        desc_output(n_out) = 'ISHMAEL: columnar-nucleated number-weighted aspect ratio'
        unit_output(n_out) = 'unitless'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,8)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'vmi3d_3'
        desc_output(n_out) = 'ISHMAEL: aggregate mass-weighted fall speeds'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,9)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'di3d_3'
        desc_output(n_out) = 'ISHMAEL: aggregate mass-weighted maximum diameter'
        unit_output(n_out) = 'm'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,10)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'rhopo3d_3'
        desc_output(n_out) = 'ISHMAEL: aggregate mass-weighted effective density'
        unit_output(n_out) = 'kg/m3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,11)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

  !.............................................

    IF( imoist.eq.1 .and. numq.gt.1 )THEN
    if(output_q      .eq.1)then
      if( ptype.eq.55 )then
        n_out = n_out + 1
        name_output(n_out) = 'phii3d_3'
        desc_output(n_out) = 'ISHMAEL: aggregate number-weighted aspect ratio'
        unit_output(n_out) = 'unitless'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = p3o(i,j,k,12)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
    ENDIF

!  end ISHMAEL diagnostics
!--------------------------------------------------------------------------------

  !.............................................

    if(output_dbz    .eq.1 .and. qd_dbz.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'dbz     '
      desc_output(n_out) = 'reflectivity'
      unit_output(n_out) = 'dBZ'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      cmpr_output(n_out) = .true.

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_dbz)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_buoyancy.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'buoyancy'
      desc_output(n_out) = 'buoyancy'
      unit_output(n_out) = 'm/s^2'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,nn)
          do k=1,maxk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = g*tha(i,j,k)/th0(i,j,k)
            enddo
            enddo
            IF(imoist.eq.1)THEN
              do j=1,nj
              do i=1,ni
                dum1(i,j,k) = dum1(i,j,k)+g*repsm1*(qa(i,j,k,nqv)-qv0(i,j,k))
              enddo
              enddo
              IF(nql1.ge.1)THEN
              do nn=nql1,nql2
                do j=1,nj
                do i=1,ni
                  dum1(i,j,k) = dum1(i,j,k)-g*qa(i,j,k,nn)
                enddo
                enddo
              enddo
              ENDIF
              IF(iice.eq.1)THEN
              do nn=nqs1,nqs2
                do j=1,nj
                do i=1,ni
                  dum1(i,j,k) = dum1(i,j,k)-g*qa(i,j,k,nn)
                enddo
                enddo
              enddo
              ENDIF
            ENDIF
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)

    endif

  !.............................................

    if(output_uinterp.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'uinterp '
      desc_output(n_out) = 'u interpolated to scalar points (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = 0.5*(ua(i,j,k)+ua(i+1,j,k))
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_vinterp.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'vinterp '
      desc_output(n_out) = 'v interpolated to scalar points (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = 0.5*(va(i,j,k)+va(i,j+1,k))
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_winterp.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'winterp '
      desc_output(n_out) = 'w interpolated to scalar points'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = 0.5*(wa(i,j,k)+wa(i,j,k+1))
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    ! radial and tangential velocity:
    doit = .false.

    docyl:  &
    IF( doit )THEN
    if( doazimavg .or. do_adapt_move )then

              ! dum1 = u at s pts
              ! dum2 = v at s pts
              ! dum3 = ur at s pts
              ! dum4 = vt at s pts

        do k=1,maxk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.5*(ua(i,j,k)+ua(i+1,j,k)) + umove
          dum2(i,j,k) = 0.5*(va(i,j,k)+va(i,j+1,k)) + vmove
          if( abs(crr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(cangle(i,j))+dum1(i,j,k)*cos(cangle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(cangle(i,j))-dum1(i,j,k)*sin(cangle(i,j))
          endif
        enddo
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

      n_out = n_out + 1
      name_output(n_out) = 'urad'
      desc_output(n_out) = 'radial velocity (grnd-rel)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(dum3,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)

  !.............................................

      n_out = n_out + 1
      name_output(n_out) = 'vtan'
      desc_output(n_out) = 'tangential velocity (grnd-rel)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(dum4,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)

  !.............................................

      ucrit = -2.0

      ! inflow layer depth:
      do j=1,nj
      do i=1,ni
        if( dum3(i,j,1).lt.ucrit )then
          k = 1
          do while( dum3(i,j,k).lt.ucrit .and. k.le.nk )
            k = k+1
          enddo
          dum2d(i,j) = zh(1,1,k-1)+(zh(1,1,k)-zh(1,1,k-1))  &
                                  *(    ucrit-dum3(i,j,k-1))  &
                                  /(dum3(i,j,k)-dum3(i,j,k-1))
        else
          dum2d(i,j) = 0.0
        endif
      enddo
      enddo

      n_out = n_out + 1
      name_output(n_out) = 'zinfl'
      desc_output(n_out) = 'depth of inflow layer (u < -2 m/s)'
      unit_output(n_out) = 'm'
      grid_output(n_out) = '2'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   write2d(dum2d,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif
    ENDIF  docyl

  !.............................................

    IF( dotimeavg .and. utav.ge.1 .and. vtav.ge.1 )THEN
      n_out = n_out + 1
      name_output(n_out) = 'wsptavg'
      desc_output(n_out) = 'time-average horiz wind speed'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dum1 = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = sqrt( (0.5*(timavg(i,j,k,utav)+timavg(i+1,j,k,utav)))**2  &
                               +(0.5*(timavg(i,j,k,vtav)+timavg(i,j+1,k,vtav)))**2 )
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    ENDIF

  !.............................................

    if(output_vort.eq.1)then
          if(timestats.ge.1) time_write=time_write+mytime()
          call     calcvort(xh,xf,uf,vf,zh,mh,zf,mf,                                         &
                            zs,gz,gzu,gzv,rgz,rgzu,rgzv,gxu,gyv,rds,sigma,rdsf,sigmaf,       &
                            ugr,vgr,wa,dum2 ,dum3 ,dum4 ,dum1,dum5,dum6,dum8,dum7,th0,tha,rr,  &
                            ust,znt,u1,v1,s1)
    endif

  !.............................................

    if(output_vort.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'xvort   '
      desc_output(n_out) = 'horizontal vorticity (x)'
      unit_output(n_out) = '1/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = dum2(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_vort.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'yvort   '
      desc_output(n_out) = 'horizontal vorticity (y)'
      unit_output(n_out) = '1/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = dum3(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_vort.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'zvort   '
      desc_output(n_out) = 'vertical vorticity'
      unit_output(n_out) = '1/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = dum4(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pv.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'pv      '
      desc_output(n_out) = 'potential vorticity'
      unit_output(n_out) = 'K m2/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          if(timestats.ge.1) time_write=time_write+mytime()
          call     calcvort(xh,xf,uf,vf,zh,mh,zf,mf,                                         &
                            zs,gz,gzu,gzv,rgz,rgzu,rgzv,gxu,gyv,rds,sigma,rdsf,sigmaf,       &
                            ugr,vgr,wa,dum2 ,dum3 ,dum4 ,dum1,dum5,dum6,dum8,dum7,th0,tha,rr,  &
                            ust,znt,u1,v1,s1)

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = dum8(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif


  !.............................................

    if(output_basestate.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'pi0     '
      desc_output(n_out) = 'base-state nondimensional pressure'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pi0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_basestate.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'th0     '
      desc_output(n_out) = 'base-state potential temperature'
      unit_output(n_out) = 'K'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = th0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_basestate.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'prs0    '
      desc_output(n_out) = 'base-state pressure'
      unit_output(n_out) = 'Pa'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = prs0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_basestate.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'qv0     '
      desc_output(n_out) = 'base-state water vapor mixing ratio'
      unit_output(n_out) = 'kg/kg'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qv0(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pblten.eq.1 .and. use_pbl)then
      n_out = n_out + 1
      name_output(n_out) = 'qcpten  '
      desc_output(n_out) = 'pbl tendency: cloudwater mixing ratio'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qcpten(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pblten.eq.1 .and. use_pbl)then
      n_out = n_out + 1
      name_output(n_out) = 'qipten  '
      desc_output(n_out) = 'pbl tendency: cloud ice mixing ratio'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qipten(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pblten.eq.1 .and. use_pbl .and. (ipbl.eq.4.or.ipbl.eq.5) )then
      n_out = n_out + 1
      name_output(n_out) = 'qnipten'
      desc_output(n_out) = 'pbl tendency: cloud ice number concentration'
      unit_output(n_out) = '#/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qnipten(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_pblten.eq.1 .and. use_pbl .and. (ipbl.eq.4.or.ipbl.eq.5) )then
      n_out = n_out + 1
      name_output(n_out) = 'qncpten  '
      desc_output(n_out) = 'pbl tendency: cloudwater number concentration'
      unit_output(n_out) = '#/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qncpten(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(cuparam.ge.1)then
      n_out = n_out + 1
      name_output(n_out) = 'thcuten'
      desc_output(n_out) = 'thcuten'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(thcuten(ibcu,jbcu,kbcu),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(cuparam.ge.1)then
      n_out = n_out + 1
      name_output(n_out) = 'qvcuten'
      desc_output(n_out) = 'qvcuten'
      unit_output(n_out) = 'g/g/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(qvcuten(ibcu,jbcu,kbcu),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(cuparam.ge.1)then
      n_out = n_out + 1
      name_output(n_out) = 'qccuten'
      desc_output(n_out) = 'qccuten'
      unit_output(n_out) = 'g/g/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(qccuten(ibcu,jbcu,kbcu),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(cuparam.ge.1)then
      n_out = n_out + 1
      name_output(n_out) = 'qicuten'
      desc_output(n_out) = 'qicuten'
      unit_output(n_out) = 'g/g/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(qicuten(ibcu,jbcu,kbcu),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(cuparam.ge.1)then
      n_out = n_out + 1
      name_output(n_out) = 'ucuten'
      desc_output(n_out) = 'ucuten'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(ucuten(ibcu,jbcu,kbcu),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(cuparam.ge.1)then
      n_out = n_out + 1
      name_output(n_out) = 'vcuten'
      desc_output(n_out) = 'vcuten'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writes(vcuten(ibcu,jbcu,kbcu),fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_radten.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'swten   '
      desc_output(n_out) = 'temperature tendency, sw radiation'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = swten(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_radten.eq.1)then

      n_out = n_out + 1
      name_output(n_out) = 'lwten   '
      desc_output(n_out) = 'temperature tendency, lw radiation'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = lwten(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif

  !.............................................

    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
      n_out = n_out + 1
      name_output(n_out) = 'cldfra  '
      desc_output(n_out) = 'cloud fraction from radiation scheme'
      unit_output(n_out) = 'nondimensional'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = cldfra(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif

  !.............................................

    if( doeff )then
    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
        n_out = n_out + 1
        name_output(n_out) = 'effc'
        desc_output(n_out) = 'effc'
        unit_output(n_out) = 'micron'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = effc(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              ! convert to microns:
              dum1(i,j,k) = dum1(i,j,k)*1.0e6
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif
    endif

  !.............................................

    if( doeff )then
    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
        n_out = n_out + 1
        name_output(n_out) = 'effi'
        desc_output(n_out) = 'effi'
        unit_output(n_out) = 'micron'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = effi(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              ! convert to microns:
              dum1(i,j,k) = dum1(i,j,k)*1.0e6
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif
    endif

  !.............................................

    if( doeff )then
    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
        n_out = n_out + 1
        name_output(n_out) = 'effs'
        desc_output(n_out) = 'effs'
        unit_output(n_out) = 'micron'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = effs(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              ! convert to microns:
              dum1(i,j,k) = dum1(i,j,k)*1.0e6
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif
    endif

  !.............................................

    if( doeff )then
    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
    if( ptype.eq.5 )then
        n_out = n_out + 1
        name_output(n_out) = 'effr'
        desc_output(n_out) = 'effr'
        unit_output(n_out) = 'micron'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = effr(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              ! convert to microns:
              dum1(i,j,k) = dum1(i,j,k)*1.0e6
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif
    endif
    endif

  !.............................................

    if( doeff )then
    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
    if( ptype.eq.5 )then
        n_out = n_out + 1
        name_output(n_out) = 'effg'
        desc_output(n_out) = 'effg'
        unit_output(n_out) = 'micron'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = effg(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              ! convert to microns:
              dum1(i,j,k) = dum1(i,j,k)*1.0e6
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif
    endif
    endif

  !.............................................

    if( doeff )then
    if(output_radten.eq.1)then
    if( radopt.eq.1 .or. radopt.eq.2 )then
    if( ptype.eq.5 )then
        n_out = n_out + 1
        name_output(n_out) = 'effis'
        desc_output(n_out) = 'effis'
        unit_output(n_out) = 'micron'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = effis(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,maxk
            do j=1,nj
            do i=1,ni
              ! convert to microns:
              dum1(i,j,k) = dum1(i,j,k)*1.0e6
            enddo
            enddo
            enddo
          endif
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
    endif
    endif
    endif

  !.............................................

!!!    if( radopt.eq.1 .or. radopt.eq.2 )then
!!!        n_out = n_out + 1
!!!        name_output(n_out) = 'o3'
!!!        desc_output(n_out) = 'ozone'
!!!        unit_output(n_out) = 'unk'
!!!        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
!!!
!!!        call   writes(o30,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
!!!    endif

  !.............................................

    doit = .false.
    if( doit )then
    dobrz:  &
    IF( ipbl.ge.1 )THEN

        n_out = n_out + 1
        name_output(n_out) = 'brz'
        desc_output(n_out) = 'sfc-layer bulk Richardson number'
        unit_output(n_out) = 'nondimensional'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

      allocate( tmpqsfc(ib:ie,jb:je) )
      tmpqsfc = 0.0
      allocate(  govthv(ib:ie,jb:je) )
      govthv = 0.0
      allocate(    thv1(ib:ie,jb:je) )
      thv1 = 0.0
      allocate(  thvsfc(ib:ie,jb:je) )
      thvsfc = 0.0

      if( imoist.eq.1 )then
        do j=1,nj
        do i=1,ni
          tmpqsfc(i,j) = rslf(psfc(i,j),tsk(i,j))
        enddo
        enddo
        !...
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = qa(i,j,k,nqv)
        enddo
        enddo
        enddo
      else
        do j=1,nj
        do i=1,ni
          tmpqsfc(i,j) = 0.0
        enddo
        enddo
        dum2 = 0.0
      endif

      do j=1,nj
      do i=1,ni
        thv = (th0(i,j,1)+tha(i,j,1))*(1.0+repsm1*dum2(i,j,1))
        govthv(i,j) = g/thv
      enddo
      enddo

    hlp:  &
    DO hloop = 1 , 2

      k = 1

    if( hloop.eq.1 )then
      do j=1,nj
      do i=1,ni
        thv = (th0(i,j,1)+tha(i,j,1))*(1.0+repsm1*dum2(i,j,1))
        uavg = 0.5*(ua(i,j,k)+ua(i+1,j,k))
        vavg = 0.5*(va(i,j,k)+va(i,j+1,k))
        thvsfc(i,j) = tsk(i,j)*((p00/psfc(i,j))**rovcp)*(1.0+repsm1*tmpqsfc(i,j))
        thv1(i,j) = thv
      enddo
      enddo
    else
      do j=1,nj
      do i=1,ni
        thv = (th0(i,j,1)+tha(i,j,1))*(1.0+repsm1*dum2(i,j,1))
        uavg = 0.5*(ua(i,j,k)+ua(i+1,j,k))
        vavg = 0.5*(va(i,j,k)+va(i,j+1,k))
        thf1 = (1.0+repsm1*dum2(i,j,1))
        thf2 = repsm1*(th0(i,j,1)+tha(i,j,1))
        thvflux = max( 0.0 , thf1*thflux(i,j) + thf2*qvflux(i,j) )
        ws = ( ust(i,j)**3 + 8.0*karman*0.5*govthv(i,j)*thvflux*hpbl(i,j) )**0.3333333
        ws = min( ws , ust(i,j)*16.0 )
        ws = max( ws , ust(i,j)/5.0 )
        ws = max( ws , 0.0001 )
        gamfac = 6.8/ws
        tmpprt = max( 0.0 , min( gamfac*thflux(i,j) , 3.0 )*thf1    &
                           +min( gamfac*qvflux(i,j) , 2.0e-3 )*thf2 )
        thv1(i,j) = thv + tmpprt
      enddo
      enddo
    endif

      do k=1,nk
      do j=1,nj
      do i=1,ni
        thv = (th0(i,j,k)+tha(i,j,k))*(1.0+repsm1*dum2(i,j,k))
        uavg = 0.5*(ua(i,j,k)+ua(i+1,j,k))
        vavg = 0.5*(va(i,j,k)+va(i,j+1,k))
        dum1(i,j,k) = govthv(i,j)*zh(i,j,k)*(thv-thv1(i,j))/max( 1.0 , uavg**2 + vavg**2 )
      enddo
      enddo
      enddo

    ENDDO  hlp

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)

      deallocate( tmpqsfc )
      deallocate(  govthv )
      deallocate(    thv1 )
      deallocate(  thvsfc )

    ENDIF  dobrz
    endif

  !.............................................

    ! arbitrary output (out3d array)

    out3dcheck:  &
    IF( nout3d.ge.1 .and. ie3d.gt.1 .and. je3d.gt.1 .and. ke3d.gt.1 )THEN

      do n=1,nout3d

        n_out = n_out + 1
        text1 = 'out     '
        if(n.lt.10)then
          write(text1(4:4),211) n
        elseif(n.lt.100)then
          write(text1(4:5),212) n
        elseif(n.lt.1000)then
          write(text1(4:6),213) n
        elseif(n.lt.10000)then
          write(text1(4:7),214) n
        elseif(n.lt.100000)then
          write(text1(4:8),215) n
        else
          print *,'  nout3d is too large '
          call stopcm1
        endif
        name_output(n_out) = text1
        desc_output(n_out) = '3d output'
        unit_output(n_out) = 'unknown'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk
                do j=1,nj
                do i=1,ni
                  dum1(i,j,k) = out3d(i,j,k,n)
                enddo
                enddo
                enddo

        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)

      enddo

    ENDIF  out3dcheck

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_hadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_hadv'
      if( hadvordrs.eq.3 .or. hadvordrs.eq.5 .or. hadvordrs.eq.7 .or. hadvordrs.eq.9 .or. advwenos.ge.1 )then
        desc_output(n_out) = 'pt budget: horiz advection (non-diff component)'
      else
        desc_output(n_out) = 'pot temp budget: horiz advection'
      endif
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_hadv)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_vadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_vadv'
      if( vadvordrs.eq.3 .or. vadvordrs.eq.5 .or. vadvordrs.eq.7 .or. vadvordrs.eq.9 .or. advwenos.ge.1 )then
        desc_output(n_out) = 'pt budget: vert advection (non-diff component)'
      else
        desc_output(n_out) = 'pot temp budget: vert advection'
      endif
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_vadv)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_hidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_hidiff'
      desc_output(n_out) = 'pot temp budget: horiz implicit diffusion'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_hidiff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_vidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_vidiff'
      desc_output(n_out) = 'pot temp budget: vert implicit diffusion'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_vidiff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_hediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_hediff'
      desc_output(n_out) = 'pot temp budget: horiz explicit diffusion'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_hediff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_vediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_vediff'
      desc_output(n_out) = 'pot temp budget: vert explicit diffusion'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_vediff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_hturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_hturb'
      desc_output(n_out) = 'pot temp budget: horiz parameterized turbulence'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_hturb)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_vturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_vturb'
      desc_output(n_out) = 'pot temp budget: vert parameterized turbulence'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_vturb)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_mp.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_mp'
      desc_output(n_out) = 'pot temp budget: microphysics scheme'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_mp)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_rdamp.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_rdamp'
      desc_output(n_out) = 'pot temp budget: Rayleigh damper'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_rdamp)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_nudge.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_nudge'
      desc_output(n_out) = 'pot temp budget: nudging'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_nudge)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_rad.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_rad'
      desc_output(n_out) = 'pot temp budget: radiation scheme'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_rad)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_div.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_div'
      desc_output(n_out) = 'pot temp budget: moist divergence term'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_div)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_diss.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_diss'
    if( ipbl.eq.3 .or. ipbl.eq.4 .or. ipbl.eq.5 )then
      desc_output(n_out) = 'pot temp budget: diss. heating (from PBL)'
    else
      desc_output(n_out) = 'pot temp budget: dissipative heating'
    endif
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_diss)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_pbl.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_pbl'
    if( ipbl.eq.3 .or. ipbl.eq.4 .or. ipbl.eq.5 )then
      desc_output(n_out) = 'pot tem. budget: PBL scheme (excluding diss. heating)'
    else
      desc_output(n_out) = 'pot tem. budget: PBL scheme'
    endif
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_pbl)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_lsw.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ptb_lsw'
      desc_output(n_out) = 'pot temp budget: advection by large-scale w'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_lsw)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_cond.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_cond'
        desc_output(n_out) = 'theta tendency: condensation'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_cond)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_evac.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_evac'
        desc_output(n_out) = 'theta tendency: cloudwater evaporation'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_evac)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_evar.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_evar'
        desc_output(n_out) = 'theta tendency: rainwater evaporation'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_evar)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_dep.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_dep'
        desc_output(n_out) = 'theta tendency: deposition'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_dep)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_subl.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_subl'
        desc_output(n_out) = 'theta tendency: sublimation'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_subl)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_melt.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_melt'
        desc_output(n_out) = 'theta tendency: melting'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_melt)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( td_frz.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'tt_frz'
        desc_output(n_out) = 'theta tendency: freezing'
        unit_output(n_out) = 'K/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_frz)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_thbudget.eq.1 )THEN
    if( td_efall.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'td_efall'
      desc_output(n_out) = 'temp. tendency: energy fallout terms'
      unit_output(n_out) = 'K/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tdiag(i,j,k,td_efall)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

    IF( output_fallvel.eq.1 )THEN
      if( qd_vtc.gt.0 )then
        n_out = n_out + 1
        name_output(n_out) = 'vtc     '
        desc_output(n_out) = 'terminal fall velocity: qc'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vtc)
            if( qa(i,j,k,nqc).le.qsmall ) dum1(i,j,k) = 0.0
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    ENDIF

  !.............................................

    IF( output_fallvel.eq.1 )THEN
      if( qd_vtr.gt.0 )then
        n_out = n_out + 1
        name_output(n_out) = 'vtr     '
        desc_output(n_out) = 'terminal fall velocity: qr'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vtr)
            if( qa(i,j,k,nqr).le.qsmall ) dum1(i,j,k) = 0.0
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    ENDIF

  !.............................................

    IF( output_fallvel.eq.1 )THEN
      if( qd_vts.gt.0 )then
        n_out = n_out + 1
        name_output(n_out) = 'vts     '
        desc_output(n_out) = 'terminal fall velocity: qs'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vts)
            if( qa(i,j,k,nqs).le.qsmall ) dum1(i,j,k) = 0.0
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    ENDIF

  !.............................................

    IF( output_fallvel.eq.1 )THEN
      if( qd_vtg.gt.0 )then
        n_out = n_out + 1
        name_output(n_out) = 'vtg     '
        desc_output(n_out) = 'terminal fall velocity: qg'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vtg)
            if( qa(i,j,k,nqg).le.qsmall ) dum1(i,j,k) = 0.0
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    ENDIF

  !.............................................

    IF( output_fallvel.eq.1 )THEN
      if( qd_vti.gt.0 )then
        n_out = n_out + 1
        name_output(n_out) = 'vti     '
        desc_output(n_out) = 'terminal fall velocity: qi'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vti)
            if( qa(i,j,k,nqi).le.qsmall ) dum1(i,j,k) = 0.0
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_hadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_hadv'
      if( hadvordrs.eq.3 .or. hadvordrs.eq.5 .or. hadvordrs.eq.7 .or. hadvordrs.eq.9 .or. advwenos.ge.1 )then
        desc_output(n_out) = 'qv budget: horizontal advection (non-diff component)'
      else
        desc_output(n_out) = 'qv budget: horizontal advection'
      endif
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_hadv)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_vadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_vadv'
      if( vadvordrs.eq.3 .or. vadvordrs.eq.5 .or. vadvordrs.eq.7 .or. vadvordrs.eq.9 .or. advwenos.ge.1 )then
        desc_output(n_out) = 'qv budget: vertical advection (non-diff component)'
      else
        desc_output(n_out) = 'qv budget: vertical advection'
      endif
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vadv)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF
  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_hidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_hidiff'
      desc_output(n_out) = 'qv budget: horiz implicit diffusion'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_hidiff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_vidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_vidiff'
      desc_output(n_out) = 'qv budget: vert implicit diffusion'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vidiff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_hediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_hediff'
      desc_output(n_out) = 'qv budget: horiz explicit diffusion'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_hediff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_vediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_vediff'
      desc_output(n_out) = 'qv budget: vert explicit diffusion'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vediff)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_hturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_hturb'
      desc_output(n_out) = 'qv budget: horizontal parameterized turbulence'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_hturb)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_vturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_vturb'
      desc_output(n_out) = 'qv budget: vertical parameterized turbulence'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_vturb)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF


  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_mp.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_mp'
      desc_output(n_out) = 'qv budget: microphysics scheme'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_mp)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_nudge.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_nudge'
      desc_output(n_out) = 'qv budget: nudging'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_nudge)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_pbl.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_pbl'
      desc_output(n_out) = 'qv budget: PBL scheme'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_pbl)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( qd_lsw.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'qvb_lsw'
      desc_output(n_out) = 'qv budget: advection by large-scale w'
      unit_output(n_out) = 'kg/kg/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_lsw)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( qd_cond.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'qt_cond'
        desc_output(n_out) = 'qv tendency: condensation'
        unit_output(n_out) = 'kg/kg/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_cond)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( qd_evac.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'qt_evac'
        desc_output(n_out) = 'qv tendency: cloudwater evaporation'
        unit_output(n_out) = 'kg/kg/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_evac)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( qd_evar.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'qt_evar'
        desc_output(n_out) = 'qv tendency: rainwater evaporation'
        unit_output(n_out) = 'kg/kg/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_evar)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( qd_dep.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'qt_dep'
        desc_output(n_out) = 'qv tendency: deposition'
        unit_output(n_out) = 'kg/kg/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_dep)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

  IF( output_qvbudget.eq.1 )THEN
    if( ptype.eq.5 )then
      if( qd_subl.ge.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'qt_subl'
        desc_output(n_out) = 'qv tendency: sublimation'
        unit_output(n_out) = 'kg/kg/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiag(i,j,k,qd_subl)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      endif
    endif
  ENDIF

  !.............................................

      IF( pdcomp )THEN

        n_out = n_out + 1
        name_output(n_out) = 'pipb'
        desc_output(n_out) = 'diagnosed pi-prime: buoyancy component'
        unit_output(n_out) = 'nondimensional'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pdiag(i,j,k,1)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF

  !.............................................

      IF( pdcomp )THEN
        n_out = n_out + 1
        name_output(n_out) = 'pipdl'
        desc_output(n_out) = 'diagnosed pi-prime: linear dynamic component'
        unit_output(n_out) = 'nondimensional'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pdiag(i,j,k,2)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF

  !.............................................

      IF( pdcomp )THEN
        n_out = n_out + 1
        name_output(n_out) = 'pipdn'
        desc_output(n_out) = 'diagnosed pi-prime: nonlinear dynamic component'
        unit_output(n_out) = 'nondimensional'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pdiag(i,j,k,3)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF

  !.............................................


      IF( pdcomp )THEN
        if( icor.eq.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'pipc'
        desc_output(n_out) = 'diagnosed pi-prime: Coriolis component'
        unit_output(n_out) = 'nondimensional'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pdiag(i,j,k,4)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
        endif
      ENDIF

  !.............................................

      IF( axisymm.eq.1 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'vgrad'
        desc_output(n_out) = 'gradient wind speed'
        unit_output(n_out) = 'm/s'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          dum1 = 0.0

          do k=1,maxk
          do j=1,nj
          do i=1,ni
            qv = 0.0
            ql = 0.0
            if( imoist.eq.1 )then
              qv = qa(i,j,k,nqv)
              IF(nql1.ge.1)THEN
              do nn=nql1,nql2
                ql = ql+qa(i,j,k,nn)
              enddo
              ENDIF
              if( iice.eq.1 )then
                do nn=nqs1,nqs2
                  ql = ql+qa(i,j,k,nn)
                enddo
              endif
            endif
            thv = (th0(i,j,k)+tha(i,j,k))*(1.0+reps*qv)/(1.0+qv+ql)
            ip = min( i+1 , ni )
            im = max( i-1 , 1 )
            dum1(i,j,k) = -0.5*fcor*xh(i) + sqrt( max(0.0,               &
                                0.25*fcor*fcor*xh(i)*xh(i)               &
               +xh(i)*cp*thv*(ppi(ip,j,k)-ppi(im,j,k))/(xh(ip)-xh(im))   &
                                             ) )
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF

  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qke'
        desc_output(n_out) = 'twice TKE from MYNN'
        unit_output(n_out) = 'm^2/s^2'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qke(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qwt'
        desc_output(n_out) = 'TKE vertical transport (MYNN)'
        unit_output(n_out) = 'm^2/s^3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qwt(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qshear'
        desc_output(n_out) = 'TKE Production - shear (MYNN)'
        unit_output(n_out) = 'm^2/s^3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qshear(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qbuoy'
        desc_output(n_out) = 'TKE Production - buoyancy (MYNN)'
        unit_output(n_out) = 'm^2/s^3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qbuoy(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qdiss'
        desc_output(n_out) = 'TKE dissipation (MYNN)'
        unit_output(n_out) = 'm^2/s^3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qdiss(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'dqke'
        desc_output(n_out) = 'TKE change (MYNN)'
        unit_output(n_out) = 'm^2/s^3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = dqke(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qke_adv'
        desc_output(n_out) = 'advection of twice TKE from MYNN'
        unit_output(n_out) = 'm^2/s^3'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qke_adv(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'sh3d'
        desc_output(n_out) = 'Stability function for heat (MYNN)'
        unit_output(n_out) = 'unk'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = sh3d(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'sm3d'
        desc_output(n_out) = 'Stability function for momentum (MYNN)'
        unit_output(n_out) = 'unk'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = sm3d(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'tsq'
        desc_output(n_out) = 'liquid water pottemp variance (MYNN)'
        unit_output(n_out) = 'K2'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = tsq(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qsq'
        desc_output(n_out) = 'liquid water variance (MYNN)'
        unit_output(n_out) = '(kg/kg)**2'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qsq(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF



  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'cov'
        desc_output(n_out) = 'liquid water-liquid water pottemp covariance (MYNN)'
        unit_output(n_out) = 'K kg/kg'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = cov(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF


  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 .and. icloud_bl.ge.1 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qc_bl'
        desc_output(n_out) = 'qc_bl (MYNN)'
        unit_output(n_out) = 'kg/kg'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qc_bl(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF


  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 .and. icloud_bl.ge.1 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'qi_bl'
        desc_output(n_out) = 'qi_bl (MYNN)'
        unit_output(n_out) = 'kg/kg'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = qi_bl(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF


  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 .and. icloud_bl.ge.1 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'cldfra_bl'
        desc_output(n_out) = 'cldfra_bl (MYNN)'
        unit_output(n_out) = 'unitless'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = cldfra_bl(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF


  !.............................................


      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
        n_out = n_out + 1
        name_output(n_out) = 'edmf_a'
        desc_output(n_out) = 'edmf_a (MYNN)'
        unit_output(n_out) = 'unitless'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = edmf_a(i,j,k)
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF


  !.............................................


      IF( do_ib )THEN
        n_out = n_out + 1
        name_output(n_out) = 'bndy'
        desc_output(n_out) = 'bndy'
        unit_output(n_out) = 'unitless'
        grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          dum1 = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni
            if( bndy(i,j,k) ) dum1(i,j,k) = 1.0
          enddo
          enddo
          enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
      ENDIF


  !.............................................

    doit = .false.
    if( doit )then
    IF( dotimeavg .and. dot2p )THEN
      n_out = n_out + 1
      name_output(n_out) = 'wspa'
      desc_output(n_out) = 'analytic wind speed'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 's'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,maxk
        do j=1,nj
        do i=1,ni
          ustbar = sfctimavg(i,j,1)
          zntbar = sfctimavg(i,j,2)
          molbar = sfctimavg(i,j,3)
          if( abs(molbar).gt.1.0e-6 )then
            zeta = zh(1,1,k)/molbar
            call stabil_funcs(zeta,tphim,tphih,tpsim,tpsih)
          else
            tpsim = 0.0
          endif
          dum1(i,j,k) = (ustbar/karman)*( alog(zh(1,1,k)/zntbar) - tpsim )
        enddo
        enddo
        enddo
        call   writes(dum1,fnum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out),sigma,zs,zh,dum9)
    ENDIF
    endif

  !.............................................

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! u-staggered variables:

  !.............................................

    if(output_u    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'u       '
      desc_output(n_out) = 'E-W (x) velocity (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writeu( ua ,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( dotimeavg .and. utav.ge.1 )THEN
    if(output_u    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'utavg'
      desc_output(n_out) = 'time-average u'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dumu = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = timavg(i,j,k,utav)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    if(output_upert.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'upert   '
      desc_output(n_out) = 'u perturbation (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = ua(i,j,k)-u0(i,j,k)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_basestate.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'u0      '
      desc_output(n_out) = 'base-state u (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writeu( u0 ,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_hadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_hadv'
      if( hadvordrv.eq.3 .or. hadvordrv.eq.5 .or. hadvordrv.eq.7 .or. hadvordrv.eq.9 .or. advwenov.ge.1 )then
        desc_output(n_out) = 'u budget: horizontal advection (non-diff component)'
      else
        desc_output(n_out) = 'u budget: horizontal advection'
      endif
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_hadv)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_vadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_vadv'
      if( vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1 )then
        desc_output(n_out) = 'u budget: vertical advection (non-diff component)'
      else
        desc_output(n_out) = 'u budget: vertical advection'
      endif
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_vadv)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF
  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_hidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_hidiff'
      desc_output(n_out) = 'u budget: horiz implicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_hidiff)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_vidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_vidiff'
      desc_output(n_out) = 'u budget: vert implicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_vidiff)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_hediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_hediff'
      desc_output(n_out) = 'u budget: horiz explicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_hediff)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_vediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_vediff'
      desc_output(n_out) = 'u budget: vert explicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_vediff)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_hturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_hturb'
      desc_output(n_out) = 'u budget: horizontal parameterized turbulence'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_hturb)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_vturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_vturb'
      desc_output(n_out) = 'u budget: vertical parameterized turbulence'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_vturb)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF


  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_pgrad.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_pgrad'
      desc_output(n_out) = 'u budget: pressure gradient'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_pgrad)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_rdamp.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_rdamp'
      desc_output(n_out) = 'u budget: Rayleigh damper'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_rdamp)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_cor.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_cor'
      desc_output(n_out) = 'u budget: Coriolis acceleration'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_cor)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_cent.ge.1 )then
      if( axisymm.eq.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_cent'
      desc_output(n_out) = 'u budget: centrifugal acceleration'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_cent)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_pbl.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_pbl'
      desc_output(n_out) = 'u budget: PBL scheme'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_pbl)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_ubudget.eq.1 )THEN
      if( ud_lsw.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'ub_lsw'
      desc_output(n_out) = 'u budget: advection by large-scale w'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj
          do i=1,ni+1
            dumu(i,j,k) = udiag(i,j,k,ud_lsw)
          enddo
          enddo
          enddo
        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    ! arbitrary output (out3d array)

! KLUDGE !
  if(ptype.eq.91919)then
    out3dchecku:  &
    IF( nout3d.ge.1 .and. ie3d.gt.1 .and. je3d.gt.1 .and. ke3d.gt.1 )THEN

      do n=1,nout3d

        n_out = n_out + 1
        text1 = 'out     '
        if(n.lt.10)then
          write(text1(4:4),211) n
        elseif(n.lt.100)then
          write(text1(4:5),212) n
        elseif(n.lt.1000)then
          write(text1(4:6),213) n
        elseif(n.lt.10000)then
          write(text1(4:7),214) n
        elseif(n.lt.100000)then
          write(text1(4:8),215) n
        else
          print *,'  nout3d is too large '
          call stopcm1
        endif
        name_output(n_out) = text1
        desc_output(n_out) = '3d output'
        unit_output(n_out) = 'unknown'
        grid_output(n_out) = 'u'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk
                do j=1,nj
                do i=1,ni+1
                  dumu(i,j,k) = out3d(i,j,k,n)
                enddo
                enddo
                enddo

        call   writeu(dumu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

      enddo

    ENDIF  out3dchecku
  endif
! KLUDGE !

  !.............................................

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! v-staggered variables:

  !.............................................

    if(output_v    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'v       '
      desc_output(n_out) = 'N-S (y) velocity (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writev( va ,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( dotimeavg .and. vtav.ge.1 )THEN
    if(output_v    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'vtavg'
      desc_output(n_out) = 'time-average v'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dumv = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = timavg(i,j,k,vtav)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    if(output_vpert.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'vpert   '
      desc_output(n_out) = 'v perturbation (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = va(i,j,k)-v0(i,j,k)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_basestate.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'v0      '
      desc_output(n_out) = 'base-state v (grid-relative)'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

        call   writev( v0 ,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_hadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_hadv'
      if( hadvordrv.eq.3 .or. hadvordrv.eq.5 .or. hadvordrv.eq.7 .or. hadvordrv.eq.9 .or. advwenov.ge.1 )then
        desc_output(n_out) = 'v budget: horizontal advection (non-diff component)'
      else
        desc_output(n_out) = 'v budget: horizontal advection'
      endif
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_hadv)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_vadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_vadv'
      if( vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1 )then
        desc_output(n_out) = 'v budget: vertical advection (non-diff component)'
      else
        desc_output(n_out) = 'v budget: vertical advection'
      endif
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_vadv)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_hidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_hidiff'
      desc_output(n_out) = 'v budget: horiz implicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_hidiff)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_vidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_vidiff'
      desc_output(n_out) = 'v budget: vert implicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_vidiff)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_hediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_hediff'
      desc_output(n_out) = 'v budget: horiz explicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_hediff)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_vediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_vediff'
      desc_output(n_out) = 'v budget: vert explicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_vediff)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_hturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_hturb'
      desc_output(n_out) = 'v budget: horizontal parameterized turbulence'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_hturb)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_vturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_vturb'
      desc_output(n_out) = 'v budget: vertical parameterized turbulence'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_vturb)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_pgrad.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_pgrad'
      desc_output(n_out) = 'v budget: pressure gradient'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_pgrad)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_rdamp.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_rdamp'
      desc_output(n_out) = 'v budget: Rayleigh damper'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_rdamp)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_cor.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_cor'
      desc_output(n_out) = 'v budget: Coriolis acceleration'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_cor)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_cent.ge.1 )then
      if( axisymm.eq.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_cent'
      desc_output(n_out) = 'v budget: centrifugal acceleration'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_cent)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_pbl.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_pbl'
      desc_output(n_out) = 'v budget: PBL scheme'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_pbl)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_vbudget.eq.1 )THEN
      if( vd_lsw.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'vb_lsw'
      desc_output(n_out) = 'v budget: advection by large-scale w'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk
          do j=1,nj+1
          do i=1,ni
            dumv(i,j,k) = vdiag(i,j,k,vd_lsw)
          enddo
          enddo
          enddo
        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    ! arbitrary output (out3d array)

! KLUDGE !
  if(ptype.eq.91919)then
    out3dcheckv:  &
    IF( nout3d.ge.1 .and. ie3d.gt.1 .and. je3d.gt.1 .and. ke3d.gt.1 )THEN

      do n=1,nout3d

        n_out = n_out + 1
        text1 = 'out     '
        if(n.lt.10)then
          write(text1(4:4),211) n
        elseif(n.lt.100)then
          write(text1(4:5),212) n
        elseif(n.lt.1000)then
          write(text1(4:6),213) n
        elseif(n.lt.10000)then
          write(text1(4:7),214) n
        elseif(n.lt.100000)then
          write(text1(4:8),215) n
        else
          print *,'  nout3d is too large '
          call stopcm1
        endif
        name_output(n_out) = text1
        desc_output(n_out) = '3d output'
        unit_output(n_out) = 'unknown'
        grid_output(n_out) = 'v'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts

                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk
                do j=1,nj+1
                do i=1,ni
                  dumv(i,j,k) = out3d(i,j,k,n)
                enddo
                enddo
                enddo

        call   writev(dumv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

      enddo

    ENDIF  out3dcheckv
  endif
! KLUDGE !

  !.............................................

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! w-staggered variables:

  !.............................................

    if(output_w  .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'w       '
      desc_output(n_out) = 'vertical velocity'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew( wa ,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( dotimeavg .and. wtav.ge.1 )THEN
    if(output_w    .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'wtavg'
      desc_output(n_out) = 'time-average w'
      unit_output(n_out) = 'm/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dumw = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = timavg(i,j,k,wtav)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    if( output_tke.eq.1 .and. idoles .and. iusetke )then
      n_out = n_out + 1
      name_output(n_out) = 'tke     '
      desc_output(n_out) = 'subgrid turbulence kinetic energy'
      unit_output(n_out) = 'm^2/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(tkea,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( dotimeavg .and. etav.ge.1 )THEN
    if(output_tke  .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'tketavg'
      desc_output(n_out) = 'time-average tke'
      unit_output(n_out) = 'm^2/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = min(keta,maxk)

          dumw = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(keta,maxk)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = timavg(i,j,k,etav)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    ENDIF

  !.............................................

    if(output_km .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'kmh     '
      IF( ipbl.eq.1 .or. ipbl.eq.3 .or. ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 )THEN
        desc_output(n_out) = 'horizontal eddy viscosity for momentum (from 2D Smagorinsky scheme)'
      ELSE
        desc_output(n_out) = 'horizontal eddy viscosity for momentum'
      ENDIF
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(kmh ,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_km .eq.1)then
    if( sgsmodel.ge.1 .or. ipbl.eq.2 )then
      n_out = n_out + 1
      name_output(n_out) = 'kmv     '
      desc_output(n_out) = 'vertical eddy viscosity for momentum'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(kmv ,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

    endif
    endif

  !.............................................

    if(output_kh .eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'khh     '
      IF( ipbl.eq.1 .or. ipbl.eq.3 .or. ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 )THEN
        desc_output(n_out) = 'horizontal eddy diffusivity for scalars (from 2D Smgorinsky scheme)'
      ELSE
        desc_output(n_out) = 'horizontal eddy diffusivity for scalars'
      ENDIF
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(khh ,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_kh .eq.1)then
    if( sgsmodel.ge.1 .or. ipbl.eq.2 )then
      n_out = n_out + 1
      name_output(n_out) = 'khv     '
      desc_output(n_out) = 'vertical eddy diffusivity for scalars'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(khv ,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
    endif

  !.............................................

  doit = .false.
  if( doit )then
    if( idoles )then
      n_out = n_out + 1
      name_output(n_out) = 'cme'
      desc_output(n_out) = 'Variable C_m in subgrid TKE scheme'
      unit_output(n_out) = 'unitless'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(cme,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif
  endif

  !.............................................

    if( ipbl.eq.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'xkzh'
      desc_output(n_out) = 'eddy diffusivity for heat (from YSU)'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = xkzh(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

!    if( ipbl.eq.1 )then
!      n_out = n_out + 1
!      name_output(n_out) = 'xkzq'
!      desc_output(n_out) = 'eddy diffusivity for moisture (from YSU)'
!      unit_output(n_out) = 'm^2/s'
!      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
!      levs_output(n_out) = maxk+1
!
!          !$omp parallel do default(shared)  &
!          !$omp private(i,j,k)
!          do k=1,maxk+1
!          do j=1,nj
!          do i=1,ni
!            dumw(i,j,k) = xkzq(i,j,k)
!          enddo
!          enddo
!          enddo
!        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
!    endif

  !.............................................

    if( ipbl.eq.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'xkzm'
      desc_output(n_out) = 'eddy viscosity (from YSU)'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = xkzm(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.3 )then
      n_out = n_out + 1
      name_output(n_out) = 'dkt3d'
      desc_output(n_out) = 'Thermal Diffusivity (from GFSEDMF)'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = xkzh(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.3 )then
      n_out = n_out + 1
      name_output(n_out) = 'dku3d'
      desc_output(n_out) = 'Momentum Diffusivity (from GFSEDMF)'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = xkzm(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 )then
      n_out = n_out + 1
      name_output(n_out) = 'exch_h'
      desc_output(n_out) = 'Thermal Diffusivity (from PBL)'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = xkzh(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.4 .or. ipbl.eq.5 .or. ipbl.eq.6 )then
      n_out = n_out + 1
      name_output(n_out) = 'exch_m'
      desc_output(n_out) = 'Momentum Diffusivity (from PBL)'
      unit_output(n_out) = 'm^2/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = xkzm(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.4 .or. ipbl.eq.5 )then
      n_out = n_out + 1
      name_output(n_out) = 'el_pbl'
      desc_output(n_out) = 'Length scale from PBL'
      unit_output(n_out) = 'm'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = el_pbl(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.6 )then
      n_out = n_out + 1
      name_output(n_out) = 'tke_myj'
      desc_output(n_out) = 'tke_myj'
      unit_output(n_out) = 'm^2/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = tke_myj(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if( ipbl.eq.6 )then
      n_out = n_out + 1
      name_output(n_out) = 'el_myj'
      desc_output(n_out) = 'el_myj'
      unit_output(n_out) = 'unk'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = el_myj(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_dissten.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'lenscl'
      desc_output(n_out) = 'lenscl'
      unit_output(n_out) = 'm'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(lenscl,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_dissten.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'dissten '
      desc_output(n_out) = 'dissipation rate'
      unit_output(n_out) = 'm^2/s^3'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(dissten,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_nm.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'nm      '
      desc_output(n_out) = 'squared Brunt-Vaisala frequency'
      unit_output(n_out) = '1/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew( nm ,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_def.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'defv    '
      desc_output(n_out) = 'vertical deformation'
      unit_output(n_out) = '1/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(defv,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_def.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'defh    '
      desc_output(n_out) = 'horizontal deformation'
      unit_output(n_out) = '1/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(defh,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    if(output_def.eq.1)then
      n_out = n_out + 1
      name_output(n_out) = 'epst'
      desc_output(n_out) = 'epst'
      unit_output(n_out) = '1/s^2'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

        call   writew(epst,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
    endif

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_hadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_hadv'
      if( hadvordrv.eq.3 .or. hadvordrv.eq.5 .or. hadvordrv.eq.7 .or. hadvordrv.eq.9 .or. advwenov.ge.1 )then
        desc_output(n_out) = 'w budget: horizontal advection (non-diff component)'
      else
        desc_output(n_out) = 'w budget: horizontal advection'
      endif
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_hadv)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_vadv.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_vadv'
      if( vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1 )then
        desc_output(n_out) = 'w budget: vertical advection (non-diff component)'
      else
        desc_output(n_out) = 'w budget: vertical advection'
      endif
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_vadv)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_hidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_hidiff'
      desc_output(n_out) = 'w budget: horiz implicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_hidiff)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_vidiff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_vidiff'
      desc_output(n_out) = 'w budget: vert implicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_vidiff)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_hediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_hediff'
      desc_output(n_out) = 'w budget: horiz explicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_hediff)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_vediff.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_vediff'
      desc_output(n_out) = 'w budget: vert explicit diffusion'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_vediff)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_hturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_hturb'
      desc_output(n_out) = 'w budget: horizontal parameterized turbulence'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_hturb)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_vturb.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_vturb'
      desc_output(n_out) = 'w budget: vertical parameterized turbulence'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_vturb)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_pgrad.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_pgrad'
      desc_output(n_out) = 'w budget: pressure gradient'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_pgrad)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_rdamp.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_rdamp'
      desc_output(n_out) = 'w budget: Rayleigh damper'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_rdamp)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

    IF( output_wbudget.eq.1 )THEN
      if( wd_buoy.ge.1 )then
      n_out = n_out + 1
      name_output(n_out) = 'wb_buoy'
      desc_output(n_out) = 'w budget: buoyancy'
      unit_output(n_out) = 'm/s/s'
      grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
      levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,maxk+1
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = wdiag(i,j,k,wd_buoy)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      endif
    ENDIF

  !.............................................

      IF( pdcomp )THEN

        n_out = n_out + 1
        name_output(n_out) = 'pgradb'
        desc_output(n_out) = 'vert pres grad: buoyancy component'
        unit_output(n_out) = 'm/s/s'
        grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
          do i=1,ni
            dumw(i,j,1) = 0.0
            dumw(i,j,nk+1) = 0.0
          enddo
          enddo

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,1)-pdiag(i,j,k-1,1))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      ENDIF

  !.............................................


      IF( pdcomp )THEN
        n_out = n_out + 1
        name_output(n_out) = 'pgraddl'
        desc_output(n_out) = 'vert pres grad: linear dynamic component'
        unit_output(n_out) = 'm/s/s'
        grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,2)-pdiag(i,j,k-1,2))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      ENDIF

  !.............................................

      IF( pdcomp )THEN
        n_out = n_out + 1
        name_output(n_out) = 'pgraddn'
        desc_output(n_out) = 'vert pres grad: nonlinear dynamic component'
        unit_output(n_out) = 'm/s/s'
        grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,3)-pdiag(i,j,k-1,3))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      ENDIF

  !.............................................

      IF( pdcomp )THEN
        if( icor.eq.1 )then
        n_out = n_out + 1
        name_output(n_out) = 'pgradc'
        desc_output(n_out) = 'vert pres grad: Coriolis component'
        unit_output(n_out) = 'm/s/s'
        grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        levs_output(n_out) = maxk+1

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,min(nk,maxk+1)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,4)-pdiag(i,j,k-1,4))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
        endif
      ENDIF

  !.............................................

      IF( (cm1setup.eq.1.or.cm1setup.eq.4) .and. sgsmodel.eq.4 .and. dotimeavg )THEN
        n_out = n_out + 1
        name_output(n_out) = 'kmw'
        desc_output(n_out) = 'kmw'
        unit_output(n_out) = 'm2/s'
        grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        levs_output(n_out) = min(ntwk,maxk+1)

          dumw = 0.0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,min(ntwk,maxk+1)
          do j=1,nj
          do i=1,ni
            dumw(i,j,k) = kmwk(i,j,k)
          enddo
          enddo
          enddo
        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))
      ENDIF

  !.............................................

    ! arbitrary output (out3d array)

! KLUDGE !
  if(ptype.eq.91919)then
    out3dcheckw:  &
    IF( nout3d.ge.1 .and. ie3d.gt.1 .and. je3d.gt.1 .and. ke3d.gt.1 )THEN

      do n=1,nout3d

        n_out = n_out + 1
        text1 = 'out     '
        if(n.lt.10)then
          write(text1(4:4),211) n
        elseif(n.lt.100)then
          write(text1(4:5),212) n
        elseif(n.lt.1000)then
          write(text1(4:6),213) n
        elseif(n.lt.10000)then
          write(text1(4:7),214) n
        elseif(n.lt.100000)then
          write(text1(4:8),215) n
        else
          print *,'  nout3d is too large '
          call stopcm1
        endif
        name_output(n_out) = text1
        desc_output(n_out) = '3d output'
        unit_output(n_out) = 'unknown'
        grid_output(n_out) = 'w'     ! s=scalar pts (3d) ; u=u pts (3d) ; v=v pts (3d) ; w=w pts (3d) ; 2=2d scalar pts
        levs_output(n_out) = maxk+1

                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,maxk+1
                do j=1,nj
                do i=1,ni
                  dumw(i,j,k) = out3d(i,j,k,n)
                enddo
                enddo
                enddo

        call   writew(dumw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output(n_out),levs_output(n_out))

      enddo

    ENDIF  out3dcheckw
  endif
! KLUDGE !

  !.............................................

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    ENDDO  bigloop


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!---------------------------------------------------------------
#ifdef MPI
      !  limit to "nlim" writes at a time:
      IF( output_filetype.eq.3 )THEN
      IF( numprocs.gt.nlim )THEN
        doit = .true.
        IF( myid+nlim .le. (numprocs-1) )THEN
          call MPI_ISEND(doit,1,mpi_logical,myid+nlim,999999,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ENDIF
      ENDIF
      ENDIF
#endif
!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Done Writing Data to File '
      if(dowr) write(outfile,*)

    IF(output_format.eq.1)THEN
      if( opens ) close(unit=fnum)
      if( openu ) close(unit=unum)
      if( openv ) close(unit=vnum)
      if( openw ) close(unit=wnum)
#ifdef NETCDF
    ELSEIF( output_format.eq.2 )THEN
      if( opens )then
        if( myid.eq.0 ) print *,'  calling nf90_close ... '
        call disp_err( nf90_close(ncid) , .true. )
        if( myid.eq.0 ) print *,'  ... done '
      endif
#endif
    ENDIF


      sout2d = 0
      sout3d = 0
      u_out = 0
      v_out = 0
      w_out = 0

      do n=1,n_out
        if( grid_output(n).eq.'2' ) sout2d = sout2d+1
        if( grid_output(n).eq.'s' ) sout3d = sout3d+1
        if( grid_output(n).eq.'u' ) u_out = u_out+1
        if( grid_output(n).eq.'v' ) v_out = v_out+1
        if( grid_output(n).eq.'w' ) w_out = w_out+1
      enddo

      s_out = sout2d+sout3d

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  sout2d = ',sout2d
      if(dowr) write(outfile,*) '  sout3d = ',sout3d
      if(dowr) write(outfile,*) '  n_out  = ',n_out
      if(dowr) write(outfile,*) '  s_out  = ',s_out
      if(dowr) write(outfile,*) '  u_out  = ',u_out
      if(dowr) write(outfile,*) '  v_out  = ',v_out
      if(dowr) write(outfile,*) '  w_out  = ',w_out
      if(dowr) write(outfile,*) '  z_out  = ',z_out

      if(dowr) write(outfile,*)

      if( output_format.eq.1 )then

        ! write GrADS descriptor files:
        call write_outputctl(xh,xf,yh,yf,xfref,yfref,sigma,sigmaf,name_output,desc_output,unit_output,grid_output,levs_output,nwrite)

      endif


!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!#endif
      if( myid.eq.0 ) print *,'  ... leaving writeout '
#ifdef MPI
      if(timestats.ge.1)then
        ! this is needed for proper accounting of timing:
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
#endif

      end subroutine writeout


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    ! writeo:
    subroutine writeo(numi,numj,numk1,numk2,nxr,nyr,var,aname,             &
                      myi1p,myi2p,myj1p,myj2p,                             &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fileunit,  &
                      ncid,time_index,output_format,output_filetype,       &
                      dat1,dat2,dat3,reqt,ppnode,d3n,d3t,                  &
                      mynode,nodeleader,nodes,d2i,d2j,d3i,d3j)
#ifdef MPI
    use mpi
#endif
#ifdef NETCDF
    use netcdf
    use writeout_nc_module , only : disp_err
#endif
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing to disk.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode,d3n,d3t,d2i,d2j,d3i,d3j
    real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character(len=*), intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey,fileunit
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(inout) :: irec
    integer, intent(in) :: ncid
    integer, intent(in) :: time_index,output_format,output_filetype
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in) :: mynode,nodeleader,nodes

    integer :: i,j,k
#ifdef NETCDF
    integer :: varid,status
#endif

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !-----------------------------------------------------------------------------
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if( myid.eq.0 ) print *,fileunit,nout,aname

#ifndef MPI

  ! single processor writeo:
      call       writeosp(numi,numj,nxr,nyr,numk1,numk2,ngxy,output_format,fileunit,irec,var,ncid,aname,time_index,dat2(1,1))

#else
  ! MPI section:

  IF(output_filetype.eq.1.or.output_filetype.eq.2)THEN

      IF(myid.ne.nodeleader)THEN
        call     writeocomm1(numi,numj,numk1,numk2,ngxy,d3i,d3j,nodeleader,dat1,var)
      ELSE
        IF(myid.ne.0)THEN
          call   writeocomm2(numi,numj,numk1,numk2,ngxy,d3i,d3j,d3n,d3t,nodeleader,myid,ppnode,reqt,dat1,dat3,var)
        ELSE
          call   writeocomm3(numi,numj,numk1,numk2,ngxy,d2i,d2j,d3i,d3j,d3n,d3t,nodeleader,myid,mynode,ppnode,nodes,numprocs,ni,nj,nxr,nyr,output_format,fileunit,irec,myi1p,myi2p,myj1p,myj2p,reqt,dat1,dat2,dat3,var,ncid,aname,time_index)
        ENDIF
      ENDIF

!!!#ifdef MPI
!!!    ! can help with memory:
!!!    call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!#endif


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !-----   output_filetype = 3   ----------------------------------------------!
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !  this section wites one output file per MPI process:
  !  (for MPI runs only)

  ELSEIF(output_filetype.eq.3)THEN
    IF( output_format.eq.1 )THEN
      ! grads format:
      DO k=numk1,numk2
        write(fileunit,rec=irec) ((var(i,j,k),i=1,numi),j=1,numj)
        irec=irec+1
      ENDDO
#ifdef NETCDF
    ELSEIF( output_format.eq.2 )THEN
      ! netcdf format:
      status = nf90_inq_varid(ncid,aname,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1c in writeo, aname = ',aname
        print *,nf90_strerror(status)
        call stopcm1
      endif
      print *,'  ncid,aname,varid = ',ncid,aname,varid,var(1,1,1)
      DO k=numk1,numk2
        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j)=var(i,j,k)
        enddo
        enddo
        if(numk1.eq.numk2)then
          status = nf90_put_var(ncid,varid,dat1,(/1,1,time_index/),(/numi,numj,1/))
          if(status.ne.nf90_noerr)then
            print *,'  Error2c in writeo, aname = ',aname
            print *,'  ncid,varid,time_index = ',ncid,varid,time_index
            print *,nf90_strerror(status)
            call stopcm1
          endif
        else
          status = nf90_put_var(ncid,varid,dat1,(/1,1,k,time_index/),(/numi,numj,1,1/))
          if(status.ne.nf90_noerr)then
            print *,'  Error3c in writeo, aname = ',aname
            print *,'  ncid,varid,time_index = ',ncid,varid,time_index
            print *,nf90_strerror(status)
            call stopcm1
          endif
        endif
      ENDDO
#endif
    ENDIF
  ENDIF
#endif

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !-----------------------------------------------------------------------------
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      end subroutine writeo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifndef MPI


      subroutine writeosp(numi,numj,nxr,nyr,numk1,numk2,ngxy,output_format,fileunit,irec,var,ncid,aname,time_index,dat2)
#ifdef NETCDF
      use netcdf
      use writeout_nc_module , only : disp_err
#endif
      implicit none

      integer, intent(in) :: numi,numj,nxr,nyr,numk1,numk2,ngxy,output_format,fileunit
      integer, intent(inout) :: irec
      real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: aname
      integer, intent(in) :: time_index
      real, intent(inout), dimension(nxr,nyr) :: dat2

      integer :: i,j,k
#ifdef NETCDF
      integer :: varid,status
#endif
#ifdef DP
      real*4, dimension(nxr,nyr) :: vout
#endif

      ! writeo for single processor:

#ifdef NETCDF
      if( output_format.eq.2 )then
        status = nf90_inq_varid(ncid,aname,varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error1a in writeo, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
        print *,'  ncid,aname,varid = ',ncid,aname,varid,var(1,1,1)
      endif
#endif

      IF(output_format.eq.1)THEN

        ! ----- grads format -----
        DO k=numk1,numk2
#ifdef DP
          ! double precision:
          do j=1,nyr
          do i=1,nxr
            vout(i,j) = var(i,j,k)
          enddo
          enddo
          write(fileunit,rec=irec) ((vout(i,j),i=1,nxr),j=1,nyr)
#else
          ! normal:
          write(fileunit,rec=irec) ((var(i,j,k),i=1,nxr),j=1,nyr)
#endif
          irec = irec+1
        ENDDO

#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN

        ! ----- netcdf format -----
        if(numk1.eq.numk2)then
          do j=1,nyr
          do i=1,nxr
            dat2(i,j) = var(i,j,1)
          enddo
          enddo
          status = nf90_put_var(ncid,varid,dat2,(/1,1,time_index/),(/nxr,nyr,1/))
          if(status.ne.nf90_noerr)then
            print *,'  Error2 in writeo, aname = ',aname
            print *,'  ncid,varid,time_index = ',ncid,varid,time_index
            print *,nf90_strerror(status)
            call stopcm1
          endif
        else
          DO k=numk1,numk2
            do j=1,nyr
            do i=1,nxr
              dat2(i,j) = var(i,j,k)
            enddo
            enddo
            status = nf90_put_var(ncid,varid,dat2,(/1,1,k,time_index/),(/nxr,nyr,1,1/))
            if(status.ne.nf90_noerr)then
              print *,'  Error3 in writeo, aname = ',aname
              print *,'  ncid,varid,time_index = ',ncid,varid,time_index
              print *,nf90_strerror(status)
              call stopcm1
            endif
          ENDDO
        endif
#endif

      ENDIF


      end subroutine writeosp


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


      subroutine writeocomm1(numi,numj,numk1,numk2,ngxy,d3i,d3j,nodeleader,dat1,var)
      use mpi
      implicit none

      integer, intent(in) :: numi,numj,numk1,numk2,ngxy,d3i,d3j,nodeleader
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

      integer :: i,j,k,tag,reqs,ierr

      ! ordinary processor ... send data to nodeleader:

      tag = 1

      DO k=numk1,numk2

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j)=var(i,j,k)
        enddo
        enddo

        call MPI_ISEND(dat1(1,1),d3i*d3j,MPI_REAL,nodeleader,tag,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)

        tag = tag+2

      ENDDO

      ! DONE, ordinary processors

      end subroutine writeocomm1


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


      subroutine writeocomm2(numi,numj,numk1,numk2,ngxy,d3i,d3j,d3n,d3t,nodeleader,myid,ppnode,reqt,dat1,dat3,var)
      use mpi
      implicit none

      integer, intent(in) :: numi,numj,numk1,numk2,ngxy,d3i,d3j,d3n,d3t,nodeleader,myid,ppnode
      integer, intent(inout), dimension(d3t) :: reqt
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
      real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var

      integer :: i,j,k,tag,reqs,ierr,proc
      integer, dimension(mpi_status_size,ppnode-1) :: status

      ! begin nodeleader section (not proc 0):

      tag = 1

      DO k=numk1,numk2

        ! start receives from all other processors on a node:
        do proc=myid+1,myid+(ppnode-1)
          call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
        enddo

        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat3(i,j,myid)=var(i,j,k)
        enddo
        enddo

        ! wait for receives to finish:
        call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status,ierr)

        ! send data to processor 0:
        call MPI_ISEND(dat3(1,1,myid),d3i*d3j*ppnode,MPI_REAL,0,tag+1,MPI_COMM_WORLD,reqs,ierr)

        ! wait for send to finish:
        call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)

        tag = tag+2

      ENDDO

      ! end nodeleader section (not proc0):

      end subroutine writeocomm2


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPI


      subroutine writeocomm3(numi,numj,numk1,numk2,ngxy,d2i,d2j,d3i,d3j,d3n,d3t,nodeleader,myid,mynode,ppnode,nodes,numprocs,ni,nj,nxr,nyr,output_format,fileunit,irec,myi1p,myi2p,myj1p,myj2p,reqt,dat1,dat2,dat3,var,ncid,aname,time_index)
      use mpi
#ifdef NETCDF
      use netcdf
#endif
      implicit none

      integer, intent(in) :: numi,numj,numk1,numk2,ngxy,d2i,d2j,d3i,d3j,d3n,d3t,nodeleader,myid,mynode,ppnode,nodes,numprocs,ni,nj,nxr,nyr,output_format,fileunit
      integer, intent(inout) :: irec
      integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
      integer, intent(inout), dimension(d3t) :: reqt
      real, intent(inout), dimension(d3i,d3j) :: dat1
      real, intent(inout), dimension(d2i,d2j) :: dat2
      real, intent(inout), dimension(d3i,d3j,0:d3n-1) :: dat3
      real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: aname
      integer, intent(in) :: time_index

      integer :: i,j,k,tag,reqs,ierr,proc
      integer :: index2,fooi,fooj,nn,nnn,ntot,n1,n2
      integer :: index,nitmp,njtmp
#ifdef DP
      real*4, dimension(nxr,nyr) :: vout
#endif
#ifdef NETCDF
      integer :: varid,status
#endif

      ! begin proc 0:

      tag = 1

      ! start receives from all other processors on a node:
      do proc=1,(ppnode-1)
        call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc),ierr)
      enddo

      ! start receives from other nodeleaders:
      do nn = 1,(nodes-1)
        proc = nn*ppnode
        call MPI_IRECV(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+1,MPI_COMM_WORLD,reqt(ppnode-1+nn),ierr)
      enddo

#ifdef NETCDF
      if( output_format.eq.2 )then
        ! while we are waiting, get varid  (k=numk1 only)
        status = nf90_inq_varid(ncid,aname,varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error1b in writeo, aname = ',aname
          print *,nf90_strerror(status)
          call stopcm1
        endif
      endif
#endif

      DO k=numk1,numk2

        ! my data:
        !$omp parallel do default(shared)  &
        !$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat2(i,j)=var(i,j,k)
        enddo
        enddo

          ! wait for data to arrive:
          ntot = ppnode-1 + nodes-1 
          do nn=1,ntot
            call mpi_waitany(ntot,reqt(1:ntot),index,MPI_STATUS_IGNORE,ierr)
            if( index.le.(ppnode-1) )then
              ! data from ordinary procs on node:
              proc = index
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

          ! ready to write ... but first, start receives for next level:
          IF( k.lt.numk2 )THEN
            do proc=1,(ppnode-1)
              call MPI_IRECV(dat3(1,1,proc),d3i*d3j,MPI_REAL,proc,tag+2,MPI_COMM_WORLD,reqt(proc),ierr)
            enddo
            do nn = 1,(nodes-1)
              proc = nn*ppnode
              call MPI_IRECV(dat3(1,1,proc),d3i*d3j*ppnode,MPI_REAL,proc,tag+3,MPI_COMM_WORLD,reqt(ppnode-1+nn),ierr)
            enddo
          ENDIF

      !-------------------- write data --------------------!
        IF(output_format.eq.1)THEN
          ! ----- grads format -----
          write(fileunit,rec=irec) ((dat2(i,j),i=1,nxr),j=1,nyr)
          irec = irec+1
#ifdef NETCDF
        ELSEIF(output_format.eq.2)THEN
          ! ----- netcdf format -----
          if(numk1.eq.numk2)then
            status = nf90_put_var(ncid,varid,dat2,(/1,1,time_index/),(/nxr,nyr,1/))
          else
            status = nf90_put_var(ncid,varid,dat2,(/1,1,k,time_index/),(/nxr,nyr,1,1/))
          endif
          if(status.ne.nf90_noerr)then
            print *,'  Error2 in writeo, aname = ',aname
            print *,'  ncid,varid,time_index = ',ncid,varid,time_index
            print *,nf90_strerror(status)
            call stopcm1
          endif
#endif
        ENDIF
      !-------------------- end write data --------------------!

        tag = tag+2

      ENDDO

      ! end proc0:

      end subroutine writeocomm3


#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_outputctl(xh,xf,yh,yf,xfref,yfref,sigma,sigmaf,name_output,desc_output,unit_output,grid_output,levs_output,nwrite)
      use input
      use constants , only : grads_undef
      implicit none

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      character(len=60), intent(in), dimension(maxvars) :: desc_output
      character(len=40), intent(in), dimension(maxvars) :: name_output,unit_output
      character(len=1),  intent(in), dimension(maxvars) :: grid_output
      integer,           intent(in), dimension(maxvars) :: levs_output
      integer, intent(in) ::  nwrite

      integer :: i,j,k,n,nn,n1,n2,ctl,ctlmax
      character(len=maxstring) :: dstring
      character(len=12) :: a12
      logical :: doit

      !----------------------------------------------------------------
      ! This subroutine writes the GrADS descriptor file for 3d output
      !----------------------------------------------------------------

    idcheck:  &
    IF( myid.eq.0 )THEN

      do i=1,maxstring
        string(i:i) = ' '
        dstring(i:i) = ' '
      enddo

      ctlmax = 4
      if( output_interp.eq.1 ) ctlmax = 5

      ctlloop:  &
      DO ctl = 1 , ctlmax

        doit = .false.

        IF(     ctl.eq.1 )THEN
          if( s_out.ge.1 )then
          string = 'cm1out_s.ctl'
          if(dowr) write(outfile,*) string
          if(output_filetype.eq.1)then
            dstring = 'cm1out_s.dat'
          elseif(output_filetype.ge.2)then
            dstring = 'cm1out_00%y4_s.dat'
          endif
          doit = .true.
          endif
        ELSEIF( ctl.eq.2 )THEN
          if( u_out.ge.1 )then
          string = 'cm1out_u.ctl'
          if(dowr) write(outfile,*) string
          if(output_filetype.eq.1)then
            dstring = 'cm1out_u.dat'
          elseif(output_filetype.ge.2)then
            dstring = 'cm1out_00%y4_u.dat'
          endif
          doit = .true.
          endif
        ELSEIF( ctl.eq.3 )THEN
          if( v_out.ge.1 )then
          string = 'cm1out_v.ctl'
          if(dowr) write(outfile,*) string
          if(output_filetype.eq.1)then
            dstring = 'cm1out_v.dat'
          elseif(output_filetype.ge.2)then
            dstring = 'cm1out_00%y4_v.dat'
          endif
          doit = .true.
          endif
        ELSEIF( ctl.eq.4 )THEN
          if( w_out.ge.1 )then
          string = 'cm1out_w.ctl'
          if(dowr) write(outfile,*) string
          if(output_filetype.eq.1)then
            dstring = 'cm1out_w.dat'
          elseif(output_filetype.ge.2)then
            dstring = 'cm1out_00%y4_w.dat'
          endif
          doit = .true.
          endif
        ELSEIF( ctl.eq.5 )THEN
          if( s_out.ge.1 )then
          string = 'cm1out_i.ctl'
          if(dowr) write(outfile,*) string
          if(output_filetype.eq.1)then
            dstring = 'cm1out_i.dat'
          elseif(output_filetype.ge.2)then
            dstring = 'cm1out_00%y4_i.dat'
          endif
          doit = .true.
          endif
        ELSE
          print *,'  98371 '
          call stopcm1
        ENDIF

        dowrite:  &
        IF( doit )THEN

          open(unit=50,file=string,status='unknown')

          write(50,201) dstring
          if(output_filetype.ge.2) write(50,221)
          if( outunits.eq.2 )then
            write(50,212) trim(cm1version)
          else
            write(50,202) trim(cm1version)
          endif
          write(50,203) grads_undef

          IF( ctl.eq.2 )THEN
            ! u staggering:
            if(stretch_x.ge.1)then
              write(50,214) nx+1
              do i=1,nx+1
                write(50,217) outunitconv*xfref(i)
              enddo
            else
              write(50,204) nx+1,outunitconv*xf(1),outunitconv*dx
            endif
          ELSE
            ! s staggering:
            if(stretch_x.ge.1)then
              write(50,214) nx
              do i=1,nx
                write(50,217) outunitconv*( 0.5*(xfref(i)+xfref(i+1)) )
              enddo
            else
              write(50,204) nx,outunitconv*xh(1),outunitconv*dx
            endif
          ENDIF

          IF( ctl.eq.3 )THEN
            ! v staggering:
            if(stretch_y.ge.1)then
              write(50,215) ny+1
              do j=1,ny+1
                write(50,217) outunitconv*yfref(j)
              enddo
            else
              write(50,205) ny+1,outunitconv*yf(1),outunitconv*dy
            endif
          ELSE
            ! s staggering:
            if(stretch_y.ge.1)then
              write(50,215) ny
              do j=1,ny
                write(50,217) outunitconv*( 0.5*(yfref(j)+yfref(j+1)) )
              enddo
            else
              write(50,205) ny,outunitconv*yh(1),outunitconv*dy
            endif
          ENDIF

          IF( ctl.eq.4 )THEN
            ! w staggering:
            if(stretch_z.eq.0)then
              write(50,206) maxk+1,0.0,outunitconv*dz
            else
              write(50,216) maxk+1
              do k=1,maxk+1
                write(50,217) outunitconv*sigmaf(k)
              enddo
            endif
          ELSE
            ! s staggering:
            if(stretch_z.eq.0)then
              write(50,206) maxk,outunitconv*sigma(1),outunitconv*dz
            else
              write(50,216) maxk
              do k=1,maxk
                write(50,217) outunitconv*sigma(k)
              enddo
            endif
          ENDIF


              write(50,227) nwrite


          IF( ctl.eq.1 .or. ctl.eq.5 )THEN
            ! scalars:
            write(50,208) s_out
            n1 = 1
            n2 = s_out
          ELSEIF( ctl.eq.2 )THEN
            ! u vars:
            write(50,208) u_out
            n1 = s_out+1
            n2 = s_out+u_out
          ELSEIF( ctl.eq.3 )THEN
            ! v vars:
            write(50,208) v_out
            n1 = s_out+u_out+1
            n2 = s_out+u_out+v_out
          ELSEIF( ctl.eq.4 )THEN
            ! w vars:
            write(50,208) w_out
            n1 = s_out+u_out+v_out+1
            n2 = s_out+u_out+v_out+w_out
          ENDIF

          do n = n1,n2
            a12 = '            '
            nn = len(trim(unit_output(n)))
            nn = min( nn , 10 )
            write(a12(2:11),314) unit_output(n)
            write(a12(1:1),301 )       '('
            write(a12(nn+2:nn+2),301 ) ')'
            ! account for both 2d and 3d output files:
            if(     grid_output(n).eq.'2' )then
              write(50,209) name_output(n),   0  ,desc_output(n),a12
            elseif( grid_output(n).eq.'s' .or. grid_output(n).eq.'u' .or. grid_output(n).eq.'v' )then
              write(50,209) name_output(n),levs_output(n),desc_output(n),a12
            elseif( grid_output(n).eq.'w' )then
              write(50,209) name_output(n),levs_output(n),desc_output(n),a12
            else
              print *,'  98371 '
              call stopcm1
            endif
          enddo

          write(50,210)
          close(unit=50)

        ENDIF  dowrite

      ENDDO  ctlloop

    ENDIF  idcheck

301   format(a1)
314   format(a10)

201   format('dset ^',a)
202   format('title CM1 output, using version ',a,'; units of x,y,z are km; time is generic, see cm1out_metadata for actual times')
212   format('title CM1 output, using version ',a,'; units of x,y,z are meters; time is generic, see cm1out_metadata for actual times')
221   format('options template')
203   format('undef ',f10.1)
!204   format('xdef ',i6,' linear ',f13.6,1x,f13.6)
204   format('xdef ',i6,' linear ',es14.6e2,1x,es14.6e2)
214   format('xdef ',i6,' levels ')
!205   format('ydef ',i6,' linear ',f13.6,1x,f13.6)
205   format('ydef ',i6,' linear ',es14.6e2,1x,es14.6e2)
215   format('ydef ',i6,' levels ')
!206   format('zdef ',i6,' linear ',f13.6,1x,f13.6)
206   format('zdef ',i6,' linear ',es14.6e2,1x,es14.6e2)
216   format('zdef ',i6,' levels ')
!217   format(2x,f13.6)
217   format(2x,es14.6e2)
227   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
208   format('vars ',i6)
209   format(a12,1x,i6,' 99 ',a60,1x,a12)
210   format('endvars')

      end subroutine write_outputctl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_grads_metadata(nwrite,mrec,nstep,mtime,dt,adt,ndt,domainlocx,domainlocy)
      use input
      use constants , only : grads_undef
      implicit none

      integer, intent(in) :: nwrite
      integer, intent(inout) :: mrec
      integer, intent(in) :: nstep
      double precision, intent(in) :: mtime
      real, intent(in) :: dt
      double precision, intent(in) :: adt
      integer, intent(in) :: ndt
      double precision, intent(in) :: domainlocx,domainlocy

      integer :: i,n,nvar
      real :: var
      character(len=80) :: a1,a2
      character(len=80), dimension(varmax) :: varname,vardesc

      nvar = 0

    IF( myid.eq.0 )THEN

      open(unit=61,file='cm1out_metadata.dat',status='unknown',  &
           form='unformatted',access='direct',recl=4)

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'mtime'
      vardesc(nvar) = 'model time (seconds since beginning of simulation)'

      write(61,rec=mrec) sngl(mtime)
      mrec = mrec+1

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'nstep'
      vardesc(nvar) = 'time steps (since beginning of simulation)'

      write(61,rec=mrec) float(nstep)
      mrec = mrec+1

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'dt'
      vardesc(nvar) = 'current time step (seconds)'

      write(61,rec=mrec) dt
      mrec = mrec+1

    !---------------------------

    IF( adapt_dt .eq. 1 )THEN

      nvar = nvar+1
      varname(nvar) = 'adt'
      vardesc(nvar) = 'average time step (seconds) since last call to statpack'

      var = adt/float(max(1,ndt))
      write(61,rec=mrec) var
      mrec = mrec+1

    ENDIF

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'nwrite'
      vardesc(nvar) = 'writeout number'

      write(61,rec=mrec) float(nwrite)
      mrec = mrec+1

    !---------------------------

    IF( imove.ge.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'umove'
      vardesc(nvar) = 'umove (m/s)'

      write(61,rec=mrec) umove
      mrec = mrec+1

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'vmove'
      vardesc(nvar) = 'vmove (m/s)'

      write(61,rec=mrec) vmove
      mrec = mrec+1

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'domainlocx'
      vardesc(nvar) = 'x location of (center of) domain (m)'

      write(61,rec=mrec) sngl(domainlocx)
      mrec = mrec+1

    !---------------------------

      nvar = nvar+1
      varname(nvar) = 'domainlocy'
      vardesc(nvar) = 'y location of (center of) domain (m)'

      write(61,rec=mrec) sngl(domainlocy)
      mrec = mrec+1

    ENDIF

    !---------------------------

    IF( radopt.ge.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'solcon'
      vardesc(nvar) = 'solar constant (W/m^2)'

      write(61,rec=mrec) rad_solcon
      mrec = mrec+1

    ENDIF

    !---------------------------

    IF( radopt.ge.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'jday'
      vardesc(nvar) = 'Julian day (for radiation calculations)'

      write(61,rec=mrec) rad_jday
      mrec = mrec+1

    ENDIF

    !---------------------------

    IF( radopt.eq.2 )THEN

      nvar = nvar+1
      varname(nvar) = 'zenith_angle'
      vardesc(nvar) = 'solar zenith angle  (radians)'

      write(61,rec=mrec) rad_zenangle
      mrec = mrec+1

    ENDIF

    !---------------------------

    IF( radopt.eq.2 )THEN

      nvar = nvar+1
      varname(nvar) = 'declin_angle'
      vardesc(nvar) = 'solar declination (radians)'

      write(61,rec=mrec) rad_declin
      mrec = mrec+1

    ENDIF

    !---------------------------

    IF( radopt.eq.2 )THEN

      nvar = nvar+1
      varname(nvar) = 'hrang'
      vardesc(nvar) = 'solar hour angle (radians)'

      write(61,rec=mrec) rad_hrang
      mrec = mrec+1

    ENDIF

    !---------------------------

      close(unit=61)

      open(unit=50,file='cm1out_metadata.ctl')
      write(50,201)
      write(50,202) trim(cm1version)
      write(50,203) grads_undef
      write(50,204)
      write(50,205)
      write(50,206)
      write(50,227) nwrite
      write(50,208) nvar
      do n=1,nvar
        do i=1,80
          a1(i:i) = ' '
          a2(i:i) = ' '
        enddo
        a1 = varname(n)
        a2 = vardesc(n)
        write(50,209) a1(1:12),0,a2(1:60)
      enddo
      write(50,210)
      close(unit=50)

201   format('dset ^cm1out_metadata.dat')
202   format('title CM1 metadata for 3d output files; using version ',a)
203   format('undef ',f10.1)
204   format('xdef 1 linear 0 1')
205   format('ydef 1 linear 0 1')
206   format('zdef 1 linear 0 1')
227   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
208   format('vars ',i6)
209   format(a12,1x,i6,' 99 ',a60,1x,a12)
210   format('endvars')

    ENDIF

    !---------------------------

      end subroutine write_grads_metadata


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine open_file(type,fnum,nwrite)
    use input
    implicit none

    character(len=1), intent(in) :: type
    integer, intent(in) :: fnum,nwrite

    integer :: i

    do i=1,maxstring
      string(i:i) = ' '
    enddo

100 format(i6.6)
102 format(i6.6)

  !---------------------------------------------------------------------

    checktype:  &
    IF(     type .eq. 's' )THEN

      sof1:  &
      IF(output_format.eq.1)THEN
        ! grads output:

!beh1
        IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN

          IF( myid.eq.nodeleader )THEN

          if( output_filetype.eq.1 )then

            ! one file:
            if(fnum.eq.51)then
              string = 'cm1out_s.dat'
            elseif(fnum.eq.71)then
              string = 'cm1out_i.dat'
            endif

          elseif( output_filetype.eq.2 )then
!beh2

            ! one file per output time:
            if(fnum.eq.51)then
              string = 'cm1out_XXXXXX_s.dat'
            elseif(fnum.eq.71)then
              string = 'cm1out_XXXXXX_i.dat'
            endif
            write(string(8:13),102) nwrite

          endif

          if(dowr) write(outfile,*) '  Opening ',trim(string)
          open(unit=fnum,file=string,form='unformatted',access='direct',   &
               recl=(nx*ny*4),status='unknown')
          opens = .true.
!beh3

          ENDIF

        ELSEIF( output_filetype.eq.3 )THEN

          if(fnum.eq.51)then
            string = 'cm1out_XXXXXX_YYYYYY_s.dat  '
          elseif(fnum.eq.71)then
            string = 'cm1out_XXXXXX_YYYYYY_i.dat  '
          endif
          write(string( 8:13),100) myid
          write(string(15:20),100) nwrite

          if(dowr) write(outfile,*) '  myid,string=',myid,'   ',string
          open(unit=fnum,file=string,                  &
               form='unformatted',access='direct',   &
               recl=(ni*nj*4),status='unknown')
          opens = .true.

        ELSE

          print *,'  unrecognized value for output_filetype in open_file_s '
          call stopcm1

        ENDIF

      ENDIF  sof1

  !---------------------------------------------------------------------

    ELSEIF( type .eq. 'u' )THEN

      uof1:  &
      IF(output_format.eq.1)THEN
        ! grads output:

!beh1
        IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN

          IF( myid.eq.nodeleader )THEN

          if( output_filetype.eq.1 )then

            ! one file:
            string = 'cm1out_u.dat'





          elseif( output_filetype.eq.2 )then
!beh2

            ! one file per output time:
            string = 'cm1out_XXXXXX_u.dat'
            write(string(8:13),102) nwrite





          endif

          if(dowr) write(outfile,*) '  Opening ',trim(string)
          open(unit=fnum,file=string,form='unformatted',access='direct',   &
               recl=((nx+1)*ny*4),status='unknown')
          openu = .true.
!beh3

          ENDIF

        ELSEIF( output_filetype.eq.3 )THEN


            string = 'cm1out_XXXXXX_YYYYYY_u.dat  '



          write(string( 8:13),100) myid
          write(string(15:20),100) nwrite

          if(dowr) write(outfile,*) '  myid,string=',myid,'   ',string
          open(unit=fnum,file=string,                  &
               form='unformatted',access='direct',   &
               recl=((ni+1)*nj*4),status='unknown')
          openu = .true.

        ELSE

          print *,'  unrecognized value for output_filetype in open_file_u '
          call stopcm1

        ENDIF

      ENDIF  uof1

  !---------------------------------------------------------------------

    ELSEIF( type .eq. 'v' )THEN

      vof1:  &
      IF(output_format.eq.1)THEN
        ! grads output:

!beh1
        IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN

          IF( myid.eq.nodeleader )THEN

          if( output_filetype.eq.1 )then

            ! one file:
            string = 'cm1out_v.dat'





          elseif( output_filetype.eq.2 )then
!beh2

            ! one file per output time:
            string = 'cm1out_XXXXXX_v.dat'
            write(string(8:13),102) nwrite





          endif

          if(dowr) write(outfile,*) '  Opening ',trim(string)
          open(unit=fnum,file=string,form='unformatted',access='direct',   &
               recl=(nx*(ny+1)*4),status='unknown')
          openv = .true.
!beh3

          ENDIF

        ELSEIF( output_filetype.eq.3 )THEN


            string = 'cm1out_XXXXXX_YYYYYY_v.dat  '



          write(string( 8:13),100) myid
          write(string(15:20),100) nwrite

          if(dowr) write(outfile,*) '  myid,string=',myid,'   ',string
          open(unit=fnum,file=string,                  &
               form='unformatted',access='direct',   &
               recl=(ni*(nj+1)*4),status='unknown')
          openv = .true.

        ELSE

          print *,'  unrecognized value for output_filetype in open_file_v '
          call stopcm1

        ENDIF

      ENDIF  vof1

  !---------------------------------------------------------------------

    ELSEIF( type .eq. 'w' )THEN

      wof1:  &
      IF(output_format.eq.1)THEN
        ! grads output:

!beh1
        IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN

          IF( myid.eq.nodeleader )THEN

          if( output_filetype.eq.1 )then

            ! one file:
            string = 'cm1out_w.dat'





          elseif( output_filetype.eq.2 )then
!beh2

            ! one file per output time:
            string = 'cm1out_XXXXXX_w.dat'
            write(string(8:13),102) nwrite





          endif

          if(dowr) write(outfile,*) '  Opening ',trim(string)
          open(unit=fnum,file=string,form='unformatted',access='direct',   &
               recl=(nx*ny*4),status='unknown')
          openw = .true.
!beh3

          ENDIF

        ELSEIF( output_filetype.eq.3 )THEN


            string = 'cm1out_XXXXXX_YYYYYY_w.dat  '



          write(string( 8:13),100) myid
          write(string(15:20),100) nwrite

          if(dowr) write(outfile,*) '  myid,string=',myid,'   ',string
          open(unit=fnum,file=string,                  &
               form='unformatted',access='direct',   &
               recl=(ni*nj*4),status='unknown')
          openw = .true.

        ELSE

          print *,'  unrecognized value for output_filetype in open_file_w '
          call stopcm1

        ENDIF

      ENDIF  wof1

  !---------------------------------------------------------------------

    ELSE

      print *,'  Unknown type in open_file '
      call stopcm1

    ENDIF  checktype

    end subroutine open_file


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine write2d(  var2d     ,snum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output,levs_output)
    use input
    implicit none

    real, intent(in), dimension(ib:ie,jb:je) :: var2d
    integer, intent(in) :: snum,nwrite,nloop
    integer, intent(inout) :: srec
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(inout) :: ncid
    integer, intent(in) :: time_index
    character(len=40), intent(in) :: name_output
    integer, intent(in) :: levs_output

    nout = n_out  ! kludge
    IF( .not. opens .and. output_format.eq.1 )THEN
      call open_file('s',snum,nwrite)
    ENDIF

    if( output_format.eq.1 .or. ( output_format.eq.2 .and. nloop.eq.2 ) )   &
    call writeo(ni,nj,1,1,nx,ny,var2d(ib,jb),trim(name_output),             &
                myi1p,myi2p,myj1p,myj2p,                                    &
                ni,nj,ngxy,myid,numprocs,nodex,nodey,srec,snum,             &
                ncid,time_index,output_format,output_filetype,              &
                dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,        &
                mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

    end subroutine write2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writes(vars,snum,srec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output,levs_output,sigma,zs,zh,dum9)
    use input
    use misclibs, only : zinterp
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: vars
    integer, intent(in) :: snum,nwrite,nloop
    integer, intent(inout) :: srec
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(inout) :: ncid
    integer, intent(in) :: time_index
    character(len=40), intent(in) :: name_output
    real, intent(in), dimension(kb:ke) :: sigma
    real, intent(in), dimension(ib:ie,jb:je) :: zs
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum9
    integer, intent(in) :: levs_output

    nout = n_out  ! kludge
    IF( .not. opens .and. output_format.eq.1 )THEN
      call open_file('s',snum,nwrite)
    ENDIF

    if( dointerp .and. trim(name_output).ne.'zhval' ) call zinterp(sigma,zs,zh,vars,dum9)

    if( output_format.eq.1 .or. ( output_format.eq.2 .and. nloop.eq.2 ) )   &
    call writeo(ni,nj,1,levs_output,nx,ny,vars(ib,jb,1),trim(name_output),         &
                myi1p,myi2p,myj1p,myj2p,                                    &
                ni,nj,ngxy,myid,numprocs,nodex,nodey,srec,snum,             &
                ncid,time_index,output_format,output_filetype,              &
                dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,        &
                mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)

    end subroutine writes


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writeu(varu,unum,urec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output,levs_output)
    use input
    implicit none

    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: varu
    integer, intent(in) :: unum,nwrite,nloop
    integer, intent(inout) :: urec
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(inout) :: ncid
    integer, intent(in) :: time_index
    character(len=40), intent(in) :: name_output
    integer, intent(in) :: levs_output

    nout = n_out  ! kludge
  if( .not. dointerp )then
    IF( .not. openu .and. output_format.eq.1 )THEN
      call open_file('u',unum,nwrite)
    ENDIF

    if( output_format.eq.1 .or. ( output_format.eq.2 .and. nloop.eq.2 ) )   &
    call writeo(ni+1,nj,1,levs_output,nx+1,ny,varu(ib,jb,1),trim(name_output),     &
                myi1p,myi2p,myj1p,myj2p,                                    &
                ni,nj,ngxy,myid,numprocs,nodex,nodey,urec,unum,             &
                ncid,time_index,output_format,output_filetype,              &
                dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,        &
                mynode,nodeleader,nodes,d2iu,d2ju,d3iu,d3ju)
  endif

    end subroutine writeu


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writev(varv,vnum,vrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output,levs_output)
    use input
    implicit none

    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: varv
    integer, intent(in) :: vnum,nwrite,nloop
    integer, intent(inout) :: vrec
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(inout) :: ncid
    integer, intent(in) :: time_index
    character(len=40), intent(in) :: name_output
    integer, intent(in) :: levs_output

    nout = n_out  ! kludge
  if( .not. dointerp )then
    IF( .not. openv .and. output_format.eq.1 )THEN
      call open_file('v',vnum,nwrite)
    ENDIF

    if( output_format.eq.1 .or. ( output_format.eq.2 .and. nloop.eq.2 ) )   &
    call writeo(ni,nj+1,1,levs_output,nx,ny+1,varv(ib,jb,1),trim(name_output),     &
                myi1p,myi2p,myj1p,myj2p,                                    &
                ni,nj,ngxy,myid,numprocs,nodex,nodey,vrec,vnum,             &
                ncid,time_index,output_format,output_filetype,              &
                dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,        &
                mynode,nodeleader,nodes,d2iv,d2jv,d3iv,d3jv)
  endif

    end subroutine writev


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writew(varw,wnum,wrec,nwrite,nloop,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p,ncid,time_index,name_output,levs_output)
    use input
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: varw
    integer, intent(in) :: wnum,nwrite,nloop
    integer, intent(inout) :: wrec
    real, intent(inout), dimension(d3i,d3j) :: dat1
    real, intent(inout), dimension(d2i,d2j) :: dat2
    real, intent(inout), dimension(d3i,d3j,d3n) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(inout) :: ncid
    integer, intent(in) :: time_index
    character(len=40), intent(in) :: name_output
    integer, intent(in) :: levs_output

    nout = n_out  ! kludge
  if( .not. dointerp )then
    IF( .not. openw .and. output_format.eq.1 )THEN
      call open_file('w',wnum,nwrite)
    ENDIF

    if( output_format.eq.1 .or. ( output_format.eq.2 .and. nloop.eq.2 ) )   &
    call writeo(ni,nj,1,levs_output,nx,ny,varw(ib,jb,1),trim(name_output),       &
                myi1p,myi2p,myj1p,myj2p,                                    &
                ni,nj,ngxy,myid,numprocs,nodex,nodey,wrec,wnum,             &
                ncid,time_index,output_format,output_filetype,              &
                dat1(1,1),dat2(1,1),dat3(1,1,1),reqt,ppnode,d3n,d3t,        &
                mynode,nodeleader,nodes,d2is,d2js,d3is,d3js)
  endif

    end subroutine writew


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  END MODULE writeout_module
