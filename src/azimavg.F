  MODULE azimavg_module

  implicit none

  private
  public :: azimavg,getpsmth,getpcent,update_adapt_move,getrangle

  integer, parameter :: varmax = 1000

  integer, parameter :: fnums = 67
  integer, parameter :: fnumw = 68

  integer :: ncid,wloop

  CONTAINS

    subroutine azimavg(nstep,mtime,nwritea,arecs,arecw,qname,dt,dosfcflx,      &
                   icrs,icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy, &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,            &
                   yh,vh,rvh,yf,vf,rvf,                                        &
                   xfref,yfref,xhref,yhref,rds,sigma,rdsf,sigmaf,              &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,                   &
                   rho0s,pi0s,prs0s,rth0s,                                     &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,u0,v0,                  &
                   qi0,rr0,rf0,rrf0,                                           &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,        &
                   tsk,znt,ust,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh,cm0, &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                    &
                   divx,rho,rf,prs,t11,t12,t13,t22,t23,t33,                    &
                   rru,u3d,ugr ,utmp ,                                         &
                   rrv,v3d,vgr ,vtmp ,                                         &
                   rrw,w3d,wten,wtmp ,                                         &
                   pp3d,ppten,sten,th3d,thv ,thten,thten1,                     &
                   q3d,qten,kmh,kmv,khh,khv,cme,csm,ce1,ce2,tkea,tke3d,tketen, &
                   nm,defv,defh,lenscl,dissten,                                &
                   thpten,qvpten,qcpten,qipten,upten,vpten,xkzh,xkzq,xkzm,     &
                   rain,raincv,u10,v10,s10,br,brcr,hpbl,prate,                 &
                   swten,lwten,cldfra,qke,qke_adv,el_pbl,edmf_a,               &
                   tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out2d,out3d,getdbz,getvt, &
                   effc,effi,effs,effr,effg,effis,                             &
                   lwupt,lwdnt,lwupb,lwdnb,                                    &
                   swupt,swdnt,swupb,swdnb,lwcf,swcf,coszen,                   &
                    ir, rr, angle,crref,kmwk,                                  &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,              &
                   dat1     ,dat2     ,dat3       ,                            &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,flag)
        ! end_azimavg
    use input
    use constants
    use cm1libs , only : rslf,rsif
#ifdef NETCDF
      use writeout_nc_module, only : writeazim_nc
#endif
#ifdef MPI
    use mpi
#endif
    implicit none

    integer, intent(in) :: nstep
    double precision, intent(in) :: mtime
    integer, intent(inout) :: nwritea,arecs,arecw
    character(len=3), intent(in), dimension(maxq) :: qname
    real, intent(inout) :: dt
    logical, intent(in) :: dosfcflx
    integer, intent(in) :: icrs
    integer, intent(inout) :: icenter,jcenter
    real,    intent(inout) :: xcenter,ycenter
    double precision, intent(in) :: domainlocx,domainlocy
    real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
    real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
    real, intent(in), dimension(jb:je) :: yh,vh,rvh
    real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
    real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
    real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
    real, intent(in), dimension(kb:ke) :: rds,sigma
    real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,c1,c2
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
    real, intent(in), dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0
    real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
    real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: qi0,rr0,rf0,rrf0
    real, intent(in), dimension(ib:ie,jb:je) :: zs
    real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
    real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
    real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,thflux,qvflux,cd,ch,cq,  &
                                                u1,v1,s1,xland,psfc,tlh,cm0
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rf,prs
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
    real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,ugr,utmp
    real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d,vgr,vtmp
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d,wten,wtmp
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: pp3d,ppten,sten
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th3d,thv,thten,thten1
    real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d,qten
    real, intent(in),    dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
    real, intent(in),    dimension(ibc:iec,jbc:jec,kbc:kec) :: cme,csm,ce1,ce2
    real, intent(in),    dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,lenscl,dissten
    real, intent(in), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten
    real, intent(in), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: xkzh,xkzq,xkzm
    real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain,raincv
    real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,br,brcr,hpbl
    real, intent(in), dimension(ib:ie,jb:je) :: prate
    real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,cldfra
    real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: qke,qke_adv,el_pbl,edmf_a
    real, intent(in) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
    real, intent(in) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
    real, intent(in) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag

    real, intent(inout), dimension(ib2d:ie2d,jb2d:je2d,nout2d) :: out2d
    real, intent(in) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
    logical, intent(in) :: getdbz,getvt
    real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: effc,effi,effs,effr,effg,effis
    real, intent(in), dimension(ibr:ier,jbr:jer) :: lwupt,lwdnt,lwupb,lwdnb
    real, intent(in), dimension(ibr:ier,jbr:jer) :: swupt,swdnt,swupb,swdnb
    real, intent(in), dimension(ibr:ier,jbr:jer) :: lwcf,swcf,coszen
    integer, intent(inout), dimension(ib:ie,jb:je) :: ir
    real, intent(inout), dimension(ib:ie,jb:je) :: rr,angle
    real, intent(in), dimension(icrs) :: crref
    real, intent(in), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: kmwk
    real, intent(in), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten
    real, intent(inout), dimension(d3is,d3js) :: dat1
    real, intent(inout), dimension(d2is,d2js) :: dat2
    real, intent(inout), dimension(d3is,d3js,0:d3n-1) :: dat3
    real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
    real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
    logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag

    !---------------------------------------------------------------------------

    integer :: i,j,k,n,nn,nq,nvar,nfile,ntmp,wloopmax
    integer :: ip,im
    real :: rddr
    character(len=80), dimension(varmax) :: varname,vardesc,varunit
    integer, dimension(varmax) :: varlvls
    character(len=1),  dimension(varmax) :: vargrid
    character(len=80) :: a1,a2
    character(len=16) :: a16

    integer, dimension(icrs) :: navg
    real, dimension(icrs) :: ravg
    double precision, dimension(icrs) :: area,savg2d
    double precision, dimension(icrs,nk+1) :: savg
    real, dimension(icrs,nk) :: uavg,vavg,wavg,wsavg,pavg,rhoavg,thavg,qvavg,mavg,thvavg,ppiavg
    real :: pisfc,prssfc,thsfc,qvsfc,ssfc,pd,p1,p2,p3,tx,qx,px,tlcl,ee,ql,qi,qvs,rrtmp
    real :: ff,qvsl,qvsi,qt,ucrit
    real :: xsave,ysave
    logical :: doit,dovt

    if( myid.eq.0 ) print *
    if( myid.eq.0 ) print *,' ..... begin azimavg code ..... '
    if( myid.eq.0 ) print *

    rddr = 1.0/ddr

    ! by default, assume 3d and scalar levels:
    varlvls = nk
    vargrid = 's'

    savg2d = 0.0
    savg = 0.0

      IF( imove.eq.1 )THEN
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ! get ground-relative winds:
          ugr(i,j,k) = u3d(i,j,k)+umove
          vgr(i,j,k) = v3d(i,j,k)+vmove
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ugr(i,j,k) = u3d(i,j,k)
          vgr(i,j,k) = v3d(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

    !---------------------------------------------------------------------------
    !  Get center point:

          call getcenter(nstep,mtime,rddr,xh,yh,zh,xfref,yfref,ugr,utmp,vgr,vtmp,pp3d,  &
                         nwritea,icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy)

          call   getrangle(xh,yh,xcenter,ycenter, ir, rr, angle,icrs,crref,xhref,xfref)

    !---------------------------------------------------------------------------
    ! open file:

      if(myid.eq.0) print *,'  nwritea = ',nwritea

      IF( output_format.eq.1 )THEN

        ! grads-format

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_azimavg_XXXXXX_s.dat'
        write(string(16:21),151) nwritea
        if( myid.eq.0 )then
          print *,string
          open(unit=fnums,file=string,form='unformatted',access='direct',recl=4*icrs)
        endif

        string = 'cm1out_azimavg_XXXXXX_w.dat'
        write(string(16:21),151) nwritea
        if( myid.eq.0 )then
          print *,string
          open(unit=fnumw,file=string,form='unformatted',access='direct',recl=4*icrs)
        endif

151     format(i6.6)

        arecs = 1
        arecw = 1

        wloopmax = 1

      ELSEIF( output_format.eq.2 )THEN

        ! netcdf-format
#ifdef NETCDF

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_azimavg_XXXXXX.nc'
        write(string(16:21),151) nwritea
        if( myid.eq.0 )then
          print *,string
        endif

        wloopmax = 2

#endif
      ENDIF

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      navg = 0.0
      area = 0.0

      do j=1,nj
      do i=1,ni
        if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
          navg(ir(i,j)) = navg(ir(i,j)) + 1
          area(ir(i,j)) = area(ir(i,j)) + dx*dy*ruh(i)*rvh(j)
        endif
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,navg,icrs,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,area,icrs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      ravg = 0.0

      do i=1,icrs
!!!        ravg(i) = 0.5*ddr + (i-1)*ddr
        ravg(i) = xhref( nx/2 + i )
      enddo

  !-----------------------------------------------------------------------------

  bigwloop:  &
  DO wloop = 1 , wloopmax

      if( myid.eq.0 ) print *,'  wloop = ',wloop

      nvar = 0

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  2d variables:

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'mtime'
      vardesc(nvar) = 'model time (seconds since beginning of simulation)'
      varunit(nvar) = 's'
      varlvls(nvar) = 0

      do i=1,icrs
        savg2d(i) = mtime
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'navg'
      vardesc(nvar) = 'number of data points in each bin'
      varunit(nvar) = 'count'
      varlvls(nvar) = 0

      do i=1,icrs
        savg2d(i) = navg(i)
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ravg'
      vardesc(nvar) = 'avg radius'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = rr(i,j)
      enddo
      enddo

      call     getavg2d(dum1(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ravg2'
      vardesc(nvar) = 'actual radius'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      do i=1,icrs
        savg2d(i) = ravg(i)
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'area'
      vardesc(nvar) = 'surface area in each averaging ring'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      do i=1,icrs
        savg2d(i) = area(i)
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'rain'
      vardesc(nvar) = 'accumulated surface rainfall'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call     getavg2d(rain(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. nrain.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'rain2'
      vardesc(nvar) = 'translated surface rainfall'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call     getavg2d(rain(ib,jb,2),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. cuparam.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'raincv'
      vardesc(nvar) = 'accumulated surface rainfall (convective, param.)'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call     getavg2d(raincv(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. nrain.eq.2 .and. cuparam.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'raincv2'
      vardesc(nvar) = 'translated surface rainfall (convective, param.)'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call     getavg2d(raincv(ib,jb,2),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'prate'
      vardesc(nvar) = 'surface precipitation rate'
      varunit(nvar) = 'kg/m2/s'
      varlvls(nvar) = 0

      call     getavg2d(prate,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. nqv.ge.1 )THEN

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
            do i=1,ni
              dum1(i,j,1) = 0.0
            enddo
              do k=1,nk
              do i=1,ni
                                                                                ! 1000 kg/m3
                dum1(i,j,1) = dum1(i,j,1) + rho(i,j,k)*q3d(i,j,k,nqv)*dz*rmh(i,j,k)/1000.0
              enddo
              enddo
          enddo

      nvar = nvar+1
      varname(nvar) = 'pwat'
      vardesc(nvar) = 'precipitable water'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'cd'
      vardesc(nvar) = 'surface drag coefficient'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(cd,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ch'
      vardesc(nvar) = 'surface sensible heat exchange coeff'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(ch,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'cq'
      vardesc(nvar) = 'surface latent heat exchange coeff'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(cq,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'znt'
      vardesc(nvar) = 'surface roughness length'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call     getavg2d(znt,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ust'
      vardesc(nvar) = 'surface friction velocity'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(ust,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        do j=1,nj
        do i=1,ni
          if( abs(rr(i,j)).lt.1.0e-4 )then
            dum3(i,j,1)=0.0
            dum4(i,j,1)=0.0
          else
            dum3(i,j,1)=v10(i,j)*sin(angle(i,j))+u10(i,j)*cos(angle(i,j))
            dum4(i,j,1)=v10(i,j)*cos(angle(i,j))-u10(i,j)*sin(angle(i,j))
          endif
        enddo
        enddo

      nvar = nvar+1
      varname(nvar) = 'u10'
    if( imove.eq.1 )then
      vardesc(nvar) = 'diagnostic 10m radial velocity (ground-rel.)'
    else
      vardesc(nvar) = 'diagnostic 10m radial velocity'
    endif
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(dum3(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'v10'
    if( imove.eq.1 )then
      vardesc(nvar) = 'diagnostic 10m tangential velocity (ground-rel.)'
    else
      vardesc(nvar) = 'diagnostic 10m tangential velocity'
    endif
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(dum4(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 's10'
    if( imove.eq.1 )then
      vardesc(nvar) = 'diagnostic 10m horiz wind speed (ground-rel.)'
    else
      vardesc(nvar) = 'diagnostic 10m horiz wind speed'
    endif
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(s10,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 's10max'
    if( imove.eq.1 )then
      vardesc(nvar) = 'max diagnostic 10m horiz wind speed (ground-rel.)'
    else
      vardesc(nvar) = 'max diagnostic 10m horiz wind speed'
    endif
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getmax2d(s10,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

    IF( ipbl.eq.1 .or. ipbl.eq.3 )THEN

      nvar = nvar+1
      varname(nvar) = 'br'
      vardesc(nvar) = 'bulk Richardson number in surface layer'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(br,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( ipbl.ge.1 .or. sfcmodel.ge.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'hpbl'
    if( ipbl.eq.1 )then
      vardesc(nvar) = 'PBL height (from YSU scheme)'
    elseif( ipbl.eq.3 )then
      vardesc(nvar) = 'PBL height (from GFSEDMF scheme)'
    elseif( ipbl.eq.4 .or. ipbl.eq.5 )then
      vardesc(nvar) = 'PBL height (from MYNN scheme)'
    else
      vardesc(nvar) = 'estimated PBL height (based on bulk Ri)'
    endif
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call     getavg2d(hpbl,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'thflux'
      vardesc(nvar) = 'surface potential temperature flux'
      varunit(nvar) = 'K m/s'
      varlvls(nvar) = 0

      call     getavg2d(thflux,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'qvflux'
      vardesc(nvar) = 'surface water vapor flux'
      varunit(nvar) = 'g/g m/s'
      varlvls(nvar) = 0

      call     getavg2d(qvflux,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      do j=1,nj
      do i=1,ni
        p1 = pi0(1,1,1)+pp3d(i,j,1)
        p2 = pi0(1,1,2)+pp3d(i,j,2)
        p3 = pi0(1,1,3)+pp3d(i,j,3)
        pisfc=cgs1*p1+cgs2*p2+cgs3*p3
        prssfc = psfc(i,j)
        thsfc = tsk(i,j)/pisfc
        qvsfc = rslf(prssfc,tsk(i,j))
        pd = prssfc/(1.0+qvsfc*reps)
        ssfc = cp*alog(tsk(i,j))   &
             + 2555000.0*qvsfc/tsk(i,j)   &
             - rd*alog(pd)
        dum1(i,j,1) = pisfc
        dum1(i,j,2) = prssfc
        dum1(i,j,3) = thsfc
        dum1(i,j,4) = qvsfc
        dum1(i,j,5) = ssfc
      enddo
      enddo

      nvar = nvar+1
      varname(nvar) = 'pisfc'
      vardesc(nvar) = 'surface nondimensional pressure'
      varunit(nvar) = 'nondimensional'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,1),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'psfc'
      vardesc(nvar) = 'surface pressure'
      varunit(nvar) = 'Pa'
      varlvls(nvar) = 0

      call     getavg2d(psfc,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'tsk'
      vardesc(nvar) = 'surface (ocean/land) temperature'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(tsk,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'thsfc'
      vardesc(nvar) = 'surface (ocean/land) potential temperature'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,3),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'qvsfc'
      vardesc(nvar) = 'surface (ocean/land) water vapor mixing ratio'
      varunit(nvar) = 'g/g'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,4),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'ssfc'
      vardesc(nvar) = 'surface (ocean/land) moist entropy'
      varunit(nvar) = 'J/kg/K'
      varlvls(nvar) = 0

      call     getavg2d(dum1(ib,jb,5),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

    if( idoles .and. cm1setup.eq.4 )then
      nvar = nvar+1
      varname(nvar) = 'cm0'
      vardesc(nvar) = 'cm0'
      varunit(nvar) = 'unitless'
      varlvls(nvar) = 0

      call     getavg2d(cm0,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    endif

      !c-c-c-c-c-c-c-c-c-c

    IF( (cm1setup.eq.2.or.cm1setup.eq.4) .and. horizturb.eq.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'tlh'
      vardesc(nvar) = 'horiz lengthscale for turbulence scheme'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call     getavg2d(tlh,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    if( radopt.ge.2 )then

      nvar = nvar+1
      varname(nvar) = 'lwupt'
      vardesc(nvar) = 'lw flux, upward, top of atmosphere (OLR)'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call     getavg2d(lwupt(ibr,jbr),icrs,savg2d,ir,rr,navg,area,ruh,rvh)
      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

    endif

      !c-c-c-c-c-c-c-c-c-c

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! 3d vars:  averages

              ! dum1 = u at s pts
              ! dum2 = v at s pts
              ! dum3 = ur at s pts
              ! dum4 = vt at s pts
              ! dum5 = w at s pts
              ! dum6 = wsp at s pts

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.5*(ugr(i,j,k)+ugr(i+1,j,k))
          dum2(i,j,k) = 0.5*(vgr(i,j,k)+vgr(i,j+1,k))
          if( abs(rr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(angle(i,j))+dum1(i,j,k)*cos(angle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(angle(i,j))-dum1(i,j,k)*sin(angle(i,j))
          endif
          dum5(i,j,k) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
          dum6(i,j,k) = sqrt( dum1(i,j,k)**2 + dum2(i,j,k)**2 )
        enddo
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'u'
        vardesc(nvar) = 'radial velocity'
        varunit(nvar) = 'm/s'

        call getavg3d(dum3,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save uavg:
        do k=1,nk
        do i=1,icrs
          uavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'v'
        vardesc(nvar) = 'tangential velocity'
        varunit(nvar) = 'm/s'

        call getavg3d(dum4,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save vavg:
        do k=1,nk
        do i=1,icrs
          vavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'umax'
        vardesc(nvar) = 'maximum radial wind speed'
        varunit(nvar) = 'm/s'

        call getmax3d(dum3,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'umin'
        vardesc(nvar) = 'minimum radial wind speed'
        varunit(nvar) = 'm/s'

        call getmin3d(dum3,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vmax'
        vardesc(nvar) = 'maximum tangential wind speed'
        varunit(nvar) = 'm/s'

        call getmax3d(dum4,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vmin'
        vardesc(nvar) = 'minimum tangential wind speed'
        varunit(nvar) = 'm/s'

        call getmin3d(dum4,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'zinfl'
      vardesc(nvar) = 'depth of inflow layer (u < -2 m/s)'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      ucrit = -2.0

      savg2d = 0.0

      do i=1,icrs
        if( uavg(i,1).lt.ucrit )then
          k = 1
          do while( uavg(i,k).lt.ucrit .and. k.lt.nk )
            k = k+1
          enddo
          savg2d(i) = zh(1,1,k-1)+(zh(1,1,k)-zh(1,1,k-1))  &
                                  *(    ucrit-uavg(i,k-1))  &
                                  /(uavg(i,k)-uavg(i,k-1))
        else
          savg2d(i) = 0.0
        endif
      enddo

      call write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'frac_in'
        vardesc(nvar) = 'fraction of area with inflow'
        varunit(nvar) = 'unitless'

        call getfrac3d(dum3,icrs,savg2d,savg,ir,rr,navg,area,ruh,rvh,2)

        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'frac_out'
        vardesc(nvar) = 'fraction of area with outflow'
        varunit(nvar) = 'unitless'

        call getfrac3d(dum3,icrs,savg2d,savg,ir,rr,navg,area,ruh,rvh,1)

        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'm'
        vardesc(nvar) = 'angular momentum'
        varunit(nvar) = 'm^2/s'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = rr(i,j)*dum4(i,j,k) + 0.5*fcor*rr(i,j)*rr(i,j)
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save mavg:
        do k=1,nk
        do i=1,icrs
          mavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'wsp'
        vardesc(nvar) = 'horizontal wind speed'
        varunit(nvar) = 'm/s'

        call getavg3d(dum6,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'wspmax'
        vardesc(nvar) = 'maximum horizontal wind speed'
        varunit(nvar) = 'm/s'

        call getmax3d(dum6,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'wspmin'
        vardesc(nvar) = 'minimum horizontal wind speed'
        varunit(nvar) = 'm/s'

        call getmin3d(dum6,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! w interpolated to scalar levels:
        nvar = nvar+1
        varname(nvar) = 'winterp'
        vardesc(nvar) = 'vertical velocity'
        varunit(nvar) = 'm/s'

        call getavg3d(dum5,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save wavg:
        do k=1,nk
        do i=1,icrs
          wavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

      IF( idoles .and. iusetke )THEN

        ! tke interpolated to s levels:
        nvar = nvar+1
        varname(nvar) = 'tkeinterp'
        vardesc(nvar) = 'subgrid tke (interp to s levels)'
        varunit(nvar) = 'm2/s2'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = 0.5*(tke3d(i,j,k)+tke3d(i,j,k+1))
        enddo
        enddo
        enddo

        call getavg3d(dum7(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

        ! w at staggered level:
        nvar = nvar+1
        varname(nvar) = 'w'
        vardesc(nvar) = 'vertical velocity'
        varunit(nvar) = 'm/s'
        vargrid(nvar) = 'w'   ! w levels

        call getavg3d(w3d(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save wsavg:
        do k=1,nk
        do i=1,icrs
          wsavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        ! w at staggered level:
        nvar = nvar+1
        varname(nvar) = 'tmfu'
        vardesc(nvar) = 'total upward mass flux'
        varunit(nvar) = 'kg/s'
        vargrid(nvar) = 's'   ! w levels

      if( imoist.eq.0 )then
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = rho(i,j,k)*max(0.0,0.5*(w3d(i,j,k)+w3d(i,j,k+1)))*ruh(i)*rvh(j)*(dx*dy)
        enddo
        enddo
        enddo
      else
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = rho(i,j,k)*(1.0+q3d(i,j,k,nqv))*max(0.0,0.5*(w3d(i,j,k)+w3d(i,j,k+1)))*ruh(i)*rvh(j)*(dx*dy)
        enddo
        enddo
        enddo
      endif

        call getsum3d(dum1(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! w at staggered level:
        nvar = nvar+1
        varname(nvar) = 'tmfd'
        vardesc(nvar) = 'total downward mass flux'
        varunit(nvar) = 'kg/s'
        vargrid(nvar) = 's'   ! w levels

      if( imoist.eq.0 )then
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = rho(i,j,k)*min(0.0,0.5*(w3d(i,j,k)+w3d(i,j,k+1)))*ruh(i)*rvh(j)*(dx*dy)
        enddo
        enddo
        enddo
      else
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = rho(i,j,k)*(1.0+q3d(i,j,k,nqv))*min(0.0,0.5*(w3d(i,j,k)+w3d(i,j,k+1)))*ruh(i)*rvh(j)*(dx*dy)
        enddo
        enddo
        enddo
      endif

        call getsum3d(dum1(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! w at staggered level:
        nvar = nvar+1
        varname(nvar) = 'wmax'
        vardesc(nvar) = 'maximum vertical velocity'
        varunit(nvar) = 'm/s'
        vargrid(nvar) = 'w'   ! w levels

        call getmax3d(w3d(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! w at staggered level:
        nvar = nvar+1
        varname(nvar) = 'wmin'
        vardesc(nvar) = 'minimum vertical velocity'
        varunit(nvar) = 'm/s'
        vargrid(nvar) = 'w'   ! w levels

        call getmin3d(w3d(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! fraction of ascending air
        nvar = nvar+1
        varname(nvar) = 'fracup'
        vardesc(nvar) = 'fraction of area that is ascending'
        varunit(nvar) = 'unitless'
        vargrid(nvar) = 'w'   ! w levels

        call getfrac3d(w3d(ib,jb,kb),icrs,savg2d,savg,ir,rr,navg,area,ruh,rvh,1)

        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        ! fraction of descending air
        nvar = nvar+1
        varname(nvar) = 'fracdn'
        vardesc(nvar) = 'fraction of area that is descending'
        varunit(nvar) = 'unitless'
        vargrid(nvar) = 'w'   ! w levels

        call getfrac3d(w3d(ib,jb,kb),icrs,savg2d,savg,ir,rr,navg,area,ruh,rvh,2)

        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'rho'
        vardesc(nvar) = 'dry air density'
        varunit(nvar) = 'kg/m3'

        call getavg3d(rho,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save rhoavg:
        do k=1,nk
        do i=1,icrs
          rhoavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'prs'
        vardesc(nvar) = 'pressure'
        varunit(nvar) = 'Pa'

        call getavg3d(prs,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save pavg:
        do k=1,nk
        do i=1,icrs
          pavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'ppi'
        vardesc(nvar) = 'nondimensional pressure perturbation'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pp3d,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save ppiavg:
        do k=1,nk
        do i=1,icrs
          ppiavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'th'
        vardesc(nvar) = 'potential temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = th0(i,j,k)+th3d(i,j,k)
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save thavg:
        do k=1,nk
        do i=1,icrs
          thavg(i,k) = savg(i,k)
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 't'
        vardesc(nvar) = 'temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'rh'
        vardesc(nvar) = 'relative humidity wrt liquid'
        varunit(nvar) = 'nondimensional'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = q3d(i,j,k,nqv)/rslf(prs(i,j,k),dum1(i,j,k))
        enddo
        enddo
        enddo

        call getavg3d(dum2,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'rh95frac'
        vardesc(nvar) = 'fraction of points that have rh >= 95'
        varunit(nvar) = 'nondimensional'

        dum7 = 0.0
        dum8 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( dum2(i,j,k).ge.0.95 ) dum7(i,j,k) = 1.0
            if( dum2(i,j,k).ge.0.99 ) dum8(i,j,k) = 1.0
          enddo
          enddo
          enddo

        call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'rh99frac'
        vardesc(nvar) = 'fraction of points that have rh >= 99'
        varunit(nvar) = 'nondimensional'

        call getavg3d(dum8,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. iice.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'rhi'
        vardesc(nvar) = 'relative humidity wrt ice'
        varunit(nvar) = 'nondimensional'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = q3d(i,j,k,nqv)/rsif(prs(i,j,k),dum1(i,j,k))
        enddo
        enddo
        enddo

        call getavg3d(dum2,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qv'
        vardesc(nvar) = 'water vapor mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqv),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save qvavg:
        do k=1,nk
        do i=1,icrs
          qvavg(i,k) = savg(i,k)
        enddo
        enddo

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqc.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qc'
        vardesc(nvar) = 'cloud water mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqc),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqr.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qr'
        vardesc(nvar) = 'rain water mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqi.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qi'
        vardesc(nvar) = 'cloud ice mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqi),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqs.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qs'
        vardesc(nvar) = 'snow mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqs),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. nqg.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'qg'
        vardesc(nvar) = 'graupel/hail mixing ratio'
        varunit(nvar) = 'g/g'

        call getavg3d(q3d(ibm,jbm,kbm,nqg),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. output_dbz.eq.1 .and. qd_dbz.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'dbz'
        vardesc(nvar) = 'reflectivity'
        varunit(nvar) = 'dBZ'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_dbz),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      doit = .true.

      dorturb:  &
      if( doit )then
        IF( cm1setup.eq.4 )THEN

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'upup'
          vardesc(nvar) = "<u'u'>, resolved rad vel var"
          varunit(nvar) = 'm2/s2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (dum3(i,j,k)-uavg(ir(i,j),k))**2
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'vpvp'
          vardesc(nvar) = "<v'v'>, resolved tang vel var"
          varunit(nvar) = 'm2/s2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (dum4(i,j,k)-vavg(ir(i,j),k))**2
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'wpwp'
          vardesc(nvar) = "<w'w'>, resolved tang vel var"
          varunit(nvar) = 'm2/s2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (dum5(i,j,k)-wavg(ir(i,j),k))**2
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'upvp'
          vardesc(nvar) = "<u'v'>"
          varunit(nvar) = 'm2/s2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (dum3(i,j,k)-uavg(ir(i,j),k))*(dum4(i,j,k)-vavg(ir(i,j),k))
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'upwp'
          vardesc(nvar) = "<u'w'>"
          varunit(nvar) = 'm2/s2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (dum3(i,j,k)-uavg(ir(i,j),k))*(dum5(i,j,k)-wavg(ir(i,j),k))
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'vpwp'
          vardesc(nvar) = "<v'w'>"
          varunit(nvar) = 'm2/s2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (dum4(i,j,k)-vavg(ir(i,j),k))*(dum5(i,j,k)-wavg(ir(i,j),k))
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'wptp'
          vardesc(nvar) = "<w't'>"
          varunit(nvar) = 'm/s K'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = ((th0(i,j,k)-thavg(ir(i,j),k))+th3d(i,j,k))*(dum5(i,j,k)-wavg(ir(i,j),k))
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        IF(imoist.eq.1)THEN

          nvar = nvar+1
          varname(nvar) = 'wpqp'
          vardesc(nvar) = "<w'q'>"
          varunit(nvar) = 'm/s kg/kg'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (q3d(i,j,k,nqv)-qvavg(ir(i,j),k))*(dum5(i,j,k)-wavg(ir(i,j),k))
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ENDIF

      !c-c-c-c-c-c-c-c-c-c

          nvar = nvar+1
          varname(nvar) = 'tptp'
          vardesc(nvar) = "<t't'>"
          varunit(nvar) = 'K^2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = ((th0(i,j,k)-thavg(ir(i,j),k))+th3d(i,j,k))**2
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        IF(imoist.eq.1)THEN

          nvar = nvar+1
          varname(nvar) = 'qpqp'
          vardesc(nvar) = "<q'q'>"
          varunit(nvar) = 'kg^2/kg^2'

          dum7 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
              dum7(i,j,k) = (q3d(i,j,k,nqv)-qvavg(ir(i,j),k))**2
            endif
          enddo
          enddo
          enddo

          call getavg3d(dum7,icrs,savg,ir,rr,navg,area,ruh,rvh)
          call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ENDIF

      !c-c-c-c-c-c-c-c-c-c

        ENDIF
      endif  dorturb

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !c-c-c-c-c-c-c-c-c-c

      dovt = .true.

      IF( imoist.eq.1 .and. dovt .and. qd_vtc.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'vtc'
        vardesc(nvar) = 'vtc'
        varunit(nvar) = 'm/s'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_vtc),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. dovt .and. qd_vtr.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'vtr'
        vardesc(nvar) = 'vtr'
        varunit(nvar) = 'm/s'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_vtr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. dovt .and. qd_vti.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'vti'
        vardesc(nvar) = 'vti'
        varunit(nvar) = 'm/s'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_vti),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. dovt .and. qd_vts.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'vts'
        vardesc(nvar) = 'vts'
        varunit(nvar) = 'm/s'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_vts),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. dovt .and. qd_vtg.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'vtg'
        vardesc(nvar) = 'vtg'
        varunit(nvar) = 'm/s'

        call getavg3d(qdiag(ibdq,jbdq,kbdq,qd_vtg),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. (ptype.eq.3.or.ptype.eq.5) .and. radopt.eq.2 )THEN

        nvar = nvar+1
        varname(nvar) = 'effc'
        vardesc(nvar) = 'effc'
        varunit(nvar) = 'm'

        call getavg3d(effc(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. (ptype.eq.3.or.ptype.eq.5) .and. radopt.eq.2 )THEN

        nvar = nvar+1
        varname(nvar) = 'effi'
        vardesc(nvar) = 'effi'
        varunit(nvar) = 'm'

        call getavg3d(effi(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. (ptype.eq.3.or.ptype.eq.5) .and. radopt.eq.2 )THEN

        nvar = nvar+1
        varname(nvar) = 'effr'
        vardesc(nvar) = 'effr'
        varunit(nvar) = 'm'

        call getavg3d(effr(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. (ptype.eq.3.or.ptype.eq.5) .and. radopt.eq.2 )THEN

        nvar = nvar+1
        varname(nvar) = 'effg'
        vardesc(nvar) = 'effg'
        varunit(nvar) = 'm'

        call getavg3d(effg(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. (ptype.eq.3.or.ptype.eq.5) .and. radopt.eq.2 )THEN

        nvar = nvar+1
        varname(nvar) = 'effs'
        vardesc(nvar) = 'effs'
        varunit(nvar) = 'm'

        call getavg3d(effs(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. (ptype.eq.3.or.ptype.eq.5) .and. radopt.eq.2 )THEN

        nvar = nvar+1
        varname(nvar) = 'effis'
        vardesc(nvar) = 'effis'
        varunit(nvar) = 'm'

        call getavg3d(effis(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c


      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'thv'
        vardesc(nvar) = 'virtual potential temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          ql = 0.0
          do nq=nql1,nql2
            ql=ql+q3d(i,j,k,nq)
          enddo
          qi = 0.0
          if(iice.eq.1)then
          do nq=nqs1,nqs2
            qi=qi+q3d(i,j,k,nq)
          enddo
          endif
          dum2(i,j,k) = (th0(i,j,k)+th3d(i,j,k))*(1.0+reps*q3d(i,j,k,nqv)) &
                                                /(1.0+q3d(i,j,k,nqv)+ql+qi)
        enddo
        enddo
        enddo

        call getavg3d(dum2,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        ! save thvavg:
        do k=1,nk
        do i=1,icrs
          thvavg(i,k) = savg(i,k)
        enddo
        enddo

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'the'
        vardesc(nvar) = 'equivalent potential temperature'
        varunit(nvar) = 'K'

        do k=1,nk
        do j=1,nj
        do i=1,ni
          tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
          qx = q3d(i,j,k,nqv)
          px = p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          if(q3d(i,j,k,nqc).ge.clwsat)then
            tlcl=tx
          else
            ee=0.01*px*qx/(eps+qx)
            tlcl=55.0+2840.0/(3.5*alog(tx)-log(1.0e-20+ee)-4.805)
          endif
          dum1(i,j,k) = tx*((p00/px)**(0.2854*(1.0-0.28*qx)))   &
                          *exp(((3376.0/tlcl)-2.54)*qx*(1.0+0.81*qx))
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'satfrac'
        vardesc(nvar) = 'fraction of points that are saturated'
        varunit(nvar) = 'nondimensional'

        dum1 = 0.0

        IF( iice.eq.0 )THEN
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
            qvs = rslf(prs(i,j,k),tx)
            ql = max( q3d(i,j,k,nqv)+q3d(i,j,k,nqc) - qvs , 0.0 )
            if( ql.gt.1.0e-10 ) dum1(i,j,k) = 1.0
          enddo
          enddo
          enddo
        ELSE
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
            qt = q3d(i,j,k,nqv)+q3d(i,j,k,nqc)+q3d(i,j,k,nqi)
            if( tx.ge.273.15 )then
              qvs = rslf(prs(i,j,k),tx)
              ql = max( qt-qvs , 0.0 )
              if( ql.gt.1.0e-10 ) dum1(i,j,k) = 1.0
            elseif( tx.le.233.15 )then
              qvs = rsif(prs(i,j,k),tx)
              qi = max( qt-qvs , 0.0 )
              if( qi.gt.1.0e-10 ) dum1(i,j,k) = 1.0
            else
              ff = (tx-233.15)/(273.15-233.15)
              qvsl = rslf(prs(i,j,k),tx)
              qvsi = rsif(prs(i,j,k),tx)
              qvs = ff*qvsl + (1.0-ff)*qvsi
              qi = max( 0.0 , (1.0-ff)*(qt-qvs) )
              ql = max( 0.0 , ff*(qt-qvs) )
              if( (ql+qi).gt.1.0e-10 ) dum1(i,j,k) = 1.0
            endif
          enddo
          enddo
          enddo
        ENDIF

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 .and. qd_dbz.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'dbz20frac'
        vardesc(nvar) = 'fraction of points that have dbz >= 20'
        varunit(nvar) = 'nondimensional'

        dum1 = 0.0
        dum2 = 0.0

          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( qdiag(i,j,k,qd_dbz).ge.20.0 ) dum1(i,j,k) = 1.0
            if( qdiag(i,j,k,qd_dbz).ge.40.0 ) dum2(i,j,k) = 1.0
          enddo
          enddo
          enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( imoist.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'dbz40frac'
        vardesc(nvar) = 'fraction of points that have dbz >= 40'
        varunit(nvar) = 'nondimensional'

        call getavg3d(dum2,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      radout:  &
      IF( radopt.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'swten'
        vardesc(nvar) = 'temperature tendency, sw radiation'
        varunit(nvar) = 'K/s'

        call getavg3d(swten(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !------------

        nvar = nvar+1
        varname(nvar) = 'lwten'
        vardesc(nvar) = 'temperature tendency, lw radiation'
        varunit(nvar) = 'K/s'

        call getavg3d(lwten(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !------------

        nvar = nvar+1
        varname(nvar) = 'cldfra'
        vardesc(nvar) = 'cloud fraction from radiation scheme'
        varunit(nvar) = 'nondimensional'

        call getavg3d(cldfra(ibr,jbr,kbr),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF  radout

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vgrad'
        vardesc(nvar) = 'gradient wind speed (from azim-avg fields)'
        varunit(nvar) = 'm/s'

        do k=1,nk
        do i=1,icrs
          ip = min( i+1 , icrs )
          im = max( i-1 , 1 )
          
          savg(i,k) = -0.5*fcor*ravg(i) + sqrt( max(0.0,               &
                            0.25*fcor*fcor*ravg(i)*ravg(i)               &
           +ravg(i)*cp*thvavg(i,k)*(ppiavg(ip,k)-ppiavg(im,k))/(ravg(ip)-ravg(im))   &
                                           ) )
        enddo
        enddo

        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

!      IF( idoles .and. iusetke )THEN
!
!        nvar = nvar+1
!        varname(nvar) = 'cme'
!        vardesc(nvar) = 'parameter in subgrid TKE scheme'
!        varunit(nvar) = 'unitless'
!        vargrid(nvar) = 'w'
!
!        call getavg3d(cme(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
!        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
!
!      ENDIF

      !c-c-c-c-c-c-c-c-c-c

!      IF( idoles .and. iusetke )THEN
!
!        nvar = nvar+1
!        varname(nvar) = 'csm'
!        vardesc(nvar) = 'parameter in subgrid TKE scheme'
!        varunit(nvar) = 'unitless'
!        vargrid(nvar) = 'w'
!
!        call getavg3d(csm(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
!        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
!
!      ENDIF

      !c-c-c-c-c-c-c-c-c-c

!      IF( idoles .and. iusetke )THEN
!
!        nvar = nvar+1
!        varname(nvar) = 'ce1'
!        vardesc(nvar) = 'parameter in subgrid TKE scheme'
!        varunit(nvar) = 'unitless'
!        vargrid(nvar) = 'w'
!
!        call getavg3d(ce1(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
!        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
!
!      ENDIF

      !c-c-c-c-c-c-c-c-c-c

!      IF( idoles .and. iusetke )THEN
!
!        nvar = nvar+1
!        varname(nvar) = 'ce2'
!        vardesc(nvar) = 'parameter in subgrid TKE scheme'
!        varunit(nvar) = 'unitless'
!        vargrid(nvar) = 'w'
!
!        call getavg3d(ce2(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
!        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
!
!      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( idoles .and. iusetke )THEN

        nvar = nvar+1
        varname(nvar) = 'tke'
        vardesc(nvar) = 'subgrid turbulence kinetic energy'
        varunit(nvar) = 'm^2/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(tke3d(ibt,jbt,kbt),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekm.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'kmh'
        vardesc(nvar) = 'eddy viscosity, horiz direction'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(kmh(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekm.eq.1 )THEN
      if( sgsmodel.ge.1 .or. ipbl.eq.2 )then

        nvar = nvar+1
        varname(nvar) = 'kmv'
        vardesc(nvar) = 'eddy viscosity, vert direction'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(kmv(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif
      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekh.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'khh'
        vardesc(nvar) = 'eddy diffusivity, horiz direction'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(khh(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( iusekh.eq.1 )THEN
      if( sgsmodel.ge.1 .or. ipbl.eq.2 )then

        nvar = nvar+1
        varname(nvar) = 'khv'
        vardesc(nvar) = 'eddy diffusivity, vert direction'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(khv(ibc,jbc,kbc),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif
      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      if( sgsmodel.eq.4 )then

        nvar = nvar+1
        varname(nvar) = 'kmw'
        vardesc(nvar) = 'kmw'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        dum1 = 0.0

        do k=1,ntwk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = kmwk(i,j,k)
        enddo
        enddo
        enddo

        call getavg3d(dum1,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif

      !c-c-c-c-c-c-c-c-c-c

      IF( ipbl.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'xkzh'
        vardesc(nvar) = 'eddy diffusivity for heat (from YSU)'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(xkzh(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

!        nvar = nvar+1
!        varname(nvar) = 'xkzq'
!        vardesc(nvar) = 'eddy diffusivity for moisture (from YSU)'
!        varunit(nvar) = 'm^2/s'
!        vargrid(nvar) = 'w'
!
!        call getavg3d(xkzq(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
!        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'xkzm'
        vardesc(nvar) = 'eddy viscosity (from YSU)'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(xkzm(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( ipbl.eq.3 )THEN

        nvar = nvar+1
        varname(nvar) = 'dkt3d'
        vardesc(nvar) = 'Thermal Diffusivity (from GFSEDMF)'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(xkzh(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'dku3d'
        vardesc(nvar) = 'Momentum Diffusivity (from GFSEDMF)'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(xkzm(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN

        nvar = nvar+1
        varname(nvar) = 'qke'
        vardesc(nvar) = 'twice TKE from MYNN'
        varunit(nvar) = 'm2/s2'
        vargrid(nvar) = 's'

        call getavg3d(qke(ibmynn,jbmynn,kbmynn),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'qke_adv'
        vardesc(nvar) = 'advection of twice TKE from MYNN'
        varunit(nvar) = 'm2/s2'
        vargrid(nvar) = 's'

        call getavg3d(qke_adv(ibmynn,jbmynn,kbmynn),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'exch_h'
        vardesc(nvar) = 'Thermal Diffusivity (from MYNN)'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(xkzh(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'exch_m'
        vardesc(nvar) = 'Momentum Diffusivity (from MYNN)'
        varunit(nvar) = 'm^2/s'
        vargrid(nvar) = 'w'

        call getavg3d(xkzm(ibb,jbb,kbb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'el_pbl'
        vardesc(nvar) = 'length scale in MYNN'
        varunit(nvar) = 'm'
        vargrid(nvar) = 'w'

        call getavg3d(el_pbl(ibmynn,jbmynn,kbmynn),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'edmf_a'
        vardesc(nvar) = 'edmf_a'
        varunit(nvar) = 'n.d.'
        vargrid(nvar) = 's'

        call getavg3d(edmf_a(ibmynn,jbmynn,kbmynn),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

        !--------!

      if( sgsmodel.ge.1 )then
        nvar = nvar+1
        varname(nvar) = 'lenscl'
        vardesc(nvar) = 'length scale in subgrid turbulence model'
        varunit(nvar) = 'm'
        vargrid(nvar) = 'w'

        call getavg3d(lenscl(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
      endif

        !--------!

        nvar = nvar+1
        varname(nvar) = 'dissten'
        vardesc(nvar) = 'dissipation rate'
        varunit(nvar) = 'm2/s3'
        vargrid(nvar) = 'w'

        call getavg3d(dissten(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

      IF( sgsmodel.ge.1 .or. output_nm.eq.1 .or. ipbl.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'nm'
        vardesc(nvar) = 'squared Brunt-Vaisala frequency'
        varunit(nvar) = '1/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(nm(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      sgshdiag:  &
      IF( sgsmodel.ge.1 .and. ud_hturb.ge.1 .and. vd_hturb.ge.1 )THEN

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = udiag(i,j,k,ud_hturb)
          dum2(i,j,k) = vdiag(i,j,k,vd_hturb)
          if( abs(rr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(angle(i,j))+dum1(i,j,k)*cos(angle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(angle(i,j))-dum1(i,j,k)*sin(angle(i,j))
          endif
        enddo
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'urhturb'
        vardesc(nvar) = 'sgsmodel tendency: radial velocity'
        varunit(nvar) = 'm/s/s'

        call getavg3d(dum3,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vthturb'
        vardesc(nvar) = 'sgsmodel tendency: tangential velocity'
        varunit(nvar) = 'm/s/s'

        call getavg3d(dum4,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF  sgshdiag

      !c-c-c-c-c-c-c-c-c-c

      sgsvdiag:  &
      IF( sgsmodel.ge.1 .and. ud_vturb.ge.1 .and. vd_vturb.ge.1 )THEN

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = udiag(i,j,k,ud_vturb)
          dum2(i,j,k) = vdiag(i,j,k,vd_vturb)
          if( abs(rr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(angle(i,j))+dum1(i,j,k)*cos(angle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(angle(i,j))-dum1(i,j,k)*sin(angle(i,j))
          endif
        enddo
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'urvturb'
        vardesc(nvar) = 'sgsmodel tendency: radial velocity'
        varunit(nvar) = 'm/s/s'

        call getavg3d(dum3,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vtvturb'
        vardesc(nvar) = 'sgsmodel tendency: tangential velocity'
        varunit(nvar) = 'm/s/s'

        call getavg3d(dum4,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF  sgsvdiag

      !c-c-c-c-c-c-c-c-c-c

      pblt:  &
      IF( use_pbl )THEN

        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = upten(i,j,k)
          dum2(i,j,k) = vpten(i,j,k)
          if( abs(rr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(angle(i,j))+dum1(i,j,k)*cos(angle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(angle(i,j))-dum1(i,j,k)*sin(angle(i,j))
          endif
        enddo
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'urpten'
        vardesc(nvar) = 'PBL tendency: radial velocity'
        varunit(nvar) = 'm/s/s'

        call getavg3d(dum3,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      !c-c-c-c-c-c-c-c-c-c

        nvar = nvar+1
        varname(nvar) = 'vtpten'
        vardesc(nvar) = 'PBL tendency: tangential velocity'
        varunit(nvar) = 'm/s/s'

        call getavg3d(dum4,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF  pblt

      !c-c-c-c-c-c-c-c-c-c

        !--------!

      IF( cm1setup.ge.1 .or. output_def.eq.1 .or. ipbl.ge.1 .or. horizturb.eq.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'defv'
        vardesc(nvar) = 'vertical deformation'
        varunit(nvar) = '1/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(defv(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'defh'
        vardesc(nvar) = 'horizontal deformation'
        varunit(nvar) = '1/s^2'
        vargrid(nvar) = 'w'

        call getavg3d(defh(ib,jb,kb),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      ENDIF

        !--------!

      !c-c-c-c-c-c-c-c-c-c

      IF( cuparam.ge.1 )THEN

        !--------!

        nvar = nvar+1
        varname(nvar) = 'thcuten'
        vardesc(nvar) = 'thcuten'
        varunit(nvar) = 'K/s'
        vargrid(nvar) = 's'

        call getavg3d(thcuten(ibcu,jbcu,kbcu),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'qvcuten'
        vardesc(nvar) = 'qvcuten'
        varunit(nvar) = 'g/g/s'
        vargrid(nvar) = 's'

        call getavg3d(qvcuten(ibcu,jbcu,kbcu),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'qccuten'
        vardesc(nvar) = 'qccuten'
        varunit(nvar) = 'g/g/s'
        vargrid(nvar) = 's'

        call getavg3d(qccuten(ibcu,jbcu,kbcu),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'qicuten'
        vardesc(nvar) = 'qicuten'
        varunit(nvar) = 'g/g/s'
        vargrid(nvar) = 's'

        call getavg3d(qicuten(ibcu,jbcu,kbcu),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'ucuten'
        vardesc(nvar) = 'ucuten'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 's'

        call getavg3d(ucuten(ibcu,jbcu,kbcu),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

        nvar = nvar+1
        varname(nvar) = 'vcuten'
        vardesc(nvar) = 'vcuten'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 's'

        call getavg3d(vcuten(ibcu,jbcu,kbcu),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !--------!

      ENDIF

      !c-c-c-c-c-c-c-c-c-c


      IF( pdcomp )THEN

        nvar = nvar+1
        varname(nvar) = 'pipb'
        vardesc(nvar) = 'diagnosed pi-prime: buoyancy component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,1),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pipdl'
        vardesc(nvar) = 'diagnosed pi-prime: linear dynamic component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,2),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pipdn'
        vardesc(nvar) = 'diagnosed pi-prime: nonlinear dynamic component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,3),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

      if( icor.eq.1 )then
        nvar = nvar+1
        varname(nvar) = 'pipc'
        vardesc(nvar) = 'diagnosed pi-prime: Coriolis component'
        varunit(nvar) = 'nondimensional'

        call getavg3d(pdiag(ib,jb,kb,4),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
      endif

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      if( wd_buoy.ge.1 )then

        nvar = nvar+1
        varname(nvar) = 'wb_buoy'
        vardesc(nvar) = 'w budget: buoyancy'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

        call getavg3d(wdiag(ib,jb,kb,wd_buoy),icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif

      !c-c-c-c-c-c-c-c-c-c

      IF( pdcomp )THEN

        nvar = nvar+1
        varname(nvar) = 'pgradb'
        vardesc(nvar) = 'vert pres grad: buoyancy component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do j=1,nj
          do i=1,ni
            wten(i,j,1) = 0.0
            wten(i,j,nk+1) = 0.0
          enddo
          enddo

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,1)-pdiag(i,j,k-1,1))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pgraddl'
        vardesc(nvar) = 'vert pres grad: linear dynamic component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,2)-pdiag(i,j,k-1,2))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

        nvar = nvar+1
        varname(nvar) = 'pgraddn'
        vardesc(nvar) = 'vert pres grad: nonlinear dynamic component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,3)-pdiag(i,j,k-1,3))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

        !-----

      if( icor.eq.1 )then
        nvar = nvar+1
        varname(nvar) = 'pgradc'
        vardesc(nvar) = 'vert pres grad: Coriolis component'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,4)-pdiag(i,j,k-1,4))*rdz*mf(i,j,k)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
      endif

        !-----

      ENDIF

      !c-c-c-c-c-c-c-c-c-c

      IF( pdcomp .and. wd_buoy.ge.1 )THEN

        nvar = nvar+1
        varname(nvar) = 'btot'
        vardesc(nvar) = 'total buoyant acceleration'
        varunit(nvar) = 'm/s/s'
        vargrid(nvar) = 'w'

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            wten(i,j,k) = -cp*(c2(i,j,k)*thv0(i,j,k)+c1(i,j,k)*thv0(i,j,k-1))  &
                             *(pdiag(i,j,k,1)-pdiag(i,j,k-1,1))*rdz*mf(i,j,k)  &
                          +wdiag(i,j,k,wd_buoy)
          enddo
          enddo
          enddo

        call getavg3d(wten,icrs,savg,ir,rr,navg,area,ruh,rvh)
        call write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)

      endif

      !c-c-c-c-c-c-c-c-c-c

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! close output file:

    id0:  IF( myid.eq.0 )THEN

  !-----------------------------------------------------------------------------

      of3:  &
      IF( output_format.eq.1 )THEN


        ! grads-format

        ! close binary output file:
        close(unit=fnums)
        close(unit=fnumw)

        ! write descriptor file:
        do nfile=1,2
          do i=1,maxstring
            string(i:i) = ' '
          enddo
          if( nfile.eq.1 )then
            string = 'cm1out_azimavg_s.ctl'
          elseif( nfile.eq.2 )then
            string = 'cm1out_azimavg_w.ctl'
          endif
          open(unit=66,file=string)
          if( nfile.eq.1 )then
            string = 'cm1out_azimavg_00%y4_s.dat'
          elseif( nfile.eq.2 )then
            string = 'cm1out_azimavg_00%y4_w.dat'
          endif
          write(66,101) string
          write(66,102)
          if( nfile.eq.1 )then
            if( outunits.eq.2 )then
              write(66,113) trim(cm1version)
            else
              write(66,103) trim(cm1version)
            endif
          elseif( nfile.eq.2 )then
            if( outunits.eq.2 )then
              write(66,213) trim(cm1version)
            else
              write(66,203) trim(cm1version)
            endif
          endif
          write(66,104) grads_undef
!!!          write(66,105) icrs,outunitconv*0.5*ddr,outunitconv*ddr
          write(66,115) icrs
          do i=1,icrs
            write(66,217) outunitconv*ravg(i)
          enddo
          write(66,106)
          if( nfile.eq.1 )then
            write(66,107) nk
            do k=1,nk
              write(66,217) outunitconv*zh(1,1,k)
            enddo
          elseif( nfile.eq.2 )then
            write(66,107) nk+1
            do k=1,nk+1
              write(66,217) outunitconv*zf(1,1,k)
            enddo
          endif
          write(66,108) nwritea
          ntmp = 0
          if( nfile.eq.1 )then
            do n=1,nvar
              if( vargrid(n).eq.'s' ) ntmp = ntmp+1
            enddo
          else
            do n=1,nvar
              if( vargrid(n).eq.'w' ) ntmp = ntmp+1
            enddo
          endif
          write(66,109) ntmp
          do n=1,nvar
            doit = .false.
            if( nfile.eq.1 .and. vargrid(n).eq.'s' ) doit = .true.
            if( nfile.eq.2 .and. vargrid(n).eq.'w' ) doit = .true.
            if( doit )then
              a1 = varname(n)
              a2 = vardesc(n)
              !---
              a16 = '                '
              nn = len(trim(varunit(n)))
              write(a16(2:15),314) varunit(n)
              write(a16(1:1),301 )       '('
              write(a16(nn+2:nn+2),301 ) ')'
              !---
              if(     nfile.eq.1 )then
                write(66,110) a1(1:12),varlvls(n),a2(1:40),a16
              elseif( nfile.eq.2 )then
                write(66,110) a1(1:12),(varlvls(n)+1),a2(1:40),a16
              endif
            endif
          enddo
          write(66,111)
          close(unit=66)
        enddo

 301    format(a1)
 314    format(a14)

 101    format('dset ^',a)
 102    format('options template')
 103    format('title CM1 azimuthally averaged output on scalar levels, using version ',a,'; x dim is radius (km), z dim is height ASL (km); time is generic, see variable mtime for actual times')
 203    format('title CM1 azimuthally averaged output on w levels, using version ',a,'; x dim is radius (km), z dim is height ASL (km); time is generic, see variable mtime for actual times')
 113    format('title CM1 azimuthally averaged output on scalar levels, using version ',a,'; x dim is radius (meters), z dim is height ASL (meters); time is generic, see variable mtime for actual times')
 213    format('title CM1 azimuthally averaged output on w levels, using version ',a,'; x dim is radius (meters), z dim is height ASL (meters); time is generic, see variable mtime for actual times')
 104    format('undef ',f10.1)
! 105    format('xdef ',i6,' linear ',f13.6,1x,f13.6)
 105    format('xdef ',i6,' linear ',es14.6e2,1x,es14.6e2)
 115    format('xdef ',i6,' levels')
 106    format('ydef 1 linear 0 1')
 107    format('zdef ',i6,' levels')
! 217    format(2x,f13.6)
 217    format(2x,es14.6e2)
 108    format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
 109    format('vars ',i6)
 110    format(a12,2x,i6,' 99 ',a40,1x,a16)
 111    format('endvars')

  !-----------------------------------------------------------------------------

#ifdef NETCDF
      ELSEIF( output_format.eq.2 )THEN  of3

        call     writeazim_nc(wloop,icrs,nvar,varmax,varname,vardesc,varunit,varlvls,vargrid,ncid,mtime,zh,zf,dat2(1,1),dum1(ib,jb,kb),dum2(ib,jb,kb))

#endif
      ENDIF  of3

    ENDIF  id0

    if( myid.eq.0 ) print *
    if( myid.eq.0 ) print *,' ..... end azimavg code ..... '
    if( myid.eq.0 ) print *

  ENDDO  bigwloop

    end subroutine azimavg


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getavg2d(s2d,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je) :: s2d
    integer, intent(in) :: icrs
    double precision, intent(inout), dimension(icrs) :: savg2d
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh

    integer :: i,j
    real :: ds

    ds = dx*dy

    savg2d = 0.0

    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg2d(ir(i,j)) = savg2d(ir(i,j)) + s2d(i,j)*(ds*(ruh(i)*rvh(j)))
      endif
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg2d,icrs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    do i=1,icrs
      if( navg(i).eq.0 )then
        savg2d(i) = grads_undef
      else
        savg2d(i) = savg2d(i)/area(i)
      endif
    enddo
   
    end subroutine getavg2d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getavg3d(s,icrs,savg,ir,rr,navg,area,ruh,rvh)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s
    integer, intent(in) :: icrs
    double precision, dimension(icrs,nk+1) :: savg
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh

    integer :: i,j,k
    real :: ds

    ds = dx*dy

    savg = 0.0

    do k=1,nk+1
    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg(ir(i,j),k) = savg(ir(i,j),k) + s(i,j,k)*(ds*(ruh(i)*rvh(j)))
      endif
    enddo
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg,icrs*(nk+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    do k=1,nk+1
    do i=1,icrs
      if( navg(i).eq.0 )then
        savg(i,k) = grads_undef
      else
        savg(i,k) = savg(i,k)/area(i)
      endif
    enddo
    enddo

    end subroutine getavg3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getsum3d(s,icrs,savg,ir,rr,navg,area,ruh,rvh)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s
    integer, intent(in) :: icrs
    double precision, dimension(icrs,nk+1) :: savg
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh

    integer :: i,j,k

    savg = 0.0

    do k=1,nk+1
    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg(ir(i,j),k) = savg(ir(i,j),k) + s(i,j,k)
      endif
    enddo
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg,icrs*(nk+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    do k=1,nk+1
    do i=1,icrs
      if( navg(i).eq.0 )then
        savg(i,k) = grads_undef
      endif
    enddo
    enddo

    end subroutine getsum3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getmax3d(s,icrs,savg,ir,rr,navg,area,ruh,rvh)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s
    integer, intent(in) :: icrs
    double precision, dimension(icrs,nk+1) :: savg
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh

    integer :: i,j,k

    savg = -1.0e30

    do k=1,nk+1
    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg(ir(i,j),k) = max( savg(ir(i,j),k) , s(i,j,k) )
      endif
    enddo
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg,icrs*(nk+1),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

    do k=1,nk+1
    do i=1,icrs
      if( navg(i).eq.0 )then
        savg(i,k) = grads_undef
      endif
    enddo
    enddo

    end subroutine getmax3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getmax2d(s2d,icrs,savg2d,ir,rr,navg,area,ruh,rvh)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je) :: s2d
    integer, intent(in) :: icrs
    double precision, dimension(icrs) :: savg2d
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh

    integer :: i,j,k

    savg2d = -1.0e30

    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg2d(ir(i,j)) = max( savg2d(ir(i,j)) , s2d(i,j) )
      endif
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg2d,icrs,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

    do i=1,icrs
      if( navg(i).eq.0 )then
        savg2d(i) = grads_undef
      endif
    enddo

    end subroutine getmax2d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getmin3d(s,icrs,savg,ir,rr,navg,area,ruh,rvh)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s
    integer, intent(in) :: icrs
    double precision, dimension(icrs,nk+1) :: savg
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh

    integer :: i,j,k

    savg = 1.0e30

    do k=1,nk+1
    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg(ir(i,j),k) = min( savg(ir(i,j),k) , s(i,j,k) )
      endif
    enddo
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg,icrs*(nk+1),MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
#endif

    do k=1,nk+1
    do i=1,icrs
      if( navg(i).eq.0 )then
        savg(i,k) = grads_undef
      endif
    enddo
    enddo

    end subroutine getmin3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getfrac3d(s,icrs,savg2d,savg,ir,rr,navg,area,ruh,rvh,op)
    use input
    use constants , only : grads_undef
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: s
    integer, intent(in) :: icrs
    double precision, dimension(icrs) :: savg2d
    double precision, dimension(icrs,nk+1) :: savg
    integer, intent(in), dimension(ib:ie,jb:je) :: ir
    real, intent(in), dimension(ib:ie,jb:je) :: rr
    integer, intent(in), dimension(icrs) :: navg
    double precision, intent(in), dimension(icrs) :: area
    real, intent(in), dimension(ib:ie) :: ruh
    real, intent(in), dimension(jb:je) :: rvh
    integer, intent(in) :: op     !  1 = use values >= 0
                                  !  2 = use values < 0

    integer :: i,j,k

    savg2d = 0.0
    savg = 0.0

    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
        savg2d(ir(i,j)) = savg2d(ir(i,j)) + dx*dy*ruh(i)*rvh(j)
      endif
    enddo
    enddo

    do k=1,nk+1
    do j=1,nj
    do i=1,ni
      if( ir(i,j).ge.1 .and. ir(i,j).le.icrs )then
      if( op.eq.1 .and. s(i,j,k).ge.0.0 )then
        savg(ir(i,j),k) = savg(ir(i,j),k) + dx*dy*ruh(i)*rvh(j)
      endif
      if( op.eq.2 .and. s(i,j,k).lt.0.0 )then
        savg(ir(i,j),k) = savg(ir(i,j),k) + dx*dy*ruh(i)*rvh(j)
      endif
      endif
    enddo
    enddo
    enddo

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg2d,icrs     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,savg,icrs*(nk+1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

    do k=1,nk+1
    do i=1,icrs
      if( navg(i).eq.0 )then
        savg(i,k) = grads_undef
      else
        savg(i,k) = savg(i,k)/max(1.0e-20,savg2d(i))
      endif
    enddo
    enddo

    end subroutine getfrac3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine write2d(icrs,savg2d,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    use input
#ifdef NETCDF
    use writeout_nc_module, only : disp_err
    use netcdf
#endif
    implicit none

    integer, intent(in) :: icrs
    double precision, intent(in), dimension(icrs) :: savg2d
    integer, intent(inout) :: arecs,arecw,nvar
    character(len=80), intent(in), dimension(varmax) :: varname,vardesc,varunit
    character(len=1), intent(in), dimension(varmax) :: vargrid

    integer :: i,varid,status
    real, dimension(1,1,1,1) :: var
      
  if( output_format.eq.1 .or. ( output_format.eq.2 .and. wloop.eq.2 ) )then

    if( myid.eq.0 )then

      print *,nvar,trim(varname(nvar))

      if( output_format.eq.1 )then

        ! grads-format file:

        ! note:  all 2d vars go into scalar output file
        write(fnums,rec=arecs) (sngl(savg2d(i)),i=1,icrs)
        arecs = arecs+1

#ifdef NETCDF
      elseif( output_format.eq.2 )then

        status = nf90_inq_varid(ncid,trim(varname(nvar)),varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error in write3d, varname = ',varname(nvar)
          print *,nf90_strerror(status)
          call stopcm1
        endif

        do i=1,icrs
          var = savg2d(i)
          call disp_err( nf90_put_var(ncid,varid,var,(/i,1,1/),(/1,1,1/)) , .true. )
        enddo
#endif

      endif

    endif

  endif

    end subroutine write2d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine write3d(icrs,savg,arecs,arecw,nvar,varname,vardesc,varunit,vargrid)
    use input
#ifdef NETCDF
    use writeout_nc_module, only : disp_err
    use netcdf
#endif
    implicit none

    integer, intent(in) :: icrs
    double precision, intent(in), dimension(icrs,nk+1) :: savg
    integer, intent(inout) :: arecs,arecw,nvar
    character(len=80), intent(in), dimension(varmax) :: varname,vardesc,varunit
    character(len=1), intent(in), dimension(varmax) :: vargrid

    integer :: i,k,kmax,varid,status
    real, dimension(1,1,1,1) :: var
      
  if( output_format.eq.1 .or. ( output_format.eq.2 .and. wloop.eq.2 ) )then

    if( myid.eq.0 )then

      print *,nvar,vargrid(nvar),' ',trim(varname(nvar))

      if( output_format.eq.1 )then
        ! grads-format file:

      if( vargrid(nvar).eq.'s' )then
        do k=1,nk
          write(fnums,rec=arecs) (sngl(savg(i,k)),i=1,icrs)
          arecs = arecs+1
        enddo
      elseif( vargrid(nvar).eq.'w' )then
        do k=1,nk+1
          write(fnumw,rec=arecw) (sngl(savg(i,k)),i=1,icrs)
          arecw = arecw+1
        enddo
      else
        print *,' 67833 '
        call stopcm1
      endif

#ifdef NETCDF
      elseif( output_format.eq.2 )then

        status = nf90_inq_varid(ncid,trim(varname(nvar)),varid)
        if(status.ne.nf90_noerr)then
          print *,'  Error in write3d, varname = ',varname(nvar)
          print *,nf90_strerror(status)
          call stopcm1
        endif

        if( vargrid(nvar).eq.'s' )then
          kmax = nk
        elseif( vargrid(nvar).eq.'w' )then
          kmax = nk+1
        else
          print *,'  67844 '
        endif

        do k=1,kmax
        do i=1,icrs
          var = savg(i,k)
          call disp_err( nf90_put_var(ncid,varid,var,(/i,1,k,1/),(/1,1,1,1/)) , .true. )
        enddo
        enddo
#endif
      endif

    endif

  endif


    end subroutine write3d


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine getcenter(nstep,mtime,rddr,xh,yh,zh,xfref,yfref,ugr,utmp,vgr,vtmp,pp3d,  &
                         nwritea,icenter,jcenter,xcenter,ycenter,domainlocx,domainlocy)
    use input
    use constants
#ifdef NETCDF
      use writeout_nc_module, only : writecenter_nc
#endif
#ifdef MPI
      use comm_module
      use mpi
#endif
    implicit none

      !----------------------------------------------------------------
      !     center-finding code:
      ! (note: for axisymmetric simulations, this code isn't used)
      !----------------------------------------------------------------

      integer, intent(in) :: nstep
      double precision, intent(in) :: mtime
      real, intent(in) :: rddr
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in   ), dimension(ib:ie+1,jb:je,kb:ke) :: ugr
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: utmp
      real, intent(in   ), dimension(ib:ie,jb:je+1,kb:ke) :: vgr
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vtmp
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: pp3d
      integer, intent(inout) :: nwritea,icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      double precision, intent(in) :: domainlocx,domainlocy

    integer :: i,j,k,ii,jj,ir,proc,jfoo,arecs,arecw,irec1,nr,kmax,nloop
    integer :: ictest,jctest,imin,imax,jmin,jmax,ninc
    real :: rr,vmax,mostmax,mostx,mosty,mostr,mostz,rmax,zmax,u,v,angle
    real :: xctest,yctest

    integer :: ipmin,jpmin
    real :: ppmin

    logical, parameter :: doit = .true.

#ifdef MPI
    integer status(MPI_STATUS_SIZE)
#endif

      integer, parameter :: centertype  =  1    ! 1 = gridpoint that produces the 
                                                !     maximum azimuthally averaged 
                                                !     tangential velocity
                                                ! 
                                                ! 2 = minimum perturbation pressure
                                                !     at lowest model level
                                                ! 

             ! assume rmw <= rcrit
             ! assume zmw <= zcrit
      real, parameter :: rcrit = 150000.0       ! max possible rmw (m)
      real, parameter :: zcrit =   4000.0       ! max possible zmw (m)
                                                ! (Note:  code uses last known location 
                                                !  of storm center as a starting point)

      real, parameter :: search_radius  =  30000.0   ! to reduce cost, search only within 
                                                     ! this radius (m) from previous known 
                                                     ! position

      integer, dimension(:), allocatable :: navg
      double precision, dimension(:), allocatable :: ravg
      double precision, dimension(:,:), allocatable :: davg
#ifdef MPI
      integer :: loc
      real, dimension(2) :: valmin
#endif

      integer, parameter :: nvars = 10

!-----------------------------------------------------------------------

  docenter:  IF(doit)THEN

    if(myid.eq.0) print *,'  Entering getcenter .... '


  centype:  &
  IF( centertype.eq.1 .or. centertype.eq.2 )THEN


    !
    ! 180620: For centype = 1 (max azim-avg tang. vel.):
    !          two-pass approach.  1st pass finds min surface pressure
    !                              2nd pass finds point with maximum azimth. avg. V
    !  (added to account for fast moving storms)
    !

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
    ! Pass 1:  find min surface pressure

#ifdef MPI
      call comm_2d_corner(pp3d(ib,jb,1))
#endif

      ppmin = 1.0e30

      do j=1,nj+1
      do i=1,ni+1
        if( pp3d(i,j,1).lt.ppmin )then
          ! nondim. pressure, lowest model level, vorticity points:
          ppmin = 0.25*( (pp3d(i-1,j-1,1)+pp3d(i  ,j  ,1)) &
                        +(pp3d(i-1,j  ,1)+pp3d(i  ,j-1,1)) )
          ipmin = i
          jpmin = j
        endif
      enddo
      enddo

#ifdef MPI
      valmin(1) = ppmin
      valmin(2) = myid
      call mpi_allreduce(mpi_in_place,valmin,1,MPI_2REAL,MPI_MINLOC,MPI_COMM_WORLD,ierr)

      ppmin = valmin(1)
      loc = nint(valmin(2))
      ipmin = ipmin + (myi1-1)
      jpmin = jpmin + (myj1-1)
      call MPI_BCAST(ipmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jpmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
#endif

      if( myid.eq.0 )then
      print *,'  myid,ppmin,ipmin,jpmin = ',myid,ppmin,ipmin,jpmin
      print *,'    xfref,yfref = ',xfref(ipmin),yfref(jpmin)
      endif

      xcenter = xfref(ipmin)
      ycenter = yfref(jpmin)

      icenter = ipmin
      jcenter = jpmin

    !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
    ! Pass 2:  find point with max azimuthally averaged tangential velocity


      nr = int( rcrit*rddr + 1 )

      kmax = 0
      do while( zh(1,1,kmax).le.zcrit .and. kmax.lt.nk )
        kmax = kmax+1
      enddo

    ! print xc,yc from last time:
    if(myid.eq.0) print *,'  Previous values: '
    if(myid.eq.0) print *,'    xcenter,ycenter = ',xcenter,ycenter
    if(myid.eq.0) print *,'    icenter,jcenter = ',icenter,jcenter
    if(myid.eq.0) print *,'  nr,kmax           = ',nr,kmax

    allocate( navg(nr) )
    navg = 0
    allocate( ravg(nr) )
    ravg = 0
    allocate( davg(nr,kmax) )
    davg = 0.0

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    do k=1,kmax
    do j=1,nj
    do i=1,ni
      utmp(i,j,k) = 0.5*( ugr(i,j,k) + ugr(i+1,j,k) )
      vtmp(i,j,k) = 0.5*( vgr(i,j,k) + vgr(i,j+1,k) )
    enddo
    enddo
    enddo

      mostmax = 0.0
      mostx = 0.0
      mosty = 0.0
      mostr = 0.0
      mostz = 0.0
      ictest = nx/2 + 1
      jctest = ny/2 + 1

    ct2:  &
    IF( centertype.eq.2 )THEN

      imin = icenter
      imax = icenter

      jmin = jcenter
      jmax = jcenter

      ninc  =  1

    ELSE  ct2


      imin = 1
      imax = nx

      do i=(nx+1),1,-1
        if( -(xfref(i)-xcenter).le.search_radius ) imin = i
      enddo

      do i=1,nx+1
        if(  (xfref(i)-xcenter).le.search_radius ) imax = i
      enddo

      jmin = 1
      jmax = ny

      do j=(ny+1),1,-1
        if( -(yfref(j)-ycenter).le.search_radius ) jmin = j
      enddo

      do j=1,ny+1
        if(  (yfref(j)-ycenter).le.search_radius ) jmax = j
      enddo

      ! to reduce cost, skip some grid points:
      !   (set ninc=1 to search every gridpoint)
      !   (set ninc=2 to search every other gridpoint)
      !   (set ninc=3 to search every third gridpoint)
      ninc  =  2
!!!      ninc  =  4

!    if( nstep.eq.0 )then
!      ! search entire domain:
!      imin = 1
!      imax = nx
!      jmin = 1
!      jmax = ny
!      ! use double the increment, to reduce cost
!      ninc = ninc*2
!    endif

    ENDIF  ct2

    if(myid.eq.0)then
      print *,'    search box: '
      print *,'      imin,imax,ninc = ',imin,imax,ninc
      print *,'      jmin,jmax,ninc = ',jmin,jmax,ninc
      print *,'      xmin,xmax,diff = ',xfref(imin),xfref(imax),xfref(imax)-xfref(imin)
      print *,'      ymin,ymax,diff = ',yfref(jmin),yfref(jmax),yfref(jmax)-yfref(jmin)
    endif

!    if( imin.lt.1 .or. imax.gt.nx .or.  &
!        jmin.lt.1 .or. jmax.gt.ny )then
!      print *
!      print *,'  too close to lateral boundary in subroutine getcenter '
!      print *
!      print *,'  56091 '
!      call stopcm1
!    endif

      jloop:  do jj=jmin,jmax,ninc
      iloop:  do ii=imin,imax,ninc
        navg = 0
        ravg = 0
        davg = 0.0
        xctest = xfref(ii)
        yctest = yfref(jj)
        if( nstep.eq.0 )then
          ! search entire domain:
          xcenter = xfref(ii)
          ycenter = yfref(jj)
        endif
        checkit: IF( sqrt( (xfref(ii)-xcenter)**2 + (yfref(jj)-ycenter)**2 ) .le. search_radius )THEN
          ! Calculate azim-avg vt:
          do j=1,nj
          do i=1,ni
            rr = sqrt( (xh(i)-xctest)**2 + (yh(j)-yctest)**2 )
            if( rr.lt.rcrit )then
              ir = 1 + int( rr*rddr )
              navg(ir)=navg(ir)+1
              if( abs(rr).lt.1.0e-4 )then
                u=0.0
                v=0.0
              else
                if( (xh(i)-xctest).ge.0.0 )then
                  angle =      asin( (yh(j)-yctest)/rr )
                else
                  angle = pi - asin( (yh(j)-yctest)/rr )
                endif
                do k=1,kmax
                  v=vtmp(i,j,k)*cos(angle)-utmp(i,j,k)*sin(angle)
                  davg(ir,k)=davg(ir,k)+v
                enddo
              endif
            endif
          enddo
          enddo
!-----------------------------------------------------------------------
!  communicate:
#ifdef MPI
        if( myid.eq.0 )then
          call MPI_REDUCE(mpi_in_place,davg,nr*kmax,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
          call MPI_REDUCE(mpi_in_place,navg,nr     ,mpi_integer         ,mpi_sum,0,mpi_comm_world,ierr)
        else
          call MPI_REDUCE(  davg      ,davg,nr*kmax,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
          call MPI_REDUCE(  navg      ,navg,nr     ,mpi_integer         ,mpi_sum,0,mpi_comm_world,ierr)
        endif
#endif
!-----------------------------------------------------------------------
        IF(myid.eq.0)THEN
          ! Search for maximum value for this center point:
          vmax = 0.0
          do i=1,nr
            ravg(i) = 1.0/dble(max(1,navg(i)))
          enddo
          do k=1,kmax
          do i=1,nr
            davg(i,k) = davg(i,k)*ravg(i)
            if( davg(i,k).gt.vmax )then
              vmax = davg(i,k)
              rmax = (i-1)*ddr + 0.5*ddr
              zmax = zh(1,1,k)
            endif
          enddo
          enddo
          if( vmax.gt.mostmax )then
            mostmax = vmax
            mostx = xctest
            mosty = yctest
            mostr = rmax
            mostz = zmax
            ictest = ii
            jctest = jj
          endif
        ENDIF
        ENDIF  checkit
      enddo  iloop
      enddo  jloop

      IF(myid.eq.0)THEN
        xcenter = mostx
        ycenter = mosty
        icenter = ictest
        jcenter = jctest
      ENDIF

  ELSE

    print *,'  87235:  unknown valye of centertype '
    call stopcm1

  ENDIF  centype

    if(myid.eq.0) print *,'  New values: '
    if(myid.eq.0) print *,'    xcenter,ycenter = ',xcenter,ycenter
    if(myid.eq.0) print *,'    icenter,jcenter = ',icenter,jcenter
    if(myid.eq.0) print *,'    vmax,rmax,zmax  = ',mostmax,mostr,mostz


  ELSE  docenter

    ! just use center of domain:
    icenter = nx/2 + 1
    jcenter = ny/2 + 1
    xcenter = minx + 0.5*(maxx-minx)
    ycenter = miny + 0.5*(maxy-miny)

  ENDIF  docenter


!-----------------------------------------------------------------------
!  more comms:

#ifdef MPI
    call MPI_BCAST(icenter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(jcenter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(xcenter,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ycenter,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

!-----------------------------------------------------------------------

  myid0:  &
  IF(myid.eq.0)THEN

    of1:  &
    IF( output_format.eq.1 )THEN
      ! grads format:

      print *,'  writing azimavg stats file'

        irec1 = 1 + nvars*(nwritea-1)

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_azimavg_stats.dat'
        print *,string
        open(unit=67,file=string,form='unformatted',access='direct',recl=4,status='unknown')

        write(67,rec=irec1) sngl(mtime)
        irec1=irec1+1
        write(67,rec=irec1) float(icenter)
        irec1=irec1+1
        write(67,rec=irec1) float(jcenter)
        irec1=irec1+1
        write(67,rec=irec1) xcenter
        irec1=irec1+1
        write(67,rec=irec1) ycenter
        irec1=irec1+1
        write(67,rec=irec1) mostmax
        irec1=irec1+1
        write(67,rec=irec1) mostr
        irec1=irec1+1
        write(67,rec=irec1) mostz
        irec1=irec1+1
        write(67,rec=irec1) sngl(domainlocx)
        irec1=irec1+1
        write(67,rec=irec1) sngl(domainlocy)
        irec1=irec1+1

        close(unit=67)


    !------------------------
    !  GrADS descriptor file:

      do i=1,maxstring
        string(i:i) = ' '
      enddo
      string = 'cm1out_azimavg_stats.ctl'
      open(unit=50,file=string,status='unknown')
      string = 'cm1out_azimavg_stats.dat'
      write(50,231) string
231   format('dset ^',a)
      write(50,202) trim(cm1version)
      write(50,203) grads_undef
      write(50,204)
      write(50,205)
      write(50,206)
      write(50,207) nwritea
      write(50,208) nvars
      write(50,209) 'mtime       ','model time (seconds since beginning of simulation)'
      write(50,209) 'icenter     ','i value of center point                           '
      write(50,209) 'jcenter     ','j value of center point                           '
      write(50,209) 'xcenter     ','x (m) value of center point                       '
      write(50,209) 'ycenter     ','y (m) value of center point                       '
      write(50,209) 'vmax        ','max azimuthally averaged tangential velocity (m/s)'
      write(50,209) 'rmax        ','radius of Vmax (m)                                '
      write(50,209) 'zmax        ','height (m ASL) of Vmax                            '
      write(50,209) 'domainlocx  ','x location of (center of) domain (m)              '
      write(50,209) 'domainlocy  ','y location of (center of) domain (m)              '
      write(50,210)
      close(unit=50)
202   format('title stats and info for CM1 azimuthally averaged output; using version ',a,'; time is generic, see variable mtime for actual times')
203   format('undef ',f10.1)
204   format('xdef 1 linear 0.0 1.0')
205   format('ydef 1 linear 0.0 1.0')
206   format('zdef 1 linear 0.0 1.0')
207   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
208   format('vars ',i4)
209   format(a12,' 0 99 ',a50)
210   format('endvars')

!-----------------------------------------------------------------------

#ifdef NETCDF
    ELSEIF( output_format.eq.2 )THEN  of1
      ! netcdf format:

    ! 210203: bug fix
      call       writecenter_nc(nwritea,ncid,mtime,icenter,jcenter,xcenter,ycenter,mostmax,mostr,mostz)

!-----------------------------------------------------------------------
#endif

    ENDIF  of1

  ENDIF  myid0

!-----------------------------------------------------------------------

    deallocate( navg )
    deallocate( ravg )
    deallocate( davg )

    if(myid.eq.0) print *,'  .... leaving getcenter '

    return
    end subroutine getcenter


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getpsmth(xh,yh,sigma,pp3d,dum1,xfref,yfref,icenter,jcenter,xcenter,ycenter,psmth)

      use input
      use constants
      use bc_module
#ifdef MPI
      use comm_module
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pp3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      integer, intent(inout) :: icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      real, intent(inout), dimension(ib:ie,jb:je) :: psmth

      integer :: i,j,k,n,ipmin,jpmin,kval
      real :: rmin,minval,zref
      real, dimension(cmp,jmp) :: west,newwest,east,neweast
      real, dimension(imp,cmp) :: south,newsouth,north,newnorth
      integer, dimension(8) :: reqs

#ifdef MPI
      integer :: loc
      real, dimension(2) :: valmin
#endif

    !---------------------------------------

    ! 210228: specify level to use
    ! (default is surface)

      zref = 0.0    ! level (ASL) to use

      minval = 1.0e30

      if( zref.gt.smeps )then
        do k=1,nk
          if( abs(sigma(k)-zref).lt.minval )then
            minval = abs(sigma(k)-zref)
            kval = k
          endif
        enddo
      else
        kval = 1
      endif

!!!      kval = nint(var2)

      if(myid.eq.0) print *,'  kval,zval = ',kval,sigma(kval)

    !---------------------------------------

      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = pp3d(i,j,kval)
      enddo
      enddo

      ! by default, run del^2 diffusion 100 times:
      do n=1,100

        call bc2d(dum1(ib,jb,1))
#ifdef MPI
        call comm_2d_start(dum1(ib,jb,1),west,newwest,east,neweast,   &
                                         south,newsouth,north,newnorth,reqs)
        call comm_2dew_end(dum1(ib,jb,1),west,newwest,east,neweast,reqs)
        call comm_2dns_end(dum1(ib,jb,1),south,newsouth,north,newnorth,reqs)
        call bcs2_2d(dum1(ib,jb,1))
        call comm_2d_corner(dum1(ib,jb,1))
#endif

        do j=1,nj
        do i=1,ni
          dum1(i,j,2) = 0.11111111*( dum1(i,j,1) + (                 &
                              ( (dum1(i+1,j  ,1)+dum1(i-1,j  ,1))    &
                               +(dum1(i  ,j+1,1)+dum1(i  ,j-1,1)) )  &
                             +( (dum1(i+1,j+1,1)+dum1(i-1,j+1,1))    &
                               +(dum1(i+1,j-1,1)+dum1(i-1,j-1,1)) ) ) )
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = dum1(i,j,2)
          psmth(i,j) = dum1(i,j,2)
        enddo
        enddo

      enddo

        call bc2d(dum1(ib,jb,1))
#ifdef MPI
        call comm_2d_start(dum1(ib,jb,1),west,newwest,east,neweast,   &
                                         south,newsouth,north,newnorth,reqs)
        call comm_2dew_end(dum1(ib,jb,1),west,newwest,east,neweast,reqs)
        call comm_2dns_end(dum1(ib,jb,1),south,newsouth,north,newnorth,reqs)
        call bcs2_2d(dum1(ib,jb,1))
        call comm_2d_corner(dum1(ib,jb,1))
#endif

    !---------------------------------------
    !  interpolate to zeta points, and find min value:

      do j=1,nj+1
      do i=1,ni+1
        dum1(i,j,2) = 0.25*( (dum1(i,j,1)+dum1(i-1,j-1,1)) &
                            +(dum1(i-1,j,1)+dum1(i,j-1,1)) )
      enddo
      enddo

      rmin = 1.0e30

      do j=1,nj
      do i=1,ni
        if( dum1(i,j,2).lt.rmin )then
          rmin = dum1(i,j,2)
          ipmin = i
          jpmin = j
        endif
      enddo
      enddo

#ifdef MPI
      valmin(1) = rmin
      valmin(2) = myid
      call mpi_allreduce(mpi_in_place,valmin,1,MPI_2REAL,MPI_MINLOC,MPI_COMM_WORLD,ierr)

      rmin = valmin(1)
      loc = nint(valmin(2))
      ipmin = ipmin + (myi-1)*ni
      jpmin = jpmin + (myj-1)*nj
      call MPI_BCAST(ipmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jpmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
#endif

      if( myid.eq.0 )then
      print *,'  myid,rmin,ipmin,jpmin = ',myid,rmin,ipmin,jpmin
      print *,'    xf,yf = ',xfref(ipmin),yfref(jpmin)
      endif

      xcenter = xfref(ipmin)
      ycenter = yfref(jpmin)

      icenter = ipmin
      jcenter = jpmin

      end subroutine getpsmth


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getpcent(xh,yh,xfref,yfref,pp3d,icenter,jcenter,xcenter,ycenter,psmth)

      use input
      use constants
      use bc_module
#ifdef MPI
      use comm_module
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pp3d
      integer, intent(inout) :: icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      real, intent(in), dimension(ib:ie,jb:je) :: psmth

      integer :: i,j
      real :: ppmin,ppmax,pcent,pminx,pmaxx,pminy,pmaxy
      real :: aa,minval
      double precision :: cx,cy,ca

      integer, parameter :: type = 2    ! 1 = center of pcent
                                        ! 2 = centroid of pcent

      real, parameter :: radmax = 300000.0   ! max radius (m) over which pcent is calculated
                                             ! (for type=2 only)

  !---------------------------------------------------------------------

      ppmin =  1.0e30
      ppmax = -1.0e30

      do j=1,nj
      do i=1,ni
        ppmin = min(ppmin,pp3d(i,j,1))
        ppmax = max(ppmax,pp3d(i,j,1))
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(mpi_in_place,ppmin,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mpi_in_place,ppmax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

      pcent = 0.5*(ppmin+ppmax)

  !---------------------------------------------------------------------

    IF( type.eq.1 )THEN

      pminx =  1.0e30
      pmaxx = -1.0e30

      pminy =  1.0e30
      pmaxy = -1.0e30

      do j=1,nj
      do i=1,ni
        if( pp3d(i,j,1).lt.pcent )then
          pminx = min(pminx,xh(i))
          pmaxx = max(pmaxx,xh(i))
          pminy = min(pminy,yh(j))
          pmaxy = max(pmaxy,yh(j))
        endif
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(mpi_in_place,pminx,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mpi_in_place,pmaxx,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mpi_in_place,pminy,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mpi_in_place,pmaxy,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

      xcenter = 0.5*(pminx+pmaxx)
      ycenter = 0.5*(pminy+pmaxy)

    ENDIF

  !---------------------------------------------------------------------

    IF( type.eq.2 )THEN

      if(myid.eq.0) print *,'  ... centroid ... '

!!!      do j=1,nj
!!!      do i=1,ni
!!!        out2d(i,j,1) = 0.0
!!!      enddo
!!!      enddo

      ca = 0.0
      cx = 0.0
      cy = 0.0

      do j=1,nj
      do i=1,ni
!!!        if( pp3d(i,j,1).lt.pcent )then
        if( sqrt( (xh(i)-centerx)**2 + (yh(j)-centery)**2 ) .lt. radmax )then
!!!          aa = (dx*ruh(i))*(dy*rvh(j))
          aa = psmth(i,j)
          ca = ca + aa
          cx = cx + aa*xh(i)
          cy = cy + aa*yh(j)
!!!          out2d(i,j,1) = 1.0
        endif
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(mpi_in_place,ca,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mpi_in_place,cx,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mpi_in_place,cy,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      if( abs(ca).gt.smeps )then
        cx = cx/ca
        cy = cy/ca
      else
        stop 19999
        cx = centerx
        cy = centery
      endif

      minval =  1.0e30
      do i=1,nx+1
        if( abs(xfref(i)-cx).lt.minval )then
          minval = abs(xfref(i)-cx)
          icenter = i
        endif
      enddo

      minval =  1.0e30
      do j=1,ny+1
        if( abs(yfref(j)-cy).lt.minval )then
          minval = abs(yfref(j)-cy)
          jcenter = j
        endif
      enddo

      if(myid.eq.0)then
        print *,'  cx,cy = ',cx,cy
        print *,'  xf,yf = ',xfref(icenter),yfref(jcenter)
        print *,'  i,j   = ',icenter,jcenter
      endif

      xcenter = cx
      ycenter = cy

    ENDIF

  !---------------------------------------------------------------------

      end subroutine getpcent


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine change_uvmove(u0,ua,u3d,v0,va,v3d,oldumove,oldvmove)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,v3d
      real, intent(in) :: oldumove,oldvmove

      integer :: i,j,k

      do k=1,nk
        do j=jb,je
        do i=ib,ie+1
          u0(i,j,k) = u0(i,j,k) + (oldumove-umove)
          ua(i,j,k) = ua(i,j,k) + (oldumove-umove)
          u3d(i,j,k) = u3d(i,j,k) + (oldumove-umove)
        enddo
        enddo
        do j=jb,je+1
        do i=ib,ie
          v0(i,j,k) = v0(i,j,k) + (oldvmove-vmove)
          va(i,j,k) = va(i,j,k) + (oldvmove-vmove)
          v3d(i,j,k) = v3d(i,j,k) + (oldvmove-vmove)
        enddo
        enddo
      enddo

      end subroutine change_uvmove


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine update_adapt_move(nstep,mtime,xh,yh,xfref,yfref,ug,vg,psmth,     &
                                   dum1,pp3d,u0,ua,u3d,v0,va,v3d,mvrec,nwritemv,  &
                                   icenter,jcenter,xcenter,ycenter,               &
                                   domainlocx,domainlocy)
      use input
      use constants
#ifdef NETCDF
      use writeout_nc_module, only : writemove_nc
#endif
      implicit none

      integer, intent(in) :: nstep
      double precision, intent(in) :: mtime
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref
      real, intent(inout), dimension(kb:ke) :: ug,vg
      real, intent(inout), dimension(ib:ie,jb:je) :: psmth
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,pp3d
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,v3d
      integer, intent(inout) :: mvrec
      integer, intent(in)    :: nwritemv
      integer, intent(inout) :: icenter,jcenter
      real, intent(inout) :: xcenter,ycenter
      double precision, intent(in) :: domainlocx,domainlocy

      integer :: k
      real :: oldumove,oldvmove,term2

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !cc   Adaptive domain-moving option for TCs  ccccccccccccccccccccccccccc
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ! Experimental (seems to work well, though)
        !
        ! Uses a smoothed field of pressure at the lowest model level to determine
        ! the position of a TC.

        ! skip if this is first time step
        if( nstep.ge.1 )then

          ! 220313: NOTE: moved getpsmth and getpcent to cm1.F

          ! attempt to move toward center of domain by next adaptmove time:
          oldumove = umove
          oldvmove = vmove

          ! 220313: increased timescale from 1.0 to 1.5 times adapt_move_frq
!!!          umove = -(centerx-xcenter)/(1.5*adapt_move_frq)
!!!          vmove = -(centery-ycenter)/(1.5*adapt_move_frq)

          term2 = 0.5

          umove = oldumove -term2*(centerx-xcenter)/(adapt_move_frq)  &
                           +(1.0-term2)*(umovem2-umovem1)

          vmove = oldvmove -term2*(centery-ycenter)/(adapt_move_frq)  &
                           +(1.0-term2)*(vmovem2-vmovem1)

          umovem2 = umovem1
          vmovem2 = vmovem1

          umovem1 = oldumove
          vmovem1 = oldvmove

          if( myid.eq.0 ) print *,'  adjusting umove,vmove:',umove,vmove

          call change_uvmove(u0,ua,u3d,v0,va,v3d,oldumove,oldvmove)

            do k=1,nk
              ug(k) = ug(k) + (oldumove-umove)
              vg(k) = vg(k) + (oldvmove-vmove)
            enddo

        endif

  !-------------------------------------------------------
  ! writeout:

    mid0:  &
    IF( myid.eq.0 )THEN

      of2:  &
      IF( output_format.eq.1 )THEN

          open(unit=82,file='cm1out_move.dat',form='unformatted',access='direct',recl=4)

          write(82,rec=mvrec) sngl(mtime)
          mvrec = mvrec+1
          write(82,rec=mvrec) xcenter
          mvrec = mvrec+1
          write(82,rec=mvrec) ycenter
          mvrec = mvrec+1
          write(82,rec=mvrec) real(icenter)
          mvrec = mvrec+1
          write(82,rec=mvrec) real(jcenter)
          mvrec = mvrec+1
          write(82,rec=mvrec) umove
          mvrec = mvrec+1
          write(82,rec=mvrec) vmove
          mvrec = mvrec+1
          write(82,rec=mvrec) sngl(domainlocx)
          mvrec = mvrec+1
          write(82,rec=mvrec) sngl(domainlocy)
          mvrec = mvrec+1

          close(unit=82)

          open(unit=50,file='cm1out_move.ctl',status='unknown')
          write(50,231)
231       format('dset ^cm1out_move.dat')
          write(50,202) trim(cm1version)
          write(50,203) grads_undef
          write(50,204)
          write(50,205)
          write(50,206)
          write(50,207) nwritemv
          write(50,208) 9
          write(50,209) 'mtime       ','model time (seconds since beginning of simulation)'
          write(50,209) 'xcenter     ','x (m) value of center point                       '
          write(50,209) 'ycenter     ','y (m) value of center point                       '
          write(50,209) 'icenter     ','i value of center point                           '
          write(50,209) 'jcenter     ','j value of center point                           '
          write(50,209) 'umove       ','umove (m/s)                                       '
          write(50,209) 'vmove       ','vmove (m/s)                                       '
          write(50,209) 'domainlocx  ','x location of (center of) domain (m)              '
          write(50,209) 'domainlocy  ','y location of (center of) domain (m)              '
          write(50,210)
          close(unit=50)
202       format('title CM1 adaptive move info, using version ',a,'; time is generic, see variable mtime for actual times')
203       format('undef ',f10.1)
204       format('xdef 1 linear 0.0 1.0')
205       format('ydef 1 linear 0.0 1.0')
206       format('zdef 1 linear 0.0 1.0')
207       format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
208       format('vars ',i4)
209       format(a12,' 0 99 ',a50)
210       format('endvars')


#ifdef NETCDF
      ELSEIF( output_format.eq.2 )THEN  of2

        call     writemove_nc(nwritemv,ncid,mtime,icenter,jcenter,xcenter,ycenter,umove ,vmove )

#endif
      ENDIF  of2

    ENDIF  mid0

      end subroutine update_adapt_move


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getrangle(xh,yh,xcenter,ycenter, ir, rr, angle,icrs,crref,xhref,xfref)
      use input
      use constants
      implicit none

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(inout) :: xcenter,ycenter
      integer, intent(inout), dimension(ib:ie,jb:je) :: ir
      real, intent(inout), dimension(ib:ie,jb:je) :: rr,angle
      integer, intent(in) :: icrs
      real, intent(in), dimension(icrs) :: crref
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xhref,xfref

      integer :: i,j,nn,ii,iimin
      real :: xsave,ysave,rddr,rrtmp
      real :: diff,diffmin

      xsave = xcenter
      ysave = ycenter

      ir = 0
      rr = 1.0e30
      angle = 0.0
      rddr = 1.0/ddr

    do nn = 1 , 1

      if( nn.eq.1 )then
        xcenter = xsave
        ycenter = ysave
      endif
      if( nn.eq.2 )then
        xcenter = xsave + (maxx-minx)
        ycenter = ysave + (maxy-miny)
      endif
      if( nn.eq.3 )then
        xcenter = xsave - (maxx-minx)
        ycenter = ysave - (maxy-miny)
      endif
      if( nn.eq.4 )then
        xcenter = xsave + (maxx-minx)
        ycenter = ysave - (maxy-miny)
      endif
      if( nn.eq.5 )then
        xcenter = xsave - (maxx-minx)
        ycenter = ysave + (maxy-miny)
      endif
      if( nn.eq.6 )then
        xcenter = xsave + (maxx-minx)
        ycenter = ysave
      endif
      if( nn.eq.7 )then
        xcenter = xsave - (maxx-minx)
        ycenter = ysave
      endif
      if( nn.eq.8 )then
        xcenter = xsave
        ycenter = ysave + (maxy-miny)
      endif
      if( nn.eq.9 )then
        xcenter = xsave
        ycenter = ysave - (maxy-miny)
      endif

      do j=1,nj
      do i=1,ni
        rrtmp = sqrt( (xh(i)-xcenter)**2 + (yh(j)-ycenter)**2 )
        if( rrtmp.lt.rr(i,j) )then
          rr(i,j) = rrtmp
          if( (xh(i)-xcenter).ge.0.0 )then
            angle(i,j) =      asin( (yh(j)-ycenter)/max(smeps,rr(i,j)) )
          else
            angle(i,j) = pi - asin( (yh(j)-ycenter)/max(smeps,rr(i,j)) )
          endif
        endif
      enddo
      enddo

    enddo

      xcenter = xsave
      ycenter = ysave

      do j=1,nj
      do i=1,ni
        rrtmp = sqrt( (xh(i)-xcenter)**2 + (yh(j)-ycenter)**2 )
        diffmin = 1.0e30
        do ii=1,icrs
          diff = rrtmp-crref(ii)
          if( abs(diff).lt.diffmin )then
            diffmin = abs(diff)
            iimin = ii
          endif
        enddo
!!!        ir(i,j) = iimin
        if( rr(i,j).lt.xfref(nx/2+icrs+1) )   ir(i,j) = iimin
      enddo
      enddo

      end subroutine getrangle


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  END MODULE azimavg_module
