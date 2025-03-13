  MODULE hifrq_module

  implicit none

  private
  public :: writeout_hifrq

  integer, parameter :: fnum = 67
  integer, parameter :: varmax = 10000

  integer :: ncid,wloop

  CONTAINS

!-------------------------------------------------------------------------------
!
!  hifrq:  code to write some output at a relatively high frequency
!          (compared to standard 3d output)
!
!-------------------------------------------------------------------------------

    subroutine writeout_hifrq(                                                 &
                   nstep,mtime,nwriteh,qname,dt,dosfcflx,                      &
                   icrs,icenter,jcenter,xcenter,ycenter,                       &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,            &
                   yh,vh,rvh,yf,vf,rvf,                                        &
                   xfref,xhref,yfref,yhref,rds,sigma,rdsf,sigmaf,              &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,                   &
                   rho0s,pi0s,prs0s,rth0s,                                     &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,u0,v0,                  &
                   qi0,rr0,rf0,rrf0,                                           &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,        &
                   tsk,znt,ust,stau,thflux,qvflux,cd,ch,cq,u1,v1,s1,xland,psfc,tlh,psmth, &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                    &
                   divx,rho,rr,rf,prs,t11,t12,t13,t22,t23,t33,                 &
                   rru,u3d,ugr ,utmp ,                                         &
                   rrv,v3d,vgr ,vtmp ,                                         &
                   rrw,w3d,wten,wtmp ,                                         &
                   pp3d,ppten,sten,th3d,sadv,thten,thten1,                     &
                   q3d,qten,kmh,kmv,khh,khv,tkea,tke3d,tketen,                 &
                   nm,defv,defh,dissten,                                       &
                   thpten,qvpten,qcpten,qipten,upten,vpten,xkzh,xkzq,xkzm,     &
                   rain,raincv,u10,v10,s10,br,brcr,hpbl,t2,q2,hfx,qfx,prate,cm0, &
                   swten,lwten,cldfra,qke,                                     &
                   lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc,        &
                   swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,        &
                   lwcf,swcf,wprof,dumk1,dumk2,                                &
                   tdiag,qdiag,udiag,vdiag,wdiag,pdiag,out2d,out3d,getdbz,getvt, &
                   cir,crr,cangle,                                             &
                   myi1p,myi2p,myj1p,myj2p,kbdy,bndy,                          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,flag,               &
                   dat1,dat2,dat3,reqt,reqs_s)
        ! end hifrq
    use input
    use constants
    use cm1libs , only : rslf,rsif
!!!    use adv_module , only : advs
    use bc_module
    use ib_module
#ifdef NETCDF
    use writeout_nc_module , only : writehifrq_nc
#endif
#ifdef MPI
    use comm_module
    use mpi
#endif
    implicit none

    integer, intent(in) :: nstep
    double precision, intent(in) :: mtime
    integer, intent(inout) :: nwriteh
    character(len=3), intent(in), dimension(maxq) :: qname
    real, intent(inout) :: dt
    logical, intent(in) :: dosfcflx
    integer, intent(in) :: icrs
    integer, intent(inout) :: icenter,jcenter
    real,    intent(inout) :: xcenter,ycenter
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
    real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,stau,thflux,qvflux,cd,ch,cq,  &
                                                u1,v1,s1,xland,psfc,tlh,psmth
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rr,rf,prs
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
    real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,ugr,utmp
    real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d,vgr,vtmp
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d,wten,wtmp
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: pp3d,ppten,sten
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th3d,sadv,thten,thten1
    real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d,qten
    real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
    real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,dissten
    real, intent(in), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten
    real, intent(in), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: xkzh,xkzq,xkzm
    real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain,raincv
    real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,br,brcr,hpbl,t2,q2,hfx,qfx
    real, intent(in), dimension(ib:ie,jb:je) :: prate,cm0
    real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,cldfra
    real, intent(in), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: qke
    real, intent(in), dimension(ibr:ier,jbr:jer) :: lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc
    real, intent(in), dimension(ibr:ier,jbr:jer) :: swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc
    real, intent(inout), dimension(ibr:ier,jbr:jer) :: lwcf,swcf
    real, intent(in),    dimension(kb:ke) :: wprof
    double precision, intent(inout), dimension(kb:ke) :: dumk1,dumk2
    real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
    real, intent(in) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
    real, intent(in) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
    real, intent(in) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag
    real, intent(inout), dimension(ib2d:ie2d,jb2d:je2d,nout2d) :: out2d
    real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
    logical, intent(in) :: getdbz,getvt
    integer, intent(in), dimension(ib:ie,jb:je) :: cir
    real, intent(in), dimension(ib:ie,jb:je) :: crr,cangle
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p
    integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
    real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
    real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
    logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
    real, intent(inout), dimension(d3is,d3js) :: dat1
    real, intent(inout), dimension(d2is,d2js) :: dat2
    real, intent(inout), dimension(d3is,d3js,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, dimension(rmp) :: reqs_s

    !---------------------------------------------------------------------------

    integer :: nvar,orec,i,j,k,n,nn,nlevels,diffit,wloopmax
    real :: delz,delzmin,pint,plast,qvs,tx,havg,ucrit
    logical :: doit,domse
    double precision :: sum,weps,bfoo

    integer, dimension(nk) :: klev
    real, dimension(nk) :: zlev

    character(len=80), dimension(varmax) :: varname,vardesc,varunit
    integer, dimension(varmax) :: varlvls
    character(len=80) :: a1,a2
    character(len=16) :: a16

    !---------------------------------------------------------------------------

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  Entering writeout_hifrq ... '
      if( myid.eq.0 ) print *

      nlevels = 5
      zlev(1) =     0.0
      zlev(2) =   500.0
      zlev(3) =  1000.0
      zlev(4) =  2500.0
      zlev(5) =  5000.0

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  levels for hifrq output: '

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'    ... requested levels ...'
      do n=1,nlevels
        if(myid.eq.0) print *,'    n,z = ',n,zlev(n)
      enddo

      ! find level that is closest to the requested level:
      do n=1,nlevels
        if( zlev(n).le.maxz )then
          delzmin = 1.0e30
          ! search from top, so that "ties" go to higher level
          do k=nk,1,-1
            delz = abs(sigma(k)-zlev(n))
            if( delz.lt.delzmin )then
              delzmin = delz
              klev(n) = k
            endif
          enddo
        else
          klev(n) = nk+1
        endif
      enddo

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'    ... first pass ...'
      do n=1,nlevels
        if(myid.eq.0) print *,'    k,z = ',klev(n),sigma(klev(n))
      enddo

      ! make sure the same level is not selected more than once:
      do n=2,nlevels
        klev(n) = max( klev(n) , klev(n-1)+1 )
      enddo

      ! check that levels do not exceed nk:
      n = 1
      do while( n .le. nlevels )
        if( sigma(n).ge.maxz .or. klev(n).gt.nk )then
          nlevels = n-1
        endif
        n = n+1
      enddo

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'    ... after checks ...'

      do n=1,nlevels
        if(myid.eq.0) print *,'    k,z = ',klev(n),sigma(klev(n))
      enddo

      if( nlevels.le.0 )then
        if(myid.eq.0)then
        print *
        print *,'  nlevels = 0,  no levels found '
        print *
        print *,'   stopping model in hifrq.F .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

    !---------------------------------------------------------------------------
    ! open file:

      if( myid.eq.0 ) print *,'  nwriteh = ',nwriteh

      IF( output_format.eq.1 )THEN

        ! grads-format

        string = 'cm1out_hf_XXXXXX.dat'

        write(string(11:16),201) nwriteh
201     format(i6.6)

        if( myid.eq.0 )then
          print *,string
          close(unit=fnum)
          open(unit=fnum,file=string,form='unformatted',access='direct',recl=4*nx*ny)
        endif

        orec = 1

        wloopmax = 1

      ELSEIF( output_format.eq.2 )THEN

        ! netcdf-format
#ifdef NETCDF

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_hf_XXXXXX.nc'
        write(string(11:16),201) nwriteh
        if( myid.eq.0 ) print *,string

        wloopmax = 2

#endif
      ELSE

        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  Invalid value for output_format in domaindiag '
        print *
        call stopcm1

      ENDIF

    !---------------------------------------------------------------------------

  bigwloop:  &
  DO wloop = 1 , wloopmax

      if( myid.eq.0 ) print *,'  wloop = ',wloop

      nvar = 0

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'mtime'
      vardesc(nvar) = 'model time (seconds since beginning of simulation)'
      varunit(nvar) = 's'
      varlvls(nvar) = 0

      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = mtime
      enddo
      enddo

      call writehf(1,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'rain'
      vardesc(nvar) = 'accumulated surface rainfall'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call writehf(1,rain(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. nrain.ge.2 )THEN

      nvar = nvar+1
      varname(nvar) = 'rain2'
      vardesc(nvar) = 'accumulated surface rainfall (translated)'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call writehf(1,rain(ib,jb,2),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. cuparam.ge.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'raincv'
      vardesc(nvar) = 'accumulated surface rainfall (convective, param.)'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call writehf(1,raincv(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. cuparam.ge.1 .and. nrain.ge.2 )THEN

      nvar = nvar+1
      varname(nvar) = 'raincv2'
      vardesc(nvar) = 'translated accumulated surface rainfall (convective, param.)'
      varunit(nvar) = 'cm'
      varlvls(nvar) = 0

      call writehf(1,raincv(ib,jb,2),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 )THEN

      nvar = nvar+1
      varname(nvar) = 'prate'
      vardesc(nvar) = 'surface precipitation rate'
      varunit(nvar) = 'kg/m2/s'
      varlvls(nvar) = 0

      call writehf(1,prate,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'thflux'
      vardesc(nvar) = 'surface potential temperature flux'
      varunit(nvar) = 'K m/s'
      varlvls(nvar) = 0

      call writehf(1,thflux,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'qvflux'
      vardesc(nvar) = 'surface water vapor mixing ratio flux'
      varunit(nvar) = 'g/g m/s'
      varlvls(nvar) = 0

      call writehf(1,qvflux,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'hfx'
      vardesc(nvar) = 'surface sensible heat flux'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,hfx,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'qfx'
      vardesc(nvar) = 'surface latent heat flux'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,qfx,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'u10'
    if( imove.eq.1 )then
      vardesc(nvar) = 'u comp. of 10m wind speed (ground-rel.)'
    else
      vardesc(nvar) = 'u comp. of 10m wind speed'
    endif
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call writehf(1,u10,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'v10'
    if( imove.eq.1 )then
      vardesc(nvar) = 'v comp. of 10m wind speed (ground-rel.)'
    else
      vardesc(nvar) = 'v comp. of 10m wind speed'
    endif
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call writehf(1,v10,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

  !---------------------------------------------------------------------

    IF( sfcmodel.ge.1 )THEN

    ! radial and tangential velocity:
    doit = .true.

    docylsfc:  &
    if( doit )then
    if( doazimavg .or. do_adapt_move )then

              ! dum1 = u at s pts
              ! dum2 = v at s pts
              ! dum3 = ur at s pts
              ! dum4 = vt at s pts

        k = 1
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = u10(i,j)
          dum2(i,j,k) = v10(i,j)
          if( abs(crr(i,j)).lt.1.0e-4 )then
            dum3(i,j,k)=0.0
            dum4(i,j,k)=0.0
          else
            dum3(i,j,k)=dum2(i,j,k)*sin(cangle(i,j))+dum1(i,j,k)*cos(cangle(i,j))
            dum4(i,j,k)=dum2(i,j,k)*cos(cangle(i,j))-dum1(i,j,k)*sin(cangle(i,j))
          endif
        enddo
        enddo

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'urad10'
      vardesc(nvar) = 'rad. comp. of 10m wind speed (ground-rel.)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call writehf(1,dum3(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'vtan10'
      vardesc(nvar) = 'tang. comp. of 10m wind speed (ground-rel.)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call writehf(1,dum4(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

    endif
    endif  docylsfc

    ENDIF

  !---------------------------------------------------------------------

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 't2'
      vardesc(nvar) = 'diagnostic 2m temperature'
      varunit(nvar) = 'K'
      varlvls(nvar) = 0

      call writehf(1,t2,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'q2'
      vardesc(nvar) = 'diagnostic 2m mixing ratio'
      varunit(nvar) = 'kg/kg'
      varlvls(nvar) = 0

      call writehf(1,q2,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'ust'
      vardesc(nvar) = 'friction velocity (u star)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

      call writehf(1,ust,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'stau'
      vardesc(nvar) = 'surface stress'
      varunit(nvar) = 'm2/s2'
      varlvls(nvar) = 0

      call writehf(1,stau,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'znt'
      vardesc(nvar) = 'surface roughness length'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call writehf(1,znt,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sfcmodel.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'psfc'
      vardesc(nvar) = 'surface pressure'
      varunit(nvar) = 'Pa'
      varlvls(nvar) = 0

      call writehf(1,psfc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    if( do_adapt_move )then
      nvar = nvar+1
      varname(nvar) = 'psmth'
      vardesc(nvar) = 'psmth'
      varunit(nvar) = 'Pa'
      varlvls(nvar) = 0

      call writehf(1,psmth,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    endif

      !c-c-c-c-c-c-c-c-c-c

    IF( output_dbz.eq.1 .and. qd_dbz.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'cref'
      vardesc(nvar) = 'composite reflectivity'
      varunit(nvar) = 'dbz'
      varlvls(nvar) = 0

      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = -1000.0
      enddo
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = max( dum1(i,j,1) , qdiag(i,j,k,qd_dbz) )
      enddo
      enddo
      enddo

      call writehf(1,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 )THEN
      ! water paths, precipitable water:

      do j=1,nj
      do i=1,ni
        dum1(i,j,1) = 0.0
        dum1(i,j,2) = 0.0
        dum1(i,j,3) = 0.0
        dum1(i,j,4) = 0.0
      enddo
      enddo

      do k=1,nk
        if( nqv.ge.1 )then
          do j=1,nj
          do i=1,ni
            dum1(i,j,1) = dum1(i,j,1) + rho(i,j,k)*q3d(i,j,k,nqv)*dz*rmh(i,j,k)
            qvs = rslf(prs(i,j,k),(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
            dum1(i,j,4) = dum1(i,j,4) + rho(i,j,k)*qvs*dz*rmh(i,j,k)
          enddo
          enddo
        endif
        if( nqc.ge.1 )then
          do j=1,nj
          do i=1,ni
            dum1(i,j,2) = dum1(i,j,2) + rho(i,j,k)*q3d(i,j,k,nqc)*dz*rmh(i,j,k)
          enddo
          enddo
        endif
        if( nqr.ge.1 )then
          do j=1,nj
          do i=1,ni
            dum1(i,j,2) = dum1(i,j,2) + rho(i,j,k)*q3d(i,j,k,nqr)*dz*rmh(i,j,k)
          enddo
          enddo
        endif
        if( nqs.ge.1 )then
          do j=1,nj
          do i=1,ni
            dum1(i,j,3) = dum1(i,j,3) + rho(i,j,k)*q3d(i,j,k,nqs)*dz*rmh(i,j,k)
          enddo
          enddo
        endif
        if( nqg.ge.1 )then
          do j=1,nj
          do i=1,ni
            dum1(i,j,3) = dum1(i,j,3) + rho(i,j,k)*q3d(i,j,k,nqg)*dz*rmh(i,j,k)
          enddo
          enddo
        endif
        if( nqi.ge.1 )then
          do j=1,nj
          do i=1,ni
            dum1(i,j,3) = dum1(i,j,3) + rho(i,j,k)*q3d(i,j,k,nqi)*dz*rmh(i,j,k)
          enddo
          enddo
        endif
      enddo

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'prw'
      vardesc(nvar) = 'water vapor path'
      varunit(nvar) = 'kg/m^2'
      varlvls(nvar) = 0

      call writehf(1,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'clwvi'
      vardesc(nvar) = 'condensed water path'
      varunit(nvar) = 'kg/m^2'
      varlvls(nvar) = 0

      call writehf(1,dum1(ib,jb,2),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'clivi'
      vardesc(nvar) = 'ice water path'
      varunit(nvar) = 'kg/m^2'
      varlvls(nvar) = 0

      call writehf(1,dum1(ib,jb,3),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'spwr'
      vardesc(nvar) = 'saturated water vapor path'
      varunit(nvar) = 'kg/m^2'
      varlvls(nvar) = 0

      call writehf(1,dum1(ib,jb,4),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

    if( idoles .and. cm1setup.eq.4 )then
      nvar = nvar+1
      varname(nvar) = 'cm0'
      vardesc(nvar) = 'cm0'
      varunit(nvar) = 'unitless'
      varlvls(nvar) = 0

      call writehf(1,cm0,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    endif

      !c-c-c-c-c-c-c-c-c-c

    IF( (cm1setup.eq.2.or.cm1setup.eq.4) .and. horizturb.eq.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'tlh'
      vardesc(nvar) = 'horiz lengthscale for turbulence scheme'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call writehf(1,tlh,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwupt'
      vardesc(nvar) = 'lw flux, upward, top of atmosphere (OLR)'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwupt,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwdnt'
      vardesc(nvar) = 'lw flux, downward, top of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwdnt,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwupb'
      vardesc(nvar) = 'lw flux, upward, bottom of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwupb,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwdnb'
      vardesc(nvar) = 'lw flux, downward, bottom of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwdnb,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swupt'
      vardesc(nvar) = 'sw flux, upward, top of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swupt,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swdnt'
      vardesc(nvar) = 'sw flux, downward, top of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swdnt,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swupb'
      vardesc(nvar) = 'sw flux, upward, bottom of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swupb,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swdnb'
      vardesc(nvar) = 'sw flux, downward, bottom of atmosphere'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swdnb,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwuptc'
      vardesc(nvar) = 'lw flux, upward, top of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwuptc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwdntc'
      vardesc(nvar) = 'lw flux, downward, top of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwdntc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwupbc'
      vardesc(nvar) = 'lw flux, upward, bottom of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwupbc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'lwdnbc'
      vardesc(nvar) = 'lw flux, downward, bottom of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,lwdnbc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swuptc'
      vardesc(nvar) = 'sw flux, upward, top of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swuptc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swdntc'
      vardesc(nvar) = 'sw flux, downward, top of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swdntc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swupbc'
      vardesc(nvar) = 'sw flux, upward, bottom of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swupbc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( radopt.eq.2 )THEN
      nvar = nvar+1
      varname(nvar) = 'swdnbc'
      vardesc(nvar) = 'sw flux, downward, bottom of atmosphere - clear sky'
      varunit(nvar) = 'W/m^2'
      varlvls(nvar) = 0

      call writehf(1,swdnbc,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    doit = .false.
    IF( doit )THEN
      nvar = nvar+1
      varname(nvar) = 'wa500'
      vardesc(nvar) = 'vertical velocity at 500 mb'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = 0

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,pint,plast)
          do j=1,nj
          do i=1,ni
            pint = 1.0e30
            k = 1
            do while( pint.gt.50000.0 .and. k.lt.nk )
              plast = pint
              k = k + 1
              pint = 0.5*(prs(i,j,k-1)+prs(i,j,k))
            enddo
            dum1(i,j,1) = w3d(i,j,k-1)+(w3d(i,j,k)-w3d(i,j,k-1))  &
                                      *(50000.0-plast)  &
                                      /(pint-plast)
          enddo
          enddo

      call writehf(1,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    domse = .false.

    IF( domse .and. imoist.eq.1 )THEN

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        dum2(i,j,1) = 0.0
      enddo
      enddo

      ! dum8 stores qice:
      if( iice.eq.1 )then
        dum8 = 0.0
      else
        do k=1,nk
          do j=1,nj
          do i=1,ni
            dum8(i,j,k) = 0.0
          enddo
          enddo
          do n=nqs1,nqs2
          do j=1,nj
          do i=1,ni
            dum8(i,j,k) = dum8(i,j,k)+q3d(i,j,k,n)
          enddo
          enddo
          enddo
        enddo
      endif

      do k=1,nk
      do j=1,nj
      do i=1,ni
        tx = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
        sadv(i,j,k) = cp*tx + g*zh(i,j,k)            &
                    +(lv1-lv2*tx)*q3d(i,j,k,nqv)     &
                    -((ls1-ls2*tx)-(lv1-lv2*tx))*dum8(i,j,k)
        dum2(i,j,1) = dum2(i,j,1) + rho(i,j,k)*sadv(i,j,k)*dz*rmh(i,j,k)
      enddo
      enddo
      enddo

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'fmse'
      vardesc(nvar) = 'mass-wgted vert intgrl of frzn moist static energy'
      varunit(nvar) = 'J/m^2'
      varlvls(nvar) = 0

      call writehf(1,dum2(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      sum = 0.0

      do j=1,nj
      do i=1,ni
        sum = sum + dum2(i,j,1)
      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      havg = ( sum/dble(nx*ny) )

      do j=1,nj
      do i=1,ni
        dum2(i,j,2) = dum2(i,j,1)-havg
      enddo
      enddo

      nvar = nvar+1
      varname(nvar) = 'fmsepert'
      vardesc(nvar) = 'perturbation fmse'
      varunit(nvar) = 'J/m^2'
      varlvls(nvar) = 0

      call writehf(1,dum2(ib,jb,2),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      call bcs(sadv)
#ifdef MPI
      call comm_3s_start(sadv,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(  sadv,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
#endif

    if( td_hadv.ge.1 .and. td_vadv.ge.1 .and. ptype.eq.91919 )then

      stop 111

        diffit = 0
        bfoo = 0.0
        weps = cp*10.0*epsilon

!        call advs( 3 ,1,0,bfoo,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,           &
!                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,divx,        &
!                   rru,rrv,rrw,th3d,sadv,thten,0,0,dt   ,weps,                             &
!                   flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,rr,diffit, &
!                   .true.,ibdt,iedt,jbdt,jedt,kbdt,kedt,ntdiag,tdiag,td_hadv,td_vadv,td_lsw, &
!                   td_hidiff,td_vidiff,td_hediff,wprof,dumk1,dumk2,hadvordrs,vadvordrs,kbdy,bndy,out3d)

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        dum3(i,j,1) = 0.0
        dum4(i,j,1) = 0.0
      enddo
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum3(i,j,1) = dum3(i,j,1) + tdiag(i,j,k,td_hadv)*rho(i,j,k)*dz*rmh(i,j,k)
        dum4(i,j,1) = dum4(i,j,1) + tdiag(i,j,k,td_vadv)*rho(i,j,k)*dz*rmh(i,j,k)
      enddo
      enddo
      enddo

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'hadvfmse'
      vardesc(nvar) = 'vert intgrl of horiz advection of fmse'
      varunit(nvar) = 'J/m^2/s'
      varlvls(nvar) = 0

      call writehf(1,dum3(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'vadvfmse'
      vardesc(nvar) = 'vert intgrl of vert advection of fmse'
      varunit(nvar) = 'J/m^2/s'
      varlvls(nvar) = 0

      call writehf(1,dum4(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

    endif

      !c-c-c-c-c-c-c-c-c-c

    ENDIF

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'th'
      vardesc(nvar) = 'potential temperature (K)'
      varunit(nvar) = 'K'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = th0(i,j,k)+th3d(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'prs'
      vardesc(nvar) = 'pressure'
      varunit(nvar) = 'Pa'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = prs(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'uinterp '
      vardesc(nvar) = 'u interpolated to scalar points (grid-relative)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = 0.5*(u3d(i,j,k)+u3d(i+1,j,k))
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'vinterp '
      vardesc(nvar) = 'v interpolated to scalar points (grid-relative)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = 0.5*(v3d(i,j,k)+v3d(i,j+1,k))
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    ! radial and tangential velocity:
    doit = .true.

    docyl:  &
    if( doit )then

      if(myid.eq.0) print *,'  doazimavg,do_adapt_move = ',doazimavg,do_adapt_move

    if( doazimavg .or. do_adapt_move )then

              ! dum1 = u at s pts
              ! dum2 = v at s pts
              ! dum3 = ur at s pts
              ! dum4 = vt at s pts

        do k=1,klev(nlevels)
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.5*(u3d(i,j,k)+u3d(i+1,j,k)) + umove
          dum2(i,j,k) = 0.5*(v3d(i,j,k)+v3d(i,j+1,k)) + vmove
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

      nvar = nvar+1
      varname(nvar) = 'urad'
      vardesc(nvar) = 'radial velocity (grnd-rel)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = dum3(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'vtan'
      vardesc(nvar) = 'tangential velocity (grnd-rel)'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = dum4(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

      ucrit = -2.0

      ! inflow layer depth:
      do j=1,nj
      do i=1,ni
        if( dum3(i,j,1).lt.ucrit )then
          k = 1
          do while( dum3(i,j,k).lt.ucrit .and. k.le.nk )
            k = k+1
          enddo
          dum1(i,j,1) = zh(1,1,k-1)+(zh(1,1,k)-zh(1,1,k-1))  &
                                  *(    ucrit-dum3(i,j,k-1))  &
                                  /(dum3(i,j,k)-dum3(i,j,k-1))
        else
          dum1(i,j,1) = 0.0
        endif
      enddo
      enddo

      nvar = nvar+1
      varname(nvar) = 'zinfl'
      vardesc(nvar) = 'depth of inflow layer (u < -2 m/s)'
      varunit(nvar) = 'm'
      varlvls(nvar) = 0

      call writehf(1,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

    endif

    ENDIF  docyl

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      !c-c-c-c-c-c-c-c-c-c

      nvar = nvar+1
      varname(nvar) = 'winterp '
      vardesc(nvar) = 'w interpolated to scalar points'
      varunit(nvar) = 'm/s'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)

      !c-c-c-c-c-c-c-c-c-c

    IF( imoist.eq.1 .and. nqv.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'qv'
      vardesc(nvar) = 'water vapor mixing ratio'
      varunit(nvar) = 'kg/kg'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = q3d(i,j,k,nqv)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( output_dbz.eq.1 .and. qd_dbz.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'dbz'
      vardesc(nvar) = 'reflectivity'
      varunit(nvar) = 'dbz'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = klev(n)
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = qdiag(i,j,k,qd_dbz)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sgsmodel.eq.1 .or. sgsmodel.eq.4 )THEN
      nvar = nvar+1
      varname(nvar) = 'tke'
      vardesc(nvar) = 'tke'
      varunit(nvar) = 'm2/s2'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = max( 2 , klev(n) )
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = tke3d(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( sgsmodel.ge.1 .or. horizturb.ge.1 )THEN
      nvar = nvar+1
      varname(nvar) = 'kmh'
      vardesc(nvar) = 'kmh'
      varunit(nvar) = 'm2/s'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = max( 2 , klev(n) )
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = kmh(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c

    IF( ipbl.eq.4 .or. ipbl.eq.5 )THEN
      nvar = nvar+1
      varname(nvar) = 'qke'
      vardesc(nvar) = 'twice TKE from MYNN'
      varunit(nvar) = 'm^2/s^2'
      varlvls(nvar) = nlevels

      do n=1,nlevels
        k = max( 1 , klev(n) )
        do j=1,nj
        do i=1,ni
          dum1(i,j,n) = qke(i,j,k)
        enddo
        enddo
      enddo

      call writehf(nlevels,dum1(ib,jb,1),dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    ENDIF

      !c-c-c-c-c-c-c-c-c-c


    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! close output file:

    id0:  IF( myid.eq.0 )THEN

!-------------------------------------------------------------------
!  grads format

    of2:  &
    IF( output_format.eq.1 )THEN

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        ! grads-format
        ! (by default, for now)

        ! close binary output file:
        close(unit=fnum)

        ! write descriptor file:
          string = 'cm1out_hf.ctl'
          open(unit=66,file=string)
          string = 'cm1out_hf_00%y4.dat'
          write(66,101) string
          write(66,102)
          if( outunits.eq.2 )then
            write(66,113) trim(cm1version)
          else
            write(66,103) trim(cm1version)
          endif
          write(66,104) grads_undef
            if(stretch_x.ge.1)then
              write(66,214) nx
              do i=1,nx
                write(66,217) outunitconv*( 0.5*(xfref(i)+xfref(i+1)) )
              enddo
            else
              write(66,204) nx,outunitconv*xh(1),outunitconv*dx
            endif
            if(stretch_y.ge.1)then
              write(66,215) ny
              do j=1,ny
                write(66,217) outunitconv*( 0.5*(yfref(j)+yfref(j+1)) )
              enddo
            else
              write(66,205) ny,outunitconv*yh(1),outunitconv*dy
            endif
            write(66,107) nlevels
            do k=1,nlevels
              write(66,217) outunitconv*sigma(klev(k))
            enddo
          write(66,108) nwriteh
          write(66,109) nvar
          do n=1,nvar
              a1 = varname(n)
              a2 = vardesc(n)
              !---
              a16 = '                '
              nn = len(trim(varunit(n)))
              write(a16(2:15),314) varunit(n)
              write(a16(1:1),301 )       '('
              write(a16(nn+2:nn+2),301 ) ')'
              !---
                write(66,110) a1(1:12),varlvls(n),a2(1:40),a16
          enddo
          write(66,111)
          close(unit=66)

 301    format(a1)
 314    format(a14)

 101    format('dset ^',a)
 102    format('options template')
 103    format('title CM1 high-frequency output, using version ',a,'; units of x,y,z are km; time is generic, see variable mtime for actual times')
 113    format('title CM1 high-frequency output, using version ',a,'; units of x,y,z are meters; time is generic, see variable mtime for actual times')
 104    format('undef ',f10.1)
! 105    format('xdef ',i6,' linear ',f13.6,1x,f13.6)
 105    format('xdef ',i6,' linear ',es14.6e2,1x,es14.6e2)
 106    format('ydef ',i6,' linear 0 1')
 107    format('zdef ',i6,' levels')
! 217    format(2x,f13.6)
 217    format(2x,es14.6e2)
 108    format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
 109    format('vars ',i6)
 110    format(a12,2x,i6,' 99 ',a40,1x,a16)
 111    format('endvars')

! 204    format('xdef ',i6,' linear ',f13.6,1x,f13.6)
! 205    format('ydef ',i6,' linear ',f13.6,1x,f13.6)
 204    format('xdef ',i6,' linear ',es14.6e2,1x,es14.6e2)
 205    format('ydef ',i6,' linear ',es14.6e2,1x,es14.6e2)

 214    format('xdef ',i6,' levels ')
 215    format('ydef ',i6,' levels ')

!-------------------------------------------------------------------
! netcdf format

#ifdef NETCDF
    ELSEIF( output_format.eq.2 )THEN  of2

      call       writehifrq_nc(wloop,nlevels,nvar,varmax,varname,vardesc,varunit,varlvls,ncid,xhref,yhref,sigma,klev,zlev,mtime,dat2(1,1),dat2(1,2),dum1(ib,jb,kb))

#endif
    ENDIF  of2

    ENDIF  id0

  ENDDO  bigwloop

      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  ... Leaving writeout_hifrq '
      if( myid.eq.0 ) print *

    end subroutine writeout_hifrq


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine writehf(nout,var,dat1,dat2,dat3,reqt,orec,nvar,varname,vardesc,varunit,myi1p,myi2p,myj1p,myj2p)
    use input
#ifdef NETCDF
    use writeout_nc_module, only : disp_err
    use netcdf
#endif
    implicit none

    integer, intent(in) :: nout
    real, intent(in), dimension(ib:ie,jb:je,nout) :: var
    real, intent(inout), dimension(d3is,d3js) :: dat1
    real, intent(inout), dimension(d2is,d2js) :: dat2
    real, intent(inout), dimension(d3is,d3js,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(inout) :: orec
    integer, intent(in) :: nvar
    character(len=80), intent(in), dimension(varmax) :: varname,vardesc,varunit
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

    integer :: i,j,k,varid,status

  if( output_format.eq.1 .or. ( output_format.eq.2 .and. wloop.eq.2 ) )then

    if(myid.eq.0) print *,nvar,trim(varname(nvar))

    DO k = 1 , nout

      call gather(var(ib,jb,k),dat1(1,1),dat2(1,1),dat3(1,1,0),reqt,myi1p,myi2p,myj1p,myj2p)

      if( myid.eq.0 )then

        if( output_format.eq.1 )then
          ! grads-format file:
          write(fnum,rec=orec) ((dat2(i,j),i=1,nx),j=1,ny)
          orec=orec+1
#ifdef NETCDF
        elseif( output_format.eq.2 )then
          status = nf90_inq_varid(ncid,trim(varname(nvar)),varid)
          if(status.ne.nf90_noerr)then
            print *,'  Error in write3d, varname = ',varname(nvar)
            print *,nf90_strerror(status)
            call stopcm1
          endif
          call disp_err( nf90_put_var(ncid,varid,dat2,(/1,1,k,1/),(/nx,ny,1,1/)) , .true. )
#endif
        endif

      endif

    ENDDO

  endif

    end subroutine writehf


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine gather(var,dat1,dat2,dat3,reqt,myi1p,myi2p,myj1p,myj2p)
    use input
#ifdef MPI
    use mpi
#endif
    implicit none

    real, intent(in), dimension(ib:ie,jb:je) :: var
    real, intent(inout), dimension(d3is,d3js) :: dat1
    real, intent(inout), dimension(d2is,d2js) :: dat2
    real, intent(inout), dimension(d3is,d3js,0:d3n-1) :: dat3
    integer, intent(inout), dimension(d3t) :: reqt
    integer, intent(in), dimension(numprocs) :: myi1p,myi2p,myj1p,myj2p

    !----

    integer :: i,j
#ifdef MPI
    integer :: fooi,fooj,nn,nnn,n1,n2,ntot,index2,nitmp,njtmp
    integer :: index,tag,proc,reqs
    integer, dimension(mpi_status_size,ppnode-1) :: status1
#endif

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) 
    !-------------------------------------------------------------------

#ifdef MPI
    tag = 1111
#endif

#ifndef MPI
      !-------------------- non-MPI section --------------------!
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        dat2(i,j)=var(i,j)
      enddo
      enddo
#else
      !-------------------- MPI section --------------------!

      IF(myid.ne.nodeleader)THEN
        ! ordinary processor ... send data to nodeleader:
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dat1(i,j)=var(i,j)
        enddo
        enddo
        call MPI_ISEND(dat1(1,1),d3is*d3js,MPI_REAL,nodeleader,tag,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
        ! DONE, ordinary processors
      ELSE
        ! begin nodeleader section:
          ! start receives from all other processors on a node:
          do proc=myid+1,myid+(ppnode-1)
            call MPI_IRECV(dat3(1,1,proc),d3is*d3js,MPI_REAL,proc,tag,MPI_COMM_WORLD,reqt(proc-myid),ierr)
          enddo
        IF(myid.ne.0)THEN
          ! nodeleader, not proc 0:
!$omp parallel do default(shared)  &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            dat3(i,j,myid)=var(i,j)
          enddo
          enddo
          ! wait for receives to finish:
          call mpi_waitall(ppnode-1,reqt(1:ppnode-1),status1,ierr)
          ! send data to processor 0:
          call MPI_ISEND(dat3(1,1,myid),d3is*d3js*ppnode,MPI_REAL,0,tag+1,MPI_COMM_WORLD,reqs,ierr)
          ! wait for send to finish:
          call MPI_WAIT(reqs,MPI_STATUS_IGNORE,ierr)
          ! DONE, nodeleader (not proc 0)
        ELSE
          ! proc 0:
            ! start receives from other nodeleaders:
            do nn = 1,(nodes-1)
              if( nn.le.mynode )then
                proc = (nn-1)*ppnode
              else
                proc = nn*ppnode
              endif
              call MPI_IRECV(dat3(1,1,proc),d3is*d3js*ppnode,MPI_REAL,proc,tag+1,MPI_COMM_WORLD,reqt(ppnode-1+nn),ierr)
            enddo
          ! my data:
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,nj
            do i=1,ni
              dat2(i,j)=var(i,j)
            enddo
            enddo
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
          ! DONE, proc 0
        ENDIF
        !-------------------- end MPI section --------------------!
      ENDIF
      ! seems to help with memory:
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif


      end subroutine gather


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  END MODULE hifrq_module
