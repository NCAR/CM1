  MODULE statpack_module

  implicit none

  private
  public :: statpack,setup_stat_vars

  CONTAINS

      subroutine statpack(mtime,nrec,ndt,dt,dtlast,rtime,adt,acfl,cloudvar,   &
                          qname,budname,qbudget,asq,bsq,                      &
                          name_stat,desc_stat,unit_stat,                      &
                          xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,zf,mf,   &
                          zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,                 &
                          rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                  &
                          dum1,dum2,dum3,dum4,dum5,rho  ,prs,                 &
                          ua,va,wa,ppi,tha,qa,vq  ,kmh,kmv,khh,khv,tkea,qke,  &
                          tke_myj,xkzh,xkzq,xkzm,                             &
                          pta,u10,v10,hpbl,prate,reset,nstatout,restarted)

      use input
      use constants
      use maxmin_module
      use misclibs
      use cm1libs , only : rslf,rsif
#ifdef NETCDF
      use writeout_nc_module, only : writestat_nc
#endif
#ifdef MPI
      use mpi
#endif
      implicit none

      double precision, intent(in) :: mtime
      integer, intent(inout) :: nrec
      integer :: ndt
      real :: dt,dtlast,rtime
      double precision :: adt,acfl
      logical, dimension(maxq) :: cloudvar
      character(len=3), dimension(maxq) :: qname
      character(len=6), dimension(maxq) :: budname
      double precision, dimension(nbudget) :: qbudget
      double precision, dimension(numq) :: asq,bsq
      character(len=40), intent(in), dimension(maxvars) :: name_stat,desc_stat,unit_stat
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: rgzu,rgzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, dimension(stat_out) :: rstat
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,thv0,th0,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,rho,prs
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,va
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,vq
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: qke
      real, intent(in), dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj
      real, intent(in), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: xkzh,xkzq,xkzm
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(in), dimension(ibl:iel,jbl:jel) :: u10,v10,hpbl
      real, intent(in), dimension(ib:ie,jb:je) :: prate
      logical, intent(inout) :: reset
      integer, intent(in) :: nstatout
      logical, intent(in) :: restarted

!-----------------------------------------------------------------------

      integer :: i,j,k,n,nstat,nloop,nkm,kmin,kmax
      character(len=6) :: text1,text2
      real :: qvs,zlev,r1,r2
#ifdef MPI
      double precision, dimension(nbudget) :: cfoo
      double precision, dimension(numq) :: afoo,bfoo
#endif

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

  dostats:  &
  IF( stat_out.gt.0 )THEN

      nstat = 1

      rstat(nstat) = mtime

    IF( adapt_dt.eq.1 )THEN
      nstat = nstat+1
      rstat(nstat) = sngl(  adt/float(max(1,ndt)) )
      if( .not. restarted )then
        acfl         = sngl( acfl/float(max(1,ndt)) )
        reset = .true.
      endif
    ENDIF

      if(stat_w.eq.1) call maxmin(ni,nj,nk+1,wa,nstat,rstat,kmin,kmax,'WMAX  ','WMIN  ')

      if( stat_w.eq.1 .and. (.not. terrain_flag) )then
        nstat = nstat + 1
        rstat(nstat) = zf(1,1,kmax)
        nstat = nstat + 1
        rstat(nstat) = zf(1,1,kmin)
      endif

      wlevs:  &
      if( stat_wlevs.eq.1 .and. maxz.ge.10000.0 )then
        wloop:  do nloop = 1 , 5
          if(     nloop.eq.1 )then
            text1 = 'WMX0.5'
            text2 = 'WMN0.5'
            zlev = 500.0
          elseif( nloop.eq.2 )then
            text1 = 'WMX1  '
            text2 = 'WMN1  '
            zlev = 1000.0
          elseif( nloop.eq.3 )then
            text1 = 'WMX2.5'
            text2 = 'WMN2.5'
            zlev = 2500.0
          elseif( nloop.eq.4 )then
            text1 = 'WMX5  '
            text2 = 'WMN5  '
            zlev = 5000.0
          elseif( nloop.eq.5 )then
            text1 = 'WMX10 '
            text2 = 'WMN10 '
            zlev = 10000.0
          endif
          if( zlev.lt.sigmaf(nk+1) )then
            nkm = nk+1
            IF(.not.terrain_flag)THEN
              ! without terrain:
              do while( sigmaf(nkm).gt.zlev .and. nkm.gt.1 )
                nkm = nkm-1
              enddo
              ! dum1(i,j,1) = wa(i,j,nkm)+(wa(i,j,nkm+1)-wa(i,j,nkm))  &
              !                          *(         zlev-sigmaf(nkm))  &
              !                          /(sigmaf(nkm+1)-sigmaf(nkm))
              r2 = (zlev-sigmaf(nkm))/(sigmaf(nkm+1)-sigmaf(nkm))
              r1 = 1.0-r2
              do j=1,nj
              do i=1,ni
                dum1(i,j,1) = r1*wa(i,j,nkm)+r2*wa(i,j,nkm+1)
              enddo
              enddo
            ELSE
              ! with terrain:
              do j=1,nj
              do i=1,ni
                nkm = nk+1
                do while( zf(i,j,nkm)-zs(i,j).gt.zlev .and. nkm.ge.1 )
                  nkm = nkm-1
                enddo
                r2 = (zlev-(zf(i,j,nkm)-zs(i,j)))/(zf(i,j,nkm+1)-zf(i,j,nkm))
                r1 = 1.0-r2
                dum1(i,j,1) = r1*wa(i,j,nkm)+r2*wa(i,j,nkm+1)
              enddo
              enddo
            ENDIF
          else
            do j=1,nj
            do i=1,ni
              dum1(i,j,1) = 0.0
            enddo
            enddo
          endif
          call maxmin2d(ni,nj,dum1(ib,jb,1),nstat,rstat,text1,text2)
        enddo  wloop
      endif  wlevs

      if(stat_u.eq.1)then
        call maxmin(ni+1,nj,nk,ua,nstat,rstat,kmin,kmax,'UMAX  ','UMIN  ')
        call maxmin2d(ni+1,nj,ua(ib,jb,1),nstat,rstat,'SUMAX ','SUMIN ')
      endif
      if(stat_v.eq.1)then
        call maxmin(ni,nj+1,nk,va,nstat,rstat,kmin,kmax,'VMAX  ','VMIN  ')
!!!      if(myid.eq.0) print *,'  umax:',rstat(nstat)+rstat(nstat-1),rstat(nstat-4)+rstat(nstat-5),rstat(nstat-1)-rstat(nstat-5)
        call maxmin2d(ni,nj+1,va(ib,jb,1),nstat,rstat,'SVMAX ','SVMIN ')
      endif
      if(stat_rmw.eq.1)then
        call getrmw(nstat,rstat,xh,zh,ua,va)
      endif
 
      if(stat_pipert.eq.1) call maxmin(ni,nj,nk,ppi,nstat,rstat,kmin,kmax,'PPIMAX','PPIMIN')

      if(stat_prspert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k)=prs(i,j,k)-p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo
        call maxmin(ni,nj,nk,dum2,nstat,rstat,kmin,kmax,'PPMAX ','PPMIN ')
      endif

      if(stat_thpert.eq.1)then
        call maxmin(ni,nj,nk,tha,nstat,rstat,kmin,kmax,'THPMAX','THPMIN')
        call maxmin2d(ni,nj,tha(ib,jb,1),nstat,rstat,'STHPMX','STHPMN')
      endif

      if(imoist.eq.1.and.stat_q.eq.1)then
        do n=1,numq
          text1='MAX   '
          text2='MIN   '
          write(text1(4:6),121) qname(n)
          write(text2(4:6),121) qname(n)
121       format(a3)
          call maxmin(ni,nj,nk,qa(ib,jb,kb,n),nstat,rstat,kmin,kmax,text1,text2)
        enddo
        call maxmin2d(ni,nj,prate(ib,jb),nstat,rstat,'PRATMX','PRATMN')
      endif

        if(stat_tke.eq.1)then
          if( idoles .and. sgsmodel.eq.1 )then
            call maxmin(ni,nj,nk+1,tkea,nstat,rstat,kmin,kmax,'TKEMAX','TKEMIN')
          endif
          if( idopbl )then
            if( ipbl.eq.4 .or. ipbl.eq.5 )then
              call maxmin(ni,nj,nk,qke,nstat,rstat,kmin,kmax,'QKEMAX','QKEMIN')
            endif
            if( ipbl.eq.6 )then
              call maxmin(ni,nj,nk+1,tke_myj,nstat,rstat,kmin,kmax,'TKEMAX','TKEMIN')
            endif
          endif
        endif

        if(stat_km.eq.1) call maxmin(ni,nj,nk+1,kmh,nstat,rstat,kmin,kmax,'KMHMAX','KMHMIN')
      if( sgsmodel.ge.1 .or. ipbl.eq.2 )then
        if(stat_km.eq.1) call maxmin(ni,nj,nk+1,kmv,nstat,rstat,kmin,kmax,'KMVMAX','KMVMIN')
      endif
        if(stat_kh.eq.1) call maxmin(ni,nj,nk+1,khh,nstat,rstat,kmin,kmax,'KHHMAX','KHHMIN')
      if( sgsmodel.ge.1 .or. ipbl.eq.2 )then
        if(stat_kh.eq.1) call maxmin(ni,nj,nk+1,khv,nstat,rstat,kmin,kmax,'KHVMAX','KHVMIN')
      endif

    if( ipbl.eq.1 )then
    if( stat_kh.eq.1 .or. stat_km.eq.1 )then
        call maxmin(ni,nj,nk+1,xkzh,nstat,rstat,kmin,kmax,'XKZHMX','XKZHMN')
!!!        call maxmin(ni,nj,nk+1,xkzq,nstat,rstat,kmin,kmax,'XKZQMX','XKZQMN')
        call maxmin(ni,nj,nk+1,xkzm,nstat,rstat,kmin,kmax,'XKZMMX','XKZMMN')
    endif
    endif

    if( ipbl.eq.3 )then
    if( stat_kh.eq.1 .or. stat_km.eq.1 )then
        call maxmin(ni,nj,nk+1,xkzh,nstat,rstat,kmin,kmax,'DKTMAX','DKTMIN')
        call maxmin(ni,nj,nk+1,xkzm,nstat,rstat,kmin,kmax,'DKUMAX','DKUMIN')
    endif
    endif

      if(stat_div.eq.1)then
      IF(axisymm.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum5(i,j,k)=                                                     &
              0.5*( (rho0(i,j,k)+rho0(i+1,j,k))*ua(i+1,j,k)                &
                   -(rho0(i,j,k)+rho0(i-1,j,k))*ua(i  ,j,k) )*rdx*uh(i)    &
             +0.5*( (rho0(i,j,k)+rho0(i,j+1,k))*va(i,j+1,k)                &
                   -(rho0(i,j,k)+rho0(i,j-1,k))*va(i,j  ,k) )*rdy*vh(j)    &
             +0.5*( (rho0(i,j,k)+rho0(i,j,k+1))*wa(i,j,k+1)                &
                   -(rho0(i,j,k)+rho0(i,j,k-1))*wa(i,j,k  ) )*rdz*mh(i,j,k)
        enddo
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum5(i,j,k)=                                                     &
              rho0(1,1,k)*( xf(i+1)*ua(i+1,j,k)                            &
                           -xf(i  )*ua(i  ,j,k) )*rdx*uh(i)*rxh(i)         &
             +0.5*( (rho0(i,j,k)+rho0(i,j,k+1))*wa(i,j,k+1)                &
                   -(rho0(i,j,k)+rho0(i,j,k-1))*wa(i,j,k  ) )*rdz*mh(i,j,k)
        enddo
        enddo
        enddo
      ENDIF
        call maxmin(ni,nj,nk,dum5,nstat,rstat,kmin,kmax,'DIVMAX','DIVMIN')
      endif

      IF(imoist.eq.1)THEN

        if(stat_rh.eq.1 .or. stat_the.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvs)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qvs=rslf( prs(i,j,k) , (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k)) )
            dum2(i,j,k)=qa(i,j,k,nqv)*(1.0+qvs*reps)    &
                       /(qvs*(1.0+qa(i,j,k,nqv)*reps))
          enddo
          enddo
          enddo
        endif

        if(stat_rh.eq.1)then
          call maxmin(ni,nj,nk,dum2,nstat,rstat,kmin,kmax,'RHMAX ','RHMIN ')
        endif

        if(iice.eq.1 .and. stat_rhi.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k,qvs)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qvs=rsif( prs(i,j,k) , (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k)) )
            dum3(i,j,k)=qa(i,j,k,nqv)*(1.0+qvs*reps)    &
                       /(qvs*(1.0+qa(i,j,k,nqv)*reps))
          enddo
          enddo
          enddo
          call maxmin(ni,nj,nk,dum3,nstat,rstat,kmin,kmax,'RHIMAX','RHIMIN')
        endif

      ENDIF

        if(iptra.eq.1)then
          do n=1,npt
            text1='MXPT  '
            text2='MNPT  '
            if( n.le.9 )then
              write(text1(5:5),122) n
              write(text2(5:5),122) n
122           format(i1)
            elseif( n.le.99 )then
              write(text1(5:6),123) n
              write(text2(5:6),123) n
123           format(i2)
            else
              write(text1(4:6),124) n
              write(text2(4:6),124) n
124           format(i3)
            endif
            call maxmin(ni,nj,nk,pta(ib,jb,kb,n),nstat,rstat,kmin,kmax,text1,text2)
          enddo
        endif

      IF(imoist.eq.1)THEN

        if(stat_the.eq.1)then
          call calcthe(zh,pi0,th0,dum4,dum2,prs,ppi,tha,qa)
          call maxmin(ni,nj,nk,dum4,nstat,rstat,kmin,kmax,'THEMAX','THEMIN')
          call maxmin2d(ni,nj,dum4(ib,jb,1),nstat,rstat,'STHEMX','STHEMN')
        endif

        if(stat_cloud.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=0.0
          enddo
          enddo
          enddo
          do n=1,numq
            if(cloudvar(n))then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum1(i,j,k)=dum1(i,j,k)+qa(i,j,k,n)
              enddo
              enddo
              enddo
            endif
          enddo
          call cloud(nstat,rstat,zh,dum1)
        endif
      ENDIF

      if(stat_sfcprs.eq.1)then
        call maxmin2d(ni,nj,prs(ib,jb,1),nstat,rstat,'SFPMAX','SFPMIN')
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = cgs1*prs(i,j,1)+cgs2*prs(i,j,2)+cgs3*prs(i,j,3)
        enddo
        enddo
        call maxmin2d(ni,nj,dum1(ib,jb,1),nstat,rstat,'PSFCMX','PSFCMN')
      endif

      if(stat_wsp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=sqrt( (umove+0.5*(ua(i,j,k)+ua(i+1,j,k)))**2     &
                           +(vmove+0.5*(va(i,j,k)+va(i,j+1,k)))**2 )
        enddo
        enddo
        enddo
        call maxmin(ni,nj,nk,dum1,nstat,rstat,kmin,kmax,'WSPMAX','WSPMIN')
        if( .not. terrain_flag )then
          nstat = nstat + 1
          rstat(nstat) = zh(1,1,kmax)
          nstat = nstat + 1
          rstat(nstat) = zh(1,1,kmin)
        endif
        call maxmin2d(ni,nj,dum1(ib,jb,1),nstat,rstat,'SWSPMX','SWSPMN')
      IF(bbc.eq.3)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1)=sqrt( u10(i,j)**2 + v10(i,j)**2 )
        enddo
        enddo
        call maxmin2d(ni,nj,dum1(ib,jb,1),nstat,rstat,'10MWMX','10MWMN')
      ENDIF
      endif

      IF( bbc.eq.3 .and. sfcmodel.ge.1 )THEN
        call maxmin2d(ni,nj,hpbl(ibl,jbl),nstat,rstat,'HPBLMX','HPBLMN')
      ENDIF

      if(stat_cfl.eq.1) call calccfl(nstat,rstat,dt,acfl,uh,vh,mh,ua,va,wa,1)

      if(stat_cfl.eq.1.and.(sgsmodel.ge.1.or.ipbl.eq.2.or.horizturb.eq.1)) call calcksmax(nstat,rstat,dt,uh,vh,mf,kmh,kmv,khh,khv)

      if(stat_vort.eq.1) call vertvort(nstat,rstat,xh,xf,uf,vf,zh,zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,dum1,dum2,ua,va)

      if(stat_tmass.eq.1) call calcmass(nstat,rstat,ruh,rvh,rmh,rho)

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k)=0.0
        dum2(i,j,k)=0.0
        dum3(i,j,k)=0.0
      enddo
      enddo
      enddo
 
      IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=qa(i,j,k,nqv)
        enddo
        enddo
        enddo

        call getqli(qa,dum2,dum3)

        if(stat_tmois.eq.1)then
          call totmois(nstat,rstat,qbudget(budrain),ruh,rvh,rmh,dum1,dum2,dum3,rho)
        endif

        if(stat_qmass.eq.1)then
          do n=1,numq
            IF( (n.eq.nqv) .or.                                 &
                (n.ge.nql1.and.n.le.nql2) .or.                  &
                (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
              text1='   MAS'
              write(text1(1:3),121) qname(n)
              call totq(nstat,rstat,ruh,rvh,rmh,qa(ib,jb,kb,n),rho,text1)
            ENDIF
          enddo
        endif

      ENDIF

        if(imoist.eq.1)then
          if(ptype.eq.1.or.ptype.eq.2)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=vq(i,j,k,3)
            enddo
            enddo
            enddo
          elseif(ptype.eq.6)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=vq(i,j,k,2)
            enddo
            enddo
            enddo
          else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum4(i,j,k)=0.0
            enddo
            enddo
            enddo
          endif
        else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum4(i,j,k)=0.0
          enddo
          enddo
          enddo
        endif
 
      if(stat_tenerg.eq.1)then
        call calcener(nstat,rstat,ruh,rvh,zh,rmh,pi0,th0,rho,ua,va,wa,ppi,tha,    &
                      dum1,dum2,dum3,dum4)
      endif

      if(stat_mo.eq.1)then
        call calcmoe(nstat,rstat,ruh,rvh,rmh,rho,ua,va,wa,dum1,dum2,dum3,dum4)
      endif

      if(stat_tmf.eq.1) call tmf(nstat,rstat,ruh,rvh,rho,wa)

!----------

      IF(imoist.eq.1 .and. stat_pcn.eq.1)THEN
      if(myid.eq.0)then
100     format(2x,a6,':',1x,e13.6)
        do n=1,nbudget
          write(6,100) budname(n),qbudget(n)
          nstat = nstat + 1
          rstat(nstat) = qbudget(n)
        enddo
      endif
      ENDIF

      IF(imoist.eq.1 .and. stat_qsrc.eq.1)THEN
      if(myid.eq.0)then
        do n=1,numq
          text1='as    '
          write(text1(3:5),121) qname(n)
          write(6,100) text1,asq(n)
          nstat = nstat + 1
          rstat(nstat) = asq(n)
        enddo
        do n=1,numq
          text1='bs    '
          write(text1(3:5),121) qname(n)
          write(6,100) text1,bsq(n)
          nstat = nstat + 1
          rstat(nstat) = bsq(n)
        enddo
      endif
      ENDIF

  IF(myid.eq.0)THEN

    if( nstat.ne.stat_out )then
      print *,'  nstat,stat_out = ',nstat,stat_out
      stop 12998
    endif

!-----------------------------------------------------------------------
!  writeitout:  GrADS format

    IF(output_format.eq.1)THEN

      if( stat_out.gt.0 )then
        ! write GrADS descriptor file:
        call write_statsctl(name_stat,desc_stat,unit_stat,nstatout)
      endif

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string = 'cm1out_stats.dat'

      open(unit=60,file=string,form='unformatted',access='direct',   &
           recl=4*nstat,status='unknown')
#ifdef DP
      write(60,rec=nrec) (sngl(rstat(n)),n=1,nstat)
#else
      write(60,rec=nrec) (rstat(n),n=1,nstat)
#endif
      close(unit=60)

!-----------------------------------------------------------------------
!  writeitout:  netcdf format

#ifdef NETCDF
    ELSEIF(output_format.eq.2)THEN

      call writestat_nc(nrec,rtime,nstat,rstat,qname,budname,name_stat,desc_stat,unit_stat)

#endif

!-----------------------------------------------------------------------

    ENDIF

  ENDIF

  ENDIF  dostats

      end subroutine statpack


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setup_stat_vars(name_stat,desc_stat,unit_stat,  &
                                 qname,qunit,budname)
      use input
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_stat,desc_stat,unit_stat
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=20), intent(in), dimension(maxq) :: qunit
      character(len=6), intent(in), dimension(maxq) :: budname

      integer :: n
      character(len=8) text1
      character(len=30) text2
      character(len=50) fname

!-----------------------------------------------------------------------
!        Define all the variables in a stats output file:

    stat_out = 0

      stat_out = stat_out+1
      name_stat(stat_out) = 'mtime'
      desc_stat(stat_out) = 'model time (seconds since beginning of simulation)'
      unit_stat(stat_out) = 's'

    if( adapt_dt.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'dt'
      desc_stat(stat_out) = 'average timestep dt'
      unit_stat(stat_out) = 's'
    endif

    if(stat_w.eq.1)then
      stat_out = stat_out+1
      name_stat(stat_out) = 'wmax'
      desc_stat(stat_out) = 'max vertical velocity'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmin'
      desc_stat(stat_out) = 'min vertical velocity'
      unit_stat(stat_out) = 'm/s'

    if( .not. terrain_flag )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'zwmax'
      desc_stat(stat_out) = 'level of max vertical velocity'
      unit_stat(stat_out) = 'm AGL'

      stat_out = stat_out+1
      name_stat(stat_out) = 'zwmin'
      desc_stat(stat_out) = 'level of min vertical velocity'
      unit_stat(stat_out) = 'm AGL'
    endif
    endif

    if( stat_wlevs.eq.1 .and. maxz.ge.10000.0 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'wmax500'
      desc_stat(stat_out) = 'max vertical velocity at 500 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmin500'
      desc_stat(stat_out) = 'min vertical velocity at 500 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmax1000'
      desc_stat(stat_out) = 'max vertical velocity at 1000 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmin1000'
      desc_stat(stat_out) = 'min vertical velocity at 1000 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmax2500'
      desc_stat(stat_out) = 'max vertical velocity at 2500 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmin2500'
      desc_stat(stat_out) = 'min vertical velocity at 2500 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmax5000'
      desc_stat(stat_out) = 'max vertical velocity at 5000 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmin5000'
      desc_stat(stat_out) = 'min vertical velocity at 5000 m AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmax10k'
      desc_stat(stat_out) = 'max vertical velocity at 10 km AGL'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wmin10k'
      desc_stat(stat_out) = 'min vertical velocity at 10 km AGL'
      unit_stat(stat_out) = 'm/s'
    endif

    if(stat_u.eq.1)then
      stat_out = stat_out+1
      name_stat(stat_out) = 'umax'
      desc_stat(stat_out) = 'max u velocity (grid-rel)'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'umin'
      desc_stat(stat_out) = 'min u velocity (grid-rel)'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sumax'
      desc_stat(stat_out) = 'max u velocity at lowst mod lev (grid-rel)'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sumin'
      desc_stat(stat_out) = 'min u velocity at lowst mod lev (grid-rel)'
      unit_stat(stat_out) = 'm/s'
    endif

    if(stat_v.eq.1)then
      stat_out = stat_out+1
      name_stat(stat_out) = 'vmax'
      desc_stat(stat_out) = 'max v velocity (grid-rel)'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'vmin'
      desc_stat(stat_out) = 'min v velocity (grid-rel)'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'svmax'
      desc_stat(stat_out) = 'max v velocity at lowst mod lev (grid-rel)'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'svmin'
      desc_stat(stat_out) = 'min v velocity at lowst mod lev (grid-rel)'
      unit_stat(stat_out) = 'm/s'
    endif

    if( stat_rmw.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'rmw'
      desc_stat(stat_out) = 'radius of maximum horizontal wind speed'
      unit_stat(stat_out) = 'm'

      stat_out = stat_out+1
      name_stat(stat_out) = 'zmw'
      desc_stat(stat_out) = 'height of maximum horizontal wind speed'
      unit_stat(stat_out) = 'm'
    endif

    if( stat_pipert.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'ppimax'
      desc_stat(stat_out) = 'max nondimensional pressure perturbation'
      unit_stat(stat_out) = 'nondimensional'

      stat_out = stat_out+1
      name_stat(stat_out) = 'ppimin'
      desc_stat(stat_out) = 'min nondimensional pressure perturbation'
      unit_stat(stat_out) = 'nondimensional'
    endif

    if( stat_prspert.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'ppmax'
      desc_stat(stat_out) = 'max pressure perturbation'
      unit_stat(stat_out) = 'Pa'

      stat_out = stat_out+1
      name_stat(stat_out) = 'ppmin'
      desc_stat(stat_out) = 'min pressure perturbation'
      unit_stat(stat_out) = 'Pa'
    endif

    if( stat_thpert.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'thpmax'
      desc_stat(stat_out) = 'max potential temperature perturbation'
      unit_stat(stat_out) = 'K'

      stat_out = stat_out+1
      name_stat(stat_out) = 'thpmin'
      desc_stat(stat_out) = 'min potential temperature perturbation'
      unit_stat(stat_out) = 'K'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sthpmax'
      desc_stat(stat_out) = 'max pot. temp. pert. at lowest model level'
      unit_stat(stat_out) = 'K'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sthpmin'
      desc_stat(stat_out) = 'min pot. temp. pert. at lowest model level'
      unit_stat(stat_out) = 'K'
    endif

    if( stat_q.eq.1 )then
      do n=1,numq
        text1='max     '
        text2='max                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = qunit(n)

        text1='min     '
        text2='min                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = qunit(n)
      enddo

      stat_out = stat_out+1
      name_stat(stat_out) = 'pratemax'
      desc_stat(stat_out) = 'maximum surface precipitation rate'
      unit_stat(stat_out) = 'kg/m2/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'pratemin'
      desc_stat(stat_out) = 'minimum surface precipitation rate'
      unit_stat(stat_out) = 'kg/m2/s'

    endif

  IF(stat_tke.eq.1)THEN
    if( idoles .and. sgsmodel.eq.1 )then

      stat_out = stat_out+1
      name_stat(stat_out) = 'tkemax'
      desc_stat(stat_out) = 'max subgrid tke'
      unit_stat(stat_out) = 'm^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'tkemin'
      desc_stat(stat_out) = 'min subgrid tke'
      unit_stat(stat_out) = 'm^2/s^2'

    endif
    if( idopbl )then
    if( ipbl.eq.4 .or. ipbl.eq.5 )then

      stat_out = stat_out+1
      name_stat(stat_out) = 'qkemax'
      desc_stat(stat_out) = 'max qke'
      unit_stat(stat_out) = 'm^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'qkemin'
      desc_stat(stat_out) = 'min qke'
      unit_stat(stat_out) = 'm^2/s^2'

    endif
    if( ipbl.eq.6 )then

      stat_out = stat_out+1
      name_stat(stat_out) = 'tkemax'
      desc_stat(stat_out) = 'max tke'
      unit_stat(stat_out) = 'm^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'tkemin'
      desc_stat(stat_out) = 'min tke'
      unit_stat(stat_out) = 'm^2/s^2'

    endif
    endif
  ENDIF

    if( stat_km.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'kmhmax'
      desc_stat(stat_out) = 'max horiz eddy viscosity for momentum'
      unit_stat(stat_out) = 'm^2/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'kmhmin'
      desc_stat(stat_out) = 'min horiz eddy viscosity for momentum'
      unit_stat(stat_out) = 'm^2/s'

    if( sgsmodel.ge.1 .or. ipbl.eq.2 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'kmvmax'
      desc_stat(stat_out) = 'max vert eddy viscosity for momentum'
      unit_stat(stat_out) = 'm^2/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'kmvmin'
      desc_stat(stat_out) = 'min vert eddy viscosity for momentum'
      unit_stat(stat_out) = 'm^2/s'
    endif
    endif

    if( stat_kh.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'khhmax'
      desc_stat(stat_out) = 'max horiz eddy diffusivity for scalars'
      unit_stat(stat_out) = 'm^2/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'khhmin'
      desc_stat(stat_out) = 'min horiz eddy diffusivity for scalars'
      unit_stat(stat_out) = 'm^2/s'

    if( sgsmodel.ge.1 .or. ipbl.eq.2 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'khvmax'
      desc_stat(stat_out) = 'max vert eddy diffusivity for scalars'
      unit_stat(stat_out) = 'm^2/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'khvmin'
      desc_stat(stat_out) = 'min vert eddy diffusivity for scalars'
      unit_stat(stat_out) = 'm^2/s'
    endif
    endif

    if( ipbl.eq.1 )then
    if( stat_kh.eq.1 .or. stat_km.eq.1 )then
      !--------
      stat_out = stat_out+1
      name_stat(stat_out) = 'xkzhmax'
      desc_stat(stat_out) = 'max eddy diffusivity for heat (from YSU)'
      unit_stat(stat_out) = 'm^2/s'
      stat_out = stat_out+1
      name_stat(stat_out) = 'xkzhmin'
      desc_stat(stat_out) = 'min eddy diffusivity for heat (from YSU)'
      unit_stat(stat_out) = 'm^2/s'
      !--------
!      stat_out = stat_out+1
!      name_stat(stat_out) = 'xkzqmax'
!      desc_stat(stat_out) = 'max eddy diffusivity for moisture (from YSU)'
!      unit_stat(stat_out) = 'm^2/s'
!      stat_out = stat_out+1
!      name_stat(stat_out) = 'xkzqmin'
!      desc_stat(stat_out) = 'min eddy diffusivity for moisture (from YSU)'
!      unit_stat(stat_out) = 'm^2/s'
      !--------
      stat_out = stat_out+1
      name_stat(stat_out) = 'xkzmmax'
      desc_stat(stat_out) = 'max eddy viscosity (from YSU)'
      unit_stat(stat_out) = 'm^2/s'
      stat_out = stat_out+1
      name_stat(stat_out) = 'xkzmmin'
      desc_stat(stat_out) = 'min eddy viscosity (from YSU)'
      unit_stat(stat_out) = 'm^2/s'
      !--------
    endif
    endif

    if( ipbl.eq.3 )then
    if( stat_kh.eq.1 .or. stat_km.eq.1 )then
      !--------
      stat_out = stat_out+1
      name_stat(stat_out) = 'dktmax'
      desc_stat(stat_out) = 'max thermal diffusivity (from GFSEDMF)'
      unit_stat(stat_out) = 'm^2/s'
      stat_out = stat_out+1
      name_stat(stat_out) = 'dktmin'
      desc_stat(stat_out) = 'min thermal diffusivity (from GFSEDMF)'
      unit_stat(stat_out) = 'm^2/s'
      !--------
      stat_out = stat_out+1
      name_stat(stat_out) = 'dkumax'
      desc_stat(stat_out) = 'max momentum diffusivity (from GFSEDMF)'
      unit_stat(stat_out) = 'm^2/s'
      stat_out = stat_out+1
      name_stat(stat_out) = 'dkumin'
      desc_stat(stat_out) = 'min momentum diffusivity (from GFSEDMF)'
      unit_stat(stat_out) = 'm^2/s'
      !--------
    endif
    endif

    if( stat_div.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'divmax'
      desc_stat(stat_out) = 'max 3d divergence'
      unit_stat(stat_out) = '1/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'divmin'
      desc_stat(stat_out) = 'min 3d divergence'
      unit_stat(stat_out) = '1/s'
    endif

    if( stat_rh.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'rhmax'
      desc_stat(stat_out) = 'max relative humidity wrt liquid'
      unit_stat(stat_out) = 'nondimensional'

      stat_out = stat_out+1
      name_stat(stat_out) = 'rhmin'
      desc_stat(stat_out) = 'min relative humidity wrt liquid'
      unit_stat(stat_out) = 'nondimensional'
    endif

    if( stat_rhi.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'rhimax'
      desc_stat(stat_out) = 'min relative humidity wrt ice'
      unit_stat(stat_out) = 'nondimensional'

      stat_out = stat_out+1
      name_stat(stat_out) = 'rhimin'
      desc_stat(stat_out) = 'min relative humidity wrt ice'
      unit_stat(stat_out) = 'nondimensional'
    endif

    if( iptra.eq.1 )then
      do n=1,npt
        text1='maxpt   '
        text2='max pt                        '
        if( n.le.9 )then
          write(text1(6:6),157) n
          write(text2(7:7),157) n
        elseif( n.le.99 )then
          write(text1(6:7),257) n
          write(text2(7:8),257) n
        else
          write(text1(6:8),258) n
          write(text2(7:9),258) n
        endif

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = 'kg/kg'

        text1='minpt   '
        text2='min pt                        '
        if( n.le.9 )then
          write(text1(6:6),157) n
          write(text2(7:7),157) n
        elseif( n.le.99 )then
          write(text1(6:7),257) n
          write(text2(7:8),257) n
        else
          write(text1(6:8),258) n
          write(text2(7:9),258) n
        endif
157     format(i1)
257     format(i2)
258     format(i3)

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = 'kg/kg'
      enddo
    endif

    if( stat_the.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'themax'
      desc_stat(stat_out) = 'max theta-e below 10 km'
      unit_stat(stat_out) = 'K'

      stat_out = stat_out+1
      name_stat(stat_out) = 'themin'
      desc_stat(stat_out) = 'min theta-e below 10 km'
      unit_stat(stat_out) = 'K'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sthemax'
      desc_stat(stat_out) = 'max theta-e at lowest model level'
      unit_stat(stat_out) = 'K'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sthemin'
      desc_stat(stat_out) = 'min theta-e at lowest model level'
      unit_stat(stat_out) = 'K'
    endif

    if( stat_cloud.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'qctop'
      desc_stat(stat_out) = 'max cloud top height'
      unit_stat(stat_out) = 'm'

      stat_out = stat_out+1
      name_stat(stat_out) = 'qcbot'
      desc_stat(stat_out) = 'min cloud base height'
      unit_stat(stat_out) = 'm'
    endif

    if( stat_sfcprs.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'sprsmax'
      desc_stat(stat_out) = 'max pressure at lowest model level'
      unit_stat(stat_out) = 'Pa'

      stat_out = stat_out+1
      name_stat(stat_out) = 'sprsmin'
      desc_stat(stat_out) = 'min pressure at lowest model level'
      unit_stat(stat_out) = 'Pa'

      stat_out = stat_out+1
      name_stat(stat_out) = 'psfcmax'
      desc_stat(stat_out) = 'max pressure at surface'
      unit_stat(stat_out) = 'Pa'

      stat_out = stat_out+1
      name_stat(stat_out) = 'psfcmin'
      desc_stat(stat_out) = 'min pressure at surface'
      unit_stat(stat_out) = 'Pa'
    endif

    if( stat_wsp.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'wspmax'
      desc_stat(stat_out) = 'max grid-rel horiz wind speed'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'wspmin'
      desc_stat(stat_out) = 'min grid-rel horiz wind speed'
      unit_stat(stat_out) = 'm/s'

    if( .not. terrain_flag )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'zwspmax'
      desc_stat(stat_out) = 'level of max grid-rel horiz wind speed'
      unit_stat(stat_out) = 'm AGL'

      stat_out = stat_out+1
      name_stat(stat_out) = 'zwspmin'
      desc_stat(stat_out) = 'level of min grid-rel horiz wind speed'
      unit_stat(stat_out) = 'm AGL'
    endif

      stat_out = stat_out+1
      name_stat(stat_out) = 'swspmax'
      desc_stat(stat_out) = 'max grid-rel horiz wind speed at l.m.l.'
      unit_stat(stat_out) = 'm/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'swspmin'
      desc_stat(stat_out) = 'min grid-rel horiz wind speed at l.m.l.'
      unit_stat(stat_out) = 'm/s'

      IF(bbc.eq.3)THEN
      if( imove.ne.1 )then
        stat_out = stat_out+1
        name_stat(stat_out) = 'wsp10max'
        desc_stat(stat_out) = 'max horiz wind speed at 10m AGL'
        unit_stat(stat_out) = 'm/s'

        stat_out = stat_out+1
        name_stat(stat_out) = 'wsp10min'
        desc_stat(stat_out) = 'min horiz wind speed at 10m AGL'
        unit_stat(stat_out) = 'm/s'
      else
        stat_out = stat_out+1
        name_stat(stat_out) = 'wsp10max'
        desc_stat(stat_out) = 'max ground-rel. horiz wind speed at 10m AGL'
        unit_stat(stat_out) = 'm/s'

        stat_out = stat_out+1
        name_stat(stat_out) = 'wsp10min'
        desc_stat(stat_out) = 'min ground-rel. horiz wind speed at 10m AGL'
        unit_stat(stat_out) = 'm/s'
      endif
      ENDIF
    endif

    IF( bbc.eq.3 .and. sfcmodel.ge.1 )THEN
        stat_out = stat_out+1
        name_stat(stat_out) = 'hpblmax'
        desc_stat(stat_out) = 'max diagnosed pbl depth'
        unit_stat(stat_out) = 'm'

        stat_out = stat_out+1
        name_stat(stat_out) = 'hpblmin'
        desc_stat(stat_out) = 'min diagnosed pbl depth'
        unit_stat(stat_out) = 'm'
    ENDIF

    if( stat_cfl.eq.1 )then

      IF( adapt_dt.eq.1 )THEN
        stat_out = stat_out+1
        name_stat(stat_out) = 'cflmax'
        desc_stat(stat_out) = 'max Courant number (average)'
        unit_stat(stat_out) = 'nondimensional'
      ELSE
        stat_out = stat_out+1
        name_stat(stat_out) = 'cflmax'
        desc_stat(stat_out) = 'max Courant number'
        unit_stat(stat_out) = 'nondimensional'
      ENDIF

      IF( sgsmodel.ge.1 .or. ipbl.eq.2 .or. horizturb.eq.1 )THEN
        stat_out = stat_out+1
        name_stat(stat_out) = 'kshmax'
        desc_stat(stat_out) = 'max horiz K stability factor'
        unit_stat(stat_out) = 'nondimensional'

        stat_out = stat_out+1
        name_stat(stat_out) = 'ksvmax'
        desc_stat(stat_out) = 'max vert K stability factor'
        unit_stat(stat_out) = 'nondimensional'
      ENDIF

    endif

    if( stat_vort.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'vortsfc'
      desc_stat(stat_out) = 'max vert vorticity at lowest model level'
      unit_stat(stat_out) = '1/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'vort1km'
      desc_stat(stat_out) = 'max vert vorticity at 1 km AGL'
      unit_stat(stat_out) = '1/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'vort2km'
      desc_stat(stat_out) = 'max vert vorticity at 2 km AGL'
      unit_stat(stat_out) = '1/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'vort3km'
      desc_stat(stat_out) = 'max vert vorticity at 3 km AGL'
      unit_stat(stat_out) = '1/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'vort4km'
      desc_stat(stat_out) = 'max vert vorticity at 4 km AGL'
      unit_stat(stat_out) = '1/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'vort5km'
      desc_stat(stat_out) = 'max vert vorticity at 5 km AGL'
      unit_stat(stat_out) = '1/s'
    endif

    if( stat_tmass.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'tmass'
      desc_stat(stat_out) = 'total mass of dry air'
      unit_stat(stat_out) = 'kg'
    endif

    if( stat_tmois.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'tmois'
      desc_stat(stat_out) = 'total mass of moisture'
      unit_stat(stat_out) = 'kg'
    endif

    if(stat_qmass  .eq.1)then
      do n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          text1='mass    '
          text2='total mass of                 '
          write(text1( 5: 7),156) qname(n)
          write(text2(15:17),156) qname(n)

          stat_out = stat_out+1
          name_stat(stat_out) = text1
          desc_stat(stat_out) = text2
          unit_stat(stat_out) = 'kg'
        ENDIF
      enddo
    endif

    if( stat_tenerg.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'ek'
      desc_stat(stat_out) = 'total kinetic energy'
      unit_stat(stat_out) = 'kg m^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'ei'
      desc_stat(stat_out) = 'total internal energy'
      unit_stat(stat_out) = 'kg m^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'ep'
      desc_stat(stat_out) = 'total potential energy'
      unit_stat(stat_out) = 'kg m^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'le'
      desc_stat(stat_out) = 'total latent energy'
      unit_stat(stat_out) = 'kg m^2/s^2'

      stat_out = stat_out+1
      name_stat(stat_out) = 'et'
      desc_stat(stat_out) = 'total energy'
      unit_stat(stat_out) = 'kg m^2/s^2'
    endif

    if( stat_mo.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'tmu'
      desc_stat(stat_out) = 'total E-W momentum'
      unit_stat(stat_out) = 'kg m/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'tmv'
      desc_stat(stat_out) = 'total N-S momentum'
      unit_stat(stat_out) = 'kg m/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'tmw'
      desc_stat(stat_out) = 'total vertical momentum'
      unit_stat(stat_out) = 'kg m/s'
    endif

    if( stat_tmf.eq.1 )then
      stat_out = stat_out+1
      name_stat(stat_out) = 'tmfu'
      desc_stat(stat_out) = 'total upward dry-air mass flux'
      unit_stat(stat_out) = 'kg m/s'

      stat_out = stat_out+1
      name_stat(stat_out) = 'tmfd'
      desc_stat(stat_out) = 'total downward dry-air mass flux'
      unit_stat(stat_out) = 'kg m/s'
    endif

    if(stat_pcn    .eq.1)then
      do n=1,nbudget
        text1='        '
        text2='                              '
        write(text1(1:6),158) budname(n)
        write(text2(1:6),158) budname(n)
158     format(a6)

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = 'unk '
      enddo
    endif

    if(stat_qsrc   .eq.1)then
      do n=1,numq
        text1='as      '
        text2='artificial source of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = 'kg'
      enddo

      do n=1,numq
        text1='bs      '
        text2='bndry source/sink of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)

        stat_out = stat_out+1
        name_stat(stat_out) = text1
        desc_stat(stat_out) = text2
        unit_stat(stat_out) = 'kg'
      enddo
    endif

156   format(a3)

!-----------------------------------------------------------------------

      end subroutine setup_stat_vars


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_statsctl(name_stat,desc_stat,unit_stat,nstatout)
      use input
      use constants , only : grads_undef
      implicit none

      !---------------------------------------------------------------
      ! This subroutine writes the GrADS descriptor file for stats
      !---------------------------------------------------------------

      character(len=40), intent(in), dimension(maxvars) :: name_stat,desc_stat,unit_stat
      integer, intent(in) :: nstatout

      integer :: i,n,nn
      character(len=16) :: a16

    idcheck:  &
    IF( myid.eq.0 )THEN

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string = 'cm1out_stats.ctl'
      if(dowr) write(outfile,*) string
      open(unit=50,file=string,status='unknown')

      string = 'cm1out_stats.dat'

      write(50,301) string
      write(50,302) trim(cm1version)
      write(50,303) grads_undef
      write(50,304)
      write(50,305)
      write(50,306)
      write(50,307) nstatout
      write(50,308) stat_out

      DO n = 1 , stat_out
        a16 = '                '
        nn = len(trim(unit_stat(n)))
        write(a16(2:15),214) unit_stat(n)
        write(a16(1:1),201 )       '('
        write(a16(nn+2:nn+2),201 ) ')'
        write(50,309) name_stat(n),desc_stat(n),a16
      ENDDO

      write(50,310)

      close(unit=50)

    ENDIF  idcheck

201   format(a1)
214   format(a14)

156   format(a3)
301   format('dset ^',a)
302   format('title CM1 stats output, using version ',a,'; time is generic, see variable mtime for actual times')
303   format('undef ',f10.1)
304   format('xdef 1 linear 0 1')
305   format('ydef 1 linear 0 1')
306   format('zdef 1 linear 0 1')
307   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
308   format('vars ',i6)
309   format(a12,' 1 99 ',a40,1x,a16)
310   format('endvars')

      end subroutine write_statsctl

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE statpack_module
