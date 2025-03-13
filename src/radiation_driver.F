  MODULE radiation_module

  implicit none

  private
  public :: radiation_driver

  CONTAINS

      subroutine radiation_driver(mtime,radtim,dt,rbufsz,xh,yh,xf,yf,zf,rmh,c1,c2,     &
                   swten,lwten,swtenc,lwtenc,cldfra,o30,thraten,                       &
                   radsw,rnflx,radswnet,radlwin,dsr,olr,rad2d,                         &
                   effc,effi,effs,effr,effg,effis,                                     &
                   lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc,                &
                   swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,                &
                   lwcf,swcf,coszr,                                                    &
                   xice,xsnow,xlat,xlong,coszen,swddir,swddni,swddif,hrang,            &
                   cldfra1_flag,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11, &
                   prs0,pi0,th0,prs,ppi,tha,rho,qa,                                    &
                   rth0s,prs0s,rho0s,tsk,albd,glw,gsw,emiss,xland,nstep,               &
                   qc_bl,qi_bl,cldfra_bl)

      use input
      use constants
      use bc_module
      use radtrns3d_module, only : radtrns,setradwrk,zenangl,nrad2d,       &
                                   nrsirbm,nrsirdf,nrsuvbm,nrsuvdf,ncosz,  &
                                   ncosss,nfdirir,nfdifir,nfdirpar,nfdifpar
      use module_ra_rrtmg_lw , only : rrtmg_lwrad
      use module_ra_rrtmg_sw , only : rrtmg_swrad
      implicit none

      double precision, intent(in) :: mtime
      double precision, intent(inout) :: radtim
      real, intent(in) :: dt
      integer, intent(in) :: rbufsz
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh,c1,c2
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,swtenc,lwtenc,cldfra,o30,thraten
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin,dsr,olr
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: effc,effi,effs,effr,effg,effis
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: lwcf,swcf,coszr
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: xice,xsnow,xlat,xlong,coszen,swddir,swddni,swddif,hrang
      integer, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: cldfra1_flag
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs0,pi0,th0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs,ppi,tha,rho
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ib:ie,jb:je) :: rho0s,prs0s,rth0s
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,xland
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: albd,glw,gsw,emiss
      integer, intent(in) :: nstep
      real, intent(in), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: qc_bl,qi_bl,cldfra_bl

      !-----------------------------------------------------------------

      integer :: i,j,k
      real :: saltitude,sazimuth,zen,rtime
      real :: albedo,albedoz,tema,temb,frac_snowcover
      logical :: doirrad,dosorad

      real, dimension(2) :: x1
      real, dimension(2) :: y1

      ! 1d arrays for radiation scheme:
      real, dimension(rbufsz) :: radbuf
      real, dimension(nkr) :: swtmp,lwtmp
      real, dimension(nkr) :: tem1,tem2,tem3,tem4,tem5,   &
                              tem6,tem7,tem8,tem9,tem10,   &
                              tem11,tem12,tem13,tem14,tem15,   &
                              tem16,tem17
      real, dimension(nkr) :: teffc,teffi,teffs,teffr,teffg,teffis
      real, dimension(nkr) :: ptprt,pprt,qv,qc,qr,qi,qs,qh,cvr,   &
                              ptbar,pbar,appi,rhostr,zpp,o31

      !-------------
      !  for rrtmg:

      real, parameter :: dpd = 360.0/365.0   ! degrees per day for earth's
                                             ! orbital position (deg/day)

      real, parameter :: degrad = 3.1415926535897932384626433/180.0   ! conversion factor for
                                             ! degrees to radians (pi/180.) (rad/deg)

      logical :: f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,warm_rain,is_CAMMGMP_used
      integer :: has_reqc,has_reqi,has_reqs,icloud,mp_physics,sf_surface_physics,yr,julday,no_src_types,aer_opt,o3input
      real :: julian,xtime,declin,solcon,gmt,tfoo,rcvm,radt_min,co2out

      real, dimension(:,:,:,:), pointer :: tauaer_sw=>null(), ssaaer_sw=>null(), asyaer_sw=>null()

      real, dimension(ibr:ier,jbr:jer) :: swdownc, swddnic, swddirc

      integer, parameter :: cldovrlp = 2
      integer, parameter :: calc_clean_atm_diag = 0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   begin radiation  ccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      dtrad = max( dtrad , dt )

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'Entering RADIATION_DRIVER '


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !-----     Goddard scheme     ----!

      rad_opt:  &
      IF( radopt.eq.1 )THEN

        ! just to be sure:
        call setradwrk(nir,njr,nkr)

        ! time at beginning of timestep:
        rtime=sngl(mtime)

          i = 1
          j = 1
          rtime=sngl(mtime+dt)
          if(dowr) write(outfile,*) '  Calculating radiation tendency:'
          if(timestats.ge.1) time_rad=time_rad+mytime()
          call bcs(prs)
          CALL zenangl( ni,nj,      zf(1,1,1),    &
                rad2d(1,1,ncosz), rad2d(1,1,ncosss), radsw,              &
                dum1(1,1,1),dum1(1,1,2),dum1(1,1,3),dum1(1,1,4),        &
                dum2(1,1,1),dum2(1,1,2),dum2(1,1,3),dum2(1,1,4),        &
                saltitude,sazimuth,dx,dy,dt,rtime,                     &
                ctrlat,ctrlon,year,month,day,hour,minute,second,jday )
          if(myid.eq.0)then
            print *,'    solar zenith angle  (degrees) = ',   &
                                   acos(rad2d(ni,nj,ncosz))*degdpi
            print *,'    solar azimuth angle (degrees) = ',sazimuth*degdpi
          endif
!-----------------------------------------------------------------------
!
!  Calculate surface albedo which is dependent on solar zenith angle
!  and soil moisture. Set the albedo for different types of solar
!  flux to be same.
!
!    rsirbm   Solar IR surface albedo for beam radiation
!    rsirdf   Solar IR surface albedo for diffuse radiation
!    rsuvbm   Solar UV surface albedo for beam radiation
!    rsuvdf   Solar UV surface albedo for diffuse radiation
!
!-----------------------------------------------------------------------
!

          radbuf = 0.0
          tem1 = 0.0
          tem2 = 0.0
          tem3 = 0.0
          tem4 = 0.0
          tem5 = 0.0
          tem6 = 0.0
          tem7 = 0.0
          tem8 = 0.0
          tem9 = 0.0
          tem10 = 0.0
          tem11 = 0.0
          tem12 = 0.0
          tem13 = 0.0
          tem14 = 0.0
          tem15 = 0.0
          tem16 = 0.0
          tem17 = 0.0
          ptprt = 0.0
          pprt = 0.0
          qv = 0.0
          qc = 0.0
          qr = 0.0
          qi = 0.0
          qs = 0.0
          qh = 0.0
          cvr = 0.0
          ptbar = 0.0
          pbar = 0.0
          appi = 0.0
          rhostr = 0.0
          zpp = 0.0
          o31 = 0.0

!$omp parallel do default(shared)  &
!$omp private(i,j,albedo,albedoz,frac_snowcover,tema)
  DO j=1,nj
    DO i=1,ni

      ! let's just use MM5/WRF value, instead:
      albedo = albd(i,j)

      ! arps code for albedo:
      ! (not sure I trust this.....)

!      albedoz = 0.01 * ( EXP( 0.003286         & ! zenith dependent albedo
!          * SQRT( ( ACOS(rad2d(i,j,ncosz))*rad2deg ) ** 3 ) ) - 1.0 )
!
!      IF ( soilmodel == 0 ) THEN             ! soil type not defined
!!!!        stop 12321
!        tema = 0
!      ELSE
!        tema = qsoil(i,j,1)/wsat(soiltyp(i,j))
!      END IF
!
!      frac_snowcover = MIN(snowdpth(i,j)/snowdepth_crit, 1.0)
!
!      IF ( tema > 0.5 ) THEN
!        albedo = albedoz + (1.-frac_snowcover)*0.14                     &
!                         + frac_snowcover*snow_albedo
!      ELSE
!        albedo = albedoz + (1.-frac_snowcover)*(0.31 - 0.34 * tema)     &
!                         + frac_snowcover*snow_albedo
!      END IF
!        albedo = albedoz

      rad2d(i,j,nrsirbm) = albedo
      rad2d(i,j,nrsirdf) = albedo
      rad2d(i,j,nrsuvbm) = albedo
      rad2d(i,j,nrsuvdf) = albedo

    END DO
  END DO
          ! big OpenMP parallelization loop:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,ptprt,pprt,qv,qc,qr,qi,qs,qh,cvr,appi,o31,        &
!$omp tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,        &
!$omp tem11,tem12,tem13,tem14,tem15,tem16,tem17,radbuf,swtmp,lwtmp,   &
!$omp doirrad,dosorad,zpp,ptbar,pbar,rhostr,x1,y1,  &
!$omp teffc,teffi,teffs,teffr,teffg,teffis)
        do j=1,nj
        do i=1,ni
          swtmp = 0.0
          lwtmp = 0.0
          do k=1,nk+2
            ptprt(k) =  tha(i,j,k-1)
             pprt(k) =  prs(i,j,k-1) - prs0(i,j,k-1)
               qv(k) =   qa(i,j,k-1,nqv)
               qc(k) =   0.0
               if( nqc.ge.1 ) qc(k) = qa(i,j,k-1,nqc)
               qr(k) =   0.0
               if( nqr.ge.1 ) qr(k) = qa(i,j,k-1,nqr)
               qi(k) =   0.0
               if( nqi.ge.1 ) qi(k) = qa(i,j,k-1,nqi)
               qs(k) =   0.0
               if( nqs.ge.1 ) qs(k) = qa(i,j,k-1,nqs)
               qh(k) =   0.0
               if( nqg.ge.1 ) qh(k) = qa(i,j,k-1,nqg)
              cvr(k) = cv+cvv*qv(k)+cpl*(qc(k)+qr(k))+cpi*(qi(k)+qs(k)+qh(k))
             appi(k) =  pi0(i,j,k-1) + ppi(i,j,k-1)
              o31(k) =  o30(i,j,k-1)
            teffc(k) = effc(i,j,k-1)
            teffi(k) = effi(i,j,k-1)
            teffs(k) = effs(i,j,k-1)
            teffr(k) = effr(i,j,k-1)
            teffg(k) = effg(i,j,k-1)
           teffis(k) = effis(i,j,k-1)
          enddo
          ptprt(1) = ptprt(2)
           pprt(1) =  pprt(2)
          ptprt(nk+2) = ptprt(nk+1)
           pprt(nk+2) =  pprt(nk+1)
          x1(1) = xf(i)
          x1(2) = xf(i+1)
          y1(1) = yf(j)
          y1(2) = yf(j+1)
          do k=1,nk+3
            zpp(k) =   zf(i,j,k-1)
          enddo
          do k=1,nk+2
            ptbar(k) =  th0(i,j,k-1)
             pbar(k) = prs0(i,j,k-1)
           rhostr(k) =  rho(i,j,k-1)
          enddo
            ptbar(1) = rth0s(i,j)**(-1)
             pbar(1) = prs0s(i,j)
           rhostr(1) = rho0s(i,j)
            doirrad = .true.
            dosorad = .true.
          CALL radtrns(nir,njr,nkr, rbufsz, 0,myid,dx,dy,            &
                 ib,ie,jb,je,kb,ke,xh,yh,prs0s(i,j),olr(i,j),dsr(i,j),  &  ! MS add olr,dsr
                 ptprt,pprt,qv,qc,qr,qi,qs,qh,cvr,                      &
                 ptbar,pbar,appi,o31,rhostr, tsk(i,j), zpp ,                                 &
                 radsw(i,j),rnflx(i,j),radswnet(i,j),radlwin(i,j), rad2d(i,j,ncosss),            &
                 rad2d(i,j,nrsirbm),rad2d(i,j,nrsirdf),rad2d(i,j,nrsuvbm),                       &
                 rad2d(i,j,nrsuvdf), rad2d(i,j,ncosz),sazimuth,                                  &
                 rad2d(i,j,nfdirir),rad2d(i,j,nfdifir),rad2d(i,j,nfdirpar),rad2d(i,j,nfdifpar),  &
                 tem1, tem2, tem3, tem4, tem5,                &
                 tem6, tem7, tem8, tem9, tem10,               &
                 tem11,tem12,tem13,tem14,tem15,tem16,         &
                 radbuf(1), tem17,swtmp,lwtmp,doirrad,dosorad, &
                 teffc,teffi,teffs,teffr,teffg,teffis,        &
                 cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,ptype,g,cp,eqtset)
          do k=1,nk
            swten(i,j,k) = swtmp(k+1)
            lwten(i,j,k) = lwtmp(k+1)
            cldfra(i,j,k) = tem5(nk+3-k)
          enddo
        enddo
        enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !-----     RRTMG scheme     ----!

      ELSEIF( radopt.eq.2 )THEN  rad_opt

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do j=1,nj
          do k=1,nk
          do i=1,ni
            dum4(i,j,k) = (pi0(i,j,k)+ppi(i,j,k))
            dum3(i,j,k) = (th0(i,j,k)+tha(i,j,k))*dum4(i,j,k)
            dum1(i,j,k) = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
            dum2(i,j,k) = c1(i,j,k)*dum3(i,j,k-1)+c2(i,j,k)*dum3(i,j,k)
            dum7(i,j,k) = 0.0
            dum8(i,j,k) = 0.0
            dum9(i,j,k) = 0.0
            dum10(i,j,k) = 0.0
            dum11(i,j,k) = 0.0
          enddo
          enddo
          enddo

          ! dum1 = p8w
          ! dum2 = t8w
          ! dum3 = t      =  t3d
          ! dum4 = pi3d
          ! dum5 = dz8w
          ! prs           =  p3d

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum5(i,j,k) = dz*rmh(i,j,k)
          enddo
          enddo
          enddo

          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            ! surface:
            dum1(i,j,1) = cgs1*prs(i,j,1)+cgs2*prs(i,j,2)+cgs3*prs(i,j,3)
!!!            dum2(i,j,1) = cgs1*dum3(i,j,1)+cgs2*dum3(i,j,2)+cgs3*dum3(i,j,3)
            dum2(i,j,1) = tsk(i,j)   ! cm1r21.1: use surface temp, instead of near-sfc air temp
            ! top of model:
            dum1(i,j,nk+1)= cgt1*prs(i,j,nk)+cgt2*prs(i,j,nk-1)+cgt3*prs(i,j,nk-2)
            dum2(i,j,nk+1)= cgt1*dum3(i,j,nk)+cgt2*dum3(i,j,nk-1)+cgt3*dum3(i,j,nk-2)
          enddo
          enddo

          f_qr = .false.
          f_qs = .false.
          f_qg = .false.
          f_qi = .false.
          f_qc = .false.

          ! store qr in dum11:
          if( nqr.ge.1 )then
            f_qr = .true.
            !$omp parallel do default(shared)   &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum11(i,j,k) = qa(i,j,k,nqr)
            enddo
            enddo
            enddo
          endif

          ! store qs in dum7:
          if( nqs.ge.1 )then
            f_qs = .true.
            !$omp parallel do default(shared)   &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum7(i,j,k) = qa(i,j,k,nqs)
            enddo
            enddo
            enddo
          endif

          ! store qg in dum8:
          if( nqg.ge.1 )then
            f_qg = .true.
            !$omp parallel do default(shared)   &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum8(i,j,k) = qa(i,j,k,nqg)
            enddo
            enddo
            enddo
          endif

          ! store qi in dum9:
          if( nqi.ge.1 )then
            f_qi = .true.
            !$omp parallel do default(shared)   &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum9(i,j,k) = qa(i,j,k,nqi)
            enddo
            enddo
            enddo
            if( nqi2.ge.1 )then
              !$omp parallel do default(shared)   &
              !$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum9(i,j,k) = dum9(i,j,k)+qa(i,j,k,nqi2)
              enddo
              enddo
              enddo
            endif
            if( nqi3.ge.1 )then
              !$omp parallel do default(shared)   &
              !$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                dum9(i,j,k) = dum9(i,j,k)+qa(i,j,k,nqi3)
              enddo
              enddo
              enddo
            endif
          endif

          ! store qc in dum10:
          if( nqc.ge.1 )then
            f_qc = .true.
            !$omp parallel do default(shared)   &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum10(i,j,k) = qa(i,j,k,nqc)
            enddo
            enddo
            enddo
          endif

          ! all microphysics schemes have qv:
          f_qv = .true.

          warm_rain = .false.
          is_CAMMGMP_used = .false.

          has_reqc = 0
          has_reqi = 0
          has_reqs = 0

          if( ptype.eq.3 )then
            ! thompson scheme
            mp_physics = 8
            has_reqc = 1
            has_reqi = 1
            has_reqs = 1
          endif
          if( ptype.eq.5 )then
            ! morrison scheme = 10
            ! use eff arrays from Morrison scheme:
            mp_physics = 10
            has_reqc = 1
            has_reqi = 1
            has_reqs = 1
          endif
          if( ptype.eq.7 )then
            !  wsm6 = 6
            ! use eff arrays from WSM6 scheme:
            mp_physics = 6
            has_reqc = 1
            has_reqi = 1
            has_reqs = 1
          endif
          if( ptype.eq.26 .or. ptype.eq.27 )then
            ! NSSL scheme
            mp_physics = 18
            has_reqc = 1
            has_reqi = 1
            has_reqs = 1
          endif
          if( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 )then
            ! P3 microphysics
            ! use eff arrays from P3 scheme:
            mp_physics = ptype
            has_reqc = 1
            has_reqi = 1
            has_reqs = 0
            f_qs = .false.
            f_qg = .false.
          endif
          if( ptype.eq.55 )then
            ! ISHMAEL scheme = 55
            ! use eff arrays from Jensen_ISHMAEL scheme:
            mp_physics = 55
            has_reqc = 1
            has_reqi = 1
            has_reqs = 0
            f_qs = .false.
            f_qg = .false.
          endif

          icloud = 1
          o3input = 0

          julian = jday - 1.0 + hour/24.0 + minute/(60.0*24.0) + (second+mtime)/(60.0*60.0*24.0)
          julday = jday
          yr = year
          gmt = hour      ! Greenwich Mean Time Hour of model start (hour)

          xtime = sngl( mtime/60.0d0 )  ! time since simulation start (min)

          xlong = ctrlon
          xlat  = ctrlat

     CALL radconst(XTIME,DECLIN,SOLCON,JULIAN,               &
                   DEGRAD,DPD                                )

!!!    print *,'  xtime,declin,solcon = ',xtime,declin,solcon
!!!    print *,'  julian,degrad,dpd = ',julian,degrad,dpd

    tfoo = ( mtime + 0.5*dtrad )/60.0

    call calc_coszen(ib,ie,jb,je,1,ni,1,nj,  &
                      julian,tfoo,gmt, &
                      declin,degrad,xlong,xlat,coszen,hrang)

        IF( testcase.eq.8 )THEN
          ! Bretherton et al 2005 (fixed values):
          solcon = 650.83
          coszen = cos(50.5*pi/180.0)
!!!          ! rcemip:
!!!          solcon = 551.58
!!!          coszen = cos(42.05*pi/180.0)
!!!          albd = 0.07
        ENDIF

        if( nstep.le.1 .and. myid.eq.0 )then
          print *,'  solcon = ',solcon
          print *,'  coszen = ',coszen(1,1),coszen(ni,nj)
          print *,'  albd   = ',albd(1,1),albd(ni,nj)
        endif

        ! save info:
        rad_solcon = solcon
        rad_jday = julian
        rad_zenangle = acos(coszen(1,1))
        rad_declin = declin
        rad_hrang = hrang(1,1)

!!!    print *,'  xtime,gmt,its,ite = ',xtime+dtrad,gmt,1,ni
!!!    print *,'  xlong  = ',xlong(1,1),xlong(ni/2,1),xlong(ni,1)
!!!    print *,'  xlat   = ',xlat(1,1),xlat(ni/2,1),xlat(ni,1)
!!!    print *,'  coszen = ',coszen(1,1),coszen(ni/2,1),coszen(ni,1)
!!!    print *,'  hrang  = ',hrang(1,1),hrang(ni/2,1),hrang(ni,1)

             !----------------------------------------------------------------

            !$omp parallel do default(shared)  &
            !$omp private(i,j)
             do j=1,nj
             do i=1,ni
               GSW(I,J)=0.0
               GLW(I,J)=0.0
             enddo
             enddo

      ! NOTE: qs is stored in dum7
      !       qg is stored in dum8
      !       qi is stored in dum9
      !       qc is stored in dum10

        CALL cal_cldfra1(CLDFRA=CLDFRA,                                             &
                  QV=qa(ib,jb,kb,nqv),QC=dum10,QI=dum9,QS=dum7,                     &
                   F_QV=F_QV,F_QC=F_QC,F_QI=F_QI,F_QS=F_QS,t_phy=dum3,p_phy=prs,    &
                   mp_physics=mp_physics,cldfra1_flag=cldfra1_flag,                 &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,          &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,            &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk )

  !-----------------------------------------------------------------------------------------

    ! cm1r20.2:  couple mynn subgrid-scale clouds to radiation
    !   (code is from WRFV4.2.1 radiation_driver)

           ! MYNN PBL only:
        if( idopbl )then
        if( ipbl.eq.4 .or. ipbl.eq.5 )then
        if( imoist.eq.1 )then
           if( myid.eq.0 ) print *,'  ... icloud_bl in radiation_driver ... '
           IF ( icloud_bl > 0 ) THEN
              IF (nstep .NE. 1) THEN
                 DO k = 1,nk
                 DO j = 1,nj
                 DO i = 1,ni
                    CLDFRA(i,j,k)=CLDFRA_BL(i,j,k)
                 ENDDO
                 ENDDO
                 ENDDO
              ENDIF

              DO k = 1,nk
              DO j = 1,nj
              DO i = 1,ni
                 IF (dum10(i,j,k) < 1.E-6 .AND. CLDFRA_BL(i,j,k) > 0.001) THEN
                     dum10(i,j,k)=dum10(i,j,k) + QC_BL(i,j,k)*CLDFRA_BL(i,j,k)
                 ENDIF
                 IF (dum9(i,j,k) < 1.E-8 .AND. CLDFRA_BL(i,j,k) > 0.001) THEN
                    dum9(i,j,k)=dum9(i,j,k) + QI_BL(i,j,k)*CLDFRA_BL(i,j,k)
                 ENDIF
              ENDDO
              ENDDO
              ENDDO

           ENDIF
        endif
        endif
        endif

  !-----------------------------------------------------------------------------------------

      ! NOTE: qs is stored in dum7
      !       qg is stored in dum8
      !       qi is stored in dum9
      !       qc is stored in dum10

             if(myid.eq.0) print *,'  rrtmg_lwrad '
             CALL RRTMG_LWRAD(                                      &
                  RTHRATENLW=lwten,lwtenc=lwtenc,                   &
                  LWUPT=LWUPT,LWUPTC=LWUPTC,                        &
                  LWDNT=LWDNT,LWDNTC=LWDNTC,                        &
                  LWUPB=LWUPB,LWUPBC=LWUPBC,                        &
                  LWDNB=LWDNB,LWDNBC=LWDNBC,                        &
                  GLW=GLW,OLR=dum6(ib,jb,1),LWCF=LWCF,              &
                  EMISS=EMISS,                                      &
                  P8W=dum1,P3D=prs,PI3D=dum4,DZ8W=dum5,TSK=tsk,T3D=dum3,    &
                  T8W=dum2,RHO3D=rho,R=rd,G=G,                      &
                  ICLOUD=icloud,WARM_RAIN=warm_rain,                &
                  CLDFRA3D=CLDFRA,                                  &
                  cldovrlp=cldovrlp,                                &
                  IS_CAMMGMP_USED=is_cammgmp_used,                  &
                  XLAND=XLAND,XICE=XICE,SNOW=XSNOW,                 &
                  QV3D=qa(ib,jb,kb,nqv),QC3D=dum10           ,QR3D=dum11           ,     &
                  QI3D=dum9            ,QS3D=dum7            ,QG3D=dum8            ,     &
                  O3INPUT=O3INPUT,O33D=o30,                         &
                  F_QV=F_QV,F_QC=F_QC,F_QR=F_QR,                    &
                  F_QI=F_QI,F_QS=F_QS,F_QG=F_QG,                    &
                  RE_CLOUD=effc,RE_ICE=effi,RE_SNOW=effs,  & ! G. Thompson
                  has_reqc=has_reqc,has_reqi=has_reqi,has_reqs=has_reqs, & ! G. Thompson
                  YR=YR,JULIAN=JULIAN,                              &
                  calc_clean_atm_diag=calc_clean_atm_diag,          &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,          &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,            &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk ,            &
                  mp_physics=mp_physics,co2out=co2out               )


            !$omp parallel do default(shared)  &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              OLR(I,J)=dum6(i,j,1)
            enddo
            enddo

            if( nstep.le.1 .and. myid.eq.0 ) print *,'  co2 = ',co2out


             !-------------------------

             no_src_types = 6

             sf_surface_physics = 1

             aer_opt = 0

             radt_min = dtrad/60.0

      ! NOTE: qs is stored in dum7
      !       qg is stored in dum8
      !       qi is stored in dum9
      !       qc is stored in dum10

             if(myid.eq.0) print *,'  rrtmg_swrad '
             CALL RRTMG_SWRAD(                                         &
                     RTHRATENSW=swten,swtenc=swtenc,                   &
                     SWUPT=SWUPT,SWUPTC=SWUPTC,                        &
                     SWDNT=SWDNT,SWDNTC=SWDNTC,                        &
                     SWUPB=SWUPB,SWUPBC=SWUPBC,                        &
                     SWDNB=SWDNB,SWDNBC=SWDNBC,                        &
                     SWCF=SWCF,GSW=GSW,                                &
                     XTIME=XTIME,GMT=GMT,XLAT=XLAT,XLONG=XLONG,        &
                     RADT=radt_min,DEGRAD=DEGRAD,DECLIN=DECLIN,        &
                     COSZR=COSZR,JULDAY=JULDAY,SOLCON=SOLCON,          &
                     ALBEDO=ALBD,t3d=dum3,t8w=dum2,TSK=TSK,            &
                     p3d=prs,p8w=dum1,pi3d=dum4,rho3d=rho,             &
                     dz8w=dum5,                                        &
                     CLDFRA3D=CLDFRA,                                  &
                     IS_CAMMGMP_USED=is_cammgmp_used,                  &
                     R=rd,G=G,              &
                     ICLOUD=icloud,WARM_RAIN=warm_rain,                &
                     cldovrlp=cldovrlp,                                &
                     XLAND=XLAND,XICE=XICE,SNOW=XSNOW,                 &
                  QV3D=qa(ib,jb,kb,nqv),QC3D=dum10           ,QR3D=dum11           ,     &
                  QI3D=dum9            ,QS3D=dum7            ,QG3D=dum8            ,     &
                  O3INPUT=O3INPUT,O33D=o30,                         &
                     AER_OPT=AER_OPT,               &
                     no_src=no_src_types,  &
                     SF_SURFACE_PHYSICS=sf_surface_physics,            &  !Zhenxin ssib sw_phy   (06/2010)
                     F_QV=f_qv,F_QC=f_qc,F_QR=f_qr,                    &
                     F_QI=f_qi,F_QS=f_qs,F_QG=f_qg,                    &
                     RE_CLOUD=effc,RE_ICE=effi,RE_SNOW=effs,  & ! G. Thompson
                     has_reqc=has_reqc,has_reqi=has_reqi,has_reqs=has_reqs, & ! G. Thompson
                     calc_clean_atm_diag=calc_clean_atm_diag,          &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,          &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,            &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk ,            &
                     tauaer3d_sw=tauaer_sw,                             & ! jararias 2013/11
                     ssaaer3d_sw=ssaaer_sw,                             & ! jararias 2013/11
                     asyaer3d_sw=asyaer_sw,                             & ! jararias 2013/11
                     swddir=swddir,swddni=swddni,swddif=swddif,         & ! jararias 2013/08/10
                     swdownc=swdownc, swddnic=swddnic, swddirc=swddirc, &
                     xcoszen=coszen,yr=yr,julian=julian,mp_physics=mp_physics,co2out=co2out ) ! jararias 2013/08/14

            if( nstep.le.1 .and. myid.eq.0 ) print *,'  co2 = ',co2out

        IF( eqtset.eq.2 )THEN
          ! for Bryan-Fritsch equation set:
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,rcvm)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! convert to cm1 temp tendency:
            rcvm = cp/( cv+cvv*max(0.0,qa(i,j,k,nqv))                               &
                          +cpl*max(0.0,(dum10(i,j,k)+dum11(i,j,k)))                 &
                          +cpi*max(0.0,(dum9(i,j,k)+dum7(i,j,k)+dum8(i,j,k))) )
            lwten(i,j,k) = lwten(i,j,k)*rcvm
            swten(i,j,k) = swten(i,j,k)*rcvm
            lwtenc(i,j,k) = lwtenc(i,j,k)*rcvm
            swtenc(i,j,k) = swtenc(i,j,k)*rcvm
          enddo
          enddo
          enddo
        ENDIF

      ENDIF  rad_opt

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! cm1r21.1: added thraten array

      IF( eqtset.eq.1 )THEN
        ! traditional equation set (pressure is constant):
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ! cm1r17:  swten and lwten now store TEMPERATURE tendencies:
          thraten(i,j,k) = (swten(i,j,k)+lwten(i,j,k))/(pi0(i,j,k)+ppi(i,j,k))
        enddo
        enddo
        enddo
      ELSEIF( eqtset.eq.2 )THEN
        ! Bryan-Fritsch equation set (mass is constant):
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ! cm1r17:  swten and lwten now store TEMPERATURE tendencies:
          thraten(i,j,k) = (swten(i,j,k)+lwten(i,j,k))*cvdcp/(pi0(i,j,k)+ppi(i,j,k))
        enddo
        enddo
        enddo
      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   end radiation  ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(myid.eq.0) print *,'Leaving RADIATION_DRIVER '
      if(myid.eq.0) print *

      end subroutine radiation_driver



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!----------------------------------------------------------------------!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


!---------------------------------------------------------------------
!BOP
! !IROUTINE: radconst - compute radiation terms
! !INTERFAC:
   SUBROUTINE radconst(XTIME,DECLIN,SOLCON,JULIAN,                   &
                       DEGRAD,DPD                                    )
!---------------------------------------------------------------------
   USE module_wrf_error
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   REAL, INTENT(IN   )      ::       DEGRAD,DPD,XTIME,JULIAN
   REAL, INTENT(OUT  )      ::       DECLIN,SOLCON
   REAL                     ::       OBECL,SINOB,SXLONG,ARG,  &
                                     DECDEG,DJUL,RJUL,ECCFAC
!
! !DESCRIPTION:
! Compute terms used in radiation physics 
!EOP

! for short wave radiation

   DECLIN=0.
   SOLCON=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.
        
   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)
        
!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)
   SXLONG=SXLONG*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)
   DECDEG=DECLIN/DEGRAD
!----SOLAR CONSTANT ECCENTRICITY FACTOR (PALTRIDGE AND PLATT 1976)
   DJUL=JULIAN*360./365.
   RJUL=DJUL*DEGRAD
   ECCFAC=1.000110+0.034221*COS(RJUL)+0.001280*SIN(RJUL)+0.000719*  &
          COS(2*RJUL)+0.000077*SIN(2*RJUL)
   SOLCON=1370.*ECCFAC
   
   END SUBROUTINE radconst


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!----------------------------------------------------------------------!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


   SUBROUTINE calc_coszen(ims,ime,jms,jme,its,ite,jts,jte,  &
                          julian,xtime,gmt, &
                          declin,degrad,xlon,xlat,coszen,hrang)
       ! Added Equation of Time correction : jararias, 2013/08/10
       implicit none
       integer, intent(in) :: ims,ime,jms,jme,its,ite,jts,jte
       real, intent(in)    :: julian,declin,xtime,gmt,degrad
       real, dimension(ims:ime,jms:jme), intent(in)    :: xlat,xlon
       real, dimension(ims:ime,jms:jme), intent(inout) :: coszen,hrang

       integer :: i,j
       real    :: da,eot,xt24,tloctm,xxlat

       da=6.2831853071795862*(julian-1)/365.
       eot=(0.000075+0.001868*cos(da)-0.032077*sin(da) &
            -0.014615*cos(2*da)-0.04089*sin(2*da))*(229.18)
       xt24=mod(xtime,1440.)+eot

     !$omp parallel do default(shared)  &
     !$omp private(i,j,tloctm,xxlat)
       do j=jts,jte
          do i=its,ite
             tloctm=gmt+xt24/60.+xlon(i,j)/15.
             hrang(i,j)=15.*(tloctm-12.)*degrad
             xxlat=xlat(i,j)*degrad
             coszen(i,j)=sin(xxlat)*sin(declin) &
                        +cos(xxlat)*cos(declin) *cos(hrang(i,j))
             coszen(i, j) = min (max (coszen(i, j), -1.0), 1.0)
          enddo
       enddo

   END SUBROUTINE calc_coszen


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!----------------------------------------------------------------------!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


! !IROUTINE: cal_cldfra1 - Compute cloud fraction
! !INTERFACE:
! cal_cldfra_xr - Compute cloud fraction.
! Code adapted from that in module_ra_gfdleta.F in WRF_v2.0.3 by James Done
!!
!!---  Cloud fraction parameterization follows Xu and Randall (JAS), 1996
!!     (see Hong et al., 1998)
!!     (modified by Ferrier, Feb '02)
!
   SUBROUTINE cal_cldfra1(CLDFRA, QV, QC, QI, QS,                    &
                         F_QV, F_QC, F_QI, F_QS, t_phy, p_phy,       &
                         F_ICE_PHY,F_RAIN_PHY,                       &
                         mp_physics, cldfra1_flag,                   &
          ids,ide, jds,jde, kds,kde,                                 &
          ims,ime, jms,jme, kms,kme,                                 &
          its,ite, jts,jte, kts,kte                                  )
     USE module_state_description, ONLY : KFCUPSCHEME, KFETASCHEME       !wig, CuP 4-Fb-2008 !BSINGH - For WRFCuP scheme

#if (HWRF == 1)
   USE module_state_description, ONLY : FER_MP_HIRES, FER_MP_HIRES_ADVECT, ETAMP_HWRF 
#else
   USE module_state_description, ONLY : FER_MP_HIRES, FER_MP_HIRES_ADVECT
#endif
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

!
   INTEGER, DIMENSION( ims:ime, jms:jme , kms:kme ), INTENT(OUT  ) :: cldfra1_flag
   REAL, DIMENSION( ims:ime, jms:jme , kms:kme ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, jms:jme , kms:kme ), INTENT(IN   ) ::    &
                                                                 QV, &
                                                                 QI, &
                                                                 QC, &
                                                                 QS, &
                                                              t_phy, &
                                                              p_phy
!                                                              p_phy, &
!                                                          F_ICE_PHY, &
!                                                         F_RAIN_PHY

   REAL, DIMENSION( ims:ime, jms:jme , kms:kme ),                     &
         OPTIONAL,                                                   &
         INTENT(IN   ) ::                                            &
                                                          F_ICE_PHY, &
                                                         F_RAIN_PHY
   LOGICAL,OPTIONAL,INTENT(IN) :: F_QC,F_QI,F_QV,F_QS
   INTEGER :: mp_physics

!  REAL thresh
   INTEGER:: i,j,k
   REAL    :: RHUM, tc, esw, esi, weight, qvsw, qvsi, qvs_weight, QIMID, QWMID, QCLD, DENOM, ARG, SUBSAT

   REAL    ,PARAMETER :: ALPHA0=100., GAMMA=0.49, QCLDMIN=1.E-12,    &
                                        PEXP=0.25, RHGRID=1.0
   REAL    , PARAMETER ::  SVP1=0.6112
   REAL    , PARAMETER ::  SVP2=17.67
   REAL    , PARAMETER ::  SVPI2=21.8745584
   REAL    , PARAMETER ::  SVP3=29.65
   REAL    , PARAMETER ::  SVPI3=7.66
   REAL    , PARAMETER ::  SVPT0=273.15
   REAL    , PARAMETER ::  r_d = 287.04
   REAL    , PARAMETER ::  r_v = 461.5
   REAL    , PARAMETER ::  ep_2=r_d/r_v
! !DESCRIPTION:
! Compute cloud fraction from input ice and cloud water fields
! if provided.
!
! Whether QI or QC is active or not is determined from the indices of
! the fields into the 4D scalar arrays in WRF. These indices are 
! P_QI and P_QC, respectively, and they are passed in to the routine
! to enable testing to see if QI and QC represent active fields in
! the moisture 4D scalar array carried by WRF.
! 
! If a field is active its index will have a value greater than or
! equal to PARAM_FIRST_SCALAR, which is also an input argument to 
! this routine.
!EOP


!-----------------------------------------------------------------------
!---  COMPUTE GRID-SCALE CLOUD COVER FOR RADIATION
!     (modified by Ferrier, Feb '02)
!
!---  Cloud fraction parameterization follows Randall, 1994
!     (see Hong et al., 1998)
!-----------------------------------------------------------------------
! Note: ep_2=287./461.6 Rd/Rv
! Note: R_D=287.

! Alternative calculation for critical RH for grid saturation
!     RHGRID=0.90+.08*((100.-DX)/95.)**.5

! Calculate saturation mixing ratio weighted according to the fractions of
! water and ice.
! Following:
! Murray, F.W. 1966. ``On the computation of Saturation Vapor Pressure''  J. Appl. Meteor.  6 p.204
!    es (in mb) = 6.1078 . exp[ a . (T-273.16)/ (T-b) ]
!
!       over ice        over water
! a =   21.8745584      17.2693882
! b =   7.66            35.86

!---------------------------------------------------------------------

  !$omp parallel do default(shared)  &
  !$omp private(i,j,k,tc,esw,esi,qvsw,qvsi,qcld,weight,qimid,qwmid,qvs_weight,rhum,subsat,denom,arg)
    DO k = kts,kte
    DO j = jts,jte
    DO i = its,ite
      tc         = t_phy(i,j,k) - SVPT0
      esw     = 1000.0 * SVP1 * EXP( SVP2  * tc / ( t_phy(i,j,k) - SVP3  ) )
      esi     = 1000.0 * SVP1 * EXP( SVPI2 * tc / ( t_phy(i,j,k) - SVPI3 ) )
      QVSW = EP_2 * esw / ( p_phy(i,j,k) - esw )
      QVSI = EP_2 * esi / ( p_phy(i,j,k) - esi )

      ifouter: IF ( PRESENT(F_QI) .and. PRESENT(F_QC) .and. PRESENT(F_QS) ) THEN

! mji - For MP options 2, 4, 6, 7, 8, etc. (qc = liquid, qi = ice, qs = snow)
         IF ( F_QI .and. F_QC .and. F_QS) THEN
            QCLD = QI(i,j,k)+QC(i,j,k)+QS(i,j,k)
            IF (QCLD .LT. QCLDMIN) THEN
               weight = 0.
            ELSE
               weight = (QI(i,j,k)+QS(i,j,k)) / QCLD
            ENDIF
         ENDIF

! for P3, mp option 50 or 51
         IF ( F_QI .and. F_QC .and. .not. F_QS) THEN
            QCLD = QI(i,j,k)+QC(i,j,k)
            IF (QCLD .LT. QCLDMIN) THEN
               weight = 0.
            ELSE
               weight = (QI(i,j,k)) / QCLD
            ENDIF
         ENDIF

! mji - For MP options 1 and 3, (qc only)
!  For MP=1, qc = liquid, for MP=3, qc = liquid or ice depending on temperature
         IF ( F_QC .and. .not. F_QI .and. .not. F_QS ) THEN
            QCLD = QC(i,j,k)
            IF (QCLD .LT. QCLDMIN) THEN
               weight = 0.
            ELSE
               if (t_phy(i,j,k) .gt. 273.15) weight = 0.
               if (t_phy(i,j,k) .le. 273.15) weight = 1.
            ENDIF
         ENDIF

! mji - For MP option 5; (qc = liquid, qs = ice)
         IF ( F_QC .and. .not. F_QI .and. F_QS .and. PRESENT(F_ICE_PHY) ) THEN

! Mixing ratios of cloud water & total ice (cloud ice + snow).
! Mixing ratios of rain are not considered in this scheme.
! F_ICE is fraction of ice
! F_RAIN is fraction of rain

           QIMID = QS(i,j,k)
           QWMID = QC(i,j,k)
! old method
!           QIMID = QC(i,j,k)*F_ICE_PHY(i,j,k)
!           QWMID = (QC(i,j,k)-QIMID)*(1.-F_RAIN_PHY(i,j,k))
!
!--- Total "cloud" mixing ratio, QCLD.  Rain is not part of cloud,
!    only cloud water + cloud ice + snow
!
           QCLD=QWMID+QIMID
           IF (QCLD .LT. QCLDMIN) THEN
              weight = 0.
           ELSE
              weight = F_ICE_PHY(i,j,k)
           ENDIF
         ENDIF
!BSF - For HWRF MP option; (qc = liquid, qi = cloud ice+snow)
!         IF ( F_QC .and. F_QI .and. .not. F_QS ) THEN
#if (HWRF == 1)
         IF ( mp_physics .eq. FER_MP_HIRES .or. &
	      mp_physics .eq. FER_MP_HIRES_ADVECT .or. &
	      mp_physics .eq. ETAMP_HWRF) THEN
#else
         IF ( mp_physics .eq. FER_MP_HIRES .or. &
              mp_physics==fer_mp_hires_advect) THEN
#endif
           QIMID = QI(i,j,k)     !- total ice (cloud ice + snow)
           QWMID = QC(i,j,k)     !- cloud water
           QCLD=QWMID+QIMID      !- cloud water + total ice
           IF (QCLD .LT. QCLDMIN) THEN
              weight = 0.
           ELSE
              weight = QIMID/QCLD
              if (tc<-40.) weight=1.
           ENDIF
         ENDIF

      ELSE
         CLDFRA(i,j,k)=0.

      ENDIF ifouter !  IF ( F_QI .and. F_QC .and. F_QS)


      QVS_WEIGHT = (1-weight)*QVSW + weight*QVSI
      RHUM=QV(i,j,k)/QVS_WEIGHT   !--- Relative humidity
!
!--- Determine cloud fraction (modified from original algorithm)
!
      cldfra1_flag(i,j,k) = 0
      IF (QCLD .LT. QCLDMIN) THEN
!
!--- Assume zero cloud fraction if there is no cloud mixing ratio
!
        CLDFRA(i,j,k)=0.
        cldfra1_flag(i,j,k) = 1
      ELSEIF(RHUM.GE.RHGRID)THEN
!
!--- Assume cloud fraction of unity if near saturation and the cloud
!    mixing ratio is at or above the minimum threshold
!
        CLDFRA(i,j,k)=1.
        cldfra1_flag(i,j,k) = 2
      ELSE
         cldfra1_flag(i,j,k) = 3
!
!--- Adaptation of original algorithm (Randall, 1994; Zhao, 1995)
!    modified based on assumed grid-scale saturation at RH=RHgrid.
!
        SUBSAT=MAX(1.E-10,RHGRID*QVS_WEIGHT-QV(i,j,k))
        DENOM=(SUBSAT)**GAMMA
        ARG=MAX(-6.9, -ALPHA0*QCLD/DENOM)    ! <-- EXP(-6.9)=.001
! prevent negative values  (new)
        RHUM=MAX(1.E-10, RHUM)
        CLDFRA(i,j,k)=(RHUM/RHGRID)**PEXP*(1.-EXP(ARG))
!!              ARG=-1000*QCLD/(RHUM-RHGRID)
!!              ARG=MAX(ARG, ARGMIN)
!!              CLDFRA(i,j,k)=(RHUM/RHGRID)*(1.-EXP(ARG))
        IF (CLDFRA(i,j,k) .LT. .01) CLDFRA(i,j,k)=0.
           
     ENDIF          !--- End IF (QCLD .LT. QCLDMIN) ...     
    ENDDO          !--- End DO i
    ENDDO          !--- End DO k
    ENDDO          !--- End DO j

   END SUBROUTINE cal_cldfra1


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!----------------------------------------------------------------------!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


  END MODULE radiation_module
