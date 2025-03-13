  MODULE turb_module

  implicit none

  private
  public :: sfc_and_turb,getepst,getepsd,tkebc,calcnm,calcdef,turbsmag

  real :: turbfac,t2pfac

  CONTAINS

      subroutine sfc_and_turb(getsfc,getpbl,nstep,dt,dosfcflx,cloudvar,qbudget,    &
                   avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp, &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,  &
                   pi0s,rth0s,pi0,rho0,prs0,thv0,th0,qv0,rr0,rf0,rrf0, &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv,        &
                   tsk,thflux,qvflux,cd,ch,cq,u1,v1,s1,t1,tlh,f2d,ustt,ut,vt,st,cm0,   &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,     &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   m11,m12,m13,m22,m23,m33,                          &
                   u0,rru,ua ,ugr ,fsu  ,                            &
                   v0,rrv,va ,vgr ,fsv  ,                            &
                   rrw,wa ,dumw1,dumw2,                              &
                   ppi ,ppten,sten,sadv,                             &
                   tha ,thten,thten1,thterm,qa ,                     &
                   kmh,kmv,khh,khv,cme,csm,ce1,ce2,tkea,tke3d,       &
                   nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2,radsw,radswnet,radlwin, &
                   thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,qrpten,qspten,qgpten, &
                   xkzh,xkzq,xkzm,                                   &
                   tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl,       &
                   qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,            &
                   edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,          &
                   sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D,                  &
                   vdfg,maxmf,nupdraft,ktop_plume,                           &
                   lwten,swten,thraten,tdiag,raincv,pratecv,                 &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,thften,qvften, &
                   lu_index,kpbl2d,psfc,u10,v10,s10,hfx,qfx,         &
                   xland,znt,rznt,ust,stau,tst,qst,z0t,z0q,          &
                   hpbl,wspd,phim,phih,psim,psih,psiq,gz1oz0,br,brcr, &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                    &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                 &
                   CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,    &
                   gsw,glw,chklowq,capg,snowc,snowh,qcg,dsxy,wstar,delta,prkpp,fm,fh,  &
                   charn,msang,scurx,scury,zkmax,cd_out,ch_out,wscale,wscaleu, &
                   zntmp,mznt,swspd,smois,taux,tauy,hpbl2d,evap2d,heat2d,  &
                   mixht,akhs,akms,ct,snow,sice,thz0,qz0,uz0,vz0,u10e,v10e,th10,q10,tshltr,qshltr,pshltr,z0base,zntmyj,lowlyr,ivgtyp,tke_myj,el_myj,tmp2d,tmp1d,  &
                   num_soil_layers,slab_zs,slab_dzs,tslb,tmn,        &
                   tml,t0ml,hml,h0ml,huml,hvml,tmoml,                &
                   cavg,uavg,vavg,savg,gamk,gamwall,kmw,ufw,vfw,u1b,v1b,l2p,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,u2pt,v2pt,kmwk,ufwk,vfwk, &
                   bndy,kbdy,timavg,sfctimavg,                       &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_x,reqs_y, &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,reqk,             &
                   iamsat,out2d,out3d,rtime,update_sfc,dotbud,dotdwrite,restarted,  &
                   dowriteout,doazimwrite,dostat)
      ! end_sfc_and_turb
      use input
      use constants
      use bc_module
      use comm_module
      use misclibs , only : calcksquick
      use sfcphys_module
      use module_sf_sfclay
      use module_sf_slab
      use module_sf_oml
      use module_bl_gfsedmf , only : bl_gfsedmf
      use module_sf_gfdl , only : sf_gfdl
      use cm1libs , only : rslf,rsif
      use module_sf_mynn , only : sfclay_mynn
      use module_bl_mynn_wrapper, only : mynnedmf_wrapper_run
      use module_bl_myjpbl , only : myjpbl
      use module_sf_myjsfc , only : myjsfc
      use turbtend_module , only : turbsz,turbuz,turbvz
      use turbnba_module
      use ib_module

      ! MMM phys:
      use mmm_physics_wrapper_module, only : ysu_wrapper,sf_sfclayrev_wrapper,cu_ntiedtke_wrapper

#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
! Arrays and variables passed into solve

      logical, intent(in) :: getsfc,getpbl
      integer, intent(in) :: nstep
      real, intent(inout) :: dt
      logical, intent(in) :: dosfcflx
      logical, intent(in), dimension(maxq) :: cloudvar
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,qv0,rr0,rf0,rrf0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,rznt,zntmp,ust,stau,tst,qst,z0t,z0q,thflux,qvflux,  &
                                                     cd,ch,cq,u1,v1,s1,t1,psfc,tlh,ustt,ut,vt,st,cm0
      real, intent(in),    dimension(ib:ie,jb:je) :: xland,f2d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,divx,rho,rr,rf,prs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(inout), dimension(ibnba:ienba,jbnba:jenba,kbnba:kenba) :: m11,m12,m13,m22,m23,m33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,ua,ugr,fsu
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,va,vgr,fsv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,dumw1,dumw2
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,ppten,sten,sadv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,thten,thten1,thterm
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: cme,csm,ce1,ce2
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2
      real, intent(inout), dimension(ni,nj) :: radsw,radswnet,radlwin
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,qrpten,qspten,qgpten
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: xkzh,xkzq,xkzm
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn,kbmynn:kemynn) :: tsq,qsq,cov,sh3d,sm3d,el_pbl,qc_bl,qi_bl,cldfra_bl, &
           qWT,qSHEAR,qBUOY,qDISS,dqke,qke_adv,qke,qke3d,edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,  &
           sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D
      real, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: vdfg,maxmf
      integer, intent(inout), dimension(ibmynn:iemynn,jbmynn:jemynn) :: nupdraft,ktop_plume
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,thraten
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: raincv
      real, intent(inout), dimension(ib:ie,jb:je) :: pratecv
      real, intent(inout), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,thften,qvften
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index,kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,s10,hfx,qfx, &
                                      hpbl,wspd,phim,phih,psim,psih,psiq,gz1oz0,br,brcr, &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      gsw,glw,chklowq,capg,snowc,snowh,qcg,dsxy,wstar,delta,prkpp,fm,fh
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: charn,msang,scurx,scury,zkmax,cd_out,ch_out,wscale,wscaleu
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: mznt,swspd,smois,taux,tauy,hpbl2d,evap2d,heat2d
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: mixht,akhs,akms,ct,snow,sice,thz0,qz0,uz0,vz0,u10e,v10e,th10,q10,tshltr,qshltr,pshltr,z0base,zntmyj
      integer, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj) :: lowlyr,ivgtyp
      real, intent(inout), dimension(ibmyj:iemyj,jbmyj:jemyj,kbmyj:kemyj) :: tke_myj,el_myj
      real, intent(inout), dimension(ibtmp:ietmp,kbtmp:ketmp  ,ntmp2d ) :: tmp2d
      real, intent(inout), dimension(ibtmp:ietmp,ntmp1d) :: tmp1d
      integer, intent(in) :: num_soil_layers
      real, intent(in), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      double precision, intent(inout), dimension(kb:ke,3+numq) :: cavg
      real, intent(inout), dimension(kb:ke) :: uavg,vavg,savg,gamk,gamwall,kmw,ufw,vfw,u1b,v1b,l2p,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt,v2pt
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: kmwk,ufwk,vfwk
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_x,reqs_y
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(jmp,kmp) :: vw1,vw2,ve1,ve2
      real, intent(inout), dimension(imp,kmp) :: vs1,vs2,vn1,vn2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      integer, intent(inout) :: reqk
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: iamsat
      real, intent(inout), dimension(ib2d:ie2d,jb2d:je2d,nout2d) :: out2d
      real, intent(inout), dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy

      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      real, intent(in), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg

      real, intent(in) :: rtime
      logical, intent(in) :: update_sfc,dotbud,dotdwrite,restarted,dowriteout,doazimwrite,dostat

!-----------------------------------------------------------------------

      integer :: i,j,k,kk,k2,n,kval
      integer :: isfflx,ifsnow
      logical :: ysu_topdown_pblmix
      real :: ep1,ep2,rovg,dtmin,dz1
      real :: p1000mb,eomeg,stbolt,tem,tem1,tem2,tem3
      real :: pisfc,thgb,tskv,thx,thvx,dthvdz,govrth,za

      logical :: flag_qi
      integer :: p_qi,p_first_scalar
      logical :: disheat,usemyjpbl
      real :: gfs_alpha,var_ric,coef_ric_l,coef_ric_s
      real :: qx
      real :: dtdz
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,rcoef

      real, parameter  ::  oml_relaxation_time  =  -1.0

      ! hwrf vars:
      integer :: ntsflg,ens_random_seed,icoef_sf,iwavecpl
      real :: sfenth,ens_Cdamp
      logical :: lcurr_sf,pert_cd
      integer, dimension(ibl:iel,jbl:jel) :: isltyp
      real :: ens_pblamp
      logical :: pert_pbl

      ! mynn:
      logical :: FLAG_QC,FLAG_QNC,FLAG_QNI,FLAG_QNWFA,FLAG_QNIFA
      logical :: restart,cycling

      ! mmm physics:
      character(len=80) :: errmsg
      integer :: errflg
      real :: rdt
      integer, dimension(ibtmp:ietmp,1) :: int1d

!-----------------------------------------------------------------------

      turbfac = 1.0
      t2pfac = 1.0

    facheck:  &
    IF( idoles .and. sgsmodel.ge.1 )THEN


      if( ramp_sgs.ge.1 )then
        ! ramp up sgs turbulence near beginning of simulation:
        turbfac = (rtime-0.0)/ramp_time
        turbfac = max(0.0,turbfac)
        turbfac = min(1.0,turbfac)
        if(myid.eq.0.and.turbfac.lt.0.999) print *,'  turbfac = ',turbfac

        ! update turbulence model "constants":
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          if( cm0(i,j).gt.cmemin )then
            cme(i,j,k) = max( 0.001 , turbfac*cm0(i,j) )
            ce1(i,j,k) = cme(i,j,k) * c_l * c_l * ( 1.0 / ri_c - 1.0 )
            ce2(i,j,k) = max( 0.0 , cme(i,j,k) * pi * pi - ce1(i,j,k) )
            csm(i,j,k) = ( cme(i,j,k) * cme(i,j,k) * cme(i,j,k) / ( ce1(i,j,k) + ce2(i,j,k) ) )**0.25   ! Smagorinsky constant
          endif
        enddo
        enddo
        enddo
      endif

      if( dot2p )then
        ! ramp up 2-part turbulence model near beginning of simulation:
        t2pfac = (rtime-1800.0)/(1800.0)
        t2pfac = max(0.0,t2pfac)
        t2pfac = min(1.0,t2pfac)
        if(myid.eq.0.and.t2pfac.lt.0.999) print *,'  t2pfac = ',t2pfac
      endif

    ENDIF  facheck


      IF( cm1setup.ge.1 .or. ipbl.ge.1 .or. horizturb.eq.1 .or. idiss.eq.1 .or. output_dissten.eq.1 )THEN

        ! cm1r17:  dissten is defined on w (full) levels:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dissten(i,j,k)=0.0
        enddo
        enddo
        enddo

      ENDIF

      do j=1,nj
      do i=1,ni
        psfc(i,j) = cgs1*prs(i,j,1)+cgs2*prs(i,j,2)+cgs3*prs(i,j,3)
      enddo
      enddo

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
        IF( pertflx.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=jb,je
          do i=ib,ie
            ugr(i,j,k) = ua(i,j,k)-u0(i,j,k)
            vgr(i,j,k) = va(i,j,k)-v0(i,j,k)
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
      ENDIF


      IF( (sfcmodel.eq.2) .or. (sfcmodel.eq.3) .or. (sfcmodel.eq.4) .or. (sfcmodel.eq.6) .or. (sfcmodel.eq.7) .or. use_pbl .or. (oceanmodel.eq.2) )THEN

        ! variables for wrf physics:

        if( sfcmodel.eq.6 )then
          k2 = 2
          if( myid.eq.0 .and. nstep.eq.1 ) print *,'  k2 = ',k2
        else
          k2 = 1
        endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        DO j=1,nj
          do k=1,k2
          do i=1,ni
            ! use ground-relative winds:
            dum1(i,j,k)=0.5*(ugr(i,j,k)+ugr(i+1,j,k))
            dum2(i,j,k)=0.5*(vgr(i,j,k)+vgr(i,j+1,k))
            dum3(i,j,k)=th0(i,j,k)+tha(i,j,k)
            dum7(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
            dum4(i,j,k)=dum3(i,j,k)*dum7(i,j,k)
            thten1(i,j,k)=dum7(i,j,k)
          enddo
          enddo
          if( sfcmodel.eq.6 )then
          if( imoist.eq.1 )then
            do k=1,k2
            do i=1,ni
              sten(i,j,k) = qa(i,j,k,nqc)
            enddo
            enddo
          else
            do k=1,k2
            do i=1,ni
              sten(i,j,k) = 0.0
            enddo
            enddo
          endif
          endif
          do k=1,2
          do i=1,ni
            dum5(i,j,k) = dz*rmh(i,j,k)
          enddo
          enddo
          k = 2
          do i=1,ni
            dum6(i,j,k) = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
          enddo
          ! surface:
          do i=1,ni
            dum6(i,j,1) = psfc(i,j)
          enddo
        ENDDO

        ! dum1 = u at scalars
        ! dum2 = v at scalars
        ! dum3 = th
        ! dum4 = t
        ! dum5 = dz between full levels (m) (ie, at s levels)
        ! dum6 = p3di (p8w) (pressure at w levels)
        ! dum7 = pi

        isfflx = 1
        p1000mb      = 100000.
        EOMEG=7.2921E-5
        STBOLT=5.67051E-8
        ep1 = rv/rd - 1.0
        ep2 = rd/rv
        rovg = rd/g

        ! note:  for RRTMG (radopt=2) gsw and glw arrays are already filled !
        IF( radopt.eq.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=radswnet(i,j)
            glw(i,j)=radlwin(i,j)
          enddo
          enddo
        ELSEIF( radopt.eq.0 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=0.0
            glw(i,j)=0.0
          enddo
          enddo
        ENDIF

      ENDIF


      !-------------------
      ! now, divx stores qv
      !      dum7 stores ql+qi

      ! GHB, 210521:
      IF(imoist.eq.1)THEN
        !$omp parallel do default(shared) private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni
            divx(i,j,k) = qa(i,j,k,nqv)
            dum7(i,j,k) = 0.0
          enddo
          enddo
        ENDDO
      IF( nql1.gt.0 )THEN
        do n=nql1,nql2
        !$omp parallel do default(shared) private(i,j,k)
        DO k=1,nk
          do j=1,nj
          do i=1,ni
            dum7(i,j,k)=dum7(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
        ENDDO
        enddo
      ENDIF
        IF(iice.eq.1)THEN
          do n=nqs1,nqs2
          !$omp parallel do default(shared) private(i,j,k)
          DO k=1,nk
          do j=1,nj
            do i=1,ni
              dum7(i,j,k)=dum7(i,j,k)+qa(i,j,k,n)
            enddo
          enddo
          ENDDO
          enddo
        ENDIF
      ELSE
          !$omp parallel do default(shared) private(i,j,k,n)
          DO k=1,nk
            do j=1,nj
            do i=1,ni
              divx(i,j,k) = 0.0
              dum7(i,j,k) = 0.0
            enddo
            enddo
          ENDDO
      ENDIF

!-----------------------------------------------------------------------
!  hwrf parameters:

      if( sfcmodel.eq.4 )then
        ! do land sfc temperature prediction if ntsflg=1
        ntsflg = 1
        ! sea spray parameter:
        sfenth = 0.0
        ! randomly perturb Cd:
        pert_cd = .false.
        ens_random_seed = 99
        ens_cdamp = 0.2
        ! Option for exchange coefficients in the surface flux scheme:
        icoef_sf = 6
        ! Option for activate coupling to sea surface wave model:
        iwavecpl = 0
        ! Option to include ocean currents in the surface flux calculations:
        lcurr_sf = .false.
        ! water:
        isltyp = 14
      endif

      if( ipbl.eq.3 )then
        p_qi = nqi
        p_first_scalar = 1
        disheat = .true.
        gfs_alpha = -1.0
        !  Flag for using variable Ric
        var_ric = 1.0
        coef_ric_l = 0.16
        coef_ric_s = 0.25
        ens_pblamp = 0.2
        pert_pbl = .false.
      endif

!-----------------------------------------------------------------------

    if(timestats.ge.1) time_turb=time_turb+mytime()

    dosfc:  IF( getsfc )THEN

      bbc3:  IF( bbc.eq.3 )THEN

        !-------------------------------
        ! u1 is u at k=kval (lowest model level by default)
        ! v1 is v at k=kval
        ! s1 is horizontal wind speed at k=kval
        ! (all defined at the scalar point of the staggered grid)
        ! for pertflx=1, account for domain (i.e., surface) motion 
        !                in calculation of wind speed

        kval = 1

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            ! use ground-relative winds:
            u1(i,j) = 0.5*( ugr(i,j,kval) + ugr(i+1,j,kval) )
            v1(i,j) = 0.5*( vgr(i,j,kval) + vgr(i,j+1,kval) )
            s1(i,j) = sqrt(u1(i,j)**2+v1(i,j)**2)
            t1(i,j) = th0(i,j,kval)+tha(i,j,kval)
            ppten(i,j,1) = zh(i,j,kval)
          enddo
          enddo

          IF( terrain_flag )THEN
            !$omp parallel do default(shared)   &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              ppten(i,j,1) = zh(i,j,kval)-zs(i,j)
            enddo
            enddo
          ENDIF

        !-------------------------------
        ! NOTE:
        ! divx stores qv
        ! dum7 stores ql+qi
        ! ppten(k=1) stores height of first model level above surface

        IF( testcase.ge.1 .and. testcase.le.7 )THEN

          ! max gradient method (for simple LES simulations only):

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k,dtdz)
          DO j=1,nj

            do i=1,ni
              dum8(i,j,1) = -1.0e30
            enddo

            do k=nk,2,-1
            do i=1,ni
              dtdz = ( (th0(i,j,k  )+tha(i,j,k  ))*(1.0+repsm1*divx(i,j,k  )-dum7(i,j,k  ))  &
                      -(th0(i,j,k-1)+tha(i,j,k-1))*(1.0+repsm1*divx(i,j,k-1)-dum7(i,j,k-1)) )*rdz*mf(i,j,k)
              if( dtdz .ge. dum8(i,j,1) )then
                dum8(i,j,1) = dtdz
                hpbl(i,j) = zf(i,j,k)
              endif
            enddo
            enddo

          ENDDO

        ELSE

          ! bulk-Ri method:

          IF( ipbl.eq.2 )THEN
            ! (note: for ipbl=1,3,4,5 pbl depth is calculated within the PBL subroutine)
            call gethpbl2(psfc,qsfc,thflux,qvflux,ust,tsk,zh,th0,tha,divx,ugr,vgr,dum8(ib,jb,1),dum8(ib,jb,2),dum8(ib,jb,3),dum8(ib,jb,4),hpbl,thten)
          ENDIF

        ENDIF


        IF( terrain_flag )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            hpbl(i,j) = hpbl(i,j)-zs(i,j)
          enddo
          enddo
         ENDIF

        IF( use_avg_sfc .and. ( sfcmodel.eq.1 .or. sfcmodel.eq.5 ) )THEN
            call getavgsfc(u1,v1,s1,t1,psfc,divx(ib,jb,1),avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp,ugr,vgr)
        ENDIF

        IF( sfcmodel.eq.1 )THEN

          call getcecd(xh,u1,v1,s1,ppten(ib,jb,1),u10,v10,s10,xland,znt,ust,cd,ch,cq,avgsfcu,avgsfcv,avgsfcs,avgsfct)

          if( tbc.eq.3 )then
            kk = nk+1-kval
            do j=1,nj
            do i=1,ni
              ut(i,j) = 0.5*( ugr(i,j,kk) + ugr(i+1,j,kk) )
              vt(i,j) = 0.5*( vgr(i,j,kk) + vgr(i,j+1,kk) )
              st(i,j) = sqrt(ut(i,j)**2+vt(i,j)**2)
            enddo
            enddo
            if( set_znt.eq.1 )then
              do j=1,nj
              do i=1,ni
                ustt(i,j) = max( st(i,j)*karman/alog((zf(i,j,nk+1)-zh(i,j,kk))/cnst_znt) , 1.0e-6 )
              enddo
              enddo
            elseif( set_ust.eq.1 )then
              do j=1,nj
              do i=1,ni
                ustt(i,j) = cnst_ust
              enddo
              enddo
            else
              print *,'  98713 '
              call stopcm1
            endif
          endif

          if(isfcflx.eq.1)then
            call sfcflux(dt,ruh,xf,rvh,pi0s,ch,cq,pi0,thv0,th0,rf0,tsk,thflux,qvflux,mavail, &
                         rho,rf,u1,v1,s1,ppi,tha,divx, &
                         qbudget(8),psfc,u10,v10,s10,qsfc,znt,rtime)
          endif

          ! get sfc diagnostics needed by pbl scheme:
          call sfcdiags(tsk,thflux,qvflux,cd,ch,cq,u1,v1,s1,wspd,        &
                        xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,gz1oz0,  &
                        psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2,fm,fh,  &
                        zs,ppten(ib,jb,1),pi0s,pi0,th0,ppi,tha,rho,rf,divx)
          if( use_pbl .and. ipbl.ne.2 )then
            !$omp parallel do default(shared)   &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              CPMM(i,j)=CP*(1.0+0.8*divx(i,j,1))                                   
              hfx(i,j) = thflux(i,j)*CPMM(i,j)*rf(i,j,1)
              qfx(i,j) = qvflux(i,j)*rf(i,j,1)
            enddo
            enddo
          endif

        ENDIF

        IF( sfcmodel.eq.5 )THEN

          call   cm1most(u1,v1,s1,t1,tst,qst,thflux,qvflux,zol,mol,rmol,  &
                         phim,phih,psim,psih,ppten(ib,jb,1),        &
                         u10,v10,s10,xland,znt,rznt,ust,cd,ch,cq,   &
                         avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,avgsfct,avgsfcq,avgsfcp,rtime,     &
                         tsk,qsfc,psfc,wspd,thv0,th0,tha,rho,divx(ib,jb,1))

          ! variables needed by some PBL codes:
          ep1 = rv/rd - 1.0

          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj
          do i=1,ni
            CPMM(i,j) = CP*(1.0+0.8*divx(i,j,1))                                   
            hfx(i,j) = thflux(i,j)*CPMM(i,j)*rf(i,j,1)
            qfx(i,j) = qvflux(i,j)*rf(i,j,1)
!!!            if( testcase.ne.14 .and. testcase.ne.11 )  &
!!!            qsfc(i,j) = rslf(psfc(i,j),tsk(i,j))
            pisfc = (psfc(i,j)*rp00)**rovcp
            thgb = tsk(i,j)/pisfc
            tskv = thgb*(1.0+ep1*qsfc(i,j))
            thx = th0(i,j,1)+tha(i,j,1)
            thvx = thx*(1.0+EP1*divx(i,j,1))
            DTHVDZ = THVX-TSKV
            govrth = g/thx
            za = ppten(i,j,1)
            br(i,j) = govrth*za*DTHVDZ/(wspd(i,j)**2)
            gz1oz0(i,j) = alog(za*rznt(i,j))
            fm(i,j) = GZ1OZ0(i,j)-PSIM(i,j)
            fh(i,j) = GZ1OZ0(i,j)-PSIH(i,j)
            CHKLOWQ(i,j) = MAVAIL(i,j)
          enddo
          enddo

          if( sfcmodel.eq.7 .or. ipbl.eq.6 )then
!!!            if(myid.eq.0) print *,'  myj sfc params '
              ! for MYJ:
            do j=1,nj
            do i=1,ni
              akms(i,j) = ust(i,j)*karman/fm(i,j)
              akhs(i,j) = ust(i,j)*karman/fh(i,j)
              thz0(i,j) = tsk(i,j)/((psfc(i,j)*rp00)**rovcp)
              qz0(i,j) = qsfc(i,j)
              uz0(i,j) = 0.0
              vz0(i,j) = 0.0
              lh(i,j) = 2.5e6*qvflux(i,j)*rf(i,j,1)
            enddo
            enddo

          endif


        ENDIF


        sfcmod1:  &
        IF( (sfcmodel.eq.2) .or. (sfcmodel.eq.3) .or. (sfcmodel.eq.4) .or. (sfcmodel.eq.6) .or. (sfcmodel.eq.7) )THEN
          ! surface layer:

        sfcmod2:  &
        if( sfcmodel.eq.2 )then
          call SFCLAY(dum1,dum2,dum4,divx,prs,dum5,      &
                       CP,G,ROVCP,RD,XLV,lv1,lv2,PSFC,CHS,CHS2,CQS2,CPMM, &
                       ZNT,UST,hpbl,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                       FM,FH,                                        &
                       XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                       U10,V10,TH2,T2,Q2,rf(ib,jb,1),                &
                       GZ1OZ0,WSPD,BR,ISFFLX,dsxy,                   &
                       SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                       KARMAN,EOMEG,STBOLT,                          &
                       P1000mb,                                      &
                       1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,              &
                       ib ,ie , jb ,je , kb ,ke ,                    &
                       1  ,ni , 1  ,nj , 1  ,nk ,                    &
                       ustm,ck,cka,cd,cda,isftcflx,iz0tlnd,z0t,z0q   )

        elseif( sfcmodel.eq.3 )then

        sfclayrev_j_loop:  &
        do j=1,nj

          do k=1,2
          do i=1,ni
            tmp2d(i,k,1) = dum1(i,j,k)  ! u3d
            tmp2d(i,k,2) = dum2(i,j,k)  ! v3d
            tmp2d(i,k,3) = dum4(i,j,k)  ! t3d
            tmp2d(i,k,4) = divx(i,j,k)  ! qv3d
            tmp2d(i,k,5) = prs(i,j,k)  ! p3d
            tmp2d(i,k,6) = dum5(i,j,k)  ! dz8w
          enddo
          enddo

          do i=1,ni
            tmp1d(i,1) = dsxy(i,j)  ! dx
            tmp1d(i,2) = 0.0        ! water_depth
          enddo

          call sf_sfclayrev_wrapper(               &
                      ux=tmp2d(ibtmp,1,1),  &
                      vx=tmp2d(ibtmp,1,2),  &
                      t1d=tmp2d(ibtmp,1,3),  &
                      qv1d=tmp2d(ibtmp,1,4),  &
                      p1d=tmp2d(ibtmp,1,5),  &
                      dz8w1d=tmp2d(ibtmp,1,6),  &
                      cp=cp,  &
                      g=g,  &
                      rovcp=rovcp,  &
                      r=rd,  &
                      xlv=xlv,  &
                      psfcpa=psfc(1,j),  &
                      chs=chs(1,j),  &
                      chs2=chs2(1,j),  &
                      cqs2=cqs2(1,j),  &
                      cpm=cpmm(1,j),      &
                      pblh=hpbl(1,j),  &
                      rmol=rmol(1,j), &
                      znt=znt(1,j),  &
                      ust=ust(1,j),  &
                      mavail=mavail(1,j),  &
                      zol=zol(1,j),  &
                      mol=mol(1,j),  &
                      regime=regime(1,j),  &
                      psim=psim(1,j),  &
                      psih=psih(1,j), &
                      fm=fm(1,j),  &
                      fh=fh(1,j),                                        &
                      xland=xland(1,j),  &
                      hfx=hfx(1,j),  &
                      qfx=qfx(1,j),  &
                      tsk=tsk(1,j),  &
                      u10=u10(1,j),  &
                      v10=v10(1,j),  &
                      th2=th2(1,j),  &
                      t2=t2(1,j),  &
                      q2=q2(1,j),                            &
                      flhc=flhc(1,j),  &
                      flqc=flqc(1,j),  &
                      qgh=qgh(1,j),  &
                      qsfc=qsfc(1,j),  &
                      lh=lh(1,j),  &
                      gz1oz0=gz1oz0(1,j),  &
                      wspd=wspd(1,j),  &
                      br=br(1,j),  &
                      isfflx=.true.,  &
                      dx=tmp1d(ibtmp,1),                     &
                      svp1=svp1,  &
                      svp2=svp2,  &
                      svp3=svp3,  &
                      svpt0=svpt0,  &
                      ep1=ep1,  &
                      ep2=ep2,                 &
                      karman=karman,  &
                      p1000mb=p1000mb,                                &
                      shalwater_z0=.false.,  &
                      water_depth=tmp1d(ibtmp,2),     &
                      isftcflx=isftcflx,  &
                      iz0tlnd=iz0tlnd,          &
                      scm_force_flux=.false.,  &
                      ustm=ustm(1,j),  &
                      ck=ck(1,j),  &
                      cka=cka(1,j),  &
                      cd=cd(1,j),  &
                      cda=cda(1,j),  &
                      its=1,  &
                      ite=ni,  &
                      errmsg=errmsg,  &
                      errflg=errflg)
        enddo  sfclayrev_j_loop

        elseif( sfcmodel.eq.4 )then
          cd_out = 0.0
          ch_out = 0.0
          mol = 0.0
          zol = 0.0
          out2d = 0.0
          do j=1,nj
          do i=1,ni
            zkmax(i,j) = zh(i,j,1)
          enddo
          enddo
          call SF_GFDL(U3D=dum1,V3D=dum2,T3D=dum4,QV3D=divx,P3D=prs,                           &
                     CP=cp,ROVCP=rovcp,R=rd,XLV=xlv,PSFC=psfc,                                 &
                     CHS=chs,CHS2=chs2,CQS2=cqs2, CPM=cpmm,                                    &
                     DT=dt, SMOIS=smois,num_soil_layers=1,ISLTYP=ISLTYP,ZNT=znt,               &
                     MZNT=mznt,                                                                &
                     UST=ust,PSIM=psim,PSIH=psih,                                              &
                     XLAND=xland,HFX=hfx,QFX=qfx,TAUX=taux,TAUY=tauy,LH=lh,GSW=gsw,GLW=glw,    &
                     TSK=tsk,FLHC=flhc,FLQC=flqc,                                              &
                     QGH=qgh,QSFC=qsfc,U10=u10,V10=v10,                                        &
                     ICOEF_SF=icoef_sf,IWAVECPL=iwavecpl,LCURR_SF=lcurr_sf,                    &
                     CHARN=charn,MSANG=msang,SCURX=scurx, SCURY=scury,                         &
                     pert_Cd=pert_cd, ens_random_seed=ens_random_seed, ens_Cdamp=ens_cdamp,    &
                     GZ1OZ0=gz1oz0,WSPD=swspd,BR=br,ZKMAX=zkmax, ISFFLX=isfflx,                &
                     EP1=ep1,EP2=ep2,KARMAN=karman,NTSFLG=ntsflg,SFENTH=sfenth,                &
                     Cd_out=cd_out,Ch_out=ch_out,mol_out=mol,zol_out=zol,z0t_out=z0t,          &
                     fm_out=fm,fh_out=fh,                                                      &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,                     &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,                       &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk  )
          do j=1,nj
          do i=1,ni
            ust(i,j) = max(ust(i,j),1.0e-8)
            wspd(i,j) = swspd(i,j)
          enddo
          enddo

        elseif( sfcmodel.eq.6 )then
          ! 6 = MYNN sfclayer


        mynn_sf_j_loop:  &
        DO j=1,nj

          do k=1,2
          do i=1,ni
            tmp2d(i,k, 1) = dum1(i,j,k)  ! u3d
            tmp2d(i,k, 2) = dum2(i,j,k)  ! v3d
            tmp2d(i,k, 6) = dum5(i,j,k)  ! dz8w
          enddo
          enddo

          k = 1
          do i=1,ni
            tmp2d(i,k, 3) = dum4(i,j,k)  ! t3d
            tmp2d(i,k, 4) = divx(i,j,k)  ! qv3d
            tmp2d(i,k, 5) =  prs(i,j,k)  ! p3d
            tmp2d(i,k, 7) = dum3(i,j,k)  ! th3d
            tmp2d(i,k, 8) =  rho(i,j,k)  ! rho3d
          enddo

          call SFCLAY_mynn(                                                                          &
                     U3D=tmp2d(ibtmp,kbtmp,1),                                                     &
                     V3D=tmp2d(ibtmp,kbtmp,2),                                                     &
                     T3D=tmp2d(ibtmp,kbtmp,3),                                                     &
                     QV3D=tmp2d(ibtmp,kbtmp,4),                                                    &
                     P3D=tmp2d(ibtmp,kbtmp,5),                                                     &
                     dz8w=tmp2d(ibtmp,kbtmp,6),                                                    &
                     th3d=tmp2d(ibtmp,kbtmp,7),                                                    &
                     rho3d=tmp2d(ibtmp,kbtmp,8),                                                   &
                     CP=cp,G=g,ROVCP=rovcp,R=rd,XLV=xlv,PSFCPA=psfc(1,j),                            &
                     CHS=chs(1,j),CHS2=chs2(1,j),CQS2=cqs2(1,j),CPM=cpmm(1,j),                       &
                     ZNT=znt(1,j),UST=ust(1,j),PBLH=hpbl(1,j),MAVAIL=mavail(1,j),                    &
                     ZOL=zol(1,j),MOL=mol(1,j),REGIME=regime(1,j),PSIM=psim(1,j),PSIH=psih(1,j),     &
                     XLAND=xland(1,j),HFX=hfx(1,j),QFX=qfx(1,j),LH=lh(1,j),TSK=tsk(1,j),             &
                     FLHC=flhc(1,j),FLQC=flqc(1,j),QGH=qgh(1,j),QSFC=qsfc(1,j),RMOL=rmol(1,j),       &
                     U10=u10(1,j),V10=v10(1,j),TH2=th2(1,j),T2=t2(1,j),Q2=q2(1,j),SNOWH=snowh(1,j),  &
                     GZ1OZ0=gz1oz0(1,j),WSPD=wspd(1,j),BR=br(1,j),ISFFLX=isfflx,DX=dx,               &
                     SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0,EP1=ep1,EP2=ep2,                      &
                     KARMAN=karman,itimestep=nstep,                                                  &
                     ch=ch(1,j),                                                                     &
                     qcg=qcg(1,j),                                                                   &
                     spp_pbl=0,                                                                      &
                     ids=1  ,ide=ni+1 , jds=1  ,jde=2    , kds=1  ,kde=nk+1 ,                        &
                     ims=1  ,ime=ni   , jms=1  ,jme=1    , kms=1  ,kme=nk ,                          &
                     its=1  ,ite=ni   , jts=1  ,jte=1    , kts=1  ,kte=nk  ,                         &
                     ustm=ustm(1,j),ck=ck(1,j),cka=cka(1,j),cd=cd(1,j),cda=cda(1,j),                 &
                     isftcflx=isftcflx,iz0tlnd=iz0tlnd)
        ENDDO  mynn_sf_j_loop


        elseif( sfcmodel.eq.7 )then
          ! 7 = MYJ sfclayer

          if( ipbl.eq.6 )then
            usemyjpbl = .true.
          else
            usemyjpbl = .false.
          endif

        do j=1,nj

          do k=1,nk
          do i=1,ni
            tmp2d(i,k, 1) = dz*rmh(i,j,k)  ! dz
            tmp2d(i,k, 2) = prs(i,j,k)     ! pmid
            tmp2d(i,k, 4) = (th0(i,j,k)+tha(i,j,k))  ! th
            tmp2d(i,k, 5) = (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k))  ! t
            tmp2d(i,k, 6) = 0.0
            tmp2d(i,k, 7) = 0.0
            tmp2d(i,k, 8) = 0.5*(ugr(i,j,k)+ugr(i+1,j,k))  ! u
            tmp2d(i,k, 9) = 0.5*(vgr(i,j,k)+vgr(i,j+1,k))  ! v
            tmp2d(i,k,10) = tke_myj(i,j,k)  ! tke_myj
          enddo
          enddo
          ! pint:
          do k=2,nk
          do i=1,ni   ! cm1r21.1 bug fix (this loop was neglected before)
            tmp2d(i,k,3) = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
          enddo
          enddo
          do i=1,ni
            tmp2d(i,1,3) = psfc(i,j)
            tmp2d(i,nk+1,3) = cgt1*prs(i,j,nk)+cgt2*prs(i,j,nk-1)+cgt3*prs(i,j,nk-2)
          enddo
          if( imoist.eq.1 .and. nqv.ge.1 )then
            do k=1,nk
            do i=1,ni
              tmp2d(i,k, 6) = qa(i,j,k,nqv)
            enddo
            enddo
          endif
          if( imoist.eq.1 .and. nqc.ge.1 )then
            do k=1,nk
            do i=1,ni
              tmp2d(i,k, 7) = qa(i,j,k,nqc)
            enddo
            enddo
          endif

          call   MYJSFC(ITIMESTEP=nstep,  &
                        HT=zs(1,j),  &
                        DZ=tmp2d(ibtmp,kbtmp,1),                               &
                        PMID=tmp2d(ibtmp,kbtmp,2),  &
                        PINT=tmp2d(ibtmp,kbtmp,3),  &
                        TH=tmp2d(ibtmp,kbtmp,4),  &
                        T=tmp2d(ibtmp,kbtmp,5),  &
                        QV=tmp2d(ibtmp,kbtmp,6),  &
                        QC=tmp2d(ibtmp,kbtmp,7),  &
                        U=tmp2d(ibtmp,kbtmp,8),  &
                        V=tmp2d(ibtmp,kbtmp,9),  &
                        Q2=tmp2d(ibtmp,kbtmp,10)                    &
                       ,TSK=tsk(1,j),QSFC=qsfc(1,j),THZ0=thz0(1,j),QZ0=qz0(1,j),UZ0=uz0(1,j),VZ0=vz0(1,j)                      &
                       ,LOWLYR=lowlyr(1,j),XLAND=xland(1,j),IVGTYP=ivgtyp(1,j),ISURBAN=0,IZ0TLND=iz0tlnd            &
                       ,USTAR=ust(1,j),ZNT=zntmyj(1,j),Z0BASE=z0base(1,j),PBLH=hpbl(1,j),MAVAIL=mavail(1,j),RMOL=rmol(1,j)              &
                       ,AKHS=akhs(1,j),AKMS=akms(1,j)                                      &
                       ,RIB=br(1,j)                                             &
                       ,CHS=chs(1,j),CHS2=chs2(1,j),CQS2=cqs2(1,j),HFX=hfx(1,j),QFX=qfx(1,j),FLX_LH=lh(1,j),FLHC=flhc(1,j),FLQC=flqc(1,j)         &
                       ,QGH=qgh(1,j),CPM=cpmm(1,j),CT=ct(1,j)                                     &
                       ,U10=u10(1,j),V10=v10(1,j),T02=t2(1,j),TH02=th2(1,j),TSHLTR=tshltr(1,j),TH10=th10(1,j),Q02=q2(1,j),QSHLTR=qshltr(1,j),Q10=q10(1,j),PSHLTR=pshltr(1,j)          &
                       ,P1000mb=p1000mb,U10E=u10e(1,j),V10E=v10e(1,j),                             &
                        psim=psim(1,j),psih=psih(1,j),fm=fm(1,j),fhh=fh(1,j),z0out=znt(1,j),usemyjpbl=usemyjpbl,rhosfc=rf(1,j,1),isftcflx=isftcflx,     &
                     ids=1  ,ide=ni+1 , jds= 1 ,jde= 2   , kds=1  ,kde=nk+1 ,                                   &
                     ims=1  ,ime=ni   , jms= 1 ,jme= 1   , kms=1  ,kme=nk+1 ,                                     &
                     its=1  ,ite=ni   , jts= 1 ,jte= 1   , kts=1  ,kte=nk  )
          do i=1,ni
            br(i,j) = max( -10.0 , br(i,j) )
            br(i,j) = min(  10.0 , br(i,j) )
            gz1oz0(I,J)=ALOG(ppten(i,j,1)/znt(i,j))
          enddo
        enddo

        endif  sfcmod2


          ifsnow = 0
          dtmin = dt/60.0

        IF( update_sfc .and. testcase.ne.11 )THEN
          ! slab scheme (MM5/WRF):
          call SLAB(dum4,divx,prs,FLHC,FLQC,                        &
                       PSFC,XLAND,TMN,HFX,QFX,LH,TSK,QSFC,CHKLOWQ,  &
                       GSW,GLW,CAPG,THC,SNOWC,EMISS,MAVAIL,         &
                       DT,ROVCP,XLV,lv1,lv2,DTMIN,IFSNOW,           &
                       SVP1,SVP2,SVP3,SVPT0,EP2,                    &
                       KARMAN,EOMEG,STBOLT,                         &
                       TSLB,slab_ZS,slab_DZS,num_soil_layers, .true. ,       &
                       P1000mb,                                     &
                         1, ni+1,   1, nj+1,   1, nk+1,             &
                        ib, ie,  jb, je,  kb, ke,                   &
                         1, ni,   1, nj,   1, nk                    )
        ELSE
        if( sfcmodel.ne.7 )then
          ! dont update tsk, but diagnose qsfc:
          !  ?  ... should this be land pts only? !
          do j=1,nj
          do i=1,ni
            qx = divx(i,j,1)
            if ( FLQC(i,j) .ne. 0.) then
               QSFC(i,j)=QX+QFX(i,j)/FLQC(i,j)
            else
               QSFC(i,j) = QX
            end if
            CHKLOWQ(i,j)=MAVAIL(i,j)
          enddo
          enddo
        endif
        ENDIF

          ! put WRF parameters into CM1 arrays:
          if( sfcmodel.eq.2 .or. sfcmodel.eq.3 .or. sfcmodel.eq.6 .or. sfcmodel.eq.7 )then
            !$omp parallel do default(shared)   &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              ch(i,j) = chs2(i,j)
              cq(i,j) = cqs2(i,j)
              s10(i,j) = sqrt( u10(i,j)**2 + v10(i,j)**2 )
            enddo
            enddo
          endif
          if( sfcmodel.eq.4 )then
            do j=1,nj
            do i=1,ni
              cd(i,j) = cd_out(i,j)
              ch(i,j) = ch_out(i,j)
              cq(i,j) = ch_out(i,j)
              s10(i,j) = sqrt( u10(i,j)**2 + v10(i,j)**2 )
              if( abs(mol(i,j)).le.smeps )then
                rmol(i,j) = sign( 1.0e20 , mol(i,j) )
              else
                rmol(i,j) = 1.0/mol(i,j)
              endif
            enddo
            enddo
          endif
          if( sfcmodel.eq.6 )then
            !$omp parallel do default(shared)   &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              ch(i,j) = ck(i,j)
              cq(i,j) = ck(i,j)
            enddo
            enddo
          endif
          if( sfcmodel.eq.7 )then
            !$omp parallel do default(shared)   &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              ustm(i,j) = ust(i,j)
              ! from wrf surface_driver:
              wspd(i,j) = MAX(SQRT(dum1(i,j,1)**2+dum2(i,j,1)**2),0.001)
!!!              ch(i,j) = chs(i,j)
!!!              cq(i,j) = akhs(i,j)
!!!              cd(i,j) = akms(i,j)
              if( abs(rmol(i,j)).le.smeps )then
                mol(i,j) = sign( 1.0e20 , rmol(i,j) )
              else
                mol(i,j) = 1.0/rmol(i,j)
              endif
              ! wspd:
              tem = max( 0.1 , ust(i,j)*rkarman*fm(i,j) )
              wspd(i,j) = max( wspd(i,j) , tem )
!!!              cd(i,j) = ust(i,j)*ust(i,j)/max(1.0e-20,tem*tem)
              cd(i,j) = ust(i,j)*ust(i,j)/max(1.0e-20,s10(i,j)*s10(i,j))
            enddo
            enddo
          endif
          IF( dosfcflx .or. output_sfcflx.eq.1 )THEN
            !$omp parallel do default(shared)   &
            !$omp private(i,j)
            do j=1,nj
            do i=1,ni
              thflux(i,j) = hfx(i,j)/(CPMM(i,j)*rf(i,j,1))
              qvflux(i,j) = qfx(i,j)/rf(i,j,1)
            enddo
            enddo
          ENDIF

        ENDIF  sfcmod1

      ENDIF  bbc3
      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

    ELSE

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ... skipping sfc stuff ... '
      if(dowr) write(outfile,*)

    ENDIF  dosfc

!---------------------------------------------
! bc/comms (very important):

  bbc3b:  &
  IF( bbc.eq.3 )THEN
    !-------------!
    call bc2d(ust)
#ifdef MPI
    call comm_1s2d_start(ust,uw31(1,1,1),uw32(1,1,1),ue31(1,1,1),ue32(1,1,1),  &
                             us31(1,1,1),us32(1,1,1),un31(1,1,1),un32(1,1,1),reqs_s)
#endif
    call bc2d(u1)
#ifdef MPI
    call comm_1s2d_start(u1 ,uw31(1,1,2),uw32(1,1,2),ue31(1,1,2),ue32(1,1,2),  &
                             us31(1,1,2),us32(1,1,2),un31(1,1,2),un32(1,1,2),reqs_u)
#endif
    call bc2d(v1)
#ifdef MPI
    call comm_1s2d_start(v1 ,uw31(1,1,3),uw32(1,1,3),ue31(1,1,3),ue32(1,1,3),  &
                             us31(1,1,3),us32(1,1,3),un31(1,1,3),un32(1,1,3),reqs_v)
#endif
    call bc2d(s1)
#ifdef MPI
    call comm_1s2d_start(s1 ,uw31(1,1,4),uw32(1,1,4),ue31(1,1,4),ue32(1,1,4),  &
                             us31(1,1,4),us32(1,1,4),un31(1,1,4),un32(1,1,4),reqs_w)
#endif
    call bc2d(znt)
#ifdef MPI
    call comm_1s2d_start(znt,uw31(1,1,5),uw32(1,1,5),ue31(1,1,5),ue32(1,1,5),  &
                             us31(1,1,5),us32(1,1,5),un31(1,1,5),un32(1,1,5),reqs_p)
#endif
  if( sfcmodel.eq.4 )then
    call bc2d(mznt)
#ifdef MPI
    call comm_1s2d_start(mznt,uw31(1,1,6),uw32(1,1,6),ue31(1,1,6),ue32(1,1,6),  &
                              us31(1,1,6),us32(1,1,6),un31(1,1,6),un32(1,1,6),reqs_x)
#endif
  endif
    call bc2d(wspd)
#ifdef MPI
    call comm_1s2d_start(wspd,uw31(1,1,7),uw32(1,1,7),ue31(1,1,7),ue32(1,1,7),  &
                              us31(1,1,7),us32(1,1,7),un31(1,1,7),un32(1,1,7),reqs_y)
#endif
    !-------------!

#ifdef MPI
    !-------------!
    call comm_1s2d_end(ust,uw31(1,1,1),uw32(1,1,1),ue31(1,1,1),ue32(1,1,1),  &
                           us31(1,1,1),us32(1,1,1),un31(1,1,1),un32(1,1,1),reqs_s)
    call bcs2_2d(ust)

    call comm_1s2d_end(u1 ,uw31(1,1,2),uw32(1,1,2),ue31(1,1,2),ue32(1,1,2),  &
                           us31(1,1,2),us32(1,1,2),un31(1,1,2),un32(1,1,2),reqs_u)
    call bcs2_2d(u1 )

    call comm_1s2d_end(v1 ,uw31(1,1,3),uw32(1,1,3),ue31(1,1,3),ue32(1,1,3),  &
                           us31(1,1,3),us32(1,1,3),un31(1,1,3),un32(1,1,3),reqs_v)
    call bcs2_2d(v1 )

    call comm_1s2d_end(s1 ,uw31(1,1,4),uw32(1,1,4),ue31(1,1,4),ue32(1,1,4),  &
                           us31(1,1,4),us32(1,1,4),un31(1,1,4),un32(1,1,4),reqs_w)
    call bcs2_2d(s1 )

    call comm_1s2d_end(znt,uw31(1,1,5),uw32(1,1,5),ue31(1,1,5),ue32(1,1,5),  &
                           us31(1,1,5),us32(1,1,5),un31(1,1,5),un32(1,1,5),reqs_p)
    call bcs2_2d(znt)

  if( sfcmodel.eq.4 )then
    call comm_1s2d_end(mznt,uw31(1,1,6),uw32(1,1,6),ue31(1,1,6),ue32(1,1,6),  &
                            us31(1,1,6),us32(1,1,6),un31(1,1,6),un32(1,1,6),reqs_x)
    call bcs2_2d(mznt)
  endif

    call comm_1s2d_end(wspd,uw31(1,1,7),uw32(1,1,7),ue31(1,1,7),ue32(1,1,7),  &
                            us31(1,1,7),us32(1,1,7),un31(1,1,7),un32(1,1,7),reqs_y)
    call bcs2_2d(wspd)
    !-------------!

    !-------------!
    call comm_2d_corner(ust)
    call comm_2d_corner(u1)
    call comm_2d_corner(v1)
    call comm_2d_corner(s1)
    call comm_2d_corner(znt)
  if( sfcmodel.eq.4 )then
    call comm_2d_corner(mznt)
  endif
    call comm_2d_corner(wspd)
    !-------------!
#endif

    ! account for annoying GFDL naming convention:
    IF( sfcmodel.eq.4 )THEN
      do j=0,nj+1
      do i=0,ni+1
        zntmp(i,j) = mznt(i,j)
      enddo
      enddo
    ELSE
      do j=0,nj+1
      do i=0,ni+1
        zntmp(i,j) = znt(i,j)
      enddo
      enddo
    ENDIF
    if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

  ENDIF  bbc3b

!-------------------------------------------------------------------

  tbc3b:  &
  IF( tbc.eq.3 )THEN
    !-------------!
    call bc2d(ustt)
#ifdef MPI
    call comm_1s2d_start(ustt,uw31(1,1,1),uw32(1,1,1),ue31(1,1,1),ue32(1,1,1),  &
                              us31(1,1,1),us32(1,1,1),un31(1,1,1),un32(1,1,1),reqs_s)
#endif
    call bc2d(ut)
#ifdef MPI
    call comm_1s2d_start(ut ,uw31(1,1,2),uw32(1,1,2),ue31(1,1,2),ue32(1,1,2),  &
                             us31(1,1,2),us32(1,1,2),un31(1,1,2),un32(1,1,2),reqs_u)
#endif
    call bc2d(vt)
#ifdef MPI
    call comm_1s2d_start(vt ,uw31(1,1,3),uw32(1,1,3),ue31(1,1,3),ue32(1,1,3),  &
                             us31(1,1,3),us32(1,1,3),un31(1,1,3),un32(1,1,3),reqs_v)
#endif
    call bc2d(st)
#ifdef MPI
    call comm_1s2d_start(st ,uw31(1,1,4),uw32(1,1,4),ue31(1,1,4),ue32(1,1,4),  &
                             us31(1,1,4),us32(1,1,4),un31(1,1,4),un32(1,1,4),reqs_w)
#endif
    !-------------!

#ifdef MPI
    !-------------!
    call comm_1s2d_end(ustt,uw31(1,1,1),uw32(1,1,1),ue31(1,1,1),ue32(1,1,1),  &
                            us31(1,1,1),us32(1,1,1),un31(1,1,1),un32(1,1,1),reqs_s)
    call bcs2_2d(ustt)

    call comm_1s2d_end(ut ,uw31(1,1,2),uw32(1,1,2),ue31(1,1,2),ue32(1,1,2),  &
                           us31(1,1,2),us32(1,1,2),un31(1,1,2),un32(1,1,2),reqs_u)
    call bcs2_2d(ut )

    call comm_1s2d_end(vt ,uw31(1,1,3),uw32(1,1,3),ue31(1,1,3),ue32(1,1,3),  &
                           us31(1,1,3),us32(1,1,3),un31(1,1,3),un32(1,1,3),reqs_v)
    call bcs2_2d(vt )

    call comm_1s2d_end(st ,uw31(1,1,4),uw32(1,1,4),ue31(1,1,4),ue32(1,1,4),  &
                           us31(1,1,4),us32(1,1,4),un31(1,1,4),un32(1,1,4),reqs_w)
    call bcs2_2d(st )

    !-------------!

    !-------------!
    call comm_2d_corner(ustt)
    call comm_2d_corner(ut)
    call comm_2d_corner(vt)
    call comm_2d_corner(st)
    !-------------!
#endif

  ENDIF  tbc3b

!-------------------------------------------------------------------
! simple ocean mixed layer model based Pollard, Rhines and Thompson (1973)
!   (from WRF)

    IF(oceanmodel.eq.2)THEN
    IF( update_sfc )THEN
      if( getsfc )then

        CALL oceanml(tml,t0ml,hml,h0ml,huml,hvml,ust,dum1,dum2, &
                     tmoml,f2d,g,oml_gamma,                     &
                    OML_RELAXATION_TIME,                        &
                     xland,hfx,lh,tsk,gsw,glw,emiss,            &
                     dt,STBOLT,                                 &
                       1, ni+1,   1, nj+1,   1, nk+1,           &
                      ib, ie,  jb, je,  kb, ke,                 &
                       1, ni,   1, nj,   1, nk                  )

        if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      endif
    ENDIF
    ENDIF

!---------------------------------------------

      IF( sgsmodel.ge.1 .or. output_nm.eq.1 .or. ipbl.ge.1 )THEN
        ! squared Brunt-Vaisala frequency:
        iamsat = .false.
        call calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,nm,dum1,dum2,dum3,dum4,dum5,dum6,   &
                    prs,ppi,tha,qa,iamsat)
      ENDIF

      IF( cm1setup.ge.1 .or. output_def.eq.1 .or. ipbl.ge.1 .or. horizturb.eq.1 )THEN
        ! deformation:
        call calcdef(    rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                     xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,mh,c1,c2,mf,defv,defh,  &
                     dum1,dum2,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy,rho,rr,rf)

      ENDIF

!--------------------------------------
!  LES subgrid models:

    les_sgs:  &
    IF( idoles )THEN

      sgsoption:  &
      IF( sgsmodel.eq.1 .or. sgsmodel.eq.3 .or. sgsmodel.eq.4 .or. sgsmodel.eq.5 )THEN

         ! Deardorff TKE:
        call     tkekm(nstep,rtime,dt,ruh,rvh,rmh,zf,mf,rmf,zntmp,ust,ustt,rf, &
                         nm,defv,defh,dum1,dum2,dum3,dum4  ,dum5   ,           &
                         kmh,kmv,khh,khv,cme,csm,ce1,ce2,tkea,lenscl,dissten,out3d, &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                   &
                         kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),       &
                         ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),       &
                         kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),       &
                         ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2),       &
                         kw1(1,1,3),kw2(1,1,3),ke1(1,1,3),ke2(1,1,3),       &
                         ks1(1,1,3),ks2(1,1,3),kn1(1,1,3),kn2(1,1,3),       &
                         kw1(1,1,4),kw2(1,1,4),ke1(1,1,4),ke2(1,1,4),       &
                         ks1(1,1,4),ks2(1,1,4),kn1(1,1,4),kn2(1,1,4))

        if( do_ib )then
          do j=0,nj+2
          do i=0,ni+2
            if( kbdy(i,j).gt.1 )then
              tkea(i,j,1) = 0.0
               kmv(i,j,1) = 0.0
               kmh(i,j,1) = 0.0
               khv(i,j,1) = 0.0
               khh(i,j,1) = 0.0
            endif
          enddo
          enddo
        endif

        call     tkebc(tkea,kmh,kmv,khh,khv,zntmp,ust,ustt,          &
                       kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),  &
                       ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),  &
                       kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),  &
                       ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))

        do j=0,nj+1
        do i=0,ni+1
          tke3d(i,j,1) = tkea(i,j,1)
          tke3d(i,j,nk+1) = tkea(i,j,nk+1)
        enddo
        enddo

        IF( sgsmodel.eq.5 )THEN
          ! Nonlinear Backscatter and Anisotropy (NBA) model:
          ! TKE version:
          call   turbnba(nstep,uh,ruh,uf,ruf,vh,rvh,vf,rvf,mh,rmh,mf,rmf,zf,c1,c2,rho,rf,zntmp,ust,cm0,  &
                         dum1,dum2,dum3,dum4,dum5,dum6,            &
                         dum7,dum8,thten,lenscl,ppten ,dum9,cme,   &   ! cm1r21.1 bug fix: replace divx with dum9
                         m11,m12,m13,m22,m23,m33,                  &
                         ua ,va ,wa ,tkea ,nm,                     &
                         kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,reqs_s,   &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        ENDIF

      ELSEIF(sgsmodel.eq.2)THEN

         ! Smagorinsky:
        call     turbsmag(nstep,dt,ruh,rvh,rmh,mf,rmf,th0,rf,              &
                          nm,defv,defh,dum4,dum5  ,thten1,zf,zntmp,ust,csm, &
                          kmh,kmv,khh,khv,lenscl,dissten,                  &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                          kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                          ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                          kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                          ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))

      ELSEIF( sgsmodel.eq.6 )THEN

        ! Nonlinear Backscatter and Anisotropy (NBA) model:
        ! Smagorinsky version:
        call     turbnba2(nstep,uh,ruh,uf,ruf,vh,rvh,vf,rvf,mh,rmh,mf,rmf,zf,c1,c2,rho,rf,zntmp,ust,cm0,  &
                         dum1,dum2,dum3,dum4,dum5,dum6,            &
                         dum7,dum8,thten,lenscl,ppten ,dum9,cme,   &   ! cm1r21.1 bug fix: replace divx with dum9
                         m11,m12,m13,m22,m23,m33,                  &
                         ua ,va ,wa ,tkea ,nm,                     &
                         kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,reqs_s,   &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)

      ELSEIF( sgsmodel.eq.0 )THEN

        ! for LES modeling without a subgrid turbulence model:
        ! (sometimes called implicit les, or ILES)

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=2,nk+1
        do j=0,nj+1
        do i=0,ni+1
          t11(i,j,k) = 0.0
          t22(i,j,k) = 0.0
          t33(i,j,k) = 0.0
          t12(i,j,k) = 0.0
          t13(i,j,k) = 0.0
          t23(i,j,k) = 0.0
        enddo
        enddo
        enddo

      ELSE

        print *,'  72383 '
        call stopcm1

      ENDIF  sgsoption

      !-------------

      ift2p:  &
      IF( dot2p )THEN
        ! Use 2-part subgrid turbulence model:

        rstrt:  &
        if( restarted )then

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... skipping t2p stuff ... '
          if(dowr) write(outfile,*)

        else

          if( sgsmodel.eq.3 )then
            call  t2psmm(dt,rtime,xf,rxf,c1,c2,                                                &
                         zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,                                  &
                         uavg,vavg,savg,l2p,kmw,ufw,vfw,gamk,gamwall,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,                                     &
                         dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                              &
                         ust,zntmp,stau,mol,rho,rr,rf,kmv,kmh,                                 &
                         ugr ,vgr ,uf,vf,arf1,arf2,                                            &
                         u1,v1,s1,u1b,v1b,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,           &
                         reqs_s,uw31,uw32,ue31,ue32,us31,us32,un31,un32)
          elseif( sgsmodel.eq.4 )then
          if( t2p_avg.eq.1 )then
            call getgamk(dt,rtime,xf,rxf,c1,c2,zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,            &
                         uavg,vavg,savg,l2p,kmw,gamwall,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,dum1,dum2,dum3,dum4,                 &
                         ust,zntmp,stau,rho,rr,rf,kmh,kmv,                                     &
                         ugr ,fsu  ,u2pt,vgr ,fsv  ,v2pt,uf,vf,arf1,arf2,                      &
                         u1,v1,s1,u1b,v1b,                                                     &
                         avgsfcu,avgsfcv,avgsfcs,                                              &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
            call t2pcode(dt,rtime,xf,rxf,c1,c2,zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,            &
                         uavg,vavg,savg,l2p,kmw,ufw,vfw,gamwall,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,dum1,dum2,dum3,dum4,                 &
                         ust,zntmp,stau,mol,rho,rr,rf,kmh,kmv,                                 &
                         ugr ,fsu  ,u2pt,vgr ,fsv  ,v2pt,uf,vf,arf1,arf2,                      &
                         u1,v1,s1,u1b,v1b,                                                     &
                         avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,                            &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
            do k=1,ntwk
            do j=0,nj+1
            do i=0,ni+1
              ufwk(i,j,k) = ufw(k)
              vfwk(i,j,k) = vfw(k)
              kmwk(i,j,k) = kmw(k)
            enddo
            enddo
            enddo
          elseif( t2p_avg.eq.2 )then
            call t2pcodetavg(dt,rtime,xf,rxf,c1,c2,zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,            &
                         uavg,vavg,savg,l2p,kmw,ufw,vfw,gamwall,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,dum1,dum2,dum3,dum4,                 &
                         ust,zntmp,stau,mol,s10,rho,rr,rf,kmh,kmv,                             &
                         ugr ,fsu  ,u2pt,vgr ,fsv  ,v2pt,uf,vf,arf1,arf2,kmwk,ufwk,vfwk,       &
                         u1,v1,s1,u1b,v1b,                                                     &
                         avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,                            &
                         timavg,sfctimavg,                                                     &
                         uw31,uw32,ue31,ue32,us31,us32,un31,un32,                              &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)
          else
            stop 92879
          endif
          else
            stop 54321
          endif

          if( cm1setup.eq.4 )then
            do k=1,ntwk
            do j=0,nj+1
            do i=0,ni+1
              if( cm0(i,j).le.cmemin )then
                ! zero-out 2-part vars outside of LES subdomain:
                ufwk(i,j,k) = 0.0
                vfwk(i,j,k) = 0.0
                kmwk(i,j,k) = 0.0
              endif
            enddo
            enddo
            enddo
          endif

        endif  rstrt

      ENDIF  ift2p

      !-------------

    ENDIF  les_sgs

!--------------------------------------
!  Simple Smagorinsky-like turbulence scheme for horiz direction:

      IF( (cm1setup.eq.2.or.cm1setup.eq.4) .and. horizturb.eq.1 )THEN
        call turbparam_horiz(nstep,zf,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                      nm,defv,defh,dum4,kmh,khh,cm0,dissten,out3d,zs,zntmp,ust,xland,psfc,tlh, &
                      nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                      kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                      ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                      kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                      ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))
      ENDIF

      IF( (cm1setup.eq.2.or.cm1setup.eq.4) .and. ipbl.eq.2 )THEN
        call turbparam_vert(nstep,zf,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                      nm,defv,defh,dum4,kmv,khv,cm0,dissten,out3d,zs,zntmp,ust,xland, &
                      nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                      kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                      ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                      kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                      ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))
      ENDIF

!--------------------------------------
!   DNS (constant viscosity)

      IF( cm1setup.eq.3 )THEN

        tem = viscosity/pr_num

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=0,nj+1
        do i=0,ni+1
          kmh(i,j,k) = viscosity
          kmv(i,j,k) = viscosity
          khh(i,j,k) = tem
          khv(i,j,k) = tem
        enddo
        enddo
        enddo

      ENDIF

!--------------------------------------
!  Check for numerical stability:

      if( cm1setup.ge.1 .and. adapt_dt.eq.1 )then
        call calcksquick(dt,uh,vh,mf,kmh,kmv,khh,khv,reqk)
      endif

!--------------------------------------

    tau:  &
    IF( idoles .or. ipbl.eq.2 .or. horizturb.eq.1 .or. cm1setup.eq.3 )THEN

      if( do_ib )then
#ifdef MPI
        call getcorneru(ua,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcu2(ua)
        call getcornerv(va,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcv2(va)
        call getcornerw(wa,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcw2(wa)
#endif
      endif

        !  now, get turbulent stresses:
      IF( sgsmodel.eq.3 )THEN
        ! SMM94 two-part model:

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=0,nj+1
        do i=0,ni+1
          dumw1(i,j,k) = gamk(k)*kmh(i,j,k)
          dumw2(i,j,k) = gamk(k)*kmv(i,j,k)
        enddo
        enddo
        enddo
        call     gettau(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,ustt,ut,vt,st,rf,mf,dum1,dum2,dum3,dum4,dum5,dum6, &
                        dumw1,dumw2,t11,t12,t13,t22,t23,t33,ua,ugr,va,vgr,wa,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,bndy,kbdy)

      ELSE

        ! all other sgsmodels:
        call     gettau(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,ustt,ut,vt,st,rf,mf,dum1,dum2,dum3,dum4,dum5,dum6, &
                        kmh,kmv,t11,t12,t13,t22,t23,t33,ua,ugr,va,vgr,wa,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,bndy,kbdy)

      ENDIF

      if( do_ib )then
        if(     sgsmodel.le.4 )then
          call drag_obstacles(xh,yh,zh,zf,rho,rf,dum1,dum2,dum3,dum4,dum5,dum6,t11,t12,t13,t22,t23,t33,ua,va,wa,kbdy)
        elseif( sgsmodel.eq.5 .or. sgsmodel.eq.6 )then
          call drag_obstacles(xh,yh,zh,zf,rho,rf,dum1,dum2,dum3,dum4,dum5,dum6,m11,m12,m13,m22,m23,m33,ua,va,wa,kbdy)
        endif
      endif

      IF( sgsmodel.eq.5 .or. sgsmodel.eq.6 )THEN
        do j=1,nj+1
        do i=1,ni+1
          m13(i,j,nk+1) = m13(i,j,nk)
          m23(i,j,nk+1) = m23(i,j,nk)
        enddo
        enddo
      ENDIF

      IF( bbc.eq.3 .and. ( sgsmodel.eq.5 .or. sgsmodel.eq.6 ) )THEN
          ! for NBA scheme, set m13 and m23 at boundaries:
        do j=1,nj+1
        do i=1,ni+1
          m13(i,j,1) = t13(i,j,1)
          m23(i,j,1) = t23(i,j,1)
        enddo
        enddo
      ENDIF
      IF( tbc.eq.3 .and. ( sgsmodel.eq.5 .or. sgsmodel.eq.6 ) )THEN
          ! for NBA scheme, set m13 and m23 at boundaries:
        do j=1,nj+1
        do i=1,ni+1
          m13(i,j,nk+1) = t13(i,j,nk+1)
          m23(i,j,nk+1) = t23(i,j,nk+1)
        enddo
        enddo
      ENDIF

    ENDIF  tau

    IF( ( idoles .and. sgsmodel.ge.1 ) .or. (cm1setup.eq.3.and.dotdwrite)  )THEN
    IF( dodomaindiag .and. dotdwrite )then

      if( sgsmodel.ge.1 .or. sgsmodel.le.4 )then
              ! use tij arrays:
            call getepst(xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,ua ,va ,wa ,  &
                         t11,t12,t13,t22,t23,t33,rf,                         &
                         dum1,dum2,epst)
      elseif( sgsmodel.eq.5 .or. sgsmodel.eq.6 )then
              ! use mij arrays:
            call getepst(xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,ua ,va ,wa ,  &
                         m11,m12,m13,m22,m23,m33,rf,                         &
                         dum1,dum2,epst)
      else
        print *,'  93872 '
        call stopcm1
      endif

      call getepsd(xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,ua ,va ,wa ,  &
                   dum1,dum2,epsd1,epsd2,rho0,rf0,rrf0,rru,rrv,rrw,dt)

    ENDIF
    ENDIF


!-------------------------------------------------------------------
!  dissip rate for DNS:

        IF( cm1setup.eq.3 )THEN

          rcoef = 1.0/viscosity

          if( axisymm.eq.0 )then

            ! Cartesian grid:

            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
            DO j=1,nj
              do k=2,nk
              do i=1,ni
                tmp11=( c1(i,j,k)*t11(i,j,k-1)**2 + c2(i,j,k)*t11(i,j,k)**2 )
                tmp22=( c1(i,j,k)*t22(i,j,k-1)**2 + c2(i,j,k)*t22(i,j,k)**2 )
                tmp33=( c1(i,j,k)*t33(i,j,k-1)**2 + c2(i,j,k)*t33(i,j,k)**2 )
                tmp12=0.25*( c1(i,j,k)*( ( t12(i,j  ,k-1)**2 + t12(i+1,j+1,k-1)**2 )     &
                                       + ( t12(i,j+1,k-1)**2 + t12(i+1,j  ,k-1)**2 ) )   &
                            +c2(i,j,k)*( ( t12(i,j  ,k  )**2 + t12(i+1,j+1,k  )**2 )     &
                                       + ( t12(i,j+1,k  )**2 + t12(i+1,j  ,k  )**2 ) ) ) 
                tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )
                tmp23=0.5*( t23(i,j,k)**2 + t23(i,j+1,k)**2 )
                dissten(i,j,k)= rcoef*( ( 2.0*( tmp33 ) + ( tmp13 + tmp23 )               &
                                         +2.0*( tmp11 + tmp22 ) + tmp12 )/(rf(i,j,k)**2) )
              enddo
              enddo
            ENDDO

          elseif( axisymm.eq.1 )then

            ! axisymmetric grid:

            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,tmp11,tmp33,tmp12,tmp13,tmp23)
            do k=2,nk
            do j=1,nj
            do i=1,ni
              tmp11=( c1(1,1,k)*t11(i,j,k-1)**2 + c2(1,1,k)*t11(i,j,k)**2 )
              tmp33=( c1(1,1,k)*t33(i,j,k-1)**2 + c2(1,1,k)*t33(i,j,k)**2 )
              tmp12=0.5*( c1(1,1,k)*( t12(i,j  ,k-1)**2 + t12(i+1,j  ,k-1)**2 ) &
                         +c2(1,1,k)*( t12(i,j  ,k  )**2 + t12(i+1,j  ,k  )**2 ) ) 
              tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )
              tmp23=      t23(i,j,k)**2
              dissten(i,j,k)= rcoef*( ( 2.0*( tmp33 ) + ( tmp13 + tmp23 )               &
                                       +2.0*( tmp11 ) + tmp12 )/(rf(i,j,k)**2) )
            enddo
            enddo
            enddo

          endif

        ENDIF

      if(timestats.ge.1) time_turb=time_turb+mytime()

!-------------------------------------------------------------------
!  Variables for PBL and/or CuParam

      get_vars:  &
      IF( idopbl .or. cuparam.ge.1 )THEN

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO k=1,nk
          ! store qi in dum8:
          if( nqi.ge.1 )then
            do j=1,nj
            do i=1,ni
              dum8(i,j,k) = qa(i,j,k,nqi)
            enddo
            enddo
          else
            do j=1,nj
            do i=1,ni
              dum8(i,j,k) = 0.0
            enddo
            enddo
          endif
        ENDDO

        ! here, ppten stores qc:

        if( nqc.ge.1 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          DO k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k) = qa(i,j,k,nqc)
          enddo
          enddo
          ENDDO
        else
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          DO k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k) = 0.0
          enddo
          enddo
          ENDDO
        endif

        if( radopt.ge.1 .or. testcase.eq.4 .or. testcase.eq.5 )then
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum9(i,j,k) = thraten(i,j,k)
          enddo
          enddo
          enddo
        else
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum9(i,j,k) = 0.0
          enddo
          enddo
          enddo
        endif

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        DO j=1,nj
          do k=1,nk
          do i=1,ni
            ! ground-relative winds:
            dum1(i,j,k)=0.5*(ugr(i,j,k)+ugr(i+1,j,k))
            dum2(i,j,k)=0.5*(vgr(i,j,k)+vgr(i,j+1,k))
            dum3(i,j,k)=th0(i,j,k)+tha(i,j,k)
            dum7(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
            dum4(i,j,k)=dum3(i,j,k)*dum7(i,j,k)
            dum5(i,j,k) = dz*rmh(i,j,k)
            dum6(i,j,k) = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
          enddo
          enddo
          do i=1,ni
            ! surface:
            dum6(i,j,1) = psfc(i,j)
            ! top of model:
            dum6(i,j,nk+1)= cgt1*prs(i,j,nk)+cgt2*prs(i,j,nk-1)+cgt3*prs(i,j,nk-2)
          enddo
        ENDDO

        ! dum1  = u (grnd-rel) at scalars
        ! dum2  = v (grnd-rel) at scalars
        ! dum3  = th
        ! dum4  = t
        ! dum5  = dz  (at scalar levels)
        ! dum6  = p3di (p8w, ie at w levels)
        ! dum7  = pi
        ! dum8  = qi
        ! dum9  = total radiative theta tendency
        ! ppten = qc
        ! divx  = qv

        if( sfcmodel.eq.4 )then
          do j=1,nj
          do i=1,ni
            wspd(i,j) = swspd(i,j)
          enddo
          enddo
        endif

        if( iice.eq.1 )then
          flag_qi = .true.
        else
          flag_qi = .false.
        endif

        if(timestats.ge.1) time_misc=time_misc+mytime()

      ENDIF  get_vars

!-------------------------------------------------------------------
!  PBL schemes:

    doingpbl:  &
    IF( idopbl .and. ipbl.ge.1 )THEN

    dopbl:  &
    IF( getpbl )THEN

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      pblopt:  &
      if( use_pbl )then

      doysu:  &
      if(ipbl.eq.1)then
        ! YSU PBL:

        ysu_topdown_pblmix = .false.

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k,tmp2d)
        ysu_j_loop:  &
        DO j=1,nj

        tmp2d = 0.0

        do k=1,nk
        do i=1,ni
          tmp2d(i,k, 1) =  dum1(i,j,k)  ! u
          tmp2d(i,k, 2) =  dum2(i,j,k)  ! v
          tmp2d(i,k, 3) =  dum4(i,j,k)  ! T
          tmp2d(i,k, 4) =  divx(i,j,k)  ! qv
          tmp2d(i,k, 5) = ppten(i,j,k)  ! qc
          tmp2d(i,k, 6) =  dum8(i,j,k)  ! qi
          tmp2d(i,k, 7) =   prs(i,j,k)  ! p3d
          tmp2d(i,k, 9) =  dum7(i,j,k)  ! pi3d
          tmp2d(i,k,10) =  dum5(i,j,k)  ! dz8w
          tmp2d(i,k,11) =  dum9(i,j,k)  ! rthraten
        enddo
        enddo

        do k=1,nk+1
        do i=1,ni
          tmp2d(i,k,8) =  dum6(i,j,k)  ! p3di
        enddo
        enddo

        tmp1d = 0.0

        do i=1,ni
          tmp1d(i,1) =  1.0              ! ctopo
          tmp1d(i,2) =  1.0              ! ctopo2
          tmp1d(i,3) =  0.0              ! uoce
          tmp1d(i,4) =  0.0              ! voce
          tmp1d(i,5) =  0.0              ! frcurb
        enddo

      ! cm1r21.1: using ysu from MMM-physics !

        call ysu_wrapper(                                          &
                  ux=tmp2d(ibtmp,kbtmp,1),                        &
                  vx=tmp2d(ibtmp,kbtmp,2),                        &
                  tx=tmp2d(ibtmp,kbtmp,3),                        &
                  qvx=tmp2d(ibtmp,kbtmp,4),                       &
                  qcx=tmp2d(ibtmp,kbtmp,5),                       &
                  qix=tmp2d(ibtmp,kbtmp,6),                       &
                  nmix=1,                                         &
                  qmix=tmp2d(ibtmp,kbtmp,23),                     &
                  p2d=tmp2d(ibtmp,kbtmp,7),                       &
                  p2di=tmp2d(ibtmp,kbtmp,8),                      &
                  pi2d=tmp2d(ibtmp,kbtmp,9),                      &
                  f_qc=.true.,                                    &
                  f_qi=.true.,                                    &
                  utnp=tmp2d(ibtmp,kbtmp,12),                     &
                  vtnp=tmp2d(ibtmp,kbtmp,13),                     &
                  ttnp=tmp2d(ibtmp,kbtmp,14),                     &
                  qvtnp=tmp2d(ibtmp,kbtmp,15),                    &
                  qctnp=tmp2d(ibtmp,kbtmp,16),                    &
                  qitnp=tmp2d(ibtmp,kbtmp,17),                    &
                  qmixtnp=tmp2d(ibtmp,kbtmp,24),                  &
                  cp=cp,                                          &
                  g=g,                                            &
                  rovcp=rovcp,                                    &
                  rd=rd,                                          &
                  rovg=rovg,                                      &
                  ep1=ep1,                                        &
                  ep2=ep2,                                        &
                  karman=karman,                                  &
                  xlv=xlv,                                        &
                  rv=rv,                                          &
                  dz8w2d=tmp2d(ibtmp,kbtmp,10),                   &
                  psfcpa=psfc(1,j),                               &
                  znt=zntmp(1,j),                                 &
                  ust=ust(1,j),                                   &
                  hpbl=hpbl(1,j),                                 &
                  psim=fm(1,j),                                   &
                  psih=fh(1,j),                                   &
                  xland=xland(1,j),                               &
                  hfx=hfx(1,j),                                   &
                  qfx=qfx(1,j),                                   &
                  wspd=wspd(1,j),                                 &
                  br=br(1,j),                                     &
                  dt=dt,                                          &
                  kpbl1d=kpbl2d(1,j),                             &
                  exch_hx=tmp2d(ibtmp,kbtmp,18),                  &
                  exch_mx=tmp2d(ibtmp,kbtmp,19),                  &
                  wstar=wstar(1,j),                               &
                  delta=delta(1,j),                               &
                  u10=u10(1,j),                                   &
                  v10=v10(1,j),                                   &
                  uox=tmp1d(ibtmp,3),                             &
                  vox=tmp1d(ibtmp,4),                             &
                  rthraten=tmp2d(ibtmp,kbtmp,11),                 &
                  ysu_topdown_pblmix=ysu_topdown_pblmix,          &
                  ctopo=tmp1d(ibtmp,1),                           &
                  ctopo2=tmp1d(ibtmp,2),                          &
                  a_u=tmp2d(ibtmp,kbtmp,31), &
                  a_v=tmp2d(ibtmp,kbtmp,32), &
                  a_t=tmp2d(ibtmp,kbtmp,33), &
                  a_q=tmp2d(ibtmp,kbtmp,34), &
                  a_e=tmp2d(ibtmp,kbtmp,35), &
                  b_u=tmp2d(ibtmp,kbtmp,36), &
                  b_v=tmp2d(ibtmp,kbtmp,37), &
                  b_t=tmp2d(ibtmp,kbtmp,38), &
                  b_q=tmp2d(ibtmp,kbtmp,39), &
                  b_e=tmp2d(ibtmp,kbtmp,40), &
                  sfk=tmp2d(ibtmp,kbtmp,41), &
                  vlk=tmp2d(ibtmp,kbtmp,42), &
                  dlu=tmp2d(ibtmp,kbtmp,43), &
                  dlg=tmp2d(ibtmp,kbtmp,44), &
                  frcurb=tmp1d(ibtmp,5), &
                  flag_bep = .false. ,                            &
                  its=1, ite=ni,                                  &
                  kts=1, kte=nk,                                  &
                  errmsg=errmsg ,                                 &
                  errflg=errflg )

        do k=1,nk
        do i=1,ni
          upten(i,j,k) = tmp2d(i,k,12)
          vpten(i,j,k) = tmp2d(i,k,13)
          thpten(i,j,k) = tmp2d(i,k,14)
          qvpten(i,j,k) = tmp2d(i,k,15)
          qcpten(i,j,k) = tmp2d(i,k,16)
          qipten(i,j,k) = tmp2d(i,k,17)
        enddo
        enddo

        IF( dotbud .and. td_diss.ge.1 )THEN
          do k=1,nk
          do i=1,ni
            tdiag(i,j,k,td_diss) = tmp2d(i,k,22)
          enddo
          enddo
        ENDIF

        IF( output_km.eq.1 .or. output_kh.eq.1 .or. dodomaindiag .or. dostat )THEN
          do k=1,nk
          do i=1,ni
            xkzh(i,j,k) = tmp2d(i,k,18)
            xkzm(i,j,k) = tmp2d(i,k,19)
          enddo
          enddo
        ENDIF

        ENDDO  ysu_j_loop

        if( cm1setup.eq.4 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! zero-out diffusivities in LES subdomain:
            if( cm0(i,j).gt.cmemin )then
              xkzh(i,j,k) = 0.0
!!!              xkzq(i,j,k) = 0.0
              xkzm(i,j,k) = 0.0
            endif
          enddo
          enddo
          enddo
        endif

      endif  doysu

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      dobr09:  &
      if(ipbl.eq.2)then
        ! Simple Smagorinsky-like turbulence scheme for vert. direction:
        ! aka, Bryan-Rotunno-09 PBL:

        if( doimpl.eq.1 )then
          !  Arrays for vimpl turbs:
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thten(i,j,k)  = khv(i,j,k  )*mf(i,j,k  )*rf(i,j,k  )*mh(i,j,k)*rr(i,j,k)
            thten1(i,j,k) = khv(i,j,k+1)*mf(i,j,k+1)*rf(i,j,k+1)*mh(i,j,k)*rr(i,j,k)
          enddo
          enddo
          enddo
        endif

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          sadv(i,j,k) = (th0(i,j,k)-th0r)+tha(i,j,k)
        enddo
        enddo
        enddo

        call     turbsz(1,dt,dosfcflx,thflux,mh,mf,  &
                       thpten,dum1,dum2,dum3,rho,rr,rf,sadv,khv,thten,thten1,1)
        if( imoist.eq.1 )  &
        call     turbsz(1,dt,dosfcflx,qvflux,mh,mf,  &
                       qvpten,dum1,dum2,dum3,rho,rr,rf,qa(ib,jb,kb,nqv),khv,thten,thten1,2)
        if( imoist.eq.1 .and. nqc.ge.1 )  &
        call     turbsz(0,dt,dosfcflx,qvflux,mh,mf,  &
                       qcpten,dum1,dum2,dum3,rho,rr,rf,qa(ib,jb,kb,nqc),khv,thten,thten1,0)
        if( imoist.eq.1 .and. nqi.ge.1 )  &
        call     turbsz(0,dt,dosfcflx,qvflux,mh,mf,  &
                       qipten,dum1,dum2,dum3,rho,rr,rf,qa(ib,jb,kb,nqi),khv,thten,thten1,0)
        call     turbuz(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                       upten,dum1,dum2,dum3,thten,thten1,ua,wa,t11,t12,t13,t22,kmv, &
                       kmw,ufw,u1b,u2pt,ufwk)
        call     turbvz(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                       vpten,dum1,dum2,dum3,thten,thten1,va,wa,t12,t22,t23,kmv, &
                       kmw,vfw,v1b,v2pt,vfwk)

      endif  dobr09

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      dogfsedmf:  &
      if( ipbl.eq.3 )then
        ! HWRF PBL  (aka, GFS-EDMF)

        ep1 = rv/rd - 1.0
        ep2 = rd/rv
        rovg = rd/g

        if( sfcmodel.ne.4 )then
          do j=1,nj
          do i=1,ni
            mznt(i,j) = znt(i,j)
          enddo
          enddo
        endif

        call  BL_GFSEDMF(U3D=dum1,V3D=dum2,TH3D=dum3,T3D=dum4,QV3D=divx,            &
                  QC3D=ppten,QI3D=dum8,P3D=prs,PI3D=dum7,                           &
                  RUBLTEN=upten,RVBLTEN=vpten,RTHBLTEN=thpten,                      &
                  RQVBLTEN=qvpten,RQCBLTEN=qcpten,RQIBLTEN=qipten,                  & 
                  CP=cp,G=g,ROVCP=rovcp,R=rd,ROVG=rovg,                             &
                  P_QI=p_qi,P_FIRST_SCALAR=p_first_scalar,                          &
                  dz8w=dum5,z=zf(ib,jb,kb),PSFC=psfc,                               &
                  UST=ust,PBL=hpbl,PSIM=psim,PSIH=psih,                             &
                  HFX=hfx,QFX=qfx,TSK=tsk,GZ1OZ0=gz1oz0,WSPD=wspd,BR=br,            &
                  DT=dt,KPBL2D=kpbl2d,EP1=ep1,KARMAN=karman,                        &
                  DISHEAT=disheat,                                                  &
                  RTHRATEN=dum9,                                                    &
                  HPBL2D=hpbl2d, EVAP2D=evap2d, HEAT2D=heat2d,                      &
                  U10=u10,V10=v10,ZNT=mznt,                                         &
                  DKU3D=xkzm,DKT3D=xkzh,                                            & 
                  VAR_RIC=var_ric,coef_ric_l=coef_ric_l,coef_ric_s=coef_ric_s,      &
                  alpha=gfs_alpha,xland=xland, pert_pbl=pert_pbl,                   &
                  ens_random_seed=ens_random_seed, ens_pblamp=ens_pblamp,           &
                  brcr=brcr,wscale=wscale,wscaleu=wscaleu,prkpp=prkpp,              &
                  stau=stau,dissten=dissten(ib,jb,kb),dheat=thten,                  &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,          &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,            &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk              )
        IF( dotbud .and. td_diss.ge.1 )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dissipative heating rate:
            tdiag(i,j,k,td_diss) = thten(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
        IF( dotbud .and. td_pbl.ge.1 )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! for budget, subtract dissip heating from PBL tendency:
            tdiag(i,j,k,td_pbl) = thpten(i,j,k)-thten(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

      endif  dogfsedmf

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      domynn:  &
      if( ipbl.eq.4 .or. ipbl.eq.5 )then
      if( .not. restarted )then  ! do not call mynn on a restart
        ! MYNN PBL:

        IF( imoist.eq.1 )THEN
          if( nnci.ge.1 )then
            do k=1,nk
            do j=1,nj
            do i=1,ni
              sadv(i,j,k) = qa(i,j,k,nnci)
            enddo
            enddo
            enddo
          else
            do k=1,nk
            do j=1,nj
            do i=1,ni
              sadv(i,j,k) = 0.0
            enddo
            enddo
            enddo
          endif
          if( nncc.ge.1 )then
            do k=1,nk
            do j=1,nj
            do i=1,ni
              thterm(i,j,k) = qa(i,j,k,nncc)
            enddo
            enddo
            enddo
          else
            do k=1,nk
            do j=1,nj
            do i=1,ni
              thterm(i,j,k) = 0.0
            enddo
            enddo
            enddo
          endif
        ELSE
          do k=1,nk
          do j=1,nj
          do i=1,ni
            sadv(i,j,k) = 0.0
            thterm(i,j,k) = 0.0
          enddo
          enddo
          enddo
        ENDIF

        FLAG_QC = .true.
        FLAG_QI = .true.
        if( ptype.eq.3 .or. ptype.eq.5 )then
          FLAG_QNI = .true.
        else
          FLAG_QNI = .false.
        endif
        FLAG_QNC = .false.
        FLAG_QNWFA = .false.
        FLAG_QNIFA = .false.

      mynn_j_loop:  &
      DO j=1,nj

        tmp1d = 0.0
        tmp2d = 0.0

        IF( nstep.eq.0 )THEN
          initflag=1
          if( myid.eq.0 .and. j.eq.1 ) print *,'  initflag = ',initflag
        ELSE
          initflag=0
        ENDIF

        if( restarted )then
          restart = .true.
        else
          restart = .false.
        endif

        cycling = .false.

        ! dum1  = u (grnd-rel) at scalars
        ! dum2  = v (grnd-rel) at scalars
        ! dum3  = th
        ! dum4  = t
        ! dum5  = dz  (at scalar levels)
        ! dum6  = p3di (p8w, ie at w levels)
        ! dum7  = pi
        ! dum8  = qi
        ! dum9  = total radiative theta tendency
        ! ppten = qc
        ! divx  = qv
        ! sadv  = nnci
        ! thterm = nncc

        do k=1,nk
        do i=1,ni
          tmp2d(i,k, 1) = dum5(i,j,k)    ! dz
          tmp2d(i,k, 2) = dum1(i,j,k)    ! u
          tmp2d(i,k, 3) = dum2(i,j,k)    ! v
          tmp2d(i,k, 4) = 0.5*(wa(i,j,k)+wa(i,j,k+1))  ! w
          tmp2d(i,k, 5) = dum3(i,j,k)    ! th
          tmp2d(i,k, 6) = divx(i,j,k)    ! qv
          tmp2d(i,k, 7) = ppten(i,j,k)   ! qc
          tmp2d(i,k, 8) = dum8(i,j,k)    ! qi
          tmp2d(i,k, 9) = 0.0            ! qnc
          tmp2d(i,k,10) = 0.0            ! qni
          tmp2d(i,k,11) = 0.0            ! qnwfa
          tmp2d(i,k,12) = 0.0            ! qnifa
          tmp2d(i,k,13) = 0.0            ! qnbca
          tmp2d(i,k,14) = prs(i,j,k)     ! p
          tmp2d(i,k,15) = dum7(i,j,k)    ! exner
          tmp2d(i,k,16) = rho(i,j,k)     ! rho
          tmp2d(i,k,17) = dum4(i,j,k)    ! t3d
          tmp2d(i,k,18) = qke(i,j,k)
          tmp2d(i,k,19) = qke3d(i,j,k)   ! qke updated after advection
          tmp2d(i,k,20) = sh3d(i,j,k)
          tmp2d(i,k,21) = sm3d(i,j,k)
          tmp2d(i,k,22) = tsq(i,j,k)
          tmp2d(i,k,23) = qsq(i,j,k)
          tmp2d(i,k,24) = cov(i,j,k)
          tmp2d(i,k,25) = 0.0            ! uten
          tmp2d(i,k,26) = 0.0            ! vten
          tmp2d(i,k,27) = 0.0            ! thten
          tmp2d(i,k,28) = 0.0            ! qvten
          tmp2d(i,k,29) = 0.0            ! qcten
          tmp2d(i,k,30) = 0.0            ! qiten
          tmp2d(i,k,31) = 0.0            ! qncten
          tmp2d(i,k,32) = 0.0            ! qniten
          tmp2d(i,k,33) = 0.0            ! qnwfablten
          tmp2d(i,k,34) = 0.0            ! qnifablten
          tmp2d(i,k,35) = 0.0            ! qnbcablten
          tmp2d(i,k,36) = 0.0            ! exch_h
          tmp2d(i,k,37) = 0.0            ! exch_m
          tmp2d(i,k,38) = el_pbl(i,j,k)
          tmp2d(i,k,39) = dqke(i,j,k)
          tmp2d(i,k,40) = qwt(i,j,k)
          tmp2d(i,k,41) = qshear(i,j,k)
          tmp2d(i,k,42) = qbuoy(i,j,k)
          tmp2d(i,k,43) = qdiss(i,j,k)
          tmp2d(i,k,44) = qc_bl(i,j,k)
          tmp2d(i,k,45) = qi_bl(i,j,k)
          tmp2d(i,k,46) = cldfra_bl(i,j,k)
          tmp2d(i,k,47) = edmf_a(i,j,k)
          tmp2d(i,k,48) = edmf_w(i,j,k)
          tmp2d(i,k,49) = edmf_qt(i,j,k)
          tmp2d(i,k,50) = edmf_thl(i,j,k)
          tmp2d(i,k,51) = edmf_ent(i,j,k)
          tmp2d(i,k,52) = edmf_qc(i,j,k)
          tmp2d(i,k,53) = sub_thl3d(i,j,k)
          tmp2d(i,k,54) = sub_sqv3d(i,j,k)
          tmp2d(i,k,55) = det_thl3d(i,j,k)
          tmp2d(i,k,56) = det_sqv3d(i,j,k)
          tmp2d(i,k,57) = dum9(i,j,k)    ! rthraten
          tmp2d(i,k,58) = 0.0            ! pattern_spp_pbl
        enddo
        enddo

        do i=1,ni
          tmp1d(i, 1) = 0.0    ! uoce
          tmp1d(i, 2) = 0.0    ! voce
        enddo

      call mynnedmf_wrapper_run(                 &
        initflag=initflag,  &
        restart=restart,  &
        cycling=cycling,                &
        delt=dt,  &
        dz=tmp2d(ibtmp,kbtmp, 1),  &
        dxc=dx,  &
        znt=znt(1,j),                         &
        u=tmp2d(ibtmp,kbtmp, 2),  &
        v=tmp2d(ibtmp,kbtmp, 3),  &
        w=tmp2d(ibtmp,kbtmp, 4),  &
        th=tmp2d(ibtmp,kbtmp, 5),                                &
        qv=tmp2d(ibtmp,kbtmp, 6),  &
        qc=tmp2d(ibtmp,kbtmp, 7),  &
        qi=tmp2d(ibtmp,kbtmp, 8),  &
        qnc=tmp2d(ibtmp,kbtmp, 9),  &
        qni=tmp2d(ibtmp,kbtmp,10),  &
        qnwfa=tmp2d(ibtmp,kbtmp,11),  &
        qnifa=tmp2d(ibtmp,kbtmp,12),  &
        qnbca=tmp2d(ibtmp,kbtmp,13),      &
!     &  ozone,                                   &
        p=tmp2d(ibtmp,kbtmp,14),  &
        exner=tmp2d(ibtmp,kbtmp,15),  &
        rho=tmp2d(ibtmp,kbtmp,16),  &
        t3d=tmp2d(ibtmp,kbtmp,17),                         &
        xland=xland(1,j),  &
        ts=tsk(1,j),  &
        qsfc=qsfc(1,j),  &
        ps=psfc(1,j),                        &
        ust=ust(1,j),  &
        ch=ch(1,j),  &
        hfx=hfx(1,j),  &
        qfx=qfx(1,j),  &
        rmol=rmol(1,j),  &
        wspd=wspd(1,j),                &
        uoce=tmp1d(1,1),  &
        voce=tmp1d(1,2),                               &
        qke=tmp2d(ibtmp,kbtmp,18),  &
        qke_adv=tmp2d(ibtmp,kbtmp,19),  &
        sh3d=tmp2d(ibtmp,kbtmp,20),  &
        sm3d=tmp2d(ibtmp,kbtmp,21),                   &
!--- chem/smoke
!#if (WRF_CHEM == 1)
!        mix_chem,chem3d,vd3d,nchem,kdvel,        &
!        ndvel,num_vert_mix,                      &
!!     &  frp_mean,emis_ant_no,enh_mix,            & !to be included soon
!#endif
!--- end chem/smoke 
        Tsq=tmp2d(ibtmp,kbtmp,22),  &
        Qsq=tmp2d(ibtmp,kbtmp,23),  &
        Cov=tmp2d(ibtmp,kbtmp,24),                             &
        rublten=tmp2d(ibtmp,kbtmp,25),  &
        rvblten=tmp2d(ibtmp,kbtmp,26),  &
        rthblten=tmp2d(ibtmp,kbtmp,27),                &
        rqvblten=tmp2d(ibtmp,kbtmp,28),  &
        rqcblten=tmp2d(ibtmp,kbtmp,29),  &
        rqiblten=tmp2d(ibtmp,kbtmp,30),              &
        rqncblten=tmp2d(ibtmp,kbtmp,31),  &
        rqniblten=tmp2d(ibtmp,kbtmp,32),                     &
        rqnwfablten=tmp2d(ibtmp,kbtmp,33),  &
        rqnifablten=tmp2d(ibtmp,kbtmp,34),  &
        rqnbcablten=tmp2d(ibtmp,kbtmp,35),     &
!     &  ro3blten,                                &
        exch_h=tmp2d(ibtmp,kbtmp,36),  &
        exch_m=tmp2d(ibtmp,kbtmp,37),  &
        pblh=hpbl(1,j),  &
        kpbl=kpbl2d(1,j),  &
        el_pbl=tmp2d(ibtmp,kbtmp,38),          &
        dqke=tmp2d(ibtmp,kbtmp,39),  &
        qwt=tmp2d(ibtmp,kbtmp,40),  &
        qshear=tmp2d(ibtmp,kbtmp,41),  &
        qbuoy=tmp2d(ibtmp,kbtmp,42),  &
        qdiss=tmp2d(ibtmp,kbtmp,43),             &
        qc_bl=tmp2d(ibtmp,kbtmp,44),  &
        qi_bl=tmp2d(ibtmp,kbtmp,45),  &
        cldfra_bl=tmp2d(ibtmp,kbtmp,46),                   &
        edmf_a=tmp2d(ibtmp,kbtmp,47),  &
        edmf_w=tmp2d(ibtmp,kbtmp,48),  &
        edmf_qt=tmp2d(ibtmp,kbtmp,49),                   &
        edmf_thl=tmp2d(ibtmp,kbtmp,50),  &
        edmf_ent=tmp2d(ibtmp,kbtmp,51),  &
        edmf_qc=tmp2d(ibtmp,kbtmp,52),               &
        sub_thl3d=tmp2d(ibtmp,kbtmp,53),  &
        sub_sqv3d=tmp2d(ibtmp,kbtmp,54),                     &
        det_thl3d=tmp2d(ibtmp,kbtmp,55),  &
        det_sqv3d=tmp2d(ibtmp,kbtmp,56),                     &
        nupdraft=nupdraft(1,j),  &
        maxMF=maxmf(1,j),  &
        ktop_plume=ktop_plume(1,j),               &
        rthraten=tmp2d(ibtmp,kbtmp,57),                                &
        dheat3d=tmp2d(ibtmp,kbtmp,59),         &
        tke_budget=bl_mynn_tkebudget,         &
        bl_mynn_tkeadvect=bl_mynn_tkeadvect,   &
        bl_mynn_cloudpdf=bl_mynn_cloudpdf,   &
        bl_mynn_mixlength=bl_mynn_mixlength,   &
        icloud_bl=icloud_bl,          &
        bl_mynn_edmf=bl_mynn_edmf,        &
        bl_mynn_edmf_mom=bl_mynn_edmf_mom,   &
        bl_mynn_edmf_tke=bl_mynn_edmf_tke,    &
        bl_mynn_cloudmix=bl_mynn_cloudmix,   &
        bl_mynn_mixqt=bl_mynn_mixqt,       &
        bl_mynn_output=bl_mynn_output,     &
        bl_mynn_closure=bl_mynn_closure,     &
        bl_mynn_mixscalars=bl_mynn_mixscalars,                      &
        spp_pbl=spp_pbl,  &
        pattern_spp_pbl=tmp2d(ibtmp,kbtmp,58),                 &
        flag_qc=flag_qc,  &
        flag_qi=flag_qi,                         &
        flag_qnc=flag_qnc,  &
        flag_qni=flag_qni,                       &
        flag_qnwfa=.false.,  &
        flag_qnifa=.false.,  &
        flag_qnbca=.false.,        &
        ids=1,ide=ni+1,jds=1,jde=2,kds=1,kde=nk+1,   &
        ims=1,ime=ni  ,jms=1,jme=1,kms=1,kme=nk  ,   &
        its=1,ite=ni  ,jts=1,jte=1,kts=1,kte=nk      )

        do k=1,nk
        do i=1,ni
          qke(i,j,k)       = tmp2d(i,k,18)
          qke_adv(i,j,k)   = tmp2d(i,k,19)
          sh3d(i,j,k)      = tmp2d(i,k,20)
          sm3d(i,j,k)      = tmp2d(i,k,21)
          tsq(i,j,k)       = tmp2d(i,k,22)
          qsq(i,j,k)       = tmp2d(i,k,23)
          cov(i,j,k)       = tmp2d(i,k,24)
          upten(i,j,k)     = tmp2d(i,k,25)
          vpten(i,j,k)     = tmp2d(i,k,26)
          thpten(i,j,k)    = tmp2d(i,k,27)
          qvpten(i,j,k)    = tmp2d(i,k,28)
          qcpten(i,j,k)    = tmp2d(i,k,29)
          qipten(i,j,k)    = tmp2d(i,k,30)
          xkzh(i,j,k)      = tmp2d(i,k,36)
          xkzm(i,j,k)      = tmp2d(i,k,37)
          el_pbl(i,j,k)    = tmp2d(i,k,38)
          dqke(i,j,k)      = tmp2d(i,k,39)
          qwt(i,j,k)       = tmp2d(i,k,40)
          qshear(i,j,k)    = tmp2d(i,k,41)
          qbuoy(i,j,k)     = tmp2d(i,k,42)
          qdiss(i,j,k)     = tmp2d(i,k,43)
          qc_bl(i,j,k)     = tmp2d(i,k,44)
          qi_bl(i,j,k)     = tmp2d(i,k,45)
          cldfra_bl(i,j,k) = tmp2d(i,k,46)
          edmf_a(i,j,k)    = tmp2d(i,k,47)
          edmf_w(i,j,k)    = tmp2d(i,k,48)
          edmf_qt(i,j,k)   = tmp2d(i,k,49)
          edmf_thl(i,j,k)  = tmp2d(i,k,50)
          edmf_ent(i,j,k)  = tmp2d(i,k,51)
          edmf_qc(i,j,k)   = tmp2d(i,k,52)
          sub_thl3d(i,j,k) = tmp2d(i,k,53)
          sub_sqv3d(i,j,k) = tmp2d(i,k,54)
          det_thl3d(i,j,k) = tmp2d(i,k,55)
          det_sqv3d(i,j,k) = tmp2d(i,k,56)
        enddo
        enddo

      ENDDO  mynn_j_loop

        if( cm1setup.eq.4 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! zero-out qke in LES subdomain:
            if( cm0(i,j).gt.cmemin )then
              qke(i,j,k) = 0.0
              xkzh(i,j,k) = 0.0
              xkzm(i,j,k) = 0.0
            endif
          enddo
          enddo
          enddo
        endif

      ! Diagnostics !
        IF( dotbud .and. td_diss.ge.1 )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dissipative heating rate:
            tdiag(i,j,k,td_diss) = tmp2d(i,k,59)
          enddo
          enddo
          enddo
        ENDIF
        IF( dotbud .and. td_pbl.ge.1 )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! for budget, subtract dissip heating from PBL tendency:
            tdiag(i,j,k,td_pbl) = thpten(i,j,k)-tmp2d(i,k,59)
          enddo
          enddo
          enddo
        ENDIF
      ! Diagnostics !

      else
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  ... skipping mynn ... '
        if(myid.eq.0) print *
      endif
      endif  domynn

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      domyj:  &
      if( ipbl.eq.6 )then
        ! MYJ PBL:

        ! dum1 = u at scalars
        ! dum2 = v at scalars
        ! dum3 = th
        ! dum4 = t
        ! dum5 = dz
        ! dum6 = p3di (p8w)
        ! dum7 = pi

        ! divx  = qv
        ! ppten = qc
        ! dum8  = qi

        ! thten1 = exch_h
        ! thten  = exch_m


      myjloop:  &
      do j=1,nj
        tmp2d = 0.0
        do k=1,nk
        do i=1,ni
          tmp2d(i,k, 1) = 0.0
          tmp2d(i,k, 2) = 0.0
          tmp2d(i,k, 3) = 0.0
          tmp2d(i,k, 4) = tke_myj(i,j,k)
          tmp2d(i,k, 5) = el_myj(i,j,k)
          tmp2d(i,k, 6) = dum5(i,j,k)
          tmp2d(i,k, 7) = prs(i,j,k)
          tmp2d(i,k, 9) = dum3(i,j,k)
          tmp2d(i,k,10) = dum4(i,j,k)
          tmp2d(i,k,11) = dum7(i,j,k)
          tmp2d(i,k,12) = divx(i,j,k)
          tmp2d(i,k,13) = ppten(i,j,k)
          tmp2d(i,k,14) = dum8(i,j,k)
          tmp2d(i,k,15) = dum1(i,j,k)
          tmp2d(i,k,16) = dum2(i,j,k)
          tmp2d(i,k,17) = rho(i,j,k)
          tmp2d(i,k,19) = upten(i,j,k)
          tmp2d(i,k,20) = vpten(i,j,k)
          tmp2d(i,k,21) = thpten(i,j,k)
          tmp2d(i,k,22) = qvpten(i,j,k)
          tmp2d(i,k,23) = qcpten(i,j,k)
          tmp2d(i,k,24) = qipten(i,j,k)
          tmp2d(i,k,25) = qspten(i,j,k)
          tmp2d(i,k,26) = qrpten(i,j,k)
          tmp2d(i,k,27) = qgpten(i,j,k)
        enddo
        enddo
        do k=1,nk+1
        do i=1,ni
          tmp2d(i,k, 8) = dum6(i,j,k)
          tmp2d(i,k,18) = 0.0
          tmp2d(i,k,28) = 0.0
        enddo
        enddo
        call     MYJPBL(dt=DT,STEPBL=1,HT=zs(1,j),DZ=tmp2d(ibtmp,kbtmp,6)                     &
                       ,PMID=tmp2d(ibtmp,kbtmp,7),PINT=tmp2d(ibtmp,kbtmp,8),TH=tmp2d(ibtmp,kbtmp,9),T=tmp2d(ibtmp,kbtmp,10),EXNER=tmp2d(ibtmp,kbtmp,11),QV=tmp2d(ibtmp,kbtmp,12),QCW=tmp2d(ibtmp,kbtmp,13),QCI=tmp2d(ibtmp,kbtmp,14),QCS=tmp2d(ibtmp,kbtmp,1),QCR=tmp2d(ibtmp,kbtmp,2),QCG=tmp2d(ibtmp,kbtmp,3)    &
                       ,U=tmp2d(ibtmp,kbtmp,15),V=tmp2d(ibtmp,kbtmp,16),RHO=tmp2d(ibtmp,kbtmp,17),TSK=tsk(1,j),QSFC=qsfc(1,j),CHKLOWQ=chklowq(1,j),THZ0=thz0(1,j),QZ0=qz0(1,j),UZ0=uz0(1,j),VZ0=vz0(1,j)      &
                       ,LOWLYR=lowlyr(1,j),XLAND=xland(1,j),SICE=sice(1,j),SNOW=snow(1,j)                         &
                       ,TKE_MYJ=tmp2d(ibtmp,kbtmp,4),EXCH_H=tmp2d(ibtmp,kbtmp,18),EXCH_M=tmp2d(ibtmp,kbtmp,28),USTAR=ust(1,j),ZNT=znt(1,j),EL_MYJ=tmp2d(ibtmp,kbtmp,5),PBLH=hpbl(1,j),KPBL=kpbl2d(1,j),CT=ct(1,j)   &
                       ,AKHS=akhs(1,j),AKMS=akms(1,j),ELFLX=lh(1,j),MIXHT=mixht(1,j)                          &
                       ,RUBLTEN=tmp2d(ibtmp,kbtmp,19),RVBLTEN=tmp2d(ibtmp,kbtmp,20),RTHBLTEN=tmp2d(ibtmp,kbtmp,21),RQVBLTEN=tmp2d(ibtmp,kbtmp,22),RQCBLTEN=tmp2d(ibtmp,kbtmp,23)     &
                       ,RQIBLTEN=tmp2d(ibtmp,kbtmp,24),RQSBLTEN=tmp2d(ibtmp,kbtmp,25),RQRBLTEN=tmp2d(ibtmp,kbtmp,26),RQGBLTEN=tmp2d(ibtmp,kbtmp,27),           &
                  ids=1  ,ide=ni+1 , jds=1  ,jde=2    , kds=1  ,kde=nk+1 ,          &
                  ims=1 ,ime=ni      , jms=1  ,jme=1    , kms=1, kme=nk+1  ,            &
                  its=1  ,ite=ni   , jts=1  ,jte=1    , kts=1  ,kte=nk )
        do k=1,nk
        do i=1,ni
          tke_myj(i,j,k) = tmp2d(i,k,4)
          el_myj(i,j,k) = tmp2d(i,k,5)
          xkzh(i,j,k+1) = tmp2d(i,k,18)
          xkzm(i,j,k+1) = tmp2d(i,k,28)
          upten(i,j,k) = tmp2d(i,k,19)
          vpten(i,j,k) = tmp2d(i,k,20)
          thpten(i,j,k) = tmp2d(i,k,21)
          qvpten(i,j,k) = tmp2d(i,k,22)
          qcpten(i,j,k) = tmp2d(i,k,23)
          qipten(i,j,k) = tmp2d(i,k,24)
          qspten(i,j,k) = tmp2d(i,k,25)
          qrpten(i,j,k) = tmp2d(i,k,26)
          qgpten(i,j,k) = tmp2d(i,k,27)
        enddo
        enddo
      enddo  myjloop

        if( cm1setup.eq.4 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! zero-out tke in LES subdomain:
            if( cm0(i,j).gt.cmemin )then
              tke_myj(i,j,k) = 0.0
              xkzh(i,j,k) = 0.0
              xkzm(i,j,k) = 0.0
            endif
          enddo
          enddo
          enddo
        endif

      endif  domyj

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      if( cm1setup.eq.4 )then
        ! zero-out PBL tendencies in LES subdomain:

        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          if( cm0(i,j).gt.cmemin )then
             thpten(i,j,k) = 0.0
             qvpten(i,j,k) = 0.0
             qcpten(i,j,k) = 0.0
             qipten(i,j,k) = 0.0
              upten(i,j,k) = 0.0
              vpten(i,j,k) = 0.0
            qnipten(i,j,k) = 0.0
            qncpten(i,j,k) = 0.0
             qrpten(i,j,k) = 0.0
             qspten(i,j,k) = 0.0
             qgpten(i,j,k) = 0.0
          endif
        enddo
        enddo
        enddo
      endif

      endif  pblopt

    ELSE

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ... skipping pbl calcs ... '
      if(dowr) write(outfile,*)

    ENDIF  dopbl

        if(timestats.ge.1) time_pbl=time_pbl+mytime()


      ysudiss:  &
      if( ipbl.eq.1 )then
        ! Dissipative heating from ysu scheme:

      IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,tem1,tem2)
        DO j=1,nj

          do i=1,ni
            ! assume t13,t23 are zero at top of domain:
            !   thten = t13
            !   thten1 = t23
            thten(i,j,nk+1) = 0.0
            thten1(i,j,nk+1) = 0.0
          enddo

          do k=nk,2,-1
          do i=1,ni
            tem1 = dz*rmh(i,j,k)
            tem2 = rdz*mf(i,j,k)
            ! get stresses from u,v tendencies:
            thten(i,j,k) = thten(i,j,k+1)-upten(i,j,k)*tem1
            thten1(i,j,k) = thten1(i,j,k+1)-vpten(i,j,k)*tem1
            ! NOTE:  dissten is defined at w points:
            !  dum1 = u
            !  dum2 = v
            dissten(i,j,k) = dissten(i,j,k)                                 &
                          +( thten(i,j,k)*( dum1(i,j,k)-dum1(i,j,k-1) )*tem2 &
                            +thten1(i,j,k)*( dum2(i,j,k)-dum2(i,j,k-1) )*tem2 )
          enddo
          enddo

        ENDDO

      ENDIF
      endif  ysudiss

        if(timestats.ge.1) time_pbl=time_pbl+mytime()

      !c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c!

      if( use_pbl )then

        call bcs(upten)
#ifdef MPI
        call comm_1s_start(upten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
#endif
        call bcs(vpten)
#ifdef MPI
        call comm_1s_start(vpten,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_p)
#endif
        IF( axisymm.eq.1 )THEN
          ! cm1r19 bug fix:
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
            upten(0,1,k) = -upten(1,1,k)
            upten(ni+1,1,k) = upten(ni,1,k)
            do i=0,ni+1
              vpten(i,0,k) = vpten(i,1,k)
              vpten(i,2,k) = vpten(i,1,k)
              vpten(i,3,k) = vpten(i,1,k)
            enddo
          enddo
        ENDIF
#ifdef MPI
        call comm_1s_end(upten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
        call comm_1s_end(vpten,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_p)
#endif
      endif

    ENDIF  doingpbl

!-----------------------------------------------------------------------------
!  Convective parameterization:

    usingcp1:  &
    IF( cuparam.eq.1 )THEN
    if( .not. restarted )then  ! do not call tiedtke cu param on a restart

      ! new Tiedtke cu param (from WRFv4.4)

      rdt = 1.0/dt

      ntiedtke_jloop:  &
      DO j=1,nj

        if( use_pbl .and. ipbl.ge.1 )then
          ! add pbl tendency to advective tendency:
          do k=1,nk
          do i=1,ni
            qvften(i,j,k)=qvften(i,j,k)+qvpten(i,j,k)
            thften(i,j,k)=thften(i,j,k)+thpten(i,j,k)
          enddo
          enddo
        endif

        if( radopt.ge.1 .or. testcase.eq.4 .or. testcase.eq.5 )then
          ! add radiative tendency to advective tendency:
          do k=1,nk
          do i=1,ni
            thften(i,j,k)=thften(i,j,k)+dum9(i,j,k)
          enddo
          enddo
        endif

    ! from wrf module_cumulus_driver:
!           RTHFTEN(i,k,j)=(RTHFTEN(i,k,j)+RTHRATEN(i,k,j) &
!                               +RTHBLTEN(i,k,j))*pi(i,k,j)     ! GHB: so, appears to be a temp. tend.
!           RQVFTEN(i,k,j)=RQVFTEN(i,k,j)+RQVBLTEN(i,k,j)

        do k=1,nk
          kk = nk+1-k
          do i=1,ni
            tmp2d(i,kk, 1) = 0.5*(ugr(i,j,k)+ugr(i+1,j,k))  ! u3d
            tmp2d(i,kk, 2) = 0.5*(vgr(i,j,k)+vgr(i,j+1,k))  ! v3d
            tmp2d(i,kk, 3) = dum4(i,j,k)                    ! t3d
            tmp2d(i,kk, 4) = divx(i,j,k)                    ! qv3d
            tmp2d(i,kk, 5) = ppten(i,j,k)                   ! qc3d
            tmp2d(i,kk, 6) = dum8(i,j,k)                    ! qi3d
            tmp2d(i,kk, 7) = qvften(i,j,k)                  ! qvften (total advective + PBL moisture tendency)
            tmp2d(i,kk, 8) = thften(i,j,k)*dum7(i,j,k)      ! thften (total advective + PBL + radiative temperature tendency)
            tmp2d(i,kk, 9) = zh(i,j,k)                      ! 
            tmp2d(i,kk,11) = -0.5*g*rho(i,j,k)*(wa(i,j,k)+wa(i,j,k+1))
            tmp2d(i,kk,12) = prs(i,j,k)
            tmp2d(i,kk,15) = tmp2d(i,kk,1)
            tmp2d(i,kk,16) = tmp2d(i,kk,2)
            tmp2d(i,kk,17) = tmp2d(i,kk,3)
            tmp2d(i,kk,18) = tmp2d(i,kk,4)
            tmp2d(i,kk,19) = tmp2d(i,kk,5)
            tmp2d(i,kk,20) = tmp2d(i,kk,6)
          enddo
        enddo

        do k=1,nk+1
          kk = nk+2-k
          do i=1,ni
            tmp2d(i,kk,10) = zf(i,j,k)
            tmp2d(i,kk,13) = dum6(i,j,k)
          enddo
        enddo

        do i=1,ni
          int1d(i,1) = (abs(xland(i,j)-2.))
          tmp1d(i,1) = 0.0
        enddo

        IF( axisymm.eq.1 .or. ny.eq.1 )THEN
          do i=1,ni
            tmp1d(i,2) = dx*ruh(i)
          enddo
        ELSEIF( nx.eq.1 )THEN
          do i=1,ni
            tmp1d(i,2) = dy*rvh(j)
          enddo
        ELSE
          do i=1,ni
            tmp1d(i,2) = max( dx*ruh(i) , dy*rvh(j) )
          enddo
        ENDIF

        call     cu_ntiedtke_wrapper(pu=tmp2d(ibtmp,kbtmp,1),  &  ! u
                                 pv=tmp2d(ibtmp,kbtmp,2),  &  ! v
                                 pt=tmp2d(ibtmp,kbtmp,3),  &  ! T
                                 pqv=tmp2d(ibtmp,kbtmp,4),  &  ! qv
                                 pqc=tmp2d(ibtmp,kbtmp,5),  &  ! qc
                                 pqi=tmp2d(ibtmp,kbtmp,6),  &  ! qi
                                 pqvf=tmp2d(ibtmp,kbtmp,7),  &  ! qv tendency due to dynamics (kg/kg/s)
                                 ptf=tmp2d(ibtmp,kbtmp,8),  &  ! T tendency due to dynamics (K/s)
                                 poz=tmp2d(ibtmp,kbtmp,9),  &  ! 
                                 pzz=tmp2d(ibtmp,kbtmp,10),  &  ! 
                                 pomg=tmp2d(ibtmp,kbtmp,11), &  ! lagrangian_tendency_of_air_pressure (Pa/s)
                                 pap=tmp2d(ibtmp,kbtmp,12),  &
                                 paph=tmp2d(ibtmp,kbtmp,13),  &
                                 evap=qfx(1,j),  &  ! kinematic surface upward latent heat flux (kg/kg m/s)
                                 hfx=hfx(1,j),  &  ! kinematic surface upward sensible heat flux (K m/s)
                                 zprecc=tmp1d(ibtmp,1),  &  ! deep convective rainfall amount on physics timestep (m)
                                 lndj=int1d(ibtmp,1),  &
                                 lq=ni,  &
                                 km=nk,  &
                                 km1=nkp1,  &
                                 dt=dt,  &
                                 dx=tmp1d(ibtmp,2),  &
                                 errmsg=errmsg,  &
                                 errflg=errflg)

        do k=1,nk
          kk = nk+1-k
          do i=1,ni
             ucuten(i,j,k) = (tmp2d(i,kk,1)-tmp2d(i,kk,15))*rdt
             vcuten(i,j,k) = (tmp2d(i,kk,2)-tmp2d(i,kk,16))*rdt
            thcuten(i,j,k) = (tmp2d(i,kk,3)-tmp2d(i,kk,17))*rdt
            qvcuten(i,j,k) = (tmp2d(i,kk,4)-tmp2d(i,kk,18))*rdt
            qccuten(i,j,k) = (tmp2d(i,kk,5)-tmp2d(i,kk,19))*rdt
            qicuten(i,j,k) = (tmp2d(i,kk,6)-tmp2d(i,kk,20))*rdt
          enddo
        enddo

        IF( cm1setup.eq.4 )THEN
          ! do not apply cuparam in fine-mesh part of domain:
          do n=1,nrain
          do i=1,ni
            if( abs( tmp1d(i,2) - sqrt( (dx_inner)**2 + (dy_inner)**2 ) ) .gt. 1.0 )then
              raincv(i,j,n) = raincv(i,j,n) + tmp1d(i,1)
            endif
          enddo
          enddo
          ! zero-out Cu param tendencies inside fine mesh:
          do k=1,nk
          do i=1,ni
            if( abs( tmp1d(i,2) - sqrt( (dx_inner)**2 + (dy_inner)**2 ) ) .le. 1.0 )then 
              pratecv(i,j) = 0.0
              thcuten(i,j,k) = 0.0
              qvcuten(i,j,k) = 0.0
              qccuten(i,j,k) = 0.0
              qicuten(i,j,k) = 0.0
               ucuten(i,j,k) = 0.0
               vcuten(i,j,k) = 0.0
            endif
          enddo
          enddo
        ELSE
          do n=1,nrain
          do i=1,ni
            raincv(i,j,n) = raincv(i,j,n) + tmp1d(i,1)
          enddo
          enddo
        ENDIF


      ENDDO  ntiedtke_jloop

      if(timestats.ge.1) time_cuprm=time_cuprm+mytime()

      call bcs(ucuten)
      call bcs(vcuten)
#ifdef MPI
      call comm_1s_start(ucuten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
      call comm_1s_start(vcuten,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_p)
      call comm_1s_end(ucuten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
      call comm_1s_end(vcuten,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_p)
#endif

    else
      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  ... skipping tiedtke ... '
      if(myid.eq.0) print *
    endif
    ENDIF  usingcp1

!-------------------------------------------------------------------



      end subroutine sfc_and_turb

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine tkekm(nstep,rtime,dt,ruh,rvh,rmh,zf,mf,rmf,zntmp,ust,ustt,rf, &
                         nm,defv,defh,tk  ,dum2,lenh,grdscl,rgrdscl,           &
                         kmh,kmv,khh,khv,cme,csm,ce1,ce2,tkea,lenscl,dissten,out3d, &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                   &
                         khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,   &
                         khdw1,khdw2,khde1,khde2,khds1,khds2,khdn1,khdn2,   &
                         kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2,   &
                         kvdw1,kvdw2,kvde1,kvde2,kvds1,kvds2,kvdn1,kvdn2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

      integer, intent(in) :: nstep
      real, intent(in) :: rtime,dt
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, intent(in),dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(in), dimension(ib:ie,jb:je) :: zntmp,ust,ustt
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tk,dum2,lenh,grdscl,rgrdscl
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(in),    dimension(ibc:iec,jbc:jec,kbc:kec) :: cme,csm,ce1,ce2
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: lenscl,dissten
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: khdw1,khdw2,khde1,khde2
      real, intent(inout), dimension(imp,kmt) :: khds1,khds2,khdn1,khdn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2
      real, intent(inout), dimension(jmp,kmt) :: kvdw1,kvdw2,kvde1,kvde2
      real, intent(inout), dimension(imp,kmt) :: kvds1,kvds2,kvdn1,kvdn2

!----------------------------------------

      integer :: i,j,k
      real :: prinv,tem1,tem2

      real, parameter :: tke_min         =  1.0e-10
      real, parameter :: nm_min          =  1.0e-6
      real, parameter :: small_len_frac  =  0.001

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
      integer reqs_khd(8)
      integer reqs_kvd(8)
#endif

!------------------------------------------------------------------
!  Get length scales:

    !  get grid scale
    IF(tconfig.eq.1)THEN

      ! single length scale:  appropriate if dx,dy are nearly the same as dz

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,prinv)
      do k=2,nk
!!!        rcs = 1.0/max(1.0e-10,csz(k))
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=( ((dx*ruh(i))*(dy*rvh(j)))*(dz*rmf(i,j,k)) )**0.33333333
        ! cm1r17:  wall condition near surface
        ! cm1r20.1: revisit near-sfc grid scale at a later date
!!!        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                  &
!!!                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j))*rcs)**2)  &
!!!                               ) )
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)

      ! Get turbulence length scale

        if( tkea(i,j,k).le.tke_min )then
          ! 170718:
          tk(i,j,k) = tke_min
          lenscl(i,j,k) = small_len_frac*grdscl(i,j,k)
        else
          tk(i,j,k)=tkea(i,j,k)
          lenscl(i,j,k)=grdscl(i,j,k)
          if(nm(i,j,k).gt.nm_min)then
            lenscl(i,j,k)=c_l*sqrt(tk(i,j,k)/nm(i,j,k))
            lenscl(i,j,k)=min(lenscl(i,j,k),grdscl(i,j,k))
            lenscl(i,j,k)=max(lenscl(i,j,k),small_len_frac*grdscl(i,j,k))
          endif 
        endif

    !  Dissipation

        dissten(i,j,k) = dissten(i,j,k)                           &
                     +(ce1(i,j,k)+ce2(i,j,k)*lenscl(i,j,k)*rgrdscl(i,j,k))    &
                     *tk(i,j,k)*sqrt(tk(i,j,k))/lenscl(i,j,k)

    !  Get km, kh

        kmh(i,j,k)=cme(i,j,k)*sqrt(tk(i,j,k))*lenscl(i,j,k)
        kmv(i,j,k)=kmh(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.nm_min)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=khh(i,j,k)

      enddo
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

      ! two length scales:  one for horizontal, one for vertical

      !$omp parallel do default(shared)   &
      !$omp private(i,j,k,prinv)
      do k=2,nk
!!!        rcs = 1.0/max(1.0e-10,csz(k))
      do j=1,nj
      do i=1,ni

        lenh(i,j,k)=sqrt( (dx*ruh(i))*(dy*rvh(j)) )
        grdscl(i,j,k)=dz*rmf(i,j,k)
        ! cm1r17:  wall condition near surface
        ! cm1r20.1: revisit near-sfc grid scale at a later date
!!!        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                  &
!!!                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j))*rcs)**2)  &
!!!                               ) )
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)

      ! Get turbulence length scale

        if( tkea(i,j,k).le.tke_min )then
          ! 170718:
          tk(i,j,k) = tke_min
          lenscl(i,j,k) = small_len_frac*grdscl(i,j,k)
        else
          tk(i,j,k)=tkea(i,j,k)
          lenscl(i,j,k)=grdscl(i,j,k)
          if(nm(i,j,k).gt.nm_min)then
            lenscl(i,j,k)=c_l*sqrt(tk(i,j,k)/nm(i,j,k))
            lenscl(i,j,k)=min(lenscl(i,j,k),grdscl(i,j,k))
            lenscl(i,j,k)=max(lenscl(i,j,k),small_len_frac*grdscl(i,j,k))
          endif 
        endif

    !  Dissipation

        dissten(i,j,k) = dissten(i,j,k)                           &
                     +(ce1(i,j,k)+ce2(i,j,k)*lenscl(i,j,k)*rgrdscl(i,j,k))    &
                     *tk(i,j,k)*sqrt(tk(i,j,k))/lenscl(i,j,k)

    !  Get km, kh

        kmh(i,j,k)=cme(i,j,k)*sqrt(tk(i,j,k))*lenh(i,j,k)
        kmv(i,j,k)=cme(i,j,k)*sqrt(tk(i,j,k))*lenscl(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.nm_min)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=kmv(i,j,k)*prinv

      enddo
      enddo
      enddo

    ENDIF


    if( nstep.eq.0 .and. myid.eq.0 )then
      print *
      print *,'  zf,grdscl:'
      i = 1
      j = 1
      do k=2,nk
        print *,k,(zf(i,j,k)-zf(i,j,1)),grdscl(i,j,k)
      enddo
    endif

      if(timestats.ge.1) time_turb=time_turb+mytime()

!------------------------------------------------------------
! Set values at boundaries, start comms:

      call bcw(kmh,1)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
#endif

      call bcw(kmv,1)
#ifdef MPI
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
#endif

      call bcw(khh,1)
#ifdef MPI
      call comm_1t_start(khh,khdw1,khdw2,khde1,khde2,   &
                             khds1,khds2,khdn1,khdn2,reqs_khd)
#endif

      call bcw(khv,1)
#ifdef MPI
      call comm_1t_start(khv,kvdw1,kvdw2,kvde1,kvde2,   &
                             kvds1,kvds2,kvdn1,kvdn2,reqs_kvd)
#endif

!--------------------------------------------------------------
!  Finish comms:
#ifdef MPI
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call bct2(kmh)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call bct2(kmv)
      call comm_1t_end(khh,khdw1,khdw2,khde1,khde2,   &
                           khds1,khds2,khdn1,khdn2,reqs_khd)
      call comm_1t_end(khv,kvdw1,kvdw2,kvde1,kvde2,   &
                           kvds1,kvds2,kvdn1,kvdn2,reqs_kvd)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
#endif

!--------------------------------------------------------------
!  finished
      
      return
      end subroutine tkekm


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine tkebc(tkea,kmh,kmv,khh,khv,zntmp,ust,ustt,              &
                       khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,  &
                       kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(in), dimension(ib:ie,jb:je) :: zntmp,ust,ustt
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer :: i,j
#ifdef MPI
      integer, dimension(8) :: reqs_khc,reqs_kvc
#endif

!--------------------------------------------------------------
!  lower surface:

      IF( bbc.eq.2 .or. bbc.eq.3 )THEN
        ! no slip (fix/check this later)
        ! semi-slip

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          ! Pope, pg 288:
!!!          tkea(i,j,1) = 3.33*ust(i,j)*ust(i,j)
          tkea(i,j,1) = 0.0
          kmh(i,j,1) = karman*zntmp(i,j)*ust(i,j)
        enddo
        enddo

      ELSE
        ! free slip

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          tkea(i,j,1) = 0.0
          kmh(i,j,1) = 0.0
        enddo
        enddo

      ENDIF

        if(timestats.ge.1) time_turb=time_turb+mytime()

        !-----
        call bc2d(tkea(ibt,jbt,1))
#ifdef MPI
        call comm_1s2d_start(tkea(ibt,jbt,1),kvcw1(1,1),kvcw2(1,1),kvce1(1,1),kvce2(1,1),   &
                                             kvcs1(1,1),kvcs2(1,1),kvcn1(1,1),kvcn2(1,1),reqs_kvc)
#endif
        call bc2d(kmh(ibc,jbc,1))
#ifdef MPI
        call comm_1s2d_start(kmh(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                            khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        !-----
        call comm_1s2d_end(tkea(ibt,jbt,1),kvcw1(1,1),kvcw2(1,1),kvce1(1,1),kvce2(1,1),   &
                                           kvcs1(1,1),kvcs2(1,1),kvcn1(1,1),kvcn2(1,1),reqs_kvc)
        call bcs2_2d(tkea(ibt,jbt,1))
        call comm_1s2d_end(kmh(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                          khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call bcs2_2d(kmh(ibc,jbc,1))
        !-----
        call comm_2d_corner(tkea(ibt,jbt,1))
        call comm_2d_corner(kmh(ibc,jbc,1))
        !-----
#endif

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          kmv(i,j,1) = kmh(i,j,1)
          khv(i,j,1) = kmv(i,j,1)*(khv(i,j,2)/(1.0e-10+kmv(i,j,2)))
          khh(i,j,1) = kmh(i,j,1)*(khh(i,j,2)/(1.0e-10+kmh(i,j,2)))
        enddo
        enddo

        if(timestats.ge.1) time_turb=time_turb+mytime()

!--------------------------------------------------------------
!  cm1r20: upper surface

      IF( tbc.eq.3 )THEN

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,nj
        do i=1,ni
          ! Pope, pg 288:
!!!          tkea(i,j,nk+1) = 3.33*ustt(i,j)*ustt(i,j)
          tkea(i,j,nk+1) = 0.0
          kmh(i,j,nk+1) = karman*zntmp(i,j)*ustt(i,j)
        enddo
        enddo

        if(timestats.ge.1) time_turb=time_turb+mytime()

        !-----
        call bc2d(tkea(ibt,jbt,nk+1))
#ifdef MPI
        call comm_1s2d_start(tkea(ibt,jbt,nk+1),kvcw1(1,1),kvcw2(1,1),kvce1(1,1),kvce2(1,1),   &
                                                kvcs1(1,1),kvcs2(1,1),kvcn1(1,1),kvcn2(1,1),reqs_kvc)
#endif
        call bc2d(kmh(ibc,jbc,nk+1))
#ifdef MPI
        call comm_1s2d_start(kmh(ibc,jbc,nk+1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                               khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        !-----
        call comm_1s2d_end(tkea(ibt,jbt,nk+1),kvcw1(1,1),kvcw2(1,1),kvce1(1,1),kvce2(1,1),   &
                                              kvcs1(1,1),kvcs2(1,1),kvcn1(1,1),kvcn2(1,1),reqs_kvc)
        call bcs2_2d(tkea(ibt,jbt,nk+1))
        call comm_1s2d_end(kmh(ibc,jbc,nk+1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                             khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call bcs2_2d(kmh(ibc,jbc,nk+1))
        !-----
        call comm_2d_corner(tkea(ibt,jbt,nk+1))
        call comm_2d_corner(kmh(ibc,jbc,nk+1))
        !-----
#endif

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          kmv(i,j,nk+1) = kmh(i,j,nk+1)
          khv(i,j,nk+1) = kmv(i,j,nk+1)*(khv(i,j,nk)/(1.0e-10+kmv(i,j,nk)))
          khh(i,j,nk+1) = kmh(i,j,nk+1)*(khh(i,j,nk)/(1.0e-10+kmh(i,j,nk)))
        enddo
        enddo

        if(timestats.ge.1) time_turb=time_turb+mytime()

      ENDIF


      return
      end subroutine tkebc


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbsmag(nstep,dt,ruh,rvh,rmh,mf,rmf,th0,rf,              &
                          nm,defv,defh,dum4,grdscl,lenh  ,zf,zntmp,ust,csm, &
                          kmh,kmv,khh,khv,lenscl,dissten,                  &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2, &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

      integer, intent(in) :: nstep
      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum4,grdscl,lenh
      real, intent(in), dimension(ib:ie,jb:je) :: zntmp,ust
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: csm
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: lenscl,dissten
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: tem,temx,temy

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

!!!      real, parameter :: cs      = 0.18
!!!      real, parameter :: csinv   = 1.0/cs
      real, parameter :: prandtl = 1.0/3.00
      real, parameter :: prinv   = 1.0/prandtl
      real, parameter :: dmin    = 1.0e-10

!-----------------------------------------------------------------------

      temx = 0.125*dx*dx/dt
      temy = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    jloop:  DO j=1,nj

    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz

      do k=2,nk
      do i=1,ni
        grdscl(i,j,k)=( ((dx*ruh(i))*(dy*rvh(j)))*(dz*rmf(i,j,k)) )**0.33333333
           ! cm1r20.1: revisit near-sfc grid scale at a later date
!!!        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                    &
!!!                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j))*csinv)**2)  &
!!!                               ) )
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical

      do i=1,ni
        tem=sqrt( (dx*ruh(i))*(dy*rvh(j)) )
        do k=2,nk
          lenh(i,j,k)=tem
        enddo
      enddo

      do k=2,nk
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
           ! cm1r20.1: revisit near-sfc grid scale at a later date
!!!        grdscl(i,j,k) = sqrt(1.0/( 1.0/(grdscl(i,j,k)**2)                                    &
!!!                                  +1.0/((karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j))*csinv)**2)  &
!!!                               ) )
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------

    IF(tconfig.eq.1)THEN

      do k=2,nk
      do i=1,ni
        kmh(i,j,k)=((csm(i,j,k)*grdscl(i,j,k))**2)     &
                 *sqrt( max(defv(i,j,k)+defh(i,j,k)-nm(i,j,k)*prinv,dmin) )
!!!        kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i)   &
!!!                                     , temy*rvh(j)*rvh(j) )
        kmv(i,j,k)=kmh(i,j,k)
        lenscl(i,j,k) = grdscl(i,j,k)
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

      do k=2,nk
      do i=1,ni
        kmh(i,j,k)=((csm(i,j,k)*lenh(i,j,k))**2)     &
                 *sqrt( max(defh(i,j,k),dmin) )
!!!        kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i)   &
!!!                                     , temy*rvh(j)*rvh(j) )
        kmv(i,j,k)=((csm(i,j,k)*grdscl(i,j,k))**2)     &
                 *sqrt( max(defv(i,j,k)-nm(i,j,k)*prinv,dmin) )
        lenscl(i,j,k) = grdscl(i,j,k)
      enddo
      enddo

    ENDIF

    ENDDO  jloop

    if( nstep.eq.0 .and. myid.eq.0 )then
      print *
!!!      print *,'  cs,csinv = ',cs,csinv
      print *,'  zf,grdscl:'
      i = 1
      j = 1
      do k=2,nk
        print *,k,(zf(i,j,k)-zf(i,j,1)),grdscl(i,j,k)
      enddo
    endif

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bcw(kmh,1)
      call bcw(kmv,1)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(kmh)
      call bct2(kmv)
#endif

!--------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        ! limit for numerical stability:
!!!        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i)   &
!!!                                     , temy*rvh(j)*rvh(j) )
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

!!!    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
    IF( tconfig.eq.1 )THEN
      do k=2,nk
      do i=1,ni
      if( csm(i,j,k).gt.cmemin )then
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/((csm(i,j,k)*grdscl(i,j,k))**4)
      endif
      enddo
      enddo
    ELSEIF( tconfig.eq.2 )THEN
      do k=2,nk
      do i=1,ni
      if( csm(i,j,k).gt.cmemin )then
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/((csm(i,j,k)*grdscl(i,j,k))**4)    &
                                        + (kmh(i,j,k)**3)/((csm(i,j,k)*lenh(i,j,k))**4)
      endif
      enddo
      enddo
    ENDIF
    ENDIF
!!!    ENDIF

    enddo

!--------------------------------------------------------------
!  cm1r18: surface

      IF( bbc.eq.3 )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          kmv(i,j,1) = karman*zntmp(i,j)*ust(i,j)
        enddo
        enddo

        call bc2d(kmv(ibc,jbc,1))
#ifdef MPI
        call comm_1s2d_start(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                            khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_1s2d_end(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                          khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_2d_corner(kmv(ibc,jbc,1))
        call bcs2_2d(kmv(ibc,jbc,1))
#endif

      IF( tconfig.eq.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          kmh(i,j,1) = kmv(i,j,1)
          khv(i,j,1) = kmv(i,j,1)*prinv
          khh(i,j,1) = khv(i,j,1)
        enddo
        enddo
      ELSEIF( tconfig.eq.2 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          khv(i,j,1) = kmv(i,j,1)*prinv
        enddo
        enddo
      ENDIF

      ENDIF

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine turbsmag


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbparam_horiz(nstep,zf,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                          nm,defv,defh,lvz,kmh,khh,cm0,dissten,out3d,zs,zntmp,ust,xland,psfc,tlh,  &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                         &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,         &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in) :: dt
      logical, intent(in) :: dosfcflx
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0
      real, intent(in), dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: lvz
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,khh
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: dissten
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      real, intent(in), dimension(ib:ie,jb:je) :: zs,zntmp,ust,xland,psfc
      real, intent(inout), dimension(ib:ie,jb:je) :: tlh
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: tem,tem1,temx,temy

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

      real, parameter :: prandtl = 1.0
      real, parameter :: prinv   = 1.0/prandtl
      real, parameter :: dmin    = 1.0e-10

!--------------------------------------------------------------
!  Smagorinsky-type scheme for parameterized turbulence:
!--------------------------------------------------------------
!  Interior:

  lhcheck:  &
  IF( l_h.gt.1.0e-12 .or. lhref1.gt.1.0e-12 .or. lhref2.gt.1.0e-12 )THEN

    if(ny.eq.1)then
      temx =  0.250*dx*dx/dt
      temy = 1000.0*dy*dy/dt
    elseif(nx.eq.1)then
      temx = 1000.0*dx*dx/dt
      temy =  0.250*dy*dy/dt
    else
      temx =  0.125*dx*dx/dt
      temy =  0.125*dy*dy/dt
    endif

    ! cm1r18:
    ! Over water, make tlh a function of surface pressure.
    !   (designed for hurricanes)
    ! Over land, simply set to tlh to l_h.
!$omp parallel do default(shared)   &
!$omp private(i,j)
    do j=1,nj
    do i=1,ni
      IF( (xland(i,j).gt.1.5) .and. (zs(i,j).lt.1.0) )THEN
        ! over water (sea level only):
        tlh(i,j) = lhref2+(lhref1-lhref2)   &
                      *(psfc(i,j)-90000.0)  &
                      /( 101500.0-90000.0)
      ELSE
        ! all other cases:
        tlh(i,j) = l_h
      ENDIF
    enddo
    enddo

  IF( cm1setup.eq.4 )THEN
    ! get km only in mesoscale-model part of domain:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
    if( cm0(i,j).le.cmemin )then
      kmh(i,j,k)=(tlh(i,j)**2)*sqrt( max(defh(i,j,k),dmin) )
!!!      kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
    endif
    enddo
    enddo
    enddo

  ELSE

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      kmh(i,j,k)=(tlh(i,j)**2)*sqrt( max(defh(i,j,k),dmin) )
!!!      kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
    enddo
    enddo
    enddo

  ENDIF

!--------------------------------------------------------------
! boundary conditions:

      if(timestats.ge.1) time_turb=time_turb+mytime()

      call bcw(kmh,1)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(kmh)
#endif

        ! Extrapolate:
        do j=0,nj+1
        do i=0,ni+1
          kmh(i,j,1) = max( 0.0 , 2.0*kmh(i,j,2) - kmh(i,j,3) )
        enddo
        enddo

!--------------------------------------------------------------
!  calculate Kh
!  and also limit horizontal coeffs for numerical stability:

  IF( cm1setup.eq.4 )THEN
    ! get km only in mesoscale-model part of domain:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem1)
    do j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
      if( cm0(i,j).le.cmemin )then
        khh(i,j,k)=kmh(i,j,k)*prinv
!!!        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
      endif
      enddo
      enddo

    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
      do k=2,nk
      do i=1,ni
      if( cm0(i,j).le.cmemin )then
        dissten(i,j,k) = dissten(i,j,k) + (kmh(i,j,k)**3)/(tlh(i,j)**4)
      endif
      enddo
      enddo
    ENDIF
    ENDIF

    enddo

  ELSE

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem1)
    do j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
!!!        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
      enddo
      enddo

    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
      do k=2,nk
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k) + (kmh(i,j,k)**3)/(tlh(i,j)**4)
      enddo
      enddo
    ENDIF
    ENDIF

    enddo

  ENDIF

!--------------------------------------------------------------

  ENDIF  lhcheck

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine turbparam_horiz


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbparam_vert(nstep,zf,dt,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,rf, &
                          nm,defv,defh,lvz,kmv,khv,cm0,dissten,out3d,zs,zntmp,ust,xland,  &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                         &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,         &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      use input
      use constants
      use bc_module
      use comm_module
      implicit none

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in) :: dt
      logical, intent(in) :: dosfcflx
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: th0
      real, intent(in), dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: lvz
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv,khv
      real, intent(in), dimension(ib:ie,jb:je) :: cm0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: dissten
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
      real, intent(in), dimension(ib:ie,jb:je) :: zs,zntmp,ust,xland
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: rlinf,tem,tem1

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

      real, parameter :: prandtl = 1.0
      real, parameter :: prinv   = 1.0/prandtl
      real, parameter :: dmin    = 1.0e-10

!--------------------------------------------------------------
!  Smagorinsky-type scheme for parameterized turbulence:
!--------------------------------------------------------------
!  Interior:

  linfcheck:  &
  IF( l_inf.gt.1.0e-12 )THEN

!!!    tem = 1.0/(1.0e-6+l_inf)
    rlinf = (1.0e-6+l_inf)**(-2)

  IF( cm1setup.eq.4 )THEN
    ! get km only in mesoscale-model part of domain:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
    if( cm0(i,j).le.cmemin )then
      lvz(i,j,k)=sqrt( ( rlinf + (karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j)))**(-2) )**(-1) )
      kmv(i,j,k)=(lvz(i,j,k)**2)*sqrt( max(defv(i,j,k)-nm(i,j,k)*prinv,dmin) )
    endif
    enddo
    enddo
    enddo

  ELSE

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      lvz(i,j,k)=sqrt( ( rlinf + (karman*((zf(i,j,k)-zf(i,j,1))+zntmp(i,j)))**(-2) )**(-1) )
      kmv(i,j,k)=(lvz(i,j,k)**2)*sqrt( max(defv(i,j,k)-nm(i,j,k)*prinv,dmin) )
    enddo
    enddo
    enddo

  ENDIF

!--------------------------------------------------------------
! boundary conditions:

      if(timestats.ge.1) time_turb=time_turb+mytime()

      call bcw(kmv,1)
#ifdef MPI
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call bct2(kmv)
#endif

!--------------------------------------------------------------
!  calculate Kh

  IF( cm1setup.eq.4 )THEN
    ! get km only in mesoscale-model part of domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
    do j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
      if( cm0(i,j).le.cmemin )then
        khv(i,j,k)=kmv(i,j,k)*prinv
      endif
      enddo
      enddo

    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
      do k=2,nk
      do i=1,ni
      if( cm0(i,j).le.cmemin )then
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/(lvz(i,j,k)**4)
      endif
      enddo
      enddo
    ENDIF
    ENDIF

    enddo

  ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1)
    do j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

    IF( idiss.eq.1 .or. output_dissten.eq.1 )THEN
    IF( j.ge.1 .and. j.le.nj )THEN
      do k=2,nk
      do i=1,ni
        dissten(i,j,k) = dissten(i,j,k) + (kmv(i,j,k)**3)/(lvz(i,j,k)**4)
      enddo
      enddo
    ENDIF
    ENDIF

    enddo

  ENDIF

!--------------------------------------------------------------
!  cm1r18: surface

      IF( bbc.eq.3 )THEN

        do j=1,nj
        do i=1,ni
          kmv(i,j,1) = karman*zntmp(i,j)*ust(i,j)
        enddo
        enddo

        call bc2d(kmv(ibc,jbc,1))
#ifdef MPI
        call comm_1s2d_start(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                            khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_1s2d_end(kmv(ibc,jbc,1),khcw1(1,1),khcw2(1,1),khce1(1,1),khce2(1,1),   &
                                          khcs1(1,1),khcs2(1,1),khcn1(1,1),khcn2(1,1),reqs_khc)
        call comm_2d_corner(kmv(ibc,jbc,1))
        call bcs2_2d(kmv(ibc,jbc,1))
#endif

        do j=0,nj+1
        do i=0,ni+1
          khv(i,j,1) = kmv(i,j,1)*prinv
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------

  ENDIF linfcheck

      if(nstep.eq.1.and.myid.eq.0)then
        print *,'  k,zf,lvz:  znt = ',zntmp(1,1)
        do k=2,nk
          print *,k,(zf(1,1,k)-zf(1,1,1)),lvz(1,1,k)
        enddo
      endif

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine turbparam_vert


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      ! gettau
      subroutine gettau(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,ustt,ut,vt,st,rf,mf,dum1,dum2,dum3,dum4,dum5,dum6, &
                        kmh,kmv,t11,t12,t13,t22,t23,t33,ua,ugr,va,vgr,wa,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,bndy,kbdy)
      use input
      use constants
      use ib_module
      implicit none
      
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: ust,u1,v1,s1,ustt,ut,vt,st
      real, intent(inout), dimension(ib:ie,jb:je) :: stau
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua,ugr
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va,vgr
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      double precision, intent(in) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
        
      integer i,j,k,kk
      real :: tem
      logical :: doit

!----------------------------------------------------------------------
!
!  This subroutine calculates the subgrid stress terms.
!
!    t_ij  =  2 * rho * K * S_ij
!
!  NOTE:  upon entering this subroutine, the t_ij arrays must already 
!         contain rho * S_ij  (see calcdef subroutine)
!
!  Since cm1r18, surface stress (ie, surface drag) is incorporated into
!  the stress arrays here.
!
!  Note:  Turbulent viscosities are defined on w points.
!
!  Note:  For axisymmetric simulations, t11 and t12 herein are 
!         actually not stresses:  the actual stresses are
!         combined in a convienent form for the sake of flux-form
!         calculations in the turbu and turbv subroutines.
!         Also note that t22 is never calculated.
!         So, if you need the actual stress components for something, 
!         beware that you will need to re-calculate t11,t12,t22.
!
!----------------------------------------------------------------------

  IF(axisymm.eq.0)THEN

    ! Cartesian grid:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        !  2.0 * 0.5 = 1.0
        tem = (kmh(i,j,k)+kmh(i,j,k+1))
        t11(i,j,k)=t11(i,j,k)*tem
        t22(i,j,k)=t22(i,j,k)*tem
        t33(i,j,k)=t33(i,j,k)*tem
        !  2.0 * 0.125 = 0.25
        t12(i,j,k)=t12(i,j,k)*0.25                                            &
     *( ( (kmh(i-1,j-1,k  )+kmh(i,j,k  ))+(kmh(i-1,j,k  )+kmh(i,j-1,k  )) )   &
       +( (kmh(i-1,j-1,k+1)+kmh(i,j,k+1))+(kmh(i-1,j,k+1)+kmh(i,j-1,k+1)) ) )
      enddo
      enddo
    ENDDO

          !-----
          ! lateral boundary conditions:
          if(wbc.eq.3.and.ibw.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(j,k)
            do k=1,nk
            do j=1,nj+1
              t12(1,j,k) = t12(2,j,k)
            enddo
            enddo
          endif
          if(ebc.eq.3.and.ibe.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(j,k)
            do k=1,nk
            do j=1,nj+1
              t12(ni+1,j,k) = t12(ni,j,k)
            enddo
            enddo
          endif
          !-----
          !-----
          if(sbc.eq.3.and.ibs.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(i,k)
            do k=1,nk
            do i=1,ni+1
              t12(i,1,k) = t12(i,2,k)
            enddo
            enddo
          endif
          if(nbc.eq.3.and.ibn.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(i,k)
            do k=1,nk
            do i=1,ni+1
              t12(i,nj+1,k) = t12(i,nj,k)
            enddo
            enddo
          endif
          !-----

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=2,nk
      do j=1,nj+1
      do i=1,ni+1
        !  2.0 x 0.5 = 1.0
        t13(i,j,k)=t13(i,j,k)*( kmv(i-1,j,k)+kmv(i,j,k) )
        t23(i,j,k)=t23(i,j,k)*( kmv(i,j-1,k)+kmv(i,j,k) )
      enddo
      enddo
    enddo
            !-----
            ! lateral boundary conditions:
            if(wbc.eq.3.and.ibw.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(j,k)
              do k=2,nk
              do j=1,nj
                t13(1,j,k) = t13(2,j,k)
              enddo
              enddo
            endif
            if(ebc.eq.3.and.ibe.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(j,k)
              do k=2,nk
              do j=1,nj
                t13(ni+1,j,k) = t13(ni,j,k)
              enddo
              enddo
            endif
            !-----
            !-----
            if(sbc.eq.3.and.ibs.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(i,k)
              do k=2,nk
              do i=1,ni
                t23(i,1,k) = t23(i,2,k)
              enddo
              enddo
            endif
            if(nbc.eq.3.and.ibn.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(i,k)
              do k=2,nk
              do i=1,ni
                t23(i,nj+1,k) = t23(i,nj,k)
              enddo
              enddo
            endif
            !-----

!------------------------------------

  ELSE

    ! axisymmetric grid:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    DO k=1,nk
      do j=1,nj
      do i=1,ni+1
        !  2.0 * 0.5 = 1.0
        tem = (kmh(i,j,k)+kmh(i,j,k+1))
        t11(i,j,k)=t11(i,j,k)*tem
        t33(i,j,k)=t33(i,j,k)*tem
        !  2.0 * 0.25  =  0.5
        t12(i,j,k)=0.5*t12(i,j,k)*( arf2(i)*(kmh(i  ,j,k+1)+kmh(i  ,j,k)) &
                                   +arf1(i)*(kmh(i-1,j,k+1)+kmh(i-1,j,k)) )
      enddo
      enddo
    ENDDO

          !-----
          ! lateral boundary conditions:
          j = 1
          if(wbc.eq.3)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(k)
            do k=1,nk
!!!              t12(1,j,k) = t12(2,j,k)
              t12(1,j,k) = 0.0
            enddo
          endif
          if(ebc.eq.3)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(k)
            do k=1,nk
              t12(ni+1,j,k) = t12(ni,j,k)
            enddo
          endif
          !-----

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    DO k=2,nk
      do j=1,nj
      do i=1,ni+1
        !  2.0 * 0.5  =  1.0
        t13(i,j,k)=t13(i,j,k)*( arf1(i)*kmv(i-1,j,k)+arf2(i)*kmv(i,j,k) )
        t23(i,j,k)=2.0*t23(i,j,k)*kmv(i,j,k)
      enddo
      enddo
    ENDDO

  ENDIF

!------------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=0,nj+1
            t11(0,j,k) = t11(1,j,k)
          enddo
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=0,nj+1
            t11(ni+1,j,k) = t11(ni,j,k)
          enddo
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=0,ni+1
            t22(i,0,k) = t22(i,1,k)
          enddo
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=0,ni+1
            t22(i,nj+1,k) = t22(i,nj,k)
          enddo
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj+1
            t12(1,j,k) = t12(2,j,k)
          enddo
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj+1
            t12(ni+1,j,k) = t12(ni,j,k)
          enddo
          enddo
        ENDIF
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni+1
            t12(i,1,k) = t12(i,2,k)
          enddo
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni+1
            t12(i,nj+1,k) = t12(i,nj,k)
          enddo
          enddo
        ENDIF
        !-----
        ! corner points:
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            t12(1,1,k) = t12(2,2,k)
          enddo
        ENDIF
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            t12(ni+1,1,k) = t12(ni,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            t12(1,nj+1,k) = t12(2,nj,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            t12(ni+1,nj+1,k) = t12(ni,nj,k)
          enddo
        ENDIF
        !-----
    ENDIF

!--------------------------------------------------------------
!  lower boundary conditions

            call sfcstress(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,rf,mf,dum1,dum2,  &
                        kmh,kmv,t13,t23,ua,ugr,va,vgr,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv)

!--------------------------------------------------------------
!  upper boundary conditions

    IF(tbc.eq.1)THEN
      ! free slip:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,nk+1)=t13(i,j,nk)
        t23(i,j,nk+1)=t23(i,j,nk)
      enddo
      enddo

    ELSEIF(tbc.eq.2)THEN
      ! no slip:

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,nk+1) = -2.0*ugr(i,j,nk)*rdz*0.5*( mf(i-1,j,nk+1)+ mf(i,j,nk+1))  &
                                              *0.5*( rf(i-1,j,nk+1)+ rf(i,j,nk+1))  &
                                              *0.5*(kmv(i-1,j,nk  )+kmv(i,j,nk  ))
          t23(i,j,nk+1) = -2.0*vgr(i,j,nk)*rdz*0.5*( mf(i,j-1,nk+1)+ mf(i,j,nk+1))  &
                                              *0.5*( rf(i,j-1,nk+1)+ rf(i,j,nk+1))  &
                                              *0.5*(kmv(i,j-1,nk  )+kmv(i,j,nk  ))
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,nk+1) = -2.0*ugr(i,j,nk)*rdz*mf(1,1,nk+1)                                       &
                                              *0.5*(arf1(i)*rf(i-1,j,nk+1)+arf2(i)*rf(i,j,nk+1))  &
                                              *0.5*(arf1(i)*kmv(i-1,j,nk)+arf2(i)*kmv(i,j,nk))
          t23(i,j,nk+1) = -2.0*vgr(i,j,nk)*rdz*mf(i,j,nk+1)  &
                                              *rf(i,j,nk+1)  &
                                              *kmv(i,j,nk)
        enddo
        enddo

      ENDIF

    ELSEIF(tbc.eq.3)THEN
      !--------------------------------------------------------!
      !-------- (this is where "drag" is set for tbc=3) -------!

      !$omp parallel do default(shared)   &
      !$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,nk+1) = -0.25*( (ustt(i-1,j)**2)*(ut(i-1,j)/max(st(i-1,j),0.01))    &
                               +(ustt(i  ,j)**2)*(ut(i  ,j)/max(st(i  ,j),0.01)) )  &
                             *( rf(i-1,j,nk+1)+rf(i,j,nk+1) )
        t23(i,j,nk+1) = -0.25*( (ustt(i,j-1)**2)*(vt(i,j-1)/max(st(i,j-1),0.01))    &
                               +(ustt(i,j  )**2)*(vt(i,j  )/max(st(i,j  ),0.01)) )  &
                             *( rf(i,j-1,nk+1)+rf(i,j,nk+1) )
      enddo
      enddo

    ENDIF

!--------------------------------------------------------------

    IF( axisymm.eq.1 )THEN
      ! lateral boundary condition:
      !$omp parallel do default(shared)   &
      !$omp private(k)
      do k=0,nk+1
        t13(1,1,k)=0.0
      enddo
    ENDIF

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()
 
      return
      end subroutine gettau
      ! gettau


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine sfcstress(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,rf,mf,dum1,dum2,  &
                        kmh,kmv,t13,t23,ua,ugr,va,vgr,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv)
      use input
      use constants
      implicit none
      
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: ust,u1,v1,s1
      real, intent(inout), dimension(ib:ie,jb:je) :: stau
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua,ugr
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va,vgr
      double precision, intent(in) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv
        
      integer :: i,j,k
      real :: tem1,tem2,s8u,s8v

    IF(bbc.eq.1)THEN
      ! free slip:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,1)=t13(i,j,2)
        t23(i,j,1)=t23(i,j,2)
      enddo
      enddo

    ELSEIF(bbc.eq.2)THEN
      ! no slip:

      IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,1) = 2.0*ugr(i,j,1)*rdz*0.5*( mf(i-1,j,1)+ mf(i,j,1))  &
                                         *0.5*( rf(i-1,j,1)+ rf(i,j,1))  &
                                         *0.5*(kmv(i-1,j,2)+kmv(i,j,2))
          t23(i,j,1) = 2.0*vgr(i,j,1)*rdz*0.5*( mf(i,j-1,1)+ mf(i,j,1))  &
                                         *0.5*( rf(i,j-1,1)+ rf(i,j,1))  &
                                         *0.5*(kmv(i,j-1,2)+kmv(i,j,2))
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,1)=2.0*ugr(i,j,1)*rdz*mf(1,1,1)                                    &
                                       *0.5*(arf1(i)*rf(i-1,j,1)+arf2(i)*rf(i,j,1))  &
                                       *0.5*( arf1(i)*kmv(i-1,j,2)+arf2(i)*kmv(i,j,2) )
          t23(i,j,1)=2.0*vgr(i,j,1)*rdz*mf(i,j,1)  &
                                       *rf(i,j,1)  &
                                       *kmv(i,j,2)
        enddo
        enddo

      ENDIF

    ELSEIF(bbc.eq.3)THEN
      !--------------------------------------------------------!
      !--------  surface stress for semi-slip lower bc --------!
      !-------- (this is where "drag" is set for bbc=3) -------!

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

        avgsfc:  &
        IF( ( .not. use_avg_sfc ) .or. ( sfcmodel.eq.2 .or. sfcmodel.eq.3 ) )THEN

          !$omp parallel do default(shared)   &
          !$omp private(i,j)
          do j=1,nj+1
          do i=1,ni+1
            t13(i,j,1) = 0.25*( (ust(i-1,j)**2)*(u1(i-1,j)/max(s1(i-1,j),0.01))    &
                               +(ust(i  ,j)**2)*(u1(i  ,j)/max(s1(i  ,j),0.01)) )  &
                             *( rf(i-1,j,1)+rf(i,j,1) )
            t23(i,j,1) = 0.25*( (ust(i,j-1)**2)*(v1(i,j-1)/max(s1(i,j-1),0.01))    &
                               +(ust(i,j  )**2)*(v1(i,j  )/max(s1(i,j  ),0.01)) )  &
                             *( rf(i,j-1,1)+rf(i,j,1) )
          enddo
          enddo

        ELSE  avgsfc

          ! Moeng (1984):

!!!          do j=0,nj+1
!!!          do i=0,ni+1
!!!              dum1(i,j,1) = ust(i,j)*ust(i,j)        &
!!!                           *( s1(i,j)*avgsfcu           &
!!!                             +avgsfcs*(u1(i,j)-avgsfcu) )  &
!!!                           /( max(0.0001,avgsfcs*sqrt( avgsfcu**2 + avgsfcv**2 )) )
!!!              dum2(i,j,1) = ust(i,j)*ust(i,j)        &
!!!                           *( s1(i,j)*avgsfcv           &
!!!                             +avgsfcs*(v1(i,j)-avgsfcv) )  &
!!!                           /( max(0.0001,avgsfcs*sqrt( avgsfcu**2 + avgsfcv**2 )) )
!!!          enddo
!!!          enddo
!!!          do j=1,nj+1
!!!          do i=1,ni+1
!!!            t13(i,j,1) = 0.25*( dum1(i-1,j,1)+dum1(i,j,1) )  &
!!!                             *( rf(i-1,j,1)+rf(i,j,1) )
!!!            t23(i,j,1) = 0.25*( dum2(i,j-1,1)+dum2(i,j,1) )  &
!!!                             *( rf(i,j-1,1)+rf(i,j,1) )
!!!          enddo
!!!          enddo

          ! cm1r20.3: calculate stresses directly at staggered grid points

          tem1 = ust(1,1)*ust(1,1)/( max(1.0e-8,avgsfcsu*sqrt( avgsfcu**2 + avgsfcv**2 )) )
          tem2 = ust(1,1)*ust(1,1)/( max(1.0e-8,avgsfcsv*sqrt( avgsfcu**2 + avgsfcv**2 )) )

          !$omp parallel do default(shared)   &
          !$omp private(i,j,s8u,s8v)
          do j=1,nj+1
          do i=1,ni+1
            s8u = sqrt( ugr(i,j,1)**2   &
                       + ( 0.25*( (vgr(i  ,j,1)+vgr(i  ,j+1,1)) &
                                 +(vgr(i-1,j,1)+vgr(i-1,j+1,1)) ) )**2 )
            s8v = sqrt( vgr(i,j,1)**2   &
                       + ( 0.25*( (ugr(i,j  ,1)+ugr(i+1,j  ,1)) &
                                 +(ugr(i,j-1,1)+ugr(i+1,j-1,1)) ) )**2 )
            t13(i,j,1) = 0.5*(rf(i-1,j,1)+rf(i,j,1))  &
                       *tem1*( s8u*avgsfcu + avgsfcsu*(ugr(i,j,1)-avgsfcu) )
            t23(i,j,1) = 0.5*(rf(i,j-1,1)+rf(i,j,1))  &
                       *tem2*( s8v*avgsfcv + avgsfcsv*(vgr(i,j,1)-avgsfcv) )
          enddo
          enddo

        ENDIF  avgsfc

      ELSE

        ! axisymmetric grid:
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,1) = 0.25*( arf1(i)*(ust(i-1,j)**2)*(u1(i-1,j)/max(s1(i-1,j),0.01))    &
                           +arf2(i)*(ust(i  ,j)**2)*(u1(i  ,j)/max(s1(i  ,j),0.01)) )  &
                         *(arf1(i)*rf(i-1,j,1)+arf2(i)*rf(i,j,1))
        t23(i,j,1) = rf(i,j,1)*(ust(i,j)**2)*(v1(i,j)/max(s1(i,j),0.01))
      enddo
      enddo

      ENDIF

      IF( testcase.eq.7 )THEN

        ! Replace normal sfc stress ... use lowest-model-level wind speed
        ! (VanZanten et al 2011, JAMES)

        !$omp parallel do default(shared)   &
        !$omp private(i,j)
        do j=1,nj+1
        do i=1,ni+1
          t13(i,j,1) =     ( 0.001229 * 0.5*(s1(i-1,j)+s1(i,j)) * ugr(i,j,1) )  &
                      *0.5*( rf(i-1,j,1)+rf(i,j,1) )
          t23(i,j,1) =     ( 0.001229 * 0.5*(s1(i,j-1)+s1(i,j)) * vgr(i,j,1) )  &
                      *0.5*( rf(i,j-1,1)+rf(i,j,1) )
        enddo
        enddo

      ENDIF

    ENDIF

    !$omp parallel do default(shared)   &
    !$omp private(i,j)
    do j=1,nj
    do i=1,ni
      stau(i,j) = sqrt( (0.5*(t13(i,j,1)+t13(i+1,j,1)))**2  &
                       +(0.5*(t23(i,j,1)+t23(i,j+1,1)))**2  &
                      )/rf(i,j,1)
    enddo
    enddo

      end subroutine sfcstress


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      ! calcdef
      subroutine calcdef(rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                     xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,mh,c1,c2,mf,defv,defh,  &
                     dum1,dum2,ua,va,wa,s11,s12,s13,s22,s23,s33,gx,gy,rho,rr,rf)
      use input
      use constants
      implicit none

      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: defv,defh
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: s11,s12,s13,s22,s23,s33
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gy
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
        
      integer :: i,j,k
      real :: r1,r2,r3,r4
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,rrf
      real :: temz

!----------------------------------------------------------------------
!
!  This subroutine calculates the strain rate terms
!
!    S_ij  =  0.5 * ( d(u_i)/d(x_j) + d(u_j)/d(x_i) )
!
!  (note: multiplied by density herein)
!  and then uses these variables to calculate deformation.
!
!  Note:
!  Since cm1r18, surface stress (ie, surface drag) is no longer 
!  calculated in this subroutine.  See gettau subroutine instead.
!
!  Note:  For axisymmetric simulations, s11 and s12 herein are 
!         actually not rate-of-strain components:  the actual 
!         components have been combined mathematically in a 
!         way to be consistent with the flux-form calculations 
!         in the turbu and turbv subroutines.
!         Also note that s22 is never calculated.
!         So, if you need the actual strain components for something, 
!         beware that you will need to re-calculate s11,s12,s22.
!
!----------------------------------------------------------------------

  IF(.not.terrain_flag)THEN

  IF( axisymm.eq.0 )THEN
    ! Cartesian without terrain:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,temz)
    DO k=1,nk

      temz = rdz*mh(1,1,k)

      do j=0,nj+1
      do i=0,ni+1 
        s11(i,j,k)=rho(i,j,k)*(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
        s22(i,j,k)=rho(i,j,k)*(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
        s33(i,j,k)=rho(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*temz
        s12(i,j,k)=0.5*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                        +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) ) &
              *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo       

    ENDDO

          !-----
          ! lateral boundary conditions:
          if(wbc.eq.3.and.ibw.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(j,k)
            do k=1,nk
            do j=1,nj+1
              s12(1,j,k) = s12(2,j,k)
            enddo
            enddo
          elseif(wbc.eq.4.and.ibw.eq.1)then
            ! no slip b.c.
            i = 1
            !$omp parallel do default(shared)   &
            !$omp private(j,k)
            do k=1,nk
            do j=1,nj+1
              s12(1,j,k) = 2.0*va(1,j,k)*rdx*uf(1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
            enddo
          endif
          if(ebc.eq.3.and.ibe.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(j,k)
            do k=1,nk
            do j=1,nj+1
              s12(ni+1,j,k) = s12(ni,j,k)
            enddo
            enddo
          elseif(ebc.eq.4.and.ibe.eq.1)then
            ! no slip b.c.
            i = ni+1
            !$omp parallel do default(shared)   &
            !$omp private(j,k)
            do k=1,nk
            do j=1,nj+1
              s12(ni+1,j,k) = -2.0*va(ni,j,k)*rdx*uf(ni+1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
            enddo
          endif
          !-----
          !-----
          if(sbc.eq.3.and.ibs.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(i,k)
            do k=1,nk
            do i=1,ni+1
              s12(i,1,k) = s12(i,2,k)
            enddo
            enddo
          elseif(sbc.eq.4.and.ibs.eq.1)then
            ! no slip b.c.
            j = 1
            !$omp parallel do default(shared)   &
            !$omp private(i,k)
            do k=1,nk
            do i=1,ni+1
              s12(i,1,k) = 2.0*ua(i,1,k)*rdy*vf(1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
            enddo
          endif
          if(nbc.eq.3.and.ibn.eq.1)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(i,k)
            do k=1,nk
            do i=1,ni+1
              s12(i,nj+1,k) = s12(i,nj,k)
            enddo
            enddo
          elseif(nbc.eq.4.and.ibn.eq.1)then
            ! no slip b.c.
            j = nj+1
            !$omp parallel do default(shared)   &
            !$omp private(i,k)
            do k=1,nk
            do i=1,ni+1
              s12(i,nj+1,k) = -2.0*ua(i,nj,k)*rdy*vf(nj+1)   &
                   *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
            enddo
            enddo
          endif
          !-----

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    DO k=2,nk

      do j=1,nj+1
      do i=1,ni+1 
        s13(i,j,k)=0.5*( (wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                        +(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))  &
                       )*0.5*( rf(i-1,j,k)+rf(i,j,k) )
        s23(i,j,k)=0.5*( (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                        +(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))  &
                       )*0.5*( rf(i,j-1,k)+rf(i,j,k) )
      enddo
      enddo       

    ENDDO

            !-----
            ! lateral boundary conditions:
            if(wbc.eq.3.and.ibw.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(j,k)
              do k=2,nk
              do j=1,nj
                s13(1,j,k) = s13(2,j,k)
              enddo
              enddo
            elseif(wbc.eq.4.and.ibw.eq.1)then
              ! no slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(j,k)
              do k=2,nk
              do j=1,nj
                s13(1,j,k) = 2.0*wa(1,j,k)*rdx*uf(1)
              enddo
              enddo
            endif
            if(ebc.eq.3.and.ibe.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(j,k)
              do k=2,nk
              do j=1,nj
                s13(ni+1,j,k) = s13(ni,j,k)
              enddo
              enddo
            elseif(ebc.eq.4.and.ibe.eq.1)then
              ! no slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(j,k)
              do k=2,nk
              do j=1,nj
                s13(ni+1,j,k) = -2.0*wa(ni,j,k)*rdx*uf(ni+1)
              enddo
              enddo
            endif
            !-----

            !-----
            if(sbc.eq.3.and.ibs.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(i,k)
              do k=2,nk
              do i=1,ni
                s23(i,1,k) = s23(i,2,k)
              enddo
              enddo
            elseif(sbc.eq.4.and.ibs.eq.1)then
              ! no slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(i,k)
              do k=2,nk
              do i=1,ni
                s23(i,1,k) = 2.0*wa(i,1,k)*rdy*vf(1)
              enddo
              enddo
            endif
            if(nbc.eq.3.and.ibn.eq.1)then
              ! free slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(i,k)
              do k=2,nk
              do i=1,ni
                s23(i,nj+1,k) = s23(i,nj,k)
              enddo
              enddo
            elseif(nbc.eq.4.and.ibn.eq.1)then
              ! no slip b.c.
              !$omp parallel do default(shared)   &
              !$omp private(i,k)
              do k=2,nk
              do i=1,ni
                s23(i,nj+1,k) = -2.0*wa(i,nj,k)*rdy*vf(nj+1)
              enddo
              enddo
            endif
            !-----

!-------------------------------------------------------------------------------

  ELSE
    ! axisymmetric:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    DO k=1,nk

      do j=1,nj
      do i=1,ni+1
        s11(i,j,k)=rho(i,j,k)*(ua(i+1,j,k)*arf1(i+1)-ua(i,j,k)*arf2(i))*rdx*uh(i)
        s33(i,j,k)=rho(i,j,k)*(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)
        !  0.5 * 0.5  =  0.25
        s12(i,j,k)=0.25*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))   &
                       *(arh1(i)*va(i,j,k)-arh2(i-1)*va(i-1,j,k))*rdx*uf(i)
      enddo
      enddo

    ENDDO

          !-----
          ! lateral boundary conditions:
          j = 1
          if(wbc.eq.3)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(k)
            do k=1,nk
!!!              s12(1,j,k) = s12(2,j,k)
              s12(1,j,k) = 0.0
            enddo
          elseif(wbc.eq.4)then
            ! no slip b.c.
            i = 1
            !$omp parallel do default(shared)   &
            !$omp private(k)
            do k=1,nk
              s12(1,j,k) = 2.0*va(1,j,k)*rdx*uf(1)   &
                        *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
            enddo
          endif
          if(ebc.eq.3)then
            ! free slip b.c.
            !$omp parallel do default(shared)   &
            !$omp private(k)
            do k=1,nk
              s12(ni+1,j,k) = s12(ni,j,k)
            enddo
          elseif(ebc.eq.4)then
            ! no slip b.c.
            i = ni+1
            !$omp parallel do default(shared)   &
            !$omp private(k)
            do k=1,nk
              s12(ni+1,j,k) = -2.0*va(ni,j,k)*rdx*uf(ni+1)   &
                        *0.5*(arf1(i)*rho(i-1,j,k)+arf2(i)*rho(i,j,k))
            enddo
          endif
          !-----

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    DO k=2,nk

      do j=1,nj
      do i=1,ni+1
        !  0.5 * 0.5  =  0.25
        s13(i,j,k)=0.25*(arf1(i)*rf(i-1,j,k)+arf2(i)*rf(i,j,k))  &
                       *( (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)  &
                         +(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i) )
        s23(i,j,k)=0.5*rf(i,j,k)*(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
      enddo
      enddo

    ENDDO

  ENDIF

!-------------------------------------------------------------------------------
!  Cartesian with terrain:

  ELSE

    ! dum1 stores u at w-pts:
    ! dum2 stores v at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
    do j=0,nj+2
      ! lowest model level:
      do i=0,ni+2
        dum1(i,j,1) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
        dum2(i,j,1) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo

      ! upper-most model level:
      do i=0,ni+2
        dum1(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
        dum2(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo

      ! interior:
      do k=2,nk
      r2 = (sigmaf(k)-sigma(k-1))*rds(k)
      r1 = 1.0-r2
      do i=0,ni+2
        dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
        dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
      enddo
      enddo
    enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1 
        s11(i,j,k)=gz(i,j)*(ua(i+1,j,k)*rgzu(i+1,j)-ua(i,j,k)*rgzu(i,j))*rdx*uh(i) &
                  +( gx(i,j,k+1)*(dum1(i,j,k+1)+dum1(i+1,j,k+1))      &
                    -gx(i,j,k  )*(dum1(i,j,k  )+dum1(i+1,j,k  ))      &
                   )*0.5*rdsf(k)
        s11(i,j,k)=s11(i,j,k)*rho(i,j,k)
        s22(i,j,k)=gz(i,j)*(va(i,j+1,k)*rgzv(i,j+1)-va(i,j,k)*rgzv(i,j))*rdy*vh(j) &
                  +( gy(i,j,k+1)*(dum2(i,j,k+1)+dum2(i,j+1,k+1))      &
                    -gy(i,j,k  )*(dum2(i,j,k  )+dum2(i,j+1,k  ))      &
                   )*0.5*rdsf(k)
        s22(i,j,k)=s22(i,j,k)*rho(i,j,k)
        s33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdsf(k)*gz(i,j)
        s33(i,j,k)=s33(i,j,k)*rho(i,j,k)
      enddo
      enddo
      do j=1,nj+1 
      do i=1,ni+1
        r1 = 0.25*( ( rho(i-1,j-1,k)*gz(i-1,j-1)   &
                     +rho(i  ,j  ,k)*gz(i  ,j  ) ) &
                   +( rho(i-1,j  ,k)*gz(i-1,j  )   &
                     +rho(i  ,j-1,k)*gz(i  ,j-1) ) )
        s12(i,j,k)=0.5*(                                                         &
                   ( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
                    +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                          -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                         )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )        &
                  +( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
                    +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                          -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                         )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )    )
      enddo
      enddo       
    ENDDO

    ! now, dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
      enddo
      enddo
    ENDDO
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=2,nk
      do j=1,nj
      do i=1,ni+1
        s13(i,j,k)=0.5*(                                                              &
                   (ua(i,j,k)-ua(i,j,k-1))*rds(k)                                     &
                  +(wa(i,j,k)*rgz(i,j)-wa(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)              &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i-1,j,k  ))       &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i-1,j,k-1)) )     &
                             *(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)                         )
        s13(i,j,k)=s13(i,j,k)*0.5*( gz(i-1,j)*rf(i-1,j,k)+gz(i,j)*rf(i,j,k) )
      enddo
      enddo
      do j=1,nj+1   
      do i=1,ni
        s23(i,j,k)=0.5*(                                                              &
                   (va(i,j,k)-va(i,j,k-1))*rds(k)                                     &
                  +(wa(i,j,k)*rgz(i,j)-wa(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)              &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i,j-1,k  ))       &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i,j-1,k-1)) )     &
                             *(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)                         )
        s23(i,j,k)=s23(i,j,k)*0.5*( gz(i,j-1)*rf(i,j-1,k)+gz(i,j)*rf(i,j,k) )
      enddo
      enddo
    ENDDO

  ENDIF

!  end of calculations for terrain
!-------------------------------------------------------------------------------
!  open boundary conditions:

    IF( wbc.eq.2 .or. ebc.eq.2 .or. sbc.eq.2 .or. nbc.eq.2 )THEN
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=0,nj+1
            s11(0,j,k) = s11(1,j,k)
          enddo
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=0,nj+1
            s11(ni+1,j,k) = s11(ni,j,k)
          enddo
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=0,ni+1
            s22(i,0,k) = s22(i,1,k)
          enddo
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=0,ni+1
            s22(i,nj+1,k) = s22(i,nj,k)
          enddo
          enddo
        ENDIF
        !-----
        !ccccc
        !-----
        IF( wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj+1
            s12(1,j,k) = s12(2,j,k)
          enddo
          enddo
        ENDIF
        IF( ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(j,k)
          do k=1,nk
          do j=1,nj+1
            s12(ni+1,j,k) = s12(ni,j,k)
          enddo
          enddo
        ENDIF
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni+1
            s12(i,1,k) = s12(i,2,k)
          enddo
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(i,k)
          do k=1,nk
          do i=1,ni+1
            s12(i,nj+1,k) = s12(i,nj,k)
          enddo
          enddo
        ENDIF
        !-----
        ! corner points:
        !-----
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            s12(1,1,k) = s12(2,2,k)
          enddo
        ENDIF
        IF( sbc.eq.2 .and. ibs.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            s12(ni+1,1,k) = s12(ni,2,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            wbc.eq.2 .and. ibw.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            s12(1,nj+1,k) = s12(2,nj,k)
          enddo
        ENDIF
        IF( nbc.eq.2 .and. ibn.eq.1 .and. &
            ebc.eq.2 .and. ibe.eq.1 )THEN
          !$omp parallel do default(shared)   &
          !$omp private(k)
          do k=1,nk
            s12(ni+1,nj+1,k) = s12(ni,nj,k)
          enddo
        ENDIF
        !-----
    ENDIF

!----------------------------------------------------------------------
!  if l_h or l_v is zero, set appropriate terms to zero:
!    (just to be sure)

!    IF( horizturb.eq.1 .or. ipbl.eq.2 )THEN
!    IF( l_h*lhref1*lhref2.lt.1.0e-12 )THEN
!!$omp parallel do default(shared)   &
!!$omp private(i,j,k)
!      do k=0,nk+1
!      do j=0,nj+1
!      do i=0,ni+1
!        s11(i,j,k) = 0.0
!        s12(i,j,k) = 0.0
!        s33(i,j,k) = 0.0
!        s22(i,j,k) = 0.0
!      enddo
!      enddo
!      enddo
!    ENDIF
!    ENDIF

!    IF( horizturb.eq.1 .or. ipbl.eq.2 )THEN
!    IF( l_inf.lt.tsmall )THEN
!!$omp parallel do default(shared)   &
!!$omp private(i,j,k)
!      do k=0,nk+1
!      do j=0,nj+1
!      do i=0,ni+1
!        s13(i,j,k) = 0.0
!        s23(i,j,k) = 0.0
!      enddo
!      enddo
!      enddo
!    ENDIF
!    ENDIF

!--------------------------------------------------------------

    IF( axisymm.eq.1 )THEN
      ! lateral boundary condition:
!$omp parallel do default(shared)   &
!$omp private(k)
      do k=0,nk+1
        s13(1,1,k)=0.0
      enddo
    ENDIF

!----------------------------------------------------------------------
!  calculate deformation:
!  Note:  deformation is defined at w points.

    IF(axisymm.eq.0)THEN
      ! Cartesian domain:

      ! Def = 2.0 * S_ij * S_ij
      !
      !     = 2.0 * (  S11*S11 + S12*S12 + S13*S13 
      !              + S21*S21 + S22*S22 + S23*S23 
      !              + S31*S31 + S32*S32 + S33*S33 )
      !
      !     =   2.0*( S11*S11 + S22*S22 + S33*S33 )
      !       + 4.0*( S12*S12 + S13*S13 + S23*S23 )

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,rrf)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=( c1(i,j,k)*s11(i,j,k-1)**2 + c2(i,j,k)*s11(i,j,k)**2 )
        tmp22=( c1(i,j,k)*s22(i,j,k-1)**2 + c2(i,j,k)*s22(i,j,k)**2 )
        tmp33=( c1(i,j,k)*s33(i,j,k-1)**2 + c2(i,j,k)*s33(i,j,k)**2 )

        tmp12=0.25*( c1(i,j,k)*( ( s12(i,j  ,k-1)**2 + s12(i+1,j+1,k-1)**2 )     &
                               + ( s12(i,j+1,k-1)**2 + s12(i+1,j  ,k-1)**2 ) )   &
                    +c2(i,j,k)*( ( s12(i,j  ,k  )**2 + s12(i+1,j+1,k  )**2 )     &
                               + ( s12(i,j+1,k  )**2 + s12(i+1,j  ,k  )**2 ) ) )

        tmp13=0.5*( s13(i,j,k)**2 + s13(i+1,j,k)**2 )

        tmp23=0.5*( s23(i,j,k)**2 + s23(i,j+1,k)**2 )

        rrf = 1.0/(rf(i,j,k)**2)

        defv(i,j,k)= 4.0*( tmp13 + tmp23 )*rrf

        defh(i,j,k) = ( 2.0*( ( tmp11 + tmp22 ) + tmp33 ) + 4.0*tmp12 )*rrf

      enddo
      enddo
      enddo

!--------------------------------------------
    ELSE
      ! axisymmetric domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp33,tmp12,tmp13,tmp23,rrf,r1,r2,r3,r4)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=( c1(1,1,k)*(s11(i,j,k-1)**2) + c2(1,1,k)*(s11(i,j,k)**2) )
        tmp33=( c1(1,1,k)*(s33(i,j,k-1)**2) + c2(1,1,k)*(s33(i,j,k)**2) )

        tmp12=0.5*(  c1(1,1,k)*( s12(i,j  ,k-1)**2 + s12(i+1,j  ,k-1)**2 )     &
                   + c2(1,1,k)*( s12(i,j  ,k  )**2 + s12(i+1,j  ,k  )**2 ) )

        tmp13=0.5*( s13(i,j,k)**2 + s13(i+1,j,k)**2 )

        tmp23=      s23(i,j,k)**2

        rrf = 1.0/(rf(i,j,k)**2)

        defv(i,j,k)= 4.0*( tmp13 + tmp23 )*rrf

        defh(i,j,k) = ( 2.0*( tmp11 + tmp33 ) + 4.0*tmp12 )*rrf

      enddo
      enddo
      enddo

    ENDIF  ! endif for axisymm


!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine calcdef
      ! calcdef


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,nm,t,qt,thv,cloud,rh,qvci,   &
                        prs,pp,th,qa,iamsat)
      use input
      use constants
      use cm1libs , only : rslf,rsif
      implicit none

      logical, intent(in), dimension(maxq) :: cloudvar
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t,qt,thv,cloud,rh,qvci
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: prs,pp,th
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: iamsat

      integer :: i,j,k,n
      real :: pavg,tavg,qtavg,qvs,lh,cpml,gamma,drdt
      real :: nml,nmi,ff,qvsl,qvsi,qvciavg,qv,ql,qi

      real, parameter :: nmsat = 1.0e-10

!----------------------------------------------------------------------
!  Dry nm

    IF(imoist.eq.0)then

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        nm(i,j,k)=alog( (th0(i,j,k)+th(i,j,k))/(th0(i,j,k-1)+th(i,j,k-1)) ) &
                    *g*rdz*mf(i,j,k)
      enddo
      enddo
      enddo

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo
      enddo

!-----------------------------------------------------------------------
!  Moist nm

    ELSE

    if( ptype.eq.0 )then
      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        qvci(i,j,k) = qa(i,j,k,nqv)
        qt(i,j,k) = 0.0
      enddo
      enddo
      enddo
    else
      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        qvci(i,j,k) = qa(i,j,k,nqv) + qa(i,j,k,nqc)
        qt(i,j,k) = 0.0
      enddo
      enddo
      enddo
    endif

      IF( iice.eq.1 )THEN
      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        qvci(i,j,k) = qvci(i,j,k) + qa(i,j,k,nqi)
      enddo
      enddo
      enddo
      ENDIF

      DO n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
          enddo
        ENDIF
      ENDDO

      k = 1

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        t(i,j,k)=(th0(i,j,k)+th(i,j,k))*(pi0(i,j,k)+pp(i,j,k))
        thv(i,j,k)=(th0(i,j,k)+th(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                         /(1.0+qt(i,j,k))
      enddo
      enddo

    icecheck:  &
    IF( iice.eq.1 )THEN

      !$omp parallel do default(shared)   &
      !$omp private(i,j,pavg,tavg,qtavg,qvs,lh,cpml,gamma,drdt,nml,nmi,ff,qvsl,qvsi,qvciavg,qv,ql,qi)
      do j=1,nj
      do k=2,nk
      do i=1,ni
        t(i,j,k)=(th0(i,j,k)+th(i,j,k))*(pi0(i,j,k)+pp(i,j,k))
        thv(i,j,k)=(th0(i,j,k)+th(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                         /(1.0+qt(i,j,k))
        ! subsaturated formulation:
        nm(i,j,k)=g*alog(thv(i,j,k)/thv(i,j,k-1))*rdz*mf(i,j,k)

        ! saturated formulation (if necessary):
        pavg = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
        tavg =   c1(i,j,k)*t(i,j,k-1)+  c2(i,j,k)*t(i,j,k)
        qtavg=  c1(i,j,k)*qt(i,j,k-1)+ c2(i,j,k)*qt(i,j,k)
        qvciavg = c1(i,j,k)*qvci(i,j,k-1)+ c2(i,j,k)*qvci(i,j,k)

        if( tavg.ge.273.15 )then
          qvs = rslf(pavg,tavg)
          ql = max( qvciavg - qvs , 0.0 )
          if( ql.gt.nmsat )then
            iamsat(i,j,k) = .true.
            qv = max( qvciavg - ql , 0.0 )
            cpml=cp+cpv*qv+cpl*ql
            lh=lv1-lv2*tavg
            drdt=17.67*(273.15-29.65)*qvs/((tavg-29.65)**2)
            gamma=g*(1.0+qtavg)*(1.0+lh*qvs/(rd*tavg))/(cpml+lh*drdt)
            nm(i,j,k)=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                              +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
                           -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
          endif
        elseif( tavg.le.233.15 )then
          qvs = rsif(pavg,tavg)
          qi = max( qvciavg - qvs , 0.0 )
          if( qi.gt.nmsat )then
            iamsat(i,j,k) = .true.
            qv = max( qvciavg - qi , 0.0 )
            cpml=cp+cpv*qv+cpi*qi
            lh=ls1-ls2*tavg
            drdt=21.8745584*(273.15-7.66)*qvs/((tavg-7.66)**2)
            gamma=g*(1.0+qtavg)*(1.0+lh*qvs/(rd*tavg))/(cpml+lh*drdt)
            nm(i,j,k)=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                              +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
                           -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
          endif
        else
          ff = (tavg-233.15)/(273.15-233.15)
          qvsl = rslf(pavg,tavg)
          qvsi = rsif(pavg,tavg)
          qvs = ff*qvsl + (1.0-ff)*qvsi
          qi = max( 0.0 , (1.0-ff)*(qvciavg-qvs) )
          ql = max( 0.0 , ff*(qvciavg-qvs) )
          if( (ql+qi).gt.nmsat )then
            iamsat(i,j,k) = .true.
            qv = max( 0.0 , qvciavg - ql - qi )
            cpml=cp+cpv*qv+cpl*ql+cpi*qi

            lh=lv1-lv2*tavg
            drdt=17.67*(273.15-29.65)*qvsl/((tavg-29.65)**2)
            gamma=g*(1.0+qtavg)*(1.0+lh*qvsl/(rd*tavg))/(cpml+lh*drdt)
            nml=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                              +gamma/tavg )*(1.0+tavg*drdt/(eps+qvsl))   &
                           -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )

            lh=ls1-ls2*tavg
            drdt=21.8745584*(273.15-7.66)*qvsi/((tavg-7.66)**2)
            gamma=g*(1.0+qtavg)*(1.0+lh*qvsi/(rd*tavg))/(cpml+lh*drdt)
            nmi=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                              +gamma/tavg )*(1.0+tavg*drdt/(eps+qvsi))   &
                           -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )

            nm(i,j,k) = ff*nml + (1.0-ff)*nmi
          endif
        endif
      enddo
      enddo
      enddo

    ELSE  icecheck

      ! liquid only:
      !$omp parallel do default(shared)   &
      !$omp private(i,j,pavg,tavg,qtavg,qvciavg,qvs,ql,qv,cpml,lh,drdt,gamma)
      do j=1,nj
      do k=2,nk
      do i=1,ni
        t(i,j,k)=(th0(i,j,k)+th(i,j,k))*(pi0(i,j,k)+pp(i,j,k))
        thv(i,j,k)=(th0(i,j,k)+th(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                         /(1.0+qt(i,j,k))
        ! subsaturated formulation:
        nm(i,j,k)=g*alog(thv(i,j,k)/thv(i,j,k-1))*rdz*mf(i,j,k)

        ! saturated formulation (if necessary):
        pavg = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
        tavg =   c1(i,j,k)*t(i,j,k-1)+  c2(i,j,k)*t(i,j,k)
        qtavg=  c1(i,j,k)*qt(i,j,k-1)+ c2(i,j,k)*qt(i,j,k)
        qvciavg = c1(i,j,k)*qvci(i,j,k-1)+ c2(i,j,k)*qvci(i,j,k)

          qvs = rslf(pavg,tavg)
          ql = max( qvciavg - qvs , 0.0 )
          if( ql.gt.nmsat )then
            iamsat(i,j,k) = .true.
            qv = max( qvciavg - ql , 0.0 )
            cpml=cp+cpv*qv+cpl*ql
            lh=lv1-lv2*tavg
            drdt=17.67*(273.15-29.65)*qvs/((tavg-29.65)**2)
            gamma=g*(1.0+qtavg)*(1.0+lh*qvs/(rd*tavg))/(cpml+lh*drdt)
            nm(i,j,k)=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                              +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
                           -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
          endif
      enddo
      enddo
      enddo

    ENDIF  icecheck

      !$omp parallel do default(shared)  &
      !$omp private(i,j)
      do j=1,nj
      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo
      enddo

    ENDIF    ! endif for imoist

!----------------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end subroutine calcnm


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! t2psmm-1
      subroutine  t2psmm(dt,rtime,xf,rxf,c1,c2,                                   &
                         zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,                     &
                         uavg,vavg,savg,l2p,kmw,ufw,vfw,gamk,gamwall,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,                        &
                         dum1,dum2,s13 ,s23 ,up  ,vp  ,fru ,frv ,                 &
                         ust,zntmp,stau,mol,rho,rr,rf,kmv,kmh,                    &
                         ugr ,vgr ,uf,vf,arf1,arf2,                               &
                         u1,v1,s1,u1b,v1b,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,  &
                         reqs_s,uw31,uw32,ue31,ue32,us31,us32,un31,un32)

      use input
      use constants , only : karman
      use sfcphys_module , only : stabil_funcs
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in) :: dt,rtime
      real, intent(in), dimension(ib:ie+1) :: xf,rxf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,zh,mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf0,rr0,rho0,rrf0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(kb:ke) :: uavg,vavg,savg,l2p,kmw,ufw,vfw,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4
      real, intent(in), dimension(kb:ke) :: gamwall
      double precision, intent(inout), dimension(kb:ke,3+numq) :: cavg
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ugr
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vgr
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,s13,s23,up,vp,fru,frv
      real, intent(in), dimension(ib:ie,jb:je) :: ust,zntmp
      real, intent(inout), dimension(ib:ie,jb:je) :: stau
      real, intent(in), dimension(ibl:iel,jbl:jel) :: mol
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv,kmh
      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie+1) :: arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,s1
      real, intent(inout), dimension(kb:ke) :: u1b,v1b
      double precision, intent(inout) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv
      integer, intent(inout), dimension(rmp) :: reqs_s
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32

      integer :: i,j,k,n
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,wbar,rratio,fstar,tem,temx,temy
      real :: d2pm1,d2pm2,d2pm3,d2pm4,d2pm5,ueo,veo,dsdz,rdsdz,dsdzu,dsdzv,dudz,dvdz,s13b,s23b
      double precision :: temd,teme,savgt,t13avg,t23avg,d13avg,d23avg,s13avg,s23avg,z13avg,z23avg,  &
                          ustbar,zntbar,molbar,kavgt,rustbar,ubar,vbar,ufrc,vfrc,ravg,rfavg,rflast,spavg
      real :: tpsim,tpsih,tphim,tphih,zeta

      real, dimension(kb:ke) :: wspa
      double precision, dimension(kb:ke) :: shravg

    !-------------------------------------------

    doingt2p:  &
    if( t2pfac.ge.0.001 )then

        kmw = 0.0

    !-------------------------------------------
    !  Get domain average profiles:

      uavg = 0.0
      vavg = 0.0

      do k=1,nk
        !----
        cavg(k,1) = 0.0
        cavg(k,2) = 0.0
        !----
        do j=1,nj
        do i=1,ni
          cavg(k,1) = cavg(k,1) + ugr(i,j,k)
          cavg(k,2) = cavg(k,2) + vgr(i,j,k)
        enddo
        enddo
        !----
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,cavg(1,1),nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,cavg(1,2),nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      temd = 1.0/dble(nx*ny)

      do k=1,nk
        uavg(k)  = cavg(k,1)*temd
        vavg(k)  = cavg(k,2)*temd
        savg(k)  = max( sqrt( uavg(k)**2 + vavg(k)**2 ) , 1.0e-10 )
        u1b(k) = uavg(k)
        v1b(k) = vavg(k)
      enddo

    !-------------------------------------------
    !  surface vars:

            call sfcstress(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,rf,mf,dum1,dum2,  &
                        kmh,kmv,s13,s23,ua,ugr,va,vgr,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv)

        teme = 1.0/dble(nx*ny)

      IF( use_avg_sfc )THEN

        ustbar = ust(1,1)
        zntbar = zntmp(1,1)
        molbar = mol(1,1)

      ELSE

        ustbar = 0.0
        zntbar = 0.0
        molbar = 0.0

        do j=1,nj
        do i=1,ni
          ustbar = ustbar + ust(i,j)
          zntbar = zntbar + zntmp(i,j)
          molbar = molbar + mol(i,j)
        enddo
        enddo

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,ustbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zntbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,molbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        ustbar = ustbar*teme
        zntbar = zntbar*teme
        molbar = molbar*teme

      ENDIF

!!!        if(myid.eq.0) print *,'  zntbar,ustbar,ust = ',zntbar,ustbar,ust(1,1)

        do k=1,(ntwk+1)
!!!          wspa(k) = (ustbar/karman)*alog( sngl(zh(1,1,k)/zntbar) )
          if( abs(molbar).gt.1.0e-6 )then
            zeta = zh(1,1,k)/molbar
            call stabil_funcs(zeta,tphim,tphih,tpsim,tpsih)
          else
            tpsim = 0.0
          endif
          wspa(k) = (ustbar/karman)*( alog(zh(1,1,k)/sngl(zntbar)) - tpsim )
        enddo

    !-------------------------------------------
    !  Get new value of gamma:

    do k=2,ntwk

      spavg = 0.0
      s13b = 0.5*(uavg(k)-uavg(k-1))*rdz*mf(1,1,k)
      s23b = 0.5*(vavg(k)-vavg(k-1))*rdz*mf(1,1,k)

      do j=1,nj
      do i=1,ni

        tmp11=( c1(i,j,k)*(t11(i,j,k-1)**2) + c2(i,j,k)*(t11(i,j,k)**2) )
        tmp22=( c1(i,j,k)*(t22(i,j,k-1)**2) + c2(i,j,k)*(t22(i,j,k)**2) )
        tmp33=( c1(i,j,k)*(t33(i,j,k-1)**2) + c2(i,j,k)*(t33(i,j,k)**2) )

        tmp12=0.25*( c1(i,j,k)*( ( (t12(i  ,j  ,k-1)**2)     &
                                 + (t12(i+1,j+1,k-1)**2) )   &
                               + ( (t12(i  ,j+1,k-1)**2)     &
                                 + (t12(i+1,j  ,k-1)**2) ) ) &
                    +c2(i,j,k)*( ( (t12(i  ,j  ,k  )**2)     &
                                 + (t12(i+1,j+1,k  )**2) )   &
                               + ( (t12(i  ,j+1,k  )**2)     &
                                 + (t12(i+1,j  ,k  )**2) ) ) )

        tmp13=0.5*( (t13(i,j,k)-s13b)**2 + (t13(i+1,j,k)-s13b)**2 )

        tmp23=0.5*( (t23(i,j,k)-s23b)**2 + (t23(i,j+1,k)-s23b)**2 )

        spavg = spavg + ( 4.0*( tmp12 + ( tmp13 + tmp23 ) ) + 2.0*( ( tmp11 + tmp22 ) + tmp33 ) )/(rf(i,j,k)**2)

      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,spavg,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      spavg = sqrt( spavg*teme )

      shravg(k) = sqrt( 4.0*( s13b**2 + s23b**2 ) )

      gamk(k) = spavg/( smeps + (spavg+gamwall(k)*shravg(k)) )

      s2p(k) = spavg
      s2b(k) = shravg(k)

    enddo

    gamk(1) = gamk(2)

    !-------------------------------------------


      rustbar = 1.0/max( 0.0001 , ustbar)

      if( abs(molbar).gt.1.0e-6 )then
        zeta = zf(1,1,2)/molbar
        call stabil_funcs(zeta,tphim,tphih,tpsim,tpsih)
      else
        tpsim = 1.0
      endif

      kloops94:  &
      DO k = 2 , 2

        ubar = 0.0
        vbar = 0.0
        t13avg = 0.0
        t23avg = 0.0
        d13avg = 0.0
        d23avg = 0.0
        s13avg = 0.0
        s23avg = 0.0
        kavgt = 0.0

        temx = 0.0
        temy = 0.0

        ! Note:  here, we assume wavg = 0

      IF( k.eq.2 )THEN
        ! Get avg u,v at w level 2:
        ! (use same algorithms as in adv_routines)
        do j=1,nj
        do i=1,ni
          !-----
          wbar = 0.5*(wa(i,j,k)+wa(i-1,j,k))
          IF( wbar.ge.0.0 )THEN
            up(i,j,k) = ( c1(1,1,k)*ua(i,j,k-1) &
                         +c2(1,1,k)*ua(i,j,k  ) )
          ELSE
            up(i,j,k) = flx3(ua(i,j,k+1),ua(i,j,k),ua(i,j,k-1))
          ENDIF
          ubar = ubar+up(i,j,k)
          !-----
          wbar = 0.5*(wa(i,j,k)+wa(i,j-1,k))
          IF( wbar.ge.0.0 )THEN
            vp(i,j,k) = ( c1(1,1,k)*va(i,j,k-1) &
                         +c2(1,1,k)*va(i,j,k  ) )
          ELSE
            vp(i,j,k) = flx3(va(i,j,k+1),va(i,j,k),va(i,j,k-1))
          ENDIF
          vbar = vbar+vp(i,j,k)
          !-----
        enddo
        enddo
      ELSE
        stop 87632
      ENDIF

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,ubar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,vbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        ! avg u,v at w level 2:
        ubar = ubar*teme
        vbar = vbar*teme

      IF( k.eq.2 )THEN
        ! Get <u'w'> and <v'w'> at w level 2:
        ! (use same algorithms as in adv_routines)
        do j=1,nj
        do i=1,ni
          !-----
!          ueo = ( c1(1,1,k)*ua(i,j,k-1)+c2(1,1,k)*ua(i,j,k  ) )
!          wbar = rf0(1,1,k)*0.5*(wa(i,j,k)+wa(i-1,j,k))
!          fru(i,j,k) = wbar*(ueo-ubar)
!          IF( wbar.lt.0.0 )THEN
!            fru(i,j,k) = fru(i,j,k) + wbar*(up(i,j,k)-ueo)
!          ENDIF
!          !-----
!          veo = ( c1(1,1,k)*va(i,j,k-1)+c2(1,1,k)*va(i,j,k  ) )
!          wbar = rf0(1,1,k)*0.5*(wa(i,j,k)+wa(i,j-1,k))
!          frv(i,j,k) = wbar*(veo-vbar)
!          IF( wbar.lt.0.0 )THEN
!            frv(i,j,k) = frv(i,j,k) + wbar*(vp(i,j,k)-veo)
!          ENDIF
          !-----
          fru(i,j,k) = 0.5*(wa(i,j,k)+wa(i-1,j,k))*(up(i,j,k)-ubar)
          frv(i,j,k) = 0.5*(wa(i,j,k)+wa(i,j-1,k))*(vp(i,j,k)-vbar)
          !-----
          t13avg = t13avg + fru(i,j,k)
          t23avg = t23avg + frv(i,j,k)
        enddo
        enddo
      ELSE
        stop 32987
      ENDIF

        do j=1,nj
        do i=1,ni
          kavgt = kavgt + gamk(k)*kmv(i,j,k)
!!!          kavgt = kavgt + kmv(i,j,k)
        enddo
        enddo

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,t13avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t23avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,kavgt ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        t13avg = t13avg*teme
        t23avg = t23avg*teme
        kavgt = kavgt*teme

        !-----
        ! Get kmw:

        rdsdz = 1.0/max( 1.0e-20 , (wspa(k)-wspa(k-1))*rdz*mf(1,1,k) )
!!!        rdsdz = karman*zf(1,1,k)/(ustbar*tphim)
!!!        rdsdz = karman*zf(1,1,k)/ustbar

        d2pm1 = -kavgt
        d2pm2 = -rdsdz*sqrt( t13avg**2 + t23avg**2 )
        d2pm3 =  0.0
        d2pm4 =  0.0

!!!        kmw(k) = t2pfac*max( 0.0 , ustbar*karman*zf(1,1,2) + d2pm1 + d2pm2 )
!!!        kmw(k) = t2pfac*max( 0.0 , ustbar*ustbar*rdsdz + d2pm1 + d2pm2 )
        ! use discretized analytic shear:
        kmw(k) = t2pfac*max( 0.0 , ustbar*ustbar*rdsdz + d2pm1 + d2pm2 )

        kmw(k) = min( kmw(k) , ustbar*karman*zf(1,1,2) )

        t2pm1(k) = d2pm2
        t2pm2(k) = d2pm1
!!!        t2pm3(k) = ustbar*karman*zf(1,1,2)
        t2pm3(k) = ustbar*ustbar*rdsdz

      ENDDO  kloops94

        l2p(1) = zntbar
        l2p(2) = sqrt( kmw(2)*karman*zf(1,1,2)*rustbar )

        do k=3,ntwk
          l2p(k) = l2p(2)
        enddo

        do k=3,ntwk
!!!          kmw(k) = gamwall(k)*l2p(k)*l2p(k)*shravg(k)
!!!          kmw(k) = l2p(k)*l2p(k)*shravg(k)
          ! use analytic shear::
          kmw(k) = gamwall(k)*l2p(k)*l2p(k)*ustbar/(karman*zf(1,1,k))
        enddo

        do k=2,ntwk
          ufw(k) = kmw(k)*(u1b(k)-u1b(k-1))*rdz*mf(1,1,k)
          vfw(k) = kmw(k)*(v1b(k)-v1b(k-1))*rdz*mf(1,1,k)
        enddo

    !-------------------------------------------

    endif  doingt2p

      if(timestats.ge.1) time_turb=time_turb+mytime()

      end subroutine t2psmm
! t2psmm-2


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! t2pcode-1
      subroutine t2pcode(dt,rtime,xf,rxf,c1,c2,zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,            &
                         uavg,vavg,savg,l2p,kmw,ufw,vfw,gamwall,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,dum1,dum2,s13 ,s23 ,                 &
                         ust,zntmp,stau,mol,rho,rr,rf,kmh,kmv,                                 &
                         ugr ,fsu  ,u2pt,vgr ,fsv  ,v2pt,uf,vf,arf1,arf2,                      &
                         u1,v1,s1,u1b,v1b,                                                     &
                         avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,                            &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)

      use input
      use constants , only : karman
      use sfcphys_module , only : stabil_funcs
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in) :: dt,rtime
      real, intent(in), dimension(ib:ie+1) :: xf,rxf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,zh,mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf0,rr0,rho0,rrf0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(kb:ke) :: uavg,vavg,savg,l2p,kmw,ufw,vfw,t2pm1,t2pm2,t2pm3,t2pm4
      real, intent(in), dimension(kb:ke) :: gamwall
      real, intent(inout), dimension(kb:ke) :: gamk,s2p,s2b
      double precision, intent(inout), dimension(kb:ke,3+numq) :: cavg
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ugr,fsu
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vgr,fsv
      real, intent(in), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt,v2pt
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,s13,s23
      real, intent(in), dimension(ib:ie,jb:je) :: ust,zntmp
      real, intent(inout), dimension(ib:ie,jb:je) :: stau
      real, intent(in), dimension(ibl:iel,jbl:jel) :: mol
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie+1) :: arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,s1
      real, intent(inout), dimension(kb:ke) :: u1b,v1b
      double precision, intent(in) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      integer, intent(inout), dimension(rmp) :: reqs_s

      integer :: i,j,k,n
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,wbar,rratio,fstar,tem,temx,temy
      real :: d2pm1,d2pm2,d2pm3,d2pm4,d2pm5,ueo,veo,dsdz,dsdzu,dsdzv,dudz,dvdz,up,vp,s13b,s23b
      double precision :: temd,teme,savgt,t13avg,t23avg,  &
                          ustbar,zntbar,molbar,kavgt,rustbar,ubar,vbar,ufrc,vfrc,ravg,rfavg,rflast,shravg,spavg
      real :: tpsim,tpsih,tphim,tphih,zeta

      real, dimension(kb:ke) :: wspa

    !-------------------------------------------

      kmw = 0.0
      ufw = 0.0
      vfw = 0.0
      l2p = 0.0
      wspa = 0.0

    doingt2p:  &
    if( t2pfac.ge.0.001 )then

      if( myid.eq.0 ) print *,'  t2pfac = ',t2pfac

    !-------------------------------------------

        call     sfcstress(xf,rxf,arf1,arf2,ust,stau,u1,v1,s1,rf,mf,dum1,dum2,  &
                        kmh,kmv,s13,s23,ua,ugr,va,vgr,avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv)

    !-------------------------------------------
    !  surface vars:

        teme = 1.0/dble(nx*ny)

      IF( use_avg_sfc )THEN

        ustbar = ust(1,1)
        zntbar = zntmp(1,1)
        molbar = mol(1,1)

      ELSE

        ustbar = 0.0
        zntbar = 0.0
        molbar = 0.0

        do j=1,nj
        do i=1,ni
          ustbar = ustbar + ust(i,j)
          zntbar = zntbar + zntmp(i,j)
          molbar = molbar + mol(i,j)
        enddo
        enddo

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,ustbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zntbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,molbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        ustbar = ustbar*teme
        zntbar = zntbar*teme
        molbar = molbar*teme

      ENDIF

        do k=1,(ntwk+1)
!!!          wspa(k) = (ustbar/karman)*alog( sngl(zh(1,1,k)/zntbar) )
          if( abs(molbar).gt.1.0e-6 )then
            zeta = zh(1,1,k)/molbar
            call stabil_funcs(zeta,tphim,tphih,tpsim,tpsih)
          else
            tpsim = 0.0
          endif
          wspa(k) = (ustbar/karman)*( alog(zh(1,1,k)/sngl(zntbar)) - tpsim )
        enddo

    !-------------------------------------------

      ! new formulation (accounts for all forcing terms in u,v equations)

      d2pm2 =  0.0

      kloop:  &
      DO k=1,(ntwk-1)

        t13avg = 0.0
        t23avg = 0.0
        rfavg = 0.0

        do j=1,nj
        do i=1,ni
          t13avg = t13avg + u2pt(i,j,k)
          t23avg = t23avg + v2pt(i,j,k)
          rfavg = rfavg + rf(i,j,k+1)*rr(i,j,k)
        enddo
        enddo

#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,t13avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,t23avg,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,rfavg ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        t13avg = t13avg*teme
        t23avg = t23avg*teme
        rfavg = rfavg*teme

    !-------------------------------------------
    !  Get kmw:

        d2pm1 =  ( uavg(k)*t13avg + vavg(k)*t23avg )/savg(k)
        d2pm3 =  ( uavg(k)*ufw(k) + vavg(k)*vfw(k) )*rfavg*rdz*mh(1,1,k)/savg(k)

        dsdz = (wspa(k+1)-wspa(k))*rdz*mf(1,1,k+1)
        tem = t2pfac*gamwall(k+1)/( rdz*mh(1,1,k)*max(1.0e-12,dsdz)*rfavg )

        ! 180903:
        kmw(k+1) = max( 0.0 , ( 0.0 - d2pm1 + d2pm3 )*tem )

        ufw(k+1) = kmw(k+1)*(uavg(k+1)-uavg(k))*rdz*mf(1,1,k+1)
        vfw(k+1) = kmw(k+1)*(vavg(k+1)-vavg(k))*rdz*mf(1,1,k+1)

        t2pm1(k+1) = -d2pm1*tem
        t2pm3(k+1) =  d2pm3*tem

      ENDDO  kloop

        ! diagnostic only:
        l2p(1) = zntbar
        do k=2,ntwk
          l2p(k) = sqrt( kmw(k)*karman*(zf(1,1,k)-zntbar)/max(0.0001,ustbar) )
        enddo

    !-------------------------------------------

    endif  doingt2p

      if(timestats.ge.1) time_turb=time_turb+mytime()

      end subroutine t2pcode
! t2pcode-2


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! t2pcode-1
      subroutine t2pcodetavg(dt,rtime,xf,rxf,c1,c2,zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,            &
                         uavg,vavg,savx,l2p,kmw,ufw,vfw,gamwall,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,dum1,dum2,s13 ,wspa,                 &
                         ust,zntmp,stau,mol,s10,rho,rr,rf,kmh,kmv,                             &
                         ugr ,fsu  ,u2pt,vgr ,fsv  ,v2pt,uf,vf,arf1,arf2,kmwk,ufwk,vfwk,       &
                         u1,v1,s1,u1b,v1b,                                                     &
                         avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv,                            &
                         timavg,sfctimavg,                                                     &
                         uw31,uw32,ue31,ue32,us31,us32,un31,un32,                              &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)

      use input
      use constants , only : karman
      use sfcphys_module , only : stabil_funcs
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in) :: dt,rtime
      real, intent(in), dimension(ib:ie+1) :: xf,rxf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,zh,mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf0,rr0,rho0,rrf0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(kb:ke) :: uavg,vavg,savx,l2p,kmw,ufw,vfw,t2pm1,t2pm2,t2pm3,t2pm4
      real, intent(in), dimension(kb:ke) :: gamwall
      real, intent(inout), dimension(kb:ke) :: gamk,s2p,s2b
      double precision, intent(inout), dimension(kb:ke,3+numq) :: cavg
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ugr,fsu
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vgr,fsv
      real, intent(in), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt,v2pt
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: kmwk,ufwk,vfwk
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,s13,wspa
      real, intent(in), dimension(ib:ie,jb:je) :: ust,zntmp,s10
      real, intent(inout), dimension(ib:ie,jb:je) :: stau
      real, intent(in), dimension(ibl:iel,jbl:jel) :: mol
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie+1) :: arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,s1
      real, intent(inout), dimension(kb:ke) :: u1b,v1b
      double precision, intent(in) :: avgsfcu,avgsfcv,avgsfcs,avgsfcsu,avgsfcsv
      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      real, intent(in), dimension(ibta:ieta,jbta:jeta,nsfctavr) :: sfctimavg
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      integer, intent(inout), dimension(rmp) :: reqs_s

      integer :: i,j,k,n
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,wbar,rratio,fstar,tem,temx,temy,savg
      real :: d2pm1,d2pm2,d2pm3,d2pm4,d2pm5,ueo,veo,dsdz,dsdzu,dsdzv,dudz,dvdz,up,vp,s13b,s23b
      real :: uu,vv,uup1,vvp1
      double precision :: temd,teme,savgt,t13avg,t23avg,d13avg,d23avg,s13avg,s23avg,z13avg,z23avg,  &
                          ustbar,zntbar,molbar,kavgt,rustbar,ubar,vbar,ufrc,vfrc,ravg,rfavg,rflast,shravg,spavg
      real :: tpsim,tpsih,tphim,tphih,zeta

#ifdef MPI
      integer, dimension(rmp) :: req1,req2
#endif

!!!      real, dimension(kb:ke) :: wspa

    !-------------------------------------------

      kmwk = 0.0
      ufwk = 0.0
      vfwk = 0.0
!!!      wspa = 0.0

    doingt2p:  &
    if( t2pfac.ge.0.001 )then

      if( myid.eq.0 ) print *,'  t2pfac = ',t2pfac

        do k=1,(ntwk+1)
        do j=1,nj
        do i=1,ni
          ustbar = sfctimavg(i,j,1)
          zntbar = sfctimavg(i,j,2)
          molbar = sfctimavg(i,j,3)
          if( abs(molbar).gt.1.0e-6 )then
            zeta = zh(i,j,k)/molbar
            call stabil_funcs(zeta,tphim,tphih,tpsim,tpsih)
          else
            tpsim = 0.0
          endif
          wspa(i,j,k) = (ustbar/karman)*( alog(zh(i,j,k)/sngl(zntbar)) - tpsim )
        enddo
        enddo
        enddo

    !-------------------------------------------

      ! new formulation (accounts for all forcing terms in u,v equations)

    kloop:  &
    DO k=1,(ntwk-1)

      do j=1,nj
      do i=1,ni

    !-------------------------------------------
    !  Get kmw:

        uu = 0.5*(timavg(i,j,k,utav)+timavg(i+1,j,k,utav))
        vv = 0.5*(timavg(i,j,k,vtav)+timavg(i,j+1,k,vtav))

        savg = 1.0/max( sqrt( uu**2 + vv**2 ) , 1.0e-10 )

        d2pm1 =  0.5*( uu*(timavg(i,j,k,uutav)+timavg(i+1,j,k,uutav))  &
                      +vv*(timavg(i,j,k,vvtav)+timavg(i,j+1,k,vvtav)) )*savg
        d2pm3 =  ( uu*ufwk(i,j,k) + vv*vfwk(i,j,k) )*rf(i,j,k)*rdz*mh(i,j,k)*savg

        dsdz = (wspa(i,j,k+1)-wspa(i,j,k))*rdz*mf(i,j,k+1)
        tem = t2pfac*gamwall(k+1)/( rdz*mh(i,j,k)*max(1.0e-12,dsdz)*rf(i,j,k+1)*rr(i,j,k) )

        ! 180903:
        kmwk(i,j,k+1) = max( 0.0 , ( 0.0 - d2pm1 + d2pm3 )*tem )

        uup1 = 0.5*(timavg(i,j,k+1,utav)+timavg(i+1,j,k+1,utav))
        vvp1 = 0.5*(timavg(i,j,k+1,vtav)+timavg(i,j+1,k+1,vtav))

        ufwk(i,j,k+1) = kmwk(i,j,k+1)*(uup1-uu)*rdz*mf(i,j,k+1)
        vfwk(i,j,k+1) = kmwk(i,j,k+1)*(vvp1-vv)*rdz*mf(i,j,k+1)

        ! GHB, 220824: set kmw to zero for very low wind speeds
        if( s10(i,j).lt.0.1 )then
          kmwk(i,j,k+1) = 0.0
          ufwk(i,j,k+1) = 0.0
          vfwk(i,j,k+1) = 0.0
        endif

      enddo
      enddo

      call bc2d(ufwk(ib2pt,jb2pt,k+1))
#ifdef MPI
      call comm_1s2d_start(ufwk(ib2pt,jb2pt,k+1),uw31(1,1,k),uw32(1,1,k),ue31(1,1,k),ue32(1,1,k),  &
                                                 us31(1,1,k),us32(1,1,k),un31(1,1,k),un32(1,1,k),req1)
#endif
      call bc2d(vfwk(ib2pt,jb2pt,k+1))
#ifdef MPI
      call comm_1s2d_start(vfwk(ib2pt,jb2pt,k+1),sw31(1,1,k),sw32(1,1,k),se31(1,1,k),se32(1,1,k),  &
                                                 ss31(1,1,k),ss32(1,1,k),sn31(1,1,k),sn32(1,1,k),req2)
      call comm_1s2d_end(ufwk(ib2pt,jb2pt,k+1),  uw31(1,1,k),uw32(1,1,k),ue31(1,1,k),ue32(1,1,k),  &
                                                 us31(1,1,k),us32(1,1,k),un31(1,1,k),un32(1,1,k),req1)
      call comm_1s2d_end(vfwk(ib2pt,jb2pt,k+1),  sw31(1,1,k),sw32(1,1,k),se31(1,1,k),se32(1,1,k),  &
                                                 ss31(1,1,k),ss32(1,1,k),sn31(1,1,k),sn32(1,1,k),req2)
#endif

    ENDDO  kloop

    !-------------------------------------------

    endif  doingt2p

      if(timestats.ge.1) time_turb=time_turb+mytime()

      end subroutine t2pcodetavg
! t2pcode-2


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine getgamk(dt,rtime,xf,rxf,c1,c2,zh,mh,zf,mf,rf0,rr0,rho0,rrf0,u0,v0,            &
                         uavg,vavg,savg,l2p,kmw,gamwall,gamk,s2p,s2b,t2pm1,t2pm2,t2pm3,t2pm4,cavg,  &
                         ua,va,wa,t11,t12,t13,t22,t23,t33,dum1,dum2,s13 ,s23 ,                 &
                         ust,zntmp,stau,rho,rr,rf,kmh,kmv,                                     &
                         ugr ,fsu  ,u2pt,vgr ,fsv  ,v2pt,uf,vf,arf1,arf2,                      &
                         u1,v1,s1,u1b,v1b,                                                     &
                         avgsfcu,avgsfcv,avgsfcs,                                              &
                         sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,reqs_s)

      use input
      use constants , only : karman
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in) :: dt,rtime
      real, intent(in), dimension(ib:ie+1) :: xf,rxf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2,zh,mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rf0,rr0,rho0,rrf0
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(kb:ke) :: uavg,vavg,savg,l2p,kmw,t2pm1,t2pm2,t2pm3,t2pm4
      real, intent(in), dimension(kb:ke) :: gamwall
      real, intent(inout), dimension(kb:ke) :: gamk,s2p,s2b
      double precision, intent(inout), dimension(kb:ke,3+numq) :: cavg
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ugr,fsu
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vgr,fsv
      real, intent(in), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt,v2pt
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,s13,s23
      real, intent(in), dimension(ib:ie,jb:je) :: ust,zntmp
      real, intent(inout), dimension(ib:ie,jb:je) :: stau
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr,rf
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, intent(in), dimension(ib:ie+1) :: uf
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie+1) :: arf1,arf2
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,s1
      real, intent(inout), dimension(kb:ke) :: u1b,v1b
      double precision, intent(in) :: avgsfcu,avgsfcv,avgsfcs
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      integer, intent(inout), dimension(rmp) :: reqs_s

      integer :: i,j,k,n
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23,wbar,rratio,fstar,tem,temx,temy
      real :: d2pm1,d2pm2,d2pm3,d2pm4,d2pm5,ueo,veo,dsdz,dsdzu,dsdzv,dudz,dvdz,up,vp,s13b,s23b
      double precision :: temd,teme,savgt,t13avg,t23avg,d13avg,d23avg,s13avg,s23avg,z13avg,z23avg,  &
                          ustbar,zntbar,kavgt,rustbar,ubar,vbar,ufrc,vfrc,ravg,rfavg,rflast,shravg,spavg

      real, dimension(kb:ke) :: wspa


    !-------------------------------------------

      u1b = 0.0
      v1b = 0.0
      savg = 0.0

    !-------------------------------------------
    !  Get domain average profiles:

      uavg = 0.0
      vavg = 0.0

      do k=1,nk
        !----
        cavg(k,1) = 0.0
        cavg(k,2) = 0.0
        !----
        do j=1,nj
        do i=1,ni
          cavg(k,1) = cavg(k,1) + ugr(i,j,k)
          cavg(k,2) = cavg(k,2) + vgr(i,j,k)
        enddo
        enddo
        !----
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,cavg(1,1),nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,cavg(1,2),nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      temd = 1.0/dble(nx*ny)

      do k=1,nk
        uavg(k)  = cavg(k,1)*temd
        vavg(k)  = cavg(k,2)*temd
        savg(k)  = max( sqrt( uavg(k)**2 + vavg(k)**2 ) , 1.0e-10 )
        u1b(k) = uavg(k)
        v1b(k) = vavg(k)
      enddo

      teme = 1.0/dble(nx*ny)

    !-------------------------------------------
    !  Get new value of gamma:

    doinggamk:  &
    if( t2pfac.ge.0.001 )then

    do k=2,ntwk

      spavg = 0.0
      s13b = 0.5*(uavg(k)-uavg(k-1))*rdz*mf(1,1,k)
      s23b = 0.5*(vavg(k)-vavg(k-1))*rdz*mf(1,1,k)

      do j=1,nj
      do i=1,ni

        tmp11=( c1(i,j,k)*(t11(i,j,k-1)**2) + c2(i,j,k)*(t11(i,j,k)**2) )
        tmp22=( c1(i,j,k)*(t22(i,j,k-1)**2) + c2(i,j,k)*(t22(i,j,k)**2) )
        tmp33=( c1(i,j,k)*(t33(i,j,k-1)**2) + c2(i,j,k)*(t33(i,j,k)**2) )

        tmp12=0.25*( c1(i,j,k)*( ( (t12(i  ,j  ,k-1)**2)     &
                                 + (t12(i+1,j+1,k-1)**2) )   &
                               + ( (t12(i  ,j+1,k-1)**2)     &
                                 + (t12(i+1,j  ,k-1)**2) ) ) &
                    +c2(i,j,k)*( ( (t12(i  ,j  ,k  )**2)     &
                                 + (t12(i+1,j+1,k  )**2) )   &
                               + ( (t12(i  ,j+1,k  )**2)     &
                                 + (t12(i+1,j  ,k  )**2) ) ) )

        tmp13=0.5*( (t13(i,j,k)-s13b)**2 + (t13(i+1,j,k)-s13b)**2 )

        tmp23=0.5*( (t23(i,j,k)-s23b)**2 + (t23(i,j+1,k)-s23b)**2 )

        spavg = spavg + ( 4.0*( tmp12 + ( tmp13 + tmp23 ) ) + 2.0*( ( tmp11 + tmp22 ) + tmp33 ) )/(rf(i,j,k)**2)

      enddo
      enddo

#ifdef MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE,spavg,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

      spavg = sqrt( spavg*teme )

      shravg = sqrt( 4.0*( s13b**2 + s23b**2 ) )

      gamk(k) = spavg/( smeps + (spavg+gamwall(k)*shravg) )

      s2p(k) = spavg
      s2b(k) = shravg

    enddo

    gamk(1) = gamk(2)

    endif  doinggamk

      if(timestats.ge.1) time_turb=time_turb+mytime()

      end subroutine getgamk

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx3
      real function flx3(s1,s2,s3)
      implicit none

      real, intent(in) :: s1,s2,s3

      ! 3rd-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx3 = (  (-1.0/6.0)*s1  &
               +( 5.0/6.0)*s2  &
               +( 2.0/6.0)*s3  )

      end function flx3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx4
      real function flx4(s1,s2,s3,s4)
      implicit none

      real, intent(in) :: s1,s2,s3,s4

      ! 4th-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx4 = (  (7.0/12.0)*(s3+s2)  &
               -(1.0/12.0)*(s4+s1)  )

      end function flx4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx5
      real function flx5(s1,s2,s3,s4,s5)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5

      ! 5th-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx5 = (  (  2.0/60.0)*s1  &
               +(-13.0/60.0)*s2  &
               +( 47.0/60.0)*s3  &
               +( 27.0/60.0)*s4  &
               +( -3.0/60.0)*s5  )

      end function flx5

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !dir$ attributes forceinline :: flx6
      real function flx6(s1,s2,s3,s4,s5,s6)
      implicit none

      real, intent(in) :: s1,s2,s3,s4,s5,s6

      ! 6th-order flux (eg, Wicker and Skamarock, 2002, MWR)

      flx6 = (  (37.0/60.0)*(s4+s3)  &
               +(-8.0/60.0)*(s5+s2)  &
               +( 1.0/60.0)*(s6+s1)  )

      end function flx6


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getepsd(xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,ua ,va ,wa ,  &
                         dum1,dum2,epsd1,epsd2,rho0,rf0,rrf0,rru,rrv,rrw,dt)
      use input
      use constants
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: epsd1,epsd2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,rf0,rrf0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in) :: dt

      integer :: i,j,k
      real :: ubar,vbar,wbar
      real :: tem
      real :: coef
      logical :: doit

      real, parameter :: onedtwelve = 1.0/12.0
      real, parameter :: onedsixty  = 1.0/60.0

    IF(.not.terrain_flag)THEN
      ! without terrain:

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k,tem)
      DO k=1,nk
        tem = rho0(1,1,k)
        do j=0,nj+2
        do i=0,ni+2
          rru(i,j,k)=tem*ua(i,j,k)
          rrv(i,j,k)=tem*va(i,j,k)
        enddo
        enddo
        IF(k.eq.1)THEN
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,   1) = 0.0
            rrw(i,j,nk+1) = 0.0
          enddo
          enddo
        ELSE
          tem = rf0(1,1,k)
          do j=0,nj+1
          do i=0,ni+1
            rrw(i,j,k)=tem*wa(i,j,k)
          enddo
          enddo
        ENDIF
      ENDDO

      doit = .true.

    ELSE
      ! with terrain:

      epsd1 = 0.0
      epsd2 = 0.0

      doit = .false.

    ENDIF

  docalc:  &
  IF( doit )THEN

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  acheck:  IF( axisymm.eq.0 )THEN

!----------
!  u-x:
!  v-y:

  IF( hadvordrv.eq.5 )THEN
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,ubar,vbar)
    do k=1,nk
    do j=1,nj
    do i=1,ni
      ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
      dum1(i,j,k) = (ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)*(  &
                      ubar*sign(1.0,ubar)               &
             *( 10.0*(ua(i+1,j,k)-ua(i  ,j,k))          &
                -5.0*(ua(i+2,j,k)-ua(i-1,j,k))          &
                    +(ua(i+3,j,k)-ua(i-2,j,k)) )*onedsixty )
      vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
      dum2(i,j,k) = (va(i,j+1,k)-va(i,j,k))*rdy*vh(j)*(  &
                      vbar*sign(1.0,vbar)               &
             *( 10.0*(va(i,j+1,k)-va(i,j  ,k))          &
                -5.0*(va(i,j+2,k)-va(i,j-1,k))          &
                    +(va(i,j+3,k)-va(i,j-2,k)) )*onedsixty )
    enddo
    enddo
    enddo
  ELSEIF( idiff.ge.1 .and. difforder.eq.6 )THEN
    coef = kdiff6/64.0/dt
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=1,nk
    do j=1,nj
    do i=1,ni
      dum1(i,j,k) = coef*( 10.0*(ua(i+1,j,k)-ua(i  ,j,k))     &
                           -5.0*(ua(i+2,j,k)-ua(i-1,j,k))     &
                               +(ua(i+3,j,k)-ua(i-2,j,k)) )   &
                        *(ua(i+1,j,k)-ua(i,j,k))
      dum2(i,j,k) = coef*( 10.0*(va(i,j+1,k)-va(i,j  ,k))     &
                           -5.0*(va(i,j+2,k)-va(i,j-1,k))     &
                               +(va(i,j+3,k)-va(i,j-2,k)) )   &
                        *(va(i,j+1,k)-va(i,j,k))
    enddo
    enddo
    enddo
  ELSE
    dum1 = 0.0
    dum2 = 0.0
  ENDIF

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          = ( c1(1,1,k)*(dum1(i,j,k-1)+dum2(i,j,k-1)) &
                      +c2(1,1,k)*(dum1(i,j,k  )+dum2(i,j,k  )) )
      epsd1(i,j,k) = max( 0.0 , tem )
      epsd2(i,j,k) = min( 0.0 , tem )
    enddo
    enddo
    enddo

!----------
!  u-y:
!  v-x:

  IF( hadvordrv.eq.5 )THEN
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,ubar,vbar)
    do k=1,nk
    do j=1,nj+1
    do i=1,ni+1
      vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
      ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
      dum1(i,j,k) = (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)*(  &
                      vbar*sign(1.0,vbar)               &
             *( 10.0*(ua(i,j  ,k)-ua(i,j-1,k))          &
                -5.0*(ua(i,j+1,k)-ua(i,j-2,k))          &
                    +(ua(i,j+2,k)-ua(i,j-3,k)) )*onedsixty ) &
                  + (va(i,j,k)-va(i-1,j,k))*rdx*uf(i)*(  &
                      ubar*sign(1.0,ubar)               &
             *( 10.0*(va(i  ,j,k)-va(i-1,j,k))          &
                -5.0*(va(i+1,j,k)-va(i-2,j,k))          &
                    +(va(i+2,j,k)-va(i-3,j,k)) )*onedsixty )
    enddo
    enddo
    enddo
  ELSEIF( idiff.ge.1 .and. difforder.eq.6 )THEN
    coef = kdiff6/64.0/dt
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=1,nk
    do j=1,nj+1
    do i=1,ni+1
      dum1(i,j,k) = coef*(ua(i,j,k)-ua(i,j-1,k))*(            &
                         ( 10.0*(ua(i,j  ,k)-ua(i,j-1,k))     &
                           -5.0*(ua(i,j+1,k)-ua(i,j-2,k))     &
                               +(ua(i,j+2,k)-ua(i,j-3,k)) ) ) &
                  + coef*(va(i,j,k)-va(i-1,j,k))*(            &
                         ( 10.0*(va(i  ,j,k)-va(i-1,j,k))     &
                           -5.0*(va(i+1,j,k)-va(i-2,j,k))     &
                               +(va(i+2,j,k)-va(i-3,j,k)) ) )
    enddo
    enddo
    enddo
  ELSE
    dum1 = 0.0
    dum2 = 0.0
  ENDIF

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          =                0.25*(                                                   &
         ( c2(1,1,k)*((dum1(i,j,k  )+dum1(i+1,j+1,k  ))+(dum1(i,j+1,k  )+dum1(i+1,j,k  )))   &
          +c1(1,1,k)*((dum1(i,j,k-1)+dum1(i+1,j+1,k-1))+(dum1(i,j+1,k-1)+dum1(i+1,j,k-1))) ) )
      epsd1(i,j,k) = epsd1(i,j,k) + max( 0.0 , tem )
      epsd2(i,j,k) = epsd2(i,j,k) + min( 0.0 , tem )
    enddo
    enddo
    enddo

!----------
!  w-x:
!  w-y:

  IF( hadvordrv.eq.5 )THEN
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,ubar,vbar)
    do k=2,nk
    do j=1,nj+1
    do i=1,ni+1
      ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
      dum1(i,j,k) = (wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)*(  &
                      ubar*sign(1.0,ubar)               &
             *( 10.0*(wa(i  ,j,k)-wa(i-1,j,k))          &
                -5.0*(wa(i+1,j,k)-wa(i-2,j,k))          &
                    +(wa(i+2,j,k)-wa(i-3,j,k)) )*onedsixty )
      vbar = c2(1,1,k)*rrv(i,j,k)+c1(1,1,k)*rrv(i,j,k-1)
      dum2(i,j,k) = (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)*(  &
                      vbar*sign(1.0,vbar)               &
             *( 10.0*(wa(i,j  ,k)-wa(i,j-1,k))          &
                -5.0*(wa(i,j+1,k)-wa(i,j-2,k))          &
                    +(wa(i,j+2,k)-wa(i,j-3,k)) )*onedsixty )
    enddo
    enddo
    enddo
  ELSEIF( idiff.ge.1 .and. difforder.eq.6 )THEN
    coef = kdiff6/64.0/dt
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k)
    do k=2,nk
    do j=1,nj+1
    do i=1,ni+1
      dum1(i,j,k) = coef*(wa(i,j,k)-wa(i-1,j,k))*(    &
                   ( 10.0*(wa(i  ,j,k)-wa(i-1,j,k))   &
                     -5.0*(wa(i+1,j,k)-wa(i-2,j,k))   &
                         +(wa(i+2,j,k)-wa(i-3,j,k)) ) )
      dum2(i,j,k) = coef*(wa(i,j,k)-wa(i,j-1,k))*(    &
                   ( 10.0*(wa(i,j  ,k)-wa(i,j-1,k))   &
                     -5.0*(wa(i,j+1,k)-wa(i,j-2,k))   &
                         +(wa(i,j+2,k)-wa(i,j-3,k)) ) )
    enddo
    enddo
    enddo
  ELSE
    dum1 = 0.0
    dum2 = 0.0
  ENDIF

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem = 0.5*( (dum1(i,j,k)+dum1(i+1,j,k)) &
                 +(dum2(i,j,k)+dum2(i,j+1,k)) )
      epsd1(i,j,k) = epsd1(i,j,k) + max( 0.0 , tem )
      epsd2(i,j,k) = epsd2(i,j,k) + min( 0.0 , tem )
    enddo
    enddo
    enddo

!----------

  ELSE

    stop 23232

  ENDIF    acheck

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  vertical terms:

  dum1 = 0.0
  dum2 = 0.0

!----------
! u-z:

  IF( vadvordrv.eq.5 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,wbar)
    do k=4,nk-2
    do j=1,nj
    do i=1,ni+1
      wbar = 0.5*(rrw(i-1,j,k)+rrw(i,j,k))
      if( wbar.ge.0.0 )then
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(ua(i,j,k-3),ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1))    &
             -flx6(ua(i,j,k-3),ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1),ua(i,j,k+2))  )
      else
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(ua(i,j,k+2),ua(i,j,k+1),ua(i,j,k  ),ua(i,j,k-1),ua(i,j,k-2))    &
             -flx6(ua(i,j,k-3),ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1),ua(i,j,k+2))  )
      endif
    enddo
    enddo
    enddo

    k = 3
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni+1
      wbar = 0.5*(rrw(i-1,j,k)+rrw(i,j,k))
      if( wbar.ge.0.0 )then
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ))   &
             -flx4(ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1))  )
      else
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(ua(i,j,k+2),ua(i,j,k+1),ua(i,j,k  ),ua(i,j,k-1),ua(i,j,k-2))    &
             -flx4(ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1))  )
      endif
    enddo
    enddo

    k = nk-1
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni+1
      wbar = 0.5*(rrw(i-1,j,k)+rrw(i,j,k))
      if( wbar.ge.0.0 )then
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(ua(i,j,k-3),ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1))    &
             -flx4(ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1))  )
      else
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(ua(i,j,k+1),ua(i,j,k),ua(i,j,k-1))    &
             -flx4(ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ),ua(i,j,k+1))  )
      endif
    enddo
    enddo

    k = 2
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni+1
      wbar = 0.5*(rrw(i-1,j,k)+rrw(i,j,k))
      if( wbar.lt.0.0 )then
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(ua(i,j,k+1),ua(i,j,k),ua(i,j,k-1))   &
             -( c1(1,1,k)*ua(i,j,k-1)+c2(1,1,k)*ua(i,j,k  ) ) )
      else
        dum1(i,j,k) = 0.0
      endif
    enddo
    enddo

    k = nk
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni+1
      wbar = 0.5*(rrw(i-1,j,k)+rrw(i,j,k))
      if( wbar.gt.0.0 )then
        dum1(i,j,k) = dum1(i,j,k) + (ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(ua(i,j,k-2),ua(i,j,k-1),ua(i,j,k  ))  &
             -( c1(1,1,k)*ua(i,j,k-1)+c2(1,1,k)*ua(i,j,k  ) ) )
      else
        dum1(i,j,k) = 0.0
      endif
    enddo
    enddo

  ELSE

    dum1 = 0.0

  ENDIF

!----------
! v-z:

  IF( vadvordrv.eq.5 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,wbar)
    do k=4,nk-2
    do j=1,nj+1
    do i=1,ni
      wbar = 0.5*(rrw(i,j-1,k)+rrw(i,j,k))
      if( wbar.ge.0.0 )then
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(va(i,j,k-3),va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1))    &
             -flx6(va(i,j,k-3),va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1),va(i,j,k+2))  )
      else
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(va(i,j,k+2),va(i,j,k+1),va(i,j,k  ),va(i,j,k-1),va(i,j,k-2))    &
             -flx6(va(i,j,k-3),va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1),va(i,j,k+2))  )
      endif
    enddo
    enddo
    enddo

    k = 3
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj+1
    do i=1,ni
      wbar = 0.5*(rrw(i,j-1,k)+rrw(i,j,k))
      if( wbar.ge.0.0 )then
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(va(i,j,k-2),va(i,j,k-1),va(i,j,k  ))   &
             -flx4(va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1))  )
      else
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(va(i,j,k+2),va(i,j,k+1),va(i,j,k  ),va(i,j,k-1),va(i,j,k-2))    &
             -flx4(va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1))  )
      endif
    enddo
    enddo

    k = nk-1
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj+1
    do i=1,ni
      wbar = 0.5*(rrw(i,j-1,k)+rrw(i,j,k))
      if( wbar.ge.0.0 )then
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx5(va(i,j,k-3),va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1))    &
             -flx4(va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1))  )
      else
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(va(i,j,k+1),va(i,j,k),va(i,j,k-1))    &
             -flx4(va(i,j,k-2),va(i,j,k-1),va(i,j,k  ),va(i,j,k+1))  )
      endif
    enddo
    enddo

    k = 2
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj+1
    do i=1,ni
      wbar = 0.5*(rrw(i,j-1,k)+rrw(i,j,k))
      if( wbar.lt.0.0 )then
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(va(i,j,k+1),va(i,j,k),va(i,j,k-1))   &
             -( c1(1,1,k)*va(i,j,k-1)+c2(1,1,k)*va(i,j,k  ) ) )
      else
        dum2(i,j,k) = 0.0
      endif
    enddo
    enddo

    k = nk
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj+1
    do i=1,ni
      wbar = 0.5*(rrw(i,j-1,k)+rrw(i,j,k))
      if( wbar.gt.0.0 )then
        dum2(i,j,k) = dum2(i,j,k) + (va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)*wbar*(   &
              flx3(va(i,j,k-2),va(i,j,k-1),va(i,j,k  ))  &
             -( c1(1,1,k)*va(i,j,k-1)+c2(1,1,k)*va(i,j,k  ) ) )
      else
        dum2(i,j,k) = 0.0
      endif
    enddo
    enddo

  ELSE

    dum2 = 0.0

  ENDIF

!----------

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          = -0.5*( (dum1(i,j,k)+dum1(i+1,j,k)) &
                           +(dum2(i,j,k)+dum2(i,j+1,k)) )
      epsd1(i,j,k) = epsd1(i,j,k) + max( 0.0 , tem )
      epsd2(i,j,k) = epsd2(i,j,k) + min( 0.0 , tem )
    enddo
    enddo
    enddo

!----------
! w-z:

  dum1 = 0.0

  acheck2:  IF( axisymm.eq.0 )THEN

  IF( vadvordrv.eq.5 )THEN

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,wbar)
    do k=3,nk-2
    do j=1,nj
    do i=1,ni
      wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
      if( wbar.ge.0.0 )then
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
            flx5(wa(i,j,k-2),wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2))  &
           -flx6(wa(i,j,k-2),wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2),wa(i,j,k+3)) )
      else
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
            flx5(wa(i,j,k+3),wa(i,j,k+2),wa(i,j,k+1),wa(i,j,k  ),wa(i,j,k-1))  &
           -flx6(wa(i,j,k-2),wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2),wa(i,j,k+3)) )
      endif
    enddo
    enddo
    enddo

    k = 2
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni
      wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
      if( wbar.ge.0.0 )then
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
            flx3(wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1))  &
           -flx4(wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2)) )
      else
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
            flx5(wa(i,j,k+3),wa(i,j,k+2),wa(i,j,k+1),wa(i,j,k  ),wa(i,j,k-1))  &
           -flx4(wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2)) )
      endif
    enddo
    enddo

    k = nk-1
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni
      wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
      if( wbar.ge.0.0 )then
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
            flx5(wa(i,j,k-2),wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2))  &
           -flx4(wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2)) )
      else
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
            flx3(wa(i,j,k+2),wa(i,j,k+1),wa(i,j,k  ))  &
           -flx4(wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1),wa(i,j,k+2)) )
      endif
    enddo
    enddo

    k = 1
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni
      wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
      if( wbar.lt.0.0 )then
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
                         flx3(wa(i,j,k+2),wa(i,j,k+1),wa(i,j,k  )) &
                        -0.5*(wa(i,j,k)+wa(i,j,k+1)) )
      else
        dum1(i,j,k) = 0.0
      endif
    enddo
    enddo

    k = nk
    !$omp parallel do default(shared)   &
    !$omp private(i,j,wbar)
    do j=1,nj
    do i=1,ni
      wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
      if( wbar.gt.0.0 )then
        dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*wbar*(  &
                         flx3(wa(i,j,k-1),wa(i,j,k  ),wa(i,j,k+1)) &
                        -0.5*(wa(i,j,k)+wa(i,j,k+1)) )
      else
        dum1(i,j,k) = 0.0
      endif
    enddo
    enddo

  ELSE

    dum1 = 0.0

  ENDIF

  ELSE

    dum1 = 0.0

  ENDIF   acheck2

!----------

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          =               -(c1(1,1,k)*dum1(i,j,k-1)+c2(1,1,k)*dum1(i,j,k))
      epsd1(i,j,k) = epsd1(i,j,k) + max( 0.0 , tem )
      epsd2(i,j,k) = epsd2(i,j,k) + min( 0.0 , tem )
      !-----
      epsd1(i,j,k) = epsd1(i,j,k)*rrf0(1,1,k)
      epsd2(i,j,k) = epsd2(i,j,k)*rrf0(1,1,k)
    enddo
    enddo
    enddo

!--------------------------------------------------------------
!  finished

  ENDIF  docalc

      if(timestats.ge.1) time_diag=time_diag+mytime()

      end subroutine getepsd


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getepst(xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,ua ,va ,wa ,  &
                         t11,t12,t13,t22,t23,t33,rf,                         &
                         dum1,dum2,epst,epstp,epstn)
      use input
      use constants
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in),    dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33,rf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: epst
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke), optional :: epstp,epstn

      integer :: i,j,k
      real :: rrf,tem

!--------------------------------------------------------------
!  horizontal terms, Cartesian grid:

!----------
!  u-x:
!  v-y:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=1,nk
    do j=1,nj
    do i=1,ni
      dum1(i,j,k) = (ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)*t11(i,j,k)
      dum2(i,j,k) = (va(i,j+1,k)-va(i,j,k))*rdy*vh(j)*t22(i,j,k)
    enddo
    enddo
    enddo

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          = ( c1(1,1,k)*(dum1(i,j,k-1)+dum2(i,j,k-1)) &
                      +c2(1,1,k)*(dum1(i,j,k  )+dum2(i,j,k  )) )
      epst(i,j,k) = tem
    enddo
    enddo
    enddo

    if( present(epstp) .and. present(epstn) )then
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          = ( c1(1,1,k)*(dum1(i,j,k-1)+dum2(i,j,k-1)) &
                      +c2(1,1,k)*(dum1(i,j,k  )+dum2(i,j,k  )) )
      epstp(i,j,k) = max( 0.0 , tem )
      epstn(i,j,k) = min( 0.0 , tem )
    enddo
    enddo
    enddo
    endif

!----------
!  u-y:
!  v-x:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=1,nk
    do j=1,nj+1
    do i=1,ni+1
      dum1(i,j,k) = t12(i,j,k)*( (ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j) &
                                +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i) )
    enddo
    enddo
    enddo

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          =                0.25*(                                                   &
         ( c2(1,1,k)*((dum1(i,j,k  )+dum1(i+1,j+1,k  ))+(dum1(i,j+1,k  )+dum1(i+1,j,k  )))   &
          +c1(1,1,k)*((dum1(i,j,k-1)+dum1(i+1,j+1,k-1))+(dum1(i,j+1,k-1)+dum1(i+1,j,k-1))) ) )
      epst(i,j,k) = epst(i,j,k) + tem
    enddo
    enddo
    enddo

    if( present(epstp) .and. present(epstn) )then
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          =                0.25*(                                                   &
         ( c2(1,1,k)*((dum1(i,j,k  )+dum1(i+1,j+1,k  ))+(dum1(i,j+1,k  )+dum1(i+1,j,k  )))   &
          +c1(1,1,k)*((dum1(i,j,k-1)+dum1(i+1,j+1,k-1))+(dum1(i,j+1,k-1)+dum1(i+1,j,k-1))) ) )
      epstp(i,j,k) = epstp(i,j,k) + max( 0.0 , tem )
      epstn(i,j,k) = epstn(i,j,k) + min( 0.0 , tem )
    enddo
    enddo
    enddo
    endif

!----------
!  w-x:
!  w-y:
!  u-z:
!  v-z:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj+1
    do i=1,ni+1
      dum1(i,j,k) = ( (wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)     &
                     +(ua(i,j,k)-ua(i,j,k-1))*rdz*mf(1,1,k) )*t13(i,j,k)
      dum2(i,j,k) = ( (wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)     &
                     +(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k) )*t23(i,j,k)
    enddo
    enddo
    enddo

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          = 0.5*( (dum1(i,j,k)+dum1(i+1,j,k)) &
                          +(dum2(i,j,k)+dum2(i,j+1,k)) )
      epst(i,j,k) = epst(i,j,k) + tem
    enddo
    enddo
    enddo

    if( present(epstp) .and. present(epstn) )then
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          = 0.5*( (dum1(i,j,k)+dum1(i+1,j,k)) &
                          +(dum2(i,j,k)+dum2(i,j+1,k)) )
      epstp(i,j,k) = epstp(i,j,k) + max( 0.0 , tem )
      epstn(i,j,k) = epstn(i,j,k) + min( 0.0 , tem )
    enddo
    enddo
    enddo
    endif

!----------
! w-z:

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem)
    do k=1,nk
    do j=1,nj
    do i=1,ni
      dum1(i,j,k) = (wa(i,j,k+1)-wa(i,j,k))*rdz*mh(1,1,k)*t33(i,j,k)
    enddo
    enddo
    enddo

    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem,rrf)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          =                (c1(1,1,k)*dum1(i,j,k-1)+c2(1,1,k)*dum1(i,j,k))
      epst(i,j,k) = epst(i,j,k) + tem
      !-------------------------
      rrf = 1.0/rf(i,j,k)
      epst(i,j,k) = epst(i,j,k)*rrf
    enddo
    enddo
    enddo

    if( present(epstp) .and. present(epstn) )then
    !$omp parallel do default(shared)   &
    !$omp private(i,j,k,tem,rrf)
    do k=2,nk
    do j=1,nj
    do i=1,ni
      tem          =                (c1(1,1,k)*dum1(i,j,k-1)+c2(1,1,k)*dum1(i,j,k))
      epstp(i,j,k) = epstp(i,j,k) + max( 0.0 , tem )
      epstn(i,j,k) = epstn(i,j,k) + min( 0.0 , tem )
      !-------------------------
      rrf = 1.0/rf(i,j,k)
      epstp(i,j,k) = epstp(i,j,k)*rrf
      epstn(i,j,k) = epstn(i,j,k)*rrf
    enddo
    enddo
    enddo
    endif

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()

      end subroutine getepst


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE turb_module
