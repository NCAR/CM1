  MODULE mmm_physics_wrapper_module

  use ccpp_kind_types,only: kind_phys
  implicit none

  private
  public :: ysu_wrapper,sf_sfclayrev_wrapper,cu_ntiedtke_wrapper,mp_wsm6_wrapper,mp_wsm6_effectRad_wrapper

  CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine ysu_wrapper(ux,vx,tx,qvx,qcx,qix,nmix,qmix,p2d,p2di,pi2d,     &
                         f_qc,f_qi,                                        &
                         utnp,vtnp,ttnp,qvtnp,qctnp,qitnp,qmixtnp,         &
                         cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,         &
                         dz8w2d,psfcpa,                                    &
                         znt,ust,hpbl,dusfc,dvsfc,dtsfc,dqsfc,psim,psih,   &
                         xland,hfx,qfx,wspd,br,                            &
                         dt,kpbl1d,                                        &
                         exch_hx,exch_mx,                                  &
                         wstar,delta,                                      &
                         u10,v10,                                          &
                         uox,vox,                                          &
                         rthraten,                                         &
                         ysu_topdown_pblmix,                               &
                         ctopo,ctopo2,                                     &
                         a_u,a_v,a_t,a_q,a_e,                              &
                         b_u,b_v,b_t,b_q,b_e,                              &
                         sfk,vlk,dlu,dlg,frcurb,                           &
                         flag_bep,                                         &
                         its,ite,kts,kte,                                  &
                         errmsg,errflg                                     &
                        )
  use bl_ysu , only : bl_ysu_run
   implicit none
!
   integer,  intent(in   )   ::     its,ite,kts,kte

   logical,  intent(in)      ::     ysu_topdown_pblmix
!
   integer,  intent(in)      ::     nmix
!
   real(kind=kind_phys),     intent(in   )   ::     dt,cp,g,rovcp,rovg,rd,xlv,rv
!
   real(kind=kind_phys),     intent(in )     ::     ep1,ep2,karman
!
   logical,  intent(in )     ::     f_qc, f_qi
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte )                    , &
             intent(in)      ::                                        dz8w2d, &
                                                                         pi2d
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte )                    , &
             intent(in   )   ::                                            tx, &
                                                                          qvx, &
                                                                          qcx, &
                                                                          qix
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte,nmix )               , &
             intent(in   )   ::                                          qmix
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte )                    , &
             intent(out  )   ::                                          utnp, &
                                                                         vtnp, &
                                                                         ttnp, &
                                                                        qvtnp, &
                                                                        qctnp, &
                                                                        qitnp
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte,nmix )               , &
             intent(out  )   ::                                       qmixtnp
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte+1 )                  , &
             intent(in   )   ::                                          p2di
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte )                    , &
             intent(in   )   ::                                           p2d
!
   real(kind=kind_phys),     dimension( its:ite )                               , &
             intent(out  )   ::                                          hpbl
!
   real(kind=kind_phys),     dimension( its:ite )                               , &
             intent(out  ), optional ::                                 dusfc, &
                                                                        dvsfc, &
                                                                        dtsfc, &
                                                                        dqsfc
!
   real(kind=kind_phys),     dimension( its:ite )                               , &
             intent(in   )   ::                                           ust, &
                                                                          znt
   real(kind=kind_phys),     dimension( its:ite )                               , &
             intent(in   )   ::                                         xland, &
                                                                          hfx, &
                                                                          qfx
!
   real(kind=kind_phys),     dimension( its:ite ), intent(in   )    ::      wspd
   real(kind=kind_phys),     dimension( its:ite ), intent(in   )    ::        br
!
   real(kind=kind_phys),     dimension( its:ite ), intent(in   )    ::      psim, &
                                                                         psih
!
   real(kind=kind_phys),     dimension( its:ite ), intent(in   )    ::    psfcpa
   integer,  dimension( its:ite ), intent(out  )   ::                     kpbl1d
!
   real(kind=kind_phys),     dimension( its:ite,kts:kte )                    , &
             intent(in   )   ::                                            ux, &
                                                                           vx, &
                                                                      rthraten
   real(kind=kind_phys),     dimension( its:ite )                               , &
             optional                                                        , &
             intent(in   )   ::                                         ctopo, &
                                                                       ctopo2
!
   logical,  intent(in   )   ::                                      flag_bep
   real(kind=kind_phys),     dimension( its:ite,kts:kte )                             , &
             optional                                                        , &
             intent(in   )   ::                                           a_u, &
                                                                          a_v, &
                                                                          a_t, &
                                                                          a_q, &
                                                                          a_e, &
                                                                          b_u, &
                                                                          b_v, &
                                                                          b_t, &
                                                                          b_q, &
                                                                          b_e, &
                                                                          sfk, &
                                                                          vlk, &
                                                                          dlu, &
                                                                          dlg
   real(kind=kind_phys),     dimension( its:ite )                               , &
             optional                                                        , &
             intent(in   )   ::                                        frcurb
!
   character(len=*), intent(out)   ::                                  errmsg
   integer,          intent(out)   ::                                  errflg

!
   real(kind=kind_phys),    dimension( its:ite, kts:kte+1 )                    , &
            intent(out  )   ::                                        exch_hx, &
                                                                      exch_mx
!
   real(kind=kind_phys),    dimension( its:ite )                             , &
            intent(inout)    ::                                           u10, &
                                                                          v10
   real(kind=kind_phys),    dimension( its:ite ), optional                   , &
            intent(in  )    ::                                            uox, &
                                                                          vox
   real(kind=kind_phys), dimension( its:ite ), intent(out) ::           wstar, &
                                                                        delta
!
   integer :: kme

   kme = kte+1

!    if(myid.eq.0) print *,'  ysu_wrapper '

     call     bl_ysu_run(ux=ux, &
                         vx=vx, &
                         tx=tx, &
                         qvx=qvx, &
                         qcx=qcx, &
                         qix=qix, &
                         nmix=nmix, &
                         qmix=qmix, &
                         p2d=p2d, &
                         p2di=p2di, &
                         pi2d=pi2d,     &
                         f_qc=f_qc, &
                         f_qi=f_qi,                                        &
                         utnp=utnp, &
                         vtnp=vtnp, &
                         ttnp=ttnp, &
                         qvtnp=qvtnp, &
                         qctnp=qctnp, &
                         qitnp=qitnp, &
                         qmixtnp=qmixtnp,         &
                         cp=cp, &
                         g=g, &
                         rovcp=rovcp, &
                         rd=rd, &
                         rovg=rovg, &
                         ep1=ep1, &
                         ep2=ep2, &
                         karman=karman, &
                         xlv=xlv, &
                         rv=rv,         &
                         dz8w2d=dz8w2d, &
                         psfcpa=psfcpa,                                    &
                         znt=znt, &
                         ust=ust, &
                         hpbl=hpbl, &
!                         dusfc=dusfc, &
!                         dvsfc=dvsfc, &
!                         dtsfc=dtsfc, &
!                         dqsfc=dqsfc, &
                         psim=psim, &
                         psih=psih,   &
                         xland=xland, &
                         hfx=hfx, &
                         qfx=qfx, &
                         wspd=wspd, &
                         br=br,                            &
                         dt=dt, &
                         kpbl1d=kpbl1d,                                        &
                         exch_hx=exch_hx, &
                         exch_mx=exch_mx,                                  &
                         wstar=wstar, &
                         delta=delta,                                      &
                         u10=u10, &
                         v10=v10,                                          &
!                         uox,vox,                                          &
                         rthraten=rthraten,                                         &
                         ysu_topdown_pblmix=ysu_topdown_pblmix,                               &
!                         ctopo,ctopo2,                                     &
!                         a_u,a_v,a_t,a_q,a_e,                              &
!                         b_u,b_v,b_t,b_q,b_e,                              &
!                         sfk,vlk,dlu,dlg,frcurb,                           &
                         flag_bep=flag_bep,                                         &
                         its=its, &
                         ite=ite, &
                         kte=kte, &
                         kme=kme,                                  &
                         errmsg=errmsg, &
                         errflg=errflg                                     &
                        )

    end subroutine ysu_wrapper

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine sf_sfclayrev_wrapper(ux,vx,t1d,qv1d,p1d,dz8w1d,                &
                             cp,g,rovcp,r,xlv,psfcpa,chs,chs2,cqs2,    &
                             cpm,pblh,rmol,znt,ust,mavail,zol,mol,     &
                             regime,psim,psih,fm,fh,                   &
                             xland,hfx,qfx,tsk,                        &
                             u10,v10,th2,t2,q2,flhc,flqc,qgh,          &
                             qsfc,lh,gz1oz0,wspd,br,isfflx,dx,         &
                             svp1,svp2,svp3,svpt0,ep1,ep2,             &
                             karman,p1000mb,                           &
                             shalwater_z0,water_depth,                 &
                             isftcflx,iz0tlnd,scm_force_flux,          &
                             ustm,ck,cka,cd,cda,                       &
                             its,ite,errmsg,errflg                     &
                            )
      use sf_sfclayrev, only: sf_sfclayrev_run
  implicit none

!=================================================================================================================

!--- input arguments:
 logical,intent(in):: isfflx
 logical,intent(in):: shalwater_z0
 logical,intent(in),optional:: scm_force_flux

 integer,intent(in):: its,ite
 integer,intent(in),optional:: isftcflx, iz0tlnd

 real(kind=kind_phys),intent(in):: svp1,svp2,svp3,svpt0
 real(kind=kind_phys),intent(in):: ep1,ep2,karman
 real(kind=kind_phys),intent(in):: p1000mb
 real(kind=kind_phys),intent(in):: cp,g,rovcp,r,xlv

 real(kind=kind_phys),intent(in),dimension(its:ite):: &
    mavail,     &
    pblh,       &
    psfcpa,     &
    tsk,        &
    xland,      &
    water_depth

 real(kind=kind_phys),intent(in),dimension(its:ite):: &
    dx,         &
    dz8w1d,     &    
    ux,         &
    vx,         &
    qv1d,       &
    p1d,        &
    t1d

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

 real(kind=kind_phys),intent(out),dimension(its:ite):: &
    lh,         &
    u10,        &
    v10,        &
    th2,        &
    t2,         &
    q2

 real(kind=kind_phys),intent(out),dimension(its:ite),optional:: &
    ck,         &
    cka,        &
    cd,         &
    cda

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite):: &
    regime,     &
    hfx,        &
    qfx,        &
    qsfc,       &
    mol,        &
    rmol,       &
    gz1oz0,     &
    wspd,       &
    br,         &
    psim,       &
    psih,       &
    fm,         &
    fh,         &
    znt,        &
    zol,        &
    ust,        &
    cpm,        &
    chs2,       &
    cqs2,       &
    chs,        &
    flhc,       &
    flqc,       &
    qgh

 real(kind=kind_phys),intent(inout),dimension(its:ite),optional:: &
    ustm

!  if(myid.eq.0) print *,'  sf_sfclayrev_wrapper '
   call     sf_sfclayrev_run(ux,vx,t1d,qv1d,p1d,dz8w1d,                &
                             cp,g,rovcp,r,xlv,psfcpa,chs,chs2,cqs2,    &
                             cpm,pblh,rmol,znt,ust,mavail,zol,mol,     &
                             regime,psim,psih,fm,fh,                   &
                             xland,hfx,qfx,tsk,                        &
                             u10,v10,th2,t2,q2,flhc,flqc,qgh,          &
                             qsfc,lh,gz1oz0,wspd,br,isfflx,dx,         &
                             svp1,svp2,svp3,svpt0,ep1,ep2,             &
                             karman,p1000mb,                           &
                             shalwater_z0,water_depth,                 &
                             isftcflx,iz0tlnd,scm_force_flux,          &
                             ustm,ck,cka,cd,cda,                       &
                             its,ite,errmsg,errflg                     &
                            )

 end subroutine sf_sfclayrev_wrapper

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cu_ntiedtke_wrapper(pu,pv,pt,pqv,pqc,pqi,pqvf,ptf,poz,pzz,pomg, &
     &         pap,paph,evap,hfx,zprecc,lndj,lq,km,km1,dt,dx,errmsg,errflg)
      use cu_ntiedtke, only : cu_ntiedtke_run
      implicit none

      integer,intent(in):: lq,km,km1
      integer,intent(in),dimension(lq):: lndj

      real(kind=kind_phys),intent(in):: dt
      real(kind=kind_phys),intent(in),dimension(lq):: dx
      real(kind=kind_phys),intent(in),dimension(lq):: evap,hfx
      real(kind=kind_phys),intent(in),dimension(lq,km):: pqvf,ptf
      real(kind=kind_phys),intent(in),dimension(lq,km):: poz,pomg,pap
      real(kind=kind_phys),intent(in),dimension(lq,km1):: pzz,paph

!--- inout arguments:
      real(kind=kind_phys),intent(inout),dimension(lq):: zprecc
      real(kind=kind_phys),intent(inout),dimension(lq,km):: pu,pv,pt,pqv,pqc,pqi

!--- output arguments:
      character(len=*),intent(out):: errmsg
      integer,intent(out):: errflg

!      if(myid.eq.0) print *,'  ntiedtke_wrapper '
      call       cu_ntiedtke_run(pu=pu,          &
                                 pv=pv,          &
                                 pt=pt,          &
                                 pqv=pqv,        &
                                 pqc=pqc,        &
                                 pqi=pqi,        &
                                 pqvf=pqvf,      &
                                 ptf=ptf,        &
                                 poz=poz,        &
                                 pzz=pzz,        &
                                 pomg=pomg,      &
                                 pap=pap,        &
                                 paph=paph,      &
                                 evap=evap,      &
                                 hfx=hfx,        &
                                 zprecc=zprecc,  &
                                 lndj=lndj,      &
                                 lq=lq,          &
                                 km=km,          &
                                 km1=km1,        &
                                 dt=dt,          &
                                 dx=dx,          &
                                 errmsg=errmsg,  &
                                 errflg=errflg)

      end subroutine cu_ntiedtke_wrapper

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine mp_wsm6_wrapper(t,q,qc,qi,qr,qs,qg,den,p,delz,delt,   &
                        g,cpd,cpv,rd,rv,t0c,ep1,ep2,qmin,xls, &
                        xlv0,xlf0,den0,denr,cliq,cice,psat,   &
                        rain,rainncv,sr,snow,snowncv,graupel, &
                        graupelncv,rainprod2d,evapprod2d,     &
                        its,ite,kts,kte,errmsg,errflg         &
                       )
 use ccpp_kind_types,only: kind_phys
    use mp_wsm6 , only : mp_wsm6_run
    implicit none

 integer,intent(in):: its,ite,kts,kte

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte)::              &
                                                             den, &
                                                               p, &
                                                            delz
 real(kind=kind_phys),intent(in)::                                &
                                                            delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             xls, &
                                                            xlv0, &
                                                            xlf0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr

!inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte)::           &
                                                               t
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte)::           &
                                                               q, &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs, &
                                                              qg
 real(kind=kind_phys),intent(inout),dimension(its:ite)::             &
                                                            rain, &
                                                         rainncv, &
                                                              sr

 real(kind=kind_phys),intent(inout),dimension(its:ite),optional::    &
                                                            snow, &
                                                         snowncv

 real(kind=kind_phys),intent(inout),dimension(its:ite),optional::    &
                                                         graupel, &
                                                      graupelncv

 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte),optional::  &
                                                      rainprod2d, &
                                                      evapprod2d

!output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!  if(myid.eq.0) print *,'  mp_wsm6_wrapper '
   call     mp_wsm6_run(t=t,                    &
                        q=q,                    &
                        qc=qc,                  &
                        qi=qi,                  &
                        qr=qr,                  &
                        qs=qs,                  &
                        qg=qg,                  &
                        den=den,                &
                        p=p,                    &
                        delz=delz,              &
                        delt=delt,              &
                        g=g,                    &
                        cpd=cpd,                &
                        cpv=cpv,                &
                        rd=rd,                  &
                        rv=rv,                  &
                        t0c=t0c,                &
                        ep1=ep1,                &
                        ep2=ep2,                &
                        qmin=qmin,              &
                        xls=xls,                &
                        xlv0=xlv0,              &
                        xlf0=xlf0,              &
                        den0=den0,              &
                        denr=denr,              &
                        cliq=cliq,              &
                        cice=cice,              &
                        psat=psat,              &
                        rain=rain,              &
                        rainncv=rainncv,        &
                        sr=sr,                  &
                        snow=snow,              &
!                        snowncv=snowncv,        &
                        graupel=graupel,        &
!                        graupelncv=graupelncv,  &
!                        rainprod2d=rainprod2d,  &
!                        evapprod2d=evapprod2d,  &
                        its=its,                &
                        ite=ite,                &
                        kts=kts,                &
                        kte=kte,                &
                        errmsg=errmsg,          &
                        errflg=errflg           &
                       )

 end subroutine mp_wsm6_wrapper

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine mp_wsm6_effectRad_wrapper(do_microp_re,t,qc,qi,qs,rho,qmin,t0c,re_qc_bg,re_qi_bg,re_qs_bg, &
                                  re_qc_max,re_qi_max,re_qs_max,re_qc,re_qi,re_qs,its,ite,kts,kte, &
                                  errmsg,errflg)
 use ccpp_kind_types,only: kind_phys
    use mp_wsm6_effectrad, only : mp_wsm6_effectRad_run
    implicit none

!..Sub arguments
 logical,intent(in):: do_microp_re
 integer,intent(in):: its,ite,kts,kte
 real(kind=kind_phys),intent(in):: qmin
 real(kind=kind_phys),intent(in):: t0c
 real(kind=kind_phys),intent(in):: re_qc_bg,re_qi_bg,re_qs_bg
 real(kind=kind_phys),intent(in):: re_qc_max,re_qi_max,re_qs_max
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(in)::  t
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(in)::  qc
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(in)::  qi
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(in)::  qs
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(in)::  rho
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(inout):: re_qc
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(inout):: re_qi
 real(kind=kind_phys),dimension(its:ite,kts:kte),intent(inout):: re_qs

!...Output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!  if(myid.eq.0) print *,'  mp_wsm6_effectRad_wrapper '
    call    mp_wsm6_effectRad_run(do_microp_re=do_microp_re,  &
                                  t=t,                        &
                                  qc=qc,                      &
                                  qi=qi,                      &
                                  qs=qs,                      &
                                  rho=rho,                    &
                                  qmin=qmin,                  &
                                  t0c=t0c,                    &
                                  re_qc_bg=re_qc_bg,          &
                                  re_qi_bg=re_qi_bg,          &
                                  re_qs_bg=re_qs_bg,          &
                                  re_qc_max=re_qc_max,        &
                                  re_qi_max=re_qi_max,        &
                                  re_qs_max=re_qs_max,        &
                                  re_qc=re_qc,                &
                                  re_qi=re_qi,                &
                                  re_qs=re_qs,                &
                                  its=its,                    &
                                  ite=ite,                    &
                                  kts=kts,                    &
                                  kte=kte,                    &
                                  errmsg=errmsg,              &
                                  errflg=errflg)

 end subroutine mp_wsm6_effectRad_wrapper

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE mmm_physics_wrapper_module
