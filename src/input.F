  MODULE input

  implicit none

  public

      integer, parameter :: maxq = 100   ! maximum possible number of
                                         ! q variables

      integer, parameter :: maxvars = 10000   ! maximum possible number of
                                              ! output variables

      integer, parameter :: maxstring = 120

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                                       cccccccccccccccccccccccccccccc
!cc   Do not change anything below here   cccccccccccccccccccccccccccccc
!cc                                       cccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

      integer ierr

!-----------------------------------

      logical terrain_flag,procfiles,dowr,                                    &
              patchsws,patchsww,patchses,patchsee,                            &
              patchnwn,patchnww,patchnen,patchnee,                            &
              p2tchsws,p2tchsww,p2tchses,p2tchsee,                            &
              p2tchnwn,p2tchnww,p2tchnen,p2tchnee,stopit,                     &
              restart_file_theta,restart_file_dbz,restart_file_th0,           &
              restart_file_prs0,restart_file_pi0,restart_file_rho0,           &
              restart_file_qv0,restart_file_u0,restart_file_v0,               &
              restart_file_zs,restart_file_zh,restart_file_zf,                &
              restart_file_diags,restart_use_theta,restart_reset_frqtim,      &
              dolsw,dodomaindiag,doazimavg,dohifrq,dohturb,dovturb,pdcomp,    &
              use_pbl,use_avg_sfc,dot2p,iusetke,do_adapt_move,                &
              bl_mynn_tkeadvect,idoles,idopbl,dotimeavg,interp_on_restart,    &
              do_recycle_w,do_recycle_s,do_recycle_e,do_recycle_n,            &
              do_lsnudge,do_lsnudge_u,do_lsnudge_v,do_lsnudge_th,do_lsnudge_qv, &
              do_ib

!-----------------------------------------------------------------------

      integer nx,ny,nz,nodex,nodey,ppnode,timeformat,timestats,outunits,      &
              ni,nj,nk,nkp1,ngxy,ngz,                                         &
              ib,ie,jb,je,kb,ke,                                              &
              ibm,iem,jbm,jem,kbm,kem,                                        &
              ibi,iei,jbi,jei,kbi,kei,iice,idm,idmplus,                       &
              ibc,iec,jbc,jec,kbc,kec,                                        &
              ibt,iet,jbt,jet,kbt,ket,                                        &
              ibp,iep,jbp,jep,kbp,kep,                                        &
              itb,ite,jtb,jte,ktb,kte,                                        &
              ipb,ipe,jpb,jpe,kpb,kpe,                                        &
              ibr,ier,jbr,jer,kbr,ker,nir,njr,nkr,                            &
              ibb,ieb,jbb,jeb,kbb,keb,                                        &
              ibp3,iep3,jbp3,jep3,kbp3,kep3,np3a,np3o,                        &
              ibmynn,iemynn,jbmynn,jemynn,kbmynn,kemynn,                      &
              ibmyj ,iemyj ,jbmyj ,jemyj ,kbmyj ,kemyj ,                      &
              ibta  ,ieta  ,jbta  ,jeta  ,kbta  ,keta  ,                      &
              ibcu  ,iecu  ,jbcu  ,jecu  ,kbcu  ,kecu  ,                      &
              ibtmp ,ietmp ,jbtmp ,jetmp ,kbtmp ,ketmp ,ntmp2d, ntmp1d,       &
              ibdt,iedt,jbdt,jedt,kbdt,kedt,ntdiag,                           &
              ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,                           &
              ibdv,iedv,jbdv,jedv,kbdv,kedv,nudiag,nvdiag,nwdiag,             &
              ibdk,iedk,jbdk,jedk,kbdk,kedk,nkdiag,                           &
              ibdp,iedp,jbdp,jedp,kbdp,kedp,npdiag,                           &
              ib2d,ie2d,jb2d,je2d,nout2d,                                     &
              ib3d,ie3d,jb3d,je3d,kb3d,ke3d,nout3d,                           &
              ibph,ieph,jbph,jeph,kbph,keph,                                  &
              ibl,iel,jbl,jel,                                                &
              ibb2,ibe2,jbb2,jbe2,kbb2,kbe2,                                  &
              d2i,d2is,d2iu,d2iv,                                             &
              d2j,d2js,d2ju,d2jv,                                             &
              d3i,d3is,d3iu,d3iv,                                             &
              d3j,d3js,d3ju,d3jv,                                             &
              d3n,d3t,mynode,nodeleader,nodes,                                &
              ibzvd,iezvd,jbzvd,jezvd,kbzvd,kezvd,nqzvd,                      &
              imirror,jmirror,                                                &
              imp,jmp,kmp,kmt,rmp,cmp,nrain,                                  &
              numq,nqv,nqc,nqr,nqi,nqs,nqg,nnci,nncc,nqi2,nqi3,               &
              prx,pry,prz,pru,prv,prw,prth,prt,prprs,prtime,                  &
              prpt1,prpt2,prqv,prq1,prq2,prnc1,prnc2,prkm,prkh,prtke,         &
              prdbz,prb,prvpg,przv,prrho,prqsl,prqsi,prznt,prust,przs,prsig,  &
              nql1,nql2,nqs1,nqs2,nnc1,nnc2,nvl1,nvl2,nzl1,nzl2,              &
              nbudget,budrain,cm1setup,testcase,                              &
              adapt_dt,irst,rstnum,iconly,                                    &
              hadvordrs,vadvordrs,hadvordrv,vadvordrv,pdscheme,apmasscon,     &
              advwenos,advwenov,weno_order,                                   &
              idiff,mdiff,difforder,imoist,sgsmodel,horizturb,                &
              tconfig,bcturbs,doimpl,                                         &
              irdamp,hrdamp,psolver,nsound,ptype,ihail,iautoc,cuparam,        &
              icor,betaplane,lspgrad,eqtset,idiss,efall,rterm,                &
              wbc,ebc,sbc,nbc,bbc,tbc,irbc,roflux,nudgeobc,                   &
              isnd,iwnd,itern,iinit,                                          &
              irandp,ibalance,iorigin,axisymm,imove,iptra,npt,pdtra,          &
              iprcl,nparcels,                                                 &
              stretch_x,stretch_y,stretch_z,                                  &
              bc_temp,ibw,ibe,ibs,ibn,npvals,                                 &
              outfile,myid,numprocs,myi,myj,nf,nu,nv,nw,                      &
              mywest,myeast,mysouth,mynorth,mysw,mynw,myne,myse,              &
              myi1,myi2,myj1,myj2,nimax,njmax,                                &
              cs1we,cs1sn,ct1we,ct1sn,cv1we,cu1sn,cw1we,cw1sn,cs2we,cs2sn,    &
              cs3we,cs3sn,ct3we,ct3sn,cv3we,cu3sn,cw3we,cw3sn,cs3weq,cs3snq,  &
              output_format,output_filetype,output_interp,                    &
              restart_format,restart_filetype,                                &
              output_rain,output_sws,output_svs,output_sps,output_srs,        &
              output_sgs,output_sus,output_shs,output_coldpool,output_zs,     &
              output_psfc,                                                    &
              output_basestate,output_sfcflx,output_sfcparams,output_sfcdiags,&
              output_zh,output_th,output_thpert,output_prs,output_prspert,    &
              output_pi,output_pipert,output_rho,output_rhopert,output_tke,   &
              output_km,output_kh,                                            &
              output_qv,output_qvpert,output_q,output_dbz,output_buoyancy,    &
              output_u,output_upert,output_uinterp,                           &
              output_v,output_vpert,output_vinterp,output_w,output_winterp,   &
              output_vort,output_pv,output_uh,output_pblten,                  &
              output_dissten,output_fallvel,output_nm,output_def,             &
              output_radten,output_cape,output_cin,output_lcl,output_lfc,     &
              output_pwat,output_lwp,                                         &
              output_thbudget,output_qvbudget,                                &
              output_ubudget,output_vbudget,output_wbudget,output_pdcomp,     &
              prcl_th,prcl_t,prcl_prs,prcl_ptra,prcl_q,prcl_nc,               &
              prcl_km,prcl_kh,prcl_tke,prcl_dbz,prcl_b,prcl_vpg,prcl_vort,    &
              prcl_rho,prcl_qsat,prcl_sfc,                                    &
              n_out,s_out,u_out,v_out,w_out,z_out,sout2d,sout3d,              &
              stat_w,stat_wlevs,stat_u,stat_v,stat_rmw,                       &
              stat_pipert,stat_prspert,stat_thpert,stat_q,                    &
              stat_tke,stat_km,stat_kh,stat_div,stat_rh,stat_rhi,stat_the,    &
              stat_cloud,stat_sfcprs,stat_wsp,stat_cfl,stat_vort,             &
              stat_tmass,stat_tmois,stat_qmass,stat_tenerg,stat_mo,stat_tmf,  &
              stat_pcn,stat_qsrc,stat_out,prcl_out,                           &
              radopt,year,month,day,hour,minute,second,jday,                  &
              isfcflx,sfcmodel,oceanmodel,ipbl,initsfc,lu0,season,            &
              cecd,pertflx,isftcflx,iz0tlnd,convinit,wnudge,maxk,             &
              set_znt,set_flx,set_ust,ramp_sgs,t2p_avg,                       &
              qd_dbz,qd_vtc,qd_vtr,qd_vts,qd_vtg,qd_vti,                      &
              td_hadv,td_vadv,td_hturb,td_vturb,td_mp,td_rdamp,               &
              td_rad,td_div,td_diss,td_pbl,td_lsw,td_efall,                   &
              td_cond,td_evac,td_evar,td_dep,td_subl,td_melt,td_frz,          &
              td_nudge,                                                       &
              qd_hadv,qd_vadv,qd_hturb,qd_vturb,qd_mp,qd_pbl,qd_lsw,          &
              qd_cond,qd_evac,qd_evar,qd_dep,qd_subl,qd_nudge,                &
              ud_hadv,ud_vadv,ud_hturb,ud_vturb,ud_pgrad,ud_rdamp,ud_nudge,   &
              ud_pbl,ud_cor,ud_cent,ud_lsw,                                   &
              vd_hadv,vd_vadv,vd_hturb,vd_vturb,vd_pgrad,vd_rdamp,vd_nudge,   &
              vd_pbl,vd_cor,vd_cent,vd_lsw,                                   &
              wd_hadv,wd_vadv,wd_hturb,wd_vturb,wd_pgrad,wd_rdamp,wd_buoy,    &
              kd_adv,kd_turb,                                                 &
              td_hidiff,td_vidiff,td_hediff,td_vediff,                        &
              qd_hidiff,qd_vidiff,qd_hediff,qd_vediff,                        &
              ud_hidiff,ud_vidiff,ud_hediff,ud_vediff,                        &
              vd_hidiff,vd_vidiff,vd_hediff,vd_vediff,                        &
              wd_hidiff,wd_vidiff,wd_hediff,wd_vediff,                        &
              iusekm,iusekh,nzeta,                                            &
              ib2pt,ie2pt,jb2pt,je2pt,kb2pt,ke2pt,ntwk,                       &
              ibnba,ienba,jbnba,jenba,kbnba,kenba,                            &
              bl_mynn_cloudpdf,icloud_bl,grav_settling,bl_mynn_tkebudget,     &
              bl_mynn_mixlength,bl_mynn_edmf,                                 &
              bl_mynn_edmf_mom,bl_mynn_edmf_tke,bl_mynn_mixscalars,           &
              bl_mynn_cloudmix,bl_mynn_mixqt,initflag,spp_pbl,bl_mynn_output, &
              nutk,nvtk,nwtk,                                                 &
              ntim,ntavr,utav,vtav,wtav,ttav,qtav,etav,uutav,vvtav,           &
              nsfctavr,les_subdomain_shape,ib_init

!-----------------------------------------------------------------------

      real dx,dy,dz,dtl,timax,run_time,                                       &
           kdiff2,kdiff6,fcor,kdiv,alph,rdalpha,zd,xhd,alphobc,umove,vmove,   &
           umovem1,umovem2,vmovem1,vmovem2,                                   &
           v_t,l_h,lhref1,lhref2,l_inf,ndcnst,nt_c,cnstce,cnstcd,             &
           dx_inner,dx_outer,nos_x_len,tot_x_len,                             &
           dy_inner,dy_outer,nos_y_len,tot_y_len,                             &
           ztop,str_bot,str_top,dz_bot,dz_top,                                &
           ptc_top,ptc_bot,viscosity,pr_num,                                  &
           rdx,rdy,rdz,rdx2,rdy2,rdz2,rdx4,rdy4,rdz4,                         &
           minx,maxx,miny,maxy,maxz,zt,rzt,pmin,                              &
           sfctheta,thec_mb,qt_mb,smeps,tsmall,qsmall,cflmax,ksmax,           &
           var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,                &
           var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,       &
           dtrad,ctrlat,ctrlon,                                               &
           tsk0,tmn0,xland0,oml_hml0,oml_gamma,                               &
           dmax,zdeep,lamx,lamy,xcent,ycent,aconv,convtime,                   &
           xc_uforce,xr_uforce,zr_uforce,alpha_uforce,t1_uforce,t2_uforce,    &
           xc_wnudge,xr_wnudge,zr_wnudge,alpha_wnudge,t1_wnudge,t2_wnudge,    &
           yc_wnudge,yr_wnudge,zc_wnudge,wmax_wnudge,                         &
           rxrwnudge,ryrwnudge,rzrwnudge,                                     &
           min_dx,min_dy,min_dz,max_dx,max_dy,max_dz,                         &
           cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,                                     &
           dgs1,dgs2,dgs3,dgt1,dgt2,dgt3,                                     &
           wbe1,wbe2,wbe3,wte1,wte2,wte3,                                     &
           csound,cstar,hurr_vg,hurr_rad,hurr_vgpl,hurr_rotate,hurr_angle,    &
           cnst_znt,cnst_shflx,cnst_lhflx,cnst_ust,ramp_time,                 &
           ddr,rlen,centerx,centery,                                          &
           gam_cp,gam_dr,                                                     &
           rad_solcon,rad_jday,rad_zenangle,rad_declin,rad_hrang,             &
           adapt_move_frq,                                                    &
           base_pbot,base_ptop,base_thbot,base_thtop,base_qvbot,base_qvtop,   &
           base_tbot,base_ttop,base_pibot,base_pitop,                         &
           les_subdomain_xlen,les_subdomain_ylen,                             &
           les_subdomain_dlen,les_subdomain_trnslen,rtim,                     &
           lsnudge_tau,lsnudge_start,lsnudge_end,lsnudge_ramp_time,           &
           top_cd,side_cd,outunitconv,                                        &
           recycle_width_dx,recycle_depth_m,                                  &
           recycle_cap_loc_m,recycle_inj_loc_m,bl_mynn_closure

!-----------------------------------------------------------------------

      double precision ::                                                     &
           tapfrq,rstfrq,statfrq,prclfrq,diagfrq,azimavgfrq,hifrqfrq

!-----------------------------------------------------------------------

      character(len=maxstring) :: string
      character(len=12) :: aunit
      character(len=12) :: cm1version
      real :: cm1rversion

!-----------------------------------------------------------------------

      integer, parameter :: sp = kind(1.0)
      integer, parameter :: dp = max( sp , selected_real_kind(2*precision(1.0_sp)) )
      integer, parameter :: qp = max( dp , selected_real_kind(2*precision(1.0_dp)) )

      real(kind=sp) :: pi_sp
      real(kind=dp) :: pi_dp
      real(kind=qp) :: pi_qp

!-----------------------------------------------------------------------

      logical :: nssl_3moment = .false.
      logical :: nssl_density_on = .true.

      integer :: output_n_as_m3 = 0   ! convert number conc from #/kg to #/m3
      integer :: output_uv_ground = 0  ! output u,v, as ground relative

!-----------------------------------------------------------------------
!  timestats:

      real clock_rate,time_sound,time_buoyan,time_turb,              &
           time_diffu,time_microphy,time_stat,time_cflq,time_bc,     &
           time_misc,time_integ,time_rdamp,time_divx,time_write,     &
           time_ttend,time_cor,time_fall,time_satadj,                &
           time_sfcphys,time_dbz,time_last,                          &
           time_advs,time_advu,time_advv,time_advw,                  &
           time_mpu1,time_mpv1,time_mpw1,time_mpp1,                  &
           time_mpu2,time_mpv2,time_mpw2,time_mpp2,                  &
           time_mps1,time_mps3,time_mpq1,time_mptk1,                 &
           time_mps2,time_mps4,time_mpq2,time_mptk2,time_mpb,        &
           time_parcels,time_rad,time_pbl,time_swath,time_pdef,      &
           time_prsrho,time_restart,time_poiss,time_diag,            &
           time_azimavg,time_hifrq,time_ercyl,time_tavg,time_cuprm
      integer count_last

!-----------------------------------------------------------------------
!  Begin namelists:

      namelist /param0/ nx,ny,nz,ppnode,                                &
                        timeformat,timestats,terrain_flag,procfiles,    &
                        outunits
      namelist /param1/ dx,dy,dz,dtl,timax,run_time,                    &
                        tapfrq,rstfrq,statfrq,prclfrq
      namelist /param2/                                                 &
          cm1setup,testcase,adapt_dt,irst,rstnum,iconly,                &
          hadvordrs,vadvordrs,hadvordrv,vadvordrv,                      &
          advwenos,advwenov,weno_order,                                 &
          apmasscon,idiff,mdiff,difforder,imoist,                       &
          ipbl,sgsmodel,tconfig,bcturbs,horizturb,doimpl,               &
          irdamp,hrdamp,psolver,ptype,ihail,iautoc,cuparam,             &
          icor,betaplane,lspgrad,eqtset,idiss,efall,rterm,              &
          wbc,ebc,sbc,nbc,bbc,tbc,irbc,roflux,nudgeobc,                 &
          isnd,iwnd,itern,iinit,                                        &
          irandp,ibalance,iorigin,axisymm,imove,iptra,npt,pdtra,        &
          iprcl,nparcels
      namelist /param3/ kdiff2,kdiff6,fcor,kdiv,alph,rdalpha,zd,xhd,alphobc, &
                        umove,vmove,v_t,l_h,lhref1,lhref2,l_inf,ndcnst,nt_c, &
                        csound,cstar
      namelist /param4/ stretch_x,dx_inner,dx_outer,nos_x_len,tot_x_len
      namelist /param5/ stretch_y,dy_inner,dy_outer,nos_y_len,tot_y_len
      namelist /param6/ stretch_z,ztop,str_bot,str_top,dz_bot,dz_top
      namelist /param7/ bc_temp,ptc_top,ptc_bot,viscosity,pr_num
      namelist /param8/ var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,   &
                        var11,var12,var13,var14,var15,var16,var17,var18,var19,var20
      namelist /param9/                                                       &
              output_format,output_filetype,                                  &
              output_interp,                                                  &
              output_rain,output_sws,output_svs,output_sps,output_srs,        &
              output_sgs,output_sus,output_shs,output_coldpool,               &
              output_sfcflx,output_sfcparams,output_sfcdiags,                 &
              output_psfc,output_zs,output_zh,output_basestate,               &
              output_th,output_thpert,output_prs,output_prspert,              &
              output_pi,output_pipert,output_rho,output_rhopert,output_tke,   &
              output_km,output_kh,                                            &
              output_qv,output_qvpert,output_q,output_dbz,output_buoyancy,    &
              output_u,output_upert,output_uinterp,                           &
              output_v,output_vpert,output_vinterp,output_w,output_winterp,   &
              output_vort,output_pv,output_uh,output_pblten,                  &
              output_dissten,output_fallvel,output_nm,output_def,             &
              output_radten,output_cape,output_cin,output_lcl,output_lfc,     &
              output_pwat,output_lwp,                                         &
              output_thbudget,output_qvbudget,                                &
              output_ubudget,output_vbudget,output_wbudget,output_pdcomp
      namelist /param16/                                                      &
              restart_format,restart_filetype,restart_reset_frqtim,           &
              restart_file_theta,restart_file_dbz,restart_file_th0,           &
              restart_file_prs0,restart_file_pi0,restart_file_rho0,           &
              restart_file_qv0,restart_file_u0,restart_file_v0,               &
              restart_file_zs,restart_file_zh,restart_file_zf,                &
              restart_file_diags,restart_use_theta
      namelist /param10/                                                      &
              stat_w,stat_wlevs,stat_u,stat_v,stat_rmw,                       &
              stat_pipert,stat_prspert,stat_thpert,stat_q,                    &
              stat_tke,stat_km,stat_kh,stat_div,stat_rh,stat_rhi,stat_the,    &
              stat_cloud,stat_sfcprs,stat_wsp,stat_cfl,stat_vort,             &
              stat_tmass,stat_tmois,stat_qmass,stat_tenerg,stat_mo,stat_tmf,  &
              stat_pcn,stat_qsrc
      namelist /param11/                                                      &
              radopt,dtrad,ctrlat,ctrlon,year,month,day,hour,minute,second
      namelist /param12/                                                      &
              isfcflx,sfcmodel,oceanmodel,initsfc,                            &
              tsk0,tmn0,xland0,lu0,season,                                    &
              cecd,pertflx,cnstce,cnstcd,                                     &
              isftcflx,iz0tlnd,oml_hml0,oml_gamma,                            &
              set_flx,cnst_shflx,cnst_lhflx,                                  &
              set_znt,cnst_znt,set_ust,cnst_ust,                              &
              ramp_sgs,ramp_time,t2p_avg    
      namelist /param13/                                                      &
              prcl_th,prcl_t,prcl_prs,prcl_ptra,prcl_q,prcl_nc,               &
              prcl_km,prcl_kh,prcl_tke,prcl_dbz,prcl_b,prcl_vpg,prcl_vort,    &
              prcl_rho,prcl_qsat,prcl_sfc
      namelist /param14/                                                      &
              dodomaindiag,diagfrq
      namelist /param15/                                                      &
              doazimavg,azimavgfrq,rlen,do_adapt_move,adapt_move_frq
      namelist /param17/                                                      &
              les_subdomain_shape,les_subdomain_xlen,les_subdomain_ylen,      &
              les_subdomain_dlen,les_subdomain_trnslen
      namelist /param18/                                                      &
              do_recycle_w,do_recycle_s,do_recycle_e,do_recycle_n,            &
              recycle_width_dx,recycle_depth_m,                               &
              recycle_cap_loc_m,recycle_inj_loc_m
      namelist /param19/                                                         &
              do_lsnudge,do_lsnudge_u,do_lsnudge_v,do_lsnudge_th,do_lsnudge_qv,  &
              lsnudge_tau,lsnudge_start,lsnudge_end,lsnudge_ramp_time
      namelist /param20/                                                      &
              do_ib,ib_init,top_cd,side_cd
      namelist /param21/                                                      &
              hurr_vg,hurr_rad,hurr_vgpl,hurr_rotate

!  End namelists:

    CONTAINS

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real function mytime()
      implicit none

      integer count,rate,max
      real time_current,rcount

      !  Platform-independent timer

      call system_clock(count,rate,max)
      if( count.lt.count_last )then
        ! simple kludge ... do nothing
        ! fix some other day   (GHB, 101018)
!!!        rcount = float(count+max)
!!!        time_current=rcount*clock_rate
!!!        mytime=time_current-time_last
!!!        rcount = float(count)
!!!        time_current=rcount*clock_rate
        rcount = float(count)
        time_current=rcount*clock_rate
        mytime=0.0
      else
        rcount = float(count)
        time_current=rcount*clock_rate
        mytime=time_current-time_last
      endif
      time_last=time_current
      count_last=count

      end function mytime

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE input
