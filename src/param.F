  MODULE param_module

  implicit none

  private
  public :: param,wenocheck

  CONTAINS

      subroutine param(dt,dtlast,stattim,taptim,rsttim,radtim,prcltim,  &
                       cloudvar,rhovar,qmag,qname,qunit,budname,        &
                       xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf, &
                       yh,vh,rvh,yf,vf,rvf,xfref,xhref,yfref,yhref,     &
                       rds,sigma,rdsf,sigmaf,tauh,taus,                 &
                       zh,mh,rmh,cc1,cc2,tauf,zf,mf,rmf,f2d,dtu0,dtv0,  &
                       gamk,gamwall,                                    &
                       zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,  &
                       reqs_u,reqs_v,reqs_s,reqs_p,                     &
                       nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                       n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,         &
                       sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,         &
                       uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                       vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                       ww31,ww32,we31,we32,ws31,ws32,wn31,wn32)
      use input
      use constants
      use index_module
      use init_terrain_module
      use bc_module
      use comm_module
      use ib_module
      use eddy_recycle
      use lsnudge_module
#ifdef MPI
      use mpi
#endif
      use module_mp_thompson , only : thompson_init
      use module_mp_graupel
      use mp_wsm6 , only : mp_wsm6_init
      use microphy_p3, only : p3_init
      use module_mp_jensen_ishmael, only : jensen_ishmael_init
      use module_mp_nssl_2mom, only: nssl_2mom_init, rho_qr, cnor, rho_qs,        &
                        cnos, rho_qh, rho_qhl, cnoh, ccn, irenuc, infall, alphah, alphahl, &
                        imurain, icdx, icdxhl, dfrz, hldnmn, iferwisventr,        &
                        iehw, iehlw,ehw0, ehlw0, dmrauto, inucopt, ioldlimiter, ihlcnh
      use goddard_module, only : consat,consat2
      use lfoice_module, only : lfoice_init
      use ccpp_kind_types,only: kind_phys
      implicit none

      real :: dt,dtlast
      double precision :: stattim,taptim,rsttim,radtim,prcltim
      logical, dimension(maxq) :: cloudvar,rhovar
      real, intent(inout), dimension(maxq) :: qmag
      character(len=3), dimension(maxq) :: qname
      character(len=20), dimension(maxq) :: qunit
      character(len=6), dimension(maxq) :: budname
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: yf,vf,rvf
      real, dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
      real, dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
      real, dimension(kb:ke) :: rds,sigma
      real, dimension(kb:ke+1) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,cc1,cc2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, intent(inout), dimension(ib:ie,jb:je) :: f2d
      real, intent(inout), dimension(ib:ie,jb:je) :: dtu0,dtv0
      real, intent(inout), dimension(kb:ke) :: gamk,gamwall
      real, dimension(ib:ie,jb:je) :: zs
      real, dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s,reqs_p
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32

!-----------------------------------------------------------------------

      integer i,j,k,n,m,nn,kst,ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3,vadv,p3stat
      integer ival,jval,ii,jj,ltmp
      integer iterrain
      integer :: inum
      real :: var
      real :: zfw1d(kb:ke+1), zfs1d(kb:ke+1)
      double precision :: gzc(kb:ke+1), gze(kb:ke+1)
      integer :: nbndlyr    = 0
      real    :: rtop       = 1.0     ! upper level stretch factor
      real    :: ztopstr    = 100000. ! height to start upper level stretching
      real    :: dzmaxtop   = 700.    ! max upper level dz
#ifdef MPI
      character(len=50) :: fname
#endif
      real :: c1,c2,nominal_dx,nominal_dy,nominal_dz,z1,z2,z3,mult
      real :: x1,x2,y1,y2,phi,beta0,max_z_out,test_len
      logical :: getvt,dothis,hrdamp_west,hrdamp_east,hrdamp_south,hrdamp_north
      double precision, dimension(:), allocatable :: xfdp,yfdp

      integer, parameter :: bigm = 4   ! highest deriv
      integer, parameter :: bign = 3   ! number of grid points minus 1
      double precision :: x0,b1,b2,b3
      double precision, dimension(0:bign) :: alpha
      ! delta(n,m,nu):
      double precision, dimension(0:bign,-1:bigm,0:bign) :: delta

      character(len=60) :: param_mp,param_sgs,param_pbl,param_sfclay,param_rad

      logical :: doit,getfall
      logical :: nssl_hail_on = .false.
      character(len=40) :: errmsg
      integer :: errflg
      real(kind=kind_phys) :: pi_phys

#ifdef MPI
      integer :: ntmp1,ntmp2,ntmp3,ntmp4,reqs,reqs1,reqs2,reqs3,reqs4
      integer, dimension(MPI_STATUS_SIZE) :: status
#endif

!-----------------------------------------------------------------------
!  for nssl microphysics:

      NAMELIST /nssl2mom_params/            &
                        rho_qr,         &
                        cnor,           &
                        rho_qs,         &
                        cnos,           &
                        rho_qh, rho_qhl,&
                        cnoh,           &
                        ccn,            &
                        irenuc,         &
                        infall,         &
                        alphah,         &
                        alphahl,        &
                        imurain,        &
                        icdx,           &
                        icdxhl,         &
                        dfrz,           &
                        hldnmn,         &
                        iferwisventr,   &
                        iehw,iehlw,     &
                        ehw0,ehlw0,     &
                        inucopt,        &
                        dmrauto,        &
                        ioldlimiter,    &
                        ihlcnh

!--------------------------------------------------------------

#ifdef MPI
      IF(procfiles)THEN
        fname='procXXXXXX.print.out'
        write(fname(5:10),100) myid
100     format(i6.6)
        open(unit=10,file=fname,status='unknown')
      ENDIF
#endif

      if(dowr) write(outfile,*) 'Inside PARAM'

#ifndef MPI
!--------------------------------------------------------------

      if(nodex.ne.1 .or. nodey.ne.1)then
        print *
        print *,'  For non-MPI runs, nodex and nodey must be = 1 !'
        print *
        call stopcm1
      endif

#endif
!--------------------------------------------------------------

#ifdef MPI
      myid0:  &
      if(myid.eq.0)then
#endif

      open(unit=20,file='namelist.input',form='formatted',status='old',    &
           access='sequential',err=8000)

      read(20,nml=param1,end=701)
701   continue
      rewind(20)

      read(20,nml=param2,end=702)
702   continue
      rewind(20)

      read(20,nml=param3,end=703)
703   continue
      rewind(20)

      read(20,nml=param11,end=711)
711   continue
      rewind(20)

      read(20,nml=param12,end=712)
712   continue
      rewind(20)

      read(20,nml=param4,end=704)
704   continue
      rewind(20)

      read(20,nml=param5,end=705)
705   continue
      rewind(20)

      read(20,nml=param6,end=706)
706   continue
      rewind(20)

      read(20,nml=param7,end=707)
707   continue
      rewind(20)

      read(20,nml=param8,end=708)
708   continue
      rewind(20)

      read(20,nml=param9,end=709)
709   continue
      rewind(20)

      read(20,nml=param16,end=716)
716   continue
      rewind(20)

      read(20,nml=param10,end=710)
710   continue
      rewind(20)

      if( iprcl.eq.1 )then
        read(20,nml=param13,end=713)
713     continue
        rewind(20)
      endif

      read(20,nml=param14,end=714)
714   continue
      rewind(20)

      read(20,nml=param15,end=715)
715   continue
      rewind(20)

      IF ( ptype.eq.26 .or. ptype.eq.27 .or. ptype.eq.28 ) THEN
         read(20,nml=nssl2mom_params,end=751)
751      continue
         rewind(20)
      ENDIF

      read(20,nml=param17,end=717)
717   continue
      rewind(20)

      read(20,nml=param18,end=718)
718   continue
      rewind(20)

      read(20,nml=param19,end=719)
719   continue
      rewind(20)

      read(20,nml=param20,end=720)
720   continue
      rewind(20)

      read(20,nml=param21,end=721)
721   continue
      rewind(20)

      close(unit=20)

#ifdef MPI
      endif  myid0

      call MPI_BCAST(dx    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dy    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dz    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtl   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(timax ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(run_time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tapfrq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rstfrq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(statfrq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prclfrq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(cm1setup ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(testcase,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adapt_dt ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irst     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rstnum   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iconly   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hadvordrs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vadvordrs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hadvordrv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vadvordrv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(advwenos ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(advwenov ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(weno_order,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(apmasscon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idiff    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mdiff    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(difforder,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(imoist   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ipbl     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sgsmodel ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tconfig  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bcturbs  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(horizturb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(doimpl   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irdamp   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hrdamp   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(psolver  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptype    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ihail    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iautoc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cuparam  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icor     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(betaplane,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lspgrad  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(eqtset   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idiss    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(efall    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rterm    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ebc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irbc     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(roflux   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nudgeobc ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(isnd     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iwnd     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(itern    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iinit    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irandp   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibalance ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iorigin  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(axisymm  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(imove    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iptra    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npt      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pdtra    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iprcl    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nparcels ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(kdiff2 ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kdiff6 ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fcor   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kdiv   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(alph   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rdalpha,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(zd     ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xhd    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(alphobc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umove  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmove  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(v_t    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(l_h    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lhref1 ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lhref2 ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(l_inf  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ndcnst ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nt_c   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(csound ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cstar  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stretch_x,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dx_inner ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dx_outer ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nos_x_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tot_x_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stretch_y,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dy_inner ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dy_outer ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nos_y_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tot_y_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stretch_z,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ztop     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(str_bot  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(str_top  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dz_bot   ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dz_top   ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(bc_temp  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptc_top  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptc_bot  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(viscosity,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pr_num   ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(var1     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var2     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var3     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var4     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var5     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var6     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var7     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var8     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var9     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var10    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var11    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var12    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var13    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var14    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var15    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var16    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var17    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var18    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var19    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var20    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(output_format  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_filetype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_interp ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_rain   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sws    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_svs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sps    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_srs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sgs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sus    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_shs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_coldpool,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sfcflx ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sfcparams,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sfcdiags,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_psfc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_zs     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_zh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_basestate,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_th     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_thpert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_prs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_prspert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pi     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pipert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_rho    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_rhopert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_tke    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_km     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_kh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_qv     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_qvpert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_q      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_n_as_m3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_dbz    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_buoyancy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_u      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_upert  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_uinterp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_v      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vpert  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vinterp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_w      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_winterp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vort   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pv     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_uh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pblten ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_dissten,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_fallvel ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_nm     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_def    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_radten ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_cape   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_cin    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_lcl    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_lfc    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pwat   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_lwp    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_thbudget,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_qvbudget,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_ubudget,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vbudget,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_wbudget,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pdcomp ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(restart_format      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_filetype    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_theta  ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_dbz    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_th0    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_prs0   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_pi0    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_rho0   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_qv0    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_u0     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_v0     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_zs     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_zh     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_zf     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_file_diags  ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_use_theta   ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart_reset_frqtim,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stat_w      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_wlevs  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_u      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_v      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_rmw    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_pipert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_prspert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_thpert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_q      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tke    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_km     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_kh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_div    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_rh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_rhi    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_the    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_cloud  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_sfcprs ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_wsp    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_cfl    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_vort   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tmass  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tmois  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_qmass  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tenerg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_mo     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tmf    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_pcn    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_qsrc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(prcl_th     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_t      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_prs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_ptra   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_q      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_nc     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_km     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_kh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_tke    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_dbz    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_b      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_vpg    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_vort   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_rho    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_qsat   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcl_sfc    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(radopt   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtrad  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ctrlat ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ctrlon ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(year   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(month  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(day    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hour   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(minute ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(second ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(isfcflx   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sfcmodel  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oceanmodel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(initsfc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tsk0      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tmn0      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xland0    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lu0       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(season    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cecd      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pertflx   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnstce    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnstcd    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(isftcflx  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iz0tlnd   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oml_hml0  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oml_gamma ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(set_flx   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(set_znt   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(set_ust   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnst_shflx,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnst_lhflx,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnst_znt  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnst_ust  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ramp_sgs  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ramp_time ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(t2p_avg   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(dodomaindiag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(diagfrq   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(doazimavg     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(azimavgfrq    ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rlen          ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_adapt_move ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adapt_move_frq,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(les_subdomain_shape   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(les_subdomain_xlen    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(les_subdomain_ylen    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(les_subdomain_dlen    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(les_subdomain_trnslen ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(do_recycle_w,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_recycle_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_recycle_e,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_recycle_n,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recycle_width_dx ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recycle_depth_m  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recycle_cap_loc_m,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recycle_inj_loc_m,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(do_lsnudge       ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_lsnudge_u     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_lsnudge_v     ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_lsnudge_th    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(do_lsnudge_qv    ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_tau      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_start    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_end      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lsnudge_ramp_time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST( do_ib           ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( ib_init         ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( top_cd          ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( side_cd         ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST( hurr_vg         ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( hurr_rad        ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( hurr_vgpl       ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( hurr_rotate     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      IF ( ptype.eq.26 .or. ptype.eq.27 .or. ptype.eq.28 ) THEN
      call MPI_BCAST(ccn       ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho_qr    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnor      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho_qs    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnos      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho_qh    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho_qhl   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnoh      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(alphah    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(alphahl   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dfrz      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hldnmn    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irenuc    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(infall    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icdx      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icdxhl    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(imurain   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iferwisventr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iehw      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iehlw     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ehw0      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ehlw0     ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(inucopt   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dmrauto   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ioldlimiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ihlcnh    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nssl_3moment,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nssl_density_on,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ENDIF
#endif

!-----------------------------------------------------------------------
!  Set constants:

      call set_constants( testcase )

!-----------------------------------------------------------------------
!  Some "dummy" checks:

      dt = dtl
      dtlast = dtl
      eqtset = max( eqtset , 1 )
      eqtset = min( eqtset , 2 )
      axisymm = max( axisymm , 0 )
      axisymm = min( axisymm , 1 )
      imoist = max( imoist , 0 )
      imoist = min( imoist , 1 )
      if(imoist.ne.1) efall=0

      icor = max( icor , 0 )
      icor = min( icor , 1 )

      if(icor.eq.0) betaplane = 0
      betaplane = max( betaplane , 0 )
      betaplane = min( betaplane , 1 )

      ! Cannot use beta plane for axisymmetric simulations, 
      ! or 2d (x,z) simulations:
      if( axisymm.eq.1 .or. ny.le.3 )  betaplane = 0

      if(imove.eq.0) umove=0.0
      if(imove.eq.0) vmove=0.0
      if(axisymm.eq.1) do_adapt_move=.false.
      if( do_adapt_move )then
        ! set to zero by default
        umove = 0.0
        vmove = 0.0
        umovem1 = umove
        umovem2 = umove
        vmovem1 = vmove
        vmovem2 = vmove
        imove = 1
        if( axisymm.eq.1 )then
          ! turn off for axisymmetric runs:
          imove = 0
        endif
      endif
      irst = max( irst , 0 )
      irst = min( irst , 1 )
      if( cm1setup.eq.2 ) tconfig = 2
      if( psolver.eq.4 .or. psolver.eq.5 ) roflux = 1
      if( cstar.lt.0.1 ) cstar = 30.0

      if( irdamp.eq.0 .and. hrdamp.eq.0 ) rdalpha = 0.0
      if( irdamp.ge.1 .and. hrdamp.ge.1 ) hrdamp = irdamp

      pdcomp = .false.
      output_pdcomp = max(0,min(1,output_pdcomp))
      if( output_pdcomp.eq.1 ) pdcomp = .true.

      ! cm1r19: hard-wire pdscheme
      pdscheme = 1

      call wenocheck

      if( axisymm.eq.1 )then
        dodomaindiag = .false.
        doazimavg = .false.
        testcase = 0
      endif

      IF( do_lsnudge )THEN
        ! for diagnostic purposes:
        if( .not. dodomaindiag )then
          dodomaindiag = .true.
          diagfrq = tapfrq
        endif
      ENDIF

      ramp_time = max( ramp_time , 1.0e-10 )

      if( cm1setup.eq.1 .or. cm1setup.eq.4 )  idoles = .true.
      if( cm1setup.eq.2 .or. cm1setup.eq.4 )  idopbl = .true.

      !------------

      IF( nk.lt.8 )THEN
        if(myid.eq.0)then
        print *
        print *,'  nk   = ',nk
        print *
        print *,'  nk must be >= 8 '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

      !------------

      vadv = max( hadvordrv , vadvordrv )

      IF( (vadv.eq.3.or.vadv.eq.4) .and. nk.lt.4 )THEN
        if(myid.eq.0)then
        print *
        print *,'  nk   = ',nk
        print *,'  vadv = ',vadv
        print *
        print *,'  nk must be >= 4 for vadv=3,4 '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      IF( (vadv.eq.5.or.vadv.eq.6) .and. nk.lt.6 )THEN
        if(myid.eq.0)then
        print *
        print *,'  nk   = ',nk
        print *,'  vadv = ',vadv
        print *
        print *,'  nk must be >= 6 for vadv=5,6 '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      IF( (vadv.eq.7.or.vadv.eq.8) .and. nk.lt.8 )THEN
        if(myid.eq.0)then
        print *
        print *,'  nk   = ',nk
        print *,'  vadv = ',vadv
        print *
        print *,'  nk must be >= 8 for vadv=7,8 '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      IF( (vadv.eq.9.or.vadv.eq.10) .and. nk.lt.10 )THEN
        if(myid.eq.0)then
        print *
        print *,'  nk   = ',nk
        print *,'  vadv = ',vadv
        print *
        print *,'  nk must be >= 10 for vadv=9,10 '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

!------------------------------------------------
!  begin non-fatal checks:

      !-----
      IF( isnd.eq.7 )THEN
        IF( iwnd.gt.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  isnd = ',isnd
        print *
        print *,'  but iwnd = ',iwnd
        print *
        print *,'  setting iwnd to 0; using wind profile from input_sounding '
        print *
        print *,'  (if you want to use iwnd >= 1 with a thermodynamic sounding '
        print *,'   from input_sounding, then set isnd = 17)'
        print *
        print *,'  -------------------------------- '
        endif
        iwnd = 0
        ENDIF
      ENDIF
      !-----
      IF( sfcmodel.ge.1 )THEN
        IF( oceanmodel.le.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  but oceanmodel = ',oceanmodel
        print *
        print *,'  setting oceanmodel to 1  (just in case its needed) '
        print *
        print *,'  -------------------------------- '
        endif
        oceanmodel = 1
        ENDIF
      ENDIF
      !-----
      IF( (sfcmodel.ge.1) .or. (oceanmodel.eq.2) .or. (ipbl.ge.1) )then
        IF( bbc.ne.3 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  oceanmodel = ',oceanmodel
        print *,'  ipbl       = ',ipbl
        print *
        print *,'  at least one of these options requires bbc = 3 '
        print *,'  ... so, setting bbc to 3 '
        print *
        print *,'  -------------------------------- '
        endif
        bbc = 3
        ENDIF
      ENDIF
      !-----
      IF( irst.eq.1 .and. iinit.ne.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  irst       = ',irst
        print *,'  iinit      = ',iinit
        print *
        print *,'  This is a restart. '
        print *,'  so, setting iinit to 0 '
        print *
        print *,'  -------------------------------- '
        endif
        iinit = 0
      ENDIF
      !-----
      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4.or.sfcmodel.eq.5.or.sfcmodel.eq.6.or.sfcmodel.eq.7) .and. isfcflx.eq.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  isfcflx    = ',isfcflx
        print *
        print *,'  sfcmodel=2,3,4,5,6,7 requires isfcflx=1 '
        print *,'  so, setting isfcflx to 1 '
        print *
        print *,'  -------------------------------- '
        endif
        isfcflx = 1
      ENDIF
      !-----
      IF( (psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6.or.psolver.eq.7) .and. eqtset.eq.2 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  psolver    = ',psolver
        print *,'  eqtset     = ',eqtset
        print *
        print *,'  psolver=4,5,6,7 requires eqtset=1 '
        print *,'  ... setting eqtset to 1 ... '
        print *
        print *,'  -------------------------------- '
        endif
        eqtset = 1
      ENDIF
      !-----
      IF( (ptype.eq.4.or.ptype.eq.7.or.ptype.eq.50.or.ptype.eq.51.or.ptype.eq.52.or.ptype.eq.53.or.ptype.eq.55) .and. eqtset.eq.2 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  ptype      = ',ptype
        print *,'  eqtset     = ',eqtset
        print *
        print *,'  ptype=4,7,50,51,52,53,55 requires eqtset=1 '
        print *
        print *,'  ... setting eqtset to 1 ... '
        print *
        print *,'  -------------------------------- '
        endif
        eqtset = 1
      ENDIF
      !-----
      IF( (psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6.or.psolver.eq.7) .and. apmasscon.eq.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  psolver    = ',psolver
        print *,'  apmasscon  = ',apmasscon
        print *
        print *,'  psolver=4,5,6,7 requires apmasscon=0 '
        print *,'  ... setting apmasscon to 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        apmasscon = 0
      ENDIF
      !-----
      IF( restart_use_theta )THEN
      IF( .not. restart_file_theta )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  restart_use_theta = ',restart_use_theta
        print *
        print *,'  ... setting restart_file_theta to true ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_file_theta = .true.
      ENDIF
      ENDIF
      !-----
#ifndef MPI
      IF( restart_filetype.ge.3 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  Single processor run: '
        print *
        print *,'  restart_filetype >= 3 not available '
        print *
        print *,'  ... setting restart_filetype to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_filetype = 2
      ENDIF
#endif
      !-----
      IF( restart_format.eq.1 .and. restart_filetype.eq.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  binary-format restart file: '
        print *
        print *,'  restart_filetype = 1 not available '
        print *
        print *,'  ... setting restart_filetype to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_filetype = 2
      ENDIF
      !-----
#ifdef NETCDF
      IF( restart_format.eq.2 .and. restart_filetype.ge.3 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  netcdf-format restart file: '
        print *
        print *,'  restart_filetype >= 3 not available '
        print *
        print *,'  ... setting restart_filetype to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        restart_filetype = 2
      ENDIF
#endif
      !-----
      IF( idopbl .and. (ipbl.eq.1.or.ipbl.eq.3.or.ipbl.eq.4.or.ipbl.eq.5.or.ipbl.eq.6) .and. abs(l_inf).gt.1.0e-6 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  ipbl = 1,3,4,5,6  require  l_inf = 0 '
        print *
        print *,'  ... setting l_inf to 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        l_inf = 0.0
      ENDIF
      !-----
      IF( testcase.eq.1 .or. testcase.eq.2 )THEN
        if( imoist.ne.0 )then
          if(myid.eq.0)then
          print *,'  -------------------------------- '
          print *
          print *,'  testcase = ',testcase
          print *
          print *,'  ... setting imoist to 0 ... '
          print *
          print *,'  -------------------------------- '
          endif
          imoist = 0
        endif
      ENDIF
      !-----
      IF( testcase.eq.3 .or. testcase.eq.4 .or. testcase.eq.5 .or. testcase.eq.7 .or. testcase.eq.8 )THEN
        if( imoist.ne.1 )then
          if(myid.eq.0)then
          print *,'  -------------------------------- '
          print *
          print *,'  testcase = ',testcase
          print *
          print *,'  ... setting imoist to 1 ... '
          print *
          print *,'  -------------------------------- '
          endif
          imoist = 1
        endif
      ENDIF
      !-----
      IF( ( .not. idoles ) .and. sgsmodel.ge.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  idoles = ',idoles
        print *
        print *,'  ... setting sgsmodel to 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        sgsmodel = 0
      ENDIF
      !-----
      IF( cm1setup.ge.1 .and. idiff.ge.1 .and. difforder.eq.2 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  cm1setup  = ',cm1setup
        print *,'  idiff     = ',idiff
        print *,'  difforder = ',difforder
        print *
        print *,'  Cannot use difforder = 2 with cm1setup >= 1 '
        print *
        print *,'  ... setting idiff = 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        idiff = 0
      ENDIF
      !-----
      IF( (ipbl.eq.3.or.ipbl.eq.4.or.ipbl.eq.5) .and. idiss.ge.1 .and. cm1setup.ne.4  )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  ipbl      = ',ipbl
        print *,'  idiss     = ',idiss
        print *
        print *,'  dissipative heating is calculated within ipbl=3,4,5 '
        print *
        print *,'  ... setting idiss = 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        idiss = 0
      ENDIF
      !-----
      IF( ( sgsmodel.eq.5 .or. sgsmodel.eq.6 ) .and. doimpl.ge.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  sgsmodel  = ',sgsmodel
        print *,'  doimpl    = ',doimpl
        print *
        print *,'  cannot use doimpl=1 with sgsmodel=5,6 (for now) '
        print *
        print *,'  ... setting doimpl = 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        doimpl = 0
      ENDIF
      !-----
      IF( do_lsnudge .and. (do_lsnudge_u .or. do_lsnudge_v .or. do_lsnudge_th) .and. irdamp.eq.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  do_lsnudge       =  ',do_lsnudge
        print *,'  do_lsnudge_u     =  ',do_lsnudge_u
        print *,'  do_lsnudge_v     =  ',do_lsnudge_v
        print *,'  do_lsnudge_th    =  ',do_lsnudge_th
        print *,'  irdamp           =  ',irdamp
        print *
        print *,'  when using lsnudge, irdamp = 2 is recommended '
        print *
        print *,'  ... setting irdamp to 2 ... '
        print *
        print *,'  -------------------------------- '
        endif
        irdamp = 2
      ENDIF
      !-----
      IF( testcase.eq.11 .and. imoist.eq.0 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  testcase         =  ',testcase
        print *,'  imoist           =  ',imoist
        print *
        print *,'  imoist = 1 is required for testcase = 11 '
        print *
        print *,'  ... setting imoist to 1 ... '
        print *
        print *,'  -------------------------------- '
        endif
        imoist = 1
      ENDIF
      !-----

!  end non-fatal checks:
!------------------------------------------------
!  begin fatal checks  (ie, model stops)

      !-----
      IF( cm1setup.lt.0 .or. cm1setup.gt.4 )THEN
        if(myid.eq.0)then
        print *
        print *,'  cm1setup = ',cm1setup
        print *
        print *,'  cm1setup must be either 0, 1, 2, 3, 4'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( hadvordrs.lt.2 .or. hadvordrs.gt.10 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  hadvordrs = ',hadvordrs
        print *
        print *,'  This value is invalid (must be between 2 and 10) '
        print *
        print *,'  -------------------------------- '
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( hadvordrv.lt.2 .or. hadvordrv.gt.10 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  hadvordrv = ',hadvordrv
        print *
        print *,'  This value is invalid (must be between 2 and 10) '
        print *
        print *,'  -------------------------------- '
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( vadvordrs.lt.2 .or. vadvordrs.gt.10 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  vadvordrs = ',vadvordrs
        print *
        print *,'  This value is invalid (must be between 2 and 10) '
        print *
        print *,'  -------------------------------- '
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( vadvordrv.lt.2 .or. vadvordrv.gt.10 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  vadvordrv = ',vadvordrv
        print *
        print *,'  This value is invalid (must be between 2 and 10) '
        print *
        print *,'  -------------------------------- '
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( advwenos.ge.1 .or. advwenov.ge.1 )THEN
      IF( weno_order.ne.3 .and. weno_order.ne.5 .and. weno_order.ne.7 .and. weno_order.ne.9 )THEN
        if(myid.eq.0)then
        print *
        print *,'  weno_order = ',weno_order
        print *
        print *,'  invalid value for weno_order '
        print *
        print *,'  (weno_order must be 3,5,7,9)'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF( psolver.lt.1 .or. psolver.gt.7 )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver  = ',psolver
        print *
        print *,'  invalid value for psolver '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( psolver.eq.1 .and. adapt_dt.eq.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver  = ',psolver
        print *,'  adapt_dt = ',adapt_dt
        print *
        print *,'  Cannot use adapt_dt with psolver=1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
!      IF( sfcmodel.eq.1 .or. sfcmodel.eq.5 )THEN
!      IF( set_znt.eq.1 .and. set_ust.eq.1 )THEN
!        if(myid.eq.0)then
!        print *
!        print *,'  set_znt = ',set_znt
!        print *,'  set_ust = ',set_ust
!        print *
!        print *,'  cannot use set_znt=1 and set_ust=1 at the same time '
!        print *,'  (one or the other must be zero to zero) '
!        print *
!        print *,'   stopping model .... '
!        print *
!        endif
!#ifdef MPI
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
!        call stopcm1
!      ENDIF
!      ENDIF
      !-----
      IF( sfcmodel.eq.1 .or. sfcmodel.eq.5 )THEN
      IF( set_znt.gt.1 .or. set_ust.gt.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  set_znt = ',set_znt
        print *,'  set_ust = ',set_ust
        print *
        print *,'  set_znt and set_ust cannot be greater than 1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF(sfcmodel.eq.1)THEN
      IF(cecd.lt.1.or.cecd.gt.3)THEN
        if(myid.eq.0)then
        print *
        print *,'  cecd  = ',cecd
        print *
        print *,'  cecd must be 1,2,3'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      !-----
      IF( tbc.eq.3 .and. sfcmodel.ne.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  tbc       = ',tbc
        print *,'  sfcmodel  = ',sfcmodel
        print *
        print *,'  tbc = 3  requires  sfcmodel = 1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( idoles )THEN
      IF( sgsmodel.lt.0 .or. sgsmodel.gt.6 )THEN
        if(myid.eq.0)then
        print *
        print *,'  idoles   = ',idoles
        print *,'  cm1setup = ',cm1setup
        print *,'  sgsmodel  = ',sgsmodel
        print *
        print *,'  for idoles, sgsmodel must be either 0,1,2,3,4,5,6 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF( idoles )THEN
      IF( sgsmodel.ge.3 .and. sgsmodel.le.6 )THEN
      IF( terrain_flag )THEN
        if(myid.eq.0)then
        print *
        print *,'  idoles        = ',idoles
        print *,'  sgsmodel      = ',sgsmodel
        print *,'  terrain_flag  = ',terrain_flag
        print *
        print *,'  cannot use terrain with sgsmodel = 3,4,5,6 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      ENDIF
      !-----
      IF( horizturb.ge.1 )THEN
      IF( .not. idopbl )THEN
        if(myid.eq.0)then
        print *
        print *,'  horizturb = ',horizturb
        print *,'  idopbl    = ',idopbl
        print *
        print *,'  horizturb is only appropriate for idopbl '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF( idopbl )THEN
      IF( ipbl.le.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  idopbl   = ',idopbl
        print *,'  ipbl     = ',ipbl
        print *
        print *,'  idopbl requires ipbl >= 1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF( ipbl.ge.1 )THEN
      IF( .not. idopbl )THEN
        if(myid.eq.0)then
        print *
        print *,'  ipbl     = ',ipbl
        print *,'  idopbl   = ',idopbl
        print *,'  cm1setup = ',cm1setup
        print *
        print *,'  ipbl >= 1  requires cm1setup=2,4 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF(bcturbs.lt.1.or.bcturbs.gt.2)THEN
        if(myid.eq.0)then
        print *
        print *,'  bcturbs = ',bcturbs
        print *
        print *,'  bcturbs must be 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(imoist.eq.1 .and. ( isnd.eq.2 .or. isnd.eq.3 ) )THEN
        if(myid.eq.0)then
        print *
        print *,'  imoist = ',imoist
        print *,'  isnd   = ',isnd
        print *
        print *,'  For this value of isnd, imoist must be 0'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( bbc.lt.1 .or. bbc.gt.3 )THEN
        if(myid.eq.0)then
        print *
        print *,'  bbc = ',bbc
        print *
        print *,'  bbc must be 1, 2, or 3'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( tbc.lt.1 .or. tbc.gt.3 )THEN
        if(myid.eq.0)then
        print *
        print *,'  tbc = ',tbc
        print *
        print *,'  tbc must be 1,2,3'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(cm1setup.eq.3 .and. (bc_temp.le.0 .or. bc_temp.ge.3))THEN
        if(myid.eq.0)then
        print *
        print *,'  cm1setup = ',cm1setup
        print *,'  bc_temp  = ',bc_temp
        print *
        print *,'  for cm1setup=3, bc_temp must be either 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(ihail.lt.0.or.ihail.gt.1)THEN
        if(myid.eq.0)then
        print *
        print *,'  ihail   = ',ihail
        print *
        print *,'  ihail must be 0 or 1'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(imoist.eq.1.and.output_dbz.eq.1.and.ptype.ne.2.and.ptype.ne.3.and.ptype.ne.5.and.ptype.ne.7  &
          .and. (.not. ptype.ge.26))then
        if(myid.eq.0)then
        print *
        print *,'  ptype      = ',ptype
        print *,'  output_dbz = ',output_dbz
        print *
        print *,'  output_dbz is only available for ptype=2,3,5,7,26,27,28,50-55'
        print *
        endif
        IF(ptype.eq.4)THEN
          print *,'   stopping model .... '
          print *
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        ELSE
          output_dbz = 0
        ENDIF
      ENDIF
      !-----
      IF( imoist.eq.1 .and. eqtset.ge.2 .and. (ptype.eq.4) )THEN
        if(myid.eq.0)then
        print *
        print *,'  eqtset  = ',eqtset
        print *,'  ptype   = ',ptype
        print *
        print *,'  eqtset = 2 is not available for ptype = 4'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(imoist.eq.1 .and. efall.eq.1)THEN
      IF(ptype.ne.1.and.ptype.ne.2.and.ptype.ne.5.and.ptype.ne.6)THEN
        if(myid.eq.0)then
        print *
        print *,'  efall   = ',efall
        print *,'  ptype   = ',ptype
        print *
        print *,'  efall = 1 is only supported with ptype = 1,2,5,6'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF((imoist.eq.1).and.(ptype.eq.4).and.terrain_flag)THEN
        if(myid.eq.0)then
        print *
        print *,'  ptype   = ',ptype
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  ptype = 4 does not work with terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF((sfcmodel.eq.5).and.terrain_flag)THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel     = ',sfcmodel
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  sfcmodel=5 does not work with terrain (yet) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( idopbl .and. ipbl.eq.6 .and. sfcmodel.ne.7 )THEN
        if(myid.eq.0)then
        print *
        print *,'  ipbl         = ',ipbl
        print *,'  sfcmodel     = ',sfcmodel
        print *
        print *,'  ipbl=6 requires sfcmodel=7  (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(imoist.eq.1 .and. ptype.eq.5 .and. ndcnst.le.1.0e-6)THEN
        if(myid.eq.0)then
        print *
        print *,'  imoist  = ',imoist
        print *,'  ptype   = ',ptype
        print *,'  ndcnst  = ',ndcnst
        print *
        print *,'  ndcnst is too small.  Please enter an appropriate value '
        print *,'  in the param3 section of namelist.input '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(imoist.eq.1 .and. ptype.eq.3 .and. nt_c.le.1.0e-6)THEN
        if(myid.eq.0)then
        print *
        print *,'  imoist  = ',imoist
        print *,'  ptype   = ',ptype
        print *,'  nt_c    = ',nt_c   
        print *
        print *,'  nt_c is too small.  Please enter an appropriate value '
        print *,'  in the param3 section of namelist.input '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( terrain_flag .and. (irdamp.eq.2 .or. hrdamp.eq.2) )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  irdamp       = ',irdamp
        print *,'  hrdamp       = ',hrdamp
        print *
        print *,'  cannot use irdamp=2 or hrdamp=2 with terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(terrain_flag .and. (psolver.eq.4.or.psolver.eq.5) )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  psolver      = ',psolver
        print *
        print *,'  for psolver = 4,5 terrain_flag must be .false.'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( terrain_flag .and. ( ibalance.eq.2 .or. pdcomp ) )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  ibalance     = ',ibalance
        print *,'  pdcomp       = ',pdcomp
        print *
        print *,'  for ibalance.eq.2 and pdcomp, terrain_flag must be .false.'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( axisymm.eq.1 )THEN
      IF( psolver.eq.4 .or. psolver.eq.5 .or. ibalance.eq.2 .or. pdcomp )THEN
        if(myid.eq.0)then
        print *
        print *,'  axisymm      = ',axisymm
        print *
        print *,'  psolver      = ',psolver
        print *,'  ibalance     = ',ibalance
        print *,'  pdcomp       = ',pdcomp
        print *
        print *,'  axisymm cannot be used with psolver=4,5 or ibalance=2 or pdcomp=T '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF( stretch_x.ge.1 .or. stretch_y.ge.1 )THEN
      IF( psolver.eq.4 .or. psolver.eq.5 .or. ibalance.eq.2 .or. pdcomp )THEN
        if(myid.eq.0)then
        print *
        print *,'  stretch_x    = ',stretch_x
        print *,'  stretch_y    = ',stretch_y
        print *
        print *,'  psolver      = ',psolver
        print *,'  ibalance     = ',ibalance
        print *,'  pdcomp       = ',pdcomp
        print *
        print *,'  stretched horiz grid cannot be used with psolver=4,5 or ibalance=2 or pdcomp=T '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
      IF(iinit.eq.6)THEN
        if(myid.eq.0)then
        print *
        print *,'  iinit        = ',iinit
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (output_format.le.0) .or. (output_format.ge.3) )THEN
        if(myid.eq.0)then
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  only output_format = 1,2 are currently supported'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (axisymm.eq.1) .and. (iorigin.ne.1) )THEN
        if(myid.eq.0)then
        print *
        print *,'  iorigin = ',iorigin
        print *
        print *,'  axisymm=1 requires iorigin=1'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (axisymm.eq.1) .and. (imove.ne.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  imove = ',imove
        print *
        print *,'  axisymm=1 requires imove=0'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( terrain_flag .and. (imove.ne.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  imove = ',imove
        print *
        print *,'  imove must be 0 when using terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( imove.eq.1 .and. abs(umove).ge.1.0e-10 )THEN
      if( wbc.eq.3 .or. wbc.eq.4 .or. ebc.eq.3 .or. ebc.eq.4 )then
        if(myid.eq.0)then
        print *
        print *,'  imove = ',imove
        print *,'  umove = ',umove
        print *
        print *,'  wbc   = ',wbc
        print *,'  ebc   = ',ebc
        print *
        print *,'  Cannot move domain in east/west direction when using rigid-wall '
        print *,'  boundary conditions on east/west boundaries '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      ENDIF
      !-----
      IF( imove.eq.1 .and. abs(vmove).ge.1.0e-10 )THEN
      if( sbc.eq.3 .or. sbc.eq.4 .or. nbc.eq.3 .or. nbc.eq.4 )then
        if(myid.eq.0)then
        print *
        print *,'  imove = ',imove
        print *,'  vmove = ',vmove
        print *
        print *,'  sbc   = ',sbc
        print *,'  nbc   = ',nbc
        print *
        print *,'  Cannot move domain in south/north direction when using rigid-wall '
        print *,'  boundary conditions on south/north boundaries '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      ENDIF
      !-----
      IF( imove.eq.1 .and. initsfc.ne.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  imove   = ',imove
        print *,'  initsfc = ',initsfc
        print *
        print *,'  For imove=1, initsfc must be 1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (axisymm.eq.1) .and. terrain_flag )THEN
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  axisymm=1 cannot be used with terrain '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( icor.eq.0 ) fcor = 0.0
      !-----
      IF( (axisymm.eq.1) .and. (wbc.ne.3) )THEN
        if(myid.eq.0)then
        print *
        print *,'  wbc = ',wbc
        print *
        print *,'  axisymm=1 requires wbc=3 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (axisymm.eq.1) .and. ( (sbc.ne.1).or.(nbc.ne.1) ) )THEN
        if(myid.eq.0)then
        print *
        print *,'  sbc = ',sbc
        print *,'  nbc = ',nbc
        print *
        print *,'  axisymm=1 requires sbc=nbc=1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (axisymm.eq.1).and.(ny.gt.1) )THEN
        if(myid.eq.0)then
        print *
        print *,'  ny = ',ny
        print *
        print *,'  axisymm=1 requires ny=1'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( axisymm.eq.1 .and. idoles )THEN
        if(myid.eq.0)then
        print *
        print *,'  cm1setup = ',cm1setup
        print *
        print *,'  axisymm=1 cannot be used with cm1setup=1,4 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( axisymm.eq.1.and.(psolver.eq.1.or.psolver.eq.4.or.psolver.eq.5) )THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver    = ',psolver
        print *
        print *,'  axisymm=1 is only available with psolver=2,3,6'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( axisymm.eq.1 .and. idiff.eq.1 .and. difforder.eq.6 )THEN
        idiff = 0
        difforder = 0
      ENDIF
      !-----
      IF( (bbc.eq.3) .and. (sfcmodel.le.0) )THEN
        if(myid.eq.0)then
        print *
        print *,'  bbc      = ',bbc
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  bbc=3 requires a setting for sfcmodel '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( isfcflx.ne.0 )THEN
      IF( sfcmodel.lt.1 .or. sfcmodel.gt.7 )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *
        print *,'  sfcmodel must be 1,2,3,4,5,6,7 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      ENDIF
      !-----
!      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4.or.oceanmodel.eq.2).and.imove.ne.0 )THEN
!        if(myid.eq.0)then
!        print *
!        print *,'  sfcmodel    = ',sfcmodel
!        print *,'  oceanmodel  = ',sfcmodel
!        print *,'  imove       = ',imove
!        print *
!        print *,'  domain translation is now allowed with sfcmodel = 2,3,4  and/or oceanmodel = 2 '
!        print *
!        print *,'   stopping model .... '
!        print *
!        endif
!#ifdef MPI
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
!        call stopcm1
!      ENDIF
      !-----
      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4.or.sfcmodel.eq.6.or.sfcmodel.eq.7).and.(season.le.0.or.season.ge.3) )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *,'  season   = ',season
        print *
        print *,'  season must have a value of 1 or 2 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( pertflx.eq.1 .and. sfcmodel.ge.2 )THEN
        if(myid.eq.0)then
        print *
        print *,'  pertflx  = ',pertflx
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  pertflx can only be used with sfcmodel = 1  '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( pertflx.eq.1 .and. imove.ge.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  pertflx  = ',pertflx
        print *,'  imove    = ',imove
        print *
        print *,'  pertflx can only be used with imove=0 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( sfcmodel.eq.1 .and. oceanmodel.ne.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  oceanmodel = ',oceanmodel
        print *
        print *,'  sfcmodel = 1 requires oceanmodel = 1 '
        print *,'  (oceanmodel = 2 requires sfcmodel = 2 ) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( ( testcase.eq.4 .or. testcase.eq.5 ) .and. radopt.ge.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  testcase = ',radopt
        print *,'  radopt    = ',radopt
        print *
        print *,'  for testcase = 4,5  radopt must be 0 '
        print *,' (simple radiative tendencies are handled in config_simple_phys subroutine)'
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( cuparam.lt.0 .or. cuparam.gt.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  cuparam  = ',cuparam
        print *
        print *,'  cuparam must be 0 or 1 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( radopt.lt.0 .or. radopt.gt.2 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *
        print *,'  radopt must be 0,1,2 '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (radopt.eq.1.or.radopt.eq.2) .and. imoist.eq.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  imoist   = ',imoist
        print *
        print *,'  radopt=1,2 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
!      IF( ipbl.eq.1 .and. imoist.eq.0 )THEN
!        if(myid.eq.0)then
!        print *
!        print *,'  ipbl     = ',ipbl
!        print *,'  imoist   = ',imoist
!        print *
!        print *,'  ipbl=1 requires imoist=1 (for now) '
!        print *
!        print *,'   stopping model .... '
!        print *
!        endif
!#ifdef MPI
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
!        call stopcm1
!      ENDIF
      !-----
      IF( (radopt.eq.1.or.radopt.eq.2) .and. rterm.eq.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  rterm    = ',rterm
        print *
        print *,'  cannot use radopt and rterm at the same time '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( radopt.eq.2 .and. ( ptype.ne.3 .and. ptype.ne.5 .and. ptype.ne.7 .and. ptype.ne.50 .and.  &
                              ptype.ne.51 .and. ptype.ne.52 .and. ptype.ne.53 .and.&
                              ptype.ne.55 .and. ptype.ne.26 .and. ptype.ne.27)  )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  ptype    = ',ptype
        print *
        print *,'  radopt=2 requires ptype=3,5,7,26,27,50,51,52,53,55 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (radopt.eq.1.or.radopt.eq.2) .and. (ptype.eq.1.or.ptype.eq.6)  )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  ptype    = ',ptype
        print *
        print *,'  radopt=1,2 requires an ice microphysics scheme (for now) '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF( (radopt.eq.1.or.radopt.eq.2) .and. sfcmodel.eq.0 )THEN
        if(myid.eq.0)then
        print *
        print *,'  radopt   = ',radopt
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  radopt=1,2 requires a surface model '
        print *
        print *,'   stopping model .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
!      IF( (sfcmodel.eq.2.or.sfcmodel.eq.3.or.sfcmodel.eq.4) .and. imoist.eq.0 )THEN
!        if(myid.eq.0)then
!        print *
!        print *,'  sfcmodel = ',sfcmodel
!        print *,'  imoist   = ',imoist
!        print *
!        print *,'  sfcmodel=2,3,4 requires imoist=1 (for now) '
!        print *
!        print *,'   stopping model .... '
!        print *
!        endif
!#ifdef MPI
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
!        call stopcm1
!      ENDIF
      !-----
#ifndef NETCDF
      !-----
      IF(output_format.eq.2)THEN
        if(myid.eq.0)then
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  You have requested netcdf output, but you have not'
        print *,'  compiled the code with netcdf capability.  Modify the'
        print *,'  Makefile, clean, and recompile'
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
      IF(restart_format.eq.2)THEN
        if(myid.eq.0)then
        print *
        print *,'  restart_format = ',restart_format
        print *
        print *,'  You have requested netcdf restart files, but you have not'
        print *,'  compiled the code with netcdf capability.  Modify the'
        print *,'  Makefile, clean, and recompile'
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF
      !-----
#endif
!#ifdef NETCDF
!      !-----
!      IF( terrain_flag .and. output_interp.ne.0 .and. output_format.eq.2 )THEN
!        if(myid.eq.0)then
!        print *
!        print *,'  output_interp = ',output_interp
!        print *
!        print *,'  output_interp=1 is not currently available for netcdf output'
!        print *
!        endif
!#ifdef MPI
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
!        call stopcm1
!      ENDIF
!      !-----
!#endif
#ifdef MPI
      !-----
      IF(psolver.eq.4.or.psolver.eq.5)THEN
        if(myid.eq.0)then
        print *
        print *,'  psolver = ',psolver
        print *
        print *,'  psolver = 4 and 5 are not supported in MPI mode'
        print *
        print *,'   stopping model .... '
        print *
        endif
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call stopcm1
      ENDIF
      !-----
      IF( output_pdcomp.ge.1 )THEN
        if(myid.eq.0)then
        print *
        print *,'  output_pdcomp = ',output_pdcomp
        print *
        print *,'  pressure decomposition output is not supported with MPI '
        print *
        print *,'   stopping model .... '
        print *
        endif
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call stopcm1
      ENDIF
      !-----
      IF(axisymm.eq.1)THEN
        if(myid.eq.0)then
        print *
        print *,'  axisymm = ',axisymm
        print *
        print *,'  axisymm is not supported in MPI mode'
        print *
        print *,'   stopping model .... '
        print *
        endif
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call stopcm1
      ENDIF
      !-----
      IF(ny.lt.3)THEN
        if(myid.eq.0)then
        print *
        print *,'  ny = ',ny
        print *
        print *,'  ny must be .ge. 3 for  MPI runs'
        print *
        print *,'   stopping model .... '
        print *
        endif
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call stopcm1
      ENDIF
      !-----
#endif


!--------------------------------------------------------------
!  Check domain size (MPI only)
#ifdef MPI

!  cm1r20:  nodex and nodey are determined automatically
!        if((nodex*nodey).ne.numprocs)then
!          if(myid.eq.0)then
!          print *
!          print *,'  WARNING!!! '
!          print *,'  nodes does not equal numprocs!'
!          print *,'  nodex,nodey,nodes=',nodex,nodey,nodex*nodey
!          print *,'  numprocs=',numprocs
!          print *
!          endif
!          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!          call stopcm1
!        endif

!  cm1r20:  no longer required
!        if(mod(nx,nodex).ne.0)then
!          if(myid.eq.0)then
!          print *
!          print *,'  nx does not divide exactly by nodex! '
!          print *,'  nx,nodex,mod(nx,nodex)=',nx,nodex,mod(nx,nodex)
!          print *
!          endif
!          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!          call stopcm1
!        endif

!  cm1r20:  no longer required
!        if(mod(ny,nodey).ne.0)then
!          if(myid.eq.0)then
!          print *
!          print *,'  ny does not divide exactly by nodey! '
!          print *,'  ny,nodey,mod(ny,nodey)=',ny,nodey,mod(ny,nodey)
!          print *
!          endif
!          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!          call stopcm1
!        endif

        IF( (ni.lt.3).or.(nj.lt.3) )THEN
          if(myid.eq.0)then
          print *
          print *,'  myid = ',myid
          print *,'  ni = ',ni
          print *,'  nj = ',nj
          print *,'  both ni and nj must be >= 3 '
          print *
          endif
          ! bug fix, cm1r20.3: do not wait for all processors
!!!          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call stopcm1
        ENDIF

        ppnode = min( ppnode , nodex*nodey )

        if(mod(nodex*nodey,ppnode).ne.0)then
          if(myid.eq.0)then
          print *
          print *,'  nodex*nodey does not divide exactly by ppnode! '
          print *,'  nodex*nodey,ppnode=',nodex*nodey,ppnode
          print *
          endif
          ! bug fix, cm1r20.3: do not wait for all processors
!!!          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call stopcm1
        endif

!!!      if(myid.eq.0)then
!!!
!!!        print *
!!!        print *,'  Everything is cool!'
!!!        print *,'  ni,nj=',ni,nj
!!!        print *
!!!
!!!      endif

#endif
!--------------------------------------------------------------
!  Check that lateral bc combinations make sense:

      if(ebc.eq.1 .and. wbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      if(wbc.eq.1 .and. ebc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      if(nbc.eq.1 .and. sbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      if(sbc.eq.1 .and. nbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

!  end fatal checks  (ie, model stops)
!--------------------------------------------------------------
!  Some basic checks:

      ! for passive tracers:
      iptra    = max(0,min(1,iptra))
      if(iptra.eq.1)then
        npt      = max(1,npt)
      else
        npt      = 1
      endif

      ! for parcels:
      nparcels = max(1,nparcels)

!-----

      if(stretch_z.lt.1) ztop = dz*float(nk)
      IF ( stretch_z == 2 ) dz = ztop/float(nk) ! nk is the number of scalar levels

      IF( advwenos.lt.0 .or. advwenos.gt.2 )THEN
        print *
        print *,'  advwenos = ',advwenos
        print *
        print *,'  unrecognized value for advwenos '
        print *
        print *,'   stopping model .... '
        print *
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

      !-----

      IF( advwenov.lt.0 .or. advwenov.gt.2 )THEN
        print *
        print *,'  advwenov = ',advwenov
        print *
        print *,'  unrecognized value for advwenov '
        print *
        print *,'   stopping model .... '
        print *
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Domain dimensions: '
      if(dowr) write(outfile,*) 'nx    =',nx
      if(dowr) write(outfile,*) 'ny    =',ny
      if(dowr) write(outfile,*) 'nz    =',nz

#ifdef MPI
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  MPI info: '
      if(dowr) write(outfile,*) 'nodex    =',nodex
      if(dowr) write(outfile,*) 'nodey    =',nodey
      if(dowr) write(outfile,*) 'numprocs =',numprocs
      if(dowr) write(outfile,*) 'ppnode   =',ppnode
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  sub-domain dimensions: '

      doit = .false.
#ifdef MPI
      doit = .false.
#endif

    IF( .not. doit )THEN
      if(dowr) write(outfile,*) 'ni    =',ni
      if(dowr) write(outfile,*) 'nj    =',nj
      if(dowr) write(outfile,*) 'nk    =',nk
#ifdef MPI
    ELSE
      if(myid.eq.0)then
        do n=0,(numprocs-1)
          if(n.eq.0)then
            ntmp1 = myid
            ntmp2 = ni
            ntmp3 = nj
            ntmp4 = nk
          else
            call MPI_IRECV(ntmp1,1,mpi_integer,n,9991,MPI_COMM_WORLD,reqs1,ierr)
            call MPI_IRECV(ntmp2,1,mpi_integer,n,9992,MPI_COMM_WORLD,reqs2,ierr)
            call MPI_IRECV(ntmp3,1,mpi_integer,n,9993,MPI_COMM_WORLD,reqs3,ierr)
            call MPI_IRECV(ntmp4,1,mpi_integer,n,9994,MPI_COMM_WORLD,reqs4,ierr)
            call MPI_WAIT(reqs1,status,ierr)
            call MPI_WAIT(reqs2,status,ierr)
            call MPI_WAIT(reqs3,status,ierr)
            call MPI_WAIT(reqs4,status,ierr)
          endif
          print *,'myid,ni,nj,nk = ',ntmp1,ntmp2,ntmp3,ntmp4
        enddo
      else
        call MPI_ISEND(myid,1,mpi_integer,0,9991,MPI_COMM_WORLD,reqs1,ierr)
        call MPI_ISEND( ni ,1,mpi_integer,0,9992,MPI_COMM_WORLD,reqs2,ierr)
        call MPI_ISEND( nj ,1,mpi_integer,0,9993,MPI_COMM_WORLD,reqs3,ierr)
        call MPI_ISEND( nk ,1,mpi_integer,0,9994,MPI_COMM_WORLD,reqs4,ierr)
        call MPI_WAIT(reqs1,status,ierr)
        call MPI_WAIT(reqs2,status,ierr)
        call MPI_WAIT(reqs3,status,ierr)
        call MPI_WAIT(reqs4,status,ierr)
      endif
#endif
    ENDIF
      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      ! default text:
      param_mp     = 'no microphysics (dry)'
      param_sgs    = 'no LES subgrid turbulence model'
      param_pbl    = 'no PBL model'
      param_sfclay = 'no surface layer param.'
      param_rad    = 'no lw/sw radiation '

      if( imoist.eq.1 )then
        ! microphysics:
        if( ptype.eq.0  ) param_mp = 'water vapor only (no microphysics)'
        if( ptype.eq.1  ) param_mp = 'Kessler'
        if( ptype.eq.2  ) param_mp = 'NASA-GSFC version of LFO'
        if( ptype.eq.3  ) param_mp = 'Thompson'
        if( ptype.eq.4  ) param_mp = 'Gilmore/Straka/Rasmussen version of LFO scheme'
        if( ptype.eq.5  ) param_mp = 'Morrison double moment'
        if( ptype.eq.6  ) param_mp = 'Rotunno-Emanuel water-only'
        if( ptype.eq.7  ) param_mp = 'WSM6'
        if( ptype.eq.26 ) param_mp = 'NSSL 2-moment scheme (graupel-only, no hail)'
        if( ptype.eq.27 ) param_mp = 'NSSL 2-moment scheme (graupel and hail)'
        if( ptype.eq.28 ) param_mp = 'NSSL single-moment scheme (graupel-only)'
        if( ptype.eq.50 ) param_mp = 'Predicted Particle Property (P3)'
        if( ptype.eq.51 ) param_mp = 'Predicted Particle Property (P3)'
        if( ptype.eq.52 ) param_mp = 'Predicted Particle Property (P3)'
        if( ptype.eq.53 ) param_mp = 'Predicted Particle Property (P3)'
        if( ptype.eq.55 ) param_mp = 'Jensen_ISHMAEL'
      endif
      if( idoles )then
        ! LES subgrid:
        if( sgsmodel.eq.1 ) param_sgs = 'Deardorff TKE scheme'
        if( sgsmodel.eq.2 ) param_sgs = 'Smagorinsky scheme'
        if( sgsmodel.eq.3 ) param_sgs = 'Deardorff TKE scheme + S94 two-part model'
        if( sgsmodel.eq.4 ) param_sgs = 'Deardorff TKE scheme + B20 two-part model'
        if( sgsmodel.eq.5 ) param_sgs = 'NBA (Deardorff-type TKE version)'
        if( sgsmodel.eq.6 ) param_sgs = 'NBA (Smagorinksy-type version)'
      endif
      if( idopbl )then
        ! PBL:
        if( ipbl.eq.1 ) param_pbl = 'YSU (Yonsei University)'
        if( ipbl.eq.2 ) param_pbl = 'CM1 simple scheme (Bryan-Rotunno, Louis-type scheme)'
        if( ipbl.eq.3 ) param_pbl = 'GFS-EDMF (as configured in HWRF)'
        if( ipbl.eq.4 ) param_pbl = 'MYNN (Mellor-Yamada Nakanishi-Niino) level 2.5'
        if( ipbl.eq.5 ) param_pbl = 'MYNN (Mellor-Yamada Nakanishi-Niino) level 3'
        if( ipbl.eq.6 ) param_pbl = 'MYJ (Mellor-Yamada-Janjic)'
      endif
      if( bbc.eq.3 )then
        ! surface layer:
        if( sfcmodel.eq.1 ) param_sfclay = 'CM1 simple scheme (neutral stability)'
        if( sfcmodel.eq.2 ) param_sfclay = 'old surface-layer scheme from WRF/MM5'
        if( sfcmodel.eq.3 ) param_sfclay = 'revised surface-layer scheme from WRF'
        if( sfcmodel.eq.4 ) param_sfclay = 'GFDL surface layer'
        if( sfcmodel.eq.5 ) param_sfclay = 'simple Monin-Obukhov Similarity Theory (MOST)'
        if( sfcmodel.eq.6 ) param_sfclay = 'MYNN (Mellor-Yamada Nakanishi-Niino) surface layer'
        if( sfcmodel.eq.7 ) param_sfclay = 'MYJ (Mellor-Yamada-Janjic) surface layer'
      endif
      if( radopt.gt.0 )then
        ! radiation:
        if( radopt.eq.1 ) param_rad = 'NASA-GSFC'
        if( radopt.eq.2 ) param_rad = 'RRTMG'
      endif

      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171)
      if(cm1setup.eq.0 )then
        if(dowr) write(outfile,172) '  cm1setup = ',cm1setup,': no subgrid turbulence, no explicit diffusion (Euler equations)'
      endif
      if(cm1setup.eq.1 )then
        if(dowr) write(outfile,172) '  cm1setup = ',cm1setup,': Large-Eddy Simulation (LES) '
      endif
      if(cm1setup.eq.2 )then
        if(dowr) write(outfile,172) '  cm1setup = ',cm1setup,': mesoscale modeling with PBL parameterization '
      endif
      if(cm1setup.eq.3 )then
        if(dowr) write(outfile,172) '  cm1setup = ',cm1setup,': Direct Numerical Simulation (DNS) '
      endif
      if(cm1setup.eq.4 )then
        if(dowr) write(outfile,172) '  cm1setup = ',cm1setup,': Large-Eddy Simulation (LES) embedded within Mesoscale Model (with PBL parameterization) '
      endif
      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171)
      if( psolver.eq.1 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': compressible, explicit p-grad, without time-splitting'
      endif
      if( psolver.eq.2 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': compressible, explicit p-grad, with time-splitting'
      endif
      if( psolver.eq.3 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': compressible, vertically implicit p-grad, with time-splitting'
      endif
      if( psolver.eq.4 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': deep anelastic solver'
      endif
      if( psolver.eq.5 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': incompressible solver'
      endif
      if( psolver.eq.6 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': compressible-Boussinesq, explicit p-grad, with time-splitting'
      endif
      if( psolver.eq.7 )then
        if(dowr) write(outfile,172) '  psolver = ',psolver,': modified compressible, explicit p-grad, with time-splitting'
      endif
      if(dowr) write(outfile,171)

      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171) '  Using these physical parameterizations: '
      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171) '    microphysics             :  ',param_mp
      if(dowr) write(outfile,171) '    LES subgrid turbulence   :  ',param_sgs
      if(dowr) write(outfile,171) '    planetary boundary layer :  ',param_pbl
      if(dowr) write(outfile,171) '    surface layer            :  ',param_sfclay
      if(dowr) write(outfile,171) '    lw/sw radiation          :  ',param_rad
    if( cuparam.ge.1 )then
      if(dowr) write(outfile,171) '    convection param.        :  ','new Tiedtke'
    endif
      if(dowr) write(outfile,171)
      if(dowr) write(outfile,171)

171   format(a,a)
172   format(a,1x,i1,1x,a)

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'dx        =',dx
      if(dowr) write(outfile,*) 'dy        =',dy
      if(dowr) write(outfile,*) 'dz        =',dz
      if(dowr) write(outfile,*) 'dtl       =',dtl
      if(dowr) write(outfile,*) 'timax     =',timax
      if(dowr) write(outfile,*) 'run_time  =',run_time
      if(dowr) write(outfile,*) 'tapfrq    =',tapfrq
      if(dowr) write(outfile,*) 'rstfrq    =',rstfrq
      if(dowr) write(outfile,*) 'statfrq   =',statfrq
      if(dowr) write(outfile,*) 'prclfrq   =',prclfrq
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'cm1setup  =',cm1setup
      if(dowr) write(outfile,*) 'testcase  =',testcase
      if(dowr) write(outfile,*) 'adapt_dt  =',adapt_dt
      if(dowr) write(outfile,*) 'irst      =',irst
      if(dowr) write(outfile,*) 'rstnum    =',rstnum
      if(dowr) write(outfile,*) 'iconly    =',iconly
      if(dowr) write(outfile,*) 'hadvordrs =',hadvordrs
      if(dowr) write(outfile,*) 'vadvordrs =',vadvordrs
      if(dowr) write(outfile,*) 'hadvordrv =',hadvordrv
      if(dowr) write(outfile,*) 'vadvordrv =',vadvordrv
      if(dowr) write(outfile,*) 'advwenos  =',advwenos
      if(dowr) write(outfile,*) 'advwenov  =',advwenov
      if(dowr) write(outfile,*) 'weno_order=',weno_order
      if(dowr) write(outfile,*) 'apmasscon =',apmasscon
      if(dowr) write(outfile,*) 'idiff     =',idiff
      if(dowr) write(outfile,*) 'mdiff     =',mdiff
      if(dowr) write(outfile,*) 'difforder =',difforder
      if(dowr) write(outfile,*) 'imoist    =',imoist
      if(dowr) write(outfile,*) 'ipbl      =',ipbl
      if(dowr) write(outfile,*) 'sgsmodel  =',sgsmodel
      if(dowr) write(outfile,*) 'tconfig   =',tconfig
      if(dowr) write(outfile,*) 'bcturbs   =',bcturbs
      if(dowr) write(outfile,*) 'horizturb =',horizturb
      if(dowr) write(outfile,*) 'doimpl    =',doimpl
      if(dowr) write(outfile,*) 'irdamp    =',irdamp
      if(dowr) write(outfile,*) 'hrdamp    =',hrdamp
      if(dowr) write(outfile,*) 'psolver   =',psolver
      if(dowr) write(outfile,*) 'ptype     =',ptype
      if(dowr) write(outfile,*) 'ihail     =',ihail
      if(dowr) write(outfile,*) 'iautoc    =',iautoc
      if(dowr) write(outfile,*) 'cuparam   =',cuparam
      if(dowr) write(outfile,*) 'icor      =',icor
      if(dowr) write(outfile,*) 'betaplane =',betaplane
      if(dowr) write(outfile,*) 'lspgrad   =',lspgrad
      if(dowr) write(outfile,*) 'eqtset    =',eqtset
      if(dowr) write(outfile,*) 'idiss     =',idiss
      if(dowr) write(outfile,*) 'efall     =',efall
      if(dowr) write(outfile,*) 'rterm     =',rterm
      if(dowr) write(outfile,*) 'wbc       =',wbc
      if(dowr) write(outfile,*) 'ebc       =',ebc
      if(dowr) write(outfile,*) 'sbc       =',sbc
      if(dowr) write(outfile,*) 'nbc       =',nbc
      if(dowr) write(outfile,*) 'bbc       =',bbc
      if(dowr) write(outfile,*) 'tbc       =',tbc
      if(dowr) write(outfile,*) 'irbc      =',irbc
      if(dowr) write(outfile,*) 'roflux    =',roflux
      if(dowr) write(outfile,*) 'nudgeobc  =',nudgeobc
      if(dowr) write(outfile,*) 'isnd      =',isnd
      if(dowr) write(outfile,*) 'iwnd      =',iwnd
      if(dowr) write(outfile,*) 'itern     =',itern
      if(dowr) write(outfile,*) 'iinit     =',iinit
      if(dowr) write(outfile,*) 'irandp    =',irandp
      if(dowr) write(outfile,*) 'ibalance  =',ibalance
      if(dowr) write(outfile,*) 'iorigin   =',iorigin
      if(dowr) write(outfile,*) 'axisymm   =',axisymm
      if(dowr) write(outfile,*) 'imove     =',imove
      if(dowr) write(outfile,*) 'iptra     =',iptra
      if(dowr) write(outfile,*) 'npt       =',npt
      if(dowr) write(outfile,*) 'pdtra     =',pdtra
      if(dowr) write(outfile,*) 'iprcl     =',iprcl
      if(dowr) write(outfile,*) 'nparcels  =',nparcels
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'kdiff2    =',kdiff2
      if(dowr) write(outfile,*) 'kdiff6    =',kdiff6
      if(dowr) write(outfile,*) 'fcor      =',fcor
      if(dowr) write(outfile,*) 'kdiv      =',kdiv
      if(dowr) write(outfile,*) 'alph      =',alph
      if(dowr) write(outfile,*) 'rdalpha   =',rdalpha
      if(dowr) write(outfile,*) 'zd        =',zd
      if(dowr) write(outfile,*) 'xhd       =',xhd
      if(dowr) write(outfile,*) 'alphobc   =',alphobc
      if(dowr) write(outfile,*) 'umove     =',umove
      if(dowr) write(outfile,*) 'vmove     =',vmove
      if(dowr) write(outfile,*) 'v_t       =',v_t
      if(dowr) write(outfile,*) 'l_h       =',l_h
      if(dowr) write(outfile,*) 'lhref1    =',lhref1
      if(dowr) write(outfile,*) 'lhref2    =',lhref2
      if(dowr) write(outfile,*) 'l_inf     =',l_inf
      if(dowr) write(outfile,*) 'ndcnst    =',ndcnst
      if(dowr) write(outfile,*) 'nt_c      =',nt_c
      if(dowr) write(outfile,*) 'csound    =',csound
      if(dowr) write(outfile,*) 'cstar     =',cstar
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'radopt    =',radopt
      if(dowr) write(outfile,*) 'dtrad     =',dtrad
      if(dowr) write(outfile,*) 'ctrlat    =',ctrlat
      if(dowr) write(outfile,*) 'ctrlon    =',ctrlon
      if(dowr) write(outfile,*) 'year      =',year
      if(dowr) write(outfile,*) 'month     =',month
      if(dowr) write(outfile,*) 'day       =',day
      if(dowr) write(outfile,*) 'hour      =',hour
      if(dowr) write(outfile,*) 'minute    =',minute
      if(dowr) write(outfile,*) 'second    =',second
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'isfcflx   =',isfcflx
      if(dowr) write(outfile,*) 'sfcmodel  =',sfcmodel
      if(dowr) write(outfile,*) 'oceanmodel=',oceanmodel
      if(dowr) write(outfile,*) 'initsfc   =',initsfc
      if(dowr) write(outfile,*) 'tsk0      =',tsk0
      if(dowr) write(outfile,*) 'tmn0      =',tmn0
      if(dowr) write(outfile,*) 'xland0    =',xland0
      if(dowr) write(outfile,*) 'lu0       =',lu0
      if(dowr) write(outfile,*) 'season    =',season
      if(dowr) write(outfile,*) 'cecd      =',cecd
      if(dowr) write(outfile,*) 'pertflx   =',pertflx
      if(dowr) write(outfile,*) 'cnstce    =',cnstce
      if(dowr) write(outfile,*) 'cnstcd    =',cnstcd
      if(dowr) write(outfile,*) 'isftcflx  =',isftcflx
      if(dowr) write(outfile,*) 'iz0tlnd   =',iz0tlnd
      if(dowr) write(outfile,*) 'oml_hml0  =',oml_hml0
      if(dowr) write(outfile,*) 'oml_gamma =',oml_gamma
      if(dowr) write(outfile,*) 'set_flx   =',set_flx
      if(dowr) write(outfile,*) 'cnst_shflx=',cnst_shflx
      if(dowr) write(outfile,*) 'cnst_lhflx=',cnst_lhflx
      if(dowr) write(outfile,*) 'set_znt   =',set_znt
      if(dowr) write(outfile,*) 'cnst_znt  =',cnst_znt
      if(dowr) write(outfile,*) 'set_ust   =',set_ust
      if(dowr) write(outfile,*) 'cnst_ust  =',cnst_ust
      if(dowr) write(outfile,*) 'ramp_sgs  =',ramp_sgs
      if(dowr) write(outfile,*) 'ramp_time =',ramp_time
      if(dowr) write(outfile,*) 't2p_avg   =',t2p_avg
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_x =',stretch_x
      if(dowr) write(outfile,*) 'dx_inner  =',dx_inner
      if(dowr) write(outfile,*) 'dx_outer  =',dx_outer
      if(dowr) write(outfile,*) 'nos_x_len =',nos_x_len
      if(dowr) write(outfile,*) 'tot_x_len =',tot_x_len
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_y =',stretch_y
      if(dowr) write(outfile,*) 'dy_inner  =',dy_inner
      if(dowr) write(outfile,*) 'dy_outer  =',dy_outer
      if(dowr) write(outfile,*) 'nos_y_len =',nos_y_len
      if(dowr) write(outfile,*) 'tot_y_len =',tot_y_len
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_z =',stretch_z
      if(dowr) write(outfile,*) 'ztop      =',ztop
      if(dowr) write(outfile,*) 'str_bot   =',str_bot
      if(dowr) write(outfile,*) 'str_top   =',str_top
      if(dowr) write(outfile,*) 'dz_bot    =',dz_bot
      if(dowr) write(outfile,*) 'dz_top    =',dz_top
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'bc_temp   =',bc_temp
      if(dowr) write(outfile,*) 'ptc_top   =',ptc_top
      if(dowr) write(outfile,*) 'ptc_bot   =',ptc_bot
      if(dowr) write(outfile,*) 'viscosity =',viscosity
      if(dowr) write(outfile,*) 'pr_num    =',pr_num
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'var1      =',var1
      if(dowr) write(outfile,*) 'var2      =',var2
      if(dowr) write(outfile,*) 'var3      =',var3
      if(dowr) write(outfile,*) 'var4      =',var4
      if(dowr) write(outfile,*) 'var5      =',var5
      if(dowr) write(outfile,*) 'var6      =',var6
      if(dowr) write(outfile,*) 'var7      =',var7
      if(dowr) write(outfile,*) 'var8      =',var8
      if(dowr) write(outfile,*) 'var9      =',var9
      if(dowr) write(outfile,*) 'var10     =',var10
      if(dowr) write(outfile,*) 'var11     =',var11
      if(dowr) write(outfile,*) 'var12     =',var12
      if(dowr) write(outfile,*) 'var13     =',var13
      if(dowr) write(outfile,*) 'var14     =',var14
      if(dowr) write(outfile,*) 'var15     =',var15
      if(dowr) write(outfile,*) 'var16     =',var16
      if(dowr) write(outfile,*) 'var17     =',var17
      if(dowr) write(outfile,*) 'var18     =',var18
      if(dowr) write(outfile,*) 'var19     =',var19
      if(dowr) write(outfile,*) 'var20     =',var20
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'hurr_vg     = ',hurr_vg
      if(dowr) write(outfile,*) 'hurr_rad    = ',hurr_rad
      if(dowr) write(outfile,*) 'hurr_vgpl   = ',hurr_vgpl
      if(dowr) write(outfile,*) 'hurr_rotate = ',hurr_rotate
      if(dowr) write(outfile,*)


      IF ( ptype >= 26 .and. dowr ) THEN
        write(outfile,NML=nssl2mom_params)
!        write(outfile,*) 'alphah    =',alphah
!        write(outfile,*) 'alphahl   =',alphahl
!        write(outfile,*) 'dfrz      =',dfrz
!        write(outfile,*) 'hldnmn    =',hldnmn
!        write(outfile,*) 'imurain   =',imurain
!        write(outfile,*) 'ccn       =',ccn
!        write(outfile,*) 'icdx      =',icdx
!        write(outfile,*) 'icdxhl    =',icdxhl
!        write(outfile,*) 'iferwisventr =',iferwisventr
!        write(outfile,*) 'iehw      =',iehw
!        write(outfile,*) 'iehlw     =',iehlw
!        write(outfile,*) 'ehw0      =',ehw0
!        write(outfile,*) 'ehlw0     =',ehlw0
!        write(outfile,*) 'dmrauto   =',dmrauto
!        write(outfile,*) 'ioldlimiter=',ioldlimiter
      ENDIF

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  sp,dp,qp = ',sp,dp,qp

      pi_sp = 4.0_sp*atan(1.0_sp)
      pi_dp = 4.0_dp*atan(1.0_dp)
      pi_qp = 4.0_qp*atan(1.0_qp)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  pi_sp = ',pi_sp
      if(dowr) write(outfile,*) '  pi_dp = ',pi_dp
      if(dowr) write(outfile,*) '  pi_qp = ',pi_qp
      if(dowr) write(outfile,*)

      pi_phys = 4.0_kind_phys*atan(1.0_kind_phys)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  kind_phys = ',kind_phys
      if(dowr) write(outfile,*) '  pi_phys   = ',pi_phys
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  some parameters related to subgrid turbulence:

      dohturb = .false.
      dovturb = .false.

      if( idoles .or. cm1setup.eq.3 )then
        ! LES or DNS
        dohturb = .true.
        dovturb = .true.
      endif

      if( cm1setup.eq.2 .and. horizturb.eq.1 )then
        ! horizontal Smagorinsky
        dohturb = .true.
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  dohturb,dovturb = ',dohturb,dovturb

    !----------
    !  use average surface conditions to determine ust,mol,etc:
    !
    !  For now, use this for a few well-established test cases. 
    !  In future: make this a namelist variable

    if( sfcmodel.eq.1 .or. sfcmodel.eq.5 )then
      if( testcase.eq.1 .or. testcase.eq.2 .or. testcase.eq.9 .or. testcase.eq.11 .or. testcase.eq.14 )then
        use_avg_sfc = .true.
      else
        use_avg_sfc = .false.
      endif
    endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  use_avg_sfc = ',use_avg_sfc


!--------------------------------------------------------------
!  Configuration for simulations with moisture
!

      !--- begin: define defaults (please do not change) ---------
      iice     = 0
      idm      = 0
      idmplus  = 0
      numq     = 1
      nqv      = 1
      nql1     = -2
      nql2     = -2
      nqs1     = -2
      nqs2     = -2
      nnc1     = -2
      nnc2     = -2
      nzl1     = -2
      nzl2     = -2
      nvl1     = -2
      nvl2     = -2
      nbudget  = 10
      budrain  = 1
      cloudvar = .false.
      rhovar   = .false.
      ibtmp=1
      ietmp=1
      jbtmp=1
      jetmp=1
      kbtmp=1
      ketmp=1
      ntmp2d=1
      ntmp1d=1
      !--- end: define defaults ----------------------------------

      IF(imoist.eq.1)THEN

!-----------------------------------------------------------------------
!-------   BEGIN:  modify stuff below here -----------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        IF(ptype.eq.0)THEN        ! water vapor only

          numq = 1
          nqv = 1
          cloudvar(1) = .false.
          qname(1) = 'qv '
          qunit(1) = 'kg/kg'
          qmag( 1) = 1.0e-2

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.1)THEN        ! Kessler scheme

          numq = 3    ! there are 3 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  ... using Kessler microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
        ELSEIF((ptype.eq.2).or.(ptype.eq.4))THEN    ! Goddard-LFO or 
                                                    ! GSR-LFO scheme

          iice = 1    ! this means that ptype=2,4 are ice schemes

          numq = 6    ! there are 6 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'
          qunit(4) = 'kg/kg'
          qunit(5) = 'kg/kg'
          qunit(6) = 'kg/kg'

          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2
          qmag( 4) = 1.0e-2
          qmag( 5) = 1.0e-2
          qmag( 6) = 1.0e-2

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the Goddard or GSR LFO scheme -----

          if(ptype.eq.2)THEN

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) 'Calling CONSAT'
            if(dowr) write(outfile,*)

            call consat
            call consat2(dt)

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) ' ----------------------------------------------- '
            if(dowr) write(outfile,*) '  ... using Goddard LFO microphysics scheme ... '
            if(dowr) write(outfile,*) '         numq   = ',numq
            if(dowr) write(outfile,*) '         ihail  = ',ihail
            if(dowr) write(outfile,*) ' ----------------------------------------------- '
            if(dowr) write(outfile,*)

          endif

          if(ptype.eq.4)then

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) 'Calling lfoice_init'
            if(dowr) write(outfile,*)

            call lfoice_init(dt)

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) ' ----------------------------------------------- '
            if(dowr) write(outfile,*) '  ... using GSR LFO microphysics scheme ... '
            if(dowr) write(outfile,*) '         numq   = ',numq
            if(dowr) write(outfile,*) ' ----------------------------------------------- '
            if(dowr) write(outfile,*)

          endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.3)THEN    ! Thompson scheme

          iice = 1    ! this means that ptype=3 is an ice scheme
          idm  = 1    ! this means that ptype=3 has at least one double moment

          numq = 8    ! there are 8 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 8    ! the last number concentration var is the eighth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.
          cloudvar(7) = .false.
          cloudvar(8) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '
          qname(7) = 'nci'
          qname(8) = 'ncr'

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'
          qunit(4) = 'kg/kg'
          qunit(5) = 'kg/kg'
          qunit(6) = 'kg/kg'
          qunit(7) = '#/kg'
          qunit(8) = '#/kg'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2
          qmag( 4) = 1.0e-2
          qmag( 5) = 1.0e-2
          qmag( 6) = 1.0e-2
          qmag( 7) = 1.0e6
          qmag( 8) = 1.0e6

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the Thompson scheme -----
             ! zh is needed first .... see below !


!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.5)THEN    ! Morrison scheme

          !----- initialize the Morrison scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling GRAUPEL_INIT'

          call graupel_init(ihail,inum,ndcnst)

          if(dowr) write(outfile,*) 'Returned from GRAUPEL_INIT'
          if(dowr) write(outfile,*)

          !------------------------------------------

          iice = 1    ! this means that ptype=5 is an ice scheme
          idm  = 1    ! this means that ptype=5 has at least one double moment

        if(inum.eq.1)then
          ! constant cloud-drop concentration

          numq = 10   ! there are 10 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 10   ! the last number concentration var is the tenth array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'nci'
          qname( 8) = 'ncs'
          qname( 9) = 'ncr'
          qname(10) = 'ncg'

          qunit( 1) = 'kg/kg'
          qunit( 2) = 'kg/kg'
          qunit( 3) = 'kg/kg'
          qunit( 4) = 'kg/kg'
          qunit( 5) = 'kg/kg'
          qunit( 6) = 'kg/kg'
          qunit( 7) = '#/kg'
          qunit( 8) = '#/kg'
          qunit( 9) = '#/kg'
          qunit(10) = '#/kg'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2
          qmag( 4) = 1.0e-2
          qmag( 5) = 1.0e-2
          qmag( 6) = 1.0e-2
          qmag( 7) = 1.0e6
          qmag( 8) = 1.0e6
          qmag( 9) = 1.0e6
          qmag(10) = 1.0e6

        elseif(inum.eq.0)then
          ! cloud-droplet concentration is a predicted variable

          numq = 11   ! there are 11 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 11   ! the last number concentration var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'nci'
          qname( 8) = 'ncs'
          qname( 9) = 'ncr'
          qname(10) = 'ncg'
          qname(11) = 'ncc'

          qunit( 1) = 'kg/kg'
          qunit( 2) = 'kg/kg'
          qunit( 3) = 'kg/kg'
          qunit( 4) = 'kg/kg'
          qunit( 5) = 'kg/kg'
          qunit( 6) = 'kg/kg'
          qunit( 7) = '#/kg'
          qunit( 8) = '#/kg'
          qunit( 9) = '#/kg'
          qunit(10) = '#/kg'
          qunit(11) = '#/kg'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2
          qmag( 4) = 1.0e-2
          qmag( 5) = 1.0e-2
          qmag( 6) = 1.0e-2
          qmag( 7) = 1.0e5
          qmag( 8) = 1.0e5
          qmag( 9) = 1.0e5
          qmag(10) = 1.0e5
          qmag(11) = 1.0e5

        else

          print *,'  unrecognized value for inum '
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1

        endif

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  ... using Morrison microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) '         ihail  = ',ihail
        if(inum.eq.1)then
          if(dowr) write(outfile,*) '         assuming constant cloud droplet concentration' 
          if(dowr) write(outfile,*) '         ndcnst = ',ndcnst
        else
          if(dowr) write(outfile,*) '         predicting cloud droplet concentration'
        endif
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.6)THEN        ! Rotunno-Emanuel scheme

          numq = 2    ! there are 2 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 2    ! the last liquid variable is the second array

          cloudvar(1) = .false.
          cloudvar(2) = .true.

          qname(1) = 'qv '
          qname(2) = 'ql '

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'

          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  ... using Rotunno-Emanuel microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.7)THEN    !  WSM6

          ibtmp  = max( ibtmp  ,  1 )
          ietmp  = max( ietmp  , ni )
          kbtmp  = max( kbtmp  ,  1 )
          ketmp  = max( ketmp  , nk )
          ntmp2d = max( ntmp2d , 15 )
          ntmp1d = max( ntmp1d ,  7 )

          iice = 1    ! this means that ptype=7 is ice scheme

          numq = 6    ! there are 6 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'
          qunit(4) = 'kg/kg'
          qunit(5) = 'kg/kg'
          qunit(6) = 'kg/kg'

          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2
          qmag( 4) = 1.0e-2
          qmag( 5) = 1.0e-2
          qmag( 6) = 1.0e-2

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

            ! Note: using some values from wrf's module_model_constants
            if(myid.eq.0) print *
            if(myid.eq.0) print *,'  calling mp_wsm6_init ... '
            call mp_wsm6_init(den0=1.28,       &
                              denr=1000.0,     &
                              dens=100.0,      &
                              cl=4190.0,       &
                              cpv=cpv,         &
                              hail_opt=ihail,  &
                              errmsg=errmsg,   &
                              errflg=errflg)
            if(myid.eq.0) print *,'  ... done '
            if(myid.eq.0) print *

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF ( ptype .eq. 26 .or. ptype .eq. 27) THEN    ! NSSL variable density scheme (with no hail category)
!        ELSEIF ( ptype .eq. 26 ) THEN    ! NSSL variable density scheme (with no hail category)

          lv = 1
          lc = 2
          lr = 3
          li = 4
          ls = 5
          lgh = 6
          lhl = 0
          lccn = 7
          lvh = 0

          IF ( ptype .eq. 26 ) THEN
            nssl_hail_on = .false.
            lhab = lgh
          ELSE
            nssl_hail_on = .true.
            lhl = 7
            lhab = lhl
          ENDIF
          
          iice = 1    ! this means that ptype=26 is an ice scheme
          idm  = 1    ! this means that ptype=26 has at least one double moment
          idmplus = 1

          ! inucopt is not used yet
          i = 0
          IF ( inucopt == 1 .or. inucopt == 3 ) THEN
           i = i + 1
          ENDIF
          IF ( inucopt == 2 .or. inucopt == 3 ) THEN
           i = i + 1
          ENDIF
          
          j = 0
          IF ( nssl_3moment ) THEN
           j = 2 ! three extra variables for 3moment
          ENDIF
          
          k = 0
          IF ( nssl_density_on ) THEN
            k = 1
          ENDIF

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = lhab
          ltmp = nqs2+1
          nnc1 = ltmp    ! the first number concentration var (lccn)

          lccn = nnc1
          
          lnc = nnc1+1+i
          lnr = nnc1+2+i
          lni = nnc1+3+i
          lns = nnc1+4+i
          lnh = nnc1+5+i

          IF ( nssl_hail_on ) THEN
            lnhl = nnc1+6+i
            ltmp = lnhl
          ELSE
            ltmp = lnh
          ENDIF
          
          nnc2 = ltmp

          
          IF ( nssl_density_on ) THEN
            ltmp = ltmp + 1
            lvh  = ltmp
            nvl1 = ltmp ! the first particle volume var
            IF ( nssl_hail_on ) THEN
              ltmp = ltmp + 1
              lvhl = ltmp
            ENDIF
            nvl2 = ltmp ! the last particle volume var
          ENDIF
          
          IF ( nssl_3moment ) THEN
            lzr = ltmp+1
            nzl1 = lzr    ! the first reflectivity var (rain)
            lzh = ltmp+2
            ltmp = lzh
            
            IF ( nssl_hail_on ) THEN
              ltmp = ltmp + 1
              lzhl = ltmp
            ENDIF
            
            nzl2 = ltmp   ! the last reflectivity var
          ENDIF

          numq = ltmp   ! number of variables


          cloudvar(lv) = .false.
          cloudvar(lc) = .true.
          cloudvar(lr) = .false.
          cloudvar(li) = .true.
          cloudvar(ls) = .false.
          cloudvar(lgh) = .false.
          cloudvar( lhab:numq) = .false.


          qname(lv) = 'qv '
          qname(lc) = 'qc '
          qname(lr) = 'qr '
          qname(li) = 'qi '
          qname(ls) = 'qs '
          qname(lgh) = 'qg '
          IF ( nssl_hail_on ) qname(lhl) = 'qhl'
          qname(lccn) = 'ccn' ! CCN concentration

          i = 0
          lccna = 1
          lcina = 1
          IF ( inucopt == 1 .or. inucopt == 3 ) THEN
           i = i + 1
          qname(lccn+i) = 'cna' ! activated CCN concentration
          lccna = lccn+i
          f_ccna = .true.
          ENDIF
          IF ( inucopt == 2 .or. inucopt == 3 ) THEN
           i = i + 1
          qname(lccn+i) = 'cia' ! activated CIN concentration
          lcina = lccn+i
          f_cina = .true.
          ENDIF
          
          qname(lnc) = 'ccw' ! droplet conc
          qname(lnr) = 'crw' ! rain conc
          qname(lni) = 'cci' ! ice crystal conc
          qname(lns) = 'csw' ! snow conc
          qname(lnh) = 'chw' ! graupel conc
          IF ( nssl_hail_on ) qname(lnhl) = 'chl' ! graupel conc
          
          IF ( nssl_density_on ) THEN
           qname(lvh) = 'vhw' ! graupel volume
           IF ( nssl_hail_on ) qname(lvhl) = 'vhl' ! hail volume
          ENDIF
          
          IF ( nssl_3moment ) THEN
            qname(lzr) = 'zrw' ! rain reflectivity moment
            qname(lzh) = 'zhw' ! graupel reflectivity moment
            IF ( nssl_hail_on ) qname(lzhl) = 'zhl' ! hail reflectivity
          ENDIF

          qunit(lv:lhab) = 'kg/kg'

          qunit(nnc1:nnc2) = '#/kg'
          IF ( nssl_density_on ) THEN
            qunit(nvl1:nvl2) = 'm^3/kg'
          ENDIF
          IF ( nssl_3moment ) THEN
            qunit(nzl1:nzl2) = 'Z/m^3/kg'
          ENDIF

          rhovar(lv:lhab) = .false.
          rhovar(lhab+1:numq) = .true.

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 0.01
          qmag( 2) = 0.01
          qmag( 3) = 0.01
          qmag( 4) = 0.01
          qmag( 5) = 0.01
          qmag( 6) = 0.01
          qmag(nnc1:nnc2) = 1.e8
          qmag( lnc) = 1.0e10
          qmag( lni) = 3.0e9
          IF ( nssl_3moment ) qmag(nzl1:nzl2) = 1.0e-12
          IF ( nssl_density_on ) qmag(nvl1:nvl2) = 1.0e-5
          qmag(:) = 1.e-12


          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the NSSL scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling nssl_2mom_init'
          if(dowr) write(outfile,*)

         IF ( nssl_3moment ) THEN
!           CALL nssl_2mom_init(ipctmp=8,mixphase=0,ihvol=-1,eqtset_tmp=eqtset,myrank=myid,mpiroot=0)
           CALL nssl_2mom_init(ipctmp=8,mixphase=0,nssl_hail_on=nssl_hail_on,nssl_density_on=nssl_density_on, &
                     eqtset_tmp=eqtset,myrank=myid,mpiroot=0)
         ELSE
!           CALL nssl_2mom_init(ipctmp=5,mixphase=0,ihvol=-1,eqtset_tmp=eqtset,myrank=myid,mpiroot=0)
           CALL nssl_2mom_init(ipctmp=5,mixphase=0,nssl_hail_on=nssl_hail_on,nssl_density_on=nssl_density_on, &
                     eqtset_tmp=eqtset,myrank=myid,mpiroot=0)
         ENDIF

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Returned from nssl_2mom_init'
          if(dowr) write(outfile,*)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          IF ( nssl_3moment ) THEN
          if(dowr) write(outfile,*) '  using NSSL 3-moment scheme  '
          ELSE
          if(dowr) write(outfile,*) '  using NSSL 2-moment scheme  '
          ENDIF

          IF ( nssl_hail_on ) THEN
          if(dowr) write(outfile,*) '  (with hail ) '
          ELSE
          if(dowr) write(outfile,*) '  (graupel-only, no hail ) '
          ENDIF
          IF ( nssl_density_on ) THEN
          if(dowr) write(outfile,*) '  (graupel/hail density predicted ) '
          ENDIF
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)


!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF ( ptype .eq. 28 ) THEN    ! single moment NSSL scheme (without hail)

          iice = 1    ! this means that ptype=28 is an ice scheme
          idm  = 0    ! this means that ptype=28 has at least one double moment
          idmplus = 1

          numq = 6   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = 6    ! the last solid variable

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '

          qunit( 1) = 'kg/kg'
          qunit( 2) = 'kg/kg'
          qunit( 3) = 'kg/kg'
          qunit( 4) = 'kg/kg'
          qunit( 5) = 'kg/kg'
          qunit( 6) = 'kg/kg'

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag(1) = 0.01
          qmag(2) = 0.01
          qmag(3) = 0.01
          qmag(4) = 0.01
          qmag(5) = 0.01
          qmag(6) = 0.01


          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the NSSL scheme -----

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Calling graupel_init'
!          if(dowr) write(outfile,*)

         CALL nssl_2mom_init(ipctmp=0,mixphase=0,ihvol=-1,eqtset_tmp=eqtset)


!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Returned from graupel_init'
!          if(dowr) write(outfile,*)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  using NSSL 1-moment scheme  '
          if(dowr) write(outfile,*) '  (graupel only; fixed graupel density) '
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 )THEN    ! P3

          !----- initialize the P3 scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling P3_INIT'

          if(     ptype.eq.50 )then
            call p3_init(lookup_file_dir='.',  &
                         nCat=1,               &
                         trplMomI=.false.,     &
                         model='WRF',          &
                         stat=p3stat,          &
                         abort_on_err=.false.)
          elseif( ptype.eq.51 )then
            call p3_init(lookup_file_dir='.',  &
                         nCat=1,               &
                         trplMomI=.false.,     &
                         model='WRF',          &
                         stat=p3stat,          &
                         abort_on_err=.false.)
          elseif( ptype.eq.52 )then
            call p3_init(lookup_file_dir='.',  &
                         nCat=2,               &
                         trplMomI=.false.,     &
                         model='WRF',          &
                         stat=p3stat,          &
                         abort_on_err=.false.)
          elseif( ptype.eq.53 )then
            call p3_init(lookup_file_dir='.',  &
                         nCat=1,               &
                         trplMomI=.true.,      &
                         model='WRF',          &
                         stat=p3stat,          &
                         abort_on_err=.false.)
          endif

          if(dowr) write(outfile,*) 'Returned from P3_INIT'
          if(dowr) write(outfile,*)

          if( p3stat.ne.0 )then
            print *,'  there was an error in p3_init '
            call stopcm1
          endif

          !------------------------------------------

          iice = 1    ! this scheme has ice microphysics
          idm  = 1    ! this scheme has at least one double moment


      !cccccccccccccccccccccccccccccccccccccccccccccccccccc!
        IF( ptype.eq.50 )THEN

          numq = 8    ! there are 8 q variables

          nqv  = 1    ! qv
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = 4    ! the last solid variable

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.
          cloudvar(7) = .false.
          cloudvar(8) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qni'
          qname(6) = 'qnr'
          qname(7) = 'qir'
          qname(8) = 'qib'

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'
          qunit(4) = 'kg/kg'
          qunit(5) = '#/kg'
          qunit(6) = '#/kg'
          qunit(7) = 'kg/kg'
          qunit(8) = 'm^-3 kg^-1'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag(1) = 0.01
          qmag(2) = 0.01
          qmag(3) = 0.01
          qmag(4) = 0.01
          qmag(5) = 1.0e5
          qmag(6) = 1.0e5
          qmag(7) = 0.01
          qmag(8) = 1.0e-5

          ! for p3 arrays:
          np3a = 22
          np3o =  3

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc!
        ELSEIF( ptype.eq.51 )THEN

          numq = 9    ! there are 9 q variables

          nqv  = 1    ! qv
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = 4    ! the last solid variable

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.
          cloudvar(7) = .false.
          cloudvar(8) = .false.
          cloudvar(9) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qni'
          qname(6) = 'qnr'
          qname(7) = 'qir'
          qname(8) = 'qib'
          qname(9) = 'qnc'

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'
          qunit(4) = 'kg/kg'
          qunit(5) = '#/kg'
          qunit(6) = '#/kg'
          qunit(7) = 'kg/kg'
          qunit(8) = 'm^-3 kg^-1'
          qunit(9) = '#/kg'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag(1) = 0.01
          qmag(2) = 0.01
          qmag(3) = 0.01
          qmag(4) = 0.01
          qmag(5) = 1.0e5
          qmag(6) = 1.0e5
          qmag(7) = 0.01
          qmag(8) = 1.0e-5
          qmag(9) = 1.0e8

          ! for p3 arrays:
          np3a = 23
          np3o =  3

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc!
        ELSEIF( ptype.eq.52 )THEN

          numq = 13   ! there are 13 q variables

          nqv  = 1    ! qv
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = 5    ! the last solid variable

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .true.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qi2'
          qname( 6) = 'nc '
          qname( 7) = 'nr '
          qname( 8) = 'ni '
          qname( 9) = 'ri '
          qname(10) = 'bi '
          qname(11) = 'ni2'
          qname(12) = 'ri2'
          qname(13) = 'bi2'

          qunit( 1) = 'kg/kg'
          qunit( 2) = 'kg/kg'
          qunit( 3) = 'kg/kg'
          qunit( 4) = 'kg/kg'
          qunit( 5) = 'kg/kg'
          qunit( 6) = '#/kg'
          qunit( 7) = '#/kg'
          qunit( 8) = '#/kg'
          qunit( 9) = 'kg/kg'
          qunit(10) = 'm^-3 kg^-1'
          qunit(11) = '#/kg'
          qunit(12) = 'kg/kg'
          qunit(13) = 'm^-3 kg^-1'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 0.01
          qmag( 2) = 0.01
          qmag( 3) = 0.01
          qmag( 4) = 0.01
          qmag( 5) = 0.01
          qmag( 6) = 1.0e5
          qmag( 7) = 1.0e5
          qmag( 8) = 1.0e5
          qmag( 9) = 0.01
          qmag(10) = 1.0e-5
          qmag(11) = 1.0e5
          qmag(12) = 0.01
          qmag(13) = 1.0e-5

          ! for p3 arrays:
          np3a = 30
          np3o =  6

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc!

        ELSEIF( ptype.eq.53 )THEN

          numq = 10   ! there are 9 q variables

          nqv  = 1    ! qv
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = 4    ! the last solid variable

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.
          cloudvar(7) = .false.
          cloudvar(8) = .false.
          cloudvar(9) = .false.
          cloudvar(10) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qni'
          qname(6) = 'qnr'
          qname(7) = 'qir'
          qname(8) = 'qib'
          qname(9) = 'qnc'
          qname(10) = 'qzi'

          qunit(1) = 'kg/kg'
          qunit(2) = 'kg/kg'
          qunit(3) = 'kg/kg'
          qunit(4) = 'kg/kg'
          qunit(5) = '#/kg'
          qunit(6) = '#/kg'
          qunit(7) = 'kg/kg'
          qunit(8) = 'm^3 kg^-1'
          qunit(9) = '#/kg'
          qunit(10) = 'm(6) kg(-1)'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag(1) = 0.01
          qmag(2) = 0.01
          qmag(3) = 0.01
          qmag(4) = 0.01
          qmag(5) = 1.0e5
          qmag(6) = 1.0e5
          qmag(7) = 0.01
          qmag(8) = 1.0e-5
          qmag(9) = 1.0e8
          qmag(10) = 1.0e-3

          ! for p3 arrays:
          np3a = 24
          np3o =  3

        ENDIF

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc!

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  ... using P3 microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.55)THEN        ! Jensen_ISHMAEL

          !----- initialize the Jensen_AHAB scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling jensen_ahab_init '

          call jensen_ishmael_init(myid)

          if(dowr) write(outfile,*) 'Returned from jensen_ahab_init '
          if(dowr) write(outfile,*)

          !------------------------------------------

          iice = 1    ! this means that ptype=8 is an ice scheme
          idm  = 1    ! this means that ptype=8 has at least one double moment

          numq = 16   ! there are 10 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable
          nql2 = 3    ! the last liquid variable
          nqs1 = 4    ! the first solid variable
          nqs2 = 6    ! the last solid variable
          nnc1 = 7    ! the first number concentration var
          nnc2 = 10   ! the last number concentration var
          nvl1 = 11   ! the first particle volume var
          nvl2 = 16   ! the last particle volume var

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .true.
          cloudvar( 6) = .true.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.
          cloudvar(14) = .false.
          cloudvar(15) = .false.
          cloudvar(16) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi1'
          qname( 5) = 'qi2'
          qname( 6) = 'qi3'
          qname( 7) = 'nr '
          qname( 8) = 'ni1'
          qname( 9) = 'ni2'
          qname(10) = 'ni3'
          qname(11) = 'ai1'
          qname(12) = 'ai2'
          qname(13) = 'ai3'
          qname(14) = 'ci1'
          qname(15) = 'ci2'
          qname(16) = 'ci3'

          qunit( 1) = 'kg/kg'
          qunit( 2) = 'kg/kg'
          qunit( 3) = 'kg/kg'
          qunit( 4) = 'kg/kg'
          qunit( 5) = 'kg/kg'
          qunit( 6) = 'kg/kg'
          qunit( 7) = '#/kg'
          qunit( 8) = '#/kg'
          qunit( 9) = '#/kg'
          qunit(10) = '#/kg'
          qunit(11) = 'm^3/kg'
          qunit(12) = 'm^3/kg'
          qunit(13) = 'm^3/kg'
          qunit(14) = 'm^3/kg'
          qunit(15) = 'm^3/kg'
          qunit(16) = 'm^3/kg'

          ! likely maximum value (order-of-magnitude)
          ! (needed for monotonic advection schemes)
          qmag( 1) = 1.0e-2
          qmag( 2) = 1.0e-2
          qmag( 3) = 1.0e-2
          qmag( 4) = 1.0e-2
          qmag( 5) = 1.0e-2
          qmag( 6) = 1.0e-2
          qmag( 7) = 1.0e6
          qmag( 8) = 1.0e6
          qmag( 9) = 1.0e6
          qmag(10) = 1.0e6
          qmag(11) = 1.0e-8
          qmag(12) = 1.0e-8
          qmag(13) = 1.0e-8
          qmag(14) = 1.0e-8
          qmag(15) = 1.0e-8
          qmag(16) = 1.0e-8

          ! use p3 arrays for diagnostic output:
          np3a =  1
          np3o = 12

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  ... using Jensen_AHAB microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)


!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  insert new ptype here

!!!        ELSEIF(ptype.eq.xxx)THEN    ! new microphysics scheme

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

        ELSE

          IF(myid.eq.0)THEN
            print *
            print *,'  ptype = ',ptype
            print *
            print *,'  Unrecognized value for ptype '
            print *
            print *,'  ... stopping cm1 ... '
            print *
          ENDIF

#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1

        ENDIF    ! endif for ptype

      ENDIF    ! endif for imoist=1

!-----------------------------------------------------------------------
!-------   END:  modify stuff above here -------------------------------
!-----------------------------------------------------------------------

      IF( (radopt.eq.1.or.radopt.eq.2) .and. iice.ne.1 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  iice     = ',iice
        print *
        print *,'  radopt=1,2 requires an ice microphysics scheme '
        print *
        print *,'   stopping model .... '
        print *
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

!-----------------------------------------------------------------------

      nqc = 0
      nqr = 0
      nqi = 0
      nqs = 0
      nqg = 0
      nnci = 0
      nncc = 0
      nqi2 = 0
      nqi3 = 0

      do n=1,numq
        if( qname(n).eq.'qc ' .or. qname(n).eq.'ql ' ) nqc = n
        if( qname(n).eq.'qr ' ) nqr = n
        if( qname(n).eq.'qi ' .or. qname(n).eq.'qi1' ) nqi = n
        if( qname(n).eq.'qs ' ) nqs = n
        if( qname(n).eq.'qg ' ) nqg = n
        if( qname(n).eq.'nci' ) nnci = n
        if( qname(n).eq.'ncc' ) nncc = n
        if( qname(n).eq.'qi2' ) nqi2 = n
        if( qname(n).eq.'qi3' ) nqi3 = n
      enddo

      if(numq .gt. maxq)then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  WARNING!   numq > maxq'
        if(dowr) write(outfile,*) '  You need to increase maxq in input.incl and recompile'
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Stopping model ....'
        if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'iice      =',iice
      if(dowr) write(outfile,*) 'idm       =',idm
      if(dowr) write(outfile,*) 'idmplus   =',idmplus
      if(dowr) write(outfile,*) 'numq      =',numq
      if(dowr) write(outfile,*) 'nqv       =',nqv
      if(dowr) write(outfile,*) 'nqc       =',nqc
      if(dowr) write(outfile,*) 'nqr       =',nqr
      if(dowr) write(outfile,*) 'nqi       =',nqi
      if(dowr) write(outfile,*) 'nqi2      =',nqi2
      if(dowr) write(outfile,*) 'nqi3      =',nqi3
      if(dowr) write(outfile,*) 'nqs       =',nqs
      if(dowr) write(outfile,*) 'nqg       =',nqg
      if(dowr) write(outfile,*) 'nnci      =',nnci
      if(dowr) write(outfile,*) 'nncc      =',nncc
      if(dowr) write(outfile,*) 'nql1      =',nql1
      if(dowr) write(outfile,*) 'nql2      =',nql2
      if(dowr) write(outfile,*) 'nqs1      =',nqs1
      if(dowr) write(outfile,*) 'nqs2      =',nqs2
      if(dowr) write(outfile,*) 'nnc1      =',nnc1
      if(dowr) write(outfile,*) 'nnc2      =',nnc2
      if(dowr) write(outfile,*) 'nvl1      =',nvl1
      if(dowr) write(outfile,*) 'nvl2      =',nvl2
      if(dowr) write(outfile,*) 'nzl1      =',nzl1
      if(dowr) write(outfile,*) 'nzl2      =',nzl2
      if(dowr) write(outfile,*)

    if( imoist.eq.1 )then
      if(dowr) write(outfile,174)
174   format('   n  qname    qunit                 qmag ')
      do n=1,numq
        if(dowr) write(outfile,173) n,qname(n),qunit(n),qmag(n)
173     format(2x,i3,2x,a3,2x,a20,2x,es12.4)
      enddo
      if(dowr) write(outfile,*)
    endif

!--------------------------------------------------------------

      iterrain = 0
      if(terrain_flag) iterrain = 1

      output_interp  = max(0,min(1,output_interp))*iterrain
      output_rain    = max(0,min(1,output_rain))
      output_sws     = max(0,min(1,output_sws))
      output_svs     = max(0,min(1,output_svs))
      output_sps     = max(0,min(1,output_sps))
      output_srs     = max(0,min(1,output_srs))
      output_sgs     = max(0,min(1,output_sgs))
      output_sus     = max(0,min(1,output_sus))
      output_shs     = max(0,min(1,output_shs))
      output_coldpool= max(0,min(1,output_coldpool))
      output_sfcflx  = max(0,min(1,output_sfcflx))
      output_sfcparams = max(0,min(1,output_sfcparams))
      output_sfcdiags = max(0,min(1,output_sfcdiags))
      output_psfc    = max(0,min(1,output_psfc))
      output_zs      = max(0,min(1,output_zs))*iterrain
      output_zh      = max(0,min(1,output_zh))
      output_basestate = max(0,min(1,output_basestate))
      output_th      = max(0,min(1,output_th))
      output_thpert  = max(0,min(1,output_thpert))
      output_prs     = max(0,min(1,output_prs))
      output_prspert = max(0,min(1,output_prspert))
      output_pi      = max(0,min(1,output_pi))
      output_pipert  = max(0,min(1,output_pipert))
      output_rho     = max(0,min(1,output_rho))
      output_rhopert = max(0,min(1,output_rhopert))
      output_tke     = max(0,min(1,output_tke))
      output_km      = max(0,min(1,output_km))
      output_kh      = max(0,min(1,output_kh))
      output_qv      = max(0,min(1,output_qv))
      output_qvpert  = max(0,min(1,output_qvpert))
      output_q       = max(0,min(1,output_q))
      output_dbz     = max(0,min(1,output_dbz))
      output_buoyancy= max(0,min(1,output_buoyancy))
      output_u       = max(0,min(1,output_u))
      output_upert   = max(0,min(1,output_upert))
      output_uinterp = max(0,min(1,output_uinterp))
      output_v       = max(0,min(1,output_v))
      output_vpert   = max(0,min(1,output_vpert))
      output_vinterp = max(0,min(1,output_vinterp))
      output_w       = max(0,min(1,output_w))
      output_winterp = max(0,min(1,output_winterp))
      output_vort    = max(0,min(1,output_vort))
      output_pv      = max(0,min(1,output_pv))
      output_uh      = max(0,min(1,output_uh))
      output_pblten  = max(0,min(1,output_pblten))
      output_dissten = max(0,min(1,output_dissten))
      output_fallvel  = max(0,min(1,output_fallvel))
      output_nm      = max(0,min(1,output_nm))
      output_def     = max(0,min(1,output_def))
      output_radten  = max(0,min(1,output_radten))
      output_cape = max(0,min(1,output_cape))
      output_cin = max(0,min(1,output_cin))
      output_lcl = max(0,min(1,output_lcl))
      output_lfc = max(0,min(1,output_lfc))
      output_lwp = max(0,min(1,output_lwp))
      output_pwat = max(0,min(1,output_pwat))
      output_thbudget = max(0,min(1,output_thbudget))
      output_qvbudget = max(0,min(1,output_qvbudget))
      output_ubudget = max(0,min(1,output_ubudget))
      output_vbudget = max(0,min(1,output_vbudget))
      output_wbudget = max(0,min(1,output_wbudget))
      prcl_dbz       = max(0,min(1,prcl_dbz))

      nrain = 1
      if(imove.eq.1) nrain = 2

      if(dowr) write(outfile,*) 'nrain     =',nrain
      if(dowr) write(outfile,*)

!!!      if( restart_file_dbz .or. prcl_dbz.ge.1 )then
!!!        if( ptype.ne.4 )then
!!!          output_dbz = 1
!!!        endif
!!!      endif
!!!      if( ptype.eq.3 .or. ptype.eq.5 )then
!!!        output_dbz = 1
!!!      endif
      if(imoist.eq.0)then
        output_rain=0
        output_srs=0
        output_sgs=0
        output_qv=0
        output_qvpert=0
        output_q=0
        output_dbz=0
        output_cape=0
        output_cin=0
        output_lcl=0
        output_lfc=0
        output_lwp=0
        output_pwat=0
        output_qvbudget=0
      endif
      if( nqr.le.0 ) output_srs = 0
      if( nqg.le.0 ) output_sgs = 0
      if( sgsmodel.le.0 .or. sgsmodel.eq.2 .or. sgsmodel.eq.6 )then
        output_tke=0
      endif
      if( sgsmodel.le.0 .and. ipbl.eq.0 )then
        output_km=0
        output_kh=0
      endif
      if(ipbl.lt.1)then
        output_pblten=0
      endif
      if(radopt.eq.0)then
        output_radten=0
      endif
      if(radopt.ge.1)then
        output_radten=1
      endif
      if( bbc.ne.3 .and. isfcflx.eq.0 .and. radopt.eq.0 )then
!!!        output_sfcflx = 0
        output_sfcparams = 0
        output_sfcdiags = 0
      endif
      if(terrain_flag)then
        output_zs = 1
        output_zh = 1
      endif
      getfall = .false.
      if( ptype.eq.3 .or. ptype.eq.5 ) getfall = .true.
      if( imoist.eq.0 .or. (.not.getfall) )then
        output_fallvel = 0
      endif

      if( outunits.eq.2 )then
        ! units of distance in output file are meters
        outunitconv = 1.0
        aunit = 'meters'
      else
        ! units of distance in output file are km
        outunitconv = 0.001
        aunit = 'km'
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'output_format    =',output_format
      if(dowr) write(outfile,*) 'output_filetype  =',output_filetype
      if(dowr) write(outfile,*) 'output_interp    =',output_interp
      if(dowr) write(outfile,*) 'output_rain      =',output_rain
      if(dowr) write(outfile,*) 'output_sws       =',output_sws
      if(dowr) write(outfile,*) 'output_svs       =',output_svs
      if(dowr) write(outfile,*) 'output_sps       =',output_sps
      if(dowr) write(outfile,*) 'output_srs       =',output_srs
      if(dowr) write(outfile,*) 'output_sgs       =',output_sgs
      if(dowr) write(outfile,*) 'output_sus       =',output_sus
      if(dowr) write(outfile,*) 'output_shs       =',output_shs
      if(dowr) write(outfile,*) 'output_coldpool  =',output_coldpool
      if(dowr) write(outfile,*) 'output_sfcflx    =',output_sfcflx
      if(dowr) write(outfile,*) 'output_sfcparams =',output_sfcparams
      if(dowr) write(outfile,*) 'output_sfcdiags  =',output_sfcdiags
      if(dowr) write(outfile,*) 'output_psfc      =',output_psfc
      if(dowr) write(outfile,*) 'output_zs        =',output_zs
      if(dowr) write(outfile,*) 'output_zh        =',output_zh
      if(dowr) write(outfile,*) 'output_basestate =',output_basestate
      if(dowr) write(outfile,*) 'output_th        =',output_th
      if(dowr) write(outfile,*) 'output_thpert    =',output_thpert
      if(dowr) write(outfile,*) 'output_prs       =',output_prs
      if(dowr) write(outfile,*) 'output_prspert   =',output_prspert
      if(dowr) write(outfile,*) 'output_pi        =',output_pi
      if(dowr) write(outfile,*) 'output_pipert    =',output_pipert
      if(dowr) write(outfile,*) 'output_rho       =',output_rho
      if(dowr) write(outfile,*) 'output_rhopert   =',output_rhopert
      if(dowr) write(outfile,*) 'output_tke       =',output_tke
      if(dowr) write(outfile,*) 'output_km        =',output_km
      if(dowr) write(outfile,*) 'output_kh        =',output_kh
      if(dowr) write(outfile,*) 'output_qv        =',output_qv
      if(dowr) write(outfile,*) 'output_qvpert    =',output_qvpert
      if(dowr) write(outfile,*) 'output_q         =',output_q
      if(dowr) write(outfile,*) 'output_dbz       =',output_dbz
      if(dowr) write(outfile,*) 'output_buoyancy  =',output_buoyancy
      if(dowr) write(outfile,*) 'output_u         =',output_u
      if(dowr) write(outfile,*) 'output_upert     =',output_upert
      if(dowr) write(outfile,*) 'output_uinterp   =',output_uinterp
      if(dowr) write(outfile,*) 'output_v         =',output_v
      if(dowr) write(outfile,*) 'output_vpert     =',output_vpert
      if(dowr) write(outfile,*) 'output_vinterp   =',output_vinterp
      if(dowr) write(outfile,*) 'output_w         =',output_w
      if(dowr) write(outfile,*) 'output_winterp   =',output_winterp
      if(dowr) write(outfile,*) 'output_vort      =',output_vort
      if(dowr) write(outfile,*) 'output_pv        =',output_pv
      if(dowr) write(outfile,*) 'output_uh        =',output_uh
      if(dowr) write(outfile,*) 'output_pblten    =',output_pblten
      if(dowr) write(outfile,*) 'output_dissten   =',output_dissten
      if(dowr) write(outfile,*) 'output_fallvel   =',output_fallvel
      if(dowr) write(outfile,*) 'output_nm        =',output_nm
      if(dowr) write(outfile,*) 'output_def       =',output_def
      if(dowr) write(outfile,*) 'output_radten    =',output_radten
      if(dowr) write(outfile,*) 'output_cape      =',output_cape
      if(dowr) write(outfile,*) 'output_cin       =',output_cin
      if(dowr) write(outfile,*) 'output_lcl       =',output_lcl
      if(dowr) write(outfile,*) 'output_lfc       =',output_lfc
      if(dowr) write(outfile,*) 'output_lwp       =',output_lwp
      if(dowr) write(outfile,*) 'output_pwat      =',output_pwat
      if(dowr) write(outfile,*) 'output_thbudget  =',output_thbudget
      if(dowr) write(outfile,*) 'output_qvbudget  =',output_qvbudget
      if(dowr) write(outfile,*) 'output_ubudget   =',output_ubudget
      if(dowr) write(outfile,*) 'output_vbudget   =',output_vbudget
      if(dowr) write(outfile,*) 'output_wbudget   =',output_wbudget
      if(dowr) write(outfile,*) 'output_pdcomp    =',output_pdcomp
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'restart_format       = ',restart_format
      if(dowr) write(outfile,*) 'restart_filetype     = ',restart_filetype
      if(dowr) write(outfile,*) 'restart_file_theta   = ',restart_file_theta
      if(dowr) write(outfile,*) 'restart_file_dbz     = ',restart_file_dbz
      if(dowr) write(outfile,*) 'restart_file_th0     = ',restart_file_th0
      if(dowr) write(outfile,*) 'restart_file_prs0    = ',restart_file_prs0
      if(dowr) write(outfile,*) 'restart_file_pi0     = ',restart_file_pi0
      if(dowr) write(outfile,*) 'restart_file_rho0    = ',restart_file_rho0
      if(dowr) write(outfile,*) 'restart_file_qv0     = ',restart_file_qv0
      if(dowr) write(outfile,*) 'restart_file_u0      = ',restart_file_u0
      if(dowr) write(outfile,*) 'restart_file_v0      = ',restart_file_v0
      if(dowr) write(outfile,*) 'restart_file_zs      = ',restart_file_zs
      if(dowr) write(outfile,*) 'restart_file_zh      = ',restart_file_zh
      if(dowr) write(outfile,*) 'restart_file_zf      = ',restart_file_zf
      if(dowr) write(outfile,*) 'restart_file_diags   = ',restart_file_diags
      if(dowr) write(outfile,*) 'restart_use_theta    = ',restart_use_theta
      if(dowr) write(outfile,*) 'restart_reset_frqtim = ',restart_reset_frqtim
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Large-scale vertical velocity:

      dolsw = .false.

      IF( testcase.eq.3 ) dolsw = .true.
      IF( testcase.eq.4 ) dolsw = .true.
      IF( testcase.eq.5 ) dolsw = .true.
      IF( testcase.eq.6 ) dolsw = .true.
      IF( testcase.eq.7 ) dolsw = .true.

      if( dolsw )then
        if( terrain_flag  .or.  axisymm.eq.1  )then
          print *
          print *,'  Cannot use large-scale vertical velocity with terrain or axisymmetric model '
          print *
          print *,'  98724 '
          print *
          call stopcm1
        endif
      endif

      IF( dolsw )THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  dolsw      = ',dolsw
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  ... Including large-scale vertical velocity ... '
        if(dowr) write(outfile,*)
      ENDIF

!--------------------------------------------------------------
!  Define dimensions for allocatable arrays

      ibb=1
      ieb=1
      jbb=1
      jeb=1
      kbb=1
      keb=1

      use_pbl = .false.

      ibmynn=1
      iemynn=1
      jbmynn=1
      jemynn=1
      kbmynn=1
      kemynn=1

      if(imoist.eq.1)then
        ibm=ib
        iem=ie
        jbm=jb
        jem=je
        kbm=kb
        kem=ke
        if(ptype.ge.26)then
          ibzvd=ib
          iezvd=ie
          jbzvd=jb
          jezvd=je
          kbzvd=kb
          kezvd=ke
          nqzvd = numq + 1
        else
          ibzvd=1
          iezvd=1
          jbzvd=1
          jezvd=1
          kbzvd=1
          kezvd=1
          nqzvd=1
        endif
      else
        ibm=1
        iem=1
        jbm=1
        jem=1
        kbm=1
        kem=1
        ibzvd=1
        iezvd=1
        jbzvd=1
        jezvd=1
        kbzvd=1
        kezvd=1
        nqzvd=1
      endif

      if(iice.eq.1)then
        ibi=ib
        iei=ie
        jbi=jb
        jei=je
        kbi=kb
        kei=ke
      else
        ibi=1
        iei=1
        jbi=1
        jei=1
        kbi=1
        kei=1
      endif

      if( radopt.ge.1 .or. testcase.eq.4 .or. testcase.eq.5 )then
        ibr=ib
        ier=ie
        jbr=jb
        jer=je
        kbr=kb
        ker=ke
      else
        ibr=1
        ier=1
        jbr=1
        jer=1
        kbr=1
        ker=1
      endif

      if( (cm1setup.eq.2.or.cm1setup.eq.4) .and. ipbl.ge.1 )then
        ! for boundary-layer schemes:
        ibb=ib
        ieb=ie
        jbb=jb
        jeb=je
        kbb=kb
        keb=ke
        use_pbl = .true.
      endif

      if( ( idopbl .and. ( ipbl.eq.4 .or. ipbl.eq.5 ) ) .or. sfcmodel.eq.6 )then
        ! for MYNN boundary-layer:
        ibmynn=ib
        iemynn=ie
        jbmynn=jb
        jemynn=je
        kbmynn=kb
        kemynn=ke
        iusetke = .true.
        ibtmp=max( ibtmp , 1 )
        ietmp=max( ietmp , ni )
        kbtmp=max( kbtmp , 1 )
        ketmp=max( ketmp , nk+1 )
        ntmp2d=max( ntmp2d , 59 )
        ntmp1d=max( ntmp1d ,  2 )
      endif

      if( idopbl .and. ipbl.eq.1 )then
        ! YSU
        ibtmp = max( ibtmp , 1 )
        ietmp = max( ietmp , ni )
        kbtmp = max( kbtmp , 1 )
        ketmp = max( ketmp , nk+1 )
        ntmp2d = max( ntmp2d , 44 )
        ntmp1d = max( ntmp1d ,  5 )
      endif

      !-----------------------------------------
      if( sfcmodel.eq.6 .or. ipbl.eq.4 .or. ipbl.eq.5 )then
        bl_mynn_cloudpdf   =   2      ! option to switch to different cloud PDFs to represent subgrid clouds
        icloud_bl          =   1      ! option to couple the subgrid-scale clouds from the PBL scheme (MYNN only) to radiation schemes
        spp_pbl            =   0      ! Perturb parameters of MYNN PBL scheme
        grav_settling      =   0      ! gravitational settling of fog/cloud droplets
        bl_mynn_tkebudget  =   1      ! adds MYNN tke budget terms to output
        bl_mynn_mixlength  =   1      ! option to change mixing length formulation in MYNN
        bl_mynn_edmf       =   1      ! 1 activates mass-flux scheme in MYNN
        bl_mynn_edmf_mom   =   1      ! activates momentum transport in MYNN mass-flux scheme
        bl_mynn_edmf_tke   =   0      ! activates TKE transport in MYNN mass-flux scheme
        bl_mynn_cloudmix   =   1      ! 1 activates mixing of qc and qi in MYNN
        bl_mynn_mixscalars =   0      ! activate mixing of scalars (qnx, qnxfa) in MYNN
        bl_mynn_mixqt      =   0      ! 0:mix moisture species separate,1: mix total water
        bl_mynn_output     =   1      ! 1:Allocate and output extra 3D arrays
        bl_mynn_tkeadvect  =  .true.  ! tke advection?  (T=yes)
        if( ipbl.eq.4 )then
          bl_mynn_closure    =  2.50_sp   ! flag to determine the closure level for the mynn
        elseif( ipbl.eq.5 )then
          bl_mynn_closure    =   3.00_sp  ! flag to determine the closure level for the mynn
        endif
        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  bl_mynn_closure = ',bl_mynn_closure
        if(myid.eq.0) print *
      endif
      !-----------------------------------------

      if( ipbl.eq.6 .or. sfcmodel.eq.7 )then
        ! for MYJ:
        ibmyj=ib
        iemyj=ie
        jbmyj=jb
        jemyj=je
        kbmyj=kb
        kemyj=ke+1

        ibtmp=max( ibtmp , 1 )
        ietmp=max( ietmp , ni )
        kbtmp=max( kbtmp , 1 )
        ketmp=max( ketmp , nk+1 )
        ntmp2d=max( ntmp2d , 28 )
      else
        ibmyj=1
        iemyj=1
        jbmyj=1
        jemyj=1
        kbmyj=1
        kemyj=1
      endif

      if( (sfcmodel.ge.1) .or. (oceanmodel.eq.2) .or. (ipbl.ge.1) .or. (bbc.eq.3) )then
        ibl=ib
        iel=ie
        jbl=jb
        jel=je
      else
        ibl=1
        iel=1
        jbl=1
        jel=1
      endif

      if( sfcmodel.eq.3 )then
        ibtmp = max( ibtmp , 1 )
        ietmp = max( ietmp , ni )
        kbtmp = max( kbtmp , 1 )
        ketmp = max( ketmp , nk )
        ntmp2d = max( ntmp2d ,  6 )
        ntmp1d = max( ntmp1d ,  3 )
      endif

      if( idoles .or. idopbl .or. ipbl.ge.1 .or. cm1setup.eq.3 )then
        iusekm = 1
        iusekh = 1
        ibc=ib
        iec=ie
        jbc=jb
        jec=je
        kbc=kb
        kec=ke+1
      else
        iusekm = 0
        iusekh = 0
        ibc=1
        iec=1
        jbc=1
        jec=1
        kbc=1
        kec=1
      endif

      if( idoles .and. ( sgsmodel.eq.1 .or. sgsmodel.eq.3 .or. sgsmodel.eq.4 .or. sgsmodel.eq.5 ) )then
        iusetke = .true.
        ibt=ib
        iet=ie
        jbt=jb
        jet=je
        kbt=kb
        ket=ke+1
      else
        ibt=1
        iet=1
        jbt=1
        jet=1
        kbt=1
        ket=1
      endif

      if( sgsmodel.eq.5 .or. sgsmodel.eq.6 )then
        ibnba=ib
        ienba=ie
        jbnba=jb
        jenba=je
        kbnba=kb
        kenba=ke
      else
        ibnba=1
        ienba=1
        jbnba=1
        jenba=1
        kbnba=1
        kenba=1
      endif

      if( cuparam.ge.1 )then
        ibcu=ib
        iecu=ie
        jbcu=jb
        jecu=je
        kbcu=kb
        kecu=ke
        ibtmp = max( ibtmp ,    1 )
        ietmp = max( ietmp ,   ni )
        kbtmp = max( kbtmp ,    1 )
        ketmp = max( ketmp , nk+1 )
        ntmp2d = max( ntmp2d , 20 )
        ntmp1d = max( ntmp1d ,  3 )
      else
        ibcu=1
        iecu=1
        jbcu=1
        jecu=1
        kbcu=1
        kecu=1
      endif

      if(iptra.eq.1)then
        ibp=ib
        iep=ie
        jbp=jb
        jep=je
        kbp=kb
        kep=ke
      else
        ibp=1
        iep=1
        jbp=1
        jep=1
        kbp=1
        kep=1
      endif

      if( psolver.eq.4 .or. psolver.eq.5 .or. ibalance.eq.2 .or. pdcomp )then

        imirror = 0
        jmirror = 0

        ipb=1
        ipe=ni

        jpb=1
        jpe=nj

        if( (wbc.eq.2.or.wbc.eq.3).or.(ebc.eq.2.or.ebc.eq.3) )then

          imirror = 1
          ipe = ni*2

        endif

        if( (sbc.eq.2.or.sbc.eq.3).or.(nbc.eq.2.or.nbc.eq.3) )then

          jmirror = 1
          jpe = nj*2

        endif

        kpb=0
        kpe=nk+1

      else

        ipb=1
        ipe=1
        jpb=1
        jpe=1
        kpb=1
        kpe=1

      endif

      if( psolver.eq.4 .or. psolver.eq.5 .or. psolver.eq.6 )then
        ibph=ib
        ieph=ie
        jbph=jb
        jeph=je
        kbph=kb
        keph=ke
      else
        ibph=1
        ieph=1
        jbph=1
        jeph=1
        kbph=1
        keph=1
      endif

      ! array indices for P3 microphysics:
      ! (use these arrays for ISHMAEL diag output, too)
      if( imoist.eq.1 .and. ( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 .or. ptype.eq.55 ) )then
        ibp3 = 1
        iep3 = ni
        jbp3 = 1
        jep3 = nj
        kbp3 = 1
        kep3 = nk
      else
        ibp3 = 1
        iep3 = 1
        jbp3 = 1
        jep3 = 1
        kbp3 = 1
        kep3 = 1
      endif

      !-------------------------------------------------------------!
      !  Diagnostic arrays and vars:

      ibdt = 1
      iedt = 1
      jbdt = 1
      jedt = 1
      kbdt = 1
      kedt = 1

      ibdq = 1
      iedq = 1
      jbdq = 1
      jedq = 1
      kbdq = 1
      kedq = 1

      ibdv = 1
      iedv = 1
      jbdv = 1
      jedv = 1
      kbdv = 1
      kedv = 1

      ntdiag   = 0
      td_hadv  = 0
      td_vadv  = 0
      td_hturb = 0
      td_vturb = 0
      td_hidiff = 0
      td_vidiff = 0
      td_hediff = 0
      td_vediff = 0
      td_mp    = 0
      td_rdamp = 0
      td_nudge = 0
      td_rad   = 0
      td_div   = 0
      td_diss  = 0
      td_pbl   = 0
      td_lsw   = 0
      td_efall = 0

      td_cond  = 0
      td_evac  = 0
      td_evar  = 0
      td_dep   = 0
      td_subl  = 0
      td_melt  = 0
      td_frz   = 0

      IF( output_thbudget.eq.1 )THEN
        ibdt = ib
        iedt = ie
        jbdt = jb
        jedt = je
        kbdt = kb
        kedt = ke

        ntdiag = ntdiag + 1
        td_hadv = ntdiag

        ntdiag = ntdiag + 1
        td_vadv = ntdiag

        if( cm1setup.ge.1 )then
          if(dohturb)then
            ntdiag = ntdiag + 1
            td_hturb = ntdiag
          endif
          if(dovturb)then
            ntdiag = ntdiag + 1
            td_vturb = ntdiag
          endif
        endif

        if( hadvordrs.eq.3 .or. hadvordrs.eq.5 .or. hadvordrs.eq.7 .or. hadvordrs.eq.9 .or. advwenos.ge.1 )then
          ntdiag = ntdiag + 1
          td_hidiff = ntdiag
        endif

        if( vadvordrs.eq.3 .or. vadvordrs.eq.5 .or. vadvordrs.eq.7 .or. vadvordrs.eq.9 .or. advwenos.ge.1 )then
          ntdiag = ntdiag + 1
          td_vidiff = ntdiag
        endif

        if( idiff.eq.1 )then
          ntdiag = ntdiag + 1
          td_hediff = ntdiag

          if( difforder.eq.2 )then
            ntdiag = ntdiag + 1
            td_vediff = ntdiag
          endif
        endif

        if( imoist.ge.1 .and. ptype.ge.1 )then
          ntdiag = ntdiag + 1
          td_mp = ntdiag
        endif

        if( irdamp.ge.1 )then
          ntdiag = ntdiag + 1
          td_rdamp = ntdiag
        endif

        if( testcase.eq.10 .or. do_lsnudge )then
          ntdiag = ntdiag + 1
          td_nudge = ntdiag
        endif

        if( rterm.eq.1 .or. radopt.ge.1 )then
          ntdiag = ntdiag + 1
          td_rad = ntdiag
        endif

        if( imoist.eq.1 .and. eqtset.eq.2 )then
          ntdiag = ntdiag + 1
          td_div = ntdiag
        endif

        if( idiss.eq.1 .or. ipbl.eq.3 .or. ipbl.eq.4 .or. ipbl.eq.5 )then
          ntdiag = ntdiag + 1
          td_diss = ntdiag
        endif

        if( use_pbl )then
          ntdiag = ntdiag + 1
          td_pbl = ntdiag
        endif

        if( dolsw )then
          ntdiag = ntdiag + 1
          td_lsw = ntdiag
        endif

        if( efall.eq.1 )then
          ntdiag = ntdiag + 1
          td_efall = ntdiag
        endif

        if( imoist.eq.1 .and. ptype.eq.5 )then
          ! Morrison microphysics:
          ntdiag = ntdiag + 1
          td_cond = ntdiag

          ntdiag = ntdiag + 1
          td_evac = ntdiag

          ntdiag = ntdiag + 1
          td_evar = ntdiag

          ntdiag = ntdiag + 1
          td_dep = ntdiag

          ntdiag = ntdiag + 1
          td_subl = ntdiag

          ntdiag = ntdiag + 1
          td_melt = ntdiag

          ntdiag = ntdiag + 1
          td_frz = ntdiag
        endif
      ENDIF

      nqdiag   = 0
      qd_dbz   = 0
      qd_vtc   = 0
      qd_vtr   = 0
      qd_vts   = 0
      qd_vtg   = 0
      qd_vti   = 0
      qd_hadv  = 0
      qd_vadv  = 0
      qd_hturb = 0
      qd_vturb = 0
      qd_hidiff = 0
      qd_vidiff = 0
      qd_hediff = 0
      qd_vediff = 0
      qd_mp    = 0
      qd_nudge = 0
      qd_pbl   = 0
      qd_lsw   = 0

      qd_cond  = 0
      qd_evac  = 0
      qd_evar  = 0
      qd_dep   = 0
      qd_subl  = 0

      ! dbz (new for cm1r19)
      IF( output_dbz.eq.1 .or. restart_file_dbz .or. (imoist.eq.1.and.ptype.eq.3) )THEN
        ibdq = ib
        iedq = ie
        jbdq = jb
        jedq = je
        kbdq = kb
        kedq = ke
        nqdiag = nqdiag+1
        qd_dbz = nqdiag
      ENDIF

      getvt = .false.
      IF( efall.eq.1 ) getvt = .true.
      IF( output_fallvel.eq.1 ) getvt = .true.
      IF( imoist.eq.1 .and. dodomaindiag .and. ( ptype.eq.3 .or. ptype.eq.5 ) ) getvt = .true.
      IF( imoist.eq.1 .and. dodomaindiag .and. testcase.eq.5 ) getvt = .true.

      IF( getvt )THEN
        if( imoist.eq.1 .and. ( ptype.eq.3 .or. ptype.eq.5 ) )then
          ! Thompson and Morrison schemes only (for now):
          ibdq = ib
          iedq = ie
          jbdq = jb
          jedq = je
          kbdq = kb
          kedq = ke

          nqdiag = nqdiag+1
          qd_vtc = nqdiag

          nqdiag = nqdiag+1
          qd_vtr = nqdiag

          nqdiag = nqdiag+1
          qd_vts = nqdiag

          nqdiag = nqdiag+1
          qd_vtg = nqdiag

          nqdiag = nqdiag+1
          qd_vti = nqdiag

        endif
      ENDIF
      IF( output_qvbudget.eq.1 .and. imoist.eq.1 )THEN
        ibdq = ib
        iedq = ie
        jbdq = jb
        jedq = je
        kbdq = kb
        kedq = ke

        nqdiag = nqdiag + 1
        qd_hadv = nqdiag

        nqdiag = nqdiag + 1
        qd_vadv = nqdiag

        if( cm1setup.ge.1 .or. ipbl.ge.1 )then
          if(dohturb)then
            nqdiag = nqdiag + 1
            qd_hturb = nqdiag
          endif
          if(dovturb)then
            nqdiag = nqdiag + 1
            qd_vturb = nqdiag
          endif
        endif

        if( hadvordrs.eq.3 .or. hadvordrs.eq.5 .or. hadvordrs.eq.7 .or. hadvordrs.eq.9 .or. advwenos.ge.1 )then
          nqdiag = nqdiag + 1
          qd_hidiff = nqdiag
        endif

        if( vadvordrs.eq.3 .or. vadvordrs.eq.5 .or. vadvordrs.eq.7 .or. vadvordrs.eq.9 .or. advwenos.ge.1 )then
          nqdiag = nqdiag + 1
          qd_vidiff = nqdiag
        endif

        if( idiff.eq.1 )then
          nqdiag = nqdiag + 1
          qd_hediff = nqdiag

          if( difforder.eq.2 )then
            nqdiag = nqdiag + 1
            qd_vediff = nqdiag
          endif
        endif

        if( ptype.ge.1 )then
          nqdiag = nqdiag + 1
          qd_mp = nqdiag
        endif

        if( testcase.eq.10 .or. do_lsnudge )then
          nqdiag = nqdiag + 1
          qd_nudge = nqdiag
        endif

        if( use_pbl )then
          nqdiag = nqdiag + 1
          qd_pbl = nqdiag
        endif

        if( dolsw )then
          nqdiag = nqdiag + 1
          qd_lsw = nqdiag
        endif

        if( ptype.eq.5 )then
          ! Morrison microphysics:
          nqdiag = nqdiag + 1
          qd_cond = nqdiag

          nqdiag = nqdiag + 1
          qd_evac = nqdiag

          nqdiag = nqdiag + 1
          qd_evar = nqdiag

          nqdiag = nqdiag + 1
          qd_dep = nqdiag

          nqdiag = nqdiag + 1
          qd_subl = nqdiag
        endif
      ENDIF

      nudiag     = 0
      ud_hadv    = 0
      ud_vadv    = 0
      ud_hturb   = 0
      ud_vturb   = 0
      ud_hidiff  = 0
      ud_vidiff  = 0
      ud_hediff  = 0
      ud_vediff  = 0
      ud_pgrad   = 0
      ud_rdamp   = 0
      ud_nudge   = 0
      ud_cor     = 0
      ud_cent    = 0
      ud_pbl     = 0
      ud_lsw     = 0
      nutk       = 0

      IF( output_ubudget.eq.1 )THEN
        ! 170721: to more easily calculate rtke budget, save budget terms

        ibdv = ib
        iedv = ie
        jbdv = jb
        jedv = je
        kbdv = kb
        kedv = ke

        nudiag = nudiag+1
        ud_hadv = nudiag

        nudiag = nudiag+1
        ud_vadv = nudiag

        if( cm1setup.ge.1 .or. ipbl.ge.1 )then
          if(dohturb)then
            nudiag = nudiag+1
            ud_hturb = nudiag
          endif
          if(dovturb)then
            nudiag = nudiag+1
            ud_vturb = nudiag
          endif
        endif

        if( hadvordrv.eq.3 .or. hadvordrv.eq.5 .or. hadvordrv.eq.7 .or. hadvordrv.eq.9 .or. advwenov.ge.1 )then
          nudiag = nudiag + 1
          ud_hidiff = nudiag
        endif

        if( vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1 )then
          nudiag = nudiag + 1
          ud_vidiff = nudiag
        endif

        if( idiff.ge.1 )then
          nudiag = nudiag + 1
          ud_hediff = nudiag

          if( difforder.eq.2 )then
            nudiag = nudiag + 1
            ud_vediff = nudiag
          endif
        endif

        nudiag = nudiag+1
        ud_pgrad = nudiag

        if( irdamp.ge.1 .or. hrdamp.ge.1 )then
          nudiag = nudiag+1
          ud_rdamp = nudiag
        endif

        if( do_lsnudge )then
          nudiag = nudiag+1
          ud_nudge = nudiag
        endif

        if( icor.eq.1 )then
          nudiag = nudiag+1
          ud_cor = nudiag
        endif

        if( axisymm.eq.1 )then
          nudiag = nudiag+1
          ud_cent = nudiag
        endif

        if( use_pbl )then
          nudiag = nudiag + 1
          ud_pbl = nudiag
        endif

        if( dolsw )then
          nudiag = nudiag + 1
          ud_lsw = nudiag
        endif

        if( dodomaindiag )then
          nudiag = nudiag + 1
          nutk = nudiag
        endif

      ENDIF

      nvdiag     = 0
      vd_hadv    = 0
      vd_vadv    = 0
      vd_hturb   = 0
      vd_vturb   = 0
      vd_hidiff  = 0
      vd_vidiff  = 0
      vd_hediff  = 0
      vd_vediff  = 0
      vd_pgrad   = 0
      vd_rdamp   = 0
      vd_nudge   = 0
      vd_cor     = 0
      vd_cent    = 0
      vd_pbl     = 0
      vd_lsw     = 0
      nvtk       = 0

      IF( output_vbudget.eq.1 )THEN
        ! 170721: to more easily calculate rtke budget, save budget terms

        ibdv = ib
        iedv = ie
        jbdv = jb
        jedv = je
        kbdv = kb
        kedv = ke

        nvdiag = nvdiag+1
        vd_hadv = nvdiag

        nvdiag = nvdiag+1
        vd_vadv = nvdiag

        if( cm1setup.ge.1 .or. ipbl.ge.1 )then
          if(dohturb)then
            nvdiag = nvdiag+1
            vd_hturb = nvdiag
          endif
          if(dovturb)then
            nvdiag = nvdiag+1
            vd_vturb = nvdiag
          endif
        endif

        if( hadvordrv.eq.3 .or. hadvordrv.eq.5 .or. hadvordrv.eq.7 .or. hadvordrv.eq.9 .or. advwenov.ge.1 )then
          nvdiag = nvdiag + 1
          vd_hidiff = nvdiag
        endif

        if( vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1 )then
          nvdiag = nvdiag + 1
          vd_vidiff = nvdiag
        endif

        if( idiff.ge.1 )then
          nvdiag = nvdiag + 1
          vd_hediff = nvdiag

          if( difforder.eq.2 )then
            nvdiag = nvdiag + 1
            vd_vediff = nvdiag
          endif
        endif

        nvdiag = nvdiag+1
        vd_pgrad = nvdiag

        if( irdamp.ge.1 .or. hrdamp.ge.1 )then
          nvdiag = nvdiag+1
          vd_rdamp = nvdiag
        endif

        if( do_lsnudge )then
          nvdiag = nvdiag+1
          vd_nudge = nvdiag
        endif

        if( icor.eq.1 )then
          nvdiag = nvdiag+1
          vd_cor = nvdiag
        endif

        if( axisymm.eq.1 )then
          nvdiag = nvdiag+1
          vd_cent = nvdiag
        endif

        if( use_pbl )then
          nvdiag = nvdiag + 1
          vd_pbl = nvdiag
        endif

        if( dolsw )then
          nvdiag = nvdiag + 1
          vd_lsw = nvdiag
        endif

        if( dodomaindiag )then
          nvdiag = nvdiag + 1
          nvtk = nvdiag
        endif

      ENDIF


      nwdiag     = 0
      wd_hadv    = 0
      wd_vadv    = 0
      wd_hturb   = 0
      wd_vturb   = 0
      wd_hidiff  = 0
      wd_vidiff  = 0
      wd_hediff  = 0
      wd_vediff  = 0
      wd_pgrad   = 0
      wd_rdamp   = 0
      wd_buoy    = 0
      nwtk       = 0

      IF( output_wbudget.eq.1 )THEN
        ! 170721: to more easily calculate rtke budget, save budget terms

        ibdv = ib
        iedv = ie
        jbdv = jb
        jedv = je
        kbdv = kb
        kedv = ke

        nwdiag = nwdiag+1
        wd_hadv = nwdiag

        nwdiag = nwdiag+1
        wd_vadv = nwdiag

        if( cm1setup.ge.1 .or. ipbl.ge.1 )then
          if(dohturb)then
            nwdiag = nwdiag+1
            wd_hturb = nwdiag
          endif
          if(dovturb)then
            nwdiag = nwdiag+1
            wd_vturb = nwdiag
          endif
        endif

        if( hadvordrv.eq.3 .or. hadvordrv.eq.5 .or. hadvordrv.eq.7 .or. hadvordrv.eq.9 .or. advwenov.ge.1 )then
          nwdiag = nwdiag + 1
          wd_hidiff = nwdiag
        endif

        if( vadvordrv.eq.3 .or. vadvordrv.eq.5 .or. vadvordrv.eq.7 .or. vadvordrv.eq.9 .or. advwenov.ge.1 )then
          nwdiag = nwdiag + 1
          wd_vidiff = nwdiag
        endif

        if( idiff.ge.1 )then
          nwdiag = nwdiag + 1
          wd_hediff = nwdiag

          if( difforder.eq.2 )then
            nwdiag = nwdiag + 1
            wd_vediff = nwdiag
          endif
        endif

        nwdiag = nwdiag+1
        wd_pgrad = nwdiag

        if( irdamp.ge.1 .or. hrdamp.ge.1 )then
          nwdiag = nwdiag+1
          wd_rdamp = nwdiag
        endif

        nwdiag = nwdiag+1
        wd_buoy = nwdiag

        if( dodomaindiag )then
          nwdiag = nwdiag + 1
          nwtk = nwdiag
        endif

      ENDIF

      !-----

      ! subgrid tke:
      ibdk = 1
      iedk = 1
      jbdk = 1
      jedk = 1
      kbdk = 1
      kedk = 1

      nkdiag   = 0
      kd_adv   = 0
      kd_turb  = 0

      ! cm1r21.1: add idoles
      IF( (idopbl.or.idoles) .and. dodomaindiag .and. iusetke )THEN

        ibdk = ib
        iedk = ie
        jbdk = jb
        jedk = je
        kbdk = kb
        kedk = ke

        nkdiag = nkdiag + 1
        kd_adv = nkdiag

        nkdiag = nkdiag + 1
        kd_turb = nkdiag

      ENDIF

      !-----

      ! for nondim pressure:
      ibdp = 1
      iedp = 1
      jbdp = 1
      jedp = 1
      kbdp = 1
      kedp = 1

      npdiag   = 0

      IF( pdcomp )THEN

        ibdp = ib
        iedp = ie
        jbdp = jb
        jedp = je
        kbdp = kb
        kedp = ke

        npdiag = 3

        if( icor.eq.1 ) npdiag = npdiag+1

      ENDIF

      !-----

      ntdiag = max( ntdiag , 1 )
      nqdiag = max( nqdiag , 1 )
      nudiag = max( nudiag , 1 )
      nvdiag = max( nvdiag , 1 )
      nwdiag = max( nwdiag , 1 )
      nkdiag = max( nkdiag , 1 )
      npdiag = max( npdiag , 1 )

      !-------------------------------------------------------------!

      ib2d = 1
      ie2d = 1
      jb2d = 1
      je2d = 1

      IF( nout2d .gt. 0 )THEN
        ib2d = ib
        ie2d = ie
        jb2d = jb
        je2d = je
      ENDIF

      nout2d = max( nout2d , 1 )

      !-------------------------------------------------------------!

      ib3d = 1
      ie3d = 1
      jb3d = 1
      je3d = 1
      kb3d = 1
      ke3d = 1

      IF( nout3d .gt. 0 )THEN
        ib3d = ib
        ie3d = ie
        jb3d = jb
        je3d = je
        kb3d = kb
        ke3d = ke
      ENDIF

      nout3d = max( nout3d , 1 )

      !-------------------------------------------------------------!

#ifdef MPI
      mynode = int( float(myid)/float(ppnode) )
      nodeleader = mynode * ppnode
      nodes = nodex * nodey / ppnode
#else
      mynode = 0
      nodeleader = 0
      nodes = 1
      ppnode = 1
#endif

      d2i  = 1
      d2is = 1
      d2iu = 1
      d2iv = 1

      d2j  = 1
      d2js = 1
      d2ju = 1
      d2jv = 1

      d3i  = 1
      d3is = 1
      d3iu = 1
      d3iv = 1

      d3j  = 1
      d3js = 1
      d3ju = 1
      d3jv = 1

      d3n = 1
      d3t = 1

    IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN

      IF( myid.eq.nodeleader )THEN

        d2i   = (ni+1)*nodex
        d2is  = nx
        d2iu  = nx+1
        d2iv  = nx

        d2j   = (nj+1)*nodey
        d2js  = ny
        d2ju  = ny
        d2jv  = ny+1

        d3n = max( numprocs , 1 )

        d3t = max( ppnode-1 + nodes-1 , 1 )

      ENDIF

    ELSEIF( output_filetype.eq.3 )THEN

        d2i   = ni+1
        d2is  = ni
        d2iu  = ni+1
        d2iv  = ni

        d2j   = nj+1
        d2js  = nj
        d2ju  = nj
        d2jv  = nj+1

        d3n = 1

        d3t = 1

    ELSE

        if( myid.eq.0 )then
        print *
        print *,'  output_filetype = ',output_filetype
        print *
        print *,'  output_filetype must be 1,2,3 '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1

    ENDIF

    IF( rstfrq.gt.0.0 .and. rstfrq.le.timax )THEN
    IF( restart_filetype.eq.1 .or. restart_filetype.eq.2 .or. restart_filetype.eq.3 )THEN

      IF( myid.eq.nodeleader )THEN

        d2i   = (ni+1)*nodex
        d2is  = nx
        d2iu  = nx+1
        d2iv  = nx

        d2j   = (nj+1)*nodey
        d2js  = ny
        d2ju  = ny
        d2jv  = ny+1

        d3n = max( numprocs , 1 )

        d3t = max( ppnode-1 + nodes-1 , 1 )

      ENDIF

    ENDIF
    ENDIF


        d3i   = nimax+1
        d3is  = nimax
        d3iu  = nimax+1
        d3iv  = nimax

        d3j   = njmax+1
        d3js  = njmax
        d3ju  = njmax
        d3jv  = njmax+1


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  myid,mynode             = ',myid,mynode
      if(dowr) write(outfile,*) '  nodeleader,nodes,ppnode = ',nodeleader,nodes,ppnode
      if(dowr) write(outfile,*) '  nodex,nodey             = ',nodex,nodey
      if(dowr) write(outfile,*) '  nimax,njmax             = ',nimax,njmax
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d2i   = ',d2i
      if(dowr) write(outfile,*) '  d2is  = ',d2is
      if(dowr) write(outfile,*) '  d2iu  = ',d2iu
      if(dowr) write(outfile,*) '  d2iv  = ',d2iv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d2j   = ',d2j
      if(dowr) write(outfile,*) '  d2js  = ',d2js
      if(dowr) write(outfile,*) '  d2ju  = ',d2ju
      if(dowr) write(outfile,*) '  d2jv  = ',d2jv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d3i   = ',d3i
      if(dowr) write(outfile,*) '  d3is  = ',d3is
      if(dowr) write(outfile,*) '  d3iu  = ',d3iu
      if(dowr) write(outfile,*) '  d3iv  = ',d3iv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d3j   = ',d3j
      if(dowr) write(outfile,*) '  d3js  = ',d3js
      if(dowr) write(outfile,*) '  d3ju  = ',d3ju
      if(dowr) write(outfile,*) '  d3jv  = ',d3jv
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  d3n   = ',d3n
      if(dowr) write(outfile,*) '  d3t   = ',d3t
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  cm1r19:  for thrd array:

      ibb2 = 1
      ibe2 = 1
      jbb2 = 1
      jbe2 = 1
      kbb2 = 1
      kbe2 = 1

!!!      ibb2 = ib
!!!      ibe2 = ie
!!!      jbb2 = jb
!!!      jbe2 = je
!!!      kbb2 = kb
!!!      kbe2 = ke
!!!
!!!      if(dowr) write(outfile,*) '  ibb2  = ',ibb2
!!!      if(dowr) write(outfile,*) '  ibe2  = ',ibe2
!!!      if(dowr) write(outfile,*) '  jbb2  = ',jbb2
!!!      if(dowr) write(outfile,*) '  jbe2  = ',jbe2
!!!      if(dowr) write(outfile,*) '  kbb2  = ',kbb2
!!!      if(dowr) write(outfile,*) '  kbe2  = ',kbe2
!!!      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      stat_w       = max(0,min(1,stat_w))
      stat_wlevs   = max(0,min(1,stat_wlevs))
      stat_u       = max(0,min(1,stat_u))
      stat_v       = max(0,min(1,stat_v))
      stat_rmw     = max(0,min(1,stat_rmw))
      IF(axisymm.ne.1) stat_rmw = 0
      stat_pipert  = max(0,min(1,stat_pipert))
      stat_prspert = max(0,min(1,stat_prspert))
      stat_thpert  = max(0,min(1,stat_thpert))
      stat_q       = max(0,min(1,stat_q))
      stat_tke     = max(0,min(1,stat_tke))
      stat_km      = max(0,min(1,stat_km))
      stat_kh      = max(0,min(1,stat_kh))
      stat_div     = max(0,min(1,stat_div))
      stat_rh      = max(0,min(1,stat_rh))
      stat_rhi     = max(0,min(1,stat_rhi))
      stat_the     = max(0,min(1,stat_the))
      stat_cloud   = max(0,min(1,stat_cloud))
      stat_sfcprs  = max(0,min(1,stat_sfcprs))
      stat_wsp     = max(0,min(1,stat_wsp))
      stat_cfl     = max(0,min(1,stat_cfl))
      stat_vort    = max(0,min(1,stat_vort))
      stat_tmass   = max(0,min(1,stat_tmass))
      stat_tmois   = max(0,min(1,stat_tmois))
      stat_qmass   = max(0,min(1,stat_qmass))
      stat_tenerg  = max(0,min(1,stat_tenerg))
      stat_mo      = max(0,min(1,stat_mo))
      stat_tmf     = max(0,min(1,stat_tmf))
      stat_pcn     = max(0,min(1,stat_pcn))
      stat_qsrc    = max(0,min(1,stat_qsrc))


      if(imoist.eq.0)then
        stat_q=0
        stat_rh=0
        stat_rhi=0
        stat_the=0
        stat_cloud=0
        stat_tmois=0
        stat_qmass=0
        stat_pcn=0
        stat_qsrc=0
      endif
      if(iice.eq.0)then
        stat_rhi=0
      endif

      if( iusekm.eq.0 )then
        stat_km = 0
      endif
      if( iusekh.eq.0 )then
        stat_kh = 0
      endif
      if( .not. iusetke )then
        stat_tke = 0
      endif


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stat_w       = ',stat_w
      if(dowr) write(outfile,*) 'stat_wlevs   = ',stat_wlevs
      if(dowr) write(outfile,*) 'stat_u       = ',stat_u
      if(dowr) write(outfile,*) 'stat_v       = ',stat_v
      if(dowr) write(outfile,*) 'stat_rmw     = ',stat_rmw
      if(dowr) write(outfile,*) 'stat_pipert  = ',stat_pipert
      if(dowr) write(outfile,*) 'stat_prspert = ',stat_prspert
      if(dowr) write(outfile,*) 'stat_thpert  = ',stat_thpert
      if(dowr) write(outfile,*) 'stat_q       = ',stat_q
      if(dowr) write(outfile,*) 'stat_tke     = ',stat_tke
      if(dowr) write(outfile,*) 'stat_km      = ',stat_km
      if(dowr) write(outfile,*) 'stat_kh      = ',stat_kh
      if(dowr) write(outfile,*) 'stat_div     = ',stat_div
      if(dowr) write(outfile,*) 'stat_rh      = ',stat_rh
      if(dowr) write(outfile,*) 'stat_rhi     = ',stat_rhi
      if(dowr) write(outfile,*) 'stat_the     = ',stat_the
      if(dowr) write(outfile,*) 'stat_cloud   = ',stat_cloud
      if(dowr) write(outfile,*) 'stat_sfcprs  = ',stat_sfcprs
      if(dowr) write(outfile,*) 'stat_wsp     = ',stat_wsp
      if(dowr) write(outfile,*) 'stat_cfl     = ',stat_cfl
      if(dowr) write(outfile,*) 'stat_vort    = ',stat_vort
      if(dowr) write(outfile,*) 'stat_tmass   = ',stat_tmass
      if(dowr) write(outfile,*) 'stat_tmois   = ',stat_tmois
      if(dowr) write(outfile,*) 'stat_qmass   = ',stat_qmass
      if(dowr) write(outfile,*) 'stat_tenerg  = ',stat_tenerg
      if(dowr) write(outfile,*) 'stat_mo      = ',stat_mo
      if(dowr) write(outfile,*) 'stat_tmf     = ',stat_tmf
      if(dowr) write(outfile,*) 'stat_pcn     = ',stat_pcn
      if(dowr) write(outfile,*) 'stat_qsrc    = ',stat_qsrc
      if(dowr) write(outfile,*)


      prcl_th       = max(0,min(1,prcl_th))
      prcl_t        = max(0,min(1,prcl_t))
      prcl_prs      = max(0,min(1,prcl_prs))
      prcl_ptra     = max(0,min(1,prcl_ptra))
      prcl_q        = max(0,min(1,prcl_q))
      prcl_nc       = max(0,min(1,prcl_nc))
      prcl_km       = max(0,min(1,prcl_km))
      prcl_kh       = max(0,min(1,prcl_kh))
      prcl_tke      = max(0,min(1,prcl_tke))
      prcl_dbz      = max(0,min(1,prcl_dbz))
      prcl_b        = max(0,min(1,prcl_b))
      prcl_vpg      = max(0,min(1,prcl_vpg))
      prcl_vort     = max(0,min(1,prcl_vort))
      prcl_rho      = max(0,min(1,prcl_rho))
      prcl_qsat     = max(0,min(1,prcl_qsat))
      prcl_sfc      = max(0,min(1,prcl_sfc))


      if( iptra.eq.0 )then
        prcl_ptra = 0
      endif
      if( imoist.eq.0 )then
        prcl_q  = 0
        prcl_nc = 0
        prcl_dbz = 0
        prcl_qsat = 0
      else
        if( ptype.eq.1 .or. ptype.eq.4 )then
          prcl_dbz = 0
        endif
        if( idm.eq.0 )then
          prcl_nc = 0
        endif
      endif
      if( sgsmodel.le.0 .or. ipbl.le.0 )then
        prcl_km = 0
        prcl_kh = 0
      endif
      if( .not. iusetke )then
        prcl_tke = 0
      endif
      if( bbc.ne.3 )then
        prcl_sfc = 0
      endif


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'prcl_th      = ',prcl_th
      if(dowr) write(outfile,*) 'prcl_t       = ',prcl_t
      if(dowr) write(outfile,*) 'prcl_prs     = ',prcl_prs
      if(dowr) write(outfile,*) 'prcl_ptra    = ',prcl_ptra
      if(dowr) write(outfile,*) 'prcl_q       = ',prcl_q
      if(dowr) write(outfile,*) 'prcl_nc      = ',prcl_nc
      if(dowr) write(outfile,*) 'prcl_km      = ',prcl_km
      if(dowr) write(outfile,*) 'prcl_kh      = ',prcl_kh
      if(dowr) write(outfile,*) 'prcl_tke     = ',prcl_tke
      if(dowr) write(outfile,*) 'prcl_dbz     = ',prcl_dbz
      if(dowr) write(outfile,*) 'prcl_b       = ',prcl_b
      if(dowr) write(outfile,*) 'prcl_vpg     = ',prcl_vpg
      if(dowr) write(outfile,*) 'prcl_vort    = ',prcl_vort
      if(dowr) write(outfile,*) 'prcl_rho     = ',prcl_rho
      if(dowr) write(outfile,*) 'prcl_qsat    = ',prcl_qsat
      if(dowr) write(outfile,*) 'prcl_sfc     = ',prcl_sfc
      if(dowr) write(outfile,*)

!------------------------------------------------------------------

    IF( dodomaindiag )THEN
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'dodomaindiag = ',dodomaindiag
      if(dowr) write(outfile,*) 'diagfrq      = ',diagfrq
      if(dowr) write(outfile,*)
    ENDIF

!------------------------------------------------------------------
!  "constants" in subgrid turbulence schemes:

      IF( dowr .and. ( sgsmodel.ge.1 ) )THEN
        write(outfile,*)
        write(outfile,*) '  c_m,c_l,ri_c  = ',c_m,c_l,ri_c
        write(outfile,*) '  c_e1,c_e2,sum = ',c_e1,c_e2,c_e1+c_e2
        write(outfile,*) '  c_s,rcs       = ',c_s,rcs
      ENDIF

!--------------------------------------------------------------

      rdx=1.0/dx
      rdy=1.0/dy
      rdz=1.0/dz
      rdx2=1.0/(2.0*dx)
      rdy2=1.0/(2.0*dy)
      rdz2=1.0/(2.0*dz)
      rdx4=1.0/(4.0*dx)
      rdy4=1.0/(4.0*dy)
      rdz4=1.0/(4.0*dz)

      thec_mb=0.0
      qt_mb=0.0

      rsttim=rstfrq

      ! do these at initial time (t=0)
      stattim=0.0
      taptim=0.0
      prcltim=0.0
      radtim=0.0

      if( iprcl.le.0 ) prcltim = 1.0d60
      if( radopt.le.0 ) radtim = 1.0d60

!--------------------------------------------------------------

      npvals = 1

      prx = 0
      pry = 0
      prz = 0
      pru = 0
      prv = 0
      prw = 0
      prtime = 0
      prth = 0
      prt = 0
      prprs = 0
      prpt1 = 0
      prpt2 = 0
      prqv = 0
      prq1 = 0
      prq2 = 0
      prnc1 = 0
      prnc2 = 0
      prkm = 0
      prkh = 0
      prtke = 0
      prdbz = 0
      prb = 0
      prvpg = 0
      przv = 0
      prrho = 0
      prqsl = 0
      prqsi = 0
      prznt = 0
      prust = 0
      przs = 0
      prsig = 0

      ! for parcels:
      if(iprcl.eq.1)then

        ! 7 basic variables for all simulations:
        ! (x,y,z,u,v,w,t)
        !  1 2 3 4 5 6 7
        npvals = 7

        prx = 1
        pry = 2
        prz = 3
        pru = 4
        prv = 5
        prw = 6
        prtime = 7

        if( prcl_th.eq.1 )then
          npvals = npvals+1
          prth = npvals
        endif
        if( prcl_t.eq.1 )then
          npvals = npvals+1
          prt = npvals
        endif
        if( prcl_prs.eq.1 )then
          npvals = npvals+1
          prprs = npvals
        endif

        ! passive tracers:
        if( prcl_ptra.eq.1 )then
          prpt1 = npvals+1
          npvals = npvals + npt
          prpt2 = npvals
        endif

        ! moisture variables:
        if( prcl_q.eq.1 )then
          !---
          npvals = npvals+1
          prqv = npvals
          prq1 = npvals+1
          !---
          npvals = npvals+(nql2-nql1+1)
          if( iice.eq.1 )then
            npvals = npvals+(nqs2-nqs1+1)
          endif
          !---
          prq2 = npvals
        endif
        if( prcl_nc.eq.1 .and. idm.eq.1 )then
          prnc1 = npvals+1
          npvals = npvals+(nnc2-nnc1+1)
          prnc2 = npvals
        endif


        ! turbulence parameters:
        if( cm1setup.ge.1 .or. ipbl.ge.1 )then
          if( prcl_km.eq.1 )then
            prkm = npvals+1
            npvals = npvals+2
          endif
          if( prcl_kh.eq.1 )then
            prkh = npvals+1
            npvals = npvals+2
          endif
          if( prcl_tke.eq.1 )then
            npvals = npvals+1
            prtke = npvals
          endif
        endif

        if( prcl_dbz.eq.1 )then
          npvals = npvals+1
          prdbz = npvals
        endif
        if( prcl_b.eq.1 )then
          npvals = npvals+1
          prb = npvals
        endif
        if( prcl_vpg.eq.1 )then
          npvals = npvals+1
          prvpg = npvals
        endif
        if( prcl_vort.eq.1 )then
        if( .not. terrain_flag )then
          npvals = npvals+1
          przv = npvals
        endif
        endif
        if( prcl_rho.eq.1 )then
          npvals = npvals+1
          prrho = npvals
        endif
        if( prcl_qsat.eq.1 )then
          npvals = npvals+1
          prqsl = npvals
          if( iice.eq.1 )then
            npvals = npvals+1
            prqsi = npvals
          endif
        endif
        if( prcl_sfc.eq.1 )then
          npvals = npvals+1
          prznt = npvals
          npvals = npvals+1
          prust = npvals
        endif
        if( terrain_flag )then
          npvals = npvals+1
          przs = npvals
          npvals = npvals+1
          prsig = npvals
        endif

      else

        nparcels = 1

      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  nparcels = ',nparcels
      if(dowr) write(outfile,*) '  npvals   = ',npvals
      if(dowr) write(outfile,*) '    prx      = ',prx
      if(dowr) write(outfile,*) '    pry      = ',pry
      if(dowr) write(outfile,*) '    prz      = ',prz
      if(dowr) write(outfile,*) '    pru      = ',pru
      if(dowr) write(outfile,*) '    prv      = ',prv
      if(dowr) write(outfile,*) '    prw      = ',prw
      if(dowr) write(outfile,*) '    prtime   = ',prtime
      if(dowr) write(outfile,*) '    prth     = ',prth
      if(dowr) write(outfile,*) '    prt      = ',prt
      if(dowr) write(outfile,*) '    prprs    = ',prprs
      if(dowr) write(outfile,*) '    prpt1    = ',prpt1
      if(dowr) write(outfile,*) '    prpt2    = ',prpt2
      if(dowr) write(outfile,*) '    prqv     = ',prqv
      if(dowr) write(outfile,*) '    prq1     = ',prq1
      if(dowr) write(outfile,*) '    prq2     = ',prq2
      if(dowr) write(outfile,*) '    prnc1    = ',prnc1
      if(dowr) write(outfile,*) '    prnc2    = ',prnc2
      if(dowr) write(outfile,*) '    prkm     = ',prkm
      if(dowr) write(outfile,*) '    prkh     = ',prkh
      if(dowr) write(outfile,*) '    prtke    = ',prtke
      if(dowr) write(outfile,*) '    prdbz    = ',prdbz
      if(dowr) write(outfile,*) '    prb      = ',prb
      if(dowr) write(outfile,*) '    prvpg    = ',prvpg
      if(dowr) write(outfile,*) '    przv     = ',przv
      if(dowr) write(outfile,*) '    prrho    = ',prrho
      if(dowr) write(outfile,*) '    prqsl    = ',prqsl
      if(dowr) write(outfile,*) '    prqsi    = ',prqsi
      if(dowr) write(outfile,*) '    prznt    = ',prznt
      if(dowr) write(outfile,*) '    prust    = ',prust
      if(dowr) write(outfile,*) '    przs     = ',przs
      if(dowr) write(outfile,*) '    prsig    = ',prsig
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Get identity

      ibw=0
      ibe=0
      ibs=0
      ibn=0

      patchsws = .false.
      patchsww = .false.
      patchses = .false.
      patchsee = .false.
      patchnwn = .false.
      patchnww = .false.
      patchnen = .false.
      patchnee = .false.

      p2tchsws = .false.
      p2tchsww = .false.
      p2tchses = .false.
      p2tchsee = .false.
      p2tchnwn = .false.
      p2tchnww = .false.
      p2tchnen = .false.
      p2tchnee = .false.

#ifdef MPI

      if(dowr) write(outfile,*) '  myi,myj=',myi,myj
      if(dowr) write(outfile,*)

      mynorth = nabor(myi,   myj+1, nodex, nodey)
      mysouth = nabor(myi,   myj-1, nodex, nodey)
      myeast  = nabor(myi+1, myj,   nodex, nodey)
      mywest  = nabor(myi-1, myj,   nodex, nodey)

      if(dowr) write(outfile,*) '  mywest  =',mywest
      if(dowr) write(outfile,*) '  myeast  =',myeast
      if(dowr) write(outfile,*) '  mysouth =',mysouth
      if(dowr) write(outfile,*) '  mynorth =',mynorth
      if(dowr) write(outfile,*)

      mysw = nabor(myi-1, myj-1,   nodex, nodey)
      mynw = nabor(myi-1, myj+1,   nodex, nodey)
      myne = nabor(myi+1, myj+1,   nodex, nodey)
      myse = nabor(myi+1, myj-1,   nodex, nodey)

      if(dowr) write(outfile,*) '  mysw    =',mysw
      if(dowr) write(outfile,*) '  mynw    =',mynw
      if(dowr) write(outfile,*) '  myne    =',myne
      if(dowr) write(outfile,*) '  myse    =',myse

      cs1we = (nj)*(nk)
      cs1sn = (ni)*(nk)
      ct1we = (nj)*(nk+1)
      ct1sn = (ni)*(nk+1)
      cv1we = (nj+1)*(nk)
      cu1sn = (ni+1)*(nk)
      cw1we = (nj)*(nk-1)
      cw1sn = (ni)*(nk-1)
      cs2we = 2*(nj)*(nk)
      cs2sn = (ni)*2*(nk)
      cs3we = cmp*(nj)*(nk)
      cs3sn = (ni)*cmp*(nk)
      ct3we = cmp*(nj)*(nk+1)
      ct3sn = (ni)*cmp*(nk+1)
      cv3we = cmp*(nj+1)*(nk)
      cu3sn = cmp*(ni+1)*(nk)
      cw3we = cmp*(nj)*(nk-1)
      cw3sn = (ni)*cmp*(nk-1)
      cs3weq = cmp*(nj)*(nk)*numq
      cs3snq = (ni)*cmp*(nk)*numq

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  cs1we   =',cs1we
      if(dowr) write(outfile,*) '  cs1sn   =',cs1sn
      if(dowr) write(outfile,*) '  ct1we   =',ct1we
      if(dowr) write(outfile,*) '  ct1sn   =',ct1sn
      if(dowr) write(outfile,*) '  cv1we   =',cv1we
      if(dowr) write(outfile,*) '  cu1sn   =',cu1sn
      if(dowr) write(outfile,*) '  cw1we   =',cw1we
      if(dowr) write(outfile,*) '  cw1sn   =',cw1sn
      if(dowr) write(outfile,*) '  cs2we   =',cs2we
      if(dowr) write(outfile,*) '  cs2sn   =',cs2sn
      if(dowr) write(outfile,*) '  cs3we   =',cs3we
      if(dowr) write(outfile,*) '  cs3sn   =',cs3sn
      if(dowr) write(outfile,*) '  ct3we   =',ct3we
      if(dowr) write(outfile,*) '  ct3sn   =',ct3sn
      if(dowr) write(outfile,*) '  cv3we   =',cv3we
      if(dowr) write(outfile,*) '  cu3sn   =',cu3sn
      if(dowr) write(outfile,*) '  cw3we   =',cw3we
      if(dowr) write(outfile,*) '  cw3sn   =',cw3sn
      if(dowr) write(outfile,*) '  cs3weq  =',cs3weq
      if(dowr) write(outfile,*) '  cs3snq  =',cs3snq
      if(dowr) write(outfile,*)
#endif

  ! define grid, assuming no stretching 
  ! (we will correct the grid later if stretch_x or stretch_y >= 1)

    allocate( xfdp(1-ngxy:nx+ngxy+1) )
    allocate( yfdp(1-ngxy:ny+ngxy+1) )

    IF(iorigin.eq.1)THEN

      do i=1-ngxy,nx+ngxy+1
        xfdp(i)=dble(dx)*(i-1)
      enddo

      do j=1-ngxy,ny+ngxy+1
        yfdp(j)=dble(dy)*(j-1)
      enddo

    ELSEIF(iorigin.eq.2)THEN

      do i=1-ngxy,nx+ngxy+1
        xfdp(i)=dble(dx)*(i-1)-0.5d0*dble(dx)*nx
      enddo

      do j=1-ngxy,ny+ngxy+1
        yfdp(j)=dble(dy)*(j-1)-0.5d0*dble(dy)*ny
      enddo

    ELSE

      print *,'  invalid option for iorigin'
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      call stopcm1

    ENDIF


    ! cm1r20:
    if( wbc.ne.1 .and. myi1.eq. 1 ) ibw=1
    if( ebc.ne.1 .and. myi2.eq.nx ) ibe=1
    if( sbc.ne.1 .and. myj1.eq. 1 ) ibs=1
    if( nbc.ne.1 .and. myj2.eq.ny ) ibn=1

      do i=ib,ie+1
        xf(i)=xfdp(i+myi1-1)
      enddo
      do i=ib,ie
        xh(i)=0.5d0*(xfdp(i+myi1-1)+xfdp(i+myi1))
      enddo

      do j=jb,je+1
        yf(j)=yfdp(j+myj1-1)
      enddo
      do j=jb,je
        yh(j)=0.5d0*(yfdp(j+myj1-1)+yfdp(j+myj1))
      enddo


#ifndef MPI
      if(wbc.eq.2)then
        ibw=1
      endif

      if(ebc.eq.2)then
        ibe=1
      endif

      if(sbc.eq.2)then
        ibs=1
      endif

      if(nbc.eq.2)then
        ibn=1
      endif
#endif

!--------------------------------------------------------------

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'g     =',g
      if(dowr) write(outfile,*) 'to    =',to
      if(dowr) write(outfile,*) 'rd    =',rd
      if(dowr) write(outfile,*) 'rv    =',rv
      if(dowr) write(outfile,*) 'cp    =',cp
      if(dowr) write(outfile,*) 'cv    =',cv
      if(dowr) write(outfile,*) 'cpv   =',cpv
      if(dowr) write(outfile,*) 'cvv   =',cvv
      if(dowr) write(outfile,*) 'p00   =',p00
      if(dowr) write(outfile,*) 'rp00  =',rp00
      if(dowr) write(outfile,*) 'th0r  =',th0r
      if(dowr) write(outfile,*) 'rcp   =',rcp
      if(dowr) write(outfile,*) 'pi    =',pi

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'cpdcv =',cpdcv
      if(dowr) write(outfile,*) 'cvdcp =',cvdcp
      if(dowr) write(outfile,*) 'rovcp =',rovcp
      if(dowr) write(outfile,*) 'rddcv =',rddcv
      if(dowr) write(outfile,*) 'cvdrd =',cvdrd
      if(dowr) write(outfile,*) 'cpdrd =',cpdrd
      if(dowr) write(outfile,*) 'eps   =',eps
      if(dowr) write(outfile,*) 'reps  =',reps
      if(dowr) write(outfile,*) 'repsm1=',repsm1
      if(dowr) write(outfile,*) 'cpt   =',cpt
      if(dowr) write(outfile,*) 'cvt   =',cvt
      if(dowr) write(outfile,*) 'pnum  =',pnum
      if(dowr) write(outfile,*) 'xlv   =',xlv
      if(dowr) write(outfile,*) 'xls   =',xls
      if(dowr) write(outfile,*) 'xlf   =',xlf
      if(dowr) write(outfile,*) 'lvdcp =',lvdcp
      if(dowr) write(outfile,*) 'condc =',condc
      if(dowr) write(outfile,*) 'cpl   =',cpl
      if(dowr) write(outfile,*) 'cpi   =',cpi
      if(dowr) write(outfile,*) 'lv1   =',lv1
      if(dowr) write(outfile,*) 'lv2   =',lv2
      if(dowr) write(outfile,*) 'ls1   =',ls1
      if(dowr) write(outfile,*) 'ls2   =',ls2
      if(dowr) write(outfile,*) 'karman=',karman

      if( psolver.eq.6 .or. psolver.eq.7 )then
        !----------------
        ! speed of sound for compressible-Boussinesq equations (psolver=6)
        ! -or- speed of sound for modified compressible equations (psolver=7)

        if( csound.lt.0.1 )then
          csound = 300.0
        endif

!        if( csound.lt.1.0 )then
!        if(myid.eq.0)then
!        print *
!        print *
!        print *,'  psolver  =  ',psolver
!        print *,'  csound   =  ',csound
!        print *
!        print *,'  csound is too small '
!        print *
!        print *,'   stopping model .... '
!        print *
!        print *
!        endif
!#ifdef MPI
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!#endif
!        call stopcm1
!        endif

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) 'csound=',csound
        if(dowr) write(outfile,*)

      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'timeformat   =',timeformat
      if(dowr) write(outfile,*) 'timestats    =',timestats
      if(dowr) write(outfile,*) 'terrain_flag =',terrain_flag
      if(dowr) write(outfile,*) 'procfiles    =',procfiles
      if(dowr) write(outfile,*) 'outunits     =',outunits
      if(dowr) write(outfile,*) 'outunitconv  =',outunitconv


      if(dowr) write(outfile,*)
 
      if(dowr) write(outfile,130) 'ib,ibm,ibi,ibc,ibt=',ib,ibm,ibi,ibc,ibt
      if(dowr) write(outfile,130) 'ie,iem,iei,iec,iet=',ie,iem,iei,iec,iet
      if(dowr) write(outfile,130) 'jb,jbm,jbi,jbc,jbt=',jb,jbm,jbi,jbc,jbt
      if(dowr) write(outfile,130) 'je,jem,jei,jec,jet=',je,jem,jei,jec,jet
      if(dowr) write(outfile,130) 'kb,kbm,kbi,kbc,kbt=',kb,kbm,kbi,kbc,kbt
      if(dowr) write(outfile,130) 'ke,kem,kei,kec,ket=',ke,kem,kei,kec,ket


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  idoles    = ',idoles
      if(dowr) write(outfile,*) '  idopbl    = ',idopbl
      if(dowr) write(outfile,*) '  iusetke   = ',iusetke
      if(dowr) write(outfile,*) '  iusekm    = ',iusekm
      if(dowr) write(outfile,*) '  iusekh    = ',iusekh
      if(dowr) write(outfile,*)

130   format(1x,a19,5(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'imirror,jmirror = ',imirror,jmirror
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,131) 'ibp,itb,ipb,ibr,ibb=',ibp,itb,ipb,ibr,ibb
      if(dowr) write(outfile,131) 'iep,ite,ipe,ier,ieb=',iep,ite,ipe,ier,ieb
      if(dowr) write(outfile,131) 'jbp,jtb,jpb,jbr,jbb=',jbp,jtb,jpb,jbr,jbb
      if(dowr) write(outfile,131) 'jep,jte,jpe,jer,jeb=',jep,jte,jpe,jer,jeb
      if(dowr) write(outfile,131) 'kbp,ktb,kpb,kbr,kbb=',kbp,ktb,kpb,kbr,kbb
      if(dowr) write(outfile,131) 'kep,kte,kpe,ker,keb=',kep,kte,kpe,ker,keb

131   format(1x,a20,5(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,132) 'ibl,ibph,ibp3     =',ibl,ibph,ibp3
      if(dowr) write(outfile,132) 'iel,ieph,iep3     =',iel,ieph,iep3
      if(dowr) write(outfile,132) 'jbl,jbph,jbp3     =',jbl,jbph,jbp3
      if(dowr) write(outfile,132) 'jel,jeph,jep3     =',jel,jeph,jep3
      if(dowr) write(outfile,136) '    kbph,kbp3     ='    ,kbph,kbp3
      if(dowr) write(outfile,136) '    keph,kep3     ='    ,keph,kep3

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'np3a,np3o = ',np3a,np3o

132   format(1x,a19,   3(4x,i5))
136   format(1x,a19,9x,2(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,138) 'ibnba             =',ibnba
      if(dowr) write(outfile,138) 'ienba             =',ienba
      if(dowr) write(outfile,138) 'jbnba             =',jbnba
      if(dowr) write(outfile,138) 'jenba             =',jenba
      if(dowr) write(outfile,138) 'kbnba             =',kbnba
      if(dowr) write(outfile,138) 'kenba             =',kenba

138   format(1x,a19,   1(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,133) imp,jmp,kmp,kmt,rmp,cmp

133   format(' imp,jmp,kmp,kmt,rmp,cmp =',6(2x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ngxy,ngz = ',ngxy,ngz

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,135) 'ibdt,ibdq,ibdv    =',ibdt,ibdq,ibdv
      if(dowr) write(outfile,135) 'iedt,iedq,iedv    =',iedt,iedq,iedv
      if(dowr) write(outfile,135) 'jbdt,jbdq,jbdv    =',jbdt,jbdq,jbdv
      if(dowr) write(outfile,135) 'jedt,jedq,jedv    =',jedt,jedq,jedv
      if(dowr) write(outfile,135) 'kbdt,kbdq,kbdv    =',kbdt,kbdq,kbdv
      if(dowr) write(outfile,135) 'kedt,kedq,kedv    =',kedt,kedq,kedv

135   format(1x,a19,3(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,135) 'ibdk,ibdp,ibmynn  =',ibdk,ibdp,ibmynn
      if(dowr) write(outfile,135) 'iedk,iedp,iemynn  =',iedk,iedp,iemynn
      if(dowr) write(outfile,135) 'jbdk,jbdp,jbmynn  =',jbdk,jbdp,jbmynn
      if(dowr) write(outfile,135) 'jedk,jedp,jemynn  =',jedk,jedp,jemynn
      if(dowr) write(outfile,135) 'kbdk,kbdp,kbmynn  =',kbdk,kbdp,kbmynn
      if(dowr) write(outfile,135) 'kedk,kedp,kemynn  =',kedk,kedp,kemynn

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,137) 'ibmyj,ibtmp,ibcu  =',ibmyj,ibtmp,ibcu
      if(dowr) write(outfile,137) 'iemyj,ietmp,iecu  =',iemyj,ietmp,iecu
      if(dowr) write(outfile,137) 'jbmyj,jbtmp,jbcu  =',jbmyj,jbtmp,jbcu
      if(dowr) write(outfile,137) 'jemyj,jetmp,jecu  =',jemyj,jetmp,jecu
      if(dowr) write(outfile,137) 'kbmyj,kbtmp,kbcu  =',kbmyj,kbtmp,kbcu
      if(dowr) write(outfile,137) 'kemyj,ketmp,kecu  =',kemyj,ketmp,kecu
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,139) 'ntmp2d            =',ntmp2d
      if(dowr) write(outfile,139) 'ntmp1d            =',ntmp1d

137   format(1x,a19,3(4x,i5))
139   format(1x,a19,1(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  use_pbl = ',use_pbl

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) ' ib3d,ie3d,ib2d,ie2d = ',ib3d,ie3d,ib2d,ie2d
      if(dowr) write(outfile,*) ' jb3d,je3d,jb2d,je2d = ',jb3d,je3d,jb2d,je2d
      if(dowr) write(outfile,*) ' kb3d,ke3d           = ',kb3d,ke3d
      if(dowr) write(outfile,*) ' nout3d,nout2d       = ',nout3d,nout2d

    !----------

      if( dowr )then
        if( ntdiag.ge.1 )then
                                write(outfile,*)
                                write(outfile,*) '  ntdiag    = ',ntdiag
          if( td_hadv   .ge.1 ) write(outfile,*) '  td_hadv   = ',td_hadv
          if( td_vadv   .ge.1 ) write(outfile,*) '  td_vadv   = ',td_vadv
          if( td_hturb  .ge.1 ) write(outfile,*) '  td_hturb  = ',td_hturb
          if( td_vturb  .ge.1 ) write(outfile,*) '  td_vturb  = ',td_vturb
          if( td_hidiff .ge.1 ) write(outfile,*) '  td_hidiff = ',td_hidiff
          if( td_vidiff .ge.1 ) write(outfile,*) '  td_vidiff = ',td_vidiff
          if( td_hediff .ge.1 ) write(outfile,*) '  td_hediff = ',td_hediff
          if( td_vediff .ge.1 ) write(outfile,*) '  td_vediff = ',td_vediff
          if( td_mp     .ge.1 ) write(outfile,*) '  td_mp     = ',td_mp
          if( td_rdamp  .ge.1 ) write(outfile,*) '  td_rdamp  = ',td_rdamp
          if( td_nudge  .ge.1 ) write(outfile,*) '  td_nudge  = ',td_nudge
          if( td_rad    .ge.1 ) write(outfile,*) '  td_rad    = ',td_rad
          if( td_div    .ge.1 ) write(outfile,*) '  td_div    = ',td_div
          if( td_diss   .ge.1 ) write(outfile,*) '  td_diss   = ',td_diss
          if( td_pbl    .ge.1 ) write(outfile,*) '  td_pbl    = ',td_pbl
          if( td_lsw    .ge.1 ) write(outfile,*) '  td_lsw    = ',td_lsw
          if( td_efall  .ge.1 ) write(outfile,*) '  td_efall  = ',td_efall
          if( ptype.eq.5 )then
            if( ntdiag.ge.2 )     write(outfile,*)
            if( td_cond   .ge.1 ) write(outfile,*) '  td_cond   = ',td_cond
            if( td_evac   .ge.1 ) write(outfile,*) '  td_evac   = ',td_evac
            if( td_evar   .ge.1 ) write(outfile,*) '  td_evar   = ',td_evar
            if( td_dep    .ge.1 ) write(outfile,*) '  td_dep    = ',td_dep
            if( td_subl   .ge.1 ) write(outfile,*) '  td_subl   = ',td_subl
            if( td_melt   .ge.1 ) write(outfile,*) '  td_melt   = ',td_melt
            if( td_frz    .ge.1 ) write(outfile,*) '  td_frz    = ',td_frz
          endif
        endif

        if( nqdiag.ge.1 )then
                                write(outfile,*)
                                write(outfile,*) '  nqdiag    = ',nqdiag
          if( qd_dbz    .ge.1 ) write(outfile,*) '  qd_dbz    = ',qd_dbz
          if( qd_vtc    .ge.1 ) write(outfile,*) '  qd_vtc    = ',qd_vtc
          if( qd_vtr    .ge.1 ) write(outfile,*) '  qd_vtr    = ',qd_vtr
          if( qd_vts    .ge.1 ) write(outfile,*) '  qd_vts    = ',qd_vts
          if( qd_vtg    .ge.1 ) write(outfile,*) '  qd_vtg    = ',qd_vtg
          if( qd_vti    .ge.1 ) write(outfile,*) '  qd_vti    = ',qd_vti
          if( qd_hadv   .ge.1 ) write(outfile,*) '  qd_hadv   = ',qd_hadv
          if( qd_vadv   .ge.1 ) write(outfile,*) '  qd_vadv   = ',qd_vadv
          if( qd_hturb  .ge.1 ) write(outfile,*) '  qd_hturb  = ',qd_hturb
          if( qd_vturb  .ge.1 ) write(outfile,*) '  qd_vturb  = ',qd_vturb
          if( qd_hidiff .ge.1 ) write(outfile,*) '  qd_hidiff = ',qd_hidiff
          if( qd_vidiff .ge.1 ) write(outfile,*) '  qd_vidiff = ',qd_vidiff
          if( qd_hediff .ge.1 ) write(outfile,*) '  qd_hediff = ',qd_hediff
          if( qd_vediff .ge.1 ) write(outfile,*) '  qd_vediff = ',qd_vediff
          if( qd_mp     .ge.1 ) write(outfile,*) '  qd_mp     = ',qd_mp
          if( qd_nudge  .ge.1 ) write(outfile,*) '  qd_nudge  = ',qd_nudge
          if( qd_pbl    .ge.1 ) write(outfile,*) '  qd_pbl    = ',qd_pbl
          if( qd_lsw    .ge.1 ) write(outfile,*) '  qd_lsw    = ',qd_lsw
          if( ptype.eq.5 )then
            if( qd_cond   .ge.1 ) write(outfile,*) '  qd_cond   = ',qd_cond
            if( qd_evac   .ge.1 ) write(outfile,*) '  qd_evac   = ',qd_evac
            if( qd_evar   .ge.1 ) write(outfile,*) '  qd_evar   = ',qd_evar
            if( qd_dep    .ge.1 ) write(outfile,*) '  qd_dep    = ',qd_dep
            if( qd_subl   .ge.1 ) write(outfile,*) '  qd_subl   = ',qd_subl
          endif
        endif

        if( nudiag.ge.1 )then
                           write(outfile,*)
                           write(outfile,*) '  nudiag    = ',nudiag
          if( ud_hadv .ge.1 ) write(outfile,*) '  ud_hadv   = ',ud_hadv
          if( ud_vadv .ge.1 ) write(outfile,*) '  ud_vadv   = ',ud_vadv
          if( ud_hturb.ge.1 ) write(outfile,*) '  ud_hturb  = ',ud_hturb
          if( ud_vturb.ge.1 ) write(outfile,*) '  ud_vturb  = ',ud_vturb
          if( ud_hidiff.ge.1) write(outfile,*) '  ud_hidiff = ',ud_hidiff
          if( ud_vidiff.ge.1) write(outfile,*) '  ud_vidiff = ',ud_vidiff
          if( ud_hediff.ge.1) write(outfile,*) '  ud_hediff = ',ud_hediff
          if( ud_vediff.ge.1) write(outfile,*) '  ud_vediff = ',ud_vediff
          if( ud_pgrad.ge.1 ) write(outfile,*) '  ud_pgrad  = ',ud_pgrad
          if( ud_rdamp.ge.1 ) write(outfile,*) '  ud_rdamp  = ',ud_rdamp
          if( ud_nudge.ge.1 ) write(outfile,*) '  ud_nudge  = ',ud_nudge
          if( ud_cor  .ge.1 ) write(outfile,*) '  ud_cor    = ',ud_cor
          if( ud_cent .ge.1 ) write(outfile,*) '  ud_cent   = ',ud_cent
          if( ud_pbl  .ge.1 ) write(outfile,*) '  ud_pbl    = ',ud_pbl
          if( ud_lsw  .ge.1 ) write(outfile,*) '  ud_lsw    = ',ud_lsw
          if( nutk    .ge.1 ) write(outfile,*) '  nutk      = ',nutk
        endif

        if( nvdiag.ge.1 )then
                           write(outfile,*)
                           write(outfile,*) '  nvdiag    = ',nvdiag
          if( vd_hadv .ge.1 ) write(outfile,*) '  vd_hadv   = ',vd_hadv
          if( vd_vadv .ge.1 ) write(outfile,*) '  vd_vadv   = ',vd_vadv
          if( vd_hturb.ge.1 ) write(outfile,*) '  vd_hturb  = ',vd_hturb
          if( vd_vturb.ge.1 ) write(outfile,*) '  vd_vturb  = ',vd_vturb
          if( vd_hidiff.ge.1) write(outfile,*) '  vd_hidiff = ',vd_hidiff
          if( vd_vidiff.ge.1) write(outfile,*) '  vd_vidiff = ',vd_vidiff
          if( vd_hediff.ge.1) write(outfile,*) '  vd_hediff = ',vd_hediff
          if( vd_vediff.ge.1) write(outfile,*) '  vd_vediff = ',vd_vediff
          if( vd_pgrad.ge.1 ) write(outfile,*) '  vd_pgrad  = ',vd_pgrad
          if( vd_rdamp.ge.1 ) write(outfile,*) '  vd_rdamp  = ',vd_rdamp
          if( vd_nudge.ge.1 ) write(outfile,*) '  vd_nudge  = ',vd_nudge
          if( vd_cor  .ge.1 ) write(outfile,*) '  vd_cor    = ',vd_cor
          if( vd_cent .ge.1 ) write(outfile,*) '  vd_cent   = ',vd_cent
          if( vd_pbl  .ge.1 ) write(outfile,*) '  vd_pbl    = ',vd_pbl
          if( vd_lsw  .ge.1 ) write(outfile,*) '  vd_lsw    = ',vd_lsw
          if( nvtk    .ge.1 ) write(outfile,*) '  nvtk      = ',nvtk
        endif

        if( nwdiag.ge.1 )then
                           write(outfile,*)
                           write(outfile,*) '  nwdiag    = ',nwdiag
          if( wd_hadv .ge.1 ) write(outfile,*) '  wd_hadv   = ',wd_hadv
          if( wd_vadv .ge.1 ) write(outfile,*) '  wd_vadv   = ',wd_vadv
          if( wd_hturb.ge.1 ) write(outfile,*) '  wd_hturb  = ',wd_hturb
          if( wd_vturb.ge.1 ) write(outfile,*) '  wd_vturb  = ',wd_vturb
          if( wd_hidiff.ge.1) write(outfile,*) '  wd_hidiff = ',wd_hidiff
          if( wd_vidiff.ge.1) write(outfile,*) '  wd_vidiff = ',wd_vidiff
          if( wd_hediff.ge.1) write(outfile,*) '  wd_hediff = ',wd_hediff
          if( wd_vediff.ge.1) write(outfile,*) '  wd_vediff = ',wd_vediff
          if( wd_pgrad.ge.1 ) write(outfile,*) '  wd_pgrad  = ',wd_pgrad
          if( wd_rdamp.ge.1 ) write(outfile,*) '  wd_rdamp  = ',wd_rdamp
          if( wd_buoy .ge.1 ) write(outfile,*) '  wd_buoy   = ',wd_buoy
          if( nwtk    .ge.1 ) write(outfile,*) '  nwtk      = ',nwtk
        endif

        if( nkdiag.ge.1 )then
                           write(outfile,*)
                           write(outfile,*) '  nkdiag    = ',nkdiag
          if( kd_adv  .ge.1  ) write(outfile,*) '  kd_adv    = ',kd_adv
          if( kd_turb .ge.1  ) write(outfile,*) '  kd_turb   = ',kd_turb
        endif

        if( npdiag.ge.1 )then
                           write(outfile,*)
                           write(outfile,*) '  npdiag    = ',npdiag
        endif

      endif

!----------

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'rdx    =',rdx
      if(dowr) write(outfile,*) 'rdy    =',rdy
      if(dowr) write(outfile,*) 'rdz    =',rdz
      if(dowr) write(outfile,*) 'rdx2   =',rdx2
      if(dowr) write(outfile,*) 'rdy2   =',rdy2
      if(dowr) write(outfile,*) 'rdz2   =',rdz2
      if(dowr) write(outfile,*) 'rdx4   =',rdx4
      if(dowr) write(outfile,*) 'rdy4   =',rdy4
      if(dowr) write(outfile,*) 'rdz4   =',rdz4
      if(dowr) write(outfile,*) 'govtwo =',govtwo
      if(dowr) write(outfile,*) 'clwsat =',clwsat
      if(dowr) write(outfile,*) 'smeps  =',smeps
      if(dowr) write(outfile,*) 'tsmall =',tsmall
      if(dowr) write(outfile,*) 'qsmall =',qsmall
      if(dowr) write(outfile,*) 'csmax  =',csmax
      if(dowr) write(outfile,*) 'epsilon=',epsilon
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'nrkmax      =',nrkmax
      if(dowr) write(outfile,*) 'grads_undef =',grads_undef

      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      do i=ib,ie
        uh(i)=1.0
      enddo

      do i=ib,ie+1
        uf(i)=1.0
      enddo

      strx:  IF(stretch_x.ge.1)THEN

!!!        ibw=0
!!!        ibe=0

        ni1 = 0
        ni2 = 0
        ni3 = 0

!-----------------------------------------------------------------------
!  Begin specify xfdp

        nominal_dx = 0.5*( dx_inner + dx_outer )

      IF(stretch_x.eq.1)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 1 ... stretching on both west and east sides of domain:'
        if(dowr) write(outfile,*)
        ni1 = nint( (tot_x_len-nos_x_len)*0.5/nominal_dx )
        ni2 = nint( nos_x_len/dx_inner )
        ni3 = ni1
        if(dowr) write(outfile,*) '  ni1,ni2,ni3 = ',(tot_x_len-nos_x_len)*0.5/nominal_dx,   &
                         nos_x_len/dx_inner,(tot_x_len-nos_x_len)*0.5/nominal_dx
        if(dowr) write(outfile,*) '    (note:  ni1,ni2,ni3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_x.eq.2)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 2 ... stretching on east side of domain only:'
        if(dowr) write(outfile,*)
        ni1 = 0
        ni2 = nint( nos_x_len/dx_inner )
        ni3 = nint( (tot_x_len-nos_x_len)/nominal_dx )
        if(dowr) write(outfile,*) '  ni1,ni2,ni3 = ',0.0,nos_x_len/dx_inner,(tot_x_len-nos_x_len)/nominal_dx
        if(dowr) write(outfile,*) '    (note:  ni1,ni2,ni3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_x.eq.3)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 3 ... grid is specified in file input_grid_x '
        if(dowr) write(outfile,*)
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x must be 1,2,3'
        if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

      if( stretch_x.eq.1 .or. stretch_x.eq.2 )then
        c2=(nominal_dx-dx_inner)/(nominal_dx*nominal_dx*float(ni3-1))
        c1=(dx_inner/nominal_dx)-c2*nominal_dx

        if(dowr) write(outfile,*) '  nominal_dx  = ',nominal_dx
        if(dowr) write(outfile,*) '  c1,c2       = ',c1,c2
        if(dowr) write(outfile,*)
      endif

        ! Test to see if stretched-grid settings all make sense:
      IF(stretch_x.eq.1)THEN
        test_len = dx_inner*ni2 + (ni1+ni3)*0.5*(dx_inner+dx_outer)
        if( nx.ne.(ni1+ni2+ni3) .or. ni1.lt.0 .or. ni2.lt.0 .or. ni3.lt.0 .or.  &
              abs(tot_x_len-test_len).gt.1.0e-6*tot_x_len )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nx = ',nx
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       ni1,ni2,ni3 = ',ni1,ni2,ni3
          if(dowr) write(outfile,*) '       ni1+ni2+ni3 = ',ni1+ni2+ni3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings ...'
          if(dowr) write(outfile,*) '       dx_inner  = ',dx_inner
          if(dowr) write(outfile,*) '       dx_outer  = ',dx_outer
          if(dowr) write(outfile,*) '       nos_x_len = ',nos_x_len
          if(dowr) write(outfile,*) '       tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be nx = ',(nos_x_len/dx_inner)+(tot_x_len-nos_x_len)/(0.5*(dx_inner+dx_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if ni1,ni2,ni3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*) '  test_len  = ',test_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
      ELSEIF(stretch_x.eq.2)THEN
        test_len = dx_inner*ni2 + ni3*0.5*(dx_inner+dx_outer)
        if( nx.ne.(ni1+ni2+ni3) .or. ni1.lt.0 .or. ni2.lt.0 .or. ni3.lt.0 .or.  &
              abs(tot_x_len-test_len).gt.1.0e-6*tot_x_len )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nx = ',nx
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       ni1,ni2,ni3 = ',ni1,ni2,ni3
          if(dowr) write(outfile,*) '       ni1+ni2+ni3 = ',ni1+ni2+ni3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value for these settings ...'
          if(dowr) write(outfile,*) '       dx_inner  = ',dx_inner
          if(dowr) write(outfile,*) '       dx_outer  = ',dx_outer
          if(dowr) write(outfile,*) '       nos_x_len = ',nos_x_len
          if(dowr) write(outfile,*) '       tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be nx = ',(nos_x_len/dx_inner)+(tot_x_len-nos_x_len)/(0.5*(dx_inner+dx_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if ni1,ni2,ni3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*) '  test_len  = ',test_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
      ENDIF

        mult = 0.0
        if(iorigin.eq.2) mult = 0.5

      IF(stretch_x.eq.1)THEN

        do i=ni1+1,ni1+ni2+1
            xfdp(i)=ni1*nominal_dx+(i-ni1-1)*dx_inner - mult*tot_x_len
        enddo
        do i=ni1+ni2+2,ni1+ni2+ni3+(ngxy+1)
            xfdp(i)=ni1*nominal_dx+(ni1+ni2+1-ni1-1)*dble(dx_inner)   &
                 +(c1+c2*dble(i-1-ni1-ni2)*nominal_dx)   &
                 *dble(i-1-ni1-ni2)*nominal_dx - mult*tot_x_len
        enddo
        do i=(1-ngxy),ni1
            xfdp(i)=ni1*nominal_dx+(ni1+1-ni1-1)*dble(dx_inner)    &
                 -(c1+c2*dble(ni1+1-i)*nominal_dx)   &
                 *dble(ni1+1-i)*nominal_dx - mult*tot_x_len
        enddo

      ELSEIF(stretch_x.eq.2)THEN

        do i=ni1+1,ni1+ni2+1
            xfdp(i)=ni1*nominal_dx+(i-ni1-1)*dx_inner - mult*tot_x_len
        enddo
        do i=ni1+ni2+2,ni1+ni2+ni3+ngxy
            xfdp(i)=ni1*nominal_dx+(ni1+ni2+1-ni1-1)*dble(dx_inner)   &
                 +(c1+c2*dble(i-1-ni1-ni2)*nominal_dx)   &
                 *dble(i-1-ni1-ni2)*nominal_dx - mult*tot_x_len
        enddo
        do i=(1-ngxy),ni1
            xfdp(i)=ni1*nominal_dx+(ni1+1-ni1-1)*dble(dx_inner)    &
                 -(c1+c2*dble(ni1+1-i)*nominal_dx)   &
                 *dble(ni1+1-i)*nominal_dx - mult*tot_x_len
        enddo

      ELSEIF( stretch_x.eq.3 )THEN

        ! New for cm1r19:
        ! User specifies scalar points in a file called "input_grid_x"
        ! (note:  specify x location of points in meters)

        if( myid.eq.0 )then

          open(unit=40,file='input_grid_x',status='old',err=9015)

          do i=1,nx
            read(40,*,err=9013) xhref(i)
            print *,'    i,xhref = ',i,xhref(i)
          enddo

          close(unit=40)

        if( axisymm.ne.1 )then
          ! work outward from center of domain :
          xfdp(nx/2+1) = 0.5*(xhref(nx/2)+xhref(nx/2+1))

          do i=(nx/2+1),(nx)
            xfdp(i+1) = 2.0*xhref(i) - xfdp(i)
          enddo
          do i=(nx/2),(1),(-1)
            xfdp(i) = 2.0*xhref(i) - xfdp(i+1)
          enddo
        else
          ! axisymm = 1
          xfdp(1) = 0.0
          do i=1,nx
            xfdp(i+1) = 2.0*xhref(i) - xfdp(i)
          enddo
        endif

        endif

#ifdef MPI
        call MPI_BCAST(xfdp(1),(nx+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

      ENDIF

!!!        if( xf(ib).lt.0.0  .and. wbc.ne.1 ) ibw=1
!!!        if( xf(ie).gt.maxx .and. ebc.ne.1 ) ibe=1

        IF(stretch_x.eq.1)THEN
          do i=1,ngxy
            xfdp(1-i)=xfdp(1)-i*dx_outer
          enddo
        ELSEIF(stretch_x.eq.2)THEN
          do i=1,ngxy
            xfdp(1-i)=xfdp(1)-i*dx_inner
          enddo
        ENDIF

          do i=1,ngxy
            xfdp(nx+1+i)=xfdp(nx+1)+i*dx_outer
          enddo

        if( stretch_x.eq.3 )then
          do i=1,ngxy
            xfdp(1-i) = xfdp(1)-i*(xfdp(2)-xfdp(1))
            xfdp(nx+1+i) = xfdp(nx+1)+i*(xfdp(nx+1)-xfdp(nx))
          enddo
        endif

!  End specify xfdp
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for xfref
!  here.  Do not change anything below here!
!
!  Note:  xfref stores the location of the staggered u points for
!  the entire domain [from x=(1-ngxy) to x=(nx+ngxy+1)] (note: this includes
!  the boundary points that extend 3 gridpoints beyond the
!  computational domain.
!
!-----------------------------------------------------------------------

      ENDIF  strx


        do i=ib,ie+1
          xf(i)=xfdp(i+myi1-1)
        enddo

        do i=ib,ie
          xh(i)=0.5d0*(xfdp(i+myi1)+xfdp(i+myi1-1))
          uh(i)=dble(dx)/(xfdp(i+myi1)-xfdp(i+myi1-1))
        enddo

        arh1 = 1.0
        arh2 = 1.0
        arf1 = 1.0
        arf2 = 1.0

      IF(axisymm.eq.1)THEN

        print *
        do i=ib,ie
          arh1(i) = ( xfdp(i  )/( 0.5d0*(xfdp(i+1)+xfdp(i)) ) )
          arh2(i) = ( xfdp(i+1)/( 0.5d0*(xfdp(i+1)+xfdp(i)) ) )
          print *,'  arh1,arh2 = ',i,arh1(i),arh2(i),0.5*(arh1(i)+arh2(i))
        enddo
        print *
        print *
        do i=ib+1,ie
          if( abs(xfdp(i)).le.smeps )then
            arf1(i) = 1.0
            arf2(i) = 1.0
          else
            arf1(i) = ( 0.5d0*(xfdp(i-1)+xfdp(i)) / xfdp(i) )
            arf2(i) = ( 0.5d0*(xfdp(i+1)+xfdp(i)) / xfdp(i) )
          endif
          print *,'  arf1,arf2 = ',i,arf1(i),arf2(i),0.5*(arf1(i)+arf2(i))
        enddo
        print *

      ENDIF

        do i=ib+1,ie
          uf(i)=dble(dx)/( 0.5d0*(xfdp(i+myi1-1)+xfdp(i+myi1)) &
                          -0.5d0*(xfdp(i+myi1-1)+xfdp(i+myi1-2)) )
        enddo

        if(ibw.eq.1)then
          do i=1-ngxy,0
            uf(i)=uf(1)
          enddo
        endif

        if(ibe.eq.1)then
          do i=ni+2,ni+1+ngxy
            uf(i)=uf(ni+1)
          enddo
        endif

      do i=ib,ie
        rxh(i)=1.0/(smeps+xh(i))
        ruh(i)=1.0/uh(i)
      enddo

      do i=ib,ie+1
        rxf(i)=1.0/(smeps+xf(i))
        ruf(i)=1.0/uf(i)
      enddo

      do i=1-ngxy,nx+ngxy+1
        xfref(i) = xfdp(i)
      enddo

      ! 190310:  dummy check
      do i=1,ni
        if( (xf(i+1)-xf(i)).lt.smeps )then
          print *
          print *,'  Error:  dx <= 0 '
          print *
          call stopcm1
        endif
      enddo

      minx = xfref(1)
      maxx = xfref(nx+1)
      centerx  =  minx + 0.5*(maxx-minx)

      do i = 1-ngxy , nx+ngxy
        xhref(i) = 0.5*(xfdp(i)+xfdp(i+1))
      enddo

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'x:'
      if(dowr) write(outfile,124)
      if(dowr) write(outfile,125)
#ifdef MPI
    ! for MPI runs without procfiles, print the entire domain info:
    IF(.not.procfiles)THEN
      do i=1-ngxy,1-1
        if(dowr) write(outfile,122) i,xfref(i),xhref(i),xfref(i+1)-xfref(i),dx/(0.5*(xfref(1)+xfref(2))-0.5*(xfref(0)+xfref(1))),dx/(xfref(i+1)-xfref(i)),'   x'
      enddo
      do i=1,nx
        if(dowr) write(outfile,122) i,xfref(i),xhref(i),xfref(i+1)-xfref(i),dx/(0.5*(xfref(i)+xfref(i+1))-0.5*(xfref(i-1)+xfref(i))),dx/(xfref(i+1)-xfref(i)),'    '
      enddo
      do i=nx+1,nx+ngxy
        if(dowr) write(outfile,122) i,xfref(i),xhref(i),xfref(i+1)-xfref(i),dx/(0.5*(xfref(nx+1)+xfref(nx+2))-0.5*(xfref(nx)+xfref(nx+1))),dx/(xfref(i+1)-xfref(i)),'   x'
      enddo
      if(dowr) write(outfile,123) nx+1+ngxy,xfref(nx+1+ngxy),dx/(0.5*(xfref(nx+1)+xfref(nx+2))-0.5*(xfref(nx)+xfref(nx+1)))
    ELSE
#endif
      do i=1-ngxy,0
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'   x'
      enddo
      do i=1,ni
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'    '
      enddo
      do i=ni+1,ni+ngxy
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'   x'
      enddo
      if(dowr) write(outfile,123) ie+1,xf(ie+1),uf(ie+1)
#ifdef MPI
    ENDIF
#endif
      if(dowr) write(outfile,*)

122   format(3x,i5,3x,f11.2,3x,f11.2,3x,f9.2,3x,f8.4,3x,f8.4,a4)
123   format(3x,i5,3x,f11.2,29x,f8.4)
124   format('      i         xf (m)       xh (m)     dx (m)     uf         uh')
125   format(' ---------------------------------------------------------------')

!--------------------------------------------------------------


!--------------------------------------------------------------

      do j=jb,je
        vh(j)=1.0
      enddo

      do j=jb,je+1
        vf(j)=1.0
      enddo

      stry:  IF(stretch_y.ge.1)THEN

!!!        ibs=0
!!!        ibn=0

        nj1 = 0
        nj2 = 0
        nj3 = 0

!-----------------------------------------------------------------------
!  Begin specify yfdp

        nominal_dy = 0.5*( dy_inner + dy_outer )

      IF(stretch_y.eq.1)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 1 ... stretching on both south and north sides of domain:'
        if(dowr) write(outfile,*)
        nj1 = nint( (tot_y_len-nos_y_len)*0.5/nominal_dy )
        nj2 = nint( nos_y_len/dy_inner )
        nj3 = nj1
        if(dowr) write(outfile,*) '  nj1,nj2,nj3 = ',(tot_y_len-nos_y_len)*0.5/nominal_dy,   &
                         nos_y_len/dy_inner,(tot_y_len-nos_y_len)*0.5/nominal_dy
        if(dowr) write(outfile,*) '    (note:  nj1,nj2,nj3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_y.eq.2)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 2 ... stretching on north side of domain only:'
        if(dowr) write(outfile,*)
        nj1 = 0
        nj2 = nint( nos_y_len/dy_inner )
        nj3 = nint( (tot_y_len-nos_y_len)/nominal_dy )
        if(dowr) write(outfile,*) '  nj1,nj2,nj3 = ',0.0,nos_y_len/dy_inner,(tot_y_len-nos_y_len)/nominal_dy
        if(dowr) write(outfile,*) '    (note:  nj1,nj2,nj3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_y.eq.3)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 3 ... grid is specified in file input_grid_y '
        if(dowr) write(outfile,*)
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y must be 1,2,3'
        if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      ENDIF

      if( stretch_y.eq.1 .or. stretch_y.eq.2 )then
        c2=(nominal_dy-dy_inner)/(nominal_dy*nominal_dy*float(nj3-1))
        c1=(dy_inner/nominal_dy)-c2*nominal_dy

        if(dowr) write(outfile,*) '  nominal_dy  = ',nominal_dy
        if(dowr) write(outfile,*) '  c1,c2       = ',c1,c2
        if(dowr) write(outfile,*)
      endif

        ! Test to see if stretched-grid settings all make sense:
      IF(stretch_y.eq.1)THEN
        test_len = dy_inner*nj2 + (nj1+nj3)*0.5*(dy_inner+dy_outer)
        if( ny.ne.(nj1+nj2+nj3) .or. nj1.lt.0 .or. nj2.lt.0 .or. nj3.lt.0 .or.  &
              abs(tot_y_len-test_len).gt.1.0e-6*tot_y_len )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of ny = ',ny
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       nj1,nj2,nj3 = ',nj1,nj2,nj3
          if(dowr) write(outfile,*) '       nj1+nj2+nj3 = ',nj1+nj2+nj3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings ...'
          if(dowr) write(outfile,*) '       dy_inner  = ',dy_inner
          if(dowr) write(outfile,*) '       dy_outer  = ',dy_outer
          if(dowr) write(outfile,*) '       nos_y_len = ',nos_y_len
          if(dowr) write(outfile,*) '       tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be ny = ',(nos_y_len/dy_inner)+(tot_y_len-nos_y_len)/(0.5*(dy_inner+dy_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if nj1,nj2,nj3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*) '  test_len  = ',test_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
      ELSEIF(stretch_y.eq.2)THEN
        if( ny.ne.(nj1+nj2+nj3) .or. nj1.lt.0 .or. nj2.lt.0 .or. nj3.lt.0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of ny = ',ny
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       nj1,nj2,nj3 = ',nj1,nj2,nj3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value for these settings ...'
          if(dowr) write(outfile,*) '       dy_inner  = ',dy_inner
          if(dowr) write(outfile,*) '       dy_outer  = ',dy_outer
          if(dowr) write(outfile,*) '       nos_y_len = ',nos_y_len
          if(dowr) write(outfile,*) '       tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be ny = ',(nos_y_len/dy_inner)+(tot_y_len-nos_y_len)/(0.5*(dy_inner+dy_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if nj1,nj2,nj3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif

      ELSEIF( stretch_y.eq.3 )THEN

        ! New for cm1r19:
        ! User specifies scalar points in a file called "input_grid_y"
        ! (note:  specify y location of points in meters)

        if( myid.eq.0 )then

          open(unit=40,file='input_grid_y',status='old',err=9016)

          do j=1,ny
            read(40,*,err=9014) yhref(j)
            print *,'    j,yhref = ',j,yhref(j)
          enddo

          close(unit=40)

          ! work outward from center of domain :
          yfdp(ny/2+1) = 0.5*(yhref(ny/2)+yhref(ny/2+1))

          do j=(ny/2+1),(ny)
            yfdp(j+1) = 2.0*yhref(j) - yfdp(j)
          enddo
          do j=(ny/2),(1),(-1)
            yfdp(j) = 2.0*yhref(j) - yfdp(j+1)
          enddo

        endif

#ifdef MPI
        call MPI_BCAST(yfdp(1),(ny+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

      ENDIF

        mult = 0.0
        if(iorigin.eq.2) mult = 0.5

      IF(stretch_y.eq.1)THEN

        do j=nj1+1,nj1+nj2+1
            yfdp(j)=nj1*nominal_dy+(j-nj1-1)*dy_inner - mult*tot_y_len
        enddo
        do j=nj1+nj2+2,nj1+nj2+nj3+(ngxy+1)
            yfdp(j)=nj1*nominal_dy+(nj1+nj2+1-nj1-1)*dble(dy_inner)   &
                 +(c1+c2*dble(j-1-nj1-nj2)*nominal_dy)   &
                 *dble(j-1-nj1-nj2)*nominal_dy - mult*tot_y_len
        enddo
        do j=(1-ngxy),nj1
            yfdp(j)=nj1*nominal_dy+(nj1+1-nj1-1)*dble(dy_inner)    &
                 -(c1+c2*dble(nj1+1-j)*nominal_dy)   &
                 *dble(nj1+1-j)*nominal_dy - mult*tot_y_len
        enddo

      ELSEIF(stretch_y.eq.2)THEN

        do j=nj1+1,nj1+nj2+1
            yfdp(j)=nj1*nominal_dy+(j-nj1-1)*dy_inner - mult*tot_y_len
        enddo
        do j=nj1+nj2+2,nj1+nj2+nj3+ngxy
            yfdp(j)=nj1*nominal_dy+(nj1+nj2+1-nj1-1)*dble(dy_inner)   &
                 +(c1+c2*dble(j-1-nj1-nj2)*nominal_dy)   &
                 *dble(j-1-nj1-nj2)*nominal_dy - mult*tot_y_len
        enddo
        do j=(1-ngxy),nj1
            yfdp(j)=nj1*nominal_dy+(nj1+1-nj1-1)*dble(dy_inner)    &
                 -(c1+c2*dble(nj1+1-j)*nominal_dy)   &
                 *dble(nj1+1-j)*nominal_dy - mult*tot_y_len
        enddo

      ENDIF

!!!        if( yf(jb).lt.0.0  .and. sbc.ne.1 ) ibs=1
!!!        if( yf(je).gt.maxy .and. nbc.ne.1 ) ibn=1

        IF(stretch_y.eq.1)THEN
          do j=1,ngxy
            yfdp(1-j)=yfdp(1)-j*dy_outer
          enddo
        ELSEIF(stretch_y.eq.2)THEN
          do j=1,ngxy
            yfdp(1-j)=yfdp(1)-j*dy_inner
          enddo
        ENDIF

          do j=1,ngxy
            yfdp(ny+1+j)=yfdp(ny+1)+j*dy_outer
          enddo

        if( stretch_y.eq.3 )then
          do j=1,ngxy
            yfdp(1-j) = yfdp(1)-j*(yfdp(2)-yfdp(1))
            yfdp(ny+1+j) = yfdp(ny+1)+j*(yfdp(ny+1)-yfdp(ny))
          enddo
        endif

!  End specify yfdp
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for yfref
!  here.  Do not change anything below here!
!
!  Note:  yfref stores the location of the staggered v points for
!  the entire domain [from y=(1-ngxy) to y=(ny+ngxy+1)] (note: this includes
!  the boundary points that extend 3 gridpoints beyond the
!  computational domain.
!
!-----------------------------------------------------------------------

      ENDIF  stry

        do j=jb,je+1
          yf(j)=yfdp(j+myj1-1)
        enddo

        do j=jb,je
          yh(j)=0.5d0*(yfdp(j+myj1)+yfdp(j+myj1-1))
          vh(j)=dble(dy)/(yfdp(j+myj1)-yfdp(j+myj1-1))
        enddo

        do j=jb+1,je
          vf(j)=dble(dy)/( 0.5d0*(yfdp(j+myj1-1)+yfdp(j+myj1)) &
                          -0.5d0*(yfdp(j+myj1-1)+yfdp(j+myj1-2)) )
        enddo

        if(ibs.eq.1)then
          do j=1-ngxy,0
            vf(j)=vf(1)
          enddo
        endif

        if(ibn.eq.1)then
          do j=nj+2,nj+1+ngxy
            vf(j)=vf(nj+1)
          enddo
        endif

      do j=jb,je
        rvh(j)=1.0/vh(j)
      enddo

      do j=jb,je+1
        rvf(j)=1.0/vf(j)
      enddo

      do j=1-ngxy,ny+ngxy+1
        yfref(j) = yfdp(j)
      enddo

      ! 190310:  dummy check
      do j=1,nj
        if( (yf(j+1)-yf(j)).lt.smeps )then
          print *
          print *,'  Error:  dy <= 0 '
          print *
          call stopcm1
        endif
      enddo

      miny = yfref(1)
      maxy = yfref(ny+1)
      centery  =  miny + 0.5*(maxy-miny)

      do j = 1-ngxy , ny+ngxy
        yhref(j) = 0.5*(yfdp(j)+yfdp(j+1))
      enddo

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'y:'
      if(dowr) write(outfile,134)
      if(dowr) write(outfile,125)
#ifdef MPI
    ! for MPI runs without procfiles, print the entire domain info:
    IF(.not.procfiles)THEN
      do j=1-ngxy,1-1
        if(dowr) write(outfile,122) j,yfref(j),yhref(j),yfref(j+1)-yfref(j),dy/(0.5*(yfref(1)+yfref(2))-0.5*(yfref(0)+yfref(1))),dy/(yfref(j+1)-yfref(j)),'   x'
      enddo
      do j=1,ny
        if(dowr) write(outfile,122) j,yfref(j),yhref(j),yfref(j+1)-yfref(j),dy/(0.5*(yfref(j)+yfref(j+1))-0.5*(yfref(j-1)+yfref(j))),dy/(yfref(j+1)-yfref(j)),'    '
      enddo
      do j=ny+1,ny+ngxy
        if(dowr) write(outfile,122) j,yfref(j),yhref(j),yfref(j+1)-yfref(j),dy/(0.5*(yfref(ny+1)+yfref(ny+2))-0.5*(yfref(ny)+yfref(ny+1))),dy/(yfref(j+1)-yfref(j)),'   x'
      enddo
      if(dowr) write(outfile,123) ny+1+ngxy,yfref(ny+1+ngxy),dy/(0.5*(yfref(ny+1)+yfref(ny+2))-0.5*(yfref(ny)+yfref(ny+1)))
    ELSE
#endif
      do j=1-ngxy,0
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'   x'
      enddo
      do j=1,nj
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'    '
      enddo
      do j=nj+1,nj+ngxy
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'   x'
      enddo
      if(dowr) write(outfile,123) je+1,yf(je+1),vf(je+1)
#ifdef MPI
    ENDIF
#endif
      if(dowr) write(outfile,*)

134   format('      j         yf (m)       yh (m)     dy (m)     vf         vh')

!--------------------------------------------------------------

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        zf(i,j,k)=dz*(k-1)
        mf(i,j,k)=1.0
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        zh(i,j,k)=0.5*(zf(i,j,k)+zf(i,j,k+1))
        mh(i,j,k)=1.0
      enddo
      enddo
      enddo

      do k=kb,ke+1
        sigmaf(k)=zf(1,1,k)
      enddo
      do k=kb,ke
        sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
      enddo

    IF(stretch_z.ge.1)THEN

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

      strz:  IF ( stretch_z == 1 ) THEN

        nominal_dz = 0.5*(dz_bot+dz_top)

        nk1 = nint( str_bot/dz_bot )
        nk3 = nint( (ztop-str_top)/dz_top )
        nk2 = nk-(nk1+nk3)

        ! dummy checks:
        if(dowr) write(outfile,*) '  bot: ',nk1*dz_bot,str_bot,nk1*dz_bot-str_bot
        if( abs(nk1*dz_bot-str_bot).gt.0.01 )then
          if(dowr) write(outfile,*) '  depth of bottom layer does not exactly divide by dz_bot! '
          if(dowr) write(outfile,*) '  nk1*dz_bot = ',nk1*dz_bot
          if(dowr) write(outfile,*) '  str_bot    = ',str_bot
          if(dowr) write(outfile,*) '  diff       = ',nk1*dz_bot-str_bot
          if(dowr) write(outfile,*) '  stopping cm1 ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
        if(dowr) write(outfile,*) '  mid: ',nk2*nominal_dz,(str_top-str_bot),nk2*nominal_dz-(str_top-str_bot)
        if( abs(nk2*nominal_dz-(str_top-str_bot)).ge.0.01 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) amod(str_top-str_bot,nominal_dz),1.0e-6*(str_top-str_bot)
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  depth of middle layer does not exactly divide by nominal_dz! '
          if(dowr) write(outfile,*) '  nk2*nominal_dz  = ',nk2*nominal_dz
          if(dowr) write(outfile,*) '  str_top-str_bot = ',str_top-str_bot
          if(dowr) write(outfile,*) '  diff            = ',nk2*nominal_dz-(str_top-str_bot)
          if(dowr) write(outfile,*) '  stopping cm1 ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
        if(dowr) write(outfile,*) '  top: ',nk3*dz_top,(ztop-str_top),nk3*dz_top-(ztop-str_top)
        if( abs(nk3*dz_top-(ztop-str_top)).gt.0.01 )then
          if(dowr) write(outfile,*) '  depth of top layer does not exactly divide by dz_top! '
          if(dowr) write(outfile,*) '  nk3*dz_top   = ',nk3*dz_top
          if(dowr) write(outfile,*) '  ztop-str_top = ',ztop-str_top
          if(dowr) write(outfile,*) '  diff         = ',nk3*dz_top-(ztop-str_top)
          if(dowr) write(outfile,*) '  stopping cm1 ... '
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif
        if( (nk1+nk2+nk3)-nk .ne. 0 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  (nk1+nk2+nk3) does not equal nk '
          if(dowr) write(outfile,*) '  (nk1+nk2+nk3) = ',(nk1+nk2+nk3)
          if(dowr) write(outfile,*) '   nk           = ',nk
          if(dowr) write(outfile,*)
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif

        nominal_dz=(str_top-str_bot)/nk2

        c2=(nominal_dz-dz_bot)/(nominal_dz*nominal_dz*float(nk2-1))
        c1=(dz_bot/nominal_dz)-c2*nominal_dz

        ! Test to see if stretched-grid settings all make sense:
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  actual-nk,test-nk:',float(nk),float(nk1+nk3)+(str_top-str_bot)/(0.5*(dz_bot+dz_top))
        test_len = dz_bot*nk1 + 0.5*(dz_bot+dz_top)*nk2 + dz_top*nk3
        if( nz.ne.(nk1+nk2+nk3) .or. nk1.lt.0 .or. nk2.lt.0 .or. nk3.lt.0 .or.  &
              abs(ztop-test_len).gt.1.0e-6*ztop )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nz = ',nz
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '       nk1,nk2,nk3 = ',nk1,nk2,nk3
          if(dowr) write(outfile,*) '       nk1+nk2+nk3 = ',nk1+nk2+nk3
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings:'
          if(dowr) write(outfile,*) '       ztop      = ',ztop
          if(dowr) write(outfile,*) '       str_bot   = ',str_bot
          if(dowr) write(outfile,*) '       str_top   = ',str_top
          if(dowr) write(outfile,*) '       dz_bot    = ',dz_bot
          if(dowr) write(outfile,*) '       dz_top    = ',dz_top
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  would be nz = ',nk1+nk3+(str_top-str_bot)/(0.5*(dz_bot+dz_top))
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ztop     = ',ztop
          if(dowr) write(outfile,*) '  test_len = ',test_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  nk1,nk2,nk3,ntot=',nk1,nk2,nk3,(nk1+nk2+nk3)
        if(dowr) write(outfile,*) '  nominal_dz =',nominal_dz
        if(dowr) write(outfile,*) '  c1,c2 = ',c1,c2
        if(dowr) write(outfile,*)

      do j=jb,je
      do i=ib,ie

        do k=1,nk1+1
          zf(i,j,k)=(k-1)*dz_bot
        enddo
        do k=(nk1+1),(nk1+nk2+1)
          zf(i,j,k)=zf(i,j,nk1+1)+(c1+c2*float(k-1-nk1)*nominal_dz)   &
                         *float(k-1-nk1)*nominal_dz
        enddo
        do k=(nk1+nk2+2),(nk1+nk2+nk3+1)
          zf(i,j,k)=zf(i,j,k-1)+dz_top
        enddo

      enddo
      enddo

!!!      if(terrain_flag)then

        do k=1,nk1+1
          sigmaf(k)=(k-1)*dz_bot
        enddo
        do k=(nk1+1),(nk1+nk2+1)
          sigmaf(k)=sigmaf(nk1+1)+(c1+c2*float(k-1-nk1)*nominal_dz)   &
                         *float(k-1-nk1)*nominal_dz
        enddo
        do k=(nk1+nk2+2),(nk1+nk2+nk3+1)
          sigmaf(k)=sigmaf(k-1)+dz_top
        enddo

        sigmaf(0)=-sigmaf(2)
        sigmaf(nk+2)=sigmaf(nk+1)+(sigmaf(nk+1)-sigmaf(nk))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ELSEIF ( stretch_z .eq. 2 ) THEN ! geometric stretching from COMMAS

! nz is ke (number of w points)
! kb = 0
! so ng = 1 and zfw1d(-ng+1:nz+ng)
! SUBROUTINE ZGRID(dz, dz_stretch, nz, nbndlyr, gzc, gze, ng, dzmax,ztopstr,rtop,dzmaxtop)
       
       nbndlyr = Max(0,Int( str_bot/dz_bot + 0.01) - 1)
       CALL ZGRID(dz, dz_bot, ke, nbndlyr, gzc, gze, 1, dz_top,ztopstr,rtop,dzmaxtop,myid)
       
         IF ( myid == 0 ) THEN
           DO k = 1,ke
             write(6,*) 'k,gzc,gze = ',k,gzc(k),gze(k)
           ENDDO
         ENDIF
         gze(0) = -gze(2)
         gze(ke+1) = gze(ke) + (gze(ke) - gze(ke-1))

        DO k = kb,ke+1
         DO j = jb,je
          DO i = ib,ie
           zf(i,j,k) = gze(k)
          ENDDO
         ENDDO
        ENDDO

        if(terrain_flag)then
          write(0,*) 'terrain not yet compatible with stretch_z == 2'
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        ENDIF

        do k=kb,ke+1
          sigmaf(k)=zf(1,1,k)
        enddo
        do k=kb,ke
          sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
        enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ELSEIF ( stretch_z .eq. 3 ) THEN

        ! New for cm1r19:
        ! User specifies scalar levels in a file called "input_grid_z"
        ! (note:  specify levels in meters ASL)
        ! See run/config_files/les_StratoCuDrizzle/input_grid_z for an example.

        ! NOTE:  this section of code expects the model heights on HALF LEVELS 
        !        (aka scalar levels), eg, 0.5*dz, 1.5*dz, 2.5*dz, etc 

        if( myid.eq.0 )then

          open(unit=40,file='input_grid_z',status='old',err=9011)

          sigmaf(1) = 0.0

          do k=1,nk
            print *,'  k = ',k
            read(40,*,err=9012) sigma(k)
            print *,'    sigma = ',sigma(k)
            sigmaf(k+1) = 2.0*sigma(k) - sigmaf(k)
          enddo

          close(unit=40)

        endif

#ifdef MPI
        call MPI_BCAST(sigma(1) ,nk  ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(sigmaf(1),nk+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

        ztop = sigmaf(nk+1)

        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  ztop = ',ztop
        if(myid.eq.0) print *

        do k=1,nk+1
        do j=jb,je
        do i=ib,ie
          zf(i,j,k) = sigmaf(k)
        enddo
        enddo
        enddo

        do j=jb,je
        do i=ib,ie

          zf(i,j,0)=-zf(i,j,2)
          zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))

          do k=0,nk+1
            zh(i,j,k)=0.5*(zf(i,j,k+1)+zf(i,j,k))
          enddo
          zh(i,j,0)=-zh(i,j,1)
          zh(i,j,nk+1)=zh(i,j,nk)+2.0*(zf(i,j,nk+1)-zh(i,j,nk))

        enddo
        enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      ELSEIF ( stretch_z .eq. 4 ) THEN

        ! New for cm1r19.8
        ! User specifies w levels in a file called "input_grid_z"
        ! (note:  specify levels in meters ASL)

        ! NOTE:  this section of code expects the model heights on FULL LEVELS 
        !        (aka w levels), eg, 0.0*dz, 1.0*dz, 2.0*dz, 3.0*dz, etc 

        if( myid.eq.0 )then

          open(unit=40,file='input_grid_z',status='old',err=9011)

          do k=1,nk+1
            print *,'  k = ',k
            read(40,*,err=9012) sigmaf(k)
            print *,'    sigmaf = ',sigmaf(k)
          enddo

          do k=1,nk
            sigma(k) = 0.5*(sigmaf(k)+sigmaf(k+1))
          enddo

          close(unit=40)


        endif

#ifdef MPI
        call MPI_BCAST(sigma(1) ,nk  ,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(sigmaf(1),nk+1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

        ztop = sigmaf(nk+1)

        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  ztop = ',ztop
        if(myid.eq.0) print *

        do k=1,nk+1
        do j=jb,je
        do i=ib,ie
          zf(i,j,k) = sigmaf(k)
        enddo
        enddo
        enddo

        do j=jb,je
        do i=ib,ie

          zf(i,j,0)=-zf(i,j,2)
          zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))

          do k=0,nk+1
            zh(i,j,k)=0.5*(zf(i,j,k+1)+zf(i,j,k))
          enddo
          zh(i,j,0)=-zh(i,j,1)
          zh(i,j,nk+1)=zh(i,j,nk)+2.0*(zf(i,j,nk+1)-zh(i,j,nk))

        enddo
        enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      ENDIF  strz


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for zf
!  here.  Do not change anything below here!
!
!  Note:  zf stores the location of the staggered w points. 
!
!  Note:  if you are using terrain, you need to also specify the nominal 
!  locations of the zf points in the sigmaf array.
!
!-----------------------------------------------------------------------

      do j=jb,je
      do i=ib,ie

        zf(i,j,0)=-zf(i,j,2)
        zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))

        do k=0,nk+1
          zh(i,j,k)=0.5*(zf(i,j,k+1)+zf(i,j,k))
          mh(i,j,k)=dz/(zf(i,j,k+1)-zf(i,j,k))
        enddo
        zh(i,j,0)=-zh(i,j,1)
        zh(i,j,nk+1)=zh(i,j,nk)+2.0*(zf(i,j,nk+1)-zh(i,j,nk))

        do k=1,nk+1
          mf(i,j,k)=dz/(zh(i,j,k)-zh(i,j,k-1))
        enddo
        mf(i,j,0)=mf(i,j,1)
        mf(i,j,nk+2)=mf(i,j,nk+1)

      enddo
      enddo

    ENDIF

! end vertical stretching section
!-----------------------------------------------------------------------

      do k=kb,ke
        sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
      enddo

      maxz = sigmaf(nk+1)

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rmh(i,j,k)=1.0/mh(i,j,k)
      enddo
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rmf(i,j,k)=1.0/mf(i,j,k)
      enddo
      enddo
      enddo

      ! 190310:  dummy check
      do k=1,nk
      do j=1,nj
      do i=1,ni
        if( (zf(i,j,k+1)-zf(i,j,k)).lt.smeps )then
          print *
          print *,'  k,zf,zfp1,dz = ',k,zf(i,j,k),zf(i,j,k+1),zf(i,j,k+1)-zf(i,j,k)
          print *
          print *,'  Error:  dz <= 0 '
          print *
          call stopcm1
        endif
      enddo
      enddo
      enddo

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'model heights:'
      if(dowr) write(outfile,104)
104   format('     k       zf (m)     zh (m)     dz (m)     mf         mh')
      if(dowr) write(outfile,105)
105   format(' ---------------------------------------------------------------')
      do k=1,nk
        if(dowr) write(outfile,102) k,zf(1,1,k),zh(1,1,k),zf(1,1,k+1)-zf(1,1,k),mf(1,1,k),mh(1,1,k)
102     format(3x,i4,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.4,3x,f8.4)
      enddo
      if(dowr) write(outfile,103) nk+1,zf(1,1,nk+1),mf(1,1,nk+1)
103   format(3x,i4,3x,f8.2,25x,f8.4)
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------

#ifdef MPI
      if( ibw.eq.1 .and. ibs.eq.0 ) patchsww = .true.
      if( ibw.eq.1 .and. ibn.eq.0 ) patchnww = .true.
      if( ibe.eq.1 .and. ibs.eq.0 ) patchsee = .true.
      if( ibe.eq.1 .and. ibn.eq.0 ) patchnee = .true.
      if( ibs.eq.1 .and. ibw.eq.0 ) patchsws = .true.
      if( ibs.eq.1 .and. ibe.eq.0 ) patchses = .true.
      if( ibn.eq.1 .and. ibw.eq.0 ) patchnwn = .true.
      if( ibn.eq.1 .and. ibe.eq.0 ) patchnen = .true.

      if(dowr) write(outfile,*) '  patchsww =',patchsww
      if(dowr) write(outfile,*) '  patchnww =',patchnww
      if(dowr) write(outfile,*) '  patchsee =',patchsee
      if(dowr) write(outfile,*) '  patchnee =',patchnee
      if(dowr) write(outfile,*) '  patchsws =',patchsws
      if(dowr) write(outfile,*) '  patchses =',patchses
      if(dowr) write(outfile,*) '  patchnwn =',patchnwn
      if(dowr) write(outfile,*) '  patchnen =',patchnen
      if(dowr) write(outfile,*)

      if( ibw.eq.1 .and. ibs.eq.1 ) p2tchsww = .true.
      if( ibw.eq.1 .and. ibn.eq.1 ) p2tchnww = .true.
      if( ibe.eq.1 .and. ibs.eq.1 ) p2tchsee = .true.
      if( ibe.eq.1 .and. ibn.eq.1 ) p2tchnee = .true.
      if( ibs.eq.1 .and. ibw.eq.1 ) p2tchsws = .true.
      if( ibs.eq.1 .and. ibe.eq.1 ) p2tchses = .true.
      if( ibn.eq.1 .and. ibw.eq.1 ) p2tchnwn = .true.
      if( ibn.eq.1 .and. ibe.eq.1 ) p2tchnen = .true.

      if(dowr) write(outfile,*) '  p2tchsww =',p2tchsww
      if(dowr) write(outfile,*) '  p2tchnww =',p2tchnww
      if(dowr) write(outfile,*) '  p2tchsee =',p2tchsee
      if(dowr) write(outfile,*) '  p2tchnee =',p2tchnee
      if(dowr) write(outfile,*) '  p2tchsws =',p2tchsws
      if(dowr) write(outfile,*) '  p2tchses =',p2tchses
      if(dowr) write(outfile,*) '  p2tchnwn =',p2tchnwn
      if(dowr) write(outfile,*) '  p2tchnen =',p2tchnen
      if(dowr) write(outfile,*)
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                  BEGIN TERRAIN !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      rds  = 1.0
      rdsf = 1.0

      zs=0.0

      gz=1.0
      rgz=1.0
      gzu=1.0
      rgzu=1.0
      gzv=1.0
      rgzv=1.0
      dzdx=0.0
      dzdy=0.0

      gx=0.0
      gxu=0.0
      gy=0.0
      gyv=0.0

      IF(terrain_flag)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Terrain included!'
        if(dowr) write(outfile,*)

        ! moved this section of code to init_terrain in cm1r15:
        call init_terrain(xh,uh,xf,uf,yh,vh,yf,vf,rds,sigma,rdsf,sigmaf,  &
                          zh,zf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,         &
                          dzdx,dzdy,gx,gxu,gy,gyv,                   &
                          reqs_u,reqs_v,reqs_s,reqs_p,               &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,           &
                          sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,   &
                          uw31,uw32,ue31,ue32,us31,us32,un31,un32,   &
                          vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,   &
                          ww31(1,1,1),ww32(1,1,1),we31(1,1,1),we32(1,1,1), &
                          ws31(1,1,1),ws32(1,1,1),wn31(1,1,1),wn32(1,1,1))

      ENDIF

!--------------------------------------------------------------
!  180212:  immersed boundary stuff


      call ib_setup


      if( myid.eq.0 ) print *
      if( myid.eq.0 ) print *,'  do_ib   = ',do_ib
      if( myid.eq.0 ) print *,'  ib_init = ',ib_init
      if( myid.eq.0 ) print *,'  top_cd  = ',top_cd
      if( myid.eq.0 ) print *,'  side_cd = ',side_cd


      !-----
      if( terrain_flag .and. do_ib )then
        if(myid.eq.0)then
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  do_ib        = ',do_ib
        print *
        print *,'  cannot use terrain_flag and do_ib at the same time (for now ... will be fixed later)'
        print *
        print *,'  stopping cm1 .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      !-----
      IF( do_ib .and. doimpl.eq.1 )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  do_ib     = ',do_ib
        print *,'  doimpl    = ',doimpl
        print *
        print *,'  cannot use doimpl=1 with do_ib (for now) '
        print *
        print *,'  ... setting doimpl = 0 ... '
        print *
        print *,'  -------------------------------- '
        endif
        doimpl = 0
      ENDIF
      !-----
      IF( do_ib .and. ( hadvordrs.gt.5 .or. vadvordrs.gt.5 .or.  &
                        hadvordrv.gt.5 .or. vadvordrv.gt.5 .or.  &
                        weno_order.gt.5 ) )THEN
        if(myid.eq.0)then
        print *,'  -------------------------------- '
        print *
        print *,'  cannot use immersed boundaries with advorder > 5 '
        print *
        print *,'  ... setting hadvordrs  = 5 ... '
        print *,'  ... setting vadvordrs  = 5 ... '
        print *,'  ... setting hadvordrv  = 5 ... '
        print *,'  ... setting vadvordrv  = 5 ... '
        print *,'  ... setting weno_order = 5 ... '
        print *
        print *,'  -------------------------------- '
        endif
        hadvordrs = 5
        vadvordrs = 5
        hadvordrv = 5
        vadvordrv = 5
        weno_order = 5
      ENDIF
      !-----
      if( do_ib .and. ( sgsmodel.eq.3 .or. sgsmodel.eq.4 ) )then
        if(myid.eq.0)then
        print *
        print *,'  sgsmodel     = ',sgsmodel
        print *,'  do_ib        = ',do_ib
        print *
        print *,'  cannot use sgsmodel = 3,4 with do_ib '
        print *
        print *,'  stopping cm1 .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      !-----
      if( do_ib .and. ( psolver.eq.3 .or. psolver.eq.4 .or. psolver.eq.5 ) )then
        if(myid.eq.0)then
        print *
        print *,'  psolver      = ',psolver
        print *,'  do_ib        = ',do_ib
        print *
        print *,'  cannot use immersed boundaries with psolver = 3,4,5 (for now) '
        print *
        print *,'  stopping cm1 .... '
        print *
        endif
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        call stopcm1
      endif
      !-----

      if( do_ib )then

        ibib=ib
        ieib=ie
        jbib=jb
        jeib=je
        kbib=kb
        keib=ke

      else

        ibib=1
        ieib=1
        jbib=1
        jeib=1
        kbib=1
        keib=1

      endif

      if( myid.eq.0 )then
        print *
        print *,'    ibib,ieib = ',ibib,ieib
        print *,'    jbib,jeib = ',jbib,jeib
        print *,'    kbib,keib = ',kbib,keib
        print *
      endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                  END   TERRAIN !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        mh(i,j,k)=dz/(zf(i,j,k+1)-zf(i,j,k))
        rmh(i,j,k)=1.0/mh(i,j,k)
      enddo
      enddo
      enddo

      do k=kb+1,ke
      do j=jb,je
      do i=ib,ie
        mf(i,j,k)=dz/(zh(i,j,k)-zh(i,j,k-1))
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        mf(i,j,0)=mf(i,j,1)
        mf(i,j,nk+2)=mf(i,j,nk+1)
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rmf(i,j,k)=1.0/mf(i,j,k)
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  minx    = ',minx
      if(dowr) write(outfile,*) '  centerx = ',centerx
      if(dowr) write(outfile,*) '  maxx    = ',maxx
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  miny    = ',miny
      if(dowr) write(outfile,*) '  centery = ',centery
      if(dowr) write(outfile,*) '  maxy    = ',maxy
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  maxz = ',maxz
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) '  ibw =',ibw
      if(dowr) write(outfile,*) '  ibe =',ibe
      if(dowr) write(outfile,*) '  ibs =',ibs
      if(dowr) write(outfile,*) '  ibn =',ibn
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!  Get min/max dx,dy,dz on grid
!  (needed for adapt_dt ... but interesting to report, nontheless)

      min_dx = 1.0e20
      min_dy = 1.0e20
      min_dz = 1.0e20

      max_dx = 0.0
      max_dy = 0.0
      max_dz = 0.0

      do i=1,ni
        min_dx = min( min_dx , xf(i+1)-xf(i) )
        max_dx = max( max_dx , xf(i+1)-xf(i) )
      enddo

      do j=1,nj
        min_dy = min( min_dy , yf(j+1)-yf(j) )
        max_dy = max( max_dy , yf(j+1)-yf(j) )
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        min_dz = min( min_dz , zf(i,j,k+1)-zf(i,j,k) )
        max_dz = max( max_dz , zf(i,j,k+1)-zf(i,j,k) )
      enddo
      enddo
      enddo

#ifdef MPI
      var=0.0
      call MPI_ALLREDUCE(min_dx,var,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_dx=var
      var=0.0
      call MPI_ALLREDUCE(min_dy,var,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_dy=var
      var=0.0
      call MPI_ALLREDUCE(min_dz,var,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_dz=var
      var=0.0
      call MPI_ALLREDUCE(max_dx,var,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_dx=var
      var=0.0
      call MPI_ALLREDUCE(max_dy,var,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_dy=var
      var=0.0
      call MPI_ALLREDUCE(max_dz,var,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_dz=var
#endif

      if(dowr) write(outfile,*) '  min_dx = ',min_dx
      if(dowr) write(outfile,*) '  max_dx = ',max_dx
      if(dowr) write(outfile,*) '  min_dy = ',min_dy
      if(dowr) write(outfile,*) '  max_dy = ',max_dy
      if(dowr) write(outfile,*) '  min_dz = ',min_dz
      if(dowr) write(outfile,*) '  max_dz = ',max_dz
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  max level for 3d output files:

      ! default: entire domain
      maxk = nk

      doit = .false.
      IF( doit )THEN

        max_z_out = 8000.0   ! max depth of 3d writeout

        do k=1,nk+1
          if( sigma(k).le.max_z_out ) maxk = k
        enddo

      ENDIF

      ! check:
      maxk = min( maxk , nk )
      maxk = max( maxk , 1 )

      if(dowr) write(outfile,*) '  maxk,maxzh,maxzf = ',maxk,sigma(maxk),sigmaf(maxk+1)
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  new (cm1r16) arrays for vertical interpolation:

      do k=1,nk+1
      do j=jb,je
      do i=ib,ie
        cc2(i,j,k)=(zf(i,j,k)-zh(i,j,k-1))/(zh(i,j,k)-zh(i,j,k-1))
      enddo
      enddo
      enddo

      call bcs(cc2)
#ifdef MPI
      call comm_all_s(cc2,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,  &
                          n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,reqs_s)
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  k,c1,c2,zhm1,zf,zh:'
      do k=1,nk+1
      do j=jb,je
      do i=ib,ie
        cc1(i,j,k)=1.0-cc2(i,j,k)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,141) k,cc1(i,j,k),cc2(i,j,k),zh(i,j,k-1),zf(i,j,k),zh(i,j,k)
141     format(3x,i4,2(3x,f7.4),3(3x,f8.2))
      enddo
      enddo
      enddo
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Specify coefficient for Rayleigh damper in vertical

      if( (irdamp.ge.1).and.(zd.lt.maxz) )then

        IF( zd.lt.(0.5*maxz) )THEN
          if(myid.eq.0)then
          print *
          print *,'  Warning:  with these settings, Rayleigh damping would  '
          print *,'  be applied over MORE than half the domain '
          print *
          print *,'  zd,maxz = ',zd,maxz
          print *
          print *,'   stopping model .... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        ENDIF

        do j=jb,je
        do i=ib,ie
          do k=1,nk
            if(zh(i,j,k).gt.zd)then
            tauh(i,j,k)=0.5*(1.0-cos(pi*(zh(i,j,k)-zd)/(zf(i,j,nk+1)-zd)))
            taus(i,j,k)=tauh(i,j,k)
            endif
          enddo
          enddo
        enddo
 
        do j=jb,je
        do i=ib,ie
          do k=1,nk+1
            if(zf(i,j,k).gt.zd)then
            tauf(i,j,k)=0.5*(1.0-cos(pi*(zf(i,j,k)-zd)/(zf(i,j,nk+1)-zd)))
            endif
          enddo
          enddo
        enddo

      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,201)
201   format(25x,'  ------ zf,tauf,   zh,tauh -----')
      do k=1,nk
        if(dowr) write(outfile,*) k,zf(1,1,k),tauf(1,1,k),zh(1,1,k),tauh(1,1,k)
      enddo
      if(dowr) write(outfile,*) nk+1,zf(1,1,nk+1),tauf(1,1,nk+1)
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Rayleigh damping near lateral boundaries:

      IF(hrdamp.ge.1)THEN

        IF( nx.gt.1 )THEN
        IF( xhd.gt.(0.5*(maxx-minx)) )THEN
          if(myid.eq.0)then
          print *
          print *,'  Warning:  with these settings, Rayleigh damping would  '
          print *,'  be applied over MORE than half the domain '
          print *
          print *,'  xhd,minx,maxx = ',xhd,minx,maxx
          print *
          print *,'   stopping model .... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        ENDIF
        ENDIF
        IF( ny.gt.1 )THEN
        IF( xhd.gt.(0.5*(maxy-miny)) )THEN
          if(myid.eq.0)then
          print *
          print *,'  Warning:  with these settings, Rayleigh damping would  '
          print *,'  be applied over MORE than half the domain '
          print *
          print *,'  xhd,miny,maxy = ',xhd,miny,maxy
          print *
          print *,'   stopping model .... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        ENDIF
        ENDIF

        ! 170803:  If you do not want/need the horizontal Rayleigh damper
        !          on certain sides of the domain, then set the 
        !          appropriate variable below to .false.

        hrdamp_west   =  .true.
        hrdamp_east   =  .true.
        hrdamp_south  =  .true.
        hrdamp_north  =  .true.

        ! cm1r19.6:  hrdamp acts only on w

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ! skip this section of code for 2d simulations:
          IF(nx.gt.1)THEN
          IF( hrdamp_west )THEN
            ! west boundary:
            IF( axisymm.ne.1 )THEN
              x1 = (xhd+minx)-xh(i)
              if( x1.gt.0.0 )then
                tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*x1/xhd)) )
              endif
            ENDIF
          ENDIF
          IF( hrdamp_east )THEN
            ! east boundary:
            x2 = xh(i)-(maxx-xhd)
            if( x2.gt.0.0 )then
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*x2/xhd)) )
            endif
          ENDIF
          ENDIF
          ! skip this section of code for 2d simulations:
          IF(ny.gt.1)THEN
          IF( hrdamp_south )THEN
            ! south boundary:
            y1 = (xhd+miny)-yh(j)
            if( y1.gt.0.0 )then
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*y1/xhd)) )
            endif
          ENDIF
          IF( hrdamp_north )THEN
            ! north boundary:
            y2 = yh(j)-(maxy-xhd)
            if( y2.gt.0.0 )then
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*y2/xhd)) )
            endif
          ENDIF
          ENDIF
        enddo
        enddo
        enddo

        IF( nx.gt.1 )THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ------ tauf for horizontal Rayleigh damping (W-to-E) -----'
          j = nint(0.5*float(nj))
          j = max( j , 1 )
          j = min( j , nj )
          if(dowr) write(outfile,*) '  i,xh,tauf:'
          do i=0,ni+1
            if(dowr) write(outfile,*) i,xh(i),tauf(i,j,1)
          enddo
          if(dowr) write(outfile,*)
        ENDIF

        IF( ny.gt.1 )THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ------ tauf for horizontal Rayleigh damping (S-to-N) -----'
          i = nint(0.5*float(ni))
          i = max( i , 1 )
          i = min( i , ni )
          if(dowr) write(outfile,*) '  j,yh,tauf:'
          do j=0,nj+1
            if(dowr) write(outfile,*) j,yh(j),tauf(i,j,1)
          enddo
          if(dowr) write(outfile,*)
        ENDIF

      ENDIF

!--------------------------------------------------------------
!  vertically implicit turbulent diffusion:

      ! Set vialpha:
      !      0.0 = forward-in-time (unstable if K dt / (dz^2) > 0.5)
      !      0.5 = centered-in-time (Crank-Nicholson) (stable but oscillatory)
      !      1.0 = backward-in-time (stable)
!      vialpha = 1.0

      ! Do not change this:
!      vibeta  = 1.0 - vialpha

!      NOTE:  these are now set in constants.F file

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  vialpha,vibeta = ',vialpha,vibeta
        if(dowr) write(outfile,*)

!--------------------------------------------------------------

        dot2p = .false.

        ntwk = 1
        ib2pt = 1
        ie2pt = 1
        jb2pt = 1
        je2pt = 1
        kb2pt = 1
        ke2pt = 1

        dotimeavg = .false.
        ntim = 1
        rtim = 0.0
        ntavr = 0
        nsfctavr = 1

        ibta = 1
        ieta = 1
        jbta = 1
        jeta = 1
        kbta = 1
        keta = 1

        utav = 0
        vtav = 0
        wtav = 0
        ttav = 0
        qtav = 0
        etav = 0
        uutav = 0
        vvtav = 0

!--------------------------------------------------------------
!  200412:  2-part turbulence model

        if( .not. idoles ) dot2p = .false.

        if( idoles .and. ( sgsmodel.eq.3 .or. sgsmodel.eq.4 ) ) dot2p = .true.

        IF( idoles .and. dot2p )THEN

          gam_cp     =    150.0        ! centerpoint of transition zone (m ASL)

          gam_dr     =     25.0        ! decay rate parameter (m)

          if( testcase.eq.11 )then
            ! SAS:
            gam_cp = 75.0
            gam_dr = 12.5
          endif

          do k=1,nk+1
            gamwall(k) = 1.0-0.5*( 1.0+tanh( (zf(1,1,k)-gam_cp)/gam_dr ) )
            if( gamwall(k).gt.0.01 ) ntwk = k
!!!            if( gamwall(k).gt.0.999 ) gamwall(k) = 1.0
            gamk(k) = 1.0
!!!            gamk(k) = 0.4 + 0.6*tanh( zf(1,1,k)/(0.5*gam_cp) )
!!!            if(    gamk(k).gt.0.999 )    gamk(k) = 1.0
          enddo

          ntwk = max( 2 , ntwk )

          ib2pt = ib
          ie2pt = ie
          jb2pt = jb
          je2pt = je
          kb2pt = 1
          ke2pt = ntwk

        !................

          if( t2p_avg.le.0 .or. t2p_avg.ge.3 )then
            print *,'  invalid value for t2p_avg: ',t2p_avg
            call stopcm1
          endif

          dotimeavg = .false.

          if( t2p_avg.eq.2 )then
            dotimeavg = .true.
          endif

        !................

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Including two-part sgs turbulence model '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '      ntwk,ztwk   = ',ntwk,zf(1,1,ntwk)
          if(dowr) write(outfile,*) '      ib2pt,ie2pt = ',ib2pt,ie2pt
          if(dowr) write(outfile,*) '      jb2pt,je2pt = ',jb2pt,je2pt
          if(dowr) write(outfile,*) '      kb2pt,ke2pt = ',kb2pt,ke2pt
          if(dowr) write(outfile,*)
          if( t2p_avg.eq.1 )then
            if(dowr) write(outfile,*) '    using horizontal average '
          elseif( t2p_avg.eq.2 )then
            if(dowr) write(outfile,*) '    using time average '
          endif
          if(dowr) write(outfile,*)

          if(dowr) write(outfile,*) '      k,zf,gamwall,gamk: '
          do k=1,nk+1
            if(dowr) write(outfile,*) k,zf(1,1,k),gamwall(k),gamk(k)
          enddo
          if(dowr) write(outfile,*)

        ELSE

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '      dot2p          = ',dot2p
          if(dowr) write(outfile,*) '      ntwk           = ',ntwk
          if(dowr) write(outfile,*) '      ib2pt,ie2pt    = ',ib2pt,ie2pt
          if(dowr) write(outfile,*) '      jb2pt,je2pt    = ',jb2pt,je2pt
          if(dowr) write(outfile,*) '      kb2pt,ke2pt    = ',kb2pt,ke2pt
          if(dowr) write(outfile,*)

        ENDIF

!--------------------------------------------------------------
!  170803:
!  Coriolis array:

      ! By default, assume an f-plane:
      f2d = fcor


      if( betaplane.eq.1 )then

        ! get central latitude from Coriolis parameter:
        phi = asin( fcor/(2.0*omega) )

        ! get gradient of f at that latitude:
        beta0 = 2.0*omega*cos(phi)/earth_radius
        if( myid.eq.0 ) print *,'  phi(degr),omega,beta0 = ',phi*180.0/pi,omega,beta0

        ! simple "midlatitude" beta-plane formulation:
        do j=jb,je
        do i=ib,ie
          f2d(i,j) = fcor + beta0*(yh(j)-centery)
        enddo
        enddo

      endif

    if( lspgrad.eq.3 )then
      ! for lspgrad = 3
      hurr_angle = dble(hurr_rotate) * ( pi_dp / 180.0d0 )
    endif

!--------------------------------------------------------------

      dt = dtl
      dtlast = dt

      deallocate( xfdp )
      deallocate( yfdp )

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Get 2nd-order extrapolation coefficients (Fornberg 1988)   ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  DO nn=1,4

    IF( nn.eq.1 )THEN
      x0 = sigmaf(1)/(sigmaf(2)-sigmaf(1))
      alpha(0) = sigma(1)/(sigmaf(2)-sigmaf(1))
      alpha(1) = sigma(2)/(sigmaf(2)-sigmaf(1))
      alpha(2) = sigma(3)/(sigmaf(2)-sigmaf(1))
      alpha(3) = sigma(4)/(sigmaf(2)-sigmaf(1))
    ELSEIF( nn.eq.2 )THEN
      x0 = sigmaf(nk+1)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(0) = sigma(nk  )/(sigmaf(nk+1)-sigmaf(nk))
      alpha(1) = sigma(nk-1)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(2) = sigma(nk-2)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(3) = sigma(nk-3)/(sigmaf(nk+1)-sigmaf(nk))
    ELSEIF( nn.eq.3 )THEN
      x0 = sigmaf(1)/(sigmaf(2)-sigmaf(1))
      alpha(0) = sigmaf(2)/(sigmaf(2)-sigmaf(1))
      alpha(1) = sigmaf(3)/(sigmaf(2)-sigmaf(1))
      alpha(2) = sigmaf(4)/(sigmaf(2)-sigmaf(1))
      alpha(3) = sigmaf(5)/(sigmaf(2)-sigmaf(1))
    ELSEIF( nn.eq.4 )THEN
      x0 = sigmaf(nk+1)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(0) = sigmaf(nk  )/(sigmaf(nk+1)-sigmaf(nk))
      alpha(1) = sigmaf(nk-1)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(2) = sigmaf(nk-2)/(sigmaf(nk+1)-sigmaf(nk))
      alpha(3) = sigmaf(nk-3)/(sigmaf(nk+1)-sigmaf(nk))
    ENDIF

      delta = 0.0

      delta(0,0,0) = 1.0
      b1 = 1.0

      do n = 1,bign
        b2 = 1.0
        do nu = 0,n-1
          b3 = alpha(n)-alpha(nu)
          b2 = b2*b3
          if( n.le.bigm ) delta(n-1,n,nu) = 0.0
          do m = 0,min(n,bigm)
            delta(n,m,nu) = ( (alpha(n)-x0)*delta(n-1,m,nu) - m*delta(n-1,m-1,nu) )/b3
          enddo
        enddo
        do m = 0,min(n,bigm)
          delta(n,m,n) = (b1/b2)*( m*delta(n-1,m-1,n-1) - (alpha(n-1)-x0)*delta(n-1,m,n-1) )
        enddo
        b1 = b2
      enddo

    IF( nn.eq.1 )THEN
      cgs1 = delta(2,0,0)
      cgs2 = delta(2,0,1)
      cgs3 = delta(2,0,2)
      var = cgs1*sigma(1)+cgs2*sigma(2)+cgs3*sigma(3)
      dgs1 = delta(2,1,0)
      dgs2 = delta(2,1,1)
      dgs3 = delta(2,1,2)
    ELSEIF( nn.eq.2 )THEN
      cgt1 = delta(2,0,0)
      cgt2 = delta(2,0,1)
      cgt3 = delta(2,0,2)
      var = cgt1*sigma(nk)+cgt2*sigma(nk-1)+cgt3*sigma(nk-2)
      dgt1 = delta(2,1,0)
      dgt2 = delta(2,1,1)
      dgt3 = delta(2,1,2)
    ELSEIF( nn.eq.3 )THEN
      wbe1 = delta(2,0,0)
      wbe2 = delta(2,0,1)
      wbe3 = delta(2,0,2)
      var = wbe1*sigmaf(2)+wbe2*sigmaf(3)+wbe3*sigmaf(4)
    ELSEIF( nn.eq.4 )THEN
      wte1 = delta(2,0,0)
      wte2 = delta(2,0,1)
      wte3 = delta(2,0,2)
      var = wte1*sigmaf(nk)+wte2*sigmaf(nk-1)+wte3*sigmaf(nk-2)
    ENDIF

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ---------------------------------- '
      if(dowr) write(outfile,*) '  nn = ',nn
      if(dowr) write(outfile,*) '  x0,alpha,delta,sum,predict,maxz:'
      if(dowr) write(outfile,*) sngl(x0)
      if(dowr) write(outfile,*) sngl(alpha(0)),sngl(alpha(1)),sngl(alpha(2)),sngl(alpha(3))
      if(dowr) write(outfile,*) sngl(delta(2,0,0)),sngl(delta(2,0,1)),sngl(delta(2,0,2)),sngl(delta(2,0,3))
      if(dowr) write(outfile,*) sngl(delta(2,0,0)+delta(2,0,1)+delta(2,0,2)+delta(2,0,3)),var,maxz
      if(dowr) write(outfile,*)

      IF( nn.eq.1 )then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  cgs1,cgs2,cgs3 = ',cgs1,cgs2,cgs3
        if(dowr) write(outfile,*) '  sum            = ',cgs1+cgs2+cgs3
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  dgs1,dgs2,dgs3 = ',dgs1,dgs2,dgs3
        if(dowr) write(outfile,*) '  sum            = ',dgs1+dgs2+dgs3
        if(dowr) write(outfile,*)
      ELSEIF( nn.eq.2 )THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  cgt1,cgt2,cgt3 = ',cgt1,cgt2,cgt3
        if(dowr) write(outfile,*) '  sum            = ',cgt1+cgt2+cgt3
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  dgt1,dgt2,dgt3 = ',dgt1,dgt2,dgt3
        if(dowr) write(outfile,*) '  sum            = ',dgt1+dgt2+dgt3
        if(dowr) write(outfile,*)
      ELSEIF( nn.eq.3 )THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  wbe1,wbe2,wbe3 = ',wbe1,wbe2,wbe3
        if(dowr) write(outfile,*) '  sum            = ',wbe1+wbe2+wbe3
        if(dowr) write(outfile,*)
      ELSEIF( nn.eq.4 )THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  wte1,wte2,wte3 = ',wte1,wte2,wte3
        if(dowr) write(outfile,*) '  sum            = ',wte1+wte2+wte3
        if(dowr) write(outfile,*)
      ENDIF

      if(dowr) write(outfile,*) '  ---------------------------------- '

  ENDDO

!--------------------------------------------------------------
!  cm1r18:  Set ghost points for zh
!  [Note:  since cm1r17, the array index (i,j,0) means the surface]
!     (upper/lower ghost points are used by parcel subroutines only)

    DO j=jb,je
    DO i=ib,ie
      zh(i,j,0) = zf(i,j,1)
      zh(i,j,nk+1) = zf(i,j,nk+1)
    ENDDO
    ENDDO

!--------------------------------------------------------------
!  moved here so zh is available:

      IF( ptype.eq.3 )THEN

          !----- initialize the Thompson scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling thompson_init'
          if(dowr) write(outfile,*) '(this can take several minutes ... please be patient)'

          call   thompson_init(hgt=zh, dowr=dowr, nt_c_in=nt_c,              &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,   &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,     &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk )

          if(dowr) write(outfile,*) 'Done with thompson_init'
          if(dowr) write(outfile,*)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*) '  ... using Thompson microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '         assuming constant cloud droplet concentration' 
          if(dowr) write(outfile,*) '         Nt_c   = ',Nt_c
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ----------------------------------------------- '
          if(dowr) write(outfile,*)

      ENDIF

!--------------------------------------------------------------

      dothis = .false.


      IF( dothis .and. myid.eq.0 )THEN

        ! write out grids to text files:

        open(unit=21,file='input_grid_x')
        do i=1,nx
          write(21,*) 0.5*(xfref(i)+xfref(i+1))
        enddo
        close(unit=21)

        open(unit=21,file='input_grid_y')
        do j=1,ny
          write(21,*) 0.5*(yfref(j)+yfref(j+1))
        enddo
        close(unit=21)

        open(unit=21,file='input_grid_z')
        do k=1,nk
          write(21,*) zh(1,1,k)
        enddo
        close(unit=21)

      ENDIF

!--------------------------------------------------------------

      dtu0 = 0.0

      do j=1,nj
      do i=1,ni+1
        dtu0(i,j) = 1.0
      enddo
      enddo

      if( ibw.eq.1 )then
        do j=1,nj
          dtu0(1,j) = 0.0
        enddo
      endif
      if( ibe.eq.1 )then
        do j=1,nj
          dtu0(ni+1,j) = 0.0
        enddo
      endif

    !.........................

      dtv0 = 0.0

      do j=1,nj+1
      do i=1,ni
        dtv0(i,j) = 1.0
      enddo
      enddo

      if( ibs.eq.1 )then
        do i=1,ni
          dtv0(i,1) = 0.0
        enddo
      endif
      if( ibn.eq.1 )then
        do i=1,ni
          dtv0(i,nj+1) = 0.0
        enddo
      endif

!------------------------------------------------------------------

    IF( doazimavg )THEN
      !----------
      ! cm1r20.1: set ddr equal to minimum grid spacing
      ddr = 1.0e30
      ddr = min( ddr , min_dx )
      ddr = min( ddr , min_dy )
      !----------
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'doazimavg    = ',doazimavg
      if(dowr) write(outfile,*) 'azimavgfrq   = ',azimavgfrq
      if(dowr) write(outfile,*) 'ddr          = ',ddr
      if(dowr) write(outfile,*) 'rlen         = ',rlen
      if(dowr) write(outfile,*)
    ENDIF

!----------------------------------------------------------------------

    IF( cm1setup.eq.4 )THEN

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'les_subdomain_shape   = ',les_subdomain_shape
      if(dowr) write(outfile,*) 'les_subdomain_xlen    = ',les_subdomain_xlen
      if(dowr) write(outfile,*) 'les_subdomain_ylen    = ',les_subdomain_ylen
      if(dowr) write(outfile,*) 'les_subdomain_dlen    = ',les_subdomain_dlen
      if(dowr) write(outfile,*) 'les_subdomain_trnslen = ',les_subdomain_trnslen
      if(dowr) write(outfile,*)

    ENDIF

!----------------------------------------------------------------------
!  Eddy recycling code:


      call eddy_recycle_setup


      urecy = 1
      vrecy = 2
      wrecy = 3
      xrecy = 4
      nrecy = 4

      tscale0 = 10.0*dtl


      do_recycle = .false.
      if( do_recycle_w ) do_recycle = .true.
      if( do_recycle_e ) do_recycle = .true.
      if( do_recycle_s ) do_recycle = .true.
      if( do_recycle_n ) do_recycle = .true.

      irecywe = 1
      jrecywe = 1
      irecysn = 1
      jrecysn = 1
      krecy = 1
      nrecy = 1

      urecy = 0
      vrecy = 0
      wrecy = 0
      trecy = 0
      qrecy = 0
      erecy = 0
      xrecy = 0

      ! dummy checks:
      if( nx.le.3 .and. ny.le.3 )then
        ! do recycling for single-column modeling
        do_recycle   = .false.
        do_recycle_w = .false.
        do_recycle_e = .false.
        do_recycle_s = .false.
        do_recycle_n = .false.
      endif
      if( nx.le.3 )then
        ! do recycling on w/e for very small domains in w-e direction
        do_recycle_w = .false.
        do_recycle_e = .false.
      endif
      if( ny.le.3 )then
        ! do recycling on s/n for very small domains in s-n direction
        do_recycle_s = .false.
        do_recycle_n = .false.
      endif


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  do_recycle   = ',do_recycle
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  do_recycle_w = ',do_recycle_w
      if(dowr) write(outfile,*) '  do_recycle_s = ',do_recycle_s
      if(dowr) write(outfile,*) '  do_recycle_e = ',do_recycle_e
      if(dowr) write(outfile,*) '  do_recycle_n = ',do_recycle_n
      if(dowr) write(outfile,*)

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ifrecycle:  &
      IF( do_recycle )THEN

        ! dummy check:
        if( cm1setup.eq.2 )then
          if(myid.eq.0)then
          print *
          print *,'  cm1setup = ',cm1setup
          print *
          print *,'  cm1setup cannot be 2 when using eddy recycling '
          print *
          print *,'   stopping model .... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif


        recy_width  =   min_dx * recycle_width_dx
        recy_depth  =   recycle_depth_m


        if( (recycle_inj_loc_m + recy_width) .ge. recycle_cap_loc_m )then
          if(myid.eq.0)then
          print *
          print *,'  recycle_cap_loc_m needs to be larger '
          print *
          print *,'  recycle_inj_loc_m + recy_width = ',recycle_inj_loc_m + recy_width
          print *,'  recycle_cap_loc_m              = ',recycle_cap_loc_m
          print *
          print *,'   stopping model .... '
          print *
          endif
#ifdef MPI
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
          call stopcm1
        endif


        IF( cm1setup.eq.4 )THEN

          ! LES within mesoscale model:
          ! injection/capture regions are based on fine-mesh part of domain

          recy_cap_w  =  centerx - 0.5*les_subdomain_xlen + recycle_cap_loc_m
          recy_cap_e  =  centerx + 0.5*les_subdomain_xlen - recycle_cap_loc_m
          recy_cap_s  =  centerx - 0.5*les_subdomain_ylen + recycle_cap_loc_m
          recy_cap_n  =  centerx + 0.5*les_subdomain_ylen - recycle_cap_loc_m

          recy_inj_w  =  centerx - 0.5*les_subdomain_xlen + recycle_inj_loc_m
          recy_inj_e  =  centerx + 0.5*les_subdomain_xlen - recycle_inj_loc_m
          recy_inj_s  =  centerx - 0.5*les_subdomain_ylen + recycle_inj_loc_m
          recy_inj_n  =  centerx + 0.5*les_subdomain_ylen - recycle_inj_loc_m

        ELSE

          ! typical simulation:

          if(myid.eq.0) print *
          if(myid.eq.0) print *,'  minx,maxx = ',minx,maxx
          if(myid.eq.0) print *,'  miny,maxy = ',miny,maxy
          if(myid.eq.0) print *,'  maxz      = ',maxz
          if(myid.eq.0) print *

          recy_cap_w     =  minx + recycle_cap_loc_m
          recy_cap_e     =  maxx - recycle_cap_loc_m
          recy_cap_s     =  miny + recycle_cap_loc_m
          recy_cap_n     =  maxy - recycle_cap_loc_m

          recy_inj_w     =  minx + recycle_inj_loc_m
          recy_inj_e     =  maxx - recycle_inj_loc_m
          recy_inj_s     =  miny + recycle_inj_loc_m
          recy_inj_n     =  maxy - recycle_inj_loc_m

        ENDIF

        if(dowr) write(outfile,*) '  recy_width    = ',recy_width
        if(dowr) write(outfile,*) '  recy_depth    = ',recy_depth
        if(dowr) write(outfile,*)
        if(dowr.and.do_recycle_w) write(outfile,*) '  recy_cap_w    = ',recy_cap_w
        if(dowr.and.do_recycle_s) write(outfile,*) '  recy_cap_s    = ',recy_cap_s
        if(dowr.and.do_recycle_e) write(outfile,*) '  recy_cap_e    = ',recy_cap_e
        if(dowr.and.do_recycle_n) write(outfile,*) '  recy_cap_n    = ',recy_cap_n
        if(dowr) write(outfile,*)
        if(dowr.and.do_recycle_w) write(outfile,*) '  recy_inj_w    = ',recy_inj_w
        if(dowr.and.do_recycle_s) write(outfile,*) '  recy_inj_s    = ',recy_inj_s
        if(dowr.and.do_recycle_e) write(outfile,*) '  recy_inj_e    = ',recy_inj_e
        if(dowr.and.do_recycle_n) write(outfile,*) '  recy_inj_n    = ',recy_inj_n
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  tscale0       = ',tscale0
        if(dowr) write(outfile,*)

        ! 3 components of velocity + theta:
        urecy = 1
        vrecy = 2
        wrecy = 3
        trecy = 4
        nrecy = 4

        if( imoist.eq.1 )then
          nrecy = nrecy+1
          qrecy = nrecy
        endif
        if( iusetke )then
          nrecy = nrecy+1
          erecy = nrecy
        endif

        nrecy = nrecy+1
        xrecy = nrecy

        if(dowr) write(outfile,*) '  nrecy = ',nrecy
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '      urecy = ',urecy
        if(dowr) write(outfile,*) '      vrecy = ',vrecy
        if(dowr) write(outfile,*) '      wrecy = ',wrecy
        if(dowr) write(outfile,*) '      trecy = ',trecy
        if(dowr) write(outfile,*) '      qrecy = ',qrecy
        if(dowr) write(outfile,*) '      erecy = ',erecy
        if(dowr) write(outfile,*) '      xrecy = ',xrecy
        if(dowr) write(outfile,*)

        !-----------------------------------------------------

        recy_depth  =  min( recy_depth , maxz )

        do k=1,nk
          if( zh(1,1,k).le.recy_depth ) krecy = k
        enddo

        ! currently uses time-average for reference state
        dotimeavg = .true.
        call setup_timeavg(nlevels=krecy)

          if( utav.le.0 )then
            ntavr = ntavr+1
            utav = ntavr
          endif
          if( vtav.le.0 )then
            ntavr = ntavr+1
            vtav = ntavr
          endif
          if( wtav.le.0 )then
            ntavr = ntavr+1
            wtav = ntavr
          endif
          if( ttav.le.0 )then
            ntavr = ntavr+1
            ttav = ntavr
          endif
          if( imoist.eq.1 .and. qtav.le.0 )then
            ntavr = ntavr+1
            qtav = ntavr
          endif
          if( idoles .and. iusetke .and. etav.le.0 )then
            ntavr = ntavr+1
            etav = ntavr
          endif

        !-----------------------------------------------------
        ! get indices for eddy recycling arrays:

          call eddy_recycling_indices(xh,xf,yh,yf,xhref,xfref,yhref,yfref)

        !-----------------------------------------------------

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '      irecywe = ',irecywe
        if(dowr) write(outfile,*) '      jrecywe = ',jrecywe
        if(dowr) write(outfile,*) '      irecysn = ',irecysn
        if(dowr) write(outfile,*) '      jrecysn = ',jrecysn
        if(dowr) write(outfile,*) '      krecy   = ',krecy
        if(dowr) write(outfile,*) '      nrecy   = ',nrecy
        if(dowr) write(outfile,*)

      ENDIF  ifrecycle

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!--------------------------------------------------------------

      ! for time averaging:
      IF( dotimeavg )THEN

        call setup_timeavg(nlevels=ntwk)

        nsfctavr = 3

        if( utav.le.0 )then
          ntavr = ntavr+1
          utav = ntavr
        endif

        if( vtav.le.0 )then
          ntavr = ntavr+1
          vtav = ntavr
        endif

        if( dot2p .and. uutav.le.0 )then
          ntavr = ntavr+1
          uutav = ntavr
        endif

        if( dot2p .and. vvtav.le.0 )then
          ntavr = ntavr+1
          vvtav = ntavr
        endif

      ENDIF

      ntavr = max( ntavr , 1 )
      nsfctavr = max( nsfctavr , 1 )

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  dotimeavg = ',dotimeavg
      if(dowr) write(outfile,*)

      IF( dotimeavg )THEN
        if(dowr) write(outfile,*) '  ... using time averages ... '
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '      ibta,ieta      = ',ibta,ieta
        if(dowr) write(outfile,*) '      jbta,jeta      = ',jbta,jeta
        if(dowr) write(outfile,*) '      kbta,keta      = ',kbta,keta
        if(dowr) write(outfile,*) '      ntavr,nsfctavr = ',ntavr,nsfctavr
        if(dowr) write(outfile,*) '      ntim,rtim      = ',ntim,rtim,ntim*rtim
        if(dowr) write(outfile,*) '      utav           = ',utav
        if(dowr) write(outfile,*) '      vtav           = ',vtav
        if(dowr) write(outfile,*) '      wtav           = ',wtav
        if(dowr) write(outfile,*) '      ttav           = ',ttav
        if(dowr) write(outfile,*) '      qtav           = ',qtav
        if(dowr) write(outfile,*) '      etav           = ',etav
        if(dowr) write(outfile,*) '      uutav          = ',uutav
        if(dowr) write(outfile,*) '      vvtav          = ',vvtav
        if(dowr) write(outfile,*)
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '      ibta,ieta      = ',ibta,ieta
        if(dowr) write(outfile,*) '      jbta,jeta      = ',jbta,jeta
        if(dowr) write(outfile,*) '      kbta,keta      = ',kbta,keta
        if(dowr) write(outfile,*) '      ntavr,nsfctavr = ',ntavr,nsfctavr
        if(dowr) write(outfile,*)
      ENDIF

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  do_lsnudge         = ',do_lsnudge
      if(dowr) write(outfile,*)
      IF( do_lsnudge )THEN
        if(dowr) write(outfile,*) '  do_lsnudge_u       = ',do_lsnudge_u
        if(dowr) write(outfile,*) '  do_lsnudge_v       = ',do_lsnudge_v
        if(dowr) write(outfile,*) '  do_lsnudge_th      = ',do_lsnudge_th
        if(dowr) write(outfile,*) '  do_lsnudge_qv      = ',do_lsnudge_qv
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  lsnudge_tau        = ',lsnudge_tau
        if(dowr) write(outfile,*) '  lsnudge_start      = ',lsnudge_start
        if(dowr) write(outfile,*) '  lsnudge_end        = ',lsnudge_end
        if(dowr) write(outfile,*) '  lsnudge_ramp_time  = ',lsnudge_ramp_time
        if(dowr) write(outfile,*)
      ENDIF

!--------------------------------------------------------------

    IF( myid.eq.0 )THEN
      open(unit=21,file='cm1_config.txt')
      write(21,*)
      write(21,*) '  CM1 version: ',trim(cm1version)
      write(21,*)
      write(21,*) '  nx                        = ',nx
      write(21,*) '  ny                        = ',ny
      write(21,*) '  nz                        = ',nz
      write(21,*)
      write(21,*) '  imoist                    = ',imoist
      write(21,*) '  sgsmodel                  = ',sgsmodel
      write(21,*) '  tconfig                   = ',tconfig
      write(21,*) '  bcturbs                   = ',bcturbs
      write(21,*) '  ptype                     = ',ptype
      IF ( ptype.eq.26 .or. ptype.eq.27  ) THEN
      write(21,*) '  nssl_3moment              = ',nssl_3moment
      ENDIF
      write(21,*) '  wbc                       = ',wbc
      write(21,*) '  ebc                       = ',ebc
      write(21,*) '  sbc                       = ',sbc
      write(21,*) '  nbc                       = ',nbc
      write(21,*) '  bbc                       = ',bbc
      write(21,*) '  tbc                       = ',tbc
      write(21,*) '  iorigin                   = ',iorigin
      write(21,*) '  axisymm                   = ',axisymm
      write(21,*) '  iptra                     = ',iptra
      write(21,*) '  npt                       = ',npt
      write(21,*)
      write(21,*) '  fcor                      = ',fcor
      write(21,*)
      write(21,*) '  radopt                    = ',radopt
    if( radopt.ge.1 )then
      write(21,*) '  dtrad                     = ',dtrad
      write(21,*) '  ctrlat                    = ',ctrlat
      write(21,*) '  ctrlon                    = ',ctrlon
      write(21,*) '  year                      = ',year
      write(21,*) '  month                     = ',month
      write(21,*) '  day                       = ',day
      write(21,*) '  hour                      = ',hour
      write(21,*) '  minute                    = ',minute
      write(21,*) '  second                    = ',second
    endif
      write(21,*)
      write(21,*) '  sfcmodel                  = ',sfcmodel
      write(21,*) '  oceanmodel                = ',oceanmodel
      write(21,*) '  ipbl                      = ',ipbl
      write(21,*)
      write(21,*) '  iice                      = ',iice
      write(21,*) '  idm                       = ',idm
      write(21,*) '  idmplus                   = ',idmplus
      write(21,*) '  numq                      = ',numq
      write(21,*) '  nql1                      = ',nql1
      write(21,*) '  nql2                      = ',nql2
      write(21,*) '  nqs1                      = ',nqs1
      write(21,*) '  nqs2                      = ',nqs2
      write(21,*) '  nnc1                      = ',nnc1
      write(21,*) '  nnc2                      = ',nnc2
      write(21,*) '  nzl1                      = ',nzl1
      write(21,*) '  nzl2                      = ',nzl2
      write(21,*) '  nvl1                      = ',nvl1
      write(21,*) '  nvl2                      = ',nvl2
      write(21,*)
      write(21,*) '  c_m                       = ',c_m
      write(21,*) '  c_e1                      = ',c_e1
      write(21,*) '  c_e2                      = ',c_e2
      write(21,*) '  c_s                       = ',c_s
      write(21,*)
      write(21,*) '  cgs1                      = ',cgs1
      write(21,*) '  cgs2                      = ',cgs2
      write(21,*) '  cgs3                      = ',cgs3
      write(21,*)
      write(21,*) '  dgs1                      = ',dgs1
      write(21,*) '  dgs2                      = ',dgs2
      write(21,*) '  dgs3                      = ',dgs3
      write(21,*)
      write(21,*) '  cgt1                      = ',cgt1
      write(21,*) '  cgt2                      = ',cgt2
      write(21,*) '  cgt3                      = ',cgt3
      write(21,*)
      write(21,*) '  dgt1                      = ',dgt1
      write(21,*) '  dgt2                      = ',dgt2
      write(21,*) '  dgt3                      = ',dgt3
      write(21,*)
      close(unit=21)
    ENDIF

!--------------------------------------------------------------

      if(dowr) write(outfile,*) 'Leaving PARAM'

      return

    !--------------------------------------------------------------

8000  print *
      print *,'  8000: error opening namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8001  print *
      print *,'  8001: error reading param1 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8002  print *
      print *,'  8002: error reading param2 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8003  print *
      print *,'  8003: error reading param3 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8004  print *
      print *,'  8004: error reading param4 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8005  print *
      print *,'  8005: error reading param5 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8006  print *
      print *,'  8006: error reading param6 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8007  print *
      print *,'  8007: error reading param7 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8008  print *
      print *,'  8008: error reading param8 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8009  print *
      print *,'  8009: error reading param9 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8010  print *
      print *,'  8010: error reading param10 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8011  print *
      print *,'  8011: error reading param11 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8012  print *
      print *,'  8012: error reading param12 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8013  print *
      print *,'  8013: error reading param13 section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

8051  print *
      print *,'  8051: error reading nssl2mom_params section of namelist.input '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9011  print *
      print *,'  9011: error opening input_grid_z file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9012  print *
      print *,'  9012: error reading input_grid_z file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9013  print *
      print *,'  9013: error reading input_grid_x file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9014  print *
      print *,'  9014: error reading input_grid_y file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9015  print *
      print *,'  9015: error opening input_grid_x file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

9016  print *
      print *,'  9016: error opening input_grid_y file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

    !--------------------------------------------------------------

      end subroutine param


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine wenocheck
      use input
      implicit none

      IF( advwenos.ge.1 )THEN
        if( weno_order.eq.3 )then
          hadvordrs = 3
          vadvordrs = 3
        endif
        if( weno_order.eq.5 )then
          hadvordrs = 5
          vadvordrs = 5
        endif
        if( weno_order.eq.7 )then
          hadvordrs = 7
          vadvordrs = 7
        endif
        if( weno_order.eq.9 )then
          hadvordrs = 9
          vadvordrs = 9
        endif
      ENDIF
      IF( advwenov.ge.1 )THEN
        if( weno_order.eq.3 )then
          hadvordrv = 3
          vadvordrv = 3
        endif
        if( weno_order.eq.5 )then
          hadvordrv = 5
          vadvordrv = 5
        endif
        if( weno_order.eq.7 )then
          hadvordrv = 7
          vadvordrv = 7
        endif
        if( weno_order.eq.9 )then
          hadvordrv = 9
          vadvordrv = 9
        endif
      ENDIF

      end subroutine wenocheck


!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE ZGRID  <<<<<<<<<<<<<<<<<<<<<<<<<<< !
!
!-------------------------------------------------------------------------------
! Creates the vertical grid using a geometric stretch 
!
! Option - set nbndlyr > 0 input deck.  A suggested value is ~ nz/4 to create a 
! layer having constant resolution near the surface.
!
!-------------------------------------------------------------------------------

 SUBROUTINE ZGRID(dz, dz_stretch, nz, nbndlyr, gzc, gze, ng, dzmax,ztopstr,rtop,dzmaxtop,myid)

   implicit none


   real dz, dz_stretch
   real  dzmax,ztopstr,rtop,dzmaxtop
   integer, intent(in) :: nbndlyr, nz, ng, myid
   double precision gzc(-ng+1:nz+ng), gze(-ng+1:nz+ng)

!   real dzmax
!   parameter( dzmax = 700. )
   integer n, k
   double precision stretch, zx, xmid, fmid, ztop

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

! Use Newton interation to find grid coefficients

   ztop    = dz * (nz-1)
   zx      = 1.0d0
   stretch = 1.0d0

   IF( dz .gt. dz_stretch ) THEN

    DO n = 1,50

     IF( abs(zx) .gt. 1.0e-12 ) THEN
      zx   = zx * 0.5
      xmid = stretch + zx
      fmid = ZHEIGHT(dz_stretch,xmid,nz-1,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop) - ztop
!      write(6,*) 'Stretch: ',n,zx,xmid,fmid,fmid + ztop
      IF( fmid .le. 0.0 ) stretch = xmid
      IF ( fmid .eq. 0.0d0 ) EXIT
     ENDIF

    ENDDO

   ENDIF

   write(6,*)
   IF( stretch .gt. 1.1 ) THEN
    IF ( myid == 0 ) THEN
    write(6,*) 'STRETCH FAC TOO BIG! - NUMERICAL ERRORS WILL BE LARGE'
    write(6,*) 'STRETCH FAC  = ',stretch
    write(6,*) 'INCREASE NZ in namelist or increase dz_bot'
    ENDIF
    call stopcm1
   ELSE
    IF ( myid == 0 ) THEN
    write(6,*) 'ZGRID:  STRETCH FAC  = ',stretch
    write(6,*) 'ZGRID:  DOMAIN  HGT  = ',ZHEIGHT(dz_stretch,stretch,nz-1,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
    ENDIF
   ENDIF

   gze(1) = 0.0
   DO k = 1,nz-1
    gze(k+1) = ZHEIGHT(dz_stretch,stretch,k,dzmax,nbndlyr,ztopstr,rtop,dzmaxtop)
    gzc(k)   = 0.5 * ( gze(k) + gze(k+1) )
!    write(6,*) 'ZGRID: gze,gzc = ',k,gze(k+1),gzc(k)
   ENDDO

   gzc(nz) = 2.*gzc(nz-1) - gzc(nz-2) 


  END SUBROUTINE ZGRID

!--------------------------------------------------------------------------
! FUNCTION ZHEIGHT:  Computes the height of a geometrically stretched grid
!                    with a few wrinkles:  It can have a layer of constant
!                    dz at the bottom 'n1' layers thick, it also limits
!                    the size of dz at the top of the model to be 'dzmax'.

 REAL FUNCTION ZHEIGHT(dzbot,r,nz,dzmax,n1,zctop,ztopr,dzmax2)

  implicit none
  integer nz, n1, k, k2
  integer n2
  real dzbot
  double precision r
  double precision sum
  double precision dznew, dzmaxdp, dzmax2dp
  real dzmax
  real zctop  ! height for upper level stretch
  real ztopr  ! upper level stretch factor
  real dzmax2 ! maximum upper dz
  real dzm

  sum = 0.0d0
  dzmaxdp = dzmax
  dzmax2dp = dzmax2
  
!  zctop = 10000.
!  ztopr = 1.09
!  dzmax2 = 1000. ! 2*dzmax
  
  n2 = 0

  DO k = 1,nz

   IF( k .le. n1 ) THEN
    dznew=dzbot
   ELSE
    k2=k-n1
    dznew = Min(dzbot * r**(k2-1),dzmaxdp)
      IF ( sum .ge. zctop ) THEN
       IF ( n2 .eq. 0 ) dzm = Min(dznew,dzmaxdp)
       n2 = n2 + 1
       dznew = Min(dzm * ztopr**n2, dzmax2)
      ENDIF
   ENDIF
   sum = sum + dznew

  ENDDO

  ZHEIGHT = sum

 RETURN
 END FUNCTION ZHEIGHT

!-------------------------------------------------------------------------------

      subroutine setup_timeavg(nlevels)
      use input
      implicit none

      integer, intent(in) :: nlevels


      !  To reduce memory usage, the time-averaging code is designed as a series 
      !  of smaller time averages.
      !
      !  Total length of time averge is ntim * rtim.  
      !  (so, ntim = 5 with rtime = 60.0s  yields 300-s average


      ! number of time-avgeraging levels:
      ntim  =   5

      ! length of time averaging in each level (seconds):
      rtim  =  60.0


      ibta = ib
      ieta = ie
      jbta = jb
      jeta = je
      kbta = 1
      keta = max( keta , nlevels )

      end subroutine setup_timeavg

!-------------------------------------------------------------------------------

  END MODULE param_module
