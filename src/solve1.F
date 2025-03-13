  MODULE solve1_module

        ! solve1: pre-RK terms
        ! (tendencies held fixed during Runge-Kutta integration and acoustic sub-stepping)

!-----------------------------------------------------------------------------
!
!  CM1 Numerical Model, Release 21.1  (cm1r21.1)
!  24 March 2024
!  https://www2.mmm.ucar.edu/people/bryan/cm1/
!
!  (c)2024 - University Corporation for Atmospheric Research 
!
!-----------------------------------------------------------------------------
!  Quick Index:
!    ua/u3d     = velocity in x-direction (m/s)  (grid-relative)
!    va/v3d     = velocity in y-direction (m/s)  (grid-relative)
!       Note: when imove=1, ground-relative winds are umove+ua, umove+u3d,
!                                                     vmove+va, vmove+v3d.
!    wa/w3d     = velocity in z-direction (m/s)
!    tha/th3d   = perturbation potential temperature (K)
!    ppi/pp3d   = perturbation nondimensional pressure ("Exner function")
!    qa/q3d     = mixing ratios of moisture (kg/kg)
!    tkea/tke3d = SUBGRID turbulence kinetic energy (m^2/s^2)
!    kmh/kmv    = turbulent diffusion coefficients for momentum (m^2/s)
!    khh/khv    = turbulent diffusion coefficients for scalars (m^2/s)
!                 (h = horizontal, v = vertical)
!    prs        = pressure (Pa)
!    rho        = density (kg/m^3)
!
!    th0,pi0,prs0,etc = base-state arrays
!
!    xh         = x (m) at scalar points
!    xf         = x (m) at u points
!    yh         = y (m) at scalar points
!    yf         = y (m) at v points
!    zh         = z (m above sea level) of scalar points (aka, "half levels")
!    zf         = z (m above sea level) of w points (aka, "full levels")
!
!    For the axisymmetric model (axisymm=1), xh and xf are radius (m).
!
!  See "The governing equations for CM1" for more details:
!        https://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!-----------------------------------------------------------------------------
!  Some notes:
!
!  - Upon entering solve, the arrays ending in "a" (eg, ua,wa,tha,qa,etc)
!    are equivalent to the arrays ending in "3d" (eg, u3d,w3d,th3d,q3d,etc).
!  - The purpose of solve is to update the variables from time "t" to time
!    "t+dt".  Values at time "t+dt" are stored in the "3d" arrays.
!  - The "ghost zones" (boundaries beyond the computational subdomain) are
!    filled out completely (3 rows/columns) for the "3d" arrays.  To save 
!    unnecessary computations, starting with cm1r15 the "ghost zones" of 
!    the "a" arrays are only filled out to 1 row/column.  Hence, if you 
!    need to do calculations that use a large stencil, you must use the 
!    "3d" arrays (not the "a") arrays.
!  - Arrays named "ten" store tendencies.  Those ending "ten1" store
!    pre-RK tendencies that are calculated once and then held fixed during
!    the RK (Runge-Kutta) sub-steps. 
!  - CM1 uses a low-storage three-step Runge-Kutta scheme.  See Wicker
!    and Skamarock (2002, MWR, p 2088) for more information.
!  - CM1 uses a staggered C grid.  Hence, u arrays have one more grid point
!    in the i direction, v arrays have one more grid point in the j 
!    direction, and w arrays have one more grid point in the k direction
!    (compared to scalar arrays).
!  - CM1 assumes the subgrid turbulence parameters (tke,km,kh) are located
!    at the w points. 
!-----------------------------------------------------------------------------

  implicit none

  private
  public :: solve1

  CONTAINS

      subroutine solve1(nstep,num_soil_layers,                       &
                   dt,dtlast,adtlast,mtime,dbldt,                    &
                   dosfcflx,bud,bud2,qbudget,asq,bsq,                &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,                            &
                   tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,         &
                   wprof,ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,           &
                   uavg,vavg,thavg,pavg,ulspg,vlspg,qavg,cavg,       &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,              &
                   qi0,rr0,rf0,rrf0,thrd,                            &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv, &
                   znt,zntmp,ust,thflux,qvflux,cm0,                  &
                   radbcw,radbce,radbcs,radbcn,                      &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,     &
                   divx,rho,rr,rf,prs,                               &
                   t11,t12,t13,t22,t23,t33,                          &
                   m11,m12,m13,m22,m23,m33,                          &
                   u0,ua,u3d,uten1,                                  &
                   v0,va,v3d,vten1,                                  &
                   rrw,wa,w3d,wten,wten1,                            &
                   ppi,pp3d,ppten,sten,sadv,                         &
                   tha,th3d,thten1,                                  &
                   qpten,qtten,qvten,qcten,qa,q3d,qten,              &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                   nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2,     &
                   thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten,o30,zir,  &
                   swten,lwten,thraten,                              &
                   thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten,    &
                   pta,pt3d,ptten,                                   &
                   gamk,gamwall,kmw,ufw,vfw,u1b,v1b,u2pt,v2pt,ufwk,vfwk, &
                   tdiag,qdiag,udiag,vdiag,wdiag,kdiag,pdiag,        &
                   out2d,out3d,                                      &
                   recy_cap,recy_inj,recywe,recysn,                  &
                   lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,        &
                   timavg,                                           &
                   dowriteout,dorad,getdbz,getvt,dotdwrite,          &
                   dotbud,doqbud,doubud,dovbud,dowbud,donudge,       &
                   doazimwrite,dorestart)
        ! end_solve1
      use input
      use constants
      use bc_module
      use diff2_module
      use turbtend_module
      use misclibs
      use simple_phys_module, only : testcase_simple_phys,get_avg_uvt,get_avg_uvtq
      use eddy_recycle
      use lsnudge_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
! Arrays and variables passed into solve

      integer, intent(in) :: nstep
      integer, intent(in) :: num_soil_layers
      real, intent(inout) :: dt,dtlast
      double precision, intent(in)    :: adtlast
      double precision, intent(in   ) :: mtime
      double precision, intent(inout) :: dbldt
      logical, intent(in) :: dosfcflx
      double precision, intent(inout), dimension(nk) :: bud
      double precision, intent(inout), dimension(nj) :: bud2
      double precision, intent(inout), dimension(nbudget) :: qbudget
      double precision, intent(inout), dimension(numq) :: asq,bsq
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, intent(in),    dimension(kb:ke) :: wprof
      real, intent(inout), dimension(kb:ke) :: ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,  &
                                               uavg,vavg,thavg,pavg,ulspg,vlspg
      real, intent(inout), dimension(kb:ke,numq) :: qavg
      double precision, intent(inout), dimension(kb:ke,3+numq) :: cavg
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: qi0,rr0,rf0,rrf0
      real, intent(in), dimension(ibb2:ibe2,jbb2:jbe2,kbb2:kbe2) :: thrd
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(ib:ie,jb:je) :: znt,zntmp,ust,thflux,qvflux
      real, intent(in),    dimension(ib:ie,jb:je) :: cm0
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rr,rf,prs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(ibnba:ienba,jbnba:jenba,kbnba:kenba) :: m11,m12,m13,m22,m23,m33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua,u3d,uten1
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va,v3d,vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten,wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten,sten,sadv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten1
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,lenscl,dissten,epst,epsd1,epsd2
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten,qnipten,qncpten
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      real, intent(inout), dimension(ibr:ier,jbr:jer) :: zir
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,thraten
      real, intent(in), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: thcuten,qvcuten,qccuten,qicuten,ucuten,vcuten
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d,ptten
      real, intent(in), dimension(kb:ke) :: gamk,gamwall,kmw,ufw,vfw,u1b,v1b
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt,v2pt
      real, intent(in),    dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: ufwk,vfwk
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
      real, intent(inout) , dimension(ibdk:iedk,jbdk:jedk,kbdk:kedk,nkdiag) :: kdiag
      real, intent(inout) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag
      real, intent(inout) , dimension(ib2d:ie2d,jb2d:je2d,nout2d) :: out2d
      real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d

      integer, intent(inout), dimension(ib:ie,jb:je) :: recy_cap,recy_inj
      real, intent(inout), dimension(irecywe,jrecywe,krecy,nrecy) :: recywe
      real, intent(inout), dimension(irecysn,jrecysn,krecy,nrecy) :: recysn

      real, intent(inout), dimension(lsnudge_maxlevels,lsnudge_nmax) :: lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv

      real, intent(in), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg

      logical, intent(in) :: dowriteout,dorad,dotdwrite,doazimwrite,dorestart
      logical, intent(inout) :: getdbz,getvt,dotbud,doqbud,doubud,dovbud,dowbud,donudge

!-----------------------------------------------------------------------
! Arrays and variables defined inside solve

      integer :: i,j,k,n,nrk,bflag,pdef,diffit,k1,k2,ii,ip
      integer :: has_reqc,has_reqi,has_reqs,do_radar_ref,ke_diag

      real :: delqv,delpi,delth,delt,fac,epsd,dheat,dz1,xs
      real :: foo1,foo2
      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad
      real :: r1,r2,tnew,pnew,pinew,thnew,qvnew
      real :: gamm,aiu,fr
      real :: qv,qli,cpli,cpm,cvm,qmax
      real :: tn,qn,nudgefac,taunudge,timenudge

      real :: umod,uref,zref,oldval,newval

      double precision :: weps,afoo,bfoo,p0,p2

      logical :: get_time_avg,tqnudge,dopf
      logical :: do_adapt_lspg

#ifdef MPI
      integer :: reqp1,reqp2
#endif

!--------------------------------------------------------------------

#ifdef MPI
      nf=0
      nu=0
      nv=0
      nw=0
#endif

      afoo=0.0d0
      bfoo=0.0d0


      if( testcase.eq.10 )then
        tqnudge = .true.
      else
        tqnudge = .false.
      endif


      timenudge = mtime+dt
      donudge = .false.
      lsnudgefac = 0.0

      if( do_lsnudge )then
        if( (mtime+dt).gt.lsnudge_start .and. (mtime+dt).lt.lsnudge_end )then

          donudge = .true.
          lsnudgefac = 1.0

          if( donudge .and. timenudge.gt.lsnudge_time2(lsnudge_count) )then

            ! update which profile we are nudging towards:
            lsnudge_count = lsnudge_count + 1

            if( myid.eq.0 ) print *,'  lsnudge times:  time1,timenudge,time2 = ',lsnudge_time1(lsnudge_count),timenudge,lsnudge_time2(lsnudge_count)
          endif

          ! ramp-up nudging starting at lsnudge_start:
          if( lsnudge_ramp_time.lt.smeps )then
            lsnudgefac = 1.0
          else
            lsnudgefac = (timenudge-lsnudge_time1(lsnudge_count))/lsnudge_ramp_time
            lsnudgefac = max( lsnudgefac , 0.0 )
            lsnudgefac = min( lsnudgefac , 1.0 )
          endif
          if( lsnudgefac.gt.0.001 .and. lsnudgefac.lt.0.999 .and. myid.eq.0 ) print *,'  lsnudgefac = ',lsnudgefac

        endif
      endif

      if(timestats.ge.1) time_misc=time_misc+mytime()

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   subgrid turbulence schemes  cccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( idoles .and. sgsmodel.eq.3 )THEN

        ! wten = kmh
        ! rrw  = kmv

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=0,nj+1
        do i=0,ni+1
          wten(i,j,k) = gamk(k)*kmh(i,j,k)
           rrw(i,j,k) = gamk(k)*kmv(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------------
!  get RHS for tke scheme:

      IF( idoles .and. iusetke )THEN

        !  Buoyancy, Dissipation, and mean-gradient (ie, shear) terms:
        !    cm1r21.1: switch back to shear-prod formulation from cm1r19 and earlier
        IF(tconfig.eq.1)THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            tketen(i,j,k) = -khv(i,j,k)*nm(i,j,k)  &
                            -dissten(i,j,k)  &
                            +kmv(i,j,k)*max(0.0,(defv(i,j,k)+defh(i,j,k)))
          enddo
          enddo
          enddo
        ELSEIF(tconfig.eq.2)THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=2,nk
          do j=1,nj
          do i=1,ni
            tketen(i,j,k) = -khv(i,j,k)*nm(i,j,k)  &
                            -dissten(i,j,k)  &
                            +kmv(i,j,k)*max(0.0,defv(i,j,k))   &
                            +kmh(i,j,k)*max(0.0,defh(i,j,k))
          enddo
          enddo
          enddo
        ENDIF
        if(timestats.ge.1) time_turb=time_turb+mytime()

        if( dotdwrite .and. kd_turb.ge.1 )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            kdiag(i,j,k,kd_turb) = tketen(i,j,k)
          enddo
          enddo
          enddo
        endif

      if( sgsmodel.eq.3 )then
        call turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,  &
                   rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,        &
                   dum1,dum2,dum3,dum4,dum5,sten,tkea,tketen,wten,rrw) ! kmh = wten, kmv = rrw
      else
        call turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,  &
                   rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,        &
                   dum1,dum2,dum3,dum4,dum5,sten,tkea,tketen,kmh,kmv)
      endif

        if( dotdwrite .and. kd_turb.ge.1 )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            kdiag(i,j,k,kd_turb) = tketen(i,j,k)-kdiag(i,j,k,kd_turb)
          enddo
          enddo
          enddo
        endif

      ENDIF


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Pre-RK calculations   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    IF( irdamp.eq.2 .or. hrdamp.eq.2 .or. tqnudge .or. do_lsnudge )THEN

      if( imoist.eq.1 )then
        call     get_avg_uvtq(uavg,vavg,thavg,qavg,cavg,th0,ua,va,tha,qa,ruh,ruf,rvh,rvf)
      else
        call     get_avg_uvt(uavg,vavg,thavg,cavg,th0,ua,va,tha,ruh,ruf,rvh,rvf)
      endif

    ENDIF

!--------------------------------------------------------------------
!  radbc
 
      if(irbc.eq.1)then

        if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
        if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)

      endif

!--------------------------------------------------------------------
!  Rayleigh damping and set tendency arrays:

      IF( irdamp.eq.1 .or. hrdamp.eq.1 )THEN
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          uten1(i,j,k) = -rdalpha*0.5*(tauh(i-1,j,k)+tauh(i,j,k))*(ua(i,j,k)-u0(i,j,k))
          vten1(i,j,k) = -rdalpha*0.5*(tauh(i,j-1,k)+tauh(i,j,k))*(va(i,j,k)-v0(i,j,k))
          wten1(i,j,k) = -rdalpha*tauf(i,j,k)*wa(i,j,k)
          thten1(i,j,k) = -rdalpha*taus(i,j,k)*tha(i,j,k)
        enddo
        enddo
        enddo
      ELSEIF( irdamp.eq.2 .or. hrdamp.eq.2 )THEN
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          uten1(i,j,k) = -rdalpha*0.5*(tauh(i-1,j,k)+tauh(i,j,k))*(ua(i,j,k)-uavg(k))
          vten1(i,j,k) = -rdalpha*0.5*(tauh(i,j-1,k)+tauh(i,j,k))*(va(i,j,k)-vavg(k))
          wten1(i,j,k) = -rdalpha*tauf(i,j,k)*wa(i,j,k)
          thten1(i,j,k) = -rdalpha*taus(i,j,k)*((th0(i,j,k)-thavg(k))+tha(i,j,k))
        enddo
        enddo
        enddo
      ELSE
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          uten1(i,j,k) = 0.0
          vten1(i,j,k) = 0.0
          wten1(i,j,k) = 0.0
          thten1(i,j,k) = 0.0
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      IF( imoist.eq.1 )THEN
        ! moisture:
        do n=1,numq
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
      ENDIF

      IF(iptra.eq.1)THEN
        ! passive tracers:
        do n=1,npt
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ptten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
        enddo
      ENDIF

    ! Diagnostics (infrequent):

      if( doubud .and. ud_rdamp.ge.1 )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          udiag(i,j,k,ud_rdamp) = uten1(i,j,k)
        enddo
        enddo
        enddo
      endif
      if( dovbud .and. vd_rdamp.ge.1 )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vdiag(i,j,k,vd_rdamp) = vten1(i,j,k)
        enddo
        enddo
        enddo
      endif
      if( dowbud .and. wd_rdamp.ge.1 )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,wd_rdamp) = wten1(i,j,k)
        enddo
        enddo
        enddo
      endif
      if( dotbud .and. td_rdamp.ge.1 )then
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          tdiag(i,j,k,td_rdamp) = thten1(i,j,k)
        enddo
        enddo
        enddo
      endif
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

!--------------------------------------------------------------------
!  large-scale pressure gradient:

      IF( lspgrad.ge.1 )THEN

        ! Include a large-scale pressure gradient:

          do_adapt_lspg = .false.

          IF( do_adapt_lspg .and. lspgrad.eq.4 )THEN
!!!            stop 3333

          if( mtime.gt.  10.0 .and. mod(mtime,  1.0).lt.1.0e-8 )then

            call get_avg_uvt(uavg,vavg,thavg,cavg,th0,ua,va,tha,ruh,ruf,rvh,rvf)

            ! 200529: modify lspg to match reference wind
            if( myid.eq.0 )then

!!!              zref = 240.0
!!!              uref = 16.7303

              if( testcase.eq.12 )then
                ! Martinuzzi and Tropea (1993, JFE) wind tunnel test case:
                zref =  0.15
                uref =  3.0
              endif

              k = 1
              do while( zh(1,1,k).lt.zref )
                k = k+1
              enddo

              k1 = k-1
              k2 = k

              umod = uavg(k1)+(uavg(k2)-uavg(k1))*(zref-zh(1,1,k2))/(zh(1,1,k1)-zh(1,1,k2))

              oldval = ulspg(1)
              newval = oldval*uref/umod

              print *,'  Get lspg:  ',mtime,mod(mtime,1.0)
              print *,'     old u,pg',umod,oldval
              print *,'     new u,pg',uref,newval
            endif

#ifdef MPI
            call MPI_BCAST( newval  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
#endif

            ulspg = newval
   
          endif
          ENDIF

        if(     lspgrad.eq.1 )then

          ulspg = 0.0
          vlspg = 0.0

          !---------------------------------------------------------------!
          ! Large-scale pressure gradient based on geostropic balance,
          ! using base-state wind profiles:
          !---------------------------------------------------------------!
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=1,ni+1
            ! 170728 bug fix:  account for grid staggering:
            ! 180618:  when imove=1, need to add vmove to base state
            uten1(i,j,k) = uten1(i,j,k)-fcor*( 0.25*( (v0(i  ,j,k)+v0(i  ,j+1,k))   &
                                                     +(v0(i-1,j,k)+v0(i-1,j+1,k)) ) + vmove )
            vten1(i,j,k) = vten1(i,j,k)+fcor*( 0.25*( (u0(i,j  ,k)+u0(i+1,j  ,k))   &
                                                     +(u0(i,j-1,k)+u0(i+1,j-1,k)) ) + umove )
          enddo
          enddo
          enddo

        elseif( lspgrad.eq.2 )then

          !---------------------------------------------------------------!
          ! Large-scale pressure gradient based on geostropic balance,
          ! using ug,vg arrays:
          !---------------------------------------------------------------!
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
            ulspg(k) = -fcor*(vg(k)+vmove)
            vlspg(k) =  fcor*(ug(k)+umove)
            do j=1,nj+1
            do i=1,ni+1
              uten1(i,j,k) = uten1(i,j,k)+ulspg(k)
              vten1(i,j,k) = vten1(i,j,k)+vlspg(k)
            enddo
            enddo
          enddo

        elseif( lspgrad.eq.3 )then

          !---------------------------------------------------------------!
          ! Large-scale pressure gradient based on gradient-wind balance,
          ! (Bryan et al, 2017, BLM, eqn 10)
          ! using base-state wind profiles:
          !---------------------------------------------------------------!

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,fr)
          do k=1,nk

            ! B17 baseline:
            fr = -fcor*vg(k)-vg(k)*vg(k)/hurr_rad

            ! rotate into Cartesian coordinates:
            ulspg(k) = fr*cos(-hurr_angle)
            vlspg(k) = fr*sin(-hurr_angle)

            do j=1,nj+1
            do i=1,ni+1
              uten1(i,j,k) = uten1(i,j,k)+ulspg(k)
              vten1(i,j,k) = vten1(i,j,k)+vlspg(k)
            enddo
            enddo

          enddo

        elseif( lspgrad.eq.4 )then

          !---------------------------------------------------------------!
          ! Specified large-scale pressure gradient:
          ! (user must set ulspg and vlspg in base.F)
          !---------------------------------------------------------------!
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
            do j=1,nj+1
            do i=1,ni+1
              uten1(i,j,k) = uten1(i,j,k)+ulspg(k)
              vten1(i,j,k) = vten1(i,j,k)+vlspg(k)
            enddo
            enddo
          enddo

        endif
        if(timestats.ge.1) time_misc=time_misc+mytime()

      ENDIF

!--------------------------------------------------------------------
!  2nd-order diffusion, fixed viscosity (eg, DNS and idealized test):

      IF( difforder.eq.2 )THEN

        idiffge1:  &
        if(idiff.ge.1)then

          call diff2u(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,uten1,ust,rho,rr,rf,divx,t11,t12,t13, &
                      doubud,udiag)
          call diff2v(xh,arh1,arh2,uh,rxf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,vten1,ust,rho,rr,rf,divx,t22,t12,t23, &
                      dovbud,vdiag)
          call diff2w(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,wten1,rho,rr,rf,divx,t33,t13,t23,  &
                      dowbud,wdiag)

        endif  idiffge1

        idiffeq1:  &
        if(idiff.eq.1)then

          call diff2s(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                      dum1,dum2,dum3,dum4,tha,thten1,rho,rr,rf,  &
                      dotbud,ibdt,iedt,jbdt,jedt,kbdt,kedt,ntdiag,tdiag,td_hediff,td_vediff)

          IF( imoist.eq.1 )THEN
            ! moisture:
            do n=1,numq
              if( n.eq.nqv )then
                call diff2s(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                            dum1,dum2,dum3,dum4,qa(ib,jb,kb,n),qten(ib,jb,kb,n),rho,rr,rf,  &
                            doqbud,ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,qdiag,qd_hediff,qd_vediff)
              else
                call diff2s(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                            dum1,dum2,dum3,dum4,qa(ib,jb,kb,n),qten(ib,jb,kb,n),rho,rr,rf,  &
                            .false.,ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,qdiag,1,1)
              endif
            enddo
          ENDIF

          IF(iptra.eq.1)THEN
            ! passive tracers:
            do n=1,npt
              call diff2s(rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,mh,mf,  &
                          dum1,dum2,dum3,dum4,pta(ib,jb,kb,n),ptten(ib,jb,kb,n),rho,rr,rf,  &
                          .false.,ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,qdiag,1,1)
            enddo
          ENDIF

        endif  idiffeq1

      ENDIF

!--------------------------------------------------------------------
!  tendencies specific to THETA:


      !----- cvm (if needed) -----!

      IF( eqtset.eq.2 .and. imoist.eq.1 .and. (idiss.eq.1.or.rterm.eq.1) )THEN
        ! for energy-conserving moist thermodynamics:
        ! store cvm in dum1:
        ! store ql  in dum2:
        ! store qi  in dum3:

        IF( nql1.ge.1 )THEN

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum2(i,j,k)=qa(i,j,k,nql1)
          enddo
          enddo
          enddo

          do n=nql1+1,nql2
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=dum2(i,j,k)+qa(i,j,k,n)
            enddo
            enddo
            enddo
          enddo

        ELSE

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum2(i,j,k)=0.0
          enddo
          enddo
          enddo

        ENDIF

        IF(iice.eq.1)THEN

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum3(i,j,k)=qa(i,j,k,nqs1)
          enddo
          enddo
          enddo

          do n=nqs1+1,nqs2
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=dum3(i,j,k)+qa(i,j,k,n)
            enddo
            enddo
            enddo
          enddo

        ELSE

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum3(i,j,k)=0.0
          enddo
          enddo
          enddo

        ENDIF

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=cv+cvv*qa(i,j,k,nqv)+cpl*dum2(i,j,k)+cpi*dum3(i,j,k)
        enddo
        enddo
        ENDDO

      ELSE

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,n)
        DO k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=cv
        enddo
        enddo
        ENDDO

      ENDIF

      !----- store appropriate rho for budget calculations in dum2 -----!

      IF(axisymm.eq.1)THEN
       ! for axisymmetric grid:
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
        enddo
        enddo
        enddo
      ELSE
       ! for Cartesian grid:
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

      !-------------------------------------------------------------

      !  budget calculations:
      if(dosfcflx.and.imoist.eq.1)then
        tem0 = dt*dx*dy*dz
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,delpi,delth,delqv,delt,n)
        do j=1,nj
        bud2(j) = 0.0d0
        do i=1,ni
          k = 1
          delth = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*thflux(i,j)
          delqv = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*qvflux(i,j)
          delpi = rddcv*(pi0(i,j,1)+ppi(i,j,1))*(           &
                                delqv/(eps+qa(i,j,1,nqv))   &
                               +delth/(th0(i,j,1)+tha(i,j,1))  )
          delt = (pi0(i,j,k)+ppi(i,j,k))*delth   &
                +(th0(i,j,k)+tha(i,j,k))*delpi
          bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*(        &
                  cv*delt                                                   &
                + cvv*qa(i,j,k,nqv)*delt                                    &
                + cvv*(pi0(i,j,k)+ppi(i,j,k))*(th0(i,j,k)+tha(i,j,k))*delqv &
                + g*zh(i,j,k)*delqv   )
        if( nql1.ge.1 )then
          do n=nql1,nql2
            bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpl*qa(i,j,k,n)*delt
          enddo
        endif
          if(iice.eq.1)then
            do n=nqs1,nqs2
              bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpi*qa(i,j,k,n)*delt
            enddo
          endif
        enddo
        enddo
        do j=1,nj
          qbudget(9) = qbudget(9) + tem0*bud2(j)
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
      endif

      !---- Dissipative heating term:

      IF( idiss.eq.1 )THEN
        IF( dotbud .and. td_diss.ge.1 )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_diss) = thten1(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
        ! note:  dissten array stores epsilon (dissipation rate) at w points
        if( bbc.eq.3 )then
          k1 = 2
        else
          k1 = 1
        endif
        if(imoist.eq.1.and.eqtset.eq.2)then
          ! moist, new equations:
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,epsd,dheat)
          do k=k1,nk
          do j=1,nj
          do i=1,ni
            epsd = 0.5*(dissten(i,j,k)+dissten(i,j,k+1))
            dheat=epsd/( cpdcv*dum1(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) )
            thten1(i,j,k)=thten1(i,j,k)+dheat
          enddo
          enddo
          enddo
          if( bbc.eq.3 )then
            k = 1
            !$omp parallel do default(shared)  &
            !$omp private(i,j,dz1,epsd,dheat)
            do j=1,nj
            do i=1,ni
              dz1 = zf(i,j,2)-zf(i,j,1)
              epsd = (ust(i,j)**3)*alog((dz1+zntmp(i,j))/zntmp(i,j))/(karman*dz1)
              dheat=epsd/( cpdcv*dum1(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) )
              thten1(i,j,k)=thten1(i,j,k)+dheat
            enddo
            enddo
          endif
        else
          ! traditional cloud-modeling equations (also dry equations):
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,epsd,dheat)
          do k=k1,nk
          do j=1,nj
          do i=1,ni
            epsd = 0.5*(dissten(i,j,k)+dissten(i,j,k+1))
            dheat=epsd/( cp*(pi0(i,j,k)+ppi(i,j,k)) )
            thten1(i,j,k)=thten1(i,j,k)+dheat
          enddo
          enddo
          enddo
          if( bbc.eq.3 )then
            k = 1
            !$omp parallel do default(shared)  &
            !$omp private(i,j,dz1,epsd,dheat)
            do j=1,nj
            do i=1,ni
              dz1 = zf(i,j,2)-zf(i,j,1)
              epsd = (ust(i,j)**3)*alog((dz1+zntmp(i,j))/zntmp(i,j))/(karman*dz1)
              dheat=epsd/( cp*(pi0(i,j,k)+ppi(i,j,k)) )
              thten1(i,j,k)=thten1(i,j,k)+dheat
            enddo
            enddo
          endif
        endif
        IF( dotbud .and. td_diss.ge.1 )THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_diss) = thten1(i,j,k)-tdiag(i,j,k,td_diss)
          enddo
          enddo
          enddo
        ENDIF
      ENDIF

      !---- Rotunno-Emanuel "radiation" term
      !---- (currently capped at 2 K/day ... see RE87 p 546)

      IF(rterm.eq.1)THEN
        if( dotbud .and. td_rad.ge.1 )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_rad) = thten1(i,j,k)
          enddo
          enddo
          enddo
        endif
        tem0 = dt*dx*dy*dz
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,thrad,prad)
        do k=1,nk
        bud(k)=0.0d0
        do j=1,nj
        do i=1,ni
          ! NOTE:  thrad is a POTENTIAL TEMPERATURE tendency
          thrad = -tha(i,j,k)/(12.0*3600.0)
          if( tha(i,j,k).gt. 1.0 ) thrad = -1.0/(12.0*3600.0)
          if( tha(i,j,k).lt.-1.0 ) thrad =  1.0/(12.0*3600.0)
          thten1(i,j,k)=thten1(i,j,k)+thrad
          ! associated pressure tendency:
          prad = (pi0(i,j,k)+ppi(i,j,k))*rddcv*thrad/(th0(i,j,k)+tha(i,j,k))
          ! budget:
          bud(k) = bud(k) + dum1(i,j,k)*dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*( &
                            thrad*(pi0(i,j,k)+ppi(i,j,k))    &
                           + prad*(th0(i,j,k)+tha(i,j,k)) )
        enddo
        enddo
        enddo
        do k=1,nk
          qbudget(10) = qbudget(10) + tem0*bud(k)
        enddo
        if( dotbud .and. td_rad.ge.1 )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_rad) = thten1(i,j,k)-tdiag(i,j,k,td_rad)
          enddo
          enddo
          enddo
        endif
      ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()

      IF( radopt.ge.1 )THEN
        ! Radiation tendencies:

        ! Notes:
        ! TEMPERATURE tendencies from radiation scheme are stored in lwten and swten
        ! since cm1r21.1: potential temp tendencies are stored in thraten

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k) = thten1(i,j,k) + thraten(i,j,k)
        enddo
        enddo
        enddo

        if( dotbud .and. td_rad.ge.1 )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_rad) = thraten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(timestats.ge.1) time_rad=time_rad+mytime()
      ENDIF

      IF( use_pbl )THEN
        ! PBL tendencies:

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k) = thten1(i,j,k) + thpten(i,j,k)
        enddo
        enddo
        enddo
        if( dotbud .and. td_pbl.ge.1 )then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            tdiag(i,j,k,td_pbl) = thpten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(timestats.ge.1) time_pbl=time_pbl+mytime()

      ENDIF


!--------------------------------------------------------------------
!  subgrid turbulence tendencies:

      doturb:  &
      IF( dohturb .or. dovturb )THEN

        if( sgsmodel.eq.5 .or. sgsmodel.eq.6 )then
          ! NBA sgs turbulence models:

          call turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                     zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                     dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,ua,uten1,wa,m11,m12,m13,m22,kmv,cm0, &
                     kmw,ufw,u1b,u2pt,ufwk,doubud,udiag)
          call turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                     zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                     dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,va,vten1,wa,m12,m22,m23,kmv,cm0, &
                     kmw,vfw,v1b,v2pt,vfwk,dovbud,vdiag)
          call turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                     dum1,dum2,dum3,dum4,dum5,dum6,wa,wten1,m13,m23,m33,m22,kmh,cm0,  &
                     dowbud,wdiag)

        elseif( sgsmodel.eq.3 )then
          ! SMM94 two-part model:
          ! kmh = wten, kmv = rrw

          call turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                     zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                     dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,ua,uten1,wa,t11,t12,t13,t22,rrw,cm0, &
                     kmw,ufw,u1b,u2pt,ufwk,doubud,udiag)
          call turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                     zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                     dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,va,vten1,wa,t12,t22,t23,rrw,cm0, &
                     kmw,vfw,v1b,v2pt,vfwk,dovbud,vdiag)
          call turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                     dum1,dum2,dum3,dum4,dum5,dum6,wa,wten1,t13,t23,t33,t22,wten,cm0,  &
                     dowbud,wdiag)

        else
          ! all other sgs turbulence models:

          call turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                     zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                     dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,ua,uten1,wa,t11,t12,t13,t22,kmv,cm0, &
                     kmw,ufw,u1b,u2pt,ufwk,doubud,udiag)
          call turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                     zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                     dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,va,vten1,wa,t12,t22,t23,kmv,cm0, &
                     kmw,vfw,v1b,v2pt,vfwk,dovbud,vdiag)
          call turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                     dum1,dum2,dum3,dum4,dum5,dum6,wa,wten1,t13,t23,t33,t22,kmh,cm0,  &
                     dowbud,wdiag)

        endif

        if( ( doimpl.eq.1 .and. dovturb ) .or. ( cm1setup.eq.2 .and. ipbl.eq.2 ) )then
          !  Arrays for vimpl turbs:
          !    NOTE:  do not change dum7,dum8 until all calls to turbs are finished
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum7(i,j,k) = khv(i,j,k  )*mf(i,j,k  )*rf(i,j,k  )*mh(i,j,k)*rr(i,j,k)
            dum8(i,j,k) = khv(i,j,k+1)*mf(i,j,k+1)*rf(i,j,k+1)*mh(i,j,k)*rr(i,j,k)
          enddo
          enddo
          enddo
        endif

        ! cm1r18: subtract th0r from theta (as in advection scheme)
        !         (reduces roundoff error)
        IF(.not.terrain_flag)THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,tem)
          do k=1,nk
          tem = th0(1,1,k)-th0r
          do j=0,nj+1
          do i=0,ni+1
            sadv(i,j,k)=tem+tha(i,j,k)
          enddo
          enddo
          enddo
        ELSE
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,tem)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            sadv(i,j,k)=(th0(i,j,k)-th0r)+tha(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
        call turbs(1,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,thflux,   &
                   rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                   dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,sadv,thten1,khh,khv,cm0,dum7,dum8, &
                   dotbud,ibdt,iedt,jbdt,jedt,kbdt,kedt,ntdiag,tdiag,td_hturb,td_vturb,1)

        IF( imoist.eq.1 )THEN
          ! moisture:
          do n=1,numq
            if( n.eq.nqv )then
              call turbs(1,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,qvflux,   &
                         rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                         dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,cm0,dum7,dum8, &
                         doqbud,ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,qdiag,qd_hturb,qd_vturb,2)
            else
              call turbs(0,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,qvflux,   &
                         rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                         dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,cm0,dum7,dum8, &
                         .false.,ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,qdiag,1,1,0)
            endif
          enddo
        ENDIF

        IF(iptra.eq.1)THEN
          ! passive tracers:
          do n=1,npt
            call turbs(0,dt,dosfcflx,xh,rxh,arh1,arh2,uh,xf,arf1,arf2,uf,vh,vf,qvflux,   &
                       rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                       dum1,dum2,dum3,dum4,dum5,sten,rho,rr,rf,pta(ib,jb,kb,n),ptten(ib,jb,kb,n),khh,khv,cm0,dum7,dum8, &
                       .false.,ibdq,iedq,jbdq,jedq,kbdq,kedq,nqdiag,qdiag,1,1,0)
          enddo
        ENDIF

      ENDIF  doturb

!--------------------------------------------------------------------
!  U-equation

      IF( iinit.eq.10 .and. mtime.lt.t2_uforce )THEN
        ! u-forcing for squall-line initialization:
        ! (Morrison et al, 2015, JAS, pg 315)
        gamm = 1.0
        if(mtime.ge.t1_uforce)THEN
          gamm = 1.0+(0.0-1.0)*(mtime-t1_uforce)/(t2_uforce-t1_uforce)
        endif
        if(myid.eq.0) print *,'  mtime,gamm = ',mtime,gamm
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,aiu)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          if( abs(xf(i)-xc_uforce).lt.xr_uforce .and. abs(zf(i,j,k)-zs(i,j)).lt.zr_uforce )then
            aiu = alpha_uforce*cos(0.5*pi*(xf(i)-xc_uforce)/xr_uforce)   &
                              *((cosh(2.5*(zf(i,j,k)-zs(i,j))/zr_uforce))**(-2))
            uten1(i,j,k)=uten1(i,j,k)+gamm*aiu
          endif
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if( use_pbl )then
        if( axisymm.eq.1 )then
          ip = 2
        else
          ip = 1
        endif
        !$omp parallel do default(shared)   &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=ip,ni+1
          uten1(i,j,k) = uten1(i,j,k) + 0.5*( upten(i-1,j,k)+ upten(i,j,k))
        enddo
        enddo
        enddo
        IF( doubud .and. ud_pbl.ge.1 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=ip,ni+1
            udiag(i,j,k,ud_pbl) = 0.5*( upten(i-1,j,k)+ upten(i,j,k))
          enddo
          enddo
          enddo
        ENDIF
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif

!--------------------------------------------------------------------
!  V-equation

      if( use_pbl )then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten1(i,j,k) = vten1(i,j,k) + 0.5*( vpten(i,j-1,k)+ vpten(i,j,k))
        enddo
        enddo
        enddo
        IF( dovbud .and. vd_pbl.ge.1 )then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj+1
          do i=1,ni
            vdiag(i,j,k,vd_pbl) = 0.5*( vpten(i,j-1,k)+ vpten(i,j,k))
          enddo
          enddo
          enddo
        ENDIF
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif
 
!-------------------------------------------------------------------
!  Moisture

      if(imoist.eq.1)then
            if( use_pbl )then
              !$omp parallel do default(shared)   &
              !$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                qten(i,j,k,nqv) = qten(i,j,k,nqv) + qvpten(i,j,k)
              enddo
              enddo
              enddo
              if( doqbud .and. qd_pbl.ge.1 )then
                !$omp parallel do default(shared)   &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_pbl) = qvpten(i,j,k)
                enddo
                enddo
                enddo
              endif
              if(timestats.ge.1) time_pbl=time_pbl+mytime()
            endif
        IF( use_pbl )THEN
        if(nqc.ge.1)then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqc) = qten(i,j,k,nqc) + qcpten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(nqi.ge.1)then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqi) = qten(i,j,k,nqi) + qipten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(nnci.ge.1)then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nnci) = qten(i,j,k,nnci) + qnipten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(nncc.ge.1)then
          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nncc) = qten(i,j,k,nncc) + qncpten(i,j,k)
          enddo
          enddo
          enddo
        endif
          if(timestats.ge.1) time_pbl=time_pbl+mytime()
        ENDIF
      endif

!--------------------------------------------------------------------

    IF( tqnudge )THEN

      if( imoist.eq.0 .and. td_nudge.ge.1   &
          .and. dotbud )then
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          tdiag(i,j,k,td_nudge) = thten1(i,j,k)
        enddo
        enddo
        enddo
      elseif( imoist.eq.1 .and. td_nudge.ge.1   &
          .and. dotbud .and. doqbud .and. qd_nudge.ge.1 )then
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          tdiag(i,j,k,td_nudge) = thten1(i,j,k)
          qdiag(i,j,k,qd_nudge) = qten(i,j,k,nqv)
        enddo
        enddo
        enddo
      endif

      ! ramp-up nudging at start of simulation:
      nudgefac = (mtime-0.0)/(900.0-0.0)
      nudgefac = max(0.0,nudgefac)
      nudgefac = min(1.0,nudgefac)

      ! nudging time scale (seconds):
      taunudge = 300.0

      if( abs(nudgefac).gt.0.001 )then
        ! make sure subsidence is off:
        dolsw = .false.
        if( myid.eq.0 ) print *,'  nudging with nudgefac,taunudge = ',nudgefac,taunudge
        tem = nudgefac/taunudge
        if( imoist.eq.0 )then
          do k=1,nk
            tn = -(thavg(k)-th0(1,1,k))*tem
            do j=1,nj+1
            do i=1,ni+1
              thten1(i,j,k) = thten1(i,j,k) + tn
            enddo
            enddo
          enddo
        else
          do k=1,nk
            tn = -(thavg(k)-th0(1,1,k))*tem
            qn = -(qavg(k,nqv)-qv0(1,1,k))*tem
            do j=1,nj+1
            do i=1,ni+1
              thten1(i,j,k) = thten1(i,j,k) + tn
              qten(i,j,k,nqv) = qten(i,j,k,nqv) + qn
            enddo
            enddo
          enddo
        endif
      else
        if(myid.eq.0) print *,'  no nudging (yet) '
      endif

      if( imoist.eq.0 .and. td_nudge.ge.1   &
          .and. dotbud )then
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          tdiag(i,j,k,td_nudge) = thten1(i,j,k)-tdiag(i,j,k,td_nudge)
        enddo
        enddo
        enddo
      elseif( imoist.eq.1 .and. td_nudge.ge.1   &
          .and. dotbud .and. doqbud .and. qd_nudge.ge.1 )then
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          tdiag(i,j,k,td_nudge) = thten1(i,j,k)-tdiag(i,j,k,td_nudge)
          qdiag(i,j,k,qd_nudge) = qten(i,j,k,nqv)-qdiag(i,j,k,qd_nudge)
        enddo
        enddo
        enddo
      endif

      if( terrain_flag ) stop 82701

      if(timestats.ge.1) time_misc=time_misc+mytime()

    ENDIF

!--------------------------------------------------------------------
!  TENDENCIES for pre-configured cases:
!    (new for cm1r19)

      IF( testcase.ge.1 )THEN

        call     testcase_simple_phys(mh,rho0,rr0,rf0,th0,u0,v0,     &
                   zh,zf,dum1,dum2,dum3,dum4,dum5,dum6,              &
                   ufrc,vfrc,thfrc,qvfrc,ug,vg,dvdr,                 &
                   uavg,vavg,thavg,qavg,cavg,                        &
                   ua,va,tha,qa,uten1,vten1,thten1,qten,             &
                   o30 ,thraten,zir,ruh,ruf,rvh,rvf,mtime)
        if(timestats.ge.1) time_misc=time_misc+mytime()

      ENDIF

!-----------------------------------------------------------------------------
!  Convective parameterization:

    IF( cuparam.ge.1 )THEN

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni+1
        thten1(i,j,k) = thten1(i,j,k) + thcuten(i,j,k)
          uten1(i,j,k) = uten1(i,j,k) + 0.5*(ucuten(i-1,j,k)+ucuten(i,j,k))
          vten1(i,j,k) = vten1(i,j,k) + 0.5*(vcuten(i,j-1,k)+vcuten(i,j,k))
      enddo
      enddo
      enddo

      IF(imoist.eq.1)THEN
        if(nqv.ge.1)then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = qten(i,j,k,nqv) + qvcuten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(nqc.ge.1)then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqc) = qten(i,j,k,nqc) + qccuten(i,j,k)
          enddo
          enddo
          enddo
        endif
        if(nqi.ge.1)then
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqi) = qten(i,j,k,nqi) + qicuten(i,j,k)
          enddo
          enddo
          enddo
        endif
      ENDIF

    ENDIF

!-------------------------------------------------------------------
!   Eddy recycler

      doeddyrec:  &
      IF( do_recycle )THEN

!!!        out3d = 0.0

        if( dowriteout )then
          recy_cap = 0.0
          recy_inj = 0.0
        endif

        if( do_recycle_w )then
          call do_eddy_recyw(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recywe,out3d,recy_cap,recy_inj,timavg,adtlast)
        endif
        if( do_recycle_e )then
          call do_eddy_recye(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recywe,out3d,recy_cap,recy_inj,timavg,adtlast)
        endif
        if( do_recycle_s )then
          call do_eddy_recys(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recysn,out3d,recy_cap,recy_inj,timavg,adtlast)
        endif
        if( do_recycle_n )then
          call do_eddy_recyn(dt,xh,xf,yh,yf,zh,zf,u0,v0,dum1,dum2,cm0,  &
                             u3d,v3d,w3d,th3d,q3d,tke3d,                &
                             uten1,vten1,wten1,thten1,qten,tketen,      &
                             recysn,out3d,recy_cap,recy_inj,timavg,adtlast)
        endif

        if(timestats.ge.1) time_ercyl=time_ercyl+mytime()

      ENDIF  doeddyrec

 
!--------------------------------------------------------------------


    IF( do_lsnudge .and. timenudge.gt.lsnudge_start )THEN

        call     get_avg_uvtq(uavg,vavg,thavg,qavg,cavg,th0,ua,va,tha,qa,ruh,ruf,rvh,rvf)

      if( donudge .and. do_lsnudge_u )then
        if( myid.eq.0 ) print *,'  applying large-scale nudging of u, lsnudge_count = ',lsnudge_count
        rdt = 1.0/dt
        do k=1,nk
          tem1 = -lsnudgefac*( (umove+uavg(k))-lsnudge_u(k,lsnudge_count) )/(lsnudge_tau)
          ! apply nudging tendency:
          do j=1,nj
          do i=1,ni+1
            uten1(i,j,k) = uten1(i,j,k) + tem1
          enddo
          enddo
          if( ud_nudge.ge.1 .and. dotdwrite )then
            do j=1,nj
            do i=1,ni+1
              udiag(i,j,k,ud_nudge) = tem1
            enddo
            enddo
          endif
        enddo
      endif
      if( donudge .and. do_lsnudge_v )then
        if( myid.eq.0 ) print *,'  applying large-scale nudging of v, lsnudge_count = ',lsnudge_count
        rdt = 1.0/dt
        do k=1,nk
          tem2 = -lsnudgefac*( (vmove+vavg(k))-lsnudge_v(k,lsnudge_count) )/(lsnudge_tau)
          ! apply nudging tendency:
          do j=1,nj+1
          do i=1,ni
            vten1(i,j,k) = vten1(i,j,k) + tem2
          enddo
          enddo
          if( vd_nudge.ge.1 .and. dotdwrite )then
            do j=1,nj+1
            do i=1,ni
              vdiag(i,j,k,vd_nudge) = tem2
            enddo
            enddo
          endif
        enddo
      endif
      if( donudge .and. do_lsnudge_th )then
        if( myid.eq.0 ) print *,'  applying large-scale nudging of theta, lsnudge_count = ',lsnudge_count
        rdt = 1.0/dt
        do k=1,nk
          tem1 = -lsnudgefac*( thavg(k)-lsnudge_th(k,lsnudge_count) )/(lsnudge_tau)
          ! apply nudging tendency:
          do j=1,nj
          do i=1,ni
            thten1(i,j,k) = thten1(i,j,k) + tem1
          enddo
          enddo
          if( td_nudge.ge.1 .and. dotdwrite )then
            do j=1,nj
            do i=1,ni
              tdiag(i,j,k,td_nudge) = tem1
            enddo
            enddo
          endif
        enddo
      endif
      if( donudge .and. do_lsnudge_qv .and. imoist.eq.1 )then
        if( myid.eq.0 ) print *,'  applying large-scale nudging of qv, lsnudge_count = ',lsnudge_count
        rdt = 1.0/dt
        do k=1,nk
          tem1 = -lsnudgefac*( qavg(k,nqv)-lsnudge_qv(k,lsnudge_count) )/(lsnudge_tau)
          ! apply nudging tendency:
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = qten(i,j,k,nqv) + tem1
          enddo
          enddo
          if( qd_nudge.ge.1 .and. dotdwrite )then
            do j=1,nj
            do i=1,ni
              qdiag(i,j,k,qd_nudge) = tem1
            enddo
            enddo
          endif
        enddo
      endif
      if(timestats.ge.1) time_misc=time_misc+mytime()
    ENDIF


!-------------------------------------------------------------------
!    NOTE:  now ok to change dum7,dum8
!-------------------------------------------------------------------
!  contribution to pressure tendency from potential temperature:
!  (for mass conservation)
!  plus, some other stuff:

    ps6:  &
    IF( psolver.le.3 )THEN

      IF(eqtset.eq.1)THEN
        ! traditional cloud modeling:
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppten(i,j,k)=0.0
        enddo
        enddo
        enddo
      ELSE
        ! mass-conserving pressure eqt:  different sections for moist/dry cases:
        rdt = 1.0/dt
        tem = 0.0001*tsmall
        IF(imoist.eq.1)THEN
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,pinew,thnew,qvnew)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! cm1r17:
            ! note:  nothing in pre-RK section should modify rho
            IF( abs(dt*thten1(i,j,k)).gt.tem .or.  &
                abs(dt*qten(i,j,k,nqv)).gt.qsmall )THEN
              thnew = tha(i,j,k)+dt*thten1(i,j,k)
              qvnew = qa(i,j,k,nqv)+dt*qten(i,j,k,nqv)
              pinew = (rho(i,j,k)*(th0(i,j,k)+thnew)*(rd+rv*qvnew)*rp00)**rddcv - pi0(i,j,k)
              ppten(i,j,k) = (pinew-ppi(i,j,k))*rdt
            ELSE
              ppten(i,j,k) = 0.0
            ENDIF
          enddo
          enddo
          enddo
          if( ptype.ge.1 )then
            ! use diabatic tendencies from last timestep as a good estimate:
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              ppten(i,j,k)=ppten(i,j,k)+qpten(i,j,k)
              thten1(i,j,k)=thten1(i,j,k)+qtten(i,j,k)
              qten(i,j,k,nqv)=qten(i,j,k,nqv)+qvten(i,j,k)
              qten(i,j,k,nqc)=qten(i,j,k,nqc)+qcten(i,j,k)
            enddo
            enddo
            enddo
          endif
        ELSE
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,pinew,thnew)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            !-----
            ! cm1r17:
            ! note:  nothing in pre-RK section should modify rho
            IF( abs(dt*thten1(i,j,k)).gt.tem )THEN
              thnew = tha(i,j,k)+dt*thten1(i,j,k)
              pinew = (rho(i,j,k)*(th0(i,j,k)+thnew)*rd*rp00)**rddcv - pi0(i,j,k)
              ppten(i,j,k) = (pinew-ppi(i,j,k))*rdt
            ELSE
              ppten(i,j,k)=0.0
            ENDIF
            !-----
          enddo
          enddo
          enddo
        ENDIF  ! endif for moist/dry
      ENDIF    ! endif for eqtset 1/2

    ENDIF  ps6

        if(timestats.ge.1) time_integ=time_integ+mytime()

!--------------------------------------------------------------------

    ! NOTE:

    ! cm1r20.1:  - moved microphysics to mp_driver (which is called from cm1.F)
    !            - moved message passing and equate to solve_finish.F


      end subroutine solve1


  END MODULE solve1_module
