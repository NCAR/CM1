  MODULE solve3_module

        ! solve3: finish primary dynamical solver
        ! (finish comms, update arrays, move surface)

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
  public :: solve3

  CONTAINS

      subroutine solve3(nstep,num_soil_layers,                       &
                   dt,dtlast,mtime,dbldt,                            &
                   xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,  &
                   yh,vh,rvh,yf,vf,rvf,                              &
                   rds,sigma,rdsf,sigmaf,                            &
                   zh,mh,rmh,c1,c2,zf,mf,rmf,wprof,ug,vg,            &
                   pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0,              &
                   qi0,rr0,rf0,rrf0,                                 &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,            &
                   rain,raincv,sws,svs,sps,srs,sgs,sus,shs,          &
                   dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,          &
                   divx,rho,rr,rf,prs,                               &
                   u0,rru,ua,u3d,uten,uten1,                         &
                   v0,rrv,va,v3d,vten,vten1,                         &
                   rrw,wa,w3d,wten,wten1,                            &
                   ppi,pp3d,ppx,phi1,                                &
                   tha,th3d,qpten,qtten,qvten,qcten,qa,q3d,tkea,tke3d, &
                   pta,pt3d,pdata,u2pt,v2pt,                         &
                   tsk,znt,zntmp,ust,f2d,hfx,qfx,qsfc,tslb,tmn,      &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_p2,reqs_p3, &
                   reqs_x,reqs_y,reqs_z,reqs_tk,reqs_q,reqs_t,       &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,          &
                   ww1,ww2,we1,we2,ws1,ws2,wn1,wn2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   p2w1,p2w2,p2e1,p2e2,p2s1,p2s2,p2n1,p2n2,          &
                   p3w1,p3w2,p3e1,p3e2,p3s1,p3s2,p3n1,p3n2,          &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   zw1,zw2,ze1,ze2,zs1,zs2,zn1,zn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                   ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,          &
                   qw31,qw32,qe31,qe32,qs31,qs32,qn31,qn32,          &
                   tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,                  &
                   tw1,tw2,te1,te2,ts1,ts2,tn1,tn2,                  &
                   udiag,vdiag,wdiag,pdiag,                          &
                   lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv,        &
                   dowriteout,dorad,getdbz,getvt,dotdwrite,          &
                   dotbud,doqbud,doubud,dovbud,dowbud,donudge,       &
                   doazimwrite,dorestart)
        ! end_solve3
      use input
      use constants
      use bc_module
      use comm_module
      use adv_routines, only : movesfc
      use misclibs
      use parcel_module, only : parcel_driver
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
      double precision, intent(in   ) :: mtime
      double precision, intent(inout) :: dbldt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(in),    dimension(kb:ke) :: wprof
      real, intent(inout), dimension(kb:ke) :: ug,vg
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,rth0,qv0,qc0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: qi0,rr0,rf0,rrf0
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,raincv,sws,svs,sps,srs,sgs,sus,shs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: divx,rho,rr,rf,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0,rru,ua,u3d,uten,uten1
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0,rrv,va,v3d,vten,vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten,wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppx
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(ib2pt:ie2pt,jb2pt:je2pt,kb2pt:ke2pt) :: u2pt,v2pt
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,zntmp,ust
      real, intent(in),    dimension(ib:ie,jb:je) :: f2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: hfx,qfx,qsfc
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_p2,reqs_p3,reqs_x,reqs_y,reqs_z,reqs_tk
      integer, intent(inout), dimension(rmp,numq) :: reqs_q
      integer, intent(inout), dimension(rmp,npt) :: reqs_t
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, intent(inout), dimension(jmp,kmp-1) :: ww1,ww2,we1,we2
      real, intent(inout), dimension(imp,kmp-1) :: ws1,ws2,wn1,wn2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(jmp,kmp) :: p2w1,p2w2,p2e1,p2e2
      real, intent(inout), dimension(imp,kmp) :: p2s1,p2s2,p2n1,p2n2
      real, intent(inout), dimension(jmp,kmp) :: p3w1,p3w2,p3e1,p3e2
      real, intent(inout), dimension(imp,kmp) :: p3s1,p3s2,p3n1,p3n2
      real, intent(inout), dimension(jmp,kmp) :: vw1,vw2,ve1,ve2
      real, intent(inout), dimension(imp,kmp) :: vs1,vs2,vn1,vn2
      real, intent(inout), dimension(jmp,kmp) :: zw1,zw2,ze1,ze2
      real, intent(inout), dimension(imp,kmp) :: zs1,zs2,zn1,zn2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: rw31,rw32,re31,re32
      real, intent(inout), dimension(imp,cmp,kmp)   :: rs31,rs32,rn31,rn32
      real, intent(inout), dimension(cmp,jmp,kmp,numq) :: qw31,qw32,qe31,qe32
      real, intent(inout), dimension(imp,cmp,kmp,numq) :: qs31,qs32,qn31,qn32
      real, intent(inout), dimension(cmp,jmp,kmt)   :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt)   :: tks1,tks2,tkn1,tkn2
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      real, intent(inout), dimension(cmp,jmp,kmp,npt) :: tw1,tw2,te1,te2
      real, intent(inout), dimension(imp,cmp,kmp,npt) :: ts1,ts2,tn1,tn2
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nudiag) :: udiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nvdiag) :: vdiag
      real, intent(inout) , dimension(ibdv:iedv,jbdv:jedv,kbdv:kedv,nwdiag) :: wdiag
      real, intent(inout) , dimension(ibdp:iedp,jbdp:jedp,kbdp:kedp,npdiag) :: pdiag
      real, intent(inout), dimension(lsnudge_maxlevels,lsnudge_nmax) :: lsnudge_u,lsnudge_v,lsnudge_th,lsnudge_qv
      logical, intent(in) :: dowriteout,dorad,dotdwrite,doazimwrite,dorestart
      logical, intent(inout) :: getdbz,getvt,dotbud,doqbud,doubud,dovbud,dowbud,donudge

!-----------------------------------------------------------------------
! Arrays and variables defined inside solve

      integer :: i,j,k,n,nrk,bflag,pdef,diffit,k1

      real :: delqv,delpi,delth,delt,fac,epsd,dheat,dz1,xs
      real :: foo1,foo2
      real :: dttmp,rtime,rdt,tem,tem0,tem1,tem2,thrad,prad
      real :: r1,r2,tnew,pnew,pinew,thnew,qvnew
      real :: gamm,aiu
      real :: qv,qli,cpli,cpm,cvm,qmax
      real :: tn,qn,taunudge

      double precision :: weps,afoo,bfoo,p0,p2

      logical :: get_time_avg,dopf

#ifdef MPI
      integer :: reqp1,reqp2
#endif

!--------------------------------------------------------------------


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Begin:  message passing

          call bcp(pp3d)
#ifdef MPI
          call comm_1s_start(pp3d,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_x)
#endif

          call bcs(th3d)
#ifdef MPI
          call comm_3s_start(th3d,rw31,rw32,re31,re32,   &
                                  rs31,rs32,rn31,rn32,reqs_y)
#endif

      IF( imoist.eq.1 )THEN
          DO n=1,numq
            call bcs(q3d(ib,jb,kb,n))
#ifdef MPI
            call comm_3s_start(q3d(ib,jb,kb,n)  &
                       ,qw31(1,1,1,n),qw32(1,1,1,n),qe31(1,1,1,n),qe32(1,1,1,n)     &
                       ,qs31(1,1,1,n),qs32(1,1,1,n),qn31(1,1,1,n),qn32(1,1,1,n)     &
                       ,reqs_q(1,n) )
#endif
          ENDDO
      ENDIF


      IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.6 .or. psolver.eq.7 )THEN
        ! 180212: moved this out of sound,sounde,soundcb

        if( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.7 )then

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppx(i,j,k)=pp3d(i,j,k)+ppx(i,j,k)
          enddo
          enddo
          enddo

        elseif( psolver.eq.6 )then

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppx(i,j,k)=phi1(i,j,k)+ppx(i,j,k)
          enddo
          enddo
          enddo

        endif

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcp(ppx)
#ifdef MPI
        call comm_1p_start(ppx,zw1,zw2,ze1,ze2,zs1,zs2,zn1,zn2,reqs_z)
#endif

      ENDIF


      call bcp(rho)
#ifdef MPI
      call comm_1s_start(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
#endif


!-----------------------------------------------------------------

!Done:  message passing
!-----------------------------------------------------------------
!  cm1r17:  diabatic tendencies for next timestep:

        IF( imoist.eq.1 .and. eqtset.eq.2 .and. ptype.ge.1 )THEN
          ! get diabatic tendencies (will be used in next timestep):
          rdt = 1.0/dt
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qpten(i,j,k)=(pp3d(i,j,k)-qpten(i,j,k))*rdt
            qtten(i,j,k)=(th3d(i,j,k)-qtten(i,j,k))*rdt
            qvten(i,j,k)=(q3d(i,j,k,nqv)-qvten(i,j,k))*rdt
            qcten(i,j,k)=(q3d(i,j,k,nqc)-qcten(i,j,k))*rdt
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
        ENDIF

!-----------------------------------------------------------------
!  Equate the two arrays

      IF( dot2p )THEN
        rdt = 1.0/dt
        ! 200413:  all tendencies except turbz term:
        do k=1,ntwk
        do j=1,nj+1
        do i=1,ni+1
          ! 211207:  all tendencies except 2pt term:
          u2pt(i,j,k) = (u3d(i,j,k)-ua(i,j,k))*rdt - u2pt(i,j,k)
          v2pt(i,j,k) = (v3d(i,j,k)-va(i,j,k))*rdt - v2pt(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_turb=time_turb+mytime()
      ENDIF

      if( do_lsnudge .and. donudge .and. do_lsnudge_u )then
        if( myid.eq.0 ) print *,'  applying large-scale nudging to base-state and reference u profiles '
        !-----
        ! update ug,u0:
        tem = lsnudgefac*dt/lsnudge_tau
        do k=1,nk
          ug(k) = ug(k)-tem*(umove+ug(k)-lsnudge_u(k,lsnudge_count))
          u0(1,1,k) = u0(1,1,k)-tem*(umove+u0(1,1,k)-lsnudge_u(k,lsnudge_count))
        enddo
        !-----
        do k=1,nk
        do j=jb,je
        do i=ib,ie+1
          u0(i,j,k) = u0(1,1,k)
        enddo
        enddo
        enddo
        do j=jb,je
        do i=ib,ie+1
          u0(i,j,0)    = cgs1*u0(i,j,1)+cgs2*u0(i,j,2)+cgs3*u0(i,j,3)
          u0(i,j,nk+1) = cgt1*u0(i,j,nk)+cgt2*u0(i,j,nk-1)+cgt3*u0(i,j,nk-2)
        enddo
        enddo
        !-----
        if(timestats.ge.1) time_misc=time_misc+mytime()
      endif
      if( do_lsnudge .and. donudge .and. do_lsnudge_v )then
        if( myid.eq.0 ) print *,'  applying large-scale nudging to base-state and reference v profiles '
        !-----
        ! update vg,v0:
        tem = lsnudgefac*dt/lsnudge_tau
        do k=1,nk
          vg(k) = vg(k)-tem*(vmove+vg(k)-lsnudge_v(k,lsnudge_count))
          v0(1,1,k) = v0(1,1,k)-tem*(vmove+v0(1,1,k)-lsnudge_v(k,lsnudge_count))
        enddo
        !-----
        do k=1,nk
        do j=jb,je+1
        do i=ib,ie
          v0(i,j,k) = v0(1,1,k)
        enddo
        enddo
        enddo
        do j=jb,je+1
        do i=ib,ie
          v0(i,j,0)    = cgs1*v0(i,j,1)+cgs2*v0(i,j,2)+cgs3*v0(i,j,3)
          v0(i,j,nk+1) = cgt1*v0(i,j,nk)+cgt2*v0(i,j,nk-1)+cgt3*v0(i,j,nk-2)
        enddo
        enddo
        !-----
        if(timestats.ge.1) time_misc=time_misc+mytime()
      endif

#ifdef MPI
      call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                           us31,us32,un31,un32,reqs_u)
      call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                           vs31,vs32,vn31,vn32,reqs_v)
      call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                           ws31,ws32,wn31,wn32,reqs_w)
#endif

#ifdef MPI
      if(   sgsmodel.eq.5 .or. sgsmodel.eq.6 )then
        call getcorneru(u3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcu2(u3d)
        call getcornerv(v3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcv2(v3d)
        call getcornerw(w3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcw2(w3d)
      endif
#endif

      if(terrain_flag)then
        call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
        call bc2d(w3d(ib,jb,1))
      endif

!----------
!  comms for parcels:

      IF(iprcl.eq.1)THEN
        ! cm1r18:  use velocities averaged over small time steps (for psolver=2,3,6,7)
        IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.6 .or. psolver.eq.7 )THEN
          ! 180713:  now saved in ten1 arrays
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            do i=1,ni+1
              rru(i,j,k)=uten1(i,j,k)
            enddo
            enddo
            IF(axisymm.eq.0)THEN
              ! Cartesian grid:
              do j=1,nj+1
              do i=1,ni
                rrv(i,j,k)=vten1(i,j,k)
              enddo
              enddo
            ENDIF
            IF(k.gt.1)THEN
              do j=1,nj
              do i=1,ni
                rrw(i,j,k)=wten1(i,j,k)
              enddo
              enddo
            ENDIF
          enddo
        ELSE
          ! psolver=1,4,5:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            do i=1,ni+1
              rru(i,j,k)=u3d(i,j,k)
            enddo
            enddo
            IF(axisymm.eq.0)THEN
              ! Cartesian grid:
              do j=1,nj+1
              do i=1,ni
                rrv(i,j,k)=v3d(i,j,k)
              enddo
              enddo
            ENDIF
            IF(k.gt.1)THEN
              do j=1,nj
              do i=1,ni
                rrw(i,j,k)=w3d(i,j,k)
              enddo
              enddo
            ENDIF
          enddo
        ENDIF
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
        ! bc/comms:
        call bcu(rru)
#ifdef MPI
        call comm_3u_start(rru,uw31,uw32,ue31,ue32,us31,us32,un31,un32,reqs_u)
#endif
        call bcv(rrv)
#ifdef MPI
        call comm_3v_start(rrv,vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,reqs_v)
#endif
        call bcw(rrw,1)
#ifdef MPI
        call comm_3w_start(rrw,ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,reqs_w)
#endif
      ENDIF

!----------
!  Diagnostics (infrequent):

      IF( doubud .and. nutk.ge.1 )THEN
        if( myid.eq.0 ) print *,'  save utk for domaindiag '
        rdt = 1.0/dt
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          udiag(i,j,k,nutk) = 0.5*(ua(i,j,k)+u3d(i,j,k))
        enddo
        enddo
        enddo
      ENDIF
      IF( dovbud .and. nvtk.ge.1 )THEN
        if( myid.eq.0 ) print *,'  save vtk for domaindiag '
        rdt = 1.0/dt
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vdiag(i,j,k,nvtk) = 0.5*(va(i,j,k)+v3d(i,j,k))
        enddo
        enddo
        enddo
      ENDIF
      IF( dowbud .and. nwtk.ge.1 )THEN
        if( myid.eq.0 ) print *,'  save wtk for domaindiag '
        rdt = 1.0/dt
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wdiag(i,j,k,nwtk) = 0.5*(wa(i,j,k)+w3d(i,j,k))
        enddo
        enddo
        enddo
      ENDIF

!----------

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=0,nj+2
      do i=0,ni+2
        ua(i,j,k)=u3d(i,j,k)
        va(i,j,k)=v3d(i,j,k)
        wa(i,j,k)=w3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

!----------

      if( idoles .and. iusetke )then
#ifdef MPI
        call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                               tks1,tks2,tkn1,tkn2,reqs_tk)
#endif
        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          tkea(i,j,k)=tke3d(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()
      endif

!----------

      if(iptra.eq.1)then
        do n=1,npt
#ifdef MPI
          call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                reqs_t(1,n))
#endif
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            pta(i,j,k,n)=pt3d(i,j,k,n)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_integ=time_integ+mytime()
        enddo
      endif

!----------

#ifdef MPI
      call comm_1s_end(pp3d,vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,reqs_x)
      call comm_3s_end(th3d,rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,reqs_y)
#endif

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
      do k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        ppi(i,j,k)=pp3d(i,j,k)
        tha(i,j,k)=th3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

!----------

      if(imoist.eq.1)then
        DO n=1,numq
#ifdef MPI
          call comm_3s_end(q3d(ib,jb,kb,n)  &
                     ,qw31(1,1,1,n),qw32(1,1,1,n),qe31(1,1,1,n),qe32(1,1,1,n)     &
                     ,qs31(1,1,1,n),qs32(1,1,1,n),qn31(1,1,1,n),qn32(1,1,1,n)     &
                     ,reqs_q(1,n) )
#endif
          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            qa(i,j,k,n)=q3d(i,j,k,n)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_integ=time_integ+mytime()
        ENDDO
      endif

!----------

#ifdef MPI
        IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.6 .or. psolver.eq.7 )THEN
          call comm_1p_end1(ppx,zw1,zw2,ze1,ze2,zs1,zs2,zn1,zn2,reqs_z)
          call comm_1p_end2(reqs_z)
        ENDIF
#endif

#ifdef MPI
        call comm_1s_end(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_s)
        call getcorner(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
        call bcs2(rho)
#endif

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
        do k=2,nk
        do j=0,nj+1
        do i=0,ni+1
          rr(i,j,k) = 1.0/rho(i,j,k)
          rf(i,j,k) = (c1(i,j,k)*rho(i,j,k-1)+c2(i,j,k)*rho(i,j,k))
        enddo
        enddo
        enddo

        ! meh 1 !
      !$omp parallel do default(shared)  &
      !$omp private(i,j)
        do j=0,nj+1
        do i=0,ni+1
          rr(i,j,1) = 1.0/rho(i,j,1)
          ! cm1r17, 2nd-order extrapolation:
          rf(i,j,1) = cgs1*rho(i,j,1)+cgs2*rho(i,j,2)+cgs3*rho(i,j,3)
          rf(i,j,nk+1) = cgt1*rho(i,j,nk)+cgt2*rho(i,j,nk-1)+cgt3*rho(i,j,nk-2)
        enddo
        enddo
        if(timestats.ge.1) time_prsrho=time_prsrho+mytime()

!-----------------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  Update parcel locations  ccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF(iprcl.eq.1)THEN
        !  get corner info, ghost zone data, etc:
        !  (may not parallelize correctly if this is not done)
#ifdef MPI
        call comm_3u_end( rru,uw31,uw32,ue31,ue32,us31,us32,un31,un32,reqs_u)
        call comm_3v_end( rrv,vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,reqs_v)
        call comm_3w_end( rrw,ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,reqs_w)
        call getcorneru3(rru,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                             s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
        call getcornerv3(rrv,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                             s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
        call getcornerw3(rrw,n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2)
        call bcu2(rru)
        call bcv2(rrv)
        call bcw2(rrw)
#endif

        prltrn:  &
        if(terrain_flag)then
          ! 180713:  get sigma-dot

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k)
          DO k=1,nk
            do j=0,nj+1
            do i=0,ni+2
              rru(i,j,k)=rru(i,j,k)*rgzu(i,j)
            enddo
            enddo
            do j=0,nj+2
            do i=0,ni+1
              rrv(i,j,k)=rrv(i,j,k)*rgzv(i,j)
            enddo
            enddo
          ENDDO

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,r1,r2)
          DO k=1,nk
            IF(k.eq.1)THEN
              do j=0,nj+1
              do i=0,ni+1
                rrw(i,j,   1) = 0.0
                rrw(i,j,nk+1) = 0.0
              enddo
              enddo
            ELSE
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              r1 = 0.5*r1
              r2 = 0.5*r2
              do j=0,nj+1
              do i=0,ni+1
                rrw(i,j,k)=rrw(i,j,k)                              &
                          +( ( r2*(rru(i,j,k  )+rru(i+1,j,k  ))               &
                              +r1*(rru(i,j,k-1)+rru(i+1,j,k-1)) )*dzdx(i,j)   &
                            +( r2*(rrv(i,j,k  )+rrv(i,j+1,k  ))               &
                              +r1*(rrv(i,j,k-1)+rrv(i,j+1,k-1)) )*dzdy(i,j)   &
                           )*(sigmaf(k)-zt)*gz(i,j)*rzt
              enddo
              enddo
            ENDIF
          ENDDO

          !$omp parallel do default(shared)  &
          !$omp private(i,j,k,r1,r2)
          do k=1,nk
          do j=0,nj+2
          do i=0,ni+2
            rru(i,j,k) = rru(i,j,k)*gzu(i,j)
            rrv(i,j,k) = rrv(i,j,k)*gzv(i,j)
            rrw(i,j,k) = rrw(i,j,k)*gz(i,j)
          enddo
          enddo
          enddo

        endif  prltrn

        call     parcel_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,zs,    &
                               sigma,sigmaf,zntmp,rho,rru,rrv,rrw,pdata)
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
      ENDIF


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   All done   cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      ! cm1r19.6  (preliminary ... needs more testing)
      IF( imove.eq.1 )THEN
      IF( (sfcmodel.eq.2) .or. (sfcmodel.eq.3) .or. (sfcmodel.eq.4) )THEN

        weps = 0.0001*epsilon
        call movesfc(0.0,dt,weps,uh,vh,znt(ib,jb),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))

        weps = 1.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,ust(ib,jb),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))

        weps = 100.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,hfx(ibl,jbl),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))

        weps = 1.0e-6*epsilon
        call movesfc(0.0,dt,weps,uh,vh,qfx(ibl,jbl),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
        call movesfc(0.0,dt,weps,uh,vh,qsfc(ibl,jbl),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))

        weps = 200.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,tsk(ib,jb),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
        call movesfc(0.0,dt,weps,uh,vh,tmn(ibl,jbl),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
        do n=1,num_soil_layers
        call movesfc(0.0,dt,weps,uh,vh,tslb(ibl,jbl,n),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
        enddo

      ENDIF
      ENDIF


!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!      if(timestats.ge.1) time_mpb=time_mpb+mytime()
!!!#endif

!--------------------------------------------------------------------
!  Calculate surface "swaths."  Move surface (if necessary). 
!--------------------------------------------------------------------

    IF( output_rain.eq.1 )THEN

      if(imove.eq.1.and.imoist.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,rain(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

    IF( output_rain.eq.1 .and. cuparam.ge.1 )THEN

      if(imove.eq.1.and.imoist.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,raincv(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!--------------------------------------------------------------------
! Maximum horizontal wind speed at lowest model level: 
! (include domain movement in calculation)

    IF( output_sws.eq.1 )THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = sqrt( (umove+0.5*(ua(i,j,1)+ua(i+1,j,1)))**2    &
                   +(vmove+0.5*(va(i,j,1)+va(i,j+1,1)))**2 ) 
        do n=1,nrain
          sws(i,j,n)=max(sws(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,sws(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!--------------------------------------------------------------------
!  Maximum vertical vorticity at lowest model level:

  IF( output_svs.eq.1 )THEN

  IF(axisymm.eq.0)THEN
    IF(.not.terrain_flag)THEN
      ! Cartesian grid, without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)-va(i-1,j,1))*rdx*uf(i)   &
             -(ua(i,j,1)-ua(i,j-1,1))*rdy*vf(j)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
    ELSE
      ! Cartesian grid, with terrain:
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
        ! interior:
        do k=2,2
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+2
          dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
          dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
        enddo
        enddo
      enddo
      k = 1
!$omp parallel do default(shared)  &
!$omp private(i,j,n,r1,tem)
      do j=1,nj+1
      do i=1,ni+1
        r1 = zt/(zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))
        tem = ( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
               +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                     -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                    )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )        &
             -( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
               +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                     -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                    )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
    ENDIF
  ELSE
      ! Axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)*xh(i)-va(i-1,j,1)*xh(i-1))*rdx*uf(i)*rxf(i)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
  ENDIF
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 1.0*epsilon
        call movesfc(-1000.0,dt,weps,uh,vh,svs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

  ENDIF

!--------------------------------------------------------------------
!  Minimum pressure perturbation at lowest model level:

  IF( output_sps.eq.1 )THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = prs(i,j,1)-prs0(i,j,1)
        do n=1,nrain
          sps(i,j,n)=min(sps(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 1000.0*epsilon
        call movesfc(-200000.0,dt,weps,uh,vh,sps(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

  ENDIF

!--------------------------------------------------------------------
!  Maximum rainwater mixing ratio (qr) at lowest model level:

  IF( output_srs.eq.1 )THEN

    IF(imoist.eq.1.and.nqr.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqr)
        do n=1,nrain
          srs(i,j,n)=max(srs(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(0.0,dt,weps,uh,vh,srs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif
    ENDIF

  ENDIF

!--------------------------------------------------------------------
!  Maximum graupel/hail mixing ratio (qg) at lowest model level:

  IF( output_sgs.eq.1 )THEN

    IF(imoist.eq.1.and.nqg.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqg)
        do n=1,nrain
          sgs(i,j,n)=max(sgs(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(0.0,dt,weps,uh,vh,sgs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif
    ENDIF

  ENDIF

!--------------------------------------------------------------------

  IF( output_sus.eq.1 )THEN

      ! get height AGL:
      if( .not. terrain_flag )then
        ! without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)
          wten(i,j,k) = zf(i,j,k)
        enddo
        enddo
        enddo
      else
        ! get height AGL:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)-zs(i,j)
          wten(i,j,k) = zf(i,j,k)-zs(i,j)
        enddo
        enddo
        enddo
      endif

!--------------------------------------------------------------------
!  Maximum updraft velocity (w) at 5 km AGL:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tem)
      do j=1,nj
      do i=1,ni
        k = 2
        ! wten is height AGL:
        do while( wten(i,j,k).lt.5000.0 .and. k.lt.nk )
          k = k + 1
        enddo
        tem = w3d(i,j,k)
        do n=1,nrain
          sus(i,j,n)=max(sus(i,j,n),tem)
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 10.0*epsilon
        call movesfc(-1000.0,dt,weps,uh,vh,sus(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!--------------------------------------------------------------------
!  Maximum integrated updraft helicity:

    IF( output_shs.eq.1 )THEN

      ! dum3 is zh (agl), wten is zf (agl)
      call calcuh(uf,vf,dum3,wten,ua,va,wa,dum1(ib,jb,1),dum2,dum5,dum6, &
                  zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf)
!$omp parallel do default(shared)  &
!$omp private(i,j,n)
      do j=1,nj
      do i=1,ni
        do n=1,nrain
          shs(i,j,n)=max(shs(i,j,n),dum1(i,j,1))
        enddo
      enddo
      enddo
      if(timestats.ge.1) time_swath=time_swath+mytime()

      if(imove.eq.1)then
        weps = 100.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,shs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!  Done with "swaths"
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!--------------------------------------------------------------------
      ! all done

      end subroutine solve3


  END MODULE solve3_module
