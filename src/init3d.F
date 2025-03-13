  MODULE init3d_module

  implicit none

  private
  public :: init3d,getset

  CONTAINS

!-----------------------------------------------------------------------------
!
!  subroutine INIT3D:  Specify initial conditions for CM1. 
!
!    Note:  init3d assumes that the base state arrays (th0,pi0,rho0,etc)
!           have already been specified.
!
!    Note:  In general, the code below assumes that the 3d initial state
!           for CM1 is simply the base-state conditions (from BASE) plus
!           specified perturbations (defined in INIT3D).  Hence, most 
!           code below calculates perturbation potential temp., pressure,
!           etc  (tha,ppi arrays).
!
!           See the documentation near the top of solve.F and/or Table 1
!           in the document "The governing equations for CM1"
!           (http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf)
!           for more details about variables and arrays in CM1
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
!        http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!-----------------------------------------------------------------------------

      subroutine init3d(xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                        xfref,xhref,yfref,yhref,sigma,c1,c2,gz,zs,        &
                        zh,mh,rmh,zf,mf,rmf,rho0s,pi0s,prs0s,             &
                        pi0,prs0,rho0,thv0,th0,rth0,qv0,                  &
                        u0,v0,qc0,qi0,rr0,rf0,rrf0,                       &
                        rain,sws,svs,sps,srs,sgs,sus,shs,                 &
                        thflux,qvflux,cd,ch,cq,f2d,                       &
                        dum1,dum2,dum3,dum4,divx,rho,prs,                 &
                        rru,ua,u3d,uten,uten1,rrv,va,v3d,vten,vten1,      &
                        rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,        &
                        tha,th3d,thten,thten1,qa,q3d,qten,                &
                        kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                        pta,pt3d,ptten,                                   &
                        pdata,cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans)

      use input
      use constants
      use misclibs
      use cm1libs , only : rslf,rsif
      use bc_module
      use module_mp_nssl_2mom, only: ccn, lccn, irenuc
      use poiss_module
      use parcel_module , only : getparcelzs
#ifdef MPI
      use mpi
#endif

      implicit none
 
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf,rxf,uf,ruf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(1-ngxy:nx+ngxy+1) :: xfref,xhref
      real, intent(in), dimension(1-ngxy:ny+ngxy+1) :: yfref,yhref
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,rho0,thv0,th0,rth0,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: qc0,qi0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,cd,ch,cq,f2d
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4
      real, dimension(ib:ie,jb:je,kb:ke) :: divx,rho,prs
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,rru,ua,u3d,uten,uten1
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,rrv,va,v3d,vten,vten1
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten,wten1
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten,sten
      real, dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thten1
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d,ptten
      real, dimension(nparcels,npvals) :: pdata
      real, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, dimension(kpb:kpe) :: cfa,cfc,ad1,ad2
      complex, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,lgbth,lgbph
      complex, dimension(ipb:ipe,jpb:jpe) :: rhs,trans
 
!-----------------------------------------------------------------------

      integer i,j,k,kk,l,n,nn,nbub,nloop
      integer ic,jc,ifoo,jfoo
      real ric,rjc
      real xc,yc,zc,bhrad,bvrad,bptpert,beta,tmp,zdep
      real thvnew(nk),pinew(nk)
      real thl,ql,qt,th1,t1,ql2,rm,cpm,v1,v2,th2

      real, dimension(:), allocatable :: rref,vt
      real, dimension(:,:), allocatable :: vref,piref,pref,thref,thvref,qvref,rhref
      real :: rmax,vmax
      double precision :: frac,angle
      real :: r0,zdd,dd2,dd1,rr,diff,xref,yref,xsave,ysave,rrmax
      real :: mult,nominal_dx
      integer :: ival,ni1,ni2,ni3
      integer :: i1,i2,ii,jj,nref
      real :: r0_re87,rmax_re87,vmax_re87

      real rmin,foo1,foo2,umax,umin,vmin
      real :: rand,rand2,amplitude
      integer, dimension(:), allocatable :: sand
      double precision :: dpi

      logical :: setppi,maintain_rh

      real :: rm1,rm2,rm3,rdc,w2,w3,v3

      logical :: rhmods
      real :: rhmax,rmid,rscale,zmid,zscale,rh0,rh1
      real :: tnew,fliq,fice,qsi,qsw,qvs

      logical :: dobsmod
      real :: th_sfc,pi_sfc,pifoo,prsfoo,thfoo,t0
      real :: prs_sfc,t_sfc,hs
      integer :: tctype
      real :: tem
      double precision :: tem1,tem2,temd
      logical :: doit


      ! Note:  By default, CM1 uses the same set of psuedorandom numbers 
      !        every time it is called.  To use truly random perturbations 
      !        (ie, different sets of random perturbations every time CM1 
      !        is used), then change this varliable to .true.
      logical, parameter :: use_truly_random_pert  =  .false.

!--------------------------

      if(dowr) write(outfile,*) 'Inside INIT3D'
      if(dowr) write(outfile,*)

      convinit = 0
      wnudge = 0
      setppi = .true.
      maintain_rh = .false.

!------------------------------------------------------------------
!  Initialize surface swath arrays:

      do n=1,nrain
      do j=jb,je
      do i=ib,ie
        ! these are all positive-definite, so set initial value to zero:
        rain(i,j,n)=0.0
        sws(i,j,n)=0.0
        srs(i,j,n)=0.0
        sgs(i,j,n)=0.0
        shs(i,j,n)=0.0
        ! for sps, we want to get a MINIMUM value at the surface, so...
        ! set sps to an absurdly large number:
        sps(i,j,n)=200000.0
        ! svs and sus can be negative or positive, 
        ! but we want to get a MAXIMUM value, so...
        ! set svs and sus to an absurdly low (negative) number:
        svs(i,j,n)=-1000.0
        sus(i,j,n)=-1000.0
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!  Set winds to base-state values:

      do k=kb,ke
      do j=jb,je
      do i=ib,ie+1
        ua(i,j,k)=u0(i,j,k)
        u3d(i,j,k)=u0(i,j,k)
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je+1
      do i=ib,ie
        va(i,j,k)=v0(i,j,k)
        v3d(i,j,k)=v0(i,j,k)
      enddo
      enddo
      enddo


      IF( idoles .and. iusetke .and. ( testcase.eq.3 .or. testcase.eq.4 .or.  testcase.eq.5 .or. testcase.eq.9 .or. testcase.eq.11 ) )THEN
          IF( testcase.eq.3 )THEN
            ! shallow Cu case (Siebesma at al, 2003, JAS)
            do k=1,nk+1
            do j=jbt,jet
            do i=ibt,iet
              tkea(i,j,k)=max(0.0,1.0-zf(i,j,k)/3000.0)
              tke3d(i,j,k)=tkea(i,j,k)
            enddo
            enddo
            enddo
          ELSEIF( testcase.eq.9 .or. testcase.eq.11 )THEN
            ! stable boundary layer (Beare et al. 2006, BLM)
            do k=1,nk+1
            do j=jbt,jet
            do i=ibt,iet
              if( zh(i,j,k).lt.255.0 )then
                tkea(i,j,k)=0.4*((1.0-zf(i,j,k)/255.0)**3)
              else
                tkea(i,j,k)=0.0
              endif
              tke3d(i,j,k)=tkea(i,j,k)
            enddo
            enddo
            enddo
          ELSE
            tkea = 1.0
            tke3d = 1.0
          ENDIF
          call bcw(tkea,1)
          call bcw(tke3d,1)
      ENDIF


!-----------------------------------------------------------------------

    IF ( ( ptype .ge. 26 .and. ptype.le.29 ) .and. lccn .gt. 0 ) THEN
! initialize CCN concentrations as constant mixing ratios througout the domain
      do k=kbm,kem
      do j=jbm,jem
      do i=ibm,iem
       IF ( irenuc < 5 ) THEN ! number of unactivated ccn
       qa (i,j,k,lccn-1) = ccn/1.225
       q3d(i,j,k,lccn-1) = ccn/1.225
       ELSE ! treat as number of activated ccn
       qa (i,j,k,lccn-1) = 0.0
       q3d(i,j,k,lccn-1) = 0.0
       ENDIF
      enddo
      enddo
      enddo
    ENDIF

!-----------------------------------------------------------------------
!  Set qv to base state value:

    im1:  &
    IF(imoist.eq.1)THEN

      do k=kbm,kem
      do j=jbm,jem
      do i=ibm,iem
        qa(i,j,k,nqv)=qv0(i,j,k)
      enddo
      enddo
      enddo

!---- This is here to ensure that certain idealized cases work ----

    pt1:  &
    IF( ptype.ge.1 )THEN

      IF( (isnd.eq.4 .or. isnd.eq.9 .or. isnd.eq.10 .or. isnd.eq.11 .or. isnd.eq.15) )THEN

        do k=kbm,kem
        do j=jbm,jem
        do i=ibm,iem
          qa(i,j,k,nqc)=qc0(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

      IF( (isnd.eq.4 .or. isnd.eq.9 .or. isnd.eq.10) .and. iice.eq.1 )THEN

        do k=kbm,kem
        do j=jbm,jem
        do i=ibm,iem
          qa(i,j,k,nqi)=qi0(i,j,k)
        enddo
        enddo
        enddo

      ENDIF

    ENDIF  pt1

    ENDIF  im1

!-----

    IF(iptra.eq.1)THEN
      ! define concentrations for passive fluid tracers here:
      do n=1,npt
      do k=kbp,kep
      do j=jbp,jep
      do i=ibp,iep
        if(n.eq.1)then
          pta(i,j,k,n)=0.0
          if(zh(i,j,k).lt.3000.0) pta(i,j,k,n)=0.001
        endif
        if(n.eq.2)then
          pta(i,j,k,n)=0.0
          if(zh(i,j,k).gt.3000.0.and.zh(i,j,k).lt.6000.0) pta(i,j,k,n)=0.001
        endif
        if(n.eq.3)then
          pta(i,j,k,n)=0.0
          if(zh(i,j,k).gt.6000.0.and.zh(i,j,k).lt.9000.0) pta(i,j,k,n)=0.001
        endif
      enddo
      enddo
      enddo
      enddo
    ENDIF

!-----
!  parcel info:

      IF(iprcl.eq.1)THEN
        ! define initial locations of parcels here:
        !   pdata(*,prx) = x location (m)
        !   pdata(*,pry) = y location (m)
        !   pdata(*,prz) = z location (m ASL)

        ! initialize to really small number (so we can use the allreduce command below)
        do n=1,nparcels
          pdata(n,prx) = -1.0e30
          pdata(n,pry) = -1.0e30
          pdata(n,prz) = -1.0e30
        enddo

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Parcels ! '
        if(dowr) write(outfile,*) '  npvals,nparcels = ',npvals,nparcels
        if(dowr) write(outfile,*) '  Initial parcel locations (x,y,z):'
        n = 0
        do k=1,10
        do j=1,60
        do i=1,60
          n = n + 1
          if(n.gt.nparcels)then
            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) ' You are trying to define too many parcels'
            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) ' Increase the value of nparcels in namelist.input'
            if(dowr) write(outfile,*)
            call stopcm1
          endif
          pdata(n,prx) = minx + 2000.0*(i-1)
          pdata(n,pry) = miny + 2000.0*(j-1)
          pdata(n,prz) = zh(1,1,1) + 1000.0*(k-1)
!!!          if(dowr) write(outfile,*) n,pdata(n,prx),pdata(n,pry),pdata(n,prz)
        enddo
        enddo
        enddo
        if(dowr) write(outfile,*)

#ifdef MPI
        ! this ensures that every processor has all parcel locations:
        call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,1),3*nparcels,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

        IF(axisymm.eq.1.or.ny.eq.1)THEN
          ! 170719,  for 2d setup (x,z), fix all y values:
          DO n=1,nparcels
            pdata(n,pry) = 0.0
          ENDDO
        ENDIF
        IF(nx.eq.1)THEN
          ! 170719,  for 2d setup (y,z), fix all x values:
          DO n=1,nparcels
            pdata(n,prx) = 0.0
          ENDDO
        ENDIF

      ENDIF

!-----------------------------------------------------------------------
!  initialize random number generator:

      IF( use_truly_random_pert )THEN

        ! randomly reinitializes the pseudorandom number generator:
        call reinit_random_seed

      ELSE

        ! generate same set of pseudorandom numbers every time:
        ! (default for cm1)
        !-----------------------------------!
        !----- don't change this code) -----!
        ! initialize the random number generator
        call random_seed(size=k)
        k = max(2,k)
        if(dowr) write(outfile,*) '  seed_size = ',k
        allocate( sand(k) )
        do n=1,k
          sand(n) = nint( 2.0e9*(2.0*float(n-1)/float(k-1)-1.0) )
        enddo
        call random_seed(put=sand(1:k))
        call random_number(rand)
        if(dowr) write(outfile,*) '  rand-1 = ',rand
        deallocate( sand )
        !----- don't change this code) -----!
        !-----------------------------------!

      ENDIF

!-----------------------------------------------------------------------
!  iinit = 1
!  Warm bubble
!  reference:

      IF(iinit.eq.1)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Warm bubble'
        if(dowr) write(outfile,*)

        ric     =  centerx  ! center of bubble in x-direction (m)
        rjc     =  centery  ! center of bubble in y-direction (m)
        zc      =   1400.0  ! height of center of bubble above ground (m)
        bhrad   =  10000.0  ! horizontal radius of bubble (m)
        bvrad   =   1400.0  ! vertical radius of bubble (m)
        bptpert =      1.0  ! max potential temp perturbation (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .false.

        do k=1,nk
        do j=1,nj
        do i=1,ni
          beta=sqrt(                             &
                    ((xh(i)-ric)/bhrad)**2       &
                   +((yh(j)-rjc)/bhrad)**2       &
                   +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            tha(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)
          else
            tha(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iinit = 2
!  Cold pool (dam break style)
!  reference:  

      ELSEIF(iinit.eq.2)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Cold pool .... periodic in N-S'
        if(dowr) write(outfile,*)

        ric      =   centerx   ! eastern edge of cold pool
        zdep     =    2500.0   ! depth of cold pool (m)
        bptpert  =      -6.0   ! max temp perturbation at sfc (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .true.

        do k=1,nk
        do j=1,nj
        do i=1,ni
          if( (xh(i).le.ric).and.(zh(i,j,k).lt.zdep) )then
            tha(i,j,k)=bptpert*(zdep-zh(i,j,k))/zdep
          else
            tha(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iinit = 3
!  Line of warm bubbles
!  reference:  

      ELSEIF(iinit.eq.3)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Line of warm bubbles'
        if(dowr) write(outfile,*)
 
        nbub    =      3     ! number of warm bubbles
        ric     =  30000.0   ! center of bubble in x-direction (m)
        zc      =   1400.0   ! height of center of bubble above ground (m)
        bhrad   =  10000.0   ! horizontal radius of bubble (m)
        bvrad   =   1400.0   ! vertical radius of bubble (m)
        bptpert =      2.0   ! max potential temp perturbation (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .false.

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          tha(i,j,k)=0.0
        enddo
        enddo
        enddo

        do n=1,nbub

          if(n.eq.1) rjc=  3000.0
          if(n.eq.2) rjc= 33000.0
          if(n.eq.3) rjc= 63000.0

          if(dowr) write(outfile,*) '  ric,rjc=',n,ric,rjc
 
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            beta=sqrt(                        &
                    ((xh(i)-ric)/bhrad)**2    &
                   +((yh(j)-rjc)/bhrad)**2    &
                   +((zh(i,j,k)-zc)/bvrad)**2)
            if(beta.lt.1.0)then
              tha(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)
            else
              tha(i,j,k)=max(0.0,tha(i,j,k))
            endif
          enddo
          enddo
          enddo

        enddo


!-----------------------------------------------------------------------
!  iinit = 4
!  moist bubble for moist benchmark
!  reference:  Bryan and Fritsch, 2002, MWR, 130, 2917-2928.

      ELSEIF(iinit.eq.4)THEN

        ! parameters for dry counterpart bubble

        ric      =      0.0       ! x-location of bubble center (m)
        zc       =   2000.0       ! z-location of bubble center (m)
        bhrad    =   2000.0       ! horizontal radius of bubble (m)
        bvrad    =   2000.0       ! vertical radius of bubble (m)
        bptpert  =      2.0       ! maximum potential temp. pert. (K)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          beta=sqrt( ((xh(i)-ric)/bhrad)**2    &
                    +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            dum1(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)
          else
            dum1(i,j,k)=0.
          endif
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          tha(i,j,k)=0.
          ppi(i,j,k)=0.
          dum2(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
        enddo
        enddo
        enddo

        do nn=1,30
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            qa(i,j,k,nqv)=dum2(i,j,k)*rslf(prs0(i,j,k),(th0(i,j,k)+tha(i,j,k))*pi0(i,j,k))
          enddo
          enddo
          enddo

          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            qa(i,j,k,nqc)=max(qt_mb-qa(i,j,k,nqv),0.0)
          enddo
          enddo
          enddo

          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            tha(i,j,k)=( (dum1(i,j,k)/300.)+(1.0+qt_mb)/(1.0+qa(i,j,k,nqv)) )  &
               *thv0(i,j,k)*(1.0+qa(i,j,k,nqv))/(1.0+reps*qa(i,j,k,nqv)) - th0(i,j,k)
            if(abs(tha(i,j,k)).lt.1.e-4) tha(i,j,k)=0.
          enddo
          enddo
          enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qa(i,j,k,nqv)=rslf(prs0(i,j,k),(th0(i,j,k)+tha(i,j,k))*pi0(i,j,k))
          qa(i,j,k,nqc)=max(qt_mb-qa(i,j,k,nqv),0.0)
        enddo
        enddo
        enddo

!-----------------------------------------------------------------
!  iinit = 5
!  density current sim

      ELSEIF(iinit.eq.5)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Cold pool (elipse, following Straka)'
        if(dowr) write(outfile,*)

        ric     =     0.0
        rjc     =     0.0
        zc      =  3000.0
        bhrad   =  4000.0
        bvrad   =  2000.0
        bptpert =   -15.0

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          beta=sqrt(                           &
                     ((xh(i)-ric)/bhrad)**2    &
!!!                    +((yh(j)-rjc)/bhrad)**2    &
                    +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            dum1(i,j,k)=bptpert*(cos(pi*beta)+1.0)*0.5
          else
            dum1(i,j,k)=0.0
          endif
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          tmp=(th0(i,j,k)*pi0(i,j,k))+dum1(i,j,k)
          tha(i,j,k)=tmp/pi0(i,j,k)-th0(i,j,k)
          if(abs(tha(i,j,k)).lt.1.e-4) tha(i,j,k)=0.0
          ppi(i,j,k)=0.0
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  Analytic tropical cyclone-like vortex

      ELSEIF(iinit.eq.7)THEN

        ! for 3d runs, place TC in center of domain, by default:
        xref  =  centerx
        yref  =  centery

        ! TC center could be placed elsewhere.  Use this code:
!!!        xref  =   200000.0    ! x location of TC center (m)
!!!        yref  =  -100000.0    ! y location of TC center (m)


        ! type of analytic vortex:
        !    1 = Rotunno and Emanuel (1987)
        !    2 = modified Rankine vortex  ( this is now default, starting in cm1r20.1 )
        tctype  =  2



        !      Default analytic vortex since cm1r20.1:
        ! for tctype=2, parameters for modified Rankine vortex:
        vmax =   15.0         ! V max (m/s)
        rm1  =  75000.0       ! 1st radius (m): location of Vmax
        rm2  = 200000.0       ! 2nd radius (m): max radius for specified radial decay rate
        rm3  = 500000.0       ! 3rd radius (m): location where V=0
        rdc  =  -0.35         ! radial decay rate (nondimensional) between rm1 and rm2
                              !   ( where v = Vmax (r/rm1)^rdc )


        ! for tctype=1, parameters for Rotunno-Emanuel (1987) analytic vortex:
        r0_re87     =   412500.0
        rmax_re87   =    82500.0
        vmax_re87   =       15.0


        ! for all tctypes:
        zdd         =    15000.0     ! height (ASL) where V=0


        !  170502, experimental code for moist initial vortex:
        rhmods  =   .false.    !  use this code?
        rhmax   =    0.98      !  max value of relative humidity in core
        rmid    =  300000.0    !  radius (m) for middle of tanh function in r
        rscale  =  100000.0    !  radial scale (m) for tanh function in r
        zmid    =   18000.0    !  height (m) for middle of tanh function in z
        zscale  =    2000.0    !  vertical scale (m) for tanh function in z


      !----------------------------------------------------------------!
      ! do not change anything below here !

        IF(axisymm.eq.1)THEN
          ! axisymmetric grid ... use actual grid:

          nref = nx
          xref = 0.0
          yref = 0.0

          allocate( rref(nref) )
          rref=0.0

          do i=1,nref
            rref(i) = 0.5*(xfref(i)+xfref(i+1))
          enddo
          rrmax = rref(nref)

        ELSE
          ! 3d Cartesian grid:

          nref = nx/2 + 3

          allocate( rref(nref) )
          rref=0.0

          do i=1,nref
            rref(i) = xhref( nx/2 + i )
          enddo

          rrmax = rref(nref)

          if( myid.eq.0 ) print *,'  xref,yref     = ',xref,yref
          if( myid.eq.0 ) print *,'  rrmax,nref = ',rrmax,nref

        ENDIF


        allocate(    vt(nref) )
        allocate(  vref(nref,0:nk+1))
        allocate( piref(nref,0:nk+1))
        allocate(  pref(nref,0:nk+1))
        allocate( thref(nref,0:nk+1))
        allocate(thvref(nref,0:nk+1))
        allocate( qvref(nref,0:nk+1))
        allocate( rhref(nref,0:nk+1))

            vt=0.0
          vref=0.0
         piref=0.0
          pref=0.0
         thref=0.0
        thvref=0.0
         qvref=0.0
         rhref=0.0

        IF(ibalance.ne.0)THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' Please use ibalance = 0 with iinit=7'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ... stopping inside init3d ... '
          if(dowr) write(outfile,*)
          call stopcm1
        ENDIF
        IF(terrain_flag)THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' iinit=7 is not setup for use with terrain'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' ... stopping inside init3d ... '
          if(dowr) write(outfile,*)
          call stopcm1
        ENDIF


      IF( tctype.eq.1 )THEN
        ! Rotunno-Emanuel (1987) vortex:
        dd2 = 2.0 * rmax_re87 / ( r0_re87 + rmax_re87 )
        do i=1,nref
          if(rref(i).lt.r0_re87)then
            dd1 = 2.0 * rmax_re87 / ( rref(i) + rmax_re87 )
            vt(i) = sqrt( vmax_re87**2 * (rref(i)/rmax_re87)**2     &
            * ( dd1 ** 3 - dd2 ** 3 ) + 0.25*fcor*fcor*rref(i)*rref(i) )   &
                    - 0.5 * fcor * rref(i)
          else
            vt(i) = 0.0
          endif
        enddo
      ELSEIF( tctype.eq.2 )THEN
        ! modified Rankine vortex:
        do i=1,nref
          if( rref(i).lt.rm1 )then
            vt(i) = vmax*rref(i)/rm1
          elseif( rref(i).lt.rm2 )then
            vt(i) = vmax*( (rref(i)/rm1)**rdc )
          elseif( rref(i).lt.rm3 )then
            v2 = vmax*( (rref(i)/rm1)**rdc )
            v3 = vmax*0.5*(1.0-(rref(i)-rm2)/(rm3-rm2))
            w3 = (rref(i)-rm2)/(rm3-rm2)
            w2 = 1.0-w3
            vt(i) = w2*v2+w3*v3
          else
            vt(i) = 0.0
          endif
        enddo
      ENDIF


        if( myid.eq.0 )then
          print *
          print *,'  tangential velocity profile:  r (km), v (m/s):'
          do i=1,nref
            print *,i,rref(i),vt(i)
          enddo
          print *
        endif


        ! just distribute linearly from surface to zdd
        ! (as in Rotunno-Emanuel 1987)
        do k=1,nk
        do i=1,nref
          if(zh(1,1,k).lt.zdd)then
            vref(i,k) = vt(i) * (zdd-zh(1,1,k))/(zdd-0.0)
          else
            vref(i,k) = 0.0
          endif
        enddo
        enddo

        do k=1,nk
          ! store environmental relative humidity (RH) in dum2 array:
          dum2(1,1,k) = qv0(1,1,k)/(rslf(prs0(1,1,k),th0(1,1,k)*pi0(1,1,k)))
          do i=1,nref
            rhref(i,k) = dum2(1,1,k)
            thref(i,k) = 0.0
            piref(i,k) = 0.0
            pref(i,k) = prs0(1,1,k)
          enddo
        enddo

        IF( rhmods )THEN
          ! increase RH in the inner core region (see parameters above)
          do k=1,nk
          do i=1,nref
            rhref(i,k) = max( dum2(1,1,k) , dum2(1,1,k) + (rhmax-dum2(1,1,k))*0.5*( 1.0+tanh( -(rref(i)-rmid)/rscale ) )  &
                                                                             *0.5*( 1.0+tanh( -(zh(1,1,k)-zmid)/zscale ) )  )
          enddo
          enddo
        ENDIF

      ! need to iterate for qv to converge:
      DO nloop=1,20

        ! get pressure (pref) and theta-v (thvref) from pi-pert (piref) and theta-pert (thref):
        do k=1,nk
        do i=1,nref
          pref(i,k) = p00*((pi0(1,1,k)+piref(i,k))**cpdrd)
          if(imoist.eq.1)then
            tnew=(th0(1,1,k)+thref(i,k))*(pi0(1,1,k)+piref(i,k))
            if( iice.eq.1 )then
              fliq=max(min((tnew-233.15)/(273.15-233.15),1.0),0.0)
              fice=1.0-fliq
              qsw=0.0
              if(tnew.gt.233.15)then
                qsw=fliq*rslf(pref(i,k),tnew)
              endif
              qsi=0.0
              if(tnew.lt.273.15)then
                qsi=fice*rsif(pref(i,k),tnew)
              endif
              qvs=qsw+qsi
            else
              qvs=rslf(pref(i,k),tnew)
            endif
            qvref(i,k) = rhref(i,k)*qvs
          endif
          thvref(i,k)=(th0(1,1,k)+thref(i,k))*(1.0+reps*qvref(i,k))   &
                                             /(1.0+qvref(i,k))
        enddo
        enddo

        ! get pi-pert from gradient-wind equation:
        do k=1,nk
          piref(nref,k)=0.0
          do i=nref,2,-1
            piref(i-1,k) = piref(i,k)                                       &
         + (rref(i-1)-rref(i))/(cp*0.5*(thvref(i-1,k)+thvref(i,k))) * 0.5 * &
             ( vref(i  ,k)*vref(i  ,k)/rref(i)                              &
              +vref(i-1,k)*vref(i-1,k)/rref(i-1)                            &
               + fcor * ( vref(i,k) + vref(i-1,k) ) )
          enddo
        enddo

        do i=1,nref
          piref(i,   0) = piref(i, 1)
          piref(i,nk+1) = piref(i,nk)
        enddo

        ! get theta-pert from hydrostatic equation:
        do k=2,nk
        do i=1,nref
          thref(i,k) = 0.5*( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))*rdz*mf(1,1,k+1)     &
                            +cp*0.5*(thvref(i,k)+thvref(i,k-1))*(piref(i,k)-piref(i,k-1))*rdz*mf(1,1,k) )   &
                          *thv0(1,1,k)/g
          thref(i,k)=(thv0(1,1,k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+reps*qvref(i,k))-th0(1,1,k)
        enddo
        enddo

        k=1
        do i=1,nref
          thref(i,k) = ( cp*0.5*(thvref(i,k)+thvref(i,k+1))*(piref(i,k+1)-piref(i,k))*rdz*mf(1,1,k+1) )   &
                          *thv0(1,1,k)/g
          thref(i,k)=(thv0(1,1,k)+thref(i,k))*(1.0+qvref(i,k))/(1.0+reps*qvref(i,k))-th0(1,1,k)
        enddo

        if(dowr) write(outfile,*) nloop,thref(1,1),qvref(1,1),piref(1,1)

      ENDDO   ! enddo for iteration


        tha = 0.0
        ppi = 0.0
        ua = 0.0
        va = 0.0

        xsave = xref
        ysave = yref

        grid:  &
        IF(axisymm.eq.1)THEN
          ! for axisymmetric model, we are done:

          do k=1,nk
          do j=0,nj+2
          do i=1,ni
             va(i,j,k) =  vref(i,k)
            ppi(i,j,k) = piref(i,k)
            tha(i,j,k) = thref(i,k)
            if(imoist.eq.1) qa(i,j,k,nqv) = qvref(i,k)
          enddo
          enddo
          enddo

        ELSE  grid
          ! for 3d model, we need to interpolate:

          if( myid.eq.0 ) print *
          if( myid.eq.0 ) print *,'  Interpolate to Cartesian grid: '

      ! loop through all possible neighbor domains:
      nnloop:  &
      do nn = 1 , 9

        doit = .false.

        if(     nn.eq.1 )then
          doit = .true.
          xref = xsave
          yref = ysave
        elseif( nn.eq.2 )then
          if( ebc.eq.1 .and. wbc.eq.1 .and. sbc.eq.1 .and. nbc.eq.1 )then
          doit = .true.
          xref = xsave + (maxx-minx)
          yref = ysave + (maxy-miny)
          endif
        elseif( nn.eq.3 )then
          if( ebc.eq.1 .and. wbc.eq.1 .and. sbc.eq.1 .and. nbc.eq.1 )then
          doit = .true.
          xref = xsave - (maxx-minx)
          yref = ysave - (maxy-miny)
          endif
        elseif( nn.eq.4 )then
          if( ebc.eq.1 .and. wbc.eq.1 .and. sbc.eq.1 .and. nbc.eq.1 )then
          doit = .true.
          xref = xsave + (maxx-minx)
          yref = ysave - (maxy-miny)
          endif
        elseif( nn.eq.5 )then
          if( ebc.eq.1 .and. wbc.eq.1 .and. sbc.eq.1 .and. nbc.eq.1 )then
          doit = .true.
          xref = xsave - (maxx-minx)
          yref = ysave + (maxy-miny)
          endif
        elseif( nn.eq.6 )then
          if( ebc.eq.1 .and. wbc.eq.1 )then
          doit = .true.
          xref = xsave + (maxx-minx)
          yref = ysave
          endif
        elseif( nn.eq.7 )then
          if( ebc.eq.1 .and. wbc.eq.1 )then
          doit = .true.
          xref = xsave - (maxx-minx)
          yref = ysave
          endif
        elseif( nn.eq.8 )then
          if( sbc.eq.1 .and. nbc.eq.1 )then
          doit = .true.
          xref = xsave
          yref = ysave + (maxy-miny)
          endif
        elseif( nn.eq.9 )then
          if( sbc.eq.1 .and. nbc.eq.1 )then
          doit = .true.
          xref = xsave
          yref = ysave - (maxy-miny)
          endif
        endif

        IF( doit )THEN

          ! interpolate to Cartesian grid:
          if( myid.eq.0 ) print *,'  nn,xref,yref = ',nn,xref,yref

          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            ! scalar points:
            tem1 = xh(i)-xref
            tem2 = yh(j)-yref
            rr = dsqrt( tem1**2 + tem2**2 )
          if( rr.lt.rrmax )then
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xref
                write(6,*) ' rref    = ',rref(ii-1),rref(ii-1)-rr
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            ppi(i,j,k) = piref(i1,k)+(piref(i2,k)-piref(i1,k))*frac
            tha(i,j,k) = thref(i1,k)+(thref(i2,k)-thref(i1,k))*frac
            if(imoist.eq.1) qa(i,j,k,nqv) = qvref(i1,k)+(qvref(i2,k)-qvref(i1,k))*frac
          else
            if(nn.eq.1)then
              tha(i,j,k) = thref(nref,k)
              ppi(i,j,k) = piref(nref,k)
              qa(i,j,k,nqv) = qvref(nref,k)
            endif
          endif
          enddo
          enddo
          enddo

          do k=1,nk
          do j=0,nj+1
          do i=0,ni+2
            ! u points:
            tem1 = xf(i)-xref
            tem2 = yh(j)-yref
            rr = dsqrt( tem1**2 + tem2**2 )
          if( rr.lt.tsmall )then
            ua(i,j,k) = 0.0
          elseif( rr.lt.rrmax )then
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xref
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            temd = -( dsin(datan2(tem2,tem1)) )
            ua(i,j,k) = ua(i,j,k) + ( vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*temd
          endif
          enddo
          enddo
          enddo

          do k=1,nk
          do j=0,nj+2
          do i=0,ni+1
            ! v points:
            tem1 = xh(i)-xref
            tem2 = yf(j)-yref
            rr = dsqrt( tem1**2 + tem2**2 )
          if( rr.lt.tsmall )then
            va(i,j,k) = 0.0
          elseif( rr.lt.rrmax )then
            ! need to account for grid stretching.  Do simple search:
            diff = -1.0e20
            ii = 0
            do while( diff.lt.0.0 )
              ii = ii + 1
              if( ii.gt.nref )then
                write(6,*)
                write(6,*) ' ii,nref = ',ii,nref
                write(6,*) ' rr      = ',rr,xref
                write(6,*)
                call stopcm1
              endif
              diff = rref(ii)-rr
            enddo
            if( abs(rr-rref(ii)).lt.tsmall .and. ii.eq.1 ) ii = 2
            i2 = ii
            i1 = i2-1
            frac = (      rr-rref(i1))   &
                  /(rref(i2)-rref(i1))
            temd =  ( dcos(datan2(tem2,tem1)) )
            va(i,j,k) = va(i,j,k) + ( vref(i1,k)+( vref(i2,k)- vref(i1,k))*frac )*temd
          endif
          enddo
          enddo
          enddo


!!!          print *
!!!          print *,'  symmtest:'
!!!          j = nj/2 + 5
!!!          k = 1
!!!          do i=1,nref
!!!            print *,i,j,ua(i,j,k),va(j,i,k),ua(i,j,k)+va(j,i,k)
!!!          enddo
!!!          print *

        ENDIF

      enddo  nnloop

          if( myid.eq.0 ) print *

          xref = xsave
          yref = ysave

          ! 180607
          ! final step: add base-state wind profile
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            ua(i,j,k) = ua(i,j,k)+u0(i,j,k)
            va(i,j,k) = va(i,j,k)+v0(i,j,k)
          enddo
          enddo
          enddo

        ENDIF  grid


        ! all done, set some boundary conditions, calculate final pressure, etc:

        call bcu(ua)
        call bcv(va)
        call bcs(ppi)
        call bcs(tha)

        call calcprs(pi0,prs,ppi)

        deallocate(  rref)
        deallocate(    vt)
        deallocate(  vref)
        deallocate( piref)
        deallocate(  pref)
        deallocate( thref)
        deallocate(thvref)
        deallocate( qvref)
        deallocate( rhref)

        setppi = .false.

      !-----------------------------------
      ! add random theta perts for 3d runs:
      ! (plus or minus this value in K)

      IF( nx.gt.3 .and. ny.gt.3 )THEN

        amplitude = 0.1

        do k=1,nk
          ! cm1r17: loop over entire domain
          do jj=1,ny
          do ii=1,nx
            call random_number(rand)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              ! only add perts in the warm core:
              if( tha(i,j,k).ge.0.1 )  &
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
        enddo

      ENDIF

!-----------------------------------------------------------------------
!  iinit = 8
!  Line thermal with random small-amplitude perturbations

      ELSEIF(iinit.eq.8)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Warm bubble'
        if(dowr) write(outfile,*)

        ric     =  centerx  ! center of bubble in x-direction (m)
        zc      =   1500.0  ! height of center of bubble above ground (m)
        bhrad   =  10000.0  ! horizontal radius of bubble (m)
        bvrad   =   1500.0  ! vertical radius of bubble (m)
        bptpert =      2.0  ! max potential temp perturbation (K)

        ! By default, CM1 sets qv=constant at a constant height level for 
        ! this value of iinit.  If you would rather have rh=constant at 
        ! a constant height level, then set this to .true.
        maintain_rh = .false.

        ! amplitude of random perturbations:
        amplitude = 0.20

        do k=1,nk
        do jj=1,ny
        do ii=1,nx
          call random_number(rand)
          i = ii - myi1 + 1
          j = jj - myj1 + 1
        IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
          beta=sqrt(                             &
                    ((xh(i)-ric)/bhrad)**2       &
                   +((zh(i,j,k)-zc)/bvrad)**2)
          if(beta.lt.1.0)then
            tha(i,j,k)=bptpert*(cos(0.5*pi*beta)**2)   &
                      +amplitude*(2.0*rand-1.0)
          else
            tha(i,j,k)=0.0
          endif
        ENDIF
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  iinit = 9
!  Forced convergence
!  Reference:  Loftus et al, 2008: MWR, v. 136, pp. 2408--2421.

      ELSEIF(iinit.eq.9)THEN

        ! User-defined settings:
        Dmax     =  -1.0e-3     ! maximum divergence (s^{-1})
        zdeep    =  2000.0      ! depth (m) of forced convergence
        lamx     = 10000.0      ! Loftus et al lambda_x parameter
        lamy     = 10000.0      ! Loftus at al lambda_y parameter
        xcent    =     0.0      ! x-location (m)
        ycent    =     0.0      ! y-location (m)
        convtime =   900.0      ! time (s) at beginning of simulation over
                                ! which convergence is applied

        ! Don't change anything below here:
        convinit = 1
        IF( ny.eq.1 )THEN
          ! 2D (x-z):
          Aconv = (-0.5*Dmax)/( (1.0/(lamx**2)) )
          lamy = 1.0e20
        ELSEIF( nx.eq.1 )THEN
          ! 2D (y-z):
          Aconv = (-0.5*Dmax)/( (1.0/(lamy**2)) )
          lamx = 1.0e20
        ELSE
          ! 3D:
          Aconv = (-0.5*Dmax)/( (1.0/(lamx**2))+(1.0/(lamy**2)) )
        ENDIF

!------------------------------------------------------------------
!  iinit = 10
!  momentum (u) forcing scheme (Morrison et al, 2015, JAS, pg 315)

      ELSEIF(iinit.eq.10)THEN

        xc_uforce     =  minx + 0.5*(maxx-minx)    ! x_c (m), center point of forcing in x
        xr_uforce     =  10000.0                   ! x_r (m), radius of forcing in x
        zr_uforce     =  10000.0                   ! z_r (m), radius of forcing in z
        alpha_uforce  =    0.1                     ! alpha (m/s/s), max intensity of forcing
        t1_uforce     =  3300.0                    ! time (s) to start ramping down u-forcing
        t2_uforce     =  3600.0                    ! time (s) to turn off u-forcing

!------------------------------------------------------------------
!  iinit = 11
!  Potential-temperature perturbation for inertia-gravity wave test case.
!  Reference:  Skamarock and Klemp, 1994, MWR, 122, 2623-2630.

      ELSEIF(iinit.eq.11)THEN

        do k=1,nk
        do j=1,nj
        do i=1,ni
          !----------
          ! Skamarock-Klemp-94 nonhydrostatic-scale inertia-gravity wave test:
          tha(i,j,k)=0.01*(sin(pi*zh(i,j,k)/10000.0))   &
                         /(1.0+((xh(i)-100000.0)/5000.0)**2)
          !----------
          ! Skamarock-Klemp-94 hydrostatic-scale inertia-gravity wave test:
!!!          tha(i,j,k)=0.01*(sin(pi*zh(i,j,k)/10000.0))   &
!!!                         /(1.0+((xh(i)-0.0)/100000.0)**2)
          !----------
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  iinit = 12
!  updraft nudging scheme (Naylor and Gilmore, 2012, MWR, pgs 3699-3705)

      ELSEIF(iinit.eq.12)THEN

        xc_wnudge     =  centerx                   ! x_c (m), center point of nudging in x
        xr_wnudge     =  10000.0                   ! x_r (m), radius of nudging in x

        yc_wnudge     =  centery                   ! y_c (m), center point of nudging in y
        yr_wnudge     =  10000.0                   ! y_r (m), radius of nudging in y

        zc_wnudge     =   1500.0                   ! z_c (m), center point of nudging in z
        zr_wnudge     =   1500.0                   ! z_r (m), radius of nudging in z

        alpha_wnudge  =      0.5                   ! alpha (1/s), inverse e-folding time for nudging
        wmax_wnudge   =     10.0                   ! w_max (m/s), max value of w to nudge towards

        t1_wnudge     =    900.0                   ! time (s) to start ramping down nudging
        t2_wnudge     =   1200.0                   ! time (s) to turn off nudging


        ! Don't change anything below here:
        wnudge = 1
        rxrwnudge = 1.0/max( xr_wnudge , 1.0e-12 )
        ryrwnudge = 1.0/max( yr_wnudge , 1.0e-12 )
        rzrwnudge = 1.0/max( zr_wnudge , 1.0e-12 )

!------------------------------------------------------------------
!  iinit = 13
!

      ELSEIF(iinit.eq.13)THEN

        bvrad = 1000.0

        do k=1,nk
        do j=1,nj
        do i=1,ni
          if( zh(i,j,k).lt.bvrad )then
            tha(i,j,k)=1.0*sin(2.0*pi*xh(i)/(0.25*(maxx-minx)))  &
                          *sin(2.0*pi*yh(j)/(0.25*(maxy-miny)))  &
                          *sin(2.0*pi*zh(i,j,k)/(2.0*bvrad))
          endif
        enddo
        enddo
        enddo

!------------------------------------------------------------------

      ENDIF    ! end of iinit options

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
      if(imoist.eq.1 .and. maintain_rh)then

        !! maintain rh
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Constant rh across domain:'
        if(dowr) write(outfile,*)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          dum2(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
          qa(i,j,k,nqv)=dum2(i,j,k)*rslf(prs0(i,j,k),(th0(i,j,k)+tha(i,j,k))*pi0(i,j,k))
        enddo
        enddo
        enddo

      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Random perturbations:

      IF( irandp.eq.1 )THEN

!        do n=1,nint(var7)
!          if(myid.eq.0) print *,'  pert loop = ',n
!        do k=1,nk
!          do jj = 0,ny+1
!          do ii = 0,nx+1
!            call random_number(rand)
!          enddo
!          enddo
!        enddo
!        enddo

      IF( testcase.eq.0 )THEN
        ! not a pre-configured test case: 

        ! this is the amplitude of the theta perturbations
        ! (plus or minus this value in K)
        amplitude = 0.25

        ! random numbers added here
        ! (can be modified to only place perturbations in certain
        !  locations, but this default code simply puts them
        !  everywhere)
        do k=1,nk
          ! cm1r17: loop over entire domain
          do jj = 0,ny+2
          do ii = 0,nx+2
            call random_number(rand)
            call random_number(rand2)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
        enddo

      ELSEIF( testcase.eq.1 .or. testcase.eq.2 .or. testcase.eq.6 .or. testcase.eq.10 .or. testcase.eq.11 .or. testcase.eq.15 )THEN

        ! 1: convective boundary layer (Sullivan and Patton, 2011, JAS)
        ! 2: shear-driven boundary layer (Moeng and Sullivan, 1994, JAS)
        ! 6: hurricane boundary layer (Bryan et al, 2017, BLM)

        amplitude = 0.10
        dpi = 4.0d0*datan(1.0d0)

        do k=1,nk
        if( zh(1,1,k).le.100.0 )then
          ! cm1r17: loop over entire domain
          do jj = 0,ny+1
          do ii = 0,nx+1
            call random_number(rand)
            call random_number(rand2)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
        endif
        enddo

      ELSEIF( testcase.eq.3 )THEN

        ! shallow Cu (Siebesma et al, 2003, JAS)

        do k=1,nk
        if( zh(1,1,k).le.1600.0 )then
          ! cm1r17: loop over entire domain
          do jj = 0,ny+1
          do ii = 0,nx+1
            call random_number(rand)
            call random_number(rand2)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+0.1*(2.0*rand-1.0)
              qa(i,j,k,nqv)=qa(i,j,k,nqv)+0.001*0.025*(2.0*rand2-1.0)
            ENDIF
          enddo
          enddo
        endif
        enddo

      ELSEIF( testcase.eq.4 .or. testcase.eq.5 .or. testcase.eq.7 )THEN

        amplitude = 0.10

        do k=1,nk
        if( zh(1,1,k).le.1000.0 )then
          ! cm1r17: loop over entire domain
          do jj = 0,ny+1
          do ii = 0,nx+1
            call random_number(rand)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
        endif
        enddo


      ELSEIF( testcase.eq.8 )THEN

        ! RCEMIP:
        do k=1,5
          amplitude = 0.10 - (k-1)*0.02
          if( myid.eq.0 ) print *,'  amplitude = ',k,amplitude
          ! cm1r17: loop over entire domain
          do jj = 0,ny+1
          do ii = 0,nx+1
            call random_number(rand)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
        enddo

      ELSEIF( testcase.eq.9 )THEN

        ! stable boundary layer (Beare et al. 2006)
        amplitude = 0.10
        do k=1,nk
          if( zh(1,1,k).lt.50.0 )then
          do jj = 0,ny+1
          do ii = 0,nx+1
            call random_number(rand)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
          endif
        enddo

      ELSEIF( testcase.eq.12 .or. testcase.eq.13 )THEN
        ! velocity perturbations:

        ! this is the amplitude of the theta perturbations
        ! (plus or minus this value in K)
        amplitude = 0.5

        ! random numbers added here
        ! (can be modified to only place perturbations in certain
        !  locations, but this default code simply puts them
        !  everywhere)
        do k=1,nk
          ! dum1,dum2 are define at scalar pts:
          do j=0,nj+1
          do i=0,ni+1
            dum1(i,j,k) = u0(1,1,k)
            dum2(i,j,k) = v0(1,1,k)
          enddo
          enddo
          ! cm1r17: loop over entire domain
          do jj = 0,ny+2
          do ii = 0+5,nx+1-5
            call random_number(rand)
            call random_number(rand2)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              dum1(i,j,k)=dum1(i,j,k)+amplitude*(2.0*rand-1.0)
              dum2(i,j,k)=dum2(i,j,k)+amplitude*(2.0*rand2-1.0)
!!!              tha(i,j,k)=tha(i,j,k)+0.5*(2.0*rand2-1.0)
            ENDIF
          enddo
          enddo
          ! interpolate to u,v pts:
          do j=1,nj+1
          do i=1,ni+1
            ua(i,j,k) = 0.5*(dum1(i-1,j,k)+dum1(i,j,k))
            va(i,j,k) = 0.5*(dum2(i,j-1,k)+dum2(i,j,k))
          enddo
          enddo
        enddo

      ELSE

        amplitude = 0.10

        do k=1,nk
          do jj = 0,ny+1
          do ii = 0,nx+1
            call random_number(rand)
            i = ii - myi1 + 1
            j = jj - myj1 + 1
            ! check to see if this processor has this gridpoint:
            IF( i.ge.ib .and. i.le.ie .and. j.ge.jb .and. j.le.je )THEN
              tha(i,j,k)=tha(i,j,k)+amplitude*(2.0*rand-1.0)
            ENDIF
          enddo
          enddo
        enddo

      ENDIF

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------------------
!  arrays for elliptic solver

      IF( (ibalance.eq.2).or.(psolver.eq.4).or.(psolver.eq.5).or.(pdcomp) )THEN

        dpi = 4.0d0*datan(1.0d0)
        if(dowr) write(outfile,*) '  dpi = ',dpi

        IF(psolver.le.3)THEN
          do k=1,nk
            cfa(k)=mh(1,1,k)*mf(1,1,k  )*rf0(1,1,k  )*0.5*(thv0(1,1,k-1)+thv0(1,1,k))/(dz*dz*rho0(1,1,k)*thv0(1,1,k))
            cfc(k)=mh(1,1,k)*mf(1,1,k+1)*rf0(1,1,k+1)*0.5*(thv0(1,1,k)+thv0(1,1,k+1))/(dz*dz*rho0(1,1,k)*thv0(1,1,k))
            ad1(k) = 1.0/(cp*rho0(1,1,k)*thv0(1,1,k))
            ad2(k) = 1.0
          enddo
          cfa( 1) = 0.0
          cfc(nk) = 0.0
          do j=jpb,jpe
          do i=ipb,ipe
            do k=1,nk
              cfb(i,j,k)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                                +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                                -2.0d0)/(dx*dx) - cfa(k) - cfc(k)
            enddo
          enddo
          enddo
        ELSE
          do k=1,nk
            cfa(k)=mh(1,1,k)*mf(1,1,k  )*rf0(1,1,k  )/(dz*dz*rho0(1,1,k-1))
            cfc(k)=mh(1,1,k)*mf(1,1,k+1)*rf0(1,1,k+1)/(dz*dz*rho0(1,1,k+1))
            ad1(k) = 1.0
            ad2(k) = 1.0/rho0(1,1,k)
          enddo
          cfa( 1) = 0.0
          cfc(nk) = 0.0
          do j=jpb,jpe
          do i=ipb,ipe
            do k=2,nk-1
              cfb(i,j,k)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                                +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                                -2.0d0)/(dx*dx)                               &
                    -mh(1,1,k)*mf(1,1,k+1)*rf0(1,1,k+1)/(dz*dz*rho0(1,1,k))   &
                    -mh(1,1,k)*mf(1,1,k  )*rf0(1,1,k  )/(dz*dz*rho0(1,1,k))
            enddo
            cfb(i,j,1)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                              +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                              -2.0d0)/(dx*dx)                               &
                  -mh(1,1,1)*mf(1,1,2  )*rf0(1,1,2  )/(dz*dz*rho0(1,1,1))
            cfb(i,j,nk)=2.0d0*( dcos(2.0d0*dpi*dble(i-1)/dble(ipe))          &
                              +dcos(2.0d0*dpi*dble(j-1)/dble(jpe))          &
                              -2.0d0)/(dx*dx)                               &
                  -mh(1,1,nk)*mf(1,1,nk  )*rf0(1,1,nk  )/(dz*dz*rho0(1,1,nk))
          enddo
          enddo
        ENDIF

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-----------------------------------------------------------------
!  Get 3d pressure


    IF(setppi)THEN

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        ppi(i,j,k)=0.0
      enddo
      enddo
      enddo

      IF(ibalance.eq.1)THEN

        ! hydrostatic balance ... integrate top-down

        do j=1,nj
        do i=1,ni
          ! virtual potential temperature

          if(imoist.eq.1)then
            do k=1,nk
              qt=0.0
              do n=nql1,nql2
                qt=qt+qa(i,j,k,n)
              enddo
              if(iice.eq.1)then
                do n=nqs1,nqs2
                  qt=qt+qa(i,j,k,n)
                enddo
              endif
              thvnew(k)=(th0(i,j,k)+tha(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                               /(1.0+qa(i,j,k,nqv)+qt)
            enddo
          else
            do k=1,nk
              thvnew(k)=th0(i,j,k)+tha(i,j,k)
            enddo
          endif

          ! non-dimensional pressure
          pinew(nk)=pi0(i,j,nk)
          do k=nk-1,1,-1
            pinew(k)=pinew(k+1)+g*(zh(i,j,k+1)-zh(i,j,k))   &
                    /(cp*0.5*(thvnew(k+1)+thvnew(k)))
          enddo

          ! new pressure
          do k=1,nk
            ppi(i,j,k)=pinew(k)-pi0(i,j,k)
            if(abs(ppi(i,j,k)).lt.1.0e-6) ppi(i,j,k)=0.0
          enddo

        enddo
        enddo

      ELSEIF(ibalance.eq.2)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  ibalance = 2'
        if(dowr) write(outfile,*)

        if(stretch_x.ge.1.or.stretch_y.ge.1)then
          print *,'  this option not supported with horizontal grid stretching'
          print *,'  (yet)'
          call stopcm1
        endif

#ifdef MPI
        print *,'  This option is not (yet) supported in MPI mode'
        print *,'  (sorry)'
        print *
        call stopcm1
#endif

        ! buoyancy pressure

        ! th3d stores theta-v

        if(imoist.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qt=0.0
            do n=nql1,nql2
              qt=qt+qa(i,j,k,n)
            enddo
            if(iice.eq.1)then
              do n=nqs1,nqs2
                qt=qt+qa(i,j,k,n)
              enddo
            endif
            th3d(i,j,k)=(th0(i,j,k)+tha(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                       /(1.0+qa(i,j,k,nqv)+qt)
          enddo
          enddo
          enddo
        else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            th3d(i,j,k)=th0(i,j,k)+tha(i,j,k)
          enddo
          enddo
          enddo
        endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          th3d(i,j,0   ) = th3d(i,j,1)
          th3d(i,j,nk+1) = th3d(i,j,nk)
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum4(i,j,k)=g*( th3d(i,j,k)/thv0(i,j,k)-1.0 )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum4(i,j,0   ) = -dum4(i,j,1)
          dum4(i,j,nk+1) = -dum4(i,j,nk)
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=0.5*( dum4(i,j,k-1)+dum4(i,j,k) )
        enddo
        enddo
        enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          ppi(i,j,k)=0.0
          dum3(i,j,k)=0.0
          divx(i,j,k)=0.0
          uten(i,j,k)=0.0
          vten(i,j,k)=0.0
        enddo
        enddo
        enddo

        call     poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                       dum1,dum2,dum3,dum4,divx,ppi,uten,vten,wten, &
                       cfb,cfa,cfc,ad1,ad2,pdt,lgbth,lgbph,rhs,trans,dtl)

        IF(psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=kb,ke
          do j=jb,je
          do i=ib,ie
            ppi(i,j,k)=((prs0(1,1,k)+ppi(i,j,k)*rho0(1,1,k))*rp00)**rovcp   &
                      -pi0(1,1,k)
            pp3d(i,j,k)=ppi(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

        call bcs(ppi)

      ENDIF

    ENDIF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  change base state:
!  WARNING:  for dycore testing only 
!    (do not modify)

  dobsmod = .false.
  bsmod:  IF( dobsmod )THEN
    print *,'  33333 '
    call stopcm1

        if(imoist.eq.1)then
          print *,' 77665 '
          call stopcm1
        endif

        print *
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
!-----------------------------------------------
! isentropic:
          th_sfc   =   288.0   ! Potential temperature of atmosphere (K)
          pi_sfc   =   1.0     ! Exner function at surface
          thfoo=th_sfc
          pifoo=pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          prsfoo=p00*(pifoo**cpdrd)
!-----------------------------------------------
! isothermal:
!          prs_sfc  = p00
!          t_sfc    = 333.15
!          hs       = rd*t_sfc/g        ! scale height of atmosphere
!          prsfoo=prs_sfc*EXP(-zh(i,j,k)/hs)
!          pifoo=(prsfoo/p00)**(rd/cp)
!          thfoo=t_sfc/pifoo
!-----------------------------------------------

          tha(i,j,k) = th0(i,j,k)+tha(i,j,k) - thfoo
          ppi(i,j,k) = pi0(i,j,k)+ppi(i,j,k) - pifoo

          th0(i,j,k) = thfoo
          pi0(i,j,k) = pifoo
          prs0(i,j,k) = prsfoo

          thv0(i,j,k) = th0(i,j,k)
          t0 = th0(i,j,k)*pi0(i,j,k)
          rho0(i,j,k) = prs0(i,j,k)/(rd*t0)

        enddo
        enddo
        enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rr0(i,j,k)=1.0/rho0(i,j,k)
      enddo
      enddo
      enddo

      do k=2,nk
      do j=jb,je
      do i=ib,ie
!!!        rf0(i,j,k)=0.5*(rho0(i,j,k-1)+rho0(i,j,k))
        rf0(i,j,k)=c1(i,j,k)*rho0(i,j,k-1)+c2(i,j,k)*rho0(i,j,k)
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,1) = cgs1*rho0(i,j,1)+cgs2*rho0(i,j,2)+cgs3*rho0(i,j,3)
        rf0(i,j,0)=rf0(i,j,1)
        rho0s(i,j) = rf0(i,j,1)
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,nk+1) = cgt1*rho0(i,j,nk)+cgt2*rho0(i,j,nk-1)+cgt3*rho0(i,j,nk-2)
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rrf0(i,j,k)=1.0/rf0(i,j,k)
      enddo
      enddo
      enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          rth0(i,j,k)=1.0/th0(i,j,k)
        enddo
        enddo
        enddo

      call bcs(ppi)

  ENDIF  bsmod

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Leaving INIT3D'

      end subroutine init3d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getset(restarted,mass1,cloudvar,dt,   &
                        xh,rxh,arh1,arh2,uh,ruh,xf,rxf,arf1,arf2,uf,ruf,yh,vh,rvh,yf,vf,rvf, &
                        zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,    &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,dumk1,dumk2, &
                        pi0,th0,rho0,rf0,prs0,thv0,ust,znt,u1,v1,s1,cm0,u0,v0,rth0,rr,rf,rho,prs, &
                        dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,                &
                        ua,u3d,va,v3d,wa,w3d,ppi,pp3d,ppx,phi1,phi2,            &
                        tha,th3d,qa,q3d,cme,csm,ce1,ce2,                        &
                        tkea,tke3d,kmh,kmv,khh,khv,nm,defv,defh,lenscl,dissten, &
                        t11,t12,t13,t22,t23,t33,iamsat,                         &
                        pta,pt3d,pdata,bndy,kbdy,timavg,                        &
                        upten,vpten,ucuten,vcuten,                              &
                        reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk,             &
                        nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                        &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                        &
                        uw31,uw32,ue31,ue32,us31,us32,un31,un32,                &
                        vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,                &
                        ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,                &
                        sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,                &
                        tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,                &
                        kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2)

      use input
      use constants
      use misclibs
      use bc_module
      use comm_module
      use parcel_module , only : getparcelzs
      use ib_module
      use turb_module , only : tkebc,calcnm,calcdef,turbsmag
#ifdef MPI
      use mpi
#endif
      implicit none
 
      logical, intent(in) :: restarted
      double precision, intent(inout) :: mass1
      logical, intent(in), dimension(maxq) :: cloudvar
      real, intent(inout) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      double precision, intent(inout), dimension(kb:ke) :: dumk1,dumk2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,rf0,prs0,thv0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th0
      real, intent(inout), dimension(ib:ie,jb:je) :: ust,znt,u1,v1,s1,cm0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rth0
      real, dimension(ib:ie,jb:je,kb:ke) :: rr,rf,rho,prs
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua,u3d
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va,v3d
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppx
      real, intent(inout), dimension(ibph:ieph,jbph:jeph,kbph:keph) :: phi1,phi2
      real, dimension(ib:ie,jb:je,kb:ke) :: tha,th3d
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: cme,csm,ce1,ce2
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: nm,defv,defh,lenscl,dissten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: iamsat
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      logical, intent(in), dimension(ibib:ieib,jbib:jeib,kbib:keib) :: bndy
      integer, intent(in), dimension(ibib:ieib,jbib:jeib) :: kbdy
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: upten,vpten
      real, intent(inout), dimension(ibcu:iecu,jbcu:jecu,kbcu:kecu) :: ucuten,vcuten
      real, intent(inout), dimension(ibta:ieta,jbta:jeta,kbta:keta,ntavr) :: timavg
      integer, dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, dimension(cmp,jmp,kmt)   :: tkw1,tkw2,tke1,tke2
      real, dimension(imp,cmp,kmt)   :: tks1,tks2,tkn1,tkn2
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2

!----------
 
      integer :: i,j,k,n
      real(kind=qp) :: temq
      real :: xw,xe,ys,yn,rtmp,rmin
      logical :: doit

      if(dowr) write(outfile,*) 'Inside GETSET'
      if(dowr) write(outfile,*)

!------------------------------------------------------------------

        IF( do_ib )THEN
          call zero_out_uv(bndy,kbdy,ua ,va )
          call zero_out_w(bndy,kbdy,wa )
          if( idoles .and. iusetke ) call zero_out_w(bndy,kbdy,tkea )
        ENDIF

!------------------------------------------------------------------
!  Make sure boundary values are set properly

      call bcu(ua)
      call bcv(va)
      call bcw(wa,1)
      call bcs(ppi)
      call bcs(tha)
      if(imoist.eq.1)then
        do n=1,numq
          call bcs(qa(ibm,jbm,kbm,n))
        enddo
      endif
      if( idoles .and. iusetke )then
        call bcw(tkea,1)
      endif
      if(iptra.eq.1)then
        do n=1,npt
          call bcs(pta(ib,jb,kb,n))
        enddo
      endif

      call bcu(u0)
      call bcv(v0)
      call bcs(th0)

#ifdef MPI
!------------------------------------------------------------------

      nf=0
      nu=0
      nv=0
      nw=0

      call comm_3u_start(ua,uw31,uw32,ue31,ue32,   &
                            us31,us32,un31,un32,reqs_u)
      call comm_3u_end(ua,uw31,uw32,ue31,ue32,   &
                          us31,us32,un31,un32,reqs_u)

      call comm_3v_start(va,vw31,vw32,ve31,ve32,   &
                            vs31,vs32,vn31,vn32,reqs_v)
      call comm_3v_end(va,vw31,vw32,ve31,ve32,   &
                          vs31,vs32,vn31,vn32,reqs_v)

      call comm_3w_start(wa,ww31,ww32,we31,we32,   &
                            ws31,ws32,wn31,wn32,reqs_w)
      call comm_3w_end(wa,ww31,ww32,we31,we32,   &
                          ws31,ws32,wn31,wn32,reqs_w)

      call comm_3s_start(ppi,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(ppi,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)

      call comm_3s_start(tha,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(tha,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)

      IF(imoist.eq.1)THEN
        do n=1,numq
          call comm_3s_start(qa(ibm,jbm,kbm,n),sw31,sw32,se31,se32,   &
                                               ss31,ss32,sn31,sn32,reqs_s)
          call comm_3s_end(qa(ibm,jbm,kbm,n),sw31,sw32,se31,se32,   &
                                             ss31,ss32,sn31,sn32,reqs_s)
        enddo
      ENDIF

      IF( idoles .and. iusetke )THEN
        call comm_3t_start(tkea,tkw1,tkw2,tke1,tke2,   &
                                tks1,tks2,tkn1,tkn2,reqs_tk)
        call comm_3t_end(tkea,tkw1,tkw2,tke1,tke2,   &
                              tks1,tks2,tkn1,tkn2,reqs_tk)
      ENDIF

      IF(iptra.eq.1)THEN
        do n=1,npt
          call comm_3s_start(pta(ib,jb,kb,n),sw31,sw32,se31,se32,   &
                                             ss31,ss32,sn31,sn32,reqs_s)
          call comm_3s_end(pta(ib,jb,kb,n),sw31,sw32,se31,se32,   &
                                           ss31,ss32,sn31,sn32,reqs_s)
        enddo
      ENDIF

      call comm_3u_start(u0,uw31,uw32,ue31,ue32,   &
                            us31,us32,un31,un32,reqs_u)
      call comm_3u_end(u0,uw31,uw32,ue31,ue32,   &
                          us31,us32,un31,un32,reqs_u)

      call comm_3v_start(v0,vw31,vw32,ve31,ve32,   &
                            vs31,vs32,vn31,vn32,reqs_v)
      call comm_3v_end(v0,vw31,vw32,ve31,ve32,   &
                          vs31,vs32,vn31,vn32,reqs_v)

      call comm_3s_start(th0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(th0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

#endif
      if(terrain_flag)then
        call bcwsfc(gz,dzdx,dzdy,ua,va,wa)
        call bc2d(wa(ib,jb,1))
      endif

      if( idoles .and. iusetke )then
        do j=1,nj
        do i=1,ni
          tkea(i,j,1) = tkea(i,j,2)
        enddo
        enddo
        call bc2d(tkea(ibt,jbt,1))
#ifdef MPI
        call comm_1s2d_start(tkea(ibt,jbt,1),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),   &
                                             ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
        call comm_1s2d_end(tkea(ibt,jbt,1),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),   &
                                           ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
        call bcs2_2d(tkea(ibt,jbt,1))
        call comm_2d_corner(tkea(ibt,jbt,1))
#endif
      endif

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          rth0(i,j,k)=1.0/th0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------

      do k=kb,ke
      do j=jb,je
      do i=ib,ie+1
        u3d(i,j,k)=ua(i,j,k)
      enddo
      enddo
      enddo
 
      do k=kb,ke
      do j=jb,je+1
      do i=ib,ie
        v3d(i,j,k)=va(i,j,k)
      enddo
      enddo
      enddo
 
      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        w3d(i,j,k)=wa(i,j,k)
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        pp3d(i,j,k)=ppi(i,j,k)
        th3d(i,j,k)=tha(i,j,k)
      enddo
      enddo
      enddo

      if(imoist.eq.1)then
        do n=1,numq
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          q3d(i,j,k,n)=qa(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      endif

      if( idoles .and. iusetke )then
        do k=kbt,ket
        do j=jbt,jet
        do i=ibt,iet
          tke3d(i,j,k)=tkea(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(iptra.eq.1)then
        do n=1,npt
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          pt3d(i,j,k,n)=pta(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      endif

      if( psolver.eq.6 )then
        call bcs(phi1)
#ifdef MPI
        call comm_1s_start(phi1,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(phi1,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
#endif
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          phi2(i,j,k)=phi1(i,j,k)
        enddo
        enddo
        enddo
      endif

!------------------------------------------------------------------
!  Get stuff

  IF( ( .not. restarted ) .or.( interp_on_restart ) )THEN

    IF(psolver.eq.4.or.psolver.eq.5.or.psolver.eq.6.or.psolver.eq.7)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        rho(i,j,k)=rho0(i,j,k)
        prs(i,j,k)=prs0(i,j,k)
      enddo
      enddo
      enddo

    ELSE

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  getting prs,rho '
      if(myid.eq.0) print *

      IF( imoist.eq.1 )THEN

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          rho(i,j,k)=prs(i,j,k)                         &
             /( (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))     &
               *(rd+max(0.0,q3d(i,j,k,nqv))*rv) )
        enddo
        enddo
        enddo

      ELSE

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

  ENDIF


    IF( psolver.le.3 )THEN

        call bcs(rho)
#ifdef MPI
        call comm_1s_start(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(rho,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call bcs2(rho)
        call getcorner(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
#endif

      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          rr(i,j,k) = 1.0/rho(i,j,k)
          rf(i,j,k) = (c1(i,j,k)*rho(i,j,k-1)+c2(i,j,k)*rho(i,j,k))
        enddo
        enddo
        enddo

        call bcs(rr)
#ifdef MPI
        call comm_1s_start(rr,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(rr,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call bcs2(rr)
        call getcorner(rr,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
#endif
        call bcs(rf)
#ifdef MPI
        call comm_1s_start(rf,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(rf,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call bcs2(rf)
        call getcorner(rf,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
#endif

        ! meh 1 !
      !$omp parallel do default(shared)  &
      !$omp private(i,j,k)
        do j=0,nj+1
        do i=0,ni+1
          ! cm1r17, 2nd-order extrapolation:
          rf(i,j,1) = cgs1*rho(i,j,1)+cgs2*rho(i,j,2)+cgs3*rho(i,j,3)
          rf(i,j,nk+1) = cgt1*rho(i,j,nk)+cgt2*rho(i,j,nk-1)+cgt3*rho(i,j,nk-2)
        enddo
        enddo

    ELSE

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rho(i,j,k) = rho0(i,j,k)
        rr(i,j,k) = 1.0/rho0(i,j,k)
        rf(i,j,k) = rf0(i,j,k)
      enddo
      enddo
      enddo

    ENDIF

!------------------------------------------------------------------
!  cm1r18:  get total mass of dry air at t=0

      IF( ( .not. restarted ) .or. ( interp_on_restart ) )THEN

        if(myid.eq.0) print *
        if(myid.eq.0) print *,'  Getting mass1 '
        if(myid.eq.0) print *

        dumk1 = 0.0

        IF( axisymm.eq.0 )THEN
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dumk1(k) = dumk1(k) + rho(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)
          enddo
          enddo
          enddo
        ELSEIF( axisymm.eq.1 )THEN
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dumk1(k) = dumk1(k) + rho(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)
          enddo
          enddo
          enddo
        ELSE
          stop 2223
        ENDIF

#ifdef MPI
        call MPI_ALLREDUCE(mpi_in_place,dumk1(1),nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

        temq = 0.0

        do k=1,nk
          temq = temq + dumk1(k)
        enddo

        mass1 = temq

        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  mass1 = ',mass1
        if( myid.eq.0 ) print *

      ENDIF

!------------------------------------------------------------------
!  cm1r19.5:

    IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.6 .or. psolver.eq.7 )THEN

      IF( ( .not. restarted ) .or. ( interp_on_restart ) )THEN

        if( myid.eq.0 ) print *
        if( myid.eq.0 ) print *,'  setting initial ppx '
        if( myid.eq.0 ) print *

        IF( psolver.eq.2 .or. psolver.eq.3 .or. psolver.eq.7 )THEN

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppx(i,j,k)=ppi(i,j,k)
          enddo
          enddo
          enddo

        ELSEIF( psolver.eq.6 )THEN

          !$omp parallel do default(shared)   &
          !$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppx(i,j,k)=phi1(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

      ENDIF

      call bcs(ppx)
#ifdef MPI
      call comm_1s_start(ppx,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      call comm_1s_end(ppx,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
#endif

    ENDIF

!------------------------------------------------------------------

      if( iprcl.eq.1 .and. terrain_flag )THEN
        ! get terrain at parcel locations:
        if( prsig.ge.1 )then
          call getparcelzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,pdata)
        else
          print *,'  invalid value for prsig: ',prsig
          call stopcm1
        endif
        if( .not. restarted )then
          ! 181022:  not a restart ... initialize sigma
          do n=1,nparcels
            ! get sigma from z:
            ! (see Section 3 of "The governing equations for CM1", 
            !  http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf)
            pdata(n,prsig) = zt*(pdata(n,prz)-pdata(n,przs))/(zt-pdata(n,przs))
          enddo
        else
          ! 181022:  this is a restart ... get z
          do n=1,nparcels
            ! get z from sigma:
            ! (see Section 3 of "The governing equations for CM1", 
            !  http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf)
            pdata(n,prz) = pdata(n,przs) + pdata(n,prsig)*((zt-pdata(n,przs))*rzt)
          enddo
        endif
      endif

!------------------------------------------------------------------

      ifidoles:  &
      IF( idoles )THEN

        cm1setupval:  &
        if( cm1setup.eq.1 )then

          ! LES:
          ! cm0 is the same everywhere

          cm0 = c_m

        elseif( cm1setup.eq.4 )then  cm1setupval

          ! LES embedded within mesoscale model:
          ! cm0 is non-zero in the inner part of the domain only

          cm0 = 0.0

          IF( les_subdomain_shape.eq.1 )THEN
            ! square/rectangle:

            xw = centerx - 0.5*les_subdomain_xlen
            xe = centerx + 0.5*les_subdomain_xlen
            ys = centery - 0.5*les_subdomain_ylen
            yn = centery + 0.5*les_subdomain_ylen

            if( myid.eq.0 ) print *
            if( myid.eq.0 ) print *,'  centerx,centery = ',centerx,centery
            if( myid.eq.0 ) print *,'  xw,xe           = ',xw,xe
            if( myid.eq.0 ) print *,'  ys,yn           = ',ys,yn
            if( myid.eq.0 ) print *

            do j=jb,je
            do i=ib,ie
              if( xh(i).gt.xw .and. xh(i).lt.xe .and. yh(j).gt.ys .and. yh(j).lt.yn )then
                cm0(i,j) = c_m
                ! 220301:
                if( xh(i).lt. xw+les_subdomain_trnslen ) cm0(i,j) = min( cm0(i,j) ,  c_m*(xh(i)-xw)/les_subdomain_trnslen )
                if( xh(i).gt. xe-les_subdomain_trnslen ) cm0(i,j) = min( cm0(i,j) , -c_m*(xh(i)-xe)/les_subdomain_trnslen )
                if( yh(j).lt. ys+les_subdomain_trnslen ) cm0(i,j) = min( cm0(i,j) ,  c_m*(yh(j)-ys)/les_subdomain_trnslen )
                if( yh(j).gt. yn-les_subdomain_trnslen ) cm0(i,j) = min( cm0(i,j) , -c_m*(yh(j)-yn)/les_subdomain_trnslen )
              endif
            enddo
            enddo

          ELSE

            if( myid.eq.0 )then
              print *
              print *,'  undefined value for les_subdomain_shape '
              print *
              print *,'    24987 '
              print *
            endif

            call stopcm1

          ENDIF

        endif  cm1setupval

        ! cm0 should now be set ... get other variables:
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          if( cm0(i,j).gt.cmemin )then
            cme(i,j,k) = cm0(i,j)
            ce1(i,j,k) = cme(i,j,k) * c_l * c_l * ( 1.0 / ri_c - 1.0 )
            ce2(i,j,k) = max( 0.0 , cme(i,j,k) * pi * pi - ce1(i,j,k) )
            csm(i,j,k) = ( cme(i,j,k) * cme(i,j,k) * cme(i,j,k) / ( ce1(i,j,k) + ce2(i,j,k) ) )**0.25   ! Smagorinsky constant
          endif
        enddo
        enddo
        enddo

      ENDIF  ifidoles

!------------------------------------------------------------------
!  Get eddy diffusivities if interpolating on restart:

    doit = .true.
    IF( doit )THEN

    iorst:  &
    IF( interp_on_restart .and. ( cm1setup.eq.1 .or. cm1setup.eq.4 ) )THEN
    if( idoles .and. iusetke )then

      if(myid.eq.0) print *
      if(myid.eq.0) print *,'  Getting initial subgrid tke ... '
      if(myid.eq.0) print *

        iamsat = .false.
        call calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,nm,dum1,dum2,dum3,dum4,dum5,dum6,   &
                    prs,ppi,tha,qa,iamsat)

        call calcdef(    rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                     xh,rxh,arh1,arh2,uh,xf,rxf,arf1,arf2,uf,vh,vf,mh,c1,c2,mf,defv,defh,  &
                     dum1,dum2,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy,rho,rr,rf)


        call     turbsmag(  0  ,dt,ruh,rvh,rmh,mf,rmf,th0,rf,              &
                          nm,defv,defh,dum4,dum5  ,dum6  ,zf,znt  ,ust,csm, &
                          kmh,kmv,khh,khv,lenscl,dissten,                  &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                          kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                          ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                          kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                          ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))

      do k=2,nk
      do j=1,nj
      do i=1,ni
        if( cm0(i,j).gt.cmemin )  &
        tkea(i,j,k) = (kmv(i,j,k)/(cme(i,j,k)*lenscl(i,j,k)))**2
      enddo
      enddo
      enddo

        IF( do_ib ) call zero_out_w(bndy,kbdy,tkea )
        call bcw(tkea,1)
#ifdef MPI
        call comm_3t_start(tkea,tkw1,tkw2,tke1,tke2,   &
                                tks1,tks2,tkn1,tkn2,reqs_tk)
        call comm_3t_end(tkea,tkw1,tkw2,tke1,tke2,   &
                              tks1,tks2,tkn1,tkn2,reqs_tk)
#endif
        do j=1,nj
        do i=1,ni
          tkea(i,j,1) = tkea(i,j,2)
        enddo
        enddo
        call bc2d(tkea(ibt,jbt,1))
#ifdef MPI
        call comm_1s2d_start(tkea(ibt,jbt,1),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),   &
                                             ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
        call comm_1s2d_end(tkea(ibt,jbt,1),sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),   &
                                           ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1),reqs_s)
        call bcs2_2d(tkea(ibt,jbt,1))
        call comm_2d_corner(tkea(ibt,jbt,1))
#endif
        do k=kbt,ket
        do j=jbt,jet
        do i=ibt,iet
          tke3d(i,j,k)=tkea(i,j,k)
        enddo
        enddo
        enddo

    endif
    ENDIF  iorst
    ENDIF

!------------------------------------------------------------------

      if( use_pbl )then

        call bcs(upten)
        call bcs(vpten)
#ifdef MPI
        call comm_1s_start(upten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(upten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_start(vpten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
        call comm_1s_end(vpten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
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

      endif

!------------------------------------------------------------------

    IF( cuparam.eq.1 )THEN
      call bcs(ucuten)
      call bcs(vcuten)
#ifdef MPI
      call comm_1s_start(ucuten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      call comm_1s_end(ucuten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      call comm_1s_start(vcuten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
      call comm_1s_end(vcuten,pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p)
#endif
      ENDIF

!------------------------------------------------------------------

        tavgcode:  &
        IF( dotimeavg .and. ( (.not. restarted) .or. interp_on_restart ) )THEN
          if(myid.eq.0) print *
          if(myid.eq.0) print *,'  initializing timeavg fields '
          if(myid.eq.0) print *
          if(utav.ge.1)then
            do k=1,keta
            do j=0,nj+1
            do i=0,ni+2
              timavg(i,j,k,utav) = u3d(i,j,k)
            enddo
            enddo
            enddo
          endif
          if(vtav.ge.1)then
            do k=1,keta
            do j=0,nj+2
            do i=0,ni+1
              timavg(i,j,k,vtav) = v3d(i,j,k)
            enddo
            enddo
            enddo
          endif
          if(wtav.ge.1)then
            do k=1,keta
            do j=0,nj+1
            do i=0,ni+1
              timavg(i,j,k,wtav) = w3d(i,j,k)
            enddo
            enddo
            enddo
          endif
          if(ttav.ge.1)then
            do k=1,keta
            do j=0,nj+1
            do i=0,ni+1
              timavg(i,j,k,ttav) = th3d(i,j,k)
            enddo
            enddo
            enddo
          endif
          if(etav.ge.1)then
            do k=1,keta
            do j=0,nj+1
            do i=0,ni+1
              timavg(i,j,k,etav) = tke3d(i,j,k)
            enddo
            enddo
            enddo
          endif
          if(qtav.ge.1)then
            do k=1,keta
            do j=0,nj+1
            do i=0,ni+1
              timavg(i,j,k,qtav) = q3d(i,j,k,nqv)
            enddo
            enddo
            enddo
          endif
        ENDIF  tavgcode

!------------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Leaving GETSET'
 
      end subroutine getset


  END MODULE init3d_module
