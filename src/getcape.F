  MODULE getcape_module

  implicit none

  private
  public :: getcape

  CONTAINS

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine getcape( source, nk , p_in , t_in , q_in, cape , cin,   &
                        zlcl, zlfc, zel , psource , tsource , qvsource )
    implicit none

    integer, intent(in) :: source,nk
    real, dimension(nk), intent(in) :: p_in,t_in,q_in
    real, intent(out) :: cape,cin,psource,tsource,qvsource

!-----------------------------------------------------------------------
!
!  getcape - a fortran90 subroutine to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!   *** Modified version to calculate Lifted Condensation Level (LCL)  ***
!   *** Level of Free Convection (LFC), and Equilibrium Level (EL)     ***
!
!  Version 1.04                           Last modified:  8 October 2010
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
!               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
!
!-----------------------------------------------------------------------
!
!  Input:     nk - number of levels in the sounding (integer)
!
!           p_in - one-dimensional array of pressure (mb) (real)
!
!           t_in - one-dimensional array of temperature (C) (real)
!
!           q_in - one-dimensional array of water vapor mixing ratio (kg/kg) (real)
!
!  Output:  cape - Convective Available Potential Energy (J/kg) (real)
!
!            cin - Convective Inhibition (J/kg) (real)
!
!-----------------------------------------------------------------------
!  User options:

    real, parameter :: pinc = 200.0   ! Pressure increment (Pa)
                                      ! (smaller number yields more accurate
                                      !  results,larger number makes code 
                                      !  go faster)

!!!    integer, parameter :: source = 2    ! Source parcel:
!!!                                        ! 1 = surface
!!!                                        ! 2 = most unstable (max theta-e)
!!!                                        ! 3 = mixed-layer (specify ml_depth)

    real, parameter :: ml_depth =  500.0  ! depth (m) of mixed layer 
                                          ! for source=3

    integer, parameter :: adiabat = 1   ! Formulation of moist adiabat:
                                        ! 1 = pseudoadiabatic, liquid only
                                        ! 2 = reversible, liquid only
                                        ! 3 = pseudoadiabatic, with ice
                                        ! 4 = reversible, with ice

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    logical :: doit,ice,cloud,not_converged
    integer :: k,kmax,n,nloop,i,orec
    real, dimension(nk) :: p,t,td,pi,q,th,thv,z,pt,pb,pc,pn,ptv,ptd,pqv,pql

    real :: the,maxthe,parea,narea,lfc
    real :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
    real :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2
    real :: thlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
    double precision :: avgth,avgqv
    real :: ee,zlcl,zlfc,zel,plcl,plfc,pel,zh

!-----------------------------------------------------------------------

    real, parameter :: g     = 9.81
    real, parameter :: p00   = 100000.0
    real, parameter :: cp    = 1005.7
    real, parameter :: rd    = 287.04
    real, parameter :: rv    = 461.5
    real, parameter :: xlv   = 2501000.0
    real, parameter :: xls   = 2836017.0
    real, parameter :: t0    = 273.15
    real, parameter :: cpv   = 1875.0
    real, parameter :: cpl   = 4190.0
    real, parameter :: cpi   = 2118.636
    real, parameter :: lv1   = xlv+(cpl-cpv)*t0
    real, parameter :: lv2   = cpl-cpv
    real, parameter :: ls1   = xls+(cpi-cpv)*t0
    real, parameter :: ls2   = cpi-cpv

    real, parameter :: rp00  = 1.0/p00
    real, parameter :: eps   = rd/rv
    real, parameter :: reps  = rv/rd
    real, parameter :: rddcp = rd/cp
    real, parameter :: cpdrd = cp/rd
    real, parameter :: cpdg  = cp/g

    real, parameter :: converge = 0.0002

    integer, parameter :: debug_level =   0

!-----------------------------------------------------------------------

!---- convert p,t,td to mks units; get pi,q,th,thv ----!

    do k=1,nk
        p(k) = 100.0*p_in(k)
        t(k) = 273.15+t_in(k)
        q(k) = max( q_in(k) , 1.0e-10 )
!!!       td(k) = 273.15+td_in(k)
       ee = alog((q(k)/eps)*p(k)/100.0/(1.0+(q(k)/eps)))
       td(k) = 273.15+(243.5*ee-440.8)/(19.48-ee)
       pi(k) = (p(k)*rp00)**rddcp
!!!        q(k) = getqvl(p(k),td(k))
       th(k) = t(k)/pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
!!!      print *,k,t(k)-273.15,td(k)-273.15
    enddo

!---- get height using the hydrostatic equation ----!

    z(1) = 0.0
    do k=2,nk
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
      z(k) = z(k-1) + dz
    enddo

!----------------------------------------------------------------

!---- find source parcel ----!

  IF(source.eq.1)THEN
    ! use surface parcel
    kmax = 1

  ELSEIF(source.eq.2)THEN
    ! use most unstable parcel (max theta-e)

    IF(p(1).lt.50000.0)THEN
      ! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthx(p(1),t(1),td(1),q(1))
    ELSE
      ! find max thetae below 500 mb
      maxthe = 0.0
      do k=1,nk
        if(p(k).ge.50000.0)then
          the = getthx(p(k),t(k),td(k),q(k))
          if( the.gt.maxthe )then
            maxthe = the
            kmax = k
          endif
        endif
      enddo
    ENDIF
    if(debug_level.ge.100) print *,'  kmax,maxthe = ',kmax,maxthe

  ELSEIF(source.eq.3)THEN
    ! use mixed layer

    IF( (z(2)-z(1)).gt.ml_depth )THEN
      ! the second level is above the mixed-layer depth:  just use the
      ! lowest level

      avgth = th(1)
      avgqv = q(1)
      kmax = 1

    ELSEIF( z(nk).lt.ml_depth )THEN
      ! the top-most level is within the mixed layer:  just use the
      ! upper-most level

      avgth = th(nk)
      avgqv = q(nk)
      kmax = nk

    ELSE
      ! calculate the mixed-layer properties:

      avgth = 0.0
      avgqv = 0.0
      k = 2
      if(debug_level.ge.100) print *,'  ml_depth = ',ml_depth
      if(debug_level.ge.100) print *,'  k,z,th,q:'
      if(debug_level.ge.100) print *,1,z(1),th(1),q(1)

      do while( (z(k).le.ml_depth) .and. (k.le.nk) )

        if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

        avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
        avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))

        k = k + 1

      enddo

      th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
      qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))

      if(debug_level.ge.100) print *,999,ml_depth,th2,qv2

      avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
      avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))

      if(debug_level.ge.100) print *,k,z(k),th(k),q(k)

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      kmax = 1

    ENDIF

    if(debug_level.ge.100) print *,avgth,avgqv

  ELSE

    print *
    print *,'  Unknown value for source'
    print *
    print *,'  source = ',source
    print *
    stop

  ENDIF

!---- define parcel properties at initial location ----!
    narea = 0.0

  if( (source.eq.1).or.(source.eq.2) )then
    k    = kmax
    th2  = th(kmax)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = t(kmax)
    thv2 = thv(kmax)
    qv2  = q(kmax)
    b2   = 0.0
  elseif( source.eq.3 )then
    k    = kmax
    th2  = avgth
    qv2  = avgqv
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = th2*pi2
    b2   = g*( thv2-thv(kmax) )/thv(kmax)
  endif

    psource = p2
    tsource = t2
   qvsource = qv2

    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

    doit = .true.
    cloud = .false.
    if(adiabat.eq.1.or.adiabat.eq.2)then
      ice = .false.
    else
      ice = .true.
    endif

      the = getthx(p2,t2,t2,qv2)
      if(debug_level.ge.100) print *,'  the = ',the

      pt(k) = t2
      if( cloud )then
        ptd(k) = t2
      else
        ptd(k) = gettd(p2,t2,qv2)
      endif
      ptv(k) = t2*(1.0+reps*qv2)/(1.0+qv2)
      pb(k) = 0.0
      pqv(k) = qv2
      pql(k) = 0.0

      zlcl = -1.0
      zlfc = -1.0
      zel  = -1.0

!---- begin ascent of parcel ----!

      if(debug_level.ge.100)then
        print *,'  Start loop:'
        print *,'  p2,th2,qv2 = ',p2,th2,qv2
      endif

    do while( doit .and. (k.lt.nk) )

        k = k+1
       b1 =  b2

       dp = p(k-1)-p(k)

      if( dp.lt.pinc )then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/float(nloop)
      endif

      do n=1,nloop

         p1 =  p2
         t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2

        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp

        thlast = th1
        i = 0
        not_converged = .true.

        do while( not_converged )
          i = i + 1
          t2 = thlast*pi2
          if(ice)then
            fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min( qt , fliq*getqvl(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max( fice*(qt-qv2) , 0.0 )
          ql2 = max( qt-qv2-qi2 , 0.0 )

          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm=rd+rv*qvbar
          cpm=cp+cpv*qvbar+cpl*qlbar+cpi*qibar
          th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar)     &
                       +lhs*(qi2-qi1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*alog(p2/p1) )

          if(i.gt.90) print *,i,th2,thlast,th2-thlast
          if(i.gt.100)then
            print *
            print *,'  Error:  lack of convergence'
            print *
            print *,'  ... stopping iteration '
            print *
!!!            stop 1001
            ! 171020:
            cape = 0.0
            cin = 0.0
            psource = 0.0
            tsource = 0.0
            qvsource = 0.0
            return
          endif
          if( abs(th2-thlast).gt.converge )then
            thlast=thlast+0.3*(th2-thlast)
          else
            not_converged = .false.
          endif
        enddo

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if( ql2.ge.1.0e-10 ) cloud = .true.
        if( cloud .and. zlcl.lt.0.0 )then
           zlcl = z(k-1)+(z(k)-z(k-1))*float(n)/float(nloop)
           plcl = p(k-1)+(p(k)-p(k-1))*float(n)/float(nloop)
        endif

        IF(adiabat.eq.1.or.adiabat.eq.3)THEN
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        ELSEIF(adiabat.le.0.or.adiabat.ge.5)THEN
          print *
          print *,'  Undefined adiabat'
          print *
          stop 10000
        ENDIF

      enddo

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
        b2 = g*( thv2-thv(k) )/thv(k)
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      if( zlcl.gt.0.0 .and. zlfc.lt.0.0 .and. b2.gt.0.0 )then
        if( b1.gt.0.0 )then
          zlfc = zlcl
          plfc = plcl
        else
          zlfc = z(k-1)+(z(k)-z(k-1))*(0.0-b1)/(b2-b1)
          plfc = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        endif
      endif

      if( zlfc.gt.0.0 .and. zel.lt.0.0 .and. b2.lt.0.0 )then
        zel = z(k-1)+(z(k)-z(k-1))*(0.0-b1)/(b2-b1)
        pel = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
      endif

      the = getthx(p2,t2,t2,qv2)

      pt(k) = t2
      if( cloud )then
        ptd(k) = t2
      else
        ptd(k) = gettd(p2,t2,qv2)
      endif
      ptv(k) = t2*(1.0+reps*qv2)/(1.0+qv2)
      pb(k) = b2
      pqv(k) = qv2
      pql(k) = ql2

      ! Get contributions to CAPE and CIN:

      if( (b2.ge.0.0) .and. (b1.lt.0.0) )then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea =  0.5*b2*dz*frac
        narea = narea-0.5*b1*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
        cin  = cin  + narea
        narea = 0.0
      elseif( (b2.lt.0.0) .and. (b1.gt.0.0) )then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea =  0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
        if(debug_level.ge.200)then
          print *,'      b1,b2 = ',b1,b2
          print *,'      p1,ps,p2 = ',p(k-1),ps,p(k)
          print *,'      frac = ',frac
          print *,'      parea = ',parea
          print *,'      narea = ',narea
        endif
      elseif( b2.lt.0.0 )then
        ! still collecting negative buoyancy
        parea =  0.0
        narea = narea-0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2)
        narea =  0.0
      endif

      cape = cape + max(0.0,parea)
      pc(k) = cape

      if(debug_level.ge.200)then
        write(6,102) p2,b1,b2,cape,cin,cloud
102     format(5(f13.4),2x,l1)
      endif

      if( (p(k).le.10000.0).and.(b2.lt.0.0) )then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      endif

    enddo

!!!    print *,'  zlcl,zlfc,zel = ',zlcl,zlfc,zel
!!!    print *,'  plcl,plfc,pel = ',plcl,plfc,pel

!---- All done ----!

    end subroutine getcape

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getqvl(p,t)
    implicit none

    real :: p,t,es

    real, parameter :: eps = 287.04/461.5

    es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
    ! 171023 (fix for very cold temps):
    es = min( es , p*0.5 )
    getqvl = eps*es/(p-es)

    end function getqvl

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getqvi(p,t)
    implicit none

    real :: p,t,es

    real, parameter :: eps = 287.04/461.5

    es = 611.2*exp(21.8745584*(t-273.15)/(t-7.66))
    ! 171023 (fix for very cold temps):
    es = min( es , p*0.5 )
    getqvi = eps*es/(p-es)

    end function getqvi

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function getthx(p,t,td,q)
    implicit none

    real :: p,t,td,q
    real :: tlcl

    if( (td-t).ge.-0.1 )then
      tlcl = t
    else
      tlcl = 56.0 + ( (td-56.0)**(-1) + 0.00125*alog(t/td) )**(-1)
    endif

    getthx=t*( (100000.0/p)**(0.2854*(1.0-0.28*q)) )   &
            *exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )

    end function getthx

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real function gettd(p,t,q)
    implicit none

    real :: p,t,q

    real :: el
    real, parameter :: eps = 287.04/461.5

    el = alog((q/eps)*p/100.0/(1.0+(q/eps)))
    gettd = 273.15+(243.5*el-440.8)/(19.48-el)

    end function gettd

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

  END MODULE getcape_module
