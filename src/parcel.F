#define TROUBLE 0
#define PITER 1048576
#ifdef _B4B
#undef _B4B01F
#undef _B4B02F
#undef _B4B03F
#endif
  MODULE parcel_module

  implicit none

  private
  public :: parcel_driver,parcel_interp,parcel_write,setup_parcel_vars,getparcelzs,tri_interp,get2d

  ! If we already know that a parcel stays on the same process 
  !    in the previous time step and it is still here, we start 
  !    the location search from the previous index and limit 
  !    the search range using the parameter below;
  ! This assumes that a parcel does not move too far and does
  !    not move between processes within this subroutine?
  integer, parameter :: location_offset = 3 

  

  CONTAINS

      subroutine parcel_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,zs,    &
                               sigma,sigmaf,znt,rho,ua,va,wa,pdata)
      use input, only : ib,ie,jb,je,kb,ke,numq,ni,nj,nk,imp,jmp,kmp,rmp,cmp, &
          kmt,npvals,nparcelsLocal,npvars,cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,bbc,tbc,imove,zt,rzt, &
          umove,vmove,nqv,terrain_flag,nx,ny,axisymm,nodex,nodey,myi,myj,viscosity, &
          prx,pry,prz,prsig,pract,prvpx,prvpy,prvpz,prtp,prrp,prmult,prms,przs, &
          dx,dy,njp1,nkp1,maxz,minx,maxx,miny,maxy,timestats,mytime,ierr,time_parcels,myid
      use constants
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine updates the parcel locations
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: znt
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d
      integer :: nrkp
      real :: dt2,uu1,vv1,ww1
      real :: z0,rznt,var
      real :: sigdot,sig1,zsp,sig3d

      logical, parameter :: debug = .false.

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      ! no slip:
      if( imove.eq.1 )then
        IF(j.le.je)THEN
        do i=ib,ie+1
          ua(i,j,0) = 0.0 - umove
        enddo
        ENDIF
        do i=ib,ie
          va(i,j,0) = 0.0 - vmove
        enddo
      else
        IF(j.le.je)THEN
        do i=ib,ie+1
          ua(i,j,0) = 0.0
        enddo
        ENDIF
        do i=ib,ie
          va(i,j,0) = 0.0
        enddo
      endif
    ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = 0.0
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      IF(j.le.je)THEN
      do i=ib,ie
        wa(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, update it's location:

    dt2 = dt/2.0

    !$acc parallel loop gang vector default(present) 
    nploop:  &
    DO np=1,nparcelsLocal

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)
      if( .not. terrain_flag )then
        z3d = pdata(np,prz)
      else
        sig3d = pdata(np,prsig)
      endif

      iflag = undefined_index
      jflag = undefined_index
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit1:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit1

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                     (jflag.ge.1.and.jflag.le.nj) )THEN

      rkloop:  DO nrkp = 1,2

      IF( nrkp.eq.1 )THEN
        i=iflag
        j=jflag
      ELSE
        iflag = undefined_index
        jflag = undefined_index
        IF(nx.eq.1)THEN
          iflag = 1
        ELSE
          ! cm1r19:
          i = ni+2
          do while( iflag.lt.0 .and. i.gt.0 )
            i = i-1
            if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
              iflag = i
            endif
          enddo
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          jflag = 1
        ELSE
          do j=0,nj+1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) ) jflag=j
          enddo
          ! cm1r19:
          j = nj+2
          do while( jflag.lt.0 .and. j.gt.0 )
            j = j-1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
              jflag = j
            endif
          enddo
        ENDIF
        i=iflag
        j=jflag
      ENDIF

        IF(debug)THEN
        if( i.lt.0 .or. i.gt.(ni+1) .or. j.lt.0 .or. j.gt.(nj+1) )then
          print *,'  myid,i,j = ',myid,i,j
          print *,'  x,x1     = ',x3d,pdata(np,prx)
          print *,'  y,y1     = ',y3d,pdata(np,pry)
          do i=0,ni+1
            print *,i,abs(xh(i)-x3d),0.5*dx*ruh(i)
          enddo
          do j=0,nj+1
            print *,j,abs(yh(j)-y3d),0.5*dy*rvh(j)
          enddo
          call stopcm1
        endif
        ENDIF

        kflag = 1
        if( .not. terrain_flag )then
          do while( z3d.gt.zf(iflag,jflag,kflag+1) )
            kflag = kflag+1
          enddo
        else
          do while( sig3d.gt.sigmaf(kflag+1) )
            kflag = kflag+1
          enddo
        endif

        IF(debug)THEN
        if( kflag.le.0 .or. kflag.ge.(nk+1) )then
          print *,myid,nrkp
          print *,iflag,jflag,kflag
          print *,pdata(np,prx),pdata(np,pry),pdata(np,prz)
          print *,x3d,y3d,z3d
          print *,uval,vval,wval
          print *,zf(iflag,jflag,kflag),z3d,zf(iflag,jflag,kflag+1)
          print *,'  16667 '
          call stopcm1
        endif
        ENDIF

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        if( .not. terrain_flag )then
          if( z3d.lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
          rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        else
          if( sig3d.lt.sigma(k) )then
            k=k-1
          endif
          rz = ( sig3d-sigma(k) )/( sigma(k+1)-sigma(k) )
        endif

        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )

        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz
        ! call calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),x3d,xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        uval = tri_interp(ie+1,je,ke,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( .not. terrain_flag )then
          if( z3d.lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( sig3d.lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        ! w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        ! w2=rx*(1.0-ry)*(1.0-rz)
        ! w3=(1.0-rx)*ry*(1.0-rz)
        ! w4=(1.0-rx)*(1.0-ry)*rz
        ! w5=rx*(1.0-ry)*rz
        ! w6=(1.0-rx)*ry*rz
        ! w7=rx*ry*(1.0-rz)
        ! w8=rx*ry*rz
        call calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),y3d,yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        vval = tri_interp(ni,njp1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va)

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        rx = rxs
        ry = rys
        if( .not. terrain_flag )then
          rz = ( z3d-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )
        else
          rz = ( sig3d-sigmaf(k) )/( sigmaf(k+1)-sigmaf(k) )
        endif

        ! saveit:
        rzw = rz

        ! w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        ! w2=rx*(1.0-ry)*(1.0-rz)
        ! w3=(1.0-rx)*ry*(1.0-rz)
        ! w4=(1.0-rx)*(1.0-ry)*rz
        ! w5=rx*(1.0-ry)*rz
        ! w6=(1.0-rx)*ry*rz
        ! w7=rx*ry*(1.0-rz)
        ! w8=rx*ry*rz
        call calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),z3d,zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        wval = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa)

        if( terrain_flag )then
          sigdot =  tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa)
          zsp = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,zs)
          z3d = zsp + sig3d*((zt-zsp)*rzt)
        endif

!----------------------------------------------------------------------
!  uv for parcels below lowest model level:

      IF( bbc.eq.3 )THEN
        ! semi-slip lower boundary condition:
        if( z3d.lt.zh(1,1,1) )then
          ! re-calculate velocities if parcel is below lowest model level:
          !------
          ! u at lowest model level:
          i=iflag
          j=jflag
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          uval = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,ua(ib,jb,1))
          !------
          ! v at lowest model level:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          vval = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,va(ib,jb,1))
          !------
          ! z0:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          z0 = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt)
          !------
          ! get u,v from (neutral) log-layer equation:
          rznt = 1.0/z0
          var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
          if( imove.eq.1 )then
            uval = (uval+umove)*var - umove
            vval = (vval+vmove)*var - vmove
          else
            uval = uval*var
            vval = vval*var
          endif
        endif
      ENDIF

!-----------------------------------------------------
!  Update parcel positions:
!-----------------------------------------------------

      ! RK2 scheme:
      IF(nrkp.eq.1)THEN
        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pdata(np,prx)+dt*uval
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pdata(np,pry)+dt*vval
        ENDIF
        if( terrain_flag )then
          sig3d = pdata(np,prsig) + dt*sigdot
          sig1 = sigdot
        else
          z3d = pdata(np,prz)+dt*wval
          ww1=wval
        endif
        uu1=uval
        vv1=vval
      ELSE
        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pdata(np,prx)+dt2*(uu1+uval)
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pdata(np,pry)+dt2*(vv1+vval)
        ENDIF
        if( terrain_flag )then
          sig3d = pdata(np,prsig) + dt2*(sig1+sigdot)
          IF( sig3d.lt.0.0 )THEN
            print *,'  parcel is below surface:  np,x3d,y3d,sig3d = ',np,x3d,y3d,sig3d
            sig3d=1.0e-6
          ENDIF
          sig3d=min(sig3d,maxz)
        else
          z3d = pdata(np,prz)+dt2*(ww1+wval)
          IF( z3d.lt.0.0 )THEN
            print *,'  parcel is below surface:  np,x3d,y3d,z3d = ',np,x3d,y3d,z3d
            z3d=1.0e-6
          ENDIF
          z3d=min(z3d,maxz)
        endif
      ENDIF


      ENDDO  rkloop

!-----------------------------------------------------
!  Account for boundary conditions (if necessary)
!-----------------------------------------------------

        ! New for cm1r17:  if parcel exits domain,
        ! just assume periodic lateral boundary conditions
        ! (no matter what actual settings are for wbc,ebc,sbc,nbc)

        if(x3d.lt.minx)then
          x3d=x3d+(maxx-minx)
        endif
        if(x3d.gt.maxx)then
          x3d=x3d-(maxx-minx)
        endif

        if( (y3d.gt.maxy).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d-(maxy-miny)
        endif
        if( (y3d.lt.miny).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d+(maxy-miny)
        endif

        pdata(np,prx)=x3d
        pdata(np,pry)=y3d
        if( .not. terrain_flag )then
          pdata(np,prz)=z3d
        else
          pdata(np,prsig)=sig3d
        endif

#ifdef MPI

      ELSE

        ! set to really small number (so we can use the allreduce command below)
        pdata(np,prx) = -1.0e30
        pdata(np,pry) = -1.0e30
        if( .not. terrain_flag )then
          pdata(np,prz) = -1.0e30
        else
          pdata(np,prsig) = -1.0e30
        endif
#endif

      ENDIF  myparcel

    ENDDO  nploop
    if(timestats.ge.1) time_parcels=time_parcels+mytime()

!----------------------------------------------------------------------
!  communicate data  (for MPI runs)

#ifdef MPI
      !$acc update host(pdata)
      if( .not. terrain_flag )then
        call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,1),npvars*nparcelsLocal,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      else
        call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,1),(npvars-1)*nparcelsLocal,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,prsig),nparcelsLocal,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      endif
      !$acc update device(pdata)
#endif

!----------------------------------------------------------------------
!  get height ASL:

      if( terrain_flag )then
            call getparcelzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,pdata)
            DO np=1,nparcelsLocal
              ! get z from sigma:
              ! (see Section 3 of "The governing equations for CM1", 
              !  http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf)
              pdata(np,prz) = pdata(np,przs) + pdata(np,prsig)*((zt-pdata(np,przs))*rzt)
            ENDDO
      endif

!----------------------------------------------------------------------

      end subroutine parcel_driver


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_interp(dt,mtime,xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                               zh,mh,rmh,zf,mf,znt,ust,c1,c2,          &
                               zs,sigma,sigmaf,rds,gz,                 &
                               pi0,th0,thv0,qv0,qc0,qi0,rth0,          &
                               dum1,dum2,dum3,dum4,zv  ,qt  ,prs,rho,  &
                               dum7,dum8,buoy,vpg  ,                   &
                               u3d,v3d,w3d,pp3d,th   ,t     ,th3d,q3d, &
                               kmh,kmv,khh,khv,tke3d,pt3d,pdata,       &
                               tdiag,qdiag,                            &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,reqs_p, &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)
      use input, only : ib,ie,jb,je,kb,ke,itb,ite,jtb,jte,ibm,iem,jbm,jem, &
          kbm,kem,numq,ibc,iec,jbc,jec,kbc,kec,ibt,iet,jbt,jet,kbt,ket,npt,ierr, &
          ibdt,iedt,jbdt,jedt,kbdt,kedt,ntdiag,ibdq,iedq,jbdq,jedq,kbdq,kedq, &
          nqdiag,jmp,imp,kmp,rmp,cmp,kmt,ibp,iep,jbp,jep,kbp,kep,nparcelsLocal,npvals,npvars, &
          ni,nj,nk,nqv,nql1,nql2,nqs1,nqs2,nqv,imoist,timestats,mytime,time_phys_d2h, &
          imoist,iice,time_parceli,bbc,tbc,rdz,cgs1,cgs2,cgs3,cgt1,cgt2,cgt3,imove, &
          qd_dbz,rdx,rdy,nx,ny,myi,myj,nodex,nodey,myid,nip1,njp1,nkp1,axisymm,nnc1, &
          umove,vmove,prth,prt,prqsl,prqsi,prvpg,prb,prprs,prrho,prpt1,prqv,prq1,prnc1, &
          prkm,prkh,prtke,prdbz,przv,prx,pry,prz,prtime,prq2,prnc2,prznt,prust,pru,prv, &
          prw,prsig,time_parceli_reduce,time_lb,terrain_flag
      use constants
      use cm1libs , only : rslf,rsif
      use bc_module, only: bcu2, bcv2, bcw2
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine interpolates model information to the parcel locations
!  (diagnostic only ... not used for model integration)
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      double precision, intent(in) :: mtime
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: znt,ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(kb:ke) :: rds
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,thv0,qv0,qc0,qi0,rth0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,zv,qt,prs,rho
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum7,dum8
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d,buoy,vpg
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pp3d,th3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th,t
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: q3d
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pt3d
      real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(cmp,jmp,kmt) :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt) :: tks1,tks2,tkn1,tkn2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tem,tem1
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d,z0,rznt,var

      logical, parameter :: debug = .false.

      if(timestats.ge.1) time_phys_D2H=time_phys_D2H+mytime()

!----------------------------------------------------------------------
!  Get derived variables:

    IF(imoist.eq.1)THEN
      ! with moisture:

    !$omp parallel do default(shared) private(i,j,k)
    !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
    do k=1,nk
      do j=1,nj
      do i=1,ni
        qt(i,j,k)=q3d(i,j,k,nqv)
      enddo
      enddo
    enddo
    if(nql1.ge.1)then
      do n=nql1,nql2
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          qt(i,j,k)=qt(i,j,k)+q3d(i,j,k,n)
        enddo
        enddo
      enddo
      enddo
    endif
      IF(iice.eq.1)THEN
        do n=nqs1,nqs2
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          qt(i,j,k)=qt(i,j,k)+q3d(i,j,k,n)
        enddo
        enddo
        enddo
    enddo
      ENDIF
      IF( prth.ge.1 .or. prt.ge.1 .or. prqsl.ge.1 .or. prqsi.ge.1 .or.  prvpg.ge.1 )THEN
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          th(i,j,k) = (th0(i,j,k)+th3d(i,j,k))
          t(i,j,k) = th(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
        enddo
      ENDIF
      IF( prb.ge.1 .or. prvpg.ge.1 )THEN
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k)             &
                           +repsm1*(q3d(i,j,k,nqv)-qv0(i,j,k))  &
                           -(qt(i,j,k)-q3d(i,j,k,nqv)-qc0(i,j,k)-qi0(i,j,k))   )
        enddo
        enddo
        enddo
      ENDIF
      IF( prvpg.ge.1 )THEN
        !$omp parallel do default(shared) private(i,j,k)
        !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum8(i,j,k) = th(i,j,k)*(1.0+reps*q3d(i,j,k,nqv))/(1.0+qt(i,j,k))
        enddo
        enddo
        enddo
      ENDIF


    ELSE
      ! dry:


      IF( prth.ge.1 .or. prt.ge.1 .or. prvpg.ge.1 )THEN
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
       do k=1,nk
        do j=1,nj
        do i=1,ni
          th(i,j,k)= (th0(i,j,k)+th3d(i,j,k))
          t(i,j,k) = th(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
       enddo
      ENDIF
      IF( prb.ge.1 .or. prvpg.ge.1 )THEN
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
       do k=1,nk
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k) )
        enddo
        enddo
       enddo
      ENDIF
      IF( prvpg.ge.1 )THEN
       !$omp parallel do default(shared) private(i,j,k)
       !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
       do k=1,nk
        do j=1,nj
        do i=1,ni
          dum8(i,j,k) = th(i,j,k)
        enddo
        enddo
       enddo
      ENDIF


    ENDIF


    IF( prb.ge.1 .or. prvpg.ge.1 )THEN
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        buoy(i,j,k) = (c1(1,1,k)*dum7(i,j,k-1)+c2(1,1,k)*dum7(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parceli=time_parceli+mytime()
      !!$acc update device(buoy)
      call    prepcornert(buoy,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,1)
      !!$acc update host(buoy)

      !$omp parallel do default(shared) private(i,j)
      !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        buoy(i,j,1) = buoy(i,j,2)+(buoy(i,j,3)-buoy(i,j,2))  &
                                 *(  zf(i,j,1)-  zf(i,j,2))  &
                                 /(  zf(i,j,3)-  zf(i,j,2))
        buoy(i,j,nk+1) = buoy(i,j,nk)+(buoy(i,j,nk  )-buoy(i,j,nk-1))  &
                                     *(  zf(i,j,nk+1)-  zf(i,j,nk  ))  &
                                     /(  zf(i,j,nk  )-  zf(i,j,nk-1))
      enddo
      enddo
    ENDIF
    IF( prvpg.ge.1 )THEN
    if( .not. terrain_flag )then
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k,tem1)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tem1 = rdz*cp*mf(1,1,k)
        vpg(i,j,k) = -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))  &
                          *(c2(1,1,k)*dum8(i,j,k)+c1(1,1,k)*dum8(i,j,k-1))
      enddo
      enddo
      enddo
    else
      !$omp parallel do default(shared) private(i,j,k)
      !$acc parallel loop gang vector collapse(3) default(present) private(i,j,k,tem1)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        tem1 = rds(k)*cp
        vpg(i,j,k) = -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))*gz(i,j)  &
                          *(c2(1,1,k)*dum8(i,j,k)+c1(1,1,k)*dum8(i,j,k-1))
      enddo
      enddo
      enddo
    endif
      if(timestats.ge.1) time_parceli=time_parceli+mytime()
      !!$acc update device(vpg)
      call    prepcornert(vpg ,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,1)
      !!$acc update host(vpg)
      ! cmr18:  at top/bottom boundaries, vpg + buoy = 0
      !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
      do j=0,nj+1
      do i=0,ni+1
        vpg(i,j,1) = -buoy(i,j,1)
        vpg(i,j,nk+1) = -buoy(i,j,nk+1)
      enddo
      enddo
    ENDIF

    if(timestats.ge.1) time_parceli=time_parceli+mytime()

!----------------------------------------------------------------------
!  get corner info for MPI runs
!  (may not parallelize correctly if this is not done)

#ifdef MPI
      call sync()
      if(timestats.ge.1) time_lb=time_lb+mytime()
      call getcorneru(u3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call getcornerv(v3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call getcornerw(w3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcu2(u3d)
      call bcv2(v3d)
      call bcw2(w3d)
#endif

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v,s the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

  IF(bbc.eq.1)THEN
    !$omp parallel do default(shared) private(i,j)
    !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
    DO j=jb,je
      ! free slip ... extrapolate:
      do i=ib,ie+1
        u3d(i,j,0) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
      enddo
    ENDDO
    !$omp parallel do default(shared) private(i,j)
    !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
    DO j=jb,je+1
      do i=ib,ie
        v3d(i,j,0) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
      enddo
    ENDDO
  ELSEIF(bbc.eq.2)THEN
    ! no slip:
    if( imove.eq.1 )then
      !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
      DO j=jb,je
        do i=ib,ie+1
          u3d(i,j,0) = 0.0 - umove
        enddo
      ENDDO
      !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
      DO j=jb,je+1
        do i=ib,ie
          v3d(i,j,0) = 0.0 - vmove
        enddo
       ENDDO
    else
      !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
      DO j=jb,je
        do i=ib,ie+1
          u3d(i,j,0) = 0.0
        enddo
      ENDDO
      !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
      DO j=jb,je+1
        do i=ib,ie
          v3d(i,j,0) = 0.0
        enddo
       ENDDO
    endif
  ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations
  ENDIF

!----------

  IF(tbc.eq.1)THEN
    !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
    DO j=jb,je
      ! free slip ... extrapolate:
      do i=ib,ie+1
        u3d(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
      enddo
    ENDDO
    !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
    DO j=jb,je+1
      do i=ib,ie
        v3d(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
      enddo
    ENDDO
  ELSEIF(tbc.eq.2)THEN
    !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
    DO j=jb,je
      ! no slip:
      do i=ib,ie+1
        u3d(i,j,nk+1) = 0.0
      enddo
    ENDDO
    !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
    DO j=jb,je+1
      do i=ib,ie
        v3d(i,j,nk+1) = 0.0
      enddo
    ENDDO
  ENDIF

!----------

  !$acc parallel loop gang vector collapse(2) default(present) private(i,j)
  DO j=jb,je
    do i=ib,ie
      w3d(i,j,nk+1) = 0.0
    enddo
  ENDDO

      if(timestats.ge.1) time_parceli=time_parceli+mytime()

      if( prth.ge.1 )then
        call prepcorners(th ,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prt.ge.1 )then
        call prepcorners(t  ,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prprs.ge.1 )then
        call prepcorners(prs,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prrho.ge.1 )then
        call prepcorners(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if(prpt1.ge.1)then
        do n=1,npt
          call prepcorners(pt3d(ib,jb,kb,n),nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        enddo
      endif
      if( prqv.ge.1 )then
        call prepcorners(q3d(ib,jb,kb,nqv),nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
      endif
      if( prq1.ge.1 .or. prnc1.ge.1 )then
        do n = 1,numq
          call prepcorners(q3d(ib,jb,kb,n),nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        enddo
      endif
      !!$acc update host(q3d)
      !print *,'parcel_interp: point #8'
      !!$acc update device(kmh,kmv)
      if( prkm.ge.1 )then
        call prepcornert(kmh,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
        call prepcornert(kmv,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
      endif
      !!$acc update host(kmh,kmv)
      !print *,'parcel_interp: point #9'
      !!$acc update device(khh,khv)
      if( prkh.ge.1 )then
        call prepcornert(khh,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
        call prepcornert(khv,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
      endif
      !!$acc update host(khh,khv)
      !print *,'parcel_interp: point #10'
      !!$acc update device(tke3d)
      if( prtke.ge.1 )then
        call prepcornert(tke3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
      endif
      !!$acc update host(tke3d)
      !print *,'parcel_interp: point #11'
      !!$acc update device(qdiag)
      if( prdbz.ge.1 )then
        call prepcorners(qdiag(ibdq,jbdq,kbdq,qd_dbz),  &
                             nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      !!$acc update host(qdiag)

!----------------------------------------------------------------------

      !print *,'parcel_interp: point #12'
    IF( prqsl.ge.1 )THEN
      !$acc parallel loop gang vector collapse(3) default(present) 
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = rslf( prs(i,j,k) , t(i,j,k) )
      enddo
      enddo
      enddo
      !$acc end parallel

      if(timestats.ge.1) time_parceli=time_parceli+mytime()

      call prepcorners(dum1,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    ENDIF
    IF( prqsi.ge.1 )THEN
      !$acc parallel loop gang vector collapse(3) default(present) 
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum2(i,j,k) = rsif( prs(i,j,k) , t(i,j,k) )
      enddo
      enddo
      enddo
      !$acc end parallel

      if(timestats.ge.1) time_parceli=time_parceli+mytime()

      call prepcorners(dum2,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    ENDIF
    !  print *,'parcel_interp: point #13'

!----------------------------------------------------------------------
!  Get zvort at appropriate C-grid location:
!  (assuming no terrain)
!  cm1r18:  below lowest model level:
!           Use extrapolated velocities for bbc=1,2
!           Use log-layer equations for bbc=3 (see below)

    IF( przv.ge.1)THEN
      !$acc parallel loop gang vector collapse(3) default(present)
      do k=0,nk+1
      do j=1,nj+1
      do i=1,ni+1
        zv(i,j,k) = (v3d(i,j,k)-v3d(i-1,j,k))*rdx*uf(i)   &
                   -(u3d(i,j,k)-u3d(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo
      !$acc end parallel
    ENDIF
    !HERE print *,'parcel_interp: before nploop2'
!----------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, get interpolated info:
    !PORTME
    !$acc parallel loop gang vector default(present)
    nploop2:  &
    DO np=1,nparcelsLocal

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)
      z3d = pdata(np,prz)

      iflag = undefined_index
      jflag = undefined_index
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit2:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit2

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myprcl:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                   (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        kflag = 1
        if( .not. terrain_flag )then
          do while( pdata(np,prz).ge.zf(iflag,jflag,kflag+1) )
            kflag = kflag+1
          enddo
        else
          do while( pdata(np,prsig).ge.sigmaf(kflag+1) )
            kflag = kflag+1
          enddo
        endif

        x3d = pdata(np,prx)
        y3d = pdata(np,pry)
        z3d = pdata(np,prz)

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,pry).lt.yh(j) )then
          j=j-1
        endif
        if( .not. terrain_flag )then
          if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
          rz = ( pdata(np,prz)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        else
          if( pdata(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
          rz = ( pdata(np,prsig)-sigma(k) )/( sigma(k+1)-sigma(k) )
        endif

        rx = ( pdata(np,prx)-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( pdata(np,pry)-yh(j) )/( yh(j+1)-yh(j) )


        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),pdata(np,prx),xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pdata(np,pry),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pdata(np,prz),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        uval =  tri_interp(nip1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,u3d)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( .not. terrain_flag )then
          if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( pdata(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = ( pdata(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( pdata(np,pry)-yf(j) )/( yf(j+1)-yf(j) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pdata(np,prx),xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),pdata(np,pry),yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pdata(np,prz),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        vval = tri_interp(ni,njp1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,v3d)

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(np,pry).lt.yh(j) )then
          j=j-1
        endif

!!!        rx = ( pdata(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
!!!        ry = ( pdata(np,pry)-yh(j) )/( yh(j+1)-yh(j) )
        rx = rxs
        ry = rys
        if( .not. terrain_flag )then
          rz = ( pdata(np,prz)-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )
        else
          rz = ( pdata(np,prsig)-sigmaf(k) )/( sigmaf(k+1)-sigmaf(k) )
        endif

        ! saveit:
        rzw = rz


        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz
        ! call calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pdata(np,prx),xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pdata(np,pry),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),pdata(np,prz),zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        wval = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,w3d)
      if(prkm.ge.1)then
        pdata(np,prkm  ) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kmh)
        pdata(np,prkm+1) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kmv)
      endif
      if(prkh.ge.1)then
        pdata(np,prkh  ) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,khh)
        pdata(np,prkh+1) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,khv)
      endif
      if( prtke.ge.1 )then
        pdata(np,prtke) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,tke3d)
      endif
      if( prb.ge.1 )then
        pdata(np,prb) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,buoy)
      endif
      if( prvpg.ge.1 )then
        pdata(np,prvpg) = tri_interp(ni,nj,nkp1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,vpg)
      endif

!----------------------------------------------------------------------
!  Data on scalar points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(np,pry).lt.yh(j) )then
          j=j-1
        endif
        if( .not. terrain_flag )then
          if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( pdata(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = rxs
        ry = rys
        rz = rzs

        ! w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        ! w2=rx*(1.0-ry)*(1.0-rz)
        ! w3=(1.0-rx)*ry*(1.0-rz)
        ! w4=(1.0-rx)*(1.0-ry)*rz
        ! w5=rx*(1.0-ry)*rz
        ! w6=(1.0-rx)*ry*rz
        ! w7=rx*ry*(1.0-rz)
        ! w8=rx*ry*rz
        call calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15558: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

      if(imoist.eq.1)then
        if(prdbz.ge.1)  &
        pdata(np,prdbz) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qdiag(ibdq,jbdq,kbdq,qd_dbz))
        if(prqv.ge.1)  &
        pdata(np,prqv) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,q3d(ib,jb,kb,nqv))
        if( prq1.ge.1 .and. nql1.ge.1 )then
          do n=nql1,nql1+(prq2-prq1)
            pdata(np,prq1+(n-nql1)) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,q3d(ib,jb,kb,n))
          enddo
        endif
        if(prnc1.ge.1)then
          do n=nnc1,nnc1+(prnc2-prnc1)
            pdata(np,prnc1+(n-nnc1)) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,q3d(ib,jb,kb,n))
          enddo
        endif
        if( prqsl.ge.1 )  &
        pdata(np,prqsl) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dum1)
        if( prqsi.ge.1 )  &
        pdata(np,prqsi) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,dum2)
      endif

        if( prth.ge.1 )  &
        pdata(np,prth) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,th)
        if( prt.ge.1 )  &
        pdata(np,prt ) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,t)
        if( prprs.ge.1 )  &
        pdata(np,prprs) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,prs)
        if( prrho.ge.1 )  &
        pdata(np,prrho) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,rho)

        if(prpt1.ge.1)then
          do n=1,npt
          pdata(np,prpt1+n-1) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pt3d(ib,jb,kb,n))
          enddo
        endif

!----------------------------------------------------------------------
!  Data on zvort points

      IF( przv.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( .not. terrain_flag )then
          if( pdata(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( pdata(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = rxu
        ry = ryv
        rz = rzs

        ! w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        ! w2=rx*(1.0-ry)*(1.0-rz)
        ! w3=(1.0-rx)*ry*(1.0-rz)
        ! w4=(1.0-rx)*(1.0-ry)*rz
        ! w5=rx*(1.0-ry)*rz
        ! w6=(1.0-rx)*ry*rz
        ! w7=rx*ry*(1.0-rz)
        ! w8=rx*ry*rz
        call calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        pdata(np,przv) = tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zv)

      ENDIF

!----------------------------------------------------------------------
!  surface variables  and  uv for parcels below lowest model level:

      IF( prznt.ge.1 .or. prust.ge.1 .or. bbc.eq.3 )THEN
        i=iflag
        j=jflag
        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        z0 = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt)
        if( prznt.ge.1 ) pdata(np,prznt) = z0
        if( prust.ge.1 )  &
        pdata(np,prust) = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,ust)
      ENDIF

      IF( bbc.eq.3 )THEN
        ! semi-slip lower boundary condition:
        if( z3d.lt.zh(1,1,1) )then
          ! re-calculate velocities if parcel is below lowest model level:
          !------
          ! u at lowest model level:
          i=iflag
          j=jflag
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          uval = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,u3d(ib,jb,1))
          !------
          ! v at lowest model level:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          vval = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,v3d(ib,jb,1))
          !------
          ! get u,v from (neutral) log-layer equation:
          rznt = 1.0/z0
          var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
          uval = (uval+umove)*var
          vval = (vval+vmove)*var
          !------
          IF( przv.ge.1 )THEN
            do j=jflag-1,jflag+1
            do i=iflag  ,iflag+1
              z0 = 0.5*(znt(i-1,j)+znt(i,j))
              rznt = 1.0/z0
              dum3(i,j,1) = (u3d(i,j,1)+umove)*alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
            enddo
            enddo
            do j=jflag  ,jflag+1
            do i=iflag-1,iflag+1
              z0 = 0.5*(znt(i,j-1)+znt(i,j))
              rznt = 1.0/z0
              dum4(i,j,1) = (v3d(i,j,1)+vmove)*alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
            enddo
            enddo
            do j=jflag,jflag+1
            do i=iflag,iflag+1
              dum7(i,j,1) = (dum4(i,j,1)-dum4(i-1,j,1))*rdx*uf(i)   &
                           -(dum3(i,j,1)-dum3(i,j-1,1))*rdy*vf(j)
            enddo
            enddo
            i=iflag
            j=jflag
            pdata(np,przv) = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 1, 0, 0,dum7(ib,jb,1))
          ENDIF
        endif
      ENDIF


!----------------------------------------------------------------------

        pdata(np,pru)=uval
        pdata(np,prv)=vval
        pdata(np,prw)=wval

#ifdef MPI
      ELSE

        ! set to really small number (so we can use the allreduce command below)
        do n=npvars+1,npvals
          pdata(np,n) = neg_huge 
        enddo
#endif

      ENDIF  myprcl

    ENDDO  nploop2
    if(timestats.ge.1) time_parceli=time_parceli+mytime()

!----------------------------------------------------------------------


      end subroutine parcel_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_write(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc)
      use input, only : maxq,maxvars,nparcelsLocal,npvals,output_format, &
                        myid,dowr,maxstring,string,outfile,numprocs,ierr, &
                        nparcelsLocalActive,nparcelsInit
#ifdef MPI
      use mpi
      use comm_droplet_module, only : CollectDropCount, CollectDropInfo
#endif
#ifdef NETCDF
      use writeout_nc_module, only : writepdata_nc
#endif
      use constants, only : write_frac_droplet, i8

      implicit none

      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      character(len=6) :: prec_string
      real, intent(in), dimension(nparcelsLocal,npvals) :: pdata
      real, intent(inout), dimension(:,:) :: ploc

      integer :: i

      integer(i8) :: nparcelstot
      integer :: nparcelsout     ! this number must be able to be represented
                                 ! by an integer type to reduce memory
                                 ! pressure on I/O
#ifdef MPI
      integer :: nparcels_per_mpi(numprocs)
      real, dimension(:,:), allocatable :: rbuf
#endif

!----------------------------------------------------------------------
!  write out data

    !$acc update host(pdata)

#ifdef MPI
    ! Collect number of active droplets from each MPI rank
    call CollectDropCount ( nparcels_per_mpi, nparcelstot )

    ! Collect information of active droplets from each MPI rank
    if ( myid .eq. 0 ) then
       nparcelsout = int( max( 1.0, write_frac_droplet*nparcelstot ) )
       allocate(rbuf(nparcelsout,npvals))
    end if
    call CollectDropInfo ( pdata, rbuf, nparcels_per_mpi, nparcelsout )
#else
    nparcelstot = int(nparcelsLocalActive, i8)
    nparcelsout = nparcelsLocalActive
#endif

    IF(myid.eq.0)THEN

      IF(output_format.eq.1)THEN
        ! GrADS format:

        ! write GrADS descriptor file:
        call write_prclctl(name_prcl,desc_prcl,unit_prcl,prec)

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        !string = 'cm1out_pdata.dat'
        write(prec_string,'(i6.6)') prec
        string = 'cm1out_pdata_'//prec_string//'.dat'

        if (dowr) write(outfile,*) string
        open(unit=61,file=string,form='unformatted',status='unknown')

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  pdata prec = ',prec
        if(dowr) write(outfile,*) '  writing for total number of parcels:',nparcelstot

        write(61) nparcelsout 
#ifdef MPI
        write(61) rbuf(1:nparcelsout,1:npvals)
#else
        write(61) pdata(1:nparcelsLocalActive,1:npvals)
#endif

        close(unit=61)

#ifdef MPI     
        deallocate(rbuf)
#endif

#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN

        !call     writepdata_nc(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc(1,1))

#endif
      ENDIF
      if(dowr) write(outfile,*)

    ENDIF   ! endif for myid=0

      return
      end subroutine parcel_write


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      real function tri_interp(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s)
      !$acc routine seq
      use input, only : ngxy,ngz
      implicit none

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      !real, dimension(ib:ie,jb:je,kb:ke) :: s
      ! The following can run as many as 8400 pracels
      !real, dimension(ib:iz,jb:jz,kb:kz) :: s

      tri_interp = s(i  ,j  ,k  )*w1    &
                 + s(i+1,j  ,k  )*w2    &
                 + s(i  ,j+1,k  )*w3    &
                 + s(i  ,j  ,k+1)*w4    &
                 + s(i+1,j  ,k+1)*w5    &
                 + s(i  ,j+1,k+1)*w6    &
                 + s(i+1,j+1,k  )*w7    &
                 + s(i+1,j+1,k+1)*w8

      end function tri_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    real function get2d(i,j,x3d,y3d,xh,xf,yh,yf,xs,ys,is,js,s)
    !$acc routine seq
    use input, only : ib,ie,jb,je
    implicit none

    integer, intent(in) :: i,j
    real, intent(in) :: x3d,y3d
    real, intent(in), dimension(ib:ie) :: xh
    real, intent(in), dimension(ib:ie+1) :: xf
    real, intent(in), dimension(jb:je) :: yh
    real, intent(in), dimension(jb:je+1) :: yf

    ! 0 = scalar point
    ! 1 = velocity point
    integer, intent(in) :: xs,ys
    integer, intent(in) :: is,js

    real, intent(in), dimension(ib:ie+is,jb:je+js) :: s

    real :: wg1,wg2,wg3,wg4
    real :: x13,x23,x33,x43
    real :: w1,w2,w3,w7,rx,ry,rz

    logical, parameter :: debug = .false.

!-----------------------------------------------------------------------
      ! tri-linear interp:

      IF(xs.eq.1)THEN
        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
      ELSE
        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
      ENDIF

      IF(ys.eq.1)THEN
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
      ELSE
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
      ENDIF

        w1=(1.0-rx)*(1.0-ry)
        w2=rx*(1.0-ry)
        w3=(1.0-rx)*ry
        w7=rx*ry

      IF( debug )THEN
        if( rx.lt.-0.000001 .or. rx.gt.1.000001 .or.        &
            ry.lt.-0.000001 .or. ry.gt.1.000001 .or.        &
            (w1+w2+w3+w7).lt.0.999999 .or.  &
            (w1+w2+w3+w7).gt.1.000001       &
          )then
          print *,'  x3d,y3d     = ',x3d,y3d
          print *,'  i,j         = ',i,j
          print *,'  rx,ry       = ',rx,ry
          print *,'  w1,w2,w3,w7 = ',w1,w2,w3,w7
          print *,'  w1+w2+w3+w7 = ',w1+w2+w3+w7
          print *,' 22346 '
          call stopcm1
        endif
      ENDIF

      get2d =s(i  ,j  )*w1    &
           +s(i+1,j  )*w2    &
           +s(i  ,j+1)*w3    &
           +s(i+1,j+1)*w7

!-----------------------------------------------------------------------

    end function get2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getparcelzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,pdata)
      use input, only : ib,ie,jb,je,nparcelsLocal,npvals,ni,nj,nx,ny, &
          myi,myj,axisymm,nodex,nodey,ierr, &
          prx,pry,przs
      use constants
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata

      integer :: i,j,iflag,jflag,np
      real :: x3d,y3d

    zsnploop:  &
    DO np=1,nparcelsLocal

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)

      iflag = undefined_index
      jflag = undefined_index

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  zshaveit1:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  zshaveit1

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      zsmyparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                       (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        pdata(np,przs) = get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,zs)

      ELSE  zsmyparcel

        pdata(np,przs) = -1.0e30

      ENDIF  zsmyparcel

    ENDDO  zsnploop

      end subroutine getparcelzs


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setup_parcel_vars(name_prcl,desc_prcl,unit_prcl,qname)
      use input, only : maxq,maxvars,npt,nql1,nql2,iice,nqs1,nqs2,nnc1,nnc2, &
          prcl_out,prcl_droplet,prth,prt,prprs,prpt1,prqv,prq1,prnc1,prkm,prkh, &
          prtke,prdbz,prb,prvpg,przv,prrho,prqsl,prqsi,prznt,prust,przs,prsig
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      character(len=3), intent(in), dimension(maxq) :: qname

      integer :: n,n2
      character(len=8) :: text1
      character(len=30) :: text2

      prcl_out = 0

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'x'
      desc_prcl(prcl_out) = 'x position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'y'
      desc_prcl(prcl_out) = 'y position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'z'
      desc_prcl(prcl_out) = 'z position (above sea level)'
      unit_prcl(prcl_out) = 'm'

      if (prcl_droplet.ge.1) then
         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'vpx'
         desc_prcl(prcl_out) = 'u velocity of droplet'
         unit_prcl(prcl_out) = 'm/s'

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'vpy'
         desc_prcl(prcl_out) = 'v velocity of droplet'
         unit_prcl(prcl_out) = 'm/s'

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'vpz'
         desc_prcl(prcl_out) = 'w velocity of droplet'
         unit_prcl(prcl_out) = 'm/s'

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'rp'
         desc_prcl(prcl_out) = 'radius of droplet'
         unit_prcl(prcl_out) = 'm'

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'ms'
         desc_prcl(prcl_out) = 'solute mass in droplet'
         unit_prcl(prcl_out) = 'kg'

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'tp'
         desc_prcl(prcl_out) = 'temperature of droplet'
         unit_prcl(prcl_out) = 'K'

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'mult'
         desc_prcl(prcl_out) = 'droplet multiplicity'
         unit_prcl(prcl_out) = ' '

         prcl_out = prcl_out+1
         name_prcl(prcl_out) = 'alive'
         desc_prcl(prcl_out) = '<0 = dead, >0 = alive'
         unit_prcl(prcl_out) = ' '

      endif

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'u'
      desc_prcl(prcl_out) = 'u velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'v'
      desc_prcl(prcl_out) = 'v velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'w'
      desc_prcl(prcl_out) = 'w velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'mtime'
      desc_prcl(prcl_out) = 'model time (seconds since beginning of simulation)'
      unit_prcl(prcl_out) = 's'

      if( prth.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'th'
        desc_prcl(prcl_out) = 'potential temperature'
        unit_prcl(prcl_out) = 'K'
      endif

      if( prt.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 't'
        desc_prcl(prcl_out) = 'temperature'
        unit_prcl(prcl_out) = 'K'
      endif

      if( prprs.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'prs'
        desc_prcl(prcl_out) = 'pressure'
        unit_prcl(prcl_out) = 'Pa'
      endif

      if(prpt1.ge.1)then
        do n=1,npt
          text1='pt      '
          if(n.le.9)then
            write(text1(3:3),155) n
155         format(i1.1)
          elseif(n.le.99)then
            write(text1(3:4),154) n
154         format(i2.2)
          else
            write(text1(3:5),153) n
153         format(i3.3)
          endif

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = 'passive tracer mixing ratio'
          unit_prcl(prcl_out) = 'kg/kg'
        enddo
      endif

      if( prqv.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qv'
        desc_prcl(prcl_out) = 'water vapor mixing ratio'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prq1.ge.1 .and. nql1.ge.1 )then
        n2 = nql2
        if( iice.eq.1 ) n2 = nqs2
        do n=nql1,n2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = text2
          unit_prcl(prcl_out) = 'kg/kg'
        enddo
      endif

      if(prnc1.ge.1)then
        do n=nnc1,nnc2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = text2
          unit_prcl(prcl_out) = '1/kg'
        enddo
      endif


      if( prkm.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'kmh'
        desc_prcl(prcl_out) = 'horiz eddy viscosity for momentum'
        unit_prcl(prcl_out) = 'm^2/s'

        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'kmv'
        desc_prcl(prcl_out) = 'vert eddy viscosity for momentum'
        unit_prcl(prcl_out) = 'm^2/s'
      endif

      if( prkh.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'khh'
        desc_prcl(prcl_out) = 'horiz eddy diffusivity for scalars'
        unit_prcl(prcl_out) = 'm^2/s'

        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'khv'
        desc_prcl(prcl_out) = 'vert eddy diffusivity for scalars'
        unit_prcl(prcl_out) = 'm^2/s'
      endif

      if( prtke.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'tke'
        desc_prcl(prcl_out) = 'subgrid tke'
        unit_prcl(prcl_out) = 'm^2/s^2'
      endif

      if( prdbz.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'dbz'
        desc_prcl(prcl_out) = 'reflectivity'
        unit_prcl(prcl_out) = 'dBZ'
      endif

      if( prb.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'b'
        desc_prcl(prcl_out) = 'buoyancy'
        unit_prcl(prcl_out) = 'm/s/s'
      endif

      if( prvpg.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'vpg'
        desc_prcl(prcl_out) = 'vertical perturbation pressure gradient'
        unit_prcl(prcl_out) = 'm/s/s'
      endif

      if( przv.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'zvort'
        desc_prcl(prcl_out) = 'vertical vorticity'
        unit_prcl(prcl_out) = '1/s'
      endif

      if( prrho.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'rho'
        desc_prcl(prcl_out) = 'dry-air density'
        unit_prcl(prcl_out) = 'kg/m^3'
      endif

      if( prqsl.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qsl'
        desc_prcl(prcl_out) = 'saturation mixing ratio wrt liquid'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prqsi.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qsi'
        desc_prcl(prcl_out) = 'saturation mixing ratio wrt ice'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prznt.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'znt'
        desc_prcl(prcl_out) = 'surface roughness length'
        unit_prcl(prcl_out) = 'm'
      endif

      if( prust.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'ust'
        desc_prcl(prcl_out) = 'surface friction velocity'
        unit_prcl(prcl_out) = 'm/s'
      endif

      if( przs.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'zs'
        desc_prcl(prcl_out) = 'terrain height'
        unit_prcl(prcl_out) = 'm'
      endif

      if( prsig.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'sigma'
        desc_prcl(prcl_out) = 'sigma (nondimensional height)'
        unit_prcl(prcl_out) = 'nondimensional'
      endif


!-----------------------------------------------------------------------

      end subroutine setup_parcel_vars


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_prclctl(name_prcl,desc_prcl,unit_prcl,prec)
      use input, only : maxvars,myid,string,maxstring,dowr,cm1version, &
          nparcelsInit,prcl_out,outfile
      use constants , only : grads_undef, i8
      implicit none

      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      integer, intent(in) :: prec

      integer :: i,n,nn
      character(len=16) :: a16

      !---------------------------------------------------------------
      ! This subroutine writes the GrADS descriptor file for parcels
      !---------------------------------------------------------------

    idcheck:  &
    IF( myid.eq.0 )THEN

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string = 'cm1out_parcel.ctl'
      if(dowr) write(outfile,*) string
      open(unit=50,file=string,status='unknown')

      string = 'cm1out_pdata.dat'

      write(50,401) string
      write(50,402) trim(cm1version)
      write(50,403) grads_undef
      write(50,404) nparcelsInit
      write(50,405)
      write(50,406)
      write(50,407) prec

      write(50,408) prcl_out

      DO n = 1 , prcl_out
        a16 = '                '
        nn = len(trim(unit_prcl(n)))
        write(a16(2:15),214) unit_prcl(n)
        write(a16(1:1),201 )       '('
        write(a16(nn+2:nn+2),201 ) ')'
        write(50,409) name_prcl(n),desc_prcl(n),a16
      ENDDO

      write(50,410)

      close(unit=50)

    ENDIF  idcheck

201   format(a1)
214   format(a14)

401   format('dset ^',a)
402   format('title CM1 parcel data output, using version ',a,'; time is generic, see variable mtime for actual times')
403   format('undef ',f10.1)
404   format('xdef ',i15,' linear 1 1')
405   format('ydef          1 linear 0 1')
406   format('zdef          1 linear 0 1')
407   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
408   format('vars ',i6)
409   format(a12,' 1 99 ',a40,1x,a16)
410   format('endvars')

      end subroutine write_prclctl

      real function calcW1(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW1=(1.0-rx)*(1.0-ry)*(1.0-rz)
      end function calcW1

      real function calcW2(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW2=rx*(1.0-ry)*(1.0-rz)
      end function calcW2

      real function calcW3(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW3=(1.0-rx)*ry*(1.0-rz)
      end function calcW3

      real function calcW4(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW4=(1.0-rx)*(1.0-ry)*rz
      end function calcW4

      real function calcW5(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW5=rx*(1.0-ry)*rz
      end function calcW5

      real function calcW6(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW6=(1.0-rx)*ry*rz
      end function calcW6

      real function calcW7(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW7=rx*ry*(1.0-rz)
      end function calcW7

      real function calcW8(rx,ry,rz)
      !$acc routine vector
      real, intent(in) :: rx,ry,rz
        calcW8=rx*ry*rz
      end function calcW8

      subroutine calcWeights(w1,w2,w3,w4,w5,w6,w7,w8,rx,ry,rz)
      !$acc routine vector
        real, intent(in)  :: rx,ry,rz
        real, intent(inout) :: w1,w2,w3,w4,w5,w6,w7,w8

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz
 
      end subroutine calcWeights

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! These two subroutines are added by following John Dennis's suggestion
      !      to optimize the location search for a parcel in a process


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE parcel_module
