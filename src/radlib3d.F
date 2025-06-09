  module radlib3d_module

  implicit none

  private
  public :: coeff,terp1,expmn,deledd,sagpol

  contains
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE COEFF                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine coeff(n,x,f,w,iop,INT,wk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  implicit none

  integer :: n
  integer :: iop(2)

  real :: x(n),f(n),w(n),wk(n,4)

  integer :: i
  integer :: i1,i2, ii,ii1, in,INT
  integer :: j0,j1,j2,j3,j4
  integer :: jm,mk
  integer :: ml,nn

  real :: a12,a13,a14,a23,a24,a34, b2,y2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ii(i)=(i-1)*INT+1
  j0=1
  do i=2,n
    jm=j0
    j0=j0+INT
    wk(i,1)=x(i)-x(i-1)
    wk(i,2)=(f(j0)-f(jm))/wk(i,1)
    wk(i,3)=wk(i,1)/6.
    wk(i,1)=wk(i,1)/3.
  end do
  nn=n
  mk=iop(1)
  ml=iop(2)
!  GO TO (102,103,104,105) ,mk
!  obsolescent feature correction by WYH
  select case (mk)
  case (1)
    go to 102
  case (2)
    go to 103
  case (3)
    go to 104
  case (4)
    go to 105
  case DEFAULT
!   Do nothing
  end select
! end of correction by wyh.

  102 continue
  wk(2,2)=wk(3,2)-wk(2,2)-wk(2,3)*w(1)
  wk(2,3)=0.
  wk(2,1)=wk(2,1)+wk(3,1)
  i1=2
  nn=nn-1
  go to 106
  103 continue
  wk(1,2)=wk(2,2)-w(1)
  wk(2,2)=wk(3,2)-wk(2,2)
  wk(1,3)=0.
  wk(1,1)=wk(2,1)
  wk(2,1)=wk(2,1)+wk(3,1)
  i1=1
  go to 106
  104 continue
  y2=wk(2,2)
  b2=wk(2,1)
  wk(2,2)=wk(3,2)-wk(2,2)
  wk(2,1)=wk(3,1)+wk(2,1)
  i1=2
  nn=nn-1
  go to 106
  105 continue
  a12=x(1)-x(2)
  a13=x(1)-x(3)
  a14=x(1)-x(4)
  a23=x(2)-x(3)
  a24=x(2)-x(4)
  a34=x(3)-x(4)
  j1=1
  j2=j1+INT
  j3=j2+INT
  j4=j3+INT
  w(1)=(1./a12+1./a13+1./a14)*f(j1)-                                    &
       a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-         &
       a12*a13/(a14*a24*a34)*f(j4)
  go to 103
  106 continue
  i2=n-2
  do i=3,i2
    wk(i,2)=wk(i+1,2)-wk(i,2)
    wk(i,1)=wk(i+1,1)+wk(i,1)
  end do
  in=ii(n)
!  GO TO (108,109,110,111) ,ml
!  obsolescent feature correction by WYH
  select case (ml)
  case (1)
    go to 108
  case (2)
    go to 109
  case (3)
    go to 110
  case (4)
    go to 111
  case DEFAULT
!   Do nothing
  end select
! end of correction by WYH

  108 continue
  wk(n-1,2)=wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
  wk(n,3)=0.
  wk(n-1,1)=wk(n-1,1)+wk(n,1)
  nn=nn-1
  go to 112
  109 continue
  wk(n-1,2)=wk(n,2)-wk(n-1,2)
  wk(n,2)=-wk(n,2)+w(in)
  wk(n-1,1)=wk(n-1,1)+wk(n,1)
  wk(1,4)=0.
  go to 112
  110 continue
  wk(n-1,2)=wk(n,2)-wk(n-1,2)
  wk(n,2)=y2-wk(n,2)
  wk(n-1,1)=wk(n-1,1)+wk(n,1)
  wk(n,1)=wk(n,1)+b2
  wk(1,4)=wk(2,3)
  go to 112
  111 continue
  a12=x(n)-x(n-1)
  a13=x(n)-x(n-2)
  a14=x(n)-x(n-3)
  a23=x(n-1)-x(n-2)
  a24=x(n-1)-x(n-3)
  a34=x(n-2)-x(n-3)
  j1=in
  j2=j1-INT
  j3=j2-INT
  j4=j3-INT
  w(in)=(1./a12+1./a13+1./a14)*f(j1)-                                   &
        a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-        &
        a12*a13/(a14*a24*a34)*f(j4)
  go to 109
  112 continue
  ii1=ii(i1)
  call trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),w(ii1),INT)
!  GO TO (114,114,113,114) ,mk
!  obsolescent feature correction by WYH
  select case (mk)
  case (1,2,4)
    go to 114
  case (3)
    go to 113
  case DEFAULT
!   Do nothing
  end select
! end of correction by wyh

  113 continue
  w(1)=w(in)
  114 continue

  return
end subroutine coeff
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TERP1                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine terp1(n,x,f,w,y,INT,tab,itab)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!fpp$ expand (search)
!fpp$ expand (interp)
!!dir$ inline always interp, search
!*$*  inline routine (interp, search)
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  implicit none

  integer :: n
  integer :: INT
  integer :: itab(3)

  real :: x(n),f(n),w(n),tab(3)

  integer :: i
  real :: y
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  call search (y,x,n,i)
  call interp (n,x,f,w,y,i,INT,tab,itab)

  return
end subroutine terp1
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRIP                       ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine trip(n,a,b,c,y,z,INT)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  implicit none

  integer :: n

  real :: a(n),b(n),c(n),y(n),z(n)

  integer :: INT

  integer :: ii, ik, in, inm
  integer :: i,j,k
  integer :: nm1, nm2

  real :: bn, yn, v, den
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ii(i)=(i-1)*INT+1
  bn=b(n)
  yn=y(n)
  v=c(n)
  y(1)=y(1)/b(1)
  a(1)=a(1)/b(1)
  b(1)=c(1)/b(1)
  nm2=n-2
  do j=2,nm2
    den=b(j)-a(j)*b(j-1)
    b(j)=c(j)/den
    y(j)=(y(j)-a(j)*y(j-1))/den
    a(j)=-a(j)*a(j-1)/den
    bn=bn-v*a(j-1)
    yn=yn-v*y(j-1)
    v=-v*b(j-1)
  end do
  den=b(n-1)-a(n-1)*b(n-2)
  b(n-1)=(c(n-1)-a(n-1)*a(n-2))/den
  y(n-1)=(y(n-1)-a(n-1)*y(n-2))/den
  bn=bn-v*a(n-2)
  yn=yn-v*y(n-2)
  v=a(n)-v*b(n-2)
  nm1=n-1
  in=ii(n)
  inm=ii(nm1)
  z(in)=(yn-v*y(nm1))/(bn-v*b(nm1))
  z(inm)=y(nm1)-b(nm1)*z(in)
  do j=2,nm1
    k=n-j
    ik=ii(k)
    z(ik)=y(k)-b(k)*z(ik+INT)-a(k)*z(in)
  end do

  return
end subroutine trip
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SEARCH                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine search(xbar,x,n,i)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  implicit none

  integer :: i,n

  real :: x(n)
  real :: xbar

  integer :: k,m,nm1

  real :: b
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  b = .69314718

  if(xbar > x(2)) go to 101
  i=1
  return

  101 continue
  if(xbar < x(n-1)) go to 102
  i=n-1
  return

  102 continue
  m=INT((ALOG(FLOAT(n)))/b)
  i=2**m
  if(i >= n) i=i/2
  k=i
  nm1=n-1
  103 continue
  k=k/2
  if(xbar >= x(i)) go to 104
  i=i-k
  go to 103
  104 continue
  if(xbar <= x(i+1)) return
  i=MIN0(i+k,nm1)
  go to 103

end subroutine search
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTERP                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine interp(n,x,f,w,y,i,INT,tab,itab)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  implicit none

  integer :: n
  integer :: i, i0, ii, INT, ip

  integer :: itab(3)
  real :: x(n),f(n),w(n),tab(3)
  real :: a,b,c, fl0,flk,flp, y
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ii(i)=(i-1)*INT+1
  flk=x(i+1)-x(i)
  flp=x(i+1)-y
  fl0=y-x(i)
  i0=ii(i)
  ip=i0+INT
  if(itab(1) /= 1) go to 102
  a=(w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)
  b=(f(ip)/flk-w(ip)*flk/6.)*fl0
  c=(f(i0)/flk-w(i0)*flk/6.)*flp
  tab(1)=a+b+c
  102 if(itab(2) /= 1) go to 104
  a=(w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
  b=(f(ip)-f(i0))/flk
  c=(w(i0)-w(ip))*flk/6.
  tab(2)=a+b+c
  104 if(itab(3) /= 1) go to 106
  tab(3)=(w(i0)*flp+w(ip)*fl0)/flk

  106 return
end subroutine interp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION EXPMN                        ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

  real function expmn(fin)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate exponential for arguments in the range 0> fin > -10.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!**************************************************************************
! compute exponential for arguments in the range 0> fin > -10.
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  implicit none

  real :: one, expmin
  parameter (one=1.0, expmin=-10.0)

  real :: e1,e2,e3,e4
  parameter (e1=1.0,        e2=-2.507213E-1)
  parameter (e3=2.92732E-2, e4=-3.827800E-3)

  real :: fin, tmp
!!!  REAL :: expmn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tmp = MAX( fin, expmin )
  expmn = ((e4*tmp + e3)*tmp+e2)*tmp+e1
  expmn = expmn * expmn
  expmn = one / (expmn * expmn)

  return
  end function expmn
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DELEDD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine deledd(tau,ssc,g0,csm,rr,tt,td)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the bulk scattering properties of a single layer using
!  delta-eddington approximation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!*********************************************************************
!
!-----uses the delta-eddington approximation to compute the
!  bulk scattering properties of a single layer
!  coded following King and Harshvardhan (JAS, 1986)
!
!  inputs:
!
!  tau: the effective optical thickness
!  ssc: the effective single scattering albedo
!  g0:  the effective asymmetry factor
!  csm: the effective secant of the zenith angle
!
!  outputs:
!
!  rr: the layer reflection of the direct beam
!  tt: the layer diffuse transmission of the direct beam
!  td: the layer direct transmission of the direct beam
!
!*********************************************************************
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!*********************************************************************
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  implicit none

  real :: zero,one,two,three,four,fourth,seven,thresh
  parameter (one =1., three=3.)
  parameter (two =2., seven=7.)
  parameter (four=4., fourth=.25)
  parameter (zero=0., thresh=1.e-8)

!-----input parameters
  real :: tau,ssc,g0,csm

!-----output parameters
  real :: rr,tt,td

!-----temporary parameters

  real :: zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2,             &
       all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4

!!!  REAL, EXTERNAL :: expmn

!---------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  zth = one / csm

!  delta-eddington scaling of single scattering albedo,
!  optical thickness, and asymmetry factor,
!  K & H eqs(27-29)

  ff  = g0*g0
  xx  = one-ff*ssc
  taup= tau*xx
  sscp= ssc*(one-ff)/xx
  gp  = g0/(one+g0)

!  gamma1, gamma2, and gamma3. see table 2 and eq(26) K & H
!  ssc and gp are the d-s single scattering
!  albedo and asymmetry factor.

  xx  =  three*gp
  gm1 =  (seven - sscp*(four+xx))*fourth
  gm2 = -(one   - sscp*(four-xx))*fourth

!  akk is k as defined in eq(25) of K & H

  akk = SQRT((gm1+gm2)*(gm1-gm2))

  xx  = akk * zth
  st7 = one - xx
  st8 = one + xx
  st3 = st7 * st8

  if (ABS(st3) < thresh) then
    zth = zth + 0.001
    xx  = akk * zth
    st7 = one - xx
    st8 = one + xx
    st3 = st7 * st8
  end if

!  extinction of the direct beam transmission

!wdt
  !td  = EXP(-taup/zth)
  td  = expmn(-taup/zth)  ! faster, but less accurate for long integrations

!  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of K & H

  gm3  = (two - zth*three*gp)*fourth
  xx   = gm1 - gm2
  alf1 = gm1 - gm3 * xx
  alf2 = gm2 + gm3 * xx

!  all is last term in eq(21) of K & H
!  bll is last term in eq(22) of K & H

  xx  = akk * two
  all = (gm3 - alf2 * zth    )*xx*td
  bll = (one - gm3 + alf1*zth)*xx

  xx  = akk * gm3
  cll = (alf2 + xx) * st7
  dll = (alf2 - xx) * st8

  xx  = akk * (one-gm3)
  fll = (alf1 + xx) * st8
  ell = (alf1 - xx) * st7

!wdt
  !st2 = EXP(-akk*taup)
  st2 = expmn(-akk*taup)  ! faster, but less accurate for long integrations
  st4 = st2 * st2

  st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)

!  rr is r-hat of eq(21) of K & H
!  tt is diffuse part of t-hat of eq(22) of K & H

  rr =   ( cll-dll*st4    -all*st2)*st1
  tt = - ((fll-ell*st4)*td-bll*st2)*st1

  rr = MAX(rr,zero)
  tt = MAX(tt,zero)

  return
end subroutine deledd
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SAGPOL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

subroutine sagpol(tau,ssc,g0,rll,tll)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate transmittance and reflectance of diffuse radiation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Formatted it in accordance with the ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!*********************************************************************
!-----transmittance (tll) and reflectance (rll) of diffuse radiation
!  follows Sagan and Pollock (JGR, 1967).
!  also, eq.(31) of Lacis and Hansen (JAS, 1974).
!
!-----input parameters:
!
!   tau: the effective optical thickness
!   ssc: the effective single scattering albedo
!   g0:  the effective asymmetry factor
!
!-----output parameters:
!
!   rll: the layer reflection of diffuse radiation
!   tll: the layer transmission of diffuse radiation
!
!*********************************************************************
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  implicit none

  real :: one,three,four
  parameter (one=1., three=3., four=4.)

!-----output parameters:

  real :: tau,ssc,g0

!-----output parameters:

  real :: rll,tll

!-----temporary arrays

  real :: xx,uuu,ttt,emt,up1,um1,st1

!-----function

!!!  REAL :: expmn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xx  = one-ssc*g0
  uuu = SQRT( xx/(one-ssc))
  ttt = SQRT( xx*(one-ssc)*three )*tau
  emt = expmn(-ttt)
  up1 = uuu + one
  um1 = uuu - one
  xx  = um1*emt
  st1 = one / ((up1+xx) * (up1-xx))
  rll = up1*um1*(one-emt*emt)*st1
  tll = uuu*four*emt         *st1

  return
end subroutine sagpol


  end module radlib3d_module
