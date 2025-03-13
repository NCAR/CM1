  MODULE radlib3d_module

  implicit none

  private
  public :: coeff,terp1,expmn,deledd,sagpol

  CONTAINS
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

SUBROUTINE coeff(n,x,f,w,iop,INT,wk)
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
  IMPLICIT NONE

  INTEGER :: n
  INTEGER :: iop(2)

  REAL :: x(n),f(n),w(n),wk(n,4)

  INTEGER :: i
  INTEGER :: i1,i2, ii,ii1, in,INT
  INTEGER :: j0,j1,j2,j3,j4
  INTEGER :: jm,mk
  INTEGER :: ml,nn

  REAL :: a12,a13,a14,a23,a24,a34, b2,y2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ii(i)=(i-1)*INT+1
  j0=1
  DO i=2,n
    jm=j0
    j0=j0+INT
    wk(i,1)=x(i)-x(i-1)
    wk(i,2)=(f(j0)-f(jm))/wk(i,1)
    wk(i,3)=wk(i,1)/6.
    wk(i,1)=wk(i,1)/3.
  END DO
  nn=n
  mk=iop(1)
  ml=iop(2)
!  GO TO (102,103,104,105) ,mk
!  obsolescent feature correction by WYH
  SELECT CASE (mk)
  CASE (1)
    GO TO 102
  CASE (2)
    GO TO 103
  CASE (3)
    GO TO 104
  CASE (4)
    GO TO 105
  CASE DEFAULT
!   Do nothing
  END SELECT
! end of correction by wyh.

  102 CONTINUE
  wk(2,2)=wk(3,2)-wk(2,2)-wk(2,3)*w(1)
  wk(2,3)=0.
  wk(2,1)=wk(2,1)+wk(3,1)
  i1=2
  nn=nn-1
  GO TO 106
  103 CONTINUE
  wk(1,2)=wk(2,2)-w(1)
  wk(2,2)=wk(3,2)-wk(2,2)
  wk(1,3)=0.
  wk(1,1)=wk(2,1)
  wk(2,1)=wk(2,1)+wk(3,1)
  i1=1
  GO TO 106
  104 CONTINUE
  y2=wk(2,2)
  b2=wk(2,1)
  wk(2,2)=wk(3,2)-wk(2,2)
  wk(2,1)=wk(3,1)+wk(2,1)
  i1=2
  nn=nn-1
  GO TO 106
  105 CONTINUE
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
  GO TO 103
  106 CONTINUE
  i2=n-2
  DO i=3,i2
    wk(i,2)=wk(i+1,2)-wk(i,2)
    wk(i,1)=wk(i+1,1)+wk(i,1)
  END DO
  in=ii(n)
!  GO TO (108,109,110,111) ,ml
!  obsolescent feature correction by WYH
  SELECT CASE (ml)
  CASE (1)
    GO TO 108
  CASE (2)
    GO TO 109
  CASE (3)
    GO TO 110
  CASE (4)
    GO TO 111
  CASE DEFAULT
!   Do nothing
  END SELECT
! end of correction by WYH

  108 CONTINUE
  wk(n-1,2)=wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
  wk(n,3)=0.
  wk(n-1,1)=wk(n-1,1)+wk(n,1)
  nn=nn-1
  GO TO 112
  109 CONTINUE
  wk(n-1,2)=wk(n,2)-wk(n-1,2)
  wk(n,2)=-wk(n,2)+w(in)
  wk(n-1,1)=wk(n-1,1)+wk(n,1)
  wk(1,4)=0.
  GO TO 112
  110 CONTINUE
  wk(n-1,2)=wk(n,2)-wk(n-1,2)
  wk(n,2)=y2-wk(n,2)
  wk(n-1,1)=wk(n-1,1)+wk(n,1)
  wk(n,1)=wk(n,1)+b2
  wk(1,4)=wk(2,3)
  GO TO 112
  111 CONTINUE
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
  GO TO 109
  112 CONTINUE
  ii1=ii(i1)
  CALL trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),w(ii1),INT)
!  GO TO (114,114,113,114) ,mk
!  obsolescent feature correction by WYH
  SELECT CASE (mk)
  CASE (1,2,4)
    GO TO 114
  CASE (3)
    GO TO 113
  CASE DEFAULT
!   Do nothing
  END SELECT
! end of correction by wyh

  113 CONTINUE
  w(1)=w(in)
  114 CONTINUE

  RETURN
END SUBROUTINE coeff
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

SUBROUTINE terp1(n,x,f,w,y,INT,tab,itab)
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
  IMPLICIT NONE

  INTEGER :: n
  INTEGER :: INT
  INTEGER :: itab(3)

  REAL :: x(n),f(n),w(n),tab(3)

  INTEGER :: i
  REAL :: y
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL search (y,x,n,i)
  CALL interp (n,x,f,w,y,i,INT,tab,itab)

  RETURN
END SUBROUTINE terp1
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

SUBROUTINE trip(n,a,b,c,y,z,INT)
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
  IMPLICIT NONE

  INTEGER :: n

  REAL :: a(n),b(n),c(n),y(n),z(n)

  INTEGER :: INT

  INTEGER :: ii, ik, in, inm
  INTEGER :: i,j,k
  INTEGER :: nm1, nm2

  REAL :: bn, yn, v, den
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
  DO j=2,nm2
    den=b(j)-a(j)*b(j-1)
    b(j)=c(j)/den
    y(j)=(y(j)-a(j)*y(j-1))/den
    a(j)=-a(j)*a(j-1)/den
    bn=bn-v*a(j-1)
    yn=yn-v*y(j-1)
    v=-v*b(j-1)
  END DO
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
  DO j=2,nm1
    k=n-j
    ik=ii(k)
    z(ik)=y(k)-b(k)*z(ik+INT)-a(k)*z(in)
  END DO

  RETURN
END SUBROUTINE trip
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

SUBROUTINE search(xbar,x,n,i)
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
  IMPLICIT NONE

  INTEGER :: i,n

  REAL :: x(n)
  REAL :: xbar

  INTEGER :: k,m,nm1

  REAL :: b
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  b = .69314718

  IF(xbar > x(2)) GO TO 101
  i=1
  RETURN

  101 CONTINUE
  IF(xbar < x(n-1)) GO TO 102
  i=n-1
  RETURN

  102 CONTINUE
  m=INT((ALOG(FLOAT(n)))/b)
  i=2**m
  IF(i >= n) i=i/2
  k=i
  nm1=n-1
  103 CONTINUE
  k=k/2
  IF(xbar >= x(i)) GO TO 104
  i=i-k
  GO TO 103
  104 CONTINUE
  IF(xbar <= x(i+1)) RETURN
  i=MIN0(i+k,nm1)
  GO TO 103

END SUBROUTINE search
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

SUBROUTINE interp(n,x,f,w,y,i,INT,tab,itab)
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
  IMPLICIT NONE

  INTEGER :: n
  INTEGER :: i, i0, ii, INT, ip

  INTEGER :: itab(3)
  REAL :: x(n),f(n),w(n),tab(3)
  REAL :: a,b,c, fl0,flk,flp, y
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
  IF(itab(1) /= 1) GO TO 102
  a=(w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)
  b=(f(ip)/flk-w(ip)*flk/6.)*fl0
  c=(f(i0)/flk-w(i0)*flk/6.)*flp
  tab(1)=a+b+c
  102 IF(itab(2) /= 1) GO TO 104
  a=(w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
  b=(f(ip)-f(i0))/flk
  c=(w(i0)-w(ip))*flk/6.
  tab(2)=a+b+c
  104 IF(itab(3) /= 1) GO TO 106
  tab(3)=(w(i0)*flp+w(ip)*fl0)/flk

  106 RETURN
END SUBROUTINE interp
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

  REAL FUNCTION expmn(fin)
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
  IMPLICIT NONE

  REAL :: one, expmin
  PARAMETER (one=1.0, expmin=-10.0)

  REAL :: e1,e2,e3,e4
  PARAMETER (e1=1.0,        e2=-2.507213E-1)
  PARAMETER (e3=2.92732E-2, e4=-3.827800E-3)

  REAL :: fin, tmp
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

  RETURN
  END FUNCTION expmn
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

SUBROUTINE deledd(tau,ssc,g0,csm,rr,tt,td)
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
  IMPLICIT NONE

  REAL :: zero,one,two,three,four,fourth,seven,thresh
  PARAMETER (one =1., three=3.)
  PARAMETER (two =2., seven=7.)
  PARAMETER (four=4., fourth=.25)
  PARAMETER (zero=0., thresh=1.e-8)

!-----input parameters
  REAL :: tau,ssc,g0,csm

!-----output parameters
  REAL :: rr,tt,td

!-----temporary parameters

  REAL :: zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2,             &
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

  IF (ABS(st3) < thresh) THEN
    zth = zth + 0.001
    xx  = akk * zth
    st7 = one - xx
    st8 = one + xx
    st3 = st7 * st8
  END IF

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

  RETURN
END SUBROUTINE deledd
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

SUBROUTINE sagpol(tau,ssc,g0,rll,tll)
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
  IMPLICIT NONE

  REAL :: one,three,four
  PARAMETER (one=1., three=3., four=4.)

!-----output parameters:

  REAL :: tau,ssc,g0

!-----output parameters:

  REAL :: rll,tll

!-----temporary arrays

  REAL :: xx,uuu,ttt,emt,up1,um1,st1

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

  RETURN
END SUBROUTINE sagpol


  END MODULE radlib3d_module
