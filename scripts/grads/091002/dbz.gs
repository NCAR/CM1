* Script to calculate estimated reflectivity (dbz) from cm1 GrADS-format output
* NOTE:  settings are consistent with the Goddard LFO scheme (ptype=2)
* assuming ihail = 1.  For ihail = 0, change setting below.
* Author:  George H. Bryan, NCAR/MMM
* Last modified:  8 August 2008


ihail = 1


"define t=th*pow((prs/100000.0),(287.04/1005.7))"
"define rho=(prs/(287.0*t*(1.0+1.607790*qv)))"

eps = 1.0e-12

*--- rain ---
n0r  = 8.0e6
rhor = 1000.0
"define gamma=pow(3.14159*"n0r"*"rhor"/(rho*qr),0.25)"
"define zer=720.0*"n0r"*pow(gamma,-7.0)"
"define zer=const(zer,0.0,-u)"

*--- graupel/hail ---
if(ihail=1)
  n0g  = 2.0e4
  rhog = 900.0
else
  n0g  = 4.0e6
  rhog = 400.0
endif
"define gamma=pow(3.14159*"n0g"*"rhog"/(rho*qg),0.25)"
"define zeg=720.0*"n0g"*pow(gamma,-7.0)*pow(("rhog"/"rhor"),2.0)*0.224"
"define zeg=const(zeg,0.0,-u)"

*--- snow ---
n0s  = 1.0e8
rhos = 100.0
"define gamma=pow(3.14159*"n0s"*"rhos"/(rho*qs),0.25)"
"define zes=720.0*"n0s"*pow(gamma,-7.0)*pow(("rhos"/"rhor"),2.0)*0.224"
"define zes=const(zes,0.0,-u)"

*--- dbz ---

"define dbz=10.0*log10("eps"+(zer+zeg+zes)*1e18)"
"define dbz=const(maskout(dbz,dbz),0.,-u)"

