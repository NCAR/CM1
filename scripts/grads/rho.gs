* Script to calculate density (rho) from cm1 GrADS-format output
* Author:  George H. Bryan, NCAR/MMM
* Last modified:  8 August 2008

"define t=th*pow((prs/100000.0),(287.04/1005.7))"
"define rho=(prs/(287.0*t*(1.0+1.607790*qv)))"
