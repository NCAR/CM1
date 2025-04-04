* Script to calculate density (rho) from cm1 GrADS-format output
* Author:  George H. Bryan, NCAR/MMM
* Last modified:  1 June 2007

"define t=th*pow((prs/100000.0),(287.0/1004.0))"
"define rho=(prs/(287.0*t*(1.0+1.606272*qv)))"
