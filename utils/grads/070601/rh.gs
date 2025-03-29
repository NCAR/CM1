* Script to calculate relative humidity wrt water (rh) from 
* cm1 GrADS-format output
* Last modified:  31 May 2007

"define t=th*pow((prs/100000.0),(287.0/1004.0))"
"define qvs=380.00*exp(17.2693882-4097.8531/(t-35.86))/prs"
"define rh=100.0*qv/qvs"
