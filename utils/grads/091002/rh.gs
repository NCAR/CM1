* Script to calculate relative humidity wrt water (rh) from 
* cm1 GrADS-format output
* Author:  George H. Bryan, NCAR/MMM
* Last modified:  8 August 2008

"define t=th*pow((prs/100000.0),(287.04/1005.7))"
"define esl=611.2*exp(17.67*(t-273.15 )/(t-29.65))"
"define qvs=0.6219718*esl/(prs-esl)"
"define rh=100.0*qv/qvs"
