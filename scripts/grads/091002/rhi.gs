* Script to calculate relative humidity wrt ice (rhi) from
* cm1 GrADS-format output
* Author:  George H. Bryan, NCAR/MMM
* Last modified:  8 August 2008

"define t=th*pow((prs/100000.0),(287.04/1005.7))"
"define esi=611.2*exp(21.8745584*(t-273.15)/(t-7.66))"
"define qvsi=0.6219718*esi/(prs-esi)"
"define rhi=100.0*qv/qvsi"
