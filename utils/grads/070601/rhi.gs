* Script to calculate relative humidity wrt ice (rhi) from
* cm1 GrADS-format output
* Last modified:  31 May 2007

"define t=th*pow((prs/100000.0),(287.0/1004.0))"
"define qvsi=380.00*exp(21.87455-5807.4743/(t-7.66))/prs"
"define rhi=100.0*qv/qvsi"
