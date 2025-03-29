* Script to calculate psuedoadiabatic equivalent potential temperature
* (the) from cm1 GrADS-format output
* Last modified:  31 May 2007

"define e=0.01*prs*qv/(0.6225596+qv)+1e-20"
"define t=th*pow((prs/100000.0),(287.0/1004.0))"
"define tlcl=55.0+2840.0/(3.5*log(t)-log(e)-4.805)"
"define the=t*pow(100000.0/prs,0.2854*(1.0-0.28*qv))*exp(((3376.0/tlcl)-2.54)*qv*(1.0+0.81*qv))"
