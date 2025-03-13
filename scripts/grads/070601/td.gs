* Script to calculate dewpoint temperature (td) from cm1 GrADS-format output
* Last modified:  31 May 2007

"define el=log((qv/0.622)*prs/100.0/(1.0+(qv/0.622)))"
"define td=273.15+(243.5*el-440.8)/(19.48-el)"
