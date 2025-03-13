* Script to calculate composite reflectivity (i.e., max dbz in column)
* NOTE:  requires script "dbz.gs"
* Author:  George H. Bryan, NCAR/MMM
* Last modified:  1 June 2007

"q file"
res=sublin(result,5)
nz=subwrd(res,9)
say ' nz = ' nz

zz=1
"define cref=0*th"
while(zz <= nz)
  say ' z = ' zz
  "set z "zz
  "run dbz"
  "define cref=const(maskout(cref,cref-dbz),0.0,-u)+const(maskout(dbz,dbz-cref),0.0,-u)"
  zz=zz+1
endwhile

