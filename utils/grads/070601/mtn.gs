* Script to plot underlying terrain on cross sections from cm1 
* GrADS-format output that is interpolated to height levels 
* (i.e., when output_interp=1)
* Requires output_zh = 1
* Last modified:  31 May 2007

"q gxinfo"
res=sublin(result,1)
gxout=subwrd(res,4)

"set gxout contour"
"set clab off"
"set clevs 0"
"set ccolor 1"
"d zh"
"set gxout "gxout

"set ccolor 1"

