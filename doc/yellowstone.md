**Suggestions for using CM1 on yellowstone:**

The new NCAR supercomputing facility, [yellowstone](https://www2.cisl.ucar.edu/resources/yellowstone), is now available. This page provides some suggestions for using CM1, based on some preliminary tests.

* * *

**Compiling CM1:**

In the Makefile, I recommend using the section titled "Linux, distributed memory, Intel compiler (eg, SHARCNET's saw)." I recommend the following compiler flags:

```
OPTS = -I../include -O3 -xHost -ip -assume byterecl -fp-model precise -ftz
```

You will likely see warning messages such as the following:  

param.f90(3677): remark #8290: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+3'.

These messages can be ignored (they only affect text output and have no impact on CM1 performance).

(Note: I haven't tested any of the other compilers that are available on yellowstone. If you have the time/interest in running some performance tests using other compilers, please keep me informed. Thanks!)

Note for users of the atmospheric radition scheme (radopt=1) in CM1: (posted 7 Feb 2013) There appears to be a problem when using the intel fortran compiler on yellowstone. To prevent the problem, add this compiler flag to the Makefile: \-assume dummy\_aliases. (This compiler flag is _only_ needed if you are using "radopt = 1" in CM1.)

Note for users of geyser: (Posted 7 Feb 2013) If you want to use CM1 on geyser, you should compile CM1 on geyser. The login nodes of yellowstone use a different processor, and so code compiled on yellowstone will not be optimized for the geyser processors, and may not actually run at all.

* * *

**Running CM1:**

I recommend using 16 MPI processes per yellowstone node. Here is an example submission script that uses 512 total MPI processes:  

```csh
#!/bin/csh

#BSUB -P PXXXXXXXX             # project code
#BSUB -W 12:00                 # wall-clock time (hrs:mins)
#BSUB -x                       # exclusive use of node
#BSUB -n 512                   # number of tasks in job
#BSUB -R "span\[ptile=16\]"      # run 16 MPI tasks per node
#BSUB -J cm1run                # job name
#BSUB -o cm1run.out            # output file name
#BSUB -e cm1run.err            # error file name
#BSUB -q regular               # queue

mpirun.lsf ./cm1.exe >&! cm1.print.out
```

If you use hyper-threading (which on bluefire was called Simultaneous Multi-Threading, or SMT) then you can use up to 32 tasks per node; i.e., you can set ptile=32. So far, I have found little improvement (only a few percent speedup) when using hyper-threading on yellowstone with CM1. However, I haven't done much testing, so I would appreciate any feedback from any performance tests. Thanks!

* * *

**Performance of CM1:**

Regarding performance relative to bluefire: if you do not use SMT on bluefire, then you should see a significant speedup on yellowstone when using the same number of processores/cores: perhaps up to 50%. If you have been using SMT on bluefire, you should see only ~5-10% speedup for equivalent runs on yellowstone.

If you find that CM1 is running slower on yellowstone (than bluefire), please contact me.

* * *

_Last updated: 7 February 2013_
