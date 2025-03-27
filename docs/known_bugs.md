## CM1: Known Bugs

Valid for cm1r14 (released 23 October 2009)

* * *

**1\. Distributed memory (MPI) runs in 2D** (first posted: 14 August 2006)

CM1 does not work properly for MPI runs in 2D. MPI should be avoided for 2D simulations until further notice. (There are no problems with single-processor runs or OpenMP (shared memory) runs. There are no problems with 3D runs in MPI.)

Updated, 12 August 2008: This problem also applies to axisymmetric simulations. Do not run axisymmetric simulations with MPI.

* * *

_Last updated: 23 October 2009_
