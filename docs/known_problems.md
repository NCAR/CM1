## Known CM1 Problems and Fixes
Valid for `cm1r18` (which was released 14 August 2015)

**1. Distributed memory (MPI) runs in 2D** (first posted: 14 August 2006) (applies to ALL releases of cm1)

CM1 does not work properly for 2d runs with MPI (i.e., distributed memory). MPI should not be used for 2D simulations until further notice. \[There are no problems with single-processor runs or OpenMP (shared memory) runs. There are no problems with 3D runs with MPI.\]

Updated, 12 August 2008: This problem also applies to axisymmetric simulations. Do not run axisymmetric simulations with MPI.

* * *

**Known Problems in earlier versions of CM1**

**Known problems in `cm1r8` - `cm1r15`:** (posted 5 January 2012)

**1. There is a bug with vertically stretched grids using terrain.** This bug does not affect simulations without terrain, and it only affects simulations with terrain that use vertically stretched grids `stretch_z = 1`.

Solution: download the following files and place them in the "src" directory. The re-compile and re-run cm1.

Files for `cm1r15`:

*   [solve.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r15_fixed/solve.F)
*   [sound.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r15_fixed/sound.F)
*   [sounde.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r15_fixed/sounde.F)

Files for `cm1r14`:

*   [solve.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/solve.F)
*   [sound.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/sound.F)
*   [sounde.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/sounde.F)
*   [thompson.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/thompson.F)
*   [morrison.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/morrison.F)

Files for `cm1r13`:

*   [solve.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r13_fixed/solve.F)
*   [sound.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r13_fixed/sound.F)
*   [sounde.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r13_fixed/sounde.F)

Files for `cm1r12`:

*   [solve.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r12_fixed/solve.F)
*   [sound.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r12_fixed/sound.F)
*   [sounde.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r12_fixed/sounde.F)

To reiterate: this bug only affects runs with terrain `terrain_flag = .true.` and only when using a vertically stretched grid `stretch_z = 1`.

* * *

**Known problems in `cm1r14`:** (posted 12 January 2011)

**1\. Thompson and Morrison microphysics (only if `neweqts={1,2}`):** There is a bug with the implementation of the mass- and energy-conserving equations for the Thompson (ptype=3) and Morrison (ptype=5) microphysics schemes in `cm1r14`. There is no problem if neweqts=0 ("traditional" equation set). This bug _only_ impacts the Thompson and Morrison microphysics schemes (ptype = 3 and 5, respectively).

Solution: download the following files and place them in the "src" directory. Then re-compile and re-run. (For `cm1r14` only.)

*   [solve.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/solve.F)
*   [thompson.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/thompson.F)
*   [morrison.F](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1r14_fixed/morrison.F)

Explanation: there was an "order-of-operations" error, wherein the old pressure and the new temperature were used to calculate the new value of potential temperature. This bug caused cold pools to be too strong (by ~15% at the lowest model level) in idealized simulations of midlatitude squall lines. The magnitude of errors in other phenomena (e.g., supercells, hurricanes) is unknown at this time. Users of `cm1r14` are encouraged to download the modified code (above) and re-run their simulations.

To reiterate: this bug only affects `cm1r14`, and it only affects the Thompson and Morrison microphysics schemes if `neweqts=1` or `neweqts=2`.

* * *
_Last updated: 6 February 2012_
