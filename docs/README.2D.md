
Brief User's Guide for running two-dimensional simulations with CM1.

Last updated:  12 August 2008.

------------------------------------------------------------------------
Setup:

To run a two-dimensional (2D) simulation with x and z as the two dimensions, 
set ny = 1, sbc = 1, and nbc = 1.

To run a two-dimensional (2D) simulation with y and z as the two dimensions,
set nx = 1, wbc = 1, and ebc = 1.

NOTE:  A 2D simulation with x and z as the two dimensions will run much 
faster than an equivalent 2D simulation with y and z as the two dimensions. 
So, if you have a choice, use the x-z option. 

I have never tested a 2D simulation with x and y as the only two dimensions,
nor do I expect this to work, so you probably shouldn't try it.

------------------------------------------------------------------------
Parallelization:

OpenMP (i.e., shared memory parallelization) will work fine for 2D simulations.

MPI (i.e., distributed memory parallelization) is currently not working 
correctly for 2d simulations in cm1.  A note of caution, though:  the model 
will actually run in 2D with MPI, but the results are meaningless.
