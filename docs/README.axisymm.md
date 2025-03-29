
Brief User's Guide for running axisymmetric simulations with CM1.

Last updated:  20 March 2009.

------------------------------------------------------------------------
Setup:

To run an axisymmetric simulation with CM1, you must set ny = 1, wbc = 3,
sbc = 1, and nbc = 1.  You must also set iorigin = 1, and imove = 0.

Other limitations:

 - You can not use terrain with the axisymmetric model.
 - You must use either psolver = 2 or 3 (one of the two time-split
   compressible solvers.)
 - The only subgrid turbulence option available is iturb = 3.

All other options should work properly with the axisymmetric model. 

------------------------------------------------------------------------
How the axisymmetric solver is formulated:

In the code, "x" is radius (r):  so, "xh" is radius at the scalar points, 
and "xf" is radius at the u points.

u is radial velocity (m/s)

v is azimuthal velocity (m/s)

------------------------------------------------------------------------
Parallelization:

OpenMP (i.e., shared memory parallelization) works fine for axisymmetric
simulations.  Scaling works best if nz is divisible by the number 
of processors being used. 

MPI (i.e., distributed memory parallelization) does not work properly 
for axisymmetric simulations with CM1.  A note of caution, though:  the 
model will actually run with the axisymmetric setup with MPI, but the 
results are meaningless. 

------------------------------------------------------------------------
Efficiency:

You might ask:  Why is CM1 slower than other axisymmetric models?

There are two primary reasons why this might be the case:

    First, CM1 uses a Runge-Kutta time integration scheme, in which 
    advection is calculated three times per model timestep.  This increases
    the time spent doing advection.  However, the large timestep (dtl) can
    usually be increased, relative to the timestep used in other popular
    axisymmetric models.  For example, with CM1 I general use three times
    the timestep that I use for the Rotunno-Emanuel axisymmetric hurricane 
    model. 

    Second, CM1 is actually a three-dimensional numerical model.  The arrays 
    are not dimensioned optimally for an axisymmetric model, and there are 
    actually several redundant and unnecessary calculations in an axisymmetric 
    simulation with CM1.  I have tried to minimize the number of redundant and 
    unnecessary calculations, but there's only so much that can be done to 
    have an axisymmetric capability and yet retain the basic three-dimensional 
    design of the code. 

------------------------------------------------------------------------
Documentation:

For more information on the axisymmetric model, see:

Bryan, G. H., and R. Rotunno, 2009:  The maximum intensity of tropical 
cyclones in axisymmetric numerical model simulations.  Mon. Wea. Rev., 
137, 1770-1789, doi:10.1175/2008MWR2709.1.

------------------------------------------------------------------------

Questions, comments, suggestions:  send email to gbryan@ucar.edu

