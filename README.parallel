
Brief User's Guide for running CM1 with multiple processors

Last updated:  10 January 2011

------------------------------------------------------------------------
Background:

CM1 has been designed to use multiple processors.  There are two 
conventions for doing this:  

One is shared memory parallelization, in which every processor accesses 
the same bank of memory.  CM1 uses OpenMP directives to parallelize
the code on these architectures.  This is relatively efficient for a few
number of processors, but probably does not scale well beyond 10 processors.

Second is distributed memory parallelization, in which each processor has 
it's own piece of the domain;  only a small amount of information is shared
between processors, and this information must be explicitly passed between 
processors, when needed.  CM1 uses MPI to handle this computing paradigm.

Before choosing your configuration, obviously you want to first consider 
the hardware you will be using.  If hundreds of processors are required, 
then choose MPI.  If 2-10 processors will be used on a shared memory 
machine, then OpenMP would be a good choice.  Then, the table below will
let you know what model configurations are available for your machine. 
If something is not listed below, then you may assume that it works 
properly in both OpenMP and MPI.  (the only exceptions to this last 
statement are 2D and axisymmetric simulations:  see README.2D and 
README.axisymm for more information)

Certain model configurations are not currently parallelized.  The specific 
reason varies for each configuration, but may be related to difficulty 
achieving good parallel efficiency (e.g., MPI with psolver=4,5), or
because of lack of time/interest/need by the model developer (e.g., 
OpenMP with ptype=4). 


------------------------------------------------------------------------

                    X    = works correctly
                    -    = does not work

                      OpenMP        MPI

Psolvers:

       psolver=1:       X            X
       psolver=2:       X            X
       psolver=3:       X            X
       psolver=4:       X            -
       psolver=5:       X            -

Microphysics:

       ptype=1:         X            X
       ptype=2:         X            X
       ptype=3:         X            X
       ptype=4:         -            X
       ptype=5:         X            X
       ptype=6:         X            X

------------------------------------------------------------------------

Questions, comments, suggestions:  send email to gbryan@ucar.edu


