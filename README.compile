
Brief User's Guide for compiling CM1.

Last updated:  14 August 2015

------------------------------------------------------------------------

How to compile cm1:

1) cd into the "src" directory  (type:  "cd src")

2) edit "Makefile"

   - Choose the hardware that you will be running on
     - if using a single processor, choose a "single memory" option
     - if using multiple processors with distributed memory parallelization
       (i.e., if you want to use MPI), choose a "distributed memory" option
     - if using multiple processors with OpenMP parallelization, choose a 
       "shared memory" option

   - Once a suitable hardware configuration has been chosen, uncomment all
     lines in that section

   - If necessary, change options in that section for your particular
     architecture.  You may need to add paths (e.g., "/lib/cpp" instead 
     of simply "cpp"), and you may want to change compiler flags.

   - If you want to include the capability to write output in netcdf 
     format, then uncomment the appropriate lines near the top of the 
     Makefile.  (See comments at top of Makefile for more information.)
     Make sure the paths to your netcdf files are correct.

3) type "make" to compile the model

   - NOTE:  on some machines, you may need to use "gmake" instead of "make".

4) if successful, two files will have been created in the "run" directory:
   cm1.exe (the executable) and onefile.F (an archive of the code used to
   make this particular executable ... I highly recommend that you retain
   onefile.F with every simulation, because it makes a nice record of the 
   model code for your particular simulation, and it sometimes helps me to
   debug potential problems)

5) to cleanup the src directory, type "make clean".  You should also do this 
   when changing compiler flags.


Questions, comments, suggestions:  send email to gbryan@ucar.edu


