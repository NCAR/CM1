## Brief User's Guide for compiling CM1.

------------------------------------------------------------------------

### How to compile cm1:

1) cd into the "src" directory  (type:  "cd src")

2) To compile the CM1 codes with this Makefile, use the following command:

```
  make [argu_list]

  where the optional argument list includes:
  - USE_OPENMP: turn on OpenMP (default=false)
  - USE_MPI: turn on MPI (default=false)
  - USE_DOUBLE: turn on double precision calculation (for PGI/NVHPC only, default=false)
  - USE_OPENACC: turn on OpenACC (for PGI/NVHPC only, default=false)
  - USE_NETCDF: turn on netCDF output (default=false)
  - DEBUG: turn on DEBUG mode (default=false)

  Some example usages are:
    make USE_OPENMP=true                 # shared memory only
    make USE_OPENMP=true USE_MPI=true    # shared + distributed memory
    make USE_OPENACC=true USE_MPI=true   # distributed memory + GPU offloading

  Note that the logical values set to the argument variables are case insensitive.

  If the "FC" env variable is not set on your system, you could also specify it through:
    make FC=ifort                        # use Intel compiler
```

   - NOTE:  on some machines, you may need to use "gmake" instead of "make".

3) if successful, two files will have been created in the "run" directory:
   cm1.exe (the executable) and onefile.F (an archive of the code used to
   make this particular executable ... I highly recommend that you retain
   onefile.F with every simulation, because it makes a nice record of the 
   model code for your particular simulation, and it sometimes helps me to
   debug potential problems)

4) to cleanup the src directory, type "make clean".  You should also do this 
   when changing compiler flags.
