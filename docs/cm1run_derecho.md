```bash
#!/bin/bash

#
#  Example submission script for CM1 (using MPI) on NCAR's derecho
#
#  Last updated:  24 March 2024
#
#  Notes:
#    - In the CM1 Makefile, we recommend the section titled:
#      "multiple processors, distributed memory (MPI), Intel compiler"
#    - If present, remove the flag "-xHost" in the Makefile before compiling.
#    - The default modules on derecho are recommended for CM1.  In our testing,
#      they are:
# cm1r21.0/src> module list
# Currently Loaded Modules:
#  1) ncarenv/23.09 (S)   3) intel/2023.2.1        5) cray-mpich/8.1.27   7) netcdf/4.9.2
#  2) craype/2.7.23       4) ncarcompilers/1.0.0   6) hdf5/1.12.2
#
#    - Also, in namelist.input we recommend using ppnode = 128
#


#-------------------------------------------
# PBS stuff below this line
# (change things here)


# job name:
#PBS -N cm1run


# project code:
# (you must specify the project code to be charged when you submit a job)
#PBS -A xxxxxxxx


# below here, "select" is the number of 128-CPU nodes to use.
# note: this example uses 512 (=4*128) CPUs:
# (do not change settings for "ncpus" or "mpiprocs" or "ompthreads")
# For more info, see: https://arc.ucar.edu/knowledge_base/74317833
#
#PBS -l select=4:ncpus=128:mpiprocs=128:ompthreads=1


# maximum wall-clock time (hh:mm:ss)
#PBS -l walltime=12:00:00 


# queue:
#PBS -q main


#-------------------------------------------
# command-line stuff below this line
# (probably should not change)

# temp directory:
export TMPDIR=/glade/derecho/scratch/$USER/temp
mkdir -p $TMPDIR

# These seem to work well for CM1 runs
export PALS_PPN=128
export PALS_DEPTH=1
export PALS_CPU_BIND=depth

# run CM1
mpiexec --cpu-bind depth ./cm1.exe >& cm1.print.out
```
