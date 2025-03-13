## Brief Users' Guide to CM1

### A quick summary of how to download, compile, and run CM1

Valid for `cm1r18` and later versions

> [!NOTE]
> Knowledge of UNIX and FORTRAN is *required* to use CM1. Experience with numerical modeling is very beneficial.

* * *

#### Step 1: Download the code.

1.  Read through the [license](/LICENSE).
2.  If you agree to the terms in the license, they you may [download and use CM1](/releases).
3.  CM1 is distributed as a gzipped tar file that can be downloaded from this page.
4.  Place the file on the disk that you will be running the model. Since the output from CM1 can be large, you should place the file on a disk with a substantial amount of free space (i.e., several GB of disk space should be available).
5.  To uncompress and extract the files, type: `gunzip cm1r18.tar.gz` followed by `tar xvf cm1r18.tar`
6.  Proceed to Step 2 below.

* * *

#### Step 2: Setup the code for your study.

1.  **Edit `Makefile`:** cd into the `src` directory. In `Makefile`, select the operating system and parallelization method that is appropriate for your hardware, and uncomment all lines in that section. Compiler flags can also be set/changed in this file. If you want to use netcdf or hdf5 output, uncomment the appropriate lines at the very top of `Makefile` and set the paths to your netcdf/hdf5 distribution.
2.  **Edit `base.F` (if necessary):** Modify the base-state conditions, as appropriate. There are two sections: one for the hydrostatic pressure, temperature, and moisture sounding (see "isnd" section); and one for the initial winds (u and v components) (see "iwnd" section).
3.  **Edit `init3d.F` (if necessary):** In `init3d.F,` you can add perturbations to the base state. Several default options are available.
4.  **Edit `init_terrain.F` (if necessary):** If you are using terrain, you will have to specify the terrain via the "zs" array in the file `init_terrain.F`.
5.  **Edit `init_surface.F` (if necessary):** If you are using surface fluxes of heat/moisture/momentum, then you might have to specify the horizontal distribution of several variables in the file `init_surface.F`. See the param12 section in [`README.namelist`](/README.namelist) for more information.
6.  **Compile the code:** Type `make` within the "src" directory. On some machines, you may need to use "gmake" instead.
7.  **Edit `namelist.input`:** If the code compiled without error, cd into the "run" directory and edit "namelist.input". See [`README.namelist`](/README.namelist) for guidance. Here, you set the domain dimensions, as well as the number of processors (using "nodex" and "nodey"). See the README files in the main directory for more information.
8.  **Place the `input_sounding` file in the same directory as `cm1.exe` (if necessary):** If you are supplying an external sounding file, make sure it is called `input_sounding` and place it in the same directory as `cm1.exe`. See ["Soundings for idealized simulations"](/soundings) for more information.
9.  **Place the `LANDUSE.TBL` file in the same directory as `cm1.exe` (if necessary):** If you are using surface fluxes of heat/momentum/moisture, or if you are using the atmospheric radiation scheme, then you need to specify the surface conditions. (See the param12 section of [`README.namelist`](/README.namelist) for more information.) The `LANDUSE.TBL` file comes with cm1 in the `run` directory.

* * *

#### Step 3: Run CM1.

1.  **Run the code:** Type `./cm1.exe` to run the code. You will probably need a different command for some MPI applications ... check the documentation for your supercomputer for more information.
2.  **stout and sterr info:** To redirect all standard output and standard error information ... i.e., all the stuff that prints to screen when running cm1 ... use this command: `./cm1.exe >&! cm1.print.out &`
3.  **Output:** All output is placed in the "run" directory, by default. You can move the output files to another directory when the code is finished.

  > [!CAUTION]
  > Running the code again will overwrite the files if they are not moved.

* * *

#### For more information ...

For more information about model settings, see [`README.namelist`](/README.namelist).

For more information about running on distributed memory computers (with MPI), see [`README.parallel`](/README.parallel).

For more information about running with terrain, see [`README.terrain`](/README.terrain).

For more information about running the axisymmetric version of the model, see [`README.axisymm`](/README.axisymm).

For more information about running the model with stretched grids, see [`README.stretch`](/README.stretch).

* * *

_Last updated: 13 August 2015_
