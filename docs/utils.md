## Fortran programs that combine the tiled GrADS output from CM1 into one file: valid for `cm1r12` and later versions

(used for processing the output from CM1 when using MPI processing and when the GrADS output format is selected):

Note: these new, more efficient "combine" programs have been released with cm1r16:

*   `combine_MPI_s.F`: combines scalar (\*\_s) files
*   `combine_MPI_u.F`: combines u files
*   `combine_MPI_v.F`: combines v files
*   `combine_MPI_w.F`: combines w files

Usage: compile with a fortran compiler. Then, run the program, and answer all the questions. Note: you will need to have a copy of your GrADS descriptor file (i.e., files ending in ".ctl") in the same directory.

Output: two files are created: a new descriptor file (cm1out\_MPI.ctl) and a new GrADS data file.

## Fortran program to combine the tiled netcdf output
(same as above, but for netcdf output format)

*   `netcdf/combine.F`: fortran program (`cm1r11` ... contributed by Dan Kirshbaum, U. Reading)
*   `netcdf/combine_r12.F`: fortran program (`cm1r12` and later versions ... contributed by Lou Wicker, NSSL)
*   `netcdf/compile`: script used to compile the program
