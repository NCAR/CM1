## Sample submission script for NCAR's [derecho](https://arc.ucar.edu/knowledge_base/74317833)
   * multiple CPUs using MPI (distributed memory)
   * [cm1run_derecho](cm1run_derecho.md)


## Run a CM1 test case

- `utils/run_case`: csh script to run a CM1 test case

```csh
Usage:
   utils/run_case <case> [-f]
     <case>: Name of the test case (e.g., "dunion_MT").
     -f    : (Optional) Force overwrite of existing run directories.
```
## Convert GrADS format to netcdf

- `output_format = 1` or unformatted direct-access binary format
- Requires CDO (climate data operators)
- On NCAR's derecho, use the two following commands:

```
module load cdo
cdo -f nc4 import_binary cm1out_s.ctl cm1out_s.nc
```
- see https://code.mpimet.mpg.de/projects/cdo for more info on CDO

## [Combining output (MPI users)](combine_output.md)

- Fortran programs that combine the tiled GrADS output from CM1 into one file

## [GrADS Scripts](grads.md)

## Links

- [Grid Analysis and Display System (GrADS)](http://cola.gmu.edu/grads)
- [NetCDF](http://www.unidata.ucar.edu/software/netcdf)
- [Visualization and Analysis Platform for Ocean, Atmosphere, and Solar Researchers (VAPOR)](http://www.vapor.ucar.edu)
