## Parallel performance of CM1:

This page presents information about the performance of CM1 on distributed memory supercomputers.

Contents: | [A Strong Scaling Test](#strong) | [A Weak Scaling Test](#weak) |

* * *

<a name="strong"></a>**A Strong Scaling Test with cm1r19 on NSF NCAR's _cheyenne_**: Supercell thunderstorm simulation (posted 27 June 2017)

**System**: NSF NCAR's [cheyenne](https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne): SGI ICE XA Cluster with Intel Broadwell processors.  
**Compiler**: ifort 16.0.3  
**Code**: cm1r19.1  
**CM1 Configuration**: MPI  
**Case**: Idealized supercell thunderstorm, 2 h integration, 250 m horizontal grid spacing  
**Total domain dimensions**: 576 × 576 × 128  
**Time steps**: 2,880  
**NOTE**: This is a **STRONG SCALING** test (i.e., problem size is fixed) that includes moisture as well as Input/Output  

**Results with I/O:** (full 3d output every 15 min; 8 output times total; 28 GB total)  
  
![ch_r19.png](https://www2.mmm.ucar.edu/people/bryan/cm1/ch_r19.png)  

Comments: This test demonstrates that the CM1 solver scales reasonably well for ~10,000 cores (black line above). However, when using more than 1,000 cores, the time required to write output begins to affect parallel performance. For binary GrADS-format output (output\_format = 1) (red line) parallel performance degrades beyond roughly 4,000 cores. For netcdf output ("output\_format = 2) (blue line) parallel performance does not scale well beyound roughly 1,000 cores.

Results will vary depending on the frequency of output, the total size of output, and (especially) based on model physical schemes. For example, simulations with more expensive microphysics schemes and simulations with atmospheric radiation will require more CPU time to complete.

Recommendation: For this configuration (that is, with the Morrison microphysics scheme, the LES subgrid turbulence model, and no radiation scheme), a good formula for estimating core-hours required for a simulation on NSF NCAR's cheyenne supercomputer is:  
![Eqn1.png](https://www2.mmm.ucar.edu/people/bryan/cm1/Eqn1.png)  
where _C_ is the total number of cheyenne core-hours, _Nx_ is the number of grid points in the _x_ direction, _Ny_ is the number of grid points in the _y_ direction, _Nz_ is the number of grid points in the _z_ direction, and _Nt_ is the total number of timesteps.

For example: a 512 × 512 × 128 domain, integrated for 2 hours with a 2.5-s timestep (and thus 2,880 total time steps), would require approximately 242 cheyenne core-hours. So, assuming 144 cores (i.e., 4 nodes) are used, then roughly 1.7 wallclock hours would be needed to run this simulation on cheyenne.

Also: for large simulations that require large processor counts, users of CM1 should use binary GrADS-format output (output\_format = 1) to ensure acceptable parallel performance. If netcdf-format output is required, then conversion software can be used after the CM1 simulation is complete. For example, the Climate Data Operators (CDO) package can be used to convert GrADS-format data to netcdf format. For example, on cheyenne, type:

```
module load nco
```

followed by a command of form:

```
cdo -f nc4 import_binary cm1out_s.ctl cm1out_s.nc
```

* * *

<a name="weak"></a>**A Weak Scaling Test with cm1r19 on NSF NCAR's _cheyenne_**: Large Eddy Simulation of a convective boundary layer (posted 29 June 2017)

**System**: NSF NCAR's [cheyenne](https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne): SGI ICE XA Cluster with Intel Broadwell processors.  
**Compiler**: ifort 16.0.3  
**Code**: cm1r19.1  
**CM1 Configuration**: MPI  
**Case**: LES of convective boundary layer, 40 m horizontal grid spacing  
**Domain dimensions**: varies with number of cores (i.e., processors). Each core has 16 × 16 × 128 grid points. The largest total domain size is 3,072 × 3,072 × 128 grid points.  
**Time steps**: 3,600  
**NOTE**: This is a **WEAK SCALING** test (i.e., problem size scales with number of processors) that does not include moisture/microphysics but _does_ include Input/Output  

**Results with I/O:** (full 3d output every 15 min; 9 output times total; 18 MB per core)  
  
![ch_r19_weak](https://www2.mmm.ucar.edu/people/bryan/cm1/ch_r19_weak.png)

Comments: For this test the domain size increases as cores are added, with the goal of testing how CM1 performs for very large domains. In this case, the total run time should (ideally) remain the same as cores are added.

Results in the figure above show that the CM1 solver performs well to at least 36,000 cores (black line). For more that roughly 10,000 cores, however, the time required to write output begins to degrade parallel performance. For binary GrADS-format output (red line) parallel performance is not seriously impacted until roughly 20,000 cores. For the netcdf output (blue line), however, parallel performance is negatively impacted with only a few thousand cores.

Recommendation: For large datasets that require large processor counts, users of CM1 should use binary GrADS-format output (output\_format = 1) to ensure acceptable parallel performance. If netcdf-format output is required, then conversion software can be used after the CM1 simulation is complete. For example, the Climate Data Operators (CDO) package can be used to convert GrADS-format data to netcdf format. For example, on cheyenne, type:

```
module load nco
```

followed by a command of form:

```
cdo -f nc4 import_binary cm1out_s.ctl cm1out_s.nc
```

* * *

_Last updated: 29 June 2017_
