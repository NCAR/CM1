### Sample submission scripts for NCAR's <a href="https://arc.ucar.edu/knowledge_base/74317833">derecho</a>:

- [cm1run_derecho](cm1run_derecho):  multiple CPUs using MPI (distributed memory) (derecho)

- [run_case](run_case): Script to run a specific CM1 simulation case.

```csh
Usage:
   run_case <case> [-f]
     <case>: Name of the simulation case (e.g., "dunion_MT").
     -f    : (Optional) Force overwrite of existing run directories.
```