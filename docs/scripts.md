### Sample submission script for NCAR's <a href="https://arc.ucar.edu/knowledge_base/74317833">derecho</a>:

- `utils/cm1run_derecho`:  multiple CPUs using MPI (distributed memory) (derecho)

### Run a CM1 test case

- `utils/run_case`: csh script to run a CM1 test case

```csh
Usage:
   utils/run_case <case> [-f]
     <case>: Name of the test case (e.g., "dunion_MT").
     -f    : (Optional) Force overwrite of existing run directories.
```
