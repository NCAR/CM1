
Brief User's Guide to restarts with cm1.

Last updated:  20 April 2022 
(see new section "Changing resolution on a restart" below)

-------------------------------------------------------------------------
Background:

Restart files allow the model simulation to "start again" (i.e., restart)
from a previous point during the model simulation.  The model automatically 
saves all information that is carried from one timestep to the next in the 
cm1rst* files.  

Please note that major changes on a restart ... like changing the 
microphysics or PBL parameterization ... may not work, and could cause 
very bad things to happen.  However, minor changes in model configuration, 
like the advection scheme or time step, could probably be changed without 
any problems.  Also, in principle, users could add parcels and passive fluid 
tracers on a restart.  

Bottom line: if you are relying on restart files to do something other 
than simply continuing a simulation from before ... as it was configured 
before ... then you may have to "get your hands dirty" and modify code in 
the restart_read.F and/or restart_write.F files.  

-------------------------------------------------------------------------
To save restart files:

Simply set a value for "rstfrq" in "namelist.input".  The model will 
generate a series of restart files, named (for example)
cm1rst_000001_*.dat.  The set of numbers (000001 in the case) is the restart 
file number:  it starts at 1, and counts upward;  use this number to start a 
subsequent restart simulation (see next section).

-------------------------------------------------------------------------
To restart the model:

In "namelist.input", set irst = 1, and then specify the restart file 
number (rstnum).  For example, if a user wanted to use the file set 
cm1rst_000004_*.dat, then specify rstnum = 4 in namelist.input.  

WARNING: restarts from CM1 may overwrite files from the previous simulation. 
For example, if you first ran CM1 for 2 days, and then restarted from restart
files written after 1 day, then all files from 1-2 days from the previous 
simulation will be overwritten.  To avoid this potential problem, it is 
usually best to move all output files from a previous simulation to a 
separate subdirectory to make sure they are not lost. 

-------------------------------------------------------------------------
Changing resolution on a restart:

Starting with cm1r21.0, users can change the resolution on a restart. 
This process is started automatically when the code recognizes that the 
current grid is different from the grid saved in the restart file. 

Note: the new domain must be equal in size, or smaller than, the domain from 
the previous simulation.  In other words, you cannot make the domain larger 
in any direction on a restart. 

Note: this option only works for restart_format=1 (binary-format output files)
at the moment. 

Note: this option is new and relatively untested, so please send comments 
and questions to the email address below if you encounter problems.  

-------------------------------------------------------------------------

Questions, comments, suggestions:  send email to gbryan@ucar.edu

