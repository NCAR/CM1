
Brief User's Guide to using parcels in CM1.

Last updated:  11 October 2006.

------------------------------------------------------------------------
Background:

In cm1, parcels can be translated every time step during the model
integration.  The parcels do not interact with the simulated flow;  they
are merely passive parcels that are used to diagnose conditions in a 
parcel-relative perspective, during the simulation.  In addition, the 
code can output environmental conditions experienced by the parcel along
its trajectory.  A few examples are included, by default, but the user
can output just about anything, with a little bit of coding.  

------------------------------------------------------------------------
How to initialize the parcels:

In init3d.F, look for the section of code where the "pdata" array is 
defined.  By default, pdata(1,n) is x, pdata(2,n) is y, and pdata(3,n) 
is z (all in meters).  In the default example, 36000 parcels are 
initialized at the beginning of the model simulations:  one every 2 km
in x, one every 2 km in y, and one every 1 km in z.  The user can have 
as few as 1 parcel, or could have millions.   

It is easiest to initalize the parcels at the model start time.  However, 
to initialize parcels in the middle of a model integration, the user 
will have to add some code to let the model know this.  For example, it 
would be relatively straightforward to add an "if" statement around the 
call the "parcel_driver" in the "solve" subroutine, such that 
parcel_driver is only called after 1 h into the simulation. 

------------------------------------------------------------------------
Information along parcels:

In principle, just about any conditions along the parcel trajectory can 
be calculated and printed to file.  By default, several are included, 
such as qv, qc, qr, equivalent potential temperature, and buoyancy.  If
other variables are desired, instead of these, then users will have to 
modify the code.  It should be fairly straightforward to modify one of 
the existing variables (e.g., by replacing kh with qg, for example).

------------------------------------------------------------------------
Output:

All data is dumped to the basename_pdata.dat file, a GrADS binary file. 
Each "x" value in this file is a different parcel;  that is, "x=1" is 
parcel number 1, "x=232" is parcel number 232, etc.  With cm1r11, there 
is a new variable called "prclfrq" in "namelist.input".  This tells the 
code how often to output parcel information.  If you would like parcel 
information every time step, then set prclfrq = dtl.  However, the 
subsequent data file can become VERY large, especially if a large number 
of parcels are used.  In my experience, it is sufficient to output parcel 
information every 30-60 seconds, depending on the application.  This will 
significantly reduce the size of the output file, but will not change the 
accuracy of the calculation, because the parcels will still be advanced 
every timestep.  


------------------------------------------------------------------------

Questions, comments, suggestions:  send email to gbryan@ucar.edu


