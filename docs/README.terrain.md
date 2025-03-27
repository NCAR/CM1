
Brief User's Guide to using terrain in CM1.

Last updated:  13 January 2011


########################################################################
########################################################################
#                                                                      #
#  Guide to terrain in CM1.                                            #
#                                                                      #
########################################################################
########################################################################


Background:

The model uses the Gal-Chen and Somerville coordinate transform.  The 
model surfaces follow the terrain at the bottom of the model, and the
model top is a surface of constant height.  The model surfaces do not
move in time ... i.e., they are fixed in space.

When there is no grid stretching, the model top is simply nz x dz.
(Both variables are specified in namelist.input)

When vertical grid stretching is used, the model top iz ztop.  (ztop is
specified in the param6 section of namelist.input)

The terrain profile is specified by the array called "zs".  The user
should set zs in init_terrain.F



References:

Gal-Chen, T., and R. Somerville, 1975:  On the use of a coordinate
transformation for the solution of the Navier-Stokes equations.
J. Comput. Phys., 17, 209-228.

Durran, D. R., and J. B. Klemp, 1983:  A compressible model for the
simulation of moist mountain waves.  Mon. Wea. Rev., 111, 2341-2361.



!-----------------------------------------------------------------------


First, a brief overview:

 step 1:  Set "terrain_flag = .true." in the param0 section of 
          "namelist.input"

 step 2:  Set a value for "itern" in the param2 section of "namelist.input"
          This tells the model what analytic terrain profile to use.
          See details in "README.namelist" for the list of pre-defined
          profiles.

 step 3:  Specify the terrain profile (zs) in "init_terrain.F" (or, just 
          double check the settings if you are using one of the pre-defined
          profiles).  The zs array is the height of the topography in m.
          NOTE:  zs is specified at the scalar points (i.e., the same as
          theta, qv, and pressure).

 step 4:  If you feel so inclined, set "output_interp = 1" in the
          file "namelist.output".  This will generate an additional
          set of output files, where the data have been interpolated to
          the nominal height levels (i.e., the height levels in the
          absence of topography ... meaning if zs = 0).  NOTE:  This 
          doubles the amount of output.  ALSO NOTE:  This is a rather
          quick-and-dirty method of looking at the output.  A more
          sophisticated plotting scheme/software should be preferred
          for journal-quality figures.

 step 5:  Compile.

 step 6:  Run model.

 step 7:  View the output.  You're on your own here.  The output in the
          regular output files (cm1out_s, cm1out_u, etc.) is on the 
          terrain-following model surfaces.  The output in the interp
          file (cm1out_i) can be rather coarse in the vertical.
          (I'm working on something better ... when I get the chance.)


!-----------------------------------------------------------------------

Questions:  send email to gbryan@ucar.edu


