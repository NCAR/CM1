
Brief User's Guide to using grid stretching in CM1.

Last updated:  15 June 2017

CM1 can be used with stretched grids in all three coordinates.  An
analytic stretching function (from Wilhelmson and Chen, 1982, JAS) is 
used to make it easy for model users to create smoothly varying grid
spacing.  An explanation for how to create a stretched grid is 
provided below.

Although the analytic stretching function is hard-wired into the code,
users can actually place the gridpoints wherever they wish ... i.e., 
a different stretching function can be used, or grid points can be 
arbitrarily located.  See relevant sections of code in "param.F" for 
more information.

%----------------------------------------------------------------------


Vertical stretching is covered first.  Notes on horizontal stretching follow.

########################################################################
########################################################################
#                                                                      #
#  Guide to vertical grid stretching.                                  #
#                                                                      #
########################################################################
########################################################################


First, a review of the namelist parameters:


 stretch_z - Use vertically stretched grid spacing?
               0 = no
               1 = Wilhelmson and Chen  (1982, JAS, p. 1481)
               2 = Smooth geometric  (NEW, contributed by Ted Mansell, NOAA/NSSL)

 ztop - Total depth of the domain (i.e., the height of the top of the
        domain) (m).

 str_bot - Level where stretching begins (m).

 str_top - Level where stretching ends (m). (not used for stretch_z = 2, 
           which computes this from the other inputs)

 dz_bot - Grid spacing at (and below) str_bot (m).

 dz_top - Grid spacing at (and above) str_top (m).

------------------------------------------------------------------------
STRETCH_Z = 1
Background:

The stretching function in the model is hard-wired, to ensure that the
stretching is smoothly distributed.  For those interested in learning 
more, the scheme is based on the one presented in Wilhelmson and Chen,
1982, JAS, 39, 1466-1483.

In CM1, the stretching layer may be placed arbitrarily in space.
That is, one can have a stretching layer in the middle of the domain,
near the surface, at the top, or throughout the entire domain.

For example, if you have a 20 km deep domain, and you wish to have
stretching throughout the entire domain, then specify the following:

        ztop =     20000.0
        str_bot =      0.0
        str_top =  20000.0

As another example, if you wanted stretching from the surface to 5 km,
and then constant grid spacing above 5 km, you could use the following:

        ztop =     20000.0
        str_bot =      0.0
        str_top =   5000.0

Or, say you used:

        stretch_z =  1,
        ztop      = 20000.0,
        str_bot   =  3000.0,
        str_top   =  9000.0,
        dz_bot    =   100.0,
        dz_top    =   500.0,

This would create a grid with constant grid spacing of 100 m below 3 km,
smoothly stretching vertical grid spacing from 100 m at 3 km to 500 m at
9 km, and then constant grid spacing of 500 m from 9 km to 20 km.


------------------------------------------------------------------------
The catch:

The number of vertical levels (i.e., "nz" in "input.incl") needs to be
set carefully.  In fact, the stretching function only allows for a very
specific value for "nz", given the other inputs (like dz_bot and dz_top).

***** The value for nz MUST be: *****

  nz =  ( str_bot / dz_bot )
      + ( ztop - str_top ) / dz_top 
      + ( str_top - str_bot )/( 0.5*( dz_bot + dz_top ) )

*************************************

Furthermore, this value MUST be an integer, otherwise the code will
complain.

As an example, the case listed above requires 72 levels:

        ztop      = 20000.0,
        str_bot   =  3000.0,
        str_top   =  9000.0,
        dz_bot    =   100.0,
        dz_top    =   500.0,

The number of grid points in the three layers would be:

   depth (km)    avg. dz (m)    depth (m)     grid points needed
      0-3:          100           3000              30
      3-9:          300           6000              20
      9-20:         500          11000              22
                                             -----------------
                                total =             72

Notice that the number of levels needed for the stretching layer is
the depth of the layer divided by the average grid spacing in that
layer [which is 0.5*( dz_bot + dz_top ) ].

If you are unsure what to do, just compile and run the model:  the code will 
recommend a value for nz.

------------------------------------------------------------------------

STRETCH_Z = 2
This method takes the input ztop, str_bot, dz_bot, and dz_top to compute
(iteratively) a grid stretch factor. (Method was used in COMMAS at least as
early as Skamarock et al. 2000, JGR.)  The resulting factor gives a constant
percentage change in dz from dz_bot until either dz_top is reached (then dz
remains constant until ztop) or ztop is reached (in which case stretching 
continues to the model top). To reduce numerical errors, the routine will 
abort if the stretch factor exceeds 1.1, in which case one can try 
increasing dz_bot, dz_top, and/or nz. The stretch factor and model levels 
are printed out, so some experimentation may be needed by running the model 
for one time step and checking the resulting grid.


------------------------------------------------------------------------
Other issues:

When a vertically stretched grid is used, the parameter "dz" is technically
now irrelevent.  However, due to issues such as truncation error and
efficiency, dz should still be set to a representative value.  For example,
if the grid stretches from 100 m to 500 m, an approriate value for dz would
be anything from 100 m to 500 m.  The goal should be to make the stretching
factors ("mf" and "mh" in the printout file, just to the right of the dz
values) close to 1.0 in value.

Also, since the vertical grid spacing tends to be much smaller than the
horizontal grid spacing when modelers use vertical stretching, it is
usually much more efficient to use the vertically implicit acoustic
solver (psolver=3).



########################################################################
########################################################################
#                                                                      #
#  Guide to horizontal grid stretching.                                #
#                                                                      #
########################################################################
########################################################################


NOTE:  The user does not have to use stretching in both horizontal directions.
       One can choose to have stretching in one direction only, or both.

NOTE:  The grid stretching in x is identical in application to stretching
       in y.  Thus, only grid stretching in one direction (x) is covered here.

------------------------------------------------------------------------
First, a review of the namelist parameters:


 stretch_x - Use horizontally stretched grid in x?  (0=no, 1=yes)

 dx_inner - Smallest grid spacing, in middle of domain (m).

 dx_outer - Largest grid spacing, at edge of domain (m).

 nos_x_len - Length of the no-stretching (i.e., middle) part of domain (m).

 tot_x_len - Total length of the domain (m).


------------------------------------------------------------------------
Background:

The stretching function in the model is hard-wired, to ensure that the
stretching is smoothly distributed.  It is based on the function
presented in Wilhelmson and Chen, 1982, JAS, 39, 1466-1483.

In CM1, the middle of the domain has no stretching.  This no-stretching
zone is flanked on both sides by equal-length stretching zones.

For example, if you have a 300 km long domain, and you wish to have
no stretching in the middle 100 km, then specify the following:

        tot_x_len = 300000.0
        nos_x_len = 100000.0

The grid will look like this:


     x =              x =                x =               x =
     0 km            100 km            200 km            300 km

      .     .    .   .  . . . . . . . . . .  .   .    .     .

      |   zone with     |   zone with     |   zone with     |
      |   variable      |   constant      |   variable      |
      |     dx          |     dx          |     dx          |


The parameter dx_inner specifies the constant grid spacing in the middle
zone.  The parameter dx_outer specifies the maximum grid spacing at the
edges of the domain (in this case, at x = 0 km and x = 300 km).  The
average grid spacing in the stretching zones is 0.5*( dx_outer + dx_inner ).


------------------------------------------------------------------------
The catch:

The number of grid points (i.e., "nx" in "input.incl") needs to be
set carefully.  In fact, the stretching function only allows for a very
specific value for "nx", given the other inputs (like dx_inner and
dx_outer).

***** The value for nx MUST be: *****

  nx =   ( nos_x_len / dx_inner )
       + ( tot_x_len - nos_x_len )/( 0.5*( dx_inner + dx_outer ) )

*************************************

Furthermore, this value MUST be an integer, otherwise the code will
complain.

As an example, this case requires 300 grid points:

         stretch_x =      1,
         dx_inner  =     500.0,
         dx_outer  =    3500.0,
         nos_x_len =  100000.0,
         tot_x_len =  300000.0,

The number of grid points in the three zones would be:

    zone (km)    avg. dx (m)   length (m)     grid points needed
      0-100:       2000         100000              50
    100-200:        500         100000             200
    200-300:       2000         100000              50
                                             -----------------
                      total =   300000             300

Notice that the number of points needed for the stretching zone is
the length of the zone divided by the average grid spacing in that
zone [which is 0.5*( dx_inner + dx_outer ) ].

If you are unsure what to do, just compile and run the model:  the code will
recommend a value for nx, based on the other parameters specified in
"namelist.input".


------------------------------------------------------------------------
Other issues:

When a horizontally stretched grid is used, the parameter "dx" is technically
now irrelevent.  However, due to issues such as truncation error and
efficiency, dx should still be set to a representative value.  For example,
if the grid stretches from 100 m to 500 m, an approriate value for dx would
be anything from 100 m to 500 m.  The goal should be to make the stretching
factors ("uf" and "uh" in the printout file, just to the right of the dx
values) close to 1.0 in value.



