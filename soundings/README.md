## Soundings for idealized simulations:

- [input_sounding_trier](input_sounding_trier):  TOGA COARE squall line (Trier et al, 1996, JAS, p. 2861)
- [input_sounding_jordan_allmean](input_sounding_jordan_allmean):  West Indies annual mean (Jordan, 1958, J. Meteor., p. 91)
- [input_sounding_jordan_hurricane](input_sounding_jordan_hurricane):  West Indies mean hurricane season (Jordan, 1958, J. Meteor., p. 91)
- [input_sounding_rotunno_emanuel](input_sounding_rotunno_emanuel):  approximately moist-neutral hurricane sounding (Rotunno and Emanuel, 1987, JAS, p. 542)
- [input_sounding_dunion_MT](input_sounding_dunion_MT):  "Moist Tropical" (MT) sounding from Dunion (2011, J Clim, p. 893)
- [input_sounding_bryan_morrison](input_sounding_bryan_morrison):  VORTEX2 squall line (Bryan and Morrison, 2012, MWR)
- [input_sounding_seabreeze_test](input_sounding_seabreeze_test):  sounding for sea breeze test case

* * *

### Description of input_sounding files:
  The format is the same as that for the WRF Model. 

<pre>
  One-line header containing:   sfc pres (mb)    sfc theta (K)    sfc qv (g/kg)

   (Note1: here, "sfc" refers to near-surface atmospheric conditions. 
    Technically, this should be z = 0, but in practice is obtained from the 
    standard reporting height of 2 m AGL/ASL from observations)
   (Note2: land-surface temperature and/or sea-surface temperature (SST) are 
    specified elsewhere: see tsk0 in namelist.input and/or tsk array in 
    init_surface.F)

 Then, the following lines are:   z (m)    theta (K)   qv (g/kg)    u (m/s)    v (m/s)

   (Note3: # of levels is arbitrary)

     Index:   sfc    =  surface (technically z=0, but typically from 2 m AGL/ASL obs)
              z      =  height AGL/ASL
              pres   =  pressure
              theta  =  potential temperature
              qv     =  mixing ratio
              u      =  west-east component of velocity
              v      =  south-north component of velocity

 Note4:  For final line of input_sounding file, z (m) must be greater than the model top 
         (which is nz * dz when stretch_z=0, or ztop when stretch_z=1,  etc)

</pre>
_Last updated:  16 November 2017_
