Some background information and usage tips for the NSSL microphysics scheme.


 IMPORTANT: Best results are attained using WENO (Weighted Essentially Non-Oscillatory) scalar advection option. This helps to limit oscillations at the edges of precipitation regions (i.e., sharp gradient), which in turns helps to prevent mismatches of moments that can show up as noisy reflectivity values.
 advwenos = 2
 weno_order = 5 or higher

NOTE TO OpenMP USERS: If a segmentation fault occurs, try setting the environment variable OMP_STACKSIZE to 8M or 16M (default is 4M, where M=MB). Note that this does not increase the shell stacksize limit [use 'ulimit -a unlimited' (bash) or 'unlimit stacksize' (tcsh)]

CHANGES:
May 2023: Main default option change is for graupel/hail fall speed options (icdx, icdxhl; changed from 3 to 6, see below), and default maximum gr/hail droplet collection efficiencies (ehw0/ehlw0 changed from 0.5/0.75 to 0.9/0.9, see below)

DESCRIPTION:

The NSSL bulk microphysical parameterization scheme describes form and phase changes among a range of liquid and ice hydrometeors, as described in Mansell et al. (2010) and Mansell and Ziegler (2013). It is designed with deep (severe) convection in mind at grid spacings of up to 4 km, but can also be run at larger grid spacing as needed for nesting etc. It is also able to capture non-severe and winter weather. The scheme predicts the mass mixing ratio and number concentration of cloud droplets, raindrops, cloud ice crystals (columns), snow particles (including large crystals and aggregates), graupel, and (optionally) hail. The 3-moment option additionally predicts the 6th moments of rain, graupel, and hail which in turn predicts the PSD shape parameters (set nssl_3moment=.true.).

Basic options in param2 namelist:
 ptype = 26 ! NSSL scheme (2-moment) without hail
 ptype = 27 ! NSSL scheme (2-moment) with hail
 nssl_3moment : default .false., setting true adds 6th moment for rain, 
                graupel (i.e., 3-moment )(ptype=26/27)  and hail (ptype=27)
 nssl_density_on : default .true.; Setting to false turns off graupel/hail predicted
                   ice density and instead uses fixed (constant) ice density 
                  for graupel (rho_qh, default 500.) and hail (rho_qhl, default 800.)

 Other options are in the nssl2mom_params namelist (see also the README.namelist)

The graupel and hail particle densities are also calculated by predicting the total particle volume. The graupel category therefore emulates a range of characteristics from high-density frozen drops (includes small hail) to low-density graupel (from rimed ice crystals/snow) in its size and density spectrum. The hail category is designed to simulate larger hail sizes. Hail is only produced from higher-density large graupel that is actively riming (esp. in wet growth).

Hydrometeor size distributions are assumed to follow a gamma functional form. Microphysical processes include cloud droplet and cloud ice nucleation, condensation, deposition, evaporation, sublimation, collection–coalescence, variable-density riming, shedding, ice multiplication, cloud ice aggregation, freezing and melting, and conversions between hydrometeor categories. 

Cloud concentration nuclei (CCN) concentration is predicted as in Mansell et al. (2010)  with a bulk activation spectrum approximating small aerosols. The model tracks the number of unactivated CCN, and the local CCN concentration is depleted as droplets are activated, either at cloud base or in cloud. The CCN are subjected to advection and subgrid turbulent mixing but have no other interactions with hydrometeors; for example, scavenging by raindrops is omitted. CCN are restored by droplet evaporation and by a gradual regeneration when no hydrometeors are present (ccntimeconst). Aerosol sensitivity is enhanced by explicitly treating droplet condensation instead of using a saturation adjustment. Supersaturation (within reason) is allowed to persist in updraft with low droplet concentration.

  ccn          - Initial concentration of cloud condensation nuclei
                   0.25e+9 maritime
                   0.6e+9 "low" continental (DEFAULT)
                   1.0e+9 "med-high" continental
                   1.5e+09 - high-extreme continental CCN)
                   Larger values run a risk of unrealistically low precipitation production
                 Value sets the concentration at MSL, and an initially
                 homogeneous number mixing ratio (ccn/1.225) is assumed throughout the depth of
                 the domain.

Droplet activation option is controlled by the 'irenuc' option. Old option (2) depletes CCN from unactivated CCN field. New option (7) instead counts the number of activated CCN (nucleated droplets) with the assumption of an initial constant CCN number mixing ratio. Option 7 better handles supersaturation at low CCN (e.g., maritime) concentrations by allowing extra droplet activation at high SS.

  irenuc (new option) : 2 = ccn field is UNactivated aerosol (previous default; old droplet activation)
                        7 = ccn field is ACTVIATED aerosol (new default 2023) (new droplet activation)

Excessive size sorting (common in 2-moment schemes) is effectively controlled by an adaptive breakup method that prevents reflectivity growth by sedimentation (Mansell 2010). For 2-moment, infall=4 (default) is recommended. For 3-moment, infall only really applies to droplets, cloud ice, and snow.

Graupel -> hail conversion: The parameter ihlcnh selects the method of converting graupel (hail embryos) to the hail category. The default value is -1 for automatic setting. The original option (ihlcnh=1) is replaced by a new option (ihlcnh=3) as of May 2023. ihlcnh=3 converts from the graupel spectrum itself based on the wet growth diameter, which generally results in fewer initiated hailstones with larger diameters (and larger mean diameter at the ground).

May 2023 update introduces changes in the default options for graupel/hail fall speeds and collection efficiencies. The original fall speed options (icdx=3; icdxhl=3) from Mansell et al. (2010) are switched to the Milbrandt and Morrison (2013) fall speed curves (icdx=6; icdxhl=6). Because the fall speeds are generally a bit lower, a partially compensating increase in maximum collection efficiency is set by default: ehw0/ehlw0 increased to 0.9. One effect is somewhat reduced total precipitation and cold pool intensity for supercell storms.

  icdx         - fall speed option for graupel (was 3, now is 6)
  icdxhl       - fall speed option for hail (was 3, now is 6)
  ehw0,ehlw0   - Maximim droplet collection efficiencies for graupel (ehw0=0.75, now 0.9)
                 and hail (ehlw0=0.75, now 0.9) 

In summary, to get something closer to previous behavior, use the following:

&nssl2mom_params
  irenuc = 2
  icdx   = 3
  icdxhl = 3
  ehw0   = 0.5
  ehlw0  = 0.75
  ihlcnh = 1
/

For the 2-moment option, shape parameters of graupel and hail can be set in the 
nssl2mom_params namelist:

  alphah       - Shape parameter for graupel (0.0)
  alphahl      - Shape parameter for hail (1.0)



 Microphysics References:

 Mansell, E. R., C. L. Ziegler, and E. C. Bruning, 2010: Simulated electrification of a small
   thunderstorm with two-moment bulk microphysics. J. Atmos. Sci., 67, 171-194, doi:10. 1175/2009JAS2965.1.

  Mansell, E. R. and C. L. Ziegler, 2013: Aerosol effects on simulated storm electrification and
     precipitation in a two-moment bulk microphysics model. J. Atmos. Sci., 70 (7), 2032-2050,
     doi:10.1175/JAS-D-12-0264.1.

 Ziegler, C. L., 1985: Retrieval of thermal and microphysical variables in observed convective storms.
    Part I: Model development and preliminary testing. J. Atmos. Sci., 42, 1487-1509.

 Sedimentation reference:

 Mansell, E. R., 2010: On sedimentation and advection in multimoment bulk microphysics.
    J. Atmos. Sci., 67, 3084-3094, doi:10.1175/2010JAS3341.1.

Possible parameters to adjust:

 ccn : base cloud condensation nuclei concentration (use namelist.input value "nssl_cccn")
 alphah, alphahl : Size distribution shape parameters for graupel (h) and hail (hl)
 infall : changes sedimentation options to see effects (see below)



