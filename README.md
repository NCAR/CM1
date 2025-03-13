<img src="https://www2.mmm.ucar.edu/people/bryan/cm1/NSF-NCAR_Lockup-UCAR-Dark_102523.png" width="30%"/>

## New announcements:
* 24 March 2024:  A new version of CM1 `cm1r21.1` has been released.  [Download the code here](releases).    

* [Here is a sample submission script for NCAR's derecho supercomputer](scripts/cm1run_derecho)  **(Last updated:  24 March 2024)**

## What is CM1?

_In scientific terms_:  CM1 is a three-dimensional, non-hydrostatic, non-linear, time-dependent numerical model designed for idealized studies of atmospheric phenomena.

_In non-scientific terms_:  CM1 is a computer program used for atmospheric research.  It is designed for studies of relatively small-scale processes in the Earth's atmosphere, such as thunderstorms.  

For more information, please read these [answers to frequently asked questions](doc/faq.md) about CM1.

## Code

[Download the code here](releases)  (Most recent version:  `cm1r21.1` available since 24 March 2024)

## Documentation
  
<summary>About CM1</summary>

- [Parallel Performance](doc/pp.md)
- [The governing equations for CM1](doc/cm1_equations.pdf)  **(Last updated:  21 October 2021)**
- [Acknowledgments](doc/ACKNOWLEDGMENTS).
- [Answers to frequently asked questions](doc/faq.md) about CM1.
- [presentation on CM1 parallelization](doc/cm1_parallelization.pdf) (pdf)  (From presentation at NCSA, December 2009)

<summary>Helpful Information for New Users of CM1</summary>

- [A brief summary of how to run cm1](doc/README.md)
- [README.namelist](README.namelist) -- explains the various settings in the namelist.input file.
- [Pre-configured namelist.input files](run/config_files)
- [Soundings for idealized simulations](soundings)
- [Sample submission scripts for NCAR's supercomputers](scripts)

<details>
  <summary>Other Information about CM1</summary>

- [CHANGES in release 21.1](doc/CHANGES) -- new features, modifications, and code fixes for the newest version (24 March 2024)
  - [History of all CHANGES](doc/changes.md) (from `r2` to present)
- [Known problems and fixes](doc/known_problems.md) (last updated:  5 January 2012)
- [Instructions for adding a new microphysics scheme to CM1](doc/new_microphysics.pdf) (pdf)
- [Some useful GrADS scripts](scripts/grads)
- [Some useful programs for MPI users](programs)
</details>

<details>
  <summary>Testing and evaluation of CM1</summary>

Here are reports on some basic tests of the accuracy and capability of CM1.  (Note:  all of these tests have been completed, but I haven't had time to write up the results.  I plan to have all of these posted online in the near future.)

- [Gravity current](https://www2.mmm.ucar.edu/people/bryan/cm1/test_gravity_current)
- [Inertia-gravity waves](https://www2.mmm.ucar.edu/people/bryan/cm1/test_inertia_gravity_waves)
- Two-dimensional mountain waves
- Potential flow over a mountain in dry and moist environments
- Bryan-Fritsch moist benchmark
- Large eddy simulation of the convective boundary layer
- A comparison of axisymmetric and three-dimensional simulations of a tropical cyclone
</details>

<details>
  <summary>Research Results</summary>

**Peer-reviewed articles** that use CM1: (Please contact George Bryan if you have something to add to this list.) 
_Last updated:  March 2024_

Check out this <a href="doc/cm1journals.md">list of journals that have published articles using CM1</a>.

![cm1pubs-211217.png](https://www2.mmm.ucar.edu/people/bryan/cm1/cm1pubs-211217.png)

<ul>
  <details>
    <summary>2024</summary>
  <ul>
    <!--451-->
    <li>Bhattacharya, A, 2024: <a href="https://doi.org/10.1007/s10652-024-09972-2">Evaluation of energy consistent entrainment rate closure for cloudy updrafts</a>. <i>Environ. Fluid Mech.,</i> doi:10.1007/s10652-024-09972-2.
    <!--450-->
    <li>Ross, T. I. D., and S. Lasher-Trapp, 2024: <a href="https://doi.org/10.1175/MWR-D-23-0154.1">On CCN Effects upon Convective Cold Pool Timing and Features</a>. <i>Mon. Wea. Rev.,</i> <b>152,</b> 891–906, doi:10.1175/MWR-D-23-0154.1.
    <!--449-->
    <li>Markowski, P. M., 2024: <a href="https://doi.org/10.1175/JAS-D-23-0161.1">A New Pathway for Tornadogenesis Exposed by Numerical Simulations of Supercells in Turbulent Environments</a>. <i>J. Atmos. Sci.,</i> <b>81,</b> 481–518, doi:10.1175/JAS-D-23-0161.1.
    <!--448-->
    <li>Fu, H., and M. O’Neill, 2024: <a href="https://doi.org/10.1175/JAS-D-23-0170.1">The Small-Amplitude Dynamics of Spontaneous Tropical Cyclogenesis. Part I: Experiments with Amplified Longwave Radiative Feedback</a>. <i>J. Atmos. Sci.,</i> <b>81,</b> 381–399, doi:10.1175/JAS-D-23-0170.1.
    <!--447-->
    <li>Régibeau-Rockett, L., O. M. Pauluis, and M. E. O’Neill, 2024: <a href="https://doi.org/10.1175/JCLI-D-22-0877.1">Investigating the Relationship between Sea Surface Temperature and the Mechanical Efficiency of Tropical Cyclones</a>. <i>J. Climate,</i> <b>37,</b> 439–456, doi:10.1175/JCLI-D-22-0877.1. 
    <!--446-->
    <li>Chen, X., and F. D. Marks, 2024: <a href="https://doi.org/10.1175/JAS-D-23-0086.1">Parameterizations of Boundary Layer Mass Fluxes in High-Wind Conditions for Tropical Cyclone Simulations</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-23-0086.1. 
    <!--445-->
    <li>Powell, S. W., 2024: <a href="https://doi.org/10.1175/JAS-D-23-0065.1">Updraft Width Implications for Cumulonimbus Growth in a Moist Marine Environment</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-23-0065.1. 
    <!--444-->
    <li>Muehr, A. J., J. H. Ruppert, M. D. Flournoy, and J. M. Peters, 2024: <a href="https://doi.org/10.1175/JAS-D-23-0082.1">The Influence of Midlevel Shear and Horizontal Rotors on Supercell Updraft Dynamics</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-23-0082.1. 
    <!--443-->
    <li>Feldmann, M., R. Rotunno, U. Germann, and A. Berne, 2024: <a href="https://doi.org/10.1175/MWR-D-22-0350.1">Supercell thunderstorms in complex topography - how mountain valleys with lakes can increase occurrence frequency</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-22-0350.1. 
    <!--442-->
    <li>Dong, Y., S. Hua, B. Chen, H. Wang, and T. Hou, 2024: <a href="https://doi.org/10.1016/j.atmosres.2023.107218">Numerical simulation of a pulse hailstorm in the plateau region in southwestern China</a>. <i>Atmospheric Research</i> <b>299</b>, doi:10.1016/j.atmosres.2023.107218.
  </ul>
  </details>

  <details>
  <summary>2023</summary>
  <ul>
    <!--441-->
    <li>Jeong, J.-H., Witte, M. K., and Smalley, M., 2023: <a href="https://doi.org/10.1029/2023JD039081">Effects of wind shear and aerosol conditions on the organization of precipitating marine stratocumulus clouds</a>. <i>Journal of Geophysical Research: Atmospheres,</i> <b>128,</b> e2023JD039081.
    <!--440-->
    <li>Yu, C., B. Tang, and R. G. Fovell, 2023: <a href="https://doi.org/10.1175/JAS-D-23-0048.1">Diverging Behaviors of Simulated Tropical Cyclones in Moderate Vertical Wind Shear</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-23-0048.1. 
    <!--439-->
    <li>Goldacker, N. A., and M. D. Parker, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0253.1">Assessing the Comparative Effects of Storm-Relative Helicity Components within Right-Moving Supercell Environments</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 2805–2822, doi:10.1175/JAS-D-22-0253.1.
    <!--438-->
    <li>Morrison, H., N. Jeevanjee, D. Lecoanet, and J. M. Peters, 2023: <a href="https://doi.org/10.1175/JAS-D-23-0063.1">What controls the entrainment rate of dry buoyant thermals with varying initial aspect ratio?</a> <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-23-0063.1. 
    <!--437-->
    <li>Hutson, A., and C. Weiss, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0288.1">Using Ensemble Sensitivity Analysis to Identify Storm Characteristics Associated with Tornadogenesis in High Resolution Simulated Supercells</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-22-0288.1.
    <!--436-->
    <li>Jiang, Q., and D. T. Dawson, 2023: <a href="https://doi.org/10.1175/MWR-D-23-0050.1">The impact of surface drag on the structure and evolution of surface boundaries associated with tornadogenesis in simulated supercells</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-23-0050.1. 
    <!--435-->
    <li>Gray, K. T., and J. W. Frame, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0309.1">Investigating the Development and Characteristics of Streamwise Vorticity Currents Produced by Outflow Surges in Simulated Supercell Thunderstorms</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-22-0309.1.
    <!--434-->
    <li>Natoli, M. B., and E. D. Maloney, 2023: <a href="https://doi.org/10.1175/JCLI-D-22-0824.1">Environmental Controls on the Tropical Island Diurnal Cycle in the Context of Intraseasonal Variability</a>. <i>J. Climate,</i> <b>36,</b> 7465–7485, doi:10.1175/JCLI-D-22-0824.1. 
    <!--433-->
    <li>Rozoff, C. M., D. S. Nolan, G. H. Bryan, E. A. Hendricks, and J. Knievel, 2023: <a href="https://doi.org/10.1175/JAMC-D-23-0024.1">Large-eddy Simulations of the Tropical Cyclone Boundary Layer at Landfall in an Idealized Urban Environment</a>. <i>J. Appl. Meteor. Climatol.,</i> doi:10.1175/JAMC-D-23-0024.1. 
    <!--432-->
    <li>Yang, H., Y. Du, and J. Wei, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0212.1">Generation of Multiple Gravity Wave Couplets from Convection</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 2323–2343, doi:10.1175/JAS-D-22-0212.1. 
    <!--431-->
    <li>Chandrakar, K. K., H. Morrison, and R. A. Shaw, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0256.1">Lagrangian and Eulerian Supersaturation Statistics in Turbulent Cloudy Rayleigh–Bénard Convection: Applications for LES Subgrid Modeling</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 2261–2285, doi:10.1175/JAS-D-22-0256.1. 
    <!--430-->
    <li>Chen, J., and X. Shi, 2023: <a href="https://doi.org/10.1175/JCLI-D-23-0018.1">Quantifying Global Warming Response of the Orographic Precipitation in a Typhoon Environment with Large-Eddy Simulations</a>. <i>J. Climate,</i> <b>36,</b> 6951–6966, doi:10.1175/JCLI-D-23-0018.1. 
    <!--429-->
    <li>Peters, J. M., D. R. Chavas, C. Su, H. Morrison, and B. E. Coffer, 2023: <a href="https://doi.org/10.1175/JAS-D-23-0003.1">An Analytic Formula for Entraining CAPE in Midlatitude Storm Environments</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 2165–2186, doi:10.1175/JAS-D-23-0003.1. 
    <!--428-->
    <li>Adams-Selin, R. D., 2023: <a href="https://doi.org/10.1175/MWR-D-22-0345.1">A Three-Dimensional Hail Trajectory Clustering Technique</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 2361–2375, doi:10.1175/MWR-D-22-0345.1. 
    <!--427-->
    <li>Ma, Z., Yan, X., and Fei, J, 2023: <a href="https://doi.org/10.1029/2023GL104578">Quantifying the rightward bias extent of tropical cyclones' cold wakes</a>. <i>Geophysical Research Letters,</i> <b>50,</b> e2023GL104578, doi:10.1029/2023GL104578. 
    <!--426-->
    <li>Coffer, B. E., M. D. Parker, J. M. Peters, and A. R. Wade, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0269.1">Supercell Low-Level Mesocyclones: Origins of Inflow and Vorticity</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 2205–2232, doi:10.1175/MWR-D-22-0269.1.
    <!--425-->
    <li>Flournoy, M. D., and E. N. Rasmussen, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0069.1">The Influence of Convection Initiation Strength on Subsequent Simulated Supercell Evolution</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 2179–2203, doi:10.1175/MWR-D-22-0069.1.
    <!--424-->
    <li>Oguejiofor, C. N., C. E. Wainwright, J. E. Rudzin, and D. H. Richter, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0158.1">Onset of Tropical Cyclone Rapid Intensification: Evaluating the Response to Length Scales of Sea Surface Temperature Anomalies</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 1971–1994, doi:10.1175/JAS-D-22-0158.1.
    <!--423-->
    <li>Dai, Y., M. S. Torn, I. N. Williams, and W. D. Collins, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0214.1">Longwave Radiative Effects beyond the Initial Intensification Phase of Tropical Cyclones</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 1829–1845, doi:10.1175/JAS-D-22-0214.1.
    <!--422-->
    <li>Loeffler, S. D., M. R. Kumjian, P. M. Markowski, B. E. Coffer, and M. D. Parker, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0228.1">Investigating the Relationship between Polarimetric Radar Signatures of Hydrometeor Size Sorting and Tornadic Potential in Simulated Supercells</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 1863–1884, doi:10.1175/MWR-D-22-0228.1.
    <!--421-->
    <li>Schecter, D. A., 2023: <a href="https://doi.org/10.1175/JAS-D-22-0188.1">Intensification Rates of Tropical Cyclone–Like Vortices in a Model with Downtilt Diabatic Forcing and Oceanic Surface Drag</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 1787–1814, doi:10.1175/JAS-D-22-0188.1.
    <!--420-->
    <li>LeBel, L. J., and P. M. Markowski, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0176.1">An Analysis of the Impact of Vertical Wind Shear on Convection Initiation Using Large-Eddy Simulations: Importance of Wake Entrainment</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 1667–1688, doi:10.1175/MWR-D-22-0176.1.
    <!--419-->
    <li>Li, Y., Y. Wang, and Z. Tan, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0186.1">Is the Outflow-Layer Inertial Stability Crucial to the Energy Cycle and Development of Tropical Cyclones?</a> <i>J. Atmos. Sci.,</i> <b>80,</b> 1605–1620, doi:10.1175/JAS-D-22-0186.1.
    <!--418-->
    <li>Wang, A., Y. Pan, G. H. Bryan, and P. M. Markowski, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0060.1">Modeling Near-Surface Turbulence in Large-Eddy Simulations of a Tornado: An Application of Thin Boundary Layer Equations</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 1587–1607, doi:10.1175/MWR-D-22-0060.1.
    <!--417-->
    <li>Parker, M. D., 2023: <a href="https://doi.org/10.1175/JAS-D-22-0195.1">How Well Must Surface Vorticity Be Organized for Tornadogenesis?</a> <i>J. Atmos. Sci.,</i> <b>80,</b> 1433–1448, doi:10.1175/JAS-D-22-0195.1.
    <!--416-->
    <li>Dahl, J. M. L., and J. Fischer, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0145.1">On the Origins of Vorticity in a Simulated Tornado-Like Vortex</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 1361–1380, doi:10.1175/JAS-D-22-0145.1.
    <!--415-->
    <li>Jo, E., and S. Lasher-Trapp, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0168.1">Entrainment in a Simulated Supercell Thunderstorm. Part III: The Influence of Decreased Environmental Humidity and General Effects upon Precipitation Efficiency</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 1107–1122, doi:10.1175/JAS-D-22-0168.1.
    <!--414-->
    <li>Schumacher, R. S., S. J. Childs, and R. D. Adams-Selin, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0103.1">Intense Surface Winds from Gravity Wave Breaking in Simulations of a Destructive Macroburst</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 775–793, doi:10.1175/MWR-D-22-0103.1.
    <!--413-->
    <li>Yu, C., B. Tang, and R. G. Fovell, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0200.1">Tropical Cyclone Tilt and Precession in Moderate Shear: Precession Hiatus in a Critical Shear Regime</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 909–932, doi:10.1175/JAS-D-22-0200.1. 
    <!--412-->
    <li>Kieu, C., W. Cai, and W.-T. Fan, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0115.1">On the Existence of Low-Dimensional Chaos of the Tropical Cyclone Intensity in an Idealized Axisymmetric Simulation</a>. <i>J. Atmos. Sci.,</i> <b>80,</b> 797–811, doi:10.1175/JAS-D-22-0115.1. 
    <!--411-->
    <li>Fei, R., and Y. Wang, 2023: <a href="https://doi.org/10.1175/JAS-D-22-0014.1">How Does Horizontal Diffusion Influence the Intensification and Maximum Intensity of Numerically Simulated Tropical Cyclones?</a> <i>J. Atmos. Sci.,</i> <b>80,</b> 705–723, doi:10.1175/JAS-D-22-0014.1. 
    <!--410-->
    <li>Zhang, S., and Coauthors, 2023: <a href="https://doi.org/10.1175/MWR-D-22-0091.1">Dynamics Governing a Simulated Bow-and-Arrow-Type Mesoscale Convective System</a>. <i>Mon. Wea. Rev.,</i> <b>151,</b> 603–623, doi:10.1175/MWR-D-22-0091.1. 
    <!--409-->
    <li>Pearson, C., T. Yu, D. Bodine, S. Torres, and A. Reinhart, 2023: <a href="https://doi.org/10.1175/JTECH-D-22-0130.1">A Framework for Comparisons of Downburst Precursor Observations Using an All-Digital Phased-Array Weather Radar</a>. <i>J. Atmos. Oceanic Technol.,</i> <b>40,</b> 919–938, doi:10.1175/JTECH-D-22-0130.1.
    <!--408-->
    <li>Naylor, J., and J. P. Mulholland, 2023: <a href="https://doi.org/10.1029/2022JD037237">The impact of vertical wind shear on the outcome of interactions between squall lines and cities</a>. <i>Journal of Geophysical Research: Atmospheres,</i> <b>128,</b> e2022JD037237, doi:10.1029/2022JD037237. 
    <!--407-->
    <li>Woods, M. J., R. J. Trapp, and H. M. Mallinson, 2023:  <a href="https://doi.org/10.1029/2023GL104796">The impact of human-induced climate change on future tornado intensity as revealed through multi-scale modeling</a>. <i>Geophysical Research Letters,</i> <b>50,</b> e2023GL104796, doi:10.1029/2023GL104796. 
    <!--406-->
    <li>Zhang, D., and Z. Ma, 2023:  <a href="https://doi.org/10.1029/2023JD038580">The generalized application of a new surface pressure tendency equation in synoptic weather systems</a>. <i>Journal of Geophysical Research: Atmospheres,</i> <b>128,</b> e2023JD038580, doi:10.1029/2023JD038580. 
    <!--405-->
    <li>Finley, C. A., M. Elmore, L. Orf, and B. D. Lee, 2023: <a href="https://doi.org/10.1029/2022GL100005">Impact of the streamwise vorticity current on low-level mesocyclone development in a simulated supercell</a>. <i>Geophysical Research Letters,</i> <b>50,</b> e2022GL100005, doi:10.1029/2022GL100005. 
    <!--404-->
    <li>Dogra, G., A. Dewan, and S. Sahany, 2023: <a href="https://doi.org/10.3390/fluids8020051">Understanding Atmospheric Convection Using Large Eddy Simulation</a>. <i>Fluids,</i> <b>8</b>, 51, doi:10.3390/fluids8020051. 
    <!--403-->
    <li>Wu, F., and K. Lombardo, 2023: <a href="https://doi.org/10.1029/2023GL102825">The impact of offshore-propagating squall lines on coastal-mountain flows</a>. <i>Geophysical Research Letters,</i> <b>50,</b> e2023GL102825, doi:10.1029/2023GL102825. 
    <!--402-->
    <li>Labriola, J. D., J. A. Gibbs, and L. J. Wicker, 2023: <a href="https://doi.org/10.5194/gmd-16-1779-2023">A method for generating a quasi-linear convective system suitable for observing system simulation experiments</a>.  <i>Geoscientific Model Development</i>, <b>16,</b> 1779-1799, doi:10.5194/gmd-16-1779-2023. 
    <!--401-->
    <li>Vich, M., and R. Romero, 2023: <a href="https://doi.org/10.1016/j.atmosres.2023.106784">Exploring severe weather environments using CM1 simulations: The 29 August 2020 event in the Balearic Islands</a>.  <i>Atmospheric Research,</i> <b>290,</b> doi:10.1016/j.atmosres.2023.106784. 
    <!--400-->
    <li>Chen, J., and Chavas, D. R., 2023: <a href="https://doi.org/10.1175/JAS-D-22-0156.1">A Model for the Tropical Cyclone Wind Field Response to Idealized Landfall</a>. <i>Journal of the Atmospheric Sciences</i>, in press, doi:10.1175/JAS-D-22-0156.1. 
    <!--399-->
    <li>Fischer, J., and Dahl, J. M. L., 2023: <a href="https://doi.org/10.1175/MWR-D-22-0026.1">Supercell-External Storms and Boundaries Acting as Catalysts for Tornadogenesis</a>. <i>Monthly Weather Review,</i> <b>151,</b> 23-38, doi:10.1175/MWR-D-22-0026.1.
    <!--398-->
    <li>Peters, J. M., Coffer, B. E., Parker, M. D., Nowotarski, C. J., Mulholland, J. P., Nixon, C. J., and Allen, J. T., 2023: <a href="https://doi.org/10.1175/JAS-D-22-0114.1">Disentangling the Influences of Storm-Relative Flow and Horizontal Streamwise Vorticity on Low-Level Mesocyclones in Supercells</a>. <i>Journal of the Atmospheric Sciences,</i> <b>80,</b> 129-149, doi:10.1175/JAS-D-22-0114.1. 
    <!--397-->
    <li>Natoli, M. B., and Maloney, E. D., 2023: <a href="https://doi.org/10.1175/JAS-D-22-0045.1">The Tropical Diurnal Cycle under Varying States of the Monsoonal Background Wind</a>. <i>Journal of the Atmospheric Sciences,</i> <b>80,</b> 235-258, doi:10.1175/JAS-D-22-0045.1. 
    <!--396-->
    <li>Wang, Y., Tan, Z., and Li, Y., 2023: <a href="https://doi.org/10.1175/JAS-D-22-0135.1">Some Refinements to the Most Recent Simple Time-Dependent Theory of Tropical Cyclone Intensification and Sensitivity</a>. <i>Journal of the Atmospheric Sciences,</i> <b>80,</b> 321-335, doi:10.1175/JAS-D-22-0135.1. 
    <!--395-->
    <li>Finley, C. A., Elmore, M., Orf, L., and Lee, B. D., 2023: <a href="https://doi.org/10.1029/2022GL100005">Impact of the streamwise vorticity current on low-level mesocyclone development in a simulated supercell</a>. <i>Geophysical Research Letters,</i> <b>50,</b> e2022GL100005, doi:10.1029/2022GL100005. 
  </ul>
  </details>

  <details>
    <summary>2022</summary>
  <ul>
    <!--394-->
    <li>Boyer, C. H., and Keeler, J. M., 2022: <a href="https://doi.org/10.1175/JAMC-D-22-0017.1">Evaluation and Improvement of an Inflow-Nudging Technique for Idealized Simulations of Convective Boundary Layers</a>. <i>Journal of Applied Meteorology and Climatology,</i> <b>61,</b> 1843-1860, doi:10.1175/JAMC-D-22-0017.1. 
    <!--393-->
    <li>Coffer, B. E., and M. D. Parker, 2022:  <a href="https://doi.org/10.1121/10.0009400">Infrasound signals in simulated nontornadic and pre-tornadic supercells"</a>. <i>Journal of the Acoustical Society of America</i>, doi:10.1121/10.0009400. 
    <!--392-->
    <li>Nelson, T. C., Marquis, J., Peters, J. M., and Friedrich, K., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0226.1">Environmental Controls on Simulated Deep Moist Convection Initiation Occurring during RELAMPAGO-CACTI</a>. <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 1941-1964, doi:10.1175/JAS-D-21-0226.1. 
    <!--391-->
    <li>Chandrakar, K. K., Morrison, H., Grabowski, W. W., and Bryan, G. H., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0138.1">Comparison of Lagrangian Superdroplet and Eulerian Double-Moment Spectral Microphysics Schemes in Large-Eddy Simulations of an Isolated Cumulus Congestus Cloud</a>.  <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 1887-1910, doi:10.1175/JAS-D-21-0138.1. 
    <!--390-->
    <li>Murdzek, S. S., Richardson, Y. P., Markowski, P. M., and Kumjian, M. R., 2022: <a href="https://doi.org/10.1175/MWR-D-21-0258.1">How the Environmental Lifting Condensation Level Affects the Sensitivity of Simulated Convective Storm Cold Pools to the Microphysics Parameterization</a>. <i>Monthly Weather Review,</i> <b>150,</b> 2527-2552, doi:10.1175/MWR-D-21-0258.1. 
    <!--389-->
    <li>Weinkaemmerer, J., Ďurán, I. B., and Schmidli, J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0195.1">The Impact of Large-Scale Winds on Boundary Layer Structure, Thermally Driven Flows, and Exchange Processes over Mountainous Terrain</a>.  <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 2685-2701, doi:10.1175/JAS-D-21-0195.1. 
    <!--388-->
    <li>Hernández Pardo, L., Morrison, H., Lauritzen, P. H., and Pöhlker, M., 2022: <a href="https://doi.org/10.1175/MWR-D-22-0025.1">Impact of Advection Schemes on Tracer Interrelationships in Large-Eddy Simulations of Deep Convection</a>. <i>Monthly Weather Review,</i> <b>150,</b> 2765-2785, doi:10.1175/MWR-D-22-0025.1. 
    <!--387-->
    <li>Wang, D., Lin, Y., and Chavas, D. R., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0325.1">Tropical Cyclone Potential Size</a>. <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 3001-3025, doi:10.1175/JAS-D-21-0325.1.
    <!--386-->
    <li>Morrison, H., Jeevanjee, N., and Yano, J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0274.1">Dynamic Pressure Drag on Rising Buoyant Thermals in a Neutrally Stable Environment</a>. <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 3045-3063, 10.1175/JAS-D-21-0274.1. 
    <!--385-->
    <li>Rotunno, R., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0306.1">Supergradient Winds in Simulated Tropical Cyclones</a>. <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 2075-2086, doi:10.1175/JAS-D-21-0306.1. 
    <!--384-->
    <li>Fu, S., R. Rotunno, and H. Xue, 2022: <a href="https://doi.org/10.5194/acp-22-7727-2022">Convective updrafts near sea-breeze fronts</a>.  <i>Atmos. Chem. Phys.,</i> doi:10.5194/acp-22-7727-2022. 
    <!--383-->
    <li>Patra, M., W.-T. Fan, and C. Kieu, 2022: <a href="https://doi.org/10.3389/feart.2022.893781">Sensitivity of Tropical Cyclone Intensity Variability to Different Stochastic Parameterization Methods</a>.  <i>Frontiers in Earth Science,</i> doi:10.3389/feart.2022.893781. 
    <!--382-->
    <li>Gordon, A. E., and Homeyer, C. R., 2022:  <a href="https://doi.org/10.1029/2022JD036713">Sensitivities of cross-tropopause transport in midlatitude overshooting convection to the lower stratosphere environment</a>. <i>Journal of Geophysical Research: Atmospheres,</i> <b>127,</b> e2022JD036713, doi:10.1029/2022JD036713. 
    <!--381-->
    <li>Bickle, M., Marsham, J. H., Griffiths, S. D., Ross, A. N., and Crook, J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0025.1">The Influence of the Diurnal Cycle in Wind Shear and Thermodynamics on Squall Lines in the West African Monsoon</a>. <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 2125-2143, doi:10.1175/JAS-D-21-0025.1. 
    <!--380-->
    <li>Gowan, T. M., Steenburgh, W. J., and Minder, J. R., 2022: <a href="https://doi.org/10.1175/MWR-D-21-0314.1">Orographic Effects on Landfalling Lake-Effect Systems</a>. <i>Monthly Weather Review,</i> <b>150,</b> 2013-2031, doi:10.1175/MWR-D-21-0314.1. 
    <!--379-->
    <li>Li, Y., Wang, Y., and Tan, Z.-M., 2022: <a href="https://doi.org/10.1029/2022JD037039">Why does the initial wind profile inside the radius of maximum wind matter to tropical cyclone development?</a> <i>Journal of Geophysical Research: Atmospheres,</i> <b>127,</b> e2022JD037039, doi:10.1029/2022JD037039. 
    <!--378-->
    <li>Martinez, J., Davis, C. A., and Bell, M. M., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0302.1">Eyewall Asymmetries and Their Contributions to the Intensification of an Idealized Tropical Cyclone Translating in Uniform Flow</a>. <i>Journal of the Atmospheric Sciences,</i> <b>79,</b> 2471-2491, doi:10.1175/JAS-D-21-0302.1. 
    <!--377-->
    <li>Labriola, J.D., and Wicker, L.J., 2022: <a href="https://doi.org/10.1002/qj.4348">Creating physically coherent and spatially correlated perturbations to initialize high-resolution ensembles of simulated convection</a>. <i>Quarterly Journal of the Royal Meteorological Society,</i> doi:10.1002/qj.4348. 
    <!--376-->
    <li>Bannigan, N., Orf, L., and Savory, E, 2022: <a href="https://doi.org/10.1007/s10546-022-00739-0">Tracking the Centre of Asymmetric Vortices Using Wind Velocity Vector Data Fields</a>. <i>Boundary-Layer Meteorol,</i> doi:10.1007/s10546-022-00739-0. 
    <!--375-->
    <li>Chandrakar, K. K., Morrison, H., and Witte, M., 2022: <a href="https://doi.org/10.1029/2022GL100511">Evolution of droplet size distributions during the transition of an ultraclean stratocumulus cloud system to open cell structure: An LES investigation using Lagrangian microphysics</a>. <i>Geophysical Research Letters,</i> <b>49,</b> e2022GL100511, doi:10.1029/2022GL100511. 
    <!--374-->
    <li>Chen, X., 2022: <a href="https://doi.org/10.1029/2022MS003088">How do planetary boundary layer schemes perform in hurricane conditions: A comparison with large-eddy simulations</a>. <i>Journal of Advances in Modeling Earth Systems,</i> <b>14,</b> doi:10.1029/2022MS003088. 
    <!--373-->
    <li>Weinkaemmerer, J., Ďurán, I.B., Westerhuis, S. and Schmidli, J., 2022: <a href="https://doi.org/10.1002/qj.4372">Stratus over rolling terrain: Large-eddy simulation reference and sensitivity to grid spacing and numerics</a>. <i>Quarterly Journal of the Royal Meteorological Society,</i> doi:10.1002/qj.4372. 
    <!--372-->
    <li>Jo, E., and S. Lasher-Trapp, S., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0289.1">Entrainment in a Simulated Supercell Thunderstorm. Part II: The Influence of Vertical Wind Shear and General Effects upon Precipitation</a>.  <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0289.1. 
    <!--371-->
    <li>Lovell, L. T., and Parker, M. D., 2022: <a href="https://doi.org/10.1175/WAF-D-21-0133.1">Simulated QLCS Vortices in a High-Shear, Low-CAPE Environment</a>. <i>Weather and Forecasting,</i> doi:10.1175/WAF-D-21-0133.1. 
    <!--370-->
    <li>Peters, J. M., Morrison, H., Nelson, T. C., Marquis, J. N., Mulholland, J. P., and Nowotarski, C. J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0145.1">The Influence of Shear on Deep Convection Initiation. Part I: Theory</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0145.1. 
    <!--369-->
    <li>Peters, J. M., Morrison, H., Nelson, T. C., Marquis, J. N., Mulholland, J. P., and Nowotarski, C. J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0144.1">The Influence of Shear on Deep Convection Initiation. Part II: Simulations</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0144.1. 
    <!--368-->
    <li>Chen, X., Bryan, G. H., Hazelton, A., Marks, F. D., and Fitzpatrick, P., 2022: <a href="https://doi.org/10.1175/WAF-D-21-0168.1">Evaluation and Improvement of a TKE-Based Eddy-Diffusivity Mass-Flux (EDMF) Planetary Boundary Layer Scheme in Hurricane Conditions</a>.  <i>Weather and Forecasting,</i> doi:10.1175/WAF-D-21-0168.1 
    <!--367-->
    <li>Powell, S. W., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0155.1">Criticality in the Shallow-to-Deep Transition of Simulated Tropical Marine Convection</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0155.1.
    <!--366-->
    <li>Nardi, K. M., Zarzycki, C. M., Larson, V. E., and Bryan, G. H., 2022: <a href="https://doi.org/10.1175/MWR-D-21-0186.1">Assessing the Sensitivity of the Tropical Cyclone Boundary Layer to the Parameterization of Momentum Flux in the Community Earth System Model</a>. <i>Monthly Weather Review,</i> doi:10.1175/MWR-D-21-0186.1. 
    <!--365-->
    <li>Singh, I., Nesbitt, S. W., and Davis, C. A., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0007.1">Quasi-idealized numerical simulations of processes involved in orogenic convection initiation over the Sierras de Córdoba mountains</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0007.1. 
    <!--364-->
    <li>Singh, M. S., and M. E. O'Neill, 2022: <a href="https://doi.org/10.1103/RevModPhys.94.015001">The climate system and the second law of thermodynamics</a>. <i>Reviews of Modern Physics,</i> doi:10.1103/RevModPhys.94.015001. 
    <!--363-->
    <li>Ye, H., Ma, Z. and Fei, J., 2022: <a href="https://doi.org/10.1007/s13351-022-1120-8">Uncertainty in TC Maximum Intensity with Fixed Ratio of Surface Exchange Coefficients for Enthalpy and Momentum</a>. <i>J Meteorol Res.</i> doi:10.1007/s13351-022-1120-8. 
    <!--362-->
    <li>Schecter, D. A., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0051.1">Intensification of Tilted Tropical Cyclones over Relatively Cool and Warm Oceans in Idealized Numerical Simulations</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0051.1. 
    <!--361-->
    <li>Lefèvre, M., X. Tan, E. K. H. Lee, and R. T. Pierrehumbert, 2022: <a href="https://doi.org/10.3847/1538-4357/ac5e2d">Cloud-convection Feedback in Brown Dwarf Atmospheres</a>.  <i>The Astrophysical Journal,</i> <b>929,</b> doi:10.3847/1538-4357/ac5e2d. 
    <!--360-->
    <li>Kirshbaum, D. J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0197.1">Large-Eddy Simulations of Convection Initiation over Heterogeneous, Low Terrain</a>. <i>Journal of the Atmospheric Sciences</i>, doi:10.1175/JAS-D-21-0197.1. 
    <!--359-->
    <li>Chandrakar, K. K., Morrison, H., Grabowski, W. W., Bryan, G. H., and Shaw, R. A., 2022. <a href="https://doi.org/10.1175/JAS-D-21-0250.1">Supersaturation Variability from Scalar Mixing: Evaluation of a New Subgrid-Scale Model Using Direct Numerical Simulations of Turbulent Rayleigh–Bénard Convection</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0250.1. 
    <!--358-->
    <li>Yan, B., Y. Yuan, C. Ma, Z. Dong, H. Huang, and Z. Wang, 2022: <a href="https://doi.org/10.1016/j.jobe.2021.103738">Modeling of downburst outflows and wind pressures on a high-rise building under different terrain conditions</a>.  <i>Journal of Building Engineering,</i> doi:10.1016/j.jobe.2021.103738. 
    <!--357-->
    <li>Fovell, R. G., M. J. Brewer, and R. J. Garmong, 2022: <a href="https://doi.org/10.3390/atmos13050765">The December 2021 Marshall Fire: Predictability and Gust Forecasts from Operational Models</a>.  <i>Atmosphere,</i> doi:10.3390/atmos13050765. 
    <!--356-->
    <li>Williams, G.J., 2022: <a href="https://doi.org/10.1007/s00703-022-00900-x">Idealized simulations of the diurnal variation within the tropical cyclone boundary layer</a>. <i>Meteorol Atmos Phys,</i> <b>134,</b> doi:10.1007/s00703-022-00900-x. 
    <!--355-->
    <li>Morrison, H., P. Lawson, P., and K. K. Chandrakar, 2022: <a href="https://doi.org/10.1029/2021JD035711">Observed and bin model simulated evolution of drop size distributions in high-based cumulus congestus over the United Arab Emirates</a>. <i>Journal of Geophysical Research: Atmospheres</i>, <b>127,</b> e2021JD035711, doi:10.1029/2021JD035711. 
    <!--354-->
    <li>Stauffer, C. L., and A. A. Wing, 2022: <a href="https://doi.org/10.1029/2021MS002917">Properties, changes, and controls of deep-convecting clouds in radiative-convective equilibrium</a>. <i>Journal of Advances in Modeling Earth Systems</i>, doi:10.1029/2021MS002917. 
    <!--353-->
    <li>Persing, J., and Montgomery, M. T., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0143.1">Does the Rotating Convection Paradigm Describe Secondary Eyewall Formation in Idealized Three-dimensional Simulations?</a> <i>Journal of the Atmospheric Sciences</i>, doi:10.1175/JAS-D-21-0143.1. 
    <!--352-->
    <li>Sokol, A. B., and D. L. Hartmann, 2022: <a href="https://doi.org/10.1029/2022MS003045">Congestus mode invigoration by convective aggregation in simulations of radiative-convective equilibrium</a>. <i>Journal of Advances in Modeling Earth Systems,</i> doi:10.1029/2022MS003045. 
    <!--351-->
    <li>Fischer, J., and Dahl, J. M. L., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0181.1">Transition of Near-Ground Vorticity Dynamics During Tornadogenesis</a>. <i>Journal of the Atmospheric Sciences</i>, doi:10.1175/JAS-D-21-0181.1.
    <!--350-->
    <li>Ma, Z., and Fei, J., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0014.1">A Comparison between Moist and Dry Tropical Cyclones: The Low Effectiveness of Surface Sensible Heat Flux in Storm Intensification</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0014.1.
    <!--349-->
    <li>Peters, J. M., Mulholland, J. P., and Chavas, D. R., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0118.1">Generalized lapse rate formulas for use in entraining CAPE calculations</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0118.1. 
    <!--348-->
    <li>Morrison, H., Peters, J. M., Chandakar, K. K., and Sherwood, S. C., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0056.1">Influences of environmental relative humidity and horizontal scale of sub-cloud ascent on deep convective initiation</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0056.1. 
    <!--347-->
    <li>Lin, Y., and Kumjian, M. R., 2022: <a href="https://doi.org/10.1175/JAS-D-21-0054.1">Influences of CAPE on Hail Production in Simulated Supercell Storms</a>. <i>Journal of the Atmospheric Sciences</i>, doi:10.1175/JAS-D-21-0054.1.
  </ul>
  </details>

  <details><summary>2021</summary>
  <ul>
    <!--346-->
    <li>Wang, Y., Li, Y., and Xu, J., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0169.1">A New Time-Dependent Theory of Tropical Cyclone Intensification</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0169.1. 
    <!--345-->
    <li>Tan, X., M. Lefèvre, and R. T. Pierrehumbert, 2021: <a href="https://doi.org/10.3847/2041-8213/ac3e69">Convection modeling of pure-steam atmospheres</a>.  <i>Astrophysical Journal Letters,</i> <b>923,</b> L15, doi:10.3847/2041-8213/ac3e69.
    <!--344-->
    <li>Li, Y., Wang, Y., Lin, Y., and Wang, X., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0129.1">Why Does Rapid Contraction of the Radius of Maximum Wind Precede Rapid Intensification in Tropical Cyclones?</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-21-0129.1. 
    <!--343-->
    <li>Lasher-Trapp, S., Scott, E. L., Järvinen, E., Schnaiter, M., Waitz, F., DeMott, P. J., et al., 2021: <a href="https://doi.org/10.1029/2021JD035479">Observations and Modeling of Rime-Splintering in Southern Ocean Cumuli</a>. <i>Journal of Geophysical Research: Atmospheres,</i> <b>126,</b> e2021JD035479, doi:10.1029/2021JD035479. 
    <!--342-->
    <li>Seeley, J.T., Wordsworth, R.D., 2021: <a href="https://doi.org/10.1038/s41586-021-03919-z">Episodic deluges in simulated hothouse climates</a>. <i>Nature,</i> <b>599,</b> 74–79, doi:10.1038/s41586-021-03919-z. 
    <!--341-->
    <li>Stern, D. P., Bryan, G. H., Lee, C., and Doyle, J. D., 2021: <a href="https://doi.org/10.1175/MWR-D-21-0059.1">Estimating the Risk of Extreme Wind Gusts in Tropical Cyclones Using Idealized Large-Eddy Simulations and a Statistical-Dynamical Model</a>. <i>Monthly Weather Review,</i> doi:10.1175/MWR-D-21-0059.1. 
    <!--340-->
    <li>Drueke, S., Kirshbaum, D. J., and Kollias, P., 2021: <a href="https://doi.org/10.5194/acp-21-14039-2021">Environmental sensitivities of shallow-cumulus dilution – Part 2: Vertical wind profile</a>. <i>Atmos. Chem. Phys.,</i> <b>21,</b> 14039–14058, doi:10.5194/acp-21-14039-2021.
    <!--339-->
    <li>O'Neill, M. E., L. Orf, G. M.  Heymsfield, and K. Halbert, 2021:  <a href="https://doi.org/10.1126/science.abh3857">Hydraulic jump dynamics above supercell thunderstorms</a>. <i>Science,</i> <b>373,</b> 1248-1251, doi:10.1126/science.abh3857. 
    <!--338 ... see below-->
    <!--337-->
    <li>Groff, F. P., Adams-Selin, R. D., and Schumacher, R. S., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0208.1">Response of MCS Low-Frequency Gravity Waves to Vertical Wind Shear and Nocturnal Thermodynamic Environments</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0208.1. 
    <!--336-->
    <li>Fu, H., and O’Neill, M., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0087.1">The role of random vorticity stretching in tropical depression genesis</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0087.1. 
    <!--335-->
    <li>Li, L., and P. Chakraborty, 2021: <a href="https://doi.org/10.1103/PhysRevFluids.6.L051801">Birth of a cold core in tropical cyclones past landfall</a>. <i>Physical Review Fluids,</i> <b>6,</b> L051801, doi:10.1103/PhysRevFluids.6.L051801. 
    <!--334-->
    <li>Chen, X., and Bryan, G. H., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0088.1">Role of Advection of Parameterized Turbulence Kinetic Energy in Idealized Tropical Cyclone Simulations</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0088.1. 
    <!--333-->
    <li>Shi, X., and Fan, Y., 2021: <a href="https://doi.org/10.1029/2021MS002632">Modulation of the Bifurcation in Radiative-Convective Equilibrium by Gray-Zone Cloud and Turbulence Parameterizations</a>. <i>Journal of Advances in Modeling Earth Systems,</i> <b>13,</b> e2021MS002632, doi:10.1029/2021MS002632. 
    <!--332-->
    <li>Fan, Y., Chung, Y. T., and Shi, X., 2021: <a href="https://doi.org/10.1029/2021GL095102">The Essential Role of Cloud-Radiation Interaction in Nonrotating Convective Self-Aggregation</a>. <i>Geophysical Research Letters,</i> <b>48,</b> e2021GL095102, doi:10.1029/2021GL095102. 
    <!--331-->
    <li>Gray, K., and Frame, J., 2021: <a href="https://doi.org/10.1175/MWR-D-21-0085.1">The Impact of Midlevel Shear Orientation on the Longevity of and Downdraft Location and Tornado-like Vortex Formation within Simulated Supercells</a>.  <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-21-0085.1. 
    <!--330-->
    <li>Chen, X., Bryan, G. H., Zhang, J. A., Cione, J. J., and  Marks, F. D., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0227.1">A Framework for Simulating the Tropical-Cyclone Boundary Layer Using Large-Eddy Simulation and Its Use in Evaluating PBL Parameterizations</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0227.1.
    <!--329-->
    <li>Hiris, ZA, and Gallus WA Jr., 2021: <a href="https://doi.org/10.3390/atmos12081019">On the Relationship of Cold Pool and Bulk Shear Magnitudes on Upscale Convective Growth in the Great Plains of the United States</a>. <i>Atmosphere,</i> <b>12,</b>, doi:10.3390/atmos12081019. 
    <!--328-->
    <li>Goldacker, N. A., and Parker, M. D., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0354.1">Low-Level Updraft Intensification in Response to Environmental Wind Profiles</a>. <i>J. Atmos. Sci.,</i> <b>78,</b> 2763-2781,doi:10.1175/JAS-D-20-0354.1. 
    <!--327-->
    <li>Lasher-Trapp, S., Jo, E., Allen, L. R., Engelsen, B. N., and Trapp, R. J., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0223.1">Entrainment in a Simulated Supercell Thunderstorm. Part I: The Evolution of Different Entrainment Mechanisms and Their Dilutive Effects</a>. <i>J. Atmos. Sci.,</i> <b>78,</b> 2725-2740, doi:10.1175/JAS-D-20-0223.1. 
    <!--326-->
    <li>Davies-Jones, R., and Markowski, P. M., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0020.1">Circulation around a Constrained Curve: An Alternative Analysis Tool for Diagnosing the Origins of Tornado Rotation in Numerical Supercell Simulations</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0020.1. 
    <!--325-->
    <li>Murdzek, S. S., Markowski, P. M., Richardson, Y. P., and Kumjian, M. R., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0069.1">Should Reversible Convective Inhibition be Used when Determining the Inflow Layer of a Convective Storm?</a>  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0069.1. 
    <!--324-->
    <li>Wu, F., and Lombardo, K., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0222.1">Precipitation Enhancement in Squall Lines Moving over Mountainous Coastal Regions</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0222.1. 
    <!--323-->
    <li>Peters, J. M., and Chavas, D. R., 2021: <A href="https://doi.org/10.1175/JAS-D-20-0351.1">Evaluating the conservation of energy variables in simulations of deep moist convection</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0351.1. 
    <!--322-->
    <li>Chen, J., and Chavas, D. R., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0037.1">Can Existing Theory Predict the Response of Tropical Cyclone Intensity to Idealized Landfall?</a> <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0037.1. 
    <!--321-->
    <li>Kumjian, M. R., Lombardo, K., and Loeffler, S., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0034.1">The Evolution of Hail Production in Simulated Supercell Storms</a>. <i>J. Atmos. Sci.,</i>, doi:10.1175/JAS-D-21-0034.1. 
    <!--320-->
    <li>Wang, A., Pan, Y., and Markowski, P. M., 2021: <a href="https://doi.org/10.1175/JAS-D-21-0033.1">The Influence of WENO Schemes on Large-Eddy Simulations of a Neutral Atmospheric Boundary Layer</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-21-0033.1.
    <!--319-->
    <li>Fu, S., Rotunno, R., Chen, J., Deng, X., and Xue, H., 2021:  <a href="https://doi.org/10.5194/acp-21-9289-2021">A large-eddy simulation study of deep-convection initiation through the collision of two sea-breeze fronts</a>. <i>Atmos. Chem. Phys.,</i> <b>21,</b> 9289–9308, doi:10.5194/acp-21-9289-2021.
    <!--318-->
    <li>Chandrakar, K. K., W. W. Grabowski, H. Morrison, and G. H. Bryan, 2021: <a href="https://doi.org/10.1175/JAS-D-20-0281.1">Impact of entrainment-mixing and turbulent fluctuations on droplet size distributions in a cumulus cloud: An investigation using Lagrangian microphysics with a sub-grid-scale model.</a>  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0281.1.
    <!--317-->
    <li>Marion, G. R., and Trapp, R. J., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0164.1">Controls of Quasi-Linear Convective System Tornado Intensity</a>. <i>Journal of the Atmospheric Sciences,</i> <b>78,</b> 1189-1205, doi:10.1175/JAS-D-20-0164.1. 
    <!--316-->
    <li>Mulholland, J. P., Peters, J. M., and Morrison, H., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0299.1">How Does Vertical Wind Shear Influence Entrainment in Squall Lines?</a> <i>Journal of the Atmospheric Sciences,</i> <b>78,</b> 1931-1946, doi:10.1175/JAS-D-20-0299.1. 
    <!--315-->
    <li>Fei, R., Wang, Y., and Li, Y., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0075.1">Contribution of Vertical Advection to Supergradient Wind in Tropical Cyclone Boundary Layer: A Numerical Study</a>. <i>Journal of the Atmospheric Sciences</i>, <b>78</b>, 1057-1073, doi:10.1175/JAS-D-20-0075.1. 
    <!--314-->
    <li>Katona, B., and P. Markowski, 2021: <a href="https://doi.org/10.1175/WAF-D-20-0136.1">Assessing the Influence of Complex Terrain on Severe Convective Environments in Northeastern Alabama</a>. <i>Weather and Forecasting,</i> <b>36,</b> 1003-1029, doi:10.1175/WAF-D-20-0136.1. 
    <!--313-->
    <li>Wang, Y., Li, Y., Xu, J., Tan, Z., and Lin, Y., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0393.1">The Intensity Dependence of Tropical Cyclone Intensification Rate in a Simplified Energetically Based Dynamical System Model.</a> <i>Journal of the Atmospheric Sciences,</i> <b>78,</b> 2033-2045, doi:10.1175/JAS-D-20-0393.1. 
    <!--312-->
    <li>Wang, S. and Smith, R.K., 2021: <a href="https://doi.org/10.1002/qj.4110">Upper-level trajectories in the prototype problem for tropical cyclone intensification</a>. <i>Q J R Meteorol Soc,</i> <b>147,</b> 2978-2987, doi:10.1002/qj.4110. 
    <!--311-->
    <li>Smith, R.K., Kilroy, G. and Montgomery, M.T., 2021: <a href="https://doi.org/10.1002/qj.4133">Tropical cyclone life cycle in a three-dimensional numerical simulation</a>. <i>Q J R Meteorol Soc.,</i> <b>147,</b> 3373-3393, doi:10.1002/qj.4133. 
    <!--310-->
    <li>Mulholland, J. P., Peters, J. M., and Morrison, H., 2021. <a href="https://doi.org/10.1029/2021GL093316">How does LCL height influence deep convective updraft width?</a> <i>Geophysical Research Letters,</i> <b>48,</b> e2021GL093316, doi:10.1029/2021GL093316. 
    <!--309-->
    <li>Thompson, C. F., and D. M. Schultz, 2021: <a href="https://doi.org/10.1029/2021GL092649">The release of inertial instability near an idealized zonal jet</a>. <i>Geophysical Research Letters,</i> <b>48,</b> e2021GL092649, doi:10.1029/2021GL092649.
    <!--308-->
    <li>Helms, C. N., and L. F. Bosart, 2021: <a href="https://doi.org/10.1175/MWR-D-20-0380.1">The Impact of a Midlevel Dry Airflow Layer on Deep Convection in the Pre-Gabrielle (2013) Tropical Disturbance on 4–5 September</a>. <i>Mon. Wea. Rev.,</i> <b>149,</b> 2695-2711. doi:10.1175/MWR-D-20-0380.1.
    <!--307-->
    <li>Richter, D. H., C. Wainwright, D. P. Stern, G. H. Bryan, and D. Chavas, 2021: <a href="https://doi.org/10.1175/JAS-D-20-0390.1">Potential low bias in high-wind drag coefficient inferred from dropsonde data in hurricanes</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0390.1. 
    <!--306-->
    <li>Schueth, A., C. Weiss, and J. M. L. Dahl, 2021: <a href="https://doi.org/10.1175/MWR-D-20-0251.1">Comparing observations and simulations of the streamwise vorticity current and the forward-flank convergence boundary in a supercell storm</a>.  <i>Mon. Wea. Rev.,</i> <b>149,</b> 1651-1671, doi:10.1175/MWR-D-20-0251.1.
    <!--305-->
    <li>Parker, M. D., 2021: <a href="https://doi.org/10.1175/MWR-D-20-0263.1">Self-Organization and Maintenance of Simulated Nocturnal Convective Systems from PECAN</a>. <i>Monthly Weather Review,</i> <b>149,</b> 999-1022, doi:10.1175/MWR-D-20-0263.1. 
    <!--304-->
    <li>Kilroy, G., 2021: <a href="https://doi.org/10.1002/qj.4011">Evolution of Convective Characteristics During Tropical Cyclogenesis</a>. <i>Q J R Meteorol Soc.,</i> <b>147,</b> 2103-2123, doi:10.1002/qj.4011. 
    <!--303-->
    <li>Wainwright, C., and D. Richter, 2021: <a href="https://doi.org/10.1007/s10546-020-00599-6">Investigating the Sensitivity of Marine Fog to Physical and Microphysical Processes Using Large-Eddy Simulation</a>. <i>Boundary-Layer Meteorol</i>, doi:s10546-020-00599-6. 
    <!--302-->
    <li>Fernando, H. J. S., Gultepe, I., Dorman, C., Pardyjak, E., Wang, Q., Hoch, S. W., Richter, D., Creegan, E., Gaberšek, S., Bullock, T., Hocut, C., Chang, R., Alappattu, D., Dimitrova, R., Flagg, D., Grachev, A., Krishnamurthy, R., Singh, D. K., Lozovatsky, I., Nagare, B., Sharma, A., Wagh, S., Wainwright, C., Wroblewski, M., Yamaguchi, R., Bardoel, S., Coppersmith, R. S., Chisholm, N., Gonzalez, E., Gunawardena, N., Hyde, O., Morrison, T., Olson, A., Perelet, A., Perrie, W., Wang, S., and Wauer, B., 2021: <a href="https://doi.org/10.1175/BAMS-D-19-0070.1">C-FOG: Life of Coastal Fog</a>. <i>Bulletin of the American Meteorological Society</i>, <b>102,</b> E244-E272, doi:10.1175/BAMS-D-19-0070.1. 
    <!--301-->
    <li>Peng, K., and J. Fang, 2021: <a href="https://doi.org/10.1029/2020JD033697">Effect of the initial vortex vertical structure on early development of an axisymmetric tropical cyclone</a>. <i>Journal of Geophysical Research: Atmospheres</i>, <b>126,</b> e2020JD033697, doi:10.1029/2020JD033697.
    <!--300-->
    <li>Hu, Q., and G. Limpert, 2021: <a href="https://doi.org/10.1002/qj.3982">Lift in the Vertical Shear of Southerly Jet Embedded in a Uniform Westerly Flow</a>. <i>Q J R Meteorol Soc.</i> <b>147,</b> 1584-1605, doi:10.1002/qj.3982. 
    <!--299-->
    <li>Lane, T. P., 2021: <a href="https://doi.org/10.1029/2020GL091025">Does lower‐stratospheric shear influence the mesoscale organization of convection?</a> <i>Geophysical Research Letters,</i> <b>48,</b> e2020GL091025, doi:10.1029/2020GL091025. 
    <!--298-->
    <li>Wade, A. R., and Parker, M. D., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0117.1">Dynamics of Simulated High-Shear Low-CAPE Supercells</a> <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-20-0117.1. 
    <!--297-->
    <li>Gowan, T. M., Steenburgh, W. J., and Minder, J. R., 2021: <a href="https://doi.org/10.1175/MWR-D-20-0253.1">Downstream Evolution and Coastal-to-Inland Transition of Landfalling Lake-Effect Systems</a>, <i>Monthly Weather Review,</i> doi:10.1175/MWR-D-20-0253.1. 
    <!--296-->
    <li>Morales, A., Posselt, D. J., and Morrison, H., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0142.1">Which combinations of environmental conditions and microphysical parameter values produce a given orographic precipitation distribution?</a>. <i>Journal of the Atmospheric Sciences,</i> </b>78,</b> 619-638, doi:10.1175/JAS-D-20-0142.1. 
    <!--295-->
    <li>Chavas, D. R., and Dawson, D. T., II., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0120.1">An idealized physical model for the severe convective storm environmental sounding</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-20-0120.1. 
    <!--294-->
    <li>Morrison, H., Peters, J. M., and Sherwood, S. C, 2021: <a href="https://doi.org/10.1175/JAS-D-20-0166.1">Comparing growth rates of simulated moist and dry convective thermals</a>. <i>Journal of the Atmospheric Sciences</i> doi:10.1175/JAS-D-20-0166.1. 
    <!--293-->
    <li>Wang, D., and Lin, Y., 2021: <a href="https://doi.org/10.1175/JAS-D-20-0192.1">Potential role of irreversible moist processes in modulating tropical cyclone surface wind structure</a>. <i>Journal of the Atmospheric Sciences,</i> doi:10.1175/JAS-D-20-0192.1. 
    <!--292-->
    <li>Peters, J.. M., Morrison, H., Zhang, G.. J., and Powell, S.. W., 2021. <a href="https://doi.org/10.1029/2020MS002282">Improving the physical basis for updraft dynamics in deep convection parameterizations</a>. <i>Journal of Advances in Modeling Earth Systems,</i> <b>13,</b> e2020MS002282, doi:10.1029/2020MS002282. 
    <!--291-->
    <li>Bickle, M.E., Marsham, J.H., Ross, A.N., Rowell, D.P., Parker, D.J. and Taylor, C.M., 2021: <a href="https://doi.org/10.1002/qj.3955">Understanding mechanisms for trends in Sahelian squall lines: Roles of thermodynamics and shear</a>. <i>Q J R Meteorol Soc.,</i> doi:10.1002/qj.3955. 
    <!--290-->
    <li>Montgomery, M. T., and Persing, J., 2021: <a href="https://doi.org/10.1175/JAS-D-19-0258.1">Does Balance Dynamics Well Capture the Secondary Circulation and Spinup of a Simulated Hurricane?</a> <i>Journal of the Atmospheric Sciences,</i> <b>78,</b> 75-95, doi:10.1175/JAS-D-19-0258.1.
    <!--289-->
    <li>Alland, J. J., B. H. Tang, K. L. Corbosiero, and G. H. Bryan, 2021: <a href="https://doi.org/10.1175/JAS-D-20-0055.1">Combined effects of midlevel dry air and vertical wind shear on tropical cyclone development. Part II: Radial ventilation.</a>.  <i>J. Atmos. Sci.</i>, doi:10.1175/JAS-D-20-0055.1.
    <!--288-->
    <li>Alland, J. J., B. H. Tang, K. L. Corbosiero, and G. H. Bryan, 2021: <a href="https://doi.org/10.1175/JAS-D-20-0054.1">Combined effects of midlevel dry air and vertical wind shear on tropical cyclone development. Part I: Downdraft ventilation.</a>.  <i>J. Atmos. Sci.</i>, doi:10.1175/JAS-D-20-0054.1.
    <!--287-->
    <li>Rousseau-Rizzi, R., R. Rotunno, and G. H. Bryan, 2021:    <a href="https://doi.org/10.1175/JAS-D-20-0140.1">A thermodynamic perspective on steady-state tropical cyclones</a>.  <i>J. Atmos. Sci.</i>, doi:10.1175/JAS-D-20-0140.1.
    <!--286-->
    <li>Haghi, K. R., and D. R. Durran, 2021: <a href="https://doi.org/10.1175/JAS-D-20-0181.1">On the Dynamics of Atmospheric Bores</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0181.1.
    <!--285-->
    <li>Hartigan, J., R. A. Warren, J. S. Soderholm, and H. Richter, 2021: <a href="https://doi.org/10.1175/MWR-D-20-0069.1">Simulated Changes in Storm Morphology Associated with a Sea-Breeze Air Mass</a>.  <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-20-0069.1.
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>

  <details>
  <summary>2020</summary>
  <ul>
    <!--284-->
    <li>Drueke, S., D. J. Kirshbaum, and P. Kollias, 2020: <a href="https://doi.org/10.5194/acp-20-13217-2020">Environmental sensitivities of shallow-cumulus dilution – Part 1: Selected thermodynamic conditions</a>.  <i>Atmospheric Chemistry and Physics,</i> <b>20,</b> 13217–13239, doi:10.5194/acp-20-13217-2020. 
    <!--283-->
    <li>Peters, J. M., H. Morrison, C. J. Nowotarski, J. P. Mulholland, and R. L. Thompson, 2020: <a href="https://doi.org/10.1175/JAS-D-20-0103.1">A Formula for the Maximum Vertical Velocity in Supercell Updrafts</a>. <i>J. Atmos. Sci.,</i> <b>77,</b> 3747–3757, doi:10.1175/JAS-D-20-0103.1.
    <!--282-->
    <li>Martinez, J., Nam, C. C., and Bell, M. M., 2020: <a href="https://doi.org/10.1029/2020JD033324">On the contributions of incipient vortex circulation and environmental moisture to tropical cyclone expansion</a>. <i>Journal of Geophysical Research: Atmospheres,</i> <b>125,</b> e2020JD033324, doi:10.1029/2020JD033324.
    <!--281-->
    <li>Li, L., Chakraborty, P., 2020: <a href="https://doi.org/10.1038/s41586-020-2867-7">Slower decay of landfalling hurricanes in a warming world</a>. <i>Nature,</i> <b>587,</b> 230–234, doi:10.1038/s41586-020-2867-7. 
    <!--280-->
    <li>Ramsay, H. A., Singh, M. S., and Chavas, D. R., 2020:  <a href="https://doi.org/10.1029/2020MS002086">Response of tropical cyclone formation and intensification rates to climate warming in idealized simulations</a>. <i>Journal of Advances in Modeling Earth Systems,</i> <b>12,</b> e2020MS002086. doi:10.1029/2020MS002086. 
    <!--279-->
    <li>Naylor, J., 2020: <a href="https://doi.org/10.3390/atmos11070707">Idealized Simulations of City-Storm Interactions in a Two-Dimensional Framework</a>.  <i>Atmosphere,</i> 11, 707, doi:10.3390/atmos11070707.
    <!--278-->
    <li>Fischer, J., and J. M. L. Dahl, 2020: <a href="https://doi.org/10.1175/JAS-D-20-0126.1">The Relative Importance of Updraft and Cold Pool Characteristics on Supercell Tornadogenesis in Highly Idealized Simulations</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0126.1.
    <!--277-->
    <li>Adams-Selin, R. D., 2020: <a href="https://doi.org/10.1175/JAS-D-19-0347.1">Sensitivity of MCS Low-Frequency Gravity Waves to Microphysical Variations</a>.  <i>J. Atmos. Sci.,</i> <b>77,</b> 3461–3477, doi:10.1175/JAS-D-19-0347.1.
    <!--276-->
    <li>Wang, A., Y. Pan, and P. M. Markowski, 2020: <a href="https://doi.org/10.1175/MWR-D-20-0031.1">The influence of turbulence memory on idealized tornado simulations</a>.  <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-20-0031.1.
    <!--275-->
    <li>Lombardo, K., 2020:  <a href="https://doi.org/10.1175/JAS-D-20-0044.1">Squall Line Response to Coastal Mid-Atlantic Thermodynamic Heterogeneities</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0044.1.
    <!--274-->
    <li>Li, Y., Y. Wang, Y. Lin, and R. Fei, 2020: <a href="https://doi.org/10.1175/MWR-D-20-0141.1">Dependence of superintensity of tropical cyclones on SST in axisymmetric numerical simulations</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-20-0141.1.
    <!--273-->
    <li>Peters, J. M., H. Morrison, A. C. Varble, W. M. Hannah, and S. E. Giangrande, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0244.1">Thermal Chains and Entrainment in Cumulus Updrafts. Part II: Analysis of Idealized Simulations</a>. <i>J. Atmos. Sci.,</i> <b>77,</b> 3661–3681, 10.1175/JAS-D-19-0244.1.
    <!--272-->
    <li>Morrison, H., Peters, J. M., Varble, A. C., Hannah, W. M., and Giangrande, S. E., 2020: <a href="https://doi.org//10.1175/JAS-D-19-0243.1">Thermal Chains and Entrainment in Cumulus Updrafts. Part I: Theoretical Description,</a> <i>Journal of the Atmospheric Sciences,</i> <b>77,</b> 3637-3660, doi:/10.1175/JAS-D-19-0243.1. 
    <!--271-->
    <li>Tao, D., R. Rotunno, and M. Bell, 2020: <a href="https://doi.org/10.1175/JAS-D-20-0057.1">Lilly’s Model for Steady-State Tropical Cyclone Intensity and Structure</a>. <i>J. Atmos. Sci.,</i> <b>77,</b> 3701–3720, doi:10.1175/JAS-D-20-0057.1.
    <!--270-->
    <li>Reif, D. W., H. B. Bluestein, T. M. Weckwerth, Z. B. Wienhoff, and M. B. Chasteen, 2020: <a href="https://doi.org/10.1175/JAS-D-20-0028.1">Estimating the Maximum Vertical Velocity at the Leading Edge of a Density Current</a>.  <i>J. Atmos. Sci.,</i> <b>77,</b> 3683–3700, doi:10.1175/JAS-D-20-0028.1.
    <!--269-->
    <li>Boyer, C. H., and J. M. L. Dahl, 2020:  <a href="https://doi.org/10.1175/MWR-D-20-0082.1">The Mechanisms Responsible for Large Near-Surface Vertical Vorticity within Simulated Supercells and Quasi-Linear Storms</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-20-0082.1.
    <!--268-->
    <li>Flournoy, M. D., M. C. Coniglio, E. N. Rasmussen, J. C. Furtado, and B. E. Coffer, 2020:  <a href="https://doi.org/10.1175/MWR-D-20-0147.1">Modes of storm-scale variability and tornado potential in VORTEX2 near- and far-field tornadic environments</a>. <i>Mon. Wea. Rev.,</i>, doi:10.1175/MWR-D-20-0147.1.
    <!--267-->
    <li>Adams-Selin, R. D., 2020: <a href="https://doi.org/10.1175/JAS-D-19-0250.1">Impact of Convectively Generated Low-Frequency Gravity Waves on Evolution of Mesoscale Convective Systems</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0250.1.
    <!--266-->
    <li>Kirshbaum, D. J., 2020:  <a href="https://doi.org/10.1175/JAS-D-20-0035.1">Numerical simulations of orographic convection across multiple grey zones</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-20-0035.1.
    <!--265-->
    <li>Dahl, J. M. L., 2020: <a href="https://doi.org/10.1175/MWR-D-20-0080.1">Near-Surface Vortex Formation in Supercells from the Perspective of Vortex Patch Dynamics</a>. <i>Mon. Wea. Rev.,</i> <b>148,</b> 3533-3547, doi:10.1175/MWR-D-20-0080.1.
    <!--264-->
    <li>Nowotarski, C. J., J. M. Peters, and J. P. Mulholland, 2020: <a href="https://doi.org/10.1175/MWR-D-20-0013.1">Evaluating the Effective Inflow Layer of Simulated Supercell Updrafts</a>. <i>Mon. Wea. Rev.,</i> 148, 3507-3532, doi:10.1175/MWR-D-20-0013.1.
    <!--263-->
    <li>Steinkruger, D., P. Markowski, and G. Young, 2020:  <a href="https://doi.org/10.1175/WAF-D-19-0249.1">An Artificially Intelligent System for the Automated Issuance of Tornado Warnings in Simulated Convective Storms</a>. <i>Wea. Forecasting,</i> doi:10.1175/WAF-D-19-0249.1.
    <!--262-->
    <li>Nystrom, R. G., R. Rotunno, C. A. Davis, and F. Zhang, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0357.1">Consistent impacts of surface enthalpy and drag coefficient uncertainty between an analytical model and simulated tropical cyclone maximum intensity and storm structure</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0357.1.
    <!--261-->
    <li>Kumjian, M. R., and K. Lombardo, 2020: <a href="https://doi.org/10.1175/JAS-D-20-0016.1">A Hail Growth Trajectory Model for Exploring the Environmental Controls on Hail Size: Model Physics and Idealized Tests</a>. <i>J. Atmos. Sci.,</i> <b>77,</b> 2765-2791, doi:10.1175/JAS-D-20-0016.1. 
    <!--260-->
    <li>Hitchcock, S. M., and R. S. Schumacher, 2020: <a href="https://doi.org/10.1175/MWR-D-19-0246.1">Analysis of Back-Building Convection in Simulations with a Strong Low-Level Stable Layer</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-19-0246.1.
    <!--259-->
    <li>Li, Y., Y. Wang, and Y. Lin, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0350.1">How Much Does the Upward Advection of the Supergradient Component of Boundary Layer Wind Contribute to Tropical Cyclone Intensification and Maximum Intensity?</a> <i>J. Atmos. Sci.,</i> <b>77,</b> 2649-2664, doi:10.1175/JAS-D-19-0350.1.
    <!--258-->
    <li>Ryglicki, D. R., D. Hodyss, and G. Rainwater, 2020:  <a href="https://doi.org/10.1175/JAS-D-20-0030.1">The Tropical Cyclone as a Divergent Source in a Background Flow</a>. <i>J. Atmos. Sci.,</i>, doi:10.1175/JAS-D-20-0030.1.
    <!--257-->
    <li>Mulholland, J. P., S. W. Nesbitt, R. J. Trapp, and J. M. Peters, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0190.1">The influence of terrain on the convective environment and associated convective morphology from an idealized modeling prospective</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0190.1.
    <!--256-->
    <li>Wing, A. A., Stauffer, C. L., Becker, T., Reed, K. A., and Coauthors, 2020: <a href="https://doi.org/10.1029/2020MS002138">Clouds and Convective Self‐Aggregation in a Multi‐Model Ensemble of Radiative‐Convective Equilibrium Simulations</a>. <i>Journal of Advances in Modeling Earth Systems,</i> <b>12,</b> e2020MS002138, doi:10.1029/2020MS002138.
    <!--255-->
    <li>Montgomery, M.T., Kilroy, G., Smith, R.K. and Črnivec, N., 2020: <a href="https://doi.org/10.1002/qj.3837">Contribution of mean and eddy momentum processes to tropical cyclone intensification</a>. <i>Q J R Meteorol Soc.,</i> doi:10.1002/qj.3837. 
    <!--254-->
    <li>Wang, S., Smith, R.K. and Montgomery, M.T., 2020: <a href="https://doi.org/10.1002/qj.3856">Upper‐tropospheric inflow layers in tropical cyclones</a>. <i>Q J R Meteorol Soc.,</i> doi:10.1002/qj.3856.
    <!--253-->
    <li>Takemi, T., and S. Yamasaki, 2020. <a href="https://doi.org/10.3390/atmos11040411">Sensitivity of the Intensity and Structure of Tropical Cyclones to Tropospheric Stability Conditions</a>. <i>Atmosphere</i>, 11, 411, doi:10.3390/atmos11040411.
    <!--252-->
    <li>Peters, J.M., C.J. Nowotarski, J.P. Mulholland, and R.L. Thompson, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0355.1">The influences of effective inflow layer streamwise vorticity and storm-relative flow on supercell updraft properties</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0355.1. 
    <!--251-->
    <li>O’Neill, M.E., and D.R. Chavas, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0330.1">Inertial waves in axisymmetric tropical cyclones</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0330.1. 
    <!--250-->
    <li>Markowski, P.M., 2020: <a href="https://doi.org/10.1175/MWR-D-20-0076.1">What is the intrinsic predictability of tornadic supercell thunderstorms?</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-20-0076.1. 
    <!--249-->
    <li>Chen, J., and D.R. Chavas, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0320.1">The Transient Responses of An Axisymmetric Tropical Cyclone to Instantaneous Surface Roughening and Drying</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0320.1. 
    <!--248-->
    <li>Wang, D., and Y. Lin, 2020: <a href="https://doi.org/10.1175/JAS-D-19-0229.1">Size and Structure of Dry and Moist Reversible Tropical Cyclones</a>. <i>J. Atmos. Sci.,</i> <b>77,</b> 2091-2114, doi:10.1175/JAS-D-19-0229.1.  
    <!--247-->
    <li>Stern, D.P., J.D. Kepert, G.H. Bryan, and J.D. Doyle, 2020: <a href="https://dx.doi.org/10.1175/JAS-D-19-0191.1">Understanding Atypical Mid-Level Wind Speed Maxima in Hurricane Eyewalls</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0191.1
    <!--246-->
    <li>Tao, D., Bell, M., Rotunno, R., and van Leeuwen, P. J.,  2020:  <a href="https://dx.doi.org/10.1029/2019GL085980">Why do the maximum intensities in modeled tropical cyclones vary under the same environmental conditions?</a>.  <i>Geophysical Research Letters,</i> doi:10.1029/2019GL085980.   
    <!--245-->
    <li>Warren, R. A., Singh, M. S., and Jakob, C., 2020: <a href="https://dx.doi.org/10.1029/2019MS001734">Simulations of radiative‐convective‐dynamical equilibrium</a>. <i>Journal of Advances in Modeling Earth Systems,</i> doi:10.1029/2019MS001734. 
    <!--244-->
    <li>Kapoor, A., Ouakka, S., Arwade, S. R., Lundquist, J. K., Lackner, M. A., Myers, A. T., Worsnop, R. P., and Bryan, G. H., 2020: <a href="https://dx.doi.org/10.5194/wes-5-89-2020">Hurricane eyewall winds and structural response of wind turbines</a>.  <i>Wind Energ. Sci.,</i>  <b>5,</b>  89-104, doi:10.5194/wes-5-89-2020.
    <!--243-->
    <li>Rotunno, R., and G.H. Bryan, 2020: <a href="https://dx.doi.org/10.1175/JAS-D-19-0142.1">Numerical Simulations of Two-Layer Flow past Topography. Part II: Lee Vortices</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0142.1. 
    <!--242-->
    <li>Peters, J.M., C.J. Nowotarski, and G.L. Mullendore, 2020: <a href="https://dx.doi.org/10.1175/JAS-D-19-0316.1">Are supercells resistant to entrainment because of their rotation?</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0316.1. 
    <!--241-->
    <li>Parker, M.D., B.S. Borchardt, R.L. Miller, and C.L. Ziegler, 2020: <a href="https://dx.doi.org/10.1175/MWR-D-19-0072.1">Simulated Evolution and Severe Wind Production by the 25–26 June 2015 Nocturnal MCS from PECAN</a>. <i>Mon. Wea. Rev.,</i> <b>148,</b> 183–209, doi:10.1175/MWR-D-19-0072.1. 
    <!--240-->
    <li>Tang, SL, and Kirshbaum, DJ, 2020. <a href="https://dx.doi.org/10.1002/qj.3726">On the sensitivity of deep‐convection initiation to horizontal grid resolution</a>. <i>Q J R Meteorol Soc.,</i> doi:10.1002/qj.3726. 
    <!--239-->
    <li>Schecter, D.A. and K. Menelaou, 2020: <a href="https://dx.doi.org/10.1175/JAS-D-19-0074.1">Development of a Misaligned Tropical Cyclone</a>. <i>J. Atmos. Sci.,</i> <b>77,</b> 79–111, doi:10.1175/JAS-D-19-0074.1.
    <!--236-->
    <li>Kilroy, G. , Smith, R. K., and Montgomery, M. T., 2020:  <a href="https://dx.doi.org/10.1002/qj.3701">An idealized numerical study of tropical cyclogenesis and evolution at the Equator.</a>  <i>Q J R Meteorol Soc.</i> doi:10.1002/qj.3701. 
    <!--235-->
    <li>Cione, J.J., G.H. Bryan, R. Dobosy, J.A. Zhang, G. de Boer, A. Aksoy, J.B. Wadler, E.A. Kalina, B.A. Dahl, K. Ryan, Jonathan, Neuhaus, E. Dumas, F.D. Marks, A.M. Farber, T. Hock, and X. Chen, 2020: <a href="https://dx.doi.org/10.1175/BAMS-D-19-0169.1">Eye of the Storm: Observing Hurricanes with a Small Unmanned Aircraft System</a>. <i>Bull. Amer. Meteor. Soc.,</i> doi:10.1175/BAMS-D-19-0169.1.
    <!--234-->
    <li>Thanh, NT, Cuong, HD, Hien, NX, Kieu, C., 2020:  <a href="https://dx.doi.org/10.1002/joc.6348">Relationship between sea surface temperature and the maximum intensity of tropical cyclones affecting Vietnam's coastline</a>.  <i>Int J Climatol.,</i> doi:10.1002/joc.6348. 
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>

  <details>
  <summary>2019</summary>
  <ul>
    <!--238-->
    <li>Tang, X., Cai, Q., Fang, J., and Tan, Z.‐M.,  2019:  <a href="https://dx.doi.org/10.1029/2019JD031454">Land–Sea Contrast in the Diurnal Variation of Precipitation from Landfalling Tropical Cyclones</a>.  <i>Journal of Geophysical Research: Atmospheres,</i> doi:10.1029/2019JD031454. 
    <!--237-->
    <li>Singh, M. S., Warren, R. A., and Jakob, C., 2019:  <a href="https://dx.doi.org/10.1029/2019MS001686">A steady‐state model for the relationship between humidity, instability, and precipitation in the tropics</a>.  <i>Journal of Advances in Modeling Earth Systems,</i> doi:10.1029/2019MS001686. 
    <!--233-->
    <li>Wang, S, Smith, RK, 2019: <a href="https://dx.doi.org/10.1002/qj.3656">Consequences of regularizing the Sawyer–Eliassen equation in balance models for tropical cyclone behaviour</a>. <i>Q J R Meteorol Soc</i>, doi:10.1002/qj.3656. 
    <!--232-->
    <li>Davenport, C.E., C.L. Ziegler, and M.I. Biggerstaff, 2019: <a href="https://dx.doi.org/10.1175/MWR-D-18-0447.1">Creating a More Realistic Idealized Supercell Thunderstorm Evolution via Incorporation of Base-State Environmental Variability</a>.  <i>Mon. Wea. Rev.,</i> <b>147,</b> 4177–4198, doi:10.1175/MWR-D-18-0447.1.
    <!--231-->
    <li>Drueke, S., Kirshbaum, D. J., and Kollias, P., 2019: <a href="https://dx.doi.org/10.1029/2019JD030889">Evaluation of shallow‐cumulus entrainment rate retrievals using large‐eddy simulation</a>.  <i>Journal of Geophysical Research: Atmospheres,</i> <b>124</b> doi:10.1029/2019JD030889. 
    <!--230-->
    <li>Dawson, D.T., B. Roberts, and M. Xue, 2019: <a href="https://dx.doi.org/10.1175/MWR-D-18-0462.1">A method to control the environmental wind profile in idealized simulations of deep convection with surface friction</a>. <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-18-0462.1. 
    <!--229-->
    <li>Shi, X., R.M. Enriquez, R.L. Street, G.H. Bryan, and F.K. Chow, 2019: <a href="https://dx.doi.org/10.1175/JAS-D-18-0375.1">An Implicit Algebraic Turbulence Closure Scheme for Atmospheric Boundary Layer Simulation</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0375.1. 
    <!--228-->
    <li>Mallinson, H.M. and S.G. Lasher-Trapp, 2019: <a href="https://dx.doi.org/10.1175/MWR-D-18-0382.1">An Investigation of Hydrometeor Latent Cooling upon Convective Cold Pool Formation, Sustainment, and Properties</a>.  <i>Mon. Wea. Rev.,</i> <b>147,</b> 3205-3222, doi:10.1175/MWR-D-18-0382.1.  
    <!--227-->
    <li>Hutson, A., C. Weiss, and G. Bryan, 2019:  <a href="https://dx.doi.org/10.1175/MWR-D-18-0439.1">Using the Translation Speed and Vertical Structure of Gust Fronts to Infer Buoyancy Deficits within Thunderstorm Outflow</a>.  <i>Mon. Wea. Rev.,</i> 10.1175/MWR-D-18-0439.1.
    <!--226-->
    <li>Li, Y., Y. Wang, and Y. Lin, 2019:  <a href="https://dx.doi.org/10.1175/JAS-D-19-0076.1">Revisiting the dynamics of eyewall contraction of tropical cyclones</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0076.1.
    <!--225-->
    <li>Peters, J.M., C.J. Nowotarski, and H. Morrison, 2019: <a href="https://dx.doi.org/10.1175/JAS-D-19-0096.1">The role of vertical wind shear in modulating maximum supercell updraft velocities</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-19-0096.1.
    <!--224-->
    <li>Kirshbaum, D. J., and D. N. Straub, 2019: <a href="https://dx.doi.org/10.1002/qj.3609">Linear theory of shallow convection in deep, vertically sheared atmospheres</a>. <i>Q J R Meteorol Soc.,</i> doi:10.1002/qj.3609.
    <!--223-->
    <li>Swenson, S., B. Argrow, E. Frew, S. Borenstein, and J. Keeler, 2019: <a href="https://dx.doi.org/10.3390/s19092149">Development and Deployment of Air-Launched Drifters from Small UAS</a>. <i>Sensors,</i> <b>19,</b> doi:10.3390/s19092149.
    <!--222-->
    <li>Yao, D., Z. Meng, and M. Xue, 2019: <a href="https://dx.doi.org/10.3390/atmos10050236">Genesis, Maintenance and Demise of a Simulated Tornado and the Evolution of Its Preceding Descending Reflectivity Core (DRC)</a>.  <i>Atmosphere,</i> <b>10,</b> 236, doi:10.3390/atmos10050236.
    <!--221-->
    <li>Ryglicki, D.R., J.D. Doyle, D. Hodyss, J.H. Cossuth, Y. Jin, K.C. Viner, and J.M. Schmidt, 2019: <a href="https://dx.doi.org/10.1175/MWR-D-18-0370.1">The Unexpected Rapid Intensification of Tropical Cyclones In Moderate Vertical Wind Shear. Part III: Outflow-Environment Interaction</a>.  <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-18-0370.1. 
    <!--220-->
    <li>Zhou, W., and S.-P. Xie, 2019:  <a href="https://dx.doi.org/10.1175/JAS-D-18-0330.1">A Conceptual Spectral Plume Model for Understanding Tropical Temperature Profile and Convective Updraft Velocities</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0330.1</a>.
    <!--219-->
    <li>Rydbeck, A.V., T.G. Jensen, and M.R. Igel, 2019: <a href="https://dx.doi.org/10.1175/JAS-D-18-0303.1">Idealized Modeling of the Atmospheric Boundary Layer Response to SST Forcing in the Western Indian Ocean</a>. <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0303.1. 
    <!--218-->
    <li>Peng, K., R. Rotunno, G.H. Bryan, and J. Fang, 2019: <a href="https://dx.doi.org/10.1175/JAS-D-18-0264.1">Evolution of an Axisymmetric Tropical Cyclone before Reaching Slantwise Moist Neutrality</a>.  <i>J. Atmos. Sci.,</i> <b>76,</b> 1865-1884, doi:10.1175/JAS-D-18-0264.1. 
    <!--217-->
    <li>Montgomery, MT, Persing, J, Smith, RK, 2019:  <a href="https://dx.doi.org/10.1002/qj.3479">On the hypothesized outflow control of tropical cyclone intensification</a>.  <i>Q J R Meteorol Soc.,</i> 1-14, doi:10.1002/qj.3479. 
    <!--216-->
    <li>Colbert, M., D.J. Stensrud, P.M. Markowski, and Y.P. Richardson, 2019: <a href="https://dx.doi.org/10.1175/WAF-D-18-0175.1">Processes Associated with Convection Initiation in the North American Mesoscale Forecast System – Version 3 (NAMv3)</a>.  <i>Wea. Forecasting,</i> doi:10.1175/WAF-D-18-0175.1. 
    <!--215-->
    <li>Sherburn, K.D., and M.D. Parker, 2019: <a href="https://dx.doi.org/10.1175/MWR-D-18-0246.1">The development of severe vortices within simulated high-shear, low-CAPE convection</a>.  <i>Mon. Wea. Rev.,</i> doi:10.1175/MWR-D-18-0246.1. 
    <!--214-->
    <li>Peters, J.M., W. Hannah, and H. Morrison, 2019: <a href="https://dx.doi.org/10.1175/JAS-D-18-0296.1">The influence of vertical wind shear on moist thermals</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0296.1. 
    <!--213-->
    <li>Rousseau-Rizzi, R., and K. Emanuel, 2019: <a href="https://dx.doi.org/10.1175/JAS-D-18-0238.1">An evaluation of hurricane superintensity in axisymmetric numerical models</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0238.1. 
    <!--212-->
    <li>Brown, M., and C. Nowotarski, 2019:  <a href="https://dx.doi.org/10.1175/JAS-D-18-0216.1">The influence of lifting condensation level on low–level outflow and rotation in simulated supercell thunderstorms</a>.  <i>J.  Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0216.1. 
    <!--211-->
    <li>Morales, A., D. Posselt, H. Morrison, and F. He, 2019:  <a href="https://dx.doi.org/10.1175/JAS-D-18-0301.1">Assessing the influence of microphysical and environmental parameter perturbations on orographic precipitation</a>.  <i>J. Atmos. Sci.,</i> doi:10.1175/JAS-D-18-0301.1.  
    <!--210-->
    <li>Molinari, J., M. Rosenmayer, D. Vollaro, and S.D. Ditchek, 2019: <a href="https://dx.doi.org/10.1175/JAMC-D-18-0148.1">Turbulence variations in the upper troposphere in tropical cyclones from NOAA G-IV flight-level vertical acceleration data</a>.  <i>J. Appl. Meteor. Climatol.,</i> <b>58,</b> 569-583, doi:10.1175/JAMC-D-18-0148.1.
    <!--209-->
    <li>Tao, J., C. Wang, N. V. Chawla, H. Guo, G. Sever, and S. H. Kim, 2019: <a href="https://dx.doi.org/10.1109/TVCG.2018.2864808">Exploring Time-Varying Multivariate Volume Data Using Matrix of Isosurface Similarity Maps</a>.  <i>IEEE Transactions on Visualization and Computer Graphics</i> <b>25,</b> 1236-1245, doi:10.1109/TVCG.2018.2864808. 
    <!--208-->
    <li>Steenkamp, S. C., G. Kilroy, and R. K. Smith, 2019:  <a href="https://dx.doi.org/10.1002/qj.3529">Tropical cyclogenesis at and near the Equator</a>.  <i>Q J R Meteorol Soc.,</i> doi:10.1002/qj.3529. 
    <!--207-->
    <li>Shi, X., F. K. Chow, R. L. Street, and G. H. Bryan, 2019: <a href="https://dx.doi.org/10.1029/2018MS001446">Key elements of turbulence closures for simulating deep convection at kilometer‐scale resolution</a>. <i>Journal of Advances in Modeling Earth Systems,</i> <b>11,</b> doi:10.1029/2018MS001446. 
    <!--206-->
    <li>Bai, L., Z. Meng, Y. Huang, Y. Zhang, S. Niu, and R. Su, 2019: <a href="https://dx.doi.org/10.1029/2018JD029832">Convection Initiation Resulting from the Interaction between a Quasi‐Stationary Dryline and Intersecting Gust Fronts: A Case Study</a>. <i>J. Geophys. Res. Atm.,</i> doi:10.1029/2018JD029832. 
    <!--205-->
    <li>Cai, Q., and X. Tang, 2019. <a href="http://dx.doi.org/10.1029/2018JD029107">Effect of the eyewall cold pool on the inner rainband of a tropical cyclone</a>. <i>J. of Geophys. Res. Atm.,</i> doi:/10.1029/2018JD029107. 
    <!--204-->
    <li>Marion, G. R., and R. J. Trapp, 2019.  <a href="http://dx.doi.org/10.1029/2018JD029055">The dynamical coupling of convective updrafts, downdrafts, and cold pools in simulated supercell thunderstorms</a>.  <i>J. of Geophys. Res. Atm.,</i> <b>124,</b> 664-683, doi:10.1029/2018JD029055.
    <!--203-->
    <li>Wang, Y. C., C. Davis, and Y. Huang, 2019: <a href="http://dx.doi.org/10.1175/JAS-D-18-0219.1"> Dynamics of lower-tropospheric vorticity in idealized simulations of tropical cyclone formation</a>. <i>J. Atmos. Sci.,</i> <b>76,</b> 707-727, doi:10.1175/JAS-D-18-0097.1. 
    <!--202-->
    <li>Duran, P., and J. Molinari, 2019: <a href="http://dx.doi.org/10.1175/JAS-D-18-0097.1">Tropopause evolution in a rapidly intensifying tropical cyclone: A static stability budget analysis in an idealized, axisymmetric framework.</a>  <i>J. Atmos. Sci.,</i> <b>76,</b> 209-229, doi.org/10.1175/JAS-D-18-0097.1. 
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>

  <details>
  <summary>2018</summary>
  <ul>
    <!--201-->
    <li>Schecter, D. A., 2018:  <a href="https://dx.doi.org/10.1080/16000870.2018.1525245">On the instabilities of tropical cyclones generated by cloud resolving models</a>. <i>Tellus A: Dynamic Meteorology and Oceanography,</i> <b>70,</b> 1-30, doi:10.1080/16000870.2018.1525245.
    <!--200-->
    <li>Belik, P., B. Dahl, D. Dokken, C. K. Potvin, K. Scholz, and M. Shvartsman, 2018:  <a href="http://dx.doi.org/10.3934/Math.2018.3.365">Possible implications of self-similarity for tornadogenesis and maintenance</a>.  <i>AIMS Mathematics,</i> <b>3,</b> 365-390, doi:10.3934/Math.2018.3.365. 
    <!--199-->
    <li>Oreskovic, C., L. G. Orf, and E. Savory, 2018:  <a href="http://dx.doi.org/10.1016/j.jweia.2018.07.020">A parametric study of downbursts using a full-scale cooling source model</a>.  <i>J. Wind Engineering and Industrial Aerodynamics,</i> <b>180,</b> 168-181, doi:10.1016/j.jweia.2018.07.020. 
    <!--198-->
    <li>Trapp, R. J., G. R. Marion, and S. W. Nesbitt, 2018: <a href="http://dx.doi.org/10.1175/JAS-D-18-0276.1">Reply to "Comments on 'The Regulation of Tornado Intensity by Updraft Width'"</a>. <i>J. Atmos. Sci.,</i> <b>75,</b> 4057–4061, doi:10.1175/JAS-D-18-0276.1. 
    <!--197-->
    <li>Coffer, B. E. and P. M. Markowski, 2018: <a href="http://dx.doi.org/10.1175/JAS-D-18-0170.1">Comments on "The Regulation of Tornado Intensity by Updraft Width"</a>. <i>J. Atmos. Sci.,</i> <b>75,</b> 4049–4056, doi:10.1175/JAS-D-18-0170.1. 
    <!--196-->
    <li>Lasher-Trapp, S., S. Kumar, D.H. Moser, A.M. Blyth, J.R. French, R.C. Jackson, D.C. Leon, and D.M. Plummer, 2018:  <a href="http://dx.doi.org/10.1175/JAMC-D-18-0041.1">On different microphysical pathways to convective rainfall</a>.  <i>J. Appl. Meteor. Climatol.,</i> <b>57,</b> 2399-2417, doi:10.1175/JAMC-D-18-0041.1.
    <!--195-->
    <li>Miyamoto, Y., D.S. Nolan, and N. Sugimoto, 2018: <a href="http://dx.doi.org/10.1175/JAS-D-18-0042.1">A dynamical mechanism for secondary eyewall formation in tropical cyclones</a>.  <i>J. Atmos. Sci.,</i> <b>57,</b> 3965-3986, doi:10.1175/JAS-D-18-0042.1. 
    <!--194-->
    <li>Ryglicki, D.R., J.D. Doyle, Y. Jin, D. Hodyss, and J. Cossuth, 2018:  <a href="http://dx.doi.org/10.1175/MWR-D-18-0021.1">The unexpected rapid intensification of tropical cyclones in moderate vertical wind shear. Part II: Vortex tilt</a>.  <i>Mon. Wea. Rev.,</i> <b>146,</b> 3801-3825, doi:10.1175/MWR-D-18-0021.1. 
    <!--193-->
    <li>Stern, D.P. and G.H. Bryan, 2018: <a href="https://doi.org/10.1175/MWR-D-18-0041.1">Using simulated dropsondes to understand extreme updrafts and wind speeds in tropical cyclones</a>.  <i>Mon. Wea. Rev.,</i> <b>146,</b> 3901–3925, doi:10.1175/MWR-D-18-0041.1. 
    <!--192-->
    <li>Williams, G.J., 2018:  <a href="https://doi.org/10.1007/s00703-018-0616-3">The effects of ice microphysics on the inner core thermal structure of the hurricane boundary layer</a>.  <i>Meteorol. Atmos. Phys.,</i> doi:10.1007/s00703-018-0616-3. 
    <!--191-->
    <li>Moser, D.H. and S. Lasher-Trapp, 2018: <a href="https://doi.org/10.1175/JAMC-D-17-0363.1">Cloud spacing effects upon entrainment and rainfall along a convective line</a>. <i>J. Appl. Meteor. Climatol.,</i> <b>57,</b> 1865-1882, doi:10.1175/JAMC-D-17-0363.1. 
    <!--190-->
    <li>Dahl, N.A. and D.S. Nolan, 2018: <a href="https://doi.org/10.1175/MWR-D-17-0333.1">Using high-resolution simulations to quantify errors in radar estimates of tornado intensity</a>.  <i>Mon. Wea. Rev.,</i> <b>146,</b> 2271–2296, doi:10.1175/MWR-D-17-0333.1</a> 
    <!--189-->
    <li>Coffer, B.E. and M.D. Parker, 2018: <a href="https://doi.org/10.1175/MWR-D-18-0050.1">Is there a “tipping point” between simulated nontornadic and tornadic supercells in VORTEX2 environments?</a>. <i>Mon. Wea. Rev.,</i> <b>146,</b> 2667–2693, doi:10.1175/MWR-D-18-0050.1. 
    <!--188-->
    <li>Houston, A.L. and J.M. Keeler, 201: <a href="https://doi.org/10.1175/JTECH-D-18-0019.1">The impact of sensor response and airspeed on the representation of the convective boundary layer and airmass boundaries by small unmanned aircraft systems</a>.  <i>J. Atmos. Oceanic Technol.,</i>, <b>35,</b> 1687-1699, doi:10.1175/JTECH-D-18-0019.1.
    <!--187-->
    <li>Smith, R. K., Montgomery, M. T. and Kilroy, G., 2018:  <a href="https://doi.org/10.1002/qj.3332">The generation of kinetic energy in tropical cyclones revisited</a>.  <i>Q. J. Roy. Meteor. Soc.</i>, <b>144,</b> 2481–2490, doi:10.1002/qj.3332.  
    <!--186-->
    <li>Yao, D., H. Xue, J. Yin, J. Sun, X. Liang, and J. Guo, 2018: <a href="https://doi.org/10.1007/s13351-018-7083-0">Investigation into the formation, structure, and evolution of an EF4 tornado in east China using a high-resolution numerical simulation</a>.  <i>J. Meteor. Res.,</i> <b>32,</b> 157-171, doi:10.1007/s13351-018-7083-0. 
    <!--185-->
    <li>Oreskovic, C., E. Savory, J. Porto, and L. G. Orf, 2018: <a href="https://doi.org/10.12989/was.2018.26.3.147">Evolution and scaling of a simulated downburst-producing thunderstorm outflow</a>.  <i>Wind and Structures,</i> <b>26,</b> 147-161, doi:10.12989/was.2018.26.3.147.
    <!--184-->
    <li>Nielsen, E. R., and R. S. Schumacher, 2018: <a href="https://doi.org/10.1175/JAS-D-17-0385.1">Dynamical insights into extreme short-term precipitation associated with supercells and mesovortices</a>.  <i>J. Atmos. Sci.,</i>, <b>75,</b> 2983-3009, doi:10.1175/JAS-D-17-0385.1.
    <!--183-->
    <li>Stratman, D. R., C. Potvin, and L. Wicker, 2018: <a href="https://doi.org/10.1175/MWR-D-17-0357.1">Correcting storm displacement errors in ensembles using the Feature Alignment Technique (FAT)</a>.  <i>Mon. Wea. Rev.,</i> <b>146,</b> 2125-2145, doi:10.1175/MWR-D-17-0357.1.
    <!--182-->
    <li>Kieu, C., and D.-L. Zhang, 2018:  <a href="https://doi.org/10.1029/2018GL078070">The control of environmental stratification on the hurricane maximum potential intensity</a>.  <i>Geophys. Res. Lett.,</i> <b>45,</b> 6272-6280, doi:10.1029/2018GL078070.
    <!--181-->
    <li>Lee, M., and T. Frisius, 2018:  <a href="https://doi.org/10.1080/16000870.2018.1433433">On the role of convective available potential energy (CAPE) in tropical cyclone intensification</a>.  <i>Tellus A,</i> <b>70:1,</b> 1433433, doi:10.1080/16000870.2018.1433433.
    <!--180-->
    <li>Xu, J., and Y. Wang, 2018: <a href="https://doi.org/10.2151/jmsj.2018-014">Effect of the initial vortex structure on intensification of a numerically simulated tropical cyclone</a>. <i>J. Meteor. Soc. Japan,</i> <b>96,</b> 111-126. 
    <!--179-->
    <li>Peng, K., R. Rotunno, and G. H. Bryan, 2018: <a href="http://dx.doi.org/10.1175/JAS-D-17-0382.1">Evaluation of a time-dependent model for the intensification of tropical cyclones</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 2125-2138, doi:10.1175/JAS-D-17-0382.1.
    <!--178-->
    <li>Morales, A., H. Morrison, and D. J. Posselt, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0389.1">Orographic precipitation response to microphysical parameter perturbations for idealized moist nearly neutral flow</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 1933-1953, doi:10.1175/JAS-D-17-0389.1.
    <!--177-->
    <li>Shi, X., F. K. Chow, R. L. Street, and G. H. Bryan, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0392.1">An evaluation of LES turbulence models for scalar mixing in the stratocumulus-capped boundary layer</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 1499-1507, doi:10.1175/JAS-D-17-0392.1.
    <!--176-->
    <li>Vande Guchte, A., and J. M. L. Dahl, 2018:   <a href="http://dx.doi.org/10.1175/MWR-D-17-0190.1">Sensitivities of parcel trajectories beneath the lowest scalar model level of a Lorenz vertical grid</a>. <i>Mon. Wea. Rev.,</i> <b>146,</b> 1427-1435, doi:10.1175/MWR-D-17-0190.1.
    <!--175-->
    <li>Limpert, G. L., and A. L. Houston, 2018:  <a href="http://dx.doi.org/10.1175/MWR-D-17-0029.1">Ensemble sensitivity analysis for targeted observations of supercell thunderstorms</a>.  <i>Mon. Wea. Rev.,</i> <b>146,</b> 1705-1721, doi:10.1175/MWR-D-17-0029.1.
    <!--174-->
    <li>Rotunno, R., and G. H. Bryan, 2018: <a href="http://dx.doi.org/10.1175/JAS-D-17-0306.1">Numerical simulations of two-layer flow past topography.  Part I. The lee-side hydraulic jump</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 1231-1241, doi:10.1175/JAS-D-17-0306.1. 
    <!--173-->
    <li>Lombardo, K., and T. Kading, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0248.1">The behavior of squall lines in horizontally heterogeneous coastal environments</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 1243-1269, doi:10.1175/JAS-D-17-0248.1. 
    <!--172-->
    <li>Guarriello, F., C. J. Nowotarski, and C. C. Epifanio, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0174.1">Effects of the low-level wind profile on outflow position and near-surface vertical vorticity in simulated lupercell thunderstorms</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 731-753, doi:10.1175/JAS-D-17-0174.1. 
    <!--171-->
    <li>Morrison, H., and J. M. Peters, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0295.1">Theoretical expressions for the ascent rate of moist deep convective thermals</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 1699-1719, doi:10.1175/JAS-D-17-0295.1. 
    <!--170-->
    <li>Frank, L. R., V. L. Galinsky, L. Orf, and J. Wurman, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0117.1">Dynamic multi-scale modes of severe storm structure detected in mobile Doppler radar data by entropy field decomposition</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 709-730, doi:10.1175/JAS-D-17-0117.1.
    <!--169-->
    <li>Shi, X., H. L. Hagen, F. K. Chow, G. H. Bryan, and R. L. Street, 2018:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0162.1">Large-eddy simulation of the stratocumulus-capped boundary layer with explicit filtering and reconstruction turbulence modeling</a>.  <i>J. Atmos. Sci.,</i> <b>75,</b> 611-637, doi:10.1175/JAS-D-17-0162.1. 
    <!--164x-->
    <li>Kilroy, G., R. K. Smith, and M. T. Montgomery, 2018:  <a href="http://dx.doi.org/10.1002/qj.3187">The role of heating and cooling associated with ice processes on tropical cyclogenesis and intensification</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>144,</b> 99-114, doi:10.1002/qj.3187. 
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>
  
  <details>
  <summary>2017</summary>
  <ul>
    <!--168-->
    <li>Coffer, B. E., M. D. Parker, J. M. L. Dahl, L. J. Wicker, and A. J. Clark, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-17-0152.1">Volatility of tornadogenesis:  An ensemble of simulated nontornadic and tornadic supercells in VORTEX2 environments</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 4605-4625, doi:10.1175/MWR-D-17-0152.1.
    <!--167-->
    <li>Stern, D. P., J. L. Vigh, D. S. Nolan, and F. Zhang, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0120.1">Reply to "Comments on 'Revisiting the relationship between eyewall contraction and intensification'"</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 4275-4286, doi:10.1175/JAS-D-17-0120.1.
    <!--166-->
    <li>Trapp, R. J., G. R. Marion, and S. W. Nesbitt, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0331.1">The regulation of tornado intensity by updraft width</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 4199-4211, doi:10.1175/JAS-D-16-0331.1. 
    <!--165-->
    <li>Parker, M. D., 2017:  <a href="http://dx.doi.org/10.1175/WAF-D-17-0064.1">How much does ‘‘backing aloft’’ actually impact a supercell?</a>  <i>Wea. Forecasting,</i> <b>32,</b> 1937-1957, doi:10.1175/WAF-D-17-0064.1.
    <!--163-->
    <li>Kilroy, G., and R. K. Smith, 2017:  <a href="http://dx.doi.org/10.1002/qj.3134">The effects of initial vortex size on tropical cyclogenesis and intensification</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>143,</b> 2832-2845, doi:10.1002/qj.3134. 
    <!--162-->
    <li>Kilroy, G., M. T. Montgomery, and R. K. Smith, 2017:  <a href="http://dx.doi.org/10.1002/qj.3104">The role of boundary-layer friction on tropical cyclogenesis and subsequent intensification</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>143,</b> 2524-2536, doi:10.1002/qj.2934. 
    <!--161-->
    <li>Diao, M., G. H. Bryan, H. Morrison, and J. B. Jensen, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0356.1">Ice nucleation parameterization and relative humidity distribution in idealized squall line simulations</a>.   <i>J. Atmos. Sci.,</i> <b>74,</b> 2761-2787, doi:10.1175/JAS-D-16-0356.1.
    <!--160-->
    <li>Giles Harrison, R., G. Pretor-Pinney, G. J. Marlton, G. D. Anderson, D. J. Kirshbaum, and R. J. Hogan, 2017:  <a href = "http://dx.doi.org/10.1002/wea.2996">Asperitas – a newly identified cloud supplementary feature.</a>  <i>Weather,</i> <b>72,</b> 132-141, doi:10.1002/wea.2996.
    <!--159-->
    <li>Navarro, E. L., G. J. Hakim, and H. E. Willoughby, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0279.1">Balanced response of an axisymmetric tropical cyclone to periodic diurnal heating</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 3325-3337, doi:10.1175/JAS-D-16-0279.1.
    <!--158-->
    <li>Dahl, J. M. L., 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-17-0091.1">Tilting of horizontal shear vorticity and the development of updraft rotation in supercell thunderstorms</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 2997-3020, doi:10.1175/JAS-D-17-0091.1.
    <!--157-->
    <li>Madaus, L. E., and G. J. Hakim, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0395.1">Constraining ensemble forecasts of discrete convective initiation with surface observations</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 2597-2610, doi:10.1175/MWR-D-16-0395.1.
    <!--156-->
    <li>MacIntosh, C. W. and M. D. Parker, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0329.1">The 6 May 2010 elevated supercell during VORTEX2</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 2635-2657, doi:10.1175/MWR-D-16-0329.1.
    <!--155-->
    <li>Schecter, D. A., 2017:  <a href="http://dx.doi.org/10.1002/2016MS000777">A computational study on the nature of meso-β scale vortex coalescence in a tropical atmosphere</a>.  <i>J. Adv. Modeling Earth Systems,</i> <b>9,</b> doi:10.1002/2016MS000777.
    <!--154-->
    <li>Worsnop, R. P., J. K. Lunquist, G. H. Bryan, R. Damiani, and W. Musial, 2017:  <a href="http://dx.doi.org/10.1002/2017GL073537">Gusts and shear within hurricane eyewalls can exceed offshore wind-turbine design standards</a>.  <i>Geophys. Res. Lett.,</i> <b>44,</b> 6413-6420, doi:10.1002/2017GL073537. 
    <!--153-->
    <li>Worsnop, R. P., G. H. Bryan, J. K. Lunquist, and J. A. Zhang, 2017:  <a href="http://dx.doi.org/10.1007/s10546-017-0266-x">Using large-eddy simulations to define spectral and coherence characteristics of the hurricane boundary layer for wind-energy applications</a>.  <i>Bound.-Layer Meteor.,</i> <b>165,</b> 55-86, doi:10.1007/s10546-017-0266-x.
    <!--152-->
    <li>Lehner, M., C. D. Whiteman, and M. Dorninger, 2017:  <a href="http://dx.doi.org/10.1007/s10546-017-0232-7">Inversion build-up and cold-air outflow in a small alpine sinkhole</a>.  <i>Bound.-Layer Meteor.,</i> <b>163,</b> 497-522, doi:10.1007/s10546-017-0232-7.
    <!--151-->
    <li>Bu, Y. P., R. G. Fovell, and K. L. Corbosiero, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0231.1">The influences of boundary layer mixing and cloud-radiative forcing on tropical cyclone size</a>.   <i>J. Atmos. Sci.,</i> <b>74,</b> 1273-1292, doi:10.1175/JAS-D-16-0258.1.
    <!--150-->
    <li>Warren, R. A., H. Richter, H. A. Ramsay, S. T. Siems, and M. J. Manton, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0412.1">Impact of variations in upper-level shear on simulated supercells</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 2659-2681, doi:10.1175/MWR-D-16-0412.1.
    <!--149-->
    <li>Rotunno, R., P. M. Markowski, and G. H. Bryan, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0288.1">`Near-ground' vertical vorticity in supercell thunderstorm models</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 1757-1766, doi:10.1175/JAS-D-16-0288.1.
    <!--148-->
    <li>Dahl, N. A., D. S. Nolan, G. H. Bryan, and R. Rotunno, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0346.1">Using high-resolution simulations to quantify underestimates of tornado intensity from in situ observations</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 1963-1982, doi:10.1175/MWR-D-16-0346.1.
    <!--147-->
    <li>Nolan, D. S., N. A. Dahl, G. H. Bryan, and R. Rotunno, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0258.1">Tornado vortex structure, intensity, and surface wind gusts in large-eddy simulations with fully developed turbulence</a>.   <i>J. Atmos. Sci.,</i> <b>74,</b> 1573-1597, doi:10.1175/JAS-D-16-0258.1.
    <!--146-->
    <li>Bryan, G. H., N. A. Dahl, D. S. Nolan, and R. Rotunno, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0339.1">An eddy injection method for large-eddy simulations of tornado-like vortices</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 1937-1961, doi:10.1175/MWR-D-16-0339.1.
    <!--145-->
    <li>Markowski, P. M., and Y. P. Richardson, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0372.1">Large sensitivity of near-surface vertical vorticity development to heat sink location in idealized simulations of supercell-like storms</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 1095-1104, doi:10.1175/JAS-D-16-0372.1.
    <!--144-->
    <li>Trapp, R. J., and J. M. Woznicki, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0266.1">Convectively induced stabilizations and subsequent recovery with supercell thunderstorms during the Mesoscale Predictability Experiment (MPEX)</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 1739-1754, doi:10.1175/MWR-D-16-0266.1.
    <!--143-->
    <li>Rousseau-Rizzi, R., D. J. Kirshbaum, and M. K. Yau, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0221.1">Initiation of deep convection over an idealized mesoscale convergence line</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 835-853, doi:10.1175/JAS-D-16-0221.1.
    <!--142-->
    <li>Sever, G., and Y.-L. Lin, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0077.1">Dynamical and physical processes associated with orographic precipitation in a conditionally unstable uniform flow: variation in basic wind speed</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 449-466, doi:10.1175/JAS-D-16-0077.1.
    <!--123-->
    <li>Orf, L., R. Wilhelmson, B. Lee, C. Finley, and A. Houston, 2017:  <a href="http://dx.doi.org/10.1175/BAMS-D-15-00073.1">Evolution of a long-track violent tornado within a simulated supercell</a>.  <i>Bull. Amer. Meteor. Soc.,</i> <b>98,</b>45-68, doi:10.1175/BAMS-D-15-00073.1.
    <!--133-->
    <li>Ditchek, S. D., J. Molinari, and D. Vollaro, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0117.1">Tropical cyclone outflow layer structure and balanced response to eddy forcings</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 133-149, doi:10.1175/JAS-D-16-0117.1. 
    <!--134-->
    <li>Alfaro, D. A., 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0168.1">Low-tropospheric shear in the structure of squall lines: Impacts on latent heating under layer-lifting ascent</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 229-248, doi:10.1175/JAS-D-16-0168.1.
    <!--135-->
    <li>Coffer, B. E., and M. D. Parker, 2017:  <a href="http://dx.doi.org/10.1175/MWR-D-16-0226.1">Simulated supercells in nontornadic and tornadic VORTEX2 environments</a>.  <i>Mon. Wea. Rev.,</i> <b>145,</b> 149-180, doi:10.1175/MWR-D-16-0226.1.
    <!--137-->
    <li>Sachsperger, J.,  S. Serafin ,  V. Grubisic,  I. Stiperski, and A. Paci, 2017:  <a href="http://dx.doi.org/10.1002/qj.2915">The amplitude of lee waves on the boundary-layer inversion</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>143,</b> 27-36, doi:10.1002/qj.2915.
    <!--138-->
    <li>Kilroy, G., R. K. Smith, and M. T. Montgomery, 2017:  <a href="http://dx.doi.org/10.1002/qj.2934">A unified view of tropical cyclogenesis and intensification</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>143,</b> 450-462, doi:10.1002/qj.2934.
    <!--139-->
    <li>Betten, D. P., M. I. Biggerstaff, and L. J. Wicker, 2017:  <a href="http://dx.doi.org/10.1175/JTECH-D-16-0043.1">A trajectory mapping technique for the visualization and analysis of three-dimensional flow in supercell storms</a>.  <i>J. Atmos. Oceanic Technol.,</i> <b>34,</b> 33-49, doi:10.1175/JTECH-D-16-0043.1.
    <!--141-->
    <li>Bryan, G. H., R. P. Worsnop, J. K. Lundquist, and J. A. Zhang, 2017:  <a href="http://dx.doi.org/10.1007/s10546-016-0207-0">A simple method for simulating wind profiles in the boundary layer of tropical cyclones</a>.  <i>Bound.-Layer Meteor.,</i> <b>3,</b> 475-502, doi:10.1007/s10546-016-0207-0. 
    <!--140-->
    <li>Dennis, E. J., and M. R. Kumjian, 2017:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0066.1">The impact of vertical wind shear on hail growth in simulated supercells</a>.  <i>J. Atmos. Sci.,</i> <b>74,</b> 641-663, doi:dx.doi.org/10.1175/JAS-D-16-0066.1.
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>

  <details>
  <summary>2016</summary>
  <ul>
    <!--136-->
    <li>Singh, M. S., and P. A. O'Gorman, 2016:  <a href="http://dx.doi.org/10.1002/2016MS000673">Scaling of the entropy budget with surface temperature in radiative-convective equilibrium</a>.  <i>J. Adv. Model. Earth Syst.,</i> <b>8,</b> 1132-1150, doi:10.1002/2016MS000673. 
</a>.
    <!--132-->
    <li>Ruppert, J. H., Jr, 2016:  <a href="https://doi.org/10.1002/2016MS000713">Diurnal timescale feedbacks in the tropical cumulus regime</a>.  <i>J. Adv. Model. Earth Syst.,</i> <b>8,</b> 1483-1500, doi:10.1002/2016MS000713.
    <!--131-->
    <li>Fovell, R. G., Y. P. Bu, K. L. Corbosiero, W.-W. Tung, Y. Cao, H.-C. Kuo, L.-H. Hsu, and H. Su, 2016:  <a href="https://doi.org/10.1175/AMSMONOGRAPHS-D-15-0006.1">Influence of cloud microphysics and radiation on tropical cyclone structure and motion</a>.  <i>Multiscale Convection-Coupled Systems in the Tropics: A Tribute to Dr. Michio Yanai, Meteor. Monogr.,</i> No. 56, 11.1-11.27.
    <!--130-->
    <li>Markowski, M., 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-16-0150.1">An idealized numerical simulation investigation of the effects of surface drag on the development of near-surface vertical vorticity in supercell thunderstorms</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 4349-4385, doi:10.1175/JAS-D-16-0150.1.
    <!--129-->
    <li>Kilroy, G., and R. K. Smith, 2016:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.2895/abstract">A numerical study of deep convection in tropical cyclones</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>142,</b> 3138-3151, doi:10.1002/qj.2895.
    <!--128-->
    <li>Peters, J. M., 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-16-0016.1">The impact of effective buoyancy and dynamic pressure forcing on vertical velocities within 2 dimensional updrafts</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 4531-4551, doi:10.1175/JAS-D-16-0016.1.
    <!--127-->
    <li>Rotunno, R., and M. Lehner, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-16-0132.1">Two-layer stratified flow past a valley</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 4065-4076, doi:10.1175/JAS-D-16-0132.1.
    <!--126-->
    <li>Lehner, M., R. Rotunno, and C. D. Whiteman, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-16-0114.1">Flow regimes over a basin induced by upstream katabatic flows -- An idealized modeling study</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 3821-3842, doi:10.1175/JAS-D-16-0114.1.
    <!--125-->
    <li>Navarro, E. L., and G. J. Hakim, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0349.1">Idealized numerical modeling of the diurnal cycle of tropical cyclones</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 4189-4201, doi:10.1175/JAS-D-15-0349.1.
    <!--124-->
    <li>Ruppert Jr., J. H. and R. H. Johnson, 2016:  <a href="http://dx.doi.org/10.1002/2015MS000610">On the cumulus diurnal cycle over the tropical warm pool</a>.  <i>J. of Advances in Modeling Earth Systems,</i> <b>8,</b> 669-690, doi:10.1002/2015MS000610.
    <!--122-->
    <li>Rotunno, R., G. H. Bryan, D. S. Nolan, and N. A. Dahl, 2016:  <a href="http://dx.doi.org/10.1175/JAS-D-16-0038.1">Axisymmetric tornado simulations at high Reynolds number</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 3843-3854, doi:10.1175/JAS-D-16-0038.1.
    <!--121-->
    <li>Crnivec, N., R. K. Smith, and G. Kilroy, 2016:  <a href="http://dx.doi.org/10.1002/qj.2752">Dependence of tropical cyclone intensification rate on sea-surface temperature</a>.  <i>Q. J. Roy. Meteor. Soc.,</i> <b>142,</b> 1618-1627, doi:10.1002/qj.2752.
    <!--120-->
    <li>Dorier, M., O. Yildizc, S. Ibrahimc, A.-C. Orgeried, and G. Antoniuc.  <a href="http://dx.doi.org/10.1016/j.future.2016.03.002">On the energy footprint of I/O management in Exascale HPC systems</a>.  <i>Future Generation Computer Systems,</i> <b>62,</b> 17-28.
    <!--119-->
    <li>Sachsperger, J., S. Serafin, and V. Grubisic, 2016:  <a href="http://dx.doi.org/10.1002/qj.2746">Dynamics of rotor formation in uniformly stratified two-dimensional flow over a mountain</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> 142,</b> 1201-1212, doi:10.1002/qj.2746.
    <!--118-->
    <li>Madaus, L. E., and Gregory J. Hakim, 2016: <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0332.1">Observable Surface Anomalies Preceding Simulated Isolated Convective Initiation.</a>  <i>Mon. Wea. Rev.,</i> <b>144,</b> 2265-2284, doi: 10.1175/MWR-D-15-0332.1.
    <!--117-->
    <li>Richter, D. H., R. Bohac, and D. P. Stern, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0331.1">An assessment of the flux profile method for determining air-sea momentum and enthalpy fluxes from dropsonde data in tropical cyclones.</a>  <i>J. Atmos. Sci.,</i> <b>73,</b> 2665-2682, doi: 10.1175/JAS-D-15-0331.1.
    <!--116-->
    <li>Morrison, H., 2016:  <a href="http://journals.ametsoc.org/doi/full/10.1175/JAS-D-15-0040.1">Impacts of Updraft Size and Dimensionality on the Perturbation Pressure and Vertical Velocity in Cumulus Convection. Part I: Simple, Generalized Analytic Solutions</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 1441–1454, doi: 10.1175/JAS-D-15-0040.1.
    <!--115-->
    <li>Morrison, H., 2016:  <a href="http://journals.ametsoc.org/doi/full/10.1175/JAS-D-15-0041.1">Impacts of Updraft Size and Dimensionality on the Perturbation Pressure and Vertical Velocity in Cumulus Convection. Part II: Comparison of Theoretical and Numerical Solutions and Fully Dynamical Simulations</a>.  <i>J. Atmos. Sci.,</i> <b>73,</b> 1455–1480, doi: 10.1175/JAS-D-15-0041.1.
    <!--114-->
    <li>Markowski, P. M., and G. H. Bryan, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0439.1">LES of laminar flow in the PBL: A potential problem for convective storm simulations</a>.  <i>Mon. Wea. Rev.,</i> <b>144,</b> 1841-1850, doi:10.1175/MWR-D-15-0439.1.
    <!--113-->
    <li>Hastings, R., and Y. Richardson, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0193.1">Long-Term Morphological Changes in Simulated Supercells Following Mergers with Nascent Supercells in Directionally Varying Shear</a>.  <i>Mon. Wea. Rev.,</i> <b>144,</b> 471-499.
    <!--112-->
    <li>Nowotarski, C. J., and P. M. Markowski, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0247.1">Modifications to the Near-Storm Environment Induced by Simulated Supercell Thunderstorms</a>.  <i>Mon. Wea. Rev.,</i> <b>144,</b> 273-293.
    <!--111-->
    <li>Kepert, J. D., J. Schwendike, and H. Ramsey, 2016:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0216.1">Why is the Tropical Cyclone Boundary Layer Not “Well Mixed”?</a>  <i>J. Atmos. Sci.,</i> <b>73,</b> 957-973.
    <!--108-->
    <li>Orf, L., R. Wilhelmson, and L. Wicker, 2016:  <a href="http://dx.doi.org/10.1016/j.parco.2015.10.014">Visualization of a simulated long-track EF5 tornado embedded within a supercell thunderstorm</a>.  <i>Parallel Computing,</i> <b>55,</b> 28-34, doi:10.1016/j.parco.2015.10.014. 
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>

  <details>
  <summary>2015</summary>
  <ul>
    <!--110-->
    <li>Schecter, D. A., 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0028.1">Development and Nondevelopment of Binary Mesoscale Vortices into Tropical Cyclones in Idealized Numerical Experiments.</a>  <i>J. Atmos. Sci.</i>, <b>73,</b> 1223–1254, doi: 10.1175/JAS-D-15-0028.1.
    <!--109-->
    <li>Lehner, M., and Coauthors, 2015: <a href="http://journals.ametsoc.org/doi/full/10.1175/BAMS-D-14-00238.1">The METCRAX II Field Experiment: A Study of Downslope Windstorm-Type Flows in Arizona’s Meteor Crater></a>. <i>Bull. Amer. Meteor. Soc.,</i> <b>97,</b> 217–235, doi: 10.1175/BAMS-D-14-00238.1.
    <!--107-->
    <li>Dahl, J. M. L., 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0115.1">Near-ground rotation in simulated supercells:  On the robustness of the baroclinic mechanism</a>.  <i>Mon. Wea. Rev.,</i> <b>143,</b> 4929-4942, doi:10.1175/MWR-D-15-0115.1.
    <!--106-->
    <li>Tushaus, S. A. Tushaus, D. J. Posselt, M. M. Miglietta, R. Rotunno, and L. Delle Monache, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0036.1">Bayesian exploration of multivariate orographic precipitation sensitivity for moist stable and neutral flows</a>.  <i>Mon. Wea. Rev.,</i> <b>143,</b> 4459-4475, doi:10.1175/MWR-D-15-0036.1.
    <!--105-->
    <li>Davenport, C. E., and M. D. Parker, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-15-0072.1">Impact of environmental heterogeneity on the dynamics of a dissipating supercell thunderstorm</a>.  <i>Mon. Wea. Rev.,</i> <b>143,</b> 4244-4277, doi:10.1175/MWR-D-15-0072.1.
    <!--104-->
    <li>Shapiro, A., S. Rahimi, C. K. Potvin, and L. Orf, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0095.1">On the use of advection correction in trajectory calculations</a>.  <i>J. Atmos. Sci.,</i> <b>72,</b> 4261-4280, doi:10.1175/JAS-D-15-0095.1.
    <!--103-->
    <li>Lane, T. P., and M. W. Moncrieff, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0073.1">Long-lived mesoscale systems in a low convective inhibition environment. Part I: Upshear propagation</a>.  <i>J. Atmos.  Sci.,</i> <b>72,</b> 4297-4318, doi:10.1175/JAS-D-15-0073.1.
    <!--102-->
    <li>Moncrieff, M. W., and T. P. Lane, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0074.1">Long-lived mesoscale systems in a low convective inhibition environment. Part II: Downshear propagation</a>.  <i>J. Atmos. Sci.,</i> <b>72,</b> 4319-4336, doi:10.1175/JAS-D-15-0074.1.
    <!--101-->
    <li>Singh, M. S., and P. A. O'Gorman, 2015:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.2567/abstract">Increases in moist-convective updraught velocities with warming in radiative-convective equilibrium</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>141,</b> 2828-2838, doi:10.1002/qj.2567.
    <!--100-->
    <li>Davis, C. A., 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-15-0027.1">The formation of moist vortices and tropical cyclones in idealized simulations</a>.  <i>J. Atmos. Sci.,</i> <b>72,</b> 3499-3516, doi:10.1175/JAS-D-15-0027.1.
    <!--99-->
    <li>Cherukuru, N. W., R. Calhoun, M. Lehner, S. W. Hoch, and C. D. Whiteman, 2015:  <a href="http://remotesensing.spiedigitallibrary.org/article.aspx?articleid=2109805">Instrument configuration for dual-Doppler lidar coplanar scans: METCRAX II</a>.  <i>J. Appl. Remote Sens.,</i> <b>9(1),</b> 096090, doi:10.1117/1.JRS.9.096090.
    <!--98-->
    <li>Schumacher, R. S., 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-14-0389.1">Sensitivity of precipitation accumulation in elevated convective systems to small changes in low-level moisture</a>.  <i>J. Atmos. Sci.,</i> <b>72,</b> 2507-2524, doi:10.1175/JAS-D-14-0389.1.
    <!--97-->
    <li>Parker, M. D., and J. M. L. Dahl, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-14-00310.1">Production of near-surface vertical vorticity by idealized downdrafts</a>.  <i>Mon. Wea. Rev.,</i> <b>143,</b> 2795-2816, doi:10.1175/MWR-D-14-00310.1.
    <!--96-->
    <li>Montgomery, M. T., J. Persing, and R. K. Smith, 2015:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/2014MS000362/abstract">Putting to rest WISHE-ful misconceptions for tropical cyclone intensification</a>.  <i>Journal of Advances in Modeling Earth Systems,</i> <b>7,</b> 92-109, doi:10.1002/2014MS000362.
    <!--95-->
    <li>Coffer, B. E., and M. D. Parker, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-14-00328.1">Impacts of increasing low-level shear on supercells during the early evening transition</a>.  <i>Mon. Wea. Rev.,</i> <b>143,</b> 1945-1969, doi:10.1175/MWR-D-14-00328.1.
    <!--94-->
    <li>Ryglicki, D. A., 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-14-0180.1">An analysis of a barotropically-unstable barotropic, high-Rossby number vortex in shear</a>.  <i>J. Atmos. Sci.,</i> <b>72,</b> 2152-2177, doi:10.1175/JAS-D-14-0180.1.
    <!--93-->
    <li>Nowotarski, C. J., P. M. Markowski, Y. P. Richardson, and G. H. Bryan, 2015:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-14-00151.1">Supercell low-level mesocyclones in simulations with a sheared convective boundary layer</a>.  <i>Mon. Wea. Rev.,</i> <b>143,</b> 272-297, doi:10.1175/MWR-D-14-00151.1. 
    <!--89-->
    <li>Kilroy, G., and R. K. Smith, 2015:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.2383/abstract">Tropical cyclone convection: the effects of a vortex boundary-layer wind profile on deep convection</a>.  <i>Q. J. Roy. Meteor. Soc.,</i> <b>141,</b> 714-726, doi:10.1002/qj.2383.
  </ul>
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  <!--ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-->
  </details>

  <details>
  <summary>2014</summary>
  <ul>
    <!--92-->
    <li>Lane, T. P., and R. D. Sharman, 2014:  <a href="http://dx.doi.org/10.1002/2014GL059299">Intensity of thunderstorm-generated turbulence revealed by large-eddy simulation</a>.  <i>Geophys. Res. Letters,</i> <b>41,</b> 2221-2227, doi:10.1002/2014GL059299.
    <!--91-->
    <li>Orf, L. G., C. Oreskovic, E. Savory, and E. Kantor, 2014:  <a href="http://www.sciencedirect.com/science/article/pii/S016761051400141X">Circumferential analysis of a simulated three-dimensional downburst-producing thunderstorm outflow</a>.  <i>J. Wind Engineering and Industrial Aerodynamics,</i> <b>135,</b> 182-190, doi:10.1016/j.jweia.2014.07.004.
    <!--90-->
    <li>French, A. J., and M. D. Parker, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-13-00356.1">Numerical simulations of bow echo formation following a squall line–supercell merger</a>.  <i>Mon. Wea. Rev.,</i> <b>142,</b> 4791-4822, doi:10.1175/MWR-D-13-00356.1.
    <!--88-->
    <li>Smith, R. K., M. T. Montgomery, and J. Persing, 2014:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.2329/abstract">On steady-state tropical cyclones</a>.  <i>Q. J. Roy. Meteor. Soc.,</i> <b>140,</b> 2638-2649, doi:10.1002/qj.2329.
    <!--87-->
    <li>Naylor, J., and D. A. Schecter, 2014:  <a href="http://onlinelibrary.wiley.com/enhanced/doi/10.1002/2014MS000366">Evaluation of the impact of moist convection on the development of asymmetric inner core instabilities in simulated tropical cyclones</a>.  <i>Journal of Advances in Modeling Earth Systems,</i> doi:10.1002/2014MS000366.
    <!--86-->
    <li>Markowski, P., Y. Richardson, and G. H. Bryan, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-14-00162.1">The origins of vortex sheets in a simulated supercell thunderstorm</a>.  <i>Mon. Wea. Rev.,</i> <b>142,</b> 3944-3954, doi:10.1175/MWR-D-14-00162.1.
    <!--85-->
    <li>Nowotarski, C. J., P. M. Markowski, Y. P. Richardson, and G. H. Bryan, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-13-00349.1">Properties of a simulated convective boundary layer in an idealized supercell thunderstorm environment</a>.  <i>Mon. Wea. Rev.,</i> <b>142,</b> 3955-3976, doi:10.1175/MWR-D-13-00349.1.
    <!--84-->
    <li>Navarro, E. L., and G. J. Hakim, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-13-00099.1">Storm-centered ensemble data assimilation for tropical cyclones</a>.  <i>Mon. Wea. Rev.,</i> <b>142,</b> 2309-2320.
    <!--83-->
    <li>Singh, M. S., and P. A. O'Gorman, 2014:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/2014GL061222/abstract">Influence of microphysics on the scaling of precipitation extremes with temperature</a>.  <i>Geophys. Res. Lett.,</i> <b>41,</b> 6037-6044, doi:10.1002/2014GL061222.
    <!--82-->
    <li>Kilroy, G., R. K. Smith, and U. Wissmeier, 2014:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.2261/abstract">Tropical convection: The effects of ambient vertical and horizontal vorticity</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>140,</b> 1756-1770, doi:10.1002/qj.2261. 
    <!--81-->
    <li>Wang, W., 2014:  <a href="http://link.springer.com/article/10.1007%2Fs10546-013-9899-6">Analytically modelling mean wind and stress profiles in canopies</a>.  <i>Bound.-Layer Meteor.,</i> <b>151,</b> 239-256. 
    <!--80-->
    <li>Dahl, J. M. L., M. D. Parker, and L. J. Wicker, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0123.1">Imported and storm-generated near-ground vertical vorticity in a simulated supercell</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 3027-3051, doi:10.1175/JAS-D-13-0123.1. 
    <!--79-->
    <li>Soderholm, B., B. Ronalds, and D. J. Kirshbaum, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-13-00280.1">The evolution of convective storms initiated by an isolated mountain ridge</a>.  <i>Mon. Wea. Rev.,</i> <b>142,</b> 1430-1451.
    <!--78-->
    <li>Miglietta, M. M., and R. Rotunno, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0297.1">Numerical simulations of sheared conditionally unstable flows over a mountain ridge</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 1747-1762. 
    <!--77-->
    <li>Chavas, D. R., and K. Emanuel, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0155.1">Equilibrium tropical cyclone size in an idealized state of axisymmetric radiative-convective equilibrium</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 1663-1680, doi:10.1175/JAS-D-13-0155.1. 
    <!--76-->
    <li>Kirshbaum, D. J., and C.-C. Wang, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0287.1">Boundary-layer updrafts driven by airflow over heated terrain</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 1425-1442, doi:10.1175/JAS-D-13-0287.1.
    <!--75-->
    <li>Markowski, P. M., and Y. P. Richardson, 2014:  <a href="http://dx.doi.org/10.1175/JAS-D-13-0159.1">The influence of environmental low-level shear and cold pools on tornadogenesis: Insights from idealized simulations</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 243-275. 
    <!--74-->
    <li>Bryan, G. H., and R. Rotunno, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0157.1">Gravity currents in confined channels with environmental shear</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 1121-1142, doi:10.1175/JAS-D-13-0157.1.
    <!--73-->
    <li>Bryan, G. H., and R. Rotunno, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0156.1">The optimal state for gravity currents in shear.</a>  <i>J. Atmos. Sci.,</i> <b>71,</b> 448-468, doi:10.1175/JAS-D-13-0156.1. 
    <!--72-->
    <li>Naylor, J., and M. S. Gilmore, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0219.1">Vorticity evolution leading to tornadogenesis and tornadogenesis failure in simulated supercells.</a>  <i>J. Atmos. Sci.,</i> <b>71,</b> 1201-1217, doi:10.1175/JAS-D-13-0219.1.
    <!--70-->
    <li>Bu, Y. P., R. G. Fovell, and K. L. Corbosiero, 2014:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-13-0265.1">Influence of cloud-radiative forcing on tropical cyclone structure</a>.  <i>J. Atmos. Sci.,</i> <b>71,</b> 1644-1662, doi:10.1175/JAS-D-13-0265.1.
  </ul>
  </details>

  <details>
  <summary>2013</summary>
  <ul>
    <!--71-->
    <li>Persing, J., M. T. Montgomery, J. C. McWilliams, and R. K. Smith, 2013:  <a href="http://www.atmos-chem-phys.net/13/12299/2013/acp-13-12299-2013.html">Asymmetric and axisymmetric dynamics of tropical cyclones</a>.  <i>Atmos. Chem. Phys.,</i> <b>13,</b> 12299-12341.
    <!--69-->
    <li>Nicolae, B., and F. Cappello, 2013:  <a href="http://dx.doi.org/10.1016/j.jpdc.2013.01.013">BlobCR: Virtual disk based checkpoint-restart for HPC applications on IaaS clouds.</a>  <i>J. Parallel and Distributed Computing.,</i> <b>73,</b> 698-711, 
    <!--68-->
    <li>Kilroy, G., and R. K. Smith, 2013:  <a href="http://dx.doi.org/10.1002/qj.2022">A numerical study of rotating convection during tropical cyclogenesis</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>139,</b> 1255-1269.
    <!--67-->
    <li>Singh, M. S., and P. A. O'Gorman, 2013:  <a href="http://dx.doi.org/10.1002/grl.50796">Influence of entrainment on the thermal stratification in simulations of radiative-convective equilibrium.</a> <i>Geophys. Res. Let.,</i> <b>40,</b> doi: 10.1002/grl.50796.
    <!--66-->
    <li>Davies-Jones, R., and P. Markowski, 2013:  <a href="http://dx.doi.org/10.1175/JAS-D-12-0149.1">Lifting of ambient air by density currents in sheared environments</a>.  <i>J. Atmos. Sci.,</i> <b>70,</b> 1204-1215. 
    <!--65-->
    <li>Adams-Selin, R. D., and R. H. Johnson, 2013:  <a href="http://dx.doi.org/10.1175/MWR-D-12-00343.1">Examination of gravity waves associated with 13 March 2003 bow echo</a>.  <i>Mon. Wea. Rev.,</i> <b>141,</b> 3735-3756. 
    <!--64-->
    <li>Hakim, G. J., 2013:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-12-0188.1">The variability and predictability of axisymmetric hurricanes in statistical equilibrium</a>.  <i>J. Atmos. Sci.,</i> <b>70,</b> 993-1005. 
    <!--63-->
    <li>Brown, B. R., and G. J. Hakim, 2013:  <a href="http://dx.doi.org/10.1175/JAS-D-12-0112.1">Variability and predictability of a three-dimensional hurricane in statistical equilibrium</a>.  <i>J. Atmos. Sci.,</i> <b>70,</b> 1806-1820. 
    <!--62-->
    <li>Kirshbaum, D. J., 2013:  <a href="http://dx.doi.org/10.1175/JAS-D-12-0199.1">On thermally forced circulations over heated terrain</a>.  <i>J. Atmos. Sci.,</i> <b>70,</b> 1690-1709. 
    <!--61-->
    <li>Letkewicz, C. E., A. J. French, and M. D. Parker, 2013:  <a href="http://dx.doi.org/10.1175/MWR-D-12-00200.1">Base-state substitution:  An idealized modeling technique for approximating environmental variability</a>.  <i>Mon. Wea. Rev.,</i> <b>141,</b> 3062-3086.
    <!--60-->
    <li>Ramsay, H. A., 2013:  <a href="http://dx.doi.org/10.1175/JCLI-D-13-00195.1">The effects of imposed stratospheric cooling on the maximum intensity of tropical cyclones in axisymmetric radiative-convective equilibrium</a>.  <i>J. Climate,</i> <b>26,</b> 9977-9985. 
    <!--59-->
    <li>Bryan, G. H., 2013:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.2066/abstract">Comments on "Sensitivity of tropical-cyclone models to the surface drag coefficient"</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>139,</b> 1957-1960, doi:10.1002/qj.2066. 
  </ul>
  </details>

  <details>
  <summary>2012</summary>
  <ul>
    <!--58-->
    <li>Naylor, J., and M. S. Gilmore, 2012:  <a href="http://www.agu.org/pubs/crossref/2012/2012GL053041.shtml">Environmental factors influential to the duration and intensity of tornadoes in simulated supercells</a>.  <i>Geophys. Res. Let.,</i> <b>39,</b> L17802, doi:10.1029/2012GL053041
    <!--57-->
    <li>Naylor, J., and M. S. Gilmore, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-12-00163.1">Convective initiation in an idealized cloud model using an updraft nudging technique</a>.  <i>Mon. Wea. Rev.,</i> <b>140,</b> 3699-3705, doi:10.1175/MWR- D-12-00163.1. 
    <!--56-->
    <li>Orf, L., and E. Kantor, and E. Savory, 2012:  <a href="http://www.sciencedirect.com/science/article/pii/S016761051200044X">Simulation of a downburst-producing thunderstorm using a very high-resolution three-dimensional cloud model</a>.  <i>J. of Wind Engineering and Industrial Aerodynamics,</i> <b>104,</b> 547-557.
    <!--55-->
    <li>Kirshbaum, D. J., and A. L. M. Grant, 2012:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.1954/abstract">Invigoration of cumulus cloud fields by mesoscale ascent</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>138,</b> 2136-2150, doi:10.1002/qj.1954
    <!--54-->
    <li>Dahl, J. M. L., M. D. Parker, and L. J. Wicker, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-12-00131.1">Uncertainties in trajectory calculations within near-surface mesocyclones of simulated supercells</a>.  <i>Mon. Wea. Rev.,</i> <b>140,</b> 2959-2966, doi:10.1175/MWR-D-12-00131.1
    <!--53-->
    <li>Naylor, J., M. S. Gilmore, R. L. Thompson, R. Edwards, and R. B. Wilhelmson, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-11-00209.1">Comparison of objective supercell identification techniques using an idealized cloud model</a>.  <i>Mon. Wea. Rev.,</i> <b>140,</b> 2090-2102, doi:10.1175/MWR-D-11-00209.1
    <!--52-->
    <li>Rotunno, R., and G. H. Bryan, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-11-0204.1">Effects of parameterized diffusion on simulated hurricanes</a>.  <i>J. Atmos. Sci.,</i> <b>69,</b> 2284-2299, doi:10.1175/JAS-D-11-0204.1
    <!--51-->
    <li>Miglietta, M. M., and R. Rotunno, 2012:  <a href="http://dx.doi.org/10.1175/MWR-D-11-00253.1">Application of theory to simulations of observed cases of orographically forced convective rainfall</a>.  <i>Mon. Wea. Rev.,</i> <b>140,</b> 3039-3053, doi:10.1175/MWR-D-11-00253.1
    <!--50-->
    <li>Kang, S.-L., D. Lenschow, and P. Sullivan, 2012:  <a href="http://www.springerlink.com/content/f8362q1337107184/">Effects of mesoscale surface thermal heterogeneity on low-level horizontal wind speeds </a>.  <i>Bound. Layer Meteor.,</i> <b>143,</b> 409-432, doi:10.1007/s10546-011-9691-4.
    <!--49-->
    <li>Bryan, G. H., 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-11-00231.1">Effects of surface exchange coefficients and turbulence length scales on the intensity and structure of numerically simulated hurricanes</a>.  <i>Mon. Wea. Rev.,</i> in <b>140,</b> 1135-1143, doi:10.1175/MWR-D-11-00231.1.
    <!--48-->
    <li>Billings, J. M., and M. D. Parker, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/WAF-D-11-00056.1">Evolution and maintenance of the 22-23 June 2003 nocturnal convection during BAMEX</a>.  <i>Wea. Forecasting,</i> <b>27,</b> 279-300, doi:10.1175/WAF-D-11-00056.1
    <!--47-->
    <li>Parker, M. D., 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/JAS-D-11-058.1">Impacts of lapse rates upon low-level rotation in idealized storms</a>.  <i>J. Atmos. Sci.,</i> <b>69,</b> 538-559, doi:10.1175/JAS-D-11-058.1
    <!--46-->
    <li>Reasor, P. D., and M. D. Eastin, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-11-00080.1">Rapidly Intensifying Hurricane Guillermo (1997). Part II: Resilience in Shear</a>.  <i>Mon. Wea. Rev.,</i> <b>140,</b> 425-444, doi:10.1175/MWR-D-11-00080.1
    <!--45-->
    <li>Bryan, G. H., and H. Morrison, 2012:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-11-00046.1">Sensitivity of a simulated squall line to horizontal resolution and parameterization of microphysics</a>.  <i>Mon. Wea. Rev.,</i> <b>140,</b> 202-225, doi:10.1175/MWR-D-11-00046.1
  </ul>
  </details>
  
  <details>
  <summary>2011</summary>
  <ul>
    <!--342-->
    <li>Rotunno, R., Klemp, J., Bryan, G., and Muraki, D., 2011: <a href="https://doi.org/10.1017/jfm.2010.648">Models of non-Boussinesq lock-exchange flow</a>. <i>Journal of Fluid Mechanics,</i> </b>675,</b> 1-26, doi:10.1017/jfm.2010.648. 
    <!--44-->
    <li>Vermeire, B. C., L. G. Orf, and E. Savory, 2011:  <a href="http://www.sciencedirect.com/science/article/pii/S0167610511000481">Improved modelling of downburst outflows for wind engineering applications using a cooling source approach</a>.  <i>J. Wind Engineering and Industrial Aerodynamics,</i> <b>99,</b> 801-814, doi:10.1016/j.jweia.2011.03.003
    <!--43-->
    <li>Wissmeier, U., and R. K. Smith, 2011:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.819/abstract">Tropical cyclone convection: the effects of ambient vertical vorticity</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>137,</b> 845-857. 
    <!--42-->
    <li>Vermeire, B. C., L. G. Orf, and E. Savory, 2011:  <a href="http://dx.doi.org/10.1016/j.jweia.2011.01.019">A parametric study of downburst line near-surface outflows</a>.  <i>J. Wind Engineering and Industrial Aerodynamics</i>, <b>99,</b> 226-238, doi:10.1016/j.jweia.2011.01.019
    <!--41-->
    <li>Letkewicz, C. E., and M. D. Parker, 2011:  <a href="http://dx.doi.org/10.1175/2011MWR3635.1">Impact of environmental variations on simulated squall lines interacting with terrain</a>.  <i>Mon. Wea. Rev.,</i> <b>139,</b> 3163-3183, doi:10.1175/2011MWR3635.1
    <!--40-->
    <li>Nowotarski, C. J., P. M. Markowski, and Y. P. Richardson, 2011:  <a href="http://dx.doi.org/10.1175/MWR-D-10-05087.1">The characteristics of numerically simulated supercell storms situated over statically stable noundary layers</a>. <i>Mon. Wea. Rev.,</i> <b>139,</b> 3139-3162, doi:10.1175/MWR-D-10-05087.1
    <!--39-->
    <li>Doyle, J. D., S. Gabersek, Q. Jiang, L. Bernardet, J. M. Brown, A. Dornbrack, E. Filaus, V. Grubisic, D. J. Kirshbaum, O. Knoth, S. Koch, J. Schmidli, I. Stiperski, S. Vosper, and S. Zhong, 2011:  <a href="http://dx.doi.org/10.1175/MWR-D-10-05042.1">An intercomparison of T-REX mountain wave simulations and implications for mesoscale predictability</a>. <i>Mon. Wea. Rev.,</i> <b>139,</b> 2811-2831, doi:10.1175/MWR-D-10-05042.1
    <!--38-->
    <li>Markowski, P. M., and N. Dotzek, 2011:  <a href="http://dx.doi.org/10.1016/j.atmosres.2010.12.027">A numerical study of the effects of orography on supercells</a>. <i>Atmos. Res.,</i> <b>100,</b> 457-478, doi:10.1016/j.atmosres.2010.12.027
    <!--37-->
    <li>Kang, S.-L., and G. H. Bryan, 2011:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/MWR-D-10-05037.1">A large eddy simulation study of moist convection initiation over heterogeneous surface fluxes</a>.  <i> Mon. Wea. Rev.,</i> <b>139,</b> 2901-2917.
    <!--36-->
    <li>Hakim, G. J., 2011:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/2010JAS3644.1">The mean state of axisymmetric hurricanes in statistical equilibrium</a>.  <i>J. Atmos. Sci.,</i> <b>68,</b> 1364-1376. 
    <!--35-->
    <li>Kirshbaum, D. J., 2011:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/2010JAS3642.1">Cloud-resolving simulations of deep convection over a heated mountain</a>.  <i>J. Atmos. Sci.,</i> <b>68,</b> 361-378. 
  </ul>
  </details>

  <details>
  <summary>2010</summary>
  <ul>
    <!--34-->
    <li>Wissmeier, U., R. K. Smith, and R. Goler, 2010:  <a href="http://onlinelibrary.wiley.com/doi/10.1002/qj.691/abstract">The formation of a multicell thunderstorm behind a sea-breeze front</a>.  <i>Quart. J. Roy. Meteor. Soc.,</i> <b>136,</b> 2176-2188. 
    <!--33-->
    <li>Wang, W., and M. W. Rotach, 2010:  <a href="http://www.springerlink.com/content/g56024753568t10t/?p=9ea1e0ce66f34548a4e0c45c6ef62ac4&pi=7">Flux footprints over an undulating surface </a>.  <i>Bound.-Layer Meteor.,</i> <b>136,</b> 325-340. 
    <!--32-->
    <li>French, A. J.,, and M. D. Parker, 2010:  <a href="http://journals.ametsoc.org/doi/abs/10.1175/2010JAS3329.1">The response of simulated nocturnal convective systems to a developing low-level jet</a>.  <i>J. Atmos. Sci.,</i> <b>67,</b> 3384-3408.
    <!--31-->
    <li>Miglietta, M. M., and R. Rotunno, 2010:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2010JAS3378.1">Numerical simulations of low-CAPE flows over a mountain ridge</a>.  <i>J. Atmos. Sci.,</i> <b>67,</b> 2391-2401. 
    <!--30-->
    <li>Parker, M. D., 2010:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2009JAS3404.1">Comments on "A comparison of tropical and midlatitude thunderstorm evolution in response to wind shear</a>."  <i>J. Atmos. Sci.,</i> <b>67,</b> 1700-1707.
    <!--29-->
    <li>James, R. P., and P. M. Markowski, 2010:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2009MWR3018.1">A numerical investigation of the effects of dry air aloft on deep convection</a>.  <i>Mon. Wea. Rev.,</i> <b>138,</b> 140-161.</li>
    <!--28-->
    <li>Wang, W., 2010: <a href="http://www.springerlink.com/content/1k81k3phh2j2772x/?p=c15c4c8a54eb49afb86a707bd087b056&pi=10">The influence of topography on single-tower-based carbon flux measurements under unstable conditions: a modeling perspective</a>.  <i>Theor. Appl. Climatol.,</i> <b>99,</b> 125-138.</li>
  </ul>
  </details>

  <details>
  <summary>pre-2010</summary>
    <li>2009</li>
  <ul>
    <!--27-->
    <li>Bryan, G. H., and R. Rotunno, 2009:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2009JAS3038.1">Evaluation of an analytical model for the maximum intensity of tropical cyclones</a>.  <i>J. Atmos. Sci.,</i> <b>66,</b> 3042-3060.</li>
    <!--26-->
    <li>Kirshbaum, D. J., and R. B. Smith, 2009:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2009JAS2990.1">Orographic precipitation in the tropics: large-eddy simulations and theory</a>.  <i>J. Atmos. Sci.,</i> <b>66,</b> 2559-2578.</li>
    <!--25-->
    <li>Wissmeier, U., and R. Goler, 2009:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2009JAS2963.1">A comparison of tropical and mid-latitude thunderstorm evolution in response to wind shear</a>.  <i>J. Atmos. Sci.,</i> <b>66,</b> 2385-2401.</li>
    <!--24-->
    <li>Miglietta, M. M., and R. Rotunno, 2009:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2009JAS2902.1">Numerical simulations of conditionally unstable flows over a mountain ridge</a>.  <i>J. Atmos. Sci.,</i> <b>66,</b> 1865-1885.</li>
    <!--23-->
    <li>Bryan, G. H., and R. Rotunno, 2009:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2008MWR2709.1">The maximum intensity of tropical cyclones in axisymmetric numerical model simulations</a>.  <i>Mon. Wea. Rev.,</i> <b>137,</b> 1770-1789.</li>
    <!--22-->
    <li>Schumacher, R. S., 2009:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2008JAS2856.1">Mechanisms for quasi-stationary behavior in simulated heavy-rain-producing convective systems</a>.  <i>J. Atmos. Sci.,</i> <b>66,</b> 1543-1568.</li>
    <!--21-->
    <li>Kang, S.-L., 2009:  <a href="http://www.springerlink.com/content/g808850w31861258/?p=a06581fe52584d768513a4f28e3c60b6&pi=0">Temporal oscillations in the convective boundary layer forced by mesoscale surface heat-flux variations.</a>  <i>Bound.-Layer Meteor.,</i> <b>132,</b> 59-81.
    <!--20-->
    <li>Wang, W., 2009:  <a href="http://www.springerlink.com/content/k75jk2128t68585g/?p=185502518ce24e76b1a5acc81f318f55&pi=4">The influence of thermally-induced mesoscale circulations on turbulence statistics over an idealized urban area under a zero background wind.</a>  <i>Bound.-Layer Meteor.,</i> <b>131,</b> 403-423.
  </ul>
    
  <li>2008</li>
  <ul>
    <!--19-->
    <li>Smith, J. W., and P. R. Bannon, 2008:<a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2008MWR2343.1">A comparison of compressible and anelastic models of deep dry convection</a>.  <i>Mon. Wea. Rev.,</i> <b>136,</b> 4555-4571.</li>
    <!--18-->
    <li>Kang, S.-L., and K. J. Davis, 2008:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2008JAS2390.1">The effects of mesoscale surface heterogeneity on the fair-weather convective atmospheric boundary layer</a>.  <i>J. Atmos. Sci.,</i> <b>65,</b> 3197-3213.</li>
    <!--17-->
    <li>Schumacher, R. S., and R. H. Johnson, 2008:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2008MWR2471.1">Mesoscale processes contributing to extreme rainfall in a midlatitude warm-season flash flood</a>.  <i>Mon. Wea. Rev.,</i> <b>136,</b> 3964-3986.</li>
    <!--16-->
    <li>Wang, W., and K. J. Davis, 2008:  <a href="http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6V8W-4SWXDB5-1&_user=1518228&_coverDate=09%2F03%2F2008&_rdoc=11&_fmt=high&_orig=browse&_srch=doc-info(%23toc%235881%232008%23998519989%23696525%23FLA%23display%23Volume)&_cdi=5881&_sort=d&_docanchor=&_ct=24&_acct=C000053462&_version=1&_urlVersion=0&_userid=1518228&md5=42809e662b6b4a1f7cac4e89928d7ab5">A numerical study of the influence of a clearcut on eddy-covariance fluxes of CO2 measured above a forest</a>.  <i>Agricultural and Forest Meteor.,</i> <b>148,</b> 1488-1500.</li>
    <!--15-->
    <li>Parker, M. D., 2008:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2007JAS2507.1">Response of simulated squall lines to low-level cooling</a>.  <i>J. Atmos. Sci.,</i> <b>65,</b> 1323-1341.</li>
    <!--14-->
    <li>Edson, A. R., and P. R. Bannon, 2008:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2007JAS2278.1">Nonlinear atmospheric adjustment to momentum forcing</a>.  <i>J. Atmos. Sci.</i>, <b>65,</b> 953-969.
    <!--13-->
    <li>Bryan, G. H., and R. Rotunno, 2008:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2007JAS2443.1">Gravity currents in a deep anelastic atmosphere</a>.  <i>J. Atmos. Sci.,</i> <b>65,</b> 536-556.</li>
  </ul>

  <li>2007</li>
  <ul>
    <!--12-->
    <li>Kirshbaum, D. J., R. Rotunno, and G. H. Bryan, 2007:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2007JAS2335.1">The spacing of orographic rainbands triggered by small-scale topography</a>.  <i>J. Atmos. Sci.,</i> <b>64,</b> 4222-4245.</li>
    <!--11-->
    <li>Lin, W. E., L. G. Orf, E. Savory, and C. Novacco, 2007:  <a href="http://technopress.kaist.ac.kr/journal/container/was/was10_4_container.jsp?order=1#abs1">Proposed large-scale modelling of the transient features of a downburst outflow</a>.  <i>Wind & Structures,</i> <b>10,</b> 315-346.</li>
    <!--10-->
    <li>Kirshbaum, D. J., G. H. Bryan, R. Rotunno, and D. R. Durran, 2007:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2FJAS3924.1">The triggering of orographic rainbands by small-scale topography</a>.  <i>J. Atmos. Sci.,</i> <b>64,</b> 1530-1549.</li>
    <!--9-->
    <li>Bryan, G. H., R. Rotunno, and J. M. Fritsch, 2007:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=1.1175%2FJAS3899.1">Roll circulations in the convective region of a simulated squall line</a>. <i>J. Atmos. Sci.,</i> <b>64,</b> 1249-1266.</li>
  </ul>
  <li>2006</li>
  <ul>
    <!--8-->
    <li>Bannon, P. R., J. M. Chagnon, and R. P. James, 2006:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2FMWR3228.1">Mass conservation and the anelastic approximation</a>.  <i>Mon. Wea. Rev.,</i> <b>134,</b> 2989-3005.
    <!--7-->
    <li>Bryan, G. H., J. C. Knievel, and M. D. Parker, 2006:  <a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2FMWR3226.1">A multimodel assessment of RKW Theory's relevance to squall-line characteristics</a>.  <i>Mon. Wea. Rev.,</i> <b>134,</b> 2772-2792.</li>
    <!--6-->
    <li>James, R. P., P. M. Markowski, and J. M. Fritsch, 2006:  <a href="http://ams.allenpress.com/amsonline/?request=get-abstract&doi=10.1175%2FMWR3109.1">Bow echo sensitivity to ambient moisture and cold pool strength</a>. <i>Mon. Wea. Rev.,</i> <b>134,</b> 950-964.
  </ul>
  <li>2005</li>
  <ul>
    <!--5-->
    <li>Fanelli, P. F., and P. R. Bannon, 2005:  <a href="http://dx.doi.org/10.1175/JAS3517.1">Nonlinear atmospheric adjustment to thermal forcing</a>.  <i>J. Atmos. Sci.</i>, <b>62,</b> 4253-4272.
    <!--4-->
    <li>James, R. P., J. M. Fritsch, and P. M. Markowski, 2005:  <a href="http://dx.doi.org/10.1175/MWR3002.1">Environmental distinctions between cellular and slabular convective lines</a>.  <i>Mon. Wea. Rev.</i>, <b>133,</b> 2669-2691.
    <!--3-->
    <li>Bryan, G. H., 2005:  <a href="http://dx.doi.org/10.1175/MWR2952.1">Spurious convective organization in simulated squall lines owing to moist absolutely unstable layers</a>.  <i>Mon. Wea. Rev.</i>, <b>133,</b> 1978-1997.
  </ul>
  <li>2003</li>
  <ul>
    <!--2-->
    <li>Bryan, G. H., J. C. Wyngaard, and J. M. Fritsch, 2003:  <a href="http://dx.doi.org/10.1175/1520-0493(2003)131<2394:RRFTSO>2.0.CO;2">Resolution requirements for the simulation of deep moist convection</a>.  <i>Mon. Wea. Rev.</i>, <b>131,</b> 2394-2416. 
  </ul>
  <li>2002</li>
  <ul>
    <!--1-->
    <li>Bryan, G. H., and J. M. Fritsch, 2002:  <a href="http://dx.doi.org/10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2">A benchmark simulation for moist nonhydrostatic numerical models</a>.  <i>Mon. Wea. Rev.,</i> <b>130,</b> 2917-2928. 
  </ul>
</ul>
  </details>
</details>

<details>
<summary>Some recent conference papers that use CM1</summary>
<ul>
  <li>S. D. Ditchek, J. Molinari, R. G. Fovell, and K. L. Corbosiero, 2018:  <a href="https://ams.confex.com/ams/33HURRICANE/webprogram/Manuscript/Paper339141/HTMC33_Ditchek_ExtendedAbstract.pdf">The tropical cyclone diurnal cycle in CM1 using an ensemble approach</a>.  <i>33rd Conference on Hurricanes and Tropical Meteorology,</i> Ponte Vedra, FL, Amer. Meteor. Soc., 8C.8.</li>
  <li>Alland, J. J., B. H. Tang, and K. L. Corbosiero, 2018:  <a href="https://ams.confex.com/ams/33HURRICANE/webprogram/Manuscript/Paper339577/Alland_Tropical_2017_Extended_Abstract_FINAL.pdf">The Synergistic Effect of Midlevel Dry Air and Vertical Wind Shear on the Development of the Tropical Cyclone Secondary Circulation</a>. <i>33rd Conference on Hurricanes and Tropical Meteorology,</i> Ponte Vedra, FL, Amer. Meteor. Soc., 8C.1.</li>
  <li>Duran, P., and J. Molinari, 2018:  <a href="https://ams.confex.com/ams/33HURRICANE/webprogram/Manuscript/Paper339247/2018_Duran%26Molinari_AMStropical_extended_abstract.pdf">Upper-Tropospheric Static Stability in Tropical Cyclones: Observations and Modeling</a>. <i>33rd Conference on Hurricanes and Tropical Meteorology,</i> Ponte Vedra, FL, Amer. Meteor. Soc., 13B.1.</li>
  <li>Naylor, J., 2016:  <a href="https://ams.confex.com/ams/28SLS/webprogram/Manuscript/Paper301741/naylor16sls.pdf">Exploring the impact of storm relative helicity on the relationship between cold pools and tornadoes</a>.  <i>28th Conference on Severe Local Storms,</i> Portland, OR, Amer. Meteor. Soc., 6A.1. 
</ul>
</details>

<details>
<summary>Honors and awards</summary>
<ul>
  <li><b>2018</b>:  <u>Patrick Duran</u> of University at Albany, SUNY, won the Max Eaton Award at the 2018 Conference on Hurricanes and Tropical Meteorology for his study <a href="https://ams.confex.com/ams/33HURRICANE/webprogram/Paper339247.html">Upper-Tropospheric Static Stability in Tropical Cyclones: Observations and Modeling</a> 
  <li><b>2018</b>:  <u>Erik Nielsen</u> of Colorado State University won a Student Poster Presentation award at the 2018 Conference on Weather Analysis and Forecasting for his study <a href="https://ams.confex.com/ams/29WAF25NWP/webprogram/Paper345419.html">Dynamical Insights into Extreme Short-Term Precipitation Associated with Supercells and Mesovortices</a>.
  <li><b>2016</b>:  <u>Tristan Kading</u> of University of Connecticut won a Student Poster Presentation award at the 2016 Conference on Severe Local Storms for his study <a href="https://ams.confex.com/ams/28SLS/webprogram/Paper300885.html">Simulated Interaction of an Idealized Squall Line with a Cool Marine Atmospheric Boundary Layer</a>.
  <li><b>2014</b>:  <u>Eli Dennis</u> of Penn State won the Best Student Poster Presentation award at the 2014 Severe Local Storms Conference for his study <a href="https://ams.confex.com/ams/27SLS/webprogram/Paper255162.html">The impact of hodograph shape on hail production in simulated supercell storms</a>.
  <li><b>2012</b>:  <u>Ryan Hastings</u> of Penn State won Best Student Poster Presentation for <a href="https://ams.confex.com/ams/26SLS/webprogram/Manuscript/Paper212524/HastingsSLS2012_2.pdf">Mergers in supercell environments. Part II: Tornadogenesis potential during merger as evaluated by changes in the near-surface low-level mesocyclone</a>.
  <li><b>2012</b>:  <u>Daniel Chavas</u> of MIT won the 2012 Max Eaton Prize from the AMS for his paper, <a href="https://ams.confex.com/ams/30Hurricane/webprogram/Manuscript/Paper204770/AMS_Hurr_2012_CHAVAS_ExtAbs.pdf">Equilibrium tropical cyclone size in an idealized state of axisymmetric radiative-convective equilibrium</a>.
  <li><b>2010</b>:  <u>George Bryan and Richard Rotunno</u> were awarded the Banner I. Miller Award by the AMS for their article "<a href="http://ams.allenpress.com/perlserv/?request=get-abstract&doi=10.1175%2F2008MWR2709.1">The maximum intensity of tropical cyclones in axisymmetric numerical model simulations</a>". 
  <li><b>2009</b>:  An analysis of a supercell thunderstorm by <u>Leigh Orf</u> of Central Michigan University was selected for the cover of the <a href="http://www.casc.org/papers/CASC2010-brochure.pdf">Coalition for Academic Scientific Computation's 2010 brochure</a>.</li>
  <li><b>2008</b>:  <u>Casey Letkewicz</u> of the Department of Marine, Earth, and Atmospheric Sciences at North Carolina State University won the Best Student Poster award at the 24th Conference on Severe Local Storms for her studies of simulated mesoscale convective systems crossing the Appalachian Mountains.</li>
  <li><b>2008</b>:  <u>Chris Nowotarski</u> of the Department of Meteorology at The Pennsylvania State University won the Best Student Oral Presentation award at the 24th Conference on Severe Local Storms for his study of numerically simulated supercells  in varying low-level environmental stability.</li>
</ul>
(Please contact George Bryan if you have something to add to this list.)
</details>

<details>
<summary>Links</summary>
<ul>
  <li><a href="http://cola.gmu.edu/grads/">Grid Analysis and Display System (GrADS)</a> 
  <li><a href="http://www.unidata.ucar.edu/software/netcdf/">NetCDF</a> 
  <li><a href="http://www.vapor.ucar.edu/">Visualization and Analysis Platform for Ocean, Atmosphere, and Solar Researchers (VAPOR)</a>  
</ul>
</details>

Send comments and/or questions about this page to:

```
George H. Bryan
NSF National Center for Atmospheric Research
3090 Center Green Drive
Boulder, CO 80301, USA
```
email: gbryan at ucar dot edu
<hr>
<font color="gray">
The National Center for Atmospheric Research is sponsored by the U.S. National Science Foundation. Any opinions, findings and conclusions or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the U.S. National Science Foundation.
</font>
