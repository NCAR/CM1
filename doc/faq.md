# CM1:  Frequently Asked Questions

## Answers to Frequently Asked Questions about CM1

**Updated for `cm1r15`, 13 January 2011**

1.  What is CM1?

CM1 is a three-dimensional, time-dependent, non-hydrostatic numerical model that has been developed primarily by George Bryan at The Pennsylvania State University (PSU) (circa 2000-2002) and at the NSF National Center for Atmospheric Research (NSF NCAR) (2003-present).  CM1 is designed primarily for idealized research, particularly for deep precipitating convection (i.e., thunderstorms).  

2.  Who can use CM1?

CM1 is freely available to anybody that wants to use it. 

3.  Is there support for CM1?

Support for CM1 is very limited.  George Bryan maintains the code, and tries to answer as many questions as possible, but it sometimes takes him a long time to respond to some questions.  Stronger support is developing gradually as more users become experts with certain aspects of the modeling system.  Nevertheless, CM1 **requires** knowledge of UNIX and FORTRAN (especially fortran90).  Experience with numerical modeling is very beneficial.  

4.  How many people use CM1?

The list of CM1 users continues to grow.  Currently, CM1 is used by researchers at NSF NCAR, The Pennsylvania State University, The University of Illinois at Urbana-Champaign, Central Michigan University, North Carolina State University, Florida State University, Massachusetts Institute of Technology, University of North Dakota, University of Washington, CNR-ISAC (Italy), University of Munich (Germany), University of Reading (England), and elsewhere.  

5.  Does CM1 conserve mass?

CM1's governing equations conserve total mass in a moist atmosphere.  But, technically:  no, CM1 does not precisely conserve mass to machine accuracy because of technical numerical reasons.  However, because the pressure equation in CM1 retains all terms, it's mass errors are several orders of magnitude smaller than those from other cloud models that integrate pressure equations (e.g., MM5, ARPS, RAMS).   See Bryan and Fritsch (2002, MWR, pg 2917) and Bryan and Rotunno (2009, MWR, pg 1770) for more information.

6.  Does CM1 conserve total energy?

CM1's governing equations conserve total energy in a moist atmosphere.  But, technically:  no, CM1 does not precisely conserve total energy to machine accuracy because of technical numerical reasons.  However, CM1 retains several terms in the internal energy equation that many other numerical models neglect, such as terms associated with the heat content of hydrometeors, and dissipative heating.  To summarize:  CM1 conserves total energy much better than most other modern cloud models.  See Bryan and Fritsch (2002, MWR, pg 2917) and Bryan and Rotunno (2009, MWR, pg 1770) for more information.   

7.  Does CM1 include the "precipitation mass-sink effect"?

Yes.  (as long as neweqts=1,2)

8.  Does CM1 include dissipative heating?

Yes.  (as long as idiss=1)

9.  What makes CM1 different from other models

One main difference was explained in the four previous questions:  CM1 conserves mass and energy better than other modern cloud models (such as ARPS and RAMS for mass and energy, and WRF for energy).  Also, CM1 was designed specifically to do very-large domain simulations using high resolution;  specifically, it has comparatively little memory overhead, which allows the code to be applied to very large problems (i.e., domains of order 10<sup>9</sup> grid points).  Also, CM1 is rare (if not unique) in its ability to use different equation sets, for different applications;  for example, the model can be run using the compressible equations (with three different solvers, depending on application and desired accuracy), but it can also be used with the anelastic equations, and even the incompressible equations.  This capability allows CM1 to be used very efficiently for a broad range of problems that span many scales. 

10.  Can CM1 be used as a large eddy simulation (LES) model?

Yes.  In fact, CM1 was originally designed to do LES of deep, precipitating convection (i.e., thunderstorms).  The code is currently being configured and evaluated for LES of planetary boundary layers, including studies of convective boundary layers (by S.-L. Kang of Texas Tech University) and of marine cumulus (by D. Kirshbaum of U. Reading). 

11.  How fast is CM1 compared to other models?

CM1 has not been rigorously benchmarked against other models.  It has only been compared head-to-head on a few occasions, usually on single-processor machines.  Furthermore, these tests have always been for situations that CM1 was primarily designed for (e.g., idealized cloud modeling).  From these limited tests, we have concluded that CM1 is roughly 2 times faster and uses roughly half as much memory (RAM) as the Advanced Research WRF.  Also, CM1 is roughly 1.5-2 times faster than ARPS (when using the same time step), and uses ~75% of the memory (RAM) required by ARPS. 

12.  What does "CM1" stand for?

Literally:  Cloud Model 1.  The name CM1 is modeled after the name for MM5, which is the Fifth-Generation Pennsylvania State University / National Center for Atmospheric Research Mesoscale Model.  One could say that CM1 is the First-Generation Pennsylvania State University / National Center for Atmospheric Research Cloud Model. 
