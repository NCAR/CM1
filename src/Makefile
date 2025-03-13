
#SHELL = /bin/csh

#-----------------------------------------------------------------------------
#                      NETCDF SECTION
#  (uncomment the following four lines if you want netcdf output capability)
#   (also, make sure the paths to netcdf files are correct for your machine)
#              (NOTE: Don't change lines 3 and 4!)
#  Note:  you may need to remove -DNCFPLUS if you use an old version of netcdf
#
OUTPUTINC = -I$(NETCDF)/include
OUTPUTLIB = -L$(NETCDF)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#                     HARDWARE SECTION
#-- Choose the appropriate architecture, and uncomment all lines 
#-- in that section.  See "README.compile" for more information.
#-----------------------------------------------------------------------------
#  single processor, Intel compiler 
#      (default for NCAR's casper)
#FC   = ifort
#OPTS = -O3 -ip -assume byterecl -fp-model precise -ftz -no-fma
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#-----------------------------------------------------------------------------
#  multiple processors, distributed memory (MPI), Intel compiler 
#      (default for NCAR's derecho)
FC   = mpif90
OPTS = -O3 -assume byterecl -fp-model precise -ftz -no-fma
CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
DM   = -DMPI
#-----------------------------------------------------------------------------
#  multiple processors, shared memory (OpenMP), Intel compiler 
#      (can be used for NCAR's derecho)
#FC   = mpif90
#OPTS = -O3 -ip -assume byterecl -fp-model precise -ftz -no-fma -qopenmp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#OMP  = -DOPENMP
#-----------------------------------------------------------------------------
#  multiple processors, hybrid distributed (MPI) and shared memory (OpenMP), Intel compiler 
#      (can be used for NCAR's derecho)
#FC   = mpif90
#OPTS = -O3 -ip -assume byterecl -fp-model precise -ftz -no-fma -qopenmp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#DM   = -DMPI -DOPENMP
#-----------------------------------------------------------------------------
#  single processor, NVIDIA (formerly Portland Group) compiler
#FC   = nvfortran
#OPTS = -Mfree -O2 -Ktrap=none -Mautoinline -Minline=reshape
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#-----------------------------------------------------------------------------
#  multiple processors, distributed memory (MPI), NVIDIA (formerly Portland Group) compiler
#FC   = mpif90
#OPTS = -Mfree -O2 -Ktrap=none -Mautoinline -Minline=reshape
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#DM   = -DMPI
#-----------------------------------------------------------------------------
#  multiple processors, shared memory (OpenMP), NVIDIA (formerly Portland Group) compiler
#FC   = nvfortran
#OPTS = -Mfree -O2 -Ktrap=none -Mautoinline -Minline=reshape -mp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#OMP  = -DOPENMP
#-----------------------------------------------------------------------------
#  single processor, GNU compiler
#FC   = gfortran
#OPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#-----------------------------------------------------------------------------
#  multiple processors, distributed memory (MPI), GNU compiler
#FC = mpif90
#OPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions -fallow-argument-mismatch
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#DM = -DMPI
#-----------------------------------------------------------------------------
#  multiple processors, shared memory (OpenMP), GNU compiler
#FC   = gfortran
#OPTS = -ffree-form -ffree-line-length-none -O2 -finline-functions -fopenmp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#OMP  = -DOPENMP
#-----------------------------------------------------------------------------
#  single processor, Cray fortran compiler
#FC   = ftn
#OPTS = -O2 -h noomp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#-----------------------------------------------------------------------------
#  multiple processors, distributed memory (MPI), Cray fortran compiler
#FC   = mpif90
#OPTS = -O2 -h noomp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#DM   = -DMPI
#-----------------------------------------------------------------------------
#  multiple processors, shared memory (OpenMP), Cray fortran compiler
#FC   = ftn
#OPTS = -O2 -h omp
#CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
#DM   = -DMPI
#OMP  = -DOPENMP
#-----------------------------------------------------------------------------
#-- You shouldn't need to change anything below here
#-----------------------------------------------------------------------------

SRC   = constants.F \
	input.F \
	adv.F \
	adv_routines.F \
	anelp.F \
	azimavg.F \
	base.F \
	bc.F \
	cm1.F \
	cm1libs.F \
	comm.F \
	testcase_simple_phys.F \
	diff2.F \
	eddy_recycle.F \
	goddard.F \
	hifrq.F \
	ib_module.F \
	init3d.F \
	init_physics.F \
	init_surface.F \
	init_terrain.F \
	interp_routines.F \
	kessler.F \
	lfoice.F \
	lsnudge.F \
	maxmin.F \
	misclibs.F \
	morrison.F \
	module_mp_jensen_ishmael.F \
	module_mp_nssl_2mom.F \
	module_mp_p3.F \
	mp_driver.F \
	param.F \
	parcel.F \
	pdef.F \
	pdcomp.F \
	poiss.F \
	sfcphys.F \
	singleton.F \
	radiation_driver.F \
	solve1.F \
	solve2.F \
	solve3.F \
	sounde.F \
	sound.F \
	soundns.F \
	soundcb.F \
	statpack.F \
	stopcm1.F \
	thompson.F \
	module_mp_radar.F \
	turb.F \
	turbtend.F \
	turbnba.F \
	domaindiag.F \
	writeout.F \
	restart_write.F \
	restart_read.F \
	writeout_nc.F \
	radlib3d.F \
	irrad3d.F \
	sorad3d.F \
	radtrns3d.F \
	getcape.F \
	sfclay.F \
	slab.F \
	oml.F \
	module_gfs_machine.F \
	module_gfs_physcons.F \
	module_gfs_funcphys.F \
	module_sf_exchcoef.F \
	module_sf_gfdl.F \
	module_bl_gfsedmf.F \
	module_ra_etc.F \
	module_ra_rrtmg_lw.F \
	module_ra_rrtmg_sw.F \
	module_bl_mynn.F \
	module_bl_mynn_common.F \
	module_bl_mynn_wrapper.F \
	module_sf_mynn.F \
	module_bl_myjpbl.F \
	module_sf_myjsfc.F \
	ccpp_kind_types.F \
	bl_ysu.F \
	mp_wsm6.F \
	mp_wsm6_effectRad.F \
	mp_radar.F \
	module_libmassv.F \
	sf_sfclayrev.F \
	cu_ntiedtke.F \
	mmm_physics_wrapper.F

OBJS = $(addsuffix .o, $(basename $(SRC)))

FFLAGS  =  $(OPTS)
AR      = ar cru

.SUFFIXES:
.SUFFIXES:      .F .f90 .o

all : cm1

cm1:			$(OBJS)
			$(FC) $(OBJS) $(FFLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS) -o ../run/cm1.exe
			$(AR) onefile.F $(SRC)
			mv onefile.F ../run

.F.o:
			$(CPP) $(DM) $(DP) $(ADV) $(OUTPUTOPT) $*.F > $*.f90
			$(FC) $(FFLAGS) $(OUTPUTINC) -c $*.f90

code:
			$(AR) onefile.F $(SRC)
			mv onefile.F ../run

clean:
			rm -f *.f90 *.o *.a *.mod

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

adv.o: constants.o input.o pdef.o adv_routines.o ib_module.o
adv_routines.o: input.o constants.o pdef.o comm.o
anelp.o: constants.o input.o misclibs.o bc.o poiss.o
azimavg.o: input.o constants.o cm1libs.o writeout_nc.o comm.o bc.o
base.o: constants.o input.o bc.o comm.o goddard.o cm1libs.o getcape.o
bc.o: constants.o input.o
cm1.o: constants.o input.o param.o base.o init3d.o misclibs.o solve1.o solve2.o solve3.o pdcomp.o diff2.o turb.o statpack.o writeout.o restart_write.o restart_read.o radiation_driver.o radtrns3d.o domaindiag.o azimavg.o hifrq.o parcel.o init_physics.o init_surface.o mp_driver.o ib_module.o eddy_recycle.o lsnudge.o
cm1libs.o: input.o constants.o
comm.o: input.o bc.o
diff2.o: constants.o input.o
domaindiag.o: constants.o input.o interp_routines.o cm1libs.o getcape.o sfcphys.o turb.o lsnudge.o writeout_nc.o testcase_simple_phys.o
eddy_recycle.o: constants.o input.o
goddard.o: constants.o input.o cm1libs.o
hifrq.o: input.o constants.o cm1libs.o adv.o bc.o ib_module.o writeout_nc.o comm.o
ib_module.o: input.o constants.o bc.o comm.o
init3d.o:  constants.o input.o misclibs.o cm1libs.o bc.o comm.o module_mp_nssl_2mom.o poiss.o parcel.o ib_module.o turb.o
init_physics.o: constants.o input.o sfclay.o slab.o radtrns3d.o irrad3d.o goddard.o module_ra_rrtmg_lw.o module_ra_rrtmg_sw.o module_sf_gfdl.o module_sf_mynn.o module_sf_myjsfc.o sf_sfclayrev.o cu_ntiedtke.o
init_surface.o: constants.o input.o oml.o 
init_terrain.o: constants.o input.o bc.o comm.o adv_routines.o
interp_routines.o: constants.o input.o
irrad3d.o: radlib3d.o
kessler.o: constants.o input.o
lfoice.o: input.o
lsnudge.o: constants.o input.o
maxmin.o: input.o
misclibs.o: constants.o input.o goddard.o lfoice.o
module_mp_radar.o: module_ra_etc.o
module_ra_rrtmg_lw.o: module_ra_etc.o
module_ra_rrtmg_sw.o: module_ra_etc.o module_ra_rrtmg_lw.o
module_gfs_funcphys.o: module_gfs_machine.o module_gfs_physcons.o
module_gfs_physcons.o: module_gfs_machine.o
module_sf_gfdl.o: module_gfs_machine.o module_gfs_physcons.o module_gfs_funcphys.o module_sf_exchcoef.o
module_bl_gfsedmf.o: module_gfs_funcphys.o module_gfs_machine.o module_gfs_physcons.o module_sf_gfdl.o
module_bl_mynn.o: module_bl_mynn_common.o
module_bl_mynn_common.o: module_ra_etc.o module_gfs_machine.o
module_bl_mynn_wrapper.o: module_bl_mynn_common.o module_bl_mynn.o
module_sf_mynn.o: module_ra_etc.o
module_bl_myjpbl.o: module_ra_etc.o
module_sf_myjsfc.o: module_ra_etc.o
morrison.o: input.o constants.o
module_mp_jensen_ishmael.o: input.o module_ra_etc.o
mp_driver.o: constants.o input.o misclibs.o kessler.o goddard.o thompson.o lfoice.o morrison.o module_mp_nssl_2mom.o module_mp_p3.o module_mp_jensen_ishmael.o mmm_physics_wrapper.o mp_wsm6.o
param.o: constants.o input.o init_terrain.o bc.o comm.o thompson.o morrison.o module_mp_nssl_2mom.o goddard.o lfoice.o module_mp_p3.o module_mp_jensen_ishmael.o ib_module.o eddy_recycle.o lsnudge.o mp_wsm6.o ccpp_kind_types.o
parcel.o: constants.o input.o cm1libs.o bc.o comm.o writeout_nc.o
pdef.o: input.o bc.o comm.o
pdcomp.o: constants.o input.o adv.o poiss.o ib_module.o
poiss.o: input.o singleton.o
radiation_driver.o: constants.o input.o bc.o radtrns3d.o module_ra_etc.o module_ra_rrtmg_lw.o module_ra_rrtmg_sw.o
radtrns3d.o: irrad3d.o sorad3d.o radlib3d.o
restart_write.o: constants.o input.o writeout_nc.o lsnudge.o
restart_read.o: constants.o input.o writeout_nc.o lsnudge.o goddard.o lfoice.o restart_write.o
sfcphys.o: constants.o input.o cm1libs.o
solve1.o: constants.o input.o bc.o diff2.o turbtend.o misclibs.o testcase_simple_phys.o eddy_recycle.o lsnudge.o
solve2.o: constants.o input.o bc.o comm.o adv.o sound.o sounde.o soundns.o soundcb.o anelp.o misclibs.o module_mp_nssl_2mom.o ib_module.o
solve3.o: constants.o input.o bc.o comm.o adv_routines.o misclibs.o parcel.o lsnudge.o
sorad3d.o: radlib3d.o
sound.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
sounde.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
soundcb.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
soundns.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
statpack.o: constants.o input.o maxmin.o misclibs.o cm1libs.o writeout_nc.o
testcase_simple_phys.o: constants.o input.o
thompson.o: input.o module_mp_radar.o module_ra_etc.o
turb.o: constants.o input.o bc.o comm.o sfcphys.o sfclay.o slab.o oml.o cm1libs.o module_sf_gfdl.o module_bl_gfsedmf.o module_sf_mynn.o module_bl_mynn_wrapper.o module_bl_myjpbl.o module_sf_myjsfc.o turbnba.o misclibs.o ib_module.o turbtend.o mmm_physics_wrapper.o
turbtend.o: constants.o input.o cm1libs.o
turbnba.o: constants.o input.o bc.o comm.o
writeout.o: constants.o input.o bc.o comm.o writeout_nc.o misclibs.o getcape.o ib_module.o cm1libs.o sfcphys.o eddy_recycle.o
writeout_nc.o: constants.o input.o

bl_ysu.o: ccpp_kind_types.o
mp_wsm6.o: ccpp_kind_types.o module_libmassv.o mp_radar.o
mp_wsm6_effectRad.o: ccpp_kind_types.o mp_wsm6.o
mp_radar.o: ccpp_kind_types.o
sf_sfclayrev.o: ccpp_kind_types.o
cu_ntiedtke.o: ccpp_kind_types.o
mmm_physics_wrapper.o: ccpp_kind_types.o bl_ysu.o mp_wsm6.o mp_wsm6_effectRad.o sf_sfclayrev.o cu_ntiedtke.o

