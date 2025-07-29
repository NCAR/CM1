#-----------------------------------------------------------------------------
#
#  To compile the CM1 codes with this Makefile, use the following command:
#
#        make [argu_list]
#
#  where the optional argument list includes:
#  - USE_OPENMP: turn on OpenMP (default=false)
#  - USE_MPI: turn on MPI (default=false)
#  - USE_DOUBLE: turn on double precision calculation (for PGI/NVHPC only, default=false)
#  - USE_OPENACC: turn on OpenACC (for PGI/NVHPC only, default=false)
#  - USE_NETCDF: turn on netCDF output (default=false)
#  - DEBUG: turn on DEBUG mode (default=false)
#
#  Some example usages are:
#    make USE_OPENMP=true                 # shared memory only
#    make USE_OPENMP=true USE_MPI=true    # shared + distributed memory
#    make USE_OPENACC=true USE_MPI=true   # distributed memory + GPU offloading
#
#  Note that the logical values set to the argument variables are case insensitive.
#
#  If the "FC" env variable is not set on your system, you could also specify it through:
#    make FC=ifort                        # use Intel compiler
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#
#  Change all boolean variables to lowercase
#
#-----------------------------------------------------------------------------
USE_OPENMP_L  = $(shell echo $(USE_OPENMP) | tr A-Z a-z)
USE_MPI_L     = $(shell echo $(USE_MPI) | tr A-Z a-z)
USE_DOUBLE_L  = $(shell echo $(USE_DOUBLE) | tr A-Z a-z)
USE_OPENACC_L = $(shell echo $(USE_OPENACC) | tr A-Z a-z)
USE_NETCDF_L  = $(shell echo $(USE_NETCDF) | tr A-Z a-z)
DEBUG_L       = $(shell echo $(DEBUG) | tr A-Z a-z)
USE_CUDAMEM_L = $(shell echo $(USE_CUDAMEM) | tr A-Z a-z)
GPU_TYPE_L    = $(shell echo $(GPU_TYPE) | tr A-Z a-z)

#-----------------------------------------------------------------------------
#  initialize some options as empty
#-----------------------------------------------------------------------------
DM            =
OMP           =
DP            =
OPTS          =
CPP           =
OUTPUTINC     =
OUTPUTLIB     =
OUTPUTOPT     =
LINKOPTS      =

#-----------------------------------------------------------------------------
#                      NETCDF SECTION
#  (uncomment the following four lines if you want netcdf output capability)
#   (also, make sure the paths to netcdf files are correct for your machine)
#              (NOTE: Don't change lines 3 and 4!)
#  Note:  you may need to remove -DNCFPLUS if you use an old version of netcdf
#-----------------------------------------------------------------------------
ifeq (${USE_NETCDF_L},true)
    OUTPUTINC    += -I$(NETCDFBASE)/include
    OUTPUTLIB    += -L$(NETCDFBASE)/lib
    OUTPUTOPT    += -DNETCDF -DNCFPLUS
    LINKOPTS     += -lnetcdf -lnetcdff
endif

#-----------------------------------------------------------------------------
#                    HARDWARE SECTION
#  Do not comment/uncomment the options below manually.
#  Specify the boolean variables above to switch configurations.
#  If a compiler flag is missing, issue a PR in CM1 Repo.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#  NVHPC compiler
#-----------------------------------------------------------------------------
ifeq (${FC},nvfortran)
    OPTS         += -Mfree -Ktrap=none -Mautoinline -Minline=reshape -Kieee -Mnofma
    CPP          += cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
    ifeq (${DEBUG_L},true)
         OPTS    += -g -O0
         DM      += -D_B4B
    else
         OPTS    += -g -O2
    endif
    ifeq (${USE_DOUBLE_L},true)
         OPTS    += -r8
         DP      += -DDP
    endif
    ifeq (${USE_OPENMP_L},true)
         OPTS    += -mp
         OMP     += -DOPENMP
    endif
    ifeq (${USE_OPENACC_L},true)
        ifeq (${DEBUG_L},true)
            ifeq (${GPU_TYPE_L},a100)
                OPTS += -acc -gpu=cc80,lineinfo,nofma,autocompare,math_uniform
            else
                OPTS += -acc -gpu=cc70,lineinfo,nofma,autocompare,math_uniform
            endif
        else
            # Potentially interesting non-default compiler flags:: 
            #     -Minfo=accel    Provides verbose output for accelerator compilation phase
            #     -Mpcast         Enables the use of PCAST functionality
            ifeq (${GPU_TYPE_L},a100)
                OPTS += -acc -gpu=cc80,lineinfo,nofma,math_uniform -Minfo=accel -cudalib=cutensor 
            else
                OPTS += -acc -gpu=cc70,lineinfo,nofma,math_uniform -Mpcast -Minfo=accel
            endif
        endif
        DM       += -D_OPENACC
    endif
    ifeq (${USE_CUDAMEM_L},true)
        DM       += -D_CUDAMEM
    endif
endif

#-----------------------------------------------------------------------------
#  PGI compiler
#-----------------------------------------------------------------------------
ifeq (${FC},pgf90)
    OPTS         += -Mfree -Ktrap=none -Mautoinline -Minline=reshape -Kieee -Mnofma
    CPP          += cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
    ifeq (${DEBUG_L},true)
         OPTS    += -g -O0
         DM      += -D_B4B
    else
         OPTS    += -O2
    endif
    ifeq (${USE_DOUBLE_L},true)
         OPTS    += -r8
         DP      += -DDP
    endif
    ifeq (${USE_OPENMP_L},true)
         OPTS    += -mp
         OMP     += -DOPENMP
    endif
    ifeq (${USE_OPENACC_L},true)
        ifeq (${DEBUG_L},true)
            OPTS += -acc -gpu=cc70,lineinfo,nofma,autocompare,math_uniform
        else
            OPTS += -acc -gpu=cc70,lineinfo,nofma,math_uniform -Mpcast -Minfo=accel
        endif
        DM       += -D_OPENACC
    endif
endif

#-----------------------------------------------------------------------------
#  Untested settings for the following compilers
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#  Intel compiler
#-----------------------------------------------------------------------------
ifeq (${FC},ifort)
    # replace -xHost with -march=core-avx2 to make sure the same code builds on
    # both Intel and AMD CPUs 
    OPTS         += -g -O2 -march=core-avx2 -ip -assume byterecl -fp-model precise -ftz -init=snan,arrays -traceback -no-fma
    CPP          += cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
    ifeq (${USE_OPENMP_L},true)
        OPTS     += -qopenmp
        OMP      += -DOPENMP
    endif
endif

#-----------------------------------------------------------------------------
#  GNU compiler
#-----------------------------------------------------------------------------
ifeq (${FC},gfortran)
    OPTS         += -ffree-form -ffree-line-length-none -O2 -finline-functions
    CPP          += cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
    ifeq (${USE_OPENMP_L},true)
        OPTS     += -fopenmp
        OMP      += -DOPENMP
    endif
endif

#-----------------------------------------------------------------------------
#  Cray compiler
#-----------------------------------------------------------------------------
ifeq (${FC},ftn)
    OPTS         += -O2 -Ovector2 -Oscalar2 -Othread2
    #LINKOPTS     += -lfast_mv
    CPP          += cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
    ifeq (${USE_OPENMP_L},true)
        OPTS     += -h omp
        OMP      += -DOPENMP
    else
        OPTS     += -h noomp
    endif
    ifeq (${USE_OPENACC_L},true)
        OPTS     += -h acc
        DM       += -D_OPENACC
    endif

endif

#-----------------------------------------------------------------------------
#  Distributed memory (MPI)
#-----------------------------------------------------------------------------
ifeq (${USE_MPI_L},true)
    DM           += -DMPI
    override FC  = mpif90
endif

$(info $$DEBUG_L is [${DEBUG_L}])
$(info $$USE_MPI_L is [${USE_MPI_L}])
$(info $$FC is [${FC}])
$(info $$OPTS is [${OPTS}])
$(info $$CPP is [${CPP}])
$(info $$OMP is [${OMP}])
$(info $$DM is [${DM}])

#-----------------------------------------------------------------------------
#-- You shouldn't need to change anything below here
#-----------------------------------------------------------------------------

SRC   = cm1.F \
	constants.F \
	mersennetwister_mod.F \
	input.F \
	adv.F \
	adv_routines.F \
	anelp.F \
	azimavg.F \
	base.F \
	bc.F \
	cm1libs.F \
	comm.F \
	comm_droplet.F \
        bcast.F \
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
	droplet.F \
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
	restart.F \
	writeout_nc.F \
	radlib3d.F \
	irrad3d.F \
	sorad3d.F \
	radtrns3d.F \
	getcape.F \
	ysu.F \
	sfclay.F \
	sfclayrev.F \
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
	module_sf_mynn.F \
	module_bl_myjpbl.F \
	module_sf_myjsfc.F

#SRC_NVTX = nvtx_mod.F
SRC_NVTX = 

OBJS = $(addsuffix .o, $(basename $(SRC)))
OBJS_OPENACC = $(addsuffix .o, $(basename $(SRC_NVTX)))

FFLAGS = $(OPTS)
ifeq (${USE_OPENACC_L},true)
ifeq (${USE_CUDAMEM_L},true)
FFLAGS = $(OPTS) -L/glade/u/apps/opt/nvhpc/22.2/Linux_x86_64/22.2/cuda/11.6/lib64 -lcudart
endif
endif
AR      = ar cru

.SUFFIXES:
.SUFFIXES:      .F .f90 .o

all : cm1

#			$(FC) $(OBJS) $(FFLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS) -o ../run/cm1.exe
cm1:			$(OBJS)
			$(FC) $(LINKOPTS) $(FFLAGS) $(OBJS) $(OUTPUTINC) $(OUTPUTLIB) -o ../run/cm1.exe
			$(AR) onefile.F $(SRC)
			mv onefile.F ../run

%.f90: %.F
			$(CPP) $(DM) $(DP) $(ADV) $(OUTPUTOPT) $< > $@

$(OBJS): %.o: %.f90
			$(FC) $(FFLAGS) $(OUTPUTINC) -c $<

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
base.o: constants.o input.o bc.o comm.o goddard.o cm1libs.o getcape.o turb.o
bc.o: constants.o input.o
bcast.o: constants.o
cm1.o: constants.o input.o param.o base.o init3d.o misclibs.o solve1.o solve2.o solve3.o diff2.o turb.o statpack.o writeout.o restart.o radiation_driver.o radtrns3d.o domaindiag.o azimavg.o hifrq.o parcel.o init_physics.o init_surface.o mp_driver.o ib_module.o eddy_recycle.o lsnudge.o nvtx_mod.o droplet.o
cm1libs.o: input.o constants.o
comm.o: input.o bc.o
comm_droplet.o: input.o constants.o
diff2.o: constants.o input.o
domaindiag.o: constants.o input.o interp_routines.o cm1libs.o getcape.o sfcphys.o turb.o lsnudge.o writeout_nc.o
eddy_recycle.o: constants.o input.o
goddard.o: constants.o input.o cm1libs.o
hifrq.o: input.o constants.o cm1libs.o adv_routines.o ib_module.o writeout_nc.o adv.o comm.o
ib_module.o: input.o constants.o bc.o comm.o
init3d.o:  constants.o input.o misclibs.o cm1libs.o bc.o comm.o module_mp_nssl_2mom.o poiss.o parcel.o ib_module.o turb.o mersennetwister_mod.o droplet.o
init_physics.o: constants.o input.o sfclay.o sfclayrev.o slab.o radtrns3d.o irrad3d.o goddard.o module_ra_rrtmg_lw.o module_ra_rrtmg_sw.o module_sf_gfdl.o module_sf_mynn.o module_bl_mynn.o module_sf_myjsfc.o
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
module_bl_mynn.o: module_ra_etc.o
module_sf_mynn.o: module_ra_etc.o module_bl_mynn.o
module_bl_myjpbl.o: module_ra_etc.o
module_sf_myjsfc.o: module_ra_etc.o
morrison.o: input.o constants.o
module_mp_jensen_ishmael.o: input.o module_ra_etc.o
mp_driver.o: constants.o input.o misclibs.o kessler.o goddard.o thompson.o lfoice.o morrison.o module_mp_nssl_2mom.o module_mp_p3.o module_mp_jensen_ishmael.o
param.o: constants.o input.o init_terrain.o bc.o comm.o bcast.o thompson.o morrison.o module_mp_nssl_2mom.o goddard.o lfoice.o module_mp_p3.o module_mp_jensen_ishmael.o ib_module.o eddy_recycle.o lsnudge.o
parcel.o: constants.o input.o cm1libs.o misclibs.o bc.o comm.o writeout_nc.o comm_droplet.o
droplet.o: constants.o input.o cm1libs.o bc.o comm.o comm_droplet.o parcel.o mersennetwister_mod.o writeout_nc.o
pdef.o: input.o bc.o comm.o
pdcomp.o: constants.o input.o adv.o poiss.o ib_module.o
poiss.o: input.o singleton.o
radiation_driver.o: constants.o input.o bc.o radtrns3d.o module_ra_etc.o module_ra_rrtmg_lw.o module_ra_rrtmg_sw.o
radtrns3d.o: irrad3d.o sorad3d.o radlib3d.o
restart.o: constants.o input.o writeout_nc.o lsnudge.o goddard.o lfoice.o
sfcphys.o: constants.o input.o cm1libs.o
solve1.o: constants.o input.o bc.o comm.o diff2.o turbtend.o misclibs.o testcase_simple_phys.o ib_module.o eddy_recycle.o lsnudge.o
solve2.o: constants.o input.o bc.o comm.o adv.o adv_routines.o sound.o sounde.o soundns.o soundcb.o anelp.o misclibs.o module_mp_nssl_2mom.o ib_module.o eddy_recycle.o
solve3.o: constants.o input.o bc.o comm.o adv_routines.o misclibs.o parcel.o pdcomp.o lsnudge.o ib_module.o droplet.o maxmin.o
sorad3d.o: radlib3d.o
sound.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
sounde.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o maxmin.o
soundcb.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
soundns.o: constants.o input.o misclibs.o bc.o comm.o ib_module.o
statpack.o: constants.o input.o maxmin.o misclibs.o cm1libs.o writeout_nc.o
testcase_simple_phys.o: constants.o input.o
thompson.o: input.o module_mp_radar.o module_ra_etc.o
turb.o: constants.o input.o bc.o comm.o sfcphys.o sfclay.o sfclayrev.o slab.o oml.o ysu.o cm1libs.o module_sf_gfdl.o module_bl_gfsedmf.o module_sf_mynn.o module_bl_mynn.o module_bl_myjpbl.o module_sf_myjsfc.o turbnba.o misclibs.o ib_module.o
turbtend.o: constants.o input.o cm1libs.o
turbnba.o: constants.o input.o bc.o comm.o
writeout.o: constants.o input.o bc.o comm.o writeout_nc.o misclibs.o getcape.o ib_module.o cm1libs.o
writeout_nc.o: constants.o input.o
mersennetwister_mod.o: constants.o
input.o: constants.o
