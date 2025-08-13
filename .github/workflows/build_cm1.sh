#!/usr/bin/env bash

set -ex

#----------------------------------------------------------------------------
# environment
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source ${SCRIPTDIR}/build_common.cfg || { echo "cannot locate ${SCRIPTDIR}/build_common.cfg!!"; exit 1; }
#----------------------------------------------------------------------------


echo "building CM1 in $(pwd)"

case "${COMPILER_FAMILY}" in
    "gcc")
        compiler_target="gfortran"
        ;;
    "oneapi")
        compiler_target="intel"
        ;;
    "nvhpc")
        compiler_target="nvfortran"
        ;;
    *)
        echo "ERROR: unrecognized COMPILER_FAMILY=${COMPILER_FAMILY}"!
        exit 1
        ;;
esac


cd src
make FC=${compiler_target} USE_MPI=true USE_OPENACC=true 
cd ../run
export MPIR_CVAR_ENABLE_GPU=1
mpiexec -n 4 ./cm1.exe
#--jobs ${MAKE_J_PROCS:-$(nproc)}
