#!/bin/bash
# Run one short CM1 regression case and compare its stats with the reference.
#
# usage: run_case.sh CASE [--write]
#   CASE     supercell_morrison | supercell_nssl | les_stable
#   --write  write the run's stats as the new reference instead of comparing
#
# Environment (all optional):
#   CM1_EXE   model executable          (default: run/cm1.exe)
#   NPROCS    MPI ranks                 (default: 4)
#   MPIEXEC   MPI launcher; set it empty for a serial build (default: mpiexec)
#   WORKDIR   where run directories go  (default: regression_runs)
#   REFDIR    where reference CSVs live (default: tests/regression/reference)
#
# Exit status: 0 = match (or no reference yet), 1 = stats differ,
#              2 = setup error, or the model did not finish.

set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
top=$(cd "$here/../.." && pwd)

case_name=${1:?usage: run_case.sh CASE [--write]}
mode=${2:-}
exe=${CM1_EXE:-$top/run/cm1.exe}
nprocs=${NPROCS:-4}
mpiexec=${MPIEXEC-mpiexec}
workdir=${WORKDIR:-$top/regression_runs}/$case_name
ref=${REFDIR:-$here/reference}/$case_name.csv

# setnml NAME VALUE: change "NAME = old," in namelist.input
setnml() {
    if ! grep -q -E "^ *$1 *=" namelist.input; then
        echo "error: $1 not found in namelist.input"
        exit 2
    fi
    sed -i -E "s/^( *$1 *= *)[^,]*,/\1 $2,/" namelist.input
}

rm -rf "$workdir"
mkdir -p "$workdir"
cd "$workdir" || exit 2

case $case_name in
    supercell_morrison)
        cp "$top/run/namelist.input" .
        setnml timax 1800.0
        ;;
    supercell_nssl)
        cp "$top/run/namelist.input" .
        setnml timax 1800.0
        setnml ptype 27
        ;;
    les_stable)
        cp "$top/run/config_files/les_StableBoundLayer/namelist.input" .
        setnml nx 32
        setnml ny 32
        setnml nz 32
        setnml timax 600.0
        ;;
    *)
        echo "error: unknown case $case_name"
        exit 2
        ;;
esac
# Stats every minute; no 3-D output (not compared, and slow to write).
setnml statfrq 60.0
setnml tapfrq 99999.0

for f in LANDUSE.TBL RRTMG_LW_DATA RRTMG_SW_DATA; do
    ln -s "$top/run/$f" .
done

start=$SECONDS
status=0
if [ -n "$mpiexec" ]; then
    # $mpiexec is unquoted so a launcher with options ("mpirun --bind-to none") works.
    OMP_NUM_THREADS=1 timeout 1800 $mpiexec -n "$nprocs" "$exe" > cm1.log 2>&1 || status=$?
else
    OMP_NUM_THREADS=1 timeout 1800 "$exe" > cm1.log 2>&1 || status=$?
fi
echo "$case_name: cm1.exe exit status $status after $((SECONDS - start)) s"
if ! grep -q "Program terminated normally" cm1.log; then
    echo "error: $case_name did not finish; end of cm1.log:"
    tail -20 cm1.log
    exit 2
fi

compare=$here/compare_stats.py
if [ "$mode" = "--write" ]; then
    mkdir -p "$(dirname "$ref")"
    python3 "$compare" cm1out_stats.ctl "$ref" --write
elif [ ! -f "$ref" ]; then
    # Until a reference is committed, report and pass rather than fail every PR.
    echo "::warning::no reference for $case_name ($ref); comparison skipped"
else
    python3 "$compare" cm1out_stats.ctl "$ref"
fi
