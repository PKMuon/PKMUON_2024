#!/bin/bash

if [ $# = 0 ] || [ $# -gt 1 ]; then
    1>&2 echo "usage: $(basename "$0") <macfile>"
    exit 1
fi

MAC="$1"

NPROC=$(nproc || sysctl -n hw.logicalcpu || getconf _NPROCESSORS_ONLN)
IPROC=0
PIDS=()
ROOTFILE="$(dirname "${MAC}")/$(grep /rlt/SetFileName "${MAC}" | awk '{print $2;}')"
ROOTDIR="$(dirname "${ROOTFILE}")"
ROOTBASE="$(basename "${ROOTFILE}")"
ROOTFILES="${ROOTFILE/.root/_*.root}"

for THRESHOLD in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0; do
for ROOTFILE in $(ls -v ${ROOTFILES}); do
    if [ $IPROC = $NPROC ]; then
        wait $PIDS
        PIDS=(${PIDS[@]:1})
    else
        let IPROC+=1
    fi
    (
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export VECLIB_MAXIMUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        ROOTBASE="$(basename "${ROOTFILE}")"
        RECOFILE="${ROOTDIR}/reco_${THRESHOLD}MeV_${ROOTBASE}"
        echo ./clfv-reco.py "${ROOTFILE}" -o "${RECOFILE}" -t "${THRESHOLD}"
        ./clfv-reco.py "${ROOTFILE}" -o "${RECOFILE}" -t "${THRESHOLD}" &> "${RECOFILE}.log"
    ) &
    PIDS=(${PIDS[@]} $!)
done
wait
done
