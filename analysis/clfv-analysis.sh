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
        RECOFILE="${ROOTDIR}/reco_${ROOTBASE}"
        echo ./clfv-reco.py "${ROOTFILE}" -o "${RECOFILE}"
        ./clfv-reco.py "${ROOTFILE}" -o "${RECOFILE}" &> "${RECOFILE}.log"
    ) &
    PIDS=(${PIDS[@]} $!)
done
wait
(
    RECOFILE="${ROOTDIR}/reco_${ROOTBASE}"
    RECOFILES="${RECOFILE/.root/_*.root}"
    RECOPLOT="${RECOFILE/.root/.pdf}"
    hadd -f "${RECOFILE}" $(ls -v ${RECOFILES})
    echo ./clfv-reco-draw.py "${RECOFILE}"
    ./clfv-reco-draw.py "${RECOFILE}"
) &
wait
