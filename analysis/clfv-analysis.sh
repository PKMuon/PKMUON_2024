#!/bin/bash

if [ $# = 0 ] || [ $# -gt 1 ]; then
    1>&2 echo "usage: $(basename "$0") <macfile>"
    exit 1
fi

MAC="$1"

NPROC=$(nproc || sysctl -n hw.logicalcpu || getconf _NPROCESSORS_ONLN)
grep -q ihep.ac.cn <<< "${HOSTNAME}" && NPROC="$(bc <<< "${NPROC}*0.75/1")"
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
        RENEFILE="${ROOTDIR}/rene_${ROOTBASE}"
        echo ./clfv-reco.py "${ROOTFILE}" -o "${RECOFILE}"
        echo ./clfv-reco.py "${ROOTFILE}" -o "${RENEFILE}" -e 1.0
        [ -f "${RECOFILE}" ] || ./clfv-reco.py "${ROOTFILE}" -o "${RECOFILE}" &> "${RECOFILE}.log"
        [ -f "${RENEFILE}" ] || ./clfv-reco.py "${ROOTFILE}" -o "${RENEFILE}" -e 1.0 &> "${RENEFILE}.log"
    ) &
    PIDS=(${PIDS[@]} $!)
done
wait
(
    RECOFILE="${ROOTDIR}/reco_${ROOTBASE}"
    RECOFILES="${RECOFILE/.root/_*.root}"
    RENEFILE="${ROOTDIR}/rene_${ROOTBASE}"
    RENEFILES="${RENEFILE/.root/_*.root}"
    hadd -f "${RECOFILE}" $(ls -v ${RECOFILES}) &
    hadd -f "${RENEFILE}" $(ls -v ${RENEFILES}) &
    wait
    echo ./clfv-reco-draw.py "${RECOFILE}"
    ./clfv-reco-draw.py "${RECOFILE}"
    echo ./clfv-reco-draw.py "${RENEFILE}"
    ./clfv-reco-draw.py "${RENEFILE}"
) &
wait
