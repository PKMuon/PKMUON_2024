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
        ROOTBASE="$(basename "${ROOTFILE}")"
        ANAFILE="${ROOTDIR}/ana_${ROOTBASE}"
        POCAFILE="${ROOTDIR}/poca_${ROOTBASE}"
        echo root -l -q -b 'src/analysis.C("'"${ROOTFILE}"'", "'"${ANAFILE}"'")'
        root -l -q -b 'src/analysis.C("'"${ROOTFILE}"'", "'"${ANAFILE}"'")' &> "${ANAFILE}.log"
        echo root -l -q -b 'src/PoCA_sim.C("'"${ANAFILE}"'", "'"${POCAFILE}"'")'
        root -l -q -b 'src/PoCA_sim.C("'"${ANAFILE}"'", "'"${POCAFILE}"'")' &> "${POCAFILE}.log"
    ) &
    PIDS=(${PIDS[@]} $!)
done
wait
(
    POCAFILE="${ROOTDIR}/poca_${ROOTBASE}"
    POCAFILES="${POCAFILE/.root/_*.root}"
    POCAPLOT="${POCAFILE/.root/.pdf}"
    hadd -f "${POCAFILE}" $(ls -v ${POCAFILES})
    root -l -q -b 'src/draw_sim.C("'"${POCAFILE}"'", "'"${POCAPLOT}"'")' &
    wait
) &
wait
