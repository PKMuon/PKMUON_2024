#!/bin/bash

N=1000
NPROC=$(nproc || sysctl -n hw.logicalcpu || getconf _NPROCESSORS_ONLN)
IPROC=0
PIDS=()
mkdir -p root_file

# 循环运行Geant4模拟，传入不同的i值
for DM_MASS in 1.0e-02 1.0e-01 1.0e+00 1.0e+01 1.0e+02; do
    XS=$(python3 -c "print('%.1e' % (${DM_MASS} * 1.0e-03))")
    for I in $(seq 0 $[N-1]); do
        if [ $IPROC = $NPROC ]; then
            wait $PIDS
            PIDS=(${PIDS[@]:1})
        else
            let IPROC+=1
        fi
        MACFILE="CryMu_DM_${DM_MASS}GeV_XS_${XS}cm2_${I}.mac"
        ROOTFILE="CryMu_DM_${DM_MASS}GeV_XS_${XS}cm2_${I}.root"
        cp CryMu.mac "root_file/${MACFILE}"
        sed -i "s/CryMu\\.root/${ROOTFILE}/g" "root_file/${MACFILE}"
        sed -i "s@/scatter/DM.*@/scatter/DM ${DM_MASS} GeV ${XS} cm2@g" "root_file/${MACFILE}"

        # 运行Geant4模拟
        echo ./muPos "root_file/${MACFILE}"
        ./muPos "root_file/${MACFILE}" &> "root_file/${ROOTFILE}.log" &
        PIDS=(${PIDS[@]} $!)
    done
done
wait
