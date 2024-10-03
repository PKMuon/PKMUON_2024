#!/bin/bash

N=500
NPROC=$(nproc || sysctl -n hw.logicalcpu || getconf _NPROCESSORS_ONLN)
IPROC=0
PIDS=()
mkdir -p root_file

# 循环运行Geant4模拟，传入不同的i值
for I in $(seq 0 $[N-1]); do
    if [ $IPROC = $NPROC ]; then
        wait $PIDS
        PIDS=(${PIDS[@]:1})
    else
        let IPROC+=1
    fi
    cp mup_50GeV.mac mup_50GeV_$I.mac
    sed -i "s/mup_50GeV\\.root/mup_50GeV_${I}.root/g" mup_50GeV_$I.mac

    # 运行Geant4模拟
    echo ./muPos mup_50GeV_$I.mac
    ./muPos mup_50GeV_$I.mac &> root_file/mup_50GeV_${I}.root.log &
    PIDS=(${PIDS[@]} $!)
done
wait
