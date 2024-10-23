#!/bin/bash

N=1000
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
    cp CryMu.mac CryMu_$I.mac
    sed -i "s/CryMu\\.root/CryMu_${I}.root/g" CryMu_$I.mac

    # 运行Geant4模拟
    echo ./muPos CryMu_$I.mac
    ./muPos CryMu_$I.mac &> root_file/CryMu_${I}.root.log &
    PIDS=(${PIDS[@]} $!)
done
wait
