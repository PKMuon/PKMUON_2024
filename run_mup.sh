#!/bin/bash

if [ $# = 0 ] || [ $# -gt 2 ]; then
    1>&2 echo "usage: $(basename "$0") <macfile> [ <ntask> ]"
    exit 1
fi

MAC="$1"
N=$[$2]
[ $N = 0 ] && N=500

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
    MACI="root_file/$(basename "${MAC/.mac/_${I}.mac}")"
    sed 's@/rlt/SetFileName\(.*\)\.root@/rlt/SetFileName\1_'${I}'.root@g' "${MAC}" > "${MACI}"
    echo ./muPos "${MACI}"
    ./muPos "${MACI}" &> "${MACI}.log" &
    PIDS=(${PIDS[@]} $!)
done
wait
