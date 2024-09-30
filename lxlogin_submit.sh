#!/bin/bash

if [ $# != 1 ]; then
    echo "usage: $(basename $0) <postfix>" >&2
    exit 1
fi
POSTFIX="$1"

set -ev

mkdir -p build
pushd build
cmake ..
make -j
popd

N=500
echo "${N}" > lxlogin_run.txt
#./scan_muemumu_"${POSTFIX}".py
ls build/mup_*_mumu_*.mac | xargs basename -a | sort -V | tee -a lxlogin_run.txt

hep_sub lxlogin_run.sh -argu %{ProcId} "${POSTFIX}" \
    -n $[($(wc -l lxlogin_run.txt | egrep -o '[0-9]+')-1)*${N}] \
    -o lxlogin_run_1.log -e lxlogin_run_2.log
