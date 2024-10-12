#!/bin/bash

set -ev

mkdir -p build
pushd build
cmake ..
make -j
popd

./scan_muemumu.py
ls build/mup_*_mumu_*.mac | xargs basename -a | sort -V | tee lxlogin_run.txt

hep_sub lxlogin_run.sh -argu %{ProcId} \
    -n $[$(wc -l lxlogin_run.txt | egrep -o '[0-9]+')*100] \
    -o lxlogin_run_1.log -e lxlogin_run_2.log
