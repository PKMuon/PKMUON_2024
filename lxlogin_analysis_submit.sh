#!/bin/bash

if [ $# != 1 ]; then
    echo "usage: $(basename $0) <postfix>" >&2
    exit 1
fi
POSTFIX="$1"

set -ev

N=$(head -1 lxlogin_run.txt)

hep_sub lxlogin_analysis.sh -mem 6000 -argu %{ProcId} "${POSTFIX}" \
    -n $[($(wc -l lxlogin_run.txt | egrep -o '[0-9]+')-1)*${N}] \
    -o lxlogin_analysis_1.log -e lxlogin_analysis_2.log
