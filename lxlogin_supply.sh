#!/bin/bash

# Collect running jobs.
echo "Collecting running jobs..."
RUNNING="$(hep_q -u $USER | egrep -o 'lxlogin_run.sh [0-9]+$' | egrep -o '[0-9]+$')"
echo "Done!"

IMAC=0
for MACI in $(cat lxlogin_run.txt); do
    echo "Checking for ${MACI}..."
    for IRUN in $(seq 0 99); do
        ROOT="${MACI/.mac/_${IRUN}.root}"
        RECO="reco_${ROOT}"
        [ -f "build/root_file/${RECO}" ] && continue
        I=$[${IMAC}*100+${IRUN}]
        if egrep -q "^${I}$" <<< "${RUNNING}"; then
            echo "Running reco ${I}: build/root_file/${RECO}"
            continue
        fi
        echo "Missing reco ${I}: build/root_file/${RECO}"
        hep_sub lxlogin_run.sh -argu "${I}" -n 1 \
            -o lxlogin_run_1.log -e lxlogin_run_2.log
    done
    let IMAC+=1
done
