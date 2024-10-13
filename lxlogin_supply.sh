#!/bin/bash

# Collect running jobs.
echo "Collecting running jobs..."
RUNNING="$(hep_q -u $USER | egrep -o 'lxlogin_run.sh [0-9]+$' | egrep -o '[0-9]+$')"
echo "Done!"

# Remove empty outputs.
if [ -z "${RUNNING}" ]; then
    count-events -r -t tree build/root_file/
fi

IMAC=0
for MACI in $(cat lxlogin_run.txt); do
    COMPLETE=1
    RECOS=""
    echo "Checking for ${MACI}..."
    [ -f "build/root_file/reco_${MACI/.mac/.root}" ] && continue
    for IRUN in $(seq 0 99); do
        ROOT="${MACI/.mac/_${IRUN}.root}"
        LOG="${MACI/.mac/_${IRUN}.mac.log}"
        RECO="reco_${ROOT}"
        if [ -f "build/root_file/${RECO}" ]; then
            rm -f "build/root_file/${ROOT}" "build/root_file/${LOG}"
            RECOS="${RECOS} build/root_file/${RECO}"
            continue
        fi
        COMPLETE=0
        I=$[${IMAC}*100+${IRUN}]
        if egrep -q "^${I}$" <<< "${RUNNING}"; then
            echo "Running reco ${I}: build/root_file/${RECO}"
            continue
        fi
        echo "Missing reco ${I}: build/root_file/${RECO}"
        hep_sub lxlogin_run.sh -argu "${I}" -n 1 \
            -o lxlogin_run_1.log -e lxlogin_run_2.log
    done
    if [ "${COMPLETE}" = 1 ]; then
        echo "Generating build/root_file/reco_${MACI/.mac/.root}"
        hadd -f "build/root_file/reco_${MACI/.mac/.root}" ${RECOS} && rm ${RECOS}
    fi
    let IMAC+=1
done
