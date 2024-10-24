#!/bin/bash

if [ $# != 1 ]; then
    echo "usage: $(basename $0) <postfix>" >&2
    exit 1
fi
POSTFIX="$1"

# Collect running jobs.
echo "Collecting running jobs..."
RUNNING="$(hep_q -u $USER | egrep -o "lxlogin_run.sh [0-9]+ ${POSTFIX}$")"
echo "Done!"

# Remove empty outputs.
if [ -z "${RUNNING}" ]; then
    count-events -r -t tree build/root_file/
fi

IMAC=-1
for MACI in $(cat lxlogin_run.txt); do
    if [ "${IMAC}" = -1 ]; then
        N="${MACI}"; IMAC=0; continue
    fi
    COMPLETE=1
    RECOS=""
    RENES=""
    echo "Checking for ${MACI}..."
    if [ ! -f "build/root_file/reco_${MACI/.mac/_${POSTFIX}.root}" ]; then
        for IRUN in $(seq 0 $[${N}-1]); do
            MAC="${MACI/.mac/_${POSTFIX}_${IRUN}.mac}"
            ROOT="${MACI/.mac/_${POSTFIX}_${IRUN}.root}"
            RECO="reco_${ROOT}"
            if [ -f "build/root_file/${RECO}" ]; then
                RECOS="${RECOS} build/root_file/${RECO}"
                continue
            fi
            COMPLETE=0
            I=$[${IMAC}*${N}+${IRUN}]
            if egrep -q " ${I} " <<< "${RUNNING}"; then
                echo "Running reco ${I}: build/root_file/${RECO}"
                continue
            fi
            echo "Missing reco ${I}: build/root_file/${RECO}"
            #hep_sub lxlogin_run.sh -argu "${I}" -n 1 \
            #    -o lxlogin_run_1.log -e lxlogin_run_2.log
        done
    fi
    if [ ! -f "build/root_file/rene_${MACI/.mac/_${POSTFIX}.root}" ]; then
        for IRUN in $(seq 0 $[${N}-1]); do
            MAC="${MACI/.mac/_${POSTFIX}_${IRUN}.mac}"
            ROOT="${MACI/.mac/_${POSTFIX}_${IRUN}.root}"
            RENE="rene_${ROOT}"
            if [ -f "build/root_file/${RENE}" ]; then
                RENES="${RENES} build/root_file/${RENE}"
                continue
            fi
            COMPLETE=0
            I=$[${IMAC}*${N}+${IRUN}]
            if egrep -q " ${I} " <<< "${RUNNING}"; then
                echo "Running rene ${I}: build/root_file/${RENE}"
                continue
            fi
            echo "Missing rene ${I}: build/root_file/${RENE}"
            #hep_sub lxlogin_run.sh -argu "${I}" -n 1 \
            #    -o lxlogin_run_1.log -e lxlogin_run_2.log
        done
    fi
    if [ "${COMPLETE}" = 1 ]; then
        if [ ! -z "${RECO}" ]; then
            echo "Generating build/root_file/reco_${MACI/.mac/_${POSTFIX}.root}"
            hadd -f "build/root_file/reco_${MACI/.mac/_${POSTFIX}.root}" ${RECOS} && rm ${RECOS}
        fi
        if [ ! -z "${RENE}" ]; then
            echo "Generating build/root_file/rene_${MACI/.mac/_${POSTFIX}.root}"
            hadd -f "build/root_file/rene_${MACI/.mac/_${POSTFIX}.root}" ${RENES} && rm ${RENES}
        fi
    fi
    let IMAC+=1
done
