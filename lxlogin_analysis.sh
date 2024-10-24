#!/bin/bash

if [ $# != 2 ]; then
    echo "usage: $(basename $0) <index> <postfix>" >&2
    exit 1
fi
N=$(head -1 lxlogin_run.txt)
IMAC=$[$1 / ${N}]
IRUN=$[$1 % ${N}]
MAC="$(head -$[$IMAC+2] lxlogin_run.txt | tail -1)"
POSTFIX="$2"

set -ve

mkdir -p build
pushd build
NEVENT=0
for I in $(seq 0 $[${N}-1]); do
    ROOTI="root_file/$(basename "${MAC/.mac/_${POSTFIX}_${I}.root}")"
    [ -f "${ROOTI}" ]
    NEVT="$(grep -Eo 'Number of events processed : [0-9]+' ${ROOTI/.root/.mac.log} | grep -Eo '[0-9]+')"
    [ ! -z "${NEVT}" ]
    let NEVENT+="NEVT"
done
echo "NEVENT : ${NEVENT}"
popd

pushd analysis
ROOT="root_file/$(basename "${MAC/.mac/_${POSTFIX}_${IRUN}.root}")"
RECO="root_file/reco_$(basename "${MAC/.mac/_${POSTFIX}_${IRUN}.root}")"
RENE="root_file/rene_$(basename "${MAC/.mac/_${POSTFIX}_${IRUN}.root}")"
&> ../build/"${RECO}".log ./clfv-reco.py -n "${NEVENT}" ../build/"${ROOT}" -o ../build/"${RECO}"
&> ../build/"${RENE}".log ./clfv-reco.py -n "${NEVENT}" ../build/"${ROOT}" -o ../build/"${RENE}" -e 1.0
popd
