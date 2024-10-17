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
source layout_"${POSTFIX}".sh
mkdir -p root_file
MACI="root_file/$(basename "${MAC/.mac/_${POSTFIX}_${IRUN}.mac}")"
ROOT="${MACI/.mac/.root}"
sed "s@/rlt/SetFileName.*@/rlt/SetFileName ${ROOT}@g" "${MAC}" > "${MACI}"
./muPos "${MACI}" &> "${MACI}.log"
popd

pushd analysis
RECO="../build/root_file/reco_$(basename "${MAC/.mac/_${POSTFIX}_${IRUN}.root}")"
RENE="../build/root_file/reco_$(basename "${MAC/.mac/_${POSTFIX}_${IRUN}.root}")"
./clfv-reco.py "../build/${ROOT}" -o "${RECO}"
./clfv-reco.py "../build/${ROOT}" -o "${RENE}" -e 1.0
popd

pushd build
(
echo rm -f "${MACI}"
echo rm -f "${MACI}.log"
#echo rm -f "${ROOT}"
) | sh
popd
