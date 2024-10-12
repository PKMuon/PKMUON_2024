#!/bin/bash

if [ $# != 1 ]; then
    echo "usage: $(basename $0) <index>" >&2
    exit 1
fi
IMAC=$[$1 / 100]
IRUN=$[$1 % 100]
MAC="$(head -$[$IMAC+1] lxlogin_run.txt | tail -1)"

set -ve

pushd build
source layout_al.sh
MACI="root_file/$(basename "${MAC/.mac/_${IRUN}.mac}")"
sed 's@/rlt/SetFileName\(.*\)\.root@/rlt/SetFileName\1_'${IRUN}'.root@g' "${MAC}" > "${MACI}"
./muPos "${MACI}" &> "${MACI}.log"
popd

pushd analysis
./clfv-reco.py "../build/${MACI/.mac/.root}"
popd

pushd build
(
echo rm -f "${MACI}"
echo rm -f "${MACI}.log"
echo rm -f "${MACI/.mac/.root}"
) | sh
popd
