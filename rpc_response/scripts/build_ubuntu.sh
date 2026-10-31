#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
garfield_prefix=${GARFIELD_INSTALL:-/home/randomuon/Desktop/sim/garfield/install}
geant4_prefix=${GEANT4_INSTALL:-/home/randomuon/Desktop/sim/geant4/geant4-install}
cmake -S "$project/rpc_response" -B "$project/rpc_response/build-ubuntu" \
  -DCMAKE_BUILD_TYPE=Release -DRPC_WITH_ROOT=ON -DRPC_WITH_GARFIELD=ON \
  -DCMAKE_PREFIX_PATH="$garfield_prefix"
cmake --build "$project/rpc_response/build-ubuntu" -j4
ctest --test-dir "$project/rpc_response/build-ubuntu" --output-on-failure
# Separate build; original build/ and its user outputs are not overwritten.
cmake -S "$project" -B "$project/build-rpc" -DCMAKE_BUILD_TYPE=Release \
  -DWITH_GEANT4_UIVIS=ON -DGeant4_DIR="$geant4_prefix/lib/cmake/Geant4"
cmake --build "$project/build-rpc" -j4
