#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
source "$project/rpc_response/scripts/env_ubuntu.sh"
results="$project/rpc_response/results/comparison_20261020"
mkdir -p "$results/energies"
bin="$project/rpc_response/build-ubuntu/rpc_make_muon_standard"
for item in 0p1:100000000 0p3:300000000 1:1000000000 3:3000000000 10:10000000000 100:100000000000; do
  label=${item%%:*}; energy=${item##*:}
  prefix="$results/energies/muon_${label}GeV"
  if [[ -f "${prefix}_metadata.json" && -f "${prefix}_transport.csv" ]]; then
    echo "Existing completed sample: $label GeV"
    continue
  fi
  "$bin" "$project/rpc_response/results/rpc_transport_v1.rpc" \
    "$project/rpc_response/results/rpc_linear_reference_v1.root" \
    "$project/rpc_response/results/r134a_298K_760Torr_scan.gas" \
    "$prefix" 10000 "$energy" 20261020 > "$prefix.log" 2>&1
  tail -1 "$prefix.log"
done
