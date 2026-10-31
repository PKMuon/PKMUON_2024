#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
source "$project/rpc_response/scripts/env_ubuntu.sh"
base="$project/rpc_response/results/comparison_20261020"
results="$base/spatial_interpolation"
bin="$project/rpc_response/build-ubuntu"
mkdir -p "$results/validation_energies"
cmake --build "$bin" -j2 --target rpc_digitize rpc_export_comparison rpc_validate rpc_make_muon_standard
ctest --test-dir "$bin" --output-on-failure
for factor in 3 5 10; do
  prefix="$results/gaisser_plane_spatial_${factor}sigma"
  if [[ ! -f "$prefix.root" ]]; then
    "$bin/rpc_digitize" "$base/gaisser_plane_100k.root" "$prefix.root" 50 20261021 --spatial-multi "$factor"
  fi
  if [[ ! -f "$prefix.csv" ]]; then
    "$bin/rpc_export_comparison" --hits "$prefix.root" "$prefix.csv"
  fi
done
# Fresh seeds and intermediate energies, separate from all training libraries.
for item in 0p2:200000000:20261101 0p6:600000000:20261102 2:2000000000:20261103 6:6000000000:20261104 30:30000000000:20261105; do
  label=${item%%:*}; tail=${item#*:}; energy=${tail%%:*}; seed=${tail##*:}
  prefix="$results/validation_energies/muon_${label}GeV"
  if [[ ! -f "${prefix}_metadata.json" ]]; then
    "$bin/rpc_make_muon_standard" "$project/rpc_response/results/rpc_transport_v1.rpc" \
      "$project/rpc_response/results/rpc_linear_reference_v1.root" \
      "$project/rpc_response/results/r134a_298K_760Torr_scan.gas" \
      "$prefix" 10000 "$energy" "$seed" > "$prefix.log" 2>&1
  fi
  tail -1 "$prefix.log"
done
