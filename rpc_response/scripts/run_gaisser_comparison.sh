#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
source "$project/rpc_response/scripts/env_ubuntu.sh"
results="$project/rpc_response/results/comparison_20261020"
bin="$project/rpc_response/build-ubuntu"
mode=${1:-original-y}
if [[ "$mode" == original-y ]]; then
  tag=gaisser_100k
elif [[ "$mode" == aligned-z ]]; then
  tag=gaisser_z_100k
  export MUPOS_GAISSER_Z_AXIS=1
elif [[ "$mode" == plane ]]; then
  tag=gaisser_plane_100k
  export MUPOS_GAISSER_PLANE=1
else
  echo 'Usage: run_gaisser_comparison.sh original-y|aligned-z|plane' >&2; exit 2
fi
mkdir -p "$results"
if [[ -e "$results/$tag.root" ]]; then echo "Existing source output protected: $tag.root" >&2; exit 1; fi
macro="$results/$tag.mac"
sed "s|gaisser_100k.root|$tag.root|" "$project/rpc_response/tests/gaisser_comparison.mac" > "$macro"
cd "$project/build-rpc"
MUPOS_ROOT_SEED=20261020 /usr/bin/time -p ./muPos "$macro" > "$results/${tag}_run.log" 2>&1
"$bin/rpc_digitize" "$results/$tag.root" "$results/${tag}_baseline_50ns.root" 50 20261020
"$bin/rpc_run_legacy_analysis" "$results/$tag.root" "$results/${tag}_legacy_analysis.root" > "$results/${tag}_legacy.log" 2>&1
csv="$results/${tag}_csv"
"$bin/rpc_export_comparison" "$results/$tag.root" "$results/${tag}_baseline_50ns.root" "$csv" "$results/${tag}_legacy_analysis.root"
