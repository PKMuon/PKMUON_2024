#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
source "$project/rpc_response/scripts/env_ubuntu.sh"
bin="$project/rpc_response/build-ubuntu"
# Unique outputs allow repeated validation without overwriting user data.
validation_dir=$(mktemp -d "$project/rpc_response/results/validation.XXXXXX")
"$bin/rpc_root_check" fixture "$validation_dir/input.root"
"$bin/rpc_digitize" "$validation_dir/input.root" "$validation_dir/hits.root" 1 20261018
"$bin/rpc_root_check" hits "$validation_dir/hits.root"
"$bin/rpc_root_check" segments "$project/rpc_response/results/geant4_smoke.root"
"$bin/rpc_root_check" baseline "$project/rpc_response/results/geant4_baseline.root"
"$bin/rpc_root_check" kernel "$project/rpc_response/results/rpc_linear_reference_v1.root"
echo "Validation files: $validation_dir"
