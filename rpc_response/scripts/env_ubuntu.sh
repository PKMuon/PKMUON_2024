#!/usr/bin/env bash
# Source this file; installation paths can be overridden by the caller.
rpc_garfield_prefix=${GARFIELD_INSTALL:-/home/randomuon/Desktop/sim/garfield/install}
rpc_geant4_prefix=${GEANT4_INSTALL:-/home/randomuon/Desktop/sim/geant4/geant4-install}
case $- in *u*) rpc_restore_nounset=1; set +u ;; *) rpc_restore_nounset=0 ;; esac
source "$rpc_garfield_prefix/share/Garfield/setupGarfield.sh"
source "$rpc_geant4_prefix/bin/geant4.sh"
if [[ $rpc_restore_nounset == 1 ]]; then set -u; fi
