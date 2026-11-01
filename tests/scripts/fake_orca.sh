#!/usr/bin/env bash
# Fake ORCA for tests: never invokes the real ORCA binary.
#
# ARG 1: path to the ORCA main input (.inp).  Copied into the directory that
# holds it, under the input's basename, the output files of a fixture run --
# exactly what a real ORCA run would have produced.  Which fixture to copy is
# chosen with $FAKE_ORCA_FIXTURE (a directory holding <basename>.out/.gbw/
# .property.json/.json/...).
#
# It mirrors the real runner's contract: invoked by delfin.orca.run_orca as
# `orca <input.inp>` with cwd = the run's working directory, exiting 0 on
# success.  A test drives `run_sp_dft` with this script on PATH so the whole
# DELFIN launch path is exercised without any real ORCA.
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "fake_orca: usage: fake_orca.sh <input.inp>" >&2
  exit 2
fi

inp="$1"
dir="$(cd "$(dirname "$inp")" && pwd)"
base="$(basename "$inp" .inp)"
fix="${FAKE_ORCA_FIXTURE:?FAKE_ORCA_FIXTURE must point at a fixture run directory}"
fix_base="$(basename "$fix")"

if [ ! -d "$fix" ]; then
  echo "fake_orca: FAKE_ORCA_FIXTURE ($fix) is not a directory" >&2
  exit 2
fi

for ext in out gbw property.json json property.txt densities engrad hess opt xyz globalminimum.xyz finalensemble.xyz; do
  if [ -f "$fix/$fix_base.$ext" ]; then
    cp "$fix/$fix_base.$ext" "$dir/$base.$ext"
  fi
done

exit 0
