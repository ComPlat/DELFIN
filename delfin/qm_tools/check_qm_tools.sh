#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${ROOT}/bin"

print_header() {
  printf "\n== %s ==\n" "$1"
}

check_prog() {
  local prog="$1"
  local path="${BIN}/${prog}"

  if [[ ! -x "${path}" ]]; then
    printf "%-12s MISSING\n" "${prog}"
    return
  fi

  printf "%-12s OK      %s\n" "${prog}" "${path}"

  case "${prog}" in
    xtb)
      "${path}" --version 2>/dev/null | sed -n '1,10p'
      ;;
    crest)
      "${path}" --version 2>/dev/null | sed -n '1,12p'
      ;;
    orca)
      strings "${path}" | rg -m 1 "Program Version" || true
      ;;
    std2)
      "${path}" 2>&1 | sed -n '1,18p'
      ;;
    stda)
      "${path}" 2>&1 | sed -n '1,18p'
      ;;
    xtb4stda)
      "${path}" 2>&1 | sed -n '1,20p'
      ;;
    dftb+)
      "${path}" --version 2>&1 | sed -n '1,16p'
      ;;
    *)
      ;;
  esac
}

print_header "Activation"
printf "source %s/env.sh\n" "${ROOT}"

print_header "Programs"
for prog in xtb crest std2 stda xtb4stda dftb+; do
  check_prog "${prog}"
done
