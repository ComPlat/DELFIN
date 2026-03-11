#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${ROOT}/bin"
XTB4STDA_RUNTIME_ROOT="${XTB4STDAHOME:-${ROOT}/share/xtb4stda}"

check_xtb4stda_runtime() {
  local required=(
    ".xtb4stdarc"
    ".param_stda1.xtb"
    ".param_stda2.xtb"
  )
  local file

  for file in "${required[@]}"; do
    if [[ ! -f "${XTB4STDA_RUNTIME_ROOT}/${file}" ]]; then
      printf "%-12s MISSING %s\n" "xtb4stda-data" "${XTB4STDA_RUNTIME_ROOT}/${file}"
      return 1
    fi
  done

  printf "%-12s OK      %s\n" "xtb4stda-data" "${XTB4STDA_RUNTIME_ROOT}"
}

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
      "${path}" --version 2>/dev/null | sed -n '1,10p' || true
      ;;
    crest)
      "${path}" --version 2>/dev/null | sed -n '1,12p' || true
      ;;
    orca)
      strings "${path}" | rg -m 1 "Program Version" || true
      ;;
    std2)
      "${path}" 2>&1 | sed -n '1,18p' || true
      ;;
    stda)
      "${path}" 2>&1 | sed -n '1,18p' || true
      ;;
    xtb4stda)
      check_xtb4stda_runtime
      "${path}" 2>&1 | sed -n '1,20p' || true
      ;;
    dftb+)
      "${path}" --version 2>&1 | sed -n '1,16p' || true
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
