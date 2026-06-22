#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${ROOT}/bin"
XTB4STDA_RUNTIME_ROOT="${XTB4STDAHOME:-${ROOT}/share/xtb4stda}"
PYTHON_BIN="${PYTHON_BIN:-${PYTHON:-python}}"

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

check_python_module() {
  local label="$1"
  local module="$2"
  local python_path

  if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    printf "%-12s MISSING %s\n" "${label}" "${PYTHON_BIN}"
    return 1
  fi

  python_path="$("${PYTHON_BIN}" -c 'import sys; print(sys.executable)' 2>/dev/null || true)"
  if "${PYTHON_BIN}" -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('${module}') else 1)" >/dev/null 2>&1; then
    printf "%-12s OK      %s (%s)\n" "${label}" "${module}" "${python_path:-${PYTHON_BIN}}"
    return 0
  fi

  printf "%-12s MISSING %s in %s\n" "${label}" "${module}" "${python_path:-${PYTHON_BIN}}"
  return 1
}

print_header() {
  printf "\n== %s ==\n" "$1"
}

check_prog() {
  local prog="$1"
  local path="${BIN}/${prog}"
  local status=0

  if [[ ! -x "${path}" ]]; then
    printf "%-12s MISSING\n" "${prog}"
    return 1
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
      if ! check_xtb4stda_runtime; then
        status=1
      fi
      "${path}" 2>&1 | sed -n '1,20p' || true
      ;;
    dftb+)
      "${path}" --version 2>&1 | sed -n '1,16p' || true
      ;;
    *)
      ;;
  esac

  return "${status}"
}

failures=0

print_header "Activation"
printf "source %s/env.sh\n" "${ROOT}"

print_header "Python"
if ! check_python_module "rdkit" "rdkit"; then
  failures=$((failures + 1))
fi

print_header "Programs"
for prog in xtb crest std2 stda xtb4stda dftb+; do
  if ! check_prog "${prog}"; then
    failures=$((failures + 1))
  fi
done

exit $(( failures > 0 ))
