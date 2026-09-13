#!/usr/bin/env bash
# ============================================================================
# DELFIN installer for any Linux login node or workstation.
#
# Installs DELFIN into a Python 3.10/3.11 venv, wires up ORCA and the OpenMPI
# it was built against, installs the external tools, and writes the settings
# and the shell environment. Nothing in it belongs to one site or one person:
# no module names, no host names, no user paths. Every location follows from
# $HOME or the command line, and running it again updates and repairs.
#
# Usage:
#   bash install_delfin.sh                     standard profile
#   bash install_delfin.sh --all               everything, ML stacks included
#   bash install_delfin.sh --only crest,gxtb,ketcher
#   bash install_delfin.sh --dry-run           print the plan, change nothing
#   bash install_delfin.sh --update            pull DELFIN, update it and every
#                                              installed tool
#   bash install_delfin.sh --repair            check everything and fix what is
#                                              broken
#
# Profiles:
#   core       DELFIN with the agent, ORCA wiring, OpenMPI
#   standard   core + xtb, g-xTB, MOPAC, CREST, DFTB+, xtb4stda/sTDA, std2
#              + CENSO/anmr, cclib, morfeus, nglview, packmol + Ketcher
#   all        standard + every MLP, CSP and AI tool (several GB)
#
# Multiwfn is licensed and is only installed when named: --only multiwfn.
# Whatever is left out can be installed later -- run this again with --only,
# press Install in the dashboard's Settings, or let DELFIN fetch a tool the
# moment a calculation needs it. All of these go through delfin/installer.py,
# the one list of what DELFIN installs.
#
# Options:
#   --profile core|standard|all   (or --core, --standard, --all)
#   --only LIST       only these tools, comma-separated, on top of core
#   --prefix DIR      where ORCA, OpenMPI and a fetched Python go
#                     (default: $HOME/software)
#   --repo DIR        the DELFIN checkout (default: the one this script is
#                     in, else PREFIX/delfin, cloned when missing)
#   --orca PATH       an ORCA directory or ORCA tarball
#   --no-orca         do not look for ORCA
#   --no-openmpi      do not set up OpenMPI
#   --python PATH     a Python 3.10/3.11 to build the venv from
#   --extras LIST     pip extras for DELFIN (default: agent,docs)
#   --dry-run         print what would be done and change nothing
#   --list            print every tool that can be installed, by group
#   --update          update instead of install (--only limits it)
#   --repair          repair instead of install (--only limits it)
#   -y, --yes         do not ask before the "all" profile
#
# ORCA cannot be downloaded for you (license). Pass --orca with an installed
# ORCA directory or the tarball from https://orcaforum.kofo.mpg.de, or leave
# either in $HOME or the prefix directory, where it is found.
# ============================================================================
set -euo pipefail

log()  { printf '[delfin-install] %s\n' "$*"; }
warn() { printf '[delfin-install] WARNING: %s\n' "$*" >&2; }
die()  { printf '[delfin-install] ERROR: %s\n' "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

usage() {
  awk 'NR > 2 && /^# =+$/ { exit } NR > 2 { sub(/^# ?/, ""); print }' "${BASH_SOURCE[0]}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Options; the environment variables work as well ------------------------
PROFILE="${DELFIN_INSTALL_PROFILE:-standard}"
ONLY="${DELFIN_INSTALL_ONLY:-}"
DRY_RUN=0
LIST_ONLY=0
MODE=install
ASSUME_YES="${DELFIN_ASSUME_YES:-0}"
PREFIX="${DELFIN_PREFIX:-$HOME/software}"
DELFIN_REPO="${DELFIN_REPO:-}"
DELFIN_REPO_URL="${DELFIN_REPO_URL:-https://github.com/ComPlat/DELFIN.git}"
DELFIN_EXTRAS="${DELFIN_EXTRAS-agent,docs}"
ORCA_SOURCE="${ORCA_DIR:-${ORCA_TARBALL:-}}"
ORCA_EXPLICIT=0
[ -z "$ORCA_SOURCE" ] || ORCA_EXPLICIT=1
WANT_ORCA=1
WANT_OPENMPI=1
# The OpenMPI current ORCA releases are linked against, used when no ORCA names one.
DEFAULT_OPENMPI="${DELFIN_OPENMPI_VERSION:-4.1.8}"
OPENMPI_PREFIX="${OPENMPI_PREFIX:-}"
PYTHON_CHOICE="${DELFIN_PYTHON:-}"
DELFIN_CALC_DIR="${DELFIN_CALC_DIR:-$HOME/calc}"
DELFIN_ARCHIVE_DIR="${DELFIN_ARCHIVE_DIR:-$HOME/archive}"
MICROMAMBA_URL="${MICROMAMBA_URL:-https://micro.mamba.pm/api/micromamba/linux-64/latest}"
if [ -z "${MAKE_JOBS:-}" ]; then
  # A login node is shared with everybody else on it: eight compile jobs at most.
  MAKE_JOBS="$(nproc 2>/dev/null || echo 4)"
  if [ "$MAKE_JOBS" -gt 8 ]; then MAKE_JOBS=8; fi
fi

while [ $# -gt 0 ]; do
  case "$1" in
    --profile)   [ $# -ge 2 ] || die "--profile needs core, standard or all"; PROFILE="$2"; shift 2 ;;
    --profile=*) PROFILE="${1#*=}"; shift ;;
    --core|--standard|--all) PROFILE="${1#--}"; shift ;;
    --only)      [ $# -ge 2 ] || die "--only needs a comma-separated list of tools"; ONLY="$2"; shift 2 ;;
    --only=*)    ONLY="${1#*=}"; shift ;;
    --prefix)    [ $# -ge 2 ] || die "--prefix needs a directory"; PREFIX="$2"; shift 2 ;;
    --prefix=*)  PREFIX="${1#*=}"; shift ;;
    --repo)      [ $# -ge 2 ] || die "--repo needs a directory"; DELFIN_REPO="$2"; shift 2 ;;
    --repo=*)    DELFIN_REPO="${1#*=}"; shift ;;
    --orca)      [ $# -ge 2 ] || die "--orca needs an ORCA directory or tarball"; ORCA_SOURCE="$2"; ORCA_EXPLICIT=1; shift 2 ;;
    --orca=*)    ORCA_SOURCE="${1#*=}"; ORCA_EXPLICIT=1; shift ;;
    --no-orca)   WANT_ORCA=0; shift ;;
    --no-openmpi) WANT_OPENMPI=0; shift ;;
    --python)    [ $# -ge 2 ] || die "--python needs a path"; PYTHON_CHOICE="$2"; shift 2 ;;
    --python=*)  PYTHON_CHOICE="${1#*=}"; shift ;;
    --extras)    [ $# -ge 2 ] || die "--extras needs a list (may be empty)"; DELFIN_EXTRAS="$2"; shift 2 ;;
    --extras=*)  DELFIN_EXTRAS="${1#*=}"; shift ;;
    --dry-run)   DRY_RUN=1; shift ;;
    --list)      LIST_ONLY=1; shift ;;
    --update)    MODE=update; shift ;;
    --repair)    MODE=repair; shift ;;
    -y|--yes)    ASSUME_YES=1; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) die "unknown option: $1 (see --help)" ;;
  esac
done

absolute() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s\n' "$PWD/$1" ;; esac; }
PREFIX="$(absolute "$PREFIX")"
[ -z "$DELFIN_REPO" ] || DELFIN_REPO="$(absolute "$DELFIN_REPO")"

# ---- What can be installed: delfin/installer.py is the list ------------------
PACKAGE_PARENT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PLAN=""
DONE=(); FAILED=()
VENV_DIR=""; VENV_PY=""; BASE_PYTHON=""; ORCA_FOUND=""; OMPI_FOUND=""

# The list, read with whatever Python is at hand -- the module needs nothing
# but the standard library, so this works before DELFIN's venv exists.
catalog() {
  local py=""
  if [ -n "$VENV_PY" ] && [ -x "$VENV_PY" ]; then
    py="$VENV_PY"
  elif [ -n "$VENV_DIR" ] && [ -x "$VENV_DIR/bin/python" ]; then
    py="$VENV_DIR/bin/python"
  else
    py="$(command -v python3 2>/dev/null || true)"
  fi
  [ -n "$py" ] || die "python3 is needed to read DELFIN's list of tools"
  PYTHONPATH="$PACKAGE_PARENT${PYTHONPATH:+:$PYTHONPATH}" "$py" -m delfin.installer "$@"
}

select_tools() {
  if [ -n "$ONLY" ]; then
    PLAN="$(catalog --plan "$ONLY")" || die "--only names a tool DELFIN does not know (see --list)"
  elif [ "$MODE" = install ]; then
    case "$PROFILE" in core|standard|all) ;; *) die "unknown profile: $PROFILE (core, standard or all)" ;; esac
    PLAN="$(catalog --profile "$PROFILE")"
  fi
}

planned_tools() {
  # shellcheck disable=SC2046
  echo $(printf '%s\n' "$PLAN" | cut -s -d: -f2-)
}

confirm_all() {
  if [ "$MODE" != install ] || [ "$PROFILE" != all ] || [ -n "$ONLY" ] || [ "$ASSUME_YES" = 1 ] || [ "$DRY_RUN" = 1 ] || [ ! -t 0 ]; then
    return 0
  fi
  printf '[delfin-install] "all" includes the ML stacks (PyTorch and friends): several GB and up to an hour. Continue? [y/N] '
  local answer=""
  read -r answer || answer=""
  case "$answer" in y|Y|yes|YES) ;; *) die "stopped before changing anything" ;; esac
}

fetch() {
  if have curl; then
    curl -fL --retry 3 -o "$2" "$1"
  elif have wget; then
    wget -O "$2" "$1"
  else
    warn "neither curl nor wget is available to download $1"
    return 1
  fi
}

# ---- DELFIN checkout -----------------------------------------------------------
resolve_repo() {
  if [ -z "$DELFIN_REPO" ]; then
    if [ -f "$SCRIPT_DIR/../../pyproject.toml" ] && [ -d "$SCRIPT_DIR/../../delfin" ]; then
      DELFIN_REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
    else
      DELFIN_REPO="$PREFIX/delfin"
    fi
  fi
  VENV_DIR="$DELFIN_REPO/.venv"
}

ensure_repo() {
  if [ -f "$DELFIN_REPO/pyproject.toml" ] && grep -q 'name = "delfin-complat"' "$DELFIN_REPO/pyproject.toml"; then
    log "DELFIN checkout: $DELFIN_REPO"
    return 0
  fi
  if [ -e "$DELFIN_REPO" ] && [ -n "$(ls -A "$DELFIN_REPO" 2>/dev/null)" ]; then
    die "$DELFIN_REPO exists but is not a DELFIN checkout; pass --repo DIR"
  fi
  have git || die "git is needed to fetch DELFIN"
  if [ "$DRY_RUN" = 1 ]; then
    log "would clone $DELFIN_REPO_URL into $DELFIN_REPO"
    return 0
  fi
  log "cloning $DELFIN_REPO_URL into $DELFIN_REPO"
  mkdir -p "$(dirname "$DELFIN_REPO")"
  git clone "$DELFIN_REPO_URL" "$DELFIN_REPO"
}

# ---- Python and the venv --------------------------------------------------------
python_ok() {
  "$1" -c 'import sys; sys.exit(0 if (3, 10) <= sys.version_info[:2] < (3, 12) else 1)' >/dev/null 2>&1
}

find_micromamba() {
  local candidate
  for candidate in "${MAMBA_EXE:-}" "$(command -v micromamba 2>/dev/null || true)" \
      "$PREFIX/micromamba/bin/micromamba" "$HOME/.delfin/qm_tools/bin/micromamba" \
      "$HOME/micromamba/bin/micromamba" "$HOME/.local/bin/micromamba"; do
    if [ -n "$candidate" ] && [ -x "$candidate" ]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

# A Python of DELFIN's own, for a machine that has no 3.10/3.11 at all. No root,
# no module system: micromamba is one static binary.
bootstrap_python() {
  local env_dir="$PREFIX/delfin-python" mamba work
  if [ -x "$env_dir/bin/python" ] && python_ok "$env_dir/bin/python"; then
    BASE_PYTHON="$env_dir/bin/python"
    return 0
  fi
  if ! mamba="$(find_micromamba)"; then
    have curl || die "no Python 3.10/3.11 here and no curl to fetch one; pass --python PATH"
    log "no Python 3.10/3.11 here; fetching micromamba (about 10 MB)"
    mkdir -p "$PREFIX/micromamba/bin"
    work="$(mktemp -d "$PREFIX/.micromamba.XXXXXX")"
    if ! curl -fsSL "$MICROMAMBA_URL" | tar -xj -C "$work" bin/micromamba; then
      rm -rf "$work"
      die "micromamba could not be fetched from $MICROMAMBA_URL; pass --python PATH"
    fi
    install -m 755 "$work/bin/micromamba" "$PREFIX/micromamba/bin/micromamba"
    rm -rf "$work"
    mamba="$PREFIX/micromamba/bin/micromamba"
  fi
  log "creating Python 3.11 at $env_dir"
  MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$PREFIX/micromamba}" \
    "$mamba" create -y -p "$env_dir" -c conda-forge "python=3.11" pip
  BASE_PYTHON="$env_dir/bin/python"
}

ensure_venv() {
  local reuse=0 candidate resolved aside spec
  spec="$DELFIN_REPO${DELFIN_EXTRAS:+[$DELFIN_EXTRAS]}"
  if [ -x "$VENV_DIR/bin/python" ] && python_ok "$VENV_DIR/bin/python"; then
    reuse=1
  else
    for candidate in "$PYTHON_CHOICE" python3.11 python3.10 python3; do
      [ -n "$candidate" ] || continue
      resolved="$(command -v "$candidate" 2>/dev/null || true)"
      if [ -n "$resolved" ] && python_ok "$resolved"; then
        BASE_PYTHON="$resolved"
        break
      fi
    done
  fi

  if [ "$DRY_RUN" = 1 ]; then
    if [ "$reuse" = 1 ]; then
      log "would update the venv at $VENV_DIR ($("$VENV_DIR/bin/python" -V 2>&1))"
    elif [ -n "$BASE_PYTHON" ]; then
      log "would create a venv at $VENV_DIR from $BASE_PYTHON ($("$BASE_PYTHON" -V 2>&1))"
    else
      log "would fetch micromamba and create Python 3.11 at $PREFIX/delfin-python for the venv at $VENV_DIR"
    fi
    log "would pip install -e \"$spec\""
    return 0
  fi

  if [ "$reuse" != 1 ]; then
    [ -n "$BASE_PYTHON" ] || bootstrap_python
    if [ -e "$VENV_DIR" ]; then
      aside="$VENV_DIR.before-$(date +%Y%m%d-%H%M%S)"
      warn "$VENV_DIR is not a Python 3.10/3.11 venv; moving it aside to $aside"
      mv "$VENV_DIR" "$aside"
    fi
    log "creating the venv at $VENV_DIR from $BASE_PYTHON"
    if ! "$BASE_PYTHON" -m venv "$VENV_DIR"; then
      # Some distributions ship a python that cannot make a venv (no ensurepip).
      rm -rf "$VENV_DIR"
      warn "$BASE_PYTHON cannot create a venv; using a Python of DELFIN's own"
      BASE_PYTHON=""
      bootstrap_python
      "$BASE_PYTHON" -m venv "$VENV_DIR"
    fi
  fi
  VENV_PY="$VENV_DIR/bin/python"
  log "installing DELFIN into $VENV_DIR"
  "$VENV_PY" -m pip install --upgrade pip wheel
  "$VENV_PY" -m pip install -e "$spec"
}

# ---- ORCA and its OpenMPI -------------------------------------------------------------
# An ORCA directory, not just any "orca": /usr/bin/orca is a screen reader.
is_orca_dir() {
  [ -x "$1/orca" ] && { [ -x "$1/orca_plot" ] || [ -x "$1/orca_mdci" ]; }
}

extract_orca_tarball() {
  local tarball="$1" top
  top="$({ tar -tf "$tarball" 2>/dev/null || true; } | head -n 1 | cut -d/ -f1)"
  if [ -z "$top" ] || [ "$top" = "." ]; then
    warn "cannot read the ORCA tarball $tarball"
    return 1
  fi
  ORCA_FOUND="$PREFIX/$top"
  if is_orca_dir "$ORCA_FOUND"; then
    return 0
  fi
  if [ "$DRY_RUN" = 1 ]; then
    log "would extract $tarball into $PREFIX"
    return 0
  fi
  log "extracting ORCA from $tarball into $PREFIX"
  mkdir -p "$PREFIX"
  tar -xf "$tarball" -C "$PREFIX"
  if ! is_orca_dir "$ORCA_FOUND"; then
    warn "$tarball did not unpack into an ORCA directory"
    ORCA_FOUND=""
    return 1
  fi
}

settings_orca() {
  local py="$VENV_PY"
  if [ -z "$py" ] || [ ! -x "$py" ]; then py="$(command -v python3 2>/dev/null || true)"; fi
  [ -n "$py" ] || return 0
  "$py" - <<'PY' 2>/dev/null || true
import json
import pathlib

try:
    runtime = json.loads((pathlib.Path.home() / ".delfin_settings.json").read_text()).get("runtime") or {}
except Exception:
    runtime = {}
for value in (runtime.get("orca_base"), (runtime.get("slurm") or {}).get("orca_base"),
              (runtime.get("local") or {}).get("orca_base")):
    if value:
        print(pathlib.Path(value).expanduser())
        break
PY
}

find_orca() {
  ORCA_FOUND=""
  if [ "$WANT_ORCA" != 1 ]; then
    log "ORCA: not looked for (--no-orca)"
    return 0
  fi
  local candidate
  if [ -n "$ORCA_SOURCE" ]; then
    if [ -d "$ORCA_SOURCE" ] && is_orca_dir "$ORCA_SOURCE"; then
      ORCA_FOUND="$(cd "$ORCA_SOURCE" && pwd)"
      return 0
    elif [ -f "$ORCA_SOURCE" ]; then
      if extract_orca_tarball "$ORCA_SOURCE"; then return 0; fi
    else
      warn "$ORCA_SOURCE is neither an ORCA directory nor an ORCA tarball"
    fi
  fi
  candidate="$(settings_orca)"
  if [ -n "$candidate" ] && is_orca_dir "$candidate"; then
    ORCA_FOUND="$candidate"
    return 0
  fi
  candidate="$(command -v orca 2>/dev/null || true)"
  if [ -n "$candidate" ]; then
    candidate="$(dirname "$(readlink -f "$candidate")")"
    if is_orca_dir "$candidate"; then
      ORCA_FOUND="$candidate"
      return 0
    fi
  fi
  # Unpacked beside the other software or in HOME, newest version first.
  while IFS= read -r candidate; do
    if is_orca_dir "$candidate"; then
      ORCA_FOUND="$candidate"
      return 0
    fi
  done < <(ls -d "$PREFIX"/orca_*/ "$HOME"/orca_*/ 2>/dev/null | sed 's:/$::' | sort -V -r)
  while IFS= read -r candidate; do
    if extract_orca_tarball "$candidate"; then return 0; fi
  done < <(ls -1 "$HOME"/orca_*.tar.xz "$PREFIX"/orca_*.tar.xz 2>/dev/null | sort -V -r)
  return 0
}

# ORCA names the OpenMPI it is linked against: ..._openmpi418_... is 4.1.8.
required_openmpi() {
  local digits
  digits="$(basename "$1" | grep -o -E 'openmpi[0-9]{3,}' | head -n 1 | tr -dc '0-9')" || true
  [ -n "$digits" ] || return 0
  printf '%s.%s.%s\n' "${digits:0:1}" "${digits:1:1}" "${digits:2}"
}

mpirun_version() {
  { "$1" --version 2>/dev/null || true; } | grep -i 'open mpi' | grep -o -E '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true
}

build_openmpi() {
  local want="$1" target="$PREFIX/openmpi-$1" build="$PREFIX/.build" url
  url="${OPENMPI_URL:-https://download.open-mpi.org/release/open-mpi/v${want%.*}/openmpi-$want.tar.gz}"
  if [ "$DRY_RUN" = 1 ]; then
    log "would build OpenMPI $want into $target"
    OMPI_FOUND="$target"
    return 0
  fi
  if ! have gcc && ! have cc; then
    warn "no C compiler on PATH, so OpenMPI $want cannot be built; set OPENMPI_PREFIX to an existing $want build"
    return 1
  fi
  log "building OpenMPI $want into $target; this takes a while (log: $build/openmpi-$want.log)"
  mkdir -p "$build"
  if [ ! -f "$build/openmpi-$want.tar.gz" ]; then
    if ! fetch "$url" "$build/openmpi-$want.tar.gz.partial"; then
      rm -f "$build/openmpi-$want.tar.gz.partial"
      warn "OpenMPI $want could not be downloaded from $url"
      return 1
    fi
    mv "$build/openmpi-$want.tar.gz.partial" "$build/openmpi-$want.tar.gz"
  fi
  rm -rf "$build/openmpi-$want"
  tar -xzf "$build/openmpi-$want.tar.gz" -C "$build" || { warn "cannot unpack $build/openmpi-$want.tar.gz"; return 1; }
  if ! ( cd "$build/openmpi-$want" \
         && ./configure --prefix="$target" --enable-mpi-cxx --enable-mca-no-build=fs-gpfs --disable-oshmem \
         && make -j "$MAKE_JOBS" \
         && make install ) >"$build/openmpi-$want.log" 2>&1; then
    warn "the OpenMPI $want build failed; see $build/openmpi-$want.log"
    return 1
  fi
  rm -rf "$build/openmpi-$want"
  OMPI_FOUND="$target"
  log "OpenMPI $want installed at $target"
}

ensure_openmpi() {
  OMPI_FOUND=""
  if [ "$WANT_OPENMPI" != 1 ]; then
    log "OpenMPI: not set up (--no-openmpi)"
    return 0
  fi
  local want="" candidate version
  if [ -n "$ORCA_FOUND" ]; then want="$(required_openmpi "$ORCA_FOUND")"; fi
  if [ -z "$want" ]; then
    # Also without an ORCA to name it: an ORCA unpacked later then runs, and
    # Genarris has an mpicc to build against.
    want="$DEFAULT_OPENMPI"
  fi
  for candidate in "${OPENMPI_PREFIX:+$OPENMPI_PREFIX/bin/mpirun}" "$PREFIX/openmpi-$want/bin/mpirun" \
      "$(command -v mpirun 2>/dev/null || true)"; do
    if [ -z "$candidate" ] || [ ! -x "$candidate" ]; then continue; fi
    version="$(mpirun_version "$candidate")"
    if [ -n "$version" ] && [ "${version%.*}" = "${want%.*}" ]; then
      OMPI_FOUND="$(dirname "$(dirname "$(readlink -f "$candidate")")")"
      log "OpenMPI $version at $OMPI_FOUND (wanted: $want)"
      return 0
    fi
  done
  build_openmpi "$want"
}

# ---- DELFIN's own checkout, on --update ----------------------------------------
pull_checkout() {
  if [ ! -d "$DELFIN_REPO/.git" ]; then
    log "$DELFIN_REPO is not a git checkout, so it is not pulled"
    return 0
  fi
  if ! git -C "$DELFIN_REPO" rev-parse --abbrev-ref --symbolic-full-name '@{u}' >/dev/null 2>&1; then
    warn "$DELFIN_REPO has no upstream branch, so it is not pulled"
    return 0
  fi
  if [ -n "$(git -C "$DELFIN_REPO" status --porcelain --untracked-files=no)" ]; then
    warn "$DELFIN_REPO has local changes, so it is not pulled; commit or stash them and run --update again"
    return 0
  fi
  if [ "$DRY_RUN" = 1 ]; then
    log "would pull $DELFIN_REPO (fast-forward only)"
    return 0
  fi
  log "pulling $DELFIN_REPO (fast-forward only)"
  if ! git -C "$DELFIN_REPO" pull --ff-only; then
    warn "the pull did not fast-forward; the checkout was left as it was"
    FAILED+=("git pull")
  fi
}

# ---- Tools: one call into delfin/installer.py ---------------------------------------
run_tools() {
  local tools line group
  tools="$(planned_tools)"
  if [ "$MODE" = install ] && [ -z "$tools" ]; then
    return 0
  fi
  if [ "$DRY_RUN" = 1 ]; then
    case "$MODE" in
      install)
        while IFS= read -r line; do
          [ -n "$line" ] || continue
          group="${line%%:*}"
          if [ "$group" = ketcher ]; then
            log "would fetch Ketcher, the structure editor, into DELFIN's store"
          else
            log "would install ($group): ${line#*: }"
          fi
        done <<< "$PLAN" ;;
      update) log "would update ${tools:-every installed tool}" ;;
      repair) log "would check and repair ${tools:-every installed tool}" ;;
    esac
    return 0
  fi
  log "---- $MODE: ${tools:-every installed tool}"
  # shellcheck disable=SC2086
  if DELFIN_PYTHON="$VENV_PY" "$VENV_PY" -m delfin.installer "--$MODE" $tools; then
    DONE+=("$MODE: ${tools:-installed tools}")
  else
    FAILED+=("$MODE of tools (see the lines above)")
  fi
}

# ---- Settings and shell ---------------------------------------------------------------
write_settings_and_env() {
  if [ "$DRY_RUN" = 1 ]; then
    log "would record ORCA in ~/.delfin_settings.json and write ~/.delfin_env.sh"
    return 0
  fi
  local qm=0
  if printf '%s\n' "$PLAN" | grep -q '^qm:'; then qm=1; fi
  DELFIN_INSTALL_ORCA="$ORCA_FOUND" DELFIN_INSTALL_ORCA_EXPLICIT="$ORCA_EXPLICIT" \
  DELFIN_INSTALL_OMPI="$OMPI_FOUND" DELFIN_INSTALL_REPO="$DELFIN_REPO" DELFIN_INSTALL_QM="$qm" \
  "$VENV_PY" - <<'PY'
import os
from pathlib import Path

from delfin.qm_runtime import get_user_qm_tools_root
from delfin.runtime_setup import ensure_shell_sources_delfin_env, write_delfin_env_file
from delfin.user_settings import load_settings, save_settings

orca = os.environ.get("DELFIN_INSTALL_ORCA", "")
settings = load_settings()
runtime = settings.get("runtime") or {}
local = runtime.get("local") or {}
slurm = runtime.get("slurm") or {}

# A choice already in the settings is the user's, unless ORCA was named just now.
for section in (runtime, local, slurm):
    if orca and (os.environ.get("DELFIN_INSTALL_ORCA_EXPLICIT") == "1" or not section.get("orca_base")):
        section["orca_base"] = orca
if os.environ.get("DELFIN_INSTALL_QM") == "1" and not runtime.get("qm_tools_root"):
    runtime["qm_tools_root"] = str(get_user_qm_tools_root())
runtime["local"], runtime["slurm"] = local, slurm
settings["runtime"] = runtime
save_settings(settings)
print(f"[delfin-install] settings: {Path.home() / '.delfin_settings.json'}")

env_file = write_delfin_env_file(
    repo_dir=os.environ["DELFIN_INSTALL_REPO"],
    orca_base=runtime.get("orca_base") or "",
    qm_tools_root=runtime.get("qm_tools_root") or "",
    openmpi_prefix=os.environ.get("DELFIN_INSTALL_OMPI", ""),
)
for rc_file in ensure_shell_sources_delfin_env(env_path=env_file):
    print(f"[delfin-install] {rc_file} sources {env_file}")
PY
}

summary() {
  local item
  log "==================== summary ===================="
  if [ -n "$VENV_PY" ] && [ -x "$VENV_PY" ]; then
    log "DELFIN   $("$VENV_PY" -c 'import delfin; print(getattr(delfin, "__version__", "?"))' 2>/dev/null || echo '?') in $VENV_DIR"
  fi
  log "ORCA     ${ORCA_FOUND:-not found; pass --orca DIR|TARBALL (https://orcaforum.kofo.mpg.de)}"
  log "OpenMPI  ${OMPI_FOUND:-not set up}"
  for item in "${DONE[@]+"${DONE[@]}"}"; do log "ok       $item"; done
  for item in "${FAILED[@]+"${FAILED[@]}"}"; do log "FAILED   $item"; done
  if [ "$DRY_RUN" = 1 ]; then
    log "dry run: nothing was changed"
    return 0
  fi
  log "next: open a new shell (or source ~/.delfin_env.sh), then run delfin-voila"
}

main() {
  if [ "$LIST_ONLY" = 1 ]; then
    catalog --list
    return 0
  fi
  [ "$(uname -s)" = Linux ] || die "this installer supports Linux"
  [ "$(uname -m)" = x86_64 ] || warn "the prebuilt tools are x86_64 builds; on $(uname -m) some will not install"
  resolve_repo
  select_tools
  local scheduler="none" scope
  if have sbatch; then scheduler="SLURM"; fi
  case "$MODE" in
    install) scope="profile $PROFILE" ;;
    *) scope="$MODE" ;;
  esac
  if [ -n "$ONLY" ]; then scope="$scope, only $ONLY"; fi
  log "plan: $scope | prefix $PREFIX | checkout $DELFIN_REPO | scheduler $scheduler"
  confirm_all
  ensure_repo
  if [ "$MODE" = update ]; then pull_checkout; fi
  ensure_venv
  if [ "$DRY_RUN" != 1 ]; then mkdir -p "$DELFIN_CALC_DIR" "$DELFIN_ARCHIVE_DIR"; fi
  find_orca
  if ! ensure_openmpi; then FAILED+=("openmpi"); fi
  run_tools
  write_settings_and_env
  summary
  [ ${#FAILED[@]} -eq 0 ]
}

main
