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
#   --extras LIST     pip extras for DELFIN (default: agent,docs; add pdf for
#                     the PDF viewer, PyMuPDF, AGPL-3.0)
#   --dry-run         print what would be done and change nothing
#   --list            print every tool that can be installed, by group
#   --update          update instead of install (--only limits it)
#   --repair          repair instead of install (--only limits it)
#   -y, --yes         do not ask before the "all" profile
#
# ORCA cannot be downloaded for you (license). Pass --orca with an installed
# ORCA directory or the tarball from https://orcaforum.kofo.mpg.de, or leave
# either in $HOME or the prefix directory, where it is found.
#
# Before anything is changed, what the run needs is checked (curl or wget,
# tar, git when DELFIN has to be fetched, a writable prefix and HOME, free
# disk space, a compiler where something is built). A real run keeps its
# whole output in PREFIX/logs/; a step that fails names itself, and running
# the same command again continues from what is already in place.
# ============================================================================
# Started with sh (dash on Debian and Ubuntu) or zsh, the bash below stops at
# its first line ("Illegal option -o pipefail"); bash takes over instead.
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

log()  { printf '[delfin-install] %s\n' "$*"; }
warn() { printf '[delfin-install] WARNING: %s\n' "$*" >&2; }
die()  { printf '[delfin-install] ERROR: %s\n' "$*" >&2; EXPECTED_EXIT=1; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

# What the run is doing, for the message when it stops: a failed command under
# set -e ends the script with no word of where it was.
STEP="reading the options"
LOG_FILE=""
EXPECTED_EXIT=0
TEE_PIDS=()
# The tee processes of the log (start_log) are waited for, so nothing they
# still hold is printed after the shell prompt has come back.
end_log() {
  [ ${#TEE_PIDS[@]} -gt 0 ] || return 0
  exec >&- 2>&-
  wait "${TEE_PIDS[@]}" 2>/dev/null || true
}
on_exit() {
  local status="$1"
  if [ "$status" -ne 0 ] && [ "$EXPECTED_EXIT" != 1 ]; then
    if [ "$status" -ge 128 ]; then
      printf '[delfin-install] ERROR: stopped (signal %s) while %s.\n' "$((status - 128))" "$STEP" >&2
    else
      printf '[delfin-install] ERROR: %s failed (exit %s); the steps after it were not run.\n' "$STEP" "$status" >&2
    fi
    printf '[delfin-install] What finished is kept. Fix the cause shown above, then run the same command again: it continues from there.\n' >&2
  fi
  if [ "$status" -ne 0 ] && [ -n "$LOG_FILE" ]; then
    printf '[delfin-install] full log: %s\n' "$LOG_FILE" >&2
  fi
  end_log
}
trap 'on_exit $?' EXIT

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
# `test` is in the default set because the agent runs the suite as part
# of ordinary work: pytest lived only in `dev`, which this default does
# not select, so every site install produced an environment whose test
# tool could not start. Two packages.
DELFIN_EXTRAS="${DELFIN_EXTRAS-agent,docs,test}"
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

# Python that imports delfin runs from this script's own directory. Started
# from anywhere else, `python -m` puts the working directory first on the
# import path, and inside another DELFIN checkout that checkout's delfin is
# the one imported -- its tool directory staged, its links copied.
in_script_dir() {
  (cd "$SCRIPT_DIR" && "$@")
}

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
  PYTHONPATH="$PACKAGE_PARENT${PYTHONPATH:+:$PYTHONPATH}" in_script_dir "$py" -m delfin.installer "$@"
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

# ---- Before anything changes: what this run needs ------------------------------
# The nearest existing directory at or above $1.
existing_parent() {
  local dir="$1"
  while [ ! -e "$dir" ] && [ "$dir" != / ]; do dir="$(dirname "$dir")"; done
  printf '%s\n' "$dir"
}

# Whole GB free on the file system holding $1; empty when df cannot say. A
# quota (usual for HOME on a cluster) is not seen by df.
free_gb() {
  { df -Pk "$1" 2>/dev/null || true; } | awk 'NR == 2 && $4 ~ /^[0-9]+$/ { printf "%d\n", $4 / 1048576 }'
}

# A path as a person reads it: ~/x for one under HOME (printed, not expanded).
shown_path() {
  # shellcheck disable=SC2088
  case "$1" in "$HOME"/*) printf '~/%s\n' "${1#"$HOME"/}" ;; *) printf '%s\n' "$1" ;; esac
}

reachable() {
  if have curl; then
    curl -fsSI --max-time 20 -o /dev/null "$1" >/dev/null 2>&1
  else
    wget -q --spider --timeout=20 "$1" >/dev/null 2>&1
  fi
}

preflight() {
  local missing=() dir parent free need where tool absent=()
  if ! have curl && ! have wget; then
    missing+=("curl or wget (downloads micromamba, OpenMPI, Ketcher and the QM tools)")
  fi
  have tar || missing+=("tar (unpacks OpenMPI, ORCA and the tools)")
  if ! grep -qs 'name = "delfin-complat"' "$DELFIN_REPO/pyproject.toml" \
      || { [ "$MODE" = update ] && [ -d "$DELFIN_REPO/.git" ]; }; then
    have git || missing+=("git (fetches and updates DELFIN)")
  fi
  if [ ${#missing[@]} -gt 0 ]; then
    for tool in "${missing[@]}"; do warn "missing: $tool"; done
    die "nothing was changed. Install what is missing (on a shared machine: ask its administrators, or load it as a module), then run the same command again."
  fi

  for dir in "$PREFIX" "$HOME"; do
    parent="$(existing_parent "$dir")"
    if [ ! -d "$parent" ] || [ ! -w "$parent" ]; then
      if [ "$dir" = "$PREFIX" ]; then
        die "$PREFIX cannot be created: $parent is not a writable directory. Nothing was changed; choose another place with --prefix DIR."
      fi
      die "HOME ($HOME) is not writable; DELFIN keeps its settings and tools there. Nothing was changed."
    fi
  done

  # Measured for a fresh install (2026-11-13): core 2.5 GB, standard 4.3 GB
  # (venv 1.6, micromamba packages 1.1, ~/.delfin 0.8, pip cache 0.4).
  case "$MODE:$PROFILE" in
    install:core) need=3 ;;
    install:all) need=20 ;;
    install:*) need=6 ;;
    *) need=2 ;;
  esac
  for dir in "$PREFIX" "$HOME"; do
    where="$(existing_parent "$dir")"
    free="$(free_gb "$where")"
    if [ -n "$free" ] && [ "$free" -lt "$need" ]; then
      warn "only $free GB free on the file system of $(shown_path "$dir"); this run needs about $need GB (the venv in the checkout, micromamba packages under the prefix, tools in ~/.delfin, the pip cache in ~/.cache). On a cluster, put the checkout and --prefix into a workspace."
    fi
  done

  if [ "$WANT_OPENMPI" = 1 ] && { ! { have gcc || have cc; } || ! have make; }; then
    warn "no C compiler (gcc or cc) or no make on PATH: if OpenMPI has to be built, that step will fail. Load a compiler first (on a cluster usually a module), or pass --no-openmpi and set OPENMPI_PREFIX to an existing build later."
  fi
  # std2 has no binary release: one on PATH or from an earlier run is linked,
  # otherwise it is built.
  if printf '%s\n' "$PLAN" | grep -q -E '^qm:.*(^| )std2( |$)' \
      && ! have std2 && [ ! -x "$HOME/.delfin/qm_tools/bin/std2" ]; then
    for tool in meson ninja gfortran gcc; do have "$tool" || absent+=("$tool"); done
    if [ ${#absent[@]} -gt 0 ]; then
      warn "std2 is built from source and needs meson, ninja, gfortran and gcc; missing: ${absent[*]}. Everything else installs; std2 will be reported as failed."
    fi
  fi

  if [ "$DRY_RUN" != 1 ]; then
    where="${PIP_INDEX_URL:-https://pypi.org/simple/}"
    if ! reachable "$where"; then
      warn "$where cannot be reached. Without a connection nothing can be installed; if this machine goes out through a proxy, export https_proxy=http://HOST:PORT (and http_proxy) first."
    fi
  fi
}

# A real run keeps everything it prints, so a failure can be read afterwards.
start_log() {
  [ "$DRY_RUN" != 1 ] || return 0
  mkdir -p "$PREFIX/logs"
  LOG_FILE="$PREFIX/logs/delfin-install-$(date +%Y%m%d-%H%M%S).log"
  exec > >(tee -a "$LOG_FILE")
  TEE_PIDS+=("$!")
  exec 2> >(tee -a "$LOG_FILE" >&2)
  TEE_PIDS+=("$!")
}

# Lines in the login files that activate another venv. They run before the
# line that loads ~/.delfin_env.sh, which then keeps that venv instead of
# activating DELFIN's.
other_activations() {
  local rc line
  for rc in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.profile"; do
    [ -f "$rc" ] || continue
    while IFS= read -r line; do
      case "$line" in *"$VENV_DIR/bin/activate"*) continue ;; esac
      printf '%s:%s\n' "$(shown_path "$rc")" "$line"
    done < <(grep -n -E '^[[:space:]]*(source|\.)[[:space:]]+[^#]*bin/activate' "$rc" 2>/dev/null || true)
  done
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
  local env_dir="$PREFIX/delfin-python" mamba fetcher
  if [ -x "$env_dir/bin/python" ] && python_ok "$env_dir/bin/python"; then
    BASE_PYTHON="$env_dir/bin/python"
    return 0
  fi
  if ! mamba="$(find_micromamba)"; then
    # The shared fetcher: one pinned release per platform, held to its
    # SHA-256. Beside this script, else in the checkout resolve_repo found
    # or cloned (this script may have been downloaded on its own).
    fetcher="$SCRIPT_DIR/fetch_micromamba.sh"
    [ -f "$fetcher" ] || fetcher="$DELFIN_REPO/delfin/installers/fetch_micromamba.sh"
    [ -f "$fetcher" ] || die "no Python 3.10/3.11 here and no fetch_micromamba.sh to fetch one with; pass --python PATH"
    log "no Python 3.10/3.11 here; fetching micromamba"
    bash "$fetcher" "$PREFIX/micromamba/bin/micromamba" \
      || die "micromamba could not be fetched and verified; pass --python PATH"
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
    STEP="creating the Python for DELFIN's venv"
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
  STEP="installing DELFIN into $VENV_DIR (pip)"
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

# delfin.qm_runtime.find_orca_dir, asked of DELFIN's venv (or any python3 that
# can import this checkout's delfin). Status 1 when no Python could answer.
delfin_orca_dir() {
  local py="$VENV_PY" answer
  if [ -z "$py" ] || [ ! -x "$py" ]; then py="$(command -v python3 2>/dev/null || true)"; fi
  [ -n "$py" ] || return 1
  answer="$(PYTHONPATH="$PACKAGE_PARENT${PYTHONPATH:+:$PYTHONPATH}" in_script_dir "$py" -c \
    'import sys; from delfin.qm_runtime import find_orca_dir; print(find_orca_dir(sys.argv[1:]))' \
    "$PREFIX" 2>/dev/null)" || return 1
  printf '%s\n' "$answer" | tail -n 1
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
  # DELFIN's own answer, the one `delfin doctor` reports: settings, the tool
  # resolver (ORCA_BINARY/ORCA_PATH, PATH, system directories, modules) and
  # orca_* directories unpacked in the prefix, ~/software or ~.
  if candidate="$(delfin_orca_dir)"; then
    if [ -n "$candidate" ] && is_orca_dir "$candidate"; then
      ORCA_FOUND="$candidate"
      return 0
    fi
  else
    # Only where no Python can import DELFIN yet (a dry run on a machine
    # without one): the unpacked directories, newest version first.
    while IFS= read -r candidate; do
      if is_orca_dir "$candidate"; then
        ORCA_FOUND="$candidate"
        return 0
      fi
    done < <(ls -d "$PREFIX"/orca_*/ "$HOME"/orca_*/ 2>/dev/null | sed 's:/$::' | sort -V -r)
  fi
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

sha256_of() {
  # The digest alone on stdout; status 1 when no tool can compute one.
  if have sha256sum; then
    sha256sum "$1" | awk '{print $1}'
  elif have shasum; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif have python3; then
    python3 -I -c 'import hashlib, sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "$1"
  else
    return 1
  fi
}

# The SHA-256 open-mpi.org publishes for openmpi-<version>.tar.gz, read
# 2026-10-16 from www.open-mpi.org/software/ompi/v4.1/ and /v5.0/; the
# tarballs of 4.1.1, 4.1.6, 4.1.8 and 5.0.8 were downloaded and agree.
# A release not listed needs OPENMPI_SHA256, the digest from that page.
openmpi_sha256() {
  case "$1" in
    4.1.0) printf '%s\n' 228467c3dd15339d9b26cf26a291af3ee7c770699c5e8a1b3ad786f9ae78140a ;;
    4.1.1) printf '%s\n' d80b9219e80ea1f8bcfe5ad921bd9014285c4948c5965f4156a3831e60776444 ;;
    4.1.2) printf '%s\n' a400719b04375cd704d2ed063a50e42d268497a3dfede342986ab7a8d7e8dcf0 ;;
    4.1.3) printf '%s\n' 9c0fd1f78fc90ca9b69ae4ab704687d5544220005ccd7678bf58cc13135e67e0 ;;
    4.1.4) printf '%s\n' e166dbe876e13a50c2882e11193fecbc4362e89e6e7b6deeb69bf095c0f4fc4c ;;
    4.1.5) printf '%s\n' c018b127619d2a2a30c1931f316fc8a245926d0f5b4ebed4711f9695e7f70925 ;;
    4.1.6) printf '%s\n' 44da277b8cdc234e71c62473305a09d63f4dcca292ca40335aab7c4bf0e6a566 ;;
    4.1.7) printf '%s\n' 3956daf86e8b570e31155d92ed4b2b1699bc36455ae24a0d98ca376d3e8b4e9d ;;
    4.1.8) printf '%s\n' fb41086bbed9300baa2f3d7572491facfe5257412fa524ec5a396aa9101d5c62 ;;
    5.0.0) printf '%s\n' 4bf81fc86f562b372c40d44a09372ba1fa9b780d50d09a46ed0c4c7f09250b71 ;;
    5.0.1) printf '%s\n' 15805510d599558aed2ef43770e8a4683b9a6b361d0a91f107cb7c377cfe2bfb ;;
    5.0.2) printf '%s\n' 095ab1cddb0fa0f9e7fc211a1d33185c6727c5237d0ee55f80a7e4311e5d279c ;;
    5.0.3) printf '%s\n' 3e12160dd1d6f5cd8897bbab63a05fb9223a1e20a1726db23607aac25dc1aa19 ;;
    5.0.4) printf '%s\n' 182763628c6cfbd732f83125cc1acf001eb445af45e5a406e65923f90a69f697 ;;
    5.0.5) printf '%s\n' 5cbefa0780b84f4126743c40cdd6a334b2f0574cd7fd95050fb1ac0ddbb7f0b8 ;;
    5.0.6) printf '%s\n' 1d6dd41f6b53c00db23b40e56733d16a996fa743957ef8add8e117feffd92689 ;;
    5.0.7) printf '%s\n' 67435fdb2560f897882e69fead171d38e7f44a6dfe1d3c90506f549a943050b8 ;;
    5.0.8) printf '%s\n' f891ddf2dab3b604f2521d0569615bc1c07a1ac86a295ac543e059aecb303621 ;;
    5.0.9) printf '%s\n' be35ac3aae68609dedf56a3c6480cb1b0af9ae40b98cc83ef27bc859084a9aa8 ;;
    5.0.10) printf '%s\n' 5692cc80554a7117c99eaa725d35100edd8bbf73423a5e265ff867979192df7d ;;
    5.0.11) printf '%s\n' 77a402c64c22f1bd8a504f76fd411df23e36e168eedce6c9bcd12c8c0cd0619a ;;
    *) return 1 ;;
  esac
}

# A prefix whose ompi_info starts: its libraries load. mpirun --version alone
# answers even when a library it needs at run time is gone.
openmpi_works() {
  [ -x "$1/bin/ompi_info" ] && "$1/bin/ompi_info" --parsable >/dev/null 2>&1
}

mpirun_version() {
  { "$1" --version 2>/dev/null || true; } | grep -i 'open mpi' | grep -o -E '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true
}

build_openmpi() {
  local want="$1" target="$PREFIX/openmpi-$1" build="$PREFIX/.build" url digest got
  url="${OPENMPI_URL:-https://download.open-mpi.org/release/open-mpi/v${want%.*}/openmpi-$want.tar.gz}"
  if [ "$DRY_RUN" = 1 ]; then
    log "would build OpenMPI $want into $target"
    OMPI_FOUND="$target"
    return 0
  fi
  digest="${OPENMPI_SHA256:-}"
  if [ -z "$digest" ] && ! digest="$(openmpi_sha256 "$want")"; then
    warn "no SHA-256 is pinned for OpenMPI $want, so its source is not downloaded; set OPENMPI_SHA256 to the digest www.open-mpi.org publishes for openmpi-$want.tar.gz, or OPENMPI_PREFIX to an existing $want build"
    return 1
  fi
  if { ! have gcc && ! have cc; } || ! have make; then
    warn "no C compiler or no make on PATH, so OpenMPI $want cannot be built; load a compiler and run again, or set OPENMPI_PREFIX to an existing $want build"
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
  # Checked every time, a tarball kept from an earlier run included.
  if ! got="$(sha256_of "$build/openmpi-$want.tar.gz")"; then
    rm -f "$build/openmpi-$want.tar.gz"
    warn "no sha256sum, shasum or python3 to check the OpenMPI $want source; it is not built"
    return 1
  fi
  if [ "$got" != "$digest" ]; then
    rm -f "$build/openmpi-$want.tar.gz"
    warn "checksum mismatch for openmpi-$want.tar.gz from $url: expected $digest, got $got; the file was removed and nothing was built"
    return 1
  fi
  log "checksum verified: openmpi-$want.tar.gz"
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
  local want="" candidate version prefix
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
      prefix="$(dirname "$(dirname "$(readlink -f "$candidate")")")"
      if openmpi_works "$prefix"; then
        OMPI_FOUND="$prefix"
        log "OpenMPI $version at $OMPI_FOUND (wanted: $want)"
        return 0
      fi
      warn "OpenMPI $version at $prefix does not start (ompi_info fails); a working one is built"
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
  if DELFIN_PYTHON="$VENV_PY" in_script_dir "$VENV_PY" -m delfin.installer "--$MODE" $tools; then
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
  in_script_dir "$VENV_PY" - <<'PY'
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
  local item line others
  log "==================== summary ===================="
  if [ -n "$VENV_PY" ] && [ -x "$VENV_PY" ]; then
    log "DELFIN   $(in_script_dir "$VENV_PY" -c 'import delfin; print(getattr(delfin, "__version__", "?"))' 2>/dev/null || echo '?') in $VENV_DIR"
  fi
  log "ORCA     ${ORCA_FOUND:-not found: download it from https://orcaforum.kofo.mpg.de (free for academic use), then run this again with --orca TARBALL}"
  log "OpenMPI  ${OMPI_FOUND:-not set up}"
  for item in "${DONE[@]+"${DONE[@]}"}"; do log "ok       $item"; done
  for item in "${FAILED[@]+"${FAILED[@]}"}"; do log "FAILED   $item"; done
  if [ "$DRY_RUN" = 1 ]; then
    log "dry run: nothing was changed"
    return 0
  fi
  log "settings ~/.delfin_settings.json; environment ~/.delfin_env.sh, loaded by ~/.bashrc and ~/.profile"
  log "log      $LOG_FILE"
  others="$(other_activations)"
  if [ -n "$others" ]; then
    while IFS= read -r line; do
      warn "$line activates another Python environment, so new shells keep that one instead of DELFIN's. Remove the line, or move it below the line that sources ~/.delfin_env.sh."
    done <<< "$others"
  fi
  if [ ${#FAILED[@]} -gt 0 ]; then
    warn "what is marked FAILED above did not install; DELFIN itself is installed. Fix the cause in the lines above and run the same command again (or: bash $DELFIN_REPO/install.sh --repair)."
  fi
  cat <<EOF

Start DELFIN (a new shell loads its environment by itself):
  source ~/.delfin_env.sh        load it into this shell now
  delfin-voila                   the dashboard: open the URL it prints
  cd $(shown_path "$DELFIN_CALC_DIR") && delfin-agent      the agent in a terminal; it works in the folder it is started in
  delfin doctor                  check ORCA, xTB, OpenMPI and the rest
  (in a VS Code terminal also: source $VENV_DIR/bin/activate)

The agent needs a language model: choose one in the dashboard under Settings,
or store an API key with  delfin-agent credentials set KIT_TOOLBOX_API_KEY
(or OPENAI_API_KEY; delfin-agent credentials list names the others).

On a remote machine such as a cluster login node the dashboard listens on
127.0.0.1 only. From your own computer:
  ssh -L 8866:127.0.0.1:8866 ${USER:-$(id -un)}@$(hostname 2>/dev/null || echo HOST)
then open the URL delfin-voila printed.

Update DELFIN and its tools:  bash $DELFIN_REPO/install.sh --update
Add a tool:                   bash $DELFIN_REPO/install.sh --only crest   (what exists: --list)
EOF
}

main() {
  if [ "$LIST_ONLY" = 1 ]; then
    catalog --list
    return 0
  fi
  STEP="checking this machine"
  [ "$(uname -s)" = Linux ] || die "this installer supports Linux"
  [ "$(uname -m)" = x86_64 ] || warn "the prebuilt tools are x86_64 builds; on $(uname -m) some will not install"
  resolve_repo
  STEP="reading DELFIN's list of tools"
  select_tools
  local scheduler="none" scope
  if have sbatch; then scheduler="SLURM"; fi
  case "$MODE" in
    install) scope="profile $PROFILE" ;;
    *) scope="$MODE" ;;
  esac
  if [ -n "$ONLY" ]; then scope="$scope, only $ONLY"; fi
  log "plan: $scope | prefix $PREFIX | checkout $DELFIN_REPO | scheduler $scheduler"
  STEP="checking what this run needs"
  preflight
  confirm_all
  start_log
  STEP="fetching DELFIN"
  ensure_repo
  STEP="updating the checkout"
  if [ "$MODE" = update ]; then pull_checkout; fi
  STEP="setting up DELFIN's venv"
  ensure_venv
  if [ "$DRY_RUN" != 1 ]; then mkdir -p "$DELFIN_CALC_DIR" "$DELFIN_ARCHIVE_DIR"; fi
  STEP="looking for ORCA"
  find_orca
  STEP="setting up OpenMPI"
  if ! ensure_openmpi; then FAILED+=("openmpi"); fi
  STEP="installing the tools"
  run_tools
  STEP="writing the settings and ~/.delfin_env.sh"
  write_settings_and_env
  STEP="printing the summary"
  summary
  if [ ${#FAILED[@]} -gt 0 ]; then
    EXPECTED_EXIT=1
    exit 1
  fi
}

main
