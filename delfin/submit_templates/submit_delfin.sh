#!/usr/bin/env bash
set -euo pipefail

JOB_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$JOB_DIR"

PYTHON_BIN="${DELFIN_PYTHON_BIN:-}"
if [ -z "$PYTHON_BIN" ]; then
  if command -v python >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python)"
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3)"
  else
    echo "ERROR: python/python3 not found in PATH for DELFIN SLURM runner." >&2
    exit 1
  fi
fi

export DELFIN_JOB_ID="${DELFIN_JOB_ID:-${SLURM_JOB_ID:-0}}"
export DELFIN_JOB_NAME="${DELFIN_JOB_NAME:-${SLURM_JOB_NAME:-slurm_job}}"

exec "$PYTHON_BIN" -m delfin.dashboard.local_runner
