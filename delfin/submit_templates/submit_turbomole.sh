#!/usr/bin/env bash
set -euo pipefail

JOB_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
RUN_DIR="${TMPDIR:-$JOB_DIR}/turbomole_${SLURM_JOB_ID:-manual}"
TM_COMMAND="${TM_MODULE:-ridft}"

mkdir -p "$RUN_DIR"

cleanup() {
  rsync -a "$RUN_DIR"/ "$JOB_DIR"/ 2>/dev/null || true
  rm -rf "$RUN_DIR" 2>/dev/null || true
}

trap cleanup EXIT SIGTERM SIGINT SIGHUP

cp -a "$JOB_DIR"/. "$RUN_DIR"/ 2>/dev/null || true
cd "$RUN_DIR"

if ! command -v "$TM_COMMAND" >/dev/null 2>&1; then
  echo "ERROR: TURBOMOLE command '$TM_COMMAND' not found in PATH." >&2
  echo "Configure a site-specific submit template or load the required module before starting DELFIN." >&2
  exit 1
fi

exec "$TM_COMMAND"
