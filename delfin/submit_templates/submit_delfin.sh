#!/usr/bin/env bash
set -euo pipefail

# ========================================================================
# DELFIN SLURM Runner — with automatic stage-in/stage-out
#
# When $TMPDIR is available (SLURM node-local SSD), all job I/O is
# performed on the fast local disk instead of the shared HOME filesystem.
# Results are copied back to the original job directory on exit.
#
# Set DELFIN_STAGE_IO=0 to disable staging and run directly in HOME.
# ========================================================================

ORIGIN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
STAGE_IO="${DELFIN_STAGE_IO:-1}"

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

# ------------------------------------------------------------------
# Stage-in: copy job directory to node-local SSD
# ------------------------------------------------------------------
WORK_DIR="$ORIGIN_DIR"

if [ "$STAGE_IO" = "1" ] && [ -n "${TMPDIR:-}" ] && [ -d "${TMPDIR}" ]; then
  WORK_DIR="$TMPDIR/delfin_job"
  echo "[stage-in] Copying $ORIGIN_DIR -> $WORK_DIR"
  mkdir -p "$WORK_DIR"
  # Use rsync to handle symlinks gracefully; fall back to cp.
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --copy-links "$ORIGIN_DIR/" "$WORK_DIR/"
  else
    cp -rL "$ORIGIN_DIR/." "$WORK_DIR/"
  fi
  export DELFIN_SCRATCH="$TMPDIR"
  echo "[stage-in] Done. Working in $WORK_DIR (local SSD)"
fi

cd "$WORK_DIR"

# ------------------------------------------------------------------
# Stage-out helper: copy results back to HOME on any exit
# This runs on normal exit, errors, SIGTERM (SLURM timeout), etc.
# ------------------------------------------------------------------
_DELFIN_PID=""

_stage_out() {
  set +e  # must not abort inside trap
  local rc=${1:-$?}
  if [ "$WORK_DIR" != "$ORIGIN_DIR" ]; then
    echo ""
    echo "[stage-out] Copying results $WORK_DIR -> $ORIGIN_DIR"
    if command -v rsync >/dev/null 2>&1; then
      rsync -a --update "$WORK_DIR/" "$ORIGIN_DIR/" 2>&1 || echo "[stage-out] WARNING: rsync had errors"
    else
      cp -ru "$WORK_DIR/." "$ORIGIN_DIR/" 2>&1 || echo "[stage-out] WARNING: cp had errors"
    fi
    echo "[stage-out] Done."
  fi
  exit "$rc"
}

_handle_signal() {
  echo ""
  echo "[signal] Received termination signal — stopping DELFIN and saving results..."
  # Kill the DELFIN process group so ORCA children also stop
  if [ -n "$_DELFIN_PID" ]; then
    kill -TERM -- -"$_DELFIN_PID" 2>/dev/null || kill -TERM "$_DELFIN_PID" 2>/dev/null || true
    wait "$_DELFIN_PID" 2>/dev/null || true
  fi
  _stage_out 124  # 124 = timeout convention
}

trap '_stage_out $?' EXIT
trap '_handle_signal' SIGTERM SIGINT SIGHUP

# ------------------------------------------------------------------
# Run DELFIN
# ------------------------------------------------------------------
if [ "$WORK_DIR" = "$ORIGIN_DIR" ]; then
  # No staging — exec directly (no trap needed)
  exec "$PYTHON_BIN" -m delfin.dashboard.local_runner
else
  # Staging active — run in background so we can catch signals for stage-out
  "$PYTHON_BIN" -m delfin.dashboard.local_runner &
  _DELFIN_PID=$!
  wait "$_DELFIN_PID"
fi
