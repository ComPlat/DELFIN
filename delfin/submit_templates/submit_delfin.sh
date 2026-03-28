#!/usr/bin/env bash
set -euo pipefail

# ========================================================================
# DELFIN SLURM Runner — with automatic stage-in/stage-out
#
# When $TMPDIR is available (SLURM node-local SSD), all job I/O is
# performed on the fast local disk instead of the shared HOME filesystem.
# Results are synced back periodically and on exit.
#
# Set DELFIN_STAGE_IO=0 to disable staging and run directly in HOME.
# Set DELFIN_SYNC_INTERVAL to control periodic sync (default: 900 = 15min).
# ========================================================================

ORIGIN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
STAGE_IO="${DELFIN_STAGE_IO:-1}"
SYNC_INTERVAL="${DELFIN_SYNC_INTERVAL:-900}"

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
  # Stage-in may fail if $TMPDIR is too small for the job data.
  # Fall back to running directly on HOME in that case.
  _STAGE_OK=0
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --copy-links "$ORIGIN_DIR/" "$WORK_DIR/" && _STAGE_OK=1
  else
    cp -rL "$ORIGIN_DIR/." "$WORK_DIR/" && _STAGE_OK=1
  fi
  if [ "$_STAGE_OK" = "1" ]; then
    export DELFIN_SCRATCH="$TMPDIR"
    echo "[stage-in] Done. Working in $WORK_DIR (local SSD)"
  else
    echo "[stage-in] WARNING: copy failed (disk full?). Running directly on HOME."
    rm -rf "$WORK_DIR" 2>/dev/null || true
    WORK_DIR="$ORIGIN_DIR"
  fi
fi

cd "$WORK_DIR"

# ------------------------------------------------------------------
# Periodic background sync: incrementally copy results back every
# SYNC_INTERVAL seconds so intermediate results survive crashes.
# Excludes temporary files to minimize HOME I/O.
# ------------------------------------------------------------------
_SYNC_PID=""
RSYNC_EXCLUDES="--exclude=.orca_iso_* --exclude=*.tmp --exclude=__pycache__ --exclude=.exit_code_*"

_start_periodic_sync() {
  if [ "$WORK_DIR" = "$ORIGIN_DIR" ]; then return; fi
  if [ "$SYNC_INTERVAL" -le 0 ] 2>/dev/null; then return; fi
  (
    while true; do
      sleep "$SYNC_INTERVAL"
      rsync -a --update $RSYNC_EXCLUDES "$WORK_DIR/" "$ORIGIN_DIR/" 2>/dev/null || true
      echo "[sync] Periodic sync completed at $(date +%H:%M:%S)"
    done
  ) &
  _SYNC_PID=$!
  echo "[sync] Periodic background sync every ${SYNC_INTERVAL}s (PID $_SYNC_PID)"
}

_stop_periodic_sync() {
  if [ -n "$_SYNC_PID" ]; then
    kill "$_SYNC_PID" 2>/dev/null || true
    wait "$_SYNC_PID" 2>/dev/null || true
    _SYNC_PID=""
  fi
}

# ------------------------------------------------------------------
# Stage-out: copy results back to HOME on any exit.
# Handles normal exit, errors, SIGTERM (SLURM timeout), SIGUSR1
# (early warning before walltime), etc.
# ------------------------------------------------------------------
_DELFIN_PID=""

_stage_out() {
  set +e  # must not abort inside trap
  local rc=${1:-$?}
  _stop_periodic_sync
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
  if [ -n "$_DELFIN_PID" ]; then
    kill -TERM -- -"$_DELFIN_PID" 2>/dev/null || kill -TERM "$_DELFIN_PID" 2>/dev/null || true
    wait "$_DELFIN_PID" 2>/dev/null || true
  fi
  _stage_out 124  # 124 = timeout convention
}

# SIGUSR1: early warning before walltime (use --signal=B:USR1@300 in sbatch).
# Syncs results immediately while DELFIN keeps running — it will be
# SIGTERMed later when the actual walltime hits.
_handle_early_warning() {
  echo ""
  echo "[warning] Approaching walltime -- performing preemptive sync..."
  if [ "$WORK_DIR" != "$ORIGIN_DIR" ]; then
    rsync -a --update $RSYNC_EXCLUDES "$WORK_DIR/" "$ORIGIN_DIR/" 2>/dev/null || true
    echo "[warning] Preemptive sync done."
  fi
}

trap '_stage_out $?' EXIT
trap '_handle_signal' SIGTERM SIGINT SIGHUP
trap '_handle_early_warning' SIGUSR1

# ------------------------------------------------------------------
# Run DELFIN
# ------------------------------------------------------------------
if [ "$WORK_DIR" = "$ORIGIN_DIR" ]; then
  # No staging — exec directly (no trap needed)
  exec "$PYTHON_BIN" -m delfin.dashboard.local_runner
else
  _start_periodic_sync
  # Run in background so we can catch signals for stage-out
  "$PYTHON_BIN" -m delfin.dashboard.local_runner &
  _DELFIN_PID=$!
  wait "$_DELFIN_PID"
fi
