#!/usr/bin/env bash
#SBATCH --job-name=pre_DFT_SP
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sp-%A_%a.out
#SBATCH --error=logs/sp-%A_%a.err
set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Missing SLURM_SUBMIT_DIR}"
: "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
ROOT="$SLURM_SUBMIT_DIR"
mapfile -t JOBS < "$ROOT/jobs.sp.submission.txt"
idx="$SLURM_ARRAY_TASK_ID"
(( idx >= 0 && idx < ${#JOBS[@]} )) || { echo "Invalid SP index $idx" >&2; exit 2; }
JOBDIR="${JOBS[$idx]}"
GDIR="$(python3 - "$ROOT" "$idx" <<'PY'
import sys
from pathlib import Path
print(Path(sys.argv[1]) / 'gxtb' / f'conf_{int(sys.argv[2]):04d}')
PY
)"
[[ -s "$JOBDIR/job.inp" && ! -e "$JOBDIR/job.out" ]] || {
    echo "Missing SP input or existing output in $JOBDIR. Refusing overwrite." >&2; exit 3;
}
python3 "$ROOT/prescreen.py" _validate-gxtb \
    --xyz "$GDIR/xtbopt.xyz" --reference "$GDIR/start.xyz" --log "$GDIR/gxtb.out"
LOCK="$JOBDIR/.run_lock"
mkdir "$LOCK" || { echo "Lock exists: $LOCK" >&2; exit 4; }
printf '%s\n' "${SLURM_JOB_ID}_${idx}" > "$LOCK/jobid"
SCR="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/${USER}_presp_${SLURM_JOB_ID}_${idx}_XXXXXX")"
WATCH_PID=""
finalize() {
    local rc="$1" file
    trap - EXIT
    set +e
    if [[ -n "$WATCH_PID" ]]; then kill "$WATCH_PID" 2>/dev/null; wait "$WATCH_PID" 2>/dev/null; fi
    for file in job.out job.gbw job.property.txt; do
        if [[ -f "$SCR/$file" ]]; then
            cp -p "$SCR/$file" "$JOBDIR/$file" || rc=97
        fi
    done
    printf '%s\n' "$rc" > "$JOBDIR/exit_code.txt" || rc=98
    rm -rf -- "$SCR"
    rm -f -- "$LOCK/jobid"
    rmdir -- "$LOCK" 2>/dev/null || true
    exit "$rc"
}
trap 'finalize $?' EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
cp "$GDIR/xtbopt.xyz" "$JOBDIR/start.xyz"
cp "$JOBDIR/start.xyz" "$SCR/start.xyz"
cp "$JOBDIR/job.inp" "$SCR/job.inp"
ORCA="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["orca_binary"])' "$ROOT/settings.json")"
[[ -x "$ORCA" ]] || { echo "Missing ORCA: $ORCA" >&2; exit 5; }
[[ "${SLURM_NTASKS:-0}" == "8" ]] || { echo "ORCA needs 8 MPI tasks" >&2; exit 6; }
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMPI_MCA_btl='^openib'
(
    trap 'exit 0' TERM
    while true; do
        sleep 60 &
        wait $! || exit 0
        left="$(squeue -h -j "$SLURM_JOB_ID" -o '%L' 2>/dev/null | head -n1)"
        if [[ "$left" =~ ^([0-9]+):([0-9][0-9]):([0-9][0-9])$ ]]; then
            sec=$((10#${BASH_REMATCH[1]}*3600+10#${BASH_REMATCH[2]}*60+10#${BASH_REMATCH[3]}))
            if (( sec <= 600 )); then
                [[ -f "$SCR/job.out" ]] && cp -p "$SCR/job.out" "$JOBDIR/job.prewalltime.out" || true
                [[ -f "$SCR/job.gbw" ]] && cp -p "$SCR/job.gbw" "$JOBDIR/job.prewalltime.gbw" || true
                break
            fi
        fi
    done
) >/dev/null 2>&1 &
WATCH_PID="$!"
cd "$SCR"
# ORCA launches its own MPI processes. Do NOT prefix this with mpirun or srun.
"$ORCA" job.inp > job.out 2>&1
python3 "$ROOT/prescreen.py" _validate-sp --log "$SCR/job.out"
