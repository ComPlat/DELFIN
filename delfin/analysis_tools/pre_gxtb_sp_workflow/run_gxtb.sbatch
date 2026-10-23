#!/usr/bin/env bash
#SBATCH --job-name=pre_gxtb_opt
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gxtb-%A_%a.out
#SBATCH --error=logs/gxtb-%A_%a.err
set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Missing SLURM_SUBMIT_DIR}"
: "${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
ROOT="$SLURM_SUBMIT_DIR"
mapfile -t JOBS < "$ROOT/jobs.gxtb.submission.txt"
idx="$SLURM_ARRAY_TASK_ID"
(( idx >= 0 && idx < ${#JOBS[@]} )) || { echo "Invalid g-xTB index $idx" >&2; exit 2; }
JOBDIR="${JOBS[$idx]}"
[[ -s "$JOBDIR/start.xyz" && ! -e "$JOBDIR/xtbopt.xyz" && ! -e "$JOBDIR/gxtb.out" ]] || {
    echo "Missing input or existing output at $JOBDIR. Refusing overwrite." >&2; exit 3;
}
LOCK="$JOBDIR/.run_lock"
mkdir "$LOCK" || { echo "Lock exists: $LOCK" >&2; exit 4; }
printf '%s\n' "${SLURM_JOB_ID}_${idx}" > "$LOCK/jobid"
SCR="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/${USER}_pregxtb_${SLURM_JOB_ID}_${idx}_XXXXXX")"
WATCH_PID=""
finalize() {
    local rc="$1" file
    trap - EXIT
    set +e
    if [[ -n "$WATCH_PID" ]]; then kill "$WATCH_PID" 2>/dev/null; wait "$WATCH_PID" 2>/dev/null; fi
    for file in gxtb.out xtbopt.xyz; do
        if [[ -f "$SCR/$file" ]]; then
            cp -p "$SCR/$file" "$JOBDIR/$file" || rc=97
        fi
    done
    printf '%s\n' "$rc" > "$JOBDIR/exit_code.txt" || rc=98
    rm -rf -- "$SCR"
    rmdir -- "$LOCK" 2>/dev/null || rm -f -- "$LOCK/jobid"
    rmdir -- "$LOCK" 2>/dev/null || true
    exit "$rc"
}
trap 'finalize $?' EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
cp "$JOBDIR/start.xyz" "$SCR/start.xyz"
GXTB="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["gxtb_binary"])' "$ROOT/settings.json")"
CHARGE="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["charge"])' "$ROOT/settings.json")"
UHF="$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["multiplicity"]-1)' "$ROOT/settings.json")"
[[ -x "$GXTB" ]] || { echo "Missing g-xTB: $GXTB" >&2; exit 6; }
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# Periodic progress + a pre-walltime checkpoint of the output log (not restart data).
(
    trap 'exit 0' TERM
    while true; do
        sleep 60 &
        wait $! || exit 0
        left="$(squeue -h -j "$SLURM_JOB_ID" -o '%L' 2>/dev/null | head -n1)"
        if [[ "$left" =~ ^([0-9]+):([0-9][0-9]):([0-9][0-9])$ ]]; then
            sec=$((10#${BASH_REMATCH[1]}*3600+10#${BASH_REMATCH[2]}*60+10#${BASH_REMATCH[3]}))
            if (( sec <= 600 )); then
                [[ -f "$SCR/gxtb.out" ]] && cp -p "$SCR/gxtb.out" "$JOBDIR/gxtb.prewalltime.out" || true
                break
            fi
        fi
    done
) >/dev/null 2>&1 &
WATCH_PID="$!"
cd "$SCR"
"$GXTB" start.xyz --gxtb --opt --chrg "$CHARGE" --uhf "$UHF" > gxtb.out 2>&1
python3 "$ROOT/prescreen.py" _validate-gxtb --xyz "$SCR/xtbopt.xyz" --reference "$JOBDIR/start.xyz" --log "$SCR/gxtb.out"
