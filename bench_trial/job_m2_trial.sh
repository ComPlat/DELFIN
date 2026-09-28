#!/bin/bash
#SBATCH --job-name=m2-bench-trial
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# M2 chemistry benchmark trial: 2 tasks x 3 repeats.
# Scratch stays node-local ($TMPDIR); artifacts go to the agent
# workspace, NOT the git checkout (the suite's own hygiene check fails
# on stray files in the tree).  The bench runner runs from the worktree
# checkout so its packaged tasks and code resolve to THIS branch.
set -euo pipefail

WT=/pfs/data6/home/ka/ka_ibcs/ka_ew7404/software/delfin/.delfin/worktrees/delfin-wt-18d68b17
OUT=/home/ka/ka_ibcs/ka_ew7404/agent_workspace/m2-benchmark
mkdir -p "$OUT"

export TMPDIR=${TMPDIR:-/dev/shm}
WORK=$TMPDIR/m2bench_${SLURM_JOB_ID:-$$}
mkdir -p "$WORK"

cd "$WT"
export PYTHONPATH="$WT"

python -m delfin.agent.cli bench run \
    --model kit.deepseek-v4-flash --backend api --provider kit \
    --task chem_acetaminophen_is_a_minimum,chem_au_cyanide_is_linear \
    --repeats 3 \
    --max-tokens 8192 \
    > "$WORK/bench_stdout.log" 2>&1 || true

cp -v "$WORK/bench_stdout.log" "$OUT/"
LATEST=$(ls -t "$HOME/.delfin/benchmark_runs/"*.jsonl 2>/dev/null | head -1 || true)
if [ -n "${LATEST:-}" ]; then
  cp -v "$LATEST" "$OUT/trial_run.jsonl"
fi
echo "TRIAL DONE"
