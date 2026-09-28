#!/bin/bash
# M2 trial run 4: bench runner ON THE LOGIN NODE (the model endpoint is
# only reachable here - compute nodes got errno 111), with the chemistry
# acceptance xtb calls routed to dev_cpu_il via DELFIN_CHEM_SLURM=1
# (sbatch --wait per call).  Scratch in a temp dir; artifacts to the
# agent workspace.
set -uo pipefail

WT=/pfs/data6/home/ka/ka_ibcs/ka_ew7404/software/delfin/.delfin/worktrees/delfin-wt-18d68b17
OUT=/home/ka/ka_ibcs/ka_ew7404/agent_workspace/m2-benchmark
mkdir -p "$OUT"

WORK=$(mktemp -d /dev/shm/m2bench_login.XXXXXX)
cleanup() { python -c "import shutil; shutil.rmtree('$WORK', ignore_errors=True)"; }
trap cleanup EXIT

cd "$WT"
export PYTHONPATH="$WT"
export DELFIN_CHEM_SLURM=1

# Bench runner itself; its own scratch home lands in /tmp by default,
# which is fine on the login node.  Run inside $WORK so nothing touches
# the checkout beyond the guarded fixtures.
python -m delfin.agent.cli bench run \
    --model kit.deepseek-v4-flash --backend api --provider kit \
    --task chem_acetaminophen_is_a_minimum,chem_au_cyanide_is_linear \
    --repeats 3 \
    --max-tokens 8192 \
    > "$WORK/bench_stdout.log" 2>&1 || true

cp -v "$WORK/bench_stdout.log" "$OUT/trial4_bench_stdout.log"
LATEST=$(ls -t "$HOME/.delfin/benchmark_runs/"*.jsonl 2>/dev/null | head -1 || true)
if [ -n "${LATEST:-}" ]; then
  cp -v "$LATEST" "$OUT/trial4_run.jsonl"
fi
echo "TRIAL4 DONE"
