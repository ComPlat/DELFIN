#!/bin/bash
#SBATCH --job-name=m2-netprobe
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/ka/ka_ibcs/ka_ew7404/agent_workspace/m2-benchmark/%x-%j.out
#SBATCH --error=/home/ka/ka_ibcs/ka_ew7404/agent_workspace/m2-benchmark/%x-%j.err

# Is the model endpoint reachable from a compute node? The bench trial
# 7348199 reached the setup but died on "Connection error." for every
# task. This probe uses the SAME client stack (api_client.create_client)
# so the answer is about reachability, not about a hand-rolled request.
set -x
WT=/pfs/data6/home/ka/ka_ibcs/ka_ew7404/software/delfin/.delfin/worktrees/delfin-wt-18d68b17
cd "$WT"
export PYTHONPATH="$WT"
python - <<'PY'
import socket, traceback
print("host:", socket.gethostname())
try:
    from delfin.agent.api_client import create_client
    c = create_client(model="kit.deepseek-v4-flash", backend="api", provider="kit")
    print("client type:", type(c).__name__)
    # Minimal listing/ping via the OpenAI-compatible surface the client uses.
    for attr in ("models", "list_models"):
        if hasattr(c, attr):
            try:
                res = getattr(c, attr)
                fn = res.list if isinstance(res, type(None)) else res
                print("probe via", attr, "->", str(fn)[:200])
                break
            except Exception as e:
                print(attr, "failed:", type(e).__name__, str(e)[:150])
except Exception:
    traceback.print_exc()
PY
echo "PROBE DONE"
