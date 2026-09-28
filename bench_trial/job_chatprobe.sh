#!/bin/bash
#SBATCH --job-name=m2-chatprobe
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/ka/ka_ibcs/ka_ew7404/agent_workspace/m2-benchmark/%x-%j.out
#SBATCH --error=/home/ka/ka_ibcs/ka_ew7404/agent_workspace/m2-benchmark/%x-%j.err

# One real chat request from a compute node, through the same client
# stack the bench runner uses. Settles whether trial 7348199's
# "Connection error." was transient or the endpoint is unreachable
# from compute nodes.
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
    cl = getattr(c, "client", None) or getattr(c, "_client", None)
    print("inner client:", type(cl).__name__ if cl is not None else None)
    if cl is not None and hasattr(cl, "chat"):
        r = cl.chat.completions.create(
            model="kit.deepseek-v4-flash",
            messages=[{"role": "user", "content": "Reply with the single word: pong"}],
            max_tokens=5)
        print("REPLY:", r.choices[0].message.content)
        print("REACHABLE")
except Exception:
    traceback.print_exc()
    print("UNREACHABLE-OR-ERROR")
PY
echo "PROBE DONE"
