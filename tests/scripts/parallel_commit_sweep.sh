#!/bin/bash
# Parallel commit sweep: each commit gets its own git worktree + runs
# benchmark concurrently.  Writes per-commit JSON to /tmp/commit_sweep/.
set -u
RESULTS_DIR=/tmp/commit_sweep
mkdir -p "$RESULTS_DIR"
LOG_DIR=/tmp/commit_sweep_logs
mkdir -p "$LOG_DIR"
SRC_REPO=/home/qmchem_max/ComPlat/DELFIN

run_one_commit() {
    local sha=$1
    local out_json="$RESULTS_DIR/$sha.json"
    local log="$LOG_DIR/$sha.log"
    local wt="/tmp/sweep_wt_$sha"
    if [ -s "$out_json" ]; then
        echo "SKIP $sha (exists)" >> "$LOG_DIR/_main.log"
        return
    fi
    # Isolated worktree checkout.
    if [ ! -d "$wt" ]; then
        (cd "$SRC_REPO" && git worktree add "$wt" "$sha" 2>>"$LOG_DIR/_main.log") || {
            echo "WT_FAIL $sha" >> "$LOG_DIR/_main.log"
            return
        }
    fi
    # Run benchmark using the worktree's own smiles_converter.py.
    timeout 18000 micromamba run -n delfin python -c "
import sys, json, subprocess, pathlib, time
sys.path.insert(0, '$wt')
# Import POOL from HEAD's test module (consistent SMILES list).
sys.path.insert(0, '$SRC_REPO')
from tests.test_isomer_benchmark import SMILES_POOL

results = {}
for entry in SMILES_POOL:
    name = entry['id']; smi = entry['smiles']
    t0 = time.time()
    worker = f'''
import sys, json, time
sys.path.insert(0, '{'"$wt"'}')
t0 = time.time()
try:
    from delfin.smiles_converter import smiles_to_xyz_isomers
    for kwargs in (
        dict(apply_uff=True, deterministic=True, collapse_label_variants=True, quality_mode='normal'),
        dict(apply_uff=True, deterministic=True, collapse_label_variants=True),
        dict(apply_uff=True),
    ):
        try:
            r, err = smiles_to_xyz_isomers({smi!r}, **kwargs)
            break
        except TypeError:
            continue
    else:
        r, err = [], 'no compatible signature'
    n = len(r) if r else 0
    print('__JSON__', json.dumps({{'n_output': n, 'wall_seconds': round(time.time()-t0,1), 'error': err if err else ''}}), '__END__')
except Exception as e:
    print('__JSON__', json.dumps({{'n_output': 0, 'wall_seconds': round(time.time()-t0,1), 'error': f'{{type(e).__name__}}: {{e}}'}}), '__END__')
'''
    try:
        out = subprocess.run(['micromamba','run','-n','delfin','python','-c',worker],
                             capture_output=True, text=True, timeout=2100)
        stdout = out.stdout
        js = stdout.find('__JSON__'); je = stdout.find('__END__')
        if js >= 0 and je > js:
            results[name] = json.loads(stdout[js+8:je].strip())
        else:
            results[name] = {'n_output':0,'error':'no_marker'}
    except subprocess.TimeoutExpired:
        results[name] = {'n_output':0,'error':'timeout_35min'}
    except Exception as e:
        results[name] = {'n_output':0,'error':f'{type(e).__name__}'}
    m = results[name]
    print(f'  {name:35s} N={m.get(\"n_output\",0):3d} {m.get(\"wall_seconds\",0):.0f}s {m.get(\"error\",\"\")}', flush=True)

with open('$out_json', 'w') as fh:
    json.dump(results, fh, indent=2, sort_keys=True)
print('DONE', '$sha')
" > "$log" 2>&1
    echo "FINISHED $sha" >> "$LOG_DIR/_main.log"
}

# Launch each commit in background.
for sha in "$@"; do
    run_one_commit "$sha" &
    # Small delay between launches so worktree creation doesn't race.
    sleep 5
done
wait
echo "ALL DONE"
