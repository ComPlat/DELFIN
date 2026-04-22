"""Run the benchmark pool, each SMILES in its own subprocess, aggregate JSON."""
import json, os, pathlib, subprocess, sys

sys.path.insert(0, '/home/qmchem_max/ComPlat/DELFIN')
from tests.test_isomer_benchmark import SMILES_POOL, BASELINE_PATH, FIXTURES_DIR

WORKER = '/tmp/bench_worker.py'
with open(WORKER, 'w') as fh:
    fh.write("""import sys, json, os
sys.path.insert(0, '/home/qmchem_max/ComPlat/DELFIN')
from tests.test_isomer_benchmark import _collect_metrics
smi = sys.argv[1]
try:
    m = _collect_metrics(smi)
except Exception as e:
    m = {'error': f'{type(e).__name__}: {e}', 'n_output': 0}
print('__JSON__', json.dumps(m), '__END__')
""")

FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
results = {}
for entry in SMILES_POOL:
    name = entry['id']
    smi = entry['smiles']
    print(f'Running {name}...', flush=True)
    try:
        out = subprocess.run(
            ['micromamba', 'run', '-n', 'delfin', 'python', WORKER, smi],
            capture_output=True, text=True, timeout=900
        )
        stdout = out.stdout
        j_start = stdout.find('__JSON__')
        j_end = stdout.find('__END__')
        if j_start >= 0 and j_end > j_start:
            payload = stdout[j_start + 8:j_end].strip()
            results[name] = json.loads(payload)
        else:
            results[name] = {'error': 'no_json_marker', 'n_output': 0}
    except subprocess.TimeoutExpired:
        results[name] = {'error': 'timeout', 'n_output': 0}
    except Exception as e:
        results[name] = {'error': f'{type(e).__name__}: {e}', 'n_output': 0}
    m = results[name]
    print(f'  -> N={m.get("n_output", 0)}  qwc={m.get("quality_weighted_coverage", 0):.2f}  '
          f'{m.get("wall_seconds", 0):.0f}s  err={m.get("error", "")}', flush=True)

with open(BASELINE_PATH, 'w') as fh:
    json.dump(results, fh, indent=2, sort_keys=True)
print(f'\nWrote {BASELINE_PATH}')
