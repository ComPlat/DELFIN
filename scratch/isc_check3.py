import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
recs = json.load(open('scratch/baseline.json'))
for r in recs:
    if r['run_id'].endswith('Stefan/Emitter8/ESD'):
        c = r.get('collect_esd_results', {})
        print('error:', c.get('error'))
        print('isc keys:', list((c.get('isc') or {}).keys()))
        print('isc sample:', json.dumps((c.get('isc') or {}), sort_keys=True)[:400])
        print('fluor:', json.dumps(c.get('fluor'))[:200])
        break
