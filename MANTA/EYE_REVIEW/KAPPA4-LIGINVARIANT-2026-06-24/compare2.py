# parallel driver: run k4_worker.py per structure with timeout, aggregate composite
import json, subprocess, statistics, concurrent.futures as cf

sample = json.load(open("/tmp/rigid_sample.json"))

def one(e):
    try:
        p = subprocess.run(["python", "/tmp/k4_worker.py", e["smi"]],
                           cwd="/tmp/delfin_k4", capture_output=True, text=True, timeout=180)
        line = p.stdout.strip().splitlines()
        if not line:
            return {**e, "err": "no-output", "stderr": p.stderr[-200:]}
        d = json.loads(line[-1])
        return {**e, **d}
    except subprocess.TimeoutExpired:
        return {**e, "err": "timeout"}
    except Exception as ex:
        return {**e, "err": str(ex)[:80]}

res = []
with cf.ThreadPoolExecutor(max_workers=24) as ex:
    for r in ex.map(one, sample):
        res.append(r)
json.dump(res, open("/tmp/compare_result.json", "w"))

# ---- summary ----
ok = [r for r in res if "off" in r and "on" in r]
err = [r for r in res if "err" in r]
print(f"=== COMPARE over {len(res)} rigid-polydentate (OFF=current vs ON=cavity-rigid) ===")
print(f"completed {len(ok)}  errors {len(err)} ({[ (r['ref'],r['err']) for r in err][:8]})")
nb_off = sum(1 for r in ok if r["off"].get("n_iso", 0) > 0)
nb_on = sum(1 for r in ok if r["on"].get("n_iso", 0) > 0)
print(f"BUILDS:   OFF {nb_off}/{len(ok)}   ON {nb_on}/{len(ok)}")

def cs(tag):
    return [r[tag]["cshm"] for r in ok if isinstance(r[tag].get("cshm"), float) and r[tag]["cshm"] == r[tag]["cshm"]]
for tag in ("off", "on"):
    c = cs(tag)
    if c:
        print(f"CShM {tag.upper():3}: n={len(c)} median={statistics.median(c):.2f} mean={statistics.mean(c):.2f} "
              f"<2:{sum(1 for x in c if x<2)} <5:{sum(1 for x in c if x<5)} >10:{sum(1 for x in c if x>10)}")

pairs = [(r['off']['cshm'], r['on']['cshm']) for r in ok
         if isinstance(r['off'].get('cshm'), float) and isinstance(r['on'].get('cshm'), float)
         and r['off']['cshm'] == r['off']['cshm'] and r['on']['cshm'] == r['on']['cshm']]
if pairs:
    bett = sum(1 for o, n in pairs if n < o - 0.5); wor = sum(1 for o, n in pairs if n > o + 0.5)
    print(f"PAIRED CShM (n={len(pairs)}): ON better {bett}, ON worse {wor}, ~tie {len(pairs)-bett-wor}")
    print(f"  mean CShM OFF={statistics.mean([o for o,_ in pairs]):.2f} ON={statistics.mean([n for _,n in pairs]):.2f}")
ddp = [(r['off'].get('dd_rmsd'), r['on'].get('dd_rmsd')) for r in ok
       if isinstance(r['off'].get('dd_rmsd'), float) and isinstance(r['on'].get('dd_rmsd'), float)]
if ddp:
    bett = sum(1 for o, n in ddp if n < o - 0.02); wor = sum(1 for o, n in ddp if n > o + 0.02)
    print(f"LIGAND donor-cavity RMSD (n={len(ddp)}): ON better {bett}, ON worse {wor}; mean OFF={statistics.mean([o for o,_ in ddp]):.3f} ON={statistics.mean([n for _,n in ddp]):.3f}")
print("--- per structure: OFF(niso,CShM,dd) | ON(niso,CShM,dd) ---")
for r in sorted(ok, key=lambda x: (x['maxdent'], x['cn'])):
    def fmt(t):
        x = r[t]; c = f"{x['cshm']:.1f}" if isinstance(x.get('cshm'), float) and x['cshm'] == x['cshm'] else "-"
        dd = f"{x['dd_rmsd']:.2f}" if isinstance(x.get('dd_rmsd'), float) else "-"
        return f"n{x.get('n_iso',0)} C={c:>5} dd={dd}"
    print(f"  {r['ref']:10} {r['metal']:3} CN{r['cn']} κ{r['maxdent']}  OFF[{fmt('off')}] | ON[{fmt('on')}]")
for r in err:
    print(f"  {r['ref']:10} {r['metal']:3} CN{r['cn']} κ{r['maxdent']}  ERR={r['err']}")
