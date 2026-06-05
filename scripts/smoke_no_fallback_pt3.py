"""MISSION A5 — smoke comparison: NO_FALLBACK ON vs OFF on a small SMILES pool.

For each input SMILES we run :func:`delfin.smiles_converter.smiles_to_xyz_isomers`
twice -- once with ``DELFIN_FFFREE_NO_FALLBACK=1`` (auto-implies FF-free
dispatch) and once with ``DELFIN_FFFREE_BUILDER=1`` (fallback allowed) --
under a forensik log that records native / blocked / fallback per call.

Outputs (printed to stdout, written to --out-dir):
    * ``forensik_on.tsv``, ``forensik_off.tsv`` -- raw breadcrumb files
    * ``summary.json``                          -- per-class native/blocked/fallback
                                                   rates + coverage delta

Determinism: each (smi, run) pair is independent; no parallelism by default
to keep the breadcrumb file append-ordered.  ``--parallel`` is exposed for
larger smokes (the log file uses O_APPEND so multi-process appends are
line-atomic on POSIX).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple


def _load_pool(path: Path, limit: int = 0) -> List[Tuple[str, str]]:
    """Load ``id|smiles`` lines.  ``limit > 0`` truncates the head."""
    rows: List[Tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "|" in line:
                sid, smi = line.split("|", 1)
            else:
                sid, smi = f"row{len(rows)+1}", line
            rows.append((sid.strip(), smi.strip()))
            if limit and len(rows) >= limit:
                break
    return rows


_WORKER_SCRIPT = """
import os, sys, json, traceback
os.environ['PYTHONHASHSEED'] = '0'
{envset}
try:
    from delfin import smiles_converter as sc
    smi = {smi!r}
    res, err = sc.smiles_to_xyz_isomers(smi, max_isomers=20, num_confs=10)
    n = len(res) if res else 0
    print(json.dumps({{'sid': {sid!r}, 'smi': smi, 'n': n, 'err': err}}))
except Exception as exc:
    print(json.dumps({{'sid': {sid!r}, 'smi': {smi!r}, 'n': 0, 'err': str(exc), 'tb': traceback.format_exc()}}))
"""


def _run_one(sid: str, smi: str, mode: str, forensik_log: str,
             timeout: int = 120) -> Dict:
    """Run one SMILES in a child process with env injected for ``mode``."""
    if mode == "on":
        envset = (
            "os.environ['DELFIN_FFFREE_NO_FALLBACK'] = '1'\n"
            f"os.environ['DELFIN_FFFREE_FORENSIK_LOG'] = {forensik_log!r}\n"
        )
    elif mode == "off":
        envset = (
            "os.environ['DELFIN_FFFREE_BUILDER'] = '1'\n"
            f"os.environ['DELFIN_FFFREE_FORENSIK_LOG'] = {forensik_log!r}\n"
        )
    else:
        raise ValueError(f"unknown mode {mode!r}")
    src = _WORKER_SCRIPT.format(envset=envset, sid=sid, smi=smi)
    try:
        cp = subprocess.run(
            [sys.executable, "-c", src],
            capture_output=True, text=True, timeout=timeout,
        )
        out = cp.stdout.strip().splitlines()
        if not out:
            return {"sid": sid, "smi": smi, "n": 0,
                    "err": f"empty stdout: {cp.stderr[:500]!r}"}
        return json.loads(out[-1])
    except subprocess.TimeoutExpired:
        return {"sid": sid, "smi": smi, "n": 0, "err": "timeout"}
    except Exception as exc:
        return {"sid": sid, "smi": smi, "n": 0, "err": str(exc)}


def _classify(sid: str) -> str:
    """Heuristic class from sid prefix.  Falls back to ``other``."""
    s = sid.lower()
    if any(t in s for t in ("hapto", "ferrocene", "cp", "eta", "η", "h5", "h6")):
        return "hapto"
    if any(t in s for t in ("multi", "bimet", "trimet", "cluster")):
        return "multi"
    if any(t in s for t in ("organ", "no_metal", "ligand")):
        return "organic"
    return "sigma"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True, help="SMILES pool file")
    ap.add_argument("--limit", type=int, default=300,
                    help="Truncate pool to first N rows (0 = full)")
    ap.add_argument("--parallel", type=int, default=40,
                    help="Parallel worker processes")
    ap.add_argument("--timeout", type=int, default=120,
                    help="Per-SMILES timeout (s)")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    pool = _load_pool(Path(args.pool), args.limit)
    print(f"[smoke] {len(pool)} SMILES -- parallel={args.parallel}, "
          f"timeout={args.timeout}s", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_on = out_dir / "forensik_on.tsv"
    log_off = out_dir / "forensik_off.tsv"
    log_on.write_text("")
    log_off.write_text("")

    # Run both modes in sequence (clearer determinism, simpler attribution).
    results: Dict[str, Dict[str, List[Dict]]] = {"on": {}, "off": {}}
    for mode, log in (("on", log_on), ("off", log_off)):
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=args.parallel) as ex:
            futs = {ex.submit(_run_one, sid, smi, mode, str(log),
                              args.timeout): (sid, smi)
                    for sid, smi in pool}
            for k, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                cls = _classify(rec["sid"])
                results[mode].setdefault(cls, []).append(rec)
                if k % 25 == 0:
                    print(f"[smoke/{mode}] {k}/{len(pool)} "
                          f"({time.time()-t0:.1f}s)", flush=True)
        print(f"[smoke/{mode}] DONE in {time.time()-t0:.1f}s", flush=True)

    # Aggregate
    summary: Dict[str, Dict] = {}
    overall = {"on": Counter(), "off": Counter()}
    for mode in ("on", "off"):
        per_class = {}
        for cls, recs in results[mode].items():
            n_total = len(recs)
            n_emit = sum(1 for r in recs if r["n"] > 0)
            n_zero = n_total - n_emit
            tot_iso = sum(r["n"] for r in recs)
            per_class[cls] = {
                "n_smiles": n_total,
                "n_emit": n_emit,
                "n_zero": n_zero,
                "coverage_pct": round(100.0 * n_emit / max(1, n_total), 2),
                "n_isomers": tot_iso,
            }
            overall[mode]["n_smiles"] += n_total
            overall[mode]["n_emit"] += n_emit
            overall[mode]["n_zero"] += n_zero
            overall[mode]["n_isomers"] += tot_iso
        summary[mode] = per_class

    # Coverage delta
    delta = {}
    on_overall = overall["on"]
    off_overall = overall["off"]
    delta["coverage_pp"] = (
        round(100.0 * on_overall["n_emit"] / max(1, on_overall["n_smiles"]), 2)
        - round(100.0 * off_overall["n_emit"] / max(1, off_overall["n_smiles"]), 2)
    )
    delta["isomers_abs"] = on_overall["n_isomers"] - off_overall["n_isomers"]

    # Forensik breadcrumb aggregation
    forensik = {"on": Counter(), "off": Counter()}
    for mode, log in (("on", log_on), ("off", log_off)):
        if log.exists():
            for line in log.read_text().splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    forensik[mode][parts[1]] += 1
    summary["forensik_breadcrumb"] = {
        "on": dict(forensik["on"]),
        "off": dict(forensik["off"]),
    }
    summary["overall"] = {
        "on": dict(on_overall),
        "off": dict(off_overall),
        "delta": delta,
    }

    out_json = out_dir / "summary.json"
    out_json.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary["overall"], indent=2), flush=True)
    print(f"[smoke] forensik breadcrumb: {summary['forensik_breadcrumb']}",
          flush=True)
    print(f"[smoke] summary: {out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
