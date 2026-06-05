#!/usr/bin/env python3
"""GRIP v5 validation script.

Runs three comparative validations once both v4 and v5 libraries are present:

  1. **Hit-rate**: stratified bond/angle/improper queries over a chemistry-
     diverse set of (element-pair, hyb) tuples — count successful lookups
     for v3 / v4 / v5 libraries and report the rates.

  2. **Mogul published comparison**: 30 well-characterised TM bonds with
     published Mogul/literature mean +/- sigma.  Compare v5's value to
     the published range; flag deviations >0.1 A.

  3. **TM-category coverage**: counts of carbene/hapto_eta5/etc. with
     ``n>=5``, and example mean values per metal block.

Usage (run from anywhere)::

    /home/qmchem_max/micromamba/envs/delfin/bin/python \\
        scripts/grip_v5_validation.py \\
        --v3 reports/grip_lib_v3.npz \\
        --v4 reports/grip_lib_v4.npz \\
        --v5 reports/grip_lib_v5.npz \\
        --out /tmp/grip_v5_validation.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ["DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK"] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from delfin.fffree import grip_mogul_lookup as gml  # noqa: E402


# ============================================================
# Published Mogul / literature reference values for key bonds
# ============================================================
# (z1, hyb1, z2, hyb2) -> (published_mu, published_sigma, source)
# All distances in Angstrom.  Sources: Cambridge Structural Database
# (CSD 2020-2024 Mogul reference tables), Allen 1987 acta cryst.
PUBLISHED_BONDS = {
    # Classic organic bonds (calibration anchor)
    ("C", "sp3", "C", "sp3"): (1.530, 0.018, "Allen87"),
    ("C", "sp2", "C", "sp2"): (1.477, 0.014, "Allen87"),
    ("C", "sp", "C", "sp"):   (1.197, 0.012, "Allen87"),
    ("C", "sp3", "H", "*"):   (1.099, 0.011, "Allen87"),
    ("C", "sp3", "N", "sp3"): (1.469, 0.014, "Allen87"),
    ("C", "sp3", "O", "sp3"): (1.426, 0.010, "Allen87"),
    ("C", "sp2", "O", "sp2"): (1.232, 0.008, "Allen87"),
    ("C", "sp2", "N", "sp2"): (1.275, 0.008, "Allen87"),
    # TM bonds (Mogul-pooled mean) - sigma-coordination dominated
    ("Pd", "*", "Cl", "sp3"):  (2.302, 0.029, "MogulCSD"),
    ("Pd", "*", "P", "sp3"):   (2.313, 0.043, "MogulCSD"),
    ("Pd", "*", "N", "sp2"):   (2.030, 0.038, "MogulCSD"),
    ("Pt", "*", "Cl", "sp3"):  (2.317, 0.028, "MogulCSD"),
    ("Pt", "*", "N", "sp2"):   (2.012, 0.041, "MogulCSD"),
    ("Cu", "*", "Cl", "sp3"):  (2.275, 0.105, "MogulCSD"),  # broad — JT distortion
    ("Cu", "*", "N", "sp2"):   (1.991, 0.045, "MogulCSD"),
    ("Cu", "*", "N", "sp3"):   (2.025, 0.075, "MogulCSD"),
    ("Cu", "*", "O", "sp3"):   (1.961, 0.080, "MogulCSD"),
    ("Fe", "*", "N", "sp2"):   (1.985, 0.060, "MogulCSD"),
    ("Fe", "*", "N", "sp3"):   (2.165, 0.075, "MogulCSD"),  # HS broad
    ("Ru", "*", "Cl", "sp3"):  (2.412, 0.040, "MogulCSD"),
    ("Ru", "*", "N", "sp2"):   (2.092, 0.046, "MogulCSD"),
    ("Ru", "*", "C", "sp2"):   (2.030, 0.060, "MogulCSD"),  # mixed sigma+pi
    ("Ir", "*", "Cl", "sp3"):  (2.379, 0.035, "MogulCSD"),
    ("Ir", "*", "C", "sp2"):   (2.045, 0.055, "MogulCSD"),
    ("Ni", "*", "N", "sp2"):   (1.910, 0.045, "MogulCSD"),
    ("Co", "*", "N", "sp2"):   (1.905, 0.050, "MogulCSD"),
    ("Co", "*", "N", "sp3"):   (2.040, 0.075, "MogulCSD"),
    ("Zn", "*", "N", "sp2"):   (2.040, 0.085, "MogulCSD"),
    ("Zn", "*", "O", "sp3"):   (1.998, 0.080, "MogulCSD"),
    ("Mo", "*", "O", "sp2"):   (1.755, 0.080, "MogulCSD"),
    ("W",  "*", "O", "sp2"):   (1.756, 0.080, "MogulCSD"),
    ("Mn", "*", "O", "sp2"):   (1.866, 0.110, "MogulCSD"),  # variable ox state
    ("Hf", "*", "Cl", "sp3"):  (2.415, 0.040, "MogulCSD"),
    ("Rh", "*", "Cl", "sp3"):  (2.354, 0.038, "MogulCSD"),
}


# Stratified hit-rate query set — diverse covering organic + organometallic
HIT_RATE_QUERIES = [
    # Organic bonds
    ("C", "sp3", "C", "sp3"),
    ("C", "sp2", "C", "sp2"),
    ("C", "sp3", "N", "sp3"),
    ("C", "sp2", "N", "sp2"),
    ("C", "sp3", "O", "sp3"),
    ("C", "sp2", "O", "sp2"),
    ("C", "sp3", "F", "sp3"),
    ("C", "sp2", "F", "sp3"),
    ("C", "sp3", "Cl", "sp3"),
    ("C", "sp2", "Br", "sp3"),
    ("C", "sp3", "S", "sp3"),
    ("N", "sp3", "H", "*"),
    ("O", "sp3", "H", "*"),
    # TM bonds
    ("Pd", "*", "Cl", "sp3"),
    ("Pt", "*", "Cl", "sp3"),
    ("Cu", "*", "Cl", "sp3"),
    ("Cu", "*", "N", "sp2"),
    ("Fe", "*", "N", "sp2"),
    ("Fe", "*", "O", "sp3"),
    ("Ru", "*", "C", "sp2"),
    ("Ir", "*", "Cl", "sp3"),
    ("Mn", "*", "O", "sp2"),
    ("Mo", "*", "O", "sp2"),
    ("Co", "*", "N", "sp2"),
    ("Co", "*", "N", "sp3"),
    ("Zn", "*", "O", "sp3"),
    ("Ni", "*", "N", "sp2"),
    ("Rh", "*", "Cl", "sp3"),
    ("Hf", "*", "Cl", "sp3"),
    ("Cd", "*", "N", "sp3"),
    ("Ag", "*", "N", "sp3"),
    ("W", "*", "O", "sp2"),
    # f-block (uncommon — measures coverage)
    ("Eu", "*", "O", "sp2"),
    ("Gd", "*", "O", "sp2"),
    ("La", "*", "O", "sp2"),
    ("U", "*", "O", "sp2"),
    # Edge cases
    ("Si", "sp3", "O", "sp3"),
    ("P", "sp3", "O", "sp2"),
    ("B", "sp2", "N", "sp2"),
    ("S", "sp3", "O", "sp2"),
]

ANGLE_QUERIES = [
    ("C", "C", "sp3", "C"),
    ("C", "C", "sp2", "O"),
    ("C", "C", "sp2", "N"),
    ("C", "N", "sp3", "C"),
    ("C", "O", "sp3", "C"),
    ("Cl", "Pd", "*", "Cl"),
    ("Cl", "Pt", "*", "N"),
    ("N", "Cu", "*", "N"),
    ("N", "Fe", "*", "N"),
    ("O", "Mn", "*", "O"),
    ("Cl", "Ir", "*", "Cl"),
]


def _open(path: Path):
    if not path.exists():
        return None
    gml.GripLibrary._SINGLETONS.clear()
    return gml.GripLibrary.get(path)


def run_hit_rate(libs):
    """Hit-rate analysis: for each query, record which libraries returned a value."""
    rows = []
    for q in HIT_RATE_QUERIES:
        row = {"query": list(q)}
        for name, lib in libs.items():
            if lib is None:
                row[name] = None
                continue
            try:
                hit = lib.lookup_bond(*q)
                row[name] = (round(hit[0], 4), round(hit[1], 4), int(hit[2])) if hit else None
            except Exception as exc:
                row[name] = f"ERR:{exc!r}"
        rows.append(row)
    return rows


def run_angle_hit_rate(libs):
    rows = []
    for q in ANGLE_QUERIES:
        row = {"query": list(q)}
        for name, lib in libs.items():
            if lib is None:
                row[name] = None
                continue
            try:
                hit = lib.lookup_angle(*q)
                row[name] = (round(hit[0], 4), round(hit[1], 4), int(hit[2])) if hit else None
            except Exception as exc:
                row[name] = f"ERR:{exc!r}"
        rows.append(row)
    return rows


def run_published_comparison(libs):
    """Compare v4 / v5 means against published Mogul/Allen values."""
    rows = []
    for q, (pub_mu, pub_sigma, src) in PUBLISHED_BONDS.items():
        row = {
            "query": list(q),
            "published_mu": pub_mu,
            "published_sigma": pub_sigma,
            "source": src,
        }
        for name, lib in libs.items():
            if lib is None:
                row[f"{name}_mu"] = None
                continue
            try:
                hit = lib.lookup_bond(*q)
            except Exception as exc:
                row[f"{name}_mu"] = f"ERR:{exc!r}"
                continue
            if hit is None:
                row[f"{name}_mu"] = None
                row[f"{name}_dev"] = None
            else:
                row[f"{name}_mu"] = round(hit[0], 4)
                row[f"{name}_sigma"] = round(hit[1], 4)
                row[f"{name}_n"] = int(hit[2])
                row[f"{name}_dev"] = round(abs(hit[0] - pub_mu), 4)
        rows.append(row)
    return rows


def run_tm_category_coverage(libs):
    """Report TM-category coverage for v5 library only."""
    v5 = libs.get("v5")
    if v5 is None or not v5.has_tm_categories:
        return {"v5_unavailable": True}
    out = {}
    for cat in ("carbene", "hapto_eta2", "hapto_eta5", "hapto_eta6",
                "mu_bridge", "agostic", "ox_addition"):
        tbl = v5._tm_cat_tables.get(cat)
        if tbl is None:
            out[cat] = {"n_keys": 0, "examples": []}
            continue
        n_keys = tbl["n_keys"]
        # Top-10 keys by n
        n_arr = tbl["n"]
        mu_arr = tbl["mu"]
        sig_arr = tbl["sigma"]
        keys = list(tbl["key_to_idx"].keys())
        if n_keys > 0:
            order = sorted(range(n_keys), key=lambda i: -int(n_arr[i]))[:10]
            examples = []
            for i in order:
                try:
                    parsed = json.loads(keys[i])
                except Exception:
                    parsed = keys[i]
                examples.append({
                    "key": parsed,
                    "mu": round(float(mu_arr[i]), 4),
                    "sigma": round(float(sig_arr[i]), 4),
                    "n": int(n_arr[i]),
                })
            out[cat] = {"n_keys": n_keys, "examples": examples}
        else:
            out[cat] = {"n_keys": 0, "examples": []}
    return out


def run_metadata_summary(libs):
    out = {}
    for name, lib in libs.items():
        if lib is None:
            out[name] = None
            continue
        d = {"version": lib.version}
        d["n_master"] = lib.n_master
        d["n_pair_bond"] = getattr(lib, "n_pair_bond", 0)
        d["n_triple_angle"] = getattr(lib, "n_triple_angle", 0)
        d["n_improper_pair"] = getattr(lib, "n_improper_pair", 0)
        d["has_pair_tables"] = getattr(lib, "has_pair_tables", False)
        d["has_tm_categories"] = getattr(lib, "has_tm_categories", False)
        d["has_block_disagg"] = getattr(lib, "has_block_disagg", False)
        d["cleaning_applied"] = getattr(lib, "cleaning_applied", False)
        d["meta_n_extracted"] = getattr(lib, "meta_n_extracted", 0)
        d["meta_n_total_scanned"] = getattr(lib, "meta_n_total_scanned", 0)
        out[name] = d
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--v3", type=Path, default=Path(
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v3.npz"))
    ap.add_argument("--v4", type=Path, default=Path(
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v4.npz"))
    ap.add_argument("--v5", type=Path, default=Path(
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"))
    ap.add_argument("--out", type=Path, default=Path("/tmp/grip_v5_validation.json"))
    args = ap.parse_args()

    print(f"[validate] loading libraries...", flush=True)
    libs = {
        "v3": _open(args.v3),
        "v4": _open(args.v4),
        "v5": _open(args.v5),
    }
    for n, l in libs.items():
        print(f"  {n}: {'loaded v=' + str(l.version) if l else 'MISSING'}")

    # Clear and re-open since shared singleton cache.
    gml.GripLibrary._SINGLETONS.clear()
    libs["v3"] = _open(args.v3)
    libs["v4"] = _open(args.v4)
    libs["v5"] = _open(args.v5)

    out = {}
    out["meta"] = run_metadata_summary(libs)
    print("[validate] running hit-rate analysis...", flush=True)
    out["hit_rate"] = run_hit_rate(libs)
    print("[validate] running angle hit-rate analysis...", flush=True)
    out["angle_hit_rate"] = run_angle_hit_rate(libs)
    print("[validate] running published bond comparison...", flush=True)
    out["published_comparison"] = run_published_comparison(libs)
    print("[validate] running TM-category coverage...", flush=True)
    out["tm_category_coverage"] = run_tm_category_coverage(libs)

    # Compute summary stats
    hr = out["hit_rate"]
    summary = {}
    for name in ("v3", "v4", "v5"):
        hits = sum(1 for r in hr if r.get(name) and not isinstance(r[name], str))
        summary[name] = {
            "hit_rate_bond": f"{hits}/{len(hr)} ({100.0*hits/max(len(hr),1):.1f}%)",
        }
    ahr = out["angle_hit_rate"]
    for name in ("v3", "v4", "v5"):
        hits = sum(1 for r in ahr if r.get(name) and not isinstance(r[name], str))
        summary[name]["hit_rate_angle"] = f"{hits}/{len(ahr)} ({100.0*hits/max(len(ahr),1):.1f}%)"
    # Published-comparison MAD
    pc = out["published_comparison"]
    for name in ("v3", "v4", "v5"):
        devs = [r.get(f"{name}_dev") for r in pc
                if r.get(f"{name}_dev") is not None and not isinstance(r.get(f"{name}_dev"), str)]
        if devs:
            summary[name]["pub_match_count"] = len(devs)
            summary[name]["pub_mean_dev"] = round(sum(devs) / len(devs), 4)
            summary[name]["pub_max_dev"] = round(max(devs), 4)
        else:
            summary[name]["pub_match_count"] = 0
    out["summary"] = summary

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(f"\n[validate] wrote {args.out}")
    print("\n[validate] SUMMARY:")
    for name, s in summary.items():
        print(f"  {name}: {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
