"""Mogul v4 Forensik: quantify v3-anomalies that fall into v4-extension buckets.

For a sample of XYZ files in a pool we:
  1. Run mogul_detector_v3 in CCDC mode
  2. For each anomaly, classify it into one or more v4-categories:
       (a) X-H bond  (involves H)
       (b) D-M-D angle (3 atoms, central is metal, both outer are donors)
       (c) M-D-X angle (3 atoms, the metal is one of the outer atoms)
       (d) torsion in bimodal/trimodal GMM region (v5 says ≥2 components)
  3. Quantify how many would be re-classified by v4-extensions (i.e.
     "false anomalies" using legacy single-Gaussian thresholds when
     reality is multimodal, or H-specific thresholds).

Output:
  paper_data/mogul_v4_forensik_report.json
  paper_data/mogul_v4_forensik_summary.md
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN")
import delfin._bond_decollapse as _bd  # noqa: E402
import delfin.fffree.mogul_detector_v3 as v3  # noqa: E402

V5_PATH = "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"

METALS = set(
    [
        "Li", "Na", "K", "Rb", "Cs",
        "Be", "Mg", "Ca", "Sr", "Ba",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi",
        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu",
    ]
)


def parse_xyz(path: str):
    """Read first frame from an XYZ file. Returns (syms, P_Nx3)."""
    text = Path(path).read_text().splitlines()
    if not text:
        return [], np.zeros((0, 3))
    na = int(text[0].strip())
    syms: List[str] = []
    P = np.zeros((na, 3), dtype=float)
    for i in range(na):
        parts = text[2 + i].split()
        syms.append(parts[0])
        P[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return syms, P


class V5Lookup:
    """Wraps the v5 lib in convenient lookup helpers."""

    def __init__(self, path: str = V5_PATH):
        self.lib = np.load(path, allow_pickle=True)
        # build dicts for fast lookup
        self._pair_bond = self._build_dict(
            self.lib["pair_bond_keys"],
            self.lib["pair_bond_mu"],
            self.lib["pair_bond_sigma"],
            self.lib["pair_bond_n"],
        )
        self._triple_angle = self._build_dict(
            self.lib["triple_angle_keys"],
            self.lib["triple_angle_mu"],
            self.lib["triple_angle_sigma"],
            self.lib["triple_angle_n"],
        )
        # torsion GMM
        self._torsion_keys = self.lib["torsion_keys"]
        self._torsion_n_components = self.lib["torsion_n_components"]
        self._torsion_pi = self.lib["torsion_pi"]
        self._torsion_mu = self.lib["torsion_mu"]
        self._torsion_sigma = self.lib["torsion_sigma"]
        # build a key-only lookup for torsion: use sorted-elements as fallback
        self._torsion_simple = {}
        for i, k in enumerate(self._torsion_keys):
            try:
                arr = json.loads(k)
                elems = sorted([arr[0], arr[2], arr[4], arr[5]])
                kk = "-".join(elems)
                if kk not in self._torsion_simple:
                    self._torsion_simple[kk] = i
            except Exception:
                continue
        # TM blocks
        self._tm_pair = self._build_dict(
            self.lib["tm_pair_bond_block_keys"],
            self.lib["tm_pair_bond_block_mu"],
            self.lib["tm_pair_bond_block_sigma"],
            self.lib["tm_pair_bond_block_n"],
        )
        self._tm_triple = self._build_dict(
            self.lib["tm_triple_angle_block_keys"],
            self.lib["tm_triple_angle_block_mu"],
            self.lib["tm_triple_angle_block_sigma"],
            self.lib["tm_triple_angle_block_n"],
        )

    @staticmethod
    def _build_dict(keys, mu, sigma, n):
        out = {}
        for i, k in enumerate(keys):
            out[str(k)] = (float(mu[i]), float(sigma[i]), int(n[i]))
        return out

    def pair_bond_xh(self, h_sym: str, x_sym: str, x_hyb: str = "*"):
        """Return (mu, sigma, n) for X-H bond. h_sym is the H side; x_sym
        with hybridization x_hyb is the heavy side. Key form sorts
        lexicographically by the JSON in the lib."""
        # try both orderings
        h_hyb = "*"
        a, b = sorted([(h_sym, h_hyb), (x_sym, x_hyb)])
        k = json.dumps([a[0], a[1], b[0], b[1]], separators=(",", ":"))
        if k in self._pair_bond:
            return self._pair_bond[k]
        # try with x_hyb='*' too
        a, b = sorted([(h_sym, h_hyb), (x_sym, "*")])
        k = json.dumps([a[0], a[1], b[0], b[1]], separators=(",", ":"))
        return self._pair_bond.get(k)

    def triple_angle(self, center: str, outer1: str, outer2: str, center_hyb: str = "*"):
        """Return (mu, sigma, n) for outer1-center-outer2 angle.
        Key form: [center, center_hyb, sorted(outer1, outer2)].
        """
        wings = sorted([outer1, outer2])
        k = json.dumps([center, center_hyb, wings[0], wings[1]], separators=(",", ":"))
        if k in self._triple_angle:
            return self._triple_angle[k]
        # fallback wildcard
        k = json.dumps([center, "*", wings[0], wings[1]], separators=(",", ":"))
        return self._triple_angle.get(k)

    def tm_pair_bond(self, m: str, m_block: str, x: str, x_hyb: str = "*"):
        a, b = sorted([(m, m_block), (x, x_hyb)])
        k = json.dumps([a[0], a[1], b[0], b[1]], separators=(",", ":"))
        return self._tm_pair.get(k)

    def tm_triple_angle(self, m: str, m_block: str, outer1: str, outer2: str):
        wings = sorted([outer1, outer2])
        k = json.dumps([m, m_block, wings[0], wings[1]], separators=(",", ":"))
        return self._tm_triple.get(k)

    def torsion_gmm(self, e0: str, e1: str, e2: str, e3: str):
        """Lookup via simple sorted-element key. Returns (n_comp, pi, mu, sigma) or None."""
        kk = "-".join(sorted([e0, e1, e2, e3]))
        idx = self._torsion_simple.get(kk)
        if idx is None:
            return None
        return (
            int(self._torsion_n_components[idx]),
            self._torsion_pi[idx],
            self._torsion_mu[idx],
            self._torsion_sigma[idx],
        )


def metal_block(sym: str) -> str:
    """Rough metal block classification (matches v5 lib's choice)."""
    s = sym
    if s in ("Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"):
        return "3d"
    if s in ("Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"):
        return "4d"
    if s in ("La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"):
        return "5d"
    if s in ("Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"):
        return "4f"
    return "*"


def classify_anomaly(rec: dict, syms: Sequence[str], P: np.ndarray, v5: V5Lookup) -> dict:
    """Return a dict describing v4-categories that apply to this v3 record."""
    axis = rec["axis"]
    atom_syms = rec["atom_syms"]
    atoms = rec["atoms"]
    obs = rec["obs"]

    cat = {
        "v4_xh": False,
        "v4_dmd": False,
        "v4_mdx": False,
        "v4_torsion_multimodal": False,
        "v5_xh_sev": None,
        "v5_dmd_sev": None,
        "v5_mdx_sev": None,
        "v5_torsion_sev_best_component": None,
        "v5_torsion_n_components": None,
        # would-eliminate flag: anomaly drops below thr=2.5 if we use the
        # v4-specific Gaussian instead of the v3 fragment-pool MAD
        "v4_eliminates": False,
        "v4_new_anomaly": False,
    }

    if axis == "bond":
        a, b = atoms
        sa, sb = syms[a], syms[b]
        if sa == "H" or sb == "H":
            x_sym = sb if sa == "H" else sa
            cat["v4_xh"] = True
            res = v5.pair_bond_xh("H", x_sym)
            if res:
                mu, sigma, n = res
                if sigma > 1e-9 and n >= 15:
                    sev = abs(obs - mu) / sigma
                    cat["v5_xh_sev"] = sev
                    # severity 2.0 ≈ MAD 2.5; flag eliminates if v5 says ≤2
                    if sev < 2.0 and rec.get("sev_mad", 0) >= 2.5:
                        cat["v4_eliminates"] = True

    elif axis == "angle":
        # three atoms; center atom is the one bonded to both wings
        # in v3 emission, the "center" key is the central atom index
        center_idx = rec.get("center")
        center_sym = syms[center_idx] if center_idx is not None else None
        wings = [i for i in atoms if i != center_idx]
        if len(wings) == 2 and center_sym:
            w1_sym = syms[wings[0]]
            w2_sym = syms[wings[1]]
            center_is_metal = center_sym in METALS
            wing1_is_metal = w1_sym in METALS
            wing2_is_metal = w2_sym in METALS

            if center_is_metal:
                # D-M-D angle: center is metal; wings are donors
                cat["v4_dmd"] = True
                blk = metal_block(center_sym)
                res = v5.tm_triple_angle(center_sym, blk, w1_sym, w2_sym)
                if res is None:
                    res = v5.triple_angle(center_sym, w1_sym, w2_sym)
                if res:
                    mu, sigma, n = res
                    if sigma > 1e-9 and n >= 15:
                        sev = abs(obs - mu) / sigma
                        cat["v5_dmd_sev"] = sev
                        if sev < 2.0 and rec.get("sev_mad", 0) >= 2.5:
                            cat["v4_eliminates"] = True
            elif wing1_is_metal or wing2_is_metal:
                # M-D-X angle: metal is an outer atom
                cat["v4_mdx"] = True
                res = v5.triple_angle(center_sym, w1_sym, w2_sym)
                if res:
                    mu, sigma, n = res
                    if sigma > 1e-9 and n >= 15:
                        sev = abs(obs - mu) / sigma
                        cat["v5_mdx_sev"] = sev
                        if sev < 2.0 and rec.get("sev_mad", 0) >= 2.5:
                            cat["v4_eliminates"] = True

    elif axis == "torsion":
        # Try GMM lookup
        if len(atom_syms) == 4:
            res = v5.torsion_gmm(*atom_syms)
            if res is not None:
                n_comp, pi, mu, sigma = res
                cat["v5_torsion_n_components"] = int(n_comp)
                if n_comp >= 2:
                    cat["v4_torsion_multimodal"] = True
                    # find best component sev (signed observation; v3 uses abs(tor))
                    # try both polarities since we don't know sign
                    best_sev = 1e9
                    for sign in (1.0, -1.0):
                        ob = sign * obs
                        for c in range(n_comp):
                            if not np.isfinite(mu[c]) or not np.isfinite(sigma[c]):
                                continue
                            if sigma[c] < 1e-9:
                                continue
                            if pi[c] < 0.05:
                                continue
                            s = abs(ob - mu[c]) / sigma[c]
                            if s < best_sev:
                                best_sev = s
                    if best_sev < 1e8:
                        cat["v5_torsion_sev_best_component"] = best_sev
                        if best_sev < 2.0 and rec.get("sev_mad", 0) >= 2.5:
                            cat["v4_eliminates"] = True
    return cat


def run_one_xyz(path: str, v5: V5Lookup, idx) -> dict:
    syms, P = parse_xyz(path)
    if len(syms) < 3:
        return {"path": path, "ok": False, "reason": "too_few_atoms", "anomalies": []}
    try:
        recs = v3.detect_anomalies_v3(syms, P, mode="ccdc", index=idx)
    except Exception as e:
        return {"path": path, "ok": False, "reason": f"err:{e}", "anomalies": []}
    classified = []
    for r in recs:
        c = classify_anomaly(r, syms, P, v5)
        classified.append({**r, **c})
    return {"path": path, "ok": True, "anomalies": classified}


def aggregate(results: List[dict]) -> dict:
    total_anomalies = 0
    n_xh = 0
    n_dmd = 0
    n_mdx = 0
    n_torsion_multi = 0
    n_eliminated = 0
    n_files_ok = 0
    n_files_err = 0

    eliminated_by_cat = {"xh": 0, "dmd": 0, "mdx": 0, "torsion": 0}
    severity_change_xh = []
    severity_change_dmd = []
    severity_change_mdx = []
    severity_change_torsion = []

    for r in results:
        if not r.get("ok"):
            n_files_err += 1
            continue
        n_files_ok += 1
        for a in r["anomalies"]:
            total_anomalies += 1
            old_sev = float(a.get("sev_mad", 0))
            elim = a.get("v4_eliminates", False)
            if a.get("v4_xh"):
                n_xh += 1
                if elim:
                    eliminated_by_cat["xh"] += 1
                if a.get("v5_xh_sev") is not None:
                    severity_change_xh.append((old_sev, a["v5_xh_sev"]))
            if a.get("v4_dmd"):
                n_dmd += 1
                if elim:
                    eliminated_by_cat["dmd"] += 1
                if a.get("v5_dmd_sev") is not None:
                    severity_change_dmd.append((old_sev, a["v5_dmd_sev"]))
            if a.get("v4_mdx"):
                n_mdx += 1
                if elim:
                    eliminated_by_cat["mdx"] += 1
                if a.get("v5_mdx_sev") is not None:
                    severity_change_mdx.append((old_sev, a["v5_mdx_sev"]))
            if a.get("v4_torsion_multimodal"):
                n_torsion_multi += 1
                if elim:
                    eliminated_by_cat["torsion"] += 1
                if a.get("v5_torsion_sev_best_component") is not None:
                    severity_change_torsion.append(
                        (old_sev, a["v5_torsion_sev_best_component"])
                    )
            if elim:
                n_eliminated += 1

    def mean_sev(arr):
        if not arr:
            return None
        old = [a for a, b in arr]
        new = [b for a, b in arr]
        return {
            "n": len(arr),
            "old_mean": float(np.mean(old)),
            "new_mean": float(np.mean(new)),
            "old_median": float(np.median(old)),
            "new_median": float(np.median(new)),
        }

    return {
        "n_files_ok": n_files_ok,
        "n_files_err": n_files_err,
        "total_anomalies": total_anomalies,
        "n_xh": n_xh,
        "n_dmd": n_dmd,
        "n_mdx": n_mdx,
        "n_torsion_multi": n_torsion_multi,
        "n_v4_eliminates": n_eliminated,
        "eliminated_by_cat": eliminated_by_cat,
        "severity_change_xh": mean_sev(severity_change_xh),
        "severity_change_dmd": mean_sev(severity_change_dmd),
        "severity_change_mdx": mean_sev(severity_change_mdx),
        "severity_change_torsion": mean_sev(severity_change_torsion),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pool",
        default="/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/2792332-aromatic-symmetry-VOLLPOOL",
    )
    ap.add_argument("--sample", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out",
        default="/home/qmchem_max/ComPlat/DELFIN/paper_data/mogul_v4_forensik_report.json",
    )
    ap.add_argument("--summary", default="/home/qmchem_max/ComPlat/DELFIN/paper_data/mogul_v4_forensik_summary.md")
    args = ap.parse_args()

    os.environ.setdefault("PYTHONHASHSEED", "0")
    random.seed(args.seed)
    np.random.seed(args.seed)

    pool = Path(args.pool)
    files = sorted([p for p in pool.iterdir() if p.suffix == ".xyz"])
    if len(files) > args.sample:
        files = random.sample(files, args.sample)
    files.sort()

    print(f"Loading v5 lib from {V5_PATH}", flush=True)
    v5 = V5Lookup()
    print(f"Loading COD fragment index", flush=True)
    idx = v3._legacy_load_index()
    print(f"Fragment index size: {len(idx)}", flush=True)

    print(f"Processing {len(files)} files", flush=True)
    t0 = time.time()
    results = []
    for i, f in enumerate(files):
        results.append(run_one_xyz(str(f), v5, idx))
        if (i + 1) % 50 == 0:
            print(f"  ...{i+1}/{len(files)} elapsed={time.time()-t0:.1f}s", flush=True)

    summary = aggregate(results)
    summary["pool"] = str(pool)
    summary["sample_size"] = len(files)
    summary["elapsed_s"] = round(time.time() - t0, 1)
    summary["v5_lib"] = V5_PATH

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(
            {"summary": summary, "per_file": results}, fh, indent=2, default=str
        )

    md = ["# Mogul v4 Forensik Summary", ""]
    md.append(f"- pool: `{pool.name}`")
    md.append(f"- sample size: **{summary['sample_size']}** ({summary['n_files_ok']} ok, {summary['n_files_err']} err)")
    md.append(f"- v5 lib: `{V5_PATH}`")
    md.append(f"- elapsed: {summary['elapsed_s']} s")
    md.append("")
    md.append("## Anomaly Buckets (v3 CCDC mode @ thr=2.5)")
    md.append(f"- **total anomalies**: {summary['total_anomalies']}")
    md.append(f"- X-H bonds: {summary['n_xh']}")
    md.append(f"- D-M-D angles: {summary['n_dmd']}")
    md.append(f"- M-D-X angles: {summary['n_mdx']}")
    md.append(f"- torsion-multimodal (v5 says n_comp >= 2): {summary['n_torsion_multi']}")
    md.append(f"- **v4 would eliminate (severity drops below 2 with v4 lookup)**: {summary['n_v4_eliminates']}")
    md.append("")
    md.append("Eliminated breakdown:")
    for k, v in summary["eliminated_by_cat"].items():
        md.append(f"- {k}: {v}")
    md.append("")
    md.append("## Severity Shift (old v3 MAD-sev vs new v5 Gaussian-sev)")
    for k in ("xh", "dmd", "mdx", "torsion"):
        s = summary[f"severity_change_{k}"]
        if s:
            md.append(
                f"- {k}: n={s['n']}, old_mean={s['old_mean']:.2f} -> new_mean={s['new_mean']:.2f}, "
                f"old_med={s['old_median']:.2f} -> new_med={s['new_median']:.2f}"
            )
        else:
            md.append(f"- {k}: no records")

    Path(args.summary).write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
