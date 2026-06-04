"""scripts/grace_scaling_plot.py — Task E: scatter + fit on the GRACE
scaling profile.

Reads ``paper_data/grace_scaling_profile.csv`` and produces an ASCII
fit-quality table + scatter (no GUI required).  Also writes a TSV
summary ``paper_data/grace_scaling_fit_summary.txt`` for the SI.

If matplotlib is available, a PNG ``paper_data/grace_scaling_scatter.png``
is also written.
"""
from __future__ import annotations

import csv
import math
import os
import sys
from typing import Any, Dict, List, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
_DATA = os.path.join(_ROOT, "paper_data")


def _read_csv(path: str) -> List[Dict[str, Any]]:
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    return rows


def _fit_powerlaw(N: List[float], T: List[float]) -> Dict[str, Any]:
    pairs = [(n, t) for n, t in zip(N, T) if n > 0 and t > 0]
    if len(pairs) < 5:
        return {"a": 0.0, "b": 0.0, "c": 0.0, "r2": 0.0, "n_used": len(pairs)}
    n_arr = np.array([p[0] for p in pairs], dtype=float)
    t_arr = np.array([p[1] for p in pairs], dtype=float)
    logN = np.log(n_arr)
    logT = np.log(t_arr)
    c, log_b = np.polyfit(logN, logT, 1)
    b = math.exp(log_b)
    pred = log_b + c * logN
    ss_res = float(np.sum((logT - pred) ** 2))
    ss_tot = float(np.sum((logT - np.mean(logT)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {"a": 0.0, "b": float(b), "c": float(c), "r2": float(r2),
            "n_used": int(len(pairs))}


def main():
    in_csv = os.path.join(_DATA, "grace_scaling_profile.csv")
    if not os.path.exists(in_csv):
        sys.exit(f"ERROR: missing {in_csv}")
    rows = _read_csv(in_csv)
    n = len(rows)

    N_atoms = [float(r.get("n_heavy_atoms", 0) or 0) for r in rows]
    CN = [float(r.get("cn", 0) or 0) for r in rows]
    Nrings = [float(r.get("n_rings", 0) or 0) for r in rows]
    T = [float(r.get("wall_time_s", 0) or 0) for r in rows]
    Emit = [float(r.get("n_emitted", 0) or 0) for r in rows]

    # Restrict to runnable (n_emitted >= 1) for cn fits — cn=0 entries
    # are decompose-failures and contribute zero structure.
    N_atoms_run = [na for na, em in zip(N_atoms, Emit) if em >= 1]
    T_run = [t for t, em in zip(T, Emit) if em >= 1]
    CN_run = [c for c, em in zip(CN, Emit) if em >= 1]
    Nrings_run = [r for r, em in zip(Nrings, Emit) if em >= 1]

    fit_N = _fit_powerlaw(N_atoms, T)
    fit_N_run = _fit_powerlaw(N_atoms_run, T_run)
    fit_CN = _fit_powerlaw(CN_run, T_run)
    # Aggregate scaling on (CN^2 * (Nrings+1)) — heuristic complexity proxy.
    proxy = [c * c * (r + 1) for c, r in zip(CN_run, Nrings_run)]
    fit_proxy = _fit_powerlaw(proxy, T_run)

    # Aggregate emit-count scaling
    fit_emit_N = _fit_powerlaw(N_atoms, Emit)

    out_path = os.path.join(_DATA, "grace_scaling_fit_summary.txt")
    with open(out_path, "w") as fh:
        fh.write("# GRACE scaling fits (log-linear power law fits)\n")
        fh.write(f"# n_total = {n}  (runnable >0 sub-fits below)\n\n")

        def _dump(tag, fit):
            fh.write(f"  {tag:32s}: t = {fit['b']:.6g} * N^{fit['c']:.3f}  "
                     f"R^2={fit['r2']:.3f}  n_used={fit['n_used']}\n")

        _dump("t vs n_heavy_atoms (all)", fit_N)
        _dump("t vs n_heavy_atoms (runnable)", fit_N_run)
        _dump("t vs cn (runnable)", fit_CN)
        _dump("t vs CN^2*(nrings+1) (runnable)", fit_proxy)
        _dump("n_emitted vs n_heavy_atoms", fit_emit_N)

        # Per-family rollup
        fh.write("\n# Per-family timing rollup:\n")
        from collections import defaultdict
        by_fam = defaultdict(list)
        for r in rows:
            fam = r.get("family", "?")
            try:
                t = float(r.get("wall_time_s", 0))
            except Exception:
                t = 0
            try:
                em = float(r.get("n_emitted", 0))
            except Exception:
                em = 0
            by_fam[fam].append((t, em))
        for fam, lst in sorted(by_fam.items()):
            ts = [t for t, _ in lst]
            ems = [e for _, e in lst]
            if ts:
                fh.write(f"  {fam:8s} n={len(ts):3d}  "
                         f"t_mean={np.mean(ts):.3f}s  "
                         f"t_max={max(ts):.3f}s  "
                         f"emit_mean={np.mean(ems):.1f}  "
                         f"emit_max={int(max(ems))}\n")

        # ASCII scatter
        fh.write("\n# ASCII scatter (n_heavy_atoms vs wall_time_s)\n")
        T_runnable = [(na, t) for na, t in zip(N_atoms, T) if na > 0 and t > 0]
        if T_runnable:
            T_runnable.sort(key=lambda x: x[0])
            n_bins = 12
            bins: Dict[int, List[float]] = {}
            n_min = min(p[0] for p in T_runnable)
            n_max = max(p[0] for p in T_runnable)
            for na, tt in T_runnable:
                bn = int((na - n_min) / max(n_max - n_min, 1) * (n_bins - 1))
                bins.setdefault(bn, []).append(tt)
            fh.write("    N_atoms_bin    | mean_t(s)  | n\n")
            fh.write("    ----------------------------------\n")
            for bn in range(n_bins):
                lst = bins.get(bn, [])
                if not lst:
                    fh.write(f"    bin{bn:2d}           |   --       |  0\n")
                    continue
                lo = n_min + bn * (n_max - n_min) / max(n_bins - 1, 1)
                fh.write(f"    bin{bn:2d}  N~{lo:5.1f}  |   {np.mean(lst):5.3f}  | {len(lst):3d}\n")

    # Try matplotlib for PNG
    try:
        import matplotlib  # noqa
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        # Left: t vs n_heavy_atoms (log-log)
        N_p = [na for na, t in zip(N_atoms, T) if na > 0 and t > 0]
        T_p = [t for na, t in zip(N_atoms, T) if na > 0 and t > 0]
        axes[0].scatter(N_p, T_p, s=14, alpha=0.6)
        if fit_N["n_used"] >= 5 and N_p:
            xs = np.linspace(min(N_p), max(N_p), 50)
            ys = fit_N["b"] * np.power(xs, fit_N["c"])
            axes[0].plot(xs, ys, "r-",
                          label=f"t={fit_N['b']:.3g}*N^{fit_N['c']:.2f}  R²={fit_N['r2']:.2f}")
            axes[0].legend(loc="upper left", fontsize=8)
        axes[0].set_xscale("log"); axes[0].set_yscale("log")
        axes[0].set_xlabel("n heavy atoms")
        axes[0].set_ylabel("GRACE wall time (s)")
        axes[0].set_title("GRACE scaling: t vs molecule size")
        axes[0].grid(True, which="both", linestyle=":", linewidth=0.5)

        # Right: n_emitted vs n_heavy_atoms
        E_p = [em for na, em in zip(N_atoms, Emit) if na > 0 and em > 0]
        N_e = [na for na, em in zip(N_atoms, Emit) if na > 0 and em > 0]
        if N_e:
            axes[1].scatter(N_e, E_p, s=14, alpha=0.6, color="green")
        axes[1].set_xlabel("n heavy atoms")
        axes[1].set_ylabel("GRACE emitted conformers")
        axes[1].set_title("GRACE emission vs molecule size")
        axes[1].grid(True, which="both", linestyle=":", linewidth=0.5)

        png_path = os.path.join(_DATA, "grace_scaling_scatter.png")
        fig.tight_layout()
        fig.savefig(png_path, dpi=120)
        print(f"# PNG: {png_path}")
    except Exception as exc:
        print(f"# matplotlib unavailable / failed: {exc}")

    print(f"# Summary written: {out_path}")
    with open(out_path) as fh:
        sys.stdout.write(fh.read())


if __name__ == "__main__":
    main()
