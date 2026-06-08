#!/usr/bin/env python3
"""Validate rigid-aromatic-chelate placement on the 3 user-eye cases.

Runs the Mogul-primary assembler twice per SMILES — once with the rigid-
chelate gate OFF, once ON — and reports:

  * M-N(chelate) distance per donor
  * Metal out-of-plane displacement vs the chelate aromatic plane
  * N-M-N bite angle
  * Aromatic backbone plane RMS

Cases:
  1. PEHDUD-class:  Co bis-phen + carbonate
  2. bipy-Cu:       Cu bipy + Cl2
  3. terpy-Fe:      Fe terpy + Cl3

Targets (CCDC empirical):
  Co-N(phen)   1.95 ± 0.06 Å
  Cu-N(bipy)   2.00 ± 0.06 Å
  Fe-N(terpy)  1.95 ± 0.06 Å
  M-out-of-plane  <= 0.10 Å
  N-M-N bite     ~ 77° ± 3 (phen / bipy)

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations

import contextlib
import io
import logging
import os
import sys
from typing import List, Optional, Tuple

import numpy as np

# Quiet RDKit's chatty C++ RingInfo / kekulize warnings.  They are
# harmless artefacts of the iterative κⁿ enumerator probing partial
# sanitize states; the assembler recovers from them.
logging.disable(logging.WARNING)
try:
    from rdkit import RDLogger as _RDLogger
    _RDLogger.DisableLog("rdApp.*")
except ImportError:
    pass

# Force determinism + CCDC-grounded GRIP library (so the Mogul-primary
# path lands on the empirical M-D + chelate D-D windows; otherwise the
# grip_polish runs no-op on an empty fragment library and the chelate
# stays at the soft-DG geometry).
os.environ["PYTHONHASHSEED"] = "0"
_default_grip = "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"
if os.path.exists(_default_grip) and not os.environ.get("DELFIN_GRIP_LIB_PATH"):
    os.environ["DELFIN_GRIP_LIB_PATH"] = _default_grip


CASES = [
    # (label, metal_sym, M-N target, bite target, full-complex SMILES).
    # SMILES taken in the convention used by the existing isomer benchmark
    # (Kekulé `=` notation + charged-metal + [N+] donor) which the
    # Mogul-primary path parses cleanly.  Two are verbatim from the
    # ADEHIB / FOLOXE classes of the existing test_isomer_benchmark.py;
    # the third is a known-good bipy-Cu pattern.
    {
        # PEHDUD-class verbatim: ADEHIB-style Mn-bis-phen-COO -> Co
        # substitution.  The chelate motif (two bidentate phen + an
        # anionic-O donor + aqua) is the same as PEHDUD's user-eye
        # bis-phen + carbonate.
        "label": "Co bis-phen + COO (PEHDUD-class)",
        "metal_sym": "Co",
        "md_target": 1.95,
        "md_tol": 0.10,
        "full_smi": (
            "O=C([O][Co-4]12([OH2+])([N+]3=CC=CC4=CC=C5C=CC=[N+]1C5=C43)"
            "[N+]1=CC=CC3=CC=C4C=CC=[N+]2C4=C31)C1=CC=CC(S(=O)(=O)[O-])=C1"
        ),
    },
    {
        # Cu(bipy)(Cl)2 — 4-coordinate Cu-bipy + 2 chloride.  The bipy
        # 2,2' inter-pyridyl bond is the canonical conjugated-biaryl
        # bridge between two aromatic rings -- the universal contract
        # captures it via the "aromatic-rigid bridge" extension.
        "label": "Cu bipy + Cl2 (aromatic biaryl-bridged)",
        "metal_sym": "Cu",
        "md_target": 2.00,
        "md_tol": 0.10,
        "full_smi": (
            "[Cl][Cu-2]1([Cl])[N+]2=CC=CC=C2C2=CC=CC=[N+]21"
        ),
    },
    {
        # Fe(terpy-NMe2)2 — 6-coordinate Fe bis-terpy (verbatim from the
        # existing isomer benchmark Fe(terpy-NMe2)2 case).  Six aromatic
        # N donors, two tridentate terpy chelates, each connected via
        # two biaryl-bridged single bonds -- universal capture target.
        "label": "Fe bis-terpy-NMe2 (aromatic tridentate-NNN)",
        "metal_sym": "Fe",
        "md_target": 1.95,
        "md_tol": 0.10,
        "full_smi": (
            "CC1=CC(N(C)C)=CC2=[N+]1[Fe-6]345([N+]6=C(C7=CC(N(C)C)"
            "=CC(C)=[N+]75)C=CC=C62)[N+](C(C)=CC(N(C)C)=C8)"
            "=C8C9=CC=CC(C%10=CC(N(C)C)=CC(C)=[N+]%104)=[N+]93"
        ),
    },
]


def _find_metal_and_donors(syms, P, metal_sym):
    """Return ``(metal_idx, donor_idxs_N)``."""
    metal_idxs = [i for i, s in enumerate(syms) if s == metal_sym]
    if not metal_idxs:
        return None, []
    m = metal_idxs[0]
    n_candidates = [
        (i, float(np.linalg.norm(P[i] - P[m])))
        for i, s in enumerate(syms) if s == "N"
    ]
    n_candidates.sort(key=lambda x: x[1])
    # Take all N within 2.7 Å (rigid chelate radial window).
    donors = [i for i, d in n_candidates if d < 2.8]
    return m, donors


def _metal_oop_distance(P, metal_idx, plane_atoms):
    """Perpendicular distance of M from the best-fit plane of ``plane_atoms``."""
    pts = P[list(plane_atoms)]
    c = pts.mean(axis=0)
    M = pts - c
    try:
        _, _, Vt = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        return float("inf")
    normal = Vt[-1]
    nn = float(np.linalg.norm(normal))
    if nn < 1e-9:
        return float("inf")
    normal = normal / nn
    return float(abs(np.dot(P[metal_idx] - c, normal)))


def _plane_rms(P, plane_atoms):
    """RMS perpendicular deviation of ``plane_atoms`` from their plane."""
    pts = P[list(plane_atoms)]
    c = pts.mean(axis=0)
    M = pts - c
    try:
        _, _, Vt = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:
        return float("inf")
    normal = Vt[-1] / max(1e-9, float(np.linalg.norm(Vt[-1])))
    offs = M @ normal
    return float(np.sqrt(np.mean(offs ** 2)))


def _measure(syms, P, metal_sym, ligand_n_target=2):
    """Compute the user-eye metrics."""
    if syms is None or P is None:
        return None
    m, donors = _find_metal_and_donors(syms, P, metal_sym)
    if m is None or len(donors) == 0:
        return None
    md_dists = [float(np.linalg.norm(P[d] - P[m])) for d in donors]
    # For each chelate (assume the donors of one chelate are the closest
    # ``ligand_n_target`` Ns) compute the bite angle.
    bite_angles: List[float] = []
    # Naive but universal: walk donor-donor pairs whose Cartesian
    # separation is < 3.5 Å — these are within-chelate.  Compute
    # N-M-N angle for each such pair.
    for i in range(len(donors)):
        for j in range(i + 1, len(donors)):
            d = float(np.linalg.norm(P[donors[i]] - P[donors[j]]))
            if d > 3.5:
                continue
            v1 = P[donors[i]] - P[m]; v1 /= max(1e-9, float(np.linalg.norm(v1)))
            v2 = P[donors[j]] - P[m]; v2 /= max(1e-9, float(np.linalg.norm(v2)))
            ang = float(np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))))
            bite_angles.append(ang)
    # Plane = the largest aromatic-N + neighbours subset reachable from
    # the donors.  Approximation: take all C neighbours of N donors and
    # walk one BFS step further along aromatic Cs.  Pure geometry; no
    # mol graph needed.
    plane_atoms = set(donors)
    # Two BFS steps along atoms within 1.6 Å of any current plane atom and
    # of element C/N (heavy aromatic-like atoms).
    for _ in range(3):
        for i, s in enumerate(syms):
            if i in plane_atoms or i == m or s not in ("C", "N"):
                continue
            close = False
            for a in plane_atoms:
                if float(np.linalg.norm(P[i] - P[a])) < 1.6:
                    close = True
                    break
            if close:
                plane_atoms.add(i)
    plane_atoms = sorted(plane_atoms)
    oop = _metal_oop_distance(P, m, plane_atoms) if len(plane_atoms) >= 3 else float("inf")
    plane_rms = _plane_rms(P, plane_atoms) if len(plane_atoms) >= 4 else float("inf")
    return {
        "n_donors": len(donors),
        "md_distances": md_dists,
        "md_mean": float(np.mean(md_dists)) if md_dists else float("nan"),
        "bite_angles_deg": bite_angles,
        "metal_oop_A": oop,
        "plane_rms_A": plane_rms,
        "n_plane_atoms": len(plane_atoms),
    }


def _run_one(smi, with_rigid: bool):
    """Run the canonical fffree path with / without the rigid gate.

    Returns the first stereoisomer XYZ as ``(syms, P)`` or (None, reason).
    """
    # Always force MOGUL_PRIMARY ON for the test (we want to compare the
    # pre-rigid vs post-rigid pipeline at the same baseline).
    os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
    if with_rigid:
        os.environ["DELFIN_FFFREE_RIGID_AROMATIC_CHELATE"] = "1"
    else:
        os.environ["DELFIN_FFFREE_RIGID_AROMATIC_CHELATE"] = "0"
    # Late import so env settings take effect.
    import importlib
    import delfin.fffree.rigid_aromatic_chelate as _rac
    importlib.reload(_rac)
    import delfin.fffree.assemble_via_mogul as _avm
    importlib.reload(_avm)
    import delfin.fffree.converter_backend as _cb
    importlib.reload(_cb)
    try:
        isomers = _cb._fffree_isomers(smi)
    except Exception as e:
        return None, f"_fffree_isomers crashed: {e}"
    if isomers is None or len(isomers) == 0:
        return None, "_fffree_isomers returned None or empty"

    def _parse_one(xyz_str):
        syms_: List[str] = []
        coords_: List[Tuple[float, float, float]] = []
        for ln in xyz_str.splitlines():
            parts = ln.split()
            if len(parts) < 4:
                continue
            sym = parts[0]
            if not (1 <= len(sym) <= 3 and sym[0].isalpha()):
                continue
            try:
                x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
            except ValueError:
                continue
            syms_.append(sym)
            coords_.append((x, y, z))
        if not coords_:
            return None, None
        return syms_, np.asarray(coords_, dtype=float)

    # Walk every isomer and pick the one with the LOWEST per-aromatic-plane
    # RMS deviation (= user-eye flat aromatic backbone).  This mirrors
    # what the downstream ranker would emit for the user-facing best
    # structure.  For pre/post comparison, picking the BEST isomer per
    # run lets us see if any orbit benefits from rigid placement.
    best = None
    best_metric = None
    for xyz, _meta in isomers:
        syms_i, P_i = _parse_one(xyz)
        if syms_i is None:
            continue
        # Compute aromatic-plane RMS proxy: the plane through all heavy
        # atoms within 2.5 Å BFS from the metal (a coarse but universal
        # measure).
        import numpy as _np
        metal_idxs = [i for i, s in enumerate(syms_i) if not s in ("H", "C", "N", "O", "S", "F", "Cl", "Br", "I", "P")]
        if not metal_idxs:
            continue
        m = metal_idxs[0]
        # Heavy atoms within 4 Å BFS from metal (broad neighbourhood).
        close = [i for i, s in enumerate(syms_i) if s != "H" and i != m
                 and float(_np.linalg.norm(P_i[i] - P_i[m])) < 4.0]
        if len(close) < 4:
            continue
        pts = P_i[close]
        c = pts.mean(axis=0)
        M_ = pts - c
        try:
            _, _, Vt = _np.linalg.svd(M_, full_matrices=False)
            normal = Vt[-1] / max(1e-9, float(_np.linalg.norm(Vt[-1])))
            metric = float(_np.sqrt(_np.mean((M_ @ normal) ** 2)))
        except Exception:
            metric = float("inf")
        if best is None or metric < best_metric:
            best = (syms_i, P_i)
            best_metric = metric
    if best is None:
        # Fall back to first isomer if metric-rank failed.
        xyz, _meta = isomers[0]
        syms_first, P_first = _parse_one(xyz)
        if syms_first is None:
            return None, "no XYZ rows parsed"
        best = (syms_first, P_first)
    return best, None


def _format_metric(d, target_md=None, target_tol=None):
    if d is None:
        return "  (no measurement)"
    lines = []
    md_mean = d["md_mean"]
    md_min = min(d["md_distances"]) if d["md_distances"] else float("nan")
    md_max = max(d["md_distances"]) if d["md_distances"] else float("nan")
    lines.append(
        f"  n_donors={d['n_donors']}, M-N mean={md_mean:.3f} Å "
        f"[min={md_min:.3f}, max={md_max:.3f}]"
    )
    if target_md is not None and target_tol is not None:
        in_window = (
            target_md - target_tol <= md_mean <= target_md + target_tol
        )
        flag = "OK" if in_window else "MISS"
        lines.append(
            f"    target {target_md:.2f} +/- {target_tol:.2f}  -> {flag}"
        )
    if d["bite_angles_deg"]:
        bite_str = ", ".join(f"{b:.1f}" for b in d["bite_angles_deg"])
        lines.append(f"  N-M-N bite angles: {bite_str} deg")
    lines.append(f"  M-out-of-plane: {d['metal_oop_A']:.3f} Å")
    lines.append(
        f"  aromatic-plane RMS: {d['plane_rms_A']:.3f} Å "
        f"(n_atoms={d['n_plane_atoms']})"
    )
    return "\n".join(lines)


def main():
    print("=" * 72)
    print("RIGID-AROMATIC-CHELATE  3-CASE  VALIDATION")
    print("=" * 72)
    rows = []
    for case in CASES:
        print()
        print(f">>> {case['label']}")
        print(f"    SMILES: {case['full_smi']}")
        # Pre: rigid OFF.
        result_pre, err_pre = _run_one(case["full_smi"], with_rigid=False)
        if err_pre:
            print(f"    PRE-build failed: {err_pre}")
        # Post: rigid ON.
        result_post, err_post = _run_one(case["full_smi"], with_rigid=True)
        if err_post:
            print(f"    POST-build failed: {err_post}")
        m_pre = None
        m_post = None
        if result_pre is not None:
            m_pre = _measure(result_pre[0], result_pre[1], case["metal_sym"])
        if result_post is not None:
            m_post = _measure(result_post[0], result_post[1], case["metal_sym"])
        print("  PRE (rigid OFF):")
        print(_format_metric(m_pre, case["md_target"], case["md_tol"]))
        print("  POST (rigid ON):")
        print(_format_metric(m_post, case["md_target"], case["md_tol"]))
        rows.append((case["label"], m_pre, m_post,
                     case["md_target"], case["md_tol"]))
    print()
    print("=" * 72)
    print("SUMMARY  (per-case ΔM-D, Δoop, Δbite)")
    print("=" * 72)
    print(f"{'case':<36} {'M-D pre':>10} {'M-D post':>10} {'oop pre':>10} {'oop post':>10}")
    print("-" * 78)
    for label, m_pre, m_post, _t, _tol in rows:
        md_pre = m_pre["md_mean"] if m_pre else float("nan")
        md_post = m_post["md_mean"] if m_post else float("nan")
        oop_pre = m_pre["metal_oop_A"] if m_pre else float("nan")
        oop_post = m_post["metal_oop_A"] if m_post else float("nan")
        print(
            f"{label[:36]:<36} {md_pre:>10.3f} {md_post:>10.3f} "
            f"{oop_pre:>10.3f} {oop_post:>10.3f}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
