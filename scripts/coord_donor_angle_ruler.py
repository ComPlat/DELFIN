#!/usr/bin/env python3
"""coord_donor_angle_ruler.py — donor-internal M-D-X coordination-angle RULER.

Closes a metric blindspot: the coordination ruler measured only the D-M-D angle
AT THE METAL (180 deg "correct" for a linear CN2 centre) and was BLIND to the
DONOR-INTERNAL M-D-X angle.  Concrete eye-validated miss (YULFAO = EtHg-Se-C6H5):
the Hg centre is correctly linear (C-Hg-Se 180 deg) but the Se donor was built
LINEAR too (Hg-Se-C(phenyl) = 180 deg) when a 2-coordinate chalcogen MUST be bent
~96-104 deg (two lone pairs).  The old ruler scored that structure clean.

This ruler measures, for every coordinated BENT-CAPABLE donor D (a chalcogen
O/S/Se/Te, or a single-bonded pyramidal pnictogen N/P/As/Sb) that carries a
substituent X, the M-D-X angle and flags it when it is far from the donor's VSEPR
ideal.  A 2-coordinate chalcogen sitting near 180 deg (linear) is therefore
flagged.  Measurement only — no mutation.

CCDC-calibrated SANITY (do NOT flag genuinely-linear donors):
  * sp nitrile N (N#C), azo/diazo/azide-terminal N, aromatic in-ring N
  * terminal oxo / double-bonded sp chalcogen (M=O, =S, isothiocyanate S=C=N)
  * a donor with NO substituent (terminal halide / bare chalcogenide) — no X.
These are recognised by graph + bond order + the same aromatic-ring detector the
runtime corrector uses, so a real linear M-D-X is NOT a false positive.

Deterministic, graph/geometry-only, open-source (no CCDC data shipped).  Reuses
the parsing / adjacency / aromatic helpers from
``delfin._coord_angle_corrector`` (single source of truth for the coordination
geometry machinery).

Usage:
    python scripts/coord_donor_angle_ruler.py STRUCT.xyz [STRUCT2.xyz ...]
    python scripts/coord_donor_angle_ruler.py --json STRUCT.xyz
    cat STRUCT.xyz | python scripts/coord_donor_angle_ruler.py -
"""
from __future__ import annotations

import json
import sys
from typing import Dict, List, Optional

import numpy as np

from delfin._coord_angle_corrector import (
    _angle_deg,
    _build_geometric_adjacency,
    _COV_RADII,
    _detect_aromatic_atoms,
    _is_metal_sym,
    _parse_xyz,
)

# VSEPR ideal M-D-X angle (deg) for a bent donor that retains lone pairs.
# CCDC-sane targets: H2O 104.5 / R2O ~111; H2Se 91 / R2Se ~96-98; R2S ~98-100;
# R2Te ~95; R3N / R3P ~107 / phosphine ~100; arsine ~96.  Single source with
# the construction-side _DONOR_BEND_DEG in assemble_complex.py.
_DONOR_IDEAL_DEG: Dict[str, float] = {
    "O": 109.0, "S": 100.0, "Se": 98.0, "Te": 95.0,
    "N": 107.0, "P": 100.0, "As": 96.0, "Sb": 95.0,
}
_CHALCOGENS = frozenset(("O", "S", "Se", "Te"))
_PNICTOGENS = frozenset(("N", "P", "As", "Sb"))

# Tolerance (deg): a bent donor whose M-D-X deviates more than this from its
# VSEPR ideal is flagged.  Wide enough to absorb crystal scatter + the natural
# spread of substituent identity (ether 111 vs water 104.5) but tight enough that
# a LINEAR (180) chalcogen — ~80 deg off the ~100 ideal — is always caught.
_TOL_DEG = 30.0
# A bent donor whose M-D-X is within this of 180 deg is treated as "linearised"
# and ALWAYS flagged regardless of tolerance (the YULFAO failure mode).
_LINEAR_FLAG_DEG = 150.0


def _bond_order_to(nbrs, bond_d, syms, d_idx: int, x_idx: int) -> float:
    """Heuristic bond multiplicity from the D-X distance vs covalent-radii sum.

    Conservative: only a CLEARLY short bond counts as multiple.  A normal polar
    single bond such as Se-C(aryl) ~1.81 A or S-C ~1.78 A sits at ~0.90-0.93 of
    the (generous) covalent-radii sum, so a 0.93 cutoff would mis-call it double
    and wrongly exempt a genuine bent donor.  We require < 0.86 rsum for double
    and < 0.78 rsum for triple, which still catches real C=O / C#N / M=O while
    leaving polar single bonds as order 1."""
    d = bond_d.get((d_idx, x_idx))
    if d is None:
        return 1.0
    rsum = _COV_RADII.get(syms[d_idx], 1.5) + _COV_RADII.get(syms[x_idx], 1.5)
    if rsum <= 0:
        return 1.0
    if d < 0.78 * rsum:
        return 3.0
    if d < 0.86 * rsum:
        return 2.0
    return 1.0


def _x_is_sp_center(nbrs, bond_d, syms, d_idx: int, x_idx: int) -> bool:
    """True if the donor's substituent X is an sp centre that makes the donor
    axis genuinely LINEAR: X has exactly two heavy neighbours and the bond from
    X to its FAR neighbour (not the donor) is a triple bond (nitrile-type
    N#C-D, isocyanide, terminal C#X).  Graph + distance only."""
    x_heavy = [k for k in nbrs[x_idx]
               if k != d_idx and syms[k] != "H" and not _is_metal_sym(syms[k])]
    if len(x_heavy) != 1:
        return False
    far = x_heavy[0]
    return _bond_order_to(nbrs, bond_d, syms, x_idx, far) >= 3.0


def measure_donor_angles(syms: List[str], pts: np.ndarray) -> Dict:
    """Return {'flagged': [...], 'all': [...], 'n_donors': int, 'n_flagged': int}.

    Each entry: M/D/X indices+symbols, observed M-D-X angle, VSEPR ideal,
    deviation, whether linearised, and the reason string.  ``all`` covers every
    bent-capable coordinated donor measured; ``flagged`` is the violating subset.
    """
    n = len(syms)
    nbrs, bond_d = _build_geometric_adjacency(syms, pts)
    aromatic = _detect_aromatic_atoms(syms, pts, nbrs)
    metals = [i for i in range(n) if _is_metal_sym(syms[i])]
    measured: List[Dict] = []
    for m_idx in metals:
        donors = [j for j in nbrs[m_idx]
                  if not _is_metal_sym(syms[j]) and syms[j] != "H"]
        for d_idx in donors:
            sym = syms[d_idx]
            ideal = _DONOR_IDEAL_DEG.get(sym)
            if ideal is None:
                continue                       # not a bent-capable element
            # substituents of D other than the metal (heavy or H)
            subs = [k for k in nbrs[d_idx]
                    if k != m_idx and not _is_metal_sym(syms[k])]
            heavy_subs = [k for k in subs if syms[k] != "H"]
            if not subs:
                continue                       # bare chalcogenide/halide -> no X, not bent
            # --- CCDC-sane LINEAR exclusions (do not flag genuinely-linear) ---
            # aromatic in-ring donor: lone pair is the donor axis (pyridine etc.)
            if d_idx in aromatic:
                continue
            # multiply-bonded terminal donor (oxo M=O, =S, imido/diazo N): a
            # clearly-short double/triple bond at the donor => sp/sp2 linear-to-
            # axis (lone pairs perpendicular, donor axis = the pi bond).
            multibonded = any(
                _bond_order_to(nbrs, bond_d, syms, d_idx, k) >= 2.0
                for k in heavy_subs
            )
            if multibonded:
                continue
            # nitrile-type linear donor: the substituent is an sp centre with a
            # terminal triple bond (N#C-D, isocyanide) => the M-D-X axis is
            # genuinely linear and must NOT be bent.
            if any(_x_is_sp_center(nbrs, bond_d, syms, d_idx, k)
                   for k in heavy_subs):
                continue
            # a pnictogen needs to be PYRAMIDAL (>=1 substituent, single-bonded);
            # a single-substituent (2-coord) or two-substituent (3-coord) bent N/P
            # is fine.  >=3 substituents (4-coord ammonium/phosphonium) -> skip
            # (no lone pair, geometry fixed tetrahedral, not a bend defect).
            if len(subs) >= 3:
                continue
            # measure M-D-X for each substituent X; report the WORST (max dev).
            worst = None
            for x_idx in subs:
                ang = _angle_deg(pts[d_idx], pts[m_idx], pts[x_idx])
                if ang is None:
                    continue
                dev = abs(ang - ideal)
                if worst is None or dev > worst[1]:
                    worst = (x_idx, dev, ang)
            if worst is None:
                continue
            x_idx, dev, ang = worst
            linearised = ang >= _LINEAR_FLAG_DEG
            flagged = linearised or dev > _TOL_DEG
            reason = ""
            if linearised:
                reason = (f"{sym} donor is LINEARISED (M-D-X={ang:.1f} deg, "
                          f"should be ~{ideal:.0f}); 2-coord chalcogen/pnictogen "
                          f"must be bent (lone pairs)")
            elif flagged:
                reason = (f"M-D-X={ang:.1f} deg deviates {dev:.1f} from VSEPR "
                          f"ideal {ideal:.0f} (tol {_TOL_DEG:.0f})")
            measured.append({
                "M_idx": m_idx, "D_idx": d_idx, "X_idx": x_idx,
                "M_sym": syms[m_idx], "D_sym": sym, "X_sym": syms[x_idx],
                "n_substituents": len(subs),
                "observed_angle": round(ang, 2),
                "ideal_angle": ideal,
                "deviation": round(dev, 2),
                "linearised": linearised,
                "flagged": flagged,
                "reason": reason,
            })
    flagged = [m for m in measured if m["flagged"]]
    return {
        "all": measured,
        "flagged": flagged,
        "n_donors": len(measured),
        "n_flagged": len(flagged),
    }


def measure_file(path: str) -> Dict:
    if path == "-":
        xyz = sys.stdin.read()
    else:
        with open(path, "r") as fh:
            xyz = fh.read()
    syms, pts, _ = _parse_xyz(xyz)
    res = measure_donor_angles(syms, pts)
    res["path"] = path
    return res


def _print_human(res: Dict) -> None:
    path = res.get("path", "?")
    n = res["n_donors"]
    nf = res["n_flagged"]
    status = "CLEAN" if nf == 0 else f"FLAGGED ({nf}/{n})"
    print(f"{path}: {status}  [{n} bent-capable donor(s) measured]")
    for m in res["all"]:
        mark = "  !!" if m["flagged"] else "    "
        print(f"{mark} M{m['M_idx']}({m['M_sym']})-"
              f"D{m['D_idx']}({m['D_sym']})-X{m['X_idx']}({m['X_sym']}): "
              f"{m['observed_angle']:.1f} deg (ideal {m['ideal_angle']:.0f}, "
              f"dev {m['deviation']:.1f})"
              + (f"  -> {m['reason']}" if m["flagged"] else ""))


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    as_json = "--json" in argv
    paths = [a for a in argv if a != "--json"]
    if not paths:
        print(__doc__)
        return 2
    results = [measure_file(p) for p in paths]
    if as_json:
        print(json.dumps(results, indent=2))
    else:
        for r in results:
            _print_human(r)
    total_flagged = sum(r["n_flagged"] for r in results)
    return 1 if total_flagged else 0


if __name__ == "__main__":
    raise SystemExit(main())
