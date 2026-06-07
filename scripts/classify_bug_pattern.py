#!/usr/bin/env python3
"""classify_bug_pattern.py — map a structure to known bug-class names.

Bug classes are deterministic boolean tests over geometry + composition.
A structure can match MULTIPLE classes (e.g., both 'cl_cl_1_3_collapse' and
'donor_drift'); we return ALL matches.  KNOWN_PATTERNS is dict-of-dict for
discoverability; the actual checks are concrete predicates so test scaffolds
can fire them on synthetic coords without importing chemistry stacks.

Public API:
  classify_bug(xyz_path) -> List[str]    # ALL matching class names
  classify_atoms(atoms) -> List[str]     # same but takes parsed atoms
  KNOWN_PATTERNS -> Dict[name, {desc, severity, test_fn}]

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations
import itertools
import math
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from diagnostic_common import (  # noqa: E402
    METALS, DONOR_ELEMS, HALOGENS, D8_METALS,
    COVRAD,
    parse_xyz_first_frame, find_metal_idx, coord_sphere, covalent_bonds,
    summarise_structure, dist,
    tetrahedral_score, square_planar_score,
)


# ---------------------------------------------------------------------------
# Bug-class predicates (each takes a 'ctx' dict pre-computed for the frame)
# ---------------------------------------------------------------------------
def _bug_cl_cl_1_3_collapse(ctx: Dict) -> bool:
    """Two halide DONORS through metal closer than 2.0 A.

    1-3 = both bonded to the metal, NOT to each other (collapse through metal).
    """
    P = ctx["positions"]
    syms = ctx["syms"]
    donors = ctx["donor_indices"]
    for i, j in itertools.combinations(donors, 2):
        if syms[i] in HALOGENS and syms[j] in HALOGENS:
            if dist(P[i], P[j]) < 2.0:
                return True
    return False


def _bug_donor_donor_1_3_collapse(ctx: Dict) -> bool:
    """Any two NON-halide donors through the metal < 2.2 A (broad collapse)."""
    P = ctx["positions"]
    syms = ctx["syms"]
    donors = ctx["donor_indices"]
    for i, j in itertools.combinations(donors, 2):
        ei, ej = syms[i], syms[j]
        if ei in HALOGENS or ej in HALOGENS:
            continue
        if dist(P[i], P[j]) < 2.2:
            return True
    return False


def _bug_ch_collapse(ctx: Dict) -> bool:
    """Any C-H bond compressed < 0.9 A."""
    P = ctx["positions"]
    syms = ctx["syms"]
    for i, j in ctx["bonds"]:
        if (syms[i] == "C" and syms[j] == "H") or (syms[i] == "H" and syms[j] == "C"):
            if dist(P[i], P[j]) < 0.9:
                return True
    return False


def _bug_xh_collapse(ctx: Dict) -> bool:
    """Any X-H (X = O,N,S,P) bond compressed < 0.85 A."""
    P = ctx["positions"]
    syms = ctx["syms"]
    for i, j in ctx["bonds"]:
        a, b = syms[i], syms[j]
        if {a, b} & {"H"} and ({a, b} - {"H"}) & {"O", "N", "S", "P"}:
            if dist(P[i], P[j]) < 0.85:
                return True
    return False


def _bug_donor_drift(ctx: Dict) -> bool:
    """A donor sits 2.4-3.0 A from metal — lost first shell but still 'donor'.

    Note: the coord_sphere() in diagnostic_common already excludes most
    out-of-shell donors, so this fires only when the cutoff barely accepts
    them (e.g. carbon-with-2.40 cap, or P with 2.6-2.8 bracket).
    """
    mp = ctx["metal_pos"]
    P = ctx["positions"]
    for d_idx, d in zip(ctx["donor_indices"], ctx["donor_distances"]):
        if 2.4 < d < 3.0:
            return True
    return False


def _bug_md_too_short(ctx: Dict) -> bool:
    """Any M-D distance < 1.5 A (M+donor radius sum is ≥1.7 A for everything)."""
    for d in ctx["donor_distances"]:
        if d < 1.5:
            return True
    return False


def _bug_h_h_clash(ctx: Dict) -> bool:
    """Any non-bonded H-H pair within 1.4 A (severe steric clash)."""
    P = ctx["positions"]
    syms = ctx["syms"]
    bond_set = {(i, j) for i, j in ctx["bonds"]}
    h_idx = [i for i, s in enumerate(syms) if s == "H"]
    for a, b in itertools.combinations(h_idx, 2):
        if (a, b) in bond_set or (b, a) in bond_set:
            continue
        if dist(P[a], P[b]) < 1.4:
            return True
    return False


def _bug_atom_overlap(ctx: Dict) -> bool:
    """Any heavy-atom pair within 0.8 A (collapsed onto each other)."""
    P = ctx["positions"]
    syms = ctx["syms"]
    heavy = [i for i, s in enumerate(syms) if s != "H"]
    for a, b in itertools.combinations(heavy, 2):
        if dist(P[a], P[b]) < 0.8:
            return True
    return False


def _bug_polyhedron_mismatch_d8_t4(ctx: Dict) -> bool:
    """d^8 metal (Ni/Pd/Pt) in T-4 instead of SP-4.

    Compares two CShM-proxies: if tetrahedral score < square-planar score,
    we believe the structure prefers T-4 — which for d^8 is the wrong choice
    geometrically (real complexes pick SP-4 by LFSE).
    """
    if ctx["metal_elem"] not in D8_METALS:
        return False
    donors = ctx["donor_indices"]
    if len(donors) != 4:
        return False
    dp = ctx["donor_positions"]
    t = tetrahedral_score(ctx["metal_pos"], dp)
    s = square_planar_score(ctx["metal_pos"], dp)
    if math.isnan(t) or math.isnan(s):
        return False
    # Decisive: t is meaningfully smaller than s
    return (t + 0.05) < s


def _bug_polyhedron_mismatch_d10_sp4(ctx: Dict) -> bool:
    """d^10 (Zn/Cd/Hg) chooses SP-4 over T-4 (LFSE 0 favours T-4)."""
    if ctx["metal_elem"] not in {"Zn", "Cd", "Hg"}:
        return False
    donors = ctx["donor_indices"]
    if len(donors) != 4:
        return False
    dp = ctx["donor_positions"]
    t = tetrahedral_score(ctx["metal_pos"], dp)
    s = square_planar_score(ctx["metal_pos"], dp)
    if math.isnan(t) or math.isnan(s):
        return False
    return (s + 0.05) < t


def _bug_donor_count_too_low(ctx: Dict) -> bool:
    """Donor count < 2 — coordination collapse."""
    return ctx["cn"] < 2


def _bug_donor_count_excessive(ctx: Dict) -> bool:
    """CN > 9 — almost certainly an over-coordination from the bond rule."""
    return ctx["cn"] > 9


def _bug_metal_isolated(ctx: Dict) -> bool:
    """No bonded donor at all — first-shell collapse OR builder failure."""
    return ctx["cn"] == 0


def _bug_md_too_long_mean(ctx: Dict) -> bool:
    """Mean M-D distance > 2.8 A (typically a face-on / mis-built donor set)."""
    if not ctx["donor_distances"]:
        return False
    mean_md = sum(ctx["donor_distances"]) / len(ctx["donor_distances"])
    return mean_md > 2.8


def _bug_donor_clash_intershell(ctx: Dict) -> bool:
    """Inter-donor distance < 1.8 A but NOT halide-halide (already caught)."""
    P = ctx["positions"]
    syms = ctx["syms"]
    donors = ctx["donor_indices"]
    for i, j in itertools.combinations(donors, 2):
        if syms[i] in HALOGENS and syms[j] in HALOGENS:
            continue
        if dist(P[i], P[j]) < 1.8:
            return True
    return False


def _bug_metal_ligand_axis_linearisation(ctx: Dict) -> bool:
    """In CN >= 3 complexes, all donors collinear (∠max-min < 30°).

    Sign that a donor sub-axis collapsed into a line through the metal.
    """
    P = ctx["positions"]
    mp = ctx["metal_pos"]
    donors = ctx["donor_positions"]
    if len(donors) < 3:
        return False
    # all pairwise angles at metal
    angles = []
    for i, j in itertools.combinations(range(len(donors)), 2):
        v1 = (donors[i][0] - mp[0], donors[i][1] - mp[1], donors[i][2] - mp[2])
        v2 = (donors[j][0] - mp[0], donors[j][1] - mp[1], donors[j][2] - mp[2])
        n1 = math.sqrt(sum(c * c for c in v1)) or 1.0
        n2 = math.sqrt(sum(c * c for c in v2)) or 1.0
        cosang = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)
        cosang = max(-1.0, min(1.0, cosang))
        angles.append(math.degrees(math.acos(cosang)))
    if not angles:
        return False
    return (max(angles) - min(angles)) < 30.0


# ---------------------------------------------------------------------------
# Public dispatch
# ---------------------------------------------------------------------------
KNOWN_PATTERNS: Dict[str, Dict] = {
    "cl_cl_1_3_collapse": {
        "description": "Two halide donors 1-3 through metal < 2 Å",
        "severity": "severe",
        "test": _bug_cl_cl_1_3_collapse,
    },
    "donor_donor_1_3_collapse": {
        "description": "Two non-halide donors 1-3 through metal < 2.2 Å",
        "severity": "severe",
        "test": _bug_donor_donor_1_3_collapse,
    },
    "ch_collapse": {
        "description": "C-H bond compressed < 0.9 Å",
        "severity": "moderate",
        "test": _bug_ch_collapse,
    },
    "xh_collapse": {
        "description": "X-H (X=O/N/S/P) bond compressed < 0.85 Å",
        "severity": "moderate",
        "test": _bug_xh_collapse,
    },
    "donor_drift": {
        "description": "Donor sits 2.4-3.0 Å from metal (out of first shell)",
        "severity": "mild",
        "test": _bug_donor_drift,
    },
    "md_too_short": {
        "description": "Any M-D distance < 1.5 Å",
        "severity": "severe",
        "test": _bug_md_too_short,
    },
    "md_too_long_mean": {
        "description": "Mean M-D > 2.8 Å (face-on / mis-built donor set)",
        "severity": "moderate",
        "test": _bug_md_too_long_mean,
    },
    "h_h_clash": {
        "description": "Non-bonded H-H pair < 1.4 Å",
        "severity": "moderate",
        "test": _bug_h_h_clash,
    },
    "atom_overlap": {
        "description": "Heavy-atom pair < 0.8 Å (collapsed)",
        "severity": "severe",
        "test": _bug_atom_overlap,
    },
    "polyhedron_mismatch_d8_t4": {
        "description": "d^8 metal (Ni/Pd/Pt/Rh/Ir/Au) in T-4 instead of SP-4",
        "severity": "moderate",
        "test": _bug_polyhedron_mismatch_d8_t4,
    },
    "polyhedron_mismatch_d10_sp4": {
        "description": "d^10 (Zn/Cd/Hg) in SP-4 instead of T-4",
        "severity": "moderate",
        "test": _bug_polyhedron_mismatch_d10_sp4,
    },
    "donor_count_too_low": {
        "description": "Donor count < 2 (coordination collapse)",
        "severity": "severe",
        "test": _bug_donor_count_too_low,
    },
    "donor_count_excessive": {
        "description": "Donor count > 9 (over-coordination)",
        "severity": "moderate",
        "test": _bug_donor_count_excessive,
    },
    "metal_isolated": {
        "description": "Metal has no bonded donors",
        "severity": "severe",
        "test": _bug_metal_isolated,
    },
    "donor_clash_intershell": {
        "description": "Non-halide donor pair < 1.8 Å (inter-donor clash)",
        "severity": "severe",
        "test": _bug_donor_clash_intershell,
    },
    "metal_axis_linearisation": {
        "description": "All donor sub-angles collapse into a single axis (range < 30°)",
        "severity": "severe",
        "test": _bug_metal_ligand_axis_linearisation,
    },
}


def _ctx_from_atoms(atoms) -> Optional[Dict]:
    """Pre-compute the context dict used by every predicate."""
    mi = find_metal_idx(atoms)
    if mi is None:
        return None
    metal_elem, metal_pos = atoms[mi]
    donors, _ = coord_sphere(atoms, mi)
    return {
        "syms": [a[0] for a in atoms],
        "positions": [a[1] for a in atoms],
        "metal_idx": mi,
        "metal_elem": metal_elem,
        "metal_pos": metal_pos,
        "donor_indices": [d[0] for d in donors],
        "donor_elems": [d[1] for d in donors],
        "donor_positions": [d[2] for d in donors],
        "donor_distances": [d[3] for d in donors],
        "cn": len(donors),
        "bonds": covalent_bonds(atoms),
    }


def classify_atoms(atoms) -> List[str]:
    """Return ALL matching bug-class names for a parsed atom list."""
    if not atoms:
        return ["parse_failed"]
    ctx = _ctx_from_atoms(atoms)
    if ctx is None:
        return ["no_metal"]
    hits: List[str] = []
    for name, spec in KNOWN_PATTERNS.items():
        try:
            if spec["test"](ctx):
                hits.append(name)
        except Exception:
            # silent skip — predicates must be robust on weird inputs
            continue
    return hits


def classify_bug(xyz_path: str) -> List[str]:
    """Return ALL matching bug-class names for an XYZ file (first frame)."""
    atoms = parse_xyz_first_frame(xyz_path)
    return classify_atoms(atoms)


def severity_label(bug_class: str) -> str:
    """Severity ('severe' / 'moderate' / 'mild') of a known class."""
    spec = KNOWN_PATTERNS.get(bug_class)
    return spec["severity"] if spec else "unknown"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xyz", nargs="*", help="one or more xyz files")
    ap.add_argument("--list", action="store_true",
                    help="list known bug classes and exit")
    args = ap.parse_args()
    if args.list:
        for name, spec in KNOWN_PATTERNS.items():
            print(f"{name:35s} [{spec['severity']:8s}] {spec['description']}")
        return
    if not args.xyz:
        ap.error("no xyz files supplied")
    for p in args.xyz:
        hits = classify_bug(p)
        line = ", ".join(hits) if hits else "(clean)"
        print(f"{p}: {line}")


if __name__ == "__main__":
    main()
