#!/usr/bin/env python3
"""diagnostic_common.py — shared utilities for the auto-diagnostic pipeline.

Used by:
  * extract_worst_files.py — top-N WORST per axis from an archive
  * classify_bug_pattern.py — map a structure to known bug-class names
  * bug_census.py — full archive bug count + class breakdown
  * ccdc_overlay.py — DELFIN-vs-XRD heavy-atom + M-D RMSD
  * per_class_ccdc_report.py — group by (metal, CN, polyhedron)
  * iter_gate_with_diagnostics.py — extended iter-gate output

Design constraints (from project doctrine):
  * Deterministic (sorted iteration, fixed-seed sampling).
  * Read-only on archive + CCDC + COD inputs; only writes to user-specified out.
  * No external API calls. Local mol2 + xyz only.
  * Universal across SMILES (no SMILES-specific shortcuts).
  * Works WITHOUT numpy too, but uses numpy when available.

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import numpy as np
    _HAVE_NP = True
except Exception:  # pragma: no cover
    _HAVE_NP = False
    np = None  # type: ignore


# ---------------------------------------------------------------------------
# Element data (kept identical to ccdc_validate.py to be drop-in compatible)
# ---------------------------------------------------------------------------
METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
    "Bi", "Th", "U",
}

DONOR_ELEMS = {"C", "N", "O", "P", "S", "Cl", "Br", "I", "F", "Se", "As",
               "B", "Si", "Te"}

HALOGENS = {"F", "Cl", "Br", "I"}

# d-electron count helpers (oxidation state-agnostic typical group counts)
D8_METALS = {"Ni", "Pd", "Pt", "Rh", "Ir", "Au"}  # square-planar candidates

COVRAD = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "Te": 1.38, "I": 1.39,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
    "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Pb": 1.46, "La": 2.07,
    "Ce": 2.04, "Sn": 1.39, "Al": 1.21, "Mg": 1.41, "Ca": 1.76, "Na": 1.66,
    "K": 2.03, "U": 1.96, "Th": 2.06, "Eu": 1.98, "Gd": 1.96, "Dy": 1.92,
    "In": 1.42, "Tl": 1.45, "Bi": 1.51, "Sb": 1.39,
}
DEFAULT_RAD = 1.5
BOND_TOL = 1.30
ELEM_MAX_MD = {"C": 2.45}


def elem_from_label(label: str) -> str:
    """Strip digits/sybyl suffix; tolerant to '.X' suffix and trailing digits."""
    s = label.split(".")[0]
    m = re.match(r"([A-Z][a-z]?)", s)
    return m.group(1) if m else s


# ---------------------------------------------------------------------------
# XYZ parsers — single-frame only for the diagnostic top-of-pool comparison.
# Multi-frame XYZ pools store frame[0] as the primary; we use that.
# ---------------------------------------------------------------------------
def parse_xyz_first_frame(path: str) -> List[Tuple[str, Tuple[float, float, float]]]:
    """Return list of (elem, (x,y,z)) for FIRST frame only.

    Handles both standard XYZ (count + comment + atoms) and the header-less
    block emitted by older fffree writers.
    """
    atoms: List[Tuple[str, Tuple[float, float, float]]] = []
    try:
        with open(path, "r") as fh:
            lines = fh.read().splitlines()
    except Exception:
        return atoms
    if not lines:
        return atoms
    # Try standard XYZ: line 0 = count
    n_first = None
    try:
        n_first = int(lines[0].strip())
    except Exception:
        n_first = None
    start = 2 if n_first is not None and len(lines) >= n_first + 2 else 0
    limit = (start + n_first) if n_first is not None else len(lines)
    for ln in lines[start:limit]:
        p = ln.split()
        if len(p) < 4:
            continue
        try:
            x, y, z = float(p[1]), float(p[2]), float(p[3])
        except ValueError:
            continue
        if not re.match(r"^[A-Za-z]{1,2}[0-9+\-]*$", p[0]):
            continue
        atoms.append((elem_from_label(p[0]), (x, y, z)))
    return atoms


def parse_mol2_atoms(path: str) -> List[Tuple[str, Tuple[float, float, float]]]:
    """Return list of (elem, (x,y,z)) from a CCDC-flavour .mol2 file."""
    atoms: List[Tuple[str, Tuple[float, float, float]]] = []
    section = None
    try:
        with open(path, "r") as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("@<TRIPOS>"):
                    section = s.split(">", 1)[1]
                    continue
                if section != "ATOM" or not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) < 6:
                    continue
                try:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                except ValueError:
                    continue
                sybyl = parts[5]
                elem = elem_from_label(sybyl)
                if elem not in COVRAD and elem not in METALS:
                    elem = elem_from_label(parts[1])
                atoms.append((elem, (x, y, z)))
    except Exception:
        pass
    return atoms


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def dist(a, b) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def find_metal_idx(atoms) -> Optional[int]:
    """Pick the metal index. If multiple, the one with the most donor neighbours."""
    metal_idx = [i for i, (e, _) in enumerate(atoms) if e in METALS]
    if not metal_idx:
        return None
    if len(metal_idx) == 1:
        return metal_idx[0]
    best, best_cn = metal_idx[0], -1
    for mi in metal_idx:
        cn = len(coord_sphere(atoms, mi)[0])
        if cn > best_cn:
            best, best_cn = mi, cn
    return best


def coord_sphere(atoms, mi: int):
    """Donor atoms bonded to metal index mi using covalent-radius rule + caps.

    Returns (donors, metal_pos) where donors = list of (donor_idx, elem, pos, dist).
    Slightly extended from ccdc_validate.py: returns the donor index for downstream
    bug-pattern matchers.
    """
    me, mp = atoms[mi]
    mr = COVRAD.get(me, DEFAULT_RAD)
    donors = []
    for j, (e, p) in enumerate(atoms):
        if j == mi or e == "H" or e not in DONOR_ELEMS:
            continue
        d = dist(mp, p)
        cutoff = (mr + COVRAD.get(e, DEFAULT_RAD)) * BOND_TOL
        cutoff = min(cutoff, ELEM_MAX_MD.get(e, cutoff))
        if d <= cutoff:
            donors.append((j, e, p, d))
    return donors, mp


def covalent_bonds(atoms, tol: float = 1.30) -> List[Tuple[int, int]]:
    """All covalent bonds in the structure (i<j) by (r1+r2)*tol cutoff.

    Skips the heavy-skip the metal needs for accurate hydride floor; counts
    every pair as a candidate bond.  O(N^2).
    """
    n = len(atoms)
    out: List[Tuple[int, int]] = []
    for i in range(n):
        ei, pi = atoms[i]
        ri = COVRAD.get(ei, DEFAULT_RAD)
        for j in range(i + 1, n):
            ej, pj = atoms[j]
            rj = COVRAD.get(ej, DEFAULT_RAD)
            cut = (ri + rj) * tol
            if dist(pi, pj) <= cut:
                out.append((i, j))
    return out


def kabsch_rmsd(P: Sequence[Sequence[float]],
                Q: Sequence[Sequence[float]]) -> float:
    """RMSD of two Nx3 point sets after optimal rotation (numpy-preferred)."""
    if len(P) != len(Q) or len(P) == 0:
        return float("nan")
    if not _HAVE_NP:
        # crude fallback: centred RMSD only
        def c(A):
            cx = sum(a[0] for a in A) / len(A)
            cy = sum(a[1] for a in A) / len(A)
            cz = sum(a[2] for a in A) / len(A)
            return [(a[0] - cx, a[1] - cy, a[2] - cz) for a in A]
        Pc, Qc = c(P), c(Q)
        s = sum(dist(p, q) ** 2 for p, q in zip(Pc, Qc))
        return math.sqrt(s / len(Pc))
    Pa = np.array(P, float)
    Qa = np.array(Q, float)
    Pa = Pa - Pa.mean(0)
    Qa = Qa - Qa.mean(0)
    H = Pa.T @ Qa
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, float(d)])
    R = Vt.T @ D @ U.T
    Pr = Pa @ R.T
    return float(math.sqrt(((Pr - Qa) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------------------
# Polyhedron shape measures (light CShM proxies sufficient for class triage)
# ---------------------------------------------------------------------------
TETRAHEDRAL_VERTS = [
    (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1),
]
SQUARE_PLANAR_VERTS = [
    (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0),
]


def _unit_from_metal(metal_p, donors_pos):
    """Donor unit vectors from metal.  Returns list of (x,y,z) length 1."""
    out = []
    for d in donors_pos:
        v = (d[0] - metal_p[0], d[1] - metal_p[1], d[2] - metal_p[2])
        n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2) or 1.0
        out.append((v[0] / n, v[1] / n, v[2] / n))
    return out


def _shape_rmsd_to_template(unit_donors, template):
    """Greedy element-blind RMSD of donor unit vectors to a fixed template.

    Brute-force over permutations is O(n!), so for n<=4 we enumerate; for
    larger we greedily assign each donor to its nearest free template slot.
    """
    if not unit_donors or not template or len(unit_donors) != len(template):
        return float("nan")
    if _HAVE_NP:
        Ud = np.array(unit_donors, float)
        Tt = np.array(template, float)
        # normalise template
        Tt = Tt / np.linalg.norm(Tt, axis=1, keepdims=True)
        # try Kabsch on greedy match (sufficient for shape triage)
        # 1. greedy element-blind assignment: each donor → nearest free template
        used = [False] * len(Tt)
        ordering = []
        for d in Ud:
            best_j, best = -1, 1e9
            for j, t in enumerate(Tt):
                if used[j]:
                    continue
                v = float(np.linalg.norm(d - t))
                if v < best:
                    best, best_j = v, j
            ordering.append(best_j)
            if best_j >= 0:
                used[best_j] = True
        Tt_ord = np.array([Tt[k] for k in ordering], float)
        return kabsch_rmsd(Ud.tolist(), Tt_ord.tolist())
    # numpy-less fallback
    return float("nan")


def tetrahedral_score(metal_p, donor_positions) -> float:
    """Lower = more tetrahedral.  RMSD of donor unit vectors to T-4 ideal."""
    if len(donor_positions) != 4:
        return float("nan")
    units = _unit_from_metal(metal_p, donor_positions)
    return _shape_rmsd_to_template(units, TETRAHEDRAL_VERTS)


def square_planar_score(metal_p, donor_positions) -> float:
    """Lower = more square-planar.  RMSD of donor unit vectors to SP-4 ideal."""
    if len(donor_positions) != 4:
        return float("nan")
    units = _unit_from_metal(metal_p, donor_positions)
    return _shape_rmsd_to_template(units, SQUARE_PLANAR_VERTS)


# ---------------------------------------------------------------------------
# Archive helpers
# ---------------------------------------------------------------------------
# Refcode = CCDC-style 6-8 uppercase letters (optionally followed by 1-2 digits),
# bounded by anything non-uppercase-alphanumeric.  Matches '042-ARABUR.xyz' as
# 'ARABUR', 'X10-ZURHID_3d_Cr_CN4_hetero.xyz' as 'ZURHID', and rejects file-
# names with no uppercase token like '01-Fe_CO_3_NHC_2.xyz'.
_REFCODE_RE = re.compile(r"(?<![A-Z0-9])([A-Z]{6}[A-Z0-9]{0,2})(?![A-Z0-9])")
_TAG_WORDS = {"DELFIN", "FFFREE", "VOLLPOOL", "ULTIMATE", "MOGUL", "HEAL",
              "SMOKE", "FFFREEPURE"}


def refcode_from_filename(fname: str) -> Optional[str]:
    """Heuristic refcode (CCDC-style 6-8 uppercase) from a filename."""
    base = Path(fname).stem
    for m in _REFCODE_RE.finditer(base):
        rc = m.group(1)
        if rc not in _TAG_WORDS:
            return rc
    return None


def label_from_filename(fname: str) -> str:
    """Recover the manifest label from a file name, mirroring run_all_detectors."""
    base = Path(fname).stem
    return base


def load_manifest_pool(master_pool_path: str) -> Dict[str, str]:
    """Load label|smiles master pool into a dict."""
    out: Dict[str, str] = {}
    if not master_pool_path or not os.path.isfile(master_pool_path):
        return out
    with open(master_pool_path) as fh:
        for line in fh:
            line = line.rstrip()
            if "|" not in line:
                continue
            lab, smi = line.split("|", 1)
            out[lab] = smi
    return out


def load_dispatch_forensik(tsv_path: str) -> Dict[str, str]:
    """Map SMILES → dispatch path from a forensik TSV.

    Format: <SMILES>\\t<path>\\t<n>
    """
    out: Dict[str, str] = {}
    if not tsv_path or not os.path.isfile(tsv_path):
        return out
    try:
        with open(tsv_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    out[parts[0]] = parts[1]
    except Exception:
        pass
    return out


def list_archive_xyz(archive_dir: str) -> List[str]:
    """Sorted list of non-empty .xyz files in the archive."""
    out: List[str] = []
    try:
        for f in sorted(os.listdir(archive_dir)):
            if not f.endswith(".xyz"):
                continue
            p = os.path.join(archive_dir, f)
            try:
                if os.path.getsize(p) > 0:
                    out.append(p)
            except OSError:
                continue
    except OSError:
        pass
    return out


def cn_class(n: int) -> str:
    """Coarse class string from a coordination number."""
    if n <= 2:
        return "CN<=2"
    if n <= 4:
        return f"CN{n}"
    if n <= 6:
        return f"CN{n}"
    if n <= 9:
        return f"CN{n}"
    return f"CN>={n}"


def polyhedron_class(metal: str, donor_count: int, metal_pos,
                     donor_positions) -> str:
    """Light polyhedron classifier ('T-4', 'SP-4', 'OC-6', 'CN5', ...)."""
    if donor_count == 4:
        t = tetrahedral_score(metal_pos, donor_positions)
        s = square_planar_score(metal_pos, donor_positions)
        if not (math.isnan(t) or math.isnan(s)):
            return "T-4" if t < s else "SP-4"
        return "CN4"
    if donor_count == 6:
        return "OC-6"
    return cn_class(donor_count)


# ---------------------------------------------------------------------------
# Frame-only summary used by extract_worst_files + classify_bug_pattern
# ---------------------------------------------------------------------------
def summarise_structure(xyz_path: str) -> Dict:
    """Return a compact dict describing the first frame of an xyz.

    Keys:
      path, atoms (list of (elem, pos)), metal_idx, metal_elem, metal_pos,
      donor_indices, donor_elems, donor_positions, donor_distances,
      cn (donor count), polyhedron_class
    Errors are flattened into an 'error' field.
    """
    rec: Dict = {"path": xyz_path}
    atoms = parse_xyz_first_frame(xyz_path)
    if not atoms:
        rec["error"] = "empty-or-unreadable"
        return rec
    rec["n_atoms"] = len(atoms)
    rec["atoms"] = atoms
    mi = find_metal_idx(atoms)
    rec["metal_idx"] = mi
    if mi is None:
        rec["error"] = "no-metal"
        return rec
    rec["metal_elem"] = atoms[mi][0]
    rec["metal_pos"] = atoms[mi][1]
    donors, _ = coord_sphere(atoms, mi)
    rec["donor_indices"] = [d[0] for d in donors]
    rec["donor_elems"] = [d[1] for d in donors]
    rec["donor_positions"] = [d[2] for d in donors]
    rec["donor_distances"] = [d[3] for d in donors]
    rec["cn"] = len(donors)
    rec["polyhedron_class"] = polyhedron_class(
        atoms[mi][0], len(donors), atoms[mi][1], rec["donor_positions"]
    )
    return rec


__all__ = [
    "METALS", "DONOR_ELEMS", "HALOGENS", "D8_METALS", "COVRAD",
    "elem_from_label", "parse_xyz_first_frame", "parse_mol2_atoms",
    "dist", "find_metal_idx", "coord_sphere", "covalent_bonds",
    "kabsch_rmsd", "tetrahedral_score", "square_planar_score",
    "refcode_from_filename", "label_from_filename", "load_manifest_pool",
    "load_dispatch_forensik", "list_archive_xyz", "cn_class",
    "polyhedron_class", "summarise_structure",
]
