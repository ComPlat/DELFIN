#!/usr/bin/env python3
"""Mission G1' — XRD-RMSD comparator (heavy-atom Kabsch + bond/angle RMSD).

For each (refcode, builder) pair we compute:
  - heavy_rmsd  : Kabsch-aligned heavy-atom (no H) RMSD vs CCDC XRD-resolved
                  cleaned structure. Element-matched greedy correspondence.
  - bond_rmsd   : RMSD over heavy-heavy covalent bonds (mean over shared bonds)
  - angle_rmsd  : RMSD over heavy-heavy-heavy angles (mean over shared angles)
  - n_match     : number of element-matched heavy atom pairs used

The CCDC ground-truth comes from
``/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_cleaned_xyz/``
(counter-ions stripped, disorder resolved, single component).

Each builder's structures come from a pool archive directory containing
XYZ files whose names embed the refcode (e.g. ``D2-ABUJAJ_4d_Rh_CN5.xyz``,
``X10-ALAWUF_3d_Ni_CN3_hetero.xyz``, ``042-ARABUR.xyz``).

USAGE
    python scripts/xrd_rmsd_comparator.py \
        --pool /path/to/2792332-aromatic-symmetry-VOLLPOOL \
        --label delfin_aromatic \
        --out paper_data/xrd_rmsd_delfin_aromatic.jsonl

Per-pool emits one JSONL line per (refcode, best-isomer-within-pool):
    {"refcode": "ABEFEV", "label": "delfin_aromatic",
     "heavy_rmsd": 0.187, "bond_rmsd": 0.014, "angle_rmsd": 4.7,
     "n_match": 22, "n_emitted": 4, "metal": "Ag", "cn": 3, "block": "d"}

We report the BEST emission per refcode (lowest heavy_rmsd), reflecting
"if DELFIN/UFF emits multiple isomers/conformers, take the one closest
to XRD" — that's the standard structure-prediction metric.

Determinism: PYTHONHASHSEED=0; greedy correspondence is lex-sorted.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Constants — covalent radii for bond detection (no FF, geometry-only)
# ---------------------------------------------------------------------------
_COV: Dict[str, float] = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "Te": 1.38, "I": 1.39,
    # Metals
    "Li": 1.28, "Be": 0.96, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "K": 2.03,
    "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50,
    "Fe": 1.42, "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90,
    "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42,
    "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62,
    "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Pb": 1.46, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Sm": 1.98,
    "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89,
    "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Th": 2.06, "U": 1.96, "Sn": 1.39,
    "Sb": 1.39, "In": 1.42, "Ga": 1.22, "Ge": 1.20, "Bi": 1.48,
}
_BOND_TOL = 1.30  # standard tolerance for covalent radii sum

METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu", "Th", "U", "Sn", "Pb", "Bi", "In", "Ga", "Ge",
}


# ---------------------------------------------------------------------------
# Filename -> refcode patterns (D2-, X10-, NNN-)
# ---------------------------------------------------------------------------
_REFCODE_PATTERNS = [
    re.compile(r"^D2-([A-Z][A-Z0-9]{4,8})_"),
    re.compile(r"^X10-([A-Z][A-Z0-9]{4,8})_"),
    re.compile(r"^[0-9]{2,4}-([A-Z][A-Z0-9]{4,8})\.xyz$"),
    # Bare REFCODE.xyz (ccdc_cleaned_xyz style)
    re.compile(r"^([A-Z][A-Z0-9]{4,8})\.xyz$"),
]


def extract_refcode(filename: str) -> Optional[str]:
    """Return CCDC refcode from filename or None."""
    bn = os.path.basename(filename)
    for pat in _REFCODE_PATTERNS:
        m = pat.match(bn)
        if m:
            return m.group(1)
    return None


# ---------------------------------------------------------------------------
# XYZ parsing
# ---------------------------------------------------------------------------
def parse_xyz_all_frames(text: str) -> List[Tuple[List[str], np.ndarray]]:
    """Parse a multi-frame XYZ file, returning list of (symbols, positions).

    Standard XYZ format:
        <N>
        <comment>
        <symbol> <x> <y> <z>
        ... N rows ...
        <N>
        <comment>
        ...
    """
    frames: List[Tuple[List[str], np.ndarray]] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        # Find next atom-count header
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue
        try:
            n_atoms = int(stripped)
        except ValueError:
            i += 1
            continue
        if n_atoms <= 0 or i + 1 + n_atoms >= len(lines) + 1:
            i += 1
            continue
        # Skip comment line
        atom_start = i + 2
        atom_end = atom_start + n_atoms
        syms: List[str] = []
        pos: List[List[float]] = []
        for j in range(atom_start, min(atom_end, len(lines))):
            parts = lines[j].split()
            if len(parts) < 4:
                continue
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError:
                continue
            m = re.match(r"([A-Z][a-z]?)", parts[0])
            if not m:
                continue
            syms.append(m.group(1))
            pos.append([x, y, z])
        if len(syms) == n_atoms:
            frames.append((syms, np.asarray(pos, dtype=float)))
        i = atom_end
    return frames


def parse_xyz(text: str) -> Tuple[List[str], np.ndarray]:
    """Return (symbols, positions) of FIRST frame only (legacy single-frame)."""
    frames = parse_xyz_all_frames(text)
    if frames:
        return frames[0]
    return [], np.zeros((0, 3))


def read_xyz_file(path: Path) -> Tuple[List[str], np.ndarray]:
    return parse_xyz(path.read_text())


def read_xyz_file_all_frames(path: Path) -> List[Tuple[List[str], np.ndarray]]:
    return parse_xyz_all_frames(path.read_text())


# ---------------------------------------------------------------------------
# Bond / angle topology (heavy atoms only)
# ---------------------------------------------------------------------------
def heavy_atom_indices(syms: List[str]) -> List[int]:
    return [i for i, s in enumerate(syms) if s != "H"]


def detect_bonds_heavy(syms: List[str], P: np.ndarray,
                      heavy_idx: List[int]) -> List[Tuple[int, int]]:
    """Return bonds as pairs of indices INTO heavy_idx (i.e. heavy-only)."""
    bonds: List[Tuple[int, int]] = []
    for a in range(len(heavy_idx)):
        ia = heavy_idx[a]
        sa = syms[ia]
        ra = _COV.get(sa, 1.5)
        for b in range(a + 1, len(heavy_idx)):
            ib = heavy_idx[b]
            sb = syms[ib]
            rb = _COV.get(sb, 1.5)
            d = float(np.linalg.norm(P[ia] - P[ib]))
            cutoff = (ra + rb) * _BOND_TOL
            # cap M-C (suppress hapto over-bond)
            if (sa in METALS and sb == "C") or (sb in METALS and sa == "C"):
                cutoff = min(cutoff, 2.45)
            if d < cutoff:
                bonds.append((a, b))
    return bonds


def angles_from_bonds(bonds: List[Tuple[int, int]],
                       n_heavy: int) -> List[Tuple[int, int, int]]:
    """Enumerate angle triplets (i, j, k) with j the center vertex."""
    adj: Dict[int, List[int]] = {i: [] for i in range(n_heavy)}
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)
    out: List[Tuple[int, int, int]] = []
    for j, nbrs in adj.items():
        nbrs_sorted = sorted(nbrs)
        for ii in range(len(nbrs_sorted)):
            for kk in range(ii + 1, len(nbrs_sorted)):
                out.append((nbrs_sorted[ii], j, nbrs_sorted[kk]))
    return out


# ---------------------------------------------------------------------------
# Kabsch alignment + element-matched greedy correspondence
# ---------------------------------------------------------------------------
def _greedy_element_match(syms_a: List[str], Pa: np.ndarray,
                          syms_b: List[str], Pb: np.ndarray
                          ) -> List[Tuple[int, int]]:
    """Greedy element-typed nearest-neighbour matching (post-centroid).

    Deterministic: iterate sources in lex-sort order; pick lowest-distance
    same-element target not yet used.
    """
    used_b = [False] * len(syms_b)
    pairs: List[Tuple[int, int]] = []
    order = sorted(range(len(syms_a)), key=lambda k: (syms_a[k], k))
    for i in order:
        si = syms_a[i]
        best_j = -1
        best_d = float("inf")
        for j in range(len(syms_b)):
            if used_b[j] or syms_b[j] != si:
                continue
            d = float(np.linalg.norm(Pa[i] - Pb[j]))
            if d < best_d:
                best_d = d
                best_j = j
        if best_j >= 0:
            used_b[best_j] = True
            pairs.append((i, best_j))
    return pairs


def kabsch_align(P_src: np.ndarray, P_tgt: np.ndarray
                 ) -> Tuple[np.ndarray, float]:
    """Return (rotated P_src centered to match P_tgt centered, rmsd)."""
    cs = P_src.mean(0)
    ct = P_tgt.mean(0)
    Ps = P_src - cs
    Pt = P_tgt - ct
    H = Ps.T @ Pt
    U, S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    Pr = Ps @ R.T
    rmsd = float(math.sqrt(((Pr - Pt) ** 2).sum(1).mean()))
    return Pr, rmsd


def heavy_rmsd_paired(syms_a: List[str], Pa: np.ndarray,
                       syms_b: List[str], Pb: np.ndarray
                       ) -> Tuple[float, List[Tuple[int, int]]]:
    """Return (heavy_rmsd, pairs of heavy-atom indices into original arrays)."""
    keep_a = heavy_atom_indices(syms_a)
    keep_b = heavy_atom_indices(syms_b)
    if len(keep_a) < 3 or len(keep_b) < 3:
        return float("nan"), []
    Pa_h = Pa[keep_a]
    Pb_h = Pb[keep_b]
    syms_ah = [syms_a[i] for i in keep_a]
    syms_bh = [syms_b[i] for i in keep_b]

    # Initial centering for matching
    Pa_h_c = Pa_h - Pa_h.mean(0)
    Pb_h_c = Pb_h - Pb_h.mean(0)
    pairs_h = _greedy_element_match(syms_ah, Pa_h_c, syms_bh, Pb_h_c)
    if len(pairs_h) < 3:
        return float("nan"), []
    Pa_m = np.array([Pa_h_c[i] for i, _ in pairs_h])
    Pb_m = np.array([Pb_h_c[j] for _, j in pairs_h])
    _, rmsd = kabsch_align(Pa_m, Pb_m)
    # Translate pairs back to original-array indices
    pairs_orig = [(keep_a[i], keep_b[j]) for i, j in pairs_h]
    return rmsd, pairs_orig


# ---------------------------------------------------------------------------
# Bond / angle RMSDs over shared topology
# ---------------------------------------------------------------------------
def bond_angle_rmsd(syms_a: List[str], Pa: np.ndarray,
                    syms_b: List[str], Pb: np.ndarray,
                    pairs_orig: List[Tuple[int, int]]
                    ) -> Tuple[float, float, int, int]:
    """Compute bond_rmsd and angle_rmsd over shared bonds/angles.

    pairs_orig maps (idx_in_A, idx_in_B). We build heavy-only graphs in BOTH
    A and B, then compute bond/angle for each (i,j) pair if BOTH endpoints
    are in pairs_orig AND both pools have that bond. Same for angles
    over shared triplets.

    Returns (bond_rmsd_angstrom, angle_rmsd_degrees, n_bonds, n_angles).
    """
    if not pairs_orig:
        return float("nan"), float("nan"), 0, 0

    a2b = dict(pairs_orig)
    b2a = {b: a for a, b in pairs_orig}

    # Heavy-only graphs in each pool
    heavy_a = heavy_atom_indices(syms_a)
    heavy_b = heavy_atom_indices(syms_b)
    a_to_local = {a: k for k, a in enumerate(heavy_a)}
    b_to_local = {b: k for k, b in enumerate(heavy_b)}
    bonds_a_local = detect_bonds_heavy(syms_a, Pa, heavy_a)
    bonds_b_local = detect_bonds_heavy(syms_b, Pb, heavy_b)
    bonds_a = {(heavy_a[i], heavy_a[j]) for i, j in bonds_a_local}
    bonds_b = {(heavy_b[i], heavy_b[j]) for i, j in bonds_b_local}

    # Shared bonds: bond (ai,aj) in A is shared iff a2b[ai],a2b[aj] is a bond in B
    bond_diffs: List[float] = []
    for (ai, aj) in bonds_a:
        if ai not in a2b or aj not in a2b:
            continue
        bi, bj = a2b[ai], a2b[aj]
        if (bi, bj) in bonds_b or (bj, bi) in bonds_b:
            la = float(np.linalg.norm(Pa[ai] - Pa[aj]))
            lb = float(np.linalg.norm(Pb[bi] - Pb[bj]))
            bond_diffs.append(la - lb)
    if bond_diffs:
        bond_rmsd = float(math.sqrt(sum(d * d for d in bond_diffs)
                                    / len(bond_diffs)))
    else:
        bond_rmsd = float("nan")

    # Angles
    n_heavy_a = len(heavy_a)
    n_heavy_b = len(heavy_b)
    angles_a = angles_from_bonds(bonds_a_local, n_heavy_a)
    angles_b = set(angles_from_bonds(bonds_b_local, n_heavy_b))
    angle_diffs: List[float] = []
    for (i, j, k) in angles_a:
        # heavy_a-local i,j,k -> original atom indices in A
        ai = heavy_a[i]
        aj = heavy_a[j]
        ak = heavy_a[k]
        if ai not in a2b or aj not in a2b or ak not in a2b:
            continue
        bi_orig = a2b[ai]
        bj_orig = a2b[aj]
        bk_orig = a2b[ak]
        if (bi_orig not in b_to_local or bj_orig not in b_to_local
                or bk_orig not in b_to_local):
            continue
        bi = b_to_local[bi_orig]
        bj = b_to_local[bj_orig]
        bk = b_to_local[bk_orig]
        triplet_b = tuple(sorted([bi, bk])) + ()  # placeholder
        triplet_b = (min(bi, bk), bj, max(bi, bk))
        if triplet_b not in angles_b:
            continue
        # Compute angles in degrees
        va = Pa[ai] - Pa[aj]
        vb = Pa[ak] - Pa[aj]
        cos_a = float(np.dot(va, vb) / (np.linalg.norm(va)
                                          * np.linalg.norm(vb) + 1e-12))
        cos_a = max(-1.0, min(1.0, cos_a))
        ang_a = math.degrees(math.acos(cos_a))
        wa = Pb[bi_orig] - Pb[bj_orig]
        wb = Pb[bk_orig] - Pb[bj_orig]
        cos_b = float(np.dot(wa, wb) / (np.linalg.norm(wa)
                                          * np.linalg.norm(wb) + 1e-12))
        cos_b = max(-1.0, min(1.0, cos_b))
        ang_b = math.degrees(math.acos(cos_b))
        angle_diffs.append(ang_a - ang_b)
    if angle_diffs:
        angle_rmsd = float(math.sqrt(sum(d * d for d in angle_diffs)
                                       / len(angle_diffs)))
    else:
        angle_rmsd = float("nan")

    return bond_rmsd, angle_rmsd, len(bond_diffs), len(angle_diffs)


# ---------------------------------------------------------------------------
# Pool scanning + per-refcode best-emission selection
# ---------------------------------------------------------------------------
def group_pool_by_refcode(pool_dir: Path) -> Dict[str, List[Path]]:
    out: Dict[str, List[Path]] = {}
    if not pool_dir.is_dir():
        return out
    for fp in sorted(pool_dir.iterdir()):
        if not fp.name.endswith(".xyz"):
            continue
        rc = extract_refcode(fp.name)
        if rc:
            out.setdefault(rc, []).append(fp)
    return out


def score_one(refcode: str, candidate_paths: List[Path],
              ccdc_syms: List[str], ccdc_P: np.ndarray,
              label: str, classification: Dict[str, str]
              ) -> Optional[dict]:
    """Pick best candidate (lowest heavy_rmsd) across files AND frames.

    Each candidate file may contain multiple isomer/conformer frames.
    We score EVERY frame and pick the globally best (min heavy_rmsd)
    over all (file, frame) pairs.
    """
    best: Optional[dict] = None
    total_frames_seen = 0
    for cpath in sorted(candidate_paths):
        try:
            frames = read_xyz_file_all_frames(cpath)
        except Exception:
            continue
        if not frames:
            continue
        for frame_idx, (syms_c, Pc) in enumerate(frames):
            total_frames_seen += 1
            if len(syms_c) == 0 or len(Pc) == 0:
                continue
            try:
                heavy_r, pairs = heavy_rmsd_paired(syms_c, Pc,
                                                    ccdc_syms, ccdc_P)
            except Exception:
                continue
            if not math.isfinite(heavy_r):
                continue
            try:
                bond_r, angle_r, n_b, n_ang = bond_angle_rmsd(
                    syms_c, Pc, ccdc_syms, ccdc_P, pairs)
            except Exception:
                bond_r, angle_r, n_b, n_ang = (float("nan"),
                                                float("nan"), 0, 0)
            rec = {
                "refcode": refcode,
                "label": label,
                "candidate_file": str(cpath.name),
                "frame_idx": frame_idx,
                "heavy_rmsd": heavy_r,
                "bond_rmsd": bond_r,
                "angle_rmsd": angle_r,
                "n_match": len(pairs),
                "n_bonds_shared": n_b,
                "n_angles_shared": n_ang,
                "n_atoms_candidate": len(syms_c),
                "n_atoms_ccdc": len(ccdc_syms),
            }
            rec.update(classification)
            if best is None or rec["heavy_rmsd"] < best["heavy_rmsd"]:
                best = rec
    if best is not None:
        best["n_emitted"] = total_frames_seen
        best["n_files"] = len(candidate_paths)
    return best


# ---------------------------------------------------------------------------
# Classification from CCDC index (cn, block, metal_element, geom_label)
# ---------------------------------------------------------------------------
def load_ccdc_index(path: Path) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            rc = rec.get("refcode")
            if rc:
                out[rc] = rec
    return out


def classify_record(rec: dict) -> Dict[str, str]:
    return {
        "cn": rec.get("cn"),
        "block": rec.get("block"),
        "metal": rec.get("metal_element"),
        "n_metals": rec.get("n_metals"),
        "geom": rec.get("geom_label"),
        "is_organometallic": rec.get("is_organometallic"),
        "n_atoms_raw": rec.get("n_atoms_raw"),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True, help="Pool dir w/ XYZ files")
    ap.add_argument("--label", required=True, help="Builder label e.g. delfin_aromatic, uff, delfin_grip")
    ap.add_argument(
        "--ccdc-xyz",
        default="/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_cleaned_xyz",
        help="Dir of CCDC cleaned XRD-resolved XYZ files (one per refcode)")
    ap.add_argument(
        "--ccdc-index",
        default="/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_tmc_index_cleaned.jsonl",
        help="CCDC cleaned index JSONL (for cn/block/metal classification)")
    ap.add_argument("--out", required=True, help="Output JSONL")
    ap.add_argument("--max-refcodes", type=int, default=0,
                    help="Optional cap (0 = all)")
    args = ap.parse_args()

    pool_dir = Path(args.pool)
    ccdc_xyz_dir = Path(args.ccdc_xyz)
    ccdc_index_path = Path(args.ccdc_index)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not pool_dir.is_dir():
        print(f"ERR: pool dir not found: {pool_dir}", file=sys.stderr)
        sys.exit(2)
    if not ccdc_xyz_dir.is_dir():
        print(f"ERR: ccdc xyz dir not found: {ccdc_xyz_dir}", file=sys.stderr)
        sys.exit(2)

    pool_by_rc = group_pool_by_refcode(pool_dir)
    ccdc_idx = load_ccdc_index(ccdc_index_path)
    print(f"[xrd-rmsd] pool refcodes: {len(pool_by_rc)}", file=sys.stderr)
    print(f"[xrd-rmsd] ccdc index records: {len(ccdc_idx)}", file=sys.stderr)

    # Available CCDC refcodes on disk
    ccdc_files = {fp.stem: fp for fp in ccdc_xyz_dir.iterdir()
                  if fp.name.endswith(".xyz")}
    print(f"[xrd-rmsd] ccdc xyz files on disk: {len(ccdc_files)}",
          file=sys.stderr)

    overlap = sorted(set(pool_by_rc) & set(ccdc_files))
    print(f"[xrd-rmsd] overlap refcodes: {len(overlap)}", file=sys.stderr)
    if args.max_refcodes > 0:
        overlap = overlap[: args.max_refcodes]

    n_written = 0
    n_failed = 0
    with out_path.open("w") as outfh:
        for rc in overlap:
            ccdc_path = ccdc_files[rc]
            try:
                ccdc_syms, ccdc_P = read_xyz_file(ccdc_path)
            except Exception as exc:
                print(f"[xrd-rmsd] skip {rc} (ccdc read err): {exc}",
                      file=sys.stderr)
                n_failed += 1
                continue
            classification = classify_record(ccdc_idx.get(rc, {}))
            rec = score_one(rc, pool_by_rc[rc],
                            ccdc_syms, ccdc_P,
                            args.label, classification)
            if rec is None:
                n_failed += 1
                continue
            outfh.write(json.dumps(rec) + "\n")
            n_written += 1
    print(f"[xrd-rmsd] wrote {n_written} records (failed: {n_failed}) -> "
          f"{out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
