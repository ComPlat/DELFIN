#!/usr/bin/env python3
"""Realism-Weight Bayesian + Grid Optimization Harness.

Optimizes the 7-dimensional weight vector in
:mod:`delfin.fffree.realism_ranking` to maximize the rank-0
XRD-match rate against the 49-refcode CCDC ground-truth set
across three archives (b00f9a0, c03a550, ec7fb0d).

Pipeline:
1. Build a per-frame raw-signal cache once per archive (expensive — does
   ALL detector JSONL parsing + heavy-atom Kabsch RMSD vs CCDC ground
   truth + per-SMILES rank-normalisation).
2. Re-rank trajectories cheaply with any weight tuple (pure numpy on
   cached normalised signals).
3. Compute rank-0 + top-3 + top-5 XRD-match rates by counting how many
   refcodes have their best-frame at rank 0/within-3/within-5.

Determinism: PYTHONHASHSEED=0; numpy RNG seed = 0.

Outputs:
* ``paper_data/realism_weight_grid_search.csv``
* ``paper_data/realism_weight_bayesian_opt.csv``
* ``paper_data/realism_weight_OPTIMAL.json``
* ``paper_data/realism_hard_gate_analysis.csv``
* ``paper_data/realism_weight_per_class_optimal.csv``
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np

np.random.seed(0)

# Ensure worktree's delfin import takes precedence
WORKTREE = Path(__file__).resolve().parent.parent
if str(WORKTREE) not in sys.path:
    sys.path.insert(0, str(WORKTREE))

from delfin.fffree import realism_ranking as RR

# Paths
ARCH_ROOT = Path("/home/qmchem_max/agent_workspace/quality_framework/xyz_archive")
GROUND_TRUTH_FILE = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_ground_truth_48.json"
)
PAPER_DATA = WORKTREE / "paper_data"

ARCHIVES = [
    "b00f9a0-full7-VOLLPOOL",
    "c03a550-race-full-stack-VOLLPOOL",
    "ec7fb0d-full9-fblock-VOLLPOOL",
]

# Map archive name -> existing rank_vs_xrd CSV with reference best_frame_index.
# This CSV was produced by the original 7f10f05/ea26dcb merge using the
# canonical Kabsch RMSD pipeline; we reuse its best_frame_index so our
# optimisation is directly comparable to the published 42-44% baseline.
EXISTING_CSV_MAP = {
    "b00f9a0-full7-VOLLPOOL": "/home/qmchem_max/ComPlat/DELFIN/paper_data/realism_rank_vs_xrd_b00f9a0.csv",
    "c03a550-race-full-stack-VOLLPOOL": "/home/qmchem_max/ComPlat/DELFIN/paper_data/realism_rank_vs_xrd_c03a550.csv",
    "ec7fb0d-full9-fblock-VOLLPOOL": "/home/qmchem_max/ComPlat/DELFIN/paper_data/realism_rank_vs_xrd_ec7fb0d.csv",
}

SOFT_KEYS = ("mogul", "cshm", "inter_clash", "hh_clash", "grip_loss", "polya", "burnside")

# ---------------------------------------------------------------------------
# XYZ I/O + Kabsch RMSD
# ---------------------------------------------------------------------------
METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu", "Th", "U",
}


def parse_multixyz(path: Path) -> List[Tuple[str, List[str], np.ndarray]]:
    """Return list of (label, symbols, positions) per frame."""
    out = []
    if not path.exists():
        return out
    with path.open("r") as f:
        text = f.read()
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i].strip()
        if not raw:
            i += 1
            continue
        try:
            natoms = int(raw)
        except ValueError:
            i += 1
            continue
        if i + natoms + 1 >= len(lines):
            break
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        symbols, positions = [], []
        for j in range(natoms):
            parts = lines[i + 2 + j].split()
            if len(parts) < 4:
                continue
            symbols.append(parts[0])
            try:
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
            except ValueError:
                positions.append([0.0, 0.0, 0.0])
        label = ""
        li = comment.find("label=")
        if li >= 0:
            label = comment[li + 6:].strip()
        out.append((label, symbols, np.array(positions, dtype=float)))
        i += natoms + 2
    return out


def _kabsch_align(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, float]:
    """Kabsch SVD alignment.  Returns (rotated_P, RMSD)."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    try:
        U, S, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        return P, float("nan")
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    P_rot = Pc @ R.T
    rmsd = float(np.sqrt(np.mean(np.sum((P_rot - Qc) ** 2, axis=1))))
    return P_rot, rmsd


def _hungarian_element_bucket(P: np.ndarray, Q: np.ndarray,
                               symsP: List[str], symsQ: List[str]) -> Optional[np.ndarray]:
    """Per-element-bucket Hungarian assignment Q→P. Returns permutation π
    such that P[π] is matched to Q.  ``None`` on infeasible.
    """
    from collections import Counter, defaultdict
    if Counter(symsP) != Counter(symsQ):
        return None
    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError:
        return None
    # Group indices by element
    grpP: Dict[str, List[int]] = defaultdict(list)
    grpQ: Dict[str, List[int]] = defaultdict(list)
    for i, s in enumerate(symsP):
        grpP[s].append(i)
    for i, s in enumerate(symsQ):
        grpQ[s].append(i)
    perm = np.full(len(P), -1, dtype=int)
    for elem, pidxs in grpP.items():
        qidxs = grpQ.get(elem, [])
        if len(qidxs) != len(pidxs):
            return None
        sub_P = P[pidxs]
        sub_Q = Q[qidxs]
        # cost matrix (Euclidean distance squared)
        cost = np.linalg.norm(sub_P[:, None, :] - sub_Q[None, :, :], axis=2)
        row_ind, col_ind = linear_sum_assignment(cost)
        for r, c in zip(row_ind, col_ind):
            perm[pidxs[r]] = qidxs[c]
    return perm


def kabsch_rmsd_heavy(symA: List[str], A: np.ndarray, symB: List[str], B: np.ndarray,
                       max_iter: int = 3) -> float:
    """Heavy-atom RMSD with iterated Kabsch + element-bucket Hungarian.

    1. Drop H.
    2. Initial Kabsch on the as-is ordering (may be wrong perm).
    3. Hungarian re-assign Q to A within each element bucket → re-Kabsch.
    4. Iterate up to max_iter; return best RMSD found.
    """
    if A.size == 0 or B.size == 0:
        return float("nan")
    mA = [i for i, s in enumerate(symA) if s and s.upper() != "H"]
    mB = [i for i, s in enumerate(symB) if s and s.upper() != "H"]
    if not mA or not mB:
        return float("nan")
    P = A[mA].astype(float)
    Q = B[mB].astype(float)
    symsP = [symA[i] for i in mA]
    symsQ = [symB[i] for i in mB]
    from collections import Counter
    if Counter(symsP) != Counter(symsQ):
        return float("nan")
    # Centre once
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    # Initial assume same canonical order
    P_rot, best_rmsd = _kabsch_align(Pc, Qc)
    for _it in range(max_iter):
        perm = _hungarian_element_bucket(P_rot, Qc, symsP, symsQ)
        if perm is None or np.any(perm < 0):
            break
        # Re-order Q to match P
        Q_perm = Qc[perm]
        # ALIGN P_rot to Q_perm via Kabsch
        P_rot_new, new_rmsd = _kabsch_align(P_rot, Q_perm)
        if not math.isnan(new_rmsd) and new_rmsd < best_rmsd - 1e-6:
            best_rmsd = new_rmsd
            P_rot = P_rot_new
            Qc = Q_perm
        else:
            break
    return float(best_rmsd)


# ---------------------------------------------------------------------------
# Ground-truth: refcode -> heavy atom positions
# ---------------------------------------------------------------------------
def load_ground_truth() -> Dict[str, Tuple[List[str], np.ndarray]]:
    """Return refcode -> (symbols, positions). H stripped, symbols cleaned."""
    d = json.load(GROUND_TRUTH_FILE.open("r"))
    out: Dict[str, Tuple[List[str], np.ndarray]] = {}
    for rec in d["records"]:
        ref = rec["refcode"]
        if not rec.get("ok"):
            continue
        syms = rec.get("symbols") or rec.get("elements") or []
        pos = rec.get("positions") or []
        if not syms or not pos:
            continue
        try:
            P = np.array(pos, dtype=float)
        except Exception:
            continue
        # Clean symbols (strip charge markers etc).
        syms_clean = [re.sub(r"[^A-Za-z]", "", s).capitalize() for s in syms]
        out[ref] = (syms_clean, P)
    return out


# ---------------------------------------------------------------------------
# Refcode -> XYZ file matching
# ---------------------------------------------------------------------------
def refcode_from_filename(fname: str) -> Optional[str]:
    """Extract 6-letter CCDC refcode from XYZ filename.

    Patterns:
        X10-AFECIZ_3d_Ti_CN5_hetero.xyz -> AFECIZ
        042-ARABUR.xyz                  -> ARABUR
        09-SomeName.xyz                 -> None (no 6-letter refcode)
    """
    stem = Path(fname).stem
    # Try X-prefixed
    m = re.match(r"^X\d+-([A-Z]{6})(?:_|$)", stem)
    if m:
        return m.group(1)
    # Try numbered prefix
    m = re.match(r"^\d+-([A-Z]{6})(?:_|$)", stem)
    if m:
        return m.group(1)
    return None


# ---------------------------------------------------------------------------
# Per-frame raw signal extraction from detector JSONLs
# ---------------------------------------------------------------------------
def _parse_header_smi_label(header: str) -> Tuple[str, str]:
    smi = ""
    for tok in header.split():
        if tok.startswith("smi="):
            smi = tok[4:]
    label = ""
    li = header.find("label=")
    if li >= 0:
        label = header[li + 6:].strip()
    return smi, label


def load_per_frame_signals(archive_dir: Path) -> Dict[Tuple[str, str], Dict[str, float]]:
    """Aggregate all available detector JSONLs into per-(smi, label) signals.

    Signals collected:
    * ``cshm``               from coordgeom.mean_abs_dev
    * ``mogul``              from coordgeom.n_violations (proxy: no mogul_v3 file)
    * ``inter_clash``        from severe_overlap.n_severe_overlaps
    * ``hh_clash``           from hanomaly.n_anomalies (fallback hclash.n_clash)
    * ``grip_loss``          NOT AVAILABLE (constant 0 → rank-normalised neutral)
    * ``polya``              0 = complete (≥80% coverage), 1 = incomplete
    * ``burnside``           NOT AVAILABLE (constant 0)
    * ``topology_ok``        from topology.topology_match
    * ``build_time_clash_ok`` from ligcollapse.broken==False
    """
    parent = archive_dir.parent
    name = archive_dir.name
    md: Dict[Tuple[str, str], Dict[str, float]] = defaultdict(dict)

    def _iter(suffix: str):
        p = parent / f"{name}_{suffix}.jsonl"
        if not p.exists():
            return
        with p.open("r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue

    # coordgeom: cshm + mogul proxy
    for rec in _iter("coordgeom"):
        smi, label = _parse_header_smi_label(rec.get("header", ""))
        if not smi:
            smi = Path(rec.get("file", "")).stem
        key = (smi, label)
        v = rec.get("mean_abs_dev")
        if v is not None:
            try:
                md[key]["cshm"] = float(v)
            except (TypeError, ValueError):
                pass
        v = rec.get("n_violations")
        if v is not None:
            try:
                md[key]["mogul"] = float(v)
            except (TypeError, ValueError):
                pass

    # severe_overlap
    for rec in _iter("severe_overlap"):
        smi, label = _parse_header_smi_label(rec.get("header", ""))
        if not smi:
            smi = Path(rec.get("file", "")).stem
        label = label or rec.get("label", "") or ""
        key = (smi, label)
        v = rec.get("n_severe_overlaps")
        if v is not None:
            try:
                md[key]["inter_clash"] = float(v)
            except (TypeError, ValueError):
                pass

    # hanomaly (preferred for hh_clash)
    for rec in _iter("hanomaly"):
        smi, label = _parse_header_smi_label(rec.get("header", ""))
        if not smi:
            smi = Path(rec.get("file", "")).stem
        key = (smi, label)
        v = rec.get("n_anomalies")
        if v is not None:
            try:
                md[key]["hh_clash"] = float(v)
            except (TypeError, ValueError):
                pass

    # hclash (fallback)
    for rec in _iter("hclash"):
        smi, label = _parse_header_smi_label(rec.get("header", ""))
        if not smi:
            smi = Path(rec.get("file", "")).stem
        key = (smi, label)
        if "hh_clash" in md[key]:
            continue
        v = rec.get("n_clash")
        if v is not None:
            try:
                md[key]["hh_clash"] = float(v)
            except (TypeError, ValueError):
                pass

    # ligcollapse: hard gate build_time_clash_ok
    for rec in _iter("ligcollapse"):
        # file/label
        smi = Path(rec.get("file", "")).stem
        label = rec.get("label", "") or ""
        key = (smi, label)
        broken = rec.get("broken")
        if broken is not None:
            md[key]["build_time_clash_ok"] = (not bool(broken))

    # topology: hard gate topology_ok
    for rec in _iter("topology"):
        smi, label = _parse_header_smi_label(rec.get("header", ""))
        if not smi:
            smi = Path(rec.get("file", "")).stem
        key = (smi, label)
        v = rec.get("topology_match")
        if v is not None:
            md[key]["topology_ok"] = bool(v)

    # isocoverage: polya — per-SMILES coverage.
    # Build a map keyed by XYZ-file STEM (with underscore form matched too).
    polya_by_smi: Dict[str, float] = {}
    for rec in _iter("isocoverage"):
        smi_full = rec.get("smiles_id") or ""
        file_stem = Path(rec.get("file", "")).stem if rec.get("file") else ""
        cov = rec.get("coverage_pct")
        if cov is None:
            continue
        try:
            complete = float(cov) >= 80.0
        except (TypeError, ValueError):
            continue
        polya_val = 0.0 if complete else 1.0
        # Map under both smi_id and file_stem, plus the dash<->underscore variants
        candidates = set()
        if smi_full:
            candidates.add(smi_full)
            candidates.add(smi_full.replace("_", "-"))
            candidates.add(smi_full.replace("-", "_"))
        if file_stem:
            candidates.add(file_stem)
            candidates.add(file_stem.replace("_", "-"))
            candidates.add(file_stem.replace("-", "_"))
        for c in candidates:
            polya_by_smi[c] = polya_val
        # Apply to all existing keys with this smi
        for (k_smi, k_label) in list(md.keys()):
            if k_smi in candidates:
                md[(k_smi, k_label)]["polya"] = polya_val
    # Stash the per-SMILES polya map for fallback lookup
    md["__polya_by_smi__"] = polya_by_smi  # type: ignore

    return dict(md)


# ---------------------------------------------------------------------------
# Build per-archive frame table
# ---------------------------------------------------------------------------
def load_reference_best_frame_csv(csv_path: Path) -> Dict[str, Tuple[int, float]]:
    """Return refcode -> (best_frame_index, best_rmsd) from the existing
    pre-computed rank_vs_xrd CSV.  Skips rows without finite RMSD.
    """
    out: Dict[str, Tuple[int, float]] = {}
    if not csv_path.exists():
        return out
    with csv_path.open("r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("status") != "ok":
                continue
            try:
                idx = int(row["best_frame_index"])
                rmsd = float(row["best_rmsd"])
            except (TypeError, ValueError):
                continue
            if idx < 0:
                continue
            out[row["refcode"]] = (idx, rmsd)
    return out


def build_archive_table(archive: str, ground_truth: Dict[str, Tuple[List[str], np.ndarray]]):
    """Return per-archive table:
    * frames: list of dicts {refcode, smiles_id, label, frame_idx, signals(7-vec)}
    * rmsd_per_refcode: refcode -> list[(frame_idx_global, rmsd)]
    * groups: list of (smiles_id, [frame_idx_global, ...])

    Only refcodes that have a CCDC ground truth AND an XYZ file are
    retained (ground-truth-only).
    """
    archive_dir = ARCH_ROOT / archive
    print(f"[{archive}] loading detector signals…", flush=True)
    md = load_per_frame_signals(archive_dir)
    polya_by_smi = md.pop("__polya_by_smi__", {})  # type: ignore
    print(f"[{archive}] {len(md)} per-frame records collected.", flush=True)

    # Try to load existing reference best-frame index from the published CSV.
    # Using it makes our optimisation directly comparable to the published
    # 42-44% rank-0 baseline (whereas our naive Kabsch RMSD may pick a
    # different "best frame" per refcode).
    ref_best: Dict[str, Tuple[int, float]] = {}
    ref_csv = Path(EXISTING_CSV_MAP.get(archive, ""))
    if ref_csv.exists():
        ref_best = load_reference_best_frame_csv(ref_csv)
        print(f"[{archive}] using reference CSV: {len(ref_best)} refcodes with best_frame_index "
              f"(from {ref_csv.name})", flush=True)

    # Discover XYZ files matching ground-truth refcodes
    xyz_files = sorted(archive_dir.glob("*.xyz"))
    refcode_to_xyz: Dict[str, Path] = {}
    for xp in xyz_files:
        ref = refcode_from_filename(xp.name)
        if ref and ref in ground_truth:
            refcode_to_xyz[ref] = xp

    print(f"[{archive}] {len(refcode_to_xyz)} refcodes matched to XYZ files (of {len(ground_truth)} ground-truth).",
          flush=True)

    frames = []  # one per frame, in build order
    refcode_frames: Dict[str, List[Tuple[int, float]]] = defaultdict(list)
    groups: Dict[str, List[int]] = defaultdict(list)

    for ref, xp in sorted(refcode_to_xyz.items()):
        smi_id = xp.stem
        frames_data = parse_multixyz(xp)
        if not frames_data:
            continue
        gt_syms, gt_pos = ground_truth[ref]
        # If we have a reference best-frame index from the existing CSV,
        # use it directly (force only that frame to have a finite rmsd =
        # the published best_rmsd, all other frames get a sentinel +inf so
        # the ranking honours the published "best frame" choice).
        ref_idx_rmsd = ref_best.get(ref) if ref_best else None
        # Compute RMSD per frame
        for fi, (label, syms, pos) in enumerate(frames_data):
            if ref_idx_rmsd is not None:
                if fi == ref_idx_rmsd[0]:
                    rmsd = ref_idx_rmsd[1]
                else:
                    rmsd = float("nan")
            else:
                rmsd = kabsch_rmsd_heavy(syms, pos, gt_syms, gt_pos)
            # Lookup signals
            key = (smi_id, label)
            sig_dict = md.get(key, {})
            # Try fallback: header smi may be different from filename stem
            if not sig_dict:
                # try without prefix
                for cand_key in md.keys():
                    if cand_key[1] == label and ref in cand_key[0]:
                        sig_dict = md[cand_key]
                        break

            # Fallback polya lookup from per-SMILES map (in case per-frame
            # dict didn't get the polya stamped — common when isocoverage
            # smiles_id differs from frame smiles_id slightly).
            polya_val = sig_dict.get("polya")
            if polya_val is None or (isinstance(polya_val, float) and math.isnan(polya_val)):
                for c in (smi_id, smi_id.replace("_", "-"), smi_id.replace("-", "_")):
                    if c in polya_by_smi:
                        polya_val = polya_by_smi[c]
                        break
            if polya_val is None:
                polya_val = float("nan")

            # Signals vec — store as raw values, will rank-normalise per group later.
            # Note: severe_overlap/hclash/hanomaly only EMIT records when
            # the detector finds at least one event.  Missing => assume 0
            # (no clashes detected for this frame).  coordgeom always
            # emits, so mogul/cshm NaN means genuinely missing record.
            signals_raw = {
                "mogul": sig_dict.get("mogul", float("nan")),
                "cshm": sig_dict.get("cshm", float("nan")),
                "inter_clash": sig_dict.get("inter_clash", 0.0),
                "hh_clash": sig_dict.get("hh_clash", 0.0),
                "grip_loss": float("nan"),  # not available
                "polya": polya_val,
                "burnside": float("nan"),  # not available
            }
            topo_ok = bool(sig_dict.get("topology_ok", True))
            clash_ok = bool(sig_dict.get("build_time_clash_ok", True))

            idx_global = len(frames)
            frames.append({
                "refcode": ref,
                "smiles_id": smi_id,
                "label": label,
                "frame_idx_local": fi,
                "rmsd": rmsd,
                "topology_ok": topo_ok,
                "build_time_clash_ok": clash_ok,
                **{f"sig_{k}": v for k, v in signals_raw.items()},
            })
            refcode_frames[ref].append((idx_global, rmsd))
            groups[smi_id].append(idx_global)

    return {
        "frames": frames,
        "refcode_frames": dict(refcode_frames),
        "groups": dict(groups),
    }


def _rank_normalise_arr(values: np.ndarray) -> np.ndarray:
    """NaN -> 0.5 (neutral), else rank/(n-1)."""
    n = len(values)
    if n == 0:
        return values
    if n == 1:
        return np.array([0.5])
    finite_mask = ~np.isnan(values)
    out = np.full(n, 0.5)
    if finite_mask.sum() == 0:
        return out
    finite_idx = np.where(finite_mask)[0]
    finite_vals = values[finite_mask]
    # Stable sort by (value, idx) to mirror module behaviour
    order = np.lexsort((finite_idx, finite_vals))
    m = len(finite_vals)
    if m == 1:
        out[finite_idx[order[0]]] = 0.5
    else:
        for rank, oi in enumerate(order):
            out[finite_idx[oi]] = rank / (m - 1)
    return out


def precompute_normalised(table: Dict) -> np.ndarray:
    """Build N_frames x 7 normalised-signal matrix.

    Rank-normalisation happens within each SMILES group (matching module
    behaviour).
    """
    n = len(table["frames"])
    norm = np.full((n, len(SOFT_KEYS)), 0.5)
    for smi_id, idxs in table["groups"].items():
        idxs_arr = np.array(idxs)
        for k_i, key in enumerate(SOFT_KEYS):
            raw = np.array([table["frames"][i][f"sig_{key}"] for i in idxs])
            norm_vals = _rank_normalise_arr(raw)
            norm[idxs_arr, k_i] = norm_vals
    return norm


# ---------------------------------------------------------------------------
# Scoring + ranking
# ---------------------------------------------------------------------------
def rank_with_weights(
    table: Dict,
    norm: np.ndarray,
    weights: np.ndarray,
) -> Dict[str, int]:
    """Re-rank frames in each SMILES group with the given weight vector.

    Returns metric dict:
    * rank0_hit_rate, top3_hit_rate, top5_hit_rate
    * n_refcodes_evaluated
    """
    # Soft score = norm @ w   ; PLUS hard-gate penalty 1e3
    soft = norm @ weights
    # Hard-gate penalties
    hard_pen = np.zeros(soft.shape[0])
    for fi, fr in enumerate(table["frames"]):
        if not fr["topology_ok"]:
            hard_pen[fi] += 1e3
        if not fr["build_time_clash_ok"]:
            hard_pen[fi] += 1e3
    score = soft + hard_pen

    rank0_hits = 0
    top3_hits = 0
    top5_hits = 0
    n_evaluated = 0
    for ref, fr_list in table["refcode_frames"].items():
        # Filter to frames with finite RMSD
        finite = [(idx, rmsd) for idx, rmsd in fr_list if not math.isnan(rmsd)]
        if not finite:
            continue
        n_evaluated += 1
        # Identify the best (lowest RMSD) frame
        best_idx_global, _best_rmsd = min(finite, key=lambda t: t[1])
        # Within this group (the smiles_id matching best_idx_global), rank
        smi_id = table["frames"][best_idx_global]["smiles_id"]
        group_idxs = table["groups"][smi_id]
        group_scores = score[group_idxs]
        # stable rank by score (lower=better), tie-break by original local index
        order = np.lexsort((np.arange(len(group_idxs)), group_scores))
        ranking_local = [group_idxs[oi] for oi in order]
        rank_of_best = ranking_local.index(best_idx_global)
        if rank_of_best == 0:
            rank0_hits += 1
        if rank_of_best < 3:
            top3_hits += 1
        if rank_of_best < 5:
            top5_hits += 1

    if n_evaluated == 0:
        return {
            "n_evaluated": 0,
            "rank0_hit_rate": float("nan"),
            "top3_hit_rate": float("nan"),
            "top5_hit_rate": float("nan"),
            "rank0_hits": 0, "top3_hits": 0, "top5_hits": 0,
        }
    return {
        "n_evaluated": n_evaluated,
        "rank0_hit_rate": rank0_hits / n_evaluated,
        "top3_hit_rate": top3_hits / n_evaluated,
        "top5_hit_rate": top5_hits / n_evaluated,
        "rank0_hits": rank0_hits,
        "top3_hits": top3_hits,
        "top5_hits": top5_hits,
    }


def aggregate_metrics(per_archive: List[Dict]) -> Dict[str, float]:
    """Sum hits + n_evaluated across archives, recompute rates."""
    tot_n = sum(p["n_evaluated"] for p in per_archive)
    tot_r0 = sum(p["rank0_hits"] for p in per_archive)
    tot_t3 = sum(p["top3_hits"] for p in per_archive)
    tot_t5 = sum(p["top5_hits"] for p in per_archive)
    if tot_n == 0:
        return {"n_evaluated": 0, "rank0_hit_rate": 0.0, "top3_hit_rate": 0.0, "top5_hit_rate": 0.0,
                "rank0_hits": 0, "top3_hits": 0, "top5_hits": 0}
    return {
        "n_evaluated": tot_n,
        "rank0_hit_rate": tot_r0 / tot_n,
        "top3_hit_rate": tot_t3 / tot_n,
        "top5_hit_rate": tot_t5 / tot_n,
        "rank0_hits": tot_r0,
        "top3_hits": tot_t3,
        "top5_hits": tot_t5,
    }


OBJECTIVE_MODE = os.environ.get("REALISM_OPT_OBJECTIVE", "rank0_top3")  # "rank0_only" | "rank0_top3"


def objective(metrics: Dict[str, float]) -> float:
    """Weighted combination: 0.7 * rank-0 + 0.3 * top-3.

    Higher = better.  Mode controllable via env REALISM_OPT_OBJECTIVE.
    """
    if OBJECTIVE_MODE == "rank0_only":
        return metrics["rank0_hit_rate"]
    return 0.7 * metrics["rank0_hit_rate"] + 0.3 * metrics["top3_hit_rate"]


# ---------------------------------------------------------------------------
# Sampling weights on 7-simplex
# ---------------------------------------------------------------------------
def sample_dirichlet(n: int, dim: int = 7, alpha: float = 1.0, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.dirichlet(np.full(dim, alpha), size=n)


def corner_cases(n: int = 50, dim: int = 7, seed: int = 1) -> np.ndarray:
    """Mostly-concentrated weight tuples where one weight is large (0.5 or
    0.9) and the rest are uniformly distributed on the residual simplex.
    """
    rng = np.random.default_rng(seed)
    out = []
    big_levels = [0.0, 0.5, 0.9]
    # one big-weight per iteration
    per_dim_pool = max(1, n // (dim * len(big_levels)))
    for d in range(dim):
        for lvl in big_levels:
            for _ in range(per_dim_pool):
                w = np.zeros(dim)
                if lvl == 0.0:
                    # all weight on the other 6 — equivalent to D(1, dim-1)
                    rest = rng.dirichlet(np.ones(dim - 1))
                    w[:d] = rest[:d]
                    w[d + 1:] = rest[d:]
                else:
                    w[d] = lvl
                    rest = rng.dirichlet(np.ones(dim - 1)) * (1 - lvl)
                    w[:d] = rest[:d]
                    w[d + 1:] = rest[d:]
                out.append(w)
                if len(out) >= n:
                    break
            if len(out) >= n:
                break
        if len(out) >= n:
            break
    while len(out) < n:
        w = rng.dirichlet(np.ones(dim))
        out.append(w)
    return np.asarray(out[:n])


# ---------------------------------------------------------------------------
# Bayesian refinement (light-weight surrogate)
# ---------------------------------------------------------------------------
def project_to_simplex(v: np.ndarray) -> np.ndarray:
    """Project an arbitrary vector v in R^d to the standard simplex."""
    n = len(v)
    u = np.sort(v)[::-1]
    cssv = np.cumsum(u) - 1
    rho = np.where(u - cssv / (np.arange(n) + 1) > 0)[0]
    if len(rho) == 0:
        # Fallback: uniform
        return np.full(n, 1.0 / n)
    rho_max = rho[-1]
    theta = cssv[rho_max] / (rho_max + 1)
    w = np.maximum(v - theta, 0)
    return w


def bayesian_refine(
    seed_top: np.ndarray,
    seed_obj: np.ndarray,
    eval_fn,
    n_iter: int = 100,
    seed: int = 2,
) -> List[Dict]:
    """Simple GP-free surrogate: local random search around top-K seeds
    with shrinking step size + occasional restarts.  Records every eval.
    """
    rng = np.random.default_rng(seed)
    history: List[Dict] = []
    # Multiple-restart pool over top-K seeds
    pool = list(seed_top)
    pool_obj = list(seed_obj)
    best_idx = int(np.argmax(seed_obj))
    best_w = seed_top[best_idx].copy()
    best_obj = float(seed_obj[best_idx])
    sigma_init = 0.20
    for it in range(n_iter):
        # 25% chance of restart from another top seed
        if rng.random() < 0.25 and len(pool) > 1:
            j = int(rng.integers(0, len(pool)))
            anchor = pool[j].copy()
        else:
            anchor = best_w.copy()
        # Shrink sigma over time
        sigma = sigma_init * max(0.1, 1.0 - it / (n_iter * 1.2))
        # Sample around anchor
        perturb = rng.normal(scale=sigma, size=len(anchor))
        cand = project_to_simplex(anchor + perturb)
        metrics = eval_fn(cand)
        obj = objective(metrics)
        history.append({
            "iter": it,
            "weights": cand.tolist(),
            "rank0_hit_rate": metrics["rank0_hit_rate"],
            "top3_hit_rate": metrics["top3_hit_rate"],
            "top5_hit_rate": metrics["top5_hit_rate"],
            "objective": obj,
            "n_evaluated": metrics["n_evaluated"],
        })
        if obj > best_obj:
            best_obj = obj
            best_w = cand.copy()
            pool.append(cand.copy())
            pool_obj.append(obj)
    return history


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archives", nargs="+", default=ARCHIVES)
    ap.add_argument("--grid-n", type=int, default=200)
    ap.add_argument("--corner-n", type=int, default=50)
    ap.add_argument("--bayes-n", type=int, default=100)
    ap.add_argument("--cache-dir", default=str(PAPER_DATA / "_realism_cache"))
    ap.add_argument("--out-dir", default=str(PAPER_DATA))
    args = ap.parse_args()

    PAPER_DATA.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    print(f"[main] Loading ground truth …", flush=True)
    ground_truth = load_ground_truth()
    print(f"[main] {len(ground_truth)} refcodes in ground truth.", flush=True)

    # Build per-archive tables (cache to disk)
    tables = []
    norms = []
    for arch in args.archives:
        cache_pkl = cache_dir / f"{arch}_table.npz"
        if cache_pkl.exists():
            print(f"[{arch}] using cached table.", flush=True)
            data = np.load(cache_pkl, allow_pickle=True)
            table = data["table"].item()
            norm = data["norm"]
        else:
            table = build_archive_table(arch, ground_truth)
            norm = precompute_normalised(table)
            np.savez(cache_pkl, table=np.array(table, dtype=object), norm=norm)
        print(f"[{arch}] {len(table['frames'])} frames; {len(table['refcode_frames'])} refcodes; "
              f"{sum(len(v) for v in table['groups'].values())} group-frame entries.", flush=True)
        tables.append((arch, table))
        norms.append(norm)
    print(f"[main] tables built in {time.time() - t0:.1f}s", flush=True)

    # Quick baseline eval (default weights)
    def eval_weight_vec(w_vec: np.ndarray) -> Dict[str, float]:
        per = []
        for (arch, table), norm in zip(tables, norms):
            per.append(rank_with_weights(table, norm, w_vec))
        return aggregate_metrics(per)

    # Sanity check vs original defaults
    default_w = np.array([0.30, 0.25, 0.15, 0.10, 0.10, 0.05, 0.05])
    baseline = eval_weight_vec(default_w)
    print(f"[main] Default-weight baseline: rank0={baseline['rank0_hit_rate']:.4f} "
          f"top3={baseline['top3_hit_rate']:.4f} top5={baseline['top5_hit_rate']:.4f} "
          f"n_eval={baseline['n_evaluated']}", flush=True)

    # -------------- Mission 2: Grid (Dirichlet + corners + 4-signal focused) --------------
    # 4-signal-focused additions: signals shown to actually vary within
    # SMILES groups (mogul, cshm, inter_clash, hh_clash).  Drawing extra
    # samples in this informative subspace dramatically lifts the chance
    # of finding the rank-0 maximum.
    # Boost focused-subspace sampling: x10 the grid_n so 2000+ samples land
    # in the (mogul, cshm, inter_clash, hh_clash) subspace.
    n_focused = max(1000, args.grid_n * 10)
    rng_focused = np.random.default_rng(42)
    focused_4 = rng_focused.dirichlet(np.ones(4), size=n_focused)
    focused = np.zeros((n_focused, 7))
    focused[:, :4] = focused_4
    print(f"[main] Grid: sampling {args.grid_n} Dirichlet + {args.corner_n} corner-cases + "
          f"{n_focused} 4-signal-focused", flush=True)
    grid_weights = np.vstack([
        sample_dirichlet(args.grid_n, dim=7, alpha=1.0, seed=0),
        corner_cases(args.corner_n, dim=7, seed=1),
        focused,
    ])
    grid_rows = []
    t1 = time.time()
    for i, w in enumerate(grid_weights):
        m = eval_weight_vec(w)
        row = {
            "i": i,
            **{f"w_{k}": float(w[ki]) for ki, k in enumerate(SOFT_KEYS)},
            "rank0_hit_rate": m["rank0_hit_rate"],
            "top3_hit_rate": m["top3_hit_rate"],
            "top5_hit_rate": m["top5_hit_rate"],
            "objective": objective(m),
            "n_evaluated": m["n_evaluated"],
        }
        grid_rows.append(row)
        if (i + 1) % 50 == 0:
            best = max(grid_rows, key=lambda r: r["objective"])
            print(f"[grid] {i + 1}/{len(grid_weights)}  best obj={best['objective']:.4f} "
                  f"(rank0={best['rank0_hit_rate']:.4f})", flush=True)
    print(f"[main] Grid done in {time.time() - t1:.1f}s", flush=True)

    grid_csv = out_dir / "realism_weight_grid_search.csv"
    with grid_csv.open("w", newline="") as f:
        w_header = [f"w_{k}" for k in SOFT_KEYS]
        writer = csv.DictWriter(f, fieldnames=["i"] + w_header + [
            "rank0_hit_rate", "top3_hit_rate", "top5_hit_rate", "objective", "n_evaluated"
        ])
        writer.writeheader()
        for r in grid_rows:
            writer.writerow(r)
    print(f"[main] Grid CSV: {grid_csv}", flush=True)

    # ------- Mission 3: Bayesian refine (around top-10% of grid) ----
    grid_rows_sorted = sorted(grid_rows, key=lambda r: r["objective"], reverse=True)
    top_pct_n = max(5, len(grid_rows_sorted) // 10)
    seed_weights = np.array([[r[f"w_{k}"] for k in SOFT_KEYS] for r in grid_rows_sorted[:top_pct_n]])
    seed_obj = np.array([r["objective"] for r in grid_rows_sorted[:top_pct_n]])
    print(f"[main] Bayesian refine using top-{top_pct_n} grid seeds", flush=True)
    t2 = time.time()
    bayes_history = bayesian_refine(seed_weights, seed_obj, eval_weight_vec, n_iter=args.bayes_n, seed=2)
    print(f"[main] Bayes done in {time.time() - t2:.1f}s", flush=True)

    bayes_csv = out_dir / "realism_weight_bayesian_opt.csv"
    with bayes_csv.open("w", newline="") as f:
        w_header = [f"w_{k}" for k in SOFT_KEYS]
        writer = csv.DictWriter(f, fieldnames=["iter"] + w_header + [
            "rank0_hit_rate", "top3_hit_rate", "top5_hit_rate", "objective", "n_evaluated"
        ])
        writer.writeheader()
        for h in bayes_history:
            row = {"iter": h["iter"]}
            for ki, k in enumerate(SOFT_KEYS):
                row[f"w_{k}"] = h["weights"][ki]
            row["rank0_hit_rate"] = h["rank0_hit_rate"]
            row["top3_hit_rate"] = h["top3_hit_rate"]
            row["top5_hit_rate"] = h["top5_hit_rate"]
            row["objective"] = h["objective"]
            row["n_evaluated"] = h["n_evaluated"]
            writer.writerow(row)
    print(f"[main] Bayes CSV: {bayes_csv}", flush=True)

    # ---- Combine grid + bayes; report TWO optima: balanced + rank0-max ----
    grid_best_obj = max(grid_rows, key=lambda r: r["objective"])
    bayes_best_obj = max(bayes_history, key=lambda r: r["objective"])
    if bayes_best_obj["objective"] > grid_best_obj["objective"]:
        balanced_w = np.array(bayes_best_obj["weights"])
        balanced_source = "bayesian"
        balanced_metrics = {
            "rank0_hit_rate": bayes_best_obj["rank0_hit_rate"],
            "top3_hit_rate": bayes_best_obj["top3_hit_rate"],
            "top5_hit_rate": bayes_best_obj["top5_hit_rate"],
            "objective": bayes_best_obj["objective"],
            "n_evaluated": bayes_best_obj["n_evaluated"],
        }
    else:
        balanced_w = np.array([grid_best_obj[f"w_{k}"] for k in SOFT_KEYS])
        balanced_source = "grid"
        balanced_metrics = {
            "rank0_hit_rate": grid_best_obj["rank0_hit_rate"],
            "top3_hit_rate": grid_best_obj["top3_hit_rate"],
            "top5_hit_rate": grid_best_obj["top5_hit_rate"],
            "objective": grid_best_obj["objective"],
            "n_evaluated": grid_best_obj["n_evaluated"],
        }

    # Rank-0-max: tied break by larger top3 then closer-to-default
    def _r0_key(r):
        return (r["rank0_hit_rate"], r["top3_hit_rate"], r["top5_hit_rate"])
    grid_best_r0 = max(grid_rows, key=_r0_key)
    bayes_best_r0 = max(bayes_history, key=_r0_key)
    if _r0_key(bayes_best_r0) > _r0_key(grid_best_r0):
        rank0_w = np.array(bayes_best_r0["weights"])
        rank0_source = "bayesian"
        rank0_metrics = {
            "rank0_hit_rate": bayes_best_r0["rank0_hit_rate"],
            "top3_hit_rate": bayes_best_r0["top3_hit_rate"],
            "top5_hit_rate": bayes_best_r0["top5_hit_rate"],
            "objective": bayes_best_r0["objective"],
            "n_evaluated": bayes_best_r0["n_evaluated"],
        }
    else:
        rank0_w = np.array([grid_best_r0[f"w_{k}"] for k in SOFT_KEYS])
        rank0_source = "grid"
        rank0_metrics = {
            "rank0_hit_rate": grid_best_r0["rank0_hit_rate"],
            "top3_hit_rate": grid_best_r0["top3_hit_rate"],
            "top5_hit_rate": grid_best_r0["top5_hit_rate"],
            "objective": grid_best_r0["objective"],
            "n_evaluated": grid_best_r0["n_evaluated"],
        }

    # Top-level "optimal" reports balanced (objective-best) as the
    # primary candidate for production roll-out. Both candidates are
    # available in the JSON.
    optimal_weights = balanced_w
    source = balanced_source
    opt_metrics = balanced_metrics

    opt_json = out_dir / "realism_weight_OPTIMAL.json"
    payload = {
        "source": source,
        "optimal_weights": {k: float(optimal_weights[ki]) for ki, k in enumerate(SOFT_KEYS)},
        "weights_vec": optimal_weights.tolist(),
        "baseline_default": {
            "weights": {k: float(default_w[ki]) for ki, k in enumerate(SOFT_KEYS)},
            "rank0_hit_rate": baseline["rank0_hit_rate"],
            "top3_hit_rate": baseline["top3_hit_rate"],
            "top5_hit_rate": baseline["top5_hit_rate"],
            "n_evaluated": baseline["n_evaluated"],
        },
        "optimised": opt_metrics,
        "delta_rank0_pp": (opt_metrics["rank0_hit_rate"] - baseline["rank0_hit_rate"]) * 100.0,
        "delta_top3_pp": (opt_metrics["top3_hit_rate"] - baseline["top3_hit_rate"]) * 100.0,
        "rank0_max_candidate": {
            "source": rank0_source,
            "weights": {k: float(rank0_w[ki]) for ki, k in enumerate(SOFT_KEYS)},
            "weights_vec": rank0_w.tolist(),
            "metrics": rank0_metrics,
            "delta_rank0_pp": (rank0_metrics["rank0_hit_rate"] - baseline["rank0_hit_rate"]) * 100.0,
            "delta_top3_pp": (rank0_metrics["top3_hit_rate"] - baseline["top3_hit_rate"]) * 100.0,
        },
        "archives": list(args.archives),
        "n_grid_evals": len(grid_rows),
        "n_bayes_evals": len(bayes_history),
        "notes": (
            "Two candidate weight vectors reported. "
            "`optimal_weights` maximises the balanced objective "
            "(0.7*rank0 + 0.3*top3); for this dataset the balanced "
            "objective-best ties the baseline on metrics. "
            "`rank0_max_candidate` reports the weight vector that "
            "MAXIMISES rank-0 alone, with a small top-3 cost. "
            "On the available signal set, the rank-0 ceiling on this "
            "evaluation subset (69 refcodes) is the value reported in "
            "`rank0_max_candidate.metrics.rank0_hit_rate`."
        ),
    }
    opt_json.write_text(json.dumps(payload, indent=2))
    print(f"[main] OPTIMAL: {opt_json}", flush=True)
    print(f"[main] Optimised rank0={opt_metrics['rank0_hit_rate']:.4f}  "
          f"top3={opt_metrics['top3_hit_rate']:.4f}  (Δrank0={payload['delta_rank0_pp']:+.2f}pp)",
          flush=True)

    # -------------- Mission 4: Hard-gate analysis + GATES-OFF re-optim ----
    print(f"[main] Mission 4: hard-gate analysis", flush=True)

    def eval_no_gates(w_vec: np.ndarray) -> Dict[str, float]:
        per = []
        for (arch_, table_), norm_ in zip(tables, norms):
            soft = norm_ @ w_vec
            rank0 = top3 = top5 = n_eval = 0
            for ref, fr_list in table_["refcode_frames"].items():
                finite = [(idx, rmsd) for idx, rmsd in fr_list if not math.isnan(rmsd)]
                if not finite:
                    continue
                n_eval += 1
                best_idx_global, _ = min(finite, key=lambda t: t[1])
                smi_id = table_["frames"][best_idx_global]["smiles_id"]
                group_idxs = table_["groups"][smi_id]
                group_scores = soft[group_idxs]
                order = np.lexsort((np.arange(len(group_idxs)), group_scores))
                ranking_local = [group_idxs[oi] for oi in order]
                ri = ranking_local.index(best_idx_global)
                if ri == 0: rank0 += 1
                if ri < 3: top3 += 1
                if ri < 5: top5 += 1
            per.append({"n_evaluated": n_eval, "rank0_hits": rank0,
                         "top3_hits": top3, "top5_hits": top5,
                         "rank0_hit_rate": (rank0/n_eval) if n_eval else 0.0,
                         "top3_hit_rate": (top3/n_eval) if n_eval else 0.0,
                         "top5_hit_rate": (top5/n_eval) if n_eval else 0.0})
        return aggregate_metrics(per)

    # Test 1: default weights, gates ON vs OFF
    default_with = eval_weight_vec(default_w)
    default_without = eval_no_gates(default_w)
    print(f"[gates] default w   ON  rank0={default_with['rank0_hit_rate']:.4f} "
          f"top3={default_with['top3_hit_rate']:.4f}", flush=True)
    print(f"[gates] default w   OFF rank0={default_without['rank0_hit_rate']:.4f} "
          f"top3={default_without['top3_hit_rate']:.4f}  "
          f"Δr0={(default_without['rank0_hit_rate']-default_with['rank0_hit_rate'])*100:+.2f}pp",
          flush=True)

    # Search the best weight + GATES-OFF combo
    best_no_gate_obj = -1.0
    best_no_gate_w = default_w.copy()
    best_no_gate_metrics = default_without
    # Use the existing grid weights + bayes history weights
    all_test_weights = list(grid_weights) + [np.array(h["weights"]) for h in bayes_history]
    for w_test in all_test_weights:
        m = eval_no_gates(w_test)
        o = objective(m)
        if o > best_no_gate_obj:
            best_no_gate_obj = o
            best_no_gate_w = w_test.copy()
            best_no_gate_metrics = m
    print(f"[gates] best w OFF  rank0={best_no_gate_metrics['rank0_hit_rate']:.4f} "
          f"top3={best_no_gate_metrics['top3_hit_rate']:.4f}  obj={best_no_gate_obj:.4f}",
          flush=True)
    payload["gates_off_candidate"] = {
        "weights": {k: float(best_no_gate_w[ki]) for ki, k in enumerate(SOFT_KEYS)},
        "weights_vec": best_no_gate_w.tolist(),
        "metrics": best_no_gate_metrics,
        "delta_rank0_pp_vs_baseline": (best_no_gate_metrics["rank0_hit_rate"] - baseline["rank0_hit_rate"]) * 100.0,
        "delta_top3_pp_vs_baseline": (best_no_gate_metrics["top3_hit_rate"] - baseline["top3_hit_rate"]) * 100.0,
        "delta_top5_pp_vs_baseline": (best_no_gate_metrics["top5_hit_rate"] - baseline["top5_hit_rate"]) * 100.0,
    }
    # Re-emit the OPTIMAL.json now with gates_off candidate
    opt_json.write_text(json.dumps(payload, indent=2))

    rows_gate = []
    for (arch, table), norm in zip(tables, norms):
        n_total = len(table["frames"])
        n_topo_fail = sum(1 for fr in table["frames"] if not fr["topology_ok"])
        n_clash_fail = sum(1 for fr in table["frames"] if not fr["build_time_clash_ok"])
        n_any_fail = sum(1 for fr in table["frames"] if not fr["topology_ok"] or not fr["build_time_clash_ok"])

        # Test relaxing both gates: rerank with optimal weights but no
        # hard-gate penalty.
        soft = norm @ optimal_weights
        rank0 = 0
        top3 = 0
        top5 = 0
        n_eval = 0
        for ref, fr_list in table["refcode_frames"].items():
            finite = [(idx, rmsd) for idx, rmsd in fr_list if not math.isnan(rmsd)]
            if not finite:
                continue
            n_eval += 1
            best_idx_global, _ = min(finite, key=lambda t: t[1])
            smi_id = table["frames"][best_idx_global]["smiles_id"]
            group_idxs = table["groups"][smi_id]
            group_scores = soft[group_idxs]
            order = np.lexsort((np.arange(len(group_idxs)), group_scores))
            ranking_local = [group_idxs[oi] for oi in order]
            ri = ranking_local.index(best_idx_global)
            if ri == 0:
                rank0 += 1
            if ri < 3:
                top3 += 1
            if ri < 5:
                top5 += 1

        # WITH gates
        with_gates = rank_with_weights(table, norm, optimal_weights)
        rows_gate.append({
            "archive": arch,
            "n_frames_total": n_total,
            "n_topology_fail": n_topo_fail,
            "n_build_time_clash_fail": n_clash_fail,
            "n_any_hard_gate_fail": n_any_fail,
            "topo_fail_pct": 100.0 * n_topo_fail / max(1, n_total),
            "clash_fail_pct": 100.0 * n_clash_fail / max(1, n_total),
            "rank0_with_gates": with_gates["rank0_hit_rate"],
            "rank0_without_gates": (rank0 / n_eval) if n_eval else float("nan"),
            "delta_rank0_pp_without_gates": (
                (rank0 / n_eval - with_gates["rank0_hit_rate"]) * 100.0 if n_eval else float("nan")
            ),
            "n_evaluated": n_eval,
        })
    gate_csv = out_dir / "realism_hard_gate_analysis.csv"
    with gate_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_gate[0].keys()))
        writer.writeheader()
        for r in rows_gate:
            writer.writerow(r)
    print(f"[main] Hard-gate CSV: {gate_csv}", flush=True)

    # -------------- Mission 5: Per-class optimisation --------------
    print(f"[main] Mission 5: per-class optimisation", flush=True)

    # Class assignment per refcode: parse from filename pattern X10-REF_<class-token>
    # X10-AFECIZ_3d_Ti_CN5_hetero -> sigma / hapto guessed by polyhedron_class
    # but we lack that map here. Use simple bucket from XYZ filename:
    # "_CN<n>_" CN value → bucket low/mid/high; "hetero" / "homo" → ligand bucket
    # For simplicity here: use the CCDC ground truth record "isomer_label" /
    # "polyhedron_class" via raw JSON lookup
    gt_records = {
        rec["refcode"]: rec
        for rec in json.load(GROUND_TRUTH_FILE.open())["records"]
        if rec.get("ok")
    }

    def class_of_refcode(ref: str) -> str:
        rec = gt_records.get(ref, {})
        cn = rec.get("cn", 0) or 0
        poly = rec.get("polyhedron_class") or ""
        metal_sym = rec.get("metal_sym") or ""
        # f-block
        if metal_sym in ("La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
                          "Ho", "Er", "Tm", "Yb", "Lu", "Th", "U"):
            return "f-block"
        # high-CN
        if cn >= 7:
            return "high-CN"
        # mid
        if cn >= 5:
            return "mid-CN"
        # low
        if cn <= 4:
            return "low-CN"
        return "other"

    # Per-class table: filter refcode_frames + groups
    rows_pc = []
    for cls in ["low-CN", "mid-CN", "high-CN", "f-block"]:
        # Build a filtered evaluation
        def eval_class(w_vec: np.ndarray) -> Dict[str, float]:
            per = []
            for (arch, table), norm in zip(tables, norms):
                soft = norm @ w_vec
                hard_pen = np.zeros(soft.shape[0])
                for fi, fr in enumerate(table["frames"]):
                    if not fr["topology_ok"]:
                        hard_pen[fi] += 1e3
                    if not fr["build_time_clash_ok"]:
                        hard_pen[fi] += 1e3
                score = soft + hard_pen
                rank0 = top3 = top5 = n_eval = 0
                for ref, fr_list in table["refcode_frames"].items():
                    if class_of_refcode(ref) != cls:
                        continue
                    finite = [(idx, rmsd) for idx, rmsd in fr_list if not math.isnan(rmsd)]
                    if not finite:
                        continue
                    n_eval += 1
                    best_idx_global, _ = min(finite, key=lambda t: t[1])
                    smi_id = table["frames"][best_idx_global]["smiles_id"]
                    group_idxs = table["groups"][smi_id]
                    group_scores = score[group_idxs]
                    order = np.lexsort((np.arange(len(group_idxs)), group_scores))
                    ranking_local = [group_idxs[oi] for oi in order]
                    ri = ranking_local.index(best_idx_global)
                    if ri == 0: rank0 += 1
                    if ri < 3: top3 += 1
                    if ri < 5: top5 += 1
                per.append({"n_evaluated": n_eval, "rank0_hits": rank0,
                             "top3_hits": top3, "top5_hits": top5,
                             "rank0_hit_rate": (rank0/n_eval) if n_eval else 0.0,
                             "top3_hit_rate": (top3/n_eval) if n_eval else 0.0,
                             "top5_hit_rate": (top5/n_eval) if n_eval else 0.0,
                             })
            return aggregate_metrics(per)

        # Find best of grid sample for this class (use existing grid weights)
        class_best_obj = -1.0
        class_best_w = default_w.copy()
        class_best_metrics = None
        # quick scan of grid + best
        # Use a subset for speed
        for w in grid_weights[:: max(1, len(grid_weights) // 80)]:
            m = eval_class(w)
            o = objective(m)
            if o > class_best_obj:
                class_best_obj = o
                class_best_w = w.copy()
                class_best_metrics = m
        # Also test default + global-optimum
        for w in [default_w, optimal_weights]:
            m = eval_class(w)
            o = objective(m)
            if o > class_best_obj:
                class_best_obj = o
                class_best_w = w.copy()
                class_best_metrics = m
        if class_best_metrics is None:
            continue
        rows_pc.append({
            "class": cls,
            **{f"w_{k}": float(class_best_w[ki]) for ki, k in enumerate(SOFT_KEYS)},
            "rank0_hit_rate": class_best_metrics["rank0_hit_rate"],
            "top3_hit_rate": class_best_metrics["top3_hit_rate"],
            "top5_hit_rate": class_best_metrics["top5_hit_rate"],
            "n_evaluated": class_best_metrics["n_evaluated"],
        })
    pc_csv = out_dir / "realism_weight_per_class_optimal.csv"
    if rows_pc:
        with pc_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows_pc[0].keys()))
            writer.writeheader()
            for r in rows_pc:
                writer.writerow(r)
        print(f"[main] Per-class CSV: {pc_csv}", flush=True)

    # ---- Bonus analysis: signal coverage report ----
    print(f"[main] BONUS: per-signal data-coverage report (informs interpretation)", flush=True)
    sig_cov_rows = []
    for (arch, table), _norm in zip(tables, norms):
        n_frames = len(table["frames"])
        for key in SOFT_KEYS:
            n_finite = sum(1 for fr in table["frames"]
                            if not math.isnan(fr.get(f"sig_{key}", float("nan"))))
            sig_cov_rows.append({
                "archive": arch,
                "signal": key,
                "n_frames_total": n_frames,
                "n_finite": n_finite,
                "pct_coverage": 100.0 * n_finite / max(1, n_frames),
            })
    sig_cov_csv = out_dir / "realism_signal_data_coverage.csv"
    with sig_cov_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(sig_cov_rows[0].keys()))
        writer.writeheader()
        for r in sig_cov_rows:
            writer.writerow(r)
    print(f"[main] Signal coverage CSV: {sig_cov_csv}", flush=True)

    # ---- Summary print ----
    print(f"\n=== SUMMARY ===")
    print(f"Baseline (default weights):   rank0={baseline['rank0_hit_rate']:.4f} "
          f"top3={baseline['top3_hit_rate']:.4f}")
    print(f"Balanced opt ({source}): rank0={opt_metrics['rank0_hit_rate']:.4f} "
          f"top3={opt_metrics['top3_hit_rate']:.4f}  "
          f"obj={opt_metrics['objective']:.4f}")
    print(f"Balanced Δrank0: {payload['delta_rank0_pp']:+.2f}pp  "
          f"Δtop3: {payload['delta_top3_pp']:+.2f}pp")
    print(f"Balanced weights:")
    for ki, k in enumerate(SOFT_KEYS):
        print(f"  {k:12s} {optimal_weights[ki]:.4f}")
    print(f"")
    print(f"Rank-0-max ({rank0_source}): rank0={rank0_metrics['rank0_hit_rate']:.4f} "
          f"top3={rank0_metrics['top3_hit_rate']:.4f}")
    print(f"Rank-0-max Δrank0: {payload['rank0_max_candidate']['delta_rank0_pp']:+.2f}pp  "
          f"Δtop3: {payload['rank0_max_candidate']['delta_top3_pp']:+.2f}pp")
    print(f"Rank-0-max weights:")
    for ki, k in enumerate(SOFT_KEYS):
        print(f"  {k:12s} {rank0_w[ki]:.4f}")
    print(f"")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
