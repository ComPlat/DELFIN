"""Full-CCDC XRD isomer + conformer recall metric (production runtime).

Successor to ``xrd_recall_metric.py`` (48 hardcoded refcodes). This module
uses the **full CSD 2026.1 TMC family table** built at build-time by
``agent_workspace/quality_framework/scripts/build_per_smiles_ccdc_families.py``
against the ``ccdc_tmc_index.jsonl`` (~500k TMC entries).

Per-SMILES Definition (User-corrected 2026-06-04 17:40):
    "48 CCDC ground truth ist unnötig wir haben jetzt alle ccdc stukturen"

For each SMILES in the master pool with >=1 CCDC family member:

    isomer_recall(smiles) =
        |DELFIN-emitted-isomers ∩ CCDC-family-isomers|
        / |CCDC-family-isomers|

    conformer_recall(smiles, rmsd_threshold=0.5) =
        fraction of CCDC family refcodes for which SOME DELFIN emission
        has heavy-atom Kabsch RMSD < threshold

Global metric = mean over SMILES with >=1 CCDC match (option:
weighted by family size).

Runtime constraints:
    - NO CCDC import at import time (legal-separation doctrine).
    - Reads only the pre-built JSONL family table + per-refcode positions
      from ccdc_tmc_index.jsonl.

Env-flags:
    DELFIN_FFFREE_XRD_RECALL_TABLE_PATH : path to per_smiles_ccdc_families.jsonl
    DELFIN_FFFREE_XRD_RECALL_INDEX_PATH : path to ccdc_tmc_index.jsonl
    DELFIN_FFFREE_XRD_RECALL_RMSD_A     : default 0.5
    DELFIN_FFFREE_XRD_RECALL_WEIGHTED   : "1" = family-size-weighted mean

Default behaviour (no env flags set): module is import-safe but
``compute_*`` returns NaN unless paths are explicitly provided.

Public API:
    load_ccdc_family_table(path) -> dict[smiles] = list[match]
    load_ccdc_tmc_index(path)    -> dict[refcode] = record  (positions, symbols)
    compute_isomer_recall(emitted_by_smiles, family_table) -> float
    compute_conformer_recall(emitted_by_smiles, family_table, tmc_index,
                              rmsd_threshold=0.5) -> float
    per_smiles_recall_report(emitted_by_smiles, family_table, tmc_index)
        -> list[dict] (one per SMILES)
"""
from __future__ import annotations

import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

try:
    import numpy as np
except Exception:  # pragma: no cover
    np = None  # type: ignore


# Reuse existing classification + RMSD helpers from xrd_recall_metric (no
# circular dep — that's a sibling module that does NOT import this one).
try:
    from delfin.fffree.xrd_recall_metric import (
        classify_isomer,
        kabsch_rmsd_heavy,
        parse_xyz,
        METALS,
        DONOR_ELEMS,
    )
    _HELPERS_AVAILABLE = True
except Exception:  # noqa
    _HELPERS_AVAILABLE = False
    METALS = {
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu", "Th", "U",
    }
    DONOR_ELEMS = {"C", "N", "O", "P", "S", "Cl", "Br", "I", "F", "Se", "As", "B", "Si", "Te"}

    _COV = {
        "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
        "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20,
        "Br": 1.20, "Te": 1.38, "I": 1.39,
        "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
        "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
        "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
        "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
        "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Pb": 1.46, "La": 2.07,
        "Ce": 2.04, "Sn": 1.39,
    }
    _BOND_TOL = 1.30

    def parse_xyz(text):
        """Minimal XYZ parser fallback."""
        atoms = []
        syms, P = [], []
        for line in text.splitlines():
            p = line.split()
            if len(p) < 4:
                continue
            try:
                x, y, z = float(p[1]), float(p[2]), float(p[3])
            except ValueError:
                continue
            mm = re.match(r"([A-Z][a-z]?)", p[0])
            if not mm:
                continue
            syms.append(mm.group(1))
            P.append([x, y, z])
        return syms, (np.asarray(P, dtype=float) if (np is not None and P) else P)

    def _coord_sphere_fb(syms, P):
        if np is None or len(syms) == 0:
            return None
        metals = [i for i, s in enumerate(syms) if s in METALS]
        if not metals:
            return None
        best = None
        for mi in metals:
            mp = P[mi]
            mr = _COV.get(syms[mi], 1.5)
            donors = []
            for j, s in enumerate(syms):
                if j == mi or s == "H" or s not in DONOR_ELEMS:
                    continue
                r = _COV.get(s, 1.0)
                d = float(np.linalg.norm(P[j] - mp))
                cutoff = (mr + r) * _BOND_TOL
                if s == "C":
                    cutoff = min(cutoff, 2.45)
                if d < cutoff:
                    donors.append((j, s, d))
            donors.sort(key=lambda t: t[2])
            if best is None or len(donors) > len(best[1]):
                best = (mi, donors)
        if best is None:
            return None
        return {
            "metal_idx": best[0],
            "donor_idx": [t[0] for t in best[1]],
            "donor_syms": [t[1] for t in best[1]],
        }

    def classify_isomer(syms, P):
        """Fallback classify_isomer — returns CN<n>-<poly> string."""
        if np is None or not isinstance(P, np.ndarray):
            try:
                P = np.asarray(P, dtype=float)
            except Exception:
                return "unknown"
        ci = _coord_sphere_fb(syms, P)
        if not ci:
            return "no-metal"
        cn = len(ci["donor_syms"])
        poly = {2: "linear", 3: "trigonal", 4: "tet-or-sp", 5: "TBP-or-SPY",
                6: "OH", 7: "PBP", 8: "Cube", 9: "TTP",
                10: "BAS", 11: "edge", 12: "ico"}.get(cn, f"CN{cn}")
        return f"CN{cn}-{poly}"

    def kabsch_rmsd_heavy(syms_a, Pa, syms_b, Pb):
        """Heavy-atom Kabsch RMSD with greedy element-matched correspondence."""
        if np is None:
            return float("nan")
        Pa = np.asarray(Pa, dtype=float)
        Pb = np.asarray(Pb, dtype=float)
        keep_a = [i for i, s in enumerate(syms_a) if s != "H"]
        keep_b = [i for i, s in enumerate(syms_b) if s != "H"]
        if not keep_a or not keep_b:
            return float("nan")
        Pa_h = Pa[keep_a] - Pa[keep_a].mean(0)
        Pb_h = Pb[keep_b] - Pb[keep_b].mean(0)
        syms_ah = [syms_a[i] for i in keep_a]
        syms_bh = [syms_b[i] for i in keep_b]
        used_b = [False] * len(Pb_h)
        pairs = []
        for i in sorted(range(len(Pa_h)), key=lambda k: (syms_ah[k], k)):
            best_j, best_d = -1, 1e9
            for j in range(len(Pb_h)):
                if used_b[j] or syms_bh[j] != syms_ah[i]:
                    continue
                d = float(np.linalg.norm(Pa_h[i] - Pb_h[j]))
                if d < best_d:
                    best_d, best_j = d, j
            if best_j >= 0:
                used_b[best_j] = True
                pairs.append((i, best_j))
        if len(pairs) < 3:
            return float("nan")
        Pa_m = np.array([Pa_h[i] for i, _ in pairs])
        Pb_m = np.array([Pb_h[j] for _, j in pairs])
        Pa_m = Pa_m - Pa_m.mean(0)
        Pb_m = Pb_m - Pb_m.mean(0)
        H = Pa_m.T @ Pb_m
        U, S, Vt = np.linalg.svd(H)
        d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        Pr = Pa_m @ R.T
        return float(math.sqrt(((Pr - Pb_m) ** 2).sum(1).mean()))


ENV_TABLE_PATH = "DELFIN_FFFREE_XRD_RECALL_TABLE_PATH"
ENV_INDEX_PATH = "DELFIN_FFFREE_XRD_RECALL_INDEX_PATH"
ENV_RMSD_A = "DELFIN_FFFREE_XRD_RECALL_RMSD_A"
ENV_WEIGHTED = "DELFIN_FFFREE_XRD_RECALL_WEIGHTED"

DEFAULT_RMSD_A = 0.5


__all__ = [
    "ENV_TABLE_PATH",
    "ENV_INDEX_PATH",
    "ENV_RMSD_A",
    "ENV_WEIGHTED",
    "DEFAULT_RMSD_A",
    "load_ccdc_family_table",
    "load_ccdc_tmc_index",
    "load_master_label_to_smiles",
    "group_archive_by_smiles",
    "compute_isomer_recall",
    "compute_conformer_recall",
    "per_smiles_recall_report",
    "score_archive",
]


# ---------------------------------------------------------------------------
# Table loaders (cache by path)
# ---------------------------------------------------------------------------
_TABLE_CACHE: dict = {}
_INDEX_CACHE: dict = {}


def load_ccdc_family_table(path: Optional[str] = None) -> Dict[str, dict]:
    """Load per-SMILES CCDC family JSONL.

    Returns dict keyed by SMILES string -> {smiles_label, n_matches, matches,
    match_method}. Cached by path.
    """
    if path is None:
        path = os.environ.get(ENV_TABLE_PATH)
    if not path:
        return {}
    if path in _TABLE_CACHE:
        return _TABLE_CACHE[path]
    p = Path(path)
    if not p.exists():
        _TABLE_CACHE[path] = {}
        return {}
    table = {}
    with p.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            smi = rec.get("smiles")
            if smi:
                table[smi] = rec
    _TABLE_CACHE[path] = table
    return table


def load_ccdc_tmc_index(path: Optional[str] = None) -> Dict[str, dict]:
    """Load CCDC TMC index JSONL keyed by refcode (positions/symbols available)."""
    if path is None:
        path = os.environ.get(ENV_INDEX_PATH)
    if not path:
        return {}
    if path in _INDEX_CACHE:
        return _INDEX_CACHE[path]
    p = Path(path)
    if not p.exists():
        _INDEX_CACHE[path] = {}
        return {}
    idx = {}
    with p.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            rc = rec.get("refcode")
            if rc:
                idx[rc] = rec
    _INDEX_CACHE[path] = idx
    return idx


# ---------------------------------------------------------------------------
# Core recall computations
# ---------------------------------------------------------------------------
def _emitted_isomer_labels(xyz_paths, max_per_smiles=200):
    """Classify each XYZ to an isomer label using classify_isomer."""
    labels = set()
    if np is None:
        return labels
    count = 0
    for xp in xyz_paths:
        if count >= max_per_smiles:
            break
        try:
            text = Path(xp).read_text()
            syms, P = parse_xyz(text)
            if not syms or P is None or len(P) == 0:
                continue
            label = classify_isomer(syms, P)
            labels.add(label)
            count += 1
        except Exception:
            continue
    return labels


def _isomer_match(delfin_iso: str, ccdc_iso: str, ccdc_cn: int) -> bool:
    """Match a DELFIN-emitted isomer label to a CCDC family isomer label.

    CCDC family labels are CN-based polyhedron strings like
    ``"octahedral_or_TP"`` / ``"TBP_or_SPY"``; DELFIN classifier produces
    ``"CN6-OH"`` / ``"CN5-TBP-or-SPY"`` / etc.

    Strategy:
      1) Direct equality.
      2) DELFIN-string starts with "CN<n>-" with same n as CCDC's `cn`.
      3) Substring of polyhedron keywords (OH/TBP/SPY/TET/SP/PBP/TTP).
    """
    if not delfin_iso or not ccdc_iso:
        return False
    if delfin_iso == ccdc_iso:
        return True
    if delfin_iso.startswith(f"CN{ccdc_cn}"):
        return True
    poly_words = {
        "octahedral": ("OH", "octahedral"),
        "tetrahedral": ("tet", "tetrahedral", "Tet"),
        "squareplanar": ("SP", "sp", "squareplanar"),
        "TBP": ("TBP",),
        "SPY": ("SPY",),
        "linear": ("linear",),
        "PB7": ("PBP", "PB7"),
        "SQAP8": ("Cube", "SQAP"),
        "TTP9": ("TTP",),
        "T-shape": ("T-shape", "trigonal"),
    }
    for key, words in poly_words.items():
        if key in ccdc_iso or any(w in ccdc_iso for w in words):
            for w in words + (key,):
                if w in delfin_iso:
                    return True
    return False


def compute_isomer_recall(
    emitted_by_smiles: Dict[str, list],
    family_table: Dict[str, dict],
    weighted: bool = False,
) -> float:
    """Global isomer-recall: mean over SMILES with >=1 family member.

    emitted_by_smiles : dict[smiles] = list[Path/str] of XYZ files DELFIN emitted
    family_table       : output of load_ccdc_family_table
    weighted           : if True, weight by family size

    Returns float in [0, 1].
    """
    totals = []
    weights = []
    for smi, xyz_paths in emitted_by_smiles.items():
        fam = family_table.get(smi)
        if not fam or fam["n_matches"] == 0:
            continue
        # CCDC-side: unique (isomer_label, cn) tuples
        ccdc_iso_tuples = set((m["isomer_label"], m.get("cn", -1))
                              for m in fam["matches"])
        if not ccdc_iso_tuples:
            continue
        delfin_isos = _emitted_isomer_labels(xyz_paths)
        matched = 0
        for ciso, ccn in ccdc_iso_tuples:
            if any(_isomer_match(d, ciso, ccn) for d in delfin_isos):
                matched += 1
        recall = matched / len(ccdc_iso_tuples)
        totals.append(recall)
        weights.append(fam["n_matches"])
    if not totals:
        return float("nan")
    if weighted:
        s = sum(t * w for t, w in zip(totals, weights))
        return s / sum(weights)
    return sum(totals) / len(totals)


def _cn_skeleton(iso_label: str) -> str:
    """Extract 'CN<N>' from a full isomer label."""
    m = re.match(r"(CN\d+)", iso_label or "")
    return m.group(1) if m else ""


def compute_conformer_recall(
    emitted_by_smiles: Dict[str, list],
    family_table: Dict[str, dict],
    tmc_index: Dict[str, dict],
    rmsd_threshold: float = DEFAULT_RMSD_A,
    weighted: bool = False,
    max_emitted_per_smiles: int = 50,
) -> float:
    """Global conformer-recall: for each CCDC family member, does SOME DELFIN
    emission match it within rmsd_threshold (heavy-atom Kabsch).

    Mean over SMILES with >=1 family member of:
       (# family-refcodes with at least one DELFIN-emission RMSD < threshold)
       / (# family-refcodes)
    """
    if np is None:
        return float("nan")
    totals = []
    weights = []
    for smi, xyz_paths in emitted_by_smiles.items():
        fam = family_table.get(smi)
        if not fam or fam["n_matches"] == 0:
            continue
        # Pre-parse emitted XYZs once per SMILES
        emitted_parsed = []
        for xp in xyz_paths[:max_emitted_per_smiles]:
            try:
                syms, P = parse_xyz(Path(xp).read_text())
                if syms and P is not None and len(P) > 0:
                    emitted_parsed.append((syms, P))
            except Exception:
                continue
        if not emitted_parsed:
            totals.append(0.0)
            weights.append(fam["n_matches"])
            continue
        n_family = len(fam["matches"])
        n_hit = 0
        for fam_rec in fam["matches"]:
            rc = fam_rec["refcode"]
            ccdc_rec = tmc_index.get(rc)
            if not ccdc_rec:
                continue
            ref_syms = ccdc_rec["symbols"]
            ref_P = np.asarray(ccdc_rec["positions"], dtype=float)
            best = float("inf")
            for syms, P in emitted_parsed:
                r = kabsch_rmsd_heavy(syms, P, ref_syms, ref_P)
                if math.isfinite(r) and r < best:
                    best = r
                    if best < rmsd_threshold:
                        break
            if math.isfinite(best) and best < rmsd_threshold:
                n_hit += 1
        recall = n_hit / max(n_family, 1)
        totals.append(recall)
        weights.append(n_family)
    if not totals:
        return float("nan")
    if weighted:
        s = sum(t * w for t, w in zip(totals, weights))
        return s / sum(weights)
    return sum(totals) / len(totals)


def per_smiles_recall_report(
    emitted_by_smiles: Dict[str, list],
    family_table: Dict[str, dict],
    tmc_index: Dict[str, dict],
    rmsd_threshold: float = DEFAULT_RMSD_A,
    max_emitted_per_smiles: int = 50,
) -> List[dict]:
    """Detailed per-SMILES report.

    Returns list of dicts (one per SMILES) with isomer_recall, conformer_recall,
    n_family, n_emitted, best_rmsd, ccdc_isomers_present, delfin_isomers_emitted.
    """
    if np is None:
        return []
    out = []
    for smi, xyz_paths in emitted_by_smiles.items():
        fam = family_table.get(smi)
        n_family = fam["n_matches"] if fam else 0
        if n_family == 0:
            out.append({
                "smiles": smi,
                "smiles_label": fam.get("smiles_label") if fam else None,
                "n_family": 0,
                "n_emitted": len(xyz_paths),
                "isomer_recall": float("nan"),
                "conformer_recall": float("nan"),
                "best_rmsd": float("nan"),
            })
            continue

        # Emit-side parse + classify
        emitted_parsed = []
        delfin_isos = set()
        for xp in xyz_paths[:max_emitted_per_smiles]:
            try:
                syms, P = parse_xyz(Path(xp).read_text())
                if syms and P is not None and len(P) > 0:
                    emitted_parsed.append((syms, P))
                    delfin_isos.add(classify_isomer(syms, P))
            except Exception:
                continue

        ccdc_iso_tuples = set((m["isomer_label"], m.get("cn", -1))
                              for m in fam["matches"])
        # isomer recall
        iso_matched = 0
        for ciso, ccn in ccdc_iso_tuples:
            if any(_isomer_match(d, ciso, ccn) for d in delfin_isos):
                iso_matched += 1
        iso_recall = iso_matched / max(len(ccdc_iso_tuples), 1)
        ccdc_isos = set(t[0] for t in ccdc_iso_tuples)

        # conformer recall
        n_hit = 0
        best_rmsd_global = float("inf")
        for fam_rec in fam["matches"]:
            rc = fam_rec["refcode"]
            ccdc_rec = tmc_index.get(rc)
            if not ccdc_rec:
                continue
            ref_syms = ccdc_rec["symbols"]
            ref_P = np.asarray(ccdc_rec["positions"], dtype=float)
            best = float("inf")
            for syms, P in emitted_parsed:
                r = kabsch_rmsd_heavy(syms, P, ref_syms, ref_P)
                if math.isfinite(r) and r < best:
                    best = r
            if math.isfinite(best):
                if best < best_rmsd_global:
                    best_rmsd_global = best
                if best < rmsd_threshold:
                    n_hit += 1
        conf_recall = n_hit / n_family

        out.append({
            "smiles": smi,
            "smiles_label": fam.get("smiles_label"),
            "n_family": n_family,
            "n_emitted": len(xyz_paths),
            "isomer_recall": iso_recall,
            "conformer_recall": conf_recall,
            "best_rmsd": (float(best_rmsd_global)
                          if math.isfinite(best_rmsd_global)
                          else float("nan")),
            "ccdc_isomer_classes": sorted(ccdc_isos)[:20],
            "delfin_isomers_emitted": sorted(delfin_isos)[:20],
        })
    return out


# ---------------------------------------------------------------------------
# Archive-scoring helper
# ---------------------------------------------------------------------------
def _smiles_from_filename(name: str) -> Optional[str]:
    """Try to extract SMILES from a pool XYZ filename if encoded; else None.

    Many DELFIN archives encode SMILES via JSON sidecar — when absent, this
    returns None and the caller must pass a mapping.
    """
    return None


def load_master_label_to_smiles(path: str) -> Dict[str, str]:
    """Load master SMILES pool (label|smiles per line) into a label->smiles map.

    Returns lookup-friendly variants: original label, label-with-special-chars-
    replaced-by-underscore (matches DELFIN's filename sanitization).
    """
    out: Dict[str, str] = {}
    p = Path(path)
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "|" not in line:
            continue
        label, smi = line.split("|", 1)
        label = label.strip()
        smi = smi.strip()
        out[label] = smi
        # Also sanitized version (paren/brackets/etc → underscore)
        sanitized = re.sub(r"[^A-Za-z0-9_\-]", "_", label)
        out[sanitized] = smi
    return out


def _smiles_from_xyz_header(xyz_path) -> Optional[str]:
    """Parse 'smi=<label>' or 'smiles=<value>' tokens from line 2 (comment) of XYZ."""
    try:
        with open(xyz_path) as fh:
            fh.readline()  # natoms
            comment = fh.readline()
    except Exception:
        return None
    if not comment:
        return None
    m = re.search(r"smi(?:les)?=([^\s]+(?:\([^)]*\))*[^\s]*)", comment)
    if m:
        return m.group(1)
    return None


def group_archive_by_smiles(archive_dir, smiles_mapping=None,
                            master_label_map: Optional[Dict[str, str]] = None
                            ) -> Dict[str, list]:
    """Group XYZ files in an archive by SMILES.

    Resolution order per file:
      1) Explicit smiles_mapping[stem] / [name]
      2) Sidecar JSON file <stem>.json with 'smiles' key
      3) XYZ comment line 'smi=<label>' resolved via master_label_map
      4) Filename stem itself resolved via master_label_map
    """
    archive = Path(archive_dir)
    by_smi: dict = defaultdict(list)
    for f in archive.glob("*.xyz"):
        smi = None
        if smiles_mapping:
            smi = smiles_mapping.get(f.stem) or smiles_mapping.get(f.name)
        if smi is None:
            sidecar = f.with_suffix(".json")
            if sidecar.exists():
                try:
                    sj = json.loads(sidecar.read_text())
                    smi = sj.get("smiles") or sj.get("input_smiles")
                except Exception:
                    pass
        if smi is None and master_label_map:
            # XYZ comment header
            label_in_header = _smiles_from_xyz_header(f)
            if label_in_header:
                smi = (master_label_map.get(label_in_header)
                       or master_label_map.get(re.sub(r"[^A-Za-z0-9_\-]", "_",
                                                       label_in_header)))
            if smi is None:
                # Filename stem
                smi = (master_label_map.get(f.stem)
                       or master_label_map.get(f.name))
        if smi is None:
            continue
        by_smi[smi].append(f)
    return by_smi


def score_archive(
    archive_dir,
    table_path: Optional[str] = None,
    index_path: Optional[str] = None,
    rmsd_threshold: Optional[float] = None,
    weighted: Optional[bool] = None,
    smiles_mapping=None,
    master_pool_path: Optional[str] = None,
) -> dict:
    """End-to-end scoring of a pool archive against the full-CCDC family table.

    Returns dict with global isomer_recall, conformer_recall, per-SMILES report
    and meta.
    """
    table = load_ccdc_family_table(table_path)
    index = load_ccdc_tmc_index(index_path)
    if not table:
        return {
            "error": ("no family table loaded — set DELFIN_FFFREE_XRD_RECALL_TABLE_PATH "
                      "or pass table_path")
        }
    if rmsd_threshold is None:
        rmsd_threshold = float(os.environ.get(ENV_RMSD_A, DEFAULT_RMSD_A))
    if weighted is None:
        weighted = bool(int(os.environ.get(ENV_WEIGHTED, "0")))

    master_map = None
    if master_pool_path is None:
        master_pool_path = os.environ.get(
            "DELFIN_FFFREE_XRD_RECALL_MASTER_POOL",
            "/home/qmchem_max/agent_workspace/quality_framework/pools/"
            "smiles_master_v3_plus.txt",
        )
    if master_pool_path:
        master_map = load_master_label_to_smiles(master_pool_path)
    by_smi = group_archive_by_smiles(archive_dir,
                                     smiles_mapping=smiles_mapping,
                                     master_label_map=master_map)
    if not by_smi:
        return {
            "archive": str(archive_dir),
            "n_smiles_with_emissions": 0,
            "note": "no SMILES mapping resolved — provide smiles_mapping or sidecar JSON",
        }

    iso = compute_isomer_recall(by_smi, table, weighted=weighted)
    conf = compute_conformer_recall(
        by_smi, table, index, rmsd_threshold=rmsd_threshold, weighted=weighted,
    )
    n_with_family = sum(1 for smi in by_smi if table.get(smi, {}).get("n_matches", 0) > 0)
    return {
        "archive": str(archive_dir),
        "isomer_recall": iso,
        "conformer_recall": conf,
        "n_smiles_total_in_archive": len(by_smi),
        "n_smiles_with_ccdc_family": n_with_family,
        "rmsd_threshold_A": rmsd_threshold,
        "weighted": weighted,
        "n_family_total": sum(table[smi]["n_matches"]
                              for smi in by_smi if smi in table),
    }
