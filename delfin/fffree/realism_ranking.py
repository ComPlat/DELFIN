"""XYZ realism ranking — realistic-first output sort.

Purpose
-------
DFT/xTB downstream consumers typically take the *first* emitted frame of a
multi-trajectory XYZ as their starting point.  This module composes ALL
available realism signals (Mogul anomaly count, aggregate CShM, inter- and
intra-ligand clash counts, GRIP final loss, Pólya completeness, Burnside
coverage) into ONE composite score per emitted XYZ frame and produces a
deterministic ranking with the most realistic frame at rank-0.

Design
------
1. Hard-filter gates (topology + build-time clash) come first — failing
   structures are pushed to the end of the rank list (never excluded —
   that would change emission count and break downstream invariants).
2. Soft signals are rank-normalised within the per-SMILES group so a single
   pathological outlier does not dominate the score.
3. Composite ``realism = Σ wᵢ · normalised_signalᵢ``; lower = more realistic.
4. Stable, deterministic sort: ``(realism_score, frame_index)``.

Determinism
-----------
All sort keys derived from immutable record fields.  ``PYTHONHASHSEED=0``
respected.  Two runs over the same input → byte-identical ``rank_xyz_group``
output (verified by ``tests/test_realism_ranking.py``).

Env-flags
---------
* ``DELFIN_FFFREE_REALISM_SORT=1``  — master enable (default OFF,
  byte-identical at module import time and inside pipeline)
* ``DELFIN_FFFREE_REALISM_WEIGHT_MOGUL=0.30``  — per-signal weight overrides
* ``DELFIN_FFFREE_REALISM_WEIGHT_CSHM=0.20``
* ``DELFIN_FFFREE_REALISM_WEIGHT_INTER_CLASH=0.15``
* ``DELFIN_FFFREE_REALISM_WEIGHT_HH_CLASH=0.10``
* ``DELFIN_FFFREE_REALISM_WEIGHT_GRIP_LOSS=0.10``
* ``DELFIN_FFFREE_REALISM_WEIGHT_POLYA=0.05``
* ``DELFIN_FFFREE_REALISM_WEIGHT_BURNSIDE=0.05``

API
---
* :func:`compute_realism_score`  — pure function, frame-level score
* :func:`rank_xyz_group`         — deterministic sort of a SMILES-group
* :func:`rank_archive`           — batch over an archive directory
* :func:`realism_sort_active`    — env-gate check
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Env-gate
# ---------------------------------------------------------------------------
MASTER_ENV = "DELFIN_FFFREE_REALISM_SORT"

# Default weights (sum = 1.0); see module docstring.
# The numeric mix follows the user-supplied 2026-06-04 directive (Mogul
# 0.30 / CShM 0.20 / inter-clash 0.15 / H-H 0.10 / GRIP-loss 0.10 / Pólya
# 0.05 / Burnside 0.05).  CShM is bumped from the indicative 0.20 to 0.25
# so the soft weights — the only ones contributing to the score — sum to
# 1.0 (hard gates carry no weight; they apply as binary post-filter).
_DEFAULT_WEIGHTS: Dict[str, float] = {
    "mogul": 0.30,
    "cshm": 0.25,
    "inter_clash": 0.15,
    "hh_clash": 0.10,
    "grip_loss": 0.10,
    "polya": 0.05,
    "burnside": 0.05,
}

# Hard-filter gates push failing frames to the back.  They do not get a
# weight — they are applied AFTER soft scoring.
_HARD_GATES = ("topology_ok", "build_time_clash_ok")


def realism_sort_active() -> bool:
    """Return True when the master env-flag is set to ``"1"``."""
    return os.environ.get(MASTER_ENV, "") == "1"


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "")
    if not raw:
        return float(default)
    try:
        return float(raw)
    except ValueError:
        return float(default)


def get_weights() -> Dict[str, float]:
    """Return the active weight map, honouring per-signal env overrides.

    Always returns a copy.  The returned mapping is normalised so the sum
    equals 1.0 — keeps the score scale stable when callers tweak weights
    one at a time.
    """
    w: Dict[str, float] = {}
    for key, default in _DEFAULT_WEIGHTS.items():
        env_name = f"DELFIN_FFFREE_REALISM_WEIGHT_{key.upper()}"
        w[key] = max(0.0, _env_float(env_name, default))
    total = sum(w.values())
    if total <= 0.0:
        # Fall back to defaults silently if all weights zeroed out.
        return dict(_DEFAULT_WEIGHTS)
    return {k: v / total for k, v in w.items()}


# ---------------------------------------------------------------------------
# Signal extraction
# ---------------------------------------------------------------------------
def _safe_float(x, default: float = float("nan")) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return default
    if math.isnan(v) or math.isinf(v):
        return default
    return v


def extract_signals(metadata: Dict) -> Dict[str, float]:
    """Pull the seven soft-signal numbers + two hard-gate booleans from a
    per-frame metadata record.

    Accepts a flat dict that may contain any subset of the canonical keys:

    * ``mogul_anomaly_count``  (lower = better)
    * ``cshm``                  (lower = better)
    * ``inter_ligand_clash_count``  (lower = better)
    * ``hh_clash_count``       (lower = better)
    * ``grip_final_loss``      (lower = better)
    * ``polya_complete``       (bool / 0-1: True = better)
    * ``burnside_coverage``    (0-1: higher = better; internally inverted)
    * ``topology_ok``          (bool — hard gate)
    * ``build_time_clash_ok``  (bool — hard gate)

    Missing values yield ``nan`` — handled by the normaliser as
    "neutral" (median rank).  Hard gates default to ``True`` (pass) when
    absent so legacy archives without these fields are not penalised.

    Returns a dict with the seven soft signals (already oriented so smaller
    = more realistic — burnside is inverted as ``1 - x``) plus the two
    hard-gate booleans.
    """
    out: Dict[str, float] = {}
    out["mogul"] = _safe_float(metadata.get("mogul_anomaly_count"))
    out["cshm"] = _safe_float(metadata.get("cshm"))
    out["inter_clash"] = _safe_float(metadata.get("inter_ligand_clash_count"))
    out["hh_clash"] = _safe_float(metadata.get("hh_clash_count"))
    out["grip_loss"] = _safe_float(metadata.get("grip_final_loss"))

    polya = metadata.get("polya_complete")
    if polya is None:
        out["polya"] = float("nan")
    else:
        # True / "full" / nonzero → 0 (good); False → 1 (bad).
        if isinstance(polya, str):
            polya_val = 1.0 if polya.lower() in ("true", "full", "complete", "1") else 0.0
        else:
            polya_val = 1.0 if bool(polya) else 0.0
        out["polya"] = 1.0 - polya_val  # invert: 0 = best

    burnside = metadata.get("burnside_coverage")
    if burnside is None:
        out["burnside"] = float("nan")
    else:
        bv = _safe_float(burnside)
        if math.isnan(bv):
            out["burnside"] = float("nan")
        else:
            # 0-1 coverage; invert so lower = better.
            out["burnside"] = 1.0 - max(0.0, min(1.0, bv))

    # Hard gates default True (pass) when absent.
    out["topology_ok"] = bool(metadata.get("topology_ok", True))
    out["build_time_clash_ok"] = bool(metadata.get("build_time_clash_ok", True))

    return out


def _rank_normalise(values: Sequence[float]) -> List[float]:
    """Map ``values`` to rank-percentile in ``[0, 1]`` (lower = lower rank).

    NaN values get the median rank (0.5) — treated as "neutral" so a
    missing signal does not penalise or favour the frame.  Ties broken
    by stable sort on original index.
    """
    n = len(values)
    if n == 0:
        return []
    if n == 1:
        return [0.5]
    # Separate finite from NaN; rank finites only.
    finite_with_idx = [(i, v) for i, v in enumerate(values) if not math.isnan(v)]
    if not finite_with_idx:
        return [0.5] * n
    # Sort finite by value (stable on idx).
    finite_with_idx.sort(key=lambda t: (t[1], t[0]))
    rank_map: Dict[int, float] = {}
    m = len(finite_with_idx)
    if m == 1:
        rank_map[finite_with_idx[0][0]] = 0.5
    else:
        for rank, (i, _v) in enumerate(finite_with_idx):
            rank_map[i] = rank / (m - 1)
    return [rank_map.get(i, 0.5) for i in range(n)]


def compute_realism_score(
    metadata: Dict,
    *,
    weights: Optional[Dict[str, float]] = None,
    group_signals: Optional[Dict[str, List[float]]] = None,
    frame_index: int = 0,
) -> float:
    """Compute the composite realism score for a single frame.

    Lower score = more realistic.  Hard-gate failures add a large constant
    penalty (1e3) so failing frames always rank after passing frames.

    Parameters
    ----------
    metadata
        Per-frame metric dict (see :func:`extract_signals`).
    weights
        Optional override of weight map; defaults to :func:`get_weights`.
    group_signals
        If supplied, used for rank-normalisation against the group.  When
        ``None`` the score is computed in "isolated" mode where each
        signal is mapped to ``[0,1]`` via a sigmoid; this is a degraded
        ranking — prefer :func:`rank_xyz_group` whenever possible.
    frame_index
        Index of this frame inside ``group_signals`` — required iff
        ``group_signals`` is provided.

    Returns
    -------
    float
        Composite score.  ``+inf`` only when the input is structurally
        broken (no usable signals AND hard-gate failure).
    """
    w = weights if weights is not None else get_weights()
    signals = extract_signals(metadata)

    penalty = 0.0
    for gate in _HARD_GATES:
        if not signals.get(gate, True):
            penalty += 1e3

    soft_keys = list(_DEFAULT_WEIGHTS.keys())
    score = 0.0
    if group_signals is not None:
        for key in soft_keys:
            vec = group_signals.get(key, [])
            if not vec or frame_index >= len(vec):
                # Fall through with neutral rank (0.5) so weight not lost.
                norm = 0.5
            else:
                # group_signals[key] is already pre-normalised by caller.
                norm = vec[frame_index]
            score += w.get(key, 0.0) * norm
    else:
        # Isolated mode — use a stable monotone squash for absent group
        # context.  Used only when caller cannot supply group context.
        for key in soft_keys:
            raw = signals.get(key, float("nan"))
            if math.isnan(raw):
                norm = 0.5
            else:
                norm = 1.0 / (1.0 + math.exp(-raw / 5.0))
            score += w.get(key, 0.0) * norm

    return float(score + penalty)


def rank_xyz_group(
    frame_metadata: Sequence[Dict],
    *,
    weights: Optional[Dict[str, float]] = None,
) -> List[Tuple[int, float]]:
    """Rank frames of a single SMILES-group by composite realism (best first).

    Parameters
    ----------
    frame_metadata
        Ordered list of per-frame metadata dicts.  Order corresponds to
        the original frame index in the multi-trajectory XYZ.
    weights
        Optional weight override.

    Returns
    -------
    list of (orig_frame_index, score)
        Sorted ascending by score (lowest = most realistic = rank-0).
        Ties broken deterministically by original frame index.

    Edge cases
    ----------
    * empty input → ``[]``
    * single frame → ``[(0, score)]``
    * all-missing signals → preserves original order (alphabetical /
      insertion order)
    """
    n = len(frame_metadata)
    if n == 0:
        return []
    if n == 1:
        s = compute_realism_score(frame_metadata[0], weights=weights)
        return [(0, s)]

    w = weights if weights is not None else get_weights()
    soft_keys = list(_DEFAULT_WEIGHTS.keys())

    # Collect raw signals per frame.
    raw_per_key: Dict[str, List[float]] = {k: [] for k in soft_keys}
    hard_flags: List[bool] = []
    for md in frame_metadata:
        sig = extract_signals(md)
        for k in soft_keys:
            raw_per_key[k].append(sig.get(k, float("nan")))
        ok = all(sig.get(g, True) for g in _HARD_GATES)
        hard_flags.append(ok)

    # Rank-normalise each signal within the group.
    norm_per_key: Dict[str, List[float]] = {
        k: _rank_normalise(raw_per_key[k]) for k in soft_keys
    }

    # Compute composite per frame.
    out: List[Tuple[int, float]] = []
    for i in range(n):
        soft = sum(w.get(k, 0.0) * norm_per_key[k][i] for k in soft_keys)
        penalty = 0.0 if hard_flags[i] else 1e3
        out.append((i, float(soft + penalty)))

    # Stable deterministic sort: score then original frame index.
    out.sort(key=lambda t: (t[1], t[0]))
    return out


# ---------------------------------------------------------------------------
# JSONL aggregation — pull per-frame metric records from run_all_detectors
# output into one frame_metadata dict per (SMILES, frame_index).
# ---------------------------------------------------------------------------
def _split_header(header: str) -> Tuple[str, str]:
    """Extract ``(smiles_id, label)`` from a detector-JSONL ``header`` field.

    Header convention: ``commit=<X> smi=<SID> label=<LABEL>`` — fields
    separated by single spaces.  Missing fields yield empty strings.
    """
    smi = ""
    label = ""
    for tok in header.split():
        if tok.startswith("smi="):
            smi = tok[4:]
        elif tok.startswith("label="):
            label = tok[6:]
    # When label spans multiple words the simple split above misses the
    # tail.  Recover by slicing after "label=".
    li = header.find("label=")
    if li >= 0:
        label = header[li + 6 :]
    return smi, label


# Map JSONL filename suffix → canonical signal key (and how to extract the
# numeric value from the record).  This adapter layer means new detector
# JSONLs can be incorporated by adding one line.
_JSONL_ADAPTERS: Dict[str, Tuple[str, str]] = {
    # suffix → (signal_key, record_field)
    "mogul_v3": ("mogul_anomaly_count", "n_anomalies"),
    "coordgeom": ("cshm", "mean_abs_dev"),
    "severe_overlap": ("inter_ligand_clash_count", "n_severe_overlaps"),
    "hclash": ("hh_clash_count", "n_clash"),
    "hanomaly": ("hh_clash_count", "n_anomaly"),
    "ligcollapse": ("build_time_clash_ok", "ok"),
    "isocoverage": ("polya_complete", "coverage_pct"),
}


def _iter_jsonl(path: Path) -> Iterable[Dict]:
    if not path.exists():
        return
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def _sanitise_smiles_id(smi: str) -> str:
    """Apply the standard ``re.sub(r"[^A-Za-z0-9_\\-]", "_", sid)`` used by
    ``pool_evaluator`` so XYZ-filename stems can be matched against the
    raw SMILES ids carried in detector-JSONL records.
    """
    import re
    return re.sub(r"[^A-Za-z0-9_\-]", "_", smi)[:120]


def load_metadata_from_archive(
    archive_dir: Path,
    metric_jsonls: Sequence[Path],
) -> Dict[Tuple[str, str], Dict]:
    """Aggregate detector-JSONL records into a per-(SMILES, frame-label) dict.

    Returns
    -------
    dict
        Keyed by ``(smiles_id, label)`` — the same identifiers used in the
        XYZ-archive trajectory.  Each value is a flat metadata dict
        suitable for :func:`extract_signals`.

    Notes
    -----
    Robust against missing JSONLs / missing fields — those simply do not
    contribute a signal.  The downstream normaliser treats absences as
    neutral (rank 0.5).
    """
    metadata: Dict[Tuple[str, str], Dict] = defaultdict(dict)
    for p in metric_jsonls:
        suffix = None
        for s in _JSONL_ADAPTERS:
            if p.name.endswith(f"_{s}.jsonl"):
                suffix = s
                break
        if suffix is None:
            continue
        signal_key, field = _JSONL_ADAPTERS[suffix]
        for rec in _iter_jsonl(p):
            smi = rec.get("smiles_id") or ""
            label = ""
            if "header" in rec:
                smi, label = _split_header(rec["header"])
            if not smi:
                # Fall back to file stem (e.g. "01-Fe_CO_3_NHC_2.xyz").
                f = rec.get("file", "")
                if f:
                    smi = Path(f).stem
            if not smi:
                continue
            # Index under BOTH the raw smi and its sanitised form so
            # callers can look up via either representation.
            smi_sanitised = _sanitise_smiles_id(smi)
            key = (smi, label)
            val = rec.get(field)
            # Special-case polya_complete from coverage_pct (≥80 = complete).
            if signal_key == "polya_complete" and val is not None:
                try:
                    val = float(val) >= 80.0
                except (TypeError, ValueError):
                    val = None
            # Special-case build_time_clash_ok from ligcollapse 'ok' flag.
            if signal_key == "build_time_clash_ok" and val is None:
                # ligcollapse uses 'ok' or 'collapse_count'
                cc = rec.get("collapse_count")
                if cc is not None:
                    val = (cc == 0)
            if val is not None:
                metadata[key][signal_key] = val
                # Also store under the sanitised SMILES key when it
                # differs (e.g. "01-Fe(CO)3(NHC)2" vs file stem
                # "01-Fe_CO_3_NHC_2").
                if smi_sanitised != smi:
                    metadata[(smi_sanitised, label)][signal_key] = val
    return dict(metadata)


# ---------------------------------------------------------------------------
# Multi-frame XYZ I/O
# ---------------------------------------------------------------------------
def _read_multixyz(path: Path) -> List[Tuple[List[str], str]]:
    """Read a multi-frame XYZ.

    Returns a list of ``(lines, label)`` per frame where ``lines`` is the
    full block (atom-count line + comment + atoms) and ``label`` is the
    value of the ``label=`` token in the comment (empty string if absent).
    """
    out: List[Tuple[List[str], str]] = []
    if not path.exists():
        return out
    with path.open("r") as f:
        text = f.read()
    if not text.strip():
        return out
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
        block = lines[i : i + natoms + 2]
        comment = block[1] if len(block) > 1 else ""
        # Extract label
        label = ""
        li = comment.find("label=")
        if li >= 0:
            label = comment[li + 6 :].strip()
        out.append((block, label))
        i += natoms + 2
    return out


def _write_multixyz(path: Path, frames: Sequence[Tuple[List[str], str]]) -> None:
    """Write frames in order to ``path``."""
    with path.open("w") as f:
        for i, (block, _label) in enumerate(frames):
            for j, line in enumerate(block):
                f.write(line)
                if j < len(block) - 1 or i < len(frames) - 1:
                    f.write("\n")
        f.write("\n")


# ---------------------------------------------------------------------------
# Archive-level batch ranking
# ---------------------------------------------------------------------------
def rank_archive(
    archive_dir: Path,
    output_jsonl: Path,
    *,
    metric_jsonls: Optional[Sequence[Path]] = None,
    weights: Optional[Dict[str, float]] = None,
    rewrite_xyz: bool = False,
) -> Dict[str, int]:
    """Rank all SMILES in ``archive_dir`` by realism, write a sidecar JSONL.

    Parameters
    ----------
    archive_dir
        Directory containing one multi-frame XYZ per SMILES.
    output_jsonl
        Sidecar JSONL — one record per SMILES, listing frames in
        realism-rank order with the score and the original frame label.
    metric_jsonls
        Detector JSONLs to pull signals from.  When ``None`` the function
        auto-discovers files named ``<archive_name>_<suffix>.jsonl`` in
        the archive directory's PARENT (the standard layout).
    weights
        Optional weight override.
    rewrite_xyz
        When True, rewrite the XYZ files with frames in rank order (and
        store the original ordering in a ``.orig_order.json`` sidecar so
        the operation is reversible).  Default False — we prefer the
        sidecar-JSONL approach so emission count and existing tooling
        continue to work.

    Returns
    -------
    dict
        ``{"n_files": <int>, "n_frames_total": <int>, "n_ranked": <int>}``
    """
    archive_dir = Path(archive_dir)
    output_jsonl = Path(output_jsonl)
    output_jsonl.parent.mkdir(parents=True, exist_ok=True)

    # Auto-discover metric JSONLs (sibling files with the archive name as
    # prefix).
    if metric_jsonls is None:
        parent = archive_dir.parent
        prefix = archive_dir.name + "_"
        metric_jsonls = sorted(parent.glob(f"{prefix}*.jsonl"))

    metadata_map = load_metadata_from_archive(archive_dir, list(metric_jsonls))

    xyz_files = sorted(archive_dir.glob("*.xyz"))
    n_frames_total = 0
    n_ranked = 0
    with output_jsonl.open("w") as fout:
        for xp in xyz_files:
            frames = _read_multixyz(xp)
            if not frames:
                continue
            n_frames_total += len(frames)

            # Build per-frame metadata list — match on (smiles_id, label).
            smiles_id = xp.stem
            # Strip leading "NN-" / "NNN-" prefix used in many archives.
            # We keep both candidates and try both.
            candidates = {smiles_id}
            if "-" in smiles_id and smiles_id.split("-", 1)[0].isdigit():
                candidates.add(smiles_id.split("-", 1)[1])
            # Also try replacing underscores with original SMILES chars
            # (best-effort — fall back to label-only match below).
            frame_md: List[Dict] = []
            for _block, label in frames:
                md: Dict = {}
                for cand in candidates:
                    key = (cand, label)
                    if key in metadata_map:
                        md = dict(metadata_map[key])
                        break
                frame_md.append(md)

            ranking = rank_xyz_group(frame_md, weights=weights)
            n_ranked += len(ranking)

            record = {
                "smiles_id": smiles_id,
                "xyz_file": str(xp),
                "n_frames": len(frames),
                "weights": weights if weights is not None else get_weights(),
                "ranking": [
                    {
                        "rank": r,
                        "orig_frame_index": int(orig),
                        "label": frames[orig][1],
                        "score": float(score),
                    }
                    for r, (orig, score) in enumerate(ranking)
                ],
            }
            fout.write(json.dumps(record) + "\n")

            if rewrite_xyz and ranking:
                # Save original order before rewriting.
                orig_path = xp.with_suffix(".orig_order.json")
                if not orig_path.exists():
                    orig_path.write_text(
                        json.dumps(
                            {"orig_labels": [lbl for _b, lbl in frames]},
                            indent=2,
                        )
                    )
                reordered = [frames[orig] for orig, _ in ranking]
                _write_multixyz(xp, reordered)

    return {
        "n_files": len(xyz_files),
        "n_frames_total": n_frames_total,
        "n_ranked": n_ranked,
    }


__all__ = [
    "MASTER_ENV",
    "realism_sort_active",
    "get_weights",
    "extract_signals",
    "compute_realism_score",
    "rank_xyz_group",
    "load_metadata_from_archive",
    "rank_archive",
]
