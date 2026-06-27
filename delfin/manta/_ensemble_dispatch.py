"""Phase 1.5 ensemble dispatcher — combine K conversion paths per SMILES.

Each path is a configuration of DELFIN env-flags that activates a specific
FPCFD champion-port mechanism (Iter-8.x or future iters).  Running multiple
paths in parallel and pooling their candidate frames lets the ensemble
recover per-class wins without committing any single path as HEAD's
default.  No path's env-flag is required to be on at HEAD level: the
ensemble adds them at runtime.

v1 → v2 changelog
-----------------
v2 introduces:

    * Per-frame composite scorer (``_ensemble_scoring.score_frame``):
      topology fidelity + M-D invariant + clashes + tridentate planarity +
      bond-length outliers.  Self-contained — no quality_framework runtime
      dependency.

    * Top-K per fingerprint (default K=2): frames are grouped by canonical-
      form polyhedron+donor-element fingerprint (mirrors
      find_isomer_coverage_named); within each group only the K best-scoring
      frames survive.  Result: each isomer-name receives its top-K conformers
      from whichever path produced them.

    * Class-aware path subset: dispatcher picks the relevant subset based on
      ``smiles_converter._classify_complex_class(mol)``.  Saves time and
      reduces noise.  ``paths_mode="all"`` overrides for forensik runs.

    * Adaptive-memory hook (Phase 4 prep): every dispatch call writes a JSONL
      line to ``~/agent_workspace/quality_framework/results/ensemble_dispatch_log.jsonl``
      capturing class, paths_run, frames_per_path, scores, winning paths and
      n_unique_fingerprints.  Phase 4 KG learner consumes this log later.

    * Mutable path catalogue: ``register_path()`` lets a future self-loop add
      new paths at runtime without restarting the process.

v2 preserves bit-exact backward-compat with v1:
    * Old positional/keyword callers still work; the new fields default to
      values that produce v1-equivalent output (same Tuple shape; score is
      now populated instead of None but the field position is preserved).
    * Path env-overrides are still applied via try/finally snapshot/restore.
    * Sequential in-process execution (smiles_to_xyz_isomers uses module
      globals — multi-threading is NOT safe).

API
---
    from delfin.manta._ensemble_dispatch import ensemble_smiles_to_xyz_isomers

    out, errors = ensemble_smiles_to_xyz_isomers(
        smiles="...",
        paths=None,                      # default: class-aware subset
        max_isomers_per_path=5,
        max_isomers_total=20,
        top_k_per_fingerprint=2,         # NEW v2
        score_weights=None,              # NEW v2 — overrides DEFAULT_WEIGHTS
        log_dispatch=True,               # NEW v2 — JSONL log on/off
    )
    # out = [(xyz, label, score, source_path), ...]
    # errors = [str, ...]                — per-path failure messages

Catalogue extensibility
-----------------------
Add a new path at runtime via::

    from delfin.manta._ensemble_dispatch import register_path, Path
    register_path(Path("iter9.0a-newport", {"DELFIN_ITER9_NEWFLAG": "1"}))

The next ``ensemble_smiles_to_xyz_isomers(..., paths_mode="all")`` call will
include it; the class-aware default mode requires also tagging
``classes={"sigma", ...}`` so the dispatcher knows when to dispatch it.
"""
from __future__ import annotations

import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path as FsPath
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------- #
# Path catalogue
# ---------------------------------------------------------------------- #

@dataclass(frozen=True)
class Path:
    """One conversion path = name + env-var overrides applied at runtime.

    Attributes
    ----------
    name
        Unique identifier (used as ``source_path`` provenance tag in output).
    env_overrides
        Mapping of env-var → string value applied around the converter call.
        Empty dict = "use HEAD baseline as-is".
    classes
        Frozenset of chemistry classes for which this path is relevant
        (members of ``{"no_metal","sigma","hapto","multi_sigma","multi_hapto"}``).
        ``frozenset()`` = "all classes" (default for HEAD baseline).
    description
        Human-readable note for forensik / log readers.
    requires_implementation
        ``True`` if the path is a placeholder for a future env-flag-only
        port (e.g. 5687b3d champion has not yet been ported to a runtime
        env-flag).  Such paths are SKIPPED at dispatch time but kept in the
        catalogue so the self-loop knows what mechanisms remain to be
        implemented.
    """

    name: str
    env_overrides: Dict[str, str] = field(default_factory=dict)
    classes: FrozenSet[str] = field(default_factory=frozenset)
    description: str = ""
    requires_implementation: bool = False


# Default catalogue.  Order matters: HEAD baseline first so it wins ties on
# fingerprint grouping (deterministic behaviour matches user expectation that
# HEAD is preferred unless an iter-port produces measurably better geometry).
DEFAULT_PATHS: List[Path] = [
    Path(
        name="head-default",
        env_overrides={},
        classes=frozenset(),  # all classes
        description="Current HEAD with all defaults — no env overrides.",
    ),
    Path(
        name="iter8.4abc-on",
        env_overrides={
            "DELFIN_SIGMA_PORT_123A130_ITER8": "1",
            "DELFIN_SIGMA_SKIP_PUCKER_ITER8":  "1",
            "DELFIN_SIGMA_TIGHT_THRESHOLD_ITER8": "1",
        },
        classes=frozenset({"sigma"}),
        description="Iter-8.4 sigma port (123a130 chelate-cap restoration).",
    ),
    Path(
        name="iter8.5a-hapto-reject",
        env_overrides={
            "DELFIN_ITER85_HARD_REJECT_CLASSES": "hapto,multi_hapto,multi_sigma",
        },
        classes=frozenset({"hapto", "multi_hapto", "multi_sigma"}),
        description="Iter-8.5a hard-reject filter for hapto-class extras tail.",
    ),
    Path(
        name="iter8.6a-hapto-port-off",
        env_overrides={"DELFIN_HAPTO_123A_PORT": "0"},
        classes=frozenset({"hapto", "multi_hapto"}),
        description="Iter-8.6a — disable 123a port, restore hapto-default σ-guard.",
    ),
    Path(
        name="iter8.5a-pump-skip-broad",
        env_overrides={
            "DELFIN_ITER85_PUMP_SKIP_CLASSES": "hapto,multi_hapto,multi_sigma",
        },
        classes=frozenset({"hapto", "multi_hapto", "multi_sigma"}),
        description="Iter-8.5a — broaden pump-skip across all hapto-likely classes.",
    ),
    Path(
        name="iter8.3-multisigma",
        env_overrides={"DELFIN_MULTISIGMA_PORT_5B3E0D2_ITER8": "1"},
        classes=frozenset({"multi_sigma"}),
        description="Iter-8.3 multi-sigma 5b3e0d2 port.",
    ),
    Path(
        name="iter8.1-multihapto-filter",
        env_overrides={"DELFIN_ITER81_FILTER": "1"},
        classes=frozenset({"multi_hapto", "hapto"}),
        description="Iter-8.1 multi-hapto safe-fallback filter (default-on).",
    ),
    # ----- Champion placeholder slots (future iters fill these) ----- #
    Path(
        name="champion-123a130",
        env_overrides={},
        classes=frozenset({"multi_sigma", "sigma"}),
        description="123a130 champion — env-flag-only port pending.",
        requires_implementation=True,
    ),
    Path(
        name="champion-5687b3d",
        env_overrides={},
        classes=frozenset({"hapto", "multi_hapto"}),
        description="5687b3d champion (H realism + topology) — port pending.",
        requires_implementation=True,
    ),
    Path(
        name="champion-e6761e4",
        env_overrides={},
        classes=frozenset({"multi_sigma", "sigma"}),
        description="e6761e4 champion (M-D invariant + sp3-C donor) — port pending.",
        requires_implementation=True,
    ),
    Path(
        name="champion-cf1d480",
        env_overrides={},
        classes=frozenset({"hapto", "multi_hapto"}),
        description="cf1d480 champion (hapto extras) — port pending.",
        requires_implementation=True,
    ),
    Path(
        name="champion-5b3e0d2",
        env_overrides={},
        classes=frozenset({"multi_sigma"}),
        description="5b3e0d2 champion (multi-sigma geom) — port pending.",
        requires_implementation=True,
    ),
    Path(
        name="champion-81f8a1f",
        env_overrides={},
        classes=frozenset({"hapto"}),
        description="81f8a1f champion (H realism) — port pending.",
        requires_implementation=True,
    ),
]


# Class-aware default selection.  Updated whenever ``register_path`` is called
# so the self-loop sees runtime additions.  Maps class → path-name list.
_CLASS_PATH_MAP: Dict[str, List[str]] = {
    "no_metal":     ["head-default"],
    "sigma":        ["head-default", "iter8.4abc-on"],
    "hapto":        ["head-default", "iter8.6a-hapto-port-off",
                     "iter8.5a-hapto-reject", "iter8.5a-pump-skip-broad"],
    "multi_sigma":  ["head-default", "iter8.3-multisigma",
                     "iter8.5a-hapto-reject"],
    "multi_hapto":  ["head-default", "iter8.6a-hapto-port-off",
                     "iter8.1-multihapto-filter",
                     "iter8.5a-pump-skip-broad"],
}


def register_path(path: Path,
                  class_path_map_update: Optional[Dict[str, List[str]]] = None
                  ) -> None:
    """Append a path to ``DEFAULT_PATHS`` and (optionally) tag it for class-
    aware dispatch.  Self-loop / runtime-mutation safe.

    ``class_path_map_update`` is a ``{class: [path_name, ...], ...}`` mapping
    that is merged into ``_CLASS_PATH_MAP`` (existing class-lists are EXTENDED
    not replaced).  If omitted, the path is reachable only via
    ``paths_mode="all"`` or by passing it explicitly.
    """
    # Reject duplicate names
    if any(p.name == path.name for p in DEFAULT_PATHS):
        raise ValueError(f"Path with name '{path.name}' already registered")
    DEFAULT_PATHS.append(path)
    if class_path_map_update:
        for cls, names in class_path_map_update.items():
            existing = _CLASS_PATH_MAP.setdefault(cls, [])
            for n in names:
                if n not in existing:
                    existing.append(n)


def select_paths_for_class(
    chemistry_class: str,
    include_placeholders: bool = False,
) -> List[Path]:
    """Return the catalogue subset relevant for a given chemistry class.

    Falls back to ``[head-default]`` when the class is unknown.  Skips paths
    flagged ``requires_implementation`` unless ``include_placeholders=True``.
    """
    name_list = _CLASS_PATH_MAP.get(chemistry_class, ["head-default"])
    by_name = {p.name: p for p in DEFAULT_PATHS}
    out: List[Path] = []
    for n in name_list:
        p = by_name.get(n)
        if p is None:
            continue
        if p.requires_implementation and not include_placeholders:
            continue
        out.append(p)
    if not out:
        # Defensive: head-default must always be reachable.
        head = by_name.get("head-default")
        if head is not None:
            out.append(head)
    return out


# ---------------------------------------------------------------------- #
# Helpers — fingerprinting + label canonicalisation
# ---------------------------------------------------------------------- #

def _isomer_label_root(label: str) -> str:
    """Extract the isomer-fingerprint-equivalent label root (mirrors v1)."""
    if not label:
        return ""
    label = label.strip()
    label = re.sub(r"-conf\d+$", "", label)
    if " " in label:
        prefix, types = label.split(" ", 1)
        toks = sorted(t.strip() for t in types.split("|") if t.strip())
        return prefix.strip() + " " + "|".join(toks)
    return label


# ---------------------------------------------------------------------- #
# Adaptive-memory log
# ---------------------------------------------------------------------- #

_DISPATCH_LOG_PATH = (
    FsPath.home() / "agent_workspace" / "quality_framework" / "results"
    / "ensemble_dispatch_log.jsonl"
)


def _write_dispatch_log(record: Dict) -> None:
    """Best-effort JSONL append.  Silently swallows any IO error so the
    converter is never broken by missing-directory / permission issues.
    Phase 4 KG learner is the sole consumer; not on the hot path.
    """
    try:
        _DISPATCH_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with _DISPATCH_LOG_PATH.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(record, separators=(",", ":")) + "\n")
    except Exception:
        pass


# ---------------------------------------------------------------------- #
# Main entry
# ---------------------------------------------------------------------- #

def ensemble_smiles_to_xyz_isomers(
    smiles: str,
    paths: Optional[List[Path]] = None,
    max_isomers_per_path: int = 5,
    max_isomers_total: int = 20,
    top_k_per_fingerprint: int = 2,
    score_weights: Optional[Dict[str, float]] = None,
    log_dispatch: bool = True,
    quality_mode: Optional[str] = "fast",
    deduplicate: bool = True,
    paths_mode: str = "class_aware",  # "class_aware" | "all" | "explicit"
    include_placeholders: bool = False,
    **converter_kwargs,
) -> Tuple[List[Tuple[str, str, float, str]], List[str]]:
    """Run K conversion paths and return their pooled top-K-per-fingerprint frames.

    Returns
    -------
    output : List[Tuple[xyz, label, score, source_path]]
        ``score`` is the composite cost (lower = better) computed by
        ``_ensemble_scoring.score_frame``.  When the scorer cannot run
        (RDKit missing, XYZ unparseable etc.) the score falls back to
        ``float('inf')``; such frames remain in the output but always lose
        the top-K cut.
    errors : List[str]
        One entry per path that failed (converter exception or non-empty
        ``error`` return).  Empty list when all paths succeeded.

    Parameters
    ----------
    paths
        Explicit path list.  When ``None`` and ``paths_mode='class_aware'``,
        the dispatcher classifies the SMILES and picks the relevant subset.
    paths_mode
        ``"class_aware"`` (default) — pick subset via class.  ``"all"`` —
        run every implemented path.  ``"explicit"`` — require ``paths`` to
        be non-None.
    include_placeholders
        Include ``requires_implementation=True`` paths.  No-op for now (their
        env-overrides are empty so they behave like head-default), useful
        once those are filled with real flag combinations.
    top_k_per_fingerprint
        Within each polyhedron-canonical-form fingerprint, keep at most this
        many best-scoring frames.  ``1`` reproduces strict deduplication.
    score_weights
        Per-metric weight overrides for ``score_frame``.  Missing keys fall
        back to ``DEFAULT_WEIGHTS``.
    log_dispatch
        Append a JSONL record to the adaptive-memory log.  Set ``False`` for
        unit-test isolation.

    Backward-compat
    ---------------
    Calling with v1 args ``(smiles, paths=DEFAULT_PATHS, max_isomers_per_path,
    max_isomers_total)`` still works.  The default ``top_k_per_fingerprint=2``
    is wider than v1's strict ``1`` — set ``top_k_per_fingerprint=1`` for
    bit-equivalent v1 behaviour (one frame per fingerprint).  All other
    arguments default to v1-compatible behaviour.
    """
    t0 = time.time()

    # Lazy imports to avoid pulling RDKit / scoring at module load.
    from delfin.smiles_converter import (
        smiles_to_xyz_isomers,
        _classify_complex_class,
    )
    from delfin.manta._ensemble_scoring import (
        score_frame,
        fingerprint as _polyhedron_fingerprint,
        build_smiles_pair_counter,
        DEFAULT_WEIGHTS,
    )
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        Chem = None  # type: ignore

    # Resolve weights (merge user → defaults so partial overrides work)
    weights = dict(DEFAULT_WEIGHTS)
    if score_weights:
        weights.update(score_weights)

    # Resolve mol + class once per call (mol can be None on parse failure)
    mol = None
    chem_class = "no_metal"
    if Chem is not None:
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                try:
                    mol.UpdatePropertyCache(strict=False)
                except Exception:
                    pass
                chem_class = _classify_complex_class(mol)
        except Exception:
            mol = None
    smi_pair_counter = build_smiles_pair_counter(mol)

    # Resolve path list per mode
    if paths is None:
        if paths_mode == "all":
            paths = [p for p in DEFAULT_PATHS
                     if include_placeholders or not p.requires_implementation]
        elif paths_mode == "explicit":
            raise ValueError(
                "paths_mode='explicit' requires `paths` to be non-None"
            )
        else:  # class_aware
            paths = select_paths_for_class(chem_class,
                                           include_placeholders=include_placeholders)

    # Per-path candidate accumulator
    per_path_frames: Dict[str, List[Tuple[str, str, float, Dict[str, float]]]] = {}
    errors: List[str] = []
    per_path_n_emitted: Dict[str, int] = {}

    for path in paths:
        per_path_frames[path.name] = []
        per_path_n_emitted[path.name] = 0
        if path.requires_implementation:
            # Placeholder — skip silently.  Not an error.
            continue

        # Snapshot env, apply overrides, restore on exit.
        snapshot = {k: os.environ.get(k) for k in path.env_overrides}
        try:
            for k, v in path.env_overrides.items():
                os.environ[k] = v
            try:
                results, err = smiles_to_xyz_isomers(
                    smiles,
                    max_isomers=max_isomers_per_path,
                    quality_mode=quality_mode,
                    **converter_kwargs,
                )
            except Exception as exc:
                errors.append(f"path={path.name} exception: {exc!r}")
                continue
            if err:
                errors.append(f"path={path.name} converter_err: {err}")
                # err does not always mean empty results — keep going if results exist
                if not results:
                    continue
            for xyz, label in (results or []):
                # Score the frame
                try:
                    score, _breakdown = score_frame(
                        xyz, smiles, mol=mol,
                        weights=weights,
                        smi_pair_counter=smi_pair_counter,
                    )
                except Exception as exc:
                    score = float("inf")
                    _breakdown = {}
                per_path_frames[path.name].append(
                    (xyz, label, score, _breakdown)
                )
                per_path_n_emitted[path.name] += 1
        finally:
            for k, prev in snapshot.items():
                if prev is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = prev

    # ----- Group by fingerprint, keep top-K per group ----- #
    # Build flat list of (xyz, label, score, source_path, breakdown)
    flat: List[Tuple[str, str, float, str, Dict[str, float]]] = []
    for src, frames in per_path_frames.items():
        for xyz, label, score, breakdown in frames:
            flat.append((xyz, label, score, src, breakdown))

    # Group by fingerprint
    grouped: Dict[str, List[Tuple[str, str, float, str, Dict[str, float]]]] = {}
    for entry in flat:
        xyz = entry[0]
        # Combine polyhedron fingerprint with isomer-label-root so different
        # isomers at the same polyhedron (e.g. cis vs trans on OH) don't
        # collapse — fingerprint is geometry+donor-elements only.
        try:
            fp_geom = _polyhedron_fingerprint(xyz)
        except Exception:
            fp_geom = "fp:err"
        fp_label = _isomer_label_root(entry[1]) if deduplicate else entry[1]
        fp = (fp_geom, fp_label)
        # Use string key for JSON serialisability + dict-hash speed
        fp_key = f"{fp_geom}::{fp_label}"
        grouped.setdefault(fp_key, []).append(entry)

    # Within each group, sort ascending by score, keep top-K.
    output: List[Tuple[str, str, float, str]] = []
    winning_paths: Dict[str, int] = {}
    n_unique_fingerprints = len(grouped)
    for fp_key, entries in grouped.items():
        entries.sort(key=lambda e: (e[2], e[3]))  # score asc, name tiebreak
        kept = entries[: max(1, int(top_k_per_fingerprint))]
        for xyz, label, score, src, _bd in kept:
            output.append((xyz, label, score, src))
            winning_paths[src] = winning_paths.get(src, 0) + 1
            if len(output) >= max_isomers_total:
                break
        if len(output) >= max_isomers_total:
            break

    # Sort final output ascending by score for caller convenience (best first).
    output.sort(key=lambda e: e[2])
    output = output[:max_isomers_total]

    # Recompute winning_paths after the global cap (output may have been
    # truncated post-grouping above).  Cheap O(n).
    winning_paths = {}
    for _, _, _, src in output:
        winning_paths[src] = winning_paths.get(src, 0) + 1

    # ----- Adaptive-memory JSONL log ----- #
    if log_dispatch:
        # Aggregate per-path score stats (min/median/max)
        scores_per_path: Dict[str, Dict[str, float]] = {}
        for src, frames in per_path_frames.items():
            if not frames:
                scores_per_path[src] = {"n": 0}
                continue
            sc = sorted(s for _x, _l, s, _b in frames if s != float("inf"))
            if not sc:
                scores_per_path[src] = {"n": len(frames), "min": None,
                                         "median": None, "max": None}
                continue
            mid = sc[len(sc) // 2]
            scores_per_path[src] = {
                "n": len(frames),
                "min": round(sc[0], 3),
                "median": round(mid, 3),
                "max": round(sc[-1], 3),
            }
        record = {
            "ts": round(time.time(), 3),
            "smiles": smiles,
            "class": chem_class,
            "paths_run": [p.name for p in paths],
            "frames_per_path": dict(per_path_n_emitted),
            "scores_per_path": scores_per_path,
            "winning_paths": dict(winning_paths),
            "n_unique_fingerprints": n_unique_fingerprints,
            "n_output": len(output),
            "wall_s": round(time.time() - t0, 3),
            "n_errors": len(errors),
        }
        _write_dispatch_log(record)

    return output, errors


# ---------------------------------------------------------------------- #
# Convenience summary (preserved from v1)
# ---------------------------------------------------------------------- #

def ensemble_summary(
    output: Sequence[Tuple[str, str, float, str]],
) -> Dict[str, int]:
    """Per-path frame count summary for logging / dashboard display."""
    counts: Dict[str, int] = {}
    for _, _, _, src in output:
        counts[src] = counts.get(src, 0) + 1
    return counts
