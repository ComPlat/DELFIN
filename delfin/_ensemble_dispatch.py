"""Phase 1.5 ensemble dispatcher — combine K conversion paths per SMILES.

Each path is a configuration of DELFIN env-flags that activates a
specific FPCFD champion-port mechanism (Iter-8.x).  Running multiple
paths in parallel and pooling their candidate frames lets the ensemble
recover per-class wins without committing any single path as HEAD's
default.  No path's env-flag is required to be on at HEAD level: the
ensemble adds them at runtime.

v1 scope (this module):
    • Sequential in-process path execution (env-flags set + reset
      around each smiles_to_xyz_isomers call)
    • Candidate-pool aggregation — concatenate frames from all paths
    • Fingerprint deduplication on heavy-atom signature
    • Provenance tagging — every output frame carries the path name
      that produced it

v2 scope (future):
    • Per-frame detector scoring + class-aware composite score
    • Top-K selection per coordination-isomer fingerprint
    • MACE refinement integration (Phase 2)
    • CCDC RMSD scoring (Phase 3)

Why in-process and not subprocess?
    DELFIN reads env-flags via _delfin_env_int at conversion time, not
    at module-load time, so changing os.environ between path invocations
    is sufficient.  In-process is ~5× faster than subprocess (no Python
    cold-start) and shares the RDKit / numpy caches.  Multi-threading is
    NOT safe (smiles_to_xyz_isomers uses module globals, e.g.
    _ITER84_SIGMA_CAPS_OVERRIDE) — call paths sequentially.

API
---
    from delfin._ensemble_dispatch import ensemble_smiles_to_xyz_isomers

    out = ensemble_smiles_to_xyz_isomers(
        smiles="...",
        paths=DEFAULT_PATHS,           # or pick subset
        max_isomers_per_path=20,       # cap per-path frames
        max_isomers_total=50,          # cap on output (post-dedup)
    )
    # out = [(xyz, label, score, path_name), ...]
    # score is None in v1; v2 will populate composite detector score

Path catalogue
--------------
    DEFAULT_PATHS — current FPCFD-port set (head + 4 iter-8.x paths).
    Add new paths by appending Path(name, env_overrides) tuples.
"""

from __future__ import annotations

import hashlib
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class Path:
    """One conversion path = name + env-var overrides applied at runtime."""

    name: str
    env_overrides: Dict[str, str] = field(default_factory=dict)


# Default catalogue of paths.  Order matters: HEAD baseline first, then
# Iter-8.4 sigma-class ports, then Iter-8.5 e6761e4 hapto / multi-class
# port.  Phase 1.5 v1 runs them all sequentially; v2 will use a
# class-aware dispatcher to select a subset based on classify(mol).
DEFAULT_PATHS: List[Path] = [
    Path("head", {}),
    Path(
        "iter8.4_sigma_full",
        {
            "DELFIN_SIGMA_PORT_123A130_ITER8": "1",
            "DELFIN_SIGMA_SKIP_PUCKER_ITER8": "1",
            "DELFIN_SIGMA_TIGHT_THRESHOLD_ITER8": "1",
        },
    ),
    Path(
        "iter8.3_multisigma",
        {"DELFIN_MULTISIGMA_PORT_5B3E0D2_ITER8": "1"},
    ),
    Path(
        "iter8.5a_hapto_classes",
        {"DELFIN_ITER85_HARD_REJECT_CLASSES": "hapto,multi_hapto,multi_sigma"},
    ),
    Path(
        "iter8.1_multihapto_filter",
        {"DELFIN_ITER81_FILTER": "1"},  # default-on, listed for completeness
    ),
]


def _xyz_signature(xyz: str) -> str:
    """Heavy-atom-position fingerprint for deduplication.

    Strips comments and H atoms, normalises numeric formatting, and
    SHA-1 hashes the result.  Same heavy-atom geometry across paths
    yields the same signature regardless of path-specific labels.
    """
    if not xyz:
        return ""
    lines = xyz.strip().splitlines()
    if len(lines) < 3:
        return ""
    try:
        n = int(lines[0].strip())
    except Exception:
        return ""
    body: List[str] = []
    for ln in lines[2 : 2 + n]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        if sym == "H":
            continue
        # Round to 3 decimals to allow tiny numerical jitter
        try:
            x = round(float(parts[1]), 3)
            y = round(float(parts[2]), 3)
            z = round(float(parts[3]), 3)
        except Exception:
            continue
        body.append(f"{sym} {x:.3f} {y:.3f} {z:.3f}")
    body.sort()
    h = hashlib.sha1("\n".join(body).encode("utf-8")).hexdigest()
    return h[:16]


def _isomer_label_root(label: str) -> str:
    """Extract the isomer-fingerprint-equivalent label root.

    DELFIN labels look like "trans-OH N1+N1|Br0+Br0|N1+N1" or
    "all-trans-conf2"; the conf-suffix and ordering of the "|"-separated
    type tuples can vary across paths.  This function returns a
    canonical form for cross-path matching.
    """
    if not label:
        return ""
    label = label.strip()
    # strip trailing -confN
    label = re.sub(r"-conf\d+$", "", label)
    # Split prefix (everything up to first space) from type-list
    if " " in label:
        prefix, types = label.split(" ", 1)
        # sort the |-separated tokens
        toks = sorted(t.strip() for t in types.split("|") if t.strip())
        return prefix.strip() + " " + "|".join(toks)
    return label


def ensemble_smiles_to_xyz_isomers(
    smiles: str,
    paths: Optional[List[Path]] = None,
    max_isomers_per_path: int = 30,
    max_isomers_total: int = 50,
    quality_mode: Optional[str] = "fast",
    deduplicate: bool = True,
    **converter_kwargs,
) -> Tuple[List[Tuple[str, str, Optional[float], str]], List[str]]:
    """Run K conversion paths and pool their candidate frames.

    Returns (output, errors) where:
        output = [(xyz, label, score, source_path), ...]
            score is None in v1 (no detector scoring at runtime).
            source_path is the Path.name that produced this frame.
        errors = list of human-readable error strings, one per path
            that failed.  Empty list if all paths succeeded.

    Deduplication runs on (heavy-atom signature, isomer label root): a
    given coordination-isomer / geometry combination is emitted at most
    once even if multiple paths produced equivalent frames.  The first
    path to produce a given fingerprint wins (deterministic order from
    ``paths`` argument).
    """
    if paths is None:
        paths = DEFAULT_PATHS

    # Lazy import to avoid pulling RDKit at module load if the caller
    # doesn't actually use the converter.
    from delfin.smiles_converter import smiles_to_xyz_isomers

    pooled: List[Tuple[str, str, Optional[float], str]] = []
    errors: List[str] = []
    seen_fingerprints: set = set()

    for path in paths:
        # Snapshot env, apply overrides, restore on exit.  ``finally``
        # block guarantees clean state even if the converter raises.
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
                continue
            for xyz, label in results or []:
                if deduplicate:
                    sig = _xyz_signature(xyz)
                    root = _isomer_label_root(label)
                    fp = (sig, root)
                    if fp in seen_fingerprints:
                        continue
                    seen_fingerprints.add(fp)
                pooled.append((xyz, label, None, path.name))
                if len(pooled) >= max_isomers_total:
                    break
        finally:
            # Restore env
            for k, prev in snapshot.items():
                if prev is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = prev
        if len(pooled) >= max_isomers_total:
            break

    return pooled, errors


def ensemble_summary(
    output: List[Tuple[str, str, Optional[float], str]],
) -> Dict[str, int]:
    """Per-path frame count summary for logging / dashboard display."""
    counts: Dict[str, int] = {}
    for _, _, _, src in output:
        counts[src] = counts.get(src, 0) + 1
    return counts
