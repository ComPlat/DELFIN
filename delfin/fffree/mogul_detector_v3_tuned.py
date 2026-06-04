"""Mogul detector v3 with per-class auto-tuned thresholds.

This module is the Mission-1 output: a thin filter over
``delfin.fffree.mogul_detector_v3.detect_anomalies_v3`` that applies a
per-class severity threshold loaded from the CSV produced by
``scripts/mogul_v3_threshold_optimization.py``.

The lookup key per record is

    ``f"{axis}::{'-'.join(sorted(atom_syms))}"``

i.e. axis + sorted element symbols of the participating atoms. This is
the same key produced by the threshold-tuning script's ``_class_key``.

For a record whose class is in the tuned table, the threshold is the
optimal value for that class. For records whose class is NOT in the
table (under-represented during tuning), the global CCDC threshold is
used (``mogul_detector_v3.CCDC_MAD_THRESHOLD`` = 2.5).

Default state
-------------
The module is **default-OFF**. Set ``DELFIN_MOGUL_V3_TUNED=1`` to enable
per-class filtering. With the flag unset and no explicit ``mode`` arg,
the underlying v3 detector behaves exactly as before (byte-identical).

CSV format expected (auto-loaded from
``paper_data/mogul_v3_threshold_optimization.csv`` or env override):

    axis,class,n_observations,...,optimal_threshold,f1,...

Only ``axis``, ``class`` and ``optimal_threshold`` columns are required.

Public API
----------
``detect_anomalies_v3_tuned(syms, P, *, table_path=None)``
    Same return-type as ``detect_anomalies_v3`` but with per-class
    thresholds applied.

``load_threshold_table(path=None)``
    Returns ``{class_key: float threshold}`` dict.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from . import mogul_detector_v3 as _v3

__all__ = [
    "ENV_FLAG",
    "DEFAULT_TABLE_PATH",
    "detect_anomalies_v3_tuned",
    "load_threshold_table",
    "class_key_for_record",
]

ENV_FLAG = "DELFIN_MOGUL_V3_TUNED"
DEFAULT_TABLE_PATH = os.environ.get(
    "DELFIN_MOGUL_V3_TUNED_TABLE",
    "/home/qmchem_max/ComPlat/DELFIN/paper_data/mogul_v3_threshold_optimization.csv",
)

# module-level cache
_TABLE_CACHE: Optional[Dict[str, float]] = None
_TABLE_PATH_CACHED: Optional[str] = None


def class_key_for_record(rec: dict) -> str:
    """Canonical class key for a v3 detector record.

    axis + '::' + sorted atom_syms joined by '-'. Mirrors the key used by
    the tuning script.
    """
    axis = rec.get("axis", "")
    if axis in ("bond", "pooled_bond"):
        ax2 = "bonds"
    elif axis == "angle":
        ax2 = "angles"
    elif axis in ("torsion", "improper"):
        ax2 = "torsions"
    else:
        ax2 = axis
    atom_syms = sorted(rec.get("atom_syms") or [])
    return f"{ax2}::" + "-".join(atom_syms)


def load_threshold_table(path: Optional[str] = None,
                         force_reload: bool = False) -> Dict[str, float]:
    """Load per-class optimal threshold map from CSV.

    Returns an empty dict if the file is missing — callers fall back
    to global CCDC threshold (2.5).
    """
    global _TABLE_CACHE, _TABLE_PATH_CACHED
    p = path or DEFAULT_TABLE_PATH
    if not force_reload and _TABLE_CACHE is not None and _TABLE_PATH_CACHED == p:
        return _TABLE_CACHE
    out: Dict[str, float] = {}
    pp = Path(p)
    if pp.exists():
        try:
            with pp.open() as fh:
                r = csv.DictReader(fh)
                for row in r:
                    axis = (row.get("axis") or "").strip()
                    cls = (row.get("class") or "").strip()
                    thr_raw = row.get("optimal_threshold") or ""
                    try:
                        thr = float(thr_raw)
                    except Exception:
                        continue
                    if not axis or not cls:
                        continue
                    # the CSV "class" field is already axis::elems but we
                    # accept both fully-qualified and bare-elems form
                    key = cls if "::" in cls else f"{axis}::{cls}"
                    out[key] = thr
        except Exception:
            out = {}
    _TABLE_CACHE = out
    _TABLE_PATH_CACHED = p
    return out


def detect_anomalies_v3_tuned(
    syms: Sequence[str],
    P,
    *,
    table_path: Optional[str] = None,
    index: Optional[dict] = None,
) -> List[Dict]:
    """v3 detector with per-class tuned thresholds.

    The underlying detector is always run in ``mode='ccdc'`` (lenient
    threshold 2.5), and the result is post-filtered by the per-class
    table. Records whose class is NOT in the table pass through the
    default 2.5 threshold (i.e. behave exactly like v3 in CCDC mode).
    """
    if os.environ.get(ENV_FLAG, "") != "1" and table_path is None:
        # default OFF — pass through (use whatever the underlying detector
        # selected via DELFIN_MOGUL_V3_DETECTOR or mode default)
        return _v3.detect_anomalies_v3(syms, P, index=index)

    # collect all CCDC-mode candidates (threshold 2.5)
    raw = _v3.detect_anomalies_v3(syms, P, mode="ccdc", index=index)
    table = load_threshold_table(table_path)
    if not table:
        return raw

    kept: List[Dict] = []
    for rec in raw:
        key = class_key_for_record(rec)
        thr = table.get(key, _v3.CCDC_MAD_THRESHOLD)
        sev = float(rec.get("sev_mad", 0.0) or 0.0)
        if sev >= thr:
            rec2 = dict(rec)
            rec2["threshold_class"] = thr
            rec2["tuned"] = True
            kept.append(rec2)
    return kept
