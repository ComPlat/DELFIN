"""Scalar calibration values that MANTA does not ship.

A few default-off construction paths were calibrated on structure statistics that are
not redistributed with this repository.  Their numbers are therefore not in the code:
such a path reads them from a tab-separated file named by an environment variable, and
without that file it has no value and falls back to the behaviour documented at its
call site (a covalent-radii, tetrahedral or regular-polygon reference, or simply the
behaviour of its switch being off).  The shipped configuration sets none of these
variables, so it never reads such a file.  See delfin/manta/PROVENANCE.md.

File format: one ``key<TAB>value`` row per entry.  Further columns are ignored, blank
lines and lines starting with ``#`` are skipped, and so is any row whose value column is
not a number (a header row, for example).
"""
from __future__ import annotations

import os
from typing import Dict

_CACHE: Dict[str, Dict[str, float]] = {}


def scalar_table(env_var: str) -> Dict[str, float]:
    """``key -> value`` from the file named by ``env_var``; empty when unset or unreadable.

    Cached per path, so a build reads each file once.
    """
    path = os.environ.get(env_var, "")
    if not path:
        return {}
    tbl = _CACHE.get(path)
    if tbl is not None:
        return tbl
    tbl = {}
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                try:
                    tbl[fields[0].strip()] = float(fields[1])
                except ValueError:
                    continue
    except OSError:
        tbl = {}
    _CACHE[path] = tbl
    return tbl
