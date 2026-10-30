"""Shared loader for the ORCA key* word dictionaries.

Each share ships one JSON file `delfin/orca_knowledge/dictionary/<share>.json`
with records of the shape:

    {"keyword": str, "aliases": [str], "category": str, "purpose": str,
     "requires": str|None, "excludes": str|None, "pairs_with": str|None,
     "defaults": str|None, "typical_use": str|None, "literature": str|None,
     "manual": {"page": str, "section": str}, "parser": str,
     "opi": str|None, "status": "grounded"|"draft"}

This module loads every *.json in this directory at import time and exposes:

    explain(word)   -> the record (dict) for a keyword/alias, case-blind, or None
    records(category=None) -> list of all records (filtered by category if given)

It is intentionally dependency-free: it must run on a plain CI interpreter and
must not import OPI or delfin internals.
"""
import json
import os

_DATADIR = os.path.dirname(os.path.abspath(__file__))
_LOADED = False
# keyword (lowercased) -> full record; every canonical keyword and each of its
# aliases maps to the same record.
_INDEX = {}


def _load():
    global _INDEX, _LOADED
    if _LOADED:
        return
    for name in sorted(os.listdir(_DATADIR)):
        if not name.endswith(".json"):
            continue
        path = os.path.join(_DATADIR, name)
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, ValueError):
            # A malformed share must not break lookup for the others; skip it.
            continue
        if isinstance(data, dict):
            for rec in data.values():
                if not isinstance(rec, dict) or "keyword" not in rec:
                    continue
                rec.setdefault("aliases", [])
                kw = str(rec["keyword"])
                _INDEX[kw.lower()] = rec
                for alias in rec["aliases"] or []:
                    _INDEX[str(alias).lower()] = rec
    _LOADED = True


def explain(word):
    """Return the record dict for a keyword or alias (case-blind), else None."""
    if word is None:
        return None
    _load()
    key = str(word).strip().lower()
    if key in _INDEX:
        return _INDEX[key]
    # hyphen/space tolerant lookup as a last resort: "no split ri j" -> "nosplitrij"
    collapsed = key.replace("-", "").replace(" ", "")
    return _INDEX.get(collapsed)


def records(category=None):
    """Return every record, optionally filtered to one category (case-blind)."""
    _load()
    seen = {}
    for rec in _INDEX.values():
        seen[id(rec)] = rec
    out = list(seen.values())
    if category is not None:
        c = str(category).lower()
        out = [r for r in out if str(r.get("category", "")).lower() == c]
    out.sort(key=lambda r: str(r.get("keyword", "")))
    return out
