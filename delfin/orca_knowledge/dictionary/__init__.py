"""ORCA keyword dictionary: every keyword and block option DELFIN may write.

Shares (swarm agents) each contribute a ``<share>.json`` in this directory;
this loader merges them into one record set. A record answers "what does
this keyword do, when is it used, what does it require or exclude, and how
is it used in the manual" so DELFIN never writes a keyword it cannot explain.

Public API:
    ``explain(word)``   case-blind, alias-resolved; the record dict or None.
    ``records(category=None)``   every record, or only records of one category.

The merge is dynamic: the loader re-reads the directory each call, so a
share file added later is picked up without a code change. With no share
files present the dictionaries are simply empty (``explain`` -> None,
``records`` -> []); nothing breaks.
"""
from __future__ import annotations

import json
import os
import re
from functools import lru_cache

__all__ = ("explain", "records")

#: Every stand-alone JSON in this directory is a share file.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))


@lru_cache(maxsize=1)
def _share_files() -> tuple[str, ...]:
    """The names of the share JSON files, sorted for determinism."""
    here = {name for name in os.listdir(_THIS_DIR)
            if name.endswith(".json")
            and not name.startswith("__")}
    return tuple(sorted(here))


def _load() -> tuple[list[dict], dict[str, "list[dict]"], set[str]]:
    """(all_records, alias_index lower->records, categories).

    A record's ``keyword`` and each of its ``aliases`` are lower-cased and
    indexed. A word may legitimately resolve to more than one record (the
    same spelling as a standalone keyword and as an alias) so the index
    holds a list and ``explain`` returns the keyword match in preference
    to an alias match.
    """
    all_records: list[dict] = []
    aliases: dict[str, list[dict]] = {}
    categories: set[str] = set()
    for filename in _share_files():
        try:
            with open(os.path.join(_THIS_DIR, filename),
                      encoding="utf-8") as fh:
                records_for_share = json.load(fh)
        except (OSError, ValueError):
            # A half-written share file must not break the loader.
            continue
        if not isinstance(records_for_share, list):
            continue
        for record in records_for_share:
            if not isinstance(record, dict):
                continue
            all_records.append(record)
            cats = record.get("category")
            if isinstance(cats, str) and cats:
                categories.add(cats)
            keys = {record.get("keyword")}
            for alias in record.get("aliases", []):
                if isinstance(alias, str) and alias:
                    keys.add(alias)
            for key in keys:
                if isinstance(key, str) and key:
                    aliases.setdefault(key.lower(), []).append(record)
    return all_records, aliases, categories


def records(category=None) -> list[dict]:
    """Every record, optionally narrowed to one category."""
    all_records, _aliases, _categories = _load()
    if category is None:
        return list(all_records)
    return [r for r in all_records if r.get("category") == category]


def explain(word: str):
    """The record for ``word`` (case-blind, aliases resolved) or None."""
    if not isinstance(word, str):
        return None
    _all, aliases, _categories = _load()
    key = word.lower()
    candidates = aliases.get(key)
    if not candidates:
        # hyphen/underscore/case tolerant fallback: the manual writes the
        # same keyword with dashes (``NoSplit-RI-J``) that the input allows
        # without them (``NoSplitRIJ``).
        collapsed = re.sub(r"[-_ ]", "", key)
        for akey, recs in aliases.items():
            if re.sub(r"[-_ ]", "", akey) == collapsed:
                candidates = recs
                break
    if not candidates:
        return None
    # Prefer the record whose exact keyword (not just an alias) is ``word``.
    for record in candidates:
        if isinstance(record.get("keyword"), str) and record["keyword"].lower() == key:
            return record
    return candidates[0]
