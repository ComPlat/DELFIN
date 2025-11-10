"""Helpers for OCCUPIER sequence profiles."""
from __future__ import annotations

import copy
from typing import Any, Dict, List, Optional


def _coerce_int(value: Any, fallback: int = 0) -> int:
    """Best-effort conversion to int with fallback."""
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    text = str(value).strip() if value is not None else ""
    if not text:
        return fallback
    try:
        return int(text)
    except ValueError:
        return fallback


def _build_sequence_cache(blocks: Optional[List[Dict[str, Any]]]) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    """Convert raw block definitions into a delta -> sequence map."""
    cache: Dict[int, Dict[str, List[Dict[str, Any]]]] = {}
    if not blocks:
        return cache

    for block in blocks:
        deltas = block.get("deltas") or []
        even_seq = block.get("even_seq")
        odd_seq = block.get("odd_seq")
        if not deltas:
            continue
        entry: Dict[str, List[Dict[str, Any]]] = {}
        if isinstance(even_seq, list):
            entry["even_seq"] = even_seq
        if isinstance(odd_seq, list):
            entry["odd_seq"] = odd_seq
        if not entry:
            continue
        for delta in deltas:
            cache[_coerce_int(delta)] = entry
    return cache


def _ensure_sequence_cache(config: Dict[str, Any]) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    cache = config.get("_occupier_sequence_cache")
    if isinstance(cache, dict):
        return cache
    blocks = config.get("_occupier_sequence_blocks", [])
    cache = _build_sequence_cache(blocks)
    config["_occupier_sequence_cache"] = cache
    return cache


def resolve_sequences_for_delta(config: Dict[str, Any], delta: int) -> Dict[str, List[Dict[str, Any]]]:
    """Return copies of even/odd sequences for a specific charge delta."""
    cache = _ensure_sequence_cache(config)
    bundle: Dict[str, List[Dict[str, Any]]] = {}
    entry = cache.get(delta)
    if entry:
        if "even_seq" in entry:
            bundle["even_seq"] = copy.deepcopy(entry["even_seq"])
        if "odd_seq" in entry:
            bundle["odd_seq"] = copy.deepcopy(entry["odd_seq"])
        if bundle:
            return bundle

    # Fall back to global sequences if no block matched
    global_even = config.get("even_seq")
    global_odd = config.get("odd_seq")
    if isinstance(global_even, list):
        bundle["even_seq"] = copy.deepcopy(global_even)
    if isinstance(global_odd, list):
        bundle["odd_seq"] = copy.deepcopy(global_odd)
    return bundle


def parse_species_delta(value: Any) -> int:
    """Convert OCCUPIER_species_delta field to int with default 0."""
    return _coerce_int(value, fallback=0)
