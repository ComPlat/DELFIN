"""Helpers for OCCUPIER sequence profiles."""
from __future__ import annotations

import copy
import os
import re
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.common.logging import get_logger

from .occupier_auto import resolve_auto_sequence_bundle

_DELTA_MARKER = ".delfin_occ_delta"
_FOLDER_PATTERN = re.compile(r"(ox|red)(?:_step)?_(\d+)", re.IGNORECASE)
logger = get_logger(__name__)


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
            parsed = _coerce_int(delta)
            if parsed not in cache:
                cache[parsed] = entry
    return cache


def _ensure_sequence_cache(config: Dict[str, Any]) -> Dict[int, Dict[str, List[Dict[str, Any]]]]:
    cache = config.get("_occupier_sequence_cache")
    if isinstance(cache, dict):
        return cache
    blocks = config.get("_occupier_sequence_blocks", [])
    cache = _build_sequence_cache(blocks)
    config["_occupier_sequence_cache"] = cache
    return cache


def _resolve_manual_sequences(config: Dict[str, Any], delta: int) -> Dict[str, List[Dict[str, Any]]]:
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


def resolve_sequences_for_delta(config: Dict[str, Any], delta: int) -> Dict[str, List[Dict[str, Any]]]:
    """Return copies of even/odd sequences for a specific charge delta."""
    method = str(config.get("OCCUPIER_method", "manually") or "manually").strip().lower()
    if method == "auto":
        auto_bundle = resolve_auto_sequence_bundle(delta)
        if auto_bundle:
            missing_even = "even_seq" not in auto_bundle
            missing_odd = "odd_seq" not in auto_bundle
            if missing_even or missing_odd:
                manual_bundle = _resolve_manual_sequences(config, delta)
                if missing_even and "even_seq" in manual_bundle:
                    auto_bundle["even_seq"] = manual_bundle["even_seq"]
                if missing_odd and "odd_seq" in manual_bundle:
                    auto_bundle["odd_seq"] = manual_bundle["odd_seq"]
            return auto_bundle
        logger.debug("[occupier_sequences] No auto sequence rule for delta=%s; falling back to CONTROL definitions.", delta)

    return _resolve_manual_sequences(config, delta)


def parse_species_delta(value: Any) -> int:
    """Convert value to int with default 0."""
    return _coerce_int(value, fallback=0)


def write_species_delta_marker(folder: Path, delta: int) -> None:
    """Persist the species delta inside a stage folder (for later reuse)."""
    try:
        marker = folder / _DELTA_MARKER
        marker.write_text(f"{delta}\n", encoding="utf-8")
    except Exception:
        # Non-fatal best-effort write
        pass


def read_species_delta_marker(folder: Path) -> Optional[int]:
    """Read the stored species delta if available."""
    marker = folder / _DELTA_MARKER
    if not marker.exists():
        return None
    try:
        value = marker.read_text(encoding="utf-8").strip()
        if not value:
            return None
        return _coerce_int(value)
    except Exception:
        return None


def _infer_delta_from_name(name: str) -> Optional[int]:
    lowered = name.lower()
    if lowered.startswith("initial"):
        return 0
    match = _FOLDER_PATTERN.search(lowered)
    if match:
        kind = match.group(1).lower()
        magnitude = _coerce_int(match.group(2), fallback=0)
        if magnitude == 0:
            return 0
        sign = 1 if kind.startswith("ox") else -1
        return sign * magnitude
    return None


def infer_species_delta(folder: Optional[Path] = None, default: int = 0) -> int:
    """Infer the species delta from env markers, files, or folder names."""
    env_value = os.environ.get("DELFIN_OCCUPIER_DELTA")
    if env_value not in (None, ""):
        try:
            return _coerce_int(env_value, fallback=default)
        except Exception:
            pass

    target_folder = folder or Path.cwd()
    marker_value = read_species_delta_marker(target_folder)
    if marker_value is not None:
        return marker_value

    guessed = _infer_delta_from_name(target_folder.name)
    if guessed is not None:
        return guessed

    return default


def append_sequence_overrides(control_path: Path, bundle: Dict[str, List[Dict[str, Any]]]) -> None:
    """Append auto-generated sequence blocks to CONTROL.txt for the stage."""
    if not bundle:
        return

    def _format_value(val: Any) -> str:
        if isinstance(val, str):
            return json.dumps(val, ensure_ascii=False)
        return str(val)

    def _format_entry(entry: Dict[str, Any]) -> str:
        ordered_keys = ["index", "m", "BS", "from"]
        remaining = [k for k in entry.keys() if k not in ordered_keys]
        parts: List[str] = []
        for key in ordered_keys + remaining:
            if key in entry:
                parts.append(f'"{key}": {_format_value(entry[key])}')
        return "{" + ", ".join(parts) + "}"

    lines: List[str] = [
        "\n# AUTO sequence overrides (generated by DELFIN)\n",
        "# The following blocks override previous even/odd_seq definitions for this stage.\n",
    ]

    sections = (
        ("even_seq", "even electron number (auto):"),
        ("odd_seq", "odd electron number (auto):"),
    )

    for key, heading in sections:
        seq = bundle.get(key)
        if not seq:
            continue
        lines.append(f"{heading}\n")
        lines.append(f"{key} = [\n")
        for entry in seq:
            lines.append(f"  {_format_entry(entry)},\n")
        lines.append("]\n")

    with control_path.open("a", encoding="utf-8") as fh:
        fh.writelines(lines)
