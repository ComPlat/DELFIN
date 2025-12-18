"""Helpers for OCCUPIER sequence profiles."""
from __future__ import annotations

import copy
import os
import re
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.common.logging import get_logger
from delfin.common.control_validator import _as_occupier_tree

from .occupier_auto import (
    resolve_auto_sequence_bundle,
    _parity_token,
)

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


def _extract_custom_tree_sequences(config: Dict[str, Any]) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    blocks = config.get("_occupier_sequence_blocks")
    target_block = None
    if isinstance(blocks, list) and blocks:
        for block in blocks:
            try:
                deltas = block.get("deltas") or []
            except AttributeError:
                continue
            if any(_coerce_int(delta) == 0 for delta in deltas):
                target_block = block
                break
        if target_block is None:
            target_block = blocks[0]
        if target_block:
            even_seq = copy.deepcopy(target_block.get("even_seq")) if isinstance(target_block.get("even_seq"), list) else []
            odd_seq = copy.deepcopy(target_block.get("odd_seq")) if isinstance(target_block.get("odd_seq"), list) else []
            return even_seq, odd_seq

    global_even = config.get("even_seq")
    global_odd = config.get("odd_seq")
    even_seq = copy.deepcopy(global_even) if isinstance(global_even, list) else []
    odd_seq = copy.deepcopy(global_odd) if isinstance(global_odd, list) else []
    return even_seq, odd_seq


def _get_custom_baseline_sequences(config: Dict[str, Any]) -> Optional[Dict[int, Dict[str, Any]]]:
    """Extract user-defined baseline sequences for own mode (no tree building).

    Returns a minimal dataset with only baseline sequences:
    {0: {"baseline": {"even": [...], "odd": [...]}}}
    """
    even_seq, odd_seq = _extract_custom_tree_sequences(config)
    if not even_seq and not odd_seq:
        return None

    # Return only baseline (no tree building needed for rule-based mode)
    return {
        0: {
            "baseline": {
                "even": even_seq,
                "odd": odd_seq,
            }
        }
    }


def _resolve_manual_sequences(config: Dict[str, Any], delta: int) -> Dict[str, List[Dict[str, Any]]]:
    method_token = str(config.get("OCCUPIER_method", "auto") or "auto").strip().lower()
    if method_token == "auto":
        return {}
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


def _infer_parity_hint(config: Dict[str, Any], delta: int) -> Optional[str]:
    neutral = config.get("_neutral_electrons")
    try:
        base_charge = _coerce_int(config.get("charge"), fallback=0)
    except Exception:
        base_charge = 0
    if neutral is not None:
        try:
            electrons = int(neutral) - (base_charge + delta)
            return "even" if electrons % 2 == 0 else "odd"
        except Exception:
            pass
    mult_guess = config.get("_multiplicity_guess")
    if mult_guess is not None:
        try:
            guess = int(mult_guess)
            base_even = (guess % 2 == 1)
            actual_even = base_even if (delta % 2 == 0) else not base_even
            return "even" if actual_even else "odd"
        except Exception:
            return None
    return None


def resolve_sequences_for_delta(config: Dict[str, Any], delta: int,
                                parity_hint: Optional[str] = None) -> Dict[str, List[Dict[str, Any]]]:
    """Return copies of even/odd sequences for a specific charge delta."""
    manual_bundle = _resolve_manual_sequences(config, delta)
    parity_target = parity_hint or _infer_parity_hint(config, delta)
    method = str(config.get("OCCUPIER_method", "auto") or "auto").strip().lower()
    if method == "auto":
        tree_mode = _as_occupier_tree(config.get("OCCUPIER_tree", "deep"))
        custom_tree = None
        if tree_mode == "own":
            custom_tree = _get_custom_baseline_sequences(config)
            if not custom_tree:
                logger.warning("OCCUPIER_tree=own requested but no custom sequences found; falling back to manual sequences.")
                return manual_bundle
        auto_bundle = resolve_auto_sequence_bundle(
            delta,
            tree_mode=tree_mode,
            parity_hint=parity_target or parity_hint,
            custom_dataset=custom_tree,
            config=config,
        )
        # Check if auto_bundle has actual sequences (not just empty dict)
        has_auto_sequences = auto_bundle and any(
            key in auto_bundle and auto_bundle[key]
            for key in ("even_seq", "odd_seq", "initial_sequence")
        )
        if has_auto_sequences:
            # Auto bundle found - augment with manual sequences if needed
            missing_even = "even_seq" not in auto_bundle
            missing_odd = "odd_seq" not in auto_bundle
            if missing_even and "even_seq" in manual_bundle:
                auto_bundle["even_seq"] = manual_bundle["even_seq"]
            if missing_odd and "odd_seq" in manual_bundle:
                auto_bundle["odd_seq"] = manual_bundle["odd_seq"]
            bundle = auto_bundle
        else:
            # No auto bundle or empty - fall back to manual sequences
            bundle = manual_bundle

        # For AUTO mode, always return both sequences (DELFIN chooses at runtime)
        return bundle
    else:
        bundle = manual_bundle

    # For manual mode, filter to single parity if hint provided
    if parity_target:
        key = f"{_parity_token(parity_target)}_seq"
        seq = bundle.get(key)
        if seq:
            return {key: seq}

    return bundle


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


def _strip_section(text: str, start_marker: str, end_markers: tuple[str, ...]) -> str:
    """Remove the text block starting at start_marker up to the first end marker."""
    start = text.find(start_marker)
    if start == -1:
        return text
    end = len(text)
    search_from = start + len(start_marker)
    for marker in end_markers:
        idx = text.find(marker, search_from)
        if idx != -1 and idx < end:
            end = idx
    return text[:start] + text[end:]


def _strip_infos_and_esd_sections(text: str) -> str:
    """Remove everything from INFOS: to the end of ESD MODULE section (including all dashes)."""
    # Find INFOS: marker
    infos_start = text.find("INFOS:")
    if infos_start == -1:
        return text

    # Look backwards to find the dash line before INFOS:
    line_start = text.rfind("\n", 0, infos_start)
    if line_start != -1:
        potential_dash_line = text[line_start+1:infos_start].strip()
        if potential_dash_line.startswith("-"):
            infos_start = line_start + 1

    # Find the last dash line after ESD MODULE (this closes the entire INFOS+ESD block)
    # We look for the third occurrence of dashes after INFOS:
    search_from = infos_start
    dash_count = 0
    end = -1

    while dash_count < 3:
        idx = text.find("-------------------------------------------------", search_from)
        if idx == -1:
            break
        dash_count += 1
        search_from = idx + 10
        if dash_count == 3:
            end = idx

    if end != -1:
        # Include the dash line and newline
        end_of_line = text.find("\n", end)
        if end_of_line != -1:
            end = end_of_line + 1
        else:
            end = end + len("-------------------------------------------------")
        return text[:infos_start] + text[end:]

    return text


def remove_existing_sequence_blocks(
    control_path: Path,
    force: bool = False,
    *,
    persist: bool = True,
    preserve_auto_sequences: bool = False,
) -> str:
    """Remove OCCUPIER sequence/INFO blocks before writing auto overrides.

    Args:
        control_path: Path to CONTROL.txt file
        force: If True, this is a copied CONTROL file - remove based on method.
               If False, this is main CONTROL - only remove OCCUPIER_sequence_profiles if method=auto.
    """
    try:
        original = control_path.read_text(encoding="utf-8")
    except OSError:
        return ""

    updated = original

    lower = original.lower()
    own_mode = "occupier_tree=own" in lower or "occupier_tree=dee4" in lower

    # In own-mode keep sequences in the main CONTROL (force=False) but strip
    # them from the OCCUPIER copies (force=True) so auto-tree handles them.
    if own_mode and not force:
        return original

    if preserve_auto_sequences:
        return original
    is_auto = "OCCUPIER_method=auto" in original

    if force:
        # For copied CONTROL files
        if is_auto:
            # Auto mode: remove all template sections (OCCUPIER_sequence_profiles, INFOS, ESD MODULE)
            updated = _strip_section(
                updated,
                "OCCUPIER_sequence_profiles:",
                ("INFOS:", "# AUTO sequence overrides"),
            )
            updated = _strip_infos_and_esd_sections(updated)
            updated = _strip_section(
                updated,
                "# AUTO sequence overrides",
                (),
            )
        else:
            # Manually mode: keep OCCUPIER_sequence_profiles, remove only INFOS and ESD MODULE
            updated = _strip_infos_and_esd_sections(updated)
    else:
        # For main CONTROL file: only remove OCCUPIER_sequence_profiles if method=auto
        if not is_auto:
            return original

        # Remove OCCUPIER_sequence_profiles section (but keep INFOS)
        updated = _strip_section(
            updated,
            "OCCUPIER_sequence_profiles:",
            ("INFOS:",),
        )
        # Remove old AUTO sequence overrides if present
        updated = _strip_section(
            updated,
            "# AUTO sequence overrides",
            ("ESD MODULE", "INFOS:"),
        )

    if updated != original and persist:
        control_path.write_text(updated.rstrip() + "\n", encoding="utf-8")

    return updated if updated != original else original
