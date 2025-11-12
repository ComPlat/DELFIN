"""Auto-selection rules for OCCUPIER sequences."""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.common.logging import get_logger
from delfin.deep_auto_tree import DEEP_AUTO_SETTINGS

logger = get_logger(__name__)

STATE_FILENAME = ".delfin_occ_auto_state.json"


def _seq(entries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [dict(item) for item in entries]

# `- 0
#    |- baseline
#    |   |- even : [3]
#    |   `- odd  : [3]
#    `- branches
#        |- even
#        |   |- FoB 1
#        |   |   |- +1 : [3]
#        |   |   |- -1 : [4]
#        |   |   |- +2 : [3]
#        |   |   |- -2 : [5]
#        |   |   |- +3 : [3]
#        |   |   `- -3 : [6]
#        |   |- FoB 2
#        |   |   |- +1 : [3]
#        |   |   |- -1 : [3]
#        |   |   |- +2 : [3]
#        |   |   |- -2 : [6]
#        |   |   |- +3 : [3]
#        |   |   `- -3 : [6]
#        |   `- FoB 3
#        |       |- +1 : [3]
#        |       |- -1 : [3]
#        |       |- +2 : [3]
#        |       |- -2 : [6]
#        |       |- +3 : [3]
#        |       `- -3 : [6]
#        `- odd
#            |- FoB 1
#            |   |- +1 : [3]
#            |   |- -1 : [3]
#            |   |- +2 : [3]
#            |   |- -2 : [4]
#            |   |- +3 : [3]
#            |   `- -3 : [3]
#            |- FoB 2
#            |   |- +1 : [3]
#            |   |- -1 : [3]
#            |   |- +2 : [3]
#            |   |- -2 : [6]
#            |   |- +3 : [3]
#            |   `- -3 : [6]
#            `- FoB 3
#                |- +1 : [3]
#                |- -1 : [3]
#                |- +2 : [3]
#                |- -2 : [6]
#                |- +3 : [3]
#                `- -3 : [6]


AUTO_SETTINGS_FLAT: Dict[int, Dict[str, Any]] = {
    0: {
        "baseline": {
            "even": _seq([
                {"index": 1, "m": 1, "BS": "", "from": 0},
                {"index": 2, "m": 3, "BS": "", "from": 0},
                {"index": 3, "m": 5, "BS": "", "from": 0},
            ]),
            "odd": _seq([
                {"index": 1, "m": 2, "BS": "", "from": 0},
                {"index": 2, "m": 4, "BS": "", "from": 0},
                {"index": 3, "m": 6, "BS": "", "from": 0},
            ]),
        },
        "branches": {
            "even": {
                1: {
                    +1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 5, "BS": "", "from": 3},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ]),
                },
                2: {
                    +1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "3,1", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
                3: {
                    +1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "4,1", "from": 3},
                        {"index": 4, "m": 6, "BS": "", "from": 3},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 3, "BS": "3,1", "from": 3},
                        {"index": 5, "m": 3, "BS": "4,2", "from": 3},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,3", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
            },
            "odd": {
                1: {
                    +1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 3, "BS": "", "from": 0},
                        {"index": 4, "m": 5, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 4, "BS": "", "from": 0},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "1,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "3,1", "from": 1},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                },
                2: {
                    +1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "3,1", "from": 1},
                        {"index": 4, "m": 6, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "3,2", "from": 1},
                        {"index": 3, "m": 1, "BS": "2,1", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 5, "BS": "", "from": 0},
                        {"index": 6, "m": 5, "BS": "4,1", "from": 4},
                    ]),
                    +3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 1, "BS": "2,1", "from": 1},
                        {"index": 3, "m": 1, "BS": "3,2", "from": 1},
                        {"index": 4, "m": 3, "BS": "", "from": 0},
                        {"index": 5, "m": 3, "BS": "4,1", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
                    ]),
                },
                3: {
                    +1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -1: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 1},
                        {"index": 3, "m": 5, "BS": "5,1", "from": 1},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 4, "BS": "4,1", "from": 4},
                        {"index": 4, "m": 4, "BS": "5,2", "from": 4},
                        {"index": 5, "m": 6, "BS": "", "from": 0},

                    ]),
                    +3: _seq([
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 5, "BS": "", "from": 0},
                    ]),
                    -3: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "3,1", "from": 1},
                        {"index": 3, "m": 2, "BS": "3,2", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,2", "from": 4},
                        {"index": 5, "m": 4, "BS": "5,3", "from": 4},
                        {"index": 6, "m": 6, "BS": "", "from": 0},
                    ]),
                },
            },
        },
    },
}

AUTO_SETTINGS: Dict[int, Dict[str, Any]] = DEEP_AUTO_SETTINGS


def _resolve_root(base_dir: Optional[Path] = None) -> Path:
    if base_dir is not None:
        return Path(base_dir)
    env_root = os.environ.get("DELFIN_OCC_ROOT")
    if env_root:
        return Path(env_root)
    cwd = Path.cwd()
    if cwd.name.endswith("_OCCUPIER"):
        return cwd.parent
    return cwd


def _state_path(root: Path) -> Path:
    return root / STATE_FILENAME


def _load_state(root: Path) -> Dict[str, Any]:
    path = _state_path(root)
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
            if isinstance(data, dict):
                return data
    except Exception as exc:  # noqa: BLE001
        logger.debug("Failed to read auto-state %s: %s", path, exc)
    return {}


def _save_state(state: Dict[str, Any], root: Path) -> None:
    path = _state_path(root)
    try:
        with path.open("w", encoding="utf-8") as fh:
            json.dump(state, fh, indent=2, sort_keys=True)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to write auto-state %s: %s", path, exc)


def _parity_token(parity: str) -> str:
    return "even" if parity.lower().startswith("even") else "odd"


def infer_parity_from_m(value: Any, fallback: Optional[str] = None) -> Optional[str]:
    try:
        m_val = int(value)
    except (TypeError, ValueError):
        return fallback
    return "even" if m_val % 2 == 1 else "odd"


def record_auto_preference(parity: str, preferred_index: Optional[int], delta: int,
                           *, root: Optional[Path] = None) -> None:
    """Remember which FoB index won for a baseline (anchor) delta."""
    if preferred_index is None:
        return
    anchor = AUTO_SETTINGS.get(delta)
    if not anchor:
        return
    root_path = _resolve_root(root)
    state = _load_state(root_path)
    entry = state.setdefault(str(delta), {})
    entry[_parity_token(parity)] = int(preferred_index)
    _save_state(state, root_path)
    logger.debug("[occupier_auto] Recorded preferred index %s for parity=%s anchor=%s at %s",
                 preferred_index, parity, delta, root_path)


def _preferred_index_for_anchor(anchor: int, parity: str, root: Path) -> Optional[int]:
    token = _parity_token(parity)
    state = _load_state(root)
    anchor_state = state.get(str(anchor))
    if not isinstance(anchor_state, dict):
        return None
    value = anchor_state.get(token)
    if value is None:
        return None
    try:
        return int(value)
    except Exception:  # noqa: BLE001
        return None


def _extract_sequence_from_tree_node(tree_node: Dict[str, Any], preferred_sub_index: int = 1) -> Optional[List[Dict[str, Any]]]:
    """Extract sequence from tree node structure.

    Tree nodes have structure: {sub_index: {"seq": [...], "branches": {...}}}
    We prefer the given sub_index, or take the first available.
    """
    if not isinstance(tree_node, dict):
        return None

    # Try preferred sub-index first
    if preferred_sub_index in tree_node:
        node_data = tree_node[preferred_sub_index]
        if isinstance(node_data, dict) and "seq" in node_data:
            return node_data["seq"]

    # Fallback: take first available sub-index
    for sub_idx in sorted(tree_node.keys()):
        node_data = tree_node.get(sub_idx)
        if isinstance(node_data, dict) and "seq" in node_data:
            return node_data["seq"]

    return None


_TREE_DATASETS = {
    "deep": AUTO_SETTINGS,
    "tree": AUTO_SETTINGS,
    "flat": AUTO_SETTINGS_FLAT,
    "flatt": AUTO_SETTINGS_FLAT,
}


def resolve_auto_sequence_bundle(delta: int, *, root: Optional[Path] = None,
                                 tree_mode: str = "deep") -> Dict[str, List[Dict[str, Any]]]:
    """Return auto-managed sequences for the requested delta (if available)."""
    normalized_mode = str(tree_mode or "deep").strip().lower()
    settings_source = _TREE_DATASETS.get(normalized_mode, AUTO_SETTINGS)
    root_path = _resolve_root(root)
    bundle: Dict[str, List[Dict[str, Any]]] = {}

    def _storage_key(seq: List[Dict[str, Any]], fallback_parity: str) -> str:
        m_candidate = None
        for entry in seq:
            if isinstance(entry, dict) and entry.get("m") is not None:
                m_candidate = entry.get("m")
                break
        parity = infer_parity_from_m(m_candidate, fallback_parity)
        parity = parity if parity in ("even", "odd") else fallback_parity
        return f"{parity}_seq"

    def _copy_sequence(seq: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return copy.deepcopy(seq)

    for anchor, settings in settings_source.items():
        offset = delta - anchor
        baseline = settings.get("baseline", {})
        branches = settings.get("branches", {})

        if offset == 0:
            for parity in ("even", "odd"):
                seq = baseline.get(parity)
                if seq:
                    key = _storage_key(seq, parity)
                    bundle[key] = _copy_sequence(seq)
            if bundle:
                return bundle
            continue

        for parity in ("even", "odd"):
            parity_branches = branches.get(parity, {})
            if not parity_branches:
                continue

            pref_index = _preferred_index_for_anchor(anchor, parity, root_path)
            ordered_indices: list[int] = []
            if pref_index is not None and pref_index in parity_branches:
                ordered_indices.append(pref_index)
            # Append remaining candidates in ascending order for determinism
            for candidate in sorted(parity_branches.keys()):
                if candidate not in ordered_indices:
                    ordered_indices.append(candidate)

            for branch_index in ordered_indices:
                delta_node = parity_branches.get(branch_index, {}).get(offset)
                if not delta_node:
                    continue

                # Extract sequence from tree node structure
                if isinstance(delta_node, list):
                    seq = delta_node
                else:
                    seq = _extract_sequence_from_tree_node(delta_node, preferred_sub_index=1)
                if not seq:
                    continue

                key = _storage_key(seq, parity)
                bundle[key] = _copy_sequence(seq)
                if pref_index is None:
                    logger.debug(
                        "[occupier_auto] Using fallback branch %s for parity=%s, anchor=%s, offset=%s",
                        branch_index,
                        parity,
                        anchor,
                        offset,
                    )
                break

        if bundle:
            return bundle

    return {}
