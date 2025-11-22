"""Auto-selection rules for OCCUPIER sequences."""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.common.logging import get_logger
from delfin.deep2_auto_tree import DEEP2_AUTO_SETTINGS
from delfin.deep3_auto_tree import DEEP3_AUTO_SETTINGS
from delfin.deep_auto_tree import DEEP_AUTO_SETTINGS

logger = get_logger(__name__)

STATE_FILENAME = ".delfin_occ_auto_state.json"
PURE_WINDOW = int(os.environ.get("OWN_TREE_PURE_WINDOW", "1"))


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
                        {"index": 4, "m": 2, "BS": "4,3", "from": 4},
                        {"index": 5, "m": 4, "BS": "", "from": 0},
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
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 3, "BS": "3,1", "from": 1},
                        {"index": 4, "m": 5, "BS": "", "from": 0},
                    ]),
                    +2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 4, "BS": "", "from": 0},
                        {"index": 3, "m": 6, "BS": "", "from": 0},
                    ]),
                    -2: _seq([
                        {"index": 1, "m": 2, "BS": "", "from": 0},
                        {"index": 2, "m": 2, "BS": "3,2", "from": 1},
                        {"index": 3, "m": 2, "BS": "2,1", "from": 1},
                        {"index": 4, "m": 4, "BS": "", "from": 0},
                        {"index": 5, "m": 4, "BS": "4,1", "from": 4},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
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
                        {"index": 1, "m": 1, "BS": "", "from": 0},
                        {"index": 2, "m": 3, "BS": "", "from": 0},
                        {"index": 3, "m": 3, "BS": "3,1", "from": 2},
                        {"index": 4, "m": 3, "BS": "4,2", "from": 2},
                        {"index": 5, "m": 3, "BS": "5,3", "from": 2},
                        {"index": 6, "m": 5, "BS": "", "from": 0},
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


def _parity_with_distance(parity: str, distance: int) -> str:
    """Return parity flipped `distance` times."""
    token = _parity_token(parity)
    if distance % 2 == 0:
        return token
    return "odd" if token == "even" else "even"


def infer_parity_from_m(value: Any, fallback: Optional[str] = None) -> Optional[str]:
    try:
        m_val = int(value)
    except (TypeError, ValueError):
        return fallback
    return "even" if m_val % 2 == 1 else "odd"


def record_auto_preference(parity: str, preferred_index: Optional[int], delta: int,
                           *, m_value: Optional[int] = None, bs_value: Optional[str] = None,
                           root: Optional[Path] = None) -> None:
    """Remember which FoB index won for a given delta, including m and BS values.

    This stores preferences for all delta values with full winner info (index, m, BS)
    to enable rule-based sequence generation for own mode.
    """
    if preferred_index is None:
        return
    root_path = _resolve_root(root)
    state = _load_state(root_path)
    entry = state.setdefault(str(delta), {})

    # Store full winner info for own mode rule-based generation
    winner_info = {"index": int(preferred_index)}
    if m_value is not None:
        winner_info["m"] = int(m_value)
    if bs_value is not None:
        winner_info["BS"] = str(bs_value).strip()

    entry[_parity_token(parity)] = winner_info
    _save_state(state, root_path)
    logger.debug("[occupier_auto] Recorded winner index=%s m=%s BS=%s for parity=%s delta=%s at %s",
                 preferred_index, m_value, bs_value, parity, delta, root_path)


def _preferred_index_from_state(state: Dict[str, Any], delta: int,
                                parity: Optional[str] = None) -> Optional[int]:
    """Get preferred index from state (backwards compatible with old int format)."""
    entry = state.get(str(delta))
    if not isinstance(entry, dict):
        return None
    tokens: List[str]
    if parity is None:
        tokens = ["even", "odd"]
    else:
        token = _parity_token(parity)
        other = "odd" if token == "even" else "even"
        tokens = [token, other]
    for token in tokens:
        value = entry.get(token)
        if value is None:
            continue
        try:
            # New format: dict with {index, m, BS}
            if isinstance(value, dict):
                return int(value.get("index", 0))
            # Old format: just int
            return int(value)
        except Exception:  # noqa: BLE001
            return None
    return None


def _get_winner_info_from_state(state: Dict[str, Any], delta: int,
                                 parity: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Get full winner info (index, m, BS) from state.

    Args:
        state: The auto state dictionary
        delta: The delta value to look up
        parity: The parity to look up. If None, returns first found (even, then odd).
                If specified, ONLY returns that parity (no fallback to other).
    """
    entry = state.get(str(delta))
    if not isinstance(entry, dict):
        return None

    if parity is None:
        # No parity specified: try both (even first, then odd)
        tokens = ["even", "odd"]
    else:
        # Parity specified: ONLY return that parity (no fallback)
        tokens = [_parity_token(parity)]

    for token in tokens:
        value = entry.get(token)
        if value is None:
            continue
        if isinstance(value, dict):
            return value
        # Old format fallback
        try:
            return {"index": int(value), "m": None, "BS": None}
        except Exception:  # noqa: BLE001
            return None
    return None


def _collect_preference_chain(anchor: int, target: int, parity: str,
                              state: Dict[str, Any]) -> List[int]:
    """Collect recorded FoB preferences along the path from anchor to target."""
    if anchor == target:
        return []
    direction = 1 if target > anchor else -1
    steps = abs(target - anchor)
    chain: List[int] = []
    for step in range(1, steps + 1):
        delta_value = anchor + direction * step
        remaining = steps - step
        parity_token = _parity_with_distance(parity, remaining)
        pref = _preferred_index_from_state(state, delta_value, parity_token)
        if pref is not None:
            chain.append(pref)
    return chain


def _ordered_candidates(source: Dict[int, Any], preferred: Optional[int]) -> List[int]:
    """Return sorted candidate keys, honoring a preferred index (int/str tolerant)."""
    if not source:
        return []

    def _normalize(key: Any) -> Any:
        try:
            return int(key)
        except Exception:  # noqa: BLE001
            return key

    ordered: List[Any] = []
    if preferred is not None:
        for token in (preferred, str(preferred)):
            if token in source and token not in ordered:
                ordered.append(token)

    numeric: List[tuple[int, Any]] = []
    non_numeric: List[tuple[str, Any]] = []
    for key in source.keys():
        norm = _normalize(key)
        if isinstance(norm, int):
            numeric.append((norm, key))
        else:
            non_numeric.append((str(norm), key))

    for _, key in sorted(numeric, key=lambda pair: pair[0]):
        if key not in ordered:
            ordered.append(key)
    for _, key in sorted(non_numeric, key=lambda pair: pair[0]):
        if key not in ordered:
            ordered.append(key)

    return ordered


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


def _navigate_recursive_tree(node: Dict[str, Any], remaining_steps: int, direction: int,
                             preferred_chain: Optional[List[int]] = None) -> Optional[List[Dict[str, Any]]]:
    """Navigate recursive tree structure to find sequence at target depth.

    For deep, structure is: node['branches'][direction][sub]['branches'][direction][sub]...
    remaining_steps=0 means we're at the target node
    """
    if remaining_steps == 0:
        # We're at the target node, extract its sequence
        if "seq" in node:
            return node.get("seq")
        return None

    # Navigate one level deeper
    inner_branches = node.get("branches", {})
    if not inner_branches:
        return None

    # Branch keys may be stored as ints (in-memory) or strings (JSON)
    next_level = inner_branches.get(direction)
    if next_level is None:
        next_level = inner_branches.get(str(direction))
    if next_level is None:
        return None

    # Get the next level (dict of sub-indices)
    preferred_sub = None
    tail: List[int] = []
    if preferred_chain:
        preferred_sub = preferred_chain[0]
        tail = preferred_chain[1:]

    for sub_idx in _ordered_candidates(next_level, preferred_sub):
        result = _navigate_recursive_tree(next_level[sub_idx], remaining_steps - 1, direction, tail)
        if result:
            return result
    return None


def _sanitize_sequence_entries(seq: Optional[List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """Return a cleaned, index-ordered copy of the provided sequence."""
    if not isinstance(seq, list):
        return []
    sanitized: List[Dict[str, Any]] = []
    for raw in seq:
        if not isinstance(raw, dict):
            continue
        try:
            m_val = int(raw.get("m"))
        except (TypeError, ValueError):
            continue
        try:
            idx_val = int(raw.get("index", len(sanitized) + 1))
        except (TypeError, ValueError):
            idx_val = len(sanitized) + 1
        bs_val = str(raw.get("BS", "") or "").strip()
        try:
            from_val = int(raw.get("from", 0))
        except (TypeError, ValueError):
            from_val = 0
        sanitized.append({"index": idx_val, "m": m_val, "BS": bs_val, "from": from_val})

    sanitized.sort(key=lambda entry: entry["index"])
    for idx, entry in enumerate(sanitized, start=1):
        entry["index"] = idx
    return sanitized


def _extract_pure_values(seq: List[Dict[str, Any]]) -> List[int]:
    """Return ordered list of pure multiplicities from a sequence."""
    values: List[int] = []
    seen: set[int] = set()
    for entry in seq:
        if entry.get("BS"):
            continue
        m_val = entry.get("m")
        if isinstance(m_val, int) and m_val not in seen:
            values.append(m_val)
            seen.add(m_val)
    return values


def _extract_bs_pairs(seq: List[Dict[str, Any]]) -> List[Tuple[int, int]]:
    """Parse BS indices as (M, N) tuples."""
    pairs: List[Tuple[int, int]] = []
    seen: set[Tuple[int, int]] = set()
    for entry in seq:
        raw = str(entry.get("BS", "") or "").strip()
        if not raw or "," not in raw:
            continue
        try:
            m_val, n_val = [int(token.strip()) for token in raw.split(",", 1)]
        except Exception:  # noqa: BLE001
            continue
        pair = (m_val, n_val)
        if pair not in seen:
            pairs.append(pair)
            seen.add(pair)
    return pairs


class _CustomTreeBuilder:
    """Build adaptive tree datasets from user-supplied baseline sequences."""

    def __init__(self, even_seq: Optional[List[Dict[str, Any]]], odd_seq: Optional[List[Dict[str, Any]]],
                 pure_window: Optional[int] = None, progressive_from: bool = False):
        sanitized = {
            "even": _sanitize_sequence_entries(even_seq),
            "odd": _sanitize_sequence_entries(odd_seq),
        }
        if pure_window is None:
            inferred = PURE_WINDOW
        else:
            inferred = pure_window
        try:
            self.pure_window = max(0, int(inferred))
        except Exception:
            self.pure_window = max(0, PURE_WINDOW)
        self.progressive_from = bool(progressive_from)
        self.pure_map = {
            "even": _extract_pure_values(sanitized["even"]),
            "odd": _extract_pure_values(sanitized["odd"]),
        }

        # If a parity is empty, auto-generate standard m-values
        # even → odd m (1,3,5), odd → even m (2,4,6)
        if not self.pure_map["even"] and self.pure_map["odd"]:
            self.pure_map["even"] = self._windowed_baseline_from_other(self.pure_map["odd"], "even")
        elif not self.pure_map["odd"] and self.pure_map["even"]:
            self.pure_map["odd"] = self._windowed_baseline_from_other(self.pure_map["even"], "odd")

        self.bs_map = {
            "even": _extract_bs_pairs(sanitized["even"]),
            "odd": _extract_bs_pairs(sanitized["odd"]),
        }
        self.min_m = {parity: min(values) if values else 1 for parity, values in self.pure_map.items()}
        # If parity is empty, use 7 as max_m (not 1) to allow BS evolution in nested levels
        self.max_m = {parity: max(values) if values else 8 for parity, values in self.pure_map.items()}
        # Baseline uses only user-defined sequences (no auto-BS)
        # BS are added dynamically in branches based on winners
        self.baseline = {
            "even": self._generate_baseline_seq("even", include_initial_bs=True),
            "odd": self._generate_baseline_seq("odd", include_initial_bs=True),
        }

    def _generate_baseline_bs_candidates(self, parity: str) -> List[Tuple[int, int]]:
        """Generate BS(m-1,1) for all pure m values to enable reduction path BS initiation."""
        pure_values = self.pure_map.get(parity, [])
        bs_candidates: List[Tuple[int, int]] = []
        for m_val in pure_values:
            target = m_val - 1
            if target >= 1:
                bs_candidates.append((m_val, 1))
        return bs_candidates

    @staticmethod
    def _parity_matches(parity: str, m_val: int) -> bool:
        if m_val <= 0:
            return False
        if parity == "even":
            return m_val % 2 == 1
        return m_val % 2 == 0

    @staticmethod
    def _default_pure_values(parity: str) -> List[int]:
        return [1, 3, 5] if parity == "even" else [2, 4, 6]

    def _align_m_for_bs(self, M: int, N: int, parity: str) -> int:
        """Derive an aligned m value for BS(M,N) using m = (M - N) + 1."""
        raw = max(1, M - N + 1)
        limit_min = max(1, N)
        limit_max = 8
        if raw < limit_min:
            raw = limit_min
        if raw > limit_max:
            raw = limit_max

        def _adjust(value: int, target_parity: str) -> int:
            if self._parity_matches(target_parity, value):
                return value
            if value + 1 <= limit_max:
                candidate = value + 1
                if self._parity_matches(target_parity, candidate):
                    return candidate
            if value - 1 >= limit_min:
                candidate = value - 1
                if self._parity_matches(target_parity, candidate):
                    return candidate
            return value

        aligned = _adjust(raw, parity)
        if aligned < limit_min:
            aligned = limit_min
        if aligned > limit_max:
            aligned = limit_max
        return aligned

    def _generate_baseline_seq(self, parity: str, *, include_initial_bs: bool,
                               add_bs: Optional[List[Tuple[int, int]]] = None) -> List[Dict[str, Any]]:
        seq: List[Dict[str, Any]] = []
        pure_values = self.pure_map.get(parity, [])

        for m_val in pure_values:
            idx = len(seq) + 1
            seq.append({"index": idx, "m": m_val, "BS": "", "from": (idx - 1 if self.progressive_from else 0)})

        bs_pool: List[Tuple[int, int]] = []
        if include_initial_bs:
            bs_pool.extend(self.bs_map.get(parity, []))
        if add_bs:
            for pair in add_bs:
                if pair not in bs_pool:
                    bs_pool.append(pair)

        for M, N in bs_pool:
            m_bs = self._align_m_for_bs(M, N, parity)
            insert_idx = 0
            for i, entry in enumerate(seq):
                if entry["m"] <= m_bs:
                    insert_idx = i + 1
            seq.insert(insert_idx, {"index": insert_idx + 1, "m": m_bs, "BS": f"{M},{N}", "from": 0})

        pure_index_map: Dict[int, int] = {}
        for idx, entry in enumerate(seq, start=1):
            entry["index"] = idx
            if entry["BS"]:
                continue
            entry["from"] = idx - 1 if self.progressive_from else 0
            pure_index_map[entry["m"]] = idx

        for entry in seq:
            if entry["BS"]:
                entry["from"] = pure_index_map.get(entry["m"], 0)

        return seq

    def _windowed_baseline_from_other(self, source_values: List[int], target_parity: str) -> List[int]:
        """Derive baseline ladder for missing parity using window around existing values."""
        if not source_values:
            return self._default_pure_values(target_parity)

        window = getattr(self, "pure_window", 0)
        if window <= 0:
            return self._default_pure_values(target_parity)

        candidates: List[int] = []
        limit_min, limit_max = 1, 8
        for base in source_values:
            for delta in range(-window, window + 1):
                candidate = base + delta
                if candidate < limit_min or candidate > limit_max:
                    continue
                if not self._parity_matches(target_parity, candidate):
                    continue
                if candidate not in candidates:
                    candidates.append(candidate)
        if not candidates:
            return self._default_pure_values(target_parity)
        candidates.sort()
        return candidates

    def _pure_window_sequence(self, parity: str, center_m: int) -> List[Dict[str, Any]]:
        """Return pure-only sequence centered around the winning m."""
        window = getattr(self, "pure_window", 1)
        baseline = self.pure_map.get(parity) or self._default_pure_values(parity)
        if window <= 0:
            values = list(baseline)
        else:
            limit_min = min(baseline) if baseline else 1
            limit_max = max(baseline) if baseline else 8
            limit_min = min(limit_min, 1)
            limit_max = max(limit_max, 8)
            values: List[int] = []
            for delta in range(-window, window + 1):
                candidate = center_m + delta
                if candidate < limit_min or candidate > limit_max:
                    continue
                if not self._parity_matches(parity, candidate):
                    continue
                if candidate not in values:
                    values.append(candidate)
            if not values:
                values = list(baseline)
        values.sort()
        seq: List[Dict[str, Any]] = []
        for idx, m_val in enumerate(values, start=1):
            seq.append({"index": idx, "m": m_val, "BS": "", "from": (idx - 1 if self.progressive_from else 0)})
        return seq

    def _inject_bs_entries(self, seq: List[Dict[str, Any]], parity: str,
                           add_bs: Optional[List[Tuple[int, int]]]) -> List[Dict[str, Any]]:
        if not add_bs:
            return seq

        for M, N in add_bs:
            aligned = self._align_m_for_bs(M, N, parity)
            pure_idx = next((i for i, entry in enumerate(seq)
                             if entry["m"] == aligned and not entry.get("BS")), None)
            if pure_idx is None:
                continue  # window excludes this multiplicity
            insert_idx = pure_idx + 1
            while insert_idx < len(seq) and seq[insert_idx]["m"] == aligned and seq[insert_idx].get("BS"):
                insert_idx += 1
            bs_entry = {"index": 0, "m": aligned, "BS": f"{M},{N}", "from": 0}
            seq.insert(insert_idx, bs_entry)
            # reindex will happen later

        for idx, entry in enumerate(seq, start=1):
            entry["index"] = idx
            if entry.get("BS"):
                pure_idx = next((item["index"] for item in seq
                                 if item["m"] == entry["m"] and not item.get("BS")), 0)
                entry["from"] = pure_idx
        return seq

    def _evolve_bs(self, parity: str, M: int, N: int) -> List[Tuple[int, int]]:
        options: List[Tuple[int, int]] = []
        max_m = self.max_m.get(parity, max(M, N))
        min_m = self.min_m.get(parity, 1)

        if M + 1 <= max_m:
            options.append((M + 1, N))
        if N + 1 <= max_m and N + 1 <= M:
            options.append((M, N + 1))
        if M - 1 >= N and M - 1 >= min_m:
            options.append((M - 1, N))
        if N - 1 >= 1:
            options.append((M, N - 1))
        return options

    def _next_bs_candidates(self, parity: str, prev_m: int, prev_bs: str) -> List[Tuple[int, int]]:
        """Generate BS candidates based on previous winner.

        Rules:
        - If BS won: test BS(M±1,N) and BS(M,N±1)
        - If pure m won: test BS(m-1,1) to initiate BS
        """
        if prev_bs:
            # BS evolution: expand/reduce
            try:
                m_val, n_val = [int(token) for token in prev_bs.split(",", 1)]
            except Exception:
                return []
            return self._evolve_bs(parity, m_val, n_val)

        # Pure state won: initiate BS with BS(m-1,1)
        target = prev_m - 1
        if target >= 1:
            return [(target, 1)]
        return []

    def _generate_reduction_sequence(self, parity: str, prev_m: int, prev_bs: str) -> List[Dict[str, Any]]:
        """Generate reduction sequence - can include BS evolution."""
        seq = self._pure_window_sequence(parity, prev_m)
        # Always generate BS candidates (both for pure and BS winners)
        add_bs = self._next_bs_candidates(parity, prev_m, prev_bs)
        if add_bs:
            seq = self._inject_bs_entries(seq, parity, add_bs)
        return seq

    # _generate_oxidation_sequence, _build_recursive_branches, and build_tree removed
    # These were only used for tree building, which is no longer needed for own mode


def generate_sequence_from_winner_rules(
    even_seq: Optional[List[Dict[str, Any]]],
    odd_seq: Optional[List[Dict[str, Any]]],
    prev_m: int,
    prev_bs: str,
    target_parity: str,
    *,
    pure_window: Optional[int] = None,
    progressive_from: bool = False,
    include_bs: bool = True,
) -> List[Dict[str, Any]]:
    """Generate next sequence based on winner using rules (no tree needed).

    Rules:
    - Pure states: windowed around prev_m
    - BS candidates (if include_bs=True):
      - If BS won: evolve BS (M±1,N), (M,N±1)
      - If pure won: initiate BS with BS(m-1,1)
    """
    builder = _CustomTreeBuilder(even_seq, odd_seq, pure_window=pure_window,
                                 progressive_from=progressive_from)

    # Generate pure states
    seq = builder._pure_window_sequence(target_parity, prev_m)

    # Add BS candidates if requested
    if include_bs:
        add_bs = builder._next_bs_candidates(target_parity, prev_m, prev_bs)
        if add_bs:
            seq = builder._inject_bs_entries(seq, target_parity, add_bs)

    # Update progressive_from values after BS injection (if needed)
    if progressive_from:
        # Build pure_index_map (m -> index in OCCUPIER.txt)
        pure_index_map: dict[int, int] = {}
        for entry in seq:
            if not entry.get("BS"):
                pure_index_map[entry["m"]] = entry["index"]

        # Get sorted pure m values
        pure_m_values = sorted([m for m in pure_index_map.keys()])

        # Update from values for pure states (progressive)
        for entry in seq:
            if not entry.get("BS"):
                # Find position of this m in sorted pure m values
                m = entry["m"]
                pos = pure_m_values.index(m)
                if pos == 0:
                    entry["from"] = 0  # First pure state
                else:
                    # Point to previous pure state's index
                    prev_m = pure_m_values[pos - 1]
                    entry["from"] = pure_index_map[prev_m]
            else:
                # BS entries point to their pure m state
                entry["from"] = pure_index_map.get(entry["m"], 0)

    return seq


# build_custom_auto_tree and persist_custom_tree removed - no longer needed for own mode
# Own mode now uses rule-based generation instead of pre-built trees


_TREE_DATASETS = {
    "flat": AUTO_SETTINGS_FLAT,
    "deep2": DEEP2_AUTO_SETTINGS,
    "deep3": DEEP3_AUTO_SETTINGS,
    "deep": DEEP_AUTO_SETTINGS,
}


def _resolve_own_mode_sequences(
    delta: int,
    custom_dataset: Dict[int, Dict[str, Any]],
    *,
    root: Optional[Path] = None,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """Rule-based sequence generation for own mode (no tree navigation).

    For delta=0: Return baseline sequences
    For delta!=0: Generate sequences based on previous winner using rules
    """
    root_path = _resolve_root(root)
    state_cache = _load_state(root_path)
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

    # Get baseline configuration
    baseline_config = custom_dataset.get(0, {}).get("baseline", {})
    even_baseline = baseline_config.get("even", [])
    odd_baseline = baseline_config.get("odd", [])

    if not even_baseline and not odd_baseline:
        return {}

    # Extract config parameters
    pure_window = None
    progressive_from = False
    bs_in_oxidation = False
    if config:
        raw_window = config.get("OWN_TREE_PURE_WINDOW")
        if raw_window not in (None, ""):
            try:
                pure_window = int(str(raw_window).strip())
            except Exception:
                pass
        raw_progressive = str(config.get("OWN_progressive_from", "no")).strip().lower()
        progressive_from = raw_progressive in ("yes", "true", "1", "on")
        raw_bs_ox = str(config.get("OWN_BS_IN_OXIDATION", "no")).strip().lower()
        bs_in_oxidation = raw_bs_ox in ("yes", "true", "1", "on")

    # Delta 0: Return baseline
    if delta == 0:
        for parity in ("even", "odd"):
            seq = baseline_config.get(parity)
            if seq:
                key = _storage_key(seq, parity)
                bundle[key] = copy.deepcopy(seq)
        return bundle

    # Delta != 0: Generate from previous winner
    # Only generate the opposite parity of the previous winner
    direction = +1 if delta > 0 else -1
    prev_delta = delta - direction

    # Determine if BS should be included
    # Oxidation (delta > 0): use bs_in_oxidation setting (default no)
    # Reduction (delta < 0): always include BS (default yes)
    is_oxidation = delta > 0
    include_bs = bs_in_oxidation if is_oxidation else True

    # First, try to find which parity won at previous delta
    winner_info_even = _get_winner_info_from_state(state_cache, prev_delta, "even")
    winner_info_odd = _get_winner_info_from_state(state_cache, prev_delta, "odd")

    # Determine which parity won (prefer one with actual m value)
    if winner_info_even and winner_info_even.get("m") is not None:
        source_parity = "even"
        winner_info = winner_info_even
    elif winner_info_odd and winner_info_odd.get("m") is not None:
        source_parity = "odd"
        winner_info = winner_info_odd
    else:
        # No winner found - generate both parities as fallback
        logger.debug("[own_mode] No winner found for delta=%s, generating both parities", prev_delta)
        for target_parity in ("even", "odd"):
            source_parity_fb = "odd" if target_parity == "even" else "even"
            source_baseline = baseline_config.get(source_parity_fb, [])
            if source_baseline:
                mid_idx = len(source_baseline) // 2
                prev_m = source_baseline[mid_idx].get("m", 3)
                prev_bs = ""
            else:
                prev_m = 3 if source_parity_fb == "even" else 4
                prev_bs = ""

            seq = generate_sequence_from_winner_rules(
                even_baseline,
                odd_baseline,
                prev_m,
                prev_bs,
                target_parity,
                pure_window=pure_window,
                progressive_from=progressive_from,
                include_bs=include_bs,
            )
            if seq:
                key = _storage_key(seq, target_parity)
                bundle[key] = seq
        return bundle

    # Generate only the opposite parity of the winner
    target_parity = "odd" if source_parity == "even" else "even"
    prev_m = winner_info["m"]
    prev_bs = winner_info.get("BS", "")

    logger.debug(
        "[own_mode] Winner from delta=%s was parity=%s (m=%s, BS=%s), generating parity=%s",
        prev_delta, source_parity, prev_m, prev_bs, target_parity
    )

    # Generate sequence using rules
    seq = generate_sequence_from_winner_rules(
        even_baseline,
        odd_baseline,
        prev_m,
        prev_bs,
        target_parity,
        pure_window=pure_window,
        progressive_from=progressive_from,
        include_bs=include_bs,
    )

    if seq:
        key = _storage_key(seq, target_parity)
        logger.debug(
            "[own_mode] Generated sequence for delta=%s target_parity=%s: %d entries, m-values=%s, storage_key=%s",
            delta, target_parity, len(seq), [e.get("m") for e in seq], key
        )
        bundle[key] = seq

    return bundle


def resolve_auto_sequence_bundle(delta: int, *, root: Optional[Path] = None,
                                 tree_mode: str = "deep",
                                 parity_hint: Optional[str] = None,
                                 custom_dataset: Optional[Dict[int, Dict[str, Any]]] = None,
                                 config: Optional[Dict[str, Any]] = None) -> Dict[str, List[Dict[str, Any]]]:
    """Return auto-managed sequences for the requested delta (if available)."""
    normalized_mode = str(tree_mode or "deep").strip().lower()

    # Special handling for own mode: rule-based generation (no tree navigation)
    if normalized_mode == "own":
        if not custom_dataset:
            return {}
        return _resolve_own_mode_sequences(
            delta, custom_dataset, root=root, config=config
        )

    settings_source = _TREE_DATASETS.get(normalized_mode, AUTO_SETTINGS)
    root_path = _resolve_root(root)
    state_cache = _load_state(root_path)
    bundle: Dict[str, List[Dict[str, Any]]] = {}
    requested_parities: tuple[str, ...]
    # For AUTO trees, always provide both parities (ORCA chooses at runtime)
    # parity_hint is only used for manual/legacy modes
    requested_parities = ("even", "odd")

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

    # Check if this is a recursive tree (deep has recursive structure)
    is_recursive_tree = normalized_mode == "deep"

    oxidation_baseline = None
    if normalized_mode == "deep":
        base_settings = settings_source.get(0)
        if base_settings:
            oxidation_baseline = base_settings.get("baseline", {})

    for anchor_key, settings in settings_source.items():
        try:
            anchor = int(anchor_key)
        except Exception:  # noqa: BLE001
            logger.debug("[occupier_auto] Skipping non-numeric anchor key: %s", anchor_key)
            continue
        offset = delta - anchor
        baseline = settings.get("baseline", {})
        branches = settings.get("branches", {})

        if offset == 0:
            for parity in requested_parities:
                seq = baseline.get(parity)
                if seq:
                    key = _storage_key(seq, parity)
                    bundle[key] = _copy_sequence(seq)
            if bundle:
                return bundle
            continue

        use_simplified_oxidation = normalized_mode in {"deep", "own"} and offset > 0 and oxidation_baseline

        for target_parity in requested_parities:
            # For recursive trees: the tree structure encodes parity flips at each depth
            # Navigate from the baseline branch that leads to the target parity after N flips
            if normalized_mode in {"deep", "own"}:
                # The tree alternates parity at each depth level
                # To reach target parity at depth N, start from the branch that flips correctly
                # depth 0 (baseline) → depth 1 (flip once) → depth 2 (flip twice)
                # For even depth: same as baseline parity
                # For odd depth: opposite of baseline parity
                depth = abs(offset)
                if depth % 2 == 0:
                    # Even depth: same parity as target
                    source_parity = target_parity
                else:
                    # Odd depth: opposite parity from target
                    source_parity = "odd" if target_parity == "even" else "even"
            else:
                source_parity = target_parity

            parity_branches = branches.get(source_parity, {})
            if not parity_branches:
                if use_simplified_oxidation:
                    seq = oxidation_baseline.get(target_parity)
                    if seq:
                        key = _storage_key(seq, target_parity)
                        bundle[key] = _copy_sequence(seq)
                        return bundle
                continue

            # Collect preference chain for ALL tree modes
            # Use source_parity to find which FoB won at the previous step
            # Follow recorded preferences from anchor to target delta
            preference_chain = _collect_preference_chain(anchor, delta, source_parity, state_cache)
            anchor_preference = _preferred_index_from_state(state_cache, anchor, source_parity)
            preferred_branch = anchor_preference
            ordered_indices = _ordered_candidates(parity_branches, preferred_branch)

            for branch_index in ordered_indices:
                if is_recursive_tree:
                    # For recursive trees (deep), navigate through the tree
                    # offset can be ±1, ±2, ±3
                    # direction is +1 or -1, depth is abs(offset)
                    direction = +1 if offset > 0 else -1
                    depth = abs(offset)

                    branch_info = parity_branches.get(branch_index, {})
                    # Branch dictionaries can have int keys in-memory and string keys when loaded from JSON
                    first_level = branch_info.get(direction)
                    if not first_level:
                        first_level = branch_info.get(str(direction))
                    if not first_level:
                        continue

                    sub_chain: List[int] = []
                    if preference_chain and (preferred_branch is None or branch_index == preferred_branch):
                        sub_chain = preference_chain

                    if depth == 1:
                        pref_sub = sub_chain[0] if sub_chain else 1
                        seq = _extract_sequence_from_tree_node(first_level, preferred_sub_index=pref_sub)
                    else:
                        ordered_first = _ordered_candidates(first_level, sub_chain[0] if sub_chain else None)
                        seq = None
                        for sub_idx in ordered_first:
                            next_node = first_level.get(sub_idx)
                            if not next_node:
                                continue
                            # If there are no deeper preferences, keep following the same sub_idx path
                            chain_tail = sub_chain[1:] if sub_chain else [sub_idx]
                            seq = _navigate_recursive_tree(next_node, depth - 1, direction, chain_tail)
                            if seq:
                                break
                elif use_simplified_oxidation:
                    seq = oxidation_baseline.get(target_parity)
                    if not seq:
                        continue
                else:
                    # Non-recursive trees (flat, deep, deep2)
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

                key = _storage_key(seq, target_parity)
                bundle[key] = _copy_sequence(seq)
                if preferred_branch is None:
                    logger.debug(
                        "[occupier_auto] Using fallback branch %s for target_parity=%s (source=%s), anchor=%s, offset=%s",
                        branch_index,
                        target_parity,
                        source_parity,
                        anchor,
                        offset,
                    )
                break

        if bundle:
            return bundle

    return {}
