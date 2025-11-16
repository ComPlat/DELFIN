"""Auto-selection rules for OCCUPIER sequences."""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.common.logging import get_logger
from delfin.deep_auto_tree import DEEP_AUTO_SETTINGS
from delfin.deep2_auto_tree import DEEP2_AUTO_SETTINGS
from delfin.deep3_auto_tree import DEEP3_AUTO_SETTINGS
from delfin.deep4_auto_tree import DEEP4_AUTO_SETTINGS
from delfin.deep5_auto_tree import DEEP5_AUTO_SETTINGS
from delfin.deep6_auto_tree import DEEP6_AUTO_SETTINGS

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


def _preferred_index_from_state(state: Dict[str, Any], delta: int,
                                parity: Optional[str] = None) -> Optional[int]:
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
            return int(value)
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
    for step in range(steps):
        delta_value = anchor + direction * step
        remaining = steps - step
        parity_token = _parity_with_distance(parity, remaining)
        pref = _preferred_index_from_state(state, delta_value, parity_token)
        if pref is None:
            break
        chain.append(pref)
    return chain


def _ordered_candidates(source: Dict[int, Any], preferred: Optional[int]) -> List[int]:
    ordered: List[int] = []
    if preferred is not None and preferred in source:
        ordered.append(preferred)
    for candidate in sorted(source.keys()):
        if candidate not in ordered:
            ordered.append(candidate)
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

    For deep3/deep4, structure is: node['branches'][direction][sub]['branches'][direction][sub]...
    remaining_steps=0 means we're at the target node
    """
    if remaining_steps == 0:
        # We're at the target node, extract its sequence
        if "seq" in node:
            return node.get("seq")
        return None

    # Navigate one level deeper
    inner_branches = node.get("branches", {})
    if not inner_branches or direction not in inner_branches:
        return None

    # Get the next level (dict of sub-indices)
    next_level = inner_branches[direction]
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


_TREE_DATASETS = {
    "deep": AUTO_SETTINGS,
    "flat": AUTO_SETTINGS_FLAT,
    "deep2": DEEP2_AUTO_SETTINGS,
    "deep3": DEEP3_AUTO_SETTINGS,
    "deep4": DEEP4_AUTO_SETTINGS,
    "deep5": DEEP5_AUTO_SETTINGS,
    "deep6": DEEP6_AUTO_SETTINGS,
}


def resolve_auto_sequence_bundle(delta: int, *, root: Optional[Path] = None,
                                 tree_mode: str = "deep",
                                 parity_hint: Optional[str] = None) -> Dict[str, List[Dict[str, Any]]]:
    """Return auto-managed sequences for the requested delta (if available)."""
    normalized_mode = str(tree_mode or "deep").strip().lower()
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

    # Check if this is a recursive tree (deep3/deep4/deep5)
    is_recursive_tree = normalized_mode in {"deep3", "deep4", "deep5"}

    oxidation_baseline = None
    if normalized_mode == "deep6":
        base_settings = settings_source.get(0)
        if base_settings:
            oxidation_baseline = base_settings.get("baseline", {})

    for anchor, settings in settings_source.items():
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

        use_simplified_oxidation = normalized_mode == "deep6" and offset > 0 and oxidation_baseline

        for parity in requested_parities:
            parity_branches = branches.get(parity, {})
            if not parity_branches:
                if use_simplified_oxidation:
                    seq = oxidation_baseline.get(parity)
                    if seq:
                        key = _storage_key(seq, parity)
                        bundle[key] = _copy_sequence(seq)
                        return bundle
                continue

            if normalized_mode == "deep6":
                preferred_branch = _preferred_index_from_state(state_cache, anchor, parity)
                preference_chain: List[int] = [preferred_branch] if preferred_branch else []
            else:
                preference_chain = _collect_preference_chain(anchor, delta, parity, state_cache)
                preferred_branch = preference_chain[0] if preference_chain else None
            ordered_indices = _ordered_candidates(parity_branches, preferred_branch)

            for branch_index in ordered_indices:
                if is_recursive_tree:
                    # For recursive trees (deep3/deep4), navigate through the tree
                    # offset can be ±1, ±2, ±3
                    # direction is +1 or -1, depth is abs(offset)
                    direction = +1 if offset > 0 else -1
                    depth = abs(offset)

                    branch_info = parity_branches.get(branch_index, {})
                    first_level = branch_info.get(direction)
                    if not first_level:
                        continue

                    sub_chain: List[int] = []
                    if preferred_branch is not None and branch_index == preferred_branch:
                        sub_chain = preference_chain[1:]

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
                            chain_tail = sub_chain[1:] if sub_chain else []
                            seq = _navigate_recursive_tree(next_node, depth - 1, direction, chain_tail)
                            if seq:
                                break
                elif use_simplified_oxidation:
                    seq = oxidation_baseline.get(parity)
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

                key = _storage_key(seq, parity)
                bundle[key] = _copy_sequence(seq)
                if preferred_branch is None:
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
