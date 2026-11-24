"""Auto-selection rules for OCCUPIER sequences."""
from __future__ import annotations

import copy
import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from delfin.common.logging import get_logger

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
    """The recorded preferences, or none where the file is absent OR
    unreadable.

    Lenient on purpose: a reader (``_preferred_index_from_state`` and the
    state_cache read) asks what is known, and a damaged store must not crash
    the step that consumes it. It does not touch the file. A WRITER may not
    build on this answer -- see :func:`_load_state_for_write`.
    """
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


def _load_state_for_write(root: Path) -> Dict[str, Any]:
    """What the writer may modify, or raise: a damaged store is kept, not
    written over.

    A WRITER's read, distinct from :func:`_load_state`, which may answer a
    damaged file with "nothing is recorded" -- true for a reader and
    destructive for a writer that would rebuild the file from it and drop
    every preference it held. On damage the file is moved to the dated
    ``.unreadable-`` quarantine of ``agent/state_file.py`` (bytes kept) and
    :class:`agent.state_file.StateFileUnreadable` is raised: a store that
    cannot be read cannot be safely replaced, so the write is refused
    rather than proceeding to empty it. The store is gone, so the next
    record writes a fresh one safely.
    """
    from delfin.agent import state_file as _sf
    data, moved = _sf.read_json_for_write(_state_path(root), default=None)
    if moved is not None:
        raise _sf.StateFileUnreadable(
            f"{_state_path(root)} cannot be read; it was kept at {moved}. "
            f"Refusing to write over it -- open it and copy out what "
            f"you need, or delete it to start fresh.")
    return data if isinstance(data, dict) else {}


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
                           populated: Optional[List[Dict[str, Any]]] = None,
                           root: Optional[Path] = None) -> None:
    """Remember which configuration won for a given delta, and which coexist.

    The winner drives the next redox step. `populated` carries every
    configuration whose Boltzmann weight cleared BOLTZMANN_PARENT_THRESHOLD,
    including the winner: where a molecule is a thermal mixture rather than a
    single state, each partner is a starting point the next step has to be
    derived from. Older state files carry no such list and fall back to the
    winner alone.
    """
    if preferred_index is None:
        return
    root_path = _resolve_root(root)
    # A writer: a damaged store is kept (quarantine) and the record
    # refused with StateFileUnreadable, not emptied by writing a store
    # that holds only this delta's preference.
    state = _load_state_for_write(root_path)
    entry = state.setdefault(str(delta), {})

    # Store full winner info for own mode rule-based generation
    winner_info: Dict[str, Any] = {"index": int(preferred_index)}
    if m_value is not None:
        winner_info["m"] = int(m_value)
    if bs_value is not None:
        winner_info["BS"] = str(bs_value).strip()

    if populated:
        kept = [p for p in populated if isinstance(p, dict) and p.get("m") is not None]
        if len(kept) > 1:
            winner_info["populated"] = kept
            logger.info(
                "[occupier_auto] delta=%s is a thermal mixture; deriving the next step "
                "from %d configurations: %s",
                delta, len(kept),
                ", ".join(
                    f"m={p.get('m')}{' BS ' + str(p['BS']) if p.get('BS') else ''}"
                    f" ({100 * float(p.get('weight', 0.0)):.0f}%)"
                    for p in kept
                ),
            )

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


#: Minimum Boltzmann population for a configuration to seed the next redox
#: step. At 298.15 K, kT is 0.593 kcal/mol; 5 % against the majority is a free
#: energy gap of about 1.7 kcal/mol, roughly 3 kT. Below that a state is not
#: meaningfully present in the equilibrium and following it would only cost
#: calculations; at or above it the substance genuinely is a mixture, and the
#: oxidised or reduced species can form from either partner.
BOLTZMANN_PARENT_THRESHOLD = 0.05


def _entry_key(entry: Dict[str, Any]) -> Tuple[Any, str]:
    """Identity of a configuration: its multiplicity and broken-symmetry label."""
    return entry.get("m"), str(entry.get("BS") or "").strip()


def _remap_from(value: Any, local: Dict[int, int]) -> Any:
    """Rewrite a `from` reference onto merged indices.

    `from` names the configuration whose converged orbitals this one starts
    from, which is what makes a broken-symmetry solution findable at all
    instead of collapsing back to the pure state. Merging sequences renumbers
    everything, so these have to follow. 0 means "start from the input
    geometry" and stays 0.
    """
    if isinstance(value, (list, tuple)):
        return [local.get(int(v), 0) for v in value if str(v).lstrip("-").isdigit()]
    try:
        index = int(value)
    except (TypeError, ValueError):
        return 0
    return local.get(index, 0) if index else 0


def merge_sequences(sequences: List[List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """Combine per-parent sequences into one, without computing anything twice.

    Two thermally populated parents usually propose overlapping daughters —
    both may want the same pure multiplicity — and running that twice would
    cost an optimisation for nothing. Configurations are therefore identified
    by (multiplicity, BS label): the first occurrence is kept, later ones fold
    into it, and any `from` pointing at a folded entry is redirected to the
    survivor.
    """
    merged: List[Dict[str, Any]] = []
    by_key: Dict[Tuple[Any, str], int] = {}

    for sequence in sequences:
        local: Dict[int, int] = {}
        added: List[Tuple[Dict[str, Any], Any]] = []
        for entry in sequence:
            try:
                old_index = int(entry["index"])
            except (KeyError, TypeError, ValueError):
                continue
            key = _entry_key(entry)
            if key in by_key:
                local[old_index] = by_key[key]
                continue
            new_entry = dict(entry)
            new_entry["index"] = len(merged) + 1
            merged.append(new_entry)
            by_key[key] = new_entry["index"]
            local[old_index] = new_entry["index"]
            added.append((new_entry, entry.get("from", 0)))
        for new_entry, original_from in added:
            new_entry["from"] = _remap_from(original_from, local)

    return merged


def populated_parents(winner_info: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """The configurations the next redox step should be derived from.

    Normally one: the state that won. When several were populated enough to
    coexist, all of them — a substance sitting as a 50:50 mixture of two spin
    states can be oxidised out of either, and the daughter states reachable
    from one are not the daughter states reachable from the other. Taking only
    the lower parent quietly assumes a purity the molecule does not have.
    """
    if not winner_info:
        return []
    populated = winner_info.get("populated")
    if isinstance(populated, list) and populated:
        usable = [p for p in populated if isinstance(p, dict) and p.get("m") is not None]
        if usable:
            return usable
    return [winner_info] if winner_info.get("m") is not None else []


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

    # One parent normally, several when the previous step left a thermal
    # mixture. Each is carried through the same rules and the results are
    # merged, so a daughter both parents reach is still computed once.
    parents = populated_parents(winner_info)
    sequences: List[List[Dict[str, Any]]] = []
    for parent in parents:
        parent_seq = generate_sequence_from_winner_rules(
            even_baseline,
            odd_baseline,
            parent["m"],
            parent.get("BS", "") or "",
            target_parity,
            pure_window=pure_window,
            progressive_from=progressive_from,
            include_bs=include_bs,
        )
        if parent_seq:
            sequences.append(parent_seq)

    seq = merge_sequences(sequences)

    if len(parents) > 1:
        logger.info(
            "[own_mode] delta=%s derived from %d populated configurations (%s) → "
            "%d candidates for parity=%s",
            delta, len(parents),
            ", ".join(
                f"m={p['m']}{' BS ' + str(p['BS']) if p.get('BS') else ''}" for p in parents
            ),
            len(seq), target_parity,
        )
    else:
        logger.debug(
            "[own_mode] Winner from delta=%s was parity=%s (m=%s, BS=%s), generating parity=%s",
            prev_delta, source_parity, winner_info.get("m"), winner_info.get("BS", ""),
            target_parity,
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
                                 parity_hint: Optional[str] = None,
                                 custom_dataset: Optional[Dict[int, Dict[str, Any]]] = None,
                                 config: Optional[Dict[str, Any]] = None) -> Dict[str, List[Dict[str, Any]]]:
    """Return auto-managed sequences for the requested delta (if available).

    Sequences are derived from rules: the configurations the previous step left
    populated, put through the pure-window and broken-symmetry rules, merged.
    The pre-built decision trees this used to be able to navigate instead are
    gone — they were the same decision enumerated ahead of time, and no
    calculation in the run archive ever selected one.
    """
    if not custom_dataset:
        return {}
    return _resolve_own_mode_sequences(delta, custom_dataset, root=root, config=config)
# ---------------------------------------------------------------------
# Spin states from the metal and its d-electron count.
#
# MEASURED FIRST (rule 9). The run archive holds 1906 DELFIN_Data.json
# records, 619 of them with OCCUPIER scans, 1316 OCCUPIER steps in all
# (initial 302, red_step_1 301, red_step_2 301, ox_step_1 294, red_step_3
# 128). Population: DELFIN OCCUPIER runs under the two archive roots, as
# they stood when this was written.
#
# What the archive DOES carry, per metal and row (steps / winner was the
# lowest scanned multiplicity / the highest / an intermediate one):
#   Fe|3d  654 / 126 / 99 / 429      Co|3d  311 / 192 /  8 / 111
#   Ni|3d  301 / 260 /  3 /  38      Mn|3d   45 /   0 /  8 /  37
#   Cu|3d    5 /   3 /  0 /   2      Zn|3d   10 /  10 /  0 /   0
# Read: where the winner is decided at all, it is the EXTREME of the
# scanned set far more often than an intermediate multiplicity -- Ni 260
# of 301, Co 192 of 311, Zn 10 of 10. That is why the proposal below
# offers the two ligand-field extremes and no middle states: a scan that
# omits an extreme cannot find a winner there, and the archive shows the
# extremes are where the winner usually is.
#
# What the archive does NOT carry: the metal's oxidation state, and the
# ligand set in a parseable form. It records each step's NET charge, which
# includes the ligands, so d-count = group - OS is NOT derivable from it.
# Deriving OS from the net charge would invent the number the rule rests
# on, so this module does not do it, and no row below cites an archive
# d-count. The rows cite ligand-field theory, which is what DELFIN's own
# occupier skill already teaches.
#
# This block ADDS a proposal. It never changes the sequence chosen when no
# metal is given: ``propose_spin_states_for_metal`` returns None and every
# caller keeps today's parity rule unchanged.
# ---------------------------------------------------------------------

# Group number per transition metal: d-count = group - oxidation state.
_METAL_GROUP: Dict[str, int] = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10,
    "Cu": 11, "Zn": 12, "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8,
    "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12, "Hf": 4, "Ta": 5, "W": 6, "Re": 7,
    "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
}

# 4d and 5d centres: Delta_o is large, so low spin is the reference state
# and a high-spin-only scan misses it.
_SECOND_AND_THIRD_ROW = frozenset({
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
})

_EVIDENCE_LIGAND_FIELD = (
    "ligand-field theory: a 4d/5d centre has a large Delta_o, so the "
    "low-spin occupancy is the reference and the high-spin one is the "
    "comparison (occupier skill, orca-spectra-spin-occupier)"
)
_EVIDENCE_FIRST_ROW = (
    "ligand-field theory: for a first-row centre both the low-spin and the "
    "high-spin occupancy must be computed -- Delta_o is not predictable "
    "from the d-count alone, so neither may be assumed "
    "(occupier skill, orca-spectra-spin-occupier); the archive's winners "
    "are the extremes of the scanned set in 260/301 Ni, 192/311 Co steps"
)


def _unpaired_in_orbitals(electrons: int, orbitals: int) -> int:
    """Unpaired electrons when ``electrons`` fill ``orbitals`` orbitals."""
    if electrons <= orbitals:
        return electrons
    return max(0, 2 * orbitals - electrons)


def _multiplicities_for_dcount(dcount: int) -> List[int]:
    """The multiplicity extremes worth computing for a d^n centre.

    Octahedral occupancy, the textbook construction: the low-spin state
    fills t2g (3 orbitals, 6 electrons) before eg (2 orbitals); the
    high-spin state fills all five d orbitals singly before pairing.
    Sorted ascending, and de-duplicated -- for d0..d3, d9..d10 both
    extremes coincide, which is correct, not a missing entry.

    The parity of the d-count is a hard constraint that this obeys by
    construction: an odd d-count cannot have M=1. d5 gives [2, 6], never
    a multiplicity-1 entry.

    d8 is the one row the octahedral construction gets wrong, so it is
    stated explicitly. Octahedral d8 is t2g^6 eg^2 (2 unpaired, M=3) --
    but a d8 centre is very often SQUARE PLANAR, where the four lower
    orbitals fill and the dx2-y2 pair is empty: diamagnetic, M=1. That is
    the classic Ni(II) case, and square planar is the geometry a d8
    complex usually adopts. A proposal of [3] alone would never compute
    the diamagnetic state, so the scan could not find it. Both go in.

    Evidence, measured from the run archive. The population is the 301
    OCCUPIER steps whose system contains nickel ("Ni|3d" in
    scratch/population.json), and their winner multiplicities are
    M=1 111, M=2 149, M=3 35, M=4 4, M=5 2. That row is per ELEMENT, and
    an element is not a d-count: 153 of those winners sit at an even
    multiplicity, which a d8 centre cannot have (8 d electrons leave an
    even number unpaired, so the multiplicity is odd). Restricting to the
    148 d8-CONSISTENT steps, M=1 won 111 (75 %) and M=3 won 35. So the
    closed-shell singlet is a real winner for a nickel centre, and it
    must be in the scanned set -- but the number that says so is 111 of
    148, not 111 of 301, and the d-count itself is inferred from nickel,
    never recorded by the archive.
    """
    n = max(0, min(10, int(dcount)))
    if n == 8:
        return [1, 3]
    in_t2g = min(n, 6)
    in_eg = max(0, n - 6)
    n_unpaired_ls = _unpaired_in_orbitals(in_t2g, 3) + _unpaired_in_orbitals(in_eg, 2)
    n_unpaired_hs = n if n <= 5 else 10 - n
    return sorted({n_unpaired_ls + 1, n_unpaired_hs + 1})


def propose_spin_states_for_metal(
    metal: Optional[str],
    oxidation_state: Optional[int],
    *,
    dcount: Optional[int] = None,
) -> Optional[Dict[str, Any]]:
    """Multiplicities and BS states to compute for one metal centre.

    Returns ``None`` when the metal or its oxidation state cannot be
    determined -- the caller then keeps today's parity rule unchanged.
    Returns ``{metal, oxidation_state, dcount, multiplicities, first,
    broken_symmetry, reason, evidence}`` otherwise.

    ``first`` is the multiplicity to try first: for a 4d/5d centre that is
    the low-spin one (large Delta_o), for a first-row centre the high-spin
    one (the ligand-field strength is what decides, so the scan is what
    answers it, not a guess). ``broken_symmetry`` is True where an
    antiferromagnetically coupled pair is worth starting: an even d-count
    with at least two unpaired electrons, where the low-spin closed shell
    and the BS open shell are different wave functions for the same
    multiplicity.

    This proposes only. It is never written into a spec: a caller's own
    sequence outranks it, and the scan itself decides the winner.
    """
    if not metal or oxidation_state is None:
        return None
    symbol = str(metal).strip()
    symbol = symbol[:1].upper() + symbol[1:].lower()
    group = _METAL_GROUP.get(symbol)
    if group is None:
        return None
    n = group - int(oxidation_state) if dcount is None else int(dcount)
    if not 0 <= n <= 10:
        return None
    heavy = symbol in _SECOND_AND_THIRD_ROW
    mults = _multiplicities_for_dcount(n)
    # Which one to try first. A heavy centre's low-spin occupancy is the
    # reference (large Delta_o), so M=1 leads. d8 is the other row where
    # the answer is not "high spin first": square-planar d8 is
    # diamagnetic, and the archive agrees -- of the 148 d8-consistent
    # nickel steps, M=1 won 111 and M=3 won 35. Keeping the 3d default
    # here would have the proposal say "M=3 first" while its own evidence
    # says the singlet usually wins. Every other 3d centre is high-spin
    # first, because Delta_o is not predictable from a d-count.
    if heavy or n == 8:
        first = mults[0]
    else:
        first = mults[-1]
    # Broken symmetry is proposed for d6 only, and the reason is narrow:
    # it is the even d-count whose low-spin state is a CLOSED shell
    # (t2g^6, M=1) while its high-spin state carries 4 unpaired electrons
    # (M=5), so the closed-shell singlet and a broken-symmetry open-shell
    # singlet at the same multiplicity are genuinely different wave
    # functions and only the BS calculation tells them apart. It is NOT
    # proposed from the d-count for a dinuclear or mixed-valence complex:
    # there the number of magnetic centres decides it, not the metal --
    # the skill's `check_contamination` branch is what inserts BS then.
    bs = n == 6
    row = "4d/5d" if heavy else "3d"
    # The reason must name the SAME state `first` names, or the proposal
    # contradicts itself: the earlier text said "which one wins is what
    # the scan decides" even where `first` already picks one, and for d8
    # that read as M=3-first while `first` was 1.
    if heavy:
        lead = "try M=1 first: a heavy centre's large Delta_o makes the low-spin occupancy the reference"
    elif n == 8:
        lead = "try M=1 first: square-planar d8 is diamagnetic, and the archive's d8-consistent nickel steps put M=1 ahead of M=3"
    else:
        lead = (f"try M={first} first, but which one wins is what the scan "
                "decides: a first-row centre's Delta_o is not predictable "
                "from the d-count alone")
    # d0..d3 and d9..d10 collapse to ONE extreme (both coincide); naming it
    # twice reads as "compute M=2 and M=2", which is not a comparison.
    scan = (f"compute M={mults[0]}" if len(set(mults)) == 1
            else f"compute M={mults[0]} and M={mults[-1]}")
    reason = (
        f"{symbol}({int(oxidation_state):+d}) is d{n} on a {row} centre: "
        f"{scan}; {lead}"
        + (", plus a broken-symmetry open-shell singlet" if bs else "")
    )
    return {
        "metal": symbol,
        "oxidation_state": int(oxidation_state),
        "dcount": n,
        "row": row,
        "multiplicities": mults,
        "first": first,
        "broken_symmetry": bs,
        "reason": reason,
        "evidence": [_EVIDENCE_LIGAND_FIELD if heavy else _EVIDENCE_FIRST_ROW],
        "fallback": False,
    }


# Element symbols of the metals above, longest first so a two-letter
# symbol is preferred over a one-letter one at the same position.
_METAL_NAMES: Dict[str, str] = {
    "iron": "Fe", "cobalt": "Co", "nickel": "Ni", "manganese": "Mn",
    "copper": "Cu", "chromium": "Cr", "vanadium": "V", "titanium": "Ti",
    "ruthenium": "Ru", "rhodium": "Rh", "palladium": "Pd",
    "osmium": "Os", "iridium": "Ir", "platinum": "Pt", "molybdenum": "Mo",
    "tungsten": "W", "zinc": "Zn", "scandium": "Sc", "yttrium": "Y",
    "zirconium": "Zr", "niobium": "Nb", "silver": "Ag", "gold": "Au",
    "cadmium": "Cd", "mercury": "Hg", "hafnium": "Hf", "tantalum": "Ta",
    "rhenium": "Re", "technetium": "Tc",
}
_ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7,
          "VIII": 8}

# A metal symbol must stand alone as a token: not preceded by a letter and
# not followed by a lowercase letter. Without this, "a co-solvent study"
# and "a Ti(IV) centre in a co-factor" are read as cobalt and titanium
# complexes -- measured on the alternative implementation before this one.
_SYMBOL_ALT = "|".join(sorted(_METAL_GROUP, key=len, reverse=True))
_OXIDE_SUFFIX = r"(?![A-Za-z])"
_METAL_OX_RE = re.compile(
    r"(?<![A-Za-z])(" + _SYMBOL_ALT + r")\s*\(\s*(?:([IVX]+)|([0-9])\+?)\s*\)"
    + _OXIDE_SUFFIX
    + r"|(?<![A-Za-z])(" + _SYMBOL_ALT + r")\s*([0-9])\+" + _OXIDE_SUFFIX,
)


def parse_metal_request(text: str) -> Optional[Dict[str, Any]]:
    """Read a metal centre and its oxidation state out of a request.

    Only an explicitly written centre is honoured -- ``Fe(III)``,
    ``Fe(3+)``, ``Fe3+``, ``iron(III)``. A request that names no metal, or
    names one without an oxidation state, returns ``None`` so the parity
    rule stays in charge (it is the fallback, not this).
    """
    if not text:
        return None
    body = str(text)
    for word, symbol in _METAL_NAMES.items():
        m = re.search(
            r"(?<![A-Za-z])" + re.escape(word)
            + r"\s*\(\s*(?:([IVX]+)|([0-9])\+?)\s*\)" + _OXIDE_SUFFIX,
            body, re.IGNORECASE)
        if m:
            os_val = _ROMAN.get(m.group(1).upper()) if m.group(1) else int(m.group(2))
            if os_val is None:
                return None
            return {"metal": symbol, "oxidation_state": os_val}
    m = _METAL_OX_RE.search(body)
    if m:
        if m.group(2):
            os_val = _ROMAN.get(m.group(2).upper())
        elif m.group(3):
            os_val = int(m.group(3))
        else:
            os_val = int(m.group(5))
        if os_val is None:
            return None
        return {"metal": m.group(1) or m.group(4), "oxidation_state": os_val}
    return None


def propose_sequence_for_request(text: str) -> Optional[Dict[str, Any]]:
    """The spin-state proposal for a natural-language request, or None.

    Deliberately narrow: a proposal only when the request states BOTH a
    metal and its oxidation state. Everything else -- a closed-shell
    organic, a request naming a metal without a charge -- falls through to
    the parity rule, unchanged.
    """
    parsed = parse_metal_request(text)
    if not parsed:
        return None
    return propose_spin_states_for_metal(parsed["metal"],
                                         parsed["oxidation_state"])
