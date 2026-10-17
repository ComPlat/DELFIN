"""prompt_head: cross-turn view of the request head.

prompt_loader.compose_sections already orders a single request head
correctly (stable -> route -> volatile -> anchor). What is missing is a
CHECK over two consecutive turns: did anything that changes between
turns sneak in BEFORE the stable prefix ends, so a provider prefix
cache misses? This module is pure: it takes the section rows of two
recorded turns and reports the shared prefix, and it checks the safety
ordering invariants (fixed layer order, critical_anchor last, the
person's newest instruction present, no fenced tool output inside a
stable-layer row). It reads nothing, writes nothing, never raises.
"""

from __future__ import annotations

from typing import Any, Iterable

# Layer numbers must mirror prompt_loader's LAYER_* constants
# (stable=0, route=1, volatile=2, anchor=3); the ordering contract is
# "layers in nondecreasing order, anchor last".
MAX_LAYER = 3

# Fenced/synthetic markers a row's CONTENT may carry when tool output
# leaked into it (untrusted.wrap style). Such text inside a stable row
# means hostile material reached the cacheable head.
_HOSTILE_MARKERS = (
    "[data]", "[UNTRUSTED EXTERNAL CONTENT", "UNTRUSTED",
    "instruction-shaped", "not quoted: the output carries")


def _row_get(row: dict, key: str, default: Any = None) -> Any:
    if not isinstance(row, dict):
        return default
    value = row.get(key, default)
    return value if value is not None else default


def stable_prefix_len(
    turn_a: Iterable[dict], turn_b: Iterable[dict]
) -> dict:
    """How many leading rows/bytes the two recorded turns share.

    The shared prefix is what a provider prefix cache can hit. Rows
    with layer < 2 (stable, route) are EXPECTED to be inside it; a
    change at one of those rows is the finding, so it is named.
    """
    a = list(turn_a or [])
    b = list(turn_b or [])
    rows = 0
    bytes_ = 0
    first_change_row: int | None = None
    first_change_section: str | None = None
    for i in range(min(len(a), len(b))):
        ra, rb = a[i], b[i]
        if (_row_get(ra, "name") == _row_get(rb, "name")
                and _row_get(ra, "content") == _row_get(rb, "content")):
            rows += 1
            bytes_ += len(str(_row_get(ra, "content") or ""))
        else:
            first_change_row = i
            first_change_section = str(_row_get(ra, "name") or "?")
            break
    else:
        if len(a) != len(b):
            first_change_row = min(len(a), len(b))
            first_change_section = (
                str(_row_get((a if len(a) < len(b) else b)[
                    first_change_row], "name") or "?"))
    expected = [r for r in a[:rows] if _row_get(r, "layer", 0) is not None
                and int(_row_get(r, "layer", 0)) < 2]
    # "all stable" means: no STABLE (layer 0) row sits after the shared
    # prefix — route rows may legitimately change per turn (cycle
    # counters); a stable row after the first change is the finding.
    tail_unstable = all(
        _row_get(r, "layer", 0) is not None
        and int(_row_get(r, "layer", 0)) >= 1 for r in a[rows:])
    return {
        "rows": rows,
        "bytes": bytes_,
        "rows_all_stable": tail_unstable,
        "first_change_row": first_change_row,
        "first_change_section": first_change_section,
        "stable_rows_in_prefix": len(expected),
    }


def check_ordering(
    rows: Iterable[dict], require_person_instruction: bool = False
) -> list[str]:
    """Safety invariants of one request head. Returns human-readable
    problems (empty list = compliant); never raises."""
    rows = list(rows or [])
    problems: list[str] = []
    last_layer = 0
    for i, r in enumerate(rows):
        layer = _row_get(r, "layer", 0)
        try:
            layer = int(layer)
        except (TypeError, ValueError):
            problems.append(f"row {i}: non-numeric layer")
            continue
        if layer < last_layer:
            problems.append(
                f"row {i} ({_row_get(r, 'name', '?')}): layer {layer} "
                f"after layer {last_layer} — fixed order broken")
        last_layer = layer
    if rows:
        last = rows[-1]
        if str(_row_get(last, "name", "")) != "critical_anchor":
            problems.append(
                "last row is not critical_anchor "
                f"(got {str(_row_get(last, 'name', '?'))!r}) — "
                "reordering dropped the safety anchor from its "
                "recency position")
    if require_person_instruction:
        found = newest_person_instruction(rows)
        if not found:
            problems.append(
                "no person instruction found in the head — "
                "the person's newest instruction must not be "
                "compacted away")
    return problems


def newest_person_instruction(rows: Iterable[dict]) -> str:
    """The content of the newest person_instruction row, or ""."""
    result = ""
    for r in rows or []:
        if (str(_row_get(r, "name", "")) == "person_instruction"):
            text = str(_row_get(r, "content", "") or "")
            if text.strip():
                result = text.strip()
    return result


def hostile_reach(rows: Iterable[dict]) -> list[str]:
    """Stable-layer rows whose content carries fenced/synthetic tool
    output markers — hostile tool output must not reach the cacheable
    head. Volatile rows are expected to hold live state and are exempt."""
    problems: list[str] = []
    for i, r in enumerate(rows or []):
        layer = _row_get(r, "layer", 0)
        try:
            layer = int(layer)
        except (TypeError, ValueError):
            continue
        if layer >= 2:
            continue  # volatile/anchor legitimately carry live state
        content = str(_row_get(r, "content", "") or "")
        for marker in _HOSTILE_MARKERS:
            if marker in content:
                problems.append(
                    f"row {i} ({str(_row_get(r, 'name', '?'))!r}): "
                    f"stable row carries fenced/synthetic marker "
                    f"{marker!r} — hostile material in the cacheable "
                    "head")
                break
    return problems
