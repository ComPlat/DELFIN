"""Mutual waits: find the cycles in what sessions wait on.

Wave 17 measured it: a builder and its QS waited on each other for an
hour; idle sessions did not notice files written for them. This module is
the pure half of the answer (the emitting half lives in
session_health.watch_events): given every session's registered waits, find
the cycles — A waits on B AND B waits on A, and longer rings too — with
the sessions and what each waits for.

Input rows are wait dicts extended with the session that registered them
(``session`` key). Only waits whose ``on_whom`` names a SESSION form
cycles: a wait on "operator", "person" or "ci" cannot deadlock with
another session — a person answers. The match is on the session NAME in
``value``, the same address the message mailboxes use.
"""

from __future__ import annotations

#: A waiting_on-session edge only counts when on_whom says "session"
#: (an empty on_whom is a legacy row from before the field existed; the
#: value names a session either way, so it counts too — a legacy wait on
#: a session is still a wait on a session).
_SESSION_WHOM = ("", "session")


def _name(value) -> str:
    """The wait value as a one-line session name (hostile input must not
    crash the finder or grow it: newlines collapse, whitespace strips)."""
    text = " ".join(str(value or "").split())
    return text


def find(waits: list) -> list[list[dict]]:
    """The wait cycles in *waits*, each as the list of its edge rows.

    A cycle is a ring A -> B -> ... -> A where each arrow is one
    session's wait naming the next session (message_from waits on a
    session name, or any kind whose on_whom is "session"). Each cycle is
    reported once — the finder collects edges first, then finds rings in
    the edge set, so duplicate rows change nothing.

    Returns [] when no cycle exists (a chain is not a cycle).
    """
    edges: dict[str, set[str]] = {}
    rows: dict[tuple[str, str], dict] = {}
    for w in waits or []:
        if not isinstance(w, dict):
            continue
        session = _name(w.get("session"))
        target = _name(w.get("value"))
        whom = str(w.get("on_whom") or "").strip().lower()
        if not session or not target or session == target:
            continue
        if whom not in _SESSION_WHOM:
            continue  # operator/person/ci waits cannot deadlock sessions
        edges.setdefault(session, set()).add(target)
        rows.setdefault((session, target), w)

    cycles: list[list[dict]] = []
    seen: set[frozenset[str]] = set()
    for start in sorted(edges):
        # ring search from start back to start (DFS, bounded by the
        # edge set — every ring through start is found once)
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            for nxt in sorted(edges.get(node, ())):
                if nxt == start:
                    ring = frozenset(path)
                    if ring not in seen:
                        seen.add(ring)
                        cycles.append([rows[(path[i], path[(i + 1) % len(path)])]
                                       for i in range(len(path))])
                    continue
                if nxt in path or nxt not in edges:
                    continue
                stack.append((nxt, path + [nxt]))
    return cycles


def describe(cycles: list) -> list[str]:
    """One line per cycle for the display: who waits on whom, and what
    for. Markup-free by construction (names are whitespace-normalised,
    for_what lines come from register()'s one-line cleaner)."""
    lines = []
    for cyc in cycles or []:
        parts = []
        for w in cyc:
            what = str(w.get("for_what") or w.get("kind") or "")
            parts.append(f"{_name(w.get('session'))} waits on "
                         f"{_name(w.get('value'))}: {what}")
        lines.append("; ".join(parts))
    return lines
