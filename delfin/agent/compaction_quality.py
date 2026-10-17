"""What a compaction lost, measured against a FIXED keep-list.

KG3 phase 2. ``compare_loss`` (compaction_log.py:119) checks four
categories and is log-only; this module completes the keep-list the
brief names — the person's goal, their latest instruction, constraints,
refusals with their reasons, open decisions, waits, the last test
result — as a PURE report over (history before, state block + summary
after). It is advisory only: it must never veto a compaction (the
context-length cliff outranks keeping — in_turn_compaction.py
"budget tiny -> still compacts").

An item is LOST when it is in the history before the compaction but in
neither the surviving state block nor the summary after it. Detection
is marker-based and deterministic: identical inputs give an identical
report. Machine turns ([Command results] and the other synthetic
prefixes) are inspected for the categories that live in tool output
(refusals, waits, test results); the person's goal and instructions
come from person and [System] rows only, so hostile tool output can
never fabricate them.
"""

from __future__ import annotations

from typing import Any, Iterable

# Machine-turn prefixes: the same set working_state.py treats as
# synthetic ([Command results], [System] is NOT synthetic there — it is
# a person/harness instruction — see working_state._machine_texts).
_MACHINE_PREFIXES = (
    "[Command results]",
    "[Message from the session",
    "[Error",
)

_CATEGORIES = ("goal", "instruction", "constraints", "refusals",
               "decisions", "waits", "tests")


def _text(msg: Any) -> str:
    content = msg.get("content", "") if isinstance(msg, dict) else ""
    if isinstance(content, str):
        return content
    try:
        return " ".join(str(b.get("text", ""))
                        for b in content if isinstance(b, dict))
    except Exception:
        return ""


def _norm(s: str) -> str:
    return " ".join(str(s).split())


def _is_machine(text: str) -> bool:
    return text.lstrip().startswith(_MACHINE_PREFIXES)


def _extractors() -> dict:
    """Marker -> category mapping. Markers are the concrete spellings
    DELFIN itself produces (refusals: 'Refused: ... Reason: ...';
    waits: 'waiting for' in a watch/schedule line; test results: a
    passed/failed count). Deterministic substring markers only."""
    return {
        "constraints": ["Constraint:", "constraint:", "read-only",
                        "never touch", "do not", "Do not", "must not"],
        "refusals": ["Refused:", "refused:", "Reason:"],
        "decisions": ["Open decision:", "open decision:",
                      "decided:", "Decided:"],
        "waits": ["waiting for", "Waiting for", "wait for",
                  "blocked on", "Blocked on"],
        "tests": ["passed", "failed", "passed,"],
    }


def quality_report(
    history_before: Iterable[dict],
    after_state_block: str = "",
    after_summary_text: str = "",
) -> dict:
    """Pure report: for each keep-list category, the items present in
    ``history_before`` but absent from both ``after_state_block`` and
    ``after_summary_text``. Never raises on odd inputs; never returns a
    'do not compact' verdict (advisory only)."""
    after = _norm(after_state_block) + "\n" + _norm(after_summary_text)
    after = _norm(after)
    report: dict[str, dict] = {}
    person: list[str] = []      # person rows (non-machine user rows)
    system_rows: list[str] = [] # [System] instruction rows
    machine: list[str] = []     # synthetic tool-result rows
    for msg in history_before or []:
        if not isinstance(msg, dict):
            continue
        role = msg.get("role")
        text = _norm(_text(msg))
        if not text:
            continue
        if role == "user":
            if text.startswith("[System]"):
                system_rows.append(text)
            elif _is_machine(text):
                machine.append(text)
            else:
                person.append(text)
        elif role == "assistant":
            machine.append(text)

    markers = _extractors()
    # goal: the person's newest non-machine user row.
    goals = [_norm(person[-1])[:300]] if person else []
    report["goal"] = _check("goal", goals, after)
    # instruction: [System] rows.
    report["instruction"] = _check(
        "instruction", [t[:300] for t in system_rows], after)
    for cat in ("constraints", "refusals", "decisions", "waits",
                "tests"):
        hits = _scan(person + machine, markers[cat])
        report[cat] = _check(cat, hits, after)
    return report


def _check(cat: str, items: list[str], after: str) -> dict:
    present = [i for i in items if _overlap(i, after)]
    lost = [i for i in items if i not in present]
    return {"category": cat, "items": items, "kept": present,
            "lost": lost}


def _scan(texts: list[str], markers: list[str]) -> list[str]:
    hits: list[str] = []
    for text in texts:
        for m in markers:
            idx = text.find(m)
            while idx != -1:
                frag = text[max(0, idx - 80): idx + 220].strip()
                if frag and frag not in hits:
                    hits.append(_norm(frag))
                idx = text.find(m, idx + 1)
    return hits[:20]


def _overlap(item: str, after: str) -> bool:
    """An item counts as kept when a meaningful chunk of it appears in
    the surviving text (not a single marker word). A synthetic machine
    prefix ([Command results], ...) is not part of the kept content —
    working_state blocks carry the bare output text — so it is stripped
    before matching."""
    item = _norm(item)
    for prefix in ("[Command results]", "[Error", "[Message from"):
        if item.startswith(prefix):
            item = _norm(item[len(prefix):])
    if not item:
        return True
    words = [w for w in item.split() if len(w) >= 4]
    if not words:
        return True
    kept = sum(1 for w in words if w in after)
    return kept >= max(1, int(len(words) * 0.6))


def category_names() -> tuple:
    return _CATEGORIES
