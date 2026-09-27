"""Learn a reusable skill from one finished DELFIN session.

Hermes-style trigger: a session that earned a skill proposal produces at
most ONE proposal — and only when it also carries evidence (a green test
run, a verified calculation/job, or a verify_recipe). A session without
evidence never becomes a rule; an unsafe success stays a suggestion at
most (skill_proposals.propose runs the safety check).

Message forms handled here are the ones DELFIN transcripts actually
contain: ``chat_messages`` dicts with roles user/system/thinking/assistant
and ``tool`` messages whose content embeds an HTML tool chip
(``<span class="tool-name">NAME</span>``) and a result block
(``<details><summary> -> {json}</details>``). Parsing is tolerant: any
message dict with a ``name``/``tool``/``tool_name`` key or an embedded
tool chip counts as a tool call.
"""

from __future__ import annotations

import re

_TOOL_CHIP = re.compile(r'tool-name">([^<]+)<')
_RESULT_BLOCK = re.compile(r"<summary>\s*(?:&rarr;|->|→)\s*(.*?)</details>",
                           re.S)

# Tool calls inside one working stretch (between user turns) that mark
# deep, non-obvious work.
_DEEP_STRETCH_CALLS = 5


def _content(msg: dict) -> str:
    """Best-effort text of one message; never raises."""
    try:
        c = msg.get("content")
        return c if isinstance(c, str) else str(c or "")
    except Exception:
        return ""


def _tool_name(msg: dict) -> str | None:
    """Name of the tool this message represents, or None."""
    for key in ("name", "tool", "tool_name"):
        v = msg.get(key)
        if isinstance(v, str) and v:
            return v
    if msg.get("role") in ("tool", "function"):
        m = _TOOL_CHIP.search(_content(msg))
        if m:
            return m.group(1)
        return "tool"
    m = _TOOL_CHIP.search(_content(msg))
    if m and msg.get("role") in (None, "tool"):
        return m.group(1)
    return None


def _result_text(msg: dict) -> str:
    m = _RESULT_BLOCK.search(_content(msg))
    return m.group(1) if m else ""


def _is_tool_failed(msg: dict, name: str | None) -> bool:
    """Did this tool result report failure? Tolerant, best-effort."""
    text = (_content(msg) + " " + _result_text(msg)).lower()
    m = re.search(r'"exit_code":\s*(-?\d+)', text)
    if m and m.group(1) not in ("0", ""):
        return True
    if re.search(r'"(status|state)":\s*"(error|failed|denied)"', text):
        return True
    if re.search(r'\b(failed|error)\b', text) and "passed" not in text:
        # "1 failed" style summaries; a summary that also says "passed"
        # (e.g. "12 passed, 1 failed") is handled above only via counts —
        # keep it conservative: count-based check first.
        if not re.search(r'"failed":\s*0\b', text):
            return True
    if "traceback" in text or "modulenotfounderror" in text:
        return True
    if name and name.lower() in ("read", "edit", "write") and "refused" in text:
        return True
    return False


def _is_tool_green(msg: dict, name: str | None) -> bool:
    """Did this tool result report clear success?"""
    text = (_content(msg) + " " + _result_text(msg)).lower()
    m = re.search(r'"exit_code":\s*(\d+)', text)
    if m and m.group(1) == "0":
        return True
    if re.search(r'"failed":\s*0\b.*"passed":\s*[1-9]', text):
        return True
    if re.search(r'"(?:passed|updated|created|written)"\s*:\s*[1-9]', text):
        return True
    if name and name.lower() == "toolrunner" and '"failed": 0' in text:
        return True
    return False


def _looks_like_correction(text: str) -> bool:
    """Does a user message read as a course correction?"""
    t = text.lower()
    markers = (
        "no, ", "nein", "falsch", "wrong", "not what i meant",
        "i meant", "stattdessen", "instead", "revert", "undo",
        "please revert", "undo that", "rollback", "that's not",
        "das ist falsch", "das war falsch", "ich meinte",
        "nicht so", "anders", "stop, ", "halt,",
    )
    return any(m in t for m in markers)


def qualifies(messages) -> list[str]:
    """Named reasons why this session earned a skill proposal.

    Empty list means: no proposal. One entry per trigger, each a short
    named English reason. Never raises on odd input.
    """
    try:
        msgs = [m for m in (messages or []) if isinstance(m, dict)]
    except Exception:
        return []
    if len(msgs) < 3:
        return []

    reasons: list[str] = []
    user_msgs = [m for m in msgs if m.get("role") == "user"]

    # 1) Deep work: >= 5 tool calls in one working stretch (between two
    #    user turns, or from the first user turn to the end).
    stretches: list[list[dict]] = [[]]
    for m in msgs:
        if m.get("role") == "user" and user_msgs and m is not user_msgs[0]:
            stretches.append([])
        stretches[-1].append(m)
    for stretch in stretches:
        n = sum(1 for m in stretch if _tool_name(m))
        if n >= _DEEP_STRETCH_CALLS:
            reasons.append(
                f"deep-work: {n} tool calls in one working stretch")
            break

    # 2) Error followed by recovery: a tool result reporting failure,
    #    then a later tool result of the same kind going green.
    seen_failed: set[str] = set()
    for m in msgs:
        name = _tool_name(m)
        if not name:
            continue
        key = name.lower()
        if _is_tool_failed(m, name):
            seen_failed.add(key)
        elif key in seen_failed and _is_tool_green(m, name):
            reasons.append(
                f"error-recovery: a {name} call failed, then a later "
                f"{name} call succeeded")
            seen_failed.discard(key)

    # 3) A course correction by the user.
    for m in user_msgs:
        if _looks_like_correction(_content(m)):
            reasons.append("user-correction: the user corrected the "
                           "agent's course mid-session")
            break

    # 4) A non-obvious route: several different tools used (>= 4 distinct)
    #    hints the solution needed exploration beyond the obvious path.
    distinct = {(_tool_name(m) or "").lower() for m in msgs} - {""}
    if len(distinct) >= 4:
        reasons.append(
            f"non-obvious-route: {len(distinct)} distinct tools used "
            f"to reach the result")

    return reasons
