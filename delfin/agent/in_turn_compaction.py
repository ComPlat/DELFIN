"""Deterministic compaction BETWEEN two rounds of one tool-call turn.

Pure module: takes the round history, returns the compacted history and
a record. The wiring (api_client's round loop) decides when to call it;
nothing here reads settings, makes model calls or touches disk.

Contract (fixed with the operator 2026-09-26):

    compact(messages, budget_chars, *, keep_last_rounds=2,
            turn_rows_base=0) -> (new_messages, record)

* ``record`` is ``None`` when nothing needed compacting (under budget).
* NEVER tear open a tool_call/tool_result pair: a round is compacted as
  a whole — both its assistant call message and its tool result go, or
  both stay.
* User goals and everything before ``turn_rows_base`` (the session
  prefix shared by every request — the prefix cache) stay verbatim.
* Compacted rounds are replaced by ONE summary message built
  deterministically from the working-state extractors (open facts:
  denials, test outcomes, touched files) plus a digest line per
  compacted tool result, so the model keeps the WHAT of earlier rounds.
"""

from __future__ import annotations

from .compaction_log import compare_loss

_SUMMARY_HEADER = "[Conversation summary — earlier rounds of this turn compacted]"

# The tool-result digest is capped so a pathological result cannot blow
# the summary; the round's identity (tool + call id) always survives.
_DIGEST_CHARS = 160


def _is_goal(msg: dict) -> bool:
    """A real user goal: user role, not a synthetic machine turn."""
    if msg.get("role") != "user":
        return False
    content = msg.get("content", "")
    if not isinstance(content, str):
        return False
    return not content.lstrip().startswith(
        ("[Command results]", "[Verify]", "[System]",
         "[Conversation summary", "[Working state"))


def _tool_chars(messages: list[dict]) -> int:
    return sum(len(str(m.get("content", "")))
               for m in messages if m.get("role") == "tool")


def _rounds(messages: list[dict], start: int) -> list[tuple[int, int]]:
    """(call_idx, result_idx) pairs of complete tool rounds at or after
    ``start``. A call whose result is missing is not a round — it stays
    untouched rather than being torn open."""
    by_call: dict[str, int] = {}
    for i in range(start, len(messages)):
        m = messages[i]
        if m.get("role") == "assistant":
            for tc in m.get("tool_calls") or []:
                cid = str(tc.get("id") or "")
                if cid:
                    by_call[cid] = i
    rounds: list[tuple[int, int]] = []
    for i in range(start, len(messages)):
        m = messages[i]
        if m.get("role") == "tool":
            j = by_call.get(str(m.get("tool_call_id") or ""))
            if j is not None:
                rounds.append((j, i))
    return rounds


def _digest_line(msg: dict) -> str:
    content = str(msg.get("content", ""))
    head = " ".join(content.split())[:_DIGEST_CHARS]
    return f"  round result {msg.get('tool_call_id', '?')}: {head}"


def compact(
    messages: list[dict],
    budget_chars: int,
    *,
    keep_last_rounds: int = 2,
    turn_rows_base: int = 0,
) -> tuple[list[dict], dict | None]:
    """Compact old rounds of one turn; see the module docstring."""
    messages = list(messages or [])
    start = max(0, int(turn_rows_base or 0))
    rounds = _rounds(messages, start)
    if len(rounds) <= max(0, int(keep_last_rounds)):
        # Nothing older than the kept tail — a no-op, not a shred.
        return messages, None
    if _tool_chars(messages) <= int(budget_chars):
        return messages, None

    # ONE deep cut: decide up front how many of the OLDEST rounds must
    # go to bring the remaining tool output under the budget (the same
    # one-cut-per-crossing rule as the elision's _ELIDE_TO_FRACTION —
    # a cut rewrites old messages, and piecemeal cuts break the prefix
    # cache for nothing). ``keep_last_rounds`` is a preference, not a
    # suicide pact: when even dropping every unprotected round cannot
    # reach the budget, protected rounds go too — OLDEST first — because
    # the alternative is the turn dying of context_length_exceeded (the
    # measured 2026-09-26 gap). The LAST round's open tool-call pair is
    # never touched: the protocol requires its answer.
    protected_pairs = rounds[-int(keep_last_rounds):] if keep_last_rounds else []
    last_pair = rounds[-1] if rounds else None
    chars = _tool_chars(messages)
    doomed: list[tuple[int, int]] = []
    for pair in rounds:
        if chars <= int(budget_chars):
            break
        if pair == last_pair:
            break
        doomed.append(pair)
        chars -= len(str(messages[pair[1]].get("content", "")))
    del protected_pairs

    if not doomed:
        return messages, None

    drop_idx: set[int] = set()
    for call_idx, result_idx in doomed:
        drop_idx.add(call_idx)
        drop_idx.add(result_idx)
    replaced = [messages[i] for i in sorted(drop_idx)]
    kept = [m for i, m in enumerate(messages) if i not in drop_idx]

    # Deterministic summary: working-state facts + one digest per round.
    texts = [str(m.get("content", "")) for m in replaced]
    state_block = ""
    try:
        from .working_state import build_working_state_block
        state_block = build_working_state_block(replaced, session_id="")
    except Exception:
        state_block = ""
    lines = [f"{_SUMMARY_HEADER}"]
    if state_block:
        lines.append(state_block.rstrip())
    lines.append("What the compacted rounds did (digests):")
    for m in replaced:
        if m.get("role") == "tool":
            lines.append(_digest_line(m))
    summary_msg = {"role": "user", "content": "\n".join(lines)}

    # Splice the summary in at the position of the first removed round
    # so the kept prefix and the kept rounds keep their relative order.
    first_drop = min(drop_idx)
    pos = sum(1 for i, _m in enumerate(messages)
              if i < first_drop and i not in drop_idx)
    new_messages = kept[:pos] + [summary_msg] + kept[pos:]

    # The determinism contract: same input, same output. The record
    # carries stable numbers; compare_loss names what was NOT carried.
    record = {
        "kind": "in_turn",
        "messages_compacted": len(replaced),
        "rounds_compacted": len(doomed),
        "budget_chars": int(budget_chars),
        "tool_chars_after": _tool_chars(new_messages),
    }
    try:
        record["lost"] = compare_loss(
            texts, state_block=state_block, summary_text=summary_msg["content"])
    except Exception:
        pass
    return new_messages, record
