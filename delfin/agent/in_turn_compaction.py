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
* A request holds at most one in-turn summary row. Each compaction
  folds the summaries earlier compactions of this turn left (rows at or
  after ``turn_rows_base`` that are ``SummaryRow`` instances) into the
  new one: their working-state block is re-read by the extractors
  together with the dropped rounds, their digest lines follow the new
  ones (newest first), and the whole row is capped at
  ``_SUMMARY_MAX_CHARS`` by dropping the OLDEST digest lines, which are
  then counted on one line. The row stands where the earliest folded or
  dropped row stood.

  Measured (tests/test_zz_u11_one_summary_per_turn.py, 32 rounds): on
  6d6295af the last request held 30 summary rows (8510 chars) on the
  tool route and 15 (6965 chars) on the whole-request route, growing by
  one row per compaction without bound.

  Recognition is by type, not by text: a user steer or an assistant
  message that starts with the header is an ordinary row and is never
  folded. A marker key on the row was rejected because the rows go to
  the endpoint unchanged and a strict endpoint refuses unknown keys;
  the header prefix alone was rejected because a pasted header would be
  swallowed. A summary carried in from an earlier turn's history is a
  plain dict before ``turn_rows_base`` and stays verbatim (prefix
  cache).
"""

from __future__ import annotations

import re

from .compaction_log import compare_loss

_SUMMARY_HEADER = "[Conversation summary — earlier rounds of this turn compacted]"

# The tool-result digest is capped so a pathological result cannot blow
# the summary; the round's identity (tool + call id) always survives.
_DIGEST_CHARS = 160

# Ceiling on the one in-turn summary row (~2000 tokens at 4 chars per
# token, 1.5 % of a 131k-token window): header, a working-state block of
# at most working_state._MAX_BLOCK_CHARS (2200) and ~29 digest lines of
# <= ~200 chars. turn_history keeps the head 4000 chars of a user row in
# later turns, which is why the digests are listed newest first.
_SUMMARY_MAX_CHARS = 8000

_DIGESTS_HEADING = "What the compacted rounds did (digests, newest first):"
_DIGEST_PREFIX = "  round result "
_UNLISTED_RE = re.compile(
    r"^  \.\.\. and (\d+) older round\(s\) of this turn, not itemised$")


class SummaryRow(dict):
    """The in-turn summary row ``compact`` creates.

    A plain ``dict`` to every consumer (equality, JSON, the SDK's
    transform, ``turn_history``); the type is what identifies it to the
    next compaction of the same turn, so text that merely looks like a
    summary is never folded. ``copy.deepcopy`` and in-place content edits
    (the pending-jobs note) keep the type; ``dict(row)`` does not, and a
    row that lost it is left in place rather than folded."""


def is_summary_row(msg) -> bool:
    """True for a summary row ``compact`` built in this request."""
    return isinstance(msg, SummaryRow) and msg.get("role") == "user"


def _fold_parts(content: str) -> tuple[str, list[str], int]:
    """(working-state block, digest lines, unlisted-round count) of a
    summary built by ``compact``; anything else in it (the pending-jobs
    note, which the caller re-appends) is dropped."""
    lines = str(content or "").splitlines()
    state: list[str] = []
    digests: list[str] = []
    unlisted = 0
    phase = "state"                  # state -> digests -> tail
    for line in lines[1:]:
        if phase == "state" and line.startswith(
                "What the compacted rounds did"):
            phase = "digests"
        elif phase == "state":
            if state or line.startswith("[Working state"):
                state.append(line)
        elif phase == "digests" and line.startswith(_DIGEST_PREFIX):
            digests.append(line)
        elif phase == "digests" and _UNLISTED_RE.match(line):
            unlisted += int(_UNLISTED_RE.match(line).group(1))
        else:
            phase = "tail"           # the pending-jobs note, or anything after
    return "\n".join(state).strip(), digests, unlisted


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
    # Rounds are cut per assistant message: one message with parallel
    # tool calls is several pairs sharing its index, and dropping it while
    # one of its results stays leaves an orphaned tool row (API 400). The
    # group holding the last pair is never cut.
    groups: dict[int, list[tuple[int, int]]] = {}
    for pair in rounds:
        groups.setdefault(pair[0], []).append(pair)
    doomed: list[tuple[int, int]] = []
    for call_idx in sorted(groups):
        if chars <= int(budget_chars):
            break
        if last_pair is not None and call_idx == last_pair[0]:
            break
        for pair in groups[call_idx]:
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
    # Earlier summaries of this turn are folded into the new one.
    folded_idx = [i for i in range(start, len(messages))
                  if i not in drop_idx and is_summary_row(messages[i])]
    folded = [messages[i] for i in folded_idx]
    gone = drop_idx | set(folded_idx)
    kept = [m for i, m in enumerate(messages) if i not in gone]

    # Deterministic summary: working-state facts + one digest per round,
    # newest first, then the digests the folded summaries carried.
    texts = [str(m.get("content", "")) for m in folded + replaced]
    old_states: list[dict] = []
    old_digests: list[str] = []
    unlisted = 0
    for m in reversed(folded):                     # newest summary first
        st, dg, un = _fold_parts(str(m.get("content", "")))
        if st:
            old_states.insert(0, {"role": "user", "content": st})
        old_digests.extend(dg)
        unlisted += un
    state_block = ""
    try:
        from .working_state import build_working_state_block
        # A prior block is read as a machine turn ("[Working state"),
        # oldest first, so the extractors' newest-first rules hold.
        state_block = build_working_state_block(
            old_states + replaced, session_id="")
    except Exception:
        state_block = ""
    head = [f"{_SUMMARY_HEADER}"]
    if state_block:
        head.append(state_block.rstrip())
    head.append(_DIGESTS_HEADING)
    digests = [_digest_line(m) for m in reversed(replaced)
               if m.get("role") == "tool"] + old_digests
    room = _SUMMARY_MAX_CHARS - sum(len(x) + 1 for x in head) - 80
    listed: list[str] = []
    for line in digests:
        if room - (len(line) + 1) < 0:
            break
        listed.append(line)
        room -= len(line) + 1
    unlisted += len(digests) - len(listed)
    lines = head + listed
    if unlisted:
        lines.append(f"  ... and {unlisted} older round(s) of this turn, "
                     "not itemised")
    text = "\n".join(lines)
    if len(text) > _SUMMARY_MAX_CHARS:      # unreachable with the caps above
        text = text[:_SUMMARY_MAX_CHARS]
    summary_msg = SummaryRow(role="user", content=text)

    # Splice the summary in at the position of the earliest removed row
    # (a folded summary or a dropped round) so the kept prefix and the
    # kept rounds keep their relative order.
    first_drop = min(gone)
    pos = sum(1 for i, _m in enumerate(messages)
              if i < first_drop and i not in gone)
    new_messages = kept[:pos] + [summary_msg] + kept[pos:]

    # The determinism contract: same input, same output. The record
    # carries stable numbers; compare_loss names what was NOT carried.
    record = {
        "kind": "in_turn",
        "messages_compacted": len(replaced),
        "rounds_compacted": len(doomed),
        "summaries_folded": len(folded),
        "summary_chars": len(text),
        "budget_chars": int(budget_chars),
        "tool_chars_after": _tool_chars(new_messages),
    }
    try:
        record["lost"] = compare_loss(
            texts, state_block=state_block, summary_text=summary_msg["content"])
    except Exception:
        pass
    return new_messages, record
