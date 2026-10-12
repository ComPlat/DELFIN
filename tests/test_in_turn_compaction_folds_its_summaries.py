"""Each in-turn compaction folds the turn's earlier summary into its own.

Input: a turn grown one round at a time, ``compact`` called after every
round as the round loop calls it. Output asserted: at most one
``SummaryRow`` per history, bounded by ``_SUMMARY_MAX_CHARS``, newest
digest first, the oldest counted when dropped; rows that only look like
a summary (a user steer, an assistant text, an earlier turn's summary)
are never folded.
"""

from __future__ import annotations

import json

from delfin.agent import api_client
from delfin.agent import in_turn_compaction as itc
from delfin.agent.in_turn_compaction import SummaryRow, compact, is_summary_row

H = itc._SUMMARY_HEADER


def _round(i, chars=4000, parallel=1):
    ids = [f"call-{i}" + (f"-{k}" if parallel > 1 else "")
           for k in range(parallel)]
    rows = [{"role": "assistant", "content": f"round {i}",
             "tool_calls": [{"id": cid, "type": "function",
                             "function": {"name": "bash",
                                          "arguments": "{}"}}
                            for cid in ids]}]
    for cid in ids:
        rows.append({"role": "tool", "tool_call_id": cid,
                     "content": f"result of {cid} " + "x" * chars})
    return rows


def _grow(n, prefix=None, base=None, between=None, parallel=1):
    msgs = list(prefix or [{"role": "user", "content": "goal"}])
    base = len(msgs) if base is None else base
    for i in range(n):
        msgs.extend(_round(i, parallel=parallel))
        if between is not None:
            msgs.extend(between(i))
        msgs, _rec = compact(msgs, 9000, turn_rows_base=base)
    return msgs


def _own(msgs):
    return [m for m in msgs if is_summary_row(m)]


def _orphans(msgs):
    calls, answered, bad = set(), set(), 0
    for m in msgs:
        for tc in m.get("tool_calls") or []:
            calls.add(tc["id"])
        if m.get("role") == "tool":
            if m["tool_call_id"] not in calls:
                bad += 1
            answered.add(m["tool_call_id"])
    return bad + len(calls - answered)


def test_forty_compactions_leave_one_summary_row():
    msgs = _grow(40)
    own = _own(msgs)
    assert len(own) == 1
    assert sum(str(m.get("content", "")).startswith(H) for m in msgs) == 1
    text = own[0]["content"]
    # Newest first: the newest compacted round precedes the oldest.
    assert text.index("call-37:") < text.index("call-0:")
    assert _orphans(msgs) == 0
    # The summary stands where the dropped rounds were: after the goal,
    # before every kept round.
    assert msgs[0]["content"] == "goal" and msgs[1] is own[0]


def test_parallel_calls_fold_without_orphans():
    msgs = _grow(25, parallel=3)
    assert len(_own(msgs)) == 1
    assert _orphans(msgs) == 0


def test_hundreds_of_rounds_stay_under_the_cap_newest_kept():
    msgs = _grow(400)
    text = _own(msgs)[0]["content"]
    assert len(text) <= itc._SUMMARY_MAX_CHARS
    kept = {m["tool_call_id"] for m in msgs if m.get("role") == "tool"}
    newest = max(i for i in range(400) if f"call-{i}" not in kept)
    assert f"round result call-{newest}:" in text
    assert "call-0:" not in text
    listed = text.count(itc._DIGEST_PREFIX)
    unlisted = int(itc._UNLISTED_RE.match(text.splitlines()[-1]).group(1))
    assert listed + unlisted == newest + 1, "a dropped round went uncounted"


def test_a_steer_that_looks_like_a_summary_is_not_swallowed():
    steer = H + "\nplease also check the docs folder"
    msgs = _grow(12, between=lambda i: (
        [{"role": "user", "content": steer}] if i == 3 else []))
    assert sum(m.get("content") == steer for m in msgs) == 1
    assert len(_own(msgs)) == 1


def test_an_assistant_text_that_looks_like_a_summary_is_not_swallowed():
    fake = {"role": "assistant", "content": H + "\nmy own notes"}
    msgs = _grow(12, between=lambda i: [dict(fake)] if i == 2 else [])
    assert sum(m.get("content") == fake["content"] for m in msgs) == 1


def test_an_earlier_turns_summary_stays_verbatim():
    old_plain = {"role": "user", "content": H + "\nfrom turn 1"}
    old_typed = SummaryRow(role="user", content=H + "\nalso turn 1")
    prefix = [{"role": "user", "content": "goal 1"}, old_plain, old_typed,
              {"role": "assistant", "content": "ok"},
              {"role": "user", "content": "goal 2"}]
    msgs = _grow(10, prefix=prefix)
    assert msgs[:5] == prefix and msgs[1] is old_plain and msgs[2] is old_typed
    assert len([m for m in msgs[5:] if is_summary_row(m)]) == 1


def test_the_working_state_of_a_folded_summary_is_carried():
    denial = "bash: permission denied for rm -rf build (do not retry)"
    first = SummaryRow(role="user", content=(
        H + "\n[Working state — session recap, rebuilt at each compaction]\n"
        "Recent denials (do not retry):\n  " + denial + "\n"
        + itc._DIGESTS_HEADING + "\n  round result call-old: x"))
    msgs = _grow(6, prefix=[{"role": "user", "content": "goal"}, first],
                 base=1)
    text = _own(msgs)[0]["content"]
    assert denial in text
    assert "round result call-old: x" in text
    assert all(m is not first for m in msgs) and len(_own(msgs)) == 1


def test_a_pending_note_is_not_folded_and_lands_once():
    msgs = _grow(4)
    note = api_client._PENDING_JOBS_MARK + ": slurm job 4711."
    api_client._append_pending_jobs_to_summary(msgs, note)
    msgs.extend(_round(99))
    msgs, rec = compact(msgs, 9000, turn_rows_base=1)
    assert rec and rec["summaries_folded"] == 1
    assert "4711" not in _own(msgs)[0]["content"]
    api_client._append_pending_jobs_to_summary(msgs, note)
    assert json.dumps(msgs).count("4711") == 1


def test_a_steer_after_the_summary_does_not_get_the_note():
    msgs = _grow(4)
    steer = {"role": "user", "content": H + " (pasted)"}
    msgs.append(steer)
    api_client._append_pending_jobs_to_summary(
        msgs, api_client._PENDING_JOBS_MARK + ": slurm job 7.")
    assert steer["content"] == H + " (pasted)"
    assert "slurm job 7" in _own(msgs)[0]["content"]


def test_the_summary_row_is_a_plain_row_on_the_wire():
    row = _own(_grow(4))[0]
    assert set(row) == {"role", "content"}
    assert json.loads(json.dumps(row)) == dict(row)
    from delfin.agent.turn_history import compact_turn_rows
    out = compact_turn_rows([row])
    assert type(out[0]) is dict and out[0]["content"].startswith(H)
