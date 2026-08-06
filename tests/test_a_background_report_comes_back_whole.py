"""A backgrounded report came back as a stub, under an id that never existed.

**The report was cut to 1.6 kB.** Stored session messages are trimmed to
head 1200 + tail 400 characters, which is right for the resume recap it
was written for -- a resumed subagent needs the gist, not the transcript.
But the session store is the ONLY route by which a BACKGROUND report
reaches its parent: ``subagent_result`` reads the last assistant message
out of it. So the longer the delegate's report, the more of it was
thrown away, and the parent got a stub with a trim marker in the middle
where the findings used to be. A foreground run of the same work
returned the whole thing.

**The id could not be collected.** Backgrounding reserves a fresh id up
front and tells the parent to collect with it. ``run_subagent`` then
picks ``resume_from`` over the supplied id whenever a prior session was
loaded, so a backgrounded RESUME stored its result under the resume id
and handed back the reserved one. The parent polled an id that never
existed and was told "unknown" forever, for work that had completed.

The report is now stored in full, separately from the trimmed
conversation, and the id handed back is the one the result will be
stored under.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import subagents as sa


@pytest.fixture
def sessions(tmp_path, monkeypatch):
    d = tmp_path / "subagent_sessions"
    d.mkdir()
    monkeypatch.setattr(sa, "_SESSIONS_DIR", d)
    return d


_LONG = "finding %03d: the parser drops the trailing field\n" * 200


def _save(sa_id: str, report: str) -> None:
    sa._save_subagent_session(
        sa_id,
        subagent_type="explore",
        description="d",
        messages=[{"role": "user", "content": "go"},
                  {"role": "assistant", "content": report}],
        interactions=[],
    )


# ---------------------------------------------------------------------------
# The whole report survives the trip
# ---------------------------------------------------------------------------

def test_a_long_report_is_not_cut_to_a_stub(sessions):
    _save("bg1", _LONG)
    out = sa.get_subagent_result("bg1")
    assert out["final_text"] == _LONG


def test_the_trim_marker_does_not_appear_in_the_report(sessions):
    _save("bg1", _LONG)
    assert "trimmed for subagent session store" not in \
        sa.get_subagent_result("bg1")["final_text"]


def test_a_short_report_is_unchanged(sessions):
    _save("bg2", "done, nothing to report")
    assert sa.get_subagent_result("bg2")["final_text"] == \
        "done, nothing to report"


def test_the_conversation_is_still_trimmed(sessions):
    """The recap a resumed subagent replays is what the trim is for, and
    it must stay bounded."""
    _save("bg3", _LONG)
    rec = json.loads((sessions / "bg3.json").read_text(encoding="utf-8"))
    assert len(rec["messages"][-1]["content"]) < len(_LONG)


def test_an_enormous_report_is_still_bounded(sessions):
    huge = "a line of findings\n" * (sa._MAX_STORED_REPORT_CHARS // 19 + 400)
    _save("bg4", huge)
    got = sa.get_subagent_result("bg4")["final_text"]
    assert len(got) <= sa._MAX_STORED_REPORT_CHARS + 200
    assert "truncated" in got.lower(), "a cut report has to say it was cut"


def test_a_report_that_is_one_long_line_does_not_stall_the_parent(sessions):
    """The shared scanners are quadratic in the length of a single LINE.
    Uncapped, a 50 kB blob cost 47 seconds of the parent's turn."""
    import time
    _save("blob", "x" * 120_000)
    start = time.time()
    sa.get_subagent_result("blob")
    assert time.time() - start < 15.0


def test_a_session_without_the_field_still_reports(sessions):
    """Written before the field existed: fall back to the conversation."""
    _save("legacy", "the old way")
    rec = json.loads((sessions / "legacy.json").read_text(encoding="utf-8"))
    rec.pop("final_report", None)
    (sessions / "legacy.json").write_text(json.dumps(rec), encoding="utf-8")
    assert sa.get_subagent_result("legacy")["final_text"] == "the old way"


# ---------------------------------------------------------------------------
# The id handed back is the id that works
# ---------------------------------------------------------------------------

def test_a_resumed_run_keeps_the_resume_id():
    """The rule itself is right: a resumed session accumulates under one id."""
    assert sa.collectable_sa_id(reserved="new123", resume_from="old456") == \
        "old456"


def test_a_fresh_run_keeps_the_reserved_id():
    assert sa.collectable_sa_id(reserved="new123", resume_from="") == "new123"


def test_no_resume_and_no_reservation_is_empty():
    assert sa.collectable_sa_id(reserved="", resume_from="") == ""


def test_the_background_dispatch_reports_the_collectable_id():
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    i = src.index('"status": "started_in_background"')
    assert "collectable_sa_id" in src[i - 3000:i + 500], (
        "the reserved id is still handed back unconditionally")
