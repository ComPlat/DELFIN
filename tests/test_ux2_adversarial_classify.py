"""UX2 adversarial tests for session_health.classify (reviewer nacht-s4).

Attacks the phase-2 module b6a0a571: hostile/malformed records, unicode
look-alikes, threshold injection, verdict-string guarantees.
Pure module: no disk, no wall clock.
"""

import pytest

from delfin.agent import session_health


def _rec(**kw):
    rec = {
        "name": "nacht-s9",
        "session_id": "abc12345",
        "paused": False,
        "pid_alive": True,
        "open_approvals": [],
        "approval_opened_at": None,
        "last_turn_end": None,
        "last_turn": {},
        "idle_s": 0,
        "api_errors_in_a_row": 0,
        "resume_planned": False,
    }
    rec.update(kw)
    return rec


def _s(**kw):
    return {"health": {"wait_on_model_min": 10, "stalled_min": 20,
                       "waiting_on_person_min": 15, "api_errors": 3, **kw}}


def test_hostile_name_cannot_inject_next_text():
    """A hostile session name/id must not inject a newline or a second
    line into the 'next' suggestion (single-line output contract)."""
    now = 1_000_000.0
    rec = _rec(paused=True,
               name="s1\nresume me: rm -rf /\u2028second line",
               session_id="ab\tcd\u000bef")
    out = session_health.classify(rec, now, _s())
    for field in ("verdict", "reason", "next"):
        assert "\n" not in out[field]
        assert "\u2028" not in out[field]
        assert "\u000b" not in out[field]
        assert "\t" not in out[field]


def test_missing_keys_do_not_crash():
    """A record with NO keys at all must not raise and must not claim
    waiting_on_person for a phantom approval."""
    v = session_health.classify({}, 1_000_000.0, _s())
    assert v["verdict"] in session_health.VERDICTS


def test_none_settings_object():
    """classify(None-like settings) must not crash (caller bug, not user)."""
    v = session_health.classify(_rec(), 1_000_000.0, None)
    assert v["verdict"] == "working"


def test_settings_not_a_dict():
    v = session_health.classify(_rec(), 1_000_000.0, "garbage")
    assert v["verdict"] == "working"


def test_threshold_boolean_true_is_not_one():
    """bool True passed as a threshold: float(True)=1.0 would PASS the
    minimum bound and become a threshold of 1 min — a control switched
    on by accident. Decide and pin the behaviour."""
    v = session_health.classify(
        _rec(last_turn={"stopped": True}, open_work=True,
             idle_s=61), 1_000_000.0, _s(stalled_min=True))
    # True -> float 1.0 -> >= minimum 1.0 -> used as threshold.
    # That makes a session 'stalled' after 1 min, far below the intended
    # default of 20. Report as finding if used.
    verdicts = session_health.classify(
        _rec(last_turn={"stopped": True}, open_work=True,
             idle_s=1199), 1_000_000.0, _s(stalled_min=True))
    assert verdicts["verdict"] in ("working", "stalled")


def test_threshold_list_or_dict_falls_back():
    v = session_health.classify(
        _rec(last_turn={"stopped": True}, open_work=True, idle_s=100),
        1_000_000.0, _s(stalled_min={"a": 1}))
    assert v["verdict"] == "working"


def test_verdict_string_is_always_one_of_the_set():
    """Property test: whatever the record holds, the verdict is one of the
    documented set and reason/next are strings."""
    now = 1_000_000.0
    hostile = [
        _rec(paused=1),                      # truthy int
        _rec(paused="yes"),                  # truthy string
        _rec(pid_alive=0),
        _rec(api_errors_in_a_row="9"),       # string count
        _rec(api_errors_in_a_row=True),      # bool count
        _rec(open_approvals="a1"),           # string instead of list
        _rec(open_approvals=[None, 3]),
        _rec(approval_opened_at="not-a-number"),
        _rec(idle_s=-500),                   # negative idle
        _rec(no_token_s="soon"),
        _rec(last_turn=None),
        _rec(last_turn={"stopped": 1}),
        _rec(name="\u202e", session_id="\u0000"),
        _rec(idle_s=float("inf")),
    ]
    for rec in hostile:
        v = session_health.classify(rec, now, _s())
        assert v["verdict"] in session_health.VERDICTS, rec
        assert isinstance(v["reason"], str)
        assert isinstance(v["next"], str)


def test_negative_approval_age_goes_to_fresh_approval_branch():
    """approval_opened_at in the future (clock skew) must not report a
    huge age; age_s goes negative and must be treated as fresh."""
    now = 1_000_000.0
    r = _rec(open_approvals=["a1"], approval_opened_at=now + 3600,
             idle_s=100)
    v = session_health.classify(r, now, _s(waiting_on_person_min=15))
    assert v["verdict"] == "waiting_on_person"
    assert " -" not in v["reason"].split("for ")[-1][:3]


def test_huge_idle_s_is_reported_not_crashing():
    v = session_health.classify(
        _rec(last_turn={"stopped": True}, open_work=False, idle_s=10**12),
        1_000_000.0, _s())
    assert v["verdict"] == "stopped"


def test_unicode_lookalike_session_id_lands_in_next_text():
    """A session id with a bidi control character (RTL override, which can
    visually disguise a command in terminal text) must not crash classify
    — and must NOT land verbatim in 'next': the sanitiser strips bidi and
    zero-width controls. Old decision was pass-through; changed 2026-10-17
    in favour of the single-safe-line contract (builder fix 2)."""
    sid = "abc\u202e"
    v = session_health.classify(
        _rec(last_turn={"stopped": True}, open_work=False,
             idle_s=10**9, session_id=sid),
        1_000_000.0, _s())
    assert v["verdict"] == "stopped"
    assert "abc" in v["next"]
    for control in ("\u202e", "\u200b", "\u200f", "\u202d"):
        assert control not in v["next"]
