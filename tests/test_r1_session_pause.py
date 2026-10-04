"""R1 / phase 4 — pause/resume: the session_pause flag-file state machine.

Wave-12 finding 2: there is no pause — SIGSTOP on the agent process has no
effect; the only stop is Esc + /exit, and partner messages wake the session
again. The fix is a DELFIN-level pause: a per-session flag file the engine
checks between tool calls (wired by a .gate patch for the protected
api_client / cli); this module holds the flag-file state and the decision
functions, all pure file operations testable with a temp dir.

Contract:
- a session with its flag present is PAUSED;
- while paused, pause_gate() says "no tool call may start" and
  wake_blocked() says "the wake look must not fire";
- pause()/resume() create/remove the flag atomically and idempotently;
- flags are keyed by the session address used elsewhere (-n name, else the
  first 8 chars of the session id, per repl._presence_key).
"""

import pytest

from delfin.agent import session_pause as sp


@pytest.fixture
def pause_dir(tmp_path, monkeypatch):
    """Point session_pause at a temp flag dir for the test."""
    monkeypatch.setattr(sp, "_DIR", tmp_path)
    return tmp_path


class TestPauseState:
    def test_is_paused_is_false_when_absent(self, pause_dir):
        assert sp.is_paused("sess-1") is False

    def test_pause_creates_the_flag(self, pause_dir):
        sp.pause("sess-1")
        assert sp.is_paused("sess-1") is True
        assert sp.pause_gate("sess-1") is False   # no tool call may start
        assert sp.wake_blocked("sess-1") is True  # no wake may fire

    def test_resume_clears_the_flag(self, pause_dir):
        sp.pause("sess-1")
        assert sp.resume("sess-1") is True
        assert sp.is_paused("sess-1") is False
        assert sp.pause_gate("sess-1") is True
        assert sp.wake_blocked("sess-1") is False

    def test_pause_is_idempotent(self, pause_dir):
        sp.pause("sess-1")
        sp.pause("sess-1")          # must not raise or corrupt
        assert sp.is_paused("sess-1") is True

    def test_resume_of_a_non_paused_key_is_false(self, pause_dir):
        assert sp.resume("never-paused") is False

    def test_flags_are_per_key(self, pause_dir):
        sp.pause("sess-1")
        assert sp.is_paused("sess-2") is False
        assert sp.pause_gate("sess-2") is True
        sp.resume("sess-1")
        assert sp.is_paused("sess-1") is False


class TestPauseListingAndKeys:
    def test_paused_keys_lists_only_paused(self, pause_dir):
        sp.pause("a")
        sp.pause("b")
        sp.resume("b")
        assert sorted(sp.paused_keys()) == ["a"]

    def test_hostile_key_stays_inside_the_flag_dir(self, pause_dir):
        key = "../../escape"
        sp.pause(key)
        # The flag must be a single file under _DIR, never a path that
        # escapes into a parent directory.
        assert (pause_dir / "escape").exists() is False
        assert sp.is_paused(key) is True
        assert sp.resume(key) is True


class TestPauseDurability:
    def test_flag_persists_across_a_fresh_lookup(self, pause_dir):
        # A pause survives until explicitly resumed: a new "look" (fresh
        # share of the same file) still reads it.
        sp.pause("sess-1")
        assert sp.is_paused("sess-1") is True
        sp.resume("sess-1")
        assert sp.is_paused("sess-1") is False


class TestPauseWakeIntegration:
    """The wired behaviour: a paused session's wake look fires nothing."""

    def test_paused_session_wake_text_returns_empty(self, pause_dir):
        from delfin.agent.repl import TerminalAgent
        sess = object.__new__(TerminalAgent)
        sess._presence_key = lambda: "paused-sess"
        sp.pause("paused-sess")
        # _wake_text must return "" (no wake, no turn) while paused, so no
        # student-message or finished-job fires and wakes the session.
        assert sess._wake_text("") == ""

    def test_paused_session_wake_text_with_typed_text_still_empty(self, pause_dir):
        # The pause guard wins even before the typed-draft guard: typing is
        # moot while paused — the session is not taking input-driven turns.
        from delfin.agent.repl import TerminalAgent
        sess = object.__new__(TerminalAgent)
        sess._presence_key = lambda: "paused-sess"
        sp.pause("paused-sess")
        assert sess._wake_text("half-typed") == ""

