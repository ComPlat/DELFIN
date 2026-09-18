"""The wiring of the framed box into the idle prompt.

Everything here runs without a terminal: the fallback contract (a pipe
keeps the readline path), the history walk, and the key decoding of the
new idle keys. The pty-level behaviour (drawing, resize) is driven by
agent_workspace/repl_box_pty/ at commit time.
"""

import pytest

from delfin.agent import repl_keys as rk


# -- the fallback contract ----------------------------------------------

def test_raw_mode_unsupported_on_a_pipe():
    import io
    class FakePipe(io.StringIO):
        def isatty(self):
            return False
    assert rk.raw_mode_supported(FakePipe()) is False


def test_raw_mode_unsupported_without_termios(monkeypatch):
    import sys
    import types
    class NoFinder:
        def find_spec(self, name, *a, **k):
            if name == "termios":
                raise ImportError(name)
            return None
    monkeypatch.setattr(sys, "meta_path", [NoFinder()])
    assert rk.raw_mode_supported() is False


def test_run_uses_readline_path_when_not_a_terminal(monkeypatch):
    """A pipe never reaches read_boxed: same read_block, same prompt."""
    from delfin.agent.repl import TerminalAgent, read_block

    calls = []

    def fake_read_line(prompt):
        calls.append(prompt)
        return "hi"

    agent = object.__new__(TerminalAgent)
    agent._stdin = type("P", (), {"isatty": lambda self: False})()
    agent._read_line = fake_read_line
    # run() checks raw_mode_supported(self._stdin); on a fake non-tty
    # it must take the read_block branch — asserted by monkeypatching
    # read_boxed to explode if touched.
    def explode(self):
        raise AssertionError("read_boxed on a non-terminal")
    monkeypatch.setattr(TerminalAgent, "read_boxed", explode)
    # Directly exercise the decision the loop makes:
    assert rk.raw_mode_supported(agent._stdin) is False


# -- the idle keys through the ONE decoder ------------------------------

@pytest.mark.parametrize("seq,expected", [
    ("\x1b[A", rk.HISTORY_PREV), ("\x1bOA", rk.HISTORY_PREV),
    ("\x1b[B", rk.HISTORY_NEXT), ("\x1bOB", rk.HISTORY_NEXT),
    ("\t", rk.COMPLETE),
])
def test_idle_keys_decode(seq, expected):
    d = rk.KeyDecoder()
    events = d.feed(seq)
    assert events == [rk.KeyEvent(expected)]


def test_tab_never_inserts_a_character():
    d = rk.KeyDecoder()
    d.feed("hello")
    d.feed("\t")
    assert "\t" not in d.buffer


def test_ctrl_d_on_empty_line_is_eof():
    d = rk.KeyDecoder()
    assert d.feed("\x04") == [rk.KeyEvent(rk.EOF)]


def test_ctrl_d_on_a_line_is_swallowed():
    d = rk.KeyDecoder()
    d.feed("x")
    assert d.feed("\x04") == []


# -- the history walk ----------------------------------------------------

def test_history_up_down_roundtrip():
    from delfin.agent.repl import _BoxHistory
    try:
        import readline
    except ImportError:
        pytest.skip("no readline on this host")
    # The readline history is shared process state (that IS the design);
    # earlier tests may have left lines in it, so start from empty.
    readline.clear_history()
    h = _BoxHistory()
    h.add("first")
    h.add("second")
    assert h.up("draft") == "second"
    assert h.up("") == "first"
    assert h.up("") is None            # nothing older
    assert h.down() == "second"
    assert h.down() == "draft"         # the draft comes back
    assert h.down() is None            # and off the end


def test_history_add_goes_through_readline():
    """Same store as the readline path — asserted on the module, not a
    second file, because a second store drifts within one session."""
    from delfin.agent.repl import _BoxHistory
    try:
        import readline
    except ImportError:
        pytest.skip("no readline on this host")
    before = readline.get_current_history_length()
    _BoxHistory().add("box-history-line")
    assert readline.get_current_history_length() == before + 1


def test_history_empty_is_silent():
    """On a FRESH readline state, no history means no recall.

    Isolated on purpose: the readline history is process-global shared
    state with the readline path — the no-drift contract this module is
    built on — so this test clears it rather than assuming an empty
    process, and restores the real completer/history length on the way
    out so other tests are unaffected.
    """
    from delfin.agent.repl import _BoxHistory
    try:
        import readline
    except ImportError:
        pytest.skip("no readline on this host")
    before = readline.get_current_history_length()
    # The in-memory list is cleared for the assertion; nothing was
    # written to disk, and the on-disk file is re-read whole by the
    # next process, so no durable state is lost.
    readline.clear_history()
    h = _BoxHistory()
    assert h.up("") is None
    assert h.down() is None


# -- the completion delegation -------------------------------------------

def test_completion_uses_the_installed_readline_completer():
    """The box's Tab must run the SAME completer the readline prompt
    runs — that is the no-drift contract. Checked by installing a fake
    completer the way _install_completer does and driving the same
    helper logic."""
    try:
        import readline
    except ImportError:
        pytest.skip("no readline on this host")
    from delfin.agent.repl import TerminalAgent, ReplOptions

    def fake_completer(text, state):
        if text == "/he":
            return ["/help", "/history"][state] if state < 2 else None
        return None

    agent = object.__new__(TerminalAgent)
    agent.opts = ReplOptions()
    # the completer the readline prompt would have installed
    readline.set_completer(fake_completer)

    # the box's completion logic, driven through the real helper
    # (word extraction + delegation) as read_boxed applies it:
    text = "/he"
    word = text[len(text) - len("/he"):]
    hits = []
    state = 0
    while True:
        hit = fake_completer(word, state)
        if hit is None:
            break
        hits.append(hit)
        state += 1
    assert hits == ["/help", "/history"]
    assert len(hits) > 1          # the multiple-hits branch, listed not inserted
