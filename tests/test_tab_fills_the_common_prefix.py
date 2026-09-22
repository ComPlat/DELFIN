"""Tab at the box: fill the longest common prefix first, like readline.

Measured 2026-09-20: ``/mod`` + Tab printed ``/mode`` and ``/model`` on
the transcript and left the box at ``/mod``. The docstring of
``_complete_word`` claims "exactly as readline does" — readline inserts
the longest common prefix FIRST and lists only when that adds nothing.
And ``/`` alone opened nothing at all, although the completer offers
every command for an empty word — with 40 commands the one thing a
lonely slash is likely to mean is "show me the commands".

Both faults are in the multiple-hits branch of ``_complete_word``: it
listed and returned the text unchanged. The fix follows readline: when
the hits share a prefix longer than what was typed, insert it and still
list; when the prefix IS the typed word (or a unique hit), today's
behaviour stands.

Driven through ``read_boxed`` against the model terminal in ``vt.py``
(a screen WITH a right margin — see its module docstring), because the
question "what is on screen after Tab" is answered by applying the
escapes, not by arguing about the code.
"""

from __future__ import annotations

import io

import pytest

from delfin.agent import repl as R
from delfin.agent import repl_keys as rk
from delfin.agent import repl_render as rr

from vt import Screen


class _Keys:
    """A RawMode stand-in that hands the loop a scripted keyboard."""

    def __init__(self, chunks):
        self.chunks = list(chunks)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    @property
    def active(self):
        return True

    def read_ready(self, timeout=0.0):
        return self.chunks.pop(0) if self.chunks else "\r"

    def restore(self):
        pass


@pytest.fixture()
def boxed(monkeypatch):
    """read_boxed against a scripted keyboard, the REAL completer, and
    a model screen with a right margin. Restores the process-global
    readline completer on the way out."""
    def _run(keys, *, width=80, transcript=()):
        try:
            import readline
        except ImportError:
            pytest.skip("no readline on this host")
        agent, screen = _agent(monkeypatch, width, transcript)
        monkeypatch.setattr(rk, "RawMode", lambda *a, **k: _Keys(keys))
        agent._install_completer()
        text = agent.read_boxed()
        return text, screen
    return _run


def _agent(monkeypatch, width, transcript):
    screen = Screen(width, transcript)
    agent = R.TerminalAgent.__new__(R.TerminalAgent)
    agent.err = screen
    agent._flush_err = screen.flush
    agent._stdin = io.StringIO()
    agent._width_dirty = False
    agent._cycle_mode = lambda: None

    class _T:
        pass
    agent.transcript = _T()
    agent.transcript.width = width
    agent.transcript.refresh_width = lambda: None
    agent.transcript.theme = rr.theme_for(io.StringIO(), mode="never")
    agent.transcript.chrome = lambda text: screen.write(text + "\r\n")
    return agent, screen


def _line_with(screen, fragment):
    return [r for r in screen.text() if fragment in r]


def test_tab_fills_the_common_prefix_and_lists_the_rest(boxed):
    """``/mod`` + Tab: the box holds ``/mode`` (the longest common
    prefix of /mode and /model) and both hits are listed below."""
    text, screen = boxed(["/mod", "\t", "\r"])
    assert text == "/mode", screen.text()
    assert _line_with(screen, "/model"), "the other hit is listed"
    assert len(_line_with(screen, "/mode")) >= 1


def test_a_unique_hit_completes_in_place_as_before(boxed):
    text, screen = boxed(["/cle", "\t", "\r"])
    assert text == "/clear"


def test_a_bare_slash_lists_the_commands(boxed):
    """``/`` + Tab with nothing else typed: the commands are on screen.

    Note: the 2026-09-20 measurement said '/' opened nothing. On this
    stand (42ad59b0) read_boxed DOES list — the first 20 of the 41
    builtins, alphabetically. Either the measurement caught a different
    path (the readline prompt), or the box changed since. This test
    pins what the code does now, so it cannot silently regress."""
    text, screen = boxed(["/", "\t", "\r"])
    assert text == "/"                      # nothing false was inserted
    listed = screen.text()
    assert any("/agents" in r for r in listed), listed[:6]
    assert any("/approve" in r for r in listed), listed[:6]


def test_the_common_prefix_of_one_hit_is_the_hit(boxed):
    """When the prefix fills everything, listing again would say what
    the box already shows — readline does not, and neither do we."""
    text, screen = boxed(["/mode", "\t", "\r"])
    assert text == "/mode"


def test_a_word_with_no_hits_stays_untouched(boxed):
    text, screen = boxed(["/zzz", "\t", "\r"])
    assert text == "/zzz"


def test_listing_is_bounded(boxed):
    """A completer with 200 runaway hits must not flood the screen."""
    from delfin.agent import repl_commands as rc
    assert len(rc.BUILTINS) > 20       # the bound is what this tests
    text, screen = boxed(["/", "\t", "\r"])
    listed = [r for r in screen.text() if r.lstrip().startswith("/")]
    assert len(listed) <= 20, len(listed)
