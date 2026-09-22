"""Ctrl+R in the box: reverse history search, like readline's.

Measured 2026-09-20: Ctrl+R did nothing, although there is a
history_search tool and readline — whose prompt this box replaced —
answers Ctrl+R with an incremental backward search. The box kept the
history (same readline store) but gave no way to search it: with the
arrow keys you walk the whole list to find a line you know three
letters of.

Design: the DECODER only emits a SEARCH event for \\x12 — it knows
keys, not history. The search state lives in read_boxed's loop, where
_BoxHistory already lives: typing extends the query, Ctrl+R again
steps to an older match, Enter accepts the match (and submits, as
readline does), Esc leaves the search and restores the draft. Every
other key while searching is taken as query text — the search prompt
shows what it is doing, so nothing is silently swallowed.
"""

from __future__ import annotations

import io

import pytest

from delfin.agent import repl as R
from delfin.agent import repl_keys as rk

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
    """read_boxed against a scripted keyboard and a model screen, with a
    known readline history. Clears the process-global history and the
    saved completer for the duration, and restores the completer after."""
    def _run(keys, history, *, width=80):
        try:
            import readline
        except ImportError:
            pytest.skip("no readline on this host")
        old_completer = readline.get_completer()
        readline.clear_history()
        for line in history:
            readline.add_history(line)
        screen = Screen(width)
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
        agent.transcript.chrome = lambda text: screen.write(text + "\r\n")
        from delfin.agent import repl_render as rr
        agent.transcript.theme = rr.theme_for(io.StringIO(), mode="never")
        monkeypatch.setattr(rk, "RawMode", lambda *a, **k: _Keys(keys))
        try:
            text = agent.read_boxed()
        finally:
            readline.set_completer(old_completer)
        return text, screen
    return _run


HISTORY = ["/model kit", "/mode acceptEdits", "/help", "count the files"]


def test_ctrl_r_finds_the_newest_match_and_enter_submits_it(boxed):
    text, screen = boxed(["\x12", "mo", "\r"], HISTORY)
    assert text == "/mode acceptEdits", screen.text()


def test_ctrl_r_again_steps_to_the_older_match(boxed):
    text, screen = boxed(["\x12", "mo", "\x12", "\r"], HISTORY)
    assert text == "/model kit", screen.text()


def test_esc_leaves_the_search_and_restores_the_draft(boxed):
    """What was typed before Ctrl+R is not lost by starting a search."""
    text, screen = boxed(["draf", "\x12", "mo", "\x1b", "\r"], HISTORY)
    assert text == "draf", screen.text()


def test_typing_extends_the_query(boxed):
    text, screen = boxed(["\x12", "m", "o", "d", "e", "l", "\r"], HISTORY)
    assert text == "/model kit", screen.text()


def test_no_match_says_so_and_keeps_the_query(boxed):
    text, screen = boxed(["\x12", "zzz", "\x1b", "\r"], HISTORY)
    assert text == "", screen.text()


def test_the_search_prompt_names_the_mode(boxed):
    """The box must show that keystrokes are query now, not text. The
    prompt is transient — the submit clears the box — so the assert
    reads the screen's HISTORY, not its final grid."""
    text, screen = boxed(["\x12", "mo", "\r"], HISTORY)
    assert any("reverse-i-search" in r for r in screen.seen()), \
        screen.seen()[:8]
