"""The inbox outlived every session that filled it, and nothing emptied it.

`/attention` listed. That was all it did. The inbox is per user and
survives restarts by design — which is right for a request somebody
still has to answer, and wrong for one whose session ended hours ago.

Reported on 2026-09-18: seven confirmations from dashboard sessions that
had long finished were still pending, and a fresh session was told they
were waiting for it. The reminder that told it so has since learned
whose request is whose; this is the other half — the user can throw away
what they no longer want.

The module's own rule is kept: a blanket clear touches notices only, so
the gesture for tidying run notices cannot also throw away the
confirmation an agent is blocked on. Naming the kind is how you mean it.
"""

from __future__ import annotations

import pytest

from delfin.agent import repl_commands as RC


@pytest.fixture()
def inbox(monkeypatch):
    calls: list = []

    class _Att:
        ATTENTION_KINDS = ("confirm_pending", "run_finished", "run_failed")

        @staticmethod
        def render_inbox(kind=None):
            return f"inbox({kind})"

        @staticmethod
        def clear_all(kind=None, *, include_blocking=False):
            calls.append((kind, include_blocking))
            return {"ok": True, "cleared": 7 if kind else 3,
                    "kept": 0 if kind else 7}

    import delfin.agent.attention as real
    for name in ("ATTENTION_KINDS", "render_inbox", "clear_all"):
        monkeypatch.setattr(real, name, getattr(_Att, name), raising=False)
    return calls


def test_listing_still_lists(inbox):
    assert "inbox(None)" in RC._attention(None, "").output
    assert "inbox(confirm_pending)" in RC._attention(
        None, "confirm_pending").output


def test_clear_empties_the_notices(inbox):
    out = RC._attention(None, "clear").output
    assert "cleared 3 item(s)" in out
    assert inbox == [(None, False)], "a blanket clear must not be blanket"


def test_a_blanket_clear_says_what_it_kept(inbox):
    out = RC._attention(None, "clear").output
    assert "kept 7" in out
    assert "name the kind" in out, "and how to mean it"


def test_naming_the_kind_clears_what_blocks(inbox):
    out = RC._attention(None, "clear confirm_pending").output
    assert "cleared 7 item(s)" in out
    assert inbox == [("confirm_pending", True)]


def test_an_unknown_kind_says_which_exist(inbox):
    out = RC._attention(None, "clear nonsense").output
    assert "usage: /attention clear" in out
    assert "confirm_pending" in out
    assert inbox == [], "nothing was cleared on a typo"


def test_the_usage_line_mentions_clear(inbox):
    out = RC._attention(None, "nonsense").output
    assert "/attention clear" in out


def test_it_never_raises(monkeypatch):
    import delfin.agent.attention as real

    def _boom(*a, **k):
        raise RuntimeError("no inbox")
    monkeypatch.setattr(real, "clear_all", _boom, raising=False)
    out = RC._attention(None, "clear").output
    assert "unavailable" in out
