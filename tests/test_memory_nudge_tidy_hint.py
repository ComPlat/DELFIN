"""The mid-session nudge also surfaces pending tidy proposals.

Phase-1 finding (2026-09-29, verified): ``memory_tidy`` closes the
duplication gap in the fact store automatically (merge/retire proposals),
but it is reachable ONLY through a manual CLI command nobody runs —
measured 0 uses across 25 real sessions. The nudge is the natural
carrier: it already fires mid-session when work has accumulated, exactly
when near-duplicate facts were just written. When the store has pending
tidy proposals, the nudge mentions them; when it does not, the text is
unchanged. Nothing is applied automatically — retire still moves notes
aside only on explicit acceptance, and the proposal stays read-only.

Pins ``delfin.agent.memory_nudge.compose_nudge``; tests never touch the
real ``~/.delfin`` (tidy proposals are computed against a tmp store).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _make_note(store: Path, name: str, text: str, mtime: int | None = None):
    store.mkdir(parents=True, exist_ok=True)
    p = store / name
    p.write_text(
        "---\nname: n\ntype: project\n---\n\n" + text, encoding="utf-8")
    if mtime is not None:
        import os
        os.utime(p, (mtime, mtime))
    return p


def _nudge(store: Path | None = None):
    from delfin.agent import memory_nudge
    return memory_nudge.compose_nudge(store=store)


def test_nudge_without_store_or_proposals_is_the_plain_text(tmp_path):
    """No store, no proposals: the nudge is the unchanged classic text."""
    out = _nudge(store=tmp_path / "empty")
    assert "worth keeping" in out
    assert "tidy" not in out.lower()


def test_nudge_with_pending_tidy_proposal_mentions_it(tmp_path):
    """A near-duplicate pair in the store adds the /tidy pointer.

    Two notes with almost identical wording at a months-old mtime are a
    merge candidate by memory_tidy's similarity rules; the nudge must
    say that tidy proposals are pending and that /tidy changes nothing.
    """
    store = tmp_path / "notes"
    _make_note(store, "project_one.md",
               "the coordinator rebuild needs anchoring")
    _make_note(store, "project_two.md",
               "the coordinator rebuild needs anchoring again")
    out = _nudge(store=store)
    assert "worth keeping" in out          # the classic ask stays
    assert "/tidy" in out                  # the new pointer
    assert "changes nothing" in out        # the safety promise stays


def test_nudge_failure_degrades_to_plain_text(tmp_path, monkeypatch):
    """A tidy pass that explodes must not take the nudge with it."""
    from delfin.agent import memory_nudge, memory_tidy
    monkeypatch.setattr(memory_tidy, "hint",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
    out = memory_nudge.compose_nudge(store=tmp_path)
    assert "worth keeping" in out
