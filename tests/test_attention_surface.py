"""User surface of the attention inbox (render/answer/dismiss/clear).

Round-trips go through the real JSONL storage (Path.home monkeypatched),
including the consumption path the engine uses (``drain_resolved``):
an answered item must reach the model exactly once, a dismissed item
never. Every surface helper must be unable to raise.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import delfin.agent.attention as attention


@pytest.fixture(autouse=True)
def _fake_home(tmp_path, monkeypatch):
    """Isolate ~ (inbox + settings) and mute out-of-band transports."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: home)
    import delfin.agent.notify as notify
    monkeypatch.setattr(notify, "send_notification", lambda *a, **k: True)
    monkeypatch.setattr(
        notify, "send_remote_trigger",
        lambda *a, **k: notify.TriggerResult(sent=False, error="test"))
    return home


def _seed_three() -> dict[str, str]:
    ids = {}
    ids["q"] = attention.emit_attention(
        "question_pending", session_id="s1", title="Which basis set?",
        options=["def2-SVP", "def2-TZVP"])
    ids["c"] = attention.emit_attention(
        "confirm_pending", title="Confirmation required: bash",
        detail="rm -rf build")
    ids["f"] = attention.emit_attention(
        "run_finished", title="entry benzene finished")
    return ids


# ---------------------------------------------------------------------------
# render_inbox
# ---------------------------------------------------------------------------

def test_render_empty_inbox(_fake_home):
    out = attention.render_inbox()
    assert isinstance(out, str)
    assert "no pending items" in out


def test_render_groups_ids_counts_and_options(_fake_home):
    ids = _seed_three()
    out = attention.render_inbox()
    assert "3 pending" in out
    assert "2 blocking the agent" in out          # question + confirm
    for eid in ids.values():
        assert eid in out                          # ids visible (full)
    # Grouped, blocking kinds first: confirm before question before finished.
    assert (out.index("Confirmations waiting (1):")
            < out.index("Questions waiting (1):")
            < out.index("Finished runs (1):"))
    assert "options: def2-SVP | def2-TZVP" in out
    assert "/attention answer <id> <text>" in out  # usage footer


def test_render_kind_filter_and_limit(_fake_home):
    _seed_three()
    only = attention.render_inbox("run_finished")
    assert "1 pending of kind 'run_finished'" in only
    assert "Questions waiting" not in only
    capped = attention.render_inbox(limit=1)
    assert "... and 2 more" in capped


# ---------------------------------------------------------------------------
# answer round-trip through storage + engine consumption path
# ---------------------------------------------------------------------------

def test_answer_round_trip_reaches_drain_exactly_once(_fake_home):
    ids = _seed_three()
    res = attention.answer_item(ids["q"], "def2-TZVP")
    assert res["ok"] is True
    assert res["id"] == ids["q"]
    assert res["kind"] == "question_pending"
    # No longer pending, and the engine-side consumption
    # (drain_resolved, see engine._build_answered_attention_block)
    # sees the answer exactly once.
    assert ids["q"] not in {e["id"] for e in attention.list_pending()}
    drained = attention.drain_resolved("s1")
    assert [(e["id"], e["answer"]) for e in drained] == \
        [(ids["q"], "def2-TZVP")]
    assert attention.drain_resolved("s1") == []


def test_answer_accepts_unique_id_prefix(_fake_home):
    eid = attention.emit_attention("plan_pending", title="plan v2")
    res = attention.answer_item(eid[:12], "approve")
    assert res["ok"] is True and res["id"] == eid


def test_answer_error_paths_never_raise(_fake_home):
    ids = _seed_three()
    assert attention.answer_item("", "text")["ok"] is False
    assert attention.answer_item("att-nope", "text")["ok"] is False
    assert "empty answer" in attention.answer_item(ids["q"], "  ")["error"]
    # "att-" prefixes every id -> ambiguous.
    amb = attention.answer_item("att-", "text")
    assert amb["ok"] is False and "ambiguous" in amb["error"]
    # Items stayed pending through all failed attempts.
    assert len(attention.list_pending()) == 3


# ---------------------------------------------------------------------------
# dismiss / clear_all — cleared items never replay to the model
# ---------------------------------------------------------------------------

def test_dismiss_item_not_replayed(_fake_home):
    ids = _seed_three()
    res = attention.dismiss_item(ids["c"])
    assert res["ok"] is True and res["id"] == ids["c"]
    assert ids["c"] not in {e["id"] for e in attention.list_pending()}
    # acknowledged=True: drain (any session) must never see it.
    assert attention.drain_resolved("") == []
    res2 = attention.dismiss_item(ids["c"])
    assert res2["ok"] is False                     # already gone


def test_clear_all_and_kind_scoped_clear(_fake_home):
    _seed_three()
    res = attention.clear_all("run_finished")
    assert (res["ok"], res["cleared"], res["kept"]) == (True, 1, 0)
    assert len(attention.list_pending()) == 2
    # A blanket clear leaves the two blocking items alone; naming the
    # kind is the explicit gesture that clears them.
    res = attention.clear_all()
    assert (res["ok"], res["cleared"], res["kept"]) == (True, 0, 2)
    assert len(attention.list_pending()) == 2
    res = attention.clear_all(include_blocking=True)
    assert (res["ok"], res["cleared"], res["kept"]) == (True, 2, 0)
    assert attention.list_pending() == []
    assert attention.drain_resolved("s1") == []    # nothing replays
    assert attention.clear_all()["cleared"] == 0


# ---------------------------------------------------------------------------
# never-raises on a broken storage layer
# ---------------------------------------------------------------------------

def test_surface_never_raises_when_storage_breaks(_fake_home, monkeypatch):
    _seed_three()

    def _boom():
        raise RuntimeError("disk gone")

    monkeypatch.setattr(attention, "_load_events", _boom)
    out = attention.render_inbox()
    assert isinstance(out, str) and "unavailable" in out
    for res in (
        attention.answer_item("att-x", "hi"),
        attention.dismiss_item("att-x"),
        attention.clear_all(),
    ):
        assert res["ok"] is False
        assert "disk gone" in res["error"]
