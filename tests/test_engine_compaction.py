"""Tests for AgentEngine._compact_history CLI-restart behaviour.

S5 introduced the user-settings opt-out ``agent.compact_resets_cli``.
Default True (kill on compact) keeps the historical behaviour;
False keeps the subprocess alive and only mutates local messages.
"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.engine import AgentEngine


@pytest.fixture
def agent_tree(tmp_path):
    """Minimal pack tree so AgentEngine can boot.

    Mirrors the layout used in tests/test_agent_prompt_loader.py
    (pack/ for roles, pack_lite/ for modes manifest).
    """
    import textwrap

    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()

    (shared / "delfin_context.md").write_text("# DELFIN Context\nTest.")
    (shared / "work_cycle_rules.md").write_text("# Work Cycle Rules\nRule 1.")
    (shared / "universal_input_template.md").write_text("# Input Template")
    (shared / "minimal_final_verdict.md").write_text("# Verdict")
    (agents / "solo_agent.md").write_text("# Solo Agent\nYou are solo.")

    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# Mode: quick\nQuick mode.")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        recommended_default_mode: quick
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - solo_agent
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _build_engine(agent_tree, *, n_messages=20):
    """Boot an engine with a fake CLI client and a long history."""
    fake_client = MagicMock()
    fake_client.session_id = "sess-1"
    with patch("delfin.agent.engine.create_client", return_value=fake_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )
    # current_role is derived from the mode's route — solo_agent mode
    # gives us solo_agent as current_role automatically.
    # Compaction is token-driven (the legacy message-count trigger was
    # removed — it fired at ~15% window full). Shrink the window so the
    # long history below crosses the auto-compact threshold; the 12-message
    # floor still protects short conversations from compacting at all.
    engine.context_window_tokens = 10
    # Long alternating history (must exceed _COMPACTION_THRESHOLD = 12)
    for i in range(n_messages):
        role = "user" if i % 2 == 0 else "assistant"
        engine.messages.append({"role": role, "content": f"msg {i}"})
    return engine


def _settings_dir_with(setting_value, tmp_path, monkeypatch):
    """Point load_settings at a temp file with the given compact_resets_cli value."""
    cfg_dir = tmp_path / "settings"
    cfg_dir.mkdir(exist_ok=True)
    cfg_path = cfg_dir / ".delfin_settings.json"
    cfg_path.write_text(json.dumps({"agent": {"compact_resets_cli": setting_value}}))
    monkeypatch.setenv("HOME", str(cfg_dir))
    return cfg_path


def test_compact_kills_cli_when_setting_true_default(agent_tree):
    """Default behaviour: kill the CLI subprocess after compaction."""
    engine = _build_engine(agent_tree, n_messages=20)
    engine._compact_history()
    engine.client.kill.assert_called_once()


def test_compact_does_not_kill_cli_when_setting_false(
    agent_tree, tmp_path, monkeypatch,
):
    """compact_resets_cli=False keeps the subprocess alive."""
    _settings_dir_with(False, tmp_path, monkeypatch)
    engine = _build_engine(agent_tree, n_messages=20)
    engine._compact_history()
    engine.client.kill.assert_not_called()


def test_compact_local_messages_always_trimmed(agent_tree, tmp_path, monkeypatch):
    """Whether or not the CLI was killed, local messages are always compacted."""
    _settings_dir_with(False, tmp_path, monkeypatch)
    engine = _build_engine(agent_tree, n_messages=20)
    engine._compact_history()
    # Compaction → 1 summary user + 1 ack assistant + last 4 originals = 6
    assert len(engine.messages) == 6
    assert "Conversation summary" in engine.messages[0]["content"]


def test_compact_below_threshold_is_noop(agent_tree):
    """Short conversations don't compact at all."""
    engine = _build_engine(agent_tree, n_messages=4)
    engine._compact_history()
    engine.client.kill.assert_not_called()
    assert len(engine.messages) == 4


# ---------------------------------------------------------------------------
# Re-compaction fidelity: a prior summary must not be re-truncated to 400 chars
# (that silently dropped any fact sitting deep in it). Exercises the compaction
# methods directly via __new__ so no agent-pack is needed.
# ---------------------------------------------------------------------------


def _bare_engine():
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 10      # force compaction (tiny window)
    eng.auto_compact_pct = 0.95
    eng.backend = "cli"                 # extractive path (no LLM call)
    eng.client = None
    eng.last_compaction_info = {}
    eng.session_id = ""
    return eng


def test_prior_summary_survives_second_compaction():
    """A load-bearing fact sitting DEEP in the first summary must survive when
    a second compaction folds that summary back in. Regression: the prior
    summary was treated as a 400-char user goal, dropping everything past it."""
    eng = _bare_engine()
    eng.messages.append({"role": "user", "content": "ORIGINAL_GOAL: build the EOS solver"})
    for i in range(5):
        eng.messages.append({"role": "user", "content": f"step {i}: " + "work " * 18})
        eng.messages.append({"role": "assistant", "content": f"did {i}: " + "done " * 18})
    eng.messages.append({"role": "user",
                         "content": "DEEP_CONSTRAINT: bug reports group-readable to qmchem_shared only"})
    for i in range(5, 8):
        eng.messages.append({"role": "user", "content": f"step {i}: " + "work " * 18})
        eng.messages.append({"role": "assistant", "content": f"did {i}: " + "done " * 18})

    eng._compact_history()
    block1 = eng.messages[0]["content"]
    assert block1.startswith("[Conversation summary")
    # precondition: the constraint really is deep (past the old 400-char cut)
    assert "DEEP_CONSTRAINT" in block1 and block1.find("DEEP_CONSTRAINT") > 400

    for i in range(12):
        eng.messages.append({"role": "user" if i % 2 == 0 else "assistant",
                             "content": f"later {i}: " + "more " * 10})
    eng._compact_history()
    block2 = eng.messages[0]["content"]
    assert "ORIGINAL_GOAL" in block2      # original intent preserved
    assert "DEEP_CONSTRAINT" in block2    # deep fact no longer truncated away


# ---------------------------------------------------------------------------
# Pinned context regions: msg["_pinned"] protects a message from every
# destructive compaction stage (slide-trim, hard-clear, full summary), and
# the pin API validates bounds. Uses the same __new__-based bare engine as
# the re-compaction tests above.
# ---------------------------------------------------------------------------


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    """Isolate every session-store side file (elided store, archives)."""
    from delfin.agent import session_store as ss
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(ss, "_SESSIONS_DIR", tmp_path / ".delfin" / "agent_sessions")
    return tmp_path


def test_pin_api_bounds_validation():
    eng = _bare_engine()
    assert eng.pin_message(0) is False          # empty history
    eng.messages = [
        {"role": "user", "content": "a"},
        {"role": "assistant", "content": "b"},
        {"role": "user", "content": "c"},
    ]
    assert eng.pin_message(5) is False           # out of range
    assert eng.pin_message(-4) is False          # negative out of range
    assert eng.pin_message("x") is False         # malformed index
    assert eng.pinned_indices() == []
    assert eng.pin_message(1) is True
    assert eng.pin_message(-1) is True           # python-style negative
    assert eng.pinned_indices() == [1, 2]
    assert eng.unpin_message(0) is False         # was never pinned
    assert eng.unpin_message(9) is False         # out of range
    assert eng.unpin_message(1) is True
    assert eng.pinned_indices() == [2]


def test_pinned_message_survives_slide_trim_and_unpin_restores():
    eng = _bare_engine()
    eng.context_window_tokens = 1000
    pinned_body = "PINNED SPEC " + "p" * 4000
    eng.messages = [
        {"role": "user", "content": "goal"},
        {"role": "assistant", "content": pinned_body},
        {"role": "assistant", "content": "q" * 4000},
        # protected recent tail
        {"role": "user", "content": "r1"},
        {"role": "assistant", "content": "r2"},
        {"role": "user", "content": "r3"},
        {"role": "assistant", "content": "r4"},
    ]
    assert eng.pin_message(1) is True
    eng._slide_window_trim()
    assert eng.messages[1]["content"] == pinned_body          # verbatim
    assert "trimmed by sliding window" in eng.messages[2]["content"]
    # Unpin -> the message becomes trimmable again.
    assert eng.unpin_message(1) is True
    eng._slide_window_trim()
    assert "trimmed by sliding window" in eng.messages[1]["content"]


def test_pinned_message_survives_hard_clear():
    eng = _bare_engine()
    eng.context_window_tokens = 100
    pinned_body = "PINNED OUTPUT " + "z" * 4000
    eng.messages = [
        {"role": "assistant", "content": pinned_body, "_pinned": True},
        {"role": "assistant", "content": "w" * 4000},
    ]
    n = eng._hard_clear_old_tool_results(eng.messages)
    assert n == 1
    assert eng.messages[0]["content"] == pinned_body          # verbatim
    assert eng.messages[1]["content"].startswith("[cleared:")


def test_pinned_message_survives_full_summary_compaction_verbatim():
    eng = _bare_engine()
    pinned_body = "PINNED CALC SPEC: CASSCF(12,12)/def2-TZVP, S1 at 2.31 eV"
    eng.messages = [{"role": "user", "content": "ORIGINAL_GOAL: solver"}]
    for i in range(6):
        eng.messages.append({"role": "assistant", "content": f"did {i}: " + "done " * 18})
        eng.messages.append({"role": "user", "content": f"step {i}: " + "work " * 18})
    eng.messages.insert(3, {"role": "assistant", "content": pinned_body})
    assert eng.pin_message(3) is True
    eng._compact_history()
    # The summary block replaced the compactable middle ...
    assert eng.messages[0]["content"].startswith("[Conversation summary")
    # ... but the pinned message came through VERBATIM, still pinned.
    kept = [m for m in eng.messages if m.get("content") == pinned_body]
    assert len(kept) == 1 and kept[0].get("_pinned") is True
    assert eng.last_compaction_info.get("pinned_kept") == 1
    # Alternation still holds after the rebuild + sanitize.
    for a, b in zip(eng.messages, eng.messages[1:]):
        assert not (a["role"] == b["role"] and a["role"] in ("user", "assistant"))


def test_sanitize_keeps_pinned_on_same_role_merge():
    eng = _bare_engine()
    eng.messages = [
        {"role": "user", "content": "PINNED SPEC", "_pinned": True},
        {"role": "user", "content": "newer user message"},
    ]
    eng._sanitize_messages()
    contents = [m["content"] for m in eng.messages]
    assert "PINNED SPEC" in contents
    assert "newer user message" in contents
    roles = [m["role"] for m in eng.messages]
    assert roles == ["user", "assistant", "user"]   # filler restores alternation


def test_context_status_block_reports_pins():
    eng = _bare_engine()
    eng.token_usage = {"input": 0, "output": 0}
    eng.messages = [
        {"role": "user", "content": "a"},
        {"role": "assistant", "content": "b"},
    ]
    assert "pinned message(s)" not in eng._build_context_status_block()
    eng.pin_message(0)
    eng.pin_message(1)
    block = eng._build_context_status_block()
    assert "2 pinned message(s) excluded from compaction" in block


def test_wire_messages_strips_private_keys_only_when_present():
    eng = _bare_engine()
    eng.messages = [
        {"role": "user", "content": "a"},
        {"role": "assistant", "content": "b"},
    ]
    # No private keys -> identity (prompt-cache/list-binding friendly).
    assert eng._wire_messages() is eng.messages
    eng.pin_message(1)
    wire = eng._wire_messages()
    assert wire is not eng.messages
    assert [m["content"] for m in wire] == ["a", "b"]
    assert all("_pinned" not in m for m in wire)
    assert eng.messages[1]["_pinned"] is True       # source untouched


def test_pin_flag_round_trips_session_store(fake_home):
    from delfin.agent import session_store as ss
    ss.save_session(
        "sess-pin-rt",
        engine_messages=[
            {"role": "user", "content": "keep", "_pinned": True},
            {"role": "assistant", "content": "ok"},
        ],
    )
    data = ss.load_session("sess-pin-rt")
    assert data["engine_messages"][0]["_pinned"] is True
    assert data["engine_messages"][1].get("_pinned") is None


# ---------------------------------------------------------------------------
# Lossless elision: the in-place trim stages persist the FULL original to
# <sessions_dir>/<session_id>.elided.jsonl and the marker carries the ref.
# ---------------------------------------------------------------------------

_REF_RE = r"history_get\('elided:([A-Za-z0-9]+)'\)"


def test_slide_trim_writes_elided_record_with_matching_ref(fake_home):
    import re as _re
    from delfin.agent import history_search as hs
    from delfin.agent import session_store as ss
    eng = _bare_engine()
    eng.session_id = "sess-elide"
    eng.context_window_tokens = 1000
    original = "ELIDE ME " + "e" * 4000
    eng.messages = [
        {"role": "assistant", "content": original},
        {"role": "user", "content": "r1"},
        {"role": "assistant", "content": "r2"},
        {"role": "user", "content": "r3"},
        {"role": "assistant", "content": "r4"},
    ]
    assert eng._slide_window_trim() == 1
    marker = eng.messages[0]["content"]
    m = _re.search(_REF_RE, marker)
    assert m, f"marker lacks retrieval ref: {marker[:200]}"
    ref = m.group(1)
    rec = ss.load_elided_record("sess-elide", ref)
    assert rec["content"] == original
    assert rec["index"] == 0 and rec["role"] == "assistant"
    assert rec["reason"] == "sliding_window"
    # history_get resolves the marker's ref to the full original.
    got = hs.history_get("sess-elide", f"elided:{ref}", max_chars=10_000)
    assert got["text"] == original and got["source"] == "elided"


def test_hard_clear_writes_elided_record_with_matching_ref(fake_home):
    import re as _re
    from delfin.agent import history_search as hs
    eng = _bare_engine()
    eng.session_id = "sess-clear"
    eng.context_window_tokens = 100
    original = "CLEARED BODY " + "c" * 3000
    eng.messages = [{"role": "assistant", "content": original}]
    assert eng._hard_clear_old_tool_results(eng.messages) == 1
    marker = eng.messages[0]["content"]
    assert marker.startswith("[cleared:")
    m = _re.search(_REF_RE, marker)
    assert m, f"marker lacks retrieval ref: {marker[:200]}"
    got = hs.history_get("sess-clear", f"elided:{m.group(1)}", max_chars=10_000)
    assert got["text"] == original and got["reason"] == "hard_clear"


def test_elision_store_failure_never_breaks_compaction(fake_home, monkeypatch):
    from delfin.agent import session_store as ss

    def _boom(*a, **k):
        raise RuntimeError("disk full")

    monkeypatch.setattr(ss, "append_elided_record", _boom)
    eng = _bare_engine()
    eng.session_id = "sess-broken"
    eng.context_window_tokens = 1000
    eng.messages = [
        {"role": "assistant", "content": "x" * 4000},
        {"role": "user", "content": "r1"},
        {"role": "assistant", "content": "r2"},
        {"role": "user", "content": "r3"},
        {"role": "assistant", "content": "r4"},
    ]
    assert eng._slide_window_trim() == 1          # trim still happened
    marker = eng.messages[0]["content"]
    assert "trimmed by sliding window" in marker
    assert "elided:" not in marker                 # no dangling ref
    # Hard-clear path degrades the same way.
    eng2 = _bare_engine()
    eng2.session_id = "sess-broken"
    eng2.context_window_tokens = 100
    eng2.messages = [{"role": "assistant", "content": "y" * 3000}]
    assert eng2._hard_clear_old_tool_results(eng2.messages) == 1
    assert "elided:" not in eng2.messages[0]["content"]


# ---------------------------------------------------------------------------
# Capability probe is pre-warmed off the hot path: constructing the engine must
# NOT block on the (cold-cache, ~5s on KIT) /v1/models probe — it runs on a
# daemon thread so the first turn never stalls. The window upgrades from the
# 100k default once the probe lands.
# ---------------------------------------------------------------------------

def test_context_window_probe_runs_in_background(agent_tree, monkeypatch):
    import threading
    import time
    from types import SimpleNamespace
    from delfin.agent import model_capabilities as mc

    # A probe that blocks until the test releases it. If construction were
    # synchronous it would hang here; the daemon thread lets construction return.
    release = threading.Event()
    caps = SimpleNamespace(context_window=262144, max_output_tokens=8192,
                           supports_tools=True, supports_vision=False)

    def gated_resolve(*a, **k):
        release.wait(timeout=5)
        return caps

    monkeypatch.setattr(mc, "resolve", gated_resolve)
    fake_client = MagicMock()
    fake_client.session_id = "s"

    with patch("delfin.agent.engine.create_client", return_value=fake_client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)

    # Construction returned while the probe is still blocked → window untouched.
    assert engine.context_window_tokens == 100_000

    # Let the probe finish; it upgrades the window for the next turn.
    release.set()
    deadline = time.monotonic() + 3.0
    while (time.monotonic() < deadline
           and engine.context_window_tokens != 262144):
        time.sleep(0.02)
    assert engine.context_window_tokens == 262144
