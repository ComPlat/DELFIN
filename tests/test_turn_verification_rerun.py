"""Characterization: the dashboard reruns a turn that cites a path that
does not exist; the terminal does not.

The dashboard's Verify-Enforcement block (the inline scan in
delfin/dashboard/tab_agent.py) forces one correction turn when the answer
cites a path the repo does not contain (``scan_for_ungrounded_code_claims``
returns a ``kind == "nonexistent"`` flag) and the turn never read a file or
consulted docs. The terminal layer drives a single ``stream_response`` and
does no post-answer scan of its own: the engine's claim-grounding guard
(``_enforce_claim_grounding``) uses ``scan_for_ungrounded_location_claims``,
which does not flag a nonexistent path, so a terminal turn that cites one
is answered once and never corrected.

This file pins that divergence (dashboard reruns, terminal does not) before
the verification decision moves into the engine. The terminal expectation
is intentional: today's absence is the behaviour the step is meant to
change, and it flips to the dashboard's count once the engine owns the
decision.
"""

from __future__ import annotations

import json
import textwrap

import pytest
from unittest.mock import MagicMock, patch

from tests.agent_terminal_probes import fake_model  # noqa: E402

# A path that does not exist in any repo.  ``verify_guard`` swallows paths
# that look like they might exist in a sibling root, so the name is
# deliberately inside a root that is not on disk.
NONEXISTENT = "delfin/src/nonexistent_file.py"
# A bare path citation with NO line number: the engine's location scanner
# only flags ``file_line`` claims, so this form slips past the terminal's
# guard while the dashboard's code scanner still flags it ``nonexistent``.
CITING = f"See `{NONEXISTENT}` for the loader."


@pytest.fixture
def agent_tree(tmp_path):
    """A minimal agent pack, the shape the engine tests construct."""
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


# --------------------------------------------------------------------------
# terminal path -- today: no rerun
# --------------------------------------------------------------------------

def _client(replies):
    """A fake backend that returns ``replies`` in order and counts calls."""
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*a, **k):
        i = calls["n"]
        calls["n"] += 1
        yield from _text_delta(replies[min(i, len(replies) - 1)])

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


def _text_delta(text):
    from delfin.agent.api_client import StreamEvent
    yield StreamEvent(type="text_delta", text=text)
    yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)


def _term_engine(client, agent_tree):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def test_terminal_does_not_rerun_a_nonexistent_path(agent_tree):
    """Today the engine's grounding guard does not treat a nonexistent
    citation as a correction trigger: one stream, one call, no [Verify]
    feedback appended. This is the behaviour the step is meant to change."""
    client = _client([CITING])
    engine = _term_engine(client, agent_tree)
    out = engine.stream_response("What is the fix?")
    assert client.stream_message.call_count == 1, (
        "terminal must not force a correction turn for a nonexistent path")
    assert "[Verify]" not in out
    feedback = [m for m in engine.messages
                if m.get("_harness") and str(m.get("content", "")).startswith("[Verify]")]
    assert feedback == [], "no verification rerun in the terminal today"


# --------------------------------------------------------------------------
# dashboard path -- today: one rerun
# --------------------------------------------------------------------------

def _walk(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "children", ()) or ())


def test_dashboard_reruns_a_nonexistent_path(tmp_path, monkeypatch):
    """The dashboard's inline Verify-Enforcement block forces one
    correction turn when the answer cites a path that does not exist.  The
    appended [Verify] feedback is a harness row on the engine's message
    list, and the corrected round is streamed again."""
    import shutil
    import ipywidgets as widgets  # noqa: F401
    from delfin import user_settings
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    MARK = "run-scenario"

    _which = shutil.which
    monkeypatch.setattr(shutil, "which", lambda name, *a, **k: (
        None if name in ("claude", "codex") else _which(name, *a, **k)))
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    ws = tmp_path / "ws"
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir()

    # Step 0 answers citing a nonexistent path; any later round (the
    # forced correction) says it has nothing more to add.
    def _plan(messages):
        starts = [i for i, m in enumerate(messages)
                  if m.get("role") == "user" and MARK in fake_model._text_of(m)]
        done = sum(1 for m in messages[starts[-1] if starts else 0:]
                   if m.get("role") == "assistant")
        if done >= 1:
            return {"text": "Nothing more."}
        return {"text": CITING}

    monkeypatch.setattr(fake_model, "plan", _plan)
    srv, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
    monkeypatch.setenv("OLLAMA_HOST", f"http://127.0.0.1:{port}")
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")

    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "ollama", "model": fake_model.MODEL, "mode": "solo",
        "permission_profile": "all_free", "probe_endpoint": False}}))

    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office", agent_dir=ws)
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(ws)
    tab, refs = tab_agent.create_tab(ctx)
    try:
        state = refs["state"]
        drop = {w.description: w for w in _walk(tab)
                if isinstance(w, widgets.Dropdown)}
        drop["Mode:"].value = "solo"
        inp = next(w for w in _walk(tab)
                   if isinstance(w, widgets.Textarea)
                   and "delfin-agent-input" in w._dom_classes)
        send = next(w for w in _walk(tab)
                    if isinstance(w, widgets.Button)
                    and w.description == "Send")
        inp.value = f"{MARK} please"
        send.click()
        thread = state.get("_turn_thread")
        assert thread is not None, "the send started no turn"
        thread.join(120)
        assert not thread.is_alive(), "the turn did not end"

        engine = state.get("engine")
        assert engine is not None, "the turn created no engine"

        verify = [m["content"] for m in engine.messages
                  if m.get("_harness")
                  and m.get("role") == "user"
                  and str(m["content"]).startswith("[Verify]")]
        # Today the dashboard's inline block injects exactly one [Verify]
        # citation correction and re-streams the turn (a second model call).
        assert len(verify) == 1, f"dashboard verification reruns: {len(verify)}"
        assert NONEXISTENT in verify[0]
    finally:
        refs["shutdown"](save=False)
        srv.shutdown()
