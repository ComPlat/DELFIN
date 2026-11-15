"""Live run 2 (operator, 2026-11-16): the offer must not depend on the
turn checkpoint, and `delfin chain list` must show workspace chains.

Run 2's crash happened BEFORE the first turn checkpoint was written, so
consume_crash_recovery_note returned None and the restarted agent saw
no offer at all (offer command lines: 0) — the exact case JC2 exists to
prevent.  Red on the pre-fix code (f00b2101c session_store/cli_chain,
same defect as the operator's 675a66374).
"""
import json
import time

import pytest

from delfin.tools import chain_state


@pytest.fixture()
def home(tmp_path, monkeypatch):
    """Redirect BOTH per-user stores (index via Path.home, session
    files via _SESSIONS_DIR — the latter is an import-time constant, so
    the Path.home patch alone does not move it)."""
    home = tmp_path / "home"
    home.mkdir(parents=True)
    monkeypatch.setattr("pathlib.Path.home", lambda: home)
    import delfin.agent.session_store as ss
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", home / ".delfin" / "agent_sessions")
    return home


def _make_session(home, session_id, workspace):
    """A saved session record with a workspace (partial-note shape)."""
    import delfin.agent.session_store as ss

    ss._SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
    now = time.time()
    rec = {
        "session_id": session_id,
        "workspace": str(workspace),
        "created_at": now,
        "updated_at": now + 10,
        "messages": [],
        "turn_count": 1,
    }
    (ss._SESSIONS_DIR / f"{session_id}.json").write_text(json.dumps(rec))


def _run_chain(monkeypatch, run_dir, finish=False):
    """A real chain round trip in run_dir (only run_step stood in)."""
    from delfin.tools import pipeline as pl

    def fake_run_step(spec, work_dir, **kw):
        return {"artifacts": {"gbw": str(work_dir / "s.gbw")}}

    monkeypatch.setattr(pl, "run_step", fake_run_step)
    pipe = pl.Pipeline("jc2-live2")
    pipe.work_dir = run_dir
    pipe.add("xtb_opt", {})
    pipe.add("orca_sp", {"charge": 0, "method": "PBE0"})
    if finish:
        monkeypatch.setattr(pl, "run_step", fake_run_step)
    run_dir.mkdir(parents=True, exist_ok=True)
    try:
        pipe.run(work_dir=run_dir)
    except Exception:
        pass
    return run_dir


def _offer(home, monkeypatch, session_id, cwd):
    monkeypatch.chdir(cwd)
    import delfin.agent.session_store as ss

    return ss.consume_crash_recovery_note(session_id)


def _unfinished_chain(monkeypatch, run_dir):
    """Chain with one pending step (crash before the step started)."""

    run_dir.mkdir(parents=True, exist_ok=True)
    chain = chain_state.ChainState(name="jc2-live2", work_dir=run_dir)
    chain.steps.append(chain_state.ChainStep(
        name="s1", tool="orca_opt", state="submitted"))
    chain.steps.append(chain_state.ChainStep(
        name="s2", tool="orca_sp", state="pending"))
    chain.save()
    return run_dir


def test_offer_without_turn_checkpoint(home, monkeypatch):
    session_id = "sess-no-ck"
    ws = home / "ws"
    _make_session(home, session_id, ws)
    _unfinished_chain(monkeypatch, ws / "chain1")
    note = _offer(home, monkeypatch, session_id, ws)
    assert note is not None, "no offer without a turn checkpoint"
    assert "chain1" in note
    assert "delfin chain resume" in note


def test_with_checkpoint_exactly_one_offer(home, monkeypatch):
    import delfin.agent.session_store as ss

    session_id = "sess-ck"
    ws = home / "ws2"
    _make_session(home, session_id, ws)
    _unfinished_chain(monkeypatch, ws / "chain1")
    ss.save_turn_checkpoint(
        session_id,
        {"user_message": "run the chain", "tool_calls": 1,
         "partial_response": "", "ts": time.time() + 100},
    )
    note = _offer(home, monkeypatch, session_id, ws)
    assert note is not None
    assert note.count("delfin chain resume") == 1, note


def test_finished_chain_second_resume_no_offer(home, monkeypatch):
    session_id = "sess-done"
    ws = home / "ws3"
    _make_session(home, session_id, ws)
    wd = _unfinished_chain(monkeypatch, ws / "chain1")
    chain = chain_state.ChainState.load(wd / "chain_state.json")
    for s in chain.steps:
        # mark finished through the real state machine: a stored
        # finished step without evidence is downgraded to pending on
        # load (torn-state rule), so the fixture must carry evidence.
        s.state = "finished"
        s.evidence = {"tool_reported_success": True}
    chain.save()
    note = _offer(home, monkeypatch, session_id, wd)
    assert note is None or "chain1" not in note


def test_second_resume_of_unfinished_session_offers_again(home, monkeypatch):
    """Pinned contract (QS ruling 8, s8 2026-11-16): once per RESUME,
    not once per chain — no consumed-marker is stored, so a second
    resume of a still-unfinished session rebuilds the offer (the
    finished-chain filter stops stale repeats)."""
    session_id = "sess-again"
    ws = home / "ws5"
    _make_session(home, session_id, ws)
    _unfinished_chain(monkeypatch, ws / "chain1")
    first = _offer(home, monkeypatch, session_id, ws)
    assert first is not None and "chain1" in first
    second = _offer(home, monkeypatch, session_id, ws)
    assert second is not None and "chain1" in second


def test_normal_turn_of_session_without_chains_no_chain_text(home, monkeypatch):
    """Pins the operator's 'no chain text at all for chainless sessions'
    on the CHECKPOINT path too: a normal turn's checkpoint of a session
    with no chains yields a partial-note without any chain text."""
    import delfin.agent.session_store as ss

    session_id = "sess-plain"
    ws = home / "ws6"
    ws.mkdir()
    _make_session(home, session_id, ws)
    ss.save_turn_checkpoint(
        session_id,
        {"user_message": "hello", "tool_calls": 1,
         "partial_response": "all done", "ts": time.time() + 100},
    )
    note = _offer(home, monkeypatch, session_id, ws)
    assert note is not None
    assert "delfin chain resume" not in note
    assert "chain" not in note


def test_stale_checkpoint_still_offers(home, monkeypatch):
    """Operator run-3 review, defect 1: a stale checkpoint (a full save
    landed after it) with an unfinished chain must return the OFFER —
    not None — and without any "[recovered]" text."""
    session_id = "sess-stale"
    ws = home / "ws7"
    _make_session(home, session_id, ws)
    _unfinished_chain(monkeypatch, ws / "chain1")
    # A checkpoint OLDER than the session's last save is stale.
    import delfin.agent.session_store as ss

    ss.save_turn_checkpoint(
        session_id,
        {"user_message": "run the chain", "tool_calls": 2,
         "partial_response": "", "ts": 1.0},
    )
    note = _offer(home, monkeypatch, session_id, ws)
    assert note is not None, "stale checkpoint dropped the offer"
    assert "chain1" in note
    assert "[recovered]" not in note


def test_chain_list_shows_workspace_chain(home, monkeypatch, capsys):
    ws = home / "ws4"
    ws.mkdir()
    _unfinished_chain(monkeypatch, ws / "chain1")
    monkeypatch.chdir(ws)
    import delfin.cli_chain as cc

    rc = cc.main(["list"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "chain1" in out
    assert "not indexed" in out
