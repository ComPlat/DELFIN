"""JC2 target 2: the chain index — red test first.

A chain created by a real Pipeline.run must be recorded in one small
per-user index under the DELFIN state dir, carrying (name, work dir,
created_at, agent session id when given); ``delfin chain list`` shows
the indexed chains with their current state re-read from each chain's
OWN state file.  Control: red on bfec5344c (no index module, no list
subcommand, no agent_session_id kwarg).
"""

from __future__ import annotations

import json
import pathlib
from pathlib import Path

import pytest


@pytest.fixture()
def home(tmp_path, monkeypatch):
    """A redirected HOME so the per-user index never touches the real one."""
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path / "home")
    return tmp_path / "home"


def _fake_run_step(step_name, geometry=None, work_dir=None, **kwargs):
    from delfin.tools.pipeline import StepResult, StepStatus

    out = Path(work_dir) / f"{step_name}.out"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("fake run output\n")
    return StepResult(
        step_name=step_name, status=StepStatus.SUCCESS,
        work_dir=Path(work_dir),
    )


def _run_chain(wd, session_id=None):
    from delfin.tools.pipeline import Pipeline

    geom = Path(wd) / "start.xyz"
    geom.parent.mkdir(parents=True, exist_ok=True)
    geom.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", charge=0)
    pipe.run(work_dir=wd, geometry=geom, agent_session_id=session_id)
    return wd


def test_new_chain_is_recorded_with_session_id(home):
    from delfin.tools import chain_index

    wd = _run_chain(home / "runs" / "a", session_id="sess-abc")
    entries = chain_index.load()
    assert len(entries) == 1
    e = entries[0]
    assert e["name"] == "demo"
    assert e["work_dir"] == str(wd)
    assert e["agent_session_id"] == "sess-abc"
    assert isinstance(e["created_at"], float)


def test_run_without_session_id_indexes_nothing_named(home):
    """No agent session → the entry carries no session id (never an
    invented one), and the run completes (wave-15b: nothing new is
    refused)."""
    from delfin.tools import chain_index

    _run_chain(home / "runs" / "plain")
    entries = chain_index.load()
    assert len(entries) == 1
    assert entries[0]["agent_session_id"] is None


def test_index_filters_by_session(home):
    from delfin.tools import chain_index

    _run_chain(home / "runs" / "a", session_id="sess-1")
    _run_chain(home / "runs" / "b", session_id="sess-2")
    mine = chain_index.load(agent_session_id="sess-1")
    assert len(mine) == 1 and mine[0]["name"] == "demo"


def test_load_is_per_call_redirected_home(home):
    """The index path resolves per call under the redirected home — no
    import-time constant (the state-paths watcher rule)."""
    from delfin.tools import chain_index

    _run_chain(home / "runs" / "a", session_id="s")
    assert chain_index.index_path() == home / ".delfin" / "chain_index.jsonl"
    assert chain_index.index_path().is_file()


def test_torn_line_is_skipped_with_log(home, caplog):
    from delfin.tools import chain_index

    p = chain_index.index_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    good = {"name": "g", "work_dir": "/no/where", "agent_session_id": None,
            "created_at": 1.0}
    p.write_text(json.dumps(good) + "\n{torn\n", encoding="utf-8")
    import logging
    with caplog.at_level(logging.WARNING):
        entries = chain_index.load()
    assert len(entries) == 1
    assert any("torn" in r.message for r in caplog.records)


def test_bounded_prunes_finished_and_gone(home):
    """Old entries whose chain is finished or vanished are pruned on the
    next record; the index stays a pointer, never a truth store."""
    from delfin.tools import chain_index
    from delfin.tools import chain_state

    finished_wd = home / "runs" / "fin"
    chain_state.ChainState.create(name="fin", work_dir=finished_wd)
    chain = chain_state.ChainState.load(finished_wd / "chain_state.json")
    chain.add_step("s0", "xtb_opt")
    chain.steps[0].mark_running()
    chain.steps[0].mark_finished({"tool_reported_success": True})
    chain.save()
    chain_index.record("fin", finished_wd)
    chain_index.record("gone", home / "runs" / "vanished")
    _run_chain(home / "runs" / "new", session_id="s")
    names = [e["name"] for e in chain_index.load()]
    assert "gone" not in names, "vanished chains are pruned"
    assert "fin" not in names, "finished chains are pruned"
    assert "demo" in names


def test_chain_list_cli_shows_state_from_each_chain(home, capsys):
    """`delfin chain list` is read-only and re-reads each chain's own
    state file — the index is a pointer."""
    from delfin.cli_chain import main
    from delfin.tools import chain_state

    wd = _run_chain(home / "runs" / "a", session_id="sess-1")
    chain = chain_state.ChainState.load(wd / "chain_state.json")
    assert "pending" in [s.state for s in chain.steps]

    rc = main(["list"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "chain 'demo'" in out
    assert str(wd) in out
    assert "pending" in out

    # Empty index → one line, not silence.
    p = chain_index_index(home)
    p.unlink()
    rc = main(["list"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "no chains recorded" in out


def chain_index_index(home):
    from delfin.tools import chain_index
    return chain_index.index_path()
