"""A long unattended run is bounded by what it declared, and a ceiling it
meets is said, saved and resumable.

Input: an engine with the unattended defaults (``run_limits``), the model
it runs on (priced, or a KIT endpoint that bills nothing), a batch-job
environment, an auto-mode run or loop started in that session (its own
iterations / wall / tokens budget), and a ceiling reached during a turn.
Output: the ceilings in force, the loop's start line, the auto-mode stop
report and the attention items.

Semantics:

* the USD default applies only where spend is measured at a published
  rate; on a KIT endpoint (no USD) the time and turn ceilings bound it;
* inside a batch job (SLURM/PBS/LSF) the job's own time limit bounds the
  wall clock, so the wall-clock default is not added;
* an auto-mode run or loop declares a bounded budget of its own; the
  defaulted session ceilings are widened to cover it (the USD default is
  lifted: the loop's token budget bounds it) and its start line says so;
  a ceiling the user set stays as set;
* a session ceiling reached while auto mode runs stops auto mode with the
  setting and the resume command in its stop report (loop log, attention
  inbox); a ceiling reached inside a turn raises one attention item
  naming the resume command (``--session`` plus the flag that raises it).
"""
from __future__ import annotations

import time
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent import auto_mode as am
from delfin.agent import cli as agent_cli
from delfin.agent import engine as E
from delfin.agent import run_limits as RL


@pytest.fixture(autouse=True)
def _no_batch(monkeypatch):
    for var in ("SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID"):
        monkeypatch.delenv(var, raising=False)


def _settings(agent=None):
    return patch("delfin.user_settings.load_settings",
                 lambda *a, **k: {"agent": dict(agent or {})})


def _engine(model="gpt-5.4", provider="openai"):
    eng = E.AgentEngine.__new__(E.AgentEngine)
    eng.client = SimpleNamespace(
        model=model,
        _permissions=SimpleNamespace(mode="acceptEdits", auto_veto=None))
    eng.provider = provider
    eng.session_id = "sess-long-1"
    eng.repo_dir = "/tmp/ws"
    eng.cost_usd = 0.0
    eng._run_started_at = time.time()
    eng._budget_attention_levels = set()
    eng.token_usage = {"input": 0, "output": 0}
    eng.messages = []
    return eng


def test_the_usd_default_only_where_spend_is_measured():
    priced, kit = _engine(), _engine("kit.glm-4.6", "kit")
    with _settings():
        RL.apply_unattended(priced)
        notes = RL.apply_unattended(kit)
    assert priced.run_budget_usd == RL.UNATTENDED_DEFAULTS["usd"]
    assert not getattr(kit, "run_budget_usd", 0)
    assert kit.run_budget_s == RL.UNATTENDED_DEFAULTS["seconds"]
    assert kit.run_budget_turns == RL.UNATTENDED_DEFAULTS["turns"]
    assert "$" not in " ".join(notes)


def test_a_batch_job_keeps_its_own_time_limit(monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "123456")
    eng = _engine("kit.glm-4.6", "kit")
    with _settings():
        notes = RL.apply_unattended(eng)
    assert not getattr(eng, "run_budget_s", 0)
    assert eng.run_budget_turns == RL.UNATTENDED_DEFAULTS["turns"]
    assert "batch job" in " ".join(notes)


def _ctl(ws, notes):
    return am.AutoController(
        ws, session_id="sess-long-1", run_check=lambda e, c: (1, ""),
        notify=lambda text, kind="": notes.append((kind, text)),
        settings={"idle_seconds": 60, "wall": 24 * 3600,
                  "max_iterations": 300, "tokens": 5_000_000})


@pytest.fixture()
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    w = tmp_path / "ws"
    w.mkdir()
    return w


def test_auto_mode_widens_the_defaulted_ceilings_to_its_own_budget(ws):
    eng = _engine()
    with _settings():
        RL.apply_unattended(eng)
        line = _ctl(ws, []).enable(eng)
        usd, secs = eng._run_budget()
        turns = eng._run_budget_turns()
    assert secs >= 24 * 3600 and turns >= 300 and usd == 0.0
    assert "ceiling" in line and "auto mode" in line


def test_a_ceiling_the_user_set_is_not_widened(ws):
    eng = _engine()
    eng.run_budget_s = 3600.0                     # --max-run-seconds 3600
    with _settings():
        RL.apply_unattended(eng)
        _ctl(ws, []).enable(eng)
        assert eng._run_budget()[1] == 3600.0


def test_a_session_ceiling_stops_auto_mode_with_the_resume_command(ws):
    eng = _engine()
    eng.run_budget_s = 3600.0
    notes = []
    ctl = _ctl(ws, notes)
    with _settings():
        RL.apply_unattended(eng)
        ctl.enable(eng)
        eng._run_started_at = time.time() - 4000
        ctl.last_turn_end = 0
        assert ctl.next_prompt(eng, work=[]) in ("", None)
    assert not ctl.enabled
    report = ctl.stop_report
    assert "agent.run_budget_s" in report and "--session sess-long-1" in report
    assert "--max-run-seconds" in report
    assert any(kind == "run_finished" and "run_budget_s" in text
               for kind, text in notes)


def test_a_ceiling_met_inside_a_turn_raises_one_attention_item():
    eng = _engine()
    eng.run_budget_s = 10.0
    eng._run_started_at = time.time() - 20
    seen = []
    with _settings(), patch("delfin.agent.attention.emit_attention",
                            lambda kind, **kw: seen.append((kind, kw))):
        assert eng._run_ceiling_hit(0.0) == "run_budget_s"
        assert eng._run_ceiling_hit(0.0) == "run_budget_s"
    assert len(seen) == 1 and seen[0][0] == "budget_warning"
    detail = seen[0][1]["detail"]
    assert "--session sess-long-1" in detail and "--max-run-seconds" in detail


def test_the_turn_ceiling_has_a_flag():
    args = agent_cli.build_parser().parse_args(["chat", "--max-run-turns", "900"])
    eng = _engine()
    agent_cli._apply_run_budget(eng, args)
    assert eng.run_budget_turns == 900
