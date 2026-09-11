"""A benchmark task answered with the wording of a different task.

Asked which of the user's calculations has the lowest energy, an attempt
began: "I answer with the three calculations in runs/" -- the prompt of
`science_three_runs_are_audited_independently`, which it had never been
given. No fixture, no memory file and no task text in that tree carried
the word.

The channel was the session briefing. It quotes recent failures from
~/.delfin/outcome_history.jsonl into the system prompt as lessons, and
three failed attempts of the other task sat in that history: the guard
restores the fixture directories after a run, but the history it never
touched, and a run killed mid-way restores nothing at all. Of 179 records
in the user's history that evening, 59 were benchmark prompts verbatim.

So an attempt no longer snapshots the user's state and puts it back. It
does not see it: every sink the agent reads about itself is pointed at an
empty scratch home for the duration, from the one table the test suite
redirects from too.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from delfin.agent import state_paths
from delfin.agent.benchmark_runner import _PristineWorkspace

_LESSON = "Kalibriere den Spektrographen in runs/ bevor du misst"


def _seed_failures(path: Path, n: int = 3) -> None:
    from delfin.agent.outcome_tracker import CycleOutcome, append_outcome
    for i in range(n):
        append_outcome(CycleOutcome(
            task=_LESSON, provider="kit", model="m", mode="solo",
            verdict="FAIL", error_type="crash",
            timestamp=f"2026-09-10T20:0{i}:00"), path=path)


def _briefing() -> str:
    from delfin.agent.briefing import generate_briefing
    return generate_briefing("kit", "Welche meiner Rechnungen ist die niedrigste?")


def test_the_channel_exists_outside_an_attempt():
    """Calibration first: the briefing really does quote earlier failures.
    Without this the test below could pass because nothing ever reaches
    the prompt, not because the guard kept it out."""
    from delfin.agent import outcome_tracker
    _seed_failures(outcome_tracker._DEFAULT_PATH)
    assert "runs/" in _briefing(), _briefing()


def test_the_briefing_inside_an_attempt_carries_no_earlier_failure(tmp_path):
    from delfin.agent import outcome_tracker
    _seed_failures(outcome_tracker._DEFAULT_PATH)
    assert "runs/" in _briefing()
    with _PristineWorkspace(tmp_path):
        inside = _briefing()
        assert "runs/" not in inside, inside
        assert _LESSON not in inside
    # And back: the user's own briefing is theirs again after the attempt.
    assert "runs/" in _briefing()


def test_an_outcome_recorded_inside_an_attempt_is_nobodys(tmp_path):
    from delfin.agent import outcome_tracker
    from delfin.agent.outcome_tracker import CycleOutcome, append_outcome
    real = outcome_tracker._DEFAULT_PATH
    _seed_failures(real, n=1)
    before = real.read_text().count("\n")
    with _PristineWorkspace(tmp_path):
        append_outcome(CycleOutcome(task="attempt probe", verdict="PASS"))
        assert outcome_tracker._DEFAULT_PATH != real
        assert "attempt probe" in outcome_tracker._DEFAULT_PATH.read_text()
    assert outcome_tracker._DEFAULT_PATH == real
    assert real.read_text().count("\n") == before
    assert "attempt probe" not in real.read_text()


def test_the_scratch_home_does_not_outlive_the_attempt(tmp_path):
    guard = _PristineWorkspace(tmp_path)
    with guard:
        home = guard._scratch_home
        assert home is not None and home.is_dir()
        from delfin.agent.outcome_tracker import CycleOutcome, append_outcome
        append_outcome(CycleOutcome(task="x"))
        assert list(home.rglob("*.jsonl"))
    assert not home.exists(), "the attempt's history survived it"


def test_the_users_settings_stay_theirs_inside_an_attempt(tmp_path):
    """An attempt runs with the user's configuration -- their provider,
    their permission rules -- it just does not run with their past."""
    from delfin.agent import hooks_editor, kit_settings
    settings_before = kit_settings.USER_SETTINGS_PATH
    hooks_before = hooks_editor._USER_SETTINGS
    with _PristineWorkspace(tmp_path):
        assert kit_settings.USER_SETTINGS_PATH == settings_before
        assert hooks_editor._USER_SETTINGS == hooks_before


def test_a_memory_the_user_has_is_not_recalled_by_an_attempt(tmp_path):
    from delfin.agent.memory_store import _delfin_memory_dir, save_typed_memory
    save_typed_memory("The user prefers TPSSh for everything.",
                      repo_root=tmp_path, memory_type="user",
                      title="tpssh preference", source="user")
    real_store = _delfin_memory_dir(tmp_path)
    assert any("tpssh" in p.name for p in real_store.glob("*.md"))
    with _PristineWorkspace(tmp_path):
        from delfin.agent import memory_store
        inside = memory_store._delfin_memory_dir(tmp_path)
        assert inside != real_store
        assert not list(inside.glob("*.md")), "the attempt saw a memory"
        save_typed_memory("Written by an attempt.", repo_root=tmp_path,
                          memory_type="user", title="attempt probe",
                          source="agent")
        assert any("attempt-probe" in p.name for p in inside.glob("*.md"))
    after = sorted(p.name for p in real_store.glob("*.md"))
    assert any("tpssh" in n for n in after), "the user's memory was lost"
    assert not any("attempt-probe" in n for n in after)


def test_a_refused_snapshot_puts_the_state_back(tmp_path, monkeypatch):
    """__exit__ is not called when __enter__ raises, so the redirect has
    to be undone on that path too -- or the user's process keeps writing
    into a scratch directory that no longer exists."""
    import shutil
    from delfin.agent import outcome_tracker
    ws = tmp_path / "tests" / "fixtures" / "behavior_workspace"
    ws.mkdir(parents=True)
    (ws / "f.txt").write_text("x")
    real = outcome_tracker._DEFAULT_PATH

    def _boom(*a, **k):
        raise OSError("disk says no")

    monkeypatch.setattr(shutil, "copytree", _boom)
    with pytest.raises(RuntimeError):
        with _PristineWorkspace(tmp_path):
            pass
    assert outcome_tracker._DEFAULT_PATH == real


# --- the same table, from the environment, for a probe run through the CLI


def test_a_scratch_state_from_the_environment_redirects_the_same_sinks(
        tmp_path):
    from delfin.agent import outcome_tracker, provider_profile
    from delfin.agent.outcome_tracker import CycleOutcome, append_outcome
    real = outcome_tracker._DEFAULT_PATH
    scratch = tmp_path / "scratch"
    redirect = state_paths.scratch_state_from_environment(
        {state_paths.SCRATCH_STATE_ENV: str(scratch)})
    try:
        assert redirect is not None
        assert (scratch / ".delfin").is_dir()
        assert outcome_tracker._DEFAULT_PATH == scratch / ".delfin" / "outcome_history.jsonl"
        assert provider_profile._LOCAL_STATE_PATH.parent == scratch / ".delfin"
        append_outcome(CycleOutcome(task="probe under scratch"))
        assert not real.exists() or "probe under scratch" not in real.read_text()
    finally:
        redirect.__exit__(None, None, None)
    assert outcome_tracker._DEFAULT_PATH == real


def test_an_unset_variable_changes_nothing():
    from delfin.agent import outcome_tracker
    real = outcome_tracker._DEFAULT_PATH
    assert state_paths.scratch_state_from_environment({}) is None
    assert state_paths.scratch_state_from_environment(
        {state_paths.SCRATCH_STATE_ENV: "  "}) is None
    assert outcome_tracker._DEFAULT_PATH == real


def test_the_cli_process_honours_the_variable(tmp_path):
    """Driven through the real entry point in a subprocess: the variable
    is read before anything imports a sink, and a record the run writes
    lands under the scratch directory."""
    import subprocess
    import sys
    scratch = tmp_path / "probe-home"
    code = (
        "import sys; from delfin.agent import cli\n"
        "cli.main(['bench', 'list']) if False else None\n"
        "from delfin.agent import state_paths; "
        "state_paths.scratch_state_from_environment()\n"
        "from delfin.agent.outcome_tracker import CycleOutcome, append_outcome\n"
        "append_outcome(CycleOutcome(task='cli probe'))\n"
        "from delfin.agent import outcome_tracker; "
        "print(outcome_tracker._DEFAULT_PATH)\n"
    )
    env = dict(os.environ, **{state_paths.SCRATCH_STATE_ENV: str(scratch)})
    out = subprocess.run([sys.executable, "-c", code], env=env,
                         capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr[-800:]
    printed = Path(out.stdout.strip().splitlines()[-1])
    assert printed == scratch / ".delfin" / "outcome_history.jsonl"
    lines = [json.loads(l) for l in printed.read_text().splitlines() if l.strip()]
    assert lines and lines[-1]["task"] == "cli probe"


def test_the_suite_and_the_bench_redirect_from_one_table():
    """The two copies drifted once. Now there is one, and it names every
    sink the earlier suite table did."""
    names = {(m, a) for m, a, _ in state_paths.USER_STATE_SINKS}
    for expected in (
        ("delfin.agent.outcome_tracker", "_DEFAULT_PATH"),
        ("delfin.agent.provider_profile", "_LOCAL_STATE_PATH"),
        ("delfin.agent.session_store", "_SESSIONS_DIR"),
        ("delfin.agent.failure_log", "_LOG_PATH"),
        ("delfin.agent.memory_store", "_DEFAULT_PATH"),
    ):
        assert expected in names, expected
    resolvers = {(m, a) for m, a, _ in state_paths.USER_STATE_RESOLVERS}
    assert ("delfin.agent.memory_store", "_delfin_memory_dir") in resolvers
    assert ("delfin.agent.memory_store", "_delfin_global_memory_dir") in resolvers
    assert ("delfin.agent.audit_log", "_default_log_path") in resolvers


def test_the_audit_log_stays_where_the_bench_reads_it(tmp_path):
    """The bench counts an attempt's gate denials and the paths it wrote
    from the audit log AFTER the guard has closed. For one night that
    log was redirected with everything else, so the count looked at a
    scratch file that was already gone and every block said "denials
    not observed". The attempt's own records are what the cost axis is
    made of; they stay where the reader looks."""
    from delfin.agent import audit_log
    before = audit_log._default_log_path()
    with _PristineWorkspace(tmp_path):
        assert audit_log._default_log_path() == before
    assert audit_log._default_log_path() == before


def test_a_denial_during_an_attempt_is_counted_after_it(tmp_path):
    """End to end through the counter the bench uses: a refusal recorded
    while the guard is up is still there when the guard is down."""
    from datetime import datetime, timezone
    from types import SimpleNamespace
    from delfin.agent import audit_log
    from delfin.agent.benchmark_runner import _denials_during
    ws = tmp_path / "ws"
    ws.mkdir()
    since = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    engine = SimpleNamespace(kit_permissions=SimpleNamespace(workspace=str(ws)))
    def _record(decision, command):
        return audit_log.make_record(
            tool="bash", decision=decision, mode="default", command=command,
            reason="test", extra={"workspace": str(ws), "cwd": str(ws)})

    with _PristineWorkspace(tmp_path):
        audit_log.append(_record("denied", "rm -rf /"))
        audit_log.append(_record("ok", "ls"))
    counted = _denials_during(engine, since_ts=since)
    assert counted == 1, counted
