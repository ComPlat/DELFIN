"""SP1 phase 3 red control: the CLI entry point.

Red on the current tree: delfin/agent/explore_pilot.py has no
main()/argparse block yet, so `python -m delfin.agent.explore_pilot`
does nothing. All cases run the CLI through run_pilot-path fakes — no
model call, no real engine.
"""
import subprocess
import sys
import types

from delfin.agent import explore_pilot as ep

_MOD = "delfin.agent.explore_pilot"


def _run_cli(*args):
    """Run the module CLI in a subprocess (the real entry, not an inner
    function); stdout carries the report JSON."""
    return subprocess.run(
        [sys.executable, "-m", _MOD, *args],
        capture_output=True, text=True, timeout=120)


def test_cli_dry_run_runs_both_arms_on_fakes():
    r = _run_cli("--model", "fake", "--dry-run")
    assert r.returncode == 0, r.stderr[-800:]
    assert '"explore_wins"' in r.stdout or "explore_wins" in r.stdout


def test_cli_dry_run_needs_no_model_key_or_network():
    """--dry-run must complete with fakes: exit 0 even though no real
    model endpoint exists in the test environment."""
    r = _run_cli("--model", "fake", "--dry-run", "--tasks",
                 "delfin/agent/pack/benchmark/pilot_tasks_sp1.yaml")
    assert r.returncode == 0, r.stderr[-800:]


def test_cli_model_is_required():
    r = _run_cli("--dry-run")
    assert r.returncode != 0


def test_cli_without_dry_run_needs_a_tasks_file():
    r = _run_cli("--model", "fake")
    assert r.returncode != 0, \
        "a live run without an explicit task file must be refused, " \
        "not silently run the production benchmark"


def test_cli_live_path_uses_explicit_tasks_file(tmp_path, monkeypatch):
    """The live path (--dry-run absent) goes through run_pilot with the
    adapter built from the named tasks file — verified with fakes via
    the module functions the CLI calls, in-process."""
    calls = []

    def fake_run_pilot(tasks, *, runner, root="."):
        calls.append([t.task_id for t in tasks])
        return [ep.TaskRow(task_id=t.task_id) for t in tasks]

    monkeypatch.setattr(ep, "run_pilot", fake_run_pilot)
    monkeypatch.setattr(ep, "make_pilot_runner", lambda **kw: lambda t, arm="plain": {})
    rc = ep.main(["--model", "fake", "--tasks",
                  "delfin/agent/pack/benchmark/pilot_tasks_sp1.yaml"])
    assert rc == 0
    assert calls and len(calls[0]) == 6
