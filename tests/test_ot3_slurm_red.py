"""OT3 Part B: red tests for the test-on-slurm defects found 01.11.

All four assertions test the RENDERED script / the real CLI wiring, never
the scheduler:

* B1: the rendered script must contain no PYTHON_PLACEHOLDER and must
  reference the session's own interpreter (sys.executable by default).
  RED before the fix: build_job renders the literal placeholder and
  cmd_test_on_slurm never passes a python path — the job exited 127.
* B2: the rendered OUTPUT_DIR must be the ABSOLUTE run dir the CLI hands
  to submit(), never the relative default "logs" (which resolved against
  the node-local $LOCAL after `cd` and was lost). RED before the fix.
* B3: a run dir on a filesystem the compute node cannot see (e.g. a
  login-node-only scratchpad) must be refused with a clear message BEFORE
  submission. RED before the fix: silently accepted.
* B4: the script pins `#SBATCH --export=ALL` exactly — never an explicit
  variable list (already true; pinned here so it stays true).
"""

import os
import sys

import pytest

from delfin.agent import slurm_tests


@pytest.fixture()
def repo_dir(tmp_path):
    """A fake repository dir (only its path crosses into the script)."""
    repo = tmp_path / "repo"
    (repo / "tests").mkdir(parents=True)
    return repo


def _args(repo, tmp_path, run_dir="."):
    """An argparse.Namespace shaped like the test-on-slurm CLI call."""
    import argparse
    return argparse.Namespace(
        repo=str(repo), ref="main", partition="cpu", minutes=10,
        tests=["tests/test_x.py"], run_dir=run_dir)


@pytest.fixture()
def monkeypatch_home(monkeypatch):
    """Marker fixture: tests that touch HOME use a fake HOME."""


# ---------------------------------------------------------------- B1

def test_cli_default_render_has_no_placeholder(repo_dir, monkeypatch,
                                               tmp_path):
    captured = {}

    def fake_submit(job_text, run_dir, **kw):
        captured["text"] = job_text
        captured["run_dir"] = run_dir
        return "Submitted batch job 4242"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    monkeypatch.setattr(slurm_tests, "_node_visible", lambda p, r: True)
    rc = slurm_tests.cmd_test_on_slurm(_args(repo_dir, tmp_path))
    assert rc == 0
    assert "PYTHON_PLACEHOLDER" not in captured["text"]


def test_cli_default_render_runs_sys_executable(repo_dir, monkeypatch,
                                                tmp_path):
    captured = {}

    def fake_submit(job_text, run_dir, **kw):
        captured["text"] = job_text
        captured["run_dir"] = run_dir
        return "Submitted batch job 4243"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    monkeypatch.setattr(slurm_tests, "_node_visible", lambda p, r: True)
    rc = slurm_tests.cmd_test_on_slurm(_args(repo_dir, tmp_path))
    assert rc == 0
    assert sys.executable in captured["text"]


# ---------------------------------------------------------------- B2

def test_cli_wiring_output_dir_matches_run_dir(repo_dir, monkeypatch,
                                               tmp_path):
    captured = {}

    def fake_submit(job_text, run_dir, **kw):
        captured["text"] = job_text
        captured["run_dir"] = run_dir
        return "Submitted batch job 4244"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    monkeypatch.setattr(slurm_tests, "_node_visible", lambda p, r: True)
    run_dir = tmp_path / "run"
    rc = slurm_tests.cmd_test_on_slurm(
        _args(repo_dir, tmp_path, run_dir=str(run_dir)))
    assert rc == 0
    text = captured["text"]
    # The rendered OUTPUT_DIR is absolute and equals the dir submit ran in.
    assert f'OUTPUT_DIR="{run_dir.resolve()}"' in text
    assert 'OUTPUT_DIR="logs"' not in text


def test_cli_wiring_output_dir_not_relative_logs(repo_dir, monkeypatch,
                                                 tmp_path):
    """Even with the CLI default run-dir (".") the rendered OUTPUT_DIR is
    the run dir's absolute path, never a relative "logs"."""
    captured = {}

    def fake_submit(job_text, run_dir, **kw):
        captured["text"] = job_text
        captured["run_dir"] = run_dir
        return "Submitted batch job 4245"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    monkeypatch.setattr(slurm_tests, "_node_visible", lambda p, r: True)
    rc = slurm_tests.cmd_test_on_slurm(_args(repo_dir, tmp_path))
    assert rc == 0
    assert 'OUTPUT_DIR="logs"' not in captured["text"]
    assert 'OUTPUT_DIR="' + os.path.join(os.getcwd(), "logs") + '"' \
        not in captured["text"]


# ---------------------------------------------------------------- B3

def test_run_dir_on_local_only_path_refused(repo_dir, monkeypatch,
                                            tmp_path, capsys):
    """A run dir on a node-invisible filesystem is refused BEFORE submit.
    Fail-closed contract (QS, phase-3 round): the product function RAISES —
    a caller that ignores rc=2 must never submit a job whose logs the node
    cannot write back. The CLI layer converts this to the promised exit 2
    (asserted by test_cli_refuses_node_local_run_dir_with_exit_2)."""
    called = []

    def fake_submit(job_text, run_dir, **kw):
        called.append(run_dir)
        return "Submitted batch job 4246"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    local_only = "/scratch/only-login-node-sees/run"
    with pytest.raises(ValueError) as raised:
        slurm_tests.cmd_test_on_slurm(
            _args(repo_dir, tmp_path, run_dir=local_only))
    called_msg = f"submit ran {called}"
    assert called == [], f"submit must not run for a node-invisible run dir: {called_msg}"
    assert "shared" in str(raised.value) or local_only in str(raised.value)


def test_cli_refuses_node_local_run_dir_with_exit_2(repo_dir, monkeypatch,
                                                    tmp_path, capsys):
    """CLI promise kept: the ValueError from the product refusal reaches a
    CLI user as the promised exit 2 with the message — no traceback.
    Real call path: build_parser().parse_args([...]) -> set_defaults(func)
    handler (the way test_t6_cli_registered.py drives the subcommand)."""
    cli_mod = pytest.importorskip("delfin.agent.cli")
    args = cli_mod.build_parser().parse_args(
        ["test-on-slurm", str(repo_dir), "main", "cpu", "10",
         "tests/test_x.py", "--run-dir",
         "/scratch/only-login-node-sees/run"])
    rc = args.func(args)
    assert rc == 2, f"expected the promised exit 2, got {rc}"
    msg = capsys.readouterr().out
    assert ("shared" in msg
            or "/scratch/only-login-node-sees/run" in msg)


def test_run_dir_under_repo_is_accepted(repo_dir, monkeypatch, tmp_path):
    """The check does not refuse a legitimate shared location: a run dir
    under the repo (a shared filesystem) submits."""
    called = []

    def fake_submit(job_text, run_dir, **kw):
        called.append(run_dir)
        return "Submitted batch job 4247"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    run_dir = repo_dir / "slurm_runs"
    rc = slurm_tests.cmd_test_on_slurm(
        _args(repo_dir, tmp_path, run_dir=str(run_dir)))
    assert rc == 0
    assert called, "a run dir under the repo must submit"


def test_run_dir_under_home_is_accepted(repo_dir, monkeypatch, tmp_path,
                                       monkeypatch_home):
    monkeypatch.setenv("HOME", str(tmp_path / "fakehome"))
    assert slurm_tests._node_visible(
        tmp_path / "fakehome" / "runs", str(repo_dir)) is True


# ---------------------------------------------------------------- B4

def test_export_is_exactly_all(repo_dir, tmp_path):
    text = slurm_tests.build_job(str(repo_dir), "main",
                                 ["tests/test_x.py"], "cpu", 10)
    assert "#SBATCH --export=ALL" in text


def test_no_explicit_sensitive_exports(repo_dir, tmp_path):
    text = slurm_tests.build_job(str(repo_dir), "main",
                                 ["tests/test_x.py"], "cpu", 10)
    assert "--export=ALL," not in text
