"""OT3 adversarial tests (reviewer nacht-s2) for delfin-agent test-on-slurm.

These attack the RENDERED script and the CLI wiring, not the inner helpers:
every test goes through build_job / cmd_test_on_slurm, the functions the
product calls. They must hold on the rendered text with no scheduler, no
network and no live SLURM involved.

Checks:
1. No interpreter placeholder survives rendering when the CLI is used
   normally (the CLI wiring cli.py:6029-6046 passes no --python flag, so
   sys.executable is the only sane default).
2. OUTPUT_DIR is absolute and equals the run dir passed to submit, so logs
   survive the ``cd "$LOCAL"`` on the compute node.
3. A run dir on a login-node-only path is refused with a clear message.
4. The SBATCH header carries --export=ALL exactly: no explicit variable
   exports, no --export=ALL,VARS spelling.
"""

import sys

from delfin.agent import slurm_tests


def _render(**kw):
    args = dict(repo="/repo", ref="abc123",
                test_paths=("tests/test_x.py",),
                partition="cpu", minutes=10)
    args.update(kw)
    return slurm_tests.build_job(**args)


def _run_cli(monkeypatch, run_dir_arg, tmp_path):
    """Drive the REAL entry point cmd_test_on_slurm with submit() captured.

    Returns (captured job_text, captured run_dir). No scheduler is touched.
    """
    captured = {}

    def fake_submit(text, run_dir):
        captured["text"] = text
        captured["run_dir"] = run_dir
        return "12345"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)

    class Args:
        repo = "/repo"
        ref = "abc123"
        partition = "cpu"
        minutes = 10
        tests = ["tests/test_x.py"]
        run_dir = run_dir_arg

    rc = slurm_tests.cmd_test_on_slurm(Args())
    assert rc == 0, rc
    assert "text" in captured and "run_dir" in captured
    return captured["text"], captured["run_dir"]


def test_cli_default_render_has_no_placeholder():
    # The CLI (delfin/agent/cli.py:6029-6046) never passes an interpreter;
    # the rendered job must still run. sys.executable is the documented
    # fallback the operator's finding demands.
    text = _render(python=None)
    assert "PYTHON_PLACEHOLDER" not in text, (
        "rendered script still carries the interpreter placeholder")


def test_cli_default_render_runs_sys_executable():
    text = _render(python=None)
    assert sys.executable in text


def test_cli_wiring_output_dir_matches_run_dir(tmp_path, monkeypatch):
    # The REAL call path: cmd_test_on_slurm renders the job AND hands the
    # run dir to submit. The rendered OUTPUT_DIR must be exactly the run dir
    # submit gets, absolute, or the logs land node-local after cd "$LOCAL"
    # (slurm_tests.py:231-235) and are lost.
    run_dir = str(tmp_path / "ot3run")
    text, got_run_dir = _run_cli(monkeypatch, run_dir, tmp_path)
    assert got_run_dir == run_dir
    assert 'OUTPUT_DIR="' + run_dir + '"' in text, (
        "rendered OUTPUT_DIR does not equal the run dir handed to submit")


def test_cli_wiring_output_dir_not_relative_logs(tmp_path, monkeypatch):
    # The product must never render the relative default "logs": on the
    # compute node it resolves against $LOCAL and the logs vanish.
    run_dir = str(tmp_path / "ot3run2")
    text, _ = _run_cli(monkeypatch, run_dir, tmp_path)
    assert 'OUTPUT_DIR="logs"' not in text


def test_run_dir_on_local_only_path_refused():
    # A path like /scratch (login-node scratchpad) is invisible to the
    # compute node; build_job must refuse it before submission.
    try:
        _render(output_dir="/scratch/someone/ot3")
    except Exception as exc:
        assert "/scratch" in str(exc) or "shared" in str(exc).lower(), exc
    else:
        raise AssertionError("local-only run dir accepted, logs will be lost")


def test_export_is_exactly_all():
    text = _render()
    for line in text.splitlines():
        if line.startswith("#SBATCH --export"):
            assert line == "#SBATCH --export=ALL", line


def test_no_explicit_sensitive_exports():
    text = _render()
    for line in text.splitlines():
        if "export" in line and line.startswith("export "):
            # only the documented PYTHONNOUSERSITE / PYTHONHASHSEED unset-
            # then-set block is allowed to export
            assert "PYTHONNOUSERSITE" in line or "PYTHONHASHSEED" in line, line
