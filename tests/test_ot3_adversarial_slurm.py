"""OT3 adversarial tests (reviewer nacht-s2) for delfin-agent test-on-slurm.

These attack the product call path, not inner helpers: every B-finding test
drives cmd_test_on_slurm — the function the CLI actually calls — with
submit() captured, so no scheduler, no network and no live SLURM is touched.
Run dirs sit under a fake repository directory, one of the two locations
_node_visible() accepts as shared.

Checks:
1. The rendered job carries the session interpreter, never
   PYTHON_PLACEHOLDER (the CLI wiring cli.py:6029-6046 passes no --python
   flag, so sys.executable is the only sane default).
2. OUTPUT_DIR is absolute and equals the run dir handed to submit, so the
   logs survive the ``cd "$LOCAL"`` on the compute node.
3. A run dir on a login-node-only path (/scratch) is refused with a clear
   message BEFORE anything is submitted.
4. The SBATCH header carries --export=ALL exactly: no explicit variable
   exports, no --export=ALL,VARS spelling.
"""

import sys
from pathlib import Path

import pytest

from delfin.agent import slurm_tests


def _render(**kw):
    args = dict(repo="/repo", ref="abc123",
                test_paths=("tests/test_x.py",),
                partition="cpu", minutes=10)
    args.update(kw)
    return slurm_tests.build_job(**args)


def _run_cli(monkeypatch, run_dir_arg, repo):
    """Drive the REAL entry point cmd_test_on_slurm with submit() captured.

    Returns (captured job_text, captured run_dir). No scheduler is touched.
    """
    captured = {}

    def fake_submit(text, run_dir):
        captured["text"] = text
        captured["run_dir"] = run_dir
        return "12345"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)

    import types
    args = types.SimpleNamespace(repo=repo, ref="abc123", partition="cpu",
                                 minutes=10, tests=["tests/test_x.py"],
                                 run_dir=run_dir_arg)
    rc = slurm_tests.cmd_test_on_slurm(args)
    assert rc == 0, rc
    assert "text" in captured and "run_dir" in captured
    return captured["text"], captured["run_dir"]


def test_cli_default_render_has_no_placeholder(tmp_path, monkeypatch):
    # The CLI never passes an interpreter; the rendered job must still run.
    repo = tmp_path / "repo"
    run_dir = str(repo / "ot3run")
    text, _ = _run_cli(monkeypatch, run_dir, str(repo))
    assert "PYTHON_PLACEHOLDER" not in text, (
        "rendered script still carries the interpreter placeholder")


def test_cli_default_render_runs_sys_executable(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    run_dir = str(repo / "ot3run")
    text, _ = _run_cli(monkeypatch, run_dir, str(repo))
    assert sys.executable in text


def test_cli_wiring_output_dir_matches_run_dir(tmp_path, monkeypatch):
    # cmd_test_on_slurm renders the job AND hands the run dir to submit.
    # The rendered OUTPUT_DIR must be exactly that run dir, absolute, or the
    # logs land node-local after cd "$LOCAL" and are lost.
    repo = tmp_path / "repo"
    run_dir = str(repo / "ot3run")
    text, got_run_dir = _run_cli(monkeypatch, run_dir, str(repo))
    resolved = str(Path(run_dir).resolve())
    # submit gets the run dir (as Path or str — same directory either way)
    assert str(Path(got_run_dir)) == resolved
    assert 'OUTPUT_DIR="' + resolved + '"' in text, (
        "rendered OUTPUT_DIR does not equal the run dir handed to submit")


def test_cli_wiring_output_dir_not_relative_logs(tmp_path, monkeypatch):
    # The product must never render the relative default "logs": on the
    # compute node it resolves against $LOCAL and the logs vanish.
    repo = tmp_path / "repo"
    run_dir = str(repo / "ot3run2")
    text, _ = _run_cli(monkeypatch, run_dir, str(repo))
    assert 'OUTPUT_DIR="logs"' not in text
    assert 'OUTPUT_DIR="' + str(Path(run_dir).resolve()) + '"' in text


def test_run_dir_on_local_only_path_refused(tmp_path, monkeypatch, capsys):
    # A path like /scratch (login-node scratchpad) is invisible to the
    # compute node; the command must refuse it BEFORE any submission.
    submitted = []

    def fake_submit(text, run_dir):
        submitted.append((text, run_dir))
        return "99999"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)

    import types
    args = types.SimpleNamespace(
        repo=str(tmp_path / "repo"),
        ref="abc123",
        partition="cpu",
        minutes=10,
        tests=["tests/test_x.py"],
        run_dir="/scratch/someone/ot3",
    )

    # Final contract (QS ruling, builder 8b85829a6): the refusal is
    # fail-closed in the product — cmd_test_on_slurm RAISES ValueError and
    # submit is never reached — while the promised exit 2 lives at the CLI
    # wrapper (cli.py:6029-6039), which catches and prints. Assert both
    # layers of the contract.
    with pytest.raises(ValueError) as raised:
        slurm_tests.cmd_test_on_slurm(args)
    assert "/scratch/someone/ot3" in str(raised.value)
    assert "shared" in str(raised.value)
    assert submitted == [], "submit was called despite the refusal"


def test_cli_wrapper_returns_2_on_refusal(tmp_path, monkeypatch, capsys):
    # The CLI layer (cli.py:6029-6039) must turn the fail-closed ValueError
    # into the promised exit 2 with the message on stdout — no traceback.
    # Driven through the REAL wiring: the handler cli.py actually registers
    # (t6.set_defaults), not cmd_test_on_slurm directly.
    submitted = []

    def fake_submit(text, run_dir):
        submitted.append((text, run_dir))
        return "99999"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)

    import types
    args = types.SimpleNamespace(
        repo=str(tmp_path / "repo"),
        ref="abc123",
        partition="cpu",
        minutes=10,
        tests=["tests/test_x.py"],
        run_dir="/scratch/someone/ot3",
    )
    # The registered handler is a closure inside the parser builder; get it
    # the way the product does: build the parser and read set func off the
    # parsed args namespace.
    from delfin.agent import cli as cli_mod
    parser = cli_mod.build_parser()
    parsed = parser.parse_args(
        ["test-on-slurm", str(tmp_path / "repo"), "abc123", "cpu", "10",
         "tests/test_x.py", "--run-dir", "/scratch/someone/ot3"])
    rc = parsed.func(parsed)
    assert rc == 2, rc
    assert submitted == [], "CLI wrapper submitted despite the refusal"
    captured = capsys.readouterr()
    assert "/scratch/someone/ot3" in captured.out
    assert "shared" in captured.out


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
