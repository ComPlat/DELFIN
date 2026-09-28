"""A chemistry task needs a workspace that exists before setup runs.

``run_setup`` launches the setup script with ``cwd=workspace``; a
workspace path that does not exist makes Popen fail with a bare
FileNotFoundError before the script's first line - on the cluster this
read as "setup script could not run: [Errno 2] No such file or
directory" while the script itself was fine.

The chemistry class follows the science pattern: a packaged fixture
directory the setup script seeds, registered in ``workspace_for`` so
the runner hands the task a real, existing folder.
"""

from pathlib import Path

from delfin.agent.benchmark_runner import workspace_for


def test_chemistry_class_gets_a_workspace():
    ws = workspace_for(Path("."), task_class="chemistry")
    assert ws is not None, (
        "chemistry has no workspace_for mapping; run_setup is started "
        "with a cwd that does not exist")
    assert ws.is_dir(), f"{ws} is mapped but not packaged"
    # The setup script seeds chem/ inside the workspace, so the fixture
    # README is the thing that tells the model what is there.
    assert (ws / "README.md").exists() or any(ws.iterdir()), (
        "fixture dir is empty - nothing to anchor the task")


def test_setup_runs_in_the_chemistry_workspace():
    import subprocess
    import sys
    from delfin.agent.benchmark_runner import run_setup
    ws = workspace_for(Path("."), task_class="chemistry")
    # A scratch copy: the real fixture must not be dirtied by the test.
    import tempfile
    with tempfile.TemporaryDirectory(dir=str(Path(".").resolve())) as td:
        target = Path(td) / "ws"
        target.mkdir()
        ok, out = run_setup("chem_start_geometries.py", target,
                            root=Path(".").resolve())
        assert ok, f"setup failed: {out}"
        assert (target / "chem" / "acetaminophen" / "start.xyz").is_file()
