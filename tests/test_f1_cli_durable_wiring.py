"""Wiring of ``delfin-agent durable`` (package F1).

Red control for the CLI consumer of
``delfin.agent.service.durable_host``: before the wiring the module's
``render_durable_job`` and heartbeat helpers have no public caller at all
(grep: only durable_host.py defines them). Every path here goes through
the CLI's own ``build_parser`` -- the real public path -- and never
submits anything to SLURM: rendering and status are offline commands.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.cli import build_parser


def _run(args, monkeypatch, tmp_path):
    """Parse + dispatch one CLI invocation; returns its exit code.

    The session store is redirected per call (state_paths style) so the
    status action reads a throwaway heartbeat, never a real home.
    """
    import delfin.agent.session_store as session_store

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    parser = build_parser()
    ns = parser.parse_args(args)
    return ns.func(ns)


def test_render_prints_a_durable_job_script(monkeypatch, capsys):
    """`durable render` prints the rendered SLURM job to stdout."""
    rc = _run(["durable", "render", "--walltime", "06:00:00"],
              monkeypatch, "/nonexistent-f1-wiring")
    out = capsys.readouterr().out
    assert rc == 0
    assert "#SBATCH --export=ALL" in out
    assert "#SBATCH --time=06:00:00" in out
    assert "-m delfin.agent.service.supervisor" in out


def test_render_refuses_a_bad_walltime(monkeypatch, capsys):
    """An injectable or malformed walltime is refused, not rendered."""
    rc = _run(["durable", "render", "--walltime", "not-a-time\n#SBATCH --x"],
              monkeypatch, "/nonexistent-f1-wiring")
    assert rc != 0
    assert "time=" not in capsys.readouterr().out


def test_status_reports_no_supervisor_when_heartbeat_absent(
        monkeypatch, tmp_path, capsys):
    """`durable status` on a store without a heartbeat says so plainly."""
    rc = _run(["durable", "status"], monkeypatch, tmp_path)
    out = capsys.readouterr().out
    assert rc == 0
    assert "none" in out.lower() or "no supervisor" in out.lower()


def test_status_sees_a_written_heartbeat(monkeypatch, tmp_path, capsys):
    """A heartbeat the real writer wrote is reported fresh with its host.

    Real writer: ``durable_host.write_supervisor_heartbeat`` -- not a dict
    typed by the test.
    """
    from delfin.agent.service import durable_host
    import delfin.agent.session_store as session_store

    # Redirect the store BEFORE the real writer runs: every path (writer
    # and reader) must see the throwaway directory, never the real home.
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", str(tmp_path))
    durable_host.write_supervisor_heartbeat(
        ["alpha", "beta"], pid=4242, host="compute-01", root=str(tmp_path))
    rc = _run(["durable", "status"], monkeypatch, tmp_path)
    out = capsys.readouterr().out
    assert rc == 0
    assert "fresh" in out.lower()
    assert "compute-01" in out
    assert "alpha" in out
