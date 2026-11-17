"""Reviewer s2 adversarial battery: `delfin-agent durable` CLI wiring
(builder commit e96154c92, durable_host.py cmd_durable / render / heartbeat).

Attacks: SLURM-directive injection through partition/qos, a corrupt
heartbeat file, and a heartbeat whose ts lies in the future (clock skew) —
a future ts keeps heartbeat_fresh() True for as long as the skew lasts,
which is exactly the second-supervisor refusal Target 2 must not allow to
be spoofed or wedged.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent.cli import build_parser
from delfin.agent.service import durable_host


def _run(args, monkeypatch):
    """One CLI invocation through build_parser — the real public path."""
    import delfin.agent.session_store as session_store

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", "/nonexistent-s2")
    parser = build_parser()
    ns = parser.parse_args(args)
    return ns.func(ns)


class TestRenderInjection:
    """Every #SBATCH value must survive the allow-list, not just walltime."""

    def test_partition_newline_refused(self, monkeypatch, capsys):
        rc = _run(["durable", "render",
                   "--partition", "standard\n#SBATCH --nodelist=evil"],
                  monkeypatch)
        assert rc != 0
        assert "nodelist" not in capsys.readouterr().out

    def test_qos_attached_value_refused(self, monkeypatch, capsys):
        rc = _run(["durable", "render", "--qos", "normal --exclusive"],
                  monkeypatch)
        assert rc != 0
        assert "--exclusive" not in capsys.readouterr().out

    def test_signal_minutes_out_of_range_refused(self, monkeypatch, capsys):
        rc = _run(["durable", "render", "--signal-minutes", "600"],
                  monkeypatch)
        assert rc != 0
        assert "signal=" not in capsys.readouterr().out


class TestStatusRobustness:
    """status must degrade to a plain message, never crash."""

    def test_corrupt_heartbeat_is_plain(self, monkeypatch, tmp_path, capsys):
        import delfin.agent.session_store as session_store

        monkeypatch.setattr(session_store, "_SESSIONS_DIR", str(tmp_path))
        (tmp_path / durable_host.HEARTBEAT_FILE).write_text(
            "{not json at all", encoding="utf-8")
        rc = _run(["durable", "status"], monkeypatch)
        assert rc == 0
        out = capsys.readouterr().out
        assert "no supervisor heartbeat" in out

    def test_future_ts_is_stale(self, monkeypatch, tmp_path, capsys):
        """Clock-skew attack: ts one hour in the future must NOT read as
        fresh — otherwise a dead-but-ahead host wedges the second
        supervisor for the whole skew (Target 2)."""
        import delfin.agent.session_store as session_store

        monkeypatch.setattr(session_store, "_SESSIONS_DIR", str(tmp_path))
        durable_host.write_supervisor_heartbeat(
            ["a"], pid=1, host="h", root=str(tmp_path))
        # Force the written record's ts one hour into the future.
        p = tmp_path / durable_host.HEARTBEAT_FILE
        rec = json.loads(p.read_text(encoding="utf-8"))
        rec["ts"] = time.time() + 3600
        p.write_text(json.dumps(rec), encoding="utf-8")
        assert durable_host.heartbeat_fresh() is False
