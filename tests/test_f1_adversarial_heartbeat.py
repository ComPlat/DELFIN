"""Adversarial contract tests for Target 2 (durable host) — reviewer s2.

The brief: a supervisor run as a SLURM job keeps a heartbeat file on the
shared filesystem and a SECOND supervisor for the same user refuses to
start while the heartbeat is fresh. QS lens: freshness decided on the
heartbeat record's clock, never wall-clock sleeps.

On main 04bf9703b no supervisor heartbeat exists at all (no
delfin/agent/service/durable_*.py; grep for a supervisor heartbeat
module: none). These tests define the contract the builder's module
must satisfy; they are RED here and turn green with the Target-2 hash.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest


def _fresh_record(ts: float) -> dict:
    return {"pid": 4242, "ts": ts, "host": "node-1",
            "sessions": ["sess-a", "sess-b"]}


class TestHeartbeatContract:
    """Contract for a module providing supervisor heartbeat state."""

    def _module(self):
        """Import the builder's heartbeat surface; skip-free probe.

        Returns the module or None. Red state on main: there is none.
        """
        try:
            from delfin.agent.service import durable_host
            return durable_host
        except ImportError:
            return None

    def test_heartbeat_module_exists_on_branch(self):
        assert self._module() is not None, (
            "Target 2: no heartbeat module yet (durable_host) — the "
            "second-supervisor refusal cannot exist")

    def _hb(self, tmp_path, ts: float) -> Path:
        # The heartbeat lives on the SHARED filesystem (HOME/.delfin), not
        # in a per-process rundir: a second supervisor on another node must
        # see it. Contract: resolved per call from the state redirect.
        from delfin.agent import session_store
        p = Path(session_store._SESSIONS_DIR) / "supervisor_heartbeat.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(_fresh_record(ts)), encoding="utf-8")
        return p

    def test_fresh_heartbeat_refuses_a_second_supervisor(self, tmp_path):
        mod = self._module()
        if mod is None:
            pytest.skip("no heartbeat module on this branch (red state)")
        self._hb(tmp_path, time.time())
        assert mod.heartbeat_fresh(str(tmp_path)) is True

    def test_stale_heartbeat_allows_the_start(self, tmp_path):
        mod = self._module()
        if mod is None:
            pytest.skip("no heartbeat module on this branch (red state)")
        self._hb(tmp_path, time.time() - 7200)
        assert mod.heartbeat_fresh(str(tmp_path)) is False

    def test_corrupt_heartbeat_is_not_fresh(self, tmp_path):
        mod = self._module()
        if mod is None:
            pytest.skip("no heartbeat module on this branch (red state)")
        from delfin.agent import session_store
        p = Path(session_store._SESSIONS_DIR) / "supervisor_heartbeat.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("{not json", encoding="utf-8")
        assert mod.heartbeat_fresh(str(tmp_path)) is False

    def test_absent_heartbeat_is_not_fresh(self, tmp_path):
        mod = self._module()
        if mod is None:
            pytest.skip("no heartbeat module on this branch (red state)")
        assert mod.heartbeat_fresh(str(tmp_path)) is False

    def test_freshness_is_decided_on_the_record_clock_not_sleep(self, tmp_path):
        """QS lens: the fresh/stale boundary moves with the record's ts,
        with no sleep in the test."""
        mod = self._module()
        if mod is None:
            pytest.skip("no heartbeat module on this branch (red state)")
        boundary = time.time() - mod.FRESH_WINDOW_S
        self._hb(tmp_path, boundary + 1.0)
        assert mod.heartbeat_fresh(str(tmp_path)) is True
        self._hb(tmp_path, boundary - 1.0)
        assert mod.heartbeat_fresh(str(tmp_path)) is False
