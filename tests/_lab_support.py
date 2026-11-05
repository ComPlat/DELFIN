"""Shared set-up for the lab-control tests: a scratch lab, the real executor.

Every test gets its own HOME, lab directory and person, the demo manifest
(two simulated modules) written and pinned, and no agent-session marker.
Tool calls go through ``api_client._DocToolExecutor.execute`` -- the path a
model's call takes -- never straight into the kernel.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


class FakeClock:
    def __init__(self, t0: float):
        self.t = float(t0)

    def __call__(self) -> float:
        return self.t

    def sleep(self, s: float) -> None:
        self.t += float(s)


@pytest.fixture
def lab(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "lab"))
    monkeypatch.setenv("DELFIN_PERSON", "alice")
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.delenv("DELFIN_LAB_MANIFEST", raising=False)
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()
    from delfin.lab import manifest, smoke, store
    from delfin.lab.safety import pin
    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(smoke.demo_manifest_text(), encoding="utf-8")
    pin(manifest.load(path), by="alice")
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def call(ws: Path, tool: str, args: dict, *, mode: str = "bypassPermissions",
         session: str = "s1") -> dict:
    from delfin.agent import api_client as A
    perms = A.KitToolPermissions(workspace=ws)
    perms.mode = mode
    perms.confirm_callback = lambda *a, **k: True
    perms.task_session_id = session
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    return json.loads(executor.execute(tool, args, perms))


def sim(device: str) -> dict:
    from delfin.lab.drivers import simulator
    return simulator.raw(device).get("values", {})


def prepare_heating(ws, session="s1"):
    """Stirrer, hood and hotplate power on: what heating above 40 degC needs."""
    for dev in ("stirrer", "fume_hood", "hotplate"):
        out = call(ws, "lab_set", {"device": dev, "channel": "power",
                                   "value": "on"}, session=session)
        assert out.get("ok"), out


def approve_and_set(ws, args, session="s1", person=None):
    from delfin.lab import approvals
    first = call(ws, "lab_set", args, session=session)
    assert first.get("code") == "approval_required", first
    approvals.grant(first["approval_request"], approver=person)
    return call(ws, "lab_set", {**args, "approval": first["approval_request"]},
                session=session)
