"""The lab tools are wired, gated and pinned like every other tool.

* the five schemas are in the catalogue and reach ``delfin.lab.tools``;
* they are advertised only on an installation with a lab manifest;
* the e-stop runs in plan mode, a write does not;
* the security-event kinds the lab records are labelled;
* the person's actions (approve, clear, estop --clear, manifest pin, node
  keygen) are refused to an agent session as command text, as a process,
  and the lab directory is not writable through the write tools (T12);
* the deploy units are hardened and the smoke check prints ok.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from _lab_support import call, lab_fixture, sim  # noqa: F401

REPO = Path(__file__).resolve().parents[1]
LAB_TOOLS = {"lab_devices", "lab_read", "lab_set", "lab_plan_check",
             "lab_estop"}


def test_the_five_tools_are_in_the_catalogue():
    from delfin.agent import api_client as A
    names = {t["function"]["name"] for t in A._DOC_TOOLS_OPENAI}
    assert LAB_TOOLS <= names
    from delfin.lab.tools import TOOL_NAMES
    assert TOOL_NAMES == LAB_TOOLS


def test_a_tool_call_reaches_the_kernel(lab):
    out = call(lab, "lab_devices", {})
    ids = {d["id"] for d in out["devices"]}
    assert {"hotplate", "uvvis", "mains_bypass"} <= ids
    hp = next(d for d in out["devices"] if d["id"] == "hotplate")
    assert hp["writes"]["setpoint"]["level"] == "hazardous"
    out = call(lab, "lab_read", {"device": "hotplate_probe"})
    assert out["readings"]["temperature"]["unit"] == "degC"


def test_advertised_only_where_a_manifest_exists(lab, tmp_path, monkeypatch):
    from delfin.agent.api_client import (ToolSurfaceContext, _lab_configured,
                                         tool_unavailable_reason)
    assert _lab_configured() is True
    ctx = ToolSurfaceContext(has_lab=_lab_configured())
    assert tool_unavailable_reason("lab_set", ctx) is None
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "empty"))
    assert _lab_configured() is False
    ctx = ToolSurfaceContext(has_lab=_lab_configured())
    why = tool_unavailable_reason("lab_set", ctx)
    assert why and "manifest" in why
    assert tool_unavailable_reason("read_file", ctx) is None
    # the live surface passes the real answer (one stat of one file)
    import inspect
    from delfin.agent import api_client as A
    assert "has_lab=_lab_configured()" in inspect.getsource(A)
    out = call(lab, "lab_set", {"device": "lamp", "channel": "power",
                                "value": "on"})
    assert "no lab manifest" in out["error"]


def test_estop_runs_in_plan_mode_and_a_write_does_not(lab):
    out = call(lab, "lab_set", {"device": "lamp", "channel": "power",
                                "value": "on"}, mode="plan")
    assert not out.get("ok") and sim("lamp").get("power", "off") == "off"
    assert call(lab, "lab_plan_check", {"steps": [
        {"device": "lamp", "channel": "power", "value": "on"}]},
        mode="plan")["ok"]
    assert call(lab, "lab_estop", {}, mode="plan")["scope"] == "global"


def test_lab_event_kinds_are_labelled():
    from delfin.agent import security_events
    from delfin.lab.audit import LAB_EVENT_KINDS
    assert LAB_EVENT_KINDS == {
        "lab_refused", "lab_safe_state", "lab_estop", "lab_estop_cleared",
        "lab_deadman", "lab_readback_mismatch", "lab_approval_refused",
        "lab_manifest_refused"}
    assert LAB_EVENT_KINDS <= security_events.known_kinds()


@pytest.mark.parametrize("cmd", [
    "delfin-agent lab approve lab-1a2b3c4d",
    "delfin-agent lab clear hotplate",
    "delfin-agent lab estop --clear",
    "delfin-agent lab manifest pin",
    "delfin-agent lab node keygen edge_synth",
    "python -m delfin.agent.cli lab approve lab-1a2b3c4d",
])
def test_the_agent_cannot_approve_or_clear(cmd):
    from delfin.agent import approval_guard
    assert approval_guard.operator_command(cmd), cmd


@pytest.mark.parametrize("cmd", [
    "delfin-agent lab estop", "delfin-agent lab devices",
    "delfin-agent lab pending", "delfin-agent lab manifest check",
])
def test_reading_and_stopping_stay_allowed(cmd):
    from delfin.agent import approval_guard
    assert approval_guard.operator_command(cmd) is None, cmd


def test_the_shell_tool_refuses_an_approval(lab):
    out = call(lab, "bash", {"command": "delfin-agent lab approve lab-1a2b3c4d"})
    text = str(out)
    assert "blocked" in text or "operator" in text, out


def test_an_agent_process_cannot_approve(lab, monkeypatch):
    prep = call(lab, "lab_set", {"device": "pump", "channel": "flow",
                                 "value": 1, "unit": "mL/min"})
    # pump needs the hood: the request is filed only once the interlock holds
    assert prep["code"] == "interlock"
    assert call(lab, "lab_set", {"device": "fume_hood", "channel": "power",
                                 "value": "on"})["ok"]
    rid = call(lab, "lab_set", {"device": "pump", "channel": "flow",
                                "value": 1, "unit": "mL/min"})[
        "approval_request"]
    env = dict(os.environ, DELFIN_AGENT_SESSION="4242")
    res = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "lab",
                          "approve", rid], cwd=str(REPO), env=env,
                         capture_output=True, text=True, timeout=120)
    assert res.returncode != 0 and "person" in res.stderr
    res = subprocess.run([sys.executable, "-m", "delfin.agent.cli", "lab",
                          "approve", rid], cwd=str(REPO),
                         env={k: v for k, v in os.environ.items()
                              if k != "DELFIN_AGENT_SESSION"},
                         capture_output=True, text=True, timeout=120)
    assert res.returncode == 0, res.stderr               # control: a person


def test_the_lab_directory_is_not_writable_through_write_tools(lab):
    from delfin.lab import store
    target = store.lab_dir() / "manifest.pin.json"
    out = call(lab, "write_file", {"path": str(target), "content": "{}"})
    assert "error" in out, out
    assert '"checksum"' in target.read_text()


def test_smoke_prints_ok(tmp_path):
    env = {k: v for k, v in os.environ.items()
           if k not in ("DELFIN_AGENT_SESSION", "DELFIN_LAB_DIR",
                        "DELFIN_LAB_MANIFEST")}
    env["HOME"] = str(tmp_path)
    res = subprocess.run([sys.executable, "-m", "delfin.lab.smoke"],
                         cwd=str(REPO), env=env, capture_output=True,
                         text=True, timeout=120)
    assert res.stdout.strip() == "smoke lab ok", res.stdout + res.stderr


def test_the_pi_units_are_hardened_and_installed_by_the_deploy():
    units = REPO / "scripts" / "pi" / "lab"
    for name in ("delfin-lab-estop.service", "delfin-lab-watchdog.service"):
        text = (units / name).read_text()
        for line in ("NoNewPrivileges=yes", "ProtectSystem=strict",
                     "RestrictAddressFamilies=AF_UNIX", "PrivateTmp=yes",
                     "Restart=always", "ReadWritePaths="):
            assert line in text, (name, line)
        assert "lab daemon --role" in text
    install = (REPO / "scripts" / "pi" / "pi_install.sh").read_text()
    assert "install_lab_units" in install
    assert "delfin.lab.smoke" in install
