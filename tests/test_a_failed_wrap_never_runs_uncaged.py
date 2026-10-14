"""A bubblewrap that fails when the command is started never leads to the
command running outside the cage, and it withdraws the caged admission.

``_bwrap_functional`` probes the wrap shape once per process and caches
the answer. The real wrap can still fail later (a nested sandbox answers
"Creating new namespace failed: Device or resource busy"; a bind source
that vanished). Two properties are pinned here:

(a) the executor starts the bwrap argv and nothing else: no fallback to a
    plain shell, so the marker file the command would write is absent,
    and the result says the sandbox could not be applied (not a bare
    exit code with bwrap's stderr);
(b) after one such failure in the process, ``_caged_command_may_run_unasked``
    answers False, so the next off-list command is asked about instead of
    being admitted unasked and failing the same way unseen.

The failing bwrap is a stand-in placed first on PATH; the isolation
decision is pinned to "bwrap", so the real argv builder and the real
executor run on every host.

The case "Landlock reported unavailable: an off-list harmless command in
default mode is asked" is pinned in
tests/test_a_caged_command_runs_without_asking.py
(test_control_without_landlock_read_support_asks at the gate,
test_without_landlock_the_command_is_asked_about through the executor).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


_BUSY = "bwrap: Creating new namespace failed: Device or resource busy"


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    d.mkdir()
    (d / "Makefile").write_text("all:\n\techo made\n")
    # A bwrap that fails before running anything, as a nested sandbox does.
    fake = tmp_path / "fakebin"
    fake.mkdir()
    calls = tmp_path / "bwrap-calls"
    (fake / "bwrap").write_text(
        "#!/bin/sh\n"
        f"echo called >> '{calls}'\n"
        f"echo '{_BUSY}' >&2\n"
        "exit 126\n")
    os.chmod(fake / "bwrap", 0o755)
    monkeypatch.setenv("PATH", f"{fake}{os.pathsep}{os.environ.get('PATH', '')}")
    # The isolation decision says bubblewrap holds commands here; the
    # process cage and Landlock read confinement are reported working, so
    # the caged admission is open before the failure.
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(A, "_BWRAP_RUNTIME_FAILURE", None, raising=False)
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "proxy")
    return d


_MARKER_CMD = 'python3 -c "open(\'marker.txt\',\'w\').write(\'ran\')"'


def _run(ws: Path, cmd: str, mode: str, attended: bool = False) -> dict:
    asked: list[str] = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=_cb if attended else None)
    out = json.loads(_doc_executor.execute("bash", {"command": cmd}, perms))
    out["_asked"] = asked
    return out


@pytest.mark.parametrize("mode", ["default", "bypassPermissions"])
def test_a_failed_wrap_does_not_run_the_command(ws, mode):
    """default: admitted by the caged rule; bypass: wrapped as usual.
    Either way the stand-in bwrap was started, failed, and the command
    did not run anywhere else."""
    out = _run(ws, _MARKER_CMD, mode)
    assert (ws.parent / "bwrap-calls").is_file(), out  # the wrap was attempted
    assert not (ws / "marker.txt").exists(), out
    assert "sandbox could not be applied" in out.get("error", ""), out
    assert "Device or resource busy" in json.dumps(out)


def test_after_a_failed_wrap_the_caged_admission_is_withdrawn(ws):
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    assert A._caged_command_may_run_unasked("make", perms, "default")
    _run(ws, _MARKER_CMD, "default")
    assert not A._caged_command_may_run_unasked("make", perms, "default")
    # The next off-list command is asked about, not admitted unseen.
    out = _run(ws, "make", "default", attended=True)
    assert out["_asked"] and out["_asked"][0].startswith("$ make"), out
    # Unattended, the old refusal names the auto-allow list.
    out = _run(ws, "make", "default")
    assert "not on the auto-allow list" in out.get("error", ""), out


def test_a_command_that_merely_fails_withdraws_nothing(ws, monkeypatch):
    """Control: a failing command under a working wrap is the command's
    own failure and leaves the admission as it was."""
    ok = ws.parent / "okbin"
    ok.mkdir()
    (ok / "bwrap").write_text(
        "#!/bin/sh\necho 'make: *** no rule' >&2\nexit 2\n")
    os.chmod(ok / "bwrap", 0o755)
    monkeypatch.setenv("PATH", f"{ok}{os.pathsep}{os.environ['PATH']}")
    out = _run(ws, "make nothing", "default")
    assert out.get("exit_code") == 2 and "error" not in out, out
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    assert A._caged_command_may_run_unasked("make", perms, "default")
