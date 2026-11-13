"""bypassPermissions runs a command only where a layer confines it.

Bypass asks nothing per command, so nothing but the cage stands between a
command and the account's files. Pinned here:

* where the native shell command is built (``_bash_isolation_argv``): a
  bypass session with no layer, or with isolation switched off, gets a
  refusal instead of a plain shell (the start check lowers such a start to
  default; this holds for any other route into bypass, and for a layer
  that stops working mid-run);
* the CLI backend: its shell runs outside DELFIN's process, and in bypass
  it was started with every check off. Its gate hook now gives a Bash call
  the whole native gate and then DELFIN's cage, as the command text the
  CLI runs instead (``updatedInput``); the cage's network goes through the
  egress proxy of the DELFIN process, whose socket the state carries.

Measured with the real CLI against a scripted endpoint (2026-11-13): on
origin/main a command that builds an outside path at run time wrote there
and a force push reached git; on this branch the cage refused the write
(read-only file system) and the deny list the push.
"""

from __future__ import annotations

import io
import json
import shlex

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli_gate, egress_proxy, kit_settings


def _host(monkeypatch, *, bwrap=False, landlock=False, seatbelt=False):
    monkeypatch.setattr(A.shutil, "which",
                        lambda name, *a, **k: "/usr/bin/bwrap" if (
                            name == "bwrap" and bwrap) else None)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: bwrap)
    monkeypatch.setattr(A, "_landlock_functional", lambda: landlock)
    monkeypatch.setattr(A, "_seatbelt_functional", lambda: seatbelt)
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)


def _refused(argv) -> str:
    return argv[4] if argv[:2] == ["/bin/sh", "-c"] and argv[3:4] == ["refused"] else ""


@pytest.fixture(autouse=True)
def _fresh():
    A.set_bash_isolation_override("")
    yield
    A.set_bash_isolation_override("")


# -- native: where the command is built --------------------------------------

def test_bypass_without_a_layer_is_refused_where_the_command_is_built(
        monkeypatch, tmp_path):
    _host(monkeypatch)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    why = _refused(A._bash_isolation_argv("touch x", tmp_path, perms, mode="auto"))
    assert "bypassPermissions" in why and "confines it" in why


def test_bypass_with_isolation_off_is_refused(monkeypatch, tmp_path):
    _host(monkeypatch, bwrap=True)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    assert "bypassPermissions" in _refused(
        A._bash_isolation_argv("touch x", tmp_path, perms, mode="off"))


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
@pytest.mark.parametrize("iso", ["auto", "off"])
def test_the_modes_that_ask_keep_their_shell_without_a_layer(
        monkeypatch, tmp_path, mode, iso):
    """Unchanged: there a person reads the command before it runs."""
    _host(monkeypatch)
    perms = A.KitToolPermissions(workspace=tmp_path, mode=mode)
    assert not _refused(A._bash_isolation_argv("touch x", tmp_path, perms, mode=iso))


# -- the CLI backend -----------------------------------------------------------

@pytest.fixture
def scene(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        home / ".delfin" / "settings.json")
    monkeypatch.setattr(egress_proxy, "network_mode", lambda: "none")
    ws = tmp_path / "ws"
    ws.mkdir()

    def build(mode="bypassPermissions"):
        token = cli_gate.new_token()
        cli_gate.write_state(token, workspace=str(ws), mode=mode, can_ask=False)
        return token

    return build, ws


def _ev(cmd, cwd):
    return {"tool_name": "Bash", "tool_input": {"command": cmd, "description": "d"},
            "cwd": str(cwd)}


def test_a_bypass_command_runs_in_the_cage(scene, monkeypatch):
    build, ws = scene
    built = []

    def cage(cmd, run_cwd, perms, mode=None, **kw):
        built.append((cmd, str(run_cwd), perms.mode, mode))
        return ["bwrap", "--ro-bind", "/", "/", "--", "/bin/bash", "-c", cmd]

    monkeypatch.setattr(A, "_bash_isolation_argv", cage)
    refusal, caged = cli_gate.decide_call(build(), _ev("touch made.txt", ws))
    assert refusal is None
    assert built == [("touch made.txt", str(ws), "bypassPermissions", "bwrap")]
    assert shlex.split(caged)[-1] == "touch made.txt"
    assert shlex.split(caged)[0] == "bwrap"


def test_without_a_layer_the_bypass_command_is_refused(scene, monkeypatch):
    build, ws = scene
    _host(monkeypatch)
    refusal, caged = cli_gate.decide_call(build(), _ev("touch made.txt", ws))
    assert caged is None
    assert refusal and "refused" in refusal


def test_the_native_gate_still_refuses_in_bypass(scene, monkeypatch):
    build, ws = scene
    monkeypatch.setattr(A, "_bash_isolation_argv",
                        lambda cmd, *a, **k: ["bwrap", "--", "/bin/bash", "-c", cmd])
    refusal, caged = cli_gate.decide_call(build(), _ev("git push --force origin HEAD", ws))
    assert caged is None and "deny-pattern" in refusal


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "plan"])
def test_other_modes_are_not_rewritten(scene, monkeypatch, mode):
    build, ws = scene
    monkeypatch.setattr(A, "_bash_isolation_argv",
                        lambda *a, **k: pytest.fail("caged outside bypass"))
    _refusal, caged = cli_gate.decide_call(build(mode), _ev("ls", ws))
    assert caged is None


def test_the_hook_hands_the_cli_the_caged_command(scene, monkeypatch, capsys):
    build, ws = scene
    monkeypatch.setattr(A, "_bash_isolation_argv",
                        lambda cmd, *a, **k: ["bwrap", "--", "/bin/bash", "-c", cmd])
    token = build()
    monkeypatch.setattr("sys.stdin", io.StringIO(json.dumps(_ev("touch made.txt", ws))))
    assert cli_gate.main(["cli_gate", token]) == 0
    out = json.loads(capsys.readouterr().out)["hookSpecificOutput"]
    assert out["permissionDecision"] == "allow"
    assert out["updatedInput"]["command"] == "bwrap -- /bin/bash -c 'touch made.txt'"
    assert out["updatedInput"]["description"] == "d"


def test_bypass_state_carries_this_process_s_proxy(scene, monkeypatch):
    build, ws = scene
    monkeypatch.setattr(egress_proxy, "network_mode", lambda: "proxy")

    class _Proxy:
        socket_path = "/run/proxy/sock"

        def token_file_for(self, perms):
            return "/run/proxy/token-1"

    monkeypatch.setattr(egress_proxy, "get_proxy", lambda: _Proxy())
    state = json.loads(cli_gate._state_path(build()).read_text())
    assert state["proxy_socket"] == "/run/proxy/sock"
    assert state["proxy_token_file"] == "/run/proxy/token-1"
    state = json.loads(cli_gate._state_path(build("default")).read_text())
    assert "proxy_socket" not in state


def test_an_external_proxy_is_what_the_cage_binds(monkeypatch):
    monkeypatch.setattr(egress_proxy, "_PROXY", None)
    egress_proxy.use_external("/run/p/sock", "/run/p/token")
    proxy = egress_proxy.get_proxy()
    assert proxy.socket_path == "/run/p/sock"
    assert proxy.token_file_for(object()) == "/run/p/token"
