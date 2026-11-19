"""`!cmd` works on the CLI backend, under that backend's own boundary.

``engine.run_gated_bash`` -- what every surface's `!` calls -- read the
policy from ``client._permissions``. Only the OpenAI-compatible client
carries one, so on the CLI backend, the ANTHROPIC API backend and the
Codex backend the shell escape answered

    this backend has no permission gate, so it has no shell

and did so on the FACTORY DEFAULT, which is the CLI backend. Measured
2026-11-19: `delfin agent chat` with no --provider reports "provider
claude", and `!` is dead there on the terminal and on the agent tab
alike.

The sentence was also wrong about the CLI backend. That backend does
gate the model's own tool calls -- through ``cli_gate``, a PreToolUse
hook whose policy lives in a state file keyed by the client's token. The
policy was there; nothing read it.

It is fetched now rather than rebuilt, so the shell escape is bound by
exactly what the model's tool calls are bound by. It cannot be looser:
the command additionally goes through DELFIN's own bash executor, whose
deny-list and secret scan the hook never applies to the external
process.

The API and Codex backends have no DELFIN gate at all and still refuse
-- there the sentence is true, and it now names the backend.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import cli_gate
from delfin.agent.engine import AgentEngine


@pytest.fixture
def own_home(monkeypatch, tmp_path):
    """cli_gate's state lives under the home directory; never the real one."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    return home


class _CLIClient:
    """The shape that matters: a gate token, no ``_permissions``."""

    def __init__(self, token, workspace):
        self._gate_token = token
        self._workspace = workspace
        self.writes = 0

    def _write_gate_state(self):
        self.writes += 1
        cli_gate.write_state(self._gate_token, workspace=str(self._workspace),
                             mode="default", can_ask=False)


def _engine(client):
    eng = AgentEngine.__new__(AgentEngine)
    eng.client = client
    return eng


def test_the_command_actually_runs(own_home, tmp_path):
    workspace = tmp_path / "ws"
    workspace.mkdir()
    eng = _engine(_CLIClient(cli_gate.new_token(), workspace))

    out = json.loads(eng.run_gated_bash("echo hi"))

    assert out.get("exit_code") == 0, out
    assert "hi" in str(out.get("stdout", "")), out


def test_the_state_is_written_when_no_callback_ever_set_one(own_home,
                                                            tmp_path):
    """``_write_gate_state`` runs when the confirm callback is set. A
    session that never set one had no state file, and the escape would
    have found no workspace and refused."""
    workspace = tmp_path / "ws"
    workspace.mkdir()
    client = _CLIClient(cli_gate.new_token(), workspace)

    _engine(client).run_gated_bash("echo hi")

    assert client.writes == 1


def test_the_boundary_is_the_one_the_hook_enforces(own_home, tmp_path):
    """Not a fresh policy built here: the same object the hook builds for
    the model's own tool calls, so the two cannot drift apart."""
    workspace = tmp_path / "ws"
    workspace.mkdir()
    token = cli_gate.new_token()
    cli_gate.write_state(token, workspace=str(workspace), mode="default",
                         can_ask=False)

    perms = cli_gate.permissions_for(token)

    assert perms is not None
    assert Path(perms.workspace) == workspace.resolve() or (
        Path(perms.workspace) == workspace)


def test_a_token_with_no_state_grants_nothing(own_home):
    """Fail closed. An absent or unreadable state file records no
    workspace, so there is no boundary -- which is a refusal, never a
    command run without one."""
    assert cli_gate.permissions_for(cli_gate.new_token()) is None
    assert cli_gate.permissions_for("") is None


def test_an_unreadable_state_file_is_not_read_as_no_restriction(own_home):
    token = cli_gate.new_token()
    path = cli_gate._state_path(token)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{ this is not json", encoding="utf-8")

    assert cli_gate.permissions_for(token) is None


def test_a_backend_with_no_gate_at_all_still_refuses_and_names_itself(
        own_home):
    """The API and Codex backends gate nothing of DELFIN's. The refusal
    stands; it just stops blaming an unnamed 'this backend'."""
    class _NoGate:
        pass

    eng = _engine(_NoGate())

    with pytest.raises(RuntimeError) as caught:
        eng.run_gated_bash("echo hi")

    assert "_NoGate" in str(caught.value), str(caught.value)
    assert "no permission gate" in str(caught.value)


def test_the_deny_list_still_applies_on_this_route(own_home, tmp_path):
    """The route must not be a way around what the agent may not do: the
    command goes through DELFIN's bash executor, deny-list included."""
    workspace = tmp_path / "ws"
    workspace.mkdir()
    eng = _engine(_CLIClient(cli_gate.new_token(), workspace))

    out = json.loads(eng.run_gated_bash("rm -rf /"))

    assert out.get("error"), out
    assert out.get("exit_code") != 0
