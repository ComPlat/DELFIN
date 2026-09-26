"""What a continued session forgets: the refusals.

A session can be continued with ``delfin-agent chat -r <session_id>``. The
engine carries a declared list of fields across that boundary
(``AgentEngine._SESSION_FIELDS``). The permission object
(``KitToolPermissions`` on ``client._permissions``) is NOT among them, so a
refusal -- the human's "no" for THIS session -- dies with the process:

observed on the night of 26 Sep 2026: a session was continued with -r after
a refusal and asked for the SAME denied file again.

This file holds the red control: build an engine whose permissions carry a
denial, export it, restore into a fresh engine, and show the denial is gone.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.engine import AgentEngine


def _make_engine_with_denial(tmp_path: Path) -> AgentEngine:
    """An engine whose client carries one denied path, one denied action
    and a refusal with a reason -- the state after the user said no."""
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.refusal_memory import Refusal, RefusalMemory

    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())

    class _Client:
        pass

    perms = KitToolPermissions(workspace=tmp_path)
    perms.denied_paths = {str(tmp_path / "secret.env")}
    perms.denied_actions = {"bash:rm -rf /": "rm -rf /"}
    perms.refusal_memory = RefusalMemory()
    perms.refusal_memory.record(Refusal(
        tool="read_file", target=str(tmp_path / "secret.env"),
        reason="operator said no", time="2026-09-26T00:00:00Z"))
    eng.client = _Client()
    eng.client._permissions = perms
    return eng


def test_a_denied_path_survives_a_resume(tmp_path):
    """The denial is the human's statement about this session; it must
    survive ``chat -r``. Red on the previous commit: restore_state never
    touched the permissions object, so the resumed session asked again."""
    eng = _make_engine_with_denial(tmp_path)
    data = eng.export_state()

    # The resumed engine is the same machine, same workspace, but a FRESH
    # process: its client rebuilt its permissions from the launch arguments.
    resumed = _make_engine_with_denial(tmp_path)
    resumed.client._permissions.denied_paths = set()  # fresh process
    resumed.client._permissions.denied_actions = {}
    resumed.client._permissions.refusal_memory = None
    resumed.restore_state(data)

    assert str(tmp_path / "secret.env") in \
        resumed.kit_permissions.denied_paths


def test_a_denied_action_survives_a_resume(tmp_path):
    eng = _make_engine_with_denial(tmp_path)
    data = eng.export_state()

    resumed = _make_engine_with_denial(tmp_path)
    resumed.client._permissions.denied_paths = set()
    resumed.client._permissions.denied_actions = {}
    resumed.client._permissions.refusal_memory = None
    resumed.restore_state(data)

    assert resumed.kit_permissions.denied_actions == {"bash:rm -rf /":
                                                      "rm -rf /"}


def test_the_refusal_memory_survives_a_resume(tmp_path):
    eng = _make_engine_with_denial(tmp_path)
    data = eng.export_state()

    resumed = _make_engine_with_denial(tmp_path)
    resumed.client._permissions.denied_paths = set()
    resumed.client._permissions.denied_actions = {}
    resumed.client._permissions.refusal_memory = None
    resumed.restore_state(data)

    mem = resumed.kit_permissions.refusal_memory
    assert mem is not None
    entries = mem.to_dict().get("entries", [])
    assert any("operator said no" in str(e) for e in entries)

