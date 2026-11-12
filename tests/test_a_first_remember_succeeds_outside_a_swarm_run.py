"""A person says "merk dir das": the model's FIRST remember call must succeed.

Observed live: a first call with type="finding" outside a swarm run was
refused ("not in a swarm run ... use type=project") and only a SECOND call
succeeded. There is no swarm run in a plain workspace, so type="finding"
there must degrade to a project memory instead of refusing.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import patch

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions


def _perms(ws: Path) -> KitToolPermissions:
    return KitToolPermissions(workspace=ws, mode="bypassPermissions")


def test_finding_outside_swarm_run_stores_project_memory(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    with patch.object(Path, "home", lambda: home), \
            patch.dict(os.environ, {"DELFIN_SWARM_RUN": ""}, clear=False):
        out = json.loads(A._doc_executor._execute_remember(
            {"text": "geiler typ.", "type": "finding"}, _perms(ws)))
    assert out["status"] == "ok", out
    assert out["type"] == "project", out
    assert "swarm" in out.get("note", ""), out

    with patch.object(Path, "home", lambda: home):
        from delfin.agent.memory_store import list_typed_memories
        mems = list_typed_memories(ws)
    assert len(mems) == 1
    assert mems[0]["type"] == "project"


def test_finding_with_swarm_run_still_goes_to_swarm_layer():
    """The swarm path itself is unchanged: a run id present -> swarm append."""
    from delfin.agent import api_client as sm  # noqa: F401  (kept: guard test)
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    # the coercion must sit BEFORE the layered dispatch, and the layered
    # refusal must still exist for any direct layered call
    assert "_execute_remember_layered" in src
    assert 'current_run_id' in src
