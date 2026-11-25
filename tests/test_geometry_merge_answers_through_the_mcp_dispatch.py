"""``geometry_merge`` answers, instead of raising on its own success.

Input: two xyz fragments, through ``execute_mcp_call`` -- the one MCP dispatch
the agent loop and Code Mode both take (``api_client.execute_mcp_call``), not
the handler directly.
Output: the concatenated structure, with the atoms of the earlier fragment
first.

What was wrong: ``delfin/tools/geometry_edit.py:712`` put a ``range`` object
into the ``change`` dict (``"added": range(len(a), len(out))``). The
concatenation itself succeeded; the answer then died in ``json.dumps`` inside
``delfin/ops_server/server.py`` (``_dumps``), so the tool answered
``{"error": "TypeError: Object of type range is not JSON serializable"}`` for
EVERY call -- two fragments, ten, a single atom. The one test that covered
``merge_fragments`` (``tests/test_geometry_edit.py:284``) calls the helper and
reads ``change`` nowhere, so a change value that cannot be serialised stayed
invisible: it is in the return type, not in the structure.

The other ``range`` in that file (``geometry_edit.py:655``, the mirror op's
``atoms_moved``) is already wrapped in ``list(...)``; this one was not, and the
same op's contract -- ``{op, atoms_moved, ...}`` as JSON -- says the value is a
list of atom indices.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

_WATER = "3\nwater\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n"
_METHANE = ("5\nmethane\nC 0 2 0\nH 0.6 2 0\nH 0 2.6 0\nH 0 2 0.6\n"
            "H -0.6 2 0\n")


def _call(tmp_path, tool: str, args: dict) -> dict:
    """One MCP tool call through the product's own dispatch."""
    from delfin.agent import api_client as A
    from delfin.agent import mcp_client as M

    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    os.environ["HOME"] = str(tmp_path / "home")
    (tmp_path / "home").mkdir(exist_ok=True)
    os.environ.pop("DELFIN_AGENT_SESSION", None)
    M.reset_registry(ws)
    M.get_registry(ws)

    perms = A.KitToolPermissions(workspace=ws)
    perms.mode = "bypassPermissions"
    perms.confirm_callback = lambda *a, **k: True
    perms.task_session_id = "s1"
    executor = A._DocToolExecutor()
    raw = A.execute_mcp_call(f"mcp__delfin-ops__{tool}", args, perms, executor)

    from delfin.agent import untrusted
    out = json.loads(untrusted.unwrap(raw))
    assert isinstance(out, dict), out
    return out


def test_merge_answers_with_the_concatenated_structure(tmp_path):
    out = _call(tmp_path, "geometry_merge", {"fragments": [_WATER, _METHANE]})
    assert "error" not in out, out
    assert out["xyz"].split("\n")[0].strip() == "8", out["xyz"][:80]
    # atoms of the earlier fragment first, unmoved
    body = out["xyz"].splitlines()[2:5]
    assert body[0].startswith("O"), body
    assert body[1].startswith("H"), body


def test_the_change_names_the_added_atoms_as_a_list(tmp_path):
    """The value ``json`` can carry: the atom indices the merge appended."""
    out = _call(tmp_path, "geometry_merge", {"fragments": [_WATER, _METHANE]})
    added = out["change"]["added"]
    assert isinstance(added, list), out["change"]
    assert added == [3, 4, 5, 6, 7], out["change"]


def test_control_an_empty_fragment_list_is_still_its_own_refusal(tmp_path):
    """Control: the fix is not a blanket. No fragments was never the bug."""
    out = _call(tmp_path, "geometry_merge", {"fragments": []})
    assert "error" in out, out
