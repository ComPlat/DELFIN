"""An approved apply_patch to a protected file reaches THAT file -- and never
silently writes to another path via an absolute path in a diff header.

Regression for the 2026-10-29 incident: an agent's apply_patch adding a class
to delfin/agent/api_client.py was approved, yet the class was absent from that
file afterwards and the agent re-sent the same change as edit_file (which
landed). Cause reproduced here: _apply_python joins the diff's target path
onto the workspace with ``Path(workspace) / rel``; pathlib discards the
workspace prefix when ``rel`` is ABSOLUTE, so an absolute path in the ``+++``
header makes the applier write to whatever absolute target is named -- outside
the workspace sandbox and outside the file the gate asked the user about --
while reporting status ``ok`` and journaling it under the RELATIVE protected
path.

Contract pinned here:
  * with an approving confirm callback, a patch whose ``+++`` header names the
    protected file RELATIVELY is applied to the workspace file;
  * a patch whose ``+++`` header names an ABSOLUTE path outside the workspace
    must NOT write there: it is refused with an error that says why, or is
    normalized and applied to the relative path -- never silently writing
    outside and reporting ok.

Gate path under test: _execute_apply_patch (api_client) ->
_run_permission_gate("apply_patch", ...) -> _gate_write_path per diff file ->
matches_path_protected -> _ask_confirm -> _pa.apply_patch -> _apply_python.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor

ORIG = (
    "class First(object):\n"
    "    pass\n"
)


def _perms(root, sid="s-protected"):
    p = KitToolPermissions(workspace=str(root), mode="acceptEdits")
    p.task_session_id = sid
    return p


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "delfin").mkdir()
    (tmp_path / "delfin" / "agent").mkdir()
    tgt = tmp_path / "delfin" / "agent" / "api_client.py"
    tgt.write_text(ORIG)
    return tmp_path, tgt


def _approving_callback(name, args, preview):
    """A fake orchestrator that approves every self-modification request."""
    assert name == "apply_patch", name
    return True


def _apply(perms, diff):
    out = json.loads(_doc_executor.execute("apply_patch", {"diff": diff}, perms))
    return out


def _patch(status, header):
    return (
        f"--- a/delfin/agent/api_client.py\n"
        f"{header}\n"
        "@@ -1 +1,2 @@\n"
        " class First(object):\n"
        "+class Second(object):\n"
        "     pass\n"
    )


def test_approved_patch_to_protected_file_lands(ws):
    root, tgt = ws
    perms = _perms(root)
    perms.confirm_callback = _approving_callback
    out = _apply(perms, _patch("--- a/delfin/agent/api_client.py",
                               "+++ b/delfin/agent/api_client.py"))
    if out.get("status") == "ok":
        assert "class Second" in tgt.read_text(), (
            f"approved apply_patch reported ok but protected file unchanged; "
            f"full result: {out}")
    else:
        assert "error" in out, out


def test_approved_patch_with_a_prefix_on_plus_side_lands(ws):
    """The incident form: a header written ``+++ a/<rel>`` (a/ on the plus
    side) previously made the applier target a/<rel> and the protected file
    never changed; it must now land on the protected relative path."""
    root, tgt = ws
    perms = _perms(root)
    perms.confirm_callback = _approving_callback
    out = _apply(perms, _patch("--- a/delfin/agent/api_client.py",
                               "+++ a/delfin/agent/api_client.py"))
    if out.get("status") == "ok":
        assert "class Second" in tgt.read_text(), (
            f"approved apply_patch (+++ a/<rel>) reported ok but protected "
            f"file unchanged; full result: {out}")
    else:
        assert "error" in out, out


def test_absolute_path_in_header_does_not_write_outside(ws, tmp_path):
    root, tgt = ws
    perms = _perms(root)
    perms.confirm_callback = _approving_callback
    # A DIFFERENT location than the workspace protected file: the applier must
    # not silently write here just because the header names an absolute path.
    outside = tmp_path / "outside" / "elsewhere.py"
    outside.parent.mkdir()
    outside.write_text(ORIG)
    out = _apply(
        perms,
        _patch("--- a/delfin/agent/api_client.py",
               f"+++ {outside}"))
    # The outside file the header named must be untouched.
    assert "class Second" not in outside.read_text(), (
        f"apply_patch wrote OUTSIDE the workspace to the absolute path named "
        f"in the header; full result: {out}")
