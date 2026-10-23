"""Grants do NOT survive a resume -- and that is a promise.

The confirmation dialog grants a read or a push "for the rest of this
session. Nothing is saved -- a later session asks again". Refusals are
carried across `chat -r` (they are the human's "no", and asking the
same question again after a "no" is a retry loop). Grants are the
opposite direction: silently inheriting a permission the user gave a
previous process would widen what this one may do, without a dialog.

Whether a resumed session counts as "the same session" for grants is
not ours to decide: in doubt DELFIN asks once too often, which is the
safe side. These tests pin that boundary so a future change must argue
with the dialog text, not with inertia.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.engine import AgentEngine


def _engine_with_grants(tmp_path: Path) -> AgentEngine:
    """An engine whose permissions hold a read grant and a push grant."""
    from delfin.agent.api_client import KitToolPermissions

    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())

    class _Client:
        pass

    perms = KitToolPermissions(workspace=tmp_path)
    perms.session_read_dirs = (Path("/etc"),)
    perms.push_grants = {"origin/main": 1}
    eng.client = _Client()
    eng.client._permissions = perms
    return eng


def _fresh_engine(tmp_path: Path) -> AgentEngine:
    """The resumed process: same machine, client rebuilt from launch
    arguments -- no grants."""
    eng = _engine_with_grants(tmp_path)
    perms = eng.kit_permissions
    perms.session_read_dirs = ()
    perms.push_grants = {}
    return eng


def test_a_read_grant_does_not_survive_a_resume(tmp_path):
    """The dialog says "for the rest of this session. Nothing is saved".
    Carrying a read grant into a resumed process would let this one read
    a directory another process was granted -- no dialog, no decision by
    the user of THIS conversation."""
    data = _engine_with_grants(tmp_path).export_state()

    resumed = _fresh_engine(tmp_path)
    resumed.restore_state(data)

    assert resumed.kit_permissions.session_read_dirs == ()


def test_a_push_grant_does_not_survive_a_resume(tmp_path):
    """One request grants one push. A carried grant would give the
    resumed session a push nobody in this conversation approved."""
    data = _engine_with_grants(tmp_path).export_state()

    resumed = _fresh_engine(tmp_path)
    resumed.restore_state(data)

    assert resumed.kit_permissions.push_grants == {}


def test_the_export_carries_no_grant_data(tmp_path):
    """The "refusals" blob must contain denials only. A grant that leaks
    into the export is a refusal away from being restored -- and this
    test is the cheap place to notice."""
    data = _engine_with_grants(tmp_path).export_state()
    blob = data.get("refusals", {})

    assert "session_read_dirs" not in blob
    assert "push_grants" not in blob
    assert set(blob.get("denied_paths", [])) == set()
