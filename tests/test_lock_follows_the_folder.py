"""The lock is carried by the folder as well as by the role.

The role alone was wrong in one direction. Switching mode mid-session
restamped the agent role on the next turn while the engine — and with it
the workspace — was kept: the session became scope-locked while still
pointing at the previous folder. Containment aimed at a place nobody
meant, and an office session wrote documents into a code checkout.

Three carriers now, ORed, and removing any one cannot unlock a session.
That redundancy is deliberate for a containment boundary and is the
intended end state, not a migration step.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A
from delfin.dashboard.tab_agent import _mode_workspace_differs


@pytest.fixture
def office_root(monkeypatch):
    root = pathlib.Path(tempfile.mkdtemp(prefix="OFFICE_")).resolve()
    monkeypatch.setattr(A, "_LOCKED_ROOTS_CACHE", (root,))
    return root


def _locked(ws, role=""):
    return A.KitToolPermissions(workspace=ws, agent_role=role).scope_locked


# ---------------------------------------------------------------------------
# The workspace carrier
# ---------------------------------------------------------------------------

def test_the_office_folder_is_locked_without_any_role(office_root):
    assert _locked(office_root) is True


def test_a_subfolder_of_it_is_locked_too(office_root):
    sub = office_root / "Reisekosten"
    sub.mkdir()
    assert _locked(sub) is True


def test_an_unrelated_folder_is_not_locked(office_root):
    other = pathlib.Path(tempfile.mkdtemp(prefix="code_"))
    assert _locked(other) is False


def test_the_role_still_locks_on_its_own(office_root):
    """Carrier two: removing the workspace term must not unlock anything."""
    other = pathlib.Path(tempfile.mkdtemp(prefix="code_"))
    assert _locked(other, "office_agent") is True


def test_the_explicit_flag_still_locks(office_root):
    """Carrier three."""
    other = pathlib.Path(tempfile.mkdtemp(prefix="code_"))
    perms = A.KitToolPermissions(workspace=other)
    perms.lock_workspace = True
    assert perms.scope_locked is True


def test_a_sibling_with_the_same_prefix_is_not_locked(office_root):
    """`/x/office_old` must not match because it starts with `/x/office`."""
    sibling = office_root.parent / (office_root.name + "_old")
    sibling.mkdir(exist_ok=True)
    assert _locked(sibling) is False


# ---------------------------------------------------------------------------
# Failure must never unlock
# ---------------------------------------------------------------------------

def test_unreadable_settings_still_lock_the_default_folder(monkeypatch):
    """"No office dir configured" must not read as "nothing is locked"."""
    monkeypatch.setattr(A, "_LOCKED_ROOTS_CACHE", None)
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda *a, **kw: (_ for _ in ()).throw(OSError("unreadable")))
    roots = A._locked_workspace_roots()
    assert (pathlib.Path.home() / "office").resolve() in roots


def test_a_missing_workspace_is_not_locked_and_does_not_raise(office_root):
    assert A._workspace_is_locked(None) is False


# ---------------------------------------------------------------------------
# A subagent inherits by recomputation, not by copy
# ---------------------------------------------------------------------------

def test_a_child_relocated_into_the_folder_is_locked(office_root):
    """dataclasses.replace re-runs __post_init__, so the child is judged on
    where it actually is — which a copied flag would not do."""
    import dataclasses

    parent = A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="code_")))
    assert parent.scope_locked is False
    child = dataclasses.replace(parent, workspace=office_root)
    assert child.scope_locked is True


# ---------------------------------------------------------------------------
# The other half: the engine must follow the folder
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("old,new,expected", [
    ("solo", "office", True),
    ("office", "solo", True),
    ("dashboard", "office", True),
    ("solo", "dashboard", True),
    ("solo", "pipeline", False),
    ("solo", "solo", False),
])
def test_only_a_folder_changing_switch_rebuilds(old, new, expected):
    """Rebuilding on every switch would discard a live session for nothing."""
    assert _mode_workspace_differs(old, new) is expected


def test_the_handler_drops_the_engine(monkeypatch):
    source = pathlib.Path(
        __import__("delfin.dashboard.tab_agent", fromlist=["x"]).__file__
    ).read_text(encoding="utf-8")
    head = source[source.index("def _on_mode_change"):]
    head = head[:head.index("# Update mode description label")]
    assert "_mode_workspace_differs(old_mode, new_mode)" in head
    assert 'state["engine"] = None' in head
    assert 'not state["streaming"]' in head, (
        "dropping the engine mid-turn would orphan the running turn")
