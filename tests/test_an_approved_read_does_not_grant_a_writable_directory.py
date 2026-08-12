"""The dialog said "this single read". The click granted a writable root.

Approving an outside-workspace read ran ``add_extra_dir`` on the file's
PARENT directory, and ``all_workspace_roots`` / ``find_root_for`` treat an
extra dir as writable. So one approved read turned a directory into
somewhere the agent may create, overwrite and delete files — while the
preview said the opposite:

    "Approving this read does NOT add the directory to the permanent
     workspace list — it only allows this single read."

``add_extra_dir`` also had no forbidden-root check, though the constructor
does and refuses $HOME by name ("it is your home directory ... which would
let the agent write anywhere"). Approving a read of ``~/notes.txt`` handed
the session exactly that root at runtime.

A read approval is now a read grant: reachable for reads, refused for
writes, never persisted — and the dialog says so.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def home(tmp_path_factory, monkeypatch):
    """A fake $HOME, sibling to the workspaces so the floor stays meaningful."""
    h = tmp_path_factory.mktemp("fake_home")
    monkeypatch.setattr(Path, "home", lambda: h)
    return h


class _Approve:
    def __init__(self, answer=True):
        self.answer = answer
        self.last_timed_out = False
        self.previews: list[str] = []

    def callback(self, name, args, preview):
        self.previews.append(preview)
        return self.answer


@pytest.fixture
def ex():
    return A._DocToolExecutor.__new__(A._DocToolExecutor)


def _perms(ws, broker):
    return A.KitToolPermissions(
        workspace=ws, mode="default", confirm_callback=broker.callback)


# ---------------------------------------------------------------------------
# What the approval grants
# ---------------------------------------------------------------------------

def test_an_approved_read_opens_the_directory_for_reading_only(
        ex, tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "notes"
    outside.mkdir()
    (outside / "a.txt").write_text("data", encoding="utf-8")
    (outside / "b.txt").write_text("more", encoding="utf-8")

    broker = _Approve(answer=True)
    perms = _perms(ws, broker)

    assert ex._check_read_access(perms, outside / "a.txt") is None
    # ...the sibling file no longer raises a second dialog...
    assert ex._check_read_access(perms, outside / "b.txt") is None
    assert len(broker.previews) == 1
    # ...and the directory did NOT become writable.
    assert outside not in perms.all_workspace_roots()
    assert perms.find_root_for(outside / "a.txt") is None
    assert perms.find_readable_root_for(outside / "a.txt") == outside
    assert perms.is_read_only_path(outside / "a.txt") is True


def test_a_write_into_that_directory_is_refused(ex, tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "notes"
    outside.mkdir()
    (outside / "a.txt").write_text("data", encoding="utf-8")

    broker = _Approve(answer=True)
    perms = _perms(ws, broker)
    ex._check_read_access(perms, outside / "a.txt")

    err = ex._gate_write_path(
        str(outside / "a.txt"), perms, "write_file",
        {"path": str(outside / "a.txt"), "content": "x"})
    assert err is not None and "READING only" in err

    # ...even with a read baseline, so the stale-write guard is not what
    # is doing the work here.
    perms.read_tracker[str((outside / "a.txt").resolve())] = (
        (outside / "a.txt").stat().st_mtime)
    shell = ex._gate_bash_write_targets(
        f"echo x > {outside / 'a.txt'}", {}, perms)
    assert shell is not None and "READING only" in shell
    assert (outside / "a.txt").read_text(encoding="utf-8") == "data"


def test_approving_a_read_in_home_does_not_make_home_a_root(ex, tmp_path, home):
    """The escalation the constructor exists to prevent, reached at runtime."""
    ws = tmp_path / "ws"
    ws.mkdir()
    (home / "notes.txt").write_text("private", encoding="utf-8")

    broker = _Approve(answer=True)
    perms = _perms(ws, broker)
    assert ex._check_read_access(perms, home / "notes.txt") is None   # read ok

    assert home not in perms.all_workspace_roots()
    assert home not in perms.session_read_dirs
    assert perms.find_root_for(home / "notes.txt") is None


def test_add_extra_dir_applies_the_same_floor_as_the_constructor(tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = A.KitToolPermissions(workspace=ws)
    with pytest.raises(ValueError):
        perms.add_extra_dir(home)
    with pytest.raises(ValueError):
        perms.add_extra_dir("/")
    assert perms.extra_workspace_dirs == ()
    # A project directory is still perfectly addable.
    project = home / "agent_workspace"
    project.mkdir()
    assert perms.add_extra_dir(project) == project.resolve()


@pytest.mark.parametrize("directory,grantable", [
    ("/tmp/project", True),
    ("/etc", False),
    ("/etc/ssl", False),
    ("/", False),
    ("/usr/lib", False),
    ("/var/log", False),
])
def test_a_single_read_never_opens_a_system_directory(directory, grantable):
    assert A._is_grantable_read_dir(directory) is grantable


def test_a_key_directory_is_never_opened(home):
    (home / ".ssh").mkdir()
    assert A._is_grantable_read_dir(home / ".ssh") is False
    assert A._is_grantable_read_dir(home) is False


def test_the_dialog_says_what_is_granted(ex, tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "notes"
    outside.mkdir()
    (outside / "a.txt").write_text("data", encoding="utf-8")

    broker = _Approve(answer=True)
    ex._check_read_access(_perms(ws, broker), outside / "a.txt")
    preview = broker.previews[0]
    assert "READING" in preview
    assert str(outside) in preview
    assert "NO write access" in preview
    # The old wording promised the opposite of what the click did.
    assert "only allows this single read" not in preview


def test_a_refused_read_grants_nothing(ex, tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "notes"
    outside.mkdir()
    (outside / "a.txt").write_text("data", encoding="utf-8")

    broker = _Approve(answer=False)
    perms = _perms(ws, broker)
    assert ex._check_read_access(perms, outside / "a.txt") is not None
    assert perms.session_read_dirs == ()
    assert perms.all_workspace_roots() == (ws.resolve(),)


def test_a_locked_scope_grants_nothing(tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = A.KitToolPermissions(workspace=ws, agent_role="office_agent")
    with pytest.raises(ValueError):
        perms.add_session_read_dir(tmp_path)
    assert perms.all_readable_roots() == (ws.resolve(),)


# ---------------------------------------------------------------------------
# The confirm panel no longer hands out a writable directory
# ---------------------------------------------------------------------------

def test_the_broker_does_not_grant_a_directory_on_a_plain_allow(tmp_path):
    """The broker's session grant was the writable one. The read gate does
    the (read-only) grant itself now, so the broker hands out nothing."""
    import threading
    import time

    widgets = pytest.importorskip("ipywidgets")
    from delfin.agent.kit_confirm import KitConfirmBroker

    broker = KitConfirmBroker(default_timeout_s=5.0)
    calls: list[tuple] = []
    broker.set_persist_callback(lambda kind, value: (calls.append((kind, value)), (True, "ok"))[1])

    f = tmp_path / "data" / "mod.py"
    f.parent.mkdir(parents=True)
    f.write_text("x", encoding="utf-8")

    out = {}
    t = threading.Thread(
        target=lambda: out.update(
            ok=broker.callback("read_file", {"path": str(f)}, "")),
        daemon=True)
    t.start()
    for _ in range(200):
        if broker._pending:
            break
        time.sleep(0.01)
    req = list(broker._pending)[0]
    row = broker._build_request_row(req, widgets)
    row.children[-1].children[0].click()             # "Allow (once)"
    t.join(timeout=5)

    assert out["ok"] is True
    assert calls == []                               # no directory handed out


def test_the_permanent_button_says_it_grants_writes(tmp_path):
    widgets = pytest.importorskip("ipywidgets")
    from delfin.agent.kit_confirm import KitConfirmBroker, _ConfirmRequest

    broker = KitConfirmBroker()
    broker.set_persist_callback(lambda kind, value: (True, "ok"))
    f = tmp_path / "docs" / "readme.md"
    f.parent.mkdir(parents=True)
    f.write_text("x", encoding="utf-8")
    req = _ConfirmRequest(seq=1, tool_name="read_file",
                          args={"path": str(f)}, preview="")
    row = broker._build_request_row(req, widgets)
    status = row.children[1].value
    assert "WRITABLE" in status
