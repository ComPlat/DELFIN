"""The messenger's directory in a workspace is read-only to the agent,
except the one place a file is offered for download.

Field report (Pi, an owner TUI session): the model, not told it was on
the phone, wrote ``.delfin-messenger/<hash>/outbox/forum.txt`` by hand.
``.delfin-messenger/<session hash>/`` holds the phone's uploads
(``inbox``: data the agent reads) and the session's ``outbox`` (what the
phone can fetch with ``datei NAME``; ``forum.txt`` there is the explicit
note to other agents).

Output, through the real write gate (write_file end to end, and the
gate every shell write target passes):

* a write into ``inbox``, into another session's directory, into a
  sub-directory of the outbox or to the directory itself: refused, with
  the outbox path of THIS session as the pointer;
* control: a plain file in this session's own outbox is written;
* control: reading a file the phone sent (inbox) still works.
"""
from __future__ import annotations

import hashlib
import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor

SID = "s-office"
OWN = hashlib.sha256(SID.encode()).hexdigest()[:24]
OTHER = hashlib.sha256(b"s-other").hexdigest()[:24]


@pytest.fixture
def ws(tmp_path):
    for h in (OWN, OTHER):
        for box in ("inbox", "outbox"):
            (tmp_path / ".delfin-messenger" / h / box).mkdir(parents=True)
    (tmp_path / ".delfin-messenger" / OWN / "inbox" / "photo.txt").write_text("from the phone")
    perms = KitToolPermissions(workspace=str(tmp_path), mode="acceptEdits")
    perms.task_session_id = SID
    perms.confirm_callback = lambda *a: True
    return tmp_path, perms


def _write(perms, path):
    """write_file's result: a JSON error, or the plain success text."""
    out = _doc_executor.execute("write_file", {"path": path, "content": "x\n"}, perms)
    try:
        return json.loads(out)
    except ValueError:
        return {"text": out}


@pytest.mark.parametrize("path", [
    f".delfin-messenger/{OWN}/inbox/photo.txt",
    f".delfin-messenger/{OWN}/inbox/new.txt",
    f".delfin-messenger/{OTHER}/outbox/forum.txt",
    f".delfin-messenger/{OWN}/outbox/sub/x.txt",
    f".delfin-messenger/{OWN}/notes.txt",
    ".delfin-messenger/x.txt",
])
def test_a_write_into_the_internals_is_refused_with_the_outbox_named(ws, path):
    root, perms = ws
    out = _write(perms, path)
    assert "error" in out, out
    assert f".delfin-messenger/{OWN}/outbox" in out["error"]
    assert (root / path).read_text() == "from the phone" if path.endswith("photo.txt") \
        else not (root / path).exists()


def test_every_shell_write_target_meets_the_same_refusal(ws):
    _root, perms = ws
    why = _doc_executor._gate_write_path(
        f".delfin-messenger/{OWN}/inbox/new.txt", perms, "bash",
        {"command": f"echo x > .delfin-messenger/{OWN}/inbox/new.txt"})
    assert why and "outbox" in why


def test_a_file_offered_in_the_own_outbox_is_written(ws):
    root, perms = ws
    out = _write(perms, f".delfin-messenger/{OWN}/outbox/report.txt")
    assert "error" not in out, out
    assert (root / ".delfin-messenger" / OWN / "outbox" / "report.txt").read_text() == "x\n"


def test_a_file_the_phone_sent_is_read(ws):
    _root, perms = ws
    out = _doc_executor.execute("read_file", {"path": f".delfin-messenger/{OWN}/inbox/photo.txt"}, perms)
    assert "from the phone" in out
