"""A damaged live-registry entry is kept, not overwritten.

Input: the per-subagent live registry entry
``~/.delfin/subagent_running/<sa_id>.json``, truncated the way a kill
mid-write or a full disk truncates it. Output: the damaged file is moved
aside with its bytes intact, the new entry is written, and the new entry
names where the old one went. Semantics: the writer of the entry
(``_running_update``) reads it back to preserve the reserving session's
``owner_session``; a lenient read there built the new file from an empty
answer and the damaged record's owner was gone without a message.

Rejected: refusing the write (the run then has no live entry at all) and
salvaging a truncated prefix (a truncated JSON prefix is not a document).
Both are rejected in ``agent/state_file.py``, which supplies the
mechanism used here.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import subagents as sa


@pytest.fixture
def running_dir(tmp_path, monkeypatch):
    d = tmp_path / "subagent_running"
    d.mkdir()
    monkeypatch.setattr(sa, "_RUNNING_DIR", d)
    return d


def test_a_damaged_running_entry_is_kept_not_overwritten(running_dir):
    f = running_dir / "sa-1.json"
    f.write_text('{"type": "explore", "description": " resear',
                 encoding="utf-8")

    sa._running_update("sa-1", {"type": "explore", "description": "new",
                                "started_at": time.time()})

    aside = sorted(q for q in running_dir.iterdir() if "unreadable" in q.name)
    assert aside, "damaged entry was not moved aside"
    assert aside[0].read_text(encoding="utf-8").startswith('{"type"')
    fresh = json.loads((running_dir / "sa-1.json").read_text(encoding="utf-8"))
    assert fresh["description"] == "new"
    assert "unreadable" in str(fresh.get("_previous_unreadable"))


def test_an_undamaged_entry_still_carries_the_owner(running_dir):
    """The normal path is unchanged: the reserving session's owner survives."""
    f = running_dir / "sa-2.json"
    f.write_text(json.dumps({"type": "explore", "owner_session": "s1"}),
                 encoding="utf-8")

    sa._running_update("sa-2", {"type": "explore", "description": "d",
                                "started_at": time.time()})

    fresh = json.loads(f.read_text(encoding="utf-8"))
    assert fresh["owner_session"] == "s1"
    assert "_previous_unreadable" not in fresh


def test_an_absent_entry_writes_without_a_note(running_dir):
    """No file, nothing to preserve: no quarantine note in the entry."""
    sa._running_update("sa-3", {"type": "explore", "description": "d",
                                "started_at": time.time()})
    fresh = json.loads((running_dir / "sa-3.json").read_text(encoding="utf-8"))
    assert fresh["description"] == "d"
    assert "_previous_unreadable" not in fresh
