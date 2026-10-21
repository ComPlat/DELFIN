"""KG1 empty-send red control (reviewer finding #3). Not committed as a
green test: it is RED until the send() refusal is in. The operator or a
granted write applies the fix; then this file proves it and joins the
next commit as the control.

delfin/agent/session_messages.py must be committed at cb40378f..f378d803
state for this to be the right red baseline.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_messages as msgs


@pytest.fixture
def inbox(tmp_path, monkeypatch):
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "session_inbox")
    return tmp_path / "session_inbox"


def test_send_with_an_empty_address_is_refused(inbox):
    with pytest.raises(ValueError):
        msgs.send("", "nowhere")


def test_send_with_a_whitespace_address_is_refused(inbox):
    with pytest.raises(ValueError):
        msgs.send("   ", "nowhere")


def test_no_stemless_file_is_created(inbox):
    with pytest.raises(ValueError):
        msgs.send("", "nowhere")
    assert inbox.exists() is False or not [
        p for p in inbox.glob("*.jsonl")
        if p.name == ".jsonl"
    ], "no stem-less mailbox file may be created"
