"""QS smoke test for package T3 coordination review tooling.

This is a QS-owned sanity test (not part of the T3 build itself): it verifies
that the module under review imports and exposes the interface the T3 phases
build on, and that the gate toolchain can run a test in this worktree. It pins
nothing transient -- only that `session_messages` is importable and its inbox
key function produces a safe file name.
"""


def test_session_messages_module_imports_and_inbox_key():
    from delfin.agent import session_messages as msgs

    # The inbox path is derived from a key the same way for every session;
    # an unsanitised char must not escape the key prefix (security-relevant
    # for the ownership of a message file).
    p = msgs._inbox("runde2-s8")
    assert p.name.endswith(".jsonl")
    assert not any(c in "/\\" for c in p.stem)


def test_send_message_carries_the_documented_fields(tmp_path, monkeypatch):
    from delfin.agent import session_messages as msgs
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "inbox")

    msg = msgs.send("peer", "hello", from_key="me", from_title="sender")
    assert msg["to"] == "peer"
    assert msg["from"] == "me"
    assert msg["text"] == "hello"
    assert msg["sent_at"] > 0

    taken = msgs.take("peer")
    assert len(taken) == 1
    assert taken[0]["text"] == "hello"
