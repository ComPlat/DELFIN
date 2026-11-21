"""The question "which host is this" is answered in one place.

Input: nothing (the machine the interpreter runs on). Output: the host
identity ``session_store._this_host`` returns, with the socket's answer and
with an absent socket (its absence is simulated by replacing
``socket.gethostname``, not by skipping the test).

Semantics: the session store's host identity is the product's single answer;
the attention feed keeps no second one (``attention_feed.this_host``, no
caller, is deleted). Rejected alternative: making the feed delegate to the
session store -- that would keep an unread name alive and the feed carries no
host field that could use it (its item kinds name id, kind, group, blocking,
sessionId, createdAt, sessionOpen, sessionSaved).

The invariant test fails at the baseline (the second definition is still
there); the identity tests pass at both arms -- the never-worse guard.
"""
from __future__ import annotations

import socket

from delfin.agent import session_store as ss
from delfin.agent.service import attention_feed


def test_the_attention_feed_keeps_no_second_host_identity():
    assert not hasattr(attention_feed, "this_host")


def test_the_host_identity_answers_when_socket_cannot(monkeypatch):
    monkeypatch.setattr(socket, "gethostname", lambda: "node-a")
    assert ss._this_host() == "node-a"

    def boom():
        raise OSError("no hostname here")

    monkeypatch.setattr(socket, "gethostname", boom)
    assert ss._this_host() == ""