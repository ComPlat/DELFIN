"""The bridge never takes its own uid as the gateway's.

Input: bridge.configured_peer_uid, which reads DELFIN_MESSENGER_PEER_UID.
The bridge admits any process whose SO_PEERCRED uid equals it. The
provisioning gives the gateway a user of its own; a value set by hand to
the service's own uid (a login node without root) admitted every process
of that user -- the agent's own shell commands among them, which could
then approve the agent's requests through the phone channel.
"""
from __future__ import annotations

import os

from delfin.agent.service import bridge


def test_the_services_own_uid_is_not_a_peer(monkeypatch):
    monkeypatch.setenv(bridge.PEER_UID_ENV, str(os.getuid()))
    assert bridge.configured_peer_uid() is None


def test_another_uid_is_the_peer(monkeypatch):
    other = os.getuid() + 1
    monkeypatch.setenv(bridge.PEER_UID_ENV, str(other))
    assert bridge.configured_peer_uid() == other
