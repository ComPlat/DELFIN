"""Every chat names the machine it comes from.

One phone may pair with several installations (a Pi, a server, a login
node). Each is its own SimpleX contact; without a name in the chats they
were three times "DELFIN" and three forums "DELFIN Forum".

Input: the gateway's group names (delfin/messenger/gateway.py) and the
provisioner (scripts/pi/pi_messenger.sh).

- ``DELFIN_MESSENGER_LABEL`` (set by the provisioner: ``--name=LABEL``,
  else the short host name, reduced to [A-Za-z0-9_-]) is part of the bot's
  display name, ``DELFIN-<label>``, of every session chat and of the forum.
- ``name NAME`` renames a session; its chat keeps the machine label.
- Without a label the names are as before.
"""
from __future__ import annotations

import asyncio

import pytest

from delfin.messenger import store
from delfin.messenger.gateway import Rooms

from test_a_session_can_be_started_from_the_phone import _Bridge, _Daemon
from test_the_messenger_provisioning_covers_every_pinned_server import run_script


class _Renaming(_Daemon):
    def __init__(self):
        super().__init__()
        self.renamed = []

    async def rename_group(self, chat, name):
        self.renamed.append((chat, name))


@pytest.fixture
def rooms(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")

    def make(label):
        if label is None:
            monkeypatch.delenv("DELFIN_MESSENGER_LABEL", raising=False)
        else:
            monkeypatch.setenv("DELFIN_MESSENGER_LABEL", label)
        daemon = _Renaming()
        return Rooms(daemon, _Bridge()), daemon
    return make


def test_session_chat_and_forum_carry_the_label(rooms):
    room, daemon = rooms("Pi")
    asyncio.run(room.ensure_room("s-new"))
    asyncio.run(room.ensure_room())
    assert daemon.groups[0].startswith("DELFIN-Pi Butan ")
    assert daemon.groups[1] == "DELFIN Forum Pi"


def test_without_a_label_the_names_are_as_before(rooms):
    room, daemon = rooms(None)
    asyncio.run(room.ensure_room("s-new"))
    asyncio.run(room.ensure_room())
    assert daemon.groups[0].startswith("DELFIN Butan ")
    assert daemon.groups[1] == "DELFIN Forum"


def test_a_renamed_session_keeps_the_machine_in_its_chat_name(rooms):
    room, daemon = rooms("Pi")
    chat = asyncio.run(room.ensure_room("s-new"))
    out = asyncio.run(room.dispatch("name Konformere", store.rooms()[chat], chat=chat))
    assert "Konformere" in out
    assert daemon.renamed[-1][1].startswith("DELFIN-Pi Konformere")


def _units(tmp_path):
    return {p.name: p.read_text() for p in (tmp_path / "units").glob("*.service")}


def test_the_provisioner_names_bot_and_gateway(tmp_path):
    assert run_script(tmp_path, args=["--name=Pi 5!"]).returncode == 0
    units = _units(tmp_path)
    assert "--create-bot-display-name DELFIN-Pi5" in units["simplex-chat-daemon.service"]
    assert "Environment=DELFIN_MESSENGER_LABEL=Pi5" in units["delfin-messenger-gateway.service"]


def test_the_host_name_is_the_default_label(tmp_path):
    assert run_script(tmp_path, env_extra={"MESSENGER_HOSTNAME": "login-node.kit.edu"}).returncode == 0
    units = _units(tmp_path)
    assert "--create-bot-display-name DELFIN-login-node" in units["simplex-chat-daemon.service"]
    assert "Environment=DELFIN_MESSENGER_LABEL=login-node" in units["delfin-messenger-gateway.service"]
