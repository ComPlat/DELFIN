"""`pi_messenger.sh --lang=en`: the phone is told everything in English.

Input: DELFIN_MESSENGER_LANG=en in the gateway's environment. Before, every
text the gateway sends was German only.

Output: every literal the gateway and the pairing pass to ``lang.t`` has an
English entry (a German text without one would stay German silently);
templates are translated before they are filled, so what the agent or the
phone wrote is never changed; a refusal the bridge sends as a whole string
is translated whole. German stays the default, and both languages' words
are understood either way.
"""
from __future__ import annotations

import ast
import asyncio
from pathlib import Path

import pytest

from delfin.messenger import lang, store
from delfin.messenger.gateway import HELP, Rooms

ROOT = Path(__file__).resolve().parents[1] / "delfin" / "messenger"


def _literals_passed_to_t(path: Path):
    for node in ast.walk(ast.parse(path.read_text())):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "t" and node.args):
            arg = node.args[0]
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                yield node.lineno, arg.value


@pytest.mark.parametrize("name", ["gateway.py", "__main__.py"])
def test_every_translated_literal_has_an_english_entry(name):
    missing = [(line, text[:60]) for line, text in _literals_passed_to_t(ROOT / name)
               if text not in lang.EN]
    assert not missing


def test_the_constants_translated_where_used_have_entries():
    from delfin.messenger import gateway
    for text in (gateway.HELP, gateway.UPDATE, gateway.Rooms._CHANGE_HINT,
                 *gateway._PLAIN.values()):
        assert text in lang.EN, text[:60]


def test_english_templates_keep_their_placeholders():
    import string
    for german, english in lang.EN.items():
        fields = lambda s: sorted(f for _, f, _, _ in string.Formatter().parse(s) if f)
        assert fields(german) == fields(english), german[:60]


class _Daemon:
    def __init__(self):
        self.sent = []

    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
        return {}


class _Bridge:
    def ask(self, command):
        if command["cmd"] == "sessions":
            return {"sessions": [{"id": "s1", "name": "t", "state": "idle",
                                  "messengerEnabled": True}]}
        if command.get("action") == "status":
            return {"ok": True, "title": "t", "state": "idle", "messengerEnabled": True}
        return {"error": "Diese Sitzung ist nicht verfügbar. Kein anderer Agent wurde angesprochen."}


@pytest.fixture
def rooms(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    store.put_room("#10", "s1", "t", owner="7", ready=True)
    daemon = _Daemon()
    return Rooms(daemon, _Bridge()), daemon


def _say(rooms, text, mid):
    asyncio.run(rooms.handle({"chat": "#10", "text": text, "sender": "7", "id": mid}))


def test_english_on_request(rooms, monkeypatch):
    room, daemon = rooms
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    _say(room, "/help", "1")
    assert daemon.sent[-1][1] == lang.EN[HELP]
    _say(room, "calculate this", "2")
    assert "This session is not available" in daemon.sent[-1][1]
    _say(room, "notify questions", "3")
    assert daemon.sent[-1][1].startswith("Automatic notifications: fragen.")


def test_german_by_default(rooms, monkeypatch):
    room, daemon = rooms
    monkeypatch.delenv("DELFIN_MESSENGER_LANG", raising=False)
    _say(room, "/help", "1")
    assert daemon.sent[-1][1] == HELP


def test_what_the_agent_wrote_is_never_translated(rooms, monkeypatch):
    room, _ = rooms
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    item = {"kind": "confirm", "action": "shell", "tool": "bash",
            "text": "echo 'Nichts offen.'", "purpose": "Nichts offen.", "options": []}
    shown = room.pending_text([item])
    assert "echo 'Nichts offen.'" in shown
    assert "The agent's statement (unchecked): Nichts offen." in shown
