"""With no phone paired the gateway waits; it never exits with an error for that.

Field report (Pi): with nothing paired the gateway exited 1 ("not paired;
run ... pair first") and systemd restarted it every 5 s, 140+ times.

Output:

* ``python -m delfin.messenger run`` with an empty pairing store serves (it
  must: people pair through it), says once that no phone is paired yet and
  how to pair, and returns 0 when its daemon connection ends;
* the gateway looks at the pairing store every 30 s and says once when the
  first phone is paired (serving) -- not on every look (control: nothing
  is said again while the state stays the same).
"""
from __future__ import annotations

from delfin.messenger import __main__ as cli
from delfin.messenger import store
from delfin.messenger.gateway import Rooms


class _Daemon:
    def __init__(self, *args, **kwargs):
        pass

    async def __aenter__(self):
        return self

    async def __aexit__(self, *args):
        pass

    async def subscribe(self):
        return {}

    async def events(self):
        if False:
            yield {}


def test_an_empty_store_waits_and_never_exits_non_zero(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "pairing"))
    monkeypatch.setattr(cli.simplex, "Daemon", _Daemon)

    class _Bridge:
        def ask(self, command):
            return {"ok": True, "requests": []}

    monkeypatch.setattr(cli, "BridgeClient", lambda *a, **k: _Bridge(), raising=False)
    assert cli._run() == 0
    out = capsys.readouterr()
    text = out.out + out.err
    assert text.count("no phone paired yet") == 1
    assert "not paired" not in text


def test_the_store_is_looked_at_every_30_s_and_a_first_phone_is_announced(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "pairing"))
    rooms = Rooms(_Daemon(), bridge=object())
    said = []
    assert rooms.watch_pairing(now=0.0, say=said.append) is True
    assert rooms.watch_pairing(now=10.0, say=said.append) is False      # within 30 s: no look
    store.save("7")
    assert rooms.watch_pairing(now=20.0, say=said.append) is False
    assert rooms.watch_pairing(now=31.0, say=said.append) is True       # looked: paired now
    assert rooms.watch_pairing(now=62.0, say=said.append) is True
    assert len(said) == 2
    assert "no phone paired yet" in said[0] and "serving" in said[1]
