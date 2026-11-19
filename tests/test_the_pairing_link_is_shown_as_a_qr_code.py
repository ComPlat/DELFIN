"""The pairing invitation is shown as a QR code, on every host.

Input: ``python -m delfin.messenger pair`` (delfin/messenger/__main__.py).
The invitation is a long simplex: link that had to be carried from an SSH
terminal to the phone by hand. It is drawn in the terminal for the SimpleX
app's scanner by delfin/messenger/qr.py -- no tool needed, and no program
on PATH is ever run for it (before: qrencode when installed, else the link
alone) -- and the link is printed after it.
"""
from __future__ import annotations

import os

from delfin.messenger import __main__ as cli


class _Daemon:
    invitation_id = "conn-1"

    def __init__(self, *args, **kwargs):
        pass

    async def __aenter__(self):
        return self

    async def __aexit__(self, *args):
        pass

    async def subscribe(self):
        pass

    async def pending_invitation(self):
        return "simplex:/invitation#/?v=2&smp=x;y"

    async def events(self):
        yield {"type": "contactConnected",
               "contact": {"contactId": 9, "activeConn": {"connId": "conn-1"}}}

    async def contact_code(self, contact_id):
        return None

    async def send_text(self, contact_id, text):
        pass


def _pair(tmp_path, monkeypatch, with_qrencode):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "pairing"))
    monkeypatch.setattr(cli.simplex, "Daemon", _Daemon)
    monkeypatch.setattr(cli.simplex, "new_contact_id",
                        lambda e: str(e["contact"]["contactId"]))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "qrencode.args"
    if with_qrencode:
        tool = bin_dir / "qrencode"
        tool.write_text(f'#!/bin/sh\nfor a in "$@"; do printf "%s\\n" "$a"; done > "{log}"\n')
        tool.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_dir))
    assert cli._pair() == 0
    return log


def test_the_code_is_drawn_without_any_tool(tmp_path, monkeypatch, capsys):
    _pair(tmp_path, monkeypatch, with_qrencode=False)
    out = capsys.readouterr().out
    assert "▄" in out or "▀" in out
    assert out.index("█") < out.index("simplex:/invitation#/?v=2&smp=x;y"), "the link comes last"
    assert "qrencode" not in out


def test_a_qrencode_on_path_is_never_run(tmp_path, monkeypatch, capsys):
    log = _pair(tmp_path, monkeypatch, with_qrencode=True)
    assert not os.path.exists(log)
    assert "simplex:/invitation#/?v=2&smp=x;y" in capsys.readouterr().out
