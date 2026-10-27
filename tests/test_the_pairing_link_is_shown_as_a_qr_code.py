"""The pairing invitation is shown as a QR code when qrencode is installed.

Input: ``python -m delfin.messenger pair`` (delfin/messenger/__main__.py).
The invitation is a long simplex: link that had to be carried from an SSH
terminal to the phone by hand. With ``qrencode`` on PATH it is also drawn
in the terminal (``-t ansiutf8``), for the SimpleX app's scanner; without
it the link alone is printed, as before. The link goes to qrencode as an
argument, never through a shell.
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


def test_a_qr_code_is_drawn_when_qrencode_is_there(tmp_path, monkeypatch, capsys):
    log = _pair(tmp_path, monkeypatch, with_qrencode=True)
    assert log.read_text().splitlines() == [
        "-t", "ansiutf8", "simplex:/invitation#/?v=2&smp=x;y"]
    assert "simplex:/invitation#/?v=2&smp=x;y" in capsys.readouterr().out


def test_without_qrencode_the_link_alone_is_printed(tmp_path, monkeypatch, capsys):
    log = _pair(tmp_path, monkeypatch, with_qrencode=False)
    assert not os.path.exists(log)
    out = capsys.readouterr().out
    assert "simplex:/invitation#/?v=2&smp=x;y" in out
    assert "qrencode" in out
