"""Speech recognition runs in an interpreter of its own, isolated.

Input: delfin.messenger.media.transcribe. faster-whisper 1.2.1 needs numpy
2.x; DELFIN in the gateway's venv keeps numpy below 2, so the two cannot
share one environment. The transcriber is a separate process already; it
now starts the interpreter named by DELFIN_WHISPER_PYTHON (default: this
one) with -I and a self-contained script, so nothing from DELFIN or the
current directory is imported there.
"""
from __future__ import annotations

import json
import sys

import pytest

from delfin.messenger import media


def _fake_python(tmp_path, transcript="Ja."):
    log = tmp_path / "argv.json"
    fake = tmp_path / "python-whisper"
    fake.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\n"
        f"json.dump(sys.argv[1:], open({str(log)!r}, 'w'))\n"
        f"print(json.dumps({transcript!r}))\n")
    fake.chmod(0o755)
    return fake, log


def test_the_named_interpreter_runs_isolated(tmp_path, monkeypatch):
    model = tmp_path / "model"
    model.mkdir()
    fake, log = _fake_python(tmp_path)
    monkeypatch.setenv("DELFIN_WHISPER_MODEL", str(model))
    monkeypatch.setenv("DELFIN_WHISPER_PYTHON", str(fake))
    assert media.transcribe(tmp_path / "note.ogg") == "Ja."
    argv = json.loads(log.read_text())
    assert argv[0] == "-I" and argv[1] == "-c"
    assert "delfin" not in argv[2]
    assert argv[3:] == [str(tmp_path / "note.ogg"), str(model)]


def test_a_missing_interpreter_is_a_plain_error(tmp_path, monkeypatch):
    model = tmp_path / "model"
    model.mkdir()
    monkeypatch.setenv("DELFIN_WHISPER_MODEL", str(model))
    monkeypatch.setenv("DELFIN_WHISPER_PYTHON", str(tmp_path / "absent"))
    with pytest.raises(ValueError):
        media.transcribe(tmp_path / "note.ogg")
