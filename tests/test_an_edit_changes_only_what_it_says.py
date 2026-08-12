"""Editing one line rewrote the whole file, and a failed write emptied it.

Every write tool used ``Path.read_text`` / ``Path.write_text``. Three
consequences, all reproduced before this was written:

* A CRLF file came back LF. Changing one word in a Windows-authored
  ``.csv`` rewrote **every** line ending, while the diff shown to the
  user was the single line that was meant to change::

      b"alpha\\r\\nbeta\\r\\ngamma\\r\\n"  ->  b"alpha\\nBETA\\ngamma\\n"

* ``errors="replace"`` made every non-UTF-8 byte a U+FFFD, and the next
  write made it permanent — on a line the diff displayed as unchanged
  context::

      b"caf\\xe9 = 1\\nx = 2\\n"  ->  b"caf\\xef\\xbf\\xbd = 1\\nx = 3\\n"

* ``write_text`` truncates first. A 2000-byte file whose write raised
  partway through was left at **zero bytes**, and the pre-image — held
  only in a local variable — was discarded by the error return, so the
  undo journal had nothing either. The one case where the journal is
  the last copy is the case it did not cover.

The internal stores already did this properly: the change journal
writes through a temporary file with ``os.replace`` and ``newline=""``.
The user's own files did not get the same care.

What is asserted here is the whole contract in one place: a write
changes the bytes it says it changes and nothing else, and a write that
fails changes nothing at all.
"""

from __future__ import annotations

import os
import stat

import pytest

from delfin.agent import text_files as TF


# ---------------------------------------------------------------------------
# Line endings survive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,newline", [
    (b"alpha\r\nbeta\r\ngamma\r\n", "\r\n"),
    (b"alpha\nbeta\ngamma\n", "\n"),
    (b"alpha\rbeta\rgamma\r", "\r"),
])
def test_the_file_keeps_its_own_line_endings(tmp_path, raw, newline):
    p = tmp_path / "data.csv"
    p.write_bytes(raw)
    shape = TF.read_text_file(p)
    assert shape.newline == newline
    TF.write_text_file(p, shape.text.replace("beta", "BETA"), like=shape)
    assert p.read_bytes() == raw.replace(b"beta", b"BETA")


def test_the_text_handed_out_is_normalised(tmp_path):
    """The model matches its edit strings against LF, whatever the file is."""
    p = tmp_path / "data.csv"
    p.write_bytes(b"alpha\r\nbeta\r\n")
    assert TF.read_text_file(p).text == "alpha\nbeta\n"


def test_a_new_file_is_written_as_given(tmp_path):
    p = tmp_path / "new.txt"
    TF.write_text_file(p, "a\nb\n", like=None)
    assert p.read_bytes() == b"a\nb\n"


# ---------------------------------------------------------------------------
# Bytes that are not UTF-8 survive
# ---------------------------------------------------------------------------

def test_a_latin1_byte_is_not_replaced(tmp_path):
    p = tmp_path / "paper.tex"
    p.write_bytes(b"caf\xe9 = 1\nx = 2\n")
    shape = TF.read_text_file(p)
    assert "café" in shape.text, "the byte has to survive the read"
    TF.write_text_file(p, shape.text.replace("x = 2", "x = 3"), like=shape)
    assert p.read_bytes() == b"caf\xe9 = 1\nx = 3\n"


def test_utf8_is_preferred_over_a_lookalike(tmp_path):
    p = tmp_path / "u.txt"
    p.write_bytes("café\n".encode("utf-8"))
    shape = TF.read_text_file(p)
    assert shape.encoding == "utf-8"
    TF.write_text_file(p, shape.text, like=shape)
    assert p.read_bytes() == "café\n".encode("utf-8")


def test_a_byte_order_mark_survives(tmp_path):
    p = tmp_path / "excel-export.csv"
    p.write_bytes(b"\xef\xbb\xbfBeleg;Betrag\r\nR-001;12,50\r\n")
    shape = TF.read_text_file(p)
    assert shape.bom is True
    assert shape.text.startswith("Beleg"), "the model must not see the mark"
    TF.write_text_file(p, shape.text.replace("12,50", "13,50"), like=shape)
    assert p.read_bytes().startswith(b"\xef\xbb\xbf")
    assert b"13,50" in p.read_bytes()


def test_changing_the_encoding_is_reported(tmp_path):
    """Silently re-encoding a file is a change nobody asked for."""
    p = tmp_path / "old.tex"
    p.write_bytes(b"caf\xe9\n")
    shape = TF.read_text_file(p)
    notes = TF.write_text_file(p, shape.text + "λ\n", like=shape)
    assert notes and "UTF-8" in notes[0]
    assert p.read_bytes().decode("utf-8") == "café\nλ\n"


# ---------------------------------------------------------------------------
# A write that fails changes nothing
# ---------------------------------------------------------------------------

def test_a_failed_write_leaves_the_file_intact(tmp_path, monkeypatch):
    p = tmp_path / "important.txt"
    p.write_bytes(b"D" * 2000)

    real = os.fdopen

    def _boom(fd, *a, **kw):
        fh = real(fd, *a, **kw)

        class _Failing:
            def __enter__(self_inner):
                return self_inner

            def __exit__(self_inner, *exc):
                fh.close()
                return False

            def write(self_inner, _data):
                raise OSError("no space left on device")

        return _Failing()

    monkeypatch.setattr(os, "fdopen", _boom)
    with pytest.raises(OSError):
        TF.write_text_file(p, "replacement\n", like=TF.read_text_file(p))
    assert p.read_bytes() == b"D" * 2000


def test_a_failed_write_leaves_no_litter(tmp_path, monkeypatch):
    p = tmp_path / "important.txt"
    p.write_bytes(b"x\n")
    before = {q.name for q in tmp_path.iterdir()}
    monkeypatch.setattr(
        os, "replace",
        lambda *_a, **_k: (_ for _ in ()).throw(OSError("cross-device")))
    with pytest.raises(OSError):
        TF.write_text_file(p, "y\n")
    assert {q.name for q in tmp_path.iterdir()} == before


def test_the_replacement_is_atomic(tmp_path):
    """Never a moment where the file exists and is half written."""
    p = tmp_path / "f.txt"
    p.write_bytes(b"old\n")
    seen: list[bytes] = []
    real_replace = os.replace

    def _watch(src, dst):
        seen.append(p.read_bytes())
        return real_replace(src, dst)

    import unittest.mock as _mock
    with _mock.patch.object(os, "replace", _watch):
        TF.write_text_file(p, "new content that is longer\n")
    assert seen == [b"old\n"], "the old content was still whole until the swap"
    assert p.read_bytes() == b"new content that is longer\n"


def test_the_permission_bits_are_carried_over(tmp_path):
    p = tmp_path / "secret.env"
    p.write_bytes(b"TOKEN=1\n")
    p.chmod(0o600)
    TF.write_text_file(p, "TOKEN=2\n", like=TF.read_text_file(p))
    assert stat.S_IMODE(p.stat().st_mode) == 0o600


# ---------------------------------------------------------------------------
# ...and the write tools actually use it
# ---------------------------------------------------------------------------

def test_the_write_tools_no_longer_call_write_text():
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    for fn in ("_execute_write_file", "_execute_edit_file",
               "_execute_multi_edit"):
        i = src.index(f"def {fn}(")
        j = src.index("\n    def ", i + 10)
        body = src[i:j]
        assert ".write_text(" not in body, (
            f"{fn} still writes through Path.write_text, which translates "
            f"line endings and truncates before it writes")
        assert "_text_files.write_text_file(" in body
