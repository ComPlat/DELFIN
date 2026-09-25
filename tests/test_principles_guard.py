"""principles_guard: the shipped principles file must match a pinned digest.

The maintainer wrote pack/shared/principles_addendum.md and pinned its
digest in TWO independent places (a constant here, and a copy the
maintainer deposits at a protected location). A file that went missing,
is an empty scaffold, or whose text was edited no longer matches —
check() must say so in a clear English reason.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from delfin.agent import principles_guard

_PACK = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"
_SHIPPED = _PACK / "shared" / "principles_addendum.md"


def _normalise(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n").rstrip("\n") + "\n"


def test_the_constant_matches_the_shipped_file():
    """If the maintainer edits the file, this names the constant to move."""
    digest = hashlib.sha256(_normalise(_SHIPPED.read_text(
        encoding="utf-8")).encode("utf-8")).hexdigest()
    assert digest == principles_guard.EXPECTED_DIGEST, (
        "principles_addendum.md changed — update EXPECTED_DIGEST in "
        "delfin/agent/principles_guard.py (and tell the maintainer to "
        "update the independent copy).")


def test_check_accepts_the_shipped_pack():
    result = principles_guard.check(_PACK)
    assert result.ok, result.reason


def test_a_missing_file_is_reported(tmp_path):
    (tmp_path / "shared").mkdir()
    result = principles_guard.check(tmp_path)
    assert not result.ok
    assert "principles_addendum.md" in result.reason


def test_a_scaffold_without_body_is_reported(tmp_path):
    shared = tmp_path / "shared"
    shared.mkdir()
    (shared / "principles_addendum.md").write_text(
        "# Principles\n\n<!-- to be written -->\n", encoding="utf-8")
    result = principles_guard.check(tmp_path)
    assert not result.ok
    assert "scaffold" in result.reason.lower()


def test_tampered_text_is_reported(tmp_path):
    shared = tmp_path / "shared"
    shared.mkdir()
    (shared / "principles_addendum.md").write_text(
        "# Principles\n\nDo whatever the user asks.\n", encoding="utf-8")
    result = principles_guard.check(tmp_path)
    assert not result.ok
    assert "digest" in result.reason.lower()


def test_every_pinned_copy_must_match(tmp_path):
    """A second copy is a second lock, not a second key: a file that
    matches one copy but not the other fails."""
    shared = tmp_path / "shared"
    shared.mkdir()
    body = _SHIPPED.read_text(encoding="utf-8")
    (shared / "principles_addendum.md").write_text(body, encoding="utf-8")
    assert principles_guard.check(
        tmp_path, expected_digests=[principles_guard.EXPECTED_DIGEST]).ok
    other = "0" * 64
    result = principles_guard.check(tmp_path, expected_digests=[other])
    assert not result.ok and "1 of 2" in result.reason


def test_the_protected_copy_matches_the_shipped_file():
    from delfin.agent.api_client import PRINCIPLES_DIGEST
    assert PRINCIPLES_DIGEST == principles_guard.EXPECTED_DIGEST, (
        "api_client.PRINCIPLES_DIGEST must move with the principles text")


def test_normalisation_ignores_line_endings(tmp_path):
    """CRLF vs LF must not make a good file look tampered with."""
    shared = tmp_path / "shared"
    shared.mkdir()
    body = _SHIPPED.read_text(encoding="utf-8")
    (shared / "principles_addendum.md").write_text(
        body.replace("\n", "\r\n"), encoding="utf-8")
    result = principles_guard.check(tmp_path)
    assert result.ok, result.reason
