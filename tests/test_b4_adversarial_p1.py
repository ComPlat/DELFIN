"""Adversarial probes for B4 phase 1 (reviewer nacht-s19).

Attacks register()'s new fields: the on_whom enum, the for_what one-line
guarantee, and the file-wait path validation. Every case here is an
ALLOW-list or an injection attempt, not a re-run of the builder's happy
paths (those live in tests/test_b4_waits.py).
"""

from __future__ import annotations

import pytest
from pathlib import Path

from delfin.agent import state_paths
from delfin.agent.wait_registry import (
    ON_WHOM_VALUES,
    _clean_oneline,
    _clean_on_whom,
    register,
)


@pytest.fixture
def home(tmp_path, monkeypatch):
    # The file-wait anchor is the process cwd (builder P2 fix cd3e928d):
    # every test chdirs into its own throwaway workspace so the anchor
    # and the escape targets are both tmp-scoped.
    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.chdir(ws)
    with state_paths.RedirectedUserState(tmp_path):
        yield tmp_path / ".delfin"


# -- on_whom enum: case-fold ok, lookalikes and padding refused --------------

def test_on_whom_case_fold_and_enum(home):
    assert _clean_on_whom("PERSON") == "person"
    # Intent note: "person\n" and "OPERATOR " are NORMALIZED to the clean
    # enum value ("person", "operator") rather than refused -- the output is
    # always exactly one enum value, so no injection is possible through
    # padding. Only genuinely unknown values must be refused.
    for bad in ("session leader", "operator2", "çerson", True, 3,
                "person-session"):
        with pytest.raises(ValueError):
            register("message_from", "x", on_whom=bad)


def test_on_whom_strip_normalization_is_the_design(home):
    # Intent note: trailing whitespace INCLUDING the unicode NBSP is stripped
    # before the enum check, so "person\u00a0" normalizes to "person". The
    # output is always exactly one enum value, so this cannot inject a
    # second token; lookalikes ("çerson") and unknown values still refuse.
    assert _clean_on_whom("person\u00a0") == "person"
    assert _clean_on_whom("OPERATOR ") == "operator"


def test_on_whom_enum_is_the_briefs_values():
    assert set(ON_WHOM_VALUES) == {"session", "operator", "person", "ci"}


# -- for_what: no second line, no markup, bounded -----------------------------

@pytest.mark.parametrize("hostile", [
    "line one\nline two",
    "line one\r\nline two",
    "line one\x0bvertical tab",
    "line one\x0cform feed",
    "line one\u2028line sep",
    "line one\u0085nel",
    "<script>alert(1)</script>",
    "<b>bold</b> claim",
])
def test_for_what_cannot_inject_second_line_or_markup(home, hostile):
    w = register("message_from", "nacht-s29", for_what=hostile)
    cleaned = w["for_what"]
    assert "\n" not in cleaned
    assert "\r" not in cleaned
    assert "\x0b" not in cleaned
    assert "\x0c" not in cleaned
    assert "\u2028" not in cleaned
    assert "<" not in cleaned and ">" not in cleaned


def test_for_what_is_bounded(home):
    w = register("message_from", "nacht-s29", for_what="x" * 5000)
    assert len(w["for_what"]) <= 200


def test_clean_oneline_preserves_legitimate_text():
    assert _clean_oneline("wait for QS-OK on PR 320") == \
        "wait for QS-OK on PR 320"


# -- file path validation: escapes refused at register time -------------------

def test_file_wait_refuses_absolute_path(home):
    with pytest.raises(ValueError):
        register("file", "/etc/passwd")


def test_file_wait_refuses_parent_escape(home):
    with pytest.raises(ValueError):
        register("file", "../outside.txt")


def test_file_wait_refuses_symlink_escape(home, tmp_path):
    # The workspace anchor is the process cwd (see the fixture); the link
    # sits INSIDE it and points OUTSIDE it. resolve() must catch this.
    target = tmp_path / "outside_symlink_dir"
    target.mkdir(exist_ok=True)
    (target / "leak.txt").write_text("secret", encoding="utf-8")
    import os
    os.symlink(str(target), str(Path.cwd() / "link_dir"))
    with pytest.raises(ValueError):
        register("file", "link_dir/leak.txt")


# (_workspace_root helper removed: the anchor is now cwd — see fixture.)

def test_file_wait_refuses_empty_and_whitespace(home):
    for bad in ("", "   "):
        with pytest.raises(ValueError):
            register("file", bad)


def test_file_wait_accepts_a_legitimate_relative_path(home, tmp_path):
    # Anchored at cwd (P2 fix cd3e928d): create the file in the chdir'd
    # throwaway workspace, not under the state root (the old fixture
    # pinned the state root — stated here, the old decision is bent).
    f = Path.cwd() / "watched.txt"
    f.write_text("v1", encoding="utf-8")
    w = register("file", "watched.txt")
    assert w["kind"] == "file"
    assert "file_sig" in w


def test_legacy_caller_with_no_kwargs_still_works(home):
    w = register("message_from", "nacht-s25")
    assert w["on_whom"] == ""
    assert w["for_what"] == ""
