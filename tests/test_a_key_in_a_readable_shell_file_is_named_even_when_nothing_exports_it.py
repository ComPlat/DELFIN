"""A key in a shell file others can read is named, even unexported.

The start-up hygiene check asked the ENVIRONMENT which provider keys are
present and then searched the shell files for those names. A key that
never reaches the environment was therefore never searched for -- and
that is the worse case, because it sits in the file permanently instead
of only for the life of a shell.

Seen on a cluster, 2026-09-20: a group-readable ``~/.bashrc`` carried

    alias codex-kit='OPENAI_API_KEY="sk-..." OPENAI_BASE_URL=... codex'

for months. ``exported_in_shell_files`` finds that line in one call; the
automatic remediation never asked it to, because ``OPENAI_API_KEY`` is
only set for the duration of the alias and so is absent from every
environment DELFIN ever inspects.

What the check may and may not do here is not symmetric:

  it names the file and line     so the user can act
  it says who else can read it   that is the risk, not the mere presence
  it does NOT edit the line      the value is not in the store, and a
                                 tool that removes the only copy of a
                                 key is worse than the leak
  a file only the owner can      read is not reported: the key is where
                                 the user put it, and nobody else sees it
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import credentials as cred


ALIAS_LINE = (
    "alias codex-kit='OPENAI_API_KEY=\"sk-live-abcdef123456\" "
    "OPENAI_BASE_URL=\"https://example.invalid/api/v1\" codex'\n"
)


@pytest.fixture
def home(a_path_others_can_walk):
    # Genuinely reachable from the root: "others can read it" is
    # the premise of these tests, and pytest's own base directory
    # is 0700, so under it nobody can read anything. Before the
    # check asked about the PATH, the file's mode bits made that
    # premise true by fiat.
    tmp_path = a_path_others_can_walk
    rc = tmp_path / ".bashrc"
    rc.write_text("# my shell\n" + ALIAS_LINE + "echo hello\n", encoding="utf-8")
    os.chmod(rc, 0o644)
    return tmp_path


def _at_risk(home, **kw):
    return cred.keys_in_readable_shell_files(home=home, **kw)


def test_the_unexported_key_is_found(home):
    rows = _at_risk(home)
    assert [r["name"] for r in rows] == ["OPENAI_API_KEY"]


def test_it_names_the_file_and_the_line(home):
    row = _at_risk(home)[0]
    assert row["file"].endswith(".bashrc")
    assert row["line_no"] == 2


def test_it_says_who_else_can_read_it(home):
    row = _at_risk(home)[0]
    assert row["mode"] == "0644"
    assert row["readable_by"] == "group, others"


def test_a_file_only_the_owner_can_read_is_not_reported(home):
    os.chmod(home / ".bashrc", 0o600)
    assert _at_risk(home) == []


def test_group_readable_alone_is_enough(home):
    os.chmod(home / ".bashrc", 0o640)
    row = _at_risk(home)[0]
    assert row["readable_by"] == "group"


def test_the_line_is_left_alone(home):
    _at_risk(home)
    assert ALIAS_LINE.strip() in (home / ".bashrc").read_text()


def test_a_commented_line_is_not_reported(home):
    rc = home / ".bashrc"
    rc.write_text("# " + ALIAS_LINE, encoding="utf-8")
    os.chmod(rc, 0o644)
    assert _at_risk(home) == []


def test_an_empty_home_is_silent(tmp_path):
    assert _at_risk(tmp_path) == []


def test_the_reported_line_carries_the_path_and_the_advice(home, monkeypatch):
    from delfin.agent import process_guard

    monkeypatch.setattr(cred, "_DEFAULT_PATH",
                        home / ".delfin" / "credentials.json")
    monkeypatch.setattr(os.path, "expanduser",
                        lambda p: p.replace("~", str(home), 1))
    line = process_guard.put_exported_keys_away((), home=home)
    assert ".bashrc" in line
    assert "OPENAI_API_KEY" in line
    assert "sk-live-abcdef123456" not in line, "the value never appears"
