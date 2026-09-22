"""A shell line that only READS a key is not an export of it.

Seen on a cluster, 2026-09-20: the key had been moved out of a
group-readable ``~/.bashrc`` into a 0600 file, and the user's alias was
rewritten to read it from there:

    alias codex-kit='OPENAI_API_KEY=$KIT_TOOLBOX_API_KEY OPENAI_BASE_URL=... codex'

The next agent start commented that line out, with "disabled by DELFIN:
KIT_TOOLBOX_API_KEY now lives in the store" above it, and the user's
alias stopped working. The line holds no key and exports none; it names
one. ``exported_in_shell_files`` took any line that contained the name
and an ``=`` anywhere for an export -- and the automatic remediation acts
on what that function returns.

A line counts when the name is SET there: the left side of an
assignment, an ``export``/``declare -x``/``typeset -x`` of it, a csh
``setenv``, a pam_env entry. And the readable-file check reports a line
only when a value is actually written in it -- a ``$REFERENCE`` or a
``$(command)`` puts no key in the file.
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import credentials as cred

NAME = "KIT_TOOLBOX_API_KEY"

READS_ONLY = [
    "alias codex-kit='OPENAI_API_KEY=$KIT_TOOLBOX_API_KEY "
    "OPENAI_BASE_URL=\"https://example.invalid/api/v1\" codex'",
    'echo "$KIT_TOOLBOX_API_KEY"',
    '[ -n "$KIT_TOOLBOX_API_KEY" ] && have_key=yes',
    'curl -H "Authorization: Bearer ${KIT_TOOLBOX_API_KEY}" url=https://x',
    'codex-kit() { OPENAI_API_KEY="$(. "$HOME/.kit_env" && '
    'printf %s "$KIT_TOOLBOX_API_KEY")" codex "$@"; }',
]

SETS = [
    "export KIT_TOOLBOX_API_KEY=abc123",
    'KIT_TOOLBOX_API_KEY="abc123"',
    "  export FOO=1 KIT_TOOLBOX_API_KEY=abc123",
    "declare -x KIT_TOOLBOX_API_KEY=abc123",
    "export KIT_TOOLBOX_API_KEY",
    "setenv KIT_TOOLBOX_API_KEY abc123",
    "KIT_TOOLBOX_API_KEY DEFAULT=abc123",
    'export KIT_TOOLBOX_API_KEY="$(cat ~/.kit_key)"',
]


def _home(tmp_path, line: str):
    rc = tmp_path / ".bashrc"
    rc.write_text("# shell\n" + line + "\n", encoding="utf-8")
    os.chmod(rc, 0o644)
    return tmp_path, rc


@pytest.mark.parametrize("line", READS_ONLY)
def test_a_line_that_reads_the_key_is_not_found(tmp_path, line):
    home, _rc = _home(tmp_path, line)
    assert cred.exported_in_shell_files(NAME, home=home) == []


@pytest.mark.parametrize("line", READS_ONLY)
def test_a_line_that_reads_the_key_is_not_edited(tmp_path, line):
    home, rc = _home(tmp_path, line)
    before = rc.read_text(encoding="utf-8")
    assert cred.comment_out_exports(NAME, home=home) == []
    assert rc.read_text(encoding="utf-8") == before


@pytest.mark.parametrize("line", SETS)
def test_a_line_that_sets_the_key_is_still_found(tmp_path, line):
    home, _rc = _home(tmp_path, line)
    assert [n for _p, n, _l in cred.exported_in_shell_files(NAME, home=home)] == [2]


def test_the_remediation_never_edits_an_alias(tmp_path):
    """An assignment inside an alias lives for one command; it exports
    nothing into the shell, so commenting it out stops no export -- it
    only breaks the alias. Same incident, from the other key's side."""
    home, rc = _home(tmp_path, READS_ONLY[0])
    before = rc.read_text(encoding="utf-8")
    assert cred.comment_out_exports("OPENAI_API_KEY", home=home) == []
    assert rc.read_text(encoding="utf-8") == before


def test_the_incident_alias_is_not_a_key_in_a_readable_file(tmp_path):
    """The value it sets OPENAI_API_KEY to is a reference: no key in the file."""
    home, _rc = _home(tmp_path, READS_ONLY[0])
    assert cred.keys_in_readable_shell_files(("OPENAI_API_KEY", NAME),
                                             home=home) == []


def test_a_literal_key_in_an_alias_is_still_reported(tmp_path):
    home, _rc = _home(
        tmp_path,
        "alias codex-kit='OPENAI_API_KEY=\"sk-live-abcdef123456\" codex'")
    rows = cred.keys_in_readable_shell_files(("OPENAI_API_KEY",), home=home)
    assert [(r["name"], r["line_no"]) for r in rows] == [("OPENAI_API_KEY", 2)]


def test_an_export_from_a_command_is_not_a_key_in_the_file(tmp_path):
    """It still exports the key -- the remediation still finds it above --
    but nothing in the file is the key, so the readable-file check is quiet."""
    home, _rc = _home(tmp_path, 'export KIT_TOOLBOX_API_KEY="$(cat ~/.kit_key)"')
    assert cred.keys_in_readable_shell_files((NAME,), home=home) == []
