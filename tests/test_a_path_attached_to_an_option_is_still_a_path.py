"""A path written onto its option (``-f/etc/shadow``) is seen by the path
scanners, exactly as the same path written as its own word.

An option value may be attached to a short option with no space, which is
how every POSIX short option takes a value. The scanners skipped every
token starting with ``-`` unless it held an ``=``, so ``grep -f /etc/x``
was gated and ``grep -f/etc/x`` -- the same read -- was not. Measured on
2026-11-13 against the real functions.

The coarse scanner additionally read ``-f./pats.txt`` as the absolute
path ``/pats.txt``: its regex starts at the first ``/`` in the token, and
a workspace-relative value became a path outside the workspace that does
not exist. Both are pinned here.
"""

from pathlib import Path

import pytest

import delfin.agent.api_client as api_client

WS = Path("/srv/work/project")


@pytest.mark.parametrize("command, value", [
    ("grep -f/etc/nope -r x .", "/etc/nope"),
    ("grep -f~/notes -r x .", "~/notes"),
    ("sort -T/srv/other/tmp a.txt", "/srv/other/tmp"),
    ("rg -f/etc/nope x .", "/etc/nope"),
])
def test_an_attached_option_value_is_read_like_its_own_word(command, value):
    """The attached spelling reports the same path as the spaced one."""
    spaced = command.replace("-f" + value, "-f " + value).replace(
        "-T" + value, "-T " + value)
    attached = api_client._bash_outside_reads(command, WS)
    expanded = str(Path(value).expanduser())
    assert any(expanded in entry or value in entry for entry in attached), \
        f"{command!r} -> {attached!r}"
    # The spaced spelling is the behaviour this follows; it must not change.
    assert any(expanded in entry or value in entry
               for entry in api_client._bash_outside_reads(spaced, WS))


def test_the_coarse_scanner_sees_an_attached_path_outside():
    """Under a locked scope the question is whether a path outside the
    workspace is named at all. ``/etc`` exists on every machine the tests
    run on; the file itself need not."""
    found = api_client._bash_paths_outside("grep -f/etc/nope -r x .", WS)
    assert "/etc/nope" in found, found


def test_a_relative_attached_value_is_not_an_absolute_path():
    """``-f./pats.txt`` names a file in the workspace, not ``/pats.txt``."""
    found = api_client._bash_paths_outside("grep -f./pats.txt -r x .", WS)
    assert "/pats.txt" not in found, found


@pytest.mark.parametrize("command", [
    "grep -rn foo delfin",
    "grep -A3 -B2 foo delfin",
    "ls -la",
    "sort -t: -k2 -nr a.txt",
    "rg -n -g '*.py' foo delfin",
])
def test_an_option_without_a_path_adds_nothing(command):
    """Flags that carry no path stay invisible to both scanners: the fix
    must not turn every short option into a path."""
    assert api_client._bash_paths_outside(command, WS) == []
    outside = api_client._bash_outside_reads(command, WS)
    assert not any(entry.startswith("/etc") or entry.startswith("/srv/other")
                   for entry in outside), outside
