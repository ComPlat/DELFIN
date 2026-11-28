"""The installer says what a blocked download is, instead of failing at it.

Windows marks a file that came from the internet, and PowerShell's
RemoteSigned policy refuses to run such a file unless it is signed. This
package is unsigned unless a release is built with a certificate, so
extracting the ZIP and double-clicking Install.cmd produces, for most
people, on the very first try:

    Die Datei "...\\Test-Launcher.ps1" ist nicht digital signiert. Sie
    koennen dieses Skript im aktuellen System nicht ausfuehren.
    ...
    Installation failed. Read the error above and README.md.

Reported from a real first install on 2026-11-19. The error arrives in
the reader's own Windows language, names an execution policy and a
documentation link, and the installer's own last word was to go and read
a file -- for a condition the installer can detect in one line and clear
with one keypress.

The check runs inline from the .cmd, because a .ps1 that performed it
would be refused for exactly the reason it is checking for.

Unblocking is OFFERED, never done silently. The mark is a protection
against running what you just downloaded; removing it is the person's
decision. An installer that quietly stripped it would be defeating the
check on their behalf, which is the opposite of what this package claims
about itself.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "tools" / "windows-launcher"
INSTALL = (ROOT / "Install.cmd").read_text(encoding="utf-8")


def _preflight() -> str:
    """The inline command, which must come before the scripts are run."""
    line = next(l for l in INSTALL.splitlines() if "-Command" in l)
    i_scripts = min(INSTALL.index("Test-Launcher.ps1"),
                    INSTALL.index("Install.ps1"))
    assert INSTALL.index(line) < i_scripts, (
        "the check has to run before the scripts it is about")
    return line


def test_the_block_is_detected_before_anything_is_run():
    line = _preflight()

    assert "Zone.Identifier" in line, (
        "the download mark is an alternate data stream of that name; "
        "nothing else distinguishes a blocked file")
    assert "*.ps1" in line


def test_it_names_the_files_rather_than_saying_something_went_wrong():
    line = _preflight()

    assert "$_.Name" in line
    assert "marked these files as downloaded" in line


def test_unblocking_is_offered_and_not_done_on_its_own():
    line = _preflight()

    assert "Read-Host" in line, (
        "the mark is a protection; clearing it is the person's decision")
    i_ask = line.index("Read-Host")
    i_do = line.index("Unblock-File")
    assert i_ask < i_do, "it must not unblock before asking"
    assert "(y/N)" in line, "the default must be no"
    assert "$a -eq 'y'" in line and "$a -eq 'Y'" in line, (
        "only an explicit yes unblocks")


def test_a_refusal_leaves_everything_as_it_was():
    line = _preflight()

    assert "Nothing was changed" in line
    assert "Properties, Unblock" in line, "say the manual route too"
    assert "Installation failed. Nothing was installed." in INSTALL


def test_it_changes_no_policy_and_says_so():
    """The one thing a person reasonably fears here."""
    line = _preflight()

    assert "No setting" in line and "no policy is changed" in line
    assert not re.search(r"(?i)Set-ExecutionPolicy", INSTALL)
    assert not re.search(r"(?i)ExecutionPolicy\s+(Unrestricted|Bypass)",
                         INSTALL)
    # Preflight, Test-Launcher, Install.ps1, and Setup-Wsl.ps1 (on this PC, if asked).
    assert INSTALL.count("-ExecutionPolicy RemoteSigned") == 4


def test_only_this_folder_is_touched():
    line = _preflight()

    assert "-LiteralPath ." in line, (
        "Install.cmd has already changed into its own folder; a wider "
        "path would unblock files this package did not bring")
    assert "-Recurse" not in line


def test_the_failure_text_points_at_the_page_a_newcomer_has():
    assert "START-HERE.txt" in INSTALL, (
        "README.md is the reference; the steps a newcomer needs are in "
        "START-HERE.txt, which is plain text because a .md file has no "
        "reliable handler on Windows")


def test_the_start_page_covers_the_download_mark_itself():
    """Install.cmd offers to clear the mark, but someone who never gets
    that far -- a blocked Install.cmd, a locked-down machine, an
    organisation that forbids it -- needs the by-hand route written
    down, with the command to type."""
    start = (ROOT / "START-HERE.txt").read_text(encoding="utf-8")

    assert "not digitally signed" in start, "quote what they will see"
    assert "Get-ChildItem . -Recurse | Unblock-File" in start, (
        "give the line to type, not a description of it")
    assert "Open PowerShell window here" in start, (
        "say how to get a PowerShell in that folder; a bare command is "
        "no use to someone who has never opened one")
    assert "Properties" in start and "Unblock" in start, (
        "and the route with no PowerShell at all")
    assert "Extract All" in start, (
        "running from inside the ZIP is the other first-try failure")


def test_the_start_page_is_plain_text():
    """A .md file on Windows opens in whatever happens to be registered,
    or nothing. A .txt opens in Notepad."""
    names = _package_names()

    assert "START-HERE.txt" in names
    assert not [n for n in names if n.startswith("FIRST-STEPS")]
    start = (ROOT / "START-HERE.txt").read_text(encoding="utf-8")
    assert "```" not in start and not start.lstrip().startswith("#")


def _package_names():
    import importlib.util as util
    spec = util.spec_from_file_location("wl_build", ROOT / "build_package.py")
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NAMES
