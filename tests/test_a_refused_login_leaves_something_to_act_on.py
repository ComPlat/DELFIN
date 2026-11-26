"""A refused first login leaves a file to send, not a puzzle.

The launcher registers a PC during its first PASSWORD login: OpenSSH asks,
the person types, and remote_bootstrap.sh appends the key. Where the
server takes no passwords the circle does not close -- the key is added
once you are in, and you get in with the key -- and a new person simply
stops there.

Breaking it means getting one line off their machine, and that is where a
lay user is lost: the key sits in a dotted folder, and the route to it was
a console command printing a file. So the launcher now writes the PUBLIC
half to their DESKTOP as a file they can attach to a mail, and says so in
plain words. The private key never leaves the PC.

The second half of the same problem is the first login itself, when
passwords DO work: OpenSSH asks up to three times -- one prompt per
connection -- which reads as a refusal, and Ctrl+V in that console sends a
control character so a CORRECT password is refused. Both were written down
in a paragraph among six other lines, after the point where people stop
reading. They are now said before the first prompt.

Text-level, like every test for this package: DELFIN.ps1 is PowerShell and
there is no Windows here to run it on.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "tools" / "windows-launcher"
APP = (ROOT / "DELFIN.ps1").read_text(encoding="utf-8")


def _before_the_first_prompt() -> str:
    """What is on screen before the dashboard connection asks anything.

    Anchored on the two remote-command builders rather than on any text
    this change added, so the test cannot pass by quoting itself: the
    dashboard block runs after the terminal one, and its first SSH call
    is built from the script version taken just before it.
    """
    i_terminal = APP.index("Get-RemoteCommand $p 'terminal'")
    i_dashboard = APP.index("Get-RemoteCommand $p 'dashboard'")
    i_first_ssh = APP.rindex("$version = Get-ScriptVersion", 0, i_dashboard)
    return APP[i_terminal:i_first_ssh]


def _dashboard_branch() -> str:
    """From the dashboard connection to the end of the file."""
    return APP[APP.index("Get-RemoteCommand $p 'dashboard'"):]


def test_the_two_traps_are_named_before_the_first_prompt():
    before = _before_the_first_prompt()
    assert "up to three times" in before, (
        "being asked once per connection reads as a refusal")
    assert "RIGHT-CLICK" in before and "Ctrl+V does not work" in before, (
        "Ctrl+V in this console sends a control character, so a correct "
        "password is refused; that has to be said before it happens")
    assert "TYPE the password" in before


def test_a_refused_login_writes_the_key_to_the_desktop():
    branch = _dashboard_branch()

    assert "GetFolderPath('Desktop')" in branch
    assert "DELFIN-key-for-the-server-owner.txt" in branch
    assert "$loginKey.Public" in branch, "the public half, by name"


def test_only_the_public_half_is_written():
    """The file is meant to be mailed. Writing the private key into it
    would send the credential itself to whoever receives the mail."""
    branch = _dashboard_branch()
    start = branch.index("DELFIN-key-for-the-server-owner.txt")
    block = branch[start:start + 1800]

    assert "$loginKey.Key" not in block, (
        "that is the PRIVATE key file; only .Public may be written")
    assert "not your password" in block


def test_it_only_fires_when_the_dashboard_never_opened():
    """A disconnect at the end of a good session is not a refused login,
    and dropping this file on the desktop every time would make it noise
    that nobody reads when it matters."""
    branch = _dashboard_branch()
    guard = next(l for l in branch.splitlines()
                 if "GetFolderPath('Desktop')" in l
                 or ("ExitCode -ne 0" in l and "$opened" in l))

    assert "-not $opened" in guard, guard
    assert "$sshProcess.ExitCode -ne 0" in guard, guard


def test_it_says_what_the_other_side_has_to_do():
    """A key with no instruction is a key that sits in someone's inbox."""
    branch = _dashboard_branch()

    assert "delfin-people --add-key" in branch, (
        "name the one command, so the server's owner need not know it")


def test_a_desktop_that_cannot_be_written_still_shows_the_key():
    """A redirected or locked-down desktop must not swallow the one thing
    the person needs."""
    branch = _dashboard_branch()
    start = branch.index("DELFIN-key-for-the-server-owner.txt")
    block = branch[start:start + 2200]

    assert "} catch {" in block
    tail = block[block.index("} catch {"):]
    assert "$loginKey.Public" in tail, "print it when it cannot be saved"
