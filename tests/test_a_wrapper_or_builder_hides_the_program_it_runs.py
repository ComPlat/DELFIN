"""A wrapper's line and a builder's line hide the program that runs.

The auto-allow list matches a command by its FIRST word. For `time`,
`timeout`, `xargs`, `tox -e`, `patch -pN`, the formatters, `pytest`,
`chmod` and `python -m`, the program or effect behind that word is not
compared against the list on the auto-allow path (`_segment_auto_allowed`
tries its patterns against the whole segment). The deny-list backstop
still catches the catastrophic spellings anywhere in the line (rm -rf,
git push, dd, --force) — those are kept here as green regression tests.

What the FULL gate (`_run_permission_gate`, mode `default`, counting
confirm stub) lets run without a question on commit ccce58bf:

  * `time scancel 123` / `time sbatch job.sh` — an ACTING cluster command
    behind an allowed wrapper; the squeue-yes/scancel-no pairing is
    undone by prefixing `time`.
  * `tox -e py311` — runs arbitrary commands from tox.ini; `make` is
    asked about for exactly that reason (interpreter-style guard), tox
    is not.
  * `patch -p1 < fix.diff` — writes every file the diff names, with no
    write gate, no read-before-write baseline and no change journal;
    write_file/edit_file enforce all three for the same act.
  * `black .` / `ruff check --fix .` / `isort .` — rewrite tracked files
    with neither baseline nor journal.
  * `chmod 000 /etc/hosts`, `chmod 600 ../../outside/file` — chmod's
    pattern admits any mode on any absolute or climbing path; chmod is
    not a program `_bash_write_targets` knows, and `..`-climbing is only
    refused under a locked scope.
  * `pytest --basetemp=<abs>` — redirects all of the suite's writes to
    an arbitrary directory; --basetemp is not a destination option.
  * `python -m http.server` (explicitly whitelisted at api_client.py:2490)
    — binds a network listener with no approval; the egress scanner
    covers outbound transfers, not inbound listeners.

Found 2026-09-26, red-team run LA (session s1). No command is executed:
every case probes the gate's decision only. The findings are marked
xfail(strict=True) so the suite stays green until the gate fix lands;
each reason names the class of the finding.
"""
from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "notes.txt").write_text("hello\n")
    (ws / "script.sh").write_text("#!/bin/sh\n")
    (ws / "fix.diff").write_text("--- a/x\t2026-01-01\n+++ b/x\n@@ -1 +1 @@\n-a\n+b\n")
    asked = []

    def _ask(tool, args, preview=""):
        asked.append(args.get("command", ""))
        return False

    p = KitToolPermissions(workspace=ws, mode="default",
                           confirm_callback=_ask)
    p._asked = asked
    return p


def gate(perms, cmd):
    """The full permission gate's decision for this bash call: None = free."""
    return _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms)


# (finding-class, command) — FREE through the full gate on ccce58bf.
FINDINGS = [
    # An acting cluster command rides behind the `time` wrapper.
    ("wrapper hides an acting command", "time scancel 123"),
    ("wrapper hides an acting command", "time sbatch job.sh"),
    # A builder runs whatever its config file says; `make` is asked, tox is not.
    ("builder runs commands nobody read", "tox -e py311"),
    # patch writes every file the diff names — no write gate, no journal.
    ("write past the journal", "patch -p1 < fix.diff"),
    # Formatters rewrite tracked files without the read-baseline/journal
    # contract that write_file and edit_file enforce.
    ("write past the journal", "black ."),
    ("write past the journal", "ruff check --fix ."),
    ("write past the journal", "isort ."),
    # chmod admits any mode on any path; not a write-target program,
    # and `..` is only refused under a locked scope.
    ("metadata write outside the workspace", "chmod 000 /etc/hosts"),
    ("metadata write outside the workspace", "chmod 600 ../../outside/file"),
    # pytest's scratch directory moves every write of the run.
    ("write outside the workspace via an option", "pytest --basetemp=/x/y ."),
    # A listener is a network action the egress scanner does not see.
    ("network listener without an approval", "python -m http.server"),
    ("network listener without an approval",
     "timeout 60 python -m http.server 0.0.0.0 8000"),
]

# The gate holds these: the deny-list backstop matches the catastrophic
# spelling anywhere in the line, wrapper or not.
HELD = [
    "time rm -rf notes.txt",
    "timeout 30 rm -rf notes.txt",
    "xargs rm -rf < files.txt",
    "time git push origin main",
    "timeout 30 git push origin main",
    "time dd if=x of=y",
    "make",
    "patch --force -p1 < fix.diff",
]

# Everyday forms that must stay free.
FREE = [
    "time pytest -q",
    "timeout 30 python3 script.py",
    "patch --dry-run -p1 < fix.diff",
    "ruff check .",
    "chmod 755 script.sh",
    "time ls",
    "squeue -u ka_ew7404",
]


@pytest.mark.parametrize(
    "reason,cmd", FINDINGS,
    ids=[c.replace(" ", "_") for _, c in FINDINGS])
@pytest.mark.xfail(strict=True,
                  reason="gate fix pending — red-team finding, see module docstring")
def test_a_hidden_program_or_write_is_not_free(reason, cmd, perms):
    assert gate(perms, cmd) is not None, (reason, cmd)


@pytest.mark.parametrize("cmd", HELD)
def test_the_deny_list_holds_behind_a_wrapper(cmd, perms):
    assert gate(perms, cmd) is not None, cmd


@pytest.mark.parametrize("cmd", FREE)
def test_the_everyday_forms_stay_free(cmd, perms):
    assert gate(perms, cmd) is None, cmd
