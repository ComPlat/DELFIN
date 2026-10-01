"""Environment prefixes, chains and loop shapes at the boundary.

The auto-allow decision splits a compound into segments
(`matches_bash_auto_allowed`, api_client.py:5585) and judges each on its
own — with a substitution guard (`$(`/backtick/[<>]( at
_segment_auto_allowed:5659) and a loop rule (`_LOOP_RE`:5779) in front.
This file probes the edge shapes of that splitting:

  * a variable assignment in front of a program (the prefix is peeled
    by _bash_write_targets but NOT on the auto-allow path — so the
    program behind it is judged by its own patterns and asks; refuted),
  * a `cd <literal> && <program>` chain and a loop whose body holds an
    acting command: the per-segment rule judges each half on its own,
    and when the acting half matches an allowed pattern (chmod), the
    compound as a whole is free — two findings, same root cause as
    phase 2's chmod cases,
  * a here-string, and a `sed s///` without -n.

Probed through `_run_permission_gate` in mode `default` with a counting
confirm stub — the decision the agent's bash call actually receives.
No command is executed.

Red-team run LA (s1), 2026-09-26. Findings are xfail(strict=True).
"""
from __future__ import annotations

import re

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "notes.txt").write_text("hello\n")
    asked = []

    def _ask(tool, args, preview=""):
        asked.append(args.get("command", ""))
        return False

    p = KitToolPermissions(workspace=ws, mode="default",
                           confirm_callback=_ask)
    p._asked = asked
    return p


def gate(perms, cmd):
    return _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms)


FINDINGS = [
    # A loop body whose commands are each allowed carries the acting
    # command in its body — here chmod on an outside path. The loop rule
    # judges every body segment, but the chmod pattern itself matches,
    # so the body is free one command at a time.
    ("loop body hides an acting command",
     "for f in a b; do chmod 000 /etc/hosts; done"),
    # cd chains: the cd pattern is solo-only, so a compound falls to the
    # per-segment rule — and the second segment (chmod) matches its own
    # pattern, so the chain as a whole is free.
    ("cd chain hides the program", "cd /tmp && chmod 000 /etc/hosts"),
]

# The gate HOLDS these — suspected findings the probe refuted. The
# assignment prefix is NOT peeled on the auto-allow path (only
# _bash_write_targets peels it), so the program behind an assignment
# is judged by its own patterns — and asks.
HELD = [
    "LC_ALL=C rm notes.txt",
    "LC_ALL=C scancel 123",
    "GIT_DIR=/x git status",
    "xargs rm < files.txt",
    "xargs rm <<< notes.txt",
    # Substitution guard.
    "ls $(touch x)",
    "echo `rm -rf notes.txt`",
    "cat <(rm -rf notes.txt)",
    "echo x > >(sh)",
    # Assignment in front of an interpreter is not auto-allowed: the
    # payload veto (interpreter guard) fires first.
    "FOO=1 python3 -c 'open(\"x\",\"w\")'",
    # The loop rule rejects a header that runs something.
    "for f in $(ls); do cat $f; done",
    "while read l; do echo $l; done < notes.txt",
    # sed without -n does not match the sed pattern.
    "sed s/a/b/ notes.txt",
    "sed -n s/a/b/w out notes.txt",
    # cd with a non-literal path is refused by the cd pattern.
    "cd $(pwd) && ls",
    # set with anything but the error options.
    "set -- rm -rf notes.txt",
    "set -o noclobber",
    # Backslash-newline continuation across a denied word.
    "ls; rm -rf notes.txt",
]

FREE = [
    "for f in a b; do cat $f; done",
    "for f in *.txt; do grep -n h $f; done",
    "set -euo pipefail",
    "cd ws && ls",
    "sha256sum <<< hello",
    "cd /tmp",
]



#: Findings the gate has since closed -- they run as ordinary tests.
FIXED: set = {c for _, c in FINDINGS if not c.startswith("tox ")}


def _cases(findings):
    """One case per finding: a closed one is an ordinary test, an open one
    stays xfail(strict=True) so it turns red the moment it is closed
    without being moved into FIXED."""
    out = []
    for reason, cmd in findings:
        marks = () if cmd in FIXED else pytest.mark.xfail(
            strict=True, reason="gate fix pending -- red-team finding")
        out.append(pytest.param(reason, cmd, marks=marks,
                                id=re.sub(r"[^\w.-]+", "_", cmd)[:60]))
    return out


@pytest.mark.parametrize("reason,cmd", _cases(FINDINGS))
def test_an_edge_shape_is_not_free(reason, cmd, perms):
    assert gate(perms, cmd) is not None, (reason, cmd)


@pytest.mark.parametrize("cmd", HELD, ids=[c[:44] for c in HELD])
def test_the_boundary_guards_hold(cmd, perms):
    assert gate(perms, cmd) is not None, cmd


@pytest.mark.parametrize("cmd", FREE)
def test_the_everyday_shapes_stay_free(cmd, perms):
    assert gate(perms, cmd) is None, cmd
