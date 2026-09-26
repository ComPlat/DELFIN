"""A refused bash command is recognised however it is quoted (LJ phase 4, landed).

Landed: the operator implemented tier 1 in api_client
(_denied_action_key), and these tests pin it.

The gap (pinned as current behaviour in
test_how_far_a_refusal_reaches.py::test_bash_refusal_locks_only_the_exact_command):
a refused command is keyed on its whitespace-collapsed text only
(api_client.py:5441-5442), so the same command with quote characters
around its arguments is a different key and raises a second dialog.

The proposal, deliberately narrow:

* Tier 1 (this file): strip quote characters from the command before
  the key is built, the same normalisation _bash_denied_path already
  applies when SCANNING a command (api_client.py:13609 -- "Keep quoted
  CONTENT, drop only the quote characters"). A refusal should recognise
  what it refused by the same reading the scanner uses. This closes
  `cp 'engine.py' '/tmp/x'` against a refused `cp engine.py /tmp/x`
  without locking anything path-wide.

* Tier 2 (rejected, documented): keying on program-sequence + write
  targets, so a differently-ARGUMENTED variant of the same act stays
  refused. Rejected on the wave-5 evidence itself: the operator refused
  `cp engine.py /tmp/x && git checkout abc -- engine.py` BECAUSE the
  backup went to the private /tmp of the cage. Re-asking the same act
  with a different sink is not a workaround -- it is the corrected ask
  the user's reason was asking for. A gate that swallows it re-creates
  a refusal broader than the decision it records.

The write-target side of the assignment ("the write targets of a
refused command apply only to that same command") is already the
implemented state: a denied bash command locks no paths at all -- pinned
in test_intent_a_denied_bash_command_does_not_lock_its_targets.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A


def _perms(tmp_path, **kw):
    return A.KitToolPermissions(workspace=tmp_path, **kw)


def _executor():
    return A._DocToolExecutor()


def _deny(name, args, preview):
    return False


def test_a_quoted_re_spelling_of_a_refused_command_stays_refused(tmp_path):
    ex = _executor()
    perms = _perms(tmp_path, mode="default", confirm_callback=_deny)
    cmd = "cp engine.py /tmp/x && git checkout abc123 -- engine.py"
    first = ex._run_permission_gate("bash", {"command": cmd}, perms)
    assert first is not None and "denied" in first

    quoted = "cp 'engine.py' '/tmp/x' && git checkout abc123 -- engine.py"
    again = ex._run_permission_gate("bash", {"command": quoted}, perms)
    assert again is not None and "already refused" in again, again


def test_a_double_quoted_re_spelling_stays_refused(tmp_path):
    ex = _executor()
    perms = _perms(tmp_path, mode="default", confirm_callback=_deny)
    # A command that asks (a plain read like `sed -n` runs free and is
    # never refused in the first place).
    ex._run_permission_gate(
        "bash", {"command": "rm notes.txt old.txt"}, perms)

    again = ex._run_permission_gate(
        "bash", {"command": 'rm "notes.txt" "old.txt"'}, perms)
    assert again is not None and "already refused" in again, again
