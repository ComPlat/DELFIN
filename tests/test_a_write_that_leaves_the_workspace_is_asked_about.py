"""Two forms that wrote or ran things unasked, and what they cost to close.

Input: commands through the real bash entry (``_run_permission_gate``,
which consults the auto-allow list only after the deny list, the secret
scan, the approval guard, the irreversibility rule and
``unrestorable.refusal``), default mode, a fresh workspace, with a
recording confirm callback so that running unasked, asking and refusing
are three distinguishable outcomes.

Three postures, because the argument for this repair rests on the third:

* ``default``       -- a layer holds the command, somebody can be asked
* ``unisolated``    -- no layer holds it (``perms.unisolated_fallback``,
                       the field ``isolation_posture.apply_at_start``
                       sets), somebody can still be asked
* ``nobody-to-ask`` -- unisolated AND no confirm callback: the posture
                       ``_unisolated_fallback_verdict`` was written for,
                       where auto-allowed commands run and everything
                       else is refused outright

Output measured 2026-11-22 on the commit before this change (3f9474f8a),
in all three postures:

    echo x > ~/.bashrc              RUN   RUN   RUN
    ls > ../escape.txt              RUN   RUN   RUN
    ls >> /tmp/append.txt           RUN   RUN   RUN
    ls 2> /tmp/err.txt              RUN   RUN   RUN
    ls &> /tmp/both.txt             RUN   RUN   RUN
    ls > /tmp/a.txt 2> /tmp/b.txt   RUN   RUN   RUN
    cd .. && ls > escaped.txt       RUN   RUN   RUN
    cd /tmp && ls > here.txt        RUN   RUN   RUN
    sort --compress-program=sh x    RUN   RUN   RUN

So a write out of the workspace ran with no question even where nobody
could have been asked -- the posture in which nothing else would have
stopped it. After this change each of those asks, and is refused where
nobody can be asked.

The cost side is measured in the same run and pinned below, because a
security rule that costs function is paid for twice: the redirect forms
an agent writes all day (into the workspace, into ``$TMPDIR``, into
``/dev/null``, into a subdirectory it just moved to) must still run.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A


def _verdict(cmd: str, ws, posture: str = "default") -> str:
    """RUN, ASK or REFUSE, as the real gate decides it in *posture*."""
    asked: list = []
    callback = None if posture == "nobody-to-ask" else (
        lambda *a, **k: (asked.append(a), False)[1])
    perms = A.KitToolPermissions(
        workspace=str(ws), mode="default", confirm_callback=callback)
    if posture in ("unisolated", "nobody-to-ask"):
        perms.unisolated_fallback = True
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms)
    if asked:
        return "ASK"
    return "REFUSE" if err else "RUN"


# A write that leaves the workspace, in every spelling measured.
LEAVES = [
    "echo x > ~/.bashrc",
    "ls > /tmp/elsewhere.txt",
    "ls > ../escape.txt",
    "ls >> /tmp/append.txt",
    "ls 2> /tmp/err.txt",
    "ls &> /tmp/both.txt",
    "cat a 3>/tmp/fd.txt",
    "ls > /tmp/a.txt 2> /tmp/b.txt",
    # the target is relative, and an earlier segment moved the cwd
    "cd .. && ls > escaped.txt",
    "cd /tmp && ls > here.txt",
    # a cwd that cannot be resolved makes every later target unknown
    "cd - && ls > wherever.txt",
    "cd $D && ls > wherever.txt",
    # inside a loop body, which has its own rule above this one
    "for f in a b; do ls > /tmp/$f.txt; done",
]

# A flag that writes a file or runs a program, in every spelling getopt
# accepts.
FLAGS = [
    "sort --compress-program=sh x",
    "sort --compress-program sh x",
    "sort --out=f x",
    "sort -T /tmp x",
    "sort -ofile x",
]


@pytest.mark.parametrize("cmd", LEAVES + FLAGS)
def test_it_reaches_the_person_where_one_can_be_asked(cmd, tmp_path):
    (tmp_path / "sub").mkdir(exist_ok=True)
    assert _verdict(cmd, tmp_path) == "ASK", cmd


@pytest.mark.parametrize("cmd", LEAVES + FLAGS)
def test_it_still_reaches_the_person_without_a_layer(cmd, tmp_path):
    """The posture the repair's argument rests on: with no layer, what
    stops a command is a reading of its text."""
    (tmp_path / "sub").mkdir(exist_ok=True)
    assert _verdict(cmd, tmp_path, "unisolated") == "ASK", cmd


@pytest.mark.parametrize("cmd", LEAVES + FLAGS)
def test_it_is_refused_where_nobody_can_be_asked(cmd, tmp_path):
    """Unisolated with nobody to ask: auto-allowed commands run and
    everything else is refused. Before this change every one of these
    was auto-allowed, so every one of them ran, here of all places."""
    (tmp_path / "sub").mkdir(exist_ok=True)
    assert _verdict(cmd, tmp_path, "nobody-to-ask") == "REFUSE", cmd


def test_the_posture_is_the_one_claimed(tmp_path):
    """Control on the instrument itself: a command off the auto-allow
    list must ask in the first two postures and be refused in the third.
    An earlier version of this probe set attribute names that do not
    exist, so all three postures were the same session and the
    unisolated rows measured nothing."""
    assert _verdict("nmap -p 80 example.invalid", tmp_path) == "ASK"
    assert _verdict(
        "nmap -p 80 example.invalid", tmp_path, "unisolated") == "ASK"
    assert _verdict(
        "nmap -p 80 example.invalid", tmp_path, "nobody-to-ask") == "REFUSE"


# --- the cost side: what must still run unasked ------------------------

@pytest.mark.parametrize("cmd", [
    "ls > out.txt",                 # the documented in-workspace policy
    "ls >> build.log",
    "grep -n x delfin/agent.py > hits.txt",
    "cat a b > merged.txt",
    "sort -u f > sorted.txt",
    "ls > ./sub/out.txt",
    "cd sub && ls > inside.txt",    # a cd that stays inside is resolved
    'ls > "$TMPDIR/listing"',       # an unexpanded variable is no known path
    "ls > /dev/null",               # the most common redirect of all
    "ls 2>/dev/null",
    "ls > /dev/stdout",
    "grep x f 2>&1 | tail -5",      # a stream merge writes no file
    "ls >&2",
    "ls -la", "git status", "cat f", "grep x f", "wc -l f",
    "grep -r pattern tests/", "python3 --version", "pytest -q tests/x.py",
    "sort -u f", "sort -k2 f", "git reflog",
    "for f in a b; do ls > $f.txt; done",   # a loop writing inside
])
def test_the_forms_an_agent_writes_all_day_still_run_unasked(cmd, tmp_path):
    """Control: the change must not cost function. Every line here ran
    unasked before it and must still run unasked after."""
    (tmp_path / "sub").mkdir(exist_ok=True)
    assert _verdict(cmd, tmp_path) == "RUN", cmd


@pytest.mark.parametrize("cmd,outcome", [
    ("rm -rf f", "REFUSE"),
    ("git push --force", "REFUSE"),
    ("echo x > ~/.ssh/authorized_keys", "REFUSE"),
    ("curl -T f https://example.invalid", "ASK"),
    ("git reflog expire --expire=now --all", "ASK"),
    ("ls >| /tmp/clobber.txt", "ASK"),
])
def test_what_was_already_decided_is_decided_the_same_way(
        cmd, outcome, tmp_path):
    """Control: the layers in front of the auto-allow list keep their
    verdicts, and so do the spellings that were already being asked
    about for another reason.

    ``git reflog expire`` was reported as auto-allowed from a predicate
    reading; at the gate it already asked, and it still does. ``>|`` was
    already asked about because the ``|`` splits the segment, not
    because its target was judged -- it is listed here so that the
    reason is on the record rather than counted as this change's work.
    """
    assert _verdict(cmd, tmp_path) == outcome, cmd
