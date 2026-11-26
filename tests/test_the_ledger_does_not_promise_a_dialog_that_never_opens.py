"""The ledger must keep saying the flag class runs unasked, because it does.

Input: ``tests/fixtures/unmounted_capabilities.txt`` read as text, and the
product's own caged-admission decision measured with the isolation decision
pinned (so the answer does not depend on what this host's kernel offers).

Output measured 2026-11-25 on ``e72a4331c`` through ``_run_permission_gate``:

    sort --compress-program=sh a.txt   cage held   -> RUN, 0 dialogs, admitted
                                       no layer    -> ASK attended, REFUSE headless
    ls -la, cat a.txt, grep -r         cage held   -> RUN

The base fixture already says this class "RUN, no question", and two real
product runs confirmed it: ``sort --compress-program=./p.sh a.txt --`` ran on
both base and branch, attended, and the script it named was executed. So the
base prose was RIGHT about the flag class, and a package must not "close" it
on the strength of the no-cage half.

What this test exists for: an earlier revision of this very branch corrected
the flag paragraph to "The flag class ASKS attended" -- true only with no
layer confining the command, and false in the state a real session is in. The
assertions below are what that error trips, so the next edit cannot repeat it.
The caged rule admits a command whose text names nothing on
``_CAGED_UNASKED_EXCLUDED_WORDS``; ``sh`` is not on that list.

Why this is not "the gate is broken": where the cage holds the command the run
is confined (writes outside the workspace, reads outside the roots and
protected files are all pinned by
``test_a_caged_command_runs_without_asking.py``). What must not be claimed is
a *dialog* that never opens.
"""
from __future__ import annotations

import pytest

import delfin.agent.api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor

BASELINE = "tests/fixtures/unmounted_capabilities.txt"
FLAG = "sort --compress-program=sh a.txt"


def _block(text: str, row: str) -> str:
    """The contiguous run of comment lines that opens with ``# <row>``."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(f"# {row}"):
            block = [line]
            for nxt in lines[i + 1:]:
                if not nxt.startswith("#"):
                    break
                block.append(nxt)
            return "\n".join(block)
    return ""


@pytest.fixture
def fixture_text() -> str:
    with open(BASELINE, encoding="utf-8") as fh:
        return fh.read()


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    (d / "a.txt").write_text("banana\napple\ncherry\n")
    return d


def _gate(cmd: str, *, caged: bool, attended: bool, ws):
    """(verdict, dialogs, admissions) of the native bash gate, isolation pinned."""
    orig = (A._bash_isolation_decision, A._process_cage_enabled,
            A._landlock_functional)
    kind = "bwrap" if caged else "none"
    A._bash_isolation_decision = lambda perms, mode=None: (kind, "auto")
    A._process_cage_enabled = lambda: caged
    A._landlock_functional = lambda: caged
    asked: list[str] = []
    try:
        def _cb(tool, args, preview):
            asked.append(preview)
            return False

        perms = KitToolPermissions(
            mode="default", workspace=str(ws),
            confirm_callback=_cb if attended else None)
        admitted: list[bool] = []
        verdict = _doc_executor._run_permission_gate(
            "bash", {"command": cmd}, perms, caged_admission=admitted)
        return verdict, asked, admitted
    finally:
        (A._bash_isolation_decision, A._process_cage_enabled,
         A._landlock_functional) = orig


def test_the_flag_class_runs_unasked_where_a_cage_holds_it(ws):
    """The measured fact the ledger states: no dialog, admitted.

    Failing input: a build where the caged rule stops admitting a command
    naming ``sh`` -- that would make the class safe by asking, and this test
    (with the prose below) must then be revisited, not silently kept.
    """
    verdict, asked, admitted = _gate(FLAG, caged=True, attended=True, ws=ws)
    assert verdict is None, (
        f"expected the caged rule to admit it, got a refusal: {verdict!r}")
    assert asked == [], (
        f"a dialog was raised ({asked!r}); if this is now true, the class is "
        "guarded and the ledger paragraph should be rewritten to say so")
    assert admitted == [True], f"not admitted by the caged rule: {admitted!r}"
    # The reading controls the package promised keep running, so the
    # assertion above cannot be satisfied by admitting nothing at all.
    for cmd in ("ls -la", "cat a.txt", "grep -r banana ."):
        v, _a, _ad = _gate(cmd, caged=True, attended=False, ws=ws)
        assert v is None, f"control {cmd!r} was refused: {v!r}"


def test_without_a_layer_it_asks_and_then_refuses(ws):
    """The no-cage half, pinned so a correction does not overstate either way.

    This is the half that is true, and that a careless edit could quote as if
    it described every session.
    """
    verdict, asked, admitted = _gate(FLAG, caged=False, attended=True, ws=ws)
    assert verdict is not None and asked, (
        "attended with no layer must ask: "
        f"verdict={verdict!r} dialogs={asked!r}")
    verdict, asked, admitted = _gate(FLAG, caged=False, attended=False, ws=ws)
    assert verdict is not None, "headless with no layer must be refused"
    assert admitted == [], f"admitted without a cage: {admitted!r}"


def test_the_fixture_does_not_claim_a_dialog_that_never_opens(fixture_text):
    """The ledger must not describe the flag class as asked about.

    Failing input: an earlier revision of this branch wrote "The flag class
    ASKS attended" here, which one no-cage measurement suggested and the real
    runs refuted.
    """
    block = _block(fixture_text, "is_safe")
    assert block, f"no is_safe block in {BASELINE}"

    # Judge the flag paragraph itself, not the whole block: "asked about"
    # appears in unrelated prose (a secret path is asked about), which is
    # true and has nothing to do with the flag class.
    flagged = [
        ln for ln in block.splitlines()
        if "compress-program" in ln or "flag!" in ln
    ]
    assert flagged, "the block no longer mentions the flag class at all"
    para = "\n".join(flagged)
    for false_claim in ("ASKS attended", "it asks", "a dialog attended"):
        assert false_claim not in para, (
            f"the flag paragraph claims {false_claim!r}, which the "
            "measurement above contradicts where a cage holds the command "
            "(admitted, no dialog)")
    assert "RUN, no question" in block or "runs unasked" in block, (
        "the block must still record that this class runs unasked -- that is "
        "the open defect the row exists for")
    flat = " ".join(block.replace("# ", " ").split())
    assert "sh` is not on that list" in flat or "sh is not on that list" in flat, (
        "the block must name WHY, so a reader can check it: the caged rule's "
        "exclusion list does not carry `sh`")
