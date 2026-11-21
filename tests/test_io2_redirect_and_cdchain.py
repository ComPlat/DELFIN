"""IO2: redirect forms and chained ``cd`` must be judged like the shell.

Three defects found by reviewer nacht-s2 (round 10, NOT APPROVED):

1. A **here-string** is the batch script but was never inspected:
   ``sbatch <<< "du -sh /pfs/..."`` was allowed. Real bash delivers the
   here-string word as stdin (``bash <<< 'echo ran'`` executes it), so
   for ``sbatch`` it is the script text -- exactly like a here-document
   body, and it must be judged as text.

2. A **here-document whose delimiter is not a bare word** was never
   read: ``sbatch <<$TAG`` / body / ``TAG``. Bash does NOT expand the
   delimiter: the terminator is the literal text ``$TAG``. The gate's
   tag pattern did not match it, so no body was found and the submit
   was allowed. Delimiting what cannot be delimited must fail closed --
   the remainder is the body.

3. A **chained relative ``cd``** over-refused a clean submit:
   ``cd A && cd relB && sbatch clean.sh`` was refused with "not found
   from the cd target 'A'", because the second ``cd`` (behind a
   separator, with a leading space) was not recognised and the pending
   dir stayed ``A``. The safe direction, but it blocks a legitimate
   submit whose script is right there.

Every "must refuse" case below is paired with a "must still run" case:
fixing an over-refusal must never open a fail-open, and closing a
fail-open must never block a clean command.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "du -sh /pfs/data6/workspace/\n"
CLEAN = "sbatch --version\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


# --------------------------------------------------------------------------
# 1. here-strings
# --------------------------------------------------------------------------

def test_a_here_string_script_is_refused(lustre):
    """`sbatch <<< '<walking text>'` must not reach sbatch.

    INTENT UNCHANGED (the walking script must not be submitted); the
    DECISION changed in round 14. Rounds 11-13 judged the here-string word
    as script text and refused only when that text walked. The QS then
    ruled that sbatch's stdin channel is refused outright, because a
    script that arrives that way cannot be identified afterwards -- so the
    verdict is still `refuse`, now with the stdin-channel reason.
    """
    verdict, text, may_run = _v(f'sbatch <<< "{WALK.strip()}"')
    assert verdict == "refuse", (
        f"FAIL-OPEN: a here-string is the batch script and was not refused "
        f"(verdict={verdict!r})")
    assert "job.sh" in text, (
        f"the refusal must name the file route, not only say no: {text[:160]}")
    assert may_run is False


def test_a_clean_here_string_is_also_refused(lustre):
    """Ruling (a) refuses the CHANNEL, so a clean body refuses too.

    This REVERSES the rounds 11-13 decision (a clean here-string ran,
    because its text was judged and was clean). QS ruling, round 14: the
    point is that no artefact is left to inspect, not that this particular
    text is harmless. Pinned so nobody relaxes it back.
    """
    verdict, text, may_run = _v('sbatch <<< "sbatch --version"')
    assert verdict == "refuse", (
        f"a clean here-string must still refuse under ruling (a); got "
        f"{verdict!r}")
    assert "job.sh" in text
    assert may_run is False


# --------------------------------------------------------------------------
# 2. here-documents whose delimiter is not a bare word
# --------------------------------------------------------------------------

def test_a_variable_tagged_heredoc_body_is_refused(lustre):
    """`<<$TAG` is a literal delimiter in bash; its body is the script.

    INTENT UNCHANGED (this body must never reach sbatch); the reason is
    the stdin channel itself under ruling (a), round 14.
    """
    cmd = "sbatch <<$TAG\n" + WALK + "$TAG\n"
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN: a heredoc tagged with a variable was not refused "
        f"(verdict={verdict!r})")
    assert "job.sh" in text


def test_an_unterminated_variable_tagged_body_is_refused(lustre):
    """A body whose terminator never matches runs to the end in bash."""
    cmd = "sbatch <<$TAG\n" + WALK + "TAG\n"
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN: an unterminated variable-tagged body was not refused "
        f"(verdict={verdict!r})")
    assert "job.sh" in text


def test_a_literal_tagged_heredoc_still_refuses(lustre):
    """Control: the plain form must keep working (no regression)."""
    verdict, text, may_run = _v("sbatch <<TAG\n" + WALK + "TAG\n")
    assert verdict == "refuse"
    assert "job.sh" in text


def test_a_clean_variable_tagged_body_is_refused_too(lustre):
    """Ruling (a): a CLEAN heredoc body is refused as well.

    REVERSES the rounds 11-13 decision (clean body ran). The channel, not
    the text, is what ruling (a) refuses.
    """
    verdict, text, may_run = _v("sbatch <<$TAG\n" + CLEAN + "$TAG\n")
    assert verdict == "refuse", (
        f"a clean heredoc body must refuse under ruling (a); got "
        f"{verdict!r}")
    assert "job.sh" in text


# --------------------------------------------------------------------------
# 3. chained cd must not over-refuse
# --------------------------------------------------------------------------

def _chain(tmp_path: Path) -> Path:
    """A/relB/ containing a walking and a clean submit script."""
    sub = tmp_path / "A" / "relB"
    sub.mkdir(parents=True)
    (sub / "walk.sh").write_text("#!/bin/bash\n" + WALK, encoding="utf-8")
    (sub / "clean.sh").write_text("#!/bin/bash\n" + CLEAN, encoding="utf-8")
    return tmp_path / "A"


def test_a_clean_submit_behind_a_chained_relative_cd_runs(lustre, tmp_path):
    """`cd A && cd relB && sbatch clean.sh` -- the script IS right there."""
    a = _chain(tmp_path)
    verdict, text, may_run = _v(f"cd {a} && cd relB && sbatch clean.sh")
    assert verdict != "refuse", (
        f"OVER-REFUSAL: a clean submit behind a chained relative cd was "
        f"refused: {text[:200]}")


def test_a_walking_submit_behind_a_chained_relative_cd_still_refuses(
        lustre, tmp_path):
    """Control: the same chain with a walking script MUST still refuse.

    This is the case that stops the over-refusal fix from becoming a
    fail-open: the second cd must reach the script, not bypass it.
    """
    a = _chain(tmp_path)
    verdict, text, may_run = _v(f"cd {a} && cd relB && sbatch walk.sh")
    assert verdict == "refuse", (
        f"FAIL-OPEN: a walking submit behind a chained relative cd was "
        f"allowed (verdict={verdict!r})")
    assert "tree_walk_on_lustre" in text
    assert "not found from the cd target" not in text
