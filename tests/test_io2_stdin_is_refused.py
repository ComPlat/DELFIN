"""IO2 round 14: an sbatch fed its script through stdin is REFUSED.

QS ruling on the operator's proposal (a), 2026-11-21: for `sbatch` the
stdin IS the job script, and a script that arrives through a pipe, a `<`
redirect, a here-document or a here-string cannot be identified as an
artefact afterwards. So the bash-tool route refuses it outright and tells
the person what to do instead: write the job script to a file and submit
that file, which IS checked.

This replaces the earlier decision (rounds 11-13) to reconstruct what the
producer wrote and judge that. That reconstruction is deleted with this
change: three rounds each closed one spelling and a fourth found another,
which is the signature of the wrong shape rather than a missing case.

Pinned here, in BOTH directions:
  * every stdin channel is refused -- and the refusal READS as advice, not
    as a wall: it must name the file route;
  * a named existing script is still checked as before (walking refuses,
    clean runs), with realistic option sets in every row;
  * `--wrap TEXT` keeps working, because its text is judgeable;
  * commands that are not a submit are untouched (`--version`, a pipe
    somewhere else, the word sbatch in ordinary text).
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho hello\n"

OPTION_SETS = [
    "",
    "--test-only",
    "--test-only -p dev_cpu",
    "--test-only -p dev_cpu -t 5",
    "-J probe -p dev_cpu",
]


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


def _opts(o):
    return f" {o}" if o else ""


# --------------------------------------------------------------------------
# refused: every stdin channel, CLEAN script included
# --------------------------------------------------------------------------

@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_clean_piped_script_is_refused(lustre, tmp_path, opts):
    """`cat clean.sh | sbatch` is refused: the file could be named instead."""
    script = tmp_path / "clean.sh"
    script.write_text(CLEAN, encoding="utf-8")
    verdict, text, may_run = _v(f"cat {script} | sbatch{_opts(opts)}")
    assert verdict == "refuse", (
        f"A CLEAN stdin submit must be refused under ruling (a); it was "
        f"{verdict!r} with options {opts!r}")
    assert may_run is False
    assert "job.sh" in text or "file" in text, (
        f"the refusal must tell the person to submit a FILE; it says: "
        f"{text[:200]}")


@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_clean_redirected_script_is_refused(lustre, tmp_path, opts):
    script = tmp_path / "clean.sh"
    script.write_text(CLEAN, encoding="utf-8")
    verdict, text, may_run = _v(f"sbatch{_opts(opts)} < {script}")
    assert verdict == "refuse", (
        f"`sbatch < file` must be refused; it was {verdict!r} with "
        f"options {opts!r}")
    assert may_run is False


@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_clean_heredoc_script_is_refused(lustre, opts):
    cmd = f"sbatch{_opts(opts)} <<'EOF'\n{CLEAN}EOF"
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"`sbatch <<EOF` must be refused; it was {verdict!r} with "
        f"options {opts!r}")
    assert may_run is False


@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_clean_here_string_script_is_refused(lustre, opts):
    cmd = f"sbatch{_opts(opts)} <<< '{CLEAN.strip()}'"
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"`sbatch <<< 'text'` must be refused; it was {verdict!r} with "
        f"options {opts!r}")
    assert may_run is False


@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_walking_stdin_script_is_refused(lustre, tmp_path, opts):
    """The original defect, still refused -- the reason must not regress."""
    script = tmp_path / "walk.sh"
    script.write_text(WALK, encoding="utf-8")
    verdict, text, may_run = _v(f"cat {script} | sbatch{_opts(opts)}")
    assert verdict == "refuse"
    assert may_run is False


def test_a_bare_sbatch_is_refused(lustre):
    """`sbatch` with nothing at all reads its script from stdin."""
    verdict, text, may_run = _v("sbatch")
    assert verdict == "refuse", (
        f"a bare `sbatch` reads stdin and must be refused; got {verdict!r}")
    assert may_run is False


# --------------------------------------------------------------------------
# still working: a named file, a wrap, and things that are not submits
# --------------------------------------------------------------------------

@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu -t 5"],
                         ids=["no_options", "with_options"])
def test_a_named_walking_file_is_still_refused(lustre, tmp_path, opts):
    script = tmp_path / "walk.sh"
    script.write_text(WALK, encoding="utf-8")
    verdict, text, may_run = _v(f"sbatch {opts} {script}".strip())
    assert verdict == "refuse"
    assert "tree_walk_on_lustre" in text
    assert may_run is False


@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu -t 5",
                                  "-J probe"],
                         ids=["no_options", "with_options", "with_job_name"])
def test_a_named_clean_file_still_runs(lustre, tmp_path, opts):
    script = tmp_path / "clean.sh"
    script.write_text(CLEAN, encoding="utf-8")
    verdict, text, may_run = _v(f"sbatch {opts} {script}".strip())
    assert verdict == "allow", (
        f"OVER-REFUSAL: a clean NAMED script was refused: {text[:200]}")
    assert may_run is True


@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu"],
                         ids=["no_options", "with_options"])
def test_a_clean_wrap_still_runs(lustre, opts):
    """`--wrap TEXT` is judgeable text, so it stays allowed."""
    verdict, text, may_run = _v(f"sbatch {opts} --wrap 'echo hello'".strip())
    assert verdict == "allow", (
        f"OVER-REFUSAL: a clean --wrap was refused: {text[:200]}")
    assert may_run is True


@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu"],
                         ids=["no_options", "with_options"])
def test_a_walking_wrap_is_still_refused(lustre, opts):
    verdict, text, may_run = _v(
        f"sbatch {opts} --wrap 'du -sh /pfs/data6/workspace/'".strip())
    assert verdict == "refuse"
    assert may_run is False


def test_sbatch_version_is_not_a_submit(lustre):
    verdict, text, may_run = _v("sbatch --version")
    assert verdict == "allow", (
        f"`sbatch --version` submits nothing and must stay allowed: "
        f"{text[:150]}")
    assert may_run is True


def test_an_ordinary_pipe_is_untouched(lustre):
    verdict, text, may_run = _v("du -sh /pfs/data6/workspace/ | head -1")
    assert verdict == "allow", (
        f"OVER-REFUSAL: an ordinary pipe was refused: {text[:150]}")


def test_a_clean_file_after_a_pipe_is_still_judged_by_the_file(lustre,
                                                               tmp_path):
    """`echo x | sbatch clean.sh` names the script as an ARGUMENT."""
    script = tmp_path / "clean.sh"
    script.write_text(CLEAN, encoding="utf-8")
    verdict, text, may_run = _v(f"echo hello | sbatch {script}")
    assert verdict == "allow", (
        f"OVER-REFUSAL: a named clean script behind a pipe was refused: "
        f"{text[:200]}")
    assert may_run is True
