"""IO2: the informational-flag carve-out must fire ONLY when no script is named.

Operator + QS ruling (round 24, #859 held): the carve-out at
hpc_io_check.py:1847 read

    if info_flag and not existing_pos:
        return "allow", "", True

`existing_pos` lists positional words that NAME AN EXISTING FILE. A word the
shell would EXPAND never names one, so a pattern beside a flag fell straight
into the allow:

    sbatch <dir>/wal*.sh --version   -> allow   <- walk.sh would run
    sbatch <dir>/walk.sh  --version  -> refuse  (correct)

That is the same class as round 22: the rule asked "does this word name a
file" when the question is "could this word BE the script". A word that
expands cannot be resolved by the gate, so it must not satisfy a condition
whose purpose is "there is no script here".

THE MEASURED PREMISE (operator run, Slurm 25.11.8, harmless script, --hold,
jobs cancelled at once) -- a flag AFTER the script name is the script's own
argument and the script IS submitted:

    sbatch clean.sh --version   -> "Submitted batch job"
    sbatch clean.sh --help      -> "Submitted batch job"
    sbatch clean.sh --usage     -> "Submitted batch job"
    sbatch --version clean.sh   -> version printed only, no job

So a flag does not make a named script stop being submitted. The carve-out
therefore exists for ONE shape only -- an sbatch invocation that names no
script at all -- and any positional word, literal or expanding, must switch
it off and let the normal rules judge.

WHY THE ROWS USE REAL OPTION SETS: the round-22 fail-open survived because
every test wrote the script word first. A person types
`sbatch -p dev_cpu -t 5 wal*.sh`; the option VALUES are non-option tokens, so
tests here carry option values on purpose.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho hello\n"

SHELL_PATTERN = "is a shell pattern"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


@pytest.fixture
def scripts(tmp_path):
    d = tmp_path / "jobs"
    d.mkdir()
    (d / "walk.sh").write_text(WALK, encoding="utf-8")
    (d / "clean.sh").write_text(CLEAN, encoding="utf-8")
    return d


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


# --- direction 1: a script that EXPANDS beside a flag is refused ------------
#
# Every row is a shape that reached "allow" at f4d3e2984. The flag is
# informational in each, which is what made the carve-out fire.

FLAG_ROWS = [
    ("glob_version", "wal*.sh", "--version"),
    ("glob_help", "wal*.sh", "--help"),
    ("glob_usage", "wal*.sh", "--usage"),
    ("glob_show_valid", "wal*.sh", "--show-valid"),
    ("question_version", "wal?.sh", "--version"),
    ("bracket_help", "wal[k].sh", "--help"),
]


@pytest.mark.parametrize(
    "name,word,flag", FLAG_ROWS, ids=[c[0] for c in FLAG_ROWS]
)
def test_an_expanding_script_beside_a_flag_is_refused(lustre, scripts, name,
                                                      word, flag):
    """The carve-out must not outrank the expansion rule.

    `--show-valid` is in the flag list because it is in the carve-out's
    tuple; s8 could not confirm it against sbatch's documentation, so it is
    exercised as a carve-out member and nothing is claimed about what sbatch
    does with it.
    """
    verdict, text, may_run = _v(f"sbatch {scripts / word} {flag}")
    assert verdict == "refuse", (
        f"FAIL-OPEN [{name}]: sbatch {word!r} {flag} was allowed; bash "
        f"expands the word before sbatch sees it")
    assert may_run is False
    assert SHELL_PATTERN in text.lower(), (
        f"[{name}] the reason must be the pattern rule, not another "
        f"refusal: {text[:200]}")


def test_a_tilde_script_beside_a_flag_is_refused(lustre, tmp_path):
    """`~/walk.sh --version` — the tilde expands, so it names no file here."""
    home = tmp_path / "home"
    home.mkdir()
    (home / "walk.sh").write_text(WALK, encoding="utf-8")
    env = dict(ENV, HOME=str(home))
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        "sbatch ~/walk.sh --version", env=env)
    assert verdict == "refuse", (
        "FAIL-OPEN [tilde_flag]: allowed; bash expands ~ before sbatch "
        "sees the word")
    assert may_run is False


# --- direction 2: the carve-out's own case must keep working ---------------


def test_a_flag_with_no_script_at_all_is_allowed(lustre):
    """`sbatch --version` is how a person checks the tool; nothing to judge.

    This is the ONLY shape the carve-out is for. If the fix removed it
    wholesale instead of narrowing it, this row goes red and the gate starts
    refusing an ordinary harmless command.
    """
    for cmd in ("sbatch --version", "sbatch --help", "sbatch --usage"):
        verdict, text, may_run = _v(cmd)
        assert verdict == "allow", (
            f"OVER-REFUSAL: {cmd!r} was refused ({text[:160]}); it names no "
            f"script and submits nothing")
        assert may_run is True


def test_a_clean_script_beside_a_flag_still_runs(lustre, scripts):
    """The measured table's allow rows: a flag after a clean script submits
    that clean script, which DELFIN must keep permitting."""
    for cmd in (
        f"sbatch {scripts / 'clean.sh'} --version",
        f"sbatch {scripts / 'clean.sh'} --help",
        f"sbatch -p x {scripts / 'clean.sh'} --version",
        f"sbatch -p x -t 5 -J job {scripts / 'clean.sh'} --version",
    ):
        verdict, text, may_run = _v(cmd)
        assert verdict == "allow", (
            f"OVER-REFUSAL: {cmd!r} was refused ({text[:160]}); the script "
            f"is clean and the flag does not stop sbatch submitting it")


def test_a_clean_script_beside_a_flag_with_options_still_runs(lustre,
                                                              scripts):
    """The operator's third row, spelled out: option VALUES must not be
    mistaken for the script, and must not switch the carve-out back on in a
    way that skips the clean script."""
    verdict, text, may_run = _v(
        f"sbatch -p dev_cpu -t 5 {scripts / 'clean.sh'} --version")
    assert verdict == "allow", f"OVER-REFUSAL: {text[:200]}"
    assert may_run is True
