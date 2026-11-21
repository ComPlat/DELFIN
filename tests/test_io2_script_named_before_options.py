"""IO2: a named script is judged even when a flag or an option follows it.

TWO findings, one root: the gate decides which word is the script by
looking for an OPTION-shaped token, so anything that makes sbatch's own
grammar disagree hides the script.

A. s2's F-R16 (`sbatch <walk.sh> --version`). hpc_io_check.py:1694-1697
   returned ``allow`` on ANY informational token BEFORE the script was
   examined. Measured at 43285ce79:

     sbatch <abs>/walk.sh --version   -> allow    <- the walk was never read
     sbatch <abs>/walk.sh --help      -> allow
     sbatch --version <abs>/walk.sh   -> allow
     sbatch <abs>/walk.sh             -> refuse

   s2's premise, "POSIX getopt stops at the first non-option, so the
   script IS that file", was the only reading under which this is a
   fail-open; the other failure mode is over-refusal. It is no longer a
   premise and no longer my guess -- the operator MEASURED it on a cluster
   (Slurm 25.11.8, harmless script, `--hold`, jobs cancelled at once):

     sbatch clean.sh --version  -> "Submitted batch job"   script submitted
     sbatch clean.sh --help     -> "Submitted batch job"   script submitted
     sbatch clean.sh --usage    -> "Submitted batch job"   script submitted
     sbatch --version clean.sh  -> version printed only,   no job

   So a flag AFTER the script name is that script's own argument and the
   script IS submitted: allowing it is a fail-open, and refusing is not
   over-reach. A flag BEFORE the script submits nothing, so refusing there
   is a deliberate, harmless asymmetry (ruling of the QS lead nacht-s8:
   nobody submits a script and asks for the version in one breath, and a
   position-keyed allow is one parser difference away from a walk). This
   is why the carve-out below keeps a flag-first command refused instead
   of allowing it: the shape is judged by whether a POSITIONAL NAMES AN
   EXISTING FILE, not by where the flag sits.

   An informational flag with NO script named still allows (`sbatch
   --version` is a version query and the ordinary way a person checks the
   tool); that case is pinned by its own test.

B. The option-VALUE case, found while probing A and not named by anyone:

     sbatch <abs>/walk.sh -J name     -> allow    <- the walk was never read

   The value `name` becomes positional[-1], the script test opens it, gets
   OSError, and the deliberate `except OSError` pass (:1652-1657) allows
   "so sbatch can refuse it". That reasoning is sound for a token that is
   genuinely the script; it is wrong when the script was already NAMED
   before the option. The first non-option word is the script (POSIX
   getopt), and it is judgeable.

Both cases are refused with the reason that names the walk, not with a
generic stdin message -- otherwise a future reader cannot tell a real
finding from the missing-name pass.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho done > /dev/null\n"


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


# --- A: an informational flag must not stop the script being judged ---------

INFO_AFTER = [
    ("version", "--version"),
    ("help", "--help"),
    ("usage", "--usage"),
    ("show_valid", "--show-valid"),
]


@pytest.mark.parametrize("name,flag", INFO_AFTER, ids=[c[0] for c in INFO_AFTER])
def test_a_flag_after_a_walking_script_does_not_hide_it(
        lustre, scripts, name, flag):
    verdict, text, may_run = _v(f"sbatch {scripts / 'walk.sh'} {flag}")
    assert verdict == "refuse", (
        f"FAIL-OPEN [{name}]: `sbatch walk.sh {flag}` was allowed; the "
        f"walking script was never read because of the informational "
        f"early return")
    assert may_run is False
    assert "tree walk" in text.lower() or "load" in text.lower(), (
        f"[{name}] the refusal must name the walk, not a generic reason: "
        f"{text[:250]}")


@pytest.mark.parametrize("name,flag", INFO_AFTER, ids=[c[0] for c in INFO_AFTER])
def test_a_flag_before_a_walking_script_does_not_hide_it(
        lustre, scripts, name, flag):
    """The MEASURED asymmetry, pinned deliberately.

    The operator measured (`sbatch --version clean.sh` -> version only, no
    job), so allowing this shape would be defensible. The QS lead ruled
    REFUSE anyway: the refusal costs nobody anything, and a position-keyed
    allow is one parser difference away from a walk. Pinned so a later
    "simplification" of the carve-out to `flag present -> allow` fails here
    rather than in production. All four spellings are rows, so a field
    that special-cases --version and --help cannot pass.
    """
    verdict, text, may_run = _v(f"sbatch {flag} {scripts / 'walk.sh'}")
    assert verdict == "refuse", (
        f"[{name}] `sbatch {flag} walk.sh` must stay refused by ruling "
        f"(it submits nothing, but the refusal is free): got {verdict!r}")
    assert may_run is False


def test_an_informational_flag_alone_still_allows(lustre):
    """The ordinary way a person checks the tool. Refusing this would be
    an over-refusal with no tree walk anywhere in sight."""
    for flag in ("--version", "--help", "--usage", "--show-valid"):
        verdict, text, may_run = _v(f"sbatch {flag}")
        assert (verdict, may_run) == ("allow", True), (
            f"`sbatch {flag}` is a version/help query, not a submit: "
            f"{verdict!r} / {text[:200]}")


def test_a_clean_script_with_a_flag_still_allows(lustre, scripts):
    """No over-refusal on the clean case: the flag is not itself a finding."""
    verdict, text, may_run = _v(f"sbatch {scripts / 'clean.sh'} --version")
    assert (verdict, may_run) == ("allow", True), text


# --- B: an option's value must not become the script ------------------------

OPTION_FORMS = [
    ("space_form", "-J name"),
    ("partition", "-p dev_cpu"),
    ("time", "-t 5"),
    ("long_form", "--job-name=name"),
]


@pytest.mark.parametrize("name,opt", OPTION_FORMS,
                         ids=[c[0] for c in OPTION_FORMS])
def test_an_option_value_after_the_script_does_not_hide_it(
        lustre, scripts, name, opt):
    verdict, text, may_run = _v(f"sbatch {scripts / 'walk.sh'} {opt}")
    assert verdict == "refuse", (
        f"FAIL-OPEN [{name}]: `sbatch walk.sh {opt}` was allowed; the "
        f"option's value became the script and the real script was never "
        f"read")
    assert may_run is False


def test_a_clean_script_with_an_option_still_allows(lustre, scripts):
    verdict, text, may_run = _v(f"sbatch {scripts / 'clean.sh'} -J name")
    assert (verdict, may_run) == ("allow", True), text
