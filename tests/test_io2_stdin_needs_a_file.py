"""IO2: an sbatch with no script FILE argument must not be fed a script
through stdin that nobody judged -- whatever its option set.

Wave 33, round 13. The operator's live run found the round-12 fix
incomplete: stdin-fed scripts were judged only when `sbatch` carried NO
options. Reproduced here first, against the real gate:

    cat walk.sh | sbatch --test-only -p dev_cpu -t 5   -> allow  (FAIL-OPEN)
    sbatch --test-only -p dev_cpu -t 5 < walk.sh       -> allow  (FAIL-OPEN)
    sbatch < walk.sh                                   -> allow  (FAIL-OPEN)
    cat walk.sh | sbatch                               -> refuse (judged)

Mechanism, and it is the same "one spelling" class as rounds 10-12:
`_sbatch_reads_stdin` skips every token that starts with `-` as an option,
so in `sbatch --test-only -p dev_cpu -t 5` the token `dev_cpu` is the first
token that does not start with `-`, and it is therefore taken to be a
script argument -- so the segment "names a script", the stdin source is
discarded and nothing is judged. For `< walk.sh` the bare `<` is skipped by
the same loop and `walk.sh` is again taken as a positional argument.

This test pins ONLY the failure direction both candidate designs share: a
script reaching sbatch through stdin, with the script text walking Lustre,
must not run. It deliberately does NOT yet pin what happens to a CLEAN
script fed through stdin -- whether that is refused outright (the
operator's proposed fail-closed rule) or judged and allowed (a correct
stdin-source parser) is a pending ruling, and pinning either now would
decide it by accident.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"

# Realistic option sets: a submit line a user actually types. `--wrap` is
# absent on purpose -- a wrap IS the script text and is judged separately.
OPTION_SETS = [
    "",
    "--test-only",
    "--test-only -p dev_cpu",
    "--test-only -p dev_cpu -t 5",
    "-p dev_cpu --test-only -t 5 -J probe",
    "-J probe",
    "--parsable -p dev_cpu",
]


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


def _walking_script(tmp_path):
    script = tmp_path / "walk.sh"
    script.write_text(WALK, encoding="utf-8")
    return script


@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_piped_walking_script_is_refused_with_any_option_set(
        lustre, tmp_path, opts):
    """`cat walk.sh | sbatch [options]` must never reach sbatch unchecked."""
    script = _walking_script(tmp_path)
    verdict, text, may_run = _v(f"cat {script} | sbatch {opts}".strip())
    assert verdict == "refuse", (
        f"FAIL-OPEN: a piped walking script got through with options "
        f"{opts!r} (verdict={verdict!r})")
    assert may_run is False


@pytest.mark.parametrize("opts", OPTION_SETS, ids=[o or "no_options"
                                                   for o in OPTION_SETS])
def test_a_redirected_walking_script_is_refused_with_any_option_set(
        lustre, tmp_path, opts):
    """`sbatch [options] < walk.sh` is the same script, different spelling."""
    script = _walking_script(tmp_path)
    verdict, text, may_run = _v(f"sbatch {opts} < {script}".strip())
    assert verdict == "refuse", (
        f"FAIL-OPEN: a redirected walking script got through with options "
        f"{opts!r} (verdict={verdict!r})")
    assert may_run is False


@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu -t 5"],
                         ids=["no_options", "with_options"])
def test_a_walking_script_in_a_here_string_is_refused(lustre, opts):
    """The ANSI-C here-string form the live run used."""
    cmd = ("sbatch " + (opts + " " if opts else "")
           + "<<< $'" + WALK.strip() + "'")
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN: here-string script through stdin was allowed "
        f"(verdict={verdict!r}) with options {opts!r}")
    assert may_run is False


# --------------------------------------------------------------------------
# must still run: a named FILE is checked as before, and a clean one runs
# --------------------------------------------------------------------------

@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu -t 5"],
                         ids=["no_options", "with_options"])
def test_a_named_walking_file_is_still_refused(lustre, tmp_path, opts):
    script = _walking_script(tmp_path)
    verdict, text, may_run = _v(f"sbatch {opts} {script}".strip())
    assert verdict == "refuse"
    assert "tree_walk_on_lustre" in text
    assert may_run is False


@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu -t 5"],
                         ids=["no_options", "with_options"])
def test_a_named_clean_file_still_runs(lustre, tmp_path, opts):
    script = tmp_path / "clean.sh"
    script.write_text("#!/bin/bash\necho hello\n", encoding="utf-8")
    verdict, text, may_run = _v(f"sbatch {opts} {script}".strip())
    assert verdict == "allow", (
        f"OVER-REFUSAL: a clean named script was refused: {text[:120]}")
    assert may_run is True


@pytest.mark.parametrize("opts", ["", "--test-only -p dev_cpu -t 5"],
                         ids=["no_options", "with_options"])
def test_a_clean_wrap_still_runs(lustre, opts):
    """`--wrap TEXT` is script text and must keep working."""
    verdict, text, may_run = _v(
        f"sbatch {opts} --wrap 'echo hello'".strip())
    assert verdict == "allow", (
        f"OVER-REFUSAL: a clean --wrap was refused: {text[:120]}")
    assert may_run is True
