"""IO2: an expansion anywhere in sbatch's argv is refused, not only in the
word the gate guessed was the script.

FOUND ON A LIVE RUN by the operator, after three rounds of my own tests were
green. Through the real bash tool and directly:

    sbatch <abs>/wal*.sh                              -> refuse   (correct)
    sbatch --test-only -p dev_cpu -t 5 <abs>/wal*.sh   -> ALLOW    <- FAIL-OPEN
    sbatch --test-only -p dev_cpu -t 5 <abs>/cle*.sh   -> ALLOW

THE ROOT, and it is the round-13 class one layer deeper. The gate elects ONE
word as "the script": the first positional that names an existing file, else
the first positional. An option's VALUE is a non-option token, so with
`-p dev_cpu -t 5` the positional list starts with `dev_cpu`, `existing_pos`
is empty (a glob names no file), and the script word becomes **`dev_cpu`**.
The expansion check then compares the glob word against `dev_cpu`, never
matches, and allows. Measured at `7ad808610`, and the same root lets the
other families through too:

    glob   + options -> allow      glob alone       -> refuse
    tilde  + options -> allow      (tilde alone was fixed at 2978fcf16)
    brace  + options -> allow
    plain walking script + options -> refuse   (the file rule still worked)
    clean script + options         -> allow    (must stay allow)

WHY THE FIX IS "EVERY WORD", not "the right word". Any attempt to elect the
script correctly has to model sbatch's option grammar -- which options take
a value, and which values are attached -- and that is the enumeration trap
this gate has already been bitten by four times (round 13: `-o`/`-e`
attached values; round 19: `-J name`). The structural fact needs no grammar:
if ANY word in sbatch's argv is one the SHELL would expand, then the argument
list sbatch receives is not the one on the command line, and the gate cannot
name the file it will run. So it refuses.

THE COST, stated because it is real: an unquoted glob in a word that is NOT
the script is also refused -- `sbatch -o out*.log job.sh`, or `-J name*`.
Bash would expand those too, so the argument list genuinely differs from
what was typed, and the refusal is cheap to obey (quote it). It is the same
trade the gate already makes elsewhere, and it is pinned below so the cost is
visible rather than discovered.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho done > /dev/null\n"

SHELL_EXPANSION = "pattern or expansion"

# Realistic submit lines: the sets a person actually types, including the
# ones the operator used on the live run.
OPTION_SETS = [
    "",
    "--test-only",
    "-p dev_cpu",
    "-t 5",
    "-J probe",
    "-p dev_cpu -t 5",
    "--test-only -p dev_cpu -t 5",
    "--test-only -p dev_cpu -t 5 -J probe",
    "--parsable -p dev_cpu",
]


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


@pytest.fixture
def jobs(tmp_path):
    d = tmp_path / "jobs"
    d.mkdir()
    (d / "walk.sh").write_text(WALK, encoding="utf-8")
    (d / "clean.sh").write_text(CLEAN, encoding="utf-8")
    return d


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


# --- direction 1: an expansion must be refused whatever the options are -----

@pytest.mark.parametrize("opts", OPTION_SETS,
                         ids=[o or "no_options" for o in OPTION_SETS])
def test_a_glob_is_refused_with_any_realistic_option_set(lustre, jobs, opts):
    """The operator's live rows, plus every option set a person types."""
    word = f"{jobs}/wal*.sh"
    cmd = f"sbatch {opts} {word}".strip()
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN: `{cmd}` was allowed; an option's value became the "
        f"script word and the glob was never judged")
    assert may_run is False
    assert SHELL_EXPANSION in text.lower(), (
        f"must refuse as an EXPANSION: {text[:250]}")


@pytest.mark.parametrize("opts", OPTION_SETS,
                         ids=[o or "no_options" for o in OPTION_SETS])
def test_a_clean_glob_is_refused_too(lustre, jobs, opts):
    """A pattern is never a named file, even when every match is clean."""
    cmd = f"sbatch {opts} {jobs}/cle*.sh".strip()
    verdict, _, may_run = _v(cmd)
    assert verdict == "refuse", f"a pattern is not a named file: {cmd}"
    assert may_run is False


@pytest.mark.parametrize("opts", OPTION_SETS,
                         ids=[o or "no_options" for o in OPTION_SETS])
def test_every_other_family_is_refused_with_options_too(lustre, jobs, opts):
    """The same root hid tilde and brace as well -- measured, not assumed."""
    for name, word in (("tilde", "~/walk.sh"),
                       ("brace", f"{jobs}/{{walk,clean}}.sh")):
        cmd = f"sbatch {opts} {word}".strip()
        verdict, text, _ = _v(cmd)
        assert verdict == "refuse", (
            f"FAIL-OPEN [{name}]: `{cmd}` was allowed")
        assert SHELL_EXPANSION in text.lower(), (
            f"[{name}] must refuse as an EXPANSION: {text[:250]}")


def test_an_expansion_in_an_option_value_is_refused_and_the_cost_is_visible(
        lustre, jobs):
    """The documented cost, pinned so it cannot be discovered by surprise.

    `-o out*.log` is not the script, but bash expands that word too, so the
    argument list sbatch receives differs from what was typed. Refused.
    """
    verdict, text, _ = _v(f"sbatch -o out*.log {jobs / 'clean.sh'}")
    assert verdict == "refuse", (
        "an unquoted glob in an option value is also expanded by the shell "
        "and must be refused; quote it to submit")
    assert SHELL_EXPANSION in text.lower(), text[:250]


# --- direction 2: no over-refusal on the ordinary work ---------------------

@pytest.mark.parametrize("opts", OPTION_SETS,
                         ids=[o or "no_options" for o in OPTION_SETS])
def test_a_clean_named_script_still_runs_with_any_option_set(
        lustre, jobs, opts):
    """The whole point: ordinary submits must keep working. Every option set
    here is one a person types, and none carries an expansion."""
    cmd = f"sbatch {opts} {jobs / 'clean.sh'}".strip()
    verdict, text, may_run = _v(cmd)
    assert (verdict, may_run) == ("allow", True), f"{cmd}: {text[:250]}"


@pytest.mark.parametrize("opts", OPTION_SETS,
                         ids=[o or "no_options" for o in OPTION_SETS])
def test_a_walking_named_script_still_refuses_for_the_walk(
        lustre, jobs, opts):
    """With options, the plain walking script must still refuse -- and for
    the WALK reason, so the expansion rule cannot be what passes this."""
    cmd = f"sbatch {opts} {jobs / 'walk.sh'}".strip()
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", f"{cmd} must refuse: {text[:200]}"
    assert may_run is False
    assert SHELL_EXPANSION not in text.lower(), (
        f"a plain script refuses for its tree walk, not as an expansion: "
        f"{text[:250]}")


def test_a_quoted_script_name_with_a_space_still_runs(lustre, tmp_path):
    """Quoting is the escape hatch the refusal tells people to use, so it
    must actually work end to end."""
    d = tmp_path / "spaced dir"
    d.mkdir()
    (d / "clean.sh").write_text(CLEAN, encoding="utf-8")
    verdict, text, may_run = _v(f'sbatch -p dev_cpu "{d}/clean.sh"')
    assert (verdict, may_run) == ("allow", True), text[:250]
