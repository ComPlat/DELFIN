"""IO2: a script argument the SHELL would expand is refused.

s8's round-18 ruling (binding, and it supersedes the narrower glob rule I
committed at d08d67ceb): a word the shell expands before sbatch sees it is
refused, because the script that runs is not the word on the command line
and the gate cannot name it. One structural rule, not a fifth spelling --
and explicitly NO simulation of the expansion, since a gate that
re-implements globbing is the same enumeration trap in a new costume.

MEASURED at the previous commit (d08d67ceb), oracle = a real bash with the
word spliced UNQUOTED into `printf '%s\\n' WORD`, so the oracle itself is
checked (my first version passed the word via "$@", which quotes it and
therefore never expands -- that bug made four fail-opens look clean):

  word                       bash hands over              gate said
  /x/wal*.sh                 /x/walk.sh                   refuse (fixed in d08d67ceb)
  ~/walk.sh                  /home/<user>/walk.sh         ALLOW  <- fail-open
  /x/{walk,clean}.sh         /x/walk.sh /x/clean.sh       ALLOW  <- fail-open
  $(echo /x/walk.sh)         /x/walk.sh                   ALLOW  <- fail-open
  "/x/wal*.sh"               "/x/wal*.sh" (literal)       allow  <- correct
  /x/walk.sh                 /x/walk.sh                   refuse <- correct

THE BOUNDARIES ARE MEASURED TOO, because a trigger that is too broad
refuses ordinary submits. tests/test_zz_io2_expand_bounds.py, real bash:

  /x/{a}.sh    -> /x/{a}.sh    (literal)   so `{` alone must NOT trigger
  /x/{}        -> /x/{}        (literal)
  /x/{1..3}.sh -> 1,2,3        (range)     so `..` inside braces triggers
  /x/a{b,c}d.sh-> two words               so brace fires MID-word
  /x/~/a.sh    -> literal                  so only a LEADING `~` triggers
  a~b.sh       -> literal
  "$(echo hi)" -> hi                       so `$(` triggers EVEN QUOTED
  '$(echo hi)' -> literal                  single quotes protect it
  \\$(echo hi) -> literal                  an escape protects it

The last one is why command substitution is not covered by the quoting
exemption the other three get; assuming it was would have been a fail-open
inside the fix for a fail-open.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho done > /dev/null\n"

# A distinctive substring of the refusal clause. It must not be a bare word
# that pytest's tmp_path can supply: tmp_path embeds the TEST NAME, and an
# assertion on a word appearing in the name is vacuously true (measured in
# the phase-1 control, where `assert "pattern" in text` went green because
# the path contained "..._pattern_...").
SHELL_EXPANSION = "pattern or expansion"


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


# --- direction 1: every expansion family the shell applies is refused --------

EXPANSIONS = [
    ("glob_star", "wal*.sh", "glob"),
    ("glob_question", "wal?.sh", "glob"),
    ("glob_bracket", "wal[k].sh", "glob"),
    ("brace_list", "{walk,clean}.sh", "brace"),
    ("brace_range", "{walk,clean}.s{h,h}", "brace"),
    ("command_sub", "$(echo walk.sh)", "command"),
    ("backtick", "`echo walk.sh`", "command"),
]


@pytest.mark.parametrize("name,word,kind", EXPANSIONS,
                         ids=[c[0] for c in EXPANSIONS])
def test_a_word_the_shell_expands_is_refused(lustre, jobs, name, word, kind):
    """Each row names the file that would run in `jobs/`."""
    verdict, text, may_run = _v(f"sbatch {jobs / word}")
    assert verdict == "refuse", (
        f"FAIL-OPEN [{name}]: `sbatch {word}` was allowed; bash expands the "
        f"word and the gate never saw the file that runs")
    assert may_run is False
    assert SHELL_EXPANSION in text.lower(), (
        f"[{name}] the refusal must be the expansion rule: {text[:250]}")
    assert kind in text.lower(), (
        f"[{name}] the refusal must name the kind ({kind}): {text[:250]}")


def test_a_leading_tilde_is_refused(lustre):
    """`~/walk.sh` expands to $HOME/walk.sh before sbatch sees it."""
    verdict, text, may_run = _v("sbatch ~/walk.sh")
    assert verdict == "refuse", (
        "FAIL-OPEN: a leading tilde was allowed; bash expands it to "
        "$HOME/walk.sh first")
    assert may_run is False
    assert "tilde" in text.lower(), f"the kind must be named: {text[:250]}"


def test_a_quoted_command_substitution_is_still_refused(lustre, jobs):
    """The boundary that separates this family: bash expands `$(...)`
    inside DOUBLE quotes too, so quoting must not rescue it (measured)."""
    verdict, text, may_run = _v(f'sbatch "{jobs}/$(echo walk.sh)"')
    assert verdict == "refuse", (
        "FAIL-OPEN: command substitution expands inside double quotes, so a "
        "quoted word is not literal here")
    assert may_run is False
    assert "command" in text.lower(), f"the kind must be named: {text[:250]}"


def test_a_relative_expansion_after_a_cd_is_refused(lustre, jobs):
    """The pattern check must not depend on the process cwd."""
    verdict, text, may_run = _v(f"cd {jobs} && sbatch wal*.sh")
    assert verdict == "refuse", "a relative pattern after cd is still a pattern"
    assert may_run is False
    assert SHELL_EXPANSION in text.lower(), (
        "the refusal must be the EXPANSION rule, not the cd path's "
        f"'not found' branch that would fire anyway: {text[:250]}")


# --- direction 2: these must NOT be refused as expansions -------------------


def test_a_quoted_glob_is_literal_and_goes_to_the_file_judgement(
        lustre, tmp_path):
    """Measured: bash hands `"/x/wal*.sh"` over unchanged. A file of that
    literal name walks, so the FILE rule refuses it -- not the glob rule."""
    d = tmp_path / "literal"
    d.mkdir()
    (d / "walk*.sh").write_text(WALK, encoding="utf-8")
    verdict, text, may_run = _v(f'sbatch "{d}/walk*.sh"')
    assert verdict == "refuse", (
        f"the literal file 'walk*.sh' walks; got {verdict!r} / {text[:200]}")
    assert may_run is False
    assert "tree" in text.lower() or "load" in text.lower(), (
        f"a quoted glob is judged as a FILE, not as an expansion: "
        f"{text[:250]}")


@pytest.mark.parametrize("word", ["{a}.sh", "{}", "a~b.sh"],
                         ids=["single_brace", "empty_brace", "mid_tilde"])
def test_a_word_the_shell_does_not_expand_is_not_refused(lustre, jobs, word):
    """Measured boundaries: a bare `{x}` and a mid-word `~` stay literal,
    so a trigger that fired on them would over-refuse an ordinary submit.
    These name no file, so the existing missing-name decision applies; what
    this pins is that the EXPANSION rule stays out of it."""
    verdict, text, _ = _v(f"sbatch {jobs / word}")
    assert SHELL_EXPANSION not in text.lower(), (
        f"[{word}] the shell does not expand this word, so the expansion "
        f"rule must not fire: {text[:250]}")
    assert verdict != "refuse" or "tree" not in text.lower()


def test_a_plain_walking_script_still_refuses_for_the_walk(lustre, jobs):
    """The expansion rule must not become the reason a plain script
    refuses -- otherwise it could pass this file by shadowing the real
    finding."""
    verdict, text, may_run = _v(f"sbatch {jobs / 'walk.sh'}")
    assert verdict == "refuse"
    assert may_run is False
    assert SHELL_EXPANSION not in text.lower(), (
        f"a plain script must refuse for its tree walk: {text[:250]}")


def test_a_plain_clean_script_is_allowed(lustre, jobs):
    """No over-refusal on the ordinary case, with and without options."""
    for cmd in (f"sbatch {jobs / 'clean.sh'}",
                f"sbatch -p dev_cpu -t 5 {jobs / 'clean.sh'}",
                f"sbatch {jobs / 'clean.sh'} --test-only"):
        verdict, text, may_run = _v(cmd)
        assert (verdict, may_run) == ("allow", True), f"{cmd}: {text[:200]}"


# --- s8's point (2): the FIRST positional is the script ---------------------

def test_several_positionals_are_judged_on_the_first(lustre, jobs):
    """`sbatch walk.sh clean.sh`: the FIRST positional is the script (POSIX
    getopt, and sbatch's own `sbatch [options] script [args]` usage), so the
    walking file must be judged even though a clean one follows it.

    Measured before the fix: verdict='allow' -- the gate read the LAST
    positional, found clean.sh, and never looked at walk.sh.
    """
    verdict, text, may_run = _v(
        f"sbatch {jobs / 'walk.sh'} {jobs / 'clean.sh'}")
    assert verdict == "refuse", (
        "FAIL-OPEN: the FIRST positional is the script and it walks; the "
        "clean second positional must not hide it")
    assert may_run is False


def test_a_clean_first_positional_with_an_argument_allows(lustre, jobs):
    """The other direction: a script argument after the script is that
    script's argument, not a second script to judge."""
    verdict, text, may_run = _v(
        f"sbatch {jobs / 'clean.sh'} {jobs / 'walk.sh'}")
    assert (verdict, may_run) == ("allow", True), text
