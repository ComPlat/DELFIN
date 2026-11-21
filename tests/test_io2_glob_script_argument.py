"""IO2: a glob in the script argument is not a script the gate can name.

QS lead nacht-s8's finding (round 18, `651a3c15d` era): `sbatch wal*.sh`
was ALLOWED. The word is a PATTERN, so ``_gate_names_existing_file``
(hpc_io_check.py:1297) finds no file for it, and the deliberate
missing-name pass (:1609-1612 -> the narrow "not found" rule) let it
through. Bash expands the pattern BEFORE sbatch sees an argument, so the
file that actually runs is whichever match sorts first -- and nothing
inspected it. Measured at `43285ce79` (the probe that reproduced it,
`tests/test_zz_io2_glob_probe.py`):

    sbatch <abs>/wal*.sh    verdict='allow'     <- walking script would run
    sbatch <abs>/walk.sh    verdict='refuse'

RULING (s8, boundary not spelling): a script argument whose glob
metacharacter is UNQUOTED is refused -- the gate cannot name the file
bash will choose, so it cannot judge it. The metacharacter set is `*`,
`?`, `[`; quoting is closed under POSIX quoting, so `"wal*.sh"` names one
literal file and stays on the existing path.

WHY THIS FILE DOES NOT TRUST shlex's TOKEN LIST: the quotes are stripped
by then. Measured (`tests/test_zz_io2_quote_survival.py`):

    shlex.split('sbatch "/a/wal*.sh"')  -> ['sbatch', '/a/wal*.sh']
    shlex.split('sbatch /a/wal*.sh')    -> ['sbatch', '/a/wal*.sh']

Those are the same token, so the quoted/unquoted decision has to come
from the RAW text (:1596-1598's `last[:1] in ("'", '"')` only sees a
quote in the `ValueError` fallback at :1524-1527, where `str.split()`
keeps them). This file pins both directions of that decision:
unquoted -> refused, quoted/escaped -> NOT the glob refusal.

WHY THE ASSERTIONS READ `SHELL_PATTERN` AND NOT `"pattern"`: pytest's
`tmp_path` embeds the TEST NAME, and this file's test names contain the
word "pattern" -- so `assert "pattern" in text` was satisfied by the
path inside the failure message and went green for the wrong reason
(measured: the cd case below was one of the six green in the red run).
`SHELL_PATTERN` is a distinctive substring of the refusal clause that no
path can contain.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho hello\n"

# A substring of the refusal clause that no temporary path can contain.
SHELL_PATTERN = "is a shell pattern"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


@pytest.fixture
def scripts(tmp_path):
    """walk.sh and clean.sh in one directory, nothing else."""
    d = tmp_path / "jobs"
    d.mkdir()
    (d / "walk.sh").write_text(WALK, encoding="utf-8")
    (d / "clean.sh").write_text(CLEAN, encoding="utf-8")
    return d


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


# --- direction 1: an unquoted pattern must be refused -----------------------

PATTERNS = [
    ("star", "wal*.sh"),
    ("question", "wal?.sh"),
    ("bracket", "wal[k].sh"),
    ("star_clean_target", "cle*.sh"),
    ("star_prefix", "*.sh"),
]


@pytest.mark.parametrize("name,word", PATTERNS, ids=[c[0] for c in PATTERNS])
def test_an_unquoted_pattern_is_refused(lustre, scripts, name, word):
    """The gate cannot name the file bash will pick, so it refuses.

    `cle*.sh` is in the list on purpose: a pattern is NEVER a named file,
    even when every match happens to be clean. Judging one guessed match
    would be a green that says nothing about the others.
    """
    verdict, text, may_run = _v(f"sbatch {scripts / word}")
    assert verdict == "refuse", (
        f"FAIL-OPEN [{name}]: sbatch {word!r} was allowed; bash expands the "
        f"pattern and the gate never saw the file that runs")
    assert may_run is False
    assert "job.sh" in text, (
        f"[{name}] the refusal must name the file route: {text[:200]}")
    assert SHELL_PATTERN in text.lower(), (
        f"[{name}] the reason must be the PATTERN rule: {text[:200]}")


def test_the_refusal_names_the_pattern(lustre, scripts):
    """s8's ruling: the refusal carries a clause naming the pattern, so a
    person is told which word to replace -- not just that something
    failed."""
    verdict, text, _ = _v(f"sbatch {scripts / 'wal*.sh'}")
    assert verdict == "refuse"
    assert "wal*.sh" in text, f"the pattern must be named: {text[:250]}"
    assert SHELL_PATTERN in text.lower(), (
        f"the clause must be the pattern rule: {text[:250]}")


def test_a_relative_pattern_resolves_against_the_cd_target(lustre, scripts):
    """`cd <dir> && sbatch wal*.sh`: bash expands in the cd directory, so
    the pattern check must not depend on the process cwd either.

    The ASSERTION IS THE REASON, not the verdict. Without the glob rule
    this command already refuses -- through the cd path's "not found from
    the cd target" branch (:1638-1644), because `wal*.sh` is no file. A
    verdict-only assertion would therefore be green before the fix and
    prove nothing (measured: it was one of the 6 green in the red run).
    """
    verdict, text, may_run = _v(f"cd {scripts} && sbatch wal*.sh")
    assert verdict == "refuse", "a relative pattern after cd is still a pattern"
    assert may_run is False
    assert SHELL_PATTERN in text.lower(), (
        "the refusal must be the PATTERN rule, not the cd-path "
        f"'not found' branch that would fire anyway: {text[:250]}")


# --- direction 2: these must NOT be refused as patterns ---------------------

def test_a_plain_walking_script_still_refuses_for_the_walk(
        lustre, scripts):
    """The glob rule must not become the reason a plain script refuses --
    otherwise the glob change could pass this file by shadowing the real
    finding."""
    verdict, text, may_run = _v(f"sbatch {scripts / 'walk.sh'}")
    assert verdict == "refuse"
    assert may_run is False
    assert SHELL_PATTERN not in text.lower(), (
        f"a plain script must refuse for its tree walk, not as a pattern: "
        f"{text[:250]}")


def test_a_plain_clean_script_is_allowed(lustre, scripts):
    """No over-refusal: the ordinary case is untouched."""
    verdict, text, may_run = _v(f"sbatch -p dev_cpu {scripts / 'clean.sh'}")
    assert (verdict, may_run) == ("allow", True), text


def test_an_escaped_star_names_a_literal_file_and_is_judged(
        lustre, tmp_path):
    """`walk\\*.sh` is a LITERAL name to bash (measured: `printf '%s\\n'
    /a/\\*.sh` prints `/a/*.sh`). So when such a file exists, the word is
    not a pattern -- and the file it names is judged like any other."""
    d = tmp_path / "literal"
    d.mkdir()
    (d / "walk*.sh").write_text(WALK, encoding="utf-8")
    verdict, text, may_run = _v(f"sbatch {d}/walk\\*.sh")
    assert verdict == "refuse", (
        "an escaped metacharacter names a real walking file, which must be "
        f"judged; got {verdict!r} / {text[:200]}")
    assert may_run is False
    assert SHELL_PATTERN not in text.lower(), (
        f"an escaped star is not a pattern: {text[:250]}")


def test_a_quoted_pattern_with_a_literal_file_of_that_name_is_judged(
        lustre, tmp_path):
    """s8's ruling, the quoted half: `"walk*.sh"` names ONE literal file.
    When it exists it is judged -- and it walks, so it refuses."""
    d = tmp_path / "quoted"
    d.mkdir()
    (d / "walk*.sh").write_text(WALK, encoding="utf-8")
    verdict, text, may_run = _v(f'sbatch "{d}/walk*.sh"')
    assert verdict == "refuse", (
        f'the literal file "walk*.sh" walks; got {verdict!r} / {text[:200]}')
    assert may_run is False


def test_a_quoted_pattern_is_not_refused_as_a_pattern(lustre, scripts):
    """The other quoted half: with NO file of that literal name, bash hands
    sbatch the literal string and sbatch fails on its own ("No such
    file"). That is not a fail-open, and s8's ruling keeps it on the
    existing path -- so the GLOB refusal must not fire.

    Asserted on the text, not the verdict, because the existing path's
    answer for a missing name is a decision of its own (rounds 6-9) and is
    pinned by its own test; this file only says the glob rule stays out of
    it."""
    verdict, text, _ = _v(f'sbatch "{scripts}/wal*.sh"')
    assert SHELL_PATTERN not in text.lower(), (
        f"a quoted word is a literal name, not a pattern: {text[:250]}")
