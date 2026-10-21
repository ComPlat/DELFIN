"""A tracked file may not carry a personal account name.

The repository is public. A commit-message guard has enforced this for
message text since it was written -- and nothing read the files. That is
how `bench_trial/`, four SLURM scripts hardcoding one account's home
paths, reached a public repository as "evidence": the finding they
documented was already in the commit message, and what the files added
was the account name and the cluster's home layout.

The shape is checked, not a list of names. A cluster account here is
the prefix `ka_`, two letters and four digits, so a guard on that shape
lets a neutral example through -- `squeue -u someuser`,
`/pfs/data6/home/xy/xy_group/xy_user/...` -- and still refuses a real
one. Tests that need a path to talk about keep working; they just stop
naming a person.

This file holds no account of that shape either: the samples the guard
must refuse are assembled at run time from parts that do not match on
their own. The file is therefore scanned like every other tracked file;
an exemption by path would let the one file that talks about accounts
carry one unseen.

Input: the text of every file `git ls-files` lists (binary suffixes
skipped). Output: the files that match, with up to three matches each.
Rejected: a list of real names (it would publish the names it guards);
a marker comment that opts a file out (any file could use it).
"""

from __future__ import annotations

import pathlib
import re
import subprocess

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]

#: This institution's accounts: the "ka_" prefix, then two letters and
#: four digits. Narrowed from a general two-letter prefix after that
#: version flagged `yy_gr0042` and `zz_ab1234` -- the deliberately
#: invented stand-ins inside the commit-message guard, which exist so a
#: test can talk about the rule without breaking it. A guard that refuses
#: the safe example teaches people to work around the guard.
_ACCOUNT_RE = re.compile(r"\bka_[a-z]{2}[0-9]{4}\b")

#: A real home under the cluster's parallel filesystem, with an account
#: in it. The neutral example paths this suite uses do not match.
_HOME_RE = re.compile(r"/pfs/data\d*/home/[a-z]{2}/[a-z]{2}_ka_[a-z]{2}[0-9]{4}\b"
                      r"|/pfs/data\d*/home/[a-z]{2}/[a-z]{2}_[a-z]+/ka_[a-z]{2}[0-9]{4}\b")

#: Extensions this guard does not try to read as text. Deliberately not
#: exhaustive: an unreadable file simply yields "" below, so a missing
#: entry costs nothing. ".npz" is left out on purpose -- the licence
#: guard reads that string in any file as a reference to CCDC data, and
#: arguing with another guard over a list that does not need the entry is
#: the wrong trade.
_BINARY = {".png", ".jpg", ".jpeg", ".gif", ".pdf", ".gbw", ".ico", ".woff",
           ".woff2", ".zip", ".gz", ".xlsx", ".docx", ".so"}


def _synthetic_account() -> str:
    """An account of the guarded shape that no person holds: letters
    "zz" and digits "0000". Joined here so the literal never appears in
    the file text."""
    return "".join(("ka", "_", "zz", "0000"))


def _tracked() -> list[pathlib.Path]:
    out = subprocess.run(["git", "ls-files"], cwd=str(_ROOT),
                         capture_output=True, text=True, timeout=120)
    return [_ROOT / line for line in out.stdout.splitlines() if line.strip()]


def _readable(path: pathlib.Path) -> str:
    if path.suffix.lower() in _BINARY:
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def test_git_is_available_or_this_checks_nothing():
    assert _tracked(), "git ls-files returned nothing; this guard is blind"


def test_this_file_is_scanned_too():
    """No self-exemption: the guard's own file is in the scanned set."""
    me = pathlib.Path(__file__).resolve()
    assert me in {p.resolve() for p in _tracked()}


def test_no_tracked_file_names_a_personal_account():
    offenders: dict = {}
    for path in _tracked():
        text = _readable(path)
        if not text:
            continue
        found = set(_ACCOUNT_RE.findall(text))
        if found:
            offenders[str(path.relative_to(_ROOT))] = sorted(found)[:3]
    assert not offenders, (
        "a public repository must not name people: "
        + "; ".join(f"{f} ({', '.join(n)})" for f, n in
                    sorted(offenders.items())[:8]))


def test_no_tracked_file_carries_a_real_cluster_home():
    offenders = []
    for path in _tracked():
        text = _readable(path)
        if text and _HOME_RE.search(text):
            offenders.append(str(path.relative_to(_ROOT)))
    assert not offenders, (
        "a home path with an account in it: " + ", ".join(offenders[:8]))


# -- the guard's own shape -------------------------------------------------

@pytest.mark.parametrize("template", [
    "squeue -u {acct}",
    "/pfs/data6/home/ka/ka_group/{acct}/software",
    "reported by {acct}",
])
def test_the_guard_would_catch_what_was_removed(template):
    sample = template.format(acct=_synthetic_account())
    assert _ACCOUNT_RE.search(sample), sample


def test_the_home_guard_would_catch_an_account_home():
    sample = f"/pfs/data6/home/ka/ka_group/{_synthetic_account()}/software"
    assert _HOME_RE.search(sample), sample


@pytest.mark.parametrize("sample", [
    "squeue -u someuser",
    "/pfs/data6/home/xy/xy_group/xy_user/software",
    "/home/<grp>/archive/AGENT_BUGS",
    "delfin-agent credentials set KIT_TOOLBOX_API_KEY",
])
def test_a_neutral_stand_in_passes(sample):
    assert not _ACCOUNT_RE.findall(sample), (
        f"{sample} is refused, so a test cannot name a path")
