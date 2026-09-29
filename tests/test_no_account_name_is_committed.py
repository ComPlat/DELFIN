"""A tracked file may not carry a personal account name.

The repository is public. A commit-message guard has enforced this for
message text since it was written -- and nothing read the files. That is
how `bench_trial/`, four SLURM scripts hardcoding one account's home
paths, reached a public repository as "evidence": the finding they
documented was already in the commit message, and what the files added
was the account name and the cluster's home layout.

The shape is checked, not a list of names. A cluster account here reads
`ka_xx0000`, so a guard on that shape lets a neutral example through --
`squeue -u someuser`, `/pfs/data6/home/xy/xy_group/xy_user/...` -- and
still refuses a real one. Tests that need a path to talk about keep
working; they just stop naming a person.

One exemption, named rather than pattern-matched: the shared bug archive
under `/home/qmchem_all`. It is a group account, not a person, and the
path is functional -- the archive is where it is, and blanking it would
break the feature that writes there. An exemption that is argued for is
different from one that is silently wide.
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

#: Named, with a reason. A group account for the shared bug archive; the
#: path is what makes the archive work.
_ALLOWED = ("qmchem_all",)

#: Extensions this guard does not try to read as text. Deliberately not
#: exhaustive: an unreadable file simply yields "" below, so a missing
#: entry costs nothing. ".npz" is left out on purpose -- the licence
#: guard reads that string in any file as a reference to CCDC data, and
#: arguing with another guard over a list that does not need the entry is
#: the wrong trade.
_BINARY = {".png", ".jpg", ".jpeg", ".gif", ".pdf", ".gbw", ".ico", ".woff",
           ".woff2", ".zip", ".gz", ".xlsx", ".docx", ".so"}


#: This file quotes what it refuses, so it would refuse itself. Excluded
#: by path rather than by some marker in the text: a marker would be a
#: way for any file to opt out.
_SELF = pathlib.Path(__file__).resolve()


def _tracked() -> list[pathlib.Path]:
    out = subprocess.run(["git", "ls-files"], cwd=str(_ROOT),
                         capture_output=True, text=True, timeout=120)
    return [p for p in (_ROOT / line for line in out.stdout.splitlines()
                        if line.strip())
            if p.resolve() != _SELF]


def _readable(path: pathlib.Path) -> str:
    if path.suffix.lower() in _BINARY:
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def test_git_is_available_or_this_checks_nothing():
    assert _tracked(), "git ls-files returned nothing; this guard is blind"


def test_no_tracked_file_names_a_personal_account():
    offenders: dict = {}
    for path in _tracked():
        text = _readable(path)
        if not text:
            continue
        found = {m for m in _ACCOUNT_RE.findall(text)
                 if not any(a in m for a in _ALLOWED)}
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

@pytest.mark.parametrize("sample", [
    "squeue -u ka_ew7404",
    "/pfs/data6/home/ka/ka_ibcs/ka_ew7404/software",
    "reported by ka_xn0397",
])
def test_the_guard_would_catch_what_was_removed(sample):
    assert _ACCOUNT_RE.search(sample), sample


@pytest.mark.parametrize("sample", [
    "squeue -u someuser",
    "/pfs/data6/home/xy/xy_group/xy_user/software",
    "/home/qmchem_all/archive/AGENT_BUGS",
    "delfin-agent credentials set KIT_TOOLBOX_API_KEY",
])
def test_a_neutral_stand_in_passes(sample):
    found = {m for m in _ACCOUNT_RE.findall(sample)
             if not any(a in m for a in _ALLOWED)}
    assert not found, f"{sample} is refused, so a test cannot name a path"


def test_the_exemption_is_a_group_not_a_person():
    """A list that grows quietly stops being a guard. One entry, and the
    docstring says why it is there."""
    assert _ALLOWED == ("qmchem_all",)
