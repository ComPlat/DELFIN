"""No tracked file carries a colleague's name, a site host name, a site
address or a real home path.

Data minimisation (GDPR Art. 5(1)(c)) for a public repository: a field
report is cited by its date and its technical content, never by who
filed it; a test that needs a host, an address or a home uses a
stand-in. `test_no_account_name_is_committed.py` covers the account
shape; this file covers the rest.

Rules, each over the text of every file `git ls-files` lists:

- ``node``: the site's node names (``uc3n`` and digits) -- login and
  compute nodes alike. Stand-ins: ``login01``, ``node024``.
- ``site_ip``: an address in the institution's /16 block (the two
  leading octets are in `_PATTERNS`).
  Stand-in: RFC 5737 documentation addresses (192.0.2.x).
- ``institution_mail``: an institutional mail of the account shape
  (two letters, four digits, at a kit.edu domain).
- ``session_scratch``: a session scratch path that encodes a home
  directory (``/tmp/<tool>-<uid>/-home-<name>-...``).
- ``home``: a top-level ``/home/<name>`` whose name is not one of the
  stand-ins in `_HOME_STAND_INS`. A new test picks one of those.
- ``name``: a word whose salted hash is in `_NAME_HASHES` -- colleagues'
  first names and the site's group-account names. Ambiguous names (one
  is also a month abbreviation) count only in the possessive form. The
  names are not written here; the hashes keep them out of the text, they
  do not keep them secret from someone who guesses a first name.

Allow-list: the name rule skips `_ATTRIBUTION_FILES` (the citation and
the README's reference list name authors on purpose). Nothing else is
exempt, this file included.

Rejected: a plain deny-list of names (it would publish what it guards);
an allow-list of files per rule (it grows quietly); a generic IPv4 rule
(public resolvers and documentation addresses are legitimate in tests).
Measured cost: one pass over the tracked tree (64 MB of text), 2.9 s.
"""

from __future__ import annotations

import hashlib
import pathlib
import re
import subprocess

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]

_BINARY = {".png", ".jpg", ".jpeg", ".gif", ".pdf", ".gbw", ".ico", ".woff",
           ".woff2", ".zip", ".gz", ".xlsx", ".docx", ".so"}

#: rule -> (literal every match contains, pattern). The literal is a
#: substring test that skips the regex on most files.
_PATTERNS: dict[str, tuple[str, re.Pattern]] = {
    "node": ("uc3n", re.compile(r"\buc3n\d+\b")),
    "site_ip": ("141.52.", re.compile(r"\b141\.52\.\d{1,3}\.\d{1,3}\b")),
    "institution_mail": ("kit.edu", re.compile(
        r"\b[a-z]{2}\d{4}@(?:[\w-]+\.)*kit\.edu\b")),
    "session_scratch": ("-home-", re.compile(r"/tmp/[A-Za-z]+-\d+/-home-")),
}

_HOME_RE = re.compile(r"(?<![\w.~-])/home/([A-Za-z_][A-Za-z0-9_.-]*)")

#: The stand-ins tests already use for a home directory. A new one is
#: added here with the test that needs it, never a real account.
_HOME_STAND_INS = frozenset({
    "a", "ab", "alice", "dev", "exampleuser", "fixtureacct", "foo", "grp",
    "jdoe42", "ka", "localuser", "me", "other", "other_user", "someone",
    "someone-else", "u", "u12345", "user", "x", "xy_fake",
})

_SALT = "delfin-personal-data-guard-v1:"

#: sha256(_SALT + lowercase word)[:24]. Matched against every word.
_NAME_HASHES = frozenset({
    "752c064907be1ba4af099fb0",
    "3566d595b83c569789798b5a",
    "47d34e0840ad03bef016336b",
    "35663ac3f1833477a6ea6a75",
    "8502466158414eef14dffcbc",
    "bd0cf1ee584bdcd6aa6ea7a9",
    "7bb2ff9daa151d6b87959d5a",
})

#: Matched only in the possessive (``<name>'s``): the word is also an
#: ordinary token elsewhere.
_POSSESSIVE_HASHES = frozenset({
    "9bf62ba94cbebda278dd5d17",
})

_ATTRIBUTION_FILES = frozenset({"CITATION.cff", "README.md"})

_WORD_RE = re.compile(r"[A-Za-z]+")
_POSSESSIVE_RE = re.compile(r"\b([A-Za-z]+)'s\b")


def _h(word: str) -> str:
    return hashlib.sha256((_SALT + word.lower()).encode()).hexdigest()[:24]


def _findings(text: str, *, names: bool = True,
              name_hashes: frozenset = _NAME_HASHES,
              possessive_hashes: frozenset = _POSSESSIVE_HASHES,
              _cache: dict = {}) -> list[tuple[str, str]]:
    """(rule, matched text) for every rule that fires on ``text``. A name
    is reported by its hash, so a failure message does not repeat it."""
    out: list[tuple[str, str]] = []
    for rule, (literal, rx) in _PATTERNS.items():
        if literal in text:
            out += [(rule, m.group(0)) for m in rx.finditer(text)]
    for m in (_HOME_RE.finditer(text) if "/home/" in text else ()):
        name = m.group(1).rstrip(".-")
        if name not in _HOME_STAND_INS:
            out.append(("home", m.group(0)))
    if names:
        for key in {w.lower() for w in _WORD_RE.findall(text)}:
            # The word, and the word less a trailing "s" (a plural or a
            # possessive written without the apostrophe).
            for form in (key, key[:-1]) if key.endswith("s") else (key,):
                h = _cache.get(form)
                if h is None:
                    h = _cache[form] = _h(form)
                if h in name_hashes:
                    out.append(("name", h))
        for word in (set(_POSSESSIVE_RE.findall(text))
                     if "'s" in text else ()):
            h = _h(word)
            if h in possessive_hashes:
                out.append(("name", h))
    return out


def _tracked() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=str(_ROOT),
                         capture_output=True, text=True, timeout=120)
    return [line for line in out.stdout.splitlines() if line.strip()]


def _readable(rel: str) -> str:
    path = _ROOT / rel
    if path.suffix.lower() in _BINARY:
        return ""
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def test_git_is_available_or_this_checks_nothing():
    files = _tracked()
    assert files, "git ls-files returned nothing; this guard is blind"
    assert "tests/test_no_personal_data_is_committed.py" in files


def test_no_tracked_file_carries_personal_data():
    offenders: dict[str, list] = {}
    for rel in _tracked():
        text = _readable(rel)
        if not text:
            continue
        found = _findings(text, names=rel not in _ATTRIBUTION_FILES)
        if found:
            offenders[rel] = sorted(set(found))[:3]
    assert not offenders, (
        "personal or site data in a public file: "
        + "; ".join(f"{f} {v}" for f, v in sorted(offenders.items())[:10]))


# -- the guard's own shape -------------------------------------------------
# Every refused sample is assembled from parts, so this file stays clean.

@pytest.mark.parametrize("rule, sample", [
    ("node", "ssh " + "uc3" + "n990" + " && tmux attach"),
    ("node", "NODELIST " + "uc3" + "n024"),
    ("site_ip", "rate limit exceeded for " + "141.52" + ".0.1"),
    ("institution_mail", "zz" + "0000" + "@" + "kit.edu"),
    ("session_scratch", "/tmp/" + "tool-1000/" + "-home-" + "someone-x/"),
    ("home", "/home/" + "rosalind" + "/calc"),
])
def test_each_rule_refuses_a_planted_sample(rule, sample):
    assert rule in {r for r, _ in _findings(sample)}, sample


def test_a_name_hash_refuses_the_name_and_its_variants():
    planted = frozenset({_h("zzplanted")})
    for sample in ("Zzplanted's report", "a zzplanted_style case",
                   "TADFs-Zzplanted_3", "finds_zzplanteds_case"):
        assert _findings(sample, name_hashes=planted), sample


def test_a_possessive_hash_refuses_only_the_possessive():
    planted = frozenset({_h("zzmonth")})
    assert _findings("Zzmonth's runs", possessive_hashes=planted)
    assert not _findings('["Zzmonth", 100]', possessive_hashes=planted)


@pytest.mark.parametrize("sample", [
    "ssh login01 && tmux attach -t delfin",
    "rate limit exceeded for 192.0.2.10",
    "https://ki-toolbox.scc.kit.edu/api/v1",
    "/tmp/session-1000/-srv-work-projects/scratchpad",
    "/home/user/projects and /home/<grp>/archive and ~/home/x",
    "/pfs/data6/home/xy/xy_group/xy_user",
    "a field report (2026-06-13)",
])
def test_a_stand_in_passes(sample):
    assert not _findings(sample), sample


def test_the_allow_lists_stay_small():
    """A list that grows quietly stops being a guard."""
    assert _ATTRIBUTION_FILES == {"CITATION.cff", "README.md"}
    assert len(_HOME_STAND_INS) <= 25
