# SPDX-License-Identifier: MIT
"""Every credential name the product reads is one the scrubber knows.

``credentials._WELL_KNOWN_KEYS`` was written for ``list_credentials``, to
show env-var-sourced keys when nothing is on disk. It later gained a second
consumer, ``web_tools._known_secret_values``, which is the exact-value
source for two security paths: ``output_guard.scrub_secrets`` on every tool
result before it enters the transcript, and ``outbound_secret_reason`` on a
URL or search query before it leaves the machine.

A credential in the ``~/.delfin`` store is covered whatever its name --
``held_secret_values`` enumerates the store. One supplied through the
ENVIRONMENT is covered only if its name is in that list.

Measured 2026-11-21 before the list was extended: with the value in the
environment, three names were scrubbed and refused egress and four were
neither, and all four are read by the product and documented as "env var or
~/.delfin credential file", i.e. offered as an equivalent route while only
one of them was protected.

Input: the credential-name call sites in ``delfin/``, read from the source.
Output: the names, compared against the list. Deriving them from the call
sites rather than restating them is the point -- a restated list is the
thing that drifted.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

from delfin.agent.credentials import _WELL_KNOWN_KEYS

_ROOT = Path(__file__).resolve().parents[1]

#: How a credential is read: the shared loader, and the modules' own
#: one-line wrappers around it.
_CALL_RE = re.compile(
    r"(?:load_credential|_credential|_load_cred|_load_cred_anth)"
    r"\(\s*[\"']([A-Z][A-Z0-9_]+)[\"']")

#: A provider table: a lowercase provider name mapping to a secret-shaped
#: environment name is a credential this product reads for that provider.
_TABLE_RE = re.compile(
    r"[\"'][a-z][a-z0-9_]*[\"']\s*:\s*"
    r"[\"']([A-Z][A-Z0-9_]*(?:_API_KEY|_TOKEN))[\"']")

#: WHAT THIS SCAN DOES NOT SEE, named rather than implied: a name assembled
#: at run time (``f"{vendor}_API_KEY"``), one read straight from
#: ``os.environ`` without going through a credential helper, and one that
#: reaches the process from a settings file. Those routes exist, so this
#: file is a ratchet on the two routes it can read and not a proof that the
#: list is complete. Widening it to every secret-shaped string literal in
#: ``delfin/`` was considered and not done here: it needs a judgement per
#: name (``JUPYTER_TOKEN`` and ``DELFIN_VOILA_TOKEN`` are credentials worth
#: scrubbing, ``HF_TOKEN_PATH`` is a path and ``DEFAULT_RUN_BUDGET_TOKENS``
#: is a number), and an exemption list nobody has reasoned through is worse
#: than a narrower claim that holds.

#: A name is a SECRET when it is shaped like one. An identifier that merely
#: travels beside a secret is not: GOOGLE_CSE_CX names a search engine, and
#: scrubbing it would redact an ordinary id out of every tool result.
_SECRET_SUFFIXES = ("_API_KEY", "_TOKEN", "_SECRET", "_PASSWORD", "_KEY")
_NOT_SECRET = frozenset({"GOOGLE_CSE_CX"})


def _read_names() -> dict[str, list[str]]:
    """name -> the files that read it as a credential."""
    found: dict[str, list[str]] = {}
    for path in sorted((_ROOT / "delfin").rglob("*.py")):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for name in set(_CALL_RE.findall(text)) | set(
                _TABLE_RE.findall(text)):
            found.setdefault(name, []).append(
                str(path.relative_to(_ROOT)))
    return found


def _secret_names() -> dict[str, list[str]]:
    return {n: f for n, f in _read_names().items()
            if n not in _NOT_SECRET
            and any(n.endswith(s) for s in _SECRET_SUFFIXES)}


def test_the_scan_finds_the_call_sites_at_all():
    """A scan that found nothing would pass every assertion below it."""
    names = _secret_names()
    # Eight names, measured 2026-11-21 over the two routes this scan
    # reads, asserted as a floor so a scan that stops finding call sites
    # fails here instead of passing every assertion below it. GITHUB_TOKEN
    # and GH_TOKEN are deliberately NOT expected: they are on the scrub
    # list and the product reads them by neither route, which is the
    # limitation named above rather than a hole -- they are covered
    # because the list names them.
    assert len(names) >= 8, sorted(names)
    for expected in ("KIT_TOOLBOX_API_KEY", "BRAVE_SEARCH_API_KEY",
                     "OPENAI_API_KEY", "SEMANTIC_SCHOLAR_API_KEY"):
        assert expected in names, (expected, sorted(names))


def test_every_credential_name_read_is_on_the_scrub_list():
    names = _secret_names()
    missing = {n: f for n, f in sorted(names.items())
               if n not in set(_WELL_KNOWN_KEYS)}
    assert not missing, (
        "read as a credential but not in _WELL_KNOWN_KEYS, so its value is "
        "neither scrubbed from a tool result nor refused egress when it "
        "comes from the environment: "
        + "; ".join(f"{n} ({', '.join(f)})" for n, f in missing.items()))


@pytest.mark.parametrize("name", sorted(_WELL_KNOWN_KEYS))
def test_a_listed_name_is_scrubbed_and_refused_egress(name, monkeypatch,
                                                      tmp_path):
    """The list is only a claim; this drives the two paths that consume it.

    Both directions of the same value: out of a tool result, and out of an
    outbound URL. A stand-in value, never a real credential, and HOME is
    redirected so no real store is read.
    """
    from delfin.agent.output_guard import scrub_secrets
    from delfin.agent.web_tools import (_known_secret_values,
                                        outbound_secret_reason)

    fake = "zzz-not-a-real-credential-0123456789abcdef"
    monkeypatch.setenv("HOME", str(tmp_path))
    for key in list(os.environ):
        if key.endswith(("_API_KEY", "_TOKEN")):
            monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv(name, fake)

    assert fake in _known_secret_values(), f"{name} is not a held value"
    assert fake not in scrub_secrets(f"the tool printed key={fake}"), (
        f"{name}'s value survived into a tool result")
    assert outbound_secret_reason(f"https://x.invalid/?d={fake}"), (
        f"{name}'s value was allowed to leave in a URL")


def test_an_identifier_beside_a_secret_is_not_scrubbed():
    """Why _NOT_SECRET exists, pinned so it cannot grow quietly.

    Every name it holds must be one the product really reads -- an entry
    for a name nobody reads would be an exemption with nothing behind it.
    """
    read = _read_names()
    for name in _NOT_SECRET:
        assert name in read, (
            f"{name} is exempted from the scrub list but nothing reads it")
        assert name not in set(_WELL_KNOWN_KEYS), (
            f"{name} is both exempted and listed")
