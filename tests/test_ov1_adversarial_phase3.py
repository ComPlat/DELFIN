"""Adversarial tests, OV1 phase 3 — probe fallback + hash cache
(reviewer nacht-s14).

Attacks check_all(use_probe=True) on 0c7c5d297+33fd55d83+a18a35e37:
- does the cache key really make identical text hit the cache (no
  second probe)?
- does an UNAVAILABLE probe leave the static verdict standing (no new
  error, no swallowed static error)?
- do None/empty/huge inputs behave (no crash, no probe)?
- is clear_probe_cache() real (a re-run probes again)?
The probe itself is monkeypatched — no real ORCA is ever run.
"""
import hashlib

import pytest

import delfin.common.orca_input_checks as oic
from delfin.common.orca_input_checks import (
    check_all, clear_probe_cache, _PROBE_CACHE, _probe_cache_key,
)

GEOM = "* xyz 0 1\nC 0.0 0.0 0.0\nH 1.089 0.0 0.0\n*\n"
# Warning-only: %scan start 0.5 A vs geometry bond 1.90 A ->
# scan_start_far_from_geometry (orca_input_checks.py:354), the input
# the builder gave as the reachable warning shape (verified:
# static=[warning], probe fires).
WARNING_ONLY = (
    "! PBE0 def2-SVP Opt\n"
    "%maxcore 2000\n"
    "%pal nprocs 2 end\n"
    "* xyz 0 1\n"
    "H 0.0 0.0 0.0\n"
    "H 0.0 0.0 1.90\n"
    "*\n"
    "%scan\n"
    "{B 0 1 0.5 2.0 10}S\n"
    "end\n"
)


@pytest.fixture(autouse=True)
def _fresh_cache(monkeypatch):
    clear_probe_cache()
    yield
    clear_probe_cache()


def test_clean_and_errored_inputs_never_probe(monkeypatch):
    calls = []
    monkeypatch.setattr(oic, "probe_input",
                        lambda text, **kw: calls.append(text) or {})
    assert check_all("", use_probe=True) == []
    assert calls == []
    # An input that already carries a static error is never probed:
    bad = check_all("! PBE0 def2-SVP NOFREQ Opt\n" + GEOM, use_probe=True)
    assert any(i.severity == "error" for i in bad)
    assert calls == []


def test_identical_text_probed_once(monkeypatch):
    calls = []
    monkeypatch.setattr(oic, "probe_input",
                        lambda text, **kw: calls.append(text) or
                        {"verdict": "parsed", "messages": []})
    check_all(WARNING_ONLY, use_probe=True)
    check_all(WARNING_ONLY, use_probe=True)
    check_all(WARNING_ONLY, use_probe=True)
    assert calls == [WARNING_ONLY] and len(_PROBE_CACHE) == 1


def test_cache_key_is_the_input_hash():
    assert _probe_cache_key("abc") == hashlib.sha256(b"abc").hexdigest()
    assert _probe_cache_key("a\n") != _probe_cache_key("ab\n")


def test_unavailable_probe_adds_no_error(monkeypatch):
    monkeypatch.setattr(oic, "probe_input",
                        lambda text, **kw: {"verdict": "unavailable",
                                            "messages": []})
    issues = check_all(WARNING_ONLY, use_probe=True)
    assert not [i for i in issues if i.severity == "error"]
    # but the static warning is still there:
    assert any(i.severity == "warning" for i in issues)


def test_rejected_probe_is_an_error(monkeypatch):
    monkeypatch.setattr(oic, "probe_input",
                        lambda text, **kw: {"verdict": "rejected",
                                            "messages": ["ORCA says: bad"]})
    issues = check_all(WARNING_ONLY, use_probe=True)
    errs = [i for i in issues if i.severity == "error"
            and i.code == "orca_parser_rejects"]
    assert errs and "bad" in errs[0].message


def test_huge_input_not_probed_static_decides(monkeypatch):
    """A 1 MB static-clean input never reaches the probe (cost guard)."""
    calls = []
    monkeypatch.setattr(oic, "probe_input",
                        lambda text, **kw: calls.append(len(text)) or {})
    text = "! PBE0 def2-SVP Opt\n" + GEOM + ("%pal nprocs 1 end\n" * 50000)
    assert check_all(text, use_probe=True) == []
    assert calls == []


def test_clear_probe_cache_allows_a_second_probe(monkeypatch):
    calls = []
    monkeypatch.setattr(oic, "probe_input",
                        lambda text, **kw: calls.append(text) or
                        {"verdict": "parsed", "messages": []})
    check_all(WARNING_ONLY, use_probe=True)
    clear_probe_cache()
    check_all(WARNING_ONLY, use_probe=True)
    assert calls == [WARNING_ONLY, WARNING_ONLY]
