"""OV1 phase 3: probe_input fallback with an input-hash cache.

When the static checks cannot decide (only warning-severity findings
remain, e.g. an option inside a catalog-gap block), ORCA's own parser
(decided by delfin.orca_knowledge.probe.probe_input) gives the verdict.
The verdict is cached by input hash so one input is never probed twice.

The tests never run a real ORCA: probe_input is monkeypatched with a
fake that counts calls and returns recorded verdicts.
"""
import hashlib

import pytest

from delfin.common import orca_input_checks as oic


def _input(gap_option=False):
    text = "! PBE0 def2-SVP Opt Freq\n%maxcore 2000\n%pal nprocs 2 end\n"
    if gap_option:
        text += "%scf\n  convergence tight\n  BogusOpt 1\nend\n"
    return text


@pytest.fixture()
def probe_calls(monkeypatch):
    """Replace probe_input with a counting fake; returns (calls, list)."""
    calls = []

    def fake_probe_input(text, look_for=None, timeout_s=None, run_dir=None):
        calls.append(text)
        return {"input": text, "verdict": "rejected",
                "messages": ["ORCA: the input is not valid."],
                "look_for": [], "timed_out": False}

    monkeypatch.setattr(oic, "probe_input", fake_probe_input)
    oic.clear_probe_cache()
    return calls


def test_gap_block_warning_triggers_probe(probe_calls):
    """An undecidable input (gap-block option warning) goes to ORCA."""
    issues = oic.check_all(_input(gap_option=True), use_probe=True)
    assert len(probe_calls) == 1
    assert any(i.code == "orca_parser_rejects" for i in issues), issues


def test_decidable_input_never_probed(probe_calls):
    """A clean, fully decidable input is not sent to ORCA."""
    issues = oic.check_all(_input(), use_probe=True)
    assert probe_calls == []
    assert not any(i.code == "orca_parser_rejects" for i in issues)


def test_same_input_probed_once_cached(probe_calls):
    """The verdict is cached by input hash: a second run does not probe."""
    oic.check_all(_input(gap_option=True), use_probe=True)
    oic.check_all(_input(gap_option=True), use_probe=True)
    oic.check_all(_input(gap_option=True), use_probe=True)
    assert len(probe_calls) == 1, probe_calls


def test_different_input_probed_again(probe_calls):
    """A different input has a different hash and is probed again."""
    oic.check_all(_input(gap_option=True), use_probe=True)
    oic.check_all(_input(gap_option=True) + "%pal nprocs 4 end\n",
                  use_probe=True)
    assert len(probe_calls) == 2, probe_calls


def test_unavailable_verdict_is_no_error(probe_calls, monkeypatch):
    """When ORCA is unavailable the probe says so and no error is added."""
    def unavailable(text, look_for=None, timeout_s=None, run_dir=None):
        return {"input": text, "verdict": "unavailable", "messages": [],
                "look_for": [], "timed_out": False}

    monkeypatch.setattr(oic, "probe_input", unavailable)
    issues = oic.check_all(_input(gap_option=True), use_probe=True)
    assert not any(i.code == "orca_parser_rejects" for i in issues)
    assert any(i.code == "unknown_block_option" for i in issues)


def test_parsed_verdict_no_error(probe_calls, monkeypatch):
    """ORCA parsing the input cleanly means no parser error is added."""
    def parsed(text, look_for=None, timeout_s=None, run_dir=None):
        return {"input": text, "verdict": "parsed", "messages": [],
                "look_for": [], "timed_out": False}

    monkeypatch.setattr(oic, "probe_input", parsed)
    issues = oic.check_all(_input(gap_option=True), use_probe=True)
    assert not any(i.code == "orca_parser_rejects" for i in issues)


def test_cache_key_is_input_hash():
    """The cache is keyed by a hash of the input text itself."""
    key = oic._probe_cache_key(_input(gap_option=True))
    assert key == hashlib.sha256(_input(gap_option=True).encode()
                                 ).hexdigest()


def test_public_path_does_not_probe_by_default(probe_calls):
    """QS ruling: the probe stays opt-in on the real public path.

    validate_orca_input calls check_all WITHOUT the probe (comment at
    delfin/api.py documents the why: a default-on probe would spawn a
    real ORCA process on every submit, an expensive step needs a
    stated reason). Nothing may silently depend on the probe through
    the default path.
    """
    from delfin.api import validate_orca_input
    validate_orca_input(_input(gap_option=True))
    assert probe_calls == []
