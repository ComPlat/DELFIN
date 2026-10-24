"""OR1 phase-3 fix: barrier_band must not pass on a band CITATION.

Red tests (reviewer F-REV-4) against unchanged helper: a prose that
merely cites the literature band ("12-20 kcal/mol", en-dash,
"between 12 and 20") must NOT grade as a computed barrier.
"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def C():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "w20_or1_common",
        Path("delfin/agent/pack/benchmark/accept/w20_or1_common.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_natural_barrier_phrasings_pass(C):
    """F-REV-6: correct, artifact-backed answers must NOT be graded FAIL."""
    assert C.barrier_band(
        "The barrier to isomerisation is 48.6 kcal/mol.", 45, 53)
    assert C.barrier_band(
        "The activation energy for the exchange is 15.0 kcal/mol.", 12, 20)
    assert C.barrier_band(
        "barrier height relative to staggered ethane: 2.9 kcal/mol", 2.3, 3.5)
    # The word window must be big enough that "barrier to <long path>
    # is <num>" still anchors: F-REV-6 restated by s22 (long-path case).
    assert C.barrier_band(
        "The barrier to the isomerisation reaction is roughly "
        "48.6 kcal/mol", 45, 53)


def test_cited_band_endpoints_do_not_pass(C):
    assert not C.barrier_band(
        "Barriers reported between 12 and 20 kcal/mol (literature).", 12, 20)


def test_endash_range_citation_does_not_pass(C):
    assert not C.barrier_band("12\u201320 kcal/mol (lit.)", 12, 20)


def test_ascii_range_citation_does_not_pass(C):
    assert not C.barrier_band("12-20 kcal/mol (lit.)", 12, 20)


def test_computed_barrier_statement_passes(C):
    assert C.barrier_band(
        "The computed barrier is 17.3 kcal/mol at the stated level.", 12, 20)


def test_computed_barrier_outside_band_fails(C):
    assert not C.barrier_band(
        "The computed barrier is 33.0 kcal/mol.", 12, 20)


def test_computed_barrier_requires_freq_artifact(C, tmp_path):
    """The artifact-gated wrapper: no frequency output -> no pass."""
    (tmp_path / "report.md").write_text(
        "The computed barrier is 17.3 kcal/mol.", encoding="utf-8")
    assert not C.barrier_with_artifact(tmp_path, "17.3 kcal/mol text", 12, 20)
    (tmp_path / "calc.out").write_text(
        "***imaginary mode***\n1 -2100.00 0.0000\n", encoding="utf-8")
    assert C.barrier_with_artifact(
        tmp_path, "The computed barrier is 17.3 kcal/mol.", 12, 20)
