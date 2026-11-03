"""EF1 phase-2 reviewer adversarial tests for cluster_efficiency.advise().

Written by nacht-s2 (reviewer); every failure is a finding to nacht-s1.
Recorded/injected inputs only — no live cluster, no machine dependence.
"""
import pytest

from delfin.agent import cluster_efficiency as ce


def _patch(monkeypatch, pal=16, maxcore=3000, parts=("cpu",)):
    monkeypatch.setattr(ce, "_auto_configure",
                        lambda config: {**config, "PAL": pal,
                                        "maxcore": maxcore})
    monkeypatch.setattr(ce, "_configured_partitions",
                        lambda explicit=None, profile=None: parts)


def test_advise_rejects_negative_system_size(monkeypatch):
    _patch(monkeypatch)
    adv = ce.advise("opt", -5)
    # A negative atom count is nonsense: it must not be echoed back as
    # data, it must be flagged.
    assert any("system_size" in n or "negative" in n for n in adv["notes"]), adv


def test_advise_rejects_zero_system_size(monkeypatch):
    _patch(monkeypatch)
    adv = ce.advise("opt", 0)
    assert any("system_size" in n or "zero" in n for n in adv["notes"]), adv


def test_advise_rejects_empty_step(monkeypatch):
    _patch(monkeypatch)
    adv = ce.advise("", 12)
    assert adv["notes"], "empty step silently sized"


def test_advise_step_case_insensitive_or_flagged(monkeypatch):
    _patch(monkeypatch)
    # DELFIN's CONTROL values are lower-case; "OPT" must either be
    # recognised or flagged — never silently treated as known.
    adv = ce.advise("OPT", 12)
    assert any("unknown" in n.lower() for n in adv["notes"]), adv


def test_advise_real_delfin_workflow_names_not_flagged_unknown(monkeypatch):
    # DELFIN's own workflow vocabulary (delfin/tools/templates/
    # __init__.py:39-50, CONTROL method=) is module names like esd,
    # occupier, imag — not the hand-made KNOWN_STEPS list. A real
    # workflow name must not be reported as unknown.
    _patch(monkeypatch)
    for step in ("esd", "occupier", "guppy", "classic"):
        adv = ce.advise(step, 12)
        assert not any("unknown step" in n for n in adv["notes"]), (step, adv)


def test_advise_empty_explicit_partitions_is_not_join_of_nothing(monkeypatch):
    _patch(monkeypatch, parts=())
    adv = ce.advise("opt", 12, explicit_partitions=())
    assert adv["partitions"] is None, adv


def test_advise_partitions_string_never_has_stray_comma(monkeypatch):
    _patch(monkeypatch, parts=("cpu", ""))
    adv = ce.advise("opt", 12, explicit_partitions=("cpu", ""))
    assert adv["partitions"] is None or adv["partitions"] == "cpu", adv


def test_advise_huge_system_size_is_accepted_with_note(monkeypatch):
    _patch(monkeypatch)
    adv = ce.advise("opt", 10**6)
    assert adv["system_size_atoms"] == 10**6
    assert "PAL" in adv
