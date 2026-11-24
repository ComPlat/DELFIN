"""Caller-switch tests: the remaining delfin.energies path callers read
through delfin.orca_read.

The switch is behaviour-identical where text is read (baseline recorder over
204 runs, 0 gibbs mismatches; see tests/test_energies_orca_read_switch.py),
because orca_read's shadow check returns the TEXT value on any mismatch.
The old parser stays as the fallback when orca_read has no value, so a run
without OPI answers exactly as before.
"""
from pathlib import Path
from unittest import mock

from delfin import orca_read


FIX = Path(__file__).resolve().parent / "fixtures" / "orca611" / \
    "butane_dlpno_ccsdt_pvtz"


def _fake_reading(value, source="text"):
    return orca_read.Reading(value=value, unit="Eh", source=source,
                             note="test")


def test_result_critic_energy_via_orca_read():
    """agent.result_critic._check_energy reads the final energy via
    orca_read.read."""
    from delfin.agent import result_critic

    out = FIX / "butane_dlpno_ccsdt_pvtz.out"
    with mock.patch.object(orca_read, "read",
                           return_value=_fake_reading(-6.6), ) as mock_read:
        critiques = result_critic._check_energy(out)
    mock_read.assert_called_once()
    args = mock_read.call_args.args
    assert args[2] == orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY
    assert critiques[0].level == "ok"
    assert critiques[0].message.endswith("Eh.")


def test_result_critic_falls_back_when_reading_is_none():
    """The text parser answers when orca_read has no value — a run without
    OPI is unchanged."""
    from delfin.agent import result_critic

    out = FIX / "butane_dlpno_ccsdt_pvtz.out"
    with mock.patch.object(orca_read, "read",
                           return_value=_fake_reading(None)):
        critiques = result_critic._check_energy(out)
    assert critiques, "the critic still answers with the text value"
    assert any("Eh" in c.message for c in critiques)


def test_claim_locator_reads_via_orca_read():
    """agent.claim_locator verifies energy claims through orca_read.read."""
    import delfin.agent.claim_locator as cl

    with mock.patch.object(orca_read, "read",
                          side_effect=lambda wd, bn, q, **kw:
                          _fake_reading({"final_single_point_energy":
                                        -6.6,
                                        "gibbs": -197.5,
                                        "zpe": 0.2}.get(q.value))) as mock_read:
        e = cl._ENERGY_LABELS["FINAL SINGLE POINT ENERGY"](str(FIX /
                                "butane_dlpno_ccsdt_pvtz.out"))
    assert e == -6.6
    assert mock_read.call_count == 1


def test_run_results_numbers_via_orca_read(tmp_path):
    """dashboard.run_results._numbers_from reads E/G/ZPE via orca_read."""
    from unittest import mock as _m
    import delfin.dashboard.run_results as rr

    out = FIX / "butane_dlpno_ccsdt_pvtz.out"
    with _m.patch.object(orca_read, "read",
                         side_effect=lambda wd, bn, q, **kw:
                         _fake_reading({"final_single_point_energy":
                                       -6.6,
                                       "gibbs": -197.5,
                                       "zpe": 0.2}.get(q.value))) as mock_read:
        numbers = rr._numbers_from(out)
    quantities = [c.args[2] for c in mock_read.call_args_list
                  if c.args and len(c.args) >= 3]
    assert orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY in quantities
    assert orca_read.Quantity.GIBBS in quantities
    assert orca_read.Quantity.ZPE in quantities
    assert numbers["energy"] == -6.6
    assert numbers["free_energy"] == -197.5
    assert numbers["zpe"] == 0.2
