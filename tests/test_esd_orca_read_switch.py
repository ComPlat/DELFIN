"""Call-site switch tests: the ESD readers go through delfin.orca_read.

The switch is proven behaviour-identical by the baseline recorder
(scratch/baseline.json, 209 runs: the four tests/fixtures/orca611
fixtures plus archived ESD runs, 0 mismatches).
"""
from pathlib import Path
from unittest import mock

from delfin import esd_input_generator, esd_results, orca_read


FIX = Path(__file__).resolve().parent / "fixtures" / "orca611" / \
    "butane_dlpno_ccsdt_pvtz"


def _fake_reading(value):
    return orca_read.Reading(value=value, unit="Eh", source="text",
                             note="test")


def test_collect_esd_results_energies_via_orca_read(tmp_path):
    """collect_esd_results reads fspe/zpe/gibbs through orca_read.read."""
    out = tmp_path / "S0.out"
    out.write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
    with mock.patch.object(orca_read, "read",
                           return_value=_fake_reading(-100.0)) as mock_read:
        esd_results.collect_esd_results(tmp_path, ["S0"], [], [])
    read_quantities = [c.args[2] for c in mock_read.call_args_list
                       if c.args and len(c.args) >= 3]
    assert orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY in read_quantities
    assert orca_read.Quantity.ZPE in read_quantities
    assert orca_read.Quantity.GIBBS in read_quantities
    # every call addresses the state output in the ESD dir
    for call in mock_read.call_args_list:
        assert call.args[0] == str(tmp_path)
        assert call.args[1] == "S0"


def test_collect_esd_results_falls_back_when_orca_read_none(tmp_path):
    """A None reading (failed/old run) falls back to the text find_*."""
    out = tmp_path / "S0.out"
    out.write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "FINAL SINGLE POINT ENERGY      -76.123456789012 Eh\n"
        "Zero point energy             0.012345678901 Eh\n"
        "Final Gibbs free energy       -76.098765432109 Eh\n")
    with mock.patch.object(orca_read, "read",
                           return_value=orca_read.Reading(
                               value=None, unit="", source="none",
                               note="no value")) as mock_read:
        summary = esd_results.collect_esd_results(tmp_path, ["S0"], [], [])
    assert mock_read.called
    assert summary.states["S0"].fspe == -76.123456789012
    assert summary.states["S0"].zpe == 0.012345678901
    assert summary.states["S0"].gibbs == -76.098765432109


def test_calculate_dele_cm1_via_orca_read(tmp_path):
    """calculate_dele_cm1 reads the two state energies via orca_read."""
    (tmp_path / "S0.out").write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
    (tmp_path / "S1.out").write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
    values = {"S0": -76.0, "S1": -75.7}
    with mock.patch.object(orca_read, "read",
                           side_effect=lambda wd, bn, q, **kw:
                           _fake_reading(values.get(bn.upper()))) as mock_read:
        dele = esd_input_generator.calculate_dele_cm1(
            str(tmp_path / "S0.out"), str(tmp_path / "S1.out"))
    assert mock_read.called
    read_quantities = [c.args[2] for c in mock_read.call_args_list
                       if c.args and len(c.args) >= 3]
    assert orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY in read_quantities
    # DELE = E(S0) - E(S1) = -0.3 Eh = -65842.3889 cm^-1
    assert abs(dele - (-0.3 * 219474.63)) < 0.1


def test_reader_text_esd_sites_equal_old_parsers(tmp_path, monkeypatch):
    """With orca.reader=text the switched ESD sites read the old way."""
    monkeypatch.setattr(orca_read, "_reader_mode", lambda: "text")
    from delfin.energies import find_electronic_energy, find_ZPE, \
        find_gibbs_energy

    out_src = tmp_path / "S0.out"
    out_src.write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "FINAL SINGLE POINT ENERGY      -76.123456789012 Eh\n"
        "Zero point energy             0.012345678901 Eh\n"
        "Final Gibbs free energy       -76.098765432109 Eh\n")
    summary = esd_results.collect_esd_results(tmp_path, ["S0"], [], [])
    st = summary.states["S0"]
    assert st.fspe == find_electronic_energy(str(out_src))
    assert st.zpe == find_ZPE(str(out_src))
    assert st.gibbs == find_gibbs_energy(str(out_src))
    dele = esd_input_generator.calculate_dele_cm1(
        str(out_src), str(out_src))
    assert dele == 0.0
