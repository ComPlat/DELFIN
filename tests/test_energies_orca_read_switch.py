"""Call-site switch tests: energy readers go through delfin.orca_read.

The switch is proven behaviour-identical by the baseline recorder
(scratch/baseline_full.json, 204 runs, 0 mismatches) and, after fixing
orca_read's multi-job FINAL selection (whole-file read on multi-job
outputs), by the final-job comparison over 199 archived runs with 0
mismatches.
"""
from pathlib import Path
from unittest import mock

from delfin import api, energies
from delfin import orca_read


FIX = Path(__file__).resolve().parent / "fixtures" / "orca611" / \
    "butane_dlpno_ccsdt_pvtz"


def _fake_reading(value):
    return orca_read.Reading(value=value, unit="Eh", source="text",
                             note="test")


def test_parse_orca_output_gibbs_and_zpe_via_orca_read():
    """api.parse_orca_output reads gibbs and zpe through orca_read.read."""
    out = FIX / "butane_dlpno_ccsdt_pvtz.out"
    with mock.patch.object(orca_read, "read",
                           side_effect=lambda wd, bn, q, **kw:
                           _fake_reading({"gibbs": -197.0, "zpe": 0.1}.get(
                               q.value))) as mock_read:
        result = api.parse_orca_output(str(out))
    calls = mock_read.call_args_list
    read_quantities = [c.args[2] for c in calls
                       if c.args and len(c.args) >= 3]
    assert orca_read.Quantity.GIBBS in read_quantities
    assert orca_read.Quantity.ZPE in read_quantities
    assert orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY in \
        read_quantities
    assert result.gibbs_free_energy == -197.0
    assert result.zpe == 0.1


def test_cli_report_gibbs_via_orca_read():
    """cli_report._extract_energy_for_state reads gibbs via orca_read."""
    import delfin.cli_report as cli_report
    out = FIX / "butane_dlpno_ccsdt_pvtz.out"
    energies = {}
    with mock.patch.object(orca_read, "read",
                           return_value=_fake_reading(-197.5)) as mock_read:
        cli_report._extract_energy_for_state("0", str(out), energies)
    mock_read.assert_called_once()
    assert energies["0"] == -197.5


def test_collect_gibbs_from_dir_via_orca_read(tmp_path):
    """energies.collect_gibbs_energies_from_dir reads gibbs via orca_read."""
    (tmp_path / "initial.out").write_text("# SYNTHETIC test fixture\n")
    with mock.patch.object(orca_read, "read",
                           return_value=_fake_reading(-100.0)) as mock_read:
        energies.collect_gibbs_energies_from_dir(tmp_path)
    mock_read.assert_called_once()
    args = mock_read.call_args.args
    assert args[0] == str(tmp_path)
    assert args[1] == "initial"
    assert args[2] == orca_read.Quantity.GIBBS
