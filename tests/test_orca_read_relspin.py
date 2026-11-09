"""Tests for ``delfin.orca_read`` reading relativistic/spin-property results.

``read_zfs_parameters`` reads the zero-field-splitting parameters (axial
``D``, rhombic ``E``, cm^-1) from the ``D_Tensor`` record ORCA writes into
the property-JSON file (manual 6.1, property_file.txt).  ``read_soc_states``
reads the SOC-corrected transition list from the ``.out`` text (the
"SOC CORRECTED ... SPECTRUM" block, manual rocis.txt; the same format the
``%tddft soc True`` run prints).

Both are tested against synthetic ORCA-format text written in the test --
never real ORCA output.  Runs under the plain CI interpreter (no OPI).
"""
import json

import pytest

from delfin.orca_read import read_soc_states, read_zfs_parameters


class TestReadZfsParameters:
    def test_reads_d_e_and_the_ratio(self, tmp_path):
        prop = {
            "D_Tensor": {
                "D": 1.23,
                "E": 0.456,
                "d_raw": [1.0, 2.0, 3.0],
                "method": "DFT",
            }
        }
        (tmp_path / "job.property.json").write_text(json.dumps(prop))
        zfs = read_zfs_parameters(tmp_path, "job")
        assert zfs is not None
        assert zfs.d_cm1 == pytest.approx(1.23)
        assert zfs.e_cm1 == pytest.approx(0.456)
        assert zfs.e_over_d == pytest.approx(0.456 / 1.23)

    def test_finds_the_record_in_a_nested_property_json(self, tmp_path):
        prop = {
            "calculation_info": {"gamess": "n"},
            "results": {"D_Tensor": {"D": 0.5, "E": 0.1}},
        }
        (tmp_path / "nested.property.json").write_text(json.dumps(prop))
        zfs = read_zfs_parameters(tmp_path, "nested")
        assert zfs is not None
        assert zfs.d_cm1 == pytest.approx(0.5)
        assert zfs.e_cm1 == pytest.approx(0.1)
        assert zfs.e_over_d == pytest.approx(0.2)

    def test_missing_file_returns_none(self, tmp_path):
        assert read_zfs_parameters(tmp_path, "absent") is None

    def test_property_json_without_the_record_returns_none(self, tmp_path):
        (tmp_path / "empty.property.json").write_text(json.dumps({}))
        assert read_zfs_parameters(tmp_path, "empty") is None

    def test_zero_d_returns_none_ratio(self, tmp_path):
        prop = {"D_Tensor": {"D": 0.0, "E": 0.0}}
        (tmp_path / "zero.property.json").write_text(json.dumps(prop))
        zfs = read_zfs_parameters(tmp_path, "zero")
        assert zfs is not None
        assert zfs.d_cm1 == 0.0
        assert zfs.e_over_d is None


class TestReadSocStates:
    @staticmethod
    def _out(path, body):
        text = (
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "\n" + body
        )
        path.write_text(text)
        return path

    def test_reads_the_soc_corrected_spectrum(self, tmp_path):
        self._out(
            tmp_path / "soc.out",
            "SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC"
            " DIPOLE MOMENTS\n"
            "State           1:    12345.67 cm**-1         1.530 eV\n"
            "State           2:    23456.78 cm**-1         2.909 eV\n",
        )
        states = read_soc_states(tmp_path, "soc")
        assert len(states) == 2
        assert states[0].states == 1
        assert states[0].e_cm1 == pytest.approx(12345.67)
        assert states[0].e_eV == pytest.approx(1.530)
        assert states[1].e_cm1 == pytest.approx(23456.78)
        assert states[1].e_eV == pytest.approx(2.909)

    def test_empty_file_returns_no_states(self, tmp_path):
        self._out(tmp_path / "none.out", "nothing here\n")
        assert read_soc_states(tmp_path, "none") == []

    def test_missing_file_returns_no_states(self, tmp_path):
        assert read_soc_states(tmp_path, "absent") == []
