"""Tests for ``delfin.orca_read`` reading relativistic/spin-property results.

``read_zfs_parameters`` reads the zero-field-splitting parameters (axial
``D``, rhombic ``E``, cm^-1) from the ``D_Tensor`` record ORCA writes into
the property-JSON file (manual 6.1, property_file.txt).  ``read_soc_states``
reads the SOC-corrected transition list from the ``.out`` text (the
"SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
table a ``%tddft soc true`` run prints; the layout below follows a real
ORCA 6.1.1 run, the numbers are invented).

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

    def test_reads_the_orca_6_1_1_layout(self, tmp_path):
        # Real ORCA 6.1.1 nesting: Geometries[i].D_Tensor is a list of
        # records (numbers invented).
        prop = {
            "Calculation_Status": {"Status": "NORMAL TERMINATION"},
            "Geometries": [{
                "DFT_Energy": [{"Method": "DFT", "Mult": 3}],
                "D_Tensor": [{
                    "D": 2.5, "E": 0.25, "Irrep": 0, "Level": "Relaxed density",
                    "Method": "SCF", "Mult": 3, "State": -1,
                    "d_eigenvalues": [[-1.0], [-0.5], [1.5]],
                    "d_raw": [[-1.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, 1.5]],
                }],
            }],
        }
        (tmp_path / "o2.property.json").write_text(json.dumps(prop))
        zfs = read_zfs_parameters(tmp_path, "o2")
        assert zfs is not None
        assert zfs.d_cm1 == pytest.approx(2.5)
        assert zfs.e_cm1 == pytest.approx(0.25)
        assert zfs.e_over_d == pytest.approx(0.1)

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


# Synthetic SOC table in the layout of a real ORCA 6.1.1 ``%tddft soc true``
# run; every number is invented.
SOC_TABLE = (
    "-" * 104 + "\n"
    "      SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
    + "-" * 104 + "\n"
    "      Transition         Energy     Energy  Wavelength fosc(D2)      D2       |DX|      |DY|      |DZ|\n"
    "                          (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)\n"
    + "-" * 104 + "\n"
    "  0-1.0A  ->  1-3.0A    2.500000   20164.0   495.9   0.000000010   0.00000   0.00000   0.00010   0.00000\n"
    "  0-1.0A  ->  2-1.0A    4.160000   33552.9   298.0   0.012345678   0.10000   0.00000   0.31623   0.00000\n"
    "  0-1.0A  -> 10-3.0A    6.250000   50409.6   198.4   0.000000000   0.00000   0.00000   0.00000   0.00000\n"
    "\n"
    + "-" * 104 + "\n"
)


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
        # ORCA 6.1.1 layout: <i>-<mult><irrep> -> <f>-<mult><irrep>, then
        # eV, cm-1, nm, fosc(D2), D2, |DX|, |DY|, |DZ|.
        self._out(tmp_path / "soc.out", SOC_TABLE)
        states = read_soc_states(tmp_path, "soc")
        assert len(states) == 3
        assert states[0].states == 1
        assert states[0].e_eV == pytest.approx(2.500000)
        assert states[0].e_cm1 == pytest.approx(20164.0)
        assert states[2].states == 10
        assert states[2].e_eV == pytest.approx(6.250000)
        assert states[2].e_cm1 == pytest.approx(50409.6)

    def test_reads_fosc_wavelength_and_spin_label(self, tmp_path):
        self._out(tmp_path / "soc.out", SOC_TABLE)
        states = read_soc_states(tmp_path, "soc")
        assert states[0].mult == pytest.approx(3.0)
        assert states[1].mult == pytest.approx(1.0)
        assert states[1].fosc == pytest.approx(0.012345678)
        assert states[1].wavelength_nm == pytest.approx(298.0)
        assert states[0].initial_state == 0

    def test_reads_only_the_electric_dipole_table(self, tmp_path):
        # The velocity-gauge table that follows must not be appended, and the
        # spin-free spectrum before it must not be read at all.
        body = (
            "     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
            "  0-1A  ->  1-3A    9.900000   79848.0   125.2   0.0   0.0   0.0   0.0   0.0\n\n"
            + SOC_TABLE
            + "      SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\n"
            "  0-1.0A  ->  1-3.0A    2.500000   20164.0   495.9   0.0   0.0   0.0   0.0   0.0\n"
        )
        self._out(tmp_path / "soc.out", body)
        states = read_soc_states(tmp_path, "soc")
        assert [s.states for s in states] == [1, 2, 10]

    def test_empty_file_returns_no_states(self, tmp_path):
        self._out(tmp_path / "none.out", "nothing here\n")
        assert read_soc_states(tmp_path, "none") == []

    def test_missing_file_returns_no_states(self, tmp_path):
        assert read_soc_states(tmp_path, "absent") == []
