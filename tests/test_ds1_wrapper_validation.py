"""DS1 R1/R2: wrappers refuse what they cannot run instead of guessing.

A Multiwfn method name that is not on the menu must not silently run the
first menu entry while the result still names the requested method, and an
atom index outside the geometry must be refused with its value named.
"""
import pytest

from delfin.analysis_tools import morfeus_wrapper, multiwfn_wrapper


def test_unknown_multiwfn_method_is_refused_not_mapped_to_entry_one(tmp_path):
    wfn = tmp_path / "x.molden"
    wfn.write_text("")
    with pytest.raises(ValueError, match="mulikken"):
        multiwfn_wrapper.population_analysis(wfn, method="mulikken")
    with pytest.raises(ValueError, match="must be a string"):
        multiwfn_wrapper.bond_order_analysis(wfn, method=None)


def _water(tmp_path):
    xyz = tmp_path / "w.xyz"
    xyz.write_text("3\n\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    return xyz


def test_atom_index_outside_the_geometry_is_refused_with_its_value(tmp_path):
    pytest.importorskip("morfeus")
    xyz = _water(tmp_path)
    with pytest.raises(ValueError, match="'metal_index'=3 is out of range"):
        morfeus_wrapper.cone_angle(xyz, 3)
    with pytest.raises(ValueError, match="'atom1_index'=-1 is negative"):
        morfeus_wrapper.sterimol(xyz, -1, 1)


def test_steric_report_records_the_validation_error_once(tmp_path):
    pytest.importorskip("morfeus")
    report = morfeus_wrapper.full_steric_report(_water(tmp_path), 7)
    assert "out of range" in report["validation_error"]
    assert "buried_volume" not in report and "cone_angle" not in report
