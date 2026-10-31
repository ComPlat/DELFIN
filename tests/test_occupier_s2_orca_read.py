"""OCCUPIER <S**2> deviation reading now goes through delfin.orca_read.

The occupier readers (delfin/occupier.py `_check_contamination_trigger` and
delfin/reporting/occupier_reports.py `get_dev`) read the <S**2> deviation of a
finished OCCUPIER run. They were switched from calling
`delfin.parser.extract_last_uhf_deviation` directly to
`delfin.orca_read.read(..., Quantity.S2_DEVIATION)`.

This must not change behaviour. The occupier flow runs ORCA into `output*.out`
and never places an OPI `*.property.json` beside it, so orca_read resolves to
DELFIN's established text reader (extract_last_uhf_deviation) and returns the
identical value. These tests pin that identity:
  * `orca.reader = "text"` restores the old path exactly (requirement);
  * in the default (auto) mode, absent a property.json, the value is the same
    text read (open-shell float and closed-shell None).

Each test copies the fixture .out/.inp into tmp_path WITHOUT the OPI
property.json, matching how the occupier flow presents runs to the readers.
"""
from pathlib import Path

from delfin.occupier import _orca_s2_deviation as occ_dev
from delfin.parser import extract_last_uhf_deviation
from delfin.reporting.occupier_reports import _orca_s2_deviation as rep_dev

FIXTURES = Path(__file__).parent / "fixtures" / "orca611"


def _copy_without_property_json(name: str, tmp_path: Path) -> Path:
    """Copy a fixture .out/.inp (minus the OPI JSON) into tmp_path."""
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    for suffix in (".out", ".inp"):
        src = FIXTURES / name / f"{name}{suffix}"
        if src.is_file():
            (dst / src.name).write_bytes(src.read_bytes())
    return dst


def test_text_mode_restores_the_old_path(tmp_path, monkeypatch):
    """orca.reader=text: `_orca_s2_deviation` == extract_last_uhf_deviation."""
    from delfin import orca_read
    monkeypatch.setattr(orca_read, "_reader_mode", lambda: "text")
    for name in ("h2o_optfreq", "h2o_tddft", "o2_triplet", "scf_not_converged"):
        d = _copy_without_property_json(name, tmp_path)
        out = d / f"{name}.out"
        expected = extract_last_uhf_deviation(str(out))
        assert occ_dev(out) == expected, f"occupier._orca_s2_deviation {name}"
        assert rep_dev(str(out)) == expected, f"occupier_reports._orca_s2_deviation {name}"


def test_default_mode_without_property_json_is_identical(tmp_path):
    """Auto mode, no property.json (occupier flow): value unchanged."""
    for name in ("h2o_optfreq", "h2o_tddft", "o2_triplet", "scf_not_converged"):
        d = _copy_without_property_json(name, tmp_path)
        out = d / f"{name}.out"
        expected = extract_last_uhf_deviation(str(out))
        assert occ_dev(out) == expected, f"occupier._orca_s2_deviation {name}"
        assert rep_dev(str(out)) == expected, f"occupier_reports._orca_s2_deviation {name}"


def test_open_shell_value_preserved(tmp_path):
    """The meaningful open-shell deviation is unchanged (float identity)."""
    name = "o2_triplet"
    d = _copy_without_property_json(name, tmp_path)
    out = d / f"{name}.out"
    expected = extract_last_uhf_deviation(str(out))
    assert expected is not None  # the fixture is open shell
    assert occ_dev(out) == expected
    assert rep_dev(str(out)) == expected


def test_closed_shell_stays_none(tmp_path):
    """Closed-shell runs keep the historical None (no UHF block)."""
    for name in ("h2o_optfreq", "h2o_tddft"):
        d = _copy_without_property_json(name, tmp_path)
        out = d / f"{name}.out"
        assert extract_last_uhf_deviation(str(out)) is None  # sanity
        assert occ_dev(out) is None
        assert rep_dev(str(out)) is None


def test_failed_run_stays_none(tmp_path):
    name = "scf_not_converged"
    d = _copy_without_property_json(name, tmp_path)
    out = d / f"{name}.out"
    assert occ_dev(out) is None
    assert rep_dev(str(out)) is None
