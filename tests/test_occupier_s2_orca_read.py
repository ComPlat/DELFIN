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

_MARK = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"

# Hand-written in ORCA's text layout (LICENSE_GUARD.md: real ORCA output is
# not committed). The open-shell run prints two contamination blocks, so the
# reader must return the LAST deviation.
_UHF = ("----------------------\n UHF SPIN CONTAMINATION\n"
        "----------------------\n\n"
        "Expectation value of <S**2>     :     {s2:.6f}\n"
        "Ideal value S*(S+1) for S=1.0   :     2.000000\n"
        "Deviation                       :     {dev:.6f}\n\n")
_DONE = "FINAL SINGLE POINT ENERGY      {e:.9f}\n" \
        "                             ****ORCA TERMINATED NORMALLY****\n"
_RUNS = {
    "h2o_optfreq": _DONE.format(e=-76.4),
    "h2o_tddft": _DONE.format(e=-76.3),
    "o2_triplet": (_UHF.format(s2=2.0123, dev=0.0123)
                   + _UHF.format(s2=2.0071, dev=0.0071)
                   + _DONE.format(e=-150.2)),
    "scf_not_converged": "SCF NOT CONVERGED AFTER 125 CYCLES\n"
                         "ORCA finished by error termination in SCF\n",
}


def _copy_without_property_json(name: str, tmp_path: Path) -> Path:
    """Write the synthetic run ``name`` as ``<name>/<name>.out`` (no OPI JSON)."""
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    (dst / f"{name}.out").write_text(_MARK + _RUNS[name])
    return dst


def test_the_synthetic_open_shell_run_reads_its_last_deviation(tmp_path):
    out = _copy_without_property_json("o2_triplet", tmp_path) / "o2_triplet.out"
    assert extract_last_uhf_deviation(str(out)) == 0.0071


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
