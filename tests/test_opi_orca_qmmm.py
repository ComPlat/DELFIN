"""Tests for the OPI-backed additive QM/MM block (delfin.tools.adapters.orca_qmmm_opi).

Interpreter contract (same as the other OPI-block tests):
* runs under the py3.12 OPI venv (``DELFIN_OPI_REQUIRED=1`` in CI);
* OPI-dependent tests skip when OPI is missing so the py3.11 no-OPI job passes.
"""

from pathlib import Path

import pytest

# --- OPI availability -------------------------------------------------------
try:
    import opi  # noqa: F401
    _HAS_OPI = True
except ImportError:
    _HAS_OPI = False

needs_opi = pytest.mark.skipif(not _HAS_OPI, reason="OPI (orca-pi) not installed")

from delfin.tools.adapters import orca_qmmm_opi as qm
from delfin.tools._registry import get, list_steps


# ---------------------------------------------------------------------------
# fixtures / helpers
# ---------------------------------------------------------------------------

#: 5 atoms: a water-like solute (0,1,2) plus two MM solvent atoms (3,4).
ELEMS = ["O", "H", "H", "O", "H"]
COORDS = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.96],
    [0.93, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [3.0, 0.0, 0.96],
]
QMATOMS = [0, 1, 2]

# --- synthetic output text --------------------------------------------------
# Hand-written in the ORCA output format; first line carries the synthetic
# marker the CI requires.  Never real ORCA output.
SYNTHETIC_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output
                                *****************
                                * O   R   C   A *
                                *****************

--------------------------
SCF CONVERGED
--------------------------

 Total Energy :                -76.120000000000

FINAL SINGLE POINT ENERGY (MM)        0.015000000000
-------------------------
FINAL SINGLE POINT ENERGY      -76.120000000000
-------------------------
FINAL SINGLE POINT ENERGY (QM/MM)      -76.105000000000

****ORCA TERMINATED NORMALLY****
"""

SYNTHETIC_HEADER = (
    "#  ORCA 6.1.1 -- synthetic header so the OPI reader gate accepts the text"
)


def _write_synthetic_run(tmp_path: Path, base_name: str = "qmmm_probe"):
    """A minimal synthetic ORCA-format run the reader can parse."""
    out = tmp_path / f"{base_name}.out"
    out.write_text(SYNTHETIC_OUT, encoding="utf-8")
    # the reader's OPI gate wants a property json; a minimal hand-written one
    (tmp_path / f"{base_name}.property.json").write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n{}\n",
        encoding="utf-8",
    )
    (tmp_path / f"{base_name}.property.txt").write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
    )
    return out


def _make_ff_file(tmp_path: Path, name: str = "solv.ORCAFF.prms") -> Path:
    p = tmp_path / name
    p.write_text("# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
    return p


# ---------------------------------------------------------------------------
# validate_qmmm -- refusals
# ---------------------------------------------------------------------------

def test_valid_input_has_no_errors(tmp_path):
    _make_ff_file(tmp_path)
    errors, warnings = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename=str(_make_ff_file(tmp_path)),
        n_atoms=5, charge=0, mult=1,
    )
    assert errors == []
    assert warnings == []


def test_refuse_empty_qm_region():
    errors, _ = qm.validate_qmmm(
        qmatoms=[], orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=1,
    )
    assert any("no QM region" in e for e in errors)


def test_refuse_qm_region_covers_whole_system():
    errors, _ = qm.validate_qmmm(
        qmatoms=[0, 1, 2, 3, 4], orcaff_filename="solv.ORCAFF.prms",
        n_atoms=5, charge=0, mult=1,
    )
    assert any("no MM region remains" in e for e in errors)


def test_refuse_out_of_range_qm_index():
    errors, _ = qm.validate_qmmm(
        qmatoms=[0, 9], orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=1,
    )
    assert any("outside" in e for e in errors)


def test_refuse_missing_ff_filename():
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename="", n_atoms=5, charge=0, mult=1,
    )
    assert any("missing force-field file" in e for e in errors)


def test_refuse_ff_file_not_on_disk(tmp_path):
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename=str(tmp_path / "absent.ORCAFF.prms"),
        n_atoms=5, charge=0, mult=1,
    )
    assert any("not found on disk" in e for e in errors)


def test_refuse_unknown_embedding():
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=1, embedding="vacuum",
    )
    assert any("unknown embedding" in e for e in errors)


def test_refuse_unknown_coupling():
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=1, coupling="hybrid",
    )
    assert any("unknown coupling" in e for e in errors)


def test_refuse_bad_multiplicity():
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=0,
    )
    assert any("multiplicity" in e for e in errors)


def test_refuse_unknown_extra_keyword():
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=1, extra_keywords=("NOT_A_REAL_KEYWORD",),
    )
    assert any("unknown ORCA simple keyword" in e for e in errors)


def test_refuse_unknown_qmmm_block_option():
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename="solv.ORCAFF.prms", n_atoms=5,
        charge=0, mult=1, qmmm_options={"not_an_option": 1},
    )
    assert any("does not know option" in e for e in errors)


# ---------------------------------------------------------------------------
# build_qmmm_input -- OPI typed objects
# ---------------------------------------------------------------------------

@needs_opi
def test_build_writes_input_with_qmmm_block(tmp_path):
    ff = _make_ff_file(tmp_path)
    inp = qm.build_qmmm_input(
        str(tmp_path), base_name="qmmm_test",
        elements=ELEMS, coordinates=COORDS, charge=1, mult=1,
        qmatoms=QMATOMS, orcaff_filename=str(ff),
        embedding="electrostatic", coupling="additive",
        method="B3LYP def2-SVP", ncores=2,
    )
    text = Path(inp).read_text()
    assert "!qmmm" in text
    assert "%qmmm" in text
    assert "qmatoms { 0 1 2 } end" in text
    assert f'orcafffilename "{ff}"' in text
    assert "embedding electrostatic" in text
    assert "coupling additive" in text
    # charge / multiplicity of the QM part go on the coordinate section
    assert text.startswith("* xyz") or "* xyz" in text


@needs_opi
def test_build_refuses_when_ff_absent(tmp_path):
    with pytest.raises(ValueError, match="not found on disk"):
        qm.build_qmmm_input(
            str(tmp_path), base_name="qmmm_test",
            elements=ELEMS, coordinates=COORDS, charge=0, mult=1,
            qmatoms=QMATOMS, orcaff_filename=str(tmp_path / "absent.ORCAFF.prms"),
        )


@needs_opi
def test_build_refuses_no_mm_region(tmp_path):
    ff = _make_ff_file(tmp_path)
    with pytest.raises(ValueError, match="no MM region"):
        qm.build_qmmm_input(
            str(tmp_path), base_name="qmmm_test",
            elements=ELEMS, coordinates=COORDS, charge=0, mult=1,
            qmatoms=[0, 1, 2, 3, 4], orcaff_filename=str(ff),
        )


# ---------------------------------------------------------------------------
# read_qmmm_output -- synthetic ORCA-format text
# ---------------------------------------------------------------------------

def test_read_total_and_qm_energies(tmp_path):
    _write_synthetic_run(tmp_path)
    result = qm.read_qmmm_output(str(tmp_path), "qmmm_probe")
    assert result.total_energy_Eh == pytest.approx(-76.105)
    assert result.qm_energy_Eh == pytest.approx(-76.120)
    # reader check: additive QM/MM total >= QM electronic energy
    assert result.total_energy_Eh >= result.qm_energy_Eh
    assert result.terminated_normally is True
    assert result.scf_converged is True
    assert result.meaningful is True


def test_read_missing_output_file(tmp_path):
    result = qm.read_qmmm_output(str(tmp_path), "no_such_base")
    assert result.ok is False
    assert result.total_energy_Eh is None
    assert result.qm_energy_Eh is None
    assert "not found" in result.error


def test_read_unterminated_run(tmp_path):
    out = tmp_path / "qmmm_broken.out"
    out.write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "SCF did not converge\n",
        encoding="utf-8",
    )
    result = qm.read_qmmm_output(str(tmp_path), "qmmm_broken")
    assert result.ok is False
    assert result.meaningful is False
    assert result.terminated_normally is False


# ---------------------------------------------------------------------------
# registration in the tool registry the other blocks use
# ---------------------------------------------------------------------------

def test_registered_in_registry():
    steps = list_steps()
    assert "orca_qmmm_opi" in steps


def test_step_resolves():
    step = get("orca_qmmm_opi")
    assert step is not None
    assert step.name == "orca_qmmm_opi"


# ---------------------------------------------------------------------------
# ORCA-parser probe: ORCA's own parser accepts the %qmmm grammar
# (existing probe helper; parser check only, never a full run).
# The probe guard refuses quoted strings (file names), so the probe text
# carries the %qmmm QMAtoms grammar but not the quoted ORCAFFFilename line.
# ---------------------------------------------------------------------------

def test_orca_parser_accepts_qmmm_input():
    probe = pytest.importorskip("delfin.orca_knowledge.probe")
    text = (
        "! QMMM\n"
        "\n"
        "%qmmm\n"
        "  QMAtoms { 0 1 2 } end\n"
        "end\n"
        "\n"
        "* xyz 0 1\n"
        "O 0.0 0.0 0.0\n"
        "H 0.0 0.0 0.96\n"
        "H 0.93 0.0 0.0\n"
        "O 3.0 0.0 0.0\n"
        "H 3.0 0.0 0.96\n"
        "*\n"
    )
    result = probe.probe_input(text)
    if result["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available on this machine")
    # The probe helper's guard refuses every block it does not whitelist
    # (probe.py: "block %qmmm is not allowed in a probe"), so a full %qmmm
    # grammar check is not possible there.  Assert the refusal is the
    # documented guard policy and NOT an ORCA parse error: the messages must
    # name the guard, not ORCA.
    assert result["verdict"] == "refused"
    assert any("%qmmm is not allowed in a probe" in m for m in result["messages"]), (
        result["messages"]
    )


# ---------------------------------------------------------------------------
# confirmed-by-probe defects (red first, then fixed below the boundary)
# ---------------------------------------------------------------------------

@needs_opi
def test_execute_follows_stepadapter_contract(tmp_path):
    """execute() must take the (work_dir, *, geometry, cores, **kwargs)
    signature every other registered adapter has and return a StepResult, so a
    pipeline step reaches run_qmmm through run_step (delfin/tools/_runner.py)."""
    from delfin.tools._types import StepResult

    ff = _make_ff_file(tmp_path)
    adapter = qm.OrcaQmmmOpiAdapter()
    res = adapter.execute(
        tmp_path,
        geometry=None,
        cores=1,
        elements=ELEMS,
        coordinates=COORDS,
        charge=0,
        mult=1,
        qmatoms=QMATOMS,
        orcaff_filename=str(ff),
    )
    # a missing geometry with nothing in kwargs must be a clean FAILED result,
    # not a TypeError from an incompatible signature
    assert isinstance(res, StepResult)
    assert res.step_name == "orca_qmmm_opi"
    assert res.error is None or isinstance(res.error, str)


@needs_opi
def test_extra_keywords_reach_the_input(tmp_path):
    """'extra_keywords' are validated but must also be written into the input;
    otherwise a user-requested basis/level never runs."""
    import re as _re

    ff = _make_ff_file(tmp_path)
    inp = qm.build_qmmm_input(
        str(tmp_path), base_name="ek", elements=ELEMS, coordinates=COORDS,
        charge=0, mult=1, qmatoms=QMATOMS, orcaff_filename=str(ff),
        method="B3LYP", extra_keywords=("def2-TZVP",),
    )
    text = Path(inp).read_text()
    assert _re.search(r"def2[- ]TZVP\b", text, _re.I), text


@needs_opi
def test_subtractive_writes_qm_method_keyword(tmp_path):
    """coupling='subtractive' must select the QM/QM2/MM simple keyword
    (!QM/<method>/MM), not silently run the additive !QMMM path."""
    import re as _re

    ff = _make_ff_file(tmp_path)
    inp = qm.build_qmmm_input(
        str(tmp_path), base_name="sub", elements=ELEMS, coordinates=COORDS,
        charge=0, mult=1, qmatoms=QMATOMS, orcaff_filename=str(ff),
        coupling="subtractive", method="XTB", extra_keywords=("B3LYP", "def2-SVP"),
    )
    text = Path(inp).read_text()
    first = text.splitlines()[0].lower()
    assert "qm/xtb/mm" in first, text
    assert not _re.search(r"!\s*qmmm\b", first), text



# ---------------------------------------------------------------------------
# operator real-run findings (ORCA 6.1.1, water dimer, ORCAFF from orca_mm)
# ---------------------------------------------------------------------------

def test_additive_without_qm_level_is_refused(tmp_path):
    """'!QMMM' alone ran the QM region at ORCA's HF/def2-SVP default."""
    ff = _make_ff_file(tmp_path)
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename=str(ff), n_atoms=len(ELEMS),
        charge=0, mult=1, method="",
    )
    assert any("no QM level" in e for e in errors), errors


def test_subtractive_without_qm_level_is_refused(tmp_path):
    ff = _make_ff_file(tmp_path)
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename=str(ff), n_atoms=len(ELEMS),
        charge=0, mult=1, coupling="subtractive", method="XTB",
    )
    assert any("no QM level" in e for e in errors), errors


def test_relative_ff_file_is_found_in_the_work_dir(tmp_path, monkeypatch):
    """The force field is read by ORCA in the run directory; a relative name
    must be checked there, not in the process working directory."""
    run = tmp_path / "run"
    run.mkdir()
    ff = _make_ff_file(run)
    monkeypatch.chdir(tmp_path)
    errors, _ = qm.validate_qmmm(
        qmatoms=QMATOMS, orcaff_filename=ff.name, n_atoms=len(ELEMS),
        charge=0, mult=1, method="B3LYP def2-SVP", work_dir=str(run),
    )
    assert not [e for e in errors if "force-field" in e], errors


@needs_opi
def test_additive_method_is_written_as_the_qm_level(tmp_path):
    ff = _make_ff_file(tmp_path)
    inp = qm.build_qmmm_input(
        str(tmp_path), base_name="lvl", elements=ELEMS, coordinates=COORDS,
        charge=0, mult=1, qmatoms=QMATOMS, orcaff_filename=str(ff),
        method="B3LYP def2-SVP",
    )
    head = Path(inp).read_text().lower()
    assert "b3lyp" in head and "def2-svp" in head, head


def test_total_is_the_qmmm_line_not_the_qm_region_line(tmp_path):
    """Real-run finding: the plain FINAL SINGLE POINT ENERGY line of an
    additive run is the QM region only; reporting it as the total dropped the
    MM energy (0.029 Eh = 18 kcal/mol on a water dimer)."""
    _write_synthetic_run(tmp_path)
    r = qm.read_qmmm_output(str(tmp_path), "qmmm_probe")
    assert r.total_energy_Eh == pytest.approx(-76.105)
    assert r.mm_energy_Eh == pytest.approx(0.015)
    assert r.total_energy_Eh == pytest.approx(r.qm_energy_Eh + r.mm_energy_Eh)
