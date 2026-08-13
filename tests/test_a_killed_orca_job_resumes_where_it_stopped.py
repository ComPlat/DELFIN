"""A re-submitted job must not repeat work the previous attempt finished.

Across the run archive, 129 TIMEOUT and 24 CANCELLED jobs had computed for
hours and left nothing a restart could pick up. The expensive part is almost
always a geometry optimisation: one archived job spent ~30 h converging an
excited state and then died in the frequency step, and resubmitting it started
again from the original SMILES-derived geometry.

ORCA offers exactly two restart handles, and these tests pin both:

* ``<base>.xyz`` — rewritten after every accepted optimisation cycle.
* ``<base>.hess`` under ``%freq restart true`` — verified against ORCA 6.1.1,
  which then recomputes none of the displacements. ORCA writes that file only
  when the Hessian is complete, so a numerical frequency run killed midway
  genuinely cannot be resumed; that is an ORCA limit, not a policy choice.

Nothing here may change the calculation. Same keywords, same method, same
basis — only the starting coordinates move to where ORCA itself would have
continued, and only in the isolated copy of the input.
"""

from pathlib import Path

import pytest

from delfin.orca import _resume_partial_orca_run


_INLINE = """! wB97X RKS def2-TZVP D4 RIJCOSX def2/J CPCM OPT numFREQ MOREAD
%base "S1"
%moinp "S0.gbw"
%pal nprocs 20 end
%maxcore 5400

%tddft
  nroots 15
  iroot 1
end

* xyz 0 1
  O    0.000000    0.000000    0.117300
  H    0.000000    0.757200   -0.469200
  H    0.000000   -0.757200   -0.469200
*
"""

_RESUMED_XYZ = """3
Coordinates from ORCA-job S1
  O    0.100000    0.200000    0.300000
  H    0.400000    0.500000    0.600000
  H    0.700000    0.800000    0.900000
"""


def _partial_run(tmp_path: Path, *, inp: str = _INLINE, out: str = "died mid-run\n"):
    """A job directory whose previous attempt was killed (no ORCA end marker)."""
    inp_path = tmp_path / "S1.inp"
    out_path = tmp_path / "S1.out"
    inp_path.write_text(inp, encoding="utf-8")
    out_path.write_text(out, encoding="utf-8")
    return inp_path, out_path


def _resume(tmp_path: Path, inp_path: Path, out_path: Path):
    return _resume_partial_orca_run(
        inp_path.read_text(encoding="utf-8"), inp_path, tmp_path, out_path
    )


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def test_the_last_accepted_geometry_replaces_the_starting_one(tmp_path):
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert any("geometry continued from S1.xyz" in n for n in notes), notes
    assert "0.100000    0.200000    0.300000" in content
    assert "0.117300" not in content, "original geometry must be gone"


def test_resuming_changes_nothing_but_the_coordinates(tmp_path):
    """The method line, the blocks and the charge/multiplicity are untouched —
    a resumed job must still be the same calculation."""
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, _ = _resume(tmp_path, inp_path, out_path)

    assert "! wB97X RKS def2-TZVP D4 RIJCOSX def2/J CPCM OPT numFREQ MOREAD" in content
    assert '%moinp "S0.gbw"' in content
    assert "nroots 15" in content
    assert "* xyz 0 1" in content
    assert content.count("* xyz") == 1


def test_the_original_input_file_is_never_rewritten(tmp_path):
    """The edit belongs to the isolated copy; the job directory keeps the
    input exactly as the workflow wrote it."""
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    before = inp_path.read_text(encoding="utf-8")
    _resume(tmp_path, inp_path, out_path)

    assert inp_path.read_text(encoding="utf-8") == before


def test_a_different_molecule_in_the_same_directory_is_not_resumed(tmp_path):
    """The element sequence is the gate. A reused job directory must not
    silently inherit somebody else's geometry."""
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.xyz").write_text(
        "3\nother molecule\n"
        "  N    0.100000    0.200000    0.300000\n"
        "  H    0.400000    0.500000    0.600000\n"
        "  H    0.700000    0.800000    0.900000\n",
        encoding="utf-8",
    )

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert notes == []
    assert "0.117300" in content


def test_a_finished_run_is_left_to_smart_recalc(tmp_path):
    inp_path, out_path = _partial_run(
        tmp_path, out="                 ****ORCA TERMINATED NORMALLY****\n"
    )
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert notes == []
    assert content == inp_path.read_text(encoding="utf-8")


def test_a_first_attempt_is_not_treated_as_a_resume(tmp_path):
    """No previous output means nothing to continue from — even if a stale
    .xyz happens to be lying around."""
    inp_path = tmp_path / "S1.inp"
    inp_path.write_text(_INLINE, encoding="utf-8")
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, notes = _resume_partial_orca_run(
        _INLINE, inp_path, tmp_path, tmp_path / "S1.out"
    )

    assert notes == []
    assert content == _INLINE


def test_a_job_without_opt_keeps_its_geometry(tmp_path):
    """A single point restarted from a half-optimised geometry would report an
    energy for coordinates nobody asked about."""
    single_point = _INLINE.replace(" OPT numFREQ", "")
    inp_path, out_path = _partial_run(tmp_path, inp=single_point)
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert notes == []
    assert "0.117300" in content


@pytest.mark.parametrize(
    "failure",
    [
        "ORCA finished by error termination in SCF\n",
        "aborting the run\n",
        "Error : ORCA finished by error termination in NUMFREQ\n",
    ],
)
def test_a_crashed_run_does_not_resume_from_its_own_bad_geometry(tmp_path, failure):
    """Being killed and having failed are different. Continuing from the last
    geometry of a run that blew up on that geometry reproduces the failure and
    burns the recovery attempts doing it."""
    inp_path, out_path = _partial_run(tmp_path, out=failure)
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert notes == []
    assert "0.117300" in content, "the input geometry must stand"


def test_a_crashed_run_may_still_reuse_a_finished_hessian(tmp_path):
    """ORCA validates the Hessian file itself, and a Hessian that completed
    before the failure is still a Hessian."""
    inp_path, out_path = _partial_run(
        tmp_path, out="ORCA finished by error termination in SCF\n"
    )
    (tmp_path / "S1.hess").write_text("$orca_hessian_file\n", encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert any("Hessian restart enabled" in n for n in notes), notes
    assert "restart true" in content


def test_an_xyzfile_reference_is_repointed(tmp_path):
    inp = """! HF-3c OPT
%base "S1"
* xyzfile 0 1 start.xyz
"""
    inp_path, out_path = _partial_run(tmp_path, inp=inp)
    (tmp_path / "start.xyz").write_text(
        "3\nstart\n"
        "  O    0.000000    0.000000    0.117300\n"
        "  H    0.000000    0.757200   -0.469200\n"
        "  H    0.000000   -0.757200   -0.469200\n",
        encoding="utf-8",
    )
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert any("geometry continued" in n for n in notes), notes
    assert str(tmp_path / "S1.xyz") in content


# --------------------------------------------------------------------------
# Hessian
# --------------------------------------------------------------------------

def test_an_existing_hessian_is_restarted_not_recomputed(tmp_path):
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.hess").write_text("$orca_hessian_file\n", encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert any("Hessian restart enabled" in n for n in notes), notes
    assert "%freq" in content and "restart true" in content


def test_the_hessian_travels_into_the_working_directory(tmp_path):
    """Regression, found by running it against ORCA 6.1.1 rather than by
    reading the code: ORCA derives the restart file from %base and looks for
    it in its working directory. Running isolated that directory is empty, and
    ORCA then prints "<<<RESTART CALCULATION>>>" followed by "No Restart data
    found - Restarting from scratch" — a restart that reports success and
    recomputes every displacement anyway."""
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.hess").write_text("$orca_hessian_file\n", encoding="utf-8")
    work = tmp_path / ".orca_iso_S1"
    work.mkdir()

    content, notes = _resume_partial_orca_run(
        inp_path.read_text(encoding="utf-8"), inp_path, tmp_path, out_path, work_dir=work
    )

    assert any("Hessian restart enabled" in n for n in notes), notes
    assert (work / "S1.hess").is_file(), "restart flag without the file is a lie"
    assert "restart true" in content


def test_no_working_directory_means_no_staging(tmp_path):
    """Callers that run in place already have the file next to the input."""
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.hess").write_text("$orca_hessian_file\n", encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert any("Hessian restart enabled" in n for n in notes), notes
    assert "restart true" in content


def test_the_restart_flag_joins_an_existing_freq_block(tmp_path):
    inp = _INLINE.replace("%tddft", "%freq\n  Temp 298.15\nend\n\n%tddft")
    inp_path, out_path = _partial_run(tmp_path, inp=inp)
    (tmp_path / "S1.hess").write_text("$orca_hessian_file\n", encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert any("Hessian restart enabled" in n for n in notes), notes
    assert content.count("%freq") == 1, "must not add a second %freq block"
    assert "Temp 298.15" in content
    assert "restart true" in content


def test_no_hessian_file_means_no_restart_flag(tmp_path):
    """ORCA writes no partial Hessian, so a frequency run killed midway has
    nothing to restart from and must simply run again."""
    inp_path, out_path = _partial_run(tmp_path)

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert "restart true" not in content
    assert not any("Hessian" in n for n in notes)


def test_a_restart_flag_is_never_added_twice(tmp_path):
    inp = _INLINE.replace("%tddft", "%freq\n  restart true\nend\n\n%tddft")
    inp_path, out_path = _partial_run(tmp_path, inp=inp)
    (tmp_path / "S1.hess").write_text("$orca_hessian_file\n", encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert content.count("restart true") == 1
    assert not any("Hessian" in n for n in notes)


# --------------------------------------------------------------------------
# robustness — resuming may never be the reason a job fails
# --------------------------------------------------------------------------

def test_only_the_first_job_of_a_compound_input_is_touched(tmp_path):
    inp = _INLINE + "\n$new_job\n! HF-3c OPT\n%base \"S1_check\"\n* xyz 0 1\n  O 9.0 9.0 9.0\n*\n"
    inp_path, out_path = _partial_run(tmp_path, inp=inp)
    (tmp_path / "S1.xyz").write_text(_RESUMED_XYZ, encoding="utf-8")

    content, _ = _resume(tmp_path, inp_path, out_path)

    assert "O 9.0 9.0 9.0" in content, "the second job must be left alone"


@pytest.mark.parametrize(
    "broken",
    ["", "not a number\n", "3\nonly a comment\n", "99\ncomment\n  O 0.0 0.0 0.0\n"],
)
def test_an_unusable_restart_geometry_is_ignored(tmp_path, broken):
    inp_path, out_path = _partial_run(tmp_path)
    (tmp_path / "S1.xyz").write_text(broken, encoding="utf-8")

    content, notes = _resume(tmp_path, inp_path, out_path)

    assert notes == []
    assert content == inp_path.read_text(encoding="utf-8")


def test_an_unreadable_job_directory_does_not_raise(tmp_path):
    """Called from the ORCA runner: it must degrade to a normal run, never
    take the job down with it."""
    inp_path, out_path = _partial_run(tmp_path)

    content, notes = _resume_partial_orca_run(
        _INLINE, inp_path, tmp_path / "does-not-exist", out_path
    )

    assert notes == []
    assert content == _INLINE
