"""Parity: OPI IMAG inputs == classic IMAG inputs.

``delfin.tools.imag_opi`` reuses the IMAG writer's input transforms BY IMPORT
(``delfin.imag._candidate_job``, ``_split_first_job``, ``_with_base``,
``_with_geometry`` — delfin/imag.py:222-295, the same functions the classic
flow calls at imag.py:410/419 and :551-554) and re-emits the transformed job
with OPI's input classes (``opi.input``).  The acceptance criterion is parity:
for the same source input and the same displacement, the OPI file must be
structurally equivalent to the classic one — same simple keywords
(order/case-insensitive), same blocks with the same options and values, same
charge/multiplicity/coordinates (per-atom ``NewGTO`` included), same
``%pal``/``%maxcore``, same file names.  The comparison is
``delfin.tools.occupier_opi.normalise_inp`` / ``inputs_equivalent`` — the
normalisation the OCCUPIER parity test uses (tests/test_opi_occupier_parity.py),
reused so the contract cannot drift between the writers.

The classic module stays the workflow's driver; this module only re-emits
inputs (no ORCA launch).
"""

from __future__ import annotations

import pytest

# The main py3.11 CI jobs run all tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from delfin.tools.imag_opi import (  # noqa: E402
    write_imag_candidate_opi_input,
    write_imag_rerun_opi_input,
)
from delfin.tools.occupier_opi import (  # noqa: E402
    inputs_equivalent,
    normalise_inp,
)

import delfin.imag as classic  # noqa: E402


# A small metal-complex-like first job built from ideal geometry (never a
# CSD/crystal structure): Fe with N ligands plus H capping atoms, per-atom
# NewGTO on the metal, %base present, an appended $new_job with an xyzfile
# reference — the shapes the classic writer must survive.
_STEP_INPUT = """! PBE0 def2-SVP D4 RIJCOSX def2/J TightSCF OPT FREQ
%scf
  maxiter 150
end
%tddft
  IRoot 1
end
%base "initial"
%maxcore 1000
%pal nprocs 4 end
* xyz 0 1
Fe  0.00000000  0.00000000  0.00000000 NewGTO "def2-TZVP" end
N   2.00000000  0.00000000  0.00000000
N  -2.00000000  0.00000000  0.00000000
N   0.00000000  2.00000000  0.00000000
N   0.00000000 -2.00000000  0.00000000
H   2.50000000  1.00000000  0.00000000
H  -2.50000000 -1.00000000  0.00000000
H   0.00000000  0.00000000  2.50000000
*
$new_job
! PBE0 def2-SVP D4 RIJCOSX def2/J
%pal nprocs 4 end
* xyzfile 0 1 initial.xyz
"""

# The classic single-point candidate input shape: optimisation and frequency
# taken off the '!' line, %base renamed, coordinates replaced (inline rows
# keep their tails, so the metal's NewGTO survives), appended jobs untouched
# by the candidate writer (only the first job is written).
_CAND_SOURCE = """! PBE0 def2-SVP D4 RIJCOSX def2/J TightSCF OPT FREQ
%scf
  maxiter 150
end
%tddft
  IRoot 1
end
%base "initial"
%maxcore 1000
%pal nprocs 4 end
* xyz 0 1
Fe  0.00000000  0.00000000  0.00000000 NewGTO "def2-TZVP" end
N   2.00000000  0.00000000  0.00000000
N  -2.00000000  0.00000000  0.00000000
N   0.00000000  2.00000000  0.00000000
N   0.00000000 -2.00000000  0.00000000
H   2.50000000  1.00000000  0.00000000
H  -2.50000000 -1.00000000  0.00000000
H   0.00000000  0.00000000  2.50000000
*
"""

_DISPLACED = [
    ("Fe", 0.11100000, 0.00000000, 0.00000000),
    ("N", 2.00000000, 0.00000000, 0.00000000),
    ("N", -2.00000000, 0.00000000, 0.00000000),
    ("N", 0.00000000, 2.00000000, 0.00000000),
    ("N", 0.00000000, -2.00000000, 0.00000000),
    ("H", 2.50000000, 1.00000000, 0.00000000),
    ("H", -2.50000000, -1.00000000, 0.00000000),
    ("H", 0.00000000, 0.00000000, 2.50000000),
]


def _classic_candidate(atoms, base, *, optimise=False, source=_CAND_SOURCE):
    return classic._candidate_job(source, atoms, base, optimise=optimise)


# ----------------------------------------------------------------- matrix --
def _control_matrix():
    """The parity matrix.  The source-input shapes are the ones the
    repository's IMAG flow transforms: the flat PBE0/TightSCF/OPT/FREQ job of
    examples' IMAG CONTROLs (e.g. examples/DELFIN_Paper_Calcualtions/
    FeTPP_classic/CONTROL.txt, IMAG=yes), per-atom NewGTO for the metal basis
    the OCCUPIER examples use (examples/DELFIN_Paper_Calcualtions/
    FeTPP_OCCUPIER/CONTROL.txt: main_basisset=def2-SVP,
    metal_basisset=def2-TZVP), %scf maxiter (maxiter_occupier=200 there),
    an ESD-style %tddft IRoot block, and the appended $new_job with an
    xyzfile geometry reference."""
    return [
        # name, source text, atoms, base, optimise
        ("pbe0_flat", _CAND_SOURCE, _DISPLACED, "initial_imag_pos_a1", False),
        ("pbe0_opt", _CAND_SOURCE, _DISPLACED, "initial_imag_pos_a1", True),
        ("no_base", _CAND_SOURCE.replace('%base "initial"\n', ""),
         _DISPLACED, "initial", False),
        ("compact_pal", _CAND_SOURCE.replace("%pal nprocs 4 end", ""),
         _DISPLACED, "initial", False),
    ]


@pytest.mark.parametrize("name,source,atoms,base,optimise", _control_matrix())
def test_imag_candidate_parity(tmp_path, name, source, atoms, base, optimise):
    classic_text = _classic_candidate(atoms, base, optimise=optimise, source=source)
    opi_path = write_imag_candidate_opi_input(
        work_dir=tmp_path / "opi",
        source_text=source,
        atoms=atoms,
        base=base,
        optimise=optimise,
    )
    opi_text = opi_path.read_text()
    diffs = inputs_equivalent(
        normalise_inp(classic_text),
        normalise_inp(opi_text),
    )
    assert diffs == [], (
        f"IMAG candidate ({name}): OPI input differs from classic input\n"
        + "\n".join(diffs)
        + f"\n--- classic ---\n{classic_text}"
        + f"\n--- opi ---\n{opi_text}"
    )


def test_imag_rerun_parity(tmp_path):
    """The re-optimisation input: first job at the chosen geometry, its %base
    forced to the step's base name, appended jobs carried verbatim
    (imag.py:551-554)."""
    classic_text = classic._with_base(
        classic._with_geometry(classic._split_first_job(_STEP_INPUT)[0], _DISPLACED),
        "initial",
    ) + classic._split_first_job(_STEP_INPUT)[1]
    opi_path = write_imag_rerun_opi_input(
        work_dir=tmp_path / "opi",
        source_text=_STEP_INPUT,
        atoms=_DISPLACED,
        base="initial",
    )
    diffs = inputs_equivalent(
        normalise_inp(classic_text),
        normalise_inp(opi_path.read_text()),
    )
    assert diffs == [], (
        "IMAG rerun: OPI input differs from classic input\n"
        + "\n".join(diffs)
        + f"\n--- classic ---\n{classic_text}"
        + f"\n--- opi ---\n{opi_path.read_text()}"
    )


def test_imag_candidate_names_follow_the_classic(tmp_path):
    """Candidate names are the classic ones: <base>.inp where base is
    initial_imag_<side>_a<n> (imag.py:416, 419)."""
    p = write_imag_candidate_opi_input(
        work_dir=tmp_path,
        source_text=_CAND_SOURCE,
        atoms=_DISPLACED,
        base="initial_imag_neg_a2",
        optimise=False,
    )
    assert p.name == "initial_imag_neg_a2.inp"


def test_imag_unknown_scf_option_is_refused_not_dropped(tmp_path):
    """An %scf option OPI's BlockScf cannot represent must raise, never be
    silently dropped (occupier parity precedent)."""
    with pytest.raises(ValueError):
        write_imag_candidate_opi_input(
            work_dir=tmp_path,
            source_text=_CAND_SOURCE.replace(
                "  maxiter 150", "  NotAnOption 7"),
            atoms=_DISPLACED,
            base="initial",
            optimise=False,
        )
