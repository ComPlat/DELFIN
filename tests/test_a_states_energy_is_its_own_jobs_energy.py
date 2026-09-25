"""An output that holds two calculations must be read where it speaks of itself.

DELFIN writes ``T1.inp`` as two calculations in one file: the triplet is
optimised under ``%base "T1"``, then a closed-shell TD-DFT check job is
appended under ``%base "T1_TDDFT"``.  The last ``FINAL SINGLE POINT ENERGY``
in ``T1.out`` therefore belongs to the check job -- the ground state at the
triplet geometry -- and not to the triplet.

Read that way, an archived phosphorescence run came out at DELE = 23162 cm-1
where the adiabatic difference is 13667 cm-1: 1.18 eV, written straight into
the ESD input.  ``T1.property.txt`` had the right number the whole time,
because ORCA gives each ``%base`` its own property file -- the manual requires
those labels for ``$new_job`` precisely so the jobs can be told apart.  So the
energy of ``<name>.out`` is the energy of the job whose base is ``<name>``,
which is what ORCA's own file naming says.

``S0.out`` had this fixed once before, by a special case keyed on the literal
filename; these tests keep that behaviour while the rule now covers every
state.  Nothing here touches real calculations: the outputs are written out
below in the shape ORCA prints them.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.energies import find_electronic_energy

OPT_ENERGY = -1538.720536103183          # the optimised triplet, job 1
CHECK_ENERGY = -1538.677272446181        # the closed-shell check job, job 2
GROUND_ENERGY = -1538.782807769000       # the optimised ground state


def _job_body(energies, *, optimised=False):
    """The part ORCA prints per calculation, ending in its final energy.

    With *optimised* the markers an optimisation leaves are included, which
    is what tells a reader a cycle count means anything.
    """
    lines = []
    for step, energy in enumerate(energies, start=1):
        if optimised:
            lines.append(f'         GEOMETRY OPTIMIZATION CYCLE {step}')
        lines += [
            f'Total Energy       :    {energy:18.10f} Eh',
            f'FINAL SINGLE POINT ENERGY     {energy:.12f}',
            '',
        ]
    if optimised:
        lines += ['       *** THE OPTIMIZATION HAS CONVERGED ***', '']
    return '\n'.join(lines)


def _output(echo: str, jobs, *, optimised=False) -> str:
    """An ORCA output: the echoed input, then one section per job.

    ``jobs`` is a list of energy lists, one per ``$new_job``.  A single entry
    produces an output with no JOB NUMBER banners at all, which is what ORCA
    prints when the input holds one calculation.  With *optimised* the first
    job carries an optimisation's markers and the rest do not, which is the
    shape DELFIN writes: an optimisation, then a single-point check job.
    """
    head = ['',
            '                                 INPUT FILE',
            '================================================================================',
            'NAME = job.inp']
    head += [f'|{n + 1:3}> {line}' for n, line in enumerate(echo.strip('\n').split('\n'))]
    head += ['                          ****END OF INPUT****', '']

    body = []
    if len(jobs) == 1:
        body.append(_job_body(jobs[0], optimised=optimised))
    else:
        for number, energies in enumerate(jobs, start=1):
            body += [
                f'                 $$$$$$$$$$$$$$$$  JOB NUMBER  {number} $$$$$$$$$$$$$$',
                '',
                _job_body(energies, optimised=optimised and number == 1),
            ]
    body.append('                             ****ORCA TERMINATED NORMALLY****')
    return '\n'.join(head + body)


COMPOUND_T1_ECHO = """
! PBE0 UKS def2-TZVP D4 RIJCOSX def2/J CPCM OPT MOREAD
%base "T1"
%moinp "S0.gbw"

* xyz 0 3
  O 0.0 0.0 0.0
*

$new_job
! PBE0 RKS def2-TZVP D4 RIJCOSX def2/J CPCM
%base "T1_TDDFT"

%tddft
  nroots 15
  triplets true
end

* xyzfile 0 1 T1.xyz
"""

COMPOUND_S0_ECHO = COMPOUND_T1_ECHO.replace('"T1"', '"S0"') \
                                   .replace('"T1_TDDFT"', '"S0_TDDFT"') \
                                   .replace('T1.xyz', 'S0.xyz') \
                                   .replace('xyz 0 3', 'xyz 0 1')


def _write(tmp_path: Path, name: str, text: str) -> str:
    path = tmp_path / name
    path.write_text(text, encoding='utf-8')
    return str(path)


# ---------------------------------------------------------------------------
# the job the file is named after
# ---------------------------------------------------------------------------

def test_a_states_energy_is_the_energy_of_its_own_job(tmp_path):
    out = _write(tmp_path, 'T1.out',
                 _output(COMPOUND_T1_ECHO, [[OPT_ENERGY - 1e-4, OPT_ENERGY],
                                            [CHECK_ENERGY]]))

    assert find_electronic_energy(out) == pytest.approx(OPT_ENERGY, abs=1e-9)


def test_the_ground_state_output_is_read_as_it_always_was(tmp_path):
    """S0.out was already protected by a special case; it must not regress."""
    out = _write(tmp_path, 'S0.out',
                 _output(COMPOUND_S0_ECHO, [[GROUND_ENERGY - 1e-4, GROUND_ENERGY],
                                            [CHECK_ENERGY]]))

    assert find_electronic_energy(out) == pytest.approx(GROUND_ENERGY, abs=1e-9)


def test_a_state_written_by_the_hybrid_path_is_covered_too(tmp_path):
    """The old special case matched the string "S0.out" and nothing else."""
    echo = COMPOUND_T1_ECHO.replace('"T1"', '"T2_second"')
    out = _write(tmp_path, 'T2_second.out',
                 _output(echo, [[OPT_ENERGY], [CHECK_ENERGY]]))

    assert find_electronic_energy(out) == pytest.approx(OPT_ENERGY, abs=1e-9)


# ---------------------------------------------------------------------------
# and everything else reads as before
# ---------------------------------------------------------------------------

def test_a_single_job_output_is_read_end_to_end(tmp_path):
    echo = '! PBE0 RKS def2-TZVP OPT numFREQ\n%base "T1"\n* xyz 0 1\n  O 0.0 0.0 0.0\n*'
    out = _write(tmp_path, 'T1.out', _output(echo, [[OPT_ENERGY - 1e-3, OPT_ENERGY]]))

    assert find_electronic_energy(out) == pytest.approx(OPT_ENERGY, abs=1e-9)


def test_a_second_job_without_a_base_of_its_own_still_counts(tmp_path):
    """Naming no base keeps the previous one, so both jobs are this file's."""
    echo = COMPOUND_T1_ECHO.replace('%base "T1_TDDFT"\n', '')
    out = _write(tmp_path, 'T1.out', _output(echo, [[OPT_ENERGY], [CHECK_ENERGY]]))

    assert find_electronic_energy(out) == pytest.approx(CHECK_ENERGY, abs=1e-9)


def test_an_output_no_job_claims_is_read_as_before(tmp_path):
    """A renamed file cannot be attributed, so the previous reading stands."""
    out = _write(tmp_path, 'irgendwas.out',
                 _output(COMPOUND_T1_ECHO, [[OPT_ENERGY], [CHECK_ENERGY]]))

    assert find_electronic_energy(out) == pytest.approx(CHECK_ENERGY, abs=1e-9)


def test_a_missing_file_is_still_no_energy(tmp_path):
    assert find_electronic_energy(str(tmp_path / 'nicht_da.out')) is None


# ---------------------------------------------------------------------------
# what it was all for: the adiabatic difference in the ESD input
# ---------------------------------------------------------------------------

def test_the_adiabatic_difference_is_between_the_two_optimised_states(tmp_path):
    from delfin.esd_input_generator import HARTREE_TO_CM1, calculate_dele_cm1

    t1 = _write(tmp_path, 'T1.out',
                _output(COMPOUND_T1_ECHO, [[OPT_ENERGY], [CHECK_ENERGY]]))
    s0 = _write(tmp_path, 'S0.out',
                _output(COMPOUND_S0_ECHO, [[GROUND_ENERGY], [CHECK_ENERGY]]))

    dele = calculate_dele_cm1(t1, s0)

    assert dele == pytest.approx((OPT_ENERGY - GROUND_ENERGY) * HARTREE_TO_CM1, abs=0.5)
    wrong = (CHECK_ENERGY - GROUND_ENERGY) * HARTREE_TO_CM1
    assert abs(dele - wrong) > 9000, 'the check job is back in the DELE'


def test_a_trajectory_counts_the_cycles_of_its_own_job_only(tmp_path):
    """The check job's single point was read as one more optimisation cycle.

    Measured on an archived T1.out: 10 cycles reported where the optimisation
    ran 9, and a final energy 1.18 eV above the triplet's.
    """
    from delfin.api import extract_optimization_trajectory

    walked = [OPT_ENERGY - 1e-2, OPT_ENERGY - 1e-3, OPT_ENERGY]
    (tmp_path / 'T1.out').write_text(
        _output(COMPOUND_T1_ECHO, [walked, [CHECK_ENERGY]], optimised=True),
        encoding='utf-8')

    result = extract_optimization_trajectory(str(tmp_path))

    assert result.output_file == 'T1.out'
    assert result.error is None, result.error
    assert result.n_cycles == len(walked), 'the check job was counted as a cycle'
    assert result.final_energy_eh == pytest.approx(OPT_ENERGY, abs=1e-9)
    assert result.converged is True


def test_a_trajectory_of_a_single_job_output_is_unchanged(tmp_path):
    from delfin.api import extract_optimization_trajectory

    echo = '! PBE0 RKS def2-TZVP OPT\n%base "S0"\n* xyz 0 1\n  O 0.0 0.0 0.0\n*'
    walked = [GROUND_ENERGY - 1e-2, GROUND_ENERGY]
    (tmp_path / 'S0.out').write_text(_output(echo, [walked], optimised=True),
                                     encoding='utf-8')

    result = extract_optimization_trajectory(str(tmp_path))

    assert result.n_cycles == len(walked)
    assert result.final_energy_eh == pytest.approx(GROUND_ENERGY, abs=1e-9)


def test_the_state_input_really_is_written_as_two_jobs(tmp_path):
    """The reader exists because the writer appends a check job -- keep them paired."""
    from delfin.esd_input_generator import create_state_input

    (tmp_path / 'S0.xyz').write_text('1\n\nO 0.0 0.0 0.0\n', encoding='utf-8')
    (tmp_path / 'T1.xyz').write_text('1\n\nO 0.0 0.0 0.0\n', encoding='utf-8')

    written = create_state_input(
        state='T1', esd_dir=tmp_path, charge=0, solvent='', metals=[],
        main_basisset='def2-TZVP', metal_basisset='def2-TZVP',
        config={'functional': 'PBE0'},
    )
    text = (tmp_path / Path(written).name).read_text(encoding='utf-8')

    assert '$new_job' in text
    assert '%base "T1"' in text
    assert '%base "T1_TDDFT"' in text, (
        'if the check job stops being appended, this reader has nothing to do')
