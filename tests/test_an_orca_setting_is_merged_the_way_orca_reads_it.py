"""A setting from CONTROL merged into an ORCA job the way ORCA itself reads the job.

The rules are the ORCA 6.1.1 manual's (section 2.1) and what ORCA 6.1.1 did
when each was tried:

* ``! TightSCF VeryTightSCF`` runs VeryTight in either order, and two
  functionals are an input error -- so a keyword replaces the member of its
  family the job has instead of being appended after it.
* ``%scf maxiter 300 tole 1e-9 end`` sets both variables, and a repeated block
  is read variable by variable with the later value winning -- so a merge
  sees each variable, and folds a repeated block before changing it.
* ``SOSCF`` alone on a line in ``%scf`` opens a sub-block.
* A bare ``%cpcm`` block switches CPCM on (the energy moved by 8.9 mEh),
  while ``%scf``, ``%geom``, ``%freq``, ``%method``, ``%basis``, ``%output``,
  ``%mp2``, ``%shark`` and ``%symmetry`` alone start nothing.

Every job of the archive reads back byte for byte when nothing is merged
(36 531 of 36 539; the rest are %compound jobs and one broken file, both
refused).
"""

from __future__ import annotations

import pytest

from delfin.common import orca_input as oi

# The shape of an OCCUPIER broken-symmetry input in the archive: two one-line
# %scf blocks and a metal basis on the metal's coordinate line.
OCCUPIER_BS = """! PBE0 def2-SVP D4 RIJCOSX def2/J CPCM(dmf) OPT FREQ PModel MOREAD TightSCF
%moinp "input_red_step_1_OCCUPIER.gbw"
%maxcore 6000
%pal nprocs 18 end
%scf maxiter 125 end
%scf BrokenSym 1,1 end
%freq
  Temp 298.15
end
%output
print[p_mos] 1
print[p_basis] 2
end
* xyz -1 1
C          -1.11736974265122      3.57845588429044      0.53176128123210
Co         -0.02912657244695     -0.34055506228113     -0.00413799049569   NewGTO "def2-TZVP" end
N           1.18159758310338     -1.73604082305673      0.43982649453832
*
"""

ESD_TWO_JOBS = """! CAM-B3LYP UKS def2-TZVPP D4 RIJCOSX def2/J CPCM(chcl3) OPT MOREAD
%base "T1"
%moinp "S0.gbw"
%pal nprocs 20 end
%maxcore 6000
* xyz 0 3
C  0.0 0.0 0.0
O  0.0 0.0 1.2
*

$new_job
! CAM-B3LYP RKS def2-TZVPP D4 RIJCOSX def2/J CPCM(chcl3)
%base "T1_TDDFT"
%pal nprocs 20 end
%maxcore 6000

%tddft
  nroots 15
  maxdim 7
  tda FALSE
  triplets true
end

* xyzfile 0 1 T1.xyz
"""


def _block(text, name):
    items = oi.parse_job(text)
    return {st.key: " ".join(st.lines[0].split()[1:]) for it in items
            if it.kind == "block" and it.name == name for st in it.block.statements if st.key}


# ------------------------------------------------------------------ reading

@pytest.mark.parametrize("text", [OCCUPIER_BS, *oi.split_jobs(ESD_TWO_JOBS), """! PBE def2-SVP
%geom
  Constraints
    { B 0 1 C }
  end
  maxiter 50
end
%scf
  SOSCF Start 0.1 MaxIt 5 end
end
%basis
  NewGTO Co "def2-TZVP" end
end
%plots
  ElDens("x.cube");
end
* xyz 0 1
H 0 0 0
H 0 0 0.74 *
"""])
def test_a_job_reads_back_exactly_when_nothing_is_merged(text):
    assert oi.render_job(oi.parse_job(text)) == text


def test_the_jobs_of_a_file_are_kept_apart_and_whole():
    jobs = oi.split_jobs(ESD_TWO_JOBS)

    assert len(jobs) == 2
    assert "".join(jobs) == ESD_TWO_JOBS
    assert jobs[1].startswith("$new_job")


def test_a_line_with_two_assignments_is_two_variables():
    # measured: ORCA sets MaxIter 300 and TolE 1e-9 from one line
    merged, notes = oi.apply_to_job("! PBE def2-SVP\n%scf maxiter 125 tole 1e-8 end\n* xyz 0 1\nH 0 0 0\n*\n",
                                    additions=["%scf tole 1e-9 end"])

    assert _block(merged, "scf") == {"maxiter": "125", "tole": "1e-9"}
    assert notes == []


def test_what_cannot_be_taken_apart_is_refused_not_guessed():
    for text in ("%compound\n  New_Step\n  ! B3LYP\n  Step_End\nend\n",
                 "! PBE\n%geom maxiter 10\n%pal nprocs 2 end\n* xyz 0 1\nH 0 0 0\n*\n"):
        with pytest.raises(oi.OrcaInputError):
            oi.parse_job(text)


# ------------------------------------------------------------------ blocks

def test_a_variable_replaces_the_jobs_own_and_a_new_one_is_added():
    merged, _ = oi.apply_to_job(OCCUPIER_BS, additions=["%scf maxiter 400 end", "%geom\n  maxiter 250\nend"])

    assert _block(merged, "scf") == {"maxiter": "400", "brokensym": "1,1"}
    assert _block(merged, "geom") == {"maxiter": "250"}
    # the geometry and everything else is where it was
    assert 'Co         -0.02912657244695     -0.34055506228113     -0.00413799049569   NewGTO "def2-TZVP" end' in merged
    assert merged.index("%geom") < merged.index("* xyz")


def test_a_repeated_block_is_folded_as_orca_reads_it_before_it_is_changed():
    text = "! PBE\n%scf maxiter 125 end\n%scf BrokenSym 1,1\n maxiter 200\nend\n* xyz 0 1\nH 0 0 0\n*\n"
    merged, _ = oi.apply_to_job(text, additions=["%scf tole 1e-9 end"])

    assert merged.count("%scf") == 1
    # the later block's maxiter is what ORCA used; it stays the value
    assert _block(merged, "scf") == {"maxiter": "200", "brokensym": "1,1", "tole": "1e-9"}


def test_merging_twice_changes_nothing():
    once, _ = oi.apply_to_job(OCCUPIER_BS, keywords=["VeryTightSCF", "DEFGRID3"],
                              additions=["%scf maxiter 400 end"])
    twice, _ = oi.apply_to_job(once, keywords=["VeryTightSCF", "DEFGRID3"], additions=["%scf maxiter 400 end"])

    assert twice == once


def test_a_bare_keyword_is_never_written_as_a_value():
    # measured: 'SOSCF' alone opens a sub-block and the next lines are read inside it
    items = oi.parse_job("! PBE\n* xyz 0 1\nH 0 0 0\n*\n")
    oi.set_block_values(items, "scf", {"CNVSOSCF": True, "MaxIter": 800})

    assert _block(oi.render_job(items), "scf") == {"cnvsoscf": "true", "maxiter": "800"}
    with pytest.raises(oi.OrcaInputError):
        oi.set_block_values(items, "scf", {"SOSCFStart": ""})


# ---------------------------------------------------------------- keywords

def test_a_keyword_takes_the_place_of_its_family_member():
    merged, notes = oi.apply_to_job(OCCUPIER_BS, keywords=["VeryTightSCF", "TightOpt", "B3LYP", "D3BJ"])
    bang = merged.splitlines()[0].split()

    assert bang[:4] == ["!", "B3LYP", "def2-SVP", "D3BJ"]
    assert "TightOpt" in bang and "OPT" not in bang
    assert "VeryTightSCF" in bang and "TightSCF" not in bang
    assert notes == []


def test_a_run_type_is_never_added_to_a_job_that_does_not_run_it():
    single_point = "! PBE0 def2-SVP TightSCF\n* xyz 0 1\nH 0 0 0\n*\n"
    merged, notes = oi.apply_to_job(single_point, keywords=["TightOpt", "NumFreq"])

    assert merged == single_point
    assert any("no optimisation level" in n for n in notes)
    assert any("no frequency" in n for n in notes)


def test_for_many_jobs_a_method_keyword_only_replaces():
    xtb = "! XTB2 OPT\n%pal nprocs 4 end\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"
    merged, notes = oi.apply_to_job(xtb, keywords=["B3LYP", "D4", "VeryTightSCF"], many_jobs=True)

    assert merged.splitlines()[0] == "! XTB2 OPT VeryTightSCF"
    assert any("B3LYP not added" in n for n in notes)
    assert any("D4 not added" in n for n in notes)


def test_a_guess_is_not_set_on_a_job_that_reads_the_previous_orbitals():
    merged, notes = oi.apply_to_job(OCCUPIER_BS, keywords=["Hueckel"])

    assert merged == OCCUPIER_BS
    assert any("MORead" in n for n in notes)


# ------------------------------------------------------- what is DELFIN's

@pytest.mark.parametrize("addition, word", [
    ("%pal nprocs 4 end", "%pal"),
    ("%maxcore 9000", "%maxcore"),
    ('%base "other"', "%base"),
    ('%moinp "x.gbw"', "%moinp"),
    ("%tddft iroot 2 end", "iroot"),
    ("%foo bar 1 end", "not an ORCA input block"),
])
def test_what_delfin_decides_per_job_is_left_alone(addition, word):
    merged, notes = oi.apply_to_job(OCCUPIER_BS, additions=[addition])

    assert merged == OCCUPIER_BS
    assert any(word in n for n in notes)


@pytest.mark.parametrize("token", ["PAL8", "MORead", "OptTS", "SP", "ESD(FLUOR)"])
def test_a_keyword_that_decides_what_a_job_is_is_refused(token):
    merged, notes = oi.apply_to_job(OCCUPIER_BS, keywords=[token])

    assert merged == OCCUPIER_BS
    assert notes


def test_for_many_jobs_tddft_roots_have_one_spelling():
    merged, notes = oi.apply_to_job(ESD_TWO_JOBS.split("$new_job")[1], additions=["%tddft nroots 20 end"],
                                    many_jobs=True)

    assert "nroots 15" in merged
    assert any("TDDFT_nroots" in n for n in notes)


def test_for_many_jobs_a_block_that_starts_a_calculation_only_joins_jobs_that_run_it():
    gas_opt = "! PBE0 def2-SVP OPT\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"
    merged, _ = oi.apply_to_job(gas_opt, many_jobs=True, additions=[
        "%tddft\n  maxiter 300\nend", "%cpcm\n  epsilon 80\nend", "%scf\n  maxiter 300\nend"])

    assert "%tddft" not in merged
    assert "%cpcm" not in merged          # a bare %cpcm would switch solvation on
    assert "%scf" in merged
    solvated, _ = oi.apply_to_job(OCCUPIER_BS, many_jobs=True, additions=["%cpcm\n  epsilon 80\nend"])
    assert "%cpcm" in solvated            # CPCM(dmf) runs; the block only tunes it


# -------------------------------------------------------------- coordinates

def test_new_coordinates_keep_the_metal_basis_on_its_atom():
    items = oi.parse_job(OCCUPIER_BS)
    oi.with_coordinates(items, [("C", 1.0, 2.0, 3.0), ("Co", 0.1, 0.2, 0.3), ("N", -1.0, -2.0, -3.0)])
    text = oi.render_job(items)

    co = next(line for line in text.splitlines() if line.split()[:1] == ["Co"])
    assert co.split()[1:4] == ["0.10000000", "0.20000000", "0.30000000"]
    assert co.endswith('NewGTO "def2-TZVP" end')


def test_coordinates_for_other_atoms_are_refused():
    items = oi.parse_job(OCCUPIER_BS)
    with pytest.raises(oi.OrcaInputError):
        oi.with_coordinates(items, [("C", 0, 0, 0), ("Ni", 0, 0, 0), ("N", 0, 0, 0)])


def test_an_xyzfile_reference_becomes_inline_coordinates():
    items = oi.parse_job(oi.split_jobs(ESD_TWO_JOBS)[1])
    oi.with_coordinates(items, [("C", 0.0, 0.0, 0.0), ("O", 0.0, 0.0, 1.21)])
    text = oi.render_job(items)

    assert "* xyzfile" not in text
    assert "* xyz 0 1" in text and text.rstrip().endswith("*")
