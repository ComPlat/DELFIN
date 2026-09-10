"""CONTROL's keyword:<job>= and additions:<job>= reach the jobs they name, merged as ORCA reads them.

Before, a keyword was appended to the '!' line and a block pasted in right
after it.  ORCA does not take the last of two SCF convergence keywords, it
takes the tighter (measured on 6.1.1), and DELFIN's own %scf further down
set MaxIter again -- so an override could add, but not change.  Now it
replaces what the job has, and the dashboard's check reads the entries with
the same reader the run uses.

The archive holds no CONTROL file with a real override (590 have the empty
template lines), so the template lines must keep validating, and nothing
that validated before may become an error.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from delfin import orca
from delfin.common import orca_overrides
from delfin.config import get_orca_override_hints, validate_control_text
from delfin.define import TEMPLATE

OPT = """! PBE0 def2-SVP D4 RIJCOSX def2/J CPCM(dmf) OPT FREQ TightSCF
%maxcore 6000
%pal nprocs 18 end
%scf maxiter 125 end
* xyz 0 1
H 0 0 0
H 0 0 0.74
*
"""

ESD = """! PBE0 def2-SVP D4 RIJCOSX def2/J TightSCF OPT
%base "S1"
%pal nprocs 8 end
%tddft
  nroots 15
  iroot 1
end
* xyz 0 1
H 0 0 0
H 0 0 0.74
*

$new_job
! PBE0 def2-SVP D4 RIJCOSX def2/J
%base "S1_TDDFT"
%pal nprocs 8 end
%tddft
  nroots 15
end
* xyzfile 0 1 S1.xyz
"""


@pytest.fixture(autouse=True)
def _fresh_caches():
    orca._CONTROL_OVERRIDE_CACHE.clear()
    orca._CONTROL_FILE_LOCATION_CACHE.clear()
    yield
    orca._CONTROL_OVERRIDE_CACHE.clear()
    orca._CONTROL_FILE_LOCATION_CACHE.clear()


def _run(tmp_path: Path, control: str, name: str, text: str, folder: str = "") -> str:
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "CONTROL.txt").write_text(control)
    where = tmp_path / folder if folder else tmp_path
    where.mkdir(parents=True, exist_ok=True)
    inp = where / name
    inp.write_text(text)
    orca._apply_control_overrides_to_input(inp, where)
    return inp.read_text()


def test_the_template_lines_change_nothing(tmp_path):
    assert _run(tmp_path, "keyword:basename=[]\nadditions:basename=[]\n", "initial.inp", OPT) == OPT


def test_a_keyword_replaces_the_one_of_its_kind(tmp_path):
    out = _run(tmp_path, "keyword:initial=[VeryTightSCF DEFGRID3]\n", "initial.inp", OPT)

    assert out.splitlines()[0] == "! PBE0 def2-SVP D4 RIJCOSX def2/J CPCM(dmf) OPT FREQ VeryTightSCF DEFGRID3"


def test_a_block_setting_changes_the_jobs_own_value(tmp_path):
    out = _run(tmp_path, 'additions:initial=["%scf maxiter 500 end"]\n', "initial.inp", OPT)

    assert out.count("%scf") == 1
    assert "maxiter 500" in out and "maxiter 125" not in out


def test_orca_input_can_be_written_between_brackets_over_several_lines(tmp_path):
    control = "additions:initial=[\n%scf\n  BrokenSym 1,1\n  maxiter 400\nend\n]\n"
    out = _run(tmp_path, control, "initial.inp", OPT)

    assert "BrokenSym 1,1" in out and "maxiter 400" in out
    assert validate_control_text(TEMPLATE + control) == validate_control_text(TEMPLATE)


def test_all_patterns_and_a_name_reach_a_job_most_specific_last(tmp_path):
    control = (
        'additions:all=["%scf maxiter 300 end"]\n'
        'additions:S*=["%scf maxiter 400 end"]\n'
        'additions:S1_TDDFT=["%scf maxiter 500 end"]\n'
        "keyword:all=[DEFGRID3]\n"
    )
    out = _run(tmp_path, control, "S1.inp", ESD, folder="ESD")
    first, second = out.split("$new_job")

    assert "maxiter 400" in first            # S1: all, then S*
    assert "maxiter 500" in second           # S1_TDDFT: all, S*, then its own name
    assert "DEFGRID3" in first.splitlines()[0] and "DEFGRID3" in second.split("\n")[1]


def test_the_esd_modules_roots_and_resources_stay_its_own(tmp_path, caplog):
    control = 'additions:S1=["%tddft iroot 3 end", "%pal nprocs 2 end"]\nkeyword:S1=[PAL4]\n'
    with caplog.at_level(logging.WARNING, logger="delfin.orca"):
        out = _run(tmp_path, control, "S1.inp", ESD, folder="ESD")

    assert out == ESD
    assert "iroot" in caplog.text and "%pal" in caplog.text and "PAL4" in caplog.text


def test_an_occupier_folder_names_every_run_in_it(tmp_path):
    out = _run(tmp_path / "a", "keyword:initial_OCCUPIER=[DEFGRID3]\n", "input2.inp", OPT, folder="initial_OCCUPIER")
    untouched = _run(tmp_path / "b", "keyword:initial=[DEFGRID3]\n", "input3.inp", OPT, folder="initial_OCCUPIER")

    assert "DEFGRID3" in out.splitlines()[0]
    assert untouched == OPT


def test_a_rerun_of_the_same_input_is_unchanged(tmp_path):
    control = 'keyword:all=[VeryTightSCF]\nadditions:initial=["%scf maxiter 500 end"]\n'
    once = _run(tmp_path, control, "initial.inp", OPT)
    inp = tmp_path / "initial.inp"
    orca._apply_control_overrides_to_input(inp, tmp_path)

    assert inp.read_text() == once


def test_a_recovery_retry_keeps_what_recovery_changed(tmp_path):
    # the retry is written from the input the override is already in; merging
    # again would set MaxIter back below what the recovery raised it to
    retry = OPT.replace("%scf maxiter 125 end", "%scf\n  maxiter 800\nend")
    out = _run(tmp_path, 'additions:initial=["%scf maxiter 500 end"]\n', "initial.retry2.inp", retry)

    assert out == retry


# ------------------------------------------------------------------ checking

def test_the_template_validates_as_before():
    assert not [e for e in validate_control_text(TEMPLATE) if "override" in e.lower()]
    assert get_orca_override_hints(TEMPLATE) == []


def test_a_name_no_job_has_is_a_hint_not_an_error():
    control = TEMPLATE + 'additions:inital=["%scf maxiter 500 end"]\n'

    assert validate_control_text(control) == validate_control_text(TEMPLATE)
    assert any("inital" in h for h in get_orca_override_hints(control))


def test_input_orca_cannot_read_is_an_error():
    errors = validate_control_text(TEMPLATE + 'additions:initial=["%scf maxiter 500"]\n')

    assert any("additions:initial" in e and "not closed" in e for e in errors)


def test_a_comma_inside_orca_input_is_part_of_it():
    # the CONTROL reader splits unbracketed values at commas; this one is ORCA's
    control = TEMPLATE + "additions:initial=%scf BrokenSym 1,1 end\n"

    assert validate_control_text(control) == validate_control_text(TEMPLATE)
    _, additions, _ = orca_overrides.parse_override_text(control)
    assert additions["initial"] == ["%scf BrokenSym 1,1 end"]


@pytest.mark.parametrize("name", ["initial", "ox_step_2", "red_step_1_OCCUPIER", "input7", "S0", "T3",
                                  "S1_TDDFT", "S1_second_deltaSCF", "S1_T1_ISC_msp1", "T1_S0_PHOSP_iroot2",
                                  "S0_IP", "E_n_cation", "XTB_GOAT", "all", "S*", "*_ISC*"])
def test_every_name_delfin_writes_is_known(name):
    hints = get_orca_override_hints(TEMPLATE + f"keyword:{name}=[DEFGRID3]\n")

    assert not [h for h in hints if "writes no job" in h]


def test_tddft_added_to_an_optimisation_is_named_for_what_it_does(tmp_path, caplog):
    # ORCA optimises root IROOT when %tddft is in an optimisation; the archive's
    # corrupted T1 retry was exactly that, by accident
    with caplog.at_level(logging.WARNING, logger="delfin.orca"):
        out = _run(tmp_path, 'additions:initial=["%tddft nroots 10 end"]\n', "initial.inp", OPT)

    assert "%tddft" in out
    assert "excited state" in caplog.text
    assert any("excited" in h for h in get_orca_override_hints(TEMPLATE + 'additions:initial=["%tddft nroots 10 end"]\n'))
    assert not any("excited" in h for h in get_orca_override_hints(TEMPLATE + 'additions:S1=["%tddft maxiter 300 end"]\n'))
