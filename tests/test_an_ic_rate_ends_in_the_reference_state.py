"""ICs are the ones ORCA's ESD(IC) can compute, each with its own root.

ORCA 6.1.1 manual 5.5.5: an ESD(IC) rate runs from root IRoot into the
reference state, and "the final state is always S0 and this cannot be
changed" -- "iroot 1 # Change to 2 for S2-S0 IC rate".  DELFIN had it the
other way round: the validator accepted only transitions into S1 or T1, so
S1>S0 was refused and the template's S2>S1 was run -- as the S2->S0 coupling
on the S1 geometry with the S1 Hessian as the ground state.  And the pipeline
defaulted a global IROOT to 1, which every IC input read before its own
transition, so S2>S0 would have been computed as S1>S0.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.common.control_validator import unsupported_ic_reason
from delfin.config import set_control_value, validate_control_text
from delfin.define import TEMPLATE

_XYZ = "3\nwater\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.467\nH 0.000 -0.757 -0.467\n"


def _control(**keys) -> str:
    text = TEMPLATE
    base = {"charge": "0", "solvent": "water", "method": "classic",
            "ESD_modus": "TDDFT", "ESD_T1_opt": "uks"}
    for key, value in {**base, **keys}.items():
        text = set_control_value(text, key, value)
    return text


@pytest.mark.parametrize("transition", ["S1>S0", "S2>S0", "S5>S0", "T2>T1", "T3>T1"])
def test_orca_can_compute(transition):
    assert unsupported_ic_reason(transition) is None


@pytest.mark.parametrize("transition, fragment", [
    ("S2>S1", "did you mean S2>S0"),
    ("S3>S2", "did you mean S3>S0"),
    ("T3>T2", "did you mean T3>T1"),
    ("S1>T1", "ISC"),
    ("S0>S0", "does not start above"),
    ("T1>T1", "does not start above"),
    ("S2-S0", "not a transition"),
])
def test_orca_cannot_compute(transition, fragment):
    assert fragment in unsupported_ic_reason(transition)


def test_the_validator_refuses_what_orca_cannot_do_when_esd_is_on():
    errors = validate_control_text(_control(ESD_modul="yes", ICs="[S2>S1]"))
    assert any("did you mean S2>S0" in e for e in errors), errors
    assert validate_control_text(_control(ESD_modul="yes", ICs="[S1>S0,S2>S0,T2>T1]")) == []


def test_an_old_file_with_esd_off_still_reads():
    # the old template shipped ICs=[S2>S1]; with ESD off it is inert
    assert validate_control_text(_control(ESD_modul="no", ICs="[S2>S1]")) == []


def test_the_template_asks_for_an_ic_orca_can_compute():
    line = re.search(r"(?m)^ICs=\[(.*)\]$", TEMPLATE).group(1)
    for transition in line.split(","):
        assert unsupported_ic_reason(transition) is None, transition


def _ic_input(tmp_path: Path, pair: str, **config) -> str:
    from delfin.esd_input_generator import create_ic_input

    esd = tmp_path / "ESD"
    esd.mkdir(exist_ok=True)
    for stem in ("S0", "S1", "S2", "S3", "T1", "T2", "T3"):
        (esd / f"{stem}.xyz").write_text(_XYZ)
    path = create_ic_input(pair, esd, 0, "water", [], "def2-SVP", "def2-TZVP",
                           {"TDDFT_nroots": 15, **config})
    return Path(path).read_text()


@pytest.mark.parametrize("pair, iroot, mult, geometry_from, gs_hess, es_hess", [
    ("S1>S0", 1, 1, "S0", "S0.hess", "S1.hess"),
    ("S2>S0", 2, 1, "S0", "S0.hess", "S2.hess"),
    ("T2>T1", 1, 3, "T1", "T1.hess", "T2.hess"),
    ("T3>T1", 2, 3, "T1", "T1.hess", "T3.hess"),
])
def test_each_ic_input_names_its_own_root(tmp_path, pair, iroot, mult, geometry_from, gs_hess, es_hess):
    # IROOT=1 is what the pipeline used to put into every config; it must not win
    text = _ic_input(tmp_path, pair, IROOT="1")
    assert re.search(rf"(?m)^\s*iroot {iroot}$", text), text
    assert f"* xyz 0 {mult}" in text
    assert f'GSHESSIAN       "{gs_hess}"' in text
    assert f'ESHESSIAN       "{es_hess}"' in text
    assert "ESD(IC)" in text and "nacme TRUE" in text


def test_the_scheduler_runs_what_the_validator_accepts(tmp_path, monkeypatch):
    import delfin.esd_module as esd_module

    class _Manager:
        total_cores = 8
        max_jobs = 2

        def __init__(self):
            self.jobs = []
            self._completed = set()

        def add_job(self, job):
            self.jobs.append(job)

    manager = _Manager()
    esd_module._populate_ic_jobs(manager, ["S1>S0", "S2>S0", "T2>T1", "S2>S1"], tmp_path,
                                 0, "water", [], "def2-SVP", "def2-TZVP", {})
    scheduled = {job.job_id: job.dependencies for job in manager.jobs}
    assert scheduled == {
        "esd_ic_S1_S0": {"esd_S1", "esd_S0"},
        "esd_ic_S2_S0": {"esd_S2", "esd_S0"},
        "esd_ic_T2_T1": {"esd_T2", "esd_T1"},
    }
