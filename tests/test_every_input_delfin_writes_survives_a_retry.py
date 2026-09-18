"""Every input DELFIN writes can be taken apart and written back by the recovery.

The recovery rewrites a failed job's input.  Measured on 543 archived inputs
with every strategy and attempt, 7 519 of 14 350 retries were wrong before the
writer was rebuilt, and ORCA refused them: a one-line ``%scf`` came back
nested, a metal's ``NewGTO`` was dropped when the coordinates were updated, a
second job's ``%base`` ended up on the first.

That audit needed the archive.  This one needs nothing: it takes the inputs
DELFIN writes today (orca_input_corpus) and puts every strategy of the
recovery to each of them, then asks of the result what ORCA would ask.  A
writer that starts writing something the recovery cannot take apart fails
here, in the commit that changes it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.common import orca_input as oi
from delfin.orca_recovery import OrcaErrorType as E
from delfin.orca_recovery import OrcaInputModifier, RecoveryStrategy
from orca_input_corpus import build_corpus

#: Every error the recovery writes a retry for, with the attempts it escalates
#: through.  INPUT_ERROR is not here: ORCA refused the input itself, and no
#: retry is written for it (that is its own test).
PLAN = [(error, attempt)
        for error, attempts in (
            (E.SCF_NO_CONVERGENCE, (1, 2, 3)),
            (E.LEANSCF_NOT_CONVERGED, (1, 2, 3, 4)),
            (E.TRAH_SEGFAULT, (1, 2)),
            (E.DIIS_ERROR, (1, 2, 3)),
            (E.GEOMETRY_NOT_CONVERGED, (1, 2, 3)),
            (E.MPI_CRASH, (1, 2, 3)),
            (E.CIS_FAILURE, (1, 2, 3)),
            (E.TDDFT_ROOT_COLLAPSE, (1, 2, 3)),
            (E.FREQUENCY_FAILURE, (1, 2)),
            (E.MEMORY_ERROR, (1, 2)),
            (E.ESD_RATE_UNPHYSICAL, (1, 2)),
            (E.ESD_WINDOW_TRUNCATED, (1,)),   # handing the window back is the whole fix
            (E.TRANSIENT_SYSTEM_ERROR, (1, 2)),
        )
        for attempt in attempts]

_ATOM = re.compile(r"^\s*([A-Z][a-z]?)\s+-?\d")


@pytest.fixture(scope="module")
def corpus(tmp_path_factory) -> dict:
    return build_corpus(tmp_path_factory.mktemp("corpus"))


def _atoms(job_text: str) -> list:
    return [m.group(1) for line in job_text.splitlines() if (m := _ATOM.match(line))]


def _scf_values(job_text: str) -> dict:
    return {st.key: " ".join(st.lines[0].split()[1:]) for item in oi.parse_job(job_text)
            if item.kind == "block" and item.name == "scf" for st in item.block.statements if st.key}


def _retry(tmp_path: Path, name: str, text: str, error, attempt: int) -> str:
    job = tmp_path / Path(name).name
    job.parent.mkdir(parents=True, exist_ok=True)
    job.write_text(text, encoding="utf-8")
    out = job.with_suffix(".out")
    out.write_text("\n $$$$$$$$$$$$$$$$  JOB NUMBER  1 $$$$$$$$$$$$$$\n", encoding="utf-8")
    strategy = RecoveryStrategy(error, attempt, {})
    strategy.output_file = out
    written = OrcaInputModifier(job, {}).apply_recovery(strategy)
    assert written != job, f"no retry was written for {name} after {error.value}"
    return written.read_text(encoding="utf-8")


@pytest.mark.parametrize("error, attempt", PLAN, ids=lambda v: getattr(v, "value", v))
def test_every_strategy_writes_an_input_orca_can_read_for_every_job(tmp_path, corpus, error, attempt):
    for name, text in corpus.items():
        retry = _retry(tmp_path / f"{error.value}_{attempt}", name, text, error, attempt)
        original_jobs = oi.split_jobs(text)
        jobs = oi.split_jobs(retry)

        assert len(jobs) == len(original_jobs), f"{name}: a job was lost or added"
        for job, original in zip(jobs, original_jobs):
            items = oi.parse_job(job)              # raises on anything ORCA could not read
            own = oi.setting_value(oi.parse_job(original), "base")
            assert oi.setting_value(items, "base") == (own or f'"{Path(name).stem}"'), \
                f"{name}: the retry's products would be written under another name"
            assert oi.has_block(items, "tddft") == oi.has_block(oi.parse_job(original), "tddft"), name
            assert _atoms(job) == _atoms(original), f"{name}: the geometry changed"
            for key, value in _scf_values(job).items():
                assert not key.startswith("sub "), f"{name}: a bare word in %scf opens a sub-block"
                if key.lower() == "convergence":
                    assert value.split()[0] in ("Tight", "VeryTight", "Strong", "Medium"), \
                        f"{name}: {value} is not an ORCA convergence value"


@pytest.mark.parametrize("error, attempt", [(E.SCF_NO_CONVERGENCE, 2), (E.GEOMETRY_NOT_CONVERGED, 1)],
                         ids=("scf", "geometry"))
def test_the_metals_own_basis_survives_every_retry(tmp_path, corpus, error, attempt):
    for name in ("red_step_1.inp", "red_step_1_OCCUPIER/input.inp", "red_step_1_OCCUPIER/input3.inp"):
        retry = _retry(tmp_path / name.replace("/", "_"), name, corpus[name], error, attempt)

        metal = [line for line in retry.splitlines() if line.strip().startswith("Co ")]
        assert metal and 'NewGTO "def2-TZVP" end' in metal[0], name


def test_broken_symmetry_is_still_there_after_a_retry(tmp_path, corpus):
    retry = _retry(tmp_path, "input3.inp", corpus["red_step_1_OCCUPIER/input3.inp"],
                   E.SCF_NO_CONVERGENCE, 3)

    assert "BrokenSym 3,2" in retry
    assert "APMethod 2" in retry


def test_a_rate_job_whose_window_was_cut_gets_orcas_own_window_back(tmp_path, corpus):
    retry = _retry(tmp_path, "S1_T1_ISC_ms0.inp", corpus["ESD/S1_T1_ISC_ms0.inp"],
                   E.ESD_WINDOW_TRUNCATED, 1)
    esd = [item for item in oi.parse_job(retry) if item.kind == "block" and item.name == "esd"]

    assert esd, "the ESD block is gone"
    assert oi.block_value(oi.parse_job(retry), "esd", "maxtime") is None
    assert oi.block_value(oi.parse_job(retry), "esd", "npoints") is None
