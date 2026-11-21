"""FL1: the slow job asks the FF-free builder why it declined.

Input:  the slow-tests workflow and the trace flag
        ``delfin.manta.converter_backend`` reads.
Output: whether a decline in CI can be diagnosed from the log of the run
        that failed.
Semantics: ``_fffree_isomers`` declines a system by returning None and
        already NAMES the reason (``_scope_no`` emits
        ``[FFREE_SCOPE] <REASON>``), but only while the trace flag is set.
        Unset, a failing assertion can say nothing but
        ``assert isinstance(None, list)``.

Measured: six tests in three files failed in BOTH the 2026-11-20 and the
2026-11-21 slow runs and pass on a developer host -- 2 of 2 is not a
flake, it is an environment difference, and the reason the builder gave
was discarded. The flag writes to stderr, which pytest captures per test
and prints in the failure report.

The pair is pinned because this instrument has already been lost once:
``_ff_trace_write``'s own docstring records that it wrote to stderr while
the caller discarded stderr, so "the project's most important diagnostic
instrument wrote into the void for a week".

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_fl1_a_declined_build_names_its_reason.py -q
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SLOW_WF = REPO_ROOT / ".github" / "workflows" / "slow-tests.yml"
BACKEND = REPO_ROOT / "delfin" / "manta" / "converter_backend.py"

FLAG = "DELFIN_FFFREE_ISO_TRACE"


def _pytest_step_env() -> dict:
    """The env of the workflow step that runs pytest."""
    yaml = pytest.importorskip("yaml")
    spec = yaml.safe_load(SLOW_WF.read_text())
    for job in spec["jobs"].values():
        for step in job.get("steps") or []:
            if "pytest" in str(step.get("run") or ""):
                return dict(step.get("env") or {})
    raise AssertionError("no step in the slow workflow runs pytest")


def test_the_slow_job_turns_the_scope_trace_on():
    assert _pytest_step_env().get(FLAG) == "1", (
        "a decline in the slow job reports only 'assert isinstance(None, "
        "list)' while this is unset")


def test_the_flag_the_workflow_sets_is_the_flag_the_code_reads():
    """The pair, not either half. A rename on one side turns the instrument
    silently back off, which is how it was lost for a week before."""
    read = set(re.findall(r'environ\.get\(\s*"(DELFIN_FFFREE_[A-Z_]*TRACE)"',
                          BACKEND.read_text()))
    assert FLAG in read, sorted(read)
    assert FLAG in _pytest_step_env()


def test_the_reason_reaches_stderr_and_names_the_decline(capfd, monkeypatch):
    """Drives the real builder on a system it declines and asserts the
    reason is emitted -- the trace is the product's, not the test's."""
    monkeypatch.setenv(FLAG, "1")
    from delfin.manta import converter_backend

    assert converter_backend._fffree_isomers("c1ccccc1") is None
    err = capfd.readouterr().err
    assert "[FFREE_SCOPE]" in err, err[-400:]
    assert re.search(r"\[FFREE_SCOPE\] [A-Z_]{3,}", err), err[-400:]


def test_the_trace_does_not_change_what_the_builder_returns():
    """A diagnostic that alters the measurement is worse than none. The same
    input is built with the trace off and on and must agree."""
    from delfin.manta import converter_backend

    smi = "[Pt](Cl)(Cl)(N)N"
    os.environ.pop(FLAG, None)
    off = converter_backend._fffree_isomers(smi)
    os.environ[FLAG] = "1"
    try:
        on = converter_backend._fffree_isomers(smi)
    finally:
        os.environ.pop(FLAG, None)
    assert (off is None) == (on is None)
    if off is not None:
        assert len(off) == len(on)
        assert [lab for _x, lab in off] == [lab for _x, lab in on]
