"""The recovery termination flag must come from delfin.orca_read.

The task switches DELFIN's reading of normal termination onto
``delfin.orca_read`` (Quantity.TERMINATED_NORMALLY) wherever the outcome
provably does not change. This test pins the wiring: ``analyze_output``
decides "succeeded" from the TERMINATED_NORMALLY reading delivered by
``orca_read``, not from its own inline marker scan.

It is RED before the call site is switched (before the fix,
``analyze_output`` ignores ``orca_read`` and returns an error from its own
scan) and GREEN after. The value ``orca_read`` reports is chosen to
differ from what the fixture's own tail would say, so the test cannot
pass by accident: when ``orca_read`` says the run terminated normally,
``analyze_output`` must treat it as success; when ``orca_read`` says it
did not, the error path runs.
"""
from pathlib import Path

from delfin import orca_read
from delfin.orca_recovery import OrcaErrorType, OrcaErrorDetector

# ORCA 6.1.1 SCF that did not converge (the committed fixture's wording).
_SCF_NOT_CONVERGED = """
[file orca_leanscf/orca_leanscf.cpp, line 306]: Error (ORCA_LEANSCF): unfortunately,
the SCF has not converged. There may be a way out but we have to stop here
ORCA finished by error termination in LEANSCF
  .... aborting the run
"""


def _out(tmp_path: Path, body: str, name: str = "job.out") -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


def test_analyze_output_takes_termination_from_orca_read(tmp_path, monkeypatch):
    """The no-error answer follows orca_read's TERMINATED_NORMALLY.

    Before the switch, analyze_output scans the fixture's own tail and
    returns LEANSCF_NOT_CONVERGED. The patch forces orca_read to report
    normal termination; after the switch the run is therefore a success.
    """
    out = _out(tmp_path, _SCF_NOT_CONVERGED)

    def fake_read(work_dir, base_name, quantity):
        assert quantity is orca_read.Quantity.TERMINATED_NORMALLY
        return orca_read.Reading(True, "", "text", "test")

    monkeypatch.setattr(orca_read, "read", fake_read)
    assert OrcaErrorDetector.analyze_output(out) is None


def test_analyze_output_falls_through_when_orca_read_says_not_terminated(
    tmp_path, monkeypatch
):
    """When orca_read reports no normal termination the error path runs."""
    out = _out(tmp_path, _SCF_NOT_CONVERGED)

    def fake_read(work_dir, base_name, quantity):
        assert quantity is orca_read.Quantity.TERMINATED_NORMALLY
        return orca_read.Reading(False, "", "text", "test")

    monkeypatch.setattr(orca_read, "read", fake_read)
    assert OrcaErrorDetector.analyze_output(out) is OrcaErrorType.LEANSCF_NOT_CONVERGED
