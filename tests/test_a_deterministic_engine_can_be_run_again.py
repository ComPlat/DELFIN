"""Replay means re-execution, and for these engines it is available.

``run_replay`` re-derives a run's metadata -- commit, prompt digest,
engine versions, contract flags -- and compares it. That is worth
having, and it is not what this field calls replay: ReproZip, Nextflow
``-resume``, Snakemake and AiiDA all re-run the work. A reviewer who
reads "replayable" and finds no engine was ever started has caught the
claim overstating itself. The command is now ``verify``, with ``replay``
kept as an alias, and THIS module is the replay.

ORCA and xtb are deterministic for a fixed binary and a fixed input,
which is the advantage a language-model system does not have: a model's
output cannot be re-derived and a calculation's can. So the strongest
check available is the plain one -- run the input again, compare the
energy.

WHAT IS ASSERTED:

  * **the original workspace is never written to.** A verification that
    overwrites the evidence destroys what it was asked about. Checked by
    hashing the directory before and after.
  * **a timeout is not a disagreement**, and neither is a missing
    engine. Both are UNAVAILABLE.
  * **a version difference is reported and left unjudged.** A newer
    program giving a different number is not a reproducibility failure
    of the run, and calling it one would make every upgrade look like a
    defect.
  * **nothing compared is not a pass.** ``reproduced()`` is None, never
    True, when no step could be compared.

The engine here is a stub script the test writes, so these run on a
machine with no quantum chemistry program installed.
"""
from __future__ import annotations

import hashlib
import stat
from pathlib import Path

import pytest

from delfin.agent import run_rerun

_E1 = -76.12345678
_E2 = -76.12349999        # 4.3e-5 Ha away: outside the default tolerance


def _stub(tmp_path, *, energy=_E1, version="6.1.0",
          emit_energy=True) -> Path:
    """An executable that writes an ORCA-shaped output beside its input.

    A stub and not a mock: ``rerun_step`` starts a process, copies files
    and reads what lands on disk, and a mock would test none of that.
    """
    body = [
        "#!/bin/sh",
        'name="${1%.inp}"',
    ]
    body.append(f'echo "           Program Version {version}  -  RELEASE  -"'
                f' > "$name.out"')
    if emit_energy:
        body.append(f'echo "FINAL SINGLE POINT ENERGY   {energy:.8f}"'
                    f' >> "$name.out"')
    body.append('echo "                   ****ORCA TERMINATED NORMALLY****"'
                ' >> "$name.out"')
    path = tmp_path / "stub_engine.sh"
    path.write_text("\n".join(body) + "\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP)
    return path


def _workspace(tmp_path, *, energy=_E1, version="6.1.0") -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "job.inp").write_text("! PBE0 def2-SVP\n* xyzfile 0 1 mol.xyz\n",
                                encoding="utf-8")
    (ws / "mol.xyz").write_text("1\n\nH 0.0 0.0 0.0\n", encoding="utf-8")
    (ws / "job.out").write_text(
        f"           Program Version {version}  -  RELEASE  -\n"
        f"FINAL SINGLE POINT ENERGY   {energy:.8f}\n"
        "                   ****ORCA TERMINATED NORMALLY****\n",
        encoding="utf-8")
    return ws


def _fingerprint(root: Path) -> str:
    h = hashlib.sha256()
    for p in sorted(root.rglob("*")):
        h.update(str(p.relative_to(root)).encode())
        if p.is_file():
            h.update(p.read_bytes())
    return h.hexdigest()


# --- the comparison ----------------------------------------------------

def test_the_same_input_giving_the_same_energy_reproduces(tmp_path):
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path)))
    v = rep.verdicts[0]
    assert v.verdict == run_rerun.REPRODUCED, v.note
    assert v.rerun == pytest.approx(_E1, abs=1e-9)
    assert v.deviation == pytest.approx(0.0, abs=1e-9)
    assert rep.reproduced() is True


def test_a_different_energy_from_the_same_input_is_a_difference(tmp_path):
    """For a deterministic engine this is a real finding, not noise."""
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path, energy=_E2)))
    v = rep.verdicts[0]
    assert v.verdict == run_rerun.DIFFERS
    assert v.deviation > v.tolerance
    assert rep.reproduced() is False


def test_the_tolerance_is_stated_in_every_verdict(tmp_path):
    """A reader must never have to look up what counted as agreement."""
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path)),
                          tolerance=1e-3)
    assert rep.verdicts[0].tolerance == 1e-3


def test_a_loose_enough_tolerance_accepts_the_difference(tmp_path):
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path, energy=_E2)),
                          tolerance=1e-3)
    assert rep.verdicts[0].verdict == run_rerun.REPRODUCED


# --- the original is evidence and stays untouched ----------------------

def test_the_workspace_is_not_written_to(tmp_path):
    """A verification that overwrites what it verifies destroys the
    thing it was asked about."""
    ws = _workspace(tmp_path)
    before = _fingerprint(ws)
    run_rerun.rerun(ws, executable=str(_stub(tmp_path)))
    assert _fingerprint(ws) == before


def test_a_stale_output_cannot_be_read_as_the_rerun_result(tmp_path):
    """The previous .out is not copied into the scratch directory: if it
    were, an engine that produced nothing would look like agreement."""
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(
        ws, executable=str(_stub(tmp_path, emit_energy=False)))
    assert rep.verdicts[0].verdict == run_rerun.UNAVAILABLE
    assert "no final energy" in rep.verdicts[0].note


# --- what is not a disagreement ----------------------------------------

def test_no_engine_is_unavailable_and_not_a_pass(tmp_path):
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(tmp_path / "does_not_exist"))
    v = rep.verdicts[0]
    assert v.verdict == run_rerun.UNAVAILABLE
    assert rep.reproduced() is None


def test_a_timeout_is_unavailable_and_not_a_difference(tmp_path,
                                                      monkeypatch):
    """A calculation that did not finish in the time allowed has not
    disagreed with anything.

    The timeout is raised rather than waited for: the branch under test
    is OUR handling of it, and a test that sleeps a real child process
    would need a shell the machine may not have. A skip runs nowhere,
    which this repository's machine ratchet refuses for exactly that
    reason.
    """
    import subprocess as _sp

    def _timeout(*_a, **_k):
        raise _sp.TimeoutExpired(cmd="stub", timeout=0.5)

    monkeypatch.setattr(run_rerun.subprocess, "run", _timeout)
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path)),
                          timeout_s=0.5)
    v = rep.verdicts[0]
    assert v.verdict == run_rerun.UNAVAILABLE
    assert "did not finish" in v.note
    assert v.seconds == pytest.approx(0.5)
    assert rep.reproduced() is None


def test_a_version_difference_is_reported_and_left_unjudged(tmp_path):
    """A newer program giving a different number is not a
    reproducibility failure of the RUN."""
    ws = _workspace(tmp_path, version="6.0.1")
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path, energy=_E2,
                                                   version="6.1.0")))
    v = rep.verdicts[0]
    assert v.verdict == run_rerun.VERSION
    assert v.recorded_version == "6.0.1" and v.rerun_version == "6.1.0"
    assert v.recorded is not None and v.rerun is not None, (
        "both numbers must be reported even though no verdict is claimed")
    assert rep.reproduced() is None


def test_an_original_without_an_energy_has_nothing_to_compare(tmp_path):
    ws = _workspace(tmp_path)
    (ws / "job.out").write_text("           Program Version 6.1.0\n",
                                encoding="utf-8")
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path)))
    assert rep.verdicts[0].verdict == run_rerun.UNAVAILABLE
    assert "nothing to compare" in rep.verdicts[0].note


# --- choosing what to run ----------------------------------------------

def test_a_step_can_be_named(tmp_path):
    ws = _workspace(tmp_path)
    (ws / "other.inp").write_text("! HF\n", encoding="utf-8")
    (ws / "other.out").write_text(
        "FINAL SINGLE POINT ENERGY   -1.00000000\n", encoding="utf-8")
    rep = run_rerun.rerun(ws, only="other", executable=str(_stub(tmp_path)))
    assert len(rep.verdicts) == 1
    assert rep.verdicts[0].step == "other"


def test_a_cap_reports_what_it_left_out(tmp_path):
    """A caller that asked for a check must not get hours of compute by
    surprise, and must be told what was skipped."""
    ws = _workspace(tmp_path)
    for name in ("b", "c"):
        (ws / f"{name}.inp").write_text("! HF\n", encoding="utf-8")
        (ws / f"{name}.out").write_text(
            "FINAL SINGLE POINT ENERGY   -1.00000000\n", encoding="utf-8")
    rep = run_rerun.rerun(ws, max_steps=1, executable=str(_stub(tmp_path)))
    assert len(rep.verdicts) == 1
    assert any("3 steps found and 1 run" in n for n in rep.notes)


def test_an_input_with_no_output_is_not_a_step(tmp_path):
    """It was never run here, so there is nothing to reproduce."""
    ws = _workspace(tmp_path)
    (ws / "planned.inp").write_text("! HF\n", encoding="utf-8")
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path)))
    assert [v.step for v in rep.verdicts] == ["job"]


def test_an_empty_workspace_compares_nothing_and_says_so(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    rep = run_rerun.rerun(empty, executable=str(_stub(tmp_path)))
    assert rep.reproduced() is None
    assert any("nothing can be run again" in n for n in rep.notes)


# --- the energy parser -------------------------------------------------

def test_the_last_energy_wins(tmp_path):
    """An optimisation prints one per cycle; the result is the last."""
    text = ("FINAL SINGLE POINT ENERGY   -1.00000000\n"
            "FINAL SINGLE POINT ENERGY   -2.00000000\n")
    assert run_rerun.final_energy(text) == pytest.approx(-2.0)


def test_no_energy_is_none_rather_than_zero(tmp_path):
    assert run_rerun.final_energy("nothing here") is None
    assert run_rerun.final_energy("") is None


def test_the_report_is_json_serialisable(tmp_path):
    import json
    ws = _workspace(tmp_path)
    rep = run_rerun.rerun(ws, executable=str(_stub(tmp_path)))
    json.dumps(rep.as_dict())
    assert rep.as_dict()["counts"][run_rerun.REPRODUCED] == 1
