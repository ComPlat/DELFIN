"""An ORCA run whose MPI cannot open a socket says why and what works.

Input: an ORCA-driven step (orca_opt, xtb_opt through ORCA) run with more
than one process inside an agent's command sandbox, which allows no network
socket.

Before: OpenMPI stopped at once ("opal_ifinit: socket() failed with
errno=13", "No network interfaces were found for out-of-band
communications") and the step said only "ORCA calculation failed" / "xTB
optimization failed" -- the agent guessed at a missing xtb and a wrong core
count (seen 2026-11-14).
Now the step's error names the cause and the two ways that work (cores=1 in
the sandbox, or a job through the scheduler). Any other failure keeps its
old message. The output below is a hand-written stand-in, not ORCA output.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin import orca
from delfin.tools.adapters import orca as orca_adapter
from delfin.tools.adapters import xtb as xtb_adapter

XYZ = "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.757 0.587\nH 0.0 -0.757 0.587\n"

NO_SOCKET = (
    "[node:00007] opal_ifinit: socket() failed with errno=13\n"
    "No network interfaces were found for out-of-band communications. We require\n"
    "at least one available network for out-of-band messaging.\n"
    "ORCA finished by error termination in orca_util\n"
)
OTHER = "SCF NOT CONVERGED AFTER 125 CYCLES\nORCA finished by error termination in SCF\n"


@pytest.fixture
def fake_orca(monkeypatch):
    """A stand-in for the ORCA run: writes the given text as the output."""
    box = {"text": NO_SOCKET}

    def run(inp, out, timeout=None, working_dir=None, **kw):
        Path(out).write_text(box["text"])
        return False

    monkeypatch.setattr(orca, "run_orca", run)
    return box


def _geom(tmp_path):
    g = tmp_path / "in.xyz"
    g.write_text(XYZ)
    return g


def test_orca_opt_names_the_mpi_cause(tmp_path, fake_orca):
    wd = tmp_path / "opt"
    wd.mkdir()
    r = orca_adapter.OrcaOptAdapter().execute(
        wd, geometry=_geom(tmp_path), cores=2, charge=0, bang="! PBE0 def2-SVP")
    assert "MPI" in r.error and "cores=1" in r.error and "scheduler" in r.error


def test_xtb_opt_through_orca_names_it_too(tmp_path, fake_orca):
    wd = tmp_path / "xtb"
    wd.mkdir()
    r = xtb_adapter.XtbOptAdapter().execute(wd, geometry=_geom(tmp_path), cores=2, charge=0)
    assert "MPI" in r.error and "cores=1" in r.error


def test_another_failure_keeps_its_message(tmp_path, fake_orca):
    fake_orca["text"] = OTHER
    wd = tmp_path / "opt"
    wd.mkdir()
    r = orca_adapter.OrcaOptAdapter().execute(
        wd, geometry=_geom(tmp_path), cores=2, charge=0, bang="! PBE0 def2-SVP")
    assert r.error == "ORCA calculation failed"


def test_a_missing_output_is_no_note(tmp_path):
    assert orca.mpi_cannot_start_note(tmp_path / "nothing.out") is None
