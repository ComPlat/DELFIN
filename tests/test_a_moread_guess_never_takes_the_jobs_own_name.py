"""An ORCA guess file never carries the job's own name.

Input: the ORCA step adapters (orca_sp after orca_opt in a Pipeline), with
``moread`` pointing at the upstream step's ``calc.gbw``.

Before: the guess was copied into the new work dir as ``calc.gbw`` and the
input said ``%moinp "calc.gbw"``. The new job is also called ``calc``; ORCA
writes ``calc.gbw`` itself when it starts, so the guess was gone before it
was read: "No orbitals were found in the gbw file", error termination in
GUESS. Found by a live three-step chain (xtb_opt -> orca_opt -> orca_sp):
every orca_sp after an orca_opt in a Pipeline failed this way.
Now a guess whose name equals the job's name is copied in as
``<calc_name>_guess.gbw`` and the input reads that; any other name is kept.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.tools.adapters import orca as orca_adapter

XYZ = "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.757 0.587\nH 0.0 -0.757 0.587\n"


@pytest.fixture
def no_orca(monkeypatch):
    """A stand-in for the ORCA run: records the input it was given."""
    seen = {}

    def fake_run_orca(inp, out, timeout=None, working_dir=None):
        seen["inp"] = Path(inp).read_text()
        seen["dir"] = Path(working_dir)
        return False

    import delfin.orca
    monkeypatch.setattr(delfin.orca, "run_orca", fake_run_orca)
    return seen


def _run(tmp_path, moread_name, **kw):
    up = tmp_path / "01_orca_opt"
    up.mkdir()
    (up / moread_name).write_bytes(b"ORBITALS-OF-THE-UPSTREAM-STEP")
    geom = tmp_path / "in.xyz"
    geom.write_text(XYZ)
    wd = tmp_path / "02_orca_sp"
    wd.mkdir()
    orca_adapter.OrcaSpAdapter().execute(
        wd, geometry=geom, cores=1, charge=0, bang="! PBE0 def2-SVP",
        moread=str(up / moread_name), **kw)
    return wd


def _moinp(inp: str) -> str:
    m = re.search(r'%moinp\s+"([^"]+)"', inp)
    assert m, inp
    return m.group(1)


def test_a_guess_named_like_the_job_is_copied_under_another_name(tmp_path, no_orca):
    wd = _run(tmp_path, "calc.gbw")
    name = _moinp(no_orca["inp"])
    assert name != "calc.gbw"
    assert name == "calc_guess.gbw"
    assert (wd / name).read_bytes() == b"ORBITALS-OF-THE-UPSTREAM-STEP"
    assert "MOREAD" in no_orca["inp"]


def test_the_rule_follows_a_custom_job_name(tmp_path, no_orca):
    wd = _run(tmp_path, "sp.gbw", base_name="sp")
    assert _moinp(no_orca["inp"]) == "sp_guess.gbw"
    assert (wd / "sp_guess.gbw").is_file()


def test_a_guess_with_another_name_keeps_it(tmp_path, no_orca):
    wd = _run(tmp_path, "previous.gbw")
    assert _moinp(no_orca["inp"]) == "previous.gbw"
    assert (wd / "previous.gbw").is_file()
