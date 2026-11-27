"""A forced recalc never destroys a good output when the rerun fails.

``--recalc --occupier-override`` recomputes a job that already finished.  The
output it replaces is the record of that job.  Measured at the product entry
point before the fix: the good ``initial.out`` was deleted, the fake ORCA left
39 bytes of junk where it had been, and ``delfin.cli.main`` still returned 0.

The rerun is forced by ``_recalc_force_outputs``, which the recalc wrapper and
the engine's precomplete check both consult.  Deleting the ``.out`` therefore
never made the rerun happen -- it only made the loss possible.

Runs the real ``delfin.cli.main`` against a fixture workspace, with a stand-in
``orca`` on PATH: no ORCA is installed or executed.
"""
from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin import cli
from delfin.config import set_control_value
from delfin.define import TEMPLATE

GOOD_OUTPUT = ("* O R C A *\nProgram Version 6.1.1\n" + "SCF done\n" * 50
               + "                             ****ORCA TERMINATED NORMALLY****\n")

# A finished OCCUPIER folder: the entry the override moves away from, and one
# to override to.  The last two lines are what read_occupier_file parses.
OCCUPIER_TXT = """Method: PBE0 def2-SVP
Charge: 0
---------------------
FINAL SINGLE POINT ENERGY (1)   = -1.000 (H) <-- PREFERRED VALUE
multiplicity 1
--------------------------------------------------------------
FINAL SINGLE POINT ENERGY (2)   = -0.900 (H)
multiplicity 1, BrokenSym 1,1
--------------------------------------------------------------

LOWEST FINAL SINGLE POINT ENERGY: -1.000 (H)

(Selection: tolerance, APmethod 2, dev_max none, precision=3, epsilon=0.0005)
(Preferred Index: 1)
(Electron number: is_even)
"""


def _workspace(root: Path) -> Path:
    """A calc folder with one finished stage whose output an override forces out."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "input.txt").write_text("1\nH2\nH 0.0 0.0 0.0\n", encoding="utf-8")
    (root / "initial.xyz").write_text("1\nH2\nH 0.0 0.0 0.0\n", encoding="utf-8")

    control = TEMPLATE
    for key, value in (("charge", "0"), ("solvent", "water"), ("method", "OCCUPIER"),
                       ("implicit_solvation_model", "CPCM"), ("PAL", "1"),
                       ("maxcore", "1000"), ("calc_initial", "yes"),
                       ("oxidation_steps", ""), ("reduction_steps", ""), ("smiles", "")):
        control = set_control_value(control, key, value)
    (root / "CONTROL.txt").write_text(control, encoding="utf-8")

    occ = root / "initial_OCCUPIER"
    occ.mkdir(exist_ok=True)
    (occ / "OCCUPIER.txt").write_text(OCCUPIER_TXT, encoding="utf-8")
    (occ / "input.xyz").write_text("1\nH2\nH 0.0 0.0 0.0\n", encoding="utf-8")

    (root / "initial.inp").write_text("! PBE0 def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    (root / "initial.out").write_text(GOOD_OUTPUT, encoding="utf-8")
    return root


@pytest.fixture
def a_failing_orca(tmp_path, monkeypatch):
    """A stand-in ORCA that writes junk and exits 1 -- a failed rerun."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "orca"
    fake.write_text("#!/usr/bin/env bash\necho 'orca: no terminator here'\nexit 1\n",
                    encoding="utf-8")
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ["PATH"])
    # The resolver caches what it found; without this the stand-in is not seen.
    from delfin import qm_runtime
    qm_runtime.clear_resolver_cache()
    yield fake
    qm_runtime.clear_resolver_cache()


def test_a_failed_forced_rerun_keeps_the_output_of_the_finished_job(tmp_path, a_failing_orca):
    """Green with the fix; RED at origin/main, where the output is deleted first."""
    ws = _workspace(tmp_path / "ws")
    here = os.getcwd()
    os.chdir(ws)
    try:
        cli.main(["--recalc", "--occupier-override", "initial_OCCUPIER=2", str(ws)])
    finally:
        os.chdir(here)

    out = ws / "initial.out"
    assert out.is_file(), (
        "the finished output was destroyed by a rerun that failed; the previous "
        "output must be restored")
    assert out.read_text(encoding="utf-8") == GOOD_OUTPUT, (
        "what the failed rerun wrote must not become the job's output")
