"""The OPI gate asks for a file OPI reads, and reads only the half it needs.

Two defects measured on this tree (2026-11-26), both through the reader the
pipeline and the reports call, ``delfin.orca_read.read``:

1. **A gate step on a file OPI never opens.** ``_opi_ready`` refused a run
   without ``<base>.property.txt``. OPI's property reader consumes only
   ``<base>.property.json``; ``.property.txt`` appears in OPI solely as the
   *input* of ``orca_2json`` inside
   ``opi.execution.core.Runner.create_property_json``, which DELFIN's reader
   never calls. Instrumenting ``open`` over a full ``Output.parse()`` counted
   **0** opens of ``.property.txt``. So a run with a valid property JSON was
   refused OPI for a file nobody would read: ``read`` answered ``text`` (or
   ``none``) with the note "property txt missing".

2. **A parse that shells out for a ``.gbw``.** ``Output.parse()`` defaults to
   ``read_gbw_json=True`` and then runs ``create_missing_gbw_json``, i.e.
   ``orca_2json``, whenever a ``<base>.gbw`` sits beside the output. Every
   real DELFIN calculation folder keeps its ``.gbw``, so on any host where
   OPI cannot find ``orca_2json`` the parse raised ``FileNotFoundError`` and
   the reader silently fell back to the text parser — from a call that only
   wanted the property JSON, which was there all along.

Both tests fail on the tree before the fix, each for the defect it names:
test A because the gate refuses, test B because the parse raises. The control
test passes on both trees (it satisfies the old gate and has no ``.gbw``), so
a red result cannot be mistaken for the fix.

These tests are machine-independent by construction: the ORCA stand-in they
need is written into ``tmp_path`` and pointed at through ``OPI_ORCA``, which
is the first place OPI looks (``opi.execution.core.Runner.set_orca_path``).
Nothing here asks what the host has installed.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

pytest.importorskip("opi", reason="OPI (orca-pi) not installed")

from delfin import orca_read  # noqa: E402
from delfin.orca_read import Quantity, read  # noqa: E402

TOL = 1e-6

#: The values the property JSON and the .out both carry (hand-written
#: fixture data, the same triple the committed synth_h2o fixture holds).
FINAL_E = -76.32126902150691
GIBBS = -76.3177312027805
ZPE = 0.021202376295716842

OUT_TEXT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output
                        Program Version 6.1.1  -  RELEASE   -

               *           SCF CONVERGED AFTER  12 CYCLES          *

                       ------------------
                       THERMOCHEMISTRY AT 298.15 K
                       ------------------
Temperature         ... 298.15 K
Zero point energy               ...   0.02120238 Eh    13.30 kcal/mol
Total thermal energy            ...  -76.29723134 Eh
Total Enthalpy                  ...  -76.29628713 Eh
Final entropy term              ...    0.02144407 Eh    13.46 kcal/mol
Final Gibbs free energy         ...  -76.31773120 Eh
       FINAL SINGLE POINT ENERGY     -76.32126902150691
****ORCA TERMINATED NORMALLY****
"""

PROPERTY_JSON = {
    "Calculation_Info": {"Charge": 0, "Mult": 1, "progName": "ORCA",
                        "version": "6.1.1"},
    "Calculation_Status": {"Status": "NORMAL TERMINATION"},
    "Geometries": [{
        "Single_Point_Data": {"Converged": True, "FinalEnergy": FINAL_E},
        "thermochemistry_energies": [{
            "zpe": ZPE,
            "innerenergyu": -76.29723134165184,
            "enthalpyh": -76.29628713260958,
            "entropys": 0.0214440701709188,
            "freeenergyg": GIBBS,
        }],
    }],
}


def _run(tmp_path: Path, *, name: str = "calc") -> Path:
    """A run folder with the output, the property JSON and the input."""
    d = tmp_path / name
    d.mkdir()
    (d / f"{name}.out").write_text(OUT_TEXT)
    (d / f"{name}.property.json").write_text(json.dumps(PROPERTY_JSON))
    (d / f"{name}.inp").write_text("! PBE0 def2-SVP\n\n* xyz 0 1\n"
                                   "O 0.0 0.0 0.117\nH 0.0 0.757 -0.469\n"
                                   "H 0.0 -0.757 -0.469\n*\n")
    return d


@pytest.fixture()
def orca_stand_in(tmp_path, monkeypatch):
    """An ORCA ``bin`` folder that answers, and records being asked.

    ``OPI_ORCA`` is the first place OPI looks for ORCA, so this replaces the
    host's installation rather than depending on one. ``orca_2json`` writes a
    trace, which is how a test proves it was consulted.
    """
    bin_dir = tmp_path / "orca-stand-in"
    bin_dir.mkdir()
    trace = tmp_path / "orca_2json.calls"
    stand_in = bin_dir / "orca_2json"
    stand_in.write_text(
        "#!/bin/sh\n"
        f'echo "$@" >> "{trace}"\n'
        "exit 1\n")
    stand_in.chmod(0o755)
    (bin_dir / "orca").write_text("#!/bin/sh\nexit 0\n")
    (bin_dir / "orca").chmod(0o755)
    monkeypatch.setenv("OPI_ORCA", str(bin_dir))
    return trace


def test_a_property_json_alone_is_enough_for_opi(tmp_path, orca_stand_in):
    """Test A: no ``.property.txt`` -- OPI still reads the run.

    Fails on the unfixed tree: the gate refuses with "property txt missing"
    and ``read`` answers from the text parser although OPI answers correctly.
    """
    d = _run(tmp_path)
    assert not (d / "calc.property.txt").exists()

    reading = read(d, "calc", Quantity.GIBBS)

    assert reading.source == "opi", reading.note
    assert reading.value == pytest.approx(GIBBS, abs=TOL)
    assert "property txt" not in reading.note.lower()


def test_a_gbw_beside_the_output_does_not_make_opi_shell_out(
        tmp_path, orca_stand_in):
    """Test B: a ``.gbw`` must not route the read through ``orca_2json``.

    Fails on the unfixed tree: ``parse()`` calls ``orca_2json`` because of the
    ``.gbw``, the stand-in exits 1, and the read falls back to text. The
    trace assertion then has a file to read, i.e. it can fail — which is the
    point of it.
    """
    d = _run(tmp_path)
    (d / "calc.property.txt").write_text("ORCA property file\n")
    (d / "calc.gbw").write_bytes(b"#GBW\0stand-in for a real gbw\n")

    reading = read(d, "calc", Quantity.FINAL_SINGLE_POINT_ENERGY)

    assert not orca_stand_in.exists(), (
        "the reader ran orca_2json for a .gbw; it only needs the property JSON")
    assert reading.source == "opi", reading.note
    assert reading.value == pytest.approx(FINAL_E, abs=TOL)


def test_control_a_run_the_old_gate_accepted_still_reads_through_opi(tmp_path):
    """Control: the shape the old gate accepted reads identically.

    Passes on both trees — it carries the ``.property.txt`` the old gate
    wanted and no ``.gbw`` — so a red result above is the fix's own.
    """
    d = _run(tmp_path, name="plain")
    (d / "plain.property.txt").write_text("ORCA property file\n")

    reading = read(d, "plain", Quantity.ZPE)

    assert reading.source == "opi", reading.note
    assert reading.value == pytest.approx(ZPE, abs=TOL)


def test_a_missing_property_json_still_falls_back_to_text(tmp_path):
    """The gate step that IS about OPI stays: no JSON, no OPI.

    Without the property JSON OPI has nothing to read, and ``read`` must
    answer from the text parser with the text value — unchanged behaviour.
    """
    d = _run(tmp_path, name="nojson")
    (d / "nojson.property.json").unlink()

    reading = read(d, "nojson", Quantity.GIBBS)

    assert reading.source == "text"
    assert reading.value == pytest.approx(-76.3177312, abs=TOL)
    assert "property json" in reading.note.lower()
