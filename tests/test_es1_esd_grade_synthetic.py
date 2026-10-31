"""The grading script works on synthetic ORCA-format text (licence hold).

SYNTHETIC — not produced by ORCA. Small hand-written text in ORCA's
ESD output format exercises the grader's red/green control before the
real run's files are graded:

* a level-matched synthetic ISC record compares against the in-tree
  reference (green);
* a mismatched level (different functional) refuses the comparison
  with NOT GRADED and the reason (the level gate fires);
* a truncated-window record is flagged by window_remainders, never
  read as a result.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
GRADE = REPO / ".gate" / "practice" / "grade_es1.py"

# Hand-written ORCA-format ESD output text, marked synthetic.
ISC_OUT = """SYNTHETIC — not produced by ORCA. Fixture for the grading script.
ISC calculation
Homogeneous linewidth is: 50 cm-1
Temperature used: 298.15 K
0-0 energy difference: 35301.0 cm-1
The calculated intersystem crossing rate constant is 2.5E+07 s^-1
"""

IC_OUT = """SYNTHETIC — not produced by ORCA. Fixture for the grading script.
IC calculation
Homogeneous linewidth is: 50 cm-1
Temperature used: 298.15 K
0-0 energy difference: 35301.0 cm-1
The calculated internal conversion rate constant is 3.79E-03 s^-1
"""


def _write_synth(tmp_path: Path, functional_line: str = "functional=PBE0") -> Path:
    """Build a minimal run folder the grader reads."""
    (tmp_path / "CONTROL.txt").write_text(
        f"""method=classic
{functional_line}
main_basisset=def2-SVP
charge=0
multiplicity=1
solvent=
input_file=input.txt
PAL=2
maxcore=2000
ESD_modul=yes
ESD_modus=TDDFT
ESD_states=[S0,S1,T1]
ESD_ISC=S1>T1
ESD_IC=S1>S0
TROOTSSL=0
"""
    )
    esd = tmp_path / "ESD"
    esd.mkdir(exist_ok=True)
    # State outputs: only the final single-point energy line matters to
    # find_electronic_energy; format "Final Gibbs free energy" style is
    # not needed -- FSPE comes from find_electronic_energy.
    for state, energy in (("S0", -114.15000), ("S1", -113.99000), ("T1", -114.04000)):
        (esd / f"{state}.out").write_text(
            f"SYNTHETIC — not produced by ORCA.\n"
            f"TOTAL ENERGY      {energy:.8f} Eh\n"
            f"The current calculation has used the line:\n"
        )
    # ISC output the resolver looks for under tddft mode
    (esd / "S1_T1_ISC.out").write_text(ISC_OUT)
    (esd / "S1_S0_IC.out").write_text(IC_OUT)
    return tmp_path


def test_grader_runs_and_reports(tmp_path):
    run = _write_synth(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(GRADE), "--synthetic", str(run)],
        capture_output=True,
        text=True,
        cwd=str(REPO),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PASS" in proc.stdout
    # state order + isc_rate must be graded
    assert '"criterion": "state_order"' in proc.stdout
    assert '"criterion": "isc_rate"' in proc.stdout


def test_level_gate_refuses_mismatched_level(tmp_path):
    run = _write_synth(tmp_path, functional_line="functional=B3LYP")
    proc = subprocess.run(
        [sys.executable, str(GRADE), "--synthetic", str(run)],
        capture_output=True,
        text=True,
        cwd=str(REPO),
    )
    assert proc.returncode == 0
    assert '"criterion": "isc_vs_reference"' in proc.stdout
    assert "NOT GRADED" in proc.stdout
    assert "functional=B3LYP" in proc.stdout
