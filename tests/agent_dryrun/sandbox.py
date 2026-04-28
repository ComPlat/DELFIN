"""Filesystem sandbox for agent dry-run tests.

Goal: the agent never touches the real DELFIN repo or the user's calc/
directories. Each test gets a freshly built temporary tree containing:

- ``calc/`` with a few fake calculation folders (BP86 / PBE0 / B3LYP)
  containing synthetic ORCA outputs that the parsers can consume.
- ``archive/`` (empty, writable target for move_to_archive).
- ``agent_workspace/`` (writable scratch for plots/scripts).
- A minimal CONTROL.txt at the root.

The runner is invoked with ``cwd=<sandbox>`` and ``--add-dir`` set to
the sandbox only. ``--permission-mode plan`` already blocks mutations,
the sandbox stops accidental reads of sensitive paths AND gives the
agent a deterministic data surface.
"""

from __future__ import annotations

import os
import shutil
import textwrap
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from tempfile import mkdtemp


_FAKE_ORCA_OUT = textwrap.dedent("""\
|  1> ! {functional} def2-TZVP Opt Freq
Number of atoms                                 ...    12

  *** SCF CONVERGED ***
FINAL SINGLE POINT ENERGY    {spe}

------------------
ORBITAL ENERGIES
------------------

  NO   OCC          E(Eh)            E(eV)
   0   2.0000     -19.246153      -523.7146
   1   2.0000      -1.005211       -27.3535
   2   2.0000      -0.502350       -13.6697
   3   2.0000      -0.412341       -11.2200
   4   2.0000      -0.150000        -4.0816
   5   0.0000       0.052345         1.4244
   6   0.0000       0.412345        11.2200
   7   0.0000       0.812345        22.1051

------------------
VIBRATIONAL FREQUENCIES
-----------------------
   0:    -45.20 cm**-1 ***imaginary mode***
   1:    250.00 cm**-1
   2:    800.00 cm**-1
   3:   1450.00 cm**-1

THERMOCHEMISTRY

Temperature        ...    298.15 K
Pressure           ...      1.00 atm
Zero point energy                ...      0.123456 Eh
Total Enthalpy                   ...   -123.444444 Eh
Final entropy term               ...     -0.022222 Eh
G-E(el)                          ...      0.090000 Eh
Final Gibbs free energy          ...   {gibbs} Eh

Number of imaginary frequencies   ...     {n_imag}

ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-------------------------------------------------------------------
States    Energy (eV)   Energy (cm-1)   Wavelength (nm)   fosc  D**2  X     Y     Z
-------------------------------------------------------------------
  0-1A  ->  1-3A    2.600831     20977.1     476.7     0.000000    0.0    0.0  0.0  0.0
  0-1A  ->  2-3A    3.105210     25048.8     399.2     0.054321    0.5    0.1  0.2  0.4
  0-1A  ->  3-3A    4.005210     32299.1     309.6     0.123456    1.2    0.0  0.5  1.0

DIPOLE MOMENT
-------------
Total Dipole Moment    :     1.20000      -0.30000       0.40000

Magnitude (a.u.)       :      1.29614
Magnitude (Debye)      :      3.29447

                          *** OPTIMIZATION RUN DONE ***

TOTAL RUN TIME: 0 days 0 hours 12 minutes 34 seconds 567 msec
""")


_FAKE_CONTROL = textwrap.dedent("""\
# DELFIN CONTROL.txt (sandbox stub)
functional = BP86
main_basisset = def2-SVP
metal_basisset = def2-TZVP
disp_corr = D3BJ
solvent =
solvation_model =
freq_type = analytical
geom_opt = TightOpt
PAL = 8
maxcore = 2000
charge = 0
multiplicity = 1
""")


@dataclass
class Sandbox:
    """A self-contained DELFIN project tree the agent can poke at safely."""
    root: Path
    calc: Path
    archive: Path
    workspace: Path
    literature: Path | None = None  # papers/ for ad-hoc PDF reads

    def add_calc_folder(
        self,
        name: str,
        *,
        functional: str = "BP86",
        gibbs: float = -100.0,
        spe: float = -100.5,
        n_imag: int = 1,
    ) -> Path:
        folder = self.calc / name
        folder.mkdir(parents=True, exist_ok=True)
        out = _FAKE_ORCA_OUT.format(
            functional=functional,
            spe=spe,
            gibbs=gibbs,
            n_imag=n_imag,
        )
        (folder / "calc.out").write_text(out, encoding="utf-8")
        (folder / "initial.xyz").write_text(
            "3\nstub xyz\nC 0 0 0\nO 1 0 0\nH -1 0 0\n",
            encoding="utf-8",
        )
        (folder / "CONTROL.txt").write_text(
            _FAKE_CONTROL.replace("BP86", functional),
            encoding="utf-8",
        )
        return folder


def build_default_sandbox() -> Sandbox:
    """Create a tmp DELFIN tree with a tiny calc/ + archive/ + workspace.

    The tree is *not* auto-cleaned — wrap usage in :func:`sandboxed` (a
    context manager) for tests, or call :meth:`Sandbox.cleanup` manually.
    """
    root = Path(mkdtemp(prefix="delfin_dryrun_"))
    calc = root / "calc"
    archive = root / "archive"
    workspace = root / "agent_workspace"
    literature = root / "literature"
    calc.mkdir()
    archive.mkdir()
    workspace.mkdir()
    literature.mkdir()
    (root / "CONTROL.txt").write_text(_FAKE_CONTROL, encoding="utf-8")
    # A minimal valid PDF (1 page, says "DELFIN dry-run sandbox") so the
    # ad-hoc PDF read test isn't trying to read a non-existent file.
    pdf_bytes = (
        b"%PDF-1.4\n"
        b"1 0 obj<</Type/Catalog/Pages 2 0 R>>endobj\n"
        b"2 0 obj<</Type/Pages/Kids[3 0 R]/Count 1>>endobj\n"
        b"3 0 obj<</Type/Page/Parent 2 0 R/MediaBox[0 0 612 792]/Contents 4 0 R"
        b"/Resources<</Font<</F1 5 0 R>>>>>>endobj\n"
        b"4 0 obj<</Length 53>>stream\n"
        b"BT /F1 12 Tf 50 700 Td (DELFIN dry-run sandbox PDF) Tj ET\n"
        b"endstream endobj\n"
        b"5 0 obj<</Type/Font/Subtype/Type1/BaseFont/Helvetica>>endobj\n"
        b"xref\n0 6\n0000000000 65535 f \n"
        b"0000000010 00000 n \n0000000053 00000 n \n0000000098 00000 n \n"
        b"0000000179 00000 n \n0000000272 00000 n \n"
        b"trailer<</Size 6/Root 1 0 R>>\nstartxref\n333\n%%EOF\n"
    )
    (literature / "example.pdf").write_bytes(pdf_bytes)
    sb = Sandbox(
        root=root, calc=calc, archive=archive, workspace=workspace,
        literature=literature,
    )
    # Three fake calcs spanning functionals — perfect for compare tests.
    sb.add_calc_folder("test_a", functional="BP86",  gibbs=-100.0, spe=-100.5, n_imag=1)
    sb.add_calc_folder("bp86",   functional="BP86",  gibbs=-150.0, spe=-150.7, n_imag=0)
    sb.add_calc_folder("pbe0",   functional="PBE0",  gibbs=-150.5, spe=-151.2, n_imag=0)
    sb.add_calc_folder("b3lyp",  functional="B3LYP", gibbs=-149.8, spe=-150.3, n_imag=2)
    return sb


def cleanup(sb: Sandbox) -> None:
    if sb.root.exists():
        shutil.rmtree(sb.root, ignore_errors=True)


@contextmanager
def sandboxed():
    """Context manager: build a sandbox, yield it, tear it down on exit."""
    sb = build_default_sandbox()
    try:
        yield sb
    finally:
        cleanup(sb)
