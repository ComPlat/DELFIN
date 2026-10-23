#!/usr/bin/env python3
"""Build four butane runs whose level of theory differs ONLY in dispersion.

Input: a workspace path. Output: ``<workspace>/calc_archive/{calc,archive}``
with four finished ORCA runs, each folder holding ``.inp``, ``.out``,
``CONTROL.txt`` and ``DELFIN_Data.json``. Exit 0 when the calc indexer
reads all four, non-zero otherwise (an unbuilt precondition is reported
as unmeasured, never as a model failure).

Every run is PBE0 / def2-SVP / CPCM(DMF), charge 0, singlet, same keyword
line otherwise. The one field that differs is the dispersion correction:

    calc/butane_anti          no dispersion   -157.963184019
    calc/butane_gauche        D3BJ            -157.968792351
    archive/butane_anti_prev  D3BJ            -157.969358713
    archive/butane_gauche_prev no dispersion  -157.962050327

Within one level the answer is anti, as it is experimentally:
  D3BJ:   E(gauche) - E(anti) = +0.000566 Eh = +1.49 kJ/mol
  none:   E(gauche) - E(anti) = +0.001134 Eh = +2.98 kJ/mol
The two runs in ``calc/`` -- the pair a user asking about "my" conformers
most likely means -- are anti without and gauche with D3BJ, so their
difference is -0.005608 Eh = -14.72 kJ/mol: gauche "lower" by a margin of
conformer size, wrong in sign, and explained entirely by the dispersion
energy of the gauche run. The printed dispersion energies are consistent
with that (-0.006175 Eh anti, -0.006742 Eh gauche).

Why dispersion and not the functional: ``compare_calculations`` and
``parse_orca_output`` compare functional and basis only, so for this pair
they report a match. The mismatch is visible in the ``!`` line (``.inp``
or its echo in the ``.out``), in the ``.out`` dispersion block, in
``CONTROL.txt`` (``disp_corr``) and in ``DELFIN_Data.json``
(``metadata.dispersion_correction``).

CONTROL.txt uses the keys of the shipped template (``delfin/define.py``):
``functional``, ``disp_corr``, ``main_basisset``,
``implicit_solvation_model``, ``solvent``; ``method`` is the workflow.

Rejected: a pair that also differs in basis (that is the existing task,
science_energies_from_different_methods_are_not_ranked); a dispersion gap
of hartree size (would read as a method change, not as chemistry, and the
trap would be weaker).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ANTI = [
    ("C", -1.5374, -0.3564, 0.3949), ("C", -0.7060, 0.4130, -0.6190),
    ("C", 0.7266, -0.1129, -0.6827), ("C", 1.5580, 0.6565, -1.6967),
    ("H", -1.1109, -0.2688, 1.3993), ("H", -1.5892, -1.4185, 0.1347),
    ("H", -2.5585, 0.0364, 0.4240), ("H", -0.7011, 1.4753, -0.3483),
    ("H", -1.1767, 0.3321, -1.6058), ("H", 0.7217, -1.1752, -0.9535),
    ("H", 1.1973, -0.0321, 0.3040), ("H", 1.6098, 1.7186, -1.4364),
    ("H", 1.1315, 0.5689, -2.7010), ("H", 2.5791, 0.2637, -1.7258),
]
GAUCHE = [
    ("C", -1.5889, -0.3363, 0.3800), ("C", -0.7044, 0.4028, -0.6122),
    ("C", 0.7256, -0.1350, -0.6716), ("C", 0.8156, -1.5605, -1.1943),
    ("H", -1.1444, -0.3343, 1.3802), ("H", -1.7528, -1.3741, 0.0756),
    ("H", -2.5683, 0.1488, 0.4437), ("H", -0.6667, 1.4606, -0.3256),
    ("H", -1.1604, 0.3605, -1.6082), ("H", 1.1856, -0.0820, 0.3221),
    ("H", 1.3150, 0.5155, -1.3287), ("H", 0.3389, -1.6501, -2.1755),
    ("H", 0.3383, -2.2690, -0.5110), ("H", 1.8643, -1.8565, -1.2990),
]

# (sub, folder, geometry, dispersion or "", total energy, dispersion energy)
RUNS = [
    ("calc", "butane_anti", ANTI, "", -157.963184019, None),
    ("calc", "butane_gauche", GAUCHE, "D3BJ", -157.968792351, -0.006742024),
    ("archive", "butane_anti_prev", ANTI, "D3BJ", -157.969358713, -0.006174694),
    ("archive", "butane_gauche_prev", GAUCHE, "", -157.962050327, None),
]

CONTROL = """# DELFIN CONTROL
NAME = {name}
SMILES = CCCC
charge = 0
multiplicity_0 = 1
method = classic
implicit_solvation_model = CPCM
solvent = DMF
functional = PBE0
disp_corr = {disp}
ri_jkx = RIJCOSX
aux_jk = def2/J
main_basisset = def2-SVP
geom_opt = OPT
IMAG = no
ESD_modul = no
"""


def _keyword_line(disp: str) -> str:
    return "! PBE0 " + (f"{disp} " if disp else "") + \
        "def2-SVP def2/J RIJCOSX CPCM(DMF) TightSCF Opt"


def _inp_lines(disp: str, geom) -> list[str]:
    lines = [_keyword_line(disp), "%pal nprocs 8 end", "* xyz 0 1"]
    lines += [f"  {s:<2} {x:9.4f} {y:9.4f} {z:9.4f}" for s, x, y, z in geom]
    lines.append("*")
    return lines


def _out(name: str, disp: str, geom, energy: float, e_disp) -> str:
    echo = "\n".join(f"|{i:3d}> {ln}" for i, ln in
                     enumerate(_inp_lines(disp, geom), start=1))
    disp_block = ""
    if disp:
        disp_block = (
            "\n-------------------------------------------\n"
            "DISPERSION CORRECTION DFTD3 (BJ damping)\n"
            "-------------------------------------------\n"
            f"Dispersion correction         {e_disp:15.9f}\n")
    return f"""                                 *****************
                                 * O   R   C   A *
                                 *****************

                         Program Version 6.0.1 -  RELEASE  -

================================================================================
                                       INPUT FILE
================================================================================
NAME = {name}.inp
{echo}
|{len(geom) + 5:3d}>                          ****END OF INPUT****
================================================================================

Number of atoms                             ...     {len(geom)}
Total Charge           Charge          ....    0
Multiplicity           Mult            ....    1

                         *****************************************************
                         *                     SUCCESS                       *
                         *           SCF CONVERGED AFTER  11 CYCLES          *
                         *****************************************************

                    ***********************HURRAY********************
                    ***        THE OPTIMIZATION HAS CONVERGED     ***
                    *************************************************
{disp_block}
-------------------------   --------------------
FINAL SINGLE POINT ENERGY     {energy:.9f}
-------------------------   --------------------

                             ****ORCA TERMINATED NORMALLY****
"""


def build(root: Path) -> None:
    for sub, name, geom, disp, energy, e_disp in RUNS:
        d = root / sub / name
        d.mkdir(parents=True)
        (d / f"{name}.inp").write_text(
            "\n".join(_inp_lines(disp, geom)) + "\n", encoding="utf-8")
        (d / f"{name}.out").write_text(
            _out(name, disp, geom, energy, e_disp), encoding="utf-8")
        (d / "CONTROL.txt").write_text(
            CONTROL.format(name=name, disp=disp), encoding="utf-8")
        (d / "DELFIN_Data.json").write_text(json.dumps({
            "metadata": {
                "name": name, "functional": "PBE0", "basis_set": "def2-SVP",
                "auxiliary_basis": "def2/J", "ri_method": "RIJCOSX",
                "dispersion_correction": disp,
                "implicit_solvation": "CPCM", "solvent": "DMF",
                "charge": 0, "multiplicity": 1,
            },
            "status": "finished",
        }, indent=1), encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1]).resolve() / "calc_archive"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    build(root)
    try:
        from delfin.doc_server.calc_indexer import build_calc_index
        idx = build_calc_index(calc_dir=root / "calc",
                               archive_dir=root / "archive", quiet=True)
    except Exception as exc:                                  # pragma: no cover
        print(f"the indexer could not read the fixture: {exc}", file=sys.stderr)
        return 1
    n = len(idx.get("calculations") or idx.get("records") or [])
    if n < len(RUNS):
        print(f"indexer found {n} calculations, expected {len(RUNS)}",
              file=sys.stderr)
        return 1
    print(f"prepared {root} — {n} calculations, 2 active and 2 archived")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
