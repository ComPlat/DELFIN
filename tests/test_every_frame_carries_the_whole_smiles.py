"""Every MANTA frame has the sum formula of the SMILES it was built from.

Two defects looked like one: ``[Pt](Cl)(Cl)(N)N`` came out with N atoms that had
lost their hydrogens on some construction paths and kept them on others (frames of
5 and of 9 atoms in one manifold), and ``[Pt](Cl)(Cl)(N)N.Cl`` therefore mixed 7-
and 11-atom frames -- which read as a lost counter-ion but was the same hydrogen
loss; the ``Cl`` was in every frame.  The root: an unbracketed atom's hydrogens are
implicit, i.e. recomputed from its valence, and preparation paths that mark every
metal-bound atom ``NoImplicit`` or cleave the M-L bond change what that recomputation
gives.  The fix pins the RDKit count into the SMILES once, at the public entry
(``_pin_metal_neighbour_hydrogens``: ``N`` -> ``[NH2]``).

The expected formula is RDKit's own reading of the SMILES: heavy atoms as written,
hydrogens as ``GetTotalNumHs`` assigns them (explicit counts for bracket atoms,
valence-derived counts for unbracketed ones).  Each frame is built through the
real ``delfin-manta`` in a fresh interpreter.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]

_CASES = (
    "[Pt](Cl)(Cl)(N)N",
    "[Pt](Cl)(Cl)([NH3])[NH3]",
    "[Pt](Cl)(Cl)(N)N.Cl",
    "[Pt](Cl)(Cl)([NH3])[NH3].[Cl-]",
)


def _expected_formula(smiles: str) -> Counter:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    formula = Counter(a.GetSymbol() for a in mol.GetAtoms())
    formula["H"] += sum(a.GetTotalNumHs() for a in mol.GetAtoms())
    return formula


def _frames_from_cli(smiles: str, out: Path) -> list:
    env = {k: v for k, v in os.environ.items() if not k.startswith("DELFIN_")}
    env["PYTHONPATH"] = str(_ROOT) + (os.pathsep + env["PYTHONPATH"]
                                      if env.get("PYTHONPATH") else "")
    proc = subprocess.run(
        [sys.executable, "-m", "delfin.cli_manta", smiles, "-o", str(out), "-q"],
        env=env, capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, proc.stderr[-3000:]
    manifest = json.loads((out / "manifest.json").read_text())
    frames = []
    for item in manifest["isomers"]:
        lines = (out / item["file"]).read_text().splitlines()[2:]
        frames.append((item["label"],
                       Counter(ln.split()[0] for ln in lines if ln.strip())))
    return frames


@pytest.mark.parametrize("smiles", _CASES)
def test_every_frame_has_the_formula_of_the_smiles(smiles, tmp_path):
    want = _expected_formula(smiles)
    frames = _frames_from_cli(smiles, tmp_path / "out")
    assert frames, "nothing was built"
    wrong = [(label, dict(got)) for label, got in frames if got != want]
    assert not wrong, f"{smiles}: expected {dict(want)}, got {wrong[:5]}"


def test_the_expected_formulas_are_what_rdkit_reads():
    """Pinned so the expectation itself cannot drift: unbracketed N on Pt is NH2."""
    assert _expected_formula("[Pt](Cl)(Cl)(N)N") == Counter(
        {"Pt": 1, "Cl": 2, "N": 2, "H": 4})
    assert _expected_formula("[Pt](Cl)(Cl)([NH3])[NH3].[Cl-]") == Counter(
        {"Pt": 1, "Cl": 3, "N": 2, "H": 6})


def test_the_pin_rewrites_only_unbracketed_metal_neighbours_with_hydrogens():
    from delfin.smiles_converter import _pin_metal_neighbour_hydrogens as pin

    assert pin("[Pt](Cl)(Cl)(N)N") == "[Pt](Cl)(Cl)([NH2])[NH2]"
    assert pin("[Pt](Cl)(Cl)(N)N.Cl") == "[Pt](Cl)(Cl)([NH2])[NH2].Cl"
    assert pin("O[Cu]") == "[OH][Cu]"          # one implicit H, pinned
    assert pin("CO[Cu]") == "CO[Cu]"           # O valence full, nothing to pin
    # already explicit, or no H to lose, or no metal: the identical string
    for same in ("[Pt](Cl)(Cl)([NH3])[NH3]", "Cl[Co+3](Cl)([NH3])([NH3])([NH3])[NH3]",
                 "CCO", "c1ccccc1", "[Cl][Cd-3]([Cl])[N+]1=CC=CC=C1"):
        assert pin(same) == same
