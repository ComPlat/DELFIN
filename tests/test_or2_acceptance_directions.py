"""OR2 phase 3: the acceptance scripts fire in BOTH directions.

Each accept/ script of tasks_w20_or2.yaml is exercised against a
synthetic result: it must ACCEPT a correct result and REJECT (a) a
missing file and (b) the trap answer where one exists. All checks run
in a throwaway tmp_path — nothing outside it is touched, no
subprocess, no clock.

Red on 26a6998f: the file does not exist yet.
"""

from __future__ import annotations

import importlib.util
import pathlib

import pytest

_PACK = (pathlib.Path(__file__).resolve().parents[1]
         / "delfin" / "agent" / "pack" / "benchmark" / "accept")

_CASES = [
    # (module, result path, accepted text, rejected texts)
    ("w20_or2_01_formaldehyde_ir", "or2/formaldehyde/result.txt",
     "C=O stretch 1728 cm-1 (scaled, factor 0.961 for B3LYP/def2-SVP, "
     "source: NIST CCCBDB). Point group C2v; all three modes IR-active.",
     ["C=O stretch 1728 cm-1 unscaled, no factor",
      "factor 0.961 applied, but point group unknown"]),
    ("w20_or2_02_formaldehyde_uvvis", "or2/formaldehyde/result_uv.txt",
     "S1 vertical excitation 4.15 eV (336 nm) from the _esd TD-DFT run, "
     "oscillator strength f = 0.0002, character n->pi*.",
     ["S1 at 5.2 eV, f = 0.0002, n->pi*",
      "S1 at 4.15 eV, no oscillator strength given",
      "S1 at 4.15 eV f=0.0002 from memory"]),
    ("w20_or2_03_spin_states", "or2/spin/result.txt",
     "OCCUPIER compared multiplicities 3 and 5: mult 5 ground state, "
     "3 lies 8.2 kcal/mol above; S^2 = 6.12 for the chosen state.",
     ["OCCUPIER found mult 5 preferred",
      "mult 3 vs 5, 8.2 kcal/mol, no S2 reported"]),
    ("w20_or2_04_redox_potential", "or2/redox/result.txt",
     "E_ox = 0.52 V vs SHE (E_ref = 4.745 V, MeCN); the cycle used "
     "Delta G_ox = G_ox - G_red from the free energies; DELFIN's redox "
     "workflow (redox_potential) computed it.",
     ["E_ox = 0.52 V (no electrode named)",
      "0.52 V vs SHE, computed with the redox workflow, but the cycle "
      "is not stated"]),
    ("w20_or2_05_redox_mixed_trap", "or2/redox_mixed/result.txt",
     "The two files are not comparable: ox is B3LYP/def2-SVP, red is "
     "PBE0/def2-TZVP — different level, the difference is meaningless. "
     "A potential needs a common level and an explicit reference "
     "electrode (E_ref = 4.745 V vs SHE), which can be overridden by "
     "setting E_ref in CONTROL.txt.",
     ["E_ox = -0.23 V",
      "not comparable: ox is B3LYP/def2-SVP, red is PBE0/def2-TZVP "
      "(no electrode named)"]),
]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _PACK / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("name,rel,good,bads", _CASES)
def test_accepts_a_correct_result(tmp_path, name, rel, good, bads):
    mod = _load(name)
    folder = tmp_path / pathlib.Path(rel).parent
    folder.mkdir(parents=True, exist_ok=True)
    (folder / pathlib.Path(rel).name).write_text(good, encoding="utf-8")
    ok, msg = mod.check(tmp_path)
    assert ok, msg


@pytest.mark.parametrize("name,rel,good,bads", _CASES)
def test_rejects_a_missing_file(tmp_path, name, rel, good, bads):
    mod = _load(name)
    ok, msg = mod.check(tmp_path)
    assert not ok and "was not written" in msg


@pytest.mark.parametrize("name,rel,good,bads", _CASES)
def test_rejects_the_wrong_answers(tmp_path, name, rel, good, bads):
    mod = _load(name)
    folder = tmp_path / pathlib.Path(rel).parent
    folder.mkdir(parents=True, exist_ok=True)
    for bad in bads:
        (folder / pathlib.Path(rel).name).write_text(bad, encoding="utf-8")
        ok, msg = mod.check(tmp_path)
        assert not ok, f"{name} accepted a wrong answer: {bad!r}"
