"""XT1 F7/F9: the solvent flag emits the VALIDATED model and a solvent name.

F7 (nacht-s17): line 150 built f"--{solvent_model}" from the RAW
parameter, so "ALPB" wrote the invalid flag "--ALPB" and "" wrote a
bare "--". The lowercased, validated ``model`` (:109) is what belongs
on the command line.

F9 (nacht-s17): a flag-shaped solvent NAME collides with the model
flag -- solvent="--alpb" produced [..., "--alpb", "--alpb"] where xtb
consumes its own duplicate as the value. Flag-shaped names are
refused.
"""

import pytest

from delfin.tools.adapters.xtb_native import build_xtb_cmd


def test_uppercase_model_writes_lowercase_flag():
    cmd = build_xtb_cmd("mol.xyz", solvent="water", solvent_model="ALPB")
    assert "--alpb" in cmd
    assert cmd[cmd.index("--alpb") + 1] == "water"


def test_empty_model_falls_back_to_alpb():
    cmd = build_xtb_cmd("mol.xyz", solvent="water", solvent_model="")
    assert "--alpb" in cmd


def test_gbsa_selects_gbsa_flag():
    cmd = build_xtb_cmd("mol.xyz", solvent="water", solvent_model="gbsa")
    assert "--gbsa" in cmd
    assert "--alpb" not in cmd


def test_solvent_name_passes_through_verbatim():
    cmd = build_xtb_cmd("mol.xyz", solvent="CH2Cl2")
    assert cmd[cmd.index("--alpb") + 1] == "CH2Cl2"


def test_flag_shaped_solvent_name_is_refused():
    with pytest.raises(ValueError, match="solvent NAME"):
        build_xtb_cmd("mol.xyz", solvent="--alpb")
    with pytest.raises(ValueError, match="solvent NAME"):
        build_xtb_cmd("mol.xyz", solvent="-water")
