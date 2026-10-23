"""Adversarial tests, solvent surface of build_xtb_cmd (XT1 review).

Reviewer nacht-s17, package XT1. Attacks the implicit-solvation path in
``delfin/tools/adapters/xtb_native.py``: model restriction (ALPB/GBSA),
protected flags around solvation, and pass-through of the solvent name.
"""

import pytest

from delfin.tools.adapters.xtb_native import build_xtb_cmd


class TestSolventModel:
    def test_alpb_default_wiring(self):
        cmd = build_xtb_cmd("mol.xyz", solvent="water")
        assert cmd[cmd.index("--alpb") + 1] == "water"

    def test_model_case_insensitive(self):
        # The fix lowercases the model (:109); uppercase must behave like
        # the lowercase form, not be refused or passed through raw.
        cmd = build_xtb_cmd("mol.xyz", solvent="water", solvent_model="ALPB")
        assert "--alpb" in cmd

    def test_unknown_model_refused(self):
        # F3 fix: CPCM and friends are ORCA-only, refused with a pointer.
        with pytest.raises(ValueError):
            build_xtb_cmd("mol.xyz", solvent="water", solvent_model="cpcm")

    def test_empty_model_falls_back_to_alpb(self):
        cmd = build_xtb_cmd("mol.xyz", solvent="water", solvent_model="")
        assert "--alpb" in cmd

    def test_gbsa_model_selects_gbsa_flag(self):
        cmd = build_xtb_cmd("mol.xyz", solvent="water", solvent_model="gbsa")
        assert "--gbsa" in cmd
        assert "--alpb" not in cmd


class TestSolventExtraArgs:
    def test_protected_solvent_flag_blocked(self):
        # Protected set (:117) covers --alpb/--gbsa in token form.
        with pytest.raises((ValueError, TypeError)):
            build_xtb_cmd("mol.xyz", extra_args=("--alpb", "water"))

    def test_assignment_form_of_solvent_flag(self):
        # F6 shape: the `=`-attached assignment --alpb=water must be
        # refused exactly like the token form, or the model choice in the
        # command line can be injected despite the validation.
        with pytest.raises((ValueError, TypeError)):
            build_xtb_cmd("mol.xyz", extra_args=("--alpb=water",))

    def test_solvent_name_value_pass_through(self):
        # The solvent NAME (xtb list, e.g. acetone, ch2cl2) is validated
        # by xtb itself at runtime; the adapter passing it through is the
        # documented contract (:8). Pinned as characterization.
        cmd = build_xtb_cmd("mol.xyz", solvent="CH2Cl2")
        assert cmd[cmd.index("--alpb") + 1] == "CH2Cl2"

    def test_solvent_without_value_is_injectable(self):
        # A bare "--alpb" as the solvent NAME would collide with the
        # model flag: cmd becomes [..., "--alpb", "--alpb"] and the flag
        # consumes the next argument as its value. Must be refused or at
        # minimum not silently produce a flag-as-value command.
        cmd = build_xtb_cmd("mol.xyz", solvent="--alpb")
        values = cmd[cmd.index("--alpb") + 1] if "--alpb" in cmd else None
        assert values != "--alpb" if values is not None else True
