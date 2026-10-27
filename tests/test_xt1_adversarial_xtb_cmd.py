"""Adversarial tests for the native xTB adapter command builder (XT1 review).

Reviewer nacht-s17, package XT1. These tests attack ``build_xtb_cmd`` and the
charge/spin wiring in ``delfin/tools/adapters/xtb_native.py`` — the function
every native xtb adapter calls (xtb_native.py:159) before running the binary.
They run against the builder branch agent/s10-xt1b unchanged; every failure
here is a finding for nacht-s10, not a patch by us.
"""

import pytest

from delfin.tools.adapters.xtb_native import build_xtb_cmd


class TestMethodSelection:
    """Unknown or malformed method names must not silently fall back."""

    @pytest.mark.parametrize("method,expected_flags", [
        ("gfn0", ["--gfn", "0"]),
        ("gfn1", ["--gfn", "1"]),
        ("gfn2", ["--gfn", "2"]),
        ("gfnff", ["--gfnff"]),
    ])
    def test_known_methods_map_to_flags(self, method, expected_flags):
        cmd = build_xtb_cmd("mol.xyz", method=method, charge=0)
        assert expected_flags == cmd[2:2 + len(expected_flags)]

    def test_unknown_method_is_rejected_not_silently_gfn2(self):
        # A typo like "gfn3" or "GFN2" (case) must not be quietly mapped onto
        # GFN2: the run would produce plausible numbers at a different level.
        with pytest.raises((ValueError, KeyError)):
            build_xtb_cmd("mol.xyz", method="gfn3", charge=0)

    def test_method_none_is_rejected(self):
        with pytest.raises((ValueError, KeyError, TypeError)):
            build_xtb_cmd("mol.xyz", method=None, charge=0)


class TestChargeSpinWiring:
    """--uhf comes from multiplicity (xtb_native.py:161); invalid values must
    fail loudly, not produce an invalid command line."""

    def test_doublet_gives_uhf_1(self):
        cmd = build_xtb_cmd("mol.xyz", charge=0, uhf=1)
        assert ["--uhf", "1"] == cmd[cmd.index("--uhf"):cmd.index("--uhf") + 2]

    def test_uhf_negative_is_rejected(self):
        # mult=0 or mult=-3 would produce --uhf -1 / --uhf -4: xtb cannot
        # interpret a negative unpaired-electron count.
        with pytest.raises((ValueError, AssertionError)):
            build_xtb_cmd("mol.xyz", charge=0, uhf=-1)

    def test_uhf_absurdly_large_is_rejected(self):
        # A multiplicity like 1000 (uhf=999) exceeds any real electron count.
        with pytest.raises((ValueError, AssertionError)):
            build_xtb_cmd("mol.xyz", charge=0, uhf=999)


class TestSolventModelWiring:
    """Solvent/model values end up as --<model> <solvent>; an unknown model
    name would build a flag xtb does not know."""

    def test_alpb_water(self):
        cmd = build_xtb_cmd("mol.xyz", charge=0, solvent="water", solvent_model="alpb")
        assert ["--alpb", "water"] == cmd[cmd.index("--alpb"):cmd.index("--alpb") + 2]

    def test_gbsa_acetonitrile(self):
        cmd = build_xtb_cmd("mol.xyz", charge=0, solvent="acetonitrile", solvent_model="gbsa")
        assert ["--gbsa", "acetonitrile"] == cmd[cmd.index("--gbsa"):cmd.index("--gbsa") + 2]

    def test_unknown_solvent_model_is_rejected(self):
        # e.g. "cpcm" is an ORCA model, not an xtb one — passing it must fail
        # at command-build time, not after a wasted binary run.
        with pytest.raises((ValueError, KeyError, AssertionError)):
            build_xtb_cmd("mol.xyz", charge=0, solvent="water", solvent_model="cpcm")

    def test_solvent_model_without_solvent_is_pointless(self):
        # --alpb with no solvent name produces a bare flag xtb rejects;
        # the builder should either reject or drop the model flag.
        cmd = build_xtb_cmd("mol.xyz", charge=0, solvent_model="alpb")
        assert "--alpb" not in cmd or any(a.startswith("-") is False and a not in ("--alpb",) for a in cmd[cmd.index("--alpb") + 1:]) if "--alpb" in cmd else True


class TestExtraArgsSurface:
    """extra_args is an escape hatch appended verbatim to the command line
    (xtb_native.py:99). It must stay a token list — a string with spaces must
    not silently split into flags, and it must not be able to override the
    validated --chrg/--uhf."""

    def test_extra_args_string_is_rejected_not_silently_split(self):
        with pytest.raises((ValueError, TypeError)):
            build_xtb_cmd("mol.xyz", charge=0, extra_args="--gfn 1")

    def test_extra_args_cannot_override_chrg(self):
        # Rewritten along its intent for the builder's fix (37235359): the
        # old behaviour appended the duplicate "--chrg -2" verbatim (xtb
        # takes the last value), so the test asserted the override never
        # happened. The fix now refuses the protected flag outright; the
        # intent — extra_args cannot override the validated charge — is
        # measured by the refusal, which is strictly stronger.
        with pytest.raises((ValueError, TypeError)):
            build_xtb_cmd("mol.xyz", charge=1, extra_args=("--chrg", "-2"))
        # And a near-miss alias must not slip past either.
        with pytest.raises((ValueError, TypeError)):
            build_xtb_cmd("mol.xyz", charge=1, extra_args=("--chrg=-2",))
