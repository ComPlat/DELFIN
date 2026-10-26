"""Adversarial tests for the XT1 accept engine and uhf trap (phase 3).

Reviewer nacht-s17, package XT1. Attacks the recomputation acceptance
scripts in delfin/agent/pack/benchmark/accept/:
  - w20_xt1_minimum.py (SYSTEMS table, band semantics, folder lookup)
  - w20_xt1_uhf_trap.py (the doublet-band + expected_signals trap)
The tests import the scripts as modules and exercise main() through
their real argv contract with throwaway workspaces.
"""

import importlib.util
from pathlib import Path

import pytest

_SPEC_MIN = importlib.util.spec_from_file_location(
    "w20_xt1_minimum",
    Path("delfin/agent/pack/benchmark/accept/w20_xt1_minimum.py"))
_mod_min = importlib.util.module_from_spec(_SPEC_MIN)
_SPEC_MIN.loader.exec_module(_mod_min)

_SPEC_TRAP = importlib.util.spec_from_file_location(
    "w20_xt1_uhf_trap",
    Path("delfin/agent/pack/benchmark/accept/w20_xt1_uhf_trap.py"))
_mod_trap = importlib.util.module_from_spec(_SPEC_TRAP)
_SPEC_TRAP.loader.exec_module(_mod_trap)


class TestSystemTable:
    def test_unknown_key_refused_before_any_workspace_io(self, tmp_path):
        # An unknown mol_key must exit with the usage error (rc 2) BEFORE
        # any filesystem read of the workspace is attempted.
        assert _mod_min.main(["prog", str(tmp_path), "no_such_mol"]) == 2

    def test_every_system_key_has_a_task(self):
        # The SYSTEMS table grades the five non-trap systems exactly; the
        # trap task (methoxy_radical) is graded by w20_xt1_uhf_trap.py.
        tasks = ("water_cation", "oh_anion", "butanone", "ammonium",
                 "methylamine_radical", "methoxy_radical")
        assert set(_mod_min.SYSTEMS) == set(tasks)
        # The trap recompute (uhf "2") is handled by the dedicated trap
        # script, not by the SYSTEMS table; the table's methoxy entry is
        # the doublet physics (uhf 1), consistent with the trap's bands.
        assert _mod_min.SYSTEMS["methoxy_radical"][2] == 1

    def test_uhf_matches_open_shell_semantics(self):
        # Radicals (13/15 valence electrons) carry uhf 1; closed shells
        # carry uhf 0; the water cation (odd electrons) uhf 1. A wrong
        # pairing would grade the wrong spin state as "the physics".
        assert _mod_min.SYSTEMS["water_cation"][2] == 1
        assert _mod_min.SYSTEMS["oh_anion"][2] == 0
        assert _mod_min.SYSTEMS["butanone"][2] == 0
        assert _mod_min.SYSTEMS["ammonium"][2] == 0
        for k in ("methylamine_radical", "methoxy_radical"):
            assert _mod_min.SYSTEMS[k][2] == 1

    def test_imag_band_honesty(self):
        # The radicals grade {0,1} imag freqs (documented reproducible
        # saddle physics); the closed shells grade exactly 0. A band
        # widened to accept ANY count (e.g. (0, 99)) would hide physics.
        for k, v in _mod_min.SYSTEMS.items():
            lo, hi = v[5]
            if k in ("methylamine_radical", "methoxy_radical"):
                assert (lo, hi) == (0, 1)
            else:
                assert (lo, hi) == (0, 0)
            assert hi <= 1

    def test_energy_tolerances_are_tight(self):
        # Bands must stay tight enough that a wrong method (gfn1, gfnff,
        # vacuum) cannot pass: worst-case band is +-0.010 Eh = 6 kcal/mol;
        # a GFN1 energy of the same molecule differs by several tenths.
        for v in _mod_min.SYSTEMS.values():
            assert v[4] <= 0.010
            assert v[7] <= 0.5


class TestFolderLookup:
    def test_missing_folder_rc1(self, tmp_path):
        assert _mod_min.main(["prog", str(tmp_path), "butanone"]) == 1

    def test_missing_folder_rc1_trap(self, tmp_path):
        assert _mod_trap.main(["prog", str(tmp_path)]) == 1

    def test_start_xyz_fallback_used(self, tmp_path, monkeypatch):
        # No xtbopt.xyz, one xyz present: the fallback picks the xyz --
        # the engine must not silently read a nonexistent start.xyz.
        # xtb is a stand-in that records the structure it was handed and
        # leaves no result: the test is about which file is read, and must
        # not depend on an xtb being installed (a scheduled CI run without
        # one ended in SystemExit 3 from the resolver).
        d = tmp_path / "xtb" / "butanone"
        d.mkdir(parents=True)
        start = "13\nbutanone\nC 0 0 0\n"
        (d / "opt1.xyz").write_text(start)
        seen = tmp_path / "seen.xyz"
        fake = tmp_path / "xtb-standin"
        fake.write_text(f"#!/bin/sh\ncat \"$1\" > '{seen}'\nexit 1\n")
        fake.chmod(0o755)
        monkeypatch.setattr(_mod_min, "_xtb", lambda: str(fake))
        assert _mod_min.main(["prog", str(tmp_path), "butanone"]) == 1
        assert seen.read_text() == start

    def test_workspace_is_resolved_not_created(self, tmp_path):
        # A nonexistent workspace path must be refused, not mkdir'd.
        ghost = tmp_path / "ghost"
        rc = _mod_min.main(["prog", str(ghost), "butanone"])
        assert rc == 1
        assert not ghost.exists()
