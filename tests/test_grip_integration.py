"""Phase 4 — GRIP × fffree integration tests.

These tests cover the env-gated GRIP hook added in
``delfin/fffree/assemble_complex.py::assemble_from_config``.  They exercise
the production entry point ``_fffree_isomers`` (the same function the SMILES
converter calls) so the GRIP path is verified end-to-end exactly as the
voll-pool will hit it.

Verified properties:

1. Default-OFF byte-identical to HEAD (no env / ``DELFIN_FFFREE_GRIP``
   unset / explicitly ``0`` / any non-"1" value all produce identical XYZ).
2. ``DELFIN_FFFREE_GRIP=1`` runs end-to-end without raising and yields
   finite XYZ for a small panel of homoleptic / mixed-ligand complexes.
3. The M-D invariant is preserved within the SPEC §11 ±0.05 Å envelope
   when GRIP is active (read off the emitted XYZ).
4. ``grip_polish`` returning NaN coords cleanly falls back (never emits
   non-finite values).
5. ``grip_polish`` returning ``None`` cleanly falls back without crashing.
6. ``GripLibrary.get_default()`` is a process-wide singleton and emits a
   WARN log when the library file is missing.

The suite runs in a few seconds — all SMILES are small / monodentate
complexes that build via the standard fffree pipeline.
"""
from __future__ import annotations

import os

# Strict determinism — must be set before numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import logging  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402


# ---------------------------------------------------------------------------
# Test panel — every SMILES must build under DELFIN_FFFREE_BUILDER=1 + OFF.
# ---------------------------------------------------------------------------
PANEL = [
    "[Pd](Cl)(Cl)(c1ccncc1)c1ccncc1",         # Pd(py)2Cl2 (SP-4)
    "[Pt](Cl)(Cl)(c1ccncc1)c1ccncc1",         # Pt(py)2Cl2 (SP-4)
    "[Ni](Cl)(Cl)(c1ccncc1)(c1ccncc1)(c1ccncc1)c1ccncc1",  # Ni(py)4Cl2 (OC-6)
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _clear_grip_env():
    for k in (
        "DELFIN_FFFREE_GRIP",
        "DELFIN_GRIP_LIB_PATH",
    ):
        os.environ.pop(k, None)


def _force_builder_on():
    """fffree path requires DELFIN_FFFREE_BUILDER=1."""
    os.environ.setdefault("DELFIN_FFFREE_BUILDER", "1")


def _isomers(smiles: str, max_isomers: int = 2):
    """Production entry point — same as smiles_converter calls."""
    from delfin.fffree.converter_backend import _fffree_isomers
    return _fffree_isomers(smiles, max_isomers=max_isomers)


def _xyz_canonical(result) -> str:
    """Concatenate all isomer XYZs (label + body) into a single canonical
    string suitable for byte-equality comparison.  6-decimal precision is
    the same precision the SMILES converter emits to disk."""
    if result is None:
        return "BUILD_FAILED"
    parts = []
    for xyz, label in result:
        parts.append(f"### {label} ###")
        parts.append(xyz)
    return "\n".join(parts)


def _parse_xyz(xyz: str):
    """Parse an XYZ string into (syms, P)."""
    syms = []
    coords = []
    for ln in xyz.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        syms.append(parts[0])
        coords.append((x, y, z))
    return syms, np.asarray(coords, dtype=np.float64)


# ---------------------------------------------------------------------------
# 1. Default-OFF byte-identical
# ---------------------------------------------------------------------------
class TestDefaultOffByteIdentical:
    """The GRIP hook must be a strict no-op when its env flag is unset
    or set to anything other than "1" (SPEC §11)."""

    @pytest.mark.parametrize("smiles", PANEL)
    def test_unset_equals_explicit_zero(self, smiles):
        _force_builder_on()
        _clear_grip_env()
        baseline = _xyz_canonical(_isomers(smiles))

        os.environ["DELFIN_FFFREE_GRIP"] = "0"
        try:
            zero = _xyz_canonical(_isomers(smiles))
        finally:
            _clear_grip_env()

        assert baseline == zero, (
            f"unset vs explicit-0 differ for {smiles}; "
            f"len(unset)={len(baseline)} len(zero)={len(zero)}"
        )

    def test_arbitrary_non_one_values_are_off(self):
        _force_builder_on()
        _clear_grip_env()
        baseline = _xyz_canonical(_isomers(PANEL[0]))
        for val in ("", "0", "false", "off", "FALSE", "no", "2"):
            os.environ["DELFIN_FFFREE_GRIP"] = val
            try:
                got = _xyz_canonical(_isomers(PANEL[0]))
            finally:
                _clear_grip_env()
            assert got == baseline, (
                f"DELFIN_FFFREE_GRIP={val!r} must be OFF (byte-identical) "
                f"— SPEC §11"
            )


# ---------------------------------------------------------------------------
# 2. GRIP-ON runs end-to-end without error
# ---------------------------------------------------------------------------
class TestGripOnNoError:
    """With ``DELFIN_FFFREE_GRIP=1`` the build must complete and emit
    finite, well-shaped coordinates."""

    @pytest.mark.parametrize("smiles", PANEL)
    def test_grip_on_runs(self, smiles):
        _force_builder_on()
        os.environ["DELFIN_FFFREE_GRIP"] = "1"
        try:
            res = _isomers(smiles)
        finally:
            _clear_grip_env()
        assert res is not None, f"GRIP-on build returned None for {smiles}"
        assert len(res) >= 1
        for xyz, label in res:
            assert isinstance(xyz, str) and len(xyz) > 0
            syms, P = _parse_xyz(xyz)
            assert P.size > 0
            assert P.shape[1] == 3
            assert np.all(np.isfinite(P)), (
                f"non-finite XYZ in isomer {label} for {smiles}"
            )


# ---------------------------------------------------------------------------
# 3. M-D invariant preserved under GRIP
# ---------------------------------------------------------------------------
class TestGripMDInvariant:
    """Defence-in-depth.  Read M-D distances off the emitted XYZ for
    GRIP-OFF and GRIP-ON; any drift > 0.05 Å indicates a safety-net
    failure (the hook itself, grip_polish, AND the builder all have
    redundant M-D guards — they must agree)."""

    @pytest.mark.parametrize("smiles", PANEL)
    def test_md_within_tolerance(self, smiles):
        _force_builder_on()
        _clear_grip_env()
        off = _isomers(smiles)
        if off is None:
            pytest.skip("baseline build returned None")

        os.environ["DELFIN_FFFREE_GRIP"] = "1"
        try:
            on = _isomers(smiles)
        finally:
            _clear_grip_env()
        assert on is not None
        assert len(off) == len(on), "isomer count must be GRIP-invariant"

        for (xyz_off, lbl_off), (xyz_on, lbl_on) in zip(off, on):
            syms_off, P_off = _parse_xyz(xyz_off)
            syms_on, P_on = _parse_xyz(xyz_on)
            assert syms_off == syms_on, (
                f"atom symbol order changed for {smiles} "
                f"({lbl_off} vs {lbl_on})"
            )
            # Metal at index 0; donors == direct neighbours within a
            # generous cutoff (we don't have donor indices here, so we
            # take the nearest len_off neighbours instead — sufficient
            # for a M-D drift check, which is per-atom not per-donor).
            d_off = np.linalg.norm(P_off - P_off[0], axis=1)
            d_on = np.linalg.norm(P_on - P_on[0], axis=1)
            drift = np.abs(d_off - d_on)
            # First-shell atoms — anything within 3.0 Å of the metal.
            shell = d_off < 3.0
            shell[0] = False  # exclude metal itself
            if shell.any():
                max_drift = float(drift[shell].max())
                assert max_drift < 0.05 + 1e-6, (
                    f"first-shell M-X drift {max_drift:.4f} Å > 0.05 "
                    f"(SPEC §11) for {smiles} {lbl_off}"
                )


# ---------------------------------------------------------------------------
# 4. NaN-safety — contrived grip_polish returning NaN coords
# ---------------------------------------------------------------------------
class TestGripNanSafety:
    def test_nan_polish_falls_back(self, monkeypatch):
        _force_builder_on()
        _clear_grip_env()
        baseline = _xyz_canonical(_isomers(PANEL[0]))

        from delfin.fffree import grip_polish as _gp_mod

        def _nan_polish(P0, mol, metal, donors, geom="", mogul_lib=None, **kw):
            nan_R = np.asarray(P0, dtype=np.float64).copy()
            nan_R[:] = np.nan
            return nan_R

        monkeypatch.setattr(_gp_mod, "grip_polish", _nan_polish)
        os.environ["DELFIN_FFFREE_GRIP"] = "1"
        try:
            res = _isomers(PANEL[0])
        finally:
            _clear_grip_env()
        assert res is not None
        for xyz, _ in res:
            _, P = _parse_xyz(xyz)
            assert np.all(np.isfinite(P)), (
                "NaN from grip_polish must not propagate to emitted XYZ"
            )
        # And the entire emitted set should equal the GRIP-OFF baseline:
        # the hook fully fell back to P, so downstream is unchanged.
        assert _xyz_canonical(res) == baseline


# ---------------------------------------------------------------------------
# 5. None-safety — contrived grip_polish returning None
# ---------------------------------------------------------------------------
class TestGripNoneSafety:
    def test_none_polish_falls_back(self, monkeypatch):
        _force_builder_on()
        _clear_grip_env()
        baseline = _xyz_canonical(_isomers(PANEL[0]))

        from delfin.fffree import grip_polish as _gp_mod

        def _none_polish(P0, mol, metal, donors, geom="", mogul_lib=None, **kw):
            return None

        monkeypatch.setattr(_gp_mod, "grip_polish", _none_polish)
        os.environ["DELFIN_FFFREE_GRIP"] = "1"
        try:
            res = _isomers(PANEL[0])
        finally:
            _clear_grip_env()
        assert res is not None
        assert _xyz_canonical(res) == baseline


# ---------------------------------------------------------------------------
# 6. GripLibrary.get_default() singleton + missing-path semantics
# ---------------------------------------------------------------------------
class TestGripLibrarySingleton:
    def test_singleton_identity(self):
        _clear_grip_env()
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        a = GripLibrary.get_default()
        b = GripLibrary.get_default()
        if a is None:
            assert b is None
        else:
            assert a is b, "get_default() must be process-wide singleton"

    def test_missing_path_returns_none_quietly(self, tmp_path, caplog):
        from delfin.fffree.grip_mogul_lookup import GripLibrary

        bogus = tmp_path / "does_not_exist.npz"
        os.environ["DELFIN_GRIP_LIB_PATH"] = str(bogus)
        GripLibrary._MISSING_WARNED = False  # reset the warned-once flag
        try:
            with caplog.at_level(
                logging.WARNING, logger="delfin.fffree.grip_mogul_lookup"
            ):
                lib = GripLibrary.get_default()
        finally:
            _clear_grip_env()
            GripLibrary._MISSING_WARNED = False
        assert lib is None
        assert any(
            "not found" in r.getMessage().lower() for r in caplog.records
        ), "missing library must emit a WARN log line"

    def test_env_override_path(self, tmp_path, caplog):
        """When DELFIN_GRIP_LIB_PATH points to a real lib, get_default()
        uses it.  We don't have a second .npz handy, so we just check
        that pointing at the default-pinned path explicitly yields a
        library identical to the default lookup."""
        from delfin.fffree.grip_mogul_lookup import (
            GripLibrary,
            DEFAULT_LIB_PATH,
        )
        if not DEFAULT_LIB_PATH.exists():
            pytest.skip("default library not present in this checkout")

        _clear_grip_env()
        default = GripLibrary.get_default()
        os.environ["DELFIN_GRIP_LIB_PATH"] = str(DEFAULT_LIB_PATH)
        try:
            via_env = GripLibrary.get_default()
        finally:
            _clear_grip_env()
        assert via_env is default, (
            "env-override pointing at the default path must hit the cache"
        )
