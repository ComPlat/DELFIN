"""Unit tests for delfin.classify + delfin.ensemble_router.

Phase 2 of Hybrid-Path (nature_project/15_HYBRID_PATH_FINAL.md).
"""
from __future__ import annotations

import pytest

from delfin.classify import (
    BLOCK_3D, BLOCK_4D, BLOCK_5D, BLOCK_LN, BLOCK_NONE,
    COORD_SIGMA, COORD_HAPTO, COORD_MULTI_SIGMA, COORD_MULTI_HAPTO,
    COORD_NO_METAL, COORD_CLUSTER,
    ClassFeatures, classify_smiles,
)
from delfin.ensemble_router import (
    ROUTING_TABLE, lookup_specialist_id, route, route_or_fallback,
    routing_summary,
)


# ----------------------------------------------------------------------
# Classifier — coord_class detection
# ----------------------------------------------------------------------

class TestCoordClass:
    def test_sigma_mono_fecl4(self):
        f = classify_smiles("[Fe-2]([Cl])([Cl])([Cl])[Cl]")
        assert f.coord_class == COORD_SIGMA
        assert f.n_metals == 1

    def test_sigma_pt_amine(self):
        f = classify_smiles("[Cl][Pt](Cl)([NH3+])[NH3+]")
        assert f.coord_class == COORD_SIGMA

    def test_no_metal_organic(self):
        f = classify_smiles("CC(=O)O")  # acetic acid
        assert f.coord_class == COORD_NO_METAL
        assert f.n_metals == 0

    def test_multi_sigma_two_metals(self):
        f = classify_smiles("[Mn]([C]#O)([C]#O)[Mn]([C]#O)([C]#O)")
        assert f.coord_class == COORD_MULTI_SIGMA
        assert f.n_metals == 2

    def test_cluster_three_metals(self):
        f = classify_smiles("[Fe]1[Fe]2[Fe]3[S]12 [S]23")
        assert f.coord_class == COORD_CLUSTER
        assert f.n_metals == 3


# ----------------------------------------------------------------------
# Classifier — metal_block detection
# ----------------------------------------------------------------------

class TestMetalBlock:
    @pytest.mark.parametrize("smi,expected_block", [
        ("[Fe][Cl]", BLOCK_3D),
        ("[Cu][Cl]", BLOCK_3D),
        ("[Mo](C#O)(C#O)C#O", BLOCK_4D),
        ("[Ru]([Cl])[Cl]", BLOCK_4D),
        ("[Pt](Cl)Cl", BLOCK_5D),
        ("[W](C#O)(C#O)C#O", BLOCK_5D),
        ("[Eu](Cl)(Cl)Cl", BLOCK_LN),
        ("CC(=O)O", BLOCK_NONE),
    ])
    def test_block_detection(self, smi, expected_block):
        f = classify_smiles(smi)
        assert f.metal_block == expected_block, f"{smi}: got {f.metal_block}"


# ----------------------------------------------------------------------
# Determinism + purity
# ----------------------------------------------------------------------

class TestDeterminism:
    def test_same_smiles_same_features(self):
        smi = "[Fe-2]([Cl])([Cl])([Cl])[Cl]"
        f1 = classify_smiles(smi)
        f2 = classify_smiles(smi)
        assert f1 == f2

    def test_immutable(self):
        f = classify_smiles("[Fe][Cl]")
        with pytest.raises(Exception):
            f.coord_class = "hapto"  # type: ignore[misc]


# ----------------------------------------------------------------------
# Router — lookup_specialist_id
# ----------------------------------------------------------------------

class TestRouterLookup:
    def test_multi_sigma_dominated_by_123a130(self):
        for block in [BLOCK_3D, BLOCK_4D, BLOCK_5D, "p"]:
            f = ClassFeatures(
                coord_class=COORD_MULTI_SIGMA, metal_block=block,
                metals=("Fe","Fe"), n_metals=2, has_hapto=False,
                has_aromatic=False, has_bridge=False,
            )
            spec = lookup_specialist_id(f)
            assert spec == "multi_sigma_123a130", f"{block}: got {spec}"

    def test_sigma_3d_routes_to_1e7eefe_prev(self):
        f = ClassFeatures(COORD_SIGMA, BLOCK_3D, ("Fe",), 1,
                          False, False, False)
        assert lookup_specialist_id(f) == "sigma_1e7eefe_prev"

    def test_sigma_4d_routes_to_229e5dc(self):
        f = ClassFeatures(COORD_SIGMA, BLOCK_4D, ("Mo",), 1,
                          False, False, False)
        assert lookup_specialist_id(f) == "sigma_229e5dc"

    def test_unrouted_returns_none(self):
        f = ClassFeatures(COORD_SIGMA, "s", ("Na",), 1,
                          False, False, False)
        # s-block not in routing table
        assert lookup_specialist_id(f) is None


# ----------------------------------------------------------------------
# Router — env-flag-gated route_or_fallback
# ----------------------------------------------------------------------

class TestRouterEnvFlag:
    def _fallback(self, smiles, **kw):
        return f"FALLBACK-{smiles}", None

    def test_mode_0_always_fallback(self, monkeypatch):
        monkeypatch.delenv("DELFIN_ENSEMBLE_ROUTER", raising=False)
        xyz, err = route_or_fallback(
            "[Fe-2]([Cl])([Cl])([Cl])[Cl]", self._fallback,
        )
        assert xyz.startswith("FALLBACK-")

    def test_mode_1_falls_back_when_no_route(self, monkeypatch):
        monkeypatch.setenv("DELFIN_ENSEMBLE_ROUTER", "1")
        # s-block sodium not in routing table → must fall back
        xyz, err = route_or_fallback("[Na+]", self._fallback)
        assert xyz is not None
        assert xyz.startswith("FALLBACK-")

    def test_mode_2_strict_no_fallback(self, monkeypatch):
        monkeypatch.setenv("DELFIN_ENSEMBLE_ROUTER", "2")
        # No specialists registered → strict mode returns error
        xyz, err = route_or_fallback("[Na+]", self._fallback)
        assert xyz is None
        assert err is not None
        assert "ensemble_router" in err


# ----------------------------------------------------------------------
# Routing summary diagnostic
# ----------------------------------------------------------------------

class TestRoutingSummary:
    def test_summary_includes_all_routes(self):
        s = routing_summary()
        assert s["n_routes_defined"] == len(ROUTING_TABLE)
        assert len(s["routes"]) == len(ROUTING_TABLE)

    def test_summary_includes_unregistered_marker(self):
        s = routing_summary()
        # Before any specialist registration, all routes are unregistered
        for r in s["routes"]:
            assert "registered" in r


# ----------------------------------------------------------------------
# End-to-end: gold-standard SMILES through classifier
# ----------------------------------------------------------------------

class TestChampionSpecialist:
    """Tests for ChampionSpecialist + registry (Phase 2B)."""

    def test_registry_register_all(self):
        from delfin.class_modules.registry import (
            register_all_specialists, _INITIALIZED,
        )
        # Idempotent
        n1 = register_all_specialists()
        n2 = register_all_specialists()
        # First call registers (or 0 if already done), second always 0
        assert n2 == 0

    def test_specialist_protocol_fields(self):
        from delfin.class_modules.registry import register_all_specialists
        from delfin.class_modules import get_specialist
        register_all_specialists()
        spec = get_specialist("multi_sigma_123a130")
        assert spec is not None
        assert spec.id == "multi_sigma_123a130"
        assert spec.source_commit == "123a130"
        assert hasattr(spec, "convert")

    def test_specialist_is_available_when_worktree_exists(self):
        from delfin.class_modules.registry import register_all_specialists
        from delfin.class_modules import get_specialist
        register_all_specialists()
        spec = get_specialist("multi_sigma_123a130")
        # Worktree at /home/qmchem_max/agent_workspace/commit_sweep/worktrees/123a130
        # should exist
        assert spec.is_available()

    def test_parse_driver_output_well_formed(self):
        from delfin.class_modules.champion_specialist import ChampionSpecialist
        out = (
            "---XYZ-START---\n"
            "Fe 0.0 0.0 0.0\n"
            "Cl 2.3 0.0 0.0\n"
            "---XYZ-END---\n"
            "---ERR---\n"
            ""
        )
        xyz, err = ChampionSpecialist._parse_driver_output(out)
        assert xyz is not None
        assert "Fe" in xyz
        assert err is None

    def test_parse_driver_output_error_path(self):
        from delfin.class_modules.champion_specialist import ChampionSpecialist
        out = (
            "---XYZ-START---\n"
            "\n"
            "---XYZ-END---\n"
            "---ERR---\n"
            "RDKit failed: invalid SMILES"
        )
        xyz, err = ChampionSpecialist._parse_driver_output(out)
        assert xyz is None
        assert err is not None
        assert "invalid" in err

    def test_archive_readback_finds_existing_xyz(self, tmp_path):
        from delfin.class_modules.champion_specialist import ChampionSpecialist
        # Set up fake archive structure
        archive_base = tmp_path
        archive_sub = archive_base / "test_archive"
        archive_sub.mkdir()
        (archive_sub / "MOL1.xyz").write_text("3\nheader\nFe 0 0 0\nCl 2 0 0\nCl -2 0 0\n")

        spec = ChampionSpecialist(
            id="test",
            source_commit="abc1234",
            archive_subdir="test_archive",
            archive_base=archive_base,
        )
        cached = spec._archive_lookup("MOL1")
        assert cached is not None
        assert "Fe" in cached

    def test_archive_readback_returns_none_for_missing(self, tmp_path):
        from delfin.class_modules.champion_specialist import ChampionSpecialist
        spec = ChampionSpecialist(
            id="test",
            source_commit="abc1234",
            archive_subdir="missing_archive",
            archive_base=tmp_path,
        )
        assert spec._archive_lookup("NONE") is None
        assert spec._archive_lookup("") is None
        assert spec._archive_lookup(None) is None


class TestEndToEndRouting:
    """End-to-end: router + specialist + legacy fallback."""

    def test_default_off_uses_fallback(self, monkeypatch):
        from delfin.ensemble_router import route_or_fallback
        from delfin.smiles_converter import smiles_to_xyz
        monkeypatch.delenv("DELFIN_ENSEMBLE_ROUTER", raising=False)
        xyz_router, _ = route_or_fallback(
            "[Fe-2]([Cl])([Cl])([Cl])[Cl]", smiles_to_xyz,
        )
        xyz_legacy, _ = smiles_to_xyz("[Fe-2]([Cl])([Cl])([Cl])[Cl]")
        assert xyz_router == xyz_legacy  # byte-identical when mode=0


class TestGoldStandard:
    @pytest.mark.parametrize("smi,exp_coord,exp_block", [
        ("[Fe-2]([Cl])([Cl])([Cl])[Cl]",        COORD_SIGMA,       BLOCK_3D),
        ("[Au-]([Cl])([Cl])([Cl])[Cl]",         COORD_SIGMA,       BLOCK_5D),
        ("[Co-2]([Br])([Br])([Br])[Br]",        COORD_SIGMA,       BLOCK_3D),
        ("[Cd-2]([Cl])([Cl])([Cl])[Cl]",        COORD_SIGMA,       BLOCK_4D),
        ("[Mo+2]([S-])([S-])([S-])[S-]",        COORD_SIGMA,       BLOCK_4D),
        ("[Cr-3]([Cl])([Cl])([Cl])([Cl])([Cl])[Cl]", COORD_SIGMA,  BLOCK_3D),
    ])
    def test_tier1_gold_standard_classify(self, smi, exp_coord, exp_block):
        f = classify_smiles(smi)
        assert f.coord_class == exp_coord, f"{smi}: coord got {f.coord_class}"
        assert f.metal_block == exp_block, f"{smi}: block got {f.metal_block}"
