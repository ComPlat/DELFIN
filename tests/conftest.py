"""Shared test fixtures.

Redirect the subagent state files (live registry, telemetry, finished
sessions) away from the real ``~/.delfin`` so test runs never leave
artifacts that show up in the user's dashboard (live panel, ``/agents``
listing, stats).
"""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _isolate_subagent_state(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_RUNNING_DIR",
                        tmp_path / "subagent_running")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH",
                        tmp_path / "subagent_telemetry.jsonl")
    monkeypatch.setattr(sa, "_SESSIONS_DIR",
                        tmp_path / "subagent_sessions")
    yield


# Secret API keys delfin resolves, each from the env OR ~/.delfin/credentials.json.
_SECRET_KEY_VARS = ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY")


@pytest.fixture(autouse=True)
def _ci_parity_no_secrets(tmp_path, monkeypatch):
    """Make the local test environment match CI's: NO API key reachable from
    either the environment or the ``~/.delfin/credentials.json`` store.

    CI runs key-free and fully mocked (see ``.github/workflows/ci.yml`` —
    "no secrets are required"). A test that quietly relies on an ambient key
    therefore passes on a developer box that has one, but fails in CI — the
    exact "lokal grün, CI rot" trap. Stripping the keys here surfaces that
    whole class of failure on the developer's own machine, before the push.

    A test that genuinely needs a key still sets it itself (this autouse
    fixture runs first; the test's own ``monkeypatch.setenv`` then wins).
    """
    for var in _SECRET_KEY_VARS:
        monkeypatch.delenv(var, raising=False)
    try:
        from delfin.agent import credentials as _cred
        monkeypatch.setattr(_cred, "_DEFAULT_PATH",
                            tmp_path / "no_credentials.json")
    except Exception:
        pass
    yield


# ---------------------------------------------------------------------------
# CI tiering (2026-06-27). The fast gate (`pytest -m "not slow"`) is the merge
# check and must stay ~minutes. Two kinds of test are auto-marked `slow` and
# moved to the non-blocking slow-tests job:
#   (a) heavy MANTA build/integration tests (construct full complexes; the
#       slowest single test is ~130s — together they push the suite >60 min),
#   (b) ORDER-DEPENDENT tests that only pass in the full serial context (they
#       rely on module-level state another test sets up — sys.modules-driven
#       scans, runtime-mutated registries); in any subset they false-fail.
# Real fix for (b) = test isolation (tracked debt); until then they run in the
# (non-blocking) slow job so they neither block nor false-fail the gate.
# Add new build-heavy / order-dependent test files to this set.
# ---------------------------------------------------------------------------
_SLOW_TEST_FILES = frozenset({
    # (a) build/integration
    "test_all_features_contracts.py", "test_classify_and_router.py",
    "test_crest_censo.py", "test_decompose_kekulize_split.py",
    "test_embed_determinism.py", "test_equatorial_square_kappa4.py",
    "test_eta_ring_hydrogens.py", "test_fffree_backbone_reembed.py",
    "test_fffree_chelate_backbone.py", "test_fffree_conformer_complete.py",
    "test_fffree_conformer_coverage.py", "test_fffree_conformer_seating.py",
    "test_fffree_diatomic_orient.py", "test_fffree_hapto_recall.py",
    "test_fffree_interlig_vdw_gate.py", "test_fffree_joint_declash.py",
    "test_fffree_metalloid_donor.py", "test_fffree_nhc_carbene.py",
    "test_highcn_coverage.py", "test_isomer_benchmark.py",
    "test_isomer_enum_determinism.py", "test_ligand_rigid_seating.py",
    "test_metal_isomer_determinism.py", "test_multihapto_patches.py",
    "test_multi_sigma_path_v2.py", "test_platform.py",
    "test_ring_bounds_embed.py", "test_smiles_converter_regressions.py",
    "test_tool_contracts.py", "test_uff_soft_determinism.py",
    "test_user_smiles_suite.py",
    "test_welle5l_t3b_daqiwaz_pyridine_coverage.py",
    # (b) order-dependent (pass only in the full serial context)
    "test_no_regression_undefined_names.py", "test_post_optimizer.py",
    "test_theorem_d_asymmetric.py", "test_welle5p_c_hotfix.py",
})


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: heavy or order-dependent test, excluded from the fast CI gate",
    )


def pytest_collection_modifyitems(config, items):
    for item in items:
        fname = item.nodeid.split("::", 1)[0].rsplit("/", 1)[-1]
        if fname in _SLOW_TEST_FILES:
            item.add_marker(pytest.mark.slow)
