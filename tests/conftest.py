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


@pytest.fixture(autouse=True)
def _isolate_user_wide_memory(tmp_path, monkeypatch):
    """Keep the user-scoped stores out of the real ``~/.delfin``.

    ``tmp_path`` bounds the PROJECT store, because that one is derived from
    the repo root a test passes in. The two USER-scoped stores are not: they
    resolve through ``Path.home()`` and ignore whatever root the test used.

    Both were observed writing into the real home directory during an
    ordinary run. The office registry had collected 178 entries, 176 of them
    reprs of mock objects a test had passed where a path was expected — and
    that registry is one of the carriers of the folder lock, so its contents
    decide which directory counts as an office workspace. The repository
    checkout itself was among them, which made a full-suite run treat DELFIN
    as an office folder and fail one memory test that passed in isolation.

    The redirect steps aside for a test that has pointed ``Path.home`` at a
    directory of its own: those tests were already isolated, by the more
    direct route, and silently overriding them would move the store out from
    under their assertions.
    """
    import pathlib

    from delfin.agent import memory_store as ms
    real_home = pathlib.Path.home()
    home = tmp_path / "user_home_delfin"

    def _redirected(name: str) -> pathlib.Path:
        current = pathlib.Path.home()
        if current != real_home:
            return current / ".delfin" / name
        return home / name

    monkeypatch.setattr(
        ms, "_delfin_global_memory_dir", lambda: _redirected("memory"))
    monkeypatch.setattr(
        ms, "_office_ws_file",
        lambda: _redirected("office_workspaces.json"))
    # The registry caches its parsed contents in a module global, so a stale
    # cache from an earlier test would outlive the redirect.
    monkeypatch.setattr(ms, "_office_ws_cache", None)
    yield
    monkeypatch.setattr(ms, "_office_ws_cache", None)


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


@pytest.fixture(autouse=True)
def _isolate_agent_telemetry(monkeypatch, tmp_path):
    """Redirect the agent's per-user telemetry sinks into the test tmp dir.

    Engine fixtures with mocked clients still execute the real
    turn-metrics/tool-trace/failure-log writers at turn end; without this
    every test run leaves mock records in the user's real ~/.delfin and
    pollutes the eval loop's health statistics. Tests that monkeypatch
    these constants themselves simply override this baseline."""
    from delfin.agent import failure_log as _fl
    from delfin.agent import tool_trace as _tt
    from delfin.agent import turn_metrics as _tm
    monkeypatch.setattr(_tm, "_DIR", tmp_path / "turn_metrics")
    monkeypatch.setattr(_tt, "_DIR", tmp_path / "tool_trace")
    monkeypatch.setattr(_fl, "_LOG_PATH", tmp_path / "failure_log.jsonl")
