"""Phase B: central MLIP model-cache environment handling.

Covers the priority order, writability checks, the read-only-home error
and the three integration points (factory, sbatch rendering, launcher).
All tests run in temporary directories; no model weights, no network.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.mlp_tools.env_setup import (
    ensure_mlp_cache_env,
    resolve_mlp_cache_dir,
)

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "delfin" / "mlp_tools" / "bin" / "delfin-mlp-extopt"


def test_priority_env_var_wins(tmp_path, monkeypatch):
    monkeypatch.setenv("AIMNET_CACHE_DIR", str(tmp_path / "a"))
    monkeypatch.setenv("DELFIN_MLP_CACHE_DIR", str(tmp_path / "b"))
    monkeypatch.delenv("HOME")  # make the home fallback unreachable
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "nonexistent-home")
    cache, source = resolve_mlp_cache_dir(os.environ)
    assert (cache, source) == (str(tmp_path / "a"), "AIMNET_CACHE_DIR")


def test_priority_delfin_key_over_home_fallback(tmp_path, monkeypatch):
    monkeypatch.delenv("AIMNET_CACHE_DIR", raising=False)
    monkeypatch.setenv("DELFIN_MLP_CACHE_DIR", str(tmp_path / "delfin-cache"))
    cache, source = resolve_mlp_cache_dir(os.environ)
    assert source == "DELFIN_MLP_CACHE_DIR"


def test_default_is_home_cache_aimnet(monkeypatch):
    monkeypatch.delenv("AIMNET_CACHE_DIR", raising=False)
    monkeypatch.delenv("DELFIN_MLP_CACHE_DIR", raising=False)
    monkeypatch.setattr(Path, "home", lambda: Path("/nonexistent-home"))
    cache, source = resolve_mlp_cache_dir({})
    assert source == "default"
    assert cache == "/nonexistent-home/.cache/aimnet"


def test_ensure_creates_and_exports(tmp_path, monkeypatch):
    monkeypatch.delenv("AIMNET_CACHE_DIR", raising=False)
    target = tmp_path / "fresh-cache"
    cache = ensure_mlp_cache_env({"AIMNET_CACHE_DIR": str(target)})
    assert cache == str(target)
    assert target.is_dir()


def test_readonly_cache_gives_facing_error(tmp_path):
    ro = tmp_path / "ro-cache"
    ro.mkdir()
    ro.chmod(0o500)
    if os.access(ro, os.W_OK):  # root ignores modes — skip honestly
        pytest.skip("running as a user that ignores read-only modes")
    with pytest.raises(RuntimeError, match="AIMNET_CACHE_DIR"):
        ensure_mlp_cache_env({"AIMNET_CACHE_DIR": str(ro)})


def test_cache_path_is_a_file_not_a_dir(tmp_path):
    victim = tmp_path / "not-a-dir"
    victim.write_text("x")
    with pytest.raises(RuntimeError, match="not a directory"):
        ensure_mlp_cache_env({"AIMNET_CACHE_DIR": str(victim)})


def test_uncreatable_cache_names_the_override(tmp_path):
    blocked = tmp_path / "blocked" / "deeper"
    (tmp_path / "blocked").write_text("file in the way")
    with pytest.raises(RuntimeError, match="AIMNET_CACHE_DIR"):
        ensure_mlp_cache_env({"AIMNET_CACHE_DIR": str(blocked)})


def test_sbatch_renders_the_resolved_cache(monkeypatch, tmp_path):
    from delfin.mlp_tools.cluster_dispatch import (
        MLPClusterConfig,
        render_mlp_sbatch,
    )

    cfg = MLPClusterConfig(gpu_partitions="gpu", gres="gpu:1")
    monkeypatch.setenv("AIMNET_CACHE_DIR", str(tmp_path / "job-cache"))
    script = render_mlp_sbatch(
        cfg, job_dir=tmp_path, xyz=tmp_path / "s.xyz", backend="aimnet2"
    )
    assert f"AIMNET_CACHE_DIR={tmp_path / 'job-cache'}" in script


def test_launcher_sets_the_home_fallback_only(tmp_path, monkeypatch):
    """The launcher exports the home fallback when nothing is configured
    and never overrides an explicit value."""
    monkeypatch.delenv("AIMNET_CACHE_DIR", raising=False)
    text = LAUNCHER.read_text(encoding="utf-8")
    assert "mlip-envs/aimnet-cache" in text
    assert '[ -z "${AIMNET_CACHE_DIR:-}" ]' in text
    mode = LAUNCHER.stat().st_mode
    assert mode & stat.S_IXUSR


def test_launcher_bash_priority_matches_python(tmp_path):
    """bin/delfin-mlp-extopt must honour DELFIN_MLP_CACHE_DIR before the
    home fallback -- priority parity with env_setup.py (Phase B2)."""
    launcher = REPO / "delfin" / "mlp_tools" / "bin" / "delfin-mlp-extopt"
    text = launcher.read_text(encoding="utf-8")
    # ordering: the DELFIN key branch must appear before the home dir test
    delfin_key_at = text.index("DELFIN_MLP_CACHE_DIR")
    home_at = text.index("mlip-envs/aimnet-cache")
    assert delfin_key_at < home_at
    # and the whole block must sit behind the explicit-override guard
    guard = text.index('if [ -z "${AIMNET_CACHE_DIR:-}" ]; then')
    assert guard < delfin_key_at


def test_installer_copies_the_repo_launcher(tmp_path):
    """install_mlp_tools.sh must copy the checked-in launcher (single
    source of truth) instead of a diverging template (Phase B2)."""
    launcher = REPO / "delfin" / "mlp_tools" / "bin" / "delfin-mlp-extopt"
    installer_text = (
        REPO / "delfin" / "mlp_tools" / "install_mlp_tools.sh"
    ).read_text(encoding="utf-8")
    assert "cp \"${_src_launcher}\"" in installer_text
    # the generated fallback template must carry the same cache priority
    assert installer_text.count("DELFIN_MLP_CACHE_DIR") >= 2
