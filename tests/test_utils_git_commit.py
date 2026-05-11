"""Tests for get_git_commit_info() resolution chain.

Compute nodes often run DELFIN from a pip-installed copy without a .git
tree, so the git fallback returns None and the run banner reports
"Git commit: unknown". The env-var override (set by the SLURM wrapper
at submit time) and the baked _commit.txt fallback fix that.
"""
from pathlib import Path

from delfin.utils import get_git_commit_info


def test_env_var_override_wins(monkeypatch):
    monkeypatch.setenv("DELFIN_GIT_COMMIT", "abc1234")
    assert get_git_commit_info() == "abc1234"


def test_env_var_override_with_dirty_marker(monkeypatch):
    monkeypatch.setenv("DELFIN_GIT_COMMIT", "abc1234-dirty")
    assert get_git_commit_info() == "abc1234-dirty"


def test_env_var_empty_falls_through(monkeypatch):
    """Empty env var must not short-circuit the resolution chain."""
    monkeypatch.setenv("DELFIN_GIT_COMMIT", "")
    result = get_git_commit_info()
    assert result is None or len(result) >= 7


def test_commit_file_fallback(monkeypatch, tmp_path):
    monkeypatch.delenv("DELFIN_GIT_COMMIT", raising=False)
    fake_commit_file = Path(__file__).parent.parent / "delfin" / "_commit.txt"
    if fake_commit_file.exists():
        return
    fake_commit_file.write_text("baked99\n")
    try:
        monkeypatch.setattr(
            "delfin.utils.subprocess.run",
            lambda *a, **kw: (_ for _ in ()).throw(FileNotFoundError("no git")),
        )
        assert get_git_commit_info() == "baked99"
    finally:
        fake_commit_file.unlink()
