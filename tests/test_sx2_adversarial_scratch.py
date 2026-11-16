"""Adversarial tests for package SX2 (reviewer nacht-s2).

Target 2 of package SX2: a calculation started inside the agent's command
sandbox (the "cage") puts its ORCA scratch on node-local storage that is
writable THERE (the session scratch, exported as TMPDIR), never on Lustre
and never in the workspace. Outside a cage nothing changes.

These tests use stand-in environments only: they never touch a real cage,
a real /scratch tree or the real home. They call the single decision
point, ``delfin.orca._ensure_orca_scratch_dir``, through its public name
and reset its module-level memo between cases.

Status on main @ a7e22dcda: RED — the cage's TMPDIR is not among
``base_candidates`` (orca.py:387-393), so with DELFIN_SCRATCH or
SLURM_TMPDIR present the scratch lands on the shared/Lustre path even
inside a cage where it is read-only or on shared storage.
"""

import os
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin import orca  # noqa: E402


@pytest.fixture
def fresh_scratch_decision(monkeypatch):
    """Reset the memo and a clean environment around each decision."""
    monkeypatch.setattr(orca, "_RUN_SCRATCH_DIR", None)
    for key in ("ORCA_SCRDIR", "ORCA_TMPDIR", "DELFIN_SCRATCH",
                "SLURM_TMPDIR", "SLURM_JOB_ID", "TMPDIR", "TMP", "TEMP"):
        monkeypatch.delenv(key, raising=False)
    return monkeypatch


def _assert_run_dir_under(result, prefix: str, tmp_path: Path,
                          in_fallback: bool = False) -> None:
    assert result is not None
    result_str = str(result)
    assert result_str.startswith(prefix), (
        f"scratch {result_str} is not under {prefix}")
    assert Path(result_str).exists(), f"{result_str} was not created"
    if in_fallback:
        # Only the last-resort path uses the per-user default name.
        assert "delfin_orca_scratch" in result_str


class TestScratchInsideCage:
    """Inside the cage the session scratch (TMPDIR) must win."""

    def test_tmpdir_writable_wins_over_delfin_scratch(self, fresh_scratch_decision, tmp_path):
        """DELFIN_SCRATCH here stands for a shared/Lustre path; the cage's
        TMPDIR stands for the node-local session scratch. In a cage the
        writable node-local one must be chosen."""
        monkeypatch = fresh_scratch_decision
        lustre_like = tmp_path / "lustre_workspace_scratch"
        lustre_like.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("DELFIN_SCRATCH", str(lustre_like))
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_tmpdir_wins_over_readonly_slurm_tmpdir(self, fresh_scratch_decision, tmp_path):
        """The measured JC failure mode: SLURM_TMPDIR is read-only inside
        the cage. The writable session scratch must be chosen without
        relying only on the OSError fallthrough of mkdir-time."""
        monkeypatch = fresh_scratch_decision
        readonly_base = tmp_path / "slurm_tmpdir_readonly"
        readonly_base.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("SLURM_TMPDIR", str(readonly_base))
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_explicit_orca_scrdir_still_wins_in_cage(self, fresh_scratch_decision, tmp_path):
        """ORCA_SCRDIR is an explicit operator choice and stays first, even
        in a cage — the node-local preference must not override an
        explicit, higher-priority candidate."""
        monkeypatch = fresh_scratch_decision
        explicit = tmp_path / "explicit_orca_scrdir"
        explicit.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("ORCA_SCRDIR", str(explicit))
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(explicit), tmp_path)


class TestScratchOutsideCage:
    """Outside a cage nothing changes: the historical order holds."""

    def test_no_tmpdir_env_delfin_scratch_still_used(self, fresh_scratch_decision, tmp_path):
        monkeypatch = fresh_scratch_decision
        delfin_scratch = tmp_path / "delfin_scratch"
        delfin_scratch.mkdir()
        monkeypatch.setenv("DELFIN_SCRATCH", str(delfin_scratch))
        monkeypatch.delenv("TMPDIR", raising=False)
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(delfin_scratch), tmp_path)

    def test_no_cage_variables_at_all_falls_back_to_tempdir(self, fresh_scratch_decision, tmp_path, monkeypatch):
        """Historical outside-cage behavior: with nothing else set, the
        system temp dir (as TMPDIR names it) is the base. The conftest
        pins ``tempfile.tempdir`` to the suite's base temp; undo that pin
        so the test measures the product's env reading, not the suite's."""
        monkeypatch = fresh_scratch_decision
        system_temp = tmp_path / "system_temp"
        system_temp.mkdir()
        monkeypatch.setenv("TMPDIR", str(system_temp))
        monkeypatch.setattr(tempfile, "tempdir", str(system_temp))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(system_temp), tmp_path, in_fallback=True)

    def test_decision_is_memoized(self, fresh_scratch_decision, tmp_path):
        """The choice is made once (single decision point, orca.py:375-377)."""
        monkeypatch = fresh_scratch_decision
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        first = orca._ensure_orca_scratch_dir()
        monkeypatch.setenv("TMPDIR", str(tmp_path / "other_place"))
        second = orca._ensure_orca_scratch_dir()
        assert first == second


class TestAdversarialEnvironmentShapes:
    """Abbreviated and adversarial environment values."""

    def test_empty_string_candidates_are_skipped(self, fresh_scratch_decision, tmp_path):
        monkeypatch = fresh_scratch_decision
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("DELFIN_SCRATCH", "")
        monkeypatch.setenv("ORCA_SCRDIR", "")
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_unicode_path_in_tmpdir(self, fresh_scratch_decision, tmp_path):
        monkeypatch = fresh_scratch_decision
        session_scratch = tmp_path / "sitzung_scratch_ünïcode"
        session_scratch.mkdir()
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_tmpdir_with_trailing_slash(self, fresh_scratch_decision, tmp_path):
        monkeypatch = fresh_scratch_decision
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("TMPDIR", str(session_scratch) + os.sep)
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)
