"""Adversarial tests for package SX2 (reviewer nacht-s2).

Target 2 of package SX2: a calculation started inside the agent's command
sandbox (the "cage") puts its ORCA scratch on node-local storage that is
writable THERE (the session scratch), never on Lustre and never in the
workspace. Outside a cage nothing changes.

Real cage environment shape (the shape the product actually exports,
api_client._session_scratch_dir + the phase-1 px patch to
_scrubbed_bash_env): TMPDIR points at the session scratch AND the cage
exports DELFIN_SESSION_SCRATCH, gated on the cage marker
DELFIN_AGENT_SESSION. These tests use that stand-in shape; they never
touch a real cage, a real /scratch tree or the real home. They call the
single decision point, ``delfin.orca._ensure_orca_scratch_dir``, through
its public name and reset its module-level memo between cases.

History: the first version of this file (commit 321c41d44) tested the
hypothesis that a generic TMPDIR rule was the fix. The builder chose an
explicit DELFIN_SESSION_SCRATCH variable instead, because the gate
itself runs with the cage marker set and a generic TMPDIR rule flipped
three neighbouring scratch tests. The tests were rewritten to the real
shape; the old TMPDIR-only behaviour is covered by TestTmpdirOnlyLegacy.
"""

import os
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin import orca  # noqa: E402

CAGE_KEYS = ("ORCA_SCRDIR", "ORCA_TMPDIR", "DELFIN_SCRATCH",
             "SLURM_TMPDIR", "SLURM_JOB_ID", "TMPDIR", "TMP", "TEMP",
             "DELFIN_SESSION_SCRATCH", "DELFIN_AGENT_SESSION")


@pytest.fixture
def fresh_scratch_decision(monkeypatch):
    """Reset the memo and a clean environment around each decision."""
    monkeypatch.setattr(orca, "_RUN_SCRATCH_DIR", None)
    for key in CAGE_KEYS:
        monkeypatch.delenv(key, raising=False)
    return monkeypatch


@pytest.fixture
def in_cage(fresh_scratch_decision):
    """Stand-in for the real cage environment: marker + session scratch,
    TMPDIR alongside (what _scrubbed_bash_env exports after the px patch)."""
    monkeypatch = fresh_scratch_decision
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "cage")

    def set_session_scratch(path: Path) -> None:
        monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(path))
        monkeypatch.setenv("TMPDIR", str(path))

    return set_session_scratch


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
    """Inside the cage the session scratch must win over derived SLURM paths."""

    def test_session_scratch_wins_over_readonly_slurm_tmpdir(self, in_cage, tmp_path, monkeypatch):
        """The measured JC failure mode: SLURM_TMPDIR is read-only inside
        the cage. The writable session scratch must be chosen without
        relying only on the OSError fallthrough of mkdir-time."""
        readonly_base = tmp_path / "slurm_tmpdir_readonly"
        readonly_base.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        monkeypatch.setenv("SLURM_TMPDIR", str(readonly_base))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_session_scratch_wins_over_readonly_delfin_scratch(self, in_cage, tmp_path, monkeypatch):
        """A DELFIN_SCRATCH that cannot take the run dir (read-only mount,
        the Lustre case the cluster operators asked us to avoid) must not
        keep the calculation from a writable node-local scratch."""
        readonly_lustre = tmp_path / "lustre_scratch"
        readonly_lustre.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        monkeypatch.setenv("DELFIN_SCRATCH", str(readonly_lustre))
        readonly_lustre.chmod(0o555)
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_session_scratch_beats_slurm_tmpdir_when_delfin_scratch_unset(self, in_cage, tmp_path, monkeypatch):
        monkeypatch.setenv("SLURM_TMPDIR", str(tmp_path / "job_tmp"))
        (tmp_path / "job_tmp").mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_session_scratch_without_marker_is_ignored(self, fresh_scratch_decision, tmp_path, monkeypatch):
        """Outside a cage the variable does not exist; a stale one must
        not change the classic order (the gate runs WITH the marker, so
        the gating must read the marker, not the variable)."""
        delfin_scratch = tmp_path / "delfin_scratch"
        delfin_scratch.mkdir()
        stale = tmp_path / "stale_session_scratch"
        stale.mkdir()
        monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(stale))
        monkeypatch.setenv("TMPDIR", str(stale))
        monkeypatch.setenv("DELFIN_SCRATCH", str(delfin_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(delfin_scratch), tmp_path)

    def test_marker_with_whitespace_and_empty_variable(self, in_cage, tmp_path, monkeypatch):
        """Adversarial env shapes: the marker is set but the scratch
        variable is empty -- the empty candidate must be skipped, not
        crash or produce an empty base path."""
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        monkeypatch.setenv("DELFIN_AGENT_SESSION", "  cage  ")
        monkeypatch.setenv("DELFIN_SESSION_SCRATCH", "")
        monkeypatch.setenv("TMPDIR", str(session_scratch))
        result = orca._ensure_orca_scratch_dir()
        assert result is not None
        assert str(result), "empty result path"


class TestScratchOutsideCage:
    """Outside a cage nothing changes: the historical order holds."""

    def test_no_tmpdir_env_delfin_scratch_still_used(self, fresh_scratch_decision, tmp_path):
        monkeypatch = fresh_scratch_decision
        delfin_scratch = tmp_path / "delfin_scratch"
        delfin_scratch.mkdir()
        monkeypatch.setenv("DELFIN_SCRATCH", str(delfin_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(delfin_scratch), tmp_path)

    def test_no_cage_variables_at_all_falls_back_to_tempdir(self, fresh_scratch_decision, tmp_path, monkeypatch):
        """Historical outside-cage behavior: with nothing else set, the
        system temp dir (as TMPDIR names it) is the base. The conftest
        pins ``tempfile.tempdir`` to the suite's base temp; undo that pin
        so the test measures the product's env reading, not the suite's."""
        system_temp = tmp_path / "system_temp"
        system_temp.mkdir()
        monkeypatch.setenv("TMPDIR", str(system_temp))
        monkeypatch.setattr(tempfile, "tempdir", str(system_temp))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(system_temp), tmp_path, in_fallback=True)

    def test_decision_is_memoized(self, in_cage, tmp_path, monkeypatch):
        """The choice is made once (single decision point, orca.py:375-377)."""
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        first = orca._ensure_orca_scratch_dir()
        monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(tmp_path / "other_place"))
        second = orca._ensure_orca_scratch_dir()
        assert first == second


class TestTmpdirOnlyLegacy:
    """The original F2 finding, kept under the real variable: when the
    cage marks the session only through TMPDIR (no DELFIN_SESSION_SCRATCH,
    as before the px patch lands) nothing regresses -- SLURM_TMPDIR is
    still preferred, unchanged from main."""

    def test_tmpdir_alone_does_not_flip_the_order(self, fresh_scratch_decision, tmp_path, monkeypatch):
        job_tmp = tmp_path / "job_tmp"
        job_tmp.mkdir()
        monkeypatch.setenv("SLURM_TMPDIR", str(job_tmp))
        monkeypatch.setenv("TMPDIR", str(tmp_path / "session_scratch"))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(job_tmp), tmp_path)


class TestExplicitChoicesKeepPriority:
    """Explicit operator choices stay above the derived session scratch
    (the builder's design, contested for the writable-Lustre case --
    QS question recorded in .gate/HANDOFF.md)."""

    def test_explicit_orca_scrdir_still_wins_in_cage(self, in_cage, tmp_path, monkeypatch):
        explicit = tmp_path / "explicit_orca_scrdir"
        explicit.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        monkeypatch.setenv("ORCA_SCRDIR", str(explicit))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(explicit), tmp_path)

    def test_writable_delfin_scratch_keeps_first_place(self, in_cage, tmp_path, monkeypatch):
        """The QS-disputed point: a WRITABLE explicit DELFIN_SCRATCH (here
        standing for an operator-chosen Lustre path) still wins over the
        session scratch. Pinned as the builder's design; if QS rules that
        node-local must always win, this test is rewritten along its
        intent, and the reason goes into the commit message."""
        delfin_scratch = tmp_path / "explicit_delfin_scratch"
        delfin_scratch.mkdir()
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        monkeypatch.setenv("DELFIN_SCRATCH", str(delfin_scratch))
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(delfin_scratch), tmp_path)


class TestAdversarialEnvironmentShapes:
    """Abbreviated and adversarial environment values on the real variable."""

    def test_empty_string_candidates_are_skipped(self, in_cage, tmp_path, monkeypatch):
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        monkeypatch.setenv("DELFIN_SCRATCH", "")
        monkeypatch.setenv("ORCA_SCRDIR", "")
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_unicode_path_in_session_scratch(self, in_cage, tmp_path):
        session_scratch = tmp_path / "sitzung_scratch_ünïcode"
        session_scratch.mkdir()
        in_cage(session_scratch)
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_session_scratch_with_trailing_slash(self, in_cage, tmp_path, monkeypatch):
        session_scratch = tmp_path / "session_scratch"
        session_scratch.mkdir()
        in_cage(session_scratch)
        monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(session_scratch) + os.sep)
        result = orca._ensure_orca_scratch_dir()
        _assert_run_dir_under(result, str(session_scratch), tmp_path)

    def test_session_scratch_path_that_does_not_exist_yet(self, in_cage, tmp_path, monkeypatch):
        """A session scratch that is not there (fresh cage before the
        first command created it) must not crash the decision."""
        session_scratch = tmp_path / "not_yet_created"
        in_cage(session_scratch)
        monkeypatch.delenv("TMPDIR", raising=False)
        result = orca._ensure_orca_scratch_dir()
        assert result is not None
