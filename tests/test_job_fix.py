"""Tests for the agent's failure-assessment + recalc-based fix bridge.

Safety-critical contract:
- SLURM OOM kill (no ORCA error in the .out) -> bounded CONTROL maxcore
  bump; OOM WITH an ORCA error -> defer to DELFIN's in-run recovery (#3).
- venv activation failure -> plain recalc resubmit (env race), retry-
  budgeted (#6); chemical / ORCA-internal -> defer; infra -> human.
- A proposed change NEVER touches chemistry/convergence keywords.
- One auto-retry per (job, signature); the second escalates.
- apply goes through DELFIN's recalc path (mode=delfin-recalc-classic),
  edits CONTROL.txt with a backup, never overwrites the old run.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import job_fix as jf


_CONTROL = """input_file=input.txt
charge=0
PAL=40
maxcore=9000
functional=B3LYP
"""


def _folder(tmp_path: Path, control: str = _CONTROL) -> Path:
    (tmp_path / "CONTROL.txt").write_text(control)
    return tmp_path


# --- CONTROL parsing + bump --------------------------------------------------

def test_parse_pal_maxcore_fallback(monkeypatch):
    # Force the light fallback (no RDKit import) and read the values.
    import builtins
    real_import = builtins.__import__

    def _no_dashboard(name, *a, **k):
        if name == "delfin.dashboard.input_processing":
            raise ImportError("blocked for test")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", _no_dashboard)
    assert jf.parse_pal_maxcore(_CONTROL) == (40, 9000)


def test_bump_maxcore_is_1_5x_only_that_line():
    new_text, old, new = jf.bump_maxcore(_CONTROL)
    assert old == "maxcore=9000"
    assert new == "maxcore=13500"          # ceil(9000 * 1.5)
    assert "maxcore=13500" in new_text
    assert "functional=B3LYP" in new_text  # nothing else touched
    assert "PAL=40" in new_text


def test_bump_maxcore_absent_returns_none():
    assert jf.bump_maxcore("PAL=40\ncharge=0\n") is None


# --- the chemistry guard (the line that must never be crossed) --------------

def test_chemistry_guard_blocks_convergence_keywords():
    assert jf._is_chemistry_safe("maxcore=9000", "maxcore=13500")
    assert not jf._is_chemistry_safe("! TightSCF", "! LooseSCF")
    assert not jf._is_chemistry_safe("TolE 1e-8", "TolE 1e-5")
    assert not jf._is_chemistry_safe("functional=B3LYP", "functional=BP86")


# --- routing ----------------------------------------------------------------

def test_oom_without_orca_error_prepares_memory_fix(tmp_path):
    # No .out files -> no ORCA error -> a pure SLURM kill the agent owns.
    a = jf.assess("12345", str(_folder(tmp_path)), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.SLURM_RESOURCE
    assert a.one_click is True
    assert a.fix.kind == "memory"
    assert a.fix.control_new == "maxcore=13500"


def test_oom_with_orca_error_defers_to_delfin(tmp_path, monkeypatch):
    # #3 disambiguation: ORCA itself reported an error -> in-run recovery
    # owns it, the agent must NOT bump maxcore.
    f = _folder(tmp_path)
    monkeypatch.setattr(jf, "_detect_orca_error", lambda folder: "memory_error")
    a = jf.assess("12345", str(f), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.DELFIN_AUTORECOVERS
    assert a.one_click is False


def test_timelimit_prepares_walltime_fix(tmp_path):
    a = jf.assess("7", str(_folder(tmp_path)), ["slurm-timelimit"],
                  attempts_path=tmp_path / "att.json")
    assert a.one_click is True
    assert a.fix.kind == "walltime"
    assert a.fix.new_time_limit == "48:00:00"


def test_venv_failure_is_retry_then_escalates(tmp_path):
    ap = tmp_path / "att.json"
    a1 = jf.assess("9", str(_folder(tmp_path)), ["venv-activation-failed"],
                   attempts_path=ap)
    assert a1.fix_class == jf.FixClass.ENV_RETRY
    assert a1.one_click is True and a1.fix.kind == "retry"
    jf.register_attempt("9", "venv-activation-failed", ap)
    a2 = jf.assess("9", str(_folder(tmp_path)), ["venv-activation-failed"],
                   attempts_path=ap)
    assert a2.escalate is True and a2.one_click is False


def test_scf_defers_to_delfin(tmp_path):
    a = jf.assess("9", str(tmp_path), ["scf-not-converged"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.DELFIN_AUTORECOVERS
    assert a.one_click is False
    assert "never loosen convergence" in a.recommendation.lower()


def test_infra_is_human_only(tmp_path):
    a = jf.assess("3", str(tmp_path), ["disk-quota"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.INFRA
    assert a.one_click is False


def test_oom_without_maxcore_line_is_advisory(tmp_path):
    _folder(tmp_path, "PAL=40\ncharge=0\n")
    a = jf.assess("5", str(tmp_path), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.SLURM_RESOURCE
    assert a.one_click is False and a.fix is None


def test_second_identical_oom_escalates(tmp_path):
    ap = tmp_path / "att.json"
    jf.register_attempt("12345", "out-of-memory", ap)
    a = jf.assess("12345", str(_folder(tmp_path)), ["out-of-memory"],
                  attempts_path=ap)
    assert a.escalate is True and a.one_click is False


# --- apply via DELFIN's recalc path -----------------------------------------

def test_apply_memory_edits_control_keeps_backup_and_recalcs(tmp_path, monkeypatch):
    monkeypatch.setattr(jf, "parse_pal_maxcore", lambda txt: (40, 13500))
    f = _folder(tmp_path)
    a = jf.assess("12345", str(f), ["out-of-memory"], attempts_path=tmp_path / "a.json")
    seen = {}
    def fake_submit(**kw):
        seen.update(kw)
        return SimpleNamespace(returncode=0, stdout="Submitted batch job 99999")
    res = jf.apply_via_recalc(a.fix, "12345", str(f), submit_delfin_fn=fake_submit,
                              attempts_path=tmp_path / "a.json")
    assert res["ok"] is True and res["new_job_id"] == "99999"
    # CONTROL.txt was bumped in place ...
    assert "maxcore=13500" in (f / "CONTROL.txt").read_text()
    # ... a backup preserving the original was kept (nothing destroyed) ...
    backups = list(f.glob("CONTROL.txt.bak-*"))
    assert backups and "maxcore=9000" in backups[0].read_text()
    # ... and it resubmitted through DELFIN's recalc path, not a raw sbatch.
    assert seen["mode"] == "delfin-recalc-classic"
    assert seen["pal"] == 40 and seen["maxcore"] == 13500
    assert jf.attempts_for("12345", "out-of-memory", tmp_path / "a.json") == 1


def test_apply_walltime_passes_time_limit_no_control_edit(tmp_path, monkeypatch):
    monkeypatch.setattr(jf, "parse_pal_maxcore", lambda txt: (40, 9000))
    f = _folder(tmp_path)
    fix = jf.Fix(kind="walltime", summary_line="x", new_time_limit="48:00:00")
    seen = {}
    res = jf.apply_via_recalc(fix, "7", str(f),
                              submit_delfin_fn=lambda **kw: seen.update(kw) or
                              SimpleNamespace(returncode=0, stdout="batch job 5"))
    assert res["ok"] is True
    assert seen["time_limit"] == "48:00:00"
    assert "maxcore=9000" in (f / "CONTROL.txt").read_text()   # untouched
    assert not list(f.glob("CONTROL.txt.bak-*"))               # no edit, no backup


def test_apply_refuses_when_control_changed(tmp_path, monkeypatch):
    monkeypatch.setattr(jf, "parse_pal_maxcore", lambda txt: (40, 9000))
    f = _folder(tmp_path)
    a = jf.assess("12345", str(f), ["out-of-memory"], attempts_path=tmp_path / "a.json")
    (f / "CONTROL.txt").write_text("PAL=40\nmaxcore=99999\n")  # changed since prep
    res = jf.apply_via_recalc(a.fix, "12345", str(f),
                              submit_delfin_fn=lambda **kw: SimpleNamespace(returncode=0, stdout=""),
                              attempts_path=tmp_path / "a.json")
    assert res["ok"] is False and "changed" in res["error"]


def test_apply_reports_recalc_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(jf, "parse_pal_maxcore", lambda txt: (40, 9000))
    f = _folder(tmp_path)
    fix = jf.Fix(kind="retry", summary_line="resubmit")
    res = jf.apply_via_recalc(fix, "9", str(f),
                              submit_delfin_fn=lambda **kw: SimpleNamespace(
                                  returncode=1, stdout="", stderr="queue full"))
    assert res["ok"] is False and "queue full" in res["error"]
