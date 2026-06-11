"""Tests for the agent's failure-assessment + bounded-fix bridge.

Safety-critical contract:
- SLURM resource kills (oom / timelimit) -> a bounded #SBATCH bump the
  agent may prepare; ORCA-internal errors -> defer to DELFIN's in-run
  recovery (never duplicate); infra -> human only.
- A proposed change NEVER touches chemistry/convergence keywords.
- One auto-retry per (job, signature); the second escalates.
- apply keeps a backup and never destroys the original.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import job_fix as jf


_SUBMIT = """#!/bin/bash
#SBATCH --job-name=MarcDy
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=8
module load orca
orca input.inp > MarcDy.out
"""


def _folder(tmp_path: Path, submit: str = _SUBMIT) -> Path:
    (tmp_path / "submit.sh").write_text(submit)
    return tmp_path


# --- resource-bump math -----------------------------------------------------

def test_mem_bump_is_1_5x_same_unit():
    assert jf._bump_mem_token("16G") == "24G"
    assert jf._bump_mem_token("16000M") == "24000M"
    assert jf._bump_mem_token("10gb") == "15gb"


def test_time_bump_doubles_all_formats():
    assert jf._bump_time_token("12:00:00") == "24:00:00"
    assert jf._bump_time_token("01:30:00") == "03:00:00"
    assert jf._bump_time_token("30:00") == "01:00:00"        # MM:SS
    assert jf._bump_time_token("1-00:00:00") == "2-00:00:00"  # D-HH:MM:SS


def test_directive_bump_reads_real_lines(tmp_path):
    f = _folder(tmp_path)
    text = (f / "submit.sh").read_text()
    old, new = jf._bump_directive(text, "memory")
    assert old.strip() == "#SBATCH --mem=16G"
    assert new.strip() == "#SBATCH --mem=24G"


# --- the chemistry guard (the line that must never be crossed) --------------

def test_chemistry_guard_blocks_convergence_keywords():
    assert jf._is_chemistry_safe("#SBATCH --mem=16G", "#SBATCH --mem=24G")
    assert not jf._is_chemistry_safe("! TightSCF", "! LooseSCF")
    assert not jf._is_chemistry_safe("TolE 1e-8", "TolE 1e-5")
    assert not jf._is_chemistry_safe("%method functional B3LYP", "x")


# --- routing ----------------------------------------------------------------

def test_oom_kill_prepares_one_click_memory_fix(tmp_path):
    a = jf.assess("12345", str(_folder(tmp_path)), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.SLURM_RESOURCE
    assert a.one_click is True
    assert a.proposal and a.proposal.kind == "memory"
    assert "--mem=24G" in a.proposal.new_line


def test_timelimit_prepares_walltime_fix(tmp_path):
    a = jf.assess("7", str(_folder(tmp_path)), ["slurm-timelimit"],
                  attempts_path=tmp_path / "att.json")
    assert a.one_click is True
    assert a.proposal.kind == "walltime"
    assert "--time=24:00:00" in a.proposal.new_line


def test_orca_internal_defers_to_delfin_no_one_click(tmp_path):
    a = jf.assess("9", str(tmp_path), ["scf-not-converged"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.DELFIN_AUTORECOVERS
    assert a.one_click is False
    assert "in-run" in a.recommendation.lower() or "in run" in a.recommendation.lower()
    assert "never loosen convergence" in a.recommendation.lower()


def test_infra_is_human_only(tmp_path):
    a = jf.assess("3", str(tmp_path), ["disk-quota"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.INFRA
    assert a.one_click is False


def test_missing_submit_script_is_advisory_not_one_click(tmp_path):
    # oom kill but no submit script present -> advise, don't offer a click
    a = jf.assess("5", str(tmp_path), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    assert a.fix_class == jf.FixClass.SLURM_RESOURCE
    assert a.one_click is False
    assert a.proposal is None


# --- retry budget -----------------------------------------------------------

def test_second_identical_failure_escalates(tmp_path):
    ap = tmp_path / "att.json"
    jf.register_attempt("12345", "out-of-memory", ap)   # already retried once
    a = jf.assess("12345", str(_folder(tmp_path)), ["out-of-memory"],
                  attempts_path=ap)
    assert a.escalate is True
    assert a.one_click is False


# --- apply (confirm-gated) --------------------------------------------------

def test_apply_edits_keeps_backup_and_resubmits(tmp_path):
    f = _folder(tmp_path)
    a = jf.assess("12345", str(f), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    calls = {}
    def fake_submit(p: Path) -> str:
        calls["path"] = p
        return "Submitted batch job 99999"
    res = jf.apply_resource_fix(a.proposal, "12345", submit_fn=fake_submit,
                                attempts_path=tmp_path / "att.json")
    assert res["ok"] is True
    assert res["new_job_id"] == "99999"
    # original directive replaced
    assert "--mem=24G" in (f / "submit.sh").read_text()
    # a backup copy preserving the original was kept (nothing destroyed)
    backups = list(f.glob("submit.sh.bak-*"))
    assert backups and "--mem=16G" in backups[0].read_text()
    # the retry was registered so it can't loop forever
    assert jf.attempts_for("12345", "out-of-memory", tmp_path / "att.json") == 1


def test_apply_refuses_when_script_changed(tmp_path):
    f = _folder(tmp_path)
    a = jf.assess("12345", str(f), ["out-of-memory"],
                  attempts_path=tmp_path / "att.json")
    (f / "submit.sh").write_text("#!/bin/bash\n# totally different now\n")
    res = jf.apply_resource_fix(a.proposal, "12345",
                                submit_fn=lambda p: "x",
                                attempts_path=tmp_path / "att.json")
    assert res["ok"] is False
    assert "changed" in res["error"]
