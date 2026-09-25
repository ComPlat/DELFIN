"""Scan orchestration: run every check, keep a baseline, report diffs.

The report directory is the only place watchpost ever writes. The first
run stores a baseline (authorized keys, binary names, git hooks); later
diff runs report only what is new against it.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from . import checks
from .model import Finding

BASELINE_NAME = "baseline.json"


def collect(home: Path, git_repos: list[Path] | None = None) -> list[Finding]:
    """Every check, read-only, against a home given by the caller."""
    findings: list[Finding] = []
    findings += checks.check_ssh(home)
    findings += checks.check_persistence(home)
    findings += checks.check_crontab()
    findings += checks.parse_last_output(
        checks.run_read_only_command(["last", "-F"]))
    findings += checks.parse_ss_output(
        checks.run_read_only_command(["ss", "-tpn"]) or "",
        baseline_hosts=set())
    findings += checks.check_user_binaries(
        home, baseline_names=set())
    findings += checks.check_credentials(home)
    findings += checks.check_git(git_repos or [], baseline_hooks=set())
    findings += checks.check_delfin_audit(home)
    return findings


def _baseline_path(report: Path) -> Path:
    return report / BASELINE_NAME


def load_baseline(report: Path) -> dict | None:
    p = _baseline_path(report)
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_report(report: Path, findings: list[Finding], mode: str) -> Path:
    """Write the report file. The only write in the whole package."""
    report.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    payload = {
        "mode": mode,
        "generated": stamp,
        "findings": [f.as_dict() for f in findings],
    }
    p = report / f"watchpost-{stamp}.json"
    with p.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return p


def run_scan(home: Path, report: Path, mode: str = "report",
             git_repos: list[Path] | None = None) -> list[Finding]:
    """One scan. mode: report (all findings, refresh baseline),
    diff (only new against the stored baseline), baseline (reset it).
    """
    baseline = load_baseline(report) or {}
    if mode == "baseline":
        baseline = {}
    findings = collect(home, git_repos)
    if mode == "diff":
        prior = [Finding(**d) for d in baseline.get("findings", [])]
        findings = diff_findings(findings, prior)
        _write_report(report, findings, mode)
    else:
        _write_report(report, findings, mode)
        _store_baseline(report, home, findings, git_repos or [])
    return findings


def diff_findings(current: list[Finding],
                  seen: list[Finding]) -> list[Finding]:
    seen_keys = {f.key for f in seen}
    return [f for f in current if f.key not in seen_keys]


def _store_baseline(report: Path, home: Path, findings: list[Finding],
                    git_repos: list[Path]) -> None:
    payload = {
        "authorized_keys": _authorized_keys(home),
        "binaries": sorted(checks.user_binary_names(home)),
        "git_hooks": sorted(
            name for repo in git_repos
            for name in checks.git_hook_names(repo)),
        "findings": [f.as_dict() for f in findings],
    }
    report.mkdir(parents=True, exist_ok=True)
    with _baseline_path(report).open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


def _authorized_keys(home: Path) -> list[str]:
    text = checks._read_text(home / ".ssh" / "authorized_keys")
    if text is None:
        return []
    return [line for line in text.splitlines() if line.strip()]
