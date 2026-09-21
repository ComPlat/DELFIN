"""Tests for the ``delfin doctor`` CLI subcommand (installation self-check)."""

from __future__ import annotations

import json

import pytest

from delfin import cli as delfin_cli
from delfin.doctor import CheckResult


def _results(*rows):
    return [CheckResult(*row) for row in rows]


def test_doctor_prints_one_line_per_check_and_fix_hints(capsys, monkeypatch):
    monkeypatch.setattr(
        "delfin.doctor.run_all",
        lambda scratch_dir=None: _results(
            ("orca", "ok", "orca 6.0.0 at /opt/orca", ""),
            ("slurm", "missing", "sinfo not found", "install SLURM"),
            ("scratch_dir", "broken", "nope: EACCES", "chmod it"),
        ),
    )
    rc = delfin_cli._run_doctor_subcommand([])
    out = capsys.readouterr().out
    assert "OK       orca" in out
    assert "MISSING  slurm" in out
    assert "BROKEN   scratch_dir" in out
    # fix hints only under missing/broken rows
    assert "fix: install SLURM" in out
    assert "fix: chmod it" in out
    # summary line
    assert "1 ok, 1 missing, 1 broken (3 checks)" in out
    assert rc == 1  # one broken -> exit code 1


def test_doctor_exit_code_zero_when_only_missing(capsys, monkeypatch):
    monkeypatch.setattr(
        "delfin.doctor.run_all",
        lambda scratch_dir=None: _results(
            ("orca", "ok", "found", ""),
            ("slurm", "missing", "not found", "optional"),
        ),
    )
    rc = delfin_cli._run_doctor_subcommand([])
    assert rc == 0
    out = capsys.readouterr().out
    assert "2 ok" not in out
    assert "1 ok, 1 missing, 0 broken (2 checks)" in out


def test_doctor_json_outputs_raw_records(capsys, monkeypatch):
    monkeypatch.setattr(
        "delfin.doctor.run_all",
        lambda scratch_dir=None: _results(
            ("orca", "ok", "found", ""),
            ("xtb", "broken", "crashed", "reinstall"),
        ),
    )
    rc = delfin_cli._run_doctor_subcommand(["--json"])
    out = capsys.readouterr().out
    data = json.loads(out)
    assert data == [
        {"name": "orca", "status": "ok", "detail": "found", "fix_hint": ""},
        {"name": "xtb", "status": "broken", "detail": "crashed",
         "fix_hint": "reinstall"},
    ]
    assert rc == 1


def test_doctor_passes_scratch_dir_through(monkeypatch):
    seen = {}

    def fake_run_all(scratch_dir=None):
        seen["scratch_dir"] = scratch_dir
        return _results(("scratch_dir", "ok", "fine", ""))

    monkeypatch.setattr("delfin.doctor.run_all", fake_run_all)
    rc = delfin_cli._run_doctor_subcommand(["--scratch", "/tmp/x"])
    assert seen["scratch_dir"] == "/tmp/x"
    assert rc == 0


def test_main_dispatches_doctor(monkeypatch, capsys):
    monkeypatch.setattr(
        "delfin.doctor.run_all",
        lambda scratch_dir=None: _results(("orca", "ok", "found", "")),
    )
    rc = delfin_cli.main(["doctor"])
    assert rc == 0
    assert "OK       orca" in capsys.readouterr().out
