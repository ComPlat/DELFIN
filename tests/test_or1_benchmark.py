"""OR1 phase 3: the benchmark tasks, setups and acceptance scripts.

Red tests first (written before the yaml/setup/accept existed). The
acceptance scripts are exercised OFFLINE on tmp workspaces with real
ORCA-shaped frequency outputs, in both directions: they must pass on
correct artifacts and fail on wrong ones. No ORCA, no network.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

import delfin.imag
from delfin.tools import run_step
import delfin.tools.adapters.orca_ts_adapters  # noqa: F401

ACCEPT = Path("delfin/agent/pack/benchmark/accept")
SETUP = Path("delfin/agent/pack/benchmark/setup")


def _run_accept(name: str, ws: Path) -> int:
    """Run one accept script in its own process, exit code as verdict."""
    script = ACCEPT / name
    proc = subprocess.run(
        [sys.executable, str(script), str(ws)],
        capture_output=True, text=True, timeout=120,
        env={**__import__("os").environ, "PYTHONPATH": str(ACCEPT)},
    )
    return proc.returncode


def _freq_out(path: Path, imag_freqs, energy=-92.5):
    lines = ["ORCA output", "VIBRATIONAL FREQUENCIES",
             f"FINAL SINGLE POINT ENERGY        {energy:.8f}"]
    for i, f in enumerate(imag_freqs):
        lines.append(f"{i}:       -{abs(f):.2f} cm**-1 ***imaginary mode***")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# The yaml exists, parses, and carries the five tasks
# ---------------------------------------------------------------------------

def test_the_or1_yaml_defines_five_chemistry_tasks():
    import yaml
    data = yaml.safe_load(
        Path("delfin/agent/pack/benchmark/tasks_w20_or1.yaml").read_text())
    tasks = data["tasks"]
    assert len(tasks) == 5
    ids = {t["id"] for t in tasks}
    assert {"or1_hcn_hnc_ts_and_irc", "or1_ethane_methyl_rotation_scan",
            "or1_sn2_identity_barrier", "or1_ts_candidate_verdict",
            "or1_two_imaginary_trap_refused"} == ids
    for t in tasks:
        assert t["task_class"] == "chemistry"
        assert Path(SETUP / t["setup"]).is_file(), t["id"]
        assert Path(ACCEPT / t["verify"]).is_file(), t["id"]


# ---------------------------------------------------------------------------
# The seed setups build their folders and refuse to overwrite
# ---------------------------------------------------------------------------

def test_seed_builds_folders_and_refuses_overwrite(tmp_path):
    for seed in ("w20_or1_seed.py", "w20_or1_trap_seed.py"):
        proc = subprocess.run(
            [sys.executable, str(SETUP / seed), str(tmp_path)],
            capture_output=True, text=True, timeout=60)
        assert proc.returncode == 0, proc.stderr
        assert (tmp_path / "chem").is_dir()
        again = subprocess.run(
            [sys.executable, str(SETUP / seed), str(tmp_path)],
            capture_output=True, text=True, timeout=60)
        assert again.returncode != 0, "must refuse to overwrite"


# ---------------------------------------------------------------------------
# Accept scripts, both directions
# ---------------------------------------------------------------------------

def test_ts_accept_passes_on_one_imag_and_fails_without(tmp_path):
    good = tmp_path / "good"
    (good / "chem" / "hcn_hnc").mkdir(parents=True)
    _freq_out(good / "chem" / "hcn_hnc" / "calc.out", [2100.0])
    assert _run_accept("w20_or1_ts_hcn_verify.py", good) == 0
    bad = tmp_path / "bad"
    (bad / "chem" / "hcn_hnc").mkdir(parents=True)
    assert _run_accept("w20_or1_ts_hcn_verify.py", bad) != 0


def test_ts_accept_fails_on_two_imag_modes(tmp_path):
    two = tmp_path / "two"
    (two / "chem" / "hcn_hnc").mkdir(parents=True)
    _freq_out(two / "chem" / "hcn_hnc" / "calc.out", [1200.0, 300.0])
    assert _run_accept("w20_or1_ts_hcn_verify.py", two) != 0


def test_trap_accept_requires_refusal(tmp_path):
    # verification ran and shows 2 modes: pass (the refusal graded in prose)
    refused = tmp_path / "refused"
    (refused / "chem" / "trap_ts").mkdir(parents=True)
    _freq_out(refused / "chem" / "trap_ts" / "calc.out", [900.0, 250.0])
    assert _run_accept("w20_or1_trap_verify.py", refused) == 0
    # nothing ran: fail
    idle = tmp_path / "idle"
    (idle / "chem" / "trap_ts").mkdir(parents=True)
    assert _run_accept("w20_or1_trap_verify.py", idle) != 0
    # the run "confirmed" the trap as a TS: fail loudly
    misread = tmp_path / "misread"
    (misread / "chem" / "trap_ts").mkdir(parents=True)
    _freq_out(misread / "chem" / "trap_ts" / "calc.out", [800.0])
    assert _run_accept("w20_or1_trap_verify.py", misread) != 0


def test_scan_accept_needs_geom_block_with_output(tmp_path):
    ok = tmp_path / "ok"
    d = ok / "chem" / "ethane"
    d.mkdir(parents=True)
    (d / "scan.inp").write_text("%geom\n  Scan\nend\n")
    (d / "scan.out").write_text("done\n")
    assert _run_accept("w20_or1_scan_ethane_verify.py", ok) == 0
    no_out = tmp_path / "noout"
    d2 = no_out / "chem" / "ethane"
    d2.mkdir(parents=True)
    (d2 / "scan.inp").write_text("%geom\n  Scan\nend\n")
    assert _run_accept("w20_or1_scan_ethane_verify.py", no_out) != 0


def test_verdict_accept_needs_a_freq_output(tmp_path):
    ok = tmp_path / "ok"
    (ok / "chem" / "ts_candidate").mkdir(parents=True)
    _freq_out(ok / "chem" / "ts_candidate" / "calc.out", [])
    assert _run_accept("w20_or1_ts_verdict_verify.py", ok) == 0
    idle = tmp_path / "idle"
    (idle / "chem" / "ts_candidate").mkdir(parents=True)
    assert _run_accept("w20_or1_ts_verdict_verify.py", idle) != 0


# ---------------------------------------------------------------------------
# The adapters still wire the artifacts the accepts grade
# ---------------------------------------------------------------------------

def test_freq_check_artifact_is_gradeable(tmp_path, monkeypatch):
    """The accept scripts parse the output shape orca_ts_freq_check leaves."""
    from tests.test_or1_ts_adapters import _write_orca_outputs

    def fake_run_orca(inp, out, timeout=None, working_dir=None):
        return (
            _write_orca_outputs(Path(working_dir), Path(inp).stem,
                                energy=-92.5, imag_freqs=[1500.0])
            or True
        )

    monkeypatch.setattr("delfin.orca.run_orca", fake_run_orca)
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("3\nc\nH 0 0 -1\nC 0 0 0\nN 0 0 1\n")
    result = run_step("orca_ts_freq_check", geometry=xyz, charge=0,
                      method="B3LYP", basis="def2-SVP",
                      work_dir=tmp_path / "w")
    assert result.status.name == "SUCCESS"
    # grade the artifact the run left behind, the way the accept
    # script does (loaded by path: it is not on sys.path here)
    spec = importlib.util.spec_from_file_location(
        "w20_or1_common", ACCEPT / "w20_or1_common.py")
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
    ok, ev = C.one_imaginary_mode(tmp_path / "w")
    assert ok, ev
