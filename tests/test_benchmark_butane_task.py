"""The butane rotation benchmark task (wave 28, bt1) exists and grades
the artifact tree it asks for.

The task: a state-of-the-art n-butane rotational profile end to end —
relaxed dihedral scan (OPI, one coordinate), anti/gauche minima with 0
imaginary modes, two eclipsed TS each with EXACTLY one imaginary mode,
the barrier and the anti-gauche gap in cited literature bands, a
multi-frame-xyz rotation movie, and the level of theory stated. The
verify script re-derives the verdict from the artifact, in the chem_*
pattern.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import yaml

import pytest

REPO = Path(__file__).resolve().parents[1]
PACK = REPO / "delfin" / "agent" / "pack" / "benchmark"
YAML = PACK / "tasks_w28_bt1.yaml"
VERIFY = PACK / "accept" / "w28_bt1_verify.py"
SEED = PACK / "setup" / "w28_bt1_seed.py"

SYN = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"

# --- hand-built fixtures ---------------------------------------------------

BUTANE = "10\nbutane\n" + "".join(
    f"{el} {x:.4f} {y:.4f} {z:.4f}\n" for el, x, y, z in [
        ("C", 0.0, 0.0, 0.0), ("C", 0.0, 0.0, 1.526),
        ("C", 1.415, 0.0, 2.034), ("C", 2.830, 0.0, 1.012),
    ] + [("H", 0.8, 0.8, -0.4)] * 6)

SCAN_INP = ("! r2SCAN-3c TightOpt TightSCF\n%geom\n"
            "Scan\n  { D 0 1 2 3 180.0 180.0 24 }\nend\nend\n"
            "* xyz 0 1\nC 0 0 0\n...\n*\n")


def _freq_out(name: str, imag: int, *, barriers: str = "") -> str:
    head = SYN + f"! {name}\nVIBRATIONAL FREQUENCIES\n"
    modes = ""
    for i in range(imag):
        # ORCA prints the marker and value on ONE line; the parser
        # (delfin.imag.collect_imaginary_modes) matches exactly this.
        modes += f"{6 + i}: -21.86 cm**-1 ***imaginary mode***\n"
    if not imag:
        modes = "no imaginary modes (this is a minimum)\n"
    modes += "Mode 7 84.72 cm**-1\n"
    return (head + modes + barriers
            + "OPTIMIZATION RUN DONE\n")


def _scan_out() -> str:
    return (SYN + "! r2SCAN-3c relaxed scan\n"
            "FINAL SINGLE POINT ENERGY  -157.41000\n")


def _movie(frames: int = 24) -> str:
    return "".join(BUTANE for _ in range(frames))


def _summary(barrier: str = "barrier is 5.10 kcal/mol",
             gap: str = "gauche is 0.89 kcal/mol above anti") -> str:
    return f"Level of theory: r2SCAN-3c. {barrier}. {gap}.\n"


def _good_tree(tmp: Path) -> Path:
    folder = tmp / "chem" / "butane"
    folder.mkdir(parents=True)
    (folder / "start.xyz").write_text(BUTANE)
    (folder / "scan.inp").write_text(SCAN_INP)
    (folder / "scan.out").write_text(_scan_out())
    (folder / "movie.xyz").write_text(_movie())
    (folder / "anti_opt.out").write_text(_freq_out("anti", 0))
    (folder / "gauche_opt.out").write_text(_freq_out("gauche", 0))
    (folder / "ts1_opt.out").write_text(_freq_out("ts1", 1))
    (folder / "ts2_opt.out").write_text(_freq_out("ts2", 1))
    (folder / "summary.md").write_text(_summary())
    return folder


def _run_verify(ws: Path) -> int:
    return subprocess.run(
        [sys.executable, str(VERIFY), str(ws)],
        capture_output=True, text=True).returncode


def _load_verify():
    spec = importlib.util.spec_from_file_location("w28_bt1_verify", VERIFY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- task definition -------------------------------------------------------


def test_task_yaml_loads_with_rubric_and_citations():
    doc = yaml.safe_load(YAML.read_text())
    tasks = [t for t in doc["tasks"]
             if t["id"] == "w28_bt1_butane_rotation_profile"]
    assert tasks, "the butane task is defined"
    task = tasks[0]
    assert task["setup"] == "w28_bt1_seed.py"
    assert task["verify"] == "w28_bt1_verify.py"
    # k >= 3 repeats: the runner aggregates replicates; the yaml carries
    # the per-task fields the runner needs.
    assert task["max_tool_calls"] and task["max_duration_s"]
    # Literature bands, cited in the file header.
    header = YAML.read_text()
    assert "0.6-1.15" in header and "4.5-5.8" in header
    assert "0.89" in header  # gauche-anti gap reference value


def test_setup_seeds_butane(tmp_path):
    out = subprocess.run([sys.executable, str(SEED), str(tmp_path)],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    folder = tmp_path / "chem" / "butane"
    assert (folder / "start.xyz").is_file()
    xyz = (folder / "start.xyz").read_text()
    assert xyz.split("\n")[0].strip() == "10"
    assert not (folder / "start.xyz").read_text().startswith("# ORCA")


def test_verify_passes_on_a_good_artifact_tree(tmp_path):
    _good_tree(tmp_path)
    assert _run_verify(tmp_path) == 0


# --- each rubric line individually fails on a broken tree ------------------


@pytest.mark.parametrize("breaker", [
    "no_scan_input", "no_scan_output", "no_movie", "one_minimum",
    "one_ts", "two_imag_ts", "barrier_out_of_band", "gap_out_of_band",
    "no_level",
])
def test_verify_fails_per_rubric_line(tmp_path, breaker):
    _good_tree(tmp_path)
    folder = tmp_path / "chem" / "butane"
    if breaker == "no_scan_input":
        (folder / "scan.inp").unlink()
    elif breaker == "no_scan_output":
        (folder / "scan.out").unlink()
    elif breaker == "no_movie":
        (folder / "movie.xyz").unlink()
    elif breaker == "one_minimum":
        (folder / "gauche_opt.out").unlink()
    elif breaker == "one_ts":
        (folder / "ts2_opt.out").unlink()
    elif breaker == "two_imag_ts":
        (folder / "ts2_opt.out").write_text(_freq_out("ts2", 2))
    elif breaker == "barrier_out_of_band":
        (folder / "summary.md").write_text(
            _summary(barrier="barrier is 12.0 kcal/mol"))
    elif breaker == "gap_out_of_band":
        (folder / "summary.md").write_text(
            _summary(gap="gauche is 2.50 kcal/mol above anti"))
    elif breaker == "no_level":
        (folder / "summary.md").write_text(
            "barrier is 5.10 kcal/mol. gauche is 0.89 kcal/mol.\n")
    assert _run_verify(tmp_path) != 0, breaker


def test_two_imaginary_ts_is_refused_not_counted_as_ts(tmp_path):
    _good_tree(tmp_path)
    folder = tmp_path / "chem" / "butane"
    (folder / "ts2_opt.out").write_text(_freq_out("ts2", 2))
    out = subprocess.run([sys.executable, str(VERIFY), str(tmp_path)],
                         capture_output=True, text=True)
    assert "no structure with 2+ imaginary modes" in out.stdout
    assert out.returncode != 0
