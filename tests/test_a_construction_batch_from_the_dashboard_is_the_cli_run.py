"""The Submit tab's construction batch is ``delfin cluster``: same arguments, same run directory.

Synthetic input only (invented IDs, textbook SMILES).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

from conftest import child_env
from delfin.dashboard import construction_batch as cb

LIST = ("cisplatin;[Cl][Pt-2]([Cl])([NH3+])[NH3+]\n"
        "hexaaqua_fe;[OH2+][Fe-3]([OH2+])([OH2+])([OH2+])([OH2+])[OH2+]\n"
        "cisplatin_again;[Cl][Pt-2]([Cl])([NH3+])[NH3+]\n")


def _run_dir_bytes(root: Path) -> dict:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


def _cb_test_panel(tmp_path, backend=None, batch_text=""):
    ctx = SimpleNamespace(calc_dir=tmp_path / "dash", backend=backend)
    batch = SimpleNamespace(value=batch_text)
    acc = cb.create_construction_batch_panel(ctx, batch)
    return acc, acc.construction_batch_widgets


def test_the_dashboard_and_the_cli_prepare_the_same_run_directory(tmp_path):
    inp = tmp_path / "list.txt"
    inp.write_text(LIST)
    acc, w = _cb_test_panel(tmp_path)
    w["run_name"].value = "r1"
    w["list_path"].value = str(inp)
    w["shard_size"].value = 2
    w["timeout"].value = 600
    w["speed"].value = "1.3"
    w["repeat"].value = 1
    w["buttons"]["Prepare"].click()
    dash = tmp_path / "dash" / cb.CB_RUNS_SUBDIR / "r1"
    assert (dash / "manifest.json").exists()

    cli = tmp_path / "cli" / "r1"
    done = subprocess.run([sys.executable, "-m", "delfin", "cluster", "prepare", "--tool", "manta",
                           "--input", str(inp), "--run-dir", str(cli), "--shard-size", "2",
                           "--timeout", "600", "--speed-factor", "1.3", "--repeat", "1"],
                          capture_output=True, text=True, env=child_env(tmp_path))
    assert done.returncode == 0, done.stderr

    a, b = _run_dir_bytes(dash), _run_dir_bytes(cli)
    assert sorted(a) == sorted(b)
    for name in a:
        if name == "manifest.json":
            ma, mb = json.loads(a[name]), json.loads(b[name])
            assert ma.pop("created") and mb.pop("created")      # the clock, nothing else
            assert ma == mb
        else:
            assert a[name] == b[name], name
    # and the panel shows the command that reproduces it
    shown = w["command"].value
    assert "delfin cluster prepare --tool manta" in shown and "--speed-factor 1.3" in shown


def test_the_batch_field_becomes_the_list_when_no_file_is_given(tmp_path):
    acc, w = _cb_test_panel(tmp_path, batch_text="Ni_1;[Cl][Ni-2]([Cl])([NH3+])[NH3+];charge=0\n\n"
                                                 "Pt_1;[Cl][Pt-2]([Cl])([NH3+])[NH3+]\n")
    w["run_name"].value = "from_field"
    w["buttons"]["Prepare"].click()
    root = tmp_path / "dash" / cb.CB_RUNS_SUBDIR
    assert (root / "from_field.input.txt").read_text() == (
        "Ni_1;[Cl][Ni-2]([Cl])([NH3+])[NH3+]\nPt_1;[Cl][Pt-2]([Cl])([NH3+])[NH3+]\n")
    man = json.loads((root / "from_field" / "manifest.json").read_text())
    assert man["n_systems"] == 2


def test_every_field_reaches_the_cli_argument_it_names():
    f = {"tool": "mace", "mode": "extended", "input": "in.txt", "run_dir": "RUN",
         "select": " sel.txt ", "specs": "", "tool_python": "/env/bin/python", "shard_size": 100,
         "timeout": 3600, "speed_factor": "1.0", "repeat": 0}
    assert cb.construction_batch_argv("prepare", f) == [
        "prepare", "--tool", "mace", "--input", "in.txt", "--run-dir", "RUN", "--mode", "extended",
        "--select", "sel.txt", "--tool-python", "/env/bin/python", "--shard-size", "100",
        "--timeout", "3600"]
    assert cb.construction_batch_argv("slurm", {"run_dir": "RUN", "set": "repeat", "throttle": 10,
                                                "time_limit": "24:00:00",
                                                "setup": ["module load x"]}) == [
        "slurm", "RUN", "--submit", "--set", "repeat", "--throttle", "10", "--time", "24:00:00",
        "--setup", "module load x"]
    from delfin.cluster_bench.cli import cbatch_parser

    for action in ("prepare", "slurm", "status", "collect", "repeat-stats"):
        cbatch_parser().parse_args(cb.construction_batch_argv(action, dict(f, set="main")))


def test_a_refusal_of_the_cli_is_shown_not_raised(tmp_path):
    rc, text = cb.construction_batch_call(["prepare", "--tool", "manta", "--input",
                                           str(tmp_path / "missing.txt"), "--run-dir",
                                           str(tmp_path / "run")])
    assert rc != 0 and not (tmp_path / "run").exists()


def test_on_slurm_submit_runs_the_cli_slurm_step_with_the_sites_modules(tmp_path, monkeypatch):
    inp = tmp_path / "list.txt"
    inp.write_text(LIST)
    backend = SimpleNamespace(backend_name="SlurmJobBackend", slurm_profile="site",
                              _PROFILE_ENV={"site": {"DELFIN_MODULES": "devel/python"}})
    monkeypatch.delenv("DELFIN_MODULES", raising=False)
    acc, w = _cb_test_panel(tmp_path, backend=backend)
    w["run_name"].value = "r2"
    w["list_path"].value = str(inp)
    w["repeat"].value = 1
    w["buttons"]["Prepare"].click()
    calls = []
    monkeypatch.setattr(cb, "construction_batch_call", lambda argv: (calls.append(argv) or (0, "")))
    w["buttons"]["Submit"].click()
    run = str(tmp_path / "dash" / cb.CB_RUNS_SUBDIR / "r2")
    assert calls == [["slurm", run, "--submit", "--setup", "module load devel/python"],
                     ["slurm", run, "--submit", "--set", "repeat", "--setup",
                      "module load devel/python"]]


def test_without_slurm_the_shards_are_built_here_one_after_the_other(tmp_path, monkeypatch):
    inp = tmp_path / "list.txt"
    inp.write_text(LIST)
    acc, w = _cb_test_panel(tmp_path, backend=SimpleNamespace(backend_name="LocalJobBackend"))
    w["run_name"].value = "r3"
    w["list_path"].value = str(inp)
    w["shard_size"].value = 1
    w["buttons"]["Prepare"].click()
    ran = []
    monkeypatch.setattr(cb.subprocess, "run",
                        lambda cmd, **kw: ran.append(cmd) or SimpleNamespace(returncode=0))
    threads = []
    real = cb.construction_batch_run_locally
    monkeypatch.setattr(cb, "construction_batch_run_locally",
                        lambda *a: threads.append(real(*a)) or threads[-1])
    w["buttons"]["Submit"].click()
    threads[0].join(30)
    run = str(tmp_path / "dash" / cb.CB_RUNS_SUBDIR / "r3")
    assert [c[3:] for c in ran] == [["run-shard", run, "--shard", "0", "--set", "main"],
                                    ["run-shard", run, "--shard", "1", "--set", "main"]]
