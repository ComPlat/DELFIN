"""The GOAT -> g-xTB -> ORCA single-point prescreen stages and submits.

``delfin/analysis_tools/pre_gxtb_sp_workflow/prescreen.py`` is a standalone
script (standard library only, copied next to the run). It is driven here
as it is run on a cluster, as a process, with fake ``xtb``, ``orca`` and
``sbatch`` programs on PATH:

* ``prepare`` stages one folder per conformer, keeps every frame and its
  GOAT index, and leaves the source folder as it was;
* it refuses an existing run folder and an ensemble whose atom order
  changes between frames;
* ``submit`` queues the g-xTB array, the SP array ``afterok`` on it and
  the collector ``afterany`` on that, with the account when one is given,
  and refuses a second submission.

Nothing here judges the chemistry or a cluster's MPI behaviour.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = (Path(__file__).resolve().parents[1] / "delfin" / "analysis_tools"
          / "pre_gxtb_sp_workflow" / "prescreen.py")

ENSEMBLE = """3
-40.123456789 frame zero
O 0.000000 0.000000 0.000000
H 0.957000 0.000000 0.000000
H -0.240000 0.927000 0.000000
3
-40.120000000 frame one
O 0.000000 0.000000 0.100000
H 0.950000 0.000000 0.100000
H -0.230000 0.930000 0.100000
"""


def _module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("prescreen", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _parse_checks(tmp_path: Path) -> None:
    """Where the Slurm/bash workflow cannot run (Windows): the part that
    is plain Python -- reading the ensemble -- is checked directly."""
    mod = _module()
    src = tmp_path / "e.xyz"
    src.write_text(ENSEMBLE)
    frames = mod.read_ensemble(src)
    assert len(frames) == 2 and frames[0][4] == pytest.approx(-40.123456789)
    src.write_text(ENSEMBLE.replace("O 0.000000 0.000000 0.100000",
                                    "N 0.000000 0.000000 0.100000", 1))
    with pytest.raises(ValueError, match="atom order differs"):
        mod.read_ensemble(src)


def _program(folder: Path, name: str, body: str) -> Path:
    """A fake program: a Python script, so it runs wherever Python does."""
    exe = folder / name
    exe.write_text(f"#!{sys.executable}\n{body}")
    exe.chmod(0o755)
    return exe


@pytest.fixture
def cluster(tmp_path, monkeypatch):
    if os.name == "nt":
        # The workflow is a Slurm/bash one; Windows has neither. What a
        # Windows host can check is the staging without programs.
        return None
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    xtb = _program(bin_dir, "xtb", "import sys\nprint('usage: xtb --gxtb ...')\n")
    orca = _program(bin_dir, "orca", "print('orca')\n")
    log = tmp_path / "sbatch.log"
    _program(bin_dir, "sbatch", (
        "import sys, pathlib\n"
        f"log = pathlib.Path({str(log)!r})\n"
        "n = len(log.read_text().splitlines()) if log.exists() else 0\n"
        "with log.open('a') as f:\n"
        "    f.write(' '.join(sys.argv[1:]) + '\\n')\n"
        "print(f'{1000 + n};cluster')\n"))
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.delenv("SBATCH_ACCOUNT", raising=False)
    return {"xtb": xtb, "orca": orca, "log": log}


def _run(*args, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], cwd=cwd,
                          capture_output=True, text=True, timeout=120)


def _ensemble(tmp_path: Path, text: str = ENSEMBLE) -> Path:
    src = tmp_path / "goat" / "molecule.finalensemble.xyz"
    src.parent.mkdir()
    src.write_text(text)
    return src


def test_prepare_stages_every_frame_and_leaves_the_source_alone(tmp_path, cluster):
    if cluster is None:
        _parse_checks(tmp_path)
        return
    src = _ensemble(tmp_path)
    before = sorted(p.name for p in src.parent.iterdir())
    run = tmp_path / "run"
    out = _run("prepare", "--ensemble", str(src), "--run-dir", str(run),
               "--gxtb", str(cluster["xtb"]), "--orca", str(cluster["orca"]),
               "--charge", "0", "--multiplicity", "1", cwd=tmp_path)
    assert out.returncode == 0, out.stdout + out.stderr
    assert sorted(p.name for p in src.parent.iterdir()) == before
    settings = json.loads((run / "settings.json").read_text())
    assert settings["frames"] == 2 and settings["atoms"] == 3
    assert settings["sp_optimize"] is False and settings["account"] is None
    rows = list(csv.DictReader((run / "manifest.csv").open()))
    assert [r["goat_index"] for r in rows] == ["0", "1"]
    assert rows[0]["goat_Eh"].startswith("-40.123456789")
    for i in (0, 1):
        assert (run / "gxtb" / f"conf_{i:04d}" / "start.xyz").is_file()
        inp = (run / "sp" / f"conf_{i:04d}" / "job.inp").read_text()
        assert "* xyzfile 0 1 start.xyz" in inp and " Opt" not in inp
    assert (run / "original.finalensemble.xyz").read_text() == ENSEMBLE
    # A second prepare into the same folder is refused, not merged.
    again = _run("prepare", "--ensemble", str(src), "--run-dir", str(run),
                 "--gxtb", str(cluster["xtb"]), "--orca", str(cluster["orca"]),
                 cwd=tmp_path)
    assert again.returncode != 0 and "Refusing to overwrite" in again.stderr


def test_a_changed_atom_order_is_refused(tmp_path, cluster):
    if cluster is None:
        _parse_checks(tmp_path)
        return
    swapped = ENSEMBLE.replace("O 0.000000 0.000000 0.100000",
                               "H 0.000000 0.000000 0.100000", 1)
    swapped = swapped.replace("H 0.950000 0.000000 0.100000",
                              "O 0.950000 0.000000 0.100000", 1)
    src = _ensemble(tmp_path, swapped)
    out = _run("prepare", "--ensemble", str(src), "--run-dir",
               str(tmp_path / "run"), "--gxtb", str(cluster["xtb"]),
               "--orca", str(cluster["orca"]), cwd=tmp_path)
    assert out.returncode != 0 and "atom order differs" in out.stderr
    assert not (tmp_path / "run").exists()


@pytest.mark.parametrize("account", ["", "my-group"])
def test_submit_chains_three_stages_once(tmp_path, cluster, account):
    if cluster is None:
        _parse_checks(tmp_path)
        return
    src = _ensemble(tmp_path)
    run = tmp_path / "run"
    args = ["prepare", "--ensemble", str(src), "--run-dir", str(run),
            "--gxtb", str(cluster["xtb"]), "--orca", str(cluster["orca"]),
            "--submit"]
    if account:
        args += ["--account", account]
    out = _run(*args, cwd=tmp_path)
    assert out.returncode == 0, out.stdout + out.stderr
    calls = cluster["log"].read_text().splitlines()
    assert len(calls) == 3
    acc = [f"--account={account}"] if account else []
    assert calls[0].split() == ["--parsable", *acc, "--array=0-1", "run_gxtb.sbatch"]
    assert calls[1].split() == ["--parsable", *acc, "--array=0-1",
                                "--dependency=afterok:1000", "run_sp.sbatch"]
    assert calls[2].split() == ["--parsable", *acc,
                                "--dependency=afterany:1001", "collect.sbatch"]
    ids = json.loads((run / "submission.json").read_text())
    assert ids == {"gxtb_array": "1000", "sp_array": "1001", "collect_job": "1002"}
    # The shipped scripts name no account of their own.
    for name in ("run_gxtb.sbatch", "run_sp.sbatch", "collect.sbatch"):
        assert "--account" not in (run / name).read_text()
    # Submitting the same run again is refused: no duplicate arrays.
    again = _run("submit", "--run-dir", str(run), cwd=tmp_path)
    assert again.returncode != 0 and "Refusing to submit" in again.stderr
    assert len(cluster["log"].read_text().splitlines()) == 3
