"""Jobs submitted from a terminal get what dashboard jobs get.

``delfin-step --slurm``, ``delfin-pipeline --slurm`` and the tools runtime
wrote their own batch scripts with no partition -- refused on a cluster that
has no default partition, bwUniCluster among them -- and the runtime's script
had no time limit and no memory either, so the partition defaults applied:
ten minutes and 2000 MB per CPU.

No real queue is touched: sbatch and its --test-only are faked.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pytest

from delfin import slurm_submit


@pytest.fixture
def fake_sbatch(monkeypatch):
    """sbatch that fits every partition in `fits` and records what it was given."""
    state = {"fits": {"cpu", "cpu_il"}, "calls": []}

    class Done:
        def __init__(self, returncode=0, stdout="", stderr=""):
            self.returncode, self.stdout, self.stderr = returncode, stdout, stderr

    def run(cmd, *args, **kwargs):
        cmd = [str(part) for part in cmd]
        state["calls"].append(cmd)
        if "--test-only" in cmd:
            part = next(c.split("=", 1)[1] for c in cmd if c.startswith("--partition="))
            if part in state["fits"]:
                return Done(stderr=f"sbatch: Job 1 to start at 2026-09-16T21:20:09 in partition {part}")
            return Done(1, stderr="sbatch: error: allocation failure: Requested node configuration is not available")
        return Done(stdout="Submitted batch job 4711\n")

    monkeypatch.setattr("subprocess.run", run)
    monkeypatch.setattr("shutil.which", lambda name=None: "/usr/bin/sbatch")
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda *a, **k: {})
    monkeypatch.setenv("DELFIN_SLURM_PARTITIONS", "cpu,cpu_il")
    return state


def _submitted(state):
    return [c for c in state["calls"] if "--test-only" not in c][-1]


def test_the_candidates_come_from_settings_then_environment_then_site(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda *a, **k: {"runtime": {"slurm": {"partitions": "cpu_il"}}})
    monkeypatch.setenv("DELFIN_SLURM_PARTITIONS", "cpu")
    assert slurm_submit.configured_partitions(profile="bwunicluster3") == ("cpu_il",)

    monkeypatch.setattr("delfin.user_settings.load_settings", lambda *a, **k: {})
    assert slurm_submit.configured_partitions(profile="bwunicluster3") == ("cpu",)

    monkeypatch.delenv("DELFIN_SLURM_PARTITIONS")
    assert slurm_submit.configured_partitions(profile="bwunicluster3") == ("cpu", "cpu_il")
    assert slurm_submit.configured_partitions(profile="") == ()


def test_a_script_that_names_its_partition_keeps_it(fake_sbatch, tmp_path):
    script = tmp_path / "job.sh"
    script.write_text("#!/bin/bash\n#SBATCH --partition=dev_cpu\nhostname\n")

    assert slurm_submit.sbatch_command("sbatch", script) == ["sbatch", str(script)]
    assert not [c for c in fake_sbatch["calls"] if "--test-only" in c]


def test_the_module_does_not_import_the_dashboard():
    """Four seconds and ipywidgets, for a command-line submit."""
    import subprocess as real_subprocess

    code = ("import sys, delfin.slurm_submit; "
            "print('delfin.dashboard' in sys.modules, 'ipywidgets' in sys.modules)")
    done = real_subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert done.stdout.strip() == "False False", done.stderr


def test_delfin_step_gets_a_partition_a_time_and_its_memory(fake_sbatch, tmp_path):
    from delfin import cli_step

    fake_sbatch["fits"] = {"cpu"}
    args = argparse.Namespace(step_name="orca_sp", geometry=None, work_dir=str(tmp_path),
                              slurm_time="48h")

    assert cli_step._submit_step_slurm(args, 40, {"maxcore": 6000}) == 0

    script = (tmp_path / "delfin-orca_sp_slurm.sh").read_text()
    assert "#SBATCH --threads-per-core=1" in script, "40 cores, not 20 cores and their threads"
    assert "#SBATCH --time=2-00:00:00" in script
    assert "#SBATCH --mem=240000M" in script
    assert "#SBATCH --nodes=1" in script
    assert "--partition=cpu" in _submitted(fake_sbatch)


def test_delfin_step_refuses_a_time_it_cannot_read(fake_sbatch, tmp_path, capsys):
    from delfin import cli_step

    args = argparse.Namespace(step_name="orca_sp", geometry=None, work_dir=str(tmp_path),
                              slurm_time="two days")
    assert cli_step._submit_step_slurm(args, 4, {}) == 1
    assert "invalid --slurm-time" in capsys.readouterr().err
    assert not fake_sbatch["calls"], "nothing reached SLURM"


def test_delfin_pipeline_keeps_a_yaml_partition_and_chooses_one_otherwise(fake_sbatch, tmp_path):
    pytest.importorskip("yaml")
    from delfin import cli_pipeline

    args = argparse.Namespace(cores="8", geometry=None, work_dir=str(tmp_path), param=[])

    named = tmp_path / "named.yaml"
    named.write_text("name: named\nslurm:\n  partition: dev_cpu\n  time: '00:30:00'\nsteps: []\n")
    assert cli_pipeline._submit_slurm(str(named), args) == 0
    assert not any(part.startswith("--partition=") for part in _submitted(fake_sbatch))

    free = tmp_path / "free.yaml"
    free.write_text("name: free\nslurm:\n  time: '12:00:00'\nsteps: []\n")
    assert cli_pipeline._submit_slurm(str(free), args) == 0
    assert "--partition=cpu,cpu_il" in _submitted(fake_sbatch)
    assert "#SBATCH --threads-per-core=1" in (tmp_path / "free_slurm.sh").read_text()


def test_a_runtime_cluster_run_states_its_time_memory_and_partition(fake_sbatch, tmp_path):
    from delfin.tools._runtime import RunStore, Runtime

    runtime = Runtime(RunStore(tmp_path / "store"))
    handle = runtime.submit_application("any_app", cores=40, maxcore=6000, backend="slurm",
                                        work_dir=tmp_path / "wd", slurm_time="1d12h")

    script = (tmp_path / "wd" / "submit.sh").read_text()
    assert "#SBATCH --threads-per-core=1" in script
    assert "#SBATCH --time=1-12:00:00" in script
    assert "#SBATCH --mem=240000M" in script
    assert "#SBATCH --nodes=1" in script
    assert "--partition=cpu,cpu_il" in _submitted(fake_sbatch)
    assert runtime.store.get(handle.id).metrics.get("slurm_job_id") == "4711"


def test_a_runtime_cluster_run_without_maxcore_still_gets_a_time_limit(fake_sbatch, tmp_path):
    from delfin.tools._runtime import RunStore, Runtime

    Runtime(RunStore(tmp_path / "store")).submit_application(
        "any_app", cores=4, backend="slurm", work_dir=tmp_path / "wd")

    script = (tmp_path / "wd" / "submit.sh").read_text()
    assert "#SBATCH --time=24:00:00" in script, "not the ten-minute partition default"
    assert "--mem=" not in script
