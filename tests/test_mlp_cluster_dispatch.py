"""Run MLP jobs where they belong: GPU on the cluster, CPU here.

cluster_dispatch (delfin/mlp_tools/cluster_dispatch.py) renders the sbatch
script a GPU run needs, decides cluster vs. local, and caps concurrent
sbatch MLP jobs at 4.  A batch of N structures is ONE job, not N.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.mlp_tools import cluster_dispatch as cd
from delfin.mlp_tools.cluster_dispatch import (
    MAX_CONCURRENT_MLP_JOBS,
    MLPClusterConfig,
    assert_can_submit,
    decide_dispatch,
    plan_batch,
    render_mlp_sbatch,
)

WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"


def _xyz(tmp_path, name="mol.xyz"):
    p = tmp_path / name
    p.write_text(WATER)
    return p


@pytest.fixture
def cluster_env(monkeypatch):
    """A cluster as settings + env describe it: sbatch on PATH, GPU partition named."""
    monkeypatch.setattr(cd, "slurm_available", lambda: True)
    return MLPClusterConfig.from_sources(
        settings={"slurm": {"gpu_partitions": "gpu_a100", "partition": "cpu"}},
        env={},
    )


# ---------------------------------------------------------------- sbatch text

def test_sbatch_text_requests_gpu_and_time(tmp_path, cluster_env):
    script = render_mlp_sbatch(
        cluster_env, job_dir=tmp_path, xyz=_xyz(tmp_path),
        backend="mace_off", task="optimize", job_name="mlp_opt",
    )
    assert "#!/bin/bash" in script
    assert "--gres=gpu:1" in script
    assert "--partition=gpu_a100" in script
    assert "--time=" in script
    assert "--cpus-per-task=" in script
    assert "--mem=" in script
    assert "--job-name=mlp_opt" in script


def test_sbatch_text_activates_env_from_config_not_hardcoded_paths(tmp_path, cluster_env):
    cfg = MLPClusterConfig.from_sources(
        settings={
            "slurm": {"gpu_partitions": "gpu_a100"},
            "mlp": {"env_activation": "source /opt/mlp/env.sh"},
        },
        env={},
    )
    script = render_mlp_sbatch(
        cfg, job_dir=tmp_path, xyz=_xyz(tmp_path),
        backend="ani2x", task="single_point", job_name="mlp_sp",
    )
    assert "source /opt/mlp/env.sh" in script
    # no home paths, no user names baked into the module
    assert "/home/" not in script


def test_sbatch_text_caps_threads_even_on_gpu(tmp_path, cluster_env):
    script = render_mlp_sbatch(
        cluster_env, job_dir=tmp_path, xyz=_xyz(tmp_path),
        backend="ani2x", task="optimize", job_name="mlp_opt", cpus_per_task=2,
    )
    assert "OMP_NUM_THREADS=2" in script
    assert "MKL_NUM_THREADS=2" in script
    assert "--cpus-per-task=2" in script


def test_sbatch_text_runs_the_batch_runner_for_rank(tmp_path, cluster_env):
    script = render_mlp_sbatch(
        cluster_env, job_dir=tmp_path, xyz=_xyz(tmp_path), batch=_xyz(tmp_path, "ens.xyz"),
        backend="ani2x", task="rank", job_name="mlp_rank",
    )
    assert "delfin.dashboard.local_runner" in script
    assert "DELFIN_MLP_TASK=rank" in script
    assert "DELFIN_MLP_BATCH_FILE" in script


def test_dry_run_writes_script_and_does_not_submit(tmp_path, cluster_env, monkeypatch):
    """submit=False (the default) never runs sbatch."""
    called = []
    monkeypatch.setattr(
        cd.subprocess, "run",
        lambda *a, **k: called.append(a)
        or type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})())
    out = cd.write_and_submit(
        cluster_env, job_dir=tmp_path, xyz=_xyz(tmp_path),
        backend="ani2x", task="optimize", job_name="mlp_opt",
    )
    assert called == [], "dry run must not execute sbatch"
    assert out["submitted"] is False
    assert Path(out["script"]).is_file()


# ------------------------------------------------------- cluster/local choice

def test_cluster_when_slurm_and_gpu_configured(cluster_env):
    target, reason = decide_dispatch(cluster_env)
    assert target == "cluster"
    assert "gpu" in reason.lower()


def test_fallback_local_when_no_cluster(monkeypatch):
    monkeypatch.setattr(cd, "slurm_available", lambda: False)
    cfg = MLPClusterConfig.from_sources(settings={"slurm": {"gpu_partitions": "gpu"}}, env={})
    target, reason = decide_dispatch(cfg)
    assert target == "local"
    assert "slurm" in reason.lower() or "cluster" in reason.lower()


def test_fallback_local_when_cluster_has_no_gpu_partition(monkeypatch):
    monkeypatch.setattr(cd, "slurm_available", lambda: True)
    cfg = MLPClusterConfig.from_sources(settings={"slurm": {"gpu_partitions": ""}}, env={})
    target, reason = decide_dispatch(cfg)
    assert target == "local"
    assert "partition" in reason.lower()


def test_fallback_local_respects_explicit_cpu_capping():
    """The local fallback runs with <= 2 threads and a time cap."""
    plan = cd.local_plan(threads=2, minutes=30)
    assert plan["threads"] <= 2
    assert plan["minutes"] <= 30
    assert plan["threads_env"]["OMP_NUM_THREADS"] == "2"


def test_force_local_overrides_a_configured_cluster(cluster_env):
    target, _ = decide_dispatch(cluster_env, force_local=True)
    assert target == "local"


# ------------------------------------------------------------------- job cap

def test_job_cap_is_four():
    assert MAX_CONCURRENT_MLP_JOBS == 4


def test_cap_refuses_the_fifth_job():
    for i in range(MAX_CONCURRENT_MLP_JOBS):
        assert_can_submit(active=i)
    with pytest.raises(cd.JobCapExceeded):
        assert_can_submit(active=MAX_CONCURRENT_MLP_JOBS)


def test_cap_message_names_the_limit():
    with pytest.raises(cd.JobCapExceeded) as exc:
        assert_can_submit(active=4)
    assert "4" in str(exc.value)


# ------------------------------------------------- batch = one job, not N

def test_plan_batch_of_n_is_one_job(tmp_path):
    ens = tmp_path / "ens.xyz"
    ens.write_text("".join(WATER for _ in range(7)))  # 7 frames = 7 structures
    n = plan_batch(ens)
    assert n == 7
    jobs = cd.jobs_for_batch(ens)
    assert len(jobs) == 1, "a batch of N structures goes as ONE job"


# ------------------------------------------------------- result schema (rank)

def test_local_rank_result_schema_matches_adapter(tmp_path, monkeypatch):
    """The result file carries the same keys as a local adapter run:
    ranking [{index, energy_eV}], n_structures, backend, status, task."""
    ens = tmp_path / "ens.xyz"
    ens.write_text("".join(WATER for _ in range(3)))
    run = tmp_path / "run"
    run.mkdir()
    monkeypatch.chdir(run)
    monkeypatch.setenv("DELFIN_MLP_TASK", "rank")
    monkeypatch.setenv("DELFIN_MLP_BATCH_FILE", str(ens))
    monkeypatch.setenv("DELFIN_MLP_BACKEND", "ani2x")
    monkeypatch.setenv("DELFIN_CHARGE", "0")
    monkeypatch.setenv("DELFIN_MULT", "1")
    monkeypatch.setenv("DELFIN_PAL", "2")

    fake_rank = {"ranking": [(0, -1.5), (2, -1.4), (1, -1.2)], "n_structures": 3,
                 "backend": "ani2x"}
    monkeypatch.setattr("delfin.mlp_tools.calculators.rank_structures",
                        lambda structures, **kw: fake_rank["ranking"])
    from delfin.dashboard import local_runner
    assert local_runner._run_mlp_job(xyz_file=str(_xyz(tmp_path)), delfin_pal="2") == 0
    data = json.loads((run / "mlp_result.json").read_text())
    assert data["task"] == "rank"
    assert data["n_structures"] == 3
    assert data["backend"] == "ani2x"
    assert data["ranking"] == [{"index": 0, "energy_eV": -1.5},
                               {"index": 2, "energy_eV": -1.4},
                               {"index": 1, "energy_eV": -1.2}]
    assert data["status"] == "success"


# ------------------------------------------------- old default path unchanged

def test_old_default_path_is_unchanged_without_new_config(tmp_path, monkeypatch):
    """No DELFIN_MLP_CLUSTER settings and no runtime.slurm GPU partitions ->
    the dispatch decision is 'local' and a single-point run behaves exactly
    as before (adapters' execute() called with device='cpu')."""
    monkeypatch.setattr(cd, "slurm_available", lambda: False)
    cfg = MLPClusterConfig.from_sources(settings={}, env={})
    target, _ = decide_dispatch(cfg)
    assert target == "local"

    from delfin.tools.adapters import mlp as mlp_adapter
    calls = {}
    monkeypatch.setattr("delfin.mlp_tools.calculators.single_point",
                        lambda atoms, **kw: calls.update(kw) or {"energy": -1.0,
                                                                "forces": None,
                                                                "backend": kw.get("backend"),
                                                                "n_atoms": len(atoms)})
    mlp_adapter_module = mlp_adapter
    # the adapter's execute path goes through run_step; call single_point path via adapter
    from delfin.tools._types import StepStatus
    result = mlp_adapter_module.MlpSinglePointAdapter().execute(
        tmp_path, geometry=_xyz(tmp_path), backend="ani2x", device="cpu")
    assert result.status == StepStatus.SUCCESS
    assert calls.get("device") == "cpu"
