"""C6 regression: ORCA/Slurm CPU counts, CUDA fallback and sbatch success."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.mlp_tools.extopt_submit import (
    ExtOptJobRequest,
    plan_extopt_job,
    submit_extopt_job,
)


_WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.75, 0.58, 0.0], [-0.75, 0.58, 0.0]])


def _req(tmp_path: Path, **changes) -> ExtOptJobRequest:
    data = dict(run_type="opt", elements=_WATER[0], coordinates=_WATER[1],
                provider="orca_mlips", mlips_backend="aimnet2", device="cuda",
                work_dir=tmp_path, job_name="c6_opt")
    data.update(changes)
    return ExtOptJobRequest(**data)


def _cfg(gpu: bool = True) -> MLPClusterConfig:
    return MLPClusterConfig(gpu_partitions="gpu_a100" if gpu else "",
                            gres="gpu:1", cpus_per_task=8)


def test_non_goat_input_pal_matches_sbatch(tmp_path):
    plan = plan_extopt_job(_req(tmp_path), cfg=_cfg())
    inp = plan.input.read_text().lower()
    script = plan.script.read_text()
    assert plan.mode == "gpu" and plan.pal == 8
    assert "nprocs 8" in inp
    assert "--device cuda" in inp
    assert "#SBATCH --cpus-per-task=8" in script
    assert 'ORCA_BIN="$(command -v orca)"' in script
    assert '"$ORCA_BIN" c6_opt.inp > c6_opt.out' in script


def test_cuda_without_partition_is_refused_not_silently_cpu(tmp_path):
    plan = plan_extopt_job(_req(tmp_path, allow_local=True), cfg=_cfg(False))
    assert plan.mode == "refuse"
    assert "CUDA" in plan.reason
    assert not plan.script.exists()


def test_sbatch_stderr_is_never_marked_submitted(tmp_path, monkeypatch):
    from delfin import scheduler_client

    def fail(kind, argv, runner=None):
        return "sbatch: error: invalid account or partition"

    monkeypatch.setattr(scheduler_client, "query", fail)
    with pytest.raises(RuntimeError, match="sbatch failed"):
        submit_extopt_job(_req(tmp_path), cfg=_cfg(), submit=True)
