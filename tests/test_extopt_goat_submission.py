"""GOAT with an orca-mlips MLIP as one GPU job (delfin/mlp_tools/extopt_submit.py).

Contract:
- one sbatch job, one node: ORCA's GOAT workers share ``%pal nprocs``
  and the MLIP side is ONE orca-mlips server per job, stopped by the
  EXIT trap at job end;
- ``--cpus-per-task`` covers the workers: max(cpus_per_task, workers),
  workers = goat_max_workers, else %goat nworkers, else 8 ("divided
  among NWorkers (8 by default)", manual 4.10.1);
- ``goat_max_workers`` caps the workers (multiple of 4, Table 4.9);
- ORCA is started by its full path (parallel ORCA runs need it);
- GOAT needs orca-mlips (the delfin-mlp-extopt wrapper reloads the model
  on every ProgExt call) and a GPU partition: no local CPU fallback;
- the plan names the files the job will produce;
- dry run by default, real submission only through delfin.scheduler_client.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.mlp_tools.extopt_submit import (
    ExtOptJobRequest,
    plan_extopt_job,
    submit_extopt_job,
)

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


def _cfg(cpus=2, gpu="gpu_a100"):
    return MLPClusterConfig(
        gpu_partitions=gpu, gres="gpu:1", cpus_per_task=cpus, mem_mb=32000,
        time_limit="12:00:00", env_activation="source /opt/mlip-venv/bin/activate",
    )


def _goat(tmp_path, **kw):
    base = dict(run_type="goat", elements=WATER[0], coordinates=WATER[1],
                provider="orca_mlips", mlips_backend="uma", device="cuda",
                job_name="goat_job", work_dir=tmp_path)
    base.update(kw)
    return ExtOptJobRequest(**base)


def test_goat_script_one_node_workers_and_trap(tmp_path):
    plan = plan_extopt_job(_goat(tmp_path), cfg=_cfg(cpus=2))
    assert plan.mode == "gpu"
    assert plan.pal == 8  # 8 default workers > cpus_per_task 2
    script = plan.script.read_text()
    assert "#SBATCH --nodes=1" in script
    assert "#SBATCH --ntasks=1" in script
    assert "#SBATCH --cpus-per-task=8" in script
    assert "#SBATCH --partition=gpu_a100" in script
    assert "#SBATCH --gres=gpu:1" in script
    assert "trap 'orca-mlips-uma --stop-server 2>/dev/null || true' EXIT" in script
    assert script.count("--stop-server") == 1
    # ORCA by full path, once; no second ORCA, no server started by hand
    assert 'ORCA_BIN="$(command -v orca)"' in script
    assert script.count('"$ORCA_BIN" goat_job.inp > goat_job.out') == 1
    assert "unverified" in script.lower()
    inp = plan.input.read_text().lower()
    assert "goat" in inp and "extopt" in inp
    assert "nprocs 8" in inp


def test_cpus_per_task_wins_when_larger(tmp_path):
    plan = plan_extopt_job(_goat(tmp_path), cfg=_cfg(cpus=16))
    assert plan.pal == 16
    assert "#SBATCH --cpus-per-task=16" in plan.script.read_text()


def test_goat_max_workers_caps_and_is_written(tmp_path):
    plan = plan_extopt_job(
        _goat(tmp_path, goat_max_workers=4, goat_options={"nworkers": 16}),
        cfg=_cfg(cpus=2),
    )
    assert plan.pal == 4
    inp = plan.input.read_text().lower()
    assert "nworkers 4" in inp and "nworkers 16" not in inp


@pytest.mark.parametrize("bad", [0, 6, -4])
def test_goat_max_workers_must_be_a_multiple_of_four(tmp_path, bad):
    with pytest.raises(ValueError, match="multiple of 4"):
        plan_extopt_job(_goat(tmp_path, goat_max_workers=bad), cfg=_cfg())


def test_goat_refused_without_gpu_partition(tmp_path):
    plan = plan_extopt_job(_goat(tmp_path, allow_local=True), cfg=_cfg(gpu=""))
    assert plan.mode == "refuse"
    assert "GPU" in plan.reason
    assert not plan.script.exists()
    with pytest.raises(RuntimeError, match="GPU"):
        submit_extopt_job(_goat(tmp_path), cfg=_cfg(gpu=""), submit=True)


def test_goat_refused_with_the_delfin_wrapper(tmp_path):
    plan = plan_extopt_job(_goat(tmp_path, provider="delfin", backend="ani2x"),
                           cfg=_cfg())
    assert plan.mode == "refuse"
    assert "orca_mlips" in plan.reason


def test_plan_names_the_goat_outputs(tmp_path):
    plan = plan_extopt_job(_goat(tmp_path), cfg=_cfg())
    names = [Path(p).name for p in plan.produces]
    assert "goat_job.finalensemble.xyz" in names
    assert "goat_job.globalminimum.xyz" in names
    assert "goat_job.out" in names


def test_dry_run_does_not_touch_the_scheduler(tmp_path, monkeypatch):
    import delfin.scheduler_client as sc

    calls = []
    monkeypatch.setattr(sc, "query", lambda kind, argv, runner=None: calls.append(argv))
    out = submit_extopt_job(_goat(tmp_path), cfg=_cfg(), submit=False)
    assert out["submitted"] is False and calls == []
    assert any(p.endswith("goat_job.finalensemble.xyz") for p in out["produces"])


def test_submission_goes_through_scheduler_client(tmp_path, monkeypatch):
    import delfin.scheduler_client as sc

    calls = []
    monkeypatch.setattr(
        sc, "query",
        lambda kind, argv, runner=None: calls.append((kind, list(argv))) or "Submitted batch job 7",
    )
    out = submit_extopt_job(_goat(tmp_path), cfg=_cfg(), submit=True)
    assert out["submitted"] is True
    assert calls == [("sbatch", ["sbatch", str(tmp_path / "goat_job.sbatch")])]


def test_opt_script_unchanged_by_goat_lines(tmp_path):
    """An orca-mlips opt job keeps the bare `orca` call and no GOAT note."""
    req = _goat(tmp_path, run_type="opt")
    plan = plan_extopt_job(req, cfg=_cfg(cpus=8))
    script = plan.script.read_text()
    assert "orca goat_job.inp > goat_job.out" in script
    assert "ORCA_BIN" not in script and "GOAT" not in script
