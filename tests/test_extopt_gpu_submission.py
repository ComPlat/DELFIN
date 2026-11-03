"""GPU-aware submission of an ORCA+ExtOpt job (delfin/mlp_tools/extopt_submit.py).

Contract (task om1 item 2):
- ORCA process and the MLIP server must run on the SAME node: one sbatch
  job, no partition split, --gres for the GPU.
- GPU partition configured -> sbatch with --partition from
  runtime.slurm.gpu_partitions, --gres from config, CPU count for %pal,
  env (venv/conda of the MLIP stack) activated in the job script; server
  lifetime bound to the job (orca-mlips scopes auto-servers per parent
  ORCA process; a non-server local wrapper gets a trap-based stop).
- No GPU partition -> local CPU run within the host limits (<= 2 cores,
  < 30 min) or a refusal with the reason.
- All submission through delfin.scheduler_client (the direct-call guard).
- Child env via delfin.agent.lifeline.child_env.
- Dry run (default): prints the input + job script without submitting.
"""

from __future__ import annotations

from pathlib import Path

from delfin.mlp_tools.extopt_submit import (
    LOCAL_CPU_LIMITS,
    ExtOptJobRequest,
    render_extopt_job_script,
    plan_extopt_job,
    submit_extopt_job,
)

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


def _gpu_cfg():
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig

    return MLPClusterConfig(
        gpu_partitions="gpu_a100,gpu", gres="gpu:1", cpus_per_task=8,
        mem_mb=32000, time_limit="04:00:00",
        env_activation="source /opt/mlip-venv/bin/activate",
    )


def _cpu_cfg():
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig

    return MLPClusterConfig(gpu_partitions="", partition="cpu", gres="gpu:1")


def _request(work_dir, **kw):
    base = dict(
        run_type="opt", elements=WATER[0], coordinates=WATER[1],
        backend="mace_off", provider="delfin",
        work_dir=work_dir,
    )
    base.update(kw)
    return ExtOptJobRequest(**base)


# ------------------------------------------------------------------ script


def test_gpu_script_carries_partition_gres_and_env(tmp_path):
    script = render_extopt_job_script(
        _request(tmp_path), cfg=_gpu_cfg(), pal=8, provider_cmd="delfin-mlp-extopt"
    )
    assert "#SBATCH --partition=gpu_a100,gpu" in script
    assert "#SBATCH --gres=gpu:1" in script
    assert "#SBATCH --nodes=1" in script
    assert "source /opt/mlip-venv/bin/activate" in script
    assert "--nodes=1" in script  # same-node guarantee for ORCA + MLIP server


def test_gpu_script_stops_server_at_job_end(tmp_path):
    script = render_extopt_job_script(
        _request(tmp_path), cfg=_gpu_cfg(), pal=8, provider_cmd="orca-mlips-mace"
    )
    # orca-mlips stops its auto-server with the parent; our wrapper gets a trap.
    assert "orca-mlips-mace --stop-server" in script
    assert "trap" in script


def test_local_cpu_limits_declared():
    assert LOCAL_CPU_LIMITS.cores == 2
    assert LOCAL_CPU_LIMITS.minutes == 30


def test_no_gpu_partition_cpu_config_plans_local(tmp_path):
    req = _request(tmp_path, allow_local=True)
    plan = plan_extopt_job(req, cfg=_cpu_cfg())
    assert plan.mode == "local"
    assert plan.pal <= LOCAL_CPU_LIMITS.cores


def test_no_gpu_partition_no_slurm_asks_local_or_refuses(tmp_path):
    req = _request(tmp_path, allow_local=True)
    plan = plan_extopt_job(req, cfg=_cpu_cfg())
    assert plan.mode in ("local", "refuse")


def test_large_system_refused_locally(tmp_path):
    elements = ["C"] * 60
    coords = [[float(i), 0.0, 0.0] for i in range(60)]
    req = _request(
        tmp_path, elements=elements, coordinates=coords,
        allow_local=False,
    )
    plan = plan_extopt_job(req, cfg=_cpu_cfg())
    assert plan.mode == "refuse"
    assert "2" in plan.reason or "local" in plan.reason.lower()


# ------------------------------------------------------------------ dry run


def test_dry_run_prints_input_and_script_without_submitting(tmp_path, monkeypatch):
    calls = []
    import delfin.scheduler_client as sc

    monkeypatch.setattr(
        sc, "query", lambda kind, argv, runner=None: calls.append(argv) or "Submitted batch job 1"
    )
    req = _request(tmp_path)
    result = submit_extopt_job(req, cfg=_gpu_cfg(), submit=False)
    assert result["submitted"] is False
    assert calls == []
    text = Path(result["script"]).read_text()
    assert "#SBATCH" in text
    inp = Path(result["input"]).read_text()
    assert "extopt" in inp.lower()  # OPI writes '!extopt', the fallback 'ExtOpt'


def test_real_submission_goes_through_scheduler_client(tmp_path, monkeypatch):
    calls = []
    import delfin.scheduler_client as sc

    monkeypatch.setattr(
        sc, "query", lambda kind, argv, runner=None: calls.append(argv) or "Submitted batch job 1"
    )
    req = _request(tmp_path)
    result = submit_extopt_job(req, cfg=_gpu_cfg(), submit=True)
    assert result["submitted"] is True
    assert len(calls) == 1
    assert calls[0][0] == "sbatch"
    assert str(Path(result["script"])) in calls[0][1] or calls[0][1][-1].endswith(".sbatch")


def test_provider_orca_mlips_builds_mlips_input(tmp_path):
    """provider='orca_mlips' routes to orca_mlips.build_orca_mlips_input."""
    req = _request(tmp_path, provider="orca_mlips", mlips_backend="mace")
    plan = plan_extopt_job(req, cfg=_gpu_cfg())
    inp = Path(plan.input).read_text()
    assert 'progext "orca-mlips-mace"' in inp.lower()


def test_provider_default_is_delfin_wrapper(tmp_path):
    req = _request(tmp_path)
    plan = plan_extopt_job(req, cfg=_gpu_cfg())
    inp = Path(plan.input).read_text()
    assert 'progext "delfin-mlp-extopt"' in inp.lower()
