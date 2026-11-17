"""Model jobs on the cluster (task hf2).

Given a registered model (hf1's registry shape; a Protocol + stub lives in
delfin.mlp_tools.model_jobs because hf1 is not merged yet) and a task (e.g.
transcribe these audio files with faster-whisper; run an inference script),
DELFIN renders an sbatch script that:

* requests the partition/GPU from the SITE'S REAL partitions, read with
  scontrol at run time (delfin.slurm_submit.discover_gpu_partitions), never a
  hard-coded host; the user may pin a partition;
* uses a per-job venv created once (name hashed from the pinned
  requirements), sets HF_HOME / the model cache to the already-fetched model
  and switches the hub to OFFLINE (no download inside the job unless allowed
  -- the third-party mirror stays refused);
* writes outputs to the job folder and a run record (inputs, model+revision,
  versions, timings).

model_job_plan is a DRY RUN (renders the script + resources, submits
nothing); submission goes through the existing gate
(delfin.scheduler_client + the concurrent-job cap in cluster_dispatch).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools import cluster_dispatch as cd
from delfin.mlp_tools import model_jobs as mj
from delfin.mlp_tools.cluster_dispatch import MAX_CONCURRENT_MLP_JOBS


#: Fake ``scontrol show partition --oneliner``: a GPU partition and a CPU one.
FAKE_SCONTROL = (
    "PartitionName=gpu_a100 State=UP TRES=cpu=64,mem=1000000M,gres/gpu=8 "
    "Nodes=n01\n"
    "PartitionName=cpu State=UP TRES=cpu=192,mem=2000000M "
    "Nodes=c01\n"
    "PartitionName=dev_test State=UP TRES=cpu=8,gres/gpu=1 Nodes=d01\n"
)


def _model(*, model_id="Systran/faster-whisper-large-v3-turbo",
           revision="abc123", local_path=None) -> dict:
    return mj.model_record(
        model_id=model_id,
        revision=revision,
        local_path=local_path or Path("/models/faster-whisper-large-v3-turbo"),
        kind="huggingface",
    )


def _req(tmp_path, *, model=None, gpu=True, partition=None, **kw) -> mj.ModelJobRequest:
    model = model if model is not None else _model()
    audio = tmp_path / "in"
    audio.mkdir(exist_ok=True)
    (audio / "a.mp3").write_text("x")
    (audio / "b.mp3").write_text("x")
    defaults = dict(
        model=model,
        kind="transcribe",
        inputs=[audio / "a.mp3", audio / "b.mp3"],
        job_name="whisper_job",
        gpu=gpu,
        partition=partition,
        out_dir=tmp_path / "out",
        requirements="faster-whisper==1.2.1\nctranslate2==4.5.0\n",
    )
    defaults.update(kw)
    return mj.ModelJobRequest(**defaults)


# ---------------------------------------------------------- partition choice


def test_gpu_partition_comes_from_the_site_scontrol(monkeypatch):
    """The GPU partition is the site's real one, read via scontrol."""
    pick = mj.choose_partition(gpu=True, scontrol_output=FAKE_SCONTROL)
    assert pick.partition == "gpu_a100"
    assert "scontrol" not in pick.reason or "gpu" in pick.partition


def test_cpu_partition_when_no_gpu_wanted(monkeypatch):
    pick = mj.choose_partition(gpu=False, scontrol_output=FAKE_SCONTROL)
    assert pick.partition == "cpu"


def test_dev_partitions_are_skipped(monkeypatch):
    pick = mj.choose_partition(gpu=True, scontrol_output=FAKE_SCONTROL)
    assert pick.partition != "dev_test"


def test_pinned_partition_wins_over_discovery(tmp_path):
    req = _req(tmp_path, partition="gpu_h100_pinned")
    pick = mj.choose_partition(gpu=True, pinned=req.partition,
                               scontrol_output=FAKE_SCONTROL)
    assert pick.partition == "gpu_h100_pinned"


def test_no_partition_when_scontrol_unreachable(monkeypatch):
    pick = mj.choose_partition(gpu=True, scontrol_output="")
    assert pick.partition is None  # no hard-coded site fallback


# -------------------------------------------------------------- sbatch text


def test_sbatch_requests_gpu_and_gres(tmp_path):
    req = _req(tmp_path)
    script = mj.render_model_sbatch(req, partition="gpu_a100")
    assert "#SBATCH --partition=gpu_a100" in script
    assert "#SBATCH --gres=gpu:1" in script


def test_cpu_sbatch_has_no_gres(tmp_path):
    req = _req(tmp_path, gpu=False)
    script = mj.render_model_sbatch(req, partition="cpu")
    assert "#SBATCH --partition=cpu" in script
    assert "--gres" not in script


def test_offline_mode_and_pinned_model_path(tmp_path):
    """The job never downloads: hub is offline, the model path is pinned."""
    req = _req(tmp_path)
    script = mj.render_model_sbatch(req, partition="gpu_a100")
    model = req.model
    assert "HF_HUB_OFFLINE=1" in script
    assert "TRANSFORMERS_OFFLINE=1" in script
    assert "HF_ENDPOINT=" in script  # the third-party mirror is overridden away
    assert model.revision in script and str(model.local_path) in script


def test_venv_name_hashed_from_requirements(tmp_path):
    req = _req(tmp_path)
    v1 = mj.venv_name(req.requirements)
    req2 = _req(tmp_path, requirements="faster-whisper==2.0.0\n")
    v2 = mj.venv_name(req2.requirements)
    assert len(v1) == 12 and v1.isalnum()
    assert v1 != v2, "a different pinned requirement set hashes to a different venv"
    script = mj.render_model_sbatch(req, partition="gpu_a100")
    assert f"venv-{v1}" in script


def test_run_record_is_written_to_the_job_folder(tmp_path):
    req = _req(tmp_path)
    script = mj.render_model_sbatch(req, partition="gpu_a100")
    assert "run_record.json" in script
    assert "-m delfin.dashboard.local_runner" not in script  # not an ORCA/MLP stencil


# --------------------------------------------------------------- dry run


def test_plan_is_a_dry_run_never_submits(tmp_path, monkeypatch):
    calls = []
    import delfin.scheduler_client as sc

    monkeypatch.setattr(sc, "query",
                        lambda kind, argv, runner=None: calls.append(argv) or "x")
    req = _req(tmp_path)
    plan = mj.model_job_plan(req, scontrol_output=FAKE_SCONTROL)
    assert calls == [], "a dry run must not touch the scheduler"
    assert plan["partition"] == "gpu_a100"
    assert "#SBATCH" in plan["script_text"]
    assert plan["submitted"] is False


def test_plan_refuses_when_no_partition_is_found(tmp_path):
    req = _req(tmp_path)
    with pytest.raises(mj.NoSuitablePartition):
        mj.model_job_plan(req, scontrol_output="")


# ------------------------------------------------------------- submission


def test_submission_goes_through_the_existing_gate(tmp_path, monkeypatch):
    """Real submission runs through scheduler_client + the job cap."""
    import delfin.scheduler_client as sc

    calls = []
    monkeypatch.setattr(
        sc, "query", lambda kind, argv, runner=None: calls.append(argv) or "Submitted batch job 17"
    )
    req = _req(tmp_path)
    result = mj.submit_model_job(req, scontrol_output=FAKE_SCONTROL, submit=True,
                                 active_jobs=0)
    assert result["submitted"] is True
    assert len(calls) == 1
    assert calls[0][0] == "sbatch"


def test_dry_run_default_submits_nothing(tmp_path, monkeypatch):
    import delfin.scheduler_client as sc

    calls = []
    monkeypatch.setattr(sc, "query",
                        lambda kind, argv, runner=None: calls.append(argv) or "x")
    req = _req(tmp_path)
    result = mj.submit_model_job(req, scontrol_output=FAKE_SCONTROL)
    assert result["submitted"] is False and calls == []


def test_submission_respects_the_concurrent_job_cap(tmp_path):
    with pytest.raises(cd.JobCapExceeded):
        mj._check_cap(active=MAX_CONCURRENT_MLP_JOBS)


# ---------------------------------------------- real-run findings (dev_cpu job)


def test_the_script_runs_under_bash_with_its_arguments(tmp_path):
    """A real job died at once: the run record's folder stood on its own
    line after the heredoc (run as a command, sys.argv[1] missing), and the
    completion stamp called a bare 'python'.  The arguments now come before
    the heredoc and both blocks use the job's venv."""
    import subprocess
    req = _req(tmp_path, gpu=False, partition="cpu")
    text = mj.render_model_sbatch(req, partition="cpu")
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line == "PYEOF":
            nxt = lines[i + 1] if i + 1 < len(lines) else ""
            assert not nxt.startswith('"'), f"argument after a heredoc: {nxt}"
    assert not any(line.startswith("python -") for line in lines)
    script = tmp_path / "s.sh"
    script.write_text(text)
    assert subprocess.run(["bash", "-n", str(script)]).returncode == 0


def test_logs_go_to_the_job_folder(tmp_path):
    req = _req(tmp_path, gpu=False, partition="cpu")
    text = mj.render_model_sbatch(req, partition="cpu")
    out = tmp_path / "out"
    assert f"#SBATCH --output={out}/whisper_job_%j.out" in text
    assert f"#SBATCH --error={out}/whisper_job_%j.err" in text


def test_the_cap_counts_only_model_jobs(monkeypatch, tmp_path):
    """A user with 21 MANTA build tasks queued could not submit one model
    job: every squeue line counted.  Only names under MODEL_JOB_PREFIX do."""
    from delfin import scheduler_client as sc
    text = mj.render_model_sbatch(_req(tmp_path, gpu=False, partition="cpu"),
                                  partition="cpu")
    assert f"#SBATCH --job-name={mj.MODEL_JOB_PREFIX}whisper_job" in text
    names = "\n".join(["manta_freeze2_cluster"] * 21 + ["model_whisper_job"])
    monkeypatch.setattr(sc, "query", lambda kind, argv, runner=None: names)
    assert mj.count_active_model_jobs() == 1
