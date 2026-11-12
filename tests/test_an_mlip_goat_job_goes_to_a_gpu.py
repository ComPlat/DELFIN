"""A DELFIN job whose CONTROL.txt asks for an MLIP GOAT gets a GPU by itself.

``global_opt_method=mlip:<backend>`` (device cuda, the default) runs GOAT
inside the DELFIN job on that job's GPU; the job stops at the GOAT step on a
node without one.  The live ethanol run needed ``--gres=gpu:1`` and a GPU
partition added by hand.  Now ``submit_delfin`` -- the dashboard's, the
agent's and the API's submission -- sends such a job to the configured GPU
partitions with one GPU: those too short for the time limit left out with
one scontrol query, the rest asked with sbatch --test-only like every
multi-partition submit; PAL, maxcore and the time limit untouched; and a
refusal, not a CPU job, when no GPU partition can take it.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.dashboard.backend_slurm import SlurmJobBackend

GOAT_MLIP = ("PAL=8\nmaxcore=2000\nglobal_optimizer=GOAT\n"
             "global_opt_method=mlip:aimnet2\nglobal_opt_mlip_device=cuda\n")
LIMITS = {"dev_gpu": "00:30:00", "gpu_long": "3-00:00:00", "gpu_mid": "2-00:00:00"}


def _answer(stdout="", stderr="", returncode=0):
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=stdout, stderr=stderr)


@pytest.fixture()
def backend(monkeypatch):
    """A backend whose sbatch --test-only accepts every partition and whose
    scontrol answers from LIMITS; returns (backend, sbatch calls, scontrol calls)."""
    import delfin.mlp_tools.extopt_submit as es

    for var in ("DELFIN_SLURM_PARTITIONS", "DELFIN_SLURM_GPU_PARTITIONS",
                "DELFIN_MLP_CLUSTER_GPU_PARTITIONS"):
        monkeypatch.delenv(var, raising=False)
    calls, asked = [], []

    def fake_run(cmd, *args, **kw):
        calls.append(list(cmd))
        if "--test-only" in cmd:
            part = next(c.split("=", 1)[1] for c in cmd if c.startswith("--partition="))
            return _answer(stderr=f"sbatch: Job 1 to start at 2026-11-12T16:00:00 using 8 "
                                  f"processors on nodes n1 in partition {part}")
        return _answer("Submitted batch job 5\n")

    def fake_scontrol(argv):
        asked.append(list(argv))
        return "\n".join(f"PartitionName={n} MaxTime={t} State=UP" for n, t in LIMITS.items())

    monkeypatch.setattr("delfin.dashboard.backend_slurm.subprocess.run", fake_run)
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold", lambda self, out: None)
    monkeypatch.setattr(es, "_ask_scheduler", fake_scontrol)

    def make(**kwargs):
        return SlurmJobBackend("/tmp", tool_binaries={"orca": "/opt/orca"},
                               slurm_profile="custom-cluster", **kwargs)
    return make, calls, asked


def _final(calls):
    sbatch = [c for c in calls if c and c[0] == "sbatch" and "--test-only" not in c]
    return sbatch[-1] if sbatch else None


def test_the_control_says_whether_goat_runs_an_mlip_on_a_gpu():
    from delfin.mlp_tools.goat_gpu_submit import mlip_goat_method

    assert mlip_goat_method(GOAT_MLIP) == "mlip:aimnet2"
    # cuda is the default device
    assert mlip_goat_method("global_optimizer=GOAT\nglobal_opt_method=mlip:orb\n") == "mlip:orb"
    # the older XTB_GOAT switch, without a global_optimizer choice
    assert mlip_goat_method("XTB_GOAT=yes\nglobal_opt_method=mlip:orb\n") == "mlip:orb"
    assert mlip_goat_method("global_optimizer=[GOAT|CREST]\nXTB_GOAT=yes\n"
                            "global_opt_method=mlip:orb\n") == "mlip:orb"
    # not GPU runs
    assert mlip_goat_method(GOAT_MLIP.replace("device=cuda", "device=cpu")) == ""
    assert mlip_goat_method("global_optimizer=GOAT\nglobal_opt_method=\n") == ""
    assert mlip_goat_method("global_optimizer=CREST\nglobal_opt_method=mlip:orb\n") == ""
    assert mlip_goat_method("global_opt_method=mlip:nonsense\nglobal_optimizer=GOAT\n") == ""


def test_an_mlip_goat_job_goes_to_the_gpu_partitions_with_one_gpu(backend, tmp_path):
    make, calls, asked = backend
    (tmp_path / "CONTROL.txt").write_text(GOAT_MLIP)

    result = make(gpu_partitions="dev_gpu,gpu_long,gpu_mid").submit_delfin(
        str(tmp_path), "goat", mode="delfin", time_limit="04:00:00")

    assert result.returncode == 0
    final = _final(calls)
    assert "--gres=gpu:1" in final
    # dev_gpu (MaxTime 30 min) cannot hold 04:00:00; the others both fit.
    assert "--partition=gpu_long,gpu_mid" in final
    # PAL, maxcore and the time limit as asked.
    assert "--time=04:00:00" in final and "--cpus-per-task=8" in final and "--mem=16000M" in final
    # One scontrol query for every partition.
    assert asked == [["scontrol", "show", "partition", "-o"]]
    tested = [c for c in calls if "--test-only" in c]
    assert tested and all("--gres=gpu:1" in c for c in tested)
    assert not any("--partition=dev_gpu" in c for c in tested)


def test_the_mlp_variable_names_the_gpu_partitions_too(backend, tmp_path, monkeypatch):
    make, calls, _ = backend
    (tmp_path / "CONTROL.txt").write_text(GOAT_MLIP)
    monkeypatch.setenv("DELFIN_MLP_CLUSTER_GPU_PARTITIONS", "gpu_mid")

    make().submit_delfin(str(tmp_path), "goat", mode="delfin", time_limit="01:00:00")

    final = _final(calls)
    assert "--partition=gpu_mid" in final and "--gres=gpu:1" in final


def test_no_gpu_partition_that_holds_the_time_is_a_refusal_not_a_cpu_job(backend, tmp_path):
    make, calls, _ = backend
    (tmp_path / "CONTROL.txt").write_text(GOAT_MLIP)

    result = make(gpu_partitions="dev_gpu").submit_delfin(
        str(tmp_path), "goat", mode="delfin", time_limit="04:00:00")

    assert result.returncode == 1
    assert "dev_gpu (MaxTime 00:30:00)" in result.stderr
    assert "global_opt_method=mlip:aimnet2" in result.stderr
    assert _final(calls) is None, "nothing was submitted"


def test_no_gpu_partition_configured_is_a_refusal(backend, tmp_path, monkeypatch):
    import delfin.dashboard.backend_slurm as backend_slurm

    make, calls, _ = backend
    (tmp_path / "CONTROL.txt").write_text(GOAT_MLIP)
    monkeypatch.setattr(backend_slurm, "discover_gpu_partitions", lambda: ())

    result = make().submit_delfin(str(tmp_path), "goat", mode="delfin", time_limit="00:30:00")

    assert result.returncode == 1
    assert "no GPU partition is configured" in result.stderr
    assert "no CPU fallback" in result.stderr
    assert _final(calls) is None


def test_a_classic_job_is_sent_as_before(backend, tmp_path):
    make, calls, asked = backend
    (tmp_path / "CONTROL.txt").write_text("PAL=8\nmaxcore=2000\nglobal_optimizer=GOAT\n")

    make(gpu_partitions="gpu_long").submit_delfin(
        str(tmp_path), "goat", mode="delfin", time_limit="04:00:00")

    final = _final(calls)
    assert not any(c.startswith("--gres") for c in final)
    assert not asked
