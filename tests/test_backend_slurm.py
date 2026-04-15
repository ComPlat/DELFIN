from unittest.mock import patch

from delfin.dashboard.backend_slurm import SlurmJobBackend


def test_slurm_backend_appends_known_profile_env():
    backend = SlurmJobBackend("/tmp", slurm_profile="bwunicluster3")

    env_vars = backend._append_profile_env("DELFIN_MODE=delfin")

    assert "DELFIN_MODE=delfin" in env_vars
    assert "DELFIN_MODULES=devel/python/3.11.7-gnu-14.2" in env_vars
    assert "DELFIN_STAGE_ORCA=1" in env_vars
    assert "DELFIN_STAGE_VENV=1" in env_vars
    assert "DELFIN_RUNTIME_CACHE=1" in env_vars


def test_slurm_backend_leaves_env_unchanged_for_unknown_profile():
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    env_vars = backend._append_profile_env("DELFIN_MODE=delfin")

    assert env_vars == "DELFIN_MODE=delfin"


def test_detect_profile_bwunicluster3_fqdn():
    with patch("delfin.dashboard.backend_slurm.socket.getfqdn",
               return_value="uc3-login1.scc.kit.edu"):
        assert SlurmJobBackend._detect_profile() == "bwunicluster3"


def test_detect_profile_bwunicluster3_hostname():
    with patch("delfin.dashboard.backend_slurm.socket.getfqdn",
               return_value="uc3-0042"), \
         patch.dict("os.environ", {"HOSTNAME": "uc3-0042"}):
        assert SlurmJobBackend._detect_profile() == "bwunicluster3"


def test_detect_profile_unknown_site():
    with patch("delfin.dashboard.backend_slurm.socket.getfqdn",
               return_value="compute01.other-uni.de"), \
         patch.dict("os.environ", {"HOSTNAME": "compute01"}, clear=False):
        assert SlurmJobBackend._detect_profile() == ""


def test_auto_detect_fills_empty_profile():
    with patch.object(SlurmJobBackend, "_detect_profile",
                      return_value="bwunicluster3"):
        backend = SlurmJobBackend("/tmp")
        assert backend.slurm_profile == "bwunicluster3"


def test_explicit_profile_overrides_auto_detect():
    with patch.object(SlurmJobBackend, "_detect_profile",
                      return_value="bwunicluster3"):
        backend = SlurmJobBackend("/tmp", slurm_profile="custom")
        assert backend.slurm_profile == "custom"


def _fake_completed(stdout: str = "Submitted batch job 1\n"):
    class R:
        returncode = 0
        stderr = ""
    r = R()
    r.stdout = stdout
    return r


def test_submit_guppy_batch_emits_array_and_env(tmp_path):
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")
    csv = tmp_path / "smiles.csv"
    csv.write_text("smiles,name\n[Fe+2],iron\nO,water\n", encoding="utf-8")

    with patch("delfin.dashboard.backend_slurm.subprocess.run",
               return_value=_fake_completed()) as run_mock:
        result = backend.submit_guppy_batch(
            job_dir=str(tmp_path),
            job_name="guppy_bench",
            smiles_csv=str(csv),
            array_size=5,
            array_concurrency=2,
            pal=8,
            maxcore=4000,
            start_strategy="isomers",
            max_isomers=150,
            goat_topk=3,
        )

    assert result.returncode == 0
    args = run_mock.call_args.args[0]
    cmd_str = " ".join(args)

    # Array spec 1-5%2 (1-based so SLURM_ARRAY_TASK_ID matches CSV row).
    assert "--array=1-5%2" in args
    # Env vars: mode, csv path, strategy, isomer cap, goat topk, rmsd/energy window.
    assert "DELFIN_MODE=guppy_batch" in cmd_str
    assert f"GUPPY_BATCH_CSV={csv}" in cmd_str
    assert "GUPPY_START_STRATEGY=isomers" in cmd_str
    assert "GUPPY_MAX_ISOMERS=150" in cmd_str
    assert "GUPPY_GOAT_TOPK=3" in cmd_str
    assert "GUPPY_RMSD_CUTOFF=0.3" in cmd_str
    assert "GUPPY_ENERGY_WINDOW_KCAL=25.0" in cmd_str
    # Resource derivation.
    assert "--cpus-per-task=8" in args
    assert "--mem=32000M" in args
    assert any(a.endswith("submit_delfin.sh") for a in args)
    assert args[0] == "sbatch"
    assert "--job-name=guppy_bench" in args


def test_submit_guppy_batch_without_concurrency(tmp_path):
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")
    csv = tmp_path / "b.csv"
    csv.write_text("CCO\nO\n", encoding="utf-8")

    with patch("delfin.dashboard.backend_slurm.subprocess.run",
               return_value=_fake_completed()) as run_mock:
        backend.submit_guppy_batch(
            job_dir=str(tmp_path), job_name="x", smiles_csv=str(csv),
            array_size=2,
        )
    args = run_mock.call_args.args[0]
    assert "--array=1-2" in args
    assert not any(a.startswith("--array=1-2%") for a in args)
