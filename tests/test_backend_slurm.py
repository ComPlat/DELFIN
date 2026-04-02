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
