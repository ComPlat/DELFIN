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
