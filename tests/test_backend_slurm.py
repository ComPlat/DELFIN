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


def _result(stdout="", returncode=0):
    class R:
        pass
    r = R()
    r.stdout = stdout
    r.stderr = ""
    r.returncode = returncode
    return r


def test_release_env_hold_releases_held_job(tmp_path):
    """A job SLURM held for failed user-env retrieval is auto-released."""
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    calls = []

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)
        if cmd[0] == "sbatch":
            return _result("Submitted batch job 42\n")
        if cmd[:2] == ["scontrol", "show"]:
            # First check: held with the env-retrieval reason.
            return _result(
                "JobId=42 JobState=PENDING "
                "Reason=user_env_retrieval_failed_requeued_held Requeue=1"
            )
        return _result("")  # scontrol release

    with patch("delfin.dashboard.backend_slurm.subprocess.run", side_effect=fake_run), \
         patch("delfin.dashboard.backend_slurm.time.sleep"):
        result = backend.submit_delfin(
            job_dir=str(tmp_path), job_name="x", mode="delfin",
            time_limit="24:00:00", pal=1, maxcore=600,
        )

    # Job id parsing must stay intact (last stdout token = the id).
    assert result.stdout.strip().split()[-1] == "42"
    assert ["scontrol", "release", "42"] in calls, (
        f"expected an scontrol release for the held job, got: {calls}"
    )


def test_release_env_hold_leaves_normal_job_alone(tmp_path):
    """A normally-pending job (Reason=Priority) is never released."""
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    calls = []

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)
        if cmd[0] == "sbatch":
            return _result("Submitted batch job 7\n")
        if cmd[:2] == ["scontrol", "show"]:
            return _result("JobId=7 JobState=PENDING Reason=Priority")
        return _result("")

    with patch("delfin.dashboard.backend_slurm.subprocess.run", side_effect=fake_run), \
         patch("delfin.dashboard.backend_slurm.time.sleep"):
        backend.submit_delfin(
            job_dir=str(tmp_path), job_name="x", mode="delfin",
            time_limit="24:00:00", pal=1, maxcore=600,
        )

    assert not any(c[:2] == ["scontrol", "release"] for c in calls), (
        f"must not release a normally-pending job, got: {calls}"
    )


def test_list_jobs_releases_env_hold_on_refresh():
    """A dispatch-time env-retrieval hold seen in squeue is auto-released,
    while a normally-pending job on the same refresh is left alone."""
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    squeue_out = (
        "       JOBID    PARTITION       NAME       USER  ST         TIME  "
        "NODES NODELIST(REASON)\n"
        "     5880813          cpu    stuckjob  ka_ew7404  PD         0:00  "
        "    1 (user env retrieval failed requeued held)\n"
        "     5880814          cpu   normaljob  ka_ew7404  PD         0:00  "
        "    1 (Priority)\n"
    )
    calls = []

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)
        if cmd[0] == "squeue":
            return _result(squeue_out)
        return _result("")  # scontrol release

    with patch("delfin.dashboard.backend_slurm.subprocess.run", side_effect=fake_run):
        jobs = backend.list_jobs()

    assert {j.job_id for j in jobs} == {"5880813", "5880814"}
    releases = [c for c in calls if c[:2] == ["scontrol", "release"]]
    assert releases == [["scontrol", "release", "5880813"]], (
        f"expected exactly the held job released, got: {releases}"
    )


def test_is_env_hold_matches_both_scontrol_and_squeue_forms():
    """The detector must match underscores (scontrol) and spaces (squeue)."""
    assert SlurmJobBackend._is_env_hold(
        "user_env_retrieval_failed_requeued_held")
    assert SlurmJobBackend._is_env_hold(
        "(user env retrieval failed requeued held)")
    assert not SlurmJobBackend._is_env_hold("Priority")
    assert not SlurmJobBackend._is_env_hold("(Resources)")
    assert not SlurmJobBackend._is_env_hold("JobHeldUser")


def test_submit_delfin_recalc_extracts_inline_pal_from_inp(tmp_path):
    """Recalc submit must derive pal=N from inline ``%pal nprocs N end``.

    Regression: same root cause as the ORCA Builder OOM (commit 12c6ed1).
    The recalc submit path (tab_calculations_browser._on_submit_recalc and
    backend_slurm._resolve_resources) both feed into parse_inp_resources.
    Before the fix the regex was line-anchored, missed the inline form, and
    SLURM fell back to the widget defaults (12/6000) -> --mem=72G while the
    .inp directed ORCA to nprocs=40 -> OOM in PROPINT
    (RSS_paper_85_NEB-TS_recalc_1, job 4151595).
    """
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")
    (tmp_path / "RSS_paper_85_NEB-TS_recalc_1.inp").write_text(
        "!PBE0 def2-SVP def2/J D4 RIJCOSX CPCM(DMF) NEB-TS Freq PModel\n\n"
        "%maxcore 6000\n"
        "%pal nprocs 40 end\n"
        "%scf maxiter 225 end\n"
        "*XYZFILE -2 3 scan.021.xyz\n",
        encoding="utf-8",
    )

    with patch("delfin.dashboard.backend_slurm.subprocess.run",
               return_value=_fake_completed()) as run_mock:
        backend.submit_delfin(
            job_dir=str(tmp_path),
            job_name="recalc_RSS_paper_85_NEB-TS_2",
            mode="delfin-recalc-classic",
            time_limit="48:00:00",
            # Caller passes widget-default fallback values; _resolve_resources
            # MUST override them from the .inp's inline %pal block.
            pal=12,
            maxcore=6000,
        )

    # _sbatch now also calls subprocess.run for the post-submit env-hold
    # check, so pick the sbatch invocation explicitly (not the last call).
    sbatch_calls = [c.args[0] for c in run_mock.call_args_list
                    if c.args and c.args[0] and c.args[0][0] == 'sbatch']
    args = sbatch_calls[-1]
    assert "--cpus-per-task=40" in args, (
        f"Expected pal=40 from inline %pal nprocs 40 end, got cmd: {args}"
    )
    assert "--mem=240000M" in args, (
        f"Expected mem=40*6000=240000M, got cmd: {args}"
    )


def test_submit_delfin_recalc_extracts_multiline_pal_from_inp(tmp_path):
    """Multi-line %pal block must also be extracted by the recalc path."""
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")
    (tmp_path / "job_recalc_1.inp").write_text(
        "!PBE0 def2-SVP\n"
        "%pal\n  nprocs 24\nend\n"
        "%maxcore 4000\n"
        "*xyz 0 1\nH 0 0 0\n*\n",
        encoding="utf-8",
    )

    with patch("delfin.dashboard.backend_slurm.subprocess.run",
               return_value=_fake_completed()) as run_mock:
        backend.submit_delfin(
            job_dir=str(tmp_path), job_name="x",
            mode="delfin-recalc-classic",
            time_limit="24:00:00", pal=12, maxcore=6000,
        )

    # _sbatch now also calls subprocess.run for the post-submit env-hold
    # check, so pick the sbatch invocation explicitly (not the last call).
    sbatch_calls = [c.args[0] for c in run_mock.call_args_list
                    if c.args and c.args[0] and c.args[0][0] == 'sbatch']
    args = sbatch_calls[-1]
    assert "--cpus-per-task=24" in args
    assert "--mem=96000M" in args


def test_submit_delfin_recalc_extracts_pal_keyword_shortcut(tmp_path):
    """``! PAL8`` keyword shortcut is the third valid form."""
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")
    (tmp_path / "job_recalc_1.inp").write_text(
        "! PBE0 def2-SVP PAL8\n%maxcore 3000\n*xyz 0 1\nH 0 0 0\n*\n",
        encoding="utf-8",
    )

    with patch("delfin.dashboard.backend_slurm.subprocess.run",
               return_value=_fake_completed()) as run_mock:
        backend.submit_delfin(
            job_dir=str(tmp_path), job_name="x",
            mode="delfin-recalc-classic",
            time_limit="24:00:00", pal=12, maxcore=6000,
        )

    # _sbatch now also calls subprocess.run for the post-submit env-hold
    # check, so pick the sbatch invocation explicitly (not the last call).
    sbatch_calls = [c.args[0] for c in run_mock.call_args_list
                    if c.args and c.args[0] and c.args[0][0] == 'sbatch']
    args = sbatch_calls[-1]
    assert "--cpus-per-task=8" in args
    assert "--mem=24000M" in args
