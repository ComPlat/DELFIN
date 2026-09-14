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


def test_an_unknown_site_still_runs_the_venv_from_node_local_disk(monkeypatch):
    """Staging the venv is not a bwUniCluster feature.

    Python started from a venv on a network HOME is an I/O problem on any
    cluster. A site without a profile gets the staging, and nothing that only
    one site has -- no module names, no node sizes.
    """
    for key in ("DELFIN_STAGE_VENV", "DELFIN_RUNTIME_CACHE"):
        monkeypatch.delenv(key, raising=False)
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    env_vars = backend._append_profile_env("DELFIN_MODE=delfin")

    assert env_vars.startswith("DELFIN_MODE=delfin,")
    assert "DELFIN_STAGE_VENV=1" in env_vars
    assert "DELFIN_RUNTIME_CACHE=1" in env_vars
    assert "DELFIN_MODULES" not in env_vars
    assert "DELFIN_NODE_CORES" not in env_vars


def test_a_staging_choice_the_user_exported_is_theirs(monkeypatch):
    monkeypatch.setenv("DELFIN_STAGE_VENV", "0")
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    env_vars = backend._append_profile_env("DELFIN_MODE=delfin")

    assert "DELFIN_STAGE_VENV=1" not in env_vars


def test_the_job_is_told_which_python_environment_submitted_it(monkeypatch, tmp_path):
    """The job used to look for ``software/delfin`` above its submit directory."""
    import delfin.dashboard.backend_slurm as backend_slurm

    for key in ("DELFIN_VENV", "DELFIN_REPO", "DELFIN_OMPI_HOME"):
        monkeypatch.delenv(key, raising=False)
    venv = tmp_path / "anywhere" / "env"
    monkeypatch.setattr(backend_slurm.sys, "prefix", str(venv))
    monkeypatch.setattr(backend_slurm.sys, "base_prefix", "/usr")
    monkeypatch.setattr(backend_slurm.shutil, "which", lambda name: None)

    location = SlurmJobBackend._runtime_location_env()

    assert location["DELFIN_VENV"] == str(venv)
    assert "DELFIN_OMPI_HOME" not in location


def test_a_system_mpirun_is_not_staged(monkeypatch):
    """The job copies DELFIN_OMPI_HOME to local disk; /usr is not that."""
    import delfin.dashboard.backend_slurm as backend_slurm

    monkeypatch.delenv("DELFIN_OMPI_HOME", raising=False)
    monkeypatch.setattr(backend_slurm.shutil, "which", lambda name: "/usr/bin/mpirun")

    assert "DELFIN_OMPI_HOME" not in SlurmJobBackend._runtime_location_env()


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


def test_sbatch_passes_vars_via_process_env_not_export_flag(tmp_path):
    """DELFIN vars must travel via the sbatch process environment.

    Regression: ``--export=ALL,VAR=…`` (explicit variables) makes slurmd
    fetch the user's login environment on the compute node at launch;
    on bwUniCluster3 that intermittently exceeds GetEnvTimeout (2 s) and
    the job is requeued/held with "user env retrieval failed". Verified
    by an A/B test pinned to the same node: plain/ALL-only completed,
    ALL,VARS failed repeatedly (jobs 5897358 vs 5897359/5897343).
    """
    backend = SlurmJobBackend("/tmp", slurm_profile="bwunicluster3")

    captured = {}

    def fake_run(cmd, *args, **kwargs):
        captured["cmd"] = cmd
        captured["env"] = kwargs.get("env")
        return _fake_completed("Submitted batch job 99\n")

    with patch("delfin.dashboard.backend_slurm.subprocess.run", side_effect=fake_run), \
         patch.object(SlurmJobBackend, "_release_env_hold"):
        backend.submit_delfin(
            job_dir=str(tmp_path), job_name="myjob", mode="delfin",
            time_limit="24:00:00", pal=4, maxcore=1000,
        )

    cmd = captured["cmd"]
    env = captured["env"]
    # Exactly --export=ALL — never with an appended variable list.
    assert "--export=ALL" in cmd
    assert not any(a.startswith("--export=ALL,") for a in cmd), cmd
    # The DELFIN variables (incl. profile env) must be in the process env.
    assert env is not None
    assert env.get("DELFIN_MODE") == "delfin"
    assert env.get("DELFIN_JOB_NAME") == "myjob"
    assert env.get("DELFIN_STAGE_ORCA") == "1"  # bwunicluster3 profile env
    # And the surrounding user environment must still be inherited.
    assert "PATH" in env


def test_sbatch_strips_inherited_sbatch_input_env_vars(tmp_path):
    """Inherited SBATCH_GET_USER_ENV / SBATCH_EXPORT must be stripped.

    sbatch reads SBATCH_* input environment variables as if they were CLI
    options; SBATCH_GET_USER_ENV has no command-line negation, so a shell
    that exports it would re-enable the compute-node env retrieval despite
    our --export=ALL. Other SBATCH_* vars (e.g. SBATCH_EXCLUDE) pass through.
    """
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster")

    captured = {}

    def fake_run(cmd, *args, **kwargs):
        captured["env"] = kwargs.get("env")
        return _fake_completed("Submitted batch job 100\n")

    with patch.dict("os.environ", {
            "SBATCH_GET_USER_ENV": "1",
            "SBATCH_EXPORT": "NONE",
            "SBATCH_EXCLUDE": "badnode01"}), \
         patch("delfin.dashboard.backend_slurm.subprocess.run",
               side_effect=fake_run), \
         patch.object(SlurmJobBackend, "_release_env_hold"):
        backend.submit_delfin(
            job_dir=str(tmp_path), job_name="x", mode="delfin",
            time_limit="24:00:00", pal=1, maxcore=600,
        )

    env = captured["env"]
    assert "SBATCH_GET_USER_ENV" not in env
    assert "SBATCH_EXPORT" not in env
    assert env.get("SBATCH_EXCLUDE") == "badnode01"  # benign ones survive


def test_env_vars_to_dict_parses_export_grammar():
    parse = SlurmJobBackend._env_vars_to_dict
    assert parse("A=1,B=two,C=/some/path") == {
        "A": "1", "B": "two", "C": "/some/path"}
    # First '=' splits key/value; later '=' stay in the value.
    assert parse("OPTS=a=b") == {"OPTS": "a=b"}
    # Malformed chunks are skipped, not fatal.
    assert parse("A=1,,novalue,B=2") == {"A": "1", "B": "2"}
    assert parse("") == {}


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


# ---------------------------------------------------------------------------
# why a job waits, and when it should start
# ---------------------------------------------------------------------------
from delfin.dashboard.backend_slurm import describe_pending_reason  # noqa: E402


def _answer(stdout="", stderr="", returncode=0):
    class R:
        pass
    r = R()
    r.stdout, r.stderr, r.returncode = stdout, stderr, returncode
    return r


def test_one_squeue_carries_why_a_job_waits_and_when_it_should_start():
    """The start estimate used to need a second, unthrottled squeue --start."""
    backend = SlurmJobBackend("/tmp", slurm_profile="custom-cluster",
                              tool_binaries={"orca": "/opt/orca"})
    out = (
        "JOBID|PARTITION|NAME|USER|ST|TIME|NODES|TIME_LIMIT|CPUS|MIN_MEMORY|START_TIME|REASON|NODELIST(REASON)\n"
        "900|cpu,cpu_il|TADF 1|u|PD|0:00|1|2-00:00:00|40|240000M|2026-09-16T21:20:09|Priority|(Priority)\n"
        "901|cpu|running|u|R|1:02:03|1|1-00:00:00|40|240000M|2026-09-14T08:00:00|None|uc3n024\n"
        "902|cpu|fresh|u|PD|0:00|1|1-00:00:00|40|240000M|N/A|None|(None)\n"
    )
    calls = []

    def fake_run(cmd, *args, **kwargs):
        calls.append(list(cmd))
        return _answer(out)

    with patch("delfin.dashboard.backend_slurm.subprocess.run", side_effect=fake_run):
        jobs = backend.list_jobs()
        pending = backend.get_pending_start_times()

    assert len([c for c in calls if c[0] == "squeue"]) == 1
    by_id = {j.job_id: j for j in jobs}
    assert by_id["900"].name == "TADF 1"
    assert by_id["900"].extra["cpus"] == "40"
    assert by_id["900"].extra["time_limit"] == "2-00:00:00"
    assert by_id["901"].extra["reason"] == "uc3n024"
    assert by_id["902"].extra["start_estimate"] == "", "N/A is no estimate"
    assert pending[0] == {
        "id": "900", "name": "TADF 1", "start": "2026-09-16T21:20:09",
        "reason": "Priority", "explanation": describe_pending_reason("Priority"),
        "partition": "cpu,cpu_il",
    }
    assert [p["id"] for p in pending] == ["900", "902"]


def test_a_waiting_reason_is_explained_in_words():
    assert "higher priority" in describe_pending_reason("(Priority)")
    assert "free cores" in describe_pending_reason("Resources")
    assert "maintenance" in describe_pending_reason("ReqNodeNotAvail, Reserved for maintenance")
    assert "QOS" in describe_pending_reason("QOSMaxCpuPerUserLimit")
    assert "account" in describe_pending_reason("AssocGrpCpuLimit")
    assert "releases" in describe_pending_reason("user env retrieval failed requeued held")
    assert describe_pending_reason("SomethingNew") == "SomethingNew"
    assert describe_pending_reason("") == ""


def _submitting_backend(monkeypatch, fits, **kwargs):
    """A backend whose sbatch --test-only answers from `fits(partition)`."""
    monkeypatch.delenv("DELFIN_SLURM_PARTITIONS", raising=False)
    backend = SlurmJobBackend("/tmp", tool_binaries={"orca": "/opt/orca"}, **kwargs)
    calls = []

    def fake_run(cmd, *args, **kw):
        calls.append(list(cmd))
        if "--test-only" in cmd:
            part = next(c.split("=", 1)[1] for c in cmd if c.startswith("--partition="))
            verdict = fits(part)
            if verdict is None:
                raise OSError("slurm controller unreachable")
            if verdict:
                return _answer(stderr=f"sbatch: Job 1 to start at 2026-09-16T21:20:09 using "
                                      f"80 processors on nodes n1 in partition {part}")
            return _answer(stderr="sbatch: error: allocation failure: Requested node "
                                  "configuration is not available", returncode=1)
        return _answer("Submitted batch job 5\n")

    monkeypatch.setattr("delfin.dashboard.backend_slurm.subprocess.run", fake_run)
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold", lambda self, out: None)
    return backend, calls


def _final_sbatch(calls):
    return [c for c in calls if c[0] == "sbatch" and "--test-only" not in c][-1]


def _submit(backend, tmp_path, maxcore=6000):
    return backend.submit_orca(str(tmp_path), "job", "missing.inp",
                               time_limit="2-00:00:00", pal=40, maxcore=maxcore)


def test_a_job_is_listed_for_every_partition_it_fits(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: True,
                                         slurm_profile="bwunicluster3")
    _submit(backend, tmp_path)
    assert "--partition=cpu,cpu_il" in _final_sbatch(calls)


def test_a_partition_the_request_does_not_fit_is_left_out(monkeypatch, tmp_path):
    """Listed anyway, SLURM rejects the whole job."""
    backend, calls = _submitting_backend(monkeypatch, lambda p: p == "cpu",
                                         slurm_profile="bwunicluster3")
    _submit(backend, tmp_path)
    assert "--partition=cpu" in _final_sbatch(calls)


def test_when_slurm_cannot_be_asked_the_template_decides(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: None,
                                         slurm_profile="bwunicluster3")
    result = _submit(backend, tmp_path)
    assert result.returncode == 0
    assert not any(c.startswith("--partition") for c in _final_sbatch(calls))


def test_the_answer_is_asked_once_per_request_shape(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: True,
                                         slurm_profile="bwunicluster3")
    for _ in range(3):
        _submit(backend, tmp_path)
    assert len([c for c in calls if "--test-only" in c]) == 2
    _submit(backend, tmp_path, maxcore=7000)
    assert len([c for c in calls if "--test-only" in c]) == 4, "a new shape is a new question"


def test_the_setting_wins_over_the_site_profile(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: True,
                                         slurm_profile="bwunicluster3", partitions="cpu_il")
    assert backend.partitions == ("cpu_il",)
    _submit(backend, tmp_path)
    assert "--partition=cpu_il" in _final_sbatch(calls)
    assert not [c for c in calls if "--test-only" in c], "one partition needs no question"


def test_a_site_without_a_profile_keeps_the_template_partition(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: True,
                                         slurm_profile="custom-cluster")
    assert backend.partitions == ()
    _submit(backend, tmp_path)
    assert not [c for c in calls if "--test-only" in c]
    assert not any(c.startswith("--partition") for c in _final_sbatch(calls))


def test_a_partition_given_by_the_caller_is_not_second_guessed(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: True,
                                         slurm_profile="bwunicluster3")
    backend._sbatch(str(tmp_path), "DELFIN_MODE=mlp", "01:00:00", 4, 16000, "gpu",
                    tmp_path / "submit.sh", gpu="gpu:1", partition="gpu_h100")
    assert not [c for c in calls if "--test-only" in c]
    assert "--partition=gpu_h100" in _final_sbatch(calls)


def test_every_slurm_time_format_is_understood():
    """A limit in days took the submit down with ValueError('2-00')."""
    seconds = SlurmJobBackend._time_limit_seconds
    assert seconds("2-00:00:00") == 2 * 86400
    assert seconds("1-12") == 86400 + 12 * 3600
    assert seconds("3-01:30") == 3 * 86400 + 3600 + 30 * 60
    assert seconds("48:00:00") == 48 * 3600
    assert seconds("30:00") == 30 * 60
    assert seconds("90") == 90 * 60


# ---------------------------------------------------------------------------
# what the dashboard's submit paths hand over
# ---------------------------------------------------------------------------
from delfin.dashboard.backend_slurm import normalize_time_limit  # noqa: E402


def test_an_orca_job_is_sized_by_the_input_it_runs_not_the_control_beside_it(monkeypatch, tmp_path):
    """A recalc edited down to 12 processes still reserved CONTROL's 40 cores."""
    (tmp_path / "CONTROL.txt").write_text("PAL=40\nmaxcore=6000\n")
    (tmp_path / "opt_recalc_1.inp").write_text("! PBE0 def2-SVP\n%pal nprocs 12 end\n%maxcore 3000\n")
    backend, calls = _submitting_backend(monkeypatch, lambda p: True, slurm_profile="custom-cluster")

    backend.submit_orca(str(tmp_path), "opt_recalc_1", "opt_recalc_1.inp",
                        time_limit="24:00:00", pal=40, maxcore=6000)

    final = _final_sbatch(calls)
    assert "--cpus-per-task=12" in final and "--mem=36000M" in final


def test_a_delfin_job_is_still_sized_by_its_control(monkeypatch, tmp_path):
    (tmp_path / "CONTROL.txt").write_text("PAL=24\nmaxcore=4000\n")
    (tmp_path / "leftover.inp").write_text("%pal nprocs 8 end\n%maxcore 1000\n")
    backend, calls = _submitting_backend(monkeypatch, lambda p: True, slurm_profile="custom-cluster")

    backend.submit_delfin(str(tmp_path), "run", mode="delfin", time_limit="24:00:00")

    final = _final_sbatch(calls)
    assert "--cpus-per-task=24" in final and "--mem=96000M" in final


def test_a_time_typed_the_way_people_write_it_is_understood():
    assert normalize_time_limit("48:00:00") == "48:00:00"
    assert normalize_time_limit(" 2-00:00:00 ") == "2-00:00:00"
    assert normalize_time_limit("90") == "90"
    assert normalize_time_limit("48h") == "2-00:00:00"
    assert normalize_time_limit("1d12h") == "1-12:00:00"
    assert normalize_time_limit("90min") == "01:30:00"
    assert normalize_time_limit("2 h 30 min") == "02:30:00"
    for bad in ("", "48 hours", "abc", "0", "1-2-3"):
        try:
            normalize_time_limit(bad)
        except ValueError:
            continue
        raise AssertionError(f"{bad!r} was accepted")


def test_an_unusable_time_limit_is_refused_with_a_reason_not_a_traceback(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(monkeypatch, lambda p: True, slurm_profile="bwunicluster3")

    result = backend.submit_orca(str(tmp_path), "job", "missing.inp", time_limit="48 hours")

    assert result.returncode == 1
    assert "Invalid time limit '48 hours'" in result.stderr and "48h" in result.stderr
    assert not [c for c in calls if c and c[0] == "sbatch"], "nothing reached SLURM"

    backend.submit_orca(str(tmp_path), "job", "missing.inp", time_limit="48h")
    assert "--time=2-00:00:00" in _final_sbatch(calls)


# ---------------------------------------------------------------------------
# GPU jobs
# ---------------------------------------------------------------------------
from delfin.slurm_submit import discover_gpu_partitions  # noqa: E402

_SCONTROL_ONELINER = (
    "PartitionName=cpu State=UP TRES=cpu=15360,mem=30920000M,node=80,billing=15360\n"
    "PartitionName=dev_gpu_h100 State=UP TRES=cpu=192,mem=773500M,node=1,billing=192,gres/gpu=4\n"
    "PartitionName=gpu_h100 State=UP TRES=cpu=2304,mem=9282000M,node=12,billing=2304,gres/gpu=48\n"
    "PartitionName=gpu_old State=DOWN TRES=cpu=96,mem=380000M,node=1,billing=96,gres/gpu=4\n"
    "PartitionName=gpu_a100_il State=UP TRES=cpu=1152,mem=4590000M,node=9,billing=1152,gres/gpu=36\n"
)


def test_gpu_partitions_are_found_without_sinfo():
    """bwUniCluster refuses sinfo; scontrol names the partitions with GPUs."""
    assert discover_gpu_partitions(_SCONTROL_ONELINER) == ("gpu_h100", "gpu_a100_il")


def test_a_gpu_job_is_listed_for_every_gpu_partition_it_fits(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(
        monkeypatch, lambda p: p in {"gpu_h100", "gpu_a100_il"}, slurm_profile="bwunicluster3")
    monkeypatch.delenv("DELFIN_SLURM_GPU_PARTITIONS", raising=False)

    backend.submit_mlp(str(tmp_path), "mlp_job", "mol.xyz", time_limit="24:00:00", pal=4, maxcore=4000)

    final = _final_sbatch(calls)
    assert "--gres=gpu:1" in final
    assert "--partition=gpu_h100,gpu_a100_il" in final
    asked = [c for c in calls if "--test-only" in c]
    assert asked and all("--gres=gpu:1" in c for c in asked), "a GPU partition is asked about a GPU job"
    assert not [c for c in calls if c and c[0] == "sinfo"]


def test_a_gpu_job_that_no_gpu_partition_can_run_runs_on_cpus(monkeypatch, tmp_path):
    backend, calls = _submitting_backend(
        monkeypatch, lambda p: p in {"cpu", "cpu_il"}, slurm_profile="bwunicluster3")
    monkeypatch.delenv("DELFIN_SLURM_GPU_PARTITIONS", raising=False)

    backend.submit_mlp(str(tmp_path), "mlp_job", "mol.xyz", time_limit="24:00:00", pal=4, maxcore=4000)

    final = _final_sbatch(calls)
    assert not any(c.startswith("--gres") for c in final)
    assert "--partition=cpu,cpu_il" in final


def test_the_gpu_setting_wins_and_a_site_without_one_asks_scontrol(monkeypatch, tmp_path):
    import delfin.dashboard.backend_slurm as backend_slurm

    backend, calls = _submitting_backend(monkeypatch, lambda p: True,
                                         slurm_profile="bwunicluster3", gpu_partitions="gpu_h100_il")
    backend.submit_mlp(str(tmp_path), "mlp_job", "mol.xyz", time_limit="02:00:00")
    assert "--partition=gpu_h100_il" in _final_sbatch(calls)

    monkeypatch.delenv("DELFIN_SLURM_GPU_PARTITIONS", raising=False)
    monkeypatch.setattr(backend_slurm, "discover_gpu_partitions", lambda: ("gpu_x",))
    elsewhere, more = _submitting_backend(monkeypatch, lambda p: True, slurm_profile="custom-cluster")
    elsewhere.submit_mlp(str(tmp_path), "mlp_job", "mol.xyz", time_limit="02:00:00")
    assert "--partition=gpu_x" in _final_sbatch(more) and "--gres=gpu:1" in _final_sbatch(more)
