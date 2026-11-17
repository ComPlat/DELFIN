"""The guided MLIP set-up: probe, plan, install with approval, configure, verify.

Input: a host described by fakes (``guided_setup.Host``): no Slurm, Slurm
with GPU partitions, an existing MLIP environment on disk. No command of
the plan is run for real; the runner, the approver and the scheduler are
recorders.

Output under test:
* the probe never raises and states what is absent;
* the plan picks a CUDA-capable partition, the torch build for the probed
  driver, and carries each backend's code and weight licence records;
* no install command runs without an approval for that command, and none
  runs when nobody can be asked;
* settings are written through save_settings only with write=True, after
  the diff, and survive the round trip;
* the verify job is rendered for the partition and submitted through
  scheduler_client under the MLP job cap; its result reader reports
  pass/fail with the numbers and the likely cause;
* the UMA path never reads a token and no output carries one.
"""

from __future__ import annotations

import builtins
import json
import os
import stat
import sys
from pathlib import Path

import pytest

from delfin import scheduler_client
from delfin.mlp_tools import guided_setup as gs
from delfin.mlp_tools import setup_cli

SCONTROL = "\n".join([
    "PartitionName=cpu State=UP MaxTime=3-00:00:00 TRES=cpu=4000,mem=16T,node=50",
    "PartitionName=gpu_a100 State=UP MaxTime=2-00:00:00 "
    "TRES=cpu=512,mem=4T,node=8,billing=512,gres/gpu:a100=32",
    "PartitionName=gpu_mi300 State=UP MaxTime=2-00:00:00 TRES=cpu=96,gres/gpu=8",
    "PartitionName=dev_gpu State=UP MaxTime=00:30:00 TRES=cpu=8,gres/gpu=1",
    "PartitionName=gpu_down State=DOWN TRES=cpu=8,gres/gpu=4",
])
SINFO = "cpu|(null)\ngpu_a100|gpu:a100:4\ngpu_mi300|gpu:mi300:4\n"


class Recorder:
    def __init__(self, answers=None):
        self.calls = []
        self.answers = answers or {}

    def __call__(self, argv, *a, **k):
        self.calls.append(list(argv))
        for key, val in self.answers.items():
            if key in " ".join(argv):
                return val
        return None


def _host(tmp_path, *, slurm=False, gpu=False, conda_envs=(), settings=None,
          scheduler=None, env=None):
    bins = {}
    if slurm:
        bins.update({c: f"/fake/bin/{c}" for c in ("sbatch", "squeue", "scontrol", "sinfo")})
    if gpu:
        bins["nvidia-smi"] = "/fake/bin/nvidia-smi"
    bins["micromamba"] = "/fake/bin/micromamba"
    runs = []

    def run(argv, timeout=20):
        runs.append(list(argv))
        if argv[0] == "nvidia-smi" and len(argv) == 1:
            return 0, "| NVIDIA-SMI 550.54  Driver Version: 550.54.15  CUDA Version: 12.4 |", ""
        if argv[0] == "nvidia-smi":
            return 0, "NVIDIA A100-SXM4-80GB\n", ""
        if argv[1:3] == ["env", "list"]:
            return 0, json.dumps({"envs": [str(e) for e in conda_envs]}), ""
        return None, "", "not runnable"

    host = gs.Host(
        which=lambda name: bins.get(name),
        run=run,
        scheduler=scheduler or Recorder({"scontrol": SCONTROL, "sinfo": SINFO}),
        modules=lambda patterns: (["orca/6.1.1"] if "orca" in patterns else ["cuda/12.4"]),
        orca_installs=lambda: [],
        settings=lambda: settings if settings is not None else {"runtime": {"slurm": {}}},
        current_interpreter=lambda: {"available": False, "progext_names": []},
        env=env if env is not None else {},
    )
    host.runs = runs
    return host


def _env_on_disk(tmp_path, name="mlipenv", backend="uma", torch="2.5.1+cu124"):
    prefix = tmp_path / "envs" / name
    site = prefix / "lib" / "python3.11" / "site-packages"
    site.mkdir(parents=True)
    (site / "orca_mlips-0.3.0.dist-info").mkdir()
    (site / f"torch-{torch}.dist-info").mkdir()
    (site / "fairchem_core-2.3.0.dist-info").mkdir()
    (prefix / "bin").mkdir()
    exe = prefix / "bin" / f"orca-mlips-{backend}"
    exe.write_text("#!/bin/sh\n")
    return prefix


# ------------------------------------------------------------------ probe


def test_probe_on_a_host_without_slurm_or_gpu_says_so(tmp_path):
    host = _host(tmp_path)
    p = gs.probe(host)
    assert p["slurm"]["available"] is False
    assert "no Slurm" in p["slurm"]["reason"]
    assert host.scheduler.calls == []          # nothing asked of a scheduler that is not there
    assert p["cuda"]["nvidia_smi"] is False
    assert not any(r[0] == "nvidia-smi" for r in host.runs)
    text = "\n".join(p["summary"])
    assert "no Slurm" in text and "No NVIDIA GPU visible" in text
    assert p["orca"]["ready"] is True        # the orca/6.1.1 module


def test_probe_lists_gpu_partitions_with_gres_through_the_scheduler_client(tmp_path):
    host = _host(tmp_path, slurm=True, gpu=True)
    p = gs.probe(host)
    names = [r["name"] for r in p["slurm"]["gpu_partitions"]]
    assert names == ["gpu_a100", "gpu_mi300"]       # no cpu, no dev_, no DOWN
    a100 = p["slurm"]["gpu_partitions"][0]
    assert a100["gpus"] == 32 and a100["gpu_type"] == "a100"
    assert a100["gres"] == ["gpu:a100:4"] and a100["cuda_capable"] is True
    assert p["slurm"]["gpu_partitions"][1]["cuda_capable"] is False
    assert [c[0] for c in host.scheduler.calls] == ["scontrol", "sinfo"]
    assert p["cuda"]["driver_cuda"] == "12.4"
    assert p["cuda"]["gpus"] == ["NVIDIA A100-SXM4-80GB"]


def test_a_throttled_scheduler_is_reported_not_raised(tmp_path):
    host = _host(tmp_path, slurm=True,
                 scheduler=lambda argv: scheduler_client.THROTTLED)
    p = gs.probe(host)
    assert p["slurm"]["gpu_partitions"] == []
    assert "throttled" in p["slurm"]["reason"]


def test_a_broken_part_becomes_a_note(tmp_path):
    host = _host(tmp_path)

    def boom():
        raise OSError("settings unreadable")

    host.settings = boom
    p = gs.probe(host)
    assert any("settings" in n for n in p["notes"])
    assert p["summary"]


def test_probe_finds_an_existing_env_and_its_progext(tmp_path):
    prefix = _env_on_disk(tmp_path)
    host = _host(tmp_path, slurm=True, conda_envs=[prefix])
    p = gs.probe(host)
    [env] = p["mlip_envs"]
    assert env["progext"] == ["orca-mlips-uma"]
    assert env["packages"]["orca_mlips"] == "0.3.0"
    assert env["torch_build"] == "cu124" and env["python"] == "3.11"
    plan = gs.propose(p, backend="uma")
    assert plan["env_exists"] and plan["env_prefix"] == str(prefix)
    assert "already exists" in plan["note"]


def test_env_overrides_are_named(tmp_path):
    host = _host(tmp_path, env={"DELFIN_MLP_CLUSTER_GRES": "gpu:2"})
    p = gs.probe(host)
    assert p["settings"]["env_overrides"] == ["DELFIN_MLP_CLUSTER_GRES"]
    assert p["settings"]["effective"]["gres"] == "gpu:2"


# ------------------------------------------------------------------- plan


@pytest.mark.parametrize("driver,gpus,tag", [
    ("12.4", [], "cu124"), ("12.2", [], "cu121"), ("11.8", [], "cu118"),
    ("12.9", [], "cu128"), ("13.0", [], "cu128"), ("", [], "cu126"),
    ("", ["NVIDIA B200"], "cu128"), ("11.4", [], "cu118"),
])
def test_the_torch_build_follows_the_driver(driver, gpus, tag):
    assert gs.cuda_build_for(driver, gpus)[0] == tag


@pytest.mark.parametrize("backend", gs.ORCA_MLIPS_BACKENDS)
def test_the_plan_commands_per_backend(tmp_path, backend):
    p = gs.probe(_host(tmp_path, slurm=True, gpu=True))
    plan = gs.propose(p, backend=backend)
    assert plan["partition"] == "gpu_a100"
    assert plan["cuda_build"] == "cu124"
    cmds = [c["argv"] for c in plan["commands"] if c["kind"] == "run"]
    assert cmds[0][:5] == ["micromamba", "create", "-y", "-n", gs.DEFAULT_ENV_NAME]
    assert "python=3.11" in cmds[0]
    assert cmds[1][-2:] == ["--index-url", "https://download.pytorch.org/whl/cu124"]
    assert cmds[2][-1] == f"orca-mlips[{backend}]"
    assert plan["remove_command"]["argv"][:4] == ["micromamba", "env", "remove", "-y"]
    assert plan["blockers"] == []


def test_licence_records_per_backend():
    uma = gs.licence_notes("uma")
    assert uma["gated"] and "hf auth login" in " ".join(uma["lines"])
    mace = gs.licence_notes("mace")
    assert mace["non_commercial"] and "non-commercial" in " ".join(mace["lines"])
    assert gs.licence_notes("orb")["weights"]["name"] == "Apache-2.0"
    assert gs.licence_notes("aimnet2")["weights"]["name"] == "MIT"
    for b in gs.ORCA_MLIPS_BACKENDS:
        assert gs.licence_notes(b)["package_code"]["name"] == "MIT"
    assert gs.default_backend()[0] == "orb"


def test_uma_adds_a_manual_login_step_that_is_never_run(tmp_path):
    plan = gs.propose(gs.probe(_host(tmp_path, slurm=True)), backend="uma")
    manual = [c for c in plan["commands"] if c["kind"] == "manual"]
    assert [c["command"] for c in manual] == ["hf auth login"]
    ran = []
    res = gs.run_install(plan["commands"], approve=lambda c, p: True,
                         runner=lambda argv: ran.append(argv) or 0)
    assert res["ok"]
    assert ["hf", "auth", "login"] not in ran
    assert len(ran) == 3


def test_no_conda_tool_is_a_blocker_without_commands(tmp_path):
    host = _host(tmp_path, slurm=True)
    host.which = (lambda name: None if name in ("micromamba", "mamba", "conda")
                  else f"/fake/bin/{name}")
    plan = gs.propose(gs.probe(host))
    assert plan["commands"] == []
    assert any("micromamba" in b for b in plan["blockers"])


def test_a_mi300_only_cluster_has_no_partition_to_propose(tmp_path):
    sched = Recorder({"scontrol": "PartitionName=gpu_mi300 State=UP TRES=gres/gpu=8",
                      "sinfo": ""})
    plan = gs.propose(gs.probe(_host(tmp_path, slurm=True, scheduler=sched)))
    assert plan["partition"] == ""
    assert any("GPU partition" in b for b in plan["blockers"])


@pytest.mark.parametrize("kwargs", [
    {"backend": "ani2x"}, {"partition": "gpu;rm -rf"}, {"env_name": "x;y"},
    {"env_name": "../up"},
])
def test_refusals_in_the_plan(tmp_path, kwargs):
    p = gs.probe(_host(tmp_path, slurm=True))
    with pytest.raises(ValueError):
        gs.propose(p, **kwargs)


def test_a_prefix_with_shell_characters_is_refused():
    with pytest.raises(ValueError):
        gs.settings_values(partition="gpu", prefix="/x/$(id)", backend="orb")


# ---------------------------------------------------------------- install


def test_no_command_runs_without_its_approval():
    cmds = [gs.Command(["micromamba", "create", "-n", "e"], "a").to_dict(),
            gs.Command(["micromamba", "run", "-n", "e", "pip", "install", "x"], "b").to_dict()]
    ran, asked = [], []

    def approve(c, p):
        asked.append(c)
        return len(asked) == 1        # yes to the first, no to the second

    res = gs.run_install(cmds, approve=approve, runner=lambda a: ran.append(a) or 0)
    assert asked == [cmds[0]["command"], cmds[1]["command"]]
    assert ran == [cmds[0]["argv"]]
    assert res["stopped"] == "declined" and not res["ok"]


def test_nobody_to_ask_runs_nothing():
    cmds = [gs.Command(["micromamba", "create", "-n", "e"], "a").to_dict()]
    ran = []
    res = gs.run_install(cmds, approve=lambda c, p: None, runner=lambda a: ran.append(a) or 0)
    assert ran == [] and res["stopped"] == "no approval channel"
    with pytest.raises(ValueError):
        gs.run_install(cmds, approve=None)


def test_a_failed_step_stops_the_rest():
    cmds = [gs.Command(["a"], "").to_dict(), gs.Command(["b"], "").to_dict()]
    ran = []
    res = gs.run_install(cmds, approve=lambda c, p: True,
                         runner=lambda a: ran.append(a) or 1)
    assert ran == [["a"]] and res["stopped"] == "exit code 1"


def _cli(argv, tmp_path, *, ask=lambda q: None, host=None, runner=None):
    lines = []
    host = host or _host(tmp_path, slurm=True, gpu=True)
    rc = setup_cli.main(argv, ask=ask, out=lines.append,
                        probe_fn=lambda: gs.probe(host), runner=runner,
                        settings_path=tmp_path / "settings.json")
    return rc, "\n".join(lines)


def test_the_cli_install_refuses_without_a_terminal(tmp_path):
    ran = []
    rc, text = _cli(["install"], tmp_path, ask=lambda q: None,
                    runner=lambda a: ran.append(a) or 0)
    assert rc == 1 and ran == []
    assert "no terminal to ask on" in text


def test_the_cli_install_asks_for_each_command(tmp_path):
    asked, ran = [], []

    def ask(q):
        asked.append(q)
        return True

    rc, _text = _cli(["install", "--backend", "aimnet2"], tmp_path, ask=ask,
                     runner=lambda a: ran.append(a) or 0)
    assert rc == 0
    assert len(asked) == len(ran) == 3
    assert all("run: micromamba" in q for q in asked)


def test_the_cli_remove_is_asked_and_refused_without_a_terminal(tmp_path):
    ran = []
    rc, text = _cli(["remove"], tmp_path, runner=lambda a: ran.append(a) or 0)
    assert rc == 1 and ran == [] and "env remove" in text
    rc, _ = _cli(["remove"], tmp_path, ask=lambda q: True,
                 runner=lambda a: ran.append(a) or 0)
    assert rc == 0 and ran == [["micromamba", "env", "remove", "-y", "-n", gs.DEFAULT_ENV_NAME]]


def test_no_yes_flag_on_install_or_remove():
    parser = setup_cli.build_parser()
    for cmd in (["install", "--yes"], ["remove", "--yes"]):
        with pytest.raises(SystemExit):
            parser.parse_args(cmd)


# -------------------------------------------------------------- configure


def test_configure_shows_the_diff_and_writes_only_with_write(tmp_path, monkeypatch):
    from delfin import user_settings

    path = tmp_path / "settings.json"
    saved = []
    real = user_settings.save_settings
    monkeypatch.setattr(user_settings, "save_settings",
                        lambda s, p=None: saved.append(p) or real(s, p))
    values = gs.settings_values(partition="gpu_a100", prefix=str(tmp_path / "env"),
                                backend="orb")
    preview = gs.configure(values, settings_path=path)
    assert not preview["written"] and saved == []
    keys = {d["key"] for d in preview["diff"]}
    assert {"runtime.slurm.gpu_partitions", "runtime.mlp.gres",
            "runtime.mlp.env_activation", "runtime.mlp.cpus_per_task",
            "runtime.mlp.time_limit", "runtime.mlp.mem_mb"} <= keys
    res = gs.configure(values, write=True, settings_path=path)
    assert res["written"] and res["not_kept"] == [] and saved == [path]
    data = json.loads(path.read_text())
    assert data["runtime"]["slurm"]["gpu_partitions"] == "gpu_a100"
    assert data["runtime"]["mlp"]["backend"] == "orb"
    # A second identical configure is no change and no write.
    again = gs.configure(values, write=True, settings_path=path)
    assert again["diff"] == [] and not again["written"]
    cfg = gs.cluster_config(user_settings.load_settings(path), env={})
    assert cfg.gpu_partitions == "gpu_a100" and cfg.env_activation.startswith("export PATH=")


@pytest.mark.parametrize("key,value", [
    ("runtime.slurm.gpu_partitions", "gpu a100"),
    ("runtime.mlp.gres", "gpu:1;id"),
    ("runtime.mlp.time_limit", "forever"),
    ("runtime.mlp.mem_mb", -1),
    ("runtime.mlp.env_activation", "a\nb"),
])
def test_configure_refuses_bad_values(tmp_path, key, value):
    with pytest.raises(ValueError):
        gs.configure({key: value}, settings_path=tmp_path / "s.json")


def test_the_cli_configure_writes_only_after_a_yes(tmp_path):
    path = tmp_path / "settings.json"
    rc, text = _cli(["configure", "--write"], tmp_path, ask=lambda q: None)
    assert rc == 1 and "pass --yes" in text and not path.exists()
    rc, text = _cli(["configure", "--write"], tmp_path, ask=lambda q: True)
    assert rc == 0 and "runtime.mlp.gres" in text
    assert json.loads(path.read_text())["runtime"]["slurm"]["gpu_partitions"] == "gpu_a100"


def test_runtime_mlp_survives_save_and_a_dashboard_payload(tmp_path):
    from delfin import user_settings as us

    path = tmp_path / "s.json"
    us.save_settings({"runtime": {"mlp": {"gres": "gpu:1", "mem_mb": "8000", "x": 1}}}, path)
    assert us.load_settings(path)["runtime"]["mlp"] == {"gres": "gpu:1", "mem_mb": 8000}
    merged = us.carry_runtime_sections({"mlp": {"gres": "gpu:1"}}, {"backend": "auto"})
    assert merged == {"backend": "auto", "mlp": {"gres": "gpu:1"}}
    # A file without runtime.mlp gets none.
    us.save_settings({"runtime": {}}, path)
    assert "mlp" not in us.load_settings(path)["runtime"]


# ----------------------------------------------------------------- verify


def _cfg(partition="gpu_a100"):
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    return MLPClusterConfig(gpu_partitions=partition, gres="gpu:1",
                            env_activation='export PATH="/envs/m/bin:$PATH"')


def test_the_verify_job_is_rendered_for_the_partition(tmp_path):
    prep = gs.prepare_verify_job(_cfg(), backend="orb", job_dir=tmp_path / "v")
    script = Path(prep["script"]).read_text()
    for needle in ("--partition=gpu_a100", "--gres=gpu:1", "--nodes=1",
                   "--job-name=mlp_verify", "--time=00:20:00",
                   'export PATH="/envs/m/bin:$PATH"', "mlip_verify.py --progext orca-mlips-orb"):
        assert needle in script, needle
    assert "--device cuda" in prep["ext_params"]
    assert Path(prep["input"]).is_file() and "orca-mlips-orb" in Path(prep["input"]).read_text()
    assert (tmp_path / "v" / "mlip_verify.py").is_file()
    assert prep["submitted"] is False and prep["submit_command"].startswith("sbatch ")


def test_no_verify_job_without_a_gpu_partition(tmp_path):
    with pytest.raises(ValueError):
        gs.prepare_verify_job(_cfg(""), backend="orb", job_dir=tmp_path / "v")
    assert not (tmp_path / "v").exists()


def test_the_verify_submission_goes_through_the_scheduler_client(tmp_path):
    scheduler_client.reset()
    prep = gs.prepare_verify_job(_cfg(), backend="orb", job_dir=tmp_path / "v")
    seen = []
    res = gs.submit_verify_job(prep, active_jobs=0,
                               runner=lambda argv: seen.append(argv) or "Submitted batch job 4242\n")
    assert seen == [["sbatch", prep["script"]]]
    assert res["submitted"] and res["job_id"] == "4242"
    scheduler_client.reset()


def test_the_job_cap_holds_for_the_verify_job(tmp_path):
    from delfin.mlp_tools.cluster_dispatch import JobCapExceeded

    scheduler_client.reset()
    prep = gs.prepare_verify_job(_cfg(), backend="orb", job_dir=tmp_path / "v")
    seen = []
    with pytest.raises(JobCapExceeded):
        gs.submit_verify_job(prep, active_jobs=4, runner=lambda a: seen.append(a) or "x")
    assert seen == []


def test_a_refused_submission_names_the_partition(tmp_path):
    scheduler_client.reset()
    prep = gs.prepare_verify_job(_cfg(), backend="orb", job_dir=tmp_path / "v")
    res = gs.submit_verify_job(
        prep, active_jobs=0,
        runner=lambda a: "sbatch: error: invalid partition specified: gpu_a100")
    assert not res["submitted"]
    assert [c["cause"] for c in res["likely_causes"]] == ["partition"]
    scheduler_client.reset()


def test_the_cli_verify_is_a_dry_run_unless_submitted_with_a_yes(tmp_path, monkeypatch):
    monkeypatch.setattr(gs, "cluster_config", lambda *a, **k: _cfg())
    calls = []
    monkeypatch.setattr(gs, "submit_verify_job",
                        lambda prep, **k: calls.append(prep) or {"submitted": True, "job_id": "1",
                                                                  "job_dir": prep["job_dir"]})
    rc, text = _cli(["verify", "--backend", "orb", "--job-dir", str(tmp_path / "v")], tmp_path)
    assert rc == 0 and calls == [] and "dry run" in text
    rc, text = _cli(["verify", "--backend", "orb", "--job-dir", str(tmp_path / "v"),
                     "--submit"], tmp_path, ask=lambda q: None)
    assert rc == 1 and calls == []
    rc, text = _cli(["verify", "--backend", "orb", "--job-dir", str(tmp_path / "v"),
                     "--submit"], tmp_path, ask=lambda q: True)
    assert rc == 0 and len(calls) == 1


def _executable(path: Path, body: str) -> Path:
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_the_verify_script_measures_cold_and_warm_calls_and_the_opt(tmp_path):
    """The stdlib verify script, run against a fake ProgExt and a fake ORCA."""
    job = tmp_path / "v"
    prep = gs.prepare_verify_job(_cfg(), backend="orb", job_dir=job)
    progext = _executable(tmp_path / "fake-progext", f"""#!{sys.executable}
import sys, time
if "--stop-server" in sys.argv:
    sys.exit(0)
inp = sys.argv[1]
time.sleep(0.05)
with open(inp.replace(".extinp.tmp", ".engrad"), "w") as fh:
    fh.write("#\\n3\\n#\\n-76.4\\n#\\n" + "\\n".join(["0.001"] * 9) + "\\n")
""")
    orca = _executable(tmp_path / "fake-orca", f"""#!{sys.executable}
print("GEOMETRY OPTIMIZATION CYCLE 1")
print("GEOMETRY OPTIMIZATION CYCLE 2")
print("FINAL SINGLE POINT ENERGY      -76.40001")
print("THE OPTIMIZATION HAS CONVERGED")
print("****ORCA TERMINATED NORMALLY****")
""")
    import subprocess
    env = dict(os.environ, PYTHONNOUSERSITE="1")
    done = subprocess.run(
        [sys.executable, "mlip_verify.py", "--progext", str(progext),
         "--ext-params", prep["ext_params"], "--orca", str(orca)],
        cwd=job, capture_output=True, text=True, timeout=300, env=env)
    assert done.returncode == 0, done.stderr
    raw = json.loads((job / gs.VERIFY_RESULT).read_text())
    assert len(raw["gradient_calls"]) == 3
    res = gs.read_verify_result(job)
    n = res["numbers"]
    assert n["energy_eh"] == pytest.approx(-76.4)
    assert n["opt_cycles"] == 2 and n["opt_final_energy_eh"] == pytest.approx(-76.40001)
    assert n["cold_call_s"] > 0 and n["warm_call_s"] > 0
    assert n["start_and_model_load_s"] == pytest.approx(n["cold_call_s"] - n["warm_call_s"])
    # The CI interpreter has no CUDA device: the one failing check is named.
    if not raw["torch"].get("cuda_available"):
        assert res["status"] == "fail" and res["detail"] == "torch sees no CUDA device"
    assert "ExtOpt Opt" in gs.verify_summary(res)


def _result(job: Path, **over) -> None:
    job.mkdir(parents=True, exist_ok=True)
    res = {"torch": {"cuda_available": True, "import_s": 2.0, "wall_s": 2.5,
                     "device": "NVIDIA A100", "torch": "2.5.1", "cuda_build": "12.4"},
           "gradient_calls": [{"rc": 0, "energy_eh": -76.4, "n_gradient": 9, "wall_s": 20.0},
                              {"rc": 0, "energy_eh": -76.4, "n_gradient": 9, "wall_s": 0.3},
                              {"rc": 0, "energy_eh": -76.4, "n_gradient": 9, "wall_s": 0.2}],
           "timing": {"cold_call_s": 20.0, "warm_call_s": 0.25, "start_and_load_s": 19.75},
           "extopt_opt": {"rc": 0, "wall_s": 12.0, "cycles": 6, "converged": True,
                          "terminated_normally": True, "final_energy_eh": -76.41}}
    res.update(over)
    (job / gs.VERIFY_RESULT).write_text(json.dumps(res))


def test_a_passing_result_reports_the_numbers(tmp_path):
    _result(tmp_path / "v")
    res = gs.read_verify_result(tmp_path / "v")
    assert res["status"] == "pass" and res["likely_causes"] == []
    assert res["numbers"]["start_and_model_load_s"] == 19.75
    assert res["numbers"]["opt_s_per_cycle"] == 2.0
    assert "19.75 s" in gs.verify_summary(res)


@pytest.mark.parametrize("stderr,cause", [
    ("RuntimeError: CUDA error: no kernel image is available for execution", "cuda"),
    ("huggingface_hub.errors.GatedRepoError: 401 Client Error", "login"),
    ("/var/spool/slurmd/job/script: line 9: orca: command not found", "module"),
    ("orca-mlips-orb: command not found", "env"),
])
def test_a_failing_result_names_the_likely_cause(tmp_path, stderr, cause):
    _result(tmp_path / "v", gradient_calls=[{"rc": 1, "wall_s": 3.0, "stderr_tail": stderr}])
    res = gs.read_verify_result(tmp_path / "v")
    assert res["status"] == "fail"
    assert cause in [c["cause"] for c in res["likely_causes"]]


def test_no_result_yet_is_pending(tmp_path):
    (tmp_path / "v").mkdir()
    assert gs.read_verify_result(tmp_path / "v")["status"] == "pending"
    (tmp_path / "v" / "mlp_verify_1.err").write_text("sbatch: error: Invalid generic resource (gres) specification")
    res = gs.read_verify_result(tmp_path / "v")
    assert res["status"] == "fail" and res["likely_causes"][0]["cause"] == "partition"


# ----------------------------------------------------------------- tokens


def test_the_uma_path_never_touches_a_token(tmp_path, monkeypatch):
    secret = "hf_" + "Q" * 30
    for var in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN"):
        monkeypatch.setenv(var, secret)
    touched = []
    real_open = builtins.open

    def spy_open(file, *a, **k):
        touched.append(str(file))
        return real_open(file, *a, **k)

    monkeypatch.setattr(builtins, "open", spy_open)
    real_stat = os.stat

    def spy_stat(path, *a, **k):
        touched.append(str(path))
        return real_stat(path, *a, **k)

    monkeypatch.setattr(os, "stat", spy_stat)
    host = _host(tmp_path, slurm=True, gpu=True, env=dict(os.environ))
    lines = []
    setup_cli.main(["plan", "--backend", "uma", "--json"], ask=lambda q: None,
                   out=lines.append, probe_fn=lambda: gs.probe(host),
                   settings_path=tmp_path / "s.json")
    setup_cli.main(["plan", "--backend", "uma"], ask=lambda q: None,
                   out=lines.append, probe_fn=lambda: gs.probe(host),
                   settings_path=tmp_path / "s.json")
    text = "\n".join(lines)
    assert secret not in text
    assert "hf auth login" in text
    assert not [t for t in touched if "huggingface" in t.lower() or "token" in t.lower()]


def test_the_setup_code_names_no_token_source():
    root = Path(gs.__file__).resolve().parent
    for name in ("guided_setup.py", "setup_cli.py"):
        src = (root / name).read_text()
        for needle in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN", ".cache/huggingface",
                       "huggingface/token", "get_token"):
            assert needle not in src, (name, needle)


# ----------------------------------------------------------------- doctor


def test_the_doctor_line(tmp_path):
    assert gs.mlip_readiness({"runtime": {"slurm": {}}})["detail"].startswith(
        "MLIP: not set up")
    prefix = _env_on_disk(tmp_path, backend="orb")
    ready = gs.mlip_readiness({"runtime": {
        "slurm": {"gpu_partitions": "gpu_a100"},
        "mlp": {"backend": "orb", "env_prefix": str(prefix)}}})
    assert ready["ready"] and ready["detail"] == "MLIP: env/partition/backend ready"
    missing = gs.mlip_readiness({"runtime": {
        "slurm": {"gpu_partitions": ""},
        "mlp": {"backend": "uma", "env_prefix": str(prefix)}}})
    assert not missing["ready"]
    assert missing["detail"].startswith("MLIP: missing env (no orca-mlips-uma")
    assert "partition" in missing["detail"]
    from_activation = gs.mlip_readiness({"runtime": {
        "slurm": {"gpu_partitions": "g"},
        "mlp": {"backend": "orb", "env_activation": f'export PATH="{prefix}/bin:$PATH"'}}})
    assert from_activation["ready"]


def test_the_agent_doctor_reports_the_mlip_line(tmp_path, monkeypatch):
    from delfin.agent import doctor

    assert ("MLIP", "_check_mlip") in doctor._CHECK_ATTRS
    rows = doctor._check_mlip({"settings": {"runtime": {
        "slurm": {"gpu_partitions": "g"}, "mlp": {"backend": "orb"}}}})
    assert rows[0]["status"] == doctor.WARN and rows[0]["detail"] == "MLIP: missing env"
    rows = doctor._check_mlip({"settings": {"runtime": {"slurm": {}}}})
    assert rows[0]["status"] == doctor.PASS
