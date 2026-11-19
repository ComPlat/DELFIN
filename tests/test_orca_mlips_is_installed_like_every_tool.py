"""orca-mlips: installed by DELFIN's installer, found where ORCA can start it.

orca-mlips is the program ORCA's ExtOpt starts as ``ProgExt`` for an MLIP
GOAT or optimisation.  Four live GPU jobs (docs/MLIP_GOAT_GPU.md) stopped on
D1-D4 before any of them finished; these tests hold DELFIN's side of each:

* the catalog has it, in an environment of its own built by
  install_ai_tools.sh (a CUDA PyTorch DELFIN's own pins cannot hold), and
  its presence is read from the distribution metadata of that environment;
* D3 -- ``ProgExt`` names the absolute ``orca-mlips-<backend>`` of that
  environment, not a bare name a compute node's PATH does not hold;
* D4 -- ``Ext_Params`` carries only options the plugins accept (there is no
  ``--backend``: the command is the backend);
* D2 -- the ExtOpt job script sets DELFIN's OpenMPI environment, so
  ``--ntasks=1`` with ``%pal nprocs 8`` does not die on "not enough slots";
* D1 -- a GPU partition whose MaxTime cannot hold the time limit is left out,
  and none left is a refusal naming each MaxTime, before sbatch.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from delfin import installer
from delfin.common import external_builders as eb
from delfin.mlp_tools import orca_mlips as om

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])
SCRIPT = Path(installer.__file__).parent / "ai_tools" / "install_ai_tools.sh"

#: Every option orca-mlips 1.3.0 accepts (plugins/runner_orca.py
#: run_orca_plugin + _add_server_args; plugins/cli_orca.py *_add_args).
_COMMON = {"--model", "-w", "--weights-file", "--device", "--solvent", "--solvent-model",
           "--xtb-cmd", "--xtb-acc", "--xtb-workdir", "--xtb-keep-files", "--dump-hessian",
           "--list-models", "--version", "--server-socket", "--serve", "--stop-server",
           "--no-server", "--server-idle-timeout"}
_PER_PLUGIN = {
    "uma": {"--task", "--workers", "--workers-per-node", "--max-neigh", "--radius",
            "--r-edges", "--otf-graph", "--no-otf-graph"},
    "orb": {"--precision", "--compile-model", "--loader-opt", "--calc-opt"},
    "mace": {"--dtype", "--calc-opt"},
    "aimnet2": {"--calc-opt"},
}


def _fake_env(root: Path, *, with_dist: bool = True, plugins=("aimnet2", "orb")) -> Path:
    """An orca-mlips environment as install_ai_tools.sh leaves it, without its stack."""
    prefix = root / ".mamba_env" / "orca_mlips"
    (prefix / "bin").mkdir(parents=True)
    python = prefix / "bin" / "python"
    python.write_text("#!/bin/sh\n")
    python.chmod(0o755)
    site = prefix / "lib" / "python3.12" / "site-packages"
    (site / "orca_mlips_plugins").mkdir(parents=True)
    (site / "orca_mlips_plugins" / "__init__.py").write_text('"""orca-mlips plugin package."""\n')
    if with_dist:
        meta = site / "orca_mlips-1.3.0.dist-info"
        meta.mkdir()
        meta.joinpath("METADATA").write_text(
            "Metadata-Version: 2.1\nName: orca-mlips\nVersion: 1.3.0\n")
        meta.joinpath("entry_points.txt").write_text(
            "[console_scripts]\n" + "".join(
                f"orca-mlips-{p} = orca_mlips_plugins.x:main\n" for p in plugins))
    for plugin in plugins:
        exe = prefix / "bin" / f"orca-mlips-{plugin}"
        exe.write_text("#!/bin/sh\n")
        exe.chmod(0o755)
    return prefix


@pytest.fixture()
def ai_root(tmp_path, monkeypatch):
    root = tmp_path / "ai_tools"
    root.mkdir()
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(root))
    monkeypatch.delenv("DELFIN_ORCA_MLIPS_PYTHON", raising=False)
    # Nothing on this PATH may answer for orca-mlips.
    monkeypatch.setenv("PATH", str(tmp_path / "empty-bin"))
    return root


# -- the catalog -----------------------------------------------------------------

def test_the_catalog_has_it_in_an_environment_of_its_own():
    tool = installer.find("orca-mlips")
    assert tool is installer.find("orca_mlips")
    assert tool.group == "ai" and tool.switch == "INSTALL_ORCA_MLIPS"
    assert tool.own_env == "orca_mlips" and tool.modules == ("orca_mlips_plugins",)
    assert "INSTALL_ORCA_MLIPS" in installer.switches("ai")
    # Not a package for DELFIN's own environment; a non-licensed tool of "all".
    assert "orca_mlips_plugins" not in installer.packages()
    assert "orca_mlips" in installer.profile("all") and not tool.licensed
    info = eb.PROGRAMS["orca_mlips"]
    assert "orca_mlips" not in eb.TOOLS  # not a structure builder
    assert info["managed_env"] == ".mamba_env/orca_mlips"
    assert info["install"] == "python -m delfin.installer --install orca-mlips"
    assert info["dist"] == "orca-mlips" and info["python_env"] == "DELFIN_ORCA_MLIPS_PYTHON"


def test_its_installer_builds_the_environment_with_the_two_licence_clean_plugins():
    text = SCRIPT.read_text()
    assert 'INSTALL_ORCA_MLIPS="${INSTALL_ORCA_MLIPS:-${INSTALL_ALL}}"' in text
    assert "install_own_env INSTALL_ORCA_MLIPS orca-mlips orca_mlips orca_mlips_plugins" in text
    deps = next(line for line in text.splitlines() if line.startswith("ORCA_MLIPS_PIP_DEPS="))
    # aimnet2 (MIT) and orb (Apache-2.0) only: never mace (non-commercial) or uma (gated).
    assert "aimnet==" in deps and "orb-models==" in deps and "torch==" in deps
    assert "mace" not in deps and "fairchem" not in deps
    spec = next(line for line in text.splitlines() if line.startswith("ORCA_MLIPS_SPEC="))
    assert "orca-mlips==" in spec


def test_an_own_environment_is_not_held_to_delfins_pins():
    """PIP_CONSTRAINT (DELFIN's core pins, set for every family installer) is
    for DELFIN's environment; inside an environment of its own it would
    refuse what that environment was built for."""
    text = SCRIPT.read_text()
    helper = text[text.index("epic_mace_python() {"):]
    helper = helper[:helper.index("\n}\n")]
    assert "-u PIP_CONSTRAINT" in helper


def test_it_is_present_by_its_distribution_metadata(ai_root):
    tool = installer.find("orca-mlips")
    assert not installer.present(tool)
    _fake_env(ai_root, with_dist=False)
    # The package directory alone is not an installed distribution.
    assert not installer.present(tool)


def test_it_is_present_once_the_distribution_is_there(ai_root):
    _fake_env(ai_root)
    assert installer.present(installer.find("orca-mlips"))


def test_it_is_not_a_structure_builder():
    with pytest.raises(ValueError, match="unknown external builder"):
        eb.build_frames("orca_mlips", "O")


# -- detection and ProgExt (D3) ---------------------------------------------------

def test_detect_reads_the_environment_of_its_own(ai_root):
    prefix = _fake_env(ai_root)
    info = om.detect_orca_mlips()
    assert info["available"] is True and info["version"] == "1.3.0"
    assert info["progext_names"] == ["orca-mlips-aimnet2", "orca-mlips-orb"]
    assert info["progext_paths"] == {
        "aimnet2": str(prefix / "bin" / "orca-mlips-aimnet2"),
        "orb": str(prefix / "bin" / "orca-mlips-orb"),
    }
    assert info["install_hint"].startswith("python -m delfin.installer --install orca-mlips")


def test_progext_is_the_absolute_command_of_the_environment(ai_root):
    prefix = _fake_env(ai_root)
    assert om.resolve_progext("aimnet2") == str(prefix / "bin" / "orca-mlips-aimnet2")
    assert os.path.isabs(om.resolve_progext("orb"))


def test_progext_stays_the_bare_name_when_nothing_is_installed(ai_root):
    assert om.resolve_progext("aimnet2") == "orca-mlips-aimnet2"


def test_the_goat_input_names_the_absolute_command(ai_root, tmp_path):
    """D3: ORCA starts ProgExt itself; a bare name was "not found" on a GPU node."""
    prefix = _fake_env(ai_root)
    inp = om.build_orca_mlips_input(
        "goat", elements=WATER[0], coordinates=WATER[1], backend="aimnet2",
        device="cuda", ncores=4, base_name="XTB_GOAT", work_dir=tmp_path,
        xyz_name="XTB_GOAT_start.xyz",
    )
    text = inp.read_text()
    assert f'"{prefix / "bin" / "orca-mlips-aimnet2"}"' in text


def test_the_extopt_job_script_starts_the_absolute_command(ai_root, tmp_path):
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import ExtOptJobRequest, plan_extopt_job

    prefix = _fake_env(ai_root)
    req = ExtOptJobRequest(run_type="goat", elements=WATER[0], coordinates=WATER[1],
                           provider="orca_mlips", mlips_backend="aimnet2", device="cuda",
                           job_name="g", work_dir=tmp_path, goat_max_workers=4)
    plan = plan_extopt_job(req, cfg=MLPClusterConfig(gpu_partitions="gpu_x", cpus_per_task=4))
    exe = str(prefix / "bin" / "orca-mlips-aimnet2")
    assert exe in plan.input.read_text() and exe in plan.script.read_text()


# -- Ext_Params (D4) ---------------------------------------------------------------

@pytest.mark.parametrize("backend", ["uma", "orb", "mace", "aimnet2"])
def test_ext_params_use_only_options_the_plugin_accepts(backend):
    params = om.ext_params_for(
        backend, model="m1", device="cuda", solvent="water", solvent_model="alpb",
        precision="float32-high" if backend == "orb" else None,
        extra=("--dump-hessian extopt.hess",),
    )
    options = {word for word in params.split() if word.startswith("-")}
    assert "--backend" not in options
    assert options <= _COMMON | _PER_PLUGIN[backend], options - _COMMON - _PER_PLUGIN[backend]


def test_precision_is_refused_for_a_plugin_without_it():
    with pytest.raises(ValueError, match="orb plugin only"):
        om.ext_params_for("aimnet2", precision="float32")


# -- the job script (D2) and the partition time (D1) -------------------------------

def _req(tmp_path):
    from delfin.mlp_tools.extopt_submit import ExtOptJobRequest

    return ExtOptJobRequest(run_type="opt", elements=WATER[0], coordinates=WATER[1],
                            provider="orca_mlips", mlips_backend="aimnet2", device="cuda",
                            job_name="o", work_dir=tmp_path)


def test_the_job_script_sets_delfins_openmpi_environment(tmp_path):
    from delfin.dashboard.mpi_env import OMPI_ENVIRONMENT
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import render_extopt_job_script

    script = render_extopt_job_script(
        _req(tmp_path), cfg=MLPClusterConfig(gpu_partitions="gpu_x", cpus_per_task=8),
        pal=8, provider_cmd="orca-mlips-aimnet2")
    assert "#SBATCH --ntasks=1" in script
    assert "export OMPI_MCA_rmaps_base_oversubscribe=true" in script
    for key, value in OMPI_ENVIRONMENT.items():
        assert f"export {key}={value}" in script
    # Set before ORCA starts.
    assert script.index("OMPI_MCA_rmaps_base_oversubscribe") < script.index('"$ORCA_BIN"')


def _scontrol(limits, asked=None):
    """scontrol show partition -o: one line per partition, all of them at once."""
    def ask(argv):
        if asked is not None:
            asked.append(list(argv))
        return "\n".join(f"PartitionName={name} MaxTime={limit} State=UP"
                         for name, limit in limits.items())
    return ask


def test_a_partition_too_short_for_the_time_limit_is_left_out():
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import partitions_for_time

    cfg = MLPClusterConfig(gpu_partitions="dev_gpu,gpu_long", time_limit="04:00:00")
    asked = []
    ask = _scontrol({"dev_gpu": "00:30:00", "gpu_long": "2-00:00:00"}, asked)
    assert partitions_for_time(cfg, ask) == ("gpu_long", "")
    # One query for every partition: a second one in a row came back
    # throttled by the scheduler client, and that partition was kept unread.
    assert asked == [["scontrol", "show", "partition", "-o"]]
    # The time limit itself is never changed.
    assert cfg.time_limit == "04:00:00"


def test_a_partition_whose_limit_cannot_be_read_is_kept():
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import partitions_for_time

    cfg = MLPClusterConfig(gpu_partitions="unknown_gpu", time_limit="04:00:00")
    assert partitions_for_time(cfg, lambda argv: None) == ("unknown_gpu", "")


def test_no_partition_holding_the_time_limit_is_refused_before_sbatch(tmp_path, monkeypatch):
    from delfin import scheduler_client
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import submit_extopt_job

    def no_sbatch(*args, **kwargs):
        raise AssertionError("sbatch must not be reached")

    monkeypatch.setattr(scheduler_client, "query", no_sbatch)
    cfg = MLPClusterConfig(gpu_partitions="dev_gpu", time_limit="04:00:00", cpus_per_task=4)
    with pytest.raises(RuntimeError, match=r"dev_gpu \(MaxTime 00:30:00\)"):
        submit_extopt_job(_req(tmp_path), cfg=cfg, submit=True,
                          ask=_scontrol({"dev_gpu": "00:30:00"}))


def test_the_submitted_script_names_only_the_fitting_partition(tmp_path, monkeypatch):
    from delfin import scheduler_client
    from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
    from delfin.mlp_tools.extopt_submit import submit_extopt_job

    seen = []

    def fake_query(kind, argv, runner=None):
        seen.append(list(argv))
        return "Submitted batch job 1"

    monkeypatch.setattr(scheduler_client, "query", fake_query)
    cfg = MLPClusterConfig(gpu_partitions="dev_gpu,gpu_long", time_limit="04:00:00",
                           cpus_per_task=4)
    out = submit_extopt_job(_req(tmp_path), cfg=cfg, submit=True,
                            ask=_scontrol({"dev_gpu": "00:30:00", "gpu_long": "2-00:00:00"}))
    assert out["submitted"] is True and seen and seen[0][0] == "sbatch"
    assert "#SBATCH --partition=gpu_long\n" in Path(out["script"]).read_text()


# -- the checker reads the run the CONTROL path made --------------------------------

def test_check_goat_e2e_reads_the_control_paths_output_name(tmp_path):
    """global_opt_method=mlip:* (delfin/xtb_crest.py) names ORCA's output
    output_XTB_GOAT.out; the checker looked only for XTB_GOAT.out and could
    not verify a run the workflow had finished.  SYNTHETIC files, not ORCA output."""
    from delfin.mlp_tools.goat_e2e_check import check_goat_e2e

    frame = "3\n{e} converged=true\nO 0 0 {z}\nH 0.757 0.586 {z}\nH -0.757 0.586 {z}\n"
    (tmp_path / "XTB_GOAT.finalensemble.xyz").write_text(
        frame.format(e=-76.40, z=0.0) + frame.format(e=-76.39, z=0.1))
    (tmp_path / "XTB_GOAT.globalminimum.xyz").write_text(
        "3\n-76.40\nO 0 0 0\nH 0.757 0.586 0\nH -0.757 0.586 0\n")
    (tmp_path / "output_XTB_GOAT.out").write_text(
        "# SYNTHETIC -- not ORCA output\nTotal Job Time (sec): 42.0\n"
        "                             ****ORCA TERMINATED NORMALLY****\n")
    report = check_goat_e2e(tmp_path, base_name="XTB_GOAT")
    assert report["terminated_normally"] is True
    assert report["walltime_s"] == 42.0
    assert report["ok"] is True, report["findings"]


# -- the Opt / Opt+Freq block (#611), as its first live run showed it ----------------

def test_the_energy_of_an_extopt_run_is_read():
    """ORCA prints an external program's energy as "FINAL SINGLE POINT ENERGY
    (From external program)"; DELFIN's reader returned None for it."""
    from delfin.energies import electronic_energy_in_text

    text = ("# SYNTHETIC -- not ORCA output\n"
            "FINAL SINGLE POINT ENERGY (From external program)      -155.139262225400\n"
            "FINAL SINGLE POINT ENERGY (From external program)      -155.139788599700\n")
    assert electronic_energy_in_text(text) == pytest.approx(-155.1397885997)
    # The ordinary line is read as before.
    assert electronic_energy_in_text("FINAL SINGLE POINT ENERGY     -76.4012\n") == pytest.approx(-76.4012)


def test_opt_freq_writes_numerical_frequencies(tmp_path):
    """ORCA 6.1.1 ran "! Opt Freq ExtOpt" as an optimisation only (no
    frequency block, normal termination); NumFreq is the route that runs."""
    from delfin.tools.adapters.orca_mlip_opt_opi import build_mlip_opt_input

    inp = build_mlip_opt_input(tmp_path, run_type="opt_freq", elements=WATER[0],
                               coordinates=WATER[1], backend="aimnet2")
    words = set(inp.read_text().lower().replace("!", " ").split())
    assert {"opt", "numfreq", "extopt"} <= words
    assert "freq" not in words


def test_opt_freq_without_a_frequency_block_is_not_a_success(tmp_path, monkeypatch):
    from delfin import orca as orca_module
    from delfin.tools.adapters.orca_mlip_opt_opi import OrcaMlipOptOpiAdapter
    from delfin.tools._types import StepStatus
    import delfin.mlp_tools.orca_mlips as om

    monkeypatch.setattr(om, "detect_orca_mlips", lambda *a, **k: {"available": True})

    def fake_run_orca(inp, out, **kwargs):  # an optimisation that ends without frequencies
        Path(out).write_text("# SYNTHETIC -- not ORCA output\n"
                             "THE OPTIMIZATION HAS CONVERGED\n"
                             "FINAL SINGLE POINT ENERGY (From external program)   -76.40\n"
                             "****ORCA TERMINATED NORMALLY****\n")
        return True

    monkeypatch.setattr(orca_module, "run_orca", fake_run_orca)
    result = OrcaMlipOptOpiAdapter().execute(
        tmp_path, elements=WATER[0], coordinates=WATER[1], run_type="opt_freq",
        backend="aimnet2", device="cpu", base_name="w")
    assert result.status is StepStatus.FAILED
    assert "VIBRATIONAL FREQUENCIES" in (result.error or "")
    assert result.data["energy_Eh"] == pytest.approx(-76.40)


def test_check_goat_e2e_reads_the_wall_time_orca_6_prints(tmp_path):
    """ORCA 6.1.1 prints "TOTAL RUN TIME: ..." and no "Total Job Time" line."""
    from delfin.mlp_tools.goat_e2e_check import _parse_walltime

    out = tmp_path / "x.out"
    out.write_text("# SYNTHETIC -- not ORCA output\n"
                   "TOTAL RUN TIME: 0 days 0 hours 9 minutes 28 seconds 822 msec\n")
    assert _parse_walltime(out) == pytest.approx(568.822)
