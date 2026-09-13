"""One installer for every Linux login node, and nobody's paths inside it.

The installer used to be bwUniCluster's: it loaded that site's modules, wrote
that site's profile and stopped wherever there was no ``module`` command. The
tools it installs are the same everywhere; what differs between machines is
found on the machine. And what it installs is what DELFIN installs on demand
and what Settings offers -- one list, not three that drift apart.
"""

import os
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
INSTALLER = REPO / "delfin" / "installers" / "install_delfin.sh"


def _dry_run(tmp_path, *args, extra_path=None, extra_env=None):
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    path = [str(pathlib.Path(sys.executable).parent), "/usr/bin", "/bin"]
    if extra_path is not None:
        path.insert(0, str(extra_path))
    done = subprocess.run(
        ["bash", str(INSTALLER), "--dry-run", *args],
        capture_output=True, text=True, timeout=120,
        env={
            "HOME": str(home),
            "PATH": os.pathsep.join(path),
            "DELFIN_REPO": str(REPO),
            "DELFIN_PYTHON": sys.executable,
            **(extra_env or {}),
        },
    )
    return done, home


def _catalog():
    done = subprocess.run(["bash", str(INSTALLER), "--list"], capture_output=True, text=True, timeout=30)
    assert done.returncode == 0, done.stderr
    return {
        group.strip(): tools.split()
        for group, tools in (line.split(":", 1) for line in done.stdout.splitlines() if ":" in line)
    }


def _executable(path: pathlib.Path) -> None:
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(0o755)


def test_the_installer_names_no_site_and_no_person():
    text = INSTALLER.read_text(encoding="utf-8")
    for marker in ("/opt/bwhpc", "bwunicluster", "bwUni", "uc3", "scc.kit.edu",
                   "devel/python", "compiler/gnu", "module load", "/pfs/", "/home/"):
        assert marker not in text, marker


def test_a_dry_run_plans_everything_and_touches_nothing(tmp_path):
    done, home = _dry_run(tmp_path)

    assert done.returncode == 0, done.stdout + done.stderr
    assert "would install (qm): xtb gxtb mopac crest dftb+ xtb4stda std2" in done.stdout
    assert "would install (analysis): censo anmr cclib morfeus nglview packmol" in done.stdout
    assert "Ketcher" in done.stdout
    assert "dry run: nothing was changed" in done.stdout
    assert list(home.iterdir()) == [], "a dry run wrote into HOME"


def test_only_installs_what_is_named_and_refuses_what_it_does_not_know(tmp_path):
    done, _ = _dry_run(tmp_path, "--only", "crest,g-xtb,mace,ketcher")

    assert done.returncode == 0, done.stderr
    assert "would install (qm): gxtb crest" in done.stdout
    assert "would install (mlp): mace" in done.stdout
    assert "Ketcher" in done.stdout
    assert "(analysis)" not in done.stdout

    refused, _ = _dry_run(tmp_path, "--only", "orca")
    assert refused.returncode != 0
    assert "unknown tool" in refused.stderr


def test_the_all_profile_does_not_fetch_a_licensed_program(tmp_path):
    everything, _ = _dry_run(tmp_path, "--all", "--no-orca")
    assert everything.returncode == 0, everything.stderr
    assert "would install (mlp):" in everything.stdout
    assert "would install (ai):" in everything.stdout
    assert "multiwfn" not in everything.stdout.lower()

    named, _ = _dry_run(tmp_path, "--only", "multiwfn", "--no-orca")
    assert "would install (analysis): multiwfn" in named.stdout


def test_a_program_merely_called_orca_is_not_orca(tmp_path):
    """/usr/bin/orca is a screen reader on many desktops."""
    fake = tmp_path / "bin"
    fake.mkdir()
    _executable(fake / "orca")

    done, _ = _dry_run(tmp_path, "--profile", "core", extra_path=fake)

    assert "ORCA     not found" in done.stdout


def test_orca_is_found_where_it_was_unpacked_and_names_its_openmpi(tmp_path):
    prefix = tmp_path / "software"
    orca = prefix / "orca_6_1_1_linux_x86-64_shared_openmpi418_avx2"
    orca.mkdir(parents=True)
    for name in ("orca", "orca_plot"):
        _executable(orca / name)

    done, _ = _dry_run(tmp_path, "--prefix", str(prefix), "--profile", "core")

    assert done.returncode == 0, done.stderr
    assert f"ORCA     {orca}" in done.stdout
    assert "4.1.8" in done.stdout, "the OpenMPI version was not read from the ORCA name"


def test_openmpi_is_set_up_without_an_orca_as_well(tmp_path):
    """Genarris needs an mpicc, and an ORCA unpacked later needs its OpenMPI."""
    done, _ = _dry_run(tmp_path, "--profile", "core", "--no-orca")
    assert done.returncode == 0, done.stderr
    assert "would build OpenMPI 4.1.8" in done.stdout

    skipped, _ = _dry_run(tmp_path, "--profile", "core", "--no-orca", "--no-openmpi")
    assert "OpenMPI: not set up (--no-openmpi)" in skipped.stdout


def test_openmpi_is_built_the_way_orca_needs_it():
    """ORCA's shared build links against this OpenMPI, so it is configured as
    it always was; a changed flag is a different library under ORCA."""
    text = INSTALLER.read_text(encoding="utf-8")
    assert '--enable-mpi-cxx --enable-mca-no-build=fs-gpfs --disable-oshmem' in text
    assert 'DEFAULT_OPENMPI="${DELFIN_OPENMPI_VERSION:-4.1.8}"' in text


def _fake_openmpi(root, ompi_info_exit):
    (root / "bin").mkdir(parents=True)
    (root / "bin" / "mpirun").write_text('#!/bin/sh\necho "mpirun (Open MPI) 4.1.8"\n')
    (root / "bin" / "ompi_info").write_text(f"#!/bin/sh\nexit {ompi_info_exit}\n")
    for name in ("mpirun", "ompi_info"):
        (root / "bin" / name).chmod(0o755)
    return root


def test_an_openmpi_that_does_not_start_is_built_again(tmp_path):
    """mpirun --version answers even when a library it needs is gone."""
    broken = _fake_openmpi(tmp_path / "broken-openmpi", ompi_info_exit=1)
    done, _ = _dry_run(tmp_path, "--profile", "core", "--no-orca",
                       extra_env={"OPENMPI_PREFIX": str(broken)})
    assert "does not start" in done.stderr
    assert "would build OpenMPI 4.1.8" in done.stdout

    working = _fake_openmpi(tmp_path / "working-openmpi", ompi_info_exit=0)
    kept, _ = _dry_run(tmp_path, "--profile", "core", "--no-orca",
                       extra_env={"OPENMPI_PREFIX": str(working)})
    assert f"OpenMPI 4.1.8 at {working}" in kept.stdout


def test_a_program_of_another_family_is_installed_when_it_is_needed(monkeypatch):
    """Packmol, Genarris and CENSO come through the one installer, not the QM one."""
    from delfin import installer, qm_health
    from delfin.dashboard import gfn_optimize

    asked = []

    def check(name, depth="answer", **kw):
        there = bool(asked)
        return qm_health.ToolHealth(name=name, label=name, present=there, healthy=there,
                                    level="ok" if there else "absent")

    def not_the_qm_installer(**kw):
        raise AssertionError("Packmol went to the QM installer")

    monkeypatch.setattr(qm_health, "check_tool", check)
    monkeypatch.setattr(qm_health, "_TRIED", set())
    monkeypatch.setattr(gfn_optimize, "auto_install_allowed", lambda: True)
    monkeypatch.setattr(gfn_optimize, "install_xtb", not_the_qm_installer)
    monkeypatch.setattr(installer, "install", lambda requested, **kw: (
        asked.append(list(requested)) or {"ok": True, "results": []}))

    licensed = qm_health.ensure_tool("orca")
    assert licensed["ok"] is False and asked == [], "a licensed program is never fetched"

    answer = qm_health.ensure_tool("packmol")
    assert asked == [["packmol"]]
    assert answer["ok"] and answer["installed"]


def test_packmol_is_installed_before_a_run_gives_up(monkeypatch):
    import shutil

    from delfin import qm_health
    from delfin.analysis_tools import packmol_wrapper

    asked = []
    monkeypatch.setattr(shutil, "which", lambda name: None)
    monkeypatch.setattr(qm_health, "provide", lambda name, **kw: asked.append(name) or {"ok": False})

    assert packmol_wrapper._find_packmol() is None
    assert asked == ["packmol"]


def test_genarris_is_installed_before_a_run_gives_up(monkeypatch, tmp_path):
    import shutil

    import pytest

    import delfin.csp_tools as csp_tools
    from delfin import qm_health
    from delfin.csp_tools import genarris_wrapper

    run = next(obj for obj in vars(genarris_wrapper).values()
               if callable(obj) and (obj.__doc__ or "").startswith("Run Genarris via its CLI"))
    asked = []
    monkeypatch.setattr(shutil, "which", lambda name: None)
    monkeypatch.setattr(csp_tools, "get_csp_tools_root", lambda: tmp_path)
    monkeypatch.setattr(qm_health, "provide", lambda name, **kw: asked.append(name) or {"ok": False})

    with pytest.raises(FileNotFoundError, match="could not be installed"):
        run(tmp_path / "genarris.ini")
    assert asked == ["gnrs"]


def test_a_delfin_in_the_working_directory_is_not_the_one_installed(tmp_path):
    """python -m puts the working directory first on the import path.

    Started from inside another DELFIN checkout, that checkout's delfin was
    imported: its tool directory was staged and its links into somebody's home
    were copied into the new install.
    """
    elsewhere = tmp_path / "elsewhere"
    (elsewhere / "delfin").mkdir(parents=True)
    (elsewhere / "delfin" / "__init__.py").write_text("")
    (elsewhere / "delfin" / "installer.py").write_text('raise SystemExit("the wrong delfin")\n')

    done = subprocess.run(["bash", str(INSTALLER), "--list"], cwd=elsewhere,
                          capture_output=True, text=True, timeout=60)

    assert done.returncode == 0, done.stderr
    assert "qm: xtb" in done.stdout
    assert "the wrong delfin" not in done.stderr


def _shell_function(script: pathlib.Path, name: str) -> str:
    text = script.read_text(encoding="utf-8")
    return name + "() {" + text.split(name + "() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"


def test_a_version_that_cannot_be_read_does_not_end_the_qm_install(tmp_path):
    """crest 3.0.2 prints "crest 3.0.2"; under pipefail the grep ended the run."""
    script = REPO / "delfin" / "qm_tools" / "install_qm_tools.sh"
    (tmp_path / "crest").write_text("#!/bin/sh\necho ' crest 3.0.2'\n")
    (tmp_path / "dftb+").write_text("#!/bin/sh\necho 'no number here'\nexit 3\n")
    for name in ("crest", "dftb+"):
        (tmp_path / name).chmod(0o755)

    body = ("set -euo pipefail\n" + _shell_function(script, "version_of")
            + 'now="$(version_of crest)"; echo "crest=$now"\n'
            + 'now="$(version_of dftb+)"; echo "dftb=$now"\n'
            + 'now="$(version_of xtb)"; echo "xtb=$now"\n'
            + 'echo "still running"\n')
    done = subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=30,
                          env={"PATH": os.environ["PATH"], "BIN_DIR": str(tmp_path)})

    assert done.returncode == 0, done.stderr
    assert "crest=3.0.2" in done.stdout
    assert "dftb=present" in done.stdout
    assert "xtb=absent" in done.stdout
    assert "still running" in done.stdout


def test_analysis_commands_go_into_the_environment_not_beside_its_interpreter():
    """A venv made from /usr/bin/python3.11 resolved to /usr/bin: Permission denied."""
    import sysconfig

    script = REPO / "delfin" / "analysis_tools" / "install_analysis_tools.sh"
    body = _shell_function(script, "python_bin_dir") + f'python_bin_dir "{sys.executable}"\n'
    done = subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=30)

    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == sysconfig.get_path("scripts")


def test_packmol_gets_an_environment_of_its_own_and_a_link_beside_delfin(tmp_path):
    """Installed without a prefix it went into a base environment on nobody's PATH."""
    venv = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", "--without-pip", str(venv)], check=True)
    fake = tmp_path / "fakebin"
    fake.mkdir()
    (fake / "micromamba").write_text(
        '#!/bin/sh\n'
        'prefix=""; while [ $# -gt 0 ]; do [ "$1" = "-p" ] && prefix="$2"; shift; done\n'
        'mkdir -p "$prefix/bin" && printf "#!/bin/sh\\necho packmol\\n" > "$prefix/bin/packmol" '
        '&& chmod 755 "$prefix/bin/packmol"\n')
    (fake / "micromamba").chmod(0o755)
    root = tmp_path / "analysis_tools"
    root.mkdir()
    script = REPO / "delfin" / "analysis_tools" / "install_analysis_tools.sh"
    switches = {f"INSTALL_{name}": "0" for name in ("ANMR", "CCLIB", "CENSO", "MORFEUS", "MULTIWFN", "NGLVIEW")}

    done = subprocess.run(
        ["bash", str(script)], capture_output=True, text=True, timeout=120,
        env={"PATH": os.pathsep.join([str(fake), "/usr/bin", "/bin"]), "HOME": str(tmp_path),
             "DELFIN_PYTHON": str(venv / "bin" / "python"), "DELFIN_ANALYSIS_TOOLS_ROOT": str(root),
             "INSTALL_PACKMOL": "1", **switches})

    link = venv / "bin" / "packmol"
    assert link.is_symlink(), done.stdout + done.stderr
    assert os.readlink(link) == str(root / ".mamba_env" / "packmol" / "bin" / "packmol")


def test_every_tool_offered_is_one_its_own_installer_knows():
    catalog = _catalog()

    qm = (REPO / "delfin" / "qm_tools" / "install_qm_tools.sh").read_text(encoding="utf-8")
    cases = qm.split("install_one() {")[1].split("esac")[0]
    for tool in catalog["qm"]:
        assert re.search(rf"^\s*(?:[\w+-]+\|)*{re.escape(tool)}(?:\|[\w+-]+)*\)", cases, re.MULTILINE), tool

    for group, script in (("analysis", "analysis_tools/install_analysis_tools.sh"),
                          ("mlp", "mlp_tools/install_mlp_tools.sh"),
                          ("ai", "ai_tools/install_ai_tools.sh")):
        text = (REPO / "delfin" / script).read_text(encoding="utf-8")
        for tool in catalog[group]:
            switch = f"INSTALL_{tool.upper()}"
            assert f'{switch}="${{{switch}:-' in text, (group, tool)


def test_what_is_installed_on_demand_is_what_the_installer_offers():
    from delfin import installer
    from delfin.qm_health import INSTALLABLE, PACKAGES
    from delfin.tools._environment import tool_info

    catalog = _catalog()
    for name in INSTALLABLE:
        assert installer.find(name).name in catalog["qm"], name
    for name in catalog["qm"]:
        assert tool_info(name, "binary").policy == "auto", name
    for module, spec in PACKAGES.items():
        assert installer.find(spec["tool"]).switch == spec["flag"], module


def test_every_name_and_alias_stands_for_exactly_one_tool():
    from delfin import installer

    seen = {}
    for tool in installer.TOOLS:
        for name in (tool.name, *tool.aliases):
            assert name.lower() not in seen, (name, seen.get(name.lower()))
            seen[name.lower()] = tool.name
    assert installer.find("G-XTB").name == "gxtb"
    assert installer.find("orca") is None, "a licensed program is not DELFIN's to install"


def _fake_family_runs(monkeypatch, tmp_path):
    """Stage nothing and run nothing: record what would have been run."""
    from delfin import installer

    calls = []
    monkeypatch.setattr(installer, "_stage", lambda group: tmp_path / group)
    monkeypatch.setattr(installer, "present", lambda tool: True)

    def run(command, *, cwd, env, on_line, timeout):
        calls.append((command, env))
        return True, []

    monkeypatch.setattr(installer, "_run", run)
    return calls


def test_one_request_runs_each_family_installer_once_with_only_its_tools(monkeypatch, tmp_path):
    from delfin import installer

    calls = _fake_family_runs(monkeypatch, tmp_path)

    outcome = installer.install(["crest", "g-xtb", "cclib", "torchani", "mace", "crest"])

    assert outcome["ok"]
    assert [pathlib.Path(command[1]).name for command, _ in calls] == [
        "install_qm_tools.sh", "install_analysis_tools.sh", "install_mlp_tools.sh"]
    qm_command, qm_env = calls[0]
    assert qm_command[2:] == ["gxtb", "crest"]
    assert qm_env["DELFIN_QM_TOOLS_ROOT"] == str(tmp_path / "qm")
    analysis_env, mlp_env = calls[1][1], calls[2][1]
    assert {s for s in installer.switches("analysis") if analysis_env[s] == "1"} == {"INSTALL_CCLIB"}
    assert {s for s in installer.switches("mlp") if mlp_env[s] == "1"} == {"INSTALL_ANI2X", "INSTALL_MACE"}
    assert "FORCE_REINSTALL" not in mlp_env or mlp_env["FORCE_REINSTALL"] != "1"


def test_a_tool_still_missing_after_its_installer_is_reported_missing(monkeypatch, tmp_path):
    """These installers say "Packmol installation requires conda" and exit 0."""
    import sysconfig

    from delfin import installer

    calls = _fake_family_runs(monkeypatch, tmp_path)
    monkeypatch.setattr(installer, "present", lambda tool: tool.name != "std2")

    outcome = installer.install(["xtb", "std2"])

    assert outcome["ok"] is False
    assert outcome["results"][0]["missing"] == ["std2"]
    env = calls[0][1]
    assert env["INSTALL_STD2_FROM_SOURCE"] == "1", "std2 has no binary release"
    assert env["PATH"].split(os.pathsep)[0] == sysconfig.get_path("scripts"), (
        "an installer looked on the PATH for what it had just put into the venv")


def test_an_update_fetches_again_what_is_installed_and_nothing_else(monkeypatch, tmp_path):
    from delfin import installer

    calls = _fake_family_runs(monkeypatch, tmp_path)
    monkeypatch.setattr(installer, "present", lambda tool: tool.name in {"xtb", "cclib"})

    outcome = installer.update()

    assert outcome["ok"]
    assert len(calls) == 2
    qm_command, qm_env = calls[0]
    assert qm_command[2:] == ["xtb"]
    assert qm_env["FORCE_CONDA_UPDATE"] == "1" and qm_env["FORCE_REDOWNLOAD"] == "1"
    assert calls[1][1]["FORCE_REINSTALL"] == "1"
    assert calls[1][1]["INSTALL_CCLIB"] == "1" and calls[1][1]["INSTALL_CENSO"] == "0"


def test_repair_fixes_what_is_broken_and_leaves_what_works_or_is_absent(monkeypatch):
    from delfin import installer, qm_health

    healths = {
        "xtb": qm_health.ToolHealth(name="xtb", label="xtb", present=True, healthy=True, level="ok"),
        "crest": qm_health.ToolHealth(name="crest", label="CREST", present=False, level="fail",
                                      why="the link points into a deleted environment", repair="relink"),
        "std2": qm_health.ToolHealth(name="std2", label="std2", present=True, level="fail",
                                     why="it does not start"),
    }
    monkeypatch.setattr(qm_health, "check_tool",
                        lambda name, depth="answer", **kw: healths.get(name, qm_health.ToolHealth(name=name)))
    repaired, reinstalled = [], []
    monkeypatch.setattr(qm_health, "repair_tool", lambda name, action, **kw: (
        repaired.append((name, action)) or {"ok": True, "status": "works now", "lines": []}))
    monkeypatch.setattr(installer, "install", lambda requested, **kw: (
        reinstalled.append((list(requested), kw.get("update"))) or {"ok": True, "results": []}))
    monkeypatch.setattr(installer, "present", lambda tool: tool.name in healths)
    monkeypatch.setattr(installer, "_module_present", lambda module: False)

    outcome = installer.repair()

    assert repaired == [("crest", "relink")]
    assert reinstalled == [(["std2"], True)]
    assert {r["tools"][0]: r["action"] for r in outcome["results"]} == {
        "xtb": "checked", "crest": "repaired", "std2": "reinstalled"}

    reinstalled.clear()
    installer.repair(["gxtb"])
    assert reinstalled == [(["gxtb"], True)], "a tool named for repair and absent is installed"


def test_asking_for_one_package_installs_that_one_and_no_other(monkeypatch, tmp_path):
    from delfin import qm_health
    from delfin.dashboard import gfn_optimize

    calls = _fake_family_runs(monkeypatch, tmp_path)
    monkeypatch.setattr(qm_health, "package_present", lambda name: False)
    monkeypatch.setattr(gfn_optimize, "auto_install_allowed", lambda: True)

    for module, script, wanted in (
        ("mace", "mlp_tools/install_mlp_tools.sh", "INSTALL_MACE"),
        ("cclib", "analysis_tools/install_analysis_tools.sh", "INSTALL_CCLIB"),
    ):
        calls.clear()
        monkeypatch.setattr(qm_health, "_PACKAGES_TRIED", set())

        qm_health.ensure_package(module)

        assert len(calls) == 1, module
        env = calls[0][1]
        text = (REPO / "delfin" / script).read_text(encoding="utf-8")
        switches = set(re.findall(r'^(INSTALL_[A-Z0-9_]+)="\$\{', text, re.MULTILINE))
        assert switches <= set(env), (module, switches - set(env))
        assert sorted(s for s in switches if env[s] == "1") == [wanted]


def test_the_old_installer_names_still_run_the_one_installer():
    for relative in ("install.sh", "scripts/install_delfin_bwu.sh",
                     "scripts/verify_delfin_bwu.sh", "delfin/installers/install_delfin_bwu.sh"):
        text = (REPO / relative).read_text(encoding="utf-8")
        assert "install_delfin.sh" in text and "exec bash" in text, relative


def test_the_dashboard_runs_the_universal_installer():
    from delfin.runtime_setup import (
        get_packaged_bwunicluster_install_script,
        get_repo_bwunicluster_install_script,
    )

    assert get_packaged_bwunicluster_install_script().name == "install_delfin.sh"
    assert get_repo_bwunicluster_install_script(REPO).name == "install_delfin.sh"


def test_the_env_file_puts_the_openmpi_orca_needs_on_both_paths(tmp_path):
    from delfin.runtime_setup import write_delfin_env_file

    ompi = tmp_path / "openmpi-4.1.8"
    (ompi / "bin").mkdir(parents=True)
    (ompi / "lib").mkdir()
    orca = tmp_path / "orca_6_1_1"
    orca.mkdir()
    for name in ("orca", "orca_plot"):
        _executable(orca / name)

    env_file = write_delfin_env_file(orca_base=str(orca), openmpi_prefix=str(ompi),
                                     env_path=tmp_path / "delfin_env.sh")
    shown = subprocess.run(
        ["bash", "-c", f'export PATH=/usr/bin:/bin LD_LIBRARY_PATH=; source "{env_file}"; '
                       'echo "$PATH"; echo "$LD_LIBRARY_PATH"; echo "$OPENMPI_PREFIX"'],
        capture_output=True, text=True, check=True,
    ).stdout.splitlines()

    assert shown[0].startswith(f"{ompi / 'bin'}:")
    assert str(orca) in shown[0].split(":")
    assert shown[1] == f"{ompi / 'lib'}:{orca}", "an empty LD_LIBRARY_PATH left a stray colon"
    assert shown[2] == str(ompi)


def test_a_login_file_that_already_sources_the_env_file_is_left_alone(tmp_path, monkeypatch):
    from delfin.runtime_setup import ensure_shell_sources_delfin_env

    monkeypatch.setenv("HOME", str(tmp_path))
    bashrc = tmp_path / ".bashrc"
    bashrc.write_text("# DELFIN environment\nsource $HOME/.delfin_env.sh\n")

    ensure_shell_sources_delfin_env([bashrc], env_path=tmp_path / ".delfin_env.sh")

    assert bashrc.read_text().count(".delfin_env.sh") == 1
