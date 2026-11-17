"""The Pi installs every optional program its hardware runs, and says which work.

Input: scripts/pi/pi_install.sh (install_deps, the optional groups, the smoke
step), scripts/pi/deploy_pi.sh (pi_env), delfin/agent/optional_tools.py and
the doctor's ``optional chemistry`` rows.

The shell functions are cut out of pi_install.sh and run in bash against a
temporary ROOT, with uv, the venv's python and DELFIN's QM installer replaced
by stubs that record their arguments. What is pinned:

* the orca extra is installed with the Pi's numpy override; a set that does
  not install is retried without the extra, and the deploy goes on;
* the dependency stamp covers pyproject.toml, the overrides and the skipped
  extras, so the orca extra arriving is a reinstall;
* each optional group installs only while missing, only when switched on,
  with the venv as constraints, and a failure is a warning with its reason;
* deploy_pi.sh forwards DELFIN_PI_* and refuses a value with shell syntax;
* the smoke lines parse back, and the doctor shows one row per group.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PI = ROOT / "scripts" / "pi"
INSTALL = (PI / "pi_install.sh").read_text(encoding="utf-8")

FUNCTIONS = ("write_requirements", "uv_install_requirements", "install_deps",
             "extra_skipped", "deps_stamp", "deps_current", "group_specs",
             "group_enabled", "first_reason", "install_pip_group",
             "install_xtb_group", "group_status", "install_optional_groups",
             "smoke_optional")


def _function(text: str, name: str) -> str:
    start = text.index(f"\n{name}() {{\n") + 1
    return text[start:text.index("\n}\n", start) + 3]


def _assignment(name: str) -> str:
    m = re.search(rf"^{name}=.*$", INSTALL, re.M)
    assert m, name
    return m.group(0)


def _harness(tmp_path: Path, *, uv_fail: str = "", status: str = "missing",
             with_overrides: bool = True) -> Path:
    """A ROOT with src/pyproject.toml, stub uv/python/installer; returns ROOT."""
    root = tmp_path / "delfin"
    (root / "src" / "scripts" / "pi").mkdir(parents=True)
    (root / "src" / "pyproject.toml").write_text(
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"), encoding="utf-8")
    if with_overrides:
        (root / "src" / "scripts" / "pi" / "pi-overrides.txt").write_text(
            (PI / "pi-overrides.txt").read_text(encoding="utf-8"), encoding="utf-8")
    log = tmp_path / "calls.log"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    uv = bin_dir / "uv"
    # Records every call; fails a call whose arguments contain $UV_FAIL.
    uv.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        printf 'uv %s\\n' "$*" >> {log}
        for a in "$@"; do
            if [[ "$a" == -c ]]; then shift_next=1; continue; fi
            if [[ -n "${{shift_next:-}}" ]]; then printf 'constraints %s\\n' "$(wc -l < "$a")" >> {log}; shift_next=; fi
        done
        if [[ "$1 $2" == "pip freeze" ]]; then printf 'numpy==2.2.6\\n-e file:///x\\nrdkit==2025.9.6\\n'; exit 0; fi
        if [[ -n "{uv_fail}" && "$*" == *"{uv_fail}"* ]]; then
            echo "  x No solution found when resolving dependencies:" >&2
            echo "error: stub refuses {uv_fail}" >&2
            exit 1
        fi
        exit 0
        """), encoding="utf-8")
    py = bin_dir / "python"
    # The venv's python: DELFIN's probe answers $STATUS for every group;
    # the smoke step records its arguments; `import pymol` succeeds;
    # anything else is the real interpreter.
    py.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        if [[ "$1 $2" == "-m delfin.agent.optional_tools" ]]; then
            printf 'py %s\\n' "$*" >> {log}
            if [[ "$3" == list ]]; then
                for g in opi pyscf openmm torch xtb; do printf 'smoke %s %s x\\n' "$g" "{status}"; done
            fi
            exit 0
        fi
        if [[ "$1" == -c && "$2" == "import pymol" ]]; then exit 0; fi
        exec {sys.executable} "$@"
        """), encoding="utf-8")
    installer = root / "src" / "delfin" / "qm_tools" / "install_qm_tools.sh"
    installer.parent.mkdir(parents=True)
    installer.write_text(textwrap.dedent(f"""\
        printf 'qm %s XTB_SPEC=%s CREST_SPEC=%s ROOT=%s\\n' "$*" "$XTB_SPEC" "$CREST_SPEC" "$DELFIN_QM_TOOLS_ROOT" >> {log}
        mkdir -p "$DELFIN_QM_TOOLS_ROOT/bin" "$DELFIN_QM_TOOLS_ROOT/.mamba_env/xtb/bin"
        for b in xtb tblite dftd4; do printf '#!/bin/sh\\n' > "$DELFIN_QM_TOOLS_ROOT/.mamba_env/xtb/bin/$b"; chmod +x "$DELFIN_QM_TOOLS_ROOT/.mamba_env/xtb/bin/$b"; done
        ln -sfn "$DELFIN_QM_TOOLS_ROOT/.mamba_env/xtb/bin/xtb" "$DELFIN_QM_TOOLS_ROOT/bin/xtb"
        """), encoding="utf-8")
    for f in (uv, py):
        f.chmod(0o755)
    return root


def _run(tmp_path: Path, root: Path, body: str, env: dict | None = None,
         skip: str | None = None) -> subprocess.CompletedProcess:
    skip_line = (f'SKIP_EXTRAS="{skip}"' if skip is not None
                 else _assignment("SKIP_EXTRAS"))
    script = "\n".join([
        "set -euo pipefail",
        f'HOME="{tmp_path}"',
        f'ROOT="{root}"',
        f'VENV_PY="{tmp_path}/bin/python"',
        f'UV="{tmp_path}/bin/uv"',
        _assignment("SKIP_DEPS_RE"),
        skip_line,
        _assignment("PI_OVERRIDES"),
        _assignment("EXTRA_PIP"),
        _assignment("PYMOL_GIT"),
        _assignment("DEPS_STAMP"),
        _assignment("PI_XTB_SPEC"),
        _assignment("PI_CREST_SPEC"),
        _assignment("PI_QM_ROOT"),
        _assignment("PI_GROUPS"),
        "install_delfin_editable() { :; }",
        *(_function(INSTALL, n) for n in FUNCTIONS),
        f'cd "{root}"',
        body,
    ])
    full_env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
    full_env.update(env or {})
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                          env=full_env, timeout=60)


def _calls(tmp_path: Path) -> list[str]:
    log = tmp_path / "calls.log"
    return log.read_text().splitlines() if log.exists() else []


def _requirements(root: Path) -> list[str]:
    return (root / "requirements.pi.txt").read_text().split("\n")


# --- the orca extra -------------------------------------------------------

def test_the_orca_extra_is_installed_with_the_pi_override(tmp_path):
    root = _harness(tmp_path)
    proc = _run(tmp_path, root, "install_deps")
    assert proc.returncode == 0, proc.stderr
    (call,) = [c for c in _calls(tmp_path) if "-r " in c]
    assert f"--overrides {root}/src/scripts/pi/pi-overrides.txt" in call
    assert "--torch-backend cpu" in call
    assert "orca-pi==2.0.0" in _requirements(root)
    assert (root / ".deps-inputs.sha256").read_text().strip()


def test_orca_in_the_skip_list_installs_without_the_override(tmp_path):
    root = _harness(tmp_path)
    proc = _run(tmp_path, root, "install_deps", skip="ai-retro orca pymol")
    assert proc.returncode == 0, proc.stderr
    (call,) = [c for c in _calls(tmp_path) if "-r " in c]
    assert "--overrides" not in call
    assert not any(r.startswith("orca-pi") for r in _requirements(root))


def test_a_set_with_orca_that_fails_is_retried_without_it(tmp_path):
    root = _harness(tmp_path, uv_fail="--overrides")
    (root / ".deps-inputs.sha256").write_text("old\n")
    proc = _run(tmp_path, root, "install_deps; echo deploy-goes-on")
    assert proc.returncode == 0, proc.stderr
    assert "deploy-goes-on" in proc.stdout
    calls = [c for c in _calls(tmp_path) if "-r " in c]
    assert len(calls) == 2 and "--overrides" in calls[0] and "--overrides" not in calls[1]
    assert not any(r.startswith("orca-pi") for r in _requirements(root))
    assert "OPI stays missing" in proc.stderr
    # No stamp: the next deploy tries the orca extra again.
    assert not (root / ".deps-inputs.sha256").exists()


def test_a_tree_without_the_overrides_file_leaves_orca_out(tmp_path):
    root = _harness(tmp_path, with_overrides=False)
    proc = _run(tmp_path, root, "install_deps")
    assert proc.returncode == 0, proc.stderr
    (call,) = [c for c in _calls(tmp_path) if "-r " in c]
    assert "--overrides" not in call
    assert not any(r.startswith("orca-pi") for r in _requirements(root))
    assert "pi-overrides.txt missing" in proc.stderr


def test_the_stamp_follows_pyproject_overrides_and_the_skip_list(tmp_path):
    root = _harness(tmp_path)
    body = textwrap.dedent("""\
        deps_stamp > "$DEPS_STAMP"
        deps_current && echo same
        SKIP_EXTRAS="$SKIP_EXTRAS orca"; deps_current || echo skip-changed
        SKIP_EXTRAS="${SKIP_EXTRAS% orca}"; deps_current && echo back
        echo 'numpy>=2.2.6,<2.3 ' >> "$PI_OVERRIDES"; deps_current || echo overrides-changed
        """)
    proc = _run(tmp_path, root, body)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == ["same", "skip-changed", "back", "overrides-changed"]


def test_a_deploy_reinstalls_when_the_dependency_inputs_change():
    install = INSTALL[INSTALL.index("\ninstall)\n"):INSTALL.index("\nrollback)\n")]
    assert "deps_current" in install and "install_deps" in install
    assert "old_hash" not in INSTALL, "pyproject.toml alone never saw the orca extra arrive"
    rollback = INSTALL[INSTALL.index("\nrollback)\n"):INSTALL.index("\nstatus)\n")]
    assert "deps_current || install_deps" in rollback
    for section in (install, rollback):
        assert section.index("install_optional_groups") < section.index("drop_foreign_tests_pkg")
        assert "smoke_optional" in section


# --- optional groups ------------------------------------------------------

def test_every_missing_group_is_installed_pinned_against_the_venv(tmp_path):
    root = _harness(tmp_path)
    proc = _run(tmp_path, root, "install_optional_groups")
    assert proc.returncode == 0, proc.stderr
    calls = _calls(tmp_path)
    installs = [c for c in calls if c.startswith("uv pip install")]
    assert any("pyscf==2.14.0" in c for c in installs)
    assert any("openmm==8.6.1" in c for c in installs)
    assert any("mace-torch==0.3.17" in c and "torch>=2.5,<3" in c for c in installs)
    for c in installs:
        assert " -c " in c and "--torch-backend cpu" in c
    # The freeze, minus the editable line, is the constraint set.
    assert "constraints 2" in calls
    (qm,) = [c for c in calls if c.startswith("qm ")]
    assert qm.startswith("qm xtb crest ")
    assert "XTB_SPEC=xtb==6.7.1" in qm and "CREST_SPEC=crest==3.0.2" in qm
    assert f"ROOT={tmp_path}/.delfin/qm_tools" in qm
    for b in ("xtb", "tblite", "dftd4"):
        assert (tmp_path / ".delfin" / "qm_tools" / "bin" / b).exists(), b
    assert not list(root.glob(".pi-frozen.*")), "the constraints file is removed"


def test_a_group_that_is_there_is_left_alone(tmp_path):
    root = _harness(tmp_path, status="ok")
    proc = _run(tmp_path, root, "install_optional_groups")
    assert proc.returncode == 0, proc.stderr
    assert not [c for c in _calls(tmp_path) if c.startswith(("uv pip install", "qm "))]


def test_a_switch_turns_a_group_off(tmp_path):
    root = _harness(tmp_path)
    proc = _run(tmp_path, root, "install_optional_groups",
                env={"DELFIN_PI_PYSCF": "0", "DELFIN_PI_XTB": "0"})
    assert proc.returncode == 0, proc.stderr
    calls = _calls(tmp_path)
    assert not any("pyscf" in c for c in calls if c.startswith("uv pip install"))
    assert not any(c.startswith("qm ") for c in calls)
    assert any("openmm==8.6.1" in c for c in calls)
    assert "optional pyscf: off (DELFIN_PI_PYSCF=0)" in proc.stdout
    assert "optional xtb: off (DELFIN_PI_XTB=0)" in proc.stdout


def test_a_group_that_fails_warns_with_its_reason_and_the_rest_goes_on(tmp_path):
    root = _harness(tmp_path, uv_fail="openmm")
    proc = _run(tmp_path, root, "install_optional_groups; echo deploy-goes-on")
    assert proc.returncode == 0, proc.stderr
    assert "deploy-goes-on" in proc.stdout
    assert re.search(r"warning: optional openmm not installed: .*No solution found", proc.stderr)
    calls = _calls(tmp_path)
    assert any("mace-torch" in c for c in calls), "the groups after it still ran"
    assert any(c.startswith("qm ") for c in calls)


def test_the_smoke_step_names_the_groups_switched_off(tmp_path):
    root = _harness(tmp_path)
    proc = _run(tmp_path, root, "smoke_optional", env={"DELFIN_PI_OPENMM": "0"},
                skip="ai-retro orca pymol")
    assert proc.returncode == 0, proc.stderr
    (call,) = [c for c in _calls(tmp_path) if "smoke" in c]
    off = call.split("--off", 1)[1].split()
    assert sorted(off) == ["openmm", "opi"]


def test_the_shell_groups_are_the_tables_groups():
    from delfin.agent import optional_tools as ot
    shell = re.search(r'^PI_GROUPS="([^"]*)"', INSTALL, re.M).group(1).split()
    assert sorted(shell + ["opi"]) == sorted(g.name for g in ot.GROUPS)
    for g in ot.GROUPS:
        if g.name != "opi":
            assert g.switch == f"DELFIN_PI_{g.name.upper()}"


def test_the_readme_lists_every_group_and_switch():
    readme = (PI / "README.md").read_text(encoding="utf-8")
    section = readme[readme.index("## What the Pi can run"):]
    section = section[:section.index("\n## ", 5)]
    from delfin.agent import optional_tools as ot
    for g in ot.GROUPS:
        assert f"| `{g.name}` |" in section, g.name
        assert g.switch in section, g.switch
    assert "cluster" in section


# --- deploy_pi.sh forwards the switches -----------------------------------

def _pi_env(env: dict) -> subprocess.CompletedProcess:
    text = (PI / "deploy_pi.sh").read_text(encoding="utf-8")
    script = _function(text, "pi_env") + "\npi_env\n"
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                          env={"PATH": os.environ.get("PATH", "/usr/bin:/bin"), **env})


def test_the_deploy_forwards_the_pi_switches():
    proc = _pi_env({"DELFIN_PI_PYSCF": "0", "DELFIN_PI_SKIP_EXTRAS": "ai-retro orca pymol",
                    "OTHER": "1"})
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.startswith("env ")
    assert "DELFIN_PI_PYSCF='0'" in proc.stdout
    assert "DELFIN_PI_SKIP_EXTRAS='ai-retro orca pymol'" in proc.stdout
    assert "OTHER" not in proc.stdout
    assert _pi_env({}).stdout == ""


@pytest.mark.parametrize("value", ["0; rm -rf ~", "$(id)", "a'b", "x`y`"])
def test_a_switch_with_shell_syntax_is_refused(value):
    proc = _pi_env({"DELFIN_PI_PYSCF": value})
    assert proc.returncode != 0 and proc.stdout == ""
    assert "does not forward" in proc.stderr


def test_every_remote_call_carries_the_switches():
    text = (PI / "deploy_pi.sh").read_text(encoding="utf-8")
    calls = re.findall(r'ssh "\$HOST" "([^"]*pi_install\.sh[^"]*)"', text)
    assert len(calls) == 2
    for c in calls:
        assert c.startswith("${remote_env}bash ")


# --- the smoke step and the doctor ----------------------------------------

def test_a_smoke_line_parses_back():
    from delfin.agent import optional_tools as ot
    lines = [ot.format_smoke_line("pyscf", ot.OK, "water RHF/STO-3G -74.963023 Eh", 0.42),
             ot.format_smoke_line("xtb", ot.PARTIAL, "xtb 6.7.1; missing: crest"),
             ot.format_smoke_line("openmm", ot.OFF, "switched off (DELFIN_PI_OPENMM=0)"),
             "==> optional tools: smoke test", "smoke summary: 1 ok, 1 partial, 1 off"]
    parsed = ot.parse_smoke("\n".join(lines))
    assert parsed == {
        "pyscf": ("ok", "water RHF/STO-3G -74.963023 Eh (0.4 s)"),
        "xtb": ("partial", "xtb 6.7.1; missing: crest"),
        "openmm": ("off", "switched off (DELFIN_PI_OPENMM=0)"),
    }


def test_the_smoke_cli_reports_every_group_and_never_fails(monkeypatch, capsys):
    from delfin.agent import optional_tools as ot
    results = {"opi": (ot.OK, "x"), "pyscf": (ot.FAILED, "boom"), "torch": (ot.MISSING, "not installed")}

    def fake(group):
        if group.name == "xtb":
            raise RuntimeError("crashed")
        return results[group.name]
    monkeypatch.setattr(ot, "smoke", fake)
    assert ot.main(["smoke", "--off", " openmm"]) == 0
    parsed = ot.parse_smoke(capsys.readouterr().out)
    assert {k: v[0] for k, v in parsed.items()} == {
        "opi": "ok", "pyscf": "failed", "openmm": "off", "torch": "missing", "xtb": "failed"}
    assert "crashed" in parsed["xtb"][1]


def test_a_python_check_that_fails_is_failed_with_its_last_line(monkeypatch):
    from delfin.agent import optional_tools as ot
    group = ot.Group("demo", "", switch="DELFIN_PI_DEMO", install="", modules=("json",),
                     _code="raise SystemExit('BLAS symbol missing')")
    assert ot.smoke(group) == (ot.FAILED, "BLAS symbol missing")
    good = ot.Group("demo", "", switch="DELFIN_PI_DEMO", install="", modules=("json",),
                    _code='print(json.dumps({"detail": "fine"}))')
    assert ot.smoke(good) == (ot.OK, "fine")
    absent = ot.Group("demo", "", switch="DELFIN_PI_DEMO", install="",
                      modules=("delfin_no_such_module_xyz",), _code="")
    assert ot.smoke(absent)[0] == ot.MISSING


def test_the_doctor_shows_one_row_per_group(monkeypatch):
    from delfin.agent import doctor
    from delfin.agent import optional_tools as ot
    states = {"opi": (ot.OK, "opi 2.0.0"), "pyscf": (ot.MISSING, "not installed"),
              "openmm": (ot.OK, "openmm 8.6.1"), "torch": (ot.PARTIAL, "torch; missing: mace"),
              "xtb": (ot.MISSING, "not installed")}
    monkeypatch.setattr(ot, "availability", lambda g: states[g.name])
    rows = doctor._check_optional_tools({})
    assert [r["check"] for r in rows] == [f"optional: {g.name}" for g in ot.GROUPS]
    by = {r["check"].split(": ")[1]: r for r in rows}
    assert by["opi"]["status"] == doctor.PASS and by["openmm"]["status"] == doctor.PASS
    for name in ("pyscf", "torch", "xtb"):
        assert by[name]["status"] == doctor.WARN
        assert f"DELFIN_PI_{name.upper()}" in by[name]["fix"]
    assert ("optional chemistry", "_check_optional_tools") in doctor._CHECK_ATTRS


def test_availability_runs_nothing(monkeypatch):
    """The doctor's probe imports and starts nothing: torch alone takes
    seconds to import on a Pi."""
    from delfin.agent import optional_tools as ot
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: pytest.fail("ran a process"))
    before = set(sys.modules)
    for g in ot.GROUPS:
        ot.availability(g)
    assert not {"torch", "pyscf", "openmm", "mace"} & (set(sys.modules) - before)
