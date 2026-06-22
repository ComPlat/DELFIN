import subprocess
from pathlib import Path

from delfin.dashboard import browser_workflows


def test_require_tools_auto_installs_supported_analysis_tools(monkeypatch, tmp_path):
    state = {"installed": False, "extra_env": None}

    def fake_resolved(tool_name):
        if tool_name in {"crest", "xtb", "orca"}:
            return f"/usr/bin/{tool_name}"
        if state["installed"] and tool_name in {"censo", "c2anmr", "anmr"}:
            return f"/env/bin/{tool_name}"
        return ""

    def fake_installer(*, target_dir=None, extra_env=None):
        state["installed"] = True
        state["extra_env"] = dict(extra_env or {})
        return Path("/tmp/analysis_tools"), subprocess.CompletedProcess(
            args=["bash", "install_analysis_tools.sh"],
            returncode=0,
            stdout="installed",
        )

    monkeypatch.setattr(browser_workflows, "_resolved_path_or_empty", fake_resolved)
    monkeypatch.setattr(browser_workflows, "run_analysis_tools_installer", fake_installer)
    monkeypatch.setenv("DELFIN_AUTO_INSTALL_ANALYSIS_TOOLS", "1")

    resolved = browser_workflows._require_tools(
        ["crest", "censo", "c2anmr", "anmr", "orca", "xtb"],
        workdir=tmp_path,
    )

    assert resolved["censo"] == "/env/bin/censo"
    assert resolved["c2anmr"] == "/env/bin/c2anmr"
    assert resolved["anmr"] == "/env/bin/anmr"
    assert state["extra_env"]["INSTALL_CENSO"] == "1"
    assert state["extra_env"]["INSTALL_ANMR"] == "1"
    assert (tmp_path / "analysis_tools_auto_install.log").is_file()


def test_require_tools_reports_unsupported_missing_tools_after_auto_install(monkeypatch, tmp_path):
    state = {"installed": False}

    def fake_resolved(tool_name):
        if tool_name in {"xtb", "orca"}:
            return f"/usr/bin/{tool_name}"
        if state["installed"] and tool_name in {"censo", "c2anmr", "anmr"}:
            return f"/env/bin/{tool_name}"
        return ""

    def fake_installer(*, target_dir=None, extra_env=None):
        state["installed"] = True
        return Path("/tmp/analysis_tools"), subprocess.CompletedProcess(
            args=["bash", "install_analysis_tools.sh"],
            returncode=0,
            stdout="installed",
        )

    monkeypatch.setattr(browser_workflows, "_resolved_path_or_empty", fake_resolved)
    monkeypatch.setattr(browser_workflows, "run_analysis_tools_installer", fake_installer)

    try:
        browser_workflows._require_tools(
            ["crest", "censo", "c2anmr", "anmr", "orca", "xtb"],
            workdir=tmp_path,
        )
    except RuntimeError as exc:
        message = str(exc)
    else:
        raise AssertionError("Expected RuntimeError for unsupported missing tool")

    assert "crest" in message
    assert "unsupported tools such as crest, xtb, and orca" in message


def test_require_tools_can_disable_auto_install(monkeypatch, tmp_path):
    state = {"called": False}

    def fake_installer(*, target_dir=None, extra_env=None):
        state["called"] = True
        return Path("/tmp/analysis_tools"), subprocess.CompletedProcess(
            args=["bash", "install_analysis_tools.sh"],
            returncode=0,
            stdout="installed",
        )

    monkeypatch.setattr(browser_workflows, "_resolved_path_or_empty", lambda tool_name: "")
    monkeypatch.setattr(browser_workflows, "run_analysis_tools_installer", fake_installer)
    monkeypatch.setenv("DELFIN_AUTO_INSTALL_ANALYSIS_TOOLS", "0")

    try:
        browser_workflows._require_tools(["censo", "anmr"], workdir=tmp_path)
    except RuntimeError as exc:
        message = str(exc)
    else:
        raise AssertionError("Expected RuntimeError when auto-install is disabled")

    assert "censo, anmr" in message
    assert state["called"] is False


def test_require_tools_auto_upgrades_outdated_censo(monkeypatch, tmp_path):
    state = {"installed": False, "extra_env": None}

    def fake_resolved(tool_name):
        if tool_name == "censo":
            return "/env/bin/censo"
        if tool_name in {"xtb", "orca"}:
            return f"/usr/bin/{tool_name}"
        return ""

    def fake_detect(_path):
        return browser_workflows._minimum_supported_censo_version() if state["installed"] else (2, 1, 2)

    def fake_installer(*, target_dir=None, extra_env=None):
        state["installed"] = True
        state["extra_env"] = dict(extra_env or {})
        return Path("/tmp/analysis_tools"), subprocess.CompletedProcess(
            args=["bash", "install_analysis_tools.sh"],
            returncode=0,
            stdout="updated",
        )

    monkeypatch.setattr(browser_workflows, "_resolved_path_or_empty", fake_resolved)
    monkeypatch.setattr(browser_workflows, "_detect_censo_version", fake_detect)
    monkeypatch.setattr(browser_workflows, "run_analysis_tools_installer", fake_installer)

    resolved = browser_workflows._require_tools(
        ["censo", "xtb", "orca"],
        workdir=tmp_path,
    )

    assert resolved["censo"] == "/env/bin/censo"
    assert state["extra_env"]["INSTALL_CENSO"] == "1"
    assert state["extra_env"]["FORCE_REINSTALL"] == "1"
    assert state["extra_env"]["CENSO_PREFER_LATEST"] == "1"


def test_helper_launch_command_uses_shell_for_shell_script(tmp_path):
    script = tmp_path / "c2anmr"
    script.write_text("#!/bin/bash\nmkdir -p anmr\n", encoding="utf-8")
    command = browser_workflows._helper_launch_command(str(script))
    assert command == ["bash", str(script.resolve())]


def test_build_censo_command_uses_legacy_cli_for_censo_2(monkeypatch):
    monkeypatch.setattr(browser_workflows, "_detect_censo_version", lambda _path: (2, 1, 2))
    command = browser_workflows._build_censo_command(
        censo_path="/usr/bin/censo",
        ensemble_name="crest_conformers.xyz",
        charge=0,
        multiplicity=1,
        rc_name="workflow.censo2rc",
        pal=12,
        solvent="chcl3",
    )
    assert "--prescreening" not in command
    assert "-O" in command
    assert "-s" in command


def test_build_censo_command_uses_stage_flags_for_censo_3(monkeypatch):
    monkeypatch.setattr(browser_workflows, "_detect_censo_version", lambda _path: (3, 0, 6))
    command = browser_workflows._build_censo_command(
        censo_path="/usr/bin/censo",
        ensemble_name="crest_conformers.xyz",
        charge=0,
        multiplicity=1,
        rc_name="workflow.censo2rc",
        pal=12,
        solvent="chcl3",
    )
    assert "--prescreening" in command
    assert "--nmr" in command
    assert "--omp-min" in command
    assert "--solvent" in command


def test_censo_refinement_threshold_is_positive_for_censo_3(monkeypatch):
    monkeypatch.setattr(browser_workflows, "_detect_censo_version", lambda _path: (3, 0, 6))

    threshold = browser_workflows._censo_refinement_threshold("/env/bin/censo")

    assert threshold > 0.0


def test_resolve_orca_reference_result_path_falls_back_to_normal_log(tmp_path):
    log_path = tmp_path / "tms_reference.log"
    log_path.write_text("...\n****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")

    resolved = browser_workflows._resolve_orca_reference_result_path(
        workdir=tmp_path,
        stem="tms_reference",
    )

    assert resolved == log_path


def test_ensure_censo_nmr_coords_backfills_missing_coord_from_nmr_input(tmp_path):
    conf_dir = tmp_path / "4_NMR" / "CONF4" / "nmr"
    conf_dir.mkdir(parents=True)
    (conf_dir / "nmr.inp").write_text(
        "! pbe0 def2-tzvp\n"
        "* xyz 0 1\n"
        "C 0.0 0.0 0.0\n"
        "H 0.0 0.0 1.0\n"
        "*\n",
        encoding="utf-8",
    )

    browser_workflows._ensure_censo_nmr_coords(tmp_path)

    coord_text = (tmp_path / "4_NMR" / "CONF4" / "coord").read_text(encoding="utf-8")
    assert coord_text.startswith("$coord\n")
    assert "c" in coord_text
    assert "h" in coord_text


def test_ensure_anmr_coord_inputs_copies_coords_into_anmr_workspace(tmp_path):
    censo_conf = tmp_path / "4_NMR" / "CONF4"
    censo_conf.mkdir(parents=True)
    coord_text = "$coord\n  0.0 0.0 0.0 c\n$end\n"
    (censo_conf / "coord").write_text(coord_text, encoding="utf-8")

    anmr_conf = tmp_path / "anmr" / "CONF4"
    anmr_conf.mkdir(parents=True)

    browser_workflows._ensure_anmr_coord_inputs(workdir=tmp_path, anmr_dir=tmp_path / "anmr")

    assert (tmp_path / "anmr" / "coord").read_text(encoding="utf-8") == coord_text
    assert (tmp_path / "anmr" / "CONF4" / "coord").read_text(encoding="utf-8") == coord_text


def test_render_anmr_png_writes_png_from_anmr_dat(tmp_path):
    data_path = tmp_path / "anmr.dat"
    data_path.write_text(
        "-0.5 0.0\n"
        "0.0 1.0\n"
        "0.5 0.0\n",
        encoding="utf-8",
    )

    output = browser_workflows._render_anmr_png(
        data_path,
        tmp_path / "anmr_spectrum.png",
        title="ANMR Spectrum: test",
    )

    assert output is not None
    assert output.is_file()


def test_resume_detectors_identify_completed_outputs(tmp_path):
    workdir = tmp_path
    (workdir / "crest_conformers.xyz").write_text("2\nx\nH 0 0 0\nH 0 0 1\n", encoding="utf-8")
    (workdir / "anmr_nucinfo").write_text("nucinfo\n", encoding="utf-8")
    (workdir / "anmr_rotamer").write_text("rotamer\n", encoding="utf-8")
    (workdir / "censo.out").write_text("...\nCENSO all done!\n", encoding="utf-8")
    anmr_dir = workdir / "anmr"
    anmr_dir.mkdir()
    (anmr_dir / "anmr_enso").write_text("enso\n", encoding="utf-8")
    (anmr_dir / "anmr.dat").write_text("0.0 1.0\n", encoding="utf-8")
    (anmr_dir / "anmr.out").write_text("...\nAll done.\n", encoding="utf-8")
    (anmr_dir / "anmr_spectrum.png").write_text("png\n", encoding="utf-8")
    (workdir / "tms_reference.log").write_text("...\nORCA TERMINATED NORMALLY\n", encoding="utf-8")
    (anmr_dir / ".anmrrc").write_text("[refs]\n", encoding="utf-8")

    assert browser_workflows._crest_outputs_complete(workdir) is True
    assert browser_workflows._censo_outputs_complete(workdir) is True
    assert browser_workflows._c2anmr_outputs_complete(workdir) is True
    assert browser_workflows._tms_reference_complete(workdir, anmr_dir) is True
    assert browser_workflows._anmr_outputs_complete(anmr_dir) is True
    assert browser_workflows._existing_anmr_plot(anmr_dir).endswith("anmr_spectrum.png")
