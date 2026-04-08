import shutil
from pathlib import Path

from delfin.runtime_setup import (
    collect_bwunicluster_verification,
    collect_runtime_diagnostics,
    detect_local_runtime_limits,
    describe_orca_installation,
    discover_orca_installations,
    get_packaged_bwunicluster_install_script,
    get_packaged_submit_templates_dir,
    get_user_qm_tools_dir,
    resolve_backend_choice,
    resolve_orca_base,
    resolve_submit_templates_dir,
    stage_packaged_qm_tools,
)
import delfin.runtime_setup as runtime_setup


def test_resolve_backend_choice_prefers_explicit_value():
    assert resolve_backend_choice("local", "slurm", slurm_available=True) == "local"
    assert resolve_backend_choice("slurm", "local", slurm_available=False) == "slurm"


def test_resolve_backend_choice_uses_settings_or_auto_detection():
    assert resolve_backend_choice(None, "slurm", slurm_available=False) == "slurm"
    assert resolve_backend_choice(None, "auto", slurm_available=True) == "slurm"
    assert resolve_backend_choice(None, "", slurm_available=False) == "local"


def test_resolve_orca_base_prefers_backend_specific_setting(tmp_path):
    runtime_settings = {
        "orca_base": str(tmp_path / "global_orca"),
        "local": {"orca_base": str(tmp_path / "local_orca")},
        "slurm": {"orca_base": str(tmp_path / "slurm_orca")},
    }

    assert resolve_orca_base(None, runtime_settings, "local", auto_candidates=[]) == str(
        tmp_path / "local_orca"
    )
    assert resolve_orca_base(None, runtime_settings, "slurm", auto_candidates=[]) == str(
        tmp_path / "slurm_orca"
    )


def test_resolve_orca_base_has_no_hardcoded_local_fallback_by_default(monkeypatch):
    monkeypatch.setattr(runtime_setup, "discover_orca_installations", lambda search_roots=None: [])
    monkeypatch.setattr(shutil, "which", lambda name: None)
    monkeypatch.delenv("DELFIN_ORCA_BASE", raising=False)
    monkeypatch.delenv("ORCA_BINARY", raising=False)
    monkeypatch.delenv("ORCA_PATH", raising=False)
    assert resolve_orca_base(None, {}, "local", auto_candidates=[]) == ""


def test_packaged_submit_templates_dir_contains_required_scripts():
    submit_dir = get_packaged_submit_templates_dir()

    assert (submit_dir / "submit_delfin.sh").is_file()
    assert (submit_dir / "submit_turbomole.sh").is_file()


def test_resolve_submit_templates_dir_uses_packaged_fallback():
    fallback = get_packaged_submit_templates_dir()
    resolved = resolve_submit_templates_dir({}, fallback)

    assert resolved == fallback


def test_collect_runtime_diagnostics_reports_missing_orca_but_existing_templates():
    submit_dir = get_packaged_submit_templates_dir()
    diagnostics = collect_runtime_diagnostics(
        {"qm_tools_root": ""},
        backend="slurm",
        effective_orca_base="",
        submit_templates_dir=submit_dir,
    )

    by_name = {item["name"]: item for item in diagnostics}
    assert by_name["backend"]["detail"] == "slurm"
    assert by_name["slurm-templates"]["status"] == "ok"
    assert by_name["orca"]["status"] in {"ok", "missing", "system"}


def test_stage_packaged_qm_tools_copies_bundle_to_target(tmp_path):
    target = tmp_path / "user_qm_tools"

    staged = stage_packaged_qm_tools(target)

    assert staged == target
    assert (staged / "install_qm_tools.sh").is_file()
    assert (staged / "bin" / "xtb4stda").is_file()
    assert (staged / "bin" / "stda_v1.6.1").is_file()
    assert (staged / "share" / "xtb4stda" / ".xtb4stdarc").is_file()


def test_get_user_qm_tools_dir_accepts_override(tmp_path):
    target = tmp_path / "custom_qm_tools"

    assert get_user_qm_tools_dir(target) == target


def test_discover_orca_installations_finds_multiple_versions(tmp_path):
    root = tmp_path / "software"
    orca_611 = root / "orca_6_1_1"
    orca_601 = root / "orca_6_0_1"
    orca_611.mkdir(parents=True)
    orca_601.mkdir(parents=True)
    (orca_611 / "orca").write_text("", encoding="utf-8")
    (orca_601 / "orca").write_text("", encoding="utf-8")

    discovered = discover_orca_installations(search_roots=[root])

    assert str(orca_611) in discovered
    assert str(orca_601) in discovered


def test_describe_orca_installation_formats_version():
    assert describe_orca_installation("/opt/orca_6_1_1") == "ORCA 6.1.1"


def test_detect_local_runtime_limits_returns_positive_values():
    cores, ram_mb = detect_local_runtime_limits()

    assert cores >= 1
    assert ram_mb >= 1


def test_packaged_bwunicluster_installer_exists():
    installer = get_packaged_bwunicluster_install_script()

    assert installer is not None
    assert installer.is_file()


def test_collect_bwunicluster_verification_reports_expected_keys(tmp_path):
    checks = collect_bwunicluster_verification(
        repo_dir=tmp_path / "missing_repo",
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
    )

    by_name = {item["name"]: item for item in checks}
    assert "install-script" in by_name
    assert "submit-templates" in by_name
    assert "sbatch" in by_name
