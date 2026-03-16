from pathlib import Path

from delfin.runtime_setup import (
    collect_runtime_diagnostics,
    get_packaged_submit_templates_dir,
    get_user_qm_tools_dir,
    resolve_backend_choice,
    resolve_orca_base,
    resolve_submit_templates_dir,
    stage_packaged_qm_tools,
)


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
    assert by_name["orca"]["status"] in {"ok", "missing"}


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
