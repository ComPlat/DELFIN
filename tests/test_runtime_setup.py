import shutil
from pathlib import Path

import pytest

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


def _qm_tools_bundle_present() -> bool:
    """The xtb4stda / stda binaries are a large optional bundle added at
    packaging time — they are NOT in a plain source checkout (or CI). Detect
    them so this test runs in a packaged/release env but skips elsewhere
    instead of being a permanent expected-failure."""
    try:
        return (runtime_setup.get_packaged_qm_tools_dir()
                / "bin" / "xtb4stda").is_file()
    except Exception:
        return False


@pytest.mark.skipif(
    not _qm_tools_bundle_present(),
    reason="packaged qm_tools binary bundle (xtb4stda/stda) not present in this checkout",
)
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


def test_staging_the_tools_keeps_the_links_and_survives_a_broken_one(tmp_path):
    """A copied binary is not the binary.

    ``bin/xtb`` points into the environment xtb was installed in. Followed
    instead of copied it became a 6.4 MB file of its own, outside the
    environment holding the libraries it is linked against, and it did not run
    at all: "libmctc-lib.so.0: cannot open shared object file". A resolved tool
    is looked for in this directory before anywhere else, so that copy then
    beat the working xtb on the PATH -- pressing Prepare left the user with a
    tool chain that could not start.

    And a link whose target has been moved away raised out of the staging
    before the installer had run a line, which took every button with it.
    """
    import os
    import shutil

    from delfin import runtime_setup

    source = tmp_path / "packaged"
    (source / "bin").mkdir(parents=True)
    real = tmp_path / "elsewhere" / "xtb"
    real.parent.mkdir()
    real.write_text("#!/bin/sh\necho xtb\n")
    real.chmod(0o755)
    os.symlink(real, source / "bin" / "xtb")
    os.symlink(tmp_path / "gone" / "crest", source / "bin" / "crest")

    monkey = tmp_path / "staged"
    original = runtime_setup.get_packaged_qm_tools_dir
    runtime_setup.get_packaged_qm_tools_dir = lambda: source
    try:
        staged = runtime_setup.stage_packaged_qm_tools(monkey)
    finally:
        runtime_setup.get_packaged_qm_tools_dir = original

    linked = staged / "bin" / "xtb"
    assert linked.is_symlink(), "the binary was copied out of its environment"
    assert os.path.realpath(linked) == str(real)
    assert (staged / "bin" / "crest").is_symlink(), "a dead link is not fatal"
    # And the whole staging completed: the dangling one did not stop it.
    assert linked.exists()

    shutil.rmtree(staged, ignore_errors=True)
