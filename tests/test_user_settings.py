import importlib.util
import json
from pathlib import Path


_MODULE_PATH = Path(__file__).resolve().parents[1] / "delfin" / "user_settings.py"
_SPEC = importlib.util.spec_from_file_location("delfin_user_settings", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load user_settings module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

get_settings_path = _MODULE.get_settings_path
load_settings = _MODULE.load_settings
load_remote_archive_enabled = _MODULE.load_remote_archive_enabled
load_runtime_settings = _MODULE.load_runtime_settings
load_transfer_settings = _MODULE.load_transfer_settings
normalize_local_directory_setting = _MODULE.normalize_local_directory_setting
save_remote_archive_enabled = _MODULE.save_remote_archive_enabled
save_transfer_settings = _MODULE.save_transfer_settings


def test_save_and_load_transfer_settings_roundtrip(tmp_path):
    settings_path = tmp_path / "settings.json"
    saved = save_transfer_settings(
        "cluster.example.org",
        "alice",
        "/remote/archive",
        2222,
        settings_path=settings_path,
    )
    loaded = load_transfer_settings(settings_path)

    assert settings_path == get_settings_path(settings_path)
    assert saved["settings_path"] == str(settings_path)
    assert loaded["settings_path"] == str(settings_path)
    assert loaded["host"] == "cluster.example.org"
    assert loaded["user"] == "alice"
    assert loaded["remote_path"] == "/remote/archive"
    assert loaded["port"] == 2222
    assert settings_path.stat().st_mode & 0o777 == 0o600


def test_load_settings_migrates_legacy_transfer_file(tmp_path):
    settings_path = tmp_path / "settings.json"
    legacy_path = tmp_path / _MODULE.LEGACY_TRANSFER_CONFIG_NAME
    legacy_path.write_text(
        json.dumps(
            {
                "host": "cluster.example.org",
                "user": "alice",
                "remote_path": "/remote/archive",
                "port": 2222,
            }
        ),
        encoding="utf-8",
    )

    loaded = load_settings(settings_path)

    assert loaded["transfer"]["host"] == "cluster.example.org"
    assert loaded["features"]["remote_archive_enabled"] is False
    assert settings_path.exists()
    persisted = json.loads(settings_path.read_text(encoding="utf-8"))
    assert persisted["transfer"]["user"] == "alice"
    assert persisted["transfer"]["port"] == 2222


def test_load_settings_appends_missing_defaults_without_overwriting_existing(tmp_path):
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(
        json.dumps(
            {
                "transfer": {
                    "host": "cluster.example.org",
                    "user": "alice",
                    "remote_path": "/remote/archive",
                    "port": 2222,
                }
            }
        ),
        encoding="utf-8",
    )

    original_defaults = _MODULE.DEFAULT_SETTINGS
    try:
        _MODULE.DEFAULT_SETTINGS = {
            "transfer": {},
            "runtime": {
                "backend": "auto",
                "orca_base": "",
                "qm_tools_root": "",
                "local": {
                    "orca_base": "",
                    "max_cores": 384,
                    "max_ram_mb": 1_400_000,
                    "allow_oversubscribe": False,
                    "oversubscribe_factor": 1.0,
                    "allow_live_load_bypass": False,
                    "live_cpu_target_factor": 0.95,
                    "live_min_free_ram_mb": 64_000,
                },
                "slurm": {"orca_base": "", "submit_templates_dir": "", "profile": ""},
            },
            "features": {"remote_archive_enabled": False},
            "ui": {"show_hidden_files": False},
        }
        loaded = load_settings(settings_path)
    finally:
        _MODULE.DEFAULT_SETTINGS = original_defaults

    assert loaded["transfer"]["host"] == "cluster.example.org"
    assert loaded["ui"]["show_hidden_files"] is False
    assert loaded["features"]["remote_archive_enabled"] is False

    persisted = json.loads(settings_path.read_text(encoding="utf-8"))
    assert persisted["transfer"]["host"] == "cluster.example.org"
    assert persisted["ui"]["show_hidden_files"] is False


def test_remote_archive_flag_defaults_to_false_and_persists(tmp_path):
    settings_path = tmp_path / "settings.json"

    assert load_remote_archive_enabled(settings_path) is False

    saved = save_remote_archive_enabled(True, settings_path=settings_path)
    loaded = load_settings(settings_path)

    assert saved is True
    assert load_remote_archive_enabled(settings_path) is True
    assert loaded["features"]["remote_archive_enabled"] is True


def test_load_settings_normalizes_workspace_paths(tmp_path):
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(
        json.dumps(
            {
                "paths": {
                    "calculations_dir": "~/calc_custom",
                    "archive_dir": str(tmp_path / "archive_custom"),
                }
            }
        ),
        encoding="utf-8",
    )

    loaded = load_settings(settings_path)

    assert loaded["paths"]["calculations_dir"] == str(Path("~/calc_custom").expanduser())
    assert loaded["paths"]["archive_dir"] == str((tmp_path / "archive_custom").expanduser())


def test_normalize_local_directory_setting_allows_blank_and_rejects_controls():
    assert normalize_local_directory_setting("", "Calculations path") == ""
    assert normalize_local_directory_setting("  ~/calc  ", "Calculations path") == str(
        Path("~/calc").expanduser()
    )
    try:
        normalize_local_directory_setting("bad\npath", "Calculations path")
    except ValueError as exc:
        assert "unsupported control characters" in str(exc)
    else:
        raise AssertionError("Expected ValueError for control characters")


def test_load_settings_normalizes_runtime_payload(tmp_path):
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(
        json.dumps(
            {
                "runtime": {
                    "backend": "SLURM",
                    "orca_base": "~/orca_global",
                    "qm_tools_root": "~/qm_tools",
                    "local": {
                        "orca_base": "~/orca_local",
                        "max_cores": 32,
                        "max_ram_mb": 128000,
                        "allow_live_load_bypass": True,
                        "live_cpu_target_factor": 0.9,
                        "live_min_free_ram_mb": 96000,
                    },
                    "slurm": {
                        "orca_base": "~/orca_slurm",
                        "submit_templates_dir": "~/submit_templates",
                        "profile": "bwunicluster3",
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    loaded = load_settings(settings_path)

    assert loaded["runtime"]["backend"] == "slurm"
    assert loaded["runtime"]["orca_base"] == str(Path("~/orca_global").expanduser())
    assert loaded["runtime"]["qm_tools_root"] == str(Path("~/qm_tools").expanduser())
    assert loaded["runtime"]["local"]["orca_base"] == str(Path("~/orca_local").expanduser())
    assert loaded["runtime"]["local"]["max_cores"] == 32
    assert loaded["runtime"]["local"]["max_ram_mb"] == 128000
    assert loaded["runtime"]["local"]["allow_live_load_bypass"] is True
    assert loaded["runtime"]["local"]["live_cpu_target_factor"] == 0.9
    assert loaded["runtime"]["local"]["live_min_free_ram_mb"] == 96000
    assert loaded["runtime"]["slurm"]["orca_base"] == str(Path("~/orca_slurm").expanduser())
    assert loaded["runtime"]["slurm"]["submit_templates_dir"] == str(
        Path("~/submit_templates").expanduser()
    )
    assert loaded["runtime"]["slurm"]["profile"] == "bwunicluster3"


def test_load_runtime_settings_returns_defaults_when_missing(tmp_path):
    settings_path = tmp_path / "settings.json"

    runtime = load_runtime_settings(settings_path)

    assert runtime["backend"] == "auto"
    assert runtime["local"]["max_cores"] == _MODULE.DEFAULT_SETTINGS["runtime"]["local"]["max_cores"]
    assert runtime["local"]["max_ram_mb"] == _MODULE.DEFAULT_SETTINGS["runtime"]["local"]["max_ram_mb"]
    assert (
        runtime["local"]["allow_live_load_bypass"]
        == _MODULE.DEFAULT_SETTINGS["runtime"]["local"]["allow_live_load_bypass"]
    )
    assert runtime["slurm"]["submit_templates_dir"] == ""


def test_load_settings_normalizes_ui_tabs_payload(tmp_path):
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(
        json.dumps(
            {
                "ui": {
                    "tabs": {
                        "order": ["submit", "archive", "submit", "", None],
                        "hidden": ["archive", "archive", "", None],
                    }
                }
            }
        ),
        encoding="utf-8",
    )

    loaded = load_settings(settings_path)

    assert loaded["ui"]["tabs"]["order"] == ["submit", "archive"]
    assert loaded["ui"]["tabs"]["hidden"] == ["archive"]
