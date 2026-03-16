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
