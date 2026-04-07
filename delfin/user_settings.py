"""User-local DELFIN settings stored outside the git repo."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

try:
    import psutil  # type: ignore
except ImportError:
    psutil = None  # type: ignore

SETTINGS_FILE_NAME = ".delfin_settings.json"
LEGACY_TRANSFER_CONFIG_NAME = ".delfin_transfer_target.json"


def _detect_default_local_limits():
    cpu_count = os.cpu_count() or 1
    memory_mb = 1_400_000
    if psutil is not None:
        try:
            memory_mb = max(1, int(psutil.virtual_memory().total // (1024 * 1024)))
        except Exception:
            memory_mb = 1_400_000
    return max(1, int(cpu_count)), memory_mb


_DEFAULT_LOCAL_CORES, _DEFAULT_LOCAL_RAM_MB = _detect_default_local_limits()
DEFAULT_SETTINGS = {
    "transfer": {
        "ssh_control_path": "",
    },
    "paths": {},
    "runtime": {
        "backend": "auto",
        "orca_base": "",
        "qm_tools_root": "",
        "csp_tools_root": "",
        "mlp_tools_root": "",
        "tool_binaries": {},
        "local": {
            "orca_base": "",
            "max_cores": _DEFAULT_LOCAL_CORES,
            "max_ram_mb": _DEFAULT_LOCAL_RAM_MB,
            "allow_oversubscribe": False,
            "oversubscribe_factor": 1.0,
            "allow_live_load_bypass": False,
            "live_cpu_target_factor": 0.95,
            "live_min_free_ram_mb": 64_000,
        },
        "slurm": {
            "orca_base": "",
            "submit_templates_dir": "",
            "profile": "",
        },
    },
    "features": {
        "remote_archive_enabled": False,
    },
    "ui": {
        "tabs": {
            "order": [],
            "hidden": [],
        },
    },
    "agent": {
        "backend": "cli",
        "model": "sonnet",
        "max_tokens": 4096,
        "role_models": {
            "session_manager": "auto",
            "critic_agent": "auto",
            "runtime_agent": "auto",
            "reviewer_agent": "auto",
            "builder_agent": "auto",
            "test_agent": "auto",
            "solo_agent": "auto",
        },
    },
}


def utc_now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def normalize_ssh_transfer_settings(host, user, remote_path, port):
    """Validate and normalize SSH transfer settings."""
    host = str(host or "").strip()
    if not host or any(ch.isspace() for ch in host) or "@" in host:
        raise ValueError("Provide a valid SSH host or IP address.")

    user = str(user or "").strip()
    if not user or any(ch.isspace() for ch in user) or any(ch in user for ch in "@:"):
        raise ValueError("Provide a valid remote account name.")

    remote_path = str(remote_path or "").strip()
    if not remote_path:
        raise ValueError("Provide a remote destination folder.")
    if "\x00" in remote_path or "\n" in remote_path or "\r" in remote_path:
        raise ValueError("Remote destination contains unsupported control characters.")

    try:
        port = int(port)
    except Exception as exc:
        raise ValueError("SSH port must be a number.") from exc
    if not 1 <= port <= 65535:
        raise ValueError("SSH port must be between 1 and 65535.")

    return host, user, remote_path, port


def normalize_local_directory_setting(path_value, label):
    """Validate and normalize an optional local directory setting."""
    path_value = str(path_value or "").strip()
    if not path_value:
        return ""
    if "\x00" in path_value or "\n" in path_value or "\r" in path_value:
        raise ValueError(f"{label} contains unsupported control characters.")
    return str(Path(path_value).expanduser())


def normalize_choice_setting(value, label, allowed_values, default):
    normalized = str(value or "").strip().lower()
    if not normalized:
        return default
    if normalized not in allowed_values:
        allowed = ", ".join(sorted(allowed_values))
        raise ValueError(f"{label} must be one of: {allowed}.")
    return normalized


def normalize_positive_int_setting(value, label, default, minimum=1):
    if value in ("", None):
        return int(default)
    try:
        normalized = int(value)
    except Exception as exc:
        raise ValueError(f"{label} must be a whole number.") from exc
    if normalized < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return normalized


def normalize_positive_float_setting(value, label, default, minimum=1.0):
    if value in ("", None):
        return float(default)
    try:
        normalized = float(value)
    except Exception as exc:
        raise ValueError(f"{label} must be a number.") from exc
    if normalized < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return float(normalized)


def normalize_bool_setting(value, default=False):
    if value in ("", None):
        return bool(default)
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "on", "enable", "enabled"}:
        return True
    if normalized in {"0", "false", "no", "off", "disable", "disabled"}:
        return False
    return bool(default)


def normalize_string_list_setting(value, label):
    if value in (None, ""):
        return []
    if not isinstance(value, list):
        raise ValueError(f"{label} must be a JSON array.")
    normalized = []
    seen = set()
    for item in value:
        text = str(item or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        normalized.append(text)
    return normalized


def get_settings_path(base_path=None):
    return Path(base_path).expanduser() if base_path else (Path.home() / SETTINGS_FILE_NAME)


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write_json_atomic(path, payload):
    target = Path(path)
    tmp_path = target.with_suffix(target.suffix + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    try:
        tmp_path.chmod(0o600)
    except Exception:
        pass
    tmp_path.replace(target)
    try:
        target.chmod(0o600)
    except Exception:
        pass


def _normalized_settings_dict(payload):
    if payload is None:
        payload = {}
    if not isinstance(payload, dict):
        raise ValueError("Settings file must contain a JSON object.")
    normalized = dict(payload)
    transfer = normalized.get("transfer", {})
    if transfer is None:
        transfer = {}
    if not isinstance(transfer, dict):
        raise ValueError("Settings key 'transfer' must be a JSON object.")
    normalized["transfer"] = dict(transfer)
    paths = normalized.get("paths", {})
    if paths is None:
        paths = {}
    if not isinstance(paths, dict):
        raise ValueError("Settings key 'paths' must be a JSON object.")
    normalized_paths = {}
    for key, label in (
        ("calculations_dir", "Calculations path"),
        ("archive_dir", "Archive path"),
    ):
        normalized_paths[key] = normalize_local_directory_setting(
            paths.get(key, ""),
            label,
        )
    normalized["paths"] = normalized_paths
    runtime = normalized.get("runtime", {})
    if runtime is None:
        runtime = {}
    if not isinstance(runtime, dict):
        raise ValueError("Settings key 'runtime' must be a JSON object.")
    local_runtime = runtime.get("local", {})
    if local_runtime is None:
        local_runtime = {}
    if not isinstance(local_runtime, dict):
        raise ValueError("Settings key 'runtime.local' must be a JSON object.")
    slurm_runtime = runtime.get("slurm", {})
    if slurm_runtime is None:
        slurm_runtime = {}
    if not isinstance(slurm_runtime, dict):
        raise ValueError("Settings key 'runtime.slurm' must be a JSON object.")
    tool_binaries = runtime.get("tool_binaries", {})
    if tool_binaries is None:
        tool_binaries = {}
    if not isinstance(tool_binaries, dict):
        raise ValueError("Settings key 'runtime.tool_binaries' must be a JSON object.")
    normalized_tool_binaries = {}
    for key, value in tool_binaries.items():
        tool_name = str(key or "").strip().lower()
        if not tool_name:
            continue
        normalized_tool_binaries[tool_name] = normalize_local_directory_setting(
            value,
            f"{tool_name} binary",
        )
    normalized["runtime"] = {
        "backend": normalize_choice_setting(
            runtime.get("backend", DEFAULT_SETTINGS["runtime"]["backend"]),
            "Runtime backend",
            {"auto", "local", "slurm"},
            DEFAULT_SETTINGS["runtime"]["backend"],
        ),
        "orca_base": normalize_local_directory_setting(
            runtime.get("orca_base", ""),
            "Global ORCA path",
        ),
        "qm_tools_root": normalize_local_directory_setting(
            runtime.get("qm_tools_root", ""),
            "qm_tools root",
        ),
        "csp_tools_root": normalize_local_directory_setting(
            runtime.get("csp_tools_root", ""),
            "csp_tools root",
        ),
        "mlp_tools_root": normalize_local_directory_setting(
            runtime.get("mlp_tools_root", ""),
            "mlp_tools root",
        ),
        "tool_binaries": normalized_tool_binaries,
        "local": {
            "orca_base": normalize_local_directory_setting(
                local_runtime.get("orca_base", ""),
                "Local ORCA path",
            ),
            "max_cores": normalize_positive_int_setting(
                local_runtime.get("max_cores", DEFAULT_SETTINGS["runtime"]["local"]["max_cores"]),
                "Local max cores",
                DEFAULT_SETTINGS["runtime"]["local"]["max_cores"],
            ),
            "max_ram_mb": normalize_positive_int_setting(
                local_runtime.get("max_ram_mb", DEFAULT_SETTINGS["runtime"]["local"]["max_ram_mb"]),
                "Local max RAM (MB)",
                DEFAULT_SETTINGS["runtime"]["local"]["max_ram_mb"],
            ),
            "allow_oversubscribe": normalize_bool_setting(
                local_runtime.get(
                    "allow_oversubscribe",
                    DEFAULT_SETTINGS["runtime"]["local"]["allow_oversubscribe"],
                ),
                DEFAULT_SETTINGS["runtime"]["local"]["allow_oversubscribe"],
            ),
            "oversubscribe_factor": normalize_positive_float_setting(
                local_runtime.get(
                    "oversubscribe_factor",
                    DEFAULT_SETTINGS["runtime"]["local"]["oversubscribe_factor"],
                ),
                "Local oversubscribe factor",
                DEFAULT_SETTINGS["runtime"]["local"]["oversubscribe_factor"],
                minimum=1.0,
            ),
            "allow_live_load_bypass": normalize_bool_setting(
                local_runtime.get(
                    "allow_live_load_bypass",
                    DEFAULT_SETTINGS["runtime"]["local"]["allow_live_load_bypass"],
                ),
                DEFAULT_SETTINGS["runtime"]["local"]["allow_live_load_bypass"],
            ),
            "live_cpu_target_factor": normalize_positive_float_setting(
                local_runtime.get(
                    "live_cpu_target_factor",
                    DEFAULT_SETTINGS["runtime"]["local"]["live_cpu_target_factor"],
                ),
                "Local live CPU target factor",
                DEFAULT_SETTINGS["runtime"]["local"]["live_cpu_target_factor"],
                minimum=0.1,
            ),
            "live_min_free_ram_mb": normalize_positive_int_setting(
                local_runtime.get(
                    "live_min_free_ram_mb",
                    DEFAULT_SETTINGS["runtime"]["local"]["live_min_free_ram_mb"],
                ),
                "Local live minimum free RAM (MB)",
                DEFAULT_SETTINGS["runtime"]["local"]["live_min_free_ram_mb"],
                minimum=1,
            ),
        },
        "slurm": {
            "orca_base": normalize_local_directory_setting(
                slurm_runtime.get("orca_base", ""),
                "SLURM ORCA path",
            ),
            "submit_templates_dir": normalize_local_directory_setting(
                slurm_runtime.get("submit_templates_dir", ""),
                "SLURM submit templates path",
            ),
            "profile": str(slurm_runtime.get("profile", "") or "").strip(),
        },
    }
    features = normalized.get("features", {})
    if features is None:
        features = {}
    if not isinstance(features, dict):
        raise ValueError("Settings key 'features' must be a JSON object.")
    normalized_features = dict(features)
    if "remote_archive_enabled" in normalized_features:
        normalized_features["remote_archive_enabled"] = bool(normalized_features["remote_archive_enabled"])
    normalized["features"] = normalized_features
    ui = normalized.get("ui", {})
    if ui is None:
        ui = {}
    if not isinstance(ui, dict):
        raise ValueError("Settings key 'ui' must be a JSON object.")
    tabs = ui.get("tabs", {})
    if tabs is None:
        tabs = {}
    if not isinstance(tabs, dict):
        raise ValueError("Settings key 'ui.tabs' must be a JSON object.")
    default_ui = DEFAULT_SETTINGS.get("ui", {}) or {}
    default_tabs = default_ui.get("tabs", {}) or {}
    normalized_ui = dict(ui)
    normalized_ui["tabs"] = {
        "order": normalize_string_list_setting(
            tabs.get("order", default_tabs.get("order", [])),
            "Tab order",
        ),
        "hidden": normalize_string_list_setting(
            tabs.get("hidden", default_tabs.get("hidden", [])),
            "Hidden tabs",
        ),
    }
    normalized["ui"] = normalized_ui
    return normalized


def _merge_missing_defaults(payload, defaults):
    merged = dict(payload)
    for key, default_value in dict(defaults or {}).items():
        if key not in merged:
            merged[key] = default_value if not isinstance(default_value, dict) else dict(default_value)
            continue
        if isinstance(default_value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_missing_defaults(merged.get(key) or {}, default_value)
    return merged


def _legacy_transfer_path(settings_path):
    return Path(settings_path).with_name(LEGACY_TRANSFER_CONFIG_NAME)


def load_settings(settings_path=None):
    path = get_settings_path(settings_path)
    if path.exists():
        normalized = _normalized_settings_dict(_read_json(path))
        merged = _merge_missing_defaults(normalized, DEFAULT_SETTINGS)
        if merged != normalized:
            _write_json_atomic(path, merged)
        return merged

    legacy_path = _legacy_transfer_path(path)
    if legacy_path.exists():
        legacy = _read_json(legacy_path)
        host, user, remote_path, port = normalize_ssh_transfer_settings(
            legacy.get("host", ""),
            legacy.get("user", ""),
            legacy.get("remote_path", ""),
            legacy.get("port", 22),
        )
        migrated = {
            "transfer": {
                "host": host,
                "user": user,
                "remote_path": remote_path,
                "port": port,
                "updated_at": utc_now_iso(),
            }
        }
        save_settings(migrated, path)
        return _merge_missing_defaults(migrated, DEFAULT_SETTINGS)

    return _merge_missing_defaults({}, DEFAULT_SETTINGS)


def save_settings(settings, settings_path=None):
    path = get_settings_path(settings_path)
    normalized = _merge_missing_defaults(_normalized_settings_dict(settings), DEFAULT_SETTINGS)
    _write_json_atomic(path, normalized)
    return normalized


def load_transfer_settings(settings_path=None):
    path = get_settings_path(settings_path)
    settings = load_settings(path)
    payload = settings.get("transfer", {}) or {}
    if not payload:
        return None
    host, user, remote_path, port = normalize_ssh_transfer_settings(
        payload.get("host", ""),
        payload.get("user", ""),
        payload.get("remote_path", ""),
        payload.get("port", 22),
    )
    return {
        "host": host,
        "user": user,
        "remote_path": remote_path,
        "port": port,
        "settings_path": str(path),
    }


def load_runtime_settings(settings_path=None):
    path = get_settings_path(settings_path)
    settings = load_settings(path)
    return settings.get("runtime", {}) or _merge_missing_defaults({}, DEFAULT_SETTINGS["runtime"])


def save_transfer_settings(host, user, remote_path, port, settings_path=None):
    path = get_settings_path(settings_path)
    host, user, remote_path, port = normalize_ssh_transfer_settings(host, user, remote_path, port)
    settings = load_settings(path)
    settings["transfer"] = {
        "host": host,
        "user": user,
        "remote_path": remote_path,
        "port": port,
        "updated_at": utc_now_iso(),
    }
    save_settings(settings, path)
    return {
        "host": host,
        "user": user,
        "remote_path": remote_path,
        "port": port,
        "settings_path": str(path),
    }


def load_remote_archive_enabled(settings_path=None):
    path = get_settings_path(settings_path)
    try:
        settings = load_settings(path)
    except Exception:
        return bool(DEFAULT_SETTINGS.get("features", {}).get("remote_archive_enabled", False))
    features = settings.get("features", {}) or {}
    return bool(features.get("remote_archive_enabled", DEFAULT_SETTINGS["features"]["remote_archive_enabled"]))


def save_remote_archive_enabled(enabled, settings_path=None):
    path = get_settings_path(settings_path)
    settings = load_settings(path)
    features = settings.get("features", {}) or {}
    features["remote_archive_enabled"] = bool(enabled)
    settings["features"] = features
    save_settings(settings, path)
    return bool(features["remote_archive_enabled"])
