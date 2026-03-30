import importlib.util
import sys
import types
from pathlib import Path


_ROOT = Path(__file__).resolve().parents[1] / "delfin"
_MODULE_PATH = _ROOT / "dashboard" / "backend_local.py"

if "delfin" not in sys.modules:
    _PKG = types.ModuleType("delfin")
    _PKG.__path__ = [str(_ROOT)]
    _PKG.__version__ = "test"
    sys.modules["delfin"] = _PKG

if "delfin.dashboard" not in sys.modules:
    _DASHBOARD_PKG = types.ModuleType("delfin.dashboard")
    _DASHBOARD_PKG.__path__ = [str(_ROOT / "dashboard")]
    sys.modules["delfin.dashboard"] = _DASHBOARD_PKG

if "delfin.dashboard.helpers" not in sys.modules:
    _HELPERS = types.ModuleType("delfin.dashboard.helpers")

    def _parse_time_to_seconds(value):
        parts = str(value).split(":")
        if len(parts) != 3:
            raise ValueError(f"Unsupported time format: {value}")
        hours, minutes, seconds = (int(part) for part in parts)
        return hours * 3600 + minutes * 60 + seconds

    _HELPERS.parse_time_to_seconds = _parse_time_to_seconds
    sys.modules["delfin.dashboard.helpers"] = _HELPERS

_SPEC = importlib.util.spec_from_file_location("delfin.dashboard.backend_local", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load backend_local module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)


def test_local_backend_falls_back_to_python_runner_when_script_is_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(_MODULE.threading.Thread, "start", lambda self: None)

    backend = _MODULE.LocalJobBackend(run_script=tmp_path / "missing.sh")

    assert backend._has_launch_target() is True
    assert backend._build_launch_command(3600) == [
        "timeout",
        "3600",
        sys.executable,
        "-m",
        "delfin.dashboard.local_runner",
    ]


def test_local_backend_prefers_explicit_run_script_when_present(monkeypatch, tmp_path):
    monkeypatch.setattr(_MODULE.threading.Thread, "start", lambda self: None)
    run_script = tmp_path / "run_local.sh"
    run_script.write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")

    backend = _MODULE.LocalJobBackend(run_script=run_script)

    assert backend._build_launch_command(60) == [
        "timeout",
        "60",
        "bash",
        str(run_script),
    ]


def test_local_backend_can_oversubscribe_cores_when_enabled(monkeypatch, tmp_path):
    monkeypatch.setattr(_MODULE.threading.Thread, "start", lambda self: None)
    backend = _MODULE.LocalJobBackend(
        run_script=tmp_path / "missing.sh",
        max_cores=12,
        max_ram_mb=1_400_000,
        allow_oversubscribe=True,
        oversubscribe_factor=2.0,
    )

    jobs_data = {
        "jobs": [
            {"job_id": 1, "status": "RUNNING", "pal": 8, "maxcore": 1000},
            {"job_id": 2, "status": "PENDING", "pal": 8, "maxcore": 1000, "job_dir": str(tmp_path)},
        ]
    }

    monkeypatch.setattr(backend, "_load_jobs", lambda: jobs_data)
    monkeypatch.setattr(backend, "_save_jobs", lambda data: None)
    monkeypatch.setattr(backend, "_update_job_status", lambda job: job)
    started = []
    monkeypatch.setattr(backend, "_start_job", lambda job, data: started.append(job["job_id"]) or True)

    changed = backend._try_start_pending_jobs()

    assert changed is True
    assert started == [2]


def test_local_backend_can_oversubscribe_ram_when_enabled(monkeypatch, tmp_path):
    monkeypatch.setattr(_MODULE.threading.Thread, "start", lambda self: None)
    backend = _MODULE.LocalJobBackend(
        run_script=tmp_path / "missing.sh",
        max_cores=384,
        max_ram_mb=1_000,
        allow_oversubscribe=True,
        oversubscribe_factor=1.5,
    )

    jobs_data = {
        "jobs": [
            {"job_id": 1, "status": "RUNNING", "pal": 1, "maxcore": 900},
            {"job_id": 2, "status": "PENDING", "pal": 1, "maxcore": 500, "job_dir": str(tmp_path)},
        ]
    }

    monkeypatch.setattr(backend, "_load_jobs", lambda: jobs_data)
    monkeypatch.setattr(backend, "_save_jobs", lambda data: None)
    monkeypatch.setattr(backend, "_update_job_status", lambda job: job)
    started = []
    monkeypatch.setattr(backend, "_start_job", lambda job, data: started.append(job["job_id"]) or True)

    changed = backend._try_start_pending_jobs()

    assert changed is True
    assert started == [2]


def test_local_backend_oversubscribe_budget_is_integer(monkeypatch, tmp_path):
    monkeypatch.setattr(_MODULE.threading.Thread, "start", lambda self: None)
    backend = _MODULE.LocalJobBackend(
        run_script=tmp_path / "missing.sh",
        max_cores=12,
        max_ram_mb=1_400_000,
        allow_oversubscribe=True,
        oversubscribe_factor=1.1,
    )

    assert backend._core_budget() == 13
    assert backend._ram_budget() == 1_540_000
