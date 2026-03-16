"""Runtime configuration helpers for DELFIN dashboard and install scripts."""

from __future__ import annotations

import os
import shutil
from contextlib import contextmanager
from pathlib import Path

from delfin.qm_runtime import check_tools


def normalize_runtime_path(path_value: str | Path | None) -> str:
    text = str(path_value or "").strip()
    if not text:
        return ""
    return str(Path(text).expanduser())


def resolve_backend_choice(
    explicit_backend: str | None,
    configured_backend: str | None,
    *,
    slurm_available: bool,
) -> str:
    requested = str(explicit_backend or "").strip().lower()
    if requested and requested != "auto":
        return requested

    configured = str(configured_backend or "").strip().lower()
    if configured and configured != "auto":
        return configured

    return "slurm" if slurm_available else "local"


def resolve_orca_base(
    explicit_orca_base: str | None,
    runtime_settings: dict | None,
    backend: str,
    *,
    auto_candidates: list[str] | None = None,
    local_default: str = "/opt/orca",
) -> str:
    runtime_settings = runtime_settings or {}
    local_settings = runtime_settings.get("local", {}) or {}
    slurm_settings = runtime_settings.get("slurm", {}) or {}

    candidates = []
    if explicit_orca_base:
        candidates.append(explicit_orca_base)
    if backend == "local":
        candidates.append(local_settings.get("orca_base", ""))
    if backend == "slurm":
        candidates.append(slurm_settings.get("orca_base", ""))
    candidates.append(runtime_settings.get("orca_base", ""))

    for env_key in ("DELFIN_ORCA_BASE", "ORCA_BINARY", "ORCA_PATH"):
        env_value = os.environ.get(env_key)
        if env_value:
            candidates.append(env_value)

    candidates.extend(auto_candidates or [])

    which_orca = shutil.which("orca")
    if which_orca:
        candidates.append(str(Path(which_orca).resolve().parent))

    if backend == "local":
        candidates.append(local_default)

    for candidate in candidates:
        normalized = normalize_runtime_path(candidate)
        if not normalized:
            continue
        path = Path(normalized)
        if path.is_file():
            return str(path.parent)
        return str(path)

    return ""


def resolve_submit_templates_dir(runtime_settings: dict | None, fallback_dir: str | Path) -> Path:
    runtime_settings = runtime_settings or {}
    slurm_settings = runtime_settings.get("slurm", {}) or {}
    configured = normalize_runtime_path(slurm_settings.get("submit_templates_dir", ""))
    return Path(configured) if configured else Path(fallback_dir)


def get_packaged_submit_templates_dir() -> Path:
    return (Path(__file__).resolve().parent / "submit_templates").resolve()


def _prepend_path_env_once(path_entry: str) -> None:
    normalized = str(Path(path_entry).expanduser())
    current_entries = [item for item in os.environ.get("PATH", "").split(os.pathsep) if item]
    if normalized in current_entries:
        return
    os.environ["PATH"] = (
        normalized
        if not current_entries
        else f"{normalized}{os.pathsep}{os.environ.get('PATH', '')}"
    )


def apply_runtime_environment(*, qm_tools_root: str = "", orca_base: str = "") -> None:
    qm_root = normalize_runtime_path(qm_tools_root)
    if qm_root:
        os.environ["DELFIN_QM_ROOT"] = qm_root
        os.environ["DELFIN_QM_TOOLS_ROOT"] = qm_root
        bin_dir = Path(qm_root) / "bin"
        if bin_dir.exists():
            _prepend_path_env_once(str(bin_dir))
        xtb4stda_share = Path(qm_root) / "share" / "xtb4stda"
        if xtb4stda_share.exists():
            os.environ.setdefault("XTB4STDAHOME", str(xtb4stda_share))
        os.environ.setdefault("STD2HOME", qm_root)

    orca_root = normalize_runtime_path(orca_base)
    if orca_root:
        os.environ["DELFIN_ORCA_BASE"] = orca_root
        os.environ["ORCA_PATH"] = orca_root
        _prepend_path_env_once(orca_root)
        orca_plot = Path(orca_root) / "orca_plot"
        if orca_plot.exists():
            os.environ["ORCA_PLOT"] = str(orca_plot)


@contextmanager
def temporary_environment(updates: dict[str, str]):
    previous = {key: os.environ.get(key) for key in updates}
    try:
        for key, value in updates.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = str(value)
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def collect_runtime_diagnostics(
    runtime_settings: dict | None,
    *,
    backend: str,
    effective_orca_base: str = "",
    submit_templates_dir: str | Path | None = None,
) -> list[dict[str, str]]:
    runtime_settings = runtime_settings or {}
    diagnostics: list[dict[str, str]] = [
        {"name": "backend", "status": "ok", "detail": backend}
    ]

    orca_root = normalize_runtime_path(effective_orca_base)
    orca_binary = ""
    if orca_root:
        orca_candidate = Path(orca_root) / "orca"
        if orca_candidate.exists():
            orca_binary = str(orca_candidate)
    if not orca_binary:
        which_orca = shutil.which("orca")
        if which_orca:
            orca_binary = which_orca

    diagnostics.append(
        {
            "name": "orca",
            "status": "ok" if orca_binary else "missing",
            "detail": orca_binary or (orca_root or "not found"),
        }
    )

    if backend == "slurm":
        submit_dir = Path(submit_templates_dir) if submit_templates_dir else None
        sbatch_path = shutil.which("sbatch")
        diagnostics.append(
            {
                "name": "sbatch",
                "status": "ok" if sbatch_path else "missing",
                "detail": sbatch_path or "sbatch not found in PATH",
            }
        )
        if submit_dir:
            diagnostics.append(
                {
                    "name": "slurm-templates",
                    "status": (
                        "ok"
                        if (submit_dir / "submit_delfin.sh").exists()
                        and (submit_dir / "submit_turbomole.sh").exists()
                        else "missing"
                    ),
                    "detail": str(submit_dir),
                }
            )
    else:
        diagnostics.append(
            {
                "name": "local-runner",
                "status": "ok",
                "detail": "bundled Python local runner",
            }
        )

    qm_tools_root = normalize_runtime_path(runtime_settings.get("qm_tools_root", ""))
    env_updates: dict[str, str] = {}
    if qm_tools_root:
        env_updates["DELFIN_QM_ROOT"] = qm_tools_root
        env_updates["DELFIN_QM_TOOLS_ROOT"] = qm_tools_root
    if orca_root:
        env_updates["ORCA_PATH"] = orca_root
        env_updates["DELFIN_ORCA_BASE"] = orca_root

    with temporary_environment(env_updates):
        qm_results = check_tools(["xtb", "crest", "std2", "stda", "xtb4stda", "dftb+"])
    for name, resolved in qm_results:
        diagnostics.append(
            {
                "name": name,
                "status": "ok" if resolved else "missing",
                "detail": resolved.path if resolved else "not found",
            }
        )

    if qm_tools_root:
        diagnostics.append(
            {
                "name": "qm_tools_root",
                "status": "ok" if Path(qm_tools_root).exists() else "missing",
                "detail": qm_tools_root,
            }
        )

    return diagnostics
