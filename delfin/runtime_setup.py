"""Runtime configuration helpers for DELFIN dashboard and install scripts."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tarfile
from contextlib import contextmanager
from pathlib import Path

from delfin.qm_runtime import check_tools

try:
    import psutil  # type: ignore
except ImportError:
    psutil = None  # type: ignore


_ORCA_VERSION_RE = re.compile(r"orca[_-]?(\d+(?:[_\.-]\d+)*)", re.IGNORECASE)


def normalize_runtime_path(path_value: str | Path | None) -> str:
    text = str(path_value or "").strip()
    if not text:
        return ""
    return str(Path(text).expanduser())


def detect_local_runtime_limits() -> tuple[int, int]:
    cpu_count = os.cpu_count() or 1
    memory_mb = 1_400_000
    if psutil is not None:
        try:
            memory_mb = max(1, int(psutil.virtual_memory().total // (1024 * 1024)))
        except Exception:
            memory_mb = 1_400_000
    return max(1, int(cpu_count)), memory_mb


def normalize_orca_base(path_value: str | Path | None) -> str:
    normalized = normalize_runtime_path(path_value)
    if not normalized:
        return ""
    path = Path(normalized)
    if path.is_file():
        if path.name.lower().startswith("orca"):
            return str(path.parent.resolve())
        return str(path.resolve().parent)
    candidate = path / "orca"
    if candidate.is_file():
        return str(path.resolve())
    return str(path.resolve())


def describe_orca_installation(path_value: str | Path | None) -> str:
    normalized = normalize_orca_base(path_value)
    if not normalized:
        return "ORCA"
    name = Path(normalized).name
    match = _ORCA_VERSION_RE.search(name)
    if match:
        version = match.group(1).replace("_", ".").replace("-", ".")
        return f"ORCA {version}"
    return f"ORCA ({name})"


def discover_orca_installations(
    *,
    seed_candidates: list[str] | None = None,
    search_roots: list[str | Path] | None = None,
) -> list[str]:
    candidates: list[str] = []
    seen: set[str] = set()

    def _add_candidate(raw_value: str | Path | None) -> None:
        normalized = normalize_orca_base(raw_value)
        if not normalized:
            return
        resolved = str(Path(normalized).expanduser().resolve())
        if resolved in seen:
            return
        orca_binary = Path(resolved) / "orca"
        if not orca_binary.is_file():
            return
        seen.add(resolved)
        candidates.append(resolved)

    for env_key in ("DELFIN_ORCA_BASE", "ORCA_BINARY", "ORCA_PATH"):
        _add_candidate(os.environ.get(env_key))

    _add_candidate(shutil.which("orca"))

    for candidate in seed_candidates or []:
        _add_candidate(candidate)

    roots = []
    for root in search_roots or []:
        normalized = normalize_runtime_path(root)
        if normalized:
            roots.append(Path(normalized))

    for root in roots:
        if not root.exists():
            continue
        try:
            entries = list(root.iterdir())
        except Exception:
            continue
        for entry in entries:
            if not entry.is_dir():
                continue
            if "orca" not in entry.name.lower():
                continue
            _add_candidate(entry)

    return sorted(
        candidates,
        key=lambda value: (
            "6.1.1" not in describe_orca_installation(value),
            describe_orca_installation(value),
            value,
        ),
    )


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
        normalized = normalize_orca_base(candidate)
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


def get_packaged_qm_tools_dir() -> Path:
    return (Path(__file__).resolve().parent / "qm_tools").resolve()


def get_packaged_installers_dir() -> Path:
    return (Path(__file__).resolve().parent / "installers").resolve()


def get_user_qm_tools_dir(target_dir: str | Path | None = None) -> Path:
    if target_dir:
        return Path(target_dir).expanduser()
    return (Path.home() / ".delfin" / "qm_tools").expanduser()


def get_repo_submit_templates_dir(repo_dir: str | Path | None) -> Path | None:
    if not repo_dir:
        return None
    repo_path = Path(repo_dir).expanduser()
    candidate = (
        repo_path
        / "examples"
        / "example_Job_Submission_Scripts"
        / "BwUniCluster"
        / "submit_sh"
    )
    return candidate if candidate.is_dir() else None


def get_repo_qm_tools_dir(repo_dir: str | Path | None) -> Path | None:
    if not repo_dir:
        return None
    repo_path = Path(repo_dir).expanduser()
    candidate = repo_path / "delfin" / "qm_tools"
    return candidate if candidate.is_dir() else None


def get_repo_bwunicluster_install_script(repo_dir: str | Path | None) -> Path | None:
    if not repo_dir:
        return None
    repo_path = Path(repo_dir).expanduser()
    candidate = repo_path / "scripts" / "install_delfin_bwu.sh"
    return candidate if candidate.is_file() else None


def get_packaged_bwunicluster_install_script() -> Path | None:
    candidate = get_packaged_installers_dir() / "install_delfin_bwu.sh"
    return candidate if candidate.is_file() else None


def write_delfin_env_file(
    *,
    repo_dir: str | Path | None = None,
    orca_base: str = "",
    qm_tools_root: str = "",
    env_path: str | Path | None = None,
) -> Path:
    target = Path(env_path).expanduser() if env_path else (Path.home() / ".delfin_env.sh")
    repo_path = Path(repo_dir).expanduser() if repo_dir else None
    orca_root = normalize_orca_base(orca_base)
    qm_root = normalize_runtime_path(qm_tools_root)

    lines = [
        "# Auto-generated by DELFIN Runtime Setup",
    ]
    if repo_path:
        lines.append(f'export DELFIN_REPO="{repo_path}"')
    if orca_root:
        lines.extend(
            [
                f'export ORCA_DIR="{orca_root}"',
                f'export ORCA_BASE="{orca_root}"',
                f'export ORCA_PATH="{orca_root}"',
            ]
        )
        orca_plot = Path(orca_root) / "orca_plot"
        if orca_plot.exists():
            lines.append(f'export ORCA_PLOT="{orca_plot}"')
    if qm_root:
        lines.extend(
            [
                f'export DELFIN_QM_ROOT="{qm_root}"',
                f'export DELFIN_QM_TOOLS_ROOT="{qm_root}"',
                f'export STD2HOME="{qm_root}"',
            ]
        )
        xtb4stda_share = Path(qm_root) / "share" / "xtb4stda"
        if xtb4stda_share.exists():
            lines.append(f'export XTB4STDAHOME="{xtb4stda_share}"')

    path_entries = []
    if qm_root and (Path(qm_root) / "bin").is_dir():
        path_entries.append(str(Path(qm_root) / "bin"))
    if orca_root:
        path_entries.append(orca_root)
    if path_entries:
        lines.append(f'export PATH="{":".join(path_entries)}:$PATH"')

    if repo_path and (repo_path / ".venv").is_dir():
        lines.extend(
            [
                'export DELFIN_AUTO_ACTIVATE_VENV="${DELFIN_AUTO_ACTIVATE_VENV:-1}"',
                'if [ "${DELFIN_AUTO_ACTIVATE_VENV}" = "1" ] '
                '&& [ -z "${VIRTUAL_ENV:-}" ] '
                f'&& [ -d "{repo_path / ".venv"}" ] '
                '&& [ -z "${VSCODE_PID:-}" ] '
                '&& [ "${TERM_PROGRAM:-}" != "vscode" ]; then',
                f'  source "{repo_path / ".venv" / "bin" / "activate"}"',
                'fi',
            ]
        )

    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    try:
        target.chmod(0o600)
    except Exception:
        pass
    return target


def ensure_shell_sources_delfin_env(
    shell_rc_paths: list[str | Path] | None = None,
    *,
    env_path: str | Path | None = None,
) -> list[Path]:
    target_env = Path(env_path).expanduser() if env_path else (Path.home() / ".delfin_env.sh")
    source_line = f'source "{target_env}"'
    updated: list[Path] = []
    for candidate in shell_rc_paths or [Path.home() / ".bashrc", Path.home() / ".profile"]:
        path = Path(candidate).expanduser()
        existing = path.read_text(encoding="utf-8") if path.exists() else ""
        if source_line not in existing:
            if existing and not existing.endswith("\n"):
                existing += "\n"
            existing += f"\n# DELFIN environment\n{source_line}\n"
            path.write_text(existing, encoding="utf-8")
            updated.append(path)
        elif path.exists():
            updated.append(path)
    return updated


def package_repo_venv(repo_dir: str | Path | None) -> tuple[Path | None, Path | None]:
    if not repo_dir:
        return None, None
    repo_path = Path(repo_dir).expanduser()
    venv_dir = repo_path / ".venv"
    if not venv_dir.is_dir():
        return None, None

    tar_path = repo_path / "delfin_venv.tar"
    if tar_path.exists():
        tar_path.unlink()
    with tarfile.open(tar_path, "w") as tar_handle:
        tar_handle.add(venv_dir, arcname=".venv")

    cache_dir = repo_path / ".runtime_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return tar_path, cache_dir


def run_bwunicluster_installer(
    *,
    repo_dir: str | Path | None,
    orca_base: str = "",
    calc_dir: str | Path | None = None,
    archive_dir: str | Path | None = None,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    installer = get_repo_bwunicluster_install_script(repo_dir)
    repo_path = Path(repo_dir).expanduser() if repo_dir else None
    if installer is None:
        installer = get_packaged_bwunicluster_install_script()
    if installer is None:
        raise FileNotFoundError("BwUniCluster installer script not found in repo or packaged resources.")

    env = os.environ.copy()
    if repo_path is not None:
        env["DELFIN_REPO"] = str(repo_path)
    if calc_dir:
        env["DELFIN_CALC_DIR"] = str(Path(calc_dir).expanduser())
    if archive_dir:
        env["DELFIN_ARCHIVE_DIR"] = str(Path(archive_dir).expanduser())
    normalized_orca = normalize_orca_base(orca_base)
    if normalized_orca:
        env["ORCA_DIR"] = normalized_orca
    if extra_env:
        env.update({str(key): str(value) for key, value in extra_env.items()})

    shell_command = (
        'if ! command -v module >/dev/null 2>&1; then '
        'for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do '
        '[ -f "$init" ] && source "$init" && break; '
        'done; '
        'fi; '
        f'bash "{installer}"'
    )
    return subprocess.run(
        ["bash", "-lc", shell_command],
        cwd=str(repo_path or Path.home()),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )


def collect_bwunicluster_verification(
    *,
    repo_dir: str | Path | None = None,
    orca_base: str = "",
    calc_dir: str | Path | None = None,
    archive_dir: str | Path | None = None,
) -> list[dict[str, str]]:
    repo_path = Path(repo_dir).expanduser() if repo_dir else None
    calc_path = Path(calc_dir).expanduser() if calc_dir else (Path.home() / "calc")
    archive_path = Path(archive_dir).expanduser() if archive_dir else (Path.home() / "archive")
    openmpi_path = Path.home() / "software" / "openmpi-4.1.8" / "bin" / "mpirun"
    submit_dir = get_repo_submit_templates_dir(repo_path) or get_packaged_submit_templates_dir()
    install_script = get_repo_bwunicluster_install_script(repo_path) or get_packaged_bwunicluster_install_script()

    effective_orca = normalize_orca_base(orca_base)
    if not effective_orca:
        discovered = discover_orca_installations(
            search_roots=[Path.home() / "software", Path.home() / "apps", Path.home() / "local", Path("/opt")]
        )
        effective_orca = discovered[0] if discovered else ""

    return [
        {
            "name": "install-script",
            "status": "ok" if install_script and Path(install_script).is_file() else "missing",
            "detail": str(install_script) if install_script else "not found",
        },
        {
            "name": "repo",
            "status": "ok" if (repo_path and repo_path.exists()) else "missing",
            "detail": str(repo_path) if repo_path else "no repo checkout configured",
        },
        {
            "name": "venv",
            "status": "ok" if (repo_path and (repo_path / ".venv").is_dir()) else "missing",
            "detail": str((repo_path / ".venv") if repo_path else Path.home() / "software" / "delfin" / ".venv"),
        },
        {
            "name": "openmpi",
            "status": "ok" if openmpi_path.is_file() else "missing",
            "detail": str(openmpi_path),
        },
        {
            "name": "orca",
            "status": "ok" if effective_orca and (Path(effective_orca) / "orca").is_file() else "missing",
            "detail": str(Path(effective_orca) / "orca") if effective_orca else "not found",
        },
        {
            "name": "orca_plot",
            "status": "ok" if effective_orca and (Path(effective_orca) / "orca_plot").is_file() else "missing",
            "detail": str(Path(effective_orca) / "orca_plot") if effective_orca else "not found",
        },
        {
            "name": "calc",
            "status": "ok" if calc_path.exists() else "missing",
            "detail": str(calc_path),
        },
        {
            "name": "archive",
            "status": "ok" if archive_path.exists() else "missing",
            "detail": str(archive_path),
        },
        {
            "name": "submit-templates",
            "status": "ok" if (submit_dir / "submit_delfin.sh").is_file() else "missing",
            "detail": str(submit_dir),
        },
        {
            "name": "sbatch",
            "status": "ok" if shutil.which("sbatch") else "missing",
            "detail": shutil.which("sbatch") or "sbatch not found in PATH",
        },
    ]


def prepare_bwunicluster_user_setup(
    *,
    repo_dir: str | Path | None = None,
    calc_dir: str | Path | None = None,
    archive_dir: str | Path | None = None,
    orca_base: str = "",
    qm_tools_root: str = "",
    install_qm_tools: bool = True,
) -> dict[str, object]:
    detected_cores, detected_ram_mb = detect_local_runtime_limits()
    effective_calc_dir = Path(calc_dir).expanduser() if calc_dir else (Path.home() / "calc")
    effective_archive_dir = Path(archive_dir).expanduser() if archive_dir else (Path.home() / "archive")
    effective_calc_dir.mkdir(parents=True, exist_ok=True)
    effective_archive_dir.mkdir(parents=True, exist_ok=True)

    repo_path = Path(repo_dir).expanduser() if repo_dir else None
    submit_templates_dir = get_repo_submit_templates_dir(repo_path) or get_packaged_submit_templates_dir()

    effective_orca = normalize_orca_base(orca_base)
    if not effective_orca:
        discovered = discover_orca_installations(
            seed_candidates=[],
            search_roots=[Path.home() / "software", Path.home() / "apps", Path("/opt")],
        )
        effective_orca = discovered[0] if discovered else ""

    qm_root = normalize_runtime_path(qm_tools_root)
    installer_result = None
    if install_qm_tools:
        target, installer_result = run_qm_tools_installer(qm_root or None)
        qm_root = str(target)
    elif not qm_root:
        repo_qm_root = get_repo_qm_tools_dir(repo_path)
        if repo_qm_root:
            qm_root = str(repo_qm_root)
        else:
            qm_root = str(stage_packaged_qm_tools())

    tar_path, cache_dir = package_repo_venv(repo_path)
    env_file = write_delfin_env_file(
        repo_dir=repo_path,
        orca_base=effective_orca,
        qm_tools_root=qm_root,
    )
    shell_files = ensure_shell_sources_delfin_env(env_path=env_file)

    runtime_payload = {
        "backend": "slurm",
        "orca_base": effective_orca,
        "qm_tools_root": qm_root,
        "local": {
            "orca_base": effective_orca,
            "max_cores": detected_cores,
            "max_ram_mb": detected_ram_mb,
        },
        "slurm": {
            "orca_base": effective_orca,
            "submit_templates_dir": str(submit_templates_dir),
            "profile": "bwunicluster3",
        },
    }

    return {
        "runtime": runtime_payload,
        "paths": {
            "calculations_dir": str(effective_calc_dir),
            "archive_dir": str(effective_archive_dir),
        },
        "env_file": str(env_file),
        "shell_files": [str(path) for path in shell_files],
        "submit_templates_dir": str(submit_templates_dir),
        "orca_base": effective_orca,
        "qm_tools_root": qm_root,
        "venv_tarball": str(tar_path) if tar_path else "",
        "runtime_cache_dir": str(cache_dir) if cache_dir else "",
        "qm_tools_installer_output": (installer_result.stdout if installer_result else ""),
        "qm_tools_installer_returncode": (installer_result.returncode if installer_result else 0),
    }


def stage_packaged_qm_tools(target_dir: str | Path | None = None) -> Path:
    source = get_packaged_qm_tools_dir()
    if not source.is_dir():
        raise FileNotFoundError(f"Packaged qm_tools directory not found: {source}")

    target = get_user_qm_tools_dir(target_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        source,
        target,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo", "downloads"),
    )
    return target


def run_qm_tools_installer(
    target_dir: str | Path | None = None,
    *,
    extra_env: dict[str, str] | None = None,
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    target = stage_packaged_qm_tools(target_dir)
    installer = target / "install_qm_tools.sh"
    if not installer.is_file():
        raise FileNotFoundError(f"qm_tools installer not found: {installer}")

    env = os.environ.copy()
    env["DELFIN_QM_ROOT"] = str(target)
    env["DELFIN_QM_TOOLS_ROOT"] = str(target)
    if extra_env:
        env.update({str(key): str(value) for key, value in extra_env.items()})

    result = subprocess.run(
        ["bash", str(installer)],
        cwd=str(target),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return target, result


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
