"""Runtime configuration helpers for DELFIN dashboard and install scripts."""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import tarfile
from contextlib import contextmanager
from pathlib import Path

from delfin.qm_runtime import check_tools, get_csp_tools_root, get_mlp_tools_root

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
    memory_mb = 8_192
    if psutil is not None:
        try:
            memory_mb = max(1, int(psutil.virtual_memory().total // (1024 * 1024)))
        except Exception:
            memory_mb = 8_192
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

    def _version_sort_key(value: str) -> tuple:
        """Sort ORCA installations by version number, highest first."""
        desc = describe_orca_installation(value)
        match = _ORCA_VERSION_RE.search(desc)
        if match:
            parts = match.group(1).replace("_", ".").replace("-", ".").split(".")
            # Pad to 4 parts, convert to ints for proper numeric sorting
            nums = []
            for p in parts:
                try:
                    nums.append(int(p))
                except ValueError:
                    nums.append(0)
            while len(nums) < 4:
                nums.append(0)
            # Negate for descending order (highest version first)
            return tuple(-n for n in nums)
        # No version found — sort to the end
        return (0, 0, 0, 0)

    return sorted(candidates, key=_version_sort_key)


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


def get_packaged_csp_tools_dir() -> Path:
    return (Path(__file__).resolve().parent / "csp_tools").resolve()


def get_packaged_mlp_tools_dir() -> Path:
    return (Path(__file__).resolve().parent / "mlp_tools").resolve()


def get_packaged_analysis_tools_dir() -> Path:
    return (Path(__file__).resolve().parent / "analysis_tools").resolve()


def get_packaged_ai_tools_dir() -> Path:
    return (Path(__file__).resolve().parent / "ai_tools").resolve()


def get_packaged_installers_dir() -> Path:
    return (Path(__file__).resolve().parent / "installers").resolve()


def get_user_qm_tools_dir(target_dir: str | Path | None = None) -> Path:
    if target_dir:
        return Path(target_dir).expanduser()
    return (Path.home() / ".delfin" / "qm_tools").expanduser()


def get_user_csp_tools_dir(target_dir: str | Path | None = None) -> Path:
    if target_dir:
        return Path(target_dir).expanduser()
    return (Path.home() / ".delfin" / "csp_tools").expanduser()


def get_user_mlp_tools_dir(target_dir: str | Path | None = None) -> Path:
    if target_dir:
        return Path(target_dir).expanduser()
    return (Path.home() / ".delfin" / "mlp_tools").expanduser()


def get_user_analysis_tools_dir(target_dir: str | Path | None = None) -> Path:
    if target_dir:
        return Path(target_dir).expanduser()
    return (Path.home() / ".delfin" / "analysis_tools").expanduser()


def get_user_ai_tools_dir(target_dir: str | Path | None = None) -> Path:
    if target_dir:
        return Path(target_dir).expanduser()
    return (Path.home() / ".delfin" / "ai_tools").expanduser()


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

    def _sq(value: str | Path) -> str:
        """Shell-quote a value for safe embedding in a shell script."""
        return shlex.quote(str(value))

    lines = [
        "# Auto-generated by DELFIN Runtime Setup",
    ]
    if repo_path:
        lines.append(f'export DELFIN_REPO={_sq(repo_path)}')
    if orca_root:
        lines.extend(
            [
                f'export ORCA_DIR={_sq(orca_root)}',
                f'export ORCA_BASE={_sq(orca_root)}',
                f'export ORCA_PATH={_sq(orca_root)}',
            ]
        )
        orca_plot = Path(orca_root) / "orca_plot"
        if orca_plot.exists():
            lines.append(f'export ORCA_PLOT={_sq(orca_plot)}')
    if qm_root:
        lines.extend(
            [
                f'export DELFIN_QM_ROOT={_sq(qm_root)}',
                f'export DELFIN_QM_TOOLS_ROOT={_sq(qm_root)}',
                f'export STD2HOME={_sq(qm_root)}',
            ]
        )
        xtb4stda_share = Path(qm_root) / "share" / "xtb4stda"
        if xtb4stda_share.exists():
            lines.append(f'export XTB4STDAHOME={_sq(xtb4stda_share)}')

    path_entries = []
    if qm_root and (Path(qm_root) / "bin").is_dir():
        path_entries.append(str(Path(qm_root) / "bin"))
    if orca_root:
        path_entries.append(orca_root)
    if path_entries:
        joined = ":".join(path_entries)
        lines.append(f'export PATH={_sq(joined)}":$PATH"')

    if repo_path and (repo_path / ".venv").is_dir():
        venv_dir = _sq(repo_path / ".venv")
        activate = _sq(repo_path / ".venv" / "bin" / "activate")
        lines.extend(
            [
                'export DELFIN_AUTO_ACTIVATE_VENV="${DELFIN_AUTO_ACTIVATE_VENV:-1}"',
                'if [ "${DELFIN_AUTO_ACTIVATE_VENV}" = "1" ] '
                '&& [ -z "${VIRTUAL_ENV:-}" ] '
                f'&& [ -d {venv_dir} ] '
                '&& [ -z "${VSCODE_PID:-}" ] '
                '&& [ "${TERM_PROGRAM:-}" != "vscode" ]; then',
                f'  source {activate}',
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
        f'bash {shlex.quote(str(installer))}'
    )
    return subprocess.run(
        ["bash", "-lc", shell_command],
        cwd=str(repo_path or Path.home()),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=1800,
    )


def _find_openmpi_mpirun() -> Path:
    """Find mpirun in ~/software/openmpi-*/bin/ or fall back to PATH."""
    software_dir = Path.home() / "software"
    if software_dir.is_dir():
        candidates = sorted(software_dir.glob("openmpi-*/bin/mpirun"), reverse=True)
        for candidate in candidates:
            if candidate.is_file():
                return candidate
    which_mpirun = shutil.which("mpirun")
    if which_mpirun:
        return Path(which_mpirun)
    return software_dir / "openmpi" / "bin" / "mpirun"


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
    openmpi_path = _find_openmpi_mpirun()
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
        if installer_result and installer_result.returncode != 0:
            logger.warning("QM tools installer exited with code %d", installer_result.returncode)
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
    tools: list[str] | None = None,
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

    cmd = ["bash", str(installer)]
    if tools:
        cmd.extend(tools)

    result = subprocess.run(
        cmd,
        cwd=str(target),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=600,
    )
    return target, result


def stage_packaged_csp_tools(target_dir: str | Path | None = None) -> Path:
    source = get_packaged_csp_tools_dir()
    if not source.is_dir():
        raise FileNotFoundError(f"Packaged csp_tools directory not found: {source}")

    target = get_user_csp_tools_dir(target_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        source,
        target,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo", ".build"),
    )
    return target


def run_csp_tools_installer(
    target_dir: str | Path | None = None,
    *,
    extra_env: dict[str, str] | None = None,
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    target = stage_packaged_csp_tools(target_dir)
    installer = target / "install_csp_tools.sh"
    if not installer.is_file():
        raise FileNotFoundError(f"csp_tools installer not found: {installer}")

    env = os.environ.copy()
    env["DELFIN_CSP_TOOLS_ROOT"] = str(target)
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
        timeout=1200,
    )
    return target, result


def stage_packaged_mlp_tools(target_dir: str | Path | None = None) -> Path:
    source = get_packaged_mlp_tools_dir()
    if not source.is_dir():
        raise FileNotFoundError(f"Packaged mlp_tools directory not found: {source}")

    target = get_user_mlp_tools_dir(target_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        source,
        target,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
    )
    return target


def run_mlp_tools_installer(
    target_dir: str | Path | None = None,
    *,
    extra_env: dict[str, str] | None = None,
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    target = stage_packaged_mlp_tools(target_dir)
    installer = target / "install_mlp_tools.sh"
    if not installer.is_file():
        raise FileNotFoundError(f"mlp_tools installer not found: {installer}")

    env = os.environ.copy()
    env["DELFIN_MLP_TOOLS_ROOT"] = str(target)
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
        timeout=600,
    )
    return target, result


def stage_packaged_analysis_tools(target_dir: str | Path | None = None) -> Path:
    source = get_packaged_analysis_tools_dir()
    if not source.is_dir():
        raise FileNotFoundError(f"Packaged analysis_tools directory not found: {source}")

    target = get_user_analysis_tools_dir(target_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        source,
        target,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
    )
    return target


def run_analysis_tools_installer(
    target_dir: str | Path | None = None,
    *,
    extra_env: dict[str, str] | None = None,
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    target = stage_packaged_analysis_tools(target_dir)
    installer = target / "install_analysis_tools.sh"
    if not installer.is_file():
        raise FileNotFoundError(f"analysis_tools installer not found: {installer}")

    env = os.environ.copy()
    env["DELFIN_ANALYSIS_TOOLS_ROOT"] = str(target)
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
        timeout=600,
    )
    return target, result


def stage_packaged_ai_tools(target_dir: str | Path | None = None) -> Path:
    source = get_packaged_ai_tools_dir()
    if not source.is_dir():
        raise FileNotFoundError(f"Packaged ai_tools directory not found: {source}")

    target = get_user_ai_tools_dir(target_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        source,
        target,
        dirs_exist_ok=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
    )
    return target


def run_ai_tools_installer(
    target_dir: str | Path | None = None,
    *,
    extra_env: dict[str, str] | None = None,
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    target = stage_packaged_ai_tools(target_dir)
    installer = target / "install_ai_tools.sh"
    if not installer.is_file():
        raise FileNotFoundError(f"ai_tools installer not found: {installer}")

    env = os.environ.copy()
    env["DELFIN_AI_TOOLS_ROOT"] = str(target)
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
        timeout=600,
    )
    return target, result


def run_pip_install_editable() -> subprocess.CompletedProcess[str]:
    """Run ``pip install -e .`` in the DELFIN repo root using the current Python.

    Detects the Python interpreter that is running DELFIN and the repository
    root (the directory containing ``pyproject.toml``).  This ensures that
    editable-mode install happens in the correct environment (venv / conda /
    system) so that new CLI entry-points become available immediately.
    """
    import sys

    python_bin = sys.executable
    # Walk upward from this file to find the repo root (contains pyproject.toml)
    repo_root = Path(__file__).resolve().parent.parent
    pyproject = repo_root / "pyproject.toml"
    if not pyproject.is_file():
        raise FileNotFoundError(
            f"Cannot locate pyproject.toml — expected at {pyproject}. "
            "Make sure DELFIN is installed from a source checkout."
        )

    return subprocess.run(
        [python_bin, "-m", "pip", "install", "-e", "."],
        cwd=str(repo_root),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=300,
    )


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


def apply_runtime_environment(*, qm_tools_root: str = "", orca_base: str = "", csp_tools_root: str = "") -> None:
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

    csp_root = normalize_runtime_path(csp_tools_root)
    if csp_root:
        os.environ["DELFIN_CSP_TOOLS_ROOT"] = csp_root
        csp_bin = Path(csp_root) / "bin"
        if csp_bin.exists():
            _prepend_path_env_once(str(csp_bin))

    orca_root = normalize_runtime_path(orca_base)
    if orca_root:
        os.environ["DELFIN_ORCA_BASE"] = orca_root
        os.environ["ORCA_PATH"] = orca_root
        _prepend_path_env_once(orca_root)
        orca_plot = Path(orca_root) / "orca_plot"
        if orca_plot.exists():
            os.environ["ORCA_PLOT"] = str(orca_plot)

    # Suppress noisy but harmless third-party warnings in the dashboard
    os.environ.setdefault("TORCHANI_NO_WARN_EXTENSIONS", "1")


@contextmanager
def temporary_environment(updates: dict[str, str | None]):
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

    # -- CSP tools (Genarris) diagnostics ---------------------------------
    csp_tools_root = normalize_runtime_path(runtime_settings.get("csp_tools_root", ""))
    try:
        from delfin.csp_tools import genarris_available, get_genarris_version

        genarris_ok = genarris_available()
        genarris_ver = get_genarris_version() or ""
        diagnostics.append(
            {
                "name": "genarris",
                "status": "ok" if genarris_ok else "missing",
                "detail": f"gnrs v{genarris_ver}" if genarris_ok else "not installed (pip install delfin-complat[csp])",
            }
        )
        if genarris_ok:
            try:
                # Verify C extension loads
                import gnrs.cgenarris  # noqa: F401

                diagnostics.append(
                    {"name": "cgenarris", "status": "ok", "detail": "C extension OK"}
                )
            except ImportError as exc:
                diagnostics.append(
                    {"name": "cgenarris", "status": "missing", "detail": f"C extension failed: {exc}"}
                )
    except ImportError:
        pass

    csp_env_updates: dict[str, str] = {}
    if csp_tools_root:
        csp_env_updates["DELFIN_CSP_TOOLS_ROOT"] = csp_tools_root
    with temporary_environment(csp_env_updates):
        gnrs_results = check_tools(["gnrs"])
    for name, resolved in gnrs_results:
        diagnostics.append(
            {
                "name": "gnrs-cli",
                "status": "ok" if resolved else "missing",
                "detail": resolved.path if resolved else "not found",
            }
        )

    if csp_tools_root:
        diagnostics.append(
            {
                "name": "csp_tools_root",
                "status": "ok" if Path(csp_tools_root).exists() else "missing",
                "detail": csp_tools_root,
            }
        )

    # -- External QM programs (auto-detect) --------------------------------
    import shutil as _shutil_diag

    _EXT_QM_PROGRAMS = [
        # ── QM engines (ab initio / DFT / semiempirical) ──
        (["turbomole", "ridft", "dscf", "define"], "turbomole"),
        (["g16", "g09", "gaussian"], "gaussian"),
        (["nwchem"], "nwchem"),
        (["qchem"], "qchem"),
        (["gamess", "rungms"], "gamess"),
        (["molpro", "molpro.exe"], "molpro"),
        (["dalton", "dalton.x"], "dalton"),
        (["psi4"], "psi4"),
        (["cfour", "xcfour"], "cfour"),
        (["mrcc", "dmrcc"], "mrcc"),
        # ── Periodic DFT / solid state ──
        (["vasp", "vasp_std", "vasp_gam", "vasp_ncl"], "vasp"),
        (["pw.x", "pw", "quantum-espresso"], "quantum-espresso"),
        (["cp2k", "cp2k.popt", "cp2k.psmp", "cp2k.sopt"], "cp2k"),
        (["aims.x", "FHIaims", "aims"], "fhi-aims"),
        (["crystal", "crystal17", "crystal23", "Pcrystal"], "crystal"),
        (["siesta"], "siesta"),
        (["gpaw", "gpaw-python"], "gpaw"),
        (["fleur", "fleur_MPI"], "fleur"),
        (["wien2k", "run_lapw", "runsp_lapw"], "wien2k"),
        (["elk"], "elk"),
        (["abinit"], "abinit"),
        # ── Multireference / correlated ──
        (["pymolcas", "molcas", "OpenMolcas"], "openmolcas"),
        (["bagel", "BAGEL"], "bagel"),
        (["columbus", "COLUMBUS"], "columbus"),
        # ── Semiempirical ──
        (["mopac", "MOPAC", "mopac2016"], "mopac"),
        (["sparrow"], "sparrow"),
        # ── MD engines ──
        (["gmx", "gmx_mpi", "gromacs", "mdrun"], "gromacs"),
        (["lmp", "lmp_mpi", "lammps", "lmp_serial"], "lammps"),
        (["sander", "pmemd", "pmemd.cuda"], "amber"),
        (["namd2", "namd3", "namd"], "namd"),
        (["openmmrun"], "openmm"),
        # ── Visualization / editors ──
        (["vmd"], "vmd"),
        (["avogadro", "avogadro2"], "avogadro"),
        (["jmol"], "jmol"),
        (["chimera", "chimerax", "ChimeraX"], "chimerax"),
        (["iqmol", "IQmol"], "iqmol"),
        # ── Analysis / post-processing ──
        (["orca_2mkl"], "orca_2mkl"),
        (["nbo7", "nbo6", "gennbo"], "nbo"),
        (["aimall", "aimqb"], "aimall"),
        (["critic2"], "critic2"),
        (["chargemol"], "chargemol"),
        (["lobster", "LOBSTER"], "lobster"),
        (["phonopy"], "phonopy"),
        (["ase"], "ase-cli"),
    ]
    # Also check Python-only tools that may not have CLI binaries
    _EXT_PYTHON_MODULES = [
        ("openmm", "openmm"),
        ("pyscf", "pyscf"),
        ("pymatgen", "pymatgen"),
        ("qcengine", "qcengine"),
        ("MDAnalysis", "mdanalysis"),
        ("plams", "plams"),
    ]
    import importlib.util as _ilu_diag
    for mod_name, label in _EXT_PYTHON_MODULES:
        try:
            found = _ilu_diag.find_spec(mod_name) is not None
        except (ModuleNotFoundError, ValueError):
            found = False
        # Don't duplicate if already found via binary check
        if any(d["name"] == label for d in diagnostics):
            continue
        diagnostics.append(
            {
                "name": label,
                "status": "ok" if found else "missing",
                "detail": f"Python module '{mod_name}'" if found else "not installed",
            }
        )
    # -- bwHPC / HPC cluster fallback paths ----------------------------------
    # Map diagnostic label -> (env_var, base_dirs, binary_name_or_glob_pattern)
    # Used when shutil.which() fails (programs behind 'module load').
    _BWHPC_ROOTS = ["/opt/bwhpc/common/chem", "/opt/bwhpc/common/phys",
                    "/opt/bwhpc/common/bio"]
    _CLUSTER_FALLBACKS: dict[str, tuple[str, list[str], str]] = {
        "turbomole": ("TURBODIR", ["/opt/bwhpc/common/chem/turbomole"], "ridft"),
        "gaussian":  ("g16root",  ["/opt/bwhpc/common/chem/gaussian"], "g16"),
        "vasp":      ("VASP_DIR", ["/opt/bwhpc/common/chem/vasp"], "vasp_std"),
        "nwchem":    ("NWCHEM_TOP", ["/opt/bwhpc/common/chem/nwchem"], "nwchem"),
        "cp2k":      ("CP2K_DIR", ["/opt/bwhpc/common/chem/cp2k"], "cp2k.psmp"),
        "gromacs":   ("GMXBIN",   ["/opt/bwhpc/common/chem/gromacs",
                                    "/opt/bwhpc/common/bio/gromacs"], "gmx"),
        "lammps":    ("LAMMPS_DIR", ["/opt/bwhpc/common/chem/lammps"], "lmp"),
        "amber":     ("AMBERHOME", ["/opt/bwhpc/common/chem/amber",
                                     "/opt/bwhpc/common/bio/amber"], "sander"),
        "namd":      ("NAMDDIR",  ["/opt/bwhpc/common/chem/namd",
                                    "/opt/bwhpc/common/bio/namd"], "namd2"),
        "quantum-espresso": ("ESPRESSO_ROOT", ["/opt/bwhpc/common/chem/quantum-espresso",
                                               "/opt/bwhpc/common/phys/quantum-espresso"], "pw.x"),
        "abinit":    ("ABINIT_DIR", ["/opt/bwhpc/common/phys/abinit",
                                      "/opt/bwhpc/common/chem/abinit"], "abinit"),
        "molpro":    ("MOLPRO_DIR", ["/opt/bwhpc/common/chem/molpro"], "molpro"),
        "openmolcas": ("MOLCAS",   ["/opt/bwhpc/common/chem/openmolcas",
                                     "/opt/bwhpc/common/chem/molcas"], "pymolcas"),
    }

    def _find_binary_recursive(base_dir: Path, binary_name: str, max_depth: int = 5) -> str | None:
        """Walk up to *max_depth* levels under *base_dir* looking for *binary_name*."""
        if not base_dir.is_dir():
            return None
        try:
            for item in sorted(base_dir.iterdir(), reverse=True):
                if item.is_file() and item.name == binary_name and os.access(str(item), os.X_OK):
                    return str(item)
                if item.is_dir() and max_depth > 0:
                    # Prioritise bin/ directories
                    if item.name == "bin":
                        result = _find_binary_recursive(item, binary_name, max_depth - 1)
                        if result:
                            return result
            # Second pass: recurse into non-bin dirs
            for item in sorted(base_dir.iterdir(), reverse=True):
                if item.is_dir() and item.name != "bin" and max_depth > 0:
                    result = _find_binary_recursive(item, binary_name, max_depth - 1)
                    if result:
                        return result
        except PermissionError:
            pass
        return None

    for binaries, label in _EXT_QM_PROGRAMS:
        found_path = None
        for name in binaries:
            path = _shutil_diag.which(name)
            if path:
                found_path = path
                break
        # Fallback: check environment variables and known HPC cluster paths
        if not found_path and label in _CLUSTER_FALLBACKS:
            env_var, base_dirs, target_bin = _CLUSTER_FALLBACKS[label]
            # Check environment variable first
            env_val = os.environ.get(env_var, "")
            if env_val:
                result = _find_binary_recursive(Path(env_val), target_bin)
                if result:
                    found_path = result
            # Then check known cluster paths
            if not found_path:
                for bdir in base_dirs:
                    result = _find_binary_recursive(Path(bdir), target_bin)
                    if result:
                        found_path = result
                        break
        diagnostics.append(
            {
                "name": label,
                "status": "ok" if found_path else "missing",
                "detail": found_path or "not found",
            }
        )

    # -- MLP tools diagnostics --------------------------------------------
    try:
        from delfin.mlp_tools import (
            torchani_available,
            aimnet2_available,
            mace_available,
            chgnet_available,
            matgl_available,
            schnetpack_available,
            nequip_available,
            alignn_available,
            get_torchani_version,
            get_aimnet2_version,
            get_mace_version,
            get_chgnet_version,
            get_matgl_version,
            get_schnetpack_version,
            get_nequip_version,
            get_alignn_version,
        )

        for check_fn, ver_fn, label in [
            (torchani_available, get_torchani_version, "ani2x"),
            (aimnet2_available, get_aimnet2_version, "aimnet2"),
            (mace_available, get_mace_version, "mace-off"),
            (chgnet_available, get_chgnet_version, "chgnet"),
            (matgl_available, get_matgl_version, "m3gnet"),
            (schnetpack_available, get_schnetpack_version, "schnetpack"),
            (nequip_available, get_nequip_version, "nequip"),
            (alignn_available, get_alignn_version, "alignn"),
        ]:
            ok = check_fn()
            ver = ver_fn() or ""
            diagnostics.append(
                {
                    "name": label,
                    "status": "ok" if ok else "missing",
                    "detail": f"v{ver}" if ok else "not installed",
                }
            )
    except ImportError:
        pass

    # -- Analysis tools diagnostics ----------------------------------------
    try:
        from delfin.analysis_tools import (
            multiwfn_available,
            censo_available,
            morfeus_available,
            cclib_available,
            nglview_available,
            packmol_available,
            get_multiwfn_version,
            get_censo_version,
            get_morfeus_version,
            get_cclib_version,
            get_nglview_version,
            get_packmol_version,
        )

        for check_fn, ver_fn, label in [
            (multiwfn_available, get_multiwfn_version, "multiwfn"),
            (censo_available, get_censo_version, "censo"),
            (morfeus_available, get_morfeus_version, "morfeus"),
            (cclib_available, get_cclib_version, "cclib"),
            (nglview_available, get_nglview_version, "nglview"),
            (packmol_available, get_packmol_version, "packmol"),
        ]:
            ok = check_fn()
            ver = ver_fn() or ""
            diagnostics.append(
                {
                    "name": label,
                    "status": "ok" if ok else "missing",
                    "detail": f"v{ver}" if ok else "not installed",
                }
            )
    except ImportError:
        pass

    # -- AI tools diagnostics -----------------------------------------------
    try:
        from delfin.ai_tools import _TOOL_REGISTRY

        for label, category, avail_fn, ver_fn, description, install_hint in _TOOL_REGISTRY:
            ok = avail_fn()
            ver = ver_fn() or ""
            diagnostics.append(
                {
                    "name": label.lower().replace(" ", "-"),
                    "status": "ok" if ok else "missing",
                    "detail": f"v{ver}" if ok else "not installed",
                }
            )
    except ImportError:
        pass

    return diagnostics
