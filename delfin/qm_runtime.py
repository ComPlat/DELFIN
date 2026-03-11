from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from shutil import which
from typing import Iterable, Mapping, Optional, Sequence

from delfin.common.logging import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class ToolSpec:
    name: str
    env_vars: tuple[str, ...] = ()
    which_targets: tuple[str, ...] = ()
    aliases: tuple[str, ...] = ()
    prefer_qm_tools: bool = True
    locator_command: Optional[str] = None


@dataclass(frozen=True)
class ResolvedTool:
    requested_name: str
    canonical_name: str
    path: str
    source: str


_TOOL_SPECS: tuple[ToolSpec, ...] = (
    ToolSpec(
        name="orca",
        env_vars=("ORCA_BINARY", "ORCA_PATH"),
        which_targets=("orca", "orca.exe"),
        aliases=("orca.exe",),
        prefer_qm_tools=False,
        locator_command="orca_locate",
    ),
    ToolSpec(name="xtb", which_targets=("xtb",), aliases=("gfnxtb",)),
    ToolSpec(name="crest", which_targets=("crest",)),
    ToolSpec(name="std2", which_targets=("std2",)),
    ToolSpec(name="stda", which_targets=("stda",)),
    ToolSpec(name="xtb4stda", which_targets=("xtb4stda",)),
    ToolSpec(name="dftb+", which_targets=("dftb+",), aliases=("dftbplus",)),
)

_SPEC_BY_CANONICAL = {spec.name: spec for spec in _TOOL_SPECS}
_ALIAS_TO_CANONICAL = {
    alias: spec.name for spec in _TOOL_SPECS for alias in (spec.name,) + spec.aliases
}


def get_qm_tools_root() -> Path:
    env_root = os.environ.get("DELFIN_QM_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()
    return (Path(__file__).resolve().parent / "qm_tools").resolve()


def get_qm_tools_bin_dir() -> Path:
    return get_qm_tools_root() / "bin"


def canonical_tool_name(name: str) -> str:
    normalized = str(name or "").strip()
    if not normalized:
        raise ValueError("Tool name must not be empty")
    return _ALIAS_TO_CANONICAL.get(normalized, normalized)


def supported_qm_tools() -> list[str]:
    names = [spec.name for spec in _TOOL_SPECS]
    bin_dir = get_qm_tools_bin_dir()
    if bin_dir.is_dir():
        for entry in sorted(bin_dir.iterdir(), key=lambda item: item.name):
            if not entry.is_file() or not os.access(entry, os.X_OK):
                continue
            if entry.name not in names:
                names.append(entry.name)
    return names


def _validate_candidate(candidate: str) -> Optional[str]:
    if not candidate:
        return None
    expanded = Path(candidate.strip()).expanduser()
    if not expanded.is_file():
        return None
    if not os.access(expanded, os.X_OK):
        return None
    return str(expanded.resolve())


def _iter_locator_candidates(locator: str) -> Iterable[str]:
    locator_path = which(locator)
    if not locator_path:
        return
    try:
        result = subprocess.run([locator_path], check=False, capture_output=True, text=True)
    except Exception as exc:  # noqa: BLE001
        logger.debug("Failed to query %s: %s", locator, exc)
        return
    if result.returncode != 0:
        logger.debug("%s returned non-zero exit status %s", locator, result.returncode)
        return
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if stripped:
            yield stripped


def _iter_tool_candidates(spec: ToolSpec) -> Iterable[tuple[str, str]]:
    bin_dir = get_qm_tools_bin_dir()
    if spec.prefer_qm_tools:
        local_candidate = bin_dir / spec.name
        if local_candidate.exists():
            yield str(local_candidate), "qm_tools"

    generic_env_key = f"DELFIN_{spec.name.upper().replace('+', 'PLUS').replace('-', '_')}_BINARY"
    for key in (generic_env_key,) + spec.env_vars:
        value = os.environ.get(key)
        if value:
            yield value, f"env:{key}"

    for target in spec.which_targets:
        located = which(target)
        if located:
            yield located, "PATH"

    if spec.locator_command:
        for candidate in _iter_locator_candidates(spec.locator_command):
            yield candidate, spec.locator_command

    if not spec.prefer_qm_tools:
        local_candidate = bin_dir / spec.name
        if local_candidate.exists():
            yield str(local_candidate), "qm_tools"


def _generic_tool_spec(name: str) -> ToolSpec:
    return ToolSpec(name=name, which_targets=(name,))


def resolve_tool(name: str) -> Optional[ResolvedTool]:
    direct_candidate = _validate_candidate(str(name))
    if direct_candidate:
        return ResolvedTool(
            requested_name=name,
            canonical_name=Path(direct_candidate).name,
            path=direct_candidate,
            source="explicit",
        )

    canonical = canonical_tool_name(name)
    spec = _SPEC_BY_CANONICAL.get(canonical, _generic_tool_spec(canonical))
    for candidate, source in _iter_tool_candidates(spec):
        valid = _validate_candidate(candidate)
        if valid:
            return ResolvedTool(
                requested_name=name,
                canonical_name=canonical,
                path=valid,
                source=source,
            )
        logger.debug("Discarding invalid %s candidate path: %r", canonical, candidate)
    return None


def find_tool_executable(name: str) -> Optional[str]:
    resolved = resolve_tool(name)
    return resolved.path if resolved else None


def check_tools(names: Optional[Sequence[str]] = None) -> list[tuple[str, Optional[ResolvedTool]]]:
    selected = list(names) if names else supported_qm_tools()
    results: list[tuple[str, Optional[ResolvedTool]]] = []
    for name in selected:
        results.append((name, resolve_tool(name)))
    return results


def _resolve_xtb4stda_home(root: Path) -> str:
    existing = os.environ.get("XTB4STDAHOME")
    if existing:
        return existing

    default_home = root / "share" / "xtb4stda"
    if not default_home.exists():
        return str(default_home)

    default_str = str(default_home)
    if len(default_str) <= 60:
        return default_str

    short_home = Path.home() / ".delfin_xtb4stda"
    try:
        if short_home.is_dir() and (short_home / ".param_stda1.xtb").exists():
            return str(short_home)
        if not short_home.exists() or short_home.is_symlink():
            if short_home.is_symlink() or short_home.exists():
                short_home.unlink()
            short_home.symlink_to(default_home)
            return str(short_home)
    except OSError as exc:
        logger.debug("Failed to prepare short xtb4stda runtime path: %s", exc)

    return default_str


def _prepare_tool_environment(extra_env: Optional[Mapping[str, str]] = None) -> dict[str, str]:
    env = os.environ.copy()
    root = get_qm_tools_root()
    bin_dir = str(root / "bin")
    current_path = env.get("PATH", "")
    env["PATH"] = bin_dir if not current_path else f"{bin_dir}{os.pathsep}{current_path}"
    share_xtb4stda = root / "share" / "xtb4stda"
    if share_xtb4stda.exists():
        env.setdefault("XTB4STDAHOME", _resolve_xtb4stda_home(root))
    env.setdefault("STD2HOME", str(root))
    env.setdefault("OMP_STACKSIZE", "4G")
    if extra_env:
        env.update({str(k): str(v) for k, v in extra_env.items()})
    return env


def run_tool(
    name: str,
    args: Sequence[str],
    *,
    cwd: Optional[str | Path] = None,
    env: Optional[Mapping[str, str]] = None,
    check: bool = False,
    capture_output: bool = False,
    text: bool = True,
    stdout=None,
    stderr=None,
    track_process: bool = False,
) -> subprocess.CompletedProcess:
    resolved = resolve_tool(name)
    if resolved is None:
        raise FileNotFoundError(f"QM tool not found: {name}")

    cmd = [resolved.path, *[str(arg) for arg in args]]
    env_map = _prepare_tool_environment(env)
    run_kwargs = {
        "cwd": str(cwd) if cwd is not None else None,
        "env": env_map,
        "check": check,
        "text": text,
    }
    if capture_output:
        run_kwargs["capture_output"] = True
    else:
        run_kwargs["stdout"] = stdout
        run_kwargs["stderr"] = stderr

    logger.info("Running QM tool %s from %s: %s", resolved.canonical_name, resolved.source, cmd)
    if not track_process:
        return subprocess.run(cmd, **run_kwargs)

    popen_kwargs = {
        "cwd": run_kwargs["cwd"],
        "env": env_map,
        "text": text,
        "start_new_session": True,
    }
    if capture_output:
        popen_kwargs["stdout"] = subprocess.PIPE
        popen_kwargs["stderr"] = subprocess.PIPE
    else:
        popen_kwargs["stdout"] = stdout
        popen_kwargs["stderr"] = stderr

    manager = None
    registration_token = None
    try:
        from delfin.global_manager import get_global_manager

        manager = get_global_manager()
    except Exception as exc:  # noqa: BLE001
        logger.debug("Failed to acquire global manager for %s: %s", resolved.canonical_name, exc)

    process = subprocess.Popen(cmd, **popen_kwargs)
    if manager is not None:
        try:
            registration_token = manager.register_subprocess(
                process,
                label=" ".join(cmd[:2]),
                cwd=run_kwargs["cwd"],
            )
        except Exception as exc:  # noqa: BLE001
            logger.debug("Failed to register QM tool subprocess %s: %s", resolved.canonical_name, exc)

    try:
        stdout_data, stderr_data = process.communicate()
    finally:
        if manager is not None:
            try:
                manager.unregister_subprocess(registration_token)
            except Exception as exc:  # noqa: BLE001
                logger.debug("Failed to unregister QM tool subprocess %s: %s", resolved.canonical_name, exc)

    result = subprocess.CompletedProcess(
        cmd,
        process.returncode,
        stdout_data if capture_output else None,
        stderr_data if capture_output else None,
    )
    if check and process.returncode != 0:
        raise subprocess.CalledProcessError(
            process.returncode,
            cmd,
            output=result.stdout,
            stderr=result.stderr,
        )
    return result
