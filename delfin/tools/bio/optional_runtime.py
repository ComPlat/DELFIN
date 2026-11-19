"""Run a fixed installed console entry point in a chosen Python environment."""
from __future__ import annotations

from pathlib import Path
import subprocess

from delfin.tools.bio.artifacts import file_hash


_VERSION = "import importlib.metadata as m,sys; print(m.version(sys.argv[1]))"
_ENTRY = ("import importlib.metadata as m,sys; "
          "d=m.distribution(sys.argv.pop(1)); n=sys.argv.pop(1); "
          "e=next(e for e in d.entry_points if e.group=='console_scripts' and e.name==n); "
          "sys.argv[0]=n; sys.exit(e.load()())")


def probe_distribution(distribution: str, *, python: str, expected_version: str, work_dir: Path) -> dict:
    from delfin import qm_runtime
    resolved = qm_runtime.resolve_tool(python)
    if resolved is None:
        raise FileNotFoundError(f"Python environment unavailable: {python}")
    version = qm_runtime.run_tool(resolved.path, ["-I", "-c", _VERSION, distribution], cwd=work_dir,
                                  capture_output=True, timeout=30, track_process=True, allow_install=False)
    if version.returncode:
        raise ImportError(f"{distribution} unavailable in {resolved.path}: {version.stderr}")
    actual = (version.stdout or "").strip()
    if actual != expected_version:
        raise ValueError(f"{distribution} version mismatch: expected {expected_version!r}, found {actual!r}")
    return {"distribution": distribution, "version": actual, "python": resolved.path,
            "python_sha256": file_hash(Path(resolved.path))}


def run_console(distribution: str, entry: str, args: list[str], *, python: str,
                expected_version: str, work_dir: Path, timeout: float) -> dict:
    from delfin import qm_runtime
    runtime = probe_distribution(distribution, python=python, expected_version=expected_version, work_dir=work_dir)
    argv = ["-I", "-c", _ENTRY, distribution, entry, *map(str, args)]
    with (work_dir / "engine.log").open("w") as log:
        result = qm_runtime.run_tool(runtime["python"], argv, cwd=work_dir, stdout=log, stderr=subprocess.STDOUT,
                                     timeout=timeout, track_process=True, allow_install=False)
    if result.returncode:
        raise RuntimeError(f"{entry} exited {result.returncode}; see engine.log")
    return {**runtime, "entry_point": entry, "arguments": list(map(str, args))}
