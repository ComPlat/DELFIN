import os
import shutil
import signal
import subprocess
import time
from pathlib import Path
from typing import Iterable, NamedTuple, Optional, Sequence, Set

from delfin.common.logging import get_logger

DEFAULT_PATTERNS: Sequence[str] = [
    "*.cpcm",
    "*.cpcm_corr",
    "*.tmp",
    "*.tmp*",
    "*.bas*",
    "*_D0*",
    "XTB_GOAT.goat.*.out",
    # ExtOpt exchange files of an MLIP GOAT (global_opt_method=mlip:*),
    # written by ORCA 6.1.1 as <base>_EXT.extinp.tmp / <base>_EXT.xyz.
    "XTB_GOAT_EXT.*",
]

logger = get_logger(__name__)


class ProcessInfo(NamedTuple):
    pid: int
    pgid: int
    cmdline: str
    cwd: Path


_ORCA_BIN_NAMES: set[str] = {
    "orca",
    "orca.exe",
    "orca_main",
    "orca_scfgrad",
    "orca_scfgrad_mpi",
    "orca_leanscf",
    "orca_leanscf_mpi",
    "orca_plot",
    "mpirun",
    "orterun",
}

OCCUPIER_DIR_PATTERNS: Sequence[str] = (
    "initial_OCCUPIER",
    "ox_step_*_OCCUPIER",
    "red_step_*_OCCUPIER",
)

def cleanup(folder: str = ".",
            recursive: bool = False,
            dry_run: bool = False,
            patterns: Optional[Sequence[str]] = None) -> int:
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"{root} ist kein Ordner.")

    pats = list(patterns or DEFAULT_PATTERNS)

    def iter_files():
        for pat in pats:
            yield from (root.rglob(pat) if recursive else root.glob(pat))

    to_delete = sorted({p for p in iter_files() if p.is_file()})
    count = 0
    for f in to_delete:
        if dry_run:
            # don't spend anything
            pass
        else:
            try:
                f.unlink()
                count += 1
            except Exception as exc:  # noqa: BLE001
                logger.debug("Failed to delete %s during cleanup: %s", f, exc)
    return count


def cleanup_all(folder: str = ".",
                dry_run: bool = False,
                patterns: Optional[Sequence[str]] = None) -> int:
    return cleanup(folder=folder, recursive=True, dry_run=dry_run, patterns=patterns)


def _read_process_cwd(pid: int) -> Optional[Path]:
    """Return the working directory for a process if accessible."""
    try:
        link = os.readlink(f"/proc/{pid}/cwd")
    except (FileNotFoundError, PermissionError, OSError):
        return None
    try:
        return Path(link).resolve()
    except OSError:
        return None


def _collect_orca_processes(target_roots: Iterable[Path]) -> list[ProcessInfo]:
    """Gather ORCA-related processes whose cwd lives under target_roots."""
    roots: Set[Path] = {root.resolve() for root in target_roots}
    try:
        result = subprocess.run(
            ["ps", "-eo", "pid=,pgid=,comm=,args="],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to enumerate processes: %s", exc)
        return []

    if result.returncode != 0 or not result.stdout:
        return []

    processes: list[ProcessInfo] = []
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split(None, 3)
        if len(parts) < 3:
            continue
        try:
            pid = int(parts[0])
            pgid = int(parts[1])
        except ValueError:
            continue

        cmd_fragment = parts[2]
        args_fragment = parts[3] if len(parts) > 3 else ""
        cmdline = f"{cmd_fragment} {args_fragment}".strip()

        # Avoid self-termination and avoid matching flag names like "--orca"
        if not _looks_like_orca_process(cmd_fragment, args_fragment):
            continue

        cwd = _read_process_cwd(pid)
        if cwd is None:
            continue

        if not any(root == cwd or root in cwd.parents for root in roots):
            continue

        processes.append(ProcessInfo(pid=pid, pgid=pgid, cmdline=cmdline, cwd=cwd))

    return processes


def _terminate_process_groups(processes: Sequence[ProcessInfo], *, dry_run: bool = False) -> list[dict]:
    """Terminate ORCA process groups and return termination summaries."""
    terminated: list[dict] = []
    groups: dict[int, list[ProcessInfo]] = {}
    for proc in processes:
        groups.setdefault(proc.pgid, []).append(proc)

    for pgid, members in groups.items():
        summary = {
            "pgid": pgid,
            "pids": [m.pid for m in members],
            "cwd": str(members[0].cwd),
            "command": members[0].cmdline,
            "status": "pending",
        }
        if dry_run:
            summary["status"] = "dry-run"
            terminated.append(summary)
            continue

        logger.info(
            "Terminating ORCA process group %s (pids=%s, cmd=%s)",
            pgid,
            summary["pids"],
            summary["command"],
        )

        try:
            os.killpg(pgid, signal.SIGTERM)
            summary["status"] = "sigterm-sent"
        except ProcessLookupError:
            summary["status"] = "missing"
            terminated.append(summary)
            continue
        except PermissionError as exc:
            summary["status"] = f"permission-error: {exc}"
            terminated.append(summary)
            continue
        except Exception as exc:  # noqa: BLE001
            summary["status"] = f"error: {exc}"
            terminated.append(summary)
            continue

        deadline = time.time() + 5.0
        while time.time() < deadline:
            if all(not Path(f"/proc/{pid}").exists() for pid in summary["pids"]):
                summary["status"] = "terminated"
                break
            time.sleep(0.2)

        if summary["status"] != "terminated":
            try:
                os.killpg(pgid, signal.SIGKILL)
                summary["status"] = "sigkilled"
            except ProcessLookupError:
                summary["status"] = "terminated"
            except Exception as exc:  # noqa: BLE001
                summary["status"] = f"kill-error: {exc}"

        terminated.append(summary)

    return terminated


def _is_kept_stage_dir(path: Path) -> bool:
    """True for an ``*_OCCUPIER`` stage directory, which cleanup must keep.

    Those folders are resume state, not scratch. A resumed OCCUPIER stage sets
    itself up from them: ``_occupier_results_fresh``
    (delfin/workflows/pipeline.py:139) reuses a stage only while its
    ``OCCUPIER.txt`` and ``.control_fingerprint`` are there, and fails closed
    to a full rerun otherwise; ``read_occ``
    (delfin/workflows/engine/occupier.py:646) raises without ``OCCUPIER.txt``
    and its caller at :738 silently falls back to ``read_occ_from_control``, a
    different electron configuration than the stage chose; the per-stage
    ``CONTROL.txt`` and the ``.gbw``/``.xyz`` a stage restarts from live there
    too. Only the temporary files inside them are scratch, and ``cleanup``
    deletes files only, so it may still remove those.
    """
    return path.is_dir() and path.name.endswith("_OCCUPIER")


def _remove_occuper_directories(workspace: Path, *, dry_run: bool = False) -> list[str]:
    """Return the stage directories cleanup refuses to delete; never removes them.

    Kept as a named no-op so the report shape and its two printers
    (delfin/cli.py ``cleanup --orca`` and ``stop --cleanup``) keep their
    meaning. The deletion this function used to perform is the data loss the
    stage directories' docstrings describe.
    """
    kept = [
        str(entry)
        for pattern in OCCUPIER_DIR_PATTERNS
        for entry in sorted(workspace.glob(pattern))
        if _is_kept_stage_dir(entry)
    ]
    if kept and dry_run:
        logger.info("Keeping %d OCCUPIER stage directory/ies (resume state)", len(kept))
    return kept


def _remove_orca_scratch_dirs(scratch_root: Path, *, dry_run: bool = False) -> list[str]:
    removed: list[str] = []

    orca_scratch = scratch_root / ".orca_scratch"
    if not orca_scratch.exists():
        return removed

    if dry_run:
        removed.append(str(orca_scratch))
        return removed

    try:
        shutil.rmtree(orca_scratch)
        removed.append(str(orca_scratch))
        logger.info("Removed ORCA scratch directory %s", orca_scratch)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to remove scratch directory %s: %s", orca_scratch, exc)

    return removed


def cleanup_orca(workspace: str | Path = ".",
                 *,
                 scratch_root: Optional[str | Path] = None,
                 dry_run: bool = False) -> dict:
    """Terminate ORCA jobs and purge ORCA scratch, keeping OCCUPIER resume state.

    The ``*_OCCUPIER`` stage directories are kept (see
    ``_remove_occuper_directories``); what is purged is ``.orca_scratch`` and
    the ``DEFAULT_PATTERNS`` temporary files, including the copies inside a
    kept stage directory, which ``cleanup`` unlinks file by file.
    """
    workspace_path = Path(workspace).expanduser().resolve()
    scratch_path = Path(scratch_root).expanduser().resolve() if scratch_root else workspace_path

    target_roots = {workspace_path}
    if scratch_path != workspace_path:
        target_roots.add(scratch_path)

    processes = _collect_orca_processes(target_roots)
    termination_report = _terminate_process_groups(processes, dry_run=dry_run) if processes else []

    # Resume state, not scratch: reported so a caller can say what was spared.
    kept_occuper = _remove_occuper_directories(workspace_path, dry_run=dry_run)
    removed_scratch = _remove_orca_scratch_dirs(scratch_path, dry_run=dry_run)

    cleaned_files = 0
    unique_cleanup_roots = {workspace_path, scratch_path}
    if not dry_run:
        for root in unique_cleanup_roots:
            try:
                cleaned_files += cleanup_all(str(root), dry_run=False)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Cleanup failed in %s: %s", root, exc)

    return {
        "processes_found": len(processes),
        "terminated_groups": termination_report,
        "occuper_dirs_kept": kept_occuper,
        # The name this key carried while the function deleted these folders.
        # It still reports the same paths so a caller that reads it keeps
        # working, but the paths in it were kept, never removed; new code
        # should read ``occuper_dirs_kept``.
        "occuper_dirs_removed": [],
        "scratch_dirs_removed": removed_scratch,
        "files_removed": cleaned_files,
        "dry_run": dry_run,
        "workspace": str(workspace_path),
        "scratch_root": str(scratch_path),
    }
def _looks_like_orca_process(cmd_fragment: str, args_fragment: str) -> bool:
    """Return True only for actual ORCA/MPI workers, not for CLI flags like '--orca'."""
    tokens = []
    if cmd_fragment:
        tokens.append(cmd_fragment.lower())
    if args_fragment:
        tokens.extend(part.lower() for part in args_fragment.split())

    for tok in tokens:
        # Skip pure option tokens (avoids matching '--orca')
        if tok.startswith("-"):
            continue

        base = tok.rsplit("/", 1)[-1]
        if base in _ORCA_BIN_NAMES:
            return True
    return False
