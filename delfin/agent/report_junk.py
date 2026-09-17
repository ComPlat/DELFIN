"""Report workspace junk that ``git status`` hides.

Read-only diagnostic: walks a workspace tree (default: the current working
directory) and reports leftover artifacts WITHOUT relying on git at all:

- large files above a size threshold (default 100 MiB)
- ``__pycache__`` directories and ``*.pyc`` files
- QM run leftovers (``*.out``, ``*.engrad``, ``*.grad``, ``*.hess``, ORCA/xtb scratch dirs such
  as ``orca.tmp*``, ``*.tmp``, ``tmp_*``, ``xtbrestart``, ``wbo``, ``charges``)
- ``*.lock`` leftovers
- total byte size per top-level folder

Nothing is ever deleted, modified, or sent anywhere; the module only reads
file names, sizes and directory structure. Files that look like secrets
(``.env*``, ``*.key``, ``*.pem``, ``credentials*``, ``.ssh``) are skipped
entirely -- not read, not even listed.

Run as ``python -m delfin.agent.report_junk [path]``.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

__all__ = ["collect", "format_text", "main"]

LARGE_FILE_THRESHOLD = 100 * 1024 * 1024  # 100 MiB

# Files/dirs that must never be touched or reported (secret deny-list).
SECRET_NAMES = {".ssh"}
SECRET_SUFFIXES = (".env", ".key", ".pem")
SECRET_PATTERNS = ("credentials", "credential")

# ORCA / xTB run leftovers by exact name.
QM_EXACT_NAMES = {
    "NORMAL_TERMINATION",
    "xtbrestart",
    "wbo",
    "charges",
    "gfnff_topo",
    "NOT_CONVERGED",
    "ORCA_ERROR",
}
# ORCA / xTB leftovers by file suffix.
QM_SUFFIXES = (
    ".out", ".engrad", ".grad", ".hess", ".incopt", ".pc",
    ".tmp", ".prop", ".cis", ".unparsed",  # ORCA scratch parts
)
# Scratch directories by prefix/name pattern.
QM_DIR_NAMES = {"orca.tmp", "xtb_tmp", "tmp_orca"}
QM_DIR_PREFIXES = ("tmp_",)


def _is_secret(path: Path) -> bool:
    name = path.name
    if name in SECRET_NAMES or any(p in name.lower() for p in SECRET_PATTERNS):
        return True
    if name.endswith(SECRET_SUFFIXES):
        return True
    return ".env" in name and name.endswith(".env")


def _is_qm_file(path: Path) -> bool:
    name = path.name
    if name in QM_EXACT_NAMES:
        return True
    return name.endswith(QM_SUFFIXES) and not name.endswith(".pyc")


def _is_qm_dir(path: Path) -> bool:
    name = path.name
    if name in QM_DIR_NAMES or name.startswith("orca.tmp"):
        return True
    return name.startswith(QM_DIR_PREFIXES) and name != "tmp"


def collect(root: str | os.PathLike | None = None) -> dict:
    """Walk ``root`` (default: cwd) and inventory workspace junk.

    Returns a dict with keys:

    - ``root``: absolute path that was scanned
    - ``large_files``: [{path, size}] above ``large_file_threshold``
    - ``pycache_dirs``: [str] __pycache__ directories found
    - ``pyc_files``: [str] *.pyc / *.pyo files found
    - ``qm_leftovers``: [str] ORCA/xtb leftover files and scratch dirs
    - ``lock_files``: [str] *.lock files
    - ``folder_sizes``: {top-level folder name: total bytes}
    - ``total_bytes``: grand total size of all regular files seen
    - ``n_files`` / ``n_dirs``: counters
    """
    root_path = Path(root).resolve() if root is not None else Path.cwd()
    data = {
        "root": str(root_path),
        "large_file_threshold": LARGE_FILE_THRESHOLD,
        "large_files": [],
        "pycache_dirs": [],
        "pyc_files": [],
        "qm_leftovers": [],
        "lock_files": [],
        "folder_sizes": {},
        "total_bytes": 0,
        "n_files": 0,
        "n_dirs": 0,
    }
    _walk(root_path, root_path, data)
    return data


def _walk(root: Path, current: Path, data: dict) -> None:
    try:
        entries = sorted(os.scandir(current), key=lambda e: e.name)
    except OSError:
        return  # unreadable dir: report nothing rather than crash
    for entry in entries:
        path = Path(entry.path)
        if _is_secret(path):
            continue
        if entry.is_symlink():
            continue  # never follow links: size would be double-counted
        if entry.is_dir(follow_symlinks=False):
            data["n_dirs"] += 1
            rel = str(path.relative_to(root))
            if path.name == "__pycache__":
                data["pycache_dirs"].append(rel)
                continue  # contents are bytecode, no need to descend
            if _is_qm_dir(path):
                data["qm_leftovers"].append(rel + "/")
            _walk(root, path, data)
        elif entry.is_file(follow_symlinks=False):
            data["n_files"] += 1
            try:
                size = entry.stat().st_size
            except OSError:
                continue
            rel = str(path.relative_to(root))
            data["total_bytes"] += size
            top = rel.split(os.sep, 1)[0]
            data["folder_sizes"][top] = data["folder_sizes"].get(top, 0) + size
            if size > data["large_file_threshold"]:
                data["large_files"].append({"path": rel, "size": size})
            if path.suffix in (".pyc", ".pyo"):
                data["pyc_files"].append(rel)
            if path.suffix == ".lock":
                data["lock_files"].append(rel)
            if _is_qm_file(path):
                data["qm_leftovers"].append(rel)


def _human(nbytes: int) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if nbytes < 1024 or unit == "TiB":
            return f"{nbytes:.1f} {unit}" if unit != "B" else f"{nbytes} B"
        nbytes /= 1024.0
    return f"{nbytes} B"


def format_text(data: dict) -> str:
    """Render :func:`collect` output as a compact worst-first listing."""
    lines = [f"Junk report for {data['root']} "
             f"({data['n_files']} files, {data['n_dirs']} dirs, "
             f"total {_human(data['total_bytes'])})"]

    def section(title: str, items: list[str], note: str | None = None) -> None:
        if not items:
            return
        lines.append(f"\n{title} ({len(items)}):")
        for item in items[:50]:
            lines.append(f"  {item}")
        if len(items) > 50:
            lines.append(f"  ... and {len(items) - 50} more")
        if note:
            lines.append(f"  ({note})")

    section("Large files (> threshold)",
            [f"{f['path']}  {_human(f['size'])}" for f in
             sorted(data["large_files"], key=lambda f: -f["size"])])
    section("QM leftovers (ORCA/xtb)", data["qm_leftovers"],
            "scratch files from quantum-chemistry runs; safe to archive/delete after review")
    section("Lock files", data["lock_files"])
    section("__pycache__ dirs", data["pycache_dirs"])
    section("Bytecode files", data["pyc_files"])

    if data["folder_sizes"]:
        lines.append("\nSize per top-level entry (worst first):")
        for name, size in sorted(data["folder_sizes"].items(),
                                 key=lambda kv: -kv[1]):
            lines.append(f"  {name:40s} {_human(size)}")

    if not any([data["large_files"], data["qm_leftovers"], data["lock_files"],
                data["pycache_dirs"], data["pyc_files"]]):
        lines.append("\nNothing suspicious found.")
    lines.append("\nRead-only report: nothing was modified or deleted.")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point; returns 0 always (this is a report, not a check)."""
    parser = argparse.ArgumentParser(
        prog="python -m delfin.agent.report_junk",
        description="Report workspace junk that git status hides. Read-only.")
    parser.add_argument("path", nargs="?", default=None,
                        help="workspace root to scan (default: cwd)")
    args = parser.parse_args(argv)
    data = collect(args.path)
    sys.stdout.write(format_text(data) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
