#!/usr/bin/env python3
"""Shared acceptance machinery for the generalist tasks. Not an
acceptance script on its own.

A ``generalist_*_verify.py`` script builds a ``Report``, records named
checks on the repository the setup built, and ends with
``sys.exit(report.finish())``. Every check reads the END STATE: test
runs, file contents, git history. Nothing reads what the model said.

``Report.check(label, ok, required=True, weight=1)``:
  required checks decide pass/fail (exit 0 only when all hold);
  every check, required or not, counts towards the quality,
  ``QUALITY = round(100 * held weight / total weight)``, printed as
  ``QUALITY: n`` for the scorer (benchmark.acceptance_quality).

Output, one line per check (``PASS``/``FAIL`` label: detail), then the
QUALITY line, then ``OK`` or ``FAIL: <first required failure>``.

Hidden tests are run on a scratch copy of the work tree (``.git``,
caches and ``node_modules`` left out), so they neither reach nor leave
anything in the repository being judged.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BASE_FILE = "delfin-bench-base"
_MARKER = re.compile(r"^(<<<<<<< |=======\s*$|>>>>>>> )", re.M)


TOOLCHAIN_EXIT = 3


def toolchain_bin() -> str:
    """Same lookup as setup/generalist_repo.toolchain_bin."""
    raw = os.environ.get("DELFIN_TOOLCHAIN_BIN") or os.path.expanduser(
        "~/micromamba/envs/delfin-toolchains/bin")
    return raw if raw and os.path.isdir(raw) else ""


def _env() -> dict:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["GIT_TERMINAL_PROMPT"] = "0"
    tc = toolchain_bin()
    if tc:
        env["PATH"] = tc + os.pathsep + env.get("PATH", "")
    env.setdefault("GOTOOLCHAIN", "local")
    env.setdefault("GOPROXY", "off")
    env.setdefault("CARGO_NET_OFFLINE", "true")
    return env


def which(tool: str) -> str:
    return shutil.which(tool, path=_env().get("PATH")) or ""


def toolchain_missing(what: str) -> int:
    """The verdict for a host without the task's toolchain: not run.
    benchmark_runner.run_acceptance reads exit 3 with this phrase as
    "not run", never as a failure."""
    print(f"toolchain unavailable: {what}", file=sys.stderr)
    return TOOLCHAIN_EXIT


def tool_test(repo: Path, argv: list[str],
              extra: dict[str, str] | None = None,
              timeout: int = 600) -> tuple[bool, str]:
    """``argv`` (cargo test / go test) on a scratch copy with ``extra``."""
    work = _scratch(repo, extra)
    try:
        proc = subprocess.run(argv, cwd=str(work), capture_output=True,
                              text=True, timeout=timeout, env=_env())
        lines = (proc.stdout + proc.stderr).strip().splitlines()
        keep = [ln for ln in lines if ln.startswith(
            ("test result:", "ok ", "FAIL", "--- FAIL", "error"))]
        return proc.returncode == 0, " | ".join((keep or lines[-3:])[:6])
    except subprocess.TimeoutExpired:
        return False, f"{argv[0]} timed out"
    finally:
        shutil.rmtree(work.parent, ignore_errors=True)


class Report:
    def __init__(self, name: str) -> None:
        self.name = name
        self.rows: list[tuple[str, bool, bool, float, str]] = []

    def check(self, label: str, ok: bool, *, required: bool = True,
              weight: float = 1.0, detail: str = "") -> bool:
        self.rows.append((label, bool(ok), required, weight,
                          " ".join(str(detail).split())[:300]))
        return bool(ok)

    def quality(self) -> int:
        total = sum(w for _, _, _, w, _ in self.rows)
        held = sum(w for _, ok, _, w, _ in self.rows if ok)
        return int(round(100 * held / total)) if total else 0

    def finish(self) -> int:
        for label, ok, req, _, detail in self.rows:
            tag = "PASS" if ok else ("FAIL" if req else "MISS")
            print(f"{tag} {label}" + (f": {detail}" if detail else ""))
        print(f"QUALITY: {self.quality()}")
        failed = [r for r in self.rows if r[2] and not r[1]]
        if failed:
            print(f"FAIL: {failed[0][0]}")
            return 1
        print("OK")
        return 0


def workspace_repo(name: str) -> Path:
    """The repository ``name`` under argv[1]; exits 2 without an argument."""
    if len(sys.argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        sys.exit(2)
    return Path(sys.argv[1]).resolve() / name


def git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(repo), capture_output=True,
                          text=True, env=_env())


def out(repo: Path, *args: str) -> str:
    return git(repo, *args).stdout.strip()


def base(repo: Path) -> str:
    try:
        return (repo / ".git" / BASE_FILE).read_text(encoding="utf-8").strip()
    except OSError:
        return ""


def changed_paths(repo: Path, since: str) -> set[str]:
    """Every path whose work-tree content differs from ``since``: tracked
    changes (committed or not) plus untracked, non-ignored files."""
    paths = set(out(repo, "diff", "--name-only", since).splitlines())
    paths |= set(out(repo, "ls-files", "--others",
                     "--exclude-standard").splitlines())
    return {p for p in paths if p}


def blob_at(repo: Path, rev: str, rel: str) -> str | None:
    proc = git(repo, "show", f"{rev}:{rel}")
    return proc.stdout if proc.returncode == 0 else None


def read(repo: Path, rel: str) -> str:
    try:
        return (repo / rel).read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return ""


def files_with_markers(repo: Path) -> list[str]:
    hits = []
    for rel in out(repo, "ls-files").splitlines():
        if _MARKER.search(read(repo, rel)):
            hits.append(rel)
    return hits


def tree_text(repo: Path, suffixes: tuple[str, ...]) -> dict[str, str]:
    """{relpath: text} for tracked and untracked work-tree files."""
    rels = set(out(repo, "ls-files").splitlines())
    rels |= set(out(repo, "ls-files", "--others",
                    "--exclude-standard").splitlines())
    return {r: read(repo, r) for r in sorted(rels)
            if r.endswith(suffixes) and (repo / r).is_file()}


def _scratch(repo: Path, extra: dict[str, str] | None) -> Path:
    tmp = Path(tempfile.mkdtemp(prefix="generalist-acc-"))
    dst = tmp / repo.name
    shutil.copytree(repo, dst, ignore=shutil.ignore_patterns(
        ".git", "__pycache__", ".pytest_cache", "node_modules", "target"))
    for rel, text in (extra or {}).items():
        p = dst / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return dst


def pytest(repo: Path, *args: str,
           extra: dict[str, str] | None = None) -> tuple[bool, str]:
    """pytest on a scratch copy with ``extra`` files added."""
    work = _scratch(repo, extra)
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
             *args], cwd=str(work), capture_output=True, text=True,
            timeout=180, env=_env())
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        return proc.returncode == 0, " | ".join(tail)
    except subprocess.TimeoutExpired:
        return False, "pytest timed out"
    finally:
        shutil.rmtree(work.parent, ignore_errors=True)


def node_prefix() -> list[str] | None:
    node = which("node")
    if not node:
        return None
    with tempfile.TemporaryDirectory(prefix="generalist-node-") as tmp:
        probe = Path(tmp) / "p.ts"
        probe.write_text("const n: number = 1;\nconsole.log(n + 1);\n",
                         encoding="utf-8")
        for prefix in ([node], [node, "--experimental-strip-types"]):
            proc = subprocess.run([*prefix, "--no-warnings", str(probe)],
                                  capture_output=True, text=True, timeout=60)
            if proc.returncode == 0 and proc.stdout.strip() == "2":
                return prefix
    return None


def node_test(repo: Path, files: list[str],
              extra: dict[str, str] | None = None) -> tuple[bool, str]:
    prefix = node_prefix()
    if prefix is None:
        return False, "no node that runs TypeScript"
    work = _scratch(repo, extra)
    try:
        proc = subprocess.run([*prefix, "--no-warnings", "--test", *files],
                              cwd=str(work), capture_output=True, text=True,
                              timeout=180)
        lines = (proc.stdout + proc.stderr).strip().splitlines()
        summary = [ln for ln in lines if ln.startswith(("# pass", "# fail"))]
        return proc.returncode == 0, " | ".join(summary or lines[-3:])
    except subprocess.TimeoutExpired:
        return False, "node --test timed out"
    finally:
        shutil.rmtree(work.parent, ignore_errors=True)


if __name__ == "__main__":
    print("not an acceptance script on its own; see the docstring",
          file=sys.stderr)
    sys.exit(2)
