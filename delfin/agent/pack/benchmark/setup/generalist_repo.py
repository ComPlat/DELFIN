#!/usr/bin/env python3
"""Shared builder for the generalist setups. Not a setup on its own.

Each ``generalist_*.py`` setup calls ``main(build)`` with a function
``build(workspace) -> Path`` that creates one git repository under the
workspace and returns it. This module supplies:

* ``new_repo(workspace, name)``: refuses an existing directory (a
  workspace that was not restored is reported, never half-rebuilt),
  ``git init -b main`` with a repository-local identity, no signing;
* ``write`` / ``commit``: files from a ``{relpath: text}`` dict;
* ``mark_base(repo, sha)``: records the commit the acceptance script
  compares against, in ``.git/delfin-bench-base`` (untracked, outside
  the work tree);
* ``pytest_rc`` / ``node_rc``: the preconditions a setup asserts on its
  own fixture (the given test is red, the rest green, ...);
* ``node_toolchain()``: a Node that runs ``.ts`` directly (type stripping,
  Node >= 22.6 with the flag, on by default from 23.6); ``require(tool)``
  for cargo and go. Tools are looked up with ``toolchain_bin()`` first on
  PATH. A host without one makes the setup exit 3 with ``toolchain
  unavailable:``, which the runner books as unmeasured, not as a model
  failure;
* every child runs offline: GOTOOLCHAIN=local, GOPROXY=off,
  CARGO_NET_OFFLINE=true (the fixtures have no dependencies).

A failed precondition exits 1 with the reason on stderr.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable

TOOLCHAIN_EXIT = 3
BASE_FILE = "delfin-bench-base"

GITIGNORE_PY = "__pycache__/\n*.pyc\n.pytest_cache/\n*.egg-info/\n.venv/\n"
GITIGNORE_NODE = "node_modules/\n.cache/\n"


class Precondition(RuntimeError):
    """The fixture is not what the task needs; exit 1."""


class ToolchainMissing(RuntimeError):
    """The host lacks a tool the task needs; exit 3."""


def toolchain_bin() -> str:
    """The directory holding node/cargo/go for these tasks, or "".

    ``DELFIN_TOOLCHAIN_BIN`` when set (the runner exports it), else
    ``~/micromamba/envs/delfin-toolchains/bin`` when that exists. Only a
    PATH prefix: a tool missing there is still looked up on PATH."""
    raw = os.environ.get("DELFIN_TOOLCHAIN_BIN") or os.path.expanduser(
        "~/micromamba/envs/delfin-toolchains/bin")
    return raw if raw and os.path.isdir(raw) else ""


def _quiet_env() -> dict:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["GIT_TERMINAL_PROMPT"] = "0"
    tc = toolchain_bin()
    if tc:
        env["PATH"] = tc + os.pathsep + env.get("PATH", "")
    # Offline and pinned: no toolchain download, no module proxy, no
    # crates index. The fixtures have no dependencies.
    env.setdefault("GOTOOLCHAIN", "local")
    env.setdefault("GOPROXY", "off")
    env.setdefault("CARGO_NET_OFFLINE", "true")
    return env


def which(tool: str) -> str:
    return shutil.which(tool, path=_quiet_env().get("PATH")) or ""


def tool_rc(repo: Path, argv: list[str],
            timeout: int = 300) -> tuple[int, str]:
    """``argv`` (cargo/go) on a scratch copy of ``repo``."""
    work = _copy(repo)
    try:
        proc = subprocess.run(argv, cwd=str(work), capture_output=True,
                              text=True, timeout=timeout, env=_quiet_env())
        return proc.returncode, proc.stdout + proc.stderr
    finally:
        shutil.rmtree(work.parent, ignore_errors=True)


def require(tool: str) -> str:
    path = which(tool)
    if not path:
        raise ToolchainMissing(f"{tool} not on PATH")
    return path


def git(repo: Path, *args: str) -> str:
    proc = subprocess.run(["git", *args], cwd=str(repo), capture_output=True,
                          text=True, env=_quiet_env())
    if proc.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed: {proc.stderr.strip()[:300]}")
    return proc.stdout.strip()


def new_repo(workspace: str | Path, name: str) -> Path:
    repo = Path(workspace).resolve() / name
    if repo.exists():
        print(f"refusing to overwrite {repo}", file=sys.stderr)
        sys.exit(1)
    repo.mkdir(parents=True)
    git(repo, "init", "-q", "-b", "main")
    git(repo, "config", "user.email", "dev@example.invalid")
    git(repo, "config", "user.name", "Dev")
    git(repo, "config", "commit.gpgsign", "false")
    return repo


def write(repo: Path, files: dict[str, str]) -> None:
    for rel, text in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")


def commit(repo: Path, message: str) -> str:
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", message)
    return git(repo, "rev-parse", "HEAD")


def mark_base(repo: Path, sha: str = "") -> None:
    sha = sha or git(repo, "rev-parse", "HEAD")
    (repo / ".git" / BASE_FILE).write_text(sha + "\n", encoding="utf-8")


def _copy(repo: Path) -> Path:
    tmp = Path(tempfile.mkdtemp(prefix="generalist-pre-"))
    dst = tmp / repo.name
    shutil.copytree(repo, dst, ignore=shutil.ignore_patterns(
        ".git", "__pycache__", ".pytest_cache", "node_modules", "target"))
    return dst


def pytest_rc(repo: Path, *args: str) -> tuple[int, str]:
    """pytest on a scratch copy, so the fixture gets no cache files."""
    work = _copy(repo)
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
             *args], cwd=str(work), capture_output=True, text=True,
            timeout=120, env=_quiet_env())
        return proc.returncode, proc.stdout + proc.stderr
    finally:
        shutil.rmtree(work.parent, ignore_errors=True)


def node_toolchain() -> list[str]:
    """argv prefix that runs a .ts file with Node's type stripping."""
    node = require("node")
    with tempfile.TemporaryDirectory(prefix="generalist-node-") as tmp:
        probe = Path(tmp) / "probe.ts"
        probe.write_text("const n: number = 1;\nconsole.log(n + 1);\n",
                         encoding="utf-8")
        for prefix in ([node], [node, "--experimental-strip-types"]):
            proc = subprocess.run([*prefix, "--no-warnings", str(probe)],
                                  capture_output=True, text=True, timeout=60)
            if proc.returncode == 0 and proc.stdout.strip() == "2":
                return prefix
    raise ToolchainMissing(
        "node cannot run TypeScript directly (needs Node >= 22.6)")


def node_rc(repo: Path, test_files: list[str]) -> tuple[int, str]:
    prefix = node_toolchain()
    work = _copy(repo)
    try:
        proc = subprocess.run(
            [*prefix, "--no-warnings", "--test", *test_files],
            cwd=str(work), capture_output=True, text=True, timeout=120)
        return proc.returncode, proc.stdout + proc.stderr
    finally:
        shutil.rmtree(work.parent, ignore_errors=True)


def expect(label: str, ok: bool, detail: str = "") -> None:
    if not ok:
        raise Precondition(f"{label}: {detail[-600:]}")


def main(build: Callable[[str], Path]) -> None:
    if len(sys.argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        sys.exit(2)
    try:
        repo = build(sys.argv[1])
    except ToolchainMissing as exc:
        print(f"toolchain unavailable: {exc}", file=sys.stderr)
        sys.exit(TOOLCHAIN_EXIT)
    except Precondition as exc:
        print(f"precondition not met: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"prepared {repo}")
    sys.exit(0)


if __name__ == "__main__":
    print("not a setup on its own; see the module docstring", file=sys.stderr)
    sys.exit(2)
