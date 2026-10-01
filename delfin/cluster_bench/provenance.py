"""What a batch run was built with: code state, interpreters, package versions, settings.

Recorded in the manifest when the run is prepared and again in every chunk; ``run-shard`` refuses
a chunk whose code or environment differs from the manifest, and ``collect`` reports chunks that
disagree with each other.  A run whose chunks were built with different code is not one result.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from decimal import ROUND_CEILING, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKERS_DIR = HERE / "tool_workers"
REPO_ROOT = HERE.parent.parent

#: Tools, their worker and their defaults.  ``python`` = None: the DELFIN interpreter itself.
TOOLS = {
    "manta": {"worker": None, "mode": "champion", "workers": 36, "threads": 7, "mem": "72G",
              "shard_size": 500, "needs_specs": False},
    "architector": {"worker": "bench_worker.py", "mode": "full", "workers": 48, "threads": 1,
                    "mem": "80G", "shard_size": 250, "needs_specs": True},
    "molsimplify": {"worker": "bench_worker.py", "mode": "full", "workers": 48, "threads": 1,
                    "mem": "80G", "shard_size": 250, "needs_specs": True},
    "mace": {"worker": "mace_worker.py", "mode": "paper", "workers": 48, "threads": 1,
             "mem": "80G", "shard_size": 250, "needs_specs": True},
}

MODES = {"manta": ("champion", "builder"), "architector": ("full", "default"),
         "molsimplify": ("full",), "mace": ("paper", "extended")}


def cbatch_sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def cbatch_effective_timeout(base_s, speed_factor) -> int:
    """ceil(base x factor) in exact decimal arithmetic (21600 x 1.3 = 28080, not 28081)."""
    return int((Decimal(str(base_s)) * Decimal(str(speed_factor))).to_integral_value(ROUND_CEILING))


def _git(*args) -> str | None:
    try:
        done = subprocess.run(["git", "-C", str(REPO_ROOT), *args], capture_output=True,
                              text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return done.stdout if done.returncode == 0 else None


def cbatch_code_sha256() -> str:
    """sha256 over path and content of every delfin/**/*.py -- the code that builds, whether or
    not git is available on the node."""
    h = hashlib.sha256()
    root = REPO_ROOT / "delfin"
    for p in sorted(root.rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        h.update(str(p.relative_to(root)).encode() + b"\0")
        h.update(hashlib.sha256(p.read_bytes()).digest())
    return h.hexdigest()


def _head_commit_from_files() -> str | None:
    """HEAD without the git binary (compute nodes may not have it)."""
    try:
        git = REPO_ROOT / ".git"
        if git.is_file():                       # worktree: "gitdir: <path>"
            git = Path(git.read_text().split(":", 1)[1].strip())
        head = (git / "HEAD").read_text().strip()
        if not head.startswith("ref:"):
            return head
        ref = head.split(":", 1)[1].strip()
        for base in (git, Path((git / "commondir").read_text().strip()) if (git / "commondir").exists()
                     else git):
            base = base if base.is_absolute() else (git / base).resolve()
            if (base / ref).exists():
                return (base / ref).read_text().strip()
            packed = base / "packed-refs"
            if packed.exists():
                for ln in packed.read_text().splitlines():
                    if ln.endswith(" " + ref):
                        return ln.split()[0]
    except (OSError, IndexError):
        pass
    return None


def cbatch_code_state() -> dict:
    """Commit of the DELFIN checkout, whether delfin/ has uncommitted changes, and the content
    hash of delfin/**/*.py (the value a chunk must agree on)."""
    commit = (_git("rev-parse", "HEAD") or "").strip() or _head_commit_from_files()
    diff = _git("diff", "HEAD", "--", "delfin")
    untracked = _git("ls-files", "--others", "--exclude-standard", "--", "delfin")
    dirty = None if diff is None else bool(diff.strip() or (untracked or "").strip())
    return {"repo": str(REPO_ROOT), "commit": commit, "dirty": dirty,
            "code_sha256": cbatch_code_sha256()}


def cbatch_probe_env(python: str) -> dict:
    """Run tool_workers/env_probe.py with ``python``; {"error": ...} if that interpreter fails."""
    try:
        done = subprocess.run([python, str(WORKERS_DIR / "env_probe.py")], capture_output=True,
                              text=True, timeout=300, env=cbatch_tool_child_env({}))
    except (OSError, subprocess.SubprocessError) as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}
    if done.returncode != 0:
        return {"error": (done.stderr or "")[-500:]}
    try:
        return json.loads(done.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return {"error": "no JSON from env_probe: " + (done.stdout or "")[-300:]}


def cbatch_tool_child_env(base) -> dict:
    """Environment of every tool process: the caller's, without anything that could make the
    tool's interpreter load foreign modules or libraries, single-threaded, fixed hash seed."""
    import os

    env = dict(os.environ)
    env.update(base)
    for k in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP", "LD_LIBRARY_PATH"):
        env.pop(k, None)
    # PYTHONNOUSERSITE is NOT forced here: the sbatch script exports it (a user site in $HOME
    # must not leak into a cluster build), while an interpreter that keeps packages in its user
    # site (a local reference env) must see them -- otherwise every build would fail on import.
    env["PYTHONHASHSEED"] = "0"
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[v] = "1"
    return env


def cbatch_construction_env(config: str) -> dict:
    from delfin.cli_manta import construction_env

    return dict(sorted(construction_env(config, environ={}).items()))


def cbatch_provenance(tool: str, tool_python: str, mode: str) -> dict:
    """Everything a chunk must agree on with the manifest."""
    prov = {"code": cbatch_code_state(), "delfin_env": cbatch_probe_env(sys.executable)}
    if tool == "manta":
        prov["construction_env"] = cbatch_construction_env(mode)
    else:
        prov["tool_env"] = cbatch_probe_env(tool_python)
        prov["worker_sha256"] = {
            p.name: cbatch_sha256_file(p) for p in sorted(WORKERS_DIR.glob("*.py"))
            if p.name != "__init__.py"}
    if tool == "manta":
        prov["worker_sha256"] = {"manta_child.py": cbatch_sha256_file(HERE / "manta_child.py")}
    return prov


def cbatch_provenance_mismatch(want: dict, have: dict) -> list:
    """Keys whose values differ.  For the code: the content hash, and the commit where both
    sides know it (git may be missing on a compute node)."""
    out = []
    for key in sorted(set(want) | set(have)):
        a, b = want.get(key), have.get(key)
        if key == "code" and isinstance(a, dict) and isinstance(b, dict):
            same = a.get("code_sha256") == b.get("code_sha256") and (
                a.get("commit") is None or b.get("commit") is None or a.get("commit") == b.get("commit"))
            if not same:
                out.append(key)
        elif a != b:
            out.append(key)
    return out
