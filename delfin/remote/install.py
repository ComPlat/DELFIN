"""DELFIN on a cluster, one install per commit: the PC's own, published, side by side.

Input: a ``Connection`` (``connection.connect``), the git checkout the PC
runs DELFIN from, and the install sources (``clusters.install_sources``).
Output: an ``Outcome`` per call; ``CommitNotPublished`` and ``InstallFailed``
(delfin.remote.errors) when DELFIN will not or could not install.

Semantics:

* The commit is the one the PC runs: HEAD of the checkout ``delfin`` is
  imported from. Changes in its working tree do not travel; the outcome
  notes them.
* DELFIN installs only a published commit unasked: one reachable from a
  branch of an install source (``git ls-remote``, then ancestry in the local
  repository, after fetching the branch heads' objects without writing a
  ref). Reachability matters: GitHub serves any commit of a fork network
  by its hash, so "the server knows this hash" would also admit a commit
  pushed only to somebody's fork.
* The tree travels over the tunnel as ``git archive`` of that commit into
  ``<install_root>/<sha>/``, with ``.delfin-commit`` naming commit and
  source. The cluster needs no access to GitHub for DELFIN itself, and gets
  exactly the commit, not the PC's working tree.
* ``install_delfin.sh`` runs there detached (``--repo`` that directory,
  ``--prefix`` the cluster's software_dir, the standard profile), so it
  outlives the tunnel and the PC. It runs on the login node within the
  share a user has there (hpc_guard's LOGIN_CORE_BUDGET, DELFIN issue
  #640): ``INSTALL_CORES`` compiles and builds at a time, one BLAS/OpenMP
  thread, ``nice -n 19``. hpc_guard's watcher counts every process of the
  user's on the login node and stops (SIGSTOP) the busiest DELFIN one --
  the install's are, by their paths -- when the user stays over the
  budget, so the install leaves a core of it to the person's own work; an
  install stopped all the same is reported (``stopped``), not waited for
  forever. Each commit has its own venv; ORCA,
  OpenMPI and the tools are found or built once and shared. It writes
  ``~/.delfin_settings.json`` and ``~/.delfin_env.sh`` on the cluster as any
  install does.
* An install is complete when install_delfin.sh succeeded, or when its venv
  imports DELFIN and OpenMPI was set up although a tool failed (std2 is
  built with gfortran, which a login node may not have): what it lacks --
  the tools that failed, ORCA not found -- comes back in ``warnings``.
  A complete install is used as it is; a running one is waited for; a
  failed one is reported with the tail of its log and started again only
  with ``retry=True``.
* ``remove_unused`` deletes DELFIN's installs (directories with
  ``.delfin-commit``, nothing else) that are not kept and that no queued or
  running job names (``--comment=delfin-commit=<sha>``, see ``job_comment``).
* Every upload, start and removal is written to the audit log.
"""

from __future__ import annotations

import json
import subprocess
import tarfile
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from typing import Callable, Optional, Sequence

from . import clusters as _clusters
from .errors import CommitNotPublished, InstallFailed

SCRIPT = Path(__file__).with_name("scripts") / "installs.sh"

#: How long ``ensure_installed`` waits for a first install (OpenMPI is
#: compiled, the tools are fetched) before it reports the install as
#: still running; the install itself carries on.
DEFAULT_WAIT_SECONDS = 3 * 3600


def _install_cores() -> int:
    """Compiles and builds an install runs at a time on a login node: hpc_guard's
    budget for all of a user's processes there, less one core for the rest."""
    from delfin.hpc_guard import LOGIN_CORE_BUDGET
    return max(1, int(LOGIN_CORE_BUDGET) - 1)


#: See the module's semantics: one compile at a time with the budget of 2.
INSTALL_CORES = _install_cores()

#: What an install that hpc_guard stopped is told with.
STOPPED_HINT = ("a process of the install was stopped (SIGSTOP), as hpc_guard does to DELFIN's "
                "processes while a user keeps a login node over its core budget; it does "
                "not go on by itself, and the next submit of this commit starts it again")


@dataclass(frozen=True)
class Commit:
    """The commit the PC runs DELFIN from."""

    sha: str
    checkout: Path
    dirty: bool  # tracked files differ from the commit; those changes do not travel


@dataclass
class Outcome:
    """What ``ensure_installed`` found or did."""

    sha: str
    path: str        # on the cluster, "~/..." or absolute
    state: str       # "complete" (present or installed now) | "running" (still going)
    installed_now: bool = False
    source: str = ""  # the install source the commit was found on
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)  # what the install lacks


def job_comment(sha: str) -> str:
    """The ``sbatch --comment`` that marks a job as running on install ``sha``,
    so ``remove_unused`` keeps that install while the job is queued or runs."""
    return f"delfin-commit={sha}"


def pc_commit(checkout: Optional[Path] = None) -> Commit:
    """The commit of the checkout DELFIN runs from (or of ``checkout``)."""
    if checkout is None:
        import delfin
        checkout = Path(delfin.__file__).resolve().parents[1]
    checkout = Path(checkout)
    try:
        sha = _git(checkout, "rev-parse", "HEAD").strip()
        status = _git(checkout, "status", "--porcelain", "--untracked-files=no")
    except (OSError, subprocess.CalledProcessError) as exc:
        raise CommitNotPublished(
            f"DELFIN on this PC does not run from a git checkout ({checkout}), so there "
            f"is no commit to install on a cluster: {exc}") from exc
    return Commit(sha=sha, checkout=checkout, dirty=bool(status.strip()))


def published_source(commit: Commit, sources: Sequence[str]) -> Optional[str]:
    """The first install source on whose branches ``commit`` lies, or None."""
    for url in sources:
        try:
            listing = _git(commit.checkout, "ls-remote", "--heads", url, timeout=120)
        except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
            continue
        heads = sorted({line.split()[0] for line in listing.splitlines() if line.strip()})
        if commit.sha in heads:
            return url
        missing = [h for h in heads if not _has_commit(commit.checkout, h)]
        if missing:
            try:
                # Objects only: no ref and no FETCH_HEAD is written into the
                # user's repository.
                _git(commit.checkout, "fetch", "--quiet", "--no-tags",
                     "--no-write-fetch-head", url, *missing, timeout=600)
            except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
                pass
        for head in heads:
            if _has_commit(commit.checkout, head) and _is_ancestor(commit.checkout,
                                                                   commit.sha, head):
                return url
    return None


def state(conn, sha: str) -> str:
    """complete | running | stopped | failed | partial | absent, for install
    ``sha`` (stopped: a process of the running install was stopped)."""
    return _script(conn, "state", _install_dir(conn, sha)).strip() or "absent"


def ensure_installed(conn, *, commit: Optional[Commit] = None,
                     sources: Optional[Sequence[str]] = None,
                     settings: Optional[dict] = None, retry: bool = False,
                     wait: bool = True, wait_seconds: float = DEFAULT_WAIT_SECONDS,
                     poll_seconds: float = 20,
                     on_progress: Optional[Callable[[str], None]] = None) -> Outcome:
    """The PC's commit installed on ``conn``'s cluster; installing it if needed."""
    say = on_progress or (lambda message: None)
    commit = commit or pc_commit()
    path = _install_dir(conn, commit.sha)
    notes = []
    if commit.dirty:
        notes.append(f"the PC's checkout has changes beyond commit {commit.sha[:10]}; "
                     "the cluster runs the commit without them")
    current = state(conn, commit.sha)
    if current == "complete":
        return Outcome(commit.sha, path, "complete", notes=notes,
                       warnings=_warnings(conn, path))
    if current == "stopped":
        # hpc_guard stopped it, and on a login node nothing may resume it: a
        # submit the person approved starts it again, on the same cores. The
        # uploaded tree is there already.
        say(f"the install of {commit.sha[:10]} on {conn.name} was stopped (hpc_guard); "
            "starting it again")
        _script(conn, "abort", path)
        _script(conn, "start", path, conn.cluster.software_dir, str(INSTALL_CORES))
        _audit("install-restart", conn.name, commit.sha, "", path)
        notes.append("the install had been stopped by hpc_guard; it started again")
        current = "running"
    if current == "failed" and not retry:
        raise InstallFailed(f"the install of {commit.sha[:10]} on {conn.name} failed before; "
                            "fix the cause and install again with retry",
                            log=_script(conn, "log", path))
    source = ""
    if current != "running":
        if sources is None:
            sources = _clusters.install_sources(settings)
        say(f"checking that {commit.sha[:10]} is published ({', '.join(sources) or 'no source'})")
        source = published_source(commit, sources) or ""
        if not source:
            raise CommitNotPublished(
                f"commit {commit.sha[:10]} is on no branch of {', '.join(sources) or 'any source'}"
                "; push it there first, DELFIN installs only published commits on a cluster")
        say(f"uploading {commit.sha[:10]} to {conn.name}:{path}")
        _upload_tree(conn, commit, source, path)
        say(f"starting install_delfin.sh on {conn.name}")
        _script(conn, "start", path, conn.cluster.software_dir, str(INSTALL_CORES))
        _audit("install", conn.name, commit.sha, source, path)
    if not wait:
        return Outcome(commit.sha, path, "running", installed_now=True, source=source,
                       notes=notes)
    deadline = time.monotonic() + wait_seconds
    while True:
        current = state(conn, commit.sha)
        if current == "complete":
            say(f"{commit.sha[:10]} is installed on {conn.name}")
            return Outcome(commit.sha, path, "complete", installed_now=True, source=source,
                           notes=notes, warnings=_warnings(conn, path))
        if current in ("failed", "partial", "absent", "stopped"):
            why = f": {STOPPED_HINT}" if current == "stopped" else ""
            raise InstallFailed(f"the install of {commit.sha[:10]} on {conn.name} failed{why}",
                                log=_script(conn, "log", path))
        if time.monotonic() > deadline:
            return Outcome(commit.sha, path, "running", installed_now=True, source=source,
                           notes=notes + ["still installing; it carries on without the PC"])
        say(f"installing {commit.sha[:10]} on {conn.name} ...")
        time.sleep(poll_seconds)


def installs(conn) -> dict[str, str]:
    """DELFIN's installs on the cluster: commit -> state."""
    out = {}
    for line in _script(conn, "list", conn.cluster.install_root).splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[0]] = parts[1]
    return out


def log_tail(conn, sha: str) -> str:
    """The last lines of install ``sha``'s log on the cluster."""
    return _script(conn, "log", _install_dir(conn, sha))


def warnings(conn, sha: str) -> list[str]:
    """What complete install ``sha`` lacks (tools that failed, ORCA not
    found), one line each; empty for an install that is not complete."""
    return _warnings(conn, _install_dir(conn, sha))


def in_use(conn) -> set[str]:
    """Commits that a queued or running job of the user's names."""
    return {line.strip() for line in _script(conn, "inuse").splitlines() if line.strip()}


def remove_unused(conn, keep: Sequence[str] = ()) -> list[str]:
    """Remove DELFIN's installs that are neither kept nor used by a job; the
    commits removed. A running install is never removed."""
    used = in_use(conn)
    removed = []
    for sha, st in installs(conn).items():
        if st in ("running", "stopped") or any(sha.startswith(k) or k.startswith(sha)
                                               for k in keep):
            continue
        if any(sha.startswith(u) or u.startswith(sha) for u in used):
            continue
        path = _install_dir(conn, sha)
        _script(conn, "remove", path)
        _audit("remove", conn.name, sha, "", path)
        removed.append(sha)
    return removed


# --- helpers ------------------------------------------------------------------------


def _install_dir(conn, sha: str) -> str:
    return f"{conn.cluster.install_root.rstrip('/')}/{sha}"


def _warnings(conn, path: str) -> list[str]:
    return [line.strip() for line in _script(conn, "warnings", path).splitlines()
            if line.strip()]


def _relative_to_home(path: str) -> str:
    """A remote path as rsync and ssh read it: relative paths start in the
    cluster's home, and "~" in an argument is not expanded."""
    return path[2:] if path.startswith("~/") else path


def _script(conn, action: str, *args: str) -> str:
    result = conn.feed_script(SCRIPT, [action, *args], timeout=600)
    if not result.ok:
        raise InstallFailed(f"{action} on {conn.name} failed (exit {result.exit_code}): "
                            f"{(result.stderr or '').strip()}")
    return result.stdout or ""


def _upload_tree(conn, commit: Commit, source: str, path: str) -> None:
    from delfin.agent import safe_git
    made = safe_git.run(["archive", "--format=tar", commit.sha], commit.checkout,
                        timeout=300, text=False)
    if made.returncode != 0:
        raise InstallFailed(f"git archive of {commit.sha[:10]} failed: "
                            f"{(made.stderr or b'').decode('utf-8', 'replace').strip()}")
    archive = made.stdout
    with tempfile.TemporaryDirectory(prefix="delfin-install-") as tmp:
        tree = Path(tmp) / commit.sha
        tree.mkdir()
        with tarfile.open(fileobj=BytesIO(archive)) as tar:
            tar.extractall(tree, filter="data")
        (tree / ".delfin-commit").write_text(json.dumps({
            "commit": commit.sha, "source": source,
            "uploaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }) + "\n", encoding="utf-8")
        remote = _relative_to_home(path)
        parent = remote.rsplit("/", 1)[0] if "/" in remote else "."
        made = conn.run(["mkdir", "-p", parent])
        if not made.ok:
            raise InstallFailed(f"cannot create {parent} on {conn.name}: "
                                f"{(made.stderr or '').strip()}")
        result = conn.upload(f"{tree}/", f"{remote}/")
        if not result.ok:
            raise InstallFailed(f"uploading {commit.sha[:10]} to {conn.name} failed: "
                                f"{(result.stderr or '').strip()}")


# Git through safe_git: the checkout is a repository an agent's commands
# can write, and git runs programs its own configuration names (filters,
# fsmonitor, hooks, an ssh command) -- here with the user's rights, outside
# every cage. safe_git.run neutralises them and never raises.
def _git(checkout: Path, *args: str, timeout: float = 60) -> str:
    from delfin.agent import safe_git
    done = safe_git.run(list(args), checkout, timeout=timeout)
    if done.returncode != 0:
        raise subprocess.CalledProcessError(done.returncode, done.args, done.stdout, done.stderr)
    return done.stdout


def _has_commit(checkout: Path, sha: str) -> bool:
    from delfin.agent import safe_git
    return safe_git.run(["cat-file", "-e", f"{sha}^{{commit}}"], checkout).returncode == 0


def _is_ancestor(checkout: Path, sha: str, head: str) -> bool:
    from delfin.agent import safe_git
    return safe_git.run(["merge-base", "--is-ancestor", sha, head], checkout).returncode == 0


def _audit(action: str, cluster: str, sha: str, source: str, path: str) -> None:
    try:
        from delfin.agent import audit_log
        audit_log.append(audit_log.make_record(
            tool="cluster_install", decision="allow", mode=f"auto_{action}",
            path=path, reason=f"{action} DELFIN {sha[:10]} on {cluster}",
            extra={"cluster": cluster, "commit": sha, "source": source}))
    except Exception:
        pass  # the audit log never blocks the action it records (audit_log.append)


__all__ = ["Commit", "INSTALL_CORES", "Outcome", "SCRIPT", "STOPPED_HINT", "ensure_installed",
           "in_use", "installs", "job_comment", "log_tail", "pc_commit", "published_source",
           "remove_unused", "state", "warnings"]
