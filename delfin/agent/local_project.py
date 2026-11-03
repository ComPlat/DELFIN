"""A workspace folder made into a local git repository, on the user's yes.

Worktrees, branches, sub-agents in their own worktrees and worktree_merge
all need git. A person's default workspace (a plain folder such as the
``agent_workspace`` each person gets on a shared machine) has none, and
``enter_worktree`` answered git's "not a git repository".

Input: a workspace directory. Output, in three steps:

* :func:`refusal` -- why the folder may not become a repository, or None.
  Refused: the home directory and its ancestors and the system roots
  (``api_client._is_forbidden_workspace_root``, the same floor a workspace
  has), any path under a system directory or a key directory
  (``_NON_GRANTABLE_READ_DIRS``, ``_SECRET_DIR_NAMES``), a mount point
  (``sandbox._is_tree_walk_root``: a whole file system is not a project),
  a folder inside another repository (that repository is the one to use),
  a missing directory.
* :func:`plan` -- what a first commit would hold and what stays out, read
  with git's own ignore matching against a scratch repository, so the
  question put to the user names real counts before anything is created.
* :func:`init_local_project` -- ``git init``, ``main`` as the branch, the
  default ``.gitignore`` (merged into an existing one, never replacing
  it), the files kept out by size or by content added to it by name, and
  one first commit. No remote is added: adding one is the user's step.
  Every step that fails after ``git init`` rolls back: the ``.git``
  directory created here is removed and ``.gitignore`` restored, so a
  failure leaves the folder as it was.

What stays out of the commit:

* QM program output and its binary artifacts by name (``*.out``,
  ``*.gbw``, ``*.hess``, ...) and, for a file whose name says nothing, by
  content: the ORCA line signatures of ``scripts/license_guard.py``
  (``QM_OUTPUT_SIGNATURES``, pinned to that list by a test). ORCA output
  may not be redistributed (LICENSE_GUARD.md); a local repository is the
  first step to a pushed one, so it is kept out here already.
* checkpoints and model weights, archives, files above
  :data:`SIZE_LIMIT_BYTES`, nested repositories;
* secrets: ``.env``, ``credentials*``, ``*.key``, ``*.pem``, ``id_*``.

After ``git add -A`` the staged list is checked again against the same
rules; a staged file that matches one (an existing ``.gitignore`` with a
``!`` negation re-including it) aborts the commit and rolls back.

Rejected: a silent ``git init`` on the first worktree request (a
repository in a person's folder is a change they must have agreed to);
asking in every mode except bypass (the question is about the folder, not
about a tool call, so no mode answers it); adding an ``origin`` (no URL is
known, and a push target is not the agent's choice); a ``.gitignore``
size rule (git has none, so files over the limit are listed by name).
"""

from __future__ import annotations

import os
import re
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from . import safe_git

#: Files above this size stay out of the first commit (listed by name in
#: .gitignore). GitHub refuses files above 100 MB and warns above 50 MB;
#: a project's own sources are far below either.
SIZE_LIMIT_BYTES = 10 * 1024 * 1024

#: More candidate files than this is not a project folder: refused rather
#: than committed.
MAX_FILES = 20000

#: Bytes of a file read for the content check.
_SNIFF_BYTES = 64 * 1024

#: The first line of the block this module writes into .gitignore.
GITIGNORE_MARKER = "# DELFIN local project: kept out of version control"

DEFAULT_GITIGNORE = f"""{GITIGNORE_MARKER}
# QM program output and artifacts. ORCA output may not be redistributed
# (see DELFIN's LICENSE_GUARD.md): it stays out of every commit.
*.out
*.gbw
*.hess
*.densities
*.densitiesinfo
*.tmp
*.tmp.*
*.engrad
*.opt
*.property.txt
*.property.json
*.bibtex
*.cis
*.loc
*.qro
*.uno
*.unso
*.finalensemble*.xyz
*.globalminimum*.xyz
*.npz
# Checkpoints and model weights
*.chk
*.fchk
*.ckpt
*.pt
*.pth
*.safetensors
*.onnx
checkpoints/
# Archives
*.zip
*.tar
*.tar.gz
*.tgz
*.tar.bz2
*.tar.xz
*.7z
*.rar
# Secrets
.env
.env.*
credentials*
*.key
*.pem
*.p12
id_*
# Caches and DELFIN's own state
.delfin/
__pycache__/
*.pyc
.ipynb_checkpoints/
"""

#: The patterns above as regular expressions over a path's base name, used
#: to check the staged list after ``git add`` (a second, independent check:
#: the commit is aborted when one matches).
_STAGED_FORBIDDEN_RE = re.compile(
    r"(?:\.(?:out|gbw|hess|densities|densitiesinfo|engrad|opt|cis|loc|qro|"
    r"uno|unso|npz|chk|fchk|ckpt|pt|pth|safetensors|onnx|zip|tar|tgz|7z|rar|"
    r"key|pem|p12)$|\.tmp(?:\..*)?$|\.property\.(?:txt|json)$|"
    r"\.tar\.(?:gz|bz2|xz)$|^\.env(?:\..*)?$|^credentials|^id_)",
    re.IGNORECASE)

#: ORCA output recognised by its lines. Kept equal to
#: ``scripts/license_guard.QM_OUTPUT_SIGNATURES`` (a test pins the pair):
#: that script is not shipped in an installed DELFIN, this module is.
QM_OUTPUT_SIGNATURES = (
    r"^\s*\* O\s+R\s+C\s+A\s+\*$",
    r"ORCA\s+TERMINATED\s+NORMALLY",
    r"^\s*ORCA\s+terminated\s+normally",
    r"Program\s+Version\s+6\.\d+\.\d+\s+-\s+RELEASE",
    r"^\s*\* Host name:",
    r"^\s*\* Working dir\.:",
    r"^\|\s*\d+>",
    r"THERMOCHEMISTRY\s+AT\s+\d+\.\d+K",
)
_QM_OUTPUT_RE = re.compile("(" + "|".join(QM_OUTPUT_SIGNATURES) + ")",
                           re.IGNORECASE | re.MULTILINE)
QM_SYNTHETIC_MARKER = (
    "# SYNTHETIC test fixture written for DELFIN -- not ORCA output")

#: Who commits when git knows nobody (a fresh account on a shared machine).
_FALLBACK_IDENTITY = {"GIT_AUTHOR_NAME": "DELFIN local project",
                      "GIT_AUTHOR_EMAIL": "delfin@localhost",
                      "GIT_COMMITTER_NAME": "DELFIN local project",
                      "GIT_COMMITTER_EMAIL": "delfin@localhost"}

FIRST_COMMIT_MESSAGE = (
    "Local project: first commit of the workspace\n\n"
    "Made by DELFIN on the user's consent. No remote is configured; "
    "adding one is the user's step.")


class LocalProjectError(RuntimeError):
    """The folder could not be made a repository; nothing was left behind."""


@dataclass
class Plan:
    """What a first commit of *workspace* would hold."""
    workspace: Path
    files: list = field(default_factory=list)      # committed, relative
    too_large: list = field(default_factory=list)  # above SIZE_LIMIT_BYTES
    qm_output: list = field(default_factory=list)  # ORCA output by content
    nested_repos: list = field(default_factory=list)
    ignored_count: int = 0                         # kept out by .gitignore

    @property
    def kept_out(self) -> list:
        return self.too_large + self.qm_output + self.nested_repos

    def summary(self) -> str:
        parts = [f"{len(self.files)} file(s) in the first commit"]
        out = []
        if self.ignored_count:
            out.append(f"{self.ignored_count} by the .gitignore rules "
                       "(QM output, weights, archives, secrets)")
        if self.too_large:
            out.append(f"{len(self.too_large)} above "
                       f"{SIZE_LIMIT_BYTES // (1024 * 1024)} MB")
        if self.qm_output:
            out.append(f"{len(self.qm_output)} recognised as ORCA output")
        if self.nested_repos:
            out.append(f"{len(self.nested_repos)} nested repositor(y/ies)")
        if out:
            parts.append("kept out: " + "; ".join(out))
        return ", ".join(parts)


@dataclass
class InitResult:
    workspace: Path
    commit: str
    branch: str
    plan: Plan


# -- refusal ------------------------------------------------------------------

def _resolved(path) -> Optional[Path]:
    try:
        return Path(path).expanduser().resolve()
    except Exception:
        return None


def enclosing_repository(path) -> str:
    """The top level of the repository *path* lies in, or ""."""
    r = safe_git.run(["rev-parse", "--show-toplevel"], path, timeout=10)
    return r.stdout.strip() if r.returncode == 0 else ""


def refusal(path) -> Optional[str]:
    """Why *path* may not become a local repository, or None."""
    from .api_client import (_is_forbidden_workspace_root,
                             _NON_GRANTABLE_READ_DIRS, _SECRET_DIR_NAMES)
    from .sandbox import _is_tree_walk_root
    p = _resolved(path)
    if p is None or not p.is_dir():
        return f"{path} is not an existing directory"
    if _is_forbidden_workspace_root(p):
        return (f"{p} is the home directory, an ancestor of it or a system "
                "root; a repository there would take in every project and "
                "key below it. Use a project folder inside it.")
    if str(p) in _NON_GRANTABLE_READ_DIRS or (
            len(p.parts) >= 2 and "/" + p.parts[1] in _NON_GRANTABLE_READ_DIRS):
        return f"{p} is under a system directory"
    if set(p.parts) & _SECRET_DIR_NAMES:
        return f"{p} is inside a key directory"
    if _is_tree_walk_root(p):
        return (f"{p} is the root of a file system (a mount point), not a "
                "project folder")
    top = enclosing_repository(p)
    if top:
        return (f"{p} is inside the repository {top}; use that repository "
                "instead of making a second one inside it")
    return None


# -- plan ---------------------------------------------------------------------

def _gitignore_text(workspace: Path) -> tuple[str, Optional[str]]:
    """(the .gitignore to write, the existing text or None)."""
    target = workspace / ".gitignore"
    try:
        old = target.read_text(encoding="utf-8")
    except FileNotFoundError:
        return DEFAULT_GITIGNORE, None
    except (OSError, UnicodeDecodeError) as exc:
        raise LocalProjectError(f"the existing .gitignore cannot be read: "
                                f"{exc}") from None
    if GITIGNORE_MARKER in old:
        return old, old
    sep = "" if old.endswith("\n") or not old else "\n"
    return old + sep + "\n" + DEFAULT_GITIGNORE, old


def _ignore_line(rel: str) -> str:
    """A .gitignore line matching exactly *rel* from the top level."""
    escaped = re.sub(r"([*?\[\]\\!# ])", r"\\\1", rel)
    return "/" + escaped


def looks_like_qm_output(path: Path) -> bool:
    """True when the start of *path* carries an ORCA output line."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(_SNIFF_BYTES)
    except OSError:
        return False
    if b"\0" in head[:4096]:
        return False
    text = head.decode("utf-8", "replace")
    if text.startswith(QM_SYNTHETIC_MARKER):
        return False
    return bool(_QM_OUTPUT_RE.search(text))


def _list_untracked(workspace: Path, git_dir: Path,
                    exclude_file: Optional[Path]) -> tuple[list, int]:
    """(not-ignored untracked paths, count of ignored ones), NUL-split."""
    env = {"GIT_DIR": str(git_dir), "GIT_WORK_TREE": str(workspace)}
    base = ["ls-files", "-z", "--others", "--exclude-standard"]
    if exclude_file is not None:
        base.append(f"--exclude-from={exclude_file}")
    r = safe_git.run(base, workspace, env=env, timeout=60)
    if r.returncode != 0:
        raise LocalProjectError(f"git ls-files failed: {r.stderr.strip()}")
    files = [n for n in r.stdout.split("\0") if n]
    ign = safe_git.run(base + ["--ignored", "--directory"], workspace,
                       env=env, timeout=60)
    ignored = len([n for n in ign.stdout.split("\0") if n]) \
        if ign.returncode == 0 else 0
    return files, ignored


def _classify(workspace: Path, names: list, plan: Plan) -> None:
    if len(names) > MAX_FILES:
        raise LocalProjectError(
            f"{len(names)} files would go into the first commit (more than "
            f"{MAX_FILES}); that is more than a project folder. Make the "
            "repository in the project's own folder instead.")
    for rel in names:
        if rel.endswith("/"):           # an embedded repository
            plan.nested_repos.append(rel.rstrip("/"))
            continue
        full = workspace / rel
        try:
            size = full.lstat().st_size
        except OSError:
            continue
        if size > SIZE_LIMIT_BYTES:
            plan.too_large.append(rel)
        elif not full.is_symlink() and looks_like_qm_output(full):
            plan.qm_output.append(rel)
        else:
            plan.files.append(rel)


def plan(workspace) -> Plan:
    """What :func:`init_local_project` would commit, without creating
    anything in *workspace*: a scratch repository elsewhere supplies the
    index git's ignore matching needs."""
    ws = _resolved(workspace)
    if ws is None:
        raise LocalProjectError(f"{workspace} cannot be resolved")
    text, _old = _gitignore_text(ws)
    with tempfile.TemporaryDirectory(prefix="delfin-localproj-") as tmp:
        tmp_path = Path(tmp)
        init = safe_git.run(["init", "-q", str(tmp_path / "r")], tmp_path)
        if init.returncode != 0:
            raise LocalProjectError(f"git init failed: {init.stderr.strip()}")
        excl = tmp_path / "exclude"
        excl.write_text(text, encoding="utf-8")
        names, ignored = _list_untracked(ws, tmp_path / "r" / ".git", excl)
    result = Plan(workspace=ws, ignored_count=ignored)
    _classify(ws, [n for n in names if n != ".gitignore"], result)
    result.files.append(".gitignore")
    return result


# -- init ---------------------------------------------------------------------

def _git(ws: Path, *args: str, env: Optional[dict] = None,
         timeout: float = 120) -> str:
    r = safe_git.run(list(args), ws, env=env, timeout=timeout)
    if r.returncode != 0:
        raise LocalProjectError(
            f"git {args[0]} failed: {(r.stderr or r.stdout).strip()[:400]}")
    return r.stdout


def _identity_env(ws: Path) -> Optional[dict]:
    """The fallback committer when git has no user.name / user.email."""
    have = all(
        os.environ.get(var) or safe_git.out(ws, "config", "--get", key).strip()
        for var, key in (("GIT_AUTHOR_EMAIL", "user.email"),
                         ("GIT_AUTHOR_NAME", "user.name")))
    return None if have else dict(_FALLBACK_IDENTITY)


def init_local_project(workspace) -> InitResult:
    """Make *workspace* a local repository with one first commit.

    Raises :class:`LocalProjectError` (after rolling back) on any failure;
    the caller has checked :func:`refusal` and has the user's consent."""
    ws = _resolved(workspace)
    why = refusal(ws) if ws is not None else f"{workspace} cannot be resolved"
    if why:
        raise LocalProjectError(why)
    git_dir = ws / ".git"
    if git_dir.exists():
        raise LocalProjectError(f"{git_dir} already exists")
    text, old = _gitignore_text(ws)
    gitignore = ws / ".gitignore"
    try:
        _git(ws, "init", "-q")
        _git(ws, "symbolic-ref", "HEAD", "refs/heads/main")
        gitignore.write_text(text, encoding="utf-8")
        names, ignored = _list_untracked(ws, git_dir, None)
        result = Plan(workspace=ws, ignored_count=ignored)
        _classify(ws, [n for n in names if n != ".gitignore"], result)
        if result.kept_out:
            lines = ["", "# Kept out of the first commit by DELFIN "
                     f"(above {SIZE_LIMIT_BYTES // (1024 * 1024)} MB, ORCA "
                     "output, or a nested repository)"]
            lines += [_ignore_line(rel) for rel in sorted(result.kept_out)]
            gitignore.write_text(text.rstrip("\n") + "\n" + "\n".join(lines)
                                 + "\n", encoding="utf-8")
        _git(ws, "add", "-A")
        staged = [n for n in _git(ws, "diff", "--cached", "--name-only",
                                  "-z").split("\0") if n]
        kept = set(result.kept_out)
        bad = [n for n in staged
               if n in kept or _STAGED_FORBIDDEN_RE.search(n.rsplit("/", 1)[-1])]
        if bad:
            raise LocalProjectError(
                "the first commit would hold files that must stay out (an "
                "existing .gitignore re-includes them?): "
                + ", ".join(bad[:10]))
        _git(ws, "commit", "-q", "--no-verify", "--no-gpg-sign", "-m",
             FIRST_COMMIT_MESSAGE, env=_identity_env(ws))
        sha = _git(ws, "rev-parse", "HEAD").strip()
        result.files = staged
        return InitResult(workspace=ws, commit=sha, branch="main",
                          plan=result)
    except BaseException:
        shutil.rmtree(git_dir, ignore_errors=True)
        try:
            if old is None:
                gitignore.unlink(missing_ok=True)
            else:
                gitignore.write_text(old, encoding="utf-8")
        except OSError:
            pass
        raise


def consent_question(p: Plan) -> str:
    """The text the user answers, in plain words."""
    return (
        f"Make {p.workspace} a LOCAL git repository?\n"
        "Worktrees, branches and sub-agents in their own worktrees need "
        "git, and this folder has none.\n\n"
        "A yes runs: git init (branch main), writes a .gitignore that keeps "
        "QM/ORCA output, model weights, archives, large files and secrets "
        f"out, and makes one first commit: {p.summary()}.\n"
        "No remote is added and nothing leaves this machine; adding a "
        "remote later is your step.\n"
        "A no leaves the folder as it is; you are not asked again in this "
        "session.")


__all__ = ["SIZE_LIMIT_BYTES", "MAX_FILES", "DEFAULT_GITIGNORE",
           "GITIGNORE_MARKER", "QM_OUTPUT_SIGNATURES", "LocalProjectError",
           "Plan", "InitResult", "refusal", "plan", "init_local_project",
           "consent_question", "looks_like_qm_output",
           "enclosing_repository"]
