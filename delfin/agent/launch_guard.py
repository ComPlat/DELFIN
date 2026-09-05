"""Is this directory safe to start an agent in, and in what posture?

The terminal agent takes the directory you launched it from as its
workspace, which makes the launch directory a security decision rather
than a convenience. Everything here runs BEFORE an engine exists: pure
functions over a path, no model, no network, no terminal.

Two things this deliberately does not do.

It does not enforce the workspace floor. ``KitToolPermissions.__post_init__``
already refuses ``$HOME`` and the system roots, and that check stays exactly
where it is — it is the guarantee, and it holds for the dashboard and every
headless caller too. What this module adds is manners: the same refusal
reached one layer earlier, as a sentence that says why, instead of a
``ValueError`` traceback out of a constructor the user never called.

It does not widen anything. ``--add-dir`` grants are only VALIDATED here;
a grant this module rejects is one the permissions object would have
dropped silently, so the effect is a named refusal rather than a session
that quietly has less access than the flags asked for.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "LaunchFinding", "GitInfo", "LaunchReport",
    "inspect_launch_dir", "workspace_refusal_message", "resolve_posture",
    "REFUSE", "ASK", "NOTICE",
]

REFUSE = "refuse"
ASK = "ask"
NOTICE = "notice"

# Directories where a workspace is legitimate but the guarantees read
# differently: /tmp and friends are exempt from write gating entirely
# (api_client._is_ephemeral_sink), so a workspace there is one where the
# write gate is largely inert. Worth a sentence, not a refusal — a scratch
# workflow is a real workflow.
_EPHEMERAL_ROOTS: frozenset[str] = frozenset({
    "/tmp", "/var/tmp", "/dev/shm",
})

_VALID_MODES: frozenset[str] = frozenset({
    "plan", "default", "acceptEdits", "bypassPermissions",
})

# What an interactive session starts in when nothing says otherwise. Plan
# mode is read-only: the first turn in a folder the agent has never seen —
# exactly when it is most likely to be wrong about what the folder IS —
# cannot write, cannot shell, and cannot reach an unknown MCP tool.
DEFAULT_MODE = "plan"

_GIT_TIMEOUT_S = 5.0
_DIRTY_SHOWN = 5


@dataclass(frozen=True)
class LaunchFinding:
    """One thing worth saying about the launch directory.

    ``level`` is ``refuse`` (do not start), ``ask`` (start only if the user
    says so) or ``notice`` (say it and carry on).
    """

    level: str
    code: str
    message: str
    detail: str = ""


@dataclass(frozen=True)
class GitInfo:
    is_repo: bool = False
    toplevel: Path | None = None
    branch: str = ""
    dirty: tuple[str, ...] = ()
    #: A repository with no commit yet. `git init` then `delfin-agent` is
    #: an ordinary way to start, and on that repository `rev-parse HEAD`
    #: fails — which is indistinguishable from a detached HEAD unless it
    #: is asked as its own question.
    unborn: bool = False


@dataclass(frozen=True)
class LaunchReport:
    workspace: Path
    findings: tuple[LaunchFinding, ...] = ()
    git: GitInfo = field(default_factory=GitInfo)
    trust_notices: tuple[str, ...] = ()
    granted_dirs: tuple[Path, ...] = ()
    read_dirs: tuple[Path, ...] = ()

    @property
    def refused(self) -> bool:
        return any(f.level == REFUSE for f in self.findings)

    @property
    def questions(self) -> tuple[LaunchFinding, ...]:
        return tuple(f for f in self.findings if f.level == ASK)

    def render(self) -> str:
        """The findings as text, refusals first. Never a traceback."""
        order = {REFUSE: 0, ASK: 1, NOTICE: 2}
        lines: list[str] = []
        for finding in sorted(self.findings, key=lambda f: order.get(f.level, 9)):
            lines.append(finding.message)
            if finding.detail:
                lines.append(finding.detail)
            lines.append("")
        return "\n".join(lines).rstrip()


def workspace_refusal_message(path: Path | str) -> str:
    """Why the agent will not take *path* as a workspace, in prose."""
    return (
        f"delfin-agent will not run in {path}.\n\n"
        "That is your home directory or a system root. The agent's workspace "
        "is the folder you launch it in, and everything inside it is "
        "writable — here that would mean your keys, your browser profile and "
        "every other project you have.\n\n"
        "    cd ~/some-project && delfin-agent"
    )


def _git(workspace: Path, *args: str) -> str:
    try:
        proc = subprocess.run(
            ["git", *args], cwd=str(workspace), capture_output=True,
            text=True, timeout=_GIT_TIMEOUT_S, check=False,
        )
    except Exception:
        return ""
    # Only the trailing newline. `--porcelain` puts the status in columns
    # 1-2, so stripping the whole output eats the leading space of the
    # FIRST line and shifts its path by one character — which is exactly
    # the kind of defect that shows up as one wrong filename in a banner
    # and nowhere else.
    return proc.stdout.rstrip("\n") if proc.returncode == 0 else ""


def _porcelain_path(line: str) -> str:
    """The path out of one `git status --porcelain` line.

    Columns 1-2 are the index and worktree status, column 3 a space, and
    the path starts at column 4. A rename reads ``R  old -> new``; the new
    name is the one that exists.
    """
    path = line[3:]
    if " -> " in path:
        path = path.split(" -> ", 1)[1]
    return path.strip().strip('"')


#: The agent's own state directory. It is created by the agent, in the
#: workspace, usually within the first second of the first session — so
#: reporting it as work the user has not committed says "you have loose
#: changes here" about a directory the user has never touched.
_OWN_STATE_DIR = ".delfin/"


def _git_info(workspace: Path) -> GitInfo:
    top = _git(workspace, "rev-parse", "--show-toplevel").strip()
    if not top:
        return GitInfo(is_repo=False)
    # `symbolic-ref` answers on a repository with no commit yet, where
    # `rev-parse --abbrev-ref HEAD` exits 128 and yields nothing. Asking
    # only the second one made a freshly initialised repository report
    # itself as a detached HEAD, which is a different and alarming state.
    branch = _git(workspace, "symbolic-ref", "--short", "HEAD").strip()
    if not branch:
        branch = _git(workspace, "rev-parse", "--abbrev-ref", "HEAD").strip()
        if branch == "HEAD":
            branch = ""              # genuinely detached
    unborn = not _git(workspace, "rev-parse", "--verify", "HEAD").strip()
    status = _git(workspace, "status", "--porcelain")
    dirty = tuple(
        path for path in (
            _porcelain_path(line) for line in status.splitlines()
            if len(line) > 3
        )
        if path and not path.startswith(_OWN_STATE_DIR)
    )
    try:
        toplevel: Path | None = Path(top).resolve()
    except Exception:
        toplevel = None
    return GitInfo(is_repo=True, toplevel=toplevel, branch=branch,
                   dirty=dirty, unborn=unborn)


def _under_ephemeral_root(path: Path) -> str:
    for root in _EPHEMERAL_ROOTS:
        rp = Path(root)
        if path == rp or rp in path.parents:
            return root
    return ""


def _check_grant(raw: str, workspace: Path, *, writable: bool
                 ) -> tuple[Path | None, LaunchFinding | None]:
    """Validate one --add-dir / --read-dir before the engine sees it.

    KitToolPermissions drops a forbidden extra directory silently, so a
    typo'd grant becomes a session that quietly has less reach than the
    command line asked for. Refusing by name is the whole point.
    """
    from .api_client import _is_forbidden_workspace_root, _SECRET_DIR_NAMES

    flag = "--add-dir" if writable else "--read-dir"
    try:
        path = Path(raw).expanduser().resolve()
    except Exception:
        return None, LaunchFinding(
            REFUSE, "bad_grant", f"{flag} {raw}: not a usable path")
    if not path.is_dir():
        return None, LaunchFinding(
            REFUSE, "bad_grant", f"{flag} {path}: not a directory")
    if set(path.parts) & _SECRET_DIR_NAMES:
        return None, LaunchFinding(
            REFUSE, "secret_grant",
            f"{flag} {path}: that is a credential directory.")
    if writable:
        if _is_forbidden_workspace_root(path):
            return None, LaunchFinding(
                REFUSE, "bad_grant",
                f"{flag} {path}: too broad to hand over — it is your home "
                "directory or a system root.")
        if path == workspace or path in workspace.parents:
            return None, LaunchFinding(
                REFUSE, "bad_grant",
                f"{flag} {path}: that contains the workspace, so it would "
                "widen the session rather than extend it.")
    return path, None


def inspect_launch_dir(
    cwd: Path | str,
    *,
    add_dirs: tuple[str, ...] = (),
    read_dirs: tuple[str, ...] = (),
    check_trust: bool = True,
) -> LaunchReport:
    """Everything worth knowing about *cwd* before an engine is built."""
    from .api_client import _is_forbidden_workspace_root, _is_delfin_workspace

    try:
        workspace = Path(cwd).expanduser().resolve()
    except Exception:
        workspace = Path(str(cwd))

    findings: list[LaunchFinding] = []

    # 1. The floor. Refusing here is manners; the constructor refuses too.
    if _is_forbidden_workspace_root(workspace):
        return LaunchReport(
            workspace=workspace,
            findings=(LaunchFinding(
                REFUSE, "forbidden_root", workspace_refusal_message(workspace)),),
        )

    if not workspace.is_dir():
        return LaunchReport(
            workspace=workspace,
            findings=(LaunchFinding(
                REFUSE, "missing_dir",
                f"delfin-agent will not run in {workspace}: no such "
                "directory."),),
        )

    # 2. Grants, validated before the permissions object silently drops them.
    granted: list[Path] = []
    readable: list[Path] = []
    for raw in add_dirs:
        path, bad = _check_grant(raw, workspace, writable=True)
        (findings.append(bad) if bad else granted.append(path))  # type: ignore[arg-type]
    for raw in read_dirs:
        path, bad = _check_grant(raw, workspace, writable=False)
        (findings.append(bad) if bad else readable.append(path))  # type: ignore[arg-type]

    # 3. Ephemeral roots — legitimate, but the write gate is mostly inert.
    root = _under_ephemeral_root(workspace)
    if root:
        findings.append(LaunchFinding(
            NOTICE, "ephemeral_root",
            f"{workspace} is under {root}. Writes there are exempt from the "
            "path gate, so approving a write is not the check it is "
            "elsewhere."))

    # 4. Git — what can be undone, and what is already changed.
    git = _git_info(workspace)
    if not git.is_repo:
        findings.append(LaunchFinding(
            NOTICE, "not_a_git_repo",
            "Not a git repository — nothing the agent changes here can be "
            "recovered with git. `git init` first, or use /rewind."))
    else:
        if git.toplevel is not None and git.toplevel != workspace:
            findings.append(LaunchFinding(
                NOTICE, "subdir_of_repo",
                f"This is a sub-directory of {git.toplevel}. The agent sees "
                "only this folder, but its changes land in that repository."))
        if git.dirty:
            shown = ", ".join(git.dirty[:_DIRTY_SHOWN])
            more = (f" (+{len(git.dirty) - _DIRTY_SHOWN} more)"
                    if len(git.dirty) > _DIRTY_SHOWN else "")
            findings.append(LaunchFinding(
                NOTICE, "dirty_tree",
                f"{len(git.dirty)} file(s) already modified: {shown}{more}",
                detail="If the agent changes more, its work and yours will "
                       "be hard to tell apart."))

    # 5. A DELFIN checkout gets a wider shell allow-list, with nothing today
    #    telling anyone that it does.
    try:
        if _is_delfin_workspace(workspace):
            findings.append(LaunchFinding(
                NOTICE, "delfin_workspace",
                "Recognised as a DELFIN checkout: the shell auto-allow list "
                "is wider here than in an ordinary project."))
    except Exception:
        pass

    # 6. Trust. These notices have existed and had no caller — a repository
    #    that ships hooks has them withheld, and nobody was ever told.
    notices: tuple[str, ...] = ()
    if check_trust:
        try:
            from . import workspace_trust
            notices = tuple(workspace_trust.pending_notices(workspace))
        except Exception:
            notices = ()
        if notices:
            # NOTICE, not ASK, and the distinction is the whole point:
            # withholding IS the safe state. Nothing here is loaded, so
            # there is no decision to block the session on — only a fact
            # the user could not otherwise learn, and the way to change it
            # if they want to. Blocking would train people to say yes.
            findings.append(LaunchFinding(
                NOTICE, "untrusted_workspace",
                "This folder offers definitions that are being withheld "
                "because it is not trusted:",
                detail="\n".join(f"  {n}" for n in notices)
                       + "\n  Nothing from them runs. `/trust` loads them."))

    return LaunchReport(
        workspace=workspace,
        findings=tuple(findings),
        git=git,
        trust_notices=notices,
        granted_dirs=tuple(granted),
        read_dirs=tuple(readable),
    )


def resolve_posture(
    *,
    flag_mode: str = "",
    persisted_mode: str = "",
    unattended_opt_in: bool = False,
    settings_path: str = "~/.delfin/settings.json",
) -> tuple[str, str]:
    """Decide the permission mode an interactive session starts in.

    Returns ``(mode, why)`` where *why* is shown in the banner, so the user
    can always see which of the three inputs decided it.

    The one rule that is not a preference: a persisted
    ``bypassPermissions`` is NOT inherited. That setting is written once —
    often in a container, often to get through one stubborn task — and it
    then makes every future session on the machine unattended, silently.
    An interactive session downgrades to plan mode and says where the
    setting lives. Reaching bypass from a terminal takes two flags typed
    together, on purpose.
    """
    flag_mode = (flag_mode or "").strip()
    persisted_mode = (persisted_mode or "").strip()

    if flag_mode:
        if flag_mode == "bypassPermissions" and not unattended_opt_in:
            return DEFAULT_MODE, (
                "refused --permission-mode bypassPermissions without "
                "--unattended; starting in plan")
        return flag_mode, "--permission-mode"

    if persisted_mode == "bypassPermissions":
        return DEFAULT_MODE, (
            f"{settings_path} sets default_mode = \"bypassPermissions\". An "
            "interactive session will not start unattended; starting in "
            "plan. Use --permission-mode acceptEdits for this session, or "
            "--permission-mode bypassPermissions --unattended if you mean it")

    if persisted_mode in _VALID_MODES:
        return persisted_mode, settings_path

    # Not the word "default": next to "approval plan" it reads as a claim
    # about the approval mode rather than about where the mode came from,
    # which is the one thing this string exists to say.
    return DEFAULT_MODE, "nothing configured it"
