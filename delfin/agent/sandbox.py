"""Defense-in-depth runner for the dashboard agent's bash execution.

Layers (applied for every command before it is executed):

  1. Allow-list   — fail-closed for unknown first-tokens; reject explicit
                    deny-substrings even if the first-token is allowed.
  2. Sandbox      — wrap in ``bwrap`` (preferred) or ``firejail``: read-only
                    rootfs, repo-dir bind-mounted rw, tmpfs over /home and
                    /tmp, network unshared by default.
  3. Audit log    — append every command + exit + bytes to a 0600 JSONL file
                    under ``$XDG_CACHE_HOME/delfin/agent-audit.jsonl``.
  4. Approval     — caller's responsibility (UI shows full command, no
                    truncation); see ``dashboard/tab_agent.py``.

Env vars (read on every call so tests can flip them with monkeypatch):

  ``DELFIN_AGENT_SANDBOX``           one of {auto,off,allowlist,bwrap,firejail}.
                                     Default: auto. ``auto`` resolves to the
                                     first available of bwrap → firejail →
                                     allowlist-only. ``off`` disables BOTH
                                     allow-list and sandbox (dangerous; debug
                                     only).
  ``DELFIN_AGENT_SANDBOX_NETWORK``   set to ``1`` to keep network in sandbox
                                     (e.g. ``pip install``). Default: deny.
  ``DELFIN_AGENT_SANDBOX_TIMEOUT``   per-command timeout in seconds. Default 60.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Layer 1 — allow-list

# First token of each pipeline segment must be in this set OR (for the few
# commands with sub-command grammar) match an entry in ``_ALLOW_SUBCOMMAND``.
_ALLOW_FIRST_TOKEN = frozenset({
    # Read-only inspection
    "ls", "cat", "head", "tail", "less", "more", "wc", "grep", "egrep",
    "fgrep", "rg", "find", "stat", "file", "tree", "du", "df", "echo",
    "pwd", "env", "which", "whoami", "uname", "date", "id", "hostname",
    "basename", "dirname", "realpath", "readlink", "sort", "uniq", "tr",
    "cut", "awk", "sed", "diff", "cmp", "tee", "xargs",
    # Python tooling
    "python", "python3", "pytest", "pip", "ruff", "black", "mypy",
    "flake8", "isort",
    # DELFIN entry-points
    "delfin", "delfin-build", "delfin-voila", "delfin-json",
    "delfin_ESD", "delfin_IR",
    # Chemistry / QM tools commonly invoked from DELFIN
    "xtb", "orca", "crest", "obabel", "x2t", "t2x", "molden", "multiwfn",
    "censo", "anmr",
    # Process / shell built-ins that are safe
    "ps", "true", "false", "test", "[", "yes", "sleep", "wait",
    "nohup",
})

# Per-binary subcommand allow-list. Used when the binary itself is
# multi-purpose (e.g. ``git``, ``pip``) — we allow only read-only or harmless
# subcommands here. Other subcommands fall through to deny.
_ALLOW_SUBCOMMAND = {
    "git": frozenset({
        "status", "log", "diff", "show", "branch", "remote", "rev-parse",
        "rev-list", "stash", "ls-files", "ls-tree", "blame", "describe",
        "tag", "shortlog", "fsck", "cat-file", "for-each-ref", "name-rev",
        "config",  # config alone is allowed; --global is deny-substring'd below
    }),
    "pip": frozenset({"list", "show", "freeze", "check", "config"}),
}

# Always-deny substrings (case-insensitive). These win even if the
# first-token check passes — designed to catch obvious destructive or
# escalation patterns and unsafe redirects.
_DENY_SUBSTRINGS = (
    " --no-verify",
    " --force ", " --force\t",
    " push --force", " push -f ",
    " reset --hard",
    " clean -fd", " clean -ffd",
    "rm -rf /", "rm -rf ~", "rm -rf .",
    "rm -fr /", "rm -fr ~", "rm -fr .",
    "mkfs", "dd if=/", "dd of=/dev/sd",
    "> /dev/sd", "> /dev/nvme",
    ":(){ :|:& };:",
    "sudo ", "sudo\t", " su -", " su\t-",
    "chmod 777", "chmod -r 777", "chmod -r 7",
    "git config --global",
    "ssh-keygen", "ssh-add",
    "curl http",  # discourage remote-fetch + pipe-to-shell
    "wget http",
    " | sh", " | bash", " |sh", " |bash",
    "/etc/passwd", "/etc/shadow",
    "~/.ssh", "~/.aws", "~/.gnupg",
    "/.ssh/id_", "/.ssh/known_hosts",
    "/.aws/credentials",
)


# Commands that recurse through a directory tree. Harmless on a job folder,
# but rooted at a whole file system they cost one metadata RPC per file --
# hundreds of thousands of them on a populated HOME -- which degrades a shared
# machine for every user on it. A site's operations team will read that as
# abuse and block the account, which is what prompted this guard.
#
# The rule derives its roots at runtime from the current user's HOME and the
# machine's mount table, so it holds for any user on any cluster without
# knowing site-specific paths.
_TREE_WALK_COMMANDS = frozenset({"du", "ncdu", "find", "tree", "rg"})
#: Walk only when told to recurse: grep -r/-R, ls -R.
_RECURSIVE_ON_FLAG = frozenset({"grep", "egrep", "fgrep", "ls"})


def _is_tree_walk_root(resolved: Path) -> bool:
    """True if a recursive walk rooted here would sweep a whole file system.

    Two universal signals, no hardcoded site paths:
      * the path is the user's HOME or an ancestor of it, and
      * the path is a mount point (``/scratch``, ``/work``, project file
        systems, ``/`` itself) -- a job directory never is.
    """
    if str(resolved) == os.sep:
        return True
    try:
        home = Path.home().resolve()
        if resolved == home or resolved in home.parents:
            return True
    except Exception:
        pass
    try:
        return os.path.ismount(str(resolved))
    except OSError:
        return False


def _recurses(base: str, toks: list[str]) -> bool:
    """Whether this command walks a tree at all."""
    if base in _TREE_WALK_COMMANDS:
        return True
    if base not in _RECURSIVE_ON_FLAG:
        return False
    letter = "R" if base == "ls" else "rR"
    for t in toks[1:]:
        if t in ("--recursive", "--dereference-recursive"):
            return True
        if t.startswith("-") and not t.startswith("--") \
                and any(c in t[1:] for c in letter):
            return True
    return False


def _check_tree_walk(toks: list[str], cwd=None) -> tuple[bool, str]:
    """Deny recursive walks rooted at a whole file system; subdirs are fine.

    The one implementation of this rule. The approval runner of the CLI
    backend reaches it through is_allowed, and the gate of every API and
    terminal session through tree_walk_refusal: it used to live on the
    CLI path only, so the sessions that run most of the work walked past
    it (measured 2026-09-22: `find ~ -name '*.py'`, three times, and walks
    over the archive, all executed).
    """
    base = os.path.basename(toks[0])
    if not _recurses(base, toks):
        return True, "ok"

    targets = [t for t in toks[1:] if not t.startswith("-")]
    if not targets or (base == "rg" and len(targets) == 1) \
            or (base in ("grep", "egrep", "fgrep") and len(targets) == 1):
        # No path given (rg/grep: only the pattern): the walk is rooted
        # at the working directory.
        targets.append(".")

    for target in targets:
        try:
            path = Path(os.path.expanduser(target))
            if not path.is_absolute() and cwd is not None:
                path = Path(cwd) / path       # where the COMMAND runs
            resolved = path.resolve()
        except Exception:
            continue
        if _is_tree_walk_root(resolved):
            return False, (
                f"{base} would walk all of {str(resolved)!r}; on a shared file "
                f"system that overloads the metadata servers. For disk usage "
                f"use delfin.quota.home_usage(); otherwise point {base} at a "
                f"specific subdirectory."
            )
    return True, "ok"


def tree_walk_refusal(cmd: str, cwd=None) -> Optional[str]:
    """Why *cmd* may not run as a walk over a whole file system, or None.

    Checks every pipeline segment, whatever the command's other merits --
    a walk rooted at the home, one of its ancestors or a mount point is
    refused in every mode, before anything is asked. Relative paths are
    read from *cwd*, the directory the command runs in (not this
    process's), and a `cd` inside the command moves it: `cd ~ && find .`
    is a walk over the home.
    """
    try:
        segments = _split_pipeline(cmd.strip())
    except ValueError:
        segments = [s for s in re.split(r"\|\||&&|[|;&\n]", cmd) if s.strip()]
    for seg in segments:
        try:
            toks = shlex.split(seg, posix=True)
        except ValueError:
            toks = seg.split()
        while toks and "=" in toks[0] and toks[0].split("=", 1)[0].isidentifier():
            toks = toks[1:]
        # Wrappers that run the next word: the walk is that word's.
        while toks and os.path.basename(toks[0]) in ("nice", "ionice", "time",
                                                      "timeout", "command",
                                                      "exec"):
            toks = toks[1:]
            # their own options and values: -n 10, -c2, a duration like 5s
            while toks and (toks[0].startswith("-")
                            or re.fullmatch(r"\d+[smhd]?", toks[0])):
                toks = toks[1:]
        if not toks:
            continue
        if toks[0] == "cd":
            dest = toks[1] if len(toks) > 1 else "~"
            try:
                p = Path(os.path.expanduser(dest))
                if not p.is_absolute() and cwd is not None:
                    p = Path(cwd) / p
                cwd = str(p.resolve())
            except Exception:
                pass
            continue
        ok, reason = _check_tree_walk(toks, cwd)
        if not ok:
            return reason
    return None


@dataclass
class AllowResult:
    allowed: bool
    reason: str


def is_allowed(cmd: str) -> AllowResult:
    """Layer 1 — allow-list check on a shell command string.

    Splits the command into pipeline segments at ``|``, ``||``, ``&&``,
    ``;``, ``&`` and verifies each segment's first token. Also rejects on
    any always-deny substring.
    """
    cmd_stripped = cmd.strip()
    if not cmd_stripped:
        return AllowResult(False, "empty command")

    # Always-deny substring scan (case-insensitive, with surrounding spaces
    # so word-boundary patterns like " --force " don't false-positive on a
    # filename like ``--force-tag.txt``).
    haystack = " " + cmd_stripped.lower() + " "
    for needle in _DENY_SUBSTRINGS:
        if needle in haystack:
            return AllowResult(False, f"matched deny pattern: {needle.strip()!r}")

    # Pipeline-aware first-token check
    try:
        segments = _split_pipeline(cmd_stripped)
    except ValueError as exc:
        return AllowResult(False, f"unparseable command: {exc}")
    for seg in segments:
        ok, reason = _check_segment(seg)
        if not ok:
            return AllowResult(False, reason)
    return AllowResult(True, "ok")


def _split_pipeline(cmd: str) -> list[str]:
    """Split *cmd* at unquoted ``|``, ``||``, ``&&``, ``;``, ``&``."""
    toks = shlex.split(cmd, posix=True)
    out: list[list[str]] = [[]]
    for t in toks:
        if t in ("|", "||", "&&", ";", "&"):
            out.append([])
        else:
            out[-1].append(t)
    return [" ".join(seg) for seg in out if seg]


def _check_segment(seg: str) -> tuple[bool, str]:
    toks = shlex.split(seg, posix=True)
    if not toks:
        return False, "empty pipeline segment"
    head = toks[0]
    # Strip leading env-var assignments (e.g. ``FOO=1 python …``)
    while head and "=" in head and head.split("=", 1)[0].isidentifier():
        toks = toks[1:]
        if not toks:
            return False, "command was just env-vars"
        head = toks[0]
    base = os.path.basename(head)
    if base in _ALLOW_FIRST_TOKEN:
        return _check_tree_walk(toks)
    if base in _ALLOW_SUBCOMMAND:
        if len(toks) >= 2 and toks[1] in _ALLOW_SUBCOMMAND[base]:
            return True, "ok"
        return False, f"{base} subcommand not in allow-list"
    return False, f"first token not in allow-list: {base!r}"


# ---------------------------------------------------------------------------
# Layer 2 — sandbox

@dataclass
class SandboxConfig:
    mode: str           # "off" | "allowlist" | "bwrap" | "firejail"
    allow_network: bool
    timeout_s: int


def detect_config(env: Optional[dict] = None) -> SandboxConfig:
    """Read env vars and decide which sandbox mode to use."""
    e = env if env is not None else os.environ
    requested = (e.get("DELFIN_AGENT_SANDBOX") or "auto").strip().lower()
    allow_net = e.get("DELFIN_AGENT_SANDBOX_NETWORK", "0").strip() == "1"
    try:
        timeout_s = int(e.get("DELFIN_AGENT_SANDBOX_TIMEOUT", "60"))
    except ValueError:
        timeout_s = 60

    if requested == "off":
        return SandboxConfig("off", allow_net, timeout_s)
    if requested == "allowlist":
        return SandboxConfig("allowlist", allow_net, timeout_s)
    if requested == "bwrap":
        return SandboxConfig("bwrap", allow_net, timeout_s)
    if requested == "firejail":
        return SandboxConfig("firejail", allow_net, timeout_s)

    # auto (and any unknown value) — pick best available. Installed is not
    # working: a bwrap whose user namespace is refused made every approved
    # command fail.
    if shutil.which("bwrap") and _bwrap_works():
        return SandboxConfig("bwrap", allow_net, timeout_s)
    if shutil.which("firejail"):
        return SandboxConfig("firejail", allow_net, timeout_s)
    return SandboxConfig("allowlist", allow_net, timeout_s)


def _bwrap_works() -> bool:
    try:
        from .api_client import _bwrap_functional
        return bool(_bwrap_functional())
    except Exception:
        return False


def _unsandboxed_argv(cmd: str, repo_dir: Path, mode: str) -> list[str]:
    """``bash -c cmd`` for the modes without bubblewrap or firejail.

    The allowlist mode ran commands with nothing around them: python, pip,
    awk and env are on the list, and every key was in the environment. Where
    the kernel has Landlock the command now runs under it (writes only in the
    repository and a private temp directory, credential folders unreadable),
    as the agent's own shell does on such a host. "off" stays off."""
    base = ["bash", "-c", cmd]
    if mode != "allowlist":
        return base
    try:
        from .api_client import (_home_secret_paths, _landlock_functional,
                                 _private_tmp_dir)
        from . import landlock_exec as _ll
        if not _landlock_functional():
            return base
        argv = [sys.executable, "-I", str(Path(_ll.__file__).resolve()),
                "--write", str(Path(repo_dir).resolve())]
        for h in _home_secret_paths():
            argv += ["--hide", h]
        return argv + ["--tmpdir", _private_tmp_dir(), "--strict", "0", "--"] + base
    except Exception:
        return base


# Home-relative directories that are tmpfs'd inside the sandbox so the agent
# cannot read or exfiltrate them. Everything else under $HOME stays bound
# read-only so language envs (micromamba, conda, venvs) remain usable for
# legitimate commands like ``python -m pytest`` or ``delfin --recalc``.
_HOME_SECRET_DIRS = (
    ".ssh", ".aws", ".gnupg", ".kube", ".docker", ".netrc", ".pgpass",
    ".azure", ".gcp", ".config/gcloud", ".config/gh",
    # NOT .config/git: it is git's own config and ignore file, no
    # credential, and hiding it made git fatal ("cannot use
    # ~/.config/git/ignore as an exclude file"). Credentials live in
    # .git-credentials and a helper, both denied separately.
    ".git-credentials", ".npmrc", ".pypirc",
    ".anthropic", ".openai", ".claude",
    # The framework's own credential store. The other providers' folders
    # were masked and this one was not, which is the wrong way round: it
    # is the file holding the key the running agent is using.
    ".delfin/credentials.json",
)


def _bwrap_argv(cmd: str, repo_dir: Path, allow_network: bool) -> list[str]:
    """Build a bwrap argv that runs *cmd* in a read-only namespace.

    The repo dir is bound read-write; common credential directories under
    ``$HOME`` are hidden behind tmpfs mounts; ``/root`` and ``/tmp`` are
    fresh tmpfs; ``--unshare-net`` blocks outbound traffic by default.
    """
    repo = str(repo_dir.resolve())
    # Resolve $HOME fully: on HPC clusters /home is often a symlink
    # (e.g. /home -> /pfs) and bwrap cannot create mountpoints through
    # symlinked path components ("Can't mkdir parents").
    home = str(Path.home().resolve())
    args: list[str] = [
        "bwrap",
        "--ro-bind", "/", "/",
        "--dev", "/dev",
        "--proc", "/proc",
        "--tmpfs", "/tmp",
        "--tmpfs", "/root",
    ]
    # Mask each secret path individually if it exists on the host. tmpfs
    # only works on directories; for single files we mask with /dev/null
    # (read-only bind of an empty source). Missing paths are skipped.
    for rel in _HOME_SECRET_DIRS:
        p = Path(home) / rel
        try:
            is_dir = p.is_dir()
            is_file = p.is_file()
        except OSError:
            continue
        if is_dir:
            args += ["--tmpfs", str(p)]
        elif is_file:
            args += ["--ro-bind", "/dev/null", str(p)]
    args += [
        "--bind", repo, repo,
        "--chdir", repo,
        "--unshare-pid",
        "--unshare-ipc",
        "--unshare-uts",
        "--die-with-parent",
        "--new-session",
    ]
    if not allow_network:
        args.append("--unshare-net")
    args += ["--", "bash", "-c", cmd]
    return args


def _firejail_argv(cmd: str, repo_dir: Path, allow_network: bool) -> list[str]:
    """Build a firejail argv with comparable restrictions."""
    repo = str(repo_dir.resolve())
    args: list[str] = [
        "firejail",
        "--quiet",
        "--noprofile",
        "--private-tmp",
        "--private-dev",
        f"--whitelist={repo}",
        "--read-only=/etc",
        "--read-only=/usr",
        "--noroot",
        "--caps.drop=all",
    ]
    if not allow_network:
        args.append("--net=none")
    args += [f"--chdir={repo}", "bash", "-c", cmd]
    return args


# ---------------------------------------------------------------------------
# Layer 3 — audit log

def _audit_log_path() -> Path:
    cache_root = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(cache_root) / "delfin" / "agent-audit.jsonl"


def _audit(record: dict) -> None:
    """Append one JSON record to the per-user audit log (mode 0600)."""
    path = _audit_log_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        return  # audit must never break the agent
    line = json.dumps(record, ensure_ascii=False, default=str) + "\n"
    try:
        # O_APPEND is atomic for small writes; create with 0600.
        fd = os.open(
            str(path),
            os.O_WRONLY | os.O_CREAT | os.O_APPEND,
            0o600,
        )
        try:
            os.write(fd, line.encode("utf-8"))
        finally:
            os.close(fd)
    except OSError:
        return


# ---------------------------------------------------------------------------
# Public entry point

@dataclass
class RunResult:
    blocked: bool
    block_reason: Optional[str]
    returncode: int
    stdout: str
    stderr: str
    mode: str               # the sandbox mode actually used
    timed_out: bool = False
    elapsed_s: float = 0.0
    extra: dict = field(default_factory=dict)


def run_agent_command(
    cmd: str,
    repo_dir: Path,
    *,
    config: Optional[SandboxConfig] = None,
) -> RunResult:
    """Run *cmd* through the layered defense and return its result.

    The caller is responsible for Layer 4 (interactive approval) — this
    function trusts that the user has already approved *cmd* in the UI.
    """
    cfg = config or detect_config()
    t0 = time.monotonic()

    # Layer 1 — allow-list (skipped only in 'off' mode)
    if cfg.mode != "off":
        decision = is_allowed(cmd)
        if not decision.allowed:
            elapsed = time.monotonic() - t0
            _audit({
                "ts": time.time(),
                "cwd": str(repo_dir),
                "mode": cfg.mode,
                "cmd": cmd,
                "blocked": True,
                "block_reason": decision.reason,
                "elapsed_s": round(elapsed, 4),
            })
            return RunResult(
                blocked=True,
                block_reason=decision.reason,
                returncode=-1,
                stdout="",
                stderr="",
                mode=cfg.mode,
                elapsed_s=elapsed,
            )

    # Layer 2 — choose argv based on sandbox mode
    if cfg.mode == "bwrap":
        argv = _bwrap_argv(cmd, repo_dir, cfg.allow_network)
        shell = False
    elif cfg.mode == "firejail":
        argv = _firejail_argv(cmd, repo_dir, cfg.allow_network)
        shell = False
    else:
        # 'off' or 'allowlist' or 'auto-fell-back-to-allowlist'. Routed
        # through bash -c explicitly so the argv stays inspectable; the
        # allowlist mode runs under Landlock where the kernel has it.
        argv = _unsandboxed_argv(cmd, repo_dir, cfg.mode)
        shell = False

    try:
        # Every mode: the shell's scrubbed environment (no provider keys),
        # its own session and process group, no terminal, the whole group
        # ended at exit and at the timeout.
        from .api_client import _scrubbed_bash_env
        from . import contained_run as _contained
        proc = _contained.run(
            argv,
            cwd=str(repo_dir),
            env=_scrubbed_bash_env(),
            timeout=cfg.timeout_s,
        )
        elapsed = time.monotonic() - t0
        _audit({
            "ts": time.time(),
            "cwd": str(repo_dir),
            "mode": cfg.mode,
            "cmd": cmd,
            "blocked": False,
            "exit": proc.returncode,
            "stdout_bytes": len(proc.stdout),
            "stderr_bytes": len(proc.stderr),
            "elapsed_s": round(elapsed, 4),
        })
        return RunResult(
            blocked=False,
            block_reason=None,
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            mode=cfg.mode,
            elapsed_s=elapsed,
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = time.monotonic() - t0
        _audit({
            "ts": time.time(),
            "cwd": str(repo_dir),
            "mode": cfg.mode,
            "cmd": cmd,
            "blocked": False,
            "exit": -1,
            "timed_out": True,
            "elapsed_s": round(elapsed, 4),
        })
        return RunResult(
            blocked=False,
            block_reason=None,
            returncode=-1,
            stdout=(exc.stdout or b"").decode("utf-8", errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or ""),
            stderr=(exc.stderr or b"").decode("utf-8", errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or ""),
            mode=cfg.mode,
            timed_out=True,
            elapsed_s=elapsed,
        )
