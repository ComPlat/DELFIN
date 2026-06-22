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
import shlex
import shutil
import subprocess
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
        return True, "ok"
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

    # auto (and any unknown value) — pick best available
    if shutil.which("bwrap"):
        return SandboxConfig("bwrap", allow_net, timeout_s)
    if shutil.which("firejail"):
        return SandboxConfig("firejail", allow_net, timeout_s)
    return SandboxConfig("allowlist", allow_net, timeout_s)


# Home-relative directories that are tmpfs'd inside the sandbox so the agent
# cannot read or exfiltrate them. Everything else under $HOME stays bound
# read-only so language envs (micromamba, conda, venvs) remain usable for
# legitimate commands like ``python -m pytest`` or ``delfin --recalc``.
_HOME_SECRET_DIRS = (
    ".ssh", ".aws", ".gnupg", ".kube", ".docker", ".netrc", ".pgpass",
    ".azure", ".gcp", ".config/gcloud", ".config/gh", ".config/git",
    ".git-credentials", ".npmrc", ".pypirc",
    ".anthropic", ".openai", ".claude",
)


def _bwrap_argv(cmd: str, repo_dir: Path, allow_network: bool) -> list[str]:
    """Build a bwrap argv that runs *cmd* in a read-only namespace.

    The repo dir is bound read-write; common credential directories under
    ``$HOME`` are hidden behind tmpfs mounts; ``/root`` and ``/tmp`` are
    fresh tmpfs; ``--unshare-net`` blocks outbound traffic by default.
    """
    repo = str(repo_dir.resolve())
    home = str(Path.home())
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
        # 'off' or 'allowlist' or 'auto-fell-back-to-allowlist' — no sandbox.
        # We still avoid shell=True by routing through bash -c explicitly,
        # which produces equivalent behavior but makes argv inspectable.
        argv = ["bash", "-c", cmd]
        shell = False

    try:
        proc = subprocess.run(
            argv,
            shell=shell,
            cwd=str(repo_dir),
            capture_output=True,
            text=True,
            timeout=cfg.timeout_s,
            check=False,
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
