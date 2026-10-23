"""The individual checks. Every one takes the home directory as an
argument and reads only fixed, named paths beneath it -- no tree walks,
no network, no subprocess beyond the read-only commands `last`, `ss`
and `crontab` (and those only through run_read_only_command, which
tests can stub).
"""
from __future__ import annotations

import os
import re
import stat
import subprocess
from pathlib import Path

from .model import Finding

# Fixed paths a check may look at. Never a glob over HOME.
SSH_DIR_FILES = ("authorized_keys", "rc", "environment", "config",
                 "known_hosts")
RC_FILES = (".bashrc", ".bash_profile", ".profile", ".zshrc")
BIN_DIRS = (Path(".local") / "bin", Path("bin"))
CRED_FILES = (Path(".delfin") / "credentials.json", Path(".kit_env"),
              Path(".netrc"), Path(".git-credentials"))
HISTORY_FILES = (".bash_history", ".zsh_history")

# Download-and-execute and persistence patterns in start-up files.
PERSISTENCE_PATTERNS = (
    ("curl|sh", re.compile(r"curl[^|]*\|\s*(ba)?sh\b")),
    ("wget-pipe", re.compile(r"wget[^|]*-O-\s*\|")),
    ("dev-tcp", re.compile(r"/dev/tcp/")),
    ("nc-exec", re.compile(r"\b(nc|ncat|socat)\b[^#\n]*\s(-e|--exec)\b")),
    ("base64-decode-pipe", re.compile(r"base64\s+-d\s*\|")),
    ("ld-preload", re.compile(r"\bLD_PRELOAD\b")),
    ("prompt-command", re.compile(r"\bPROMPT_COMMAND\b")),
    ("alias-sudo", re.compile(r"alias\s+sudo\s*=")),
    ("alias-ssh", re.compile(r"alias\s+ssh\s*=")),
)

# Secret-shaped patterns in history files -- a hit is counted, the
# matched text is never reported.
SECRET_PATTERNS = (
    ("aws-access-key-id", re.compile(r"AKIA[0-9A-Z]{16}")),
    ("aws-secret-access-key", re.compile(
        r"AWS_SECRET_ACCESS_KEY\s*=\s*\S+")),
    ("private-key-block", re.compile(r"BEGIN (RSA|OPENSSH|EC) PRIVATE KEY")),
    ("api-token-assignment", re.compile(
        r"(API_KEY|API_TOKEN|GITHUB_TOKEN)\s*=\s*\S+")),
)

SSH_CONFIG_RISKY = (
    ("ProxyCommand", re.compile(r"^\s*ProxyCommand\s+\S+", re.M)),
    ("LocalCommand", re.compile(r"^\s*LocalCommand\s+\S+", re.M)),
    ("ForwardAgent", re.compile(
        r"^\s*ForwardAgent\s+(yes|ask)\s*$", re.M)),
)

SUSPICIOUS_PORTS = (4444, 1337, 6667)

# Reverse-shell markers, matched on the argv entries joined with
# spaces (so flags that span several argv entries, like `nc -e`, are
# still seen). `bash -i` alone is deliberately NOT a marker: without
# the socket redirection it is any interactive shell; the redirect
# carries /dev/tcp/ and trips dev-tcp.
# Reverse-shell shapes, judged on the PROGRAM and its own arguments. A
# text search over the whole command line flagged any process that merely
# carried these words -- the first live run alarmed on a DELFIN agent
# whose task text described them (night run 2026-09-25).
_SHELLS = frozenset({"sh", "bash", "zsh", "dash", "ksh", "ash"})
_NETCATS = frozenset({"nc", "ncat", "netcat"})
# A shell that dials out names host AND port: /dev/tcp/<host>/<port>. The
# bare words are text -- a shell whose script merely mentions them (a
# commit message, a grep) is no connection.
_DEV_NET_TARGET = re.compile(r"/dev/(?:tcp|udp)/[A-Za-z0-9_.:-]+/\d+")


def _reverse_shell_shapes(argv: list[str]) -> list[str]:
    """Labels of the reverse-shell shapes this argv has, from its program."""
    if not argv or not argv[0]:
        return []
    prog = os.path.basename(argv[0])
    args = argv[1:]
    found = []
    if prog in _SHELLS and any(_DEV_NET_TARGET.search(a) for a in args):
        found.append("dev-tcp")
    if prog in _NETCATS and any(a in ("-e", "--exec", "-c", "--sh-exec")
                                or (a.startswith("-e") and len(a) > 2)
                                for a in args):
        found.append("nc-exec")
    if prog == "socat" and any(a.lower().startswith(("exec:", "system:"))
                               or ",exec:" in a.lower() for a in args):
        found.append("socat-exec")
    return found


# Mining markers: `xmrig` and `--donate-level` as whole argv entries
# (so a file named report_xmrig.txt does not trip them); stratum+tcp
# as a prefix so the pool URL counts.
MINING_MARKERS = ("xmrig", "stratum+tcp", "--donate-level")

# World-writable scratch paths a running program should not live in.
TMP_EXE_PREFIXES = ("/tmp/", "/dev/shm/", "/var/tmp/")


def _readlink(path: Path) -> str | None:
    try:
        return os.readlink(path)
    except OSError:
        return None


def _read_text(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _lines(path: Path) -> list[tuple[int, str]]:
    text = _read_text(path)
    if text is None:
        return []
    return list(enumerate(text.splitlines(), start=1))


def _mode(path: Path) -> int | None:
    try:
        return stat.S_IMODE(path.stat().st_mode)
    except OSError:
        return None


def run_read_only_command(argv: list[str]) -> str | None:
    """Run one of the allowed read-only commands; None on any failure.

    Tests monkeypatch this to plant outputs without running anything.
    """
    if not argv or argv[0] not in ("last", "ss", "crontab"):
        raise ValueError(f"command not allowed: {argv!r}")
    try:
        out = subprocess.run(argv, capture_output=True, text=True,
                             timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout if out.returncode == 0 else None


def _proc_scan_uid() -> int:
    """The UID whose processes the scanner owns (its own)."""
    return os.getuid()


def check_processes(proc_root: "Path | str") -> list[Finding]:
    """The own running processes, read-only over a /proc root.

    Only /proc/<pid> of the scanning user's own UID is looked at:
    cmdline and exe of each entry, nothing else. Processes of other
    users, kernel threads (empty cmdline) and unreadable entries are
    skipped silently. Never signals or touches a process.
    """
    proc_root = Path(proc_root)
    uid = _proc_scan_uid()
    out: list[Finding] = []
    try:
        entries = list(os.scandir(proc_root))
    except OSError:
        return out
    for entry in entries:
        if not entry.name.isdigit():
            continue
        try:
            if entry.stat(follow_symlinks=False).st_uid != uid:
                continue
        except OSError:
            continue
        cmdline = _read_text(proc_root / entry.name / "cmdline")
        if not cmdline:
            # kernel threads and unreadable entries: skip silently
            continue
        argv = [a for a in cmdline.split("\0") if a != ""]
        for label in _reverse_shell_shapes(argv):
            out.append(Finding(
                "processes", "alert", f"/proc/{entry.name}", None,
                f"{label} in command line",
                "this is reverse-shell behavior; a shell that "
                "dials out is not something you started on purpose"))
        exe = _readlink(proc_root / entry.name / "exe")
        if exe is not None:
            base = exe[:-len(" (deleted)")] if exe.endswith(" (deleted)") \
                else None
            if base is not None:
                out.append(Finding(
                    "processes", "alert", f"/proc/{entry.name}", None,
                    f"deleted-exe: {base}",
                    "the running program's file has been removed; malware "
                    "often deletes itself after starting"))
            elif exe.startswith(TMP_EXE_PREFIXES):
                out.append(Finding(
                    "processes", "warn", f"/proc/{entry.name}", None,
                    f"tmp-exe: {exe}",
                    "programs run from world-writable scratch directories "
                    "are usually not part of a normal workflow"))
        for arg in argv:
            for marker in MINING_MARKERS:
                if arg == marker or (marker.endswith("+tcp")
                                     and arg.startswith(marker)):
                    out.append(Finding(
                        "processes", "warn", f"/proc/{entry.name}", None,
                        f"mining: {marker}",
                        "a known cryptocurrency-mining marker in the "
                        "process arguments"))
    return out


# --- 1. SSH ---------------------------------------------------------------

def _authorized_key_blobs(text: str) -> list[str]:
    """The key material of each authorized_keys line (no options)."""
    blobs = []
    for line in text.splitlines():
        parts = line.split()
        for i, part in enumerate(parts):
            if part in ("ssh-ed25519", "ssh-rsa", "ecdsa-sha2-nistp256",
                        "ssh-dss") and i + 1 < len(parts):
                blobs.append(parts[i + 1])
                break
    return blobs


def check_ssh(home: Path) -> list[Finding]:
    out: list[Finding] = []
    ssh = home / ".ssh"
    mode = _mode(ssh)
    if mode is not None and (mode & 0o077):
        out.append(Finding("ssh", "warn", str(ssh), None,
                           "directory mode is permissive",
                           "~/.ssh must not be readable or writable by "
                           "group or others"))
    ak = ssh / "authorized_keys"
    text = _read_text(ak)
    if text is not None:
        for n, line in enumerate(text.splitlines(), start=1):
            if not line.strip() or line.startswith("#"):
                continue
            if "command=" in line:
                out.append(Finding(
                    "ssh", "alert", str(ak), n,
                    "authorized key carries a command= option",
                    "a forced command limits this key to one program; "
                    "an attacker uses it to pin a backdoor"))
            if "from=" in line:
                out.append(Finding(
                    "ssh", "info", str(ak), n,
                    "authorized key restricted by from= option",
                    "source restriction on a key; verify it is yours"))
        blobs = _authorized_key_blobs(text)
        for blob in set(blobs):
            if blobs.count(blob) > 1:
                out.append(Finding(
                    "ssh", "warn", str(ak), None,
                    "duplicate key entry in authorized_keys",
                    "the same key is authorized more than once"))
                break
    for name in ("rc", "environment"):
        for n, line in _lines(ssh / name):
            for label, pattern in PERSISTENCE_PATTERNS:
                if pattern.search(line):
                    out.append(Finding(
                        "ssh", "alert", str(ssh / name), n,
                        f"{label} pattern in ~/.ssh/{name}",
                        "start-up hooks that ssh runs on every login"))
    cfg = _read_text(ssh / "config")
    if cfg is not None:
        for label, pattern in SSH_CONFIG_RISKY:
            for m in pattern.finditer(cfg):
                out.append(Finding(
                    "ssh", "alert" if label != "ForwardAgent" else "warn",
                    str(ssh / "config"), cfg[:m.start()].count("\n") + 1,
                    f"{label} in ssh config",
                    f"{label} runs or forwards things automatically; "
                    "confirm every host it names is known"))
    # permissions of private keys: any named key file in ~/.ssh that
    # is group/other accessible and looks like a key
    for entry in _ssh_dir_entries(ssh):
        kmode = _mode(entry)
        if kmode is not None and (kmode & 0o077) and _looks_like_key(entry):
            out.append(Finding(
                "ssh", "warn", str(entry), None,
                "private key file is too open",
                "group/other can read it; ssh refuses, but a copy "
                "left behind still leaks"))
    return out


def _ssh_dir_entries(ssh: Path) -> list[Path]:
    """Named entries only: at most one level, no tree walk."""
    names = ("identity", "id_rsa", "id_dsa", "id_ecdsa", "id_ed25519",
             "id_ecdsa_sk", "id_ed25519_sk")
    return [ssh / n for n in names]


def _looks_like_key(path: Path) -> bool:
    text = _read_text(path)
    return bool(text and text.lstrip().startswith("-----BEGIN"))


# --- 3. persistence -------------------------------------------------------

def check_persistence(home: Path) -> list[Finding]:
    out: list[Finding] = []
    for name in RC_FILES:
        for n, line in _lines(home / name):
            for label, pattern in PERSISTENCE_PATTERNS:
                if pattern.search(line):
                    out.append(Finding(
                        "persistence", "alert", str(home / name), n,
                        f"{label} pattern in {name}",
                        "shell start-up files run on every interactive "
                        "login; download-and-execute here is the classic "
                        "persistence spot"))
    return out


def check_crontab() -> list[Finding]:
    """The user crontab, only through the allowed read command."""
    cron = run_read_only_command(["crontab", "-l"])
    if cron is None:
        return [Finding("persistence", "info", "crontab -l", None,
                        "user crontab not checkable",
                        "crontab -l produced no readable output")]
    out: list[Finding] = []
    for n, line in enumerate(cron.splitlines(), start=1):
        for label, pattern in PERSISTENCE_PATTERNS:
            if pattern.search(line):
                out.append(Finding(
                    "persistence", "alert", "crontab -l", n,
                    f"{label} pattern in user crontab",
                    "a scheduled command that downloads and executes"))
    return out


# --- 6. user binaries -----------------------------------------------------

def check_user_binaries(home: Path, baseline_names: set[str]) -> list[Finding]:
    out: list[Finding] = []
    seen: set[str] = set()
    for rel in BIN_DIRS:
        d = home / rel
        if not d.is_dir():
            continue
        try:
            # one level only, names only
            for entry in sorted(os.scandir(d), key=lambda e: e.name):
                if not entry.is_file():
                    continue
                seen.add(entry.name)
                mode = _mode(Path(entry.path))
                if mode is None:
                    continue
                if entry.name not in baseline_names:
                    out.append(Finding(
                        "user-binaries", "warn",
                        str(Path(entry.path)), None,
                        f"new executable in ~/{rel.as_posix()}: "
                        f"{entry.name}",
                        "a program that appeared since the baseline; "
                        "confirm you installed it"))
                if entry.name.startswith(".") and mode & 0o111:
                    out.append(Finding(
                        "user-binaries", "warn", str(Path(entry.path)), None,
                        f"hidden executable name: {entry.name}",
                        "hidden names dodge a casual glance at the dir"))
                if mode & stat.S_ISUID:
                    out.append(Finding(
                        "user-binaries", "alert", str(Path(entry.path)),
                        None, f"setuid bit on {entry.name}",
                        "setuid user binaries are rare and a strong "
                        "escalation sign"))
                if mode & stat.S_ISGID:
                    out.append(Finding(
                        "user-binaries", "alert", str(Path(entry.path)),
                        None, f"setgid bit on {entry.name}",
                        "setgid user binaries are rare and a strong "
                        "escalation sign"))
        except OSError:
            continue
    return out


def user_binary_names(home: Path) -> set[str]:
    names: set[str] = set()
    for rel in BIN_DIRS:
        d = home / rel
        if not d.is_dir():
            continue
        try:
            for entry in os.scandir(d):
                if entry.is_file():
                    names.add(entry.name)
        except OSError:
            pass
    return names


# --- 7. credentials -------------------------------------------------------

def check_credentials(home: Path) -> list[Finding]:
    out: list[Finding] = []
    for rel in CRED_FILES:
        p = home / rel
        mode = _mode(p)
        if mode is None:
            continue
        if mode != 0o600:
            out.append(Finding(
                "credentials", "warn", str(p), None,
                f"mode {oct(mode)} instead of 0600",
                "credentials files must be readable by the owner only"))
    for rel in HISTORY_FILES:
        p = home / rel
        for n, line in _lines(p):
            for label, pattern in SECRET_PATTERNS:
                if pattern.search(line):
                    out.append(Finding(
                        "credentials", "alert", str(p), n,
                        f"{label} pattern in shell history",
                        "a secret-shaped value sits in the history file; "
                        "the value itself is not reported here"))
    return out


# --- 8. git ---------------------------------------------------------------

def git_hook_names(repo: Path) -> set[str]:
    hooks = repo / ".git" / "hooks"
    if not hooks.is_dir():
        return set()
    try:
        # one level, names only
        return {e.name for e in os.scandir(hooks) if e.is_file()}
    except OSError:
        return set()


def check_git(repos: list[Path], baseline_hooks: set[str]) -> list[Finding]:
    """Hooks in the repos the user configured, plus remote URL changes."""
    out: list[Finding] = []
    for repo in repos:
        hooks = repo / ".git" / "hooks"
        for name in sorted(git_hook_names(repo)):
            p = hooks / name
            if name in baseline_hooks:
                continue
            out.append(Finding(
                "git", "alert", str(p), None,
                f"new git hook: {name}",
                "hooks run on every fetch, commit or push into this repo"))
        cfg = _read_text(repo / ".git" / "config")
        if cfg is not None:
            for n, line in enumerate(cfg.splitlines(), start=1):
                m = re.match(r"\s*url\s*=\s*(\S+)", line)
                if m and m.group(1).startswith("http://"):
                    out.append(Finding(
                        "git", "warn", str(repo / ".git" / "config"), n,
                        "remote URL over plain http",
                        "credentials and code travel unencrypted"))
    return out


# --- 9. DELFIN's own audit log -------------------------------------------

AUDIT_LOG_NAME = re.compile(r"audit(-.*)?\.log$")


def check_delfin_audit(home: Path) -> list[Finding]:
    out: list[Finding] = []
    d = home / ".delfin"
    if not d.is_dir():
        return out
    denied = 0
    outside = 0
    blocklist = 0
    try:
        # fixed dir, one level, only audit-*.log files
        entries = [e for e in os.scandir(d)
                   if e.is_file() and AUDIT_LOG_NAME.search(e.name)]
    except OSError:
        return out
    for entry in entries:
        for n, line in _lines(Path(entry.path)):
            if '"denied"' in line or "'denied'" in line:
                denied += 1
            if "outside the workspace" in line:
                outside += 1
            if "blocklist" in line or "deny-list" in line:
                blocklist += 1
    if denied >= 20:
        out.append(Finding(
            "audit", "warn", "audit-*.log", None,
            f"{denied} denied actions pile up in the audit log",
            "a pile of denials is either a struggling legitimate task "
            "or something probing the gates"))
    if outside:
        out.append(Finding(
            "audit", "warn", "audit-*.log", None,
            f"{outside} access attempts outside the workspace",
            "the agent reached for paths it should not"))
    if blocklist:
        out.append(Finding(
            "audit", "warn", "audit-*.log", None,
            f"{blocklist} blocklisted commands attempted",
            "commands from the deny list were tried"))
    return out


# --- 2. logins ------------------------------------------------------------

_LOGIN_LINE = re.compile(
    r"^\S+\s+\S+\s+(\S+)\s+\w{3}\s+\w{3}\s+\d+\s+(\d{2}):(\d{2})")


def parse_last_output(text: str | None) -> list[Finding]:
    if text is None:
        return [Finding("logins", "info", "last", None,
                        "logins not checkable",
                        "no readable output from last(1)")]
    out: list[Finding] = []
    for n, line in enumerate(text.splitlines(), start=1):
        m = _LOGIN_LINE.match(line)
        if not m:
            continue
        host, hh, mm = m.group(1), int(m.group(2)), int(m.group(3))
        if hh < 5 or hh >= 23:
            out.append(Finding(
                "logins", "warn", "last", n,
                f"login at unusual time {hh:02d}:{mm:02d} from {host}",
                "logins in the small hours deserve a second look"))
    return out


# --- 5. network connections ----------------------------------------------

_SS_LISTEN = re.compile(r"^LISTEN\s+\S+\s+\S+\s+\S+:(\d+)")
_SS_ESTAB = re.compile(
    r"^ESTAB\s+\S+\s+\S+\s+\S+:(\d+)\s+(\S+):(\d+)\s+users:"
)


def parse_ss_output(text: str, baseline_hosts: set[str]) -> list[Finding]:
    out: list[Finding] = []
    for n, line in enumerate(text.splitlines(), start=1):
        m = _SS_LISTEN.match(line)
        if m:
            port = int(m.group(1))
            severity = "alert" if port in SUSPICIOUS_PORTS else "warn"
            out.append(Finding(
                "net", severity, "ss", n,
                f"LISTEN on port {port}",
                "a listening port of your own processes; known one?"))
            if port in SUSPICIOUS_PORTS:
                out.append(Finding(
                    "net", "alert", "ss", n,
                    f"LISTEN on suspicious port {port}",
                    f"port {port} is a classic backdoor port"))
        m = _SS_ESTAB.match(line)
        if m:
            port = int(m.group(3))
            peer = m.group(2)
            if port in SUSPICIOUS_PORTS:
                out.append(Finding(
                    "net", "alert", "ss", n,
                    f"outgoing connection to port {port}",
                    f"port {port} is a classic backdoor port"))
            elif peer not in baseline_hosts:
                out.append(Finding(
                    "net", "info", "ss", n,
                    f"connection to a peer outside the baseline: {peer}",
                    "confirm the destination is yours"))
    return out
