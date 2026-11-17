"""Check a job script for Lustre metadata load BEFORE it is submitted.

The cluster operators measured 43 million open/close calls on the shared
Lustre file systems from one Python start per system out of a Python
environment that lived on Lustre, plus a test suite run straight from a
HOME worktree (operator notes, 14.11. and 16.11.2026). Both came from job
scripts that looked harmless. The rule that fixed it, by hand: unpack the
environment per job to the node-local ``$TMPDIR``, keep scratch on the
node, copy results back once, smoke-test first. This module encodes that
rule so DELFIN checks a script before it goes to the cluster.

:func:`check` is a pure function over the script text: no file-system
walk, no subprocess. Where Lustre is mounted comes from ONE read of
``/proc/mounts`` (fs type ``lustre``), injectable as text via
:class:`MountsSource` for tests -- never a hard-coded cluster name.

Public sources for the rules (the operators' mails are private):
- bwUniCluster wiki, "Workspaces" and "$TMPDIR / scratch":
  https://wiki.bwidm.de/doku.php?id=services:bwunicluster:workspaces
  https://wiki.bwidm.de/doku.php?id=services:bwunicluster:tmpdir
- Lustre best practice: metadata (stat/open/close) load from many small
  operations, Intel HPC glossary / LLNL HPC tutorials
  https://hpc.llnl.gov/banks-jobs/running-jobs/io-best-practices
- NERSC I/O best practice on Lustre: avoid directory walks, keep temp
  on node-local scratch
  https://docs.nersc.gov/filesystems/lustre/best-practices/
"""
from __future__ import annotations

import re
from dataclasses import dataclass

__all__ = ["Finding", "MountsSource", "lustre_paths", "check", "main"]


def _is_lustre(path: str, lustres: list[str]) -> bool:
    return any(path == lp or path.startswith(lp.rstrip("/") + "/") for lp in lustres)


def _path_is_lustre_home(env: dict | None, lustres: list[str]) -> bool:
    home = str((env or {}).get("HOME", ""))
    return bool(home) and _is_lustre(home, lustres)


# Words that say "a whole tree is about to be walked". The person's
# standing rule on the shared Lustre: no tree walks (bwUniCluster wiki,
# NERSC best practice: a find/du/rsync -a over a million files is 43M
# metadata ops). Recursive only: `find -maxdepth 1` is a directory read.
_WALK_CMDS = ("find", "du", "ls", "rsync", "cp")
_WALK_FLAGS = re.compile(
    r"(?:^|\s)-{1,2}[a-zA-Z]*(?:[Rr]|recursive|archive|no-compress)[a-zA-Z]*"
    r"(?:\s|$|=)"
    # rsync -a / cp -a: the GNU/BSD archive flag, recursive by itself.
    r"|(?:^|\s)-a(?:\s|$)"
    r"|(?:(?<=^)|(?<=\s))ls\s+-[a-zA-Z]*R"
)
_MAXDEPTH = re.compile(r"-maxdepth\s+(\d+)")

# Interpreters whose start reads every .so/.pth on their path.
_INTERPRETER = re.compile(
    r"(?:^|\s|/)(?:python3?(?:\d+(?:\.\d+)*)?|pip3?|conda|pytest)\s", re.M)

_LOOP_HEADS = re.compile(r"(?:^|\n)\s*(?:for\s+\S+\s+in|while\s)[^\n]*\n", re.M)

# Scratch settings that would send ORCA/xtb/CREST temp to Lustre
# (bwUniCluster wiki: keep scratch on the node, $TMPDIR).
_SCRATCH_VARS = re.compile(
    r"\b(TMPDIR|XTB_TMPDIR|ORCA_SCRATCH)\b[^\n]*?(/[\w./-]+)"
    r"|(^\s*%scratch\s+/[\w./-]+)", re.M | re.I)


def _strip_env_prefixes(line: str) -> str:
    """Drop leading VAR=val assignments (LC_ALL=C, X=1): they are env
    prefixes, not the command."""
    while True:
        m = re.match(r"\s*[A-Za-z_][\w.]*=(?:\"[^\"]*\"|'[^']*'|\S+)\s+(.*)",
                     line)
        if not m:
            return line
        line = m.group(1)


def _join_continuations(text: str) -> str:
    """Join backslash line continuations so a command and its paths sit
    on one logical line (line numbers then follow the last continued
    line, which is where the walk fires)."""
    out = []
    buf = ""
    for line in text.splitlines():
        buf += line
        if buf.rstrip().endswith("\\"):
            buf = buf.rstrip()[:-1] + " "
            continue
        out.append(buf)
        buf = ""
    if buf:
        out.append(buf)
    return "\n".join(out)


def _expand_env_paths(line: str, env: dict | None) -> list[str]:
    """Absolute paths on the line, with $HOME / ${HOME} (and quoted
    "$HOME/...") expanded from env, and /.. resolved away."""
    env = env or {}
    paths = re.findall(r"(\"\$HOME[^\"]*\"|/[^\s\"';&|]+)", line)
    out = []
    for p in paths:
        p = p.strip("\"")
        home = str(env.get("HOME", ""))
        p = p.replace("${HOME}", home).replace("$HOME", home)
        if "/.." in p:
            from posixpath import normpath
            p = normpath(p)
        out.append(p)
    return out


def _in_loop(lines: list[str], idx: int) -> bool:
    """True when line ``idx`` sits inside a for/while loop body (an
    unterminated loop head above it). Shell and bash arrays: the
    interpreter starts once per item."""
    depth = 0
    for j in range(idx, -1, -1):
        line = lines[j]
        if re.match(r"\s*(for\s+\S+\s+in|while\s)", line):
            return True
        if re.match(r"\s*done\b", line):
            depth += 1
    return False


def _loop_lines(lines: list[str]) -> set[int]:
    """Line indexes (0-based) inside any for/while ... done block, also
    when the whole loop sits on ONE line (``for x in a; do python ...;
    done``)."""
    inside: set[int] = set()
    stack: list[int] = []
    for i, line in enumerate(lines):
        if re.match(r"\s*(for\s+\S+\s+in|while\s)", line):
            if re.search(r"\bdone\b", line):
                inside.add(i)
                continue
            stack.append(i)
        elif re.match(r"\s*done\b", line):
            if stack:
                start = stack.pop()
                for k in range(start, i + 1):
                    inside.add(k)
    return inside


_WALK_WHY = ("a recursive walk over a shared Lustre tree issues millions of "
             "open/stat/close calls and stalls the file system for every user "
             "(NERSC Lustre best practices: avoid directory walks; the "
             "person's standing rule: no tree walks)")
_WALK_INSTEAD = ("run it on the node-local copy or over one directory: "
                 "unpack to $TMPDIR first, or find -maxdepth 1; copy results "
                 "back once at the end")


def check(script_text: str, *, env: dict | None = None,
          mounts_path: MountsSource | None = None) -> list[Finding]:
    """The Lustre metadata-load findings for one job script. Pure: the
    script text, an env dict (HOME on the cluster) and one read of the
    mounts source; no walk, no subprocess. The check itself costs one
    regex pass per rule over the text -- microseconds, no I/O beyond the
    one /proc/mounts read."""
    src = mounts_path or MountsSource("/proc/mounts")
    lustres = lustre_paths(src)
    home_lustre = _path_is_lustre_home(env, lustres)
    if not lustres:
        return []
    text = _join_continuations(script_text)
    lines = text.splitlines()
    inside = _loop_lines(lines)
    findings: list[Finding] = []

    # --- Rule: interpreter / pip / pytest start on Lustre, in a loop ----
    for i, line in enumerate(lines):
        if i not in inside:
            continue
        m = _INTERPRETER.search(line)
        if not m:
            continue
        # The interpreter's code/module path (argv after the word) or,
        # when run bare, the HOME it would read from.
        tail = line[m.end():]
        arg = re.search(r"(/[^\s\"']+)", tail)
        path = arg.group(1) if arg else (str((env or {}).get("HOME", "")))
        if not _is_lustre(path, lustres) and not home_lustre:
            continue
        if re.match(r"\s*(python|pip)", m.group(0).strip()):
            findings.append(Finding(
                "lustre_interpreter_in_loop", i + 1, "warn",
                "an interpreter or pip started once per loop item on Lustre "
                "reads its environment (.pth, shared objects) from the shared "
                "file system every iteration -- the measured 43M open/close "
                "shape (bwUniCluster wiki $TMPDIR guidance)",
                "unpack the Python environment once per job to the node-local "
                "$TMPDIR (tar -x of the env into $TMPDIR; PATH=$TMPDIR/env/bin)"
                " and run it from there"))
        if "pytest" in m.group(0):
            findings.append(Finding(
                "pytest_from_lustre_checkout", i + 1, "refuse",
                "pytest run from a HOME/Lustre checkout stats the whole tree "
                "of tests and the installed environment on the shared file "
                "system (the 16.11.2026 incident)",
                "copy the checkout to the node first: "
                "git archive HEAD | tar -x -C $TMPDIR && cd $TMPDIR && pytest"))

    # --- Rule: pytest from a HOME/Lustre checkout (any line) ------------
    for i, line in enumerate(lines):
        m = _INTERPRETER.search(line)
        if m and "pytest" in m.group(0):
            tail = line[m.end():]
            arg = re.search(r"(/[^\s\"']+)", tail)
            here = re.search(r"(?:cd\s+|;\s*cd\s+)\"?(\$[\w{]+|/[^\\s\"';]+)",
                             script_text)
            # A cd onto $TMPDIR (or another node-local variable) is the
            # FIX this rule asks for: the checkout there is node-local.
            if here and here.group(1).startswith("$"):
                continue
            path = arg.group(1) if arg else (
                here.group(1) if here else str((env or {}).get("HOME", "")))
            if _is_lustre(path, lustres) or home_lustre:
                if not any(f.rule == "pytest_from_lustre_checkout" and f.line == i + 1
                           for f in findings):
                    findings.append(Finding(
                        "pytest_from_lustre_checkout", i + 1, "refuse",
                        "pytest run from a HOME/Lustre checkout stats the whole "
                        "tree of tests and the installed environment on the "
                        "shared file system (the 16.11.2026 incident)",
                        "copy the checkout to the node first: "
                        "git archive HEAD | tar -x -C $TMPDIR && cd $TMPDIR && pytest"))

    # --- Rule: recursive tree walks over Lustre -> refuse ---------------
    # du is recursive by nature; find needs -r/--recursive or -maxdepth>1;
    # ls needs -R; rsync -a / cp -r recurse by their flags. Leading
    # VAR=val words are env prefixes (LC_ALL=C du ...), stripped first.
    for i, line in enumerate(lines):
        cmd = _strip_env_prefixes(line)
        head = re.match(r"\s*(du|find|ls|rsync|cp)\b", cmd)
        if not head:
            continue
        word = head.group(1)
        recursive = (
            word == "du"
            or (word in ("rsync", "cp") and _WALK_FLAGS.search(cmd))
            or (word == "ls" and re.search(r"\s-[a-zA-Z]*R", cmd))
            or (word == "find"
                and not (lambda d: d and int(d.group(1)) <= 1)(
                    _MAXDEPTH.search(cmd))))
        if not recursive:
            continue
        depth = _MAXDEPTH.search(cmd)
        if depth and int(depth.group(1)) <= 1:
            continue
        # Absolute paths on the line, $HOME expanded from env; $TMPDIR,
        # a relative path or the job's own run directory are fine.
        paths = _expand_env_paths(cmd, env)
        if not any(_is_lustre(p, lustres) for p in paths):
            continue
        findings.append(Finding(
            "tree_walk_on_lustre", i + 1, "refuse", _WALK_WHY, _WALK_INSTEAD))

    # --- Rule: scratch on Lustre -> warn --------------------------------
    for i, line in enumerate(lines):
        m = _SCRATCH_VARS.search(line)
        if not m:
            continue
        path = m.group(2) or m.group(3)
        if path and _is_lustre(path, lustres):
            findings.append(Finding(
                "scratch_on_lustre", i + 1, "warn",
                "ORCA/xtb/CREST scratch on the shared Lustre file system "
                "writes and stats millions of small files during the run "
                "(bwUniCluster wiki: keep scratch on the node)",
                "point the scratch at the node-local $TMPDIR "
                "(TMPDIR=$TMPDIR, XTB_TMPDIR=$TMPDIR, %scratch $TMPDIR)"))
    return findings


def _format(findings: list[Finding], name: str) -> str:
    if not findings:
        return f"{name}: no Lustre metadata-load findings."
    out = [f"{name}: {len(findings)} finding(s):"]
    for f in findings:
        out.append(f"  line {f.line} [{f.severity}] {f.rule}: {f.why}")
        out.append(f"           instead: {f.instead}")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    """`python -m delfin.hpc_io_check <script>`: the same check a person
    runs by hand. Prints the findings; exit 1 when one is a refusal."""
    import sys
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print(__doc__.strip().splitlines()[0])
        print("usage: python -m delfin.hpc_io_check <script>")
        return 2 if not args else 0
    name = args[0]
    try:
        with open(name, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError as exc:
        print(f"cannot read {name}: {exc}", file=sys.stderr)
        return 2
    findings = check(text)
    print(_format(findings, name))
    return 1 if any(f.severity == "refuse" for f in findings) else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


@dataclass(frozen=True)
class Finding:
    rule: str
    line: int
    severity: str          # "refuse" | "warn"
    why: str
    instead: str


class MountsSource:
    """Where the mounts come from: a path (default /proc/mounts, read
    once) or literal text (tests)."""

    def __init__(self, source, path: str = "/proc/mounts"):
        self.text = source if isinstance(source, str) else None
        self.path = path

    def read(self) -> str:
        if self.text is not None:
            return self.text
        try:
            with open(self.path, "r", encoding="utf-8", errors="replace") as fh:
                return fh.read()
        except OSError:
            return ""


_LUSTRE_LINE = re.compile(r"^\S+\s+(\S+)\s+lustre\s", re.M)


def lustre_paths(source: MountsSource) -> list[str]:
    """The mount points of fs type ``lustre``, from one read of the
    mounts text. A fixture (tests) or ``/proc/mounts`` (default) -- never
    a hard-coded cluster name."""
    text = source.read() if isinstance(source, MountsSource) else str(source)
    return _LUSTRE_LINE.findall(text)
