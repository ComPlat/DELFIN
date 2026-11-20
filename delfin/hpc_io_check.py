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

Public sources for the rules (the operators' mails are private; each
page checked this session via DELFIN's web tools):
- bwHPC wiki, bwUniCluster 2.0/3.0 Hardware and Architecture (workspaces,
  node-local scratch): https://wiki.bwhpc.de/e/BwUniCluster2.0/Hardware_and_Architecture
- INCD "Lustre best practices": writing/accessing thousands of files
  loads the metadata servers, can take the file system offline:
  https://wiki.incd.pt/books/filesystem-user-guide/page/lustre-best-practices
- NAS HECC "Lustre Best Practices" (shared metadata contention):
  https://www.nas.nasa.gov/hecc/support/kb/lustre-best-practices_226.html
- University of Helsinki Lustre User Guide: avoid plain ``find`` and
  long ``ls`` — metadata operation load on the MDS:
  https://wiki.helsinki.fi/xwiki/bin/view/it4sci/IT%20for%20Science%20group/HPC%20Environment%20User%20Guide/Lustre%20User%20Guide/
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

__all__ = ["Finding", "MountsSource", "lustre_paths", "check",
           "check_submit", "main"]

# The mounts file the production check reads by default (L1 fix: the
# default wiring must READ this file, never treat the path as text).
# Tests monkeypatch this constant.
DEFAULT_MOUNTS = "/proc/mounts"


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

# Interpreters whose start reads every .so/.pth on their path. The word
# must be complete: it ends at whitespace OR at the end of the line --
# a bare `pytest` as the last word of a line (the commonest job-script
# spelling, wave-33 review F1) was invisible when only \s counted.
_INTERPRETER = re.compile(
    r"(?:^|\s|/)(?:python3?(?:\d+(?:\.\d+)*)?|pip3?|conda|pytest)(?=\s|$)",
    re.M)

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


def _join_continuations(text: str) -> tuple[str, list[int]]:
    """Join backslash line continuations so a command and its paths sit
    on one logical line. Returns ``(joined, source_lines)`` where
    ``source_lines[i]`` is the 0-based SOURCE line index of logical line
    ``i`` -- a report must name the line the user can open (F5, wave-33
    review): a walk on real line 5 behind continuations was reported as
    line 4."""
    out = []
    source: list[int] = []
    buf = ""
    buf_start = 0
    for idx, line in enumerate(text.splitlines()):
        if not buf:
            buf_start = idx
        buf += line
        if buf.rstrip().endswith("\\"):
            buf = buf.rstrip()[:-1] + " "
            continue
        out.append(buf)
        source.append(buf_start)
        buf = ""
    if buf:
        out.append(buf)
        source.append(buf_start)
    return "\n".join(out), source


def _strip_wrappers(cmd: str) -> str:
    """Drop wrappers that do not change what runs (wave-33 review F4 +
    W1/W2): leading `bash -c`/`sh -c` (their QUOTED argument is the
    command; trailing $0 names are dropped), `command`, `nice` with or
    without its own adjustment (-n 10, -10), `env` with or without
    assignments / -i (re-stripped inside the loop), `nohup`, and one
    level of grouping parentheses. The command word must then anchor
    the rule; before this a wrapper hid the walk from it."""
    while True:
        cmd = _strip_env_prefixes(cmd)
        m = re.match(
            r"\s*(?:bash|sh)\s+-c\s+(?:'([^']*)'|\"([^\"]*)\")(?:\s.*)?$",
            cmd)
        if m:
            cmd = m.group(1) or m.group(2) or ""
            continue
        m = re.match(r"\s*(?:command|nohup|exec)\s+(.*)", cmd)
        if m:
            cmd = m.group(1)
            continue
        # nice with or without its own arguments: nice -n 10 cmd,
        # nice -10 cmd, nice cmd.
        m = re.match(r"\s*nice(?:\s+-n\s+\S+|\s+-\d+)?\s+(.*)", cmd)
        if m and m.group(1):
            cmd = m.group(1)
            continue
        # env: env X=1 ... / env -i ... / env ...  -- assignments are
        # stripped by _strip_env_prefixes on the next loop pass.
        m = re.match(r"\s*env(?:\s+-i|\s+(?:-\w+\s+)*[A-Za-z_]\w*=(?:\S+|\"[^\"]*\"|'[^']*'))*\s+(.+)",
                     cmd)
        if m and m.group(1):
            cmd = m.group(1)
            continue
        m = re.match(r"\s*\(\s*(.*)\s*\)\s*$", cmd)
        if m and ("&&" in m.group(1) or ";" not in m.group(1)):
            # A single command or a chain inside parentheses; a lone
            # `( cmd )` in a script is one command.
            cmd = m.group(1)
            continue
        return cmd


def _collect_vars(text: str, env: dict | None) -> dict[str, str]:
    """Every variable a script assigns itself plus the caller's env --
    the map `$VAR` expansion needs so a variable cannot hide a Lustre
    path (wave-33 review F3). Script assignments win: they run."""
    merged: dict[str, str] = {}
    merged.update(env or {})
    for m in re.finditer(
            r"(?:^|\n)\s*(?:export\s+)?([A-Za-z_]\w*)=(?:\"([^\"]*)\"|'([^']*)'|(\S*))",
            text):
        merged[m.group(1)] = m.group(2) or m.group(3) or m.group(4) or ""
    return merged


def _expand_env_paths(line: str, env: dict | None) -> list[str]:
    """Absolute paths on the line, with $VAR / ${VAR} expanded from the
    script's own assignments and the caller's env (F3: only $HOME was
    expanded, so any other variable hid a Lustre path completely),
    `~/` from HOME, and /.. resolved away."""
    env = env or {}
    paths = re.findall(
        r"(~/(?:[\w./-]*)|\"\$\{?\w+\}?[^\"]*\"|\$\{?\w+\}?|/[^\s\"';&|]+)",
        line)
    out = []
    for p in paths:
        p = p.strip("\"")
        if p.startswith("~/"):
            home = str(env.get("HOME", ""))
            p = home + p[1:] if home else p
        # Expand every $VAR / ${VAR}, repeatedly (nested assignments),
        # resolving from the merged map.
        for _ in range(3):
            expanded = re.sub(r"\$\{(\w+)\}|\$(\w+)",
                              lambda m: str(env.get(m.group(1) or m.group(2),
                                                    m.group(0))),
                              p)
            if expanded == p:
                break
            p = expanded
        if "/.." in p:
            from posixpath import normpath
            p = normpath(p)
        if p and not p.startswith("$"):
            out.append(p)
    return out


def _source_is_wide_root(path: str, env: dict | None,
                         lustres: list[str]) -> bool:
    """True when the copy SOURCE is HOME, a workspace root or a mount
    point (a whole-tree walk), not one job folder (operator live-run L5:
    the stage-in of one folder is the fix; walking HOME or a workspace
    root is still the problem)."""
    home = str((env or {}).get("HOME", ""))
    wide = [home] + lustres if home else list(lustres)
    for root in wide:
        if not root:
            continue
        r = root.rstrip("/")
        if path == r:
            return True
        # One path level below a wide root is still a wide root
        # (HOME/projects, /pfs/data6/workspace) -- a JOB FOLDER sits
        # deeper (runs/water under a workspace).
        rest = path[len(r):].lstrip("/")
        if rest and "/" not in rest:
            return True
    return False


def _dest_is_node_local(cmd: str) -> bool:
    """For cp/rsync: is the LAST path argument a node-local destination
    ($TMPDIR, $TMP, /tmp, /scratch, /dev/shm, $SCRATCH)? Such a copy is
    the stage-in the rules ask for (operator live-run L5); the walk
    refusal is for source-side reads over Lustre."""
    # Strip flags so the last token is a path, not "-a" or "--archive".
    tokens = [t for t in cmd.split()
              if not (t.startswith("-") and len(t) > 1)]
    # Drop the command word itself.
    tokens = tokens[1:]
    if not tokens:
        return False
    dest = tokens[-1].strip("'\"")
    if dest.startswith("$"):
        var = re.match(r"\$\{?(\w+)\}?", dest).group(1)
        return var.upper() in ("TMPDIR", "TMP", "SCRATCH")
    if not dest.startswith("/"):
        return False
    return bool(re.match(r"^/(tmp|scratch|dev/shm)(/.*)?$", dest))


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
    src = mounts_path or MountsSource(path=DEFAULT_MOUNTS)
    lustres = lustre_paths(src)
    home_lustre = _path_is_lustre_home(env, lustres)
    if not lustres:
        return []
    text, source_lines = _join_continuations(script_text)
    lines = text.splitlines()
    inside = _loop_lines(lines)
    findings: list[Finding] = []
    vars_map = _collect_vars(script_text, env)

    # A helper the line-numbered rules use: the 1-based SOURCE line of
    # logical line ``i`` (F5: continuations collapsed the count, so a
    # walk on real line 5 was reported as line 4).

    def _src(i: int) -> int:
        return (source_lines[i] if 0 <= i < len(source_lines) else i) + 1

    def _cd_is_node_local(cd_target: str) -> bool:
        """True when the cd target is node-local: a $VAR, or a path
        that is not on Lustre either as written or expanded from the
        script's own variables."""
        if cd_target.startswith("$"):
            return True
        if _is_lustre(cd_target, lustres):
            return False
        expanded = _expand_env_paths(cd_target, vars_map)
        resolved = expanded[0] if expanded else cd_target
        return not _is_lustre(resolved, lustres)

    # --- Rule: interpreter / pip / pytest start on Lustre, in a loop ----
    for i, line in enumerate(lines):
        if i not in inside:
            continue
        m = _INTERPRETER.search(line)
        if not m:
            continue
        # L3 fix: the match may START with the leading "/" of an
        # absolute interpreter path (m.group(0) = "/python "); classify
        # by the interpreter WORD, which is always last in the match.
        interp_word = m.group(0).strip().rstrip("/").split("/")[-1].strip()
        # The interpreter's code/module path (argv after the word) or,
        # when run bare, the HOME it would read from.
        tail = line[m.end():]
        arg = re.search(r"(/[^\s\"']+)", tail)
        path = arg.group(1) if arg else (str((env or {}).get("HOME", "")))
        if not _is_lustre(path, lustres) and not home_lustre:
            continue
        if re.match(r"(python|pip)", interp_word):
            findings.append(Finding(
                "lustre_interpreter_in_loop", _src(i), "warn",
                "an interpreter or pip started once per loop item on Lustre "
                "reads its environment (.pth, shared objects) from the shared "
                "file system every iteration -- the measured 43M open/close "
                "shape (bwUniCluster wiki $TMPDIR guidance)",
                "unpack the Python environment once per job to the node-local "
                "$TMPDIR (tar -x of the env into $TMPDIR; PATH=$TMPDIR/env/bin)"
                " and run it from there"))
        if "pytest" in interp_word:
            findings.append(Finding(
                "pytest_from_lustre_checkout", _src(i), "refuse",
                "pytest run from a HOME/Lustre checkout stats the whole tree "
                "of tests and the installed environment on the shared file "
                "system (the 16.11.2026 incident)",
                "copy the checkout to the node first: "
                "git archive HEAD | tar -x -C $TMPDIR && cd $TMPDIR && pytest"))

    # --- Rule: pytest from a HOME/Lustre checkout (any line) ------------
    # L4 fix: `python -m pytest` (any interpreter path) is the same
    # shape as bare `pytest` -- the module runs the whole test tree.
    for i, line in enumerate(lines):
        m = _INTERPRETER.search(line)
        is_pytest_shape = (
            (m and "pytest" in m.group(0))
            or bool(re.search(
                r"(?:^|\s)[-\w./]*(?:python3?\b|pytest\b)[^\n]*\s-m\s+pytest\b",
                line)))
        if m and is_pytest_shape:
            tail = line[m.end():]
            arg = re.search(r"(/[^\s\"']+)", tail)
            # The cd that counts is the LAST one on the pytest's OWN
            # line (a chained `cd a && cd b` is two cds on one line:
            # the shell ends up in b), or — when the line has none —
            # the last cd on the lines before it (W3: taking any
            # earlier line's cd missed `cd /tmp && cd /pfs/...` and
            # falsely excused/accused on the chain direction).
            own = list(re.finditer(r"(?:^|&&|;)\s*cd\s+\"?([^&;\n]+)",
                                   line))
            if own:
                cd_target = own[-1].group(1).strip().strip('"\'')
                if _cd_is_node_local(cd_target):
                    continue
                # The shell really ends on a Lustre dir: refuse.
                path = cd_target
            else:
                text_before = "\n".join(lines[:i])
                here = None
                for cd_m in re.finditer(
                        r"(?:^|\n|&&|;)\s*cd\s+\"?([^&;\n]+)", text_before):
                    here = cd_m
                if here is not None:
                    cd_target = here.group(1).strip().strip('"\'')
                    if _cd_is_node_local(cd_target):
                        continue
                    path = cd_target
                else:
                    path = str(vars_map.get("HOME", ""))
            path = (_expand_env_paths(path, vars_map)[0]
                    if path.startswith("$") else path)
            if _is_lustre(path, lustres) or home_lustre or (
                    not arg and path == str(vars_map.get("HOME", ""))
                    and _is_lustre(path, lustres)):
                if not any(f.rule == "pytest_from_lustre_checkout" and f.line == _src(i)
                           for f in findings):
                    findings.append(Finding(
                        "pytest_from_lustre_checkout", _src(i), "refuse",
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
        cmd = _strip_wrappers(line)
        head = re.match(r"\s*(du|find|ls|rsync|cp)\b", cmd)
        # A single-line loop ('for i in 1 2; do du ...; done') hides the
        # command behind the loop head; also try the post-'do' segment
        # so the start-anchored head match does not defeat the rule.
        seg = None
        if not head and " do " in cmd:
            seg = cmd.split(" do ", 1)[1]
            head = re.match(r"\s*(du|find|ls|rsync|cp)\b", seg)
        if not head:
            continue
        if seg is not None:
            cmd = seg
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
        # Absolute paths on the line, every $VAR expanded from the
        # script's own assignments and env (F3), $TMPDIR, a relative
        # path or the job's own run directory are fine.
        paths = _expand_env_paths(cmd, vars_map)
        # L5 fix: a stage-in -- one job folder copied TO a node-local
        # destination -- is the fix the rule asks for, not a walk. But a
        # copy whose SOURCE is HOME, a workspace root or a mount point
        # still walks a whole tree on Lustre and stays refused.
        if word in ("cp", "rsync") and _dest_is_node_local(cmd):
            # The destination ($TMPDIR/...) expands to no absolute Lustre
            # path, so paths holds only the source side: judge the
            # source directly. HOME, a workspace root or a mount point
            # stays a refusal; one job folder deeper is a stage-in.
            if not any(_source_is_wide_root(p, env, lustres) for p in paths):
                continue
        if not any(_is_lustre(p, lustres) for p in paths):
            continue
        findings.append(Finding(
            "tree_walk_on_lustre", _src(i), "refuse", _WALK_WHY, _WALK_INSTEAD))

    # --- Rule: scratch on Lustre -> warn --------------------------------
    for i, line in enumerate(lines):
        m = _SCRATCH_VARS.search(line)
        if not m:
            continue
        path = m.group(2) or m.group(3)
        if path and _is_lustre(path, lustres):
            findings.append(Finding(
                "scratch_on_lustre", _src(i), "warn",
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


_SCRIPT_SUFFIXES = (".sh", ".bash", ".sub", ".py")


def check_submit(script, env: dict | None = None,
                 mounts_path: MountsSource | None = None
                 ) -> tuple[list[str], list[str]]:
    """The shared submit-time check: ONE helper both submit routes call.

    ``script`` is a path to a submit script (Path or str) or the script
    TEXT itself. Returns ``(refused, warned)`` -- finding strings in the
    exact form ``cluster_tools.submit`` has reported them since IO1
    (``name:line [rule] why instead: how``, cluster_tools.py:403-404), so
    the dashboard backend, the agent's bash gate and ``cluster_submit``
    cannot drift apart in what they call a refusal. ``name`` is the
    script's basename; for text passed inline it is the first word of
    the shebang or simply ``<script>``.

    ``env`` defaults to ``os.environ`` (H4): with ``None`` the HOME rule
    could not fire and a caller that forgot env would silently weaken
    the check. Pure-text callers that must not inherit the process env
    pass ``env={}`` explicitly.

    Cost: one text read (a path) plus the regex passes of :func:`check`;
    the mounts source is one /proc/mounts read per call, injectable for
    tests. Raises OSError when a script path cannot be read -- the check
    never fails open (H1: a str that names a script-shaped file but
    misses on disk raises like the Path spelling; it is never re-read as
    script text, which would be a clean verdict for a script nobody
    read). Callers may catch OSError and refuse the submit on their own
    wording; backend_slurm catches FileNotFoundError and leaves the
    missing script to sbatch's own error.
    """
    if isinstance(script, Path):
        text = script.read_text(encoding="utf-8", errors="replace")
        name = script.name
    elif (isinstance(script, str) and "\n" not in script
          and script.endswith(_SCRIPT_SUFFIXES)):
        # A script-shaped str is a PATH (H1): if it names nothing on
        # disk that is an error, never a silent reinterpretation as
        # one-line script text.
        text = Path(script).read_text(encoding="utf-8", errors="replace")
        name = Path(script).name
    elif (isinstance(script, str) and "\n" not in script
          and Path(script).is_file()):
        text = Path(script).read_text(encoding="utf-8", errors="replace")
        name = Path(script).name
    else:
        text = str(script)
        name = "<script>"
    if env is None:
        import os as _os
        env = dict(_os.environ)
    findings = check(text, env=env, mounts_path=mounts_path)

    def _fmt(f: Finding) -> str:
        return f"{name}:{f.line} [{f.rule}] {f.why} instead: {f.instead}"

    refused = [_fmt(f) for f in findings if f.severity == "refuse"]
    warned = [_fmt(f) for f in findings if f.severity == "warn"]
    return refused, warned


def sbatch_gate(cmd: str, env: dict | None = None,
                mounts_path: MountsSource | None = None
                ) -> tuple[str, str, bool]:
    """The decision hook the agent's bash gate calls for a command that
    runs sbatch (IO2 phase 3, the px in ``api_client.py``).

    Returns ``(verdict, text, may_run)``:

    - ``verdict``: ``"allow"`` when the command is not an sbatch
      command or its script is clean; ``"refuse"`` when the script
      carries a refuse finding; ``"warn"`` when it only carries warn
      findings.
    - ``text``: the finding text for refuse/warn, ``""`` for allow.
    - ``may_run``: False only for refuse -- the command must not
      execute; warn and allow run.

    An sbatch command is one whose segment starts with the word
    ``sbatch`` after any of ``;`` ``&&`` ``|`` ``env X=1`` ``exec``
    ``nice`` ``command`` ``nohup`` (the same wrapper shapes
    :func:`_strip_wrappers` handles). The script to check is the LAST
    positional argument (an ``-o/-o``-style attached token is an
    option, not the script); ``--wrap=<text>`` or ``--wrap <text>`` is
    checked as script text instead. When neither a script argument nor
    --wrap is present the verdict is allow -- sbatch reads its script
    from stdin then, and there is nothing to inspect.

    Pure over text; the same mounts-source rules as :func:`check`.
    """
    import shlex
    # The LAST segment the shell would run, stripped of wrappers: the
    # gate asks about the sbatch command the shell actually executes
    # (a `cd /tmp && sbatch ...` runs sbatch from /tmp).
    segments = re.split(r"&&|\|\||;|\|", cmd)
    target = _strip_wrappers(segments[-1] if segments else cmd)
    tokens = target.split()
    if not tokens or tokens[0] != "sbatch":
        return "allow", "", True
    if env is None:
        import os as _os
        env = dict(_os.environ)
    # --wrap=<text> or --wrap <text>: the wrapped text IS the script.
    wrap = re.search(r"--wrap(?:=(?:'([^']*)'|\"([^\"]*)\"|(\S+))|\s+('([^']*)'|\"([^\"]*)\"|(\S+)))",
                     target)
    if wrap:
        text = next(g for g in wrap.groups() if g is not None)
        refused, warned = check_submit(text, env=env,
                                       mounts_path=mounts_path)
    else:
        # The script is the last token that is not an option or an
        # option's value. A token starting with '-' is an option; a
        # leading '--' token with '=' carries its value attached. An
        # ATTACHED short-option value (`-o<value>`, reviewer W-shape)
        # is not the script either: the value after the single-letter
        # flag is stripped before the positional scan.
        stripped = []
        for t in tokens[1:]:
            m = re.match(r"^-([A-Za-z])(\S.*)$", t)
            if m and m.group(1) in "oe":  # -o outfile, -e errfile
                continue  # attached option VALUE, not a positional
            stripped.append(t)
        positional = [t for t in stripped if not t.startswith("-")]
        if not positional:
            # The whole command line was options (e.g. `sbatch
            # -o<value>` only): no batch script is named, sbatch reads
            # stdin. But an attached value that LOOKS like a script
            # path is checked conservatively -- if it exists and would
            # be refused, refuse (a wrong -o path is never a walk, so
            # this can only fire on a real script path).
            for t in tokens[1:]:
                m = re.match(r"^-([A-Za-z])(\S.*)$", t)
                if m and m.group(1) in "oe" and Path(m.group(2)).is_file():
                    script_arg = m.group(2)
                    break
            else:
                return "allow", "", True
        else:
            script_arg = positional[-1]
        try:
            refused, warned = check_submit(script_arg, env=env,
                                           mounts_path=mounts_path)
        except FileNotFoundError:
            # sbatch will fail on its own with "No such file"; the gate
            # does not add a second error for a path that names nothing.
            return "allow", "", True
        except OSError:
            # A token like `-o<missing>` was taken as the script because
            # an attached option value is indistinguishable from a bare
            # path by shape alone; if it names nothing the command goes
            # to sbatch, which refuses it. Same rule as FileNotFoundError.
            return "allow", "", True
    if refused:
        return ("refuse",
                "hpc_io_check refused: the sbatch script would load the "
                "shared Lustre file system. " + " ".join(refused)
                + " Fix the script and submit again.", False)
    if warned:
        return "warn", "hpc_io_check warn: " + " ".join(warned), True
    return "allow", "", True


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


@dataclass(frozen=True)
class Finding:
    rule: str
    line: int
    severity: str          # "refuse" | "warn"
    why: str
    instead: str


class MountsSource:
    """Where the mounts come from. Two explicit forms (operator live-run
    L1 fix): ``MountsSource(path)`` or ``MountsSource(path=...)`` reads
    that FILE; ``MountsSource(text=...)`` takes literal mounts text
    (tests). A bare string is a PATH, never content -- the L1 bug was
    the path string being stored as mounts text, which disabled the
    production check entirely."""

    def __init__(self, source=None, *, path: str | None = None,
                 text: str | None = None):
        given = [x for x in (source, path, text) if x is not None]
        if len(given) != 1:
            raise TypeError(
                "MountsSource takes exactly one of: a positional path, "
                "path=..., or text=... (L1 contract)")
        if text is not None:
            self.text = text
            self.path = None
        else:
            self.text = None
            self.path = str(source if source is not None else path)

    def read(self) -> str:
        if self.text is not None:
            return self.text
        try:
            with open(self.path, "r", encoding="utf-8", errors="replace") as fh:
                return fh.read()
        except OSError as exc:
            # QS finding (IO1, nacht-s8): an unreadable mounts file used to
            # yield "" -> no Lustre paths -> the check silently disabled
            # itself. The check may never fail open: raise instead.
            raise RuntimeError(
                f"hpc_io_check: cannot read mounts file {self.path!r}: {exc}"
            ) from exc


_LUSTRE_LINE = re.compile(r"^\S+\s+(\S+)\s+lustre\s", re.M)


def lustre_paths(source: MountsSource) -> list[str]:
    """The mount points of fs type ``lustre``, from one read of the
    mounts text. A fixture (tests) or ``/proc/mounts`` (default) -- never
    a hard-coded cluster name."""
    text = source.read()
    return _LUSTRE_LINE.findall(text)


# The CLI entry point sits at the very END of the module, below every
# definition (operator live-run L2: above them it raised NameError on
# python -m delfin.hpc_io_check).
if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
