"""Whether a calculation may start without a dialog.

Input: a shell command (or the arguments of a submission tool) that starts
a quantum-chemistry program -- ``xtb``, ``crest``, ``orca`` -- or submits a
batch job with ``sbatch``, plus the host it would run on.

Output: a :class:`Verdict`. ``calculation`` says whether the command starts
a calculation at all; ``admit`` whether it may start without asking;
``summary`` what runs where (for the security-event panel), ``reason`` why
it was not admitted (for the dialog, or the refusal when nobody can be
asked).

The rule (user decision 2026-10-22): a calculation runs without a dialog
when its setup passes the check below and its resources are moderate and
compliant with the host. Everything else keeps the dialog.

Setup check
    The input file exists and is read. Charge and spin are taken from where
    the program takes them (xtb/CREST: ``--chrg``/``--uhf``, an xcontrol
    ``$chrg``/``$spin``, the ``.CHRG``/``.UHF`` files; ORCA: the ``* xyz``
    or ``* xyzfile`` line). The electron count of the geometry must allow
    that spin: ``(electrons - unpaired) % 2 == 0``. An element the table
    does not know, a geometry format this module does not read, or an input
    that names no charge/spin line is ``not checked``, which is not a pass.
    ORCA inputs also pass ``delfin.api.validate_orca_input`` without an
    error.

Resources
    ``login_node``: Slurm's client is installed and this process is not in
    a job. Locally at most :data:`LOGIN_LOCAL_CORES` cores and
    :data:`LOGIN_LOCAL_SECONDS` of wall time (a ``timeout`` wrapper, or the
    bash tool's own limit for a foreground run). Anything bigger goes
    through ``sbatch``.

    ``slurm_job``: inside an allocation, its cores
    (``SLURM_CPUS_PER_TASK``/``SLURM_CPUS_ON_NODE``).

    ``workstation``: at most half of the cores this process may use, and
    for ORCA at most half of the memory (``pal * maxcore``).

    ``sbatch``: the requested wall time and nodes within the partition's
    limits (``scontrol show partition``), the cores the job's ORCA input
    asks for within the cores it requests, and at most
    :data:`MAX_QUEUED_JOBS` of the user's jobs pending or running with the
    new one.

Rejected: a regular expression per program on the auto-allow list (what
``xtb <file>.xyz`` had). It admits ``-P 64`` on a login node and an odd
electron count as readily as a correct run.
"""
from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

LOGIN_LOCAL_CORES = 2
LOGIN_LOCAL_SECONDS = 1800
MAX_QUEUED_JOBS = 4

LOCAL_PROGRAMS = ("xtb", "crest", "orca")
_WRAPPERS = ("time", "nice", "nohup", "setsid", "timeout")
_GEOMETRY_SUFFIXES = (".xyz", ".coord", ".sdf", ".mol", ".tmol")


@dataclass
class Verdict:
    calculation: bool = False
    admit: bool = False
    summary: str = ""
    reason: str = ""
    notes: list = field(default_factory=list)


@dataclass
class _Run:
    """One calculation command, taken apart."""
    program: str
    argv: list
    env: dict
    wall_s: Optional[float]
    detached: bool


# --------------------------------------------------------------------------
# Host
# --------------------------------------------------------------------------

def host_kind(env: Optional[dict] = None) -> str:
    """``slurm_job``, ``login_node`` or ``workstation``."""
    env = os.environ if env is None else env
    if str(env.get("SLURM_JOB_ID", "")).strip():
        return "slurm_job"
    if shutil.which("sbatch") and shutil.which("squeue"):
        return "login_node"
    return "workstation"


def _available_cpus() -> int:
    from .bash_jobs import _available_cpus as _cpus
    return max(1, int(_cpus()))


def local_core_limit(kind: str) -> int:
    if kind == "login_node":
        return LOGIN_LOCAL_CORES
    if kind == "slurm_job":
        return _available_cpus()
    return max(1, _available_cpus() // 2)


def default_threads(kind: str, env: Optional[dict] = None) -> int:
    """The thread count a program gets when its command names none.

    The same number the agent's shell environment sets
    (``_scrubbed_bash_env``): the allocation inside a job, the login-node
    limit on a login node, else at most four.
    """
    env = os.environ if env is None else env
    cap = str(env.get("SLURM_CPUS_PER_TASK", "")).strip()
    if cap.isdigit() and int(cap) > 0:
        return int(cap)
    if kind == "login_node":
        return LOGIN_LOCAL_CORES
    return max(1, min(4, os.cpu_count() or 1))


def _memory_mb() -> Optional[int]:
    try:
        with open("/proc/meminfo", encoding="ascii") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        return None
    return None


# --------------------------------------------------------------------------
# Chemistry: electrons and spin
# --------------------------------------------------------------------------

def _atomic_number(symbol: str) -> Optional[int]:
    sym = re.sub(r"[^A-Za-z]", "", symbol or "")
    if not sym or len(sym) > 3:
        return None
    sym = sym[0].upper() + sym[1:].lower()
    try:
        from rdkit import Chem
        z = int(Chem.GetPeriodicTable().GetAtomicNumber(sym))
    except Exception:
        return None
    return z if z > 0 else None


def electrons_of(symbols) -> Optional[int]:
    """Sum of atomic numbers, or None when any symbol is unknown."""
    total = 0
    n = 0
    for s in symbols:
        z = _atomic_number(s)
        if z is None:
            return None
        total += z
        n += 1
    return total if n else None


def _xyz_symbols(text: str) -> Optional[list]:
    lines = [l for l in text.splitlines()]
    if not lines:
        return None
    try:
        count = int(lines[0].split()[0])
        body = lines[2:2 + count]
    except (ValueError, IndexError):
        body = [l for l in lines if l.strip()]
        count = len(body)
    symbols = []
    for line in body:
        parts = line.split()
        if len(parts) < 4:
            return None
        symbols.append(parts[0])
    return symbols if symbols and len(symbols) == count else None


def _coord_symbols(text: str) -> Optional[list]:
    out, inside = [], False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("$coord"):
            inside = True
            continue
        if inside and s.startswith("$"):
            break
        if inside and s:
            parts = s.split()
            if len(parts) < 4:
                return None
            out.append(parts[3])
    return out or None


def _molfile_symbols(text: str) -> Optional[list]:
    lines = text.splitlines()
    if len(lines) < 4 or "V2000" not in lines[3]:
        return None
    try:
        n = int(lines[3][:3])
    except ValueError:
        return None
    out = []
    for line in lines[4:4 + n]:
        parts = line.split()
        if len(parts) < 4:
            return None
        out.append(parts[3])
    return out or None


def geometry_symbols(path: Path) -> Optional[list]:
    """Element symbols of a geometry file, or None when it cannot be read."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    name = path.name.lower()
    if name.endswith(".xyz"):
        return _xyz_symbols(text)
    if name.endswith((".sdf", ".mol")):
        return _molfile_symbols(text)
    if name.endswith((".coord", ".tmol")) or name == "coord" or "$coord" in text:
        return _coord_symbols(text)
    return None


def spin_allowed(electrons: int, unpaired: int) -> bool:
    return unpaired >= 0 and electrons >= unpaired and (
        electrons - unpaired) % 2 == 0


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------

def _duration_s(text: str) -> Optional[float]:
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([smhd]?)", text.strip())
    if not m:
        return None
    return float(m.group(1)) * {"": 1, "s": 1, "m": 60, "h": 3600,
                                "d": 86400}[m.group(2)]


def _parse_run(segment: str) -> Optional[_Run]:
    """The calculation in one shell segment, or None if it starts none."""
    text = segment.strip()
    detached = text.endswith("&") and not text.endswith("&&")
    if detached:
        text = text[:-1].strip()
    try:
        words = shlex.split(text)
    except ValueError:
        return None
    env, wall = {}, None
    while words and re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", words[0]):
        k, v = words.pop(0).split("=", 1)
        env[k] = v
    while words and words[0] in _WRAPPERS:
        w = words.pop(0)
        if w in ("nohup", "setsid"):
            detached = True
        elif w == "time":
            if words and words[0] == "-p":
                words.pop(0)
        elif w == "nice":
            if words and words[0] == "-n" and len(words) > 1:
                del words[:2]
            elif words and re.fullmatch(r"-\d+", words[0]):
                words.pop(0)
        elif w == "timeout":
            while words and words[0].startswith("-"):
                opt = words.pop(0)
                if opt in ("-s", "-k", "--signal", "--kill-after") and words:
                    words.pop(0)
            if not words:
                return None
            wall = _duration_s(words.pop(0))
    if not words:
        return None
    program = os.path.basename(words[0])
    if program not in LOCAL_PROGRAMS + ("sbatch",):
        return None
    argv = []
    skip = False
    for w in words[1:]:
        if skip:
            skip = False
            continue
        if re.fullmatch(r"\d*>>?|\d*<|&>>?|>&|<<<?", w):
            skip = True
            continue
        if re.match(r"^(\d*>>?|\d*<|&>>?)\S", w) or re.fullmatch(r"\d*>&\d+", w):
            continue
        argv.append(w)
    return _Run(program, argv, env, wall, detached)


def _flag_value(argv: list, *names: str) -> Optional[str]:
    for i, w in enumerate(argv):
        for n in names:
            if w == n and i + 1 < len(argv):
                return argv[i + 1]
            if n.startswith("--") and w.startswith(n + "="):
                return w.split("=", 1)[1]
    return None


def _int(value) -> Optional[int]:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _resolve(cwd: Path, name: str) -> Path:
    """Lexical: a directory a preceding `mkdir` makes does not exist yet."""
    p = Path(os.path.expanduser(name))
    return Path(os.path.normpath(p if p.is_absolute() else cwd / p))


# xtb and CREST options that take a value, so it is not read as the input.
_XTB_VALUED = (
    "-c", "--chrg", "-u", "--uhf", "-P", "--parallel", "-T", "--threads",
    "--gfn", "--alpb", "--gbsa", "--cosmo", "--cpcmx", "--input", "-I",
    "--namespace", "--cycles", "--acc", "--etemp", "--ewin", "--temp",
    "--len", "--mdlen", "--tstep", "--cinp", "--json", "--copy",
    "--iterations", "--vparam", "--xparam", "--rthr", "--ethr", "--bthr",
)
_XTB_LEVELS = {"crude", "sloppy", "loose", "lax", "normal", "tight",
               "vtight", "extreme"}


def _geometry_argument(argv: list, cwd: Path) -> Optional[Path]:
    """The first positional argument that is an existing geometry file."""
    i = 0
    while i < len(argv):
        w = argv[i]
        if w.startswith("-"):
            i += 2 if w in _XTB_VALUED else 1
            continue
        if w.lower() in _XTB_LEVELS:
            i += 1
            continue
        return _resolve(cwd, w)
    return None


def _first_int(path: Path) -> Optional[int]:
    try:
        words = path.read_text(encoding="utf-8", errors="replace").split()
    except OSError:
        return None
    return _int(words[0]) if words else None


def _xcontrol_charge_spin(path: Path) -> tuple:
    charge = unpaired = None
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None, None
    m = re.search(r"^\s*\$chrg\s+(-?\d+)", text, re.M)
    if m:
        charge = int(m.group(1))
    m = re.search(r"^\s*\$spin\s+(\d+)", text, re.M)
    if m:
        unpaired = int(m.group(1))
    return charge, unpaired


def _check_xtb_like(run: _Run, cwd: Path, kind: str) -> tuple:
    """(cores, problem) for an xtb or CREST command."""
    geom = _geometry_argument(run.argv, cwd)
    if geom is None:
        return 0, f"{run.program}: no input file named on the command line"
    if not geom.is_file():
        return 0, f"{run.program}: input {geom.name} does not exist"
    if geom.suffix.lower() in (".toml",):
        return 0, f"{run.program}: a TOML input is not checked here"
    symbols = geometry_symbols(geom)
    if symbols is None:
        return 0, (f"{run.program}: the geometry {geom.name} could not be "
                   "read (format not checked here)")
    electrons = electrons_of(symbols)
    if electrons is None:
        return 0, f"{run.program}: {geom.name} contains an unknown element"
    charge = _int(_flag_value(run.argv, "--chrg", "-c"))
    unpaired = _int(_flag_value(run.argv, "--uhf", "-u"))
    xc = _flag_value(run.argv, "--input", "-I")
    if xc:
        c2, u2 = _xcontrol_charge_spin(_resolve(cwd, xc))
        charge = c2 if charge is None else charge
        unpaired = u2 if unpaired is None else unpaired
    # xtb reads .CHRG and .UHF from the directory it runs in.
    if charge is None:
        charge = _first_int(cwd / ".CHRG")
    if unpaired is None:
        unpaired = _first_int(cwd / ".UHF")
    charge = 0 if charge is None else charge
    unpaired = 0 if unpaired is None else unpaired
    n = electrons - charge
    if not spin_allowed(n, unpaired):
        return 0, (f"{run.program}: {geom.name} has {n} electrons at charge "
                   f"{charge}, which cannot have {unpaired} unpaired "
                   "electron(s) (set --uhf or .UHF to match the parity)")
    cores = (_int(_flag_value(run.argv, "-P", "--parallel", "-T",
                              "--threads"))
             or _int(run.env.get("OMP_NUM_THREADS"))
             or default_threads(kind))
    return cores, ""


_ORCA_GEOM_RE = re.compile(
    r"^[ \t]*\*[ \t]*(xyz|xyzfile|int|gzmtfile|pdbfile)[ \t]+(-?\d+)[ \t]+"
    r"(\d+)[ \t]*(\S*)",
    re.I | re.M)


def orca_resources(text: str) -> tuple:
    """(cores, maxcore MB) an ORCA input asks for."""
    cores = 1
    m = re.search(r"%\s*pal\b.*?\bnprocs\s+(\d+)", text, re.I | re.S)
    if m:
        cores = int(m.group(1))
    for kw in re.findall(r"^\s*!.*$", text, re.M):
        for m2 in re.finditer(r"\bPAL(\d+)\b", kw, re.I):
            cores = max(cores, int(m2.group(1)))
    m = re.search(r"%\s*maxcore\s+(\d+)", text, re.I)
    maxcore = int(m.group(1)) if m else None
    return cores, maxcore


def check_orca_input(path: Path) -> tuple:
    """(cores, maxcore, problem) for one ORCA input file."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 0, None, f"orca: input {path.name} does not exist"
    try:
        from delfin.api import validate_orca_input
        errors = [i for i in validate_orca_input(text)
                  if getattr(i, "severity", "") == "error"]
    except Exception:
        errors = []
    if errors:
        return 0, None, ("orca: " + path.name + ": "
                         + "; ".join(str(getattr(e, "message", e))
                                     for e in errors[:3]))
    m = _ORCA_GEOM_RE.search(text)
    if not m:
        return 0, None, (f"orca: {path.name} names no '* xyz' or "
                         "'* xyzfile' line, so charge and spin are not "
                         "checked here")
    kind, charge, mult, ref = m.group(1).lower(), int(m.group(2)), \
        int(m.group(3)), m.group(4)
    if kind == "xyz":
        block = text[m.end():]
        end = re.search(r"^\s*\*\s*$", block, re.M)
        lines = (block[:end.start()] if end else "").splitlines()
        symbols = [l.split()[0] for l in lines if len(l.split()) >= 4]
    elif kind == "xyzfile":
        symbols = geometry_symbols(path.parent / ref) if ref else None
    else:
        symbols = None
    if not symbols:
        return 0, None, (f"orca: the geometry of {path.name} could not be "
                         "read (format not checked here)")
    electrons = electrons_of(symbols)
    if electrons is None:
        return 0, None, f"orca: {path.name} contains an unknown element"
    n = electrons - charge
    if mult < 1 or not spin_allowed(n, mult - 1):
        return 0, None, (f"orca: {path.name} has {n} electrons at charge "
                         f"{charge}, which cannot have multiplicity {mult}")
    cores, maxcore = orca_resources(text)
    return cores, maxcore, ""


def _orca_input_argument(run: _Run, cwd: Path) -> Optional[Path]:
    for w in run.argv:
        if not w.startswith("-"):
            return _resolve(cwd, w)
    return None


def _check_local(run: _Run, cwd: Path, kind: str, tool_wall_s: float,
                 ) -> tuple:
    """(summary, problem) for a calculation run on this host."""
    maxcore = None
    if run.program == "orca":
        inp = _orca_input_argument(run, cwd)
        if inp is None:
            return "", "orca: no input file named"
        cores, maxcore, problem = check_orca_input(inp)
        name = inp.name
    else:
        cores, problem = _check_xtb_like(run, cwd, kind)
        geom = _geometry_argument(run.argv, cwd)
        name = geom.name if geom else "?"
    if problem:
        return "", problem
    limit = local_core_limit(kind)
    where = {"login_node": "on this login node", "slurm_job":
             "in this Slurm allocation", "workstation": "on this machine"}[kind]
    if cores > limit:
        hint = (" Submit it with sbatch instead." if kind == "login_node"
                else "")
        return "", (f"{run.program} {name} asks for {cores} cores; at most "
                    f"{limit} run {where} without asking.{hint}")
    if run.program == "orca" and maxcore:
        total = _memory_mb()
        share = (total // 2) if total else None
        if kind == "login_node":
            share = min(share or 8000, LOGIN_LOCAL_CORES * 4000)
        if share and cores * maxcore > share:
            return "", (f"orca {name} asks for {cores * maxcore} MB "
                        f"({cores} x %maxcore {maxcore}); at most {share} MB "
                        f"run {where} without asking.")
    if kind == "login_node":
        wall = run.wall_s
        if wall is None and not run.detached:
            wall = tool_wall_s
        if wall is None or wall > LOGIN_LOCAL_SECONDS:
            return "", (f"{run.program} {name} has no wall-time limit of "
                        f"{LOGIN_LOCAL_SECONDS // 60} min; on a login node "
                        f"a local run needs one (`timeout "
                        f"{LOGIN_LOCAL_SECONDS} {run.program} ...`) or goes "
                        "through sbatch.")
    return f"{run.program} {name}, {cores} core(s) {where}", ""


# --------------------------------------------------------------------------
# sbatch
# --------------------------------------------------------------------------

def slurm_seconds(text: str) -> Optional[int]:
    """Slurm's time formats: M, M:S, H:M:S, D-H, D-H:M, D-H:M:S."""
    try:
        return _slurm_seconds(text)
    except (TypeError, ValueError):
        return None


def _slurm_seconds(text: str) -> Optional[int]:
    t = (text or "").strip()
    if not t or t.upper() in ("UNLIMITED", "INFINITE"):
        return None
    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        days = int(d)
        parts = [int(p) for p in t.split(":")]
        while len(parts) < 3:
            parts.append(0)
        h, m, s = parts[:3]
    else:
        parts = [int(p) for p in t.split(":")]
        if len(parts) == 1:
            h, m, s = 0, parts[0], 0
        elif len(parts) == 2:
            h, m, s = 0, parts[0], parts[1]
        else:
            h, m, s = parts[:3]
    return ((days * 24 + h) * 60 + m) * 60 + s


_SBATCH_VALUED = {"-t", "--time", "-N", "--nodes", "-n", "--ntasks",
                  "-c", "--cpus-per-task", "-p", "--partition", "-D",
                  "--chdir", "-J", "--job-name", "-o", "--output", "-e",
                  "--error", "-A", "--account", "--mem", "--mem-per-cpu",
                  "-q", "--qos", "--gres", "-C", "--constraint"}

_SBATCH_SHORT = {"-t": "--time", "-N": "--nodes", "-n": "--ntasks",
                 "-c": "--cpus-per-task", "-p": "--partition",
                 "-D": "--chdir"}


def sbatch_options(script_text: str, argv: list) -> dict:
    """``#SBATCH`` directives, overridden by command-line options."""
    opts: dict = {}

    def take(words):
        i = 0
        while i < len(words):
            w = words[i]
            key, val = w, None
            if "=" in w and w.startswith("--"):
                key, val = w.split("=", 1)
            if val is None and key in _SBATCH_VALUED and i + 1 < len(words):
                val = words[i + 1]
                i += 1
            key = _SBATCH_SHORT.get(key, key)
            if key.startswith("--"):
                opts[key] = val
            i += 1

    for line in script_text.splitlines():
        s = line.strip()
        if s.startswith("#SBATCH"):
            try:
                take(shlex.split(s[len("#SBATCH"):], comments=True))
            except ValueError:
                pass
        elif s and not s.startswith("#"):
            break
    take(list(argv))
    return opts


def _run_quiet(cmd: list, timeout: float = 10.0) -> Optional[str]:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=timeout)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout if out.returncode == 0 else None


def partition_limits(name: Optional[str]) -> Optional[dict]:
    """``scontrol show partition`` as key=value pairs, or None."""
    out = _run_quiet(["scontrol", "show", "partition", "-o"]
                     + ([name] if name else []))
    if not out:
        return None
    rows = []
    for line in out.splitlines():
        pairs = dict(p.split("=", 1) for p in line.split() if "=" in p)
        if pairs:
            rows.append(pairs)
    if name:
        return rows[0] if rows else None
    for r in rows:
        if r.get("Default", "").upper() == "YES":
            return r
    return None


def queued_jobs() -> Optional[int]:
    out = _run_quiet(["squeue", "-h", "-u", os.environ.get("USER", ""),
                      "-t", "PENDING,RUNNING", "-o", "%i"])
    if out is None:
        return None
    return len([l for l in out.splitlines() if l.strip()])


def _script_calculations(text: str, cwd: Path) -> list:
    """ORCA/xtb/CREST lines of a job script, with the directory they run in."""
    runs = []
    here = cwd
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        for seg in re.split(r"&&|;|\|\|", line):
            seg = seg.strip()
            m = re.fullmatch(r"cd\s+(\S+)", seg)
            if m:
                target = m.group(1).strip("\"'")
                if target in ("$SLURM_SUBMIT_DIR", "${SLURM_SUBMIT_DIR}"):
                    here = cwd
                elif "$" not in target:
                    here = _resolve(here, target)
                else:
                    here = None
                continue
            run = _parse_run(seg)
            if run is not None and run.program in LOCAL_PROGRAMS:
                runs.append((run, here))
    return runs


def check_sbatch(run: _Run, cwd: Path,
                 limits_of: Callable = partition_limits,
                 jobs_now: Callable = queued_jobs) -> tuple:
    """(summary, problem) for an sbatch submission."""
    # sbatch's own options come first; the first word that is no option and
    # no option's value is the script, and the words after it are its own.
    script, before, i = None, [], 0
    while i < len(run.argv):
        w = run.argv[i]
        if w.startswith("-"):
            before.append(w)
            if "=" not in w and w in _SBATCH_VALUED and i + 1 < len(run.argv):
                before.append(run.argv[i + 1])
                i += 1
            i += 1
            continue
        script = _resolve(cwd, w)
        break
    if any(w.startswith("--wrap") for w in before):
        return "", "sbatch: a --wrap command is not checked here"
    if script is None or not script.is_file():
        return "", "sbatch: no job script file named"
    text = script.read_text(encoding="utf-8", errors="replace")
    opts = sbatch_options(text, before)
    chdir = opts.get("--chdir")
    job_cwd = _resolve(cwd, chdir) if chdir and "$" not in chdir else cwd
    limits = limits_of(opts.get("--partition"))
    if limits is None:
        return "", ("sbatch: the partition's limits could not be read "
                    "(scontrol), so the request is not checked")
    want_s = slurm_seconds(opts.get("--time") or "") if opts.get(
        "--time") else None
    max_s = slurm_seconds(limits.get("MaxTime", ""))
    if want_s is not None and max_s is not None and want_s > max_s:
        return "", (f"sbatch: --time {opts['--time']} exceeds the partition "
                    f"limit {limits.get('MaxTime')}")
    nodes = _int((opts.get("--nodes") or "1").split("-")[0]) or 1
    max_nodes = _int(limits.get("MaxNodes", ""))
    if max_nodes and nodes > max_nodes:
        return "", (f"sbatch: {nodes} nodes exceed the partition limit "
                    f"{max_nodes}")
    tasks = _int(opts.get("--ntasks")) or 1
    per_task = _int(opts.get("--cpus-per-task")) or 1
    job_cores = tasks * per_task
    max_cpn = _int(limits.get("MaxCPUsPerNode", ""))
    if max_cpn and job_cores > max_cpn * nodes:
        return "", (f"sbatch: {job_cores} cores exceed the partition limit "
                    f"of {max_cpn} per node")
    calcs = _script_calculations(text, job_cwd)
    for crun, here in calcs:
        if here is None:
            return "", ("sbatch: the job script changes into a directory "
                        "named by a variable, so its inputs are not checked")
        if crun.program == "orca":
            inp = _orca_input_argument(crun, here)
            if inp is None:
                return "", "sbatch: the job's orca line names no input"
            cores, _mc, problem = check_orca_input(inp)
            if not problem and cores > job_cores:
                problem = (f"orca {inp.name} asks for {cores} cores and the "
                           f"job requests {job_cores}")
        else:
            cores, problem = _check_xtb_like(crun, here, "slurm_job")
        if problem:
            return "", "sbatch: " + problem
    if not calcs:
        return "", ("sbatch: the job script runs no orca, xtb or crest line "
                    "this check can read")
    queued = jobs_now()
    if queued is None:
        return "", "sbatch: the queue could not be read (squeue)"
    if queued + 1 > MAX_QUEUED_JOBS:
        return "", (f"sbatch: {queued} of your jobs are pending or running; "
                    f"at most {MAX_QUEUED_JOBS} are submitted without asking")
    part = limits.get("PartitionName", opts.get("--partition") or "default")
    return (f"sbatch {script.name}: {job_cores} core(s), "
            f"{opts.get('--time') or 'default time'} on {part}"), ""


# --------------------------------------------------------------------------
# The verdict
# --------------------------------------------------------------------------

def judge_command(segments: list, cwd: Path,
                  segment_allowed: Callable[[str], bool],
                  *, tool_wall_s: Optional[float] = None,
                  kind: Optional[str] = None,
                  limits_of: Callable = partition_limits,
                  jobs_now: Callable = queued_jobs) -> Verdict:
    """The verdict for a shell command given as its segments.

    ``segment_allowed`` is the gate's own auto-allow test for a segment
    that starts no calculation (``cd``, ``mkdir``...); a compound admits
    only when every other segment passes it.
    """
    kind = kind or host_kind()
    here: Optional[Path] = Path(cwd)
    summaries, problems, found = [], [], False
    for seg in segments:
        run = _parse_run(seg)
        if run is None:
            m = re.fullmatch(r"\s*cd\s+(\S+)\s*", seg)
            if m:
                target = m.group(1).strip("\"'")
                here = None if "$" in target or here is None else \
                    _resolve(here, target)
            elif not segment_allowed(seg):
                problems.append(f"'{seg.strip()[:60]}' is not a calculation "
                                "and is not auto-allowed")
            continue
        found = True
        if here is None:
            problems.append(f"{run.program}: runs in a directory named by a "
                            "variable, so its input is not checked")
            continue
        if run.program == "sbatch":
            summary, problem = check_sbatch(run, here, limits_of, jobs_now)
        else:
            summary, problem = _check_local(run, here, kind, tool_wall_s)
        if problem:
            problems.append(problem)
        else:
            summaries.append(summary)
    if not found:
        return Verdict()
    if problems:
        return Verdict(True, False, "", "; ".join(problems))
    return Verdict(True, True, "; ".join(summaries), "")


def judge_orca_submission(folder: Path, pal: int, maxcore: int,
                          time_limit: str, *, kind: Optional[str] = None,
                          limits_of: Callable = partition_limits,
                          jobs_now: Callable = queued_jobs) -> Verdict:
    """The verdict for the ops ``submit_calculation`` in ORCA mode."""
    kind = kind or host_kind()
    inputs = sorted(Path(folder).glob("*.inp"))
    if len(inputs) != 1:
        return Verdict(True, False, "", (
            f"submit_calculation: {len(inputs)} .inp files in the folder; "
            "exactly one is checked"))
    _cores, _mc, problem = check_orca_input(inputs[0])
    if problem:
        return Verdict(True, False, "", problem)
    if kind == "workstation":
        limit = local_core_limit(kind)
        if pal > limit:
            return Verdict(True, False, "", (
                f"submit_calculation: pal {pal} exceeds {limit} cores on "
                "this machine"))
        return Verdict(True, True,
                       f"{inputs[0].name}: {pal} core(s) on this machine", "")
    limits = limits_of(None)
    if limits is None:
        return Verdict(True, False, "", (
            "submit_calculation: the partition's limits could not be read"))
    want_s, max_s = slurm_seconds(time_limit), slurm_seconds(
        limits.get("MaxTime", ""))
    if want_s is not None and max_s is not None and want_s > max_s:
        return Verdict(True, False, "", (
            f"submit_calculation: time {time_limit} exceeds the partition "
            f"limit {limits.get('MaxTime')}"))
    max_cpn = _int(limits.get("MaxCPUsPerNode", ""))
    if max_cpn and pal > max_cpn:
        return Verdict(True, False, "", (
            f"submit_calculation: pal {pal} exceeds {max_cpn} per node"))
    queued = jobs_now()
    if queued is None or queued + 1 > MAX_QUEUED_JOBS:
        return Verdict(True, False, "", (
            "submit_calculation: the queue could not be read" if queued is None
            else f"submit_calculation: {queued} of your jobs are pending or "
                 f"running; at most {MAX_QUEUED_JOBS} without asking"))
    return Verdict(True, True, (f"{inputs[0].name}: {pal} core(s), "
                                f"{time_limit} via Slurm"), "")


def judge_orca_run(input_file: Path, *, kind: Optional[str] = None) -> Verdict:
    """The verdict for the ops ``run_orca_input`` (a local ORCA run)."""
    run = _Run("orca", [str(input_file)], {}, None, False)
    kind = kind or host_kind()
    # The ops tool has no wall-time limit of its own.
    summary, problem = _check_local(run, Path(input_file).parent, kind, None)
    if problem:
        return Verdict(True, False, "", problem)
    return Verdict(True, True, summary, "")
