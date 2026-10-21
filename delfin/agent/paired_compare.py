"""A paired comparison of two runs that refuses to be read as more than it is.

Two commands back it (``delfin-agent provenance`` and ``delfin-agent
compare-arms``); both only read and print, so the shell gate treats them
as read-only and no tool schema is sent with every request.

Input of :func:`provenance`: the code paths that make up the instrument,
optionally the arm's command line (its leading ``VAR=value``, ``env`` and
``export`` assignments are overlaid on this shell's environment) or the
pid of a runner that is already running (its environment, start time and
working directory are read from the process). Output: a stamp dict --
content hash of the code paths, git head and a hash of the uncommitted
diff of the checkout, platform and interpreter, the environment variables
that can change behaviour (configured prefixes plus the determinism
variables), and for a pid the code files modified after the process
started.

Input of :func:`compare`: two arm result files (JSONL or a JSON list, one
record per item), their stamps, the declared axis (environment variable
names and/or ``code``/``platform``), the item key and the fields that say
ok / defect / quality. Output: a verdict dict whose ``state`` is one of
:data:`STATES`, the absolute counts on the paired population, every
never-worse case by item, and a stratified sample for a person to look
at. No state says "land": the strongest is ``needs_human_check``.

Semantics, in the order the checks run:

* Same population: duplicate keys within an arm, or item sets that differ
  between the arms, make the pair ``skewed_pair``. A rate over a
  population that moved is not compared.
* One instrument state: platform, interpreter or code hash that differ
  outside the declared axis make the pair ``skewed_pair``.
* Axis guard: every behaviour variable whose effective value differs
  between the arms must be declared; an undeclared one is ``skewed_pair``.
  A declared variable with the same effective value in both arms (a
  switch at its default) is ``null_test_suspected``. A declared variable
  set in one arm only, with no ``--default`` to resolve it, is
  ``axis_unknown``: null -> value is a new measurement.
* Shared switches: a behaviour variable set to the same value in both
  arms and not declared is refused (``undeclared_switch``) unless the
  caller opts in (``allow_shared_switches`` or the environment variable
  :data:`ALLOW_SHARED_ENV`). A value containing "/" is a data path and an
  empty value is dropped; every other value -- numbers included -- is a
  switch.
* Reach probe: when the first ``probe`` items (key order of arm A) have
  identical content in both arms, or no item changed at all, the state is
  ``null_test_suspected``. Content is the record with volatile fields
  (time, duration, pid, host, ...) removed.
* Never-worse terms, any single case blocks (``blocked``), each listed:
  ok -> not ok, a defect field that rose on an item, a quality value that
  fell on an item, and items whose content differs between two runs of
  arm A (``rerun_a``: nondeterminism). Absolute defect totals that rise
  block as well.
* Minimum n: fewer paired items than ``min_n`` is ``insufficient_n``;
  nothing is concluded automatically, the user decides.
* Otherwise ``no_difference`` or ``needs_human_check``.

Rejected: a tool in the model's schema (about 200 tokens in every
request of every session for a check most sessions never run); a
significance test as the deciding signal (one never-worse case blocks
regardless of p); deriving the behaviour prefixes from the axis name only
(the field case was a variable of the same project under another name);
treating 0/false as "off" (whether 0 is off depends on the code's
default, so only a declared ``--default`` resolves it).

Limits: an environment read from ``/proc/<pid>/environ`` is the process's
START environment; a process that changes its own environment later is
not seen. Off Linux a pid's environment is unavailable and the stamp says
so, which makes the comparison ``axis_unknown``.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shlex
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

SCHEMA = "delfin.provenance/1"

#: Verdict states, from "no comparison" to "a person should look".
STATES = ("skewed_pair", "undeclared_switch", "null_test_suspected",
          "axis_unknown", "blocked", "insufficient_n", "no_difference",
          "needs_human_check")

#: Opt-in for behaviour variables set identically in both arms.
ALLOW_SHARED_ENV = "DELFIN_COMPARE_ALLOW_SHARED_SWITCHES"

#: Behaviour prefixes taken in every stamp besides the checkout's own
#: (``<CHECKOUT NAME>_``, which is ``DELFIN_`` in DELFIN's repository) and
#: the ones the caller names. Empty: in another project ``DELFIN_*`` are
#: the agent's own settings, not the project's behaviour.
DEFAULT_PREFIXES: tuple[str, ...] = ()

#: Variables outside any project prefix that change results or their
#: determinism.
DETERMINISM_VARS = (
    "PYTHONHASHSEED", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    "CUBLAS_WORKSPACE_CONFIG", "CUDA_VISIBLE_DEVICES", "TF_DETERMINISTIC_OPS",
    "LC_ALL", "LANG", "TZ",
)

#: Variables the agent harness sets in every shell it starts; they say
#: which session ran the arm, not how the arm behaves.
HARNESS_VARS = frozenset({"DELFIN_AGENT_SESSION", ALLOW_SHARED_ENV})

_SECRET_RE = re.compile(r"KEY|TOKEN|SECRET|PASSW|CREDENTIAL|AUTH",
                        re.IGNORECASE)
_VOLATILE_KEY_RE = re.compile(
    r"(?:^|_)(?:time|times|timestamp|ts|duration|elapsed|wall|took|pid|"
    r"host|hostname|node|date|started|ended|start|end|seconds|ms|runtime)"
    r"(?:_s|_ms|_sec)?$|_at$", re.IGNORECASE)
_SKIP_DIRS = frozenset({".git", "__pycache__", ".pytest_cache", ".mypy_cache",
                        "node_modules", ".tox", ".venv"})
_KEY_CANDIDATES = ("id", "key", "item", "item_id", "task_id", "name", "case")
_OFF_VALUES = frozenset({""})
_MAX_LISTED = 40
_START_SLACK_S = 1.0


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def _redact(name: str, value: str) -> str:
    return "<redacted>" if _SECRET_RE.search(name) else value


def checkout_prefix(path: str | os.PathLike = ".") -> str:
    """``<NAME>_`` for the git checkout holding *path* (its top directory
    name upper-cased, every non-alphanumeric run as ``_``), or ""."""
    top = _git(path, "rev-parse", "--show-toplevel")
    if not top:
        return ""
    name = re.sub(r"[^A-Za-z0-9]+", "_", Path(top).name).strip("_").upper()
    return f"{name}_" if name else ""


def behaviour_env(env: Mapping[str, str], prefixes: Iterable[str]) -> dict:
    """The variables of *env* that can change how an arm behaves: every
    name starting with one of *prefixes*, plus :data:`DETERMINISM_VARS`.
    Secret-named values are replaced by ``<redacted>``."""
    pre = tuple(p for p in prefixes if p)
    out = {}
    for name, value in env.items():
        if name in DETERMINISM_VARS or (pre and name.startswith(pre)):
            out[name] = _redact(name, str(value))
    return dict(sorted(out.items()))


_ASSIGN_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)=(.*)$", re.DOTALL)


def command_env(command: str, base: Mapping[str, str]) -> dict:
    """*base* with the assignments the arm's *command* makes before (or
    instead of) its program: leading ``VAR=value`` words, ``env [-u VAR]
    VAR=value`` and ``export VAR=value`` segments. Later segments see the
    exports of earlier ones; the assignments of a segment's own prefix
    are the ones in force for that segment's program, and the program
    segment that is last wins (it is the arm). Unparsable text leaves
    *base* unchanged."""
    env = dict(base)
    try:
        segments = re.split(r"&&|\|\||;|\n", command or "")
    except Exception:
        return env
    effective = dict(env)
    for seg in segments:
        try:
            words = shlex.split(seg)
        except ValueError:
            continue
        if not words:
            continue
        if words[0] == "export":
            for w in words[1:]:
                m = _ASSIGN_RE.match(w)
                if m:
                    env[m.group(1)] = m.group(2)
            effective = dict(env)
            continue
        if words[0] == "unset":
            for w in words[1:]:
                env.pop(w, None)
            effective = dict(env)
            continue
        local = dict(env)
        i = 0
        while i < len(words):
            w = words[i]
            m = _ASSIGN_RE.match(w)
            if m:
                local[m.group(1)] = m.group(2)
                i += 1
                continue
            if w == "env":
                i += 1
                while i < len(words):
                    w2 = words[i]
                    if w2 == "-u" and i + 1 < len(words):
                        local.pop(words[i + 1], None)
                        i += 2
                        continue
                    if w2 in ("-i", "--ignore-environment"):
                        local = {}
                        i += 1
                        continue
                    m2 = _ASSIGN_RE.match(w2)
                    if m2:
                        local[m2.group(1)] = m2.group(2)
                        i += 1
                        continue
                    break
                continue
            break
        if i < len(words):
            effective = local
        else:
            # A segment of assignments only sets shell variables, which a
            # later program does not inherit unless exported.
            effective = dict(env)
    return effective


def _git(path: str | os.PathLike, *args: str) -> str:
    from . import safe_git
    r = safe_git.run(args, path, timeout=20)
    return r.stdout.strip() if r.returncode == 0 else ""


def _iter_files(paths: Sequence[str]) -> list[Path]:
    out: list[Path] = []
    for raw in paths:
        p = Path(raw)
        if p.is_file():
            out.append(p)
        elif p.is_dir():
            for root, dirs, files in os.walk(p):
                dirs[:] = sorted(d for d in dirs if d not in _SKIP_DIRS)
                for f in sorted(files):
                    if f.endswith((".pyc", ".pyo")):
                        continue
                    out.append(Path(root) / f)
    return out


def code_hash(paths: Sequence[str]) -> dict:
    """sha256 over every file under *paths* (name and bytes, sorted), the
    file count and the newest modification time. A path that does not
    exist is listed under ``missing`` and never silently skipped."""
    h = hashlib.sha256()
    files = _iter_files(paths)
    newest = 0.0
    for f in sorted(files, key=lambda p: str(p.resolve())):
        try:
            data = f.read_bytes()
            st = f.stat()
        except OSError:
            continue
        h.update(str(f.resolve()).encode("utf-8", "replace") + b"\0")
        h.update(len(data).to_bytes(8, "big") + data)
        newest = max(newest, st.st_mtime)
    missing = [p for p in paths if not Path(p).exists()]
    return {"paths": [str(Path(p).resolve()) for p in paths],
            "hash": h.hexdigest()[:20], "files": len(files),
            "newest_mtime": newest, "missing": missing}


def _etime_seconds(text: str) -> Optional[float]:
    """``ps -o etime`` ([[dd-]hh:]mm:ss) in seconds."""
    m = re.fullmatch(r"\s*(?:(\d+)-)?(?:(\d+):)?(\d+):(\d+)\s*", text or "")
    if not m:
        return None
    d, h, mi, s = (int(x) if x else 0 for x in m.groups())
    return float(((d * 24 + h) * 60 + mi) * 60 + s)


def process_start_epoch(pid: int) -> Optional[float]:
    """When *pid* started, as epoch seconds; None when it cannot be told.
    /proc on Linux (``stop_all.process_started_at``, the reader the
    emergency stop uses), ``ps -o etime`` elsewhere."""
    from .stop_all import process_started_at
    started = process_started_at(int(pid))
    if started is not None:
        return started
    try:
        r = subprocess.run(["ps", "-o", "etime=", "-p", str(int(pid))],
                           capture_output=True, text=True, timeout=5)
    except Exception:
        return None
    secs = _etime_seconds(r.stdout) if r.returncode == 0 else None
    return None if secs is None else time.time() - secs


def _process_env(pid: int) -> Optional[dict]:
    try:
        raw = Path(f"/proc/{int(pid)}/environ").read_bytes()
    except (OSError, ValueError):
        return None
    env = {}
    for item in raw.split(b"\0"):
        if b"=" in item:
            k, _, v = item.partition(b"=")
            env[k.decode("utf-8", "replace")] = v.decode("utf-8", "replace")
    return env


def _process_cwd(pid: int) -> str:
    try:
        return os.readlink(f"/proc/{int(pid)}/cwd")
    except (OSError, ValueError):
        return ""


def _process_cmdline(pid: int) -> str:
    try:
        raw = Path(f"/proc/{int(pid)}/cmdline").read_bytes()
        return " ".join(shlex.quote(p.decode("utf-8", "replace"))
                        for p in raw.split(b"\0") if p)
    except (OSError, ValueError):
        pass
    try:
        r = subprocess.run(["ps", "-o", "command=", "-p", str(int(pid))],
                           capture_output=True, text=True, timeout=5)
        return r.stdout.strip() if r.returncode == 0 else ""
    except Exception:
        return ""


def provenance(*, code: Sequence[str] = (), command: str = "",
               pid: Optional[int] = None, prefixes: Sequence[str] = (),
               cwd: str = "", env: Optional[Mapping[str, str]] = None,
               now: Optional[float] = None) -> dict:
    """The stamp of one arm. See the module docstring."""
    where = cwd or os.getcwd()
    source = "self"
    process: dict = {}
    base_env: Optional[Mapping[str, str]] = env if env is not None else os.environ
    if pid is not None:
        penv = _process_env(pid)
        process = {"pid": int(pid), "start_epoch": process_start_epoch(pid),
                   "cmdline": _process_cmdline(pid)}
        pcwd = _process_cwd(pid)
        if pcwd:
            where = cwd or pcwd
            process["cwd"] = pcwd
        if penv is None:
            source, base_env = "unavailable", None
        else:
            source, base_env = "pid", penv
    elif command:
        source = "command"
        base_env = command_env(command, base_env or {})
    pre = list(DEFAULT_PREFIXES) + [p for p in prefixes if p]
    own = checkout_prefix(where)
    if own:
        pre.append(own)
    pre = list(dict.fromkeys(pre))
    top = _git(where, "rev-parse", "--show-toplevel")
    diff = _git(where, "diff", "HEAD") if top else ""
    stamp: dict[str, Any] = {
        "schema": SCHEMA,
        "taken_at": float(now if now is not None else time.time()),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "cwd": str(where),
        "git": {"toplevel": top,
                "head": _git(where, "rev-parse", "HEAD") if top else "",
                "dirty_hash": (hashlib.sha256(diff.encode()).hexdigest()[:20]
                               if diff else "")},
        "prefixes": pre,
        "env_source": source,
        "env": behaviour_env(base_env, pre) if base_env is not None else None,
        "command": command,
    }
    if code:
        stamp["code"] = code_hash(code)
    if process:
        stamp["process"] = process
        start = process.get("start_epoch")
        if code and start:
            # /proc/stat's btime is whole seconds, so the start is known to
            # within one second; a file written in that second is not
            # called later than the process.
            stale = [str(f) for f in _iter_files(code)
                     if _mtime(f) > float(start) + _START_SLACK_S]
            stamp["code_changed_after_start"] = stale[:_MAX_LISTED]
    return stamp


_STAMP_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}")


def stamp_dir() -> Path:
    """Where saved stamps live: ``~/.delfin/provenance``."""
    return Path.home() / ".delfin" / "provenance"


def save_stamp(name: str, stamp: Mapping) -> Path:
    """Write *stamp* once under :func:`stamp_dir` as ``<name>.json``.

    The name is one plain word (no path), and an existing name is refused:
    a saved stamp is the record of an arm's instrument at that moment, and
    re-stamping under the same name would let a later state stand in for
    it. Raises ``ValueError`` on a bad or taken name."""
    if not _STAMP_NAME_RE.fullmatch(name or "") or ".." in name:
        raise ValueError(f"a stamp name is one word of letters, digits, "
                         f"'.', '_' or '-': {name!r}")
    from .state_paths import ensure_dir, secure_file
    d = ensure_dir(stamp_dir())
    path = d / f"{name}.json"
    body = dict(stamp)
    body["sealed_name"] = name
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError:
        raise ValueError(f"stamp {name!r} exists already ({path}); a stamp "
                         "is written once -- choose a new name") from None
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(body, fh, indent=1, sort_keys=True)
    secure_file(path)
    return path


def load_stamp(ref: str) -> Optional[dict]:
    """A stamp by saved name, else by file path; None for "". A stamp read
    from the store carries ``sealed: True``; one read from any other file
    ``sealed: False`` (a file the session could have edited)."""
    if not ref:
        return None
    if _STAMP_NAME_RE.fullmatch(ref) and (stamp_dir() / f"{ref}.json").is_file():
        data = json.loads((stamp_dir() / f"{ref}.json").read_text("utf-8"))
        data["sealed"] = True
        return data
    data = json.loads(Path(ref).read_text(encoding="utf-8"))
    data["sealed"] = False
    return data


def _mtime(p: Path) -> float:
    try:
        return p.stat().st_mtime
    except OSError:
        return 0.0


# ---------------------------------------------------------------------------
# reading the arms
# ---------------------------------------------------------------------------

def load_records(path: str | os.PathLike) -> list[dict]:
    """Records of one arm: a JSON list, a JSON object with one list value
    (``{"items": [...]}``), or JSON lines. Non-dict entries are dropped."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    stripped = text.lstrip()
    data: Any = None
    if stripped.startswith(("[", "{")):
        try:
            data = json.loads(text)
        except ValueError:
            data = None
    if isinstance(data, dict):
        lists = [v for v in data.values() if isinstance(v, list)]
        data = lists[0] if len(lists) == 1 else [data]
    if isinstance(data, list):
        return [r for r in data if isinstance(r, dict)]
    out = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except ValueError:
            continue
        if isinstance(rec, dict):
            out.append(rec)
    return out


def pick_key(records: Sequence[dict], key: str = "") -> str:
    if key:
        return key
    for cand in _KEY_CANDIDATES:
        if records and all(cand in r for r in records):
            return cand
    return ""


def content_signature(record: Mapping[str, Any]) -> str:
    """sha256 of the record without its volatile fields, recursively."""
    def strip(v: Any) -> Any:
        if isinstance(v, Mapping):
            return {k: strip(x) for k, x in sorted(v.items())
                    if not _VOLATILE_KEY_RE.search(str(k))}
        if isinstance(v, list):
            return [strip(x) for x in v]
        return v
    body = json.dumps(strip(record), sort_keys=True, default=str,
                      ensure_ascii=False)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()[:16]


def _truthy(v: Any) -> bool:
    if isinstance(v, str):
        return v.strip().lower() in ("1", "true", "yes", "ok", "pass",
                                     "passed", "success")
    return bool(v)


def _number(v: Any) -> Optional[float]:
    if isinstance(v, bool):
        return float(v)
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, (list, tuple, dict)):
        return float(len(v))
    try:
        return float(str(v))
    except (TypeError, ValueError):
        return None


def _index(records: Sequence[dict], key: str) -> tuple[dict, list]:
    seen: dict = {}
    dups: list = []
    for r in records:
        k = str(r.get(key))
        if k in seen:
            dups.append(k)
        seen[k] = r
    return seen, dups


# ---------------------------------------------------------------------------
# the axis guard
# ---------------------------------------------------------------------------

def _is_switch(value: str) -> bool:
    """A behaviour switch, as opposed to a data path or nothing. Numbers
    are switches (thresholds, lengths, bands)."""
    v = str(value)
    if v in _OFF_VALUES:
        return False
    return "/" not in v


def _axis_prefixes(axis: Iterable[str]) -> list[str]:
    out = []
    for name in axis:
        if "_" in name:
            out.append(name.split("_", 1)[0] + "_")
    return out


def axis_check(stamp_a: Mapping, stamp_b: Mapping, axis: Sequence[str], *,
               defaults: Optional[Mapping[str, str]] = None,
               allow_shared_switches: bool = False) -> dict:
    """The instrument and environment comparison of two stamps.

    Returns ``{"state": <state or "">, "reasons": [...], "env_diff":
    {...}, "shared_switches": {...}}``; an empty state means the pair is
    one instrument state apart from exactly the declared axis."""
    defaults = dict(defaults or {})
    env_axis = [a for a in axis if a not in ("code", "platform")]
    reasons: list[str] = []
    state = ""

    def worst(new: str) -> None:
        nonlocal state
        if not state or STATES.index(new) < STATES.index(state):
            state = new

    if "platform" not in axis:
        for field in ("platform", "machine", "python"):
            if stamp_a.get(field) != stamp_b.get(field):
                worst("skewed_pair")
                reasons.append(f"{field} differs: {stamp_a.get(field)!r} vs "
                               f"{stamp_b.get(field)!r} and is not the axis")
    code_a = (stamp_a.get("code") or {}).get("hash")
    code_b = (stamp_b.get("code") or {}).get("hash")
    git_a = (stamp_a.get("git") or {})
    git_b = (stamp_b.get("git") or {})
    code_differs = (code_a != code_b
                    or git_a.get("head") != git_b.get("head")
                    or git_a.get("dirty_hash") != git_b.get("dirty_hash"))
    if "code" in axis:
        if not code_differs:
            worst("null_test_suspected")
            reasons.append("code is the declared axis but both arms ran "
                           "the same code (hash, head and diff equal)")
    elif code_differs:
        worst("skewed_pair")
        reasons.append("the code differs between the arms (hash "
                       f"{code_a} vs {code_b}, head "
                       f"{str(git_a.get('head'))[:10]} vs "
                       f"{str(git_b.get('head'))[:10]}) and is not the axis")
    for st, arm in ((stamp_a, "A"), (stamp_b, "B")):
        missing = (st.get("code") or {}).get("missing") or []
        if missing:
            worst("skewed_pair")
            reasons.append(f"arm {arm}'s code paths do not exist: {missing}")
        stale = st.get("code_changed_after_start") or []
        if stale:
            worst("skewed_pair")
            reasons.append(f"arm {arm}'s runner started before {len(stale)} "
                           f"of its code files changed (e.g. {stale[0]}); its "
                           "results are not the stamped code's")

    env_a, env_b = stamp_a.get("env"), stamp_b.get("env")
    if env_a is None or env_b is None:
        worst("axis_unknown")
        reasons.append("an arm's environment could not be read "
                       f"(A: {stamp_a.get('env_source')}, B: "
                       f"{stamp_b.get('env_source')}); check the axis by hand")
        return {"state": state, "reasons": reasons, "env_diff": {},
                "shared_switches": {}}
    pre = set(stamp_a.get("prefixes") or ()) | set(stamp_b.get("prefixes") or ())
    pre |= set(_axis_prefixes(env_axis))
    pre_t = tuple(p for p in pre if p)

    def relevant(name: str) -> bool:
        return (name not in HARNESS_VARS
                and (name in DETERMINISM_VARS or name.startswith(pre_t)
                     or name in env_axis))

    names = {n for n in set(env_a) | set(env_b) if relevant(n)}

    def resolved(env: Mapping, name: str) -> Optional[str]:
        v = env.get(name)
        if v is None or v == "":
            return defaults.get(name)
        return str(v)

    diff = {}
    for n in sorted(names | set(env_axis)):
        va, vb = resolved(env_a, n), resolved(env_b, n)
        if va != vb:
            diff[n] = [va, vb]
    undeclared = [n for n in diff if n not in env_axis]
    if undeclared:
        worst("skewed_pair")
        reasons.append("environment differs outside the declared axis: "
                       + ", ".join(f"{n}={diff[n][0]!r}->{diff[n][1]!r}"
                                   for n in undeclared[:8]))
    for n in env_axis:
        va, vb = resolved(env_a, n), resolved(env_b, n)
        if va == vb:
            worst("null_test_suspected")
            reasons.append(f"axis {n} has the same effective value in both "
                           f"arms ({va!r}): a switch at one value is no A/B")
        elif va is None or vb is None:
            worst("axis_unknown")
            reasons.append(f"axis {n} is unset in arm "
                           f"{'A' if va is None else 'B'} and has no declared "
                           "default: null -> value is a new measurement, not "
                           "an improvement (declare --default "
                           f"{n}=<the code's default>)")
    shared = {n: env_a[n] for n in sorted(names)
              if n not in env_axis and n not in diff
              and n not in DETERMINISM_VARS
              and n in env_a and _is_switch(env_a[n])}
    if shared and not allow_shared_switches:
        worst("undeclared_switch")
        reasons.append("behaviour switches set in both arms and not declared: "
                       + ", ".join(f"{k}={v!r}" for k, v in
                                   list(shared.items())[:8])
                       + f" (declare them, unset them, or set {ALLOW_SHARED_ENV}=1)")
    return {"state": state, "reasons": reasons, "env_diff": diff,
            "shared_switches": shared}


# ---------------------------------------------------------------------------
# the comparison
# ---------------------------------------------------------------------------

def _sample(keys: list, rows: dict, changed: set, k: int) -> list:
    """Every changed item (up to *k*), then unchanged ones round-robin
    over the (ok_a, ok_b) classes."""
    out = [x for x in keys if x in changed][:k]
    classes: dict = {}
    for x in keys:
        if x not in changed:
            classes.setdefault(rows[x], []).append(x)
    pools = [list(v) for _, v in sorted(classes.items(), key=lambda t: str(t[0]))]
    while len(out) < k and any(pools):
        for pool in pools:
            if pool and len(out) < k:
                out.append(pool.pop(0))
    return out


def compare(records_a: Sequence[dict], records_b: Sequence[dict], *,
            stamp_a: Optional[Mapping] = None,
            stamp_b: Optional[Mapping] = None,
            axis: Sequence[str] = (), key: str = "", ok: str = "",
            defects: Sequence[str] = (), quality: str = "",
            lower_is_better: bool = False, min_n: int = 100,
            probe: int = 24, sample: int = 100,
            defaults: Optional[Mapping[str, str]] = None,
            allow_shared_switches: Optional[bool] = None,
            rerun_a: Optional[Sequence[dict]] = None,
            prereg: str = "", env: Optional[Mapping[str, str]] = None) -> dict:
    """The verdict on two arms. See the module docstring."""
    env = os.environ if env is None else env
    if allow_shared_switches is None:
        allow_shared_switches = str(env.get(ALLOW_SHARED_ENV, "")).strip() in (
            "1", "true", "yes")
    out: dict[str, Any] = {"axis": list(axis), "reasons": [], "state": ""}

    def worst(new: str, why: str) -> None:
        cur = out["state"]
        if not cur or STATES.index(new) < STATES.index(cur):
            out["state"] = new
        out["reasons"].append(why)

    k = pick_key(list(records_a) + list(records_b), key)
    if not k:
        worst("skewed_pair", "no item key: name it with --key (a field "
              "every record of both arms carries)")
        return out
    out["key"] = k
    idx_a, dup_a = _index(records_a, k)
    idx_b, dup_b = _index(records_b, k)
    for dups, arm in ((dup_a, "A"), (dup_b, "B")):
        if dups:
            worst("skewed_pair", f"arm {arm} lists {len(dups)} item(s) more "
                  f"than once (e.g. {dups[0]!r}): the denominator moves")
    only_a = sorted(set(idx_a) - set(idx_b))
    only_b = sorted(set(idx_b) - set(idx_a))
    if only_a or only_b:
        worst("skewed_pair", f"the arms are different populations: "
              f"{len(only_a)} item(s) only in A {only_a[:5]}, {len(only_b)} "
              f"only in B {only_b[:5]}")
    if not axis:
        worst("axis_unknown", "no axis declared (--axis VAR, code or "
              "platform): what differs between the arms is not known")
    if stamp_a is None or stamp_b is None:
        worst("axis_unknown", "an arm has no provenance stamp "
              "(delfin-agent provenance): the instrument state is unchecked")
    else:
        ax = axis_check(stamp_a, stamp_b, axis, defaults=defaults,
                        allow_shared_switches=bool(allow_shared_switches))
        out["env_diff"] = ax["env_diff"]
        out["stamps_sealed"] = bool(stamp_a.get("sealed")) and bool(
            stamp_b.get("sealed"))
        if ax["state"]:
            for why in ax["reasons"]:
                worst(ax["state"], why)
        if prereg:
            p = Path(prereg)
            first = min(float(stamp_a.get("taken_at") or 0),
                        float(stamp_b.get("taken_at") or 0))
            if not p.is_file():
                out["prereg"] = "missing"
            elif p.stat().st_mtime > first:
                out["prereg"] = "written_after_the_first_arm_was_stamped"
            else:
                out["prereg"] = "before_measurement"
    if out["state"] == "skewed_pair":
        out["comparison"] = "none: not one population under one instrument"
        return out

    keys = [str(r.get(k)) for r in records_a if str(r.get(k)) in idx_b]
    keys = list(dict.fromkeys(keys))
    n = len(keys)
    sig_a = {x: content_signature(idx_a[x]) for x in keys}
    sig_b = {x: content_signature(idx_b[x]) for x in keys}
    changed = {x for x in keys if sig_a[x] != sig_b[x]}
    head = keys[:max(0, int(probe))]
    if n and not changed:
        worst("null_test_suspected", f"none of the {n} items differs between "
              "the arms: the axis did not reach the results")
    elif head and all(sig_a[x] == sig_b[x] for x in head):
        worst("null_test_suspected", f"the first {len(head)} items are "
              "identical in both arms: check that the axis reaches the code "
              "path (a counter that the new path fired) before reading more")

    never_worse: list[dict] = []
    totals: dict[str, Any] = {"n": n, "changed": len(changed)}
    cls: dict = {}
    if ok:
        ok_a = {x: _truthy(idx_a[x].get(ok)) for x in keys}
        ok_b = {x: _truthy(idx_b[x].get(ok)) for x in keys}
        totals["ok_a"] = sum(ok_a.values())
        totals["ok_b"] = sum(ok_b.values())
        for x in keys:
            cls[x] = (ok_a[x], ok_b[x])
            if ok_a[x] and not ok_b[x]:
                never_worse.append({"item": x, "term": f"{ok}: ok -> not ok"})
    else:
        cls = {x: (None, None) for x in keys}
    for field in defects:
        ta = tb = 0.0
        for x in keys:
            va = _number(idx_a[x].get(field)) or 0.0
            vb = _number(idx_b[x].get(field)) or 0.0
            ta += va
            tb += vb
            if vb > va:
                never_worse.append({"item": x, "term":
                                    f"{field}: {va:g} -> {vb:g}"})
        totals[f"{field}_a"] = ta
        totals[f"{field}_b"] = tb
        if tb > ta:
            out["reasons"].append(f"absolute {field} rose: {ta:g} -> {tb:g}")
    if quality:
        for x in keys:
            qa = _number(idx_a[x].get(quality))
            qb = _number(idx_b[x].get(quality))
            if qa is None or qb is None:
                if qa is not None and qb is None:
                    never_worse.append({"item": x, "term":
                                        f"{quality}: value -> none"})
                continue
            worse = qb > qa if lower_is_better else qb < qa
            if worse:
                never_worse.append({"item": x, "term":
                                    f"{quality}: {qa:g} -> {qb:g}"})
    if rerun_a is not None:
        idx_r, _ = _index(rerun_a, k)
        flips = [x for x in keys if x in idx_r
                 and content_signature(idx_r[x]) != sig_a[x]]
        totals["nondeterministic_a"] = len(flips)
        for x in flips:
            never_worse.append({"item": x, "term":
                                "arm A gave different content on a rerun"})
    out["totals"] = totals
    out["never_worse"] = never_worse[:_MAX_LISTED]
    out["never_worse_count"] = len(never_worse)
    if never_worse or any(r.startswith("absolute ") for r in out["reasons"]):
        worst("blocked", f"{len(never_worse)} never-worse case(s): explain "
              "each one; an instrument error is fixed in the instrument and "
              "both arms are read again, never by loosening a threshold")
    if n < int(min_n):
        worst("insufficient_n", f"{n} paired items, fewer than the minimum "
              f"{min_n}: nothing is concluded automatically, the user decides")
    if not out["state"]:
        if ok and totals.get("ok_a") == totals.get("ok_b") and not changed:
            out["state"] = "no_difference"
        else:
            out["state"] = "needs_human_check"
    out["sample"] = _sample(keys, cls, changed, int(sample))
    return out


def format_report(v: Mapping) -> str:
    """The verdict as text for a terminal or a tool result."""
    lines = [f"state: {v.get('state')}"]
    for why in v.get("reasons") or ():
        lines.append(f"- {why}")
    if v.get("comparison"):
        lines.append(f"comparison: {v['comparison']}")
    t = v.get("totals") or {}
    if t:
        lines.append("totals (absolute, paired items): " + ", ".join(
            f"{k}={t[k]:g}" if isinstance(t[k], float) else f"{k}={t[k]}"
            for k in t))
    if v.get("env_diff"):
        lines.append("effective env diff: " + ", ".join(
            f"{k}: {a!r} -> {b!r}" for k, (a, b) in v["env_diff"].items()))
    if v.get("prereg"):
        lines.append(f"pre-registration: {v['prereg']}")
    if v.get("stamps_sealed") is False:
        lines.append("stamps: read from editable files, not from the saved "
                     "store (provenance --save NAME)")
    nw = v.get("never_worse") or []
    if nw:
        lines.append(f"never-worse cases ({v.get('never_worse_count')}):")
        lines += [f"  {c['item']}: {c['term']}" for c in nw]
    if v.get("sample"):
        lines.append(f"look at these items in both arms ({len(v['sample'])}): "
                     + ", ".join(map(str, v["sample"][:30]))
                     + (" ..." if len(v["sample"]) > 30 else ""))
    if v.get("state") in ("needs_human_check", "no_difference"):
        lines.append("Not a landing: a person looks at the sample, and only "
                     "the user's own message approves.")
    return "\n".join(lines)


def exit_code(v: Mapping) -> int:
    """0 when a comparison stands (a person should look), 3 otherwise."""
    return 0 if v.get("state") in ("needs_human_check", "no_difference") else 3


__all__ = ["SCHEMA", "STATES", "ALLOW_SHARED_ENV", "provenance", "compare",
           "save_stamp", "load_stamp", "stamp_dir",
           "axis_check", "command_env", "behaviour_env", "code_hash",
           "content_signature", "load_records", "format_report", "exit_code",
           "checkout_prefix", "process_start_epoch"]
