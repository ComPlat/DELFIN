"""Loop definitions: a process the agent repeats, written down as a file.

Input: ``<workspace>/.delfin/loops/<name>.md`` (project) or
``~/.delfin/loops/<name>.md`` (personal; the project file wins on a name
collision). YAML front matter, then the steps as a numbered list::

    ---
    name: cluster-table            # optional, defaults to the file stem
    goal: Every finished job has its energy in results/energies.csv
    trigger: job                   # manual | idle | job | every 30m
    per_iteration: all             # all steps each iteration | one
    check: python scripts/check_table.py   # optional; exit 0 = goal met
    metric: 'missing=(\\d+)'       # optional; a number in the check output
    better: lower                  # lower | higher
    budget: {iterations: 10, wall: 4h, tokens: 2000000}
    stop: {plateau: 3, repeat_failures: 2}
    touch: [results/, scripts/]    # what an iteration may change ([] = all)
    measure: [tests/test_table.py] # files an iteration may not change
    branch: loop/cluster-table     # never the default branch
    after_user: pause              # pause | resume, once the user's turn ended
    ---
    1. Read the output folder of the job the wake names.
       run: python scripts/parse_energy.py <folder>
    2. Update its row in results/energies.csv.
       check: python scripts/check_table.py

Output: :class:`LoopSpec`, whose ``errors`` lists everything wrong with
the file (an empty list is a valid loop), plus :func:`discover`,
:func:`confirm` / :func:`is_confirmed`, :func:`write_template` and
:func:`ad_hoc`.

Semantics:

* The measure is the declared ``measure`` files, the loop file itself and
  every existing file the ``check`` command names. auto_guard refuses an
  iteration's write to any of them and the end-of-iteration diff reverts
  one that got through another way.
* A loop runs only when the user confirmed exactly this content. The
  ledger (``~/.delfin/loops_confirmed.json``) maps the file's resolved path
  to the sha256 of its bytes, so an edit after the confirmation -- by the
  agent writing its own loop file, or by a pull bringing a new one --
  withdraws it. Rejected: a ``confirmed: true`` key in the file, which the
  file's author (the agent, or a cloned repository) could write itself;
  and the workspace-trust gate of project skills, which trusts a directory
  rather than one loop's content.
* Budgets have defaults and ceilings; a loop cannot ask for an unbounded
  run. A command-line loop (``/loop 30m <prompt>``) is an ad-hoc spec that
  needs no confirmation, the user having typed it.

Nothing here runs a command, reads git or talks to a model.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

__all__ = [
    "LoopSpec", "Step", "discover", "load", "confirm", "is_confirmed",
    "ledger_path", "write_template", "ad_hoc", "parse",
    "DEFAULT_ITERATIONS", "DEFAULT_WALL_S", "DEFAULT_TOKENS",
]

DEFAULT_ITERATIONS = 10
MAX_ITERATIONS = 500
DEFAULT_WALL_S = 4 * 3600
MAX_WALL_S = 7 * 86400
DEFAULT_TOKENS = 2_000_000
MAX_TOKENS = 50_000_000
DEFAULT_PLATEAU = 3
DEFAULT_REPEAT_FAILURES = 2
MAX_STEPS = 30
TRIGGERS = ("manual", "idle", "job", "every")
#: Branch names a loop may never commit on, whatever the remote says.
FORBIDDEN_BRANCHES = frozenset({"main", "master", "trunk", "HEAD"})

_NAME_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{0,47}$")
_STEP_RE = re.compile(r"^\s{0,3}(\d{1,2})[.)]\s+(\S.*)$")
_STEP_ATTR_RE = re.compile(r"^\s+(run|check):\s*(\S.*)$")
_DURATION_RE = re.compile(r"^\s*(\d+)\s*([smhd])\s*$")


@dataclass
class Step:
    text: str
    run: str = ""
    check: str = ""


@dataclass
class LoopSpec:
    name: str
    goal: str = ""
    steps: list[Step] = field(default_factory=list)
    trigger: str = "manual"
    every_seconds: int = 0
    per_iteration: str = "all"
    check: str = ""
    metric: str = ""
    better: str = "lower"
    max_iterations: int = DEFAULT_ITERATIONS
    wall_seconds: int = DEFAULT_WALL_S
    max_tokens: int = DEFAULT_TOKENS
    plateau: int = DEFAULT_PLATEAU
    repeat_failures: int = DEFAULT_REPEAT_FAILURES
    touch: list[str] = field(default_factory=list)
    measure: list[str] = field(default_factory=list)
    branch: str = ""
    after_user: str = "pause"
    path: Optional[Path] = None
    scope: str = "project"          # project | personal | ad-hoc
    sha256: str = ""
    errors: list[str] = field(default_factory=list)

    def metric_value(self, output: str) -> Optional[float]:
        """The number the metric regex finds in *output* (its last match),
        or None when there is no metric or no match."""
        if not self.metric:
            return None
        try:
            hits = re.findall(self.metric, output or "")
        except re.error:
            return None
        if not hits:
            return None
        last = hits[-1]
        if isinstance(last, tuple):
            last = next((x for x in last if x), "")
        try:
            return float(last)
        except (TypeError, ValueError):
            return None

    def summary(self) -> str:
        trig = (f"every {self.every_seconds // 60}m" if self.trigger == "every"
                else self.trigger)
        return (f"{self.name} — {self.goal[:70]} · trigger {trig} · "
                f"{self.max_iterations} iterations max")


# ---------------------------------------------------------------------------
# Where loop files live
# ---------------------------------------------------------------------------

def project_dir(workspace: Path | str) -> Path:
    return Path(workspace) / ".delfin" / "loops"


def personal_dir() -> Path:
    return Path.home() / ".delfin" / "loops"


def ledger_path() -> Path:
    return Path.home() / ".delfin" / "loops_confirmed.json"


def _files(workspace: Path | str | None) -> dict[str, tuple[Path, str]]:
    out: dict[str, tuple[Path, str]] = {}
    dirs = [(personal_dir(), "personal")]
    if workspace:
        dirs.append((project_dir(workspace), "project"))
    for d, scope in dirs:              # later wins: the project's
        try:
            entries = sorted(d.glob("*.md"))
        except OSError:
            continue
        for p in entries:
            if p.is_file():
                out[p.stem] = (p, scope)
    return out


def discover(workspace: Path | str | None) -> list[LoopSpec]:
    """Every loop file the session can see, parsed (invalid ones included,
    with their errors), sorted by name."""
    specs = []
    for _stem, (p, scope) in sorted(_files(workspace).items()):
        specs.append(_load_path(p, scope, workspace))
    return specs


def load(name: str, workspace: Path | str | None) -> Optional[LoopSpec]:
    """The loop called *name* (file stem or front-matter name), or None."""
    name = str(name or "").strip()
    files = _files(workspace)
    if name in files:
        p, scope = files[name]
        return _load_path(p, scope, workspace)
    for spec in discover(workspace):
        if spec.name == name:
            return spec
    return None


def _load_path(p: Path, scope: str, workspace) -> LoopSpec:
    try:
        raw = p.read_bytes()
    except OSError as exc:
        return LoopSpec(name=p.stem, path=p, scope=scope,
                        errors=[f"cannot read {p}: {exc}"])
    spec = parse(raw.decode("utf-8", errors="replace"), stem=p.stem,
                 workspace=workspace, path=p)
    spec.scope = scope
    spec.sha256 = hashlib.sha256(raw).hexdigest()
    return spec


# ---------------------------------------------------------------------------
# Parsing and validation
# ---------------------------------------------------------------------------

def _front_matter(text: str) -> tuple[dict, str, str]:
    """(meta, body, error)."""
    if not text.startswith("---"):
        return {}, text, "no front matter (the file must start with ---)"
    rest = text.split("\n", 1)[1] if "\n" in text else ""
    end = rest.find("\n---")
    if rest.startswith("---"):
        end = 0
    if end < 0:
        return {}, text, "front matter is not closed with ---"
    block, body = rest[:end], rest[end + 4:].lstrip("\n")
    try:
        import yaml
        meta = yaml.safe_load(block) or {}
    except Exception as exc:
        return {}, body, f"front matter is not valid YAML: {exc}"
    if not isinstance(meta, dict):
        return {}, body, "front matter must be a mapping"
    return meta, body, ""


def _duration(value: Any) -> Optional[int]:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return int(value)
    m = _DURATION_RE.match(str(value or ""))
    if not m:
        return None
    return int(m.group(1)) * {"s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)]


def _int(value: Any) -> Optional[int]:
    if isinstance(value, bool):
        return None
    try:
        return int(str(value).replace("_", ""))
    except (TypeError, ValueError):
        return None


def _steps(body: str) -> list[Step]:
    steps: list[Step] = []
    for line in body.splitlines():
        m = _STEP_RE.match(line)
        if m:
            steps.append(Step(text=m.group(2).strip()))
            continue
        if not steps:
            continue
        a = _STEP_ATTR_RE.match(line)
        if a:
            setattr(steps[-1], a.group(1), a.group(2).strip())
        elif line.startswith((" ", "\t")) and line.strip():
            steps[-1].text += " " + line.strip()
    return steps


def _rel_paths(value: Any, key: str, errors: list[str]) -> list[str]:
    if value in (None, ""):
        return []
    items = value if isinstance(value, list) else [value]
    out = []
    for item in items:
        s = str(item or "").strip()
        if not s:
            continue
        parts = Path(s).parts
        if Path(s).is_absolute() or ".." in parts or s.startswith("~"):
            errors.append(f"{key}: {s!r} must be a path inside the workspace")
            continue
        out.append(s)
    return out


def _check_files(check: str, workspace) -> list[str]:
    """The relative paths the check command names: a word with a slash or
    a file suffix that is not an option. Named, not only existing: a check
    script the loop is about to create is the measure as much as one that
    exists."""
    if not check:
        return []
    try:
        words = shlex.split(check)
    except ValueError:
        words = check.split()
    out = []
    for w in words:
        if w.startswith("-") or "=" in w or ("/" not in w
                                             and not Path(w).suffix):
            continue
        p = Path(w)
        if p.is_absolute() or ".." in p.parts or w.startswith("~"):
            continue
        out.append(str(p))
    return out


def parse(text: str, *, stem: str = "", workspace=None,
          path: Optional[Path] = None) -> LoopSpec:
    """One loop file's text as a LoopSpec (see the module docstring)."""
    meta, body, err = _front_matter(text)
    errors: list[str] = [err] if err else []
    name = str(meta.get("name") or stem or "").strip()
    spec = LoopSpec(name=name, path=path)
    if not _NAME_RE.match(name):
        errors.append(f"name: {name!r} must be lower case letters, digits, "
                      "- or _, at most 48 characters")
    spec.goal = " ".join(str(meta.get("goal") or "").split())
    if not spec.goal:
        errors.append("goal: required (one sentence: what is true when the "
                      "loop is done)")
    spec.steps = _steps(body)
    if not spec.steps:
        errors.append("steps: none found (write them as a numbered list: "
                      "'1. ...')")
    elif len(spec.steps) > MAX_STEPS:
        errors.append(f"steps: {len(spec.steps)} > {MAX_STEPS}")

    trig = str(meta.get("trigger") or "manual").strip().lower()
    if trig.startswith("every"):
        secs = _duration(trig[5:])
        if not secs or secs < 60:
            errors.append(f"trigger: {trig!r} needs an interval of at least "
                          "1m, e.g. 'every 30m'")
        else:
            spec.every_seconds = secs
        spec.trigger = "every"
    elif trig in TRIGGERS:
        spec.trigger = trig
    else:
        errors.append(f"trigger: {trig!r} is not one of manual, idle, job, "
                      "every <N>m")

    per = str(meta.get("per_iteration") or "all").strip().lower()
    if per not in ("all", "one"):
        errors.append(f"per_iteration: {per!r} is not all or one")
    spec.per_iteration = per

    spec.check = str(meta.get("check") or "").strip()
    spec.metric = str(meta.get("metric") or "").strip()
    if spec.metric:
        try:
            if re.compile(spec.metric).groups < 1:
                errors.append("metric: the regex needs one (group) around "
                              "the number")
        except re.error as exc:
            errors.append(f"metric: not a valid regex ({exc})")
        if not spec.check:
            errors.append("metric: needs a check command whose output it reads")
    spec.better = str(meta.get("better") or "lower").strip().lower()
    if spec.better not in ("lower", "higher"):
        errors.append(f"better: {spec.better!r} is not lower or higher")

    budget = meta.get("budget") or {}
    if not isinstance(budget, dict):
        errors.append("budget: must be a mapping (iterations, wall, tokens)")
        budget = {}
    it = _int(budget.get("iterations", DEFAULT_ITERATIONS))
    if it is None or not 1 <= it <= MAX_ITERATIONS:
        errors.append(f"budget.iterations: {budget.get('iterations')!r} must "
                      f"be 1..{MAX_ITERATIONS}")
    else:
        spec.max_iterations = it
    wall = _duration(budget.get("wall", DEFAULT_WALL_S))
    if wall is None or not 60 <= wall <= MAX_WALL_S:
        errors.append(f"budget.wall: {budget.get('wall')!r} must be 1m..7d")
    else:
        spec.wall_seconds = wall
    tok = _int(budget.get("tokens", DEFAULT_TOKENS))
    if tok is None or not 1000 <= tok <= MAX_TOKENS:
        errors.append(f"budget.tokens: {budget.get('tokens')!r} must be "
                      f"1000..{MAX_TOKENS}")
    else:
        spec.max_tokens = tok

    stop = meta.get("stop") or {}
    if not isinstance(stop, dict):
        errors.append("stop: must be a mapping (plateau, repeat_failures)")
        stop = {}
    pl = _int(stop.get("plateau", DEFAULT_PLATEAU))
    if pl is None or not 1 <= pl <= 50:
        errors.append("stop.plateau: must be 1..50")
    else:
        spec.plateau = pl
    rf = _int(stop.get("repeat_failures", DEFAULT_REPEAT_FAILURES))
    if rf is None or not 1 <= rf <= 20:
        errors.append("stop.repeat_failures: must be 1..20")
    else:
        spec.repeat_failures = rf

    spec.touch = _rel_paths(meta.get("touch"), "touch", errors)
    declared = _rel_paths(meta.get("measure"), "measure", errors)
    own = []
    if path is not None and workspace:
        try:
            own = [str(Path(path).resolve().relative_to(
                Path(workspace).resolve()))]
        except (ValueError, OSError):
            own = []
    spec.measure = list(dict.fromkeys(
        declared + own + _check_files(spec.check, workspace)))

    branch = str(meta.get("branch") or f"loop/{name}").strip()
    if branch in FORBIDDEN_BRANCHES or not re.match(r"^[\w./-]{1,80}$", branch) \
            or branch.startswith(("-", "/")) or ".." in branch:
        errors.append(f"branch: {branch!r} is not allowed (a loop commits on "
                      "its own branch, never the default branch)")
    spec.branch = branch

    # The harness runs these itself, after the iteration's veto is off: a
    # command that publishes would publish unattended. Classified by the
    # same rules the veto uses.
    from .auto_guard import allowed_ssh_hosts, classify_call
    ssh_hosts = allowed_ssh_hosts()
    commands = [("check", spec.check)] + [
        (f"step {i} {k}", getattr(st, k)) for i, st in enumerate(spec.steps, 1)
        for k in ("check", "run") if getattr(st, k)]
    for label, command in commands:
        if not command:
            continue
        verdict = classify_call("bash", {"command": command},
                                workspace=str(workspace or ""),
                                ssh_hosts=ssh_hosts)
        if verdict is not None and verdict.queue:
            errors.append(f"{label}: {command!r} would publish or send "
                          f"({verdict.reason}); a loop's commands stay local")

    after = str(meta.get("after_user") or "pause").strip().lower()
    if after not in ("pause", "resume"):
        errors.append(f"after_user: {after!r} is not pause or resume")
    spec.after_user = after
    spec.errors = errors
    return spec


# ---------------------------------------------------------------------------
# Confirmation
# ---------------------------------------------------------------------------

def _ledger() -> dict:
    try:
        data = json.loads(ledger_path().read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _key(spec: LoopSpec) -> str:
    try:
        return str(Path(spec.path).resolve()) if spec.path else ""
    except OSError:
        return str(spec.path)


def is_confirmed(spec: Optional[LoopSpec]) -> bool:
    """True when the user confirmed exactly this file content (an ad-hoc
    spec, typed by the user, always is)."""
    if spec is None or spec.errors:
        return False
    if spec.scope == "ad-hoc":
        return True
    rec = _ledger().get(_key(spec))
    return bool(isinstance(rec, dict) and spec.sha256
                and rec.get("sha256") == spec.sha256)


def confirm(spec: LoopSpec) -> None:
    """Record the user's confirmation of this content. Called only from a
    command the user typed (``/loop confirm``, the run dialog)."""
    if spec is None or spec.errors:
        raise ValueError("an invalid loop cannot be confirmed: "
                         + "; ".join(spec.errors if spec else ["missing"]))
    if spec.scope == "ad-hoc":
        return
    data = _ledger()
    data[_key(spec)] = {"sha256": spec.sha256, "name": spec.name,
                        "at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    path = ledger_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=1, sort_keys=True), encoding="utf-8")
    try:
        os.chmod(tmp, 0o600)
    except OSError:
        pass
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# New loops
# ---------------------------------------------------------------------------

_TEMPLATE = """---
name: {name}
goal: {goal}
trigger: manual          # manual | idle | job | every 30m
per_iteration: all       # all | one
# check: <command>       # exit 0 = goal met; its output feeds metric
# metric: 'failed=(\\d+)'
# better: lower
budget:
  iterations: {iterations}
  wall: 4h
  tokens: 2000000
stop:
  plateau: 3
  repeat_failures: 2
touch: []                # paths an iteration may change; [] = the workspace
measure: []              # files an iteration may not change
branch: loop/{name}
after_user: pause
---
# Steps
1. {first_step}
"""


def write_template(workspace: Path | str, name: str, *, goal: str = "",
                   first_step: str = "Describe the first step.",
                   iterations: int = DEFAULT_ITERATIONS) -> Path:
    """Write a new, unconfirmed loop file; FileExistsError when one exists."""
    if not _NAME_RE.match(name or ""):
        raise ValueError(f"{name!r} is not a loop name (lower case letters, "
                         "digits, - or _)")
    path = project_dir(workspace) / f"{name}.md"
    if path.exists():
        raise FileExistsError(str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_TEMPLATE.format(
        name=name, goal=" ".join((goal or "Say what is true when done.").split()),
        iterations=int(iterations), first_step=first_step), encoding="utf-8")
    return path


def ad_hoc(interval: str, prompt: str, workspace=None) -> LoopSpec:
    """A loop typed on the command line: ``/loop 30m <prompt>``."""
    secs = _duration(interval)
    spec = LoopSpec(name="adhoc", goal=" ".join(str(prompt or "").split()),
                    steps=[Step(text=" ".join(str(prompt or "").split()))],
                    trigger="every", every_seconds=int(secs or 0),
                    scope="ad-hoc", branch="loop/adhoc")
    if not secs or secs < 60:
        spec.errors.append(f"{interval!r} is not an interval of at least 1m "
                           "(5m, 2h, 1d)")
    if not spec.goal:
        spec.errors.append("the prompt is empty")
    return spec
