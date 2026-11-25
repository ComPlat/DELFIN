"""A member's real runs: its change used the way a person would use it.

A member cannot make a real run itself: its cage strips every provider
key from the environment, and that is right. So it writes what it wants
run into ``scratch/realrun.json``, and the stage (``swarm/stages.py``)
runs it here, outside the cage, twice -- on the member's base and on its
branch HEAD -- with the same prompt and model, each in a fresh HOME and
workspace that are removed afterwards, a stop included. A test passing is
not a run: until 2026-11-22 members reported green suites for changes
whose behaviour nobody had seen.

``scratch/realrun.json``::

    {"probes": [{"name": "inline-xyz", "prompt": "...", "workspace": "empty" | "repo",
                 "turns": 20}],
     "provider": "kit", "model": "..."}      # default: the member's own

At most three probes; a name is ``[a-z0-9-]{1,40}``; a prompt at most 2000
characters; at most 40 turns. Output: ``scratch/realrun/<probe>-<side>.jsonl``
(redacted: secrets, and home paths as ``<home>``) and
``scratch/realrun/summary.md``, a table of stop reason, tool calls,
tokens, seconds and error per probe and side, then each final answer.

Each run is ``python -P -m delfin.agent.cli chat`` with ``PYTHONPATH`` on
the side's code copy, checked to import from there: ``-P`` because a
``repo`` workspace is itself a copy of the repository, which ``python
-m`` would otherwise put first on the path (all six runs of one member
were refused by that check, 2026-11-24). The key is read with DELFIN's
own loader in this process and given to the child in its environment
only; DELFIN copies it into the child's HOME, which is why that HOME is
always removed.

Did the run touch the change at all
----------------------------------
A run of the branch that never executes the changed code compares two
copies of the same behaviour and measures nothing (three probes of one
package on 2026-11-25 never reached it, and only the transcripts showed
that). So each probe and side also records:

* ``offered`` -- the tool list the model was GIVEN in its first turn: the
  ``tools`` keyword of the first request that carries one. The action
  classifier's own requests (``_classifier_once``) carry none, so they
  cannot be mistaken for a turn;
* ``called`` -- the tools the model CALLED, from the result line's
  ``tool_calls``;
* ``reached`` -- ``yes`` / ``no`` / ``unknown``: whether a function the
  branch changed was executed. The changed functions' qualified names come
  from the branch's diff against its base (``_changed_names``: the diff's
  hunks intersected with the base and branch syntax trees); execution is
  observed by ``_PROBE`` below, which the child loads through a
  ``sitecustomize`` on its ``PYTHONPATH`` (measured: ``sitecustomize``
  still runs under ``python -P``, whose ``sys.path[0]`` is the harness
  directory, not the workspace).

``reached`` is ``yes`` when a changed function's frame ran, ``no`` when the
trace ran to the end and none did, and ``unknown`` when the trace cannot
say -- the child died before writing it, the diff changed no Python
function, or the run never imported ``delfin`` at all. ``reached_reason``
says which, because a bare ``unknown`` is not evidence of anything.

Execution is watched with ``sys.setprofile``, whose callback is entered on
every Python call in the child. Measured (3.12.3,
``scratch/r5-bench-probe-overhead.py``): **358 ns per call**, **+22 %** on an
82 ms import-bound child. The path filter is NOT a speed measure -- with it
removed the overhead was the same within noise (0.99x) because entering the
callback is the cost, not deciding inside it. The filter is kept for what it
does change: the record holds DELFIN's own frames (16 of 575 on that
workload), so the trace stays small and a changed name is matched against the
code under test rather than against the whole standard library. It is not
free; a probe's 1500 s budget absorbs it, and a member whose change is worth
measuring is not timed out by 358 ns a call.

The record is written as the run goes and again when it ends, so a probe that
is killed still leaves what it saw; a probe that leaves nothing degrades to
``unknown`` rather than failing the member's whole stage.
"""
from __future__ import annotations

import ast
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from . import wave

MAX_PROBES = 3
_NAME = re.compile(r"^[a-z0-9-]{1,40}$")
#: The key each provider's client reads (``api_client.create_client``,
#: pinned by tests/test_a_wave_is_run_through_delfin.py).
PROVIDER_KEYS = {"kit": "KIT_TOOLBOX_API_KEY", "blablador": "BLABLADOR_API_KEY",
                 "openai": "OPENAI_API_KEY", "anthropic": "ANTHROPIC_API_KEY"}


def read_request(worktree: Path) -> dict:
    """The member's ``scratch/realrun.json``, checked; WaveError names what
    is wrong."""
    try:
        req = json.loads((Path(worktree) / "scratch" / "realrun.json").read_text())
    except (OSError, ValueError) as exc:
        raise wave.WaveError(f"scratch/realrun.json is missing or not JSON ({exc})")
    probes = req.get("probes") if isinstance(req, dict) else None
    if not isinstance(probes, list) or not 1 <= len(probes) <= MAX_PROBES:
        raise wave.WaveError(f"realrun.json needs 1 to {MAX_PROBES} probes")
    for p in probes:
        if not isinstance(p, dict) or not _NAME.match(str(p.get("name") or "")):
            raise wave.WaveError("a probe's name is [a-z0-9-]{1,40}")
        if not str(p.get("prompt") or "").strip() or len(p["prompt"]) > 2000:
            raise wave.WaveError(f"probe {p['name']}: the prompt is empty or over 2000 characters")
        if p.get("workspace", "") not in ("", "empty", "repo"):
            raise wave.WaveError(f"probe {p['name']}: workspace is 'empty' or 'repo'")
        p["turns"] = min(int(p.get("turns") or 30), 40)
    return req


def _redact(text: str) -> str:
    try:
        from delfin.agent.bug_report import redact_for_report
        text = redact_for_report(text)[0]
    except Exception:
        text = re.sub(r"[A-Za-z0-9_\-]{32,}", "<long>", text)
    return re.sub(r"/home/[^/\s\"'`]*", "<home>", text)


def _archive(repo: Path, ref: str, dest: Path) -> None:
    dest.mkdir(parents=True)
    from delfin.agent import safe_git
    r = safe_git.run(["archive", ref], str(repo), text=False, timeout=300)
    if r.returncode != 0:
        raise wave.WaveError(f"git archive {ref} failed: {r.stderr[-300:]!r}")
    tar = r.stdout
    subprocess.run(["tar", "-x", "-C", str(dest)], input=tar, check=True)


# -- did the run touch the change --------------------------------------------------

_HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", re.M)


def _changed_lines(diff: str) -> dict:
    """``{path: [(first_line, count), ...]}`` of the ADDED side of a diff.

    Reads a ``--unified=0`` diff, so every hunk is a run of lines the branch
    touched. ``+++ b/<path>`` names the file (``/dev/null``: it is gone in the
    branch, and the caller falls back to the base's source); a rename is told
    apart from a real edit by its own header lines.
    """
    out: dict[str, list] = {}
    path = None
    for line in diff.splitlines():
        if line.startswith("+++ "):
            target = line[4:].strip()
            if target == "/dev/null":
                path = None
            elif target.startswith("b/"):
                path = target[2:]
            else:
                path = target
        elif line.startswith("@@ ") and path:
            m = _HUNK.match(line)
            if m:
                out.setdefault(path, []).append((int(m.group(1)), int(m.group(2) or 1)))
    return out


def _def_ranges(src: str) -> dict:
    """``{qualified_name: (first_line, last_line)}`` for one Python source.

    The qualified name is built the way CPython builds ``co_qualname``, because
    that is what the child's trace reports: a class member joins with ``.``, a
    function nested in a function joins with ``.<locals>.``. Measured on
    CPython 3.12 (a nested ``leaf`` in ``Outer.Inner.deep`` reports
    ``Outer.Inner.deep.<locals>.leaf``); a rule that joined plain functions with
    ``.`` would silently never match anything and read as "the change was not
    reached". ``<module>`` is added for a file whose diff touches module level,
    which is the frame that runs when the module is imported.
    """
    out: dict[str, tuple] = {}
    try:
        tree = ast.parse(src)
    except (SyntaxError, ValueError):
        return out

    def walk(node, prefix, in_fn):
        for n in ast.iter_child_nodes(node):
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                join = ".<locals>." if in_fn else "."
                qual = n.name if prefix == "" else prefix + join + n.name
                try:
                    out[qual] = (n.lineno, n.end_lineno or n.lineno)
                except AttributeError:
                    out[qual] = (n.lineno, n.lineno)
                walk(n, qual, not isinstance(n, ast.ClassDef))
            else:
                walk(n, prefix, in_fn)

    walk(tree, "", False)
    return out


def changed_functions(repo, base: str, head: str) -> dict:
    """What the branch changed that Python can execute, and where.

    ``{"names": {"pkg.mod.Class.method": "pkg/mod.py", ...}, "files": [...],
    "note": str}``. A name is ``<module>``-suffixed when the change sits at
    module level: ``{"pkg.mod.<module>": "pkg/mod.py"}``. Files the branch
    changed that hold no Python, and Python files whose diff touches no
    statement inside a def/class, contribute nothing -- and the ``note`` says
    so, because "the diff changed no function" must not read as "the probe
    missed the change".
    """
    from delfin.agent import safe_git
    r = safe_git.run(["diff", "--unified=0", "--find-renames", f"{base}...{head}", "--", "*.py"],
                     str(repo), timeout=300)
    if r.returncode != 0:
        raise wave.WaveError(f"git diff {base}...{head} failed: {(r.stderr or '')[-300:]!r}")
    hunks = _changed_lines(r.stdout or "")
    names: dict[str, str] = {}
    for path, spans in sorted(hunks.items()):
        if not path.endswith(".py"):
            continue
        src = safe_git.run(["show", f"{head}:{path}"], str(repo), timeout=120)
        if src.returncode != 0:                      # deleted in the branch
            src = safe_git.run(["show", f"{base}:{path}"], str(repo), timeout=120)
        if src.returncode != 0:
            continue
        ranges = _def_ranges(src.stdout or "")
        module_only = False
        for first, count in spans:
            last = first + max(count - 1, 0)
            # A hunk inside a nested def names that def and every def enclosing
            # it: running the outer one does not run the inner one, so each is
            # recorded and the trace decides which actually ran.
            hit = [q for q, (a, b) in ranges.items() if a <= last and b >= first]
            for q in hit:
                names[_qualify(path, q)] = path
            if not hit:
                module_only = True
        if module_only:
            names[_qualify(path, "<module>")] = path
    note = ""
    if not names:
        note = ("the diff touches no Python function" if hunks
                else "the branch changed no .py file")
    return {"names": names, "files": sorted({p for p in hunks}), "note": note}


def _qualify(path: str, qual: str) -> str:
    """``pkg/mod.py`` + ``C.m`` -> ``pkg.mod.C.m`` (the name the trace reports)."""
    mod = path[:-3] if path.endswith(".py") else path
    parts = [p for p in mod.split("/") if p]
    if parts and parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts + [qual])


#: The child's own observer, written to ``sitecustomize.py`` next to the run
#: and put on the child's ``PYTHONPATH``. ``sitecustomize`` is imported by
#: ``site`` before anything else, which is the only hook that sees the frames
#: of ``python -m delfin.agent.cli`` itself (measured: under ``-P``
#: ``sys.path[0]`` is the harness directory, so the workspace cannot shadow
#: it, and the module still runs). It records two things into the run's own
#: temporary directory -- not the workspace, so a probe cannot edit its own
#: evidence, and not ``~/.delfin``, which is nobody's to write into but its
#: owner's:
#:
#: ``offered``  the tool names of the FIRST request that carries a tool list.
#:              DELFIN builds that list inside ``stream_message`` and hands it
#:              to the provider client as ``kwargs["tools"]``, so the request
#:              frame is the only place the whole surface exists; the action
#:              classifier's own requests carry none, so a turn is not
#:              confused with a classification. Both wire shapes are watched
#:              (``chat.completions.create`` and ``messages.stream``).
#: ``funcs``    every ``delfin`` function frame that ran, as the qualified
#:              name ``changed_functions`` produces.
#:
#: ``threading.setprofile`` is not decoration: DELFIN issues its request from
#: a thread it starts for that purpose (``_spawn_create``), and a profile set
#: with ``sys.setprofile`` alone never sees that thread (measured on 3.12 --
#: without it every probe would answer "no tools were offered").
_PROBE = r'''
"""DELFIN's real-run observer (installed by delfin/agent/swarm/realrun.py)."""
import atexit
import json
import os
import sys
import threading

OUT = os.environ.get("DELFIN_PROBE_OUT") or ""
ROOT = os.environ.get("DELFIN_PROBE_ROOT") or ""
_NAMES = os.environ.get("DELFIN_PROBE_NAMES") or ""
#: The request functions that carry the tool list. Matched on the code object's
#: name only, so the expensive part (reading ``f_locals``) happens for a
#: handful of frames, never for every call.
_WIRE = frozenset(("create", "stream"))

state = {"offered": None, "offered_from": "", "funcs": [], "seen": set(),
         "wire_seen": 0, "delfin_seen": 0, "broken": ""}


def _names():
    try:
        with open(_NAMES) as fh:
            return frozenset(json.load(fh))
    except Exception:
        return frozenset()


CHANGED = _names()


def _write():
    if not OUT:
        return
    rec = {"offered": state["offered"], "offered_from": state["offered_from"],
           "funcs": state["funcs"], "wire_seen": state["wire_seen"],
           "delfin_seen": state["delfin_seen"], "broken": state["broken"],
           "changed_seen": sorted(set(state["funcs"]) & CHANGED)}
    try:
        tmp = OUT + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(rec, fh)
        os.replace(tmp, OUT)
    except Exception:
        pass


def _module_of(filename):
    """``<root>/delfin/agent/x.py`` -> ``delfin.agent.x`` (root stripped)."""
    try:
        rel = os.path.relpath(filename, ROOT)
    except Exception:
        return ""
    if rel.startswith(".."):
        return ""
    if not rel.endswith(".py"):
        return ""
    parts = rel[:-3].split(os.sep)
    if parts and parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _profile(frame, event, arg):
    """One ``call`` event, from every thread. Runs for the child's whole life.

    Kept to a name test and, for a ``delfin`` frame, one set lookup: the cost
    is in the realrun package's PR, measured, not assumed.
    """
    if event != "call":
        return None
    try:
        code = frame.f_code
        name = code.co_name
        if name in _WIRE:
            state["wire_seen"] += 1
            if state["offered"] is None:
                tools = frame.f_locals.get("tools")
                if not tools:                     # kwargs built by the caller
                    back = frame.f_back
                    while back is not None:
                        kw = back.f_locals.get("kwargs")
                        if isinstance(kw, dict) and kw.get("tools"):
                            tools = kw["tools"]
                            break
                        back = back.f_back
                if tools:
                    got = []
                    for t in tools:
                        try:
                            got.append(t["function"]["name"])
                        except Exception:
                            try:
                                got.append(t.name)
                            except Exception:
                                pass
                    if got:
                        state["offered"] = got
                        state["offered_from"] = name
                        _write()
            return None
        fn = code.co_filename
        if not fn.startswith(ROOT):               # not DELFIN's own code
            return None
        state["delfin_seen"] += 1
        qual = code.co_qualname if hasattr(code, "co_qualname") else name
        if qual == "<module>":
            q = _module_of(fn) + ".<module>" if _module_of(fn) else ""
        else:
            mod = _module_of(fn)
            q = mod + "." + qual if mod else ""
        if q and q not in state["seen"]:
            state["seen"].add(q)
            state["funcs"].append(q)
            if q in CHANGED:                       # rare: write as it happens,
                _write()                           # so a killed run still tells
    except Exception as exc:                       # never break the run
        state["broken"] = type(exc).__name__
    return None


if OUT:
    sys.setprofile(_profile)
    threading.setprofile(_profile)
    atexit.register(_write)
'''


def _instrument(codes: Path, side: str, names: dict) -> tuple:
    """The observer directory for one side, and the child's extra environment.

    Returns ``(dir, {ENV: value})``. The directory holds ``sitecustomize.py``
    and the changed names; it goes on ``PYTHONPATH`` before the code copy,
    which stays the source of ``delfin`` (the existing import check still
    decides that). ``DELFIN_PROBE_ROOT`` is that copy: the observer calls a
    frame DELFIN's own only when its file lives under it, so a second checkout
    in the workspace (a ``repo`` probe) is not traced as the code under test.
    """
    d = codes / f"probe-{side}"
    d.mkdir(parents=True, exist_ok=True)
    (d / "sitecustomize.py").write_text(_PROBE)
    (d / "names.json").write_text(json.dumps(sorted(names)))
    trace = d / "trace.json"
    # Stale evidence is worse than none: the file is removed before the run
    # starts, so a child that dies before writing cannot be read as a run that
    # reached the change (its previous probe's trace would still be there).
    trace.unlink(missing_ok=True)
    return d, {"DELFIN_PROBE_OUT": str(trace),
               "DELFIN_PROBE_NAMES": str(d / "names.json"),
               # Resolved, because the child reports the path it imported from
               # and a symlinked temp root would match nothing: every frame
               # would read as "not DELFIN" and answer ``unknown``.
               "DELFIN_PROBE_ROOT": os.path.realpath(codes / side)}


def _read_trace(d: Path, names: dict, row: dict) -> dict:
    """``reached`` / ``reached_reason`` / ``offered`` for one probe and side."""
    path = d / "trace.json"
    if not path.is_file():
        return {"reached": "unknown", "reached_reason": "the observer wrote nothing "
                "(the child did not reach the point it records at)",
                "offered": [], "called_tools": row.get("called") or []}
    try:
        tr = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        return {"reached": "unknown", "reached_reason": f"the trace did not parse ({exc})",
                "offered": [], "called_tools": row.get("called") or []}
    offered = [str(t) for t in (tr.get("offered") or [])]
    called = row.get("called") or []
    hit = sorted(set(tr.get("funcs") or []) & set(names))
    out = {"reached": "no", "reached_reason": "", "offered": offered, "called_tools": called,
           "wire_seen": tr.get("wire_seen"), "delfin_frames": tr.get("delfin_seen")}
    if hit:
        out["reached"] = "yes"
        out["reached_reason"] = "executed " + ", ".join(hit[:4]) + (
            f" (+{len(hit) - 4} more)" if len(hit) > 4 else "")
    elif not names:
        out["reached"] = "unknown"
        out["reached_reason"] = "the branch changed no Python function to reach"
    elif not tr.get("delfin_seen"):
        out["reached"] = "unknown"
        out["reached_reason"] = "the run never executed DELFIN code, so nothing of " \
                                "the change ran (no delfin frame was entered)"
    else:
        out["reached_reason"] = (f"the trace covered {len(set(tr.get('funcs') or []))} "
                                 "delfin functions and none of the changed ones")
    return out


def _summary(log: Path) -> dict:
    tools, result, called = 0, {}, []
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith('{"type": "tool_use"'):
            tools += 1
            try:
                name = json.loads(line).get("name")
            except ValueError:
                name = None
            if name and name not in called:
                called.append(name)
        elif line.startswith('{"type": "result"'):
            try:
                result = json.loads(line)
            except ValueError:
                pass
    # The result line is the child's own account of the turn and names every
    # call it made; the emitted ``tool_use`` events are the second source, and
    # they are what a stand-in child writes. Union, in first-seen order.
    for c in result.get("tool_calls") or []:
        name = c.get("name") if isinstance(c, dict) else c
        if name and name not in called:
            called.append(name)
    return {"stop": result.get("stop_reason"), "tools": tools, "called": called,
            "in": result.get("input_tokens"), "out": result.get("output_tokens"),
            "answer": (result.get("text") or "")[-1200:], "error": (result.get("error") or "")[:300]}


def _one(probe: dict, side: str, code: Path, provider: str, model: str, key: str,
         out: Path, probe_dir: Path | None = None,
         probe_env: dict | None = None) -> dict:
    with tempfile.TemporaryDirectory(prefix="delfin-realrun-") as tmp:
        tmp = Path(tmp)
        home, ws = tmp / "home", tmp / "ws"
        home.mkdir()
        if probe.get("workspace") == "repo":
            subprocess.run(["cp", "-a", f"{code}/.", str(ws)], check=True)
        else:
            ws.mkdir()
        env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "HOME": str(home),
               "LANG": "C.UTF-8", "PYTHONPATH": str(code), "DELFIN_NO_BROWSER": "1"}
        if key and provider in PROVIDER_KEYS:
            env[PROVIDER_KEYS[provider]] = key
        for name in ("XDG_RUNTIME_DIR", "DBUS_SESSION_BUS_ADDRESS"):
            if os.environ.get(name):
                env[name] = os.environ[name]
        where = subprocess.run(
            [sys.executable, "-P", "-c", "import delfin, os; print(os.path.realpath(delfin.__file__))"],
            cwd=str(ws), env=env, capture_output=True, text=True)
        if not where.stdout.strip().startswith(os.path.realpath(code)):
            return {"error": "the run imported DELFIN from elsewhere: "
                             + _redact(where.stdout.strip() or where.stderr[-200:])}
        if probe_dir is not None:
            # The observer rides along the code copy, first on the path so its
            # sitecustomize is the one Python imports. The check above ran
            # without it on purpose: a probe that never started must not leave
            # a trace behind and read as a run that reached the change. The
            # code path is resolved for the same reason the observer's root is:
            # the child reports the file names it imported under, and a path
            # through a symlink would match no root and read as ``unknown``.
            env["PYTHONPATH"] = os.pathsep.join([str(probe_dir), os.path.realpath(code)])
            env.update(probe_env or {})
        cmd = [sys.executable, "-P", "-m", "delfin.agent.cli", "chat", "-p", probe["prompt"],
               "--cwd", str(ws), "--provider", provider, "--model", model,
               "--permission-mode", "default", "--unattended", "--isolate",
               "--max-run-seconds", "1200", "--max-run-turns", str(probe["turns"]),
               "--output-format", "stream-json"]
        log = out / f"{probe['name']}-{side}.jsonl"
        t0 = time.time()
        try:
            proc = subprocess.run(wave._scoped(cmd, wave.settings()), cwd=str(ws), env=env,
                                  capture_output=True, text=True, timeout=1500)
            text = proc.stdout + ("\n" + proc.stderr[-4000:] if proc.stderr else "")
        except subprocess.TimeoutExpired as exc:
            text = (exc.stdout or "") if isinstance(exc.stdout, str) else ""
            text += "\n(timed out after 1500 s)"
        log.write_text(_redact(text))
        row = _summary(log)
        row["s"] = round(time.time() - t0)
        return row


#: The three states the summary column can hold, spelled out so a reader
#: cannot mistake ``unknown`` for ``no``.
_REACHED_LABEL = {"yes": "yes", "no": "no", "unknown": "unknown"}


def run(run_name: str, tid: str) -> dict:
    """The member's real runs, base and branch; returns the rows and writes
    ``scratch/realrun/summary.md``.

    The changed functions come from ONE diff -- the branch against its base --
    and both sides are traced against that same set. That is deliberate: a
    ``yes`` on the base and a ``yes`` on the branch means the probe walked the
    changed code on both, which is what makes a difference between the two
    rows attributable to the change; ``no`` on both means the probe measured
    something else entirely.
    """
    st = wave.load_state(run_name)
    a = st["agents"].get(tid)
    if a is None:
        raise wave.WaveError(f"{run_name}/{tid} was never started")
    wt = Path(a["worktree"])
    req = read_request(wt)
    task = a.get("task") or {}
    provider = str(req.get("provider") or task.get("provider") or "kit")
    model = str(req.get("model") or task.get("model") or "")
    if not model:
        raise wave.WaveError("no model: name one in realrun.json or in the task")
    key = ""
    if provider in PROVIDER_KEYS:
        from delfin.agent.credentials import load_credential
        key = load_credential(PROVIDER_KEYS[provider]) or ""
    repo = Path(st.get("repo") or ".")
    out = wt / "scratch" / "realrun"
    out.mkdir(parents=True, exist_ok=True)
    sides = {"base": str(task.get("base") or "origin/main"), "branch": "HEAD"}
    changed = changed_functions(wt, sides["base"], "HEAD")
    names = changed["names"]
    rows = []
    with tempfile.TemporaryDirectory(prefix="delfin-realrun-code-") as codes:
        for side, ref in sides.items():
            _archive(wt if side == "branch" else repo, ref, Path(codes) / side)
        for probe in req["probes"]:
            for side in sides:
                # Re-written before every run: each run's trace lands in its own
                # file, so the two sides can never read each other's result.
                d, env_extra = _instrument(Path(codes), side, names)
                row = _one(probe, side, Path(codes) / side, provider, model, key, out,
                           d, env_extra)
                row.update(_read_trace(d, names, row))
                rows.append({"probe": probe["name"], "side": side, **row})
    key = ""
    reached = _reached_overall(rows)
    lines = [f"# Real runs for {run_name}/{tid}", "",
             f"{provider} / {model}; each probe on the base ({sides['base']}) and on the "
             "branch HEAD, fresh HOME and workspace each.", "",
             f"Changed functions the branch diff names: "
             + (f"{len(names)} -- " + ", ".join(sorted(names)[:8])
                + (f" (+{len(names) - 8} more)" if len(names) > 8 else "")
                if names else f"none ({changed['note'] or 'no Python change'})"), "",
             "**reached the change** is whether a function the branch changed was "
             "actually executed in that run: `yes`, `no`, or `unknown` when the run "
             "cannot say (its reason is in the section below).", "",
             "| probe | side | stop | reached the change | offered | called | in tok | "
             "out tok | s | error |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['probe']} | {r['side']} | {r.get('stop')} | "
                     f"{_REACHED_LABEL.get(r.get('reached'), 'unknown')} | "
                     f"{len(r.get('offered') or [])} | {len(r.get('called_tools') or [])} | "
                     f"{r.get('in')} | {r.get('out')} | {r.get('s')} | "
                     f"{(r.get('error') or '').replace('|', '/')[:120]} |")
    for r in rows:
        lines += ["", f"## {r['probe']} / {r['side']}: final answer (last 1200 characters)", "",
                  "```", r.get("answer") or "(none)", "```", "",
                  f"- reached the change: {_REACHED_LABEL.get(r.get('reached'), 'unknown')} "
                  f"-- {r.get('reached_reason') or 'the changed functions ran'}",
                  f"- tools offered in its first turn: "
                  + (f"{len(r['offered'])}: " + ", ".join(r["offered"]) if r.get("offered")
                     else "none recorded"),
                  f"- tools it called: " + (", ".join(r["called_tools"])
                                            if r.get("called_tools") else "none")]
    (out / "summary.md").write_text(_redact("\n".join(lines) + "\n"))
    wave.record(run_name, "realrun", tid,
                rows=[{k: r.get(k) for k in ("probe", "side", "stop", "tools", "error",
                                             "reached")} for r in rows],
                reached=reached, changed=len(names))
    return {"run": run_name, "id": tid, "rows": rows, "reached": reached,
            "changed_functions": len(names), "summary": str(out / "summary.md")}


def _reached_overall(rows: list) -> str:
    """One word for the member: did ANY probe of the branch reach the change?

    Only the branch side counts -- the base has no change to reach. ``unknown``
    outranks ``no`` when nothing reached it, because "the probe missed it" and
    "the run could not tell" are different things to say to a member.
    """
    branch = [r.get("reached") for r in rows if r.get("side") == "branch"]
    if not branch:
        return "unknown"
    if "yes" in branch:
        return "yes"
    return "unknown" if "unknown" in branch else "no"
