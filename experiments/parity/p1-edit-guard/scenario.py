"""Parity measurement: read-before-edit and edit precision, on the product path.

Input: the DELFIN tree this file sits in (the base commit for ``before.json``,
the branch for ``after.json``). Output: one row per scenario variant with
three independent verdicts, plus the exact command that produced them.

The entry point is the product's own, not a helper:
``delfin.agent.cli.main(["-p", ..., "--output-format", "stream-json"])``,
the same one ``delfin-agent chat -p`` calls. The model is the scripted
endpoint in ``tests/agent_terminal_probes/fake_model.py``, whose ``plan``
this script replaces per variant, so a live provider and a key are neither
needed nor used.

Three verdicts per variant, measured apart so a number cannot be one number:

REACHED  the code that is supposed to make this decision ran. Proven by a
         ``sys.settrace`` line-hit trace over ``delfin/agent/api_client.py``:
         the variant names source markers, this script resolves each marker
         to a line at run time and reports whether the trace hit it.
         ``reached_marker`` says which marker carried the verdict, so the
         same script measures a tree where the guard still sits in the tool
         executor and one where it moved into a helper.
USED     the product path actually invoked the tool where it should. Read
         off the ``tool_use`` events of the stream-json output, never off
         the prompt.
CORRECT  the outcome is right: a refusal that says why, an edit that landed
         exactly where it was meant, or an unchanged file. Read off the file
         on disk after the turn and off the tool result the scripted endpoint
         logged.

A correct REFUSAL counts as CORRECT. Nothing here loosens a check.

Usage::

    python experiments/parity/p1-edit-guard/scenario.py --out before.json
"""

from __future__ import annotations

import argparse
import io
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

MARK = "PARITYP1"

#: Files the trace keeps, and the markers are looked for. The refusal a
#: variant expects is not always built in the tool executor: an ambiguous or
#: absent old_string is explained by the edit engine, so a trace of the
#: executor alone would report "not reached" for a decision that was made.
TRACED_FILES = ("delfin/agent/api_client.py", "delfin/agent/edit_engine.py")


# ---------------------------------------------------------------------------
# the line-hit trace
# ---------------------------------------------------------------------------

class Trace:
    """Line numbers hit inside one file, plus the code objects entered."""

    def __init__(self, filenames) -> None:
        self.filenames = tuple(filenames)
        self.lines: dict[str, set[int]] = {name: set() for name in filenames}
        self.entered: set[str] = set()
        self._previous = None

    def _local(self, frame, event, arg):
        if event == "line":
            for name in self.filenames:
                if frame.f_code.co_filename.endswith(name):
                    self.lines[name].add(frame.f_lineno)
                    break
        return self._local

    def _global(self, frame, event, arg):
        if frame.f_code.co_filename.endswith(self.filenames):
            self.entered.add(frame.f_code.co_name)
            return self._local
        return None

    def __enter__(self) -> "Trace":
        # Threads started before this point keep their own trace function, and
        # the tool dispatch runs off the main thread. Both hooks are set, or
        # the trace reports the main thread only and a refusal that happened
        # on a worker thread reads as "not reached".
        import threading
        self._previous = sys.gettrace()
        self._previous_thread = threading.gettrace()
        sys.settrace(self._global)
        threading.settrace(self._global)
        return self

    def __exit__(self, *exc) -> None:
        import threading
        sys.settrace(self._previous)
        threading.settrace(self._previous_thread)


def _bare(name: str) -> str:
    """``mcp__kit-coding__write_file`` -> ``write_file``: the product
    registers its tools under a namespace, the scenario names the tool."""
    return str(name).rsplit("__", 1)[-1] if "__" in str(name) else str(name)


_MARKER_CACHE: dict[str, tuple] = {}
_AST_CACHE: dict[str, tuple] = {}


def _parsed(name: str):
    """(source lines, AST) of one traced file, parsed once per process."""
    if name not in _AST_CACHE:
        import ast
        try:
            source = (REPO / name).read_text(encoding="utf-8")
            _AST_CACHE[name] = (source.splitlines(), ast.parse(source))
        except (OSError, SyntaxError):
            _AST_CACHE[name] = ([], None)
    return _AST_CACHE[name]


def _marker_lines(marker: str) -> tuple:
    """``(file, line)`` for the innermost statements whose source holds it.

    Not the line the marker's text sits on: the guard's refusals are built
    inside multi-line statements, and CPython reports the STATEMENT's lines,
    so a sentence on the third line of a ``return json.dumps({...})`` is
    reported at the ``return``. Resolved through the AST instead.

    INNERMOST on purpose. ``if tracked is None:`` encloses the refusal and is
    reached on every write, guard fired or not -- a statement reached either
    way cannot answer "did this decision run". A statement is kept only when
    no statement nested inside it carries the marker, which leaves the
    ``return`` the refusal is built in.

    Read from disk once per process. Nothing edits the files during a run, so
    the numbers the trace is compared against are the code that ran.
    """
    if marker in _MARKER_CACHE:
        return _MARKER_CACHE[marker]
    import ast
    found: set[tuple] = set()
    for name in TRACED_FILES:
        lines, tree = _parsed(name)
        if tree is None or not any(marker in line for line in lines):
            continue

        def _segment(stmt) -> str:
            end = getattr(stmt, "end_lineno", None) or stmt.lineno
            return "\n".join(lines[stmt.lineno - 1:end])

        for node in ast.walk(tree):
            for stmt in ast.walk(node):
                if not isinstance(stmt, ast.stmt):
                    continue
                if marker not in _segment(stmt):
                    continue
                nested = any(m is not stmt and isinstance(m, ast.stmt)
                             and marker in _segment(m) for m in ast.walk(stmt))
                if not nested:
                    found.add((name, stmt.lineno))
    _MARKER_CACHE[marker] = tuple(sorted(found))
    return _MARKER_CACHE[marker]


# ---------------------------------------------------------------------------
# variants
# ---------------------------------------------------------------------------

def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _rewrite_newer(ws: Path) -> None:
    """Another process writes the file; its mtime moves forward."""
    target = ws / "target.py"
    _write(target, "x = 111\ny = 2\nz = 3\n")


def _rewrite_preserving_mtime(ws: Path) -> None:
    """Another process writes it and restores the old stamp (cp -p, git
    checkout of a tracked file, rsync -t, an editor that keeps mtime)."""
    target = ws / "target.py"
    stamp = target.stat()
    _write(target, "x = 111\ny = 2\n")
    shutil.copystat(target, target)  # no-op; the stamp below is the point
    os.utime(target, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))


def _rewrite_same_length_and_utime(ws: Path) -> None:
    """Same byte length, old stamp put back -- mtime and size both agree."""
    target = ws / "target.py"
    stamp = target.stat()
    _write(target, "x = 999\ny = 2\n")
    os.utime(target, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))


class Variant:
    def __init__(self, name: str, tool: str, expect: str, *,
                 steps: list, files: dict[str, str] | None = None,
                 correct, note: str = "", reached_markers: tuple = (),
                 needs_shell: bool = False) -> None:
        self.name = name
        self.tool = tool
        self.expect = expect
        self.steps = steps
        self.files = files or {"target.py": "x = 1\ny = 2\n"}
        self.correct = correct
        self.note = note
        self.reached_markers = reached_markers
        #: The variant drives a shell command. Where the host's cage refuses
        #: to start one at all, the variant measures the cage and is reported
        #: UNTESTED rather than counted as a capability failure.
        self.needs_shell = needs_shell


def _read_then_write(path: str, content: str, *, mode: str = "") -> dict:
    args = {"path": path, "content": content}
    if mode:
        args["mode"] = mode
    return {"tool": "write_file", "args": args}


def _read(path: str = "target.py") -> dict:
    return {"tool": "read_file", "args": {"path": path}}


_PATCH_REPLACE = (
    "--- a/target.py\n"
    "+++ b/target.py\n"
    "@@ -1,2 +1,2 @@\n"
    "-x = 1\n"
    "+x = 42\n"
    " y = 2\n"
)

#: The notebook a notebook_edit case works on: one markdown cell, one code
#: cell, the smallest file that is a real .ipynb.
_NOTEBOOK = json.dumps({
    "cells": [
        {"cell_type": "markdown", "metadata": {}, "source": ["# t\n"]},
        {"cell_type": "code", "metadata": {}, "source": ["print(1)\n"],
         "execution_count": None, "outputs": []},
    ],
    "metadata": {}, "nbformat": 4, "nbformat_minor": 5,
}) + "\n"


def _unchanged(ws: Path, name: str, original: str, phrase: str):
    def check(result: str) -> bool:
        return (ws / name).read_text(encoding="utf-8") == original \
            and phrase in (result or "")
    return check


def _applied(ws: Path, name: str, text: str):
    def check(result: str) -> bool:
        return (ws / name).read_text(encoding="utf-8") == text \
            and '"error"' not in (result or "")
    return check


def _build_variants(ws: Path) -> list[Variant]:
    two = "x = 1\ny = 2\n"
    return [
        # --- blind overwrite -------------------------------------------------
        Variant(
            "blind_overwrite", "write_file", "refused",
            steps=[_read_then_write("target.py", "x = 9\n")],
            correct=_unchanged(ws, "target.py", two, "refusing to overwrite"),
            note="write_file on an existing file the session never read",
            reached_markers=("refusing to overwrite existing file",),
        ),
        Variant(
            "read_then_overwrite", "write_file", "applied",
            steps=[_read(), _read_then_write("target.py", "x = 9\n")],
            correct=_applied(ws, "target.py", "x = 9\n"),
            note="the control: a read file may be overwritten",
            reached_markers=("refusing to overwrite existing file",
                             "read_tracker[str(full.resolve())] = full.stat"),
        ),
        # --- stale writes ----------------------------------------------------
        Variant(
            "stale_rewrite_newer_mtime", "write_file", "refused",
            steps=[_read(), {"before": _rewrite_newer,
                             **_read_then_write("target.py", "x = 9\n")}],
            correct=_unchanged(ws, "target.py", "x = 111\ny = 2\nz = 3\n",
                               "was modified since last read_file"),
            note="another process wrote it; the stamp moved forward",
            reached_markers=("was modified since last read_file",),
        ),
        Variant(
            "stale_rewrite_mtime_restored", "write_file", "refused",
            steps=[_read(), {"before": _rewrite_preserving_mtime,
                             **_read_then_write("target.py", "x = 9\n")}],
            correct=_unchanged(ws, "target.py", "x = 111\ny = 2\n",
                               "was modified since last read_file"),
            note="another process wrote it and put the old stamp back",
            reached_markers=("was modified since last read_file",),
        ),
        Variant(
            "stale_same_length_utime", "edit_file", "refused",
            steps=[_read(), {"before": _rewrite_same_length_and_utime,
                             "tool": "edit_file", "args": {
                                 "path": "target.py",
                                 "old_string": "x = 999",
                                 "new_string": "x = 5"}}],
            correct=_unchanged(ws, "target.py", "x = 999\ny = 2\n",
                               "was modified since last read_file"),
            note="same byte length, old stamp: size and mtime both agree",
            reached_markers=("was modified since last read_file",),
        ),
        # --- edit precision --------------------------------------------------
        Variant(
            "edit_unique", "edit_file", "applied",
            steps=[_read(), {"tool": "edit_file",
                             "args": {"path": "target.py",
                                      "old_string": "y = 2",
                                      "new_string": "y = 20"}}],
            correct=_applied(ws, "target.py", "x = 1\ny = 20\n"),
            note="one unique match lands where it was meant",
            reached_markers=("old_text.count(old_string)",),
        ),
        Variant(
            "edit_ambiguous", "edit_file", "refused",
            files={"target.py": "x = 1\nx = 1\ny = 2\n"},
            steps=[_read(), {"tool": "edit_file",
                             "args": {"path": "target.py",
                                      "old_string": "x = 1",
                                      "new_string": "x = 5"}}],
            correct=_unchanged(ws, "target.py", "x = 1\nx = 1\ny = 2\n",
                               "matches 2 times"),
            note="two matches: must name them, not pick one",
            reached_markers=("matches {count} times",),
        ),
        Variant(
            "edit_absent", "edit_file", "refused",
            steps=[_read(), {"tool": "edit_file",
                             "args": {"path": "target.py",
                                      "old_string": "y = 77",
                                      "new_string": "y = 78"}}],
            correct=_unchanged(ws, "target.py", two, "not found"),
            note="text that is not there: must say so",
            reached_markers=("not found (no exact match)",),
        ),
        Variant(
            "multi_edit_ambiguous_second", "multi_edit", "refused",
            files={"target.py": "a = 1\nb = 2\nb = 2\n"},
            steps=[_read(), {"tool": "multi_edit", "args": {
                "path": "target.py",
                "edits": [
                    {"old_string": "a = 1", "new_string": "a = 10"},
                    {"old_string": "b = 2", "new_string": "b = 20"},
                ]}}],
            correct=_unchanged(ws, "target.py", "a = 1\nb = 2\nb = 2\n",
                               "Nothing in this batch was applied"),
            note="batch is atomic: one ambiguous edit refuses the whole batch",
            reached_markers=("Nothing in this batch was",),
        ),
        # --- append ----------------------------------------------------------
        Variant(
            "append_after_read", "write_file", "applied",
            steps=[_read(), _read_then_write("target.py", "z = 3\n",
                                             mode="append")],
            correct=_applied(ws, "target.py", "x = 1\ny = 2\nz = 3\n"),
            note="append to a file that was read",
            reached_markers=("build_new_text",),
        ),
        Variant(
            "append_without_read", "write_file", "refused",
            steps=[_read_then_write("target.py", "z = 3\n", mode="append")],
            correct=_unchanged(ws, "target.py", two, "refusing to overwrite"),
            note="append is a write: it needs the same baseline",
            reached_markers=("refusing to overwrite existing file",),
        ),
        # --- the diff route --------------------------------------------------
        Variant(
            "apply_patch_unread", "apply_patch", "refused",
            steps=[{"tool": "apply_patch", "args": {"diff": _PATCH_REPLACE}}],
            correct=_unchanged(ws, "target.py", two, "read_file"),
            note="a diff is an overwrite by another route",
            # The capability being measured: a read-baseline refusal for the
            # files a diff would modify. Absent on a tree that has none, so
            # this reports "not reached" there -- which is the truth.
            reached_markers=("apply_patch needs a read baseline",),
        ),
        Variant(
            "apply_patch_after_read", "apply_patch", "applied",
            steps=[_read(), {"tool": "apply_patch",
                             "args": {"diff": _PATCH_REPLACE}}],
            correct=_applied(ws, "target.py", "x = 42\ny = 2\n"),
            note="the control: a read file may be patched",
            # A `def` line runs at import, before any trace: mark the body.
            reached_markers=("from . import patch_apply as _pa",),
        ),
        # --- working controls that must not regress --------------------------
        Variant(
            "notebook_edit_without_read", "notebook_edit", "refused",
            files={"book.ipynb": _NOTEBOOK},
            steps=[{"tool": "notebook_edit", "args": {
                "path": "book.ipynb", "mode": "replace", "cell_idx": 1,
                "source": "print(2)\n"}}],
            correct=_unchanged(ws, "book.ipynb", _NOTEBOOK,
                               "call notebook_read"),
            note="control: an unread notebook is refused",
            reached_markers=("call notebook_read on",),
        ),
        Variant(
            "shell_read_then_edit", "edit_file", "applied",
            steps=[{"tool": "bash",
                     "args": {"command": "sed -n '1,2p' target.py"}},
                   {"tool": "edit_file", "args": {
                       "path": "target.py",
                       "old_string": "y = 2", "new_string": "y = 22"}}],
            correct=_applied(ws, "target.py", "x = 1\ny = 22\n"),
            note="control: text printed through the shell counts as read",
            # Same rule: the body's first statement, not the `def` line.
            reached_markers=("tracker = getattr(perms, \"read_tracker\", None)",),
            needs_shell=True,
        ),
    ]


# ---------------------------------------------------------------------------
# one variant: one real headless turn
# ---------------------------------------------------------------------------

def _run_variant(variant: Variant, root: Path) -> dict:
    """Drive ``delfin.agent.cli.main`` for one variant; return its row."""
    from tests.agent_terminal_probes import fake_model
    from delfin.agent import api_client as _A
    from delfin.agent import cli as agent_cli
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent import provider_profile as PP

    home = root / "home"
    ws = root / "ws"
    ws.mkdir(parents=True, exist_ok=True)
    for name, text in variant.files.items():
        _write(ws / name, text)

    saved_env = {k: os.environ.get(k) for k in
                 ("HOME", "OLLAMA_HOST", "DELFIN_CLASSIFIER",
                  "DELFIN_SCRATCH_STATE", "http_proxy", "https_proxy",
                  "HTTP_PROXY", "HTTPS_PROXY")}
    saved_cwd = os.getcwd()
    os.environ.update(HOME=str(home), DELFIN_CLASSIFIER="off",
                      DELFIN_PROCESS_GUARD="off")
    os.chdir(ws)

    class _EmptyReg:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    cap = mc.ModelCapabilities(model=fake_model.MODEL, provider="ollama",
                               context_window=200_000, supports_tools=True)
    patches = [
        (mc, "resolve", lambda *a, **k: cap),
        (M, "get_registry", lambda *a, **k: _EmptyReg()),
        (agent_cli, "_save_session", lambda engine, root, **kw: "sid-1"),
        (agent_cli, "_resume_or_create", lambda engine, args: ""),
        (_A, "_resolve_max_tool_rounds", lambda *a, **k: 8),
    ]
    originals = [(mod, name, getattr(mod, name)) for mod, name, _ in patches]
    for mod, name, value in patches:
        setattr(mod, name, value)

    baseline = root / "learned_profiles.json"
    baseline.write_bytes(PP._DEFAULT_PATH.read_bytes())
    _pp_orig = PP._DEFAULT_PATH
    PP._DEFAULT_PATH = baseline

    steps = [dict(s) for s in variant.steps]

    def _plan(messages):
        # The action watcher ("safety classifier") asks the SAME endpoint a
        # second question under its own system prompt. An unreadable answer
        # makes it hold an action nobody can be asked about, which would look
        # like the guard under test refusing. Its header names it, so it is
        # scripted here, and only it: the turn's own steps are matched on the
        # marker the scenario typed.
        joined = " ".join(fake_model._text_of(m) for m in messages)
        if "safety classifier" in joined:
            return {"text": "DECISION: allow\nREASON: an ordinary step "
                            "inside the workspace for the task at hand."}
        starts = [i for i, m in enumerate(messages)
                  if m.get("role") == "user"
                  and MARK in fake_model._text_of(m)]
        done = sum(1 for m in messages[starts[-1] if starts else 0:]
                   if m.get("role") == "assistant")
        if done >= len(steps):
            return {"text": "Nothing more."}
        step = dict(steps[done])
        before = step.pop("before", None)
        if before is not None:
            before(ws)
        return step

    fake_model.plan = _plan
    model_log = root / "model.log"
    srv, port = fake_model.serve(log_path=str(model_log))
    os.environ["OLLAMA_HOST"] = f"http://127.0.0.1:{port}"

    captured: list[str] = []
    real_write = sys.stdout.write

    def _tee(text):
        captured.append(text)
        return real_write(text)

    stdin_saved = sys.stdin
    stdin = io.StringIO("")
    stdin.isatty = lambda: False
    sys.stdin = stdin
    sys.stdout.write = _tee  # type: ignore[method-assign]
    trace = Trace(TRACED_FILES)
    try:
        with trace:
            try:
                agent_cli.main(["-p", f"{MARK} do the scenario",
                                "--cwd", str(ws), "--provider", "ollama",
                                "--model", fake_model.MODEL,
                                "--permission-mode", "default",
                                "--output-format", "stream-json"])
            except SystemExit:
                pass
    finally:
        sys.stdout.write = real_write  # type: ignore[method-assign]
        sys.stdin = stdin_saved
        try:
            srv.shutdown()
        except Exception:
            pass
        for mod, name, value in originals:
            setattr(mod, name, value)
        PP._DEFAULT_PATH = _pp_orig
        os.chdir(saved_cwd)
        for k, v in saved_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    output = "".join(captured)
    events = []
    for line in output.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            events.append(json.loads(line))
        except ValueError:
            continue
    used_tools = [_bare(e.get("name") or "") for e in events
                  if e.get("type") == "tool_use"]

    results: list[str] = []
    if model_log.exists():
        for line in model_log.read_text(encoding="utf-8").splitlines():
            try:
                row = json.loads(line)
            except ValueError:
                continue
            if row.get("tool_result"):
                results.append(row["tool_result"])
    last_result = results[-1] if results else ""

    reached: dict = {}
    for marker in variant.reached_markers:
        cands = _marker_lines(marker)
        reached[marker] = any(name.endswith(tuple(TRACED_FILES))
                              and line in trace.lines.get(name, set())
                              for name, line in cands)
    hit_markers = [m for m, hit in reached.items() if hit]
    shell_blocked = variant.needs_shell and any(
        '"exit_code": 126' in r or "sandbox could not be applied" in r
        or "the cage" in r for r in results)
    if shell_blocked:
        return {
            "variant": variant.name,
            "tool_expected": variant.tool,
            "expect": variant.expect,
            "note": variant.note,
            "untested": ("this host refuses to start a shell command "
                         "(the sandbox answered instead), so the shell path "
                         "was never exercised"),
            "reached": None, "used": 1 if variant.tool in used_tools else 0,
            "correct": None, "reached_marker": "",
            "used_tools": used_tools,
            "tool_result": (last_result or "")[:400],
            "trace_entered": sorted(trace.entered)[:40],
        }

    try:
        correct = bool(variant.correct(last_result))
    except Exception as exc:            # a missing file is an outcome too
        correct = False
        last_result = f"{last_result}\n[{type(exc).__name__}: {exc}]"

    return {
        "variant": variant.name,
        "tool_expected": variant.tool,
        "expect": variant.expect,
        "note": variant.note,
        "reached": 1 if hit_markers else 0,
        "reached_marker": hit_markers[0] if hit_markers else "",
        "used": 1 if variant.tool in used_tools else 0,
        "used_tools": used_tools,
        "correct": 1 if correct else 0,
        "tool_result": (last_result or "")[:400],
        "trace_entered": sorted(trace.entered)[:40],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="before.json",
                        help="where the rows are written (JSON)")
    parser.add_argument("--only", default="",
                        help="comma-separated variant names (default: all)")
    args = parser.parse_args(argv)

    probe = Path(tempfile.mkdtemp(prefix="p1-layout-"))
    names = {v.name for v in _build_variants(probe)}
    shutil.rmtree(probe, ignore_errors=True)
    wanted = [n.strip() for n in args.only.split(",") if n.strip()]
    unknown = sorted(set(wanted) - names)
    if unknown:
        print(f"unknown variant(s): {', '.join(unknown)}", file=sys.stderr)
        return 2

    rows: list[dict] = []
    base = Path(tempfile.mkdtemp(prefix="p1-parity-"))
    # The variants' correctness closures need the workspace, so they are
    # built per run against that run's own ``ws``.
    for name in sorted(names):
        if wanted and name not in wanted:
            continue
        run_root = base / name
        run_root.mkdir(parents=True, exist_ok=True)
        ws = run_root / "ws"
        ws.mkdir(exist_ok=True)
        variants = {v.name: v for v in _build_variants(ws)}
        row = _run_variant(variants[name], run_root)
        rows.append(row)
        print(f"{row['variant']:<34} "
              + ("UNTESTED " + row["untested"] if row.get("untested") else
                 f"REACHED={row['reached']} USED={row['used']} "
                 f"CORRECT={row['correct']}"), flush=True)

    measured = [r for r in rows if r.get("reached") is not None]
    totals = {
        "variants": len(rows),
        "measured": len(measured),
        "untested": [r["variant"] for r in rows
                     if r.get("reached") is None],
        "reached": sum(r["reached"] for r in measured),
        "used": sum(r["used"] for r in measured),
        "correct": sum(r["correct"] for r in measured),
    }
    import delfin as _delfin
    from delfin.agent import api_client as _api
    payload = {
        "package": "p1-edit-guard",
        "provenance": {
            "delfin_package": _delfin.__file__,
            "api_client": _api.__file__,
            "repo": str(REPO),
        },
        "command": ("python experiments/parity/p1-edit-guard/scenario.py "
                    f"--out {args.out}"),
        "entry_point": ("delfin.agent.cli.main([\"-p\", ..., "
                        "\"--output-format\", \"stream-json\"])"),
        "model": "tests/agent_terminal_probes/fake_model.py (scripted plan)",
        "totals": totals,
        "rows": rows,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n",
                              encoding="utf-8")
    print(f"\n{totals['reached']}/{totals['measured']} REACHED  "
          f"{totals['used']}/{totals['measured']} USED  "
          f"{totals['correct']}/{totals['measured']} CORRECT"
          + (f"  ({len(totals['untested'])} UNTESTED: "
             f"{', '.join(totals['untested'])})" if totals["untested"] else ""))
    shutil.rmtree(base, ignore_errors=True)
    return 0 if totals["correct"] == totals["measured"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
