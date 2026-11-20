#!/usr/bin/env python3
"""Score one attempt at an unreachable-capability task, blind to the arm.

Input: one attempt directory named by an opaque id, a checkout of the tree
with the attempt applied, and a checkout of the same repository at the
pre-attempt sha (the control). Output: one JSON record, printed and
optionally written, holding the four checks, the route, and the cost figures.

The score is 1 only when all four checks hold, and 0 otherwise. A failure, a
time-out and a refusal all score 0 and stay in the denominator; nothing is
dropped, and the reason is recorded beside the zero.

The four checks
---------------
1. ``capability_resolved`` -- the capability is no longer in the product's
   unreachable set, re-measured on the attempt tree with the scan in
   ``ab_task_pools.py``. The outcome is a **route**, never a sum:
   ``connected`` (still defined, and something under ``delfin/`` now names
   it) or ``deleted`` (no definition of that name is left under
   ``delfin/``). ``still_unreachable`` and ``deleted_but_named`` are the two
   ways of not finishing.
2. ``tests_still_pass`` -- the consuming suites of the touched files, chosen
   from the tree **before** the attempt, run in both trees. Per test:
   passed, failed, skipped, missing. Charged to the attempt are only the
   tests that passed in the control, so a test already red at the baseline
   is not the attempt's debt; a test that passed before and skips now is
   recorded as ``turned_into_skip`` and fails the check, because a skip is
   not a pass. The control is required: without it a skip and a
   pre-existing failure cannot be told apart from a regression.
3. ``instrument_untouched`` -- ``git diff --name-only <pre>..<post>`` names
   none of: the reachability test, its baseline fixture, any
   ``scripts/ab_*.py``, any ``conftest.py``. One exception, with its own
   check: the baseline fixture may lose lines and may not gain any, because
   the reachability test fails until the connected or deleted name is taken
   out of it -- so a correct attempt is forced to edit that file. A line
   added to it, or a removed line naming anything other than this task's
   capability, scores zero.
4. ``deletion_removed_its_tests`` -- when the route is ``deleted``, no test
   file that named the capability before the attempt may still name it. A
   deleted function whose tests were left behind red is not a finished task.
   Vacuously true on the ``connected`` route.

Why the scorer cannot see the arm
---------------------------------
* It takes no arm, role or model argument. ``build_parser`` has no such
  option and ``tests/test_ab_apparatus.py`` asserts that the option strings
  contain none of the arm vocabulary.
* It opens exactly two files inside the attempt directory, both by fixed
  name: ``manifest.json`` and ``log.jsonl``. Every other file there is
  listed in ``ignored_files`` and never read, so the runner may drop
  transcripts and briefs beside them without reaching the score.
* The manifest is read against a key allowlist (``MANIFEST_KEYS``,
  ``TASK_KEYS``). An unknown key makes the scorer refuse with exit 3 and
  emit no score. A runner that tries to hand it the arm therefore fails
  loudly instead of quietly biasing a number.
* The log is **projected** onto the numeric usage fields before anything
  else looks at it (``LOG_KEYS``): prompt text, system prompts and tool
  arguments are discarded at the read, not filtered later.
* The arm lives in a mapping file of attempt id to arm that this script has
  no path to and never globs for. Joining the scores to the arms is a step
  after all attempts are scored.

A refusal under any of those rules is an apparatus defect, not an attempt
outcome: it is reported with exit 3 and no score, because scoring zero would
charge the agent for the runner's mistake.

Cost axis
---------
``input_tokens``, ``output_tokens``, ``wall_seconds`` and ``turns`` come from
the attempt's own log and manifest. A figure that is not in the log is
``null``. It is never 0: 0 tokens and an unreported token count are
different facts, and summing them would understate the cost of the arm whose
logs are thinner.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
import json
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ab_task_pools import (  # noqa: E402  (same directory, one scan only)
    _chem_file_cache,
    _definitions,
    _mentions,
    pool_roots,
)

#: Files that decide the score. An attempt that edits one scores zero.
INSTRUMENT_PATHS = (
    "tests/test_every_capability_has_a_caller.py",
    "tests/fixtures/unmounted_capabilities.txt",
)
INSTRUMENT_GLOBS = ("scripts/ab_*.py",)
INSTRUMENT_BASENAMES = ("conftest.py",)

#: The one instrument file a correct attempt is forced to change, and only
#: by removing lines (see check 3).
BASELINE_FIXTURE = "tests/fixtures/unmounted_capabilities.txt"

MANIFEST_KEYS = frozenset({
    "attempt_id", "task", "pre_sha", "post_sha", "started_at", "ended_at",
})
TASK_KEYS = frozenset({
    "name", "file", "line", "module", "pool", "tier", "task",
})
#: The only log fields that reach the score.
LOG_KEYS = frozenset({
    "ts", "input_tokens", "output_tokens", "type", "role",
})
#: Vocabulary that would reveal the arm. Checked against the manifest.
ARM_WORDS = (
    "arm", "role", "specialist", "generalist", "variant", "condition",
    "treatment", "control_group", "model", "prompt", "brief", "subagent",
    "delegate", "persona",
)

_WORDRE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


class Refusal(Exception):
    """An apparatus defect: no score is emitted."""


# --------------------------------------------------------------------------
# attempt artefacts
# --------------------------------------------------------------------------

def read_manifest(attempt_dir: Path) -> dict:
    """The attempt manifest, against the key allowlist.

    Input: the attempt directory. Output: the manifest dict. Raises
    ``Refusal`` on an unknown key or on arm vocabulary anywhere in it.
    """
    path = attempt_dir / "manifest.json"
    if not path.is_file():
        raise Refusal(f"no manifest.json in {attempt_dir}")
    doc = json.loads(path.read_text())
    if not isinstance(doc, dict):
        raise Refusal("manifest.json is not an object")
    unknown = sorted(set(doc) - MANIFEST_KEYS)
    if unknown:
        raise Refusal(f"manifest holds keys the scorer may not read: {unknown}")
    task = doc.get("task")
    if not isinstance(task, dict):
        raise Refusal("manifest.task is missing")
    unknown = sorted(set(task) - TASK_KEYS)
    if unknown:
        raise Refusal(f"manifest.task holds keys the scorer may not read: "
                      f"{unknown}")
    blob = json.dumps(doc).lower()
    leaked = sorted({w for w in ARM_WORDS
                     if re.search(rf"\b{re.escape(w)}\b", blob)})
    if leaked:
        raise Refusal(f"the manifest names the arm or the agent: {leaked}")
    for key in ("pre_sha", "post_sha"):
        if not doc.get(key):
            raise Refusal(f"manifest.{key} is missing")
    return doc


def read_cost(attempt_dir: Path) -> dict:
    """The cost figures, projected onto ``LOG_KEYS`` at the read.

    Output: ``{"input_tokens", "output_tokens", "turns", "wall_seconds",
    "log_records", "dropped_fields"}``; a figure absent from the log is
    ``None``, never 0. ``dropped_fields`` counts the log fields the
    projection discarded -- the count, so that the projection is visible
    from outside without the discarded names coming with it.
    """
    out = {"input_tokens": None, "output_tokens": None, "turns": None,
           "wall_seconds": None, "log_records": 0, "dropped_fields": 0}
    path = attempt_dir / "log.jsonl"
    if not path.is_file():
        return out
    rows = []
    dropped = 0
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(rec, dict):
            usage = rec.get("usage") if isinstance(rec.get("usage"), dict) \
                else {}
            merged = {**usage, **rec}
            projected = {k: merged[k] for k in LOG_KEYS if k in merged}
            # The count, never the names: a key name can carry the arm as
            # surely as a value can, and the count is what shows that the
            # projection happened at all.
            dropped += len(merged) - len(projected)
            rows.append(projected)
    out["log_records"] = len(rows)
    out["dropped_fields"] = dropped
    for key in ("input_tokens", "output_tokens"):
        vals = [r[key] for r in rows
                if isinstance(r.get(key), (int, float))]
        out[key] = sum(vals) if vals else None
    turns = [r for r in rows
             if r.get("type") == "turn" or r.get("role") == "assistant"]
    out["turns"] = len(turns) if turns else None
    stamps = sorted(r["ts"] for r in rows if isinstance(r.get("ts"), str))
    if len(stamps) >= 2:
        try:
            first = dt.datetime.fromisoformat(stamps[0])
            last = dt.datetime.fromisoformat(stamps[-1])
            out["wall_seconds"] = round((last - first).total_seconds(), 1)
        except ValueError:
            pass
    return out


def ignored_files(attempt_dir: Path) -> list[str]:
    """Everything in the attempt directory the scorer does not open."""
    read = {"manifest.json", "log.jsonl"}
    return sorted(p.name for p in attempt_dir.iterdir()
                  if p.name not in read)


# --------------------------------------------------------------------------
# git helpers
# --------------------------------------------------------------------------

def _git(tree: Path, *args: str, check: bool = True) -> str:
    out = subprocess.run(["git", "-C", str(tree), *args],
                         capture_output=True, text=True, timeout=300)
    if check and out.returncode != 0:
        raise Refusal(f"git {' '.join(args)} failed in {tree}: "
                      f"{out.stderr.strip()}")
    return out.stdout


def changed_files(tree: Path, pre: str, post: str) -> list[str]:
    return sorted(p for p in _git(tree, "diff", "--name-only",
                                  f"{pre}..{post}").split("\n") if p)


# --------------------------------------------------------------------------
# check 1: the capability left the unreachable set
# --------------------------------------------------------------------------

def _defined_anywhere(root: Path, name: str) -> list[str]:
    """Every file under ``delfin/`` that defines a callable called ``name``.

    Any nesting, not only module level: a function moved into a class is
    not a deletion, and a module-level-only answer would call it one.
    """
    hits = []
    for path in sorted((root / "delfin").rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                    and node.name == name:
                hits.append(f"{path.relative_to(root)}:{node.lineno}")
    return hits


def _call_sites(root: Path, name: str) -> list[str]:
    """Where ``name`` is actually *called* in the product, by syntax tree.

    Input: the tree and the capability name. Output: ``file:line`` for every
    ``ast.Call`` under ``delfin/``, outside test files, whose callee is that
    name plainly or as an attribute.

    Why it is recorded and does not decide the score: the scan counts a name
    in a string literal as a reference, because name-based dispatch is real
    in this product and a scan that ignored it called two live functions
    unreachable. The same rule means a string literal naming the capability
    is enough to move it out of the unreachable set -- and a string literal
    is something the agent can write. This number is what lets that be
    checked afterwards instead of assumed, for both arms.
    """
    out = []
    for path in sorted((root / "delfin").rglob("*.py")):
        if path.name.startswith("test_") or path.name.endswith("_test.py"):
            continue
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            hit = (isinstance(func, ast.Name) and func.id == name) or \
                  (isinstance(func, ast.Attribute) and func.attr == name)
            if hit:
                out.append(f"{path.relative_to(root)}:{node.lineno}")
    return out


def capability_route(root: Path, task: dict) -> dict:
    """Where the capability stands in the attempt tree.

    Input: the attempt tree and the task row. Output: ``{"route", "mentions",
    "call_sites", "definitions"}`` with route in ``connected``, ``deleted``,
    ``still_unreachable``, ``deleted_but_named``. ``call_sites`` is a
    diagnosis, not part of the score (see ``_call_sites``).
    """
    name = task["name"]
    defs_here = _defined_anywhere(root, name)
    pool = task.get("pool") or "infra"
    chem = _chem_file_cache(root)
    def_root, keep = pool_roots(root, chem)[pool]
    scoped = _definitions(def_root, root, keep)
    defs = {name: scoped.get(name, (task["file"], task["line"]))}
    skip = {root / BASELINE_FIXTURE}
    mentions = _mentions(root / "delfin", defs, True, root, skip)[name]
    if defs_here:
        route = "connected" if mentions else "still_unreachable"
    else:
        route = "deleted_but_named" if mentions else "deleted"
    return {"route": route, "mentions": mentions,
            "call_sites": _call_sites(root, name), "definitions": defs_here}


# --------------------------------------------------------------------------
# check 2: the consuming suites of the touched files
# --------------------------------------------------------------------------

_UNSEARCHABLE_STEMS = {"__init__", "utils", "main", "config", "cli"}


def consuming_suites(pre_tree: Path, touched: list[str],
                     limit: int = 12) -> list[str]:
    """The test files that cover ``touched``, read from the pre-attempt tree.

    Mirrors the product's own rule (``_related_test_files`` in
    ``delfin/agent/api_client.py``): the conventional test file of each
    touched module first, then the test files whose import lines name a
    touched module. Re-implemented rather than imported because that rule
    lives in a 14,000-line module whose import pulls in the agent runtime,
    and because the list has to be taken from a tree that is not the one
    this process runs in.

    Input: the control tree and the paths the attempt changed. Output: paths
    relative to the tree, most direct first, at most ``limit``.
    """
    found: list[str] = []
    stems: list[str] = []
    for rel in touched:
        p = Path(rel)
        if p.suffix != ".py":
            continue
        if p.name.startswith("test_") or p.name.endswith("_test.py"):
            if (pre_tree / rel).is_file() and rel not in found:
                found.append(rel)
            continue
        for cand in (f"tests/test_{p.stem}.py",
                     f"tests/{p.parent}/test_{p.stem}.py",
                     f"tests/{p.stem}_test.py",
                     str(p.parent / f"test_{p.stem}.py")):
            if (pre_tree / cand).is_file():
                if cand not in found:
                    found.append(cand)
                break
        if p.stem not in _UNSEARCHABLE_STEMS:
            stems.append(p.stem)
    if stems:
        names = "|".join(re.escape(s) for s in dict.fromkeys(stems))
        importer = re.compile(
            rf"^[ \t]*(?:from|import)[ \t][^\n]*\b(?:{names})\b", re.MULTILINE)
        for cand in sorted((pre_tree / "tests").rglob("test_*.py")):
            if len(found) >= limit:
                break
            rel = str(cand.relative_to(pre_tree))
            if rel in found:
                continue
            try:
                text = cand.read_text(errors="replace")
            except OSError:
                continue
            if importer.search(text[:200_000]):
                found.append(rel)
    return found[:limit]


def run_suites(tree: Path, files: list[str], python: str,
               timeout: int = 1800) -> dict:
    """Run ``files`` in ``tree`` and read the junit report.

    Output: ``{"outcomes": {nodeid: passed|failed|skipped}, "node_file":
    {nodeid: path}, "counts": {...}, "present": [...], "missing": [...],
    "returncode", "timed_out"}``. A time-out is reported, not raised: it
    scores zero like any other failure.
    """
    present = [f for f in files if (tree / f).is_file()]
    missing = [f for f in files if f not in present]
    out = {"outcomes": {}, "node_file": {},
           "counts": {"passed": 0, "failed": 0, "skipped": 0},
           "present": present, "missing": missing, "returncode": {},
           "timed_out": False}
    # One run per file. The junit report's ``file`` attribute is empty in
    # the default report family, and a node has to be attributed to its
    # file: otherwise a test file a deletion was allowed to remove cannot
    # be told from a test that quietly vanished.
    for rel in present:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "junit.xml"
            cmd = [python, "-m", "pytest", "-q", "--no-header",
                   "-p", "no:cacheprovider", f"--junit-xml={report}", rel]
            try:
                proc = subprocess.run(cmd, cwd=str(tree), capture_output=True,
                                      text=True, timeout=timeout)
                out["returncode"][rel] = proc.returncode
            except subprocess.TimeoutExpired:
                out["timed_out"] = True
                continue
            if not report.is_file():
                continue
            for case in ET.parse(report).getroot().iter("testcase"):
                node = (f"{rel}::{case.get('classname', '')}"
                        f"::{case.get('name', '')}")
                kinds = {child.tag for child in case}
                if kinds & {"failure", "error"}:
                    verdict = "failed"
                elif "skipped" in kinds:
                    verdict = "skipped"
                else:
                    verdict = "passed"
                out["outcomes"][node] = verdict
                out["node_file"][node] = rel
                out["counts"][verdict] += 1
    return out


def compare_suites(pre: dict, post: dict, allowed_missing: list[str]) -> dict:
    """What the attempt did to the tests that passed before it.

    Input: two ``run_suites`` results and the test files a deletion was
    allowed to remove. Output the four counts and the node ids behind them.
    A test in a file the deletion was allowed to remove is not counted as
    vanished: check 4 is what judges that file, and charging it here as
    well would make the only correct deletion impossible to pass.
    """
    allowed = set(allowed_missing)
    regressed, turned_skip, vanished = [], [], []
    for node, before in pre["outcomes"].items():
        if before != "passed":
            continue
        after = post["outcomes"].get(node)
        if after == "failed":
            regressed.append(node)
        elif after == "skipped":
            turned_skip.append(node)
        elif after is None and pre["node_file"].get(node) not in allowed:
            vanished.append(node)
    already_red = sorted(n for n, v in pre["outcomes"].items()
                         if v == "failed")
    unexpected_missing = sorted(set(post["missing"]) - set(allowed_missing))
    return {
        "passed_before": sum(1 for v in pre["outcomes"].values()
                             if v == "passed"),
        "counts_before": pre["counts"],
        "counts_after": post["counts"],
        "regressed": sorted(regressed),
        "turned_into_skip": sorted(turned_skip),
        "vanished": sorted(vanished),
        "already_red_before": already_red,
        "missing_files": unexpected_missing,
        "timed_out": bool(pre["timed_out"] or post["timed_out"]),
    }


# --------------------------------------------------------------------------
# check 3: the instrument
# --------------------------------------------------------------------------

def instrument_verdict(tree: Path, pre: str, post: str, touched: list[str],
                       name: str) -> dict:
    """Whether the attempt changed anything that decides its own score.

    Output: ``{"ok", "touched_instrument", "fixture_change"}``. A
    removal-only change to the baseline fixture, naming this capability and
    nothing else, is permitted and reported; anything else is a zero.
    """
    hit = []
    for rel in touched:
        p = Path(rel)
        if rel in INSTRUMENT_PATHS or p.name in INSTRUMENT_BASENAMES:
            hit.append(rel)
            continue
        if any(p.match(g) for g in INSTRUMENT_GLOBS):
            hit.append(rel)
    fixture_change = None
    if BASELINE_FIXTURE in hit:
        diff = _git(tree, "diff", "--unified=0", f"{pre}..{post}", "--",
                    BASELINE_FIXTURE)
        added = [ln[1:] for ln in diff.split("\n")
                 if ln.startswith("+") and not ln.startswith("+++")]
        removed = [ln[1:] for ln in diff.split("\n")
                   if ln.startswith("-") and not ln.startswith("---")]
        removed_names = sorted({
            parts[1] for parts in (ln.split("\t") for ln in removed)
            if len(parts) >= 2})
        # A removed header line carries no tab, so it would otherwise be
        # invisible here: the header is what says what the tiers mean.
        removed_prose = [ln for ln in removed
                         if ln.strip() and len(ln.split("\t")) < 2]
        allowed = (not added and not removed_prose
                   and removed_names == [name])
        fixture_change = {"added_lines": len(added),
                          "removed_names": removed_names,
                          "removed_prose_lines": len(removed_prose),
                          "allowed": allowed}
        if allowed:
            hit.remove(BASELINE_FIXTURE)
    return {"ok": not hit, "touched_instrument": sorted(hit),
            "fixture_change": fixture_change}


# --------------------------------------------------------------------------
# check 4: a deletion took its tests with it
# --------------------------------------------------------------------------

def tests_naming(tree: Path, name: str) -> list[str]:
    """Test files whose text names ``name``.

    Both ``tests/`` and the test files that sit beside the source under
    ``delfin/``: five of those exist, and a deletion that left one of them
    behind would otherwise pass.
    """
    out = []
    candidates = list((tree / "tests").rglob("*.py"))
    candidates += [p for p in (tree / "delfin").rglob("*.py")
                   if p.name.startswith("test_") or p.name.endswith("_test.py")]
    for path in sorted(candidates):
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        if name in _WORDRE.findall(text):
            out.append(str(path.relative_to(tree)))
    return out


def deletion_verdict(pre_tree: Path, post_tree: Path, name: str,
                     route: str) -> dict:
    """Whether a deletion removed the capability's own tests.

    Vacuous on any route but ``deleted``. Output: ``{"ok", "before",
    "left_behind"}``.
    """
    before = tests_naming(pre_tree, name)
    if route != "deleted":
        return {"ok": True, "applies": False, "before": before,
                "left_behind": []}
    after = tests_naming(post_tree, name)
    return {"ok": not after, "applies": True, "before": before,
            "left_behind": after}


# --------------------------------------------------------------------------
# the score
# --------------------------------------------------------------------------

def score_attempt(attempt_dir: Path, post_tree: Path, pre_tree: Path,
                  python: str, timeout: int = 1800) -> dict:
    """Score one attempt. Input and output as in the module docstring."""
    manifest = read_manifest(attempt_dir)
    task = manifest["task"]
    name = task["name"]
    pre_sha, post_sha = manifest["pre_sha"], manifest["post_sha"]

    touched = changed_files(post_tree, pre_sha, post_sha)
    route = capability_route(post_tree, task)
    instrument = instrument_verdict(post_tree, pre_sha, post_sha, touched,
                                    name)
    deletion = deletion_verdict(pre_tree, post_tree, name, route["route"])

    suites = consuming_suites(pre_tree, touched)
    allowed_missing = deletion["before"] if route["route"] == "deleted" else []
    pre_run = run_suites(pre_tree, suites, python, timeout)
    post_run = run_suites(post_tree, suites, python, timeout)
    tests = compare_suites(pre_run, post_run, allowed_missing)

    checks = {
        "capability_resolved": route["route"] in ("connected", "deleted"),
        "tests_still_pass": not (tests["regressed"] or
                                 tests["turned_into_skip"] or
                                 tests["vanished"] or
                                 tests["missing_files"] or
                                 tests["timed_out"]),
        "instrument_untouched": instrument["ok"],
        "deletion_removed_its_tests": deletion["ok"],
    }
    reasons = []
    if not checks["capability_resolved"]:
        reasons.append(f"route={route['route']}")
    if not checks["tests_still_pass"]:
        for key in ("regressed", "turned_into_skip", "vanished",
                    "missing_files"):
            if tests[key]:
                reasons.append(f"{key}={tests[key]}")
        if tests["timed_out"]:
            reasons.append("suite timed out")
    if not checks["instrument_untouched"]:
        reasons.append(f"instrument touched: {instrument['touched_instrument']}")
    if not checks["deletion_removed_its_tests"]:
        reasons.append(f"tests left behind: {deletion['left_behind']}")

    return {
        "attempt_id": manifest.get("attempt_id") or attempt_dir.name,
        "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(
            timespec="seconds"),
        "task": task,
        "pre_sha": pre_sha,
        "post_sha": post_sha,
        "score": 1 if all(checks.values()) else 0,
        "checks": checks,
        "reasons": reasons,
        "route": route["route"],
        "route_detail": route,
        "touched_files": touched,
        "consuming_suites": suites,
        "tests": tests,
        "instrument": instrument,
        "deletion": deletion,
        "cost": read_cost(attempt_dir),
        "ignored_files": ignored_files(attempt_dir),
    }


def build_parser() -> argparse.ArgumentParser:
    """The command line. It has no arm, role or model option by design."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--attempt-dir", required=True,
                    help="one attempt directory, named by an opaque id")
    ap.add_argument("--post-tree", required=True,
                    help="checkout with the attempt applied")
    ap.add_argument("--pre-tree", required=True,
                    help="checkout at the pre-attempt sha, the control")
    ap.add_argument("--python", default=sys.executable,
                    help="interpreter that runs the suites")
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--out", default=None, help="write the record here")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        record = score_attempt(Path(args.attempt_dir), Path(args.post_tree),
                               Path(args.pre_tree), args.python,
                               args.timeout)
    except Refusal as exc:
        print(json.dumps({"refused": str(exc),
                          "note": "apparatus defect, not an attempt outcome; "
                                  "no score emitted"}, indent=2))
        return 3
    text = json.dumps(record, indent=2, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
