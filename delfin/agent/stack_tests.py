"""Which test runner a ``run_tests`` target belongs to, and what it said.

Input: the workspace and the ``target`` string of a run_tests call.
Output: ``plan(workspace, target)`` -> a list of ``Group``s, one per
(stack, project directory): ``framework`` (pytest | cargo | go | npm |
node), ``project`` (the directory holding the stack's marker), ``targets``
(the parts of ``target`` inside it) and, for the non-Python stacks, the
shell ``command`` that runs them; ``summarise(framework, exit_code,
output)`` -> the run_tests summary (passed, failed, errors, skipped,
xfailed) parsed from the runner's own output.

Detection, per target, from the target itself and the markers between it
and the workspace root, nearest first: a ``.rs`` path or Cargo.toml ->
cargo; a ``.go`` path or go.mod -> go; a ``.ts/.js`` path or package.json
with ``scripts.test`` -> npm (node --test when package.json has no test
script); anything else -> pytest, as before, so an existing call does
not change framework.

The commands are the stack's plain test verbs (``cargo test --offline``,
``go test -v ./...``, ``npm test``, ``node --test``). They are run through
the bash tool, i.e. through its gate: inside the process cage they are
admitted unasked (caged_test_runners), outside it they are asked or
refused exactly as the same command typed into bash would be.

Rejected: running them from here with subprocess, which would be a second
execution path beside the gated shell. Pure; never raises.
"""

from __future__ import annotations

import json
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

_SUFFIX_STACK = {".rs": "cargo", ".go": "go", ".ts": "npm", ".tsx": "npm",
                 ".js": "npm", ".mjs": "npm", ".cjs": "npm"}


@dataclass
class Group:
    framework: str
    project: Path
    targets: list = field(default_factory=list)
    command: str = ""


def _inside(p: Path, ws: Path) -> bool:
    try:
        p, ws = p.resolve(), ws.resolve()
    except OSError:
        return False
    return p == ws or ws in p.parents


def _package_test_script(d: Path) -> Optional[bool]:
    """True: package.json with scripts.test; False: without; None: none."""
    pj = d / "package.json"
    if not pj.is_file():
        return None
    try:
        data = json.loads(pj.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    scripts = (data or {}).get("scripts") if isinstance(data, dict) else None
    return bool(isinstance(scripts, dict) and scripts.get("test"))


def _detect(start: Path, ws: Path, prefer: str = "") -> tuple[str, Path]:
    """(framework, project dir) for a path, walking up to the workspace."""
    cur = start if start.is_dir() else start.parent
    while _inside(cur, ws):
        found = []
        if (cur / "Cargo.toml").is_file():
            found.append("cargo")
        if (cur / "go.mod").is_file():
            found.append("go")
        pkg = _package_test_script(cur)
        if pkg is not None:
            found.append("npm" if pkg else "node")
        if any((cur / m).exists() for m in ("pyproject.toml", "pytest.ini",
                                            "setup.cfg", "conftest.py",
                                            "setup.py", "tox.ini")):
            found.append("pytest")
        if found:
            if prefer:
                for f in found:
                    if f == prefer or (prefer == "npm" and f == "node"):
                        return f, cur
            return found[0], cur
        if cur.resolve() == ws.resolve():
            break
        cur = cur.parent
    return "pytest", ws


def _command(g: Group, ws: Path) -> str:
    rels = []
    for t in g.targets:
        p = (ws / t).resolve() if not Path(t).is_absolute() else Path(t)
        try:
            rels.append(str(p.relative_to(g.project.resolve())) or ".")
        except ValueError:
            rels.append(t)
    rels = [r for r in rels if r not in ("", ".")]
    if g.framework == "cargo":
        extra = ""
        for r in rels:
            if r.startswith("tests/") and r.endswith(".rs"):
                extra += f" --test {shlex.quote(Path(r).stem)}"
        return "cargo test --offline" + extra
    if g.framework == "go":
        pkgs = []
        for r in rels:
            d = r if not r.endswith(".go") else (str(Path(r).parent) or ".")
            pkgs.append("./" + d.strip("./") + "/..." if d.strip("./")
                        else "./...")
        return "go test -v " + " ".join(dict.fromkeys(pkgs or ["./..."]))
    if g.framework == "npm":
        return "npm test"
    if g.framework == "node":
        files = [r for r in rels if re.search(r"\.[cm]?[jt]sx?$", r)]
        return "node --test" + ("".join(f" {shlex.quote(f)}" for f in files)
                                if files else "")
    return ""


def plan(workspace, target: str) -> list[Group]:
    """The runs a target stands for (module docstring)."""
    ws = Path(str(workspace))
    text = str(target or "").strip()
    try:
        parts = shlex.split(text) if text else []
    except ValueError:
        parts = [text]
    groups: dict[tuple[str, str], Group] = {}
    if not parts:
        fw, proj = _detect(ws, ws)
        groups[(fw, str(proj))] = Group(fw, proj, [])
    for part in parts:
        head = part.split("::", 1)[0]
        p = Path(head) if Path(head).is_absolute() else ws / head
        if not p.exists():
            # Not a path (a pytest -k-like id, a typo): pytest decides.
            groups.setdefault(("pytest", str(ws)), Group("pytest", ws, []))
            groups[("pytest", str(ws))].targets.append(part)
            continue
        prefer = _SUFFIX_STACK.get(p.suffix.lower(), "")
        if p.suffix == ".py":
            prefer = "pytest"
        fw, proj = _detect(p, ws, prefer)
        if p.suffix == ".py":
            fw, proj = "pytest", ws
        g = groups.setdefault((fw, str(proj)), Group(fw, proj, []))
        g.targets.append(part)
    out = list(groups.values())
    for g in out:
        if g.framework != "pytest":
            g.command = _command(g, ws)
    return out


_CARGO_RE = re.compile(r"test result: \w+\. (\d+) passed; (\d+) failed; (\d+) ignored")
_GO_LINE_RE = re.compile(r"^\s*--- (PASS|FAIL|SKIP): ", re.M)
_NODE_RE = {k: re.compile(rf"^\s*[ℹ#] {k} (\d+)", re.M)
            for k in ("tests", "pass", "fail", "skipped", "todo", "cancelled")}
_VITEST_RE = re.compile(r"^\s*Tests\s+(.*)$", re.M)
_JEST_RE = re.compile(r"^Tests:\s+(.*)$", re.M)
_MOCHA_RE = {k: re.compile(rf"^\s*(\d+) {k}\b", re.M)
             for k in ("passing", "failing", "pending")}


def summarise(framework: str, exit_code: Optional[int], output: str) -> dict:
    """passed / failed / errors / skipped / xfailed from a runner's output.

    A non-zero exit with nothing counted is one error (the build or the
    collection failed before a test ran)."""
    out = {"passed": 0, "failed": 0, "errors": 0, "skipped": 0, "xfailed": 0}
    text = str(output or "")
    if framework == "cargo":
        for p, f, i in _CARGO_RE.findall(text):
            out["passed"] += int(p)
            out["failed"] += int(f)
            out["skipped"] += int(i)
    elif framework == "go":
        for kind in _GO_LINE_RE.findall(text):
            out[{"PASS": "passed", "FAIL": "failed", "SKIP": "skipped"}[kind]] += 1
    else:
        node = {k: rx.findall(text) for k, rx in _NODE_RE.items()}
        if node["tests"]:
            out["passed"] = sum(int(x) for x in node["pass"])
            out["failed"] = sum(int(x) for x in node["fail"]) + sum(
                int(x) for x in node["cancelled"])
            out["skipped"] = sum(int(x) for x in node["skipped"]) + sum(
                int(x) for x in node["todo"])
        else:
            line = (_VITEST_RE.findall(text) or _JEST_RE.findall(text) or [""])[-1]
            if line:
                for n, word in re.findall(r"(\d+) (passed|failed|skipped|todo)",
                                          line):
                    key = {"passed": "passed", "failed": "failed"}.get(
                        word, "skipped")
                    out[key] += int(n)
            else:
                m = {k: rx.findall(text) for k, rx in _MOCHA_RE.items()}
                out["passed"] = sum(int(x) for x in m["passing"])
                out["failed"] = sum(int(x) for x in m["failing"])
                out["skipped"] = sum(int(x) for x in m["pending"])
    if exit_code not in (0, None) and not (out["passed"] or out["failed"]):
        out["errors"] = 1
    return out


def status_of(summary: dict, exit_code: Optional[int]) -> str:
    """ok / failed / nothing_ran / error, the run_tests vocabulary.
    0 passed with nothing failed is nothing_ran, never ok."""
    ran = summary.get("passed", 0) + summary.get("failed", 0) + \
        summary.get("errors", 0)
    if summary.get("failed") or (exit_code not in (0, None) and ran):
        return "failed"
    if summary.get("errors"):
        return "error"
    if not summary.get("passed"):
        return "nothing_ran"
    return "ok" if exit_code in (0, None) else "failed"


__all__ = ["Group", "plan", "summarise", "status_of"]
