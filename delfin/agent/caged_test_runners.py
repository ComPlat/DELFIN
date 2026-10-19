"""The stack test runners the process cage may run without a question.

Input: a shell command, the directory it runs in and the workspace root.
Output: ``admitted_runners(cmd, cwd, workspace)`` -> the runner segments
of the command that are exactly a verification run of the stack found in
that directory, or [] when any part of the command is something else;
``without_admitted_runners(...)`` -> the command with those segments
replaced by ``true``, for the rest of the cage rule to judge.

Consulted only by ``api_client._caged_command_may_run_unasked``, i.e. only
when bubblewrap, the process cage, Landlock read confinement and a
non-open egress policy all hold. Outside the cage nothing here applies
and the runner is asked, as before.

Admitted, when the marker of the stack is in the run's directory or one
above it inside the workspace:

* ``cargo test|check|clippy`` (Cargo.toml), flags from ``_CARGO_FLAGS``;
  anything after ``--`` goes to the test binary;
* ``go test|vet`` (go.mod), flags from ``_GO_FLAGS``, package patterns;
* ``npm test``, ``npm t``, ``npm run test`` (package.json with a
  ``scripts.test``), ``--silent``; anything after ``--``;
* ``npx vitest [run]``, ``npx tsc --noEmit`` -- only when the tool is
  installed in the project (``node_modules/.bin/<tool>``), so npx runs the
  local copy and fetches nothing;
* ``node --test`` needs nothing here: ``node`` is not an excluded word;
* an install INTO THE PROJECT (since 2026-10-19): ``npm install|i|ci|add``,
  ``npm init -y``, ``yarn [install|add]``, ``pnpm install|i|add`` with
  package specs and flags from ``_JS_INSTALL_FLAGS``; ``uv venv``, ``uv pip
  install|sync``, ``uv add|sync|lock`` with ``_UV_FLAGS``. It writes
  node_modules / .venv in the workspace and fetches from the package
  indexes the egress proxy knows -- the cage's own bounds. Measured: a
  browser game's agent could neither install nor fetch three.js and
  rebuilt a 3D renderer by hand.

Never admitted: ``cargo install|publish|add|run|update|fetch``,
``+toolchain``, ``--config``, ``-Z``; ``go get|install|mod|run``,
``-exec``, ``-toolexec``, ``-o``, ``-modfile``, ``-mod=mod`` (cargo and
go install into the home directory, outside the project); ``npm
publish|exec|run <other>``; any install with ``-g``/``--global``,
``global``, ``--prefix``, ``--location``, ``--system``, ``--target``, a
path or a URL; ``npx`` with
``-p``/``--package``/``-y``/``--yes`` or any other package; environment
assignments in the runner segment; quotes, command substitution,
backticks, redirects other than ``2>&1`` and a background ``&``.

Security: these runners execute repository code (build scripts, test
code, npm lifecycle scripts) exactly like pytest does. What bounds that
is the cage (writes confined to the roots with protected files
read-only, no visible foreign process, egress through the proxy policy,
reads confined by Landlock), not the program name; the grammar only keeps
the network and install verbs, whose effect leaves the cage's purpose,
on the dialog. Rejected: admitting by stack detection alone without the
grammar (``cargo install`` would pass), and admitting outside the cage
(the rule that security holds on every host).
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Optional

_SPLIT_RE = re.compile(r"\s*(?:&&|\|\||;|\n|\|)\s*")
# Inside one runner segment (``&&``, ``||``, ``;`` and ``|`` already split
# off, ``2>&1`` removed): quoting, substitution, redirects, background.
_UNSAFE_RE = re.compile(r"[\"'`$<>\\(){}&;]")

_CARGO_SUB = {"test", "check", "clippy"}
_CARGO_FLAGS = {
    "--offline", "--locked", "--frozen", "-q", "--quiet", "-v", "--verbose",
    "--release", "--workspace", "--all", "--all-targets", "--all-features",
    "--no-default-features", "--lib", "--bins", "--tests", "--doc",
    "--no-fail-fast", "--no-run", "--benches", "--examples",
}
_CARGO_VALUE_FLAGS = {"--features", "-F", "--package", "-p", "--test",
                      "--bin", "--example", "--bench", "-j", "--jobs",
                      "--color", "--message-format", "--manifest-path"}
_GO_SUB = {"test", "vet"}
_GO_FLAGS = {"-v", "-race", "-short", "-cover", "-failfast", "-json",
             "-benchmem", "-mod=vendor", "-mod=readonly", "-trimpath"}
_GO_VALUE_FLAGS = {"-run", "-count", "-timeout", "-bench", "-benchtime",
                   "-skip", "-parallel", "-cpu", "-tags", "-shuffle",
                   "-covermode"}
_NAME_RE = re.compile(r"^[\w.:/@*,=+-]+$")


def _find_up(start: Path, workspace: Path, name: str) -> Optional[Path]:
    try:
        start, workspace = start.resolve(), workspace.resolve()
    except OSError:
        return None
    cur = start
    while True:
        cand = cur / name
        if cand.exists():
            return cand
        if cur == workspace or workspace not in cur.parents:
            return None
        cur = cur.parent


def _inside(path: Path, workspace: Path) -> bool:
    try:
        p, w = path.resolve(), workspace.resolve()
    except OSError:
        return False
    return p == w or w in p.parents


def _flags_ok(words: list[str], plain: set, valued: set,
              value_ok=None) -> bool:
    """Every word a known flag, a known flag's value, or a plain name;
    after ``--`` everything is a plain argument to the test binary."""
    i = 0
    while i < len(words):
        w = words[i]
        if w == "--":
            return all(_NAME_RE.match(x) for x in words[i + 1:])
        if "=" in w and w.split("=", 1)[0] in valued:
            if value_ok and not value_ok(w.split("=", 1)[0],
                                         w.split("=", 1)[1]):
                return False
            i += 1
            continue
        if w in plain:
            i += 1
            continue
        if w in valued:
            if i + 1 >= len(words) or not _NAME_RE.match(words[i + 1]):
                return False
            if value_ok and not value_ok(w, words[i + 1]):
                return False
            i += 2
            continue
        if w.startswith("-") or not _NAME_RE.match(w):
            return False
        i += 1
    return True


def _cargo(words: list[str], here: Path, ws: Path) -> bool:
    if len(words) < 2 or words[1] not in _CARGO_SUB:
        return False

    def manifest_ok(flag: str, value: str) -> bool:
        if flag != "--manifest-path":
            return True
        return _inside(here / value, ws) and value.endswith("Cargo.toml")

    if not _flags_ok(words[2:], _CARGO_FLAGS, _CARGO_VALUE_FLAGS, manifest_ok):
        return False
    if "--manifest-path" in " ".join(words):
        return True
    return _find_up(here, ws, "Cargo.toml") is not None


def _go(words: list[str], here: Path, ws: Path) -> bool:
    if len(words) < 2 or words[1] not in _GO_SUB:
        return False
    rest = words[2:]
    i = 0
    while i < len(rest):
        w = rest[i]
        if w.startswith("-"):
            name = w.split("=", 1)[0]
            if w in _GO_FLAGS:
                i += 1
            elif name in _GO_VALUE_FLAGS and "=" in w:
                if not _NAME_RE.match(w.split("=", 1)[1]):
                    return False
                i += 1
            elif name in _GO_VALUE_FLAGS:
                if i + 1 >= len(rest) or not _NAME_RE.match(rest[i + 1]):
                    return False
                i += 2
            else:
                return False
        elif _NAME_RE.match(w):
            i += 1
        else:
            return False
    return _find_up(here, ws, "go.mod") is not None


def _package_json(here: Path, ws: Path) -> Optional[dict]:
    p = _find_up(here, ws, "package.json")
    if p is None:
        return None
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _npm(words: list[str], here: Path, ws: Path) -> bool:
    rest = words[1:]
    if rest[:1] in (["test"], ["t"]):
        rest = rest[1:]
    elif rest[:2] == ["run", "test"]:
        rest = rest[2:]
    else:
        return False
    rest = [w for w in rest if w not in ("--silent", "-s")]
    if rest and (rest[0] != "--" or not all(_NAME_RE.match(x)
                                            for x in rest[1:])):
        return False
    pkg = _package_json(here, ws)
    scripts = (pkg or {}).get("scripts") or {}
    return isinstance(scripts, dict) and bool(scripts.get("test"))


def _npx(words: list[str], here: Path, ws: Path) -> bool:
    rest = [w for w in words[1:] if w not in ("--no-install", "--offline")]
    if not rest:
        return False
    tool, args = rest[0], rest[1:]
    if tool == "vitest":
        if args[:1] in (["run"], ["--run"]):
            args = args[1:]
        ok = _flags_ok(args, {"--run", "--reporter=dot", "--silent",
                              "--passWithNoTests"},
                       {"-t", "--testNamePattern", "--project", "--dir",
                        "--root", "--config"})
        if not ok or any(a.startswith("--config") for a in args):
            return False
    elif tool == "tsc":
        if "--noEmit" not in args:
            return False
        if not _flags_ok(args, {"--noEmit", "--pretty", "--strict"},
                         {"-p", "--project"}):
            return False
    else:
        return False
    proj = _find_up(here, ws, "package.json")
    if proj is None:
        return False
    return (proj.parent / "node_modules" / ".bin" / tool).exists()


# A package spec as an install names it: `three`, `three@0.160.1`,
# `@types/node@^20`. Never a flag, a path or a URL (no "-" first, no "/"
# outside a scope, no ":", and not `.` or `..`, which npm and uv install as
# the directory they name).
_PKG_RE = re.compile(
    r"^(?!-)(?!\.+(?:@|$))(?:@[\w.-]+/)?[\w.-]+(?:@[\w.^~<>=*|-]+)?$")
_JS_INSTALL_FLAGS = {
    "-D", "-E", "-O", "-P", "--save", "--save-dev", "--save-exact",
    "--save-optional", "--save-peer", "--no-save", "--no-audit", "--no-fund",
    "--ignore-scripts", "--legacy-peer-deps", "--production",
    "--frozen-lockfile", "--prefer-offline", "--dev", "-W", "--silent",
    "--omit=dev", "--omit=optional", "--omit=peer",
}
_UV_FLAGS = {"--upgrade", "-U", "--dev", "--frozen", "--locked", "--seed",
             "-q", "--quiet", "--no-cache"}
_UV_VALUE_FLAGS = {"--python", "-p", "-r", "--requirement", "-e",
                   "--editable"}


def _specs_ok(words: list[str], flags: set) -> bool:
    # `yarn global add`: an install outside the project, spelled as a word.
    return all((w in flags or _PKG_RE.match(w)) and w != "global"
               for w in words)


def _npm_install(words: list[str], here: Path, ws: Path) -> bool:
    if words[1:] in (["init", "-y"], ["init", "--yes"]):
        return True
    if len(words) < 2 or words[1] not in ("install", "i", "ci", "add"):
        return False
    return _specs_ok(words[2:], _JS_INSTALL_FLAGS)


def _yarn(words: list[str], here: Path, ws: Path) -> bool:
    rest = words[1:]
    if rest[:1] in (["install"], ["add"]):
        rest = rest[1:]
    if not _specs_ok(rest, _JS_INSTALL_FLAGS):
        return False
    return _find_up(here, ws, "package.json") is not None


def _pnpm(words: list[str], here: Path, ws: Path) -> bool:
    if len(words) < 2 or words[1] not in ("install", "i", "add"):
        return False
    if not _specs_ok(words[2:], _JS_INSTALL_FLAGS):
        return False
    return _find_up(here, ws, "package.json") is not None


def _uv(words: list[str], here: Path, ws: Path) -> bool:
    rest = words[1:]
    if rest[:1] == ["venv"]:
        rest = rest[1:]
        names = [w for w in rest if not w.startswith("-")]
        if any(not re.match(r"^\.?venv[\w.-]*$", n) and not re.match(
                r"^[\d.]+$", n) for n in names):
            return False
    elif rest[:2] in (["pip", "install"], ["pip", "sync"]):
        rest = rest[2:]
    elif rest[:1] in (["add"], ["sync"], ["lock"]):
        rest = rest[1:]
    else:
        return False

    def value_ok(flag: str, value: str) -> bool:
        if flag in ("-e", "--editable"):
            return value == "."
        if flag in ("-r", "--requirement"):
            return _inside(here / value, ws) and value.endswith(".txt")
        return bool(re.match(r"^[\d.]+$", value))

    # A path (a requirements file of `uv pip sync`, a local package) only
    # inside the workspace: `uv pip install ../other` or `/abs/pkg` builds
    # and installs what lies outside it.
    return _flags_ok(rest, _UV_FLAGS, _UV_VALUE_FLAGS, value_ok) and all(
        w.startswith("-") or _PKG_RE.match(w)
        or (re.match(r"^[\w./-]+$", w) and _inside(here / w, ws))
        for w in rest)


def _npm_any(words: list[str], here: Path, ws: Path) -> bool:
    return _npm(words, here, ws) or _npm_install(words, here, ws)


_RUNNERS = {"cargo": _cargo, "go": _go, "npm": _npm_any, "npx": _npx,
            "yarn": _yarn, "pnpm": _pnpm, "uv": _uv}


def is_project_install(cmd: str) -> bool:
    """Whether ``cmd`` is, by its words alone, an install into the project
    (npm/yarn/pnpm/uv; module docstring). No marker file is checked: the
    caller asks whether it leaves the workspace, not whether it may run."""
    words = str(cmd or "").split()
    if not words or _UNSAFE_RE.search(str(cmd)) or re.search(
            r"\bglobal\b", str(cmd)):
        return False
    tmp = Path("/nonexistent-delfin-probe")
    if words[0] == "npm":
        return _npm_install(words, tmp, tmp)
    if words[0] == "yarn":
        rest = words[1:]
        if rest[:1] in (["install"], ["add"]):
            rest = rest[1:]
        return _specs_ok(rest, _JS_INSTALL_FLAGS)
    if words[0] == "pnpm":
        return (len(words) > 1 and words[1] in ("install", "i", "add")
                and _specs_ok(words[2:], _JS_INSTALL_FLAGS))
    if words[0] == "uv":
        return _uv(words, tmp, tmp)
    return False


def _segments(cmd: str) -> Optional[list[str]]:
    """The command's segments, or None when it has a background ``&``."""
    text = str(cmd or "").replace("2>&1", " ")
    if re.search(r"(?<![&|])&(?![&>])", text):
        return None
    return [s for s in _SPLIT_RE.split(text) if s.strip()]


def admitted_runners(cmd: str, cwd: str, workspace: str) -> list[str]:
    """The segments of ``cmd`` that are an admitted stack runner (module
    docstring); a ``cd`` before a runner moves its directory."""
    ws = Path(workspace or ".")
    here = Path(cwd) if cwd and os.path.isabs(cwd) else ws / (cwd or "")
    segs = _segments(cmd)
    if not segs:
        return []
    out: list[str] = []
    for seg in segs:
        words = seg.split()
        if not words:
            continue
        if words[0] == "cd" and len(words) == 2 and not _UNSAFE_RE.search(seg):
            nxt = here / words[1]
            if not _inside(nxt, ws):
                return []
            here = nxt
            continue
        fn = _RUNNERS.get(words[0])
        if fn is None:
            continue
        if _UNSAFE_RE.search(seg) or not fn(words, here, ws):
            return []
        out.append(seg)
    return out


def without_admitted_runners(cmd: str, cwd: str, workspace: str) -> Optional[str]:
    """``cmd`` with each admitted runner segment replaced by ``true``;
    None when it holds none."""
    segs = admitted_runners(cmd, cwd, workspace)
    if not segs:
        return None
    text = str(cmd)
    for seg in segs:
        text = text.replace(seg, "true", 1)
    return text


__all__ = ["admitted_runners", "without_admitted_runners",
           "is_project_install"]
