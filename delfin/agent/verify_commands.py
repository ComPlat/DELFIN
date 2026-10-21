"""Recognising a verification run in a tool call and reading its result.

One answer to three questions that the benchmark scorer, the turn-end
verification check and the test-claim scanner all ask:

* ``command_kind(cmd)``: is this shell command a test run (``"test"``),
  a static check (``"check"``) or neither (``None``)? Stack-neutral:
  pytest/unittest/tox/nox, cargo test/nextest, go test, npm/yarn/pnpm/bun
  test, npx vitest/jest/mocha, node --test, deno test, make test/check,
  ctest; checks are cargo check/clippy/build, go vet/build, tsc --noEmit,
  mypy, ruff, pyflakes, flake8, py_compile, eslint.
* ``bash_writes_file(cmd)``: does the command write a file (``sed -i``,
  a redirect into a path, ``tee``, a here-document into a file,
  ``patch``/``git apply``)? Counted as an edit, because the edited state
  is what a later test has to cover.
* ``run_outcome(tool, result)``: did the call execute, and how did it
  end? ``refused`` (a gate said no), ``error`` (the tool failed before
  running anything), ``started`` (a background job, outcome not yet
  seen), ``passed``, ``failed``, ``nothing_ran`` (0 collected / all
  skipped), ``unknown`` (no result to read).

The runner must be the command, not an argument: ``echo pytest`` is not a
test run. A command is matched at its start or after ``;``, ``&&``,
``||``, ``|`` or ``(``, after optional ``VAR=value`` assignments and a
``timeout N`` prefix.

Rejected: a separate pattern per consumer. Four existed
(benchmark._TEST_CMD_RE, api_client._TEST_BASH_CMD_RE,
verify_guard._TEST_CMD_RE, counter_test._TEST_CMD_RE) and disagreed: the
benchmark counted ``cargo test`` as a test, the evidence ledger did not,
and neither knew vitest. Pure functions, no I/O, never raise.
"""

from __future__ import annotations

import json
import re
from typing import Any, Optional

# Where a command may start inside a compound command line.
_HEAD = (r"(?:^|[;&|(\n]|&&|\|\|)\s*"
         r"(?:[A-Za-z_][A-Za-z0-9_]*=\S*\s+)*"
         r"(?:timeout\s+(?:-\S+\s+)*\S+\s+)?")

_TEST_BODY = (
    r"(?:python[\d.]*\s+-m\s+(?:pytest|unittest|tox|nox)\b"
    r"|(?:\S*/)?(?:pytest|py\.test)\b"
    r"|(?:uv\s+run|poetry\s+run|pipenv\s+run|hatch\s+run)\s+(?:python[\d.]*\s+-m\s+)?pytest\b"
    r"|tox\b|nox\b"
    r"|cargo\s+(?:\+\S+\s+)?(?:test|nextest)\b"
    r"|go\s+test\b"
    r"|(?:npm|pnpm|yarn|bun)\s+(?:run\s+)?test(?::[\w:-]+)?\b"
    r"|npm\s+t\b"
    r"|(?:npx|pnpm\s+exec|yarn\s+exec|bunx)\s+(?:--\S+\s+)*(?:vitest|jest|mocha|ava|tap)\b"
    r"|(?:\S*/)?(?:vitest|jest|mocha)\b"
    r"|node\s+(?:--[\w-]+(?:=\S+)?\s+)*--test\b"
    r"|deno\s+test\b|bun\s+test\b"
    r"|make\s+(?:-\S+\s+)*(?:test|check|tests)\b"
    r"|ctest\b"
    r")"
)

_CHECK_BODY = (
    r"(?:cargo\s+(?:\+\S+\s+)?(?:check|clippy|build)\b"
    r"|go\s+(?:vet|build)\b"
    r"|(?:npx\s+|pnpm\s+exec\s+|yarn\s+exec\s+)?(?:\S*/)?tsc\b"
    r"|python[\d.]*\s+-m\s+(?:mypy|py_compile|pyflakes|flake8|ruff|compileall)\b"
    r"|(?:mypy|pyflakes|flake8|pyright)\b"
    r"|ruff\s+check\b"
    r"|(?:npx\s+)?eslint\b"
    r")"
)

TEST_COMMAND_RE = re.compile(_HEAD + _TEST_BODY)
CHECK_COMMAND_RE = re.compile(_HEAD + _CHECK_BODY)


def command_kind(cmd: Any) -> Optional[str]:
    """``"test"``, ``"check"`` or None for a shell command string."""
    text = str(cmd or "")
    if not text.strip():
        return None
    if TEST_COMMAND_RE.search(text):
        return "test"
    if CHECK_COMMAND_RE.search(text):
        return "check"
    return None


# A shell command that writes a file. A redirect into /dev/null or into a
# descriptor (2>&1) writes nothing; a redirect into a path does.
_BASH_WRITE_RE = re.compile(
    r"\bsed\s+(?:-\w*\s+)*-i|\bperl\s+(?:-\w*\s+)*-i"
    r"|(?<![0-9&<>])>>?\s*(?!&|/dev/null\b|\s*$)['\"]?[\w./~-]"
    r"|(?:^|[;&|(]|&&)\s*tee\b|\|\s*tee\b"
    r"|\bpatch\s+-p|\bgit\s+apply\b"
)


_HEREDOC_RE = re.compile(r"<<-?\s*['\"]?(\w+)['\"]?")


def shell_text(cmd: Any) -> str:
    """The command with quoted strings and here-document bodies blanked.

    What the shell itself acts on: ``print('1,5 ->')`` inside a quoted
    argument redirects nothing, and a here-document body is data. A
    ``cat > f << 'EOF'`` keeps its first line, which holds the redirect."""
    text = str(cmd or "")
    lines = text.split("\n")
    out: list[str] = []
    end: Optional[str] = None
    for ln in lines:
        if end is not None:
            if ln.strip() == end:
                end = None
            continue
        out.append(ln)
        m = _HEREDOC_RE.search(ln)
        if m:
            end = m.group(1)
    text = "\n".join(out)
    text = re.sub(r"'[^'\n]*'", "''", text)
    return re.sub(r'"(?:\\.|[^"\\\n])*"', '""', text)


def bash_writes_file(cmd: Any) -> bool:
    """True when the shell command writes a file (see module docstring).
    Judged on ``shell_text``: a ``>`` inside a quoted string or a
    here-document body is not a redirect."""
    return bool(_BASH_WRITE_RE.search(shell_text(cmd)))


# The gate's refusals and asks. Each phrase is the head of a real refusal
# text in api_client (allow list, plan mode, read scope, approval).
_REFUSAL_MARKERS = (
    "not on the auto-allow list", "blocked:", "plan mode", "refused",
    "denied", "needs approval", "needs your approval", "not allowed",
    "command rejected", "approval", "outside the allowed",
)

_FAILED_TEXT_RE = re.compile(
    r"\b[1-9]\d*\s+(?:failed|errors?)\b|^FAILED\b|^FAIL\b|test result: FAILED"
    r"|\bpanicked at\b|^--- FAIL:|^not ok \d|ℹ fail [1-9]|\bTests?:\s+[1-9]\d*\s+failed"
    r"|^error(?:\[E\d+\])?:|AssertionError|Traceback \(most recent call last\)",
    re.M)
_NOTHING_RAN_RE = re.compile(
    r"\bno tests ran\b|\bcollected 0 items\b|\brunning 0 tests\b"
    r"|\bno test files\b|ℹ tests 0\b|\bNo tests found\b|\b0 passed\b"
    r"|\bRan 0 tests\b|\[no test files\]", re.I)
_PASSED_COUNT_RE = re.compile(
    r"\b([1-9]\d*)\s+passed\b|test result: ok\. ([1-9]\d*) passed"
    r"|ℹ pass ([1-9]\d*)|\bRan ([1-9]\d*) tests?\b|^ok\s+\S+", re.M)


_TRUNC_FIELDS_RE = re.compile(
    r'"(exit_code|status|framework|job_id|timed_out|passed|failed|errors|'
    r'skipped)"\s*:\s*("(?:[^"\\]|\\.)*"|-?\d+|true|false)')


def _as_json(result: Any) -> Optional[dict]:
    """A result's JSON object; for a head slice that no longer parses (the
    engine hands callers 2,000 characters), the scalar fields that made it
    into the slice, with ``stdout`` the raw text."""
    if isinstance(result, dict):
        return result
    text = str(result or "").strip()
    if not text.startswith("{"):
        return None
    try:
        data = json.loads(text)
    except (TypeError, ValueError):
        if text.startswith('{"error"'):
            return {"error": text[10:600]}
        fields: dict = {}
        summary: dict = {}
        for key, val in _TRUNC_FIELDS_RE.findall(text[:6000]):
            v: Any = (json.loads(val) if val[:1] in '"tf' or val.lstrip("-").isdigit()
                      else val)
            if key in ("passed", "failed", "errors", "skipped"):
                summary.setdefault(key, v)
            else:
                fields.setdefault(key, v)
        if not fields:
            return None
        if summary:
            fields["summary"] = summary
        fields["stdout"] = text
        return fields
    return data if isinstance(data, dict) else None


def _is_refusal(text: str) -> bool:
    low = text.lower()
    return any(m in low for m in _REFUSAL_MARKERS)


def run_outcome(tool: str, result: Any) -> str:
    """How a bash / bash_background / run_tests call ended (module doc).

    ``tool`` is the bare tool name. Output-based: an exit code 0 whose
    output reports a failure (``cargo test | head`` exits with head's
    status) is ``failed``; an exit code 0 with nothing collected is
    ``nothing_ran``.
    """
    text = str(result or "")
    if not text.strip():
        return "unknown"
    data = _as_json(result)
    if data is None:
        # Plain text: a refusal some backends return unwrapped, or a
        # result the caller could not pair with its call.
        if text.lstrip().lower().startswith(("error", "bash: ")) and \
                _is_refusal(text[:400]):
            return "refused"
        return "unknown"
    if "error" in data and "exit_code" not in data:
        return "refused" if _is_refusal(str(data.get("error"))[:600]) \
            else "error"
    if tool == "run_tests" or "framework" in data:
        status = str(data.get("status") or "").lower()
        summary = data.get("summary") or {}
        try:
            passed = int(summary.get("passed", 0) or 0)
            failed = int(summary.get("failed", 0) or 0) + int(
                summary.get("errors", 0) or 0)
        except (TypeError, ValueError, AttributeError):
            passed, failed = 0, 0
        if status in ("refused", "needs_approval"):
            return "refused"
        if status == "nothing_ran" or (status == "ok" and passed + failed == 0):
            return "nothing_ran"
        if failed or status in ("failed", "fail"):
            return "failed"
        if status == "ok":
            return "passed"
        if status in ("error", "timeout", "timed_out"):
            return "failed" if "exit_code" in data else "error"
        return "unknown"
    if str(data.get("status") or "") == "started" and "job_id" in data:
        return "started"
    if "exit_code" in data:
        out = (str(data.get("stdout") or "") + "\n"
               + str(data.get("stderr") or ""))
        try:
            code = int(data.get("exit_code"))
        except (TypeError, ValueError):
            code = None
        if data.get("timed_out") or code is None:
            return "failed"
        if code != 0 or _FAILED_TEXT_RE.search(out):
            return "failed"
        if _NOTHING_RAN_RE.search(out) and not _PASSED_COUNT_RE.search(out):
            return "nothing_ran"
        return "passed"
    return "unknown"


EXECUTED = frozenset({"passed", "failed", "nothing_ran", "started"})

_COUNT_PATTERNS = (
    re.compile(r"\b(\d+)\s+passed\b"),
    re.compile(r"test result: \w+\. (\d+) passed"),
    re.compile(r"ℹ pass (\d+)"),
    re.compile(r"\bTests:\s+(\d+) passed"),
    re.compile(r"\bRan (\d+) tests?\b"),
)


def passed_count(tool: str, result: Any) -> Optional[int]:
    """The number of passed tests a run reported; None when it reported
    none (go test prints no count). The largest figure wins: a summary
    line sums the per-file lines above it."""
    data = _as_json(result)
    if data is not None and (tool == "run_tests" or "framework" in data):
        try:
            return int((data.get("summary") or {}).get("passed"))
        except (TypeError, ValueError):
            return None
    text = str(result or "")
    if data is not None:
        text = str(data.get("stdout") or "") + "\n" + str(data.get("stderr") or "")
    found = [int(m) for rx in _COUNT_PATTERNS for m in rx.findall(text)]
    return max(found) if found else None


# The measure: files that decide whether a test run is green. Test files
# in each stack's naming, test directories, runner configuration and the
# benchmark's own acceptance scripts. An edit here after a red run can turn
# the run green without changing the code under test.
_MEASURE_NAME_RE = re.compile(
    r"(?:^|/)(?:test_[^/]+\.py|[^/]+_test\.(?:py|go)|conftest\.py"
    r"|[^/]+\.(?:test|spec)\.[cm]?[jt]sx?"
    r"|pytest\.ini|tox\.ini|noxfile\.py|setup\.cfg"
    r"|(?:vitest|jest|mocha|playwright)\.config\.[cm]?[jt]s"
    r"|\.mocharc(?:\.\w+)?)$")
_MEASURE_DIR_RE = re.compile(
    r"(?:^|/)(?:tests?|__tests__|spec|testdata)/"
    r"|(?:^|/)pack/benchmark/(?:accept|setup)/")


# Lines that switch a test off, in each stack's spelling.
_SKIP_MARKER_RE = re.compile(
    r"pytest\.mark\.(?:skip|xfail)|@unittest\.skip|\bskipTest\b|pytest\.skip\("
    r"|#\[ignore\]|\bt\.Skip(?:Now|f)?\(|\b(?:it|test|describe)\.(?:skip|todo|only)\b"
    r"|\bx(?:it|describe|test)\(|@Disabled\b|@Ignore\b")


def measure_edit_is_additive(old: Any, new: Any) -> bool:
    """True when an edit to a test file only ADDS: the old text survives
    unchanged inside the new one, and no skip/xfail/ignore/only marker is
    among the additions. Such an edit (a regression test added beside the
    red one) cannot have turned the red test green.

    Interim, text-level. The content-level question -- did a change weaken
    the measure (a loosened tolerance, a changed expected value, a deleted
    case) -- belongs to the general measure-edit detection; the interface
    this module expects from it is ``measure_edit_weakens(path, before,
    after) -> Optional[str]`` (a reason, or None), which would replace this
    function at its two call sites."""
    o, n = str(old or ""), str(new or "")
    if not o or o not in n:
        return False
    added = n.replace(o, "", 1)
    return not _SKIP_MARKER_RE.search(added)


def is_measure_path(path: Any) -> bool:
    """True when ``path`` is part of what decides a test verdict.

    Path-based on purpose: the caller asks "was the measure touched", and
    whether the change weakened it (a changed assertion, an added skip, a
    loosened tolerance) is the content-level question the general
    measure-edit detection answers. ``pyproject.toml`` and
    ``package.json`` are not listed: they carry the build as well, and an
    edit to them is far more often a dependency than a test switch.
    """
    p = str(path or "").replace("\\", "/").strip()
    if not p:
        return False
    return bool(_MEASURE_NAME_RE.search(p) or _MEASURE_DIR_RE.search(p))


def result_belongs_to(command: str, result: Any) -> bool:
    """Whether ``result`` is the answer to the bash call ``command``.

    The executor echoes the command in its JSON (``"command": ...``) and
    the allow-list refusal quotes its first 120 characters, so the pair
    can be recovered when the result list is not index-aligned with the
    call list (a call whose arguments were invalid JSON reaches the
    results but not the calls)."""
    cmd = str(command or "").strip()
    if not cmd:
        return False
    data = _as_json(result)
    if data is not None and isinstance(data.get("command"), str):
        return data["command"].strip() == cmd
    text = str(result or "")
    return f"'{cmd[:60]}" in text


# ---------------------------------------------------------------------------
# Relevance: did the run touch the edited file?
# ---------------------------------------------------------------------------

_SEG_SPLIT_RE = re.compile(r"\s*(?:&&|\|\||;|\n|\|)\s*")
_PATHLIKE_RE = re.compile(
    r"(?:/|^\.{1,2}$|\.(?:py|pyi|ts|tsx|js|jsx|mjs|cjs|go|rs|rb|jl|r|sh)$"
    r"|\.\.\.$)")
# Flags whose next word is a value, not a target (a filter, a marker, a
# package name).
_VALUE_FLAGS = frozenset({
    "-k", "-m", "-p", "-c", "-o", "-n", "-x", "--maxfail", "-run", "-bench",
    "-count", "-timeout", "-tags", "--features", "--package", "--test",
    "--bin", "--example", "--target", "--config", "-t", "--testNamePattern",
    "--grep", "-g", "--reporter", "--test-name-pattern", "--project",
    "--rootdir", "--durations", "-W", "-e", "--deselect", "--ignore",
})
_TEST_DIR_NAMES = frozenset({"tests", "test", "__tests__", "spec", "testdata"})


def _norm(path: str) -> str:
    p = str(path or "").replace("\\", "/").strip().strip("'\"")
    parts: list[str] = []
    for part in p.split("/"):
        if part in ("", "."):
            continue
        if part == "..":
            if parts:
                parts.pop()
            continue
        parts.append(part)
    return "/".join(parts)


def _join(base: str, rel: str) -> str:
    rel = str(rel or "")
    if rel.startswith("/"):
        return _norm(rel)
    return _norm(f"{base}/{rel}" if base else rel)


def _under(path: str, d: str) -> bool:
    return not d or path == d or path.startswith(d + "/")


def run_scope(command: Any, cwd: str = "",
              to_rel: Optional[Any] = None) -> tuple[str, list[str]]:
    """Where a test, check or program run ran, and what it named.

    Input: the shell command and its ``cwd`` (relative to the workspace,
    "" for the workspace itself). Output: ``(directory, targets)``: the
    directory after any ``cd`` segments before the runner, and the
    path-like words the runner was given (files, directories, globs,
    ``./...``), resolved against that directory. Flags and their values
    (``-k expr``, ``-run Name``) are not targets. No runner segment:
    ``(cwd, [])``. ``to_rel(abs_path) -> str`` maps an absolute path into
    the workspace (the caller knows where it is)."""
    def _j(base: str, rel: str) -> str:
        if str(rel).startswith("/") and to_rel is not None:
            try:
                return _norm(to_rel(rel))
            except Exception:
                return _norm(rel)
        return _join(base, rel)

    d = _norm(to_rel(cwd)) if (to_rel is not None and str(cwd).startswith("/")) \
        else _norm(cwd)
    for seg in _SEG_SPLIT_RE.split(str(command or "")):
        words = seg.split()
        if not words:
            continue
        if words[0] == "cd" and len(words) > 1:
            d = _j(d, words[1])
            continue
        if not (command_kind(seg) or re.match(
                r"(?:\S*/)?(?:python[\d.]*|node|deno|bun|bash|sh|ruby|perl"
                r"|Rscript|julia)\b|\./", seg.strip())):
            continue
        targets: list[str] = []
        skip = False
        for w in words[1:]:
            if skip:
                skip = False
                continue
            if w in _VALUE_FLAGS:
                skip = True
                continue
            if w.startswith("--manifest-path="):
                targets.append(_j(d, w.split("=", 1)[1].rsplit("/", 1)[0]
                                     if "/" in w else "."))
                continue
            if w.startswith("-") or "=" in w:
                continue
            if _PATHLIKE_RE.search(w):
                t = w.split("::", 1)[0]
                if t.endswith("..."):
                    t = t[:-3]
                targets.append(_j(d, t))
        return d, targets
    return d, []


def _stem(path: str) -> str:
    base = path.rsplit("/", 1)[-1]
    return base.split(".", 1)[0]


def run_touches(edited: str, directory: str, targets: list[str],
                read: Optional[Any] = None) -> bool:
    """Whether a run scoped ``(directory, targets)`` exercised ``edited``.

    All paths relative to the workspace. True when

    * the run named no target and ran in a directory containing the file
      (the whole suite of that project);
    * a target is the file, or a directory or glob containing it;
    * a target is a file whose name carries the edited module's stem
      (``test_words.py`` for ``words.py``, ``port_test.go`` for
      ``port.go``, ``split.test.ts`` for ``split.ts``; stems of 3+ chars);
    * a target is a test directory, or a glob in one, whose parent
      contains the file (that project's suite; a single named test file in
      it is not the suite);
    * ``read(path) -> str`` is given and a target file imports or names the
      edited module's stem.

    Rejected: any test anywhere in the repository. That was satisfiable by
    running an unrelated suite after the edit, which is the cheapest way
    to make "tested after the last edit" true."""
    e = _norm(edited)
    if not e:
        return False
    if not targets:
        return _under(e, _norm(directory))
    stem = _stem(e)
    for t in targets:
        tn = _norm(t.split("*", 1)[0].rstrip("/")) if "*" in t else _norm(t)
        if tn == e or _under(e, tn):
            return True
        name = t.rsplit("/", 1)[-1]
        if len(stem) >= 3 and "." in name and stem in name:
            return True
        # A test directory (or a glob in one) is that project's suite; a
        # single named test file in it is not.
        if "*" in t or "." not in name:
            parts = tn.split("/")
            for i, part in enumerate(parts):
                if part in _TEST_DIR_NAMES and _under(e, "/".join(parts[:i])):
                    return True
        if read is not None and len(stem) >= 3:
            try:
                body = str(read(tn) or "")
            except Exception:
                body = ""
            if re.search(r"\b" + re.escape(stem) + r"\b", body):
                return True
    return False


__all__ = [
    "TEST_COMMAND_RE", "CHECK_COMMAND_RE", "EXECUTED", "command_kind",
    "bash_writes_file", "run_outcome", "result_belongs_to", "is_measure_path",
    "run_scope", "run_touches", "passed_count", "shell_text",
    "measure_edit_is_additive",
]
