"""Structured test-runner wrapper.

The agent currently runs ``pytest`` via the ``bash`` tool and parses
human-readable output, which is fragile (cosmetic changes in pytest
break parsing). ``run_tests()`` shells out to pytest with the
``--json-report`` plugin if installed, falls back to ``--junitxml``
otherwise, and returns a stable JSON shape::

    {
        "framework": "pytest",
        "summary": {"passed": 12, "failed": 1, "errors": 0,
                    "skipped": 0, "duration_s": 4.21},
        "failures": [
            {"node_id": "tests/test_x.py::test_y",
             "message": "AssertionError: 1 != 2",
             "duration_s": 0.05}
        ],
        "raw_stdout_tail": "..."     // last 30 lines for context
    }

The function NEVER raises — failures (missing pytest, plugin error,
timeout, …) are reported via ``status``. Output is capped so a
runaway test suite doesn't blow up the agent's context.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


_DEFAULT_TIMEOUT_S = 300
_TAIL_LINES = 30


def _has_json_report() -> bool:
    """Detect pytest-json-report installation in the active interpreter."""
    try:
        import importlib.util
        return importlib.util.find_spec("pytest_jsonreport") is not None
    except ImportError:
        return False


def _tail(text: str, n: int = _TAIL_LINES) -> str:
    if not text:
        return ""
    lines = text.splitlines()
    if len(lines) <= n:
        return text
    return "...\n" + "\n".join(lines[-n:])


def _parse_json_report(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"summary": {}, "failures": [], "_parse_error": str(exc)}
    summary = data.get("summary") or {}
    out: dict[str, Any] = {
        "summary": {
            "passed": int(summary.get("passed", 0)),
            "failed": int(summary.get("failed", 0)),
            "errors": int(summary.get("error", 0)),
            "skipped": int(summary.get("skipped", 0)),
            "duration_s": round(float(data.get("duration", 0.0)), 3),
        },
        "failures": [],
    }
    for test in data.get("tests", []) or []:
        outcome = test.get("outcome", "")
        if outcome not in ("failed", "error"):
            continue
        # The failure message is most reliably in the longrepr / call.crash
        crash = (test.get("call") or {}).get("crash") or {}
        message = (
            crash.get("message")
            or (test.get("call") or {}).get("longrepr")
            or test.get("longrepr")
            or "<no message>"
        )
        if isinstance(message, dict):
            message = json.dumps(message)[:600]
        out["failures"].append({
            "node_id": test.get("nodeid", ""),
            "outcome": outcome,
            "message": str(message)[:600],
            "duration_s": round(float(
                (test.get("call") or {}).get("duration", 0.0)
            ), 3),
        })
    return out


def _parse_junitxml(path: Path) -> dict[str, Any]:
    try:
        tree = ET.parse(str(path))
        root = tree.getroot()
    except (OSError, ET.ParseError) as exc:
        return {"summary": {}, "failures": [], "_parse_error": str(exc)}
    # Either <testsuites> wraps multiple, or a lone <testsuite>
    suites = root.findall("testsuite") if root.tag == "testsuites" else [root]
    passed = failed = errors = skipped = 0
    duration = 0.0
    failures: list[dict[str, Any]] = []
    for s in suites:
        try:
            failed += int(s.get("failures", 0))
            errors += int(s.get("errors", 0))
            skipped += int(s.get("skipped", 0))
            tests = int(s.get("tests", 0))
            duration += float(s.get("time", 0.0))
            passed += tests - failed - errors - skipped
        except (TypeError, ValueError):
            continue
        for tc in s.findall("testcase"):
            for tag in ("failure", "error"):
                fail = tc.find(tag)
                if fail is None:
                    continue
                node_id = (
                    f"{tc.get('classname', '')}::{tc.get('name', '')}"
                    if tc.get("classname")
                    else tc.get("name", "")
                )
                msg = (fail.get("message")
                       or (fail.text or "").strip()
                       or "<no message>")
                failures.append({
                    "node_id": node_id,
                    "outcome": tag,
                    "message": msg[:600],
                    "duration_s": round(float(tc.get("time", 0.0)), 3),
                })
    return {
        "summary": {
            "passed": max(0, passed),
            "failed": failed,
            "errors": errors,
            "skipped": skipped,
            "duration_s": round(duration, 3),
        },
        "failures": failures,
    }


def run_tests(
    workspace: Path | str,
    *,
    target: str = "",
    pytest_args: list[str] | None = None,
    timeout_s: int = _DEFAULT_TIMEOUT_S,
    python: str = "",
) -> dict[str, Any]:
    """Run pytest on ``target`` (relative to ``workspace``) and return JSON.

    Parameters
    ----------
    workspace : Path
        Working directory. pytest is invoked with ``cwd=workspace``.
    target : str
        Optional test path / nodeid. Empty = run pytest's default discovery.
    pytest_args : list[str]
        Extra pytest CLI args (``-x``, ``-k pattern``, etc.).
    timeout_s : int
        Wall-clock cap. Hitting it returns ``status='timeout'``.
    python : str
        Python interpreter to use. Default = ``sys.executable``.
    """
    import sys
    workspace = Path(workspace).expanduser().resolve()
    if not workspace.is_dir():
        return {
            "status": "error",
            "framework": "pytest",
            "error": f"workspace not a directory: {workspace}",
        }
    py = python or sys.executable
    if not Path(py).is_file() and shutil.which(py) is None:
        return {
            "status": "error", "framework": "pytest",
            "error": f"python interpreter not found: {py}",
        }

    use_json = _has_json_report()
    tmpdir = Path(tempfile.mkdtemp(prefix="delfin_tests_"))
    report_path = tmpdir / ("report.json" if use_json else "report.xml")
    cmd = [py, "-m", "pytest"]
    if target:
        cmd.append(target)
    if pytest_args:
        cmd.extend(pytest_args)
    if use_json:
        cmd.extend(["--json-report", f"--json-report-file={report_path}"])
    else:
        cmd.extend([f"--junitxml={report_path}"])
    # Quiet output keeps the captured stdout small.
    cmd.append("-q")

    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            cmd, cwd=str(workspace),
            capture_output=True, text=True,
            timeout=max(5, timeout_s),
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "status": "timeout",
            "framework": "pytest",
            "elapsed_s": round(time.monotonic() - t0, 2),
            "raw_stdout_tail": _tail(exc.stdout or ""),
            "raw_stderr_tail": _tail(exc.stderr or ""),
        }
    except FileNotFoundError as exc:
        return {
            "status": "error", "framework": "pytest",
            "error": f"failed to launch pytest: {exc}",
        }
    finally:
        # The temp report file will be cleaned up below; the dir survives
        # until the process exits if Python errors before us.
        pass

    elapsed = round(time.monotonic() - t0, 2)
    parsed: dict[str, Any]
    if report_path.is_file():
        parsed = (_parse_json_report(report_path)
                  if use_json else _parse_junitxml(report_path))
    else:
        parsed = {"summary": {}, "failures": [],
                  "_parse_error": "no report file produced"}
    try:
        shutil.rmtree(tmpdir, ignore_errors=True)
    except OSError:
        pass

    summary = parsed.get("summary") or {}
    status = "ok"
    if proc.returncode == 0:
        status = "ok"
    elif summary.get("failed", 0) or summary.get("errors", 0):
        status = "failed"
    elif "_parse_error" in parsed:
        status = "error"
    else:
        status = "failed"

    return {
        "status": status,
        "framework": "pytest",
        "report_format": "json-report" if use_json else "junitxml",
        "exit_code": proc.returncode,
        "elapsed_s": elapsed,
        "summary": summary,
        "failures": parsed.get("failures", []),
        "raw_stdout_tail": _tail(proc.stdout or ""),
        "raw_stderr_tail": _tail(proc.stderr or ""),
    }


__all__ = ["run_tests"]
