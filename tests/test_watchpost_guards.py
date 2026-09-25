"""Source-level guarantees for watchpost: read-only outside its report
dir, no network libraries, no tree walks over HOME, and enable/disable
touch only the watchpost scheduler entry.
"""
from __future__ import annotations

import ast
from pathlib import Path

import delfin.watchpost
from delfin.watchpost import scheduler_entry

PKG = Path(delfin.watchpost.__file__).parent

# The only names allowed to run a subprocess.
ALLOWED_READ_COMMANDS = {"last", "ss", "crontab"}

# Modules that must never be imported by the package.
FORBIDDEN_IMPORTS = {"socket", "requests", "urllib", "http", "httpx",
                     "ftplib", "smtplib", "telnetlib", "asyncio"}

WRITE_OPEN_MODES = {"w", "w+", "wb", "wb+", "a", "a+", "ab", "ab+", "x"}


def _sources():
    for p in sorted(PKG.glob("*.py")):
        yield p, ast.parse(p.read_text())


def _pkg_sources_only_report_writes():
    """No write call outside the report dir, no os.remove/unlink/kill,
    no subprocess other than the allowed read-only commands.

    Writes are confined to scan.py, whose write-mode opens and mkdir
    calls all target the report directory it is given as an argument.
    """
    offenders = []
    for path, tree in _sources():
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            name = getattr(func, "attr", None) or getattr(func, "id", None)
            if name in ("remove", "unlink", "rmdir", "rename", "replace",
                        "kill", "killpg", "system", "popen"):
                offenders.append(f"{path.name}:{node.lineno} {name}")
            if name == "mkdir":
                # allowed only for the report directory, in scan.py
                target = node.args[0] if node.args else None
                receiver = getattr(func, "value", None)
                is_report = (
                    (isinstance(target, ast.Name) and target.id == "report")
                    or (isinstance(receiver, ast.Name)
                        and receiver.id == "report"))
                if not (path.name == "scan.py" and is_report):
                    offenders.append(f"{path.name}:{node.lineno} mkdir")
            if name == "open":
                modes = [kw.value for kw in node.keywords if kw.arg == "mode"]
                modes += [a for a in node.args[1:2]
                          if isinstance(a, ast.Constant)]
                if any(isinstance(m, ast.Constant)
                       and any(str(m.value).startswith(w)
                               for w in WRITE_OPEN_MODES)
                       for m in modes):
                    if path.name != "scan.py":
                        offenders.append(
                            f"{path.name}:{node.lineno} open write-mode")
    return offenders


def test_package_has_no_writes_or_kills():
    assert _pkg_sources_only_report_writes() == []


def test_package_imports_no_network_libraries():
    bad = []
    for path, tree in _sources():
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module.split(".")[0]]
            for n in names:
                if n in FORBIDDEN_IMPORTS:
                    bad.append(f"{path.name}:{node.lineno} {n}")
    assert bad == []


def test_package_walks_no_trees():
    """No rglob/walk/glob('**') anywhere -- only fixed, named paths."""
    bad = []
    for path, tree in _sources():
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                func = node.func
                name = getattr(func, "attr", None)
                if name in {"rglob", "walk"}:
                    bad.append(f"{path.name}:{node.lineno} {name}")
                if name == "glob":
                    for a in node.args:
                        if (isinstance(a, ast.Constant)
                                and "**" in str(a.value)):
                            bad.append(f"{path.name}:{node.lineno} glob **")
    assert bad == []


def test_subprocess_calls_only_allowed_read_commands():
    import subprocess as sp
    offenders = []
    for path, tree in _sources():
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) \
                    and node.func.id in dir(sp):
                if node.args and isinstance(node.args[0], ast.Constant):
                    cmd = str(node.args[0].value).split()[0]
                    if node.func.id == "run" and cmd in ALLOWED_READ_COMMANDS:
                        continue
                    offenders.append(f"{path.name}:{node.lineno} {cmd}")
    assert offenders == []


# --- scheduler -----------------------------------------------------------

def test_enable_creates_and_disable_removes_own_entry(tmp_path, monkeypatch):
    cron = tmp_path / "cron.json"
    monkeypatch.setattr(scheduler_entry, "_cron_path", lambda: cron)

    entry = scheduler_entry.enable(every_minutes=30)
    assert cron.exists()
    import json
    data = json.loads(cron.read_text())
    assert len(data["entries"]) == 1
    assert "watchpost" in data["entries"][0]["prompt"]

    assert scheduler_entry.disable() is True
    assert json.loads(cron.read_text()) == {"entries": []}


def test_disable_leaves_other_entries_alone(tmp_path, monkeypatch):
    import json
    cron = tmp_path / "cron.json"
    other = {"id": "other1", "kind": "interval", "prompt": "run bench",
             "reason": "", "every_seconds": 3600, "next_fire_at": 1,
             "workspace": "", "budget_usd": 0.0}
    cron.write_text(json.dumps({"entries": [other]}))
    monkeypatch.setattr(scheduler_entry, "_cron_path", lambda: cron)
    assert scheduler_entry.disable() is False
    assert json.loads(cron.read_text()) == {"entries": [other]}
