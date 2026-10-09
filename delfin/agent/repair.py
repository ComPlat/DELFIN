"""Repair — put DELFIN's settings, venv and tools back the way they should be.

The doctor report (`delfin.agent.doctor.run_doctor`) already numbers the
prerequisites DELFIN can fix itself: rows carry a prose ``fix`` plus --- where
a fixable remedy exists --- a machine-actionable ``command`` (a shell line to
run) or ``setting`` (a ``(key, value)`` pair to write into the user settings
file). ``repair`` turns those rows into ordered steps and applies them one at
a time, so a person approves each before anything runs.

Contract:

- ``plan(report)`` returns the ordered, fixable steps derived from a doctor
  report. Prose-only rows (a system package, a credential helper) are left
  out: repair never improvises a remedy out of a sentence.
- ``apply(step, approved=...)`` applies exactly one step and REFUSES unless
  ``approved`` is true. Each step records what it changed and how to undo it.
- Settings are never overwritten or deleted: before a ``setting`` step writes,
  the current settings file is MOVED to a dated backup, and the new file is
  written fresh from the backed-up content plus the change. ``undo`` for that
  step names the backup to restore.
- The settings path defaults to the real ``~/.delfin/settings.json`` but is a
  parameter, so a test can point it at a faked home and the real file is
  never read or written in a test.
"""

from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

#: Where the user-global DELFIN settings live; overridable per call for tests.
DEFAULT_SETTINGS_PATH = Path("~/.delfin/settings.json").expanduser()


# ---------------------------------------------------------------------------
# plan()
# ---------------------------------------------------------------------------

def _step(row: dict, seq: int) -> dict:
    """One fixable step out of a doctor row, with its undo stated."""
    command = str(row.get("command") or "").strip()
    setting = row.get("setting")
    if command:
        what = (f"{row['check']}: run '{command}'")
        undo = ("the command has no clean inverse; reverse it by hand if it "
                "made a change (e.g. uninstall what it installed)")
        return {
            "id": seq, "check": str(row.get("check", "")),
            "kind": "command", "what": what, "undo": undo,
            "command": command,
        }
    if (isinstance(setting, (list, tuple)) and len(setting) == 2
            and isinstance(setting[0], str) and setting[0].strip()):
        key, value = str(setting[0]), str(setting[1])
        what = (f"{row['check']}: set {key} = {value!r} in the user settings "
                f"file (the whole file is backed up to a dated name first)")
        undo = (f"restore the prior settings from the dated backup named in "
                f"the apply result")
        return {
            "id": seq, "check": str(row.get("check", "")),
            "kind": "setting", "what": what, "undo": undo,
            "setting_key": key, "setting_value": value,
        }
    return {}  # prose-only: not fixable by repair, caller skips it


def plan(report: list) -> list:
    """Ordered, fixable steps derived from a doctor report.

    Only rows with a ``command`` or ``setting`` remedy become steps. The
    order is the report's own (its stable check order), each step gets a
    stable increasing ``id`` so a caller can present and approve them one
    at a time.
    """
    steps: list[dict] = []
    seq = 0
    for row in report or []:
        if not isinstance(row, dict):
            continue
        step = _step(row, seq)
        if step:
            seq += 1
            steps.append(step)
    return steps


# ---------------------------------------------------------------------------
# apply()
# ---------------------------------------------------------------------------

def _split_nested(key: str) -> list[str]:
    return [part for part in key.split(".") if part]


def _set_nested(root: dict, dotted: str, value: Any) -> None:
    """Set ``dotted`` (e.g. ``agent.mcp_isolation``) into ``root`` in place."""
    parts = _split_nested(dotted)
    node = root
    for part in parts[:-1]:
        child = node.setdefault(part, {})
        if not isinstance(child, dict):
            child = {}
            node[part] = child
        node = child
    node[parts[-1]] = value


def _dated_backup_path(settings: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return settings.with_name(f"{settings.name}.{stamp}.bak")


def _apply_setting(step: dict, settings: Path, on_line: Callable) -> dict:
    """Settings step: move the old file to a dated backup, write the new one.

    The old file is never overwritten or deleted in place: it is renamed to
    a dated sibling, and the fresh file is written from the backed-up
    content with the requested change applied. When no file exists there is
    nothing to back up and a fresh file is written directly.
    """
    if not isinstance(settings, Path):
        settings = Path(settings)
    prior = {}
    backup_path = ""
    if settings.exists():
        try:
            prior = json.loads(settings.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            prior = {}
        if not isinstance(prior, dict):
            prior = {}
        backup_path = str(_dated_backup_path(settings))
        _say(on_line, f"moving {settings} to {backup_path}")
        os.replace(settings, Path(backup_path))
    else:
        settings.parent.mkdir(parents=True, exist_ok=True)

    updated = dict(prior)
    _set_nested(updated, step.get("setting_key", ""),
                step.get("setting_value"))
    _atomic_write_json(settings, updated)
    _say(on_line, f"wrote {settings} with {step.get('setting_key')}")
    what = step.get("what", "")
    if backup_path:
        undo = (f"restore the prior settings by moving the dated backup "
                f"{backup_path} back over {settings}")
    else:
        undo = (f"remove {settings} to restore the absence of a settings file")
    return {
        "ok": True, "kind": "setting", "check": step.get("check", ""),
        "what": what, "backup_path": backup_path,
        "settings_path": str(settings),
        "undo": undo,
    }


def _apply_command(step: dict, on_line: Callable) -> dict:
    """Command step: run the declared shell line and report its result."""
    command = step.get("command", "")
    _say(on_line, f"running: {command}")
    try:
        proc = subprocess.run(command, shell=True, text=True,
                              capture_output=True)
        status = "ok" if proc.returncode == 0 else f"exit {proc.returncode}"
        return {
            "ok": proc.returncode == 0, "kind": "command",
            "check": step.get("check", ""), "what": step.get("what", ""),
            "status": status,
            "undo": step.get("undo", ""),
            "output": (proc.stdout or "") + (proc.stderr or ""),
        }
    except Exception as exc:  # noqa: BLE001 — report, never raise
        return {"ok": False, "kind": "command", "check": step.get("check", ""),
                "what": step.get("what", ""), "status": f"error: {exc}",
                "undo": step.get("undo", "")}


def _atomic_write_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def _say(on_line, line: str) -> None:
    if on_line is not None:
        on_line(line)


def apply(step: dict, *, approved: bool = False,
          user_settings_path: str | Path | None = None,
          on_line: Callable | None = None) -> dict:
    """Apply one repair step; refuses unless explicitly approved.

    ``approved`` must be true for anything to run -- the whole point of the
    repair is a person in the loop, one approval per step.

    For a ``setting`` step the settings file is never overwritten or
    deleted: the old file is moved to a dated backup first, then the new
    file is written fresh from that content. ``user_settings_path`` points
    at the settings file (defaults to ``~/.delfin/settings.json``); pass an
    explicit path to repair a faked home in a test.
    """
    if not approved:
        raise ValueError(
            "repair step not approved: "
            f"{step.get('kind')} {step.get('check')} refuses to run "
            "without approved=True")
    settings = Path(user_settings_path) if user_settings_path else \
        DEFAULT_SETTINGS_PATH
    if step.get("kind") == "setting":
        return _apply_setting(step, settings, on_line)
    if step.get("kind") == "command":
        return _apply_command(step, on_line)
    raise ValueError(f"unknown repair step kind: {step.get('kind')!r}")
