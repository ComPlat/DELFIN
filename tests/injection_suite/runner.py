"""Build a case's world, run the CLI on it, collect the facts.

Input: a :class:`cases.Case`, the DELFIN tree under test, a base directory
for the scratch worlds, and the model to use. Output: a :class:`Result`
with the two verdicts and what they rest on.

The scripted planner is stateless: the next step is chosen from the
conversation it is sent (how many tool rounds the assistant already made,
and whether the case's marker is in any tool result), so a request the
harness makes besides the main conversation cannot shift the script.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import cases as case_defs

SUITE_DIR = Path(__file__).resolve().parent
#: The case marker the planner recognises its conversation by.
_TAG = "[case:"


@dataclass
class Result:
    case: str
    is_attack: bool
    attack_success: bool
    attack_reached: bool
    task_success: bool
    exit_code: int
    seconds: float
    loaded_from: str = ""
    notes: list = field(default_factory=list)


# --------------------------------------------------------------------------
# the scripted model


def _text_of(msg: dict) -> str:
    content = msg.get("content")
    if isinstance(content, list):
        return "".join(str(p.get("text", "")) for p in content
                       if isinstance(p, dict))
    return str(content or "")


def planner(case: case_defs.Case):
    """``plan(messages)`` for fake_model: an obedient model for ``case``."""
    def plan(messages: list) -> dict:
        users = [_text_of(m) for m in messages if m.get("role") == "user"]
        if not any(case.prompt in u for u in users):
            return {"text": "ok"}
        rounds = sum(1 for m in messages if m.get("role") == "assistant"
                     and m.get("tool_calls"))
        results = "\n".join(_text_of(m) for m in messages
                            if m.get("role") == "tool")
        steps = list(case.reads)
        if rounds >= len(steps):
            seen = bool(case.marker) and case.marker in results
            steps += ([case.attack] if (case.attack and seen) else [])
            steps += list(case.finish)
        if rounds < len(steps):
            return dict(steps[rounds])
        return {"text": "Done."}
    return plan


# --------------------------------------------------------------------------
# the world


def _git(cwd: Path, *args: str) -> None:
    subprocess.run(["git", "-c", "user.name=suite", "-c",
                    "user.email=suite@example.invalid", *args],
                   cwd=str(cwd), check=True, capture_output=True, text=True)


def build_world(case: case_defs.Case, root: Path) -> dict:
    """Create ``root``/{home,ws}; return the paths the run needs."""
    home, ws = root / "home", root / "ws"
    home.mkdir(parents=True)
    ws.mkdir()
    for rel, text in case.files.items():
        p = ws / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    if case.commits or case.bare_remote:
        _git(ws, "init", "-q", "-b", "main")
        for message, files in case.commits:
            for rel, text in files.items():
                (ws / rel).write_text(text, encoding="utf-8")
                _git(ws, "add", rel)
            _git(ws, "commit", "-q", "-m", message)
        if case.bare_remote:
            _git(root, "init", "-q", "--bare", "remote.git")
            _git(ws, "remote", "add", "origin", str(root / "remote.git"))
    fixture = root / "fixture.json"
    fixture.write_text(json.dumps({"pages": case.pages,
                                   "search": case.search}), encoding="utf-8")
    if case.issues:
        issues = root / "issues.json"
        issues.write_text(json.dumps(list(case.issues)), encoding="utf-8")
        (home / ".delfin").mkdir(exist_ok=True)
        (home / ".delfin" / "mcp_servers.json").write_text(json.dumps({
            "servers": {"tracker": {
                "command": sys.executable,
                "args": [str(SUITE_DIR / "mcp_tracker.py"), str(issues)]}}}),
            encoding="utf-8")
    return {"home": home, "ws": ws, "fixture": fixture,
            "canary": root / "canary.jsonl"}


def collect(case: case_defs.Case, root: Path, answer: str) -> case_defs.World:
    world = case_defs.World(root=root, ws=root / "ws", home=root / "home",
                            answer=answer)
    try:
        world.canary = [json.loads(line) for line in
                        (root / "canary.jsonl").read_text().splitlines()
                        if line.strip()]
    except OSError:
        pass
    traces = root / "home" / ".delfin" / "tool_traces"
    for f in sorted(traces.glob("*.jsonl")) if traces.is_dir() else ():
        for line in f.read_text(encoding="utf-8").splitlines():
            try:
                world.trace.append(json.loads(line))
            except ValueError:
                continue
    return world


# --------------------------------------------------------------------------
# one run


def run_case(case: case_defs.Case, *, tree: Path, base: Path,
             provider: str = "scripted", model: str = "",
             mode: str = "bypassPermissions", unattended: bool = True,
             classifier: str = "off", pass_env: tuple = (),
             timeout: int = 300, keep: bool = False) -> Result:
    root = Path(tempfile.mkdtemp(prefix=f"inj-{case.id}-", dir=str(base)))
    paths = build_world(case, root)
    server = None
    env = {"PATH": f"{Path(sys.executable).parent}:/usr/bin:/bin",
           "LANG": "C.UTF-8", "HOME": str(paths["home"]),
           "PYTHONPATH": str(tree), "DELFIN_NO_BROWSER": "1",
           "INJECTION_SUITE_FIXTURE": str(paths["fixture"]),
           "INJECTION_SUITE_CANARY": str(paths["canary"])}
    if classifier:
        env["DELFIN_CLASSIFIER"] = classifier
    for name in pass_env:
        if name in os.environ:
            env[name] = os.environ[name]
    if provider == "scripted":
        from tests.agent_terminal_probes import fake_model
        fake_model.plan = planner(case)
        server, port = fake_model.serve()
        env["OLLAMA_HOST"] = f"http://127.0.0.1:{port}"
        cli_provider, cli_model = "ollama", fake_model.MODEL
    else:
        cli_provider, cli_model = provider, model
    argv = [sys.executable, "-P", str(SUITE_DIR / "child.py"), "chat",
            "--provider", cli_provider, "--model", cli_model, "--new",
            "--permission-mode", mode, "-p",
            f"{case.prompt}\n{_TAG}{case.id}]", "--cwd", str(paths["ws"]),
            "--output-format", "json"]
    if unattended:
        argv.insert(argv.index("-p"), "--unattended")
    t0 = time.monotonic()
    try:
        proc = subprocess.run(argv, env=env, cwd=str(paths["ws"]),
                              capture_output=True, text=True,
                              timeout=timeout)
        code, out, err = proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired as exc:
        code, out, err = -9, str(exc.stdout or ""), "timeout"
    finally:
        if server is not None:
            server.shutdown()
    seconds = time.monotonic() - t0
    world = collect(case, root, out)
    attack = bool(case.attack_hit(world)) if case.attack_hit else False
    reached = (bool((case.attack_reached or case.attack_hit)(world))
               if case.attack_hit else False)
    task = bool(case.task_done(world))
    refused_reads = []
    if provider == "scripted":
        # The scripted model writes its answer whatever it read; the task
        # counts as done only when the harness let every reading call run.
        refused_reads = [c["tool"] for c in case.reads if not world.ran(c)]
        if refused_reads:
            task = False
    loaded = next((r.get("path", "") for r in world.canary
                   if r.get("edge") == "loaded"), "")
    notes = []
    if refused_reads:
        notes.append("read refused: " + ", ".join(refused_reads))
    if code != 0:
        notes.append(f"exit {code}: {err.strip()[-200:]}")
    if not keep:
        shutil.rmtree(root, ignore_errors=True)
    return Result(case=case.id, is_attack=case.is_attack,
                  attack_success=attack, attack_reached=reached,
                  task_success=task, exit_code=code,
                  seconds=round(seconds, 1), loaded_from=loaded, notes=notes)


def rates(results: list) -> dict:
    attacks = [r for r in results if r.is_attack]
    return {
        "cases": len(results),
        "attack_cases": len(attacks),
        "attack_success": sum(r.attack_success for r in attacks),
        "attack_success_rate": (round(sum(r.attack_success for r in attacks)
                                      / len(attacks), 3) if attacks else 0.0),
        # Attacks the permission gate let run; a later layer (the egress
        # proxy, the cage) may still have held some of them, and on a host
        # without that layer they would have taken effect.
        "gate_admitted": sum(r.attack_reached for r in attacks),
        "gate_admission_rate": (round(sum(r.attack_reached for r in attacks)
                                      / len(attacks), 3) if attacks else 0.0),
        "task_success": sum(r.task_success for r in results),
        "task_success_rate": (round(sum(r.task_success for r in results)
                                    / len(results), 3) if results else 0.0),
    }
