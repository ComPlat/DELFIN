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
"""
from __future__ import annotations

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


def _summary(log: Path) -> dict:
    tools, result = 0, {}
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith('{"type": "tool_use"'):
            tools += 1
        elif line.startswith('{"type": "result"'):
            try:
                result = json.loads(line)
            except ValueError:
                pass
    return {"stop": result.get("stop_reason"), "tools": tools,
            "in": result.get("input_tokens"), "out": result.get("output_tokens"),
            "answer": (result.get("text") or "")[-1200:], "error": (result.get("error") or "")[:300]}


def _one(probe: dict, side: str, code: Path, provider: str, model: str, key: str,
         out: Path) -> dict:
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


def run(run_name: str, tid: str) -> dict:
    """The member's real runs, base and branch; returns the rows and writes
    ``scratch/realrun/summary.md``."""
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
    rows = []
    with tempfile.TemporaryDirectory(prefix="delfin-realrun-code-") as codes:
        sides = {"base": str(task.get("base") or "origin/main"), "branch": "HEAD"}
        for side, ref in sides.items():
            _archive(wt if side == "branch" else repo, ref, Path(codes) / side)
        for probe in req["probes"]:
            for side in sides:
                row = _one(probe, side, Path(codes) / side, provider, model, key, out)
                rows.append({"probe": probe["name"], "side": side, **row})
    key = ""
    lines = [f"# Real runs for {run_name}/{tid}", "",
             f"{provider} / {model}; each probe on the base ({sides['base']}) and on the "
             "branch HEAD, fresh HOME and workspace each.", "",
             "| probe | side | stop | tool calls | in tok | out tok | s | error |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['probe']} | {r['side']} | {r.get('stop')} | {r.get('tools')} | "
                     f"{r.get('in')} | {r.get('out')} | {r.get('s')} | "
                     f"{(r.get('error') or '').replace('|', '/')[:120]} |")
    for r in rows:
        lines += ["", f"## {r['probe']} / {r['side']}: final answer (last 1200 characters)", "",
                  "```", r.get("answer") or "(none)", "```"]
    (out / "summary.md").write_text(_redact("\n".join(lines) + "\n"))
    wave.record(run_name, "realrun", tid, rows=[{k: r.get(k) for k in ("probe", "side", "stop", "tools", "error")}
                                                for r in rows])
    return {"run": run_name, "id": tid, "rows": rows, "summary": str(out / "summary.md")}
