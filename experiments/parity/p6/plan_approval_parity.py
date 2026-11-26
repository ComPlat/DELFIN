#!/usr/bin/env python3
"""Parity harness: plan -> approval -> execution, and an honest report.

Measures one capability from ONE package -- a headless (``-p``) session that
plans, is approved by a person, executes, and reports what it did -- through
the PRODUCT ENTRY POINT, never through a helper:

    python -P -m delfin.agent.cli chat --provider ollama --model fake-model
        --new --permission-mode plan -p <marker> --output-format stream-json

The model is the scripted OpenAI-compatible endpoint of
``tests/agent_terminal_probes/fake_model.py``; its ``plan`` is replaced from
HERE (that module's own docstring documents the extension), so the shared
probe file is not edited and no live model is needed. Each variant is one
child process with a private HOME (conftest.child_env), so nothing reaches
the operator's real ~/.delfin.

A person's approval is given the way a person gives it: a second child runs

    python -P -m delfin.agent.cli approvals approve|deny <id>

with the same HOME and NO ``DELFIN_AGENT_SESSION`` marker, which is what a
terminal run by a human has. It reads the one-time code from the code file the
asking run wrote for a job with no terminal (``approval_code.write_code_file``)
-- the harness never sees the code, so no check is weakened and the code never
lands in a log. The asking child is started with ``start_new_session`` for the
same reason a batch job has none: with no controlling terminal the code goes
to its code file instead of to /dev/tty.

Per variant the script records three things:

  REACHED  the plan/approval code ran          (the tool result that came
                                               back to the model proves it,
                                               from $FAKE_MODEL_LOG)
  USED     a question actually reached the     (a *.request.json appeared in
           approvals directory where one is     that HOME's approvals room)
           needed
  CORRECT  the outcome AND the report are      (the result object's
           what they should be                  stop_reason + exit code, and
                                               files_changed: nothing may be
                                               written before an approval)

A security refusal that is correct counts as CORRECT: the run an agent tries
to answer for itself must never be answered, on either tree.

Usage
    python experiments/parity/p6/plan_approval_parity.py --repo <tree> \
        --out <file.json> [--only id,id]

Exit 0 when every variant ran (a CORRECT=no row is a finding, not a crash).
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path


# --- the scripted model ----------------------------------------------------

PLAN_TEXT = ("1. Write notes.txt with the plan's result\n"
             "2. Report what was written")


def _steps_for(kind: str) -> list[dict]:
    """The scripted model's turns, one dict per request, in order.

    A variant whose real request count exceeds its steps repeats the last
    step; every list ends in a text step, so a run always terminates.
    """
    submit = {"tool": "exit_plan_mode", "args": {"plan": PLAN_TEXT}}
    write = {"tool": "write_file",
             "args": {"path": "notes.txt", "content": "planned\n"}}
    read = {"tool": "read_file", "args": {"path": "readme.txt"}}
    done = {"text": "That is the whole answer."}
    refused = {"text": "Understood, nothing more from me."}
    return {
        # plan mode refuses a write: the read-only contract holds, and the
        # turn must be reported as refused, not as finished.
        "plan_write_refused": [write, refused],
        # no --supervised: there is genuinely nobody to ask.
        "submit_no_channel": [submit, done],
        # --supervised, a person watches the approvals directory.
        "submit_supervised_approve": [submit, write, done],
        "submit_supervised_deny": [submit, refused],
        # asked, and nobody answered in time (not a rejection).
        "submit_supervised_timeout": [submit, refused],
        # a malformed plan: nothing is submitted, nothing is pending.
        "empty_plan": [{"tool": "exit_plan_mode", "args": {"plan": "   "}}, done],
        # a second submission while the first is still on screen.
        "duplicate_submit": [submit, submit, done],
        # a write attempted BEFORE the plan is submitted: 0 changes may
        # exist before an approval, whatever else the run reports.
        "plan_write_then_submit": [write, submit, refused],
        # --- controls, expected identical on both trees ---
        # exit_plan_mode outside plan mode is an error, not a gate.
        "exit_plan_in_default_mode": [submit, done],
        # default mode: a plain write is allowed and the turn is done.
        "default_write_control": [write, done],
        # an ordinary planning turn: read, answer, submit nothing.
        "plan_readonly_no_submit": [read, done],
        # the run is person-launched (so it CAN be answered); only the
        # ANSWERER carries the agent marker. Same turns as no-channel.
        "answerer_is_an_agent_refused": [submit, refused],
    }[kind]


# --- what each variant must satisfy ---------------------------------------

MODES = {
    # id: (permission mode, --supervised, confirm-timeout, answer action,
    #      keep the agent marker in the ANSWERER, expected)
    "plan_write_refused": ("plan", False, 0, "none", False,
                           {"stop": "refused", "code": 5, "before": []}),
    "submit_no_channel": ("plan", False, 0, "none", False,
                          {"stop": "awaiting_approval", "code": 11,
                           "before": []}),
    "submit_supervised_approve": ("plan", True, 25, "approve", False,
                                  {"stop": "done", "code": 0,
                                   "before": [], "after_ok": True}),
    "submit_supervised_deny": ("plan", True, 25, "deny", False,
                               {"stop": "refused", "code": 5, "before": []}),
    "submit_supervised_timeout": ("plan", True, 3, "none", False,
                                  {"stop": "awaiting_approval", "code": 11,
                                   "before": []}),
    "empty_plan": ("plan", False, 0, "none", False,
                   {"stop": "done", "code": 0, "before": []}),
    "duplicate_submit": ("plan", False, 0, "none", False,
                         {"stop": "awaiting_approval", "code": 11,
                          "before": []}),
    "plan_write_then_submit": ("plan", False, 0, "none", False,
                               {"stop": "awaiting_approval", "code": 11,
                                "before": []}),
    "exit_plan_in_default_mode": ("default", False, 0, "none", False,
                                  {"stop": "done", "code": 0}),
    "default_write_control": ("default", False, 0, "none", False,
                              {"stop": "done", "code": 0}),
    "plan_readonly_no_submit": ("plan", False, 0, "none", False,
                                {"stop": "done", "code": 0, "before": []}),
    # SECURITY CONTROL: an agent session may not answer a question for the
    # run that waits, on either tree. The asking run IS supervised, so a real
    # question exists to be answered -- otherwise this row would assert
    # nothing -- and nobody else answers it, so it expires after 3 s.
    "answerer_is_an_agent_refused": ("plan", True, 3, "agent-deny", True,
                                     {"answer_refused": True}),
}

# Variants whose point is that a person's approval channel was reached.
NEEDS_CHANNEL = ("submit_supervised_approve", "submit_supervised_deny",
                 "submit_supervised_timeout")

ALL_IDS = list(MODES)


# --- the run ---------------------------------------------------------------

def _child_env(home: Path, repo: Path, *, marker: bool) -> dict:
    """A child's environment: private HOME, no credential, no agent marker.

    Mirrors conftest.child_env. ``DELFIN_AGENT_SESSION`` is dropped unless
    ``marker`` is set, because the marker marks a process an agent started:
    its approval code is then shown nowhere and it cannot be answered from
    outside at all. The person-launched case is the one without it.
    """
    import re
    home.mkdir(parents=True, exist_ok=True)
    (home.parent / "scratch" / ".delfin").mkdir(parents=True, exist_ok=True)
    env = {}
    for k, v in os.environ.items():
        if re.search(r"(KEY|TOKEN|SECRET|PASSW|CRED|AUTH)", k, re.I):
            continue
        env[k] = v
    env["HOME"] = str(home)
    env["DELFIN_SCRATCH_STATE"] = str(home.parent / "scratch")
    env["PYTHON_KEYRING_BACKEND"] = "keyring.backends.null.Keyring"
    env["PYTHONPATH"] = str(repo)
    env["DELFIN_CLASSIFIER"] = "off"
    # The three the suite's own parent gets from conftest (996-997 set the
    # retry knobs, the autouse fixture at 1329 turns the process guard off).
    # Without the guard off, a child of a NON-pytest parent cannot connect
    # to the loopback endpoint the harness listens on: every request dies
    # as APIConnectionError and the openai backoff makes the run look hung
    # (measured here: FAIL -> PASS flips on this one variable alone). With
    # the retries scaled down, a real connection failure ends the variant
    # in seconds and is reported, not waited out.
    env["DELFIN_PROCESS_GUARD"] = "off"
    env["DELFIN_STREAM_RETRY_MAX"] = "3"
    env["DELFIN_STREAM_RETRY_DELAY_SCALE"] = "0"
    env.pop("DELFIN_LIFELINE_PID", None)
    env.pop("NO_COLOR", None)
    if not marker:
        env.pop("DELFIN_AGENT_SESSION", None)
    return env


def _answer(home: Path, repo: Path, request_id: str, verb: str,
            marker: bool) -> dict:
    """Run the person's approvals command; return what it did."""
    env = _child_env(home, repo, marker=marker)
    env["OLLAMA_HOST"] = "http://127.0.0.1:1"     # no model is needed here
    cmd = [sys.executable, "-P", "-m", "delfin.agent.cli", "approvals",
           verb, request_id]
    try:
        p = subprocess.run(cmd, cwd=str(home), env=env, capture_output=True,
                           text=True, timeout=120, stdin=subprocess.DEVNULL,
                           start_new_session=True)
        return {"rc": p.returncode,
                "out": (p.stdout + p.stderr)[-400:]}
    except subprocess.TimeoutExpired:
        return {"rc": None, "out": "answerer timed out"}


def _watch_and_answer(home: Path, repo: Path, ws: Path, mode: dict,
                      log: list, stop: threading.Event) -> None:
    """Answer every question this run publishes, as the person would.

    An approved plan flips the session out of plan mode; the next write then
    raises a CONFIRM question of its own, which a watching person answers
    too. The workspace is snapshotted at the FIRST request: that is the
    measured state "before any approval", which must hold no plan artifact.
    """
    verb = {"approve": "approve", "deny": "deny",
            "agent-deny": "deny"}[mode[3]]
    marker = mode[4]
    room = home / ".delfin" / "approvals"
    answered: set = set()
    seen_any = False
    deadline = time.time() + 120
    while not stop.is_set() and time.time() < deadline:
        for path in sorted(glob.glob(str(room / "*.request.json"))):
            try:
                rec = json.loads(Path(path).read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            rid = str(rec.get("id") or "")
            if not rid or rid in answered:
                continue
            if not seen_any:
                seen_any = True
                log.append({"seen": True, "tool": rec.get("tool", ""),
                            "id": rid,
                            "ws_at_first_request": sorted(
                                os.listdir(ws)) if ws.is_dir() else []})
            out = _answer(home, repo, rid, verb, marker=marker)
            answered.add(rid)
            log.append({"answered": verb, "rc": out["rc"],
                        "out": out["out"]})
        time.sleep(0.15)
    if not seen_any:
        log.append({"seen": False})


def run_variant(kind: str, repo: Path, plan_fn_holder: dict) -> dict:
    """One variant: one child of the real CLI, scored on REACHED/USED/CORRECT."""
    mode = MODES[kind]
    perm, supervised, confirm_timeout, answer_action, marker, expect = mode
    root = Path(os.environ.get("DELFIN_PARITY_TMP") or "/tmp") / f"p6-{kind}"
    if root.exists():
        for p in sorted(root.rglob("*"), reverse=True):
            try:
                p.rmdir() if p.is_dir() else p.unlink()
            except OSError:
                pass
    home = root / "home"
    ws = root / "ws"
    ws.mkdir(parents=True, exist_ok=True)
    (ws / "readme.txt").write_text("readable\n", encoding="utf-8")

    env = _child_env(home, repo, marker=False)
    env["OLLAMA_HOST"] = plan_fn_holder["url"]
    # The request log belongs to the SERVER process (this one), not to the
    # child: fake_model writes it from the handler. Point the shared server
    # at this variant's file rather than putting FAKE_MODEL_LOG in the
    # child's environment, where nothing would ever read it.
    log_path = root / "model.log"
    plan_fn_holder["server"].log_path = str(log_path)

    cmd = [sys.executable, "-P", "-m", "delfin.agent.cli", "chat",
           "--provider", "ollama", "--model", "fake-model", "--new",
           "--permission-mode", perm, "--cwd", str(ws),
           "--output-format", "stream-json", "-p", f"parity {kind}"]
    if supervised:
        cmd.append("--supervised")
    if confirm_timeout:
        cmd += ["--confirm-timeout", str(confirm_timeout)]

    # The scripted model for THIS variant, installed on the server thread.
    steps = _steps_for(kind)
    state = {"i": 0}

    def _plan(messages: list) -> dict:
        step = steps[min(state["i"], len(steps) - 1)]
        state["i"] += 1
        return step

    plan_fn_holder["set_plan"](_plan)

    watch_log: list = []
    stop_watch = threading.Event()
    if answer_action in ("approve", "deny", "agent-deny"):
        t = threading.Thread(target=_watch_and_answer,
                             args=(home, repo, ws, mode, watch_log, stop_watch),
                             daemon=True)
        t.start()

    started = time.time()
    try:
        p = subprocess.run(cmd, cwd=str(ws), env=env, capture_output=True,
                           text=True, timeout=300, stdin=subprocess.DEVNULL,
                           start_new_session=True)
        rc, out = p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired as exc:
        rc = None
        out = (exc.stdout or b"").decode() if isinstance(exc.stdout, bytes) \
            else (exc.stdout or "")
    if answer_action in ("approve", "deny", "agent-deny"):
        stop_watch.set()
        t.join(timeout=15)
    else:
        watch_log.append({"seen": False, "not_needed": True})

    events = []
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("{") and '"type"' in line:
            try:
                events.append(json.loads(line))
            except ValueError:
                pass
    result = next((e for e in reversed(events) if e.get("type") == "result"),
                  {})

    tool_results = []
    if log_path.is_file():
        for line in log_path.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line)
            except ValueError:
                continue
            if rec.get("tool_result"):
                tool_results.append(str(rec["tool_result"]))

    if answer_action == "agent-deny":
        # The security control: the answer attempt itself must be refused.
        reached = bool(watch_log)
        used = any(w.get("seen") for w in watch_log)
        answer_rc = next((w.get("rc") for w in watch_log if "answered" in w),
                         None)
        correct = answer_rc not in (0, None)
        detail = {"answer_rc": answer_rc,
                  "answer_out": next((w.get("out", "") for w in watch_log
                                      if "answered" in w), "")}
    else:
        reached = bool(tool_results)
        used = any(w.get("seen") for w in watch_log)
        stop = str(result.get("stop_reason") or "")
        files = [str(f) for f in (result.get("files_changed") or [])]
        approved = answer_action == "approve"
        first_seen = next((w for w in watch_log if w.get("seen")), {})
        ws_before = first_seen.get("ws_at_first_request", [])
        # "changes made before approval = 0": the workspace held no plan
        # artifact at the moment the question was on screen.
        before_ok = "notes.txt" not in ws_before
        if approved and expect.get("after_ok"):
            # and it exists only AFTER the approval: the plan executed.
            before_ok = before_ok and (ws / "notes.txt").is_file()
        want = expect.get("stop")
        correct = (stop == want and rc == expect.get("code")
                   and before_ok)
        if kind in NEEDS_CHANNEL:
            # A row whose point is that a person was asked cannot pass by
            # never having asked: on the untouched tree the supervised run
            # refuses the follow-up write on its own and the turn happens to
            # end "refused", which would look like a correct denial.
            correct = correct and used
        detail = {"stop": stop, "rc": rc, "files": files,
                  "ws_before": ws_before,
                  "n_requests": state["i"],
                  "first_result": tool_results[0][:160] if tool_results else ""}
    return {"id": kind, "reached": bool(reached), "used": bool(used),
            "correct": bool(correct), "seconds": round(time.time() - started, 1),
            "needs_channel": kind in NEEDS_CHANNEL or
            kind == "answerer_is_an_agent_refused",
            "detail": detail}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="tree under test")
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    sys.path.insert(0, str(repo))
    from tests.agent_terminal_probes import fake_model            # noqa: E402

    server, port = fake_model.serve()
    holder = {"url": f"http://127.0.0.1:{port}", "server": server,
              "set_plan": lambda fn: setattr(fake_model, "plan", fn)}
    ids = [i.strip() for i in args.only.split(",") if i.strip()] or ALL_IDS
    rows = []
    try:
        for kind in ids:
            row = run_variant(kind, repo, holder)
            rows.append(row)
            print(f"{kind:30s} reached={'Y' if row['reached'] else 'n'}"
                  f"  used={'Y' if row['used'] else 'n'}"
                  f"  correct={'Y' if row['correct'] else 'n'}"
                  f"  {row['seconds']}s  {row['detail']}", flush=True)
    finally:
        server.shutdown()
        server.server_close()

    n = len(rows)
    wanted = [r for r in rows if r["needs_channel"]]
    summary = {
        "tree": str(repo), "commit": _git(repo), "n": n,
        "reached": sum(1 for r in rows if r["reached"]),
        "used": sum(1 for r in wanted if r["used"]),
        "used_of": len(wanted),
        "correct": sum(1 for r in rows if r["correct"]),
        "rows": rows,
    }
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n",
                              encoding="utf-8")
    print(f"\nREACHED {summary['reached']}/{n}  USED {summary['used']}/"
          f"{summary['used_of']}  CORRECT {summary['correct']}/{n}"
          f" -> {args.out}", flush=True)
    return 0


def _git(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


if __name__ == "__main__":
    sys.exit(main())
