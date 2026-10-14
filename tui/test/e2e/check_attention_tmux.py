"""What needs the user, live: `delfin-agent tui` in tmux (full screen)
against a private agent service and a private home, into whose attention
inbox and job-monitor findings this check writes what other DELFIN
processes would: three parked questions (one of a saved session, one of
a session another process still holds, one of no session), a failed job
with the monitor's finding, a finished job.

The hints line names what arrived at its right, then counts it (◉ 3 · ✗ 1)
and the terminal's title says (4); ← lists it under NEEDS YOU and
NOTICES. ctrl+x dismisses the notice at once; typing + enter on the
peeked question of no session answers it in the inbox; → on a question
of a saved session resumes it (and closes the question: it is answered
there); → on one whose session is held elsewhere opens a copy; → on the
failed job starts a session about it whose context row ctrl+o shows in
full, with a suggested prompt in the box -- nothing sent.

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_attention_tmux.py
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[2]
REPO = TUI.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(TUI / "spike"))

from check_tmux import Tmux  # noqa: E402

from delfin.agent.service.client import ServiceClient  # noqa: E402

PYTHON = os.environ.get("DELFIN_PYTHON") or sys.executable

#: Run with the service's environment: the paths its modules read, and
#: what this check writes there (or holds: a session's writer lock).
HELPER = r"""
import json, sys, time
from delfin.agent import state_paths
state_paths.scratch_state_from_environment()
from delfin.agent import attention, job_monitor, session_store
what = sys.argv[1]
if what == "paths":
    print(json.dumps({"inbox": str(attention._inbox_path()), "findings": str(job_monitor._FINDINGS_PATH)}))
elif what == "append":
    path, records = sys.argv[2], json.loads(sys.argv[3])
    with open(path, "a", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")
elif what == "hold":
    session_store.acquire_session_lock(sys.argv[2])
    print("held", flush=True)
    time.sleep(600)
"""


def entry(kind: str, title: str, *, ident: str, detail: str = "", session_id: str = "", workspace: str = "",
          options=None, age: float = 0.0) -> dict:
    return {"id": ident, "kind": kind, "session_id": session_id, "title": title, "detail": detail,
            "options": options, "workspace": workspace, "host": os.uname().nodename,
            "created_at": time.time() - age, "status": "pending", "answer": None, "resolved_at": None,
            "acknowledged": False, "delivered": []}


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-attention-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtt", dir="/tmp"))
    ws = tmp / "project"
    calc = tmp / "calc"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws, calc):
        d.mkdir(parents=True)
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO),
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    child_env = {**os.environ, **env}
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(110, 40)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"
    holder: subprocess.Popen | None = None

    def helper(*args: str) -> str:
        return subprocess.run([PYTHON, "-c", HELPER, *args], env=child_env, cwd=str(tmp), check=True,
                              capture_output=True, text=True, timeout=60).stdout

    def expect(pattern: str, what: str, timeout: float = 15) -> re.Match | None:
        found = _screen_wait(tmux, pattern, timeout)
        if not found:
            problems.append(f"{what}: no {pattern!r}\n{tmux.screen()}")
        return found

    def select(text: str) -> bool:
        """Move the agent view's selection to the row that says `text`."""
        for _ in range(16):
            if any(re.match(r"^\s*›", line) and text in line for line in tmux.screen().splitlines()):
                return True
            tmux.key("Down")
            time.sleep(0.25)
        problems.append(f"no row with {text!r} to select\n{tmux.screen()}")
        return False

    def inbox() -> dict[str, dict]:
        out = {}
        for line in Path(paths["inbox"]).read_text().splitlines():
            if line.strip():
                rec = json.loads(line)
                out[rec["id"]] = rec
        return out

    def title() -> str:
        return tmux.run("display-message", "-p", "-t", tmux.name, "#{pane_title}").strip()

    try:
        paths = json.loads(helper("paths"))
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window did not come up", 60):
            raise SystemExit
        # Two sessions that ran and ended (saved in the store), as a
        # dashboard's would be; the second one's lock then held by another
        # process, as a dashboard that still has it open holds it.
        ids = {}
        with ServiceClient(sock) as other:
            for name in ("benzene-scan", "held-elsewhere"):
                sid = other.request("session.create", {"cwd": str(ws), "fake": True, "fakeDelay": 0.003,
                                                       "fakeSave": True, "name": name}, timeout=60)["sessionId"]
                mark = len(other.events)
                other.request("session.prompt", {"sessionId": sid, "text": f"hello from {name}"})
                other.wait_for(lambda e, s=sid: e.get("type") == "turn.ended" and e.get("session") == s,
                               after=mark, timeout=30)
                other.request("session.shutdown", {"sessionId": sid}, timeout=60)
                ids[name] = sid
            deadline = time.monotonic() + 30
            while time.monotonic() < deadline:
                running = {s["sessionId"] for s in other.request("sessions.list")["sessions"]}
                if not running & set(ids.values()):
                    break
                time.sleep(0.3)
        holder = subprocess.Popen([PYTHON, "-c", HELPER, "hold", ids["held-elsewhere"]], env=child_env,
                                  cwd=str(tmp), stdout=subprocess.PIPE, text=True)
        if holder.stdout.readline().strip() != "held":
            problems.append("the helper could not hold the session's lock")
        # What the dashboard, the REPL and the job monitor would have written.
        helper("append", paths["findings"], json.dumps([{
            "job_id": "4711", "folder": str(calc), "state": "OUT_OF_MEMORY", "signatures": ["oom-kill"],
            "ts": time.time(), "diagnosis_session": "", "summary": "the SCF ran out of memory",
            "fix": {"fix_class": "resource", "recommendation": "raise maxcore",
                    "fix": {"kind": "memory", "summary_line": "maxcore 3000 -> 4500",
                            "diff": "-%maxcore 3000\n+%maxcore 4500"}}}]))
        helper("append", paths["inbox"], json.dumps([
            entry("question_pending", "Question waiting: Which functional?", ident="att-qx", age=3,
                  detail="Which functional for the excited states?", options=["B3LYP", "CAM-B3LYP"],
                  session_id=ids["benzene-scan"], workspace=str(ws)),
            entry("question_pending", "Question waiting: Which basis set?", ident="att-qfree", age=2,
                  detail="Which basis set for the scan?", workspace=str(ws)),
            entry("question_pending", "Question waiting: Which solvent?", ident="att-qy", age=1,
                  detail="Which solvent model?", session_id=ids["held-elsewhere"], workspace=str(ws)),
            entry("run_failed", "Job 4711 failed (OUT_OF_MEMORY)", ident="att-failed", age=5,
                  detail="out of memory", workspace=str(calc)),
            entry("run_finished", "Job 4702 finished (COMPLETED)", ident="att-done", age=6,
                  detail="finished after 2h14m", workspace=str(calc)),
            # Written on another machine: its folder is not here.
            entry("run_finished", "Job 4800 finished (COMPLETED)", ident="att-far", age=7,
                  detail="finished on the other cluster", workspace="/nonexistent-delfin-check/calc"),
        ]))
        # Quiet, at the bottom right: what arrived names itself, then the count.
        expect(r"[◉✗] (job 4711 failed|\S+ asks)\s*$|◉ 3 · ✗ 1\s*$", "the note at the hints' right", 20)
        expect(r"agents\s+◉ 3 · ✗ 1\s*$", "the count at the hints' right", 20)
        if not re.match(r"\(4\) DELFIN", title()):
            problems.append(f"the terminal's title should count 4: {title()!r}")
        # The agent view: NEEDS YOU and NOTICES.
        tmux.key("Left")
        expect(r"◆ AGENTS\s+\d+ sessions? · 4 need you", "the agent view")
        screen = tmux.screen()
        for text in ("NEEDS YOU", "NOTICES", "job 4711 failed", "a question in project", "Which basis set",
                     "held-elsewhere", "open elsewhere", "job 4702 finished"):
            if text not in screen:
                problems.append(f"the agent view lacks {text!r}\n{screen}")
        if "NOTICES" in screen and "DONE" in screen and screen.index("DONE") > screen.index("NOTICES"):
            problems.append("NOTICES should come after DONE")
        # ctrl+x on the notice: gone at once.
        if select("job 4702 finished"):
            tmux.key("C-x")
            expect(r"dismissed: job 4702 finished", "the dismiss")
            time.sleep(0.5)
            done = inbox()["att-done"]
            if done.get("status") != "resolved" or not done.get("acknowledged"):
                problems.append(f"the notice should be dismissed in the inbox: {done}")
        # The question of no session: peeked and answered from the view.
        if select("Which basis set"):
            tmux.key("Space")
            expect(r"│ Which basis set for the scan\?", "the question's peek")
            expect(r"type the answer and press enter", "the peek's hint")
            tmux.type("def2-TZVP")
            tmux.key("Enter")
            expect(r"answered: the next session started in", "the answer's note")
            time.sleep(0.5)
            got = inbox()["att-qfree"]
            if got.get("status") != "resolved" or got.get("answer") != "def2-TZVP" or got.get("acknowledged"):
                problems.append(f"the answer should wait in the inbox for the next session: {got}")
        # → on the question of a saved session: resumed here, the question closed.
        if select("Which functional"):
            tmux.key("Right")
            # Full screen: the session's own rows from the top (no divider).
            expect(r"DELFIN .*· benzene-scan", "the resumed session's banner", 30)
            expect(r"hello from benzene-scan", "its replay", 30)
            expect(r"it asked, and nobody was there to answer -- answer it here", "what it asked", 15)
            expect(r"Which functional for the excited states\?", "the question under the replay", 15)
            time.sleep(0.5)
            got = inbox()["att-qx"]
            if got.get("status") != "resolved" or not got.get("acknowledged"):
                problems.append(f"the question should be closed once its session is open: {got}")
        # → on the question of a session another process holds: a copy.
        time.sleep(0.5)
        tmux.key("Left")
        expect(r"◆ AGENTS", "the agent view again")
        if select("Which solvent"):
            tmux.key("Right")
            expect(r"DELFIN .*· held-elsewhere \(fork\)", "the copy's banner", 30)
            expect(r"earlier in this session, as it was saved", "the copy's history notice", 30)
            expect(r"hello from held-elsewhere", "the copy's history", 30)
            expect(r"this is a copy of held-elsewhere, which is still open elsewhere", "the copy's notice", 30)
        # → on the failed job: a session about it, its context, a suggested prompt.
        time.sleep(0.5)
        tmux.key("Left")
        expect(r"◆ AGENTS", "the agent view once more")
        if select("job 4711 failed"):
            tmux.key("Right")
            expect(r"⏺ context  job 4711 failed", "the context row", 30)
            expect(r"─ job 4711 failed ─", "the session's name on the prompt box", 15)
            expect(r"⎿ \d+ lines, given to the model with the first prompt", "its summary", 15)
            expect(r"❯ Job 4711 failed with OUT_OF_MEMORY\. What went wrong, and what should we change\?",
                   "the suggested prompt", 15)
            tmux.key("C-o")
            expect(r"A fix is prepared \(not applied\): maxcore 3000 -> 4500", "the context's details", 15)
            expect(rf"Folder: {re.escape(str(calc))}", "the job's folder in the context", 15)
            tmux.key("C-o")
        # → on one whose folder is not on this machine: a session about it
        # here, and the note saying so -- after the session is in view (full
        # screen clears the screen when it shows one).
        time.sleep(0.5)
        tmux.key("C-c")
        time.sleep(0.4)
        tmux.key("Left")
        expect(r"◆ AGENTS", "the agent view a last time")
        if select("job 4800 finished"):
            tmux.key("Right")
            expect(r"⏺ context  job 4800 finished", "the far job's context row", 30)
            expect(r"/nonexistent-delfin-check/calc could not take a session", "the folder note", 15)
            expect(rf"started in {re.escape(str(ws))}", "where it started instead", 15)
        time.sleep(0.5)
        tmux.key("C-c")                       # the suggested prompt goes, unsent
        time.sleep(0.4)
        tmux.key("C-d")
        exited = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not exited or exited.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {exited and exited.group(1)}")
        with ServiceClient(sock) as cl:
            sessions = cl.request("sessions.list")["sessions"]
            about = [s for s in sessions if s.get("workspace") == str(calc)]
            if not any(s.get("title") == "job 4800 finished" and s.get("workspace") == str(ws) for s in sessions):
                problems.append(f"the far job's session should run in {ws}: {sessions}")
            if len(about) != 1:
                problems.append(f"one session should run in {calc}: {sessions}")
            else:
                if about[0].get("title") != "job 4711 failed":
                    problems.append(f"the session about the job should be named after it: {about[0]}")
                events = cl.request("session.replay", {"sessionId": about[0]["sessionId"]})["events"]
                if not any(e.get("type") == "context" for e in events):
                    problems.append("the session about the job has no context event")
                said = [e for e in events if e.get("type") == "user.message"]
                if said:
                    problems.append(f"nothing should have been sent to the session about the job: {said}")
            copies = [s for s in sessions if s["sessionId"] not in ids.values()
                      and "held-elsewhere" in str(s.get("title") or "")]
            if not copies:
                problems.append(f"a copy of held-elsewhere should run under a new id: {sessions}")
            if any(s["sessionId"] == ids["held-elsewhere"] for s in sessions):
                problems.append("the held session itself should not have been opened")
            cl.request("service.shutdown")
    except SystemExit:
        pass
    finally:
        if holder is not None:
            holder.terminate()
            try:
                holder.wait(timeout=10)
            except subprocess.TimeoutExpired:
                holder.kill()
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux attention check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: what needs the user is counted, listed, dismissed, answered and opened from the window, live")
    return 0


def _screen_wait(tmux: Tmux, pattern: str, timeout: float):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        m = re.search(pattern, tmux.screen(), re.MULTILINE)
        if m:
            return m
        time.sleep(0.1)
    return None


if __name__ == "__main__":
    sys.exit(main())
