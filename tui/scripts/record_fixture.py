"""Record a fake conversation from a real agent service, as a test fixture.

    DELFIN_PYTHON=... python3 tui/scripts/record_fixture.py

A private service (own HOME, scratch state and runtime directory) runs one
scripted session (the fake engine, no model) through the prompts below; every
event the window would receive is kept, with what changes from run to run
(times, ids, the temporary workspace) made fixed. The TUI's snapshot test
renders the fixture headless; the tmux test runs the same prompts live and
must show the same text (tui/test/fixtures/).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[1]
REPO = TUI.parent
sys.path.insert(0, str(REPO))

from delfin.agent.service.client import ServiceClient  # noqa: E402

#: What the conversation does; the tmux test types the same (fixtures/).
SCRIPT = [
    {"prompt": "hello there"},
    {"prompt": "please edit the parser", "reply": {"choice": "y"}},
    {"prompt": "ask me which solvent", "reply": {"answers": ["SMD"]}},
]
WORKSPACE = "/home/user/project"


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-fixture-"))
    runtime = Path(tempfile.mkdtemp(prefix="dfx", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws):
        d.mkdir(parents=True)
    env = dict(os.environ, HOME=str(tmp / "home"), DELFIN_SCRATCH_STATE=str(tmp / "scratch"),
               DELFIN_RUNTIME_DIR=str(runtime), PYTHONPATH=str(REPO))
    env.pop("DELFIN_LIFELINE_PID", None)
    env.pop("DELFIN_LIFELINE_TICKS", None)
    service = subprocess.Popen([sys.executable, "-m", "delfin.agent.service.supervisor",
                                "--idle-exit", "60"], env=env, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, start_new_session=True)
    sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"
    try:
        for _ in range(200):
            if sock.exists():
                break
            time.sleep(0.05)
        with ServiceClient(sock) as c:
            created = c.request("session.create", {"cwd": str(ws), "fake": True,
                                                   "fakeDelay": 0.003,
                                                   "permissionMode": "default"})
            sid = created["sessionId"]
            for step in SCRIPT:
                mark = len(c.events)
                c.request("session.prompt", {"sessionId": sid, "text": step["prompt"]})
                if "reply" in step:
                    asked, _ = c.wait_for("permission.asked", after=mark)
                    c.request("permission.reply", {"sessionId": sid,
                                                   "requestId": asked["requestId"],
                                                   **step["reply"]})
                c.wait_for(lambda e: e.get("type") == "turn.ended", after=mark)
                c.wait_for(lambda e: e.get("type") == "session.state"
                           and e.get("state") == "idle", after=mark)
            events = [e for e in c.events if e.get("session") == sid]
            c.request("service.shutdown")
        fixture = {"script": SCRIPT, "session": created["session"],
                   "events": [e for e in events if e.get("type") != "usage"]}
        text = json.dumps(fixture, indent=1, ensure_ascii=False)
        text = text.replace(str(ws), WORKSPACE).replace(sid, "s1")
        data = json.loads(text)
        for e in data["events"]:
            e["ts"] = 0
            if "durationMs" in e:
                e["durationMs"] = 400
        data["session"]["sessionId"] = "s1"
        out = TUI / "test" / "fixtures" / "fake-conversation.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n")
        print(f"wrote {out}: {len(data['events'])} events")
        return 0
    finally:
        service.terminate()
        service.wait(timeout=30)
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.rmtree(runtime, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
