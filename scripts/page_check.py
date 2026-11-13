#!/usr/bin/env python3
"""Real-page check of the dashboard's agent tab, in one command.

Starts the real dashboard (``python -m delfin.cli_voila``) against a
scripted model server, opens it in headless Chromium and records what the
agent tab shows.

Input:
    --out DIR          where record.json, screenshots and logs are written
    --tree PATH        the DELFIN checkout the dashboard runs from
                       (default: the checkout this script lives in)
    --scenario NAME    composer | service-session
    --timeout S        hard deadline for the whole check (SIGALRM)
    --server-python P  interpreter that runs the dashboard and the agent
                       service (default: $DELFIN_PAGE_CHECK_SERVER_PYTHON,
                       else this interpreter)

This script itself must run on an interpreter that has ``playwright`` and
its chromium. The dashboard interpreter needs voila and DELFIN's
dashboard dependencies; the two may be different environments.

Output: ``DIR/record.json`` (also printed), ``DIR/*.png``,
``DIR/logs/{voila,service,model}.log``. Exit code 0 when the scenario's
check passed, 1 when it ran and failed, 2 on a usage error.

Scenarios:
    composer         open the agent tab, select the scripted provider and
                     model, send one prompt, wait for the agent's reply
    service-session  start a private agent service, run one turn in a
                     session there, then open the dashboard and record
                     what the sidebar's running-sessions list shows
    sidebar          a service with a front session, a background run and
                     a saved (ended) session; the agent tab's sidebar and
                     control row photographed at 1366x768 and 1920x1080,
                     the resume field with a typed filter, and the
                     sidebar's text recorded
    keep             the keep switches measured by process: page A runs
                     two sessions, one switched to "Keep session"; the
                     page is closed; the record says whether A's kernel
                     process ended (pid, not id), which session the
                     agent service runs afterwards, and what a newly
                     opened page B lists as running. Then B switches the
                     kept session off and is closed: the session must
                     end with it. On a tree without the switches the
                     same steps run as the control, and the record says
                     which switch was missing

Isolation, each point measured as a failure before it was encoded here:
  * The dashboard's environment is built from nothing (no VSCODE_*,
    BROWSER, TERM_PROGRAM, LD_LIBRARY_PATH, no provider keys of the
    caller); a VS Code terminal's variables made the dashboard open a tab
    on the desktop, and LD_LIBRARY_PATH swapped libstdc++ under sqlite3.
  * Its own HOME, DELFIN_SCRATCH_STATE and workspace: no write reaches the
    caller's ~/.delfin.
  * DELFIN_RUNTIME_DIR is a short 0700 directory under /tmp: AF_UNIX
    paths are limited to 108 bytes, and a longer runtime path breaks the
    agent service's socket without an error on the page.
  * The agent tab is only built when a backend is configured; a
    placeholder ANTHROPIC_API_KEY makes it visible, and the page is then
    switched to the scripted ollama endpoint (OLLAMA_HOST).
  * Cleanup runs in ``finally`` and on the deadline: the dashboard's
    process group is terminated, an agent service on the private socket
    is asked to shut down, any remaining process carrying this run's
    DELFIN_RUNTIME_DIR is killed, and the runtime directory is removed.
    Only processes this check started are touched. An earlier ad-hoc
    check without this left a dashboard running for seven hours.
  * No CDP debugger: Debugger.enable/pause hangs this page.

Rejected: a pytest test with a real browser. The suite's browser tests
are skipped on most machines and take minutes under load; this is a tool
a person or an agent runs on purpose, with its record kept as evidence.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Optional

SCENARIOS = ("composer", "service-session", "sidebar", "keep")

#: The viewports the sidebar scenario photographs: a laptop and a wide
#: screen.
VIEWPORTS = ((1366, 768), (1920, 1080), (1280, 720), (2560, 1440),
             (900, 1400), (390, 844))

#: Measured at every viewport: whether the page itself scrolls (it should
#: fit the window; inner areas scroll), and every pair of visible controls
#: in the sidebar and in the session's control row that overlap.
_LAYOUT_JS = """() => {
  const vis = e => { const r = e.getBoundingClientRect();
    return e.offsetParent !== null && r.width > 0 && r.height > 0 ? r : null; };
  const doc = document.scrollingElement || document.documentElement;
  const cells = document.querySelector('#rendered_cells');
  const out = {page_scroll_y: Math.max(0, doc.scrollHeight - window.innerHeight),
               page_scroll_x: Math.max(0, doc.scrollWidth - window.innerWidth),
               cells_scroll_y: cells ? Math.max(0, cells.scrollHeight - cells.clientHeight) : null,
               overlaps: []};
  const groups = [];
  const side = document.querySelector('.delfin-sessions');
  if (side) {
    groups.push(['sidebar', Array.from(side.children)]);
    const list = side.querySelector('.delfin-session-list');
    if (list) groups.push(['list', Array.from(list.children)]);
  }
  const rows = Array.from(document.querySelectorAll('.delfin-agent-root .widget-hbox'))
    .filter(e => vis(e));
  if (rows.length) groups.push(['controls', Array.from(rows[0].children)]);
  for (const [name, kids] of groups) {
    const boxes = kids.map(k => [k, vis(k)]).filter(p => p[1]);
    for (let i = 0; i < boxes.length; i++) for (let j = i + 1; j < boxes.length; j++) {
      const a = boxes[i][1], b = boxes[j][1];
      const w = Math.min(a.right, b.right) - Math.max(a.left, b.left);
      const h = Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top);
      if (w > 1 && h > 1) out.overlaps.push(name + ': ' +
        (boxes[i][0].innerText || boxes[i][0].className).slice(0, 30) + ' / ' +
        (boxes[j][0].innerText || boxes[j][0].className).slice(0, 30));
    }
  }
  if (side) { const r = side.getBoundingClientRect();
    out.sidebar = [Math.round(r.left), Math.round(r.top), Math.round(r.width), Math.round(r.height)];
    out.sidebar_offscreen_x = r.right > window.innerWidth + 1 || r.left < -1; }
  return out;
}"""

#: The keep scenario's window grace (DELFIN_SESSION_GRACE_SECONDS): short,
#: so a closed page's kernel ends within the check, not after 90 s.
KEEP_GRACE_S = 8
#: How long the keep scenario waits for a kernel process to end, and for a
#: kept session to appear in (or leave) the agent service.
WAIT_KERNEL_END = 180
WAIT_HANDOVER = 180
TOKEN = "delfin-page-check-token-0123456789abcdef"
PLACEHOLDER_KEY = "placeholder-only-shows-the-agent-tab"
PROMPT = "hello from the page check"
SERVER_PYTHON_ENV = "DELFIN_PAGE_CHECK_SERVER_PYTHON"

#: AF_UNIX sun_path is 108 bytes on Linux; leave room for the name.
MAX_SOCKET_PATH = 100

# Waits, in seconds. Generous on purpose: on a shared host at load 300+
# the first render of the dashboard takes minutes.
WAIT_SERVER = 300
WAIT_PAGE = 600
WAIT_COMPOSER = 420
WAIT_REPLY = 180
WAIT_RUNNING = 180
REPLY_SETTLE = 5
SCREENSHOT_TIMEOUT_MS = 120_000


class Deadline(Exception):
    """The --timeout ran out, or the check was asked to stop."""


def default_tree() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="page_check.py",
        description="Real-page check of the dashboard's agent tab.")
    p.add_argument("--out", required=True, type=Path,
                   help="output directory for record.json, screenshots, logs")
    p.add_argument("--tree", type=Path, default=default_tree(),
                   help="DELFIN checkout the dashboard runs from")
    p.add_argument("--scenario", choices=SCENARIOS, default="composer")
    p.add_argument("--timeout", type=int, default=1800,
                   help="hard deadline for the whole check, seconds")
    p.add_argument("--server-python", type=Path,
                   default=Path(os.environ.get(SERVER_PYTHON_ENV)
                                or sys.executable),
                   help=f"interpreter for the dashboard (env {SERVER_PYTHON_ENV})")
    args = p.parse_args(argv)
    if args.timeout <= 0:
        p.error("--timeout must be positive")
    args.tree = args.tree.resolve()
    if not (args.tree / "delfin" / "cli_voila.py").is_file():
        p.error(f"--tree {args.tree} is not a DELFIN checkout")
    args.out = args.out.resolve()
    return args


def make_runtime_dir(base: str = "/tmp") -> Path:
    """A short, private DELFIN_RUNTIME_DIR; refused if the socket path
    inside it would not fit in an AF_UNIX address."""
    path = Path(tempfile.mkdtemp(prefix="dpc-", dir=base))
    os.chmod(path, 0o700)
    sock = service_socket(path)
    if len(str(sock).encode()) > MAX_SOCKET_PATH:
        shutil.rmtree(path, ignore_errors=True)
        raise ValueError(f"socket path too long ({len(str(sock))} bytes): {sock}")
    return path


def service_socket(runtime: Path) -> Path:
    """Where the agent service listens for this runtime dir
    (delfin/agent/service/rundir.py: socket_path)."""
    return Path(runtime) / f"delfin-{os.getuid()}" / "agent.sock"


def server_env(*, server_python: Path, tree: Path, home: Path, scratch: Path,
               runtime: Path, workspace: Path, model_port: int) -> dict:
    """The dashboard's whole environment, built from nothing.

    Nothing of the caller's environment is inherited: no VS Code or
    browser hooks, no loader path, no provider credentials."""
    return {
        "PATH": f"{Path(server_python).parent}:/usr/bin:/bin",
        "LANG": "C.UTF-8",
        "HOME": str(home),
        "PYTHONPATH": str(tree),
        "DELFIN_SCRATCH_STATE": str(scratch),
        "DELFIN_RUNTIME_DIR": str(runtime),
        "DELFIN_AGENT_WORKSPACE": str(workspace),
        "DELFIN_NO_BROWSER": "1",
        "DELFIN_CLASSIFIER": "off",
        "JUPYTER_TOKEN": TOKEN,
        # Sessions start in the workspace, not in the checkout the server
        # runs from: a second session there was given a worktree of the
        # checkout under test.
        "DELFIN_LAUNCH_CWD": str(workspace),
        "OLLAMA_HOST": f"http://127.0.0.1:{model_port}",
        "PYTHON_KEYRING_BACKEND": "keyring.backends.null.Keyring",
        "ANTHROPIC_API_KEY": PLACEHOLDER_KEY,
    }


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def _pids_with_env(name: str, value: str) -> list[int]:
    """This user's processes whose environment has name=value."""
    needle = f"{name}={value}".encode()
    found = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit() or int(entry) == os.getpid():
            continue
        try:
            if os.stat(f"/proc/{entry}").st_uid != os.getuid():
                continue
            with open(f"/proc/{entry}/environ", "rb") as fh:
                env = fh.read().split(b"\0")
        except OSError:
            continue
        if needle in env:
            found.append(int(entry))
    return found


def _in_main_thread() -> bool:
    import threading
    return threading.current_thread() is threading.main_thread()


def _signal_all(pids: list[int], sig: int) -> None:
    for pid in pids:
        try:
            os.kill(pid, sig)
        except ProcessLookupError:
            pass


class Cleanup:
    """Everything this check started, and how to end it. Idempotent."""

    def __init__(self, grace: float = 20.0) -> None:
        #: seconds between SIGTERM and SIGKILL
        self.grace = grace
        self.groups: list[subprocess.Popen] = []
        self.runtime: Path | None = None
        self.work: Path | None = None
        self.model_server = None
        self.tree: Path | None = None
        self.report: dict = {}
        self._done = False

    def add_group(self, proc: subprocess.Popen) -> None:
        self.groups.append(proc)

    def _shutdown_service(self) -> str:
        if self.runtime is None:
            return "none"
        sock = service_socket(self.runtime)
        if not sock.exists():
            return "none"
        try:
            if self.tree is not None and str(self.tree) not in sys.path:
                sys.path.insert(0, str(self.tree))
            from delfin.agent.service.client import ServiceClient
            with ServiceClient(sock, timeout=5.0, client="page-check") as c:
                c.request("service.shutdown", timeout=20)
            return "shutdown"
        except Exception as exc:  # the sweep below still ends it
            return f"shutdown failed: {type(exc).__name__}"

    def _end_group(self, proc: subprocess.Popen) -> str:
        if proc.poll() is not None:
            return "exited"
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            return "exited"
        try:
            proc.wait(timeout=self.grace)
            return "terminated"
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait(timeout=10)
            return "killed"

    def run(self) -> dict:
        if self._done:
            return self.report
        self._done = True
        signal.alarm(0)
        # A second Ctrl-C must not abandon the cleanup half way; it is
        # bounded (about a minute) and SIGKILL still ends this process.
        if _in_main_thread():
            for sig in (signal.SIGTERM, signal.SIGINT):
                signal.signal(sig, signal.SIG_IGN)
        rep = self.report
        rep["service"] = self._shutdown_service()
        rep["groups"] = [self._end_group(p) for p in self.groups]
        leftover: list[int] = []
        if self.runtime is not None:
            # Kernels and services leave the dashboard's process group;
            # they still carry this run's runtime dir in their environment.
            leftover = _pids_with_env("DELFIN_RUNTIME_DIR", str(self.runtime))
            _signal_all(leftover, signal.SIGTERM)
            deadline = time.time() + self.grace
            while time.time() < deadline and _pids_with_env(
                    "DELFIN_RUNTIME_DIR", str(self.runtime)):
                time.sleep(0.5)
            _signal_all(_pids_with_env("DELFIN_RUNTIME_DIR", str(self.runtime)),
                        signal.SIGKILL)
            shutil.rmtree(self.runtime, ignore_errors=True)
            rep["runtime_dir_removed"] = not self.runtime.exists()
        rep["leftovers_ended"] = len(leftover)
        if self.model_server is not None:
            self.model_server.shutdown()
        if self.work is not None:
            shutil.rmtree(self.work, ignore_errors=True)
        return rep


_SIGNALS = (signal.SIGALRM, signal.SIGTERM, signal.SIGINT)


def install_deadline(seconds: int) -> dict:
    """SIGALRM after ``seconds``, SIGTERM and SIGINT: each raises Deadline
    in the main thread, so the ``finally`` with the cleanup runs. Returns
    the previous handlers for :func:`restore_signals`."""
    def _raise(signum, _frame):
        raise Deadline(f"stopped by signal {signal.Signals(signum).name}")
    previous = {sig: signal.getsignal(sig) for sig in _SIGNALS}
    for sig in _SIGNALS:
        signal.signal(sig, _raise)
    signal.alarm(seconds)
    return previous


def restore_signals(previous: dict) -> None:
    signal.alarm(0)
    for sig, handler in previous.items():
        signal.signal(sig, handler)


def _git_head(tree: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(tree), "rev-parse", "HEAD"],
                              capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


def _start_service(args, env, work, cleanup: Cleanup, rec: dict, model_name: str,
                   *, extra: bool = False):
    """A private agent service with one session that ran a turn. ``extra``:
    also a background run (a session in a worktree folder) and a saved
    session that ended, for the sidebar's groups and the resume list."""
    from delfin.agent.service.client import ServiceClient
    logs = args.out / "logs"
    proc = subprocess.Popen(
        [str(args.server_python), "-m", "delfin.agent.service.supervisor",
         "--idle-exit", str(args.timeout)],
        cwd=str(work / "ws"), env=env, stdout=open(logs / "service.log", "wb"),
        stderr=subprocess.STDOUT, start_new_session=True)
    cleanup.add_group(proc)
    sock = service_socket(cleanup.runtime)
    t0 = time.time()
    while not sock.exists():
        if proc.poll() is not None:
            raise RuntimeError(f"agent service exited with {proc.returncode}")
        if time.time() - t0 > WAIT_SERVER:
            raise RuntimeError("agent service socket did not appear")
        time.sleep(0.5)
    with ServiceClient(sock, client="page-check", timeout=30.0) as c:
        sid = c.request("session.create", {
            "cwd": str(work / "ws"), "provider": "ollama", "model": model_name,
            "attach": True}, timeout=180)["sessionId"]
        start = len(c.events)
        c.request("session.prompt", {"sessionId": sid, "text": PROMPT})
        c.wait_for(lambda e: e.get("type") == "turn.ended"
                   and e.get("session") == sid, after=start, timeout=WAIT_REPLY)
        rec["service_session"] = sid
        if extra:
            background = work / "ws" / ".delfin" / "worktrees" / "delfin-wt-0a1b2c3d"
            background.mkdir(parents=True, exist_ok=True)
            for cwd, text, end in ((background, "background run: nightly check", False),
                                   (work / "ws", "saved session: earlier work", True)):
                other = c.request("session.create", {
                    "cwd": str(cwd), "provider": "ollama", "model": model_name,
                    "attach": True}, timeout=180)["sessionId"]
                start = len(c.events)
                c.request("session.prompt", {"sessionId": other, "text": text})
                c.wait_for(lambda e, _o=other: e.get("type") == "turn.ended"
                           and e.get("session") == _o, after=start, timeout=WAIT_REPLY)
                if end:
                    c.request("session.shutdown", {"sessionId": other}, timeout=60)
            time.sleep(3)
        rec["service_sessions"] = [s["sessionId"]
                                   for s in c.request("sessions.list")["sessions"]]
    rec["timings_s"]["service_turn"] = round(time.time() - t0, 1)


def _shot(page, args, name: str, rec: dict) -> None:
    path = args.out / f"{name}.png"
    try:
        page.screenshot(path=str(path), timeout=SCREENSHOT_TIMEOUT_MS)
        rec["screenshots"].append(path.name)
    except Exception as exc:
        rec.setdefault("screenshot_errors", []).append(f"{name}: {type(exc).__name__}")


def _visible(page, selector: str) -> int:
    return page.evaluate(
        "sel => Array.from(document.querySelectorAll(sel))"
        ".filter(e => e.offsetParent !== null).length", selector)


def _text(page, selector: str) -> str:
    return page.evaluate(
        "sel => Array.from(document.querySelectorAll(sel))"
        ".map(e => e.innerText).join(' | ')", selector)


def _open_agent_tab(browser, args, url: str, rec: dict, *, name: str = "",
                    viewport: tuple = (1400, 1000)):
    """A new page on the dashboard with its agent tab open and the composer
    shown. Each page is a kernel of its own. Returns (page, errors)."""
    t = rec["timings_s"]
    t0 = time.time()
    page = browser.new_page(viewport={"width": viewport[0], "height": viewport[1]})
    errors: list[str] = []
    page.on("pageerror", lambda e: errors.append(str(e)[:200]))
    # "commit", not "domcontentloaded": voila holds the document
    # open while the kernel renders, and a render that never ends
    # should still leave a screenshot of how far it came.
    page.goto(f"{url}/?token={TOKEN}", wait_until="commit",
              timeout=WAIT_PAGE * 1000)
    # Under a load of several hundred Voila sometimes sends a widget's view
    # before its model ("Error displaying widget: model not found", seen on
    # origin/main as on a branch): that page never renders, and a new one
    # (a new kernel) is the way on. Twice at most.
    deadline, reloads = time.time() + WAIT_PAGE, 0
    while True:
        if _visible(page, ".delfin-session-strip, .delfin-sessions"):
            break
        broken = "model not found" in (page.evaluate(
            "() => document.body ? document.body.innerText.slice(0, 2000) : ''") or "")
        if time.time() > deadline or (broken and reloads >= 2):
            _shot(page, args, f"not_rendered{name}", rec)
            rec["page_errors"] = errors[:10]
            raise RuntimeError("the dashboard did not render"
                               + (" (model not found)" if broken else ""))
        if broken:
            reloads += 1
            rec.setdefault("reloads", []).append(name or "page")
            time.sleep(5)
            page.goto(f"{url}/?token={TOKEN}", wait_until="commit",
                      timeout=WAIT_PAGE * 1000)
            continue
        time.sleep(3)
    t[f"page_loaded{name}"] = round(time.time() - t0, 1)
    tabs = page.evaluate("() => Array.from(document.querySelectorAll("
                         "'.lm-TabBar-tab')).map(t => t.innerText.trim())")
    rec["tabs"] = tabs
    target = next((x for x in tabs if "agent" in x.lower()), None)
    rec["agent_tab"] = target
    if target is None:
        _shot(page, args, f"no_agent_tab{name}", rec)
        raise RuntimeError("no agent tab on the page")
    # The tab is built on first open.
    page.locator(".lm-TabBar-tab").filter(has_text=target).first.click()
    try:
        # The agent tab's own composer, or a service session's view.
        page.wait_for_selector(".delfin-agent-input textarea, .delfin-live-input textarea",
                               timeout=WAIT_COMPOSER * 1000)
    except Exception:
        _shot(page, args, f"stuck{name}", rec)
        raise
    t[f"composer{name}"] = round(time.time() - t0, 1)
    time.sleep(3)
    return page, errors


def _drive_page(args, url: str, rec: dict, model_name: str, env: dict,
                cleanup: "Cleanup") -> None:
    from playwright.sync_api import sync_playwright
    t0 = time.time()
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        try:
            if args.scenario == "keep":
                _keep_check(browser, args, url, rec, model_name, env, cleanup)
                return
            viewport = VIEWPORTS[0] if args.scenario == "sidebar" else (1400, 1000)
            page, errors = _open_agent_tab(browser, args, url, rec, viewport=viewport)
            rec["composer_visible"] = _visible(page, ".delfin-agent-input textarea")
            _shot(page, args, "opened", rec)
            if args.scenario == "composer":
                _send_prompt(page, args, rec, model_name, t0)
            elif args.scenario == "sidebar":
                _photograph_sidebar(page, args, rec)
            else:
                _read_running(page, args, rec, t0)
            rec["page_errors"] = errors[:10]
        finally:
            browser.close()


#: The element a selector finds that is on screen (several sessions share
#: the page and only the active one is visible).
_VISIBLE_JS = """(sel) => Array.from(document.querySelectorAll(sel))
    .filter(e => e.offsetParent !== null && e.getBoundingClientRect().width > 0)"""


def _pick_visible_model(page, model_name: str) -> dict:
    """Select the scripted provider and model in the visible session's
    dropdowns; returns what was picked."""
    picked = {}
    for label, want in (("Provider:", "ollama"), ("Model:", model_name)):
        got = page.evaluate(
            "([label, want]) => {"
            " const boxes = Array.from(document.querySelectorAll('.widget-dropdown'))"
            "   .filter(b => b.offsetParent !== null && b.innerText.includes(label));"
            " for (const b of boxes) { const s = b.querySelector('select'); if (!s) continue;"
            "   const o = Array.from(s.options).find(o => o.text.toLowerCase().includes(want.toLowerCase()));"
            "   if (!o) return null;"
            "   if (s.value !== o.value) { s.value = o.value;"
            "     s.dispatchEvent(new Event('change', {bubbles: true})); }"
            "   return o.text; }"
            " return null; }", [label, want])
        picked[label.rstrip(":")] = got
        time.sleep(3)
    return picked


def _send_visible(page, text: str, *, wait_s: float = WAIT_REPLY) -> str:
    """Type *text* into the visible composer, send it, and return the
    visible chat's last agent reply once it has settled."""
    page.evaluate(
        "(text) => { const ta = (" + _VISIBLE_JS + ")('.delfin-agent-input textarea')[0];"
        " ta.focus(); ta.value = text;"
        " ta.dispatchEvent(new Event('input', {bubbles: true})); }", text)
    time.sleep(1)
    page.evaluate("() => { const b = (" + _VISIBLE_JS + ")('.delfin-agent-send-row button')[0];"
                  " b.click(); }")
    t_send = time.time()
    reply, settled_at = "", None
    while time.time() - t_send < wait_s:
        now = page.evaluate(
            "() => { const r = (" + _VISIBLE_JS + ")('.delfin-chat-agent');"
            " return r.length ? r[r.length - 1].innerText.trim() : ''; }")
        if now and now == reply and "received" in now:
            settled_at = settled_at or time.time()
            if time.time() - settled_at >= REPLY_SETTLE:
                break
        else:
            reply, settled_at = now, None
        time.sleep(1)
    return reply[:300]


def _sidebar_rows(page) -> list[dict]:
    """Every row of the session sidebar, in order: its text and classes."""
    return page.evaluate(
        "() => Array.from(document.querySelectorAll('.delfin-sessions "
        ".delfin-session-item, .delfin-sessions .delfin-running-item, "
        ".delfin-sessions .delfin-resume-row, .delfin-sessions .delfin-list-divider'))"
        ".filter(e => e.offsetParent !== null)"
        ".map(e => ({text: e.innerText.trim().replace(/\\s+/g, ' ').slice(0, 160),"
        " classes: Array.from(e.classList).filter(c => c.startsWith('delfin-'))}))")


def _switches(page) -> list[dict]:
    """The visible switches and toggles: header and control row."""
    return page.evaluate(
        "() => Array.from(document.querySelectorAll('.delfin-switch, "
        ".delfin-keep-session, .delfin-keep-dashboard'))"
        ".filter(e => e.offsetParent !== null)"
        ".map(e => ({text: e.innerText.trim(), title: e.title || '',"
        " on: e.classList.contains('mod-active'),"
        " classes: Array.from(e.classList).filter(c => c.startsWith('delfin-'))}))")


def _clip(page, args, name: str, rec: dict, selector: str) -> None:
    """A screenshot of the first visible element matching *selector*."""
    path = args.out / f"{name}.png"
    try:
        handle = page.evaluate_handle(
            "(sel) => (" + _VISIBLE_JS + ")(sel)[0] || null", selector)
        element = handle.as_element()
        if element is None:
            rec.setdefault("screenshot_errors", []).append(f"{name}: no {selector}")
            return
        element.screenshot(path=str(path), timeout=SCREENSHOT_TIMEOUT_MS)
        rec["screenshots"].append(path.name)
    except Exception as exc:
        rec.setdefault("screenshot_errors", []).append(f"{name}: {type(exc).__name__}")


def _photograph_sidebar(page, args, rec: dict) -> None:
    """The sidebar and the control row at both viewports, and the resume
    field with a filter typed into it."""
    shots = rec.setdefault("sidebar", {})
    for width, height in VIEWPORTS:
        tag = f"{width}x{height}"
        page.set_viewport_size({"width": width, "height": height})
        time.sleep(4)
        _shot(page, args, f"page_{tag}", rec)
        _clip(page, args, f"sidebar_{tag}", rec, ".delfin-sessions")
        _clip(page, args, f"controls_{tag}", rec, ".delfin-agent-root .widget-hbox")
        shots[tag] = {"rows": _sidebar_rows(page), "switches": _switches(page),
                      "layout": page.evaluate(_LAYOUT_JS)}
    page.set_viewport_size({"width": VIEWPORTS[0][0], "height": VIEWPORTS[0][1]})
    time.sleep(2)
    field = ".delfin-session-resume input"
    if _visible(page, field):
        page.locator(field).first.click()
        time.sleep(2)
        _shot(page, args, "resume_open", rec)
        page.locator(field).first.type("saved", delay=80)
        time.sleep(4)
        _shot(page, args, "resume_typed", rec)
        shots["resume_typed_rows"] = _sidebar_rows(page)
        shown = ("() => document.querySelectorAll('.delfin-resume-pop.delfin-resume-shown')"
                 ".length")
        shots["list_shown_typed"] = page.evaluate(shown)
        page.keyboard.press("ArrowDown")
        time.sleep(1)
        shots["active_after_arrow"] = page.evaluate(
            "() => Array.from(document.querySelectorAll('.delfin-resume-active'))"
            ".map(e => e.innerText.trim())")
        _shot(page, args, "resume_arrow", rec)
        page.keyboard.press("Escape")
        time.sleep(2)
        shots["list_shown_after_escape"] = page.evaluate(shown)
        shots["after_escape_rows"] = _sidebar_rows(page)
        # Enter opens the first offer, and the list closes.
        page.locator(field).first.click()
        page.locator(field).first.fill("")
        page.locator(field).first.type("saved", delay=80)
        time.sleep(4)
        page.keyboard.press("Enter")
        time.sleep(15)
        shots["list_shown_after_enter"] = page.evaluate(shown)
        shots["field_after_enter"] = page.locator(field).first.input_value()
        shots["rows_after_enter"] = _sidebar_rows(page)
        _shot(page, args, "resume_entered", rec)
    rec["sidebar_text"] = page.evaluate(
        "() => ((document.querySelector('.delfin-sessions')||{}).innerText||'')")[:1500]
    rec["ok"] = bool(shots.get(f"{VIEWPORTS[0][0]}x{VIEWPORTS[0][1]}", {}).get("rows"))


# -- the keep check ----------------------------------------------------------

def _kernel_pids(url: str, *, watched: bool = False) -> dict:
    """kernel id -> pid of the process running it (read from /proc by the
    connection file in its command line), for every kernel the server
    lists; ``watched``: only those a window is connected to (a page that
    did not render and was reloaded leaves one without)."""
    try:
        with urllib.request.urlopen(f"{url}/api/kernels?token={TOKEN}", timeout=10) as r:
            ids = [k["id"] for k in json.loads(r.read().decode())
                   if not watched or int(k.get("connections") or 0) > 0]
    except Exception:
        ids = []
    out = {}
    for kid in ids:
        out[kid] = _pid_of_kernel(kid)
    return out


def _pid_of_kernel(kid: str) -> int:
    needle = f"kernel-{kid}.json".encode()
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/cmdline", "rb") as fh:
                if needle in fh.read():
                    return int(entry)
        except OSError:
            continue
    return 0


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        with open(f"/proc/{pid}/stat", "rb") as fh:
            # A zombie has ended; only its parent has not read it yet.
            return fh.read().split(b") ", 1)[1][:1] != b"Z"
    except (OSError, IndexError):
        return False


def _saved_sessions(args, env: dict) -> list[dict]:
    """The session store of the dashboard's HOME, read by the dashboard's
    interpreter: id, title, workspace."""
    code = ("import json; from delfin.agent.session_store import list_sessions; "
            "print(json.dumps([{k: r.get(k) for k in ('session_id', 'title', 'workspace')}"
            " for r in list_sessions(50)]))")
    try:
        out = subprocess.run([str(args.server_python), "-c", code], env=env,
                             cwd=str(args.tree), capture_output=True, text=True,
                             timeout=120)
        return json.loads(out.stdout.strip().splitlines()[-1])
    except Exception:
        return []


def _service_sessions(runtime: Path, tree: Path) -> list[dict]:
    """The agent service's sessions (sessionId, state, keep), or [] when no
    service listens on this run's socket."""
    sock = service_socket(runtime)
    if not sock.exists():
        return []
    try:
        if str(tree) not in sys.path:
            sys.path.insert(0, str(tree))
        from delfin.agent.service.client import ServiceClient
        with ServiceClient(sock, timeout=10.0, client="page-check") as c:
            rows = c.request("sessions.list", timeout=30)["sessions"]
        return [{"sessionId": s.get("sessionId"), "state": s.get("state"),
                 "keep": s.get("keep"), "title": s.get("title")} for s in rows]
    except Exception as exc:
        return [{"error": f"{type(exc).__name__}: {exc}"[:200]}]


def _kept_records(env: dict) -> dict:
    """What the kept-session records and the hand-over log say (the
    dashboard's HOME and its scratch state): file name -> text."""
    out = {}
    for root in (env.get("HOME"), env.get("DELFIN_SCRATCH_STATE")):
        folder = Path(str(root or "")) / ".delfin" / "kept_agent_sessions"
        if not root or not folder.is_dir():
            continue
        for path in sorted(folder.iterdir()):
            try:
                out[str(path.relative_to(Path(root)))] = path.read_text(
                    encoding="utf-8", errors="replace")[-2000:]
            except OSError:
                continue
    return out


def _click_row(page, text: str) -> bool:
    """Click the sidebar row whose title contains *text*."""
    return bool(page.evaluate(
        "(text) => { const rows = Array.from(document.querySelectorAll("
        "'.delfin-sessions .delfin-session-title')).filter(e => e.offsetParent !== null"
        " && e.innerText.includes(text)); if (!rows.length) return false;"
        " rows[0].click(); return true; }", text))


def _toggle_keep(page) -> dict:
    """Click the visible "Keep session" switch; its state before and after,
    or {"present": False}."""
    sel = ".delfin-keep-session"
    if not _visible(page, sel):
        return {"present": False}
    before = page.evaluate("(sel) => (" + _VISIBLE_JS + ")(sel)[0].classList.contains('mod-active')", sel)
    page.evaluate("(sel) => (" + _VISIBLE_JS + ")(sel)[0].click()", sel)
    time.sleep(4)
    after = page.evaluate("(sel) => (" + _VISIBLE_JS + ")(sel)[0].classList.contains('mod-active')", sel)
    return {"present": True, "before": before, "after": after}


def _close_page(page) -> float:
    """Leave the page the way a closed browser tab does (pagehide, the
    page's close beacon), then close it. Returns when it was closed."""
    try:
        page.goto("about:blank", timeout=30_000)
    except Exception:
        pass
    page.close()
    return time.time()


def _wait_pid_gone(pid: int, timeout: float) -> Optional[float]:
    t0 = time.time()
    while time.time() - t0 < timeout:
        if not _pid_alive(pid):
            return round(time.time() - t0, 1)
        time.sleep(1)
    return None


def _keep_check(browser, args, url: str, rec: dict, model_name: str, env: dict,
                cleanup: "Cleanup") -> None:
    """The keep scenario (module docstring). Every fact by process or by
    the service's own list, never by a widget's word."""
    keep = rec.setdefault("keep", {})
    runtime = Path(env["DELFIN_RUNTIME_DIR"])
    # -- page A: two sessions, the first one kept
    page, errors = _open_agent_tab(browser, args, url, rec, name="_a")
    keep["picked_a"] = _pick_visible_model(page, model_name)
    keep["reply_keep"] = _send_visible(page, "keep-me: hello")
    page.locator(".delfin-session-add").first.click()
    time.sleep(2)
    # A folder of its own: a second session in the first one's folder is
    # refused (one folder, one live session).
    second = Path(env["DELFIN_AGENT_WORKSPACE"]).parent / "ws2"
    second.mkdir(exist_ok=True)
    page.locator(".delfin-path-input input").first.fill(str(second))
    time.sleep(3)
    page.locator(".delfin-session-start").first.click()
    time.sleep(8)
    keep["picked_a2"] = _pick_visible_model(page, model_name)
    keep["reply_end"] = _send_visible(page, "end-me: hello")
    keep["clicked_keep_row"] = _click_row(page, "keep-me")
    time.sleep(3)
    keep["switch_a"] = _toggle_keep(page)
    keep["switches_a"] = _switches(page)
    keep["rows_a"] = _sidebar_rows(page)
    _shot(page, args, "keep_page_a", rec)
    _clip(page, args, "keep_sidebar_a", rec, ".delfin-sessions")
    saved = _saved_sessions(args, env)
    keep["saved"] = saved
    ids = {("keep" if "keep-me" in str(s.get("title")) else
            "end" if "end-me" in str(s.get("title")) else "other"): s["session_id"]
           for s in saved}
    keep["ids"] = ids
    kernels = _kernel_pids(url, watched=True)
    keep["kernels_a"] = kernels
    pid_a = next(iter(kernels.values()), 0) if len(kernels) == 1 else 0
    keep["kernel_a_pid"] = pid_a
    keep["page_errors_a"] = errors[:10]
    # -- close A
    _close_page(page)
    keep["kernel_a_ended_after_s"] = _wait_pid_gone(pid_a, WAIT_KERNEL_END)
    t0 = time.time()
    service: list = []
    while time.time() - t0 < WAIT_HANDOVER:
        service = _service_sessions(runtime, args.tree)
        if any(s.get("sessionId") == ids.get("keep") for s in service):
            break
        time.sleep(3)
    keep["service_after_a"] = service
    keep["kept_records_after_a"] = _kept_records(env)
    keep["kept_runs_after_a"] = any(s.get("sessionId") == ids.get("keep") for s in service)
    keep["ended_is_gone_after_a"] = not any(s.get("sessionId") == ids.get("end")
                                            for s in service)
    # -- page B: is it listed as running? A dashboard server whose last
    # window closed with nothing kept stops (resume_server); a person then
    # starts it again, and so does this check.
    server = rec["_servers"][-1]
    time.sleep(5)
    keep["server_stopped_after_a"] = server.poll() is not None
    if keep["server_stopped_after_a"]:
        url = start_dashboard(args, env, cleanup, rec, log="voila_b.log")
    page, errors = _open_agent_tab(browser, args, url, rec, name="_b")
    time.sleep(8)
    keep["rows_b"] = _sidebar_rows(page)
    keep["switches_b"] = _switches(page)
    _shot(page, args, "keep_page_b", rec)
    _clip(page, args, "keep_sidebar_b", rec, ".delfin-sessions")
    keep["listed_running_in_b"] = any(
        "keep-me" in r["text"] and "delfin-state-running" in r["classes"]
        for r in keep["rows_b"])
    # -- B switches it off and closes: it must end with B
    if keep["kept_runs_after_a"]:
        keep["opened_in_b"] = _click_row(page, "keep-me")
        time.sleep(10)
        keep["switch_b"] = _toggle_keep(page)
        keep["rows_b_after_off"] = _sidebar_rows(page)
        _shot(page, args, "keep_page_b_off", rec)
    kernels_b = {k: v for k, v in _kernel_pids(url, watched=True).items()
                 if k not in kernels}
    keep["kernels_b"] = kernels_b
    pid_b = next(iter(kernels_b.values()), 0)
    keep["page_errors_b"] = errors[:10]
    _close_page(page)
    keep["kernel_b_ended_after_s"] = _wait_pid_gone(pid_b, WAIT_KERNEL_END)
    if keep["kept_runs_after_a"]:
        time.sleep(20)
        service = _service_sessions(runtime, args.tree)
        keep["service_after_b"] = service
        keep["off_session_ended_with_b"] = not any(
            s.get("sessionId") == ids.get("keep") for s in service)
    rec["ok"] = bool(keep["kernel_a_ended_after_s"] is not None
                     and keep["kept_runs_after_a"] and keep["ended_is_gone_after_a"]
                     and keep["listed_running_in_b"]
                     and keep.get("off_session_ended_with_b"))


def _send_prompt(page, args, rec, model_name: str, t0: float) -> None:
    picked = rec.setdefault("picked", {})
    for label, want in (("Provider:", "ollama"), ("Model:", model_name)):
        sel = page.locator(".widget-dropdown", has_text=label).locator("select").first
        opts = sel.evaluate("s => Array.from(s.options).map(o => o.text)")
        choice = next((o for o in opts if want.lower() in o.lower()), None)
        picked[label.rstrip(":")] = choice
        if choice:
            sel.select_option(label=choice)
            time.sleep(3)
    page.locator(".delfin-agent-input textarea").first.fill(PROMPT)
    t_send = time.time()
    page.locator(".delfin-agent-send-row button").first.click()
    # The reply streams in: take it once it has not changed for
    # REPLY_SETTLE seconds, and time its first appearance separately.
    reply, first, settled_at = "", None, None
    while time.time() - t_send < WAIT_REPLY:
        now = _text(page, ".delfin-chat-agent").strip()
        if now and first is None:
            first = time.time()
        if now and now == reply:
            settled_at = settled_at or time.time()
            if time.time() - settled_at >= REPLY_SETTLE:
                break
        else:
            reply, settled_at = now, None
        time.sleep(1)
    rec["reply_text"] = reply[:400]
    if first is not None:
        rec["timings_s"]["first_reply_text"] = round(first - t_send, 1)
    rec["timings_s"]["reply_settled"] = round(time.time() - t_send, 1)
    rec["timings_s"]["total_page"] = round(time.time() - t0, 1)
    _shot(page, args, "answered", rec)
    rec["ok"] = bool(rec.get("composer_visible")) and "received" in reply


def _read_running(page, args, rec, t0: float) -> None:
    t_wait = time.time()
    try:
        page.wait_for_selector(".delfin-running-item", timeout=WAIT_RUNNING * 1000)
    except Exception:
        pass
    rec["timings_s"]["running_list"] = round(time.time() - t_wait, 1)
    rec["running_items"] = page.evaluate(
        "() => Array.from(document.querySelectorAll('.delfin-running-item'))"
        ".map(e => e.innerText.trim().slice(0, 200))")
    rec["sidebar_text"] = page.evaluate(
        "() => ((document.querySelector('.delfin-sessions')||{}).innerText||'')"
    )[:800]
    rec["timings_s"]["total_page"] = round(time.time() - t0, 1)
    _shot(page, args, "running", rec)
    rec["ok"] = len(rec["running_items"]) >= 1


def start_dashboard(args, env: dict, cleanup: "Cleanup", rec: dict, *,
                    log: str = "voila.log") -> str:
    """Start the dashboard server; return its URL once it answers."""
    port = free_port()
    dash = subprocess.Popen(
        [str(args.server_python), "-m", "delfin.cli_voila", "--port", str(port),
         "--ip", "127.0.0.1", "--token", TOKEN, "--no-browser"],
        cwd=str(args.tree), env=env,
        stdout=open(args.out / "logs" / log, "wb"),
        stderr=subprocess.STDOUT, start_new_session=True)
    cleanup.add_group(dash)
    url = f"http://127.0.0.1:{port}"
    t_srv = time.time()
    while True:
        if dash.poll() is not None:
            raise RuntimeError(f"dashboard exited with {dash.returncode}")
        try:
            urllib.request.urlopen(f"{url}/api/kernels?token={TOKEN}", timeout=5)
            break
        except Exception:
            if time.time() - t_srv > WAIT_SERVER:
                raise RuntimeError("dashboard server did not answer")
            time.sleep(2)
    rec["timings_s"].setdefault("server_ready", round(time.time() - t_srv, 1))
    rec.setdefault("_servers", []).append(dash)
    return url


def run_check(args: argparse.Namespace) -> dict:
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "logs").mkdir(exist_ok=True)
    rec: dict = {"scenario": args.scenario, "tree_head": _git_head(args.tree),
                 "started": time.strftime("%Y-%m-%dT%H:%M:%S"), "ok": False,
                 "timings_s": {}, "screenshots": []}
    cleanup = Cleanup()
    cleanup.tree = args.tree
    t0 = time.time()
    previous = install_deadline(args.timeout)
    try:
        sys.path.insert(0, str(args.tree / "tests"))
        sys.path.insert(0, str(args.tree))
        from agent_terminal_probes import fake_model
        work = Path(tempfile.mkdtemp(prefix="work-", dir=str(args.out)))
        cleanup.work = work
        cleanup.runtime = make_runtime_dir()
        home, scratch, ws = work / "home", work / "scratch", work / "ws"
        home.mkdir()
        (scratch / ".delfin").mkdir(parents=True)
        ws.mkdir()
        (home / ".delfin_settings.json").write_text(json.dumps({"agent": {
            "provider": "ollama", "model": fake_model.MODEL,
            "permission_profile": "all_free", "probe_endpoint": False}}))
        srv, model_port = fake_model.serve(
            log_path=str(args.out / "logs" / "model.log"))
        cleanup.model_server = srv
        # Scratch state at the run's own HOME: a dashboard kernel reads
        # ~/.delfin and an agent worker honours DELFIN_SCRATCH_STATE, and two
        # different folders were two session stores -- a session saved in
        # the page was "not saved" for the service.
        env = server_env(server_python=args.server_python, tree=args.tree,
                         home=home, scratch=home, runtime=cleanup.runtime,
                         workspace=ws, model_port=model_port)
        if args.scenario == "keep":
            env["DELFIN_SESSION_GRACE_SECONDS"] = str(KEEP_GRACE_S)
        if args.scenario == "service-session":
            _start_service(args, env, work, cleanup, rec, fake_model.MODEL)
        elif args.scenario == "sidebar":
            _start_service(args, env, work, cleanup, rec, fake_model.MODEL,
                           extra=True)
        url = start_dashboard(args, env, cleanup, rec)
        _drive_page(args, url, rec, fake_model.MODEL, env, cleanup)
    except BaseException as exc:  # noqa: BLE001 -- recorded, then cleaned up
        rec["error"] = f"{type(exc).__name__}: {exc}"[:500]
        rec["ok"] = False
    finally:
        rec.pop("_servers", None)
        rec["cleanup"] = cleanup.run()
        restore_signals(previous)
        rec["timings_s"]["wall"] = round(time.time() - t0, 1)
        (args.out / "record.json").write_text(json.dumps(rec, indent=1))
    return rec


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rec = run_check(args)
    print(json.dumps(rec, indent=1))
    return 0 if rec.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
