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

SCENARIOS = ("composer", "service-session")
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


def _start_service(args, env, work, cleanup: Cleanup, rec: dict, model_name: str):
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


def _drive_page(args, url: str, rec: dict, model_name: str) -> None:
    from playwright.sync_api import sync_playwright
    t = rec["timings_s"]
    t0 = time.time()
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        try:
            page = browser.new_page(viewport={"width": 1400, "height": 1000})
            errors: list[str] = []
            page.on("pageerror", lambda e: errors.append(str(e)[:200]))
            # "commit", not "domcontentloaded": voila holds the document
            # open while the kernel renders, and a render that never ends
            # should still leave a screenshot of how far it came.
            page.goto(f"{url}/?token={TOKEN}", wait_until="commit",
                      timeout=WAIT_PAGE * 1000)
            try:
                page.wait_for_selector(".delfin-session-strip, .delfin-sessions",
                                       timeout=WAIT_PAGE * 1000)
            except Exception:
                _shot(page, args, "not_rendered", rec)
                rec["page_errors"] = errors[:10]
                raise
            t["page_loaded"] = round(time.time() - t0, 1)
            tabs = page.evaluate("() => Array.from(document.querySelectorAll("
                                 "'.lm-TabBar-tab')).map(t => t.innerText.trim())")
            rec["tabs"] = tabs
            target = next((x for x in tabs if "agent" in x.lower()), None)
            rec["agent_tab"] = target
            if target is None:
                _shot(page, args, "no_agent_tab", rec)
                raise RuntimeError("no agent tab on the page")
            # The tab is built on first open.
            page.locator(".lm-TabBar-tab").filter(has_text=target).first.click()
            try:
                page.wait_for_selector(".delfin-agent-input textarea",
                                       timeout=WAIT_COMPOSER * 1000)
            except Exception:
                _shot(page, args, "stuck", rec)
                raise
            t["composer"] = round(time.time() - t0, 1)
            time.sleep(3)
            rec["composer_visible"] = _visible(page, ".delfin-agent-input textarea")
            _shot(page, args, "opened", rec)
            if args.scenario == "composer":
                _send_prompt(page, args, rec, model_name, t0)
            else:
                _read_running(page, args, rec, t0)
            rec["page_errors"] = errors[:10]
        finally:
            browser.close()


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
        env = server_env(server_python=args.server_python, tree=args.tree,
                         home=home, scratch=scratch, runtime=cleanup.runtime,
                         workspace=ws, model_port=model_port)
        if args.scenario == "service-session":
            _start_service(args, env, work, cleanup, rec, fake_model.MODEL)
        port = free_port()
        dash = subprocess.Popen(
            [str(args.server_python), "-m", "delfin.cli_voila", "--port", str(port),
             "--ip", "127.0.0.1", "--token", TOKEN, "--no-browser"],
            cwd=str(args.tree), env=env,
            stdout=open(args.out / "logs" / "voila.log", "wb"),
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
        rec["timings_s"]["server_ready"] = round(time.time() - t_srv, 1)
        _drive_page(args, url, rec, fake_model.MODEL)
    except BaseException as exc:  # noqa: BLE001 -- recorded, then cleaned up
        rec["error"] = f"{type(exc).__name__}: {exc}"[:500]
        rec["ok"] = False
    finally:
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
