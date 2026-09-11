"""The session mechanism, driven the way a person drives it.

Every part of this was already covered by tests that pass without a
browser, and the browser found three things none of them could:

  the control was not on the page   ``build_status_strip`` existed, had
                                    tests, and ``create_dashboard`` never
                                    displayed it. Neither did it display
                                    the heartbeat field, so no page ever
                                    beat and the watchdog never armed --
                                    the default teardown did not happen
                                    either

  a kernel id is not a kernel       an unkept kernel WAS killed, and the
                                    server put a replacement in its place
                                    under the same id. Reading
                                    ``/api/kernels`` showed a survivor
                                    where there was a resurrection

  ending has to be asked for        exiting on our own looks like a crash,
                                    which is what invites the restarter --
                                    so the SERVER ends it now, on its own
                                    count of windows

So the checks here are about the process behind the id, not the id, and
they run against a real Voila server with a real chromium.

Slow on purpose: two dashboard renders and one grace period. Skipped
wherever there is no browser, no voila, or no free port, which is most
machines.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]

#: Short enough that a test can wait it out, long enough that a slow
#: render is not mistaken for a window that closed.
_GRACE = 10.0

_WATER = """3
water
O  0.000  0.000  0.117
H  0.000  0.757 -0.469
H  0.000 -0.757 -0.469
"""
_TOKEN = "delfin-browser-session-test-token"


def _requirements():
    """Everything that must be present, or the reason it is not."""
    try:
        import playwright.sync_api  # noqa: F401
    except Exception:
        return "playwright is not installed"
    try:
        import voila  # noqa: F401
    except Exception:
        return "voila is not installed"
    return ""


pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(bool(_requirements()), reason=_requirements() or "ok"),
]


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def _kernels(root: str) -> list:
    req = urllib.request.Request(f"{root}/api/kernels?token={_TOKEN}")
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp)


def _pid_for(kid: str):
    """The OS process behind a kernel id, or None.

    The id alone cannot answer whether a kernel survived: a restarted
    kernel keeps it. The connection file in the command line is the link
    back to a process.

    Read from /proc rather than ps on purpose. ps truncates its output to
    the width in COLUMNS, pytest sets COLUMNS, and the connection path
    sits past column 80 -- so the same helper answered correctly in a
    terminal and returned None under the suite.
    """
    needle = f"kernel-{kid}.json"
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/cmdline", "rb") as handle:
                cmd = handle.read().decode("utf-8", "replace")
        except OSError:
            continue
        if needle in cmd:
            return int(entry)
    return None


def _wait_for_pid(kid: str, timeout: float = 20.0):
    """The process may be a moment behind the registry entry."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        pid = _pid_for(kid)
        if pid:
            return pid
        time.sleep(0.5)
    return None


def _alive(pid) -> bool:
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    """A real dashboard server, with its own state directory."""
    port = _free_port()
    records = tmp_path_factory.mktemp("kept_sessions")
    env = dict(os.environ)
    env.update({
        "PYTHONPATH": str(_REPO),
        "JUPYTER_TOKEN": _TOKEN,
        "DELFIN_SESSION_GRACE_SECONDS": str(_GRACE),
        "DELFIN_SESSION_RECORD_DIR": str(records),
    })
    # The launcher generates its own token unless it is given one, so a
    # test that only sets JUPYTER_TOKEN polls a 403 until it gives up.
    log = tmp_path_factory.mktemp("server") / "voila.log"
    with open(log, "wb") as handle:
        proc = subprocess.Popen(
            [sys.executable, "-m", "delfin.cli_voila",
             "--port", str(port), "--ip", "127.0.0.1", "--token", _TOKEN],
            cwd=str(_REPO), env=env, stdout=handle, stderr=subprocess.STDOUT,
        )
    root = f"http://127.0.0.1:{port}"
    last = ""
    deadline = time.time() + 300
    while time.time() < deadline:
        if proc.poll() is not None:
            pytest.skip(f"the server exited while starting: "
                        f"{log.read_text()[-400:]}")
        try:
            _kernels(root)
            break
        except Exception as exc:
            last = f"{type(exc).__name__}: {exc}"
            time.sleep(2)
    else:
        proc.terminate()
        pytest.skip(f"the server did not answer in 300s ({last}); "
                    f"log: {log.read_text()[-400:]}")
    try:
        yield root, records, log
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=30)
        except Exception:
            proc.kill()


def _open_dashboard(page, root: str):
    page.goto(f"{root}/?token={_TOKEN}", wait_until="domcontentloaded",
              timeout=180_000)
    page.wait_for_selector(".delfin-session-strip", timeout=240_000)



_DIAG = """() => ({
  lib: typeof window.$3Dmol,
  chatScroll: window.__delfinChatScroll,
  agentKeys: window.__delfinAgentKeys,
  canvases: document.querySelectorAll('.submit-mol-output canvas').length,
  previewText: (document.querySelector('.submit-mol-output') || {}).innerText,
  scripts: document.querySelectorAll('.jp-OutputArea script, .widget-output script').length,
  outputs: document.querySelectorAll('.jp-OutputArea-output').length,
})"""


def _expect(page, what: str, expression: str, timeout_ms: int):
    """Wait for *expression*; on timeout say what the page looks like."""
    from playwright.sync_api import TimeoutError as _Timeout

    try:
        page.wait_for_function(f"() => {expression}", timeout=timeout_ms)
    except _Timeout:
        diag = page.evaluate(_DIAG)
        errors = [e[:300] for e in getattr(page, "_delfin_errors", [])[:8]]
        raise AssertionError(
            f"{what} did not come back with the page: {diag}; errors: {errors}")


def test_the_control_is_on_the_page(server):
    """The strip reaches the rendered page.

    It was built and not displayed once; nothing without a browser could
    see that.
    """
    from playwright.sync_api import sync_playwright

    root, _, _ = server
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        try:
            _open_dashboard(page, root)
            assert page.locator("button", has_text="Keep session").count() >= 1
            # The baseline for the resume case: a fresh page has the
            # scripts an Output widget carries, not only the bundle.
            _expect(page, "the agent tab's Output-carried script",
                    "window.__delfinChatScroll === true", 60_000)
        finally:
            browser.close()


def test_a_session_nobody_kept_ends_and_stays_ended(server):
    """The default, which is the behaviour the opt-in is an exception to."""
    from playwright.sync_api import sync_playwright

    root, _, _ = server
    before = {k["id"] for k in _kernels(root)}
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        try:
            _open_dashboard(page, root)
            new = {k["id"] for k in _kernels(root)} - before
            assert len(new) == 1, f"expected one new kernel, saw {new}"
            kid = new.pop()
            pid = _wait_for_pid(kid)
            assert pid, f"no process found for kernel {kid}"
            # Nothing to arm: the server counts this window itself.
            page.wait_for_timeout(2_000)
        finally:
            browser.close()

    deadline = time.time() + _GRACE * 6 + 60
    while time.time() < deadline and _alive(pid):
        time.sleep(1)

    assert not _alive(pid), "the kernel outlived the window that owned it"
    # The replacement is the part a kernel id cannot show: the restarter
    # reuses it, so a live id here means the teardown built a leak.
    assert _pid_for(kid) is None, "a replacement kernel took the same id"
    assert kid not in {k["id"] for k in _kernels(root)}


def test_a_kept_session_survives_the_window_and_comes_back(server):
    """The opt-in: the window goes, the kernel and its state stay, and
    the address in the strip renders from that same kernel."""
    from playwright.sync_api import sync_playwright

    root, records, log = server
    before = {k["id"] for k in _kernels(root)}
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        try:
            _open_dashboard(page, root)
            new = {k["id"] for k in _kernels(root)} - before
            assert len(new) == 1
            kid = new.pop()
            pid = _wait_for_pid(kid)
            assert pid

            # Something on the page that only script can draw: a
            # structure in the Submit preview. The widgets come back by
            # themselves; the viewer is JavaScript and has to be given
            # its library and its bootstrap again.
            box = page.locator("textarea[placeholder^='Paste XYZ']").first
            box.fill(_WATER)
            page.wait_for_selector(".submit-mol-output canvas", timeout=90_000)

            page.locator("button", has_text="Keep session").first.click()
            page.wait_for_selector(".delfin-session-dot.on", timeout=60_000)
            note = page.locator(".delfin-session-dot.on").last.get_attribute("title") or ""
            resume = page.locator(".delfin-session-link").last.get_attribute("href") or ""
        finally:
            browser.close()

        assert "Kept as" in note, f"the dot does not say it is kept: {note!r}"
        assert resume.startswith("http"), f"the dot links nowhere: {resume!r}"

        # The return address is behind the same token as the dashboard:
        # without it nobody gets in, with it only the person who started
        # the server. A name alone opens nothing.
        import re as _re
        naked = _re.sub(r"([?&])token=[^&]*&?", r"\1", resume).rstrip("?&")
        assert "token=" not in naked

        class _NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, *a, **k):
                return None
        opener = urllib.request.build_opener(_NoRedirect)
        try:
            with opener.open(naked, timeout=30) as resp:
                status, body = resp.status, resp.read(4000).decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            status, body = exc.code, ""
        assert status != 200 or "delfin-session-strip" not in body, (
            f"the return address answered without a token: HTTP {status}")

        # A record is what a later request follows back to this kernel.
        assert list(Path(records).glob("*.json")), "the session announced nothing"

        # And the terminal that runs the server is told. A print from
        # the kernel does not get there -- ipykernel forwards it to the
        # frontend -- so this reads the server's own output.
        deadline = time.time() + 15
        while time.time() < deadline and "Return:" not in log.read_text(errors="replace"):
            time.sleep(0.5)
        told = log.read_text(errors="replace")
        assert "Return:" in told and resume.split("?")[0] in told, (
            f"the server terminal was not told the address; log tail: {told[-300:]!r}"
        )

        # Well past the grace, in the same units the watchdog uses.
        time.sleep(_GRACE * 3 + 10)
        assert _alive(pid), "a kept session was torn down anyway"
        assert _pid_for(kid) == pid, "the kept kernel was replaced"

        # And the address comes back INTO that kernel: no new one appears.
        listed_before_resume = {k["id"] for k in _kernels(root)}
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        errors: list = []
        page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
        page.on("console", lambda m: errors.append(f"console.{m.type}: {m.text}")
                if m.type == "error" else None)
        page._delfin_errors = errors
        try:
            page.goto(resume, wait_until="domcontentloaded", timeout=180_000)
            page.wait_for_selector(".delfin-session-strip", timeout=240_000)
            # The control comes back armed: it is the same widget object.
            assert "Kept as" in (page.locator(
                ".delfin-session-dot.on").last.get_attribute("title") or "")
            # The scripts came back with the page: the library, the
            # startup scripts, and a viewer drawn from them. Before this
            # was built, a resumed page showed the box and no molecule.
            _expect(page, "the 3Dmol library",
                    "typeof window.$3Dmol !== 'undefined'", 60_000)
            _expect(page, "the agent tab's startup script",
                    "window.__delfinChatScroll === true", 60_000)
            _expect(page, "a viewer drawn in the Submit preview",
                    "!!document.querySelector('.submit-mol-output canvas')", 90_000)
            assert page.locator("textarea[placeholder^='Paste XYZ']").first.input_value().strip() == _WATER.strip()
        finally:
            browser.close()

        assert {k["id"] for k in _kernels(root)} == listed_before_resume, (
            "the resume address started a new kernel instead of returning "
            "to the kept one"
        )
        assert _pid_for(kid) == pid
