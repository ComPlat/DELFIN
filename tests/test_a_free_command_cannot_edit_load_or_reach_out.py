"""What runs unasked writes nothing, loads no script and reaches nothing.

Three forms that the stack #328-#332 made free, and the form each keeps:

* sed: every prefix GNU sed accepts for ``--in-place`` and ``--file``
  (``--i``, ``--fi``) and ``-i``/``-f`` in a bundle edited a file or ran a
  script file's ``e`` command unasked; ``sed -n -f`` already did on main.
* An install named ``..``, ``.`` or a path outside the workspace installs
  that directory.
* The headless page check runs on the host, outside the cage and the egress
  proxy. Measured 2026-10-19 with Chromium through Playwright 1.x: before,
  an <img>, fetch(), a WebSocket, a ``localhost`` name and a WebRTC STUN
  packet each reached a loopback listener, and a script outside the page's
  directory put a word of itself into the error returned to the model;
  after, no request and no packet arrived, and errors of the page's own
  scripts were still reported. The browser is not installed on every
  machine, so the checks below drive a stand-in for Playwright and pin
  what the measurement depended on.
"""

from __future__ import annotations

import sys
import types

import pytest

from delfin.agent import artifact_check as AC
from delfin.agent import caged_test_runners as C
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=tmp_path)


@pytest.mark.parametrize("cmd", [
    "sed --i 's/a/b/' .delfin/settings.json",
    "sed --in-pl=.bak 's/a/b/' f.txt",
    "sed --in-place 's/a/b/' f.txt",
    "sed -i 's/a/b/' f.txt",
    "sed -ni 's/a/b/' f.txt",
    "sed -es/x/y/ -i f.txt",
    "sed --fi=x.sed f.txt",
    "sed --fil x.sed f.txt",
    "sed -nf x.sed f.txt",
    "sed -n -f x.sed f.txt",
    "sed '1e id' f.txt",
    "sed 's/a/b/w out' f.txt",
])
def test_a_sed_that_edits_or_loads_a_script_is_asked(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is False, cmd


@pytest.mark.parametrize("cmd", [
    "sed 's/a/b/' f.txt",
    "sed '1d;$d' f.txt",
    "sed -n '1,20p' f.txt",
    "sed -nE '/fi/p' f.txt",
    "sed -E 's/(a)b/\\1/g' f.txt",
    "sed -e 's/a/b/' -e 's/c/d/' f.txt",
    "sed --expression='s/a/b/' f.txt",
    "sed -s -n '$p' a.txt b.txt",
    "sed --follow-symlinks -n 1p f.txt",
    "sed -n '/<script>/,/<\\/script>/p' game.html",
])
def test_a_sed_that_prints_stays_free(perms, cmd):
    assert perms.matches_bash_auto_allow(cmd) is True, cmd


@pytest.mark.parametrize("cmd", [
    "npm install ..", "npm install .", "npm install ../x",
    "uv pip install ..", "uv pip install ../other", "uv pip install /abs/pkg",
])
def test_an_install_of_a_directory_outside_is_not_a_project_install(cmd):
    assert C.is_project_install(cmd) is False, cmd


@pytest.mark.parametrize("cmd", [
    "npm install three", "npm install three@0.160.1",
    "npm i -D @types/node@^20", "uv pip install numpy",
    "uv pip install -e .", "uv pip install .",
    "uv pip sync requirements.txt", "uv venv .venv --python 3.11",
])
def test_an_install_into_the_project_stays_one(cmd):
    assert C.is_project_install(cmd) is True, cmd


class _Fake:
    """Playwright's sync API as far as _browser_errors uses it. ``requests``
    are the URLs the page asks for, in order; each goes through the route
    and then the page's request listener."""

    def __init__(self, requests, errors):
        self.requests, self.errors = requests, errors
        self.launch_args = None
        self.context_kwargs = None
        self.aborted: list[str] = []
        self.continued: list[str] = []

    def install(self, monkeypatch):
        fake = self

        class Req:
            def __init__(self, url):
                self.url = url

        class Route:
            def __init__(self, url):
                self.request = Req(url)

            def continue_(self):
                fake.continued.append(self.request.url)

            def abort(self):
                fake.aborted.append(self.request.url)

        class Page:
            def __init__(self, route_fn):
                self.route_fn, self.handlers = route_fn, {}

            def on(self, event, fn):
                self.handlers[event] = fn

            def goto(self, url, **kw):
                for u in [url, *fake.requests]:
                    self.route_fn(Route(u))
                    self.handlers["request"](Req(u))
                for e in fake.errors:
                    self.handlers["pageerror"](e)

            def wait_for_timeout(self, ms):
                pass

        class Context:
            def route(self, pattern, fn):
                self.fn = fn

            def new_page(self):
                return Page(self.fn)

        class Browser:
            def new_context(self, **kw):
                fake.context_kwargs = kw
                return Context()

            def close(self):
                pass

        class Chromium:
            def launch(self, headless=True, args=()):
                fake.launch_args = list(args)
                return Browser()

        class PW:
            chromium = Chromium()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

        mod = types.ModuleType("playwright.sync_api")
        mod.sync_playwright = PW
        monkeypatch.setitem(sys.modules, "playwright", types.ModuleType("playwright"))
        monkeypatch.setitem(sys.modules, "playwright.sync_api", mod)
        return self


def test_the_browser_resolves_nothing_and_connects_nowhere(tmp_path, monkeypatch):
    page = tmp_path / "game.html"
    page.write_text("<script>broken()</script>")
    fake = _Fake([], ["broken is not defined"]).install(monkeypatch)
    out = AC._browser_errors(page)
    assert out and "broken is not defined" in out[0]
    args = " ".join(fake.launch_args)
    assert "--host-resolver-rules=MAP * ~NOTFOUND" in args
    assert "--proxy-server=" in args
    assert "<-loopback>" in args        # the loopback goes through it too
    assert "disable_non_proxied_udp" in args
    assert fake.context_kwargs.get("offline") is True
    assert fake.context_kwargs.get("service_workers") == "block"


def test_a_page_that_reaches_out_is_refused_and_not_judged(tmp_path, monkeypatch):
    page = tmp_path / "game.html"
    page.write_text("<script>THREE.Scene()</script>")
    url = "https://cdn.jsdelivr.net/npm/three/build/three.module.js"
    fake = _Fake([url], ["THREE is not defined"]).install(monkeypatch)
    assert AC._browser_errors(page) == []
    assert fake.aborted == [url]


def test_a_file_outside_the_page_directory_is_refused(tmp_path, monkeypatch):
    sub = tmp_path / "site"
    sub.mkdir()
    page = sub / "index.html"
    page.write_text("<script src='../../secret.js'></script>")
    outside = (tmp_path.parent / "secret.js").as_uri()
    fake = _Fake([outside], ["Unexpected identifier 'nope'"]).install(monkeypatch)
    assert AC._browser_errors(page) == []
    assert fake.aborted == [outside]


def test_the_page_s_own_files_are_loaded_and_judged(tmp_path, monkeypatch):
    page = tmp_path / "index.html"
    page.write_text("<script src='lib.js'></script>")
    (tmp_path / "lib.js").write_text("missingThing()")
    own = (tmp_path / "lib.js").as_uri()
    fake = _Fake([own, "data:image/png;base64,AA=="],
                 ["missingThing is not defined"]).install(monkeypatch)
    out = AC._browser_errors(page)
    assert out and "missingThing" in out[0]
    assert fake.aborted == []
    assert own in fake.continued


@pytest.mark.parametrize("cmd,free", [
    (".venv/bin/python -m pip install numpy", True),
    ("proj/.venv/bin/python -m pip install numpy==1.26 scipy>=1", True),
    (".venv/bin/python -m pip install -r requirements.txt -U", True),
    ("a/../../venv/bin/python -m pip install numpy", False),
    (".venv/bin/python -m pip install -r ../r.txt", False),
    ("/abs/.venv/bin/python -m pip install x", False),
    (".venv/bin/python -m pip install --index-url http://x y", False),
])
def test_pip_through_the_venv_python_stays_in_the_workspace(cmd, free):
    from delfin.agent import api_client as A
    assert A._venv_pip_module_install(cmd) is free, cmd


def test_a_long_pip_line_is_judged_in_linear_time():
    """The one-regex form backtracked exponentially (code scanning
    py/redos on #332)."""
    import time
    from delfin.agent import api_client as A
    cmd = ".venv/bin/python -m pip install" + " -r a.txt" * 5000 + " ./x"
    t = time.monotonic()
    assert A._venv_pip_module_install(cmd) is False
    assert time.monotonic() - t < 1.0
