"""After login the access token leaves the address bar.

Found in the bug hunt of 11.11. (job 8543241): after loading
``.../delfin_dashboard.ipynb?token=...`` the token stayed in
``location.search`` -- in the browser history, in screenshots, in stack
traces. Jupyter Notebook and JupyterLab drop it with
``history.replaceState`` once it is consumed (the first request answers it
with the login cookie). The page's start-up script does the same now.
"""

from __future__ import annotations

import json
import shutil
import subprocess

import pytest

from delfin.dashboard import login_session

_NODE = shutil.which("node")


def test_the_startup_script_strips_the_token():
    js = login_session.keepalive_js()
    assert login_session.strip_token_js() in js
    assert "replaceState" in js
    assert "searchParams.delete('token')" in js


def _run(href):
    script = (
        "var calls = [];\n"
        f"var window = {{location: {{href: {json.dumps(href)}}},\n"
        "  history: {state: null, replaceState: function(s, t, u){ calls.push(u); }}};\n"
        "var document = {title: 'DELFIN'};\n"
        + login_session.strip_token_js()
        + "console.log(JSON.stringify(calls));\n"
    )
    proc = subprocess.run([_NODE, "-e", script], capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


@pytest.mark.skipif(_NODE is None, reason="node is not installed")
def test_only_the_token_goes_and_everything_else_stays():
    calls = _run("http://127.0.0.1:8899/voila/render/delfin_dashboard.ipynb"
                 "?token=s3cret&tab=jobs#top")
    assert calls == ["/voila/render/delfin_dashboard.ipynb?tab=jobs#top"]
    assert _run("http://127.0.0.1:8899/voila/render/x.ipynb?token=s3cret") == [
        "/voila/render/x.ipynb"]


@pytest.mark.skipif(_NODE is None, reason="node is not installed")
def test_an_address_without_a_token_is_left_alone():
    assert _run("http://127.0.0.1:8899/voila/render/x.ipynb?tab=jobs") == []
