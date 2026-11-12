"""After login the access token leaves the address bar.

Found in the bug hunt of 11.11. (job 8543241): after loading
``.../delfin_dashboard.ipynb?token=...`` the token stayed in
``location.search`` -- in the browser history, in screenshots, in stack
traces. Jupyter Notebook and JupyterLab drop it with
``history.replaceState`` once it is consumed (the page request answers it
with the login cookie). The dashboard's page template does the same now, at
the end of the page.

The first attempt put the snippet into the login start-up script
(``login_session.keepalive_js``); the browser check (job 8544581) showed
that script does not run on a fresh page load, so it lives in the page
template, which every load runs.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from delfin import cli_voila

jinja2 = pytest.importorskip("jinja2")

_NODE = shutil.which("node")
TEMPLATE_DIR = (Path(cli_voila.__file__).resolve().parent / "voila_assets" / "jupyter_data"
                / "nbconvert" / "templates" / cli_voila.PAGE_TEMPLATE)

# Voila's lab template, reduced to the blocks the dashboard template uses.
STUB_LAB = (
    "<html>{% block notebook_css %}<link rel=\"stylesheet\" href=\"static/index.css\">"
    "{% endblock notebook_css %}"
    "<body>cells</body>"
    "{%- block footer %}{% block footer_js %}{% endblock footer_js %}</html>"
    "{%- endblock footer -%}"
)


def _rendered() -> str:
    loader = jinja2.ChoiceLoader([
        jinja2.FileSystemLoader(str(TEMPLATE_DIR)),
        jinja2.DictLoader({"voila/templates/lab/index.html.j2": STUB_LAB}),
    ])
    env = jinja2.Environment(loader=loader, autoescape=True)
    return env.get_template("index.html.j2").render()


def _script() -> str:
    match = re.search(r'<script type="text/javascript" id="delfin-strip-token">(.*?)</script>',
                      _rendered(), re.S)
    assert match, "the page carries no token-stripping script"
    return match.group(1)


def test_the_page_ends_with_the_token_stripping_script():
    out = _rendered()
    script = _script()
    assert "replaceState" in script and "searchParams.delete('token')" in script
    assert out.index("delfin-strip-token") > out.index("cells")  # at the end of the page
    assert "&#39;" not in script  # not escaped by autoescape


def _run(href):
    script = (
        "var calls = [];\n"
        f"var window = {{location: {{href: {json.dumps(href)}}},\n"
        "  history: {state: null, replaceState: function(s, t, u){ calls.push(u); }}};\n"
        "var document = {title: 'DELFIN'};\n"
        + _script()
        + "\nconsole.log(JSON.stringify(calls));\n"
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
