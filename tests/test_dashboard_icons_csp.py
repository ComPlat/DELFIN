"""The dashboard page shows Font Awesome icons without asking a CDN for them.

Voila's ``lab`` template links ``fontawesome-free@^5`` (``all.min.css`` and
``v4-shims.min.css``) from cdn.jsdelivr.net: a version range without an SRI
hash. The page CSP refuses both, and the page logged two CSP violations on
every load. The same Font Awesome is bundled in Voila's ``voila-style.js``
(fonts served from the dashboard's own origin), so the dashboard's page
template drops the two links instead of the CSP letting them in.

Input: the shipped template (``delfin/voila_assets/jupyter_data/nbconvert/
templates/delfin``), ``cli_voila._page_template_args``, ``web_security.
PAGE_CSP`` and, when installed, Voila's own template and bundle. Properties
asserted:

* the template extends Voila's lab template and removes exactly the lines
  that link the Font Awesome CDN, keeping every other line of the block
  (with and without autoescaping);
* the installed Voila links Font Awesome only on such lines, and bundles
  Font Awesome 5 with its v4 shims in ``voila-style.js`` (the icons'
  source once the links are gone);
* Voila resolves the template through ``JUPYTER_PATH`` and keeps its own
  data directory (themes);
* the launcher selects the template and keeps an existing ``JUPYTER_PATH``;
* the page CSP was not widened: no CDN host or Font Awesome URL in
  ``style-src``/``font-src``, no new ``'unsafe-*'``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

jinja2 = pytest.importorskip("jinja2")

from delfin import cli_voila  # noqa: E402
from delfin.dashboard import pinned_cdn, web_security  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO / "delfin" / "voila_assets" / "jupyter_data"
TEMPLATE_DIR = DATA_ROOT / "nbconvert" / "templates" / cli_voila.PAGE_TEMPLATE
FA_CDN = "cdn.jsdelivr.net/npm/@fortawesome/"

FA_ALL = ('<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/'
          '@fortawesome/fontawesome-free@^5/css/all.min.css" type="text/css" />')
FA_SHIMS = ('<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/'
            '@fortawesome/fontawesome-free@^5/css/v4-shims.min.css" type="text/css" />')

STUB_LAB = (
    "{% block notebook_css %}"
    "<script src=\"static/voila-style.js\"></script>\n"
    "  <style>\n    .jp-mod-noOutputs { display: none; }\n  </style>\n"
    f"  {FA_ALL}\n"
    f"  {FA_SHIMS}\n"
    "<link rel=\"stylesheet\" href=\"static/index.css\">"
    "{% endblock notebook_css %}"
)


def _voila_share() -> Path | None:
    try:
        from voila.paths import collect_template_paths
    except Exception:
        return None
    for path in collect_template_paths(["voila"], "lab"):
        candidate = Path(path) / "index.html.j2"
        if candidate.is_file() and Path(path).name == "lab":
            return Path(path).parents[2]
    return None


def _render_with_stub(autoescape: bool) -> str:
    loader = jinja2.ChoiceLoader([
        jinja2.FileSystemLoader(str(TEMPLATE_DIR)),
        jinja2.DictLoader({"voila/templates/lab/index.html.j2": STUB_LAB}),
    ])
    env = jinja2.Environment(loader=loader, autoescape=autoescape)
    return env.get_template("index.html.j2").render()


def test_the_template_is_shipped_and_based_on_voilas_lab_template():
    conf = json.loads((TEMPLATE_DIR / "conf.json").read_text(encoding="utf-8"))
    assert conf["base_template"] == "lab"
    text = (TEMPLATE_DIR / "index.html.j2").read_text(encoding="utf-8")
    assert "extends 'voila/templates/lab/index.html.j2'" in text
    # Under nbconvert/, never voila/: Voila takes its data directory (the
    # themes) from the first voila/ directory on the Jupyter path.
    assert not (DATA_ROOT / "voila").exists()


@pytest.mark.parametrize("autoescape", [False, True])
def test_the_template_drops_the_cdn_links_and_keeps_every_other_line(autoescape):
    out = _render_with_stub(autoescape)
    assert FA_CDN not in out
    assert "fontawesome" not in out.lower()
    assert '<script src="static/voila-style.js"></script>' in out
    assert ".jp-mod-noOutputs { display: none; }" in out
    assert '<link rel="stylesheet" href="static/index.css">' in out
    assert "&lt;" not in out  # nothing of the parent block got escaped


def test_voila_links_font_awesome_only_on_lines_the_template_drops():
    share = _voila_share()
    if share is None:
        pytest.skip("voila is not installed")
    lab = (share / "voila" / "templates" / "lab" / "index.html.j2").read_text(encoding="utf-8")
    lines = [ln.strip() for ln in lab.splitlines() if "fontawesome" in ln.lower()]
    assert lines, "Voila's lab template no longer links Font Awesome; the template may go"
    for line in lines:
        assert FA_CDN in line and line.startswith("<link ") and line.endswith("/>"), line


def test_voila_bundles_font_awesome_5_with_the_v4_shims():
    share = _voila_share()
    if share is None:
        pytest.skip("voila is not installed")
    bundle = (share / "voila" / "templates" / "base" / "static" / "voila-style.js").read_text(
        encoding="utf-8", errors="replace")
    assert '@font-face{font-family:"Font Awesome 5 Free"' in bundle
    assert ".fa,.fas{font-weight:900}" in bundle
    # ipywidgets writes v4 names (fa fa-refresh, fa-trash-o): the shims
    assert ".fa.fa-trash-o" in bundle


def test_voila_resolves_the_template_and_keeps_its_own_data_directory(monkeypatch):
    try:
        from voila.paths import collect_template_paths
        from voila.utils import get_data_dir
    except Exception:
        pytest.skip("voila is not installed")
    before = get_data_dir()
    env: dict = {"JUPYTER_PATH": os.environ.get("JUPYTER_PATH", "")}
    assert cli_voila._page_template_args(env)
    monkeypatch.setenv("JUPYTER_PATH", env["JUPYTER_PATH"])
    paths = collect_template_paths(["voila", "nbconvert"], cli_voila.PAGE_TEMPLATE)
    existing = [p for p in paths if Path(p).is_dir()]
    assert Path(existing[0]).resolve() == TEMPLATE_DIR.resolve()
    assert any(Path(p).name == "lab" for p in paths)
    loader = jinja2.FileSystemLoader(paths)
    env_j = jinja2.Environment(loader=loader)
    assert Path(env_j.get_template("index.html.j2").filename).resolve() == (
        TEMPLATE_DIR / "index.html.j2").resolve()
    parent = env_j.get_template("voila/templates/lab/index.html.j2").filename
    assert "fontawesome" in Path(parent).read_text(encoding="utf-8").lower()
    assert get_data_dir() == before


def test_the_launcher_selects_the_template_and_keeps_jupyter_path(tmp_path):
    env: dict = {}
    assert cli_voila._page_template_args(env) == [
        f"--VoilaConfiguration.template={cli_voila.PAGE_TEMPLATE}"]
    assert Path(env["JUPYTER_PATH"]).resolve() == DATA_ROOT.resolve()

    env = {"JUPYTER_PATH": str(tmp_path)}
    cli_voila._page_template_args(env)
    cli_voila._page_template_args(env)  # twice: no duplicate entry
    parts = env["JUPYTER_PATH"].split(os.pathsep)
    assert [Path(p).resolve() for p in parts] == [DATA_ROOT.resolve(), tmp_path.resolve()]


def _directive(name: str) -> list[str]:
    for part in web_security.PAGE_CSP.split(";"):
        words = part.split()
        if words and words[0] == name:
            return words[1:]
    raise AssertionError(name)


def test_the_page_csp_was_not_widened_for_the_icons():
    assert _directive("style-src") == [
        "'self'", "'unsafe-inline'", *pinned_cdn.urls("style")]
    assert _directive("font-src") == ["'self'", "data:"]
    assert "fontawesome" not in web_security.PAGE_CSP.lower()
