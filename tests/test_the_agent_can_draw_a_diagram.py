"""A fenced ``mermaid`` block is a diagram in the dashboard chat.

Input: an answer containing a ```mermaid fence. Output: a container holding
the diagram source twice -- as text inside a ``<pre>``, which is what is on
screen until anything else happens, and in a data attribute the renderer
reads.

Semantics. The readable source is the resting state, not an error path.
Three outcomes, by construction rather than by detection:

    a vendored static/mermaid.min.js   drawn, no network
    outbound network                   drawn, fetched once per page, pinned
    neither                            the source, which was already there

Nothing has to notice a failure, because nothing had to succeed first. The
terminal REPL renders text only -- no sixel, no iTerm2, no kitty protocol --
so there the block is the source, which is the same fallback.

The library is not committed: 3.5 MB against a repository whose largest
tracked file is 797 KB. The pinned CDN is the default and a site without
outbound network drops the file in once.
"""

from __future__ import annotations

import re

import pytest

from delfin.dashboard import tab_agent as T


def _html(md: str) -> str:
    return T._md_to_html(md)


# -- what a fence becomes --------------------------------------------------

def test_a_mermaid_fence_becomes_a_diagram_container():
    out = _html("before\n\n```mermaid\ngraph TD\n  A-->B\n```\n\nafter")
    assert 'class="delfin-mermaid"' in out
    assert "data-mmd=" in out
    assert "before" in out and "after" in out


def test_the_source_is_on_screen_before_anything_renders():
    """The fallback is the starting state, so a machine with no library and
    no network shows the same text a terminal would."""
    out = _html("```mermaid\ngraph TD\n  A-->B\n```")
    assert "delfin-mermaid-src" in out
    body = re.search(r"<pre class=\"delfin-mermaid-src\"><code>(.*?)</code>",
                     out, re.S)
    assert body and "graph TD" in body.group(1)
    assert "A--&gt;B" in body.group(1), "escaped, and still readable"


def test_another_language_is_still_a_code_block():
    out = _html("```python\nx = 1\n```")
    assert "delfin-mermaid" not in out
    assert "delfin-code-wrap" in out


def test_an_empty_mermaid_fence_does_not_crash():
    out = _html("```mermaid\n```")
    assert isinstance(out, str)


# -- the source comes from a model and lands in a page ---------------------

@pytest.mark.parametrize("payload", [
    '</div><script>alert(1)</script>',
    '"><img src=x onerror=alert(1)>',
    "graph TD\n  A[<script>alert(1)</script>]-->B",
])
def test_nothing_in_the_source_becomes_markup(payload):
    """The property is that no ELEMENT is created, not that a string is
    absent: `onerror=` occurs harmlessly as escaped text, and asserting on
    the string rather than the structure is how a test passes while the
    page is unsafe -- or fails while it is fine, which is what happened
    when this was first written.
    """
    out = _html("```mermaid\n" + payload + "\n```")
    # Every tag in the output is one this function emitted.
    import re as _re
    tags = {t.lower() for t in _re.findall(r"<\s*/?\s*([a-zA-Z][\w-]*)", out)}
    assert tags <= {"div", "pre", "code"}, f"unexpected elements: {tags}"
    # And the payload survived as text, so the reader still sees it.
    assert "&lt;" in out or "&quot;" in out or "&gt;" in out


def test_the_data_attribute_cannot_be_closed_early():
    out = _html('```mermaid\ngraph TD\n  A["a quote \\" here"]-->B\n```')
    head = out.split("<pre", 1)[0]
    assert head.count('data-mmd="') == 1
    # Everything between the opening quote and the container's close is one
    # attribute value: a bare quote would end it and start an attribute.
    assert '"' not in head.split('data-mmd="', 1)[1].rstrip('">')


# -- how the library is reached -------------------------------------------

def test_the_cdn_script_is_pinned_to_a_version_and_to_its_bytes():
    """An unpinned CDN script is a supply-chain surface. With integrity the
    browser refuses a changed file instead of running it."""
    assert re.fullmatch(r"\d+\.\d+\.\d+", T._MERMAID_VERSION)
    assert T._MERMAID_VERSION in T._MERMAID_URL
    assert T._MERMAID_SRI.startswith("sha384-")
    assert len(T._MERMAID_SRI) > 40


def test_a_missing_vendored_build_is_not_an_error():
    assert T.vendored_mermaid_js() == "" or isinstance(
        T.vendored_mermaid_js(), str)


def test_a_vendored_build_is_used_when_a_site_installs_one(monkeypatch,
                                                           tmp_path):
    monkeypatch.setattr(T, "_VENDORED_MERMAID_CACHE", None)
    fake = tmp_path / "mermaid.min.js"
    fake.write_text("/* a build */", encoding="utf-8")

    class _Files:
        def joinpath(self, rel):
            assert rel.endswith("mermaid.min.js")
            return fake

    monkeypatch.setattr("importlib.resources.files", lambda pkg: _Files())
    assert "a build" in T.vendored_mermaid_js()


def test_the_library_is_not_committed():
    """Measured: 3.5 MB, against a 797 KB largest tracked file. The pinned
    CDN is the default; a site without network installs it deliberately."""
    import pathlib

    static = (pathlib.Path(T.__file__).resolve().parent / "static")
    shipped = static / "mermaid.min.js"
    assert not shipped.exists(), (
        "the build is meant to be a site's choice, not a clone's weight")


# -- the model has to be told ---------------------------------------------

def test_the_built_prompt_says_a_diagram_is_possible(monkeypatch):
    """A capability the prompt does not mention is one the model does not
    use. Checked in the COMPOSED prompt, not in the file: a role file is
    not automatically part of what a session receives."""
    from delfin import user_settings
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(user_settings, "load_settings",
                        lambda *a, **k: {"agent": {"slim_prompt": True}})
    built = PromptLoader().build_system_prompt(
        role_id="dashboard_agent", mode_id="dashboard",
        route=["dashboard_agent"], task_text="explain the architecture",
        session_key="mermaid-1")
    text = built if isinstance(built, str) else str(built)
    assert "mermaid" in text.lower(), "the diagram rule never reaches a model"


# -- the renderer is wired in ---------------------------------------------

def test_the_renderer_is_registered_with_the_page():
    import inspect

    src = inspect.getsource(T)
    assert "_mermaid_init_js" in src
    assert "ctx.add_init_js(_mermaid_init_js)" in src


def test_the_renderer_is_pinned_to_the_strict_security_level():
    """The source is model-written and rendered into the page; strict
    sanitises labels and htmlLabels:false keeps them text to begin with."""
    import inspect

    src = inspect.getsource(T)
    assert "securityLevel: 'strict'" in src
    assert "htmlLabels: false" in src
