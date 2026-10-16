"""MathJax's loader shim never stopped trying, and it cost the whole machine.

The page asks for MathJax at its default entry point, ``latest.min.js``. That
file is not MathJax: it is a shim that works out which version is newest and
then appends a script tag for it. On a host that cannot reach the CDN -- a
compute cluster, most of the time -- it never stops. Measured on the running
dashboard with a CPU profile, the hot frames were ``e.loadLatest`` and the
``appendChild`` under it.

What it cost, and what removing it gave back, measured on the same page:

    idle 5 s, main thread busy    4.895 s  ->  0.064 s
    idle 5 s, script              4.098 s  ->  0.036 s
    document                      7626 elements, 15 script tags, either way

That is 98% of the main thread held by a page that was doing nothing, which is
what a dashboard that "hinkt hinterher" feels like.

Naming the versioned bundle removes the shim: where there is a network it is
the same MathJax, where there is not it fails once and is quiet.
"""

import json
import re

from delfin import cli_voila


def _voila_command(monkeypatch):
    """The argument list the launcher would run, without running it."""
    captured = {}

    def fake_execvp(*args, **kwargs):
        captured['args'] = args
        raise SystemExit(0)

    monkeypatch.setattr(cli_voila.os, 'execvp', fake_execvp, raising=False)
    return captured


def _dewrapped(text):
    """The source with adjacent string literals joined, as Python sees them."""
    return re.sub(r"'\s*\n\s*'", '', text)


def test_the_launcher_names_a_bundle_and_not_the_loader():
    source = _dewrapped(open(cli_voila.__file__, encoding='utf-8').read())

    assert 'mathjax_url' in source, 'the launcher must say which MathJax it wants'
    assert 'latest.min.js' not in source, 'that is the shim that spins'
    # Nothing at all, not a working MathJax: window.MathJax is undefined in
    # this dashboard today and no output is typeset, so a real bundle would
    # start typesetting where there is a network -- a change in what is shown.
    # The nothing is an empty script on the dashboard's own origin, so the
    # page CSP needs no data: in script-src.
    assert '"mathjax_url": NO_TYPESETTING_ROUTE' in source
    assert 'data:text/javascript' not in source


def test_it_is_set_where_this_deployment_reads_it():
    """The dashboard runs Voila as a jupyter-server extension, and in that
    mode the Voila app's own ``mathjax_url`` trait is never consulted -- the
    handler reads the tornado settings. Setting the trait looked right and did
    nothing: the page kept the shim URL."""
    source = open(cli_voila.__file__, encoding='utf-8').read()

    assert '--ServerApp.tornado_settings=' in source
    assert '--Voila.mathjax_url=' not in source


def test_the_setting_is_valid_json_tornado_can_read():
    from delfin.dashboard import web_security

    source = open(cli_voila.__file__, encoding='utf-8').read()
    assert ('"--ServerApp.tornado_settings="\n'
            '        + json.dumps({"mathjax_url": NO_TYPESETTING_ROUTE})'
            in source)
    payload = json.loads(json.dumps(
        {"mathjax_url": web_security.NO_TYPESETTING_ROUTE}))
    assert set(payload) == {'mathjax_url'}
    # Voila appends ?config=... to whatever it is given; the route is a
    # path, so that is a query string the route ignores (the server-level
    # check: test_the_dashboard_serves_files_only_to_its_owner).
    assert payload['mathjax_url'].startswith('/')
    assert payload['mathjax_url'].endswith('.js')
