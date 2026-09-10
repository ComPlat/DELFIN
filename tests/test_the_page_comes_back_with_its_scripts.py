"""A resumed page has its widgets and none of the code around them.

The widget objects live in the kernel and come back by themselves. The
JavaScript that turns a box into a viewer does not: the bundled 3Dmol,
every tab's startup script and the editor's two lazy bootstraps were
sent once through ``ctx.run_js``, whose Output widget is cleared by the
next script -- so a session resumed into a new window rendered the
Submit tab with an empty preview, and worse, whatever script had been
sent LAST ran again. After a branch switch that script is
``window.location.reload()``.

What is pinned here:

  kept scripts are kept        keep_js remembers; plain run_js does not
  the page runs them first     resume() emits them as cell output BEFORE
                               the roots, so they have run by the time
                               an Output widget replays a viewer
  the shared output is cleared  a resume does not replay the last one-off
  hooks run after, and alone   one failing hook does not keep the others
                               -- or the page -- off the screen
  the dashboard is wired        the bundle is kept and the getter is
                               registered; the editor's bootstraps are
                               kept
"""

from __future__ import annotations

import inspect
import re

import pytest

from delfin.dashboard import session as S


@pytest.fixture(autouse=True)
def _clean():
    S._reset_for_tests()
    yield
    S._reset_for_tests()


def _ctx():
    from delfin.dashboard.context import DashboardContext

    return DashboardContext()


# ---------------------------------------------------------------------------
# What ctx keeps
# ---------------------------------------------------------------------------

def test_a_kept_script_is_kept_and_a_plain_one_is_not():
    ctx = _ctx()
    ctx.keep_js("window.__lib = 1;")
    ctx.run_js("copyToClipboard('x');")
    ctx.keep_js("window.__boot = 1;")
    assert ctx.page_scripts == ["window.__lib = 1;", "window.__boot = 1;"]


def test_the_bootstrap_joins_the_kept_scripts_in_order_each_guarded():
    """One script throwing must not skip the ones after it."""
    ctx = _ctx()
    ctx.keep_js("a();")
    ctx.keep_js("b();")
    out = ctx.resume_bootstrap_js()
    assert out.index("a();") < out.index("b();")
    assert out.count("try {") == 2 and out.count("catch (e)") == 2


def test_the_bootstrap_clears_the_shared_output():
    """The last one-off must not run again in the new window.

    Outside a kernel an Output widget captures nothing, so the one-off
    is put into it the way the frontend would hold it.
    """
    ctx = _ctx()
    ctx.js_output.outputs = ({
        "output_type": "display_data",
        "data": {"application/javascript": "window.location.reload();"},
        "metadata": {},
    },)
    ctx.resume_bootstrap_js()
    assert ctx.js_output.outputs == ()


def test_ctx_on_resume_reaches_the_session_registry():
    ctx = _ctx()
    calls = []
    ctx.on_resume(lambda: calls.append("hook"))
    assert S.resume_hooks() and S.resume_hooks()[0]() is None
    assert calls == ["hook"]


# ---------------------------------------------------------------------------
# What resume() does, in what order
# ---------------------------------------------------------------------------

def _resume_with_capture(monkeypatch):
    """Run resume() and record what it displayed, in order."""
    shown = []

    def _display(obj):
        shown.append(obj)

    import IPython.display as ipd

    monkeypatch.setattr(ipd, "display", _display)
    return shown


def test_scripts_come_before_the_roots_and_hooks_after(monkeypatch):
    import ipywidgets as widgets
    from IPython.display import Javascript

    shown = _resume_with_capture(monkeypatch)
    order = []
    root = widgets.VBox()
    S.register_root(root)
    S.register_bootstrap(lambda: (order.append("scripts"), "window.__x=1;")[1])
    S.on_resume(lambda: order.append("hook"))

    assert S.resume() is True

    assert isinstance(shown[0], Javascript) and "window.__x=1;" in shown[0].data
    assert shown[1] is root
    assert order == ["scripts", "hook"]


def test_scripts_go_ahead_of_the_root_they_were_registered_against(monkeypatch):
    """At startup the bundle sits at the top of the body, under a header
    that is already there. A resume keeps that order."""
    import ipywidgets as widgets
    from IPython.display import Javascript

    shown = _resume_with_capture(monkeypatch)
    header, body = widgets.VBox(), widgets.VBox()
    S.register_root(header, body)
    S.register_bootstrap(lambda: "boot();", before=body)
    assert S.resume() is True
    assert shown[0] is header
    assert isinstance(shown[1], Javascript) and "boot();" in shown[1].data
    assert shown[2] is body


def test_no_scripts_means_no_script_output(monkeypatch):
    import ipywidgets as widgets

    shown = _resume_with_capture(monkeypatch)
    root = widgets.VBox()
    S.register_root(root)
    S.register_bootstrap(lambda: "   ")
    assert S.resume() is True
    assert shown == [root]


def test_a_failing_hook_is_skipped_not_fatal(monkeypatch, capsys):
    import ipywidgets as widgets

    _resume_with_capture(monkeypatch)
    S.register_root(widgets.VBox())
    ran = []

    def _bad():
        raise RuntimeError("the frame is not there")

    S.on_resume(_bad)
    S.on_resume(lambda: ran.append("second"))
    assert S.resume() is True
    assert ran == ["second"]
    assert "skipped" in capsys.readouterr().out


def test_a_failing_bootstrap_still_shows_the_page(monkeypatch, capsys):
    import ipywidgets as widgets

    shown = _resume_with_capture(monkeypatch)
    root = widgets.VBox()
    S.register_root(root)

    def _boom():
        raise OSError("no scripts today")

    S.register_bootstrap(_boom)
    assert S.resume() is True
    assert shown == [root]
    assert "could not be collected" in capsys.readouterr().out


def test_a_hook_registered_twice_runs_once(monkeypatch):
    import ipywidgets as widgets

    _resume_with_capture(monkeypatch)
    S.register_root(widgets.VBox())
    ran = []
    hook = lambda: ran.append(1)                    # noqa: E731
    S.on_resume(hook)
    S.on_resume(hook)
    S.resume()
    assert ran == [1]


def test_reset_forgets_scripts_and_hooks():
    S.register_bootstrap(lambda: "x")
    S.on_resume(lambda: None)
    S._reset_for_tests()
    assert S.resume_hooks() == []
    assert S._bootstrap is None


# ---------------------------------------------------------------------------
# Wired into the dashboard, not just available
# ---------------------------------------------------------------------------

def test_the_dashboard_keeps_its_bundle_and_registers_the_getter():
    from delfin import dashboard as d

    src = inspect.getsource(d.create_dashboard)
    assert "ctx.keep_js(_calc_init)" in src, (
        "the 3Dmol + startup bundle must be a kept script"
    )
    assert "_session.register_bootstrap(ctx.resume_bootstrap_js, before=_body_root)" in src
    # After the roots are registered and before they are shown.
    assert src.index("_session.register_root(") < src.index(
        "_session.register_bootstrap(") < src.index("display(_header_root)")


def test_the_editor_keeps_its_two_bootstraps():
    from delfin.dashboard import structure_editor

    src = inspect.getsource(structure_editor)
    assert "ctx.keep_js(submit_manip_bootstrap_js())" in src
    assert "ctx.keep_js(molecule_ff_bootstrap_js())" in src


# ---------------------------------------------------------------------------
# The drawing inside the editor's frame
# ---------------------------------------------------------------------------
#
# The frame comes back as HTML and the editor inside it starts empty. The
# panel keeps the last drawing it put in or read out, and a resume loads
# it again once the frame is there.

_MOL = """
  Ketcher  1 1

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""


@pytest.fixture
def panel(monkeypatch, tmp_path):
    from delfin.dashboard import ketcher as k
    from delfin.dashboard import ketcher_panel
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(k, "app_url", lambda: "/static/ketcher/index.html")
    monkeypatch.setattr(k, "installed_version", lambda: "3.17.0")
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    sent = []
    ctx.run_js = sent.append
    built = ketcher_panel.build(ctx, scope="resume-test", folder=tmp_path / "drawings")
    built.sent_js = sent
    return built


def test_the_panel_registers_a_resume_hook(panel):
    assert S.resume_hooks(), "the panel registered nothing for a resume"


def test_nothing_drawn_means_nothing_sent_on_resume(panel):
    panel.sent_js.clear()
    for hook in S.resume_hooks():
        hook()
    assert panel.sent_js == []


def test_a_drawing_put_in_is_given_back_on_resume(panel):
    assert panel.open_text(_MOL, "water") is True
    panel.sent_js.clear()
    for hook in S.resume_hooks():
        hook()
    assert len(panel.sent_js) == 1
    assert "setMolecule" in panel.sent_js[0]
    assert "M  END" in panel.sent_js[0]


def test_a_drawing_read_out_is_the_one_given_back(panel, monkeypatch):
    """What came out of the editor last is newer than what went in."""
    from delfin.dashboard import ketcher as k

    monkeypatch.setattr(
        k, "smiles_from_drawing",
        lambda payload: {"ok": True, "smiles": "O", "status": "read", "reaction": False},
    )
    panel.open_text(_MOL, "first")
    newer = _MOL.replace("O   0", "N   0")
    panel.sync.value = f"7\nsmiles\n{newer}"
    panel.sent_js.clear()
    for hook in S.resume_hooks():
        hook()
    assert len(panel.sent_js) == 1
    assert "N   0" in panel.sent_js[0] and "O   0" not in panel.sent_js[0]


def test_keep_js_goes_through_run_js_so_a_recorder_sees_it():
    """Tests and tabs replace run_js with a one-argument recorder."""
    ctx = _ctx()
    seen = []
    ctx.run_js = seen.append
    ctx.keep_js("lib();")
    assert seen == ["lib();"]
    assert ctx.page_scripts == ["lib();"]
