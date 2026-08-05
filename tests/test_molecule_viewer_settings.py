import importlib.util
import json
from pathlib import Path


_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "molecule_viewer.py"
)
_SPEC = importlib.util.spec_from_file_location("delfin_molecule_viewer", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load molecule_viewer module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


_SETTINGS_UI_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "tab_settings.py"
).read_text(encoding="utf-8")

_ORCA_VIEWER_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "tab_orca_builder.py"
).read_text(encoding="utf-8")

_CALC_VIEWER_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "tab_calculations_browser.py"
).read_text(encoding="utf-8")

_REMOTE_VIEWER_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "tab_remote_archive.py"
).read_text(encoding="utf-8")


def test_settings_ui_exposes_and_persists_viewer_controls():
    for widget_name in (
        "viewer_representation_dropdown",
        "viewer_atom_scale_slider",
        "viewer_bond_radius_slider",
        "viewer_multiple_bonds_toggle",
        "viewer_depth_fog_toggle",
        "viewer_ambient_occlusion_toggle",
    ):
        assert f"{widget_name} = widgets." in _SETTINGS_UI_SOURCE

    for settings_key in (
        "'representation': str(viewer_representation_dropdown.value)",
        "'atom_scale': float(viewer_atom_scale_slider.value)",
        "'bond_radius': float(viewer_bond_radius_slider.value)",
        "'multiple_bonds': bool(viewer_multiple_bonds_toggle.value)",
        "'depth_fog': bool(viewer_depth_fog_toggle.value)",
        "'ambient_occlusion': bool(viewer_ambient_occlusion_toggle.value)",
    ):
        assert settings_key in _SETTINGS_UI_SOURCE
    assert "XYZ/ORCA coordinates do not" in _SETTINGS_UI_SOURCE
    assert "only affects MOL/SDF" in _SETTINGS_UI_SOURCE


def test_viewer_representation_styles_use_global_dimensions():
    ball_and_stick = _MODULE.build_molecule_view_style(
        "ball_and_stick",
        atom_scale=0.42,
        bond_radius=0.17,
        multiple_bonds=True,
    )
    assert set(ball_and_stick) == {"stick", "sphere"}
    assert ball_and_stick["sphere"]["scale"] == 0.42
    assert ball_and_stick["stick"]["radius"] == 0.17
    assert ball_and_stick["stick"]["singleBonds"] is False

    stick = _MODULE.build_molecule_view_style(
        "stick",
        atom_scale=0.9,
        bond_radius=0.23,
        multiple_bonds=False,
    )
    assert set(stick) == {"stick"}
    assert stick["stick"]["radius"] == 0.23
    assert stick["stick"]["singleBonds"] is True

    sphere = _MODULE.build_molecule_view_style("sphere", atom_scale=1.0)
    assert set(sphere) == {"sphere"}
    assert sphere["sphere"]["scale"] == 1.0

    line = _MODULE.build_molecule_view_style("line", bond_radius=0.22)
    assert set(line) == {"line"}
    assert line["line"]["linewidth"] == 2.0


def test_overlay_style_keeps_representation_and_replaces_element_colors():
    style = _MODULE.build_molecule_view_style(
        "ball_and_stick",
        atom_scale=0.36,
        bond_radius=0.15,
    )

    colored = json.loads(
        _MODULE.molecule_view_style_js(style, color="#1f5fff")
    )

    assert set(colored) == {"stick", "sphere"}
    assert colored["stick"]["color"] == "#1f5fff"
    assert colored["sphere"]["color"] == "#1f5fff"
    assert "colorscheme" not in colored["stick"]
    assert colored["stick"]["radius"] == 0.15
    assert colored["sphere"]["scale"] == 0.36


def test_quality_profiles_control_renderer_independently():
    low = _MODULE.VIEWER_QUALITY_PROFILES["low"]["viewer_config"]
    medium = _MODULE.VIEWER_QUALITY_PROFILES["medium"]["viewer_config"]
    high = _MODULE.VIEWER_QUALITY_PROFILES["high"]["viewer_config"]

    assert low["antialias"] is False
    assert medium["antialias"] is True
    assert low["upscale"] is False
    assert "upscale" not in medium
    assert "upscale" not in high
    assert "style" not in medium
    assert high["antialias"] is True
    assert "style" not in high

    assert _MODULE.build_viewer_config("high", depth_fog=True)["disableFog"] is False
    assert _MODULE.build_viewer_config("high", depth_fog=False)["disableFog"] is True
    ao = _MODULE.build_viewer_config(
        "high",
        depth_fog=True,
        ambient_occlusion=True,
    )
    assert ao["style"] == "ambientOcclusion"


def test_py3dmol_style_application_injects_renderer_config(monkeypatch):
    profile = {
        "style": {"sphere": {"colorscheme": "Jmol", "scale": 0.8}},
        "viewer_config": {
            "backgroundColor": "white",
            "antialias": False,
        },
        "viewer_config_js": '{"backgroundColor":"white","antialias":false}',
    }
    monkeypatch.setattr(_MODULE, "get_viewer_profile", lambda: profile)

    class FakeView:
        def __init__(self):
            self.startjs = (
                'viewer_UNIQUEID = $3Dmol.createViewer(el,'
                '{backgroundColor:"white"});'
            )
            self.style = None
            self.background = None

        def setStyle(self, selection, style):
            self.style = (selection, style)

        def setBackgroundColor(self, color):
            self.background = color

        def zoomTo(self):
            pass

        def center(self):
            pass

        def zoom(self, value):
            pass

        def render(self):
            pass

    view = FakeView()
    _MODULE.apply_molecule_view_style(view)

    assert '{"backgroundColor":"white","antialias":false}' in view.startjs
    assert view.style == ({}, profile["style"])
    assert view.background == "white"


def test_mouse_patch_coalesces_drag_rendering_and_disposes_old_contexts():
    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS

    assert "pendingDx += dx" in patch
    assert "pendingMoveEvent=e" in patch
    assert "firstMovePending" in patch
    assert "now-moveFrameStarted>34" in patch
    assert "requestAnimationFrame" in patch
    assert "moveFrame=requestFrame(flushTranslation)" in patch
    assert 'viewer.translate(dx, dy);\nif(typeof viewer.render' not in patch
    assert "viewer.__delfinInteracting=true" in patch
    assert "useFastInteractionStyle" in patch
    assert 'style.replace("ambientOcclusion","")' in patch
    assert "__delfinInteractionEndHandlers" in patch
    assert "window.__delfinDisposeViewer" in patch
    assert "WEBGL_lose_context" in patch
    assert "window.setInterval" not in patch


def test_mouse_patch_leaves_editor_drag_handlers_on_window_alive():
    """The Submit-tab editor binds its drag handlers on window in the capture
    phase, same as this patch. Using stopImmediatePropagation there would kill
    them and disable rubber-band selection and atom dragging."""
    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS

    window_handlers = patch.split('el.addEventListener("contextmenu"')[0]
    move_handler = window_handlers.split('var onMouseMove = function(e){')[1]
    assert 'stopImmediatePropagation' not in move_handler
    assert 'softStopEvt(e)' in move_handler

    up_handler = window_handlers.split('var onMouseUp = function(e){')[1]
    assert 'stopImmediatePropagation' not in up_handler
    assert 'if(wasDragging) softStopEvt(e);' in patch

    # softStopEvt must stop descent to 3Dmol's own listeners but nothing else.
    soft = patch.split('var softStopEvt = function(e){')[1].split('};')[0]
    assert 'e.stopPropagation()' in soft
    assert 'stopImmediatePropagation' not in soft

    # While the editor owns the drag, 3Dmol must not be fed forwarded moves.
    assert 'if(editorDragActive()) return;' in patch
    assert '&&!editorDragActive()' in patch
    assert 'window._submitManipStateByScope' in patch

    # An already-open dashboard rebinds only when the version changes.
    assert 'var PATCH_VERSION=9;' in patch


def test_py3dmol_viewers_release_the_contexts_they_leave_behind():
    """py3Dmol builds a fresh container and a fresh window.viewer_<id> on every
    render, so the caller has no previous handle to release. Browsers cap live
    WebGL contexts and evict the oldest, which blanks viewers in other tabs."""
    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS
    sweep = patch.split('window.__delfinDisposeOrphanedViewers = function(){')[1]
    sweep = sweep.split('window.__delfinPatchAllKnown3DmolViewers')[0]
    # Only viewers whose element has left the document may be released.
    assert 'document.body.contains(el)' in sweep
    assert 'window.__delfinDisposeViewer(v)' in sweep
    assert 'v.__delfinDisposed' in sweep

    submit = (
        Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_submit.py"
    ).read_text(encoding="utf-8")
    # The Submit tab keeps a scope-keyed handle, so it can release directly.
    assert '__delfinDisposeViewer(prev)' in submit


def test_deferred_recentering_never_overrides_the_user():
    """apply_molecule_view_style recenters again at +120 ms because the canvas
    is not always at its final size on the first frame. Any pan or zoom the
    user managed in that window used to be thrown away."""

    class FakeView:
        def __init__(self):
            self.startjs = '/*head*/'

        def setStyle(self, *a):
            pass

        def setBackgroundColor(self, *a):
            pass

        def zoomTo(self):
            pass

        def center(self):
            pass

        def zoom(self, *a):
            pass

        def render(self):
            pass

    view = FakeView()
    _MODULE.apply_molecule_view_style(view)
    recenter = view.startjs.split('var __delfinRecenter=function(){')[1]
    guard = 'if(viewer_UNIQUEID.__delfinUserInteracted) return;'
    assert guard in recenter
    # The guard has to come before any camera call, not after.
    body = recenter.split('};')[0]
    assert body.index(guard) < body.index('zoomTo()')

    # The flag is sticky, unlike __delfinInteracting which clears on release.
    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS
    start = patch.split('var notifyInteractionStart = function(){')[1].split('};')[0]
    assert 'viewer.__delfinUserInteracted=true;' in start
    assert start.index('__delfinUserInteracted') < start.index(
        'if(viewer.__delfinInteracting) return;')


def test_orca_label_occlusion_is_deferred_until_mouse_interaction_ends():
    assert "grid=Object.create(null)" in _ORCA_VIEWER_SOURCE
    assert "hideForInteraction" not in _ORCA_VIEWER_SOURCE
    assert "__delfinInteractionEndHandlers" in _ORCA_VIEWER_SOURCE
    assert "schedule(0)" in _ORCA_VIEWER_SOURCE
    assert "v.render=function()" not in _ORCA_VIEWER_SOURCE
    assert "raf(loop)" not in _ORCA_VIEWER_SOURCE
    assert "window.__delfinDisposeViewer(prev)" in _ORCA_VIEWER_SOURCE


def test_trajectory_playback_yields_rendering_while_viewer_is_dragged():
    for source in (_CALC_VIEWER_SOURCE, _REMOTE_VIEWER_SOURCE):
        assert "viewer && !viewer.__delfinInteracting" in source
        assert "state.get('traj_playing')" in source or 'state.get("traj_playing")' in source
        assert "window.__delfinDisposeViewer(previousViewer)" in source


def test_3dmol_is_vendored_and_satisfies_both_loader_guards():
    """py3Dmol's loader fetches jsdelivr and guards on $3Dmolpromise; the
    hand-written viewers fetch 3Dmol.org and guard on $3Dmol. The guards never
    see each other, so a page mixing tabs downloaded the library twice and the
    second copy replaced the global for every later viewer. Without outbound
    network no molecule appeared at all, and nothing said why."""
    bundle = (
        Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "static" / "3Dmol-min.js"
    )
    assert bundle.is_file()
    assert bundle.stat().st_size > 100_000
    assert bundle.with_suffix(".js.LICENSE.txt").is_file()

    js = _MODULE.vendored_3dmol_js()
    assert 'typeof __delfinGlobal.$3Dmol === "undefined"' in js
    assert "__delfinGlobal.$3Dmolpromise = __delfinGlobal.$3Dmolpromise" in js

    # It has to run before anything that creates a viewer.
    dashboard = (
        Path(__file__).resolve().parents[1] / "delfin" / "dashboard" / "__init__.py"
    ).read_text(encoding="utf-8")
    assert "vendored_3dmol_js()" in dashboard
    init_list = dashboard.split("_calc_init = ")[1].split("])")[0]
    assert init_list.index("vendored_3dmol_js()") < init_list.index(
        "RIGHT_MOUSE_TRANSLATE_PATCH_JS"
    )

    # And it has to reach the wheel, or every viewer quietly returns to the CDNs.
    pyproject = (
        Path(__file__).resolve().parents[1] / "pyproject.toml"
    ).read_text(encoding="utf-8")
    assert '"delfin.dashboard" = ["static/*.js", "static/*.txt"]' in pyproject


def test_quality_levels_render_at_different_pixel_densities():
    """Medium and High used to be pixel-identical: they differed only in
    cartoonQuality, which is dead because DELFIN never renders cartoons, and
    3Dmol's own `upscale` flag merely raises devicePixelRatio to 2. The one
    real lever is the pixel ratio the viewer is constructed with."""
    ratios = _MODULE.VIEWER_QUALITY_PIXEL_RATIO
    assert ratios["low"] < ratios["medium"] < ratios["high"]
    assert _MODULE.get_viewer_profile()["pixel_ratio"] in ratios.values()

    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS
    # 3Dmol reads window.devicePixelRatio once, at construction. The factory
    # itself cannot be wrapped -- $3Dmol.createViewer is a non-configurable
    # getter, so assigning over it silently does nothing and reports success --
    # so every call site goes through this helper instead.
    assert "window.__delfinWithPixelRatio" in patch
    assert "window.__delfinCreateViewer = function" in patch
    assert "window.__delfinViewerPixelRatio" in patch

    for name in (
        "tab_calculations_browser.py",
        "tab_remote_archive.py",
        "tab_orca_builder.py",
    ):
        source = (
            Path(__file__).resolve().parents[1] / "delfin" / "dashboard" / name
        ).read_text(encoding="utf-8")
        assert "$3Dmol.createViewer(" not in source, name
        assert "window.__delfinCreateViewer(" in source, name

    # py3Dmol builds its own construction call, so it is rerouted in startjs.
    viewer_source = _MODULE_PATH.read_text(encoding="utf-8")
    assert "if '$3Dmol.createViewer(' in view.startjs:" in viewer_source
    assert "'window.__delfinCreateViewer'," in viewer_source or (
        "'window.__delfinCreateViewer(',"
    ) in viewer_source
    # The override must always be undone, whatever the build does.
    helper = patch.split("window.__delfinWithPixelRatio = function")[1]
    assert "finally" in helper.split("};")[0] or "finally" in helper[:900]


def test_png_export_re_renders_instead_of_upscaling_the_screen():
    """Export took the on-screen canvas and blew it up 6x with drawImage --
    bilinear interpolation, so a tens-of-megabyte file that is blurrier than
    what the user is looking at."""
    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS
    assert "window.__delfinRenderViewerPng" in patch

    download = patch.split("window.__delfinDownloadViewerPng = function")[1]
    assert "__delfinRenderViewerPng(viewer, opts)" in download
    assert "__delfinCloneCanvasDataUrl" not in download

    render = patch.split("window.__delfinRenderViewerPng = function")[1].split(
        "window.__delfinDownloadViewerPng"
    )[0]
    # It renders the scene again at a higher pixel density...
    assert "window.__delfinWithPixelRatio(scale" in render
    # ...from the geometry on screen, which may be an edit that exists nowhere
    # else, and matched to the on-screen camera.
    assert "selectedAtoms({})" in render
    assert "shot.setView(viewer.getView())" in render
    # The off-screen viewer must not leak a WebGL context.
    assert "__delfinDisposeViewer(shot)" in render
    assert "removeChild(host)" in render


def test_viewer_creation_is_observable_without_patching_the_factory():
    """The Calculations browser carried a monkey-patch of $3Dmol.createViewer
    meant to kick every scoped resizer 300 ms after a viewer appeared. It never
    ran once: createViewer is a non-configurable getter, so the assignment
    fails silently in sloppy mode — no exception, so the defineProperty
    fallback inside the catch was never reached either, and the 'patched' flag
    was set regardless."""
    patch = _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS
    assert "window.__delfinOnViewerCreated" in patch
    assert "window.__delfinViewerCreatedHooks" in patch
    funnel = patch.split("window.__delfinCreateViewer = function")[1][:1400]
    assert "hooks[i](viewer, element)" in funnel

    calc = (
        Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_calculations_browser.py"
    ).read_text(encoding="utf-8")
    assert "patchCreateViewer" not in calc
    assert "_calcResizePatched" not in calc
    assert "window.__delfinOnViewerCreated(function()" in calc


def test_resize_handlers_are_coalesced():
    """They are called from viewer creation, the splitter, several timers and a
    MutationObserver watching every inline-style change in the tab, and each
    call reallocates the renderer's buffers and draws three times."""
    assert "window.__delfinCoalesce" in _MODULE.RIGHT_MOUSE_TRANSLATE_PATCH_JS
    assert "window.__delfinCoalesce(window[" in _CALC_VIEWER_SOURCE
    assert "__delfinCoalesced" in _REMOTE_VIEWER_SOURCE


def test_viewers_no_longer_give_up_while_their_tab_is_hidden():
    """80 tries at 50 ms was a hard four-second deadline, and the poll skips a
    host whose tab is not the visible one — so switching away for four seconds
    lost the molecule with no error anywhere."""
    for source in (_CALC_VIEWER_SOURCE, _REMOTE_VIEWER_SOURCE, _ORCA_VIEWER_SOURCE):
        assert "tries < 80" not in source
        assert "tries<80" not in source
    assert "tries < 400" in _CALC_VIEWER_SOURCE
    assert "tries < 40 ? 50 : 250" in _CALC_VIEWER_SOURCE
    assert "tries<400" in _ORCA_VIEWER_SOURCE


def test_vendored_bundle_does_not_depend_on_the_host_calling_convention():
    """The bundle is a UMD wrapper: its first statement is ``root = this`` and
    its last assigns ``root["3Dmol"]``. ctx.run_js hands scripts to
    JupyterLab/Voila through display(Javascript(...)), which evaluates them
    with ``this`` undefined, so the library threw

        Cannot set properties of undefined (setting '3Dmol')

    and — because that abort killed the rest of the start-up script — neither
    $3Dmolpromise nor the mouse patch that defines __delfinCreateViewer ever
    ran. Every viewer construction then failed and no molecule appeared."""
    js = _MODULE.vendored_3dmol_js()
    assert js.startswith("(function(__delfinGlobal){")
    # The library is invoked with an explicit receiver, not the host's `this`.
    assert ").call(__delfinGlobal);" in js
    assert 'typeof window !== "undefined" ? window' in js
    # And the guard the loaders check is set on that same object.
    assert "__delfinGlobal.$3Dmolpromise = __delfinGlobal.$3Dmolpromise" in js
