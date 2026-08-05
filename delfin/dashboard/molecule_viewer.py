"""3D molecule visualisation helpers using py3Dmol."""

import json

import py3Dmol
from IPython.display import HTML, clear_output, display

DEFAULT_3DMOL_BACKGROUND = 'white'
DEFAULT_3DMOL_ZOOM = 0.90
DEFAULT_3DMOL_REPRESENTATION = 'ball_and_stick'
DEFAULT_3DMOL_ATOM_SCALE = 0.28
DEFAULT_3DMOL_BOND_RADIUS = 0.11
DEFAULT_3DMOL_MULTIPLE_BONDS = True
DEFAULT_3DMOL_DEPTH_FOG = True
DEFAULT_3DMOL_AMBIENT_OCCLUSION = False

VIEWER_REPRESENTATIONS = {
    'line',
    'stick',
    'ball_and_stick',
    'sphere',
}


def build_molecule_view_style(
    representation=DEFAULT_3DMOL_REPRESENTATION,
    atom_scale=DEFAULT_3DMOL_ATOM_SCALE,
    bond_radius=DEFAULT_3DMOL_BOND_RADIUS,
    multiple_bonds=DEFAULT_3DMOL_MULTIPLE_BONDS,
    *,
    color=None,
):
    """Build one 3Dmol atom style from the global display controls.

    ``atom_scale`` scales van-der-Waals radii (1.0 is a space-filling model),
    while ``bond_radius`` is the cylinder radius in Angstrom. A fixed ``color``
    is used by overlay viewers; otherwise the standard Jmol element palette is
    retained.
    """
    representation = str(representation or '').strip().lower()
    if representation not in VIEWER_REPRESENTATIONS:
        representation = DEFAULT_3DMOL_REPRESENTATION
    atom_scale = max(0.05, min(1.50, float(atom_scale)))
    bond_radius = max(0.02, min(0.50, float(bond_radius)))
    color_spec = {'color': str(color)} if color else {'colorscheme': 'Jmol'}

    line = dict(color_spec)
    # Browser/WebGL support for wide native lines varies, but passing the
    # 3Dmol linewidth still honors the thickness control where supported.
    line['linewidth'] = round(max(1.0, bond_radius / DEFAULT_3DMOL_BOND_RADIUS), 3)
    stick = {
        **color_spec,
        'radius': bond_radius,
        'singleBonds': not bool(multiple_bonds),
        'doubleBondScaling': 0.65,
        'tripleBondScaling': 0.65,
    }
    sphere = {**color_spec, 'scale': atom_scale}

    if representation == 'line':
        return {'line': line}
    if representation == 'stick':
        return {'stick': stick}
    if representation == 'sphere':
        return {'sphere': sphere}
    return {'stick': stick, 'sphere': sphere}


def molecule_view_style_js(style, *, color=None):
    """Serialize a style for inline JavaScript, optionally recoloring it."""
    rendered = {}
    for name, options in dict(style or {}).items():
        rendered_options = dict(options or {})
        if color:
            rendered_options.pop('colorscheme', None)
            rendered_options['color'] = str(color)
        rendered[name] = rendered_options
    return json.dumps(rendered, separators=(',', ':'))


DEFAULT_3DMOL_STYLE = build_molecule_view_style()
DEFAULT_3DMOL_STYLE_JS = molecule_view_style_js(DEFAULT_3DMOL_STYLE)

_VIEWER_FIXED_WIDTH = 560
_VIEWER_FIXED_HEIGHT = 420
_VIEWER_FIXED_ZOOM = 0.90

# Container size used by tabs that embed the 3D viewer inside an
# ipywidgets/HTML wrapper (calc-browser, remote-archive, submit, fukui).
# Single source of truth so the dashboard stays visually consistent.
VIEWER_CONTAINER_HEIGHT_PX = 450
VIEWER_CONTAINER_DYNAMIC_SCALE = 0.9725

VIEWER_QUALITY_PROFILES = {
    'low': {
        'viewer_config': {
            'backgroundColor': DEFAULT_3DMOL_BACKGROUND,
            'antialias': False,
            'upscale': False,
            'cartoonQuality': 5,
        },
    },
    'medium': {
        'viewer_config': {
            'backgroundColor': DEFAULT_3DMOL_BACKGROUND,
            'antialias': True,
            'cartoonQuality': 10,
        },
    },
    'high': {
        'viewer_config': {
            'backgroundColor': DEFAULT_3DMOL_BACKGROUND,
            'antialias': True,
            'cartoonQuality': 20,
        },
    },
}


def build_viewer_config(
    quality='high',
    depth_fog=DEFAULT_3DMOL_DEPTH_FOG,
    ambient_occlusion=DEFAULT_3DMOL_AMBIENT_OCCLUSION,
):
    """Build the GLViewer renderer config from independent visual controls."""
    quality = str(quality or '').strip().lower()
    if quality not in VIEWER_QUALITY_PROFILES:
        quality = 'high'
    config = dict(VIEWER_QUALITY_PROFILES[quality]['viewer_config'])
    # 3Dmol's fog is the depth cue that blends distant atoms into the
    # background. State it explicitly instead of relying on a library default.
    config['disableFog'] = not bool(depth_fog)
    if ambient_occlusion:
        config.update({
            'style': 'ambientOcclusion',
            'strength': 0.65,
            'radius': 4.0,
        })
    return config


_VIEWER_DISABLED_PLACEHOLDER_HTML = (
    "<div style=\"width:100%;max-width:560px;padding:18px 24px;"
    "border:1px dashed #b0b6bf;background:#f6f7f9;color:#4a525c;"
    "font-family:Arial,sans-serif;font-size:13px;border-radius:6px;\">"
    "The 3D viewer is disabled in Settings."
    "</div>"
)


def get_viewer_profile():
    """Return the complete active global viewer profile.

    Representation, atom size, bond thickness, depth fog and ambient occlusion
    are independent controls. Quality controls antialiasing and supported
    geometry detail. Falls back to the historical ball-and-stick/fog defaults
    if the user settings cannot be loaded.
    """
    try:
        from delfin.user_settings import load_viewer_settings
        viewer = load_viewer_settings()
    except Exception:
        viewer = {
            'enabled': True,
            'quality': 'high',
            'representation': DEFAULT_3DMOL_REPRESENTATION,
            'atom_scale': DEFAULT_3DMOL_ATOM_SCALE,
            'bond_radius': DEFAULT_3DMOL_BOND_RADIUS,
            'multiple_bonds': DEFAULT_3DMOL_MULTIPLE_BONDS,
            'depth_fog': DEFAULT_3DMOL_DEPTH_FOG,
            'ambient_occlusion': DEFAULT_3DMOL_AMBIENT_OCCLUSION,
        }
    quality = viewer.get('quality', 'high')
    if quality not in VIEWER_QUALITY_PROFILES:
        quality = 'high'
    representation = str(
        viewer.get('representation', DEFAULT_3DMOL_REPRESENTATION)
    )
    if representation not in VIEWER_REPRESENTATIONS:
        representation = DEFAULT_3DMOL_REPRESENTATION
    try:
        atom_scale = float(viewer.get('atom_scale', DEFAULT_3DMOL_ATOM_SCALE))
    except (TypeError, ValueError):
        atom_scale = DEFAULT_3DMOL_ATOM_SCALE
    try:
        bond_radius = float(viewer.get('bond_radius', DEFAULT_3DMOL_BOND_RADIUS))
    except (TypeError, ValueError):
        bond_radius = DEFAULT_3DMOL_BOND_RADIUS
    atom_scale = max(0.05, min(1.50, atom_scale))
    bond_radius = max(0.02, min(0.50, bond_radius))
    multiple_bonds = bool(
        viewer.get('multiple_bonds', DEFAULT_3DMOL_MULTIPLE_BONDS)
    )
    depth_fog = bool(viewer.get('depth_fog', DEFAULT_3DMOL_DEPTH_FOG))
    ambient_occlusion = bool(
        viewer.get('ambient_occlusion', DEFAULT_3DMOL_AMBIENT_OCCLUSION)
    )
    style = build_molecule_view_style(
        representation,
        atom_scale,
        bond_radius,
        multiple_bonds,
    )
    viewer_config = build_viewer_config(
        quality,
        depth_fog,
        ambient_occlusion,
    )
    return {
        'enabled': bool(viewer.get('enabled', True)),
        'quality': quality,
        'representation': representation,
        'atom_scale': atom_scale,
        'bond_radius': bond_radius,
        'multiple_bonds': multiple_bonds,
        'depth_fog': depth_fog,
        'ambient_occlusion': ambient_occlusion,
        'width': _VIEWER_FIXED_WIDTH,
        'height': _VIEWER_FIXED_HEIGHT,
        'zoom': _VIEWER_FIXED_ZOOM,
        'style': style,
        'style_js': molecule_view_style_js(style),
        'viewer_config': viewer_config,
        'viewer_config_js': json.dumps(viewer_config, separators=(',', ':')),
    }


def viewer_disabled_html(width=None, height=None):
    """Return the placeholder HTML shown when the viewer is disabled."""
    return _VIEWER_DISABLED_PLACEHOLDER_HTML
_VENDORED_3DMOL_CACHE = None


def vendored_3dmol_js():
    """Return JS that installs the bundled 3Dmol build, or '' if it is missing.

    Every viewer used to pull 3Dmol over the network, and from two different
    origins: py3Dmol's generated loader fetches jsdelivr and guards on
    ``$3Dmolpromise``, while the hand-written viewers fetch 3Dmol.org and guard
    on ``$3Dmol``. The guards do not see each other, so a page mixing tabs
    downloaded the library twice and whichever landed second replaced the
    global for every later viewer. On a cluster without outbound network the
    molecule simply never appeared, with no error shown anywhere.

    Running this once at start-up satisfies *both* guards: the library is
    already defined, and ``$3Dmolpromise`` is an already-resolved promise, so
    neither loader fetches anything.
    """
    global _VENDORED_3DMOL_CACHE
    if _VENDORED_3DMOL_CACHE is not None:
        return _VENDORED_3DMOL_CACHE
    try:
        from importlib.resources import files
        source = (
            files('delfin.dashboard')
            .joinpath('static/3Dmol-min.js')
            .read_text(encoding='utf-8')
        )
    except Exception:
        _VENDORED_3DMOL_CACHE = ''
        return ''
    _VENDORED_3DMOL_CACHE = (
        'if (typeof $3Dmol === "undefined") {\n'
        + source
        + '\n}\n'
        'try { window.$3Dmolpromise = window.$3Dmolpromise '
        '|| Promise.resolve(window.$3Dmol); } catch (e) {}\n'
    )
    return _VENDORED_3DMOL_CACHE


RIGHT_MOUSE_TRANSLATE_PATCH_JS = (
    '(function(){\n'
    'var PATCH_VERSION=8;\n'
    'if(window.__delfinRightDragTranslateVersion===PATCH_VERSION) return;\n'
    'if(window.__delfinRightDragTranslateTimer){\n'
    'try{clearInterval(window.__delfinRightDragTranslateTimer);}catch(e){}\n'
    'window.__delfinRightDragTranslateTimer=null;\n'
    '}\n'
    'window.__delfinRightDragTranslateSetup = true;\n'
    'window.__delfinRightDragTranslateVersion = PATCH_VERSION;\n'
    'window.__delfinResolveViewerElement = function(viewer, fallbackEl){\n'
    'try {\n'
    'if (fallbackEl && fallbackEl.nodeType === 1) return fallbackEl;\n'
    'if (viewer && viewer.__delfinRootElement && viewer.__delfinRootElement.nodeType === 1) '
    'return viewer.__delfinRootElement;\n'
    'if (viewer && viewer.container && viewer.container.nodeType === 1) return viewer.container;\n'
    '} catch(e) {}\n'
    'return null;\n'
    '};\n'
    'window.__delfinEnableRightDragTranslate = function(viewer, rootEl){\n'
    'try {\n'
    'if(!viewer) return false;\n'
    'var el = window.__delfinResolveViewerElement(viewer, rootEl);\n'
    'if(!el) return false;\n'
    'if(viewer.__delfinRightDragTranslateBound '
    '&& viewer.__delfinRightDragTranslateVersion===PATCH_VERSION) return true;\n'
    'if(viewer.__delfinRightDragTranslateCleanup){\n'
    'try{viewer.__delfinRightDragTranslateCleanup();}catch(e){}\n'
    '}\n'
    'viewer.__delfinRightDragTranslateBound = true;\n'
    'viewer.__delfinRightDragTranslateVersion = PATCH_VERSION;\n'
    'viewer.__delfinRootElement = el;\n'
    'var dragging = false;\n'
    'var mouseActive = false;\n'
    'var lastX = 0;\n'
    'var lastY = 0;\n'
    'var pendingDx = 0;\n'
    'var pendingDy = 0;\n'
    'var pendingMoveEvent = null;\n'
    'var moveFrame = null;\n'
    'var moveFrameStarted = 0;\n'
    'var firstMovePending = false;\n'
    'var wheelTimer = 0;\n'
    'var restoreViewStyle = null;\n'
    'var requestFrame = typeof window.requestAnimationFrame==="function" '
    '?function(fn){return window.requestAnimationFrame(fn);}'
    ':function(fn){return window.setTimeout(fn,16);};\n'
    'var cancelFrame = typeof window.cancelAnimationFrame==="function" '
    '?function(id){window.cancelAnimationFrame(id);}'
    ':function(id){window.clearTimeout(id);};\n'
    'var nowMs=function(){\n'
    'return window.performance&&typeof window.performance.now==="function"'
    '?window.performance.now():Date.now();\n'
    '};\n'
    'var useFastInteractionStyle = function(){\n'
    'if(restoreViewStyle||typeof viewer.setViewStyle!=="function") return;\n'
    'var config=null;\n'
    'if(typeof viewer.getConfig==="function"){try{config=viewer.getConfig();}catch(e){}}\n'
    'var style=config&&config.style?String(config.style):"";\n'
    'if(style.indexOf("ambientOcclusion")<0) return;\n'
    'restoreViewStyle={style:style,strength:config.strength,radius:config.radius};\n'
    'try{viewer.setViewStyle({style:style.replace("ambientOcclusion","")});}catch(e){}\n'
    '};\n'
    'var restoreFullViewStyle = function(){\n'
    'if(!restoreViewStyle||typeof viewer.setViewStyle!=="function") return false;\n'
    'var style=restoreViewStyle;restoreViewStyle=null;\n'
    'try{viewer.setViewStyle(style);return true;}catch(e){return false;}\n'
    '};\n'
    'var notifyInteractionStart = function(){\n'
    # Sticky, unlike __delfinInteracting: once the user has touched this viewer,
    # deferred set-up work must not reset the camera under their hands.
    'viewer.__delfinUserInteracted=true;\n'
    'if(viewer.__delfinInteracting) return;\n'
    'viewer.__delfinInteracting=true;\n'
    'useFastInteractionStyle();\n'
    'var handlers=viewer.__delfinInteractionStartHandlers||[];\n'
    'handlers.slice().forEach(function(fn){\n'
    'if(typeof fn==="function"){try{fn();}catch(e){}}\n'
    '});\n'
    '};\n'
    'var notifyInteractionEnd = function(){\n'
    'if(mouseActive||dragging) return;\n'
    'viewer.__delfinInteracting=false;\n'
    'var styleRestored=restoreFullViewStyle();\n'
    'var handlers=viewer.__delfinInteractionEndHandlers||[];\n'
    'handlers.slice().forEach(function(fn){\n'
    'if(typeof fn==="function"){try{fn();}catch(e){}}\n'
    '});\n'
    'if(styleRestored&&handlers.length===0&&typeof viewer.render==="function"){\n'
    'try{viewer.render();}catch(e){}\n'
    '}\n'
    '};\n'
    'var translateNow = function(dx,dy){\n'
    'if(typeof viewer.translate === "function"&&(dx||dy)) {\n'
    'viewer.translate(dx, dy);\n'
    '} else if(typeof viewer.translateScene === "function"&&(dx||dy)) {\n'
    'viewer.translateScene(dx, dy);\n'
    '} else if(typeof viewer.pan === "function"&&(dx||dy)) {\n'
    'viewer.pan(dx, dy);\n'
    '}\n'
    '};\n'
    'var flushTranslation = function(){\n'
    'moveFrame=null;\n'
    'moveFrameStarted=0;\n'
    'var dx=pendingDx,dy=pendingDy;\n'
    'var moveEvent=pendingMoveEvent;\n'
    'pendingDx=0;pendingDy=0;\n'
    'pendingMoveEvent=null;\n'
    'try {\n'
    'if(dragging) {\n'
    'translateNow(dx,dy);\n'
    '} else if(!dragging&&moveEvent&&typeof viewer._handleMouseMove==="function"'
    '&&!editorDragActive()) {\n'
    'viewer._handleMouseMove(moveEvent);\n'
    '}\n'
    '} catch(_e) {}\n'
    '};\n'
    'var scheduleTranslation = function(){\n'
    'var now=nowMs();\n'
    'if(moveFrame!==null&&now-moveFrameStarted>34){\n'
    'cancelFrame(moveFrame);moveFrame=null;flushTranslation();return;\n'
    '}\n'
    'if(moveFrame===null){\n'
    'moveFrameStarted=now;moveFrame=requestFrame(flushTranslation);\n'
    '}\n'
    '};\n'
    'var stopEvt = function(e){\n'
    'if(!e) return;\n'
    'if(typeof e.preventDefault === "function") e.preventDefault();\n'
    'if(typeof e.stopPropagation === "function") e.stopPropagation();\n'
    'if(typeof e.stopImmediatePropagation === "function") e.stopImmediatePropagation();\n'
    '};\n'
    # Window-level handlers must not use stopImmediatePropagation: the editor
    # overlay binds its own drag handlers on window in the same capture phase,
    # and killing sibling listeners there would disable rubber-band selection
    # and atom dragging. stopPropagation still keeps the event from descending
    # to 3Dmol's own element listeners, which is all the throttling needs.
    'var softStopEvt = function(e){\n'
    'if(!e) return;\n'
    'if(typeof e.preventDefault === "function") e.preventDefault();\n'
    'if(typeof e.stopPropagation === "function") e.stopPropagation();\n'
    '};\n'
    'var editorDragActive = function(){\n'
    'try {\n'
    'var states = window._submitManipStateByScope;\n'
    'if(!states) return false;\n'
    'for (var k in states) {\n'
    'if(!Object.prototype.hasOwnProperty.call(states, k)) continue;\n'
    'if(states[k] && states[k].drag) return true;\n'
    '}\n'
    '} catch(e) {}\n'
    'return false;\n'
    '};\n'
    'var onContextMenu = function(e){\n'
    'stopEvt(e);\n'
    'return false;\n'
    '};\n'
    'var onMouseDown = function(e){\n'
    'if(!e) return;\n'
    'mouseActive=true;\n'
    'firstMovePending=true;\n'
    'notifyInteractionStart();\n'
    'if(e.button !== 2) return;\n'
    'dragging = true;\n'
    'lastX = e.clientX || 0;\n'
    'lastY = e.clientY || 0;\n'
    'stopEvt(e);\n'
    '};\n'
    'var onMouseMove = function(e){\n'
    'if(!mouseActive) return;\n'
    'if(!dragging){\n'
    'if(editorDragActive()) return;\n'
    'if(typeof viewer._handleMouseMove!=="function") return;\n'
    'if(firstMovePending){\n'
    'firstMovePending=false;\n'
    'try{viewer._handleMouseMove(e);}catch(_e){}\n'
    'softStopEvt(e);return;\n'
    '}\n'
    'pendingMoveEvent=e;\n'
    'scheduleTranslation();\n'
    'softStopEvt(e);\n'
    'return;\n'
    '}\n'
    'var x = e.clientX || 0;\n'
    'var y = e.clientY || 0;\n'
    'var dx = x - lastX;\n'
    'var dy = y - lastY;\n'
    'lastX = x;\n'
    'lastY = y;\n'
    'pendingDx += dx;\n'
    'pendingDy -= dy;\n'
    'if(firstMovePending){\n'
    'firstMovePending=false;\n'
    'var firstDx=pendingDx,firstDy=pendingDy;pendingDx=0;pendingDy=0;\n'
    'try{translateNow(firstDx,firstDy);}catch(_e){}\n'
    'softStopEvt(e);return;\n'
    '}\n'
    'scheduleTranslation();\n'
    'softStopEvt(e);\n'
    '};\n'
    'var onMouseUp = function(e){\n'
    'var wasDragging=dragging;\n'
    'mouseActive=false;\n'
    'firstMovePending=false;\n'
    'if(moveFrame!==null){cancelFrame(moveFrame);moveFrame=null;flushTranslation();}\n'
    'dragging = false;\n'
    'notifyInteractionEnd();\n'
    'if(wasDragging) softStopEvt(e);\n'
    '};\n'
    'var onWheel = function(){\n'
    'notifyInteractionStart();\n'
    'if(wheelTimer) window.clearTimeout(wheelTimer);\n'
    'wheelTimer=window.setTimeout(function(){\n'
    'wheelTimer=0;notifyInteractionEnd();\n'
    '},120);\n'
    '};\n'
    'var onTouchStart = function(){\n'
    'mouseActive=true;notifyInteractionStart();\n'
    '};\n'
    'var onTouchEnd = function(){\n'
    'mouseActive=false;notifyInteractionEnd();\n'
    '};\n'
    'var onBlur = function(){\n'
    'mouseActive=false;dragging=false;firstMovePending=false;notifyInteractionEnd();\n'
    '};\n'
    'el.addEventListener("contextmenu", onContextMenu, true);\n'
    'el.addEventListener("mousedown", onMouseDown, true);\n'
    'el.addEventListener("wheel", onWheel, true);\n'
    'el.addEventListener("touchstart", onTouchStart, true);\n'
    'window.addEventListener("mousemove", onMouseMove, true);\n'
    'window.addEventListener("mouseup", onMouseUp, true);\n'
    'window.addEventListener("touchend", onTouchEnd, true);\n'
    'window.addEventListener("touchcancel", onTouchEnd, true);\n'
    'window.addEventListener("blur", onBlur, true);\n'
    'viewer.__delfinRightDragTranslateCleanup = function(){\n'
    'try {\n'
    'el.removeEventListener("contextmenu", onContextMenu, true);\n'
    'el.removeEventListener("mousedown", onMouseDown, true);\n'
    'el.removeEventListener("wheel", onWheel, true);\n'
    'el.removeEventListener("touchstart", onTouchStart, true);\n'
    'window.removeEventListener("mousemove", onMouseMove, true);\n'
    'window.removeEventListener("mouseup", onMouseUp, true);\n'
    'window.removeEventListener("touchend", onTouchEnd, true);\n'
    'window.removeEventListener("touchcancel", onTouchEnd, true);\n'
    'window.removeEventListener("blur", onBlur, true);\n'
    'if(moveFrame!==null) cancelFrame(moveFrame);\n'
    'if(wheelTimer) window.clearTimeout(wheelTimer);\n'
    'restoreFullViewStyle();\n'
    '} catch(e) {}\n'
    'moveFrame=null;moveFrameStarted=0;wheelTimer=0;firstMovePending=false;\n'
    'pendingDx=0;pendingDy=0;pendingMoveEvent=null;\n'
    'viewer.__delfinRightDragTranslateBound=false;\n'
    '};\n'
    'return true;\n'
    '} catch(e) {\n'
    'return false;\n'
    '}\n'
    '};\n'
    'window.__delfinDisposeViewer = function(viewer){\n'
    'try {\n'
    'if(!viewer||viewer.__delfinDisposed) return;\n'
    'viewer.__delfinDisposed=true;\n'
    'var el=window.__delfinResolveViewerElement(viewer,null);\n'
    'if(typeof viewer.__delfinCleanup==="function"){\n'
    'try{viewer.__delfinCleanup();}catch(e){}\n'
    '}\n'
    'if(typeof viewer.__delfinRightDragTranslateCleanup==="function"){\n'
    'try{viewer.__delfinRightDragTranslateCleanup();}catch(e){}\n'
    '}\n'
    'if(viewer.divwatcher&&typeof viewer.divwatcher.disconnect==="function"){\n'
    'try{viewer.divwatcher.disconnect();}catch(e){}\n'
    '}\n'
    'if(viewer.intwatcher&&typeof viewer.intwatcher.disconnect==="function"){\n'
    'try{viewer.intwatcher.disconnect();}catch(e){}\n'
    '}\n'
    'if(typeof viewer.pauseAnimate==="function"){try{viewer.pauseAnimate();}catch(e){}}\n'
    'if(viewer.spinInterval){try{clearInterval(viewer.spinInterval);}catch(e){}}\n'
    'var canvas=null;\n'
    'if(typeof viewer.getCanvas==="function"){try{canvas=viewer.getCanvas();}catch(e){}}\n'
    'if(!canvas&&el&&el.querySelector)canvas=el.querySelector("canvas");\n'
    'if(canvas){try{canvas._3dmol_viewer=null;}catch(e){}}\n'
    'var gl=null;\n'
    'if(canvas){\n'
    'try{gl=canvas.getContext("webgl2")||canvas.getContext("webgl")'
    '||canvas.getContext("experimental-webgl");}catch(e){}\n'
    '}\n'
    'if(gl&&typeof gl.getExtension==="function"){\n'
    'var ext=gl.getExtension("WEBGL_lose_context");\n'
    'if(ext&&typeof ext.loseContext==="function") ext.loseContext();\n'
    '}\n'
    'try{\n'
    'if(Array.isArray(viewer.models))viewer.models.length=0;\n'
    'if(Array.isArray(viewer.labels))viewer.labels.length=0;\n'
    'if(Array.isArray(viewer.shapes))viewer.shapes.length=0;\n'
    'viewer.surfaces={};viewer.scene=null;viewer.modelGroup=null;viewer.rotationGroup=null;\n'
    '}catch(e){}\n'
    '} catch(e) {}\n'
    '};\n'
    # py3Dmol-backed viewers (Submit, ChemDarwin, TURBOMOLE) get a fresh
    # container and a fresh window.viewer_<id> global on every render, so there
    # is no previous handle for the caller to release. Their WebGL contexts,
    # observers and window-level listeners otherwise survive the div that owned
    # them. Browsers cap live contexts and evict the oldest, which blacks out
    # viewers the user is still working in elsewhere.
    'window.__delfinDisposeOrphanedViewers = function(){\n'
    'try {\n'
    'for (var k in window) {\n'
    'if (!Object.prototype.hasOwnProperty.call(window, k)) continue;\n'
    'if (k.indexOf("viewer_") !== 0) continue;\n'
    'var v = window[k];\n'
    'if (!v || typeof v !== "object" || v.__delfinDisposed) continue;\n'
    'var el = window.__delfinResolveViewerElement(v, null);\n'
    'if (!el) continue;\n'
    'if (document.body && document.body.contains(el)) continue;\n'
    'window.__delfinDisposeViewer(v);\n'
    'try { window[k] = null; } catch(e) {}\n'
    '}\n'
    '} catch(e) {}\n'
    '};\n'
    'window.__delfinPatchAllKnown3DmolViewers = function(){\n'
    'try {\n'
    'for (var k in window) {\n'
    'if (!Object.prototype.hasOwnProperty.call(window, k)) continue;\n'
    'if (k.indexOf("viewer_") !== 0) continue;\n'
    'var v = window[k];\n'
    'if (v && typeof v === "object") {\n'
    'var suffix = k.substring(7);\n'
    'var el = document.getElementById("3dmolviewer_" + suffix);\n'
    'window.__delfinEnableRightDragTranslate(v, el);\n'
    '}\n'
    '}\n'
    'var scopeMaps = [window._calcMolViewerByScope,window._calcTrajViewerByScope,'
    'window._remoteMolViewerByScope,window._remoteTrajViewerByScope];\n'
    'for (var i = 0; i < scopeMaps.length; i++) {\n'
    'var map = scopeMaps[i];\n'
    'if(!map || typeof map !== "object") continue;\n'
    'for (var key in map) {\n'
    'if(!Object.prototype.hasOwnProperty.call(map, key)) continue;\n'
    'var sv = map[key];\n'
    'if(sv && typeof sv === "object") window.__delfinEnableRightDragTranslate(sv, null);\n'
    '}\n'
    '}\n'
    'if(window._orcaBuildViewer) '
    'window.__delfinEnableRightDragTranslate(window._orcaBuildViewer,null);\n'
    '} catch(e) {}\n'
    '};\n'
    'window.__delfinDownloadDataUrl = function(dataUrl, filename){\n'
    'try {\n'
    'if(!dataUrl) return false;\n'
    'var link = document.createElement("a");\n'
    'link.href = dataUrl;\n'
    'link.download = filename || "viewer.png";\n'
    'document.body.appendChild(link);\n'
    'link.click();\n'
    'document.body.removeChild(link);\n'
    'return true;\n'
    '} catch(e) {\n'
    'return false;\n'
    '}\n'
    '};\n'
    'window.__delfinCanvasToPngDataUrl = function(canvas){\n'
    'try {\n'
    'if(!canvas || typeof canvas.toDataURL !== "function") return null;\n'
    'return canvas.toDataURL("image/png");\n'
    '} catch(e) {\n'
    'return null;\n'
    '}\n'
    '};\n'
    'window.__delfinCloneCanvasDataUrl = function(sourceCanvas, scale){\n'
    'try {\n'
    'if(!sourceCanvas) return null;\n'
    'var srcW = sourceCanvas.width || sourceCanvas.clientWidth || 0;\n'
    'var srcH = sourceCanvas.height || sourceCanvas.clientHeight || 0;\n'
    'if(!srcW || !srcH) return null;\n'
    'var factor = Math.max(1, parseInt(scale || 1, 10) || 1);\n'
    'var exportCanvas = document.createElement("canvas");\n'
    'exportCanvas.width = srcW * factor;\n'
    'exportCanvas.height = srcH * factor;\n'
    'var ctx = exportCanvas.getContext("2d");\n'
    'if(!ctx) return null;\n'
    'ctx.imageSmoothingEnabled = true;\n'
    'ctx.imageSmoothingQuality = "high";\n'
    'ctx.drawImage(sourceCanvas, 0, 0, exportCanvas.width, exportCanvas.height);\n'
    'return window.__delfinCanvasToPngDataUrl(exportCanvas);\n'
    '} catch(e) {\n'
    'return null;\n'
    '}\n'
    '};\n'
    'window.__delfinDownloadViewerPng = function(viewer, options){\n'
    'try {\n'
    'if(!viewer) return false;\n'
    'var opts = options || {};\n'
    'if(typeof viewer.render === "function") viewer.render();\n'
    'var el = window.__delfinResolveViewerElement(viewer, opts.element || null);\n'
    'var canvas = el ? el.querySelector("canvas") : null;\n'
    'var dataUrl = window.__delfinCloneCanvasDataUrl(canvas, opts.scale || 1)\n'
    ' || window.__delfinCanvasToPngDataUrl(canvas);\n'
    'return window.__delfinDownloadDataUrl(dataUrl, opts.filename || "viewer.png");\n'
    '} catch(e) {\n'
    'return false;\n'
    '}\n'
    '};\n'
    'var attachPromiseHook = function(){\n'
    'var p = window.$3Dmolpromise;\n'
    'if(p && typeof p.then === "function" && p !== window.__delfinMousePatchPromise){\n'
    'window.__delfinMousePatchPromise = p;\n'
    'p.then(function(){\n'
    'window.__delfinPatchAllKnown3DmolViewers();\n'
    '}).catch(function(){});\n'
    '}\n'
    '};\n'
    'var bootstrapAttempt=0;\n'
    'function bootstrap(){\n'
    'window.__delfinPatchAllKnown3DmolViewers();\n'
    'attachPromiseHook();\n'
    'bootstrapAttempt++;\n'
    'if(typeof window.$3Dmol==="undefined"&&bootstrapAttempt<20){\n'
    'window.setTimeout(bootstrap,250);\n'
    '}\n'
    '}\n'
    'bootstrap();\n'
    '})();'
)


def patch_viewer_mouse_controls_js(viewer_var='viewer', viewer_element_var='null'):
    """Return JS that installs right-drag translate and applies it to a viewer."""
    return (
        RIGHT_MOUSE_TRANSLATE_PATCH_JS
        + '\n'
        + (
            'if (window.__delfinEnableRightDragTranslate) '
            f'{{ window.__delfinEnableRightDragTranslate({viewer_var}, {viewer_element_var}); }}'
        )
    )


MEASUREMENT_BOOTSTRAP_JS = r"""
(function() {
    if (window.__delfinMeasureReady) return;
    window.__delfinMeasureReady = true;
    window._delfinMeasurePicks = window._delfinMeasurePicks || {};
    window._delfinMeasureActive = window._delfinMeasureActive || {};
    window._delfinMeasureShapes = window._delfinMeasureShapes || {};

    var VIEWER_MAPS = [
        '_remoteTrajViewerByScope', '_remoteMolViewerByScope',
        '_calcTrajViewerByScope', '_calcMolViewerByScope'
    ];
    var COLORS = ['#ffd54f', '#4fc3f7', '#81c784', '#f06292'];

    function getViewer(scopeKey) {
        for (var i = 0; i < VIEWER_MAPS.length; i++) {
            var m = window[VIEWER_MAPS[i]];
            if (m && m[scopeKey]) return m[scopeKey];
        }
        return null;
    }
    function sub(a, b){ return {x:a.x-b.x, y:a.y-b.y, z:a.z-b.z}; }
    function dot(a, b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
    function norm(a){ return Math.sqrt(dot(a,a)); }
    function cross(a, b){
        return {x:a.y*b.z-a.z*b.y, y:a.z*b.x-a.x*b.z, z:a.x*b.y-a.y*b.x};
    }
    function scale(a, s){ return {x:a.x*s, y:a.y*s, z:a.z*s}; }
    function dist(a, b){ var d = sub(a,b); return Math.sqrt(dot(d,d)); }
    function angleDeg(a, b, c){
        var u = sub(a,b), v = sub(c,b);
        var nu = norm(u), nv = norm(v);
        if (nu < 1e-9 || nv < 1e-9) return 0;
        var c1 = Math.max(-1, Math.min(1, dot(u,v)/(nu*nv)));
        return Math.acos(c1) * 180 / Math.PI;
    }
    function dihedralDeg(a, b, c, d){
        var b1 = sub(b,a), b2 = sub(c,b), b3 = sub(d,c);
        var nb2 = norm(b2);
        if (nb2 < 1e-9) return 0;
        var b2n = scale(b2, 1/nb2);
        var n1 = cross(b1, b2), n2 = cross(b2, b3);
        var m1 = cross(n1, b2n);
        var x = dot(n1, n2), y = dot(m1, n2);
        return Math.atan2(y, x) * 180 / Math.PI;
    }
    function getAtomBySerial(viewer, serial) {
        try {
            var atoms = viewer.selectedAtoms({serial: serial});
            if (atoms && atoms.length) return atoms[0];
        } catch (e) {}
        try {
            var model = viewer.getModel();
            if (model) {
                var atoms2 = model.selectedAtoms({serial: serial});
                if (atoms2 && atoms2.length) return atoms2[0];
            }
        } catch (e) {}
        return null;
    }
    function findDisplay(scopeKey) {
        var root = document.querySelector('.' + scopeKey);
        if (!root) return null;
        return root.querySelector('.delfin-xyz-measure-display');
    }
    function updateDisplay(scopeKey) {
        var el = findDisplay(scopeKey);
        if (!el) return;
        var picks = window._delfinMeasurePicks[scopeKey] || [];
        if (!picks.length) {
            el.innerHTML = '<span style="color:#888;">— click atoms (2/3/4) —</span>';
            return;
        }
        var labels = picks.map(function(p){
            var e = (p.atom.elem || p.atom.atom || '?');
            return e + (p.atom.serial != null ? p.atom.serial : '?');
        });
        var coords = picks.map(function(p){ return {x:p.atom.x, y:p.atom.y, z:p.atom.z}; });
        var lines = [];
        lines.push(
            '<div style="color:#555;font-size:0.9em;margin-bottom:2px;">[' +
            labels.join(' → ') + ']</div>'
        );
        if (picks.length === 1) {
            lines.push('<div style="color:#888;">pick more atoms…</div>');
        }
        for (var i = 0; i + 1 < picks.length; i++) {
            var d = dist(coords[i], coords[i+1]).toFixed(3);
            lines.push(
                '<div><span style="color:#1976d2;font-weight:600;">d(' +
                labels[i] + ',' + labels[i+1] + ')</span> = ' + d + ' Å</div>'
            );
        }
        for (var j = 0; j + 2 < picks.length; j++) {
            var a = angleDeg(coords[j], coords[j+1], coords[j+2]).toFixed(2);
            lines.push(
                '<div><span style="color:#2e7d32;font-weight:600;">∠(' +
                labels[j] + ',' + labels[j+1] + ',' + labels[j+2] +
                ')</span> = ' + a + '°</div>'
            );
        }
        if (picks.length >= 4) {
            var t = dihedralDeg(coords[0], coords[1], coords[2], coords[3]).toFixed(2);
            lines.push(
                '<div><span style="color:#c62828;font-weight:600;">τ(' +
                labels[0] + ',' + labels[1] + ',' + labels[2] + ',' + labels[3] +
                ')</span> = ' + t + '°</div>'
            );
        }
        el.innerHTML = lines.join('');
    }
    function redraw(viewer, scopeKey) {
        var prev = window._delfinMeasureShapes[scopeKey] || [];
        prev.forEach(function(s){ try { viewer.removeShape(s); } catch (_){} });
        var picks = window._delfinMeasurePicks[scopeKey] || [];
        var shapes = [];
        picks.forEach(function(p, i) {
            var fresh = getAtomBySerial(viewer, p.atom.serial);
            if (fresh) p.atom = fresh;
            var sh = viewer.addSphere({
                center: {x: p.atom.x, y: p.atom.y, z: p.atom.z},
                radius: 0.78,
                color: COLORS[i % COLORS.length],
                opacity: 0.45
            });
            shapes.push(sh);
        });
        window._delfinMeasureShapes[scopeKey] = shapes;
        try { viewer.render(); } catch (_){}
        updateDisplay(scopeKey);
    }
    function attach(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return false;
        try {
            viewer.setClickable({}, true, function(atom) {
                if (!window._delfinMeasureActive[scopeKey]) return;
                if (!atom || atom.serial === undefined || atom.serial === null) return;
                var picks = window._delfinMeasurePicks[scopeKey] = window._delfinMeasurePicks[scopeKey] || [];
                var found = -1;
                for (var i = 0; i < picks.length; i++) {
                    if (picks[i].atom.serial === atom.serial) { found = i; break; }
                }
                if (found >= 0) {
                    picks.splice(found, 1);
                } else {
                    if (picks.length >= 4) picks.shift();
                    picks.push({atom: atom});
                }
                redraw(viewer, scopeKey);
            });
        } catch (e) { return false; }
        return true;
    }
    function detach(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        try { viewer.setClickable({}, false, function(){}); } catch (e) {}
    }
    function clearPicks(scopeKey) {
        var viewer = getViewer(scopeKey);
        var prev = window._delfinMeasureShapes[scopeKey] || [];
        if (viewer) prev.forEach(function(s){ try { viewer.removeShape(s); } catch (_){} });
        window._delfinMeasureShapes[scopeKey] = [];
        window._delfinMeasurePicks[scopeKey] = [];
        if (viewer) { try { viewer.render(); } catch (_){} }
        updateDisplay(scopeKey);
    }
    function refresh(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        redraw(viewer, scopeKey);
    }
    function setActive(scopeKey, active) {
        window._delfinMeasureActive[scopeKey] = !!active;
        if (active) { attach(scopeKey); refresh(scopeKey); }
        else { detach(scopeKey); }
    }
    function ensureAfterRender(scopeKey) {
        if (window._delfinMeasureActive[scopeKey]) {
            attach(scopeKey);
            refresh(scopeKey);
        }
    }
    window.__delfinMeasure = {
        attach: attach,
        detach: detach,
        clear: clearPicks,
        refresh: refresh,
        setActive: setActive,
        ensureAfterRender: ensureAfterRender
    };
})();
"""


def measurement_bootstrap_js():
    """Return the one-time JS that installs window.__delfinMeasure helpers."""
    return MEASUREMENT_BOOTSTRAP_JS


SUBMIT_MANIP_BOOTSTRAP_JS = r"""
(function() {
    if (window.__delfinSubmitManipReady) return;
    window.__delfinSubmitManipReady = true;

    window._submitMolViewerByScope = window._submitMolViewerByScope || {};
    window._submitManipStateByScope = window._submitManipStateByScope || {};

    var COLORS = ['#ffd54f','#4fc3f7','#81c784','#f06292','#ba68c8','#ffb74d'];
    var PIVOT_COLOR = '#e53935';
    var UNDO_LIMIT = 50;
    var ROT_RAD_PER_PX = 0.01;
    var DRAG_THRESHOLD_PX = 3;
    var PICK_RADIUS_PX = 14;

    function getViewer(scopeKey) {
        return (window._submitMolViewerByScope || {})[scopeKey] || null;
    }
    function getState(scopeKey) {
        var s = window._submitManipStateByScope[scopeKey];
        if (!s) {
            s = {
                mode: 'off',
                scopeKey: scopeKey,
                ffActive: false,
                ffFrameMs: 16,
                picks: [],
                pivot: null,
                shapes: [],
                pivotShape: null,
                undo: [],
                overlay: null,
                viewerEl: null,
                canvas: null,
                rect: null,
                drag: null
            };
            window._submitManipStateByScope[scopeKey] = s;
        }
        return s;
    }
    function getAtoms(viewer) {
        if (!viewer) return [];
        try {
            var m = viewer.getModel();
            if (m && typeof m.selectedAtoms === 'function') {
                return m.selectedAtoms({}) || [];
            }
        } catch (e) {}
        return [];
    }
    function getRoot(scopeKey) {
        return document.querySelector('.' + scopeKey);
    }
    function getSyncInput(scopeKey) {
        var root = getRoot(scopeKey);
        if (!root) return null;
        var wrap = root.querySelector('.submit-manip-sync');
        if (!wrap) return null;
        return wrap.querySelector('input, textarea');
    }
    function getStatusEl(scopeKey) {
        var root = getRoot(scopeKey);
        if (!root) return null;
        return root.querySelector('.submit-manip-status');
    }

    function vecAdd(a,b){return {x:a.x+b.x,y:a.y+b.y,z:a.z+b.z};}
    function vecSub(a,b){return {x:a.x-b.x,y:a.y-b.y,z:a.z-b.z};}
    function vecScale(a,s){return {x:a.x*s,y:a.y*s,z:a.z*s};}
    function vecLen(a){return Math.sqrt(a.x*a.x+a.y*a.y+a.z*a.z);}
    function vecNorm(a){var l=vecLen(a); return l<1e-9?{x:0,y:0,z:0}:vecScale(a,1/l);}
    function crossV(a,b){return {x:a.y*b.z-a.z*b.y, y:a.z*b.x-a.x*b.z, z:a.x*b.y-a.y*b.x};}

    // --- Camera-space basis from rotationGroup matrix (3Dmol uses THREE.js) ---
    function getCameraBasis(viewer) {
        try {
            // Scene rotation matrix elements (column-major in three.js)
            var e = viewer.rotationGroup.matrix.elements;
            // Rows of rotation matrix = camera axes expressed in world coords
            var right = {x: e[0], y: e[4], z: e[8]};
            var up    = {x: e[1], y: e[5], z: e[9]};
            var fwd   = {x: e[2], y: e[6], z: e[10]};
            return {right: right, up: up, fwd: fwd};
        } catch (e) {
            return {right:{x:1,y:0,z:0}, up:{x:0,y:1,z:0}, fwd:{x:0,y:0,z:1}};
        }
    }
    function getPixelToWorld(viewer, canvas) {
        try {
            var camZ = Math.abs(viewer.rotationGroup.position.z || 150);
            var fov = (viewer.camera && viewer.camera.fov) ? viewer.camera.fov : 20;
            var fovRad = fov * Math.PI / 180;
            var h = canvas.clientHeight || canvas.height || 600;
            return 2 * camZ * Math.tan(fovRad / 2) / h;
        } catch (e) {
            return 0.03;
        }
    }
    function mat4Apply(m, x, y, z, w) {
        var e = m.elements;
        return {
            x: e[0]*x + e[4]*y + e[8]*z + e[12]*w,
            y: e[1]*x + e[5]*y + e[9]*z + e[13]*w,
            z: e[2]*x + e[6]*y + e[10]*z + e[14]*w,
            w: e[3]*x + e[7]*y + e[11]*z + e[15]*w
        };
    }
    function modelToScreen(viewer, canvas, p) {
        if (!viewer || !canvas) return null;
        var coord = {x: p.x, y: p.y, z: p.z};  // strip extra atom props
        // Try 3Dmol's built-in first
        if (typeof viewer.modelToScreen === 'function') {
            try {
                var r = viewer.modelToScreen(coord);
                if (Array.isArray(r)) r = r[0];
                if (r && isFinite(r.x) && isFinite(r.y)) return r;
            } catch (e) {}
        }
        // Manual projection via THREE.js-style Vector3.project if available
        try {
            var THREE_NS = window.$3Dmol || window.THREE;
            if (THREE_NS && THREE_NS.Vector3 && viewer.camera) {
                if (viewer.modelGroup) {
                    viewer.modelGroup.updateMatrixWorld && viewer.modelGroup.updateMatrixWorld(true);
                }
                var v = new THREE_NS.Vector3(coord.x, coord.y, coord.z);
                // Apply modelGroup transform (scene rotation+pan)
                if (viewer.modelGroup && typeof v.applyMatrix4 === 'function') {
                    v.applyMatrix4(viewer.modelGroup.matrixWorld);
                }
                // Project to NDC using the camera
                if (typeof v.project === 'function') {
                    v.project(viewer.camera);
                    var w = canvas.clientWidth || canvas.width || 600;
                    var h = canvas.clientHeight || canvas.height || 600;
                    var sx = (v.x + 1) * w / 2;
                    var sy = (-v.y + 1) * h / 2;
                    if (isFinite(sx) && isFinite(sy)) return {x: sx, y: sy};
                }
            }
        } catch (e) {}
        // Last-resort manual projection using matrices (force render for freshness)
        try { if (typeof viewer.render === 'function') viewer.render(); } catch (e) {}
        try {
            var mgMat = viewer.modelGroup && viewer.modelGroup.matrixWorld;
            var camInv = viewer.camera && viewer.camera.matrixWorldInverse;
            var proj = viewer.camera && viewer.camera.projectionMatrix;
            if (!mgMat || !camInv || !proj) return null;
            var world = mat4Apply(mgMat, coord.x, coord.y, coord.z, 1);
            var cam = mat4Apply(camInv, world.x, world.y, world.z, world.w);
            var clip = mat4Apply(proj, cam.x, cam.y, cam.z, cam.w);
            if (Math.abs(clip.w) < 1e-9) return null;
            var ndcX = clip.x / clip.w;
            var ndcY = clip.y / clip.w;
            var w2 = canvas.clientWidth || canvas.width || 600;
            var h2 = canvas.clientHeight || canvas.height || 600;
            return {x: (ndcX + 1) * w2 / 2, y: (-ndcY + 1) * h2 / 2};
        } catch (e) { return null; }
    }

    function getAtomBySerial(viewer, serial) {
        try {
            var atoms = getAtoms(viewer);
            for (var i = 0; i < atoms.length; i++) {
                if (atoms[i].serial === serial) return atoms[i];
            }
        } catch (e) {}
        return null;
    }

    function serializeXyz(viewer, header) {
        var atoms = getAtoms(viewer);
        if (!atoms.length) return '';
        var lines = [atoms.length.toString(), header || 'Edited in DELFIN viewer'];
        for (var i = 0; i < atoms.length; i++) {
            var a = atoms[i];
            var el = a.elem || a.atom || 'X';
            lines.push(el + ' ' +
                a.x.toFixed(6) + ' ' +
                a.y.toFixed(6) + ' ' +
                a.z.toFixed(6)
            );
        }
        return lines.join('\n');
    }

    function snapshotForUndo(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        var state = getState(scopeKey);
        var atoms = getAtoms(viewer);
        var snap = atoms.map(function(a) {
            return {serial: a.serial, x: a.x, y: a.y, z: a.z};
        });
        state.undo.push(snap);
        if (state.undo.length > UNDO_LIMIT) state.undo.shift();
    }

    function restoreFromSnapshot(scopeKey, snap) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        var atoms = getAtoms(viewer);
        var bySerial = {};
        for (var i = 0; i < atoms.length; i++) bySerial[atoms[i].serial] = atoms[i];
        for (var j = 0; j < snap.length; j++) {
            var a = bySerial[snap[j].serial];
            if (a) { a.x = snap[j].x; a.y = snap[j].y; a.z = snap[j].z; }
        }
        try { viewer.render(); } catch (e) {}
    }

    function pushXyzToPython(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        var input = getSyncInput(scopeKey);
        if (!input) return;
        var xyz = serializeXyz(viewer);
        var proto = (input.tagName === 'TEXTAREA')
            ? window.HTMLTextAreaElement.prototype
            : window.HTMLInputElement.prototype;
        var setter = Object.getOwnPropertyDescriptor(proto, 'value');
        if (setter && setter.set) setter.set.call(input, xyz);
        else input.value = xyz;
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
    }

    // Force 3Dmol to rebuild atom/bond geometry after we mutated atom.x/y/z
    // directly. Without this, viewer.render() re-uses the cached mesh and the
    // molecule appears frozen.
    function invalidateGeometry(viewer) {
        try {
            var m = viewer.getModel();
            if (m) {
                m.molObj = null;
                if (typeof m.rebuildBonds === 'function') {
                    try { m.rebuildBonds(); } catch (e) {}
                }
            }
        } catch (e) {}
    }

    function distV(a, b) {
        var dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }
    function angleV(a, b, c) {
        var ux = a.x - b.x, uy = a.y - b.y, uz = a.z - b.z;
        var vx = c.x - b.x, vy = c.y - b.y, vz = c.z - b.z;
        var nu = Math.sqrt(ux*ux + uy*uy + uz*uz);
        var nv = Math.sqrt(vx*vx + vy*vy + vz*vz);
        if (nu < 1e-9 || nv < 1e-9) return 0;
        var cs = (ux*vx + uy*vy + uz*vz) / (nu * nv);
        cs = Math.max(-1, Math.min(1, cs));
        return Math.acos(cs) * 180 / Math.PI;
    }
    function dihedralV(a, b, c, d) {
        var b1x = b.x-a.x, b1y = b.y-a.y, b1z = b.z-a.z;
        var b2x = c.x-b.x, b2y = c.y-b.y, b2z = c.z-b.z;
        var b3x = d.x-c.x, b3y = d.y-c.y, b3z = d.z-c.z;
        var nb2 = Math.sqrt(b2x*b2x + b2y*b2y + b2z*b2z);
        if (nb2 < 1e-9) return 0;
        var b2nx = b2x/nb2, b2ny = b2y/nb2, b2nz = b2z/nb2;
        var n1x = b1y*b2z - b1z*b2y, n1y = b1z*b2x - b1x*b2z, n1z = b1x*b2y - b1y*b2x;
        var n2x = b2y*b3z - b2z*b3y, n2y = b2z*b3x - b2x*b3z, n2z = b2x*b3y - b2y*b3x;
        var m1x = n1y*b2nz - n1z*b2ny;
        var m1y = n1z*b2nx - n1x*b2nz;
        var m1z = n1x*b2ny - n1y*b2nx;
        var x = n1x*n2x + n1y*n2y + n1z*n2z;
        var y = m1x*n2x + m1y*n2y + m1z*n2z;
        return Math.atan2(y, x) * 180 / Math.PI;
    }

    function ensureMeasureBox(scopeKey) {
        var state = getState(scopeKey);
        if (state.measureBox) return state.measureBox;
        if (!state.viewerEl) return null;
        var box = document.createElement('div');
        box.className = 'submit-manip-measure-box';
        box.style.position = 'absolute';
        box.style.top = '8px';
        box.style.right = '8px';
        box.style.maxWidth = '260px';
        box.style.maxHeight = 'calc(100% - 16px)';
        box.style.padding = '6px 9px';
        box.style.background = 'rgba(255,255,255,0.92)';
        box.style.border = '1px solid #cfd8dc';
        box.style.borderRadius = '4px';
        box.style.boxShadow = '0 1px 3px rgba(0,0,0,0.15)';
        box.style.font = '12px/1.35 monospace';
        box.style.color = '#37474f';
        box.style.pointerEvents = 'auto';
        box.style.overflowY = 'auto';
        box.style.zIndex = '25';
        box.style.display = 'none';
        state.viewerEl.appendChild(box);
        state.measureBox = box;
        return box;
    }
    function updateMeasureBox(scopeKey) {
        var state = getState(scopeKey);
        var viewer = getViewer(scopeKey);
        var box = ensureMeasureBox(scopeKey);
        if (!box || !viewer) return;
        var picks = state.picks || [];
        if (!picks.length) {
            box.style.display = 'none';
            box.innerHTML = '';
            return;
        }
        var byS = {};
        var all = getAtoms(viewer);
        for (var i = 0; i < all.length; i++) byS[all[i].serial] = all[i];
        var pts = [], labels = [];
        for (var j = 0; j < picks.length; j++) {
            var a = byS[picks[j].serial];
            if (!a) continue;
            pts.push({x: a.x, y: a.y, z: a.z});
            labels.push((a.elem || '?') + a.serial);
        }
        if (!pts.length) { box.style.display = 'none'; return; }
        if (pts.length > 4) { box.style.display = 'none'; box.innerHTML = ''; return; }
        var parts = ['<div style="color:#555;margin-bottom:2px;">[' + labels.join(' → ') + ']</div>'];
        for (var k = 0; k + 1 < pts.length; k++) {
            parts.push(
                '<div><span style="color:#1976d2;font-weight:600;">d(' +
                labels[k] + ',' + labels[k+1] + ')</span> = ' +
                distV(pts[k], pts[k+1]).toFixed(3) + ' Å</div>'
            );
        }
        for (var m = 0; m + 2 < pts.length; m++) {
            parts.push(
                '<div><span style="color:#2e7d32;font-weight:600;">∠(' +
                labels[m] + ',' + labels[m+1] + ',' + labels[m+2] + ')</span> = ' +
                angleV(pts[m], pts[m+1], pts[m+2]).toFixed(2) + '°</div>'
            );
        }
        for (var q = 0; q + 3 < pts.length; q++) {
            parts.push(
                '<div><span style="color:#c62828;font-weight:600;">τ(' +
                labels[q] + ',' + labels[q+1] + ',' + labels[q+2] + ',' + labels[q+3] + ')</span> = ' +
                dihedralV(pts[q], pts[q+1], pts[q+2], pts[q+3]).toFixed(2) + '°</div>'
            );
        }
        box.innerHTML = parts.join('');
        box.style.display = 'block';
    }

    function redrawHighlights(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        var state = getState(scopeKey);
        invalidateGeometry(viewer);
        // Remove previous shapes
        state.shapes.forEach(function(s) {
            try { viewer.removeShape(s); } catch (e) {}
        });
        state.shapes = [];
        if (state.pivotShape) {
            try { viewer.removeShape(state.pivotShape); } catch (e) {}
            state.pivotShape = null;
        }
        // Pick highlights
        state.picks.forEach(function(p, i) {
            var atom = getAtomBySerial(viewer, p.serial);
            if (!atom) return;
            var sh = viewer.addSphere({
                center: {x: atom.x, y: atom.y, z: atom.z},
                radius: 0.72,
                color: COLORS[i % COLORS.length],
                opacity: 0.45
            });
            state.shapes.push(sh);
        });
        // Pivot ring
        if (state.pivot) {
            var pa = getAtomBySerial(viewer, state.pivot.serial);
            if (pa) {
                state.pivotShape = viewer.addSphere({
                    center: {x: pa.x, y: pa.y, z: pa.z},
                    radius: 0.95,
                    color: PIVOT_COLOR,
                    opacity: 0.28
                });
            }
        }
        try { viewer.render(); } catch (e) {}
        updateStatus(scopeKey);
        updateMeasureBox(scopeKey);
    }

    function updateStatus(scopeKey) {
        var el = getStatusEl(scopeKey);
        if (!el) return;
        var state = getState(scopeKey);
        var n = state.picks.length;
        var pivotTxt = '';
        if (state.pivot) {
            pivotTxt = ' · pivot: <b>' +
                (state.pivot.elem || '?') + state.pivot.serial + '</b>';
        }
        var undoTxt = state.undo.length
            ? ' · <span style="color:#888;">' + state.undo.length + ' undo</span>'
            : '';
        var modeTxt = state.mode === 'select' ? 'SELECT'
                    : state.mode === 'manipulate' ? 'MANIPULATE' : '';
        var modeBadge = modeTxt
            ? '<span style="color:#1976d2;font-weight:600;">' + modeTxt + '</span> · '
            : '';
        var hint = '';
        if (state.mode === 'select' && n === 0) {
            hint = ' <span style="color:#888;font-size:0.9em;">(click atom to pick · hold <b>Shift</b>+drag = rect)</span>';
        } else if (state.mode === 'manipulate' && n === 0) {
            hint = ' <span style="color:#888;font-size:0.9em;">(pick atoms first)</span>';
        }
        el.innerHTML = modeBadge +
            '<b>' + n + '</b> atom' + (n === 1 ? '' : 's') + ' selected' +
            pivotTxt + undoTxt + hint;
    }

    // --- Pick toggle ---
    function togglePick(scopeKey, atom, additive) {
        var state = getState(scopeKey);
        if (!atom || atom.serial === undefined) return;
        var found = -1;
        for (var i = 0; i < state.picks.length; i++) {
            if (state.picks[i].serial === atom.serial) { found = i; break; }
        }
        if (!additive) {
            // plain click: toggle
            if (found >= 0) state.picks.splice(found, 1);
            else state.picks.push({serial: atom.serial, elem: atom.elem || 'X'});
        } else {
            // additive: only add if not present
            if (found < 0) state.picks.push({serial: atom.serial, elem: atom.elem || 'X'});
        }
        redrawHighlights(scopeKey);
    }

    // --- Atom picking ---
    // 3Dmol's own picker and our projection-based one used to run as two
    // independent click handlers that each called togglePick, so a click on an
    // atom selected and immediately deselected it — the most basic operation in
    // the editor did nothing whenever 3Dmol's picker worked. There is now a
    // single resolution path: 3Dmol only *records* what it hit, and one
    // deferred handler decides and toggles exactly once.
    function attachClickable(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        try {
            viewer.setClickable({}, true, function(atom) {
                var state = getState(scopeKey);
                if (state.mode !== 'select') return;
                if (!atom || atom.serial === undefined) return;
                state._nativeHit = {serial: atom.serial, elem: atom.elem || 'X'};
            });
            // Force 3Dmol to rebuild any internal pick buffers.
            try { viewer.render(); } catch (e) {}
        } catch (e) {}
        installCanvasClickFallback(scopeKey);
    }
    function detachClickable(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        try { viewer.setClickable({}, false, function(){}); } catch (e) {}
        uninstallCanvasClickFallback(scopeKey);
    }

    // Resolve the atom under a viewport position. Depth-aware, no synthetic
    // events: the old implementation dispatched a mousedown/mouseup/click
    // triple at the canvas and swapped 3Dmol's click callback around it just to
    // read back one atom.
    function probeClickAtom(scopeKey, clientX, clientY) {
        var atom = raycastAtom(scopeKey, clientX, clientY);
        if (!atom) return null;
        return {serial: atom.serial, elem: atom.elem || 'X'};
    }

    // The single click-resolution path. Runs deferred so 3Dmol's own picker has
    // already recorded its hit (it is depth-exact where it works); we fall back
    // to our projection when it reports nothing.
    function installCanvasClickFallback(scopeKey) {
        var state = getState(scopeKey);
        var canvas = state.canvas;
        if (!canvas) return;
        if (state._canvasClickHandler) {
            try { canvas.removeEventListener('click', state._canvasClickHandler, true); } catch (e) {}
        }
        var handler = function(e) {
            if (state.mode !== 'select') return;
            if (state.shiftHeld) return; // shift+drag path handles the rubber band
            if (state.drag && state.drag.movedEnough) return;
            var additive = !!(e.shiftKey || e.ctrlKey || e.metaKey);
            var cx = e.clientX, cy = e.clientY;
            state._nativeHit = null;
            window.setTimeout(function() {
                var hit = state._nativeHit;
                state._nativeHit = null;
                var atom = hit ? getAtomBySerial(getViewer(scopeKey), hit.serial) : null;
                if (!atom) atom = raycastAtom(scopeKey, cx, cy);
                if (atom) togglePick(scopeKey, atom, additive);
            }, 0);
        };
        canvas.addEventListener('click', handler, true);
        state._canvasClickHandler = handler;
    }
    function uninstallCanvasClickFallback(scopeKey) {
        var state = getState(scopeKey);
        if (state._canvasClickHandler && state.canvas) {
            try { state.canvas.removeEventListener('click', state._canvasClickHandler, true); } catch (e) {}
        }
        state._canvasClickHandler = null;
    }

    // --- Overlay + mouse handling ---
    function ensureOverlay(scopeKey) {
        var state = getState(scopeKey);
        if (state.overlay) return state.overlay;
        var el = state.viewerEl;
        if (!el) return null;
        if (window.getComputedStyle(el).position === 'static') {
            el.style.position = 'relative';
        }
        var ov = document.createElement('div');
        ov.className = 'submit-manip-overlay';
        ov.style.position = 'absolute';
        ov.style.left = '0';
        ov.style.top = '0';
        ov.style.right = '0';
        ov.style.bottom = '0';
        ov.style.pointerEvents = 'none';
        ov.style.zIndex = '20';
        el.appendChild(ov);
        state.overlay = ov;
        bindOverlayEvents(scopeKey);
        return ov;
    }

    // Project one model point to canvas pixels *and* normalised device depth.
    // Depth is what lets picking prefer the atom in front, which plain
    // screen-distance picking cannot do: in a fused ring or a crowded metal
    // coordination sphere the nearest atom on screen is regularly the one
    // hidden behind the structure.
    function projectWithDepth(viewer, canvas, p) {
        try {
            if (viewer.modelGroup && viewer.modelGroup.updateMatrixWorld) {
                viewer.modelGroup.updateMatrixWorld(true);
            }
            if (viewer.camera && viewer.camera.updateMatrixWorld) {
                viewer.camera.updateMatrixWorld(true);
            }
            var mgMat = viewer.modelGroup && viewer.modelGroup.matrixWorld;
            var camInv = viewer.camera && viewer.camera.matrixWorldInverse;
            var proj = viewer.camera && viewer.camera.projectionMatrix;
            if (!mgMat || !camInv || !proj) return null;
            var world = mat4Apply(mgMat, p.x, p.y, p.z, 1);
            var cam = mat4Apply(camInv, world.x, world.y, world.z, world.w);
            var clip = mat4Apply(proj, cam.x, cam.y, cam.z, cam.w);
            if (!isFinite(clip.w) || Math.abs(clip.w) < 1e-9) return null;
            var w = canvas.clientWidth || canvas.width || 600;
            var h = canvas.clientHeight || canvas.height || 600;
            var sx = (clip.x / clip.w + 1) * w / 2;
            var sy = (-clip.y / clip.w + 1) * h / 2;
            if (!isFinite(sx) || !isFinite(sy)) return null;
            return {x: sx, y: sy, depth: clip.z / clip.w};
        } catch (e) { return null; }
    }

    // Screen positions of every atom, in canvas pixels, with depth. One pass,
    // reused by picking and by rubber-band selection.
    function projectAllAtoms(scopeKey) {
        var viewer = getViewer(scopeKey);
        var state = getState(scopeKey);
        if (!viewer || !state.canvas) return [];
        var atoms = getAtoms(viewer);
        var out = [];
        for (var i = 0; i < atoms.length; i++) {
            var p = projectWithDepth(viewer, state.canvas, atoms[i]);
            if (!p) continue;
            out.push({atom: atoms[i], x: p.x, y: p.y, depth: p.depth});
        }
        return out;
    }

    function raycastAtom(scopeKey, clientX, clientY) {
        var state = getState(scopeKey);
        if (!state.canvas) return null;
        var rect = state.canvas.getBoundingClientRect();
        var sx = clientX - rect.left;
        var sy = clientY - rect.top;
        var projected = projectAllAtoms(scopeKey);
        var best = null, bestDepth = Infinity, bestDist = Infinity;
        var radius2 = PICK_RADIUS_PX * PICK_RADIUS_PX;
        for (var i = 0; i < projected.length; i++) {
            var q = projected[i];
            var dx = q.x - sx, dy = q.y - sy;
            var d2 = dx*dx + dy*dy;
            if (d2 > radius2) continue;
            // Nearest to the camera wins; screen distance only breaks ties
            // between atoms at effectively the same depth.
            if (q.depth < bestDepth - 1e-6 ||
                (Math.abs(q.depth - bestDepth) <= 1e-6 && d2 < bestDist)) {
                best = q.atom; bestDepth = q.depth; bestDist = d2;
            }
        }
        return best;
    }

    function beginRectDraw(scopeKey, x0, y0) {
        var state = getState(scopeKey);
        if (state.rect) try { state.overlay.removeChild(state.rect); } catch(e) {}
        var r = document.createElement('div');
        r.style.position = 'absolute';
        r.style.left = x0 + 'px';
        r.style.top = y0 + 'px';
        r.style.width = '0px';
        r.style.height = '0px';
        r.style.border = '1px dashed #1976d2';
        r.style.background = 'rgba(25,118,210,0.08)';
        r.style.pointerEvents = 'none';
        state.overlay.appendChild(r);
        state.rect = r;
    }
    function updateRect(scopeKey, x0, y0, x1, y1) {
        var state = getState(scopeKey);
        if (!state.rect) return;
        var left = Math.min(x0, x1), top = Math.min(y0, y1);
        state.rect.style.left = left + 'px';
        state.rect.style.top = top + 'px';
        state.rect.style.width = Math.abs(x1 - x0) + 'px';
        state.rect.style.height = Math.abs(y1 - y0) + 'px';
    }
    // Rubber-band selection: project every atom once and keep those inside the
    // band. The previous implementation dispatched a synthetic mousedown /
    // mouseup / click triple at every 7 px of the rectangle — thousands of
    // events for a band over a medium molecule, which froze the tab and, if
    // Shift came up before the mouse button, re-entered the toggle path and
    // scrambled the selection.
    // ``additive`` unions with the current picks (Shift/Ctrl held); otherwise
    // the band replaces the selection, which is what every other editor does.
    function finishRect(scopeKey, x0, y0, x1, y1, additive) {
        var state = getState(scopeKey);
        if (state.rect) { try { state.overlay.removeChild(state.rect); } catch(e) {} state.rect = null; }
        if (!state.canvas) return;

        var rect = state.canvas.getBoundingClientRect();
        var minX = Math.min(x0, x1) - rect.left, maxX = Math.max(x0, x1) - rect.left;
        var minY = Math.min(y0, y1) - rect.top,  maxY = Math.max(y0, y1) - rect.top;
        if (maxX - minX < 3 || maxY - minY < 3) return;

        if (!additive) state.picks = [];

        var projected = projectAllAtoms(scopeKey);
        for (var i = 0; i < projected.length; i++) {
            var q = projected[i];
            if (q.x < minX || q.x > maxX || q.y < minY || q.y > maxY) continue;
            // Behind the camera / outside the frustum in depth.
            if (q.depth < -1 || q.depth > 1) continue;
            var serial = q.atom.serial;
            var exists = state.picks.some(function(p) { return p.serial === serial; });
            if (!exists) state.picks.push({serial: serial, elem: q.atom.elem || 'X'});
        }
        updateStatus(scopeKey);
        redrawHighlights(scopeKey);
    }

    // --- internal coordinates: set a bond, angle or dihedral --------------
    // The routine edits a chemist actually asks for: "make this bond 2.10 A",
    // "open that angle to 120 deg", "turn this substituent to 60 deg". Which
    // half of the molecule moves is decided from the model's own connectivity:
    // cut the bond the coordinate turns about and move the fragment that is no
    // longer attached to the anchor. Inside a ring both halves stay connected,
    // so only the terminal atom moves and the caller is told why.
    function indexOfSerial(atoms, serial) {
        for (var i = 0; i < atoms.length; i++) {
            if (atoms[i].serial === serial) return i;
        }
        return -1;
    }
    function bondAdjacency(viewer) {
        var atoms = getAtoms(viewer);
        var adj = new Array(atoms.length);
        for (var i = 0; i < atoms.length; i++) {
            adj[i] = (atoms[i].bonds || []).slice();
        }
        return adj;
    }
    function fragmentFrom(adj, start, cutA, cutB) {
        var seen = {}, stack = [start], out = [];
        seen[start] = true;
        while (stack.length) {
            var a = stack.pop();
            out.push(a);
            var nb = adj[a] || [];
            for (var t = 0; t < nb.length; t++) {
                var b = nb[t];
                if ((a === cutA && b === cutB) || (a === cutB && b === cutA)) continue;
                if (seen[b]) continue;
                seen[b] = true;
                stack.push(b);
            }
        }
        return {atoms: out, seen: seen};
    }
    function pickedIndices(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return null;
        var state = getState(scopeKey);
        var atoms = getAtoms(viewer);
        var idx = [];
        for (var i = 0; i < state.picks.length; i++) {
            var at = indexOfSerial(atoms, state.picks[i].serial);
            if (at < 0) return null;
            idx.push(at);
        }
        return idx;
    }
    // The value the current selection describes: two atoms are a bond, three
    // an angle, four a dihedral, in the order they were picked.
    function readInternal(scopeKey) {
        var viewer = getViewer(scopeKey);
        var idx = pickedIndices(scopeKey);
        if (!viewer || !idx) return null;
        var a = getAtoms(viewer);
        if (idx.length === 2) {
            return {kind: 'bond', unit: 'A', idx: idx,
                    value: distV(a[idx[0]], a[idx[1]])};
        }
        if (idx.length === 3) {
            return {kind: 'angle', unit: 'deg', idx: idx,
                    value: angleV(a[idx[0]], a[idx[1]], a[idx[2]])};
        }
        if (idx.length === 4) {
            return {kind: 'dihedral', unit: 'deg', idx: idx,
                    value: dihedralV(a[idx[0]], a[idx[1]], a[idx[2]], a[idx[3]])};
        }
        return null;
    }
    function translateAtoms(atoms, indices, delta) {
        for (var i = 0; i < indices.length; i++) {
            var at = atoms[indices[i]];
            at.x += delta.x; at.y += delta.y; at.z += delta.z;
        }
    }
    function rotateAtomsAbout(atoms, indices, origin, axis, angle) {
        for (var i = 0; i < indices.length; i++) {
            var at = atoms[indices[i]];
            var rel = rotateAboutAxis(vecSub(at, origin), axis, angle);
            var np = vecAdd(rel, origin);
            at.x = np.x; at.y = np.y; at.z = np.z;
        }
    }
    function setInternal(scopeKey, target) {
        var viewer = getViewer(scopeKey);
        var info = readInternal(scopeKey);
        if (!viewer || !info) {
            return {ok: false, error: 'pick 2, 3 or 4 atoms first'};
        }
        if (typeof target !== 'number' || !isFinite(target)) {
            return {ok: false, error: 'not a number'};
        }
        var atoms = getAtoms(viewer);
        var adj = bondAdjacency(viewer);
        var idx = info.idx;
        var note = '';
        snapshotForUndo(scopeKey);

        if (info.kind === 'bond') {
            if (target <= 0) return {ok: false, error: 'a bond must be positive'};
            var i = idx[0], j = idx[1];
            var u = vecNorm(vecSub(atoms[j], atoms[i]));
            if (vecLen(u) < 0.5) return {ok: false, error: 'atoms coincide'};
            var frag = fragmentFrom(adj, j, i, j);
            var moving = frag.atoms;
            if (frag.seen[i]) { moving = [j]; note = 'ring: only the second atom moved'; }
            translateAtoms(atoms, moving, vecScale(u, target - info.value));
        } else if (info.kind === 'angle') {
            var i2 = idx[0], j2 = idx[1], k2 = idx[2];
            var axis = vecNorm(crossV(vecSub(atoms[i2], atoms[j2]),
                                     vecSub(atoms[k2], atoms[j2])));
            if (vecLen(axis) < 0.5) {
                return {ok: false, error: 'the three atoms are collinear'};
            }
            var frag2 = fragmentFrom(adj, k2, j2, k2);
            var moving2 = frag2.atoms;
            if (frag2.seen[j2]) { moving2 = [k2]; note = 'ring: only the third atom moved'; }
            var d2 = (target - info.value) * Math.PI / 180;
            rotateAtomsAbout(atoms, moving2, atoms[j2], axis, d2);
            if (Math.abs(angleV(atoms[i2], atoms[j2], atoms[k2]) - target) > 1e-3) {
                rotateAtomsAbout(atoms, moving2, atoms[j2], axis, -2 * d2);
            }
        } else {
            var j3 = idx[1], k3 = idx[2];
            var axis3 = vecNorm(vecSub(atoms[k3], atoms[j3]));
            if (vecLen(axis3) < 0.5) {
                return {ok: false, error: 'the central bond has no length'};
            }
            var frag3 = fragmentFrom(adj, k3, j3, k3);
            var moving3 = frag3.atoms;
            if (frag3.seen[j3]) {
                return {ok: false,
                        error: 'that dihedral turns about a ring bond'};
            }
            var d3 = (target - info.value) * Math.PI / 180;
            rotateAtomsAbout(atoms, moving3, atoms[j3], axis3, d3);
            var got = dihedralV(atoms[idx[0]], atoms[j3], atoms[k3], atoms[idx[3]]);
            if (Math.abs(((got - target + 540) % 360) - 180) > 1e-3) {
                rotateAtomsAbout(atoms, moving3, atoms[j3], axis3, -2 * d3);
            }
        }

        invalidateGeometry(viewer);
        try { viewer.render(); } catch (e) {}
        redrawHighlights(scopeKey);
        pushXyzToPython(scopeKey);
        var after = readInternal(scopeKey);
        return {ok: true, kind: info.kind, was: info.value,
                now: after ? after.value : target, note: note};
    }

    // --- live force field -------------------------------------------------
    // Avogadro relaxes the molecule while you drag: the grabbed atom follows
    // the cursor exactly and everything else settles around it. The relaxation
    // runs in the browser (window.__delfinFF) because a per-frame round trip to
    // the kernel costs 45 ms and collapses to 13 Hz under a drag. Term indices
    // are 0-based into the model's atom order, which is the XYZ order Python
    // exported the parameters from.
    function ffIndicesOf(viewer, serials) {
        var atoms = getAtoms(viewer);
        var out = [];
        for (var i = 0; i < atoms.length; i++) {
            if (serials.indexOf(atoms[i].serial) >= 0) out.push(i);
        }
        return out;
    }
    function ffReadPositions(viewer) {
        var atoms = getAtoms(viewer);
        var p = new Float64Array(3 * atoms.length);
        for (var i = 0; i < atoms.length; i++) {
            p[3*i] = atoms[i].x; p[3*i+1] = atoms[i].y; p[3*i+2] = atoms[i].z;
        }
        return p;
    }
    function ffWritePositions(viewer, pos) {
        var atoms = getAtoms(viewer);
        if (!pos || pos.length < 3 * atoms.length) return false;
        for (var i = 0; i < atoms.length; i++) {
            atoms[i].x = pos[3*i]; atoms[i].y = pos[3*i+1]; atoms[i].z = pos[3*i+2];
        }
        return true;
    }
    function ffEnabled(state) {
        return !!(state.ffActive && window.__delfinFF &&
                  window._delfinFFByScope && window._delfinFFByScope[state.scopeKey]);
    }
    function ffBeginDrag(scopeKey, targets) {
        var state = getState(scopeKey);
        if (!ffEnabled(state)) return;
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        try {
            window.__delfinFF.grab(scopeKey, ffIndicesOf(viewer, targets));
            state.ffFrameMs = 16;
        } catch (e) {}
    }
    function ffRelaxFrame(scopeKey) {
        var state = getState(scopeKey);
        if (!ffEnabled(state)) return false;
        var viewer = getViewer(scopeKey);
        if (!viewer) return false;
        var t0 = nowMs();
        try {
            var out = window.__delfinFF.step(
                scopeKey, ffReadPositions(viewer), state.ffFrameMs || 16);
            if (!out) return false;
            ffWritePositions(viewer, out);
        } catch (e) { return false; }
        // The budget the engine adapts to is a *full frame*: 3Dmol's own
        // geometry rebuild costs up to 12 ms at 400 atoms and comes out of the
        // same 33 ms, so it has to be measured together with the relaxation.
        state.ffFrameMs = nowMs() - t0;
        return true;
    }
    // Avogadro's Auto Optimization tool: the force field keeps running on a
    // timer whether or not the mouse is down, so the structure visibly settles
    // and a grabbed atom drags the relaxed molecule along with it. Avogadro 1
    // ticks at 50 ms; an animation frame is the browser's equivalent and lets
    // the engine's own 33 ms batch budget do the pacing.
    function autoOptimizeRunning(scopeKey) {
        var state = getState(scopeKey);
        return !!state.autoOpt;
    }
    function autoOptimizeTick(scopeKey) {
        var state = getState(scopeKey);
        if (!state.autoOpt) return;
        var viewer = getViewer(scopeKey);
        if (!viewer || !ffEnabled(state)) { stopAutoOptimize(scopeKey); return; }
        // While a drag is in flight its own handler already relaxes and
        // redraws; running a second batch here would double the step budget.
        if (!state.drag) {
            if (ffRelaxFrame(scopeKey)) {
                redrawHighlights(scopeKey);
                var now = nowMs();
                // The coordinate box follows at a readable rate, not per frame:
                // each push is a widget round trip.
                if (now - (state.autoPushed || 0) > 500) {
                    state.autoPushed = now;
                    pushXyzToPython(scopeKey);
                }
            }
        }
        state.autoRaf = window.requestAnimationFrame(function() {
            autoOptimizeTick(scopeKey);
        });
    }
    function startAutoOptimize(scopeKey) {
        var state = getState(scopeKey);
        if (state.autoOpt) return true;
        if (!ffEnabled(state)) return false;
        snapshotForUndo(scopeKey);
        state.autoOpt = true;
        state.autoPushed = nowMs();
        autoOptimizeTick(scopeKey);
        updateStatus(scopeKey);
        return true;
    }
    function stopAutoOptimize(scopeKey) {
        var state = getState(scopeKey);
        if (!state.autoOpt) return false;
        state.autoOpt = false;
        if (state.autoRaf) {
            try { window.cancelAnimationFrame(state.autoRaf); } catch (e) {}
            state.autoRaf = null;
        }
        // The last frames since the throttled push have to reach Python too.
        pushXyzToPython(scopeKey);
        updateStatus(scopeKey);
        return true;
    }

    function ffEndDrag(scopeKey) {
        var state = getState(scopeKey);
        if (!ffEnabled(state)) return;
        try { window.__delfinFF.release(scopeKey); } catch (e) {}
    }
    function nowMs() {
        return (window.performance && typeof window.performance.now === 'function')
            ? window.performance.now() : Date.now();
    }

    // ``serials`` is the set of atoms this drag owns — the current selection
    // when the user grabbed a selected atom, a single atom when they grabbed an
    // unselected one. Falls back to the selection so callers without a drag
    // context keep working.
    function applyTranslate(scopeKey, deltaWorld, serials) {
        var viewer = getViewer(scopeKey);
        var state = getState(scopeKey);
        if (!viewer) return;
        var targets = serials || state.picks.map(function(p) { return p.serial; });
        var atoms = getAtoms(viewer);
        var byS = {};
        for (var i = 0; i < atoms.length; i++) byS[atoms[i].serial] = atoms[i];
        targets.forEach(function(serial) {
            var a = byS[serial];
            if (!a) return;
            a.x += deltaWorld.x;
            a.y += deltaWorld.y;
            a.z += deltaWorld.z;
        });
        redrawHighlights(scopeKey);
    }

    function rotateAboutAxis(v, axis, angle) {
        // Rodrigues
        var c = Math.cos(angle), s = Math.sin(angle);
        var k = vecNorm(axis);
        var dot = k.x*v.x + k.y*v.y + k.z*v.z;
        var cross = {
            x: k.y*v.z - k.z*v.y,
            y: k.z*v.x - k.x*v.z,
            z: k.x*v.y - k.y*v.x
        };
        return {
            x: v.x*c + cross.x*s + k.x*dot*(1-c),
            y: v.y*c + cross.y*s + k.y*dot*(1-c),
            z: v.z*c + cross.z*s + k.z*dot*(1-c)
        };
    }

    function applyRotate(scopeKey, dxPx, dyPx) {
        var viewer = getViewer(scopeKey);
        var state = getState(scopeKey);
        if (!viewer || !state.pivot) {
            return;
        }
        var pivotAtom = getAtomBySerial(viewer, state.pivot.serial);
        if (!pivotAtom) {
            return;
        }
        var pivot = {x: pivotAtom.x, y: pivotAtom.y, z: pivotAtom.z};
        var basis = getCameraBasis(viewer);
        var yaw = dxPx * ROT_RAD_PER_PX;
        var pitch = dyPx * ROT_RAD_PER_PX;
        var atoms = getAtoms(viewer);
        var byS = {};
        for (var i = 0; i < atoms.length; i++) byS[atoms[i].serial] = atoms[i];
        var moved = 0;
        state.picks.forEach(function(p) {
            if (p.serial === state.pivot.serial) return;
            var a = byS[p.serial];
            if (!a) return;
            var rel = vecSub({x:a.x,y:a.y,z:a.z}, pivot);
            rel = rotateAboutAxis(rel, basis.up, yaw);
            rel = rotateAboutAxis(rel, basis.right, pitch);
            var np = vecAdd(rel, pivot);
            a.x = np.x; a.y = np.y; a.z = np.z;
            moved++;
        });
        redrawHighlights(scopeKey);
    }

    function bindOverlayEvents(scopeKey) {
        var state = getState(scopeKey);
        var ov = state.overlay;
        if (!ov) return;

        ov.addEventListener('contextmenu', function(e) {
            if (state.mode === 'manipulate') { e.preventDefault(); e.stopPropagation(); }
        });

        ov.addEventListener('mousedown', function(e) {
            if (state.mode === 'off') return;
            var rect = ov.getBoundingClientRect();
            var x = e.clientX - rect.left, y = e.clientY - rect.top;
            var atom = raycastAtom(scopeKey, e.clientX, e.clientY);

            if (state.mode === 'manipulate') {
                if (e.button === 2) {
                    // Right mouse always tries to pick a pivot atom under the cursor.
                    e.preventDefault(); e.stopPropagation();
                    var picked = probeClickAtom(scopeKey, e.clientX, e.clientY);
                    if (picked) {
                        state.pivot = {serial: picked.serial, elem: picked.elem || 'X'};
                        redrawHighlights(scopeKey);
                    }
                    // Rotation only makes sense with picks; still allow pivot change without picks.
                    if (state.pivot && state.picks.length > 0) {
                        state.drag = {
                            kind: 'rotate',
                            startX: e.clientX, startY: e.clientY,
                            lastX: e.clientX, lastY: e.clientY,
                            movedEnough: false, snapshotted: false
                        };
                    }
                    return;
                }
                if (e.button !== 0) return;

                // Avogadro's manipulate gesture: press on an atom and drag it.
                // Grabbing an atom that is part of the selection moves the whole
                // selection; grabbing any other atom moves just that atom,
                // without disturbing the selection.
                if (atom) {
                    e.preventDefault(); e.stopPropagation();
                    var inSelection = state.picks.some(function(p) {
                        return p.serial === atom.serial;
                    });
                    state.drag = {
                        kind: 'translate',
                        targets: inSelection
                            ? state.picks.map(function(p) { return p.serial; })
                            : [atom.serial],
                        grabbed: atom.serial,
                        startX: e.clientX, startY: e.clientY,
                        lastX: e.clientX, lastY: e.clientY,
                        movedEnough: false, snapshotted: false
                    };
                    ffBeginDrag(scopeKey, state.drag.targets);
                    return;
                }

                // Empty space: hand the press to 3Dmol so the camera turns.
                // Manipulate mode used to swallow every left press, so the
                // molecule could not be reoriented without leaving the mode and
                // losing the pivot workflow. Bubbling cannot do this for us —
                // 3Dmol binds its mouse handlers on the canvas, and the overlay
                // is the canvas's sibling, not its ancestor.
                if (!(e.shiftKey || e.ctrlKey || e.metaKey)) {
                    var v = getViewer(scopeKey);
                    if (v && typeof v._handleMouseDown === 'function') {
                        try { v._handleMouseDown(e); } catch (_e) {}
                    }
                    return;
                }
                if (state.picks.length > 0) {
                    // Modifier on empty space keeps the old behaviour: drag the
                    // whole selection from anywhere in the viewport.
                    e.preventDefault(); e.stopPropagation();
                    state.drag = {
                        kind: 'translate',
                        targets: state.picks.map(function(p) { return p.serial; }),
                        startX: e.clientX, startY: e.clientY,
                        lastX: e.clientX, lastY: e.clientY,
                        movedEnough: false, snapshotted: false
                    };
                    ffBeginDrag(scopeKey, state.drag.targets);
                }
                return;
            }

            // Select mode — the overlay is only interactive while Shift is
            // held, so reaching here means the user asked for a rubber band.
            if (e.button === 0) {
                e.preventDefault();
                e.stopPropagation();
                state.drag = {
                    kind: 'maybe-rect',
                    startX: e.clientX, startY: e.clientY,
                    origX: x, origY: y,
                    additive: !!(e.ctrlKey || e.metaKey),
                    movedEnough: false,
                    atomRef: atom
                };
            }
        });

        if (state._globalBound) return;
        state._globalBound = true;

        window.addEventListener('mousemove', function(e) {
            var state2 = window._submitManipStateByScope[scopeKey];
            if (!state2 || !state2.drag) return;
            var d = state2.drag;
            var dx = e.clientX - d.lastX, dy = e.clientY - d.lastY;
            var totX = e.clientX - d.startX, totY = e.clientY - d.startY;
            var totMag2 = totX*totX + totY*totY;
            if (!d.movedEnough && totMag2 > DRAG_THRESHOLD_PX*DRAG_THRESHOLD_PX) {
                d.movedEnough = true;
            }
            d.lastX = e.clientX; d.lastY = e.clientY;

            var viewer = getViewer(scopeKey);
            if (!viewer) return;

            if (d.kind === 'translate' && d.movedEnough) {
                e.preventDefault();
                if (!d.snapshotted) { snapshotForUndo(scopeKey); d.snapshotted = true; }
                var basis = getCameraBasis(viewer);
                var s = getPixelToWorld(viewer, state2.canvas);
                var delta = {
                    x: basis.right.x * dx * s - basis.up.x * dy * s,
                    y: basis.right.y * dx * s - basis.up.y * dy * s,
                    z: basis.right.z * dx * s - basis.up.z * dy * s
                };
                applyTranslate(scopeKey, delta, d.targets);
                // The grabbed atoms are already where the cursor put them;
                // the relaxation pulls everything else after them.
                if (ffRelaxFrame(scopeKey)) redrawHighlights(scopeKey);
            } else if (d.kind === 'rotate' && d.movedEnough) {
                e.preventDefault();
                if (!d.snapshotted) { snapshotForUndo(scopeKey); d.snapshotted = true; }
                applyRotate(scopeKey, dx, dy);
            } else if (d.kind === 'maybe-rect' && d.movedEnough) {
                // Lazily begin rect
                if (!state2.rect) beginRectDraw(scopeKey, d.origX, d.origY);
                var rect = state2.overlay.getBoundingClientRect();
                updateRect(scopeKey, d.origX, d.origY,
                    e.clientX - rect.left, e.clientY - rect.top);
            }
        }, true);

        window.addEventListener('mouseup', function(e) {
            var state2 = window._submitManipStateByScope[scopeKey];
            if (!state2 || !state2.drag) return;
            var d = state2.drag;
            state2.drag = null;
            if (d.kind === 'translate' || d.kind === 'rotate') {
                ffEndDrag(scopeKey);
                if (d.movedEnough) {
                    pushXyzToPython(scopeKey);
                }
            } else if (d.kind === 'maybe-rect') {
                if (d.movedEnough) {
                    // Use client (viewport) coords for hit-test — robust to
                    // any overlay/canvas offset mismatch.
                    finishRect(scopeKey,
                        d.startX, d.startY,
                        e.clientX, e.clientY,
                        d.additive);
                }
                // click-no-drag branches intentionally removed: in shift+drag
                // mode, a bare click without movement is a no-op. Non-shift
                // atom clicks are handled by 3Dmol's setClickable path.
            }
        }, true);
    }

    function setOverlayInteractive(scopeKey) {
        var state = getState(scopeKey);
        if (!state.overlay) return;
        if (state.mode === 'off') {
            state.overlay.style.pointerEvents = 'none';
            state.overlay.style.cursor = '';
        } else if (state.mode === 'select') {
            // In select mode the overlay is passthrough by default so 3Dmol
            // can handle atom clicks directly (via setClickable). Holding
            // Shift turns the overlay on for rubber-band selection.
            var active = !!state.shiftHeld;
            state.overlay.style.pointerEvents = active ? 'auto' : 'none';
            state.overlay.style.cursor = active ? 'crosshair' : '';
        } else if (state.mode === 'manipulate') {
            // Always interactive: any atom can be grabbed directly, and a press
            // on empty space is forwarded so the camera still turns.
            state.overlay.style.pointerEvents = 'auto';
            state.overlay.style.cursor = 'move';
        }
    }

    // --- Public API ---
    function onViewerReady(scopeKey, viewerEl) {
        var state = getState(scopeKey);
        // External re-render → previous picks refer to stale atoms: reset.
        state.picks = [];
        state.pivot = null;
        state.shapes = [];
        state.pivotShape = null;
        state.undo = [];
        state.viewerEl = viewerEl || state.viewerEl ||
            (getRoot(scopeKey) ? getRoot(scopeKey).querySelector('.submit-mol-output') : null);
        if (state.viewerEl) {
            state.canvas = state.viewerEl.querySelector('canvas');
        }
        // Overlay is attached fresh per render (old one is gone with the HTML)
        state.overlay = null;
        state.rect = null;
        state.drag = null;
        // The parameters were assigned for the geometry that just went away.
        if (state.autoOpt && state.autoRaf) {
            try { window.cancelAnimationFrame(state.autoRaf); } catch (e) {}
        }
        state.autoOpt = false;
        state.autoRaf = null;
        state.ffActive = false;
        state.ffInfo = null;
        state.measureBox = null;
        ensureOverlay(scopeKey);
        setOverlayInteractive(scopeKey);
        redrawHighlights(scopeKey);
        if (state.mode === 'select') attachClickable(scopeKey);
    }

    function setMode(scopeKey, mode) {
        var state = getState(scopeKey);
        state.mode = (mode === 'select' || mode === 'manipulate') ? mode : 'off';
        ensureOverlay(scopeKey);
        setOverlayInteractive(scopeKey);
        updateStatus(scopeKey);
        if (state.mode === 'select') {
            attachClickable(scopeKey);
        } else {
            detachClickable(scopeKey);
        }
    }


    function clearPicks(scopeKey) {
        var state = getState(scopeKey);
        state.picks = [];
        state.pivot = null;
        redrawHighlights(scopeKey);
    }

    function undo(scopeKey) {
        var state = getState(scopeKey);
        if (!state.undo.length) return;
        var snap = state.undo.pop();
        restoreFromSnapshot(scopeKey, snap);
        redrawHighlights(scopeKey);
        pushXyzToPython(scopeKey);
    }

    // Window-level mousedown intercept: the apply_molecule_view_style patch
    // binds right-drag-pan on the viewer element in capture phase, which
    // would swallow our right-click before the overlay sees it. Binding on
    // window catches the event earlier and lets us stop the patch.
    if (!window.__delfinSubmitManipWindowBound) {
        window.__delfinSubmitManipWindowBound = true;
        window.addEventListener('mousedown', function(e) {
            var states = window._submitManipStateByScope || {};
            for (var k in states) {
                var s = states[k];
                if (!s || !s.viewerEl) continue;
                if (!s.viewerEl.contains(e.target) && e.target !== s.viewerEl) continue;
                if (s.mode !== 'manipulate') continue;
                if (e.button !== 2) continue;  // only steal right-button
                e.preventDefault(); e.stopImmediatePropagation();
                // Forward to our normal mousedown logic by synthesising
                // a direct call (our overlay listener expects this shape).
                try {
                    var picked = probeClickAtom(k, e.clientX, e.clientY);
                    if (picked) {
                        s.pivot = {serial: picked.serial, elem: picked.elem || 'X'};
                        redrawHighlights(k);
                    }
                    if (s.pivot && s.picks.length > 0) {
                        s.drag = {
                            kind: 'rotate',
                            startX: e.clientX, startY: e.clientY,
                            lastX: e.clientX, lastY: e.clientY,
                            movedEnough: false, snapshotted: false
                        };
                    }
                } catch (_) {}
                return;
            }
        }, true);
        // Suppress context menu anywhere inside a manipulate-active viewer.
        window.addEventListener('contextmenu', function(e) {
            var states = window._submitManipStateByScope || {};
            for (var k in states) {
                var s = states[k];
                if (!s || !s.viewerEl) continue;
                if (!s.viewerEl.contains(e.target) && e.target !== s.viewerEl) continue;
                if (s.mode === 'manipulate') {
                    e.preventDefault(); e.stopImmediatePropagation();
                    return;
                }
            }
        }, true);
    }

    // Keyboard shortcuts + Shift tracking for rubber-band
    if (!window.__delfinSubmitManipKeyBound) {
        window.__delfinSubmitManipKeyBound = true;

        function propagateShift(isDown) {
            var states = window._submitManipStateByScope || {};
            Object.keys(states).forEach(function(k) {
                var s = states[k];
                if (!s) return;
                s.shiftHeld = isDown;
                if (s.mode === 'select') setOverlayInteractive(k);
            });
        }
        window.addEventListener('keydown', function(e) {
            if (e.key === 'Shift') { propagateShift(true); }
            var key = e.key || '';
            if ((e.ctrlKey || e.metaKey) && (key === 'z' || key === 'Z') && !e.shiftKey) {
                var states = window._submitManipStateByScope || {};
                var keys = Object.keys(states);
                for (var i = 0; i < keys.length; i++) {
                    var s = states[keys[i]];
                    if (s && (s.mode === 'select' || s.mode === 'manipulate') && s.undo.length) {
                        e.preventDefault();
                        undo(keys[i]);
                        break;
                    }
                }
            }
        }, true);
        window.addEventListener('keyup', function(e) {
            if (e.key === 'Shift') { propagateShift(false); }
        }, true);
        window.addEventListener('blur', function() { propagateShift(false); }, true);
    }

    // Fullscreen: on toggle, move viewer + toolbar + isomer nav + copy row into
    // a floating overlay; on exit, put them back where they were. No layout
    // changes to the default DOM, so canvas alignment stays intact.
    if (!window.__delfinSubmitFullscreenBound) {
        window.__delfinSubmitFullscreenBound = true;
        window._submitFsByScope = window._submitFsByScope || {};

        function findScope(el) {
            while (el && el.classList) {
                for (var i = 0; i < el.classList.length; i++) {
                    if (el.classList[i].indexOf('submit-scope-') === 0) {
                        return el.classList[i];
                    }
                }
                el = el.parentElement;
            }
            return null;
        }
        function resizeScopeViewer(scopeKey) {
            try {
                var viewer = (window._submitMolViewerByScope || {})[scopeKey];
                if (!viewer) return;
                [60, 250].forEach(function(delay) {
                    setTimeout(function() {
                        try {
                            if (typeof viewer.resize === 'function') viewer.resize();
                            if (typeof viewer.render === 'function') viewer.render();
                        } catch (e) {}
                    }, delay);
                });
            } catch (e) {}
        }
        function setFsIcon(btn, active) {
            if (!btn) return;
            var icon = btn.querySelector('i.fa');
            if (!icon) return;
            icon.classList.remove('fa-expand');
            icon.classList.remove('fa-compress');
            icon.classList.add(active ? 'fa-compress' : 'fa-expand');
            btn.setAttribute('title', active ? 'Exit fullscreen (Esc)' : 'Toggle fullscreen (Esc to exit)');
        }
        function enterFullscreen(scopeKey) {
            var root = document.querySelector('.' + scopeKey);
            if (!root) return;
            var selectors = [
                '.submit-fs-member-toolbar',
                '.submit-fs-member-viewer',
                '.submit-fs-member-isomer',
                '.submit-fs-member-copyrow'
            ];
            var members = [];
            for (var i = 0; i < selectors.length; i++) {
                var el = root.querySelector(selectors[i]);
                if (el) members.push(el);
            }
            if (!members.length) return;
            var overlay = document.createElement('div');
            overlay.className = 'submit-fs-overlay ' + scopeKey;
            var restore = members.map(function(el) {
                return { el: el, parent: el.parentNode, next: el.nextSibling };
            });
            members.forEach(function(el) { overlay.appendChild(el); });
            document.body.appendChild(overlay);
            window._submitFsByScope[scopeKey] = { overlay: overlay, restore: restore };
            var btn = overlay.querySelector('.submit-fullscreen-btn');
            setFsIcon(btn, true);
            resizeScopeViewer(scopeKey);
        }
        function exitFullscreen(scopeKey) {
            var entry = window._submitFsByScope[scopeKey];
            if (!entry) return;
            // Restore in reverse so each element's recorded nextSibling is
            // already back in the original parent before we insertBefore.
            for (var i = entry.restore.length - 1; i >= 0; i--) {
                var r = entry.restore[i];
                try {
                    if (r.next && r.next.parentNode === r.parent) {
                        r.parent.insertBefore(r.el, r.next);
                    } else {
                        r.parent.appendChild(r.el);
                    }
                } catch (e) {}
            }
            try { entry.overlay.parentNode.removeChild(entry.overlay); } catch (e) {}
            delete window._submitFsByScope[scopeKey];
            var root = document.querySelector('.' + scopeKey);
            var btn = root && root.querySelector('.submit-fullscreen-btn');
            setFsIcon(btn, false);
            resizeScopeViewer(scopeKey);
        }
        document.addEventListener('click', function(e) {
            var t = e.target;
            if (!t || !t.closest) return;
            var btn = t.closest('.submit-fullscreen-btn');
            if (!btn) return;
            var scopeKey = findScope(btn);
            if (!scopeKey) return;
            if (window._submitFsByScope[scopeKey]) {
                exitFullscreen(scopeKey);
            } else {
                enterFullscreen(scopeKey);
            }
        }, true);
        document.addEventListener('keydown', function(e) {
            if (e.key !== 'Escape') return;
            var keys = Object.keys(window._submitFsByScope || {});
            if (!keys.length) return;
            exitFullscreen(keys[0]);
        }, true);
    }

    // Hand the browser the force-field parameters Python assigned for the
    // current geometry, or null to switch live relaxation off again. Called
    // once when the mode is entered, never during a drag.
    function setForceField(scopeKey, terms) {
        var state = getState(scopeKey);
        if (!terms) {
            state.ffActive = false;
            state.ffInfo = null;
            if (window.__delfinFF) {
                try { window.__delfinFF.dispose(scopeKey); } catch (e) {}
            }
            updateStatus(scopeKey);
            return {ok: false, error: 'force field cleared'};
        }
        if (!window.__delfinFF) {
            return {ok: false, error: 'force-field engine not loaded'};
        }
        var result;
        try {
            result = window.__delfinFF.load(scopeKey, terms);
        } catch (e) {
            return {ok: false, error: 'force field failed to load'};
        }
        state.ffActive = !!(result && result.ok);
        state.ffInfo = result;
        updateStatus(scopeKey);
        return result;
    }

    window.__delfinSubmitManip = {
        onViewerReady: onViewerReady,
        setMode: setMode,
        clear: clearPicks,
        undo: undo,
        setForceField: setForceField,
        readInternal: readInternal,
        setInternal: setInternal,
        startAutoOptimize: startAutoOptimize,
        stopAutoOptimize: stopAutoOptimize,
        autoOptimizeRunning: autoOptimizeRunning
    };
})();
"""


def submit_manip_bootstrap_js():
    """Return one-time JS that installs window.__delfinSubmitManip helpers."""
    return SUBMIT_MANIP_BOOTSTRAP_JS


# Shared fullscreen support for the ORCA Builder, Calculations Browser, and
# Remote Archive molecule viewers.  The normal widget layout is deliberately
# left untouched: only the marked viewer/header/control members are moved into
# a body-level overlay, then restored to their exact original DOM positions.
STRUCTURE_VIEWER_FULLSCREEN_CSS = r"""
body.delfin-structure-fs-open {
    overflow: hidden !important;
}
.delfin-structure-fs-overlay,
.delfin-structure-fs-overlay * {
    box-sizing: border-box !important;
}
.delfin-structure-fs-overlay {
    position: fixed !important;
    inset: 0 !important;
    z-index: 2147483000 !important;
    display: flex !important;
    flex-direction: column !important;
    width: 100vw !important;
    height: 100vh !important;
    max-width: none !important;
    max-height: none !important;
    min-width: 0 !important;
    min-height: 0 !important;
    padding: 8px !important;
    margin: 0 !important;
    gap: 6px !important;
    overflow: hidden !important;
    background: #fff !important;
}
.delfin-structure-fs-overlay > .delfin-structure-fs-header,
.delfin-structure-fs-overlay > .delfin-structure-fs-toolbar {
    flex: 0 0 auto !important;
    width: 100% !important;
    max-width: none !important;
    min-width: 0 !important;
    margin: 0 !important;
    overflow: visible !important;
}
.delfin-structure-fs-overlay > .delfin-structure-fs-view-row {
    display: flex !important;
    flex: 1 1 0 !important;
    flex-flow: row nowrap !important;
    align-items: stretch !important;
    width: 100% !important;
    height: auto !important;
    max-width: none !important;
    min-width: 0 !important;
    min-height: 0 !important;
    margin: 0 !important;
    overflow: hidden !important;
}
.delfin-structure-fs-overlay > .delfin-structure-fs-viewer,
.delfin-structure-fs-overlay .delfin-structure-fs-view-wrap,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer {
    flex: 1 1 0 !important;
    width: 100% !important;
    height: 100% !important;
    max-width: none !important;
    max-height: none !important;
    min-width: 0 !important;
    min-height: 0 !important;
    margin: 0 !important;
    overflow: hidden !important;
}
.delfin-structure-fs-overlay .delfin-structure-fs-controls {
    flex: 0 0 clamp(300px, 28vw, 420px) !important;
    align-self: stretch !important;
    width: clamp(300px, 28vw, 420px) !important;
    min-width: 300px !important;
    max-width: 420px !important;
    max-height: 100% !important;
    overflow-y: auto !important;
    overflow-x: hidden !important;
}
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .output_area,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .output_subarea,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .output_wrapper,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .jp-OutputArea,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .jp-OutputArea-child,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .jp-OutputArea-output {
    width: 100% !important;
    height: 100% !important;
    max-width: none !important;
    max-height: none !important;
    min-width: 0 !important;
    min-height: 0 !important;
    padding: 0 !important;
    margin: 0 !important;
    overflow: hidden !important;
}
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .calc-mol-stage-wrapper,
.delfin-structure-fs-overlay .delfin-structure-fs-viewer .remote-mol-stage-wrapper {
    display: block !important;
    position: relative !important;
    width: 100% !important;
    height: 100% !important;
    max-width: none !important;
    max-height: none !important;
    min-width: 0 !important;
    min-height: 0 !important;
    margin: 0 !important;
    padding: 0 !important;
    overflow: hidden !important;
}
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="orca-mol-"],
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="orca-overlay-"],
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="mol3d_"],
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="calc_trj"],
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="remote_mol3d_"],
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="remote_trj_viewer_"],
.delfin-structure-fs-overlay .delfin-structure-fs-viewer [id^="3dmolviewer_"] {
    display: block !important;
    position: relative !important;
    width: 100% !important;
    height: 100% !important;
    max-width: none !important;
    max-height: none !important;
    min-width: 0 !important;
    min-height: 0 !important;
    margin: 0 !important;
    padding: 0 !important;
    overflow: hidden !important;
}
.delfin-structure-fs-overlay .delfin-structure-fs-viewer canvas {
    width: 100% !important;
    height: 100% !important;
    max-width: none !important;
    max-height: none !important;
}
@media (max-width: 800px) {
    .delfin-structure-fs-overlay > .delfin-structure-fs-view-row {
        flex-direction: column !important;
    }
    .delfin-structure-fs-overlay .delfin-structure-fs-controls {
        flex: 0 0 auto !important;
        align-self: stretch !important;
        width: 100% !important;
        min-width: 0 !important;
        max-width: none !important;
        max-height: 34vh !important;
    }
}
"""


STRUCTURE_VIEWER_FULLSCREEN_BOOTSTRAP_JS = r"""
(function() {
    if (window.__delfinStructureViewerFullscreenBound) return;
    window.__delfinStructureViewerFullscreenBound = true;
    window._delfinStructureFullscreen = window._delfinStructureFullscreen || null;

    function classWithPrefix(el, prefix) {
        while (el && el.classList) {
            for (var i = 0; i < el.classList.length; i++) {
                var cls = el.classList[i];
                if (cls.indexOf(prefix) === 0) return cls;
            }
            el = el.parentElement;
        }
        return null;
    }

    function moduleType(module) {
        if (module.classList.contains('orca-structure-fs-module')) return 'orca';
        if (module.classList.contains('calc-structure-fs-module')) return 'calc';
        if (module.classList.contains('remote-structure-fs-module')) return 'remote';
        return '';
    }

    function scopeFor(module, type) {
        if (type === 'orca') return classWithPrefix(module, 'orca-scope-');
        if (type === 'calc') return classWithPrefix(module, 'calc-scope-');
        if (type === 'remote') return classWithPrefix(module, 'remote-archive-scope-');
        return null;
    }

    function viewerFor(type, scopeKey) {
        if (type === 'orca') return window._orcaBuildViewer || null;
        if (type === 'calc') {
            return ((window._calcMolViewerByScope || {})[scopeKey]
                || (window._calcTrajViewerByScope || {})[scopeKey]
                || null);
        }
        if (type === 'remote') {
            return ((window._remoteMolViewerByScope || {})[scopeKey]
                || (window._remoteTrajViewerByScope || {})[scopeKey]
                || null);
        }
        return null;
    }

    function viewerView(viewer) {
        if (!viewer || typeof viewer.getView !== 'function') return null;
        try {
            var view = viewer.getView();
            return view && typeof view.slice === 'function' ? view.slice() : view;
        } catch (_e) {
            return null;
        }
    }

    function resizeViewer(entry) {
        if (!entry) return;
        [0, 80, 260].forEach(function(delay) {
            setTimeout(function() {
                var viewer = viewerFor(entry.type, entry.scopeKey);
                if (!viewer) return;
                try {
                    if (typeof viewer.resize === 'function') viewer.resize();
                    if (entry.view && typeof viewer.setView === 'function') {
                        viewer.setView(entry.view);
                    }
                    if (typeof viewer.render === 'function') viewer.render();
                } catch (_e) {}
            }, delay);
        });
    }

    function setButtonState(btn, active) {
        if (!btn) return;
        var icon = btn.querySelector('i.fa, i');
        if (icon) {
            icon.classList.remove('fa-expand');
            icon.classList.remove('fa-compress');
            icon.classList.add(active ? 'fa-compress' : 'fa-expand');
        }
        btn.setAttribute(
            'title',
            active ? 'Exit fullscreen (Esc)' : 'Toggle fullscreen (Esc to exit)'
        );
        btn.setAttribute('aria-pressed', active ? 'true' : 'false');
    }

    function exitFullscreen() {
        var entry = window._delfinStructureFullscreen;
        if (!entry) return;
        // Preserve the camera as it currently appears in fullscreen, including
        // any rotation, pan, or zoom the user applied while inspecting it.
        var currentView = viewerView(viewerFor(entry.type, entry.scopeKey));
        if (currentView) entry.view = currentView;
        for (var i = entry.restore.length - 1; i >= 0; i--) {
            var item = entry.restore[i];
            try {
                if (item.next && item.next.parentNode === item.parent) {
                    item.parent.insertBefore(item.el, item.next);
                } else {
                    item.parent.appendChild(item.el);
                }
            } catch (_e) {}
        }
        try {
            if (entry.overlay.parentNode) entry.overlay.parentNode.removeChild(entry.overlay);
        } catch (_e2) {}
        setButtonState(entry.btn, false);
        if (!entry.bodyHadOpenClass) {
            document.body.classList.remove('delfin-structure-fs-open');
        }
        window._delfinStructureFullscreen = null;
        resizeViewer(entry);
        [60, 280].forEach(function(delay) {
            setTimeout(function() {
                try { window.dispatchEvent(new Event('resize')); } catch (_e) {}
            }, delay);
        });
    }

    function enterFullscreen(btn) {
        var module = btn.closest('.delfin-structure-fs-module');
        if (!module) return;
        if (window._delfinStructureFullscreen) exitFullscreen();
        var type = moduleType(module);
        if (!type) return;
        var scopeKey = scopeFor(module, type);
        var candidates = module.querySelectorAll('.delfin-structure-fs-member');
        var members = [];
        for (var i = 0; i < candidates.length; i++) {
            var candidate = candidates[i];
            if (candidate.closest('.delfin-structure-fs-module') === module) {
                members.push(candidate);
            }
        }
        if (!members.length) return;
        var overlay = document.createElement('div');
        overlay.className = 'delfin-structure-fs-overlay delfin-structure-fs-' + type;
        if (scopeKey) overlay.classList.add(scopeKey);
        var restore = members.map(function(el) {
            return {el: el, parent: el.parentNode, next: el.nextSibling};
        });
        var entry = {
            btn: btn,
            module: module,
            overlay: overlay,
            restore: restore,
            type: type,
            scopeKey: scopeKey,
            view: viewerView(viewerFor(type, scopeKey)),
            bodyHadOpenClass: document.body.classList.contains('delfin-structure-fs-open')
        };
        members.forEach(function(el) { overlay.appendChild(el); });
        document.body.appendChild(overlay);
        document.body.classList.add('delfin-structure-fs-open');
        window._delfinStructureFullscreen = entry;
        setButtonState(btn, true);
        resizeViewer(entry);
        try { btn.focus(); } catch (_e) {}
    }

    document.addEventListener('click', function(e) {
        var target = e.target;
        if (!target || !target.closest) return;
        var btn = target.closest('.delfin-structure-fullscreen-btn');
        if (!btn) return;
        e.preventDefault();
        var entry = window._delfinStructureFullscreen;
        if (entry && entry.btn === btn) exitFullscreen();
        else enterFullscreen(btn);
    }, true);

    document.addEventListener('keydown', function(e) {
        if (e.key !== 'Escape' || !window._delfinStructureFullscreen) return;
        e.preventDefault();
        exitFullscreen();
    }, true);
})();
"""


def structure_viewer_fullscreen_css():
    """Return shared CSS for body-level 3D-viewer fullscreen overlays."""
    return STRUCTURE_VIEWER_FULLSCREEN_CSS


def structure_viewer_fullscreen_bootstrap_js():
    """Return one-time JS for ORCA/calc/remote fullscreen viewer buttons."""
    return STRUCTURE_VIEWER_FULLSCREEN_BOOTSTRAP_JS


def apply_molecule_view_style(view, zoom=DEFAULT_3DMOL_ZOOM, style=None):
    """Apply a shared ChemDarwin/MSILES-like style to a py3Dmol viewer.

    When ``style`` is None the function pulls the active representation and
    dimensions from the global viewer settings. Passing an explicit ``style``
    dict overrides only the molecule style; render-quality configuration still
    applies. Window zoom always honors the caller's ``zoom`` argument.
    """
    try:
        profile = get_viewer_profile()
    except Exception:
        profile = {
            'style': DEFAULT_3DMOL_STYLE,
            'viewer_config': {
                'backgroundColor': DEFAULT_3DMOL_BACKGROUND,
                'antialias': True,
                'disableFog': False,
            },
            'viewer_config_js': json.dumps(
                {
                    'backgroundColor': DEFAULT_3DMOL_BACKGROUND,
                    'antialias': True,
                    'disableFog': False,
                },
                separators=(',', ':'),
            ),
        }
    if style is None:
        style = profile['style']
    if hasattr(view, 'startjs'):
        # py3Dmol does not expose GLViewer's constructor config directly. The
        # view is not created until ``show()``, so replacing its stock config
        # here safely applies antialiasing/depth-shading to these viewers too.
        default_config_js = '{backgroundColor:"white"}'
        if default_config_js in view.startjs:
            view.startjs = view.startjs.replace(
                default_config_js,
                profile['viewer_config_js'],
                1,
            )
        marker = (
            'window.__delfinEnableRightDragTranslate('
            'viewer_UNIQUEID,document.getElementById("3dmolviewer_UNIQUEID"));'
        )
        if marker not in view.startjs:
            view.startjs += (
                '\n'
                + patch_viewer_mouse_controls_js(
                    'viewer_UNIQUEID',
                    'document.getElementById("3dmolviewer_UNIQUEID")',
                )
                + '\n'
                + (
                    '(function(){'
                    # Deferred recentering exists because the canvas is not
                    # always at its final size on the first frame. It must never
                    # run once the user has started interacting, or the pan,
                    # zoom or rotation they just made is silently thrown away.
                    'var __delfinRecenter=function(){'
                    'try{'
                    'if(viewer_UNIQUEID.__delfinUserInteracted) return;'
                    'viewer_UNIQUEID.zoomTo();'
                    'viewer_UNIQUEID.center();'
                    + (
                        f'viewer_UNIQUEID.zoom({zoom});'
                        if zoom is not None
                        else ''
                    )
                    + 'viewer_UNIQUEID.render();'
                    '}catch(_e){}'
                    '};'
                    'setTimeout(__delfinRecenter,0);'
                    'setTimeout(__delfinRecenter,120);'
                    'setTimeout(function(){'
                    'if(window.__delfinDisposeOrphanedViewers)'
                    'window.__delfinDisposeOrphanedViewers();'
                    '},0);'
                    '})();'
                )
                + '\n'
            )
    view.setStyle({}, style)
    view.setBackgroundColor(
        profile['viewer_config'].get('backgroundColor', DEFAULT_3DMOL_BACKGROUND)
    )
    view.zoomTo()
    view.center()
    if zoom is not None:
        view.zoom(zoom)
    view.render()
    return view


def render_xyz_in_output(output_widget, xyz_text, width=None, height=None):
    """Render an XYZ string inside an ipywidgets Output widget.

    Honors the global viewer settings (``ui.viewer.enabled`` and
    ``ui.viewer.quality``). Window size + zoom are constant across quality
    levels; the quality knob controls the representation style (line vs
    stick vs stick+sphere). Explicit ``width``/``height`` arguments still win.
    """
    profile = get_viewer_profile()
    with output_widget:
        clear_output()
        if not profile['enabled']:
            display(HTML(viewer_disabled_html()))
            return
        if not xyz_text or not xyz_text.strip():
            print('No coordinates to display.')
            return
        eff_width = width if width is not None else profile['width']
        eff_height = height if height is not None else profile['height']
        view = py3Dmol.view(width=eff_width, height=eff_height)
        view.addModel(xyz_text, 'xyz')
        apply_molecule_view_style(view, zoom=profile['zoom'], style=profile['style'])
        view.show()


def coord_to_xyz(coord_text):
    """Convert TURBOMOLE ``coord`` format to XYZ.

    Handles both ``$coord`` blocks (Bohr) and plain Cartesian lines.
    Returns an XYZ string with atom count + comment, or *None*.
    """
    lines = coord_text.strip().split('\n')
    xyz_lines = []
    bohr_to_ang = 0.529177249
    in_coord = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('$coord'):
            in_coord = True
            continue
        if stripped.startswith('$'):
            in_coord = False
            continue
        if in_coord:
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    x = float(parts[0]) * bohr_to_ang
                    y = float(parts[1]) * bohr_to_ang
                    z = float(parts[2]) * bohr_to_ang
                    xyz_lines.append(
                        f'{parts[3].capitalize()}  {x:12.6f}  {y:12.6f}  {z:12.6f}'
                    )
                except (ValueError, IndexError):
                    continue
    # Fallback: plain Cartesian lines (element x y z)
    if not xyz_lines:
        for line in lines:
            parts = line.split()
            if len(parts) >= 4 and parts[0][0].isalpha():
                try:
                    float(parts[1])
                    xyz_lines.append(line)
                except (ValueError, IndexError):
                    continue
    if not xyz_lines:
        return None
    return f'{len(xyz_lines)}\nTM Preview\n' + '\n'.join(xyz_lines)


def parse_xyz_frames(content):
    """Parse a multi-frame XYZ file into a list of ``(comment, xyz_block, n_atoms)`` tuples."""
    frames = []
    lines = content.strip().split('\n')
    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
        except (ValueError, IndexError):
            break
        if i + 1 >= len(lines):
            break
        comment = lines[i + 1].strip()
        coord_lines = []
        for j in range(n_atoms):
            if i + 2 + j < len(lines):
                coord_lines.append(lines[i + 2 + j])
        if len(coord_lines) == n_atoms:
            xyz_block = '\n'.join(coord_lines)
            frames.append((comment, xyz_block, n_atoms))
        i += 2 + n_atoms
    return frames


def strip_xyz_header(text):
    """Remove the XYZ header (atom count + comment line) if present."""
    text = text.strip()
    if not text:
        return text
    lines = text.split('\n')
    if len(lines) < 2:
        return text
    first_line = lines[0].strip()
    try:
        int(first_line)
        return '\n'.join(lines[2:]).strip()
    except ValueError:
        return text


# ---------------------------------------------------------------------------
# Fukui visualization helpers
# ---------------------------------------------------------------------------

_FUKUI_CUBE_OPTIONS = [
    ('none', None, False),
    ('ρ(N) — density(neutral)',  'density_neutral.cube', False),
    ('ρ(N+1) — density(anion)',  'density_anion.cube',   False),
    ('ρ(N−1) — density(cation)', 'density_cation.cube',  False),
    ('f⁺ — Fukui plus',          'fukui_plus.cube',      True),
    ('f⁻ — Fukui minus',         'fukui_minus.cube',     True),
    ('f⁰ — Fukui zero',          'fukui_zero.cube',      True),
]

_FUKUI_LABEL_OPTIONS = [
    ('none',          'none'),
    ('atom index',    'index'),
    ('q(N)',          'q_neutral'),
    ('q(N+1)',        'q_anion'),
    ('q(N−1)',        'q_cation'),
    ('f⁺',            'f_plus'),
    ('f⁻',            'f_minus'),
    ('f⁰',            'f_zero'),
]


def _fukui_atom_positions_from_xyz(xyz_text):
    """Return ``[(symbol, x, y, z)]`` from a (possibly headered) XYZ string."""
    text = strip_xyz_header(xyz_text or '')
    positions = []
    for raw in text.splitlines():
        parts = raw.split()
        if len(parts) < 4:
            continue
        try:
            x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
        except ValueError:
            continue
        positions.append((parts[0], x, y, z))
    return positions


def fukui_atom_labels_js(xyz_text, values, *, decimals=3, color='black'):
    """Build a JS snippet that adds floating per-atom labels to a 3Dmol viewer.

    Args:
        xyz_text: XYZ block (with or without header) used to place labels.
        values: Iterable of values matching the atom order in ``xyz_text``.
            Strings are emitted verbatim; numbers get ``decimals`` precision.
        decimals: Float precision for numeric values.
        color: Foreground color of the label text.

    Returns:
        JavaScript string. Callers paste it between ``addModel`` and ``zoomTo``.
    """
    if not values:
        return ''
    positions = _fukui_atom_positions_from_xyz(xyz_text)
    pieces = []
    for (sym, x, y, z), value in zip(positions, values):
        if isinstance(value, (int, float)):
            text = f"{value:+.{decimals}f}"
        else:
            text = str(value)
        text = text.replace("'", "\\'")
        pieces.append(
            "viewer.addLabel('" + text + "', {"
            "position:{x:" + f"{x:.4f}" + ",y:" + f"{y:.4f}" + ",z:" + f"{z:.4f}" + "},"
            "backgroundColor:'white',backgroundOpacity:0.7,"
            "fontColor:'" + color + "',fontSize:13,"
            "borderThickness:0.5,borderColor:'#666',"
            "showBackground:true,inFront:true});"
        )
    return "\n".join(pieces)


def fukui_cube_isosurface_js(cube_text, *, isoval=0.02, signed=False):
    """Build a JS snippet that paints a cube as an isosurface in a 3Dmol viewer.

    Args:
        cube_text: Raw Gaussian/ORCA cube file contents.
        isoval: Iso-value (positive). For signed cubes, ``+isoval`` (blue) and
            ``-isoval`` (red) are drawn together; for unsigned, only the
            positive lobe is drawn.
        signed: True for difference cubes (Fukui f⁺/f⁻/f⁰), False for raw
            densities.
    """
    if not cube_text:
        return ''
    import json as _json
    cube_json = _json.dumps(cube_text)
    iso = float(isoval)
    if signed:
        return (
            "viewer.addVolumetricData(" + cube_json + ",'cube',"
            "{isoval:" + f"{iso}" + ",color:'#0026ff',opacity:0.80});\n"
            "viewer.addVolumetricData(" + cube_json + ",'cube',"
            "{isoval:" + f"{-iso}" + ",color:'#b00010',opacity:0.80});"
        )
    return (
        "viewer.addVolumetricData(" + cube_json + ",'cube',"
        "{isoval:" + f"{iso}" + ",color:'#0026ff',opacity:0.65});"
    )


def build_fukui_viewer_html(
    xyz_text,
    *,
    labels=None,
    cube_text=None,
    isoval=0.02,
    cube_signed=False,
    width=None,
    height=None,
    viewer_id=None,
):
    """Return a self-contained HTML fragment with a 3Dmol viewer showing a Fukui scene.

    The HTML loads 3Dmol from CDN if not already present, adds the XYZ model,
    applies the standard DELFIN style, optionally adds per-atom labels, and
    optionally adds an isosurface from a cube.

    Honors the global viewer settings; returns a placeholder when the viewer
    is disabled and uses quality-preset dimensions/zoom unless overridden.
    """
    import json as _json
    import random as _random

    profile = get_viewer_profile()
    if not profile['enabled']:
        return viewer_disabled_html()
    eff_width = width if width is not None else profile['width']
    eff_height = height if height is not None else profile['height']
    eff_zoom = profile['zoom']

    if viewer_id is None:
        viewer_id = f"fukui_view_{_random.randint(0, 10**9)}"
    xyz_json = _json.dumps(xyz_text or '')
    label_js = fukui_atom_labels_js(xyz_text, labels) if labels else ''
    cube_js = fukui_cube_isosurface_js(cube_text, isoval=isoval, signed=cube_signed) if cube_text else ''
    style_js = profile['style_js']
    viewer_config_js = profile['viewer_config_js']
    # Reuse the same right-drag-translate patch every other DELFIN viewer
    # installs (orca-builder, calc-browser, remote-archive) so the mouse
    # behaviour stays uniform across the dashboard.
    mouse_patch_js = patch_viewer_mouse_controls_js('viewer', 'el')
    # Slightly larger than CALC_MOL_SIZE (450 px) — user-tuned for the
    # Fukui panel which packs labels at every atom plus optional cubes.
    height_px = int(eff_height) if eff_height else 560
    return f"""
<div id="{viewer_id}"
     style="width:100%;height:{height_px}px;position:relative;"></div>
<script>
(function() {{
    var CURRENT_ID = '{viewer_id}';
    // Tear down any previous Fukui viewer in this document so the new
    // render is the ONLY one painting into the host Output area. We
    // (1) drop the WebGL context, (2) remove the container DOM.
    try {{
        document.querySelectorAll('[id^="fukui_view_"]').forEach(function(node) {{
            if (node.id === CURRENT_ID) return;
            try {{
                var canvas = node.querySelector('canvas');
                if (canvas) {{
                    var gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
                    if (gl && gl.getExtension('WEBGL_lose_context')) {{
                        gl.getExtension('WEBGL_lose_context').loseContext();
                    }}
                }}
            }} catch(_e) {{}}
            try {{
                if (node.parentNode) node.parentNode.removeChild(node);
            }} catch(_e) {{}}
        }});
    }} catch(e) {{}}

    function ensureLoaded(cb) {{
        if (typeof $3Dmol !== 'undefined') {{ cb(); return; }}
        var s = document.createElement('script');
        s.src = 'https://3Dmol.org/build/3Dmol-min.js';
        s.onload = cb;
        document.head.appendChild(s);
    }}
    ensureLoaded(function() {{
        var el = document.getElementById('{viewer_id}');
        if (!el) return;
        var viewer = $3Dmol.createViewer(el, {viewer_config_js});
        viewer.addModel({xyz_json}, 'xyz');
        viewer.setStyle({{}}, {style_js});
        {label_js}
        // Centre + zoom on the ATOM model BEFORE adding volumetric data —
        // cube bounding boxes are usually much larger than the molecule
        // and would otherwise dominate the auto-fit.
        viewer.zoomTo({{model: 0}});
        viewer.center({{model: 0}});
        viewer.zoom({eff_zoom});
        {cube_js}
        viewer.render();
        {mouse_patch_js}
        // Dynamic fit-to-host: the .fukui-viewer-box ancestor knows the
        // current available width/height (set by the parent panel
        // layout); resize the inner div + 3Dmol canvas to match so the
        // viewer fills the frame as the page is resized.
        function fitToHost() {{
            try {{
                var host = el.closest('.fukui-viewer-box') || el.parentElement;
                if (host) {{
                    var w = host.clientWidth || 0;
                    var h = host.clientHeight || 0;
                    if (w >= 200) el.style.width = w + 'px';
                    if (h >= 200) el.style.height = h + 'px';
                }}
                viewer.resize();
                viewer.zoomTo({{model: 0}});
                viewer.center({{model: 0}});
                viewer.zoom({eff_zoom});
                viewer.render();
            }} catch(e) {{}}
        }}
        setTimeout(fitToHost, 40);
        setTimeout(fitToHost, 250);
        if (typeof ResizeObserver !== 'undefined') {{
            try {{
                var ro = new ResizeObserver(fitToHost);
                ro.observe(el);
                var host = el.closest('.fukui-viewer-box');
                if (host) ro.observe(host);
            }} catch(e) {{}}
        }}
        window.addEventListener('resize', fitToHost);
    }});
}})();
</script>
"""


def _build_fukui_table_html(result):
    """Return an HTML table summarizing per-atom Fukui data from ``fukui_result.json``."""
    import html as _html
    atoms = result.get('atoms') or []
    q_neutral = result.get('q_neutral') or []
    q_anion = result.get('q_anion') or []
    q_cation = result.get('q_cation') or []
    f_plus = result.get('f_plus') or []
    f_minus = result.get('f_minus') or []
    f_zero = result.get('f_zero') or []
    n = len(atoms)

    def _cell(value):
        try:
            return f"{float(value):+.3f}"
        except (TypeError, ValueError):
            return _html.escape(str(value))

    rows = []
    for i in range(n):
        rows.append(
            "<tr>"
            f"<td style='text-align:right;padding:2px 8px'>{i + 1}</td>"
            f"<td style='padding:2px 8px'>{_html.escape(str(atoms[i]))}</td>"
            f"<td style='text-align:right;padding:2px 8px;font-family:monospace'>{_cell(q_neutral[i]) if i < len(q_neutral) else ''}</td>"
            f"<td style='text-align:right;padding:2px 8px;font-family:monospace'>{_cell(q_anion[i])   if i < len(q_anion)   else ''}</td>"
            f"<td style='text-align:right;padding:2px 8px;font-family:monospace'>{_cell(q_cation[i])  if i < len(q_cation)  else ''}</td>"
            f"<td style='text-align:right;padding:2px 8px;font-family:monospace;background:#e8f4ff'>{_cell(f_plus[i])    if i < len(f_plus)    else ''}</td>"
            f"<td style='text-align:right;padding:2px 8px;font-family:monospace;background:#ffe8e8'>{_cell(f_minus[i])   if i < len(f_minus)   else ''}</td>"
            f"<td style='text-align:right;padding:2px 8px;font-family:monospace;background:#f0e8ff'>{_cell(f_zero[i])    if i < len(f_zero)    else ''}</td>"
            "</tr>"
        )

    scheme = _html.escape(str(result.get('scheme') or ''))
    settings = result.get('orca_settings') or {}
    settings_line = " · ".join(
        f"{_html.escape(str(k))}={_html.escape(str(v))}"
        for k, v in settings.items() if v not in (None, '')
    )
    geom_origin = _html.escape(str(result.get('geometry_origin') or ''))

    return f"""
<div style="font-family:Arial,sans-serif">
  <div style="margin-bottom:6px">
    <b>Fukui indices</b> · scheme: <code>{scheme}</code> · geometry: <code>{geom_origin}</code>
    <div style="color:#555;font-size:12px">{settings_line}</div>
  </div>
  <table style="border-collapse:collapse;font-size:12px">
    <thead>
      <tr style="background:#f0f0f0">
        <th style="padding:4px 8px">#</th>
        <th style="padding:4px 8px;text-align:left">atom</th>
        <th style="padding:4px 8px">q(N)</th>
        <th style="padding:4px 8px">q(N+1)</th>
        <th style="padding:4px 8px">q(N−1)</th>
        <th style="padding:4px 8px;background:#cfe6ff">f⁺</th>
        <th style="padding:4px 8px;background:#ffcccc">f⁻</th>
        <th style="padding:4px 8px;background:#dccaff">f⁰</th>
      </tr>
    </thead>
    <tbody>
      {"".join(rows)}
    </tbody>
  </table>
</div>
"""


def render_fukui_panel(workdir, output_widget=None):
    """Build a Fukui-results panel as an ``ipywidgets.VBox``.

    The panel owns its own 3D viewer (an ``HTML`` widget whose value is
    swapped on every dropdown change). Decoupling the viewer from the
    surrounding tab's shared ``Output`` widget kills the cross-render
    bleed-through where switching cubes or revisiting ``fukui_result.json``
    left an earlier 3Dmol canvas underneath the new one.

    Args:
        workdir: Path to a Fukui job directory (must contain
            ``fukui_result.json``, ``fukui_geom.xyz``, and the cube files
            produced by ``delfin-fukui``).
        output_widget: Deprecated, ignored. Kept for callers that still
            pass it positionally.

    Returns:
        ``widgets.VBox`` with the controls, the 3D-viewer HTML widget,
        and the per-atom result table. The caller displays it.
    """
    import json as _json
    from pathlib import Path as _Path

    import ipywidgets as widgets

    workdir = _Path(workdir)
    try:
        result = _json.loads((workdir / 'fukui_result.json').read_text(encoding='utf-8'))
    except Exception as exc:
        return widgets.HTML(f"<b style='color:#a00'>Failed to load fukui_result.json:</b> {exc}")

    geom_path = workdir / 'fukui_geom.xyz'
    xyz_text = geom_path.read_text(encoding='utf-8') if geom_path.exists() else ''

    cube_cache = {}

    def _load_cube(name):
        if name is None:
            return None, False
        if name in cube_cache:
            return cube_cache[name]
        for label, fname, signed in _FUKUI_CUBE_OPTIONS:
            if fname == name:
                path = workdir / fname
                if path.exists():
                    try:
                        cube_cache[name] = (path.read_text(encoding='utf-8'), signed)
                        return cube_cache[name]
                    except Exception:
                        pass
        cube_cache[name] = (None, False)
        return cube_cache[name]

    label_select = widgets.Dropdown(
        options=_FUKUI_LABEL_OPTIONS,
        value='f_plus',
        description='Labels:',
        style={'description_width': '60px'},
        layout={'width': '200px'},
    )
    cube_options = [(label, fname) for label, fname, _signed in _FUKUI_CUBE_OPTIONS]
    cube_select = widgets.Dropdown(
        options=cube_options,
        value=None,
        description='Isosurface:',
        style={'description_width': '80px'},
        layout={'width': '320px'},
    )
    iso_slider = widgets.BoundedFloatText(
        value=0.02, min=0.001, max=0.5, step=0.001,
        description='Iso (au):',
        style={'description_width': '70px'},
        layout={'width': '180px'},
    )
    decimals_slider = widgets.BoundedIntText(
        value=3, min=1, max=5, step=1,
        description='Decimals:',
        style={'description_width': '70px'},
        layout={'width': '160px'},
    )

    def _values_for(label_kind):
        if label_kind == 'none' or not label_kind:
            return None
        if label_kind == 'index':
            return [str(i + 1) for i in range(len(result.get('atoms') or []))]
        return result.get(label_kind) or None

    # The viewer is an Output widget OWNED by this panel — it must be
    # an Output (not an HTML widget) because ``widgets.HTML.value``
    # sanitises content and strips ``<script>`` tags, so 3Dmol would
    # never bootstrap. ``with output: display(HTML(...))`` keeps the
    # script tag intact while still scoping the canvas inside the panel,
    # so the canvas is removed together with the panel.
    # Viewer box sized 1:1 like the standard calc-browser 3D viewer
    # (square-ish region, fixed pixel dimensions). flex 0 0 auto pins
    # it on the left of the HBox so the controls sidebar can sit to
    # its right without re-flowing.
    # Same nominal dimensions as the standard calc-browser 3D viewer
    # (CALC_MOL_SIZE = 450 px content height + 2 px border on each side).
    # Width fixed so the table sidebar to the right always fits without
    # overflow; ``flex: 0 0 auto`` pins the viewer at its set width.
    # Own class (NOT ``calc-mol-viewer`` — that would inherit the calc
    # browser's scoping CSS and confuse its JS resize routines).
    # Viewer pinned at fixed width so the sidebar (controls + table)
    # always sits next to it without being pushed off-screen. Inside
    # this fixed box the 3Dmol canvas itself is fluid: the fitToHost
    # JS resizes it whenever the surrounding column changes size, so
    # the viewer-to-canvas relationship is still dynamic.
    viewer_box = widgets.Output(
        layout=widgets.Layout(
            width='560px',
            height='564px',
            min_width='560px',
            min_height='564px',
            margin='0',
            padding='0',
            overflow='hidden',
            border='2px solid #1976d2',
            flex='0 0 auto',
        ),
    )
    viewer_box.add_class('fukui-viewer-box')

    def _redraw(*_):
        from IPython.display import HTML as _HTML, clear_output as _clear, display as _display

        label_kind = label_select.value
        cube_name = cube_select.value
        cube_text, signed = _load_cube(cube_name)
        # Either-or: when a cube is selected, suppress labels and grey
        # out the label dropdown so only the isosurface is shown.
        # Switching the cube dropdown back to "none" re-enables labels.
        label_select.disabled = bool(cube_text)
        decimals_slider.disabled = bool(cube_text)
        iso_slider.disabled = not bool(cube_text)
        if cube_text:
            values = None
            decimals = 0
        else:
            values = _values_for(label_kind)
            decimals = decimals_slider.value if (values and label_kind != 'index') else 0
        iso = iso_slider.value
        html = build_fukui_viewer_html(
            xyz_text,
            labels=values,
            cube_text=cube_text,
            isoval=iso,
            cube_signed=signed,
        )
        if values and label_kind != 'index' and decimals != 3:
            label_js = fukui_atom_labels_js(xyz_text, values, decimals=decimals)
            html = html.replace(
                fukui_atom_labels_js(xyz_text, values, decimals=3),
                label_js,
            )
        with viewer_box:
            _clear(wait=False)
            _display(_HTML(html))

    for control in (label_select, cube_select, iso_slider, decimals_slider):
        control.observe(_redraw, names='value')

    _redraw()

    table_html = _build_fukui_table_html(result)

    # Tighter dropdowns / fields for the side-bar layout (right of viewer).
    for w in (label_select, cube_select):
        w.layout = widgets.Layout(width='240px')
        w.style = {'description_width': '70px'}
    for w in (iso_slider, decimals_slider):
        w.layout = widgets.Layout(width='180px')
        w.style = {'description_width': '70px'}

    controls_block = widgets.VBox(
        [label_select, cube_select, iso_slider, decimals_slider],
        layout=widgets.Layout(
            gap='14px',
            align_items='flex-start',
            margin='0 0 12px 0',
            width='auto',
        ),
    )
    sidebar = widgets.VBox(
        [
            controls_block,
            widgets.HTML(
                table_html,
                layout=widgets.Layout(width='auto', overflow_x='auto'),
            ),
        ],
        layout=widgets.Layout(
            gap='4px',
            align_items='flex-start',
            margin='0 0 0 12px',
            flex='1 1 auto',
            min_width='0',
            width='auto',
            overflow_x='auto',
        ),
    )

    viewer_row = widgets.HBox(
        [viewer_box, sidebar],
        layout=widgets.Layout(
            gap='12px',
            align_items='flex-start',
            justify_content='flex-start',
            flex_flow='row nowrap',
            width='100%',
        ),
    )

    return widgets.VBox(
        [viewer_row],
        layout=widgets.Layout(
            width='100%',
            align_items='flex-start',
            margin='0',
            gap='6px',
            overflow_x='auto',
        ),
    )
