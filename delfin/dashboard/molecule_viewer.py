"""3D molecule visualisation helpers using py3Dmol."""

import py3Dmol
from IPython.display import HTML, clear_output, display

DEFAULT_3DMOL_STYLE = {
    'stick': {
        'colorscheme': 'Jmol',
        'radius': 0.11,
        'singleBonds': False,
        'doubleBondScaling': 0.65,
        'tripleBondScaling': 0.65,
    },
    'sphere': {'colorscheme': 'Jmol', 'scale': 0.28},
}
DEFAULT_3DMOL_STYLE_JS = (
    '{'
    'stick:{colorscheme:"Jmol",radius:0.11,singleBonds:false,doubleBondScaling:0.65,tripleBondScaling:0.65},'
    'sphere:{colorscheme:"Jmol",scale:0.28}'
    '}'
)
DEFAULT_3DMOL_BACKGROUND = 'white'
DEFAULT_3DMOL_ZOOM = 0.90

VIEWER_QUALITY_PROFILES = {
    'low':    {'width': 360, 'height': 280, 'zoom': 0.80, 'show_labels': False},
    'medium': {'width': 480, 'height': 360, 'zoom': 0.85, 'show_labels': True},
    'high':   {'width': 560, 'height': 420, 'zoom': 0.90, 'show_labels': True},
}
_VIEWER_DISABLED_PLACEHOLDER_HTML = (
    "<div style=\"width:100%;max-width:560px;padding:18px 24px;"
    "border:1px dashed #b0b6bf;background:#f6f7f9;color:#4a525c;"
    "font-family:Arial,sans-serif;font-size:13px;border-radius:6px;\">"
    "3D-Viewer ist in den Einstellungen deaktiviert."
    "</div>"
)


def get_viewer_profile():
    """Return the active viewer profile {enabled, quality, width, height, zoom, show_labels}.

    Falls back to the high-quality preset if the user-settings module cannot be
    loaded (e.g. tests that import this module without the full DELFIN env).
    """
    try:
        from delfin.user_settings import load_viewer_settings
        viewer = load_viewer_settings()
    except Exception:
        viewer = {'enabled': True, 'quality': 'high'}
    quality = viewer.get('quality', 'high')
    preset = VIEWER_QUALITY_PROFILES.get(quality, VIEWER_QUALITY_PROFILES['high'])
    return {
        'enabled': bool(viewer.get('enabled', True)),
        'quality': quality,
        **preset,
    }


def viewer_disabled_html(width=None, height=None):
    """Return the placeholder HTML shown when the viewer is disabled."""
    return _VIEWER_DISABLED_PLACEHOLDER_HTML
RIGHT_MOUSE_TRANSLATE_PATCH_JS = (
    '(function(){\n'
    'if(window.__delfinRightDragTranslateSetup) return;\n'
    'window.__delfinRightDragTranslateSetup = true;\n'
    'window.__delfinResolveViewerElement = function(viewer, fallbackEl){\n'
    'try {\n'
    'if (fallbackEl && fallbackEl.nodeType === 1) return fallbackEl;\n'
    'if (viewer && viewer.container && viewer.container.nodeType === 1) return viewer.container;\n'
    '} catch(e) {}\n'
    'return null;\n'
    '};\n'
    'window.__delfinEnableRightDragTranslate = function(viewer, rootEl){\n'
    'try {\n'
    'if(!viewer) return false;\n'
    'var el = window.__delfinResolveViewerElement(viewer, rootEl);\n'
    'if(!el) return false;\n'
    'if(viewer.__delfinRightDragTranslateBound) return true;\n'
    'viewer.__delfinRightDragTranslateBound = true;\n'
    'var dragging = false;\n'
    'var lastX = 0;\n'
    'var lastY = 0;\n'
    'var stopEvt = function(e){\n'
    'if(!e) return;\n'
    'if(typeof e.preventDefault === "function") e.preventDefault();\n'
    'if(typeof e.stopPropagation === "function") e.stopPropagation();\n'
    'if(typeof e.stopImmediatePropagation === "function") e.stopImmediatePropagation();\n'
    '};\n'
    'var onContextMenu = function(e){\n'
    'stopEvt(e);\n'
    'return false;\n'
    '};\n'
    'var onMouseDown = function(e){\n'
    'if(!e || e.button !== 2) return;\n'
    'dragging = true;\n'
    'lastX = e.clientX || 0;\n'
    'lastY = e.clientY || 0;\n'
    'stopEvt(e);\n'
    '};\n'
    'var onMouseMove = function(e){\n'
    'if(!dragging) return;\n'
    'var x = e.clientX || 0;\n'
    'var y = e.clientY || 0;\n'
    'var dx = x - lastX;\n'
    'var dy = y - lastY;\n'
    'var tdy = -dy;\n'
    'lastX = x;\n'
    'lastY = y;\n'
    'try {\n'
    'if(typeof viewer.translate === "function") {\n'
    'viewer.translate(dx, tdy);\n'
    '} else if (typeof viewer.translateScene === "function") {\n'
    'viewer.translateScene(dx, tdy);\n'
    '} else if (typeof viewer.pan === "function") {\n'
    'viewer.pan(dx, tdy);\n'
    '}\n'
    'if(typeof viewer.render === "function") viewer.render();\n'
    '} catch(_e) {}\n'
    'stopEvt(e);\n'
    '};\n'
    'var onMouseUp = function(e){\n'
    'if(!dragging) return;\n'
    'dragging = false;\n'
    'stopEvt(e);\n'
    '};\n'
    'el.addEventListener("contextmenu", onContextMenu, true);\n'
    'el.addEventListener("mousedown", onMouseDown, true);\n'
    'window.addEventListener("mousemove", onMouseMove, true);\n'
    'window.addEventListener("mouseup", onMouseUp, true);\n'
    'viewer.__delfinRightDragTranslateCleanup = function(){\n'
    'try {\n'
    'el.removeEventListener("contextmenu", onContextMenu, true);\n'
    'el.removeEventListener("mousedown", onMouseDown, true);\n'
    'window.removeEventListener("mousemove", onMouseMove, true);\n'
    'window.removeEventListener("mouseup", onMouseUp, true);\n'
    '} catch(e) {}\n'
    '};\n'
    'return true;\n'
    '} catch(e) {\n'
    'return false;\n'
    '}\n'
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
    'var scopeMaps = [window._calcMolViewerByScope, window._calcTrajViewerByScope];\n'
    'for (var i = 0; i < scopeMaps.length; i++) {\n'
    'var map = scopeMaps[i];\n'
    'if(!map || typeof map !== "object") continue;\n'
    'for (var key in map) {\n'
    'if(!Object.prototype.hasOwnProperty.call(map, key)) continue;\n'
    'var sv = map[key];\n'
    'if(sv && typeof sv === "object") window.__delfinEnableRightDragTranslate(sv, null);\n'
    '}\n'
    '}\n'
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
    'window.__delfinPatchAllKnown3DmolViewers();\n'
    'attachPromiseHook();\n'
    'if(!window.__delfinRightDragTranslateTimer){\n'
    'window.__delfinRightDragTranslateTimer = window.setInterval(function(){\n'
    'window.__delfinPatchAllKnown3DmolViewers();\n'
    'attachPromiseHook();\n'
    '}, 250);\n'
    '}\n'
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

    function getViewer(scopeKey) {
        return (window._submitMolViewerByScope || {})[scopeKey] || null;
    }
    function getState(scopeKey) {
        var s = window._submitManipStateByScope[scopeKey];
        if (!s) {
            s = {
                mode: 'off',
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
        if (!viewer) { try { console.log('delfin push: no viewer'); } catch (e) {} return; }
        var input = getSyncInput(scopeKey);
        if (!input) { try { console.log('delfin push: no sync input'); } catch (e) {} return; }
        var xyz = serializeXyz(viewer);
        var proto = (input.tagName === 'TEXTAREA')
            ? window.HTMLTextAreaElement.prototype
            : window.HTMLInputElement.prototype;
        var setter = Object.getOwnPropertyDescriptor(proto, 'value');
        if (setter && setter.set) setter.set.call(input, xyz);
        else input.value = xyz;
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
        try {
            console.log('delfin push: sync input', input.tagName,
                'len=', xyz.length, 'firstLine=', xyz.split('\n')[0]);
        } catch (e) {}
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

    // --- 3Dmol-native atom picking ---
    function attachClickable(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) {
            try { console.log('delfin attachClickable: no viewer for', scopeKey); } catch (e) {}
            return;
        }
        try { console.log('delfin attachClickable:', scopeKey,
            'hasSetClickable=', typeof viewer.setClickable); } catch (e) {}
        try {
            viewer.setClickable({}, true, function(atom, v, ev, container) {
                try { console.log('delfin setClickable fired:', atom && atom.serial); } catch (e) {}
                var state = getState(scopeKey);
                if (state.mode !== 'select') return;
                if (state.drag && state.drag.movedEnough) return;
                if (!atom || atom.serial === undefined) return;
                var additive = !!(ev && (ev.shiftKey || ev.ctrlKey || ev.metaKey));
                togglePick(scopeKey, atom, additive);
            });
            // Force 3Dmol to rebuild any internal pick buffers.
            try { viewer.render(); } catch (e) {}
        } catch (e) {
            try { console.log('delfin setClickable error:', e.message); } catch (_) {}
        }
        installCanvasClickFallback(scopeKey);
    }
    function detachClickable(scopeKey) {
        var viewer = getViewer(scopeKey);
        if (!viewer) return;
        try { viewer.setClickable({}, false, function(){}); } catch (e) {}
        uninstallCanvasClickFallback(scopeKey);
    }

    // Dispatch a single synthetic click at (clientX, clientY), collect the atom
    // that 3Dmol's setClickable reports. Restores any prior callback afterwards.
    function probeClickAtom(scopeKey, clientX, clientY) {
        var viewer = getViewer(scopeKey);
        var state = getState(scopeKey);
        if (!viewer || !state.canvas) return null;
        var hit = null;
        try {
            viewer.setClickable({}, true, function(atom) {
                if (atom && atom.serial !== undefined) {
                    hit = {serial: atom.serial, elem: atom.elem || 'X'};
                }
            });
            try { viewer.render(); } catch (e) {}
            var opts = {
                bubbles: true, cancelable: true,
                clientX: clientX, clientY: clientY, button: 0, buttons: 1
            };
            state.canvas.dispatchEvent(new MouseEvent('mousedown', opts));
            opts.buttons = 0;
            state.canvas.dispatchEvent(new MouseEvent('mouseup', opts));
            state.canvas.dispatchEvent(new MouseEvent('click', opts));
        } catch (e) {
            try { console.log('delfin probeClickAtom error:', e.message); } catch (_) {}
        }
        // Restore the select-mode callback (even in manipulate, harmless)
        if (state.mode === 'select') {
            attachClickable(scopeKey);
        } else {
            try { viewer.setClickable({}, false, function(){}); } catch (e) {}
        }
        try { console.log('delfin probeClickAtom @', clientX, clientY, '→', hit && hit.serial); } catch (e) {}
        return hit;
    }

    // Fallback: direct click listener on canvas that uses raycastAtom.
    // Runs in addition to 3Dmol's setClickable so we pick atoms even if
    // 3Dmol's internal click detection is broken in this build.
    function installCanvasClickFallback(scopeKey) {
        var state = getState(scopeKey);
        var canvas = state.canvas;
        if (!canvas) {
            try { console.log('delfin: no canvas for fallback'); } catch (e) {}
            return;
        }
        if (state._canvasClickHandler) {
            try { canvas.removeEventListener('click', state._canvasClickHandler, true); } catch (e) {}
        }
        var handler = function(e) {
            try { console.log('delfin canvas click @', e.clientX, e.clientY); } catch (_) {}
            if (state.mode !== 'select') return;
            if (state.shiftHeld) return; // shift+drag path handles rect
            if (state.drag && state.drag.movedEnough) return;
            var atom = raycastAtom(scopeKey, e.clientX, e.clientY);
            try { console.log('delfin raycast hit:', atom && atom.serial); } catch (_) {}
            if (atom) {
                var additive = !!(e.shiftKey || e.ctrlKey || e.metaKey);
                togglePick(scopeKey, atom, additive);
            }
        };
        canvas.addEventListener('click', handler, true);
        state._canvasClickHandler = handler;
        try { console.log('delfin canvas click fallback installed'); } catch (e) {}
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

    function raycastAtom(scopeKey, clientX, clientY) {
        var viewer = getViewer(scopeKey);
        var state = getState(scopeKey);
        if (!viewer || !state.canvas) return null;
        var rect = state.canvas.getBoundingClientRect();
        var sx = clientX - rect.left;
        var sy = clientY - rect.top;
        var atoms = getAtoms(viewer);
        var best = null, bestDist = 22*22; // 22 px pick radius
        for (var i = 0; i < atoms.length; i++) {
            var p = modelToScreen(viewer, state.canvas, atoms[i]);
            if (!p) continue;
            var dx = p.x - sx, dy = p.y - sy;
            var d2 = dx*dx + dy*dy;
            if (d2 < bestDist) { bestDist = d2; best = atoms[i]; }
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
    // Rect selection by probing: dispatch synthetic click events at a grid of
    // points inside the rect. 3Dmol's setClickable callback fires for each hit
    // atom. This bypasses 3D→2D projection entirely and uses 3Dmol's own
    // internal pick buffer, which we know works (individual clicks pick atoms).
    function finishRect(scopeKey, x0, y0, x1, y1, additive) {
        var state = getState(scopeKey);
        if (state.rect) { try { state.overlay.removeChild(state.rect); } catch(e) {} state.rect = null; }
        var viewer = getViewer(scopeKey);
        if (!viewer || !state.canvas) return;
        var canvas = state.canvas;

        var minX = Math.min(x0, x1), maxX = Math.max(x0, x1);
        var minY = Math.min(y0, y1), maxY = Math.max(y0, y1);
        if (maxX - minX < 3 || maxY - minY < 3) return;

        if (!additive) state.picks = [];

        // Swap setClickable callback to a collector.
        var collected = {};
        try {
            viewer.setClickable({}, true, function(atom) {
                if (atom && atom.serial !== undefined) {
                    collected[atom.serial] = {
                        serial: atom.serial, elem: atom.elem || 'X'
                    };
                }
            });
            try { viewer.render(); } catch (e) {}
        } catch (e) {
            try { console.log('delfin rect: setClickable swap failed:', e.message); } catch (_) {}
            return;
        }

        // Step = ~half of a small atom radius on screen. 6–8px works for most
        // zoom levels; trade-off between completeness and speed.
        var step = 7;
        var points = 0;
        function dispatchClick(x, y) {
            var md = new MouseEvent('mousedown', {
                bubbles: true, cancelable: true,
                clientX: x, clientY: y, button: 0, buttons: 1
            });
            var mu = new MouseEvent('mouseup', {
                bubbles: true, cancelable: true,
                clientX: x, clientY: y, button: 0, buttons: 0
            });
            var cl = new MouseEvent('click', {
                bubbles: true, cancelable: true,
                clientX: x, clientY: y, button: 0
            });
            canvas.dispatchEvent(md);
            canvas.dispatchEvent(mu);
            canvas.dispatchEvent(cl);
        }
        for (var y = minY; y <= maxY; y += step) {
            for (var x = minX; x <= maxX; x += step) {
                points++;
                dispatchClick(x, y);
            }
        }

        // Restore normal callback
        attachClickable(scopeKey);

        // Merge collected atoms into picks
        var added = 0;
        Object.keys(collected).forEach(function(s) {
            var serial = +s;
            var exists = state.picks.some(function(q) { return q.serial === serial; });
            if (!exists) {
                state.picks.push(collected[s]);
                added++;
            }
        });

        try {
            console.log('delfin rect probe:',
                'points=', points,
                'hits=', Object.keys(collected).length,
                'added=', added);
        } catch (e) {}
        redrawHighlights(scopeKey);
    }

    function applyTranslate(scopeKey, deltaWorld) {
        var viewer = getViewer(scopeKey);
        var state = getState(scopeKey);
        if (!viewer) return;
        var atoms = getAtoms(viewer);
        var byS = {};
        for (var i = 0; i < atoms.length; i++) byS[atoms[i].serial] = atoms[i];
        state.picks.forEach(function(p) {
            var a = byS[p.serial];
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
            try { console.log('delfin rotate: no viewer/pivot'); } catch (_) {}
            return;
        }
        var pivotAtom = getAtomBySerial(viewer, state.pivot.serial);
        if (!pivotAtom) {
            try { console.log('delfin rotate: pivot atom serial', state.pivot.serial, 'not found'); } catch (_) {}
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
        try { console.log('delfin rotate: moved', moved, 'yaw=', yaw.toFixed(3), 'pitch=', pitch.toFixed(3)); } catch (_) {}
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
                e.preventDefault(); e.stopPropagation();
                try { console.log('delfin manip mousedown button=', e.button, 'picks=', state.picks.length); } catch (_) {}
                if (e.button === 2) {
                    // Right mouse always tries to pick pivot atom under cursor.
                    var picked = probeClickAtom(scopeKey, e.clientX, e.clientY);
                    try { console.log('delfin pivot probe →', picked && picked.serial); } catch (_) {}
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
                        try { console.log('delfin rotate drag started, pivot=', state.pivot.serial); } catch (_) {}
                    }
                    return;
                }
                if (state.picks.length === 0) return;
                if (e.button === 0) {
                    state.drag = {
                        kind: 'translate',
                        startX: e.clientX, startY: e.clientY,
                        lastX: e.clientX, lastY: e.clientY,
                        movedEnough: false, snapshotted: false
                    };
                }
                return;
            }

            // Select mode — we fully capture the mouse so 3Dmol doesn't rotate.
            if (e.button === 0) {
                e.preventDefault();
                e.stopPropagation();
                state.drag = {
                    kind: 'maybe-rect',
                    startX: e.clientX, startY: e.clientY,
                    origX: x, origY: y,
                    additive: !!(e.shiftKey || e.ctrlKey || e.metaKey),
                    movedEnough: false,
                    hitAtom: !!atom,
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
                applyTranslate(scopeKey, delta);
            } else if (d.kind === 'rotate' && d.movedEnough) {
                e.preventDefault();
                if (!d.snapshotted) { snapshotForUndo(scopeKey); d.snapshotted = true; }
                applyRotate(scopeKey, dx, dy);
            } else if (d.kind === 'maybe-rect' && d.movedEnough && !d.hitAtom) {
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
            state.overlay.style.pointerEvents = 'auto';
            state.overlay.style.cursor = state.picks.length ? 'move' : 'not-allowed';
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
                        try { console.log('delfin window-level rotate drag, pivot=', s.pivot.serial); } catch (_) {}
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

    window.__delfinSubmitManip = {
        onViewerReady: onViewerReady,
        setMode: setMode,
        clear: clearPicks,
        undo: undo
    };
})();
"""


def submit_manip_bootstrap_js():
    """Return one-time JS that installs window.__delfinSubmitManip helpers."""
    return SUBMIT_MANIP_BOOTSTRAP_JS


def apply_molecule_view_style(view, zoom=DEFAULT_3DMOL_ZOOM):
    """Apply a shared ChemDarwin/MSILES-like style to a py3Dmol viewer."""
    if hasattr(view, 'startjs'):
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
                    'var __delfinRecenter=function(){'
                    'try{'
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
                    '})();'
                )
                + '\n'
            )
    view.setStyle({}, DEFAULT_3DMOL_STYLE)
    view.setBackgroundColor(DEFAULT_3DMOL_BACKGROUND)
    view.zoomTo()
    view.center()
    if zoom is not None:
        view.zoom(zoom)
    view.render()
    return view


def render_xyz_in_output(output_widget, xyz_text, width=None, height=None):
    """Render an XYZ string inside an ipywidgets Output widget.

    Honors the global viewer settings (``ui.viewer.enabled`` and
    ``ui.viewer.quality``). When ``width``/``height`` are passed explicitly
    they override the quality preset.
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
        apply_molecule_view_style(view, zoom=profile['zoom'])
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
    style_js = DEFAULT_3DMOL_STYLE_JS
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
        var viewer = $3Dmol.createViewer(el, {{ backgroundColor: 'white' }});
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
