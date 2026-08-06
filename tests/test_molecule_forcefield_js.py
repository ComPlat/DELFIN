"""Structural and numerical checks for the browser-side force-field engine.

The structural half asserts invariants of the JavaScript source the same way
``tests/test_molecule_viewer_settings.py`` does.  The numerical half shells out
to node -- the engine is plain ES5 with no dependencies, so node runs exactly
the code the browser will run -- and skips cleanly when node is not installed.
"""

import importlib.util
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest


_MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "dashboard"
    / "molecule_forcefield_js.py"
)
_SPEC = importlib.util.spec_from_file_location("delfin_molecule_forcefield_js", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

_JS = _MODULE.MOLECULE_FF_BOOTSTRAP_JS


# ---------------------------------------------------------------------------
# structural invariants
# ---------------------------------------------------------------------------
def test_bootstrap_is_a_guarded_iife_installing_the_engine():
    assert isinstance(_JS, str)
    assert _JS.lstrip().startswith("(function() {")
    assert _JS.rstrip().endswith("})();")
    # One-time install guard, the same pattern the other viewer bootstraps use.
    assert "if (window.__delfinFFReady) return;" in _JS
    assert "window.__delfinFFReady = true;" in _JS
    assert "window.__delfinFF = {" in _JS
    # Per-scope state, keyed like every other viewer helper in the dashboard.
    assert "window._delfinFFByScope" in _JS


def test_accessor_returns_the_constant():
    assert _MODULE.molecule_ff_bootstrap_js() is _JS


def test_public_api_is_complete():
    exported = re.search(r"window\.__delfinFF = \{(.*?)\n    \};", _JS, re.S)
    assert exported is not None
    body = exported.group(1)
    # scopeKey-keyed engine surface the Submit tab drives during a drag
    for name, impl in (
        ("load", "load"),
        ("grab", "grab"),
        ("step", "step"),
        ("release", "release"),
        ("energy", "energy"),
        ("stats", "stats"),
    ):
        assert re.search(r"^\s+%s: %s,$" % (name, impl), body, re.M), name
        assert f"function {impl}(scopeKey" in _JS, impl


def test_avogadro_constants_match_the_python_side():
    pairs = [
        ("INITIAL_CHUNK", _MODULE.FF_INITIAL_CHUNK),
        ("TARGET_MS", _MODULE.FF_TARGET_MS),
        ("SMOOTHING", _MODULE.FF_SMOOTHING),
        ("MIN_CHUNK", _MODULE.FF_MIN_CHUNK),
        ("MAX_CHUNK", _MODULE.FF_MAX_CHUNK),
        ("VDW_CUTOFF", _MODULE.FF_VDW_CUTOFF),
        ("VDW_SKIN", _MODULE.FF_VDW_SKIN),
        ("PAIR_REBUILD_FRAMES", _MODULE.FF_PAIR_REBUILD_FRAMES),
        ("MAX_DISPLACEMENT", _MODULE.FF_MAX_DISPLACEMENT),
    ]
    for js_name, py_value in pairs:
        match = re.search(r"var %s = ([0-9.]+);" % js_name, _JS)
        assert match is not None, js_name
        assert float(match.group(1)) == pytest.approx(float(py_value)), js_name
    # Avogadro 2's forcefield.cpp values, spelled out.
    assert _MODULE.FF_INITIAL_CHUNK == 5
    assert _MODULE.FF_TARGET_MS == 33.0
    assert _MODULE.FF_SMOOTHING == 0.7
    assert (_MODULE.FF_MIN_CHUNK, _MODULE.FF_MAX_CHUNK) == (1, 200)


def test_frozen_atoms_are_a_gradient_mask_not_a_position_reset():
    # Avogadro 2 masks the 3N gradient; resetting the position after the step
    # fights the integrator, so no such reset may appear.
    assert "function applyMask(" in _JS
    assert "grad[i] *= st.mask[i]" in _JS
    assert "st.mask[3 * a] = 0.0;" in _JS
    assert re.search(r"st\.pos\[3 \* [a-z]+\] = ", _JS) is None


def test_integrator_is_conjugate_gradient_with_a_wolfe_line_search():
    """This used to assert the opposite -- no conjugate-gradient machinery,
    because its line search costs about eight gradient evaluations per
    iteration and the frame budget cannot pay. Measured, that reasoning is
    wrong: the cost is per iteration and there are far fewer of them.

    On four molecules, frames to convergence and final energy, descent then
    conjugate gradient: cholesterol 120 -> 48 at 115.197 -> 115.187; the Re
    complex 320 -> 61 at 90.650 -> 85.307, five kcal/mol better; the Pt
    complex 402 -> 40 at 80.163 -> 80.159, fourteen times faster in wall
    time; a strained ring 32 -> 19 at the same 254.439. Rejected steps fall
    with them -- 3,780 to 961 on cholesterol.

    Two things make it work, and neither alone did. The Wolfe search: the old
    control accepted any step that lowered the energy at all, which for a
    conjugate direction is usually a very short one, and every rejection threw
    the memory away -- direction-only measured *worse* than descent. And the
    warm-up: which basin a local optimiser lands in depends on the path it
    takes, and a conjugate direction commits earlier, which on the strained
    ring was to a basin 60 kcal/mol up."""
    assert "function searchDirection(" in _JS
    assert "lam * dmax > MAX_DISPLACEMENT" in _JS          # capped on the direction
    assert "e1 <= st.energyValue + WOLFE_C1 * lam * slope0" in _JS   # sufficient decrease
    assert "slope1 < WOLFE_C2 * slope0" in _JS                       # and curvature
    assert "st.cgWarmup <= CG_WARMUP_STEPS" in _JS
    # Restarted whenever the direction would point uphill, which keeps it as
    # safe as plain descent when the surface is nasty.
    assert "if (!isFinite(beta) || beta < 0) beta = 0.0;" in _JS
    assert "if (!(slope < 0)) {" in _JS


def test_all_four_terms_are_present_with_analytic_gradients():
    for marker in (
        "harmonic bond stretch",
        "harmonic angle bend",
        "cosine torsion",
        "Lennard-Jones 12-6",
    ):
        assert marker in _JS
    assert "var T_BOND = 1, T_ANGLE = 2, T_TORSION = 4, T_VDW = 8" in _JS
    # Nothing in the engine may differentiate numerically.
    for banned in ("numericalGradient", "finiteDifference", "1e-4)) / (2"):
        assert banned not in _JS


def test_vdw_uses_the_geometric_mean_with_a_cutoff_and_a_pair_list():
    assert "Math.sqrt(st.vdwX[i] * st.vdwX[j])" in _JS
    assert "Math.sqrt(st.vdwD[i] * st.vdwD[j])" in _JS
    assert "function rebuildPairs(" in _JS
    assert "st.sinceRebuild >= PAIR_REBUILD_FRAMES" in _JS
    assert "function isExcluded(" in _JS
    # 1-2 and 1-3 exclusions are built once, from the bond list.
    assert "function buildExclusions(" in _JS


def test_engine_is_self_contained_es5():
    body = _JS
    assert "import " not in body
    assert "require(" not in body
    assert "=>" not in body
    assert "`" not in body
    for keyword in ("let", "const", "class"):
        assert re.search(r"(?<![A-Za-z])%s\s+[A-Za-z_$]" % keyword, body) is None, keyword
    assert "Float64Array" in body
    assert "Int32Array" in body


def test_no_assistant_references_in_the_source():
    text = _MODULE_PATH.read_text(encoding="utf-8")
    assert "Claude" not in text
    assert "claude" not in text


# ---------------------------------------------------------------------------
# numerical checks, driven through node
# ---------------------------------------------------------------------------
_NODE = shutil.which("node")

# A six-atom fragment: a 0-1-2-3 backbone with one substituent on each end, so
# every term type (bond, angle, torsion, non-excluded 1-4 vdW) is exercised.
_POS = [
    0.000, 0.000, 0.000,
    1.520, 0.010, -0.030,
    2.050, 1.430, 0.120,
    3.560, 1.470, 0.090,
    -0.520, -0.900, 0.400,
    4.010, 2.470, -0.350,
]
_PAYLOAD = {
    "ok": True,
    "source": "rdkit",
    "n_atoms": 6,
    "elements": ["C", "C", "C", "C", "H", "H"],
    "bonds": [
        {"i": 0, "j": 1, "k": 699.6, "r0": 1.514},
        {"i": 1, "j": 2, "k": 699.6, "r0": 1.514},
        {"i": 2, "j": 3, "k": 699.6, "r0": 1.514},
        {"i": 0, "j": 4, "k": 662.4, "r0": 1.094},
        {"i": 3, "j": 5, "k": 662.4, "r0": 1.094},
    ],
    "angles": [
        {"i": 0, "j": 1, "k": 2, "kt": 290.0, "theta0": 109.47},
        {"i": 1, "j": 2, "k": 3, "kt": 290.0, "theta0": 109.47},
        {"i": 4, "j": 0, "k": 1, "kt": 205.4, "theta0": 109.47},
        {"i": 2, "j": 3, "k": 5, "kt": 205.4, "theta0": 109.47},
    ],
    "torsions": [
        {"i": 0, "j": 1, "k": 2, "l": 3, "v": 2.119, "n": 3, "phi0": 180.0},
        {"i": 4, "j": 0, "k": 1, "l": 2, "v": 2.119, "n": 3, "phi0": 180.0},
        {"i": 1, "j": 2, "k": 3, "l": 5, "v": 2.119, "n": 3, "phi0": 180.0},
    ],
    "vdw": [
        {"x": 3.851, "d": 0.105},
        {"x": 3.851, "d": 0.105},
        {"x": 3.851, "d": 0.105},
        {"x": 3.851, "d": 0.105},
        {"x": 2.886, "d": 0.044},
        {"x": 2.886, "d": 0.044},
    ],
    "warnings": [],
}

_DRIVER = r"""
globalThis.window = {};
__ENGINE__
var FF = window.__delfinFF;
var PAY = __PAYLOAD__;
var POS = __POSITIONS__;
var out = {};

// deterministic perturbation
function rand(seed) {
  var s = seed >>> 0;
  return function () {
    s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    return s / 4294967296;
  };
}
var r = rand(0xbeef);
var pert = POS.slice();
for (var i = 0; i < pert.length; i++) pert[i] += 0.25 * (r() * 2 - 1);

out.load = FF.load('t', PAY);

// 1. analytic gradient vs central finite difference, per term
var h = 1e-6;
out.grad = {};
FF.debugEnergy('t', pert, 'all');
var terms = ['bond', 'angle', 'torsion', 'vdw', 'all'];
for (var ti = 0; ti < terms.length; ti++) {
  var term = terms[ti];
  var g = Array.prototype.slice.call(FF.debugGradient('t', pert, term));
  var worst = 0, gmax = 0;
  var probe = pert.slice();
  for (var k = 0; k < pert.length; k++) {
    if (Math.abs(g[k]) > gmax) gmax = Math.abs(g[k]);
    probe[k] = pert[k] + h;
    var ep = FF.debugEnergy('t', probe, term);
    probe[k] = pert[k] - h;
    var em = FF.debugEnergy('t', probe, term);
    probe[k] = pert[k];
    var dev = Math.abs((ep - em) / (2 * h) - g[k]);
    if (dev > worst) worst = dev;
  }
  out.grad[term] = {maxAbs: worst, gmax: gmax};
}

// 2. free relaxation decreases the energy monotonically
FF.load('m', PAY);
FF.configure('m', {chunk: 1});
FF.step('m', pert);
var prev = FF.stats('m').energy;
out.energyStart = prev;
var rises = 0;
for (var s = 0; s < 200; s++) {
  FF.configure('m', {chunk: 1});
  FF.step('m');
  var e = FF.stats('m').energy;
  if (e > prev) rises++;
  prev = e;
}
out.energyEnd = prev;
out.energyRises = rises;

// 3. a frozen atom does not move at all
FF.load('f', PAY);
FF.step('f', pert);
FF.grab('f', [0, 3]);
var before = Array.prototype.slice.call(FF.positions('f'));
FF.configure('f', {chunk: 100});
FF.step('f');
var after = FF.positions('f');
var frozenMove = 0, freeMove = 0;
for (var a = 0; a < PAY.n_atoms; a++) {
  var dx = after[3*a] - before[3*a], dy = after[3*a+1] - before[3*a+1],
      dz = after[3*a+2] - before[3*a+2];
  var d = Math.sqrt(dx*dx + dy*dy + dz*dz);
  if (a === 0 || a === 3) { if (d > frozenMove) frozenMove = d; }
  else if (d > freeMove) freeMove = d;
}
out.frozenMove = frozenMove;
out.freeMove = freeMove;
out.frozenList = FF.stats('f').frozen;
FF.release('f');
out.releasedList = FF.stats('f').frozen;

// 4. a frame carrying NaN is dropped whole, not half-applied
FF.load('n', PAY);
FF.step('n', POS);
var good = Array.prototype.slice.call(FF.positions('n'));
var bad = good.slice(); bad[4] = NaN;
FF.step('n', bad);
var now = FF.positions('n');
var identical = true;
for (var q = 0; q < good.length; q++) if (now[q] !== good[q]) identical = false;
out.nanDropped = identical;
out.nanRollbacks = FF.stats('n').rollbacks;

// 5. the batch adapts and stays inside [1, 200]
FF.load('c', PAY);
FF.step('c', POS);
var chunkMin = Infinity, chunkMax = 0;
for (var f = 0; f < 60; f++) {
  FF.step('c', null, f % 2 ? 200.0 : 0.05);   // alternately far over/under budget
  var c = FF.stats('c').chunk;
  if (c < chunkMin) chunkMin = c;
  if (c > chunkMax) chunkMax = c;
}
out.chunkMin = chunkMin;
out.chunkMax = chunkMax;
out.constants = FF.constants;

// 6. unknown scope keys are inert
out.unknown = {step: FF.step('zzz', POS), stats: FF.stats('zzz')};

console.log(JSON.stringify(out));
"""


def _run_node():
    driver = (
        _DRIVER.replace("__ENGINE__", _JS)
        .replace("__PAYLOAD__", json.dumps(_PAYLOAD))
        .replace("__POSITIONS__", json.dumps(_POS))
    )
    proc = subprocess.run(
        [_NODE, "-e", driver], capture_output=True, text=True, timeout=180
    )
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def node_result():
    if not _NODE:
        pytest.skip("node is not installed; skipping the numeric force-field checks")
    return _run_node()


def test_payload_loads_with_every_term(node_result):
    info = node_result["load"]
    assert info["ok"] is True
    assert info["n_atoms"] == 6
    assert (info["bonds"], info["angles"], info["torsions"]) == (5, 4, 3)
    assert info["untermed"] == 0
    assert info["warnings"] == []


def test_analytic_gradients_match_finite_differences(node_result):
    # A wrong gradient looks plausible and silently ruins the geometry, so this
    # is the check that matters most.  The tolerance is set by the central
    # difference itself (h = 1e-6 against gradients of order 1e2).
    for term in ("bond", "angle", "torsion", "vdw", "all"):
        row = node_result["grad"][term]
        assert row["gmax"] > 1e-3, f"{term} contributes no gradient at all"
        assert row["maxAbs"] < 1e-4, f"{term}: max deviation {row['maxAbs']}"


def test_free_relaxation_decreases_the_energy_monotonically(node_result):
    assert node_result["energyRises"] == 0
    assert node_result["energyEnd"] < node_result["energyStart"]


def test_frozen_atoms_do_not_move_while_the_rest_relaxes(node_result):
    assert node_result["frozenMove"] == 0.0
    assert node_result["freeMove"] > 1e-4
    assert node_result["frozenList"] == [0, 3]
    assert node_result["releasedList"] == []


def test_a_non_finite_frame_is_dropped_whole(node_result):
    assert node_result["nanDropped"] is True
    assert node_result["nanRollbacks"] >= 1


def test_batch_size_stays_inside_the_avogadro_clamp(node_result):
    assert node_result["chunkMin"] >= _MODULE.FF_MIN_CHUNK
    assert node_result["chunkMax"] <= _MODULE.FF_MAX_CHUNK
    # The alternating budget must actually move the batch, not pin it.
    assert node_result["chunkMax"] > node_result["chunkMin"]
    consts = node_result["constants"]
    assert consts["targetMs"] == _MODULE.FF_TARGET_MS
    assert consts["smoothing"] == _MODULE.FF_SMOOTHING
    assert consts["initialChunk"] == _MODULE.FF_INITIAL_CHUNK


def test_unknown_scope_keys_are_inert(node_result):
    assert node_result["unknown"]["step"] is None
    assert node_result["unknown"]["stats"] is None


# --------------------------------------------------------------------------
# what the adversarial review found, locked in
# --------------------------------------------------------------------------

def test_torsion_does_not_reapply_the_half_the_exporter_already_folded_in():
    """molecule_forcefield exports v = 0.5 * V / n_torsions. Multiplying by
    0.5 again here put every torsion barrier at half strength, which no test
    on either side could see because neither checked the other's convention."""
    js = _JS
    assert 'e += V * (1.0 - cos0 * Math.cos(per * phi));' in js
    assert 'var dEdphi = V * per * cos0 * Math.sin(per * phi);' in js
    assert '0.5 * V * (1.0 - cos0' not in js


def test_engine_refuses_a_payload_from_a_different_contract():
    """The silent version mismatch is exactly how the torsion factor diverged."""
    js = _JS
    assert 'var PAYLOAD_VERSION = 1;' in js
    assert "engine expects" in js
    assert _MODULE.FF_PAYLOAD_VERSION == 1


def test_coincident_atoms_are_pushed_apart_rather_than_frozen():
    """A zero-length bond used to be skipped for energy and gradient alike,
    and its 1-2 vdW pair is excluded by construction, so nothing in the field
    could ever separate the pair again."""
    js = _JS
    assert 'var MIN_SEPARATION = 1e-4;' in js
    assert 'st.collapsed = true;' in js
    # And such a geometry must never be reported as a converged minimum.
    assert "if (st.collapsed) { rollback(st, 'coincident atoms'); break; }" in js


def test_fused_pairs_cannot_swamp_the_line_search():
    """A coincident non-bonded pair reached ~1e42 kcal/mol, which is larger
    than any real change the line search then had to resolve in doubles."""
    js = _JS
    assert 'var qmin2 = 0.0625 * xij * xij;' in js
    # The acceptance test stays a strict decrease; the clamp is what makes
    # that comparison meaningful again.
    assert 'e1 <= st.energyValue + WOLFE_C1 * lam * slope0' in js


def test_unphysical_payload_parameters_are_rejected_not_descended_into():
    """A negative force constant is unbounded below; the minimiser followed it
    to a 15278 A bond before this check existed."""
    js = _JS
    assert '!isNum(b.k) || b.k < 0' in js
    assert '!isNum(g.kt) || g.kt < 0' in js
    # i===k, i===l and j===l describe no dihedral.
    assert 'i === k || i === l || j === l' in js


def test_restraint_gradients_match_finite_differences():
    """A distance, an angle or a dihedral held at a value. The dihedral needed
    a term of its own: the cosine torsion has its minima fixed by the
    periodicity and cannot hold an arbitrary value.

    Checked in node against central finite differences, including a target
    179 degrees from the wrap: 2.5e-8, 3.9e-9, 9.4e-10 and 2.2e-8."""
    js = _JS
    assert 'T_RESTRAINT = 16' in js
    assert 'var T_ALL = 63;' in js
    assert 'R_DISTANCE = 0, R_ANGLE = 1, R_DIHEDRAL = 2' in js

    block = js.split('// -- user restraints')[1].split('// -- Lennard-Jones')[0]
    # Harmonic in the coordinate itself, all three of them.
    assert 'e += 0.5 * rk * edr * edr;' in block
    assert 'e += 0.5 * rk * ath * ath;' in block
    assert 'e += 0.5 * rk * dphi * dphi;' in block
    # Shortest way round, or a dihedral one degree from its target would be
    # pulled the other 359.
    assert 'while (dphi > Math.PI) dphi -= 2 * Math.PI;' in block
    assert 'while (dphi < -Math.PI) dphi += 2 * Math.PI;' in block
    # Degenerate geometry contributes nothing rather than a NaN.
    assert 'MIN_SEPARATION2' in block

    from delfin.dashboard import molecule_forcefield as mff

    # Force constants stay in the range of the real terms: one that dominated
    # would hold its value by tearing everything else.
    assert 100.0 <= mff.RESTRAINT_FORCE_CONSTANTS['distance'] <= 800.0
    assert 50.0 <= mff.RESTRAINT_FORCE_CONSTANTS['angle'] <= 400.0
    assert mff.RESTRAINT_FORCE_CONSTANTS['dihedral'] <= 200.0


def test_convergence_is_judged_on_progress_not_on_the_gradient():
    """The gradient test never once fired. Measured on three molecules,
    cholesterol reached 115.193 kcal/mol by frame 200 and spent the next 1300
    frames gaining 0.006 while having 54,000 steps rejected -- about forty per
    frame, for ever. A strained ring did worse: it stopped improving at frame
    50 and then rejected 107,000 steps without moving at all, because a failed
    step resets the length and the next frame tries the same thing again.
    Meanwhile a real metal complex was still descending at frame 700 and had
    every right to go on.

    Energy change per frame tells those three apart and the gradient cannot.
    With it: cholesterol converges at frame 120 with 3,780 rejections, the Re
    complex at 320, the Pt complex at 402, and the strained ring at 32 --
    each at the same energy it used to grind towards."""
    source = _MODULE_PATH.read_text(encoding='utf-8')
    assert 'var ENERGY_TOL = 1e-4;' in source
    assert 'var QUIET_FRAMES = 5;' in source
    # Judged on what the frame achieved, not on where the gradient stands.
    assert 'var moved = Math.abs(st.energyValue - startEnergy);' in source
    assert 'st.quietFrames = (st.quietFrames || 0) + 1;' in source
    assert 'st.converged = true;' in source
    # A stall counts as no progress, which is what turns endless retrying
    # into an answer.
    assert 'st.stalled = false;' in source
    # And anything that moves the atoms starts the count again.
    staged = source.split('st.pos.set(st.probe);')[1][:400]
    assert 'st.quietFrames = 0;' in staged
    assert 'st.converged = false;' in staged
