"""Browser-side force-field engine for live molecule manipulation.

The Submit tab lets the user grab an atom and drag it while the rest of the
molecule relaxes underneath the cursor.  The round trip from the browser into
the Python kernel and back is ~45 ms and collapses to ~13 Hz under a continuous
drag, so a per-frame relaxation cannot live in Python.  This module therefore
ships the whole relaxation loop as JavaScript: Python assigns force-field
parameters once when manipulate mode is entered, hands them over as a single
JSON payload, and is not touched again until the mouse is released.

The engine follows Avogadro's design, including its constants:

* frozen atoms are a gradient mask (Avogadro 2 ``frozenAtomMask``/``setMask``),
  never a position reset after the step -- resetting fights the integrator;
* steepest descent, not conjugate gradients, because CG's line search costs
  roughly eight gradient evaluations per iteration and Avogadro documents
  steepest descent as "the most fluid and interactive" choice;
* the batch size adapts to a wall-clock budget (Avogadro 2
  ``adaptChunkIterations``): ``chunkIterations`` starts at 5, targets 33 ms,
  is smoothed with 0.7 and clamped to [1, 200].

The budget is a *full frame* budget: the caller passes the measured time of the
previous frame, renderer included, into :js:func:`step`, because 3Dmol's own
geometry rebuild costs 0.8-12.1 ms per frame in the default stick+sphere style
and has to be paid out of the same 33 ms.

Payload contract (produced on the Python side, one JSON blob)::

    {
      "ok": true,
      "version": 1,
      "source": "rdkit-uff" | "openbabel-uff" | "geometric-fallback",
      "n_atoms": 104,
      "elements": ["C", "H", ...],                     # informational
      "bonds":    [{"i": 0, "j": 1, "k": 699.59, "r0": 1.514}, ...],
      "angles":   [{"i": 0, "j": 1, "k": 2, "kt": 290.03, "theta0": 109.47}, ...],
      "torsions": [{"i": 0, "j": 1, "k": 2, "l": 3, "v": 0.1177, "n": 3,
                    "phi0": 180.0}, ...],
      "vdw":      [{"x": 3.851, "d": 0.105}, ...],     # one entry per atom
      "warnings": []
    }

Units are the ones rdkit's ``rdForceFieldHelpers`` hands out, unchanged:
lengths in angstrom, energies in kcal/mol, ``k`` in kcal/mol/A^2, ``kt`` in
kcal/mol/rad^2, and **``theta0`` and ``phi0`` in degrees**.  ``vdw`` holds the
per-atom UFF self terms; pairs are formed in the browser with the geometric
mean combining rule (verified to be exactly what rdkit's pair call returns),
which keeps the payload at N entries instead of N^2 -- 10 KiB rather than
2.85 MiB at 392 atoms.

Functional form of each term:

* bond     ``E = 0.5 * k * (r - r0)^2``
* angle    ``E = 0.5 * kt * (theta - theta0)^2`` -- deliberately the harmonic
  surrogate, not UFF's cosine expansion. It is correct to second order at
  every ``theta0`` and the two forms agree to 0.011 kcal/mol over a whole
  molecule at the geometries a drag visits, while the exact n=1 (linear)
  form is the one that inverts its own well if the ``cos(n*theta0)`` factor
  is dropped. ``angles[].n`` is therefore informational here.
* torsion  ``E = v * (1 - cos(n*phi0) * cos(n*phi))``  (UFF's own form; ``v``
  arrives ready to use -- the exporter has already folded in UFF's 1/2 and
  the division by the number of torsions about the central bond, so applying
  a factor of 1/2 here would halve every barrier)
* vdW      ``E = D_ij * ((x_ij/r)^12 - 2*(x_ij/r)^6)``, shifted to zero at the
  cutoff, with ``x_ij = sqrt(x_i*x_j)`` and ``D_ij = sqrt(D_i*D_j)``

All four gradients are analytic; nothing in the frame loop differentiates
numerically.  ``tests/test_molecule_forcefield_js.py`` checks every one of them
against central finite differences.

A note for whoever wires this into the Submit tab: the transition-metal hazard
lives on the Python side, not here.  rdkit's UFF silently drops transition
metals from every bonded and non-bonded term, so a payload for anything that
``contains_metal()`` must come from Open Babel (see
``_optimize_xyz_openbabel`` in ``delfin/smiles_converter.py``).  This engine
relaxes exactly the terms it is given and cannot notice that a metal is
missing -- it only reports the atom as untermed in ``stats()``.
"""

# Avogadro 2 forcefield.cpp constants, mirrored here so the Python side and the
# tests can assert against the same numbers the JavaScript uses.
FF_INITIAL_CHUNK = 5
FF_TARGET_MS = 33.0
FF_SMOOTHING = 0.7
FF_MIN_CHUNK = 1
FF_MAX_CHUNK = 200

import json

# Payload contract version, mirrored from molecule_forcefield.PAYLOAD_VERSION.
FF_PAYLOAD_VERSION = 1

# Non-bonded bookkeeping.
FF_VDW_CUTOFF = 10.0
FF_VDW_SKIN = 2.0
FF_PAIR_REBUILD_FRAMES = 20

# Integrator safety rails.
FF_MAX_DISPLACEMENT = 0.10


MOLECULE_FF_BOOTSTRAP_JS = r"""
(function() {
    if (window.__delfinFFReady) return;
    window.__delfinFFReady = true;

    window._delfinFFByScope = window._delfinFFByScope || {};

    // --- Avogadro 2 (forcefield.cpp) batch-adaptation constants ---
    var INITIAL_CHUNK = 5;
    var TARGET_MS = 33.0;
    var SMOOTHING = 0.7;
    var MIN_CHUNK = 1;
    var MAX_CHUNK = 200;

    // --- non-bonded bookkeeping ---
    var VDW_CUTOFF = 10.0;      // A
    var VDW_SKIN = 2.0;         // A, Verlet skin on the neighbour list
    var PAIR_REBUILD_FRAMES = 20;

    // --- integrator safety rails ---
    var MAX_DISPLACEMENT = 0.10;   // A, hardest cap on one atom in one step
    var INITIAL_LAMBDA = 1e-4;     // A per (kcal/mol/A)
    var MAX_LAMBDA = 0.05;
    var MIN_LAMBDA = 1e-12;
    var GROW = 1.2;
    var SHRINK = 0.5;
    var GRAD_TOL = 1e-8;
    // Coincident atoms have no direction of their own. Evaluating the term at
    // this floor along a fixed axis keeps the restoring force finite and large;
    // skipping the term instead left a zero-length bond that nothing in the
    // field could ever pull apart again.
    var MIN_SEPARATION = 1e-4;     // A
    var MIN_SEPARATION2 = 1e-8;

    // Must equal molecule_forcefield.PAYLOAD_VERSION. The two sides carry
    // matching conventions for the torsion factor and the angle form; a
    // mismatch is refused rather than relaxed wrongly.
    var PAYLOAD_VERSION = 1;

    var DEG2RAD = Math.PI / 180.0;

    // Term selectors, so the drag loop can run everything while the tests and
    // the diagnostics panel can isolate one term at a time.
    var T_BOND = 1, T_ANGLE = 2, T_TORSION = 4, T_VDW = 8, T_RESTRAINT = 16;
    var T_INVERSION = 32;
    //: The hand, as a force rather than as a position.
    //
    // Everyone who does interactive manipulation well applies a force.
    // VMD/NAMD's IMD translates a drag into forces on the atoms; Narupa
    // biases the dynamics through an external field whose strength the user
    // controls; Reiher's real-time quantum chemistry hands the response back
    // as force feedback.  Nobody sets a coordinate, and that is why none of
    // them has to prevent unphysical geometry -- it cannot arise.  You pull,
    // the bond pulls back, and the atom follows exactly as far as the balance
    // allows, so tearing one means pulling harder than it holds.
    var T_TETHER = 64;
    var T_ALL = 127;
    // Restraint kinds, as sent in the payload.
    var R_DISTANCE = 0, R_ANGLE = 1, R_DIHEDRAL = 2;

    var CLOCK = (typeof performance !== 'undefined' && performance &&
                 typeof performance.now === 'function') ? performance : null;
    function now() { return CLOCK ? CLOCK.now() : Date.now(); }

    function isNum(v) { return typeof v === 'number' && isFinite(v); }

    // ------------------------------------------------------------------
    // payload -> typed arrays
    // ------------------------------------------------------------------
    function growI32(arr, need) {
        if (arr && arr.length >= need) return arr;
        var cap = arr ? arr.length : 64;
        while (cap < need) cap *= 2;
        var next = new Int32Array(cap);
        if (arr) next.set(arr);
        return next;
    }
    function growF64(arr, need) {
        if (arr && arr.length >= need) return arr;
        var cap = arr ? arr.length : 64;
        while (cap < need) cap *= 2;
        var next = new Float64Array(cap);
        if (arr) next.set(arr);
        return next;
    }

    function buildExclusions(st, bonds) {
        // 1-2 and 1-3 pairs, stored as a flat per-atom adjacency so the pair
        // build can scan a handful of ints instead of hashing.
        var n = st.n, i, j;
        var deg = new Int32Array(n);
        var nb = st.nBonds;
        for (i = 0; i < nb; i++) { deg[st.bondI[i]]++; deg[st.bondJ[i]]++; }
        var start = new Int32Array(n + 1);
        for (i = 0; i < n; i++) start[i + 1] = start[i] + deg[i];
        var adj = new Int32Array(start[n]);
        var fill = new Int32Array(n);
        for (i = 0; i < nb; i++) {
            var a = st.bondI[i], b = st.bondJ[i];
            adj[start[a] + fill[a]++] = b;
            adj[start[b] + fill[b]++] = a;
        }
        st.adjStart = start;
        st.adjList = adj;

        // Collect exclusions per atom: bonded neighbours plus their neighbours.
        var lists = new Array(n);
        for (i = 0; i < n; i++) lists[i] = [];
        for (i = 0; i < nb; i++) {
            lists[st.bondI[i]].push(st.bondJ[i]);
            lists[st.bondJ[i]].push(st.bondI[i]);
        }
        for (i = 0; i < n; i++) {
            var s = start[i], e = start[i + 1];
            for (j = s; j < e; j++) {
                var mid = adj[j];
                var ms = start[mid], me = start[mid + 1];
                for (var q = ms; q < me; q++) {
                    var far = adj[q];
                    if (far !== i) lists[i].push(far);
                }
            }
        }
        var exStart = new Int32Array(n + 1);
        var total = 0;
        for (i = 0; i < n; i++) {
            lists[i].sort(function(x, y) { return x - y; });
            var uniq = [];
            for (j = 0; j < lists[i].length; j++) {
                if (j === 0 || lists[i][j] !== lists[i][j - 1]) uniq.push(lists[i][j]);
            }
            lists[i] = uniq;
            total += uniq.length;
            exStart[i + 1] = total;
        }
        var exList = new Int32Array(total);
        var at = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < lists[i].length; j++) exList[at++] = lists[i][j];
        }
        st.exStart = exStart;
        st.exList = exList;
    }

    function isExcluded(st, i, j) {
        var s = st.exStart[i], e = st.exStart[i + 1];
        for (var q = s; q < e; q++) {
            if (st.exList[q] === j) return true;
        }
        return false;
    }

    function buildState(scopeKey, terms) {
        var warnings = [];
        if (!terms || typeof terms !== 'object') {
            return {error: 'no force-field payload'};
        }
        if (terms.ok === false) {
            return {error: 'payload reports ok=false'};
        }
        // Refuse a payload built against a different contract. Silently loading
        // one is how the torsion factor diverged between the two sides without
        // anything noticing: the molecule just relaxed wrongly.
        if (terms.version !== undefined && (terms.version | 0) !== PAYLOAD_VERSION) {
            return {error: 'force-field payload version ' + terms.version +
                           ', engine expects ' + PAYLOAD_VERSION};
        }
        var n = terms.n_atoms | 0;
        if (!(n > 0)) return {error: 'n_atoms must be positive'};
        if (terms.warnings && terms.warnings.length) {
            for (var w = 0; w < terms.warnings.length; w++) {
                warnings.push(String(terms.warnings[w]));
            }
        }

        var st = {
            scopeKey: scopeKey,
            n: n,
            source: terms.source ? String(terms.source) : 'unknown',
            collapsed: false,
            elements: terms.elements || [],
            warnings: warnings,
            chunk: INITIAL_CHUNK,
            lambda: INITIAL_LAMBDA,
            frozen: [],
            tetN: 0,
            tetA: new Int32Array(0), tetK: new Float64Array(0),
            tetX: new Float64Array(0), tetY: new Float64Array(0),
            tetZ: new Float64Array(0),
            frames: 0,
            totalSteps: 0,
            lastSteps: 0,
            lastFrameMs: 0,
            lastFfMs: 0,
            msPerStep: 0,
            rejects: 0,
            rollbacks: 0,
            pairRebuilds: 0,
            sinceRebuild: PAIR_REBUILD_FRAMES,
            stalled: false,
            converged: false,
            haveEnergy: false,
            energyValue: NaN,
            lastError: null,
            cutoff: VDW_CUTOFF,
            skin: VDW_SKIN
        };

        // -- bonds -------------------------------------------------------
        var src = terms.bonds || [];
        var nb = 0;
        var bi = new Int32Array(src.length), bj = new Int32Array(src.length);
        var bk = new Float64Array(src.length), br = new Float64Array(src.length);
        var t, i, j, k, l;
        for (t = 0; t < src.length; t++) {
            var b = src[t];
            i = b.i | 0; j = b.j | 0;
            if (i < 0 || j < 0 || i >= n || j >= n || i === j) { warnings.push('bad bond index'); continue; }
            if (!isNum(b.k) || b.k < 0 || !isNum(b.r0) || b.r0 <= 0) {
                // A negative force constant is unbounded below and the
                // minimiser will happily descend into it forever.
                warnings.push('bad bond params'); continue;
            }
            bi[nb] = i; bj[nb] = j; bk[nb] = b.k; br[nb] = b.r0; nb++;
        }
        st.nBonds = nb; st.bondI = bi; st.bondJ = bj; st.bondK = bk; st.bondR0 = br;

        // -- angles ------------------------------------------------------
        src = terms.angles || [];
        var na = 0;
        var aI = new Int32Array(src.length), aJ = new Int32Array(src.length),
            aK = new Int32Array(src.length);
        var aKt = new Float64Array(src.length), aT0 = new Float64Array(src.length);
        for (t = 0; t < src.length; t++) {
            var g = src[t];
            i = g.i | 0; j = g.j | 0; k = g.k | 0;
            if (i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n ||
                i === j || j === k || i === k) { warnings.push('bad angle index'); continue; }
            if (!isNum(g.kt) || g.kt < 0 || !isNum(g.theta0)) {
                warnings.push('bad angle params'); continue;
            }
            aI[na] = i; aJ[na] = j; aK[na] = k;
            aKt[na] = g.kt; aT0[na] = g.theta0 * DEG2RAD; na++;
        }
        st.nAngles = na; st.angI = aI; st.angJ = aJ; st.angK = aK;
        st.angKt = aKt; st.angT0 = aT0;

        // -- torsions ----------------------------------------------------
        src = terms.torsions || [];
        var nt = 0;
        var tI = new Int32Array(src.length), tJ = new Int32Array(src.length),
            tK = new Int32Array(src.length), tL = new Int32Array(src.length);
        var tV = new Float64Array(src.length), tN = new Float64Array(src.length),
            tC = new Float64Array(src.length);
        for (t = 0; t < src.length; t++) {
            var d = src[t];
            i = d.i | 0; j = d.j | 0; k = d.k | 0; l = d.l | 0;
            if (i < 0 || j < 0 || k < 0 || l < 0 ||
                i >= n || j >= n || k >= n || l >= n ||
                i === j || j === k || k === l ||
                i === k || i === l || j === l) { warnings.push('bad torsion index'); continue; }
            if (!isNum(d.v)) { warnings.push('bad torsion params'); continue; }
            var per = isNum(d.n) ? d.n : 3;
            var ph = isNum(d.phi0) ? d.phi0 : 180.0;
            tI[nt] = i; tJ[nt] = j; tK[nt] = k; tL[nt] = l;
            tV[nt] = d.v; tN[nt] = per;
            // UFF folds the phase into cos(n*phi0), which is +-1 for every
            // phase the rules produce; caching it keeps the frame loop clean.
            tC[nt] = Math.cos(per * ph * DEG2RAD);
            nt++;
        }
        st.nTors = nt; st.torI = tI; st.torJ = tJ; st.torK = tK; st.torL = tL;
        st.torV = tV; st.torN = tN; st.torCos0 = tC;

        // -- per-atom vdW ------------------------------------------------
        src = terms.vdw || [];
        var vx = new Float64Array(n), vd = new Float64Array(n);
        var untermed = 0;
        for (i = 0; i < n; i++) {
            var v = src[i];
            if (v && isNum(v.x) && isNum(v.d) && v.x > 0 && v.d >= 0) {
                vx[i] = v.x; vd[i] = v.d;
            } else {
                vx[i] = 0; vd[i] = 0; untermed++;
            }
        }
        if (untermed) warnings.push(untermed + ' atom(s) without vdW parameters');
        st.vdwX = vx; st.vdwD = vd;
        st.untermed = untermed;

        // -- working buffers --------------------------------------------
        var n3 = 3 * n;
        st.pos = new Float64Array(n3);
        st.grad = new Float64Array(n3);
        st.trial = new Float64Array(n3);
        st.tgrad = new Float64Array(n3);
        st.dir = new Float64Array(n3);
        st.prevGrad = new Float64Array(n3);
        st.havePrevGrad = false;
        st.probe = new Float64Array(n3);
        st.pgrad = new Float64Array(n3);
        st.lastGood = new Float64Array(n3);
        st.pairRef = new Float64Array(n3);
        st.mask = new Float64Array(n3);
        for (i = 0; i < n3; i++) st.mask[i] = 1.0;
        st.havePos = false;

        st.nPairs = 0;
        st.pairI = null; st.pairJ = null;
        st.pairX = null; st.pairD = null; st.pairShift = null;

        // -- inversions --------------------------------------------------
        // UFF's out-of-plane term. Without it, planarity at a three-coordinate
        // centre came only from three angles wanting 120 degrees, and that is
        // a quartic penalty, not a quadratic one: lifting a carbon one degree
        // out of its own plane cost 0.002 kcal/mol, which is nothing. It also
        // did not reach a donor at a metal at all, because those angles are
        // not typed. One entry per apex, three per centre.
        src = terms.inversions || [];
        var ni = 0;
        var invI = new Int32Array(src.length), invJ = new Int32Array(src.length),
            invK = new Int32Array(src.length), invL = new Int32Array(src.length);
        var invW = new Float64Array(src.length);
        for (t = 0; t < src.length; t++) {
            var iv = src[t];
            var ii = iv.i | 0, ij = iv.j | 0, ik = iv.k | 0, il = iv.l | 0;
            if (ii < 0 || ii >= n || ij < 0 || ij >= n ||
                ik < 0 || ik >= n || il < 0 || il >= n) {
                warnings.push('inversion index out of range'); continue;
            }
            if (!isNum(iv.k_inv) || iv.k_inv < 0) continue;
            invI[ni] = ii; invJ[ni] = ij; invK[ni] = ik; invL[ni] = il;
            invW[ni] = iv.k_inv; ni++;
        }

        // -- restraints --------------------------------------------------
        // Values the user asked to hold while everything else relaxes. Kept as
        // ordinary harmonic terms, so a set of restraints that cannot all be
        // satisfied at once settles at a compromise instead of tearing the
        // structure apart.
        src = terms.restraints || [];
        var nr = 0;
        var rKind = new Int32Array(src.length);
        var rA = new Int32Array(src.length), rB = new Int32Array(src.length),
            rC = new Int32Array(src.length), rD = new Int32Array(src.length);
        var rK = new Float64Array(src.length), rT = new Float64Array(src.length);
        for (t = 0; t < src.length; t++) {
            var r = src[t];
            var kindName = String(r.kind || '');
            var kind = kindName === 'angle' ? R_ANGLE
                     : kindName === 'dihedral' ? R_DIHEDRAL
                     : kindName === 'distance' ? R_DISTANCE : -1;
            if (kind < 0) { warnings.push('unknown restraint kind'); continue; }
            var need = kind === R_DISTANCE ? 2 : (kind === R_ANGLE ? 3 : 4);
            var list = r.atoms || [];
            if (list.length !== need) { warnings.push('bad restraint atoms'); continue; }
            var ok = true, idx = [];
            for (var q = 0; q < need; q++) {
                var one = list[q] | 0;
                if (one < 0 || one >= n || idx.indexOf(one) >= 0) { ok = false; break; }
                idx.push(one);
            }
            if (!ok) { warnings.push('bad restraint atoms'); continue; }
            if (!isNum(r.value) || !isNum(r.k) || r.k < 0) {
                warnings.push('bad restraint params'); continue;
            }
            if (kind === R_DISTANCE && r.value <= 0) {
                warnings.push('bad restraint params'); continue;
            }
            rKind[nr] = kind;
            rA[nr] = idx[0]; rB[nr] = idx[1];
            rC[nr] = need > 2 ? idx[2] : 0;
            rD[nr] = need > 3 ? idx[3] : 0;
            rK[nr] = r.k;
            rT[nr] = kind === R_DISTANCE ? r.value : r.value * DEG2RAD;
            nr++;
        }
        st.nInversions = ni;
        st.invI = invI; st.invJ = invJ; st.invK = invK; st.invL = invL;
        st.invW = invW;
        st.nRestraints = nr;
        st.resKind = rKind; st.resA = rA; st.resB = rB; st.resC = rC; st.resD = rD;
        st.resK = rK; st.resTarget = rT;

        buildExclusions(st, terms.bonds || []);
        return {state: st};
    }

    // ------------------------------------------------------------------
    // neighbour pair list
    // ------------------------------------------------------------------
    function rebuildPairs(st) {
        var n = st.n, pos = st.pos;
        var listCut = st.cutoff + st.skin;
        var listCut2 = listCut * listCut;
        var count = 0, i, j;
        var pi = st.pairI, pj = st.pairJ, px = st.pairX, pd = st.pairD, ps = st.pairShift;
        var cutInv = 1.0 / st.cutoff;
        for (i = 0; i < n; i++) {
            if (st.vdwX[i] <= 0) continue;
            var i3 = 3 * i;
            var xi = pos[i3], yi = pos[i3 + 1], zi = pos[i3 + 2];
            for (j = i + 1; j < n; j++) {
                if (st.vdwX[j] <= 0) continue;
                var j3 = 3 * j;
                var dx = pos[j3] - xi, dy = pos[j3 + 1] - yi, dz = pos[j3 + 2] - zi;
                var r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > listCut2) continue;
                if (isExcluded(st, i, j)) continue;
                pi = growI32(pi, count + 1); pj = growI32(pj, count + 1);
                px = growF64(px, count + 1); pd = growF64(pd, count + 1);
                ps = growF64(ps, count + 1);
                // Geometric-mean combining rule, verified against rdkit's own
                // pairwise GetUFFVdWParams.
                var xij = Math.sqrt(st.vdwX[i] * st.vdwX[j]);
                var dij = Math.sqrt(st.vdwD[i] * st.vdwD[j]);
                pi[count] = i; pj[count] = j;
                px[count] = xij; pd[count] = dij;
                // Shift the potential to zero at the cutoff so the energy stays
                // continuous when a pair enters or leaves the list.
                var tc2 = xij * xij * cutInv * cutInv;
                var tc6 = tc2 * tc2 * tc2;
                ps[count] = dij * (tc6 * tc6 - 2.0 * tc6);
                count++;
            }
        }
        st.pairI = pi; st.pairJ = pj; st.pairX = px; st.pairD = pd; st.pairShift = ps;
        st.nPairs = count;
        st.pairRef.set(pos);
        st.sinceRebuild = 0;
        st.pairRebuilds++;
        st.haveEnergy = false;
    }

    function pairsAreStale(st) {
        if (!st.pairI) return true;
        if (st.sinceRebuild >= PAIR_REBUILD_FRAMES) return true;
        // Verlet criterion: half the skin is the most any atom may drift
        // before a pair could have crossed the cutoff unseen.
        var half = 0.5 * st.skin, n3 = 3 * st.n;
        for (var i = 0; i < n3; i++) {
            var d = st.pos[i] - st.pairRef[i];
            if (d > half || d < -half) return true;
        }
        return false;
    }

    // ------------------------------------------------------------------
    // energy + analytic gradient
    // ------------------------------------------------------------------
    function evaluate(st, pos, grad, which) {
        var n3 = 3 * st.n, i;
        var e = 0.0;
        if (grad) { for (i = 0; i < n3; i++) grad[i] = 0.0; }
        st.collapsed = false;

        // -- harmonic bond stretch --------------------------------------
        if (which & T_BOND) {
            var nb = st.nBonds;
            for (i = 0; i < nb; i++) {
                var a = 3 * st.bondI[i], b = 3 * st.bondJ[i];
                var dx = pos[b] - pos[a], dy = pos[b + 1] - pos[a + 1],
                    dz = pos[b + 2] - pos[a + 2];
                var r2 = dx * dx + dy * dy + dz * dz;
                var r;
                if (r2 < MIN_SEPARATION2) {
                    // Coincident bonded atoms: substitute a fixed direction so
                    // the bond still pushes them apart. Dropping the term here
                    // froze the molecule permanently, because the 1-2 vdW pair
                    // is excluded by construction and no angle term survives a
                    // zero-length leg either.
                    st.collapsed = true;
                    r = MIN_SEPARATION;
                    dx = MIN_SEPARATION; dy = 0.0; dz = 0.0;
                } else {
                    r = Math.sqrt(r2);
                }
                var dr = r - st.bondR0[i];
                e += 0.5 * st.bondK[i] * dr * dr;
                if (grad) {
                    var f = st.bondK[i] * dr / r;
                    grad[a] -= f * dx; grad[a + 1] -= f * dy; grad[a + 2] -= f * dz;
                    grad[b] += f * dx; grad[b + 1] += f * dy; grad[b + 2] += f * dz;
                }
            }
        }

        // -- harmonic angle bend ----------------------------------------
        if (which & T_ANGLE) {
            var na = st.nAngles;
            for (i = 0; i < na; i++) {
                var ia = 3 * st.angI[i], ja = 3 * st.angJ[i], ka = 3 * st.angK[i];
                var ux = pos[ia] - pos[ja], uy = pos[ia + 1] - pos[ja + 1],
                    uz = pos[ia + 2] - pos[ja + 2];
                var vx = pos[ka] - pos[ja], vy = pos[ka + 1] - pos[ja + 1],
                    vz = pos[ka + 2] - pos[ja + 2];
                var ru2 = ux * ux + uy * uy + uz * uz;
                var rv2 = vx * vx + vy * vy + vz * vz;
                if (ru2 < 1e-16 || rv2 < 1e-16) continue;
                var ru = Math.sqrt(ru2), rv = Math.sqrt(rv2);
                var c = (ux * vx + uy * vy + uz * vz) / (ru * rv);
                if (c > 1.0) c = 1.0; else if (c < -1.0) c = -1.0;
                var theta = Math.acos(c);
                var dth = theta - st.angT0[i];
                e += 0.5 * st.angKt[i] * dth * dth;
                if (grad) {
                    var s = Math.sqrt(1.0 - c * c);
                    if (s < 1e-8) continue;   // collinear: dtheta/dx is singular
                    var pref = -st.angKt[i] * dth / s;
                    var inv = 1.0 / (ru * rv);
                    var gix = pref * (vx * inv - c * ux / ru2);
                    var giy = pref * (vy * inv - c * uy / ru2);
                    var giz = pref * (vz * inv - c * uz / ru2);
                    var gkx = pref * (ux * inv - c * vx / rv2);
                    var gky = pref * (uy * inv - c * vy / rv2);
                    var gkz = pref * (uz * inv - c * vz / rv2);
                    grad[ia] += gix; grad[ia + 1] += giy; grad[ia + 2] += giz;
                    grad[ka] += gkx; grad[ka + 1] += gky; grad[ka + 2] += gkz;
                    grad[ja] -= gix + gkx; grad[ja + 1] -= giy + gky;
                    grad[ja + 2] -= giz + gkz;
                }
            }
        }

        // -- cosine torsion ---------------------------------------------
        if (which & T_TORSION) {
            var nt = st.nTors;
            for (i = 0; i < nt; i++) {
                var i0 = 3 * st.torI[i], j0 = 3 * st.torJ[i],
                    k0 = 3 * st.torK[i], l0 = 3 * st.torL[i];
                var b1x = pos[j0] - pos[i0], b1y = pos[j0 + 1] - pos[i0 + 1],
                    b1z = pos[j0 + 2] - pos[i0 + 2];
                var b2x = pos[k0] - pos[j0], b2y = pos[k0 + 1] - pos[j0 + 1],
                    b2z = pos[k0 + 2] - pos[j0 + 2];
                var b3x = pos[l0] - pos[k0], b3y = pos[l0 + 1] - pos[k0 + 1],
                    b3z = pos[l0 + 2] - pos[k0 + 2];
                var ax = b1y * b2z - b1z * b2y,
                    ay = b1z * b2x - b1x * b2z,
                    az = b1x * b2y - b1y * b2x;
                var bx = b2y * b3z - b2z * b3y,
                    by = b2z * b3x - b2x * b3z,
                    bz = b2x * b3y - b2y * b3x;
                var a2 = ax * ax + ay * ay + az * az;
                var bb2 = bx * bx + by * by + bz * bz;
                var c2 = b2x * b2x + b2y * b2y + b2z * b2z;
                // Degenerate: three collinear atoms leave the dihedral plane
                // undefined, so the term simply does not contribute.
                if (a2 < 1e-12 || bb2 < 1e-12 || c2 < 1e-12) continue;
                var rb2 = Math.sqrt(c2);
                var xdot = ax * bx + ay * by + az * bz;
                var ydot = rb2 * (b1x * bx + b1y * by + b1z * bz);
                var phi = Math.atan2(ydot, xdot);
                var per = st.torN[i], cos0 = st.torCos0[i], V = st.torV[i];
                // E = v * (1 - cos(n*phi0) * cos(n*phi)). No factor of 1/2 here:
                // the exporter already folds UFF's 1/2 and the division by the
                // number of torsions about the central bond into v. Applying it
                // twice put every torsion barrier at half strength.
                e += V * (1.0 - cos0 * Math.cos(per * phi));
                if (grad) {
                    var dEdphi = V * per * cos0 * Math.sin(per * phi);
                    var fi = -rb2 / a2, fl = rb2 / bb2;
                    var dix = fi * ax, diy = fi * ay, diz = fi * az;
                    var dlx = fl * bx, dly = fl * by, dlz = fl * bz;
                    // The two central atoms take the remainder; this is the
                    // standard redistribution and the four gradients sum to
                    // zero, as a translation-invariant coordinate demands.
                    var p = (b1x * b2x + b1y * b2y + b1z * b2z) / c2;
                    var q = (b3x * b2x + b3y * b2y + b3z * b2z) / c2;
                    var djx = -(1.0 + p) * dix + q * dlx;
                    var djy = -(1.0 + p) * diy + q * dly;
                    var djz = -(1.0 + p) * diz + q * dlz;
                    var dkx = p * dix - (1.0 + q) * dlx;
                    var dky = p * diy - (1.0 + q) * dly;
                    var dkz = p * diz - (1.0 + q) * dlz;
                    grad[i0] += dEdphi * dix; grad[i0 + 1] += dEdphi * diy;
                    grad[i0 + 2] += dEdphi * diz;
                    grad[j0] += dEdphi * djx; grad[j0 + 1] += dEdphi * djy;
                    grad[j0 + 2] += dEdphi * djz;
                    grad[k0] += dEdphi * dkx; grad[k0 + 1] += dEdphi * dky;
                    grad[k0 + 2] += dEdphi * dkz;
                    grad[l0] += dEdphi * dlx; grad[l0 + 1] += dEdphi * dly;
                    grad[l0 + 2] += dEdphi * dlz;
                }
            }
        }

        // -- user restraints --------------------------------------------
        // A distance, an angle or a dihedral the user asked to hold. Harmonic
        // in the coordinate itself -- unlike the cosine torsion above, which
        // has its minima fixed by the periodicity and cannot hold an arbitrary
        // dihedral.
        if (which & T_INVERSION) {
            // E = k (1 - cos w), w the angle between the bond to the apex and
            // the plane of the other two. Zero and flat at planarity, and
            // quadratic in w just off it, which is the whole point: the angle
            // terms alone are quartic there and hold almost nothing.
            var nInv = st.nInversions;
            for (i = 0; i < nInv; i++) {
                var oc = 3 * st.invI[i], oj = 3 * st.invJ[i],
                    ok2 = 3 * st.invK[i], ol = 3 * st.invL[i];
                var okw = st.invW[i];
                // Closed form for the energy; the gradient is a central
                // difference over the twelve components involved. There are
                // only three entries per planar centre, so the cost is small
                // and the derivative cannot be got subtly wrong by hand.
                var oidx = [oc, oj, ok2, ol];
                var invEnergy = function() {
                    var ax = pos[oj] - pos[oc], ay = pos[oj + 1] - pos[oc + 1],
                        az = pos[oj + 2] - pos[oc + 2];
                    var bx = pos[ok2] - pos[oc], by = pos[ok2 + 1] - pos[oc + 1],
                        bz = pos[ok2 + 2] - pos[oc + 2];
                    var cx = pos[ol] - pos[oc], cy = pos[ol + 1] - pos[oc + 1],
                        cz = pos[ol + 2] - pos[oc + 2];
                    var nx = ay * bz - az * by, ny = az * bx - ax * bz,
                        nz = ax * by - ay * bx;
                    var nn = Math.sqrt(nx * nx + ny * ny + nz * nz);
                    var cc = Math.sqrt(cx * cx + cy * cy + cz * cz);
                    if (nn < 1e-9 || cc < 1e-9) return 0.0;
                    var sn = (nx * cx + ny * cy + nz * cz) / (nn * cc);
                    if (sn > 1.0) sn = 1.0; else if (sn < -1.0) sn = -1.0;
                    return okw * (1.0 - Math.sqrt(1.0 - sn * sn));
                };
                e += invEnergy();
                if (grad) {
                    var h = 1e-5;
                    for (var oa = 0; oa < 4; oa++) {
                        for (var ob = 0; ob < 3; ob++) {
                            var slot = oidx[oa] + ob;
                            var keep = pos[slot];
                            pos[slot] = keep + h;
                            var eUp = invEnergy();
                            pos[slot] = keep - h;
                            var eDown = invEnergy();
                            pos[slot] = keep;
                            grad[slot] += (eUp - eDown) / (2 * h);
                        }
                    }
                }
            }
        }

        if (which & T_RESTRAINT) {
            var nres = st.nRestraints;
            for (i = 0; i < nres; i++) {
                var rk = st.resK[i], rt = st.resTarget[i];
                var ra = 3 * st.resA[i], rb = 3 * st.resB[i];
                var rc = 3 * st.resC[i], rd = 3 * st.resD[i];
                var kind = st.resKind[i];

                if (kind === R_DISTANCE) {
                    var ex = pos[rb] - pos[ra], ey = pos[rb + 1] - pos[ra + 1],
                        ez = pos[rb + 2] - pos[ra + 2];
                    var er2 = ex * ex + ey * ey + ez * ez;
                    var er;
                    if (er2 < MIN_SEPARATION2) {
                        st.collapsed = true;
                        er = MIN_SEPARATION; ex = MIN_SEPARATION; ey = 0.0; ez = 0.0;
                    } else {
                        er = Math.sqrt(er2);
                    }
                    var edr = er - rt;
                    e += 0.5 * rk * edr * edr;
                    if (grad) {
                        var ef = rk * edr / er;
                        grad[ra] -= ef * ex; grad[ra + 1] -= ef * ey; grad[ra + 2] -= ef * ez;
                        grad[rb] += ef * ex; grad[rb + 1] += ef * ey; grad[rb + 2] += ef * ez;
                    }
                } else if (kind === R_ANGLE) {
                    var aux = pos[ra] - pos[rb], auy = pos[ra + 1] - pos[rb + 1],
                        auz = pos[ra + 2] - pos[rb + 2];
                    var avx = pos[rc] - pos[rb], avy = pos[rc + 1] - pos[rb + 1],
                        avz = pos[rc + 2] - pos[rb + 2];
                    var au2 = aux * aux + auy * auy + auz * auz;
                    var av2 = avx * avx + avy * avy + avz * avz;
                    if (au2 < 1e-16 || av2 < 1e-16) continue;
                    var au = Math.sqrt(au2), av = Math.sqrt(av2);
                    var ac = (aux * avx + auy * avy + auz * avz) / (au * av);
                    if (ac > 1.0) ac = 1.0; else if (ac < -1.0) ac = -1.0;
                    var ath = Math.acos(ac) - rt;
                    e += 0.5 * rk * ath * ath;
                    if (grad) {
                        var asn = Math.sqrt(1.0 - ac * ac);
                        if (asn < 1e-8) continue;
                        var apref = -rk * ath / asn;
                        var ainv = 1.0 / (au * av);
                        var agix = apref * (avx * ainv - ac * aux / au2);
                        var agiy = apref * (avy * ainv - ac * auy / au2);
                        var agiz = apref * (avz * ainv - ac * auz / au2);
                        var agkx = apref * (aux * ainv - ac * avx / av2);
                        var agky = apref * (auy * ainv - ac * avy / av2);
                        var agkz = apref * (auz * ainv - ac * avz / av2);
                        grad[ra] += agix; grad[ra + 1] += agiy; grad[ra + 2] += agiz;
                        grad[rc] += agkx; grad[rc + 1] += agky; grad[rc + 2] += agkz;
                        grad[rb] -= agix + agkx; grad[rb + 1] -= agiy + agky;
                        grad[rb + 2] -= agiz + agkz;
                    }
                } else {
                    var tb1x = pos[rb] - pos[ra], tb1y = pos[rb + 1] - pos[ra + 1],
                        tb1z = pos[rb + 2] - pos[ra + 2];
                    var tb2x = pos[rc] - pos[rb], tb2y = pos[rc + 1] - pos[rb + 1],
                        tb2z = pos[rc + 2] - pos[rb + 2];
                    var tb3x = pos[rd] - pos[rc], tb3y = pos[rd + 1] - pos[rc + 1],
                        tb3z = pos[rd + 2] - pos[rc + 2];
                    var tax = tb1y * tb2z - tb1z * tb2y,
                        tay = tb1z * tb2x - tb1x * tb2z,
                        taz = tb1x * tb2y - tb1y * tb2x;
                    var tbx = tb2y * tb3z - tb2z * tb3y,
                        tby = tb2z * tb3x - tb2x * tb3z,
                        tbz = tb2x * tb3y - tb2y * tb3x;
                    var ta2 = tax * tax + tay * tay + taz * taz;
                    var tbb2 = tbx * tbx + tby * tby + tbz * tbz;
                    var tc2 = tb2x * tb2x + tb2y * tb2y + tb2z * tb2z;
                    if (ta2 < 1e-12 || tbb2 < 1e-12 || tc2 < 1e-12) continue;
                    var trb2 = Math.sqrt(tc2);
                    var txd = tax * tbx + tay * tby + taz * tbz;
                    var tyd = trb2 * (tb1x * tbx + tb1y * tby + tb1z * tbz);
                    var tphi = Math.atan2(tyd, txd);
                    // Shortest way round: a dihedral 179 deg from its target is
                    // 1 deg away, not 359.
                    var dphi = tphi - rt;
                    while (dphi > Math.PI) dphi -= 2 * Math.PI;
                    while (dphi < -Math.PI) dphi += 2 * Math.PI;
                    e += 0.5 * rk * dphi * dphi;
                    if (grad) {
                        var tdE = rk * dphi;
                        var tfi = -trb2 / ta2, tfl = trb2 / tbb2;
                        var tdix = tfi * tax, tdiy = tfi * tay, tdiz = tfi * taz;
                        var tdlx = tfl * tbx, tdly = tfl * tby, tdlz = tfl * tbz;
                        var tp = (tb1x * tb2x + tb1y * tb2y + tb1z * tb2z) / tc2;
                        var tq = (tb3x * tb2x + tb3y * tb2y + tb3z * tb2z) / tc2;
                        var tdjx = -(1.0 + tp) * tdix + tq * tdlx;
                        var tdjy = -(1.0 + tp) * tdiy + tq * tdly;
                        var tdjz = -(1.0 + tp) * tdiz + tq * tdlz;
                        var tdkx = tp * tdix - (1.0 + tq) * tdlx;
                        var tdky = tp * tdiy - (1.0 + tq) * tdly;
                        var tdkz = tp * tdiz - (1.0 + tq) * tdlz;
                        grad[ra] += tdE * tdix; grad[ra + 1] += tdE * tdiy;
                        grad[ra + 2] += tdE * tdiz;
                        grad[rb] += tdE * tdjx; grad[rb + 1] += tdE * tdjy;
                        grad[rb + 2] += tdE * tdjz;
                        grad[rc] += tdE * tdkx; grad[rc + 1] += tdE * tdky;
                        grad[rc + 2] += tdE * tdkz;
                        grad[rd] += tdE * tdlx; grad[rd + 1] += tdE * tdly;
                        grad[rd + 2] += tdE * tdlz;
                    }
                }
            }
        }

        // -- Lennard-Jones 12-6 -----------------------------------------
        if (which & T_VDW) {
            var np = st.nPairs, cut2 = st.cutoff * st.cutoff;
            for (i = 0; i < np; i++) {
                var pa = 3 * st.pairI[i], pb = 3 * st.pairJ[i];
                var qx = pos[pb] - pos[pa], qy = pos[pb + 1] - pos[pa + 1],
                    qz = pos[pb + 2] - pos[pa + 2];
                var q2 = qx * qx + qy * qy + qz * qz;
                if (q2 > cut2) continue;
                if (q2 < MIN_SEPARATION2) {
                    // Exactly coincident: the separation vector is zero, so the
                    // gradient would be zero however large the energy grows.
                    // Give the pair a direction to be pushed apart along.
                    st.collapsed = true;
                    qx = MIN_SEPARATION; qy = 0.0; qz = 0.0;
                    q2 = MIN_SEPARATION2;
                }
                var xij = st.pairX[i], dij = st.pairD[i];
                // Below a quarter of the pair's vdW distance the two atoms are
                // fused and there is no chemistry left to model, only a need to
                // push them apart. The true 12-6 energy there reaches ~1e42,
                // which swamps every other term and blinds the line search's
                // energy comparison in double precision. Clamp the separation so
                // the repulsion stays enormous but representable.
                var qmin2 = 0.0625 * xij * xij;
                if (q2 < qmin2) q2 = qmin2;
                var t2 = xij * xij / q2;
                var t6 = t2 * t2 * t2;
                var t12 = t6 * t6;
                e += dij * (t12 - 2.0 * t6) - st.pairShift[i];
                if (grad) {
                    // (dE/dr)/r, so the direction vector needs no normalising
                    var fr = 12.0 * dij * (t6 - t12) / q2;
                    grad[pa] -= fr * qx; grad[pa + 1] -= fr * qy; grad[pa + 2] -= fr * qz;
                    grad[pb] += fr * qx; grad[pb + 1] += fr * qy; grad[pb + 2] += fr * qz;
                }
            }
        }
        if ((which & T_TETHER) && st.tetN) {
            var tg = grad ? new Float64Array(3 * st.n) : null;
            for (var tt = 0; tt < st.tetN; tt++) {
                var ta = 3 * st.tetA[tt], tk = st.tetK[tt];
                var tdx = pos[ta] - st.tetX[tt],
                    tdy = pos[ta + 1] - st.tetY[tt],
                    tdz = pos[ta + 2] - st.tetZ[tt];
                e += 0.5 * tk * (tdx * tdx + tdy * tdy + tdz * tdz);
                if (tg) {
                    tg[ta] += tk * tdx;
                    tg[ta + 1] += tk * tdy;
                    tg[ta + 2] += tk * tdz;
                }
            }
            if (tg) {
                // A pull that only deforms.
                //
                // Every other term is invariant under moving the molecule, so
                // its net force and torque are zero.  A spring is not: it has
                // a preferred place, and sliding the whole structure towards
                // it is a zero-energy mode of everything else.  So a minimiser
                // answers a pull by *translating*, never by stretching --
                // measured, at every force constant from 1 to 2000 the atom
                // reached its target in full and the bond it hangs on stayed
                // at 1.094 A.  Re-framing afterwards cannot help: a rigid move
                // changes no internal coordinate.
                //
                // Taking the net force and torque out of the pull leaves the
                // part that can only be answered by changing shape, which is
                // what the hand meant.  It is what inertia does in interactive
                // MD, where the rest of the molecule simply cannot respond in
                // one step; here there is no time, so it has to be said.
                projectOutRigid(st, pos, tg);
                for (var tq = 0; tq < 3 * st.n; tq++) grad[tq] += tg[tq];
            }
        }
        return e;
    }

    function projectOutRigid(st, pos, f) {
        // Six modes: three translations, three rotations about the centroid.
        // Gram-Schmidt, so no matrix is inverted and a linear molecule (whose
        // rotation modes are degenerate) costs nothing extra.
        var n = st.n, n3 = 3 * n, a, k, m;
        if (n < 2) return;
        var cx = 0, cy = 0, cz = 0;
        for (a = 0; a < n; a++) { cx += pos[3*a]; cy += pos[3*a+1]; cz += pos[3*a+2]; }
        cx /= n; cy /= n; cz /= n;
        var modes = [];
        for (k = 0; k < 3; k++) {
            var t = new Float64Array(n3);
            for (a = 0; a < n; a++) t[3*a + k] = 1.0;
            modes.push(t);
        }
        for (k = 0; k < 3; k++) {
            var r = new Float64Array(n3);
            for (a = 0; a < n; a++) {
                var qx = pos[3*a] - cx, qy = pos[3*a+1] - cy, qz = pos[3*a+2] - cz;
                if (k === 0) { r[3*a] = 0; r[3*a+1] = -qz; r[3*a+2] = qy; }
                else if (k === 1) { r[3*a] = qz; r[3*a+1] = 0; r[3*a+2] = -qx; }
                else { r[3*a] = -qy; r[3*a+1] = qx; r[3*a+2] = 0; }
            }
            modes.push(r);
        }
        for (m = 0; m < modes.length; m++) {
            var v = modes[m];
            for (var j = 0; j < m; j++) {
                var u = modes[j], dot = 0.0;
                for (k = 0; k < n3; k++) dot += v[k] * u[k];
                for (k = 0; k < n3; k++) v[k] -= dot * u[k];
            }
            var norm = 0.0;
            for (k = 0; k < n3; k++) norm += v[k] * v[k];
            norm = Math.sqrt(norm);
            if (norm < 1e-8) { modes[m] = null; continue; }
            for (k = 0; k < n3; k++) v[k] /= norm;
        }
        for (m = 0; m < modes.length; m++) {
            var w = modes[m];
            if (!w) continue;
            var along = 0.0;
            for (k = 0; k < n3; k++) along += f[k] * w[k];
            for (k = 0; k < n3; k++) f[k] -= along * w[k];
        }
    }

    function applyMask(st, grad) {
        if (!st.frozen.length) return;
        var n3 = 3 * st.n;
        for (var i = 0; i < n3; i++) grad[i] *= st.mask[i];
    }

    function evalMasked(st, pos, grad, which) {
        var e = evaluate(st, pos, grad, which);
        if (grad) applyMask(st, grad);
        return e;
    }

    function allFinite(arr) {
        for (var i = 0; i < arr.length; i++) {
            if (!isFinite(arr[i])) return false;
        }
        return true;
    }

    //: A frame that moves the energy less than this has stopped achieving
    //: anything, whatever the gradient still says.
    var ENERGY_TOL = 1e-4;          // kcal/mol per frame
    var QUIET_FRAMES = 5;           // how many in a row before we believe it

    var WOLFE_C1 = 1e-4;        // sufficient decrease
    var WOLFE_C2 = 0.4;         // curvature; loose, as fits an inexact search
    var LINE_TRIES = 8;
    //: Steps of plain descent before the conjugate memory is engaged. Which
    //: basin a local optimiser lands in depends on the path it takes, and a
    //: conjugate direction commits to one earlier than the gradient does. On
    //: a badly strained start that was the wrong one: a ring drawn shut
    //: settled 60 kcal/mol above where descent got to. Letting the worst of
    //: the strain out first costs a few frames and removes the only case
    //: where this was worse than what it replaces.
    var CG_WARMUP_STEPS = 200;

    // The direction: Polak-Ribiere, restarted whenever it would point uphill.
    // Steepest descent goes downhill, which is not the same as going towards
    // the minimum -- in a long narrow valley, and every stiff bond against
    // every soft torsion is one, it crosses back and forth and creeps along.
    function searchDirection(st, n3) {
        var i, beta = 0.0;
        st.cgWarmup = (st.cgWarmup || 0) + 1;
        if (st.cgWarmup <= CG_WARMUP_STEPS) st.havePrevGrad = false;
        if (st.havePrevGrad) {
            var num = 0.0, den = 0.0;
            for (i = 0; i < n3; i++) {
                num += st.grad[i] * (st.grad[i] - st.prevGrad[i]);
                den += st.prevGrad[i] * st.prevGrad[i];
            }
            if (den > 1e-16) beta = num / den;
            if (!isFinite(beta) || beta < 0) beta = 0.0;
        }
        var slope = 0.0;
        for (i = 0; i < n3; i++) {
            st.dir[i] = -st.grad[i] + beta * st.dir[i];
            slope += st.dir[i] * st.grad[i];
        }
        if (!(slope < 0)) {
            slope = 0.0;
            for (i = 0; i < n3; i++) {
                st.dir[i] = -st.grad[i];
                slope += st.dir[i] * st.grad[i];
            }
        }
        for (i = 0; i < n3; i++) st.prevGrad[i] = st.grad[i];
        st.havePrevGrad = true;
        return slope;
    }

    // ------------------------------------------------------------------
    // conjugate gradient with a Wolfe line search
    //
    // The direction alone is not enough: a first attempt that changed only
    // that measured worse than plain descent, because the old step control
    // accepted any step that lowered the energy at all -- for a conjugate
    // direction usually a very short one -- and every rejection threw the
    // memory away. So the step comes from a Wolfe line search. Sufficient
    // decrease decides whether a step is good enough to take; the curvature
    // along the direction decides whether it was too short, which is what
    // lets the length grow into the valley instead of being tuned against
    // the size of the gradient.
    //
    // Convergence is judged on progress, not on the gradient alone. Measured
    // on three molecules, the gradient test never once fired: cholesterol
    // reached 115.193 kcal/mol by frame 200 and spent the next 1300 frames
    // gaining 0.006 while having 54,000 steps rejected -- about forty per
    // frame, for ever. A strained ring did worse: it stopped improving at
    // frame 50 and then rejected 107,000 steps without moving at all, because
    // a failed step resets the length and the next frame tries the same thing
    // again. Meanwhile a real metal complex was still descending at frame 700
    // and had every right to go on.
    //
    // Energy change per frame tells those three apart, and the gradient
    // cannot: it stops the two that are finished and lets the third run.
    // ------------------------------------------------------------------
    function relax(st, k) {
        var n3 = 3 * st.n, i;
        if (!st.haveEnergy) {
            st.energyValue = evalMasked(st, st.pos, st.grad, T_ALL);
            st.haveEnergy = true;
        }
        st.converged = false;
        st.stalled = false;
        var startEnergy = st.energyValue;
        var taken = 0;
        for (var s = 0; s < k; s++) {
            if (!isFinite(st.energyValue)) { rollback(st, 'non-finite energy'); break; }
            var gmax = 0.0;
            for (i = 0; i < n3; i++) {
                var g = st.grad[i];
                if (g < 0) g = -g;
                if (g > gmax) gmax = g;
            }
            if (!isFinite(gmax)) { rollback(st, 'non-finite gradient'); break; }
            if (gmax < GRAD_TOL) {
                // A fully collapsed geometry produces a zero gradient from
                // every term and used to be reported as converged. It is the
                // opposite of converged.
                if (st.collapsed) { rollback(st, 'coincident atoms'); break; }
                st.converged = true;
                break;
            }
            var slope0 = searchDirection(st, n3);
            var dmax = 0.0;
            for (i = 0; i < n3; i++) {
                var dv = st.dir[i] < 0 ? -st.dir[i] : st.dir[i];
                if (dv > dmax) dmax = dv;
            }
            if (!isFinite(dmax) || dmax < 1e-16) { st.converged = true; break; }
            var lam = st.lambda;
            if (lam * dmax > MAX_DISPLACEMENT) lam = MAX_DISPLACEMENT / dmax;
            var accepted = false, e1 = 0.0, slope1 = 0.0;
            for (var attempt = 0; attempt < LINE_TRIES; attempt++) {
                for (i = 0; i < n3; i++) st.trial[i] = st.pos[i] + lam * st.dir[i];
                e1 = evalMasked(st, st.trial, st.tgrad, T_ALL);
                if (isFinite(e1) && e1 <= st.energyValue + WOLFE_C1 * lam * slope0) {
                    slope1 = 0.0;
                    for (i = 0; i < n3; i++) slope1 += st.dir[i] * st.tgrad[i];
                    accepted = true;
                    break;
                }
                st.rejects++;
                lam *= SHRINK;
                if (lam * dmax < MIN_LAMBDA) break;
            }
            if (accepted) {
                var swap = st.pos; st.pos = st.trial; st.trial = swap;
                swap = st.grad; st.grad = st.tgrad; st.tgrad = swap;
                st.energyValue = e1;
                // The curvature half of Wolfe: a step that left the slope
                // still steeply downhill was too short and the next may be
                // longer; one that overshot into a rising slope was too long.
                if (slope1 < WOLFE_C2 * slope0) {
                    st.lambda = Math.min(lam * GROW, MAX_LAMBDA);
                } else if (slope1 > -WOLFE_C2 * slope0) {
                    st.lambda = Math.max(lam * SHRINK, MIN_LAMBDA);
                } else {
                    st.lambda = lam;
                }
                taken++;
            } else if (st.havePrevGrad) {
                // Nothing along the conjugate direction helps. Forget it and
                // try the plain gradient before giving up on the geometry.
                st.havePrevGrad = false;
                st.lambda = INITIAL_LAMBDA;
                continue;
            } else {
                st.lambda = INITIAL_LAMBDA;
                st.stalled = true;
                break;
            }
        }
        // Did this frame achieve anything? A stall counts as no, which is
        // what turns endless retrying into an answer.
        var moved = Math.abs(st.energyValue - startEnergy);
        if (st.converged) {
            st.quietFrames = QUIET_FRAMES;
        } else if (isFinite(moved) && moved < ENERGY_TOL) {
            st.quietFrames = (st.quietFrames || 0) + 1;
            if (st.quietFrames >= QUIET_FRAMES) {
                st.converged = true;
                st.stalled = false;
            }
        } else {
            st.quietFrames = 0;
        }
        st.totalSteps += taken;
        return taken;
    }

    function rollback(st, why) {
        st.pos.set(st.lastGood);
        st.lambda = INITIAL_LAMBDA;
        st.haveEnergy = false;
        st.rollbacks++;
        st.lastError = why;
    }

    // ------------------------------------------------------------------
    // position marshalling
    // ------------------------------------------------------------------
    function readPositions(st, src, out) {
        var n = st.n, i;
        if (!src) return false;
        if (src.length === 3 * n && typeof src[0] === 'number') {
            for (i = 0; i < 3 * n; i++) out[i] = src[i];
        } else if (src.length === n) {
            for (i = 0; i < n; i++) {
                var p = src[i];
                if (!p) return false;
                if (p.length >= 3) {
                    out[3 * i] = p[0]; out[3 * i + 1] = p[1]; out[3 * i + 2] = p[2];
                } else {
                    out[3 * i] = p.x; out[3 * i + 1] = p.y; out[3 * i + 2] = p.z;
                }
            }
        } else {
            return false;
        }
        return allFinite(out);
    }

    // ------------------------------------------------------------------
    // batch adaptation (Avogadro 2 adaptChunkIterations)
    // ------------------------------------------------------------------
    function adaptChunk(st, elapsedMs, steps) {
        if (!(elapsedMs > 0) || !(steps > 0)) return;
        var perStep = elapsedMs / steps;
        st.msPerStep = perStep;
        if (!(perStep > 0)) return;
        var ideal = TARGET_MS / perStep;
        var next = SMOOTHING * st.chunk + (1.0 - SMOOTHING) * ideal;
        if (!isFinite(next)) return;
        if (next < MIN_CHUNK) next = MIN_CHUNK;
        var ceiling = st.maxChunk || MAX_CHUNK;
        if (next > ceiling) next = ceiling;
        st.chunk = next;
    }

    // ------------------------------------------------------------------
    // public API
    // ------------------------------------------------------------------
    function getState(scopeKey) {
        return window._delfinFFByScope[scopeKey] || null;
    }

    function load(scopeKey, terms) {
        var built;
        try {
            built = buildState(scopeKey, terms);
        } catch (err) {
            return {ok: false, error: 'load failed: ' + err};
        }
        if (built.error) return {ok: false, error: built.error};
        var st = built.state;
        window._delfinFFByScope[scopeKey] = st;
        return {
            ok: true,
            n_atoms: st.n,
            source: st.source,
            bonds: st.nBonds,
            angles: st.nAngles,
            torsions: st.nTors,
            untermed: st.untermed,
            warnings: st.warnings.slice(0)
        };
    }

    function tether(scopeKey, pulls) {
        // Where the hand wants each atom it is holding, and how hard.
        //
        // The whole set at once, because that is what a mouse move is: one
        // message a frame, replacing the last.  An empty list lets go.
        var st = getState(scopeKey);
        if (!st) return false;
        var src = pulls || [];
        var a = new Int32Array(src.length), k = new Float64Array(src.length);
        var x = new Float64Array(src.length), y = new Float64Array(src.length),
            z = new Float64Array(src.length);
        var m = 0;
        for (var t = 0; t < src.length; t++) {
            var one = src[t] || {};
            var at = one.atom | 0;
            if (at < 0 || at >= st.n) continue;
            if (!isFinite(one.x) || !isFinite(one.y) || !isFinite(one.z)) continue;
            if (!isFinite(one.k) || one.k < 0) continue;
            a[m] = at; x[m] = one.x; y[m] = one.y; z[m] = one.z; k[m] = one.k;
            m++;
        }
        st.tetN = m; st.tetA = a; st.tetX = x; st.tetY = y; st.tetZ = z;
        st.tetK = k;
        return true;
    }

    function grab(scopeKey, atomIndices) {
        var st = getState(scopeKey);
        if (!st) return false;
        var n3 = 3 * st.n, i;
        for (i = 0; i < n3; i++) st.mask[i] = 1.0;
        st.frozen = [];
        var list = atomIndices || [];
        if (typeof list.length !== 'number') list = [list];
        for (i = 0; i < list.length; i++) {
            var a = list[i] | 0;
            if (a < 0 || a >= st.n) continue;
            st.frozen.push(a);
            st.mask[3 * a] = 0.0;
            st.mask[3 * a + 1] = 0.0;
            st.mask[3 * a + 2] = 0.0;
        }
        // A fresh grab is a fresh interaction: restart the batch at Avogadro's
        // initial chunk and force a neighbour rebuild on the next frame.
        st.chunk = INITIAL_CHUNK;
        st.maxChunk = st.maxChunk || MAX_CHUNK;
        st.lambda = INITIAL_LAMBDA;
        st.haveEnergy = false;
        st.sinceRebuild = PAIR_REBUILD_FRAMES;
        st.stalled = false;
        st.quietFrames = 0;
        st.converged = false;
        st.lastError = null;
        return true;
    }

    function release(scopeKey) {
        var st = getState(scopeKey);
        if (!st) return false;
        var n3 = 3 * st.n;
        for (var i = 0; i < n3; i++) st.mask[i] = 1.0;
        st.frozen = [];
        st.haveEnergy = false;
        st.quietFrames = 0;
        st.havePrevGrad = false;
        st.converged = false;
        st.lambda = INITIAL_LAMBDA;
        return true;
    }

    /**
     * Run one adaptive batch and return the relaxed positions.
     *
     * positions   flat [x,y,z,...], or N triples, or N {x,y,z} objects.  Pass
     *             the geometry with the dragged atom already moved to the
     *             cursor; omit it to continue from the engine's own state.
     * frameMs     wall time of the PREVIOUS full frame, force field plus
     *             renderer.  The renderer's cost has to sit inside the 33 ms
     *             budget or the batch grows until the viewer stutters.
     *
     * Returns the Float64Array the engine owns -- it is reused every frame,
     * so copy it if you need to keep a snapshot.
     */
    function step(scopeKey, positions, frameMs) {
        var st = getState(scopeKey);
        if (!st) return null;
        try {
            if (positions) {
                // Stage into the probe buffer first: a frame that arrives with
                // a NaN in it must be dropped whole, never half-applied.
                if (!readPositions(st, positions, st.probe)) {
                    st.lastError = 'caller supplied non-finite positions';
                    st.rollbacks++;
                    return st.havePos ? st.pos : null;
                }
                st.pos.set(st.probe);
                // Someone moved the atoms: whatever was converged is not any
                // more, and the count of quiet frames means nothing now.
                st.quietFrames = 0;
                st.converged = false;
                st.haveEnergy = false;
                st.havePos = true;
            } else if (!st.havePos) {
                return null;
            }
            st.lastGood.set(st.pos);

            if (isNum(frameMs) && frameMs > 0 && st.lastSteps > 0) {
                // Adapt from the measured full frame, renderer included.
                adaptChunk(st, frameMs, st.lastSteps);
                st.lastFrameMs = frameMs;
            }
            if (pairsAreStale(st)) rebuildPairs(st);
            st.sinceRebuild++;

            var k = Math.round(st.chunk);
            var ceiling = st.maxChunk || MAX_CHUNK;
            if (k > ceiling) k = ceiling;
            if (k < MIN_CHUNK) k = MIN_CHUNK;
            if (k > MAX_CHUNK) k = MAX_CHUNK;

            var t0 = now();
            var taken = relax(st, k);
            var ffMs = now() - t0;

            st.lastFfMs = ffMs;
            st.lastSteps = taken > 0 ? taken : k;
            st.frames++;
            if (!(isNum(frameMs) && frameMs > 0)) {
                // No renderer timing offered: fall back to our own cost, which
                // under-counts the frame but at least keeps the batch bounded.
                adaptChunk(st, ffMs, st.lastSteps);
            }
            if (!allFinite(st.pos)) rollback(st, 'non-finite positions after step');
            return st.pos;
        } catch (err) {
            st.lastError = 'step failed: ' + err;
            if (st.havePos) st.pos.set(st.lastGood);
            st.rollbacks++;
            return st.havePos ? st.pos : null;
        }
    }

    function energy(scopeKey, positions) {
        var st = getState(scopeKey);
        if (!st) return NaN;
        var target = st.pos;
        if (positions) {
            if (!readPositions(st, positions, st.probe)) return NaN;
            target = st.probe;
        } else if (!st.havePos) {
            return NaN;
        }
        if (!st.pairI) {
            var keep = st.pos;
            st.pos = target;
            rebuildPairs(st);
            st.pos = keep;
        }
        return evaluate(st, target, null, T_ALL);
    }

    function stats(scopeKey) {
        var st = getState(scopeKey);
        if (!st) return null;
        return {
            n_atoms: st.n,
            source: st.source,
            frames: st.frames,
            steps: st.totalSteps,
            last_steps: st.lastSteps,
            chunk: st.chunk,
            max_chunk: st.maxChunk || MAX_CHUNK,
            ms_per_frame: st.lastFrameMs || st.lastFfMs,
            ms_force_field: st.lastFfMs,
            ms_per_step: st.msPerStep,
            energy: st.energyValue,
            frozen: st.frozen.slice(0),
            pairs: st.nPairs,
            pair_rebuilds: st.pairRebuilds,
            rejected_steps: st.rejects,
            rollbacks: st.rollbacks,
            stalled: st.stalled,
            converged: st.converged,
            untermed_atoms: st.untermed,
            last_error: st.lastError,
            warnings: st.warnings.slice(0)
        };
    }

    function positionsOf(scopeKey) {
        var st = getState(scopeKey);
        return st && st.havePos ? st.pos : null;
    }

    function configure(scopeKey, opts) {
        var st = getState(scopeKey);
        if (!st || !opts) return false;
        if (isNum(opts.cutoff) && opts.cutoff > 0) { st.cutoff = opts.cutoff; st.pairI = null; }
        if (isNum(opts.skin) && opts.skin >= 0) { st.skin = opts.skin; st.pairI = null; }
        if (isNum(opts.chunk) && opts.chunk > 0) st.chunk = opts.chunk;
        if (isNum(opts.maxChunk) && opts.maxChunk > 0) {
            st.maxChunk = Math.max(MIN_CHUNK, Math.min(MAX_CHUNK, opts.maxChunk));
            if (st.chunk > st.maxChunk) st.chunk = st.maxChunk;
        }
        return true;
    }

    function dispose(scopeKey) {
        if (window._delfinFFByScope[scopeKey]) {
            delete window._delfinFFByScope[scopeKey];
            return true;
        }
        return false;
    }

    // Diagnostics: energy/gradient for one term type at a time.  The drag loop
    // never calls these; the gradient tests and the debug panel do.
    function termBit(term) {
        if (!term || term === 'all') return T_ALL;
        if (term === 'bond') return T_BOND;
        if (term === 'angle') return T_ANGLE;
        if (term === 'torsion') return T_TORSION;
        if (term === 'vdw') return T_VDW;
        if (term === 'restraints' || term === 'restraint') return T_RESTRAINT;
        if (term === 'tether' || term === 'tethers') return T_TETHER;
        if (term === 'inversion' || term === 'inversions') return T_INVERSION;
        return T_ALL;
    }
    function debugEnergy(scopeKey, positions, term) {
        var st = getState(scopeKey);
        if (!st) return NaN;
        var target = st.pos;
        if (positions) {
            if (!readPositions(st, positions, st.probe)) return NaN;
            target = st.probe;
        }
        if (!st.pairI) rebuildPairs(st);
        return evaluate(st, target, null, termBit(term));
    }
    function debugGradient(scopeKey, positions, term) {
        var st = getState(scopeKey);
        if (!st) return null;
        var target = st.pos;
        if (positions) {
            if (!readPositions(st, positions, st.probe)) return null;
            target = st.probe;
        }
        if (!st.pairI) rebuildPairs(st);
        evaluate(st, target, st.pgrad, termBit(term));
        return st.pgrad;
    }

    window.__delfinFF = {
        load: load,
        grab: grab,
        tether: tether,
        step: step,
        release: release,
        energy: energy,
        stats: stats,
        positions: positionsOf,
        configure: configure,
        dispose: dispose,
        debugEnergy: debugEnergy,
        debugGradient: debugGradient,
        constants: {
            initialChunk: INITIAL_CHUNK,
            targetMs: TARGET_MS,
            smoothing: SMOOTHING,
            minChunk: MIN_CHUNK,
            maxChunk: MAX_CHUNK,
            cutoff: VDW_CUTOFF,
            skin: VDW_SKIN,
            pairRebuildFrames: PAIR_REBUILD_FRAMES,
            maxDisplacement: MAX_DISPLACEMENT
        }
    };
})();
"""


def molecule_ff_bootstrap_js():
    """Return the one-time JS that installs the window.__delfinFF engine.

    The same source is left behind as a string as well.  A relaxation batch
    aims at a whole frame (``TARGET_MS``), and a whole frame spent here is a
    whole frame the page does not have: measured in a browser on a 100-atom
    peptide, one batch is 31 ms of the 17 ms a display frame lasts, so every
    click, every widget update and every message from the kernel waits behind
    the physics.  The viewer therefore runs the engine in a Worker, and a
    Worker needs its program as text.  It is the same text, not a second
    engine -- there is one place where the force field is written.
    """
    return (MOLECULE_FF_BOOTSTRAP_JS
            + '\nwindow.__delfinFFSource = '
            + json.dumps(MOLECULE_FF_BOOTSTRAP_JS) + ';\n')
