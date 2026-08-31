"""Deforming a structure by pulling on it, and letting it choose how to give.

A relaxed scan drives a coordinate somebody chose.  Its weakness is not that
it is slow or coarse, it is that the choice is a guess: measured on the
sixteen-atom Diels-Alder in :func:`delfin.dashboard.gfn_optimize.paths_disagree`,
the walk's summit and the bond that actually forms are 1.2 A apart in a
coordinate nobody was driving.  You cannot drive the coordinate that matters
until you know which one it is, and knowing that is most of the question.

So: do not drive one.  Hang a force on an atom, ramp it up, and relax at every
step.  What gives way is the structure's answer rather than the user's guess --
the path of least resistance under that load.  This is the general case of the
push the scan already has, which is the same thing along the axis of one pair.

The construction is the published one for mechanochemistry, *external force
explicitly included*: a constant force **F** on atom *i* deforms the surface to

    E'(x) = E(x) - sum_i F_i . x_i

and the geometry that minimises E' is the one where the internal gradient
carries the load.  What is *reported* is E at that geometry, never E' -- the
bias is how the structure was got to, not part of what it costs, the same
distinction the restraint-free energies in ``gfn_optimize`` are careful about.

Two things have to be true or nothing converges, and both are about the fact
that a molecule in free space has six directions in which it costs nothing to
move:

* A load with a net force accelerates the molecule.  E' is then unbounded
  below -- travelling in the direction of the net force lowers it for ever --
  and a minimiser walks the structure off to infinity without deforming it at
  all.
* A load with a net torque spins it, for the same reason.

Both are removed by projecting the load onto the space orthogonal to rigid
motion at the geometry it is applied to.  What is left deforms and nothing
else, and E' is then invariant under rigid motion, bounded, and minimisable.
The part that was removed is reported rather than silently dropped: it is
usually the user meaning "pull these apart" with one arrow, and a pair is what
they wanted.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import climb as _climb
from . import gfn_optimize as _gfn


#: One Hartree per Bohr, in the kcal/mol/A the editor's forces are drawn in.
#:
#: The editor speaks in kcal/mol/A everywhere a force is shown -- the push
#: ramps from :data:`gfn_optimize.PUSH_FORCE_FROM` to ``PUSH_FORCE_TO`` in
#: them, and :data:`gfn_optimize.A_BOND_HOLDS` is what a bond withstands in
#: them -- so a load drawn by hand is in them too, and this is the only place
#: it becomes atomic units.
FORCE_IN_KCAL_PER_ANGSTROM = _climb.HARTREE_IN_KCAL / _climb.BOHR

#: How much of a load may be net force or torque before it is worth saying so.
#:
#: A share of the load's own size rather than an absolute: a gentle pull that
#: is 90 per cent net translation is doing nothing, and a hard one that is 2
#: per cent is doing what was asked.
RIGID_SHARE_WORTH_SAYING = 0.05

#: How far any one atom may travel while one level of the load settles, in A.
#:
#: A trust region, and it is about the method rather than about the chemistry:
#: L-BFGS's first step on a soft direction under a hard load can be enormous,
#: and a geometry a long way from anything chemical is where a semiempirical
#: method stops converging at all.  One level is a small step from the one
#: before it by construction, so this never binds on a walk that is going well.
MOST_A_LEVEL_MOVES = 0.8

#: And how far it may travel while the load comes off, in A.
#:
#: Wider, because this is the one step that is allowed to be large: a bond
#: that has just given lets two pieces fall away from each other, and the
#: trust region that keeps a *level* honest would stop that halfway and call
#: the result a minimum.
SETTLE_REACH = 3.0


def rigid_directions(coords: Any) -> np.ndarray:
    """The six ways a molecule can move without changing shape, orthonormal.

    Three translations and three rotations about the centre, as vectors in the
    3N space the forces live in.  A linear molecule has only two rotations and
    a single atom none at all; the ones that vanish come out as null vectors
    and are dropped, so the answer is however many there really are.

    Not mass weighted.  These are used to take the rigid part out of a *force*,
    and a force pattern that translates the molecule is one that is the same on
    every atom -- which is the unweighted vector.  Mass weighting belongs to
    normal modes, where the question is what the motion costs rather than what
    the load is.
    """
    xyz = np.asarray(coords, dtype=float).reshape(-1, 3)
    n = len(xyz)
    if n == 0:
        return np.zeros((0, 0))
    middle = xyz - xyz.mean(axis=0)
    basis = []
    for axis in range(3):
        one = np.zeros((n, 3))
        one[:, axis] = 1.0
        basis.append(one.ravel())
    for axis in range(3):
        turn = np.zeros(3)
        turn[axis] = 1.0
        basis.append(np.cross(turn, middle).ravel())
    kept: List[np.ndarray] = []
    for one in basis:
        for already in kept:
            one = one - float(one @ already) * already
        size = float(np.linalg.norm(one))
        # A rotation of a single atom about the axis it sits on, or of a linear
        # molecule about its own axis, is not a direction: it is zero, and
        # normalising it would be a division by nothing.
        if size > 1e-8:
            kept.append(one / size)
    return np.array(kept) if kept else np.zeros((0, 3 * n))


def deforming_part(forces: Any, coords: Any) -> Tuple[np.ndarray, float]:
    """The part of a load that deforms, and what share of it did not.

    Returns the projected forces and the share of the original that was rigid
    -- net push and net spin together -- as a number between 0 and 1.

    A load that is entirely rigid comes back as zeros and a share of 1: one
    arrow on one atom is exactly that, and it is worth knowing that it moves
    the molecule and bends nothing.
    """
    load = np.asarray(forces, dtype=float).ravel()
    size = float(np.linalg.norm(load))
    if size < 1e-12:
        return load.reshape(-1, 3), 0.0
    rigid = rigid_directions(coords)
    left = load.copy()
    for one in rigid:
        left = left - float(left @ one) * one
    return left.reshape(-1, 3), max(0.0, 1.0 - float(np.linalg.norm(left)) / size)


def _turned_onto(moved: Any, reference: Any) -> np.ndarray:
    """*moved* rotated and shifted to lie over *reference*, shape unchanged.

    Kabsch.  The load is made to deform and not to move, but a minimiser takes
    finite steps on a surface whose rigid invariance is exact only at the
    geometry the load was projected at -- so a long relaxation can still drift
    a little.  Putting the answer back where it started costs nothing and means
    the frames of a walk sit still on screen instead of swimming.
    """
    a = np.asarray(moved, dtype=float).reshape(-1, 3)
    b = np.asarray(reference, dtype=float).reshape(-1, 3)
    if len(a) != len(b) or len(a) < 2:
        return a
    here, there = a - a.mean(axis=0), b - b.mean(axis=0)
    u, _s, vt = np.linalg.svd(here.T @ there)
    turn = u @ vt
    if np.linalg.det(turn) < 0:
        u[:, -1] *= -1.0
        turn = u @ vt
    return here @ turn + b.mean(axis=0)


class _Loaded:
    """One structure, one engine, and a load that can be scaled."""

    def __init__(self, symbols, numbers, angstrom, method,
                 charge, uhf, solvent, cores):
        self.symbols = list(symbols)
        self.numbers = numbers
        self.bohr = np.asarray(angstrom, dtype=float) / _climb.BOHR
        # GFN-FF goes out through the command line even where the library is
        # there, for the reason Climb gives: in process it drops its topology
        # file into whatever directory the dashboard was launched from.
        fast = _climb.have_fast_gradients() and method != 'gfnff'
        maker = _climb._InProcess if fast else _climb._CommandLine
        self.engine = maker(self.numbers, self.bohr, method,
                            charge, uhf, solvent, cores)
        self.fast = bool(fast)

    def close(self) -> None:
        closer = getattr(self.engine, 'close', None)
        if callable(closer):
            try:
                closer()
            except Exception:
                pass


def relax_under_load(loaded: _Loaded, bohr: Any, forces: Any, *,
                     steps: int = 200, reach: Optional[float] = None,
                     should_stop: Optional[Callable[[], bool]] = None,
                     ) -> Dict[str, Any]:
    """Let the structure settle under a load it cannot move away from.

    *forces* are already the deforming part, in Hartree per Bohr, one row an
    atom.  Returns ``{'ok', 'bohr', 'energy', 'gradients', 'status'}`` where
    *energy* is the plain energy of the geometry it settled at -- not the
    biased one it was found by.

    L-BFGS with the analytic gradient, which is what makes this affordable: the
    engine answers energy and gradient together in one call, and a load level
    is a small displacement from the one before it, so a few dozen calls settle
    it.  Measured in :class:`delfin.dashboard.climb.Climb`, a gradient is about
    6 ms for sixteen atoms in process.
    """
    from scipy.optimize import minimize

    load = np.asarray(forces, dtype=float).ravel()
    start = np.asarray(bohr, dtype=float).ravel()
    spent = {'calls': 0, 'stopped': False, 'refused': 0}
    #: The last point the method could answer about, and its biased value.
    good: Dict[str, Any] = {'x': start.copy(), 'value': None}

    def biased(flat):
        if should_stop is not None and should_stop():
            spent['stopped'] = True
            # Nothing raises: a stop is not an error, and the geometry reached
            # so far is the answer.  Handing back a zero gradient is how a
            # minimiser is told it has arrived.
            return 0.0, np.zeros_like(flat)
        try:
            energy, gradient = loaded.engine(flat.reshape(-1, 3))
        except Exception:
            # A self-consistent field that did not converge is not a value,
            # and it is not always a mistake either: pull hard enough and the
            # bond goes, and two radicals are exactly where a closed-shell
            # method stops having an answer.  Refused as a wall rather than
            # raised -- uphill, with a gradient pointing back at the last
            # point that worked -- so the line search backs off instead of the
            # whole level ending on an exception.
            spent['refused'] += 1
            away = flat - good['x']
            floor = good['value'] if good['value'] is not None else 0.0
            return floor + 1.0 + float(away @ away), 2.0 * away
        spent['calls'] += 1
        flat_gradient = np.asarray(gradient, dtype=float).ravel()
        # E' = E - F.x, so dE'/dx = dE/dx - F.
        value = float(energy) - float(load @ flat)
        good['x'] = flat.copy()
        good['value'] = value
        return value, (flat_gradient - load)

    # A trust region, as bounds.  L-BFGS's first step on a soft direction under
    # a hard load can be enormous, and a geometry a long way from anything
    # chemical is where the method stops converging -- so the refusal above
    # would be doing all the work.  One level is a small step from the one
    # before it by construction, so a bound this wide never binds on a walk
    # that is going well and catches only the ones that are not.
    far = (MOST_A_LEVEL_MOVES if reach is None else float(reach)) / _climb.BOHR
    edges = [(float(v) - far, float(v) + far) for v in start]
    answer = minimize(biased, start, jac=True, method='L-BFGS-B',
                      bounds=edges,
                      options={'maxiter': int(steps), 'gtol': 1e-5,
                               'ftol': 1e-12})
    settled = np.asarray(answer.x, dtype=float).reshape(-1, 3)
    try:
        # The plain energy of what it settled at, which is what is reported.
        # Asked again rather than taken from the last biased call: the
        # minimiser's last evaluation is not always the point it returns.
        energy, _gradient = loaded.engine(settled)
        spent['calls'] += 1
    except Exception:
        return {'ok': False, 'bohr': settled, 'energy': None,
                'gradients': spent['calls'], 'stopped': spent['stopped'],
                'refused': spent['refused'],
                'status': 'the method had no answer for this geometry'}
    return {'ok': True, 'bohr': settled, 'energy': float(energy),
            'gradients': spent['calls'], 'stopped': spent['stopped'],
            'refused': spent['refused'],
            'status': str(getattr(answer, 'message', '') or '')}


def walk_under_load(xyz_text: str, loads: Sequence[Dict[str, Any]],
                    method: str = 'gfn2', *,
                    charge: int = 0, uhf: int = 0,
                    solvent: Optional[str] = None, cores: int = 4,
                    steps: int = 20,
                    whole: bool = False,
                    force_from: Optional[float] = None,
                    force_to: Optional[float] = None,
                    on_point: Optional[Callable[[Dict[str, Any]], None]] = None,
                    should_stop: Optional[Callable[[], bool]] = None,
                    ) -> Dict[str, Any]:
    """Ramp a load from gentle to hard, and say what the structure did.

    *loads* is what the user drew: ``[{'atom': i, 'vector': (x, y, z)}]`` with
    the vector in kcal/mol/A.  The drawing sets the directions and how hard the
    arrows pull *relative to each other*; the ramp sets how hard the whole set
    pulls, from *force_from* to *force_to* geometrically -- the same shape of
    ramp, and the same two ends, as the push the scan already has.

    Returns ``{'ok', 'points', 'rigid', 'status'}``.  Each point is
    ``{'force', 'energy', 'xyz', 'gradients'}``: the strength of the largest
    arrow at that level in kcal/mol/A, the plain energy in kcal/mol against the
    first point, and the geometry.

    The load is projected once at each level, at the geometry it is applied to,
    so what is applied deforms and nothing else -- see :func:`deforming_part`.
    Every level starts from the geometry the level before it settled at, which
    is what makes each one cheap and what makes the sequence a path rather than
    a set of unrelated answers.
    """
    found = _climb._elements(xyz_text)
    if found is None:
        return {'ok': False, 'points': [], 'rigid': 0.0,
                'status': 'There is no structure to pull on.'}
    symbols = found['symbols']
    angstrom = np.asarray(found['angstrom'], dtype=float)
    n = len(symbols)
    if n == 0:
        return {'ok': False, 'points': [], 'rigid': 0.0,
                'status': 'There is no structure to pull on.'}

    drawn = np.zeros((n, 3))
    for one in loads or ():
        index = int(one.get('atom', -1))
        if 0 <= index < n:
            drawn[index] += np.asarray(one.get('vector') or (0.0, 0.0, 0.0),
                                       dtype=float)
    biggest = float(np.max(np.linalg.norm(drawn, axis=1))) if n else 0.0
    if biggest < 1e-9:
        return {'ok': False, 'points': [], 'rigid': 0.0,
                'status': 'No arrow has been drawn, so there is nothing to '
                          'pull with.'}
    shape = drawn / biggest          # directions and relative strengths

    low = float(force_from if force_from is not None else _gfn.PUSH_FORCE_FROM)
    high = float(force_to if force_to is not None else _gfn.PUSH_FORCE_TO)
    many = max(2, int(steps))
    growth = (high / low) ** (1.0 / (many - 1))

    loaded = _Loaded(symbols, found['numbers'], angstrom, method,
                     charge, uhf, solvent, cores)
    points: List[Dict[str, Any]] = []
    rigid_share = 0.0
    was_bonded: Any = None
    gave: Optional[Dict[str, Any]] = None
    settled: Optional[Dict[str, Any]] = None
    try:
        here = loaded.bohr.copy()
        # The zero is the structure as it was handed over, unloaded.
        #
        # It used to be the first *level* -- already under the gentlest load --
        # so every energy in the walk was quoted against a structure that was
        # itself slightly bent, and the whole ramp was measured from somewhere
        # nobody had been.  Measured: an ethane pulled and released came back
        # at -0.05 kcal/mol against its own starting point, and the twentieth
        # was the reference relaxing, not the molecule.
        #
        # A single point and no relaxation: the structure the user handed over
        # is the structure they meant, and moving it to make a nicer zero is
        # answering about a molecule they did not ask about.
        try:
            first = float(loaded.engine(here)[0])
        except Exception:
            first = None
        for level in range(many):
            if should_stop is not None and should_stop():
                break
            strength = low * (growth ** level)
            # Projected at the geometry it is applied to: the rigid directions
            # turn with the molecule, so a load balanced at the start is not
            # balanced after it has bent.
            wanted = shape * strength / FORCE_IN_KCAL_PER_ANGSTROM
            forces, share = deforming_part(wanted, here)
            rigid_share = max(rigid_share, share)
            got = relax_under_load(loaded, here, forces,
                                   should_stop=should_stop)
            if not got.get('ok'):
                break
            # Back where it started, so a walk of these frames sits still on
            # screen rather than swimming: see :func:`_turned_onto`.
            here = _turned_onto(got['bohr'], loaded.bohr)
            point = {
                'force': strength,
                'energy': ((float(got['energy']) - first)
                           * _climb.HARTREE_IN_KCAL) if first is not None else 0.0,
                'xyz': _climb.xyz_document(
                    symbols, here * _climb.BOHR,
                    f'pulled at {strength:.0f} kcal/mol/A'),
                'gradients': int(got.get('gradients') or 0),
            }
            points.append(point)
            if on_point is not None:
                try:
                    on_point(point)
                except Exception:
                    pass
            if got.get('stopped'):
                break
            # What gave, and at what load.
            #
            # This is the answer, not a diagnostic.  A ramp is asked "which
            # bond goes first and how hard do I have to pull" -- the whole
            # reason for not driving a coordinate is that the first half of
            # that is not known in advance.  Measured on an ethane pulled
            # along its C-C under GFN2: it holds at 77 kcal/mol/A and is
            # broken by 113, which is where A_BOND_HOLDS says a bond gives.
            graph = _gfn.bond_graph(point['xyz'])
            if was_bonded is None:
                was_bonded = graph
            elif graph != was_bonded and gave is None:
                gave = {
                    'said': _gfn.graph_changed(was_bonded, graph, symbols),
                    'held': points[-2]['force'] if len(points) > 1 else None,
                    'broke': point['force'],
                }
            elif gave is not None and not whole:
                # One more level, and then stop.  A ramp goes from one settled
                # structure to the next: the bond gave, the pieces relaxed,
                # and that is a minimum.  Pressing again carries on from the
                # load it stopped at, which is the next one.
                #
                # Past the break the load is otherwise only drawing two pieces
                # apart -- every further level by the trust region -- and says
                # nothing new.  "Whole profile" is asking for exactly that:
                # the whole ramp, whatever gives on the way.
                break
        # And then the load comes off.
        #
        # Every level is a minimum of the *loaded* surface, which is a real
        # thing and not the thing anybody wants to keep: it is a structure
        # held out of shape by a force that is about to stop existing.  Handed
        # back, it goes to Optimise and to the queue strained, and the energy
        # quoted for it is the energy of the strain as much as of the
        # chemistry.
        #
        # So the arrows are taken off and it settles once more.  That is the
        # far end of "from one minimum to the next", and it is what the walk
        # was for: the load is how the structure was got over the hill, not
        # part of where it landed.
        if points and not (should_stop is not None and should_stop()):
            free = relax_under_load(loaded, here, np.zeros_like(here),
                                    reach=SETTLE_REACH,
                                    should_stop=should_stop)
            if free.get('ok') and free.get('energy') is not None:
                rested = _turned_onto(free['bohr'], loaded.bohr)
                text = _climb.xyz_document(
                    symbols, rested * _climb.BOHR,
                    'pulled, and settled with the load off')
                # And whether anything survived the load coming off.
                #
                # This is the question, not a detail.  Measured on an ethane
                # pulled along its C-C under GFN2: at 164 kcal/mol/A the bond
                # is 4.899 A apart and the structure is 130 kcal/mol up, which
                # reads as a broken bond and a barrier.  Released, it comes
                # back to 1.521 A and -0.1 kcal/mol -- two radicals held apart
                # by a force recombine the moment the force stops, and all but
                # a tenth of a kcal/mol of that 130 was the strain of holding
                # them.  Reported without this, a pull says a bond broke every
                # time it is pulled hard enough, which is a claim about the
                # load and not about the molecule.
                settled = {
                    'xyz': text,
                    'energy': ((float(free['energy']) - first)
                               * _climb.HARTREE_IN_KCAL),
                    'gradients': int(free.get('gradients') or 0),
                    'said': _gfn.graph_changed(
                        _gfn.bond_graph(points[0]['xyz']),
                        _gfn.bond_graph(text), symbols),
                }
    finally:
        loaded.close()

    if not points:
        return {'ok': False, 'points': [], 'rigid': rigid_share, 'gave': None,
                'settled': None,
                'status': 'The load could not be applied.'}
    return {'ok': True, 'points': points, 'rigid': rigid_share, 'gave': gave,
            'settled': settled, 'status': ''}
