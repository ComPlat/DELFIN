"""Let ORCA's optimisers drive g-xTB, which ORCA has no keyword for.

ORCA 6.1.1 knows ``XTB1``, ``XTB2`` and ``XTBFF`` and drives them through its
own bundled ``otool_xtb`` -- an xtb 6.7.1 from 2024.  g-xTB is not among them
and cannot be: the method ships as a statically linked xtb 6.7.1 of its own,
built in 2026, and an ordinary xtb accepts ``--gxtb`` and silently runs GFN2
instead.  So ``! XTB2 OPTTS`` has no g-xTB counterpart, and the most accurate
method the editor offers -- the one somebody picks when a transition state has
to be worth quoting -- had no route to a saddle at all.

It has one, and it is the interface ORCA publishes for exactly this: **ExtOpt**.
Given ``! ExtOpt``, ORCA writes ``<base>_EXT.extinp.tmp`` naming a geometry, a
charge, a multiplicity, a core count and whether it wants a gradient, calls the
program named by ``EXTOPTEXE``, and reads ``<base>_EXT.engrad`` back.  Whatever
answers is the method.  It is the same shape GRRM uses to drive xtb through
``%link=non-supported``, and this module is what answers it with g-xTB.

What it costs was measured here rather than assumed, on the sixteen-atom
Diels-Alder estimate, ORCA 6.1.1 with ``! ExtOpt OptTS``:

* ``Calc_Hess true`` on its own stops in PROPINT -- "ERROR (SHARK): Failed to
  read input file" -- because ORCA's analytic Hessian has no basis set to work
  from.  ``NumHess true`` is ORCA's own remedy and is what is asked for.
* with ``Recalc_Hess 5``, ORCA's usual cadence, it converged in **268 s** over
  sixteen cycles and four numerical Hessians.  The Hessians are nearly all of
  that: each is forty-odd separate g-xTB processes with an ORCA start-up
  apiece.
* with **one** Hessian and Bofill updates after it, it is faster and wrong: the
  two forming bonds walked from 2.524 A out to 3.03 and 3.05, the run hit its
  sixty-cycle bound, and what it left had no imaginary mode at all.  So the
  recalculation stays.

It is still not free, and the reason is not the SCF: a g-xTB gradient is a
whole process, and a process is where the cost is -- see :func:`answer`.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

__all__ = ['read_request', 'engrad_document', 'answer', 'write_hook',
           'BOHR_IN_ANGSTROM']

#: Angstrom in a Bohr, from the same place the rest of the editor takes it.
#: Imported lazily so the hook -- which ORCA starts forty-odd times over one
#: numerical Hessian -- does not pull the dashboard in behind it.
BOHR_IN_ANGSTROM = 0.52917721067

#: Element symbols by atomic number, which is all the periodic table this
#: needs: ORCA's engrad format carries atomic numbers and an xyz carries
#: symbols.  As far as g-xTB itself goes, which is radon.
_SYMBOLS = (
    '', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
    'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
    'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
    'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
    'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
    'Rn',
)
_NUMBER_OF = {name.upper(): z for z, name in enumerate(_SYMBOLS) if name}

#: What ORCA calls the file it wants back, given what it wrote.
_REQUEST_SUFFIX = '.extinp.tmp'

#: Where g-xTB is run, under the directory ORCA is working in.
#:
#: Its own directory rather than ORCA's, because xtb leaves ``xtbrestart``,
#: ``charges``, ``wbo`` and ``xtbtopo.mol`` behind wherever it runs and ORCA's
#: own files are in there.  Kept between calls rather than made fresh, because
#: ``xtbrestart`` is the previous wavefunction and starting the SCF from it is
#: measurably cheaper: on a fifty-atom structure a gradient in a kept
#: directory took 0.12 s against 0.19 in a fresh one.
#:
#: One per request name, and that is not tidiness.  ORCA computes a numerical
#: Hessian **in parallel** -- with ``%pal nprocs 8`` it runs eight displaced
#: geometries at once, each a request of its own named ``in_D00024_EXT`` and
#: the like.  Sharing one directory, they overwrite each other's ``in.xyz``
#: and ``gradient`` mid-flight: measured, a sixteen-atom Hessian of 46
#: displacements came back with numbers 24 and 27 missing, ORCA reported "the
#: Numerical calculation ISN'T COMPLETE", could not write its Hessian, and the
#: whole search ended in nine seconds having done nothing.  Named after the
#: request, two calls collide only if ORCA gave them the same name, and then
#: they are the same call.
_WORK_PREFIX = '_gxtb_'


def _find_gxtb() -> Optional[str]:
    """The g-xTB binary, through the resolver the rest of DELFIN uses.

    Imported here rather than at the top so that this module can be run as a
    hook by ORCA without importing the dashboard: it is started once per
    gradient, and one numerical Hessian is dozens of them -- forty-six for
    sixteen atoms, with translation invariance taken out.

    The environment first, which is the same precedence
    :func:`gfn_optimize.find_gxtb` uses and here is also what keeps the
    dashboard out of the loop: measured, a bare interpreter starts in 37 ms
    and one that imports ``gfn_optimize`` in 119, so being told where the
    binary is saves 80 ms on every gradient -- four seconds over one Hessian,
    and :mod:`delfin.dashboard.saddle` knows the answer already because it
    checked before starting ORCA.
    """
    told = os.environ.get('DELFIN_GXTB_BINARY')
    if told and Path(told).is_file():
        return told
    try:
        from delfin.dashboard.gfn_optimize import find_gxtb
    except Exception:
        return None
    return find_gxtb()


def read_request(path: Any) -> Dict[str, Any]:
    """What ORCA asked for, out of the file it wrote.

    Six lines, each with its meaning in a comment after a hash: the xyz to
    read, the charge, the multiplicity, how many cores it is lending, whether
    it wants a gradient, and a point-charge file.  The comments are ORCA's
    own; only what is in front of them is read.
    """
    rows: List[str] = []
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        said = line.split('#')[0].strip()
        if said:
            rows.append(said)
    if len(rows) < 5:
        raise ValueError(f'{path} is not an ORCA ExtOpt request')
    return {
        'xyz': rows[0],
        'charge': int(float(rows[1])),
        'multiplicity': int(float(rows[2])),
        'cores': max(1, int(float(rows[3]))),
        'gradient': int(float(rows[4])) != 0,
        'point_charges': rows[5] if len(rows) > 5 else '',
    }


def _read_xyz(path: Path) -> Dict[str, Any]:
    """Symbols and Angstrom coordinates, counting the body rather than the header.

    The count in the first line is not trusted, for the reason
    :func:`gfn_optimize.atom_lines` gives: a block whose header disagrees with
    its body silently loses atoms, and every other part of the run succeeds.
    """
    raw = path.read_text(encoding='utf-8').splitlines()
    body: List[List[str]] = []
    for line in raw:
        parts = line.split()
        if len(parts) < 4:
            continue
        if parts[0].upper() not in _NUMBER_OF:
            continue
        try:
            float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        body.append(parts[:4])
    if not body:
        raise ValueError(f'{path} holds no atoms this can read')
    return {
        'symbols': [row[0] for row in body],
        'angstrom': [[float(row[1]), float(row[2]), float(row[3])]
                     for row in body],
    }


def _read_turbomole_gradient(path: Path, count: int):
    """Energy and gradient out of what ``xtb --grad`` leaves behind.

    The file rather than the printed energy, which carries fewer digits.
    """
    rows = path.read_text(encoding='utf-8', errors='replace').splitlines()
    energy = None
    for row in rows:
        if 'energy =' in row:
            try:
                energy = float(row.split('energy =')[1].split()[0])
            except (IndexError, ValueError):
                pass
    body = [row for row in rows if row.strip() and not row.startswith('$')]
    wanted = body[1 + count:1 + 2 * count]
    grad = []
    for row in wanted:
        grad.append([float(word.replace('D', 'E').replace('d', 'E'))
                     for word in row.split()])
    if energy is None or len(grad) != count:
        raise RuntimeError('g-xTB wrote no gradient')
    return energy, grad


def engrad_document(symbols: Sequence[str],
                    angstrom: Sequence[Sequence[float]],
                    energy: float,
                    gradient: Sequence[Sequence[float]]) -> str:
    """The answer in ORCA's own engrad format.

    Three blocks and their headings, exactly as ORCA writes them itself: how
    many atoms, the energy in Hartree, the gradient in Hartree per Bohr one
    component to a line, and the atoms with their coordinates in Bohr.  ORCA
    checks the atom count against its own and stops if it disagrees, which is
    the one mistake this could make silently.
    """
    out = ['#', '# Number of atoms', '#', f' {len(symbols)}',
           '#', '# The current total energy in Eh', '#',
           f'  {float(energy):.12f}',
           '#', '# The current gradient in Eh/bohr', '#']
    for one in gradient:
        for value in one:
            out.append(f'  {float(value):.12f}')
    out += ['#', '# The atomic numbers and current coordinates in Bohr', '#']
    for name, place in zip(symbols, angstrom):
        out.append('%4d %18.12f %18.12f %18.12f' % (
            _NUMBER_OF[str(name).upper()],
            place[0] / BOHR_IN_ANGSTROM,
            place[1] / BOHR_IN_ANGSTROM,
            place[2] / BOHR_IN_ANGSTROM))
    return '\n'.join(out) + '\n'


def answer(request_path: Any, binary: Optional[str] = None) -> int:
    """Answer one of ORCA's requests with g-xTB.  Zero when it worked.

    A gradient is always computed, whether ORCA asked for one or not: in g-xTB
    the gradient is not what the call costs.  Measured on water, a whole
    ``--grad`` process is 0.114 s and the SCF inside it is 0.004 -- the rest is
    the parameters and the basis being set up, which every process pays again.
    That fixed cost is also why this is a saddle *button* and not the
    interactive climb next door: at sixteen atoms one g-xTB gradient is 0.19
    to 0.29 s, against the 6 ms :mod:`delfin.dashboard.climb` measured for a
    GFN2 gradient taken in process -- and a climb answers ten times a second.
    """
    request_path = Path(request_path).resolve()
    asked = read_request(request_path)
    folder = request_path.parent
    if asked['point_charges'] and asked['point_charges'] not in ('-', 'none'):
        sys.stderr.write(
            'g-xTB through ORCA has no external point charges: ORCA named '
            f'{asked["point_charges"]} and nothing here can apply it.\n')
        return 2
    binary = binary or _find_gxtb()
    if binary is None:
        sys.stderr.write(
            'g-xTB was not found. It is a build of its own, xtb-gxtb -- an '
            'ordinary xtb accepts --gxtb and silently runs GFN2 instead.\n')
        return 3

    found = _read_xyz(folder / asked['xyz'])
    count = len(found['symbols'])
    base = request_path.name[:-len(_REQUEST_SUFFIX)] \
        if request_path.name.endswith(_REQUEST_SUFFIX) else request_path.stem
    work = folder / (_WORK_PREFIX + base)
    work.mkdir(exist_ok=True)
    # The previous answer goes before the next question is asked.  The
    # directory is deliberately kept between calls for its xtbrestart, and a
    # run that ends without writing a gradient would otherwise leave the last
    # one there to be read as this one's -- the same geometry twice, and
    # nothing to say it happened.
    try:
        (work / 'gradient').unlink()
    except OSError:
        pass
    (work / 'in.xyz').write_text(
        '%d\nfor ORCA\n' % count + '\n'.join(
            '%-3s %18.12f %18.12f %18.12f' % (name, *place)
            for name, place in zip(found['symbols'], found['angstrom'])
        ) + '\n', encoding='utf-8')

    # Threads: what ORCA said it was lending, unless DELFIN has said otherwise.
    # It is set in the child's environment rather than in this process,
    # because OpenMP reads the variable when the runtime starts and not
    # afterwards -- measured on this box, an unpinned sixteen-atom gradient
    # costs 1.66 s where four threads cost 6 ms, and setting the variable
    # after the fact changed nothing at all.
    cores = os.environ.get('DELFIN_GXTB_CORES') or str(asked['cores'])
    room = dict(os.environ)
    room['OMP_NUM_THREADS'] = room['MKL_NUM_THREADS'] = str(cores)
    room.setdefault('OMP_STACKSIZE', '1G')
    order = [binary, 'in.xyz', '--gxtb', '--grad',
             '--chrg', str(asked['charge']),
             '--uhf', str(max(0, asked['multiplicity'] - 1))]
    try:
        done = subprocess.run(order, cwd=str(work), env=room,
                              capture_output=True, text=True, timeout=3600)
    except (OSError, subprocess.SubprocessError) as exc:
        sys.stderr.write(f'g-xTB did not run: {exc}\n')
        return 4
    if done.returncode != 0:
        # The last thousand characters of its own complaint, because ORCA
        # prints only an exit code and a run that failed for a reason is worth
        # more than a number.
        sys.stderr.write((done.stdout or '')[-1200:] + (done.stderr or ''))
        return 5

    energy, gradient = _read_turbomole_gradient(work / 'gradient', count)
    (folder / (base + '.engrad')).write_text(
        engrad_document(found['symbols'], found['angstrom'], energy, gradient),
        encoding='utf-8')
    return 0


def write_hook(folder: Any, name: str = 'delfin_gxtb_external') -> Path:
    """Write the executable ORCA is pointed at, and return where it is.

    ``EXTOPTEXE`` has to name something the operating system can execute, and
    what has to run is this module inside the interpreter the dashboard is
    already running -- so a two-line shell script that execs it is written
    into the run's own directory.  Naming the interpreter rather than
    ``python`` matters: a dashboard is very often in an environment that is
    not the one on the PATH.
    """
    folder = Path(folder)
    hook = folder / name
    package = str(Path(__file__).resolve().parent.parent.parent)
    hook.write_text(
        '#!/bin/sh\n'
        f'PYTHONPATH="{package}${{PYTHONPATH:+:$PYTHONPATH}}" '
        f'exec "{sys.executable}" -m delfin.dashboard.gxtb_engrad "$@"\n',
        encoding='utf-8')
    hook.chmod(0o755)
    return hook


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        sys.stderr.write('usage: gxtb_engrad <base>_EXT.extinp.tmp\n')
        return 1
    try:
        return answer(argv[0])
    except Exception as problem:            # ORCA sees only the exit code
        sys.stderr.write(f'{type(problem).__name__}: {problem}\n')
        return 6


if __name__ == '__main__':
    raise SystemExit(main())
