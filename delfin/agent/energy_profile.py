"""The group's energy-profile figure (HPES): ``make_plot`` kind ``energy_profile``.

One layout for every reaction energy diagram the agent draws, so a profile
looks the same whoever asked for it, however many steps and pathways it has.
It is the layout the group's own ``HPES*.py`` scripts drew by hand:

* every state is a flat bar at its energy, the value printed under it with
  two decimals; a computed transition state is a state like any other
  (``TS1``), with its bar and value, joined straight to its neighbours;
* consecutive states of a pathway are joined by dashes in the pathway's
  colour -- what lies between two computed states is not to scale;
* a transition state that was NOT computed is a dashed hump between its two
  neighbours, rising above both, with no bar, no value and no column of its
  own -- it says "there is a barrier here" without claiming a height;
* several pathways share one axis, coloured tab10 in order; a state that two
  pathways reach at the same energy is ONE bar split into their colours and
  labelled once;
* serif type, only the left spine, a bold title, bold state names as the x
  ticks and a bold ``ΔG (kcal/mol)`` y label.

Where the scripts depended on the data, the figure does not: the bar
thickness and the label gap are in points, not in energy units (a 7 kcal/mol
profile drew slabs and a 100 kcal/mol one hairlines), and the y range pads
the span rather than the extremes (a maximum of 0 left no room above it).
A value label that would sit on another pathway's bar at the same state goes
above its own bar instead. The figure widens with the number of steps, and
the type shrinks and the state names turn once the columns get narrow, so
two steps and eighty are both legible.

Input -- a list of pathways, each ``(label, states, values, start)``:
``states`` names the x ticks, ``values`` are the energies, ``None`` for an
uncomputed transition state. A pathway that names a state an earlier
pathway already placed is lined up on it -- a branch ``4, 5, 6`` sits under
states 4 to 6 of the main route -- so nobody has to count columns; ``start``
(default 1) places only a pathway that shares no state. The live run that
asked for this passed ``start: true`` and then ``start: 0`` for a branch,
which would have drawn it two columns off; a name cannot be miscounted.
"""

from __future__ import annotations

import re
import textwrap
from collections import Counter
from dataclasses import dataclass

#: Width of one state bar, in units of the spacing between states.
BAR_WIDTH = 0.5
#: Bar thickness and connector width, in points.
BAR_LW = 4.0
JOIN_LW = 1.5
#: Two values that print the same at two decimals are the same state.
_SAME = 0.005
#: Labels closer than this fraction of the y span collide.
_CROWDED = 0.07
#: How far an uncomputed transition state's hump rises above the higher of
#: its two neighbours, as a fraction of the y span.
_HUMP = 0.12
#: Figure size: 10 x 6 in up to ten columns, then wider per column, capped.
_HEIGHT = 6.0
_MIN_WIDTH, _MAX_WIDTH, _PER_COLUMN = 10.0, 32.0, 0.85
#: A column narrower than this (inches) gets smaller type.
_ROOMY = 0.85


@dataclass(frozen=True)
class Pathway:
    label: str
    states: tuple[str, ...]
    values: tuple[float | None, ...]
    start: int = 1
    #: Columns of the computed states, set by :func:`placed`; without them
    #: the states follow ``start`` one column each.
    cols: tuple[int, ...] | None = None

    def points(self):
        """``(position, state, value)`` of the computed states, in order.

        An uncomputed transition state (value ``None``) takes no column.
        """
        known = [(s, v) for s, v in zip(self.states, self.values)
                 if v is not None]
        cols = self.cols or range(self.start, self.start + len(known))
        return [(c, s, v) for c, (s, v) in zip(cols, known)]

    def humps(self) -> set[int]:
        """Indices ``k`` into :meth:`points` with an uncomputed transition
        state between point ``k`` and point ``k + 1``."""
        out, k = set(), -1
        for v in self.values:
            if v is None:
                out.add(k)
            else:
                k += 1
        return out


def check(pw: Pathway, where: str) -> None:
    """Raise ``ValueError`` for a pathway that cannot be drawn."""
    vals = pw.values
    if vals[0] is None or vals[-1] is None:
        raise ValueError(f"{where}: an uncomputed transition state needs a "
                         f"computed state on both sides")
    if sum(v is not None for v in vals) < 2:
        # A live run gave every state its own one-point series: two bars in
        # one column, nothing joined. States go in x of ONE series.
        raise ValueError(f"{where} has {sum(v is not None for v in vals)} "
                         f"computed state(s); a pathway joins at least two "
                         f"-- the states of one route go in x/y of ONE "
                         f"series")
    for a, b in zip(vals, vals[1:]):
        if a is None and b is None:
            raise ValueError(f"{where}: two uncomputed transition states in "
                             f"a row; put the intermediate between them")


# --------------------------------------------------------------------------
# energies from ORCA outputs

#: Hartree in the units an energy profile is drawn in (CODATA 2018).
HARTREE = {"kcal/mol": 627.5094740631, "kj/mol": 2625.4996394799,
           "ev": 27.211386245988}

_FREQ_BLOCK = "VIBRATIONAL FREQUENCIES"


def imaginary_modes(text: str) -> int | None:
    """Imaginary modes in the LAST frequency block, or ``None`` without one.

    Only the last block: a restarted or repeated frequency run prints one
    block per pass, and counting the whole file counts a TS's mode twice.
    """
    at = text.rfind(_FREQ_BLOCK)
    if at < 0:
        return None
    block = text[at:]
    end = block.find("NORMAL MODES")
    if end > 0:
        block = block[:end]
    return block.count("***imaginary mode***")


def output_energy(text: str, path: str, energy: str) -> float | None:
    """The ``gibbs`` or ``electronic`` energy (Eh) an ORCA output reports,
    read by DELFIN's own parser (:mod:`delfin.energies`)."""
    from delfin import energies
    if energy == "gibbs":
        return energies.gibbs_energy_in_text(text, path)
    return energies.electronic_energy_in_text(text, path)


_COORDS = "CARTESIAN COORDINATES (ANGSTROEM)"
_CHARGE = re.compile(r"Total\s+Charge\s+Charge\s+\.+\s+(-?\d+)")


def composition(text: str) -> Counter | None:
    """Element counts of the LAST geometry an ORCA output prints."""
    at = text.rfind(_COORDS)
    if at < 0:
        return None
    lines = text[at:].splitlines()[2:]
    counts: Counter = Counter()
    for line in lines:
        cells = line.split()
        if len(cells) != 4:
            break
        element = cells[0].rstrip(">:").capitalize()
        if not element.isalpha():
            break
        counts[element] += 1
    return counts or None


def total_charge(text: str) -> int | None:
    m = _CHARGE.search(text)
    return int(m.group(1)) if m else None


def formula(counts: Counter, charge: int | None = None) -> str:
    """``C1 Co1 O2`` in Hill order, with the charge when there is one."""
    def rank(e):      # Hill: C, then H when there is C, then alphabetical
        return (0 if e == "C" else 1 if e == "H" and "C" in counts else 2, e)
    order = sorted(counts, key=rank)
    text = " ".join(f"{e}{counts[e]}" for e in order if counts[e])
    if charge:
        text += f" ({charge:+d})"
    return text or "nothing"


def balance_difference(state: tuple[Counter, int], ref: tuple[Counter, int]) -> str:
    """What ``state`` has more (+) or less (-) than ``ref``; empty if equal."""
    (c, q), (rc, rq) = state, ref
    diff = []
    for e in sorted(set(c) | set(rc)):
        d = c[e] - rc[e]
        if d:
            diff.append(f"{d:+d} {e}")
    if q != rq:
        diff.append(f"charge {q:+d} vs {rq:+d}")
    return ", ".join(diff)


#: G and E of one output differ by the thermal correction -- under an Eh
#: even for a few hundred atoms. More than this is a unit mix-up: cclib
#: 1.8.1 hands xtb's free energy over in eV, 27 times the number in Eh.
_G_MINUS_E_MAX = 2.0

#: Lines a program prints that cclib does not hand over, read in Eh as the
#: program prints them. Measured against cclib 1.8.1 on real outputs: xtb
#: 6.7's free energy comes back in eV, Psi4 1.7's energies not at all. The
#: NWChem line is its standard SCF/DFT summary; cclib's NWChem 7 reference
#: output (a frequency restart) prints no energy, so it is unverified.
_FREE_LINES = (
    re.compile(r"TOTAL FREE ENERGY\s+(-?\d+\.\d+)\s+Eh"),               # xtb
    re.compile(r"Total G, Free enthalpy at\s+[\d.]+\s+\[K\]\s+"
               r"(-?\d+\.\d+)\s+\[Eh\]"),                               # Psi4
)
_ELECTRONIC_LINES = (
    re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh"),                    # xtb
    re.compile(r"Total Energy =\s+(-?\d+\.\d+)"),                        # Psi4
    re.compile(r"Total (?:DFT|SCF) energy =\s+(-?\d+\.\d+)"),            # NWChem
)


def _last(patterns, text: str) -> float | None:
    for pat in patterns:
        found = pat.findall(text)
        if found:
            return float(found[-1])
    return None


def is_orca(text: str) -> bool:
    return "O   R   C   A" in text[:20000]


def has_free_energy(text: str, path: str) -> bool:
    """Whether an output carries a free energy, i.e. a frequency run.

    Read the same way ``read_output`` reads it for ``energy: gibbs``: ORCA
    through DELFIN's parser, xtb and Psi4 by their printed line, anything
    else through cclib's ``freeenergy``. An output cclib cannot read has
    none; ``read_output`` then names the file when it is drawn.
    """
    if is_orca(text):
        return output_energy(text, path, "gibbs") is not None
    if _last(_FREE_LINES, text) is not None:
        return True
    try:
        import io
        import logging

        import cclib
        data = cclib.io.ccread(io.StringIO(text), quiet=True,
                               loglevel=logging.CRITICAL)
    except Exception:
        return False
    return getattr(data, "freeenergy", None) is not None


def read_output(text: str, path: str, energy: str) -> dict:
    """``{"eh", "imag", "atoms", "charge", "program"}`` of one output.

    ORCA through DELFIN's own parser (cclib 1.8.1 stops with an IndexError
    on ORCA 6 outputs); every other program through cclib, with the lines
    cclib misses or mis-scales read as the program prints them. ``eh`` is in
    Hartree; ``atoms``/``charge``/``imag`` are ``None`` where the output
    does not say. Raises ``ValueError`` naming the file when there is no
    energy to draw -- a state without its energy is never guessed.
    """
    if is_orca(text):
        eh = output_energy(text, path, energy)
        if eh is None:
            raise ValueError(
                f"{path}: no "
                + ("'Final Gibbs free energy' (was it run with Freq?)"
                   if energy == "gibbs" else "'FINAL SINGLE POINT ENERGY'"))
        return {"eh": eh, "imag": imaginary_modes(text),
                "atoms": composition(text), "charge": total_charge(text),
                "program": "ORCA"}
    return _read_with_cclib(text, path, energy)


def composite_output(elec_text: str, elec_path: str,
                     freq_text: str, freq_path: str) -> dict:
    """``{"eh", "imag", "atoms", "charge", "program"}`` of one composite state.

    The state's Gibbs energy is built from two levels of theory, per the
    conformer/rotamer protocol::

        G = E(elec output) + (G - E)(freq output)

    * ``elec_text`` / ``elec_path`` -- the high-level single point (e.g.
      DLPNO-CCSD(T)/CBS) that carries only an electronic energy;
    * ``freq_text`` / ``freq_path`` -- the low-level ``Opt Freq`` output
      (e.g. r2SCAN-3c) that carries BOTH a final electronic energy and a
      Final Gibbs free energy, so ``G - E`` of that one output is the thermal
      correction (ZPE + vib/rot/transl + the small 1 atm PV term).

    Both levels must give their number, named in the refusal. The thermal
    correction ``G - E`` is sanity-checked against the same unit-guard as a
    plain Gibbs read (ORCA 6.1 manual, "Vibrational frequencies and
    thermochemistry": free energies and SCF energies are both in Eh and the
    difference is the thermal correction, a fraction of an Eh). ``imag``,
    ``atoms`` and ``charge`` come from the frequency output, which is the one
    that defines the stationary-point character and the stoichiometry.
    """
    from delfin import energies

    e_hi = energies.electronic_energy_in_text(elec_text, elec_path)
    if e_hi is None:
        raise ValueError(
            f"{elec_path}: no 'FINAL SINGLE POINT ENERGY' (the high-level "
            f"electronic energy of the composite); draw energy: electronic "
            f"or pass the DLPNO single-point output")
    g_lo = energies.gibbs_energy_in_text(freq_text, freq_path)
    e_lo = energies.electronic_energy_in_text(freq_text, freq_path)
    if g_lo is None or e_lo is None:
        raise ValueError(
            f"{freq_path}: composite needs BOTH a Final Gibbs free energy "
            f"and a final electronic energy in one output (run it with "
            f"Freq), got g={g_lo}, e={e_lo}")
    # A pair from two different molecules would add a correction to the
    # wrong energy without any sign of it.
    pairs = ((composition(elec_text), composition(freq_text)),
             (total_charge(elec_text), total_charge(freq_text)))
    if any(a is not None and b is not None and a != b for a, b in pairs):
        raise ValueError(
            f"{elec_path} and {freq_path} are not the same species "
            f"(composition or charge differ); a composite pairs the "
            f"single point and the freq output of one structure")
    corr = g_lo - e_lo
    if abs(corr) > _G_MINUS_E_MAX:
        raise ValueError(
            f"{freq_path}: thermal correction G - E = {corr:.3f} Eh is "
            f"larger than any real one -- a unit mix-up, not a correction; "
            f"nothing is drawn from it")
    return {"eh": e_hi + (g_lo - e_lo),
            "imag": imaginary_modes(freq_text),
            "atoms": composition(freq_text),
            "charge": total_charge(freq_text),
            "program": "ORCA"}


def _read_with_cclib(text: str, path: str, energy: str) -> dict:
    try:
        import cclib
    except ImportError:
        raise ValueError(
            f"{path} is not an ORCA output; reading other programs needs "
            f"cclib, a dependency of DELFIN that is not installed in this "
            f"environment (pip install 'cclib>=1.8.1,<1.9')") from None
    import io
    import logging
    try:
        data = cclib.io.ccread(io.StringIO(text), quiet=True,
                               loglevel=logging.CRITICAL)
    except Exception as exc:
        raise ValueError(f"{path}: cclib could not read it "
                         f"({type(exc).__name__}: {exc})") from None
    if data is None:
        raise ValueError(f"{path}: not an output of ORCA or of a program "
                         f"cclib reads (Gaussian, xtb, TURBOMOLE, Q-Chem, "
                         f"Psi4, GAMESS, NWChem, Molpro, ...)")
    program = (getattr(data, "metadata", {}) or {}).get("package") or "?"

    electronic = None
    for attr in ("ccenergies", "mpenergies", "scfenergies"):
        values = getattr(data, attr, None)
        if values is not None and len(values):
            last = values[-1]
            if hasattr(last, "__len__"):        # mpenergies: one row per order
                last = last[-1]
            # cclib's own constant: the line in the output, to the digit.
            from cclib.parser.utils import convertor
            electronic = float(convertor(float(last), "eV", "hartree"))
            break
    if electronic is None:
        electronic = _last(_ELECTRONIC_LINES, text)

    if energy == "gibbs":
        free = _last(_FREE_LINES, text)
        if free is None and getattr(data, "freeenergy", None) is not None:
            free = float(data.freeenergy)
        if free is None:
            raise ValueError(f"{path} ({program}): no free energy in this "
                             f"output (was it run with frequencies?) -- or "
                             f"draw energy: electronic")
        if electronic is not None and abs(free - electronic) > _G_MINUS_E_MAX:
            raise ValueError(
                f"{path} ({program}): free energy {free:.6f} and electronic "
                f"energy {electronic:.6f} Eh are {abs(free - electronic):.1f} "
                f"Eh apart -- a unit mix-up, not a thermal correction; "
                f"nothing is drawn from it")
        eh = free
    else:
        if electronic is None:
            raise ValueError(f"{path} ({program}): no electronic energy found")
        eh = electronic

    atoms = None
    nos = getattr(data, "atomnos", None)
    if nos is not None and len(nos):
        from cclib.parser.utils import PeriodicTable
        table = PeriodicTable()
        atoms = Counter(table.element[int(n)] for n in nos)
    charge = getattr(data, "charge", None)
    freqs = getattr(data, "vibfreqs", None)
    imag = (sum(1 for v in freqs if float(v) < -1.0)
            if freqs is not None else None)
    return {"eh": eh, "imag": imag, "atoms": atoms,
            "charge": None if charge is None else int(charge),
            "program": program}


def is_ts_name(state: str) -> bool:
    """Whether the user named this state a transition state (``TS1``, ``ts``)."""
    return state.strip().upper().startswith("TS")


# --------------------------------------------------------------------------
# the group's CSV

_UNKNOWN = frozenset({"?", "ts?", "ts"})


def _delimiter(rows: list[str]) -> str:
    text = "\n".join(rows)
    if ";" in text:
        return ";"
    if "\t" in text:
        return "\t"
    return ","


def _value(cell: str, delim: str) -> float | None:
    cell = cell.strip()
    if cell.lower() in _UNKNOWN:
        return None
    if delim != ",":
        cell = cell.replace(",", ".")
    return float(cell)


def read_csv(text: str) -> list[Pathway]:
    """Pathways from the group's ``Index;deltaG`` file.

    The format the HPES scripts read: one ``state;value`` row per state, a
    header row, Excel's comma decimals (``2,62``), an optional BOM. An
    optional third column names the state (``3;-12,54;INT1``), and the value
    ``?`` marks an uncomputed transition state between its neighbours.
    Pathways are separated by an empty row (``;`` or a blank line). A file
    without separators starts a new pathway where the first state comes back
    with rows still after it; the first state as the LAST row closes a cycle
    (``1,2,3,4,1``) and stays in its pathway.

    Raises ``ValueError`` with the offending row for anything else.
    """
    lines = [ln.lstrip("﻿").rstrip("\r") for ln in text.splitlines()]
    delim = _delimiter(lines)
    rows: list[tuple[str, float | None] | None] = []
    for ln in lines:
        cells = [c.strip() for c in ln.split(delim)]
        if not any(cells):
            rows.append(None)
            continue
        if len(cells) not in (2, 3) or not cells[1]:
            raise ValueError(f"expected 'state{delim}value[{delim}name]', "
                             f"got {ln!r}")
        try:
            value = _value(cells[1], delim)
        except ValueError:
            if not rows or all(r is None for r in rows):
                continue                       # the header row
            raise ValueError(f"not a number in row {ln!r}") from None
        name = cells[2] if len(cells) == 3 and cells[2] else cells[0]
        if not name:
            raise ValueError(f"row {ln!r} names no state")
        rows.append((name, value))

    groups: list[list[tuple[str, float | None]]] = [[]]
    separated = any(r is None for r in _trimmed(rows))
    data = [r for r in rows if r is not None]
    for row in rows:
        if row is None:
            if groups[-1]:
                groups.append([])
            continue
        group = groups[-1]
        if (not separated and group and row[0] == group[0][0]
                and row is not data[-1]):
            groups.append([row])
            continue
        group.append(row)
    pathways = []
    for n, group in enumerate(g for g in groups if g):
        first = group[0][0]
        pw = Pathway(label="", states=tuple(s for s, _ in group),
                     values=tuple(v for _, v in group),
                     start=int(first) if first.isdigit() else 1)
        check(pw, f"pathway {n + 1}")
        pathways.append(pw)
    if not pathways:
        raise ValueError("no state;value rows found")
    return pathways


def csv_quantity(text: str) -> str | None:
    """The quantity the CSV's header names (``Index;deltaG`` -> ``ΔG``)."""
    for line in text.splitlines():
        line = line.lstrip("\ufeff").strip()
        if not line:
            continue
        cells = re.split(r"[;\t,]", line)
        if len(cells) < 2:
            return None
        head = cells[1].strip().lower().replace("δ", "delta")
        for prefix in ("delta", "d"):
            if head.startswith(prefix) and len(head) > len(prefix):
                head = head[len(prefix):]
                break
        return {"g": "ΔG", "e": "ΔE", "h": "ΔH"}.get(head.strip(" _"))
    return None


def _trimmed(rows):
    """``rows`` without leading and trailing empty rows."""
    lo, hi = 0, len(rows)
    while lo < hi and rows[lo] is None:
        lo += 1
    while hi > lo and rows[hi - 1] is None:
        hi -= 1
    return rows[lo:hi]


# --------------------------------------------------------------------------
# drawing


def placed(pathways: list[Pathway]) -> list[Pathway]:
    """``pathways`` with every state lined up on the states it shares.

    The first state whose name an earlier pathway placed fixes where the
    pathway begins. After that each state takes the next column -- or, if
    an earlier pathway placed its name further right, that column, so a
    shorter route to the same product ends under the product (the dashes
    then span the columns it skips). A name placed twice (a cycle back to
    ``1``) keeps its first column.
    """
    column: dict[str, int] = {}
    out = []
    for pw in pathways:
        names = [s for _, s, _ in Pathway(
            pw.label, pw.states, pw.values, 0).points()]
        start = pw.start
        for i, state in enumerate(names):
            if state in column and column[state] - i >= 0:
                start = column[state] - i
                break
        cols: list[int] = []
        for i, state in enumerate(names):
            if not cols:
                cols.append(start)
                continue
            known = column.get(state)
            cols.append(known if known is not None and known > cols[-1]
                        else cols[-1] + 1)
        pw = Pathway(pw.label, pw.states, pw.values, start, tuple(cols))
        for pos, state, _ in pw.points():
            column.setdefault(state, pos)
        out.append(pw)
    return out


def columns(pathways: list[Pathway]) -> int:
    """How many state columns the profile spans."""
    pos = [p for pw in placed(pathways) for p, _, _ in pw.points()]
    return max(pos) - min(pos) + 1


def figsize(pathways: list[Pathway]) -> tuple[float, float]:
    """10 x 6 in for a short profile, wider per column for a long one."""
    width = _PER_COLUMN * columns(pathways) + 1.5
    return (min(max(width, _MIN_WIDTH), _MAX_WIDTH), _HEIGHT)


def _colours(n: int):
    from matplotlib import colormaps
    tab10 = colormaps["tab10"].colors
    return [tab10[i % len(tab10)] for i in range(n)]


def _hump(x0, y0, x1, y1, top, n=40):
    """Two half-parabolas meeting flat at ``top`` halfway between."""
    xm = (x0 + x1) / 2
    xs, ys = [], []
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        if x <= xm:
            y = top - (top - y0) * ((x - xm) / (x0 - xm)) ** 2
        else:
            y = top - (top - y1) * ((x - xm) / (x1 - xm)) ** 2
        xs.append(x)
        ys.append(y)
    return xs, ys


def draw(fig, pathways: list[Pathway], *, ylabel: str = "ΔG",
         units: str = "kcal/mol", xlabel: str = "", title: str = "") -> None:
    """Draw ``pathways`` onto the empty figure ``fig``."""
    from matplotlib.lines import Line2D

    pathways = placed(pathways)
    ax = fig.subplots()
    colours = _colours(len(pathways))

    n_cols = columns(pathways)
    col_in = (fig.get_figwidth() - 1.5) / n_cols
    scale = min(1.0, col_in / _ROOMY)
    value_fs = max(7.0, 12 * scale)
    tick_fs = max(8.0, 14 * scale)

    # Which pathways reach each (position, value): one bar, split.
    bars: dict[tuple[int, float], list[int]] = {}
    for k, pw in enumerate(pathways):
        for pos, _, v in pw.points():
            key = next((b for b in bars
                        if b[0] == pos and abs(b[1] - v) < _SAME), (pos, v))
            if k not in bars.setdefault(key, []):
                bars[key].append(k)

    values = [v for pw in pathways for _, _, v in pw.points()]
    lo, hi = min(values), max(values)
    span = max(hi - lo, 1.0)
    rise = _HUMP * span
    tops = [max(a[2], b[2]) + rise for pw in pathways
            for k, (a, b) in enumerate(zip(pw.points(), pw.points()[1:]))
            if k in pw.humps()]
    ax.set_ylim(lo - 0.14 * span, max([hi] + tops) + 0.10 * span)

    half = BAR_WIDTH / 2
    for (pos, v), owners in bars.items():
        part = BAR_WIDTH / len(owners)
        for j, k in enumerate(owners):
            x0 = pos - half + j * part
            ax.plot([x0, x0 + part], [v, v], color=colours[k], lw=BAR_LW,
                    solid_capstyle="butt", zorder=3)

    for k, pw in enumerate(pathways):
        pts, humps = pw.points(), pw.humps()
        for i, ((p0, _, v0), (p1, _, v1)) in enumerate(zip(pts, pts[1:])):
            if i in humps:
                xs, ys = _hump(p0 + half, v0, p1 - half, v1,
                               max(v0, v1) + rise)
            else:
                xs, ys = [p0 + half, p1 - half], [v0, v1]
            ax.plot(xs, ys, ls="--", color=colours[k], lw=JOIN_LW, zorder=2)

    # One label per bar, under it -- above it when the next bar down at the
    # same state is so close that the label would sit on that bar.
    by_pos: dict[int, list[float]] = {}
    for pos, v in bars:
        by_pos.setdefault(pos, []).append(v)
    for pos, vs in by_pos.items():
        vs.sort()
        for i, v in enumerate(vs):
            crowded = i > 0 and (v - vs[i - 1]) < _CROWDED * span
            ax.annotate(_value_text(v), (pos, v),
                        xytext=(0, 7 if crowded else -7),
                        textcoords="offset points", ha="center",
                        va="bottom" if crowded else "top",
                        fontsize=value_fs, zorder=4)

    ticks: dict[int, list[str]] = {}
    for pw in pathways:
        for pos, state, _ in pw.points():
            names = ticks.setdefault(pos, [])
            if state not in names:
                names.append(state)
    positions = sorted(ticks)
    # A name wider than its column breaks at its spaces; if names still
    # touch -- measured, not estimated: bold serif is too uneven to count
    # characters -- they all turn, rather than running into each other.
    fits = max(1, int(col_in * 72 * 0.9 / (tick_fs * 0.6)))
    joined = ["/".join(ticks[p]) for p in positions]
    ax.set_xticks(positions, [textwrap.fill(n, fits, break_long_words=False)
                              for n in joined])
    ax.set_xlim(positions[0] - 0.5, positions[-1] + 0.5)
    ax.tick_params(axis="x", bottom=False, labelsize=tick_fs)
    for t in ax.get_xticklabels():
        t.set_fontweight("bold")
    if _touching(fig, ax.get_xticklabels()):
        ax.set_xticks(positions, joined)
        ax.tick_params(axis="x", labelrotation=45)
        for t in ax.get_xticklabels():
            t.set_fontweight("bold")
            t.set_ha("right")
            t.set_rotation_mode("anchor")
    ax.spines[["right", "bottom", "top"]].set_visible(False)

    ax.set_ylabel(_ylabel(ylabel, units), fontsize=14, fontweight="bold")
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=14, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=16, fontweight="bold", pad=12)
    # A legend tells pathways apart; one pathway has nothing to tell apart.
    if sum(bool(pw.label) for pw in pathways) >= 2:
        ax.legend([Line2D([], [], color=colours[k], lw=BAR_LW)
                   for k, pw in enumerate(pathways) if pw.label],
                  [pw.label for pw in pathways if pw.label],
                  frameon=False, fontsize=12)
    fig.tight_layout()


def _touching(fig, labels, gap_pt: float = 4.0) -> bool:
    """Whether two neighbouring tick labels come closer than ``gap_pt``."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    canvas = fig.canvas if isinstance(fig.canvas, FigureCanvasAgg) \
        else FigureCanvasAgg(fig)
    renderer = canvas.get_renderer()
    gap = gap_pt * fig.dpi / 72
    boxes = [t.get_window_extent(renderer) for t in labels]
    return any(a.x1 + gap > b.x0 for a, b in zip(boxes, boxes[1:]))


def _ylabel(quantity: str, units: str) -> str:
    """The y label in quantity calculus (IUPAC Green Book, ISO 80000):
    ``ΔG / kcal mol⁻¹`` -- the quantity symbol italic, Δ, subscripts and the
    unit upright, the exponent raised. Bold, as the layout sets the label.

    A label that already carries a unit (``ΔG / kJ·mol⁻¹``, ``ΔG (kcal/mol)``)
    keeps ITS unit; a live run passed one and got ``ΔG / kJ·mol⁻¹
    (kcal/mol)``, two units on one axis.
    """
    q = quantity.strip()
    for sep in (" / ", "(", "[", "/"):
        if sep in q:
            q, _, given = q.partition(sep)
            q, units = q.strip(), given.strip().rstrip(")]").strip()
            break
    q = _quantity_tex(q)
    u = _unit_tex(units) if units else ""
    return f"{q} / {u}" if u else q


def _plain(text: str) -> str:
    return text.replace("$", r"\$")


def _quantity_tex(q: str) -> str:
    """``ΔG`` -> bold upright Δ + bold italic G; ``ΔG‡``, ``ΔE_rel`` too.

    Anything that is not a symbol (``Relative energy``) stays plain text.
    """
    import re
    m = re.fullmatch(r"(Δ|\\Delta|Delta)?\s*([A-Za-z])\s*(‡|\^?\\ddagger)?"
                     r"\s*(?:_\{?([A-Za-z0-9,]+)\}?)?", q)
    if not m:
        return _plain(q)
    delta, sym, dagger, sub = m.groups()
    tex = (r"\mathbf{\Delta}" if delta else "") + rf"\mathbfit{{{sym}}}"
    if dagger:
        tex += r"^{\mathbf{\ddagger}}"
    if sub:
        tex += rf"_{{\mathbf{{{sub}}}}}"
    return f"${tex}$"


_UNITS = {
    "kcal/mol": r"kcal mol$^{\mathbf{-1}}$",
    "kj/mol": r"kJ mol$^{\mathbf{-1}}$",
    "ev": "eV",
    "eh": r"$\mathbfit{E}_{\mathbf{h}}$",
    "hartree": r"$\mathbfit{E}_{\mathbf{h}}$",
    "cm-1": r"cm$^{\mathbf{-1}}$",
}


def _unit_tex(unit: str) -> str:
    """``kcal/mol``, ``kcal mol-1``, ``kcal·mol⁻¹`` -> ``kcal mol⁻¹``."""
    key = (unit.lower().replace("·", "").replace(" ", "").replace("⁻¹", "-1")
           .replace("^-1", "-1").replace("^{-1}", "-1"))
    key = key.replace("mol-1", "/mol") if key != "cm-1" else key
    return _UNITS.get(key, _plain(unit))


def _value_text(v: float) -> str:
    """Two decimals with a real minus sign (U+2212), never ``−0.00``."""
    text = f"{v:.2f}"
    if text == "-0.00":
        return "0.00"
    return text.replace("-", "−")


#: rcParams of the layout; applied with ``rc_context`` so the caller's
#: process (the dashboard kernel) keeps its own.
RC = {"font.family": "serif", "font.size": 11, "ytick.labelsize": 11,
      "axes.unicode_minus": True, "mathtext.fontset": "dejavuserif"}
