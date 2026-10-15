"""Declarative figures to the agent workspace (the ``make_plot`` tool).

Input: a JSON spec (below). Output: a PNG or SVG written inside the
workspace, and the ``make_plot`` tool result -- the ``DELFIN_CARD:`` marker
followed by ``{"card": "plot", "path", "title", "text"}``. The result
carries no figure bytes and no markup: the dashboard reads the file from
``path`` and shows it with the chat's existing artifact preview
(``tab_agent._make_plot_card_html``); the model and the terminal receive
``text`` only (``chat_viewer.card_model_text``).

Properties:

- Deterministic: the same spec yields byte-identical output (SVG hash salt
  pinned, no date metadata).
- No process-wide state: figures are ``matplotlib.figure.Figure`` objects
  drawn without pyplot, so the calling process's backend and rcParams (the
  dashboard kernel's inline backend) are not touched.
- Bounded: at most ``MAX_POINTS`` values, ``MAX_SERIES`` series,
  ``MAX_BINS`` histogram bins, ``MAX_TEXT`` characters per label, a data
  file of at most ``MAX_DATA_FILE_BYTES``, and a figure of at most
  ``MAX_BYTES``. Non-finite values are refused.
- Nothing in the spec is executed: every key is read as a string, a number
  or a list of them; no key is passed to matplotlib as a keyword argument.

Rejected: an ``<iframe sandbox srcdoc>`` card carrying the base64 figure in
the tool result. A three-point line was a 39 kB result; the dashboard sees
the first 2000 characters of a result, so the card could never be drawn,
and the model was handed the base64 instead of the caption.

Spec (all keys optional except ``kind`` and the data):

    kind:        "line" | "scatter" | "bar" | "histogram" | "box"
                 | "energy_profile" (the group's reaction energy diagram,
                 :mod:`delfin.agent.energy_profile`)
    title, xlabel, ylabel, units: str ("units" is appended to the y label
                 and the caption)
    format:      "svg" (default) | "png"
    filename:    output name inside the workspace; default
                 ``plot_<kind>.<format>``
    bins:        int, histogram only
    label:       legend label of a single series
    # data -- one of:
    x / y:             parallel lists (line/scatter/bar; x optional)
    data:              list of numbers (histogram/box) OR list of {x, y}
    data: {"file": p}  CSV with x,y (or value) columns, or JSON, inside the
                       workspace
    series: [{"label", "x"?, "y"} | {"label", "data"}]  several series

energy_profile: one series per pathway, ``x`` the state names, ``y`` the
energies, ``null`` for a transition state that was not computed (drawn as a
hump); a state named in two pathways is one column, so a branch lines up
by its names (``start`` places only a pathway that shares none); or
``data: {"file": "<name>.csv"}`` in the group's ``Index;deltaG`` format. ``ylabel`` defaults to ``ΔG`` and
``units`` to ``kcal/mol``.
"""

from __future__ import annotations

import csv
import io
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path

#: Maximum size of one written figure (bytes); a larger one is refused
#: before it is written. Equal to the dashboard's raster preview cap.
MAX_BYTES = 2_000_000
#: Total values over all series. matplotlib's time grows with the count
#: and a data file can hold millions of rows.
MAX_POINTS = 100_000
MAX_SERIES = 20
#: An unbounded ``bins`` is a denial of service: 3,000,000 bins ran for
#: more than three CPU minutes without finishing.
MAX_BINS = 1_000
MAX_TEXT = 300
MAX_DATA_FILE_BYTES = 5_000_000
#: An ORCA output read for ``energy_profile``; a frequency run is tens of MB.
MAX_OUTPUT_BYTES = 200_000_000
#: Source lines in a caption; the dashboard reads 2000 characters of it.
_MAX_NOTES = 14

#: Prefix of a card result; shared with ``show_molecule``
#: (``delfin.dashboard.chat_viewer.CARD_MARKER``).
CARD_MARKER = "DELFIN_CARD:"

_KINDS = frozenset({"line", "scatter", "bar", "histogram", "box",
                    "energy_profile"})
_FORMATS = frozenset({"svg", "png"})
_KINDS_WITH_VALUES = frozenset({"histogram", "box"})
_SVG_HASHSALT = "delfin-chat-plots"


class SpecError(Exception):
    """The spec is malformed, incomplete, over a bound or unsupported.

    The message is shown to the model verbatim.
    """


@dataclass(frozen=True)
class PlotResult:
    """The written figure and what it shows."""

    path: Path
    caption: str
    kind: str
    series: int
    labels: tuple[str, ...]
    format: str
    n_points: int


def _is_within(child: Path, parent: Path) -> bool:
    """True when ``child`` resolves to ``parent`` or a path below it."""
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


# --------------------------------------------------------------------------
# spec parsing / data resolution


def _require_kind(spec: dict) -> str:
    kind = spec.get("kind")
    if not isinstance(kind, str) or kind not in _KINDS:
        raise SpecError(
            f"unknown kind {kind!r}: choose one of "
            f"{', '.join(sorted(_KINDS))}"
        )
    return kind


def _require_str(spec: dict, key: str, default: str = "") -> str:
    val = spec.get(key, default)
    if val is None:
        return default
    if not isinstance(val, str):
        raise SpecError(f"{key} must be a string, got {type(val).__name__}")
    if len(val) > MAX_TEXT:
        raise SpecError(f"{key} is {len(val)} characters; at most {MAX_TEXT}")
    return val


def _number(v, what: str) -> float:
    """A finite float from ``v``; bool, NaN and infinity are refused."""
    if isinstance(v, bool):
        raise SpecError(f"{what} must be numbers, got {v!r}")
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise SpecError(f"{what} must be numbers, got {v!r}") from None
    if not math.isfinite(f):
        raise SpecError(f"{what} contains {v!r}; drop or replace non-finite "
                        f"values")
    return f


def _numbers(values, what: str) -> list[float]:
    if not isinstance(values, list):
        raise SpecError(f"{what} must be a list of numbers")
    if len(values) > MAX_POINTS:
        raise SpecError(f"{what} has {len(values)} values; at most "
                        f"{MAX_POINTS}")
    return [_number(v, what) for v in values]


def _categories(xs, what: str = "x") -> list:
    """x values: a list of numbers or of category labels (str/int/float)."""
    if not isinstance(xs, list):
        raise SpecError(f"{what} must be a list")
    if len(xs) > MAX_POINTS:
        raise SpecError(f"{what} has {len(xs)} values; at most {MAX_POINTS}")
    for v in xs:
        if isinstance(v, bool) or not isinstance(v, (str, int, float)):
            raise SpecError(f"{what} entries must be numbers or strings, "
                            f"got {v!r}")
        if isinstance(v, str) and len(v) > MAX_TEXT:
            raise SpecError(f"an {what} label is longer than {MAX_TEXT} "
                            f"characters")
        if isinstance(v, float) and not math.isfinite(v):
            raise SpecError(f"{what} contains {v!r}; drop or replace "
                            f"non-finite values")
    return xs


def _read_data_file(path, root: Path, reader=None) -> tuple[str, str]:
    """``(text, suffix)`` of a CSV or JSON file inside ``root``.

    A path that resolves outside ``root`` (absolute, ``..``, symlink) is
    refused, as is a file over ``MAX_DATA_FILE_BYTES``. ``reader`` -- the
    dispatch layer's read gate (secret-deny list, read roots) -- is called
    with the resolved path and returns the text; without it the file is
    read directly.
    """
    if not isinstance(path, str) or not path:
        raise SpecError("data.file must be a path string")
    root_res = root.resolve()
    p = (root / path if not Path(path).is_absolute() else Path(path)).resolve()
    if not _is_within(p, root_res):
        raise SpecError(f"data file {path!r} resolves outside the workspace "
                        f"{root}")
    if not p.is_file():
        raise SpecError(f"data file not found: {path}")
    size = p.stat().st_size
    if size > MAX_DATA_FILE_BYTES:
        raise SpecError(f"data file {path} is {size} bytes; at most "
                        f"{MAX_DATA_FILE_BYTES}")
    suffix = p.suffix.lower()
    if suffix not in (".csv", ".json"):
        raise SpecError(f"unsupported data file type: {suffix or path}")
    if reader is None:
        def reader(q):
            with q.open(encoding="utf-8") as fh:
                return fh.read()
    return reader(p), suffix


def _resolve_values_from_file(path, root: Path, reader=None):
    """Numbers or ``{x, y}`` rows from a CSV or JSON file inside ``root``."""
    text, suffix = _read_data_file(path, root, reader)
    try:
        if suffix == ".json":
            return json.loads(text)
        rows = list(csv.DictReader(io.StringIO(text)))
    except Exception as exc:
        raise SpecError(f"cannot read data file {path}: {exc}") from exc
    if not rows:
        raise SpecError(f"data file has no rows: {path}")
    if len(rows) > MAX_POINTS:
        raise SpecError(f"data file {path} has {len(rows)} rows; at most "
                        f"{MAX_POINTS}")
    if "x" in rows[0] and "y" in rows[0]:
        return [{"x": r["x"], "y": r["y"]} for r in rows]
    if "value" in rows[0]:
        return [r["value"] for r in rows]
    raise SpecError(f"CSV must have x,y (or value) columns, has: "
                    f"{sorted(k for k in rows[0].keys() if k)}")


def _resolve_points(rows) -> tuple[list, list[float]]:
    """Rows -> (xs, ys): a list of ``{x, y}`` or a bare list of numbers."""
    if not isinstance(rows, list) or not rows:
        raise SpecError("data must be a non-empty list")
    if len(rows) > MAX_POINTS:
        raise SpecError(f"data has {len(rows)} rows; at most {MAX_POINTS}")
    if all(isinstance(r, dict) for r in rows):
        xs, ys = [], []
        for r in rows:
            if "x" not in r or "y" not in r:
                raise SpecError("each data row needs 'x' and 'y'")
            xs.append(r["x"])
            ys.append(_number(r["y"], "y"))
        # CSV cells are strings: numeric ones become numbers.
        try:
            xs = [_number(v, "x") for v in xs]
        except SpecError:
            pass
        return _categories(xs), ys
    return [], _numbers(rows, "data")


def _resolve_series(spec: dict, kind: str, root: Path, reader=None):
    """``([(label, xs, ys), ...], n_points)`` for the spec."""
    if kind in _KINDS_WITH_VALUES:
        if spec.get("series"):
            raise SpecError(f"kind {kind!r} does not take a 'series' list")
        data = spec.get("data")
        if isinstance(data, dict) and "file" in data:
            data = _resolve_values_from_file(data["file"], root, reader)
        if data is None or data == []:
            raise SpecError(f"kind {kind!r} needs 'data' (a list of numbers)")
        values = _numbers(data, f"kind {kind!r}: data")
        return [(None, [], values)], len(values)

    if spec.get("series"):
        raw = spec["series"]
        if not isinstance(raw, list):
            raise SpecError("series must be a list of objects")
        if len(raw) > MAX_SERIES:
            raise SpecError(f"{len(raw)} series; at most {MAX_SERIES}")
        series, n = [], 0
        for idx, s in enumerate(raw):
            if not isinstance(s, dict):
                raise SpecError(f"series[{idx}] must be an object")
            if s.get("data") is not None:
                xs, ys = _resolve_points(s["data"])
            else:
                ys = _numbers(s.get("y"), f"series[{idx}].y")
                xs = _categories(s.get("x") or [], f"series[{idx}].x")
            if not ys:
                raise SpecError(f"series[{idx}] has no points")
            if xs and len(xs) != len(ys):
                raise SpecError(f"series[{idx}]: x has {len(xs)} values, "
                                f"y has {len(ys)}")
            label = _require_str(s, "label", f"series {idx + 1}")
            series.append((label, xs, ys))
            n += len(ys)
        if n > MAX_POINTS:
            raise SpecError(f"{n} points over all series; at most "
                            f"{MAX_POINTS}")
        return series, n

    data = spec.get("data")
    if isinstance(data, dict) and "file" in data:
        xs, ys = _resolve_points(
            _resolve_values_from_file(data["file"], root, reader))
    elif isinstance(data, list) and data:
        xs, ys = _resolve_points(data)
    elif "y" in spec:
        ys = _numbers(spec["y"], "y")
        xs = _categories(spec.get("x") or [])
    else:
        raise SpecError(
            f"kind {kind!r} needs data: 'x'/'y', 'data', or data: {{file: ...}}"
        )
    if not ys:
        raise SpecError(f"kind {kind!r}: no points to plot")
    if xs and len(xs) != len(ys):
        raise SpecError(f"x has {len(xs)} values, y has {len(ys)}")
    return [(_require_str(spec, "label", ""), xs, ys)], len(ys)


def _read_output(path, root: Path, reader=None) -> tuple[str, Path]:
    """``(text, resolved path)`` of an ORCA output for ``energy_profile``.

    Unlike a data file it may lie outside the workspace -- outputs live in
    the calculation folders -- so the dispatch layer's read gate
    (``reader``: secret-deny list, read roots) decides, exactly as it does
    for ``read_file``. Bounded by ``MAX_OUTPUT_BYTES``.
    """
    if not isinstance(path, str) or not path.strip():
        raise SpecError("each files entry must be a path string or null")
    p = Path(path).expanduser()
    p = (p if p.is_absolute() else root / p).resolve()
    if not p.is_file():
        raise SpecError(f"output not found: {path}")
    size = p.stat().st_size
    if size > MAX_OUTPUT_BYTES:
        raise SpecError(f"output {path} is {size} bytes; at most "
                        f"{MAX_OUTPUT_BYTES}")
    if reader is None:
        with p.open(encoding="utf-8", errors="replace") as fh:
            return fh.read(), p
    return reader(p), p


def _resolve_pathways(spec: dict, root: Path, reader=None):
    """``([Pathway, ...], n_points, notes)`` for an ``energy_profile`` spec.

    ``notes`` are lines for the caption: where each energy read from an
    output came from, and every state whose frequencies contradict its name.
    """
    from delfin.agent import energy_profile as ep

    data = spec.get("data")
    for alias in ("data_file", "file"):         # what models write instead
        if data is None and isinstance(spec.get(alias), str):
            data = {"file": spec[alias]}
    if isinstance(data, dict) and "file" in data:
        text, suffix = _read_data_file(data["file"], root, reader)
        if suffix == ".csv":
            try:
                pathways = ep.read_csv(text)
            except ValueError as exc:
                raise SpecError(f"cannot read {data['file']}: {exc}") from exc
            if len(pathways) > MAX_SERIES:
                raise SpecError(f"{len(pathways)} pathways; at most "
                                f"{MAX_SERIES}")
            return pathways, sum(len(p.values) for p in pathways), []
        try:
            raw = json.loads(text)
        except Exception as exc:
            raise SpecError(f"cannot read {data['file']}: {exc}") from exc
        spec = dict(spec, series=(raw.get("series")
                                  if isinstance(raw, dict) else raw))
    raw = spec.get("series")
    if not raw:
        if "y" not in spec and "files" not in spec:
            raise SpecError(
                "kind 'energy_profile' needs 'series' (one per pathway: "
                "label, x = state names, y = energies or files = ORCA "
                "outputs), 'x'/'y', or data: {file: <csv>}")
        raw = [{k: spec[k] for k in ("x", "y", "files", "start")
                if k in spec}]
    if not isinstance(raw, list):
        raise SpecError("series must be a list of objects")
    if len(raw) > MAX_SERIES:
        raise SpecError(f"{len(raw)} pathways; at most {MAX_SERIES}")

    energy = _require_str(spec, "energy", "gibbs").strip().lower()
    if energy not in ("gibbs", "electronic"):
        raise SpecError(f"energy must be gibbs or electronic, got {energy!r}")

    # First pass: values as given, or Hartree read from the outputs.
    parsed, outputs, n = [], {}, 0
    for idx, s in enumerate(raw):
        if not isinstance(s, dict):
            raise SpecError(f"series[{idx}] must be an object")
        from_files = s.get("files") is not None
        key = "files" if from_files else "y"
        entries = s.get(key)
        if not isinstance(entries, list) or not entries:
            raise SpecError(
                f"series[{idx}] needs 'y' (energies, null for an "
                f"uncomputed transition state) or 'files' (ORCA outputs, "
                f"null likewise)")
        if len(entries) > MAX_POINTS:
            raise SpecError(f"series[{idx}].{key} has {len(entries)} "
                            f"entries; at most {MAX_POINTS}")
        if from_files:
            ys = []
            for f in entries:
                if f is None:
                    ys.append(None)
                    continue
                key = _state_key(f, idx)
                if key not in outputs:
                    outputs[key] = _state_record(key, root, reader, energy)
                ys.append(outputs[key]["eh"])
        else:
            ys = [None if v is None else _number(v, f"series[{idx}].y")
                  for v in entries]
        start = s.get("start", 1)
        if isinstance(start, bool) or not isinstance(start, int) \
                or not 0 <= start <= MAX_POINTS:
            raise SpecError(f"series[{idx}].start must be a whole number "
                            f"from 0, got {start!r}")
        xs = [str(v) for v in
              _categories(s.get("x") or [], f"series[{idx}].x")]
        computed = sum(v is not None for v in ys)
        if xs and len(xs) == computed < len(ys):
            # Names for the computed states only: the uncomputed
            # transition states have none to give.
            names = iter(xs)
            xs = ["TS?" if v is None else next(names) for v in ys]
        if xs and len(xs) != len(ys):
            hint = (" -- several outputs that make ONE state (catalyst + "
                    "free CO2) go in a nested list: files: [[\"a.out\", "
                    "\"b.out\"], \"ab.out\"]") if from_files else ""
            raise SpecError(f"series[{idx}]: x has {len(xs)} states, "
                            f"{key} has {len(ys)}{hint}")
        if not xs:
            n_seen = iter(range(start, start + len(ys)))
            xs = ["TS?" if v is None else str(next(n_seen)) for v in ys]
        parsed.append((idx, s, xs, ys, start,
                       [None if f is None else _state_key(f, idx)
                        for f in entries] if from_files else None))
        n += len(ys)
    if n > MAX_POINTS:
        raise SpecError(f"{n} states over all pathways; at most {MAX_POINTS}")

    # Second pass: outputs become energies relative to ONE reference -- the
    # first output of the first pathway that has outputs, or ``reference`` --
    # so a branch is measured from where the route started, not from itself.
    notes: list[str] = []
    ref_eh, factor = None, 1.0
    if outputs:
        units = (_require_str(spec, "units", "kcal/mol").strip()
                 if spec.get("units") is not None else "kcal/mol")
        factor = ep.HARTREE.get(units.lower().replace(" ", ""))
        if factor is None:
            raise SpecError(f"energies from outputs are drawn in kcal/mol, "
                            f"kJ/mol or eV, not {units!r}")
        ref_file = spec.get("reference")
        if ref_file is None:
            ref_file = next(f for *_, files in parsed if files
                            for f in files if f is not None)
        else:
            ref_file = _state_key(ref_file, "reference")
        if ref_file not in outputs:
            outputs[ref_file] = _state_record(ref_file, root, reader, energy)
        ref = outputs[ref_file]
        ref_eh = ref["eh"]
        notes.append(f"{energy} energies relative to {' + '.join(ref_file)} "
                     f"({ref_eh:.6f} Eh), x {factor:g} {units}/Eh")
        _check_balance(parsed, outputs, ref_file, notes)

    pathways = []
    for idx, s, xs, ys, start, files in parsed:
        if files:
            ys = [None if v is None else (v - ref_eh) * factor for v in ys]
            for state, f, v in zip(xs, files, ys):
                if f is None:
                    continue
                rec = outputs[f]
                notes.append(f"{state}: {' + '.join(f)} {rec['eh']:.6f} Eh "
                             f"-> {v:.2f}" + _frequency_note(state, rec["imag"]))
        pw = ep.Pathway(label=_require_str(s, "label", ""),
                        states=tuple(xs), values=tuple(ys), start=start)
        try:
            ep.check(pw, f"series[{idx}]")
        except ValueError as exc:
            raise SpecError(str(exc)) from None
        pathways.append(pw)
    return pathways, n, notes


def _state_key(entry, where) -> tuple[str, ...]:
    """A state's outputs as a tuple: one path, or several fragments whose
    energies add up (``["cat.out", "co2.out"]`` for catalyst + free CO2)."""
    if isinstance(entry, str):
        return (entry,)
    if (isinstance(entry, list) and entry
            and all(isinstance(f, str) and f.strip() for f in entry)):
        return tuple(entry)
    raise SpecError(f"series[{where}]: a files entry is a path, a list of "
                    f"paths (fragments of one state) or null, got {entry!r}")


def _state_record(key: tuple[str, ...], root: Path, reader,
                  energy: str) -> dict:
    """``{"eh", "imag", "atoms", "charge"}`` of one state, its fragments
    summed; a refusal naming the output that lacks its energy."""
    from collections import Counter
    from delfin.agent import energy_profile as ep

    eh, imag, atoms, charge = 0.0, 0, Counter(), 0
    for path in key:
        text, resolved = _read_output(path, root, reader)
        e = ep.output_energy(text, str(resolved), energy)
        if e is None:
            what = ("'Final Gibbs free energy' (was it run with Freq?)"
                    if energy == "gibbs" else "'FINAL SINGLE POINT ENERGY'")
            raise SpecError(f"{path}: no {what}; nothing is drawn for a "
                            f"state without its energy")
        eh += e
        n = ep.imaginary_modes(text)
        imag = None if imag is None or n is None else imag + n
        comp, q = ep.composition(text), ep.total_charge(text)
        if atoms is None or comp is None or q is None or charge is None:
            atoms = charge = None
        else:
            atoms += comp
            charge += q
    return {"eh": eh, "imag": imag, "atoms": atoms, "charge": charge}


def _check_balance(parsed, outputs, ref_key, notes) -> None:
    """Every state holds the same atoms and charge as the reference.

    A profile is one system followed through a reaction: CO2 binding to a
    metal starts at catalyst + FREE CO2, not at the bare catalyst. A state
    missing a fragment is not on the surface at all, and its energy against
    the others means nothing -- so the plot is refused, naming what is short.
    """
    from delfin.agent import energy_profile as ep

    ref = outputs[ref_key]
    if ref["atoms"] is None:
        notes.append("mass balance NOT checked: no geometry or charge in "
                     + " + ".join(ref_key))
        return
    want = (ref["atoms"], ref["charge"])
    wrong, unread = [], []
    for idx, s, xs, ys, start, files in parsed:
        for state, key in zip(xs, files or []):
            if key is None:
                continue
            rec = outputs[key]
            if rec["atoms"] is None:
                unread.append(state)
                continue
            diff = ep.balance_difference((rec["atoms"], rec["charge"]), want)
            if diff:
                wrong.append(f"{state} ({' + '.join(key)}): {diff}")
    if wrong:
        raise SpecError(
            "mass balance: every state must hold the same atoms and charge "
            f"as the reference {' + '.join(ref_key)} "
            f"[{ep.formula(*want)}]. Off: " + "; ".join(wrong) + ". Add "
            "the missing molecules as fragments of that state, e.g. files: "
            "[[\"cat.out\", \"co2.out\"], \"cat-co2.out\"].")
    notes.append("mass balance checked: every state "
                 f"{ep.formula(*want)}" + (
                     f"; NOT checked for {', '.join(unread)} (no geometry)"
                     if unread else ""))


def _frequency_note(state: str, imag) -> str:
    """What the frequencies say about a state's name, when they disagree."""
    from delfin.agent import energy_profile as ep

    if imag is None:
        return " [no frequencies: not checked]" if ep.is_ts_name(state) else ""
    if ep.is_ts_name(state) and imag != 1:
        return f" [WARNING: {imag} imaginary modes, a TS needs exactly 1]"
    if not ep.is_ts_name(state) and imag:
        return f" [WARNING: {imag} imaginary mode(s), not a minimum]"
    return ""


def _positions(xs):
    """``(positions, tick_labels)``; tick_labels is None for numeric x."""
    if all(isinstance(v, (int, float)) for v in xs):
        return [float(v) for v in xs], None
    labels = [str(v) for v in xs]
    return list(range(len(labels))), labels


def _caption(title: str, units: str, kind: str, n_points: int) -> str:
    title, units = title.strip(), units.strip()
    base = (f"{title} ({kind}, {n_points} points)" if title
            else f"{kind} plot, {n_points} points")
    return f"{base} [{units}]" if units else base


# --------------------------------------------------------------------------
# rendering


def resolve_out_path(spec: dict, out_dir: Path) -> Path:
    """The file :func:`plot` writes for ``spec`` under ``out_dir``.

    Shared by the renderer and the dispatch executor, so the path the write
    gate checks is the path written. Refuses a name outside ``out_dir``.
    """
    kind = _require_kind(spec)
    # An energy profile goes on a slide or into a paper draft: PNG at
    # 300 dpi as the group's scripts wrote it. Two live runs asked for a PNG
    # without saying format and then converted the SVG by hand.
    fmt = _require_str(spec, "format",
                       "png" if kind == "energy_profile" else "svg").lower()
    if fmt not in _FORMATS:
        raise SpecError(f"format must be png or svg, got {fmt!r}")
    filename = _require_str(spec, "filename", "").strip() or f"plot_{kind}.{fmt}"
    if not filename.lower().endswith(f".{fmt}"):
        filename = f"{filename}.{fmt}"
    out_path = (out_dir / filename).resolve()
    if not _is_within(out_path, out_dir.resolve()):
        raise SpecError(f"filename {filename!r} resolves outside the workspace "
                        f"{out_dir}")
    return out_path


def plot(spec: dict, *, out_dir: str = "", reader=None) -> PlotResult:
    """Render ``spec`` to a figure under ``out_dir`` (default: cwd).

    Raises :class:`SpecError` for a malformed, unsupported or over-bound
    spec; nothing is written in that case. The only write is the figure
    at :func:`resolve_out_path`.
    """
    from matplotlib.figure import Figure

    if not isinstance(spec, dict):
        raise SpecError("spec must be an object with a 'kind'")
    kind = _require_kind(spec)
    out_dir_path = Path(out_dir or os.getcwd())
    out_path = resolve_out_path(spec, out_dir_path)
    fmt = out_path.suffix[1:].lower()
    title = _require_str(spec, "title", "")
    units = _require_str(spec, "units", "").strip()
    xlabel = _require_str(spec, "xlabel", "").strip()
    ylabel = _require_str(spec, "ylabel", "").strip()
    bins = spec.get("bins")
    if bins is not None and (isinstance(bins, bool)
                             or not isinstance(bins, int)
                             or not 1 <= bins <= MAX_BINS):
        raise SpecError(f"bins must be an integer from 1 to {MAX_BINS}")
    if kind == "energy_profile":
        return _plot_energy_profile(spec, out_path, out_dir_path, reader,
                                    title=title, xlabel=xlabel)
    series, n_points = _resolve_series(spec, kind, out_dir_path, reader)

    fig = Figure(figsize=(6, 4))
    ax = fig.subplots()
    if kind in ("line", "scatter", "bar"):
        draw = {"line": ax.plot, "scatter": ax.scatter, "bar": ax.bar}[kind]
        for label, xs, ys in series:
            pos, ticks = _positions(xs) if xs else (list(range(len(ys))), None)
            draw(pos, ys, label=label or None)
            if ticks is not None:
                ax.set_xticks(pos, ticks)
        if any(lbl for lbl, _, _ in series):
            ax.legend()
    elif kind == "histogram":
        ax.hist(series[0][2], bins=bins if bins is not None else "auto")
    else:  # box
        ax.boxplot([series[0][2]], tick_labels=[xlabel])
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(f"{ylabel} ({units})" if units else ylabel)
    if title:
        ax.set_title(title)

    _write(fig, out_path, fmt)
    return PlotResult(
        path=out_path,
        caption=_caption(title, units, kind, n_points),
        kind=kind,
        series=len(series),
        labels=(xlabel, ylabel),
        format=fmt,
        n_points=n_points,
    )


def _plot_energy_profile(spec, out_path, out_dir_path, reader, *,
                         title, xlabel) -> PlotResult:
    from matplotlib import rc_context
    from matplotlib.figure import Figure
    from delfin.agent import energy_profile as ep

    fmt = out_path.suffix[1:].lower()
    ylabel = _require_str(spec, "ylabel", "ΔG").strip() or "ΔG"
    units = (_require_str(spec, "units", "kcal/mol").strip()
             if spec.get("units") is not None else "kcal/mol")
    pathways, n_points, notes = _resolve_pathways(spec, out_dir_path,
                                                  reader)
    if not notes:
        # Said in the result, where the model reads it at once: a live run
        # grepped the outputs, did the Hartree arithmetic itself and then
        # reported "no imaginary modes" for an adduct that had one.
        notes.append("values drawn as given: mass balance and imaginary "
                     "modes NOT checked -- if they come from ORCA outputs, "
                     "draw again with files: [...] instead of y")
    labels = [p.label for p in pathways if p.label]
    if not title and len(pathways) == 1 and labels:
        # One pathway has no legend to carry its name; a live run put the
        # requested heading there and it vanished.
        title = labels[0]
    with rc_context(ep.RC):
        fig = Figure(figsize=ep.figsize(pathways))
        ep.draw(fig, pathways, ylabel=ylabel, units=units, xlabel=xlabel,
                title=title)
        _write(fig, out_path, fmt, dpi=300)
    highest = max((v, s) for p in pathways for _, s, v in p.points())
    starts = [p.start for p in ep.placed(pathways)]
    caption = (f"{_caption(title, units, 'energy_profile', n_points)}; "
               f"{len(pathways)} pathway(s) starting at column(s) "
               f"{', '.join(map(str, starts))}, highest {highest[0]:.2f} "
               f"{units} at state {highest[1]}")
    if notes:
        # What contradicts the picture first, so a cut caption keeps it.
        notes = ([n for n in notes if "WARNING" in n or "NOT" in n]
                 + [n for n in notes if not ("WARNING" in n or "NOT" in n)])
        shown = notes[:_MAX_NOTES]
        if len(notes) > _MAX_NOTES:
            shown.append(f"... {len(notes) - _MAX_NOTES} more states")
        caption += "\n" + "\n".join(shown)
    return PlotResult(
        path=out_path,
        caption=caption,
        kind="energy_profile",
        series=len(pathways),
        labels=(xlabel, ylabel),
        format=fmt,
        n_points=n_points,
    )


def _write(fig, out_path: Path, fmt: str, dpi=None) -> None:
    """Save ``fig`` to ``out_path`` within ``MAX_BYTES``."""
    from matplotlib import rc_context

    buf = io.BytesIO()
    try:
        with rc_context({"svg.hashsalt": _SVG_HASHSALT}):
            fig.savefig(buf, format=fmt, dpi=dpi,
                        metadata={"Date": None} if fmt == "svg" else {})
    except ValueError as exc:   # a mathtext label it cannot parse, ...
        raise SpecError(f"cannot render the figure: {exc}") from exc
    data = buf.getvalue()
    if len(data) > MAX_BYTES:
        raise SpecError(
            f"figure {out_path.name} would be {len(data)} bytes, above the "
            f"{MAX_BYTES}-byte cap; reduce the data or use format png"
        )
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(data)
    except OSError as exc:
        raise SpecError(f"cannot write {out_path}: {exc}") from exc


def card_result(result: PlotResult) -> str:
    """The ``make_plot`` tool result for a written figure.

    ``DELFIN_CARD:`` + ``{"card": "plot", "path": <resolved path>,
    "title": <file name>, "text": <plain line for the model>}`` -- a few
    hundred characters, so it fits the dashboard's 2000-character view of
    a result.
    """
    return CARD_MARKER + json.dumps({
        "card": "plot",
        "path": str(result.path),
        "title": result.path.name,
        "text": f"Wrote {result.path} -- {result.caption}",
    })


def make_plot(spec: dict, *, out_dir: str = "", reader=None) -> str:
    """Render ``spec`` (:func:`plot`) and return :func:`card_result`."""
    return card_result(plot(spec, out_dir=out_dir, reader=reader))
