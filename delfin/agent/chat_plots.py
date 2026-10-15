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


def _resolve_pathways(spec: dict, root: Path, reader=None):
    """``([Pathway, ...], n_points)`` for an ``energy_profile`` spec."""
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
            return pathways, sum(len(p.values) for p in pathways)
        try:
            raw = json.loads(text)
        except Exception as exc:
            raise SpecError(f"cannot read {data['file']}: {exc}") from exc
        spec = {"series": raw.get("series") if isinstance(raw, dict) else raw}
    raw = spec.get("series")
    if not raw:
        if "y" not in spec:
            raise SpecError(
                "kind 'energy_profile' needs 'series' (one per pathway: "
                "label, x = state names, y = energies, start), 'x'/'y', or "
                "data: {file: <csv>}")
        raw = [{"x": spec.get("x"), "y": spec["y"],
                "start": spec.get("start", 1)}]
    if not isinstance(raw, list):
        raise SpecError("series must be a list of objects")
    if len(raw) > MAX_SERIES:
        raise SpecError(f"{len(raw)} pathways; at most {MAX_SERIES}")
    pathways, n = [], 0
    for idx, s in enumerate(raw):
        if not isinstance(s, dict):
            raise SpecError(f"series[{idx}] must be an object")
        raw_y = s.get("y")
        if not isinstance(raw_y, list) or not raw_y:
            raise SpecError(f"series[{idx}].y must be a non-empty list of "
                            f"energies (null for an uncomputed transition "
                            f"state)")
        if len(raw_y) > MAX_POINTS:
            raise SpecError(f"series[{idx}].y has {len(raw_y)} values; at "
                            f"most {MAX_POINTS}")
        ys = [None if v is None else _number(v, f"series[{idx}].y")
              for v in raw_y]
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
            raise SpecError(f"series[{idx}]: x has {len(xs)} states, "
                            f"y has {len(ys)}")
        if not xs:
            n_seen = iter(range(start, start + len(ys)))
            xs = ["TS?" if v is None else str(next(n_seen)) for v in ys]
        states = xs
        pw = ep.Pathway(label=_require_str(s, "label", ""),
                        states=tuple(states), values=tuple(ys), start=start)
        try:
            ep.check(pw, f"series[{idx}]")
        except ValueError as exc:
            raise SpecError(str(exc)) from None
        pathways.append(pw)
        n += len(ys)
    if n > MAX_POINTS:
        raise SpecError(f"{n} states over all pathways; at most {MAX_POINTS}")
    return pathways, n


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
    fmt = _require_str(spec, "format", "svg").lower()
    if fmt not in _FORMATS:
        raise SpecError(f"format must be png or svg, got {fmt!r}")
    kind = _require_kind(spec)
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
    pathways, n_points = _resolve_pathways(spec, out_dir_path, reader)
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
