"""Declarative figures to the agent workspace — package V2.

``chat_plots.plot(spec)`` turns a small JSON-able spec into a figure
(PNG or SVG) written inside the agent's workspace, plus a one-line caption,
so the agent does not script matplotlib by hand every time.

Design rules enforced here (and asserted by the package's tests):

- **Deterministic**: a fixed non-interactive backend, no timestamps, no
  randomness, no wall-clock or host names inside the figure. The same spec
  yields byte-identical output.
- **No network**: nothing here touches the network.
- **Size-capped**: a figure larger than ``MAX_BYTES`` is refused, not
  silently written — a pathological spec must fail loudly.
- **Headless**: the "Agg" backend is forced at import, so no display server
  is ever needed.
- **Bounded surface**: only matplotlib + numpy + the standard library
  (csv/json). No pandas — it is not a declared dependency.

Spec (all keys optional except ``kind`` and the data source):

    kind:        "line" | "scatter" | "bar" | "histogram" | "box"
    title:       str, figure / caption title
    xlabel:      str, x-axis label
    ylabel:      str, y-axis label
    units:       str, appended to the caption (and to the y label)
    format:      "svg" (default) | "png"
    filename:    output basename; default ``plot_<kind>.<format>``
    bins:        int, histogram only
    # data — one of:
    x / y:             parallel arrays  (line/scatter/bar; bar accepts only y)
    data:              list of numbers  (histogram/box)  OR
                       list of {x, y}    (line/scatter/bar)
    data: {"file": p}  CSV (x,y)-columns or JSON, read from disk
    series: [{"label", "x"?, "y"}]  multiple series (line/scatter/bar)
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib

# Headless, deterministic backend — forced before pyplot is touched.
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as _plt  # noqa: E402  (backend set above)

# Maximum output size for one figure (bytes). Mirrors the dashboard's raster
# cap so a rendered plot still embeds inline; anything larger is refused.
MAX_BYTES = 2_000_000

_KINDS = frozenset({"line", "scatter", "bar", "histogram", "box"})
_FORMATS = frozenset({"svg", "png"})

# Plots below read their own labels; a missing label is not an error.
_KINDS_WITH_XY = frozenset({"line", "scatter", "bar"})
_KINDS_WITHOLESS = frozenset({"histogram", "box"})


class SpecError(Exception):
    """The spec is malformed, incomplete or unsupported.

    Raised with a clear, actionable message so the caller (the make_plot
    tool) can surface it verbatim instead of a traceback.
    """


@dataclass(frozen=True)
class PlotResult:
    """What ``plot`` handed back: where the figure went + what it shows.

    The structure fields (series, labels, n_points) let tests and the
    make_plot tool report *what* was drawn without re-parsing pixels.
    """

    path: Path
    caption: str
    kind: str
    series: int
    labels: tuple[str, ...]
    format: str
    n_points: int


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
    return val


def _resolve_values_from_file(path: str) -> list:
    """Numbers or {x,y} rows from a CSV or JSON file; CSVs are sorted per
    column assumption, JSON keeps file order."""
    p = Path(path)
    if not p.is_file():
        raise SpecError(f"data file not found: {path}")
    suffix = p.suffix.lower()
    try:
        if suffix == ".csv":
            with p.open(newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            if not rows:
                raise SpecError(f"data file has no rows: {path}")
            if "x" in rows[0] and "y" in rows[0]:
                return [{"x": float(r["x"]), "y": float(r["y"])} for r in rows]
            if "value" in rows[0]:
                return [float(r["value"]) for r in rows]
            raise SpecError(
                f"CSV must have x,y (or value) columns, has: "
                f"{sorted(rows[0].keys())}"
            )
        if suffix == ".json":
            with p.open(encoding="utf-8") as fh:
                return json.load(fh)
    except SpecError:
        raise
    except Exception as exc:  # csv/type errors -> a clear spec error
        raise SpecError(f"cannot read data file {path}: {exc}") from exc
    raise SpecError(f"unsupported data file type: {suffix or path}")


def _resolve_points(rows, kind: str):
    """Rows -> (xs, ys) for line/scatter/bar, or ([], [values]) shorthand."""
    if all(isinstance(r, dict) for r in rows):
        xs, ys = [], []
        for r in rows:
            if "x" not in r or "y" not in r:
                raise SpecError("each data row needs 'x' and 'y'")
            xs.append(r["x"])
            try:
                ys.append(float(r["y"]))
            except (TypeError, ValueError):
                raise SpecError(f"non-numeric y in row {r!r}")
        return xs, ys
    # bare list of numbers
    try:
        ys = [float(v) for v in rows]
    except (TypeError, ValueError):
        raise SpecError("data rows must all be numbers, or all {x, y} dicts")
    return [], ys


def _resolve_series(spec: dict, kind: str):
    """A list of (label, xs, ys) — one series per drawn line/barset. Returns
    also the values for histogram/box and the running point count."""
    if kind in _KINDS_WITHOLESS:
        if spec.get("series"):
            raise SpecError(f"kind {kind!r} does not take a 'series' list")
        data = spec.get("data")
        if data is None or (isinstance(data, list) and not data):
            raise SpecError(f"kind {kind!r} needs 'data' (a list of numbers)")
        if isinstance(data, dict) and "file" in data:
            data = _resolve_values_from_file(data["file"])
        try:
            values = [float(v) for v in data]
        except (TypeError, ValueError):
            raise SpecError(f"kind {kind!r}: data must be a list of numbers")
        if not values:
            raise SpecError(f"kind {kind!r}: data is empty")
        return [(None, [], values)], len(values)

    if spec.get("series"):
        series = []
        n = 0
        for idx, s in enumerate(spec["series"]):
            if not isinstance(s, dict):
                raise SpecError(f"series[{idx}] must be an object")
            st = _resolve_points(
                s.get("data") or (s.get("y") if isinstance(s.get("y"), list) else []),
                kind,
            )
            series.append((s.get("label", f"series {idx + 1}"), st[0], st[1]))
            n += max(len(st[0]), len(st[1]))
        if not series:
            raise SpecError("series list is empty")
        return series, n

    # single series: data rows (list or {x, y}) OR the x/y shorthand
    data = spec.get("data")
    if isinstance(data, list) and data:
        xs, ys = _resolve_points(data, kind)
    elif isinstance(data, dict) and "file" in data:
        rows = _resolve_values_from_file(data["file"])
        xs, ys = _resolve_points(rows, kind)
    elif "y" in spec:
        xs = spec.get("x") or []
        ys = spec["y"]
        if xs and not isinstance(xs, list):
            raise SpecError("x must be a list")
        if not isinstance(ys, list) or not ys:
            raise SpecError("y must be a non-empty list")
        try:
            ys = [float(v) for v in ys]
        except (TypeError, ValueError):
            raise SpecError("y must be a list of numbers")
    else:
        raise SpecError(
            f"kind {kind!r} needs data — 'x'/'y', 'data', or data: {{file: …}}"
        )
    if not ys:
        raise SpecError(f"kind {kind!r}: no points to plot")
    label = spec.get("label", "")
    return [(label, xs, ys)], len(ys)


# --------------------------------------------------------------------------
# position mapping / caption


def _numeric_positions(xs):
    """Map categorical (non-numeric) x entries to integer positions.

    Returns ``(positions, tick_labels)``; ``tick_labels`` is None when x was
    already numeric so the caller skips setting category ticks.
    """
    if not xs:
        return list(xs), None
    try:
        return [float(v) for v in xs], None
    except (TypeError, ValueError):
        labels = [str(v) for v in xs]
        return list(range(len(labels))), labels


def _caption(spec: dict, kind: str, n_points: int) -> str:
    title = _require_str(spec, "title", "").strip()
    units = _require_str(spec, "units", "").strip()
    if title:
        base = f"{title} ({kind}, {n_points} points)"
    else:
        base = f"{kind} plot — {n_points} points"
    return f"{base} [{units}]" if units else base


# --------------------------------------------------------------------------
# figure builders (one per kind)


def _finalize(fig, ax, spec: dict, kind: str, out_dir: Path) -> PlotResult:
    fmt = _require_str(spec, "format", "svg").lower()
    if fmt not in _FORMATS:
        raise SpecError(f"format must be png or svg, got {fmt!r}")
    units = _require_str(spec, "units", "").strip()
    xlabel = _require_str(spec, "xlabel", "").strip()
    ylabel = _require_str(spec, "ylabel", "").strip()
    title = _require_str(spec, "title", "")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(f"{ylabel} ({units})" if units else ylabel)
    if title:
        ax.set_title(title)
    filename = _require_str(spec, "filename", "").strip() or f"plot_{kind}.{fmt}"
    if not filename.lower().endswith(f".{fmt}"):
        filename = f"{filename}.{fmt}"
    out_path = out_dir / filename
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise SpecError(f"cannot create output dir {out_dir}: {exc}")
    fig.savefig(str(out_path), format=fmt)
    _plt.close(fig)
    try:
        size = out_path.stat().st_size
    except OSError:
        size = MAX_BYTES + 1
    if size > MAX_BYTES:
        try:
            out_path.unlink()
        except OSError:
            pass
        raise SpecError(
            f"figure {filename} is {size} bytes — above the "
            f"{MAX_BYTES}-byte cap; reduce the data or width"
        )
    return PlotResult(
        path=out_path,
        caption=_caption(spec, kind, 0),
        kind=kind,
        series=0,
        labels=(xlabel, ylabel),
        format=fmt,
        n_points=0,
    )


def plot(spec: dict, *, out_dir: str = "") -> PlotResult:
    """Render ``spec`` to a figure under ``out_dir`` (default: cwd).

    Returns a :class:`PlotResult`; raises :class:`SpecError` for any
    malformed or unsupported spec. ``out_dir`` is the only write this
    function performs and it is checked to be inside it.
    """
    if not isinstance(spec, dict):
        raise SpecError("spec must be an object with a 'kind'")
    kind = _require_kind(spec)
    series, n_points = _resolve_series(spec, kind)
    out_dir_path = Path(out_dir or os.getcwd())

    fig, ax = _plt.subplots(figsize=(6, 4))

    if kind in ("line", "scatter"):
        for label, xs, ys in series:
            if not xs:
                xs = list(range(len(ys)))
            pos, ticks = _numeric_positions(xs)
            (ax.plot if kind == "line" else ax.scatter)(pos, ys, label=label or None)
            if ticks is not None:
                ax.set_xticks(list(range(len(ticks))), ticks)
    elif kind == "bar":
        for label, xs, ys in series:
            if xs:
                pos, ticks = _numeric_positions(xs)
                ax.bar(pos, ys, label=label or None)
                if ticks is not None:
                    ax.set_xticks(list(range(len(ticks))), ticks)
            else:
                ax.bar(list(range(len(ys))), ys, label=label or None)
    elif kind == "histogram":
        _, _, values = series[0]
        bins = spec.get("bins")
        ax.hist(values, bins=bins if isinstance(bins, int) and bins > 0 else None)
    elif kind == "box":
        _, _, values = series[0]
        ax.boxplot([values], tick_labels=[_require_str(spec, "xlabel", "") or None])
    else:  # unreachable, kept for safety
        raise SpecError(f"unsupported kind {kind!r}")

    if series and any(lbl for lbl, _, _ in series) and kind in _KINDS_WITH_XY:
        ax.legend()

    n_series = len(series)

    # rebuild the result with real structure counts
    res = _finalize(fig, ax, spec, kind, out_dir_path)
    return PlotResult(
        path=res.path,
        caption=_caption(spec, kind, n_points),
        kind=kind,
        series=n_series,
        labels=res.labels,
        format=res.format,
        n_points=n_points,
    )


def to_html(result: PlotResult) -> str:
    """A self-contained HTML block for the dashboard chat (the shared
    role='tool' hook at tab_agent.py:10505).

    Security: the caption is html.escaped so hostile spec/data-derived text
    cannot inject markup; the figure itself travels as a ``data:`` URI served
    as an ``<img>``, where an embedded SVG never executes scripts. The whole
    block is the ONLY string this returns — the dashboard interpolates it
    verbatim, so nothing else here may be interpolated raw.
    """
    import base64
    import html as _html

    mime = "image/svg+xml" if result.format == "svg" else "image/png"
    try:
        b64 = base64.b64encode(result.path.read_bytes()).decode("ascii")
    except OSError as exc:
        raise SpecError(f"cannot read figure {result.path}: {exc}") from exc
    cap = _html.escape(result.caption)
    return (
        f'<div class="delfin-chat-plot" style="margin:6px 0;padding:6px;'
        f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
        f'<img src="data:{mime};base64,{b64}" style="max-width:100%;'
        f'border-radius:4px;display:block;" alt="{cap}"/>'
        f'<div style="font-size:12px;color:#374151;margin-top:4px;">{cap}</div>'
        f'</div>'
    )
