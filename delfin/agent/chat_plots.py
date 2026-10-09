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

# Deterministic SVG glyph IDs: matplotlib otherwise salts the internal path
# IDs with a random value drawn from a per-call hashsalt, making byte-identical
# output impossible. Pin the salt so identical specs give identical SVGs.
matplotlib.rcParams["svg.hashsalt"] = "delfin-v2-chat-plots-deterministic"

import matplotlib.pyplot as _plt  # noqa: E402  (backend set above)

# Maximum output size for one figure (bytes). Mirrors the dashboard's raster
# cap so a rendered plot still embeds inline; anything larger is refused.
MAX_BYTES = 2_000_000

_KINDS = frozenset({"line", "scatter", "bar", "histogram", "box"})
_FORMATS = frozenset({"svg", "png"})

# Plots below read their own labels; a missing label is not an error.
_KINDS_WITH_XY = frozenset({"line", "scatter", "bar"})
_KINDS_WITHOLESS = frozenset({"histogram", "box"})


def _is_within(child: Path, parent: Path) -> bool:
    """True when ``child`` (resolved) is ``parent`` (resolved) or below it.

    Guards the workspace boundary for both the written figure and any data
    file a spec points at: nothing may escape out_dir."""
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


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


def _resolve_values_from_file(path: str, root: Path) -> list:
    """Numbers or {x,y} rows from a CSV or JSON file; CSVs are sorted per
    column assumption, JSON keeps file order.

    The file must live inside ``root`` (the workspace/out_dir): an absolute
    or ``..`` path that resolves outside it is refused — plot must never
    read a file the workspace does not own (operator data-is-never-
    instruction contract)."""
    root_res = root.resolve()
    p = (root / path if not Path(path).is_absolute() else Path(path)).resolve()
    if not _is_within(p, root_res):
        raise SpecError(f"data file {path!r} resolves outside the workspace "
                        f"{root}")
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


def _resolve_series(spec: dict, kind: str, root: Path):
    """A list of (label, xs, ys) — one series per drawn line/barset. Returns
    also the values for histogram/box and the running point count."""
    if kind in _KINDS_WITHOLESS:
        if spec.get("series"):
            raise SpecError(f"kind {kind!r} does not take a 'series' list")
        data = spec.get("data")
        if data is None or (isinstance(data, list) and not data):
            raise SpecError(f"kind {kind!r} needs 'data' (a list of numbers)")
        if isinstance(data, dict) and "file" in data:
            data = _resolve_values_from_file(data["file"], root)
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
        rows = _resolve_values_from_file(data["file"], root)
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
    out_path = (out_dir / filename).resolve()
    if not _is_within(out_path, out_dir.resolve()):
        raise SpecError(f"filename {filename!r} resolves outside the workspace "
                        f"{out_dir}")
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise SpecError(f"cannot create output dir {out_dir}: {exc}")
    _metadata = {"Date": None} if fmt == "svg" else {}
    fig.savefig(str(out_path), format=fmt, metadata=_metadata)
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
    out_dir_path = Path(out_dir or os.getcwd())
    series, n_points = _resolve_series(spec, kind, out_dir_path)

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
    """The escaped card BODY — the div+figure+caption that goes inside the
    sandboxed iframe that ``to_card`` builds. Caption is html.escaped so
    hostile spec/data-derived text cannot inject markup; the figure travels
    as a ``data:`` URI served as an ``<img>`` (an embedded SVG never executes
    scripts in an image). This body alone is never inlined into the dashboard
    raw — it must be escaped again inside ``srcdoc="..."`` (see to_card).
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


def to_card(result: PlotResult) -> str:
    """The self-contained card the dashboard inlines for make_plot, behind
    the shared ``DELFIN_CARD:`` marker (one hook with V1 show_molecule).

    Security (operator contract): the figure + caption are rendered inside an
    ``<iframe sandbox="allow-scripts" srcdoc="...">`` without
    ``allow-same-origin``, so the card runs scripts but is isolated and cannot
    read the dashboard. The whole body is escaped a second time inside the
    ``srcdoc`` attribute (quotes -> &quot;, < > -> &lt; &gt;), so hostile
    spec/data-derived text cannot terminate the attribute or emit a raw
    </script>/</iframe> that climbs out of the sandbox into the parent.
    """
    import html as _html

    srcdoc = _html.escape(to_html(result), quote=True)
    return (
        f'<iframe sandbox="allow-scripts" '
        f'srcdoc="{srcdoc}" '
        f'frameborder="0" style="width:100%;min-height:380px;'
        f'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
        f'</iframe>'
    )


# --------------------------------------------------------------------------
# the make_plot tool result (shared DELFIN_CARD: contract with package V1)


#: Marker every chat-card tool result starts with. When the shared V1 emit
#: point (delfin.dashboard.chat_viewer.tool_result) is present we delegate to
#: it so the marker + keys live in ONE place across both packages; this is the
#: conforming local fallback for branches where chat_viewer is not merged yet.
CARD_MARKER = "DELFIN_CARD:"


def _emit_card(html_card: str, text: str) -> str:
    """Prefer the shared V1+V2 emit point, fall back to a conforming local one.

    chat_viewer.tool_result is the single emit point (V1 builder, one hook for
    both packages). It is only reachable once V1's dashboard module is merged
    to main — which has not happened on this branch — so we import it lazily
    and keep a byte-identical local emit for the interim, so make_plot works
    on the current base without a hard dependency on an unmerged module.
    """
    try:
        from delfin.dashboard import chat_viewer as _cv
        return _cv.tool_result(html_card, text)
    except (ImportError, AttributeError):
        payload = json.dumps({
            "escape": "html",
            "html": html_card,
            "text": text,
        })
        return CARD_MARKER + payload


def make_plot(spec: dict, *, out_dir: str = "") -> str:
    """Render ``spec`` and return the shared ``DELFIN_CARD:`` tool result.

    The single tool entry point for package V2: draws the figure (see
    :func:`plot`), wraps it as the sandboxed iframe card (see :func:`to_card`)
    and returns the one string the dashboard's role="tool" hook inlines and a
    terminal shows as its plain ``text`` line — ``DELFIN_CARD:`` followed by
    ``{"escape":"html","html":<iframe card>,"text":<plain path + caption>}``.

    The model-facing ``text`` is deliberately plain (no HTML): the dashboard
    renders the escaped card from ``html`` (an isolated sandboxed iframe);
    the model receives only ``text`` — the workspace path and one-line
    caption — per the operator's data-is-never-instruction rule. A malformed
    spec raises :class:`SpecError` exactly as :func:`plot` does, before
    anything is written.
    """
    result = plot(spec, out_dir=out_dir)
    html_card = to_card(result)
    text = f"Wrote {result.path} — {result.caption}"
    return _emit_card(html_card, text)
