"""Declarative figures to the agent workspace (the ``make_plot`` tool).

Input: a JSON spec (below). Output: a PNG or SVG written inside the
workspace, and the ``make_plot`` tool result -- the ``DELFIN_CARD:`` marker
followed by ``{"card": "plot", "path", "title", "text"}``. The result
carries no figure bytes and no markup: the dashboard reads the file from
``path`` and shows it with the chat's existing artifact preview
(``tab_agent._make_plot_card_html``); the model and the terminal receive
``text`` and the caption lines that follow the object
(``chat_viewer.card_model_text``).

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
import re
import stat
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
#: Source lines listed in an energy profile's caption; further states are
#: counted, not listed. Bounds what the model re-reads every round.
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


def read_text(path: Path, max_bytes: int) -> str:
    """The text of the regular file ``path``, of at most ``max_bytes``.

    Decoded as UTF-8 with undecodable bytes replaced (outputs echo their
    input in whatever encoding it was written). Existence and size are
    looked at here and nowhere before: a caller with a read gate calls
    this only once the gate admitted ``path``, so a refused path is
    neither found, missed nor measured. A file larger than ``max_bytes``
    is refused before it is opened.
    """
    try:
        st = path.stat()
    except FileNotFoundError:
        raise SpecError(f"not found: {path}") from None
    except OSError as exc:
        raise SpecError(f"cannot read {path}: {exc.strerror}") from None
    if not stat.S_ISREG(st.st_mode):
        raise SpecError(f"not a regular file: {path}")
    if st.st_size > max_bytes:
        raise SpecError(f"{path} is {st.st_size} bytes; at most {max_bytes}")
    try:
        with path.open("rb") as fh:
            data = fh.read(max_bytes + 1)
    except OSError as exc:
        raise SpecError(f"cannot read {path}: {exc.strerror}") from None
    if len(data) > max_bytes:               # grew since the stat
        raise SpecError(f"{path} is over {max_bytes} bytes")
    return data.decode("utf-8", errors="replace")


def _read_data_file(path, root: Path, reader=None) -> tuple[str, str]:
    """``(text, suffix)`` of a CSV or JSON file inside ``root``.

    A path that resolves outside ``root`` (absolute, ``..``, symlink) is
    refused, as is a file over ``MAX_DATA_FILE_BYTES``. ``reader`` -- the
    dispatch layer's read gate (secret-deny list, read roots) -- is called
    as ``reader(path, max_bytes)`` with the resolved path before the file
    is looked at, and returns the text (:func:`read_text` once admitted);
    without it :func:`read_text` reads the file directly.
    """
    if not isinstance(path, str) or not path:
        raise SpecError("data.file must be a path string")
    root_res = root.resolve()
    p = (root / path if not Path(path).is_absolute() else Path(path)).resolve()
    if not _is_within(p, root_res):
        raise SpecError(f"data file {path!r} resolves outside the workspace "
                        f"{root}")
    suffix = p.suffix.lower()
    if suffix not in (".csv", ".json"):
        raise SpecError(f"unsupported data file type: {suffix or path}")
    return (reader or read_text)(p, MAX_DATA_FILE_BYTES), suffix


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


def _as_columns(raw: list) -> list:
    """Series in the column form (``x`` + ``y`` or ``files``) the profile
    is built from.

    Two other shapes come in. ``states: [{name, value | files}]`` is the
    one models write when left to themselves -- an object per state -- and
    a list of files inside a state is then plainly its fragments. And a run
    of one-state series (each state given as its own series, as live runs
    did twice) is one pathway: a pathway of one state joins nothing, so the
    reading is unambiguous. A scalar ``x``/``y`` is a one-element list.
    """
    out = []
    for idx, s in enumerate(raw):
        if not isinstance(s, dict):
            raise SpecError(f"series[{idx}] must be an object")
        s = dict(s)
        for k in ("x", "y"):
            if k in s and not isinstance(s[k], list) and s[k] is not None:
                s[k] = [s[k]]
        if isinstance(s.get("files"), str):
            s["files"] = [s["files"]]
        states = s.pop("states", None)
        if states is not None:
            s.update(_states_to_columns(states, idx))
        out.append(s)

    def one(s):
        n = len(s.get("files") or s.get("y") or [])
        return n == 1 and "states" not in s
    if len(out) > 1 and all(one(s) for s in out):
        names = [str((s.get("x") or [s.get("label") or f"{i + 1}"])[0])
                 for i, s in enumerate(out)]
        merged = {"x": names}
        if all(s.get("files") for s in out):
            merged["files"] = [s["files"][0] for s in out]
        elif not any(s.get("files") for s in out):
            merged["y"] = [s["y"][0] for s in out]
        else:
            raise SpecError("give every state either a value or files, not "
                            "a mix")
        merged["from_user"] = any(s.get("from_user") is True for s in out)
        out = [merged]
    return out


def _states_to_columns(states, idx) -> dict:
    """``[{name, value | files}]`` -> ``{x, y}`` or ``{x, files}``."""
    if not isinstance(states, list) or not states:
        raise SpecError(f"series[{idx}].states must be a non-empty list of "
                        f"{{name, value or files}}")
    names, values, files = [], [], []
    for k, st in enumerate(states):
        if not isinstance(st, dict):
            raise SpecError(f"series[{idx}].states[{k}] must be an object "
                            f"{{name, value or files}}")
        names.append(str(st.get("name", st.get("label", f"{k + 1}"))))
        f = st.get("files", st.get("file"))
        values.append(st.get("value", st.get("y")))
        files.append(f)
    if any(f is not None for f in files):
        if any(v is not None for v, f in zip(values, files) if f is None):
            raise SpecError(f"series[{idx}]: give every state either a value "
                            f"or files, not a mix")
        return {"x": names, "files": files}
    return {"x": names, "y": values}


def _refuse_hartree_as_given(parsed, spec: dict) -> None:
    """Typed-in values that are Hartree, not kcal/mol, are refused.

    A live run read xtb energies, passed their difference in Eh as ``y``
    -- -0.0334 for an adduct bound by 20.95 kcal/mol -- and ignored the
    note that typed-in values are unchecked; the figure was 627 times too
    flat. Two shapes give it away: every value under 1 (a reaction profile
    in kcal/mol that flat is not one), or every value large and within a
    few units of the others (absolute energies). A user whose numbers
    really are that small says so with ``units``.
    """
    if spec.get("units") is not None:
        return
    given = [v for *_, ys, _, files in parsed if not files
             for v in ys if v is not None]
    if len(given) < 2:
        return
    small = max(abs(v) for v in given) < 1.0 and any(given)
    absolute = (min(abs(v) for v in given) > 20.0
                and max(given) - min(given) < 2.0)
    if small or absolute:
        what = ("all under 1" if small else
                "all large and within 2 of each other -- absolute energies")
        raise SpecError(
            f"the values are {what}: they look like Hartree, not kcal/mol. "
            f"Pass the outputs as files: [...] (the tool converts them), or "
            f"convert relative to the first state x 627.509; if they really "
            f"are kcal/mol, set units: \"kcal/mol\"")


def _read_output(path, root: Path, reader=None,
                 texts: dict | None = None) -> tuple[str, Path]:
    """``(text, resolved path)`` of an output for ``energy_profile``.

    Unlike a data file it may lie outside the workspace -- outputs live in
    the calculation folders -- so the dispatch layer's read gate
    (``reader``: secret-deny list, read roots) decides, exactly as it does
    for ``read_file``. The gate runs first: ``reader(path, max_bytes)``
    looks at the file (:func:`read_text`) only once it admitted the path,
    so a refusal is the same whether the path exists or not. Bounded by
    ``MAX_OUTPUT_BYTES``. ``texts`` keeps what one call already read, so
    each output passes the gate once.
    """
    if not isinstance(path, str) or not path.strip():
        raise SpecError("each files entry must be a path string or null")
    p = Path(path).expanduser()
    p = (p if p.is_absolute() else root / p).resolve()
    if texts is not None and p in texts:
        return texts[p], p
    text = (reader or read_text)(p, MAX_OUTPUT_BYTES)
    if texts is not None:
        texts[p] = text
    return text, p


def _resolve_pathways(spec: dict, root: Path, reader=None):
    """``([Pathway, ...], n_points, notes, quantity)`` for an
    ``energy_profile`` spec; ``quantity`` is what the source says it holds
    (``ΔG`` from free energies or a ``deltaG`` CSV header), else None.

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
            return (pathways, sum(len(p.values) for p in pathways), [],
                    ep.csv_quantity(text))
        try:
            raw = json.loads(text)
        except Exception as exc:
            raise SpecError(f"cannot read {data['file']}: {exc}") from exc
        spec = dict(spec, series=(raw.get("series")
                                  if isinstance(raw, dict) else raw))
    raw = spec.get("series")
    if not raw and spec.get("states") is not None:
        raw = [{"states": spec["states"], "label": spec.get("label", "")}]
    if not raw:
        if "y" not in spec and "files" not in spec:
            raise SpecError(
                "kind 'energy_profile' needs 'states' (one pathway: "
                "[{name, value or files}]), 'series' (several: [{label, "
                "states}]), or data: {file: <csv>}")
        raw = [{k: spec[k] for k in ("x", "y", "files", "start")
                if k in spec}]
    if not isinstance(raw, list):
        raise SpecError("series must be a list of objects")
    if len(raw) > MAX_SERIES:
        raise SpecError(f"{len(raw)} pathways; at most {MAX_SERIES}")
    raw = _as_columns(raw)

    energy = _require_str(spec, "energy", "auto").strip().lower()
    if energy not in ("gibbs", "electronic", "composite", "auto"):
        raise SpecError(f"energy must be gibbs, electronic, composite or "
                        f"auto, got {energy!r}")
    corrections = (spec.get("corrections") or "").strip()
    if energy == "composite" and not corrections:
        raise SpecError("energy: \"composite\" needs a \"corrections\" "
                        "key naming the low level of theory whose thermal "
                        "correction (G - E) is added to the high-level "
                        "electronic energy, e.g. corrections: "
                        "\"r2SCAN-3c\"")
    energy_note = ""
    # Each output once through the read gate: deciding which energy the
    # outputs carry reads them before they are drawn.
    texts: dict = {}
    if energy == "auto":
        energy, energy_note = _energy_the_outputs_carry(raw, spec, root,
                                                        reader, texts)

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
                    outputs[key] = _state_record(key, root, reader,
                                                 energy, texts)
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
            outputs[ref_file] = _state_record(ref_file, root, reader, energy,
                                            texts)
        ref = outputs[ref_file]
        ref_eh = ref["eh"]
        programs = sorted({part for rec in outputs.values()
                           for part in rec["program"].split("/")})
        if len(programs) > 1:
            # Energies of two methods are not on one surface; drawn side by
            # side the gap between them would read as a reaction step.
            raise SpecError(
                f"the outputs come from different programs "
                f"({', '.join(programs)}): their energies are not on one "
                f"surface -- every state must be computed the same way")
        notes.append(f"{energy} energies relative to {' + '.join(ref_file)} "
                     f"({ref_eh:.6f} Eh), x {factor:g} {units}/Eh")
        if energy_note:
            notes.append(energy_note)
        if energy == "composite":
            notes.append("composite G = E(high-level single point) + "
                         f"G_corr({corrections}); G_corr is G - E of the "
                         "freq output")
        _check_balance(parsed, outputs, ref_file, notes)
        if energy == "gibbs":
            _standard_state_note(parsed, notes)

    typed = any(not files for *_, files in parsed)
    from_user = spec.get("from_user") is True or any(
        isinstance(r, dict) and r.get("from_user") is True for r in raw)
    if typed and not from_user:
        # Live runs read the outputs, typed the numbers in themselves and
        # then reported "no imaginary mode" for an adduct that had one -- a
        # note in the result was ignored three times. Typed-in values are
        # the user's numbers; anything read from a calculation goes through
        # files, where the tool converts and checks it.
        raise SpecError(
            "y is for numbers the USER gave you. If you read them from "
            "calculation outputs, pass the paths as files: [...] instead -- "
            "the tool converts them and checks mass balance and imaginary "
            "modes. If the user gave them, add \"from_user\": true.")
    _refuse_hartree_as_given(parsed, spec)

    pathways = []
    for idx, s, xs, ys, start, files in parsed:
        if files:
            ys = [None if v is None else (v - ref_eh) * factor for v in ys]
            for state, f, v in zip(xs, files, ys):
                if f is None:
                    continue
                rec = outputs[f]
                notes.append(f"{state}: {' + '.join(f)} ({rec['program']}) "
                             f"{rec['eh']:.6f} Eh -> {v:.2f}"
                             + _frequency_note(state, rec["imag"])
                             + (f" [{rec['opi']}]" if rec.get("opi") else ""))
        pw = ep.Pathway(label=_require_str(s, "label", ""),
                        states=tuple(xs), values=tuple(ys), start=start)
        try:
            ep.check(pw, f"series[{idx}]")
        except ValueError as exc:
            raise SpecError(str(exc)) from None
        pathways.append(pw)
    quantity = ({"gibbs": "ΔG", "composite": "ΔG",
                 "electronic": "ΔE"}[energy] if outputs else None)
    if quantity:
        _check_label_quantity(_require_str(spec, "ylabel", "").strip(),
                              quantity, energy)
    return pathways, n, notes, quantity


#: How a y label names an energy: the symbol (``ΔG``, ``\Delta G``,
#: ``deltaG``, ``dG``; the letter case-sensitive) or the word.
_LABEL_QUANTITY = {
    q: re.compile(
        r"(?:Δ|\\Delta|(?<![A-Za-z])[dD]elta|(?<![A-Za-z\\])d)"
        r"[\s{]*" + q[1] + r"(?![a-z])"
        + r"|(?i:" + words + r")")
    for q, words in (("ΔG", r"gibbs|free\s+energ"),
                     ("ΔH", r"enthalp"),
                     ("ΔE", r"electronic|single[\s-]point|\bscf\b"))
}


def _check_label_quantity(label: str, quantity: str, energy: str) -> None:
    """Refuse a y label that names another energy than the outputs give.

    ``quantity`` is fixed by what was read (``energy``: Gibbs or
    electronic). A label naming it, or naming no energy ("Relative
    energy"), is kept; one naming another (``ΔE`` over Gibbs energies) is
    refused, since the axis would misstate the data without a word.
    """
    named = {q for q, rx in _LABEL_QUANTITY.items() if rx.search(label)}
    other = sorted(named - {quantity})
    if not label or not other:
        return
    fix = (", or pass energy: \"electronic\" to draw ΔE" if "ΔE" in other
           else ", or pass energy: \"gibbs\" to draw ΔG" if "ΔG" in other
           else "; ΔH is not read from outputs")
    raise SpecError(
        f"ylabel {label!r} contradicts the outputs: it names "
        f"{' and '.join(other)}, they give {quantity} (energy: "
        f"\"{energy}\"). Leave ylabel out (the axis then reads "
        f"{quantity}){fix}")


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


def _energy_the_outputs_carry(raw, spec, root: Path, reader,
                              texts: dict) -> tuple[str, str]:
    """``(energy, note)`` for a profile drawn without an ``energy`` key.

    Every output of the call counts -- all pathways, every fragment of a
    state, and ``reference``: ``gibbs`` when each carries a free energy,
    ``electronic`` with a note saying so when none does. A mix is refused
    naming both groups: ΔG and ΔE are not one surface, and the gap between
    them would read as a reaction step. Without outputs (typed values) the
    choice is moot and ``gibbs`` is returned.
    """
    from delfin.agent import energy_profile as ep

    keys = []
    for idx, s in enumerate(raw):
        if isinstance(s, dict) and isinstance(s.get("files"), list):
            keys += [_state_key(f, idx) for f in s["files"] if f is not None]
    if keys and spec.get("reference") is not None:
        keys.append(_state_key(spec["reference"], "reference"))
    with_g: list[str] = []
    without: list[str] = []
    for path in dict.fromkeys(p for key in keys for p in key):
        text, resolved = _read_output(path, root, reader, texts)
        (with_g if ep.has_free_energy(text, str(resolved))
         else without).append(path)
    if not without:
        return "gibbs", ""
    if not with_g:
        return "electronic", (
            "no output carries a free energy (no frequency calculation): "
            "ΔE from electronic energies; frequencies on every state give ΔG")

    def _names(paths: list[str]) -> str:
        shown = ", ".join(paths[:6])
        return shown + (f" and {len(paths) - 6} more" if len(paths) > 6
                        else "")

    raise SpecError(
        f"free energies in {_names(with_g)}, none in {_names(without)}: "
        f"ΔG and ΔE are not one surface. Run frequencies on the second "
        f"group, or pass energy: \"electronic\" to draw ΔE for every state")


def _state_record(key: tuple[str, ...], root: Path, reader,
                  energy: str, texts: dict | None = None) -> dict:
    """``{"eh", "imag", "atoms", "charge"}`` of one state, its fragments
    summed; a refusal naming the output that lacks its energy."""
    from collections import Counter
    from delfin.agent import energy_profile as ep

    if energy == "composite":
        return _composite_record(key, root, reader, texts)

    eh, imag, atoms, charge, programs = 0.0, 0, Counter(), 0, []
    witnessed = []
    for path in key:
        text, resolved = _read_output(path, root, reader, texts)
        try:
            rec = ep.read_output(text, str(resolved), energy)
        except ValueError as exc:
            raise SpecError(f"{exc}; nothing is drawn for a state without "
                            f"its energy") from None
        if rec["program"] == "ORCA":
            witnessed.append(_opi_says(resolved, rec["eh"], energy, reader))
        eh += rec["eh"]
        programs.append(rec["program"])
        n = rec["imag"]
        imag = None if imag is None or n is None else imag + n
        comp, q = rec["atoms"], rec["charge"]
        if atoms is None or comp is None or q is None or charge is None:
            atoms = charge = None
        else:
            atoms += comp
            charge += q
    # One verdict for the state: any disagreement wins, and "confirmed"
    # only when every fragment was confirmed.
    opi = next((w for w in witnessed if "WARNING" in w), "")
    if not opi and witnessed and all(witnessed):
        opi = "confirmed by OPI"
    return {"eh": eh, "imag": imag, "atoms": atoms, "charge": charge,
            "program": "/".join(dict.fromkeys(programs)), "opi": opi}


def _composite_record(key: tuple[str, ...], root: Path, reader,
                      texts: dict | None = None) -> dict:
    """``{"eh", "imag", "atoms", "charge", "program"}`` of a composite state.

    ``key`` must be the two-path pair ``(elec.out, freq.out)`` a composite
    ``files`` entry produces -- the high-level single point and the low-level
    ``Opt Freq`` output. Its Gibbs energy is E(elec) + (G - E)(freq), read by
    :func:`delfin.agent.energy_profile.composite_output`. A pair of any other
    length is refused, and a missing energy names the output.
    """
    from delfin.agent import energy_profile as ep

    if len(key) != 2:
        raise SpecError(
            "energy: \"composite\" pairs the high-level single point with "
            "the low-level freq output as a two-element files entry "
            f"[elec.out, freq.out] (two outputs); got {len(key)} path(s)")
    hi, lo = key[0], key[1]
    hi_text, hi_res = _read_output(hi, root, reader, texts)
    lo_text, lo_res = _read_output(lo, root, reader, texts)
    try:
        rec = ep.composite_output(hi_text, str(hi_res),
                                  lo_text, str(lo_res))
    except ValueError as exc:
        raise SpecError(f"{exc}; nothing is drawn for a state without "
                        f"its composite energy") from None
    return {"eh": rec["eh"], "imag": rec["imag"], "atoms": rec["atoms"],
            "charge": rec["charge"], "program": rec["program"], "opi": ""}


def _opi_says(output: Path, eh: float, energy: str, reader) -> str:
    """ORCA's property.json beside ``output``, through OPI, on ``eh``.

    Used automatically when the [orca] extra is installed and the JSON is
    there (:mod:`delfin.opi_witness`); otherwise ``""`` and nothing
    changes. The JSON passes the same read gate as the output.
    """
    from delfin import opi_witness

    found = opi_witness.json_beside(output)
    if found is None or not opi_witness.available():
        return ""
    if reader is not None:
        try:
            reader(found)
        except SpecError:
            return ""
    return opi_witness.compare(eh, opi_witness.structured(found), energy)


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


def _standard_state_note(parsed, notes) -> None:
    """Say so when the number of separate molecules changes along a profile.

    Free energies from outputs refer to 1 atm; in solution the reference is
    1 M, RT ln 24.46 = 1.89 kcal/mol per molecule at 298.15 K, which does
    not cancel when molecules associate or dissociate. Not applied: the
    temperature, the phase and whether the outputs already include it are
    the caller's to state.
    """
    counts = [len(key) for *_, files in parsed for key in files or []
              if key is not None]
    if counts and min(counts) != max(counts):
        notes.append(
            f"number of molecules changes ({counts[0]} -> "
            f"{next(c for c in counts if c != counts[0])}): the 1 atm -> 1 M "
            f"standard-state correction (RT ln 24.46 = 1.89 kcal/mol per "
            f"molecule at 298.15 K) is NOT applied")


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


#: Lines with at most this many points get a marker on each one.
_MARKERS_UP_TO = 80


def _house_style(ax) -> None:
    """The look every line/scatter/bar/histogram/box plot shares: no
    top/right frame, a faint grid under the data, readable tick labels.
    Set on the axes, not through rcParams, so the caller's process keeps
    its own matplotlib defaults (see the module docstring)."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis="y", alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=10)


def resolve_out_path(spec: dict, out_dir: Path) -> Path:
    """The file :func:`plot` writes for ``spec`` under ``out_dir``.

    Shared by the renderer and the dispatch executor, so the path the write
    gate checks is the path written. Refuses a name outside ``out_dir``.
    """
    kind = _require_kind(spec)
    # An energy profile goes on a slide or into a paper draft: PNG at
    # 300 dpi as the group's scripts wrote it. Two live runs asked for a PNG
    # without saying format and then converted the SVG by hand.
    filename = _require_str(spec, "filename", "").strip()
    # A filename that ends in a format IS the format when none is given:
    # "butane_pes.svg" with the energy profile's PNG default was written as
    # "butane_pes.svg.png".
    named = Path(filename).suffix[1:].lower() if filename else ""
    default_fmt = (named if named in _FORMATS
                   else "png" if kind == "energy_profile" else "svg")
    fmt = _require_str(spec, "format", default_fmt).lower()
    if fmt not in _FORMATS:
        raise SpecError(f"format must be png or svg, got {fmt!r}")
    filename = filename or f"plot_{kind}.{fmt}"
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
    at :func:`resolve_out_path`. ``reader(path, max_bytes) -> str`` is the
    caller's read gate for data files and outputs: it refuses with a
    :class:`SpecError` or returns :func:`read_text`.
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

    fig = Figure(figsize=(7.5, 4.2), layout="constrained")
    ax = fig.subplots()
    _house_style(ax)
    if kind in ("line", "scatter", "bar"):
        draw = {"line": ax.plot, "scatter": ax.scatter, "bar": ax.bar}[kind]
        for label, xs, ys in series:
            pos, ticks = _positions(xs) if xs else (list(range(len(ys))), None)
            extra = {}
            if kind == "line":
                # Each computed point visible while there are few enough to
                # read: a smooth line over 25 calculations hides that they
                # are 25.
                extra = {"linewidth": 1.8}
                if len(ys) <= _MARKERS_UP_TO:
                    extra.update(marker="o", markersize=4)
            draw(pos, ys, label=label or None, **extra)
            if ticks is not None:
                ax.set_xticks(pos, ticks)
        if any(lbl for lbl, _, _ in series):
            ax.legend(frameon=False)
    elif kind == "histogram":
        ax.hist(series[0][2], bins=bins if bins is not None else "auto")
    else:  # box
        ax.boxplot([series[0][2]], tick_labels=[xlabel])
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=11.5)
    if ylabel:
        ax.set_ylabel(f"{ylabel} ({units})" if units else ylabel,
                      fontsize=11.5)
    if title:
        ax.set_title(title, fontsize=13, fontweight="semibold", loc="left")

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
    units = (_require_str(spec, "units", "kcal/mol").strip()
             if spec.get("units") is not None else "kcal/mol")
    pathways, n_points, notes, quantity = _resolve_pathways(
        spec, out_dir_path, reader)
    ylabel = _require_str(spec, "ylabel", "").strip() or quantity
    if not ylabel:
        # A live run drew electronic energies under the default ΔG axis.
        # Typed-in numbers do not say what they are; the caller must.
        raise SpecError(
            "say which energy the values are: ylabel \"ΔG\", \"ΔE\" or "
            "\"ΔH\" (from files the tool knows it; typed-in values do not "
            "say)")
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
        # What contradicts the picture first, then what was not checked or
        # applied, so a cut caption keeps both.
        notes = sorted(notes, key=lambda n: (0 if "WARNING" in n else
                                             1 if "NOT" in n else 2))
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
    "title": <file name>, "text": "Wrote <path> -- <first caption line>"}``,
    then, when the caption has more lines (an energy profile's source
    notes), a newline and those lines as plain text.

    The dashboard parses the JSON object from the first 2000 characters of
    a result (``chat_viewer.card_payload``). The object holds the path
    three times and one caption line (``MAX_TEXT`` bounds the title in it),
    so it fits for any realistic path; the notes grow with the number of
    states -- 11 states of ~90-character paths passed 2000 -- and travel
    after it, where the preview need not parse them. The model receives
    ``text`` and the notes (``chat_viewer.card_model_text``).

    Rejected: a larger preview slice (any fixed limit is passed by enough
    states) and fewer or shorter notes (the paths and per-state checks are
    what the model reports from).
    """
    head, _, notes = result.caption.partition("\n")
    card = CARD_MARKER + json.dumps({
        "card": "plot",
        "path": str(result.path),
        "title": result.path.name,
        "text": f"Wrote {result.path} -- {head}",
    })
    return card + "\n" + notes if notes else card


def make_plot(spec: dict, *, out_dir: str = "", reader=None) -> str:
    """Render ``spec`` (:func:`plot`) and return :func:`card_result`."""
    return card_result(plot(spec, out_dir=out_dir, reader=reader))
