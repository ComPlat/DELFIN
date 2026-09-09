"""Jupyter notebook (.ipynb) read / edit primitives.

Notebooks are JSON; ``read_file`` would dump the raw structure (verbose
and noisy) and ``edit_file`` would corrupt cell delimiters when string-
matching. This module provides cell-aware operations:

* :func:`read_cells` — return a list of ``{idx, cell_type, source,
  output_summary}`` records, with outputs collapsed to a one-line
  summary so the agent sees structure without drowning in
  base64-encoded plot images.
* :func:`apply_edit` — atomic mutation: ``replace``, ``insert_before``,
  ``insert_after``, or ``delete`` a cell at a given index.

Both functions take an already-resolved absolute path. The api_client
caller is responsible for sandbox checks (``_resolve_in_workspace``)
and Self-Modification-Guard before calling.

Notebook ``source`` can be either a list of strings or a single
string per the nbformat spec; we always normalise to a string for the
agent and always emit a list of strings on write (with trailing
newlines on all but the last line — that's the format Jupyter writes
itself, so re-saving doesn't churn the diff).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional


_VALID_MODES = {"replace", "insert_before", "insert_after", "delete"}
_VALID_CELL_TYPES = {"code", "markdown", "raw"}


@dataclass
class CellInfo:
    idx: int
    cell_type: str
    source: str
    output_summary: str  # what kind of output, and how much
    output_text: str     # the output itself, text only, capped


# What a notebook's outputs are MADE of, and what is worth carrying.
#
# A figure's PNG is the thing this module was written to keep out of the
# prompt, and it stays out. Everything else in an output is text a
# scientist opened the notebook to read: the printed numbers, the repr of
# the array, and above all the traceback of the cell that failed.
_TEXT_MIMES = ("text/plain", "text/markdown", "application/json")
_DEFAULT_OUTPUT_CHARS = 2000

# Jupyter writes tracebacks with terminal colour codes in them.
_ANSI_RE = None  # compiled lazily; see _strip_ansi


def _source_to_string(source: Any) -> str:
    """nbformat allows source as list[str] or str; normalise to str."""
    if isinstance(source, list):
        return "".join(source)
    if isinstance(source, str):
        return source
    return ""


def _string_to_source(text: str) -> list[str]:
    """Emit nbformat-style list-of-lines (each but last with trailing \\n)."""
    if text == "":
        return []
    lines = text.split("\n")
    out = [ln + "\n" for ln in lines[:-1]]
    if lines[-1] != "":
        out.append(lines[-1])
    return out


def _strip_ansi(text: str) -> str:
    global _ANSI_RE
    if _ANSI_RE is None:
        import re
        _ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
    return _ANSI_RE.sub("", text)


def _joined(value: Any) -> str:
    """nbformat writes text as list[str] or str, the same as source."""
    if isinstance(value, list):
        return "".join(v for v in value if isinstance(v, str))
    if isinstance(value, str):
        return value
    return ""


def _output_text(outputs: list[Any], max_chars: int) -> str:
    """The text a cell produced, without the bytes it produced.

    ``_summarise_outputs`` says what KIND of output a cell has, and for a
    plot that is the whole useful answer. For everything else it is not:
    the printed number, the repr of the array and the traceback of the
    cell that failed are what a scientist opens the notebook to read, and
    a summary reading ``1 output(s): error(ValueError)`` names the
    exception class and drops the message, the line and the stack.

    That made the specialised tool return strictly less of what matters
    than the generic ``read_file`` it tells the agent to use it instead
    of -- read_file would at least have dumped the traceback as JSON.

    Binary MIME data is still dropped, which was the point of the
    summary: a base64 PNG is thousands of tokens of nothing readable.
    """
    if not outputs:
        return ""
    parts: list[str] = []
    for o in outputs:
        if not isinstance(o, dict):
            continue
        kind = o.get("output_type", "?")
        if kind == "stream":
            parts.append(_joined(o.get("text", "")))
        elif kind == "error":
            # ename alone is not a diagnosis. The traceback's last frames
            # are where the fault is; its first are import machinery.
            head = f"{o.get('ename', 'Error')}: {o.get('evalue', '')}".strip()
            tb = _strip_ansi("\n".join(
                ln for ln in (o.get("traceback") or [])
                if isinstance(ln, str)))
            parts.append(f"{head}\n{tb}".strip())
        elif kind in ("execute_result", "display_data"):
            data = o.get("data") or {}
            if not isinstance(data, dict):
                continue
            for mime in _TEXT_MIMES:
                if mime in data:
                    parts.append(_joined(data[mime]))
                    break
    text = "\n".join(p for p in parts if p).strip()
    if len(text) > max_chars:
        # The TAIL, not the head. Source truncates in the middle because
        # a function's name is at the top and its return at the bottom;
        # output is the other way round. A traceback's first frames are
        # import machinery and its last is the line that failed, and a
        # loop that printed for ten minutes ends with the result.
        keep = max(0, max_chars - 60)
        if keep:
            dropped = len(text) - keep
            text = (f"... ({dropped} chars of earlier output omitted)\n"
                    + text[-keep:])
        else:
            text = f"... ({len(text)} chars of output omitted)"
    return text


def _summarise_outputs(outputs: list[Any]) -> str:
    """Collapse a cell's outputs to one short line.

    Says what KIND of output a cell has and how much of it there is. The
    output itself is carried separately by :func:`_output_text`.
    """
    if not outputs:
        return ""
    types: list[str] = []
    text_chars = 0
    for o in outputs:
        if not isinstance(o, dict):
            continue
        t = o.get("output_type", "?")
        if t == "stream":
            text_chars += sum(len(s) for s in o.get("text", []) if isinstance(s, str))
            types.append(f"stream/{o.get('name', '?')}")
        elif t in ("execute_result", "display_data"):
            # `"data": null` is out of spec and does occur -- a
            # hand-edited or truncated notebook. `.keys()` on it raised
            # out of read_cells, whose caller turns any exception into
            # "notebook read failed", so one malformed output made the
            # whole notebook unreadable.
            data = o.get("data") or {}
            mimes = list(data.keys()) if isinstance(data, dict) else []
            types.append(f"{t}({','.join(mimes[:2])})")
        elif t == "error":
            ename = o.get("ename", "Error")
            types.append(f"error({ename})")
        else:
            types.append(t)
    if text_chars:
        return f"{len(outputs)} output(s): {', '.join(types[:4])} (~{text_chars} chars)"
    return f"{len(outputs)} output(s): {', '.join(types[:4])}"


def read_cells(
    path: Path,
    max_source_chars: int = 4000,
    max_output_chars: int = _DEFAULT_OUTPUT_CHARS,
) -> list[CellInfo]:
    """Parse a notebook and return cell records.

    ``max_source_chars`` truncates each cell's source so a single
    huge cell can't blow the agent's tool-result window. The result
    keeps the head + tail with a marker in the middle.

    ``max_output_chars`` does the same for each cell's OUTPUT text --
    a loop that printed for ten minutes is exactly as capable of
    flooding the window as a long cell, and unlike source, output is
    usually most informative at its end.
    """
    with path.open("r", encoding="utf-8") as fh:
        nb = json.load(fh)
    cells = nb.get("cells", []) or []
    records: list[CellInfo] = []
    for i, cell in enumerate(cells):
        if not isinstance(cell, dict):
            continue
        source = _source_to_string(cell.get("source", ""))
        if len(source) > max_source_chars:
            head = source[: max_source_chars // 2]
            tail = source[-max_source_chars // 2:]
            omitted = len(source) - len(head) - len(tail)
            source = (
                f"{head}\n... ({omitted} chars from middle of cell omitted)\n{tail}"
            )
        outputs = cell.get("outputs", []) or []
        records.append(CellInfo(
            idx=i,
            cell_type=cell.get("cell_type", "?"),
            source=source,
            output_summary=_summarise_outputs(outputs),
            output_text=_output_text(outputs, max_output_chars),
        ))
    return records


def apply_edit(
    path: Path,
    cell_idx: int,
    mode: str,
    source: Optional[str] = None,
    cell_type: str = "code",
) -> tuple[int, int]:
    """Apply a single cell-level edit, atomically.

    Returns ``(cells_before, cells_after)`` so the caller can render a
    succinct status line. Raises ValueError on invalid input.

    The whole notebook is read, mutated, then written back via
    ``Path.write_text``. nbformat's ``execution_count`` and ``id``
    fields are preserved on replace; new cells get a random id and
    null execution_count (Jupyter regenerates these on save).
    """
    if mode not in _VALID_MODES:
        raise ValueError(
            f"mode must be one of {sorted(_VALID_MODES)}, got {mode!r}"
        )
    if cell_type not in _VALID_CELL_TYPES:
        raise ValueError(
            f"cell_type must be one of {sorted(_VALID_CELL_TYPES)}, got {cell_type!r}"
        )
    if mode != "delete" and source is None:
        raise ValueError(f"source is required for mode={mode!r}")

    with path.open("r", encoding="utf-8") as fh:
        nb = json.load(fh)
    cells = nb.get("cells", []) or []
    n_before = len(cells)

    if mode == "delete":
        if not (0 <= cell_idx < n_before):
            raise ValueError(
                f"cell_idx {cell_idx} out of range (notebook has {n_before} cells)"
            )
        cells.pop(cell_idx)
    elif mode == "replace":
        if not (0 <= cell_idx < n_before):
            raise ValueError(
                f"cell_idx {cell_idx} out of range (notebook has {n_before} cells)"
            )
        old = cells[cell_idx]
        new_cell = {
            "cell_type": cell_type,
            "metadata": old.get("metadata", {}),
            "source": _string_to_source(source or ""),
        }
        if cell_type == "code":
            new_cell["outputs"] = []
            new_cell["execution_count"] = None
        # Preserve id if present (nbformat 4.5+).
        if "id" in old:
            new_cell["id"] = old["id"]
        cells[cell_idx] = new_cell
    elif mode in ("insert_before", "insert_after"):
        if not (0 <= cell_idx <= n_before):
            raise ValueError(
                f"cell_idx {cell_idx} out of range for insert "
                f"(notebook has {n_before} cells; valid 0..{n_before})"
            )
        new_cell = {
            "cell_type": cell_type,
            "metadata": {},
            "source": _string_to_source(source or ""),
        }
        if cell_type == "code":
            new_cell["outputs"] = []
            new_cell["execution_count"] = None
        # Generate a short cell id (nbformat 4.5+).
        import secrets
        new_cell["id"] = secrets.token_hex(4)
        insert_at = cell_idx if mode == "insert_before" else cell_idx + 1
        cells.insert(insert_at, new_cell)

    nb["cells"] = cells
    # Preserve trailing newline convention.
    text = json.dumps(nb, indent=1, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")
    return n_before, len(cells)
