"""The Submit tab's GOAT-method selector: CONTROL text <-> choice, and availability.

GOAT optimises with the ``xTB_method`` unless CONTROL.txt names an orca-mlips
machine-learned potential in ``global_opt_method=mlip:<plugin>`` (validated by
delfin.common.control_validator); ``global_opt_mlip_model`` optionally names
the model the plugin loads.  Only GOAT's method changes: every later step of
the CONTROL runs as before.  An MLIP GOAT runs on the GPU of the node the
DELFIN job runs on, so the submission sends the whole job to a GPU partition
(delfin.mlp_tools.goat_gpu_submit).

The text is the single source of truth.  :func:`read_method` reads it,
:func:`write_method` / :func:`write_model` change exactly the lines they own
(and ``global_optimizer`` when an MLIP needs GOAT selected), so the selector
and the textarea can follow each other both ways.  Nothing here imports
ipywidgets, torch or an MLIP stack.
"""

from __future__ import annotations

import html
import re
import time
from typing import Mapping, Optional

from delfin.config import set_control_value

METHOD_KEY = "global_opt_method"
MODEL_KEY = "global_opt_mlip_model"
OPTIMIZER_KEY = "global_optimizer"

#: The selector's value for the classic GOAT with the xTB_method.
XTB = ""
#: The selector's value while the text holds something it cannot show.
UNKNOWN = "?"

#: value, short name, one line: what it is, licence, GPU.
METHODS = (
    (XTB, "xTB (as xTB_method)",
     "GOAT with the semi-empirical xTB_method of this CONTROL (GFN2-xTB by default); "
     "LGPL-3.0, runs on CPUs."),
    ("mlip:aimnet2", "AIMNet2",
     "AIMNet2 neural-network potential (14 main-group elements, neutral and charged "
     "molecules); weights MIT; needs a GPU node."),
    ("mlip:orb", "ORB",
     "ORB v3 conservative OMol potential from Orbital Materials (broad chemistry); "
     "weights Apache-2.0; needs a GPU node."),
    ("mlip:mace", "MACE",
     "MACE-OMOL-0 equivariant foundation potential; weights under the Academic Software "
     "License, non-commercial use only; needs a GPU node."),
    ("mlip:uma", "UMA",
     "Meta FAIR Chemistry UMA universal potential (OMol task); weights gated on "
     "HuggingFace (your own login + accepted terms); needs a GPU node."),
)

_BY_VALUE = {value: (name, line) for value, name, line in METHODS}

_LINE = r"(?m)^[ \t]*{key}[ \t]*=(.*)$"


def _values(text: str, key: str) -> list:
    return [m.group(1).strip() for m in re.finditer(_LINE.format(key=re.escape(key)), text or "")]


def _value(text: str, key: str) -> Optional[str]:
    """The value of *key* as DELFIN's parser takes it (the last line wins), or None."""
    found = _values(text, key)
    return found[-1] if found else None


def read_method(text: str) -> tuple:
    """``(selector value, raw value)`` of the CONTROL text's GOAT method.

    An absent or empty ``global_opt_method`` (and every spelling the validator
    reads as "xTB": none, no, xtb ...) is :data:`XTB`; ``mlip:<plugin>`` in any
    spelling the validator accepts is that plugin; anything else is
    :data:`UNKNOWN` with its raw text, which the selector shows as such and
    never replaces.
    """
    from delfin.common.control_validator import _as_global_opt_method

    raw = _value(text, METHOD_KEY)
    if raw is None:
        return XTB, ""
    try:
        value = _as_global_opt_method(raw)
    except ValueError:
        return UNKNOWN, raw
    return (value if value in _BY_VALUE else UNKNOWN), raw


def read_model(text: str) -> str:
    """The ``global_opt_mlip_model`` value, "" when absent."""
    return _value(text, MODEL_KEY) or ""


def _insert_after(text: str, anchors: tuple, line: str) -> str:
    """*text* with *line* on its own line after the last of the first anchor key found.

    Appended at the end when no anchor key has a line.  Nothing else changes.
    """
    for key in anchors:
        matches = list(re.finditer(_LINE.format(key=re.escape(key)), text))
        if matches:
            end = matches[-1].end()
            return text[:end] + "\n" + line + text[end:]
    body = text.rstrip("\n")
    return (body + "\n" if body else "") + line + "\n"


def _goat_selected(text: str) -> bool:
    """Whether the text runs GOAT (``global_optimizer``, else the older XTB_GOAT)."""
    from delfin.common.control_validator import _as_global_optimizer

    chosen = _value(text, OPTIMIZER_KEY)
    if chosen is not None:
        try:
            return _as_global_optimizer(chosen) == "GOAT"
        except ValueError:
            return False
    return str(_value(text, "XTB_GOAT") or "").strip().lower() in ("yes", "true", "1", "on")


def write_method(text: str, value: str) -> str:
    """*text* with GOAT's method set to *value* (:data:`XTB` or ``mlip:<plugin>``).

    Rewrites only the ``global_opt_method`` line (inserted after
    ``global_optimizer`` when there is none).  An MLIP also selects GOAT the
    way the Only-GOAT template does (``global_optimizer=GOAT``) when the text
    does not already run GOAT -- CREST runs xTB only.  xTB on a text without
    the line leaves the text as it is.
    """
    text = text or ""
    if value not in _BY_VALUE:
        raise ValueError(f"not a GOAT method of the selector: {value!r}")
    has_line = _value(text, METHOD_KEY) is not None
    if value == XTB:
        return set_control_value(text, METHOD_KEY, "") if has_line else text
    if not _goat_selected(text):
        if _value(text, OPTIMIZER_KEY) is not None:
            text = set_control_value(text, OPTIMIZER_KEY, "GOAT")
        else:
            text = _insert_after(text, ("xTB_method", "XTB_preOPT"), f"{OPTIMIZER_KEY}=GOAT")
    if has_line:
        return set_control_value(text, METHOD_KEY, value)
    return _insert_after(text, (OPTIMIZER_KEY, "XTB_GOAT", "xTB_method"), f"{METHOD_KEY}={value}")


def write_model(text: str, model: str) -> str:
    """*text* with ``global_opt_mlip_model`` set to *model* ("" = the plugin's default).

    An empty model on a text without the line leaves the text as it is.
    """
    text = text or ""
    model = str(model or "").strip()
    if _value(text, MODEL_KEY) is not None:
        return set_control_value(text, MODEL_KEY, model)
    if not model:
        return text
    return _insert_after(text, (METHOD_KEY, OPTIMIZER_KEY), f"{MODEL_KEY}={model}")


def model_error(model: str) -> str:
    """The validator's complaint about a model name, "" when it is fine."""
    from delfin.common.control_validator import _as_global_opt_mlip_model

    try:
        _as_global_opt_mlip_model(model)
    except ValueError as exc:
        return str(exc)
    return ""


# ---- availability -----------------------------------------------------------

#: How long one availability reading is kept: the selector is built with every
#: tab at start and asked again on each change; an install in between shows
#: after this many seconds (or at the next choice after them).
STATUS_TTL = 60.0
_CACHE: dict = {"at": None, "status": None}


def availability(*, refresh: bool = False, environ: Optional[Mapping[str, str]] = None) -> dict:
    """:func:`delfin.mlp_tools.orca_mlips.mlips_backend_status`, kept for :data:`STATUS_TTL`.

    Metadata only (directory listings; no torch, nothing started).  A failure
    to read it is an empty dict: the selector then shows no availability
    rather than not being built.
    """
    now = time.monotonic()
    if (not refresh and environ is None and _CACHE["status"] is not None
            and now - _CACHE["at"] < STATUS_TTL):
        return _CACHE["status"]
    try:
        from delfin.mlp_tools.orca_mlips import mlips_backend_status

        status = mlips_backend_status(environ)
    except Exception:
        status = {}
    if environ is None:
        _CACHE.update(at=now, status=status)
    return status


def _plugin(value: str) -> str:
    return value.split(":", 1)[1] if value.startswith("mlip:") else ""


def short_state(value: str, status: Mapping) -> str:
    """A few words on whether *value* can run here, for the option label."""
    plugin = _plugin(value)
    if not plugin:
        return ""
    entry = status.get(plugin)
    if not entry:
        return ""
    if not entry.get("installed"):
        return "not installed"
    if plugin == "uma" and not entry.get("hf_login"):
        return "installed, HuggingFace login needed"
    return "installed"


def options(status: Mapping, unknown_raw: Optional[str] = None) -> list:
    """``[(label, value), ...]`` for the selector, availability next to each MLIP."""
    out = []
    for value, name, _line in METHODS:
        state = short_state(value, status)
        out.append((f"{name} - {state}" if state else name, value))
    if unknown_raw is not None:
        out.append((f"Unknown: {unknown_raw or '(empty)'}", UNKNOWN))
    return out


def tooltip() -> str:
    """The selector's tooltip: every choice in one line each."""
    lines = ["Method GOAT optimises with (CONTROL key global_opt_method)."]
    lines += [f"{name}: {line}" for _value, name, line in METHODS]
    lines.append("Only GOAT changes; every later step of the CONTROL runs as before. "
                 "An MLIP GOAT sends the whole DELFIN job to a GPU partition.")
    return "\n".join(lines)


def status_html(value: str, status: Mapping, raw: str = "", model_problem: str = "") -> str:
    """The line under the selector: what the choice is and what it needs here."""
    esc = html.escape
    if value == UNKNOWN:
        return ('<span style="color:#b71c1c;">global_opt_method=' + esc(raw)
                + ' is not a method this selector knows (empty = xTB, or mlip:aimnet2, '
                  'mlip:orb, mlip:mace, mlip:uma); VALIDATE CONTROL says what is wrong. '
                  'The text is left as it is.</span>')
    name, line = _BY_VALUE.get(value, ("", ""))
    parts = [f"<b>{esc(name)}</b>: {esc(line)}"]
    plugin = _plugin(value)
    if plugin:
        entry = status.get(plugin) or {}
        if entry.get("installed"):
            parts.append(
                '<span style="color:#2e7d32;">Installed</span> (orca-mlips '
                f"{esc(str(entry.get('orca_mlips_version') or ''))}, "
                f"{esc(str(entry.get('stack')))} {esc(str(entry.get('stack_version') or ''))}).")
        elif entry:
            parts.append(
                '<span style="color:#b71c1c;">Not installed</span> '
                f"({esc(str(entry.get('missing') or ''))}). Install with "
                f"<code>{esc(str(entry.get('install') or ''))}</code>.")
        if plugin == "uma" and entry:
            if entry.get("hf_login"):
                parts.append("HuggingFace login found; access to the gated model is "
                             "checked on its first download.")
            else:
                parts.append("No HuggingFace login found: run <code>hf auth login</code> "
                             "(older huggingface_hub: <code>huggingface-cli login</code>) "
                             "in the UMA environment and accept the model terms at "
                             '<a href="https://huggingface.co/facebook/UMA" target="_blank" '
                             'rel="noopener">huggingface.co/facebook/UMA</a>.')
        lic = entry.get("licence") or {}
        if lic.get("non_commercial"):
            parts.append(f"Weights: {esc(str(lic.get('name')))} - academic use only.")
        model = entry.get("default_model")
        if model:
            parts.append(f"Model: empty = {esc(str(model))}.")
        parts.append("Only GOAT changes; later steps run as before. An MLIP GOAT sends "
                     "the whole DELFIN job to a GPU partition (--gres=gpu:1).")
    if model_problem:
        parts.append(f'<span style="color:#b71c1c;">Model: {esc(model_problem)}</span>')
    return '<span style="font-size:12px;">' + " ".join(parts) + "</span>"
