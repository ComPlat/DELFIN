"""Application keys files — a CONTROL.txt-style input file for an application.

Collect a workflow's input keys in one editable text file, adjust them, and run
it — exactly like CONTROL.txt, but generated from the application's input
contract so every key shows its default, type and (from the central key
vocabulary) its allowed values.

    delfin-app template redox_potential > my_redox.txt   # generate
    #  ... edit my_redox.txt ...
    delfin-app run my_redox.txt                          # run

Format: ``key=value`` lines; full-line ``#`` comments (so values may contain
``#``, e.g. SMILES triple bonds).  One reserved key ``application=<name>`` names
the application.  Empty values fall back to the contract default.
"""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Tuple

from delfin.tools._application import Application, get_application


def application_keyfile(name: str) -> str:
    """Generate an editable CONTROL.txt-style keys file for an application."""
    app = get_application(name)
    if app is None:
        raise ValueError(f"unknown application {name!r}")

    lines: List[str] = [
        f"# DELFIN application keys: {app.name}  (v{app.version})",
    ]
    if app.description:
        lines.append(f"# {app.description}")
    lines += [
        "# Edit the values below, then run:  delfin-app run <this-file>",
        "# '*' = required;  full-line '#' = comment;  empty value = use default.",
        "",
        f"application={app.name}",
        "",
    ]

    for p in app.inputs:
        bits: List[str] = []
        if p.required:
            bits.append("*")
        if p.description:
            bits.append(p.description)
        bits.append(f"({p.type})")
        if p.unit:
            bits.append(f"[{p.unit}]")
        if p.enum:
            allowed = list(p.enum)
            shown = ", ".join(str(v) for v in allowed[:6])
            more = f", … (+{len(allowed) - 6} more)" if len(allowed) > 6 else ""
            bits.append(f"— allowed: {shown}{more}")
        lines.append(f"# {p.name}: {' '.join(bits)}")
        default = "" if p.default is None else p.default
        lines.append(f"{p.name}={default}")
        lines.append("")

    if app.outputs:
        outs = ", ".join(o.name for o in app.outputs)
        lines.append(f"# outputs: {outs}")
    return "\n".join(lines) + "\n"


def parse_keyfile(text: str) -> Tuple[Optional[str], Dict[str, str]]:
    """Parse a keys file → (application name, {key: raw string value}).

    Full-line ``#`` comments and blanks are skipped; values keep ``#`` so SMILES
    work.  Empty values are dropped (the contract default applies).
    """
    app_name: Optional[str] = None
    values: Dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key, value = key.strip(), value.strip()
        if key == "application":
            app_name = value or None
            continue
        if value == "":
            continue
        values[key] = value
    return app_name, values


def _coerce(value: str, type_: str) -> Any:
    if type_ == "int":
        return int(value)
    if type_ == "float":
        return float(value)
    if type_ == "bool":
        return value.strip().lower() in ("1", "true", "yes", "on")
    if type_ in ("list", "dict"):
        try:
            return json.loads(value)
        except (ValueError, TypeError):
            return [v.strip() for v in value.split(",")] if type_ == "list" else value
    return value  # str / path


def coerce_to_contract(app: Application, raw_values: Dict[str, str]) -> Dict[str, Any]:
    """Coerce raw string values to the types declared by the app's inputs."""
    types = {p.name: p.type for p in app.inputs}
    out: Dict[str, Any] = {}
    for key, value in raw_values.items():
        out[key] = _coerce(value, types.get(key, "str"))
    return out


def inputs_from_keyfile(text: str) -> Tuple[str, Dict[str, Any]]:
    """Resolve a keys file to (application name, typed inputs dict)."""
    app_name, raw = parse_keyfile(text)
    if not app_name:
        raise ValueError("keys file has no 'application=<name>' line")
    app = get_application(app_name)
    if app is None:
        raise ValueError(f"unknown application {app_name!r}")
    return app_name, coerce_to_contract(app, raw)


def run_from_keyfile(path: str, *, cores: int = 1, submit: bool = False, **overrides: Any):
    """Run (or submit) an application described by a keys file.

    ``submit=False`` runs synchronously and returns an ApplicationResult;
    ``submit=True`` submits to the runtime and returns a run id.  *overrides*
    take precedence over the file's values.
    """
    from pathlib import Path

    from delfin.tools import platform

    text = Path(path).read_text(encoding="utf-8")
    app_name, inputs = inputs_from_keyfile(text)
    inputs.update(overrides)
    if submit:
        return platform.submit_application(app_name, cores=cores, **inputs)
    return platform.run_application(app_name, cores=cores, **inputs)


__all__ = [
    "application_keyfile",
    "parse_keyfile",
    "coerce_to_contract",
    "inputs_from_keyfile",
    "run_from_keyfile",
]
