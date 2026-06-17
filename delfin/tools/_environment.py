"""Environment probing + tool-installation policy for the platform.

Two jobs, both driven by the building blocks' declared ``requires_binaries`` /
``requires_python`` contracts:

1. **probe** — which tools each capability needs and whether they are present, so
   settings / the agent can show readiness and route around missing tools.
2. **install** — install the *open-source* tools DELFIN may legally ship an
   installer for, and for *license-restricted* tools (ORCA, Turbomole) emit a
   clear hint on how to obtain them and where to put them so DELFIN finds them —
   **never** attempting to install or redistribute them.

The install **policy** is per tool:

* ``auto``   — open source; handled by DELFIN's bundled ``install_qm_tools.sh``
  (binaries) or pip/conda (Python).
* ``manual`` — license-restricted or commercial; DELFIN must not install it.
  The user installs it themselves and makes it discoverable (PATH / runtime
  settings); the platform only provides guidance.
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delfin.tools._registry import list_steps


# --- tool metadata / install policy ---------------------------------------


@dataclass(frozen=True)
class ToolInfo:
    name: str
    kind: str                 # "binary" | "python"
    policy: str               # "auto" | "manual"
    source: str = ""          # official download / info URL
    install_hint: str = ""    # how to install / where to place the file
    installer: str = ""       # "qm_tools" | "pip" | ""


_KNOWN: Dict[str, ToolInfo] = {}


def _t(name: str, kind: str, policy: str, *, source: str = "",
       hint: str = "", installer: str = "") -> None:
    _KNOWN[name.lower()] = ToolInfo(name, kind, policy, source, hint, installer)


# License-restricted — DELFIN must NOT install or redistribute these.
_t(
    "orca", "binary", "manual",
    source="https://orcaforum.kit.edu",
    hint=(
        "ORCA is free for academic use after registering on the ORCA forum, but "
        "may not be redistributed — DELFIN cannot install it. Download the ORCA "
        "binaries for your OS yourself, then either add the ORCA directory to PATH "
        "or set the ORCA base in DELFIN's runtime settings; DELFIN auto-discovers "
        "installations (see runtime_setup.discover_orca_installations)."
    ),
)
for _tm in ("turbomole", "define", "x2t", "ridft", "dscf", "jobex", "aoforce"):
    _t(
        _tm, "binary", "manual",
        source="https://www.turbomole.org",
        hint=(
            "Turbomole is commercial, licensed software and cannot be installed by "
            "DELFIN. Install it from your licensed distribution and put its binaries "
            "(define, x2t, ridft, …) on PATH."
        ),
    )

# Open-source QM binaries — installable via DELFIN's bundled installer.
for _b in ("xtb", "crest", "xtb4stda", "std2", "stda", "dftb+"):
    _t(
        _b, "binary", "auto",
        source="https://github.com/grimme-lab",
        hint="Open source; install via DELFIN (bundled install_qm_tools.sh) or your package manager.",
        installer="qm_tools",
    )

# Common Python dependencies — pip/conda installable.
for _m in ("rdkit", "mendeleev", "morfeus", "cclib", "ase", "openbabel"):
    _t(_m, "python", "auto", hint=f"pip install {_m}  (or via conda)", installer="pip")


def tool_info(name: str, kind: str) -> ToolInfo:
    """Return known metadata for *name*, or a safe default.

    Unknown binaries default to ``manual`` (DELFIN will not assume it may
    auto-install an arbitrary binary); unknown Python modules default to ``auto``
    with a pip hint.
    """
    info = _KNOWN.get(name.lower())
    if info is not None:
        return info
    if kind == "python":
        return ToolInfo(name, "python", "auto", install_hint=f"pip install {name}", installer="pip")
    return ToolInfo(
        name, "binary", "manual",
        install_hint=(f"'{name}' is required but was not found on PATH. Install it and add it "
                      f"to PATH, or place it where DELFIN's qm_tools resolver looks."),
    )


# --- probing --------------------------------------------------------------


@dataclass
class RequirementStatus:
    name: str
    kind: str                 # "binary" | "python"
    available: bool
    detail: str = ""
    policy: str = "auto"
    source: str = ""
    install_hint: str = ""


@dataclass
class CapabilityReadiness:
    step_name: str
    ready: bool
    missing: List[RequirementStatus] = field(default_factory=list)


def _binary_available(name: str) -> Tuple[bool, str]:
    found = shutil.which(name)
    if found:
        return True, found
    # Best-effort: a binary staged under the packaged qm_tools bin dir.
    try:
        from delfin.runtime_setup import get_packaged_qm_tools_dir
        cand = get_packaged_qm_tools_dir() / "bin" / name
        if cand.exists():
            return True, str(cand)
    except Exception:
        pass
    return False, "not found on PATH"


def _python_available(name: str) -> Tuple[bool, str]:
    try:
        ok = importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        ok = False
    return ok, ("importable" if ok else "not importable")


def _status(name: str, kind: str) -> RequirementStatus:
    available, detail = (_binary_available(name) if kind == "binary"
                         else _python_available(name))
    if available:
        return RequirementStatus(name, kind, True, detail)
    info = tool_info(name, kind)
    return RequirementStatus(
        name, kind, False, detail,
        policy=info.policy, source=info.source, install_hint=info.install_hint,
    )


def probe() -> List[CapabilityReadiness]:
    """Readiness of every registered capability against its declared requirements."""
    out: List[CapabilityReadiness] = []
    for name, adapter in sorted(list_steps().items()):
        c = adapter.contract()
        missing: List[RequirementStatus] = []
        for b in sorted(c.requires_binaries):
            s = _status(b, "binary")
            if not s.available:
                missing.append(s)
        for m in sorted(c.requires_python):
            s = _status(m, "python")
            if not s.available:
                missing.append(s)
        out.append(CapabilityReadiness(name, ready=not missing, missing=missing))
    return out


def missing_tools() -> Dict[str, Dict[str, List[str]]]:
    """Unique missing requirements → the capabilities that need them.

    Returns ``{"binaries": {name: [caps]}, "python": {name: [caps]}}``.
    """
    binaries: Dict[str, List[str]] = {}
    python: Dict[str, List[str]] = {}
    for cr in probe():
        for r in cr.missing:
            bucket = binaries if r.kind == "binary" else python
            bucket.setdefault(r.name, []).append(cr.step_name)
    return {"binaries": binaries, "python": python}


# --- installation ---------------------------------------------------------


def installer_path() -> Optional[Path]:
    """Path to the bundled open-source QM-tools installer, if present."""
    try:
        from delfin.runtime_setup import get_packaged_qm_tools_dir
        p = get_packaged_qm_tools_dir() / "install_qm_tools.sh"
        return p if p.is_file() else None
    except Exception:
        return None


def install_plan() -> Dict[str, object]:
    """Classify missing tools into auto-installable vs. manual (license-restricted).

    Never schedules ORCA/Turbomole (or any ``manual`` tool) for installation —
    those come back under ``manual`` with a hint on how to obtain and place them.
    """
    miss = missing_tools()
    auto_binaries: List[str] = []
    manual: List[Dict[str, str]] = []
    python: List[Dict[str, str]] = []

    for name in sorted(miss["binaries"]):
        info = tool_info(name, "binary")
        if info.policy == "manual":
            manual.append({"tool": name, "source": info.source, "hint": info.install_hint})
        else:
            auto_binaries.append(name)
    for name in sorted(miss["python"]):
        info = tool_info(name, "python")
        if info.policy == "manual":
            manual.append({"tool": name, "source": info.source, "hint": info.install_hint})
        else:
            python.append({"tool": name, "hint": info.install_hint or f"pip install {name}"})

    return {
        "auto_binaries": auto_binaries,      # handled by install_qm_tools.sh
        "manual": manual,                    # license-restricted: never auto-installed
        "python": python,                    # pip/conda — surfaced as guidance
        "installer": str(installer_path() or ""),
    }


def install_tools(*, run: bool = False, timeout: Optional[float] = None) -> Dict[str, object]:
    """Install the *open-source* QM binaries DELFIN may legally install.

    With ``run=False`` (default) returns the plan without executing anything.
    With ``run=True`` runs the bundled ``install_qm_tools.sh`` for the missing
    open-source binaries only.  License-restricted tools (ORCA, Turbomole) are
    **never** installed; their guidance is always returned under ``plan.manual``.
    """
    plan = install_plan()
    out: Dict[str, object] = {"plan": plan, "executed": False}
    if not run:
        return out

    script = installer_path()
    if script is None or not plan["auto_binaries"]:
        out["note"] = "no auto-installable binaries missing, or installer not found"
        return out
    try:
        proc = subprocess.run(
            ["bash", str(script)], capture_output=True, text=True, timeout=timeout,
        )
        out["executed"] = True
        out["ok"] = proc.returncode == 0
        out["returncode"] = proc.returncode
        out["stdout"] = (proc.stdout or "")[-4000:]
        out["stderr"] = (proc.stderr or "")[-4000:]
    except Exception as exc:  # never raise out of the platform call
        out["ok"] = False
        out["error"] = str(exc)
    return out


__all__ = [
    "ToolInfo",
    "tool_info",
    "RequirementStatus",
    "CapabilityReadiness",
    "probe",
    "missing_tools",
    "installer_path",
    "install_plan",
    "install_tools",
]
