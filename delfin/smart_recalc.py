"""Smart recalculation: skip ORCA when inp+dependencies are unchanged and output is complete.

Uses a SHA-256 fingerprint of:
 - The ORCA input file content (after all overrides have been applied)
 - Referenced dependency files from the input (e.g. ``%moinp``, ``* xyzfile``,
   Hessian references in ESD blocks, and other quoted file references)
 - Optional extra dependencies passed by the caller (e.g. ``copy_files``)

The fingerprint is persisted as a ``.fprint`` sidecar file next to the input file.

Modes:
 - Smart recalc (default): ``DELFIN_SMART_RECALC=1`` (or unset)
 - Classic recalc: ``DELFIN_SMART_RECALC=0`` (skip by marker only)

On ``delfin --recalc`` (``DELFIN_RECALC=1``), ``should_skip()`` returns:
 - Smart mode: True only when fingerprint matches and output is complete.
 - Classic mode: True when output is complete.
"""

import hashlib
import logging
import os
import re
from pathlib import Path
from typing import Iterable, List, Optional

logger = logging.getLogger(__name__)

OK_MARKER = "ORCA TERMINATED NORMALLY"

# Explicit dependency patterns commonly used in ORCA/DELFIN inputs.
_MOINP_RE = re.compile(
    r'^\s*%?moinp\s+(?:"([^"]+)"|(\S+))',
    re.IGNORECASE | re.MULTILINE,
)
_XYZFILE_RE = re.compile(
    r'^\s*\*\s*xyzfile\s+\S+\s+\S+\s+(?:"([^"]+)"|(\S+))',
    re.IGNORECASE | re.MULTILINE,
)
_HESS_DIRECTIVE_RE = re.compile(
    r'^\s*(?:gshessian|eshessian|tshessian|iscishess|iscfshess|inhess(?:ian)?|inhessname)\s+'
    r'(?:"([^"]+)"|(\S+))',
    re.IGNORECASE | re.MULTILINE,
)

# Broad fallback: quoted tokens that look like file paths.
_QUOTED_TOKEN_RE = re.compile(r'"([^"\n]+)"')
_PATHLIKE_EXT_RE = re.compile(r'.+\.[A-Za-z0-9]{1,12}$')


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def recalc_enabled() -> bool:
    """Return True when the ``DELFIN_RECALC`` environment variable is set."""
    return str(os.environ.get("DELFIN_RECALC", "0")).lower() in (
        "1", "true", "yes", "on", "y"
    )


def smart_mode_enabled() -> bool:
    """Return True when smart dependency-aware recalc logic is enabled."""
    return str(os.environ.get("DELFIN_SMART_RECALC", "1")).lower() not in (
        "0", "false", "no", "off", "n"
    )


def has_ok_marker(out_path) -> bool:
    """Return True when *out_path* exists and contains ``ORCA TERMINATED NORMALLY``."""
    try:
        return OK_MARKER in Path(out_path).read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Fingerprint logic
# ---------------------------------------------------------------------------

def _looks_like_file_reference(token: str) -> bool:
    """Heuristic filter for path-like tokens found in ORCA input text."""
    value = str(token or "").strip().strip('"').strip("'").rstrip(",;")
    if not value:
        return False
    low = value.lower()
    if low in {"true", "false", "none", "null", "auto"}:
        return False
    if "/" in value or "\\" in value:
        return True
    return bool(_PATHLIKE_EXT_RE.match(value))


def _resolve_dep_path(raw_dep: str, inp_path: Path) -> Path:
    dep = Path(str(raw_dep).strip().strip('"').strip("'").rstrip(",;"))
    if not dep.is_absolute():
        dep = inp_path.parent / dep
    try:
        return dep.resolve()
    except Exception:
        return dep


def _collect_input_dependencies(inp_path: Path) -> List[Path]:
    """Return resolved absolute dependency paths referenced by the ORCA input."""
    try:
        text = inp_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return []

    raw_candidates: list[str] = []
    for match in _MOINP_RE.finditer(text):
        raw = match.group(1) or match.group(2)
        if raw:
            raw_candidates.append(raw)
    for match in _XYZFILE_RE.finditer(text):
        raw = match.group(1) or match.group(2)
        if raw:
            raw_candidates.append(raw)
    for match in _HESS_DIRECTIVE_RE.finditer(text):
        raw = match.group(1) or match.group(2)
        if raw:
            raw_candidates.append(raw)
    for match in _QUOTED_TOKEN_RE.finditer(text):
        raw = match.group(1)
        if raw:
            raw_candidates.append(raw)

    deps: list[Path] = []
    seen: set[str] = set()
    for raw in raw_candidates:
        if not _looks_like_file_reference(raw):
            continue
        dep = _resolve_dep_path(raw, inp_path)
        key = str(dep)
        if key in seen:
            continue
        seen.add(key)
        deps.append(dep)
    return deps


def _normalize_extra_dependencies(inp_path: Path, extra_deps: Optional[Iterable]) -> List[Path]:
    if not extra_deps:
        return []
    deps: list[Path] = []
    seen: set[str] = set()
    for raw in extra_deps:
        if raw is None:
            continue
        dep = Path(raw)
        if not dep.is_absolute():
            dep = inp_path.parent / dep
        try:
            dep = dep.resolve()
        except Exception:
            pass
        key = str(dep)
        if key in seen:
            continue
        seen.add(key)
        deps.append(dep)
    return deps


def compute_fingerprint(inp_path, extra_deps: Optional[Iterable] = None) -> str:
    """Return a hex-digest fingerprint for *inp_path* and its dependencies.

    The fingerprint covers:
    - The full byte content of the input file.
    - For each dependency: ``<resolved_path>:<size>:<mtime_ns>``.
      (input-referenced dependencies + optional ``extra_deps``).

    Returns an empty string when the input file cannot be read.
    """
    inp = Path(inp_path)
    h = hashlib.sha256()
    try:
        h.update(inp.read_bytes())
    except Exception:
        return ""

    deps: list[Path] = []
    deps.extend(_collect_input_dependencies(inp))
    deps.extend(_normalize_extra_dependencies(inp, extra_deps))

    # Stable ordering independent of discovery order.
    unique_sorted = sorted({str(dep): dep for dep in deps}.values(), key=lambda p: str(p))

    for dep in unique_sorted:
        dep_id = str(dep)
        try:
            st = dep.stat()
            h.update(f"{dep_id}:{st.st_size}:{st.st_mtime_ns}".encode())
        except Exception:
            # Dependency missing — include a marker so the fingerprint differs
            # from a previous run where it existed.
            h.update(f"{dep_id}:missing".encode())
    return h.hexdigest()


def _fprint_path(inp_path: Path) -> Path:
    return inp_path.with_suffix(inp_path.suffix + ".fprint")


def fingerprint_unchanged(inp_path, extra_deps: Optional[Iterable] = None) -> bool:
    """Return True when the stored ``.fprint`` sidecar matches the current fingerprint."""
    inp = Path(inp_path)
    fp_file = _fprint_path(inp)
    if not fp_file.exists():
        return False
    try:
        stored = fp_file.read_text(encoding="utf-8").strip()
    except Exception:
        return False
    current = compute_fingerprint(inp, extra_deps=extra_deps)
    return bool(current) and stored == current


def store_fingerprint(inp_path, extra_deps: Optional[Iterable] = None) -> None:
    """Persist the current fingerprint of *inp_path* as a ``.fprint`` sidecar file."""
    inp = Path(inp_path)
    fp = compute_fingerprint(inp, extra_deps=extra_deps)
    if not fp:
        return
    try:
        _fprint_path(inp).write_text(fp + "\n", encoding="utf-8")
    except Exception as exc:
        logger.debug("[smart_recalc] Could not store fingerprint for %s: %s", inp, exc)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def should_skip(inp_path, out_path, extra_deps: Optional[Iterable] = None) -> bool:
    """Return True when ORCA execution can safely be skipped.

    All three conditions must hold:
    1. ``DELFIN_RECALC`` is enabled in the environment.
    2. *out_path* exists and contains ``ORCA TERMINATED NORMALLY``.
    3. Smart mode only: the ``.fprint`` sidecar of *inp_path* matches the
       current fingerprint (i.e. neither inp content nor dependencies changed).
    """
    if not recalc_enabled():
        return False
    if not has_ok_marker(out_path):
        return False
    if not smart_mode_enabled():
        return True
    if not fingerprint_unchanged(inp_path, extra_deps=extra_deps):
        return False
    return True
