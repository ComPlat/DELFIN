"""Smart recalculation: skip ORCA when inp+dependencies are unchanged and output is complete.

Uses a SHA-256 fingerprint of:
 - The ORCA input file content (after all overrides have been applied)
 - The name, size, and modification time of each GBW file referenced via ``%moinp``

The fingerprint is persisted as a ``.fprint`` sidecar file next to the input file.
On ``delfin --recalc`` (``DELFIN_RECALC=1``), ``should_skip()`` returns *True* only when
the stored fingerprint matches the current one **and** the output contains
``ORCA TERMINATED NORMALLY``.
"""

import hashlib
import logging
import os
import re
from pathlib import Path

logger = logging.getLogger(__name__)

OK_MARKER = "ORCA TERMINATED NORMALLY"

# Matches: %moinp "file.gbw"  (case-insensitive, optional surrounding whitespace)
_MOINP_RE = re.compile(r'%moinp\s+"([^"]+)"', re.IGNORECASE)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def recalc_enabled() -> bool:
    """Return True when the ``DELFIN_RECALC`` environment variable is set."""
    return str(os.environ.get("DELFIN_RECALC", "0")).lower() in (
        "1", "true", "yes", "on", "y"
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

def _gbw_deps(inp_path: Path) -> list:
    """Return resolved absolute paths of GBW files referenced via ``%moinp``."""
    try:
        text = inp_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return []
    deps = []
    for m in _MOINP_RE.finditer(text):
        dep = Path(m.group(1))
        if not dep.is_absolute():
            dep = inp_path.parent / dep
        try:
            deps.append(dep.resolve())
        except Exception:
            deps.append(dep)
    return deps


def compute_fingerprint(inp_path) -> str:
    """Return a hex-digest fingerprint for *inp_path* and its GBW dependencies.

    The fingerprint covers:
    - The full byte content of the input file.
    - For each ``%moinp`` dependency: ``<name>:<size>:<mtime_ns>``.

    Returns an empty string when the input file cannot be read.
    """
    inp = Path(inp_path)
    h = hashlib.sha256()
    try:
        h.update(inp.read_bytes())
    except Exception:
        return ""
    for dep in _gbw_deps(inp):
        try:
            st = dep.stat()
            h.update(f"{dep.name}:{st.st_size}:{st.st_mtime_ns}".encode())
        except Exception:
            # Dependency missing — include a marker so the fingerprint differs
            # from a previous run where it existed.
            h.update(f"{dep.name}:missing".encode())
    return h.hexdigest()


def _fprint_path(inp_path: Path) -> Path:
    return inp_path.with_suffix(inp_path.suffix + ".fprint")


def fingerprint_unchanged(inp_path) -> bool:
    """Return True when the stored ``.fprint`` sidecar matches the current fingerprint."""
    inp = Path(inp_path)
    fp_file = _fprint_path(inp)
    if not fp_file.exists():
        return False
    try:
        stored = fp_file.read_text(encoding="utf-8").strip()
    except Exception:
        return False
    current = compute_fingerprint(inp)
    return bool(current) and stored == current


def store_fingerprint(inp_path) -> None:
    """Persist the current fingerprint of *inp_path* as a ``.fprint`` sidecar file."""
    inp = Path(inp_path)
    fp = compute_fingerprint(inp)
    if not fp:
        return
    try:
        _fprint_path(inp).write_text(fp + "\n", encoding="utf-8")
    except Exception as exc:
        logger.debug("[smart_recalc] Could not store fingerprint for %s: %s", inp, exc)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def should_skip(inp_path, out_path) -> bool:
    """Return True when ORCA execution can safely be skipped.

    All three conditions must hold:
    1. ``DELFIN_RECALC`` is enabled in the environment.
    2. *out_path* exists and contains ``ORCA TERMINATED NORMALLY``.
    3. The ``.fprint`` sidecar of *inp_path* matches the current fingerprint
       (i.e. neither the inp content nor any referenced GBW dependency changed).
    """
    if not recalc_enabled():
        return False
    if not has_ok_marker(out_path):
        return False
    if not fingerprint_unchanged(inp_path):
        return False
    return True
