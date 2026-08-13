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
_BASE_DIRECTIVE_RE = re.compile(
    r'^\s*%base\s+(?:"([^"]+)"|(\S+))',
    re.IGNORECASE | re.MULTILINE,
)
_ORCA_BANG_RE = re.compile(r'^\s*!', re.MULTILINE)
_ORCA_BLOCK_RE = re.compile(r'^\s*%(?:pal|scf|base|moinp)\b', re.IGNORECASE | re.MULTILINE)
_ORCA_XYZFILE_RE = re.compile(r'^\s*\*\s*xyzfile\b', re.IGNORECASE | re.MULTILINE)
_SPECIAL_REQUIRED_OUTPUTS = {
    "xtb.inp": ("XTB.xyz",),
    "xtb_goat.inp": ("XTB_GOAT.globalminimum.xyz",),
    "xtb_solvator.inp": ("XTB_SOLVATOR.solvator.xyz",),
}


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
    """Return True when *out_path* exists and contains ``ORCA TERMINATED NORMALLY``.

    Only reads the last 8 KB of the file since the marker always appears
    near the end, avoiding full reads of large (100–500 MB) output files.
    """
    try:
        p = Path(out_path)
        size = p.stat().st_size
        with p.open(encoding="utf-8", errors="replace") as f:
            if size > 8192:
                f.seek(size - 8192)
            return OK_MARKER in f.read()
    except Exception:
        return False


def _looks_like_orca_job(inp_path: Path) -> bool:
    """Heuristically detect ORCA inputs without affecting xTB marker files."""
    try:
        text = inp_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return bool(
        _ORCA_BANG_RE.search(text)
        or _ORCA_BLOCK_RE.search(text)
        or _ORCA_XYZFILE_RE.search(text)
    )


def required_orca_outputs(inp_path=None, out_path=None) -> List[Path]:
    """Return generated artifacts required for a calculation to count as complete.

    Generic ORCA jobs only require a normal `.out` terminator marker. Extra
    side products like `.gbw` are workflow-specific and should be enforced by
    callers via `required_outputs` when a downstream step truly depends on
    them. This keeps legacy completed outputs recalc-skippable.
    """
    inp = Path(inp_path) if inp_path is not None else None
    out = Path(out_path) if out_path is not None else None

    if inp is None and out is None:
        return []

    if inp is not None:
        special_outputs = _SPECIAL_REQUIRED_OUTPUTS.get(inp.name.lower())
        if special_outputs is not None:
            return [inp.parent / name for name in special_outputs]

    if inp is not None and not _looks_like_orca_job(inp):
        return []
    return []


#: How far back from the end of an output the optimiser's verdict is looked
#: for. ORCA prints it before the final energy evaluation, the property
#: section and any frequency job that follows, so it can sit well before the
#: end — measured at 1 MB in a real ESD run. Bounded because these outputs
#: reach hundreds of MB and this is consulted for every job on every recalc.
_OPT_VERDICT_TAIL_BYTES = 8 * 1024 * 1024

#: ORCA's own words when an optimisation runs out of cycles. Testing for this
#: rather than for the absence of "THE OPTIMIZATION HAS CONVERGED" keeps job
#: types that never print a convergence banner — single points, frequency-only
#: runs, compound inputs — out of it entirely.
_OPT_GAVE_UP_MARKER = "The optimization did not converge"


def optimization_gave_up(out_path) -> bool:
    """True when ORCA reported an optimisation that ran out of cycles.

    Such a run still ends with ``ORCA TERMINATED NORMALLY``, so on the marker
    alone it is indistinguishable from a finished one. Left at that, recalc
    skips it on every resubmission — the job can never finish — and dependent
    jobs go on to use a geometry that is not a stationary point.
    """
    try:
        p = Path(out_path)
        size = p.stat().st_size
        with p.open(encoding="utf-8", errors="replace") as f:
            if size > _OPT_VERDICT_TAIL_BYTES:
                f.seek(size - _OPT_VERDICT_TAIL_BYTES)
            return _OPT_GAVE_UP_MARKER in f.read()
    except Exception:
        return False


#: ORCA copies the input it was given into the head of its output, one line at
#: a time behind a `|  12> ` marker. That echo is the only self-contained
#: record of which input produced a given output.
_ECHOED_INPUT_RE = re.compile(r"^\|\s*\d+>\s?(.*)$")

#: The echo sits just after the banner. Reading the first megabyte finds it in
#: any real output without pulling a hundreds-of-MB file into memory.
_ECHO_SCAN_HEAD_BYTES = 1024 * 1024


def _normalise_input_lines(lines: Iterable[str]) -> List[str]:
    """Reduce input lines to what has to match for an output to be reusable.

    Blank lines and comments go, whitespace collapses, and anything that looks
    like a path keeps only its final component. That last part matters: ORCA is
    run in an isolated directory, and the copy it actually sees has its
    %moinp/xyzfile/hess references rewritten to absolute paths. Comparing those
    verbatim would report a difference on every single job.
    """
    out: List[str] = []
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        line = " ".join(line.split())
        line = re.sub(r'(?<=[\s"])/\S+/([^/\s"]+)', r"\1", line)
        out.append(line)
    return out


def echoed_input(out_path) -> Optional[List[str]]:
    """The input ORCA recorded in *out_path*, or None when it recorded none."""
    echoed: List[str] = []
    try:
        with Path(out_path).open(encoding="utf-8", errors="replace") as f:
            # Line by line, and stopping at ORCA's own end marker: the echo is
            # in the first pages, and nothing here should scale with the size
            # of an output that can run to hundreds of megabytes.
            budget = _ECHO_SCAN_HEAD_BYTES
            for line in f:
                budget -= len(line)
                if budget <= 0:
                    break
                m = _ECHOED_INPUT_RE.match(line.rstrip("\n"))
                if m is None:
                    continue
                body = m.group(1)
                if "END OF INPUT" in body:
                    break
                echoed.append(body)
    except OSError:
        return None
    return echoed or None


def output_matches_input(inp_path, out_path) -> Optional[bool]:
    """Whether *out_path* was produced by the input currently in *inp_path*.

    Returns None when the output carries no echo — an xTB run driven through a
    wrapper, a truncated file — because then there is nothing to compare and
    the caller has to fall back on its own judgement.

    This exists because a complete output is not evidence that it belongs to
    the input sitting next to it. Recalc used to treat it as such whenever the
    fingerprint sidecar was missing, and would then stamp the current input as
    the one that produced it. Change a functional, resubmit after a walltime
    kill, and the finished jobs are skipped and reported under settings they
    were never computed with.
    """
    echoed = echoed_input(out_path)
    if echoed is None:
        return None
    try:
        current = Path(inp_path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return None
    return _normalise_input_lines(echoed) == _normalise_input_lines(current)


def outputs_complete(inp_path, out_path, required_outputs: Optional[Iterable] = None) -> bool:
    """Return True when the main output and required generated files all exist."""
    if not has_ok_marker(out_path):
        return False

    if optimization_gave_up(out_path):
        return False

    if required_outputs is None:
        outputs = required_orca_outputs(inp_path=inp_path, out_path=out_path)
    else:
        outputs = [Path(path) for path in required_outputs]

    return all(path.exists() for path in outputs)


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

def should_skip(
    inp_path,
    out_path,
    extra_deps: Optional[Iterable] = None,
    required_outputs: Optional[Iterable] = None,
) -> bool:
    """Return True when ORCA execution can safely be skipped.

    All three conditions must hold:
    1. ``DELFIN_RECALC`` is enabled in the environment.
    2. *out_path* exists, contains ``ORCA TERMINATED NORMALLY``, and the
       required generated artifacts for dependent jobs are present.
    3. Smart mode only: the ``.fprint`` sidecar of *inp_path* matches the
       current fingerprint (i.e. neither inp content nor dependencies changed).
    """
    if not recalc_enabled():
        return False
    if not outputs_complete(inp_path, out_path, required_outputs=required_outputs):
        return False
    if not smart_mode_enabled():
        return True
    if not fingerprint_unchanged(inp_path, extra_deps=extra_deps):
        return False
    return True


def can_precomplete(
    inp_path,
    out_path,
    extra_deps: Optional[Iterable] = None,
    required_outputs: Optional[Iterable] = None,
) -> bool:
    """Return True when a recalc job can be marked complete before scheduling.

    In smart mode this mirrors ``should_skip()`` and also bootstraps a missing
    fingerprint once for legacy completed outputs.
    """
    if not recalc_enabled():
        return False

    inp = Path(inp_path) if inp_path is not None else None
    out = Path(out_path)

    if smart_mode_enabled():
        if inp is not None and should_skip(
            inp,
            out,
            extra_deps=extra_deps,
            required_outputs=required_outputs,
        ):
            return True

        if inp is None:
            return False

        if not outputs_complete(inp, out, required_outputs=required_outputs):
            return False

        fp_path = _fprint_path(inp)
        if not fp_path.exists():
            store_fingerprint(inp, extra_deps=extra_deps)
            return fp_path.exists()
        return fingerprint_unchanged(inp, extra_deps=extra_deps)

    return outputs_complete(inp, out, required_outputs=required_outputs)
