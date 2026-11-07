"""TURBOMOLE version detection and control compatibility (package TM2, phase 1).

The installed TURBOMOLE version decides the generator and MCP target. The
``turbomole_mcp`` service speaks 8.0/8.1 only, while the cluster module is
7.9. This module never assumes a version -- it detects it from, in order:

1.  the ``TURBODIR`` environment variable (a path that usually embeds the
    version, e.g. ``/opt/turbomole/7.9/``);
2.  the resolved module name (e.g. ``chem/turbomole/7.9``);
3.  a binary's own version output text, passed in by the caller/tests.

When none of these yields a version the module returns ``None`` (unknown) --
it never falls back to a guessed default. An *unknown* installed version is
kept distinct from a *known* one: the caller may not refuse, but may not
claim compatibility either.

A generated ``control`` file is checked against the active version: a
``$directive`` that is not part of DELFIN's known TURBOMOLE control
vocabulary, or that requires a newer version than is installed, is REFUSED
with a clear message -- never run, never silently dropped. The generator
version is recorded so provenance says which rule set produced the file.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

__all__ = [
    "TurbomoleVersion",
    "Violation",
    "ControlCheckResult",
    "version_from_turbodir",
    "version_from_module_name",
    "version_from_binary_output",
    "detect_installed_version",
    "mcp_target_version",
    "check_control",
    "stamp_generator_version",
    "KNOWN_DIRECTIVES",
]


# --------------------------------------------------------------------------
# Version model
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class TurbomoleVersion:
    """A detected TURBOMOLE version (major.minor) and where it came from."""

    major: int
    minor: int
    source: str
    raw: str

    def __str__(self) -> str:
        return f"{self.major}.{self.minor}"

    def __lt__(self, other: "TurbomoleVersion") -> bool:
        return (self.major, self.minor) < (other.major, other.minor)

    def __le__(self, other: "TurbomoleVersion") -> bool:
        return (self.major, self.minor) <= (other.major, other.minor)

    @property
    def is_7x_line(self) -> bool:
        """True for the 7.x line (7.0 .. 7.9)."""
        return self.major == 7


# A regex that matches the version part of a path/module/output text. Anchored
# on digits separated by '.'/'_'/'-' so unrelated numbers (copyright years,
# etc.) are not grabbed; the caller narrows the string first (TURBODIR path,
# module name, binary version line). An optional single ``v``/``V`` marker is
# allowed directly before the major digit so version strings like
# "TURBOMOLE V8.1" or a path "turbomole-v8.1" parse -- but the marker letter
# itself must not be part of a word ("ver8.1" is not a version here).
_VERSION_RE = re.compile(r"(?<![A-Za-z0-9])[Vv]?(\d{1,2})[._-](\d{1,2})(?![A-Za-z0-9])")


def _parse(raw: str, source: str) -> Optional[TurbomoleVersion]:
    m = _VERSION_RE.search(raw)
    if not m:
        return None
    return TurbomoleVersion(int(m.group(1)), int(m.group(2)), source, raw.strip())


# --------------------------------------------------------------------------
# Detection -- never assume
# --------------------------------------------------------------------------

def version_from_turbodir(env: Optional[dict] = None) -> Optional[TurbomoleVersion]:
    """Read the version from ``TURBODIR``, e.g. ``.../turbomole/7.9/``.

    Absent or unparseable ``TURBODIR`` -> ``None`` (never assumed).
    """
    env = env if env is not None else os.environ
    turbodir = env.get("TURBODIR")
    if not turbodir:
        return None
    return _parse(turbodir, "TURBODIR")


def version_from_module_name(module_name: Optional[str]) -> Optional[TurbomoleVersion]:
    """Read the version from a resolved module name, e.g. ``chem/turbomole/7.9``.

    ``None`` when the name carries no parseable version (never assumed).
    """
    if not module_name:
        return None
    return _parse(module_name, "module")


def version_from_binary_output(text: Optional[str]) -> Optional[TurbomoleVersion]:
    """Read the version from a binary's version-output text.

    The caller supplies the tool's stdout/stderr (or a test fixture). An
    absent number -> ``None`` (never assumed).
    """
    if not text:
        return None
    return _parse(text, "binary")


def detect_installed_version(
    *,
    turbodir: Optional[str] = None,
    module_name: Optional[str] = None,
    binary_output: Optional[str] = None,
) -> Optional[TurbomoleVersion]:
    """Detect the installed version from the available sources, in order.

    Order of trust: ``TURBODIR``, then the resolved module name, then the
    binary's version output. If every source is absent or unparseable this
    returns ``None`` -- the caller must treat that as *unknown*, not as any
    particular version.
    """
    if turbodir:
        from_td = version_from_turbodir({"TURBODIR": turbodir})
        if from_td:
            return from_td
    if module_name:
        from_mod = version_from_module_name(module_name)
        if from_mod:
            return from_mod
    if binary_output:
        from_bin = version_from_binary_output(binary_output)
        if from_bin:
            return from_bin
    # Fall back to the environment only when no source was given explicitly.
    if turbodir is None and module_name is None and binary_output is None:
        from_env = version_from_turbodir()
        if from_env:
            return from_env
    return None


# --------------------------------------------------------------------------
# MCP target rule
# --------------------------------------------------------------------------

def mcp_target_version(
    installed: Optional[TurbomoleVersion],
) -> tuple[Optional[str], bool]:
    """The version to ask the MCP for, given the installed version.

    The MCP speaks 8.0/8.1 only. When the installed version is below the
    8.x line (7.x, or any earlier line) we still ask the MCP for the oldest
    version it speaks (8.0) and say so explicitly -- the second element of
    the tuple is True when the MCP target differs from the installed
    version. For an *unknown* installed version the target is the
    conservative 8.0 and ``differs`` is True (we
    cannot prove compatibility, so we flag it). For 8.x we ask for the exact
    installed version and ``differs`` is False.
    """
    if installed is None or installed.major < 8:
        return "8.0", True
    return str(installed), False


# --------------------------------------------------------------------------
# Control-directive compatibility check
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Violation:
    """One reason a generated control is refused for the installed version."""

    directive: str
    installed: Optional[str]
    message: str


@dataclass(frozen=True)
class ControlCheckResult:
    ok: bool
    violations: tuple[Violation, ...]

    def __bool__(self) -> bool:
        return self.ok


# The control ``$directives`` that DELFIN's own TURBOMOLE setup writes or
# knows. Everything DELFIN emits from ``_setup_control`` lives here. A
# directive that is NOT in this set is unknown to DELFIN => REFUSED with a
# clear message rather than silently passed to the engine. This set is the
# honest boundary: DELFIN only ever generates control files out of this
# vocabulary, and an MCP answer that proposes something outside it is
# rejected, not run.
KNOWN_DIRECTIVES = frozenset(
    {
        "title",
        "operating system",
        "symmetry",
        "coord",
        "redundant",
        "atoms",
        "basis",
        "basis file",
        "charge",
        "unpaired",
        "uhf",
        "dft",
        "ricore",
        "maxcor",
        "scfconv",
        "scfiterlimit",
        "scfiterlimit0",
        "ridft",
        "rij",
        "dscf",
        "cosmo",
        "cosmo_rs",
        "disp3",
        "disp",
        "disp4",
        "idft",
        "exx",
        "tddft",
        "grad",
        "energy",
        "hessian",
        "vibrational spectrum",
        "vibrational analysis",
        "last energy",
        "actual step",
        "lock",
        "data",
        "end",
    }
)

# Directives that require a newer TURBOMOLE line than 7.x. This gated map is
# deliberately small and conservative: it only lists directives with a
# documented line requirement, and it is the *mechanism* the compatibility
# check uses. Relaxing it for a specific version is a one-line change inside
# ``_allows``.
_GATED_TO_8X = frozenset({"cosmo_rs"})  # COSMO-RS solvent model: 8.x recommended


def _allows(directive: str, installed: Optional[TurbomoleVersion]) -> bool:
    """Does ``installed`` support ``directive``?

    Unknown installed version: no gate can be proved, so a directive is
    allowed -- but the caller must not claim full compatibility (the check
    still refuses directives unknown to DELFIN). A 7.x line cannot emit the
    Cosmo-rs block. Everything else in DELFIN's vocabulary is long
    established on 7.x as well.
    """
    if installed is None:
        return True
    if directive in _GATED_TO_8X:
        return not installed.is_7x_line
    return True


def _directives_of(
    control_text: str, known_directives: Sequence[str] = ()
) -> Iterable[str]:
    """Yield the ``$directive`` names present in a control file's text.

    A directive is a ``$name`` token at the start of a line (allowing the
    ``$name`` or ``$name value`` forms TURBOMOLE uses, and multi-word names
    such as ``$operating system``). Lines inside a block (``$atoms ...`` /
    ``basis=...``) have no ``$`` prefix and are ignored.

    Multi-word names are matched greedily against ``known_directives`` (the
    longest name that is a whitespace-/``=``-aligned prefix wins); anything not
    in that set still yields its first word so an *unknown* directive remains
    visible and can be refused by the caller.
    """
    known_desc = tuple(sorted(known_directives, key=len, reverse=True))
    for line in (control_text or "").splitlines():
        stripped = line.strip()
        if not stripped.startswith("$"):
            continue
        head = stripped[1:].strip()
        if not head:
            continue
        matched = None
        for name in known_desc:
            if head == name or head.startswith(name + " ") or head.startswith(name + "="):
                matched = name
                break
        if matched is not None:
            yield matched
            continue
        token = head.split(None, 1)[0]
        if "=" in token:
            token = token.split("=", 1)[0]
        if token:
            yield token.strip()


def check_control(
    control_text: str,
    installed: Optional[TurbomoleVersion] = None,
    *,
    known_directives: Optional[Sequence[str]] = None,
) -> ControlCheckResult:
    """Refuse a generated control when a directive is not supported.

    Two refusal reasons are produced, each with a clear message:

    * a ``$directive`` that is not in DELFIN's known control vocabulary
      (a typo, or an MCP answer proposing something DELFIN cannot map) is
      REFUSED -- it would be run silently otherwise;
    * a directive that the installed version does not support (version 7.x
      trying to emit an 8.x-gated block) is REFUSED.

    An unknown installed version never silently passes: directives unknown
    to DELFIN are still refused, and the caller is expected to keep it
    distinct from a verified ``ok``.
    """
    known = frozenset(known_directives) if known_directives is not None else KNOWN_DIRECTIVES
    known_seq = tuple(known)
    violations: list[Violation] = []
    for directive in _directives_of(control_text, known_seq):
        if directive not in known:
            violations.append(
                Violation(
                    directive,
                    str(installed) if installed else None,
                    f"refusing directive '${directive}': not a known TURBOMOLE "
                    f"control directive, so it cannot be checked against the "
                    f"installed version ({installed or 'unknown'})",
                )
            )
        elif not _allows(directive, installed):
            violations.append(
                Violation(
                    directive,
                    str(installed) if installed else None,
                    f"refusing directive '${directive}': requires the 8.x line, "
                    f"but the installed TURBOMOLE is {installed}",
                )
            )
    return ControlCheckResult(not violations, tuple(violations))


# --------------------------------------------------------------------------
# Provenance stamp
# --------------------------------------------------------------------------

def stamp_generator_version(control_text: str, generator: str, version: str) -> str:
    """Append a generator-version provenance note to a control file's text.

    Returns the text with a ``$title`` comment recording which generator and
    generator version produced it, so a later reader knows which rule set
    made the file. If the text already carries the same note it is returned
    unchanged (idempotent).
    """
    marker = f"# generated by {generator} {version}"
    if marker in control_text:
        return control_text
    return control_text.rstrip() + "\n" + marker + "\n"
