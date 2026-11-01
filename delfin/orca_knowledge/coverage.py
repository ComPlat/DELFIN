"""The ORCA coverage map: what DELFIN's workflows write, for what, and what OPI can build.

``coverage.json`` (data, one record per keyword / block option) tells a builder
what a keyword DELFIN already emits is *for*, which CONTROL key sets it, which
workflows use it, and whether the pinned OPI can express it.  It is the map the
OPI building blocks consult before they claim to reproduce a DELFIN input.

Every record::

    {"key": "RIJCOSX", "kind": "simple_keyword", "purpose": "...",
     "manual": "essentialelements/integralhandling.txt", "probe": "accepted",
     "control": "ri_jkx", "opi": "opi.input.simple_keywords.approximation",
     "workflows": ["classic", "esd"], "source": "delfin/xyz_io.py:527"}

Guarantees enforced by the tests:

* every record names a manual section and carries a probe verdict;
* every literal ``%block`` / bang-line keyword a DELFIN writer emits is in the
  map -- a new keyword added to a workflow cannot appear unexplained.
"""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

COVERAGE_FILE = Path(__file__).with_name("coverage.json")

_REQUIRED = frozenset({"key", "kind", "purpose", "manual", "probe", "control",
                       "opi", "workflows", "source"})
_KINDS = frozenset({"simple_keyword", "block_option"})
_KNOWN_PROBES = frozenset({"accepted", "rejected", "parsed", "finished",
                           "timeout", "refused", "unavailable"})


class CoverageError(ValueError):
    """coverage.json is missing or violates the schema."""


@lru_cache(maxsize=1)
def load_coverage() -> dict:
    """Load and validate ``coverage.json`` once per process."""
    if not COVERAGE_FILE.exists():
        raise CoverageError(f"missing {COVERAGE_FILE.name} next to {__file__}")
    with COVERAGE_FILE.open(encoding="utf-8") as fh:
        data = json.load(fh)
    records = data.get("records")
    if not isinstance(records, dict) or not records:
        raise CoverageError("coverage.json: 'records' must be a non-empty dict")
    for key, rec in records.items():
        if not isinstance(rec, dict):
            raise CoverageError(f"coverage.json: record {key!r} is not an object")
        missing = _REQUIRED - set(rec)
        if missing:
            raise CoverageError(f"coverage.json: record {key!r} lacks fields {sorted(missing)}")
        if rec.get("kind") not in _KINDS:
            raise CoverageError(f"coverage.json: record {key!r}: unknown kind {rec.get('kind')!r}")
        if not rec.get("manual"):
            raise CoverageError(f"coverage.json: record {key!r} names no manual section")
        if rec.get("probe") not in _KNOWN_PROBES:
            raise CoverageError(
                f"coverage.json: record {key!r}: unknown probe {rec.get('probe')!r}")
    return data


def record(keyword: str) -> dict | None:
    """The coverage record for a keyword / block option, or None."""
    return load_coverage()["records"].get(keyword)


def is_covered(keyword: str) -> bool:
    """Is a keyword / block option present in the map?"""
    return keyword in load_coverage()["records"]


def covers_block(block: str) -> bool:
    """Has the map at least one record whose key opens a given ``%block``?

    ``block`` must carry the ``%`` (e.g. ``"%scf"``).  True when any record key
    starts with ``%<block>`` (so ``"%scf maxiter"`` covers ``"%scf"``); False
    when the block is unknown to DELFIN's writers entirely.  Case-insensitive
    (``%QMMM`` matches the ``%QMMM`` record).
    """
    stripped = block.strip()
    if not stripped.startswith("%"):
        return False
    prefix = stripped[1:]  # the block stem, no %
    for k in load_coverage()["records"]:
        stem = k[1:].split(" ", 1)[0] if k.startswith("%") else ""
        if stem.lower() == prefix.lower():
            return True
    return False


def covered_writers() -> frozenset[str]:
    """Literal ``%block`` and bang-line keywords the writers are known to emit."""
    return frozenset(load_coverage()["records"])


#: Writer files whose literal keywords the map must cover.  Paths are relative to
#: the repo root and are the ones the 'every emitted keyword is mapped' test greps.
WRITER_FILES = [
    "delfin/common/orca_blocks.py",
    "delfin/common/solvation.py",
    "delfin/common/tddft_settings.py",
    "delfin/xyz_io.py",
    "delfin/esd_input_generator.py",
    "delfin/workflows/engine/classic.py",
    "delfin/workflows/engine/occupier.py",
    "delfin/build_up_complex.py",
    "delfin/calculators.py",
]

#: Literal block stems DELFIN writes in these files (the ``%`` prefix is dropped
#: because it also opens non-block keywords like ``%pal``).  A writer that emits a
#: ``%<name>`` block not listed here -- or not in the map -- fails the test.
_BLOCK_NAME = re.compile(r"%([A-Za-z][A-Za-z0-9_]*)")


def emitted_blocks(source: str) -> set[str]:
    """The literal ``%<word>`` tokens a writer emits, lowercased and deduplicated.

    Deny-list based: every ``%<word>`` not proven to be prose/format-spec is
    returned, so a *new* ``%block`` a developer adds to a writer is seen and
    must be covered.  Prose tokens (``%s``, ``%d``, ``%(name)s`` format specs,
    log words like ``%dh``) are excluded so the grep is not noisy.  A writer
    that emits a ``%block`` for which the map has no record fails the suite.
    """
    return {f"%{m.group(1).lower()}"
            for m in _BLOCK_NAME.finditer(source)
            if m.group(1).lower() not in _PROSE_TOKENS}


#: ``%<word>`` tokens that are Python/log format specifiers or prose, not ORCA
#: input DELFIN emits.  A token here is a false positive of the grep; anything
#: else the writers contain must be mapped.
_PROSE_TOKENS = frozenset({
    # Python %-format specifiers (writers format strings with them)
    "s", "d", "r", "f", "x", "o", "dh", "dh:2d", "s_step", "s_step_",
})
