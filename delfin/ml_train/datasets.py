"""Dataset records, I/O and the licence / hold-out guards for delfin/ml_train.

A :class:`Record` describes one molecule (SMILES and/or an XYZ-coordinate
block) with named targets that carry units, plus per-record ``source``,
``licence`` and ``provenance``. Two hard guards are enforced here:

* **Source / licence guard** (``validate_record`` / :func:`guard_source`):
  data derived from the Cambridge Structural Database (``CSD`` / ``CCDC``)
  is refused, as is any licence outside the whitelist in
  :data:`ALLOWED_LICENCES`. COD (CC0) and our own calculations are fine.
* **Hold-out guard** (``guard_hold_out_path`` / :func:`guard_hold_out_text``):
  input files or recorded provenance that name one of the paper hold-out
  sets -- ``LOCKBOX2``, ``justus_batch/.../stage1_v2.txt``,
  ``.../heldout_stage1.txt``, ``pool_6000`` -- are refused by path
  component, so a training set can never silently train or tune on them.

Nothing here reads, opens or lists those hold-out files; the guard matches
against path components the caller passes in, and the tests use synthetic
paths only.

rdkit (a core DELFIN dependency) is imported lazily inside the functions
that need it, per DELFIN's convention; it is used for canonical InChIKey
computation in deduplication.
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.parse import quote, unquote

__all__ = [
    "ALLOWED_LICENCES",
    "HoldOutRefused",
    "LicenceRefused",
    "Record",
    "Target",
    "DedupReport",
    "guard_hold_out_path",
    "guard_hold_out_text",
    "guard_source",
    "validate_record",
    "deduplicate",
    "save_jsonl",
    "load_jsonl",
    "save_extxyz",
    "load_extxyz",
    "save",
    "load",
]


# ---------------------------------------------------------------------------
# Licence / source policy
# ---------------------------------------------------------------------------

#: Licences under which a dataset may be used. COD is CC0; our own
#: calculations are ``own-calculation`` / ``own``. Anything else is refused
#: as unknown (this is the hard rule: an unrecognised licence is not loaded).
ALLOWED_LICENCES = frozenset(
    {
        "CC0",
        "CC0-1.0",
        "CC-BY-4.0",
        "own-calculation",
        "own",
    }
)

#: Tokens that mark CSD/CCDC-derived data (refused whatever the licence).
_CSD_TOKENS = ("csd", "ccdc")


class LicenceRefused(ValueError):
    """Source is CSD/CCDC-derived, or the licence is unknown."""


class HoldOutRefused(ValueError):
    """Path or provenance names one of the paper hold-out sets."""


#: Path components that identify a paper hold-out set. Matching is by path
#: component (split on ``/`` and ``\\``), so a file whose *path* or
#: *provenance* contains any of these is refused. Keeping the tokens narrow
#: to the exact names from the rules avoids over-blocking ordinary data.
HOLD_OUT_COMPONENTS = frozenset(
    {
        "LOCKBOX2",
        "justus_batch",
        "competitors",
        "stages_v2",
        "stage1_v2.txt",
        "ANALYSE_2026_10_12_f3_haenger",
        "heldout_stage1.txt",
        "pool_6000",
    }
)


# ---------------------------------------------------------------------------
# Record types
# ---------------------------------------------------------------------------


@dataclass
class Target:
    """One named regression/classification target with its unit."""

    name: str
    value: float
    unit: Optional[str] = None

    def __post_init__(self) -> None:
        self.value = float(self.value)


@dataclass
class Record:
    """One dataset record: a molecule with targets and provenance.

    ``smiles`` and/or ``xyz`` must be present (at least one). ``xyz`` is an
    extended-XYZ coordinate block (count line, comment line, then one
    ``element x y z`` line per atom), consistent with :func:`save_extxyz`.
    ``provenance`` records where the molecule came from and is scanned by the
    hold-out guard; it may be a str or a dict of string values.
    """

    identifier: str
    source: str
    licence: str
    smiles: Optional[str] = None
    xyz: Optional[str] = None
    targets: List[Target] = field(default_factory=list)
    provenance: Optional[Any] = None

    def __post_init__(self) -> None:
        if not self.smiles and not self.xyz:
            raise ValueError(
                f"record {self.identifier!r} must carry a smiles or an xyz block"
            )


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------


def guard_source(source: str, licence: str, identifier: str = "?") -> None:
    """Refuse CSD/CCDC-derived sources and unknown licences."""
    src = (source or "").lower()
    if any(tok in src for tok in _CSD_TOKENS):
        raise LicenceRefused(
            f"record {identifier!r}: source {source!r} is CSD/CCDC-derived "
            "and is refused (no licence permits its use)"
        )
    if licence not in ALLOWED_LICENCES:
        raise LicenceRefused(
            f"record {identifier!r}: licence {licence!r} is not in the "
            f"allowed set {sorted(ALLOWED_LICENCES)}"
        )


def _components(value: Any) -> Iterable[str]:
    return (p for p in os.fspath(value).replace("\\", "/").split("/") if p)


def _refuses_hold_out(value: str) -> bool:
    return any(comp in HOLD_OUT_COMPONENTS for comp in _components(value))


def guard_hold_out_text(text: str, what: str = "value") -> None:
    """Refuse free text whose path components name a hold-out set."""
    if _refuses_hold_out(text):
        raise HoldOutRefused(f"{what} names a paper hold-out set: {text!r}")


def guard_hold_out_path(path: str) -> None:
    """Refuse an input *path* whose path components name a hold-out set.

    Called on every file before it is loaded into a trainable dataset; the
    guard only inspects the string components, it never opens the file.
    """
    if _refuses_hold_out(path):
        raise HoldOutRefused(f"path names a paper hold-out set: {path!r}")


def _provenance_text(provenance: Any) -> str:
    """Flatten a provenance dict/str/list into one component-matching string."""
    if provenance is None:
        return ""
    if isinstance(provenance, str):
        return provenance
    if isinstance(provenance, dict):
        return " ".join(str(v) for v in provenance.values())
    if isinstance(provenance, (list, tuple)):
        return " ".join(str(v) for v in provenance)
    return str(provenance)


def validate_record(record: Record) -> None:
    """Run every guard on one record; raise on the first failure."""
    guard_source(record.source, record.licence, record.identifier)
    prov = _provenance_text(record.provenance)
    if prov:
        guard_hold_out_text(prov, f"provenance of record {record.identifier!r}")
    # a record that already stores coordinates is also scanned, in case a
    # provenance-free xyz header quoted a hold-out set
    if record.xyz and _refuses_hold_out(record.xyz):
        raise HoldOutRefused(
            f"xyz header of record {record.identifier!r} names a hold-out set"
        )


# ---------------------------------------------------------------------------
# Serialisation helpers
# ---------------------------------------------------------------------------


def _record_to_dict(record: Record) -> Dict[str, Any]:
    dct = asdict(record)
    dct["targets"] = [
        {"name": t.name, "value": t.value, "unit": t.unit} for t in record.targets
    ]
    return dct


def _record_from_dict(dct: Dict[str, Any]) -> Record:
    dct = dict(dct)
    targets = dct.pop("targets", []) or []
    dct["targets"] = [
        Target(name=t["name"], value=t["value"], unit=t.get("unit"))
        for t in targets
    ]
    return Record(**dct)


# ---------------------------------------------------------------------------
# JSONL I/O
# ---------------------------------------------------------------------------


def save_jsonl(records: Iterable[Record], path: str) -> None:
    """Write records one JSON object per line. Validates each record first."""
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            validate_record(record)
            fh.write(json.dumps(_record_to_dict(record), sort_keys=True) + "\n")


def load_jsonl(path: str) -> List[Record]:
    """Load records from JSONL, running every guard on each."""
    guard_hold_out_path(path)
    records: List[Record] = []
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = _record_from_dict(json.loads(line))
            except (json.JSONDecodeError, TypeError, KeyError) as exc:
                raise ValueError(
                    f"invalid JSONL record at {path}:{lineno}: {exc}"
                ) from exc
            validate_record(record)
            records.append(record)
    return records


# ---------------------------------------------------------------------------
# extXYZ I/O
# ---------------------------------------------------------------------------

def _as_extxyz_block(record: Record) -> str:
    if not record.xyz:
        raise ValueError(
            f"record {record.identifier!r} has no xyz block; cannot write extxyz"
        )
    lines = record.xyz.splitlines()
    if not lines:
        raise ValueError(f"record {record.identifier!r} has empty xyz")
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(
            f"record {record.identifier!r} xyz has unreadable atom count"
        ) from exc
    coord = lines[2 : 2 + natoms]
    # metadata + targets as comment-line ``key=value`` tokens; values are
    # percent-encoded so that any whitespace saves verbatim instead of
    # splitting the token in two and silently truncating at the first space
    def _q(v: str) -> str:
        """Quote a metadata value; empty stays empty (unquoted)."""
        return quote(str(v), safe="") if v else ""

    parts = [
        f"identifier={_q(record.identifier)}",
        f"smiles={_q(record.smiles)}",
        f"source={_q(record.source)}",
        f"licence={_q(record.licence)}",
    ]
    for tgt in record.targets:
        parts.append(f"target.{_q(tgt.name)}={_encode_target(tgt)}")
    return f"{natoms}\n{' '.join(parts)}\n" + "\n".join(coord)


def save_extxyz(records: Iterable[Record], path: str) -> None:
    """Write records as one extended-XYZ frame each (multi-frame xyz)."""
    guard_hold_out_path(path)
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            validate_record(record)
            fh.write(_as_extxyz_block(record) + "\n")


def load_extxyz(path: str) -> List[Record]:
    """Read a multi-frame ext XYZ file written by :func:`save_extxyz`."""
    guard_hold_out_path(path)
    records: List[Record] = []
    with open(path, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    i = 0
    n = len(lines)
    while i < n:
        if not lines[i].strip():
            i += 1
            continue
        try:
            natoms = int(lines[i].strip())
        except ValueError as exc:
            raise ValueError(f"invalid extxyz natoms at frame starting line {i + 1}") from exc
        meta = lines[i + 1] if i + 1 < n else ""
        coord = lines[i + 2 : i + 2 + natoms]
        records.append(_record_from_extxyz(meta, coord))
        i += 2 + natoms
    return records


def _record_from_extxyz(meta: str, coord: List[str]) -> Record:
    kv: Dict[str, str] = {}
    targets: List[Target] = []
    for token in meta.split():
        k, sep, v = token.partition("=")
        if not sep:
            continue
        v = unquote(v)  # values were percent-encoded on write (spaces, '%')
        if k.startswith("target."):
            name = unquote(k[len("target."):])
            if not name:
                continue
            value, unit = _decode_target(v)
            if value is not None:
                targets.append(Target(name=name, value=value, unit=unit))
        else:
            kv[k] = v
    rec = Record(
        identifier=kv.get("identifier", "?"),
        source=kv.get("source", "?"),
        licence=kv.get("licence", "?"),
        smiles=kv.get("smiles") or None,
        xyz=f"{len(coord)}\n\n" + "\n".join(coord),
        targets=targets,
        provenance=None,
    )
    validate_record(rec)
    return rec


def _encode_target(t: Target) -> str:
    """``value`` in shortest exact round-trip form, then the unit, e.g. ``-6.4eV``.

    ``repr`` (not ``f"{value:g}"``) keeps every significant digit: ``:g`` at its
    default 6-significant-digit precision silently rounds, e.g. ``1234567.0 -> 1.23457e+06``,
    losing the last digits of an absolute energy. ``repr`` always round-trips via
    the reader, including its exponent form for very small/large magnitudes.
    """
    value = repr(t.value)
    return value if not t.unit else f"{value}{t.unit}"


def _decode_target(raw: str) -> Tuple[Optional[float], Optional[str]]:
    """Split ``-6.4eV`` into ``(-6.4, 'eV')``; ``(None, None)`` if not numeric.

    The value is written with ``repr(value)`` (:func:`_encode_target`), the
    shortest exact round-trip float representation, which uses exponent form
    for very small/large magnitudes (e.g. ``-5e-05Ha``, ``1.23457e+06eV``).
    The mantissa scan must therefore consume an exponent (``e``/``E`` with an
    optional sign and digits) as part of the number, or small absolute energies
    read back ~6 orders too large. A bare trailing ``e`` is not part of the
    number (a unit such as ``eV`` starts with ``e``).
    """
    i = 0
    n = len(raw)
    if i < n and raw[i] in "+-":
        i += 1
    saw_exp = False
    while i < n:
        c = raw[i]
        if c.isdigit() or c == ".":
            i += 1
        elif c in "eE" and not saw_exp:
            nxt = raw[i + 1] if i + 1 < n else ""
            if nxt in "+-" or nxt.isdigit() or nxt == ".":
                saw_exp = True
                i += 1
            else:
                break  # bare 'e': start of a unit like 'eV', not an exponent
        elif c in "+-" and saw_exp:
            i += 1
        else:
            break
    if i == 0:
        return None, None
    try:
        return float(raw[:i]), (raw[i:] or None)
    except ValueError:
        return None, None


# ---------------------------------------------------------------------------
# Generic save / load by extension
# ---------------------------------------------------------------------------


def save(records: Iterable[Record], path: str) -> None:
    """Dispatch on extension: ``.jsonl`` -> JSONL, ``.xyz``/``.extxyz`` -> ext XYZ."""
    armed = str(path).lower()
    if armed.endswith(".jsonl"):
        save_jsonl(records, path)
    elif armed.endswith((".xyz", ".extxyz")):
        save_extxyz(records, path)
    else:
        raise ValueError(
            f"cannot infer format from {path!r}; use a .jsonl or .xyz/.extxyz path"
        )


def load(path: str) -> List[Record]:
    """Load records from a file, dispatching on extension (see :func:`save`)."""
    armed = str(path).lower()
    if armed.endswith(".jsonl"):
        return load_jsonl(path)
    if armed.endswith((".xyz", ".extxyz")):
        return load_extxyz(path)
    raise ValueError(
        f"cannot infer format from {path!r}; use a .jsonl or .xyz/.extxyz path"
    )


# ---------------------------------------------------------------------------
# Deduplication by canonical InChIKey
# ---------------------------------------------------------------------------


@dataclass
class DedupEntry:
    """One dropped record and the kept record it collided with."""

    identifier: str
    duplicate_of: str
    inchikey: str


@dataclass
class DedupReport:
    """Summary of a :func:`deduplicate` pass."""

    kept_count: int
    dropped_count: int
    dropped: List[DedupEntry] = field(default_factory=list)

    @property
    def dropped_identifiers(self) -> List[str]:
        return [d.identifier for d in self.dropped]


def _inchikey_of(record: Record) -> str:
    """Canonical InChIKey of a record's molecule (lazy rdkit import)."""
    from rdkit import Chem

    if record.smiles:
        mol = Chem.MolFromSmiles(record.smiles)
        if mol is not None:
            inchi = Chem.MolToInchiKey(mol)
            if inchi:
                return inchi
    return record.identifier


def deduplicate(records: Iterable[Record]) -> Tuple[List[Record], DedupReport]:
    """Drop duplicate molecules by canonical InChIKey, keeping the first.

    Returns ``(kept, report)``. Two records with the same canonical InChIKey
    are considered the same molecule; the first occurrence is kept and later
    ones are reported with the kept record's identifier they collided with.
    Records without a parseable SMILES fall back to their own identifier as
    their key (never merged with another record).
    """
    kept: List[Record] = []
    seen: Dict[str, str] = {}  # inchikey -> kept identifier
    dropped: List[DedupEntry] = []
    for record in records:
        key = _inchikey_of(record)
        if key in seen:
            dropped.append(
                DedupEntry(
                    identifier=record.identifier,
                    duplicate_of=seen[key],
                    inchikey=key,
                )
            )
            continue
        seen[key] = record.identifier
        kept.append(record)
    report = DedupReport(
        kept_count=len(kept), dropped_count=len(dropped), dropped=dropped
    )
    return kept, report


def canonical_inchikey(smiles: str) -> Optional[str]:
    """Canonical InChIKey for a SMILES string, or None if unparseable (lazy rdkit)."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)
