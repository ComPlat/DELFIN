"""How ORCA is used correctly: knowledge cards, checked before they are used.

Input: cards in ``cards/*.json`` (each file a list of cards). Output: the
verified cards (:func:`cards`), the reasons a card is refused
(:func:`validate`), and whether an ORCA keyword exists (:func:`known_keyword`).

A card is one statement about ORCA, with the ORCA versions it holds for and
its sources::

    {"id": "scf.slowconv.open-shell-tm", "topic": "scf", "kind": "recommendation",
     "statement": "...", "orca": ">=6.0", "keywords": ["SLOWCONV"],
     "when": {...}, "recommend": {...}, "avoid": [...], "severity": "warning",
     "sources": [{"type": "manual", "version": "6.1",
                  "page": "contents/essentialelements/scf.html", "section": "..."},
                 {"type": "paper", "doi": "10.1063/5.0004608"}],
     "status": "verified", "author": "agent-07", "verified_by": "agent-11"}

Rules a card cannot pass:

* the schema: a known ``topic`` and ``kind``, a statement of 20 to 1500
  characters, an ``orca`` version range (``>=6.0``, ``>=6.0,<7``), an id
  of ``[a-z0-9._-]``;
* sources: at least one; a manual source names its version and page, a
  paper its DOI, an ORCA message its text, an archive finding its count;
* keywords: every keyword named must be one ORCA knows -- the catalog of
  the pinned OPI (``opi_keywords.json``, written by :func:`keywords_from_opi`),
  its solvents inside ``CPCM(...)``, ``SMD(...)``, ``ALPB(...)``,
  ``DDCOSMO(...)``, and a short list of aliases the manual uses, each with
  its page. A misspelt or invented keyword is refused;
* use: only ``verified`` cards are served, and a card is verified only by
  someone other than its author.

Rejected: free prose notes (nothing can check them), knowledge looked up
afresh at every run (not reproducible, and each look can misread the
manual differently), YAML (one more parser for data a JSON file holds).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .checks import applies, check_graph, findings  # noqa: F401  (re-exported)
from .coverage import (  # noqa: F401  (re-exported)
    CoverageError, covers_block, covered_writers, emitted_blocks,
    is_covered, load_coverage, record,
)

CARDS_DIR = Path(__file__).with_name("cards")
KEYWORDS_FILE = Path(__file__).with_name("opi_keywords.json")

TOPICS = frozenset({
    "general", "input", "parallel", "memory", "output", "scf", "guess", "stability",
    "basis", "ecp", "ri", "grid", "dispersion", "solvation", "relativistic", "symmetry",
    "hf", "dft", "double-hybrid", "composite", "semiempirical", "xtb", "mp2", "dlpno",
    "coupled-cluster", "casscf", "nevpt2", "caspt2", "mrci", "dmrg", "ice-ci", "autoci",
    "tddft", "eom", "roci", "excited-states", "esd", "properties", "nmr", "epr",
    "magnetism", "spectroscopy", "frequencies", "thermochemistry", "optimization",
    "transition-state", "irc", "neb", "scan", "conformers", "goat", "docker", "solvator",
    "md", "qmmm", "embedding", "broken-symmetry", "spin", "compound", "workflow",
    "error", "performance",
    # essential elements
    "counterpoise", "gcp", "fragments", "cpscf", "field", "frozen-core",
    "integrals", "coordinates", "numerical-gradients",
    # solvation models (Section 2.13, distinct models)
    "cpcm", "smd", "cosmors",
    # model chemistries
    "rpa", "adft", "casdft", "dcd-cas", "fci", "wavefunction-types", "mrcc",
    "mreom", "mm", "3c",
    # population and decomposition analyses (Section 5)
    "population", "nbo", "nocv", "eda", "led", "hfld", "adld", "fod",
    # vibrations and bandshapes
    "vibrations", "vpt2", "asa",
    # further property types
    "electric-properties", "deltascf", "xas", "xes", "mcd", "dichroism",
    "moessbauer", "spin-orbit", "aniso", "magprop",
    # excited-state methods not covered above
    "steom", "rocis", "mcrpa", "ops", "ihfsmr", "casscf-resp",
    # structure and reactivity beyond opt/TS/IRC/NEB
    "conical-intersections", "mecp", "dock",
    # multiscale
    "fmm", "crystal-qmmm",
    # utilities, workflows, environment
    "plots", "orca-2json", "property-file", "utilities", "opi",
    "installation", "gui", "troubleshooting",
    # ab initio molecular dynamics
    "aimd",
})
KINDS = frozenset({
    "requirement", "recommendation", "default", "incompatibility", "error-pattern",
    "keyword", "limit", "procedure", "check",
})
SOURCE_TYPES = frozenset({"manual", "paper", "orca-message", "orca-run", "archive", "delfin-doc"})
STATUSES = frozenset({"proposed", "verified", "rejected"})
SEVERITIES = frozenset({"error", "warning", "info"})

_ID = re.compile(r"^[a-z0-9][a-z0-9._-]{2,119}$")
_DOI = re.compile(r"^10\.\d{4,9}/\S+$")
_RANGE = re.compile(r"^(?:(?:>=|<=|>|<|==)\d+(?:\.\d+){0,2})(?:,(?:>=|<=|>|<|==)\d+(?:\.\d+){0,2})*$")
_PARAM = re.compile(r"^([a-z0-9-]+)\((.+)\)$")
_SOLVATION_FORMS = ("cpcm", "smd", "alpb", "ddcosmo", "cosmors", "cpcmc", "smd18")


def _catalog() -> dict:
    try:
        return json.loads(KEYWORDS_FILE.read_text())
    except (OSError, ValueError):
        return {"keywords": [], "solvents": [], "aliases": {}}


_CATALOG: dict | None = None


def _keywords() -> tuple[set[str], set[str], set[str], dict[str, set[str]]]:
    global _CATALOG
    if _CATALOG is None:
        _CATALOG = _catalog()
    # OPI's catalog, and the words ORCA's own parser accepted (the catalog
    # does not list every keyword ORCA knows: FIC-NEVPT2, RI-NEVPT2, ...).
    kws = {k.lower() for k in [*_CATALOG.get("keywords", []), *_CATALOG.get("orca_parser", [])]}
    solvents = {s.lower().replace(" ", "") for s in _CATALOG.get("solvents", [])}
    aliases = {k.lower() for k in (_CATALOG.get("aliases") or {})}
    # Block names -> option names. "blocks" is what OPI's typed models
    # know; "blocks_orca_parser" are "block.option" pairs ORCA's own parser
    # accepted (measured by probes, kept when the catalog is rewritten).
    blocks: dict[str, set[str]] = {}
    for block, options in (_CATALOG.get("blocks") or {}).items():
        blocks[str(block).lower()] = {str(o).lower() for o in options or []}
    for pair in _CATALOG.get("blocks_orca_parser") or []:
        block, _, option = str(pair).partition(".")
        if option:
            blocks.setdefault(block.lower(), set()).add(option.lower())
    return kws, solvents, aliases, blocks


def known_keyword(word: str) -> bool:
    """Whether ORCA knows ``word`` as a simple-input keyword (case-blind)."""
    w = str(word or "").strip().lower()
    if not w:
        return False
    kws, solvents, aliases, _blocks = _keywords()
    if w in kws or w in aliases:
        return True
    m = _PARAM.match(w)
    if m and m.group(1) in _SOLVATION_FORMS:
        return m.group(2).replace(" ", "") in solvents
    return False


def known_block_option(block: str, option: str = "") -> bool:
    """Whether ORCA knows ``option`` in ``block`` (case-blind).

    From the ``blocks`` section of the catalog (OPI's typed block models,
    written by ``refresh-keywords``) plus ``blocks_orca_parser``, the
    options a probe measured ORCA's own parser accepting. With an empty
    ``option`` it asks whether the ``block`` itself is known.
    """
    b, o = str(block or "").strip().lower(), str(option or "").strip().lower()
    if not b:
        return False
    known = _keywords()[3]
    if not o:
        return b in known
    return o in known.get(b, ())


def validate(card: Any) -> list[str]:
    """Every reason ``card`` is refused; [] when it passes."""
    if not isinstance(card, dict):
        return ["a card is a JSON object"]
    out: list[str] = []
    cid = str(card.get("id") or "")
    if not _ID.match(cid):
        out.append(f"id {cid!r}: lowercase letters, digits, '.', '_', '-' (3-120)")
    if card.get("topic") not in TOPICS:
        out.append(f"topic {card.get('topic')!r} is not one of the known topics")
    if card.get("kind") not in KINDS:
        out.append(f"kind {card.get('kind')!r} is not one of {sorted(KINDS)}")
    statement = card.get("statement")
    if not isinstance(statement, str) or not 20 <= len(statement.strip()) <= 1500:
        out.append("statement: 20 to 1500 characters")
    if not isinstance(card.get("orca"), str) or not _RANGE.match(card["orca"].replace(" ", "")):
        out.append("orca: a version range such as '>=6.0' or '>=6.0,<7'")
    if "severity" in card and card["severity"] not in SEVERITIES:
        out.append(f"severity is one of {sorted(SEVERITIES)}")
    for kw in _named_keywords(card):
        if not known_keyword(kw):
            out.append(f"keyword {kw!r} is not one ORCA knows (OPI catalog, ORCA's parser)")
    rec = card.get("recommend")
    if isinstance(rec, dict) and "blocks" in rec:
        out += blocks_problems(rec["blocks"], "recommend.blocks")
    avoid = card.get("avoid")
    if isinstance(avoid, dict):
        if "blocks" in avoid:
            out += blocks_problems(avoid["blocks"], "avoid.blocks")
        if set(avoid) - {"blocks"}:
            out.append("avoid as an object carries only 'blocks' "
                       "(advice text stays a list of strings)")
    sources = card.get("sources")
    if not isinstance(sources, list) or not sources:
        out.append("sources: at least one source")
    else:
        for i, src in enumerate(sources):
            out += [f"source {i}: {p}" for p in _source_problems(src)]
    status = card.get("status")
    if status not in STATUSES:
        out.append(f"status is one of {sorted(STATUSES)}")
    if not str(card.get("author") or "").strip():
        out.append("author: who wrote the card")
    if status == "verified":
        verifier = str(card.get("verified_by") or "").strip()
        if not verifier or verifier == str(card.get("author") or "").strip():
            out.append("verified_by: a verifier other than the author")
    out += _v2_problems(card)
    return out


def _v2_problems(card: dict) -> list[str]:
    """v2 fields (all optional): ``applies``, ``check``, ``relations``."""
    from .checks import (check_problems, condition_problems,
                         relations_problems)
    out: list[str] = []
    if "applies" in card:
        out += condition_problems(card["applies"], known_keyword,
                                  known_block_option=known_block_option)
    if "check" in card:
        out += check_problems(card["check"], known_keyword,
                              known_block_option=known_block_option)
    if "relations" in card:
        out += relations_problems(card["relations"])
    return out


# A keyword marker is "!" at the start of the text or after whitespace or
# "(" -- never after a letter (so "WARNING!" is no marker). Quoted text is
# checked like any other: an input line quoted in a card is advice too.
# The keyword itself stops at whitespace, trailing punctuation (".,;:")
# and an unbalanced closing ")" (so "(see !TightSCF)" yields TightSCF,
# not "TightSCF)").
_BANG = re.compile(r"(?:^|(?<=[\s(\"']))!\s*([A-Za-z0-9][A-Za-z0-9+\-*/(),_.']*)")
_TRAILING = re.compile(r"[.,;:\"']+$")


def blocks_problems(blocks: Any, field: str) -> list[str]:
    """Every reason a structured ``blocks`` advice is not one ORCA reads:
    block and option names must be known (catalog), the mapping must be
    ``{block: {option: value}}``. Called for ``recommend.blocks`` and
    ``avoid.blocks`` -- the audits found prose options and options put
    into the wrong block; both are refused here, with the name.
    """
    if not isinstance(blocks, dict) or not blocks:
        return [f"{field}: {blocks!r} is not {{block: {{option: value}}}}"]
    out: list[str] = []
    for block, options in blocks.items():
        name = str(block or "").strip().lower()
        if not name or not isinstance(options, dict) or not options:
            out.append(f"{field}.{block!r}: {{option: value}} pairs")
            continue
        if not known_block_option(name, "") and name not in _keywords()[3]:
            out.append(f"{field}: block {block!r} is not one ORCA knows "
                       "(OPI catalog, ORCA's parser)")
        for option, value in options.items():
            if not str(option or "").strip():
                out.append(f"{field}.{block!r}: an option name")
            elif not known_block_option(name, str(option)):
                out.append(f"{field}: option {option!r} is not one block {block!r} "
                           "knows (OPI catalog, ORCA's parser)")
            if not isinstance(value, (str, int, float, bool)) or isinstance(value, dict):
                out.append(f"{field}.{block}.{option}: {value!r} is not an input value")
    return out


def to_input(blocks: dict) -> str:
    """Render structured blocks as complete ORCA block input::

    >>> to_input({"mdci": {"TCutPNO": "3.33e-7"}})
    '%mdci\n  TCutPNO 3.33e-7\nend'

    Multiple blocks in insertion order, two spaces before each option.
    """
    out: list[str] = []
    for block, options in (blocks or {}).items():
        if not str(block).strip():
            continue
        out.append(f"%{str(block).strip().lower()}")
        for option, value in (options or {}).items():
            out.append(f"  {str(option).strip()} {value}")
        out.append("end")
    return "\n".join(out)


def _named_keywords(card: dict) -> list[str]:
    """Every keyword a card names: its ``keywords``, those of its advice
    (``recommend.keywords``, ``avoid``) and every ``!Keyword`` in its text --
    a keyword advised only in prose is checked like a listed one."""
    names = [str(k) for k in card.get("keywords") or []]
    rec = card.get("recommend")
    if isinstance(rec, dict):
        names += [str(k) for k in rec.get("keywords") or []]
    texts = [str(card.get("statement") or "")]
    for item in card.get("avoid") or []:
        texts.append(str(item))
    if isinstance(rec, dict):
        texts += [str(v) for v in rec.values() if isinstance(v, str)]
    for text in texts:
        for m in _BANG.finditer(text):
            token = m.group(1)
            # Trailing punctuation and an unbalanced closing ")" are not part
            # of the keyword: "(see !TightSCF)." names TightSCF.
            while True:
                stripped = _TRAILING.sub("", token)
                if stripped.endswith(")") and stripped.count(")") > stripped.count("("):
                    stripped = stripped[:-1]
                if stripped == token:
                    break
                token = stripped
            if token:
                names.append(token)
    return list(dict.fromkeys(n for n in names if n))


def _source_problems(src: Any) -> list[str]:
    if not isinstance(src, dict):
        return ["a source is an object"]
    kind = src.get("type")
    if kind not in SOURCE_TYPES:
        return [f"type is one of {sorted(SOURCE_TYPES)}"]
    if kind == "manual":
        out = []
        if not str(src.get("version") or "").strip():
            out.append("a manual source names its version")
        page = str(src.get("page") or "")
        if not page.endswith(".html") or ".." in page or page.startswith("/"):
            out.append("a manual source names its page (contents/....html)")
        return out
    if kind == "paper":
        return [] if _DOI.match(str(src.get("doi") or "")) else ["a paper source names its doi (10.xxxx/...)"]
    if kind == "orca-run":
        out = [f"an ORCA run source names its {f}" for f in ("version", "input", "output")
               if not str(src.get(f) or "").strip()]
        return out
    if kind == "orca-message":
        return [] if len(str(src.get("text") or "")) >= 8 else ["an ORCA message source quotes the message"]
    if kind == "archive":
        count = src.get("count")
        return [] if isinstance(count, int) and count > 0 else ["an archive source gives the count of runs"]
    if kind == "delfin-doc":
        path = str(src.get("path") or "")
        if not path:
            return ["a delfin-doc source names its path (docs/... or delfin/...)"]
        if ".." in path or path.startswith("/"):
            return ["a delfin-doc source path stays inside the repository (no '..', not absolute)"]
        if not path.startswith(("docs/", "delfin/")):
            return ["a delfin-doc source path is under docs/ or delfin/"]
        # Existence is checked by a test in the repository, not here: an
        # installed DELFIN has no docs/, and the working folder of a session
        # is not the repository.
        return []


def _files() -> list[Path]:
    return sorted(CARDS_DIR.glob("*.json")) if CARDS_DIR.is_dir() else []


def _read(path: Path) -> list:
    data = json.loads(path.read_text(encoding="utf-8"))
    return data if isinstance(data, list) else [data]


def cards(*, include_unverified: bool = False, topic: str | None = None) -> list[dict]:
    """The cards to use: verified ones (all with ``include_unverified``),
    each passing :func:`validate`."""
    out = []
    for path in _files():
        try:
            rows = _read(path)
        except (OSError, ValueError):
            continue
        for card in rows:
            if validate(card):
                continue
            if card["status"] == "rejected":
                continue
            if card["status"] != "verified" and not include_unverified:
                continue
            if topic and card["topic"] != topic:
                continue
            out.append(card)
    return out


def check_store() -> list[str]:
    """Every problem of the shipped store: unreadable files, invalid cards,
    duplicate ids."""
    problems: list[str] = []
    seen: dict[str, str] = {}
    for path in _files():
        try:
            rows = _read(path)
        except (OSError, ValueError) as exc:
            problems.append(f"{path.name}: not readable JSON ({exc})")
            continue
        for card in rows:
            cid = str(card.get("id") if isinstance(card, dict) else "")
            for p in validate(card):
                problems.append(f"{path.name} {cid}: {p}")
            if cid in seen:
                problems.append(f"{path.name} {cid}: duplicate id (also in {seen[cid]})")
            seen[cid] = path.name
    return problems


def keywords_from_opi() -> list[str]:
    """Every simple-input keyword the installed OPI knows, sorted."""
    import inspect
    import opi.input.simple_keywords as sk
    from opi.input.simple_keywords import SimpleKeyword
    found = set()
    for name in dir(sk):
        obj = getattr(sk, name)
        if inspect.isclass(obj):
            for attr in dir(obj):
                value = getattr(obj, attr)
                if isinstance(value, SimpleKeyword):
                    found.add(str(value))
    return sorted(found)


def solvents_from_opi() -> list[str]:
    from opi.input.simple_keywords import Solvent
    return sorted({str(getattr(Solvent, a)) for a in dir(Solvent) if a.isupper()})


def blocks_from_opi() -> dict[str, list[str]]:
    """Every input block OPI's typed models know, with its option names.

    Read from the pydantic fields of ``opi.input.blocks.Block*`` (the class
    attribute ``_name`` is the block name); field aliases, where OPI sets
    one, are included so a card may use either spelling.
    """
    import inspect
    import opi.input.blocks as blocks_mod
    from opi.input.blocks import Block
    out: dict[str, set[str]] = {}
    for name, obj in vars(blocks_mod).items():
        if not inspect.isclass(obj) or not issubclass(obj, Block) or obj is Block:
            continue
        pa = obj.__private_attributes__.get("_name")
        block = str(pa.default) if pa is not None else ""
        if not block:
            continue
        options = out.setdefault(block.lower(), set())
        for field in obj.model_fields.values():
            options.add(str(field.alias or field.validation_alias or "").lower())
        for field in obj.model_fields:
            options.add(field.lower())
    for options in out.values():
        options.discard("")
        options.discard("aftercoord")
    return {b: sorted(o) for b, o in sorted(out.items()) if o}


__all__ = ["validate", "cards", "check_store", "known_keyword", "known_block_option",
           "keywords_from_opi", "solvents_from_opi", "blocks_from_opi",
           "blocks_problems", "to_input",
           "TOPICS", "KINDS", "CARDS_DIR", "KEYWORDS_FILE",
           "applies", "findings", "check_graph"]
