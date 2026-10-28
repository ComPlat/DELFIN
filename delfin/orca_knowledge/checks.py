"""Executable knowledge cards (v2): conditions, checks and the card graph.

A v1 card is text a human reads. A v2 card may additionally carry:

``"applies"``
    A structured condition over a *calculation description* dict::

        {"keywords": ["RKS", "TIGHTSCF"],           # simple-input keywords, upper-case
         "blocks": {"scf": {"maxiter": 200}},       # block options actually set
         "method_family": "dft",                    # dft|hf|mp2|dlpno|casscf|...
         "job": "opt",                              # sp|opt|freq|ts|irc|neb|md|...
         "charge": 0, "mult": 1,
         "elements": ["Fe"],                        # element symbols
         "orca_version": "6.1.1"}

    Condition language (JSON, recursive)::

        {"all": [cond, ...]}          every sub-condition holds
        {"any": [cond, ...]}          at least one holds
        {"not": cond}                 the sub-condition does not hold
        {"keyword": "TightSCF"}       the simple-input keyword is set
        {"block_option": ["scf", "maxiter"]}
                                      the option is set in the block
        {"method_family": ["dft", "hf"]}
                                      the family is one of these
        {"job": ["opt", "ts"]}        the job is one of these
        {"mult": ">1"}                multiplicity compared with > >= < <= == (also "charge")
        {"elements_in": ["3d", "Fe"]} an element of the calc is in the set
                                      (period families 3d/4d/5d/4f/5f or symbols)
        {"orca": ">=6.0"}             the ORCA version satisfies the range

``"check"``
    What must hold when ``applies`` is true::

        {"requires_keyword": ["K", ...]}
        {"forbids_keyword": ["K", ...]}
        {"requires_block_option": [["casscf", "nel"], ...]}   # or "casscf.nel"
        {"forbids_combination": ["K", "scf.maxiter", ...]}    # not all at once

    Evaluated to findings with the card's severity and statement.

``"relations"``
    Links to other cards::

        {"requires": [id, ...], "conflicts_with": [id, ...],
         "supersedes": [id, ...], "see_also": [id, ...]}

All three fields are optional: every v1 card is a valid v2 card. Invalid
conditions are refused by :func:`delfin.orca_knowledge.validate` at write
time, so evaluation can be strict: :func:`applies` and the check evaluation
raise :class:`ValueError` on a malformed node instead of guessing.
"""

from __future__ import annotations

import re
from typing import Any, Iterable

#: Numeric comparison operators accepted by ``mult``/``charge`` conditions.
_OPS = (">=", "<=", "==", ">", "<")

_BLOCK_OPTION = re.compile(r"^[a-z][a-z0-9_-]*$")
_ID_SHAPE = re.compile(r"^[a-z0-9][a-z0-9._-]{2,119}$")

# Element symbol -> d/f family of the periodic table (period families only;
# main-group elements are named by symbol in conditions, not by family).
_ELEMENT_FAMILIES: dict[str, str] = {}
for _symbols, _fam in (
    (("Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"), "3d"),
    (("Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"), "4d"),
    (("Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"), "5d"),
    (("Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
      "Er", "Tm", "Yb"), "4f"),
    (("Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
      "Fm", "Md", "No", "Lr"), "5f"),
):
    for _s in _symbols:
        _ELEMENT_FAMILIES[_s.lower()] = _fam

CHECK_KEYS = ("requires_keyword", "forbids_keyword", "requires_block_option",
              "forbids_combination")
RELATION_KEYS = ("requires", "conflicts_with", "supersedes", "see_also")


# ---------------------------------------------------------------- conditions

def _version_tuple(version: str) -> tuple[int, ...]:
    parts = []
    for piece in str(version).strip().split("."):
        m = re.match(r"\d+", piece)
        if not m:
            raise ValueError(f"not a version number: {version!r}")
        parts.append(int(m.group(0)))
    return tuple(parts or [0])


def _compare_version(left: str, op: str, right: str) -> bool:
    a, b = _version_tuple(left), _version_tuple(right)
    n = max(len(a), len(b))
    a += (0,) * (n - len(a))
    b += (0,) * (n - len(b))
    return {"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b, "==": a == b}[op]


def _compare_int(left: Any, op: str, right: int) -> bool:
    if not isinstance(left, int):
        return False
    return {"<": left < right, "<=": left <= right, ">": left > right,
            ">=": left >= right, "==": left == right}[op]


def _parse_comparison(value: Any, field: str) -> tuple[str, int]:
    m = re.match(r"^(>=|<=|==|>|<)\s*(-?\d+)$", str(value).strip())
    if not m:
        raise ValueError(f"{field}: comparison must be like '>1' "
                         f"(one of {', '.join(_OPS)}, an integer)")
    return m.group(1), int(m.group(2))


def _family_of(symbol: str) -> str | None:
    return _ELEMENT_FAMILIES.get(str(symbol).strip().lower())


def _elements_match(wanted: list, present: list) -> bool:
    for w in wanted:
        w = str(w).strip()
        if w.lower() in ("3d", "4d", "5d", "4f", "5f"):
            if any(_family_of(e) == w.lower() for e in present):
                return True
        elif any(str(e).strip().lower() == w.lower() for e in present):
            return True
    return False


def condition_problems(cond: Any, known_keyword=None, _path: str = "applies") -> list[str]:
    """Every reason ``cond`` is not a well-formed condition; [] when it is.

    ``known_keyword`` (the :func:`delfin.orca_knowledge.known_keyword`
    callable) makes keywords inside conditions checked like the card's own.
    """
    if not isinstance(cond, dict):
        return [f"{_path}: a condition is a JSON object"]
    if len(cond) != 1:
        return [f"{_path}: a condition has exactly one key"]
    key, value = next(iter(cond.items()))
    if key in ("all", "any"):
        if not isinstance(value, list) or not value:
            return [f"{_path}.{key}: a non-empty list of conditions"]
        out = []
        for i, sub in enumerate(value):
            out += condition_problems(sub, known_keyword, f"{_path}.{key}[{i}]")
        return out
    if key == "not":
        return condition_problems(value, known_keyword, f"{_path}.not")
    if key == "keyword":
        if not isinstance(value, str) or not value.strip():
            return [f"{_path}.keyword: a keyword name"]
        if known_keyword and not known_keyword(value):
            return [f"{_path}.keyword {value!r} is not one ORCA knows "
                    "(OPI catalog, ORCA's parser)"]
        return []
    if key == "block_option":
        if (not isinstance(value, list) or len(value) != 2
                or not all(isinstance(p, str) and _BLOCK_OPTION.match(p) for p in value)):
            return [f"{_path}.block_option: [block, option], lowercase names"]
        return []
    if key in ("method_family", "job"):
        if not isinstance(value, list) or not value or not all(
                isinstance(v, str) and v.strip() for v in value):
            return [f"{_path}.{key}: a non-empty list of names"]
        return []
    if key in ("mult", "charge"):
        try:
            _parse_comparison(value, f"{_path}.{key}")
        except ValueError as exc:
            return [str(exc)]
        return []
    if key == "elements_in":
        if not isinstance(value, list) or not value or not all(
                isinstance(v, str) and v.strip() for v in value):
            return [f"{_path}.elements_in: a non-empty list of symbols or families"]
        bad = [v for v in value
               if str(v).strip().lower() not in ("3d", "4d", "5d", "4f", "5f")
               and len(str(v).strip()) not in (1, 2, 3)]
        if bad:
            return [f"{_path}.elements_in: {bad[0]!r} is not an element symbol "
                    "or 3d/4d/5d/4f/5f family"]
        return []
    if key == "orca":
        if not isinstance(value, str):
            return [f"{_path}.orca: a version range such as '>=6.0'"]
        expr = value.replace(" ", "")
        parts = expr.split(",")
        if not parts or not all(re.match(r"^(?:>=|<=|>|<|==)\d+(?:\.\d+){0,2}$", p)
                                for p in parts):
            return [f"{_path}.orca: a version range such as '>=6.0' or '>=6.0,<7'"]
        return []
    return [f"{_path}: unknown condition key {key!r}"]


def _eval(cond: Any, calc: dict, _path: str = "applies") -> bool:
    """Evaluate one condition node against ``calc``. Raises ValueError on a
    malformed node (validate refuses those before a card is shipped)."""
    if not isinstance(cond, dict) or len(cond) != 1:
        raise ValueError(f"{_path}: a condition is an object with one key")
    key, value = next(iter(cond.items()))
    if key == "all":
        return all(_eval(sub, calc, f"{_path}.all[{i}]") for i, sub in enumerate(value))
    if key == "any":
        return any(_eval(sub, calc, f"{_path}.any[{i}]") for i, sub in enumerate(value))
    if key == "not":
        return not _eval(value, calc, f"{_path}.not")
    if key == "keyword":
        return str(value).strip().lower() in {str(k).lower() for k in calc.get("keywords") or []}
    if key == "block_option":
        block, option = value
        block_data = calc.get("blocks") or {}
        return option in (block_data.get(block) or {})
    if key == "method_family":
        return str(calc.get("method_family") or "").lower() in {str(v).lower() for v in value}
    if key == "job":
        return str(calc.get("job") or "").lower() in {str(v).lower() for v in value}
    if key in ("mult", "charge"):
        op, number = _parse_comparison(value, f"{_path}.{key}")
        return _compare_int(calc.get(key), op, number)
    if key == "elements_in":
        return _elements_match(value, calc.get("elements") or [])
    if key == "orca":
        version = calc.get("orca_version")
        if version is None:
            # No version stated: the calc could be any version. Conservative:
            # true only when the range is not an upper bound.
            return _versionless_ok(value)
        return all(_compare_version(str(version), op, rhs)
                   for op, rhs in (re.match(r"^(>=|<=|==|>|<)(.+)$", p).groups()
                                   for p in value.replace(" ", "").split(",")))
    raise ValueError(f"{_path}: unknown condition key {key!r}")


def _versionless_ok(range_expr: str) -> bool:
    """True when the range holds for *some* version we cannot pin down -- a
    calc without a version is treated as "unknown", and a condition is true
    only if it cannot be excluded for any version (i.e. its upper bounds are
    open). Conservative: '>=6.0' -> True, '>=6.0,<7' -> True (most versions
    fit), '<6.0' -> False."""
    for op, rhs in (re.match(r"^(>=|<=|==|>|<)(.+)$", p).groups()
                    for p in range_expr.replace(" ", "").split(",")):
        if op in ("<", "<=", "=="):
            return False
    return True


def applies(card: dict, calc: dict) -> bool:
    """Whether the card's ``applies`` condition holds for ``calc``.

    A card without ``applies`` applies to every calculation (v1 behaviour).
    """
    cond = card.get("applies")
    if cond is None:
        return True
    return _eval(cond, calc)


# -------------------------------------------------------------------- checks

def check_problems(check: Any, known_keyword=None, _path: str = "check") -> list[str]:
    """Every reason ``check`` is not a well-formed check spec."""
    if not isinstance(check, dict) or not check:
        return [f"{_path}: a non-empty object"]
    out = []
    for key, value in check.items():
        if key not in CHECK_KEYS:
            out.append(f"{_path}: unknown check key {key!r} (one of {list(CHECK_KEYS)})")
            continue
        if key == "requires_block_option":
            if not _pairs(value):
                out.append(f"{_path}.{key}: [block, option], [[block, option], ...] "
                           "or \"block.option\"")
            continue
        if not isinstance(value, list) or not value or not all(
                isinstance(v, str) and v.strip() for v in value):
            out.append(f"{_path}.{key}: a non-empty list of keyword names")
            continue
        if key == "forbids_combination" and len(value) < 2:
            out.append(f"{_path}.{key}: at least two entries (a combination)")
        for i, entry in enumerate(value):
            if "." in entry and len(entry.split(".")) == 2:
                if not all(_BLOCK_OPTION.match(p) for p in entry.split(".")):
                    out.append(f"{_path}.{key}[{i}]: {entry!r} is not \"block.option\"")
                continue
            if known_keyword and not known_keyword(entry):
                out.append(f"{_path}.{key}[{i}] {entry!r} is not one ORCA knows "
                           "(OPI catalog, ORCA's parser)")
    return out


def _pairs(entries: Any) -> list[list[str]] | None:
    """``requires_block_option`` entries as (block, option) pairs.

    Both shapes are accepted: one flat pair ``["casscf", "nel"]``, a list of
    pairs ``[["casscf", "nel"], ["scf", "tol"]]``, or "block.option" strings.
    Returns None when the shape cannot be read.
    """
    if not isinstance(entries, list) or not entries:
        return None
    if len(entries) == 2 and all(isinstance(e, str) and _BLOCK_OPTION.match(e)
                                 for e in entries):
        return [list(entries)]
    out = []
    for entry in entries:
        if isinstance(entry, str) and "." in entry:
            parts = entry.split(".")
        elif isinstance(entry, list) and len(entry) == 2:
            parts = entry
        else:
            return None
        if len(parts) != 2 or not all(isinstance(p, str) and _BLOCK_OPTION.match(p)
                                      for p in parts):
            return None
        out.append(parts)
    return out


def evaluate_check(card: dict, calc: dict) -> list[dict]:
    """The findings of one card's ``check`` against ``calc``.

    Empty when the card does not apply, has no check, or the check holds.
    Each finding carries the card's id, severity, statement and sources plus
    a ``detail`` naming what failed. Malformed conditions raise ValueError.
    """
    if not applies(card, calc):
        return []
    check = card.get("check")
    if not check:
        return []
    cid = str(card.get("id") or "")
    severity = card.get("severity") if card.get("severity") in ("error", "warning", "info") \
        else "warning"
    base = {"id": cid, "severity": severity,
            "statement": str(card.get("statement") or ""),
            "sources": card.get("sources") or []}
    keywords = {str(k).lower() for k in calc.get("keywords") or []}
    blocks = calc.get("blocks") or {}
    out = []

    def finding(detail: str) -> dict:
        out.append({**base, "detail": detail})

    for kw in check.get("requires_keyword") or []:
        if kw.lower() not in keywords:
            finding(f"required keyword {kw} is not set")
    for kw in check.get("forbids_keyword") or []:
        if kw.lower() in keywords:
            finding(f"forbidden keyword {kw} is set")
    pairs = _pairs(check.get("requires_block_option") or []) or []
    for block, option in pairs:
        if option not in (blocks.get(block) or {}):
            finding(f"required block option {block}.{option} is not set")
    entries = check.get("forbids_combination") or []
    if entries:
        present = [e for e in entries if _entry_set(e, keywords, blocks)]
        if len(present) == len(entries):
            names = [e if isinstance(e, str) else ".".join(e) for e in entries]
            finding(f"forbidden combination is set: {', '.join(names)}")
    return out


def _entry_set(entry: Any, keywords: set[str], blocks: dict) -> bool:
    """Whether one ``forbids_combination`` entry is set in the calc: a
    keyword name, or a "block.option" / [block, option] reference."""
    if isinstance(entry, str) and "." not in entry:
        return entry.lower() in keywords
    parts = entry.split(".") if isinstance(entry, str) else list(entry)
    return parts[1] in (blocks.get(parts[0]) or {})


def findings(calc: dict, cards: Iterable[dict] | None = None) -> list[dict]:
    """Every finding of every applicable card against ``calc``.

    ``cards`` defaults to the shipped store (verified cards only). A card
    whose condition or check is malformed yields one finding naming the
    problem instead of being skipped silently.
    """
    if cards is None:
        from . import cards as load_cards
        cards = load_cards()
    out = []
    for card in cards:
        try:
            out += evaluate_check(card, calc)
        except ValueError as exc:
            out.append({"id": str(card.get("id") or ""),
                        "severity": "error", "statement": str(card.get("statement") or ""),
                        "sources": card.get("sources") or [],
                        "detail": f"card is malformed and cannot be evaluated: {exc}"})
    return out


# ------------------------------------------------------------------ relations

def relations_problems(relations: Any, _path: str = "relations") -> list[str]:
    """Every reason ``relations`` is not a well-formed relations object."""
    if not isinstance(relations, dict) or not relations:
        return [f"{_path}: a non-empty object"]
    out = []
    for key, value in relations.items():
        if key not in RELATION_KEYS:
            out.append(f"{_path}: unknown relation key {key!r} (one of {list(RELATION_KEYS)})")
        elif not isinstance(value, list) or not value or not all(
                isinstance(v, str) and _ID_SHAPE.match(v.strip()) for v in value):
            out.append(f"{_path}.{key}: a non-empty list of card ids")
    return out


def _flatten_conditions(cond: Any, out: dict, negated_keywords: set) -> None:
    """Collect the atomic constraints of a condition into ``out`` for the
    overlap test: method_family set, job set, mult/charge ranges, orca
    bounds, required keywords; keywords under ``not`` into
    ``negated_keywords``."""
    if not isinstance(cond, dict) or len(cond) != 1:
        return
    key, value = next(iter(cond.items()))
    if key in ("all", "any"):
        for sub in value or []:
            _flatten_conditions(sub, out, negated_keywords)
        return
    if key == "not":
        sub = value if isinstance(value, dict) else {}
        if len(sub) == 1 and next(iter(sub)) == "keyword":
            negated_keywords.add(str(next(iter(sub.values()))).strip().lower())
        return
    if key in ("method_family", "job"):
        out.setdefault(key, set()).update(str(v).lower() for v in value or [])
    elif key in ("mult", "charge"):
        try:
            op, number = _parse_comparison(value, key)
        except ValueError:
            return
        out.setdefault(key, []).append((op, number))
    elif key == "orca":
        out.setdefault("orca", []).append(value)
    elif key == "keyword":
        out.setdefault("keywords", set()).add(str(value).strip().lower())


def _range_conflicting(bounds_a: list, bounds_b: list) -> bool:
    """Whether two lists of (op, number) bounds cannot both hold."""
    def clamp(op: str, n: int) -> tuple[int, int]:
        lo, hi = (n, n) if op == "==" else ((n, 10**9) if op in (">", ">=") else (-10**9, n))
        if op == ">":
            lo += 1
        if op == "<":
            hi -= 1
        return lo, hi
    a_lo, a_hi = max((clamp(*b)[0] for b in bounds_a), default=-10**9), \
        min((clamp(*b)[1] for b in bounds_a), default=10**9)
    b_lo, b_hi = max((clamp(*b)[0] for b in bounds_b), default=-10**9), \
        min((clamp(*b)[1] for b in bounds_b), default=10**9)
    return a_lo > b_hi or b_lo > a_hi


def _disjoint(a: Any, b: Any) -> bool:
    """Whether two ``applies`` conditions provably never overlap. Only clear
    contradictions count; anything else is treated as overlapping (the graph
    check may then report a pair that in practice never co-occurs)."""
    ca: dict = {}
    cb: dict = {}
    na: set = set()
    nb: set = set()
    _flatten_conditions(a if a is not None else {}, ca, na)
    _flatten_conditions(b if b is not None else {}, cb, nb)
    for field in ("method_family", "job"):
        sa, sb = ca.get(field), cb.get(field)
        if sa and sb and not (sa & sb):
            return True
    for field in ("mult", "charge"):
        if field in ca and field in cb and _range_conflicting(ca[field], cb[field]):
            return True
    ka, kb = ca.get("keywords"), cb.get("keywords")
    # a keyword required by one condition and excluded by the other (or
    # required and negated inside the same condition) is a contradiction
    if ka and (ka & na or ka & nb):
        return True
    if kb and (kb & nb or kb & na):
        return True
    return False


def _check_of(card: dict, key: str) -> list[str]:
    check = card.get("check") or {}
    values = check.get(key) or []
    return [str(v).lower() for v in values]


def check_graph(cards: Iterable[dict] | None = None) -> list[str]:
    """Problems of the card graph: unknown ids in relations, ``requires``
    cycles, ``conflicts_with`` not mutual, and two verified cards where one
    requires and the other forbids the same keyword under overlapping
    conditions."""
    if cards is None:
        from . import cards as load_cards
        cards = load_cards(include_unverified=True)
    known: dict[str, dict] = {}
    problems: list[str] = []
    for card in cards:
        cid = str(card.get("id") or "")
        if cid:
            known[cid] = card
    # unknown ids + non-mutual conflicts_with
    for cid, card in known.items():
        for kind, targets in (card.get("relations") or {}).items():
            for target in targets or []:
                if target not in known:
                    problems.append(f"{cid}: {kind} names unknown card {target}")
        conflicts = set(card.get("relations", {}).get("conflicts_with") or [])
        for other in conflicts:
            back = set((known.get(other) or {}).get("relations", {}).get("conflicts_with") or [])
            if cid not in back:
                problems.append(f"{cid}: conflicts_with {other} is not mutual")
    # requires cycles
    WHITE, GREY, BLACK = 0, 1, 2
    colour = {cid: WHITE for cid in known}

    def walk(cid: str, stack: list[str]) -> None:
        colour[cid] = GREY
        for dep in (known[cid].get("relations") or {}).get("requires") or []:
            if dep not in known:
                continue
            if colour[dep] == GREY:
                cycle = stack[stack.index(dep):] + [dep]
                problems.append("requires cycle: " + " -> ".join(cycle))
            elif colour[dep] == WHITE:
                walk(dep, stack + [dep])
        colour[cid] = BLACK

    for cid in known:
        if colour[cid] == WHITE:
            walk(cid, [cid])
    # verified cards requiring and forbidding the same keyword under
    # overlapping conditions
    verified = [c for c in known.values() if c.get("status") == "verified"]
    for i, a in enumerate(verified):
        for b in verified[i + 1:]:
            if _disjoint(a.get("applies"), b.get("applies")):
                continue
            a_req, a_forb = set(_check_of(a, "requires_keyword")), set(_check_of(a, "forbids_keyword"))
            b_req, b_forb = set(_check_of(b, "requires_keyword")), set(_check_of(b, "forbids_keyword"))
            both = (a_req & b_forb) | (b_req & a_forb)
            if both:
                names = ", ".join(sorted(both))
                problems.append(f"{a['id']} and {b['id']}: verified cards require and "
                                f"forbid the same keyword under overlapping conditions ({names})")
    return problems


__all__ = ["applies", "findings", "check_graph", "evaluate_check",
           "condition_problems", "check_problems", "relations_problems",
           "CHECK_KEYS", "RELATION_KEYS"]
