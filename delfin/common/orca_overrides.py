"""CONTROL's ``keyword:<job>=`` and ``additions:<job>=`` entries: read, aimed and checked in one place.

The run (delfin.orca) and the dashboard's check (delfin.config) both read the
entries here, so what the check accepts is what a job gets.

``<job>`` is a job's name -- its ``%base``, else its input file name, with
``X`` and ``X_OCCUPIER`` naming the same job -- or the folder of an OCCUPIER
step (``initial_OCCUPIER``: every run in it), a pattern (``S*``,
``*_ISC*``), or ``all``.  How a value is merged into a job is
delfin.common.orca_input's.
"""

from __future__ import annotations

import ast
import fnmatch
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

from delfin.common import orca_input

__all__ = ["parse_override_text", "normalize_target", "target_aliases", "override_groups",
           "is_pattern", "findings"]

_OVERRIDE_KEY_RE = re.compile(r"^(keyword|addition|additions)\s*:\s*(.+)$", re.IGNORECASE)
_VALID_TARGET_RE = re.compile(r"^[A-Za-z0-9_.*?\[\]-]+$")

#: Names DELFIN gives jobs, by shape: redox steps and their OCCUPIER runs,
#: every ESD job (S1, S1_TDDFT, S1_second_deltaSCF, S1_T1_ISC_msp1, S0_IP,
#: T1_S0_PHOSP_iroot1, ...), the reorganisation-energy jobs, and the xTB ones.
_NAME_SHAPES = (
    re.compile(r"^(initial|(ox|red)_step_\d+)(_occupier)?$"),
    re.compile(r"^input\d*$"),
    re.compile(r"^[st]\d+(_[a-z0-9]+)*$"),
    re.compile(r"^e_n_(cation|anion)$"),
    re.compile(r"^xtb2?(_goat|_solvator)?$"),
)
_NAMED = {"basename", "genolate"}


def normalize_target(raw: str) -> str:
    text = str(raw or "").strip().strip('"').strip("'")
    text = Path(text).name
    if text.lower().endswith(".inp"):
        text = text[:-4]
    return text.strip().lower()


def target_aliases(raw: str) -> List[str]:
    """A job's names for matching: ``X`` and ``X_occupier`` are the same job."""
    norm = normalize_target(raw)
    if not norm:
        return []
    aliases = {norm}
    if norm.endswith("_occupier"):
        aliases.add(norm[:-9])
    else:
        aliases.add(norm + "_occupier")
    return [a for a in aliases if a]


def is_pattern(target: str) -> bool:
    return target == "all" or any(ch in target for ch in "*?[")


def _parse_value(lines: List[str], start_idx: int, initial_value: str) -> Tuple[Any, int]:
    """One value; a ``[`` opens a list that may run over several lines."""
    value = initial_value.strip()
    if value.startswith("["):
        buffer = value + "\n"
        depth = value.count("[") - value.count("]")
        idx = start_idx + 1
        while idx < len(lines) and depth > 0:
            line = lines[idx]
            buffer += line
            depth += line.count("[") - line.count("]")
            idx += 1
        try:
            return ast.literal_eval(buffer), idx
        except Exception:
            # ORCA input between the brackets, unquoted:
            #   additions:T1=[
            #   %scf
            #     maxiter 400
            #   end
            #   ]
            raw = buffer.strip()
            if raw.startswith("[") and raw.endswith("]"):
                inner = raw[1:-1].strip()
                return ([inner] if inner else []), idx
            return value, idx
    try:
        return ast.literal_eval(value), start_idx + 1
    except Exception:
        return value, start_idx + 1


def _flatten(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [text for item in value for text in _flatten(item)]
    text = str(value).strip()
    return [text] if text else []


def parse_override_text(text: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], List[Tuple[str, str, List[str]]]]:
    """(keyword values by target, addition values by target, every entry as (key, kind, values))."""
    lines = text.splitlines(keepends=True)
    keyword_map: Dict[str, List[str]] = {}
    addition_map: Dict[str, List[str]] = {}
    entries: List[Tuple[str, str, List[str]]] = []
    idx = 0
    while idx < len(lines):
        raw_line = lines[idx]
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or "=" not in raw_line:
            idx += 1
            continue
        key_raw, value_raw = raw_line.split("=", 1)
        key = key_raw.strip()
        match = _OVERRIDE_KEY_RE.match(key)
        if not match:
            idx += 1
            continue
        kind = "keyword" if match.group(1).lower() == "keyword" else "additions"
        parsed, next_idx = _parse_value(lines, idx, value_raw)
        values = _flatten(parsed)
        entries.append((key, kind, values))
        target = normalize_target(match.group(2))
        if target and values:
            (keyword_map if kind == "keyword" else addition_map).setdefault(target, []).extend(values)
        idx = next_idx
    return keyword_map, addition_map, entries


def override_groups(names: List[str], folder: str, keyword_map: Dict[str, List[str]],
                    addition_map: Dict[str, List[str]]) -> List[Tuple[bool, List[str], List[str]]]:
    """The values that reach a job with ``names`` in ``folder``, as (aimed at many jobs, keywords, additions).

    ``all`` first, then patterns, then the job's own name, so the more
    specific is merged later and wins.  One target's additions are one piece
    of input, as written.
    """
    groups: List[Tuple[bool, List[str], List[str]]] = []
    targets = list(dict.fromkeys(list(keyword_map) + list(addition_map)))
    folder = folder.lower()

    def collect(selected: List[str], many: bool) -> None:
        kws = [tok for t in selected for value in keyword_map.get(t, []) for tok in str(value).split()]
        adds = ["\n".join(addition_map[t]) for t in selected if addition_map.get(t)]
        if kws or adds:
            groups.append((many, kws, adds))

    collect([t for t in targets if t == "all"], True)
    collect([t for t in targets if t != "all" and is_pattern(t)
             and any(fnmatch.fnmatchcase(n, t) for n in names)], True)
    collect([t for t in targets if not is_pattern(t)
             and (t in names or (folder.endswith("_occupier") and t == folder))], False)
    return groups


def _known(target: str) -> bool:
    return is_pattern(target) or target in _NAMED or any(shape.match(target) for shape in _NAME_SHAPES)


def findings(text: str) -> Tuple[List[str], List[str]]:
    """(errors, hints) for the entries of a CONTROL text.

    An error is what no job could run with: a malformed key, or input ORCA
    cannot read.  What would merely not take effect -- a name no job has, a
    setting that is DELFIN's to make -- is a hint, so a file that ran before
    still runs.
    """
    errors: List[str] = []
    hints: List[str] = []
    _, _, entries = parse_override_text(text)
    for key, kind, values in entries:
        raw_target = _OVERRIDE_KEY_RE.match(key).group(2).strip()
        if not _VALID_TARGET_RE.match(raw_target):
            errors.append(f"Invalid ORCA override key: {key!r} (expected additions:<job>=[...] or "
                          f"keyword:<job>=[...]; <job> is a job name, a pattern such as S*, or all)")
            continue
        if not values:
            continue
        target = normalize_target(raw_target)
        many = is_pattern(target)
        if not _known(target):
            hints.append(f"{key}: DELFIN writes no job named {raw_target!r}, so this changes nothing "
                         f"(names: initial, ox_step_1, red_step_1, S0, S1, T1, S1_TDDFT, ...; "
                         f"patterns such as S* or all)")
        tokens: List[str] = []
        if kind == "keyword":
            tokens = [tok for value in values for tok in value.split()]
        else:
            try:
                settings = orca_input.parse_settings("\n".join(values))
            except orca_input.OrcaInputError as exc:
                errors.append(f"{key}: not ORCA input ({exc})")
                continue
            hints.extend(f"{key}: {reason}; it is left out"
                         for reason in orca_input.refused_settings(
                             [s for s in settings if s.kind != "simple"], many_jobs=many))
            if not many and not re.match(r"^[st][1-9]\d*$", target) and \
                    any(s.kind == "block" and s.name == "cis" for s in settings):
                hints.append(f"{key}: in a job that optimises, %tddft makes ORCA optimise the excited "
                             f"state IROOT (default 1) instead of the ground state")
            tokens = [tok for s in settings if s.kind == "simple"
                      for tok in s.lines[0].strip()[1:].split("#", 1)[0].split()]
        for token in tokens:
            reason = orca_input.refused_keyword(token)
            if reason:
                hints.append(f"{key}: {reason}; it is left out")
    return errors, list(dict.fromkeys(hints))
