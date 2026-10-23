"""An ORCA job as the ORCA manual defines it, and a setting merged into it.

The rules, from the ORCA 6.1.1 manual, section 2.1:

* ``!`` lines: "all the simple input lines are collected into a single
  string", and known keywords are looked up "in a predefined order,
  regardless of the order in the input file" -- except basis sets and
  auxiliary basis sets of one type, where the latter takes priority.  Measured
  on 6.1.1: ``TightSCF VeryTightSCF`` and ``VeryTightSCF TightSCF`` both run
  VeryTight, and two functionals are an INPUT ERROR.  So a keyword that is to
  replace another must take its place, not be appended.
* Blocks run from ``%name`` to ``end``.  "Variable assignments within blocks
  have the following general structure: VariableName Value", with an optional
  ``=``; arrays as ``Name[i] Value``; "some input block keywords open a nested
  sub-block, which must be closed with an additional end".  "If a keyword is
  duplicated, the latter value is used", but "it is not recommended to have
  multiple instances of the same block", since some blocks reset data when
  they are opened.  So a setting for a block the job already has is merged
  into that block: the same variable is replaced, a new one appended.
* Table 2.2: ``base``, ``cclib``, ``id``, ``ljcoefficients``, ``maxcore``,
  ``moinp`` and ``pointcharges`` are ``%`` settings with no ``end``.
* Table 2.1 synonyms: ``cis``/``tddft``, ``ice``/``iceci``/``cipsi``,
  ``symmetry``/``sym``.
* Comments run from ``#`` to the end of the line; input is not case
  sensitive, file names excepted.

:func:`parse_job` keeps every line it read, so :func:`render_job` gives the
same text back byte for byte when nothing was merged -- checked on the 35 294
ORCA inputs of the archive.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional, Sequence

__all__ = [
    "OrcaInputError", "Statement", "Block", "Item",
    "parse_job", "render_job", "parse_settings", "merge_settings", "split_jobs",
    "canonical_block_name", "ONE_LINE_SETTINGS", "KNOWN_BLOCKS",
    "job_keywords", "add_keywords", "remove_keywords", "replace_in_family", "keyword_family",
    "block_value", "set_block_values", "remove_block_values", "has_block",
    "setting_value", "set_setting", "remove_setting",
    "geometry_atoms", "with_coordinates", "apply_to_job",
]

#: Table 2.2 of the manual: '%' keywords with no closing 'end'.
ONE_LINE_SETTINGS = frozenset({"base", "cclib", "id", "ljcoefficients", "maxcore", "moinp", "pointcharges"})

#: Table 2.1 of the manual, synonyms mapped onto one name.
_SYNONYMS = {"tddft": "cis", "iceci": "ice", "cipsi": "ice", "sym": "symmetry"}
KNOWN_BLOCKS = frozenset({
    "autoci", "basis", "casresp", "casscf", "chelpg", "cim", "cis", "compound", "conical", "coords",
    "cosmors", "cpcm", "docker", "eda", "elprop", "eprnmr", "esd", "frag", "freq", "geom", "goat",
    "ice", "irc", "lft", "loc", "mcrpa", "md", "mdci", "mecp", "method", "mm", "mp2", "mrcc", "mrci",
    "mtr", "nbo", "ndoparas", "neb", "numgrad", "output", "pal", "paras", "plots", "qmmm", "rel",
    "rocis", "rr", "scf", "shark", "solvator", "symmetry", "vpt2", "xtb",
})

#: Sub-block openers that take arguments on their own line (the manual's
#: %basis examples); every other nested sub-block is opened by a line holding
#: only its keyword (SOSCF, Constraints, Scan, ...).
_ARGUMENT_OPENERS = frozenset({"newgto", "newauxjgto", "newauxcgto", "newauxjkgto", "newauxgto",
                               "newecp", "addgto", "addauxgto"})

#: %compound steps pair New_Step with Step_End (manual, section 8.2).
_PAIRED_CLOSERS = {"new_step": "step_end"}

_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class OrcaInputError(ValueError):
    """The text is not an ORCA job this module can take apart safely."""


def canonical_block_name(name: str) -> str:
    name = name.strip().lstrip("%").lower()
    return _SYNONYMS.get(name, name)


def _tokens(line: str) -> List[str]:
    """The line's words without its comment."""
    return line.split("#", 1)[0].split()


@dataclass
class Statement:
    """One assignment (one line) or one nested sub-block (opener to its end); key None for blank/comment."""

    key: Optional[str]
    lines: List[str]


@dataclass
class Block:
    name: str                       # canonical
    statements: List[Statement]
    raw: List[str]                  # the lines as read
    changed: bool = False
    one_line: bool = False          # written as "%name ... end" on one line


@dataclass
class Item:
    kind: str                       # "simple" | "block" | "setting" | "geometry" | "text"
    lines: List[str]
    name: str = ""                  # canonical block / setting name
    block: Optional[Block] = None
    changed: bool = False


def _statement_key(tokens: Sequence[str]) -> str:
    head = tokens[0].split("=", 1)[0].lower()
    if head in _ARGUMENT_OPENERS and len(tokens) > 1:
        return f"{head} {tokens[1].lower()}"
    return head


_ASSIGNED = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\[\d+\])?$")


def _assignments(tokens: Sequence[str]) -> List[List[str]]:
    """The variable assignments on one line, each as its tokens.

    ORCA reads ``maxiter 300 tole 1e-9`` as two assignments (measured on
    6.1.1), so a merge must see two variables there.  Only plain ``Name
    Value`` pairs (or ``Name=Value`` words) are split; anything else --
    quoted strings, lists, ``Scan``/``{ }`` syntax -- stays one statement.
    """
    tokens = list(tokens)
    if len(tokens) < 2 or any('"' in t or "'" in t for t in tokens):
        return [tokens]
    if all("=" in t and _ASSIGNED.match(t.split("=", 1)[0]) and t.split("=", 1)[1] for t in tokens):
        return [[t] for t in tokens]
    if len(tokens) >= 4 and len(tokens) % 2 == 0:
        names, values = tokens[0::2], tokens[1::2]
        if all(_ASSIGNED.match(n) and n.lower() not in _ARGUMENT_OPENERS for n in names) and \
                not any(v in ("=", "{", "}") or v.startswith("{") for v in values):
            return [[n, v] for n, v in zip(names, values)]
    return [tokens]


def _statements_of_line(line: str) -> List[Statement]:
    tokens = _tokens(line)
    parts = _assignments(tokens)
    if len(parts) == 1:
        return [Statement(_statement_key(tokens), [line])]
    return [Statement(_statement_key(part), ["  " + " ".join(part) + "\n"]) for part in parts]


def _opens_sub_block(tokens: Sequence[str]) -> bool:
    """A line that opens a nested sub-block: a bare keyword (SOSCF, Constraints) or NewGTO & co.

    A single word that is not a plain identifier -- ``ElDens("x.cube");`` in
    %plots, ``{C 44 C}`` in a constraint list -- is a statement, not an opener.
    """
    if not tokens or tokens[0].lower() in ("end", "step_end"):
        return False
    closed_here = len(tokens) > 1 and tokens[-1].lower() == "end"
    if closed_here:
        return False
    if len(tokens) == 1:
        return bool(_IDENTIFIER.match(tokens[0]))
    return tokens[0].lower() in _ARGUMENT_OPENERS


def _closes(tokens: Sequence[str], opener: str) -> bool:
    closer = _PAIRED_CLOSERS.get(opener, "end")
    return len(tokens) == 1 and tokens[0].lower() == closer


def _parse_body(lines: List[str], start: int, name: str) -> tuple[List[Statement], int]:
    """Statements from ``lines[start:]`` up to the block's own 'end'; returns (statements, index after end)."""
    statements: List[Statement] = []
    i = start
    while i < len(lines):
        tokens = _tokens(lines[i])
        if not tokens:
            statements.append(Statement(None, [lines[i]]))
            i += 1
            continue
        if len(tokens) == 1 and tokens[0].lower() == "end":
            return statements, i + 1
        if tokens[0].startswith("%"):
            # ORCA reads this as an unknown identifier of the open block too
            raise OrcaInputError(f"{tokens[0]} starts inside %{name}, which is not closed with end")
        if _opens_sub_block(tokens):
            stack, j = [tokens[0].lower()], i + 1
            while j < len(lines) and stack:
                inner = _tokens(lines[j])
                if inner and _closes(inner, stack[-1]):
                    stack.pop()
                elif _opens_sub_block(inner):
                    stack.append(inner[0].lower())
                j += 1
            depth = len(stack)
            if depth:
                raise OrcaInputError(f"%{name}: sub-block {tokens[0]} is not closed with end")
            statements.append(Statement("sub " + _statement_key(tokens), lines[i:j]))
            i = j
            continue
        statements.extend(_statements_of_line(lines[i]))
        i += 1
    raise OrcaInputError(f"%{name} is not closed with end")


def parse_job(text: str) -> List[Item]:
    """One ORCA job (no $new_job inside) as items, every line kept."""
    lines = text.splitlines(keepends=True)
    items: List[Item] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        tokens = _tokens(line)
        if not stripped or stripped.startswith("#"):
            items.append(Item("text", [line]))
            i += 1
        elif stripped.startswith("!"):
            items.append(Item("simple", [line]))
            i += 1
        elif stripped == "*":
            # a stray terminator after the coordinates; ORCA ignores it
            items.append(Item("text", [line]))
            i += 1
        elif stripped.startswith("*"):
            words = stripped[1:].split()
            kind = words[0].lower() if words else ""
            if kind.endswith("file") or kind in ("pdbfile",):
                items.append(Item("geometry", [line]))
                i += 1
            else:
                j = i + 1
                while j < len(lines) and lines[j].strip() != "*" and not lines[j].rstrip().endswith("*"):
                    j += 1
                if j >= len(lines):
                    raise OrcaInputError("coordinates are not closed with *")
                items.append(Item("geometry", lines[i:j + 1]))
                i = j + 1
        elif stripped.startswith("%"):
            name = canonical_block_name(tokens[0])
            if name == "compound":
                # the compound language (New_Step, Step_End, loops) is not merged into
                raise OrcaInputError("%compound jobs are left as they are")
            if name in ONE_LINE_SETTINGS:
                items.append(Item("setting", [line], name=name))
                i += 1
                continue
            rest = tokens[1:]
            if rest and rest[-1].lower() == "end":
                # the whole block on one line: "%pal nprocs 4 end", "%scf maxiter 300 end"
                body = " ".join(rest[:-1])
                statements = _statements_of_line("  " + body + "\n") if body else []
                items.append(Item("block", [line], name=name, block=Block(name, statements, [line], one_line=True)))
                i += 1
                continue
            first_body = []
            if rest:  # "%geom MaxIter 60" followed by more lines
                first_body = ["  " + " ".join(rest) + "\n"]
            body_lines = first_body + lines[i + 1:]
            statements, after = _parse_body(body_lines, 0, name)
            consumed = after - len(first_body)
            raw = lines[i:i + 1 + consumed]
            items.append(Item("block", raw, name=name, block=Block(name, statements, raw)))
            i += 1 + consumed
        else:
            items.append(Item("text", [line]))
            i += 1
    return items


def _render_block(block: Block) -> List[str]:
    name = block.name if block.name != "cis" else "tddft"
    if block.one_line and all(len(st.lines) == 1 for st in block.statements):
        # kept on one line as it was written ("%pal nprocs 4 end"): ORCA reads
        # several assignments on a line, and DELFIN's own editors of %pal and
        # %scf look for that form
        body = " ".join(" ".join(_tokens(st.lines[0])) for st in block.statements if st.key)
        return [f"%{name} {body} end\n" if body else f"%{name} end\n"]
    out = [f"%{name}\n"]
    for st in block.statements:
        out.extend(st.lines)
    out.append("end\n")
    return out


def render_job(items: Sequence[Item]) -> str:
    out: List[str] = []
    for item in items:
        if item.kind == "block" and item.block is not None and item.block.changed:
            out.extend(_render_block(item.block))
        else:
            out.extend(item.lines)
    return "".join(out)


def parse_settings(text: str) -> List[Item]:
    """What a user wrote to be merged: blocks, one-line settings and '!' lines, nothing else.

    Coordinates and free text are refused; a block must be closed.
    """
    items = [it for it in parse_job(text if text.endswith("\n") else text + "\n") if it.kind != "text"
             or it.lines[0].strip() and not it.lines[0].strip().startswith("#")]
    for it in items:
        if it.kind == "geometry":
            raise OrcaInputError("coordinates cannot be set this way; DELFIN writes each job's geometry")
        if it.kind == "text":
            raise OrcaInputError(f"not an ORCA block, % setting or ! line: {it.lines[0].strip()!r}")
    return items


def _normalised_line(line: str) -> str:
    return "  " + line.strip() + "\n"


def merge_settings(items: List[Item], settings: Sequence[Item], *, only_existing_blocks: bool = False) -> List[Item]:
    """``settings`` (from :func:`parse_settings`) merged into ``items``; returns the new item list.

    A block the job has gets each variable replaced or appended; one it does
    not have is added before the coordinates -- unless ``only_existing_blocks``,
    for settings aimed at many jobs at once, where a block DELFIN did not write
    (a %tddft in an optimisation) would change what the job is.  A one-line
    setting replaces the job's own.  '!' lines are not handled here.
    """
    items = list(items)

    def insert_before_geometry(new_item: Item) -> None:
        for idx, it in enumerate(items):
            if it.kind == "geometry":
                items.insert(idx, new_item)
                return
        items.append(new_item)

    for setting in settings:
        if setting.kind == "setting":
            for idx, it in enumerate(items):
                if it.kind == "setting" and it.name == setting.name:
                    items[idx] = Item("setting", [setting.lines[0].strip() + "\n"], name=setting.name, changed=True)
                    break
            else:
                insert_before_geometry(Item("setting", [setting.lines[0].strip() + "\n"], name=setting.name,
                                            changed=True))
        elif setting.kind == "block":
            existing = _one_block(items, setting.name)
            target = Item("block", [], name=setting.name, block=existing) if existing is not None else None
            if target is None:
                if only_existing_blocks:
                    continue
                block = Block(setting.name, [Statement(st.key, [_normalised_line(l) for l in st.lines]
                                                       if st.key and not st.key.startswith("sub ") else st.lines)
                                             for st in setting.block.statements if st.key], [], changed=True)
                insert_before_geometry(Item("block", [], name=setting.name, block=block, changed=True))
                continue
            block = target.block
            for st in setting.block.statements:
                if not st.key:
                    continue
                new_lines = st.lines if st.key.startswith("sub ") else [_normalised_line(st.lines[0])]
                for pos, existing in enumerate(block.statements):
                    if existing.key == st.key:
                        if [l.strip() for l in existing.lines] != [l.strip() for l in new_lines]:
                            block.statements[pos] = Statement(st.key, new_lines)
                            block.changed = True
                        break
                else:
                    block.statements.append(Statement(st.key, new_lines))
                    block.changed = True
    return items


# ---------------------------------------------------------------- '!' keywords

#: Keyword families from the manual's simple-input tables (SCF convergence
#: Table 2.11, grids 2.48, RI/COSX 2.44, dispersion 3.13, relativity 2.54,
#: guesses 2.68, run types 2.4).  Within a family ORCA does not take the last
#: one written, so a keyword replaces the family member already there.
_FAMILIES = {
    "scf convergence": {"sloppyscf", "loosescf", "normalscf", "strongscf", "tightscf", "verytightscf",
                        "extremescf", "scfconv6", "scfconv7", "scfconv8", "scfconv9", "scfconv10"},
    "integration grid": {"defgrid1", "defgrid2", "defgrid3"},
    "reference": {"rhf", "uhf", "rohf", "rks", "uks", "roks"},
    "RI approximation": {"ri", "rijcosx", "rijk", "rijonx", "rijdx", "nori"},
    "dispersion": {"d2", "d3", "d3bj", "d3zero", "d4"},
    "relativity": {"zora", "dkh", "dkh2", "x2c"},
    "guess": {"pmodel", "patom", "hueckel", "hcore", "moread"},
    "optimisation level": {"opt", "normalopt", "tightopt", "looseopt", "verytightopt", "sloppyopt",
                           "crudeopt", "copt", "zopt"},
    "frequency": {"freq", "numfreq", "anfreq"},
}
#: A job's run type is DELFIN's to decide: these may replace a member of
#: their family the job already has (TightOpt for Opt, NumFreq for Freq), but
#: are never added to a job that does not run that kind of calculation.
_RUN_TYPE_FAMILIES = {"optimisation level", "frequency"}
#: Keywords that say what method a job is.  Aimed at one job they are the
#: user's to set; aimed at many they only replace what a job already has.
_METHOD_FAMILIES = {"functional", "orbital basis", "auxiliary basis /J", "auxiliary basis /JK",
                    "auxiliary basis /C", "reference", "RI approximation", "dispersion", "relativity",
                    "implicit solvation"}
#: Never settable: the scheduler sets the cores, DELFIN chains the orbitals,
#: and these change what a job is.
_REFUSED_KEYWORDS = {"moread", "noautostart", "optts", "scants", "neb", "neb-ts", "neb-ci", "irc", "md",
                     "goat", "engrad", "numgrad", "sp", "energy"}
_REFUSED_PREFIXES = ("pal", "esd(", "compound")

_SOLVATION_PREFIXES = ("cpcm(", "smd(", "cpcmc(", "alpb(", "ddcosmo(", "cosmors(")


def _functionals_and_bases():
    try:
        from delfin.common.control_validator import ORCA_BASIS_SETS, ORCA_FUNCTIONALS
    except Exception:  # noqa: BLE001 - families still work without the curated lists
        return frozenset(), frozenset()
    squash = lambda s: re.sub(r"[-_\s]", "", s.lower())  # noqa: E731
    return frozenset(squash(f) for f in ORCA_FUNCTIONALS), frozenset(b.lower() for b in ORCA_BASIS_SETS)


_FUNCTIONALS, _BASES = _functionals_and_bases()


def keyword_family(token: str) -> Optional[str]:
    low = token.lower()
    for family, members in _FAMILIES.items():
        if low in members:
            return family
    if low.startswith(_SOLVATION_PREFIXES):
        return "implicit solvation"
    if low.endswith("/jk"):
        return "auxiliary basis /JK"
    if low.endswith("/j"):
        return "auxiliary basis /J"
    if low.endswith("/c"):
        return "auxiliary basis /C"
    if low in _BASES:
        return "orbital basis"
    if re.sub(r"[-_\s]", "", low) in _FUNCTIONALS or low.startswith("libxc("):
        return "functional"
    return None


def refused_keyword(token: str) -> Optional[str]:
    """Why ``token`` cannot be set through CONTROL, or None."""
    low = token.lower()
    if low in _REFUSED_KEYWORDS or low.startswith(_REFUSED_PREFIXES) or re.match(r"^pal\d+$", low):
        return (f"{token} is DELFIN's to set: cores come from the scheduler, orbitals from the "
                f"previous step, and the run type from the workflow")
    return None


def _split_bang(line: str) -> tuple[List[str], str]:
    body, _, comment = line.rstrip("\n").lstrip()[1:].partition("#")
    return body.split(), ("#" + comment if comment else "")


def _bang_words(items: Sequence[Item]) -> dict:
    return {idx: _split_bang(items[idx].lines[0]) for idx, it in enumerate(items) if it.kind == "simple"}


def _store_bang_words(items: List[Item], words: dict) -> None:
    for idx, (tokens_here, comment) in words.items():
        old = items[idx].lines[0]
        new_line = ("! " + " ".join(tokens_here) + (" " + comment if comment else "")).rstrip() + "\n"
        if [t.lower() for t in tokens_here] != [t.lower() for t in _split_bang(old)[0]]:
            items[idx] = Item("simple", [new_line], changed=True)


def job_keywords(items: Sequence[Item]) -> List[str]:
    """Every '!' keyword of the job, in the order written."""
    return [w for tokens_here, _ in _bang_words(items).values() for w in tokens_here]


def add_keywords(items: List[Item], tokens: Sequence[str]) -> bool:
    """``tokens`` the job does not have yet, appended to its last '!' line; True if any was added."""
    words = _bang_words(items)
    if not words:
        return False
    present = {w.lower() for tokens_here, _ in words.values() for w in tokens_here}
    last = max(words)
    added = False
    for token in tokens:
        if token.lower() not in present:
            words[last][0].append(token)
            present.add(token.lower())
            added = True
    _store_bang_words(items, words)
    return added


def remove_keywords(items: List[Item], tokens: Sequence[str]) -> bool:
    """``tokens`` (any case) taken off every '!' line; True if any was there."""
    drop = {t.lower() for t in tokens}
    words = _bang_words(items)
    removed = False
    for idx, (tokens_here, comment) in words.items():
        kept = [w for w in tokens_here if w.lower() not in drop]
        removed = removed or len(kept) != len(tokens_here)
        words[idx] = (kept, comment)
    _store_bang_words(items, words)
    return removed


def replace_in_family(items: List[Item], token: str, *, add_if_absent: bool = True) -> bool:
    """``token`` in place of the member of its family the job has (first one kept in position).

    With no member present it is appended, unless ``add_if_absent`` is False.
    A keyword with no family is only added.  True if the job changed.
    """
    family = keyword_family(token)
    words = _bang_words(items)
    if not words:
        return False
    if family is None:
        return add_keywords(items, [token]) if add_if_absent else False
    placed = False
    changed = False
    for idx in sorted(words):
        tokens_here, comment = words[idx]
        out = []
        for w in tokens_here:
            if keyword_family(w) == family:
                if not placed:
                    out.append(w if w.lower() == token.lower() else token)
                    changed = changed or w.lower() != token.lower()
                    placed = True
                else:
                    changed = True
                continue
            out.append(w)
        words[idx] = (out, comment)
    if not placed:
        if not add_if_absent:
            return False
        words[max(words)][0].append(token)
        changed = True
    _store_bang_words(items, words)
    return changed


def merge_keywords(items: List[Item], tokens: Sequence[str], *, many_jobs: bool = False) -> tuple[List[Item], List[str]]:
    """CONTROL's ``tokens`` into the job's '!' lines; (items, notes on what was left out).

    A keyword of a family takes the place of the member the job has.  One
    that the job has no member of is added -- except a run type (DELFIN
    decides whether a job optimises or computes frequencies) and, for
    settings aimed at many jobs, a keyword that defines the method: a
    functional or basis added to an xTB job, or a dispersion correction to a
    functional that has its own, would change what that job is.
    """
    items = list(items)
    notes: List[str] = []
    if not any(it.kind == "simple" for it in items):
        return items, [f"no '!' line to put {' '.join(tokens)} on"]
    for token in tokens:
        reason = refused_keyword(token)
        if reason:
            notes.append(reason)
            continue
        family = keyword_family(token)
        if family == "guess" and any(w.lower() == "moread" for w in job_keywords(items)):
            notes.append(f"{token} not set: this job starts from the orbitals of the step before it (MORead)")
            continue
        add = family not in _RUN_TYPE_FAMILIES and not (many_jobs and family in _METHOD_FAMILIES)
        if not replace_in_family(items, token, add_if_absent=add) and not add:
            if not any(keyword_family(w) == family for w in job_keywords(items)):
                why = (f"this job runs no {family} step" if family in _RUN_TYPE_FAMILIES
                       else f"this job has no {family} to replace")
                notes.append(f"{token} not added: {why}")
    return items, notes


# ------------------------------------------------------------- what is DELFIN's

#: % settings and blocks DELFIN writes per job and a user must not replace.
_REFUSED_SETTINGS = {
    "base": "each job's file names are DELFIN's; other steps read them",
    "moinp": "DELFIN hands each step the orbitals of the one before it",
    "maxcore": "memory is set from CONTROL's maxcore for every job",
}
_REFUSED_BLOCKS = {
    "pal": "cores are given to each job by the scheduler (PAL in CONTROL)",
    "coords": "DELFIN writes each job's geometry",
    "compound": "compound jobs are not merged into",
}
#: Blocks that only hold settings for what a job runs anyway; everything
#: else in Table 2.1 starts a calculation of its own when it is present.
_SETTINGS_BLOCKS = frozenset({"scf", "geom", "freq", "output", "basis", "method", "rel", "shark", "mp2",
                              "symmetry"})
#: A block whose calculation a '!' keyword already runs: %cpcm in a job
#: with CPCM(...) only tunes it.  A bare %cpcm switches CPCM on (measured on
#: 6.1.1), so it is never added to a gas-phase job.
_BLOCK_RUN_BY_KEYWORD = {"cpcm": "implicit solvation"}
#: %tddft keywords the ESD module sets per job; for settings aimed at many
#: jobs also the ones TDDFT_* keys own, so a value has one spelling.
_PER_JOB_TDDFT = {"iroot", "irootmult", "irootlist", "triplets", "sroot", "troot", "trootssl", "nacme", "etf"}
_TDDFT_OWNED = {"nroots": "TDDFT_nroots", "maxdim": "TDDFT_maxdim", "maxiter": "TDDFT_maxiter",
                "tda": "TDDFT_TDA", "followiroot": "TDDFT_followiroot", "dosoc": "TDDFT_SOC"}


def refused_settings(settings: Sequence[Item], *, many_jobs: bool) -> List[str]:
    """What in ``settings`` DELFIN will not merge, with the reason."""
    problems: List[str] = []
    for it in settings:
        if it.kind == "setting" and it.name in _REFUSED_SETTINGS:
            problems.append(f"%{it.name}: {_REFUSED_SETTINGS[it.name]}")
        elif it.kind == "block":
            if it.name in _REFUSED_BLOCKS:
                problems.append(f"%{it.name}: {_REFUSED_BLOCKS[it.name]}")
            elif it.name not in KNOWN_BLOCKS:
                problems.append(f"%{it.name} is not an ORCA input block (manual Table 2.1)")
            elif it.name == "cis":
                for st in it.block.statements:
                    key = (st.key or "").split()[0] if st.key else ""
                    if key in _PER_JOB_TDDFT:
                        problems.append(f"%tddft {key}: the ESD module sets it per job")
                    elif many_jobs and key in _TDDFT_OWNED:
                        problems.append(f"%tddft {key}: use {_TDDFT_OWNED[key]}, which reaches every TD-DFT job")
        elif it.kind == "simple":
            for token in it.lines[0].strip()[1:].split("#", 1)[0].split():
                reason = refused_keyword(token)
                if reason:
                    problems.append(reason)
    return problems


def apply_to_job(text: str, *, keywords: Sequence[str] = (), additions: Sequence[str] = (),
                 many_jobs: bool = False) -> tuple[str, List[str]]:
    """One job with CONTROL's keyword:/additions: merged in; (new text, notes on what was left out).

    A job this module cannot take apart safely comes back unchanged, with a
    note -- nothing is merged into text that is not understood.
    """
    try:
        items = parse_job(text)
    except OrcaInputError as exc:
        return text, [f"left unchanged: {exc}"]
    notes: List[str] = []
    tokens = list(keywords)
    for value in additions:
        try:
            settings = parse_settings(value)
        except OrcaInputError as exc:
            notes.append(f"not merged: {exc}")
            continue
        refused = refused_settings(settings, many_jobs=many_jobs)
        if refused:
            notes.extend(refused)
            settings = [s for s in settings if not refused_settings([s], many_jobs=many_jobs)]
        for s in settings:
            if s.kind == "simple":
                tokens.extend(s.lines[0].strip()[1:].split("#", 1)[0].split())
        merged = [s for s in settings if s.kind in ("block", "setting")]
        if not many_jobs and any(s.kind == "block" and s.name == "cis" for s in merged) \
                and not has_block(items, "cis") \
                and any(keyword_family(w) == "optimisation level" for w in job_keywords(items)):
            notes.append("%tddft added to an optimisation: ORCA now optimises the excited state IROOT "
                         "(default 1), not the ground state")
        if many_jobs:
            # a block that only holds settings may be added to any job; one
            # that starts a calculation (%tddft, %eprnmr, %casscf, ...) only
            # changes jobs that already run it
            runs = {name for name, family in _BLOCK_RUN_BY_KEYWORD.items()
                    if any(keyword_family(w) == family for w in job_keywords(items))}
            addable = _SETTINGS_BLOCKS | runs
            items = merge_settings(items, [s for s in merged if s.kind == "setting" or s.name in addable])
            items = merge_settings(items, [s for s in merged if s.kind == "block" and s.name not in addable],
                                   only_existing_blocks=True)
        else:
            items = merge_settings(items, merged)
    if tokens:
        items, more = merge_keywords(items, tokens, many_jobs=many_jobs)
        notes.extend(more)
    return render_job(items), notes


# ------------------------------------------------------------ jobs of a file

_NEW_JOB = re.compile(r"^\s*\$new_job\b", re.IGNORECASE)


def split_jobs(text: str) -> List[str]:
    """The jobs of an input file; each job after the first starts with its ``$new_job`` line.

    ``"".join(split_jobs(text)) == text``.
    """
    jobs: List[str] = []
    current: List[str] = []
    for line in text.splitlines(keepends=True):
        if _NEW_JOB.match(line):
            jobs.append("".join(current))
            current = [line]
        else:
            current.append(line)
    jobs.append("".join(current))
    return jobs


# ------------------------------------------------------------ block variables

def _put_statement(block: Block, statement: Statement) -> None:
    for pos, existing in enumerate(block.statements):
        if existing.key == statement.key:
            if [l.strip() for l in existing.lines] != [l.strip() for l in statement.lines]:
                block.statements[pos] = statement
                block.changed = True
            return
    block.statements.append(statement)
    block.changed = True


def has_block(items: Sequence[Item], name: str) -> bool:
    name = canonical_block_name(name)
    return any(it.kind == "block" and it.name == name for it in items)


def _one_block(items: List[Item], name: str) -> Optional[Block]:
    """The job's %name as one block, or None.

    ORCA reads a repeated block variable by variable, the later value winning
    (measured for %scf and %tddft on 6.1.1), and the manual advises against
    repeating a block.  So before one is changed, the later ones are folded
    into the first the way ORCA would have read them.
    """
    idxs = [i for i, it in enumerate(items) if it.kind == "block" and it.name == name]
    if not idxs:
        return None
    first = items[idxs[0]].block
    for i in idxs[1:]:
        for st in items[i].block.statements:
            if st.key:
                _put_statement(first, st)
        first.changed = True
    for i in reversed(idxs[1:]):
        del items[i]
    return first


def _value_text(value) -> str:
    if value is True:
        return "true"
    if value is False:
        return "false"
    return str(value)


def block_value(items: Sequence[Item], name: str, key: str) -> Optional[str]:
    """The value ORCA reads for ``%name key`` (the last one written), or None."""
    name, key = canonical_block_name(name), key.lower()
    found = None
    for it in items:
        if it.kind == "block" and it.name == name:
            for st in it.block.statements:
                if st.key == key:
                    words = _tokens(st.lines[0])
                    rest = " ".join(words[1:]) if "=" not in words[0] else words[0].split("=", 1)[1]
                    found = rest.lstrip("= ").strip()
    return found


def set_block_values(items: List[Item], name: str, values: dict, *, create: bool = True) -> bool:
    """Set ``%name`` variables (key -> value); the block is added before the coordinates if missing.

    A value must be a single ORCA value -- a word, a number, true/false or a
    quoted string; a bare keyword on its own line opens a sub-block in ORCA
    (manual section 2.1: ``SOSCF`` in %scf), so none is written.
    """
    name = canonical_block_name(name)
    block = _one_block(items, name)
    if block is None:
        if not create:
            return False
        block = Block(name, [], [], changed=True)
        for idx, it in enumerate(items):
            if it.kind == "geometry":
                items.insert(idx, Item("block", [], name=name, block=block, changed=True))
                break
        else:
            items.append(Item("block", [], name=name, block=block, changed=True))
    before = block.changed
    block.changed = False
    for key, value in values.items():
        text = _value_text(value).strip()
        if not text:
            raise OrcaInputError(f"%{name} {key}: no value given")
        _put_statement(block, Statement(key.lower(), [f"  {key} {text}\n"]))
    changed = block.changed
    block.changed = before or changed
    return changed


def remove_block_values(items: List[Item], name: str, keys: Sequence[str]) -> bool:
    """Take ``keys`` (and sub-blocks of that name) out of ``%name``; True if any was there."""
    name = canonical_block_name(name)
    drop = {k.lower() for k in keys}
    if not has_block(items, name):
        return False
    block = _one_block(items, name)
    kept = [st for st in block.statements
            if not (st.key and (st.key in drop or (st.key.startswith("sub ") and st.key[4:].split()[0] in drop)))]
    if len(kept) == len(block.statements):
        return False
    block.statements = kept
    block.changed = True
    return True


# ------------------------------------------------------- one-line % settings

def setting_value(items: Sequence[Item], name: str) -> Optional[str]:
    """What follows ``%name`` (e.g. the quoted file of %moinp), or None."""
    for it in items:
        if it.kind == "setting" and it.name == name:
            words = it.lines[0].split("#", 1)[0].split(None, 1)
            return words[1].strip() if len(words) > 1 else ""
    return None


def set_setting(items: List[Item], name: str, value: str) -> bool:
    line = f"%{name} {value}\n"
    for idx, it in enumerate(items):
        if it.kind == "setting" and it.name == name:
            if it.lines[0].split("#", 1)[0].split() == line.split():
                return False
            items[idx] = Item("setting", [line], name=name, changed=True)
            return True
    simple = [idx for idx, it in enumerate(items) if it.kind == "simple"]
    at = simple[-1] + 1 if simple else 0
    items.insert(at, Item("setting", [line], name=name, changed=True))
    return True


def remove_setting(items: List[Item], name: str) -> bool:
    idxs = [i for i, it in enumerate(items) if it.kind == "setting" and it.name == name]
    for i in reversed(idxs):
        del items[i]
    return bool(idxs)


# --------------------------------------------------------------- coordinates

def _element(label: str) -> str:
    match = re.match(r"[A-Za-z]{1,2}", label)
    return match.group(0).capitalize() if match else label


def geometry_atoms(items: Sequence[Item]) -> Optional[List[str]]:
    """The element of every atom of inline Cartesian coordinates, or None for any other geometry."""
    geom = next((it for it in items if it.kind == "geometry"), None)
    if geom is None or len(geom.lines) < 2:
        return None
    head = geom.lines[0].split()
    if len(head) < 2 or head[1].lower() != "xyz":
        return None
    rows = [l for l in geom.lines[1:] if l.strip() and l.strip() != "*"]
    return [_element(l.split()[0]) for l in rows]


def with_coordinates(items: List[Item], atoms: Sequence[tuple]) -> None:
    """Give the job's coordinates the positions of ``atoms`` ((element, x, y, z) each).

    Inline coordinates keep everything after x y z on each line, so a metal's
    ``NewGTO ... end`` stays with its atom.  ``* xyzfile`` becomes inline.
    Atoms that do not match the job's, one by one, are refused.
    """
    for idx, it in enumerate(items):
        if it.kind != "geometry":
            continue
        head = it.lines[0].split()
        if len(head) < 4 or head[0] != "*" or head[1].lower() not in ("xyz", "xyzfile"):
            raise OrcaInputError(f"coordinates given as {' '.join(head[:2])} are not replaced")
        if head[1].lower() == "xyzfile":
            rows = [f"  {el:<3s} {x:14.8f} {y:14.8f} {z:14.8f}\n" for el, x, y, z in atoms]
            items[idx] = Item("geometry", [f"* xyz {head[2]} {head[3]}\n", *rows, "*\n"], changed=True)
            return
        body = it.lines[1:]
        closing = body[-1] if body and body[-1].strip() == "*" else None
        rows = [l for l in (body[:-1] if closing else body) if l.strip()]
        if len(rows) != len(atoms) or any(_element(r.split()[0]) != _element(a[0]) for r, a in zip(rows, atoms)):
            raise OrcaInputError("the new positions are not for the atoms of this job")
        new_rows = []
        for row, (_, x, y, z) in zip(rows, atoms):
            parts = row.split()
            tail = " ".join(parts[4:])
            new_rows.append(f"  {parts[0]:<3s} {x:14.8f} {y:14.8f} {z:14.8f}" + (f" {tail}" if tail else "") + "\n")
        if not closing:
            new_rows[-1] = new_rows[-1].rstrip("\n").rstrip().rstrip("*").rstrip() + "\n"
            closing = "*\n"
        items[idx] = Item("geometry", [it.lines[0], *new_rows, closing], changed=True)
        return
    raise OrcaInputError("the job has no coordinates")
