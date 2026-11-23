"""Compare a run's own numbers against a reference the agent cannot author.

The execution contract (``run_evidence.contract``) answers *did this
run compute*. It cannot answer *is the number right*, and an agent can
write a confident report over a wrong number as easily as over a right
one. A reference published by someone else closes that: a model can
fabricate a report, and it cannot fabricate agreement with a value it
did not choose.

A field survey of 26 autonomous-research systems (arXiv:2608.05179)
found, under its own coding rule, that **no system of the LLM era
demonstrates an externally validated oracle inside the loop**. This is
that oracle, in the smallest honest form.

THE REFERENCE SET IS INPUT, NOT CONTENT. Nothing here embeds published
data, deliberately and for two reasons. The practical one: QUESTDB, the
obvious source for vertical excitation energies (GitHub pfloos/QUESTDB,
Zenodo DOI 10.5281/zenodo.15671385, Loos, Boggio-Pasqua, Blondel,
Lipparini and Jacquemin, J. Chem. Theory Comput. 2025, 21, 8010), is
CC-BY-SA-4.0 -- a share-alike licence, with EU database rights on top --
and vendoring it would put obligations on this repository that have
nothing to do with the comparison. The design one: a comparator that
takes its references as an argument works with ANY reference set, and a
comparator that embeds one works with exactly that one.

FOUR VERDICTS, NEVER TWO:

    agrees        the computed value is within the reference's own
                  stated tolerance
    differs       both values are known and outside it
    not_measured  the run produced no comparable value -- never an
                  agreement, for the same reason an unjudged contract
                  flag is never a pass
    ambiguous     a value exists and which transition it corresponds to
                  cannot be decided by the matching rule; reported
                  rather than guessed, because assigning a state by
                  proximity is how a benchmark scores luck

That last verdict is the honest part. Matching a computed transition to
a published one is the hard half of this comparison: state ordering
differs between methods, a dark state can sit below a bright one, and
"closest in energy" silently rewards a wrong assignment. So the matching
rule is declared per reference and named in every verdict.

Reading is done through ``delfin.api.extract_excited_states``, which
already owns the parse and already defines what "bright" means with its
own threshold. One place for one answer.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

#: The report's own version, so a reader can refuse what it cannot read.
SCHEMA = 1

AGREES = "agrees"
DIFFERS = "differs"
NOT_MEASURED = "not_measured"
AMBIGUOUS = "ambiguous"

#: How a computed transition is matched to a reference one. Declared per
#: reference and carried into the verdict, because the choice changes
#: what the comparison means.
MATCH_FIRST_BRIGHT = "first_bright"   #: the lowest-energy bright transition
MATCH_LOWEST = "lowest"               #: the lowest transition, bright or not
MATCH_INDEX = "index"                 #: the n-th transition, 1-based
_MATCH_RULES = (MATCH_FIRST_BRIGHT, MATCH_LOWEST, MATCH_INDEX)


@dataclass(frozen=True)
class Reference:
    """One published value, with what it takes to compare against it.

    ``tolerance`` is the reference's own, not a global constant: a
    coupled-cluster benchmark value and an experimental band maximum do
    not deserve the same window, and a single number hidden in this
    module would quietly decide both.
    """
    key: str                  #: the folder or molecule this is about
    quantity: str             #: "vertical_excitation_ev" is the only one so far
    value: float
    tolerance: float
    source: str               #: a citation a reader can follow, not a name
    match: str = MATCH_FIRST_BRIGHT
    index: int = 1            #: for MATCH_INDEX
    note: str = ""


@dataclass
class Verdict:
    key: str
    quantity: str
    verdict: str
    reference: Optional[float] = None
    computed: Optional[float] = None
    deviation: Optional[float] = None
    tolerance: Optional[float] = None
    match: str = ""
    source: str = ""
    why: str = ""


@dataclass
class Report:
    schema: int = SCHEMA
    workspace: str = ""
    verdicts: list[Verdict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def of(self, verdict: str) -> list[Verdict]:
        return [v for v in self.verdicts if v.verdict == verdict]

    def coverage(self) -> Optional[float]:
        """The share of references a comparable value came out for.

        The field reports this beside accuracy and for a good reason: an
        agent that answers two cases out of ten very accurately has not
        done better than one that answers nine adequately, and a single
        error figure hides which happened. None when there is nothing to
        cover, because 0/0 is not 0 %.
        """
        if not self.verdicts:
            return None
        measured = len(self.verdicts) - len(self.of(NOT_MEASURED))
        return measured / len(self.verdicts)

    def mean_absolute_deviation(self) -> Optional[float]:
        """Over the comparable values only, and never over the others.

        None when nothing was comparable -- an average over an empty set
        reported as 0.0 would read as perfect agreement.
        """
        devs = [v.deviation for v in self.verdicts
                if v.deviation is not None]
        return (sum(devs) / len(devs)) if devs else None

    def as_dict(self) -> dict:
        return {
            "schema": self.schema,
            "workspace": self.workspace,
            "coverage": self.coverage(),
            "mean_absolute_deviation_ev": self.mean_absolute_deviation(),
            "counts": {v: len(self.of(v)) for v in
                       (AGREES, DIFFERS, NOT_MEASURED, AMBIGUOUS)},
            "verdicts": [vars(v) for v in self.verdicts],
            "notes": list(self.notes),
        }


def load_references(path: Any) -> tuple[list[Reference], list[str]]:
    """Read a reference file; return what parsed and what did not.

    A row this cannot read is REPORTED, not skipped: a reference set
    silently shortened by three rows turns a coverage figure into a
    different measurement with the same name.
    """
    refs: list[Reference] = []
    problems: list[str] = []
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception as exc:                            # noqa: BLE001
        return [], [f"cannot read {path}: {type(exc).__name__}: {exc}"]
    rows = raw.get("references") if isinstance(raw, dict) else raw
    if not isinstance(rows, list):
        return [], [f"{path}: no 'references' list"]
    for i, row in enumerate(rows):
        if not isinstance(row, dict):
            problems.append(f"row {i}: not an object")
            continue
        try:
            match = str(row.get("match") or MATCH_FIRST_BRIGHT)
            if match not in _MATCH_RULES:
                problems.append(f"row {i}: unknown match rule {match!r}")
                continue
            refs.append(Reference(
                key=str(row["key"]),
                quantity=str(row.get("quantity")
                             or "vertical_excitation_ev"),
                value=float(row["value"]),
                tolerance=float(row["tolerance"]),
                source=str(row["source"]),
                match=match,
                index=int(row.get("index", 1)),
                note=str(row.get("note") or "")))
        except Exception as exc:                        # noqa: BLE001
            problems.append(f"row {i}: {type(exc).__name__}: {exc}")
    return refs, problems


def _computed(folder: Path, ref: Reference) -> tuple[Optional[float], str]:
    """The run's own value for *ref*, or None and the reason.

    Read through ``delfin.api.extract_excited_states``: it owns the
    parse and it already defines "bright" with a stated threshold, so
    this module does not get a second opinion on either.
    """
    try:
        from delfin.api import extract_excited_states
    except Exception as exc:                            # noqa: BLE001
        return None, f"the extractor is unavailable ({type(exc).__name__})"
    try:
        res = extract_excited_states(str(folder))
    except Exception as exc:                            # noqa: BLE001
        return None, f"the extractor raised {type(exc).__name__}"
    if getattr(res, "error", None):
        return None, f"no transition table: {res.error}"
    transitions = list(getattr(res, "transitions", None) or [])
    if not transitions:
        return None, "the output holds no transitions"

    if ref.match == MATCH_FIRST_BRIGHT:
        bright = getattr(res, "first_bright", None)
        if not bright:
            return None, (
                "no transition reaches the extractor's own bright "
                f"threshold ({getattr(res, 'bright_threshold_fosc', '?')})")
        return float(bright.get("energy_ev")), "the lowest bright transition"
    if ref.match == MATCH_LOWEST:
        energies = [float(getattr(t, "energy_ev", 0.0) or 0.0)
                    for t in transitions]
        return min(energies), "the lowest transition in the table"
    # MATCH_INDEX
    if ref.index < 1 or ref.index > len(transitions):
        return None, (f"the table has {len(transitions)} transitions and "
                      f"index {ref.index} was asked for")
    value = getattr(transitions[ref.index - 1], "energy_ev", None)
    if value is None:
        return None, f"transition {ref.index} carries no energy"
    return float(value), f"transition {ref.index} of the table"


def compare(workspace: Any, references: list[Reference]) -> Report:
    """Judge a run's folders against the references. Never raises.

    ``Reference.key`` is matched against a directory under *workspace*,
    which is how a run that computed two functionals into two folders is
    judged twice. A key no folder answers to is NOT_MEASURED with the
    key named -- a reference silently dropped would raise coverage by
    removing the case it failed.
    """
    report = Report(workspace=str(workspace))
    root = Path(workspace)
    if not root.is_dir():
        report.notes.append(f"the workspace is not here ({root})")
        for ref in references:
            report.verdicts.append(Verdict(
                key=ref.key, quantity=ref.quantity, verdict=NOT_MEASURED,
                reference=ref.value, tolerance=ref.tolerance,
                match=ref.match, source=ref.source,
                why="the workspace is not here"))
        return report

    for ref in references:
        folder = root / ref.key
        if not folder.is_dir():
            matches = sorted(p for p in root.rglob(ref.key) if p.is_dir())
            folder = matches[0] if matches else None
        if folder is None:
            report.verdicts.append(Verdict(
                key=ref.key, quantity=ref.quantity, verdict=NOT_MEASURED,
                reference=ref.value, tolerance=ref.tolerance,
                match=ref.match, source=ref.source,
                why=f"no folder named {ref.key!r} under the workspace"))
            continue
        value, why = _computed(folder, ref)
        if value is None:
            report.verdicts.append(Verdict(
                key=ref.key, quantity=ref.quantity, verdict=NOT_MEASURED,
                reference=ref.value, tolerance=ref.tolerance,
                match=ref.match, source=ref.source, why=why))
            continue
        deviation = abs(value - ref.value)
        report.verdicts.append(Verdict(
            key=ref.key, quantity=ref.quantity,
            verdict=AGREES if deviation <= ref.tolerance else DIFFERS,
            reference=ref.value, computed=value, deviation=deviation,
            tolerance=ref.tolerance, match=ref.match, source=ref.source,
            why=why))
    return report


def compare_file(workspace: Any, reference_file: Any) -> Report:
    """``compare`` with the references read from a file. Never raises."""
    refs, problems = load_references(reference_file)
    report = compare(workspace, refs) if refs else Report(
        workspace=str(workspace))
    for p in problems:
        report.notes.append(f"reference file: {p}")
    if not refs:
        report.notes.append(
            f"no usable reference in {reference_file} -- nothing was "
            "compared, which is not the same as agreement")
    return report


def render(report: Report) -> str:
    """One line per reference, what disagrees or is missing first."""
    order = {DIFFERS: 0, AMBIGUOUS: 1, NOT_MEASURED: 2, AGREES: 3}
    lines = [f"workspace: {report.workspace}"]
    for v in sorted(report.verdicts,
                    key=lambda v: (order.get(v.verdict, 4), v.key)):
        line = f"  {v.verdict:<13} {v.key} ({v.quantity})"
        if v.computed is not None:
            line += (f"\n                  reference {v.reference} "
                     f"computed {v.computed:.4g} "
                     f"deviation {v.deviation:.4g} "
                     f"tolerance {v.tolerance}")
        if v.match:
            line += f"\n                  matched by: {v.match} ({v.why})"
        elif v.why:
            line += f"\n                  ({v.why})"
        if v.source:
            line += f"\n                  source: {v.source}"
        lines.append(line)
    cov = report.coverage()
    mad = report.mean_absolute_deviation()
    lines.append("")
    lines.append(
        f"coverage {'n/a' if cov is None else f'{cov:.0%}'}"
        f"  mean absolute deviation "
        f"{'n/a' if mad is None else f'{mad:.4g} eV'}"
        f"  ({len(report.of(AGREES))} agree, {len(report.of(DIFFERS))} "
        f"differ, {len(report.of(NOT_MEASURED))} not measured)")
    for note in report.notes:
        lines.append(f"note: {note}")
    return "\n".join(lines)
