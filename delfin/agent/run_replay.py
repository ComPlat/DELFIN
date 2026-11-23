"""Check a stored run manifest against what can be derived again now.

A manifest is a claim that a result can be arrived at again. This is the
check of that claim, and it is the half that makes the claim worth
making: a record nobody can test is a record nobody has to be right
about.

What replay does NOT do is re-run the model. A model call is not
reproducible and costs money; what is reproducible is the derivation —
the commit, the prompt the model read, the engine versions, the files
the run left behind and the execution contract read off them. Those are
re-derived here and compared to what the manifest recorded.

**One collector, not two.** The comparison is made by building a FRESH
manifest for the same session and workspace (``run_manifest.build``) and
diffing it field by field against the stored one. A second collector
written for this module would drift from the first, and then "replay
disagrees" would mean "the two readers disagree" instead of "the run
changed".

Three verdicts per field, never two:

    match        recorded and re-derived agree
    differs      both are known and they disagree
    unavailable  the field cannot be re-derived here, and WHY

``reproduced()`` is true only when something was compared and nothing
differs AND nothing is unavailable. An unavailable field is not a pass:
that is the same rule the execution contract uses for an unjudged flag,
and it exists because the opposite reading — absence as agreement — is
what makes a record of this kind useless.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from . import run_manifest

#: The report's own version, so a reader can refuse what it cannot read.
SCHEMA = 1

#: The fields whose agreement decides whether a run can be derived again.
#: Everything else in a manifest is informational and is reported as
#: such: ``written_at`` and the activity timestamps differ between two
#: honest readings of the same run, and comparing them would make every
#: replay fail for a reason that says nothing.
COMPARED: tuple[str, ...] = (
    "code.commit",
    "code.delfin_version",
    "model.id",
    "prompt.pack_digest",
    "prompt.role",
)

MATCH = "match"
DIFFERS = "differs"
UNAVAILABLE = "unavailable"


@dataclass
class Check:
    name: str
    verdict: str
    recorded: Any = None
    now: Any = None
    note: str = ""
    #: Whether this check standing open blocks `reproduced()`.
    #:
    #: Not every unavailable field is an unanswered question about the
    #: RUN. A flag that was unjudged when the run was recorded and is
    #: unjudged now -- `no_egress` on a workspace that logged no egress
    #: events -- is CONSISTENT: the record and the re-derivation agree
    #: that it is unknown. Nothing was claimed, so nothing failed to
    #: verify, and the contract's own completeness stays the contract's
    #: question (`certified`, `unjudged` in the manifest) rather than
    #: replay's.
    #:
    #: What does block is a value that WAS recorded and cannot be
    #: re-derived: there a claim exists and cannot be checked.
    #:
    #: Measured while writing the tests: without this distinction
    #: `reproduced()` was false for every run whose contract had one
    #: unjudged flag, which is nearly all of them -- and a predicate
    #: that is always false says nothing.
    blocking: bool = True


@dataclass
class Report:
    schema: int = SCHEMA
    manifest: str = ""
    session_id: str = ""
    checks: list[Check] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    #: The REQUIRED manifest fields the stored record never held
    #: (``run_manifest.missing``). Carried, not re-decided: whether a
    #: record is complete is the manifest's question, and whether what it
    #: holds derives again is this module's. Mixing the two would make a
    #: run with an incomplete record indistinguishable from one whose
    #: result changed, and those call for different work.
    incomplete: list[str] = field(default_factory=list)

    def of(self, verdict: str) -> list[Check]:
        return [c for c in self.checks if c.verdict == verdict]

    def differs(self) -> list[Check]:
        return self.of(DIFFERS)

    def unavailable(self) -> list[Check]:
        return self.of(UNAVAILABLE)

    def matches(self) -> list[Check]:
        return self.of(MATCH)

    def open_questions(self) -> list[Check]:
        """Unavailable fields that leave a recorded claim unchecked."""
        return [c for c in self.unavailable() if c.blocking]

    def consistent_unknowns(self) -> list[Check]:
        """Unknown when recorded and unknown now: nothing was claimed."""
        return [c for c in self.unavailable() if not c.blocking]

    def reproduced(self) -> bool:
        """Something was compared, nothing differs, no claim is unchecked.

        A field that was unknown when recorded and is unknown now does
        not block: see ``Check.blocking``.
        """
        return (bool(self.checks) and not self.differs()
                and not self.open_questions())

    def as_dict(self) -> dict:
        return {
            "schema": self.schema,
            "manifest": self.manifest,
            "session_id": self.session_id,
            "reproduced": self.reproduced(),
            "counts": {MATCH: len(self.matches()),
                       DIFFERS: len(self.differs()),
                       UNAVAILABLE: len(self.unavailable()),
                       "open_questions": len(self.open_questions())},
            "checks": [vars(c) for c in self.checks],
            "incomplete": list(self.incomplete),
            "notes": list(self.notes),
        }


def _dig(node: Any, dotted: str) -> Any:
    for part in dotted.split("."):
        if not isinstance(node, dict):
            return None
        node = node.get(part)
    return node


def _compare(report: Report, name: str, recorded: Any, now: Any,
             *, unresolved_now: set[str]) -> None:
    """One field, with the reason when it cannot be judged."""
    if recorded is None:
        # Nothing was claimed, so nothing failed to verify. Still shown,
        # because a reader must see which fields the record never held.
        report.checks.append(Check(
            name, UNAVAILABLE, recorded, now,
            "the manifest does not carry this field", blocking=False))
        return
    if now is None:
        why = ("it could not be re-derived here"
               if name in unresolved_now
               else "nothing here produced a value for it")
        report.checks.append(Check(name, UNAVAILABLE, recorded, now, why))
        return
    report.checks.append(Check(
        name, MATCH if recorded == now else DIFFERS, recorded, now))


def _engines(report: Report, stored: dict, fresh: dict,
             unresolved_now: set[str]) -> None:
    """Every engine either manifest mentions, by name.

    The union, not the stored list: an engine that appeared since the run
    is not a difference in the result, but an engine that is GONE is
    exactly the thing a reader needs to see before trusting a number
    that came out of it.
    """
    rec = (_dig(stored, "environment.engines") or {})
    cur = (_dig(fresh, "environment.engines") or {})
    for name in sorted(set(rec) | set(cur)):
        r, c = rec.get(name), cur.get(name)
        label = f"environment.engines.{name}"
        if r is None and c is None:
            report.checks.append(Check(label, MATCH, None, None,
                                       "absent then and now"))
        elif r is None:
            report.checks.append(Check(
                label, MATCH, None, (c or {}).get("version"),
                "not installed when the run was recorded; present now, "
                "which cannot change that run's result"))
        elif c is None:
            report.checks.append(Check(
                label, DIFFERS, (r or {}).get("version"), None,
                "the engine that produced this result is no longer here"))
        else:
            _compare(report, f"{label}.version",
                     (r or {}).get("version"), (c or {}).get("version"),
                     unresolved_now=unresolved_now)


def _workspace(report: Report, stored: dict, workspace: Optional[Path],
               ) -> None:
    """The digest of every evidence file the manifest listed.

    A missing workspace root makes every file UNAVAILABLE, which is the
    honest reading on a different machine. A root that is here with a
    file missing is a DIFFERENCE: the evidence a claim rests on is gone.
    """
    files = _dig(stored, "workspace.files")
    if not files:
        report.checks.append(Check(
            "workspace.files", UNAVAILABLE, None, None,
            "the manifest lists no evidence files"))
        return
    root = workspace or (
        Path(_dig(stored, "workspace.root"))
        if _dig(stored, "workspace.root") else None)
    if root is None or not root.is_dir():
        report.checks.append(Check(
            "workspace.files", UNAVAILABLE, f"{len(files)} files", None,
            f"the workspace root is not here ({root})"))
        return
    from .change_journal import sha256_file
    for entry in files:
        rel = str(entry.get("path") or "")
        label = f"workspace.{rel}"
        target = root / rel
        if not target.is_file():
            report.checks.append(Check(
                label, DIFFERS, entry.get("sha256"), None,
                "the file this result was read from is gone"))
            continue
        report.checks.append(Check(
            label,
            MATCH if sha256_file(target) == entry.get("sha256") else DIFFERS,
            entry.get("sha256"), sha256_file(target)))


def _contract(report: Report, stored: dict, fresh: dict) -> None:
    """The execution contract, flag by flag.

    Each flag is its own check because they fail for different reasons,
    and "the contract differs" would hide which one. None stays None:
    a flag that could not be judged then or now is unavailable, never a
    match with itself.
    """
    rec = _dig(stored, "contract.flags")
    cur = _dig(fresh, "contract.flags")
    if not isinstance(rec, dict):
        report.checks.append(Check("contract.flags", UNAVAILABLE, None, None,
                                   "the manifest carries no contract"))
        return
    if not isinstance(cur, dict):
        report.checks.append(Check(
            "contract.flags", UNAVAILABLE, rec, None,
            "the contract cannot be re-derived here (no workspace)"))
        return
    for flag in sorted(set(rec) | set(cur)):
        r, c = rec.get(flag), cur.get(flag)
        label = f"contract.{flag}"
        if r is None and c is None:
            report.checks.append(Check(
                label, UNAVAILABLE, r, c,
                "unjudged when recorded and unjudged now, consistently",
                blocking=False))
        elif r is None:
            report.checks.append(Check(
                label, UNAVAILABLE, r, c,
                "unjudged when recorded; judged now, so the record "
                "claimed nothing here", blocking=False))
        elif c is None:
            report.checks.append(Check(
                label, UNAVAILABLE, r, c,
                "recorded as a judgement and cannot be judged now",
                blocking=True))
        else:
            report.checks.append(Check(
                label, MATCH if r == c else DIFFERS, r, c))


def replay(stored: dict, *, workspace: Any = None,
           manifest_path: str = "") -> Report:
    """Compare a stored manifest against a fresh one for the same run.

    ``workspace`` overrides the root the manifest recorded, which is how
    a run is checked after it has been moved or copied elsewhere.
    """
    report = Report(manifest=str(manifest_path or ""))
    if not isinstance(stored, dict):
        report.notes.append("the manifest is not a record this can read")
        return report
    if stored.get("schema") != run_manifest.SCHEMA:
        report.notes.append(
            f"manifest schema {stored.get('schema')!r}, this reader knows "
            f"{run_manifest.SCHEMA} -- a field may have moved, so nothing "
            "is compared rather than compared wrongly")
        return report

    sid = str(_dig(stored, "run.session_id") or "")
    report.session_id = sid
    report.incomplete = run_manifest.missing(stored)
    root = workspace
    if root is None:
        recorded_root = _dig(stored, "workspace.root")
        root = recorded_root if recorded_root else None

    fresh = run_manifest.build(
        session_id=sid, workspace=root,
        expected_steps=_dig(stored, "contract.expected_steps"),
        termination=str(_dig(stored, "run.termination") or ""),
        settings={"model": _dig(stored, "model.id"),
                  "provider": _dig(stored, "model.provider"),
                  "temperature": _dig(stored, "model.temperature"),
                  "seed": _dig(stored, "model.seed")},
        role=str(_dig(stored, "prompt.role") or ""))
    unresolved_now = set(fresh.get("unresolved") or ())

    for name in COMPARED:
        _compare(report, name, _dig(stored, name), _dig(fresh, name),
                 unresolved_now=unresolved_now)
    _engines(report, stored, fresh, unresolved_now)
    _workspace(report, stored, Path(root) if root else None)
    _contract(report, stored, fresh)

    # The fields a replay must NOT compare, said once so a reader does
    # not mistake their absence for an oversight.
    report.notes.append(
        "not compared: written_at, the activity counters and the wall "
        "clock -- two honest readings of one run differ in those, and "
        "comparing them would fail every replay for a reason that says "
        "nothing about the result")
    return report


def replay_file(path: Any, *, workspace: Any = None) -> Report:
    """Read a manifest from disk and check it. Never raises."""
    p = Path(path)
    stored = run_manifest.read(p)
    if stored is None:
        rep = Report(manifest=str(p))
        try:
            raw = json.loads(p.read_text(encoding="utf-8"))
            rep.notes.append(
                f"manifest schema {raw.get('schema')!r} is not "
                f"{run_manifest.SCHEMA}; nothing was compared")
        except Exception as exc:                        # noqa: BLE001
            rep.notes.append(f"cannot read {p}: {type(exc).__name__}")
        return rep
    return replay(stored, workspace=workspace, manifest_path=str(p))


def render(report: Report) -> str:
    """One line per check, differences and open fields first.

    Ordered that way on purpose: a reader scanning the top must see what
    is wrong or unknown without reading the matches.
    """
    lines = [f"manifest: {report.manifest or '(inline)'}"]
    if report.session_id:
        lines.append(f"session:  {report.session_id}")
    order = {DIFFERS: 0, UNAVAILABLE: 1, MATCH: 2}
    for c in sorted(report.checks, key=lambda c: (order.get(c.verdict, 3),
                                                  c.name)):
        line = f"  {c.verdict:<12} {c.name}"
        if c.verdict == DIFFERS:
            line += f"\n                 recorded: {c.recorded!r}\n" \
                    f"                 now:      {c.now!r}"
        if c.note:
            line += f"\n                 ({c.note})"
        lines.append(line)
    lines.append("")
    lines.append(f"{len(report.matches())} match, "
                 f"{len(report.differs())} differ, "
                 f"{len(report.open_questions())} recorded and unverifiable, "
                 f"{len(report.consistent_unknowns())} unknown then and now")
    lines.append("REPRODUCED" if report.reproduced()
                 else "NOT REPRODUCED -- see above")
    if report.incomplete:
        lines.append(
            "the manifest itself was incomplete when written, missing: "
            + ", ".join(report.incomplete)
            + " (a gap in the record, not a change in the result)")
    for note in report.notes:
        lines.append(f"note: {note}")
    return "\n".join(lines)
