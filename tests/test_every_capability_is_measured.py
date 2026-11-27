"""A capability nothing measures cannot be shown to improve.

The benchmark pack holds 280 tasks and three registries hold 131
capabilities, and until this file existed nothing connected them: no
count could say whether ``tool:read_document`` was measured, and the
honest answer for most of them was "nobody knows". An untested capability
is not merely unproven -- it is the one place a regression enters without
anyone noticing, because no task fails when it breaks.

WHAT IS ENUMERATED, AND FROM WHERE. The capability universe is never
written down here. It is read from the registries the product itself
serves from:

  * ``tool:<name>``     -- the advertised tool catalogue, the same list
                           the model's schema is built from
  * ``preset:<name>``   -- ``subagent_type_names()``, the live preset
                           registry (built-ins plus the pack's .md files)
  * ``workflow:<name>`` -- the workflow list the platform offers

All three are imported from the product, so a tool added tomorrow appears
in the census the moment it is advertised -- nobody has to remember to
add it. A hand-written list of capabilities would be the very defect this
file exists to catch, one level up: it would describe the product as it
was when the list was typed.

WHAT A TASK CLAIMS. Each task's ``capability:`` tags, read through the
real loader (``load_tasks``), so a tag survives only if the loader is
able to carry it -- a key the loader drops is data nothing can measure.
A task with no tag measures nothing and covers nothing.

THE RATCHET. A capability with no task must be named in
``tests/fixtures/unmeasured_capabilities.txt`` with a reason. The list
may only SHRINK, in both directions:

  * a capability with no task and no entry        -> FAIL (new hole)
  * an entry that now HAS a task                  -> FAIL (stale entry)

The second half is what stops the file from becoming a place to put
things and forget them: covering a capability forces its removal, so the
number in the file is a real count of what is unmeasured, not a count of
what was unmeasured once.

WHY THE BASELINE IS BIG AND STARTS AT 122. The obvious shortcut --
inferring a tag from a capability's name appearing in the task prompt --
was tried and rejected on evidence: it covered 7 of 128 capabilities, and
its tags were accidents of phrasing (an OCCUPIER task scored only
``workflow:OCCUPIER`` although it plainly reads and writes files). The
counts a reader can check today say the same thing without that
measurement: of the pack's 280 tasks, 11 carry a tag and 269 do not, and
those 11 name 9 of the 128 capabilities. Inferring a tag is not the same
as proving a task measures the capability, so a tag is now an author's
claim, checked for existence and accounted for, never inferred.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks

ROOT = Path(__file__).resolve().parents[1]
BASELINE = Path(__file__).resolve().parent / "fixtures" / "unmeasured_capabilities.txt"

# What a reason line must carry to be a reason. A bare "TODO" or "n/a"
# silences the census while explaining nothing, which is the same defect
# as a passing assertion that cannot fail.
_MIN_REASON_CHARS = 20


def _registry_capabilities() -> set[str]:
    """Every capability the product offers, from its own registries."""
    from delfin.agent import tool_catalog
    from delfin.agent.subagents import subagent_type_names
    from delfin.tools.workflows import list_workflows

    caps = {f"tool:{t['function']['name']}"
            for t in tool_catalog._DOC_TOOLS_OPENAI}
    caps |= {f"preset:{n}" for n in subagent_type_names()}
    caps |= {f"workflow:{getattr(w, 'name', str(w))}" for w in list_workflows()}
    return caps


def _baseline() -> dict[str, str]:
    rows: dict[str, str] = {}
    # A missing file reads as EMPTY, not as a crash. The interesting
    # statement is "121 capabilities have no task and no entry in the
    # baseline", which names what to do; a FileNotFoundError from a
    # helper would hide that behind a path and make the red look like a
    # setup accident rather than the finding.
    if not BASELINE.is_file():
        return rows
    for line in BASELINE.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cap, _sep, reason = line.partition("\t")
        rows[cap.strip()] = reason.strip()
    return rows


@pytest.fixture(scope="module")
def measured():
    caps = _registry_capabilities()
    # A scan that finds nothing must not report full coverage. The three
    # registries served 128 capabilities when this was written; a
    # collapse to a handful means the registries moved and the census is
    # now reading the wrong thing -- the failure mode where an empty scan
    # reports every capability covered.
    assert len(caps) > 100, (
        f"only {len(caps)} capabilities found from the product's "
        "registries -- the census is reading the wrong surface, and an "
        "empty universe reports nothing unmeasured")
    tasks = load_tasks()
    assert len(tasks) > 200, f"only {len(tasks)} benchmark tasks loaded"
    covered: set[str] = set()
    for task in tasks:
        covered |= set(getattr(task, "capability", ()) or ())
    return caps, tasks, covered


def test_the_registries_are_actually_read(measured):
    """The control for every other test here: a scan that looked at the
    wrong place would report no holes and pass forever."""
    caps, tasks, covered = measured
    assert any(c.startswith("tool:") for c in caps)
    assert any(c.startswith("preset:") for c in caps)
    assert any(c.startswith("workflow:") for c in caps)
    assert covered, (
        "no task carries a capability tag -- the census is measuring an "
        "empty claim set, so every assertion below it would pass")


def test_a_tag_names_a_capability_that_exists(measured):
    """A tag that names nothing is a typo that covers nothing, and the
    census would read it as coverage."""
    caps, tasks, _covered = measured
    unknown: list[str] = []
    for task in tasks:
        for tag in getattr(task, "capability", ()) or ():
            if tag not in caps:
                unknown.append(f"{task.id}: {tag}")
    assert not unknown, (
        "these tasks tag a capability that no registry offers, so the "
        "tag measures nothing:\n  " + "\n  ".join(sorted(unknown)))


def test_no_capability_is_unmeasured_without_a_reason(measured):
    """The ratchet: a hole must be named in the baseline, with a reason.

    Fails on a capability that has no task and no entry. This is the
    assertion that is red on the base commit, by construction: the base
    tags no capability at all -- it has no ``capability`` field, so it
    covers 0 of 128 -- and the baseline did not exist there.
    """
    caps, _tasks, covered = measured
    known = _baseline()
    missing = sorted(c for c in caps if c not in covered and c not in known)
    assert not missing, (
        f"{len(missing)} capabilit"
        f"{'y has' if len(missing) == 1 else 'ies have'} no benchmark task "
        "and no entry in tests/fixtures/unmeasured_capabilities.txt. Either "
        "write a task that measures it, or list it there with the reason it "
        "cannot be measured:\n  " + "\n  ".join(missing))


def test_a_covered_capability_is_not_still_listed_as_unmeasured(measured):
    """The other direction: the list may only shrink.

    Without this, the baseline becomes a place to put a capability and
    forget it -- the number would stop describing what is unmeasured and
    start describing what was unmeasured once.
    """
    _caps, tasks, covered = measured
    known = _baseline()
    stale = sorted(c for c in known if c in covered)
    assert not stale, (
        "these capabilities now have a benchmark task but are still listed "
        "as unmeasured; remove them from "
        "tests/fixtures/unmeasured_capabilities.txt:\n  "
        + "\n  ".join(stale))


def test_every_baseline_entry_carries_a_reason(measured):
    """A baseline of bare names answers "what is unmeasured" but not
    "why", which is the half that decides whether the hole is
    acceptable."""
    caps, _tasks, _covered = measured
    known = _baseline()
    thin = sorted(c for c, reason in known.items()
                  if len(reason) < _MIN_REASON_CHARS)
    assert not thin, (
        "these baseline entries carry no usable reason (fewer than "
        f"{_MIN_REASON_CHARS} characters):\n  " + "\n  ".join(thin)
        + "\nEach line is '<capability>\\t<why it has no task yet>'.")
    unknown = sorted(c for c in known if c not in caps)
    assert not unknown, (
        "the baseline names capabilities no registry offers any more; the "
        "list is measuring something that does not exist:\n  "
        + "\n  ".join(unknown))


def test_the_numbers_in_the_header_are_the_numbers_measured(measured):
    """The file's own prose is a claim, and a stale claim re-teaches the
    wrong thing.

    The counts above were once written as 265 tasks, a baseline starting at
    121 and "the base covers 7 of 128". Measured on the tree that ships them:
    270 tasks, 119 entries, and a base that tags no capability at all. The
    assertions here derive from the live registries, so they cannot go stale
    the way the prose did -- which is why the prose names them instead of the
    memory of a measurement nobody can repeat.
    """
    caps, tasks, covered = measured
    known = _baseline()
    source = Path(__file__).read_text(encoding="utf-8")
    tagged = sum(1 for t in tasks if getattr(t, "capability", ()))
    # Prose wraps, so line breaks are not information: compare on collapsed
    # whitespace. (the phrase the header uses, the number beside it) -- one
    # entry per count the header claims, either order.
    flat = " ".join(source.split())
    claimed = [
        ("tasks", len(tasks)), ("capabilities", len(caps)),
        ("STARTS AT", len(known)), ("carry a tag", tagged),
        ("do not", len(tasks) - tagged), ("name", len(covered)),
    ]
    absent = [f"{phrase} {value}" for phrase, value in claimed
              if f"{phrase} {value}" not in flat
              and f"{value} {phrase}" not in flat]
    assert not absent, (
        "the header states numbers that the live registries do not produce, "
        "or states them in prose that no longer says what it measured: "
        f"missing {absent}")

