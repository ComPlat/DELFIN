"""A capability the product never reaches is not a capability.

The dominant defect shape in this repository is a mechanism whose parts are
each tested and which is wired to nothing. The commit trailer was one: a
line builder, a quote-aware and heredoc-aware ``--trailer`` rewriter, a
setting that defaults to on, and a passing test file of its own -- and
nothing in ``delfin/`` ever called ``set_commit_trailer_model``, so no
commit any agent made carried the line. The parts were covered. The whole
was unreachable, for months, in silence.

This asks the question for the whole tree. For every module-level callable
in ``delfin/agent/``, is its name mentioned anywhere else in the product?

Input: the source tree. Output: two lists, and a baseline file to compare
them against (``tests/fixtures/unmounted_capabilities.txt``).

Semantics: the baseline only shrinks. A name not in it fails the test, so a
newly orphaned mechanism is caught on the commit that orphans it; a name
that is connected or deleted has to be removed from the file, so the list
cannot quietly keep stale entries either.

Two tiers, because they are different findings:

* ``nowhere``   -- nothing in ``delfin/`` or ``tests/`` names it. Dead code,
                   or a capability that was never connected to anything.
* ``test-only`` -- only the suite reaches it. The capability exists and is
                   proven, and nothing in the running product uses it.

What this measure is NOT: a defect count. A function reached through
``getattr``, a documented entry point for outside users, and a module
loaded by filename from data all look identical to it. It says where to
look. That is why the baseline is a file of names a person has seen rather
than a number.

Four ways of asking this were measured and three were wrong:

* scanning ``tests/`` without excluding the baseline file makes every name
  in it read as reached, because the file lists them: the ``nowhere`` tier
  came out empty and all 90 looked merely test-only. The instrument was
  measuring itself;

* counting identifier tokens only reports false orphans, because
  ``doctor.py`` dispatches its checks from a table of
  ``("credentials", "_check_credentials")`` pairs and the RPC surface is
  listed as tool-name strings -- so string literals count as references;
* reading Python files only reports the benchmark pack as orphaned, because
  ``pack/benchmark/tasks.yaml`` names its setup and accept modules by
  filename -- so every text file under ``delfin/`` is read;
* asking "does another FILE name it" calls every correctly-private helper
  an orphan: 2,968 of them. For a private name the question is whether
  anything names it at all.
"""
from __future__ import annotations

import ast
import io
import re
import tokenize
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BASELINE = Path(__file__).resolve().parent / "fixtures" / "unmounted_capabilities.txt"

#: Not source text. Reading these as text would spray random identifiers
#: into the reference count and hide real orphans.
_BINARY = frozenset({
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".ico", ".woff", ".woff2",
    ".zip", ".gz", ".so", ".pyc", ".xlsx", ".docx", ".pptx", ".odt", ".ods",
})

_WORD = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def _definitions() -> dict[str, tuple[str, int]]:
    """``{name: (file, line)}`` for every module-level callable in the agent.

    Module level only: a method is reached through its class and asking the
    same question of it would need the class's own reachability first.
    """
    out: dict[str, tuple[str, int]] = {}
    for path in sorted((ROOT / "delfin" / "agent").rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                    and not node.name.startswith("__"):
                out.setdefault(node.name,
                               (str(path.relative_to(ROOT)), node.lineno))
    return out


def _mentions(tree_root: Path, defs: dict, skip_own_definition: bool) -> dict:
    """How often each name is mentioned under ``tree_root``.

    Python files are tokenised, so a name in a comment does not count while
    a name in a string literal does -- that is name-based dispatch, which is
    a real reference. Everything else readable as text is scanned for words,
    which is how a module named from YAML is found. With
    ``skip_own_definition`` the ``def`` line itself does not count as a
    mention of the thing it defines.
    """
    hits = {name: 0 for name in defs}
    for path in sorted(tree_root.rglob("*")):
        if not path.is_file() or path.suffix.lower() in _BINARY:
            continue
        # The baseline lists every name it holds, and it lives under
        # tests/. Scanning it made all 90 of them read as "reached by the
        # suite" and the "nowhere" tier came out empty -- the instrument
        # measuring itself. Measured: 31 of the 90 are named nowhere.
        if path == BASELINE:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        rel = str(path.relative_to(ROOT))
        if path.suffix == ".py":
            try:
                for tok in tokenize.generate_tokens(io.StringIO(text).readline):
                    if tok.type == tokenize.NAME and tok.string in hits:
                        name = tok.string
                        if (skip_own_definition and rel == defs[name][0]
                                and tok.start[0] == defs[name][1]):
                            continue
                        hits[name] += 1
                    elif tok.type == tokenize.STRING:
                        for word in _WORD.findall(tok.string):
                            if word in hits:
                                hits[word] += 1
                continue
            except (tokenize.TokenError, IndentationError, SyntaxError):
                pass        # unparseable: fall through to the word scan
        for word in _WORD.findall(text):
            if word in hits:
                hits[word] += 1
    return hits


@pytest.fixture(scope="module")
def measured():
    defs = _definitions()
    assert len(defs) > 3000, (
        f"only {len(defs)} callables found -- the scan is looking at the "
        "wrong tree, and an empty scan reports no orphans")
    in_product = _mentions(ROOT / "delfin", defs, True)
    unmounted = sorted(n for n in defs if in_product[n] == 0)
    in_suite = _mentions(ROOT / "tests", defs, False)
    return defs, unmounted, in_suite


def _baseline() -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in BASELINE.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tier, name, _where = line.split("\t")
        rows[name] = tier
    return rows


def test_the_baseline_file_is_readable():
    rows = _baseline()
    assert rows, "the baseline is empty, so the test below asserts nothing"
    assert set(rows.values()) <= {"nowhere", "test-only"}


def test_no_capability_is_newly_unreachable(measured):
    """The one that catches the next commit-trailer.

    A mechanism that stops being called -- or arrives without ever being
    called -- is named here, on the commit that does it.
    """
    defs, unmounted, _in_suite = measured
    known = _baseline()
    new = [n for n in unmounted if n not in known]
    assert not new, (
        "nothing in delfin/ names these, so they cannot run: "
        + ", ".join(f"{n} ({defs[n][0]}:{defs[n][1]})" for n in sorted(new))
        + ". Connect it, delete it, or -- if it is reached in a way this "
          "scan cannot see -- add it to "
          "tests/fixtures/unmounted_capabilities.txt with the reason.")


def test_the_baseline_holds_no_name_that_is_reached_again(measured):
    """The other direction: the list may not keep stale entries.

    Without this the file would grow into a permanent excuse, and the count
    would stop meaning anything.
    """
    _defs, unmounted, _in_suite = measured
    gone = sorted(set(_baseline()) - set(unmounted))
    assert not gone, (
        "these are reached from delfin/ now; take them out of "
        f"tests/fixtures/unmounted_capabilities.txt: {gone}")


def test_the_tier_of_each_name_is_right(measured):
    """``nowhere`` and ``test-only`` are different findings, and the file
    has to say which: a capability the suite proves and nobody uses is not
    the same thing as dead code."""
    _defs, unmounted, in_suite = measured
    known = _baseline()
    wrong = []
    for name in unmounted:
        tier = "test-only" if in_suite.get(name) else "nowhere"
        if known.get(name) != tier:
            wrong.append(f"{name}: file says {known.get(name)}, measured {tier}")
    assert not wrong, wrong


def test_the_count_only_goes_down(measured):
    """A ratchet. 90 on 2026-11-20, of which 31 are named nowhere at all.

    Paid down, not raised: a change that needs the ceiling lifted is a
    change that leaves another capability unreachable.
    """
    _defs, unmounted, _in_suite = measured
    assert len(unmounted) <= 90, (
        f"{len(unmounted)} capabilities the product does not reach, against "
        "a ceiling of 90")
