"""A marker that fires once per process must not cross tests.

The suite runs every test in one process, so a module-level "already
said this" marker is shared by all of them. The first test to trip one
silences it for the rest, and a test that asserts the notice passes on
its own and fails behind a neighbour -- which reads as flakiness and
sends whoever is on duty looking in the wrong place.

Two were measured on 2026-11-19. ``_ENDPOINT_TTFT_S`` is fed by every
completed stream round and decides the watcher's wait: under three
samples the patched default, from three a computed value with a 5 s
floor. Stubbed rounds leaving 0.0006 s behind were enough, because the
floor decides, and a classifier test failed on CI as ``assert 0 == 1``
in a shard of a branch that does not touch the classifier.
``_MCP_WITHHELD_REPORTED`` records a withheld tool description once per
tool: running the test that asserts that event twice in one process
fails the second time.

The fixture in conftest empties the registered markers before each
test. This file holds it to its job and, more importantly, makes the
registry fail when the package grows a marker nobody registered -- the
defect is invisible until it bites, so the list cannot be kept by
memory.
"""

from __future__ import annotations

import ast
import importlib
import pathlib
import re

import pytest

from conftest import _REPORT_ONCE_MARKERS


def test_the_registry_is_not_empty():
    """A guard on the guard: an import that quietly produced an empty
    tuple would make every assertion below vacuous."""
    assert len(_REPORT_ONCE_MARKERS) >= 18


@pytest.mark.parametrize("module_name,attr", _REPORT_ONCE_MARKERS)
def test_each_registered_marker_exists_and_can_be_emptied(module_name, attr):
    """A renamed or deleted marker must fail here rather than leave the
    registry quietly pointing at nothing."""
    module = importlib.import_module(module_name)
    value = getattr(module, attr, None)

    assert value is not None, f"{module_name}.{attr} is gone or renamed"
    assert hasattr(value, "clear"), (
        f"{module_name}.{attr} is a {type(value).__name__}, which the "
        f"fixture cannot empty in place")


@pytest.mark.parametrize("module_name,attr", _REPORT_ONCE_MARKERS)
def test_every_registered_marker_starts_this_test_empty(module_name, attr):
    """What the fixture is for, asserted directly."""
    value = getattr(importlib.import_module(module_name), attr, None)
    if value is None:
        pytest.skip(f"{module_name} not importable here")

    assert len(value) == 0, (
        f"{module_name}.{attr} carried {len(value)} entries into this "
        f"test; the fixture did not empty it")


# --- the isolation, demonstrated rather than asserted ----------------------
#
# The two tests below run in file order. The first fills a marker, as any
# test driving a real round does; the second asserts it arrives empty. If
# the fixture were removed, the second would fail -- which is the control
# for this whole file.


def test_one_test_fills_a_marker():
    from delfin.agent import api_client as ac

    ac._MCP_WITHHELD_REPORTED.add("chem/lookup")
    ac._ENDPOINT_TTFT_S.append(0.0006)

    assert ac._MCP_WITHHELD_REPORTED, "the marker is filled on purpose here"


def test_the_next_test_does_not_inherit_it():
    from delfin.agent import api_client as ac

    assert "chem/lookup" not in ac._MCP_WITHHELD_REPORTED, (
        "the previous test's marker crossed into this one; a notice it "
        "silenced would be missing from whatever asserts it next")
    assert len(ac._ENDPOINT_TTFT_S) == 0, (
        "and the latency window is its own, so the watcher's wait is the "
        "one this test sets rather than a neighbour's measurements")


# --- the registry must keep up with the package ---------------------------


_ONCE_NAME = re.compile(
    r"^_?[A-Z0-9_]*(SAID|REPORTED|NOTED|NOTES|SENT|TRIED|SEEN|WARNED|"
    r"ANNOUNCED|TOLD)[A-Z0-9_]*$")

#: Names the pattern catches that are NOT once-per-process markers, with
#: the reason. Being here is a decision; being absent is an omission the
#: test reports.
_NOT_A_MARKER: dict[str, str] = {
    "_PACKAGES_TRIED": "qm_health install bookkeeping, not a notice; it "
                       "records what was attempted so a failed install is "
                       "not retried in a loop",
    "_TRIED": "the same shape for model fallbacks inside one turn",
    "_GXTB_VERIFIED": "a binary check whose answer does not change while "
                      "the process runs",
    "_KR_TRIED": "keyring names already probed. Emptying it per test would "
                 "send every test back to the OS keyring, which is a D-Bus "
                 "service that can block -- the cost of the reset is worse "
                 "than what it protects",
    "_NOTES": "subagents routes delegate messages by sa_id through this; "
              "the name matches the pattern, the contents are message text "
              "rather than a notice already given",
    "MODEL_NOTES": "a constant table of prompt text per model knob",
    "VIEWER_REPRESENTATIONS": "a constant set of 3Dmol representation names",
}


def _dotted(path: pathlib.Path, root: pathlib.Path) -> str:
    """``delfin/agent/api_client.py`` -> ``delfin.agent.api_client``.

    The registry in conftest is keyed by the dotted module name, and the
    guard below compares PAIRS against it, so this side must speak the
    same language. A path key made the comparison fall back to bare
    attribute names, which is the defect this function exists to remove.
    """
    rel = path.relative_to(root.parent).with_suffix("")
    return ".".join(rel.parts)


def _declared_markers() -> dict[str, set[str]]:
    """{dotted module: {name}} for once-shaped module-level containers."""
    found: dict[str, set[str]] = {}
    root = pathlib.Path(__file__).resolve().parent.parent / "delfin"
    for path in sorted(root.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in tree.body:
            targets = (node.targets if isinstance(node, ast.Assign)
                       else [node.target] if isinstance(node, ast.AnnAssign)
                       else [])
            names = [t.id for t in targets if isinstance(t, ast.Name)]
            if not names:
                continue
            value = node.value
            container = False
            if isinstance(value, (ast.Set, ast.List, ast.Dict)):
                container = True
            elif isinstance(value, ast.Call):
                fn = (getattr(value.func, "attr", None)
                      or getattr(value.func, "id", ""))
                container = fn in ("set", "list", "dict", "deque",
                                   "defaultdict", "OrderedDict")
            if not container:
                continue
            for name in names:
                if _ONCE_NAME.match(name) and name not in _NOT_A_MARKER:
                    found.setdefault(_dotted(path, root), set()).add(name)
    return found


def test_the_package_grew_no_marker_the_registry_does_not_know():
    """The list in conftest cannot be kept by memory: a new "already
    said this" set is invisible until a test two shards away fails for
    no reason of its own. This finds it at the commit that adds it.

    A name the pattern catches that is not a notice belongs in
    _NOT_A_MARKER with its reason -- a decision on the record, rather
    than a pattern loosened until it matches nothing.
    """
    # PAIRS, not bare names. Reading only the attribute asked "is this
    # name registered somewhere" instead of "is THIS module's marker
    # registered": _ANNOUNCED is registered for three modules
    # (container_exec, exchange, workspace_trust), so a fourth module
    # declaring _ANNOUNCED = set() satisfied the guard and was never
    # emptied between tests -- the exact leak this file exists to stop.
    registered = set(_REPORT_ONCE_MARKERS)
    missing: list[str] = []
    for module, names in _declared_markers().items():
        for name in sorted(names):
            if (module, name) not in registered:
                missing.append(f"{module}: {name}")

    assert not missing, (
        "once-per-process markers that no test isolates:\n  "
        + "\n  ".join(missing)
        + "\n\nAdd each to conftest._REPORT_ONCE_MARKERS, or to "
          "_NOT_A_MARKER in this file with the reason it is not one.")
