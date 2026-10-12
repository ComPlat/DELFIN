"""Green-flip probe: the H1G phase-3 wiring red test, SEEN green on the patch.

Runs the EXACT test module `tests/test_zz_h1g_wiring_red.py` — which feeds
REAL tool_trace-shaped entries `{tool:"bash", input:{command}, output, ok, ts}`
through `counter_test.durable_test_runs` and into
`task_evidence.check_completion_claim(durable_test_runs=...)` — but with the
patched shadow package (`.gate/px/h1g_pkg`: origin/main task_evidence.py +
api_client.py WITH the wiring applied, plus the F1-F5-corrected
counter_test.py) at the front of sys.path, so `from delfin.agent import ...`
resolves to the patched copies. Every red-test function that must go green is
run from the shadow; the green only surfaces if the operator's patch delivers
the `durable_test_runs=` keyword AND the F1-F5 fix makes a real trace-run map
onto the changed counter-test. This is the "seen going green once" evidence
for the operator build: the red test goes RED now (TypeError on the
un-wired signature), GREEN here against the patched+i corrected shadow.
"""
import importlib.util
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_SHADOW = str(_ROOT / ".gate" / "px" / "h1g_pkg")
_RED_TEST = _ROOT / "tests" / "test_zz_h1g_wiring_red.py"


def _load_red_under_shadow():
    """Insert the shadow at the front of sys.path, clear delfin from
    sys.modules so the import machinery re-searches and picks up the patched
    shadow copy, then import the red test module through importlib."""
    if _SHADOW not in sys.path:
        sys.path.insert(0, _SHADOW)
    for mod in list(sys.modules):
        if mod == "delfin" or mod.startswith("delfin."):
            del sys.modules[mod]
    spec = importlib.util.spec_from_file_location(
        "test_zz_h1g_wiring_red_shadow", str(_RED_TEST))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_green_flip_runs_the_red_test_against_the_shadow():
    """The same scenarios the red test asserts RED on the un-wired code pass
    GREEN on the patched shadow — the operator-applied wiring + F1-F5 fix."""
    mod = _load_red_under_shadow()
    # The red test defines seven ``test_*`` functions (grep confirms 7); the
    # green flip must flip EVERY one of them, not a subset — any of the 7
    # that fails raises here, and `ran == 7` confirms none were skipped.
    ran = 0
    for name in dir(mod):
        if name.startswith("test_") and callable(getattr(mod, name)):
            fn = getattr(mod, name)
            if getattr(fn, "__module__", None) == mod.__name__:
                fn()
                ran += 1
    assert ran == 7, (
        f"expected all 7 red scenarios to green-flip against the shadow, "
        f"got {ran}"
    )
