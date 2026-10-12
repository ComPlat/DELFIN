"""Probe: the H1G phase-3 wiring test SEEN green against the patched copy.

Runs the same scenarios as tests/test_zz_h1g_wiring_red.py, but first inserts
the shadow package at .gate/px/h1g_pkg (fresh origin/main task_evidence.py +
api_client.py WITH the wiring applied, plus the branch's counter_test) at the
front of sys.path. This is the builder-side proof the operator's patch turns
the red test green: the red test itself fails against the un-wired code
(TypeError: unexpected keyword 'durable_test_runs'); here the keyword exists
and the note surfaces the sentence.
"""
import sys
from pathlib import Path

_SHADOW = str(Path(__file__).resolve().parents[1] / ".gate" / "px" / "h1g_pkg")


def _load_shadow():
    if _SHADOW not in sys.path:
        sys.path.insert(0, _SHADOW)
    # Force a fresh import of the shadow package: clear every delfin module
    # so the import machinery re-searches sys.path and picks up .gate/px/h1g_pkg.
    for mod in list(sys.modules):
        if mod == "delfin" or mod.startswith("delfin."):
            del sys.modules[mod]
    from delfin.agent import task_evidence as te  # noqa: F401
    return te


def _changes(*paths):
    return [{"path": p, "ts": 0.0, "created": True} for p in paths]


def _green_run(target):
    return {"tool": "run_tests", "command": target, "exit_code": 0,
            "status": "ok", "passed": 1, "failed": 0, "ts": 100.0}


def _claim(te, *, runs):
    return te.check_completion_claim(
        "run the tests for the fix",
        "",
        changes=_changes("delfin/agent/thing.py", "tests/test_thing.py"),
        tests=[_green_run("tests/test_thing.py::test_x")],
        durable_test_runs=runs,
        window_start=0.0,
    )


def test_shadow_imports_patched_signature():
    te = _load_shadow()
    # The signature now takes the keyword (it did not before the patch).
    assert "durable_test_runs" in te.check_completion_claim.__code__.co_varnames


def test_green_even_with_durable_runs_surfaces_sentence():
    te = _load_shadow()
    res = _claim(te, runs=[{"target": "tests/test_thing.py::test_x",
                            "red": False, "green": True}])
    assert res["verdict"] == "verified"
    assert "never seen failing" in (res.get("note") or ""), res


def test_controlled_fix_not_perturbed_by_wiring():
    te = _load_shadow()
    res = _claim(te, runs=[{"target": "tests/test_thing.py::test_x",
                            "red": True, "green": True}])
    assert res["verdict"] == "verified"
    assert "never seen failing" not in (res.get("note") or ""), res


def test_changed_never_run_test_names_the_test():
    te = _load_shadow()
    res = _claim(te, runs=[])
    assert "changed but never run" in (res.get("note") or ""), res
