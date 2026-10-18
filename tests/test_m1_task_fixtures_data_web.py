"""M1 phase-3 fixture controls (data + web).

Same contract as tests/test_m1_task_fixtures.py: build under tmp_path via
the setup builder's build(), verdict through the accept main([name, ws])
convention (accept/l5_method_verify.py:100), correct solution passes,
wrong solution fails.
"""
import sys
from pathlib import Path

SETUP_DIR = Path(__file__).resolve().parents[1] / "delfin" / "agent" / \
    "pack" / "benchmark" / "setup"
ACCEPT_DIR = Path(__file__).resolve().parents[1] / "delfin" / "agent" / \
    "pack" / "benchmark" / "accept"


def _build(module_name: str, target: Path, variant: str) -> None:
    sys.path.insert(0, str(SETUP_DIR))
    import importlib
    try:
        mod = importlib.import_module(module_name)
        importlib.reload(mod)
        mod.build(target, variant=variant)
    finally:
        sys.path.remove(str(SETUP_DIR))


def _accept(module_name: str, target: Path, variant: str) -> bool:
    sys.path.insert(0, str(ACCEPT_DIR))
    import importlib
    try:
        mod = importlib.import_module(module_name)
        importlib.reload(mod)
        return mod.main([module_name + ".py", str(target)],
                        variant=variant) == 0
    finally:
        sys.path.remove(str(ACCEPT_DIR))


# ---------------------------------------------------------------- data

def test_data_audit_fixture_setup_and_accept(tmp_path):
    target = tmp_path / "data_audit"
    _build("w17_method_data_repo", target, variant="audit")
    assert (target / "data_notes" / "m1_brief.txt").exists()
    # wrong solution: nothing done -> accept fails
    assert _accept("w17_method_data_verify", target, "audit") is False
    # correct solution: the audit's output file, as the brief asks
    answer = target / "data_notes" / "data" / "AUDIT.txt"
    answer.parent.mkdir(parents=True, exist_ok=True)
    answer.write_text("duplicates: r1 r2\n", encoding="utf-8")
    assert _accept("w17_method_data_verify", target, "audit") is True


def test_data_convert_fixture_setup_and_accept(tmp_path):
    target = tmp_path / "data_convert"
    _build("w17_method_data_repo", target, variant="convert")
    # wrong: nothing
    assert _accept("w17_method_data_verify", target, "convert") is False
    # correct: converted.json mapping every id to its float
    out = target / "data_notes" / "data" / "converted.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text('{"a": 1234.5, "b": -0.5, "c": 10000.25, '
                   '"d": 2.7, "e": 1000.1}', encoding="utf-8")
    assert _accept("w17_method_data_verify", target, "convert") is True


# ---------------------------------------------------------------- web

def test_web_manual_fixture_setup_and_accept(tmp_path):
    target = tmp_path / "web_manual"
    _build("w17_method_web_repo", target, variant="manual")
    # wrong: nothing
    assert _accept("w17_method_web_verify", target, "manual") is False
    # correct: names --baseline, does not repeat the decoy --bessel
    answer = target / "web_notes" / "docs" / "answer.txt"
    answer.parent.mkdir(parents=True, exist_ok=True)
    answer.write_text("Use the --baseline flag from the advanced-setups "
                      "page.\n", encoding="utf-8")
    assert _accept("w17_method_web_verify", target, "manual") is True


def test_web_changelog_fixture_setup_and_accept(tmp_path):
    target = tmp_path / "web_changelog"
    _build("w17_method_web_repo", target, variant="changelog")
    # wrong: nothing
    assert _accept("w17_method_web_verify", target, "changelog") is False
    # correct: removal stated AND the strict replacement named, without
    # claiming lenient still exists
    answer = target / "web_notes" / "docs" / "answer.txt"
    answer.parent.mkdir(parents=True, exist_ok=True)
    answer.write_text("--lenient was removed in 5.0; use the strict "
                      "--exact instead.\n", encoding="utf-8")
    assert _accept("w17_method_web_verify", target, "changelog") is True
