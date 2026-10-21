"""M1 phase-3 fixture controls.

Builds each fixture under tmp_path through the setup builder functions
(NEVER run setup/ or accept/ scripts directly — operator rule after s26's
stop), then applies a scripted solution and checks the accept verdict
functions through the gate. Each fixture: correct solution passes, wrong
solution fails.
"""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SETUP_DIR = Path(__file__).resolve().parents[1] / "delfin" / "agent" / \
    "pack" / "benchmark" / "setup"
ACCEPT_DIR = Path(__file__).resolve().parents[1] / "delfin" / "agent" / \
    "pack" / "benchmark" / "accept"


def _run_builder(module_name: str, target: Path, variant: str | None = None):
    sys.path.insert(0, str(SETUP_DIR))
    import importlib
    try:
        mod = importlib.import_module(module_name)
        importlib.reload(mod)
        if variant is None:
            mod.build(target)
        else:
            mod.build(target, variant=variant)
    finally:
        sys.path.remove(str(SETUP_DIR))


def _scripted_correct_and_wrong(module_name: str, target: Path) -> bool:
    """Call the accept module through its runner contract.

    The convention (accept/l5_method_verify.py:100 runs main(sys.argv)) is
    that argv carries the script name first: main([name, workspace]).
    """
    sys.path.insert(0, str(ACCEPT_DIR))
    import importlib
    try:
        mod = importlib.import_module(module_name)
        importlib.reload(mod)
        return mod.main([module_name + ".py", str(target)]) == 0
    finally:
        sys.path.remove(str(ACCEPT_DIR))


# ---------------------------------------------------------------- code

def test_code_red_fixture_setup_and_accept(tmp_path):
    target = tmp_path / "code_red"
    _run_builder("w17_method_code_repo", target, variant="red")
    assert (target / "parse_tools" / "m1_brief.txt").exists()
    # wrong solution: nothing changed -> suite still red -> accept fails
    assert _scripted_correct_and_wrong("w17_method_code_verify",
                                       target) is False
    # correct solution: fix the parser, suite turns green
    parse = target / "parse_tools" / "parse_x.py"
    src = parse.read_text()
    assert 'raise ValueError("not implemented")' in src
    src = src.replace('raise ValueError("not implemented")',
                      'return float(cleaned.replace(",", "."))')
    parse.write_text(src)
    # The verdict channel is the accept script (the runner's contract);
    # it runs the suite itself with the workspace on PYTHONPATH.
    assert _scripted_correct_and_wrong("w17_method_code_verify",
                                       target) is True


def test_code_feature_fixture_setup_and_accept(tmp_path):
    target = tmp_path / "code_feature"
    _run_builder("w17_method_code_repo", target, variant="feature")
    ok = _scripted_correct_and_wrong("w17_method_code_verify", target)
    assert ok is False
    src = target / "parse_tools" / "parse_x.py"
    src.write_text(src.read_text() + '\n\ndef slugify(s: str) -> str:\n'
                   '    return s.lower().replace(" ", "-")\n')
    assert _scripted_correct_and_wrong("w17_method_code_verify",
                                       target) is True
