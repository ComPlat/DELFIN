"""Package A2W reviewer (nacht-s12): red tests for the operator's findings
F1 (default OFF, fail closed) and F2 (pip shape requires --only-binary=:all:).

F1 (operator, 2026-10-17): letting commands skip a check is the person's
decision -- an opt-in, not a default. Missing payload, missing block,
missing key, None, non-bool (including truthy strings and 1) must all read
as OFF. Today the module maps missing/None/non-bool -> True (fail open).

F2: --require-hashes alone only pins what the AGENT itself wrote into the
requirements file; without --only-binary=:all: pip still runs setup.py
from an sdist at install time. So the rated pip shape requires BOTH flags.
"""

import pytest

import delfin.agent.workspace_selfservice as ws


# The xfails of the review round are gone: the fix (c83d35f9) landed and
# all 8 formerly-red cases pass plain. The tests pin the fail-closed
# contract from here on.

# --- F1: fail closed ---------------------------------------------------------

def test_f1_defaults_off():
    assert ws.selfservice_defaults() == {"enabled": False}


def test_f1_missing_payload_is_off():
    assert ws.settings_from_payload({}) == {"enabled": False}
    assert ws.settings_from_payload(None) == {"enabled": False}


def test_f1_missing_key_is_off():
    assert ws.settings_from_payload({"agent": {}}) == {"enabled": False}


def test_f1_non_bool_is_off():
    for bad in (None, "yes", "true", 1, 0, [], {}, "False"):
        assert ws.settings_from_payload({"agent": {"workspace_selfservice":
                                                   {"enabled": bad}}}) == {
            "enabled": False}, bad


def test_f1_wrong_type_block_is_off():
    assert ws.settings_from_payload({"agent": {"workspace_selfservice": "on"}}
                                    ) == {"enabled": False}


def test_f1_explicit_true_still_on():
    # The opt-in path must keep working: the person sets enabled: true.
    assert ws.settings_from_payload(
        {"agent": {"workspace_selfservice": {"enabled": True}}}
    ) == {"enabled": True}


# --- F2: pip shape needs --only-binary=:all: ----------------------------------

def _venv(tmp_path):
    venv = tmp_path / ".venv"
    (venv / "bin").mkdir(parents=True)
    (tmp_path / "requirements.txt").write_text("numpy==1.26.4\n")
    return venv


def test_f2_pip_without_only_binary_is_none(tmp_path):
    venv = _venv(tmp_path)
    cmd = f"{venv}/bin/pip install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_f2_pip_without_only_binary_interpreter_form_is_none(tmp_path):
    venv = _venv(tmp_path)
    cmd = f"{venv}/bin/python -m pip install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_f2_pip_with_only_binary_is_auto(tmp_path):
    venv = _venv(tmp_path)
    cmd = (f"{venv}/bin/pip install -r requirements.txt --require-hashes "
           "--only-binary=:all:")
    assert ws.classify(cmd, str(tmp_path)) == ws.AUTO


def test_f2_attached_only_binary_is_auto(tmp_path):
    venv = _venv(tmp_path)
    cmd = (f"{venv}/bin/pip install --only-binary=:all: -r requirements.txt "
           "--require-hashes")
    assert ws.classify(cmd, str(tmp_path)) == ws.AUTO


def test_f2_only_binary_wrong_value_is_none(tmp_path):
    # A narrowed --only-binary for specific packages does NOT ban sdists.
    venv = _venv(tmp_path)
    cmd = (f"{venv}/bin/pip install -r requirements.txt --require-hashes "
           "--only-binary=numpy")
    assert ws.classify(cmd, str(tmp_path)) is None
