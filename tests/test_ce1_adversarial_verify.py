"""CE1 reviewer adversarial tests against the phase-2 verifier + fixture.

Red tests for defects in accept/w20_ce1_verify.py; green cases pin the
correct behaviour. The fixture runs for real (subprocess-free: its main()
is called directly with the workspace path).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(
        name, _REPO / "delfin" / "agent" / "pack" / "benchmark" / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_fixture = _load("w20_ce1_fixture", "setup/w20_ce1_ensemble_fixture.py")
_verify = _load("w20_ce1_verify", "accept/w20_ce1_verify.py")


@pytest.fixture()
def ws(tmp_path):
    assert _fixture.main(["fixture", str(tmp_path)]) == 0
    return tmp_path


def _accept(ws):
    return _verify.check(ws)


# --- red: wrong weights accepted because the exact values are never used ---
def test_wrong_weights_are_rejected_not_accepted(ws):
    d = ws / "w20_ce1_ensemble"
    # Structurally plausible (3 numbers, [0,1], sum 1, first largest) but
    # far from the exact Boltzmann weights at 350 K (lowest ~0.988,
    # pair ~0.006 each): the lowest frame's weight is wildly wrong.
    (d / "weights.txt").write_text("0.5\n0.3\n0.2\n", encoding="utf-8")
    ok, msg = _accept(ws)
    assert not ok, msg


# --- red: gibberish refusal mentioning "level" is accepted -----------------
def test_gibberish_refusal_is_rejected(ws):
    d = ws / "w20_ce1_ensemble"
    (d / "refusal.txt").write_text("lorem level ipsum\n", encoding="utf-8")
    ok, msg = _accept(ws)
    assert not ok, msg


# --- green: no artifact at all is a rejection ------------------------------
def test_no_artifact_is_rejected(ws):
    ok, msg = _accept(ws)
    assert not ok
    assert "none of the four artifacts" in msg


# --- green: a handoff/ with the wrong survivor count is rejected -----------
def test_wrong_handoff_count_is_rejected(ws):
    d = ws / "w20_ce1_ensemble"
    h = d / "handoff"
    h.mkdir()
    for i in range(3):  # 3 files; the 4.0 threshold leaves exactly 2
        (h / f"conf{i}.xyz").write_text(
            "3\nenergy: -12.34 Eh\nH 0.0 0.0 0.0\n", encoding="utf-8")
    ok, msg = _accept(ws)
    assert not ok
    assert "threshold leaves 2" in msg
