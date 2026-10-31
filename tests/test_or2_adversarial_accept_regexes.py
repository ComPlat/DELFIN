"""Adversarial tests for OR2 phase 1 (reviewer nacht-s26).

The acceptance scripts decide on regexes over the agent's result text.
Allow-list / bypass attacks: a bare potential reported as valid must be
rejected by the trap acceptance (w20-or2-05) no matter how it is
formatted, and the redox acceptance (w20-or2-04) must not accept a
band-fitting number that is not an oxidation potential at all.
"""

from __future__ import annotations

import importlib
import pathlib
import sys

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark" / "accept"
sys.path.insert(0, str(_PACK))

_acc04 = importlib.import_module("w20_or2_04_redox_potential")
_acc05 = importlib.import_module("w20_or2_05_redox_mixed_trap")


def _ws(tmp_path, text):
    d = tmp_path / "or2"
    (d / "redox").mkdir(parents=True)
    (d / "redox" / "result.txt").write_text(text, encoding="utf-8")
    return tmp_path


def _ws05(tmp_path, text):
    d = tmp_path / "or2" / "redox_mixed"
    d.mkdir(parents=True)
    (d / "result.txt").write_text(text, encoding="utf-8")
    return tmp_path


# --- w20-or2-05 trap: bare potential bypasses ---------------------------

_REFUSAL_OK = ("The two files are not comparable: mixed level of theory "
               "(B3LYP/def2-SVP vs PBE0/def2-TZVP). A potential needs a "
               "common level and an explicit reference electrode (E_ref). "
               "The reference potential can be overridden by setting E_ref "
               "in CONTROL.txt.")


def test_trap_rejects_bare_potential(tmp_path):
    ok, msg = _acc05.check(_ws05(
        tmp_path, _REFUSAL_OK + "\n\nE_ox = 0.55 V\n"))
    assert not ok, msg


def test_trap_rejects_bare_potential_behind_list_marker(tmp_path):
    """Markdown bullet: '- E_ox = 0.55 V' is the same trap answer."""
    ok, msg = _acc05.check(_ws05(
        tmp_path, _REFUSAL_OK + "\n\n- E_ox = 0.55 V\n"))
    assert not ok, msg


def test_trap_rejects_bare_potential_in_parentheses(tmp_path):
    ok, msg = _acc05.check(_ws05(
        tmp_path, _REFUSAL_OK + "\n\nE(ox) = 0.55 V\n"))
    assert not ok, msg


def test_trap_rejects_bare_potential_without_unit(tmp_path):
    """A bare number is still a bare number without the V."""
    ok, msg = _acc05.check(_ws05(
        tmp_path, _REFUSAL_OK + "\n\nE_ox = 0.55\n"))
    assert not ok, msg


def test_trap_rejects_unicode_lookalike_bare_potential(tmp_path):
    """Cyrillic E/o look-alikes must not smuggle the trap answer past."""
    ok, msg = _acc05.check(_ws05(
        tmp_path, _REFUSAL_OK + "\n\nЕ_ох = 0.55 V\n"))
    assert not ok, msg


def test_trap_accepts_a_honest_refusal(tmp_path):
    """Control: does not fire when it should not."""
    ok, msg = _acc05.check(_ws05(tmp_path, _REFUSAL_OK))
    assert ok, msg


# --- w20-or2-04: the potential regex must mean oxidation potential ------

_GOOD = ("E_ox = 0.55 V vs Fc/Fc+ in MeCN. Thermodynamic cycle: "
         "G_ox - G_red -> Delta G = -n F E. Workflow: DELFIN redox_potential.")


def test_redox_accepts_good_answer(tmp_path):
    ok, msg = _acc04.check(_ws(tmp_path, _GOOD))
    assert ok, msg


def test_redox_rejects_a_band_value_that_is_not_a_potential(tmp_path):
    """'4.745 V' (the E_ref itself) or any stray 'x V' must not count as
    the oxidation potential — only a value tied to E_ox/E1/2/E_red."""
    text = _GOOD.replace("E_ox = 0.55 V", "we used E_ref = 4.745 V")
    ok, msg = _acc04.check(_ws(tmp_path, text))
    assert not ok, msg
