"""QS ruling (nacht-s27, relayed by operator): the trap task's accepted
answer must state that the E_ref reference potential can be OVERRIDDEN
manually (delfin/config.py:1042-1044 — get_E_ref returns the user's
E_ref before the solvent table). An answer that names only the
solvent-table electrode value, without the override route, misses the
correction and must be rejected by the trap acceptance script.

Red on 6d5404f8: the accept script has no such requirement.
"""
from __future__ import annotations

import importlib.util
import pathlib

_PACK = (pathlib.Path(__file__).resolve().parents[1]
         / "delfin" / "agent" / "pack" / "benchmark" / "accept")

_HEADER = ("The two files are not comparable: mixed level (B3LYP/def2-SVP "
           "for ox, PBE0/def2-TZVP for red). "
           "A potential across different levels of theory is meaningless.\n")


def _load():
    spec = importlib.util.spec_from_file_location(
        "w20_or2_05_trap", _PACK / "w20_or2_05_redox_mixed_trap.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _check(workspace, text):
    ws = workspace / "or2" / "redox_mixed"
    ws.mkdir(parents=True, exist_ok=True)
    (ws / "result.txt").write_text(text, encoding="utf-8")
    return _load().check(workspace)


def test_rejects_answer_without_the_override_route(tmp_path):
    ok, msg = _check(tmp_path, _HEADER +
                     "Reference electrode: E_ref(acetonitrile) = 4.745 V "
                     "vs SHE.\n")
    assert not ok, msg


def test_accepts_answer_that_names_the_override_route(tmp_path):
    ok, msg = _check(tmp_path, _HEADER +
                     "Reference electrode: E_ref(acetonitrile) = 4.745 V "
                     "vs SHE. The E_ref value can be overridden manually "
                     "by setting E_ref in CONTROL.txt.\n")
    assert ok, msg


def test_accepts_answer_with_electrode_override_wording(tmp_path):
    ok, msg = _check(tmp_path, _HEADER +
                     "Reference electrode: E_ref = 4.745 V vs SHE "
                     "(override: set E_ref in CONTROL.txt).\n")
    assert ok, msg
