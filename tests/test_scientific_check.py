"""The scientific-correctness critic is injected into get_calc_info so the agent
can't inspect a result without seeing its red flags (DELFIN's domain moat:
don't trust a numerically-done but scientifically-wrong result)."""

from __future__ import annotations

from pathlib import Path

from delfin.agent.api_client import _doc_executor


def test_flags_a_bad_result(tmp_path):
    (tmp_path / "job.out").write_text(
        "Some ORCA output\nSCF NOT CONVERGED\n")          # no clean termination
    info = {"calc_id": "x", "path": str(tmp_path)}
    out = _doc_executor._inject_scientific_check(info)
    assert "scientific_check" in out
    assert out["scientific_check"]["worst"] == "error"
    assert "SCF" in out["scientific_check"]["report"]
    assert "do NOT" in out["scientific_check"]["note"]


def test_no_path_is_unchanged():
    info = {"calc_id": "x"}                                # no path
    assert _doc_executor._inject_scientific_check(info) == info


def test_folder_without_outputs_is_unchanged(tmp_path):
    info = {"calc_id": "x", "path": str(tmp_path)}         # empty dir, no .out
    assert "scientific_check" not in _doc_executor._inject_scientific_check(info)


def test_wired_into_get_calc_info():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert "info = self._inject_scientific_check(info)" in src
