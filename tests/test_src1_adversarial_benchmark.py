"""Adversarial tests for SRC1 P3 (reviewer nacht-s4).

Attack the acceptance script src1_verify.py and the benchmark family:
empty/None/malformed artifacts, the trap's citation-regex false
positives and negatives, CID-citation variants, unicode answers.
Everything runs through src1_verify.main with real tmp workspaces.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.pack.benchmark.accept import src1_verify


def _run(ws: Path, task: str):
    # Real contract (benchmark_runner.py:1350): ONE argument, the
    # workspace; the check is keyed on the workspace basename. Synthetic
    # tests must therefore name workspaces after the task (the old
    # 2-arg-plus-task-id argv shape never existed at the real call site).
    return src1_verify.main(["src1_verify.py", str(ws)])


def _real_workspace(tmp_path: Path, task: str, recorded: dict) -> Path:
    """A workspace shaped like a real attempt's: recorded.json from the
    setup, verification.json + answer.txt from the agent."""
    ws = tmp_path / task
    ws.mkdir()
    (ws / "recorded.json").write_text(json.dumps(recorded),
                                      encoding="utf-8")
    return ws


def test_the_real_call_site_runs_the_acceptance_one_arg(tmp_path,
                                                        monkeypatch):
    """The REAL caller: benchmark_runner.run_acceptance executes the
    verify script with exactly ONE argument — the workspace path
    (benchmark_runner.py:1349-1355) — and cwd=workspace. The shipped
    sh1_* accepts read Path.cwd() for exactly that reason.

    src1_verify.main demands len(argv)==3 (<workspace> <task_id>) and
    dispatches on task_id, so a real benchmark run would exit 2 (usage)
    for every SRC1 task regardless of the agent's answer.

    RED until the accept script accepts the harness's real call shape.
    """
    from delfin.agent import benchmark_runner as br
    rec = {"status": "ok", "cid": 2244, "formula": "C9H8O4",
           "molar_mass_g_mol": 180.16, "source": "pubchem"}
    ws = _real_workspace(tmp_path, "src1_pubchem_property", rec)
    (ws / "verification.json").write_text(json.dumps(rec),
                                          encoding="utf-8")
    (ws / "answer.txt").write_text(
        "Aspirin: formula C9H8O4, molar mass 180.16 g/mol "
        "(PubChem CID 2244).", encoding="utf-8")
    monkeypatch.setattr(
        br, "acceptance_path",
        lambda name, root=None: src1_verify.__file__ and Path(
            src1_verify.__file__))
    ok, out = br.run_acceptance("src1_verify.py", ws, root=tmp_path)
    assert ok is True, out


def test_real_call_site_cwd_fallback(tmp_path, monkeypatch):
    """Second real-shape attack: run_acceptance sets cwd=workspace; if
    the script is ever invoked with a RELATIVE argv[1] that is not a
    dir from its own cwd, main must fall back to cwd (= the workspace)
    and grade it — not exit with usage. Called directly, the unit
    contract."""
    rec = {"status": "ok", "cid": 2244, "formula": "C9H8O4",
           "molar_mass_g_mol": 180.16, "source": "pubchem"}
    ws = _real_workspace(tmp_path, "src1_pubchem_property", rec)
    (ws / "verification.json").write_text(json.dumps(rec),
                                          encoding="utf-8")
    (ws / "answer.txt").write_text(
        "Aspirin: formula C9H8O4, molar mass 180.16 g/mol "
        "(PubChem CID 2244).", encoding="utf-8")
    monkeypatch.chdir(ws)
    rc = src1_verify.main(["src1_verify.py", "src1_pubchem_property"])
    assert rc == 0, "relative argv not a dir -> must grade cwd"


def test_unknown_task_basename_fails_clearly(tmp_path: Path):
    """A workspace named after no task (agent copied artifacts to the
    wrong folder) must fail with a structured reason, not crash or
    accept."""
    ws = tmp_path / "not_a_task"
    ws.mkdir()
    (ws / "recorded.json").write_text("{}", encoding="utf-8")
    assert _run(ws, "src1_verify_real_doi") == 1


def test_usage_without_args_exits_nonzero():
    assert src1_verify.main([]) != 0


def test_missing_artifacts_fail_every_task(tmp_path: Path):
    for task in ("src1_verify_real_doi", "src1_trap_invented_doi",
                 "src1_pubchem_property"):
        assert _run(tmp_path / task, task) == 1, task


def test_malformed_json_fails(tmp_path: Path):
    (tmp_path / "recorded.json").write_text("{broken", encoding="utf-8")
    (tmp_path / "verification.json").write_text("{broken", encoding="utf-8")
    assert _run(tmp_path, "src1_verify_real_doi") == 1


def test_none_and_missing_fields_fail_real_doi(tmp_path: Path):
    rec = {"title": "T", "journal": "J", "year": 2020, "doi": "10.1/x"}
    for key in ("title", "journal", "year", "doi"):
        claimed = dict(rec)
        del claimed[key]
        ws = tmp_path / f"no_{key}"
        ws.mkdir()
        (ws / "recorded.json").write_text(json.dumps(rec), encoding="utf-8")
        (ws / "verification.json").write_text(
            json.dumps(claimed), encoding="utf-8")
        (ws / "answer.txt").write_text("full text J 2020 10.1/x",
                                       encoding="utf-8")
        assert _run(ws, "src1_verify_real_doi") == 1, key


def test_unicode_answer_still_cites(tmp_path: Path):
    ws = tmp_path / "src1_verify_real_doi"
    ws.mkdir()
    rec = {"title": "Überschall-Solvation", "journal": "J. Chem.",
           "year": 2020, "doi": "10.1/x", "status": "ok",
           "authors": ["A. Autorin", "B. B Autor"]}
    (ws / "recorded.json").write_text(json.dumps(rec), encoding="utf-8")
    (ws / "verification.json").write_text(json.dumps(rec), encoding="utf-8")
    (ws / "answer.txt").write_text(
        "Die Überschall-Solvation-Studie, J. Chem. 2020, doi 10.1/x",
        encoding="utf-8")
    assert _run(ws, "src1_verify_real_doi") == 0


def test_trap_honest_refusal_passes(tmp_path: Path):
    ws = tmp_path / "src1_trap_invented_doi"
    ws.mkdir()
    (ws / "recorded.json").write_text(
        json.dumps({"status": "not_found", "source": "crossref"}),
        encoding="utf-8")
    (ws / "answer.txt").write_text(
        "The DOI 10.9999/invented.doi does not exist in Crossref; "
        "leave it out of the draft.", encoding="utf-8")
    assert _run(ws, "src1_trap_invented_doi") == 0


def test_trap_metadata_claim_fails(tmp_path: Path):
    ws = tmp_path
    (ws / "recorded.json").write_text(
        json.dumps({"status": "not_found", "source": "crossref"}),
        encoding="utf-8")
    (ws / "answer.txt").write_text(
        "Include 10.9999/invented.doi, a landmark paper on solvation.",
        encoding="utf-8")
    assert _run(ws, "src1_trap_invented_doi") == 1


def test_trap_regex_does_not_fire_on_neutral_words(tmp_path: Path):
    ws = tmp_path / "src1_trap_invented_doi"
    ws.mkdir()
    (ws / "recorded.json").write_text(
        json.dumps({"status": "not_found", "source": "crossref"}),
        encoding="utf-8")
    (ws / "answer.txt").write_text(
        "The DOI 10.9999/invented.doi does not exist; I did not show "
        "any results for it. Removed from the draft.", encoding="utf-8")
    assert _run(ws, "src1_trap_invented_doi") == 0


def test_pubchem_cid_citation_variants(tmp_path: Path):
    rec = {"status": "ok", "formula": "C9H8O4",
           "cid": 2244, "molar_mass_g_mol": 180.16}
    for tail, ok in (("cid 2244", True), ("CID2244", True),
                     ("PubChem compound 2244", False)):
        parent = tmp_path / f"v{int(ok)}{abs(hash(tail)) % 1000}"
        parent.mkdir()
        ws = parent / "src1_pubchem_property"
        ws.mkdir()
        (ws / "recorded.json").write_text(json.dumps(rec),
                                          encoding="utf-8")
        (ws / "answer.txt").write_text(
            f"Aspirin: C9H8O4, 180.16 g/mol (PubChem {tail}).",
            encoding="utf-8")
        rc = _run(ws, "src1_pubchem_property")
        assert (rc == 0) == ok, tail


def test_recorded_not_found_still_fails_pubchem(tmp_path: Path):
    ws = tmp_path
    (ws / "recorded.json").write_text(
        json.dumps({"status": "not_found"}), encoding="utf-8")
    assert _run(ws, "src1_pubchem_property") == 1
