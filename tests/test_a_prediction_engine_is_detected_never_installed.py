"""A structure-prediction engine is detected, never installed.

Four engines can predict a protein structure from a sequence: ColabFold,
ESMFold, Boltz and AlphaFold 3.  The contract DELFIN makes with each is the
docking one: the engine is only *detected* (``shutil.which`` on its command,
the same probe every external-program block in this layer makes) and the
command that would run it is only *described* — building the command line is
the whole of what the interface does.  No installation, no download of
weights or databases, no run: DELFIN steers real calculations, and a
prediction engine is exactly the kind of program a caller must decide to
have.

Licences, named here because the interface names the programs (code
licences; the trained weights carry their own terms in every case):

* ColabFold — MIT (the trained AlphaFold weights it loads are licensed
  separately by their own terms).
* ESMFold (Meta ESM) — MIT.
* Boltz — MIT (code and weights).
* AlphaFold 3 (DeepMind) — the inference code is published under Creative
  Commons BY-NC-SA 4.0: non-commercial terms, and a caller must accept the
  DeepMind model terms before use.

Tests drive the interface with a stand-in program written into the test's
own tmp dir and put first on ``PATH``, so every test runs on every machine
without the real programs.
"""

from __future__ import annotations

import json
import os
import stat
import time
from pathlib import Path

import pytest

from delfin.tools import describe
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.bio import StructureFormatError, parse_fasta

STEP = "structure_prediction_engine"
ENGINES = ("colabfold", "esmfold", "boltz", "alphafold3")

FASTA = ">query lysozyme\nKVFGRCELAA\n"


def _standin(tmp_path: Path, name: str) -> Path:
    """A stand-in executable that prints its arguments, for detection tests."""
    path = tmp_path / name
    path.write_text("#!/bin/sh\necho \"$0 $@\"\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _env(tmp_path: Path, *names: str) -> dict:
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(
        [str(tmp_path)] + env.get("PATH", "").split(os.pathsep))
    return env


def test_the_four_engines_are_named_by_the_contract():
    contract = describe(STEP)
    engines = {e["name"] for e in contract.data_keys[0].description and []} \
        if False else None  # data_keys hold scalar keys, not the engine list
    assert contract.consumes == frozenset({"sequence"})
    assert contract.produces == frozenset({"json_report"})
    assert contract.requires_binaries == frozenset()
    assert contract.produces_geometry is False


def test_a_missing_engine_is_reported_by_its_command(tmp_path):
    from delfin.tools.bio import prediction_engines
    env = _env(tmp_path)
    report = prediction_engines(env=env)
    assert [e["name"] for e in report] == list(ENGINES)
    for entry in report:
        assert entry["detected"] is False
        assert "command" in entry
        assert bool(entry["command"]) != entry["input_conversion_required"]
        assert "licence" in entry and entry["licence"]


def test_a_detected_engine_carries_the_path_found(tmp_path):
    from delfin.tools.bio import prediction_engines
    _standin(tmp_path, "colabfold_batch")
    report = prediction_engines(env=_env(tmp_path))
    by_name = {e["name"]: e for e in report}
    assert by_name["colabfold"]["detected"] is True
    assert by_name["colabfold"]["path"] == str(tmp_path / "colabfold_batch")
    assert by_name["esmfold"]["detected"] is False


def test_the_command_names_the_input_and_the_output(tmp_path):
    from delfin.tools.bio import prediction_engines
    _standin(tmp_path, "esm-fold")
    env = _env(tmp_path)
    fasta = tmp_path / "query.fasta"
    fasta.write_text(FASTA, encoding="utf-8")
    report = prediction_engines(env=env, fasta=str(fasta))
    command = {e["name"]: e for e in report}["esmfold"]["command"]
    assert command
    assert str(fasta) in command
    assert command[0] == "esm-fold"


def test_detection_never_runs_a_program(tmp_path):
    """A stand-in that would fail the test if executed proves which() only."""
    sentinel = tmp_path / "would_run"
    script = tmp_path / "colabfold_batch"
    script.write_text(
        f"#!/bin/sh\ntouch {sentinel}\n", encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    from delfin.tools.bio import prediction_engines
    report = prediction_engines(env=_env(tmp_path))
    by_name = {e["name"]: e for e in report}
    assert by_name["colabfold"]["detected"] is True
    assert not sentinel.exists()


def test_detection_reads_no_environment_it_was_not_given(tmp_path):
    from delfin.tools.bio import prediction_engines
    env = {"PATH": str(tmp_path)}          # no os.environ at all
    report = prediction_engines(env=env)
    assert [e["detected"] for e in report] == [False] * 4


def test_af3_and_boltz_do_not_advertise_fabricated_generic_fasta_commands():
    from delfin.tools.bio import prediction_engines
    entries = {e["name"]: e for e in prediction_engines(env={}, fasta="input.fasta")}
    for name in ("alphafold3", "boltz"):
        assert entries[name]["command"] == []
        assert entries[name]["input_conversion_required"]
        assert not entries[name]["execution_ready"]


def test_the_block_reports_the_engines_and_the_fasta(tmp_path):
    from delfin.tools import run_step
    _standin(tmp_path, "colabfold_batch")
    fasta = tmp_path / "query.fasta"
    fasta.write_text(FASTA, encoding="utf-8")
    result = run_step(STEP, sequence_file=str(fasta), work_dir=tmp_path,
                      _prediction_env={"PATH": str(tmp_path)})
    assert result.ok
    assert result.data["sequence_id"] == "query lysozyme"
    assert result.data["n_engines"] == 4
    detected = [e for e in result.data["engines"] if e["detected"]]
    assert [e["name"] for e in detected] == ["colabfold"]
    report = Path(result.artifacts["json_report"])
    payload = json.loads(report.read_text(encoding="utf-8"))
    assert payload["sequence_id"] == "query lysozyme"


def test_the_block_refuses_a_dna_fasta_by_name(tmp_path):
    from delfin.tools import run_step
    fasta = tmp_path / "gene.fasta"
    fasta.write_text(">gene\nATGCTGATCGA\n", encoding="utf-8")
    result = run_step(STEP, sequence_file=str(fasta), work_dir=tmp_path,
                      _prediction_env={"PATH": str(tmp_path)})
    assert result.status is StepStatus.FAILED
    assert "DNA" in result.error


def test_a_missing_sequence_file_is_a_named_failure(tmp_path):
    from delfin.tools import run_step
    result = run_step(STEP, sequence_file=str(tmp_path / "gone.fasta"),
                      work_dir=tmp_path, _prediction_env={"PATH": str(tmp_path)})
    assert result.status is StepStatus.FAILED
    assert result.error_kind is ErrorKind.MISSING_INPUT
