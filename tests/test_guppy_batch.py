"""Batch GUPPY orchestration tests (M3)."""
from pathlib import Path

import pytest

from delfin import guppy_batch, guppy_sampling


def test_read_smiles_lines_plain_text(tmp_path: Path) -> None:
    input_file = tmp_path / "smiles.txt"
    input_file.write_text(
        "# header comment\n"
        "[Fe+2].N.N.N.N\n"
        "\n"
        "CCO ethanol\n"
        "* star comment\n"
        "O water\n",
        encoding="utf-8",
    )
    entries = guppy_sampling._read_smiles_lines(input_file)
    assert entries == [
        ("[Fe+2].N.N.N.N", None),
        ("CCO", "ethanol"),
        ("O", "water"),
    ]


def test_read_smiles_lines_csv_with_header(tmp_path: Path) -> None:
    input_file = tmp_path / "smiles.csv"
    input_file.write_text(
        "smiles,name,charge\n"
        "[Cu+2].N.N.N.N,tetraamminecopper,2\n"
        "O,water,0\n"
        "[Pt+2].[Cl-].[Cl-].N.N,cisplatin,0\n",
        encoding="utf-8",
    )
    entries = guppy_sampling._read_smiles_lines(input_file)
    # Charge column must be ignored — always derived from SMILES.
    assert [e[0] for e in entries] == [
        "[Cu+2].N.N.N.N",
        "O",
        "[Pt+2].[Cl-].[Cl-].N.N",
    ]
    assert [e[1] for e in entries] == ["tetraamminecopper", "water", "cisplatin"]


def test_read_smiles_lines_csv_without_header(tmp_path: Path) -> None:
    input_file = tmp_path / "smiles.csv"
    input_file.write_text("[Fe+2],iron\nCCO,ethanol\n", encoding="utf-8")
    entries = guppy_sampling._read_smiles_lines(input_file, name_column=1)
    assert entries == [("[Fe+2]", "iron"), ("CCO", "ethanol")]


def test_read_smiles_lines_errors_on_missing_file(tmp_path: Path) -> None:
    missing = tmp_path / "nope.txt"
    with pytest.raises(FileNotFoundError):
        guppy_sampling._read_smiles_lines(missing)


def test_read_smiles_lines_errors_on_empty(tmp_path: Path) -> None:
    empty = tmp_path / "empty.txt"
    empty.write_text("# only comments\n\n*\n", encoding="utf-8")
    with pytest.raises(ValueError):
        guppy_sampling._read_smiles_lines(empty)


def test_slug_derivation_uses_name_when_present() -> None:
    slug = guppy_batch._derive_slug("CCO", "Ethanol!!", 3)
    assert slug.startswith("0003_")
    assert "Ethanol" in slug


def test_slug_derivation_falls_back_to_hash() -> None:
    slug_a = guppy_batch._derive_slug("CCO", None, 1)
    slug_b = guppy_batch._derive_slug("CCO", None, 1)
    slug_c = guppy_batch._derive_slug("CC", None, 1)
    assert slug_a == slug_b  # deterministic
    assert slug_a != slug_c  # SMILES-dependent
    assert slug_a.startswith("0001_")


def test_run_sampling_batch_per_entry_isolation(tmp_path: Path, monkeypatch) -> None:
    """Each batch entry gets its own workdir and derives its own charge."""
    captured = []

    def fake_run_sampling(**kwargs):
        # Record charge-derivation by reading input.txt the batch wrote.
        input_text = Path(kwargs["input_file"]).read_text(encoding="utf-8").strip()
        derived = guppy_sampling._derive_charge_from_smiles(input_text)
        captured.append({
            "workdir": str(kwargs["workdir"]),
            "charge_override": kwargs["charge"],
            "derived": derived,
            "start_strategy": kwargs["start_strategy"],
        })
        # Simulate a successful run: write the best-coordination sentinel.
        out = Path(kwargs["output_file"])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("1\ncomment -1.0\nC 0.0 0.0 0.0\n", encoding="utf-8")
        out.with_name("best_coordniation.xyz").write_text(
            "1\nrun_01 -1.23\nC 0.0 0.0 0.0\n", encoding="utf-8"
        )
        out.with_name("best_coordniation_pre_goat.xyz").write_text(
            "1\nrun_01 -1.23\nC 0.0 0.0 0.0\n", encoding="utf-8"
        )
        return 0

    monkeypatch.setattr(guppy_batch, "run_sampling", fake_run_sampling)

    entries = [
        ("[Fe+2].N.N.N.N", "iron_tetraammine"),
        ("[Cu+2].N.N.N.N", "copper_tetraammine"),
        ("O", None),
    ]
    shared = guppy_batch.BatchSharedArgs(
        runs=1, pal=1, maxcore=100, parallel_jobs=1, goat_topk=0, allow_partial=True,
    )
    rc = guppy_batch.run_sampling_batch(
        entries, base_workdir=tmp_path / "GUPPY_BATCH", shared_args=shared
    )
    assert rc == 0

    # Each entry called run_sampling once with its own workdir.
    assert len(captured) == 3
    workdirs = {c["workdir"] for c in captured}
    assert len(workdirs) == 3
    # Charge always auto-derived; override never sent.
    assert all(c["charge_override"] is None for c in captured)
    # Derived charges match expectation.
    charges = [c["derived"] for c in captured]
    assert charges == [2, 2, 0]
    # Default start strategy propagates.
    assert all(c["start_strategy"] == "isomers" for c in captured)

    # Aggregate file exists and has all three entries.
    aggregate = tmp_path / "GUPPY_BATCH" / "batch_results.json"
    assert aggregate.exists()
    import json as _json
    payload = _json.loads(aggregate.read_text(encoding="utf-8"))
    assert payload["total_entries"] == 3
    assert payload["ok"] == 3
    assert {e["resolved_charge"] for e in payload["entries"]} == {0, 2}


def test_run_sampling_batch_row_mode(tmp_path: Path, monkeypatch) -> None:
    """--row picks one entry; charge still derived fresh for that entry."""
    calls = []

    def fake_run_sampling(**kwargs):
        calls.append(kwargs["workdir"])
        out = Path(kwargs["output_file"])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.with_name("best_coordniation.xyz").write_text(
            "1\nrun_01 -1.0\nC 0 0 0\n", encoding="utf-8"
        )
        out.with_name("best_coordniation_pre_goat.xyz").write_text(
            "1\nrun_01 -1.0\nC 0 0 0\n", encoding="utf-8"
        )
        return 0

    monkeypatch.setattr(guppy_batch, "run_sampling", fake_run_sampling)

    entries = [("CCO", "a"), ("[Fe+2]", "b"), ("O", "c")]
    shared = guppy_batch.BatchSharedArgs(runs=1, pal=1, maxcore=100, parallel_jobs=1)
    rc = guppy_batch.run_sampling_batch(
        entries, base_workdir=tmp_path / "GB", shared_args=shared, row=2
    )
    assert rc == 0
    assert len(calls) == 1  # only row 2 executed
