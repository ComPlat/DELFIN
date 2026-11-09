"""Tests for delfin/ml_train/datasets.py.

Synthetic records only -- no real structures, no CSD data, no real paths.
Runs on the py3.11 CI interpreter (no OPI required).
"""
import json

import pytest

from delfin.ml_train import datasets as ds


# ---------------------------------------------------------------- helpers


def _rec(**kw):
    base = {
        "identifier": "rec1",
        "source": "own-calculation",
        "licence": "CC0",
        "smiles": "CCO",  # ethanol: can be canonicalized, gives an InChIKey
    }
    base.update(kw)
    return ds.Record(**base)


def _sample_records():
    """Distinct simple molecules with a couple of targets."""
    return [
        ds.Record(
            identifier="eth",
            smiles="CCO",
            source="own-calculation",
            licence="CC0",
            targets=[ds.Target(name="homo", value=-6.4, unit="eV")],
            provenance={"origin": "own", "path": "calc/qm01"},
        ),
        ds.Record(
            identifier="meoh",
            smiles="CO",
            source="own-calculation",
            licence="CC0",
            targets=[ds.Target(name="homo", value=-6.8, unit="eV")],
        ),
        ds.Record(
            identifier="c6h6",
            smiles="c1ccccc1",
            source="cod",
            licence="CC0",
            targets=[ds.Target(name="homo", value=-6.1, unit="eV")],
        ),
    ]


# ---------------------------------------------------------------- guards


def test_validate_refuses_csd_source():
    rec = _rec(source="CSD refcode ABCDEF", licence="CC-BY-4.0")
    with pytest.raises(ds.LicenceRefused):
        ds.validate_record(rec)


def test_validate_refuses_ccdc_source():
    rec = _rec(source="CCDC deposition", licence="CC0")
    with pytest.raises(ds.LicenceRefused):
        ds.validate_record(rec)


def test_validate_refuses_unknown_licence():
    rec = _rec(source="own-calculation", licence="totally-made-up-licence")
    with pytest.raises(ds.LicenceRefused):
        ds.validate_record(rec)


def test_validate_accepts_cod_and_own():
    for source, licence in (("cod", "CC0"), ("own-calculation", "CC0")):
        record = _rec(source=source, licence=licence)
        ds.validate_record(record)  # must not raise


def test_hold_out_path_lockbox2():
    with pytest.raises(ds.HoldOutRefused):
        ds.guard_hold_out_path("agent_workspace/LOCKBOX2")


def test_hold_out_path_stage1_v2_txt():
    with pytest.raises(ds.HoldOutRefused):
        ds.guard_hold_out_path("justus_batch/competitors/stages_v2/stage1_v2.txt")


def test_hold_out_path_heldout_stage1_txt():
    with pytest.raises(ds.HoldOutRefused):
        ds.guard_hold_out_path("ANALYSE_2026_10_12_f3_haenger/heldout_stage1.txt")


def test_hold_out_path_pool_6000():
    with pytest.raises(ds.HoldOutRefused):
        ds.guard_hold_out_path("data/raw/pool_6000/part_remaining.csv")


def test_hold_out_path_ok_for_normal_file(tmp_path):
    # a benign path that merely shares a non-marker directory name must pass
    p = tmp_path / "my_data" / "CC0.csv"
    ds.guard_hold_out_path(str(p))  # must not raise


def test_record_provenance_naming_holdout_is_refused():
    rec = _rec(
        source="own-calculation",
        licence="CC0",
        provenance={"origin": "agent_workspace/LOCKBOX2"},
    )
    with pytest.raises(ds.HoldOutRefused):
        ds.validate_record(rec)


# ---------------------------------------------------------------- jsonl round trip


def test_jsonl_round_trip(tmp_path):
    records = _sample_records()
    path = tmp_path / "ds.jsonl"
    ds.save_jsonl(records, path)
    loaded = ds.load_jsonl(path)
    assert len(loaded) == len(records)
    for a, b in zip(records, loaded):
        assert a.identifier == b.identifier
        assert a.smiles == b.smiles
        assert a.source == b.source
        assert a.licence == b.licence
        assert a.targets == b.targets
        assert a.provenance == b.provenance


def test_load_jsonl_refuses_csd_record(tmp_path):
    path = tmp_path / "bad.jsonl"
    path.write_text(
        json.dumps(
            {
                "identifier": "bad",
                "source": "CSD refcode ZZZZZZ",
                "licence": "CC0",
                "smiles": "CCO",
            }
        )
        + "\n"
    )
    with pytest.raises(ds.LicenceRefused):
        ds.load_jsonl(path)


def test_load_jsonl_refuses_holdout_provenance(tmp_path):
    path = tmp_path / "hold.jsonl"
    path.write_text(
        json.dumps(
            {
                "identifier": "h",
                "source": "own-calculation",
                "licence": "CC0",
                "smiles": "CCO",
                "provenance": {"path": "agent_workspace/LOCKBOX2"},
            }
        )
        + "\n"
    )
    with pytest.raises(ds.HoldOutRefused):
        ds.load_jsonl(path)


# ---------------------------------------------------------------- extxyz round trip


def test_extxyz_round_trip(tmp_path):
    # xyz block for ethanol: N atoms then N lines
    xyz = "9\n\nC  -0.0000  0.0000  0.0000\n" + "\n".join(
        "X  0.0  0.0  0.0" for _ in range(8)
    )
    records = [
        ds.Record(
            identifier="eth",
            smiles="CCO",
            xyz=xyz,
            source="own-calculation",
            licence="CC0",
            targets=[ds.Target(name="homo", value=-6.4, unit="eV")],
        )
    ]
    path = tmp_path / "mol.xyz"
    ds.save_extxyz(records, path)
    loaded = ds.load_extxyz(path)
    assert len(loaded) == 1
    assert loaded[0].identifier == "eth"
    assert loaded[0].smiles == "CCO"
    assert loaded[0].targets[0].value == pytest.approx(-6.4)


def test_extxyz_round_trip_scientific_notation_target(tmp_path):
    # targets are written with f"{value:g}" (datasets.py:364), which emits
    # exponent form for |value| < 1e-4 or magnitude >= 1e6; the reader must
    # not lose the magnitude (bug: '5e-05Ha' read back as 5.0 with unit
    # 'e-05Ha'). Real absolute Hartree energies can be this small.
    small_abs_energy = -0.00005  # Hartree; f"{:g}" -> '-5e-05'
    xyz = "9\n\nC  -0.0000  0.0000  0.0000\n" + "\n".join(
        "X  0.0  0.0  0.0" for _ in range(8)
    )
    records = [
        ds.Record(
            identifier="eth",
            smiles="CCO",
            xyz=xyz,
            source="own-calculation",
            licence="CC0",
            targets=[ds.Target(name="e_abs", value=small_abs_energy, unit="Ha")],
        )
    ]
    path = tmp_path / "sci.xyz"
    ds.save_extxyz(records, path)
    loaded = ds.load_extxyz(path)
    assert len(loaded) == 1
    t = loaded[0].targets[0]
    assert t.value == pytest.approx(small_abs_energy), (
        f"{small_abs_energy} round-tripped to {t.value} (unit {t.unit!r})"
    )
    assert t.unit == "Ha"


def test_encode_decode_target_exponent_form():
    # low-level pair must survive exponent-form values across the :g boundary
    for v, unit in [(-5e-05, "Ha"), (-0.000123, "Ha"), (1e20, "au"),
                    (1234567.0, "kJ/mol"), (3.0e-7, "eV"), (-1234567.0, None)]:
        raw = ds._encode_target(ds.Target("x", v, unit))
        dec_v, dec_u = ds._decode_target(raw)
        assert dec_v == pytest.approx(v), f"{v!r} -> {raw!r} -> {dec_v!r}"
        assert dec_u == unit, f"{v!r} {unit!r} -> {dec_u!r}"


def test_extxyz_round_trip_metadata_with_spaces(tmp_path):
    # metadata goes into one comment line as whitespace-separated key=value
    # tokens (datasets.py:290-298); a multi-word identifier / source must not
    # be truncated at the first space on read back
    xyz = "9\n\nC  -0.0000  0.0000  0.0000\n" + "\n".join(
        "X  0.0  0.0  0.0" for _ in range(8)
    )
    records = [
        ds.Record(
            identifier="ethane conformer 1",
            smiles="CCO",
            xyz=xyz,
            source="own-calc batch A",
            licence="CC0",
            targets=[ds.Target(name="g", value=-5.5, unit="kcal/mol")],
        )
    ]
    path = tmp_path / "sp.xyz"
    ds.save_extxyz(records, path)
    loaded = ds.load_extxyz(path)
    assert len(loaded) == 1
    r = loaded[0]
    assert r.identifier == "ethane conformer 1", r.identifier
    assert r.source == "own-calc batch A", r.source
    assert r.targets[0].value == pytest.approx(-5.5)


# ---------------------------------------------------------------- dedup by InChIKey


def test_dedup_canonical_inchikey_drops_duplicate():
    # same molecule written two different ways -> same canonical InChIKey
    a = _rec(identifier="a", smiles="CCO")
    b = _rec(identifier="b", smiles="OCC")
    c = _rec(identifier="c", smiles="CO")
    kept, report = ds.deduplicate([a, b, c])
    assert len(kept) == 2  # ethanol + methanol
    assert {r.identifier for r in kept} == {"a", "c"} or "b" in {r.identifier for r in kept}
    dropped = report.dropped
    assert len(dropped) == 1
    # the dropped one reports which kept record it collided with
    assert dropped[0].duplicate_of in kept_ids(kept)


def kept_ids(kept):
    return {r.identifier for r in kept}


def test_dedup_report_counts(tmp_path):
    a = _rec(identifier="a", smiles="CCO")
    b = _rec(identifier="b", smiles="OCC")
    kept, report = ds.deduplicate([a, b])
    assert report.dropped_count == 1
    assert report.kept_count == 1
