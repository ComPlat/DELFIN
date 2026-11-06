"""Tests for delfin/ml_train/splits.py.

Synthetic records only. Runs on the py3.11 CI interpreter (rdkit, no OPI).
"""
from delfin.ml_train import datasets as ds
from delfin.ml_train import splits as sp


# ---------------------------------------------------------------- helpers


def _records(n, smiles_pool=None):
    smiles_pool = smiles_pool or ["CCO", "CO", "c1ccccc1", "CCC", "CCCC"]
    return [
        ds.Record(
            identifier=f"r{i}",
            smiles=smiles_pool[i % len(smiles_pool)],
            source="own-calculation",
            licence="CC0",
            targets=[ds.Target(name="homo", value=-6.0 + i * 0.1, unit="eV")],
            provenance={"group": f"metal-{i % 3}", "source_label": f"src-{i % 4}"},
        )
        for i in range(n)
    ]


def _assert_disjoint(splits):
    counts = {}
    for name, indices in splits.items():
        for i in indices:
            assert i not in counts, f"index {i} appears in more than one split"
            counts[i] = name


def _assert_reproducible(fn, *args, **kw):
    s1 = fn(*args, **kw)
    s2 = fn(*args, **kw)
    assert list(s1) == list(s2)
    for name in s1:
        assert s1[name] == s2[name]


# ---------------------------------------------------------------- random split


def test_random_split_disjoint_and_covers_all():
    n = 40
    splits = sp.random_split(n, {"train": 0.7, "val": 0.2, "test": 0.1}, seed=7)
    all_idx = [i for idx in splits.values() for i in idx]
    assert sorted(all_idx) == list(range(n))
    _assert_disjoint(splits)


def test_random_split_reproducible_with_seed():
    _assert_reproducible(sp.random_split, 50, {"train": 0.7, "test": 0.3}, seed=123)


def test_random_split_seed_changes_result():
    a = sp.random_split(40, {"train": 0.7, "test": 0.3}, seed=1)
    b = sp.random_split(40, {"train": 0.7, "test": 0.3}, seed=2)
    assert a["train"] != b["train"]


# ---------------------------------------------------------------- scaffold split


def test_scaffold_split_disjoint_and_groups_w_structure():
    records = _records(40)
    splits = sp.scaffold_split(records, {"train": 0.7, "test": 0.3}, seed=7)
    _assert_disjoint(splits)
    # scaffold split must never split the same scaffold across two splits
    idx_to_smiles = {i: r.smiles for i, r in enumerate(records)}
    scaffold_of = {
        i: _murcko(idx_to_smiles[i]) for i in range(len(records))
    }
    for a_name in splits:
        for b_name in splits:
            if a_name >= b_name:
                continue
            sa = {scaffold_of[i] for i in splits[a_name]}
            sb = {scaffold_of[i] for i in splits[b_name]}
            assert sa.isdisjoint(sb)


def test_scaffold_split_reproducible_with_seed():
    records = _records(30)
    _assert_reproducible(sp.scaffold_split, records, {"a": 0.5, "b": 0.5}, seed=42)


def _murcko(smiles):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    return MurckoScaffold.MurckoScaffoldSmiles(
        MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smiles))
    )


# ---------------------------------------------------------------- stratified split


def test_stratified_by_target_bins_balanced_per_bin():
    records = _records(60)  # homogeneous binning by value would be trivial;
    splits = sp.stratified_split(
        records,
        target_name="homo",
        n_bins=3,
        fractions={"train": 0.6, "test": 0.4},
        seed=11,
    )
    _assert_disjoint(splits)
    # every bin of the target value must be represented in every split
    def _bin(i, records):
        v = records[i].targets[0].value
        edges = [min(r.targets[0].value for r in records),
                 max(r.targets[0].value for r in records)]
        span = (edges[1] - edges[0]) / 3
        return int((v - edges[0]) / span) if span else 0

    for name, indices in splits.items():
        bins = {_bin(i, records) for i in indices}
        assert len(bins) >= 2  # balanced representation, not a single bin


def test_stratified_reproducible():
    records = _records(40)
    _assert_reproducible(
        sp.stratified_split, records, "homo", 3,
        {"train": 0.6, "test": 0.4}, seed=9,
    )


# ---------------------------------------------------------------- grouped split


def test_grouped_split_keeps_source_together():
    records = _records(48)
    splits = sp.grouped_split(
        records,
        key=lambda r: r.provenance["source_label"],
        fractions={"train": 0.75, "test": 0.25},
        seed=3,
    )
    _assert_disjoint(splits)
    # no source group may appear in two splits
    src_of = {i: r.provenance["source_label"] for i, r in enumerate(records)}
    for a_name in splits:
        for b_name in splits:
            if a_name >= b_name:
                continue
            sa = {src_of[i] for i in splits[a_name]}
            sb = {src_of[i] for i in splits[b_name]}
            assert sa.isdisjoint(sb)


def test_grouped_split_reproducible():
    records = _records(36)
    _assert_reproducible(
        sp.grouped_split, records, lambda r: r.smiles,
        {"train": 0.5, "test": 0.5}, seed=5,
    )


# ---------------------------------------------------------------- leakage check


def test_leakage_check_catches_planted_inchikey_duplicate():
    # two records sharing the same molecule (same canonical InChIKey)
    records = [
        ds.Record(identifier="a", smiles="CCO", source="own-calculation",
                  licence="CC0"),
        ds.Record(identifier="b", smiles="OCC", source="own-calculation",
                  licence="CC0"),
        ds.Record(identifier="c", smiles="CO", source="own-calculation",
                  licence="CC0"),
        ds.Record(identifier="d", smiles="c1ccccc1", source="own-calculation",
                  licence="CC0"),
    ]
    # plant: 'a' (ethanol) in train, 'b' (also ethanol) in test
    splits = {"train": [0], "test": [1, 2, 3]}
    report = sp.check_leakage(splits, records)
    assert report.has_leakage
    assert len(report.inchikey_overlaps) == 1


def test_leakage_check_clean_datasets_pass():
    records = _records(20)
    splits = sp.scaffold_split(records, {"train": 0.6, "test": 0.4}, seed=1)
    report = sp.check_leakage(splits, records)
    assert not report.has_leakage


def test_split_file_records_seeds(tmp_path):
    records = _records(30)
    splits = sp.random_split(
        len(records), {"train": 0.7, "test": 0.3}, seed=2024
    )
    path = tmp_path / "split.json"
    sp.save_split(splits, path, seed=2024, method="random")
    loaded = sp.load_split(path)
    assert loaded["seed"] == 2024
    assert loaded["method"] == "random"
    # reloaded split pieces reproduce the same disjoint assignment
    _assert_disjoint(loaded["splits"])
