"""An Open Babel UFF minimisation gives the same result whatever ran before it.

Open Babel hands every caller the one process-wide UFF instance, and its Setup keeps
the previous setup when the new molecule has the same atoms and bonds.  Two frames of
one complex that differ only in their coordinates then minimise differently depending
on which frame the process minimised before.  The batch UFF of the topology enumerator
runs in a pool of min(n, cpu_count, 64) worker processes, so on a machine with fewer
cores than inputs (a 4-core CI runner) the predecessor of an input in its worker is up to
the scheduler, and two builds of [Pt](Cl)(Cl)(N)N differed (CI, PR #548).

The two inputs are batch-UFF inputs of that build (trans-[PtCl2(NH3)2] with the N-H
hydrogens turned); the second minimises differently after the first than after an
unrelated molecule.
"""

from __future__ import annotations

import pytest

pytest.importorskip("openbabel")

from delfin import cli_manta  # noqa: E402

_CONSTRAINTS = {
    "fix_atoms": [1, 0, 2, 3, 4],
    "distances": [(1, 0, 2.03), (1, 2, 2.03), (1, 3, 2.3), (1, 4, 2.3)],
    "angles": [],
    "torsions": [],
    "_soft_donor_meta": {"metal_indices": [1], "donor_indices": [0, 2, 3, 4],
                         "pairs": [(1, 0), (1, 2), (1, 3), (1, 4)], "class_label": "sigma"},
}
_HEAVY = (
    "N        1.172021     1.172021     1.172021\n"
    "Pt       0.000000     0.000000     0.000000\n"
    "N       -1.172021    -1.172021     1.172021\n"
    "Cl      -1.327906     1.327906    -1.327906\n"
    "Cl       1.327906    -1.327906    -1.327906\n"
)
_BEFORE = _HEAVY + (
    "H        1.502057     2.128829     0.939266\n"
    "H        1.559935     0.933163     2.122726\n"
    "H       -0.799741    -1.762245     1.930847\n"
    "H       -2.197768    -1.235264     1.066661\n"
)
_TARGET = _HEAVY + (
    "H        1.316678     1.051966     2.190186\n"
    "H        1.721915     1.986628     0.848408\n"
    "H       -1.603408    -0.853511     2.063798\n"
    "H       -1.402801    -2.152698     0.942411\n"
)
_UNRELATED = "O 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"


def _set_the_champion_construction_env(monkeypatch):
    for key, value in cli_manta.construction_env("champion", environ={}).items():
        monkeypatch.setenv(key, value)
    monkeypatch.setenv("DELFIN_DETERMINISTIC", "1")


def test_the_champion_minimises_every_input_with_its_own_setup(monkeypatch):
    _set_the_champion_construction_env(monkeypatch)
    from delfin.manta.openbabel_optimize import _optimize_xyz_openbabel

    # Open Babel's setup is that of the first of a run of same-atoms-and-bonds inputs,
    # so each order starts after the unrelated molecule.
    _optimize_xyz_openbabel(_UNRELATED, 500, None)
    after_unrelated = _optimize_xyz_openbabel(_TARGET, 500, _CONSTRAINTS)
    _optimize_xyz_openbabel(_UNRELATED, 500, None)
    _optimize_xyz_openbabel(_BEFORE, 500, _CONSTRAINTS)
    after_sibling = _optimize_xyz_openbabel(_TARGET, 500, _CONSTRAINTS)
    assert after_sibling == after_unrelated


def test_the_switch_is_part_of_the_champion():
    assert cli_manta.construction_env("champion", environ={})["DELFIN_OB_FULL_SETUP"] == "1"


def _manifold(out):
    import json

    manifest = json.loads((out / "manifest.json").read_text())
    return [(item["label"], cli_manta._atom_lines((out / item["file"]).read_text()))
            for item in manifest["isomers"]]


@pytest.mark.slow
def test_one_uff_worker_builds_what_many_workers_build(monkeypatch, tmp_path):
    """With one batch-UFF worker every input follows the previous one in the same
    process; with many workers each input mostly gets a process of its own.  Before the
    fix the two manifolds differed (and with 4 workers nine manifolds came out of 24
    builds)."""
    for key in [k for k in __import__("os").environ if k.startswith("DELFIN_")]:
        monkeypatch.delenv(key, raising=False)
    built = []
    for workers in ("1", "64"):
        monkeypatch.setenv("DELFIN_MAX_PROCESS_WORKERS", workers)
        out = tmp_path / f"w{workers}"
        assert cli_manta.main(["[Pt](Cl)(Cl)(N)N", "-o", str(out), "-q"]) == 0
        built.append(_manifold(out))
    assert built[0] and all(lines for _label, lines in built[0])
    assert built[0] == built[1]
