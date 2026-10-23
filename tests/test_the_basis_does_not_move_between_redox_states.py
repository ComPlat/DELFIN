"""A ligand that drifts must not change which basis set it gets.

The first coordination sphere decides which atoms carry the larger metal
basis, and membership is a distance test. Metal-ligand bonds lengthen on
reduction, so a ligand can cross the cutoff between two states — and a redox
potential is the difference between two states, so the basis-set step lands
straight in the reported number.

Measured on a real run (Co complex, two DMF ligands, CONTROL scale 1.3): one
Co-O bond grew 2.11 -> 2.34 -> 2.58 A over two reductions and crossed the
2.50 A cutoff at the second one. That state was computed with def2-SVP on an
oxygen every other state had def2-TZVP on, its energy came out 0.086 Ha too
high, and the second reduction potential was reported as -4.16 V instead of
about -1.8 V.
"""

from pathlib import Path

import pytest

from delfin.xyz_io import (
    _first_sphere_indices,
    _parse_xyz_atoms,
    read_xyz_and_create_input3,
)

CONFIG = {
    "first_coordination_sphere_metal_basisset": "yes",
    "first_coordination_sphere_scale": "1.3",
    "functional": "PBE0",
    "disp_corr": "D4",
    "ri_jkx": "RIJCOSX",
    "aux_jk": "def2/J",
    "implicit_solvation_model": "CPCM",
    "geom_opt": "OPT",
    "freq_type": "FREQ",
    "initial_guess": "PModel",
    "relativity": "none",
    "maxcore": 4500,
    "PAL": 8,
    "temperature": "298.15",
    "maxiter": 125,
}

# Co-O cutoff at scale 1.3 is 1.3 * (1.26 + 0.66) = 2.496 A.
BOUND = 2.10      # both ligands clearly inside
DRIFTED = 2.58    # one ligand past the cutoff, as after two reductions


def _geometry(second_o_distance: float) -> str:
    rows = [
        ("Co", 0.0, 0.0, 0.0),
        ("O", 2.10, 0.0, 0.0),
        ("O", -second_o_distance, 0.0, 0.0),
        ("H", 2.10, 0.96, 0.0),
        ("H", -second_o_distance - 0.96, 0.0, 0.0),
    ]
    body = "\n".join(f"{el} {x:.6f} {y:.6f} {z:.6f}" for el, x, y, z in rows)
    return f"{len(rows)}\ncomment\n{body}\n"


def _write(tmp_path: Path, name: str, distance: float) -> Path:
    path = tmp_path / name
    path.write_text(_geometry(distance), encoding="utf-8")
    return path


def _newgto_atoms(inp: Path) -> list[str]:
    found, n, started = [], 0, False
    for line in inp.read_text(encoding="utf-8").splitlines():
        if line.startswith("* xyz"):
            started = True
            continue
        if started and line[:1].isupper():
            n += 1
            if "NewGTO" in line:
                found.append(f"{line.split()[0]}{n}")
    return found


def _build(src: Path, out: Path) -> list[str]:
    read_xyz_and_create_input3(
        str(src), str(out), 0, 4, "dmf", ["Co"], "def2-TZVP", "def2-SVP", dict(CONFIG), "",
    )
    return _newgto_atoms(out)


def test_the_bare_rule_does_drop_a_drifted_ligand():
    """The distance test itself — this is what used to reach the input."""
    def sphere(distance):
        lines = _geometry(distance).splitlines()[2:]
        atoms = _parse_xyz_atoms([f"{l}\n" for l in lines])
        metals = [i for i, a in enumerate(atoms) if a["elem"] == "Co"]
        return _first_sphere_indices(atoms, metals, 1.3, None)

    assert len(sphere(BOUND)) == 2      # both oxygens
    assert len(sphere(DRIFTED)) == 1    # the drifted one is gone


def test_a_drifted_ligand_keeps_the_basis_the_reference_state_gave_it(tmp_path):
    reference = _build(_write(tmp_path, "ref.xyz", BOUND), tmp_path / "ref.inp")
    drifted = _build(_write(tmp_path, "drift.xyz", DRIFTED), tmp_path / "drift.inp")

    assert reference == ["Co1", "O2", "O3"]
    assert drifted == reference, (
        "the drifted state was given a different basis than the state its "
        "energy is subtracted from"
    )


def test_a_different_molecule_is_decided_afresh(tmp_path):
    """Freezing is per composition, so explicit solvent or another system
    does not inherit a sphere that was never about it."""
    _build(_write(tmp_path, "ref.xyz", BOUND), tmp_path / "ref.inp")

    other = tmp_path / "other.xyz"
    other.write_text("2\ncomment\nCo 0.0 0.0 0.0\nO 2.10 0.0 0.0\n", encoding="utf-8")
    assert _build(other, tmp_path / "other.inp") == ["Co1", "O2"]


def test_the_first_geometry_to_ask_decides(tmp_path):
    """A run that starts at the drifted state has no reference to inherit, so
    it freezes there — consistently, for every state it then computes."""
    drifted = _build(_write(tmp_path, "drift.xyz", DRIFTED), tmp_path / "drift.inp")
    later = _build(_write(tmp_path, "ref.xyz", BOUND), tmp_path / "ref.inp")
    assert drifted == ["Co1", "O2"]
    assert later == drifted


def test_two_jobs_beside_each_other_do_not_share_a_sphere(tmp_path):
    """The sphere belongs to one run, not to whatever directory it sits in."""
    jobs = {}
    for name, distance in (("jobA", BOUND), ("jobB", DRIFTED)):
        job = tmp_path / name
        job.mkdir()
        (job / "CONTROL.txt").write_text("charge=0\n", encoding="utf-8")
        jobs[name] = _build(_write(job, "start.xyz", distance), job / "start.inp")
    assert jobs["jobA"] == ["Co1", "O2", "O3"]
    assert jobs["jobB"] == ["Co1", "O2"]


def test_an_occupier_subfolder_and_the_run_root_agree(tmp_path):
    """This is the real path: OCCUPIER optimises in a subfolder, the ligand
    drifts, and the next input is built in the run root from that geometry."""
    job = tmp_path / "job"
    (job / "red_step_2_OCCUPIER").mkdir(parents=True)
    (job / "CONTROL.txt").write_text("charge=0\n", encoding="utf-8")

    # the reference state is set up first, in the run root
    root = _build(_write(job, "initial.xyz", BOUND), job / "initial.inp")
    # OCCUPIER hands back a geometry whose ligand has drifted out
    sub = job / "red_step_2_OCCUPIER"
    candidate = _build(_write(sub, "input3.xyz", DRIFTED), sub / "input3.inp")
    # ...and the final step is built from it, back in the run root
    final = _build(_write(job, "red_step_2.xyz", DRIFTED), job / "red_step_2.inp")

    assert root == ["Co1", "O2", "O3"]
    assert candidate == root
    assert final == root
    assert (job / ".delfin_first_sphere.json").is_file()
