"""XE1 benchmark target #5: the four structure tasks load and grade end to end.

For each task the setup builds a fixture workspace (run_setup); the model's
correct action is simulated by writing the NEW xyz file the preset would
produce; the acceptance script (run_acceptance) grades the WRITTEN
ARTIFACT — the grader must exit 0 only when the written file is correct.
The no-overwrite task is also exercised in the NEGATIVE: modifying the
input must fail the acceptance.
"""
from __future__ import annotations

from pathlib import Path

from delfin.agent.benchmark_runner import run_acceptance, run_setup

_REPO = Path(__file__).resolve().parents[1]  # repo root


def _write(ws: Path, rel: str, text: str) -> None:
    p = ws / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def _read(ws: Path, rel: str) -> str:
    return (ws / rel).read_text(encoding="utf-8")


def test_renumber_setup_and_verify_grade_the_written_file(tmp_path):
    ok, out = run_setup("xe1_renumber_seed.py", tmp_path, root=_REPO)
    assert ok, out
    # Correct action: renumber product onto educt -> [C, N, O] order,
    # written to a NEW file. Use the same structures the setup placed.
    from delfin.structure_edit import renumber_to_reference

    educt = _read(tmp_path, "struc/educt.xyz")
    product = _read(tmp_path, "struc/product.xyz")
    res = renumber_to_reference(product, educt)
    _write(tmp_path, "struc/product_renumbered.xyz", res["reordered_xyz"])
    ok2, out2 = run_acceptance("xe1_renumber_verify.py", tmp_path, root=_REPO)
    assert ok2, out2


def test_no_overwrite_verify_refuses_a_modified_input(tmp_path):
    ok, out = run_setup("xe1_no_overwrite_seed.py", tmp_path, root=_REPO)
    assert ok, out
    water = _read(tmp_path, "struc/water.xyz")
    # Correct: a NEW water_fixed.xyz, input untouched.
    _write(tmp_path, "struc/water_fixed.xyz", water)
    ok2, out2 = run_acceptance("xe1_no_overwrite_verify.py", tmp_path, root=_REPO)
    assert ok2, out2
    # Negative: overwriting the input must FAIL the acceptance.
    _write(tmp_path, "struc/water.xyz", water.replace("water", "ethanol"))
    ok3, out3 = run_acceptance("xe1_no_overwrite_verify.py", tmp_path, root=_REPO)
    assert not ok3, out3


def test_fix_h_setup_and_verify_grade_the_written_file(tmp_path):
    ok, out = run_setup("xe1_fix_h_seed.py", tmp_path, root=_REPO)
    assert ok, out
    # Correct action: add the missing amido H (fix_hydrogens) -> Zn1 N2 H3,
    # written to a NEW file.
    from delfin.structure_edit import fix_hydrogens

    complex_ = _read(tmp_path, "struc/complex.xyz")
    fixed = fix_hydrogens(complex_)
    _write(tmp_path, "struc/complex_fixed.xyz", fixed["xyz"])
    ok2, out2 = run_acceptance("xe1_fix_h_verify.py", tmp_path, root=_REPO)
    assert ok2, out2


def test_split_verify_grades_the_written_fragments(tmp_path):
    ok, out = run_setup("xe1_split_seed.py", tmp_path, root=_REPO)
    assert ok, out
    # Correct: each CO molecule becomes its own NEW file.
    frag_a = """2
CO fragment A
C    -3.500    0.000    0.000
O    -4.630    0.000    0.000
"""
    frag_b = """2
CO fragment B
C     3.500    0.000    0.000
O     4.630    0.000    0.000
"""
    _write(tmp_path, "struc/fragment_a.xyz", frag_a)
    _write(tmp_path, "struc/fragment_b.xyz", frag_b)
    ok2, out2 = run_acceptance("xe1_split_verify.py", tmp_path, root=_REPO)
    assert ok2, out2

