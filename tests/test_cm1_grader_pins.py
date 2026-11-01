"""Red/green pins for the CM1 census grader (reviewer finding D).

The group-size check was substring containment: "6 B3LYP" passed because
"5" matched inside energy digits or job_05/job_15 tokens. The pin demands
the numbers ADJACENT to the method name (found red on 75d25391c by the
reviewer's probe tests/test_zz_cm1r2_graders.py, then fixed).
"""
from pathlib import Path

from delfin.agent.pack.benchmark.setup import cm1_run_tree_fixture as fx

ACCEPT = Path(__file__).resolve().parents[1] / "delfin/agent/pack/benchmark/accept"


def _grade(ws, text):
    doc = ws / "cm1" / "census.md"
    doc.write_text(text, encoding="utf-8")
    import subprocess
    import sys
    return subprocess.run([sys.executable, str(ACCEPT / "cm1_census_verify.py"),
                           str(ws)], capture_output=True, text=True)


def _correct_census(group_sizes_line):
    rows = "\n".join(
        f"| job_{i:02d} | ... | -402.51 | converged |" for i in range(1, 21))
    return ("20 folders tabulated.\n" + rows +
            "\n\n18 converged; job_04 failed: SCF; job_13 failed: walltime.\n"
            "Lowest within r2SCAN-3c: job_05 (-402.521 Eh). "
            "Lowest within B3LYP: job_07 (-401.781 Eh). "
            f"{group_sizes_line} Never compare energies across methods.\n")


def test_wrong_group_size_is_rejected(tmp_path):
    ws = tmp_path / "ws"
    fx.build(ws)
    doc = ws / "cm1" / "census.md"
    doc.write_text(_correct_census("Group sizes: 6 B3LYP."),
                   encoding="utf-8")
    r = _grade(ws, doc.read_text(encoding="utf-8"))
    assert r.returncode != 0
    assert "B3LYP group size" in r.stderr


def test_correct_group_size_adjacent_still_accepts(tmp_path):
    ws = tmp_path / "ws"
    fx.build(ws)
    doc = ws / "cm1" / "census.md"
    doc.write_text(_correct_census("Group sizes: 15 r2SCAN-3c, 5 B3LYP."),
                   encoding="utf-8")
    assert _grade(ws, doc.read_text(encoding="utf-8")).returncode == 0
