"""CAS1 phase 3b: the benchmark tasks grade on physics.

Red before tasks_w25_cas1.yaml / accept/cas1_verify.py exist. The
acceptance must re-run the module's checks on the agent's structured
deliverable (proposal.json), refuse a bad proposal, and grade the trap
task on refusal -- not on prose.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TASKS = ROOT / "delfin/agent/pack/benchmark/tasks_w25_cas1.yaml"
VERIFY = ROOT / "delfin/agent/pack/benchmark/accept/cas1_verify.py"


def _task_ids():
    data = yaml.safe_load(TASKS.read_text(encoding="utf-8"))
    return [t["id"] for t in data["tasks"]]


def test_tasks_file_loads_with_three_tasks():
    ids = _task_ids()
    assert ids == ["cas1_propose_ozone", "cas1_propose_n2_stretched",
                   "cas1_trap_oversized"]


def _run_verify(ws: Path, task_id: str) -> int:
    # The accept script imports delfin.casscf_active_space. The gate
    # provides that environment; a bare `python -m pytest` outside it
    # does not, so the repo root is put on PYTHONPATH explicitly.
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, str(VERIFY), str(ws), task_id],
        capture_output=True, text=True, env=env).returncode


def _workspace(tmp_path: Path, task_id: str, files: dict) -> Path:
    d = tmp_path / task_id
    d.mkdir(parents=True)
    for name, text in files.items():
        (d / name).write_text(text, encoding="utf-8")
    return tmp_path


def test_good_proposal_is_accepted(tmp_path):
    ws = _workspace(tmp_path, "cas1_propose_ozone", {
        "proposal.json": json.dumps({
            "n_electrons": 6, "n_orbitals": 6,
            "justification": "full pi space of the delocalised ozone "
                             "pi system, bonding and antibonding",
            "named_orbital_classes": ["pi bonding", "pi* antibonding"],
        }),
        "answer.txt": "the pi system is multireference",
    })
    assert _run_verify(ws, "cas1_propose_ozone") == 0


def test_bad_proposal_is_refused(tmp_path):
    ws = _workspace(tmp_path, "cas1_propose_n2_stretched", {
        "proposal.json": json.dumps({
            "n_electrons": 4, "n_orbitals": 2,
            "justification": "",
            "named_orbital_classes": [],
        }),
        "answer.txt": "trust me",
    })
    assert _run_verify(ws, "cas1_propose_n2_stretched") != 0


def test_pinned_occupations_would_not_validate(tmp_path):
    # A space that is valid in size but whose justification carries no
    # chemistry is refused — the grading is on the physics, not prose.
    ws = _workspace(tmp_path, "cas1_propose_ozone", {
        "proposal.json": json.dumps({
            "n_electrons": 6, "n_orbitals": 6,
            "justification": "the numbers look nice",
            "named_orbital_classes": ["everything"],
        }),
        "answer.txt": "ok",
    })
    assert _run_verify(ws, "cas1_propose_ozone") != 0


def test_trap_task_requires_refusal(tmp_path):
    ws = _workspace(tmp_path / "bad", "cas1_trap_oversized", {
        "answer.txt": "Nactel 24 and Nactorb 24 are real keywords, "
                      "the space is fine.",
    })
    assert _run_verify(ws, "cas1_trap_oversized") != 0
    ws2 = _workspace(tmp_path / "good", "cas1_trap_oversized", {
        "answer.txt": "Nactel/Nactorb do not exist; the real keys are "
                      "nel and norb (manual 3.13.3). 24 orbitals exceed "
                      "CASSCF; use RAS, ICE or DMRG instead.",
    })
    assert _run_verify(ws2, "cas1_trap_oversized") == 0
