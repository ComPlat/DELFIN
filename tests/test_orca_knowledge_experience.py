"""The ORCA experience extractor: records, healings, privacy, CLI walk."""

import json
from pathlib import Path

import pytest

from delfin.orca_knowledge.experience import (
    extract, link, system_class, _walk_records,
)

FIX = Path(__file__).parent / "fixtures" / "orca_experience"


def test_normal_run_record():
    r = extract(FIX / "normal" / "run1.inp", FIX / "normal" / "run1.out")
    assert r["orca_version"] == "6.1.1"
    assert "PBE0" in r["simple_input_keywords"]
    assert "FREQ" in r["simple_input_keywords"]
    assert r["charge"] == 0 and r["multiplicity"] == 1
    assert r["element_counts"] == {"O": 1, "H": 2}
    assert r["job_types"] == ["freq"]
    assert r["method_family"] == "dft"
    assert r["normal_termination"] is True
    assert r["scf_iterations"] == 12 and r["scf_converged"] is True
    assert r["imaginary_frequencies"] == 0
    assert r["wall_time_seconds"] == pytest.approx(123.45)
    assert r["error_signatures"] == []
    assert "dir_hash" in r and len(r["dir_hash"]) == 16


def test_scf_failure_record():
    r = extract(FIX / "scf_failure" / "run2.inp",
                FIX / "scf_failure" / "run2.out")
    assert r["normal_termination"] is False
    assert r["scf_converged"] is False and r["scf_iterations"] == 124
    assert "SCF NOT CONVERGED" in r["error_signatures"]
    assert "ORCA finished by error termination in LEANSCF" in r["error_signatures"]
    assert r["method_family"] == "dft"
    assert r["job_types"] == ["opt"]
    sc = r["system_class"]
    assert sc["metals"] == [{"element": "Fe", "count": 1}]
    assert "Fe" in sc["transition_metals"]
    assert sc["total_atoms"] == 3


def test_no_coordinates_or_paths_in_record():
    r = extract(FIX / "scf_failure" / "run2.inp",
                FIX / "scf_failure" / "run2.out")
    text = json.dumps(r)
    assert "Fe 0.0" not in text
    for forbidden in ("0.0 0.0", str(FIX), "scf_failure"):
        assert forbidden not in text


def test_link_pairs_failure_with_later_success():
    failed = extract(FIX / "scf_failure" / "run2.inp",
                     FIX / "scf_failure" / "run2.out")
    healed = extract(FIX / "slowconv_pair" / "run3.inp",
                     FIX / "slowconv_pair" / "run3.out")
    # same directory is required by link(); emulate by overriding the hash
    healed = dict(healed, dir_hash=failed["dir_hash"])
    heals = link([failed, healed])
    assert len(heals) == 1
    h = heals[0]
    assert h["error_signatures"] == failed["error_signatures"]
    assert "+SLOWCONV" in h["keyword_changes"]
    assert h["healed_run"] and h["failed_run"] and h["healed_run"] != h["failed_run"]


def test_link_ignores_two_successes():
    ok1 = extract(FIX / "normal" / "run1.inp", FIX / "normal" / "run1.out")
    assert link([ok1, dict(ok1)]) == []


def test_walk_records_walks_fixture_tree(tmp_path):
    root = tmp_path / "root"
    (root / "a").mkdir(parents=True)
    for name in ("run1.inp", "run1.out"):
        (root / "a" / name).write_bytes((FIX / "normal" / name).read_bytes())
    # a symlink pointing outside must not be followed
    (root / "b").symlink_to(FIX.parent.parent.parent)
    records = _walk_records(root)
    assert len(records) == 1
    # File names can carry molecule names: a record holds a hash, not a name.
    assert "run1" not in json.dumps(records[0]) and records[0]["run_hash"]


def test_system_class_summary():
    r = extract(FIX / "scf_failure" / "run2.inp",
                FIX / "scf_failure" / "run2.out")
    sc = system_class(r)
    assert sc["total_charge"] == 0 and sc["multiplicity"] == 2
    assert sc["total_atoms"] == 3
    assert sc["light_elements"] == ["N", "O"]

import re
from pathlib import Path



FIX = Path(__file__).parent / "fixtures" / "orca_experience"
MULTI_JOB = FIX / "multi_job"


def _record_json(record) -> str:
    return json.dumps(record)


def test_multi_job_blocks_stop_at_end_and_new_job():
    """A block ends at its 'end'; a $new_job never leaks into it."""
    r = extract(MULTI_JOB / "run4.inp", MULTI_JOB / "run4.out")
    # same-named blocks in different jobs merge their option lines; the
    # first %maxcore (job 1) must not swallow $new_job or job 2's options
    assert r["blocks"].get("TDDFT") == ["nroots 5", "nroots 3"]
    assert r["blocks"].get("MAXCORE") == ["4000"]
    assert "OPT" in r["simple_input_keywords"]


def test_multi_job_element_counts_are_real_elements_only():
    """Geometry atoms only -- never input tokens like '!', 'end', 'nroots'."""
    r = extract(MULTI_JOB / "run4.inp", MULTI_JOB / "run4.out")
    assert set(r["element_counts"]) <= {"H", "O", "N"}
    assert "$new_job" not in r["element_counts"]
    assert "end" not in r["element_counts"]
    assert "nroots" not in r["element_counts"]


def test_multi_job_no_coordinates_in_record_json():
    """Privacy: no 4+ decimal number and no geometry token anywhere."""
    r = extract(MULTI_JOB / "run4.inp", MULTI_JOB / "run4.out")
    text = _record_json(r)
    assert not re.search(r"\d\.\d{4,}", text), "coordinate-like number leaked"
    for token in ("5.57510392", "1.26479274", "1.99310029", "6.10000000",
                  "5.575103", "$new_job"):
        assert token not in text, f"geometry token {token!r} leaked"


def test_multi_job_second_job_geometry_counted():
    """Elements from EVERY job's geometry are counted."""
    r = extract(MULTI_JOB / "run4.inp", MULTI_JOB / "run4.out")
    # each job's geometry contributes: job1 H+O, job2 H+O, job3 N
    assert sum(r["element_counts"].values()) == 5
