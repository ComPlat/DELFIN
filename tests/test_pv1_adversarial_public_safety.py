"""Adversarial test (reviewer nacht-s27, package PV1): the provenance record that
`delfin.provenance_record` produces for any result must satisfy the PR-public rule
(PV1 target #4): no user names, hosts, or absolute home paths inside the record.

This is RED on `main` because `delfin.provenance_record` does not exist yet. It
turns green once the builder's record writer ships and writes a record through
its public path for a synthetic run. It is the guard that prevents the reuse-base
trap found in `delfin/cluster_bench/provenance.py:100` (`cbatch_code_state`
returns `{'repo': <ABS_PATH>}` — an absolute filesystem path embedding account
name + home path), which must be relativized/removed before it enters the record.

The module is imported unconditionally so the test fails cleanly (ImportError)
until the package exists -- a true red-before control, not a skip.
"""
from pathlib import Path

import pytest

import delfin.provenance_record as pr  # noqa: F401  (red until the package exists)


@pytest.fixture
def synthetic_run(tmp_path):
    """A minimal finished run: one result file + its inputs, under tmp_path."""
    run = tmp_path / "run_demo"
    run.mkdir()
    (run / "result.xyz").write_text("1\n\nH 0 0 0\n")
    (run / "input.inp").write_text("! B3LYP def2-SVP\n")
    return run


def _text_values(record: dict, _path: str = ""):
    """Yield every string value (including nested dict keys/values) that might
    carry a path or a name, flattened with their key path."""
    for k, v in record.items():
        keypath = f"{_path}.{k}" if _path else str(k)
        if isinstance(v, dict):
            yield from _text_values(v, keypath)
        elif isinstance(v, (str,)):
            yield keypath, v


def test_record_json_has_schema_version(synthetic_run):
    record = pr.write_record(result_dir=str(synthetic_run))
    assert isinstance(record, dict)
    assert "schema_version" in record, "record must carry a JSON schema version"


def test_record_contains_no_absolute_home_or_account_path(synthetic_run):
    home = str(Path.home())
    record = pr.write_record(result_dir=str(synthetic_run))
    for keypath, value in _text_values(record):
        # Absolute path under a home directory (any user's), or that contains the
        # current user's home path verbatim, is a PR-public violation.
        assert not (value.startswith("/")
                    and (value.startswith(home) or "/home/" in value)), (
            f"PR-public violation: '{keypath}' holds an absolute home path: {value!r}"
        )


def test_record_input_files_are_hashing_and_relative(synthetic_run):
    record = pr.write_record(result_dir=str(synthetic_run))
    inputs = record.get("input_files")
    assert inputs, "record must list input files (hash + path relative to the run)"
    for entry in inputs:
        assert not entry["path"].startswith("/"), (
            f"input path must be relative to the run, got absolute: {entry['path']!r}"
        )
        assert len(entry.get("sha256", "")) == 64, "input must carry a sha256 hash"
