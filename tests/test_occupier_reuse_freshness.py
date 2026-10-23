"""Grounding controls for OCCUPIER result reuse (j1-grounding, L2).

An existing OCCUPIER.txt is state -- it was produced against ONE
CONTROL configuration. Reusing it silently after the configuration
changed runs all downstream chemistry on a stale electronic
structure. Each test names that: reuse must compare the stamp the
results carry against the CURRENT control file, and a missing stamp
(legacy results) is unconfirmed, not evidence.
"""

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.workflows import pipeline as pl  # noqa: E402


def _write_control(root: Path, text: str) -> Path:
    p = root / "CONTROL.txt"
    p.write_text(text, encoding="utf-8")
    return p


def _make_occ(root: Path, name: str = "initial_OCCUPIER") -> Path:
    folder = root / name
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "OCCUPIER.txt").write_text("OPT 1\n", encoding="utf-8")
    return folder


def test_matching_stamp_allows_reuse(tmp_path):
    control = _write_control(tmp_path, "method = OCCUPIER\n")
    folder = _make_occ(tmp_path)
    pl._stamp_occupier_control(control, folder)
    assert pl._occupier_results_fresh(control, folder) is True


def test_changed_control_forces_rerun(tmp_path):
    """Red before the fix: an existing OCCUPIER.txt with a control that
    changed since it was produced must NOT be reused."""
    control = _write_control(tmp_path, "method = OCCUPIER\n")
    folder = _make_occ(tmp_path)
    pl._stamp_occupier_control(control, folder)
    control.write_text("method = OCCUPIER\nbasis = def2-TZVP\n", encoding="utf-8")
    assert pl._occupier_results_fresh(control, folder) is False


def test_missing_stamp_is_unconfirmed_not_reused(tmp_path):
    """Legacy results carry no stamp: unconfirmed state, so rerun."""
    control = _write_control(tmp_path, "method = OCCUPIER\n")
    folder = _make_occ(tmp_path)
    assert pl._occupier_results_fresh(control, folder) is False


def test_corrupt_stamp_is_unconfirmed(tmp_path):
    control = _write_control(tmp_path, "method = OCCUPIER\n")
    folder = _make_occ(tmp_path)
    (folder / ".control_fingerprint").write_text("garbage\n", encoding="utf-8")
    assert pl._occupier_results_fresh(control, folder) is False


def test_missing_control_file_fails_closed(tmp_path):
    folder = _make_occ(tmp_path)
    assert pl._occupier_results_fresh(tmp_path / "nope.txt", folder) is False


def test_stamp_is_deterministic_and_content_sensitive(tmp_path):
    control = _write_control(tmp_path, "a = 1\n")
    d1 = tmp_path / "d1"
    d1.mkdir()
    pl._stamp_occupier_control(control, d1)
    stamp = (d1 / ".control_fingerprint").read_text(encoding="utf-8").strip()
    expect = hashlib.sha256(control.read_bytes()).hexdigest()
    assert stamp == expect
