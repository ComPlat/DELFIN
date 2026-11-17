"""Tests for delfin.ts_goat_traj: orca_pltvib multi-frame XYZ for the
TS -> GOAT panel.

The unit tests build a small, synthetic Hessian and do not rely on
fixtures removed from main. Real orca_pltvib is tested only when its
binary is available (skipped otherwise).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

# Anchor the import to THIS checkout so a pytest run started from a
# different working directory cannot silently pick up an installed
# `delfin` that lacks the new module (observed: the auto-verifier ran
# from outside the worktree and imported /pfs/.../software/delfin).
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from delfin import ts_goat_traj  # noqa: E402  (path anchored above)

@pytest.fixture
def hessian(tmp_path):
    """Synthetic ORCA Hessian with a real 3N-mode table and known frequency."""
    from tests.test_ts_goat_reactive_center import _write_hess
    atoms = [('O', 0.0, 0.0, 0.0), ('H', 0.96, 0.0, 0.0),
             ('H', -0.3, 0.8, 0.0)]
    zero = [(0.0, 0.0, 0.0)] * len(atoms)
    mode = [(0.1, 0, 0), (0, 0, 0), (0, 0, 0)]
    hess = tmp_path / 'm.hess'
    _write_hess(hess, atoms, [zero] * 6 + [mode] + [zero] * 2,
                [100.0] * 6 + [-21.86] + [100.0] * 2)
    return hess


MULTIFRAME = """3
frame a
O 0.0 0.0 0.0
H 0.9 0.0 0.0
H -0.3 0.8 0.0
3
frame b
O 0.0 0.0 0.1
H 0.9 0.0 0.0
H -0.3 0.8 0.1
3
frame c
O 0.0 0.0 -0.1
H 0.9 0.0 0.1
H -0.3 0.8 0.0
"""


@pytest.fixture
def multiframe_xyz(tmp_path):
    p = tmp_path / 'traj.xyz'
    p.write_text(MULTIFRAME)
    return p


def test_pltvib_output_path_mapping():
    """orca_pltvib writes <hess>.v<NNN>.xyz for a 1-based mode index."""
    hess = Path('run.hess')
    assert ts_goat_traj.pltvib_output_path(hess, 6) == Path(
        'run.hess.v006.xyz')
    assert ts_goat_traj.pltvib_output_path(hess, 1) == Path(
        'run.hess.v001.xyz')


def test_parse_multiframe_xyz_keeps_atom_order(multiframe_xyz):
    frames = ts_goat_traj.parse_multiframe_xyz(multiframe_xyz)
    assert len(frames) == 3
    order = [[el for el, *_ in f] for f in frames]
    assert all(o == order[0] for o in order), (
        'element order must be identical across frames')


def test_parse_multiframe_xyz_rejects_mixed_order(tmp_path):
    lines = MULTIFRAME.replace(
        'frame c\nO 0.0 0.0 -0.1', 'frame c\nH 0.0 0.0 -0.1').splitlines()
    p = tmp_path / 'bad.xyz'
    p.write_text('\n'.join(lines))
    with pytest.raises(ValueError):
        ts_goat_traj.parse_multiframe_xyz(p)


def test_parse_multiframe_xyz_rejects_count_mismatch(tmp_path):
    lines = MULTIFRAME.splitlines()
    # third frame header says 3 but holds only 2 atom lines
    p = tmp_path / 'short.xyz'
    p.write_text('\n'.join(lines[:-2]))
    with pytest.raises(ValueError):
        ts_goat_traj.parse_multiframe_xyz(p)


def test_scaled_is_linear_and_preserves_order():
    traj = ts_goat_traj.VibrationTrajectory(
        mode_index=6, frequency_cm1=-100.0,
        frames=[
            [('O', 0.0, 0.0, 0.0), ('H', 1.0, 0.0, 0.0)],
            [('O', 0.0, 0.0, 0.5), ('H', 1.0, 0.0, 0.0)],
            [('O', 0.0, 0.0, -0.5), ('H', 1.0, 0.0, 0.0)],
        ])
    out = traj.scaled(2.0)
    # largest displacement must equal the requested amplitude
    biggest = max(abs(f[0][3]) for f in out)
    assert abs(biggest - 2.0) < 1e-9
    # order preserved
    assert [[el for el, *_ in f] for f in out] == [['O', 'H']] * 3
    # zero amplitude collapses to the equilibrium frame
    zero = traj.scaled(0.0)
    assert all(abs(f[0][3]) < 1e-12 for f in zero)


def test_frames_to_xyz_roundtrip(multiframe_xyz, tmp_path):
    frames = ts_goat_traj.parse_multiframe_xyz(multiframe_xyz)
    text = ts_goat_traj.frames_to_xyz(frames)
    out = tmp_path / 'round.xyz'
    out.write_text(text)
    again = ts_goat_traj.parse_multiframe_xyz(out)
    assert [[(el, round(x, 4), round(y, 4), round(z, 4))
             for el, x, y, z in f] for f in again] == [
        [(el, round(x, 4), round(y, 4), round(z, 4))
         for el, x, y, z in f] for f in frames]


def test_mode_outside_hessian_rejected():
    result = ts_goat_traj.generate_vibration_trajectory(
        Path('nonexistent.hess'), 1)
    assert result.trajectory is None and result.error


def test_pltvib_unavailable_fails_cleanly(hessian, tmp_path):
    """Missing orca_pltvib -> visible error, no rendered trajectory."""
    result = ts_goat_traj.generate_vibration_trajectory(
        hessian, 6, cache_dir=tmp_path / 'cache',
        pltvib='orca_pltvib-definitely-missing')
    assert result.trajectory is None
    assert result.error and 'orca_pltvib' in result.error.lower()


def _seed_cache(hessian, cache_dir):
    cache_dir.mkdir(parents=True, exist_ok=True)
    local = cache_dir / hessian.name
    shutil.copy2(hessian, local)
    xyz = ts_goat_traj.pltvib_output_path(local, 6)
    xyz.write_text(
        '3\nframe 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.3 0.8 0.0\n'
        '3\nframe 2\nO 0.1 0.0 0.0\nH 0.96 0.0 0.0\nH -0.3 0.8 0.0\n',
        encoding='utf-8')
    return xyz


def test_cache_reuse_returns_hessian_frequency(hessian, tmp_path):
    """A real cached XYZ is reused with the current Hessian's frequency."""
    cached = _seed_cache(hessian, tmp_path / 'cache')
    result = ts_goat_traj.generate_vibration_trajectory(
        hessian, 6, cache_dir=cached.parent, pltvib='not-needed')
    assert result.trajectory is not None and result.reused
    assert result.trajectory.frequency_cm1 == pytest.approx(-21.86)


def test_cache_check_runs_without_hashlib_file_digest(
        hessian, tmp_path, monkeypatch):
    """hashlib.file_digest is 3.11+; requires-python admits 3.10."""
    import hashlib
    monkeypatch.delattr(hashlib, 'file_digest', raising=False)
    cached = _seed_cache(hessian, tmp_path / 'cache')
    result = ts_goat_traj.generate_vibration_trajectory(
        hessian, 6, cache_dir=cached.parent, pltvib='not-needed')
    assert result.trajectory is not None and result.reused


def test_changed_hessian_invalidates_cached_mode(hessian, tmp_path):
    """Same-name Hessian with changed modes must discard old XYZ frames."""
    cached = _seed_cache(hessian, tmp_path / 'cache')
    assert cached.is_file()
    # Changing a frequency alone must invalidate the cached trajectory.
    text = hessian.read_text(encoding='utf-8')
    hessian.write_text(text.replace('-21.86', '-31.86'), encoding='utf-8')
    result = ts_goat_traj.generate_vibration_trajectory(
        hessian, 6, cache_dir=cached.parent, pltvib='orca_pltvib-not-installed')
    assert result.trajectory is None and result.error
    assert not cached.exists()


def test_real_pltvib_against_committed_fixture(hessian):
    """Integration gate for real ORCA, if a real compatible Hessian is present.

    Synthetic ORCA Hessians are for unit tests only and may not be accepted
    by the external ORCA utility.
    """
    pltvib = shutil.which('orca_pltvib')
    if pltvib is None:
        pytest.skip('orca_pltvib not installed on this host')
    fixture = _REPO_ROOT / (
        'tests/fixtures/orca611/sh1_cu_nh3_imagmode/'
        'sh1_cu_nh3_imagmode.hess'
    )
    if not fixture.is_file():
        pytest.skip('optional real ORCA Hessian fixture not installed')
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        hess = Path(td) / 'fixture.hess'
        shutil.copy(fixture, hess)
        result = ts_goat_traj.generate_vibration_trajectory(
            hess, 6, cache_dir=Path(td) / 'cache', pltvib=pltvib)
        assert result.error is None, result.error
        traj = result.trajectory
        assert traj.n_frames >= 2
        order = [[el for el, *_ in f] for f in traj.frames]
        assert all(o == order[0] for o in order)
