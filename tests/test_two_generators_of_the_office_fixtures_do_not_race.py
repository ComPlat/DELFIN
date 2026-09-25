"""Two generators racing on one office workspace.

The parallel suite run (SLURM 7188718) starts one pytest process per
test file, so two processes can reach for the office workbooks at the
same moment. This file checks the property the generation already has
and must keep: whatever two concurrent generators leave behind, every
workbook is either wholly old, wholly new, or -- because two builds of
one spec are byte-identical -- wholly the same. No half-written file,
no mixed zip, no reader that fails.

The write path earned that property on purpose: ``_write_workbook``
builds each file under a ``.wb-*`` temp name and ``os.replace``s it
into place (delfin/agent/benchmark_fixtures.py), and
``_make_deterministic`` strips the clock from the zip, so two builds
of one spec produce identical bytes. What the incident ACTUALLY
reported were tracked CSVs a neighbour's restore had momentarily
removed -- a different window, covered by the confirmation scan in
conftest. This file pins the generator half so the race stays closed.
"""

from __future__ import annotations

import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

from conftest import child_env

from delfin.agent import benchmark_fixtures as BF

_WORKER = r"""
import sys
from pathlib import Path
from delfin.agent import benchmark_fixtures as BF
root = Path(sys.argv[1])
written, reason = BF.ensure_office_fixtures(root)
print("OK" if not reason else "REASON:" + reason, flush=True)
"""


def _repo_like_checkout(tmp_path: Path) -> Path:
    """A minimal checkout root: the office workspace directory with the
    tracked files git carries, so the generator finds a home for the
    workbooks."""
    root = tmp_path / "repo"
    ws = root / BF.OFFICE_WS_REL
    ws.mkdir(parents=True)
    (ws / ".gitignore").write_text("*.xlsx\n.fixture-stamp\n",
                                   encoding="utf-8")
    (ws / "README.md").write_text("fixtures\n", encoding="utf-8")
    for name in ("buchungen.csv", "inventar.csv",
                 "kostenstellen_roh.csv", "rechnungen.csv"):
        (ws / name).write_text("a;b\n", encoding="utf-8")
    return root


def _spawn_two_generators(root: Path, tmp_path, child_env_fixture):
    """Two child processes building the same fixtures at the same time.

    Children, not threads: the suite runs as one process per file, so
    the race is between processes, and a thread would share the module
    state the processes do not. ``-I`` would ignore PYTHONPATH, and
    without it a child run from tmp_path cannot import the checkout's
    delfin package.
    """
    env = dict(child_env_fixture)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
    script = tmp_path / "worker.py"
    script.write_text(_WORKER, encoding="utf-8")
    procs = [
        subprocess.Popen(
            [sys.executable, str(script), str(root)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env=env)
        for _ in range(2)
    ]
    outs = [p.communicate(timeout=240) for p in procs]
    return procs, outs


@pytest.fixture
def child_env_fixture(tmp_path):
    return child_env(tmp_path)


def test_two_generators_leave_whole_workbooks(tmp_path,
                                              child_env_fixture):
    """Both children run ``ensure_office_fixtures`` from empty; both
    must finish saying OK (or a named dependency reason), and every
    workbook on disk must be a complete zip that a reader can open --
    no half-written product of the race."""
    if BF.missing_dependency_reason():
        pytest.skip(BF.missing_dependency_reason())
    root = _repo_like_checkout(tmp_path)
    procs, outs = _spawn_two_generators(root, tmp_path, child_env_fixture)
    for p, (out, _err) in zip(procs, outs):
        assert p.returncode == 0, out
        assert out.strip().endswith("OK"), out
    ws = root / BF.OFFICE_WS_REL
    for name in BF.WORKBOOK_NAMES:
        target = ws / name
        assert target.is_file(), f"{name} missing after the race"
        # A whole zip: the central directory parses and every entry
        # reads. A file a generator was still writing when the other
        # opened it would not pass this.
        with zipfile.ZipFile(target) as zf:
            assert zf.testzip() is None, f"{name} is corrupt"
            assert zf.namelist(), f"{name} has no entries"


def test_two_generators_of_one_spec_leave_identical_bytes(
        tmp_path, child_env_fixture):
    """The other half of the race property: two builds of one spec are
    byte-identical, so whichever generator won a given file, a reader
    sees the same workbook. Without that, "wholly old or wholly new"
    would still mean two different files in flight."""
    if BF.missing_dependency_reason():
        pytest.skip(BF.missing_dependency_reason())
    # Two separate roots, built independently by one process each --
    # the same comparison the race resolves by identity rather than by
    # locking.
    digests = []
    roots = []
    for i in range(2):
        root = _repo_like_checkout(tmp_path / f"repo{i}")
        roots.append(root)
        script = tmp_path / "worker.py"
        script.write_text(_WORKER, encoding="utf-8")
        env = dict(child_env_fixture)
        env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
        done = subprocess.run(
            [sys.executable, str(script), str(root)],
            capture_output=True, text=True, timeout=240,
            env=env)
        assert done.returncode == 0, done.stdout + done.stderr
        assert done.stdout.strip().endswith("OK"), done.stdout
        ws = root / BF.OFFICE_WS_REL
        digests.append({name: (ws / name).read_bytes()
                        for name in BF.WORKBOOK_NAMES})
    assert digests[0] == digests[1], (
        "two builds of one spec differ -- the deterministic write is "
        "broken and the race has no identity to fall back on")
    # ... and the stamp matches the spec both children built.
    ws0 = roots[0] / BF.OFFICE_WS_REL
    stamp = ws0 / BF.STAMP_NAME
    assert BF.fixtures_are_current(ws0), (
        "the workspace the children built does not match the spec: "
        f"files={sorted(p.name for p in ws0.iterdir())}, "
        f"stamp={stamp.read_text(encoding='utf-8').strip() if stamp.is_file() else None!r}, "
        f"spec={BF.spec_digest()!r}")


def test_no_temporary_workbook_survives_the_race(tmp_path,
                                                 child_env_fixture):
    """The atomic replace must leave no ``.wb-*`` leftovers beside the
    workbooks: a temp file the race dropped is a path the checkout
    guard would rightly report."""
    if BF.missing_dependency_reason():
        pytest.skip(BF.missing_dependency_reason())
    root = _repo_like_checkout(tmp_path)
    procs, outs = _spawn_two_generators(root, tmp_path, child_env_fixture)
    for p, (out, _err) in zip(procs, outs):
        assert p.returncode == 0, out
    leftovers = [p.name for p in (root / BF.OFFICE_WS_REL).iterdir()
                 if p.name.startswith(".wb-")]
    assert not leftovers, f"race left temp workbooks: {leftovers}"
