"""A benchmark run edited the repository it was measuring.

After a full suite the working tree held three deleted fixture files and
a MODIFIED ``tests/test_backend_local.py`` -- a block of import guarding
removed along with the comment recording why it exists.

There is already a guard: ``_PristineWorkspace`` snapshots three fixture
directories around every attempt and restores them afterwards. It works
(asserted below). It is simply not the whole story, in two ways.

Its SCOPE is three directories while the agent's workspace is the entire
checkout, so anything it touches elsewhere survives the run. That is how
a test file two directories away came to be rewritten.

And it FAILS OPEN. If the snapshot raises, ``_pairs`` is set to empty and
``__exit__`` then restores nothing, silently -- the one moment you most
need to hear about it is the one moment it says nothing.

Restoring the whole checkout automatically is not the answer: a run
would then throw away whatever the developer had uncommitted. So the
mechanism is detection. A run states plainly what it changed outside the
directories it is allowed to change, and a failed snapshot is loud.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import benchmark_runner as BR


@pytest.fixture
def repo(tmp_path):
    """A checkout-shaped tree: guarded fixture dirs plus a file elsewhere."""
    ws = tmp_path / "tests" / "fixtures" / "behavior_workspace"
    ws.mkdir(parents=True)
    (ws / "thermo.py").write_text("original\n", encoding="utf-8")
    other = tmp_path / "tests" / "test_backend_local.py"
    other.write_text("import sys\n", encoding="utf-8")
    return tmp_path, ws, other


# ---------------------------------------------------------------------------
# The guard that exists works, and must keep working
# ---------------------------------------------------------------------------

def test_a_deleted_fixture_is_restored(repo):
    root, ws, _ = repo
    with BR._PristineWorkspace(root):
        (ws / "thermo.py").unlink()
    assert (ws / "thermo.py").read_text(encoding="utf-8") == "original\n"


def test_an_edited_fixture_is_restored(repo):
    root, ws, _ = repo
    with BR._PristineWorkspace(root):
        (ws / "thermo.py").write_text("clobbered\n", encoding="utf-8")
    assert (ws / "thermo.py").read_text(encoding="utf-8") == "original\n"


def test_a_file_the_agent_invented_is_removed(repo):
    root, ws, _ = repo
    with BR._PristineWorkspace(root):
        (ws / "invented.txt").write_text("junk", encoding="utf-8")
    assert not (ws / "invented.txt").exists()


# ---------------------------------------------------------------------------
# ...and its limits are stated rather than discovered
# ---------------------------------------------------------------------------

def test_a_change_outside_the_guarded_dirs_is_reported(repo):
    """The guard cannot restore it; the run has to say it happened."""
    root, _, other = repo
    before = BR.checkout_fingerprint(root)
    other.write_text("import sys\n# the agent was here\n", encoding="utf-8")
    changed = BR.changed_outside_workspaces(root, before)
    assert any("test_backend_local.py" in p for p in changed)


def test_a_change_inside_the_guarded_dirs_is_not_reported(repo):
    """Those are restored, so naming them would be noise."""
    root, ws, _ = repo
    before = BR.checkout_fingerprint(root)
    (ws / "thermo.py").write_text("edited by a task\n", encoding="utf-8")
    assert BR.changed_outside_workspaces(root, before) == []


def test_a_deletion_outside_the_guarded_dirs_is_reported(repo):
    root, _, other = repo
    before = BR.checkout_fingerprint(root)
    other.unlink()
    changed = BR.changed_outside_workspaces(root, before)
    assert any("test_backend_local.py" in p for p in changed)


def test_a_new_file_outside_the_guarded_dirs_is_reported(repo):
    root, _, _ = repo
    before = BR.checkout_fingerprint(root)
    (root / "delfin_stray.py").write_text("x", encoding="utf-8")
    assert any("delfin_stray.py" in p for p in
               BR.changed_outside_workspaces(root, before))


def test_an_untouched_checkout_reports_nothing(repo):
    root, _, _ = repo
    before = BR.checkout_fingerprint(root)
    assert BR.changed_outside_workspaces(root, before) == []


# ---------------------------------------------------------------------------
# A guard that could not snapshot says so
# ---------------------------------------------------------------------------

def test_a_failed_snapshot_is_not_silent(repo, monkeypatch):
    """It used to set _pairs empty and restore nothing without a word."""
    root, ws, _ = repo

    def _boom(*_a, **_k):
        raise OSError("no space left on device")
    monkeypatch.setattr("shutil.copytree", _boom)

    guard = BR._PristineWorkspace(root)
    with guard:
        pass
    assert guard.failed, "a snapshot that did not happen must be visible"


def test_a_working_snapshot_does_not_claim_failure(repo):
    root, _, _ = repo
    guard = BR._PristineWorkspace(root)
    with guard:
        pass
    assert not guard.failed


def test_the_suite_checks_the_checkout():
    """Wired in, not merely available."""
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "benchmark_runner.py").read_text(encoding="utf-8")
    i = src.index("def run_suite(")
    assert "changed_outside_workspaces" in src[i:], (
        "the suite still cannot tell you it edited your checkout")
