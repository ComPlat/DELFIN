"""A snapshot that failed left the run going, unprotected and silent.

``_PristineWorkspace.__enter__`` copies the fixture workspaces so
``__exit__`` can put them back after each attempt. When the copy raises
it sets the pair list empty, and ``__exit__`` then restores nothing.
This morning that failure was made *visible* as ``failed``; nothing acted
on it, so a run whose protection could not be established carried on
looking exactly like one that was protected.

It is not hypothetical. Writing new fixture files into
``tests/fixtures/office_workspace`` while a benchmark was running made
``shutil.copytree`` race a file that appeared mid-copy. The snapshot
failed, the run continued, an attempt deleted ``buchungen.csv``, and
nothing put it back — the file was still missing when the run ended, and
two unrelated tests failed against the half-state before anyone noticed
why.

A transient race deserves a retry, not an abort. A snapshot that fails
twice does not: the guard exists to keep the user's files, and a run that
cannot keep them has no business touching them. So it retries once and
then refuses, saying which directory it could not copy.

The refusal is loud on purpose. The previous behaviour — carry on and
say nothing — is the one shape this project keeps finding and removing:
the worse the situation, the quieter the mechanism.
"""

from __future__ import annotations

import shutil

import pytest

from delfin.agent import benchmark_runner as BR


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "tests" / "fixtures" / "behavior_workspace"
    d.mkdir(parents=True)
    (d / "thermo.py").write_text("original\n", encoding="utf-8")
    return tmp_path, d


# ---------------------------------------------------------------------------
# A transient failure is retried
# ---------------------------------------------------------------------------

def test_a_snapshot_that_fails_once_is_retried(ws, monkeypatch):
    root, workspace = ws
    calls = {"n": 0}
    real = shutil.copytree

    def _flaky(src, dst, *a, **kw):
        calls["n"] += 1
        if calls["n"] == 1:
            raise OSError("file vanished mid-copy")
        return real(src, dst, *a, **kw)

    monkeypatch.setattr(shutil, "copytree", _flaky)
    with BR._PristineWorkspace(root) as guard:
        (workspace / "thermo.py").unlink()
    assert calls["n"] >= 2
    assert not guard.failed
    assert (workspace / "thermo.py").read_text(encoding="utf-8") == "original\n"


# ---------------------------------------------------------------------------
# A failure that persists stops the run
# ---------------------------------------------------------------------------

def test_a_snapshot_that_keeps_failing_refuses_to_run(ws, monkeypatch):
    root, _ = ws

    def _always(*_a, **_kw):
        raise OSError("no space left on device")

    monkeypatch.setattr(shutil, "copytree", _always)
    with pytest.raises(RuntimeError) as exc:
        with BR._PristineWorkspace(root):
            pass
    assert "snapshot" in str(exc.value).lower()


def test_the_refusal_names_the_directory(ws, monkeypatch):
    root, _ = ws
    monkeypatch.setattr(
        shutil, "copytree",
        lambda *_a, **_kw: (_ for _ in ()).throw(OSError("nope")))
    with pytest.raises(RuntimeError) as exc:
        with BR._PristineWorkspace(root):
            pass
    assert "behavior_workspace" in str(exc.value)


def test_the_refusal_says_why_it_matters(ws, monkeypatch):
    root, _ = ws
    monkeypatch.setattr(
        shutil, "copytree",
        lambda *_a, **_kw: (_ for _ in ()).throw(OSError("nope")))
    with pytest.raises(RuntimeError) as exc:
        with BR._PristineWorkspace(root):
            pass
    assert "restore" in str(exc.value).lower()


# ---------------------------------------------------------------------------
# ...without changing anything that already worked
# ---------------------------------------------------------------------------

def test_a_healthy_snapshot_restores_as_before(ws):
    root, workspace = ws
    with BR._PristineWorkspace(root) as guard:
        (workspace / "thermo.py").write_text("clobbered\n", encoding="utf-8")
        (workspace / "invented.txt").write_text("junk", encoding="utf-8")
    assert not guard.failed
    assert (workspace / "thermo.py").read_text(encoding="utf-8") == "original\n"
    assert not (workspace / "invented.txt").exists()


def test_a_root_with_no_fixture_directories_is_not_an_error(tmp_path):
    """Running the suite from somewhere else must not raise."""
    with BR._PristineWorkspace(tmp_path) as guard:
        pass
    assert not guard.failed
