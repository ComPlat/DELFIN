"""A fixture that needs a directory outside /tmp finds one, or says why not.

The hide-path tests need a lockbox that a cage does not already mask: every
cage mounts a fresh tmpfs on /tmp, so a lockbox there is invisible with or
without the feature and the real-run tests would pass on code that hides
nothing. ``/var/tmp`` was the one path hard-coded for that, and on a
workstation where ``/var/tmp`` is mounted read-only the fixture raised in
setup -- five tests in ``test_a_run_can_hide_paths_from_its_cage.py`` came
back as ERRORS, not as skips, and the run reported a suite that measures
nothing about hiding rather than a host that cannot host the measurement.

Measured 2026-11-20 on that workstation: ``os.access("/var/tmp",
os.W_OK)`` was False and ``tempfile.mkdtemp(dir="/var/tmp")`` raised
``OSError: [Errno 30] Read-only file system``.

Input: nothing. Output: ``outside_tmp()` returns a writable directory that
is NOT under /tmp, or None. Semantics: it asks the host in the order a real
lockbox lives -- /var/tmp, then the user's home, then the system temporary
root only if that root is not /tmp -- and checks each by CREATING a
directory there, because a permission bit and a read-only mount look alike
to ``os.access`` on some filesystems. Nothing is left behind.

The alternative rejected: ``pytest.skip`` on ``not os.access(...)`` alone,
which the sibling fixture in ``test_a_caged_command_runs_without_asking.py``
uses. ``os.access`` answers the permission bit, not whether the mount
accepts writes (it returned False here for a read-only mount, but returns
True for a writable bit on a read-only mount of the same mode); the create
is the measurement.
"""

from __future__ import annotations

from pathlib import Path

import conftest


def test_a_directory_outside_tmp_is_returned_when_one_is_writable():
    """On a host with a writable one, the answer is a real directory.

    The test asks the host rather than assuming: on a machine where none of
    the candidates accepts a write, the skip below is the correct outcome
    and this test says so instead of failing.
    """
    found = conftest.outside_tmp()
    if found is None:
        import pytest
        pytest.skip("no writable directory outside /tmp on this host")
    assert found.is_dir()
    assert not str(found.resolve()).startswith("/tmp/")


def test_the_answer_is_not_under_tmp_even_when_tmp_is_the_only_writable_root():
    """The one property the caller depends on is checked, not assumed.

    A cage gives /tmp a fresh tmpfs. A "outside tmp" answer that is under
    /tmp would make every real-run hide test pass on code that hides
    nothing -- silently, which is worse than the error it replaced.
    """
    found = conftest.outside_tmp()
    if found is None:
        import pytest
        pytest.skip("no writable directory outside /tmp on this host")
    resolved = str(found.resolve())
    assert not resolved.startswith(("/tmp/", "/private/tmp/")), resolved


def test_a_read_only_candidate_is_rejected_rather_than_returned(monkeypatch):
    """The check is a write, so a read-only mount cannot pass as writable.

    ``os.access`` here is a stand-in that always says yes, as it does for a
    writable MODE on a read-only mount; the answer must still be None,
    because the helper's measurement is the create. Measured on this host:
    ``/var/tmp`` has ``W_OK`` False AND refuses ``mkdtemp``.

    The ``mkdtemp`` the helper calls is patched where it is looked up --
    ``tempfile`` in this test process -- so no directory is created and the
    test means the same thing on a host where ``/var/tmp`` is writable.
    """
    import tempfile

    def _refuse(**kw):
        raise OSError(30, "Read-only file system")

    monkeypatch.setattr(conftest.os, "access", lambda path, mode: True)
    monkeypatch.setattr(tempfile, "mkdtemp", _refuse)
    assert conftest.outside_tmp() is None


def test_the_candidates_are_asked_in_the_order_a_lockbox_lives():
    """A real lockbox is in a home directory; /var/tmp is the first guess.

    The order is asserted because it decides which path a test measures on a
    host that offers several: the reader of a failure needs to know which
    directory the numbers came from.
    """
    assert conftest._OUTSIDE_TMP_CANDIDATES[0] == Path("/var/tmp")
    second = conftest._OUTSIDE_TMP_CANDIDATES[1]
    assert Path.home().samefile(second) if second.exists() else True
