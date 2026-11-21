"""A test that binds a Unix socket is not limited by the shell's TMPDIR.

Input: a pytest run started from a shell whose ``TMPDIR`` is deep (a DELFIN
agent session exports its scratch there), so ``tmp_path`` for a test is that
much longer. Output: whether a socket made for the test binds.

Semantics: ``bind`` on an ``AF_UNIX`` path is capped near 107 bytes by
``sockaddr_un``. ``TMPDIR`` is a property of the shell that started the
suite, not of the code under test, so a checkout under a long temporary root
fails tests that would bind a socket anywhere while CI -- whose ``TMPDIR`` is
``/tmp`` -- never sees it. The suite's answer is the one
``_child_runtime_dir`` and the dashboard-kernel tests already use: a socket
path is made in a short private directory, and it is asked for rather than
assumed from the depth of a name.

Measured 2026-11-20 on a workstation whose ``TMPDIR`` was
``/tmp/delfin-cmd-tmp-<id>/scratch-<id>`` (64 characters): three cases of
``test_the_bridge_admits_its_own_uid_only_for_the_gateway`` failed with
``OSError: AF_UNIX path too long`` on a 110-character path, and two of the
three passed with ``TMPDIR=/tmp``.

The alternative rejected: shortening the test's own directory name, or
pointing pytest's ``basetemp`` at ``/tmp``. The first fits one test name and
overflows on the next; the second was measured and is worse than the bug --
pytest removes its base directory at the end of the run, so it tried to
remove ``/tmp`` itself, and only the read-only mount stopped it.

The module only measures the predicate; the fix lives where a socket is
bound (``tests/_messenger_people_rig.py``).
"""
from __future__ import annotations

import os
import shutil
import socket
import tempfile
from pathlib import Path

import pytest

from conftest import child_env


def _binding_fits(where: Path) -> bool:
    """Whether a socket path a test builds under ``where`` will bind."""
    path = str(where / "b.sock")
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        server.bind(path)
        return True
    except OSError:
        return False
    finally:
        server.close()
        try:
            os.unlink(path)
        except OSError:
            pass


def test_a_deep_tmp_path_does_not_fit_a_socket():
    """The premise: pytest's ``tmp_path`` under a deep TMPDIR is too long.

    Built here rather than taken from the environment, so it holds on every
    host: the length is what the limit reacts to, not this machine's TMPDIR.
    """
    with tempfile.TemporaryDirectory(prefix="dx" + "d" * 60 + "-") as deep:
        where = Path(deep) / ("t" * 40)
        where.mkdir()
        assert len(str(where / "b.sock").encode()) > 107
        assert not _binding_fits(where)


def test_a_short_private_directory_fits_a_socket():
    """The answer the suite uses: a short private directory binds."""
    root = "/tmp" if os.path.isdir("/tmp") else tempfile.gettempdir()
    made = tempfile.mkdtemp(prefix="ds", dir=root)
    try:
        assert _binding_fits(Path(made))
    finally:
        shutil.rmtree(made, ignore_errors=True)


def test_a_child_gets_a_runtime_that_fits_a_socket(tmp_path):
    """``child_env`` hands a child a runtime directory short enough for a
    socket, whatever TMPDIR the suite was started with."""
    env = child_env(tmp_path)
    made = None
    for key in ("JUPYTER_RUNTIME_DIR", "XDG_RUNTIME_DIR"):
        value = env.get(key)
        if value:
            made = Path(value)
            break
    if made is None:
        pytest.skip("no runtime directory is handed to the child")
    assert len(str(made / "b.sock").encode()) < 108, str(made)
