"""`2>/dev/null` was read as reaching outside the workspace.

A locked-scope session refuses any command naming a path outside its
folder, and it resolves symlinks first so a link inside the folder
pointing out of it cannot be used as a door. /dev/null is a link, so
`find . -name "*.csv" 2>/dev/null` was refused as an escape.

That is the commonest idiom in shell there is. In the audit log it is 58
refusals, most of them a find or an ls with its errors silenced, and
nothing whatever is contained by refusing them: /dev/null discards what
it is given and yields nothing when read.

The exemption is a list of devices, NOT a "/dev/" prefix. A prefix would
also exempt /dev/sda1 and /dev/mem, which is precisely the escape the
check exists for. Every name on the list carries no user data in either
direction, or aliases a descriptor the process already holds.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor, _is_safe_device,
)

_WS = Path("tests/fixtures/office_workspace").resolve()


@pytest.fixture
def locked():
    if not _WS.is_dir():
        pytest.skip("the office fixture is not in this checkout")
    perms = KitToolPermissions(workspace=str(_WS), lock_workspace=True)
    perms.mode = "acceptEdits"
    perms.task_session_id = "dev-null"
    assert perms.scope_locked, "premise: this test needs a locked scope"
    return perms


def _refused(perms, command: str) -> bool:
    out = _DocToolExecutor().execute("bash", {"command": command}, perms)
    return '"error"' in out[:40]


# ---------------------------------------------------------------------------
# What must now run
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("command", [
    'find . -name "*.csv" 2>/dev/null',
    'ls nonexistent 2>/dev/null || echo none',
    'grep -r Betrag . 2>/dev/null | head -3',
    'ls . >/dev/null 2>&1 || echo none',
    'cat buchungen.csv 2>/dev/null | head -2',
])
def test_silencing_errors_is_not_an_escape(locked, command):
    assert not _refused(locked, command), command


# ---------------------------------------------------------------------------
# What must still be refused
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("command", [
    "cat /dev/sda1",
    "cat /dev/mem",
    "dd if=/dev/sda of=copy.img",
    "cat /etc/passwd",
    "cat ../rechnungen.csv",
])
def test_a_real_device_or_a_real_escape_is_still_refused(locked, command):
    assert _refused(locked, command), command


# ---------------------------------------------------------------------------
# The predicate
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", [
    "/dev/null", "/dev/zero", "/dev/urandom", "/dev/random", "/dev/full",
    "/dev/tty", "/dev/stdin", "/dev/stdout", "/dev/stderr", "/dev/fd/3",
])
def test_the_pseudo_devices_are_safe(path):
    assert _is_safe_device(path)


@pytest.mark.parametrize("path", [
    "/dev/sda", "/dev/sda1", "/dev/mem", "/dev/kmem", "/dev/nvme0n1",
    "/dev/disk/by-uuid/x", "/dev", "/dev/shm/secret",
])
def test_a_real_device_is_not_on_the_list(path):
    """A "/dev/" prefix would have exempted all of these, which is the
    escape the containment check is for."""
    assert not _is_safe_device(path)


def test_the_list_is_names_not_a_prefix():
    """Stated as a test because the tempting fix is one line and wrong."""
    from delfin.agent.api_client import _SAFE_DEVICE_PREFIXES

    assert "/dev/" not in _SAFE_DEVICE_PREFIXES
    assert all(p.startswith("/dev/") and len(p) > len("/dev/")
               for p in _SAFE_DEVICE_PREFIXES)
