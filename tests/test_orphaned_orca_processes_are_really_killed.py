"""The cleanup that removes orphaned ORCA workers has to survive its own path.

ORCA leaves MPI workers behind when it is killed, and
``_ensure_process_group_terminated`` exists to escalate SIGTERM to SIGKILL
until they are gone. On a machine without psutil that escalation raised
UnboundLocalError on its own ``time.sleep``: a second ``import time`` further
down the same function made ``time`` a local name for the whole of it. The
SIGKILL was therefore never sent, and the workers the function exists to
remove kept running and holding cores — while the caller's blanket
``except Exception`` logged it at debug level and moved on.

Found by running ruff's F823 over the tree, which is why the lint job exists.
"""

import subprocess
import sys

import pytest

import delfin.orca as orca


@pytest.fixture
def sleeper():
    """A child in its own process group, so killpg cannot reach the test."""
    proc = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(30)"],
        start_new_session=True,
    )
    yield proc
    try:
        proc.kill()
        proc.wait(timeout=5)
    except Exception:  # noqa: BLE001 - best effort teardown
        pass


def test_the_fallback_kills_the_group_when_psutil_is_absent(monkeypatch, sleeper):
    monkeypatch.setattr(orca, "psutil", None)

    orca._ensure_process_group_terminated(sleeper, grace_timeout=0.2)

    # -15 is SIGTERM, -9 SIGKILL; either means the escalation ran through.
    assert sleeper.wait(timeout=5) in (-15, -9)


def test_the_fallback_does_not_raise_on_its_own_sleep(monkeypatch, sleeper):
    """The regression itself: the call must not blow up before it escalates."""
    monkeypatch.setattr(orca, "psutil", None)

    try:
        orca._ensure_process_group_terminated(sleeper, grace_timeout=0.2)
    except UnboundLocalError as exc:  # pragma: no cover - the bug being pinned
        pytest.fail(f"cleanup raised instead of killing the group: {exc}")


def test_the_psutil_failure_path_also_escalates(monkeypatch, sleeper):
    """The second fallback, reached when psutil itself throws, shares the
    SIGTERM/sleep/SIGKILL sequence and must be just as intact."""
    class _Exploding:
        @staticmethod
        def process_iter(*_args, **_kwargs):
            raise RuntimeError("psutil unavailable mid-scan")

    monkeypatch.setattr(orca, "psutil", _Exploding)

    orca._ensure_process_group_terminated(sleeper, grace_timeout=0.2)

    assert sleeper.wait(timeout=5) in (-15, -9)


def test_an_already_dead_process_is_not_an_error(monkeypatch):
    """Cleanup runs after ORCA exits, so the group is usually gone already."""
    proc = subprocess.Popen([sys.executable, "-c", "pass"], start_new_session=True)
    proc.wait(timeout=5)
    monkeypatch.setattr(orca, "psutil", None)

    orca._ensure_process_group_terminated(proc, grace_timeout=0.1)
