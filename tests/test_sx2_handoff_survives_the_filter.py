"""SX2 phase 2b deadlock fix: the cage's own handoff must survive the
filter.

The live run on f4217992a found the cage deadlocks before /bin/bash:
landlock_exec.py:465 hands the supervisor listener with
socket.send_fds (= sendmsg), and the restricted program notified
sendmsg -- that notification goes to USER_NOTIF with nobody serving
yet. Target state, per the operator's directions + s2's recommendation:

1. sendto: notify ONLY when the destination pointer (args[4]) is
   non-NULL; a NULL destination uses the connected peer (already
   judged) and stays in-kernel ALLOW. Cheap BPF test, supervisor out
   of the no-destination hot path.
2. sendmsg/sendmmsg: back to in-kernel EACCES under a restricted
   network (measurement: OpenMPI's interface read is socket+ioctl, the
   rank traffic is TCP connect; s2's pin). The supervisor handoff stops
   using sendmsg entirely: the parent dups the listener with
   pidfd_getfd and the child no longer sends fds at all.
3. The handoff must work end to end: a caged `true` returns within a
   few seconds (real landlock_exec).
"""
import errno
import socket
import struct
import time

import pytest

from delfin.agent import socket_guard as SG


# -- target 1: sendto judged by destination presence --------------------------

def _x86_sendto():
    return SG._ARCH_DGRAM["x86_64"][1]                # 44


def _type_check(p, entry):
    """Follow a family JEQ's true edge: LD 24, AND 15, then the type JEQ."""
    t = entry + 1 + p[entry][1]
    assert p[t] == (SG._LD, 24), p[t]
    assert p[t + 1] == (SG._AND, 0xF), p[t + 1]
    return t + 2


def test_sendto_notified_only_with_a_destination_pointer():
    """sendto with dest_addr != NULL is notified; with a NULL pointer it
    falls through to allow (connected peer, already judged)."""
    p = SG._program(block_inet_datagrams=True)
    nr_sendto = _x86_sendto()
    send_idx = next(i for i, x in enumerate(p)
                    if x[0] == SG._JEQ and x[3] == nr_sendto)
    ins = p[send_idx]
    # true edge (dest != NULL): LD args[1] then JEQ 0: NULL -> allow
    # (peer judged), non-NULL -> notify; false edge (not sendto): the
    # socket() checks
    jt = send_idx + 1 + ins[1]
    assert p[jt] == (SG._LD, 0, 0, 40), p[jt]    # load the dest pointer
    jeq = p[jt + 1]
    assert jeq[0] == SG._JEQ and jeq[3] == 0, jeq
    to_null = jt + 2 + jeq[1]                    # JEQ==0 true edge
    to_dest = jt + 2 + jeq[2]                    # JEQ!=0 false edge
    assert p[to_null] == (SG._RET, 0, 0, SG._RET_ALLOW), p[to_null]
    assert p[to_dest] == (SG._RET, 0, 0, SG._RET_USER_NOTIF), p[to_dest]


def test_sendmsg_sendmmsg_back_in_kernel():
    """The msghdr destination is invisible to BPF; sendmsg/sendmmsg are
    in-kernel EACCES under a restricted network (s2's pin; the measured
    OpenMPI path needs socket+ioctl and TCP connect only)."""
    p = SG._program(block_inet_datagrams=True)
    jeq_targets = {x[3] for x in p if x[0] == SG._JEQ}
    nr_sendmsg = SG._ARCH_DGRAM["x86_64"][2]
    nr_sendmmsg = SG._ARCH_DGRAM["x86_64"][3]
    assert nr_sendmsg not in jeq_targets
    assert nr_sendmmsg not in jeq_targets
    # and their JEQs are gone from the program entirely
    for nr in (nr_sendmsg, nr_sendmmsg):
        assert not any(x[0] == SG._JEQ and x[3] == nr for x in p), nr


# -- target 2: the handoff survives its own filter ----------------------------

def test_the_handoff_path_no_longer_sends_fds():
    """landlock_exec's child branch hands the listener without sendmsg:
    no send_fds call in the source, and the parent fetches the fd with
    pidfd_getfd."""
    import inspect
    import delfin.agent.landlock_exec as LE
    src = inspect.getsource(LE)
    # the comment names the removed mechanism; assert it is gone from
    # CODE, not prose: no call site remains
    assert "send_fds(" not in src, (
        "the child still hands the listener via sendmsg -- it notifies "
        "itself to nobody and deadlocks")
    assert "pidfd_getfd" in src or "PIDFD_GETFD" in src, (
        "the parent must fetch the listener itself")


# -- target 3: end to end, the cage still opens -------------------------------

def test_a_caged_true_returns_quickly():
    """The real landlock_exec path with a restricted network: a caged
    `true` returns within seconds -- the handoff no longer notifies
    itself. Skipped where the cage cannot be nested (login node)."""
    import os
    if os.environ.get("DELFIN_AGENT_SESSION"):
        pytest.skip("nested cage on this login node")
    import delfin.agent.landlock_exec as LE
    start = time.monotonic()
    rc = LE._run_guarded(write=[], hide=[], allowed=[],
                         command=["/bin/true"], strict=True,
                         use_fs=False, net_mode="proxy",
                         proxy_socket="", allow_tcp=[],
                         deny_socket=[], read_only=[], read=[])
    assert rc == 0, f"caged true returned {rc}"
    assert time.monotonic() - start < 10, "the cage handoff still blocks"
