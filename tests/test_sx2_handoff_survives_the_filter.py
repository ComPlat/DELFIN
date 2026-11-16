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


def test_sendmsg_sendmmsg_refused_in_kernel():
    """The msghdr destination is invisible to BPF; sendmsg/sendmmsg are
    in-kernel EACCES under a restricted network (s2's pin; the measured
    OpenMPI path needs socket+ioctl and TCP connect only). Walk-asserted
    on the REAL program: an unheld nr must END at EACCES, not merely
    lack a JEQ (QS lesson: assert the walk result, not the absence)."""
    p = SG._program(block_inet_datagrams=True)
    nr_sendmsg = SG._ARCH_DGRAM["x86_64"][2]
    nr_sendmmsg = SG._ARCH_DGRAM["x86_64"][3]
    for nr in (nr_sendmsg, nr_sendmmsg):
        jeq = next((i for i, x in enumerate(p)
                    if x[0] == SG._JEQ and x[3] == nr), None)
        assert jeq is not None, (
            f"nr {nr} has no JEQ: with no JEQ it cannot be refused in-kernel"
            " (the socket JEQ's false edge falls through to the family"
            " dispatch and ends in ALLOW -- the reopened escape s2 walked)")
        end_true = jeq + 1 + p[jeq][1]
        # the TRUE edge (nr matches) must end in EACCES. The false edge
        # continues into the socket/connect family dispatch -- normal
        # routing for other nrs, asserted by the wiring tests.
        assert p[end_true][3] == SG._RET_ERRNO | errno.EACCES, (nr, p[end_true])
        # and the nr is really reachable: it sits after the sendto
        # block, before socket()'s fall-through JEQ
        socket_idx = next(i for i, x in enumerate(p)
                          if x[0] == SG._JEQ and x[3] == 41)
        assert jeq < socket_idx, (nr, jeq, socket_idx)


# -- target 2: the handoff survives its own filter ----------------------------

def test_the_handoff_protocol_is_synchronised():
    """The operator's live-run deadlock (pidfd_getfd(-1) before the child
    has installed its filter): the handoff must (1) tell the parent the
    listener's fd NUMBER, (2) wait for the parent's go before exec --
    both through a plain pipe or the socketpair, never a filtered
    syscall. Asserted on the real child/parent source."""
    import inspect
    import delfin.agent.landlock_exec as LE
    src = inspect.getsource(LE._run_guarded)
    child = src.split("if pid == 0:")[1].split("child_end.close()")[0]
    # the child sends the listener NUMBER (an int, not just a marker)
    assert "listener" in child and (".send(" in child or "os.write" in child), (
        "the child does not report the listener fd number")
    # the child waits for the parent's go AFTER reporting
    assert child.find(".send(") < child.find(".recv(") or (
        child.find("os.write") < child.find("os.read")), (
        "the child does not wait for the parent's go")
    # the parent gets the number from the child, not an init value
    parent = src.split("child_end.close()", 2)[2].split("if go:")[0]
    assert "recv" in parent or "os.read" in parent, (
        "the parent fetches nothing from the child -- pidfd_getfd(-1) "
        "again")


# -- target 3: end to end, the cage still opens -------------------------------

def test_a_caged_true_returns_quickly():
    """The real landlock_exec path with a restricted network: a caged
    `true` returns within seconds -- the handoff no longer notifies
    itself. Skipped where the cage cannot be nested (login node):
    a nested seccomp filter fails with EPERM/EACCES at install."""
    import os
    if os.environ.get("DELFIN_AGENT_SESSION"):
        pytest.skip("nested cage on this login node")
    try:
        import delfin.agent.socket_guard as SG
        SG.install_filter(block_inet_datagrams=True)
        raise RuntimeError("install_filter unexpectedly succeeded")
    except (PermissionError, OSError) as exc:
        if exc.errno in (errno.EACCES, errno.EPERM, errno.EBUSY):
            pytest.skip(f"nested seccomp not possible here ({exc})")
        raise
    except RuntimeError:
        pass
    import delfin.agent.landlock_exec as LE
    start = time.monotonic()
    rc = LE._run_guarded(write=[], hide=[], allowed=[],
                         command=["/bin/true"], strict=True,
                         use_fs=False, net_mode="proxy",
                         proxy_socket="", allow_tcp=[],
                         deny_socket=[], read_only=[], read=[])
    assert rc == 0, f"caged true returned {rc}"
    assert time.monotonic() - start < 10, "the cage handoff still blocks"
