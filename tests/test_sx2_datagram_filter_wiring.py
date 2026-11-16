"""The socket filter's datagram rules (SX2 phase 2b).

The measured live failure (JC1/JC2/JC3 runs, 14.-15.11.2026): inside the
agent cage OpenMPI's interface discovery died with
``opal_ifinit: socket() failed with errno=13`` -- a UDP datagram socket
refused AT CREATION by ``socket_guard._program(block_inet_datagrams=True)``,
before any connect() the supervisor could judge. The rule has to move to a
creation-time decision judged by the supervisor.

These tests read the assembled BPF program itself (``_assemble`` returns
tuples, so they hold on every little-endian architecture and need no real
seccomp kernel). Assertions trace the real jump chains, not labels.
"""
import errno
import socket

from delfin.agent import socket_guard as SG

_LD, _JEQ, _RET, _AND = 0x20, 0x15, 0x06, 0x54
RET_ALLOW = 0x7FFF0000            # SECCOMP_RET_ALLOW, as the kernel defines it
RET_EACCES = 0x00050000 | errno.EACCES
RET_NOTIFY = 0x7FC00000
RET_ENOSYS = 0x00050000 | errno.ENOSYS
AF_UNIX, AF_INET, AF_INET6 = socket.AF_UNIX, socket.AF_INET, socket.AF_INET6
DGRAM, STREAM = socket.SOCK_DGRAM, socket.SOCK_STREAM


def _target(program, idx, which):
    """The instruction a jump's true (jt) or false (jf) edge lands on."""
    jt, jf = program[idx][1], program[idx][2]
    return program[idx + 1 + (jt if which else jf)]


def _first(program, code, k):
    return next(i for i, ins in enumerate(program)
                if ins[0] == code and ins[3] == k)


def _type_check(program, entry):
    """From a family comparison, the type comparison the branch reaches:
    LD 24; AND 15; JEQ <type>. The true edge lands LD 24 directly
    (jt counts from entry+1)."""
    t = entry + 1 + program[entry][1]
    assert program[t] == (_LD, 0, 0, 24), program[t]
    assert program[t + 1] == (_AND, 0, 0, 15), program[t + 1]
    return t + 2


def _type_jump(program, entry):
    """The JEQ <type> instruction of a family's branch, and its index.
    The branch lands directly on the JEQ SOCK_DGRAM (the type is loaded and
    masked before the family dispatch; the two inet families share the
    JEQ's two return instructions after it)."""
    t = entry + program[entry][1]
    return t, program[t]


# -- the measured defect -----------------------------------------------------

def test_restricted_inet_datagram_creation_is_handed_to_the_supervisor():
    """The target state the wiring must reach: under a restricted network
    the AF_INET[6] SOCK_DGRAM branch is NOTIFIED (the supervisor judges
    creation) instead of an in-kernel EACCES -- today (before the wiring
    commit) the filter ends it in EACCES, the measured errno-13 death of
    OpenMPI's opal_ifinit, so this test is RED on the unfixed tree and
    green only when creation is judged. The inet families' dispatch lands
    on the inet branch's own JEQ DGRAM (LD 24; AND 15; JEQ), NOT on the
    shared UNIX-branch JEQ whose true edge must stay EACCES."""
    p = SG._program(block_inet_datagrams=True)
    for family in (AF_INET, AF_INET6):
        entry = _first(p, _JEQ, family)
        j = _type_check(p, entry)
        assert p[j][0] == _JEQ and p[j][3] == DGRAM, p[j]
        assert _target(p, j, True)[3] == RET_NOTIFY, (
            f"{family!r} datagram creation denied in-kernel "
            f"({_target(p, j, True)[3]:#x}) -- opal_ifinit dies with "
            "errno 13 before any supervisor decision")
        assert _target(p, j, False)[3] == RET_EACCES, _target(p, j, False)


# -- the invariants that must survive ----------------------------------------

def test_restricted_inet_stream_creation_stays_allowed_in_kernel():
    """The proxy path (TCP streams) is unaffected by the change."""
    p = SG._program(block_inet_datagrams=True)
    for family in (AF_INET, AF_INET6):
        type_idx, type_jump = _type_jump(p, _first(p, _JEQ, family))
        assert type_jump[0] == _JEQ and type_jump[3] == DGRAM, type_jump
        assert _target(p, type_idx, False)[3] == RET_ALLOW, _target(p, type_idx, False)


def test_datagram_unix_sockets_stay_denied_in_kernel():
    """The pre-existing rule survives: datagram UNIX sockets stay EACCES in
    the kernel (sendto carries its own address, past the supervisor)."""
    p = SG._program(block_inet_datagrams=True)
    # In the restricted program the UNIX branch shares the JEQ DGRAM with
    # the inet families (9 -> 14, 10 -> 14). The JEQ's two edges are the
    # two return instructions after it; the EACCES denial is one of them.
    idx = _type_jump(p, _first(p, _JEQ, AF_UNIX))[0]
    edges = (idx + 1 + p[idx][1], idx + 1 + p[idx][2])
    # Each edge reloads args[1], ANDs it, and jumps on SOCK_DGRAM again;
    # that jump's two ends are the RET instructions: one is the EACCES
    # denial that keeps the rule in the kernel today, one allows the
    # stream side. Both must stay reachable after the change.
    ends = set()
    for x in edges:
        j = x + 2
        ends.add(j + 1 + p[j][1])
        ends.add(j + 1 + p[j][2])
    assert any(p[o][0] == _RET and p[o][3] == RET_EACCES for o in ends), (
        [p[o] for o in ends])
    assert any(p[o][0] == _RET and p[o][3] == RET_ALLOW for o in ends), (
        [p[o] for o in ends])


def test_connect_and_uring_rules_unchanged():
    """connect() is still handed to the supervisor; io_uring is still
    refused with ENOSYS (it connects without the system call)."""
    p = SG._program(block_inet_datagrams=True)
    audit, _sc, nr_connect, _ns, nr_uring, _st, _sm = SG._arch()
    assert _target(p, _first(p, _JEQ, nr_connect), True)[3] == RET_NOTIFY
    assert _target(p, _first(p, _JEQ, nr_uring), True)[3] == RET_ENOSYS


def test_open_mode_program_keeps_its_old_shape():
    """Outside a restricted network the program is byte-for-byte today's."""
    assert SG._program(block_inet_datagrams=False) == SG._program(False)
    p = SG._program(block_inet_datagrams=False)
    # In the open program the UNIX branch has its own JEQ DGRAM (12), whose
    # true edge is the EACCES denial and whose false edge allows.
    idx = 12
    assert p[idx][0] == _JEQ and p[idx][3] == DGRAM, p[idx]
    assert _target(p, idx, True)[3] == RET_EACCES
    assert _target(p, idx, False)[3] == RET_ALLOW
