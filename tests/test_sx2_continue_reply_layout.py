"""The supervisor's CONTINUE reply puts the flag where the kernel reads it.

struct seccomp_notif_resp is (u64 id, s64 val, s32 error, u32 flags). A
reply that let the kernel carry out the notified call with CONTINUE in the
*error* field made the command's socket() fail with EINVAL: OpenMPI's
pmix_ifinit inside the agent's cage, live run 2026-11-17. With CONTINUE
the kernel requires val == 0 and error == 0.
"""
from __future__ import annotations

import struct

import pytest

sg = pytest.importorskip("delfin.agent.socket_guard")


class _Capture:
    def __init__(self):
        self.raw = b""

    def ioctl(self, fd, req, buf):
        self.raw = bytes(buf.raw[:24])
        return 0


def _reply(monkeypatch, call):
    cap = _Capture()
    monkeypatch.setattr(sg, "_libc", cap)
    sup = sg.Supervisor.__new__(sg.Supervisor)
    sup.listener = 99
    call(sup)
    return struct.unpack("<QqiI", cap.raw)


def test_continue_sets_the_flags_field_and_nothing_else(monkeypatch):
    nid, val, err, flags = _reply(
        monkeypatch, lambda s: s._datagram_socket(7, 1234, (0,) * 6, 41))
    assert (nid, val, err, flags) == (7, 0, 0, sg._NOTIF_FLAG_CONTINUE)


def test_a_refusal_carries_its_errno_and_no_flag(monkeypatch):
    import errno
    nid, val, err, flags = _reply(monkeypatch, lambda s: s._send(8, -1, -errno.EACCES))
    assert (nid, val, err, flags) == (8, -1, -errno.EACCES, 0)


def test_a_notification_is_read_at_the_kernels_offsets():
    """struct seccomp_notif: id @0, pid @8, flags @12, seccomp_data @16
    (nr @16, arch @20, instruction pointer @24, args @32). Read at 24, the
    number was the instruction pointer and socket() was never recognised."""
    raw = struct.pack("<QIIiIQ6Q", 77, 4321, 0, 41, 0xC000003E, 0x7F00DEADBEEF,
                      2, 2, 0, 0, 0, 0) + b"\0" * 8
    nid, pid, nr, args = sg.parse_notification(raw)
    assert (nid, pid, nr) == (77, 4321, 41)
    assert args[:3] == (2, 2, 0)


def test_the_socket_number_reaches_the_datagram_judge(monkeypatch):
    """The arch's own socket number dispatches to _datagram_socket."""
    seen = []
    sup = sg.Supervisor.__new__(sg.Supervisor)
    monkeypatch.setattr(sup, "_datagram_socket", lambda *a: seen.append("socket"), raising=False)
    monkeypatch.setattr(sup, "_send", lambda *a: seen.append(("send",) + a), raising=False)
    nr_socket = sg._arch()[3]
    sup._handle(5, 1, (2, 2, 0, 0, 0, 0), nr_socket)
    assert seen == ["socket"], seen


def run_filter(program, nr, args, arch=None):
    """Evaluate a classic-BPF seccomp program on a seccomp_data built with
    the kernel's layout: nr @0, arch @4, instruction_pointer @8, args[i]
    @16+8*i (64-bit, little endian). Only the opcodes the guard uses."""
    arch = sg._arch()[0] if arch is None else arch
    a6 = list(args) + [0] * (6 - len(args))
    data = struct.pack("<iIQ6Q", nr, arch, 0x401000, *a6)
    acc, pc = 0, 0
    for _ in range(len(program) + 2):
        code, jt, jf, k = program[pc]
        if code == sg._RET:
            return k
        if code == sg._LD:
            (acc,) = struct.unpack_from("<I", data, k)
            pc += 1
        elif code == sg._AND:
            acc &= k
            pc += 1
        elif code == sg._JEQ:
            pc += 1 + (jt if acc == k else jf)
        elif code == sg._JSET:
            pc += 1 + (jt if acc & k else jf)
        else:
            raise AssertionError(f"opcode {code:#x} not modelled")
    raise AssertionError("program did not return")


def test_a_walk_of_the_real_program_with_the_kernel_layout():
    """socket(AF_INET, SOCK_DGRAM) -> notify; SOCK_STREAM -> allow;
    unix datagram -> EACCES; sendmsg/sendmmsg -> EACCES; connect -> notify."""
    import errno
    import socket as so
    p = sg._program(block_inet_datagrams=True)
    a = sg._arch()
    d = sg._arch_dgram()
    eacces = sg._RET_ERRNO | errno.EACCES
    assert run_filter(p, a[3], (so.AF_INET, so.SOCK_DGRAM, 0)) == sg._RET_USER_NOTIF
    assert run_filter(p, a[3], (so.AF_INET6, so.SOCK_DGRAM, 0)) == sg._RET_USER_NOTIF
    assert run_filter(p, a[3], (so.AF_INET, so.SOCK_STREAM, 0)) == sg._RET_ALLOW
    assert run_filter(p, a[3], (so.AF_UNIX, so.SOCK_DGRAM, 0)) == eacces
    assert run_filter(p, a[2], (3, 0x1000, 16)) == sg._RET_USER_NOTIF
    assert run_filter(p, d[2], (3, 0x1000, 0)) == eacces
    assert run_filter(p, d[3], (3, 0x1000, 1, 0)) == eacces
