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
