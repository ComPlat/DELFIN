"""SX2 phase 2b fix 2: the datagram escape is judged, not opened.

Operator findings on 7c75a2796 (all three verified in code):
1. _NR_SOCKET/_NR_BIND were x86_64 clone/shutdown, not socket/bind -- the
   notified socket() fell through to the connect() judge. Target state:
   the numbers come from _arch(), never module constants.
2. bind() is never notified by the filter -- dead branch in the judge.
   Target state: bind() is notified only where an owned loopback endpoint
   must be distinguished (kept out of the filter; the judge stays for
   testability and future wiring) -- OR simpler: removed. Pinned here as
   dead-code-free.
3. SECURITY: socket(AF_INET/6, SOCK_DGRAM) allowed unconditionally, but
   sendto()/sendmsg()/sendmmsg() carry a destination address and leave
   past the supervisor (the proxy sees only connect()). Target state:
   the address-bearing datagram sends are notified and their destination
   is judged by the UDP ownership rule; a send without an address uses
   the connected peer connect() already judged.
"""
import errno
import socket
import struct

import pytest

from delfin.agent import socket_guard as SG


# -- finding 1: the arch's own numbers ---------------------------------------

def test_handle_dispatches_on_the_archs_own_numbers():
    """socket()/bind() reach the datagram judge by the numbers _arch()
    reports (x86_64 socket=41, bind=49; aarch64 socket=198, bind=200) --
    the fixed constants 56/48 were clone/shutdown."""
    audit, _sc, _nc, nr_socket, _nu, _po, _pg, *_rest = SG._arch()
    sup = SG.Supervisor(0, [])
    s = SG.Supervisor._datagram_socket
    calls = []
    sup._datagram_socket = lambda *a: calls.append(a)
    sup._handle(1, 1, (0, 0, 0), nr_socket)
    assert calls, "socket() by the arch's number did not reach the judge"
    calls.clear()
    sup._handle(1, 1, (0, 0, 0), nr=9999)          # not a judged call
    assert not calls


def test_arch_reports_socket_and_bind():
    """The architecture tables carry the datagram surface: socket() next
    to connect in _ARCHES, bind/sendto/sendmsg/sendmmsg in _ARCH_DGRAM."""
    for machine, entry in SG._ARCHES.items():
        assert entry[3] == {"x86_64": 41, "aarch64": 198}[machine], entry
        d = SG._ARCH_DGRAM[machine]
        assert d == {"x86_64": (49, 44, 45, 46),
                     "aarch64": (200, 205, 211, 212)}[machine], (entry, d)


# -- finding 2: no dead bind branch ------------------------------------------

def test_the_bind_judge_branch_is_not_dead_code():
    """The filter never notifies bind(); a judge branch that no
    notification can reach is dead code. Pinned: _handle has no bind
    dispatch (the filter's bind stays in-kernel-allowed, its datagram
    escape is judged at sendto)."""
    assert not hasattr(SG.Supervisor, "_bind"), (
        "a bind judge with no notifying filter is dead code")


# -- finding 3: the address-bearing datagram sends are judged ----------------

def test_sendto_with_a_destination_is_notified():
    """sendto with a destination address is handed to the supervisor
    under a restricted network -- a UDP datagram otherwise leaves past
    the supervisor (the proxy sees only connect()). sendmsg/sendmmsg
    are refused in-kernel since the deadlock fix (the destination is
    inside the msghdr, invisible to BPF); only sendto is notified."""
    p = SG._program(block_inet_datagrams=True)
    # sendto must have a JEQ in the restricted program
    jeqs = {i[3] for i in p if i[0] == SG._JEQ}
    from delfin.agent.socket_guard import _ARCH_DGRAM
    x86 = _ARCH_DGRAM["x86_64"]
    assert x86[1] in jeqs, "sendto not notified"


def test_a_send_without_the_arch_numbers_is_not_dispatched():
    """A non-send notification still goes to the connect judge."""
    sup = SG.Supervisor(0, [])
    calls = []
    sup._datagram_socket = lambda *a: calls.append(a)
    sup._datagram_send = lambda *a: calls.append(("send", a))
    sup._handle(1, 1, (0, 0, 0), nr=9999)
    assert not calls
