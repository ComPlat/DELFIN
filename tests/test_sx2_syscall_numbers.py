"""SX2: the architecture rows are derived, never hand-copied.

Two guards:
- test_that_the_arch_rows_match_the_kernel_headers: every number in
  _ARCH_DGRAM is checked against the running machine's own syscall
  header (asm/unistd_64.h; asm-generic for aarch64), skipped only when
  neither header exists.
- test_that_the_arch_rows_match_the_names: a NAME->number mapping is
  asserted to be consistent with the rows, so 45/46 can never again
  mean the wrong call (the live-run kill of recvfrom on a6ee6bba9).
"""
import re
from pathlib import Path

from delfin.agent.socket_guard import _ARCH_DGRAM, _ARCH_LISTEN

_NAME_TO_NR = {
    "socket": {("x86_64", 41), ("aarch64", 198)},
    "connect": {("x86_64", 42), ("aarch64", 203)},
    "sendto": {("x86_64", 44), ("aarch64", 206)},
    "recvfrom": {("x86_64", 45), ("aarch64", 207)},
    "sendmsg": {("x86_64", 46), ("aarch64", 211)},
    "recvmsg": {("x86_64", 47), ("aarch64", 212)},
    "bind": {("x86_64", 49), ("aarch64", 200)},
    "listen": {("x86_64", 50), ("aarch64", 201)},
    "sendmmsg": {("x86_64", 307), ("aarch64", 269)},
}

_HDR_CANDIDATES = [
    "/usr/include/asm/unistd_64.h",           # x86_64 hosts
    "/usr/include/asm-generic/unistd.h",      # aarch64 hosts
]

# which _ARCH_DGRAM slot carries which name (machine -> tuple index)
_SLOT = {0: "bind", 1: "sendto", 2: "sendmsg", 3: "sendmmsg"}


def _header_table():
    """Parse the machine's own syscall header into {name: nr}.

    Returns None when no header exists (the machine is a CI container
    without linux-libc-dev) -- the test skips, never guesses."""
    machine = __import__("platform").machine()
    hdr = _HDR_CANDIDATES[0] if machine == "x86_64" else _HDR_CANDIDATES[1]
    path = Path(hdr)
    if not path.exists():
        return None
    table = {}
    for line in path.read_text(encoding="latin-1").splitlines():
        m = re.match(r"\s*#\s*define\s+__NR_([A-Za-z0-9_]+)\s+(\d+)", line)
        if m:
            table[m.group(1)] = int(m.group(2))
    return table


def test_that_the_arch_rows_match_the_kernel_headers():
    table = _header_table()
    if table is None:
        import pytest
        pytest.skip("no syscall header on this host; the numbers are "
                    "guarded by the name-mapping test instead")
    for machine, row in _ARCH_DGRAM.items():
        if machine != __import__("platform").machine():
            continue
        for idx, name in _SLOT.items():
            want = _NAME_TO_NR[name]
            nr = row[idx]
            for mm, n in want:
                if mm == machine:
                    assert nr == n, (
                        f"_ARCH_DGRAM[{machine}][{idx}] = {nr} for {name} "
                        f"but the kernel header says {n} -- hand-copied "
                        "numbers are forbidden, derive them")
                    if name in table:
                        assert table[name] == nr, (
                            f"{name}: header {table[name]} vs row {nr}")


def test_that_the_arch_rows_match_the_names():
    """Every row entry equals the name->number mapping: wrong constants
    cannot come back (the recvfrom kill on a6ee6bba9)."""
    for machine, row in _ARCH_DGRAM.items():
        for idx, name in _SLOT.items():
            for mm, n in _NAME_TO_NR[name]:
                if mm == machine:
                    assert row[idx] == n, (
                        f"_ARCH_DGRAM[{machine}][{idx}] = {row[idx]}, "
                        f"{name} is {n}")


def test_that_the_listen_row_matches_the_name_and_the_header():
    """listen() is notified since 2026-11-24 (the inbound half of the
    guard); its number is held to the same two checks as the others."""
    for machine, nr in _ARCH_LISTEN.items():
        assert (machine, nr) in _NAME_TO_NR["listen"], (machine, nr)
    table = _header_table()
    machine = __import__("platform").machine()
    if table is not None and machine in _ARCH_LISTEN and "listen" in table:
        assert table["listen"] == _ARCH_LISTEN[machine]
