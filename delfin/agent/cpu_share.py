"""Run a command on a given set of CPUs.

Usage: ``python -m delfin.agent.cpu_share 2,5 <program> [args...]``.

The process binds itself (``os.sched_setaffinity``) and then replaces
itself with the program, so the program and everything it starts inherit
the binding. It stands in for ``taskset -c`` where util-linux is not
installed. On a system without the affinity call (macOS) the program runs
unbound.
"""
from __future__ import annotations

import os
import sys


def main(argv: list) -> None:
    cpus = {int(c) for c in str(argv[0]).split(",") if c.strip().isdigit()}
    program = list(argv[1:])
    setter = getattr(os, "sched_setaffinity", None)
    if setter is not None and cpus:
        try:
            setter(0, cpus)
        except OSError:
            pass                     # a CPU outside this process's set
    os.execvp(program[0], program)


if __name__ == "__main__":
    main(sys.argv[1:])
