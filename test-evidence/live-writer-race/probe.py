"""Measure what a reader sees while a writer is appending.

The audit log is appended to all day. ``append`` writes the line and then
the head file, so a reader can arrive between the two and find the log a
few bytes longer than its head file says. Before this change ``verify``
called that tampering; ``_just_grew`` accepts the gap only when the bytes
past the marker continue the chain.

Input: a log directory, a writer thread appending records, a reader
calling ``verify`` and counting reports that name a ``head`` problem.
Output: reads, not-intact reports, how many of them were ``head``, and
the module path that was measured -- so the two arms cannot be confused.

Not a test: the count depends on thread scheduling on the host, which is
why the committed test constructs the same state deterministically
instead. Run it from inside a worktree to measure that worktree.
"""
from __future__ import annotations

import collections
import pathlib
import sys
import tempfile
import threading
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from delfin.agent import audit_log as A  # noqa: E402


def measure(seconds: float = 3.0) -> tuple[int, int, collections.Counter]:
    directory = pathlib.Path(tempfile.mkdtemp())
    path = directory / "audit.log"
    stop = threading.Event()

    def writer() -> None:
        i = 0
        while not stop.is_set():
            A.append(A.make_record(tool="write_file", decision="ok",
                                   path=f"f{i}.py", session_id="s"),
                     log_path=path)
            i += 1

    thread = threading.Thread(target=writer, daemon=True)
    thread.start()
    reads = not_intact = 0
    kinds: collections.Counter = collections.Counter()
    deadline = time.monotonic() + seconds
    try:
        while time.monotonic() < deadline:
            rep = A.verify(path)
            reads += 1
            if not rep.get("ok"):
                not_intact += 1
                for problem in rep.get("problems") or []:
                    kinds[str(problem.get("kind"))] += 1
    finally:
        stop.set()
        thread.join(timeout=5)
    return reads, not_intact, kinds


if __name__ == "__main__":
    print("module:", A.__file__)
    reads, not_intact, kinds = measure()
    print(f"reads={reads} not_intact={not_intact} kinds={dict(kinds)}")
