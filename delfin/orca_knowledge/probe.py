"""Probe ORCA's own parser: is a word a simple-input keyword?

For each word a one-line input ``! WORD`` on H2 is written into a fresh
temporary directory and ORCA is started on it. If ORCA prints
``UNRECOGNIZED OR DUPLICATED KEYWORD`` the word is rejected; if parsing
finishes without that message it is accepted. Without an ORCA binary the
answer is ``unavailable``.

Only one ORCA runs at a time (``probe_all`` is strictly serial). Each run
gets its own process group; on timeout that group -- and only it -- is
killed.
"""
from __future__ import annotations

import glob
import os
import re
import shutil
import signal
import subprocess
import tempfile
from contextlib import contextmanager

#: Message ORCA's parser prints for a word it does not accept.
REJECT_MARKER = "UNRECOGNIZED OR DUPLICATED KEYWORD"
#: Printed once the input has been read: past this point the parser accepted it.
PARSED_MARKER = "INPUT FILE"
#: What a keyword may contain; anything else (a newline above all) could add
#: input lines of its own to the probe.
_WORD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9+\-*/(),_.'|]{0,79}$")

DEFAULT_TIMEOUT_S = 6

INPUT_TEMPLATE = """\
! {word}
* xyz 0 1
H 0.0 0.0 0.0
H 0.0 0.0 0.74
*
"""

# Overrides find_orca(); set by tests to a fake script. None = look it up.
ORCA_BINARY_CACHE: str | None = None


def get_orca() -> str | None:
    """The binary to use: the test/override cache, else find_orca()."""
    return ORCA_BINARY_CACHE if ORCA_BINARY_CACHE is not None else find_orca()


def find_orca() -> str | None:
    """Path of an ORCA binary: $ORCA_BINARY, else `orca` on PATH, else the
    newest /opt/orca*/orca; None when none exists."""
    env = os.environ.get("ORCA_BINARY")
    if env:
        return env
    found = shutil.which("orca")
    if found:
        return found
    candidates = sorted(glob.glob("/opt/orca*/orca"))
    return candidates[-1] if candidates else None


def _write_input(run_dir: str, word: str) -> str:
    inp = os.path.join(run_dir, "probe.inp")
    with open(inp, "w") as fh:
        fh.write(INPUT_TEMPLATE.format(word=word))
    return inp


@contextmanager
def _fresh_dir():
    """A fresh temporary directory, removed afterwards."""
    d = tempfile.mkdtemp(prefix="orca_probe_")
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _kill_group(pgid: int) -> None:
    try:
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def probe_one(word: str, run_dir: str | None = None,
              timeout_s: int = DEFAULT_TIMEOUT_S) -> str:
    """Ask ORCA's parser about one word -> accepted | rejected | unavailable."""
    if not _WORD.match(str(word or "")):
        return "rejected"
    binary = get_orca()
    if not binary or not os.access(binary, os.X_OK):
        # A missing/non-executable binary must not look like 'accepted':
        # `nice` would exit 127 with empty output otherwise.
        return "unavailable"
    with _fresh_dir() as fresh:
        if run_dir is not None:  # test hook: keep the run dir observable
            fresh = str(run_dir)
            os.makedirs(fresh, exist_ok=True)
        inp = _write_input(fresh, word)
        # Own process group so the timeout kills ORCA and its children only.
        try:
            proc = subprocess.Popen(
                ["nice", "-n", "19", binary, inp],
                cwd=fresh, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        except OSError:
            return "unavailable"
        pgid = proc.pid  # start_new_session makes pid == pgid
        output = b""
        timed_out = False
        try:
            output, _ = proc.communicate(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill_group(pgid)
            try:
                output, _ = proc.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                pass
        text = output.decode("utf-8", errors="replace")
        if REJECT_MARKER in text:
            return "rejected"
        if timed_out and PARSED_MARKER not in text:
            # Cut off before ORCA had even read the input: no verdict.
            return "unavailable"
        # Read without a refusal (a valid method may then run past the budget).
        return "accepted"


def probe_all(words: list[str], run_dir: str | None = None,
              timeout_s: int = DEFAULT_TIMEOUT_S) -> tuple[list[str], int]:
    """Probe words strictly one at a time (never more than one ORCA)."""
    results = [probe_one(w, run_dir, timeout_s) for w in words]
    rc = 2 if all(r == "unavailable" for r in results) else 0
    return results, rc
