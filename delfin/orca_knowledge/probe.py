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


# --------------------------------------------------------------------------
# Full-input probe: ORCA itself decides.
# --------------------------------------------------------------------------
#: Marks printed by ORCA's parser for input it refuses (measured against real
#: ORCA 6.1.1, see PR_NOTES.md; the forms without "INPUT ERROR" match the
#: ORCA 6.1 manual, contents/quickstartguide/troubleshooting.html 1.7.4):
REJECT_MARKERS = (
    "UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE",  # 1.7.4.1
    "INPUT ERROR",
    "Unknown identifier in",      # unknown block option, 1.7.4.2
    "Unknown identifier %",       # unknown block name (qcscan line form)
    "Unrecognized symbol in",     # missing end, 1.7.4.4
    "ERROR: expect a '$', '!', '%', '*' or '[' in the input",  # 1.7.4.5/4.6
    "expected an identifier after '%%' in the input",          # 1.7.4.6
    "ERROR: non-ASCII character(s) found",                     # 1.7.4.9
    "CANNOT OPEN FILE",           # 1.7.4.7
    "Error: Quotes expected",     # 1.7.4.6 (%moinp without quotes)
)
#: Printed once the input has been read completely (end of the input echo).
PARSED_INPUT_MARKER = "INPUT FILE"
#: Printed only when a run completed to the very end.
FINISHED_MARKER = "ORCA TERMINATED NORMALLY"

MAX_ATOMS = 6
MAX_NPROCS = 2
MAX_MAXCORE = 2000
DEFAULT_INPUT_TIMEOUT_S = 60


class ProbeGuardError(ValueError):
    """The input is not small enough to be a probe."""


#: Blocks a probe may use: settings of a calculation, nothing that names a
#: file or a program.
ALLOWED_BLOCKS = frozenset({
    "scf", "method", "basis", "mdci", "casscf", "mrci", "mp2", "tddft", "cis", "rel",
    "grid", "freq", "geom", "pal", "maxcore", "output", "cpcm", "elprop", "eprnmr",
    "loc", "rocis", "esd", "nbo", "xtb", "neb", "irc", "goat", "autoci", "ice", "cosmors",
    "symmetry", "coords", "paras", "frag", "solvator", "docker", "md", "plots", "eda",
})
#: Anything that starts a program or reads a named file.
_FORBIDDEN = re.compile(
    r"(?i)\b(?:sys_?cmd|extopt|progext|ext_params|xyzfile|gzmtfile|pdbfile|moread|moinp|"
    r"inhessname|gtoname|gbwname|hessname|readfile|includefile|external|orca_\w+)\b"
    r"|%compound")


def _guard_shape(text: str) -> None:
    if _FORBIDDEN.search(text):
        raise ProbeGuardError(
            f"{_FORBIDDEN.search(text).group(0)!r} is not allowed in a probe: it runs a "
            "program or reads a file")
    if '"' in text or "'" in text:
        raise ProbeGuardError("quoted strings are not allowed in a probe (file names)")
    for m in re.finditer(r"(?im)^[ \t]*%(\w+)", text):
        if m.group(1).lower() not in ALLOWED_BLOCKS:
            raise ProbeGuardError(f"block %{m.group(1)} is not allowed in a probe")
    for m in re.finditer(r"(?im)^[ \t]*\*[ \t]*(\w+)", text):
        if m.group(1).lower() not in ("xyz", "int", "gzmt"):
            raise ProbeGuardError(f"geometry '* {m.group(1)}' is not allowed in a probe "
                                  "(inline coordinates only)")


def guard_input(text: str) -> None:
    """Refuse inputs that are more calculation than probe.

    Checks the whole input for coordinate-like lines (element + three
    numbers), ``%pal nprocs`` and ``%maxcore`` against the probe limits.
    Counting atom lines anywhere (not per coordinate block) is deliberately
    conservative for a probe: extra atoms may sit outside the first block.
    Raises :class:`ProbeGuardError` with the violated limit.

    It is also an allow-list, because the inputs come from agents and the
    probe runs outside their cage: ORCA can start programs from an input
    (``Sys_cmd`` in compound scripts, ``ExtOpt`` with ``ProgExt``) and read
    any file it is pointed at (``*xyzfile``, ``%moinp``, ``InHessName``, a
    basis file). So a probe may hold only simple-input lines, blocks from
    ``ALLOWED_BLOCKS`` and an inline ``* xyz`` geometry -- no compound
    script, no external program, no file reference, no quoted string.
    """
    _guard_shape(text)
    natoms = sum(
        1 for line in text.splitlines()
        if re.match(r"^\s*[A-Za-z]{1,3}\s+[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+"
                    r"[-+0-9.eEdD]+\s*$", line)
        and not line.strip().startswith("#"))
    if natoms > MAX_ATOMS:
        raise ProbeGuardError(
                f"{natoms} atoms > probe limit {MAX_ATOMS} (it is a probe, "
                "not a calculation)")
    for m in re.finditer(r"(?i)%pal\s*\n?\s*nprocs\s+(\d+)", text):
        if int(m.group(1)) > MAX_NPROCS:
            raise ProbeGuardError(
                f"%pal nprocs {m.group(1)} > probe limit {MAX_NPROCS}")
    for m in re.finditer(r"(?i)%maxcore\s+(\d+)", text):
        if int(m.group(1)) > MAX_MAXCORE:
            raise ProbeGuardError(
                f"%maxcore {m.group(1)} > probe limit {MAX_MAXCORE}")


def look_for_hits(text: str, texts: list[str]) -> list[dict]:
    """For each searched text: found? plus the lines around every occurrence."""
    lines = text.splitlines()
    hits = []
    for wanted in texts:
        positions = [i for i, line in enumerate(lines) if wanted in line]
        context = []
        for i in positions[:5]:  # report the first five occurrences
            lo, hi = max(0, i - 2), min(len(lines), i + 3)
            context.append("\n".join(lines[lo:hi]))
        hits.append({"text": wanted, "found": bool(positions),
                     "count": len(positions), "context": context})
    return hits


#: Printed when ORCA stops during setup, after the input was accepted
#: (e.g. "WARNING: For DLPNO calculations we need at least one auxiliary
#: basis set" -> "Error (ORCA_MAIN): ... aborting the run"; measured with
#: real ORCA 6.1.1). The input was read without a refusal -> "parsed",
#: though the "INPUT FILE" echo itself was not reached.
ABORT_MARKER = "aborting the run"


def _verdict_for_input(text: str, timed_out: bool) -> tuple[str, list[str]]:
    """(verdict, messages) for one full-input run."""
    messages: list[str] = []
    for marker in REJECT_MARKERS:
        if marker in text:
            for line in text.splitlines():
                if marker in line or (
                        messages and line.strip().startswith("Last token")):
                    messages.append(line.strip())
    if messages:
        return "rejected", messages[:10]
    if ABORT_MARKER in text:
        return "parsed", []
    if timed_out:
        if PARSED_INPUT_MARKER in text:
            return "parsed", []  # read, but ran past the budget
        return "timeout", []
    if FINISHED_MARKER in text:
        return "finished", []
    if PARSED_INPUT_MARKER in text:
        return "parsed", []
    return "timeout", []


def probe_input(text: str, look_for: list[str] | None = None,
                timeout_s: int = DEFAULT_INPUT_TIMEOUT_S,
                run_dir: str | None = None) -> dict:
    """Run ORCA on a complete small input and let ORCA decide the verdict.

    Returns a dict: ``input`` (echo), ``verdict``
    (rejected|parsed|finished|timeout|unavailable), ``messages`` (ORCA's own
    lines when rejected), ``look_for`` (per searched text: found/count/context),
    ``timed_out``.
    """
    result: dict = {"input": text, "verdict": "unavailable", "messages": [],
                    "look_for": [], "timed_out": False}
    try:
        guard_input(text)
    except ProbeGuardError as exc:
        result["verdict"] = "refused"
        result["messages"] = [str(exc)]
        return result
    binary = get_orca()
    if not binary or not os.access(binary, os.X_OK):
        return result
    look_for = list(look_for or [])
    with _fresh_dir() as fresh:
        if run_dir is not None:  # test hook: keep the run dir observable
            fresh = str(run_dir)
            os.makedirs(fresh, exist_ok=True)
        inp = os.path.join(fresh, "probe_input.inp")
        with open(inp, "w") as fh:
            fh.write(text)
        try:
            proc = subprocess.Popen(
                ["nice", "-n", "19", binary, inp],
                cwd=fresh, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                start_new_session=True)
        except OSError:
            return result
        pgid = proc.pid
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
        text_out = output.decode("utf-8", errors="replace")
        verdict, messages = _verdict_for_input(text_out, timed_out)
        result.update(verdict=verdict, messages=messages, timed_out=timed_out,
                      look_for=look_for_hits(text_out, look_for))
    return result

