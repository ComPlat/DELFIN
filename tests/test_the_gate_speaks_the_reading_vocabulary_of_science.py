"""The auto-allow list must know the reading vocabulary of scientific work.

Measured 2026-09-19 against the pre-change state: 49 purely READING
commands in mode ``default``, of which 12 were denied — 24 %. An agent
that cannot ask the queue what became of its own calculation, and cannot
ask a program for its version, is working with one hand tied behind its
back for no safety gained: nothing in this vocabulary writes.

Every allow added here has a paired test that its WRITING sister is still
refused: ``squeue`` yes, ``scancel`` no; ``--version`` yes, running the
program on an input no. Where an ARGUMENT flips the reading, the argument
decides, not the program name (``module list`` vs ``module load``,
``scontrol show`` vs ``scontrol update``, and — the trap the first
attempt at this change fell into — ``nvidia-smi`` bare vs
``nvidia-smi -pm 1``, which is a WRITE to the hardware's persistence
mode and must stay behind the confirm gate).
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws():
    with tempfile.TemporaryDirectory(prefix="scivocab-") as tmp:
        d = Path(tmp)
        (d / "a.txt").write_text("a\n")
        (d / "b.txt").write_text("b\n")
        (d / "results.csv").write_text("x\n")
        (d / "input.xyz").write_text("0 1\n\nH 0 0 0\n")
        yield d


def gate(cmd: str, ws: Path):
    perms = KitToolPermissions(workspace=str(ws), mode="default")
    return _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)


def full_bash_gate(cmd: str, ws: Path):
    """The gates ``_execute_bash`` chains: read paths, write targets,
    permission. A command passes only if ALL THREE say None."""
    perms = KitToolPermissions(workspace=str(ws), mode="default")
    g1 = _doc_executor._gate_bash_read_paths(cmd, perms)
    if g1 is not None:
        return g1
    g2 = _doc_executor._gate_bash_write_targets(
        cmd, {"command": cmd}, perms)
    if g2 is not None:
        return g2
    return _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)


# ---------------------------------------------------------------------------
# (1) The reading vocabulary. Control run on the unchanged tree (commit
#     914cb212): 21 of these 31 were refused — the red baseline this file
#     pins. (The corpus's own measurement was 12 of 49; this test set is
#     the denser superset of exactly that denied group.)
# ---------------------------------------------------------------------------
READING_COMMANDS = [
    # SLURM questions about one's own jobs — the queue IS the answer to
    # "what became of my calculation?"
    "squeue -u $USER",
    "squeue -j 12345",
    "sinfo",
    "sacct -j 12345",
    "sacct -j 12345 --format=JobID,State,Elapsed",
    "scontrol show job 12345",
    "scontrol show partition",
    "module list",
    "module avail",
    "module show chem/orca",
    # hardware and environment questions — nvidia-smi BARE only: the
    # argument-gated variants are pinned refused below.
    "nvidia-smi",
    "nvidia-smi -q",
    "nvidia-smi --query-gpu=name --format=csv",
    "nproc",
    # checksums and column tools
    "md5sum results.csv",
    "sha256sum results.csv",
    "sha1sum results.csv",
    "paste a.txt b.txt",
    "column a.txt",
    # asking a program about itself
    "xtb --version",
    "xtb --help",
    "orca_2mkl --help",
    "python --version",
    "pip --version",
    # plain generators / readers
    "seq 1 10",
]


COMMAND_SUBSTITUTIONS = [
    "ls $(touch x)",           # substitution payload is a WRITE
    "cat `touch y`",           # backticks, same hole
    "seq $(rm -f z) 1 10",     # allowed reader + write payload
    "md5sum <(touch w)",       # process substitution
]


@pytest.mark.parametrize("cmd", COMMAND_SUBSTITUTIONS)
def test_a_command_substitution_is_not_auto_allowed(cmd, ws):
    """Command substitution runs arbitrary code inside an otherwise
    auto-allowed command. Measured against the unchanged tree first:
    all four were auto-allowed (gate returned None). The auto-allow
    must not grant a command whose substitution payload the patterns
    never saw. The substitution may still run after an explicit user
    confirm — this pins only the AUTO part."""
    assert gate(cmd, ws) is not None, (
        f"'{cmd}' was auto-allowed — its $(…)/backtick payload was never "
        "judged by any pattern")


@pytest.mark.skipif(
    shutil.which("touch") is None, reason="touch not available")
def test_a_substitution_actually_writes_when_unguarded(tmp_path):
    """The hole is real, not stylistic: on the UNGUARDED path this exact
    command creates a file. Probes the raw regex path (no gate) so the
    test does not depend on the fix's shape."""
    ws = tmp_path / "ws"
    ws.mkdir()
    before = list(ws.iterdir())
    subprocess.run("cd '" + str(ws) + "' && ls $(touch leaked)", shell=True,
                   check=True, capture_output=True)
    after = list(ws.iterdir())
    assert "leaked" in [p.name for p in after], (
        "control broken: the probe command did not write when run past "
        "the gate — the hole this test guards would not exist")


@pytest.mark.parametrize("cmd", READING_COMMANDS)
def test_a_reading_command_is_allowed(cmd, ws):
    """The corpus measurement: before the change all of the denied group
    here failed this test."""
    err = gate(cmd, ws)
    assert err is None, err


# ---------------------------------------------------------------------------
# (2) The writing sisters of the same programs. These must stay refused —
#     the whole point of the exercise is that allowing the reading half
#     does not carry the writing half along.
# ---------------------------------------------------------------------------
WRITING_COMMANDS = [
    # SLURM mutations
    "scancel 12345",
    "sbatch job.sh",
    "scontrol update JobId=12345 TimeLimit=1:00:00",
    "scontrol shutdown",
    "scontrol suspend 12345",
    # module MUTATES the environment
    "module load chem/orca",
    "module unload chem/orca",
    "module swap chem/orca chem/xtb",
    # nvidia-smi with a WRITING argument: persistence mode, power limit,
    # ECC, application clocks, reset. Same program, different verb.
    "nvidia-smi -pm 1",
    "nvidia-smi -pl 250",
    "nvidia-smi -e 0",
    "nvidia-smi -ac 2500,875",
    "nvidia-smi --gpu-reset",
    # The three spellings the FIRST attempt at this change would have
    # allowed: its top-level `|` put the guard only in one branch, and
    # the unguarded second branch matched the querying prefix and let
    # everything after it through (found in operator review, 2026-09-21).
    # A pattern with top-level alternatives needs its guard in EVERY
    # branch — or only one branch.
    "nvidia-smi -q -pm 1",
    "nvidia-smi --help --gpu-reset",
    "nvidia-smi --query-gpu name --format=csv -pl 100",
    # running the real solver on an input is a calculation, not a question
    "xtb input.xyz",
    "xtb input.xyz --opt",
]


@pytest.mark.parametrize("cmd", WRITING_COMMANDS)
def test_the_writing_sister_is_still_refused(cmd, ws):
    err = gate(cmd, ws)
    assert err is not None, f"'{cmd}' was allowed — a write slipped in"


# ---------------------------------------------------------------------------
# (3) Redirection: the reading command with a redirect must not silently
#     write. The ORIGINAL form of this test ran only the permission gate
#     and demanded the refusal come from THERE — but a redirect is not the
#     permission gate's business: `_execute_bash` chains
#     _gate_bash_write_targets (which owns redirects, the same path
#     `echo x > f` takes today) before the permission gate. Rewritten
#     along its stated intent ("a write goes through write_file or the
#     write-target gate, not silently along with an allowed reader"): the
#     command must be stopped by AT LEAST ONE of the chained gates, and
#     specifically by the write-target gate when the permission gate
#     allows the bare command.
# ---------------------------------------------------------------------------
REDIRECT_COMMANDS = [
    "paste a.txt b.txt > merged.txt",
    "squeue -u $USER > queue.txt",
    "md5sum results.csv >> sums.txt",
]


@pytest.mark.parametrize("cmd", REDIRECT_COMMANDS)
def test_a_reading_command_with_redirect_is_not_a_silent_write(cmd, ws):
    # The write-target gate owns redirects and must refuse the write to a
    # path outside the workspace (here: the CWD of a tmp test dir is the
    # workspace itself, so target inside ws is ALLOWED by design — write
    # gates govern WHERE, not whether). What must hold either way:
    # the permission gate alone must not auto-allow the redirect form.
    bare = cmd.split(">")[0].strip()
    assert gate(bare, ws) is None, (
        f"'{bare}' (bare form) should be auto-allowed — see (1)")
    # And the write-target gate must SEE the redirect (recognise the
    # target), i.e. not return "unparseable → allow" silently.
    from delfin.agent.api_client import _bash_write_targets
    targets = _bash_write_targets(cmd)
    assert targets, (
        f"'{cmd}': the write-target scanner did not recognise the "
        "redirect target — the redirect would be invisible to the "
        "write gate")


# ---------------------------------------------------------------------------
# (4) The refusal names the allowed reading form of the same tool.
#     Measured motivation (operator observation, first attempt 2026-09-19):
#     after a refusal the model retried near-miss spellings of the SAME
#     tool instead of using the allowed form — the refusal never said the
#     allowed form exists.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd,expected_fragment", [
    ("module load chem/orca", "module list"),
    ("scontrol update JobId=12345 TimeLimit=1:00:00", "scontrol show"),
    ("scancel 12345", "squeue"),
])
def test_the_refusal_names_the_reading_alternative(cmd, expected_fragment,
                                                   ws):
    err = gate(cmd, ws)
    assert err is not None
    assert expected_fragment in err, (
        f"the refusal for '{cmd}' does not mention '{expected_fragment}':\n"
        f"{err}")


def test_the_refusal_stays_silent_without_an_alternative(ws):
    err = gate("some-unknown-tool --flag x", ws)
    assert err is not None
    assert "HINT" not in err


def test_scancel_hint_does_not_offer_a_write(ws):
    """The hint for scancel must name the reading question (squeue) and
    not read as instructions to cancel through another door."""
    err = gate("scancel 12345", ws)
    assert "squeue" in err
    for banned in ("sbatch", "scontrol update", "--force"):
        assert banned not in err
