"""V22W CI follow-up (PR #308 red): the quote-aware supplement must not
ignore a cd the command performs.

Defect: `cd sub && ./tool2 > f.txt` — the gate resolves the redirect
against `sub/` (correctly, no refusal match), but _redirect_targets
supplemented the bare `f.txt`, which resolves against the workspace and
falsely matched a refusal of `./f.txt` -> a refused write that was never
refused. Conservative rule (operator, 2026-10-17): when the command
changes directory (cd/pushd in any segment), take NO quote-aware
candidate — the pass only ever ADDS candidates, so skipping is safe.
"""
import pytest

from delfin.agent import refusal_memory as R


def test_cd_command_gets_no_quote_aware_candidate():
    refusal = R.Refusal("bash", "/ws/f.txt", "outside workspace", "00:00")
    # Quoted relative write after cd: the gate says sub/f.txt, the pass
    # must NOT add workspace-level f.txt.
    args = {"command": 'cd sub && ./tool2 > "my f.txt"'}
    got = R.same_change(refusal, "bash", args, "/ws",
                        write_targets=lambda cmd: ["sub/my f.txt"])
    assert got is None


def test_without_cd_the_quoted_candidate_still_arrives():
    refusal = R.Refusal("bash", "/ws/my f.txt", "outside workspace", "00:00")
    args = {"command": './tool2 > "my f.txt"'}
    got = R.same_change(refusal, "bash", args, "/ws",
                        write_targets=lambda cmd: ["my"])
    assert got == "/ws/my f.txt"


def test_cd_without_quotes_also_skips_the_pass():
    # Even unquoted: cd changes the resolution base, so the supplement
    # is skipped entirely (conservative in both spellings).
    refusal = R.Refusal("bash", "/ws/f.txt", "outside workspace", "00:00")
    args = {"command": "cd sub && ./tool2 > f.txt"}
    got = R.same_change(refusal, "bash", args, "/ws",
                        write_targets=lambda cmd: ["sub/f.txt"])
    assert got is None


def test_pushd_also_counts_as_changing_directory():
    refusal = R.Refusal("bash", "/ws/f.txt", "outside workspace", "00:00")
    args = {"command": "pushd sub && ./tool2 > f.txt"}
    got = R.same_change(refusal, "bash", args, "/ws",
                        write_targets=lambda cmd: ["sub/f.txt"])
    assert got is None
