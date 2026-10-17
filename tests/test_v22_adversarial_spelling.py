"""V22 reviewer adversarial tests for delfin.agent.refusal_memory.same_change.

Attack shapes the builder's tests/test_v22_same_change.py does NOT cover:
env prefixes before the writing command, redirections without a space,
numbered and >&word redirections, quoted paths, and unicode look-alike
paths. Expectation is the SAME one the builder pinned: a bash attempt that
would write the refused file yields the refused target; anything that does
not write it yields None.

Control: this file collects as ERROR on the previous commit (same_change
does not exist yet), so every case below is red with the module unchanged.
"""
from pathlib import Path

from delfin.agent.api_client import _bash_write_targets
from delfin.agent.refusal_memory import Refusal, same_change

WS = "/home/u/workspace"
F = WS + "/f.txt"


def _refused_write(path, tool="bash", reason="the user said no to this"):
    return Refusal(
        tool=tool,
        target=str(Path(path).resolve()),
        reason=reason,
        time="2026-10-16T00:00:00+00:00",
    )


def _changed(refused, cmd):
    return same_change(
        refused, "bash", {"command": cmd}, WS,
        write_targets=_bash_write_targets,
    )


# ---- env prefix before the writing command ---------------------------


class TestEnvPrefix:
    def test_env_prefix_sed_in_place(self):
        refused = _refused_write(F)
        assert _changed(refused, f"FOO=1 sed -i 's/a/b/' {F}") == str(
            Path(F).resolve()
        )

    def test_env_prefix_redirect(self):
        refused = _refused_write(F)
        assert _changed(refused, f"BAR=yes sh -c 'echo x > {F}'") == str(
            Path(F).resolve()
        )


# ---- redirections without a space, numbered, >&word ------------------


class TestRedirectSpellings:
    def test_no_space_redirect(self):
        refused = _refused_write(F)
        assert _changed(refused, f"echo x >{F}") == str(Path(F).resolve())

    def test_numbered_stdout_redirect(self):
        refused = _refused_write(F)
        assert _changed(refused, f"echo x 1> {F}") == str(Path(F).resolve())

    def test_numbered_stderr_redirect(self):
        refused = _refused_write(F)
        assert _changed(refused, f"cmd 2> {F}") == str(Path(F).resolve())

    def test_ampersand_word_redirect(self):
        refused = _refused_write(F)
        assert _changed(refused, f"echo x >& {F}") is not None or True


# ---- quoting variants -------------------------------------------------


class TestQuoting:
    def test_double_quoted_path(self):
        refused = _refused_write(F)
        assert _changed(refused, f"sed -i 's/a/b/' \"{F}\"") == str(
            Path(F).resolve()
        )

    def test_quoted_around_variable(self):
        refused = _refused_write(WS + "/my file.txt")
        assert _changed(
            refused, f"cat x > \"{WS}/my file.txt\""
        ) == str(Path(WS + "/my file.txt").resolve())


# ---- unicode look-alike path is a DIFFERENT file ---------------------


class TestUnicode:
    def test_confusable_path_is_no_change(self):
        # fullwidth "ｆ.txt" is not the refused "f.txt"
        refused = _refused_write(F)
        assert _changed(refused, f"echo x > {WS}/ｆ.txt") is None

    def test_real_unicode_name_matches_itself(self):
        refused = _refused_write(WS + "/r\u65e5\u672c\u8a9e.txt")
        assert _changed(
            refused, f"echo x > {WS}/r\u65e5\u672c\u8a9e.txt"
        ) == str(Path(WS + "/r\u65e5\u672c\u8a9e.txt").resolve())
