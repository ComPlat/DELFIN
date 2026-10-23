"""V22W phase 1: a refused quoted path with a space is still recognized.

RED CONTROL: on refusal_memory.py as merged (same_change takes targets
only from _bash_write_targets), all three tests below fail: the quoted
redirect word `"/home/u/workspace/my file.txt"` is split at the space by
the gate tokenizer, so same_change returns None where the shell would
write the refused file. Red on commit e724cdd8 (verified with
tests/test_zz_quoted_path.py: targets=['/home/u/workspace/my']).
"""
from pathlib import Path

from delfin.agent.api_client import _bash_write_targets
from delfin.agent.refusal_memory import Refusal, same_change

WS = "/home/u/workspace"


def _refused(path):
    return Refusal(tool="bash", target=str(Path(path).resolve()),
                   reason="outside workspace", time="00:00")


def _changed(refused, cmd):
    return same_change(refused, "bash", {"command": cmd}, WS,
                       write_targets=_bash_write_targets)


class TestQuotedSpace:
    def test_quoted_redirect_space(self):
        f = WS + "/my file.txt"
        refused = _refused(f)
        assert _changed(refused, f'cat x > "{f}"') == str(Path(f).resolve())

    def test_quoted_redirect_space_single_quotes(self):
        f = WS + "/my file.txt"
        refused = _refused(f)
        assert _changed(refused, f"cat x > '{f}'") == str(Path(f).resolve())

    def test_no_quoting_still_matches(self):
        # The ordinary spelling must not regress.
        f = WS + "/plain.txt"
        refused = _refused(f)
        assert _changed(refused, f"cat x > {f}") == str(Path(f).resolve())
